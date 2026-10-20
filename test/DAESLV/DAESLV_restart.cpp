// DAESLV_restart -- the restart system is ONE parameterised system, not one per stop.
// =====================================================================================================
// Gate: a model with MANY stage boundaries must leave the model's dag the same size after a solve as
// before it, and give the same answer as the same model with one stage.  Rebuilding per stop was what
// added ~4 variables to FFModel::_dag per restart, so dag growth IS the instrument here.
//
//   y' = -y + u ,  0 = z - u y     with u piecewise constant over 8 elements, u_e = 1 + e/8
//
// The input changes at every boundary, so every stop must bind a DIFFERENT u into the same system.
#include <cmath>
#include <iomanip>
#include <iostream>
#include "daeslv_idas.hpp"

using namespace mc;

static int n_pass = 0, n_fail = 0;
static void check( std::string const& w, double got, double want, double tol )
{
  double const d = std::fabs( got - want );
  ( d <= tol )? ++n_pass: ++n_fail;
  std::cout << "    " << std::left << std::setw( 46 ) << w << std::right
            << " got=" << std::scientific << std::setprecision(10) << got
            << " want=" << want << " |d|=" << std::setprecision(2) << d
            << ( d <= tol? "  PASS\n": "  FAIL\n" );
}

int main()
{
  size_t const NE = 8;
  std::cout << "\n--- one parameterised restart system across " << NE << " elements\n";

  FFGraph G; DAESLV_IDAS M; M.set_dag( &G );
  FFVar t = G.add_var("t"), y = G.add_var("y(t)"), z = G.add_var("z(t)"), u = G.add_var("u");
  FFPartial OpP;  int const TI = FFDom::ALL - FFDom::LB;

  M.add_domain( t, FFDom( 0., 2., NE, FFDom::LGR, 4 ) );
  M.add_state( y, {t} );  M.add_state( z, {t} );
  M.add_input( u, {t}, FFDom::LGR, 1 );                      // piecewise constant: one level per element
  M.set_evolution_domain( t );
  M.update_ref( y, 1. );  M.update_ref( z, 1. );
  auto io = FFModel::EqnOptions( FFModel::EqnRole::INTERIOR, 0 );
  auto ii = FFModel::EqnOptions( FFModel::EqnRole::INITIAL,  0 );
  M.add_equation( OpP( y, t ) + y - u, {t}, {TI}, io );
  // FFDom::ALL, not TI: the ALGEBRAIC row must hold at t0 as well.  Declared on (0,T] it left z(0)
  // determined by nothing, and setup() said so -- "DOF: the model is NOT balanced -- rows - unknowns
  // = -1".  Every check here is about y, which IS determined, so the driver passed anyway: a gate
  // weaker than it looked, on a model the framework considered ill-posed.
  M.add_equation( z - u * y,           {t}, {FFDom::ALL}, io );
  M.add_equation( y - 1.,              {t}, {FFDom::LB}, ii );

  M.options.SOLVE.RTOL = M.options.SOLVE.ATOL = 1e-10;
  M.options.DISPLAY_LEVEL = 0;
  M.options.INITIAL.SOLVER = BASE_IDAS::Options::t_Initial::KINSOL;

  if( !M.setup() ){ std::cout << "setup failed: " << M.extract_error() << "\n"; return 1; }
  std::cout << "    stages=" << M.nx() << " states, np=" << M.np()
            << " (one level per element), ic_rows=" << M.ic_rows().size() << "\n";

  std::vector<double> P;
  for( size_t e = 0; e < NE; ++e ) P.push_back( 1. + (double)e / (double)NE );

  unsigned long const nv_before = G.nvar();
  auto const st = M.solve( P );
  unsigned long const nv_after = G.nvar();
  std::cout << "    solve -> status " << st << ", consistent initialisations = " << M.stats.numIC << "\n";
  std::cout << "    model dag variables: before=" << nv_before << " after=" << nv_after << "\n";
  auto const& xk_dbg = M.val_state();

  { double ye = 1.; std::cout << "    stage-by-stage y:  solver / exact\n";
    for( size_t e = 0; e < NE; ++e ){
      double const ue = P[e], h = 2./(double)NE;
      ye = ue + ( ye - ue ) * std::exp( -h );
      double const ys = xk_dbg.at(e+1).at(0);
      std::cout << "      after element " << e << "  u=" << std::fixed << std::setprecision(4) << ue
                << "   y=" << std::setprecision(10) << ys << "   exact=" << ye
                << "   d=" << std::scientific << std::setprecision(2) << std::fabs(ys-ye) << "\n";
    } }
  // GATE the balance, rather than leave it to a setup() warning that scrolls past under "all pass".
  // This driver ran for a week with rows - unknowns = -1 and every check still green, because every
  // check was about y and the deficit was in z.
  check( "the model is DOF-balanced", M.dof_balance().balanced? 1.: 0., 1., 0. );

  check( "the solve restarted at every boundary", (double)M.stats.numIC, (double)NE, 0. );
  check( "the model dag did NOT grow across the stops", (double)nv_after, (double)nv_before, 0. );

  // The algebraic state tracks the CURRENT level: z = u y on the last element, u = 1 + 7/8.
  auto const& xk = M.val_state();
  double const yT = xk.back()[0], zT = xk.back()[1];
  check( "z(T) = u_last * y(T)", zT, P.back() * yT, 1e-8 );

  // The CLOSED FORM, element by element: y' = -y + u with u constant on each element gives
  // y(t_{e+1}) = u_e + ( y(t_e) - u_e ) exp(-h).  This is the gate -- agreement between two correctors
  // would be weak evidence, and on this model it is also false (see below).
  { double ye = 1.;
    for( size_t e = 0; e < NE; ++e ){ double const ue = P[e]; ye = ue + ( ye - ue ) * std::exp( -2./(double)NE ); }
    check( "y(T) == the closed form, across 7 restarts", yT, ye, 1e-7 );
    check( "z(T) == u_last * the closed form",           zT, P.back() * ye, 1e-7 ); }

  // solve() must be IDEMPOTENT on a MULTI-STAGE model.  DAESLV_test1 checks idempotence on a single
  // stage, where the rows never change, so it could not see that the loop left _build_rows at the LAST
  // stage and a second solve started stage 0 with the previous run's input levels.  The error here was
  // 3.4e-02 at T, and it was reproducible -- a wrong answer, not a flaky one.
  M.solve( P );
  check( "solve() is idempotent on y(T)", M.val_state().back()[0], yT, 0. );
  check( "solve() is idempotent on z(T)", M.val_state().back()[1], zT, 0. );

  // The two correctors AGREE on this model, each from its own setup().  Switching INITIAL.SOLVER after
  // setup() is what used to differ, and that was the idempotence defect above, not the corrector:
  // IDACalcIC in YA_YDP_INIT tracks the closed form here to 5e-10.
  M.options.INITIAL.SOLVER = BASE_IDAS::Options::t_Initial::IDACALCIC;
  M.setup();
  M.solve( P );
  check( "IDACalcIC agrees with KINSOL on y(T)", M.val_state().back()[0], yT, 1e-7 );

  std::cout << "\n==========  DAESLV_restart: " << n_pass << " passed, " << n_fail << " failed -- "
            << ( n_fail? "FAILURES": "ALL PASS" ) << "  ==========\n\n";
  return n_fail? 1: 0;
}
