// DAESLV_dummy -- index reduction by DUMMY DERIVATIVES against SUBSTITUTION.
// =====================================================================================================
// REDUCE.DUMMY_DERIV keeps a high-index constraint as an interior row, adds its derivative as another,
// and mints one algebraic state per differentiation to stand for the derivative of a chosen state
// (Mattsson & Soderlind).  Substitution -- the default -- overwrites the constraint's row with its
// derivative and keeps the constraint only at the initial point, so on the interior it drifts with
// the integration error of its derivative (DAESLV_drift measures that).
//
// This driver gates the CASCADE: the index-3 pendulum, two dummies chosen at once under the
// nonsingularity test, and requires the two reductions to agree on the solve and the forward
// sensitivities, and the nonlinear constraint to hold TIGHTER under dummy derivatives.  The index-2
// case -- the vented reactor, with values, both sensitivities and the closed form tabulated -- is
// DAESLV_reduction's.  The horizon here stops before the bob crosses the axis where the static pivot
// vanishes: past it the solve stops by design (see DAESLV_IDAS's pivot diagnostic).
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
  std::cout << "    " << std::left << std::setw( 48 ) << w << std::right
            << " got=" << std::scientific << std::setprecision(8) << got << " want=" << want
            << " |d|=" << std::setprecision(2) << d << ( d <= tol? "  PASS\n": "  FAIL\n" ) << std::flush;
}
static void check_lt( std::string const& w, double got, double bound )
{
  ( got <= bound )? ++n_pass: ++n_fail;
  std::cout << "    " << std::left << std::setw( 48 ) << w << std::right
            << " got=" << std::scientific << std::setprecision(2) << got << " <= " << bound
            << ( got <= bound? "  PASS\n": "  FAIL\n" ) << std::flush;
}

struct Result { bool ok = false; size_t nx = 0; std::vector<double> f, gf, ga; double worst = 0.; };
// What the sensitivity paths do on a dummy-reduced model, as of 2026-10-25 -- the driver asserts only on
// what the method delivers, and names what it does not:
//   (index 2, the reactor: solve, fsens and asens all agree with substitution -- DAESLV_reduction)
//   TRANSITIONS              work by the ordinary index-1 rule (case C): a jumped state stays differential
//                             -- the pivot ranks it first -- and the algebraic states re-solve at tau+.
//   POINT OUTPUTS on an       seeded by projection through dF/d(y_d', y_a) at T, with the direct -F_p^T w
//   ALGEBRAIC state           term in the gradient (case C's nB(1)); formerly refused.
//   FSENS on the pendulum     works, to 2.5e-11 against substitution -- once the dummies are out of the
//                             error test.  A dummy is a derivative (dv = v' is an acceleration) and its
//                             BDF error estimate is meaningless: 1749 steps and a stalled fsens with it
//                             in, 179 steps with it out.  DAESLV excludes algebraics itself whenever the
//                             model carries dummies.
//   ASENS on the pendulum     works, to 5e-10 and 4e-8 against the forward.  It was wrong by ONE TERM of the
//                             ordinary index-1 adjoint: d/dt(F_y'^T lambda) = F_y'^T lambda' + (dF_y'/dt)^T
//                             lambda, and the second term -- zero for a constant mass matrix, every model
//                             in the corpus before this one -- was missing from the backward residual.
//                             Nothing about a dummy-reduced DAE is special for the adjoint: it is index 1.
// The flags below gate what works; set them true to exercise the gaps.
static bool const SENS_PENDULUM = true, ASENS_PENDULUM = true;

// ---------------------------------------------------------------- B: the pendulum, index 3
static Result pendulum( bool dummy, double const T )
{
  FFPartial OpP;  FFEval OpE;
  int const TI = FFDom::ALL - FFDom::LB;
  auto const INTER = FFModel::EqnRole::INTERIOR, INIT = FFModel::EqnRole::INITIAL;
  Result R;
  FFGraph G;  DAESLV_IDAS M;  M.set_dag( &G );
  FFVar t = G.add_var("t"), x = G.add_var("x(t)"), y = G.add_var("y(t)"), u = G.add_var("u(t)"),
        v = G.add_var("v(t)"), lam = G.add_var("lam(t)"), g = G.add_var("g");
  M.add_domain( t, FFDom( 0., T, 2, FFDom::LGR, 3 ) );  M.set_evolution_domain( t );
  for( auto s : { x, y, u, v, lam } ) M.add_state( s, {t} );
  M.add_input( g );
  M.update_ref( x, .6 );  M.update_ref( y, .8 );  M.update_ref( u, 0. );  M.update_ref( v, 0. );  M.update_ref( lam, 1. );
  M.add_equation( OpP( x, t ) - u,              {t}, {TI}, INTER );
  M.add_equation( OpP( y, t ) - v,              {t}, {TI}, INTER );
  M.add_equation( OpP( u, t ) + x*lam,          {t}, {TI}, INTER );
  M.add_equation( OpP( v, t ) + y*lam + g,      {t}, {TI}, INTER );
  M.add_equation( x*x + y*y - 1.,               {t}, {FFDom::ALL}, INTER );
  M.add_equation( x - 0.6,                      {t}, {FFDom::LB}, INIT );
  M.add_equation( u - 0.,                       {t}, {FFDom::LB}, INIT );
  // integral outputs: a POINT output on an algebraic state is refused by the adjoint seed by design
  // (lam is algebraic in both reductions, y in the dummy one), and the substitution-reduced index-3
  // pendulum's terminal seed is a separate, pre-existing limitation -- integrals bypass both
  FFIntegral OpI;
  M.add_output( OpI( std::vector<FFVar>{ x }, {t} )[0] );  M.add_output( OpI( std::vector<FFVar>{ lam }, {t} )[0] );
  M.options.SOLVE.RTOL = M.options.SOLVE.ATOL = 1e-10;
  M.options.DISPLAY_LEVEL = 0;  M.options.SOLVE.VERBOSE = 0;  M.options.OUTPUT.RESRECORD = 20;
  M.options.REDUCE.DUMMY_DERIV = dummy;
  if( !M.setup() ){ std::cout << "    setup failed: " << M.extract_error() << "\n"; return R; }
  R.nx = M.nx();
  std::vector<double> P{ 9.81 };
  if( M.solve( P ) != DAESLV_IDAS::NORMAL ){ std::cout << "    solve failed: " << M.extract_error() << "\n"; return R; }
  R.f = M.val_function();
  for( auto const& r : M.results_solve ) R.worst = std::max( R.worst, std::fabs( r.x[0]*r.x[0] + r.x[1]*r.x[1] - 1. ) );
  if( SENS_PENDULUM || !dummy ){
    if( M.solve_fsens( P ) != DAESLV_IDAS::NORMAL ){ std::cout << "    fsens failed: " << M.extract_error() << "\n"; return R; }
    for( auto const& row : M.val_function_gradient() ) for( double v : row ) R.gf.push_back( v );
  }
  if( ASENS_PENDULUM || !dummy ){      // the substitution-reduced pendulum's adjoint runs since the seed was generalised
    if( M.solve_asens( P ) != DAESLV_IDAS::NORMAL ){ std::cout << "    asens failed: " << M.extract_error() << "\n"; return R; }
    for( auto const& row : M.val_function_gradient() ) for( double v : row ) R.ga.push_back( v );
  }
  R.ok = true;  return R;
}

// ---------------------------------------------------------------- C: the reactor with a JUMP, dummy-reduced
// nA(0.5+) = 0.8 nA(0.5-) on the index-2 reactor.  Under dummy derivatives the stop system applies the
// ordinary index-1 rule: the jumped state nA keeps its continuity row replaced by the jump, and the
// algebraic nB RE-SOLVES from nA + nB = 1 at tau+ (substitution leaves the constraint at 0.947 there, with
// a warning).  The point output nB(1) reads an ALGEBRAIC state, which the adjoint seed PROJECTS through
// the full residual at T rather than refusing (2026-10-25).  Forward and adjoint must agree.
static void reactor_jump()
{
  FFPartial OpP;  FFEval OpE;
  int const TI = FFDom::ALL - FFDom::LB;
  auto const INTER = FFModel::EqnRole::INTERIOR, INIT = FFModel::EqnRole::INITIAL;
  FFGraph G;  DAESLV_IDAS M;  M.set_dag( &G );
  FFVar t = G.add_var("t"), nA = G.add_var("nA(t)"), nB = G.add_var("nB(t)"), q = G.add_var("q(t)"),
        k = G.add_var("k"), n0 = G.add_var("nA0");
  M.add_domain( t, FFDom( 0., 1., 2, FFDom::LGR, 3 ) );  M.set_evolution_domain( t );
  M.add_state( nA, {t} );  M.add_state( nB, {t} );  M.add_state( q, {t} );  M.add_input( k );  M.add_input( n0 );
  M.update_ref( nA, .5 );  M.update_ref( nB, .5 );  M.update_ref( q, .1 );
  M.add_equation( OpP( nA, t ) + k*nA + q*nA,    {t}, {TI}, INTER );
  M.add_equation( OpP( nB, t ) - 2.*k*nA + q*nB, {t}, {TI}, INTER );
  M.add_equation( nA + nB - 1.,                  {t}, {FFDom::ALL}, INTER );
  M.add_equation( nA - n0,                       {t}, {FFDom::LB}, INIT );
  M.add_transition( 0.8*nA, nA, t, 0.5 );
  M.add_output( OpE( nA, t, 1. ) );  M.add_output( OpE( nB, t, 1. ) );
  M.options.SOLVE.RTOL = M.options.SOLVE.ATOL = 1e-10;
  M.options.DISPLAY_LEVEL = 0;  M.options.SOLVE.VERBOSE = 0;  M.options.OUTPUT.RESRECORD = 4;
  M.options.REDUCE.DUMMY_DERIV = true;
  if( !M.setup() ){ std::cout << "    setup failed: " << M.extract_error() << "\n"; n_fail += 6; return; }
  check( "jump: states minted (nB is the dummy)", (double)M.nx(), 4., 0. );
  std::vector<double> P{ 1.5, 0.8 };
  if( M.solve( P ) != DAESLV_IDAS::NORMAL ){ std::cout << "    solve failed: " << M.extract_error() << "\n"; n_fail += 5; return; }
  double sum_after = 0.;  bool seen = false;
  for( auto const& r : M.results_solve ) if( std::fabs( r.t - 0.5 ) < 1e-12 ){ sum_after = r.x[0] + r.x[1];  seen = true; }
  check( "jump: nA + nB = 1 holds at tau+ (nB re-solved)", seen? sum_after: -1., 1., 1e-9 );
  check( "jump: nA + nB = 1 at the end", M.val_function()[0] + M.val_function()[1], 1., 1e-9 );
  std::vector<double> gf, ga;
  if( M.solve_fsens( P ) != DAESLV_IDAS::NORMAL ){ std::cout << "    fsens failed: " << M.extract_error() << "\n"; n_fail += 3; return; }
  for( auto const& row : M.val_function_gradient() ) for( double v : row ) gf.push_back( v );
  if( M.solve_asens( P ) != DAESLV_IDAS::NORMAL ){ std::cout << "    asens failed: " << M.extract_error() << "\n"; n_fail += 3; return; }
  for( auto const& row : M.val_function_gradient() ) for( double v : row ) ga.push_back( v );
  for( size_t i = 0; i < gf.size(); ++i )
    check( "jump: asens == fsens " + std::to_string(i) + ( i >= 2? " (nB(1): algebraic, projected seed)": "" ), ga[i], gf[i], 1e-7 );
}

static void compare( std::string const& tag, Result const& S, Result const& D, size_t const nx_dummy, double tol_f, double tol_g, bool const nonlinear )
{
  if( !S.ok || !D.ok ){ std::cout << "    " << tag << ": a run failed\n"; n_fail += 3; return; }
  std::cout << "    " << tag << " gradients  subst fsens:"; for( double g : S.gf ) std::cout << " " << std::setprecision(6) << g;
  std::cout << "  subst asens:"; for( double g : S.ga ) std::cout << " " << g;
  std::cout << "  dummy fsens:"; for( double g : D.gf ) std::cout << " " << g;
  std::cout << "  dummy asens:"; for( double g : D.ga ) std::cout << " " << g; std::cout << "\n";
  check( tag + ": states minted", (double)D.nx, (double)nx_dummy, 0. );
  for( size_t i = 0; i < S.f.size(); ++i )
    check( tag + ": output " + std::to_string(i), D.f[i], S.f[i], tol_f );
  for( size_t i = 0; i < S.gf.size() && i < D.gf.size(); ++i )
    check( tag + ": fsens " + std::to_string(i), D.gf[i], S.gf[i], tol_g );
  for( size_t i = 0; i < S.ga.size() && i < D.ga.size(); ++i )
    check( tag + ": asens " + std::to_string(i), D.ga[i], S.ga[i], tol_g );
  for( size_t i = 0; i < D.gf.size() && i < D.ga.size(); ++i )
    check( tag + ": dummy fsens == asens " + std::to_string(i), D.ga[i], D.gf[i], tol_g );
  if( D.gf.empty() ) std::cout << "    " << tag << ": fsens on the dummy model NOT exercised (see SENS_PENDULUM)\n";
  if( D.ga.empty() ) std::cout << "    " << tag << ": asens on the dummy model NOT exercised (see ASENS_PENDULUM)\n";
  std::cout << "    " << std::left << std::setw( 48 ) << ( tag + ": constraint residual on the path" )
            << std::scientific << std::setprecision(2) << " subst=" << S.worst << "  dummy=" << D.worst << "\n";
  if( nonlinear ) check_lt( tag + ": dummy holds the constraint TIGHTER", D.worst, S.worst );
  else{           // a LINEAR constraint is conserved exactly by BDF under either reduction: both are round-off
    check_lt( tag + ": linear constraint at round-off (subst)", S.worst, 1e-12 );
    check_lt( tag + ": linear constraint at round-off (dummy)", D.worst, 1e-12 );
  }
}

int main()
{
  std::cout << "\n--- B. pendulum, INDEX 3: the cascade -- two dummies on the chain y' -> dy, v' -> dv\n"
            << "    horizon 0.4, before the bob crosses y = 0 (t = 0.51) where the static pivot 2y is singular\n";
  compare( "pendulum", pendulum( false, 0.4 ), pendulum( true, 0.4 ), 7, 1e-8, 1e-5, /*nonlinear=*/ true );

  std::cout << "\n--- C. the reactor with a JUMP nA(0.5+) = 0.8 nA(0.5-), dummy-reduced: the index-1 stop rule, and\n"
            << "    the adjoint of a point output on the ALGEBRAIC nB(1) through the projected seed\n";
  reactor_jump();

  std::cout << "\n==========  DAESLV_dummy: " << n_pass << " passed, " << n_fail << " failed -- "
            << ( n_fail? "FAILURES": "ALL PASS" ) << "  ==========\n\n";
  return n_fail? 1: 0;
}
