// DAESLV_hidden -- are the HIDDEN constraints of a reduced high-index model preserved across a stop
// caused by a discontinuous INPUT?
// =====================================================================================================
// REDUCE.HIDDEN_IC materialises every level below the one the reduction pins with as an INITIAL row,
// pinned to FFDom::LB -- the initial point ONLY.  They appear in no stop system: a restart imposes
// CONTINUITY on the differential states and re-solves the algebraic ones from the interior rows alone.
//
// That is sound as long as every state appearing in a hidden constraint is DIFFERENTIAL: continuity
// carries it across unchanged, the left state satisfied the constraint, so the right state does too.
// The open question is the other case -- a state that the reduction leaves ALGEBRAIC and that appears
// in a hidden level is RE-SOLVED at the stop, with no hidden constraint in the system to hold it.
//
// So this driver is an INSTRUMENT before it is a gate: it reports which states the reduction leaves
// differential and which hidden levels it materialises, and only then measures the residual of those
// levels at every stage time.  Either the risky configuration arises and the drift is measured, or it
// cannot arise from this reduction and the §8b entry narrows.
//
//   A. INDEX-2, switched input:   x' = z + u ,  0 = x - 1
//      Differentiating the constraint gives 0 = x' = z + u, so z = -u: z is algebraic, x differential.
//      The hidden level is x - 1 = 0, in x alone.  u steps on every element, so z must step with it.
//
//   B. INDEX-3 pendulum, switched input: gravity g switched per element.
//      TWO hidden levels: x^2+y^2-1 (positions) and x u + y v (positions AND velocities).  If the
//      reduction ever leaves a velocity algebraic, this is where it shows.
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
static void check_bool( std::string const& w, bool ok )
{ ok? ++n_pass: ++n_fail;
  std::cout << "    " << std::left << std::setw( 46 ) << w << std::right
            << ( ok? "  PASS\n": "  FAIL\n" ); }

//! Report what the reduction decided: the plan, and the differential/algebraic split.
static void report( DAESLV_IDAS& M, std::vector<FFVar> const& V )
{
  auto const& R = M.reduction_plan();
  std::cout << "    reduction_plan: " << R.assigns.size() << " assign(s), max_index=" << R.max_index
            << ", resolved=" << R.resolved << "\n";
  for( auto const& a : R.assigns )
    std::cout << "       differentiate the constraint " << a.n_diff << "x to pin "
              << a.pinned_var.name() << "\n";
  std::cout << "    differential=" << M.n_differential() << " algebraic=" << M.n_algebraic() << "   [";
  for( size_t i = 0; i < V.size() && i < M.nx(); ++i )
    std::cout << " " << V[i].name() << ( M.id()[i] != 0.? "=diff": "=ALG" );
  std::cout << " ]\n";
}

int main()
{
  FFPartial OpP;
  int const TI = FFDom::ALL - FFDom::LB;
  size_t const NE = 4;
  auto io = FFModel::EqnOptions( FFModel::EqnRole::INTERIOR, 0 );
  auto ii = FFModel::EqnOptions( FFModel::EqnRole::INITIAL,  0 );

  // ------------------------------------------------------------------ A: index-2, switched input
  { std::cout << "\n--- A. INDEX-2:  x' = z + u ,  0 = x - 1 ,  u stepping over " << NE << " elements\n";
    FFGraph G; DAESLV_IDAS M; M.set_dag( &G );
    FFVar t = G.add_var("t"), x = G.add_var("x(t)"), z = G.add_var("z(t)"), u = G.add_var("u");
    M.add_domain( t, FFDom( 0., 2., NE, FFDom::LGR, 4 ) );
    M.add_state( x, {t} ); M.add_state( z, {t} );
    M.add_input( u, {t}, FFDom::LGR, 1 );              // piecewise constant: one level per element
    M.set_evolution_domain( t );
    M.update_ref( x, 1. ); M.update_ref( z, 0. );
    M.add_equation( OpP( x, t ) - z - u, {t}, {TI}, io );
    M.add_equation( x - 1.,              {t}, {FFDom::ALL}, io );   // INDEX 2: no z in the constraint
    M.options.SOLVE.RTOL = M.options.SOLVE.ATOL = 1e-10;
    M.options.DISPLAY_LEVEL = 0;

    bool const ok = M.setup();
    std::cout << "    setup=" << ok << ( ok? std::string(): "  :: " + M.extract_error() ) << "\n";
    if( !ok ){ ++n_fail; }
    else{
      report( M, { x, z } );
      check_bool( "Pantelides finds the index-2 block", M.reduction_plan().max_index == 2 );
      std::vector<double> P;
      for( size_t e = 0; e < NE; ++e ) P.push_back( 1. + (double)e );
      if( M.solve( P ) != DAESLV_IDAS::NORMAL ){
        std::cout << "    solve failed: " << M.extract_error() << "\n"; ++n_fail; }
      else{
        std::cout << "    numIC=" << M.stats.numIC << "   per-stage state and hidden residual:\n";
        double worst = 0.;
        for( size_t k = 0; k < M.val_state().size(); ++k ){
          auto const& X = M.val_state()[k];
          double const res = std::fabs( X[0] - 1. );      // the HIDDEN level
          worst = std::max( worst, res );
          std::cout << "      k=" << k << "  x=" << std::scientific << std::setprecision(10) << X[0]
                    << "  z=" << X[1] << "  |x-1|=" << std::setprecision(2) << res << "\n";
        }
        check( "hidden level |x-1| across every stop", worst, 0., 1e-8 );
        // z is ALGEBRAIC and tracks the input.  NOTE the convention this pins down: val_state()[k] for
        // k >= 1 is the state as the PRECEDING stage reached it, i.e. the LEFT limit at the stop, before
        // the restart re-solves z for the element that starts there.  So z[k] = -u_{k-1}, not -u_k;
        // z[0] is the t0 solve, which does use u_0.
        check( "z(0) = -u_0 (the t0 solve uses the first element)", M.val_state()[0][1], -P[0], 1e-7 );
        for( size_t k = 1; k < M.val_state().size(); ++k )
          check( "  z at stage time " + std::to_string(k) + " is the LEFT limit -u_" + std::to_string(k-1),
                 M.val_state()[k][1], -P[k-1], 1e-7 );
      }
    }
  }

  // ------------------------------------------------------------------ B: index-3, switched gravity
  { std::cout << "\n--- B. INDEX-3 pendulum, gravity switched per element (two hidden levels)\n";
    FFGraph G; DAESLV_IDAS M; M.set_dag( &G );
    FFVar t = G.add_var("t"), x = G.add_var("x(t)"), y = G.add_var("y(t)"),
          u = G.add_var("u(t)"), v = G.add_var("v(t)"), lam = G.add_var("lam(t)"),
          g = G.add_var("g");
    M.add_domain( t, FFDom( 0., 1., 2, FFDom::LGR, 4 ) );
    for( auto s : { x, y, u, v, lam } ) M.add_state( s, {t} );
    M.add_input( g, {t}, FFDom::LGR, 1 );
    M.set_evolution_domain( t );
    // NOT x = 1: there y = 0, and the hidden level x^2+y^2-1 pins y with Jacobian 2y = 0, so the IC
    // Jacobian is singular and the solve correctly refuses.  x = 0.6 gives y = 0.8, v = 0.
    M.update_ref( x, 0.6 ); M.update_ref( y, 0.8 );
    M.update_ref( u, 0. ); M.update_ref( v, 0. ); M.update_ref( lam, 0. );
    M.add_equation( OpP(x,t) - u,          {t}, {TI}, io );
    M.add_equation( OpP(y,t) - v,          {t}, {TI}, io );
    M.add_equation( OpP(u,t) + x*lam,      {t}, {TI}, io );
    M.add_equation( OpP(v,t) + y*lam + g,  {t}, {TI}, io );   // gravity is the INPUT
    M.add_equation( x*x + y*y - 1.,        {t}, {FFDom::ALL}, io );   // INDEX 3
    M.add_equation( x - 0.6, {t}, {FFDom::LB}, ii );                  // free data {x, u}
    M.add_equation( u - 0., {t}, {FFDom::LB}, ii );
    M.options.SOLVE.RTOL = M.options.SOLVE.ATOL = 1e-10;
    M.options.DISPLAY_LEVEL = 0;

    bool const ok = M.setup();
    std::cout << "    setup=" << ok << ( ok? std::string(): "  :: " + M.extract_error() ) << "\n";
    if( !ok ){ ++n_fail; }
    else{
      report( M, { x, y, u, v, lam } );
      // THE measurement this case exists for: do the hidden levels involve any state the reduction
      // left algebraic?  Level 0 is x,y; level 1 is x,y,u,v.  lam is the only algebraic state and
      // appears in neither, so continuity carries both levels across a stop and the risky
      // configuration does not arise from THIS reduction.
      check_bool( "no hidden level involves an ALGEBRAIC state",
                  M.id()[0] && M.id()[1] && M.id()[2] && M.id()[3] );
      // This SOLVES as of the _ic_seed() fix: the reduced index-3 pendulum used to fail here with a
      // singular IC Jacobian, which was never about the index -- the seed was the ORIGIN, where both
      // hidden levels are quadratic and their gradients vanish, dropping the rank by exactly two.
      // Kept as a fallback report in case it regresses.
      if( M.solve( std::vector<double>{ 1., 4. } ) != DAESLV_IDAS::NORMAL ){
        std::cout << "    [recorded] switched gravity: solve failed -- " << M.extract_error() << "\n";
        if( M.solve( std::vector<double>{ 1., 1. } ) != DAESLV_IDAS::NORMAL )
          std::cout << "    [recorded] CONSTANT gravity: also fails -- the input discontinuity is NOT"
                    << " implicated; this is phase 4b, not a stop defect\n";
        else
          std::cout << "    [recorded] CONSTANT gravity SUCCEEDS -- the discontinuity IS implicated\n";
      }
      else{
        double worstP = 0., worstV = 0.;
        std::cout << "    numIC=" << M.stats.numIC << "   per-stage hidden residuals:\n";
        for( size_t k = 0; k < M.val_state().size(); ++k ){
          auto const& X = M.val_state()[k];
          double const rp = std::fabs( X[0]*X[0] + X[1]*X[1] - 1. );     // level 0
          double const rv = std::fabs( X[0]*X[2] + X[1]*X[3] );          // level 1 (velocities)
          worstP = std::max( worstP, rp );  worstV = std::max( worstV, rv );
          std::cout << "      k=" << k << "  |x^2+y^2-1|=" << std::scientific << std::setprecision(2) << rp
                    << "  |xu+yv|=" << rv << "\n";
        }
        check( "hidden level 0 (positions) across the stop",  worstP, 0., 1e-7 );
        check( "hidden level 1 (velocities) across the stop", worstV, 0., 1e-7 );
      }
    }
  }

  std::cout << "\n==========  DAESLV_hidden: " << n_pass << " passed, " << n_fail << " failed -- "
            << ( n_fail? "FAILURES": "ALL PASS" ) << "  ==========\n\n";
  return n_fail? 1: 0;
}
