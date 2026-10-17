// DAESLV_test2 -- high index: what the integrator does with a reduced model, and where the automatic
// route currently stops.
// =====================================================================================================
// The planar pendulum in Cartesian coordinates, g = L = m = 1, released from rest at (1,0):
//     x' = u ,  y' = v ,  u' = -x lam ,  v' = -y lam - g ,  0 = x^2 + y^2 - 1
// The constraint is INDEX 3: it carries no lam at all, so d(alg)/d(lam) = 0 and no sequential
// integrator can start.  Differentiating it twice and using x^2+y^2 = 1 gives the index-1 form
//     0 = u^2 + v^2 - lam - g y .
//
//  A  HAND-reduced index-1 form.  Integrates, and the three invariants the reduced dynamics do NOT
//     impose -- the position constraint, the energy, and the velocity-level hidden constraint -- all
//     hold, which is what makes this a test of the integrator rather than of the algebra.
//  B  The index-3 form, left to FFModel's Pantelides + dummy-derivative reduction.  The reduction
//     WORKS: it plans the double differentiation and the reduced row reaches the extractor.  What stops
//     it is the INITIAL rows: a materialised hidden constraint is NONLINEAR in the state it pins, so
//     _split_in's affine factorisation gives ic[y] = (1 + y^2 - x^2)/(2y) -- self-referential, and no
//     evaluation order repairs it.  B therefore PINS the current refusal, and is the case to flip when
//     the shared IC solve lands.
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
{
  ok? ++n_pass: ++n_fail;
  std::cout << "    " << std::left << std::setw( 46 ) << w << std::right
            << ( ok? "                                        PASS\n": "                                        FAIL\n" );
}

//! Declare the pendulum.  @p index1 picks the twice-differentiated constraint; otherwise the raw one.
//! @p full_ic gives every state an initial row (possible only for the index-1 form, where the ICs are
//! independent); the index-3 form leaves the hidden constraints to supply the rest.
template <typename T>
static void build( T& M, FFGraph& G, bool const index1, std::vector<FFVar>* out=nullptr )
{
  FFVar t=G.add_var("t"), x=G.add_var("x(t)"), y=G.add_var("y(t)"),
        u=G.add_var("u(t)"), v=G.add_var("v(t)"), lam=G.add_var("lam(t)");
  if( out ) *out = { t, x, y, u, v, lam };
  FFPartial OpP; int const TI = FFDom::ALL - FFDom::LB;
  M.add_domain( t, FFDom( 0., 1., 4, FFDom::LGR, 4 ) );
  for( auto const& s : { x, y, u, v, lam } ) M.add_state( s, {t} );
  M.set_evolution_domain( t );
  M.update_ref( x, 1. ); M.update_ref( y, 0. ); M.update_ref( u, 0. );
  M.update_ref( v, 0. ); M.update_ref( lam, 0. );
  auto io=FFModel::EqnOptions(FFModel::EqnRole::INTERIOR,0), ii=FFModel::EqnOptions(FFModel::EqnRole::INITIAL,0);
  M.add_equation( OpP(x,t) - u,          {t}, {TI}, io );
  M.add_equation( OpP(y,t) - v,          {t}, {TI}, io );
  M.add_equation( OpP(u,t) + x*lam,      {t}, {TI}, io );
  M.add_equation( OpP(v,t) + y*lam + 1., {t}, {TI}, io );
  if( index1 ){
    M.add_equation( u*u + v*v - lam - y, {t}, {FFDom::ALL}, io );      // d^2/dt^2 of the constraint
    M.add_equation( x - 1., {t}, {FFDom::LB}, ii );
    M.add_equation( y - 0., {t}, {FFDom::LB}, ii );
    M.add_equation( u - 0., {t}, {FFDom::LB}, ii );
    M.add_equation( v - 0., {t}, {FFDom::LB}, ii );
  }
  else{
    M.add_equation( x*x + y*y - 1., {t}, {FFDom::ALL}, io );           // INDEX 3
    M.add_equation( x - 1., {t}, {FFDom::LB}, ii );                    // free data {x, u}
    M.add_equation( u - 0., {t}, {FFDom::LB}, ii );
  }
}

int main()
{
  // ---------------------------------------------------------------- A
  { std::cout << "\n--- A. HAND-reduced index-1 pendulum: 0 = u^2 + v^2 - lam - g y\n";
    FFGraph G; DAESLV_IDAS M; M.set_dag( &G );
    build( M, G, true );
    M.options.DISPLAY_LEVEL = 0;
    M.options.SOLVE.RTOL = M.options.SOLVE.ATOL = 1e-10;
    if( !M.setup() ){ std::cout << "    setup failed: " << M.extract_error() << "\n"; ++n_fail; }
    else if( M.solve( std::vector<double>{} ) != DAESLV_IDAS::NORMAL ){
      std::cout << "    solve failed: " << M.extract_error() << "\n"; ++n_fail; }
    else{
      auto const& X = M.val_state().back();
      std::cout << "    differential=" << M.n_differential() << " algebraic=" << M.n_algebraic()
                << "   x=" << std::fixed << std::setprecision(6) << X[0] << " y=" << X[1]
                << " lam=" << X[4] << "\n";
      // none of these is imposed on the reduced dynamics; they hold only if the integration is right
      check( "position constraint x^2 + y^2 = 1", X[0]*X[0] + X[1]*X[1], 1., 1e-7 );
      check( "energy (u^2+v^2)/2 + g y = 0", 0.5*(X[2]*X[2]+X[3]*X[3]) + X[1], 0., 1e-7 );
      check( "hidden constraint x u + y v = 0", X[0]*X[2] + X[1]*X[3], 0., 1e-7 );
      check_bool( "the bob stays below the pivot", X[1] <= 1e-9 );
    } }

  // ---------------------------------------------------------------- B
  { std::cout << "\n--- B. INDEX-3 form, left to the automatic reduction\n";
    FFGraph G; DAESLV_IDAS M; M.set_dag( &G );
    std::vector<FFVar> V;
    build( M, G, false, &V );
    M.options.DISPLAY_LEVEL = 0;
    bool const ok = M.FFModel::setup();
    auto const& F = M.dynamic_form();
    auto const& R = M.reduction_plan();
    std::cout << "    model setup=" << ok << "   index=" << F.declared_index
              << " resolved=" << F.resolved << "\n"
              << "    reduction_plan: " << R.assigns.size() << " assign(s), max_index=" << R.max_index
              << ", resolved=" << R.resolved << "\n";
    for( auto const& a : R.assigns )
      std::cout << "       differentiate the constraint " << a.n_diff << "x to pin "
                << a.pinned_var.name() << "\n";
    // the reduction itself is sound, and that is worth asserting
    check_bool( "Pantelides finds the index-3 block", R.max_index == 3 );
    check_bool( "it plans the DOUBLE differentiation", !R.assigns.empty() && R.assigns.front().n_diff == 2 );
    check_bool( "and resolves it", R.resolved );
    // ... but the materialised hidden constraints are nonlinear in the states they pin, so the INITIAL
    // rows cannot be factorised affinely.  PINNED: flip this when the shared IC solve lands.
    bool const refused = !M.setup();
    std::cout << "    DAESLV setup refused=" << refused << "\n      " << M.extract_error() << "\n";
    check_bool( "refused, naming the implicit INITIAL rows",
                refused && M.extract_error().find( "implicit in each other" ) != std::string::npos );
    // The UNDIVIDED rows are what an IC solve will work from.  They must vanish at the known consistent
    // initial point (1,0,0,0,0) -- the divided form cannot, since ic[y] = (1+y^2-x^2)/(2y) is singular there.
    auto const& rows = M.initial_rows();
    std::cout << "    initial rows kept undivided: " << rows.size()
              << "  (2 declared + the materialised hidden constraints)\n";
    check_bool( "more rows than the 2 declared", rows.size() > 2 );
    // The four rows are  x - 1 ,  u ,  x^2 + y^2 - 1  and  2x u + 2y v : the two declared ones plus the
    // materialised POSITION and VELOCITY constraints, all plain algebraic residuals in the states, with no
    // derivative operators -- the reduction has already substituted x' -> u and y' -> v.  So the IC solve
    // here is an ordinary 4x4 nonlinear system, which is what makes KINSOL over these rows the right shape.
    // They live in the WORKING model's dag, not the one the driver declared into; evaluating them against
    // the wrong graph is what threw while this was being characterised.
  }

  // ---------------------------------------------------------------- C
  { std::cout << "\n--- C. a MODEL-layer option set on the solver must reach the model\n";
    // CLASSIFY.MODE gates _build_reduction_plan/_reduce_high_index, so the reduction plan is an OBSERVABLE
    // consequence of a model option.  Once the options are restructured, the solver's object SHADOWS
    // FFModel's and only _apply_options() slices it down (ocfeslv.hpp:1383) -- a solver that forgets that
    // override silently ignores every REDUCE/CLASSIFY setting, and passes every other check in this file.
    FFGraph G; DAESLV_IDAS M; M.set_dag( &G );
    std::vector<FFVar> V;
    build( M, G, false, &V );
    M.options.DISPLAY_LEVEL = 0;
    M.options.CLASSIFY.MODE = FFModel::Options::CLASS_NONE;   // reduction must NOT run
    // EXPECTED on stderr here: "DOF: the model is NOT balanced -- rows - unknowns = -2".  That is the
    // point of the case: with the reduction suppressed the two hidden constraints are never materialised,
    // so the model IS short two rows.  Warnings are never gated by DISPLAY_LEVEL, by design, so it cannot
    // be silenced -- and should not be.
    M.FFModel::setup();
    std::cout << "    with CLASS_NONE: reduction_plan has " << M.reduction_plan().assigns.size()
              << " assign(s)\n";
    check_bool( "CLASS_NONE reaches the model: no reduction planned", M.reduction_plan().empty() ); }

  std::cout << "\n==========  DAESLV_test2: " << n_pass << " passed, " << n_fail << " failed -- "
            << ( n_fail? "FAILURES": "ALL PASS" ) << "  ==========\n\n";
  return n_fail? 1: 0;
}
