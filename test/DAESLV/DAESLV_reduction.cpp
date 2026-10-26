// DAESLV_reduction -- SUBSTITUTION against DUMMY DERIVATIVES on the index-2 vented reactor, side by side.
// =====================================================================================================
// CRONOS reduces a high-index block one of two ways.  SUBSTITUTION (the default) differentiates the
// consumed constraint and overwrites its row with the derivative, keeping the constraint itself only at
// the initial point.  DUMMY DERIVATIVES (REDUCE.DUMMY_DERIV) keep the constraint as an interior row,
// add its derivative as another, and mint one algebraic state to stand for the derivative of a chosen
// state.  DAESLV_dummy gates their agreement; this driver TABULATES it -- values, forward and adjoint
// sensitivities, the constraint residual along the path, step and Jacobian counts, wall time -- at two
// tolerances, against the closed form where one exists.
//
// The vented batch reactor of the DAESLV tutorial: A -> 2B at constant pressure, the vent rate q the
// algebraic state, 0 = nA + nB - nT the constraint; index 2, since q is absent from it.  Its constraint
// is LINEAR, and that puts the two methods in an unusual order: SUBSTITUTION holds it to ROUND-OFF,
// because BDF conserves a linear invariant exactly, while DUMMY DERIVATIVES hold it as an algebraic row,
// to Newton's tolerance -- 4e-10 at rtol 1e-6, round-off at 1e-10 -- well inside the integration
// tolerance but not exact.  The drift that dummy derivatives exist to remove shows on NONLINEAR
// constraints (DAESLV_drift, DAESLV_dummy's pendulum: 2e-9 against 3e-13); on this model the two are a
// wash on accuracy.  What the table does separate is COST: the dummy's row is algebraic and the system
// IDAS sees is better conditioned, so it takes fewer steps and Jacobians for one extra state.
//
// Every line is a check: the two reductions must agree to the integration tolerance, both must match
// the closed form, and the residual must be round-off for both.  The numbers print either way.
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include "daeslv_idas.hpp"
using namespace mc;
struct R { bool ok=false; size_t nx=0; double f[2]={0,0}; double gf[4]={0,0,0,0}, ga[4]={0,0,0,0}; double worst=0.; long steps=0, jac=0; double tsolve=0, tfs=0, tas=0; };
static R run( bool dummy, double rtol ){
  FFPartial OpP; FFEval OpE; FFIntegral OpI; int const TI = FFDom::ALL - FFDom::LB;
  auto const INTER = FFModel::EqnRole::INTERIOR, INIT = FFModel::EqnRole::INITIAL;
  R r; FFGraph G; DAESLV_IDAS M; M.set_dag( &G );
  FFVar t = G.add_var("t"), nA = G.add_var("nA(t)"), nB = G.add_var("nB(t)"), q = G.add_var("q(t)"), k = G.add_var("k"), n0 = G.add_var("nA0");
  M.add_domain( t, FFDom( 0., 1., 4, FFDom::LGR, 3 ) ); M.set_evolution_domain( t );
  M.add_state( nA, {t} ); M.add_state( nB, {t} ); M.add_state( q, {t} ); M.add_input( k ); M.add_input( n0 );
  M.update_ref( nA, .5 ); M.update_ref( nB, .5 ); M.update_ref( q, .1 );
  M.add_equation( OpP(nA,t) + k*nA + q*nA, {t}, {TI}, INTER );
  M.add_equation( OpP(nB,t) - 2.*k*nA + q*nB, {t}, {TI}, INTER );
  M.add_equation( nA + nB - 1., {t}, {FFDom::ALL}, INTER );
  M.add_equation( nA - n0, {t}, {FFDom::LB}, INIT );
  M.add_output( OpE( nA, t, 1. ) ); M.add_output( OpI( std::vector<FFVar>{ q }, {t} )[0] );
  M.options.SOLVE.RTOL = M.options.SOLVE.ATOL = rtol; M.options.DISPLAY_LEVEL = 0; M.options.SOLVE.VERBOSE = 0; M.options.OUTPUT.RESRECORD = 50;
  M.options.REDUCE.DUMMY_DERIV = dummy;
  if( !M.setup() ) return r;
  r.nx = M.nx(); std::vector<double> P{ 1.5, 0.8 };
  auto t0 = std::chrono::steady_clock::now();
  if( M.solve( P ) != DAESLV_IDAS::NORMAL ) return r;
  r.tsolve = std::chrono::duration<double,std::milli>( std::chrono::steady_clock::now() - t0 ).count();
  r.f[0] = M.val_function()[0]; r.f[1] = M.val_function()[1]; r.steps = M.stats_solve.numSteps; r.jac = M.stats_solve.numJAC;
  for( auto const& x : M.results_solve ) r.worst = std::max( r.worst, std::fabs( x.x[0] + x.x[1] - 1. ) );
  t0 = std::chrono::steady_clock::now();
  if( M.solve_fsens( P ) != DAESLV_IDAS::NORMAL ) return r;
  r.tfs = std::chrono::duration<double,std::milli>( std::chrono::steady_clock::now() - t0 ).count();
  { auto g = M.val_function_gradient(); r.gf[0]=g[0][0]; r.gf[1]=g[0][1]; r.gf[2]=g[1][0]; r.gf[3]=g[1][1]; }
  t0 = std::chrono::steady_clock::now();
  if( M.solve_asens( P ) != DAESLV_IDAS::NORMAL ) return r;
  r.tas = std::chrono::duration<double,std::milli>( std::chrono::steady_clock::now() - t0 ).count();
  { auto g = M.val_function_gradient(); r.ga[0]=g[0][0]; r.ga[1]=g[0][1]; r.ga[2]=g[1][0]; r.ga[3]=g[1][1]; }
  r.ok = true; return r; }
static int n_pass = 0, n_fail = 0;
static void note( bool ok ){ ok? ++n_pass: ++n_fail; }
int main(){
  // closed form: nA(1) = a e^-k / (1 - a e^-k), a = nA0/(1+nA0), k=1.5, nA0=0.8 (nT=1)
  double const K=1.5, N0=0.8, a = N0/(1.+N0), e = std::exp(-K), nAex = a*e/(1.-a*e);
  double const dnA_dn0 = ( 1./((1.+N0)*(1.+N0)) ) * ( e/((1.-a*e)*(1.-a*e)) );
  for( double rtol : { 1e-6, 1e-10 } ){
    R S = run( false, rtol ), D = run( true, rtol );
    if( !S.ok || !D.ok ){ std::cout << "  a run FAILED at rtol " << rtol << "\n"; n_fail += 10; continue; }
    double const tolSD = 50.*rtol, tolEX = 50.*rtol;   // agreement to the integration tolerance
    double const tolAS = 1000.*rtol;                     // the adjoint's backward integration is an order looser
    std::cout << "\n=== rtol = atol = " << std::scientific << std::setprecision(0) << rtol << "   (closed form: nA(1) = " << std::setprecision(10) << nAex << ",  d nA(1)/d nA0 = " << dnA_dn0 << ")\n";
    auto row = [&]( char const* name, double s, double d, double ex = NAN, bool adjoint = false ){
      double const tS = adjoint? tolAS: tolSD, tE = adjoint? tolAS: tolEX;
      bool ok = std::fabs(s-d) <= tS && ( std::isnan(ex) || ( std::fabs(s-ex) <= tE && std::fabs(d-ex) <= tE ) );
      note( ok );
      std::cout << "  " << std::left << std::setw(26) << name << std::right << std::setprecision(10) << std::setw(18) << s << std::setw(18) << d
                << "   |S-D| " << std::setprecision(1) << std::fabs(s-d);
      if( !std::isnan(ex) ) std::cout << "   |S-exact| " << std::fabs(s-ex) << "  |D-exact| " << std::fabs(d-ex);
      std::cout << ( ok? "  PASS\n": "  FAIL\n" ); };
    std::cout << "  " << std::left << std::setw(26) << "" << std::right << std::setw(18) << "SUBSTITUTION" << std::setw(18) << "DUMMY DERIV" << "\n";
    note( S.nx == 3 && D.nx == 4 );
    std::cout << "  " << std::left << std::setw(26) << "states" << std::right << std::setw(18) << S.nx << std::setw(18) << D.nx << ( S.nx == 3 && D.nx == 4? "  PASS\n": "  FAIL\n" );
    row( "nA(1)", S.f[0], D.f[0], nAex );  row( "vented moles", S.f[1], D.f[1] );
    row( "fsens d nA(1)/dk",   S.gf[0], D.gf[0] ); row( "fsens d nA(1)/d nA0", S.gf[2], D.gf[2], dnA_dn0 );
    row( "fsens d vent/dk",    S.gf[1], D.gf[1] ); row( "fsens d vent/d nA0",  S.gf[3], D.gf[3] );
    row( "asens d nA(1)/dk",   S.ga[0], D.ga[0], NAN, true ); row( "asens d nA(1)/d nA0", S.ga[2], D.ga[2], dnA_dn0, true );
    row( "asens d vent/dk",    S.ga[1], D.ga[1], NAN, true ); row( "asens d vent/d nA0",  S.ga[3], D.ga[3], NAN, true );
    // linear constraint: substitution at ROUND-OFF (BDF conserves it), dummy derivatives at Newton's tolerance
    bool const okres = S.worst < 1e-12 && D.worst < 1e-3*rtol;
    note( okres );
    std::cout << "  " << std::left << std::setw(26) << "max|nA+nB-1| on path" << std::right << std::scientific << std::setprecision(2) << std::setw(18) << S.worst << std::setw(18) << D.worst
              << ( okres? "  PASS (subst round-off; dummy << rtol)\n": "  FAIL\n" );
    std::cout << "  " << std::left << std::setw(26) << "steps / Jacobians" << std::right << std::setw(18) << (std::to_string(S.steps)+" / "+std::to_string(S.jac)) << std::setw(18) << (std::to_string(D.steps)+" / "+std::to_string(D.jac)) << "\n";
    std::cout << "  " << std::left << std::setw(26) << "ms: solve/fsens/asens" << std::right << std::fixed << std::setprecision(1)
              << std::setw(18) << (std::to_string((int)S.tsolve)+"/"+std::to_string((int)S.tfs)+"/"+std::to_string((int)S.tas))
              << std::setw(18) << (std::to_string((int)D.tsolve)+"/"+std::to_string((int)D.tfs)+"/"+std::to_string((int)D.tas)) << "\n";
  }
  std::cout << "\n==========  DAESLV_reduction: " << n_pass << " passed, " << n_fail << " failed -- " << ( n_fail? "FAILURES": "ALL PASS" ) << "  ==========\n\n";
  return n_fail? 1: 0; }
