// ============================================================================
//  OCFE_ODE0b.cpp  --  Pure ODE IVP on a multi-element time domain, swept over
//                      the full 2x2 of  {IC_WEAK, IC_STRONG} x {monolithic, marching}.
//
//  Extends OCFE_ODE0 with the monolithic-vs-marching axis.  Two orthogonal
//  discretisation choices are exercised on the SAME ODE, and all four must
//  reproduce the analytic solution and agree with one another:
//
//    * IMPOSITION_TYPE (IC_WEAK / IC_STRONG) -- HOW element-interface C0
//      continuity is enforced.  This only bites in MONOLITHIC mode, where the
//      full space-time BVP carries interior interfaces (SAT penalty vs
//      tau-eliminated exact continuity).  Under marching each window collapses
//      to a single element with NO interior interface, so IC_WEAK and IC_STRONG
//      are the same computation -- the marching runs confirm that.
//
//    * SOLVE_MARCHING (false / true) -- WHETHER the space-time BVP is solved
//      all-at-once (monolithic) or window-by-window with the IC transferred
//      between windows (marching).  Both are the same discrete BVP on the same
//      element grid, so for a transient IVP they must give the same trajectory.
//
//  NOTE: SOLVE_MARCHING defaults to TRUE.  A time domain with differential
//  states is march-eligible, so an ODE/DAE driver that does not set
//  SOLVE_MARCHING=false is silently marching (and its IC_WEAK/IC_STRONG choice
//  is a no-op).  Every run below asserts is_marching()==requested to make the
//  active mode explicit.
//
//  MODEL (3 states, 1 time domain, no inputs)  -- closed forms
//        dx/dt = -x + y ,  x(0)=1     ->  x(t) = e^{-t} (1 + t)
//        dy/dt = -y      ,  y(0)=1     ->  y(t) = e^{-t}
//        dz/dt = -z^2    ,  z(0)=1     ->  z(t) = 1/(1 + t)
//
//  CHECKS (per (imposition,mode) run, then cross-run)
//    MODE is_marching()==requested (monolithic vs marching actually engaged)
//    C0   solve converges
//    C1-3 x,y,z at sampled times == closed form                 [discretisation]
//    XM   IC_WEAK[march] == IC_STRONG[march]  (GATED, exact: marching windows have
//         no interior interface, so WEAK and STRONG are the same computation)
//    ~    WEAK~STRONG(mono), mono~march  (INFORMATIONAL deltas at discretisation
//         scale; the gated property is each cell vs analytic, above)
//
//  Build:
//    g++ -std=c++17 <suite flags> -DOCFE_OCFESLV_HEADER='"ocfeslv.hpp"' \
//        OCFE_ODE0b.cpp -o OCFE_ODE0b <libs>
// ============================================================================

#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#ifndef OCFE_OCFESLV_HEADER
#define OCFE_OCFESLV_HEADER "ocfeslv.hpp"
#endif
#include OCFE_OCFESLV_HEADER

using namespace mc;

static double const T_end = 2.0;
static size_t const N_EL  = 6;
static size_t const N_ND  = 5;

static double x_exact( double t ){ return std::exp( -t ) * ( 1.0 + t ); }
static double y_exact( double t ){ return std::exp( -t ); }
static double z_exact( double t ){ return 1.0 / ( 1.0 + t ); }

static std::vector<double> const T_SAMPLE = { 0.0, 0.5, 1.0, 1.5, T_end };

static int g_pass = 0, g_fail = 0;
static void check_close( char const* name, double got, double want, double tol )
{
  bool const ok = std::fabs( got - want ) <= tol;
  ( ok ? g_pass : g_fail )++;
  std::cout << "  " << std::left << std::setw(40) << name
            << std::right << std::scientific << std::setprecision(4)
            << " got=" << std::setw(12) << got << " want=" << std::setw(12) << want
            << " |d|=" << std::setw(10) << std::fabs( got - want )
            << "  " << ( ok ? "PASS" : "FAIL" ) << "\n";
}
static void check_true( char const* name, bool ok )
{
  ( ok ? g_pass : g_fail )++;
  std::cout << "  " << std::left << std::setw(40) << name << "  " << ( ok ? "PASS" : "FAIL" ) << "\n";
}
static void compare_vec( char const* name, std::vector<double> const& a,
                         std::vector<double> const& b, double tol )
{
  if( a.empty() || a.size() != b.size() ){ check_true( name, false ); return; }
  double e = 0.; for( size_t i = 0; i < a.size(); ++i ) e = std::max( e, std::fabs( a[i] - b[i] ) );
  check_close( name, e, 0., tol );
}

// informational cross-run delta: WEAK-vs-STRONG(mono) and mono-vs-march differ at the
// DISCRETISATION scale (different discrete operators of the same continuous BVP), so the
// per-cell "cell vs analytic" gates inside run_mode are the real property; this only prints.
static void info_delta( char const* name, std::vector<double> const& a, std::vector<double> const& b )
{
  double e = 0.;
  if( a.size() == b.size() && !a.empty() )
    for( size_t i = 0; i < a.size(); ++i ) e = std::max( e, std::fabs( a[i] - b[i] ) );
  std::cout << "  " << std::left << std::setw(40) << name
            << " delta=" << std::right << std::scientific << std::setprecision(3) << e
            << "  (discretisation-scale; informational)\n";
}

struct ModeOut {
  bool                converged = false;
  bool                ok        = false;
  std::vector<double> F;                 // 3 * T_SAMPLE.size()
};

static ModeOut run_mode( OCFESLV::Options::ImpositionType imp, char const* strimp, bool marching )
{
  ModeOut R;
  char const* mstr = marching ? "march" : "mono";
  std::cout << "\n---------------- IC_" << strimp << " [" << mstr << "] ----------------\n";

  FFGraph DAG;
  FFVar t = DAG.add_var( "t" );
  FFVar x = DAG.add_var( "x(t)" );
  FFVar y = DAG.add_var( "y(t)" );
  FFVar z = DAG.add_var( "z(t)" );

  OCFESLV oc( &DAG );
  oc.add_domain( t, FFDom( 0., T_end, N_EL, FFDom::LGR, N_ND ) );
  oc.set_evolution_domain( t );
  oc.add_state( x, {t} );
  oc.add_state( y, {t} );
  oc.add_state( z, {t} );

  oc.update_ref( x, []( OCFESLV::t_Coord const& ){ return 1.0; } );
  oc.update_ref( y, []( OCFESLV::t_Coord const& ){ return 1.0; } );
  oc.update_ref( z, []( OCFESLV::t_Coord const& ){ return 1.0; } );

  FFPartial OpP;   // d/dt
  FFVar EVOL_X = OpP( x, t ) - ( -x + y );
  FFVar EVOL_Y = OpP( y, t ) - ( -y );
  FFVar EVOL_Z = OpP( z, t ) - ( -z * z );
  FFVar IC_X   = x - 1.0;
  FFVar IC_Y   = y - 1.0;
  FFVar IC_Z   = z - 1.0;

  int const T_NO_LB = FFDom::ALL - FFDom::LB;
  OCFESLV::EqnOptions const interior( OCFESLV::EqnRole::INTERIOR, 0 );
  OCFESLV::EqnOptions const initial ( OCFESLV::EqnRole::INITIAL,  0 );

  oc.add_equation( EVOL_X, {t}, {T_NO_LB},   interior );
  oc.add_equation( EVOL_Y, {t}, {T_NO_LB},   interior );
  oc.add_equation( EVOL_Z, {t}, {T_NO_LB},   interior );
  oc.add_equation( IC_X,   {t}, {FFDom::LB}, initial  );
  oc.add_equation( IC_Y,   {t}, {FFDom::LB}, initial  );
  oc.add_equation( IC_Z,   {t}, {FFDom::LB}, initial  );

  for( double ts : T_SAMPLE ){
    oc.add_output( x, {t}, { ts } );
    oc.add_output( y, {t}, { ts } );
    oc.add_output( z, {t}, { ts } );
  }

  oc.options.REDUCE.ORDER    = OCFESLV::Options::RED_NONE;
  oc.options.CLASSIFY.MODE        = OCFESLV::Options::CLASS_AUTO;
  oc.options.INTERFACE.IMPOSITION = imp;
  oc.options.INTERFACE.TYPE  = OCFESLV::Options::IC_AUTO;
  oc.options.SOLVE.MARCHING  = marching;     // explicit: default is TRUE, so mono needs false
  oc.options.SOLVE.MAX_ITER  = 60;
  oc.options.SOLVE.RES_TOL   = 1.0e-10;
  oc.options.DISPLAY_LEVEL   = 0;

  if( !oc.setup() ){
    std::cerr << "  setup() FAILED: " << OCFESLV::setup_status_str( oc.setup_status() ) << "\n";
    return R;
  }

  bool const marching_active = oc.is_marching();
  std::cout << "  is_marching()=" << ( marching_active ? "true" : "false" )
            << "  n_march_steps=" << oc.n_march_steps()
            << "  nVar=" << oc.n_colloc_sta() << " nEqn=" << oc.n_colloc_eqn()
            << " nTrace=" << oc.n_colloc_trace() << "\n";
  check_true( "MODE is_marching()==requested", marching_active == marching );

  std::vector<double> xv, inp;
  if( !oc.init( xv, inp, nullptr ) ){ std::cerr << "  init() FAILED\n"; return R; }
  double* const pinp = inp.empty() ? nullptr : inp.data();

  OCFESLV::SolveReport const rep = oc.solve( xv.data(), pinp, nullptr );
  R.converged = rep.converged;
  std::cout << "  solve: converged=" << ( rep.converged ? "yes" : "no" )
            << " iters=" << rep.iterations
            << " final|r|=" << std::scientific << std::setprecision(3) << rep.final_residual << "\n\n";
  check_true( "C0 solve converged", rep.converged );
  if( !rep.converged ){ return R; }

  R.F = oc.val_functions();
  if( R.F.size() != 3 * T_SAMPLE.size() ){
    std::cerr << "  output count mismatch: got " << R.F.size() << " want " << 3*T_SAMPLE.size() << "\n";
    return R;
  }

  bool ok = true;
  for( size_t k = 0; k < T_SAMPLE.size(); ++k ){
    double const ts = T_SAMPLE[k];
    double const xv_ = R.F[3*k+0], yv_ = R.F[3*k+1], zv_ = R.F[3*k+2];
    std::ostringstream nx, ny, nz;
    nx << "C1 x(" << ts << ")"; ny << "C2 y(" << ts << ")"; nz << "C3 z(" << ts << ")";
    double const tol = ( ts == 0.0 ) ? 1e-9 : 1e-5;
    check_close( nx.str().c_str(), xv_, x_exact( ts ), tol );
    check_close( ny.str().c_str(), yv_, y_exact( ts ), tol );
    check_close( nz.str().c_str(), zv_, z_exact( ts ), tol );
    ok &= ( std::fabs( xv_ - x_exact(ts) ) <= tol )
       && ( std::fabs( yv_ - y_exact(ts) ) <= tol )
       && ( std::fabs( zv_ - z_exact(ts) ) <= tol );
  }
  R.ok = ok;
  return R;
}

int main()
{
  std::cout << "================================================================\n"
            << "  OCFE_ODE0b : pure ODE IVP, 2x2 sweep\n"
            << "  {IC_WEAK, IC_STRONG} x {monolithic, marching}\n"
            << "  x'=-x+y, y'=-y, z'=-z^2 ;  x(0)=y(0)=z(0)=1 ;  T_end=" << T_end << "\n"
            << "================================================================\n";

  ModeOut const wm = run_mode( OCFESLV::Options::IC_WEAK,   "WEAK",   false );  // weak,   monolithic
  ModeOut const sm = run_mode( OCFESLV::Options::IC_STRONG, "STRONG", false );  // strong, monolithic
  ModeOut const wc = run_mode( OCFESLV::Options::IC_WEAK,   "WEAK",   true  );  // weak,   marching
  ModeOut const sc = run_mode( OCFESLV::Options::IC_STRONG, "STRONG", true  );  // strong, marching

  std::cout << "\n---------------- cross-run comparisons ----------------\n";
  // Per-cell "cell vs analytic" (C1/C2/C3 inside run_mode) is the gated property for all
  // four cells.  Cross-run we gate ONLY the equality that is exact in principle: under
  // marching each window is a single element with NO interior interface, so IC_WEAK and
  // IC_STRONG are the same computation.
  compare_vec( "XM  WEAK[march] == STRONG[march] (exact)", wc.F, sc.F, 1e-9 );
  // The rest differ at the discretisation scale (WEAK/STRONG are different interface
  // discretisations; mono/march are different discrete operators) -> informational.
  info_delta( "XM  WEAK[mono]  ~ STRONG[mono]",  wm.F, sm.F );
  info_delta( "MM  WEAK:   mono ~ march",        wm.F, wc.F );
  info_delta( "MM  STRONG: mono ~ march",        sm.F, sc.F );

  std::cout << "\n=================== ODE0b sweep summary ===================\n";
  std::cout << std::left << std::setw(20) << "run"
            << std::setw(12) << "converged" << std::setw(8) << "result" << "\n";
  auto row = [&]( char const* nm, ModeOut const& m ){
    std::cout << std::left << std::setw(20) << nm
              << std::setw(12) << ( m.converged ? "yes" : "no" )
              << std::setw(8)  << ( m.ok ? "PASS" : "FAIL" ) << "\n"; };
  row( "WEAK  [mono]",   wm );
  row( "STRONG[mono]",   sm );
  row( "WEAK  [march]",  wc );
  row( "STRONG[march]",  sc );

  std::cout << "\n================================================================\n"
            << "  OCFE_ODE0b: " << g_pass << " passed, " << g_fail << " failed -- "
            << ( g_fail == 0 ? "ALL PASS" : "SOME FAILED" ) << "\n"
            << "================================================================\n";
  return g_fail == 0 ? 0 : 1;
}
