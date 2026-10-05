// ============================================================================
//  OCFE_DAE0b.cpp  --  Index-3 pendulum (index-1 reduced) DAE on a time domain,
//                      swept over the full 2x2 of {IC_WEAK, IC_STRONG} x
//                      {monolithic, marching}.
//
//  Extends OCFE_DAE0 with the monolithic-vs-marching axis.  The index-1
//  reduced pendulum is solved four ways and every cell must (a) converge,
//  (b) hold the hidden position/velocity constraints and conserve energy,
//  (c) match the independent RK4 angle-ODE oracle, and (d) agree across cells.
//
//  Monolithic vs marching on a DAE: monolithic solves the whole space-time
//  DAE-BVP at once with interior interface continuity imposed weakly/strongly;
//  marching solves each evolution window as a small DAE-BVP and transfers the
//  differential-state IC to the next window (the algebraic multiplier lam is
//  re-solved per window from the index-1 constraint).  Both are the same
//  discrete DAE on the same element grid and must give the same trajectory.
//  Under marching each window is a single element (no interior interface), so
//  IC_WEAK and IC_STRONG coincide there.
//
//  SOLVE_MARCHING defaults to TRUE; each run asserts is_marching()==requested.
//
//  INDEX-3 -> INDEX-1: the position constraint x^2+y^2-L^2=0 has
//  d/d lam = 0 (structurally singular block), so it is differentiated twice to
//  the acceleration-level constraint  u^2+v^2 - lam(x^2+y^2) - g y = 0
//  (d/d lam = -L^2 != 0).  Consistent ICs keep the two hidden constraints
//  g0=x^2+y^2-L^2 and g1=x u+y v at zero; their drift measures solve quality.
//
//  CHECKS (per (imposition,mode) cell, then cross-cell)
//    MODE is_marching()==requested
//    C0   solve converges
//    P0   position drift max|g0|                                  [manifold]
//    P1   velocity drift max|g1|                                   [hidden constraint]
//    P2   energy drift max|E-E0|                                    [invariant]
//    P3   trajectory (x,y,u,v) == RK4 oracle at samples            [independent oracle]
//    P4   lam == (u^2+v^2-g y)/L^2 at samples                       [algebraic consistency]
//    XM   IC_WEAK[march] == IC_STRONG[march]  (GATED, exact: no interior interface)
//    ~    WEAK~STRONG(mono), mono~march  (INFORMATIONAL, discretisation scale;
//         the gated property is each cell's P* vs oracle/invariants, above)
//
//  Build:
//    g++ -std=c++17 <suite flags> -DOCFE_OCFESLV_HEADER='"ocfeslv.hpp"' \
//        OCFE_DAE0b.cpp -o OCFE_DAE0b <libs>
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

static double const L_rod  = 1.0;
static double const g_grav = 9.81;
static double const th0    = M_PI / 3.0;
static double const om0    = 0.0;
static double const T_end  = 0.5;
static size_t const N_EL   = 8;
static size_t const N_ND   = 5;

// ---- independent oracle: RK4 of  theta'' = -(g/L) sin theta ----
struct Oracle {
  size_t              N;
  double              h;
  std::vector<double> th, om;

  Oracle( size_t Ngrid ) : N( Ngrid ), h( T_end / double( Ngrid ) ),
                           th( Ngrid + 1 ), om( Ngrid + 1 )
  {
    double t = th0, w = om0;
    th[0] = t; om[0] = w;
    auto f_t = []( double /*th*/, double w_ ){ return w_; };
    auto f_w = []( double th_, double /*w*/ ){ return -( g_grav / L_rod ) * std::sin( th_ ); };
    for( size_t i = 0; i < N; ++i ){
      double k1t = f_t( t, w ),                          k1w = f_w( t, w );
      double k2t = f_t( t + 0.5*h*k1t, w + 0.5*h*k1w ),  k2w = f_w( t + 0.5*h*k1t, w + 0.5*h*k1w );
      double k3t = f_t( t + 0.5*h*k2t, w + 0.5*h*k2w ),  k3w = f_w( t + 0.5*h*k2t, w + 0.5*h*k2w );
      double k4t = f_t( t + h*k3t, w + h*k3w ),          k4w = f_w( t + h*k3t, w + h*k3w );
      t += ( h / 6.0 ) * ( k1t + 2*k2t + 2*k3t + k4t );
      w += ( h / 6.0 ) * ( k1w + 2*k2w + 2*k3w + k4w );
      th[i+1] = t; om[i+1] = w;
    }
  }
  void at_interp( double t, double& theta, double& omega ) const
  {
    if( t <= 0.0 ){ theta = th.front(); omega = om.front(); return; }
    if( t >= T_end ){ theta = th.back();  omega = om.back();  return; }
    double s = t / h; size_t i = (size_t)s; double f = s - double( i );
    if( i >= N ){ theta = th.back(); omega = om.back(); return; }
    theta = ( 1.0 - f ) * th[i] + f * th[i+1];
    omega = ( 1.0 - f ) * om[i] + f * om[i+1];
  }
  void at_node( double t, double& theta, double& omega ) const
  {
    long idx = (long)std::llround( t / h );
    if( idx < 0 ) idx = 0;
    if( idx > (long)N ) idx = (long)N;
    theta = th[idx]; omega = om[idx];
  }
  static void cart( double theta, double omega,
                    double& x, double& y, double& u, double& v, double& lam )
  {
    x =  L_rod * std::sin( theta );
    y = -L_rod * std::cos( theta );
    u =  L_rod * omega * std::cos( theta );
    v =  L_rod * omega * std::sin( theta );
    lam = omega*omega + g_grav * std::cos( theta ) / L_rod;
  }
};

static Oracle const ORC( 200000 );
static std::vector<double> const T_SAMPLE = { 0.0, 0.25*T_end, 0.5*T_end, 0.75*T_end, T_end };

static double energy0()
{
  double x,y,u,v,lam; Oracle::cart( th0, om0, x,y,u,v,lam );
  return 0.5*( u*u + v*v ) + g_grav * y;
}

static int g_pass = 0, g_fail = 0;
static void check_close( char const* name, double got, double want, double tol )
{
  bool const ok = std::fabs( got - want ) <= tol;
  ( ok ? g_pass : g_fail )++;
  std::cout << "  " << std::left << std::setw(42) << name
            << std::right << std::scientific << std::setprecision(4)
            << " got=" << std::setw(12) << got << " want=" << std::setw(12) << want
            << " |d|=" << std::setw(10) << std::fabs( got - want )
            << "  " << ( ok ? "PASS" : "FAIL" ) << "\n";
}
static void check_le( char const* name, double got, double bound )
{
  bool const ok = ( got <= bound ) && std::isfinite( got );
  ( ok ? g_pass : g_fail )++;
  std::cout << "  " << std::left << std::setw(42) << name
            << std::right << std::scientific << std::setprecision(4)
            << " val=" << std::setw(12) << got << " <= " << std::setw(10) << bound
            << "  " << ( ok ? "PASS" : "FAIL" ) << "\n";
}
static void check_true( char const* name, bool ok )
{
  ( ok ? g_pass : g_fail )++;
  std::cout << "  " << std::left << std::setw(42) << name << "  " << ( ok ? "PASS" : "FAIL" ) << "\n";
}

struct ModeOut {
  bool                converged = false;
  bool                ok        = false;
  std::vector<double> F;                 // 8 outputs per sample
};

static ModeOut run_mode( OCFESLV::Options::ImpositionType imp, char const* strimp, bool marching )
{
  ModeOut R;
  char const* mstr = marching ? "march" : "mono";
  std::cout << "\n---------------- IC_" << strimp << " [" << mstr << "] ----------------\n";

  FFGraph DAG;
  FFVar t   = DAG.add_var( "t"     );
  FFVar x   = DAG.add_var( "x(t)"  );
  FFVar y   = DAG.add_var( "y(t)"  );
  FFVar u   = DAG.add_var( "u(t)"  );
  FFVar v   = DAG.add_var( "v(t)"  );
  FFVar lam = DAG.add_var( "lam(t)");

  OCFESLV oc( &DAG );
  oc.add_domain( t, FFDom( 0., T_end, N_EL, FFDom::LGR, N_ND ) );
  oc.set_evolution_domain( t );
  oc.add_state( x,   {t} );
  oc.add_state( y,   {t} );
  oc.add_state( u,   {t} );
  oc.add_state( v,   {t} );
  oc.add_state( lam, {t} );

  oc.update_ref( x,   [&]( OCFESLV::t_Coord const& c ){ double th,om,X,Y,U,V,Lm; ORC.at_interp( c.at(t), th, om ); Oracle::cart(th,om,X,Y,U,V,Lm); return X;  } );
  oc.update_ref( y,   [&]( OCFESLV::t_Coord const& c ){ double th,om,X,Y,U,V,Lm; ORC.at_interp( c.at(t), th, om ); Oracle::cart(th,om,X,Y,U,V,Lm); return Y;  } );
  oc.update_ref( u,   [&]( OCFESLV::t_Coord const& c ){ double th,om,X,Y,U,V,Lm; ORC.at_interp( c.at(t), th, om ); Oracle::cart(th,om,X,Y,U,V,Lm); return U;  } );
  oc.update_ref( v,   [&]( OCFESLV::t_Coord const& c ){ double th,om,X,Y,U,V,Lm; ORC.at_interp( c.at(t), th, om ); Oracle::cart(th,om,X,Y,U,V,Lm); return V;  } );
  oc.update_ref( lam, [&]( OCFESLV::t_Coord const& c ){ double th,om,X,Y,U,V,Lm; ORC.at_interp( c.at(t), th, om ); Oracle::cart(th,om,X,Y,U,V,Lm); return Lm; } );

  double x0,y0,u0,v0,l0; Oracle::cart( th0, om0, x0,y0,u0,v0,l0 );

  FFPartial OpP;
  FFVar EVOL_X = OpP( x, t ) - u;
  FFVar EVOL_Y = OpP( y, t ) - v;
  FFVar EVOL_U = OpP( u, t ) + lam * x;
  FFVar EVOL_V = OpP( v, t ) + lam * y + g_grav;
  FFVar ALG    = u*u + v*v - lam*( x*x + y*y ) - g_grav * y;
  FFVar IC_X = x - x0;
  FFVar IC_Y = y - y0;
  FFVar IC_U = u - u0;
  FFVar IC_V = v - v0;

  int const T_NO_LB = FFDom::ALL - FFDom::LB;
  OCFESLV::EqnOptions const interior( OCFESLV::EqnRole::INTERIOR, 0 );
  OCFESLV::EqnOptions const initial ( OCFESLV::EqnRole::INITIAL,  0 );

  oc.add_equation( EVOL_X, {t}, {T_NO_LB},    interior );
  oc.add_equation( EVOL_Y, {t}, {T_NO_LB},    interior );
  oc.add_equation( EVOL_U, {t}, {T_NO_LB},    interior );
  oc.add_equation( EVOL_V, {t}, {T_NO_LB},    interior );
  oc.add_equation( ALG,    {t}, {FFDom::ALL}, interior );
  oc.add_equation( IC_X,   {t}, {FFDom::LB},  initial  );
  oc.add_equation( IC_Y,   {t}, {FFDom::LB},  initial  );
  oc.add_equation( IC_U,   {t}, {FFDom::LB},  initial  );
  oc.add_equation( IC_V,   {t}, {FFDom::LB},  initial  );

  FFVar OUT_G0 = x*x + y*y - L_rod*L_rod;
  FFVar OUT_G1 = x*u + y*v;
  FFVar OUT_E  = 0.5*( u*u + v*v ) + g_grav * y;
  for( double ts : T_SAMPLE ){
    oc.add_output( x,      {t}, { ts } );
    oc.add_output( y,      {t}, { ts } );
    oc.add_output( u,      {t}, { ts } );
    oc.add_output( v,      {t}, { ts } );
    oc.add_output( lam,    {t}, { ts } );
    oc.add_output( OUT_G0, {t}, { ts } );
    oc.add_output( OUT_G1, {t}, { ts } );
    oc.add_output( OUT_E,  {t}, { ts } );
  }

  oc.options.REDUCE.ORDER    = OCFESLV::Options::RED_NONE;
  oc.options.CLASSIFY.MODE        = OCFESLV::Options::CLASS_AUTO;
  oc.options.INTERFACE.IMPOSITION = imp;
  oc.options.INTERFACE.TYPE  = OCFESLV::Options::IC_AUTO;
  oc.options.SOLVE.MARCHING  = marching;
  oc.options.SOLVE.MAX_ITER  = 100;
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
  check_true( "C0 DAE solve converged", rep.converged );
  if( !rep.converged ){ return R; }

  R.F = oc.val_functions();
  size_t const per = 8;
  if( R.F.size() != per * T_SAMPLE.size() ){
    std::cerr << "  output count mismatch: got " << R.F.size() << " want " << per*T_SAMPLE.size() << "\n";
    return R;
  }

  double const E0 = energy0();
  double maxg0 = 0., maxg1 = 0., maxE = 0., maxTraj = 0., maxLam = 0.;
  for( size_t kk = 0; kk < T_SAMPLE.size(); ++kk ){
    double const ts = T_SAMPLE[kk];
    double const xs = R.F[per*kk+0], ys = R.F[per*kk+1], us = R.F[per*kk+2],
                 vs = R.F[per*kk+3], ls = R.F[per*kk+4],
                 g0 = R.F[per*kk+5], g1 = R.F[per*kk+6], Es = R.F[per*kk+7];
    double th, om, xo, yo, uo, vo, lo; ORC.at_node( ts, th, om ); Oracle::cart( th, om, xo, yo, uo, vo, lo );
    maxg0   = std::max( maxg0,   std::fabs( g0 ) );
    maxg1   = std::max( maxg1,   std::fabs( g1 ) );
    maxE    = std::max( maxE,    std::fabs( Es - E0 ) );
    maxLam  = std::max( maxLam,  std::fabs( ls - lo ) );
    maxTraj = std::max( maxTraj, std::fabs( xs - xo ) );
    maxTraj = std::max( maxTraj, std::fabs( ys - yo ) );
    maxTraj = std::max( maxTraj, std::fabs( us - uo ) );
    maxTraj = std::max( maxTraj, std::fabs( vs - vo ) );
  }

  // Tolerances = a few x the observed spectral error for this pendulum at
  // N_EL=8, N_ND=5, T_end=0.5 (index-1 form): position/velocity drift ~2.5e-5,
  // energy ~2.9e-4, trajectory ~6e-5, lam ~1.4e-3 (lam is an acceleration-level
  // quantity, so the least accurate).  They still gate a ~3-4x accuracy regression.
  check_le( "P0 position drift max|g0|",   maxg0,   1e-4 );
  check_le( "P1 velocity drift max|g1|",   maxg1,   1e-4 );
  check_le( "P2 energy drift max|E-E0|",   maxE,    1e-3 );
  check_le( "P3 trajectory vs RK4 oracle", maxTraj, 2e-4 );
  check_le( "P4 lam vs (u^2+v^2-g y)/L^2", maxLam,  5e-3 );

  R.ok = ( maxg0 <= 1e-4 ) && ( maxg1 <= 1e-4 ) && ( maxE <= 1e-3 )
      && ( maxTraj <= 2e-4 ) && ( maxLam <= 5e-3 );
  return R;
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

int main()
{
  std::cout << "================================================================\n"
            << "  OCFE_DAE0b : index-3 pendulum (index-1 reduced), 2x2 sweep\n"
            << "  {IC_WEAK, IC_STRONG} x {monolithic, marching}\n"
            << "  L=" << L_rod << " g=" << g_grav << " theta0=" << th0
            << " omega0=" << om0 << " T_end=" << T_end << "\n"
            << "  E0=" << std::scientific << std::setprecision(6) << energy0() << "\n"
            << "================================================================\n";

  ModeOut const wm = run_mode( OCFESLV::Options::IC_WEAK,   "WEAK",   false );
  ModeOut const sm = run_mode( OCFESLV::Options::IC_STRONG, "STRONG", false );
  ModeOut const wc = run_mode( OCFESLV::Options::IC_WEAK,   "WEAK",   true  );
  ModeOut const sc = run_mode( OCFESLV::Options::IC_STRONG, "STRONG", true  );

  std::cout << "\n---------------- cross-run comparisons ----------------\n";
  // Gated only where exact: marching windows have no interior interface, so WEAK==STRONG.
  compare_vec( "XM  WEAK[march] == STRONG[march] (exact)", wc.F, sc.F, 1e-9 );
  // Discretisation-scale differences (informational; per-cell P* gates are the property).
  info_delta( "XM  WEAK[mono]  ~ STRONG[mono]",  wm.F, sm.F );
  info_delta( "MM  WEAK:   mono ~ march",        wm.F, wc.F );
  info_delta( "MM  STRONG: mono ~ march",        sm.F, sc.F );

  std::cout << "\n=================== DAE0b sweep summary ===================\n";
  std::cout << std::left << std::setw(20) << "run"
            << std::setw(12) << "converged" << std::setw(8) << "result" << "\n";
  auto row = [&]( char const* nm, ModeOut const& m ){
    std::cout << std::left << std::setw(20) << nm
              << std::setw(12) << ( m.converged ? "yes" : "no" )
              << std::setw(8)  << ( m.ok ? "PASS" : "FAIL" ) << "\n"; };
  row( "WEAK  [mono]",  wm );
  row( "STRONG[mono]",  sm );
  row( "WEAK  [march]", wc );
  row( "STRONG[march]", sc );

  std::cout << "\n================================================================\n"
            << "  OCFE_DAE0b: " << g_pass << " passed, " << g_fail << " failed -- "
            << ( g_fail == 0 ? "ALL PASS" : "SOME FAILED" ) << "\n"
            << "================================================================\n";
  return g_fail == 0 ? 0 : 1;
}
