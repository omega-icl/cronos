// OCFE_DAE3q.cpp  ---  DAE with PIECEWISE-QUADRATIC (n_node=3) time control
// ===========================================================================
// Extends OCFE_DAE3 (piecewise-LINEAR, n_node=2) to a piecewise-QUADRATIC input, to verify the
// causal evolution-interface fix PAST piecewise-linear.
//
//   differential   dc/dt = -a c + u(t)        (linear dynamics, quadratic forcing)
//   algebraic      r     = c^2                 (nonlinear output relation)
//   initial        c(0)  = c_ic                (value-continuity IC; marching transfer)
//
// u(t) is a DISTRIBUTED input control over {t} with its OWN n_node=3 (LGR) discretisation: on each
// element k it is the quadratic  u = A_k + B_k tau + C_k tau^2  (tau = t - t_k).  The A_k/B_k/C_k are
// chosen so u JUMPS at every interior element interface -> DISCONTINUOUS in the evolution direction,
// so evolution_input_may_jump = true and the causal guard is exercised.  The state carries the whole
// control ( ncd = n_el * 3 = 15 input DOFs, + 1 for c_ic ).
//
// Because u is exactly quadratic per element, the n_node=3 interpolant reproduces it EXACTLY, so the
// analytic recurrence below (integrating dc/dt = -a c + quadratic in closed form) is the exact
// solution of the DISCRETISED problem -- no reconstruction of node positions is needed.
//
// Checks, per imposition (IC_WEAK, IC_STRONG) x mode (mono, march):
//   Q0  solve converged
//   Q1  c(t_m) == analytic quadratic recurrence            [discretisation-limited]
//   Q2  r(T) == c(T)^2                                     [algebraic, exact]
//   Q3  forward-AD reduced Jacobian == adjoint             [tight]
//   Q4  forward-AD reduced Jacobian == central FD          [loose]
// Cross-cell:
//   MM  mono ~ march  ( |dFval|, |dJ| )   -- STRONG must be ~machine precision (the causal-fix goal)
// ===========================================================================
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <string>
#include <sstream>

#ifndef OCFE_OCFESLV_HEADER
#define OCFE_OCFESLV_HEADER "ocfeslv.hpp"
#endif
#include OCFE_OCFESLV_HEADER

#include "ffocfe.hpp"

using namespace mc;

// --- problem parameters ---
static double const a_decay = 1.0;     // linear decay rate
static double const c0_ic   = 0.5;     // initial condition c(0)
static size_t const Nel     = 5;       // evolution elements on [0,Nel], width 1
static double const T_end   = double(Nel);
static size_t const n_nd    = 7;       // STATE nodes (LGL) -- high order so c(t_m) resolves the recurrence
static size_t const u_nd    = 3;       // INPUT nodes (LGR) -- piecewise-QUADRATIC control

// per-element quadratic coefficients of u (discontinuous at interior interfaces by construction)
static double const Acoef[5] = { 1.0,  1.4,  0.8,  1.3,  0.7 };
static double const Bcoef[5] = { 0.6, -0.4,  0.7, -0.3,  0.5 };
static double const Ccoef[5] = { 0.3,  0.25,-0.20, 0.35,-0.15 };

static double u_fun( double t )
{
  int k = (int)std::floor( t + 1e-9 ); if( k < 0 ) k = 0; if( k > (int)Nel-1 ) k = (int)Nel-1;  // +1e-9: LB node t=k belongs to element k (returns A_k), not k-1
  double const tau = t - double(k);
  return Acoef[k] + Bcoef[k]*tau + Ccoef[k]*tau*tau;
}

// Exact solution of  dc/dt = -a c + (A + B tau + C tau^2)  on [t_k,t_k+h], marched element to element.
// Particular quadratic p(tau)=p0+p1 tau+p2 tau^2:  p2=C/a, p1=B/a-2C/a^2, p0=A/a-B/a^2+2C/a^3.
// c(h) = (c_k - p0) e^{-a h} + p0 + p1 h + p2 h^2.
static std::vector<double> recurrence()
{
  std::vector<double> cout( Nel, 0. );   // cout[k] = c at t = k+1
  double c = c0_ic, h = 1.0, a = a_decay;
  for( size_t k = 0; k < Nel; ++k ){
    double const A = Acoef[k], B = Bcoef[k], C = Ccoef[k];
    double const p2 = C/a;
    double const p1 = B/a - 2.0*C/(a*a);
    double const p0 = A/a - B/(a*a) + 2.0*C/(a*a*a);
    c = ( c - p0 )*std::exp( -a*h ) + p0 + p1*h + p2*h*h;
    cout[k] = c;
  }
  return cout;
}

static int g_pass = 0, g_fail = 0;
static void check_true( std::string const& nm, bool ok )
{ std::cout << "  " << std::left << std::setw(46) << nm << ( ok ? " PASS" : " FAIL" ) << "\n"; ok ? ++g_pass : ++g_fail; }
static void check_close( std::string const& nm, double got, double want, double tol )
{ double e = std::fabs( got - want );
  std::cout << "  " << std::left << std::setw(46) << nm << " |d|=" << std::scientific << std::setprecision(3) << e
            << " tol=" << tol << ( e <= tol ? "  PASS" : "  FAIL" ) << std::defaultfloat << "\n";
  ( e <= tol ) ? ++g_pass : ++g_fail; }

struct ModeResult {
  bool ok = false;
  size_t ncd = 0, ncf = 0;
  size_t cic_col = 0;                          // reduced-control column index of c_ic
  std::vector<double> F;                       // output values
  std::vector<std::vector<double>> Jf;         // forward-AD reduced Jacobian [ncf][ncd]
};

static ModeResult run_mode( OCFESLV::Options::ImpositionType imp, bool marching,
                            std::string const& tag )
{
  ModeResult R;

  FFGraph DAG;
  FFVar t    = DAG.add_var( "t" );
  FFVar c    = DAG.add_var( "c(t)" );      // differential state
  FFVar r    = DAG.add_var( "r(t)" );      // algebraic state  r = c^2
  FFVar u    = DAG.add_var( "u(t)" );      // distributed input control (piecewise-quadratic)
  FFVar c_ic = DAG.add_var( "c_ic" );      // IC value (marching transfer + decision)

  OCFESLV oc( &DAG );
  oc.options.INTERFACE.IMPOSITION = imp;
  oc.options.SOLVE.MARCHING  = marching;
  oc.options.DISPLAY_LEVEL   = 0;

  oc.add_domain( t, FFDom( 0., T_end, Nel, FFDom::LGL, n_nd ) );
  oc.set_evolution_domain( t );

  oc.add_state( c, {t} );
  oc.add_state( r, {t} );
  oc.add_input( u, {t}, FFDom::LGR, u_nd );          // piecewise-quadratic control, DISCONTINUOUS (default)
  oc.add_input( c_ic, c0_ic, /*is_decision=*/true ); // scalar IC transfer + decision

  oc.update_ref( u, [&]( OCFESLV::t_Coord const& cr ){ return u_fun( cr.at(t) ); } );
  oc.update_ref( c, [&]( OCFESLV::t_Coord const& ){ return c0_ic; } );
  oc.update_ref( r, [&]( OCFESLV::t_Coord const& ){ return c0_ic*c0_ic; } );

  FFPartial OpP;
  FFVar EVOL = OpP( c, t ) - ( -a_decay*c + u );
  FFVar ALG  = r - c*c;
  FFVar IC   = c - c_ic;

  int const T_NO_LB = FFDom::ALL - FFDom::LB;
  oc.add_equation( EVOL, {t}, {T_NO_LB},    OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( ALG,  {t}, {FFDom::ALL}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( IC,   {t}, {FFDom::LB},  OCFESLV::EqnOptions( OCFESLV::EqnRole::INITIAL,  0 ) );

  for( size_t m = 1; m <= Nel; ++m ) oc.add_output( c, {t}, { double(m) } );  // c(1)..c(Nel)
  oc.add_output( r, {t}, { double(Nel) } );                                    // r(T) = c(T)^2

  if( !oc.setup() ){ std::cerr << "  [" << tag << "] setup() FAILED\n"; return R; }
  oc.register_control( u );   // register u's n_el*u_nd DOFs (c_ic auto-registered via is_decision)

  std::vector<double> xv, inp;
  if( !oc.init( xv, inp, nullptr ) ){ std::cerr << "  [" << tag << "] init() FAILED\n"; return R; }

  R.ncd = oc.n_control_dof();
  R.ncf = oc.n_colloc_fct();
  R.cic_col = oc.controls().at( c_ic ).offset;   // which reduced column is the IC control
  std::vector<double> p0v; oc.encode_controls( inp.data(), p0v );

  // ---- primary solve ----
  std::vector<double> xs( xv ), is( inp );
  OCFESLV::SolveReport const rep = oc.solve( xs.data(), is.data(), nullptr );
  R.ok = rep.converged;
  R.F  = oc.val_functions(); R.F.resize( R.ncf );
  if( !rep.converged ) return R;

  // ---- forward-AD reduced Jacobian ----
  R.Jf.assign( R.ncf, std::vector<double>( R.ncd, 0. ) );
  { std::vector<double> x( xv ), i( inp );
    if( oc.solve_fsens( x.data(), i.data(), nullptr ) ){
      auto const& J = oc.sens_jacobian();
      for( size_t f = 0; f < R.ncf; ++f ) for( size_t j = 0; j < R.ncd; ++j ) R.Jf[f][j] = J[f*R.ncd + j];
    } else std::cerr << "  [" << tag << "] solve_fsens FAILED\n";
  }
  // ---- adjoint reduced Jacobian ----
  std::vector<std::vector<double>> Ja( R.ncf, std::vector<double>( R.ncd, 0. ) );
  { std::vector<double> x( xv ), i( inp );
    if( oc.solve_asens( x.data(), i.data(), nullptr ) ){
      auto const& J = oc.sens_jacobian();
      for( size_t f = 0; f < R.ncf; ++f ) for( size_t j = 0; j < R.ncd; ++j ) Ja[f][j] = J[f*R.ncd + j];
    } else std::cerr << "  [" << tag << "] solve_asens FAILED\n";
  }
  // ---- central-FD reduced Jacobian ----
  std::vector<std::vector<double>> Jfd( R.ncf, std::vector<double>( R.ncd, 0. ) );
  { double const h = 1e-6;
    for( size_t col = 0; col < R.ncd; ++col ){
      std::vector<double> pp( p0v ), pm( p0v ); pp[col] += h; pm[col] -= h;
      std::vector<double> ip( inp ), im( inp ), xp( xv ), xm( xv );
      oc.decode_controls( pp, ip.data() ); oc.solve( xp.data(), ip.data(), nullptr );
      std::vector<double> Fp = oc.val_functions(); Fp.resize( R.ncf );
      oc.decode_controls( pm, im.data() ); oc.solve( xm.data(), im.data(), nullptr );
      std::vector<double> Fm = oc.val_functions(); Fm.resize( R.ncf );
      for( size_t f = 0; f < R.ncf; ++f ) Jfd[f][col] = ( Fp[f] - Fm[f] )/( 2.0*h );
    }
  }

  // ================= checks =================
  std::cout << "  ---- " << tag << "  (ncd=" << R.ncd << " ncf=" << R.ncf
            << " marching=" << ( oc.is_marching() ? "yes" : "no" ) << ") ----\n";
  check_true( "Q0 solve converged", rep.converged );

  std::vector<double> const cref = recurrence();
  { double emax = 0.;
    for( size_t m = 0; m < Nel; ++m ) emax = std::max( emax, std::fabs( R.F[m] - cref[m] ) );
    check_close( "Q1 c(t_m) == quadratic recurrence", emax, 0., 1e-4 ); }

  check_close( "Q2 r(T) == c(T)^2", R.F[Nel], R.F[Nel-1]*R.F[Nel-1], 1e-10 );

  { double e = 0.; for( size_t f=0; f<R.ncf; ++f ) for( size_t j=0; j<R.ncd; ++j ) e = std::max( e, std::fabs( R.Jf[f][j]-Ja[f][j] ) );
    check_close( "Q3 reduced Jacobian fwd == adjoint", e, 0., 1e-9 ); }

  { double emax = 0.; size_t wf=0, wj=0;
    for( size_t f=0; f<R.ncf; ++f ) for( size_t j=0; j<R.ncd; ++j ){
      double const scale = 1e-2*std::max( 1e-4, std::fabs( R.Jf[f][j] ) ) + 1e-5;
      double const over = std::fabs( R.Jf[f][j]-Jfd[f][j] ) - scale;        // <=0 means within tol
      if( over > emax ){ emax = over; wf = f; wj = j; }
    }
    check_true( "Q4 reduced Jacobian fwd ~ central FD", emax <= 0. );
    if( emax > 0. )
      std::cout << "      worst dF" << wf << "/dp" << wj
                << ( wj==R.cic_col ? "  [== c_ic column]" : "  [u column]" )
                << "  fwd=" << std::scientific << std::setprecision(4) << R.Jf[wf][wj]
                << "  FD="  << Jfd[wf][wj] << std::defaultfloat << "\n";
  }

  return R;
}

static void cross_cell( std::string const& tag, ModeResult const& mono, ModeResult const& march )
{
  if( !mono.ok || !march.ok ){ check_true( tag + " mono~march (both converged)", false ); return; }
  double dF = 0., dJ = 0.; size_t wf=0, wj=0;
  for( size_t f=0; f<mono.ncf && f<march.ncf; ++f ) dF = std::max( dF, std::fabs( mono.F[f]-march.F[f] ) );
  for( size_t f=0; f<mono.ncf && f<march.ncf; ++f )
    for( size_t j=0; j<mono.ncd && j<march.ncd; ++j ){
      double const d = std::fabs( mono.Jf[f][j]-march.Jf[f][j] );
      if( d > dJ ){ dJ = d; wf = f; wj = j; }
    }
  std::cout << "  MM " << std::left << std::setw(8) << tag
            << " mono~march  |dFval|=" << std::scientific << std::setprecision(4) << dF
            << "  |dJ|=" << dJ << std::defaultfloat
            << "   (worst dF" << wf << "/dp" << wj
            << ( wj==mono.cic_col ? " == c_ic col" : " = u col" )
            << ": mono=" << std::scientific << std::setprecision(4) << mono.Jf[wf][wj]
            << " march=" << march.Jf[wf][wj] << std::defaultfloat << ")\n";
}

int main()
{
  std::cout << "================================================================\n"
            << "  OCFE_DAE3q : DAE with PIECEWISE-QUADRATIC (n_node=3) time control\n"
            << "  dc/dt = -a c + u(t),  r = c^2 ;  u piecewise-quadratic (3 DOFs/element)\n"
            << "  forward + adjoint sensitivity, FD-validated, mono vs march\n"
            << "================================================================\n\n";

  std::cout << "================ IC_WEAK ================\n";
  ModeResult const wMono  = run_mode( OCFESLV::Options::IC_WEAK,   false, "IC_WEAK  [mono]"  );
  ModeResult const wMarch = run_mode( OCFESLV::Options::IC_WEAK,   true,  "IC_WEAK  [march]" );

  std::cout << "\n================ IC_STRONG ================\n";
  ModeResult const sMono  = run_mode( OCFESLV::Options::IC_STRONG, false, "IC_STRONG[mono]"  );
  ModeResult const sMarch = run_mode( OCFESLV::Options::IC_STRONG, true,  "IC_STRONG[march]" );

  std::cout << "\n---------------- cross-cell comparisons ----------------\n";
  cross_cell( "WEAK",   wMono, wMarch );
  cross_cell( "STRONG", sMono, sMarch );

  std::cout << "\n============================================================\n"
            << "  RESULT: " << g_pass << " passed, " << g_fail << " failed  -- "
            << ( g_fail == 0 ? "ALL PASS" : "SOME FAILED" ) << "\n"
            << "============================================================\n";
  return g_fail ? 1 : 0;
}
