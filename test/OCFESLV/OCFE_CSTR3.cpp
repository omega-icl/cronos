// OCFE_CSTR_gic.cpp  ---  General-IC marching: framework-automated consistent init (option b, proper)
// ===========================================================================
// The user writes the GENERAL steady-state INITIAL equation directly, with cf0 a CONSTANT:
//
//   INITIAL (LB):  Dil*(cf0 - c) - r = 0        (steady state at cf0; cf0 is a literal)
//   EVOLUTION    :  dc/dt = Dil*(cf - c) - r     (dynamics at the input feed cf)
//   ALGEBRAIC    :  r = k*c/(1+K*c)
//
//   (no transfer designation needed -- the framework auto-detects c as the differential state
//    closed by the general INITIAL equation and registers the override transfer itself)
//
// The framework runs the INITIAL equation as-is at window 0 (consistent initialisation) and, via
// the residual-row provenance (_rowEqnNdx), overrides c's INITIAL row with value-continuity
// c(LB) - terminal at every k>0.  No selector, no stiffness, no user-written blend.
//
// Validated against OCFE_CSTR_march_v2 (canonical IC, precomputed c_ss) -- the two are
// mathematically identical, so this must reproduce the same transient.
//
//   G0  marched solve converged (window-0 general IC + k>0 continuity override)
//   G1  c(0) == cf0 steady state   (consistent init from the general IC, captured at window 0)
//   G2  r(0) == k c0/(1+K c0)      (algebraic recovery at the LB)
//   G3  c(T_end) ~ cf steady state (relaxation, captured at the last window)
//   G4  forward == adjoint reduced Jacobian
//   G5..G7  dF/dcf vs central FD for c(0), c(T_end), r(0)  (sensitivity THROUGH the override)
//
// NOTE: point outputs are captured in the window whose evolution range contains their t (c(0) at
// window 0), and the reduced-space sensitivity linearises the k>0 continuity rows (not the general
// IC), with the terminal coupling carried by the override RHS (forward) / adjoint (reverse).
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

static double const Dil = 1.0, krate = 2.0, Ksat = 1.0;
static double const cf0_const = 1.0;   // pre-step feed -- a CONSTANT in the INITIAL equation
static double const cf_nom = 2.0, T_end = 5.0;

static double c_steady( double u )
{ double c = 0.5*u;
  for( int it=0; it<100; ++it ){
    double f  = Dil*( u - c ) - krate*c/( 1.0 + Ksat*c );
    double df = -Dil - krate/( ( 1.0+Ksat*c )*( 1.0+Ksat*c ) );
    double dc = -f/df; c += dc; if( std::fabs(dc)<1e-14 ) break; }
  return c; }
static double r_of_c( double c ){ return krate*c/( 1.0 + Ksat*c ); }

static int g_pass = 0, g_fail = 0;
static void check_true( char const* nm, bool ok )
{ std::cout << "  " << std::left << std::setw(46) << nm << ( ok ? " PASS" : " FAIL" ) << "\n"; ok?++g_pass:++g_fail; }
static void check_close( char const* nm, double got, double want, double tol )
{ double e = std::fabs( got - want );
  std::cout << "  " << std::left << std::setw(46) << nm << " |got-want|=" << std::scientific << std::setprecision(3) << e
            << " tol=" << tol << ( e<=tol ? "  PASS" : "  FAIL" ) << std::defaultfloat << "\n"; (e<=tol)?++g_pass:++g_fail; }

int main()
{
  std::cout << "================================================================\n"
            << "  General-IC marching CSTR -- framework-automated consistent init\n"
            << "  Dil*(cf0-c)-r=0 at LB (cf0=" << cf0_const << " constant) -> dynamics at cf=" << cf_nom << "\n"
            << "================================================================\n\n";

  size_t const n_el = 6, n_nd = 5;
  FFDom::TYPE const coltype = FFDom::LGL;

  FFGraph DAG;
  FFVar t  = DAG.add_var( "t" );
  FFVar c  = DAG.add_var( "c(t)" );
  FFVar r  = DAG.add_var( "r(t)" );
  FFVar cf = DAG.add_var( "cf" );        // dynamics feed -- decision

  OCFESLV oc( &DAG );
  oc.add_domain( t, FFDom( 0., T_end, n_el, coltype, n_nd ) );
  oc.set_evolution_domain( t );          // REQUIRED for marching

  oc.add_state( c, {t} );
  oc.add_state( r, {t} );
  oc.add_input( cf, cf_nom, /*is_decision=*/true );

  double const css0 = c_steady( cf0_const );
  oc.update_ref( c, [&]( OCFESLV::t_Coord const& ){ return css0; } );
  oc.update_ref( r, [&]( OCFESLV::t_Coord const& ){ return r_of_c( css0 ); } );

  FFPartial OpP;
  FFVar EVOL = OpP( c, t ) - ( Dil*( cf - c ) - r );
  FFVar ALG  = r - krate*c/( 1.0 + Ksat*c );
  FFVar IC   = Dil*( cf0_const - c ) - r;                     // GENERAL steady-state IC, cf0 constant

  int const T_NO_LB = FFDom::ALL - FFDom::LB;
  oc.add_equation( EVOL, {t}, {T_NO_LB},    OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( ALG,  {t}, {FFDom::ALL}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( IC,   {t}, {FFDom::LB},  OCFESLV::EqnOptions( OCFESLV::EqnRole::INITIAL,  0 ) );

  // No transfer designation: the framework auto-detects that c is a marching-continuous
  // differential state (it carries d/dt in EVOL) closed by a GENERAL INITIAL equation, and
  // registers the override-path transfer itself.  The user writes only the model + evolution domain.

  oc.add_output( c, {t}, {0.0}   );   // F0 = c(0)
  oc.add_output( c, {t}, {T_end} );   // F1 = c(T_end)
  oc.add_output( r, {t}, {0.0}   );   // F2 = r(0)

  if( !oc.setup() ){ std::cerr << "  setup() FAILED\n"; return 1; }

  std::vector<double> xv, inp;
  if( !oc.init( xv, inp, nullptr ) ){ std::cerr << "  init() FAILED\n"; return 1; }
  size_t const ncd = oc.n_control_dof(), ncf = oc.n_colloc_fct();
  auto const& C = oc.controls();
  size_t const icf = C.at( cf ).offset;
  std::vector<double> p0; oc.encode_controls( inp.data(), p0 );

  std::cout << "  is_marching()=" << ( oc.is_marching() ? "true" : "false" )
            << " ;  n_march_steps=" << oc.n_march_steps()
            << " ;  ncd=" << ncd << " (cf@" << icf << ")  ncf=" << ncf << "\n\n";
  check_true( "G-pre is_marching()", oc.is_marching() );

  // ---- marched solve (window 0 general IC -> consistent init; k>0 override -> continuity) ----
  std::vector<double> xm( xv ), im( inp );
  OCFESLV::SolveReport const repM = oc.solve( xm.data(), im.data(), nullptr );
  std::vector<double> Fm = oc.val_functions();
  check_true( "G0 marched solve converged", repM.converged );
  if( !repM.converged ){ std::cerr << "  abort\n"; return 1; }

  // ---- PRIMAL consistent-initialisation, validated against analytic steady states ----
  double const c0 = Fm[0], cT = Fm[1], r0 = Fm[2];
  check_close( "G1 c(0) == cf0 steady state (consistent init)", c0, css0, 1e-8 );
  check_close( "G2 r(0) == k c0/(1+K c0)  (algebraic recovery)", r0, r_of_c(c0), 1e-9 );
  check_close( "G3 c(T_end) ~ cf steady state (relaxation)",     cT, c_steady(cf_nom), 1e-3 );

  // ---- Reduced-space sensitivity THROUGH the general-IC override (fsens/asens block Jacobian uses
  // the continuity row at k>0, terminal coupling carried by the override RHS/adjoint). ----
  std::vector<std::vector<double>> Jf( ncf, std::vector<double>( ncd, 0. ) ), Ja( ncf, std::vector<double>( ncd, 0. ) );
  { std::vector<double> x(xv), i(inp); if(!oc.solve_fsens(x.data(),i.data(),nullptr)){ std::cerr<<"  solve_fsens FAILED\n"; return 1; }
    auto const& J=oc.sens_jacobian(); for(size_t a=0;a<ncf;++a) for(size_t b=0;b<ncd;++b) Jf[a][b]=J[a*ncd+b]; }
  { std::vector<double> x(xv), i(inp); if(!oc.solve_asens(x.data(),i.data(),nullptr)){ std::cerr<<"  solve_asens FAILED\n"; return 1; }
    auto const& J=oc.sens_jacobian(); for(size_t a=0;a<ncf;++a) for(size_t b=0;b<ncd;++b) Ja[a][b]=J[a*ncd+b]; }
  { double e=0.; for(size_t a=0;a<ncf;++a) for(size_t b=0;b<ncd;++b) e=std::max(e,std::fabs(Jf[a][b]-Ja[a][b]));
    check_close( "G4 forward == adjoint reduced Jacobian", e, 0., 1e-9 ); }

  // ---- central FD on cf vs the reduced Jacobian, for every output (c(0), c(T_end), r(0)) ----
  {
    double const h=1e-6;
    std::vector<double> pp(p0), pm(p0); pp[icf]+=h; pm[icf]-=h;
    std::vector<double> ip(inp), im2(inp), xp(xv), xm2(xv);
    oc.decode_controls(pp, ip.data());  oc.solve(xp.data(), ip.data(), nullptr);  std::vector<double> Fp=oc.val_functions();
    oc.decode_controls(pm, im2.data()); oc.solve(xm2.data(),im2.data(),nullptr);  std::vector<double> Fn=oc.val_functions();
    char const* nm[3] = { "G5 dC(0)/dcf ~ FD", "G6 dC(T_end)/dcf ~ FD", "G7 dR(0)/dcf ~ FD" };
    for( size_t a=0; a<ncf && a<3; ++a ){
      double const fd = ( Fp[a]-Fn[a] )/( 2.0*h );
      double const tol = 5e-3*std::max( 1e-3, std::fabs(Jf[a][icf]) ) + 1e-6;
      check_close( nm[a], Jf[a][icf], fd, tol );
    }
  }

  std::cout << "\n============================================================\n"
            << "  RESULT: " << g_pass << " passed, " << g_fail << " failed  -- "
            << ( g_fail==0 ? "ALL PASS" : "SOME FAILED" ) << "\n"
            << "============================================================\n";
  return g_fail ? 1 : 0;
}
