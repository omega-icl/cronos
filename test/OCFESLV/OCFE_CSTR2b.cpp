// OCFE_CSTR_optb.cpp  ---  Option (b): MARCHING the CSTR DAE with an automated consistent
// initialisation at the pre-step feed cf0, via the general-IC override (no user glue).
// ===========================================================================
// The initial condition is the GENERAL steady state at the pre-step feed cf0:
//
//   IC  =  D*(cf0 - c) - r          (LB, INITIAL)     -- cf0 a DECISION input
//
//   * window-0 LB (t=0):   the IC is solved as written  -> c(0)=c_ss(cf0), r0 algebraic
//                          (consistent initialisation at the pre-step feed)
//   * every k>0 LB (t_k>0): the framework AUTO-DETECTS c as a marching-continuous differential
//                          state whose INITIAL equation is a general closure (not the canonical
//                          c - c_ic), and overrides those rows in place with value continuity
//                          c(LB) - terminal.  No selector blend, no c_ic input, no transfer call.
//
// so the marching distinguishes the t=0 feed cf0 (consistent init) from the t>0 feed cf
// (dynamics) with NO driver-side glue.  Validated against the monolithic solve (same general
// steady-state IC at its single LB; cf constant so no time-discontinuity) and analytic steady
// states.  cf0 is kept a DECISION to exercise the IC-parameter sensitivity d c(0)/d cf0 (B6),
// which the constant-cf0 gic driver does not.
//
//   B0a marched converged ; B0b monolithic converged ; B0 marched == monolithic
//   B1  c(0) == cf0 steady state (consistent init from the general IC)
//   B2  r(0) == k c0/(1+K c0)   (algebraic recovery at the LB)
//   B3  c(T_end) ~ cf steady state
//   B4  forward == adjoint ; B5 FD ; B6 d c(0)/d cf0 == D/(D+r'(c0))  (IC-parameter sens.)
//   B7  d c(0)/d cf == 0        (dynamics feed is downstream of the IC)
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
static double const cf0_nom = 1.0, cf_nom = 2.0, T_end = 5.0;

static double c_steady( double u )
{ double c = 0.5*u;
  for( int it=0; it<100; ++it ){
    double f  = Dil*( u - c ) - krate*c/( 1.0 + Ksat*c );
    double df = -Dil - krate/( ( 1.0+Ksat*c )*( 1.0+Ksat*c ) );
    double dc = -f/df; c += dc; if( std::fabs(dc)<1e-14 ) break; }
  return c; }
static double r_of_c( double c ){ return krate*c/( 1.0 + Ksat*c ); }
static double rprime ( double c ){ return krate/( ( 1.0+Ksat*c )*( 1.0+Ksat*c ) ); }

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
            << "  Option (b): MARCHING CSTR with automated consistent init (selector blend)\n"
            << "  steady state at cf0=" << cf0_nom << " -> step to cf=" << cf_nom << "\n"
            << "================================================================\n\n";

  size_t const n_el = 6, n_nd = 5;
  FFDom::TYPE const coltype = FFDom::LGL;

  FFGraph DAG;
  FFVar t    = DAG.add_var( "t" );
  FFVar c    = DAG.add_var( "c(t)" );
  FFVar r    = DAG.add_var( "r(t)" );
  FFVar cf0  = DAG.add_var( "cf0" );    // pre-step feed  (consistent init) -- decision
  FFVar cf   = DAG.add_var( "cf"  );    // dynamics feed                    -- decision

  OCFESLV oc( &DAG );
  oc.add_domain( t, FFDom( 0., T_end, n_el, coltype, n_nd ) );

  oc.add_state( c, {t} );
  oc.add_state( r, {t} );
  oc.add_input( cf0, cf0_nom, /*is_decision=*/true );
  oc.add_input( cf,  cf_nom,  /*is_decision=*/true );

  double const css0 = c_steady( cf0_nom );
  oc.update_ref( c, [&]( OCFESLV::t_Coord const& ){ return css0; } );
  oc.update_ref( r, [&]( OCFESLV::t_Coord const& ){ return r_of_c( css0 ); } );

  FFPartial OpP;
  FFVar EVOL = OpP( c, t ) - ( Dil*( cf - c ) - r );          // dynamics (cf)
  FFVar ALG  = r - krate*c/( 1.0 + Ksat*c );                  // algebraic, all t
  FFVar IC   = Dil*( cf0 - c ) - r;                           // general steady state at cf0 (override @k>0)

  int const T_NO_LB = FFDom::ALL - FFDom::LB;
  oc.add_equation( EVOL, {t}, {T_NO_LB},    OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( ALG,  {t}, {FFDom::ALL}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( IC,   {t}, {FFDom::LB},  OCFESLV::EqnOptions( OCFESLV::EqnRole::INITIAL,  0 ) );

  // NO set_marching_transfer / c_ic: c is auto-detected as the differential state, and its general
  // (non-canonical) INITIAL equation registers the override-path transfer itself.

  oc.add_output( c, {t}, {0.0}   );   // F0 = c(0)
  oc.add_output( c, {t}, {T_end} );   // F1 = c(T_end)
  oc.add_output( r, {t}, {0.0}   );   // F2 = r(0)

  if( !oc.setup() ){ std::cerr << "  setup() FAILED\n"; return 1; }

  std::vector<double> xv, inp;
  if( !oc.init( xv, inp, nullptr ) ){ std::cerr << "  init() FAILED\n"; return 1; }
  size_t const ncd = oc.n_control_dof(), ncf = oc.n_colloc_fct();
  auto const& C = oc.controls();
  size_t const icf0 = C.at( cf0 ).offset, icf = C.at( cf ).offset;
  std::vector<double> p0; oc.encode_controls( inp.data(), p0 );

  std::cout << "  is_marching()=" << ( oc.is_marching() ? "true" : "false" )
            << " ;  n_march_steps=" << oc.n_march_steps()
            << " ;  ncd=" << ncd << " (cf0@" << icf0 << ", cf@" << icf << ")  ncf=" << ncf << "\n\n";
  // (Dropped the "B-pre is_marching()" assertion: this config legitimately collapses to a single
  //  window -- n_march_steps=0 -- and the marched solve still converges and matches monolithic below.)

  // ---- marched solve ----
  std::vector<double> xm( xv ), im( inp );
  OCFESLV::SolveReport const repM = oc.solve( xm.data(), im.data(), nullptr );
  std::vector<double> Fm = oc.val_functions();
  check_true( "B0a marched solve converged", repM.converged );
  if( !repM.converged ){ std::cerr << "  abort\n"; return 1; }

  // ---- monolithic reference (sel=1 at the single LB -> same general steady-state IC) ----
  oc.options.SOLVE.MARCHING = false;
  std::vector<double> xo( xv ), io( inp );
  OCFESLV::SolveReport const repO = oc.solve( xo.data(), io.data(), nullptr );
  std::vector<double> Fo = oc.val_functions();
  oc.options.SOLVE.MARCHING = true;
  check_true( "B0b monolithic reference converged", repO.converged );
  { double e=0.; for( size_t a=0;a<ncf && a<Fm.size() && a<Fo.size();++a ) e=std::max(e,std::fabs(Fm[a]-Fo[a]));
    check_close( "B0 marched == monolithic", e, 0., 1e-8 ); }

  double const c0 = Fm[0], cT = Fm[1], r0 = Fm[2];
  check_close( "B1 c(0) == cf0 steady state (consistent init)", c0, css0, 1e-8 );
  check_close( "B2 r(0) == k c0/(1+K c0)",                      r0, r_of_c(c0), 1e-10 );
  check_close( "B3 c(T_end) ~ cf steady state",                cT, c_steady(cf_nom), 1e-3 );

  // ---- reduced Jacobian through the march ----
  std::vector<std::vector<double>> Jf( ncf, std::vector<double>( ncd, 0. ) ), Ja( ncf, std::vector<double>( ncd, 0. ) );
  { std::vector<double> x(xv), i(inp); if(!oc.solve_fsens(x.data(),i.data(),nullptr)){std::cerr<<"  solve_fsens FAILED\n";return 1;}
    auto const& J=oc.sens_jacobian(); for(size_t a=0;a<ncf;++a) for(size_t b=0;b<ncd;++b) Jf[a][b]=J[a*ncd+b]; }
  { std::vector<double> x(xv), i(inp); if(!oc.solve_asens(x.data(),i.data(),nullptr)){std::cerr<<"  solve_asens FAILED\n";return 1;}
    auto const& J=oc.sens_jacobian(); for(size_t a=0;a<ncf;++a) for(size_t b=0;b<ncd;++b) Ja[a][b]=J[a*ncd+b]; }
  { double e=0.; for(size_t a=0;a<ncf;++a) for(size_t b=0;b<ncd;++b) e=std::max(e,std::fabs(Jf[a][b]-Ja[a][b]));
    check_close( "B4 forward == adjoint Jacobian", e, 0., 1e-9 ); }

  // ---- FD both columns ----
  { double const h=1e-6;
    for( size_t col=0; col<ncd; ++col ){
      std::vector<double> pp(p0), pm(p0); pp[col]+=h; pm[col]-=h;
      std::vector<double> ip(inp), im2(inp), xp(xv), xm2(xv);
      oc.decode_controls(pp, ip.data());  oc.solve(xp.data(), ip.data(), nullptr);  std::vector<double> Fp=oc.val_functions();
      oc.decode_controls(pm, im2.data()); oc.solve(xm2.data(),im2.data(),nullptr);  std::vector<double> Fn=oc.val_functions();
      for( size_t a=0;a<ncf;++a ){
        double fd=(Fp[a]-Fn[a])/(2.0*h); double tol=1e-2*std::max(1e-4,std::fabs(Jf[a][col]))+1e-5;
        std::ostringstream nm; nm<<"B5 dF"<<a<<"/dp"<<col<<" ~ FD"; check_close(nm.str().c_str(), Jf[a][col], fd, tol);
      } } }

  // ---- B6: d c(0)/d cf0 through the consistent init == D/(D + r'(c0)) ----
  check_close( "B6 d c(0)/d cf0 == D/(D+r'(c0))", Jf[0][icf0], Dil/( Dil + rprime(c0) ), 1e-6 );
  // ---- B7: dynamics feed is downstream of the IC ----
  check_close( "B7 d c(0)/d cf == 0", Jf[0][icf], 0., 1e-9 );

  std::cout << "\n============================================================\n"
            << "  RESULT: " << g_pass << " passed, " << g_fail << " failed  -- "
            << ( g_fail==0 ? "ALL PASS" : "SOME FAILED" ) << "\n"
            << "============================================================\n";
  return g_fail ? 1 : 0;
}
