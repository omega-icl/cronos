// TRANS_psa.cpp -- OCFESLV transitions on the PSA1 model (OCFE_PSA1.cpp, declared EXACTLY as there), monolithic and
// marching (IC_STRONG, RED_FULL, 5x6 LGL, T_end = 5), with a PURGE STEP at tau = 2 (a window seam when marching):
//     c1(tau+, z) = (1 - d) c1(tau-, z),  d = 0.3        (every other state continuous; c1's IC is canonical: the data path)
// No closed form -- the oracles are internal: the identity (d = 0) == no transition; CAUSALITY (c1(tau-) unaffected by
// the jump: exact when marching, to discretisation level when monolithic); the jump c1(tau+,1/2) = (1-d) c1(tau-,1/2);
// marching == monolithic (values and gradients); forward == adjoint; dF/dd == central finite differences; and the
// physics (purging component 1 lowers the final inventory Inv1).
#include <cmath>
#include <cstdio>
#include <sstream>
#include <iostream>
#include <string>
#include <vector>
#include "ocfeslv.hpp"
using namespace mc;
static double const U_VEL = 1.0, D_ax = 0.1, F_ph = 0.5, k_ldf = 2.0;
static double const qs1 = 1.0, qs2 = 1.0, b02_L = 1.0;
static double const beta1 = 2.0, beta2 = 1.0, T0_ref = 1.0;
static double const Cp_e = 1.0, G_cv = 1.0, lam = 0.1, dH1 = 1.0, dH2 = 0.5, hw = 0.5, Tw = 1.0;
static double const c0_1 = 0.5, c0_2 = 0.5, tau_in = 0.15, T_end = 5.0;
static double const b01_nom = 3.0;

static inline double c1feed_d( double t ){ return c0_1*( 1.0 - std::exp( -t/tau_in ) ); }
static inline double c2feed_d( double t ){ return c0_2*( 1.0 - std::exp( -t/tau_in ) ); }
static inline double q1star_g( double c1g, double c2g, double b01 )
{ return qs1*b01*c1g/( 1.0 + b01*c1g + b02_L*c2g ); }
static inline double q2star_g( double c1g, double c2g, double b01 )
{ return qs2*b02_L*c2g/( 1.0 + b01*c1g + b02_L*c2g ); }

static double const TAU = 2.0, DJ = 0.3;
static int npass = 0, nfail = 0;
static void check( bool c, std::string const& w, double v = -1. ){ std::printf( "  %s  %-78s", c? "PASS": "FAIL", w.c_str() ); if( v >= 0. ) std::printf( " (%.2e)", v ); std::printf( "\n" ); c? ++npass: ++nfail; }
struct Out { bool ok = false; std::string msg; std::vector<double> F, Jf, Ja; };      // F: Inv1, Eff1, c1(tau-,.5), c1(tau+,.5)
//! kind 0: no transition; 1: the purge with d an input.  Controls (b01, d); Jacobians 4 x 2, row-major.
static Out run( int kind, double d, bool march, bool sens = true ){
  Out R;
  size_t const n_el = getenv( "PSA_NEL" )? atoi( getenv( "PSA_NEL" ) ): 5, n_nd = getenv( "PSA_NND" )? atoi( getenv( "PSA_NND" ) ): 6;   // refinement study
  FFDom::TYPE const coltype = FFDom::LGL;
  FFGraph DAG;
  FFVar t  = DAG.add_var( "t" ), z = DAG.add_var( "z" );
  FFVar c1 = DAG.add_var( "c1(t,z)" ), c2 = DAG.add_var( "c2(t,z)" ), q1 = DAG.add_var( "q1(t,z)" ), q2 = DAG.add_var( "q2(t,z)" ), T = DAG.add_var( "T(t,z)" );
  FFVar c1_ic = DAG.add_var( "c1_ic(z)" ), c2_ic = DAG.add_var( "c2_ic(z)" ), q1_ic = DAG.add_var( "q1_ic(z)" ), q2_ic = DAG.add_var( "q2_ic(z)" ), T_ic = DAG.add_var( "T_ic(z)" );
  FFVar b01 = DAG.add_var( "b01" ), pd = DAG.add_var( "d" );
  FFPartial OpP;  FFIntegral OpI;  FFEval OpE;
  FFVar c1feed = c0_1*( 1.0 - exp( -t/tau_in ) );
  FFVar c2feed = c0_2*( 1.0 - exp( -t/tau_in ) );
  FFVar b1 = b01  *exp( beta1*( 1.0/T - 1.0/T0_ref ) );
  FFVar b2 = b02_L*exp( beta2*( 1.0/T - 1.0/T0_ref ) );
  FFVar den = 1.0 + b1*c1 + b2*c2;
  FFVar q1star = qs1*b1*c1/den;
  FFVar q2star = qs2*b2*c2/den;

  FFVar CONT1 = OpP( c1, t ) + U_VEL*OpP( c1, z ) - D_ax*OpP( OpP( c1, z ), z ) + F_ph*OpP( q1, t );
  FFVar CONT2 = OpP( c2, t ) + U_VEL*OpP( c2, z ) - D_ax*OpP( OpP( c2, z ), z ) + F_ph*OpP( q2, t );
  FFVar LDF1  = OpP( q1, t ) - k_ldf*( q1star - q1 );
  FFVar LDF2  = OpP( q2, t ) - k_ldf*( q2star - q2 );
  FFVar ENE_T = Cp_e*OpP( T, t ) + G_cv*OpP( T, z ) - lam*OpP( OpP( T, z ), z )
              - F_ph*( dH1*OpP( q1, t ) + dH2*OpP( q2, t ) ) + hw*( T - Tw );

  FFVar IC_c1 = c1 - c1_ic;
  FFVar IC_c2 = c2 - c2_ic;
  FFVar IC_q1 = q1 - q1_ic;
  FFVar IC_q2 = q2 - q2_ic;
  FFVar IC_T  = T  - T_ic;

  FFVar BC_L1 = U_VEL*c1 - D_ax*OpP( c1, z ) - U_VEL*c1feed;
  FFVar BC_L2 = U_VEL*c2 - D_ax*OpP( c2, z ) - U_VEL*c2feed;
  FFVar BC_LT = G_cv*T  - lam*OpP( T,  z ) - G_cv*T0_ref;
  FFVar BC_U1 = OpP( c1, z );
  FFVar BC_U2 = OpP( c2, z );
  FFVar BC_UT = OpP( T,  z );

  FFVar Inv1 = OpI( c1 + F_ph*q1, z );   // fct 0 (terminal)
  FFVar Eff1 = OpI( U_VEL*c1, t );        // fct 1 (accumulated)

  OCFESLV oc( &DAG );
  oc.add_domain( t, FFDom( 0., T_end, n_el, coltype, n_nd ) );
  oc.add_domain( z, FFDom( 0., 1.0,   n_el, coltype, n_nd ) );
  oc.add_state ( c1, {t,z} );  oc.add_state ( c2, {t,z} );
  oc.add_state ( q1, {t,z} );  oc.add_state ( q2, {t,z} );
  oc.add_state ( T,  {t,z} );
  oc.add_input ( b01, {} );

  auto c1g = [&]( OCFESLV::t_Coord const& cr ){ return c1feed_d( cr.at(t) )*( 1.0 - 0.5*cr.at(z) ); };
  auto c2g = [&]( OCFESLV::t_Coord const& cr ){ return c2feed_d( cr.at(t) )*( 1.0 - 0.5*cr.at(z) ); };
  oc.update_ref( c1, c1g );
  oc.update_ref( c2, c2g );
  oc.update_ref( q1, [&]( OCFESLV::t_Coord const& cr ){ return q1star_g( c1g(cr), c2g(cr), b01_nom ); } );
  oc.update_ref( q2, [&]( OCFESLV::t_Coord const& cr ){ return q2star_g( c1g(cr), c2g(cr), b01_nom ); } );
  oc.update_ref( T,  [&]( OCFESLV::t_Coord const& cr ){
    return T0_ref + F_ph*( dH1*q1star_g(c1g(cr),c2g(cr),b01_nom) + dH2*q2star_g(c1g(cr),c2g(cr),b01_nom) )/Cp_e; } );

  oc.add_input ( c1_ic, {z} );  oc.add_input ( c2_ic, {z} );
  oc.add_input ( q1_ic, {z} );  oc.add_input ( q2_ic, {z} );  oc.add_input ( T_ic, {z} );
  oc.update_ref( c1_ic, [&]( OCFESLV::t_Coord const& cr ){ return 0.0;    } );
  oc.update_ref( c2_ic, [&]( OCFESLV::t_Coord const& cr ){ return 0.0;    } );
  oc.update_ref( q1_ic, [&]( OCFESLV::t_Coord const& cr ){ return 0.0;    } );
  oc.update_ref( q2_ic, [&]( OCFESLV::t_Coord const& cr ){ return 0.0;    } );
  oc.update_ref( T_ic,  [&]( OCFESLV::t_Coord const& cr ){ return T0_ref; } );

  oc.set_evolution_domain( t );

  int const T_NO_LB = FFDom::ALL - FFDom::LB;
  int const Z_INT   = FFDom::ALL - FFDom::LB - FFDom::UB;

  oc.add_equation( CONT1, {t,z}, {T_NO_LB,   Z_INT},      OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( CONT2, {t,z}, {T_NO_LB,   Z_INT},      OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( LDF1,  {t,z}, {T_NO_LB,   FFDom::ALL}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( LDF2,  {t,z}, {T_NO_LB,   FFDom::ALL}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( ENE_T, {t,z}, {T_NO_LB,   Z_INT},      OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( IC_c1, {t,z}, {FFDom::LB, FFDom::ALL}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INITIAL,  0 ) );
  oc.add_equation( IC_c2, {t,z}, {FFDom::LB, FFDom::ALL}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INITIAL,  0 ) );
  oc.add_equation( IC_q1, {t,z}, {FFDom::LB, FFDom::ALL}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INITIAL,  0 ) );
  oc.add_equation( IC_q2, {t,z}, {FFDom::LB, FFDom::ALL}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INITIAL,  0 ) );
  oc.add_equation( IC_T,  {t,z}, {FFDom::LB, FFDom::ALL}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INITIAL,  0 ) );
  oc.add_equation( BC_L1, {t,z}, {T_NO_LB,   FFDom::LB},  OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );
  oc.add_equation( BC_L2, {t,z}, {T_NO_LB,   FFDom::LB},  OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );
  oc.add_equation( BC_LT, {t,z}, {T_NO_LB,   FFDom::LB},  OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );
  oc.add_equation( BC_U1, {t,z}, {T_NO_LB,   FFDom::UB},  OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );
  oc.add_equation( BC_U2, {t,z}, {T_NO_LB,   FFDom::UB},  OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );
  oc.add_equation( BC_UT, {t,z}, {T_NO_LB,   FFDom::UB},  OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );
  oc.add_input( pd, {} );
  if( kind == 1 ) oc.add_transition( ( 1. - pd )*c1, c1, t, TAU );
  oc.add_output( Inv1, {t}, {T_end} );
  oc.add_output( Eff1, {z}, {1.0}   );
  if( getenv( "PSA_POINTFORM" ) ){                        // the same two values through POINT-form outputs (no auxiliary)
    oc.add_output( c1, {t,z}, {TAU, 0.5}, {FFDom::MINUS, FFDom::MINUS} );
    oc.add_output( c1, {t,z}, {TAU, 0.5}, {FFDom::PLUS,  FFDom::MINUS} );
  }
  else{
    oc.add_output( OpE( OpE( c1, t, TAU, FFDom::MINUS ), z, 0.5 ) );
    oc.add_output( OpE( OpE( c1, t, TAU, FFDom::PLUS  ), z, 0.5 ) );
  }
  oc.options.REDUCE.ORDER     = OCFESLV::Options::RED_FULL;
  oc.options.CLASSIFY.MODE    = OCFESLV::Options::CLASS_AUTO;
  oc.options.INTERFACE.TYPE   = OCFESLV::Options::IC_AUTO;
  oc.options.INTERFACE.IMPOSITION = OCFESLV::Options::IC_STRONG;
  oc.options.INTERFACE.SAT_SIGMA0 = 10.0;
  oc.options.DISPLAY_LEVEL    = 0;
#if defined(CRONOS__WITH_SPQR)
  oc.options.SOLVE.FACTORIZATION = OCFESLV::Options::SOLVE_SPQR;
#else
  oc.options.SOLVE.FACTORIZATION = OCFESLV::Options::SOLVE_SUPERLU;
#endif
  oc.options.SOLVE.MARCHING = march;
  if( march ) oc.options.SOLVE.WARMSTART = OCFESLV::Options::BROADCAST_IC;
  std::ostringstream os; auto* old = std::cerr.rdbuf( os.rdbuf() );
  std::vector<double> xv, inp; bool ok = oc.setup() && oc.init( xv, inp, nullptr );
  if( ok ){
    oc.set_input_values( b01, { b01_nom }, inp.data() );
    size_t const ntic = oc.get_input_values( T_ic, inp.data() ).size();
    oc.set_input_values( T_ic, std::vector<double>( ntic, T0_ref ), inp.data() );
    oc.set_input_values( pd, { d }, inp.data() );
    oc.register_control( b01 );  oc.register_control( pd );
  }
  std::vector<double> x1( xv ), x2( xv ), x3( xv );
  if( ok ) ok = oc.solve( x1.data(), inp.data(), nullptr ).converged;  if( ok ) R.F = oc.val_functions();
  if( ok && sens ){ ok = oc.solve_fsens( x2.data(), inp.data(), nullptr );  if( ok ) R.Jf = oc.sens_jacobian(); }
  if( ok && sens ){ ok = oc.solve_asens( x3.data(), inp.data(), nullptr );  if( ok ) R.Ja = oc.sens_jacobian(); }
  std::cerr.rdbuf( old );  R.msg = os.str();
  R.ok = ok && R.F.size() == 4 && ( !sens || ( R.Jf.size() == 8 && R.Ja.size() == 8 ) );
  return R; }
static std::string reason( std::string const& m ){
  for( char const* key : { "TRANSITION", "FAILED", "not converge", "structural", "**" } ){
    auto const k = m.find( key ); if( k == std::string::npos ) continue;
    auto const b = m.rfind( '\n', k ); auto const e = m.find( '\n', k );
    return m.substr( b == std::string::npos? 0: b+1, std::min<size_t>( 110, ( e == std::string::npos? m.size(): e ) - ( b == std::string::npos? 0: b+1 ) ) ); }
  return "failed (no message)"; }
int main(){
  if( getenv( "PSA_REFINE" ) ){                          // refinement study: marching vs monolithic, with and without the jump
    Out const A0 = run( 0, 0., false, false ), B0 = run( 0, 0., true, false ), A1 = run( 1, DJ, false, false ), B1 = run( 1, DJ, true, false );
    auto gap = [&]( Out const& a, Out const& b ){ double g = 0.; for( size_t j = 0; j < 2; ++j ) g = std::max( g, std::fabs( a.F[j]-b.F[j] )/std::max( 1., std::fabs( a.F[j] ) ) ); return g; };
    if( A0.ok && B0.ok && A1.ok && B1.ok ){
      std::printf( "  n_el=%2s n_nd=%2s  |mono - march|, with the jump:", getenv( "PSA_NEL" )? getenv( "PSA_NEL" ): "5", getenv( "PSA_NND" )? getenv( "PSA_NND" ): "6" );
      char const* nm[4] = { "Inv1", "Eff1", "c1(tau-)", "c1(tau+)" };
      for( size_t j = 0; j < 4; ++j ) std::printf( "  %s %.1e", nm[j], std::fabs( A1.F[j]-B1.F[j] ) );
      std::printf( "   (no jump: Inv1 %.1e, c1(tau) %.1e)\n", std::fabs( A0.F[0]-B0.F[0] ), std::fabs( A0.F[2]-B0.F[2] ) );
      std::printf( "      CAUSALITY  c1(tau-) with the jump minus c1(tau) without:  monolithic %+.2e   marching %+.2e\n", A1.F[2]-A0.F[2], B1.F[2]-B0.F[2] );
      if( getenv( "PSA_VALUES" ) ) std::printf( "      c1(tau-,1/2): mono %.10f march %.10f   c1(tau+,1/2): mono %.10f march %.10f\n", A1.F[2], B1.F[2], A1.F[3], B1.F[3] );
    }
    else std::printf( "  n_el=%s n_nd=%s   a solve failed\n", getenv( "PSA_NEL" ), getenv( "PSA_NND" ) );
    return 0;
  }
  Out J[2];
  for( int mi = 0; mi < 2; ++mi ){
    bool const march = mi == 1;
    std::string const tag = march? "marching   : ": "monolithic : ";
    Out const N = run( 0, 0., march, false ), I0 = run( 1, 0., march, false );
    J[mi] = run( 1, DJ, march );
    if( !N.ok || !I0.ok || !J[mi].ok ){ check( false, tag + "setup and solves -- " + reason( !N.ok? N.msg: !I0.ok? I0.msg: J[mi].msg ) ); continue; }
    Out const& R = J[mi];
    double ei = 0.; for( size_t j = 0; j < 4; ++j ) ei = std::max( ei, std::fabs( I0.F[j]-N.F[j] ) );
    // IDENTITY: exact when marching; monolithic to discretisation level, because an identity transition's lifted claim
    // receives AFTER tau only while the plain model's claim receives on both sides (see CAUSALITY)
    check( ei < ( march? 1e-14: 1e-4 ), tag + ( march? "IDENTITY (d = 0) == no transition (exact)": "IDENTITY (d = 0) == no transition (to discretisation level)" ), ei );
    // CAUSALITY: the value BEFORE the jump must not depend on the jump.  Marching (window by window) is causal exactly.
    // Monolithic is causal up to DISCRETISATION level: ordinary evolution continuity claims receive on both sides of a
    // face, so their multipliers couple past and future (the long-standing monolithic/marching gap); a lifted claim
    // receives after the jump only (2026-09-30).  Measured 9.2e-6 at n_nd = 6, decreasing with refinement; received
    // downstream-only EVERYWHERE, monolithic equals marching to 1e-13 (the receiver-policy study on the list).
    double const ec = std::fabs( R.F[2] - N.F[2] );
    check( ec < ( march? 1e-14: 1e-4 ), tag + ( march? "CAUSALITY: c1(tau-) unaffected by the jump (exact)"
                                                   : "CAUSALITY: c1(tau-) unaffected by the jump (to discretisation level)" ), ec );
    double const ej = std::fabs( R.F[3] - ( 1.-DJ )*R.F[2] );
    check( ej < 1e-10, tag + "the jump: c1(tau+,1/2) == (1-d) c1(tau-,1/2)", ej );
    check( R.F[0] < N.F[0], tag + "physics: purging component 1 lowers the final inventory Inv1" );
    double fa = 0.; for( size_t q = 0; q < 8; ++q ) fa = std::max( fa, std::fabs( R.Jf[q]-R.Ja[q] )/std::max( 1., std::fabs( R.Jf[q] ) ) );
    check( fa < 1e-8, tag + "FORWARD == ADJOINT (b01, d; relative)", fa );
    double const h = 1e-6;
    Out const Pp = run( 1, DJ+h, march, false ), Pm = run( 1, DJ-h, march, false );
    if( Pp.ok && Pm.ok ){
      double ef = 0.; for( size_t j = 0; j < 4; ++j ){ double const D = ( Pp.F[j]-Pm.F[j] )/( 2.*h ); ef = std::max( ef, std::fabs( R.Jf[j*2+1]-D )/std::max( 1., std::fabs( D ) ) ); }
      check( ef < 1e-6, tag + "dF/dd == central finite differences (relative)", ef );
    } else check( false, tag + "finite-difference solves" );
  }
  { Out const A = run( 0, 0., false, false ), B = run( 0, 0., true, false );   // CONTROL: the model WITHOUT a transition
    if( A.ok && B.ok ){ double eb = 0.; for( size_t j = 0; j < 4; ++j ) eb = std::max( eb, std::fabs( A.F[j]-B.F[j] )/std::max( 1., std::fabs( A.F[j] ) ) );
      std::printf( "    [control] WITHOUT a transition: marching vs monolithic values differ by %.2e (relative)\n", eb ); } }
  if( J[0].ok && J[1].ok ){
    double ev = 0., eg = 0.;
    for( size_t j = 0; j < 4; ++j ) ev = std::max( ev, std::fabs( J[0].F[j]-J[1].F[j] )/std::max( 1., std::fabs( J[0].F[j] ) ) );
    for( size_t q = 0; q < 8; ++q ) eg = std::max( eg, std::fabs( J[0].Jf[q]-J[1].Jf[q] )/std::max( 1., std::fabs( J[0].Jf[q] ) ) );
    // at the monolithic formulation's discretisation level (see CAUSALITY above): not round-off
    check( ev < 2e-4, "marching == monolithic: values (relative, discretisation level)", ev );
    check( eg < 2e-3, "marching == monolithic: gradients (relative, discretisation level)", eg );
  }
  std::printf( "\n  TRANS_psa: %d passed, %d failed -- %s\n", npass, nfail, nfail? "FAILURES": "ALL PASS" );
  return nfail? 1: 0;
}
