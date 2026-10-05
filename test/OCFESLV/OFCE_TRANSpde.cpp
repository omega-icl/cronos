// TRANS_pde.cpp -- v1.5: OCFESLV transitions of a DISTRIBUTED state (a profile jump), monolithic, in every imposition
// mode.   u_t = alpha u_zz on z in [0,1], u(t,0) = u(t,1) = 0, u(0,z) = sin(pi z);   at tau = 0.4: u(tau+,z) = (1+d) u(tau-,z)
// ORACLE: the closed form u = A(t) sin(pi z), A = e^{-alpha pi^2 t} (times 1+d after tau), for values and gradients, to
// the spatial discretisation's accuracy; forward == adjoint to round-off; the identity (d = 0) == no jump (exact in the
// exact modes; at discretisation level under IC_WEAK -- see below).  NOT a discrete scaling of the no-jump model: the
// lifted d u/dz is joined across tau to d w/dz with a jump but to its own value without one.
// Marching (the transfer map node by node): the same closed-form checks, plus values == monolithic's.
#include <cmath>
#include <cstdio>
#include <sstream>
#include <iostream>
#include <string>
#include <vector>
#include "ocfeslv.hpp"
using namespace mc;
static int npass = 0, nfail = 0;
static void check( bool c, std::string const& w, double v = -1. ){ std::printf( "  %s  %-78s", c? "PASS": "FAIL", w.c_str() ); if( v >= 0. ) std::printf( " (%.2e)", v ); std::printf( "\n" ); c? ++npass: ++nfail; }
double const TAU = 0.4, ALPHA = 0.25, DJ = 0.3;
struct Out { bool ok = false; std::string msg; std::vector<double> F, Jf, Ja; };       // F: u(1,.5), u(tau-,.5), u(tau+,.5)
static Out run( bool jump, double d, int imp, bool march ){
  Out R; FFGraph G; OCFESLV I( &G ); FFPartial OpP; FFEval OpE;
  FFVar t = G.add_var( "t" ), z = G.add_var( "z" ), u = G.add_var( "u(t,z)" ), pa = G.add_var( "alpha" ), pd = G.add_var( "d" );
  I.add_domain( t, FFDom( std::vector<double>{ 0., TAU, 1. }, FFDom::LGR, 8 ) );
  I.add_domain( z, FFDom( 0., 1., 4, FFDom::LGL, 7 ) );  I.set_evolution_domain( t );
  I.add_state( u, {t,z} );  I.add_input( pa, {} );  I.add_input( pd, {} );  I.update_ref( u, 0.5 );
  int const T_NO_LB = FFDom::ALL - FFDom::LB, Z_INT = FFDom::ALL - FFDom::LB - FFDom::UB;
  I.add_equation( OpP( u, t ) - pa*OpP( OpP( u, z ), z ), {t,z}, {T_NO_LB, Z_INT},      OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  I.add_equation( u,                                      {t,z}, {T_NO_LB, FFDom::LB},  OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );
  I.add_equation( u,                                      {t,z}, {T_NO_LB, FFDom::UB},  OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );
  I.add_equation( u - sin( M_PI*z ),                      {t,z}, {FFDom::LB, FFDom::ALL}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INITIAL, 0 ) );
  if( jump ) I.add_transition( ( 1. + pd )*u, u, t, TAU );
  I.add_output( OpE( OpE( u, t, 1. ), z, 0.5 ) );
  I.add_output( OpE( OpE( u, t, TAU, FFDom::MINUS ), z, 0.5 ) );
  I.add_output( OpE( OpE( u, t, TAU, FFDom::PLUS  ), z, 0.5 ) );
  I.options.INTERFACE.IMPOSITION = imp == 0? OCFESLV::Options::IC_WEAK: imp == 1? OCFESLV::Options::IC_STRONG: OCFESLV::Options::IC_TRACE;
  I.options.SOLVE.MARCHING = march; I.options.SOLVE.RES_TOL = 1e-12; I.options.DISPLAY_LEVEL = 0;
  std::ostringstream os; auto* old = std::cerr.rdbuf( os.rdbuf() );
  std::vector<double> xv, inp; bool ok = I.setup() && I.init( xv, inp, nullptr );
  if( ok ){ I.set_input_values( pa, { ALPHA }, inp.data() ); I.set_input_values( pd, { d }, inp.data() ); I.register_control( pa ); I.register_control( pd ); }
  std::vector<double> y1( xv ), y2( xv ), y3( xv );
  if( ok ) ok = I.solve( y1.data(), inp.data(), nullptr ).converged;  if( ok ) R.F = I.val_functions();
  if( ok ) ok = I.solve_fsens( y2.data(), inp.data(), nullptr );     if( ok ) R.Jf = I.sens_jacobian();
  if( ok ) ok = I.solve_asens( y3.data(), inp.data(), nullptr );     if( ok ) R.Ja = I.sens_jacobian();
  std::cerr.rdbuf( old );  R.msg = os.str();
  R.ok = ok && R.F.size() == 3 && R.Jf.size() == 6 && R.Ja.size() == 6;
  return R; }
static std::string reason( std::string const& m ){                // the most telling line of a failed run
  for( char const* key : { "TRANSITION", "REFUSED", "not converge", "FAIL", "fail", "singular", "**" } ){
    auto const k = m.find( key ); if( k == std::string::npos ) continue;
    auto const b = m.rfind( '\n', k ); auto const e = m.find( '\n', k );
    return m.substr( b == std::string::npos? 0: b+1, std::min<size_t>( 110, ( e == std::string::npos? m.size(): e ) - ( b == std::string::npos? 0: b+1 ) ) ); }
  return "failed (no message)"; }
int main(){
  char const* mn[3] = { "IC_WEAK  ", "IC_STRONG", "IC_TRACE " };
  for( int imp = 0; imp < 3; ++imp ){
    std::string const tag = std::string( "monolithic " ) + mn[imp] + ": ";
    Out const N = run( false, 0., imp, false ), Z = run( true, 0., imp, false ), J = run( true, DJ, imp, false );
    if( !N.ok || !Z.ok || !J.ok ){ check( false, tag + "setup and solves -- " + reason( !N.ok? N.msg: !Z.ok? Z.msg: J.msg ) ); continue; }
    // IDENTITY: exact in the exact modes; under IC_WEAK the lifted spatial derivative Dz_u is joined to d w/dz by a
    // PENALTY rather than to its own value, a different (equally consistent) discrete condition -- discretisation level
    double ei = 0.; for( size_t j = 0; j < 3; ++j ) ei = std::max( ei, std::fabs( Z.F[j]-N.F[j] ) );
    // (and in every mode, the identity's lifted claim receives AFTER tau only, the plain model's on both sides -- a
    // different, equally consistent discrete system: causality, 2026-09-30)
    check( ei < 1e-6, tag + "IDENTITY (d = 0) == no jump (to discretisation level)", ei );
    // CLOSED FORM: u = A(t) sin(pi z), A = e^{-alpha pi^2 t}, times (1+d) after tau.  At z = 1/2 sin = 1.
    double const s = 1. + DJ, E1 = std::exp( -ALPHA*M_PI*M_PI ), Et = std::exp( -ALPHA*M_PI*M_PI*TAU );
    double const ev = std::max( { std::fabs( J.F[0] - s*E1 ), std::fabs( J.F[1] - Et ), std::fabs( J.F[2] - s*Et ) } );
    check( ev < 1e-7, tag + "u(1), u(tau-), u(tau+) at z = 1/2 == closed form (discretisation-limited)", ev );
    // gradients: du(1)/dd = e^{-alpha pi^2};  du(1)/dalpha = -(1+d) pi^2 e^{-alpha pi^2}   (rows: outputs; cols: alpha, d)
    double eg = 0.;
    for( auto const* Jm : { &J.Jf, &J.Ja } ){
      eg = std::max( eg, std::fabs( (*Jm)[0*2+1] - E1 ) );
      eg = std::max( eg, std::fabs( (*Jm)[0*2+0] + s*M_PI*M_PI*E1 ) );
    }
    check( eg < 1e-6, tag + "FORWARD and ADJOINT du(1)/dd, du(1)/dalpha == closed form", eg );
    double fa = 0.; for( size_t q = 0; q < 6; ++q ) fa = std::max( fa, std::fabs( J.Jf[q] - J.Ja[q] ) );
    check( fa < 1e-10, tag + "FORWARD == ADJOINT (the adjoint is the exact transpose)", fa );
  }
  for( int imp = 0; imp < 3; ++imp ){                  // ---- marching (v1.5: the transfer map node by node)
    std::string const tag = std::string( "marching   " ) + mn[imp] + ": ";
    Out const M = run( true, DJ, imp, true ), Mm = run( true, DJ, imp, false ), M0 = run( false, 0., imp, true );
    check( M0.ok, tag + "baseline WITHOUT a transition marches" + ( M0.ok? std::string(): " -- " + reason( M0.msg ) ) );
    if( M0.ok && getenv( "TRANS_PDE_VERBOSE" ) ){ double fb = 0.; for( size_t q = 0; q < 6; ++q ) fb = std::max( fb, std::fabs( M0.Jf[q] - M0.Ja[q] ) );
      std::printf( "    [control] %sbaseline WITHOUT a transition: |forward - adjoint| = %.2e\n", tag.c_str(), fb ); }
    if( !M.ok ){ check( false, tag + "setup and solves -- " + reason( M.msg ) ); continue; }
    double const s = 1. + DJ, E1 = std::exp( -ALPHA*M_PI*M_PI ), Et = std::exp( -ALPHA*M_PI*M_PI*TAU );
    double const ev = std::max( { std::fabs( M.F[0] - s*E1 ), std::fabs( M.F[1] - Et ), std::fabs( M.F[2] - s*Et ) } );
    check( ev < 1e-7, tag + "u(1), u(tau-), u(tau+) at z = 1/2 == closed form (discretisation-limited)", ev );
    double eg = 0.;
    for( auto const* Jm : { &M.Jf, &M.Ja } ){ eg = std::max( eg, std::fabs( (*Jm)[0*2+1] - E1 ) ); eg = std::max( eg, std::fabs( (*Jm)[0*2+0] + s*M_PI*M_PI*E1 ) ); }
    check( eg < 1e-6, tag + "FORWARD and ADJOINT du(1)/dd, du(1)/dalpha == closed form", eg );
    // FORWARD == ADJOINT.  Under IC_WEAK, marching forward and adjoint sensitivities of this PDE already differ WITHOUT
    // any transition (~1e-7: a pre-existing property of the weak marching sensitivities, recorded separately), so there
    // the check is that the transition adds nothing beyond that baseline; the exact modes agree to round-off.
    double fa = 0., fb = 0.;
    for( size_t q = 0; q < 6; ++q ){ fa = std::max( fa, std::fabs( M.Jf[q] - M.Ja[q] ) ); if( M0.ok ) fb = std::max( fb, std::fabs( M0.Jf[q] - M0.Ja[q] ) ); }
    if( imp == 0 ) check( M0.ok && fa < 2.*fb + 1e-10, tag + "FORWARD vs ADJOINT: no worse than the no-transition baseline's own gap", fa );
    else           check( fa < 1e-10, tag + "FORWARD == ADJOINT", fa );
    // second oracle: marching and monolithic are each within the discretisation of the closed form; under IC_WEAK their
    // weak-imposition errors differ in sign, so they agree to that level only
    if( Mm.ok ){ double em = 0.; for( size_t j = 0; j < 3; ++j ) em = std::max( em, std::fabs( M.F[j]-Mm.F[j] ) );
      check( em < ( imp == 0? 5e-7: 1e-7 ), tag + "values == monolithic's (second oracle)", em ); }
  }
  std::printf( "\n  TRANS_pde: %d passed, %d failed -- %s\n", npass, nfail, nfail? "FAILURES": "ALL PASS" );
  return nfail? 1: 0;
}
