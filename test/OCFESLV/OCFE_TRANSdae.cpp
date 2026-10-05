// TRANS_dae.cpp -- OCFESLV transitions on an index-1 DAE, in every imposition mode, monolithic and marching:
//   dx/dt = -k x + z,  0 = z - c x,  x(0) = x0;   at tau = 0.4: x(tau+) = x(tau-) + d   (z is ALGEBRAIC)
// Closed form: x = x0 e^{-(k-c)t} before tau, (x(tau-)+d) e^{-(k-c)(t-tau)} after; z = c x.
// Checks: values; ALGEBRAIC CONSISTENCY after the jump -- z(tau+) = c x(tau+), so z jumps by c d without being told;
// forward and adjoint gradients w.r.t. (k, c, d, x0) against the closed form; and a transition on the algebraic z is
// REFUSED with the reason (it would otherwise be silently ineffective).
#include <cmath>
#include <cstdio>
#include <sstream>
#include <iostream>
#include <string>
#include <vector>
#include "ocfeslv.hpp"
using namespace mc;
static int npass = 0, nfail = 0;
static void check( bool c, std::string const& w, double v = -1. ){ std::printf( "  %s  %-76s", c? "PASS": "FAIL", w.c_str() ); if( v >= 0. ) std::printf( " (%.2e)", v ); std::printf( "\n" ); c? ++npass: ++nfail; }
double const TAU = 0.4;
static std::vector<double> closed( std::vector<double> const& P ){          // P = k, c, d, x0 -> x(1), z(tau-), z(tau+), z(1), int z
  double const k = P[0], c = P[1], d = P[2], x0 = P[3], r = k - c;
  double const xm = x0*std::exp( -r*TAU ), xp = xm + d, x1 = xp*std::exp( -r*( 1.-TAU ) );
  double const Iz = c*( x0*( 1.-std::exp( -r*TAU ) )/r + xp*( 1.-std::exp( -r*( 1.-TAU ) ) )/r );
  return { x1, c*xm, c*xp, c*x1, Iz }; }
struct Out { bool ok = false; std::string msg; std::vector<double> F, Jf, Ja; };
static Out run( int imp, bool march, std::vector<double> const& P0, bool on_z = false ){
  Out R; FFGraph G; OCFESLV I( &G ); FFPartial OpP; FFEval OpE; FFIntegral OpI;
  FFVar t = G.add_var( "t" ), x = G.add_var( "x(t)" ), z = G.add_var( "z(t)" );
  char const* pn[4] = { "k", "c", "d", "x0" }; std::vector<FFVar> p; for( auto n : pn ) p.push_back( G.add_var( n ) );
  I.add_domain( t, FFDom( std::vector<double>{ 0., TAU, 0.7, 1. }, FFDom::LGR, 8 ) ); I.set_evolution_domain( t );
  I.add_state( x, {t} ); I.add_state( z, {t} ); for( auto const& v : p ) I.add_input( v, {} );
  I.update_ref( x, 1. ); I.update_ref( z, 0.5 );
  int const T_INT = FFDom::ALL - FFDom::LB;
  I.add_equation( OpP( x, t ) + p[0]*x - z, {t}, {T_INT}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  I.add_equation( z - p[1]*x, {t}, {FFDom::ALL}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  I.add_equation( x - p[3], {t}, {FFDom::LB}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INITIAL, 0 ) );
  if( on_z ) I.add_transition( z + p[2], z, t, TAU ); else I.add_transition( x + p[2], x, t, TAU );
  I.add_output( OpE( x, t, 1. ) ); I.add_output( OpE( z, t, TAU, FFDom::MINUS ) ); I.add_output( OpE( z, t, TAU, FFDom::PLUS ) );
  I.add_output( OpE( z, t, 1. ) ); I.add_output( OpI( std::vector<FFVar>{ z }, { t } )[0] );
  I.options.INTERFACE.IMPOSITION = imp == 0? OCFESLV::Options::IC_WEAK: imp == 1? OCFESLV::Options::IC_STRONG: OCFESLV::Options::IC_TRACE;
  I.options.SOLVE.MARCHING = march; I.options.SOLVE.RES_TOL = 1e-12; I.options.DISPLAY_LEVEL = 0;
  std::ostringstream os; auto* old = std::cerr.rdbuf( os.rdbuf() );
  std::vector<double> xv, inp; bool ok = I.setup() && I.init( xv, inp, nullptr );
  if( ok ) for( size_t q = 0; q < 4; ++q ){ I.set_input_values( p[q], { P0[q] }, inp.data() ); I.register_control( p[q] ); }
  std::vector<double> y1( xv ), y2( xv ), y3( xv );
  if( ok ) ok = I.solve( y1.data(), inp.data(), nullptr ).converged;  if( ok ) R.F = I.val_functions();
  if( ok ) ok = I.solve_fsens( y2.data(), inp.data(), nullptr );     if( ok ) R.Jf = I.sens_jacobian();
  if( ok ) ok = I.solve_asens( y3.data(), inp.data(), nullptr );     if( ok ) R.Ja = I.sens_jacobian();
  std::cerr.rdbuf( old );  R.msg = os.str();
  R.ok = ok && R.F.size() == 5 && R.Jf.size() == 20 && R.Ja.size() == 20;
  return R; }
int main(){
  std::vector<double> const P0{ 1.3, 0.5, 0.7, 1.1 };
  char const* mn[3] = { "IC_WEAK  ", "IC_STRONG", "IC_TRACE " };
  for( bool march : { false, true } ) for( int imp = 0; imp < 3; ++imp ){
    std::string const tag = std::string( march? "marching   ": "monolithic " ) + mn[imp] + ": ";
    Out const R = run( imp, march, P0 );
    if( !R.ok ){ auto const k = R.msg.find( "**" ); check( false, tag + "setup and solves -- " + ( k != std::string::npos? R.msg.substr( k, 70 ): std::string( "failed" ) ) ); continue; }
    auto const F0 = closed( P0 ); double ev = 0.; for( size_t j = 0; j < 5; ++j ) ev = std::max( ev, std::fabs( R.F[j]-F0[j] ) );
    check( ev < 1e-9, tag + "values == closed form", ev );
    // z(tau+) = c x(tau+) = c ( x(tau-) + d ) and z(tau-) = c x(tau-): the algebraic state jumps by exactly c d
    double const jz = R.F[2] - R.F[1];
    check( std::fabs( jz - P0[1]*P0[2] ) < 1e-9, tag + "ALGEBRAIC z re-determined: z(tau+) - z(tau-) == c d", std::fabs( jz - P0[1]*P0[2] ) );
    double ef = 0., ea = 0.;
    for( size_t q = 0; q < 4; ++q ){ auto Pp = P0, Pm = P0; double const h = 1e-6*std::max( 1., std::fabs( P0[q] ) ); Pp[q] += h; Pm[q] -= h;
      auto const fp = closed( Pp ), fm = closed( Pm );
      for( size_t j = 0; j < 5; ++j ){ double const D = ( fp[j]-fm[j] )/( 2.*h ); ef = std::max( ef, std::fabs( R.Jf[j*4+q]-D ) ); ea = std::max( ea, std::fabs( R.Ja[j*4+q]-D ) ); } }
    check( ef < 1e-7, tag + "FORWARD gradients (k, c, d, x0) == closed form", ef );
    check( ea < 1e-7, tag + "ADJOINT gradients (k, c, d, x0) == closed form", ea );
  }
  { Out const R = run( 1, false, P0, true );
    check( !R.ok && R.msg.find( "not differentiated along the evolution direction" ) != std::string::npos,
           "a transition on the ALGEBRAIC z is REFUSED, with the reason" ); }
  std::printf( "\n  TRANS_dae: %d passed, %d failed -- %s\n", npass, nfail, nfail? "FAILURES": "ALL PASS" );
  return nfail? 1: 0;
}
