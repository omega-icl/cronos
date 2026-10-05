// TRANS_oc.cpp -- OCFESLV transitions (add_transition), in EVERY imposition mode (IC_WEAK / IC_STRONG / IC_TRACE),
// monolithic and marching (a transfer map at the window seam), against CLOSED FORMS.  Mesh {0, 0.3, 0.5, 1}: both transition times are element boundaries.
//   dx1/dt = -a x1, dx2/dt = -b x2, x1(0) = x10, x2(0) = x20
//   I  IDENTITY transition x1(0.3+) = x1(0.3-)  ==  the model WITHOUT a transition (isolates the claim replacement)
//   E  EXPLICIT: x1(0.3+) = x1(0.3-) + d ;  x2(0.5+) = x2(0.5-) + c x1(0.5-)^2      (the TRANS_multi model)
//   M  IMPLICIT (OCFESLV only): x1+ + x2+ = x1- + x2- + d ,  x1+ - x2+ = x1- - x2-  at 0.3  =>  both jump by d/2
// Values and FORWARD / ADJOINT gradients w.r.t. (a, b, c, d, x10, x20) against central differences of the closed form.
#include <cmath>
#include <cstdio>
#include <sstream>
#include <iostream>
#include <string>
#include <vector>
#include "ocfeslv.hpp"
using namespace mc;
static int npass = 0, nfail = 0;
//! the reason a case did not run: the refusal line if there is one, else the first informative line
static std::string reason( std::string const& msg ){
  for( char const* key : { "REFUSED", "no solver supports", "not supported", "** " } ){
    auto const k = msg.find( key ); if( k == std::string::npos ) continue;
    auto const b = msg.rfind( '\n', k ); auto const e = msg.find( '\n', k );
    return msg.substr( b == std::string::npos? 0: b+1, std::min<size_t>( 90, ( e == std::string::npos? msg.size(): e ) - ( b == std::string::npos? 0: b+1 ) ) ); }
  return "failed (no message)"; }
static void check( bool c, std::string const& w, double v = -1. ){ std::printf( "  %s  %-74s", c? "PASS": "FAIL", w.c_str() ); if( v >= 0. ) std::printf( " (%.2e)", v ); std::printf( "\n" ); c? ++npass: ++nfail; }
double const T1 = 0.3, T2 = 0.5;
static std::vector<double> closed( int kind, std::vector<double> const& P ){   // outputs: x1(1), x2(1), int (x1+x2)
  double const a = P[0], b = P[1], c = P[2], d = P[3], x10 = P[4], x20 = P[5];
  double j1 = 0., j2 = 0.;  double tj1 = T1, tj2 = T2;                         // jumps and their times
  auto run = [&]( double x0, double k, double tj, double jmp ){               // x(1) and int_0^1 x, one jump at tj
    double const xm = x0*std::exp( -k*tj ), xp = xm + jmp;
    return std::pair<double,double>{ xp*std::exp( -k*( 1.-tj ) ), x0*( 1.-std::exp( -k*tj ) )/k + xp*( 1.-std::exp( -k*( 1.-tj ) ) )/k }; };
  if( kind == 1 ){ j1 = d; double const x1m = x10*std::exp( -a*T1 ) + d; double const x1_T2 = x1m*std::exp( -a*( T2-T1 ) ); j2 = c*x1_T2*x1_T2; }
  if( kind == 2 ){ j1 = d/2.; j2 = d/2.; tj2 = T1; }
  auto const r1 = run( x10, a, tj1, j1 ), r2 = run( x20, b, tj2, j2 );
  return { r1.first, r2.first, r1.second + r2.second }; }
struct Out { bool ok = false; std::string msg; std::vector<double> F, Jf, Ja; };
//! kind: 0 none, 1 explicit (E), 2 implicit (M), 3 identity (I)
static Out run( int kind, int imp, bool march, std::vector<double> const& P0 ){
  Out R; FFGraph G; OCFESLV I( &G ); FFPartial OpP; FFEval OpE; FFIntegral OpI;
  FFVar t = G.add_var( "t" ), x1 = G.add_var( "x1(t)" ), x2 = G.add_var( "x2(t)" );
  char const* pn[6] = { "a", "b", "c", "d", "x10", "x20" }; std::vector<FFVar> p; for( auto n : pn ) p.push_back( G.add_var( n ) );
  I.add_domain( t, FFDom( std::vector<double>{ 0., T1, T2, 1. }, FFDom::LGR, 8 ) ); I.set_evolution_domain( t );
  I.add_state( x1, {t} ); I.add_state( x2, {t} ); for( auto const& v : p ) I.add_input( v, {} );
  I.update_ref( x1, 1. ); I.update_ref( x2, 0.5 );
  int const T_INT = FFDom::ALL - FFDom::LB;
  I.add_equation( OpP( x1, t ) + p[0]*x1, {t}, {T_INT}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  I.add_equation( OpP( x2, t ) + p[1]*x2, {t}, {T_INT}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  I.add_equation( x1 - p[4], {t}, {FFDom::LB}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INITIAL, 0 ) );
  I.add_equation( x2 - p[5], {t}, {FFDom::LB}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INITIAL, 0 ) );
  if( kind == 1 ){ I.add_transition( x1 + p[3], x1, t, T1 );  I.add_transition( x2 + p[2]*x1*x1, x2, t, T2 ); }
  if( kind == 2 )  I.add_transition( std::vector<FFVar>{ x1 + x2 + p[3], x1 - x2 }, std::vector<FFVar>{ x1 + x2, x1 - x2 }, t, T1 );
  if( kind == 3 )  I.add_transition( x1, x1, t, T1 );
  I.add_output( OpE( x1, t, 1. ) ); I.add_output( OpE( x2, t, 1. ) ); I.add_output( OpI( std::vector<FFVar>{ x1 + x2 }, { t } )[0] );
  I.options.INTERFACE.IMPOSITION = imp == 0? OCFESLV::Options::IC_WEAK: imp == 1? OCFESLV::Options::IC_STRONG: OCFESLV::Options::IC_TRACE;
  I.options.SOLVE.MARCHING = march; I.options.SOLVE.RES_TOL = 1e-12; I.options.DISPLAY_LEVEL = 0;
  std::ostringstream os; auto* old = std::cerr.rdbuf( os.rdbuf() );
  bool ok = I.setup();
  std::vector<double> xv, inp;
  if( ok ) ok = I.init( xv, inp, nullptr );
  std::cerr.rdbuf( old ); R.msg = os.str();
  if( !ok ) return R;
  for( size_t q = 0; q < 6; ++q ){ I.set_input_values( p[q], { P0[q] }, inp.data() ); I.register_control( p[q] ); }
  std::vector<double> x1v( xv ), x2v( xv ), x3v( xv );
  old = std::cerr.rdbuf( os.rdbuf() );
  ok = I.solve( x1v.data(), inp.data(), nullptr ).converged;  if( ok ) R.F = I.val_functions();
  if( ok ) ok = I.solve_fsens( x2v.data(), inp.data(), nullptr );  if( ok ) R.Jf = I.sens_jacobian();
  if( ok ) ok = I.solve_asens( x3v.data(), inp.data(), nullptr );  if( ok ) R.Ja = I.sens_jacobian();
  std::cerr.rdbuf( old ); R.msg = os.str();
  R.ok = ok && R.F.size() == 3 && R.Jf.size() == 18 && R.Ja.size() == 18;
  return R; }
int main(){
  std::vector<double> const P0{ 1.3, 0.8, 0.6, 0.7, 1.1, 0.4 };
  char const* mn[3] = { "IC_WEAK  ", "IC_STRONG", "IC_TRACE " };
  for( bool march : { false, true } ) for( int imp = 0; imp < 3; ++imp ){
    std::string const tag = std::string( march? "marching   ": "monolithic " ) + mn[imp] + ": ";
    Out const N = run( 0, imp, march, P0 ), Id = run( 3, imp, march, P0 );
    check( N.ok, tag + "baseline WITHOUT a transition solves" + ( N.ok? std::string(): " -- " + reason( N.msg ) ) );
    if( !Id.ok ) check( false, tag + "I  identity transition -- " + reason( Id.msg ) );
    else{ double w = 0.; for( size_t j = 0; j < 3; ++j ) w = std::max( w, std::fabs( Id.F[j]-N.F[j] ) );
          check( N.ok && w < 1e-10, tag + "I  identity transition == no transition (values)", w ); }
    for( int kind : { 1, 2 } ){
      Out const R = run( kind, imp, march, P0 ); std::string const k = kind == 1? "E  explicit": "M  implicit";
      if( !R.ok ){ check( false, tag + k + " -- " + reason( R.msg ) ); continue; }
      auto const F0 = closed( kind, P0 ); double ev = 0.; for( size_t j = 0; j < 3; ++j ) ev = std::max( ev, std::fabs( R.F[j]-F0[j] ) );
      if( getenv( "TRANS_OC_VERBOSE" ) ) for( size_t j = 0; j < 3; ++j ) std::printf( "      output %zu: computed %.10f  closed %.10f\n", j, R.F[j], F0[j] );
      check( ev < 1e-9, tag + k + ": values == closed form", ev );
      double ef = 0., ea = 0.;
      for( size_t q = 0; q < 6; ++q ){ auto Pp = P0, Pm = P0; double const h = 1e-6*std::max( 1., std::fabs( P0[q] ) ); Pp[q] += h; Pm[q] -= h;
        auto const fp = closed( kind, Pp ), fm = closed( kind, Pm );
        for( size_t j = 0; j < 3; ++j ){ double const D = ( fp[j]-fm[j] )/( 2.*h ); ef = std::max( ef, std::fabs( R.Jf[j*6+q]-D ) ); ea = std::max( ea, std::fabs( R.Ja[j*6+q]-D ) ); } }
      check( ef < 1e-7, tag + k + ": FORWARD gradients (6 parameters) == closed form", ef );
      check( ea < 1e-7, tag + k + ": ADJOINT gradients (6 parameters) == closed form", ea );
    }
  }
  std::printf( "\n  TRANS_oc: %d passed, %d failed -- %s\n", npass, nfail, nfail? "FAILURES": "ALL PASS" );
  return nfail? 1: 0;
}
