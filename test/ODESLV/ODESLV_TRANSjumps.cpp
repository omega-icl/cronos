// TRANS_jumps.cpp -- the jump analysis of FFModel (add_transition): per tau, the post-jump states S+ must be as many as
// the components (SQUARE) and matchable to them (STRUCTURALLY determined).  Each case must get the right verdict AND
// the right reason.  ODESLV additionally requires EXPLICIT maps -- an implicit square map passes the jump analysis and
// is then refused by ODESLV, not by the analysis (OCFESLV will accept it).
#include <cstdio>
#include <sstream>
#include <iostream>
#include <string>
#include "ffode.hpp"
using namespace mc;
static int npass = 0, nfail = 0;
static void check( bool c, std::string const& w ){ std::printf( "  %s  %s\n", c? "PASS": "FAIL", w.c_str() ); c? ++npass: ++nfail; }
struct R { bool ok; std::string msg; };
static R run( int kind ){
  FFGraph G; ODESLVS_CVODES I( &G ); FFPartial OpP; FFEval OpE;
  FFVar t = G.add_var( "t" ), x1 = G.add_var( "x1(t)" ), x2 = G.add_var( "x2(t)" ), x3 = G.add_var( "x3(t)" ), d = G.add_var( "d" );
  I.add_domain( t, FFDom( std::vector<double>{ 0., 0.5, 1. }, FFDom::LGR, 3 ) ); I.set_evolution_domain( t );
  for( auto const& x : { x1, x2, x3 } ){ I.add_state( x, {t} ); I.update_ref( x, 1. ); }
  I.add_input( d );
  int const T_INT = FFDom::ALL - FFDom::LB;
  for( auto const& x : { x1, x2, x3 } ){
    I.add_equation( OpP( x, t ) + x, {t}, {T_INT}, FFModel::EqnOptions( FFModel::EqnRole::INTERIOR, 0 ) );
    I.add_equation( x - 1., {t}, {FFDom::LB}, FFModel::EqnOptions( FFModel::EqnRole::INITIAL, 0 ) ); }
  switch( kind ){
    case 1: I.add_transition( x1 + d, x1, t, 0.3 );  I.add_transition( x2 + d, x2, t, 0.7 );  break;
    case 2: I.add_transition( x1 + d, x1, t, 0.5 );  I.add_transition( x2 + d, x2, t, 0.5 );  break;
    case 3: I.add_transition( std::vector<FFVar>{ x1 + d, x2 }, std::vector<FFVar>{ x1 + x2, x1 - x2 }, t, 0.5 );  break;
    case 4: I.add_transition( x1 + d, x1 + x2, t, 0.5 );  break;
    case 5: I.add_transition( x1 + d, x1, t, 0.5 );  I.add_transition( x1 - d, x1, t, 0.5 );  break;
    case 6: I.add_transition( std::vector<FFVar>{ x1, x1 + d, x2 + x3 }, std::vector<FFVar>{ x1, x1*x1, x2 + x3 }, t, 0.5 );  break;
  }
  I.add_output( OpE( x1 + x2 + x3, t, 1. ) );
  I.options.DISPLAY = 0;
  std::ostringstream os; auto* old = std::cerr.rdbuf( os.rdbuf() ); bool const ok = I.setup(); std::cerr.rdbuf( old );
  return { ok, os.str() + I.extract_error() }; }
int main(){
  { R r = run( 1 ); check( r.ok, "1  two explicit transitions at different tau: ACCEPTED" ); }
  { R r = run( 2 ); check( r.ok, "2  two transitions at the SAME tau, jointly square (x1, x2): ACCEPTED" ); }
  { R r = run( 3 ); check( !r.ok && r.msg.find( "EXPLICIT" ) != std::string::npos && r.msg.find( "post-jump" ) == std::string::npos,
                           "3  implicit square map: passes the jump analysis, refused by ODESLV as not EXPLICIT" ); }
  { R r = run( 4 ); check( !r.ok && r.msg.find( "1 component(s) but 2 post-jump state(s)" ) != std::string::npos,
                           "4  one component, two post-jump states: REFUSED, not square" ); }
  { R r = run( 5 ); check( !r.ok && r.msg.find( "2 component(s) but 1 post-jump state(s)" ) != std::string::npos,
                           "5  one state mapped by two transitions at one tau: REFUSED, not square" ); }
  { R r = run( 6 ); check( !r.ok && r.msg.find( "STRUCTURALLY singular" ) != std::string::npos,
                           "6  three components over three states, two involving only x1: REFUSED, structurally singular" ); }
  std::printf( "\n  TRANS_jumps: %d passed, %d failed -- %s\n", npass, nfail, nfail? "FAILURES": "ALL PASS" );
  return nfail? 1: 0;
}
