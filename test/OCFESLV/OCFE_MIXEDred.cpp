// MIXED_red.cpp -- reductions consuming the EVOLUTION direction AND a spatial one (mixed), in outputs and in
// equations, monolithic and marching, against closed forms.  FFModel defers a reduction (captures it after the solve)
// only when it consumes exactly the evolution direction (_is_deferred_reduction); a mixed one stays IN-SOLVE as an
// auxiliary, which every marching window carries but only the window holding the point can satisfy.  The output-point
// site was fixed on 2026-09-30 (O1); this driver checks the others.
// DECISION (2026-09-30): such a reduction may NOT appear in an EQUATION, in either mode (symmetry): refused, saying
// whether its operand depends on a STATE (causality) or on INPUTS only (pre-solve evaluation: an open item).  Before,
// E1 / I1 / I2 returned a silent 0 in both modes (the post-solve latch read during the solve), E3 integrated one window
// when marching, and E2 / I1-I4 failed marching setup.  OUTPUTS keep deferring correctly (O1, O2).
//   Heat equation u_t = a u_zz, u(t,0) = u(t,1) = 0, u(0,z) = sin(pi z)  ->  u = e^{-a pi^2 t} sin(pi z);  a = 1/4
//   Scalar ODE x' = -k x, x(0) = 1  ->  x(1) = e^{-k};  k = 0.7
//   O1  output  OpE( u, {t,z}, {1, 1/2} )             = e^{-a pi^2}                        (output point: the control)
//   O2  output  int int u dz dt                         = (2/pi)(1 - e^{-a pi^2})/(a pi^2)  (output integral)
//   E1  s = OpE( x, t, 1 )           (equation; s domain-less)                             -> REFUSED (STATE)
//   E2  s = OpE( u, {t,z}, {1, 1/2} ) (equation)                                          -> REFUSED (STATE)
//   E3  s = int int u dz dt           (equation)                                            -> REFUSED (STATE)
//   I1-I4  the same with INPUTS p(t), q(t,z): OpE(p,t,1), int p dt, OpE(q,{t,z},{1,1/2}), int int q dz dt -> REFUSED (INPUTS)
#include <cmath>
#include <cstdio>
#include <sstream>
#include <iostream>
#include <string>
#include <vector>
#include "ocfeslv.hpp"
using namespace mc;
static int npass = 0, nfail = 0;
static void check( bool c, std::string const& w, double v = -1. ){ std::printf( "  %s  %-74s", c? "PASS": "FAIL", w.c_str() ); if( v >= 0. ) std::printf( " (%.2e)", v ); std::printf( "\n" ); c? ++npass: ++nfail; }
double const A = 0.25, K = 0.7;
struct Out { bool ok = false; double v = 0.; std::string msg; };
static std::string reason( std::string const& m ){
  for( char const* key : { "REFUSED", "refused", "structural", "not square", "FAILED", "failed", "**" } ){
    auto const k = m.find( key ); if( k == std::string::npos ) continue;
    auto const b = m.rfind( '\n', k ); auto const e = m.find( '\n', k );
    return m.substr( b == std::string::npos? 0: b+1, std::min<size_t>( 120, ( e == std::string::npos? m.size(): e ) - ( b == std::string::npos? 0: b+1 ) ) ); }
  return "no message"; }
static Out run( std::string const& cs, bool march ){
  Out R; FFGraph G; OCFESLV I( &G ); FFPartial OpP; FFEval OpE; FFIntegral OpI;
  FFVar t = G.add_var( "t" ), z = G.add_var( "z" ), u = G.add_var( "u(t,z)" ), x = G.add_var( "x(t)" ), s = G.add_var( "s" );
  FFVar pin = G.add_var( "p(t)" ), qin = G.add_var( "q(t,z)" );
  I.add_domain( t, FFDom( std::vector<double>{ 0., 0.4, 1. }, FFDom::LGR, 8 ) );
  I.add_domain( z, FFDom( 0., 1., 4, FFDom::LGL, 7 ) );  I.set_evolution_domain( t );
  I.add_state( u, {t,z} );  I.add_state( x, {t} );
  bool const inputs = cs[0] == 'I';
  if( inputs ){ I.add_input( pin, {t} );  I.add_input( qin, {t,z} ); }  I.update_ref( u, 0.5 );  I.update_ref( x, 1. );
  int const T_NO_LB = FFDom::ALL - FFDom::LB, Z_INT = FFDom::ALL - FFDom::LB - FFDom::UB;
  I.add_equation( OpP( u, t ) - A*OpP( OpP( u, z ), z ), {t,z}, {T_NO_LB, Z_INT},       OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  I.add_equation( u,                                     {t,z}, {T_NO_LB, FFDom::LB},   OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );
  I.add_equation( u,                                     {t,z}, {T_NO_LB, FFDom::UB},   OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );
  I.add_equation( u - sin( M_PI*z ),                     {t,z}, {FFDom::LB, FFDom::ALL}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INITIAL, 0 ) );
  I.add_equation( OpP( x, t ) + K*x,                     {t},   {T_NO_LB},              OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  I.add_equation( x - 1.,                                {t},   {FFDom::LB},            OCFESLV::EqnOptions( OCFESLV::EqnRole::INITIAL, 0 ) );
  SMon<FFVar,lt_FFVar> const TZ( std::map<FFVar,unsigned,lt_FFVar>{ { t, 1u }, { z, 1u } } );
  FFVar const uPt = OpE( u, TZ, std::map<FFVar,double,lt_FFVar>{ { t, 1. }, { z, 0.5 } } );
  FFVar const uII = OpI( std::vector<FFVar>{ OpI( std::vector<FFVar>{ u }, { z } )[0] }, { t } )[0];   // folds into one {t,z} node
  if( cs == "O1" ) I.add_output( uPt );
  if( cs == "O2" ) I.add_output( uII );
  FFVar const qPt = OpE( qin, TZ, std::map<FFVar,double,lt_FFVar>{ { t, 1. }, { z, 0.5 } } );
  FFVar const qII = OpI( std::vector<FFVar>{ OpI( std::vector<FFVar>{ qin }, { z } )[0] }, { t } )[0];
  if( cs[0] == 'E' || inputs ){
    I.add_state( s, {} );  I.update_ref( s, 0.5 );
    FFVar const rhs = cs == "E1"? OpE( x, t, 1. ): cs == "E2"? uPt: cs == "E3"? uII
                    : cs == "I1"? OpE( pin, t, 1. ): cs == "I2"? OpI( std::vector<FFVar>{ pin }, { t } )[0]: cs == "I3"? qPt: qII;
    I.add_equation( s - rhs, {}, {}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
    I.add_output( s );
  }
  I.options.SOLVE.MARCHING = march;  I.options.SOLVE.RES_TOL = 1e-12;  I.options.DISPLAY_LEVEL = 0;
  std::ostringstream os; auto* old = std::cerr.rdbuf( os.rdbuf() );
  std::vector<double> xv, inp;  bool ok = false;
  try{ ok = I.setup() && I.init( xv, inp, nullptr );
       if( ok && inputs ){ for( auto const& w : { pin, qin } ){ size_t const n = I.get_input_values( w, inp.data() ).size(); I.set_input_values( w, std::vector<double>( n, 2. ), inp.data() ); } }
       ok = ok && I.solve( xv.data(), inp.data(), nullptr ).converged; }
  catch( std::exception& e ){ os << "** exception: " << e.what(); ok = false; }
  std::cerr.rdbuf( old );  R.msg = os.str();
  if( ok ){ auto const F = I.val_functions(); if( F.size() == 1 ){ R.ok = true; R.v = F[0]; } }
  return R; }
int main(){
  double const E1pt = std::exp( -A*M_PI*M_PI ), II = 2./M_PI*( 1. - E1pt )/( A*M_PI*M_PI );
  struct Case { char const* id; char const* what; double ref; };
  Case const cases[] = { { "O1", "output   OpE(u,{t,z},{1,1/2})  (output point: fixed, the control)", E1pt },
                         { "O2", "output   int int u dz dt        (output integral)", II },
                         { "E1", "equation s = OpE(x,t,1)         (single direction: reference)", std::exp( -K ) },
                         { "E2", "equation s = OpE(u,{t,z},{1,1/2}) (mixed point)", E1pt },
                         { "E3", "equation s = int int u dz dt   (mixed integral)", II },
                         { "I1", "equation s = OpE(p,t,1)         (INPUT, single direction; p = 2)", 2. },
                         { "I2", "equation s = int p dt           (INPUT, single direction)", 2. },
                         { "I3", "equation s = OpE(q,{t,z},{1,1/2}) (INPUT, mixed point; q = 2)", 2. },
                         { "I4", "equation s = int int q dz dt   (INPUT, mixed integral)", 2. } };
  for( auto const& c : cases ){
    std::printf( "\n  %s  %s   closed form %.10f\n", c.id, c.what, c.ref );
    for( bool march : { false, true } ){
      std::string const tag = std::string( c.id ) + ( march? " marching  : ": " monolithic: " );
      Out const R = run( c.id, march );
      if( c.id[0] != 'O' ){                             // equations: REFUSED, with the kind of dependence
        bool const st = c.id[0] == 'E';
        bool const refused = !R.ok && R.msg.find( "EQUATION REFUSED" ) != std::string::npos
                             && R.msg.find( st? "depends on a STATE": "of INPUTS only" ) != std::string::npos;
        check( refused, tag + ( st? "REFUSED: depends on a STATE (causality)": "REFUSED: INPUTS only (pre-solve evaluation: open item)" ) );
        continue;
      }
      if( !R.ok ){ check( false, tag + "solves -- " + reason( R.msg ) ); continue; }
      double const e = std::fabs( R.v - c.ref );
      std::printf( "      value %.10f\n", R.v );
      check( e < 1e-6, tag + "== closed form (discretisation level)", e );
    }
  }
  std::printf( "\n  MIXED_red: %d passed, %d failed -- %s\n", npass, nfail, nfail? "FAILURES": "ALL PASS" );
  return nfail? 1: 0;
}
