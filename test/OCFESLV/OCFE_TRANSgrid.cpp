// OCFE_TRANSgrid.cpp -- GATE for WORKPLAN 1.6 (2026-10-06): a transition must lie on an element boundary of the
// evolution direction -- in MONOLITHIC and in MARCHING solves alike.  Monolithic refused a tau off the grid at setup;
// marching accepted it and silently IGNORED it (its one-element windows never meet tau): x(T) came back without the
// dose, error 0.42, no message.  Checked on x' = -k x with a dose d at tau, 4 elements on [0, 1]:
//   tau on an interior boundary (0.25, 0.5, 0.75): sets up in both modes; x(T) = exact (the dose is applied);
//   tau inside an element (0.35) or at the domain end (1.0): REFUSED at setup in both modes.
#include <cstdio>
#include <cmath>
#include "ocfeslv.hpp"
using namespace mc;  typedef OCFESLV::EqnRole Role;  typedef OCFESLV::EqnOptions EO;
static int nfail = 0;
static void check( char const* what, bool ok ){ std::printf( "  %-66s %s\n", what, ok? "PASS": "FAIL" ); nfail += !ok; }
struct R { bool ok = false; double v = NAN; };
static R run( double tau, bool march ){
  double const K = .8, D = .7, T = 1.;
  FFGraph G;  FFPartial OpP;  FFEval OpE;  int const NLB = FFDom::ALL - FFDom::LB;
  FFVar t = G.add_var( "t" ), x = G.add_var( "x" ), d = G.add_var( "d" );
  OCFESLV S( &G );  S.add_domain( t, FFDom( 0., T, 4, FFDom::LGR, 5 ) );  S.set_evolution_domain( t );
  S.add_state( x, {t} );  S.update_ref( x, 1. );  S.add_input( d );  S.update_ref( d, D );
  S.add_equation( OpP( x, t ) + K * x, {t}, {NLB}, EO( Role::INTERIOR ) );
  S.add_equation( x - 1., {t}, {(int)FFDom::LB}, EO( Role::INITIAL ) );
  S.add_transition( x + d, x, t, tau );
  S.add_output( OpE( x, t, T ) );
  S.options.INTERFACE.IMPOSITION = OCFESLV::Options::IC_STRONG;  S.options.SOLVE.MARCHING = march;
  S.options.SOLVE.RES_TOL = 1e-12;  S.options.DISPLAY_LEVEL = 0;
  R r;  r.ok = S.setup();
  if( r.ok ){ std::vector<double> var, inp;  S.init( var, inp );
              auto rep = S.solve( var.data(), inp.data(), nullptr );  if( rep.converged ) r.v = S.val_functions()[0]; }
  return r;
}
int main(){
  std::printf( "OCFE_TRANSgrid -- where a transition may lie, in both solve modes\n" );
  double const K = .8, D = .7, T = 1.;  char buf[120];
  for( double tau : { 0.25, 0.5, 0.75 } ) for( bool march : { false, true } ){
    R r = run( tau, march );  double const ex = ( std::exp( -K * tau ) + D ) * std::exp( -K * ( T - tau ) );
    std::snprintf( buf, sizeof buf, "tau = %.2f (a boundary), %-10s: sets up, x(T) = exact (%.1e)", tau, march? "marching": "monolithic", std::fabs( r.v - ex ) );
    check( buf, r.ok && std::fabs( r.v - ex ) < 1e-6 );
  }
  for( double tau : { 0.35, 1.0 } ) for( bool march : { false, true } ){
    std::printf( "    (the refusal below is expected)\n" );
    R r = run( tau, march );
    std::snprintf( buf, sizeof buf, "tau = %.2f (%s), %-10s: REFUSED at setup", tau, tau < T? "inside an element": "the domain end", march? "marching": "monolithic" );
    check( buf, !r.ok );
  }
  std::printf( "  OCFE_TRANSgrid: %d failed\n", nfail );
  return nfail? 1: 0;
}
