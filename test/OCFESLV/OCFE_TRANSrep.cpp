// OCFE_TRANSrep.cpp -- reproduction of the OCFE_TRANSode failure under the row-level causal design: x' = -k x,
// a DOSE at tau, x(tau+) = x(tau-) + d, and y' = x (continuous).  x(T) = (x0 e^{-k tau} + d) e^{-k (T - tau)};
// y(T) = (x0 - x(tau-))/k + (x(tau+) - x(T))/k.  Monolithic and marching, IC_STRONG / IC_WEAK, with
// INTERFACE.EVOLUTION_CAUSAL as argument 1 (0/1); the transition sits at an element boundary (tau = 0.4, 5 elements).
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include "ocfeslv.hpp"
using namespace mc;  typedef FFModel::EqnRole Role;
int main( int argc, char** argv ){
  bool const causal = argc > 1 && std::atoi( argv[1] );
  double const K = .8, D = .7, X0 = 1., T = 1., TAU = .4;
  double const xm = X0 * std::exp( -K * TAU ), xp = xm + D, xT = xp * std::exp( -K * ( T - TAU ) );
  double const yT = ( X0 - xm ) / K + ( xp - xT ) / K;
  int nfail = 0;
  for( int imp : { 1, 0 } ) for( bool march : { false, true } ){
    FFGraph G;  FFPartial OpP;  FFEval OpE;
    FFVar t = G.add_var( "t" ), x = G.add_var( "x" ), y = G.add_var( "y" ), k = G.add_var( "k" ), d = G.add_var( "d" );
    OCFESLV S( &G );
    S.add_domain( t, FFDom( 0., T, 5, FFDom::LGR, 6 ) );  S.set_evolution_domain( t );
    S.add_state( x, {t} );  S.update_ref( x, X0 );  S.add_state( y, {t} );  S.update_ref( y, 0. );
    S.add_input( k );  S.update_ref( k, K );  S.add_input( d );  S.update_ref( d, D );
    S.add_equation( OpP( x, t ) + k * x, {t}, {FFDom::ALL - FFDom::LB}, OCFESLV::EqnOptions( Role::INTERIOR ) );
    S.add_equation( OpP( y, t ) - x, {t}, {FFDom::ALL - FFDom::LB}, OCFESLV::EqnOptions( Role::INTERIOR ) );
    S.add_equation( x - X0, {t}, {FFDom::LB}, OCFESLV::EqnOptions( Role::INITIAL ) );
    S.add_equation( y, {t}, {FFDom::LB}, OCFESLV::EqnOptions( Role::INITIAL ) );
    S.add_transition( x + d, x, t, TAU );                                  // left(tau-) = right(tau+)
    S.add_output( OpE( x, t, T ) );  S.add_output( OpE( y, t, T ) );
    S.options.INTERFACE.IMPOSITION = imp? OCFESLV::Options::IC_STRONG: OCFESLV::Options::IC_WEAK;
    S.options.INTERFACE.EVOLUTION_CAUSAL = causal;
    S.options.SOLVE.MARCHING = march;  S.options.SOLVE.RES_TOL = 1e-12;  S.options.DISPLAY_LEVEL = 0;
    bool ok = S.setup();  std::vector<double> var, inp;  if( ok ) S.init( var, inp );
    ok = ok && S.solve( var.data(), inp.data() ).converged;
    double const ex = ok? std::fabs( S.val_functions()[0] - xT ): 1., ey = ok? std::fabs( S.val_functions()[1] - yT ): 1.;
    bool const pass = ex < 1e-7 && ey < 1e-7;  nfail += !pass;
    std::printf( "causal=%d %-6s %-10s | converged=%d  |x(T) - exact| %.1e  |y(T) - exact| %.1e  %s\n", (int)causal,
                 imp? "STRONG": "WEAK", march? "marching": "monolithic", (int)ok, ex, ey, pass? "PASS": "FAIL" );
  }
  return nfail? 1: 0;
}
