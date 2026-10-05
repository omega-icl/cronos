// OCFE_TRANSgrad.cpp -- GRADIENTS across a transition (marching and monolithic, forward and adjoint) against the
// closed form: x' = -k x, a DOSE at tau, x(tau+) = x(tau-) + d.  x(T) = (x0 e^{-k tau} + d) e^{-k (T - tau)}:
//   d x(T)/dk = -tau x0 e^{-kT} - (T - tau)(x0 e^{-k tau} + d) e^{-k (T - tau)},   d x(T)/dd = e^{-k (T - tau)}.
#include <cmath>
#include <cstdio>
#include "ocfeslv.hpp"
using namespace mc;  typedef FFModel::EqnRole Role;
int main(){
  double const K = .8, D = .7, X0 = 1., T = 1., TAU = .4;
  double const xT = ( X0 * std::exp( -K * TAU ) + D ) * std::exp( -K * ( T - TAU ) );
  double const gk = -TAU * X0 * std::exp( -K * T ) - ( T - TAU ) * ( X0 * std::exp( -K * TAU ) + D ) * std::exp( -K * ( T - TAU ) );
  double const gd = std::exp( -K * ( T - TAU ) );
  int nfail = 0;
  for( bool march : { false, true } ) for( int adj : { 0, 1 } ){
    FFGraph G;  FFPartial OpP;  FFEval OpE;
    FFVar t = G.add_var( "t" ), x = G.add_var( "x" ), k = G.add_var( "k" ), d = G.add_var( "d" );
    OCFESLV S( &G );
    S.add_domain( t, FFDom( 0., T, 5, FFDom::LGR, 6 ) );  S.set_evolution_domain( t );
    S.add_state( x, {t} );  S.update_ref( x, X0 );
    S.add_input( k );  S.update_ref( k, K );  S.add_input( d );  S.update_ref( d, D );
    S.add_equation( OpP( x, t ) + k * x, {t}, {FFDom::ALL - FFDom::LB}, OCFESLV::EqnOptions( Role::INTERIOR ) );
    S.add_equation( x - X0, {t}, {FFDom::LB}, OCFESLV::EqnOptions( Role::INITIAL ) );
    S.add_transition( x + d, x, t, TAU );
    S.add_output( OpE( x, t, T ) );
    S.register_control( k );  S.register_control( d );
    S.options.INTERFACE.IMPOSITION = OCFESLV::Options::IC_STRONG;
    S.options.SOLVE.MARCHING = march;  S.options.SOLVE.RES_TOL = 1e-12;  S.options.DISPLAY_LEVEL = 0;
    bool ok = S.setup();  std::vector<double> var, inp;  if( ok ) S.init( var, inp );
    ok = ok && ( adj? S.solve_asens( var.data(), inp.data(), nullptr ): S.solve_fsens( var.data(), inp.data(), nullptr ) );
    std::vector<double> const& J = S.sens_jacobian();
    double const ev = ok? std::fabs( S.sens_functions()[0] - xT ): 1.;
    double const ek = ok && J.size() >= 2? std::fabs( J[0] - gk ): 1., ed = ok && J.size() >= 2? std::fabs( J[1] - gd ): 1.;
    bool const pass = ev < 1e-7 && ek < 1e-6 && ed < 1e-6;  nfail += !pass;
    std::printf( "%-10s %-7s | x(T) %.1e | dx/dk %+.8f (exact %+.8f) | dx/dd %+.8f (exact %+.8f)  %s\n",
                 march? "marching": "monolithic", adj? "adjoint": "forward", ev, ok? J[0]: 0., gk, ok? J[1]: 0., gd, pass? "PASS": "FAIL" );
  }
  return nfail? 1: 0;
}
