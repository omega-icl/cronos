// TRANS_multi.cpp -- ODESLV STATE DISCONTINUITIES (add_transition), against a CLOSED FORM, values + FORWARD and
// ADJOINT gradients w.r.t. every parameter (central differences of the closed form).
//   dx1/dt = -a x1,  dx2/dt = -b x2,   x1(0) = x10, x2(0) = x20,   mesh {0, 0.5, 1}
//   tau1 = 0.3 (INSIDE an element):       x1(tau1+) = x1(tau1-) + d               (x2 not mentioned: continuous)
//   tau2 = 0.5 (an ELEMENT BOUNDARY):     x2(tau2+) = x2(tau2-) + c x1(tau2-)^2   (nonlinear, cross-coupled map)
//   outputs: x1(1), x2(1), x2(tau2-), x2(tau2+), int (x1 + x2) dt
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>
#include "ffode.hpp"
using namespace mc;
static int npass = 0, nfail = 0;
static void check( bool c, std::string const& w, double v ){ std::printf( "  %s  %-62s (%.2e)\n", c? "PASS": "FAIL", w.c_str(), v ); c? ++npass: ++nfail; }
double const T1 = 0.3, T2 = 0.5;
static std::vector<double> closed( std::vector<double> const& P ){          // P = a, b, c, d, x10, x20
  double const a = P[0], b = P[1], c = P[2], d = P[3], x10 = P[4], x20 = P[5];
  double const x1m = x10*std::exp( -a*T1 ), x1p = x1m + d;                   // at tau1
  auto x1 = [&]( double t ){ return t < T1? x10*std::exp( -a*t ): x1p*std::exp( -a*( t-T1 ) ); };
  double const x2m = x20*std::exp( -b*T2 ), x2p = x2m + c*x1( T2 )*x1( T2 );
  double const I1 = x10*( 1.-std::exp( -a*T1 ) )/a + x1p*( 1.-std::exp( -a*( 1.-T1 ) ) )/a;
  double const I2 = x20*( 1.-std::exp( -b*T2 ) )/b + x2p*( 1.-std::exp( -b*( 1.-T2 ) ) )/b;
  return { x1( 1. ), x2p*std::exp( -b*( 1.-T2 ) ), x2m, x2p, I1 + I2 }; }
int main(){
  std::vector<double> const P0{ 1.3, 0.8, 0.6, 0.7, 1.1, 0.4 };
  char const* nm[5] = { "x1(1)", "x2(1)", "x2(tau2-)", "x2(tau2+)", "int (x1+x2) dt" };
  char const* pn[6] = { "a", "b", "c", "d", "x10", "x20" };
  FFGraph G; ODESLVS_CVODES I( &G ); FFPartial OpP; FFEval OpE; FFIntegral OpI;
  FFVar t = G.add_var( "t" ), x1 = G.add_var( "x1(t)" ), x2 = G.add_var( "x2(t)" );
  std::vector<FFVar> p; for( auto n : pn ) p.push_back( G.add_var( n ) );
  I.add_domain( t, FFDom( std::vector<double>{ 0., 0.5, 1. }, FFDom::LGR, 4 ) ); I.set_evolution_domain( t );
  I.add_state( x1, {t} ); I.add_state( x2, {t} ); for( auto const& v : p ) I.add_input( v );
  I.update_ref( x1, 1. ); I.update_ref( x2, 0.5 );
  int const T_INT = FFDom::ALL - FFDom::LB;
  I.add_equation( OpP( x1, t ) + p[0]*x1, {t}, {T_INT}, FFModel::EqnOptions( FFModel::EqnRole::INTERIOR, 0 ) );
  I.add_equation( OpP( x2, t ) + p[1]*x2, {t}, {T_INT}, FFModel::EqnOptions( FFModel::EqnRole::INTERIOR, 0 ) );
  I.add_equation( x1 - p[4], {t}, {FFDom::LB}, FFModel::EqnOptions( FFModel::EqnRole::INITIAL, 0 ) );
  I.add_equation( x2 - p[5], {t}, {FFDom::LB}, FFModel::EqnOptions( FFModel::EqnRole::INITIAL, 0 ) );
  I.add_transition( x1 + p[3], x1, t, T1 );
  I.add_transition( x2 + p[2]*x1*x1, x2, t, T2 );
  I.add_output( OpE( x1, t, 1. ) ); I.add_output( OpE( x2, t, 1. ) );
  I.add_output( OpE( x2, t, T2, FFDom::MINUS ) ); I.add_output( OpE( x2, t, T2, FFDom::PLUS ) );
  I.add_output( OpI( std::vector<FFVar>{ x1 + x2 }, { t } )[0] );
  I.options.DISPLAY = 0; I.options.ATOL = I.options.ATOLS = I.options.ATOLB = 1e-12; I.options.RTOL = I.options.RTOLS = I.options.RTOLB = 1e-12;
  bool ok = I.setup(); check( ok, "setup (two transitions)", 0. ); if( !ok ) return 1;
  std::vector<double> Q( I.np() ); std::vector<size_t> ix; for( size_t q = 0; q < 6; ++q ){ ix.push_back( I.parameter_index( p[q] )[0] ); Q[ix[q]] = P0[q]; }
  bool const of = I.solve_fsens( Q ) == ODESLVS_CVODES::STATUS::NORMAL; auto const F = I.val_function(); auto const Gf = I.val_function_gradient();
  bool const oa = I.solve_asens( Q ) == ODESLVS_CVODES::STATUS::NORMAL;     auto const Ga = I.val_function_gradient();
  auto const F0 = closed( P0 );
  for( size_t j = 0; j < 5; ++j ) check( of && std::fabs( F[j]-F0[j] ) < 1e-9, std::string( "value " ) + nm[j] + " == closed form", of? std::fabs( F[j]-F0[j] ): 1. );
  for( size_t q = 0; q < 6; ++q ){
    auto Pp = P0, Pm = P0; double const h = 1e-6*std::max( 1., std::fabs( P0[q] ) ); Pp[q] += h; Pm[q] -= h;
    auto const fp = closed( Pp ), fm = closed( Pm ); double ef = 0., ea = 0.;
    for( size_t j = 0; j < 5; ++j ){ double const D = ( fp[j]-fm[j] )/( 2.*h ); ef = std::max( ef, of? std::fabs( Gf[ix[q]][j]-D ): 1. ); ea = std::max( ea, oa? std::fabs( Ga[ix[q]][j]-D ): 1. ); }
    check( ef < 1e-7, std::string( "FORWARD d/d" ) + pn[q] + " of every output == closed form", ef );
    check( ea < 1e-7, std::string( "ADJOINT d/d" ) + pn[q] + " of every output == closed form", ea );
  }
  std::printf( "\n  TRANS_multi: %d passed, %d failed -- %s\n", npass, nfail, nfail? "FAILURES": "ALL PASS" );
  return nfail? 1: 0;
}
