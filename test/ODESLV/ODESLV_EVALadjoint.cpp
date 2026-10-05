// EVAL_adjoint.cpp -- ODESLV outputs with POINTWISE terms at INTERMEDIATE times: each term makes the adjoint jump by
// dF/dx there.  No state discontinuity.  Values + FORWARD and ADJOINT gradients against a CLOSED FORM (central
// differences of it).   dx1/dt = -k1 x1, x1(0) = 1;  dx2/dt = x1 - k2 x2, x2(0) = 0;  mesh {0, 0.5, 1}
//   x1 = e^{-k1 t},   x2 = ( e^{-k1 t} - e^{-k2 t} ) / ( k2 - k1 )
//   F0 = x1(0.3)^2 + 2 x2(0.7) + int x1 dt     (nonlinear pointwise term, interior-time term, integral)
//   F1 = w (x1 x2)(0.45)                        (parameter-weighted pointwise term)
//   F2 = x2(0.5) - x1(0.8)                      (element-boundary term, interior term)
//   F3 = x2(1)
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>
#include "ffode.hpp"
using namespace mc;
static int npass = 0, nfail = 0;
static void check( bool c, std::string const& w, double v ){ std::printf( "  %s  %-62s (%.2e)\n", c? "PASS": "FAIL", w.c_str(), v ); c? ++npass: ++nfail; }
static std::vector<double> closed( std::vector<double> const& P ){          // P = k1, k2, w
  double const k1 = P[0], k2 = P[1], w = P[2];
  auto x1 = [&]( double t ){ return std::exp( -k1*t ); };
  auto x2 = [&]( double t ){ return ( std::exp( -k1*t ) - std::exp( -k2*t ) )/( k2 - k1 ); };
  return { x1( 0.3 )*x1( 0.3 ) + 2.*x2( 0.7 ) + ( 1.-std::exp( -k1 ) )/k1, w*x1( 0.45 )*x2( 0.45 ), x2( 0.5 ) - x1( 0.8 ), x2( 1. ) }; }
int main(){
  std::vector<double> const P0{ 1.3, 0.6, 1.7 };
  char const* nm[4] = { "F0 = x1(.3)^2 + 2 x2(.7) + int x1", "F1 = w (x1 x2)(.45)", "F2 = x2(.5) - x1(.8)", "F3 = x2(1)" };
  char const* pn[3] = { "k1", "k2", "w" };
  FFGraph G; ODESLVS_CVODES I( &G ); FFPartial OpP; FFEval OpE; FFIntegral OpI;
  FFVar t = G.add_var( "t" ), x1 = G.add_var( "x1(t)" ), x2 = G.add_var( "x2(t)" );
  std::vector<FFVar> p; for( auto n : pn ) p.push_back( G.add_var( n ) );
  I.add_domain( t, FFDom( std::vector<double>{ 0., 0.5, 1. }, FFDom::LGR, 4 ) ); I.set_evolution_domain( t );
  I.add_state( x1, {t} ); I.add_state( x2, {t} ); for( auto const& v : p ) I.add_input( v );
  I.update_ref( x1, 1. ); I.update_ref( x2, 0. );
  int const T_INT = FFDom::ALL - FFDom::LB;
  I.add_equation( OpP( x1, t ) + p[0]*x1,      {t}, {T_INT}, FFModel::EqnOptions( FFModel::EqnRole::INTERIOR, 0 ) );
  I.add_equation( OpP( x2, t ) - x1 + p[1]*x2, {t}, {T_INT}, FFModel::EqnOptions( FFModel::EqnRole::INTERIOR, 0 ) );
  I.add_equation( x1 - 1., {t}, {FFDom::LB}, FFModel::EqnOptions( FFModel::EqnRole::INITIAL, 0 ) );
  I.add_equation( x2,      {t}, {FFDom::LB}, FFModel::EqnOptions( FFModel::EqnRole::INITIAL, 0 ) );
  I.add_output( OpE( x1*x1, t, 0.3 ) + 2.*OpE( x2, t, 0.7 ) + OpI( std::vector<FFVar>{ x1 }, { t } )[0] );
  I.add_output( p[2]*OpE( x1*x2, t, 0.45 ) );
  I.add_output( OpE( x2, t, 0.5 ) - OpE( x1, t, 0.8 ) );
  I.add_output( OpE( x2, t, 1. ) );
  I.options.DISPLAY = 0; I.options.ATOL = I.options.ATOLS = I.options.ATOLB = 1e-12; I.options.RTOL = I.options.RTOLS = I.options.RTOLB = 1e-12;
  bool ok = I.setup(); check( ok, "setup", 0. ); if( !ok ){ std::printf( "  %s\n", I.extract_error().c_str() ); return 1; }
  std::vector<double> Q( I.np() ); std::vector<size_t> ix; for( size_t q = 0; q < 3; ++q ){ ix.push_back( I.parameter_index( p[q] )[0] ); Q[ix[q]] = P0[q]; }
  bool const of = I.solve_fsens( Q ) == ODESLVS_CVODES::STATUS::NORMAL; auto const F = I.val_function(); auto const Gf = I.val_function_gradient();
  bool const oa = I.solve_asens( Q ) == ODESLVS_CVODES::STATUS::NORMAL;     auto const Ga = I.val_function_gradient();
  auto const F0 = closed( P0 );
  for( size_t j = 0; j < 4; ++j ) check( of && std::fabs( F[j]-F0[j] ) < 1e-9, std::string( "value " ) + nm[j], of? std::fabs( F[j]-F0[j] ): 1. );
  for( size_t q = 0; q < 3; ++q ){
    auto Pp = P0, Pm = P0; double const h = 1e-6*std::max( 1., std::fabs( P0[q] ) ); Pp[q] += h; Pm[q] -= h;
    auto const fp = closed( Pp ), fm = closed( Pm ); double ef = 0., ea = 0.;
    for( size_t j = 0; j < 4; ++j ){ double const D = ( fp[j]-fm[j] )/( 2.*h ); ef = std::max( ef, of? std::fabs( Gf[ix[q]][j]-D ): 1. ); ea = std::max( ea, oa? std::fabs( Ga[ix[q]][j]-D ): 1. ); }
    check( ef < 1e-7, std::string( "FORWARD d/d" ) + pn[q] + " of every output == closed form", ef );
    check( ea < 1e-7, std::string( "ADJOINT d/d" ) + pn[q] + " of every output == closed form", ea );
  }
  std::printf( "\n  EVAL_adjoint: %d passed, %d failed -- %s\n", npass, nfail, nfail? "FAILURES": "ALL PASS" );
  return nfail? 1: 0;
}
