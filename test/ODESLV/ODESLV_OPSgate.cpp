// OPS_gate.cpp -- EVERY smooth nonlinear operation of a MULTI-node input is formed at the collocation nodes (OCFESLV,
// monolithic and marching).  u piecewise LINEAR (LGL, 2 nodes/element, 0.2 < u < 0.8) on {0, 0.3, 0.6, 1}; one state
// per function:  dx_k/dt = f_k(u),  x_k(0) = 0,  so  x_k(1) = int f_k(u) dt.  ORACLE: the true integral (20-point
// Gauss-Legendre per element), with 12 collocation nodes per element so the scheme's own quadrature error is far below
// the 1e-7 tolerance for these smooth f (~1e-8 at worst, for 1/u and log u near u = 0.25); formed at the input's own 2 nodes instead, each is a trapezoidal rule, ~1e-2
// off.  (Checked with a control run on headers without the lift.)
// A new or changed nonlinear OCVar operation that skips OCVar::_nl_lift shows up here as a FAIL.
#include <cmath>
#include <cstdio>
#include <functional>
#include <vector>
#include "ocfeslv.hpp"
using namespace mc;
static int npass = 0, nfail = 0;
std::vector<double> const EB{ 0., 0.3, 0.6, 1. }, LV{ 0.25, 0.70, 0.35, 0.55, 0.78, 0.30 };   // u[e][j]
struct Fn { char const* name; std::function<FFVar( FFVar const& )> sym; std::function<double( double )> num; };
int main(){
  std::vector<Fn> F{
    { "sqr(u)",       []( FFVar const& u ){ return sqr( u ); },            []( double u ){ return u*u; } },
    { "u*u",          []( FFVar const& u ){ return u*u; },                 []( double u ){ return u*u; } },
    { "pow(u,3)",     []( FFVar const& u ){ return pow( u, 3 ); },         []( double u ){ return u*u*u; } },
    { "pow(u,1.5)",   []( FFVar const& u ){ return pow( u, 1.5 ); },       []( double u ){ return std::pow( u, 1.5 ); } },
    { "exp(u)",       []( FFVar const& u ){ return exp( u ); },            []( double u ){ return std::exp( u ); } },
    { "log(u)",       []( FFVar const& u ){ return log( u ); },            []( double u ){ return std::log( u ); } },
    { "sqrt(u)",      []( FFVar const& u ){ return sqrt( u ); },           []( double u ){ return std::sqrt( u ); } },
    { "sin(3u)",      []( FFVar const& u ){ return sin( 3.*u ); },         []( double u ){ return std::sin( 3.*u ); } },
    { "cos(3u)",      []( FFVar const& u ){ return cos( 3.*u ); },         []( double u ){ return std::cos( 3.*u ); } },
    { "tan(u)",       []( FFVar const& u ){ return tan( u ); },            []( double u ){ return std::tan( u ); } },
    { "asin(u)",      []( FFVar const& u ){ return asin( u ); },           []( double u ){ return std::asin( u ); } },
    { "acos(u)",      []( FFVar const& u ){ return acos( u ); },           []( double u ){ return std::acos( u ); } },
    { "atan(3u)",     []( FFVar const& u ){ return atan( 3.*u ); },        []( double u ){ return std::atan( 3.*u ); } },
    { "sinh(2u)",     []( FFVar const& u ){ return sinh( 2.*u ); },        []( double u ){ return std::sinh( 2.*u ); } },
    { "cosh(2u)",     []( FFVar const& u ){ return cosh( 2.*u ); },        []( double u ){ return std::cosh( 2.*u ); } },
    { "tanh(3u)",     []( FFVar const& u ){ return tanh( 3.*u ); },        []( double u ){ return std::tanh( 3.*u ); } },
    { "xlog(u)",      []( FFVar const& u ){ return xlog( u ); },           []( double u ){ return u*std::log( u ); } },
    { "erf(2u)",      []( FFVar const& u ){ return erf( 2.*u ); },         []( double u ){ return std::erf( 2.*u ); } },
    { "erfc(2u)",     []( FFVar const& u ){ return erfc( 2.*u ); },        []( double u ){ return std::erfc( 2.*u ); } },
    { "inv(u)",       []( FFVar const& u ){ return inv( u ); },            []( double u ){ return 1./u; } },
    { "1/u",          []( FFVar const& u ){ return 1./u; },                []( double u ){ return 1./u; } },
    { "u/(1+u)",      []( FFVar const& u ){ return u/( 1.+u ); },          []( double u ){ return u/( 1.+u ); } },
    { "2u+3 (LINEAR)",[]( FFVar const& u ){ return 2.*u + 3.; },           []( double u ){ return 2.*u + 3.; } } };
  static double const gx[10] = { 0.0765265211334973, 0.2277858511416451, 0.3737060887154195, 0.5108670019508271, 0.6360536807265150,
                                 0.7463319064601508, 0.8391169718222188, 0.9122344282513259, 0.9639719272779138, 0.9931285991850949 };
  static double const gw[10] = { 0.1527533871307258, 0.1491729864726037, 0.1420961093183820, 0.1316886384491766, 0.1181945319615184,
                                 0.1019301198172404, 0.0832767415767048, 0.0626720483341091, 0.0406014298003869, 0.0176140071391521 };
  auto exact = [&]( std::function<double( double )> const& f ){ double I = 0.;
    for( size_t e = 0; e < 3; ++e ) for( int k = 0; k < 10; ++k ) for( int sg : { -1, 1 } ){
      double const h = EB[e+1]-EB[e], s = 0.5*( 1. + sg*gx[k] ), u = LV[2*e]*( 1.-s ) + LV[2*e+1]*s; I += 0.5*h*gw[k]*f( u ); }
    return I; };
  for( bool march : { false, true } ){
    FFGraph G; OCFESLV I( &G ); FFVar t = G.add_var( "t" ), u = G.add_var( "u(t)" ); FFPartial OpP; FFEval OpE;
    I.add_domain( t, FFDom( EB, FFDom::LGR, 12 ) ); I.set_evolution_domain( t ); I.add_input( u, {t}, FFDom::LGL, 2 );
    int const T_INT = FFDom::ALL - FFDom::LB;
    std::vector<FFVar> X;
    for( size_t k = 0; k < F.size(); ++k ){
      X.push_back( G.add_var( "x" + std::to_string( k ) + "(t)" ) ); I.add_state( X.back(), {t} ); I.update_ref( X.back(), 0. );
      I.add_equation( OpP( X.back(), t ) - F[k].sym( u ), {t}, {T_INT}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
      I.add_equation( X.back(), {t}, {FFDom::LB}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INITIAL, 0 ) );
      I.add_output( OpE( X.back(), t, 1. ) );
    }
    I.options.SOLVE.MARCHING = march; I.options.SOLVE.RES_TOL = 1e-11; I.options.DISPLAY_LEVEL = 0;
    if( !I.setup() ){ std::printf( "  setup failed\n" ); ++nfail; continue; }
    std::vector<double> xv, inp; I.init( xv, inp, nullptr ); I.set_input_values( u, LV, inp.data() );
    bool const ok = I.solve( xv.data(), inp.data(), nullptr ).converged; auto const V = I.val_functions();
    for( size_t k = 0; k < F.size(); ++k ){
      double const ex = exact( F[k].num ), err = ok? std::fabs( V[k]-ex )/std::max( 1., std::fabs( ex ) ): 1.;
      bool const c = err < 1e-7;   // scheme quadrature error of steep f near u = 0.25 is ~1e-8; unlifted errors are >= 1e-2
      std::printf( "  %s  %-10s %-14s x(1) = %.13f   int f(u) = %.13f   rel %.2e\n", c? "PASS": "FAIL", march? "marching": "monolithic", F[k].name, ok? V[k]: 0., ex, err );
      c? ++npass: ++nfail;
    }
  }
  std::printf( "\n  OPS_gate: %d passed, %d failed -- %s\n", npass, nfail, nfail? "FAILURES": "ALL PASS" );
  return nfail? 1: 0;
}
