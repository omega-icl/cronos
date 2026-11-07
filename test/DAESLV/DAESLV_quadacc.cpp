// DAESLV_quadacc.cpp -- quadratures ACCUMULATE over the horizon (4 Nov 2026), in ODESLV and DAESLV alike.
//   x' = -k x, x(0) = X0, k a parameter; outputs f0 = int x dt and f1 = x(1) + int x dt (two reductions side by side),
//   on 1, 2 and 4 stages; values, forward and adjoint gradients against closed forms.
// Before: the quadrature was zeroed at every stage and its per-stage increments summed.  DAESLV then failed at a
// stage boundary without a state restart once the integrand was large against ATOL (x0 = 10, ATOL = 1e-12), and both
// integrators got f1 wrong by 62% on any multi-stage horizon.
// The model is also EXACTLY the case that stalled IDAS's order selection (k = 1, unit weights: the costate of f1 is
// constant, lambda = 1, while its parameter quadrature varies): unpatched IDAS decided an order decrease on the state
// norms alone (0 <= 0) and its quadrature test could not overrule it, so DAESLV's adjoint hit mxstep at order 1.  Needs
// SUNDIALS with the IDAS order-decision patch (sundials-idas-order-decision-20261104.patch) until upstream fixes it.
#include <iostream>
#include <iomanip>
#include <cmath>
#include "odeslv_cvodes.hpp"
#include "daeslv_idas.hpp"
using namespace mc;
static int nfail = 0;
static double X0G = 10., ATOLG = 1e-12;
template <typename S>
void run( char const* name, size_t NST )
{
  double const X0 = X0G, K = 1.;
  FFGraph G;  S s( &G );
  FFVar t = G.add_var("t"), x = G.add_var("x"), k = G.add_var("k");
  s.add_domain( t, FFDom( 0., 1., NST, FFDom::LGR, 3 ) );
  s.set_evolution_domain( t );
  s.add_state( x, {t} );  s.update_ref( x, X0 );
  s.add_input( k, {} );
  FFPartial OpP;  FFIntegral OpI;
  typename S::EqnOptions io( S::EqnRole::INTERIOR, 0 ), ii( S::EqnRole::INITIAL, 0 );
  s.add_equation( OpP(x,t) + k*x, {t}, {FFDom::ALL - FFDom::LB}, io );
  s.add_equation( x - X0,         {t}, {FFDom::LB},              ii );
  s.add_output( OpI( x, t ) );
  FFEval OpE;
  s.add_output( OpE( x, t, 1. ) + OpI( x, t ) );                  // two reductions side by side (not nested)
  s.options.SOLVE.RTOL = 1e-10;  s.options.SOLVE.ATOL = ATOLG;  s.options.SOLVE.VERBOSE = 0;  s.options.DISPLAY_LEVEL = 0;
  double const e = std::exp( -K ), I = X0*( 1. - e )/K, dI = X0*( e/K - ( 1. - e )/( K*K ) );
  double const f[2] = { I, X0*e + I }, df[2] = { dI, -X0*e + dI };
  if( !s.setup() ){ std::cout << name << " " << NST << " stages: setup failed -- " << s.extract_error() << "\n"; ++nfail; return; }
  std::vector<double> P{ K };
  auto line = [&]( char const* what, bool ok, double err ){
    std::cout << "  " << std::left << std::setw(7) << name << NST << " stage(s)  " << std::setw(10) << what << ( ok? "": " FAILED (" + s.extract_error() + ")" )
              << "  max rel err " << std::scientific << std::setprecision(2) << err;
    bool const adj = std::string( what ) == "adjoint";
    double const tol = adj? 1e-5: 1e-7;                   // an adjoint is as accurate as its own tolerances, ~1e-6
    std::cout << ( ok && err < tol? "  PASS": "  FAIL" ) << "\n";
    if( !( ok && err < tol ) ) ++nfail; };
  bool ok = s.solve( P ) == S::NORMAL;
  double ev = 0.;  if( ok ) for( int i = 0; i < 2; ++i ) ev = std::max( ev, std::fabs( s.val_function()[i] - f[i] )/std::fabs( f[i] ) );
  line( "value", ok, ev );
  for( int mode = 0; mode < 2; ++mode ){
    ok = ( mode? s.solve_asens( P ): s.solve_fsens( P ) ) == S::NORMAL;
    double eg = 0.;
    if( ok ){
      auto const& Gr = s.val_function_gradient();
      for( int i = 0; i < 2; ++i ){
        double const g = Gr.size() == 2? Gr[i].at( 0 ): Gr.at( 0 ).at( i );          // [output][param] or [param][output]
        eg = std::max( eg, std::fabs( g - df[i] )/std::fabs( df[i] ) ); } }
    line( mode? "adjoint": "forward", ok, eg );
  }
}
// NONLINEAR functions of horizon integrals (7 Nov 2026; refused before, Z10: one record or an affine combination):
// A = int_0^1 x dt, sqrt(A) and A/B with B = int_0^1 x^2 dt on x' = -k x, x(0) = 1 -- values and d/dk against the
// closed forms, forward and adjoint, on both solvers.  Each record a hidden function, the output finished from them.
template <typename S>
void run_nonlinear( char const* name )
{
  double const K = 0.7;
  FFGraph G;  S s( &G );  FFPartial OpP;  FFIntegral OpI;
  FFVar t = G.add_var("t"), x = G.add_var("x"), k = G.add_var("k");
  s.add_domain( t, FFDom( 0., 1., 2, FFDom::LGR, 3 ) );
  s.set_evolution_domain( t );
  s.add_state( x, {t} );  s.update_ref( x, 1. );
  s.add_input( k, {} );
  typename S::EqnOptions io( S::EqnRole::INTERIOR, 0 ), ii( S::EqnRole::INITIAL, 0 );
  s.add_equation( OpP(x,t) + k*x, {t}, {FFDom::ALL - FFDom::LB}, io );
  s.add_equation( x - 1.,         {t}, {FFDom::LB},              ii );
  FFVar A = OpI( x, t ), B = OpI( x*x, t );
  s.add_output( A );  s.add_output( sqrt( A ) );  s.add_output( A / B );
  s.options.SOLVE.RTOL = 1e-10;  s.options.SOLVE.ATOL = ATOLG;  s.options.SOLVE.VERBOSE = 0;  s.options.DISPLAY_LEVEL = 0;
  double const Ae = ( 1 - std::exp(-K) )/K, Be = ( 1 - std::exp(-2*K) )/( 2*K );
  double const dA = -( Ae - std::exp(-K) )/K, dB = -( Be - std::exp(-2*K) )/K;
  double const f[3] = { Ae, std::sqrt( Ae ), Ae/Be }, df[3] = { dA, dA/( 2*std::sqrt( Ae ) ), ( dA*Be - Ae*dB )/( Be*Be ) };
  std::vector<double> const P{ K };
  auto line = [&]( char const* what, bool ok, double err, double tol ){
    ++nfail;  if( ok && err < tol ) --nfail;
    std::cout << "  " << std::left << std::setw(7) << name << "nonlinear  " << std::setw(10) << what << ( ok? "": " FAILED (" + s.extract_error() + ")" )
              << "max rel err " << std::scientific << std::setprecision(2) << err << ( ok && err < tol? "  PASS": "  FAIL" ) << "\n"; };
  if( !s.setup() ){ line( "setup", false, 0., 0. ); return; }
  bool ok = s.solve( P ) == S::NORMAL;
  double ev = 0.;  if( ok ) for( int i = 0; i < 3; ++i ) ev = std::max( ev, std::fabs( s.val_function()[i] - f[i] )/std::fabs( f[i] ) );
  line( "value", ok && s.nf() == 3, ev, 1e-8 );
  for( int mode = 0; mode < 2; ++mode ){
    ok = ( mode? s.solve_asens( P ): s.solve_fsens( P ) ) == S::NORMAL;
    double eg = 0.;
    if( ok ){
      auto const& Gr = s.val_function_gradient();
      for( int i = 0; i < 3; ++i ){
        double const g = Gr.size() == 3? Gr[i].at( 0 ): Gr.at( 0 ).at( i );
        eg = std::max( eg, std::fabs( g - df[i] )/std::fabs( df[i] ) ); } }
    line( mode? "adjoint": "forward", ok, eg, mode? 1e-6: 1e-7 );
  }
}

int main( int argc, char* argv[] )
{
  if( argc > 1 ) X0G = std::atof( argv[1] );
  if( argc > 2 ) ATOLG = std::atof( argv[2] );
  for( size_t n : { 1, 2, 4 } ){ run<ODESLV_CVODES>( "ODESLV", n );  run<DAESLV_IDAS>( "DAESLV", n ); }
  run_nonlinear<ODESLV_CVODES>( "ODESLV" );  run_nonlinear<DAESLV_IDAS>( "DAESLV" );
  std::cout << ( nfail? "FAILURES: ": "ALL PASS" ) << ( nfail? std::to_string( nfail ): std::string() ) << "\n";
  return nfail? 1: 0;
}
