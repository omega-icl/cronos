// ODESLV_stageout.cpp -- GATE for OUTPUTS AT STAGE TIMES (WORKPLAN 1.7, 1.8, 2.2; 2026-10-06).
// ODESLV has no collocation nodes, so an output DISTRIBUTED over the evolution direction is carried at the STAGE TIMES:
// one value per kept stage time (the end of the stage; the initial state at t0), the mask selecting them.  Before, such an
// output was accepted and silently collapsed to ONE value (the final one).  Also checked here: an output AT THE INITIAL TIME
// (stage 0) had a wrong ADJOINT gradient -- the chain through the initial state multiplied by the initial VALUE instead of
// d(x0)/dp -- which is why the initial state depends on the parameter (x0 = 2 k), not on a constant.
// Model: x' = -k x, x(0) = 2 k  =>  x(t) = 2 k exp(-k t);  4 stages on [0, 1] (stage times 0, .25, .5, .75, 1).
//   rows: 0 point x(.25) | 1..5 x on ALL | 6..9 x on ALL-LB | 10 x on LB | 11 x on UB | 12..14 x^2 on ALL-LB-UB | 15 point x(1)
#include <cmath>
#include <cstdio>
#include <iostream>
#include <sstream>
#include <vector>
#include "ffode.hpp"
using namespace mc;
static int nfail = 0;
static void check( char const* what, bool ok, double err = -1. ){
  if( err >= 0. ) std::printf( "  %-70s %s (%.1e)\n", what, ok? "PASS": "FAIL", err ); else std::printf( "  %-70s %s\n", what, ok? "PASS": "FAIL" );
  nfail += !ok; }
struct Model { FFGraph G; ODESLVS_CVODES I{ &G }; FFVar t, x, k; };
static bool build( Model& m ){
  m.t = m.G.add_var( "t" ); m.x = m.G.add_var( "x(t)" ); m.k = m.G.add_var( "k" );
  FFPartial OpP; FFEval OpE;
  m.I.add_domain( m.t, FFDom( 0., 1., 4, FFDom::LGR, 3 ) ); m.I.set_evolution_domain( m.t );
  m.I.add_state( m.x, {m.t} ); m.I.add_input( m.k );  m.I.update_ref( m.x, 1. );
  int const T_INT = FFDom::ALL - FFDom::LB;
  m.I.add_equation( OpP( m.x, m.t ) + m.k*m.x, {m.t}, {T_INT}, FFModel::EqnOptions( FFModel::EqnRole::INTERIOR, 0 ) );
  m.I.add_equation( m.x - 2.*m.k,              {m.t}, {(int)FFDom::LB}, FFModel::EqnOptions( FFModel::EqnRole::INITIAL, 0 ) );
  size_t const i0 = m.I.add_output( OpE( m.x, m.t, 0.25 ) );
  size_t const i1 = m.I.add_output( m.x, {m.t}, {(int)FFDom::ALL} );
  size_t const i2 = m.I.add_output( m.x, {m.t}, {(int)( FFDom::ALL - FFDom::LB )} );
  size_t const i3 = m.I.add_output( m.x, {m.t}, {(int)FFDom::LB} );
  size_t const i4 = m.I.add_output( m.x, {m.t}, {(int)FFDom::UB} );
  size_t const i5 = m.I.add_output( m.x*m.x, {m.t}, {(int)( FFDom::ALL - FFDom::LB - FFDom::UB )} );
  size_t const i6 = m.I.add_output( OpE( m.x, m.t, 1. ) );
  check( "add_output returns the output's index: 0..6", i0 == 0 && i1 == 1 && i2 == 2 && i3 == 3 && i4 == 4 && i5 == 5 && i6 == 6 );
  m.I.options.DISPLAY = 0; m.I.options.ATOL = m.I.options.ATOLS = m.I.options.ATOLB = 1e-12; m.I.options.RTOL = m.I.options.RTOLS = m.I.options.RTOLB = 1e-12;
  std::ostringstream os; auto* old = std::cerr.rdbuf( os.rdbuf() ); bool const ok = m.I.setup(); std::cerr.rdbuf( old );
  return ok; }
int main(){
  std::printf( "ODESLV_stageout -- outputs at stage times; the adjoint at the initial time\n" );
  double const K = 0.8;  double const TS[5] = { 0., .25, .5, .75, 1. };
  auto x  = [&]( double t ){ return 2.*K*std::exp( -K*t ); };
  auto dx = [&]( double t ){ return 2.*std::exp( -K*t )*( 1. - K*t ); };                  // d x / d k
  auto x2 = [&]( double t ){ return 4.*K*K*std::exp( -2.*K*t ); };
  auto dx2= [&]( double t ){ return 8.*K*std::exp( -2.*K*t )*( 1. - K*t ); };
  Model m;  bool const ok = build( m );  check( "setup", ok );  if( !ok ){ std::printf( "  %s\n", m.I.extract_error().c_str() ); return 1; }
  // the blocks
  std::pair<size_t,size_t> const B[7] = { {0,1}, {1,5}, {6,4}, {10,1}, {11,1}, {12,3}, {15,1} };
  bool blk = true;  for( size_t i = 0; i < 7; ++i ) blk = blk && m.I.blk_fct( i ) == B[i];
  check( "blk_fct: (0,1) (1,5) (6,4) (10,1) (11,1) (12,3) (15,1) -- one row per kept stage time", blk );
  check( "blk_fct beyond the last output: the sentinel", m.I.blk_fct( 7 ).first == std::numeric_limits<size_t>::max() );
  std::vector<double> P( m.I.np() );  P[ m.I.parameter_index( m.k )[0] ] = K;
  // the exact value of every row
  std::vector<double> ex( 16 ), dex( 16 );
  ex[0] = x(.25); dex[0] = dx(.25);
  for( size_t j = 0; j < 5; ++j ){ ex[1+j] = x(TS[j]); dex[1+j] = dx(TS[j]); }
  for( size_t j = 0; j < 4; ++j ){ ex[6+j] = x(TS[j+1]); dex[6+j] = dx(TS[j+1]); }
  ex[10] = x(0.); dex[10] = dx(0.);  ex[11] = x(1.); dex[11] = dx(1.);
  for( size_t j = 0; j < 3; ++j ){ ex[12+j] = x2(TS[j+1]); dex[12+j] = dx2(TS[j+1]); }
  ex[15] = x(1.); dex[15] = dx(1.);
  bool const of = m.I.solve_fsens( P ) == ODESLVS_CVODES::STATUS::NORMAL;
  auto const F = m.I.val_function();  auto const Gf = m.I.val_function_gradient();
  // the gradient is stored [parameter][function]: here one parameter, 16 functions
  auto grad = []( std::vector<std::vector<double>> const& G, size_t row ){ return G.size() == 16? G[row][0]: G[0][row]; };
  auto nrows = []( std::vector<std::vector<double>> const& G ){ return G.size() == 1? G[0].size(): G.size(); };
  double ev = 0., egf = 0., eg0f = 0.;
  for( size_t r = 0; r < 16 && F.size() == 16; ++r ){ ev = std::max( ev, std::fabs( F[r] - ex[r] ) ); egf = std::max( egf, std::fabs( grad( Gf, r ) - dex[r] ) ); }
  eg0f = std::fabs( grad( Gf, 10 ) - dex[10] );
  check( "val_function has 16 rows (one per kept stage time)", of && F.size() == 16 );
  check( "every row = the closed form (1e-8)", of && F.size() == 16 && ev < 1e-8, ev );
  check( "forward gradient of every row = the closed form (1e-7)", of && F.size() == 16 && egf < 1e-7, egf );
  bool const oa = m.I.solve_asens( P ) == ODESLVS_CVODES::STATUS::NORMAL;
  auto const Ga = m.I.val_function_gradient();  double ega = 0., eg0a = 0.;
  for( size_t r = 0; r < 16 && oa; ++r ) ega = std::max( ega, std::fabs( grad( Ga, r ) - dex[r] ) );
  eg0a = std::fabs( grad( Ga, 10 ) - dex[10] );
  check( "adjoint gradient of every row = the closed form (1e-7)", oa && nrows( Ga ) == 16 && ega < 1e-7, ega );
  check( "adjoint gradient of the output AT THE INITIAL TIME: d(2k)/dk = 2 (was x0*df/dx = 1.6)", oa && nrows( Ga ) == 16 && eg0a < 1e-9, eg0a );
  // refusals: a distributed output must not involve an evaluation or an integral
  { FFGraph G; ODESLVS_CVODES I( &G ); FFVar t = G.add_var( "t" ), x = G.add_var( "x" ), k = G.add_var( "k" ); FFPartial OpP; FFEval OpE;
    I.add_domain( t, FFDom( 0., 1., 2, FFDom::LGR, 3 ) ); I.set_evolution_domain( t ); I.add_state( x, {t} ); I.add_input( k ); I.update_ref( x, 1. );
    I.add_equation( OpP( x, t ) + k*x, {t}, {(int)( FFDom::ALL - FFDom::LB )}, FFModel::EqnOptions( FFModel::EqnRole::INTERIOR, 0 ) );
    I.add_equation( x - 1., {t}, {(int)FFDom::LB}, FFModel::EqnOptions( FFModel::EqnRole::INITIAL, 0 ) );
    I.add_output( OpE( x, t, 0.5 ), {t}, {(int)FFDom::ALL} );
    std::ostringstream os; auto* old = std::cerr.rdbuf( os.rdbuf() ); bool const ok2 = I.setup(); std::cerr.rdbuf( old );
    check( "a distributed output of an evaluation is REFUSED at setup", !ok2 ); }
  std::printf( "  ODESLV_stageout: %d failed\n", nfail );
  return nfail? 1: 0;
}
