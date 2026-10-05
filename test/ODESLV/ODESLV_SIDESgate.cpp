// SIDES_gate.cpp -- one-sided point evaluation (FFDom::MINUS / FFDom::PLUS), ODESLV and OCFESLV (monolithic and
// marching), against EXACT answers.  dx/dt = u,  x(0) = 0,  u piecewise constant on elements {0, 0.5, 1}: u0 = 2, u1 = 5.
//   outputs: u(0.5-) = 2, u(0.5+) = 5, x(0.5-) = x(0.5+) = 1, u(0.3-) = u(0.3+) = 2 (not a boundary),
//            u(1+) = u(1-) = 5 (the upper bound: the one limit that exists), int u dt = 3.5
//   gradients w.r.t. (u0, u1), forward and adjoint: exact.
#include <cmath>
#include <cstdio>
#include <vector>
#include "ffode.hpp"
#include "ocfeslv.hpp"
using namespace mc;
static int npass = 0, nfail = 0;
static void check( bool c, char const* w, double v ){ std::printf( "  %s  %-70s (%.2e)\n", c? "PASS": "FAIL", w, v ); c? ++npass: ++nfail; }
double const U0 = 2., U1 = 5.;
// value and d/d(u0,u1) of each output, exactly
std::vector<double> const VAL{ 2., 5., 1., 1., 2., 2., 5., 5., 3.5 };
std::vector<std::vector<double>> const GRD{ {1,0}, {0,1}, {.5,0}, {.5,0}, {1,0}, {1,0}, {0,1}, {0,1}, {.5,.5} };
char const* const NAME[9] = { "u(0.5-)", "u(0.5+)", "x(0.5-)", "x(0.5+)", "u(0.3-)", "u(0.3+)", "u(1-)", "u(1+)", "int u" };
template <class M> static void outputs( M& I, FFVar const& t, FFVar const& x, FFVar const& u ){
  FFEval OpE; FFIntegral OpI;
  I.add_output( OpE( u, t, 0.5, FFDom::MINUS ) ); I.add_output( OpE( u, t, 0.5, FFDom::PLUS ) );
  I.add_output( OpE( x, t, 0.5, FFDom::MINUS ) ); I.add_output( OpE( x, t, 0.5, FFDom::PLUS ) );
  I.add_output( OpE( u, t, 0.3, FFDom::MINUS ) ); I.add_output( OpE( u, t, 0.3, FFDom::PLUS ) );
  I.add_output( OpE( u, t, 1.0, FFDom::MINUS ) ); I.add_output( OpE( u, t, 1.0, FFDom::PLUS ) );
  I.add_output( OpI( std::vector<FFVar>{ u }, { t } )[0] ); }
static void report( char const* who, bool ok, std::vector<double> const& F, std::vector<std::vector<double>> const& Gf, std::vector<std::vector<double>> const& Ga ){
  char w[140];
  if( !ok || F.size() != 9 ){ std::snprintf( w, sizeof w, "%s: solves and 9 outputs", who ); check( false, w, 0. ); return; }
  for( size_t k = 0; k < 9; ++k ){
    double e = std::fabs( F[k] - VAL[k] );
    for( size_t d = 0; d < 2; ++d ){ if( !Gf.empty() ) e = std::max( e, std::fabs( Gf[k][d] - GRD[k][d] ) ); if( !Ga.empty() ) e = std::max( e, std::fabs( Ga[k][d] - GRD[k][d] ) ); }
    std::snprintf( w, sizeof w, "%s: %-8s value %g, forward and adjoint gradients exact", who, NAME[k], VAL[k] );
    check( e < 1e-6, w, e );
  }
}
int main(){
  { // ---- ODESLV
    FFGraph G; ODESLVS_CVODES I( &G ); FFVar t = G.add_var( "t" ), x = G.add_var( "x(t)" ), u = G.add_var( "u(t)" ); FFPartial OpP;
    I.add_domain( t, FFDom( std::vector<double>{ 0., 0.5, 1. }, FFDom::LGR, 4 ) ); I.set_evolution_domain( t );
    I.add_state( x, {t} ); I.add_input( u, {t}, FFDom::LGR, 1 ); I.update_ref( x, 0. );
    int const T_INT = FFDom::ALL - FFDom::LB;
    I.add_equation( OpP( x, t ) - u, {t}, {T_INT}, FFModel::EqnOptions( FFModel::EqnRole::INTERIOR, 0 ) );
    I.add_equation( x, {t}, {FFDom::LB}, FFModel::EqnOptions( FFModel::EqnRole::INITIAL, 0 ) );
    outputs( I, t, x, u ); I.options.DISPLAY = 0; I.options.ATOL = I.options.ATOLS = I.options.ATOLB = 1e-12; I.options.RTOL = I.options.RTOLS = I.options.RTOLB = 1e-12;
    bool ok = I.setup(); std::vector<double> P( I.np() ); auto iu = I.parameter_index( u ); P[iu[0]] = U0; P[iu[1]] = U1;
    ok = ok && I.solve_fsens( P ) == ODESLVS_CVODES::STATUS::NORMAL; auto const F = I.val_function(); auto Gf0 = I.val_function_gradient();
    ok = ok && I.solve_asens( P ) == ODESLVS_CVODES::STATUS::NORMAL; auto Ga0 = I.val_function_gradient();
    std::vector<std::vector<double>> Gf( 9, std::vector<double>( 2 ) ), Ga( Gf );          // [output][u0,u1]
    for( size_t k = 0; ok && k < 9; ++k ) for( size_t d = 0; d < 2; ++d ){ Gf[k][d] = Gf0[iu[d]][k]; Ga[k][d] = Ga0[iu[d]][k]; }
    report( "ODESLV          ", ok, F, Gf, Ga );
  }
  for( bool march : { false, true } ){ // ---- OCFESLV
    FFGraph G; OCFESLV I( &G ); FFVar t = G.add_var( "t" ), x = G.add_var( "x(t)" ), u = G.add_var( "u(t)" ); FFPartial OpP;
    I.add_domain( t, FFDom( std::vector<double>{ 0., 0.5, 1. }, FFDom::LGR, 3 ) ); I.set_evolution_domain( t );
    I.add_state( x, {t} ); I.add_input( u, {t}, FFDom::LGR, 1 ); I.update_ref( x, 0. );
    int const T_INT = FFDom::ALL - FFDom::LB;
    I.add_equation( OpP( x, t ) - u, {t}, {T_INT}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
    I.add_equation( x, {t}, {FFDom::LB}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INITIAL, 0 ) );
    outputs( I, t, x, u ); I.options.SOLVE.MARCHING = march; I.options.SOLVE.RES_TOL = 1e-13; I.options.DISPLAY_LEVEL = 0;
    bool ok = I.setup(); std::vector<double> xv, inp; ok = ok && I.init( xv, inp, nullptr );
    std::vector<double> uv{ U0, U1 }; I.set_input_values( u, uv, inp.data() ); I.register_control( u );
    std::vector<double> x1( xv ), x2( xv ), x3( xv );
    ok = ok && I.solve( x1.data(), inp.data(), nullptr ).converged; auto const F = I.val_functions();
    ok = ok && I.solve_fsens( x2.data(), inp.data(), nullptr ); auto const Jf = I.sens_jacobian();
    ok = ok && I.solve_asens( x3.data(), inp.data(), nullptr ); auto const Ja = I.sens_jacobian();
    std::vector<std::vector<double>> Gf( 9, std::vector<double>( 2 ) ), Ga( Gf );          // sens_jacobian: nf x ncd, row-major
    for( size_t k = 0; ok && k < 9 && Jf.size() >= 18; ++k ) for( size_t d = 0; d < 2; ++d ){ Gf[k][d] = Jf[k*2+d]; Ga[k][d] = Ja[k*2+d]; }
    report( march? "OCFESLV marching": "OCFESLV monolith", ok && Jf.size() == 18, F, Gf, Ga );
  }
  std::printf( "\n  SIDES_gate: %d passed, %d failed -- %s\n", npass, nfail, nfail? "FAILURES": "ALL PASS" );
  return nfail? 1: 0;
}
