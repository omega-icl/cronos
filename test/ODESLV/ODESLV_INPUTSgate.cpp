// INPUTS_gate.cpp -- outputs that EVALUATE and INTEGRATE time-varying inputs, ODESLV and OCFESLV (monolithic and
// marching), against a closed form.   dx/dt = -a x + u(t), x(0) = 1;  u piecewise LINEAR (LGL, 2 nodes per element,
// DISCONTINUOUS at element boundaries) on the NON-uniform mesh {0, 0.3, 0.6, 1}.
//   F0 = u(0.45)                (inside an element)          F1 = int u^2 dt
//   F2 = int x u dt             (state and input)            F3 = x(0.6) + 2 u(0.6+)   (affine: ODESLV takes it)
//   F4 = u(0.6-) - u(0.6+)      (the jump, sides explicit)
// Oracle: closed-form x and u, integrals by 20-point Gauss-Legendre per element; gradients w.r.t. a and the 6 levels
// by central differences OF THE ORACLE.  Input-only outputs exact everywhere; state outputs: ODESLV to 1e-8, OCFESLV
// to discretisation accuracy.
// F1 guards a FIXED OCFESLV defect (2026-09-28): a nonlinear function of a multi-node input used to be formed at the
// input's own nodes and interpolated linearly -- a trapezoidal rule, ~30% off here.  Nonlinear OCVar operations now
// lift such operands onto the collocation grid first (OCVar::_nl_lift); OPS_gate checks every operation.
#include <cmath>
#include <cstdio>
#include <vector>
#include "ffode.hpp"
#include "ocfeslv.hpp"
using namespace mc;
static int npass = 0, nfail = 0;
static void check( bool c, char const* w, double v ){ std::printf( "  %s  %-66s (%.2e)\n", c? "PASS": "FAIL", w, v ); c? ++npass: ++nfail; }
std::vector<double> const EB{ 0., 0.3, 0.6, 1. };
size_t const NE = 3, NN = 2;
// the oracle: P = { a, u[e][j] element-major }
static std::vector<double> oracle( std::vector<double> const& P ){
  double const a = P[0];
  auto lv = [&]( size_t e, size_t j ){ return P[1 + e*NN + j]; };
  auto uu = [&]( size_t e, double t ){ double const s = ( t - EB[e] )/( EB[e+1]-EB[e] ); return lv(e,0)*( 1.-s ) + lv(e,1)*s; };
  // x on element e from x(te): u = al + be (t - te)
  std::vector<double> x0( NE+1 ); x0[0] = 1.;
  auto xx = [&]( size_t e, double t ){ double const te = EB[e], he = EB[e+1]-te, al = lv(e,0), be = ( lv(e,1)-lv(e,0) )/he;
    double const xp0 = al/a - be/(a*a); return xp0 + (be/a)*(t-te) + ( x0[e] - xp0 )*std::exp( -a*(t-te) ); };
  for( size_t e = 0; e < NE; ++e ) x0[e+1] = xx( e, EB[e+1] );
  // 20-point Gauss-Legendre on [-1,1]
  static double const gx[10] = { 0.0765265211334973, 0.2277858511416451, 0.3737060887154195, 0.5108670019508271, 0.6360536807265150,
                                 0.7463319064601508, 0.8391169718222188, 0.9122344282513259, 0.9639719272779138, 0.9931285991850949 };
  static double const gw[10] = { 0.1527533871307258, 0.1491729864726037, 0.1420961093183820, 0.1316886384491766, 0.1181945319615184,
                                 0.1019301198172404, 0.0832767415767048, 0.0626720483341091, 0.0406014298003869, 0.0176140071391521 };
  double I1 = 0., I2 = 0.;
  for( size_t e = 0; e < NE; ++e ) for( int k = 0; k < 10; ++k ) for( int sg : { -1, 1 } ){
    double const h = EB[e+1]-EB[e], t = EB[e] + 0.5*h*( 1. + sg*gx[k] ), w = 0.5*h*gw[k];
    I1 += w * uu(e,t)*uu(e,t);  I2 += w * xx(e,t)*uu(e,t); }
  return { uu( 1, 0.45 ), I1, I2, x0[2] + 2.*uu( 2, 0.6 ), uu( 1, 0.6 ) - uu( 2, 0.6 ) };
}
static std::vector<std::vector<double>> oracle_grad( std::vector<double> const& P ){     // [param][output]
  std::vector<std::vector<double>> D( P.size() );
  for( size_t i = 0; i < P.size(); ++i ){ auto Pp = P, Pm = P; double const h = 1e-6*std::max( 1., std::fabs( P[i] ) ); Pp[i] += h; Pm[i] -= h;
    auto const fp = oracle( Pp ), fm = oracle( Pm ); D[i].resize( fp.size() ); for( size_t k = 0; k < fp.size(); ++k ) D[i][k] = ( fp[k]-fm[k] )/( 2.*h ); }
  return D; }
template <class M> static void outputs( M& I, FFVar const& t, FFVar const& x, FFVar const& u ){
  FFEval OpE; FFIntegral OpI;
  I.add_output( OpE( u, t, 0.45 ) );
  I.add_output( OpI( std::vector<FFVar>{ u*u }, { t } )[0] );
  I.add_output( OpI( std::vector<FFVar>{ x*u }, { t } )[0] );
  I.add_output( OpE( x, t, 0.6 ) + 2.*OpE( u, t, 0.6, FFDom::PLUS ) );
  I.add_output( OpE( u, t, 0.6, FFDom::MINUS ) - OpE( u, t, 0.6, FFDom::PLUS ) ); }
static char const* const NM[5] = { "u(0.45)", "int u^2", "int x u", "x(0.6) + 2u(0.6+)", "u(0.6-) - u(0.6+)" };
static void report( char const* who, bool ok, std::vector<double> const& F, std::vector<std::vector<double>> const& Gf,
                    std::vector<std::vector<double>> const& Ga, std::vector<double> const& P, double tolx ){
  char w[140];
  if( !ok || F.size() != 5 ){ std::snprintf( w, sizeof w, "%s: solves and 5 outputs", who ); check( false, w, 0. ); return; }
  auto const F0 = oracle( P ); auto const D = oracle_grad( P );
  for( size_t k = 0; k < 5; ++k ){
    bool const has_x = ( k == 2 || k == 3 );  double const tv = has_x? tolx: 1e-10, tg = has_x? 10.*tolx: 1e-7;
    double ev = std::fabs( F[k]-F0[k] ), eg = 0.;
    for( size_t i = 0; i < P.size(); ++i ) eg = std::max( { eg, std::fabs( Gf[i][k]-D[i][k] ), std::fabs( Ga[i][k]-D[i][k] ) } );
    std::snprintf( w, sizeof w, "%s: %-18s value", who, NM[k] ); check( ev < tv, w, ev );
    std::snprintf( w, sizeof w, "%s: %-18s forward and adjoint gradients", who, NM[k] ); check( eg < tg, w, eg );
  }
}
int main(){
  std::vector<double> P{ 1.7, 0.5, 1.2, -0.4, 0.9, 2.0, 0.3 };      // a, then u[e][j]
  { // ---- ODESLV
    FFGraph G; ODESLVS_CVODES I( &G ); FFVar t = G.add_var( "t" ), x = G.add_var( "x(t)" ), a = G.add_var( "a" ), u = G.add_var( "u(t)" ); FFPartial OpP;
    I.add_domain( t, FFDom( EB, FFDom::LGR, 4 ) ); I.set_evolution_domain( t );
    I.add_state( x, {t} ); I.add_input( a ); I.add_input( u, {t}, FFDom::LGL, NN ); I.update_ref( x, 1. );
    int const T_INT = FFDom::ALL - FFDom::LB;
    I.add_equation( OpP( x, t ) + a*x - u, {t}, {T_INT}, FFModel::EqnOptions( FFModel::EqnRole::INTERIOR, 0 ) );
    I.add_equation( x - 1., {t}, {FFDom::LB}, FFModel::EqnOptions( FFModel::EqnRole::INITIAL, 0 ) );
    outputs( I, t, x, u ); I.options.DISPLAY = 0; I.options.ATOL = I.options.ATOLS = I.options.ATOLB = 1e-12; I.options.RTOL = I.options.RTOLS = I.options.RTOLB = 1e-12;
    bool ok = I.setup(); std::vector<double> Q( I.np() ); auto const ia = I.parameter_index( a )[0]; auto const iu = I.parameter_index( u );
    Q[ia] = P[0]; for( size_t k = 0; k < iu.size(); ++k ) Q[iu[k]] = P[1+k];
    ok = ok && I.solve_fsens( Q ) == ODESLVS_CVODES::STATUS::NORMAL; auto const F = I.val_function(); auto const G0 = I.val_function_gradient();
    ok = ok && I.solve_asens( Q ) == ODESLVS_CVODES::STATUS::NORMAL; auto const G1 = I.val_function_gradient();
    std::vector<std::vector<double>> Gf( P.size(), std::vector<double>( 5 ) ), Ga( Gf );
    for( size_t k = 0; ok && k < 5; ++k ){ Gf[0][k] = G0[ia][k]; Ga[0][k] = G1[ia][k]; for( size_t j = 0; j < iu.size(); ++j ){ Gf[1+j][k] = G0[iu[j]][k]; Ga[1+j][k] = G1[iu[j]][k]; } }
    report( "ODESLV          ", ok && iu.size() == 6, F, Gf, Ga, P, 1e-8 );
  }
  for( bool march : { false, true } ){ // ---- OCFESLV
    FFGraph G; OCFESLV I( &G ); FFVar t = G.add_var( "t" ), x = G.add_var( "x(t)" ), a = G.add_var( "a" ), u = G.add_var( "u(t)" ); FFPartial OpP;
    I.add_domain( t, FFDom( EB, FFDom::LGR, 7 ) ); I.set_evolution_domain( t );
    I.add_state( x, {t} ); I.add_input( a, {} ); I.add_input( u, {t}, FFDom::LGL, NN ); I.update_ref( x, 1. );
    int const T_INT = FFDom::ALL - FFDom::LB;
    I.add_equation( OpP( x, t ) + a*x - u, {t}, {T_INT}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
    I.add_equation( x - 1., {t}, {FFDom::LB}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INITIAL, 0 ) );
    outputs( I, t, x, u ); I.options.SOLVE.MARCHING = march; I.options.SOLVE.RES_TOL = 1e-13; I.options.DISPLAY_LEVEL = 0;
    bool ok = I.setup(); std::vector<double> xv, inp; ok = ok && I.init( xv, inp, nullptr );
    I.set_input_values( a, { P[0] }, inp.data() ); I.set_input_values( u, std::vector<double>( P.begin()+1, P.end() ), inp.data() );
    I.register_control( a ); I.register_control( u );
    std::vector<double> x1( xv ), x2( xv ), x3( xv );
    ok = ok && I.solve( x1.data(), inp.data(), nullptr ).converged; auto const F = I.val_functions();
    ok = ok && I.solve_fsens( x2.data(), inp.data(), nullptr ); auto const Jf = I.sens_jacobian();
    ok = ok && I.solve_asens( x3.data(), inp.data(), nullptr ); auto const Ja = I.sens_jacobian();
    size_t const nc = P.size();  ok = ok && Jf.size() == 5*nc;
    std::vector<std::vector<double>> Gf( nc, std::vector<double>( 5 ) ), Ga( Gf );          // sens_jacobian: nf x ncd, row-major
    for( size_t k = 0; ok && k < 5; ++k ) for( size_t i = 0; i < nc; ++i ){ Gf[i][k] = Jf[k*nc+i]; Ga[i][k] = Ja[k*nc+i]; }
    report( march? "OCFESLV marching": "OCFESLV monolith", ok, F, Gf, Ga, P, 1e-7 );
  }
  std::printf( "\n  INPUTS_gate: %d passed, %d failed -- %s\n", npass, nfail, nfail? "FAILURES": "ALL PASS" );
  return nfail? 1: 0;
}
