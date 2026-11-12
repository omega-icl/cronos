// DAESLV_sides -- ONE-SIDED point outputs (FFDom::MINUS / FFDom::PLUS) on DAESLV, against exact answers: values,
// forward and adjoint gradients, each one-sided output as an FFEval expression and in point form.  The counterpart of
// ODESLV_SIDESgate, with what only a DAE solver meets: an ALGEBRAIC state, re-solved at a stop, and a jump that is a
// ROW of the stop system rather than an explicit map.  At a stage time tau where the integrator stops, x(tau^-) is
// the end of the stage before and x(tau^+) the solution of the stop system; an output reads the side it names, its
// forward sensitivity the same side's, and its adjoint seeds the costate on that side of the stop's transpose.
//   (A) ODESLV_SIDESgate's model: x' = u, x(0) = 0, u piecewise constant on {0, 0.5, 1} (u0 = 2, u1 = 5)
//   (B) an algebraic state at an INPUT jump: x' = -k x, x(0) = a, z = u x:  z(0.5-) = u0 x(0.5), z(0.5+) = u1 x(0.5)
//   (B2) the same algebraic state with u constant: interior values with no change of rows
//   (C) a TRANSITION, x(tau+) = x(tau-) + d, at tau = 0.5 (an element boundary) and 0.4 (inside an element), with
//       y' = x; ODESLV on the same declaration as the reference
//   (D) an algebraic state across the transition: z = x^2 jumps with x
//   (E) an IMPLICIT transition, x(0.5+)^2 = x(0.5-) (the post-jump state the positive root)
//   (G) values at t0, on a mass-matrix row and an algebraic state
//   (F) the method of lines: the zero-flux rod with the source p e^-t and u(0.3+) = u(0.3-)/2 inside element 0 of t;
//       J = int u dz on either side: J(0.3-) = 1 + p(1 - e^-0.3), J(0.3+) = J(0.3-)/2
#include <cmath>
#include <cstdio>
#include <functional>
#include <sstream>
#include <string>
#include <vector>
#include "odeslv_cvodes.hpp"
#include "daeslv_idas.hpp"
using namespace mc;

static int npass = 0, nfail = 0;
static void check( bool c, std::string const& w, double v = -1. )
{ std::printf( "  %s  %-84s", c? "PASS": "FAIL", w.c_str() ); if( v >= 0. ) std::printf( " (%.2e)", v ); std::printf( "\n" ); c? ++npass: ++nfail; }
template <class S> static bool normal( S s ){ return s == decltype( s )::NORMAL; }
int const TI = FFDom::ALL - FFDom::LB;
FFModel::EqnOptions const IO( FFModel::EqnRole::INTERIOR, 0 ), II( FFModel::EqnRole::INITIAL, 0 );

// the values and gradients, [output] and [output][input DOF], the DOFs of the listed inputs in order
struct Res { bool ok = false; std::string why; std::vector<double> F; std::vector<std::vector<double>> Gf, Ga; };
template <class SOLVER>
static Res run( std::function<std::vector<FFVar>(SOLVER&,FFGraph&)> const& declare, std::vector<double> const& pv )
{
  Res R;  FFGraph G;  SOLVER I( &G );
  std::vector<FFVar> const in = declare( I, G );
  I.options.DISPLAY_LEVEL = 0;
  I.options.SOLVE.ATOL = I.options.SOLVE.RTOL = I.options.SOLVE.ATOLS = I.options.SOLVE.RTOLS = 1e-12;
  I.options.SOLVE.ATOLB = I.options.SOLVE.RTOLB = 1e-12;
  std::ostringstream os;  auto* old = std::cerr.rdbuf( os.rdbuf() );  bool const su = I.setup();  std::cerr.rdbuf( old );
  if( !su ){ R.why = "setup: " + I.extract_error();  return R; }
  std::vector<double> P( I.np(), 0. );  std::vector<size_t> col;
  size_t q = 0;
  for( auto const& u : in ) for( size_t j : I.parameter_index( u ) ){ P[j] = pv.at( q++ );  col.push_back( j ); }
  if( !normal( I.solve_fsens( P ) ) ){ R.why = "solve_fsens: " + I.extract_error();  return R; }
  R.F = I.val_function();  auto const Gf = I.val_function_gradient();
  if( !normal( I.solve_asens( P ) ) ){ R.why = "solve_asens: " + I.extract_error();  return R; }
  auto const Ga = I.val_function_gradient();
  R.Gf.assign( R.F.size(), std::vector<double>( col.size() ) );  R.Ga = R.Gf;
  for( size_t k = 0; k < R.F.size(); ++k ) for( size_t c = 0; c < col.size(); ++c ){ R.Gf[k][c] = Gf[col[c]][k];  R.Ga[k][c] = Ga[col[c]][k]; }
  R.ok = true;  return R;
}
// central differences of a closed form F(v), [output][input]
static std::vector<std::vector<double>> fd( std::function<std::vector<double>(std::vector<double> const&)> const& f, std::vector<double> const& v )
{
  std::vector<std::vector<double>> g( f( v ).size(), std::vector<double>( v.size() ) );
  for( size_t q = 0; q < v.size(); ++q ){ auto a = v, b = v;  double const h = 1e-6*std::max( 1., std::fabs( v[q] ) );  a[q] += h;  b[q] -= h;
    auto const fp = f( a ), fm = f( b );  for( size_t k = 0; k < fp.size(); ++k ) g[k][q] = ( fp[k]-fm[k] )/( 2.*h ); }
  return g;
}
// one line per output: value, forward and adjoint gradients against the exact ones
static void report( std::string const& who, Res const& R, std::vector<std::string> const& name,
                    std::vector<double> const& val, std::vector<std::vector<double>> const& grd, double tol = 1e-7 )
{
  if( !R.ok || R.F.size() != val.size() ){ check( false, who + ": sets up and solves -- " + ( R.ok? "output count": R.why ) );  return; }
  for( size_t k = 0; k < val.size(); ++k ){
    double ev = std::fabs( R.F[k] - val[k] ), eg = 0.;
    for( size_t c = 0; c < grd[k].size(); ++c ) eg = std::max( { eg, std::fabs( R.Gf[k][c] - grd[k][c] ), std::fabs( R.Ga[k][c] - grd[k][c] ) } );
    char w[200];  std::snprintf( w, sizeof w, "%s: %-12s value %.6f, forward and adjoint gradients exact", who.c_str(), name[k].c_str(), val[k] );
    check( ev < 1e-9 && eg < tol, w, std::max( ev, eg ) );
    if( ev < 1e-9 && eg < tol ) continue;                  // on a failure, what was got: value, then per input fwd / adj / exact
    std::printf( "        got %.10f;", R.F[k] );
    for( size_t c = 0; c < grd[k].size(); ++c ) std::printf( "  [%zu] %.6f / %.6f / %.6f", c, R.Gf[k][c], R.Ga[k][c], grd[k][c] );
    std::printf( "\n" );
  }
}

int main()
{
  FFPartial OpP;  FFEval OpE;  FFIntegral OpI;
  auto const MINUS = FFDom::MINUS, PLUS = FFDom::PLUS;

  std::printf( "\n(A) x' = u, x(0) = 0, u piecewise constant on {0, 0.5, 1}: ODESLV_SIDESgate's model\n" );
  { double const U0 = 2., U1 = 5.;
    std::vector<std::string> const nm{ "u(0.5-)", "u(0.5+)", "x(0.5-)", "x(0.5+)", "u(0.3-)", "u(0.3+)", "u(1-)", "u(1+)", "int u",
                                       "u(0.5-) pt", "u(0.5+) pt", "x(0.5-) pt", "x(0.5+) pt" };
    std::vector<double> const val{ 2., 5., 1., 1., 2., 2., 5., 5., 3.5,  2., 5., 1., 1. };
    std::vector<std::vector<double>> const grd{ {1,0}, {0,1}, {.5,0}, {.5,0}, {1,0}, {1,0}, {0,1}, {0,1}, {.5,.5},  {1,0}, {0,1}, {.5,0}, {.5,0} };
    auto const R = run<DAESLV_IDAS>( [&]( DAESLV_IDAS& I, FFGraph& G ){
      FFVar t = G.add_var( "t" ), x = G.add_var( "x(t)" ), u = G.add_var( "u(t)" );
      I.add_domain( t, FFDom( std::vector<double>{ 0., 0.5, 1. }, FFDom::LGR, 4 ) );  I.set_evolution_domain( t );
      I.add_state( x, {t} );  I.add_input( u, {t}, FFDom::LGR, 1 );  I.update_ref( x, 0. );
      I.add_equation( OpP( x, t ) - u, {t}, {TI}, IO );  I.add_equation( x, {t}, {FFDom::LB}, II );
      I.add_output( OpE( u, t, 0.5, MINUS ) );  I.add_output( OpE( u, t, 0.5, PLUS ) );
      I.add_output( OpE( x, t, 0.5, MINUS ) );  I.add_output( OpE( x, t, 0.5, PLUS ) );
      I.add_output( OpE( u, t, 0.3, MINUS ) );  I.add_output( OpE( u, t, 0.3, PLUS ) );
      I.add_output( OpE( u, t, 1.0, MINUS ) );  I.add_output( OpE( u, t, 1.0, PLUS ) );
      I.add_output( OpI( std::vector<FFVar>{ u }, { t } )[0] );
      I.add_output( u, {t}, {0.5}, {MINUS} );   I.add_output( u, {t}, {0.5}, {PLUS} );
      I.add_output( x, {t}, {0.5}, {MINUS} );   I.add_output( x, {t}, {0.5}, {PLUS} );
      return std::vector<FFVar>{ u }; }, { U0, U1 } );
    report( "(A)", R, nm, val, grd ); }

  std::printf( "\n(B) an ALGEBRAIC state at an input jump: x' = -k x, x(0) = a, z = u x, u piecewise constant (u0, u1)\n" );
  { std::vector<double> const v{ 1.3, 2.1, 0.8, 1.7 };                       // k, a, u0, u1
    auto f = []( std::vector<double> const& w ){ double const xh = w[1]*std::exp( -w[0]*0.5 ), xq = w[1]*std::exp( -w[0]*0.25 ), x1 = w[1]*std::exp( -w[0] );
      return std::vector<double>{ w[2]*xh, w[3]*xh, w[2]*xq, w[3]*x1, w[2]*xh, w[3]*xh }; };
    std::vector<std::string> const nm{ "z(0.5-)", "z(0.5+)", "z(0.25)", "z(1)", "z(0.5-) pt", "z(0.5+) pt" };
    auto const R = run<DAESLV_IDAS>( [&]( DAESLV_IDAS& I, FFGraph& G ){
      FFVar t = G.add_var( "t" ), x = G.add_var( "x(t)" ), z = G.add_var( "z(t)" ), u = G.add_var( "u(t)" ), k = G.add_var( "k" ), a = G.add_var( "a" );
      I.add_domain( t, FFDom( 0., 1., 2, FFDom::LGR, 4 ) );  I.set_evolution_domain( t );
      I.add_state( x, {t} );  I.add_state( z, {t} );  I.add_input( k );  I.add_input( a );  I.add_input( u, {t}, FFDom::LGR, 1 );
      I.update_ref( x, 1. );  I.update_ref( z, 1. );
      I.add_equation( OpP( x, t ) + k*x, {t}, {TI}, IO );  I.add_equation( z - u*x, {t}, {FFDom::ALL}, IO );  I.add_equation( x - a, {t}, {FFDom::LB}, II );
      I.add_output( OpE( z, t, 0.5, MINUS ) );  I.add_output( OpE( z, t, 0.5, PLUS ) );  I.add_output( OpE( z, t, 0.25 ) );  I.add_output( OpE( z, t, 1. ) );
      I.add_output( z, {t}, {0.5}, {MINUS} );   I.add_output( z, {t}, {0.5}, {PLUS} );
      return std::vector<FFVar>{ k, a, u }; }, v );
    report( "(B)", R, nm, f( v ), fd( f, v ) ); }

  std::printf( "\n(B2) the same with u CONSTANT: no stop changes the rows, the interior values still read the algebraic row\n" );
  { std::vector<double> const v{ 1.3, 2.1, 0.8 };                            // k, a, u
    auto f = []( std::vector<double> const& w ){
      return std::vector<double>{ w[2]*w[1]*std::exp( -w[0]*0.25 ), w[2]*w[1]*std::exp( -w[0]*0.5 ), w[2]*w[1]*std::exp( -w[0] ) }; };
    std::vector<std::string> const nm{ "z(0.25)", "z(0.5)", "z(1)" };
    auto const R = run<DAESLV_IDAS>( [&]( DAESLV_IDAS& I, FFGraph& G ){
      FFVar t = G.add_var( "t" ), x = G.add_var( "x(t)" ), z = G.add_var( "z(t)" ), u = G.add_var( "u" ), k = G.add_var( "k" ), a = G.add_var( "a" );
      I.add_domain( t, FFDom( 0., 1., 2, FFDom::LGR, 4 ) );  I.set_evolution_domain( t );
      I.add_state( x, {t} );  I.add_state( z, {t} );  I.add_input( k );  I.add_input( a );  I.add_input( u );
      I.update_ref( x, 1. );  I.update_ref( z, 1. );
      I.add_equation( OpP( x, t ) + k*x, {t}, {TI}, IO );  I.add_equation( z - u*x, {t}, {FFDom::ALL}, IO );  I.add_equation( x - a, {t}, {FFDom::LB}, II );
      I.add_output( OpE( z, t, 0.25 ) );  I.add_output( OpE( z, t, 0.5 ) );  I.add_output( OpE( z, t, 1. ) );
      return std::vector<FFVar>{ k, a, u }; }, v );
    report( "(B2)", R, nm, f( v ), fd( f, v ) ); }

  std::printf( "\n(C) a TRANSITION x(tau+) = x(tau-) + d: x' = -k x, y' = x, x(0) = x0, on and off the element boundaries\n" );
  for( double const tau : { 0.5, 0.4 } ){
    std::vector<double> const v{ 1.3, 0.7, 1.1 };                            // k, d, x0
    auto f = [tau]( std::vector<double> const& w ){ double const xm = w[2]*std::exp( -w[0]*tau ), xp = xm + w[1], xf = xp*std::exp( -w[0]*( 1.-tau ) );
      double const I = w[2]*( 1.-std::exp( -w[0]*tau ) )/w[0] + xp*( 1.-std::exp( -w[0]*( 1.-tau ) ) )/w[0];
      return std::vector<double>{ xm, xp, xf, I, xm, xp }; };
    std::vector<std::string> const nm{ "x(tau-)", "x(tau+)", "x(1)", "y(1)", "x(tau-) pt", "x(tau+) pt" };
    auto declare = [&]( FFModel& I, FFGraph& G ){
      FFVar t = G.add_var( "t" ), x = G.add_var( "x(t)" ), y = G.add_var( "y(t)" ), k = G.add_var( "k" ), d = G.add_var( "d" ), x0 = G.add_var( "x0" );
      I.add_domain( t, FFDom( 0., 1., 2, FFDom::LGR, 4 ) );  I.set_evolution_domain( t );
      I.add_state( x, {t} );  I.add_state( y, {t} );  I.add_input( k );  I.add_input( d );  I.add_input( x0 );  I.update_ref( x, 1. );  I.update_ref( y, 0. );
      I.add_equation( OpP( x, t ) + k*x, {t}, {TI}, IO );  I.add_equation( OpP( y, t ) - x, {t}, {TI}, IO );
      I.add_equation( x - x0, {t}, {FFDom::LB}, II );       I.add_equation( y, {t}, {FFDom::LB}, II );
      I.add_transition( x + d, x, t, tau );
      I.add_output( OpE( x, t, tau, MINUS ) );  I.add_output( OpE( x, t, tau, PLUS ) );  I.add_output( OpE( x, t, 1. ) );  I.add_output( OpE( y, t, 1. ) );
      I.add_output( x, {t}, {tau}, {MINUS} );   I.add_output( x, {t}, {tau}, {PLUS} );
      return std::vector<FFVar>{ k, d, x0 }; };
    std::string const tag = tau == 0.5? "(C) tau = 0.5": "(C) tau = 0.4";
    report( tag + " ODESLV", run<ODESLVS_CVODES>( [&]( ODESLVS_CVODES& I, FFGraph& G ){ return declare( I, G ); }, v ), nm, f( v ), fd( f, v ) );
    report( tag + " DAESLV", run<DAESLV_IDAS>(    [&]( DAESLV_IDAS& I,    FFGraph& G ){ return declare( I, G ); }, v ), nm, f( v ), fd( f, v ) ); }

  std::printf( "\n(D) an ALGEBRAIC state across the transition: z = x^2, x(0.4+) = x(0.4-) + d\n" );
  { std::vector<double> const v{ 1.3, 0.7, 1.1 };
    auto f = []( std::vector<double> const& w ){ double const xm = w[2]*std::exp( -w[0]*0.4 ), xp = xm + w[1], xf = xp*std::exp( -w[0]*0.6 );
      return std::vector<double>{ xm*xm, xp*xp, xf*xf }; };
    std::vector<std::string> const nm{ "z(0.4-)", "z(0.4+)", "z(1)" };
    auto const R = run<DAESLV_IDAS>( [&]( DAESLV_IDAS& I, FFGraph& G ){
      FFVar t = G.add_var( "t" ), x = G.add_var( "x(t)" ), z = G.add_var( "z(t)" ), k = G.add_var( "k" ), d = G.add_var( "d" ), x0 = G.add_var( "x0" );
      I.add_domain( t, FFDom( 0., 1., 2, FFDom::LGR, 4 ) );  I.set_evolution_domain( t );
      I.add_state( x, {t} );  I.add_state( z, {t} );  I.add_input( k );  I.add_input( d );  I.add_input( x0 );  I.update_ref( x, 1. );  I.update_ref( z, 1. );
      I.add_equation( OpP( x, t ) + k*x, {t}, {TI}, IO );  I.add_equation( z - x*x, {t}, {FFDom::ALL}, IO );  I.add_equation( x - x0, {t}, {FFDom::LB}, II );
      I.add_transition( std::vector<FFVar>{ x + d }, std::vector<FFVar>{ x }, t, 0.4, std::vector<FFVar>{ x } );
      I.add_output( OpE( z, t, 0.4, MINUS ) );  I.add_output( OpE( z, t, 0.4, PLUS ) );  I.add_output( OpE( z, t, 1. ) );
      return std::vector<FFVar>{ k, d, x0 }; }, v );
    report( "(D)", R, nm, f( v ), fd( f, v ) ); }

  std::printf( "\n(E) an IMPLICIT transition x(0.5+)^2 = x(0.5-): x' = -k x, x(0) = a\n" );
  { std::vector<double> const v{ 1.3, 2.1 };                                 // k, a
    auto f = []( std::vector<double> const& w ){ double const xm = w[1]*std::exp( -w[0]*0.5 ), xp = std::sqrt( xm ), xf = xp*std::exp( -w[0]*0.5 );
      return std::vector<double>{ xm, xp, xf }; };
    std::vector<std::string> const nm{ "x(0.5-)", "x(0.5+)", "x(1)" };
    auto const R = run<DAESLV_IDAS>( [&]( DAESLV_IDAS& I, FFGraph& G ){
      FFVar t = G.add_var( "t" ), x = G.add_var( "x(t)" ), k = G.add_var( "k" ), a = G.add_var( "a" );
      I.add_domain( t, FFDom( 0., 1., 2, FFDom::LGR, 4 ) );  I.set_evolution_domain( t );
      I.add_state( x, {t} );  I.add_input( k );  I.add_input( a );  I.update_ref( x, 1. );
      I.add_equation( OpP( x, t ) + k*x, {t}, {TI}, IO );  I.add_equation( x - a, {t}, {FFDom::LB}, II );
      I.add_transition( x, x*x, t, 0.5 );                                    // left(tau-) = right(tau+)
      I.add_output( OpE( x, t, 0.5, MINUS ) );  I.add_output( OpE( x, t, 0.5, PLUS ) );  I.add_output( OpE( x, t, 1. ) );
      return std::vector<FFVar>{ k, a }; }, v );
    report( "(E)", R, nm, f( v ), fd( f, v ) ); }

  std::printf( "\n(G) values at t0 on a MASS-MATRIX row and an algebraic state: 2 x' = -k x, x(0) = a, z = u x\n" );
  { std::vector<double> const v{ 1.3, 2.1, 0.8 };                            // k, a, u
    auto f = []( std::vector<double> const& w ){ return std::vector<double>{ w[1], w[2]*w[1], w[1]*std::exp( -0.5*w[0] ) }; };
    std::vector<std::string> const nm{ "x(0)", "z(0)", "x(1)" };
    auto const R = run<DAESLV_IDAS>( [&]( DAESLV_IDAS& I, FFGraph& G ){
      FFVar t = G.add_var( "t" ), x = G.add_var( "x(t)" ), z = G.add_var( "z(t)" ), u = G.add_var( "u" ), k = G.add_var( "k" ), a = G.add_var( "a" );
      I.add_domain( t, FFDom( 0., 1., 2, FFDom::LGR, 4 ) );  I.set_evolution_domain( t );
      I.add_state( x, {t} );  I.add_state( z, {t} );  I.add_input( k );  I.add_input( a );  I.add_input( u );
      I.update_ref( x, 1. );  I.update_ref( z, 1. );
      I.add_equation( 2.*OpP( x, t ) + k*x, {t}, {TI}, IO );  I.add_equation( z - u*x, {t}, {FFDom::ALL}, IO );  I.add_equation( x - a, {t}, {FFDom::LB}, II );
      I.add_output( OpE( x, t, 0. ) );  I.add_output( OpE( z, t, 0. ) );  I.add_output( OpE( x, t, 1. ) );
      return std::vector<FFVar>{ k, a, u }; }, v );
    report( "(G)", R, nm, f( v ), fd( f, v ) ); }

  std::printf( "\n(F) the method of lines: the zero-flux rod with the source p e^-t, u(0.3+) = u(0.3-)/2 inside element 0 of t\n" );
  { double const P = 0.7, TAU = 0.3, e = std::exp( -TAU );
    std::vector<double> const val{ 1. + P*( 1.-e ), 0.5*( 1. + P*( 1.-e ) ) };
    std::vector<std::vector<double>> const grd{ { 1.-e }, { 0.5*( 1.-e ) } };
    Res R;  FFGraph G;  DAESLV_IDAS S( &G );
    FFVar t = G.add_var("t"), z = G.add_var("z"), u = G.add_var("u"), p = G.add_var("p");
    S.add_domain( t, FFDom( 0., 1., 2, FFDom::LGR, 3 ) );  S.add_domain( z, FFDom( 0., 1., 3, FFDom::LGL, 5 ) );
    S.set_evolution_domain( t );  S.add_state( u, {t,z} );  S.update_ref( u, 1. );  S.add_input( p, {} );
    int const ZI = FFDom::ALL - FFDom::LB - FFDom::UB;
    FFModel::EqnOptions const ib( FFModel::EqnRole::BOUNDARY, 0 );
    S.add_equation( OpP(u,t) - OpP( 0.1*OpP(u,z), z ) - p*exp( -t ), {t,z}, {TI, ZI},                IO );
    S.add_equation( OpP(u,z),                                       {t,z}, {TI, FFDom::LB},         ib );
    S.add_equation( OpP(u,z),                                       {t,z}, {TI, FFDom::UB},         ib );
    S.add_equation( u - 1.,                                         {t,z}, {FFDom::LB, FFDom::ALL}, II );
    S.add_transition( { 0.5*u }, { u }, t, TAU );
    S.add_output( OpI( u, z ), {t}, {TAU}, {MINUS} );  S.add_output( OpI( u, z ), {t}, {TAU}, {PLUS} );
    S.options.INTERFACE.IMPOSITION = OCFEMOL::Options::IC_STRONG;  S.options.DISPLAY_LEVEL = 0;  S.options.SOLVE.VERBOSE = 0;
    S.options.SOLVE.RTOL = S.options.SOLVE.RTOLS = S.options.SOLVE.RTOLB = 1e-10;  S.options.SOLVE.ATOL = S.options.SOLVE.ATOLS = S.options.SOLVE.ATOLB = 1e-12;
    std::ostringstream os;  auto* old = std::cerr.rdbuf( os.rdbuf() );  bool const su = S.setup();  std::cerr.rdbuf( old );
    std::vector<FFModel::InputVal> const in{ { p, { P } } };
    if( !su ) R.why = "setup: " + S.extract_error();
    else if( !normal( S.solve_fsens( in ) ) ) R.why = "solve_fsens: " + S.extract_error();
    else{ R.F = S.val_function();  auto const Gf = S.val_function_gradient();
      if( !normal( S.solve_asens( in ) ) ) R.why = "solve_asens: " + S.extract_error();
      else{ auto const Ga = S.val_function_gradient();  R.Gf = { { Gf[0][0] }, { Gf[0][1] } };  R.Ga = { { Ga[0][0] }, { Ga[0][1] } };  R.ok = true; } }
    report( "(F)", R, { "J(0.3-)", "J(0.3+)" }, val, grd, 1e-8 ); }

  std::printf( "\n  DAESLV_sides: %d passed, %d failed -- %s\n", npass, nfail, nfail? "FAILURES": "ALL PASS" );
  return nfail? 1: 0;
}
