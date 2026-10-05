// ODESLV_nonuniform.cpp -- ODESLV on a NON-UNIFORM evolution mesh, FFDom( elem_bnd, ... ), against a CLOSED FORM.
//   dx/dt = -a x + u(t),  x(0) = 1,  u piecewise constant (LGR, 1 node) or piecewise linear (LGL, 2 nodes) per element
//   outputs: x(0.3) (inside an element: an inserted stage), x(0.45), x(tf)
//   Checks, uniform and non-uniform meshes: values, FSA and ASA gradients (w.r.t. a and every control level) against
//   the closed form (gradients by central differences OF THE CLOSED FORM).
#include <cmath>
#include <cstdio>
#include <vector>
#include "ffode.hpp"

static int npass = 0, nfail = 0;
static void check( bool c, char const* w, double v ){ std::printf( "  %s  %-62s (%.3e)\n", c? "PASS": "FAIL", w, v ); c? ++npass: ++nfail; }
std::vector<double> const TOUT{ 0.3, 0.45, 1.0 };

//! closed form: x at each time in TOUT; lv = element-major levels (nn per element), nodes at the element ends for nn=2
static std::vector<double> exact( double a, std::vector<double> const& lv, std::vector<double> const& eb, size_t nn ){
  std::vector<double> out; double x = 1.;
  auto advance = [&]( double te, double he, double al, double be, double t ){      // u = al + be (t - te) on [te, te+he]
    double const xp0 = al/a - be/(a*a);                                               // xp(t) = xp0 + (be/a)(t - te)
    return xp0 + (be/a)*(t-te) + ( x - xp0 ) * std::exp( -a*(t-te) ); };
  size_t io = 0;
  for( size_t e = 0; e+1 < eb.size(); ++e ){
    double const te = eb[e], he = eb[e+1]-eb[e];
    double const al = lv[e*nn], be = nn > 1? ( lv[e*nn+1] - lv[e*nn] ) / he: 0.;
    while( io < TOUT.size() && TOUT[io] <= eb[e+1] + 1e-14 ){ out.push_back( advance( te, he, al, be, TOUT[io] ) ); ++io; }
    x = advance( te, he, al, be, eb[e+1] );
  }
  return out;
}

int main(){
  struct Mesh { char const* name; std::vector<double> eb; };
  for( Mesh M : { Mesh{ "uniform    ", { 0., 0.25, 0.5, 0.75, 1. } }, Mesh{ "NON-uniform", { 0., 0.1, 0.45, 0.5, 1. } } } )
  for( size_t nn : { 1, 2 } ){
    mc::FFGraph G; mc::ODESLVS_CVODES I( &G );
    mc::FFVar t = G.add_var( "t" ), x = G.add_var( "x(t)" ), a = G.add_var( "a" ), u = G.add_var( "u(t)" );
    mc::FFPartial OpP; mc::FFEval OpE;
    I.add_domain( t, mc::FFDom( M.eb, mc::FFDom::LGR, 4 ) );  I.set_evolution_domain( t );
    I.add_state( x, {t} );  I.add_input( a );  I.add_input( u, {t}, nn == 1? mc::FFDom::LGR: mc::FFDom::LGL, nn );
    I.update_ref( x, 1. );
    int const T_INT = mc::FFDom::ALL - mc::FFDom::LB;
    I.add_equation( OpP( x, t ) + a*x - u, {t}, {T_INT},         mc::FFModel::EqnOptions( mc::FFModel::EqnRole::INTERIOR, 0 ) );
    I.add_equation( x - 1.,               {t}, {mc::FFDom::LB}, mc::FFModel::EqnOptions( mc::FFModel::EqnRole::INITIAL,  0 ) );
    for( double to : TOUT ) I.add_output( OpE( x, t, to ) );
    I.options.LINSOL = mc::BASE_CVODES::Options::DENSE; I.options.DISPLAY = 0; I.options.NMAX = 100000;
    I.options.ATOL = I.options.ATOLS = I.options.ATOLB = 1e-12; I.options.RTOL = I.options.RTOLS = I.options.RTOLB = 1e-12;
    if( !I.setup() ){ std::printf( "  setup failed: %s\n", I.extract_error().c_str() ); ++nfail; continue; }
    size_t const ne = M.eb.size()-1;
    std::vector<double> lv( ne*nn ); for( size_t k = 0; k < lv.size(); ++k ) lv[k] = 0.8*std::sin( 1.3*k + 0.4 ) + 0.2;
    double const av = 1.7;
    std::vector<double> P( I.np() );  P[ I.parameter_index( a )[0] ] = av;
    auto const iu = I.parameter_index( u );  for( size_t k = 0; k < iu.size(); ++k ) P[ iu[k] ] = lv[k];
    // oracle: value and d/dP by central differences of the closed form, in parameter order
    auto closed = [&]( std::vector<double> const& Q ){ std::vector<double> l( lv.size() ); for( size_t k = 0; k < l.size(); ++k ) l[k] = Q[ iu[k] ];
                                                     return exact( Q[ I.parameter_index( a )[0] ], l, M.eb, nn ); };
    auto const F0 = closed( P );
    std::vector<std::vector<double>> D( P.size() );
    for( size_t i = 0; i < P.size(); ++i ){ auto Qp = P, Qm = P; double const h = 1e-6*std::max( 1., std::fabs( P[i] ) ); Qp[i] += h; Qm[i] -= h;
      auto const fp = closed( Qp ), fm = closed( Qm ); D[i].resize( fp.size() ); for( size_t k = 0; k < fp.size(); ++k ) D[i][k] = ( fp[k]-fm[k] )/( 2.*h ); }
    bool const s1 = I.solve_fsens( P ) == mc::ODESLVS_CVODES::STATUS::NORMAL;  auto const F = I.val_function();  auto const Gf = I.val_function_gradient();
    bool const s2 = I.solve_asens( P )     == mc::ODESLVS_CVODES::STATUS::NORMAL;  auto const Ga = I.val_function_gradient();
    double wv = 0., wf = 0., wa = 0.;
    for( size_t k = 0; k < F0.size() && k < F.size(); ++k ) wv = std::max( wv, std::fabs( F[k]-F0[k] ) );
    for( size_t i = 0; i < D.size() && s1 && s2; ++i ) for( size_t k = 0; k < D[i].size(); ++k ){
      wf = std::max( wf, std::fabs( Gf[i][k]-D[i][k] ) ); wa = std::max( wa, std::fabs( Ga[i][k]-D[i][k] ) ); }
    char w[120];
    std::snprintf( w, sizeof w, "%s, %s: values == closed form", M.name, nn == 1? "piecewise-constant u": "piecewise-linear u  " ); check( s1 && s2 && F.size() == 3 && wv < 1e-8, w, wv );
    std::snprintf( w, sizeof w, "%s, %s: FSA gradient == closed form", M.name, nn == 1? "piecewise-constant u": "piecewise-linear u  " ); check( s1 && wf < 1e-6, w, wf );
    std::snprintf( w, sizeof w, "%s, %s: ASA gradient == closed form", M.name, nn == 1? "piecewise-constant u": "piecewise-linear u  " ); check( s2 && wa < 1e-6, w, wa );
  }
  std::printf( "\n  ODESLV_nonuniform: %d passed, %d failed -- %s\n", npass, nfail, nfail? "FAILURES": "ALL PASS" );
  return nfail? 1: 0;
}
