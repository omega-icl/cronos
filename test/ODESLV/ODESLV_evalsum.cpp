// Combined outputs vs the same combination of SEPARATE outputs (which already work): value, forward, adjoint.
#include <cmath>
#include <iostream>
#include <iomanip>
#include <functional>
#include "ffode.hpp"
double const t0 = 0., tf = 10., T1 = 5., P = 3.;
struct M { mc::FFGraph G; mc::ODESLVS_CVODES I{ &G }; mc::FFVar t, x0, x1, p; };
template <class F> static bool build( M& m, F outputs, size_t nstage = 2 ){
  m.t = m.G.add_var( "t" ); m.x0 = m.G.add_var( "x0(t)" ); m.x1 = m.G.add_var( "x1(t)" ); m.p = m.G.add_var( "p" );
  mc::FFPartial OpP;
  m.I.add_domain( m.t, mc::FFDom( t0, tf, nstage, mc::FFDom::LGR, 4 ) ); m.I.set_evolution_domain( m.t );
  m.I.add_state( m.x0, {m.t} ); m.I.add_state( m.x1, {m.t} ); m.I.add_input( m.p ); m.I.update_ref( m.x0, 1.2 ); m.I.update_ref( m.x1, 1.1 );
  int const T_INT = mc::FFDom::ALL - mc::FFDom::LB;
  auto const io = mc::FFModel::EqnOptions( mc::FFModel::EqnRole::INTERIOR, 0 ), ii = mc::FFModel::EqnOptions( mc::FFModel::EqnRole::INITIAL, 0 );
  m.I.add_equation( OpP( m.x0, m.t ) - m.p*m.x0*(1.-m.x1), {m.t}, {T_INT}, io ); m.I.add_equation( OpP( m.x1, m.t ) - m.p*m.x1*(m.x0-1.), {m.t}, {T_INT}, io );
  m.I.add_equation( m.x0 - 1.2, {m.t}, {mc::FFDom::LB}, ii ); m.I.add_equation( m.x1 - ( 1.1 + 0.01*m.p ), {m.t}, {mc::FFDom::LB}, ii );
  outputs( m );
  m.I.options.DISPLAY = 0; m.I.options.LINSOL = mc::BASE_CVODES::Options::DENSE;
  m.I.options.ATOL = m.I.options.ATOLS = m.I.options.ATOLB = 1e-10; m.I.options.RTOL = m.I.options.RTOLS = m.I.options.RTOLB = 1e-10;
  return m.I.setup(); }
struct R { std::vector<double> F, Gf, Ga; bool ok = false; };
static R solve( M& m ){
  R r; std::vector<mc::FFModel::InputVal> v{ { m.p, { P } } };
  r.ok = m.I.solve_fsens( v ) == 0; r.F = m.I.val_function(); for( double g : m.I.val_function_gradient()[0] ) r.Gf.push_back( g );
  r.ok = r.ok && m.I.solve_asens( v ) == 0;     for( double g : m.I.val_function_gradient()[0] ) r.Ga.push_back( g ); return r; }
int main(){
  std::cout << std::scientific << std::setprecision(6);
  mc::FFEval OpE; mc::FFIntegral OpI;
  // reference: e0 = x0^2(T1), e1 = x0^2(tf), e2 = x1(tf), q = int x1 dt   -- four SEPARATE outputs
  M ref; build( ref, [&]( M& m ){ m.I.add_output( OpE( sqr( m.x0 ), m.t, T1 ) ); m.I.add_output( OpE( sqr( m.x0 ), m.t, tf ) );
                                  m.I.add_output( OpE( m.x1, m.t, tf ) ); m.I.add_output( OpI( std::vector<mc::FFVar>{ m.x1 }, { m.t } )[0] ); } );
  R const r0 = solve( ref );
  std::cout << "  reference (4 separate outputs): solves NORMAL=" << r0.ok << "   fwd~adj on the references: ";
  double wr = 0.; for( size_t k=0;k<4;++k ) wr = std::max( wr, std::fabs( r0.Ga[k]-r0.Gf[k] )/std::fabs( r0.Gf[k] ) ); std::cout << wr << "\n";
  { // ORACLE: the stage-independent reference outputs x0^2(tf), x1(tf), int x1 against a SINGLE-stage model
    M one; build( one, [&]( M& m ){ m.I.add_output( OpE( sqr( m.x0 ), m.t, tf ) ); m.I.add_output( OpE( m.x1, m.t, tf ) );
                                   m.I.add_output( OpI( std::vector<mc::FFVar>{ m.x1 }, { m.t } )[0] ); }, 1 );
    R const r1 = solve( one );  double wv = 0., wg = 0.;
    for( size_t k = 0; k < 3; ++k ){ wv = std::max( wv, std::fabs( r0.F[k+1]-r1.F[k] )/std::fabs( r1.F[k] ) ); wg = std::max( wg, std::fabs( r0.Gf[k+1]-r1.Gf[k] )/std::fabs( r1.Gf[k] ) ); }
    std::cout << "  ORACLE 2-stage references == single-stage (value " << wv << ", gradient " << wg << ")" << ( wv < 1e-7 && wg < 1e-6? "  OK": "  BAD" ) << "\n"; }
  { // an evaluation OFF the element grid: ODESLV inserts a stage boundary there -- same value as a model whose elements
    // explicitly have a boundary at 2.5 (4 elements on [0,10])
    M off, on;
    bool const s1 = build( off, [&]( M& m ){ m.I.add_output( OpE( sqr( m.x0 ), m.t, 2.5 ) ); } );
    bool const s2 = build( on,  [&]( M& m ){ m.I.add_output( OpE( sqr( m.x0 ), m.t, 2.5 ) ); }, 4 );
    R const a = solve( off ), b = solve( on );
    double const w = ( a.ok && b.ok )? std::fabs( a.F[0]-b.F[0] )/std::fabs( b.F[0] ): 1.;
    std::cout << "  evaluation at t = 2.5, off the element grid: setup=" << s1 << s2 << "  == model with a boundary at 2.5 (rel " << w << ")" << ( s1 && s2 && w < 1e-7? "  OK": "  BAD" ) << "\n";
    M out; bool const s3 = build( out, [&]( M& m ){ m.I.add_output( OpE( sqr( m.x0 ), m.t, 12. ) ); } );
    std::cout << "  evaluation at t = 12, OUTSIDE the domain: setup=" << s3 << "  err=[" << out.I.extract_error().substr( 0, 60 ) << "]" << ( !s3? "  REFUSED as it should be": "  BAD" ) << "\n"; }
  struct C { char const* name; std::function<mc::FFVar( M& )> g; std::function<double( R const&, int )> want; };   // want(ref, 0=value 1=grad)
  auto V = [&]( R const& r, int w, int k ){ return w? r.Gf[k]: r.F[k]; };
  std::vector<C> cases{
    { "x0^2(T1) + x0^2(tf)              ", [&]( M& m ){ return OpE( sqr( m.x0 ), m.t, T1 ) + OpE( sqr( m.x0 ), m.t, tf ); },     [&]( R const& r, int w ){ return V(r,w,0) + V(r,w,1); } },
    { "x0^2(tf) + x0^2(tf)  (same record)", [&]( M& m ){ return OpE( sqr( m.x0 ), m.t, tf ) + OpE( sqr( m.x0 ), m.t, tf ); },     [&]( R const& r, int w ){ return 2.*V(r,w,1); } },
    { "2 x0^2(T1) - 3 int x1 + 5        ", [&]( M& m ){ return 2.*OpE( sqr( m.x0 ), m.t, T1 ) - 3.*OpI( std::vector<mc::FFVar>{ m.x1 }, { m.t } )[0] + 5.; },
                                                                                                                      [&]( R const& r, int w ){ return 2.*V(r,w,0) - 3.*V(r,w,3) + ( w? 0.: 5. ); } },
    { "p int x1 + x1(tf)  (param coeff) ", [&]( M& m ){ return m.p*OpI( std::vector<mc::FFVar>{ m.x1 }, { m.t } )[0] + OpE( m.x1, m.t, tf ); },
                                                                                                                      [&]( R const& r, int w ){ return w? r.F[3] + P*r.Gf[3] + r.Gf[2]: P*r.F[3] + r.F[2]; } } };
  for( auto& c : cases ){
    M m; bool const su = build( m, [&]( M& mm ){ mm.I.add_output( c.g( mm ) ); } );
    std::cout << "  " << c.name << "  setup=" << su << " nf=" << m.I.nf();
    if( !su ){ std::cout << "  err=[" << m.I.extract_error() << "]\n"; continue; }
    R const r = solve( m );
    double const ev = std::fabs( r.F[0] - c.want( r0, 0 ) ) / std::fabs( c.want( r0, 0 ) );
    double const ef = std::fabs( r.Gf[0] - c.want( r0, 1 ) ) / std::fabs( c.want( r0, 1 ) );
    double const ea = std::fabs( r.Ga[0] - r.Gf[0] ) / std::fabs( r.Gf[0] );
    std::cout << "  solves=" << r.ok << "  value~ref " << ev << "  fwd~ref " << ef << "  adj~fwd " << ea << ( r.ok && ev < 1e-8 && ef < 1e-7 && ea < 1e-6? "  OK": "  BAD" ) << "\n";
  }
  M nl; bool const su = build( nl, [&]( M& m ){ m.I.add_output( OpE( sqr( m.x0 ), m.t, T1 ) * OpE( m.x1, m.t, tf ) ); } );
  std::cout << "  x0^2(T1) * x1(tf)  (NONLINEAR)    setup=" << su << "  err=[" << nl.I.extract_error() << "]" << ( !su? "  REFUSED as it should be": "  BAD" ) << "\n";
}
