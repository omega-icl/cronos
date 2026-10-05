// OCFE_CAUSAL.cpp -- the causal time discretisation of a monolithic OCFESLV solve (always on since 2026-10-03), and init()'s
// reference fill of a time-distributed input under marching, and a step in a boundary input (defect 3) (2026-10-01).  Model: u_t = (D(u) u_z)_z,
// D = 0.1 (1 + u), z in (0,1), t in (0,T], u(t,0) = 0, D(u) u_z(t,1) = q (a constant, or a piecewise-constant input
// q(t) whose values come from init()), initial data compatible and transient; strong imposition.
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <memory>
#include <vector>
#include "ocfeslv.hpp"
using namespace mc;
typedef FFModel::EqnRole Role;
static int npass = 0, nfail = 0;
static void check( bool c, char const* w ){ std::printf( "  %s  %s\n", c? "PASS": "FAIL", w ); c? ++npass: ++nfail; }

struct Run { std::unique_ptr<FFGraph> G; std::unique_ptr<OCFESLV> S; FFVar t, z, u; bool ok = false; };
static Run solve( bool marching, size_t NT, double T, bool causal, bool q_input, bool step = false, bool trace = false ){
  Run r;  r.G.reset( new FFGraph );  FFGraph& G = *r.G;  FFPartial OpP;
  r.t = G.add_var( "t" );  r.z = G.add_var( "z" );  r.u = G.add_var( "u" );
  FFVar const &t = r.t, &z = r.z, &u = r.u;
  r.S.reset( new OCFESLV( &G ) );  OCFESLV& S = *r.S;
  S.add_domain( t, FFDom( 0., T, NT, FFDom::LGR, 3 ) );  S.add_domain( z, FFDom( 0., 1., 6, FFDom::LGL, 5 ) );
  S.set_evolution_domain( t );  S.add_state( u, {t,z} );  S.update_ref( u, 0. );
  int const TI = FFDom::ALL - FFDom::LB, ZI = FFDom::ALL - FFDom::LB - FFDom::UB;
  double const D0 = .1, beta = 1., c = 1., b = 1. + 2. * beta * c / 3.;
  double const a = ( -b + std::sqrt( b*b + 4. * beta / D0 ) ) / ( 2. * beta );
  FFVar const D = D0 * ( 1. + beta * u );
  S.add_equation( OpP( u, t ) - D * OpP( u, {{z,2}} ) - D0 * beta * sqr( OpP( u, z ) ), {t,z}, {TI,ZI}, OCFESLV::EqnOptions( Role::INTERIOR ) );
  S.add_equation( u, {t,z}, {TI,FFDom::LB}, OCFESLV::EqnOptions( Role::BOUNDARY ) );
  if( q_input ){ FFVar const q = G.add_var( "q" );  S.add_input( q, {t}, FFDom::LGR, 1 );  S.update_ref( q, 1. );
                 S.add_equation( D * OpP( u, z ) - q, {t,z}, {TI,FFDom::UB}, OCFESLV::EqnOptions( Role::BOUNDARY ) ); }
  else           S.add_equation( D * OpP( u, z ) - 1., {t,z}, {TI,FFDom::UB}, OCFESLV::EqnOptions( Role::BOUNDARY ) );
  S.add_equation( u - a * z - c * ( z - z * z * z / 3. ), {t,z}, {FFDom::LB,FFDom::ALL}, OCFESLV::EqnOptions( Role::INITIAL ) );
  S.options.INTERFACE.IMPOSITION = trace? OCFESLV::Options::IC_TRACE: OCFESLV::Options::IC_STRONG;  (void)causal;   /* the time discretisation is always causal (2026-10-03) */
  S.options.SOLVE.MARCHING = marching;  S.options.SOLVE.RES_TOL = 1e-11;  S.options.DISPLAY_LEVEL = 0;
  if( !S.setup() ) return r;
  std::vector<double> var, inp;  S.init( var, inp );                     // q's values: init()'s reference fill only
  if( q_input && step ){                                                  // a STEP heater: 1 | 0.25 at T/2
    FFVar const* qv = nullptr;  for( auto const& [w, d] : S.var_declared_input() ) if( w.name() == "q" ) qv = &w;
    std::vector<double> lev( NT );  for( size_t k = 0; k < NT; ++k ) lev[k] = ( k < NT / 2 )? 1.: .25;
    S.set_input_values( *qv, lev, inp.data() );
  }
  r.ok = S.solve( var.data(), inp.data() ).converged;
  return r;
}

static double at( Run const& r, double t, double z ){
  typename OCFESLV::t_Coord pt{ { r.t, t }, { r.z, z } };  return r.S->eval_solution( r.u, pt ); }

int main(){
  double const CAUSAL = 2.068511328630;   // u(1,1/2) at NT = 10, monolithic (the former two-way value: 2.0686363287)
  // 1. the monolithic solve gives the causal value (bit for bit)
  { Run M = solve( false, 10, 1., false, false );
    double const v = at( M, 1., .5 );
    std::printf( "  monolithic u(1,1/2) = %.12f (causal %.12f)\n", v, CAUSAL );
    check( M.ok && std::fabs( v - CAUSAL ) < 5e-11, "monolithic: the causal value (2.06851132863 at NT = 10)" ); }
  // 2. causal: u(0.1, 1/2) independent of the elements that follow
  { Run A = solve( false, 1, .1, true, false ), B = solve( false, 2, .2, true, false ), C = solve( false, 10, 1., true, false );
    double const a = at( A, .1, .5 ), b = at( B, .1, .5 ), c = at( C, .1, .5 );
    std::printf( "  monolithic u(0.1,1/2): 1 el %.12f  2 el %.12f  10 el %.12f\n", a, b, c );
    check( A.ok && B.ok && C.ok && std::fabs( a - b ) < 1e-11 && std::fabs( a - c ) < 1e-11,
           "the monolithic solve is causal (u(0.1) independent of later elements, 1e-11)" ); }
  // 3. monolithic == marching
  { Run M = solve( false, 10, 1., true, false ), K = solve( true, 10, 1., true, false );
    double const m = at( M, 1., .5 ), k = at( K, 1., .5 );
    std::printf( "  monolithic %.12f  marching %.12f\n", m, k );
    check( M.ok && K.ok && std::fabs( m - k ) < 1e-10, "monolithic == marching (1e-10)" ); }
  // 4. init(): a time-distributed input with a constant reference fills every window under marching
  { Run M = solve( false, 10, 1., false, true ), K = solve( true, 10, 1., false, true );
    double const m = at( M, 1., .5 ), k = at( K, 1., .5 );
    std::printf( "  q(t) from init(): monolithic %.12f  marching %.12f\n", m, k );
    check( M.ok && K.ok && std::fabs( m - k ) < 1e-10, "init(): q(t) = 1 in every window under marching (marching == monolithic, 1e-10)" ); }
  // 5. a STEP in a boundary input: the boundary node stays continuous in time (monolithic == marching everywhere,
  //    including the heated end at the switch -- defect 3: the boundary row used to hold the element's first node)
  { Run M = solve( false, 10, 1., false, true, true ), K = solve( true, 10, 1., false, true, true );
    double dmax = 0., jump = std::fabs( at( M, .5, 1. ) - at( M, .5 - 1e-7, 1. ) );
    for( int i = 0; i <= 40; ++i ) for( int j = 0; j <= 10; ++j )
      dmax = std::max( dmax, std::fabs( at( M, i / 40., j / 10. ) - at( K, i / 40., j / 10. ) ) );
    std::printf( "  step heater: max |mono - march| over (t,z) = %.2e | monolithic jump at (0.5, 1) = %.2e\n", dmax, jump );
    check( M.ok && K.ok && dmax < 1e-10 && jump < 1e-5,
           "step boundary input: no jump at the boundary node, monolithic == marching everywhere (1e-10)" ); }
  // 6. the same under IC_TRACE (the per-node rule holds for every imposition)
  { Run M = solve( false, 10, 1., false, true, true, true ), K = solve( true, 10, 1., false, true, true, true );
    double const jump = std::fabs( at( M, .5, 1. ) - at( M, .5 - 1e-7, 1. ) ), dm = std::fabs( at( M, 1., .5 ) - at( K, 1., .5 ) );
    std::printf( "  IC_TRACE step heater: monolithic jump at (0.5, 1) = %.2e | |mono - march| at (1, 1/2) = %.2e\n", jump, dm );
    check( M.ok && K.ok && jump < 1e-5 && dm < 1e-10, "IC_TRACE: step boundary input, no jump, monolithic == marching" ); }
  // 7. the reduced-order flux of a 2D heat problem at the second element's first time node, on and off the
  //    x interface: monolithic == marching (the auxiliaries of OCFE_PDE0 / PDE5 drifted at spatial interfaces of
  //    causal element starts: the time face was still treated as a corner of the spatial face)
  { auto heat2d = [&]( bool march, std::vector<double>& out ) -> bool {
      FFGraph G;  FFPartial OpP;
      FFVar t = G.add_var( "t" ), x = G.add_var( "x" ), y = G.add_var( "y" ), U = G.add_var( "U" );
      OCFESLV S( &G );
      S.add_domain( t, FFDom( 0., .1, 2, FFDom::LGR, 4 ) );  S.add_domain( x, FFDom( 0., 1., 2, FFDom::LGL, 5 ) );
      S.add_domain( y, FFDom( 0., 1., 2, FFDom::LGL, 5 ) );
      S.set_evolution_domain( t );  S.add_state( U, {t,x,y} );  S.update_ref( U, 0. );
      int const TI = FFDom::ALL - FFDom::LB, I = FFDom::ALL - FFDom::LB - FFDom::UB;
      S.add_equation( OpP( U, t ) - OpP( U, {{x,2}} ) - OpP( U, {{y,2}} ), {t,x,y}, {TI,I,I}, OCFESLV::EqnOptions( Role::INTERIOR ) );
      for( int b : { FFDom::LB, FFDom::UB } ){
        S.add_equation( U, {t,x,y}, {TI,b,FFDom::ALL}, OCFESLV::EqnOptions( Role::BOUNDARY ) );
        S.add_equation( U, {t,x,y}, {TI,I,b}, OCFESLV::EqnOptions( Role::BOUNDARY ) ); }
      S.add_equation( U - sin( M_PI * x ) * sin( M_PI * y ) - .3 * x * x * y * ( 1. - y ) * ( 1. - x ), {t,x,y},
                      {FFDom::LB,FFDom::ALL,FFDom::ALL}, OCFESLV::EqnOptions( Role::INITIAL ) );
      S.options.INTERFACE.IMPOSITION = OCFESLV::Options::IC_STRONG;  
      S.options.SOLVE.MARCHING = march;  S.options.SOLVE.RES_TOL = 1e-12;  S.options.DISPLAY_LEVEL = 0;
      if( !S.setup() ) return false;
      std::vector<double> var, inp;  S.init( var, inp );
      if( !S.solve( var.data(), inp.data() ).converged ) return false;
      FFVar Dx;  for( auto const& v : S.states_colloc() ) if( v.name().find( "Dx" ) != std::string::npos ) Dx = v;
      out.clear();
      for( double xv : { .25, .5, .75 } ) for( double yv : { .25, .5, .75 } ){
        typename OCFESLV::t_Coord pt{ { t, .05 }, { x, xv }, { y, yv } };
        out.push_back( S.eval_solution( Dx, pt ) ); }
      return true; };
    std::vector<double> M, K;  bool const ok = heat2d( false, M ) && heat2d( true, K );
    double d = 0.;  if( ok ) for( size_t k = 0; k < M.size(); ++k ) d = std::max( d, std::fabs( M[k] - K[k] ) );
    std::printf( "  2D heat, reduced flux Dx_U at the element start (t = 0.05): max |mono - march| = %.2e\n", d );
    check( ok && d < 1e-9, "option on: reduced-order fluxes at an element start, monolithic == marching (1e-9)" ); }
  std::printf( "\n  OCFE_CAUSAL: %d passed, %d failed -- %s\n", npass, nfail, nfail? "FAILURES": "ALL PASS" );
  return nfail? 1: 0;
}
