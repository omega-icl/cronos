// OCFE_MARCHflux.cpp -- monolithic vs marching OCFESLV on the nonlinear heat equation
//   u_t = (D(u) u_z)_z,  D(u) = D0 (1 + beta u),  z in (0,1), t in (0,1],  u(t,0) = 0,
// written expanded: u_t - D(u) u_zz - D0 beta u_z^2 = 0.  At z = 1, by argument BC:
//   0  Dirichlet       u(t,1) = 1                 (IC u0 = z)
//   1  linear flux     D0 u_z(t,1) = 1            (IC u0 = a z + c (z - z^3/3), a = 1/D0)
//   2  nonlinear flux  D(u) u_z(t,1) = 1          (IC u0 = a z + c (z - z^3/3), D0 (1 + beta (a + 2c/3)) a = 1)
//   3  as 2, the flux being a piecewise-constant INPUT q(t) = 1 (opens the causal evolution-interface gate)
// The initial data satisfy both boundary conditions (no corner singularity) and are transient.  Prints u(1,1/2)
// from both modes and, through eval_solution(), max |u_mono - u_march| over z at the end of every window.
// Usage: OCFE_MARCHflux [beta=1] [bc=2] [NT=10] [display=0] [imposition=-1 (default) | 0 IC_WEAK, 1 IC_STRONG, 2 IC_TRACE] [NZ=6] [SAT_SIGMA0=default 1] [no_output=0] [T=1] [reduce=-1 default | 0 RED_NONE] [setq=0] [noout=0]
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include "ocfeslv.hpp"
using namespace mc;
typedef FFModel::EqnRole Role;

struct Run { std::unique_ptr<FFGraph> G; std::unique_ptr<OCFESLV> S; FFVar t, z, u; bool ok = false; double out = std::nan( "" ); };

static int g_imp = -1;   // INTERFACE.IMPOSITION override (-1: default)
static size_t g_nz = 6;   // z elements
static double g_sigma0 = -1.;   // INTERFACE.SAT_SIGMA0 override (< 0: default)
static double g_T = 1.;   // horizon: NT elements over [0, T]
static int g_red = -1;   // REDUCE.ORDER override (-1 default; 0 RED_NONE)
static int g_setq = 0;   // bc 3: also write q = 1 explicitly with set_input_values (0: rely on init()'s reference fill)
static bool g_noout = false;   // no outputs (no Evtz_u_R auxiliary): read u(1,1/2) through eval_solution
static Run solve( bool marching, double beta, int bc, size_t NT, int disp, double D0 = 0.1, double c = 1. ){
  Run r;  r.G.reset( new FFGraph );  FFGraph& G = *r.G;  FFPartial OpP;  FFEval OpE;
  r.t = G.add_var( "t" );  r.z = G.add_var( "z" );  r.u = G.add_var( "u" );
  FFVar const &t = r.t, &z = r.z, &u = r.u;
  r.S.reset( new OCFESLV( &G ) );  OCFESLV& S = *r.S;
  S.add_domain( t, FFDom( 0., g_T, NT, FFDom::LGR, 3 ) );  S.add_domain( z, FFDom( 0., 1., g_nz, FFDom::LGL, 5 ) );
  S.set_evolution_domain( t );  S.add_state( u, {t,z} );  S.update_ref( u, 0. );
  int const TI = FFDom::ALL - FFDom::LB, ZI = FFDom::ALL - FFDom::LB - FFDom::UB;
  FFVar const D = D0 * ( 1. + beta * u );
  S.add_equation( OpP( u, t ) - D * OpP( u, {{z,2}} ) - D0 * beta * sqr( OpP( u, z ) ), {t,z}, {TI,ZI}, OCFESLV::EqnOptions( Role::INTERIOR ) );
  S.add_equation( u, {t,z}, {TI,FFDom::LB}, OCFESLV::EqnOptions( Role::BOUNDARY ) );
  double a = 1.;
  if( bc == 0 )      S.add_equation( u - 1., {t,z}, {TI,FFDom::UB}, OCFESLV::EqnOptions( Role::BOUNDARY ) );
  else if( bc == 1 ){ a = 1. / D0;  S.add_equation( D0 * OpP( u, z ) - 1., {t,z}, {TI,FFDom::UB}, OCFESLV::EqnOptions( Role::BOUNDARY ) ); }
  else if( bc == 3 ){ FFVar const q = G.add_var( "q" );
                      S.add_input( q, {t}, FFDom::LGR, 1 );  S.update_ref( q, 1. );
                      double const b = 1. + 2. * beta * c / 3.;
                      a = beta? ( -b + std::sqrt( b*b + 4. * beta / D0 ) ) / ( 2. * beta ): 1. / D0;
                      S.add_equation( D * OpP( u, z ) - q, {t,z}, {TI,FFDom::UB}, OCFESLV::EqnOptions( Role::BOUNDARY ) ); }
  else{               double const b = 1. + 2. * beta * c / 3.;
                      a = beta? ( -b + std::sqrt( b*b + 4. * beta / D0 ) ) / ( 2. * beta ): 1. / D0;
                      S.add_equation( D * OpP( u, z ) - 1., {t,z}, {TI,FFDom::UB}, OCFESLV::EqnOptions( Role::BOUNDARY ) ); }
  FFVar const u0 = bc == 0? z: a * z + c * ( z - z * z * z / 3. );
  S.add_equation( u - u0, {t,z}, {FFDom::LB,FFDom::ALL}, OCFESLV::EqnOptions( Role::INITIAL ) );
  if( !g_noout ) S.add_output( OpE( u, {{t,1},{z,1}}, {{t,g_T},{z,.5}} ) );
  if( g_imp >= 0 ) S.options.INTERFACE.IMPOSITION = (OCFESLV::Options::ImpositionType)g_imp;
  if( g_sigma0 > 0. ) S.options.INTERFACE.SAT_SIGMA0 = g_sigma0;
  if( g_red == 0 ) S.options.REDUCE.ORDER = FFModel::Options::RED_NONE;
  S.options.SOLVE.MARCHING = marching;  S.options.SOLVE.RES_TOL = 1e-11;  S.options.DISPLAY_LEVEL = disp;
  if( !S.setup() ) return r;
  std::vector<double> var, inp;  S.init( var, inp );
  if( bc == 3 && g_setq ){
    FFVar const* q = nullptr;  for( auto const& [w, d] : S.var_declared_input() ) if( w.name() == "q" ) q = &w;
    if( q ) S.set_input_values( *q, std::vector<double>( S.control_dofs( *q ).size()? S.control_dofs( *q ).size(): NT, 1. ), inp.data() );
    if( disp ){ auto const v = S.get_input_values( *q, inp.data() ); std::printf( "  q after set: %zu values, first %g last %g\n", v.size(), v.front(), v.back() ); }
  }
  if( bc == 3 && disp ){ FFVar const* q = nullptr;  for( auto const& [w, d] : S.var_declared_input() ) if( w.name() == "q" ) q = &w;
    auto const v = S.get_input_values( *q, inp.data() );  std::printf( "  [%s] q in inp: %zu values:", marching? "march": "mono ", v.size() );
    for( double x : v ) std::printf( " %g", x );  std::printf( "\n" ); }
  r.ok = S.solve( var.data(), inp.data() ).converged;
  if( r.ok ) r.out = g_noout? S.eval_solution( u, { { t, g_T }, { z, .5 } } ): S.val_functions()[0];
  return r;
}

int main( int argc, char** argv ){
  double const beta = argc > 1? std::atof( argv[1] ): 1.;
  int    const bc   = argc > 2? std::atoi( argv[2] ): 2;
  size_t const NT   = argc > 3? std::atoi( argv[3] ): 10;
  int    const disp = argc > 4? std::atoi( argv[4] ): 0;
  g_imp = argc > 5? std::atoi( argv[5] ): -1;
  g_nz  = argc > 6? std::atoi( argv[6] ): 6;
  g_sigma0 = argc > 7? std::atof( argv[7] ): -1.;
  g_noout = argc > 8 && std::atoi( argv[8] );
  g_T     = argc > 9? std::atof( argv[9] ): 1.;
  g_red   = argc > 10? std::atoi( argv[10] ): -1;
  g_setq  = argc > 11? std::atoi( argv[11] ): 0;
  char const* BC[] = { "Dirichlet u=1", "linear flux D0 u_z=1", "nonlinear flux D(u) u_z=1", "nonlinear flux D(u) u_z=q(t), q=1 input" };
  Run M = solve( false, beta, bc, NT, disp ), K = solve( true, beta, bc, NT, disp );
  std::printf( "beta=%g  %s  NT=%zu NZ=%zu | u(1,1/2): monolithic %.10f (%s)  marching %.10f (%s) | diff %.2e\n", beta, BC[bc], NT, g_nz,
               M.out, M.ok? "ok": "FAILED", K.out, K.ok? "ok": "FAILED", std::fabs( M.out - K.out ) );
  if( !M.ok || !K.ok ) return 1;
  { typename OCFESLV::t_Coord pm{ { M.t, .1 }, { M.z, .5 } }, pk{ { K.t, .1 }, { K.z, .5 } };
    std::printf( "  u(0.1, 1/2): monolithic %.12f  marching %.12f\n", M.S->eval_solution( M.u, pm ), K.S->eval_solution( K.u, pk ) ); }
  std::printf( "  window end t | max_z |u_mono - u_march|   (z at the 31 points of a uniform grid)\n" );
  for( size_t k = 1; k <= NT; ++k ){
    double const tk = g_T * double( k ) / NT;  double dmax = 0.;
    for( int j = 0; j <= 30; ++j ){
      typename OCFESLV::t_Coord pt{ { M.t, tk }, { M.z, j / 30. } }, pk{ { K.t, tk }, { K.z, j / 30. } };
      dmax = std::max( dmax, std::fabs( M.S->eval_solution( M.u, pt ) - K.S->eval_solution( K.u, pk ) ) );
    }
    if( k <= 3 || k == NT || NT <= 10 ) std::printf( "    %8.4f   %.3e\n", tk, dmax );
  }
  return 0;
}
