// ============================================================================
//  OCFE_Burger_fixed.cpp -- Burgers' equation on a FIXED mesh: every solve mode against every interface imposition
//  (2026-10-11).  Companion of OCFE_Burger_move.cpp (the same physics on a moving mesh).
//
//    u_t = eps u_xx - u u_x,  0 < x < 1,  u(0,t) = u(1,t) = 0,  u(x,0) = sin(2 pi x) + 0.5 sin(pi x)
//    (Problem 2 of Huang, Ren & Russell, J. Comput. Phys. 113 (1994), there with eps = 1e-4 on a moving mesh).
//
//  First-order form on the fixed mesh x = xi:  q = u_x,  r = u_xx (= q_x),
//    u_t = eps r - u q                  (interior t, interior x)
//    q - u_x = 0,  r - q_x = 0          (all t, all x)
//    u = 0 at x = 0, 1 (interior t);    u = u0(x) at t = 0
//
//  For each configuration in {monolithic, marching} x {weak, trace, strong}: setup status, convergence, time, the
//  largest seam jump of u and q (layout-free: each state evaluated just left and right of every element boundary),
//  and u on a common grid at t = tf, compared pairwise against the first configuration that converged.
//  The initial guess comes from reference functions (update_ref before setup; init() seeds from them), so it is
//  independent of the state-vector layout, which differs between solve modes.
//
//  A resolved front needs eps of order 5e-2 on the default mesh (front width ~ 2 eps / [u]); at eps = 1e-2 and below
//  the fixed mesh under-resolves it, which is the benchmark's point and the reason for OCFE_Burger_move.
//
//  USAGE  ./OCFE_Burger_fixed [--eps E] [--nel N] [--nnd n] [--nelt M] [--tf T] [--rtol R] [--maxit K]
//                             [--cases list] [--form q|r] [--link] [--out case [--outfile prefix]] [--disp d]
//    --form q (default): u_t = eps q_x - u q with q = u_x -- the minimal first-order form;  --form r: also r = q_x
//    --link  : declare the definition rows LINK (the role reduce_order() gives its auxiliaries)
//    --out   : write that case's u(t, x) at 11 times to prefix.dat with a gnuplot script prefix.gp (-> prefix.png)
//    --cases : comma-separated subset of mono-weak,mono-trace,mono-strong,march-weak,march-trace,march-strong
//  Exit status 0 iff every configuration run converged and all converged profiles agree to --agree (1e-6).
// ============================================================================
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <chrono>
#include <string>
#include <vector>
#include <iostream>
#include "ocfeslv.hpp"
using namespace mc;

static double const kPI = 3.14159265358979323846;
typedef decltype( OCFESLV::Options::IC_WEAK ) Imp;
struct Par { double eps = 5e-2, tf = 0.5, rtol = 1e-7, agree = 1e-6; size_t nel = 10, nnd = 8, nelt = 5; int maxit = 60, disp = 0;
             std::string form = "q";        // "q": q = u_x only, u_t = eps q_x - u q;  "r": also r = q_x, u_t = eps r - u q
             bool link = false;             // declare the definition rows LINK (as reduce_order() does) instead of INTERIOR
             std::string out, outfile = "burger_fixed"; };
struct Res { std::string name; bool setup = false, conv = false; int its = 0; double resid = NAN, sec = 0., ju = NAN, jq = NAN; std::vector<double> uf; };

static double u0 ( double z ){ return std::sin( 2*kPI*z ) + 0.5*std::sin( kPI*z ); }
static double u0x( double z ){ return 2*kPI*std::cos( 2*kPI*z ) + 0.5*kPI*std::cos( kPI*z ); }
static double u0xx( double z ){ return -4*kPI*kPI*std::sin( 2*kPI*z ) - 0.5*kPI*kPI*std::sin( kPI*z ); }

static Res run( Par const& P, bool mono, Imp imp, std::string const& name )
{
  Res R; R.name = name;
  auto const t0 = std::chrono::steady_clock::now();
  FFGraph DAG; FFPartial OpP;
  FFVar t = DAG.add_var( "t" ), x = DAG.add_var( "x" );
  FFVar u = DAG.add_var( "u(t,x)" ), q = DAG.add_var( "q(t,x)" ), r = DAG.add_var( "r(t,x)" );
  OCFESLV oc( &DAG );
  oc.add_domain( t, FFDom( 0., P.tf, P.nelt, FFDom::LGR, 4 ) );
  oc.add_domain( x, FFDom( 0., 1., P.nel, FFDom::LGL, P.nnd ) );
  oc.set_evolution_domain( t );
  bool const RF = ( P.form == "r" );
  oc.add_state( u, { t, x } ); oc.add_state( q, { t, x } ); if( RF ) oc.add_state( r, { t, x } );
  int const ALL = FFDom::ALL, LB = FFDom::LB, UB = FFDom::UB, T_INT = ALL - LB, X_INT = ALL - LB - UB;
  OCFESLV::EqnOptions inter( OCFESLV::EqnRole::INTERIOR ), bnd( OCFESLV::EqnRole::BOUNDARY ), ini( OCFESLV::EqnRole::INITIAL );
  OCFESLV::EqnOptions def( P.link ? OCFESLV::EqnRole::LINK : OCFESLV::EqnRole::INTERIOR );
  oc.add_equation( OpP( u, t ) - P.eps * ( RF ? r : OpP( q, x ) ) + u * q, { t, x }, { T_INT, X_INT }, inter );
  oc.add_equation( q - OpP( u, x ),                 { t, x }, { ALL, ALL },     def );
  if( RF ) oc.add_equation( r - OpP( q, x ),        { t, x }, { ALL, ALL },     def );
  oc.add_equation( u - ( sin( 2.*kPI*x ) + 0.5*sin( kPI*x ) ), { t, x }, { LB, ALL }, ini );
  oc.add_equation( u, { t, x }, { T_INT, LB }, bnd );
  oc.add_equation( u, { t, x }, { T_INT, UB }, bnd );
  // initial guess: the initial data at every time (layout-free, through the reference functions)
  oc.update_ref( u, [x]( OCFESLV::t_Coord const& c ){ return u0  ( c.at( x ) ); } );
  oc.update_ref( q, [x]( OCFESLV::t_Coord const& c ){ return u0x ( c.at( x ) ); } );
  if( RF ) oc.update_ref( r, [x]( OCFESLV::t_Coord const& c ){ return u0xx( c.at( x ) ); } );
  oc.options.INTERFACE.IMPOSITION = imp;
  oc.options.SOLVE.MARCHING = !mono;
  oc.options.SOLVE.RES_TOL = P.rtol; oc.options.SOLVE.MAX_ITER = P.maxit;
  oc.options.DISPLAY_LEVEL = P.disp;
  try { R.setup = oc.setup(); } catch( std::exception const& e ){ std::printf( "  %s: setup threw %s\n", name.c_str(), e.what() ); }
  if( R.setup ){
    std::vector<double> var, inp; oc.init( var, inp );
    OCFESLV::SolveReport rep = oc.solve( var.data(), nullptr, nullptr );
    R.conv = rep.converged; R.its = rep.iterations; R.resid = rep.final_residual;
    auto ev = [&]( FFVar const& s, double tt, double z ){ OCFESLV::t_Coord pt; pt[t] = tt; pt[x] = z;
      try { return oc.eval_colloc<double>( s, pt, var.data(), nullptr, nullptr ); } catch( ... ){ return std::nan( "" ); } };
    R.ju = R.jq = 0.;
    for( int k = 1; k <= 10; ++k ) for( size_t e = 1; e < P.nel; ++e ){
      double const tt = P.tf * k / 10., z = double( e ) / double( P.nel );
      R.ju = std::max( R.ju, std::fabs( ev( u, tt, z - 1e-10 ) - ev( u, tt, z + 1e-10 ) ) );
      R.jq = std::max( R.jq, std::fabs( ev( q, tt, z - 1e-10 ) - ev( q, tt, z + 1e-10 ) ) );
    }
    for( int s = 0; s <= 200; ++s ) R.uf.push_back( ev( u, P.tf, s / 200. ) );
    if( name == P.out ){   // the solution u(t, x) for gnuplot: 11 times, blank-line separated blocks, and a script that plots them
      std::string const dat = P.outfile + ".dat", gp = P.outfile + ".gp", png = P.outfile + ".png";
      std::FILE* f = std::fopen( dat.c_str(), "w" );
      std::fprintf( f, "# %s: Burgers, eps = %g, fixed mesh %zu x %zu LGL, %zu LGR elements on [0, %g]\n# t x u\n", name.c_str(), P.eps, P.nel, P.nnd, P.nelt, P.tf );
      for( int k = 0; k <= 10; ++k ){ double const tt = P.tf * k / 10.;
        for( int s = 0; s <= 400; ++s ) std::fprintf( f, "%.6g %.6g %.10g\n", tt, s / 400., ev( u, tt, s / 400. ) );
        std::fprintf( f, "\n\n" ); }
      std::fclose( f );
      std::FILE* g = std::fopen( gp.c_str(), "w" );
      std::fprintf( g, "set terminal pngcairo size 900,560 enhanced font ',11'\nset output '%s'\n", png.c_str() );
      std::fprintf( g, "set title 'Burgers, {/Symbol e} = %g, fixed mesh, %s'\nset xlabel 'x'\nset ylabel 'u'\nset grid\nset key outside right\nunset colorbox\nset cbrange [0:1]\n", P.eps, name.c_str() );
      std::fprintf( g, "plot for [k=0:10] '%s' index k using 2:3 with lines lw 1.5 lc palette frac k/10. title sprintf('t = %%.3g', k*%g/10.)\n", dat.c_str(), P.tf );
      std::fclose( g );
      std::printf( "  wrote %s and %s (run: gnuplot %s -> %s)\n", dat.c_str(), gp.c_str(), gp.c_str(), png.c_str() );
    }
  }
  R.sec = std::chrono::duration<double>( std::chrono::steady_clock::now() - t0 ).count();
  return R;
}

int main( int argc, char** argv )
{
  Par P; std::string cases = "mono-weak,mono-trace,mono-strong,march-weak,march-trace,march-strong";
  for( int i = 1; i < argc; ++i ){
    auto num = [&]( char const* f ){ return !std::strcmp( argv[i], f ) && i+1 < argc; };
    if( num( "--eps"   ) ){ P.eps   = std::atof( argv[++i] ); continue; }
    if( num( "--nel"   ) ){ P.nel   = std::atoi( argv[++i] ); continue; }
    if( num( "--nnd"   ) ){ P.nnd   = std::atoi( argv[++i] ); continue; }
    if( num( "--nelt"  ) ){ P.nelt  = std::atoi( argv[++i] ); continue; }
    if( num( "--tf"    ) ){ P.tf    = std::atof( argv[++i] ); continue; }
    if( num( "--rtol"  ) ){ P.rtol  = std::atof( argv[++i] ); continue; }
    if( num( "--maxit" ) ){ P.maxit = std::atoi( argv[++i] ); continue; }
    if( num( "--agree" ) ){ P.agree = std::atof( argv[++i] ); continue; }
    if( num( "--disp"  ) ){ P.disp  = std::atoi( argv[++i] ); continue; }
    if( num( "--cases" ) ){ cases   = argv[++i]; continue; }
    if( num( "--form"  ) ){ P.form  = argv[++i]; continue; }
    if( num( "--out"   ) ){ P.out   = argv[++i]; continue; }
    if( num( "--outfile" ) ){ P.outfile = argv[++i]; continue; }
    if( !std::strcmp( argv[i], "--link" ) ){ P.link = true; continue; }
    std::cerr << "unknown argument " << argv[i] << "\n"; return 2;
  }
  std::printf( "OCFE_Burger_fixed: eps=%g, x: %zu LGL elements of %zu nodes, t: %zu LGR elements on [0,%g], rtol=%g, form %s%s\n",
               P.eps, P.nel, P.nnd, P.nelt, P.tf, P.rtol, P.form.c_str(), P.link ? ", LINK rows" : "" );
  std::vector<Res> all;
  for( std::string const mode : { "mono", "march" } )
    for( auto const& [iname, imp] : std::vector<std::pair<std::string,Imp>>{
           { "weak", OCFESLV::Options::IC_WEAK }, { "trace", OCFESLV::Options::IC_TRACE }, { "strong", OCFESLV::Options::IC_STRONG } } ){
      std::string const name = mode + "-" + iname;
      if( ( "," + cases + "," ).find( "," + name + "," ) == std::string::npos ) continue;
      all.push_back( run( P, mode == "mono", imp, name ) );
      Res const& R = all.back();
      std::printf( "  %-13s setup %-3s conv %-3s its %3d  |r| %.2e  %7.1f s  seam jump u %.1e  q %.1e\n", R.name.c_str(),
                   R.setup ? "ok" : "NO", R.conv ? "yes" : "no", R.its, R.resid, R.sec, R.ju, R.jq );
      std::fflush( stdout );
    }
  // agreement of the converged final profiles.  TESTED: the same imposition across solve modes (mono vs march), and
  // trace against strong (both exact).  REPORTED only: weak against strong -- the penalty method's own consistency error.
  auto find = [&]( std::string const& n )->Res const* { for( auto const& R : all ) if( R.name == n && R.conv ) return &R; return nullptr; };
  auto diff = []( Res const& a, Res const& b ){ double d = 0.; for( size_t s_ = 0; s_ < a.uf.size() && s_ < b.uf.size(); ++s_ ) d = std::max( d, std::fabs( a.uf[s_] - b.uf[s_] ) ); return d; };
  bool ok = !all.empty();
  for( auto const& R : all ) ok = ok && R.conv;
  std::printf( "  final profiles u(x, tf), max |difference| on 201 points:\n" );
  for( std::string const imp : { "weak", "trace", "strong" } ){
    Res const *a = find( "mono-" + imp ), *b = find( "march-" + imp );
    if( a && b ){ double const d = diff( *a, *b ); std::printf( "    mono vs march, %-7s %.2e%s\n", imp.c_str(), d, d > P.agree ? "  <-- DIFFERS" : "" ); ok = ok && d <= P.agree; }
  }
  for( std::string const mode : { "mono", "march" } ){
    Res const *st = find( mode + "-strong" ), *tr = find( mode + "-trace" ), *wk = find( mode + "-weak" );
    if( st && tr ){ double const d = diff( *tr, *st ); std::printf( "    %-5s trace vs strong   %.2e%s\n", mode.c_str(), d, d > P.agree ? "  <-- DIFFERS" : "" ); ok = ok && d <= P.agree; }
    if( st && wk ) std::printf( "    %-5s weak vs strong    %.2e  (reported: the penalty's consistency error)\n", mode.c_str(), diff( *wk, *st ) );
  }
  std::printf( "%s\n", ok ? "PASS: all configurations converged; same-imposition and exact-imposition profiles agree" : "FAIL" );
  return ok ? 0 : 1;
}
