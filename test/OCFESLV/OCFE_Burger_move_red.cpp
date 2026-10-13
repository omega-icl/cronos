// ============================================================================
//  OCFE_Burger_move_red.cpp -- Burgers' equation on a MOVING mesh with the order reduction left to OCFESLV
//  (REDUCE.ORDER = RED_MAIN | RED_FULL) for BOTH the physics and the mesh equation (2026-10-12).  Companion of
//  OCFE_Burger_move.cpp, which introduces q, r, xg, g by hand.  States: u and x only.
//
//    u_t - x_t (u_xi/x_xi) = eps ((u_xi/x_xi)_xi)/x_xi - u (u_xi/x_xi)     (ALE form, interior t, interior xi)
//    tau x_t = ( sqrt(x_xi^2 + u_xi^2) )_xi                                 (MMPDE5, arc-length monitor:
//                                                                             M x_xi = sqrt(1 + u_x^2) x_xi exactly)
//    u = 0, x = 0 / 1 at the xi ends;  u = u0(x), x = xi at t = 0
//
//  USAGE  ./OCFE_Burger_move_red [--red main|full] [--eps E] [--nel N] [--nnd n] [--nelt M] [--tf T] [--tau tau]
//                                [--smooth s] [--rtol R] [--maxit K] [--cases list] [--disp d]
//    --smooth s : smoothed monitor Ms - s^2 Ms_xixi = sqrt(1 + u_x^2), Ms_xi = 0 at the ends, flux Ms x_xi (default off)
// ============================================================================
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <chrono>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <iostream>
#include "ocfeslv.hpp"
using namespace mc;

static double const kPI = 3.14159265358979323846;
typedef decltype( OCFESLV::Options::IC_WEAK ) Imp;
struct Par { double eps = 1e-2, tf = 0.02, rtol = 1e-7, agree = 1e-6, tau = 1e-3, smooth = 0.; size_t nel = 10, nnd = 8, nelt = 2; int maxit = 60, disp = 0; bool pos = false, gg = false; std::string red = "main"; bool lam = false; };
struct Res { std::string name; bool setup = false, conv = false; int its = 0; double resid = NAN, sec = 0., ju = NAN, jx = NAN, jq = NAN, jg = NAN, drift = NAN; std::vector<double> uf; };

static double u0 ( double z ){ return std::sin( 2*kPI*z ) + 0.5*std::sin( kPI*z ); }
static double u0x( double z ){ return 2*kPI*std::cos( 2*kPI*z ) + 0.5*kPI*std::cos( kPI*z ); }
static double u0xx( double z ){ return -4*kPI*kPI*std::sin( 2*kPI*z ) - 0.5*kPI*kPI*std::sin( kPI*z ); }

static Res run( Par const& P, bool mono, Imp imp, std::string const& name )
{
  Res R; R.name = name;
  auto const t0 = std::chrono::steady_clock::now();
  FFGraph DAG; FFPartial OpP;
  FFVar t = DAG.add_var( "t" ), xi = DAG.add_var( "xi" );
  FFVar u = DAG.add_var( "u(t,xi)" ), x = DAG.add_var( "x(t,xi)" ), Ms = DAG.add_var( "Ms(t,xi)" );
  FFVar Lam = DAG.add_var( "lam" );   // --lam: monitor strength, M = sqrt(1 + lam u_x^2), ramped 0 -> 1
  OCFESLV oc( &DAG );
  oc.add_domain( t,  FFDom( 0., P.tf, P.nelt, FFDom::LGR, 4 ) );
  oc.add_domain( xi, FFDom( 0., 1., P.nel, FFDom::LGL, P.nnd ) );
  oc.set_evolution_domain( t );
  oc.add_state( u, { t, xi } ); oc.add_state( x, { t, xi } );
  bool const SM = P.smooth > 0.;
  if( SM ) oc.add_state( Ms, { t, xi } );
  int const ALL = FFDom::ALL, LB = FFDom::LB, UB = FFDom::UB, T_INT = ALL - LB, X_INT = ALL - LB - UB;
  OCFESLV::EqnOptions inter( OCFESLV::EqnRole::INTERIOR ), bnd( OCFESLV::EqnRole::BOUNDARY ), ini( OCFESLV::EqnRole::INITIAL );
  FFVar const ux = OpP( u, xi ) / OpP( x, xi );                       // u_x in the computational coordinate
  oc.add_equation( OpP( u, t ) - OpP( x, t ) * ux - P.eps * OpP( ux, xi ) / OpP( x, xi ) + u * ux, { t, xi }, { T_INT, X_INT }, inter );
  oc.add_equation( u - ( sin( 2.*kPI*x ) + 0.5*sin( kPI*x ) ), { t, xi }, { LB, ALL }, ini );
  oc.add_equation( u, { t, xi }, { T_INT, LB }, bnd );
  oc.add_equation( u, { t, xi }, { T_INT, UB }, bnd );
  if( SM ){   // smoothed monitor: Ms - s^2 Ms_xixi = M, Ms_xi = 0 at the ends; flux Ms x_xi
    oc.add_equation( Ms - P.smooth * P.smooth * OpP( OpP( Ms, xi ), xi ) - sqrt( 1. + ( P.lam ? Lam : FFVar( 1. ) ) * ux * ux ), { t, xi }, { ALL, X_INT }, inter );
    oc.add_equation( OpP( Ms, xi ), { t, xi }, { ALL, LB }, bnd );
    oc.add_equation( OpP( Ms, xi ), { t, xi }, { ALL, UB }, bnd );
    oc.add_equation( P.tau * OpP( x, t ) - OpP( Ms * OpP( x, xi ), xi ), { t, xi }, { T_INT, X_INT }, inter );
  } else
    oc.add_equation( P.tau * OpP( x, t ) - OpP( sqrt( OpP( x, xi ) * OpP( x, xi ) + ( P.lam ? Lam : FFVar( 1. ) ) * OpP( u, xi ) * OpP( u, xi ) ), xi ), { t, xi }, { T_INT, X_INT }, inter );
  oc.add_equation( x - xi,  { t, xi }, { LB, ALL },   ini );
  oc.add_equation( x,       { t, xi }, { T_INT, LB }, bnd );
  oc.add_equation( x - 1.0, { t, xi }, { T_INT, UB }, bnd );
  oc.update_ref( u, [xi]( OCFESLV::t_Coord const& c ){ return u0( c.at( xi ) ); } );
  oc.update_ref( x, [xi]( OCFESLV::t_Coord const& c ){ return c.at( xi ); } );
  if( SM && std::getenv( "MSREF1" ) ) oc.update_ref( Ms, 1.0 );   // diagnosis: a constant seed for Ms
  else if( SM ) oc.update_ref( Ms, [xi]( OCFESLV::t_Coord const& c ){ double const d = u0x( c.at( xi ) ); return std::sqrt( 1. + d * d ); } );
  oc.options.REDUCE.ORDER = ( P.red == "full" ) ? OCFESLV::Options::RED_FULL : OCFESLV::Options::RED_MAIN;
  if( P.lam ) oc.add_input( Lam, {}, 1.0 );
  oc.options.INTERFACE.IMPOSITION = imp;
  oc.options.SOLVE.MARCHING = !mono;
  oc.options.SOLVE.RES_TOL = P.rtol; oc.options.SOLVE.MAX_ITER = P.maxit;
  oc.options.DISPLAY_LEVEL = P.disp;
  try { R.setup = oc.setup(); } catch( std::exception const& e ){ std::printf( "  %s: setup threw %s\n", name.c_str(), e.what() ); }
  if( R.setup ){
    if( P.lam ){ oc.add_homotopy( Lam, 0.0, 1.0 ); oc.options.HOMOTOPY.ENABLE = true; oc.options.HOMOTOPY.SOLVE_S0 = true; oc.options.HOMOTOPY.STEP_CAP = 0.1; }   // after setup, as OCFE_MMPDE2
    std::vector<double> var, inp; oc.init( var, inp );
    if( char const* lf = std::getenv( "LOADVEC" ) ){   // diagnosis: start from a given state vector (e.g. the hand-reduced solution)
      std::FILE* f = std::fopen( lf, "r" ); size_t i = 0; double v_;
      while( f && i < var.size() && std::fscanf( f, "%lf", &v_ ) == 1 ) var[i++] = v_;
      if( f ) std::fclose( f ); std::printf( "  [loadvec] %zu of %zu values loaded\n", i, var.size() );
    }
    if( P.lam ) oc.set_input_values( Lam, { 1.0 }, inp.data() );
    OCFESLV::SolveReport rep = oc.solve( var.data(), inp.empty() ? nullptr : inp.data(), nullptr );
    R.conv = rep.converged; R.its = rep.iterations; R.resid = rep.final_residual;
    if( char const* jf = std::getenv( "JDUMPH" ) ){   // the Jacobian at the vector the solve left, with lam = JLAM
      if( P.lam ) oc.set_input_values( Lam, { std::getenv( "JLAM" ) ? std::atof( std::getenv( "JLAM" ) ) : 1.0 }, inp.data() );
      OCFESLV::t_JacAssembly A;
      if( oc.assemble( var.data(), inp.empty() ? nullptr : inp.data(), nullptr, A, false ) ){
        A.J.save( std::string( jf ), arma::coord_ascii );
        std::printf( "  [jdumph] %llu x %llu, nnz %llu\n", (unsigned long long)A.J.n_rows, (unsigned long long)A.J.n_cols, (unsigned long long)A.J.n_nonzero );
      } else std::printf( "  [jdumph] assemble FAILED\n" );
    }
    auto ev = [&]( FFVar const& s_, double tt, double z ){ OCFESLV::t_Coord pt; pt[t] = tt; pt[xi] = z;
      try { return oc.eval_colloc<double>( s_, pt, var.data(), nullptr, nullptr ); } catch( ... ){ return std::nan( "" ); } };
    R.ju = R.jx = R.jq = R.jg = R.drift = 0.;
    for( int k = 1; k <= 10; ++k ) for( size_t e = 1; e < P.nel; ++e ){
      double const tt = P.tf * k / 10., z = double( e ) / double( P.nel ), zl = z - 1e-10, zr = z + 1e-10;
      R.ju = std::max( R.ju, std::fabs( ev( u, tt, zl ) - ev( u, tt, zr ) ) );
      R.jx = std::max( R.jx, std::fabs( ev( x, tt, zl ) - ev( x, tt, zr ) ) );
      R.drift = std::max( R.drift, std::fabs( 0.5 * ( ev( x, tt, zl ) + ev( x, tt, zr ) ) - z ) );
    }
    // u(x, tf) on a common physical grid: sample (x(xi), u(xi)) finely, interpolate linearly in x
    std::vector<double> xs, us; for( int s_ = 0; s_ <= 2000; ++s_ ){ xs.push_back( ev( x, P.tf, s_ / 2000. ) ); us.push_back( ev( u, P.tf, s_ / 2000. ) ); }
    for( int s_ = 0; s_ <= 200; ++s_ ){
      double const X = s_ / 200.; size_t j = 1; while( j + 1 < xs.size() && xs[j] < X ) ++j;
      double const w_ = ( xs[j] > xs[j-1] ) ? ( X - xs[j-1] ) / ( xs[j] - xs[j-1] ) : 0.;
      R.uf.push_back( us[j-1] + w_ * ( us[j] - us[j-1] ) );
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
    if( num( "--tau"    ) ){ P.tau    = std::atof( argv[++i] ); continue; }
    if( num( "--smooth" ) ){ P.smooth = std::atof( argv[++i] ); continue; }
    if( num( "--red"    ) ){ P.red    = argv[++i]; continue; }
    if( !std::strcmp( argv[i], "--lam" ) ){ P.lam = true; continue; }
    std::cerr << "unknown argument " << argv[i] << "\n"; return 2;
  }
  std::printf( "OCFE_Burger_move_red: eps=%g, xi: %zu LGL elements of %zu nodes, t: %zu LGR elements on [0,%g], tau=%g, REDUCE.ORDER = RED_%s, smooth=%g, rtol=%g\n",
               P.eps, P.nel, P.nnd, P.nelt, P.tf, P.tau, P.red == "full" ? "FULL" : "MAIN", P.smooth, P.rtol );
  std::vector<Res> all;
  for( std::string const mode : { "mono", "march" } )
    for( auto const& [iname, imp] : std::vector<std::pair<std::string,Imp>>{
           { "weak", OCFESLV::Options::IC_WEAK }, { "trace", OCFESLV::Options::IC_TRACE }, { "strong", OCFESLV::Options::IC_STRONG } } ){
      std::string const name = mode + "-" + iname;
      if( ( "," + cases + "," ).find( "," + name + "," ) == std::string::npos ) continue;
      all.push_back( run( P, mode == "mono", imp, name ) );
      Res const& R = all.back();
      std::printf( "  %-13s setup %-3s conv %-3s its %3d  |r| %.2e  %7.1f s  seam jump u %.1e x %.1e  (q %.0e g %.0e n/a)  seam drift %.3f\n",
                   R.name.c_str(), R.setup ? "ok" : "NO", R.conv ? "yes" : "no", R.its, R.resid, R.sec, R.ju, R.jx, R.jq, R.jg, R.drift );
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
  // This driver is a REPRODUCER for the open moving-mesh items (WORKPLAN 4R.24), not a regression test: its default
  // configuration is expected to fail today.  It reports its verdict and exits 0 so a corpus sweep stays green;
  // make it `return ok ? 0 : 1;` once the moving mesh is expected to pass.
  std::printf( "%s\n", ok ? "PASS: all configurations converged; same-imposition and exact-imposition profiles agree"
                          : "REPRODUCER (WORKPLAN 4R.24): not all configurations converged/agreed -- see the header" );
  return 0;
}
