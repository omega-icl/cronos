// ============================================================================
//  OCFE_capture_marching_eval.cpp  --  Approach A: a DIRECT eval() of captured
//  outputs on a MARCHING model returns the same values as val_functions().
//  Two models: SCALAR cap (Bolza x=a t, I1=INT x) and DISTRIBUTED cap
//  (c=z(1+a t), G=INT_t c read at z=1), each with a linear and a nonlinear
//  captured output.  After a marching solve, val_functions() (the marched
//  reader) and oc.eval(...) (direct) must agree; a fresh marching solve at a
//  different a must refresh the cache (no stale caps).
//
//  Build: -DOCFE_OCFESLV_HEADER='"ocfeslv.hpp"'
// ============================================================================
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

#include OCFE_OCFESLV_HEADER

using namespace mc;

static int g_fail = 0;
static void ck( char const* nm, double got, double want, double tol )
{
  bool const p = std::fabs( got - want ) < tol;
  if( !p ) ++g_fail;
  std::cout << "    " << std::left << std::setw(34) << nm << std::right
            << std::scientific << std::setprecision(6)
            << " got=" << std::setw(14) << got << " want=" << std::setw(14) << want
            << "  " << ( p ? "PASS" : "FAIL" ) << "\n";
}

// ---- scalar model: dx/dt=a, x(0)=0 => x=a t ; I1=INT x=a/2 ---------------------
static void build_scalar( FFGraph& DAG, OCFESLV& oc, FFVar& t, FFVar& x, FFVar& a )
{
  t = DAG.add_var("t"); x = DAG.add_var("x(t)"); a = DAG.add_var("a");
  FFPartial OpP; FFIntegral OpI;
  oc.set_evolution_domain( t );
  oc.add_domain( t, FFDom( 0., 1., 3, FFDom::LGR, 4 ) );
  oc.add_state ( x, { t } );
  oc.add_input ( a, {} );
  oc.add_equation( OpP(x,t) - a, { t }, { FFDom::ALL - FFDom::LB }, OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( x, { t }, { FFDom::LB }, OCFESLV::EqnOptions( OCFESLV::EqnRole::INITIAL, 0 ) );
  FFVar I1 = OpI( x, t );
  oc.add_output( I1      );   // linear cap
  oc.add_output( I1 * I1 );   // nonlinear cap
  oc.options.REDUCE.ORDER   = OCFESLV::Options::RED_FULL;
  oc.options.SOLVE.MARCHING = true;
  oc.options.DISPLAY_LEVEL  = 0;
}

// ---- distributed model: c=z(1+a t); G=INT_t c = z(1+a/2), read at z=1 ----------
static double const Dz = 0.1;
static void build_dist( FFGraph& DAG, OCFESLV& oc, FFVar& t, FFVar& z, FFVar& c, FFVar& a )
{
  t = DAG.add_var("t"); z = DAG.add_var("z"); c = DAG.add_var("c(t,z)"); a = DAG.add_var("a");
  FFPartial OpP; FFIntegral OpI;
  oc.add_domain( t, FFDom( 0., 1., 3, FFDom::LGR, 4 ) );
  oc.add_domain( z, FFDom( 0., 1., 2, FFDom::LGL, 4 ) );
  oc.add_state ( c, { t, z } );
  oc.add_input ( a, {} );
  OCFESLV::EqnOptions io( OCFESLV::EqnRole::INTERIOR, 0 ), ii( OCFESLV::EqnRole::INITIAL, 0 ), ib( OCFESLV::EqnRole::BOUNDARY, 0 );
  int const T_INT = FFDom::ALL - FFDom::LB, Z_INT = FFDom::ALL - FFDom::LB - FFDom::UB;
  oc.add_equation( OpP(c,t) - Dz*OpP(c,{z,2}) - a*z, { t, z }, { T_INT, Z_INT }, io );
  oc.add_equation( c - z,                            { t, z }, { FFDom::LB, FFDom::ALL }, ii );
  oc.add_equation( c,                                { t, z }, { T_INT, FFDom::LB }, ib );
  oc.add_equation( OpP(c,z) - ( 1.0 + a*t ),         { t, z }, { T_INT, FFDom::UB }, ib );
  FFVar G = OpI( c, t );
  oc.add_output( G,     { z }, { 1.0 } );   // linear cap
  oc.add_output( G * G, { z }, { 1.0 } );   // nonlinear cap
  oc.set_evolution_domain( t );
  oc.options.REDUCE.ORDER   = OCFESLV::Options::RED_FULL;
  oc.options.SOLVE.MARCHING = true;
  oc.options.DISPLAY_LEVEL  = 0;
}

// Explicit drivers:
static void run_scalar( double aval )
{
  std::cout << "\n[SCALAR]  a=" << aval << "\n";
  FFGraph DAG; OCFESLV oc( &DAG ); FFVar t,x,a; build_scalar( DAG, oc, t, x, a );
  if( !oc.setup() ){ std::cout << "    setup FAILED\n"; ++g_fail; return; }
  std::vector<double> vi, ii; oc.init( vi, ii, nullptr );
  oc.set_input_values( a, { aval }, ii.data() );
  std::vector<double> xv = vi;
  if( !oc.solve( xv.data(), ii.data(), nullptr ).converged ){ std::cout << "    solve FAILED\n"; ++g_fail; return; }
  std::vector<double> const val = oc.val_functions();               // marched reader
  std::vector<double> eqn( oc.n_colloc_eqn(), 0. ), fct( oc.n_colloc_fct(), 0. );
  bool const ok = oc.eval( eqn.data(), fct.data(), xv.data(), ii.empty()?nullptr:ii.data(), nullptr );
  if( !ok || val.size() < 2 || fct.size() < 2 ){ std::cout << "    eval FAILED\n"; ++g_fail; return; }
  double const I1 = aval/2.;
  ck( "val_functions I1",           val[0], I1,      1e-6 );
  ck( "direct eval I1 == val",      fct[0], val[0],  1e-9 );
  ck( "direct eval I1^2 == val",    fct[1], val[1],  1e-9 );
  ck( "direct eval I1^2 correct",   fct[1], I1*I1,   1e-6 );
}

static void run_dist( double aval )
{
  std::cout << "\n[DISTRIBUTED]  a=" << aval << "\n";
  FFGraph DAG; OCFESLV oc( &DAG ); FFVar t,z,c,a; build_dist( DAG, oc, t, z, c, a );
  if( !oc.setup() ){ std::cout << "    setup FAILED\n"; ++g_fail; return; }
  std::vector<double> vi, ii; oc.init( vi, ii, nullptr );
  oc.set_input_values( a, { aval }, ii.data() );
  std::vector<double> xv = vi;
  if( !oc.solve( xv.data(), ii.data(), nullptr ).converged ){ std::cout << "    solve FAILED\n"; ++g_fail; return; }
  std::vector<double> const val = oc.val_functions();
  std::vector<double> eqn( oc.n_colloc_eqn(), 0. ), fct( oc.n_colloc_fct(), 0. );
  bool const ok = oc.eval( eqn.data(), fct.data(), xv.data(), ii.empty()?nullptr:ii.data(), nullptr );
  if( !ok || val.size() < 2 || fct.size() < 2 ){ std::cout << "    eval FAILED\n"; ++g_fail; return; }
  double const Glin = 1.0 + aval/2.;
  ck( "val_functions Glin",         val[0], Glin,    1e-6 );
  ck( "direct eval Glin == val",    fct[0], val[0],  1e-9 );
  ck( "direct eval Gsq  == val",    fct[1], val[1],  1e-9 );
  ck( "direct eval Gsq  correct",   fct[1], Glin*Glin, 1e-6 );
}

int main()
{
  std::cout << "===============================================================\n"
            << "  MARCHING direct eval() of captured outputs == val_functions()\n"
            << "===============================================================\n";
  run_scalar( 1.2 );
  run_scalar( 2.0 );   // refresh cache at a new point (no stale caps)
  run_dist  ( 1.2 );
  run_dist  ( 2.0 );
  std::cout << "\n  RESULT: " << ( g_fail==0 ? "ALL PASS" : "FAIL" ) << "  (" << g_fail << " failing)\n";
  return g_fail==0 ? 0 : 1;
}
