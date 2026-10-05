// ============================================================================
//  OCFE_mono_eval_direct.cpp  --  isolate the sweep regression:
//  a BARE (non-captured) evolution-integral output read via a DIRECT oc.eval()
//  call in MONOLITHIC, exactly as PDE32b/PDE36b/PSA read Eff/F1_in.
//
//  x(t)=t  (dx/dt=1, x(0)=0).  Bare evolution integral INT_0^1 x dt = 1/2.
//  Case A: 1-D, scalar output add_output(INT x dt),          read via oc.eval().
//  Case B: 2-D (t,z), output add_output(INT_t x dt, {z},{1}), read via oc.eval()
//          -- x independent of z, so INT_t x |z=1 = 1/2 (matches driver shape:
//          evolution integral presented at a spatial point).
//
//  Both are LINEAR/bare so the linearity gate must leave them on fctacc/_eval_fct
//  (NOT captured).  Expected 1/2 in each.  A 0 here reproduces the driver bug in
//  isolation and localises it to _eval_fct's monolithic evolution-integral sum.
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
static void ck( char const* nm, double got, double want )
{
  bool const p = std::fabs( got - want ) < 1e-7;
  if( !p ) ++g_fail;
  std::cout << "    " << std::left << std::setw(34) << nm << std::right
            << std::scientific << std::setprecision(8)
            << " got=" << std::setw(15) << got << " want=" << std::setw(15) << want
            << "  " << ( p ? "PASS" : "FAIL" ) << "\n";
}

static double case_A()
{
  std::cout << "\n[A] 1-D bare INT x dt, monolithic, DIRECT oc.eval()\n";
  FFGraph DAG;
  FFVar t = DAG.add_var("t"), x = DAG.add_var("x(t)");
  FFPartial OpP; FFIntegral OpI;
  OCFESLV oc( &DAG );
  oc.set_evolution_domain( t );
  oc.add_domain( t, FFDom( 0., 1., 3, FFDom::LGR, 4 ) );
  oc.add_state ( x, { t } );
  oc.add_equation( OpP(x,t) - 1., { t }, { FFDom::ALL - FFDom::LB }, OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( x, { t }, { FFDom::LB }, OCFESLV::EqnOptions( OCFESLV::EqnRole::INITIAL, 0 ) );
  oc.add_output( OpI( x, t ) );                        // bare evolution integral (scalar)
  oc.options.REDUCE.ORDER = OCFESLV::Options::RED_FULL;
  oc.options.SOLVE.MARCHING = false;
  oc.options.DISPLAY_LEVEL = 0;
  if( !oc.setup() ){ std::cout << "    setup FAILED\n"; ++g_fail; return -1; }
  std::cout << oc;                                     // dump: OpI should NOT be a captured _Q input
  std::vector<double> vi, ii; oc.init( vi, ii, nullptr );
  std::vector<double> xv = vi;
  oc.solve( xv.data(), ii.empty()?nullptr:ii.data(), nullptr );
  std::vector<double> fct( oc.n_colloc_fct(), 0. ), eqn( oc.n_colloc_eqn(), 0. );
  oc.eval( eqn.data(), fct.data(), xv.data(), ii.empty()?nullptr:ii.data(), nullptr );  // DIRECT eval
  return fct.empty()? -1. : fct[0];
}

static double case_B()
{
  std::cout << "\n[B] 2-D bare INT_t x dt |z=1, monolithic, DIRECT oc.eval()\n";
  FFGraph DAG;
  FFVar t = DAG.add_var("t"), z = DAG.add_var("z"), x = DAG.add_var("x(t,z)");
  FFPartial OpP; FFIntegral OpI;
  OCFESLV oc( &DAG );
  oc.set_evolution_domain( t );
  oc.add_domain( t, FFDom( 0., 1., 3, FFDom::LGR, 4 ) );
  oc.add_domain( z, FFDom( 0., 1., 1, FFDom::LGR, 4 ) );
  oc.add_state ( x, { t, z } );
  oc.add_equation( OpP(x,t) - 1., { t, z }, { FFDom::ALL - FFDom::LB, FFDom::ALL }, OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( x, { t, z }, { FFDom::LB, FFDom::ALL }, OCFESLV::EqnOptions( OCFESLV::EqnRole::INITIAL, 0 ) );
  oc.add_output( OpI( x, t ), { z }, { 1.0 } );        // evolution integral presented at z=1 (like Eff)
  oc.options.REDUCE.ORDER = OCFESLV::Options::RED_FULL;
  oc.options.SOLVE.MARCHING = false;
  oc.options.DISPLAY_LEVEL = 0;
  if( !oc.setup() ){ std::cout << "    setup FAILED\n"; ++g_fail; return -1; }
  std::cout << oc;
  std::vector<double> vi, ii; oc.init( vi, ii, nullptr );
  std::vector<double> xv = vi;
  oc.solve( xv.data(), ii.empty()?nullptr:ii.data(), nullptr );
  std::vector<double> fct( oc.n_colloc_fct(), 0. ), eqn( oc.n_colloc_eqn(), 0. );
  oc.eval( eqn.data(), fct.data(), xv.data(), ii.empty()?nullptr:ii.data(), nullptr );  // DIRECT eval
  return fct.empty()? -1. : fct[0];
}

static double case_C()
{
  std::cout << "\n[C] 2-D bare INT_t x dt |z=1, 5 z-elements, monolithic, DIRECT oc.eval()\n";
  FFGraph DAG;
  FFVar t = DAG.add_var("t"), z = DAG.add_var("z"), x = DAG.add_var("x(t,z)");
  FFPartial OpP; FFIntegral OpI;
  OCFESLV oc( &DAG );
  oc.set_evolution_domain( t );
  oc.add_domain( t, FFDom( 0., 1., 3, FFDom::LGR, 4 ) );
  oc.add_domain( z, FFDom( 0., 1., 5, FFDom::LGR, 4 ) );   // 5 z-elements (driver shape)
  oc.add_state ( x, { t, z } );
  oc.add_equation( OpP(x,t) - 1., { t, z }, { FFDom::ALL - FFDom::LB, FFDom::ALL }, OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( x, { t, z }, { FFDom::LB, FFDom::ALL }, OCFESLV::EqnOptions( OCFESLV::EqnRole::INITIAL, 0 ) );
  oc.add_output( OpI( x, t ), { z }, { 1.0 } );            // evolution integral at z=1 (upper boundary)
  oc.options.REDUCE.ORDER = OCFESLV::Options::RED_FULL;
  oc.options.SOLVE.MARCHING = false;
  oc.options.DISPLAY_LEVEL = 0;
  if( !oc.setup() ){ std::cout << "    setup FAILED\n"; ++g_fail; return -1; }
  std::cout << oc;
  std::vector<double> vi, ii; oc.init( vi, ii, nullptr );
  std::vector<double> xv = vi;
  oc.solve( xv.data(), ii.empty()?nullptr:ii.data(), nullptr );
  std::vector<double> fct( oc.n_colloc_fct(), 0. ), eqn( oc.n_colloc_eqn(), 0. );
  oc.eval( eqn.data(), fct.data(), xv.data(), ii.empty()?nullptr:ii.data(), nullptr );  // DIRECT eval
  return fct.empty()? -1. : fct[0];
}

int main()
{
  std::cout << "===============================================================\n"
            << "  monolithic bare evolution integral via DIRECT oc.eval()\n"
            << "===============================================================\n";
  double a = case_A();  ck( "A INT x dt == 1/2", a, 0.5 );
  double b = case_B();  ck( "B INT_t x |z=1 == 1/2", b, 0.5 );
  double c = case_C();  ck( "C INT_t x |z=1 (5 z-el) == 1/2", c, 0.5 );
  std::cout << "\n  RESULT: " << ( g_fail==0 ? "ALL PASS" : "FAIL" ) << "  (" << g_fail << " failing)\n";
  return g_fail==0 ? 0 : 1;
}
