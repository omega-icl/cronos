// ============================================================================
//  OCFE_mono_capture.cpp  --  gate for the MIXED extract-then-capture form of
//  evolution-direction reductions consumed by OUTPUTS (val_functions).
//
//  Marched state:   dx/dt = 1 , x(0)=0  ->  x(t)=t   (3 evolution windows)
//  Captured atoms:  I1 = INT_0^1 x   dt = 1/2
//                   I2 = INT_0^1 x^2 dt = 1/3
//
//  Outputs (read post-march via val_functions(), in add_output order):
//    [0]  I1          = 1/2              LINEAR   (matches the retired fctacc sum)
//    [1]  I2          = 1/3              LINEAR
//    [2]  I1*I1       = 1/4              NONLINEAR  <-- the discriminating check
//    [3]  I1*I2       = 1/6              NONLINEAR (product of two captures)
//
//  Why [2]/[3] matter: the OLD per-window fctacc sum would give
//    sum_k (INT_win_k x)^2 = (1/18)^2+(1/6)^2+(5/18)^2 = 35/324 ~= 0.108  (WRONG),
//  whereas the captured form gives (INT_0^1 x)^2 = 1/4 exactly.  So a PASS on [2]
//  proves the wrapper f(.) is evaluated ONCE over the completed capture, not summed
//  per window.  [0]/[1] confirm the linear case still lands on the horizon value.
//
//  The reduction is now a captured INPUT (no aux state), consumed only by outputs
//  read post-march -- the causality contract.  There is deliberately NO
//  eval_colloc(state) read of the reduction here (a state alias would lag).
//
//  Build: -DOCFE_OCFESLV_HEADER='"ocfeslv.hpp"'
// ============================================================================
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

#include OCFE_OCFESLV_HEADER

using namespace mc;

static bool g_ok = true;
static void check( char const* nm, double got, double want, double tol )
{
  double const e = std::fabs( got - want );
  bool const p = ( e < tol );
  g_ok &= p;
  std::cout << "  " << std::left << std::setw(28) << nm << std::right
            << std::scientific << std::setprecision(8)
            << " got=" << std::setw(15) << got << " want=" << std::setw(15) << want
            << std::setprecision(2) << " err=" << std::setw(9) << e
            << "  " << ( p ? "PASS" : "FAIL" ) << "\n";
}

int main()
{
  std::cout << "===============================================================\n"
            << "  mixed capture: outputs f(INT R dt) over MONOLITHIC t (SOLVE_MARCHING=false)\n"
            << "===============================================================\n";

  FFGraph DAG;
  FFVar t  = DAG.add_var( "t" );
  FFVar x  = DAG.add_var( "x(t)" );

  FFPartial  OpP;
  FFIntegral OpI;

  OCFESLV oc( &DAG );
  oc.set_evolution_domain( t );
  oc.add_domain( t, FFDom( 0., 1., 3, FFDom::LGR, 4 ) );   // 3 evolution elements -> 3 march windows
  oc.add_state ( x,  { t } );

  FFVar dxdt = OpP( x, t );
  oc.add_equation( dxdt - 1., { t }, { FFDom::ALL - FFDom::LB },
                   OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );   // dx/dt = 1  (interior + terminal)
  oc.add_equation( x, { t }, { FFDom::LB },
                   OCFESLV::EqnOptions( OCFESLV::EqnRole::INITIAL, 0 ) );    // x(0) = 0

  FFVar I1 = OpI( x,     t );    // INT x   dt  -> captured
  FFVar I2 = OpI( x * x, t );    // INT x^2 dt  -> captured
  oc.add_output( I1 );           // [0] linear
  oc.add_output( I2 );           // [1] linear
  oc.add_output( I1 * I1 );      // [2] NONLINEAR  (discriminates capture vs fctacc)
  oc.add_output( I1 * I2 );      // [3] NONLINEAR  (product of two captures)

  oc.options.REDUCE.ORDER   = OCFESLV::Options::RED_FULL;
  oc.options.CLASSIFY.MODE       = OCFESLV::Options::CLASS_AUTO;
  oc.options.SOLVE.MARCHING = false;
  oc.options.DISPLAY_LEVEL  = 1;

  if( !oc.setup() ){ std::cerr << "ERROR: setup() failed\n"; return 2; }
  std::cout << oc;    // dump reduced system: the OpI-over-t should appear as captured INPUTS, no aux state
  std::cout << "  marching=" << ( oc.is_marching() ? "YES" : "no (monolithic)" )
            << "  march_steps=" << oc.n_march_steps() << "\n";

  std::vector<double> varInit, inpInit;
  if( !oc.init( varInit, inpInit, nullptr ) ){ std::cerr << "ERROR: init() failed\n"; return 4; }

  std::vector<double> xv = varInit;
  OCFESLV::SolveReport rep = oc.solve( xv.data(), nullptr, nullptr );
  std::cout << "  solve: converged=" << ( rep.converged ? "yes" : "no" )
            << " iters=" << rep.iterations
            << " final|r|=" << std::scientific << std::setprecision(3) << rep.final_residual << "\n";
  if( !rep.converged ){ std::cerr << "ERROR: solve did not converge\n"; return 5; }

  std::vector<double> const& F = oc.val_functions();
  std::cout << "\n  val_functions().size()=" << F.size() << "\n";
  for( size_t i = 0; i < F.size(); ++i )
    std::cout << "    F[" << i << "]=" << std::scientific << std::setprecision(8) << F[i] << "\n";

  if( F.size() < 4 ){ std::cerr << "ERROR: expected >=4 output values\n"; return 6; }

  std::cout << "\n  checks:\n";
  check( "F0  INT x  dt      == 1/2",  F[0], 0.5,        1e-7 );
  check( "F1  INT x^2 dt     == 1/3",  F[1], 1.0/3.0,    1e-7 );
  check( "F2  (INT x dt)^2   == 1/4",  F[2], 0.25,       1e-7 );
  check( "F3  I1*I2          == 1/6",  F[3], 1.0/6.0,    1e-7 );

  std::cout << "\n  RESULT: " << ( g_ok ? "ALL PASS" : "FAIL" ) << "\n";
  return g_ok ? 0 : 1;
}
