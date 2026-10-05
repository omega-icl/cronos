// ============================================================================
//  OCFE_march_capture_eval.cpp  --  gate for LATCH capture of evolution-
//  direction POINT reductions (OpEval) consumed by OUTPUTS (val_functions),
//  the OpEval analogue of OCFE_march_capture.cpp.
//
//  Marched state:  dx/dt = 1 , x(0)=0  ->  x(t)=t   (3 evolution windows)
//  Captured atoms (interior points, latched at their window):
//     E1 = x(1/2)      -> 1/2   (interior, window 1 of [1/3,2/3])
//     E2 = x(1/3)      -> 1/3   (grid node -> left-window tie-break)
//  And one MIXED output combining an OpI capture with an OpEval capture.
//
//  Outputs (read post-march via val_functions(), in add_output order):
//    [0]  x(1/2)              = 1/2      LINEAR point
//    [1]  x(1/3)              = 1/3      LINEAR point (boundary tie-break)
//    [2]  (x(1/2))^2          = 1/4      NONLINEAR point
//    [3]  x(1/2) * INT x dt   = 1/4      MIXED  (OpEval capture * OpI capture)
//
//  An interior-tau OpEval could NOT be served by an aux state w=x(tau) under
//  marching (the state is unsatisfiable in windows whose range excludes tau);
//  the LATCH capture writes cap=x(tau) at tau's window and holds it, so the
//  enclosing output f(cap) is correct post-march.  Read via val_functions();
//  deliberately no eval_colloc(state) read of the reduction.
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
            << "  latch capture: outputs f(x(tau)) over marched t\n"
            << "===============================================================\n";

  FFGraph DAG;
  FFVar t  = DAG.add_var( "t" );
  FFVar x  = DAG.add_var( "x(t)" );

  FFPartial  OpP;
  FFIntegral OpI;
  FFEval     OpEval;

  OCFESLV oc( &DAG );
  oc.set_evolution_domain( t );
  oc.add_domain( t, FFDom( 0., 1., 3, FFDom::LGR, 4 ) );   // 3 evolution elements -> 3 march windows
  oc.add_state ( x,  { t } );

  FFVar dxdt = OpP( x, t );
  oc.add_equation( dxdt - 1., { t }, { FFDom::ALL - FFDom::LB },
                   OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );   // dx/dt = 1  (interior + terminal)
  oc.add_equation( x, { t }, { FFDom::LB },
                   OCFESLV::EqnOptions( OCFESLV::EqnRole::INITIAL, 0 ) );    // x(0) = 0

  FFVar E1 = OpEval( x, t, 0.5     );    // x(1/2)  -> latched at window 1
  FFVar E2 = OpEval( x, t, 1.0/3.0 );    // x(1/3)  -> boundary, left-window tie-break
  FFVar I1 = OpI   ( x, t          );    // INT x dt -> captured (accumulated)
  oc.add_output( E1        );            // [0] linear point
  oc.add_output( E2        );            // [1] linear point (tie-break)
  oc.add_output( E1 * E1   );            // [2] NONLINEAR point
  oc.add_output( E1 * I1   );            // [3] MIXED: OpEval capture * OpI capture

  oc.options.REDUCE.ORDER   = OCFESLV::Options::RED_FULL;
  oc.options.CLASSIFY.MODE       = OCFESLV::Options::CLASS_AUTO;
  oc.options.SOLVE.MARCHING = true;
  oc.options.DISPLAY_LEVEL  = 1;

  if( !oc.setup() ){ std::cerr << "ERROR: setup() failed\n"; return 2; }
  std::cout << oc;    // dump: E1/E2 should appear as captured _L inputs, I1 as a _Q input
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
  check( "F0  x(1/2)         == 1/2",  F[0], 0.5,      1e-7 );
  check( "F1  x(1/3)         == 1/3",  F[1], 1.0/3.0,  1e-7 );
  check( "F2  (x(1/2))^2     == 1/4",  F[2], 0.25,     1e-7 );
  check( "F3  x(1/2)*INT x   == 1/4",  F[3], 0.25,     1e-7 );

  std::cout << "\n  RESULT: " << ( g_ok ? "ALL PASS" : "FAIL" ) << "\n";
  return g_ok ? 0 : 1;
}
