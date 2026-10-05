// ============================================================================
//  OCFE_march_capture_cse.cpp  --  common-subexpression check for captured
//  evolution-direction reductions.
//
//  The SAME reduction, written SEPARATELY in several places, must collapse to a
//  SINGLE captured input -- both because the DAG hash-conses an identical
//  external op (lt_FFOp: type, external-op id, operands, then FFIntegral::lt /
//  FFEval::lt on the integration/evaluation data), and because reduce_order keys
//  its substitution map on the reduction NODE id (one node -> one cap, reused by
//  every equation/output that references it).
//
//  x(t)=t marched.  Four reductions are written as FOUR SEPARATE calls:
//     I1a = INT x dt , I1b = INT x dt        (identical -> one node)
//     E1a = x(1/2)   , E1b = x(1/2)          (identical -> one node)
//  Outputs:
//     [0] I1a          = 1/2
//     [1] I1b          = 1/2   (same cap as [0])
//     [2] I1a * I1b    = 1/4   (the SAME cap squared -- nonlinear, so captured)
//     [3] E1a * E1b    = 1/4   (the SAME latch, squared)
//
//  CSE assertion: n_colloc_inp() == 2  (one _Q for INT x dt, one _L for x(1/2))
//  -- four written reductions, two captured inputs.
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
static void check_int( char const* nm, long got, long want )
{
  bool const p = ( got == want );
  g_ok &= p;
  std::cout << "  " << std::left << std::setw(28) << nm << std::right
            << " got=" << got << " want=" << want << "  " << ( p ? "PASS" : "FAIL" ) << "\n";
}

int main()
{
  std::cout << "===============================================================\n"
            << "  CSE: separately-written identical reductions -> one capture\n"
            << "===============================================================\n";

  FFGraph DAG;
  FFVar t  = DAG.add_var( "t" );
  FFVar x  = DAG.add_var( "x(t)" );

  FFPartial  OpP;
  FFIntegral OpI;
  FFEval     OpEval;

  OCFESLV oc( &DAG );
  oc.set_evolution_domain( t );
  oc.add_domain( t, FFDom( 0., 1., 3, FFDom::LGR, 4 ) );
  oc.add_state ( x,  { t } );

  FFVar dxdt = OpP( x, t );
  oc.add_equation( dxdt - 1., { t }, { FFDom::ALL - FFDom::LB },
                   OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( x, { t }, { FFDom::LB },
                   OCFESLV::EqnOptions( OCFESLV::EqnRole::INITIAL, 0 ) );

  // FOUR separate reduction calls -- two identical integrals, two identical points.
  FFVar I1a = OpI   ( x, t      );
  FFVar I1b = OpI   ( x, t      );   // identical to I1a -> same DAG node
  FFVar E1a = OpEval( x, t, 0.5 );
  FFVar E1b = OpEval( x, t, 0.5 );   // identical to E1a -> same DAG node

  oc.add_output( I1a       );        // [0]
  oc.add_output( I1b       );        // [1]  same cap as [0]
  oc.add_output( I1a * I1b );        // [2]  (same cap)^2 -- nonlinear, forces capture
  oc.add_output( E1a * E1b );        // [3]  the same latch, squared

  oc.options.REDUCE.ORDER   = OCFESLV::Options::RED_FULL;
  oc.options.CLASSIFY.MODE       = OCFESLV::Options::CLASS_AUTO;
  oc.options.SOLVE.MARCHING = true;
  oc.options.DISPLAY_LEVEL  = 1;

  if( !oc.setup() ){ std::cerr << "ERROR: setup() failed\n"; return 2; }
  std::cout << oc;    // INPUTS line should list exactly one _Q and one _L
  std::cout << "  n_colloc_inp()=" << oc.n_colloc_inp() << "\n";

  std::vector<double> varInit, inpInit;
  if( !oc.init( varInit, inpInit, nullptr ) ){ std::cerr << "ERROR: init() failed\n"; return 4; }

  std::vector<double> xv = varInit;
  OCFESLV::SolveReport rep = oc.solve( xv.data(), nullptr, nullptr );
  if( !rep.converged ){ std::cerr << "ERROR: solve did not converge\n"; return 5; }

  std::vector<double> const& F = oc.val_functions();
  std::cout << "\n  val_functions().size()=" << F.size() << "\n";
  for( size_t i = 0; i < F.size(); ++i )
    std::cout << "    F[" << i << "]=" << std::scientific << std::setprecision(8) << F[i] << "\n";
  if( F.size() < 4 ){ std::cerr << "ERROR: expected >=4 output values\n"; return 6; }

  std::cout << "\n  checks:\n";
  check_int( "n_colloc_inp (CSE dedup)", (long)oc.n_colloc_inp(), 2 );  // 4 reductions -> 2 caps
  check( "F0  I1a           == 1/2", F[0], 0.5,  1e-7 );
  check( "F1  I1b           == 1/2", F[1], 0.5,  1e-7 );
  check( "F2  I1a*I1b       == 1/4", F[2], 0.25, 1e-7 );
  check( "F3  E1a*E1b       == 1/4", F[3], 0.25, 1e-7 );

  std::cout << "\n  RESULT: " << ( g_ok ? "ALL PASS" : "FAIL" ) << "\n";
  return g_ok ? 0 : 1;
}
