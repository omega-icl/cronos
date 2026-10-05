// ============================================================================
//  OCFE_capture_nested.cpp  --  consolidated coverage for NESTED evolution-
//  direction reductions in the captured-reduction form.  Supersedes
//  OCFE_march_nested.cpp (diagnostic) and OCFE_capture_nested_refuse.cpp.
//
//  Two regimes, distinguished by whether the inner reduction lives on the SAME
//  domain as the outer (evolution) reduction:
//
//   REFUSED (same evolution direction) -- the inner reduction's value is a
//   POST-solve captured input, so the outer's IN-solve materialised operand
//   would read a window-lagged value.  setup() must reject with
//   CAPTURE_NESTED_REFUSED (was silently wrong: 0.593 vs 1, 0.139 vs 1/4).
//     [1] INT_t ( x + INT_t' x dt' ) dt      (OpI in OpI)
//     [2] INT_t ( x * x(1/2) ) dt            (OpEval in OpI)
//
//   WORKS (inner on a SEPARATE, spatial domain) -- the inner INT_z is an
//   in-solve quantity (fctacc-summed output, or an aux/user STATE), available to
//   the outer evolution reduction.  No capture-into-in-solve-state, so no refusal.
//     [3] INT_t ( INT_z x dz ) dt              = 1/2   (linear; rides fctacc)
//     [4] S = INT_z x dz (state);  INT_t S dt  = 1/2   (outer captured over a state)
//                                  (INT_t S)^2 = 1/4   (nonlinear wrapper over the capture)
//
//  x(t)   = t          (1-D cases [1],[2];  x from dx/dt=1, x(0)=0)
//  x(t,z) = t          (2-D cases [3],[4];  z in [0,1] so INT_z x dz = t)
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
static void ok_flag( char const* nm, bool pass )
{
  if( !pass ) ++g_fail;
  std::cout << "    " << std::left << std::setw(40) << nm << std::right
            << "  " << ( pass ? "PASS" : "FAIL" ) << "\n";
}
static void ok_near( char const* nm, double got, double want, double tol = 1e-7 )
{
  bool const p = ( std::fabs( got - want ) < tol );
  if( !p ) ++g_fail;
  std::cout << "    " << std::left << std::setw(40) << nm << std::right
            << std::scientific << std::setprecision(8)
            << " got=" << std::setw(15) << got << " want=" << std::setw(15) << want
            << "  " << ( p ? "PASS" : "FAIL" ) << "\n";
}

// ---- shared 1-D evolution model:  x(t)=t on 3 march windows -------------------
static void build_1d( FFGraph& DAG, OCFESLV& oc, FFVar& t, FFVar& x )
{
  t = DAG.add_var( "t" );
  x = DAG.add_var( "x(t)" );
  FFPartial OpP;
  oc.set_evolution_domain( t );
  oc.add_domain( t, FFDom( 0., 1., 3, FFDom::LGR, 4 ) );
  oc.add_state ( x, { t } );
  oc.add_equation( OpP( x, t ) - 1., { t }, { FFDom::ALL - FFDom::LB },
                   OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( x, { t }, { FFDom::LB },
                   OCFESLV::EqnOptions( OCFESLV::EqnRole::INITIAL, 0 ) );
}

// ---- shared 2-D model:  x(t,z)=t,  t evolution (3 windows), z spatial ---------
static void build_2d( FFGraph& DAG, OCFESLV& oc, FFVar& t, FFVar& z, FFVar& x )
{
  t = DAG.add_var( "t" );
  z = DAG.add_var( "z" );
  x = DAG.add_var( "x(t,z)" );
  FFPartial OpP;
  oc.set_evolution_domain( t );
  oc.add_domain( t, FFDom( 0., 1., 3, FFDom::LGR, 4 ) );
  oc.add_domain( z, FFDom( 0., 1., 1, FFDom::LGR, 4 ) );
  oc.add_state ( x, { t, z } );
  oc.add_equation( OpP( x, t ) - 1., { t, z }, { FFDom::ALL - FFDom::LB, FFDom::ALL },
                   OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( x, { t, z }, { FFDom::LB, FFDom::ALL },
                   OCFESLV::EqnOptions( OCFESLV::EqnRole::INITIAL, 0 ) );
}

static void set_opts( OCFESLV& oc )
{
  oc.options.REDUCE.ORDER   = OCFESLV::Options::RED_FULL;
  oc.options.CLASSIFY.MODE       = OCFESLV::Options::CLASS_AUTO;
  oc.options.SOLVE.MARCHING = true;
  oc.options.DISPLAY_LEVEL  = 0;
}

// ============================================================================
//  [1] REFUSE: OpI inside OpI, same evolution direction
// ============================================================================
static void test_refuse_opi_in_opi()
{
  std::cout << "\n[1] REFUSE  INT_t ( x + INT_t' x ) dt  (OpI in OpI, same direction)\n";
  FFGraph DAG; OCFESLV oc( &DAG ); FFVar t, x; build_1d( DAG, oc, t, x );
  FFIntegral OpI;
  FFVar I_in = OpI( x, t );
  oc.add_output( OpI( x + I_in, t ) );
  set_opts( oc );
  bool const s = oc.setup();
  ok_flag( "setup rejected", !s );
  ok_flag( "status == CAPTURE_NESTED_REFUSED",
           oc.setup_status() == OCFESLV::SetupStatus::CAPTURE_NESTED_REFUSED );
}

// ============================================================================
//  [2] REFUSE: OpEval inside OpI, same evolution direction
// ============================================================================
static void test_refuse_opeval_in_opi()
{
  std::cout << "\n[2] REFUSE  INT_t ( x * x(1/2) ) dt  (OpEval in OpI, same direction)\n";
  FFGraph DAG; OCFESLV oc( &DAG ); FFVar t, x; build_1d( DAG, oc, t, x );
  FFIntegral OpI; FFEval OpEval;
  FFVar E_in = OpEval( x, t, 0.5 );
  oc.add_output( OpI( x * E_in, t ) );
  set_opts( oc );
  bool const s = oc.setup();
  ok_flag( "setup rejected", !s );
  ok_flag( "status == CAPTURE_NESTED_REFUSED",
           oc.setup_status() == OCFESLV::SetupStatus::CAPTURE_NESTED_REFUSED );
}

// ============================================================================
//  [3] WORKS: INT_t ( INT_z x dz ) dt  -- two domains, linear (fctacc path)
// ============================================================================
static void test_two_domain_linear()
{
  std::cout << "\n[3] WORKS   INT_t ( INT_z x dz ) dt = 1/2  (two domains, linear)\n";
  FFGraph DAG; OCFESLV oc( &DAG ); FFVar t, z, x; build_2d( DAG, oc, t, z, x );
  FFIntegral OpI;
  oc.add_output( OpI( OpI( x, z ), t ) );        // INT_t ( INT_z x dz ) dt
  set_opts( oc );
  if( !oc.setup() ){ ok_flag( "setup ok", false ); return; }
  ok_flag( "setup ok (not refused)", true );
  std::vector<double> vi, ii; oc.init( vi, ii, nullptr );
  std::vector<double> xv = vi;
  OCFESLV::SolveReport r = oc.solve( xv.data(), nullptr, nullptr );
  ok_flag( "converged", r.converged );
  std::vector<double> const& F = oc.val_functions();
  ok_near( "INT_t INT_z x == 1/2", F.empty() ? -1. : F[0], 0.5 );
}

// ============================================================================
//  [4] WORKS: S = INT_z x (state); INT_t S dt and (INT_t S)^2 -- capture path
// ============================================================================
static void test_two_domain_capture()
{
  std::cout << "\n[4] WORKS   S=INT_z x;  INT_t S dt = 1/2,  (INT_t S)^2 = 1/4  (two domains, capture)\n";
  FFGraph DAG; OCFESLV oc( &DAG ); FFVar t, z, x; build_2d( DAG, oc, t, z, x );
  FFIntegral OpI;
  FFVar S = DAG.add_var( "S(t)" );
  oc.add_state( S, { t } );
  oc.add_equation( S - OpI( x, z ), { t }, { FFDom::ALL },     // S(t) = INT_z x dz  (spatial reduction)
                   OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  FFVar It = OpI( S, t );                                      // outer evolution integral over a STATE
  oc.add_output( It      );                                    // linear   -> 1/2
  oc.add_output( It * It );                                    // nonlinear-> 1/4  (capture, not fctacc)
  set_opts( oc );
  if( !oc.setup() ){ ok_flag( "setup ok", false ); return; }
  ok_flag( "setup ok (not refused)", true );
  std::vector<double> vi, ii; oc.init( vi, ii, nullptr );
  std::vector<double> xv = vi;
  OCFESLV::SolveReport r = oc.solve( xv.data(), nullptr, nullptr );
  ok_flag( "converged", r.converged );
  std::vector<double> const& F = oc.val_functions();
  ok_near( "INT_t S     == 1/2", F.size() < 1 ? -1. : F[0], 0.5  );
  ok_near( "(INT_t S)^2 == 1/4", F.size() < 2 ? -1. : F[1], 0.25 );
}

int main()
{
  std::cout << "===============================================================\n"
            << "  nested reductions: same-direction REFUSED, two-domain WORKS\n"
            << "===============================================================\n";

  test_refuse_opi_in_opi();
  test_refuse_opeval_in_opi();
  test_two_domain_linear();
  test_two_domain_capture();

  std::cout << "\n  RESULT: " << ( g_fail == 0 ? "ALL PASS" : "FAIL" )
            << "  (" << g_fail << " failing checks)\n";
  return g_fail == 0 ? 0 : 1;
}
