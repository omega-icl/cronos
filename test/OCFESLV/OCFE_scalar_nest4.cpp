// ============================================================================
//  OCFE_scalar_nest4.cpp  --  deeper reduction/derivative operand coverage
//
//  Exercises the paths nest3 did NOT: OpEval rebuilt around a compound
//  derivative operand (OpEvalLoc), a variable (state) coefficient, high
//  nonlinearity, and DEEP multi-pass peeling (a 2nd derivative inside a
//  nonlinear reduction operand).
//
//      a(z) = z(1-z) ,  a'(z)=1-2z ,  a''(z)=-2
//
//    t1 = OpEval( OpP(a,z)^2, z, 0.25 )      = (a'(0.25))^2 = 0.25   [ OpEvalLoc rebuild ]
//    t2 = OpI( a * OpP(a,z)^2, z )           = INT z(1-z)(1-2z)^2 dz = 1/30  [ STATE coeff x nonlinear ]
//    t3 = OpI( OpP(a,z)^4, z )               = INT (1-2z)^4 dz = 1/5  [ quartic of derivative ]
//    t4 = OpI( OpP(OpP(a,z),z)^2, z )        = INT (a'')^2 dz = 4      [ 2nd-deriv squared: deep peel ]
//
//  2026-09-18 -- A 3 x 2 MATRIX: impositions WEAK / TRACE / STRONG  x  reductions RED_FULL / RED_MAIN.
//  The driver was IC_STRONG + RED_FULL only, and that hid a defect for as long as it existed: under RED_FULL the
//  balance row a'' + 2 = 0 is rewritten as  Dpz_Dz_a + 2 = 0  -- ALGEBRAIC in the deepest auxiliary, no
//  derivative left -- so the principal symbol is empty and NO interface claim gets a balance receiver ([a] and
//  [Dz_a] land in the chain's LINK rows, [Dpz_Dz_a] only has refused rescued edges).  The two interface
//  conditions [a]=0, [a']=0 are load-bearing (a''=-2 at 5 points per element leaves one linear mode free per
//  element; the redundant rows are the PDE rows at the two interface copies).  IC_STRONG drops those rows and
//  imposes the constraints hard: right to 1e-13.  IC_WEAK penalises LINK rows, which are not the redundant ones:
//  rank 48/49 (STATE CONTENT), no convergence.  IC_TRACE puts the multipliers there too, B has rank 2 of 3, the
//  projector springs the dependent combination: t1-t3 wrong by 3-7%.  Under RED_MAIN the balance row keeps
//  Dz_a' + 2 and every mode is exact to 1e-16 (NOTES_20260918a, T-P1: the chain-aware balance receiver).
//  Until the chain-aware symbol is the default: RED_FULL x {WEAK, TRACE} are documented XFAIL and the other four
//  cells must PASS; with CRONOS_CHAIN_SYMBOL=1 (rev283) all six must PASS.
//
//  Build:
//    g++ -std=c++17 ... -DOCFE_OCFESLV_HEADER='"ocfeslv.hpp"'
//        OCFE_scalar_nest4.cpp -o OCFE_scalar_nest4  <libs>
// ============================================================================
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include OCFE_OCFESLV_HEADER

using namespace mc;

static bool g_ok = true;
static void check( char const* nm, double got, double want, double tol )
{
  double const err = std::fabs( got - want );
  bool const pass = ( err < tol );
  g_ok &= pass;
  std::cout << "  " << std::left << std::setw(32) << nm
            << std::right << std::scientific << std::setprecision(8)
            << " got=" << std::setw(15) << got << " want=" << std::setw(15) << want
            << std::setprecision(2) << "  err=" << std::setw(9) << err
            << "  " << ( pass ? "PASS" : "FAIL" ) << "\n";
}


typedef OCFESLV::Options O;
struct CellR { std::string red, imp; bool ok=false; std::string verdict; };

static CellR run_cell( O::ReductionType red, O::ImpositionType imp )
{
  CellR C; C.red = ( red==O::RED_FULL ? "RED_FULL" : "RED_MAIN" );
  C.imp = ( imp==O::IC_WEAK ? "WEAK" : imp==O::IC_TRACE ? "TRACE" : "STRONG" );
  std::cout << "\n---- " << C.red << " / " << C.imp << " ----\n";
  g_ok = true;

  FFGraph DAG;
  FFVar z  = DAG.add_var( "z" );
  FFVar a  = DAG.add_var( "a(z)" );
  FFVar t1 = DAG.add_var( "t1" );
  FFVar t2 = DAG.add_var( "t2" );
  FFVar t3 = DAG.add_var( "t3" );
  FFVar t4 = DAG.add_var( "t4" );
  FFPartial  OpP;
  FFEval     OpEval;
  FFIntegral OpI;

  OCFESLV oc( &DAG );
  oc.add_domain( z, FFDom( 0., 1., 2, FFDom::CGL, 6 ) );
  oc.add_state ( a,  { z } );
  oc.add_state ( t1, {} );
  oc.add_state ( t2, {} );
  oc.add_state ( t3, {} );
  oc.add_state ( t4, {} );

  int const INT = FFDom::ALL - FFDom::LB - FFDom::UB;

  FFVar PDE_a = OpP( OpP( a, z ), z ) + 2.;   // a'' = -2 -> a = z(1-z)
  FFVar Da  = OpP( a, z );
  FFVar D2a = OpP( OpP( a, z ), z );

  FFVar E1 = t1 - OpEval( Da * Da, z, 0.25 );
  FFVar E2 = t2 - OpI( a * ( Da * Da ), z );
  FFVar E3 = t3 - OpI( ( Da * Da ) * ( Da * Da ), z );
  FFVar E4 = t4 - OpI( D2a * D2a, z );

  oc.add_equation( PDE_a, { z }, { INT }, OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( a, { z }, { FFDom::LB }, OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );
  oc.add_equation( a, { z }, { FFDom::UB }, OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );
  oc.add_equation( E1, {}, {}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( E2, {}, {}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( E3, {}, {}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( E4, {}, {}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );

  oc.options.REDUCE.ORDER    = red;
  oc.options.CLASSIFY.MODE        = OCFESLV::Options::CLASS_AUTO;
  oc.options.INTERFACE.IMPOSITION = imp;
  oc.options.INTERFACE.TYPE  = OCFESLV::Options::IC_AUTO;
  oc.options.INTERFACE.SAT_SIGMA0      = 1.0;
  oc.options.DISPLAY_LEVEL   = 1;

  if( !oc.setup() ){ std::cerr << "ERROR: setup() failed\n"; C.verdict = "FAIL(setup)"; return C; }
  std::cout << "  setup: states=" << oc.n_colloc_sta() << " equations=" << oc.n_colloc_eqn()
            << " square=" << ( oc.n_colloc_sta() == oc.n_colloc_eqn() ? "yes" : "NO" ) << "\n";
  if( oc.n_colloc_sta() != oc.n_colloc_eqn() ){ C.verdict = "FAIL(not square)"; return C; }
  std::vector<double> varInit, inpInit;
  if( !oc.init( varInit, inpInit, nullptr ) ){ C.verdict = "FAIL(init)"; return C; }
  std::vector<double> xv = varInit;
  OCFESLV::SolveReport rep = oc.solve( xv.data(), nullptr, nullptr );
  std::cout << "  solve: converged=" << ( rep.converged ? "yes" : "no" ) << " iters=" << rep.iterations
            << " final|r|=" << std::scientific << std::setprecision(3) << rep.final_residual << "\n";
  if( !rep.converged ){ C.verdict = "FAIL(no convergence)"; return C; }
  OCFESLV::t_Coord sp;
  double const v1 = oc.eval_colloc<double>( t1, sp, xv.data(), nullptr, nullptr );
  double const v2 = oc.eval_colloc<double>( t2, sp, xv.data(), nullptr, nullptr );
  double const v3 = oc.eval_colloc<double>( t3, sp, xv.data(), nullptr, nullptr );
  double const v4 = oc.eval_colloc<double>( t4, sp, xv.data(), nullptr, nullptr );
  std::cout << "  checks:\n";
  check( "t1 OpEval(OpP(a,z)^2,z,0.25)==1/4", v1, 0.25,      1e-7 );
  check( "t2 INT a*(a')^2 dz        ==1/30", v2, 1.0/30.0,  1e-7 );
  check( "t3 INT (a')^4 dz          ==1/5",  v3, 0.2,       1e-7 );
  check( "t4 INT (a'')^2 dz         ==4",    v4, 4.0,       1e-6 );
  C.ok = g_ok; C.verdict = g_ok ? "PASS" : "FAIL";
  return C;
}

int main()
{
  std::cout << "================================================================\n"
            << "  deeper reduction/derivative operand coverage -- 3 impositions x 2 reductions\n"
            << "================================================================\n";
  O::ReductionType const REDS[2] = { O::RED_FULL, O::RED_MAIN };
  O::ImpositionType const IMPS[3] = { O::IC_WEAK, O::IC_TRACE, O::IC_STRONG };
  std::vector<CellR> cells;
  for( auto red : REDS ) for( auto imp : IMPS ) cells.push_back( run_cell( red, imp ) );
  std::cout << "\n================================ SUMMARY ================================\n";
  bool all = true;
  for( auto& C : cells ){
    // T-P1: RED_FULL under WEAK and TRACE is the documented defect (see the header); expected to fail until
    // the chain-aware balance receiver exists.  A PASS there is a finding (the defect is gone) and is reported.
    // rev283 (CRONOS_CHAIN_SYMBOL): with the chain-aware symbol the balance row keeps its derivative column and
    // these two cells must PASS; without it they are the documented defect.  The expectation follows the knob.
    static bool const chain_symbol = !( std::getenv( "CRONOS_CHAIN_SYMBOL" ) && !std::atoi( std::getenv( "CRONOS_CHAIN_SYMBOL" ) ) );   // rev284: on unless =0
    bool const xfail = !chain_symbol && ( C.red == "RED_FULL" && ( C.imp == "WEAK" || C.imp == "TRACE" ) );
    std::string v = C.ok ? ( xfail ? "PASS (unexpected: T-P1 defect absent?)" : "PASS" )
                         : ( xfail ? "XFAIL (T-P1: RED_FULL balance row algebraic in the deepest auxiliary)" : C.verdict );
    if( !C.ok && !xfail ) all = false;
    std::cout << "  " << std::left << std::setw(10) << C.red << std::setw(8) << C.imp << v << "\n";
  }
  std::cout << "\nOCFE_scalar_nest4: " << ( all ? "PASS" : "FAIL" ) << "\n";
  return all ? 0 : 1;
}
