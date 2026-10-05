// ============================================================================
//  OCFE_scalar_smoke.cpp
//
//  Purely-scalar smoke test for OCFESLV Stage 1 (native lumped algebraic states).
//
//  A standard steady elliptic BVP provides a genuine distributed state so the
//  scalar rows are exercised *alongside* distributed ones:
//
//      d2c/dz2 = -2 ,  c(0)=0 , c(1)=0     ->   c(z) = z(1-z)      [distributed]
//
//  plus two LUMPED (empty-domain) scalar algebraic states, closed by scalar-only
//  equations (no distributed coupling -- that FFIntegral case is Stage 2/section 4):
//
//      s1 - 42      = 0                    ->   s1 = 42            [scalar]
//      s2 - 3 * s1  = 0                    ->   s2 = 126           [scalar, scalar-coupled]
//
//  What this validates for Stage 1:
//    * add_state(s,{}) + add_equation(SCi,{},{}) assemble a square system
//      (the LOOP-2 empty-domain row idiom actually emits the scalar rows);
//    * classification handles lumped algebraic states coexisting with a
//      differential distributed state;
//    * the solve converges and recovers s1=42, s2=126, c(0.5)=0.25;
//    * eval_colloc() and eval_solution() both read scalar states correctly.
//
//  Build (against the Stage-1 header):
//    g++ -std=c++17 ... -DOCFE_OCFESLV_HEADER='"ocfeslv_scalar_stage1r.hpp"' \
//        OCFE_scalar_smoke.cpp -o OCFE_scalar_smoke  <libs>
// ============================================================================
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>

#ifndef OCFE_OCFESLV_HEADER
#define OCFE_OCFESLV_HEADER "ocfeslv.hpp"
#endif
#include OCFE_OCFESLV_HEADER

using namespace mc;

static bool g_ok = true;
static void check( char const* nm, double got, double want, double tol )
{
  double const err = std::fabs( got - want );
  bool const pass = ( err < tol );
  g_ok &= pass;
  std::cout << "  " << std::left << std::setw(30) << nm
            << std::right << std::scientific << std::setprecision(8)
            << " got=" << std::setw(15) << got
            << " want=" << std::setw(15) << want
            << std::setprecision(2) << "  err=" << std::setw(9) << err
            << "  " << ( pass ? "PASS" : "FAIL" ) << "\n";
}

int main()
{
  std::cout << "================================================================\n"
            << "  OCFESLV Stage 1 smoke test: lumped scalar algebraic states\n"
            << "================================================================\n";

  FFGraph DAG;
  FFVar z  = DAG.add_var( "z" );
  FFVar c  = DAG.add_var( "c(z)" );
  FFVar s1 = DAG.add_var( "s1" );      // lumped scalar
  FFVar s2 = DAG.add_var( "s2" );      // lumped scalar

  FFPartial OpP;

  OCFESLV oc( &DAG );
  oc.add_domain( z, FFDom( 0., 1., 2, FFDom::CGL, 6 ) );
  oc.add_state ( c,  { z } );          // distributed differential state
  oc.add_state ( s1, {}    );          // scalar (empty domain)
  oc.add_state ( s2, {}    );          // scalar (empty domain)

  // Distributed: c'' + 2 = 0 on the interior, c=0 at both ends.
  FFVar PDE   = OpP( OpP( c, z ), z ) + 2.;
  FFVar BC_lo = c;                     // c(0) = 0
  FFVar BC_hi = c;                     // c(1) = 0
  // Scalar closures (reference only scalars/constants -> no distributed coupling).
  FFVar SC1   = s1 - 42.;
  FFVar SC2   = s2 - 3. * s1;

  oc.add_equation( PDE,   { z }, { FFDom::ALL - FFDom::LB - FFDom::UB },
                   OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( BC_lo, { z }, { FFDom::LB },
                   OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );
  oc.add_equation( BC_hi, { z }, { FFDom::UB },
                   OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );
  oc.add_equation( SC1,   {}, {},
                   OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( SC2,   {}, {},
                   OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );

  oc.options.REDUCE.ORDER    = OCFESLV::Options::RED_FULL;
  oc.options.CLASSIFY.MODE        = OCFESLV::Options::CLASS_AUTO;
  oc.options.INTERFACE.IMPOSITION = OCFESLV::Options::IC_STRONG;
  oc.options.INTERFACE.TYPE  = OCFESLV::Options::IC_AUTO;
  oc.options.INTERFACE.SAT_SIGMA0      = 1.0;
  oc.options.DISPLAY_LEVEL   = 1;

  if( !oc.setup() ){ std::cerr << "ERROR: setup() failed\n"; return 2; }

  std::cout << "  nColloc states=" << oc.n_colloc_sta()
            << "  equations="       << oc.n_colloc_eqn()
            << "  square="          << ( oc.n_colloc_sta() == oc.n_colloc_eqn() ? "yes" : "NO" )
            << "\n";
  if( oc.n_colloc_sta() != oc.n_colloc_eqn() ){
    std::cerr << "ERROR: system is not square -- scalar row not assembled\n";
    return 3;
  }

  std::vector<double> varInit, inpInit;
  if( !oc.init( varInit, inpInit, nullptr ) ){ std::cerr << "ERROR: init() failed\n"; return 4; }

  std::vector<double> xv = varInit;
  OCFESLV::SolveReport rep = oc.solve( xv.data(), nullptr, nullptr );
  std::cout << "  [B] solve: converged=" << ( rep.converged ? "yes" : "no" )
            << " iters=" << rep.iterations
            << " final|r|=" << std::scientific << std::setprecision(3) << rep.final_residual << "\n";
  if( !rep.converged ){ std::cerr << "ERROR: solve did not converge\n"; return 5; }

  // ---- reads ----
  OCFESLV::t_Coord scalarpt;            // empty coord for a lumped state
  OCFESLV::t_Coord midpt;   midpt[z] = 0.5;

  double const s1_colloc  = oc.eval_colloc<double>( s1, scalarpt, xv.data(), nullptr, nullptr );
  double const s2_colloc  = oc.eval_colloc<double>( s2, scalarpt, xv.data(), nullptr, nullptr );
  double const c_mid      = oc.eval_colloc<double>( c,  midpt,    xv.data(), nullptr, nullptr );
  // buffer-free accessor should agree on the scalars too
  double const s1_sol     = oc.eval_solution( s1, scalarpt );
  double const s2_sol     = oc.eval_solution( s2, scalarpt );

  std::cout << "\n  [D] checks:\n";
  check( "scalar s1 == 42",          s1_colloc, 42.,  1e-9 );
  check( "scalar s2 == 3*s1 (126)",  s2_colloc, 126., 1e-9 );
  check( "distributed c(0.5)==0.25", c_mid,     0.25, 1e-8 );
  check( "eval_solution s1 == eval_colloc", s1_sol, s1_colloc, 1e-12 );
  check( "eval_solution s2 == eval_colloc", s2_sol, s2_colloc, 1e-12 );

  std::cout << "\n  RESULT: " << ( g_ok ? "ALL PASS" : "FAIL" ) << "\n";
  return g_ok ? 0 : 1;
}
