// ============================================================================
//  OCFE_AE.cpp  --  OCFESLV applied to a COMPLETELY ALGEBRAIC system
//
//  The degenerate case: NO domain, NO partial derivatives, NO integrals.  Every
//  state is a lumped (empty-domain) scalar, every equation is a pointwise scalar
//  closure.  The environment therefore holds ZERO collocation domains, and setup()
//  must assemble a plain square nonlinear system F(x)=0 that solve() drives to a
//  root with its equilibrated LM/Newton hybrid.
//
//  WHAT THIS EXERCISES
//  -------------------
//    * add_state(x,{})            -- lumped scalar states, no domain at all
//    * add_equation(F,{},{})      -- the empty-domain ("LOOP-2") scalar-row idiom
//                                    with the environment carrying no domains
//    * a square 5x5 assembly, a converged solve, and scalar reads via an EMPTY
//      t_Coord through both eval_colloc() and eval_solution().
//    * four flavours of algebraic row:
//        - linear                         F1 : x1 - 3
//        - scalar-scalar linear coupling  F2 : x2 - 2 x1
//        - transcendental nonlinear       F3 : exp(x3) - 2
//        - a genuinely coupled nonlinear  F4 : x4 x5 - 6
//          2x2 block                      F5 : x4 - x5 - 1
//
//  UNIQUE SOLUTION (reachable from the seeded guesses)
//  ---------------------------------------------------
//        x1 = 3
//        x2 = 6
//        x3 = ln 2 ~ 0.6931471805599453
//        x4 = 3 ,  x5 = 2      (positive branch of the coupled block)
//
//  The coupled block's OTHER root is (x4,x5)=(-2,-3); the origin makes its
//  Jacobian [[x5,x4],[1,-1]] singular, so x4,x5 are seeded on the positive branch.
//  The linear and transcendental rows converge from the default-zero guess.
//
//  Build (against the live header):
//    g++ -std=c++17 <suite flags> -DOCFE_OCFESLV_HEADER='"ocfeslv.hpp"' \
//        OCFE_AE.cpp -o OCFE_AE <libs>
//  Run:
//    ./OCFE_AE        # prints a per-row PASS/FAIL table; exit 0 iff all pass
// ============================================================================

#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
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
  double const err  = std::fabs( got - want );
  bool   const pass = ( err < tol );
  g_ok &= pass;
  std::cout << "  " << std::left << std::setw(34) << nm
            << std::right << std::scientific << std::setprecision(8)
            << " got=" << std::setw(15) << got
            << " want=" << std::setw(15) << want
            << std::setprecision(2) << "  err=" << std::setw(9) << err
            << "  " << ( pass ? "PASS" : "FAIL" ) << "\n";
}

int main()
{
  std::cout << "================================================================\n"
            << "  OCFE_AE : OCFESLV on a COMPLETELY ALGEBRAIC system\n"
            << "  no domain, no partial derivatives, no integrals\n"
            << "  (all states are lumped / empty-domain scalars)\n"
            << "================================================================\n";

  FFGraph DAG;
  FFVar x1 = DAG.add_var( "x1" );
  FFVar x2 = DAG.add_var( "x2" );
  FFVar x3 = DAG.add_var( "x3" );
  FFVar x4 = DAG.add_var( "x4" );
  FFVar x5 = DAG.add_var( "x5" );

  OCFESLV oc( &DAG );
  // No add_domain(): the environment carries ZERO collocation domains.
  oc.add_state( x1, {} );
  oc.add_state( x2, {} );
  oc.add_state( x3, {} );
  oc.add_state( x4, {} );
  oc.add_state( x5, {} );

  // Initial guesses.  State references default to zero; the coupled block is
  // singular at the origin, so seed x4,x5 on the positive branch.  x1,x2,x3 are
  // left at the default zero to confirm those rows converge from it.
  oc.update_ref( x4, 1.0 );
  oc.update_ref( x5, 1.0 );

  FFVar F1 = x1 - 3.;
  FFVar F2 = x2 - 2. * x1;
  FFVar F3 = exp( x3 ) - 2.;
  FFVar F4 = x4 * x5 - 6.;
  FFVar F5 = x4 - x5 - 1.;

  OCFESLV::EqnOptions alg( OCFESLV::EqnRole::INTERIOR, 0 );
  oc.add_equation( F1, {}, {}, alg );
  oc.add_equation( F2, {}, {}, alg );
  oc.add_equation( F3, {}, {}, alg );
  oc.add_equation( F4, {}, {}, alg );
  oc.add_equation( F5, {}, {}, alg );

  oc.options.REDUCE.ORDER   = OCFESLV::Options::RED_NONE;    // no derivatives to reduce
  oc.options.CLASSIFY.MODE       = OCFESLV::Options::CLASS_NONE;  // no PDE to classify
  oc.options.SOLVE.MAX_ITER = 60;
  oc.options.SOLVE.RES_TOL  = 1.0e-9;
  oc.options.DISPLAY_LEVEL  = 1;

  if( !oc.setup() ){
    std::cerr << "ERROR: setup() failed: "
              << OCFESLV::setup_status_str( oc.setup_status() ) << "\n";
    return 2;
  }

  size_t const nsta = oc.n_colloc_sta();
  size_t const neqn = oc.n_colloc_eqn();
  std::cout << "  nColloc states=" << nsta << "  equations=" << neqn
            << "  square=" << ( nsta == neqn ? "yes" : "NO" ) << "\n";
  if( nsta != neqn || nsta != 5 ){
    std::cerr << "ERROR: expected a square 5x5 algebraic system (got "
              << nsta << "x" << neqn << ")\n";
    return 3;
  }

  std::vector<double> varInit, inpInit;
  if( !oc.init( varInit, inpInit, nullptr ) ){
    std::cerr << "ERROR: init() failed\n";
    return 4;
  }
  double* const pinp = inpInit.empty() ? nullptr : inpInit.data();

  std::vector<double> xv = varInit;
  OCFESLV::SolveReport rep = oc.solve( xv.data(), pinp, nullptr );
  std::cout << "  solve: converged=" << ( rep.converged ? "yes" : "no" )
            << " iters=" << rep.iterations
            << " |r0|=" << std::scientific << std::setprecision(3) << rep.initial_residual
            << " final|r|=" << rep.final_residual << "\n";
  if( !rep.converged ){
    std::cerr << "ERROR: solve did not converge\n";
    return 5;
  }

  OCFESLV::t_Coord pt;   // empty coordinate: lumped scalar read
  double const v1 = oc.eval_colloc<double>( x1, pt, xv.data(), pinp, nullptr );
  double const v2 = oc.eval_colloc<double>( x2, pt, xv.data(), pinp, nullptr );
  double const v3 = oc.eval_colloc<double>( x3, pt, xv.data(), pinp, nullptr );
  double const v4 = oc.eval_colloc<double>( x4, pt, xv.data(), pinp, nullptr );
  double const v5 = oc.eval_colloc<double>( x5, pt, xv.data(), pinp, nullptr );

  std::cout << "\n  checks:\n";
  check( "x1 == 3",                          v1, 3.,             1e-9 );
  check( "x2 == 2 x1 (6)",                   v2, 6.,             1e-9 );
  check( "x3 == ln 2",                       v3, std::log( 2. ), 1e-9 );
  check( "x4 == 3 (coupled block)",          v4, 3.,             1e-9 );
  check( "x5 == 2 (coupled block)",          v5, 2.,             1e-9 );
  check( "eval_solution x1 == eval_colloc",  oc.eval_solution( x1, pt ), v1, 1e-12 );
  check( "eval_solution x4 == eval_colloc",  oc.eval_solution( x4, pt ), v4, 1e-12 );

  std::cout << "\n================================================================\n"
            << "  OCFE_AE (completely algebraic, zero-domain): "
            << ( g_ok ? "ALL PASS" : "SOME FAILED" ) << "\n"
            << "================================================================\n";
  return g_ok ? 0 : 1;
}
