// ============================================================================
//  OCFE_scalar_fctopeval.cpp  --  OpEval point reduction in an OUTPUT
//
//  Task (3): add_output( OpEval(c,z,z0) ) as the point-REDUCTION spelling, on
//  the same footing as the equation OpEval.  The FFEval consumes z at z0; the
//  output scan extracts it into an aux state + LINK, so the output is
//  reduction-free (F = w) and _eval_fct reads the collocated w.  Cross-checked
//  against the existing fctrec.point spelling add_output(c,{z},{z0}) (output-node
//  PRESENTATION), which must give the same value.
//
//      c'' + 2 = 0 , c(0)=0.2 , c(1)=0.5   ->  c = -z^2 + 1.3 z + 0.2
//                    ( c(0)=0.2 , c(0.5)=0.6 )
//
//      out0 = OpEval(c,z,0.5)            -> 0.6   [OpEval, interior]
//      out1 = c  @ {z},{0.5}            -> 0.6   [fctrec.point]
//      out2 = OpEval(c,z,0.0)            -> 0.2   [OpEval, boundary]
//      out3 = c  @ {z},{0.0}            -> 0.2   [fctrec.point]
//
//  Build:
//    g++ -std=c++17 ... -DOCFE_OCFESLV_HEADER='"ocfeslv_fct.hpp"'
//        OCFE_scalar_fctopeval.cpp -o OCFE_scalar_fctopeval  <libs>
// ============================================================================
#include <cmath>
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

int main()
{
  std::cout << "================================================================\n"
            << "  OpEval in an OUTPUT: add_output(OpEval(c,z,z0))\n"
            << "================================================================\n";

  FFGraph DAG;
  FFVar z = DAG.add_var( "z" );
  FFVar c = DAG.add_var( "c(z)" );

  FFPartial OpP;
  FFEval    OpEval;

  OCFESLV oc( &DAG );
  oc.add_domain( z, FFDom( 0., 1., 2, FFDom::CGL, 6 ) );
  oc.add_state ( c, { z } );

  FFVar PDE   = OpP( OpP( c, z ), z ) + 2.;   // c'' = -2
  FFVar BC_lo = c - 0.2;                        // c(0) = 0.2
  FFVar BC_hi = c - 0.5;                        // c(1) = 0.5

  oc.add_equation( PDE,   { z }, { FFDom::ALL - FFDom::LB - FFDom::UB },
                   OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( BC_lo, { z }, { FFDom::LB }, OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );
  oc.add_equation( BC_hi, { z }, { FFDom::UB }, OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );

  // out0/out2: OpEval point reduction in the output expression
  oc.add_output( OpEval( c, z, 0.5 ) );
  oc.add_output( c, { z }, { 0.5 } );          // out1: fctrec.point (presentation) cross-check
  oc.add_output( OpEval( c, z, 0.0 ) );  // out2
  oc.add_output( c, { z }, { 0.0 } );          // out3: fctrec.point cross-check

  oc.options.REDUCE.ORDER    = OCFESLV::Options::RED_FULL;
  oc.options.CLASSIFY.MODE        = OCFESLV::Options::CLASS_AUTO;
  oc.options.INTERFACE.IMPOSITION = OCFESLV::Options::IC_STRONG;
  oc.options.INTERFACE.TYPE  = OCFESLV::Options::IC_AUTO;
  oc.options.INTERFACE.SAT_SIGMA0      = 1.0;
  oc.options.DISPLAY_LEVEL   = 1;

  if( !oc.setup() ){ std::cerr << "ERROR: setup() failed\n"; return 2; }
  std::cout << "  setup: states=" << oc.n_colloc_sta()
            << " equations=" << oc.n_colloc_eqn()
            << " functions=" << oc.n_colloc_fct()
            << " square=" << ( oc.n_colloc_sta() == oc.n_colloc_eqn() ? "yes" : "NO" ) << "\n";
  if( oc.n_colloc_sta() != oc.n_colloc_eqn() ){ std::cerr << "ERROR: not square\n"; return 3; }

  std::vector<double> varInit, inpInit;
  if( !oc.init( varInit, inpInit, nullptr ) ){ std::cerr << "ERROR: init() failed\n"; return 4; }

  std::vector<double> xv = varInit;
  OCFESLV::SolveReport rep = oc.solve( xv.data(), nullptr, nullptr );
  std::cout << "  solve: converged=" << ( rep.converged ? "yes" : "no" )
            << " iters=" << rep.iterations
            << " final|r|=" << std::scientific << std::setprecision(3) << rep.final_residual << "\n";
  if( !rep.converged ){ std::cerr << "ERROR: solve did not converge\n"; return 5; }

  std::vector<double> res( oc.n_colloc_eqn(), 0. );
  std::vector<double> fct( oc.n_colloc_fct(), 0. );
  if( !oc.eval( res.data(), fct.data(), xv.data(), nullptr, nullptr ) ){
    std::cerr << "ERROR: eval() failed\n"; return 6;
  }
  if( fct.size() < 4 ){ std::cerr << "ERROR: expected >=4 outputs\n"; return 7; }

  std::cout << "\n  checks:\n";
  check( "out0 OpEval(c,z,0.5) == 0.6",   fct[0], 0.6, 1e-8 );
  check( "out1 fctrec.point   == 0.6",    fct[1], 0.6, 1e-8 );
  check( "OpEval == fctrec.point @ 0.5",  fct[0], fct[1], 1e-10 );
  check( "out2 OpEval(c,z,0.0) == 0.2",   fct[2], 0.2, 1e-8 );
  check( "out3 fctrec.point   == 0.2",    fct[3], 0.2, 1e-8 );
  check( "OpEval == fctrec.point @ 0.0",  fct[2], fct[3], 1e-10 );

  std::cout << "\n  RESULT: " << ( g_ok ? "ALL PASS" : "FAIL" ) << "\n";
  return g_ok ? 0 : 1;
}
