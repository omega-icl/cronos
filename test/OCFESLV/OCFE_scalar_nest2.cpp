// ============================================================================
//  OCFE_scalar_nest2.cpp  --  reduction/derivative nesting coverage
//
//  Two elliptic 1-D fields (spatial, no evolution direction):
//      a(z) = z(1-z)  ( INT a = 1/6 , a'(z) = 1-2z )
//      b(y) = y(1-y)  ( INT b = 1/6 )
//  Nested reductions/derivatives on the product a*b (distributed over {z,y}):
//
//    p1 = OpEval( OpI(a*b,z), y, 0.5 )   = b(0.5)/6 = 1/24   [ integral-inside-point ]
//    p2 = OpI( OpI(a*b,z), y )           = (1/6)^2  = 1/36   [ double integral ]
//
//  Both cases produce a DISTRIBUTED integral LINK (free_dom-for-distributed).
//
//  Build:
//    g++ -std=c++17 ... -DOCFE_OCFESLV_HEADER='"ocfeslv_nest7.hpp"'
//        OCFE_scalar_nest2.cpp -o OCFE_scalar_nest2  <libs>
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
  std::cout << "  " << std::left << std::setw(34) << nm
            << std::right << std::scientific << std::setprecision(8)
            << " got=" << std::setw(15) << got << " want=" << std::setw(15) << want
            << std::setprecision(2) << "  err=" << std::setw(9) << err
            << "  " << ( pass ? "PASS" : "FAIL" ) << "\n";
}

int main()
{
  std::cout << "================================================================\n"
            << "  distributed integral LINKs: OpEval(OpI), OpI(OpI)\n"
            << "================================================================\n";

  FFGraph DAG;
  FFVar z  = DAG.add_var( "z" );
  FFVar y  = DAG.add_var( "y" );
  FFVar a  = DAG.add_var( "a(z)" );
  FFVar b  = DAG.add_var( "b(y)" );
  FFVar p1 = DAG.add_var( "p1" );
  FFVar p2 = DAG.add_var( "p2" );

  FFPartial  OpP;
  FFEval     OpEval;
  FFIntegral OpI;

  OCFESLV oc( &DAG );
  oc.add_domain( z, FFDom( 0., 1., 2, FFDom::CGL, 6 ) );
  oc.add_domain( y, FFDom( 0., 1., 2, FFDom::CGL, 6 ) );
  oc.add_state ( a,  { z } );
  oc.add_state ( b,  { y } );
  oc.add_state ( p1, {} );
  oc.add_state ( p2, {} );

  int const INT = FFDom::ALL - FFDom::LB - FFDom::UB;

  FFVar PDE_a = OpP( OpP( a, z ), z ) + 2.;
  FFVar PDE_b = OpP( OpP( b, y ), y ) + 2.;

  FFVar E1 = p1 - OpEval( OpI( a * b, z ), y, 0.5 );
  FFVar E2 = p2 - OpI( OpI( a * b, z ), y );

  oc.add_equation( PDE_a, { z }, { INT }, OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( a, { z }, { FFDom::LB }, OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );
  oc.add_equation( a, { z }, { FFDom::UB }, OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );
  oc.add_equation( PDE_b, { y }, { INT }, OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( b, { y }, { FFDom::LB }, OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );
  oc.add_equation( b, { y }, { FFDom::UB }, OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );
  oc.add_equation( E1, {}, {}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( E2, {}, {}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );

  oc.options.REDUCE.ORDER    = OCFESLV::Options::RED_FULL;
  oc.options.CLASSIFY.MODE        = OCFESLV::Options::CLASS_AUTO;
  oc.options.INTERFACE.IMPOSITION = OCFESLV::Options::IC_STRONG;
  oc.options.INTERFACE.TYPE  = OCFESLV::Options::IC_AUTO;
  oc.options.INTERFACE.SAT_SIGMA0      = 1.0;
  oc.options.DISPLAY_LEVEL   = 1;

  if( !oc.setup() ){ std::cerr << "ERROR: setup() failed\n"; return 2; }
  std::cout << "  setup: states=" << oc.n_colloc_sta()
            << " equations=" << oc.n_colloc_eqn()
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

  OCFESLV::t_Coord sp;
  double const v1 = oc.eval_colloc<double>( p1, sp, xv.data(), nullptr, nullptr );
  double const v2 = oc.eval_colloc<double>( p2, sp, xv.data(), nullptr, nullptr );

  std::cout << "\n  checks:\n";
  check( "p1 OpEval(OpI(a*b,z),y,0.5)==1/24",  v1, 1.0/24.0, 1e-8 );
  check( "p2 OpI(OpI(a*b,z),y)      ==1/36",  v2, 1.0/36.0, 1e-8 );

  std::cout << "\n  RESULT: " << ( g_ok ? "ALL PASS" : "FAIL" ) << "\n";
  return g_ok ? 0 : 1;
}
