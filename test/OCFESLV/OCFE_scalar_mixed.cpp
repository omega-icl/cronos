// ============================================================================
//  OCFE_scalar_mixed.cpp  --  mixed & nonlinear reductions via reduce_order
//                             extraction (OpI + OpEval -> aux states + LINKs)
//
//  Exercises what the uniform reduce_order extraction unlocks, which neither the
//  direct free-accumulate nor the direct OpEval path can handle in one residual:
//
//      c'' + 2 = 0 , c(0)=0.2 , c(1)=0.5   ->  c = -z^2 + 1.3 z + 0.2  (c(0.5)=0.6)
//      d'' = 0     , d(0)=0   , d(1)=1     ->  d = z                   (INT d = 0.5)
//
//    MIXED (point + integral over the same z, in one equation):
//      p - OpEval(c,z,0.5) - OpI(d,z) = 0   ->  p = 0.6 + 0.5 = 1.1
//    NONLINEAR in a reduction:
//      q - OpI(d,z)^2 = 0                    ->  q = 0.5^2 = 0.25
//
//  Extraction turns each reduction into its own linear LINK
//  ( w1 - OpEval(c,z,0.5)=0 , w2 - OpI(d,z)=0 ) and the parents become
//  reduction-free ( p - w1 - w2 = 0 , q - w2^2 = 0 ).
//
//  Build:
//    g++ -std=c++17 ... -DOCFE_OCFESLV_HEADER='"ocfeslv_ro.hpp"'
//        OCFE_scalar_mixed.cpp -o OCFE_scalar_mixed  <libs>
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
  std::cout << "  " << std::left << std::setw(30) << nm
            << std::right << std::scientific << std::setprecision(8)
            << " got=" << std::setw(15) << got << " want=" << std::setw(15) << want
            << std::setprecision(2) << "  err=" << std::setw(9) << err
            << "  " << ( pass ? "PASS" : "FAIL" ) << "\n";
}

int main()
{
  std::cout << "================================================================\n"
            << "  mixed/nonlinear reductions (reduce_order extraction)\n"
            << "================================================================\n";

  FFGraph DAG;
  FFVar z = DAG.add_var( "z" );
  FFVar c = DAG.add_var( "c(z)" );
  FFVar d = DAG.add_var( "d(z)" );
  FFVar p = DAG.add_var( "p" );   // p = c(0.5) + INT d      (mixed)
  FFVar q = DAG.add_var( "q" );   // q = (INT d)^2           (nonlinear)

  FFPartial  OpP;
  FFEval     OpEval;
  FFIntegral OpI;

  OCFESLV oc( &DAG );
  oc.add_domain( z, FFDom( 0., 1., 2, FFDom::CGL, 6 ) );
  oc.add_state ( c, { z } );
  oc.add_state ( d, { z } );
  oc.add_state ( p, {}    );
  oc.add_state ( q, {}    );

  FFVar PDE_c = OpP( OpP( c, z ), z ) + 2.;   // c'' = -2
  FFVar BCc_lo = c - 0.2, BCc_hi = c - 0.5;    // c(0)=0.2, c(1)=0.5
  FFVar PDE_d = OpP( OpP( d, z ), z );         // d'' = 0
  FFVar BCd_lo = d,       BCd_hi = d - 1.;     // d(0)=0,   d(1)=1

  FFVar MIX = p - OpEval( c, z, 0.5 ) - OpI( d, z );   // p = c(0.5) + INT d
  FFVar NL  = q - OpI( d, z ) * OpI( d, z );           // q = (INT d)^2

  oc.add_equation( PDE_c,  { z }, { FFDom::ALL - FFDom::LB - FFDom::UB },
                   OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( BCc_lo, { z }, { FFDom::LB }, OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );
  oc.add_equation( BCc_hi, { z }, { FFDom::UB }, OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );
  oc.add_equation( PDE_d,  { z }, { FFDom::ALL - FFDom::LB - FFDom::UB },
                   OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( BCd_lo, { z }, { FFDom::LB }, OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );
  oc.add_equation( BCd_hi, { z }, { FFDom::UB }, OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );
  oc.add_equation( MIX, {}, {}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( NL,  {}, {}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );

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
  double const pv = oc.eval_colloc<double>( p, sp, xv.data(), nullptr, nullptr );
  double const qv = oc.eval_colloc<double>( q, sp, xv.data(), nullptr, nullptr );

  std::cout << "\n  checks:\n";
  check( "p == c(0.5) + INT d == 1.1", pv, 1.1,  1e-8 );
  check( "q == (INT d)^2     == 0.25", qv, 0.25, 1e-8 );

  std::cout << "\n  RESULT: " << ( g_ok ? "ALL PASS" : "FAIL" ) << "\n";
  return g_ok ? 0 : 1;
}
