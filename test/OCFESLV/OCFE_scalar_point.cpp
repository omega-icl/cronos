// ============================================================================
//  OCFE_scalar_point.cpp  --  point reduction via OpEval (FFEval)
//
//  A lumped scalar closed by a DISTRIBUTED field read at a POINT, written as an
//  OpEval operator in the residual.  (The former equation-record point spec has
//  been retired in favour of OpEval.)  The equation is a plain domain-less
//  scalar; the FFEval node consumes z at the target coordinate during the DAG
//  eval, so OCres.dom() comes back empty and the row-write is trivial.
//
//      distributed:  c'' + 2 = 0 , c(0)=0.2 , c(1)=0.5   ->  c = -z^2 + 1.3 z + 0.2
//                    ( c(0)=0.2 , c(0.5)=0.6 , c(1)=0.5 )
//      p_mid : p_mid - OpEval(c,z,0.5) = 0   [interior]  ->  0.6
//      p_lo  : p_lo  - OpEval(c,z,0.0) = 0   [boundary]  ->  0.2
//
//  Build:
//    g++ -std=c++17 ... -DOCFE_OCFESLV_HEADER='"ocfeslv.hpp"'
//        OCFE_scalar_point.cpp -o OCFE_scalar_point  <libs>
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
  std::cout << "  " << std::left << std::setw(26) << nm
            << std::right << std::scientific << std::setprecision(8)
            << " got=" << std::setw(15) << got << " want=" << std::setw(15) << want
            << std::setprecision(2) << "  err=" << std::setw(9) << err
            << "  " << ( pass ? "PASS" : "FAIL" ) << "\n";
}

int main()
{
  std::cout << "================================================================\n"
            << "  point (OpEval): p = OpEval(c,z,z0)   (c = -z^2 + 1.3 z + 0.2)\n"
            << "================================================================\n";

  FFGraph DAG;
  FFVar z     = DAG.add_var( "z" );
  FFVar c     = DAG.add_var( "c(z)" );
  FFVar p_mid = DAG.add_var( "p_mid" );      // interior point z=0.5
  FFVar p_lo  = DAG.add_var( "p_lo"  );      // boundary point z=0

  FFPartial OpP;
  FFEval    OpEval;

  OCFESLV oc( &DAG );
  oc.add_domain( z, FFDom( 0., 1., 2, FFDom::CGL, 6 ) );
  oc.add_state ( c,     { z } );
  oc.add_state ( p_mid, {}    );
  oc.add_state ( p_lo,  {}    );

  FFVar PDE   = OpP( OpP( c, z ), z ) + 2.;   // c'' = -2
  FFVar BC_lo = c - 0.2;                        // c(0) = 0.2
  FFVar BC_hi = c - 0.5;                        // c(1) = 0.5
  FFVar PT_M  = p_mid - OpEval( c, z, 0.5 );   // p_mid = c(0.5)
  FFVar PT_L  = p_lo  - OpEval( c, z, 0.0 );   // p_lo  = c(0.0)

  oc.add_equation( PDE,   { z }, { FFDom::ALL - FFDom::LB - FFDom::UB },
                   OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( BC_lo, { z }, { FFDom::LB },
                   OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );
  oc.add_equation( BC_hi, { z }, { FFDom::UB },
                   OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );
  oc.add_equation( PT_M, {}, {}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( PT_L, {}, {}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );

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
  OCFESLV::t_Coord midpt;  midpt[z] = 0.5;
  double const pm = oc.eval_colloc<double>( p_mid, sp,    xv.data(), nullptr, nullptr );
  double const pl = oc.eval_colloc<double>( p_lo,  sp,    xv.data(), nullptr, nullptr );
  double const cm = oc.eval_colloc<double>( c,     midpt, xv.data(), nullptr, nullptr );

  std::cout << "\n  checks:\n";
  check( "p_mid == c(0.5) == 0.6", pm, 0.6, 1e-8 );
  check( "p_lo  == c(0.0) == 0.2", pl, 0.2, 1e-8 );
  check( "c(0.5) == 0.6",          cm, 0.6, 1e-8 );

  std::cout << "\n  RESULT: " << ( g_ok ? "ALL PASS" : "FAIL" ) << "\n";
  return g_ok ? 0 : 1;
}
