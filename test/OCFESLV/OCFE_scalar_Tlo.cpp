// ============================================================================
//  OCFE_scalar_Tlo.cpp   --   section-4a acceptance test
//
//  A native lumped scalar (Tlo) closed by an FFIntegral of the DISTRIBUTED field
//  -- the exact shape of the MBC liquid-side energy balance:
//
//      distributed:  c'' + 2 = 0 , c(0)=c(1)=0        ->  c(z) = z(1-z)
//      scalar:       Tlo - INT_0^1 c dz = 0            ->  Tlo = 1/6 = 0.16666667
//
//  Exercises the section-4a machinery: _eqnEvalDom (the integrated direction z,
//  precomputed at setup) + the LOOP-2 scalar free-accumulate (sum of per-element
//  FFIntegral contributions into one scalar row).  Setup must be square (the row
//  assembles and the structural audit clears), the solve converges, and Tlo = 1/6.
//
//  Build (against the section-4a header), e.g.:
//    g++ -std=c++17 ... -DOCFE_OCFESLV_HEADER='"ocfeslv.hpp"'
//        OCFE_scalar_Tlo.cpp -o OCFE_scalar_Tlo  <libs>
// ============================================================================
#include <cmath>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
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
            << "  section-4a: scalar Tlo = INT_0^1 c dz  (c = z(1-z))\n"
            << "================================================================\n";

  FFGraph DAG;
  FFVar z   = DAG.add_var( "z" );
  FFVar c   = DAG.add_var( "c(z)" );
  FFVar Tlo = DAG.add_var( "Tlo" );      // lumped scalar

  FFPartial  OpP;
  FFIntegral OpI;

  OCFESLV oc( &DAG );
  oc.add_domain( z, FFDom( 0., 1., 2, FFDom::CGL, 6 ) );
  oc.add_state ( c,   { z } );           // distributed
  oc.add_state ( Tlo, {}    );           // scalar

  FFVar PDE   = OpP( OpP( c, z ), z ) + 2.;   // c'' = -2  ->  c = z(1-z)
  FFVar BC_lo = c;                            // c(0) = 0
  FFVar BC_hi = c;                            // c(1) = 0
  FFVar TL_EB = Tlo - OpI( c, z );            // scalar <- distributed field: Tlo = INT_0^1 c dz

  oc.add_equation( PDE,   { z }, { FFDom::ALL - FFDom::LB - FFDom::UB },
                   OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( BC_lo, { z }, { FFDom::LB },
                   OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );
  oc.add_equation( BC_hi, { z }, { FFDom::UB },
                   OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );
  oc.add_equation( TL_EB, {}, {},
                   OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );

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

  OCFESLV::t_Coord scalarpt;
  OCFESLV::t_Coord midpt;  midpt[z] = 0.5;
  double const Tlo_v = oc.eval_colloc<double>( Tlo, scalarpt, xv.data(), nullptr, nullptr );
  double const c_mid = oc.eval_colloc<double>( c,   midpt,    xv.data(), nullptr, nullptr );
  double const Tlo_s = oc.eval_solution( Tlo, scalarpt );

  std::cout << "\n  checks:\n";
  check( "Tlo == 1/6",              Tlo_v, 1.0/6.0, 1e-8 );
  check( "c(0.5) == 0.25",          c_mid, 0.25,    1e-8 );
  check( "eval_solution == colloc", Tlo_s, Tlo_v,   1e-12 );

  std::cout << "\n  RESULT: " << ( g_ok ? "ALL PASS" : "FAIL" ) << "\n";
  return g_ok ? 0 : 1;
}
