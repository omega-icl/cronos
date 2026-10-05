// ============================================================================
//  OCFE_scalar_nest3.cpp  --  reductions over COMPOUND derivative operands
//
//  Tests OpI/OpEval whose operand is a constant*derivative or a NONLINEAR
//  function of a derivative -- all resolved by the tainted-operand
//  materialization (materialise the whole operand as a state, so the reduction
//  sees a state and FAD never reaches a deriv()).
//
//      a(z) = z(1-z) , a'(z) = 1-2z ,  INT_0^1 (a')^2 dz = INT (1-2z)^2 dz = 1/3
//      u(z) = z      , u'(z) = 1    ,  INT_0^1 u' dz     = u(1)-u(0) = 1
//
//    q1 = OpI( OpP(a,z)^2, z )        = 1/3        [ nonlinear f(derivative) ]
//    q2 = OpI( 2.0*OpP(a,z)^2, z )    = 2/3        [ const * nonlinear ]
//    q3 = OpI( 3.0*OpP(u,z), z )      = 3          [ const * derivative ]
//
//  Build:
//    g++ -std=c++17 ... -DOCFE_OCFESLV_HEADER='"ocfeslv.hpp"'
//    w4 = d/dz( u*du/dz ) = 1     [ OpP chain: OpP(u*OpP(u,z),z) ]
//        OCFE_scalar_nest3.cpp -o OCFE_scalar_nest3  <libs>
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
            << "  reductions over compound derivative operands\n"
            << "================================================================\n";

  FFGraph DAG;
  FFVar z  = DAG.add_var( "z" );
  FFVar a  = DAG.add_var( "a(z)" );
  FFVar u  = DAG.add_var( "u(z)" );
  FFVar q1 = DAG.add_var( "q1" );
  FFVar q2 = DAG.add_var( "q2" );
  FFVar q3 = DAG.add_var( "q3" );
  FFVar w4 = DAG.add_var( "w4(z)" );

  FFPartial  OpP;
  FFIntegral OpI;

  OCFESLV oc( &DAG );
  oc.add_domain( z, FFDom( 0., 1., 2, FFDom::CGL, 6 ) );
  oc.add_state ( a,  { z } );
  oc.add_state ( u,  { z } );
  oc.add_state ( q1, {} );
  oc.add_state ( q2, {} );
  oc.add_state ( q3, {} );
  oc.add_state ( w4, { z } );  // w4 = d/dz( u * du/dz )  [ FFPartial chain: OpP(u*OpP(u,z),z) ]

  int const INT = FFDom::ALL - FFDom::LB - FFDom::UB;

  FFVar PDE_a = OpP( OpP( a, z ), z ) + 2.;   // a'' = -2 -> a = z(1-z)
  FFVar PDE_u = OpP( OpP( u, z ), z );         // u'' = 0  -> u = z (with u(0)=0,u(1)=1)

  FFVar Da = OpP( a, z );
  FFVar Du = OpP( u, z );
  FFVar E1 = q1 - OpI( Da * Da, z );
  FFVar E2 = q2 - OpI( 2.0 * ( Da * Da ), z );
  FFVar E3 = q3 - OpI( 3.0 * Du, z );
  FFVar E4 = w4 - OpP( u * OpP( u, z ), z );   // OpP OF a product containing a derivative

  oc.add_equation( PDE_a, { z }, { INT }, OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( a, { z }, { FFDom::LB }, OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );
  oc.add_equation( a, { z }, { FFDom::UB }, OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );
  oc.add_equation( PDE_u, { z }, { INT }, OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( u,       { z }, { FFDom::LB }, OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );  // u(0)=0
  oc.add_equation( u - 1.0, { z }, { FFDom::UB }, OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );  // u(1)=1
  oc.add_equation( E1, {}, {}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( E2, {}, {}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( E3, {}, {}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( E4, { z }, { FFDom::ALL }, OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );

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
  double const v1 = oc.eval_colloc<double>( q1, sp, xv.data(), nullptr, nullptr );
  double const v2 = oc.eval_colloc<double>( q2, sp, xv.data(), nullptr, nullptr );
  double const v3 = oc.eval_colloc<double>( q3, sp, xv.data(), nullptr, nullptr );
  OCFESLV::t_Coord zc; zc[z] = 0.5;
  double const v4 = oc.eval_colloc<double>( w4, zc, xv.data(), nullptr, nullptr );

  std::cout << "\n  checks:\n";
  check( "q1 INT (a')^2 dz    == 1/3", v1, 1.0/3.0, 1e-7 );
  check( "q2 INT 2(a')^2 dz   == 2/3", v2, 2.0/3.0, 1e-7 );
  check( "q3 INT 3 u' dz      == 3",   v3, 3.0,     1e-7 );
  check( "w4 d(u u')/dz @0.5  == 1",   v4, 1.0,     1e-7 );

  std::cout << "\n  RESULT: " << ( g_ok ? "ALL PASS" : "FAIL" ) << "\n";
  return g_ok ? 0 : 1;
}
