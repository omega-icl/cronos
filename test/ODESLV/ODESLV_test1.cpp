// test1_ff.cpp -- the original test1, with the model defined through FFModel instead of BASE_DE.
// ==============================================================================================
// The solver classes now derive from FFModel, so the IVP object IS the model: the description goes in through
// add_domain / add_state / add_input / add_equation / add_output, and setup() does both halves.
//
// WHAT MAPS DIRECTLY
//   set_time( T )        -> add_domain( t, FFDom( t0, tf, NS, ... ) ): the stages ARE the domain's elements
//   set_state / _parameter -> add_state( x, {t} ) / add_input( p )
//   set_differential     -> add_equation( dx/dt - f, {t}, {(a,b]}, INTERIOR )
//   set_initial          -> add_equation( x - x0, {t}, {LB}, INITIAL )
//   set_quadrature       -> add_output( OpI( integrand, t ) ): a deferred INTEGRAL record, and its value over the
//                           whole horizon IS test1's terminal quadrature function
//   set_function         -> add_output( OpEval( expr, t, tau ) ): a function of the states at a time point, which
//                           need NOT be a stage boundary (see the one at t = 3.7 below)
//
// WHAT DOES NOT MAP, and it is a real limit rather than an oversight
//   test1 registers FCT[k] = { {1, Q[0]} } at EVERY stage -- the RUNNING quadrature.  Writing that here as
//   OpEval( Q, t, tau ) is accepted when declared and REFUSED at setup: "nested evolution-direction reduction --
//   an inner reduction's captured (post-solve) value feeds an outer evolution".  The guard is right for a
//   collocation solve, where a captured value is post-solve and using it inside the same solve is circular, but
//   it is coarser than needed here, where the value is only REPORTED.  Until that is relaxed, the running
//   quadrature is not a model output; the terminal one is.
//
// ALSO DROPPED, by design: test1's stage-wise initial values IC[k] (state discontinuities), which the refactored
// path does not model.
// ==============================================================================================
#define SAVE_RESULTS

#include "odeslvs_cvodes.hpp"

int main()
{
  mc::FFGraph IVPDAG;
  mc::ODESLVS_CVODES IVP( &IVPDAG );          // the solver IS the model: it takes the DAG

  double const t0 = 0., tf = 10.;
  size_t const NS = 4;                        // stages -> the evolution domain's elements

  mc::FFVar t = IVPDAG.add_var( "t" );
  mc::FFVar X0 = IVPDAG.add_var( "x0(t)" ), X1 = IVPDAG.add_var( "x1(t)" );
  mc::FFVar P0 = IVPDAG.add_var( "p0" ),    P1 = IVPDAG.add_var( "p1" );
  mc::FFPartial OpP;  mc::FFIntegral OpI;  mc::FFEval OpE;

  IVP.add_domain( t, mc::FFDom( t0, tf, NS, mc::FFDom::LGR, 4 ) );
  IVP.add_state( X0, {t} );  IVP.add_state( X1, {t} );
  IVP.add_input( P0 );       IVP.add_input( P1 );
  IVP.set_evolution_domain( t );
  IVP.update_ref( X0, 1.2 ); IVP.update_ref( X1, 1.1 );

  int const T_INT = mc::FFDom::ALL - mc::FFDom::LB;
  auto const io = mc::FFModel::EqnOptions( mc::FFModel::EqnRole::INTERIOR );//, 0 );
  auto const ii = mc::FFModel::EqnOptions( mc::FFModel::EqnRole::INITIAL  );//,  0 );

  IVP.add_equation( OpP( X0, t ) - P0 * X0 * ( 1. - X1 ), {t}, {T_INT},          io );
  IVP.add_equation( OpP( X1, t ) - P0 * X1 * ( X0 - 1. ), {t}, {T_INT},          io );
  IVP.add_equation( X0 - 1.2,                             {t}, {mc::FFDom::LB},  ii );
  IVP.add_equation( X1 - ( 1.1 + 0.01 * P1 ),             {t}, {mc::FFDom::LB},  ii );

  IVP.add_output( OpI( X1, t ) );                   // the quadrature, and its terminal value
  IVP.add_output( OpE( X0 * X1, t, tf ) );          // test1's terminal state function
  //IVP.add_output( OpE( X0, t, 3.7 ) );              // a function OFF the stage partition (the addition)

  IVP.options.INTMETH   = mc::BASE_CVODES::Options::MSBDF;
  IVP.options.NLINSOL   = mc::BASE_CVODES::Options::NEWTON;
  IVP.options.LINSOL    = mc::BASE_CVODES::Options::SPARSE;//DENSE;//SPARSE;
  IVP.options.FSACORR   = mc::BASE_CVODES::Options::STAGGERED;
  IVP.options.NMAX      = 2000;
  IVP.options.DISPLAY   = 1;
  IVP.options.ATOL      = IVP.options.ATOLB = IVP.options.ATOLS = 1e-9;
  IVP.options.RTOL      = IVP.options.RTOLB = IVP.options.RTOLS = 1e-9;
  IVP.options.QERR      = IVP.options.QERRS = 1;
  IVP.options.ASACHKPT  = 2000;
#if defined( SAVE_RESULTS )
  IVP.options.RESRECORD = 100;
#endif

  if( !IVP.setup() ){                                   // one call: the model, then the local copy
    std::cerr << "setup failed: " << IVP.extract_error() << std::endl;
    return 1;
  }

  std::vector<double> p( { 2.96, 3. } );

  std::cout << "\nCONTINUOUS-TIME INTEGRATION:\n\n";
  IVP.solve( p );

  std::cout << "\nCONTINUOUS-TIME INTEGRATION WITH FORWARD SENSITIVITY ANALYSIS:\n\n";
  IVP.solve_fsens( p );

  std::cout << "\nCONTINUOUS-TIME INTEGRATION WITH ADJOINT SENSITIVITY ANALYSIS:\n\n";
  IVP.solve_asens( p );

  return 0;
}
