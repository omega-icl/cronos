// test2_ff.cpp -- test2, with the model defined through FFModel.
// =============================================================================================================
// The interesting part of test2 is the CONTROL.  The original declares NP = 2 + NS parameters -- TF, X10 and one
// level U[k] per stage -- and writes NS stage-wise right-hand sides, RHS[k] using U[k].  Here the control is ONE
// time-varying input, declared piecewise constant on the same grid:
//
//     IVP.add_input( U, { t }, mc::FFDom::LGR, 1 );     // one value per element
//
// so there is ONE pair of equations, not NS.  The extraction then does what the original did by hand: a
// piecewise-constant input becomes ONE PARAMETER PER ELEMENT, appended to the parameter vector in element order,
// and the right-hand sides become stage-wise with that element's level substituted.  The parameter vector is
// therefore [ TF, X10, U[0], U[1], U[2] ] -- the same NP = 2 + NS, in the same order, so the p0 below is the
// original's unchanged.
//
// TF and X10 are time-invariant inputs; the evolution variable t appears explicitly in the second equation, which
// the model carries as itself.
//
// A profile that is NOT piecewise constant is REFUSED with a message rather than silently left dangling: an
// integrator would need its interpolant evaluated at t inside the right-hand side, which is not implemented yet.
// =============================================================================================================
#define SAVE_RESULTS

#include "odeslvs_cvodes.hpp"

int main()
{
  mc::FFGraph IVPDAG;
  mc::ODESLVS_CVODES IVP( &IVPDAG );

  size_t const NS = 3;                       // stages -> the evolution domain's elements
  double const t0 = 0., tf = 1.;

  mc::FFVar t   = IVPDAG.add_var( "t" );
  mc::FFVar X0  = IVPDAG.add_var( "x0(t)" ), X1 = IVPDAG.add_var( "x1(t)" );
  mc::FFVar TF  = IVPDAG.add_var( "TF" ), X10 = IVPDAG.add_var( "x10" );
  mc::FFVar U   = IVPDAG.add_var( "u(t)" );
  mc::FFPartial OpP;  mc::FFEval OpEval;

  IVP.add_domain( t, mc::FFDom( t0, tf, NS, mc::FFDom::LGR, 4 ) );
  IVP.add_state( X0, {t} );  IVP.add_state( X1, {t} );
  IVP.add_input( TF );  IVP.add_input( X10 );                  // time-invariant
  IVP.add_input( U, { t }, mc::FFDom::LGR, 1 );                // piecewise constant: one level per stage
  IVP.set_evolution_domain( t );
  IVP.update_ref( X0, 0. );  IVP.update_ref( X1, 0.5 );

  int const T_INT = mc::FFDom::ALL - mc::FFDom::LB;
  auto const io = mc::FFModel::EqnOptions( mc::FFModel::EqnRole::INTERIOR, 0 );
  auto const ii = mc::FFModel::EqnOptions( mc::FFModel::EqnRole::INITIAL,  0 );

  IVP.add_equation( OpP( X0, t ) - TF * X1,                       {t}, {T_INT},          io );
  IVP.add_equation( OpP( X1, t ) - TF * ( U * X0 - 2. * X1 ) * t, {t}, {T_INT},          io );
  IVP.add_equation( X0 - 0.,                                      {t}, {mc::FFDom::LB},  ii );
  IVP.add_equation( X1 - X10,                                     {t}, {mc::FFDom::LB},  ii );

  IVP.add_output( OpEval( X0 - 1., t, tf ) );                  // test2's FCT[1] and FCT[2] at the end; FCT[0] was
  IVP.add_output( OpEval( X1,      t, tf ) );                  // the parameter TF itself, which needs no output

  IVP.options.INTMETH   = mc::BASE_CVODES::Options::MSBDF;
  IVP.options.NLINSOL   = mc::BASE_CVODES::Options::NEWTON;
  IVP.options.LINSOL    = mc::BASE_CVODES::Options::DENSE;//SPARSE;
  IVP.options.NMAX      = 0;
  IVP.options.DISPLAY   = 1;
  IVP.options.ATOL      = IVP.options.ATOLB = IVP.options.ATOLS = 1e-9;
  IVP.options.RTOL      = IVP.options.RTOLB = IVP.options.RTOLS = 1e-9;
#if defined( SAVE_RESULTS )
  IVP.options.RESRECORD = true;
#endif

  if( !IVP.setup() ){
    std::cerr << "setup failed: " << IVP.extract_error() << std::endl;
    return 1;
  }

  std::vector<double> p0( { 6., 0.5, 0.5, 0.5, 0.5 } );        // TF, X10, then one level per stage

  std::ofstream direcSTA;

  std::cout << "\nCONTINUOUS-TIME INTEGRATION:\n\n";
  IVP.solve( p0 );
#if defined( SAVE_RESULTS )
  direcSTA.open( "test2_STA.dat", std::ios_base::out );
  IVP.record( direcSTA );
#endif

  std::cout << "\nCONTINUOUS-TIME INTEGRATION WITH FORWARD SENSITIVITY ANALYSIS:\n\n";
  IVP.solve_fsens( p0 );

  std::cout << "\nCONTINUOUS-TIME INTEGRATION WITH ADJOINT SENSITIVITY ANALYSIS:\n\n";
  IVP.solve_asens( p0 );

  return 0;
}
