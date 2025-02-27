#define SAVE_RESULTS		// <- Whether to save bounds to file

#include "odeslvs_cvodes.hpp"

int main()
{
  mc::FFGraph IVPDAG;  // DAG describing the problem

  const unsigned NS = 4;  // Time stages
  double t0 = 0., tf = 10.;       // Time span
  std::vector<double> T( NS+1 );  // Time stages
  for( unsigned int i=0; i<=NS; i++ ) T[i] = t0 + i * ( tf - t0 ) / NS; 

  const unsigned NP = 2;  // Number of parameters
  std::vector<mc::FFVar> P(NP);   // Parameters
  for( unsigned int i=0; i<NP; i++ ) P[i].set( &IVPDAG );

  const unsigned NX = 2;  // Number of states
  std::vector<mc::FFVar> X(NX);   // States
  for( unsigned int i=0; i<NX; i++ ) X[i].set( &IVPDAG );

  const unsigned NQ = 1;  // Number of state quadratures
  std::vector<mc::FFVar> Q(NQ);   // State quadratures
  for( unsigned i=0; i<NQ; i++ ) Q[i].set( &IVPDAG );

  std::vector<mc::FFVar> RHS(NX); // Right-hand side function
  RHS[0] = P[0] * X[0] * ( 1. - X[1] );
  RHS[1] = P[0] * X[1] * ( X[0] - 1. );

//  std::vector<mc::FFVar> IC(NX);  // Initial value function
//  IC[0] = 1.2;
//  IC[1] = 1.1 + 0.01*P[1];

  std::vector<std::vector<mc::FFVar>> IC(NS);   // Initial value function
  IC[0].assign( { 1.2, 1.1 + 0.01*P[1] } );
  for( unsigned k=1; k<NS; k++ )
     IC[k].assign( { X[0], X[1] - 0.02*P[1] } );

  std::vector<mc::FFVar> QUAD(NQ);  // Quadrature function
  QUAD[0] = X[1];

  const unsigned NF = 2;  // Number of state functions

//  std::vector<mc::FFVar> FCT(NF);  // State functions
//  FCT[0] = X[0] * X[1];
//  FCT[1] = Q[0]; //P[0] * pow( X[0], 2 );
//  //FCT[0] = Q[0]; //P[0] * pow( X[0], 2 );

  std::vector<std::map<size_t,mc::FFVar>> FCT(NS+1);  // State functions
  for( unsigned k=1; k<NS; k++ )
    FCT[k] = { { 1, Q[0] } };
  FCT[NS] = { { 0, X[0] * X[1] }, { 1, Q[0] } };


  /////////////////////////////////////////////////////////////////////////
  // Compute ODE solutions

  mc::ODESLVS_CVODES IVP;

  IVP.options.INTMETH   = mc::BASE_CVODES::Options::MSBDF;//MSADAMS;
  IVP.options.NLINSOL   = mc::BASE_CVODES::Options::NEWTON;//FIXEDPOINT;
  IVP.options.LINSOL    = mc::BASE_CVODES::Options::SPARSE;//DENSE;//DIAG;
  IVP.options.FSACORR   = mc::BASE_CVODES::Options::STAGGERED;//STAGGERED1;//SIMULTANEOUS;
  IVP.options.NMAX      = 2000;
  IVP.options.DISPLAY   = 1;
  IVP.options.ATOL      = IVP.options.ATOLB     = IVP.options.ATOLS  = 1e-9;
  IVP.options.RTOL      = IVP.options.RTOLB     = IVP.options.RTOLS  = 1e-9;
  IVP.options.QERR      = IVP.options.QERRS     = 1;
  IVP.options.ASACHKPT  = 2000;
#if defined( SAVE_RESULTS )
  IVP.options.RESRECORD = 100;
#endif

  IVP.set_dag( &IVPDAG );
  IVP.set_time( T );
  IVP.set_state( X );
  IVP.set_parameter( P );
  IVP.set_differential( RHS );
  IVP.set_initial( IC );
  IVP.set_quadrature( QUAD, Q );
  IVP.set_function( FCT );
  IVP.setup();
  
  std::vector<double> p( { 2.96, 3. } );  // Parameter values

#if defined( SAVE_RESULTS )
  std::ofstream direcSTA, direcFSA[NP], direcASA[NF];
  char fname[50];
#endif

  std::cout << "\nCONTINUOUS-TIME INTEGRATION:\n\n";
  IVP.solve_state( p );
#if defined( SAVE_RESULTS )
  direcSTA.open( "test1_STA.dat", std::ios_base::out );
  IVP.record( direcSTA );
#endif

  std::cout << "\nCONTINUOUS-TIME INTEGRATION WITH FORWARD SENSITIVITY ANALYSIS:\n\n";
  IVP.solve_sensitivity( p );
#if defined( SAVE_RESULTS )
  direcSTA.open( "test1_STA.dat", std::ios_base::out );
  for( unsigned i=0; i<NP; ++i ){
    sprintf( fname, "test1_FSA%d.dat",i );  
    direcFSA[i].open( fname, std::ios_base::out );
  }
  IVP.record( direcSTA, direcFSA );
#endif

  std::cout << "\nCONTINUOUS-TIME INTEGRATION WITH ADJOINT SENSITIVITY ANALYSIS:\n\n";
  //for( unsigned i=0; i<1000; i++ )
  IVP.solve_adjoint( p );
#if defined( SAVE_RESULTS )
  direcSTA.open( "test1_STA.dat", std::ios_base::out );
  for( unsigned i=0; i<NF; ++i ){
    sprintf( fname, "test1_ASA%d.dat",i );  
    direcASA[i].open( fname, std::ios_base::out );
  }
  IVP.record( direcSTA, direcASA );
#endif

  return 0;
}


