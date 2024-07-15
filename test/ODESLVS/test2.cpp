#define SAVE_RESULTS		// <- Whether to save bounds to file

#include "odeslvs_cvodes.hpp"

int main()
{
  mc::FFGraph IVPDAG;  // DAG describing the problem

  const unsigned NS = 3;  // Time stages
  double t0 = 0., tf = 1.;       // Time span
  std::vector<double> TS( NS+1 );  // Time stages
  for( unsigned int i=0; i<=NS; i++ ) TS[i] = t0 + i * ( tf - t0 ) / NS; 
  mc::FFVar T; T.set( &IVPDAG );  // Time

  const unsigned NP = 2+NS;  // Number of parameters
  std::vector<mc::FFVar> P(NP);   // Parameters
  for( unsigned int i=0; i<NP; i++ ) P[i].set( &IVPDAG );
  mc::FFVar TF = P[0], X10 = P[1], *U = &P[2];

  const unsigned NX = 2;  // Number of states
  std::vector<mc::FFVar> X(NX);   // States
  for( unsigned int i=0; i<NX; i++ ) X[i].set( &IVPDAG );

  std::vector<mc::FFVar> IC(NX);  // Initial value function
  IC[0] = 0.;
  IC[1] = X10;

  std::vector<std::vector<mc::FFVar>> RHS(NS); // Right-hand side function
  for( unsigned k=0; k<NS; k++ )
    RHS[k].assign( { TF * X[1],
                     TF * ( U[k]*X[0] - 2.*X[1] ) * T } );

  const unsigned NF = 3;  // Number of state functions
  std::vector<mc::FFVar> FCT(NF);  // State functions
  FCT[0] = TF;
  FCT[1] = X[0]-1.;
  FCT[2] = X[1];


  /////////////////////////////////////////////////////////////////////////
  // Compute ODE solutions

  mc::ODESLVS_CVODES IVP;

  IVP.options.LINSOL    = mc::BASE_CVODES::Options::DIAG;//DENSE;
  IVP.options.INTMETH   = mc::BASE_CVODES::Options::MSBDF;//MSADAMS;
  IVP.options.NMAX      = 0; //20000;
  IVP.options.DISPLAY   = 1;
  IVP.options.ATOL      = IVP.options.ATOLB      = IVP.options.ATOLS  = 1e-9;
  IVP.options.RTOL      = IVP.options.RTOLB      = IVP.options.RTOLS  = 1e-9;
#if defined( SAVE_RESULTS )
  IVP.options.RESRECORD = true;
#endif

  IVP.set_dag( &IVPDAG );
  IVP.set_time( TS, &T );
  IVP.set_state( X );
  IVP.set_parameter( P );
  IVP.set_differential( RHS );
  IVP.set_initial( IC );
  IVP.set_function( FCT );
  IVP.setup();
  
  std::vector<double> p0(NP);
  p0[0] = 6.;
  p0[1] = 0.5;
  for( unsigned int k=0; k<NS; k++ )
    p0[2+k] = 0.5;

  std::ofstream direcSTA, direcFSA[NP], direcASA[NF];
  char fname[50];

  std::cout << "\nCONTINUOUS-TIME INTEGRATION:\n\n";
  IVP.solve_state( p0 );
#if defined( SAVE_RESULTS )
  direcSTA.open( "test2_STA.dat", std::ios_base::out );
  IVP.record( direcSTA );
#endif

  std::cout << "\nCONTINUOUS-TIME INTEGRATION WITH FORWARD SENSITIVITY ANALYSIS:\n\n";
  IVP.solve_sensitivity( p0 );
#if defined( SAVE_RESULTS )
  direcSTA.open( "test2_STA.dat", std::ios_base::out );
  for( unsigned i=0; i<NP; ++i ){
    sprintf( fname, "test2_FSA%d.dat",i );  
    direcFSA[i].open( fname, std::ios_base::out );
  }
  IVP.record( direcSTA, direcFSA );
#endif

  std::cout << "\nCONTINUOUS-TIME INTEGRATION WITH ADJOINT SENSITIVITY ANALYSIS:\n\n";
  IVP.solve_adjoint( p0 );
#if defined( SAVE_RESULTS )
  direcSTA.open( "test2_STA.dat", std::ios_base::out );
  for( unsigned i=0; i<NF; ++i ){
    sprintf( fname, "test2_ASA%d.dat",i );  
    direcASA[i].open( fname, std::ios_base::out );
  }
  IVP.record( direcSTA, direcASA );
#endif

  return 0;
}

