#define SAVE_RESULTS		// <- Whether to save bounds to file

#include "odeslvs_sundials.hpp"

int main()
{
  mc::FFGraph IVP;  // DAG describing the problem

  double t0 = 0., tf = 10.;     // Time span
  const unsigned NS = 1;  // Time stages

  const unsigned NP = 1;  // Number of parameters
  const unsigned NX = 2;  // Number of states
  const unsigned NQ = 1;  // Number of state quadratures
  const unsigned NF = 2;  // Number of state functions

  mc::FFVar P[NP];  // Parameters
  for( unsigned int i=0; i<NP; i++ ) P[i].set( &IVP );

  mc::FFVar X[NX];  // States
  for( unsigned int i=0; i<NX; i++ ) X[i].set( &IVP );

  mc::FFVar Q[NQ];  // State quadratures
  for( unsigned i=0; i<NQ; i++ ) Q[i].set( &IVP );

  mc::FFVar RHS[NX];  // Right-hand side function
  RHS[0] = P[0] * X[0] * ( 1. - X[1] );
  RHS[1] = P[0] * X[1] * ( X[0] - 1. );

  mc::FFVar IC[NX];   // Initial value function
  IC[0] = 1.2;
  IC[1] = 1.1 + 0.01*(P[0]-3.);
/*
  mc::FFVar IC[NX*NS];   // Initial value function
  for( unsigned k=0; k<NX*NS; k++ ) IC[k] = X[k%NX];
  IC[0] = 1.2;
  IC[1] = 1.1;// + 0.01*P[0];
  //if( NS > 1 ) IC[(NS/2)*NX+1] = X[1] - 0.5;
*/
  mc::FFVar QUAD[NQ];  // Quadrature function
  QUAD[0] = X[1];
/*
  mc::FFVar FCT[NF];  // State functions
  FCT[0] = X[0] * X[1];
  //FCT[1] = P[0] * pow( X[0], 2 );
*/
  mc::FFVar FCT[NF*NS];  // State functions
  for( unsigned k=0; k<NF*NS; k++ ) FCT[k] = 0.;
  if( NS > 1 ) FCT[((NS-1)/NF)*NF+0] = X[0] + 0.1*P[0];
  FCT[(NS-1)*NF+0] = X[0] * X[1];
  FCT[(NS-1)*NF+1] = P[0] * pow( X[0], 2 );
  for( unsigned k=0; k<NS; k++ ) FCT[k*NF+1] += Q[0];

  /////////////////////////////////////////////////////////////////////////
  // Compute ODE solutions

  mc::ODESLVS_SUNDIALS LV;

  LV.options.JACAPPROX = mc::BASE_SUNDIALS::Options::CV_DIAG;//CV_DENSE;//CV_DIAG;//
  LV.options.INTMETH   = mc::BASE_SUNDIALS::Options::MSBDF;//MSADAMS;//MSBDF;
  LV.options.NMAX      = 2000;
  LV.options.DISPLAY   = 1;
  LV.options.ATOL      = LV.options.ATOLB      = LV.options.ATOLS  = 1e-9;
  LV.options.RTOL      = LV.options.RTOLB      = LV.options.RTOLS  = 1e-9;
  LV.options.ASACHKPT  = 200;
#if defined( SAVE_RESULTS )
  LV.options.RESRECORD = 100;
#endif

  LV.set_dag( &IVP );
  LV.set_state( NX, X );
  LV.set_time( t0, tf );
  LV.set_parameter( NP, P );
  LV.set_differential( NX, RHS );
  LV.set_initial( NX, IC );
  //LV.set_initial( NS, NX, IC );
  LV.set_quadrature( NQ, QUAD, Q );
  //LV.set_function( NF, FCT );
  LV.set_function( NS, NF, FCT );

  double p[NP] = { 2.95 };  // Parameter values

  std::ofstream direcSTA, direcFSA[NP], direcASA[NF];
  char fname[50];

  std::cout << "\nCONTINUOUS-TIME INTEGRATION:\n\n";
  LV.states( p ); //, xk, f );
#if defined( SAVE_RESULTS )
  direcSTA.open( "test1_STA.dat", std::ios_base::out );
  LV.record( direcSTA );
#endif

  std::cout << "\nCONTINUOUS-TIME INTEGRATION WITH FORWARD SENSITIVITY ANALYSIS:\n\n";
  LV.states_FSA( p ); //, xk, f, xpk, fp );
#if defined( SAVE_RESULTS )
  direcSTA.open( "test1_STA.dat", std::ios_base::out );
  for( unsigned i=0; i<NP; ++i ){
    sprintf( fname, "test1_FSA%d.dat",i );  
    direcFSA[i].open( fname, std::ios_base::out );
  }
  LV.record( direcSTA, direcFSA );
#endif

  std::cout << "\nCONTINUOUS-TIME INTEGRATION WITH ADJOINT SENSITIVITY ANALYSIS:\n\n";
  //for( unsigned i=0; i<1000; i++ )
  LV.states_ASA( p ); //, xk, f, lk, fp );
#if defined( SAVE_RESULTS )
  direcSTA.open( "test1_STA.dat", std::ios_base::out );
  for( unsigned i=0; i<NF; ++i ){
    sprintf( fname, "test1_ASA%d.dat",i );  
    direcASA[i].open( fname, std::ios_base::out );
  }
  LV.record( direcSTA, direcASA );
#endif

  return 0;
}


