#define SAVE_RESULTS		// <- Whether to save bounds to file

#include "odeslvs_sundials.hpp"

int main()
{
  mc::FFGraph IVP;  // DAG describing the problem

  const unsigned int NS = 3;   // Time stages
  double tk[NS+1]; tk[0] = 0.;
  for( unsigned k=0; k<NS; k++ ) tk[k+1] = tk[k] + 1./(double)NS;
  mc::FFVar T; T.set( &IVP );  // Time

  const unsigned NP = 2+NS;    // Number of parameters
  const unsigned NX = 2;       // Number of states
  const unsigned NF = 3;       // Number of state functions

  mc::FFVar P[NP];             // Parameters
  for( unsigned int i=0; i<NP; i++ ) P[i].set( &IVP );
  mc::FFVar TF = P[0], X10 = P[1], *U = P+2;

  mc::FFVar X[NX];             // States
  for( unsigned int i=0; i<NX; i++ ) X[i].set( &IVP );

  mc::FFVar RHS[NX*NS];        // Right-hand side function
  for( unsigned k=0; k<NS; k++ ){
    RHS[NX*k+0] = TF * X[1];
    RHS[NX*k+1] = TF * ( U[k]*X[0] - 2.*X[1] ) * T;
  }

  mc::FFVar IC[NX];            // Initial value function
  IC[0] = 0.;
  IC[1] = X10;

  mc::FFVar FCT[NF*NS];        // State functions
  for( unsigned k=0; k<NF*NS; k++ ) FCT[k] = 0.;
  FCT[0] = TF;
  FCT[(NS-1)*NF+1] = X[0]-1.;
  FCT[(NS-1)*NF+2] = X[1];

  /////////////////////////////////////////////////////////////////////////
  // Compute ODE solutions

  mc::ODESLVS_SUNDIALS LV;

  LV.options.JACAPPROX = mc::BASE_SUNDIALS::Options::CV_DIAG;//CV_DENSE;
  LV.options.INTMETH   = mc::BASE_SUNDIALS::Options::MSBDF;//MSADAMS;
  LV.options.NMAX      = 0; //20000;
  LV.options.DISPLAY   = 1;
  LV.options.ATOL      = LV.options.ATOLB      = LV.options.ATOLS  = 1e-9;
  LV.options.RTOL      = LV.options.RTOLB      = LV.options.RTOLS  = 1e-9;
#if defined( SAVE_RESULTS )
  LV.options.RESRECORD = true;
#endif

  LV.set_dag( &IVP );
  LV.set_time( NS, tk, &T );
  LV.set_state( NX, X );
  LV.set_parameter( NP, P );
  LV.set_differential( NS, NX, RHS );
  LV.set_initial( NX, IC );
  LV.set_function( NS, NF, FCT );

  double p0[NP];
  p0[0] = 6.;
  p0[1] = 0.5;
  for( unsigned int is=0; is<NS; is++ )
    p0[2+is] = 0.5;

  double *xk[NS+1], f[NF], *xpk[NS+1], *lk[NS+1], fp[NF*NP];

  std::ofstream direcSTA, direcFSA[NP], direcASA[NF];
  char fname[50];

  std::cout << "\nCONTINUOUS-TIME INTEGRATION:\n\n";
  LV.states( NS, tk, p0, xk, f );
#if defined( SAVE_RESULTS )
  direcSTA.open( "test2_STA.dat", std::ios_base::out );
  LV.record( direcSTA );
#endif

  std::cout << "\nCONTINUOUS-TIME INTEGRATION WITH FORWARD SENSITIVITY ANALYSIS:\n\n";
  LV.states_FSA( NS, tk, p0, xk, f, xpk, fp );
#if defined( SAVE_RESULTS )
  direcSTA.open( "test2_STA.dat", std::ios_base::out );
  direcSEN.open( "test2_SEN.dat", std::ios_base::out );
  LV.record( direcSTA, direcSEN );
#endif

  std::cout << "\nCONTINUOUS-TIME INTEGRATION WITH ADJOINT SENSITIVITY ANALYSIS:\n\n";
  LV.states_ASA( NS, tk, p0, xk, f, lk, fp );
#if defined( SAVE_RESULTS )
  direcSTA.open( "test2_STA.dat", std::ios_base::out );
  direcADJ.open( "test2_ADJ.dat", std::ios_base::out );
  LV.record( direcSTA, direcADJ );
#endif

  // cleanup
  for( unsigned k=0; k<=NS; k++ ){
    delete[] xk[k];
    delete[] xpk[k];
    delete[] lk[k];
  }

  return 0;
}


