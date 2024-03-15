#undef SAVE_RESULTS		// <- Whether to save bounds to file

#include "odeslvs_cvodes.hpp"

int main()
{
  mc::FFGraph IVPDAG;  // DAG describing the problem

  const unsigned NS = 10;  // Time stages
  double tk[NS+1]; tk[0] = 0.;
  for( unsigned k=0; k<NS; k++ ) tk[k+1] = tk[k] + 1./(double)NS;

  const unsigned NP = NS; // Number of parameters
  const unsigned NX = 1;  // Number of states
  const unsigned NQ = 1;  // Number of state quadratures
  const unsigned NF = 2;  // Number of state functions

  mc::FFVar P[NP];  // Parameters
  for( unsigned int i=0; i<NP; i++ ) P[i].set( &IVPDAG );

  mc::FFVar X[NX];  // States
  for( unsigned int i=0; i<NX; i++ ) X[i].set( &IVPDAG );

  mc::FFVar Q[NQ];  // State quadratures
  for( unsigned i=0; i<NQ; i++ ) Q[i].set( &IVPDAG );

  mc::FFVar RHS[NX*NS];  // Right-hand side function
  for( unsigned k=0; k<NS; k++ )
    RHS[k] = P[k] - X[0];

  mc::FFVar IC[NX];   // Initial value function
  IC[0] = 1e0;

  mc::FFVar QUAD[NQ*NS];  // Quadrature function
  for( unsigned k=0; k<NS; k++ )
    QUAD[k] = 0.5 * mc::sqr( P[k] );

  mc::FFVar FCT[NF];  // State functions
  FCT[0] = Q[0];
  FCT[1] = X[0];

  /////////////////////////////////////////////////////////////////////////
  // Compute ODE solutions

  mc::ODESLVS_CVODES IVP;

  IVP.options.INTMETH   = mc::BASE_CVODES::Options::MSBDF;//MSADAMS;
  IVP.options.NLINSOL   = mc::BASE_CVODES::Options::NEWTON;//FIXEDPOINT;
  IVP.options.LINSOL    = mc::BASE_CVODES::Options::DENSE;//DIAG;
  IVP.options.FSACORR   = mc::BASE_CVODES::Options::STAGGERED;//STAGGERED1;//SIMULTANEOUS;
  IVP.options.NMAX      = 0;
  IVP.options.DISPLAY   = 1;
  IVP.options.ATOL      = IVP.options.ATOLB     = IVP.options.ATOLS  = 1e-8;
  IVP.options.RTOL      = IVP.options.RTOLB     = IVP.options.RTOLS  = 1e-7;
  IVP.options.QERR      = IVP.options.QERRS     = 1;
  IVP.options.ASACHKPT  = 2000;
#if defined( SAVE_RESULTS )
  IVP.options.RESRECORD = 20;
#endif
  
  IVP.set_dag( &IVPDAG );
  IVP.set_time( NS, tk );
  IVP.set_state( NX, X );
  IVP.set_parameter( NP, P );
  IVP.set_differential( NS, NX, RHS );
  IVP.set_initial( NX, IC );
  IVP.set_quadrature( NS, NQ, QUAD, Q );
  IVP.set_function( NF, FCT );
  IVP.setup();
  
  double p0[NP];
  for( unsigned int is=0; is<NS; is++ )
    p0[is] = -3e-1;

  std::ofstream direcSTA, direcFSA[NP], direcASA[NF];
  char fname[50];

  std::cout << "\nCONTINUOUS-TIME INTEGRATION:\n\n";
  IVP.states( p0 );
#if defined( SAVE_RESULTS )
  direcSTA.open( "test4_STA.dat", std::ios_base::out );
  IVP.record( direcSTA );
#endif

  std::cout << "\nCONTINUOUS-TIME INTEGRATION WITH FORWARD SENSITIVITY ANALYSIS:\n\n";
  IVP.states_FSA( p0 );
#if defined( SAVE_RESULTS )
  direcSTA.open( "test4_STA.dat", std::ios_base::out );
  for( unsigned i=0; i<NP; ++i ){
    sprintf( fname, "test4_FSA%d.dat",i );  
    direcFSA[i].open( fname, std::ios_base::out );
  }
  IVP.record( direcSTA, direcFSA );
#endif

  std::cout << "\nCONTINUOUS-TIME INTEGRATION WITH ADJOINT SENSITIVITY ANALYSIS:\n\n";
  IVP.states_ASA( p0 );
#if defined( SAVE_RESULTS )
  direcSTA.open( "test4_STA.dat", std::ios_base::out );
  for( unsigned i=0; i<NF; ++i ){
    sprintf( fname, "test4_ASA%d.dat",i );  
    direcASA[i].open( fname, std::ios_base::out );
  }
  IVP.record( direcSTA, direcASA );
#endif

  return 0;
}


