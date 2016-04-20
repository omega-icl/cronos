#define SAVE_RESULTS		// <- Whether to save bounds to file
#define USE_SUNDIALS		// <- whether to use SUNDIALS or GSL integrator

#ifndef USE_SUNDIALS
  #include "odeslvs_gsl.hpp"
  #include "interval.hpp"
  typedef mc::Interval I;
#else
  #include "odeslvs_sundials.hpp"
#endif

int main()
{
  mc::FFGraph IVP;  // DAG describing the problem

  double t0 = 0., tf = 6.;     // Time span
  const unsigned int NS = 1;  // Time stages
  double tk[NS+1]; tk[0] = t0;
  for( unsigned k=0; k<NS; k++ ) tk[k+1] = tk[k] + (tf-t0)/(double)NS;

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

  //mc::FFVar FCT[NF];  // State functions
  //FCT[0] = X[0] * X[1];
  ////FCT[1] = P[0] * pow( X[0], 2 );

  mc::FFVar FCT[NF*NS];  // State functions
  for( unsigned k=0; k<NF*NS; k++ ) FCT[k] = 0.;
  if( NS > 1 ) FCT[((NS-1)/NF)*NF+0] = X[0] + 0.1*P[0];
  FCT[(NS-1)*NF+0] = X[0] * X[1];
  FCT[(NS-1)*NF+1] = P[0] * pow( X[0], 2 );
  for( unsigned k=0; k<NS; k++ ) FCT[k*NF+1] += Q[0];

  /////////////////////////////////////////////////////////////////////////
  // Compute ODE solutions
#ifndef USE_SUNDIALS // GSL integrator
  mc::ODESLVS_GSL<I> LV;
  LV.options.INTMETH  = mc::ODESLV_GSL<I>::Options::MSBDF;
#else // SUNDIALS integrator
  mc::ODESLVS_SUNDIALS LV;
  LV.options.JACAPPROX = mc::ODESLV_SUNDIALS::Options::CV_DIAG;//CV_DENSE;//CV_DIAG;//
  LV.options.INTMETH   = mc::ODESLV_SUNDIALS::Options::MSBDF;//MSADAMS;//MSBDF;
  LV.options.NMAX      = 20000;
#endif

  LV.set_dag( &IVP );
  LV.set_state( NX, X );
  LV.set_parameter( NP, P );
  LV.set_differential( NX, RHS );
  LV.set_initial( NX, IC );
  //LV.set_initial( NS, NX, IC );
  LV.set_quadrature( NQ, QUAD, Q );
  //LV.set_function( NF, FCT );
  LV.set_function( NS, NF, FCT );

#if defined( SAVE_RESULTS )
  LV.options.RESRECORD = true;
#endif
  LV.options.DISPLAY   = 1;
  LV.options.ATOL      = LV.options.ATOLB      = LV.options.ATOLS  = 1e-9;
  LV.options.RTOL      = LV.options.RTOLB      = LV.options.RTOLS  = 1e-9;

  double p[NP] = { 2.95 };  // Parameter values
  double *xk[NS+1], f[NF], *xpk[NS+1], *lk[NS+1], fp[NF*NP];
  for( unsigned k=0; k<=NS; k++ ){
    xk[k]  = new double[NX];
    xpk[k] = new double[NX*NP];
    lk[k]  = new double[NX*NF];
  }
  std::ofstream direcSTA, direcSEN, direcADJ;

  std::cout << "\nCONTINUOUS-TIME INTEGRATION:\n\n";
  LV.states( NS, tk, p, xk, f );
#if defined( SAVE_RESULTS )
  direcSTA.open( "test1_STA.dat", std::ios_base::out );
  LV.record( direcSTA );
#endif

  std::cout << "\nCONTINUOUS-TIME INTEGRATION WITH FORWARD SENSITIVITY ANALYSIS:\n\n";
  LV.states_FSA( NS, tk, p, xk, f, xpk, fp );
#if defined( SAVE_RESULTS )
  direcSTA.open( "test1_STA.dat", std::ios_base::out );
  direcSEN.open( "test1_SEN.dat", std::ios_base::out );
  LV.record( direcSTA, direcSEN );
#endif

  std::cout << "\nCONTINUOUS-TIME INTEGRATION WITH ADJOINT SENSITIVITY ANALYSIS:\n\n";
  LV.states_ASA( NS, tk, p, xk, f, lk, fp );
#if defined( SAVE_RESULTS )
  direcSTA.open( "test1_STA.dat", std::ios_base::out );
  direcADJ.open( "test1_ADJ.dat", std::ios_base::out );
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


