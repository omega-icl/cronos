#define USE_THREADS		// <- Whether to use parallel threads
#define MC__BASE_CVODES_CHECK

#include "odeslv_cvodes.hpp"

void states( mc::ODESLV_CVODES<>* ivp, double*p, double*f )
{
  ivp->states( p, nullptr, f );
}

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
  // Define ODE problem

  mc::ODESLV_CVODES LV;

  LV.options.INTMETH   = mc::BASE_CVODES::Options::MSBDF;//MSADAMS;
  LV.options.NLINSOL   = mc::BASE_CVODES::Options::NEWTON;//FIXEDPOINT;
  LV.options.LINSOL    = mc::BASE_CVODES::Options::DENSE;//DENSEDQ;//DIAG;
  LV.options.NMAX      = 2000;
  LV.options.DISPLAY   = 0;
  LV.options.ATOL      = 1e-9;
  LV.options.RTOL      = 1e-9;

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

  /////////////////////////////////////////////////////////////////////////
  // Make local copies
  
  unsigned const NTH = 20;
#if defined( USE_THREADS )
  std::vector<std::thread> vth( NTH );
#endif
  std::vector<mc::ODESLV_CVODES<>> vivp( NTH );
  std::vector<double[NF]> vf( NTH );
  std::vector<double[2][NX]> vxk( NTH );
  std::vector<double> vp( NTH );
  double p0 = 3.0, dp = 0.1/(NTH-1.), pL = p0 - 0.1/2.;
  std::ofstream os;// = std::cout;

  for( unsigned ith=0; ith<NTH; ++ith ){
    vp[ith] = pL+ith*dp;  // Parameter values
    vivp[ith].set( LV );
    vivp[ith].options = LV.options;
    vivp[ith].setup();
#if defined( USE_THREADS )
    vth[ith] = std::thread( states, &vivp[ith], &vp[ith], vf[ith] );
    //vth[ith] = std::thread( &mc::ODESLV_CVODES<>::states, &vivp[ith], &vp[ith], &vxk[ith], vf[ith], os );
#else
    vivp[ith].states( &vp[ith], nullptr, vf[ith] );
#endif
  }

  // Join the threads with the main thread
  for( unsigned ith=0; ith<NTH; ++ith ){
#if defined( USE_THREADS )
    vth[ith].join();
#endif
    std::cout << vp[ith];
    for( unsigned k=0; k<NF;++k )
      std::cout << "  " << vf[ith][k];
    std::cout << std::endl;
  }

  return 0;
}


