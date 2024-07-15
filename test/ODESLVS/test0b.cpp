#define USE_THREADS		// <- Whether to use parallel threads
#define MC__BASE_CVODES_CHECK

#include "odeslv_cvodes.hpp"

void solve_state( mc::ODESLV_CVODES<>& ivp, std::vector<double>& p, std::vector<double>& f )
{
  ivp.solve_state( p );
  f = ivp.val_function();
}

int main()
{
  mc::FFGraph IVP;  // DAG describing the problem

  const unsigned NS = 2;  // Time stages
  double t0 = 0., tf = 10.;       // Time span
  std::vector<double> T( NS+1 );  // Time stages
  for( unsigned int i=0; i<=NS; i++ ) T[i] = t0 + i * ( tf - t0 ) / NS; 

  const unsigned NP = 1;  // Number of parameters
  std::vector<mc::FFVar> P(NP);   // Parameters
  for( unsigned int i=0; i<NP; i++ ) P[i].set( &IVP );

  const unsigned NX = 2;  // Number of states
  std::vector<mc::FFVar> X(NX);   // States
  for( unsigned int i=0; i<NX; i++ ) X[i].set( &IVP );

  const unsigned NQ = 1;  // Number of state quadratures
  std::vector<mc::FFVar> Q(NQ);   // State quadratures
  for( unsigned i=0; i<NQ; i++ ) Q[i].set( &IVP );

  std::vector<mc::FFVar> RHS(NX); // Right-hand side function
  RHS[0] = P[0] * X[0] * ( 1. - X[1] );
  RHS[1] = P[0] * X[1] * ( X[0] - 1. );
/*
  std::vector<mc::FFVar> IC(NX);  // Initial value function
  IC[0] = 1.2;
  IC[1] = 1.1 + 0.01*(P[0]-3.);
*/
  std::vector<std::vector<mc::FFVar>> IC(NS);   // Initial value function
  IC[0].assign( { 1.2, 1.1 + 0.01*P[0] } );
  for( unsigned k=1; k<NS; k++ )
     IC[k].assign( { X[0], X[1] - 0.5 } );

  std::vector<mc::FFVar> QUAD(NQ);  // Quadrature function
  QUAD[0] = X[1];

  const unsigned NF = 2;  // Number of state functions
/*
  std::vector<mc::FFVar> FCT(NF);  // State functions
  FCT[0] = X[0] * X[1];
  FCT[1] = P[0] * pow( X[0], 2 );
*/
  std::vector<std::vector<mc::FFVar>> FCT(NS);  // State functions
  for( unsigned k=0; k<NS-1; k++ )
    FCT[k].assign( { 0., Q[0] } );
  FCT[NS-1].assign( { X[0] * X[1], Q[0] } );


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
  LV.set_state( X );
  LV.set_time( T );//t0, tf );
  LV.set_parameter( P );
  LV.set_differential( RHS );
  LV.set_initial( IC );
  LV.set_quadrature( QUAD, Q );
  LV.set_function( FCT );

  /////////////////////////////////////////////////////////////////////////
  // Make local copies
  
  unsigned const NTH = 20;
#if defined( USE_THREADS )
  std::vector<std::thread> vth( NTH );
#endif
  std::vector<mc::ODESLV_CVODES<>> vivp( NTH );
  std::vector<std::vector<double>> vf( NTH );
  std::vector<std::vector<double>> vp( NTH );
  double p0 = 3.0, dp = 0.1/(NTH-1.), pL = p0 - 0.1/2.;
  std::ofstream os;// = std::cout;

  for( unsigned ith=0; ith<NTH; ++ith ){
    vp[ith].assign( { pL+ith*dp } );  // Parameter values
    vivp[ith].set( LV );
    vivp[ith].options = LV.options;
    vivp[ith].setup();
#if defined( USE_THREADS )
    vth[ith] = std::thread( solve_state, std::ref(vivp[ith]), std::ref(vp[ith]), std::ref(vf[ith]) );
#else
    solve_state( vivp[ith], vp[ith], vf[ith] );
#endif
  }

  // Join the threads with the main thread
  for( unsigned ith=0; ith<NTH; ++ith ){
#if defined( USE_THREADS )
    vth[ith].join();
#endif
    std::cout << vp[ith][0];
    for( unsigned k=0; k<NF;++k )
      std::cout << "  " << vf[ith][k];
    std::cout << std::endl;
  }

  return 0;
}


