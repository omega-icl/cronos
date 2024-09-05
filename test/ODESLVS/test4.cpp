#define USE_PROFIL
#include <fstream>
#include <iomanip>

#include "ffode.hpp"

////////////////////////////////////////////////////////////////////////
int main()
////////////////////////////////////////////////////////////////////////
{
  /////////////////////////////////////////////////////////////////////////
  // Define IVP-ODE

  mc::FFGraph IVPDAG;  // DAG describing the problem

  const unsigned NS = 5;  // Time stages
  double t0 = 0., tf = 1.;       // Time span
  std::vector<double> T( NS+1 );  // Time stages
  for( unsigned int i=0; i<=NS; i++ ) T[i] = t0 + i * ( tf - t0 ) / NS; 

  const unsigned NP = NS;  // Number of parameters
  std::vector<mc::FFVar> P(NP);   // Parameters
  for( unsigned int i=0; i<NP; i++ ) P[i].set( &IVPDAG );

  const unsigned NX = 1;  // Number of states
  std::vector<mc::FFVar> X(NX);   // States
  for( unsigned int i=0; i<NX; i++ ) X[i].set( &IVPDAG );

  const unsigned NQ = 1;  // Number of state quadratures
  std::vector<mc::FFVar> Q(NQ);   // State quadratures
  for( unsigned i=0; i<NQ; i++ ) Q[i].set( &IVPDAG );

  std::vector<std::vector<mc::FFVar>> RHS(NS); // Right-hand side function
  for( unsigned k=0; k<NS; k++ )
    RHS[k].assign( { P[k] - X[0] } );

  std::vector<mc::FFVar> IC(NX);  // Initial value function
  IC[0] = 1e0;

  std::vector<std::vector<mc::FFVar>> QUAD(NS); // Quadrature function
  for( unsigned k=0; k<NS; k++ )
    QUAD[k].assign( { 0.5 * mc::sqr( P[k] ) } );

  const unsigned NF = 2;  // Number of state functions
  std::vector<std::vector<mc::FFVar>> FCT(NS);  // State functions
  for( unsigned k=0; k<NS-1; k++ )
    FCT[k].assign( { 0., Q[0] } );
  FCT[NS-1].assign( { X[0], Q[0] } );


  mc::ODESLVS_CVODES IVP;
  
  IVP.options.INTMETH   = mc::BASE_CVODES::Options::MSBDF;//MSADAMS;//
  IVP.options.NLINSOL   = mc::BASE_CVODES::Options::FIXEDPOINT;//NEWTON;//
  IVP.options.LINSOL    = mc::BASE_CVODES::Options::DIAG;//DENSE;//
  IVP.options.FSACORR   = mc::BASE_CVODES::Options::STAGGERED;//STAGGERED1;//SIMULTANEOUS;
  IVP.options.NMAX      = 2000;
  IVP.options.DISPLAY   = 1;
  IVP.options.ATOL      = IVP.options.ATOLB     = IVP.options.ATOLS  = 1e-9;
  IVP.options.RTOL      = IVP.options.RTOLB     = IVP.options.RTOLS  = 1e-8;
  IVP.options.FSAERR    = IVP.options.QERR      = IVP.options.QERRS     = 1;
  IVP.options.ASACHKPT  = 2000;

  IVP.set_dag( &IVPDAG );
  IVP.set_time( T );
  IVP.set_state( X );
  IVP.set_parameter( P );
  IVP.set_differential( RHS );
  IVP.set_initial( IC );
  IVP.set_quadrature( QUAD, Q );
  IVP.set_function( FCT );

  IVP.setup();
  IVP.solve_state( std::vector<double>( NP, -1e0 ) );

  /////////////////////////////////////////////////////////////////////////
  // Define DAG

  mc::FFGraph DAG;

  mc::FFVar PP[NP];  // Parameters
  for( unsigned int i=0; i<NP; i++ ) PP[i].set( &DAG );

  mc::FFODE OpODE;
  mc::FFVar F[NF];
  for( unsigned int j=0; j<NF; j++ ) F[j] = OpODE( j, NP, PP, &IVP );
  std::cout << DAG;

  std::vector<double> dP( NP, -1e0 ), dF( NF );
  DAG.eval( NF, F, dF.data(), NP, PP, dP.data() ); 

  return 0;
}
