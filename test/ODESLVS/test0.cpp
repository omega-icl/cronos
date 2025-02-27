#define SAVE_RESULTS		// <- Whether to save bounds to file
#define MC__BASE_CVODES_CHECK

#include "odeslvs_cvodes.hpp"

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
  IC[1] = 1.1 + 0.01*P[0];
*/
  std::vector<std::vector<mc::FFVar>> IC(NS);   // Initial value function
  IC[0].assign( { 1.2, 1.1 + 0.01*P[0] } );
  for( unsigned k=1; k<NS; k++ )
     IC[k].assign( { X[0], X[1] - 0.5 } );

  std::vector<mc::FFVar> QUAD(NQ);  // Quadrature function
  QUAD[0] = X[1];

  const unsigned NF = 2;  // Number of state functions
/*
  std::map<size_t,mc::FFVar> FCT{ { 0, X[0] * X[1] }, { 1, Q[0] } };

  std::vector<mc::FFVar> FCT{ X[0] * X[1], Q[0] };

  std::vector<std::vector<mc::FFVar>> FCT(NS+1);  // State functions
  for( unsigned k=0; k<NS; k++ )
    FCT[k].assign( { 0., Q[0] } );
  FCT[NS].assign( { X[0] * X[1], Q[0] } );
*/
  std::vector<std::map<size_t,mc::FFVar>> FCT(NS+1); // State functions
  //FCT[0].insert( { 0, P[0] } );
  for( unsigned k=1; k<NS; k++ )
    FCT[k].insert( { 1, Q[0] } );
  FCT[NS].insert( { 0, X[0] * X[1] } );
  FCT[NS].insert(  { 1, Q[0] } );


  /////////////////////////////////////////////////////////////////////////
  // Compute ODE solutions

  mc::ODESLVS_CVODES LV;

  LV.options.INTMETH   = mc::BASE_CVODES::Options::MSBDF;//MSADAMS;
  LV.options.NLINSOL   = mc::BASE_CVODES::Options::NEWTON;//FIXEDPOINT;
  LV.options.LINSOL    = mc::BASE_CVODES::Options::SPARSE;//DENSE;//DENSEDQ;//DIAG;
  LV.options.NMAX      = 2000;
  LV.options.DISPLAY   = 1;
  LV.options.ATOL      = 1e-9;
  LV.options.RTOL      = 1e-9;
#if defined( SAVE_RESULTS )
  LV.options.RESRECORD = 100;
#endif

  LV.set_dag( &IVP );
  LV.set_state( X );
  LV.set_time( T );//t0, tf );
  LV.set_parameter( P );
  LV.set_differential( RHS );
  LV.set_initial( IC );
  LV.set_quadrature( QUAD, Q );
  LV.set_function( FCT );
  LV.setup();

  std::cout << "\nCONTINUOUS-TIME STATE INTEGRATION:\n\n";
  std::vector<double> p( { 2.95 } );  // Parameter values
  LV.solve_state( p );
#if defined( SAVE_RESULTS )
  std::ofstream direcSTA;
  direcSTA.open( "test0_STA.dat", std::ios_base::out );
  LV.record( direcSTA );
#endif

  std::cout << "\nCONTINUOUS-TIME STATE SENSITIVITY INTEGRATION:\n\n";
  LV.solve_sensitivity( p );
#if defined( SAVE_RESULTS )
  std::ofstream direcFSA[NP];
  char fname[50];
  direcSTA.open( "test0_STA.dat", std::ios_base::out );
  for( unsigned i=0; i<NP; ++i ){
    sprintf( fname, "test0_FSA%d.dat",i );  
    direcFSA[i].open( fname, std::ios_base::out );
  }
  LV.record( direcSTA, direcFSA );
#endif
  
  std::cout << "\nCONTINUOUS-TIME ADJOINT SENSITIVITY INTEGRATION:\n\n";
  LV.solve_adjoint( p );
#if defined( SAVE_RESULTS )
  std::ofstream direcASA[NF];
  direcSTA.open( "test0_STA.dat", std::ios_base::out );
  for( unsigned i=0; i<NF; ++i ){
    sprintf( fname, "test0_ASA%d.dat",i );  
    direcASA[i].open( fname, std::ios_base::out );
  }
  LV.record( direcSTA, direcASA );
#endif

  /////////////////////////////////////////////////////////////////////////
  // Differentiate ODE system

  mc::ODESLVS_CVODES* GradLV = LV.fdiff( P );
  GradLV->setup();
  GradLV->solve_state( p );
  delete GradLV;

  return 0;
}


