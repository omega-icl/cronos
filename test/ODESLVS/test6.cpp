#define SAVE_RESULTS		// <- Whether to save bounds to file

#include "odeslvs_cvodes.hpp"

int main()
{
  mc::FFGraph IVPDAG;  // DAG describing the problem

  double T0 = 0., TF = 20.;   // Time span
  size_t const NS = 1;      // Time stages
  std::vector<double> TS( NS+1 );  // Time stages
  for( size_t i=0; i<=NS; i++ ) TS[i] = T0 + i * ( TF - T0 ) / NS; 
  mc::FFVar T; T.set( &IVPDAG );  // Time
  
  double V0 = -30., VF = 20.; // Space range
  size_t const NM = 101; // Mesh stages -- Keep this as an odd number
  double H = (VF-V0)/(double)NM;

  size_t NP = 2;      // Number of parameters
  std::vector<mc::FFVar> P(NP);   // Parameters
  for( size_t i=0; i<NP; i++ ) P[i].set( &IVPDAG );

  size_t const NU = 2;      // Number of states
  size_t const NX = NU*NM;  // Number of discretized states (ODEs)
  std::vector<mc::FFVar> X(NX);   // States
  for( size_t i=0; i<NX; i++ ) X[i].set( &IVPDAG );
  
  // Parameters
  double    rc = 10.,     bc = 1.25; 
  mc::FFVar r  = rc*P[0], b  = bc*P[1];
  //double    pc = 2e-2;
  
  // Map u -> x
  std::vector<mc::FFVar> U(NU*(NM+2));   // Discretized states
  for( size_t i=0; i<NM; i++ ){ U[i+1] = X[i]; U[NM+3+i] = X[NM+i]; }

  // Boundary Conditions
  U[0]    = U[1];    U[NM+1]        = U[NM];          // u_1
  U[NM+2] = U[NM+3]; U[NU*(NM+2)-1] = U[NU*(NM+2)-2]; // u_2

  std::vector<mc::FFVar> IC(NX);  // Initial value function
  IC[NM/2] = IC[NM/2 + NM] = 0.5;
  for( size_t i=0     ; i<NM/2; i++ ){ IC[i] = 0.; IC[i+NM] = 1.; }
  for( size_t i=NM/2+1; i<NM  ; i++ ){ IC[i] = 1.; IC[i+NM] = 0.; }

  std::vector<mc::FFVar> RHS(NX); // Right-hand side function
  for( size_t i=1; i<NM+1; i++ ){
    RHS[i-1]    = ( U[i+1]    - 2*U[i]      + U[i-1] )   /( H*H ) + U[i]*( 1. - U[i] - r*U[i+NM+2] );
    RHS[i+NM-1] = ( U[i+NM+3] - 2*U[i+NM+2] + U[i+NM+1] )/( H*H ) - b * U[i] * U[i+NM+2];
  }

  size_t const NF = 2;      // Number of state functions
  std::vector<mc::FFVar> FCT(NF);  // State functions
  FCT[0] = FCT[1] = 0.;
  for( size_t i=1; i<NM+1; i++ ){
    FCT[0] += U[i]      / double(NM);
    FCT[1] += U[i+NM+1] / double(NM);
  }

  /////////////////////////////////////////////////////////////////////////
  // Compute ODE solutions

  mc::ODESLVS_CVODES IVP;

  IVP.options.INTMETH   = mc::BASE_CVODES::Options::MSBDF;//MSADAMS;
  IVP.options.NLINSOL   = mc::BASE_CVODES::Options::NEWTON;//FIXEDPOINT;
  IVP.options.LINSOL    = mc::BASE_CVODES::Options::SPARSE;//DENSE;//DENSEDQ;//DIAG;
  IVP.options.NMAX      = 20000;
  IVP.options.DISPLAY   = 1;
  IVP.options.ATOL      = IVP.options.ATOLB      = IVP.options.ATOLS  = 1e-9;
  IVP.options.RTOL      = IVP.options.RTOLB      = IVP.options.RTOLS  = 1e-8;
#if defined( SAVE_RESULTS )
  IVP.options.RESRECORD = 40;
#endif

  IVP.set_dag( &IVPDAG );
  IVP.set_time( TS, &T );
  IVP.set_state( X );
  IVP.set_parameter( P );
  IVP.set_differential( RHS );
  IVP.set_initial( IC );
  IVP.set_function( FCT );
  IVP.setup();

  std::vector<double> p0( { 1., 1. } );

  std::ofstream direcSTA, direcFSA[NP], direcASA[NF];
  char fname[50];

  std::cout << "\nCONTINUOUS-TIME INTEGRATION:\n\n";
  IVP.solve_state( p0 );
#if defined( SAVE_RESULTS )
  direcSTA.open( "test6_STA.dat", std::ios_base::out );
  IVP.record( direcSTA );
#endif

  std::cout << "\nCONTINUOUS-TIME INTEGRATION WITH FORWARD SENSITIVITY ANALYSIS:\n\n";
  IVP.solve_sensitivity( p0 );
#if defined( SAVE_RESULTS )
  direcSTA.open( "test6_STA.dat", std::ios_base::out );
  for( unsigned i=0; i<NP; ++i ){
    sprintf( fname, "test6_FSA%d.dat",i );  
    direcFSA[i].open( fname, std::ios_base::out );
  }
  IVP.record( direcSTA, direcFSA );
#endif

  std::cout << "\nCONTINUOUS-TIME INTEGRATION WITH ADJOINT SENSITIVITY ANALYSIS:\n\n";
  IVP.solve_adjoint( p0 );
#if defined( SAVE_RESULTS )
  direcSTA.open( "test6_STA.dat", std::ios_base::out );
  for( unsigned i=0; i<NF; ++i ){
    sprintf( fname, "test6_ASA%d.dat",i );  
    direcASA[i].open( fname, std::ios_base::out );
  }
  IVP.record( direcSTA, direcASA );
#endif

  return 0;
}

