#include "odeslvs_cvodes.hpp"

int main()
{
  /////////////////////////////////////////////////////////////////////////
  // Define IVP-ODE

  mc::FFGraph IVPDAG;  // DAG describing the problem

  double t0 = 0., tf = 10.;     // Time span
  const unsigned NS = 1;  // Time stages

  const unsigned NP = 1;  // Number of parameters
  const unsigned NX = 2;  // Number of states
  const unsigned NQ = 1;  // Number of state quadratures
  const unsigned NF = 2;  // Number of state functions

  mc::FFVar P[NP];  // Parameters
  for( unsigned int i=0; i<NP; i++ ) P[i].set( &IVPDAG );

  mc::FFVar X[NX];  // States
  for( unsigned int i=0; i<NX; i++ ) X[i].set( &IVPDAG );

  mc::FFVar Q[NQ];  // State quadratures
  for( unsigned i=0; i<NQ; i++ ) Q[i].set( &IVPDAG );

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

  mc::FFVar FCT[NF];  // State functions
  FCT[0] = X[0] * X[1];
  FCT[1] = P[0] * pow( X[0], 2 ) + Q[0];
/*
  mc::FFVar FCT[NF*NS];  // State functions
  for( unsigned k=0; k<NF*NS; k++ ) FCT[k] = 0.;
  if( NS > 1 ) FCT[((NS-1)/NF)*NF+0] = X[0] + 0.1*P[0];
  FCT[(NS-1)*NF+0] = X[0] * X[1];
  FCT[(NS-1)*NF+1] = P[0] * pow( X[0], 2 );
  for( unsigned k=0; k<NS; k++ ) FCT[k*NF+1] += Q[0];
*/

  mc::ODESLVS_CVODES IVP;

  IVP.options.INTMETH   = mc::BASE_CVODES::Options::MSBDF;//MSADAMS;
  IVP.options.NLINSOL   = mc::BASE_CVODES::Options::NEWTON;//FIXEDPOINT;
  IVP.options.LINSOL    = mc::BASE_CVODES::Options::DENSE;//DIAG;
  IVP.options.FSACORR   = mc::BASE_CVODES::Options::STAGGERED;//STAGGERED1;//SIMULTANEOUS;
  IVP.options.NMAX      = 2000;
  IVP.options.DISPLAY   = 1;
  IVP.options.ATOL      = IVP.options.ATOLB     = IVP.options.ATOLS  = 1e-9;
  IVP.options.RTOL      = IVP.options.RTOLB     = IVP.options.RTOLS  = 1e-9;
  IVP.options.QERR      = IVP.options.QERRS     = 1;
  IVP.options.ASACHKPT  = 1000;
#if defined( SAVE_RESULTS )
  LV.options.RESRECORD = 100;
#endif
  
  IVP.set_dag( &IVPDAG );
  IVP.set_state( NX, X );
  IVP.set_time( t0, tf );
  IVP.set_parameter( NP, P );
  IVP.set_differential( NX, RHS );
  IVP.set_initial( NX, IC );
  //IVP.set_initial( NS, NX, IC );
  IVP.set_quadrature( NQ, QUAD, Q );
  IVP.set_function( NF, FCT );
  //IVP.set_function( NS, NF, FCT );
  IVP.setup();

  /////////////////////////////////////////////////////////////////////////
  // Define DAG

  mc::FFGraph< mc::FFODE<0>, mc::FFGRADODE<0> > DAG;
  mc::FFVar PP[NP];  // Parameters
  for( unsigned int i=0; i<NP; i++ ) PP[i].set( &DAG );
  mc::FFODE<0> OpODE;
  mc::FFVar F[NF];
  mc::ODESLVS_CVODES<>* pIVP = &IVP;
  for( unsigned int j=0; j<NF; j++ ) F[j] = OpODE( j, NP, PP, pIVP );
  std::cout << DAG;
  
  auto F_op  = DAG.subgraph( NF, F );
  DAG.output( F_op );

  std::ofstream o_F( "test1b_F.dot", std::ios_base::out );
  DAG.dot_script( NF, F, o_F );
  o_F.close();

  // DAG evaluation in double arithmetic
  double dP[NP] = { 2.95 }, dF[NF];
  std::vector<double> dwk;
  DAG.eval( F_op, dwk, NF, F, dF, NP, PP, dP );
  for( unsigned i=0; i<NF; i++ )
    std::cout << "F[" << i << "] = " << dF[i] << std::endl;

#if defined( SAVE_RESULTS )
  std::ofstream of_state;
  of_state.open( "test1b_state.dat", std::ios_base::out );
  pIVP->record( of_state );
#endif

  // DAG evaluation in fadbad<double> arithmetic
  fadbad::F<double> FdP[NP] = { 2.95 }, FdF[NF];
  for( unsigned j=0; j<NP; j++ ) FdP[j].diff(j,NP);
  std::vector<fadbad::F<double>> Fdwk;
  DAG.eval( F_op, Fdwk, NF, F, FdF, NP, PP, FdP );
  for( unsigned i=0; i<NF; i++ ){
    std::cout << "F[" << i << "] = " << FdF[i].x() << std::endl;
    for( unsigned j=0; j<NP; j++ )
      std::cout << "dFdP[" << i << "][" << j << "] = " << FdF[i].d(j) << std::endl;
  }

  // DAG differentiation
  mc::FFVar const* dFdP = DAG.BAD( NF, F, NP, PP );
  std::cout << DAG;
  
  auto dFdP_op = DAG.subgraph( NF*NP, dFdP );
  DAG.output( dFdP_op, " dFdP" );

  std::ofstream o_dFdP( "test1b_dFdP.dot", std::ios_base::out );
  DAG.dot_script( NF*NP, dFdP, o_dFdP );
  o_dFdP.close();

  double ddFdP[NF*NP];
  DAG.eval( dFdP_op, dwk, NF*NP, dFdP, ddFdP, NP, PP, dP );
  for( unsigned i=0; i<NF; i++ )
    for( unsigned j=0; j<NP; j++ )
      std::cout << "dFdP[" << i << "][" << j << "] = " << ddFdP[i+j*NF] << std::endl;

  delete[] dFdP;
  return 0;
}


