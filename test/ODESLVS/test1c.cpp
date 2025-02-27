#define SAVE_RESULTS		// <- Whether to save bounds to file
#define CRONOS__ODESLVS_FDIFF_DEBUG

#include "ffode.hpp"

int main()
{
  /////////////////////////////////////////////////////////////////////////
  // Define IVP-ODE

  mc::FFGraph DAG;  // DAG describing the problem

  const unsigned NS = 4;  // Time stages
  double t0 = 0., tf = 10.;       // Time span
  std::vector<double> T( NS+1 );  // Time stages
  for( unsigned int i=0; i<=NS; i++ ) T[i] = t0 + i * ( tf - t0 ) / NS; 

  const unsigned NC = 1;  // Number of constants
  std::vector<mc::FFVar> C(NC);   // Parameters
  for( unsigned int i=0; i<NC; i++ ) C[i].set( &DAG );

  const unsigned NP = 2;  // Number of parameters
  std::vector<mc::FFVar> P(NP);   // Parameters
  for( unsigned int i=0; i<NP; i++ ) P[i].set( &DAG );

  const unsigned NX = 2;  // Number of states
  std::vector<mc::FFVar> X(NX);   // States
  for( unsigned int i=0; i<NX; i++ ) X[i].set( &DAG );

  const unsigned NQ = 1;  // Number of state quadratures
  std::vector<mc::FFVar> Q(NQ);   // State quadratures
  for( unsigned i=0; i<NQ; i++ ) Q[i].set( &DAG );

  std::vector<mc::FFVar> RHS(NX); // Right-hand side function
  RHS[0] = P[0] * X[0] * ( 1. - X[1] );
  RHS[1] = P[0] * X[1] * ( X[0] - 1. );

//  std::vector<mc::FFVar> IC(NX);  // Initial value function
//  IC[0] = 1.2;
//  IC[1] = 1.1 + 0.01*P[1];

  std::vector<std::vector<mc::FFVar>> IC(NS);   // Initial value function
  IC[0].assign( { C[0], 1.1 + 0.01*P[1] } );
  for( unsigned k=1; k<NS; k++ )
     IC[k].assign( { X[0], X[1] - 0.02*P[1] } );

  std::vector<mc::FFVar> QUAD(NQ);  // Quadrature function
  QUAD[0] = X[1];

  const unsigned NF = 2;  // Number of state functions

//  std::vector<mc::FFVar> FCT(NF);  // State functions
//  FCT[0] = X[0] * X[1];
//  FCT[1] = Q[0]; //P[0] * pow( X[0], 2 );
//  //FCT[0] = Q[0]; //P[0] * pow( X[0], 2 );

  std::vector<std::map<size_t,mc::FFVar>> FCT(NS+1);  // State functions
  for( unsigned k=1; k<NS; k++ )
    FCT[k] = { { 1, Q[0] } };
  FCT[NS] = { { 0, X[0] * X[1] }, { 1, Q[0] } };


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
  IVP.options.RESRECORD = 100;
#endif
  
  IVP.set_dag( &DAG );
  IVP.set_time( T );
  IVP.set_state( X );
  IVP.set_constant( C );
  IVP.set_parameter( P );
  IVP.set_differential( RHS );
  IVP.set_initial( IC );
  IVP.set_quadrature( QUAD, Q );
  IVP.set_function( FCT );
  IVP.setup();

  /////////////////////////////////////////////////////////////////////////
  // Insert IVP-ODE into DAG

  mc::FFODE OpODE;
  std::vector<mc::FFVar> F(NF);
  for( unsigned int j=0; j<NF; j++ ) F[j] = OpODE( j, NP, P.data(), NC, C.data(), &IVP, 0 );
  std::cout << DAG;

  auto F_op  = DAG.subgraph( NF, F.data() );
  DAG.output( F_op );

  std::ofstream o_F( "test1c_F.dot", std::ios_base::out );
  DAG.dot_script( NF, F.data(), o_F );
  o_F.close();

  // DAG evaluation in double arithmetic
  std::vector<double> dP( { 2.96, 3. } ), dC( { 1.2 } ), dF(NF), dwk;
/*
  DAG.eval( F_op, dwk, NF, F.data(), dF.data(), NP, P.data(), dP.data(), NC, C.data(), dC.data() );
  for( unsigned i=0; i<NF; i++ )
    std::cout << "F[" << i << "] = " << dF[i] << std::endl;

#if defined( SAVE_RESULTS )
  std::ofstream of_state;
  of_state.open( "test1c_STA.dat", std::ios_base::out );
  IVP.record( of_state );
#endif
*/
  // DAG symbolic differentiation
  mc::FFODE::options.DIFF = mc::FFODE::Options::SYM_C;
  mc::FFVar const* dFdC = DAG.FAD( NF, F.data(), NC, C.data() );
  //std::cout << DAG;

  auto dFdC_op  = DAG.subgraph( NF*NC, dFdC );
  //DAG.output( dFdC_op );

  //std::ofstream o_dFdC( "test1c_dFdC.dot", std::ios_base::out );
  //DAG.dot_script( NF*NC, dFdC, o_dFdC );
  //o_dFdC.close();

  // DAG evaluation in double arithmetic
  std::vector<double> ddFdC(NF*NC);
  DAG.eval( dFdC_op, dwk, NF*NC, dFdC, ddFdC.data(), NP, P.data(), dP.data(), NC, C.data(), dC.data() );
  for( unsigned i=0; i<NF; i++ )
    for( unsigned j=0; j<NC; j++ )
      std::cout << "dF[" << i << "]dC[" << j << "] = " << ddFdC[i*NC+j] << std::endl;

  delete[] dFdC;

  // DAG symbolic differentiation
  mc::FFODE::options.DIFF = mc::FFODE::Options::SYM_P;
  mc::FFVar const* dFdP = DAG.FAD( NF, F.data(), NP, P.data() );
  //std::cout << DAG;

  auto dFdP_op  = DAG.subgraph( NF*NP, dFdP );
  //DAG.output( dFdP_op );

  //std::ofstream o_dFdP( "test1c_dFdP.dot", std::ios_base::out );
  //DAG.dot_script( NF*NP, dFdP, o_dFdP );
  //o_dFdP.close();

  // DAG evaluation in double arithmetic
  std::vector<double> ddFdP(NF*NP);
  DAG.eval( dFdP_op, dwk, NF*NP, dFdP, ddFdP.data(), NP, P.data(), dP.data(), NC, C.data(), dC.data() );
  for( unsigned i=0; i<NF; i++ )
    for( unsigned j=0; j<NP; j++ )
      std::cout << "dF[" << i << "]dP[" << j << "] = " << ddFdP[i*NP+j] << std::endl;

  delete[] dFdP;

  // DAG symbolic differentiation
  mc::FFODE::options.DIFF = mc::FFODE::Options::SYM_PC;
  mc::FFVar const* dFdPC = DAG.FAD( NF, F.data(), NP, P.data(), NC, C.data() );
  //std::cout << DAG;

  auto dFdPC_op  = DAG.subgraph( NF*(NP+NC), dFdPC );
  //DAG.output( dFdPC_op );

  //std::ofstream o_dFdPC( "test1c_dFdPC.dot", std::ios_base::out );
  //DAG.dot_script( NF*(NP+NC), dFdPC, o_dFdPC );
  //o_dFdPC.close();

  // DAG evaluation in double arithmetic
  std::vector<double> ddFdPC(NF*(NP+NC));
  DAG.eval( dFdPC_op, dwk, NF*(NP+NC), dFdPC, ddFdPC.data(), NP, P.data(), dP.data(), NC, C.data(), dC.data() );
  for( unsigned i=0; i<NF; i++ ){
    for( unsigned j=0; j<NP; j++ )
      std::cout << "dF[" << i << "]dP[" << j << "] = " << ddFdPC[i*(NP+NC)+j] << std::endl;
    for( unsigned j=0; j<NC; j++ )
      std::cout << "dF[" << i << "]dC[" << j << "] = " << ddFdPC[i*(NP+NC)+NP+j] << std::endl;
  }
  delete[] dFdPC;

  // DAG numerical differentiation
  mc::FFODE::options.DIFF = mc::FFODE::Options::NUM_P;
  dFdP = DAG.BAD( NF, F.data(), NP, P.data() );
  //std::cout << DAG;

  dFdP_op = DAG.subgraph( NF*NP, dFdP );
  //DAG.output( dFdP_op, " dFdP" );

  //o_dFdP.open( "test1c_dFdP.dot", std::ios_base::out );
  //DAG.dot_script( NF*NP, dFdP, o_dFdP );
  //o_dFdP.close();

  //std::vector<double> ddFdP(NF*NP);
  DAG.eval( dFdP_op, dwk, NF*NP, dFdP, ddFdP.data(), NP, P.data(), dP.data(), NC, C.data(), dC.data() );
  for( unsigned i=0; i<NF; i++ )
    for( unsigned j=0; j<NP; j++ )
      std::cout << "dF[" << i << "]dP[" << j << "] = " << ddFdP[i*NP+j] << std::endl;

  delete[] dFdP;

  // DAG evaluation in fadbad<double> arithmetic
  std::vector<fadbad::F<double>> FdP( { 2.96, 3. } ), FdC( { 1.2 } ), FdF(NF), Fdwk;
  for( unsigned j=0; j<NP; j++ ) FdP[j].diff(j,NP);
  DAG.eval( F_op, Fdwk, NF, F.data(), FdF.data(), NP, P.data(), FdP.data(), NC, C.data(), FdC.data() );
  for( unsigned i=0; i<NF; i++ ){
    std::cout << "F[" << i << "] = " << FdF[i].x() << std::endl;
    for( unsigned j=0; j<NP; j++ )
      std::cout << "dF[" << i << "]dP[" << j << "] = " << FdF[i].d(j) << std::endl;
  }

  return 0;
}


