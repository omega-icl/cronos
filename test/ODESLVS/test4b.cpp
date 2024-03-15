#include "odeslvs_cvodes.hpp"

int main()
{
  /////////////////////////////////////////////////////////////////////////
  // Define IVP-ODE

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
  IVP.options.RESRECORD = 20;
  
  IVP.set_dag( &IVPDAG );
  IVP.set_time( NS, tk );
  IVP.set_state( NX, X );
  IVP.set_parameter( NP, P );
  IVP.set_differential( NS, NX, RHS );
  IVP.set_initial( NX, IC );
  IVP.set_quadrature( NS, NQ, QUAD, Q );
  IVP.set_function( NF, FCT );
  IVP.setup();

  /////////////////////////////////////////////////////////////////////////
  // Define DAG

  mc::FFGraph< mc::FFODE<0>, mc::FFGRADODE<0> > DAG;

  mc::FFVar PP[NP];  // Parameters
  for( unsigned int i=0; i<NP; i++ ) PP[i].set( &DAG );

  mc::FFODE<0> OpODE;
  mc::ODESLVS_CVODES<>* pIVP = &IVP;
  mc::FFVar F[NF];
  for( unsigned int j=0; j<NF; j++ ) F[j] = OpODE( j, NP, PP, pIVP );
  std::cout << DAG;
  
  auto F_op  = DAG.subgraph( NF, F );
  DAG.output( F_op );

  std::ofstream o_F( "test4b_F.dot", std::ios_base::out );
  DAG.dot_script( NF, F, o_F );
  o_F.close();

  // DAG evaluation in double arithmetic
  double dP[NP], dF[NF];
  for( unsigned int i=0; i<NP; i++ ) dP[i] = -3e-1;
  std::vector<double> dwk;
  DAG.eval( F_op, dwk, NF, F, dF, NP, PP, dP );
  for( unsigned i=0; i<NF; i++ )
    std::cout << "F[" << i << "] = " << dF[i] << std::endl;

  std::ofstream of_state;
  of_state.open( "test4b_state.dat", std::ios_base::out );
  pIVP->record( of_state );

  // DAG evaluation in fadbad<double> arithmetic
  fadbad::F<double> FdP[NP], FdF[NF];
  for( unsigned j=0; j<NP; j++ ){ FdP[j] = -3e-1; FdP[j].diff(j,NP); }
  std::vector<fadbad::F<double>> Fdwk;
  DAG.eval( F_op, Fdwk, NF, F, FdF, NP, PP, FdP );
  for( unsigned i=0; i<NF; i++ ){
    std::cout << "F[" << i << "] = " << FdF[i].x() << std::endl;
    for( unsigned j=0; j<NP; j++ )
      std::cout << "dFdP[" << i << "][" << j << "] = " << FdF[i].d(j) << std::endl;
  }

  // DAG forward AD
  mc::FFVar const* dFdP = DAG.FAD( NF, F, NP, PP );
  std::cout << DAG;
  
  auto dFdP_op = DAG.subgraph( NF*NP, dFdP );
  DAG.output( dFdP_op, " dFdP" );

  std::ofstream o_dFdP( "test4b_dFdP.dot", std::ios_base::out );
  DAG.dot_script( NF*NP, dFdP, o_dFdP );
  o_dFdP.close();

  double ddFdP[NF*NP];
  DAG.eval( dFdP_op, dwk, NF*NP, dFdP, ddFdP, NP, PP, dP );
  for( unsigned i=0; i<NF; i++ )
    for( unsigned j=0; j<NP; j++ )
      std::cout << "dFdP[" << i << "][" << j << "] = " << ddFdP[j+i*NP] << std::endl;

  delete[] dFdP;
  return 0;
}


