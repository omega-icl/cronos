const unsigned int NPM   = 3;	// <- Order of poynomial expansion
const unsigned int NSAMP = 50;	// <- Number of sampling points for inner approx.
#define SAVE_RESULTS            // <- Whether to save bounds to file
#define USE_CMODEL              // <- whether to use Chebyshev models or Taylor models

#include "interval.hpp"
typedef mc::Interval I;

#ifdef USE_CMODEL
  #include "cmodel.hpp"
  typedef mc::CModel<I> PM;
  typedef mc::CVar<I> PV;
#else
  #include "tmodel.hpp"
  typedef mc::TModel<I> PM;
  typedef mc::TVar<I> PV;
#endif

#include "odebnds_sundials.hpp"

int main()
{
  mc::FFGraph IVP;  // DAG describing the problem

  double t0 = 0., tf = 6.;  // Time span
  const unsigned int NS = 120;  // Time stages
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

  mc::FFVar QUAD[NQ];  // Quadrature function
  QUAD[0] = X[1];

  //mc::FFVar FCT[NF];  // State functions
  //FCT[0] = X[0] * X[1];
  //FCT[1] = P[0] * pow( X[0], 2 );

  mc::FFVar FCT[NF*NS];  // State functions
  for( unsigned k=0; k<NF*NS; k++ ) FCT[k] = 0.;
  //if( NS > 1 ) FCT[((NS-1)/NF)*NF+0] = X[0] + 0.1*P[0];
  FCT[(NS-1)*NF+0] = X[0] * X[1];
  FCT[(NS-1)*NF+1] = P[0] * pow( X[0], 2 );
  for( unsigned k=0; k<NS; k++ ) FCT[k*NF+1] += Q[0];

  I Ip[NP] = { I(2.95,3.05) };
  I *Ixk[NS+1], *Iyk[NS+1], *Ixpk[NS+1];
  for( unsigned k=0; k<=NS; k++ ){
    Ixk[k] = new I[NX];
    Iyk[k] = new I[NF*NX];
    Ixpk[k] = new I[NP*NX];
  }
  I If[NF], Ifp[NF*NP];

  PM PMEnv( NP, NPM );
  PV PMp[NP] = { PV( &PMEnv, 0, Ip[0] ) };
  PV *PMxk[NS+1], *PMyk[NS+1], *PMxpk[NS+1];
  for( unsigned k=0; k<=NS; k++ ){
    PMxk[k] = new PV[NX];
    PMyk[k] = new PV[NF*NX];
    PMxpk[k] = new PV[NP*NX];
  }
  PV PMf[NF], PMfp[NF*NP];

  /////////////////////////////////////////////////////////////////////////////
  //// DIFFERENTIAL INEQUALITIES
  mc::ODEBNDS_SUNDIALS<I,PM,PV> LV;

  LV.options.DISPLAY   = 1;
#if defined( SAVE_RESULTS )
  LV.options.RESRECORD = true;
#endif
  LV.options.ATOL      = LV.options.ATOLB      = LV.options.ATOLS  = 1e-20;
  LV.options.RTOL      = LV.options.RTOLB      = LV.options.RTOLS  = 1e-8;
  LV.options.NMAX      = 20000;
  LV.options.MAXFAIL   = 10;
  LV.options.FSAERR    = true;
  LV.options.QERRS     = true;
  //LV.options.AUTOTOLS  = true;
  LV.options.ASACHKPT  = 20000;
  LV.options.INTMETH   = mc::ODEBNDS_SUNDIALS<I,PM,PV>::Options::MSBDF;//MSADAMS;//MSBDF;
  LV.options.JACAPPROX = mc::ODEBNDS_SUNDIALS<I,PM,PV>::Options::CV_DIAG;//CV_DENSE;//CV_DIAG;
  LV.options.ORDMIT    = -2; //PMp->nord();
  LV.options.WRAPMIT   = mc::ODEBNDS_SUNDIALS<I,PM,PV>::Options::ELLIPS;//NONE;//DINEQ;
  LV.options.QERRB     = true;
  LV.options.QSCALE    = 1e-5;

  LV.options.ODESLV.NMAX = 20000;
  LV.options.ODESLV.INTMETH   = mc::ODESLVS_SUNDIALS::Options::MSBDF;//MSADAMS;//MSBDF;
  LV.options.ODESLV.JACAPPROX = mc::ODESLVS_SUNDIALS::Options::CV_DIAG;//CV_DENSE;//CV_DIAG;

  LV.set_dag( &IVP );
  LV.set_state( NX, X );
  LV.set_parameter( NP, P );
  LV.set_differential( NX, RHS );
  LV.set_initial( NX, IC );
  LV.set_quadrature( NQ, QUAD, Q );
  //LV.set_function( NF, FCT );
  LV.set_function( NS, NF, FCT );

  std::ofstream ofSTA, ofSEN;

  std::cout << "\nNON_VALIDATED INTEGRATION - INNER-APPROXIMATION OF REACHABLE SET AND ADJOINT SENSITIVITY:\n\n";
  LV.bounds_ASA( NS, tk, Ip, Ixk, If, Iyk, Ifp, NSAMP );
#if defined( SAVE_RESULTS )
  ofSTA.open( "test1_APPROX_STA.dat", std::ios_base::out );
  ofSEN.open( "test1_APPROX_ASA.dat", std::ios_base::out );
  LV.record( ofSTA, ofSEN );
  ofSTA.close();
  ofSEN.close();
#endif
  { int dum; std::cout << "PAUSED--"; std::cin >> dum; }
/*
  std::cout << "\nCONTINUOUS SET-VALUED INTEGRATION - INTERVAL ENCLOSURE OF REACHABLE SET AND ADJOINT SENSITIVITY:\n\n";
  LV.bounds_ASA( NS, tk, Ip, Ixk, If, Iyk, Ifp );
#if defined( SAVE_RESULTS )
  ofSTA.open( "test1_DINEQI_STA.dat", std::ios_base::out );
  ofSEN.open( "test1_DINEQI_ASA.dat", std::ios_base::out );
  LV.record( ofSTA, ofSEN );
  ofSTA.close();
  ofSEN.close();
#endif
*/
  std::cout << "\nCONTINUOUS SET-VALUED INTEGRATION - POLYNOMIAL MODEL ENCLOSURE OF REACHABLE SET AND ADJOINT SENSITIVITY:\n\n";
  //LV.bounds( NS, tk, PMp, PMxk, PMf );
  LV.bounds_ASA( NS, tk, PMp, PMxk, PMf, PMyk, PMfp );
#if defined( SAVE_RESULTS )
  ofSTA.open( "test1_DINEQPM_STA.dat", std::ios_base::out );
  ofSEN.open( "test1_DINEQPM_ASA.dat", std::ios_base::out );
  LV.record( ofSTA, ofSEN );
  ofSTA.close();
  ofSEN.close();
#endif
  { int dum; std::cout << "--PAUSED "; std::cin >> dum; }

  std::cout << "\nNON_VALIDATED INTEGRATION - INNER-APPROXIMATION OF REACHABLE SET AND FORWARD SENSITIVITY:\n\n";
  LV.bounds_FSA( NS, tk, Ip, Ixk, If, Ixpk, Ifp, NSAMP );
#if defined( SAVE_RESULTS )
  ofSTA.open( "test1_APPROX_STA.dat", std::ios_base::out );
  ofSEN.open( "test1_APPROX_FSA.dat", std::ios_base::out );
  LV.record( ofSTA, ofSEN );
  ofSTA.close();
  ofSEN.close();
#endif
  { int dum; std::cout << "PAUSED--"; std::cin >> dum; }
/*
  std::cout << "\nCONTINUOUS SET-VALUED INTEGRATION - INTERVAL ENCLOSURE OF REACHABLE SET AND FORWARD SENSITIVITY:\n\n";
  //LV.bounds( NS, tk, Ip, Ixk, If );
  LV.bounds_FSA( NS, tk, Ip, Ixk, If, Ixpk, Ifp );
#if defined( SAVE_RESULTS )
  ofSTA.open( "test1_DINEQI_STA.dat", std::ios_base::out );
  ofSEN.open( "test1_DINEQI_FSA.dat", std::ios_base::out );
  LV.record( ofSTA, ofSEN );
  ofSTA.close();
  ofSEN.close();
#endif
*/
  std::cout << "\nCONTINUOUS SET-VALUED INTEGRATION - POLYNOMIAL MODEL ENCLOSURE OF REACHABLE SET AND FORWARD SENSITIVITY:\n\n";
  //LV.bounds( NS, tk, PMp, PMxk, PMf );
  LV.bounds_FSA( NS, tk, PMp, PMxk, PMf, PMxpk, PMfp );
#if defined( SAVE_RESULTS )
  ofSTA.open( "test1_DINEQPM_STA.dat", std::ios_base::out );
  ofSEN.open( "test1_DINEQPM_FSA.dat", std::ios_base::out );
  LV.record( ofSTA, ofSEN );
  ofSTA.close();
  ofSEN.close();
#endif

  // Clean up
  for( unsigned k=0; k<=NS; k++ ){
    delete[] Ixk[k];
    delete[] Iyk[k];
    delete[] Ixpk[k];
    delete[] PMxk[k];
    delete[] PMyk[k];
    delete[] PMxpk[k];
  }

  return 0;
}
