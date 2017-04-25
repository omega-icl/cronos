/// LOTKA-VOLTERRA SYSTEM ///
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

  double T0 = 0., TF = 8.;      // Time span
  const unsigned int NS = 1;    // Time stages
  double TS[NS+1]; TS[0] = T0;
  for( unsigned k=0; k<NS; k++ ) TS[k+1] = TS[k] + (TF-T0)/(double)NS;

  const unsigned NP = 1;  // Number of parameters
  const unsigned NX = 2;  // Number of states
  const unsigned NQ = 1;  // Number of state quadratures
  const unsigned NF = 2;  // Number of state functions

  mc::FFVar P[NP];    // Parameters
  for( unsigned int i=0; i<NP; i++ ) P[i].set( &IVP );

  mc::FFVar X[NX];    // States
  for( unsigned int i=0; i<NX; i++ ) X[i].set( &IVP );

  mc::FFVar Q[NQ];    // State quadratures
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
  //I *Ixk[NS+1], *Iyk[NS+1], *Ixpk[NS+1];
  //for( unsigned k=0; k<=NS; k++ ){
  //  Ixk[k] = new I[NX+NQ];
  //  Iyk[k] = new I[NF*(NX+NQ)];
  //  Ixpk[k] = new I[NP*(NX+NQ)X];
  //}
  //I If[NF], Ifp[NF*NP];

  PM PMEnv( NP, NPM );
  PV PMp[NP] = { PV( &PMEnv, 0, Ip[0] ) };
  //PV *PMxk[NS+1], *PMyk[NS+1], *PMxpk[NS+1];
  //for( unsigned k=0; k<=NS; k++ ){
  //  PMxk[k] = new PV[NX+NQ];
  //  PMyk[k] = new PV[NF*(NX+NQ)];
  //  PMxpk[k] = new PV[NP*(NX+NQ)];
  //}
  //PV PMf[NF], PMfp[NF*NP];

  std::ofstream ofSTA, ofFSA[NP], ofASA[NF];
  char fname[50];

  /////////////////////////////////////////////////////////////////////////////
  //// DIFFERENTIAL INEQUALITIES
  mc::ODEBNDS_SUNDIALS<I,PM,PV> LV;

  LV.options.DISPLAY   = 1;
#if defined( SAVE_RESULTS )
  LV.options.RESRECORD = 100;
#endif
  LV.options.ATOL      = LV.options.ATOLB      = LV.options.ATOLS  = 1e-8;
  LV.options.RTOL      = LV.options.RTOLB      = LV.options.RTOLS  = 1e-8;
  LV.options.ETOL      = LV.options.ETOLB      = LV.options.ETOLS  = 1e-20;
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
  LV.options.HMIN      = 1e-12;

  LV.options.ODESLVS.NMAX = 20000;
  LV.options.ODESLVS.INTMETH   = mc::ODESLVS_SUNDIALS::Options::MSBDF;//MSADAMS;//MSBDF;
  LV.options.ODESLVS.JACAPPROX = mc::ODESLVS_SUNDIALS::Options::CV_DIAG;//CV_DENSE;//CV_DIAG;

  LV.set_dag( &IVP );
  LV.set_time( NS, TS );
  LV.set_state( NX, X );
  LV.set_parameter( NP, P );
  LV.set_differential( NX, RHS );
  LV.set_initial( NX, IC );
  LV.set_quadrature( NQ, QUAD, Q );
  //LV.set_function( NF, FCT );
  LV.set_function( NS, NF, FCT );

  std::cout << "\nNON_VALIDATED INTEGRATION - INNER-APPROXIMATION OF REACHABLE SET AND ADJOINT SENSITIVITY:\n\n";
  LV.bounds_ASA( NSAMP, Ip );//, Ixk, If, Iyk, Ifp );
#if defined( SAVE_RESULTS )
  ofSTA.open( "test1_APPROX_STA.dat", std::ios_base::out );
  for( unsigned i=0; i<NF; ++i ){
    sprintf( fname, "test1_APPROX_ASA%d.dat",i );  
    ofASA[i].open( fname, std::ios_base::out );
  }
  LV.record( ofSTA, ofASA );
  ofSTA.close();
  for( unsigned i=0; i<NF; ++i ) ofASA[i].close();
#endif

  std::cout << "\nCONTINUOUS SET-VALUED INTEGRATION - INTERVAL ENCLOSURE OF REACHABLE SET AND ADJOINT SENSITIVITY:\n\n";
  LV.bounds_ASA( Ip );//, Ixk, If, Iyk, Ifp );
#if defined( SAVE_RESULTS )
  ofSTA.open( "test1_DINEQI_STA.dat", std::ios_base::out );
  for( unsigned i=0; i<NF; ++i ){
    sprintf( fname, "test1_DINEQI_ASA%d.dat",i );  
    ofASA[i].open( fname, std::ios_base::out );
  }
  LV.record( ofSTA, ofASA );
  ofSTA.close();
  for( unsigned i=0; i<NF; ++i ) ofASA[i].close();
#endif

  std::cout << "\nCONTINUOUS SET-VALUED INTEGRATION - POLYNOMIAL MODEL ENCLOSURE OF REACHABLE SET AND ADJOINT SENSITIVITY:\n\n";
  LV.bounds_ASA( PMp );//, PMxk, PMf, PMyk, PMfp );
#if defined( SAVE_RESULTS )
  ofSTA.open( "test1_DINEQPM_STA.dat", std::ios_base::out );
  for( unsigned i=0; i<NF; ++i ){
    sprintf( fname, "test1_DINEQPM_ASA%d.dat",i );  
    ofASA[i].open( fname, std::ios_base::out );
  }
  LV.record( ofSTA, ofASA );
  ofSTA.close();
  for( unsigned i=0; i<NF; ++i ) ofASA[i].close();
#endif


  std::cout << "\nNON_VALIDATED INTEGRATION - INNER-APPROXIMATION OF REACHABLE SET AND FORWARD SENSITIVITY:\n\n";
  LV.bounds_FSA( NSAMP, Ip );//, Ixk, If, Ixpk, Ifp );
#if defined( SAVE_RESULTS )
  ofSTA.open( "test1_APPROX_STA.dat", std::ios_base::out );
  for( unsigned i=0; i<NP; ++i ){
    sprintf( fname, "test1_APPROX_FSA%d.dat",i );  
    ofFSA[i].open( fname, std::ios_base::out );
  }
  LV.record( ofSTA, ofFSA );
  ofSTA.close();
  for( unsigned i=0; i<NP; ++i ) ofFSA[i].close();
#endif

  std::cout << "\nCONTINUOUS SET-VALUED INTEGRATION - INTERVAL ENCLOSURE OF REACHABLE SET AND FORWARD SENSITIVITY:\n\n";
  LV.bounds_FSA( Ip );//, Ixk, If, Ixpk, Ifp );
#if defined( SAVE_RESULTS )
  ofSTA.open( "test1_DINEQI_STA.dat", std::ios_base::out );
  for( unsigned i=0; i<NP; ++i ){
    sprintf( fname, "test1_DINEQI_FSA%d.dat",i );  
    ofFSA[i].open( fname, std::ios_base::out );
  }
  LV.record( ofSTA, ofFSA );
  ofSTA.close();
  for( unsigned i=0; i<NP; ++i ) ofFSA[i].close();
#endif

  std::cout << "\nCONTINUOUS SET-VALUED INTEGRATION - POLYNOMIAL MODEL ENCLOSURE OF REACHABLE SET AND FORWARD SENSITIVITY:\n\n";
  LV.bounds_FSA( PMp );//, PMxk, PMf, PMxpk, PMfp );
#if defined( SAVE_RESULTS )
  ofSTA.open( "test1_DINEQPM_STA.dat", std::ios_base::out );
  for( unsigned i=0; i<NP; ++i ){
    sprintf( fname, "test1_DINEQPM_FSA%d.dat",i );  
    ofFSA[i].open( fname, std::ios_base::out );
  }
  LV.record( ofSTA, ofFSA );
  ofSTA.close();
  for( unsigned i=0; i<NP; ++i ) ofFSA[i].close();
#endif

  // Clean up
  //for( unsigned k=0; k<=NS; k++ ){
  //  delete[] Ixk[k];
  //  delete[] Iyk[k];
  //  delete[] Ixpk[k];
  //  delete[] PMxk[k];
  //  delete[] PMyk[k];
  //  delete[] PMxpk[k];
  //}

  return 0;
}
