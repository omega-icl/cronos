/// CUBIC OSCILLATOR ///
const unsigned int NPM   = 5;	// <- Order of poynomial expansion
const unsigned int NSAMP = 20;	// <- Number of sampling points for inner approx.
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

  double T0 = 0., TF = 100.;  // Time span
  const unsigned int NS = 1;  // Time stages
  double TS[NS+1]; TS[0] = T0;
  for( unsigned k=0; k<NS; k++ ) TS[k+1] = TS[k] + (TF-T0)/(double)NS;

  const unsigned NP = 2;  // Number of parameters
  const unsigned NX = 2;  // Number of states
  const unsigned NQ = 1;  // Number of state quadratures
  const unsigned NF = 1;  // Number of state functions

  mc::FFVar P[NP];  // Parameters
  for( unsigned int i=0; i<NP; i++ ) P[i].set( &IVP );

  mc::FFVar X[NX];  // States
  for( unsigned int i=0; i<NX; i++ ) X[i].set( &IVP );

  mc::FFVar Q[NQ];  // State quadratures
  for( unsigned i=0; i<NQ; i++ ) Q[i].set( &IVP );

  mc::FFVar RHS[NX];  // Right-hand side function
  const double q = 0.1;
  RHS[0] = q*X[0]*(1.-X[0]*X[0]) + X[1]*(1.-q*X[0]*X[1]);
  RHS[1] = q*X[1]*(1.-X[1]*X[1]) - X[0]*(1.+q*X[0]*X[1]) - 0.2*X[1];

  mc::FFVar IC[NX];   // Initial value function
  IC[0] = P[0];
  IC[1] = P[1];

  mc::FFVar QUAD[NQ];  // Quadrature function
  QUAD[0] = X[0]*X[0];

  mc::FFVar FCT[NF*NS];  // State functions
  for( unsigned k=0; k<NS; k++ ) FCT[k*NF] = Q[0];

  // Parameter and state bounds
  I Ip[NP] = { I(1.5,2.5), I(-0.5,0.5) };
  //I *Ixk[NS+1], *Iyk[NS+1], *Ixpk[NS+1];
  //for( unsigned k=0; k<=NS; k++ ){
  //  Ixk[k] = new I[NX];
  //  Iyk[k] = new I[NF*NX];
  //  Ixpk[k] = new I[NP*NX];
  //}
  //I If[NF], Ifp[NF*NP];

  PM PMEnv( NP, NPM );
  PV PMp[NP] = { PV( &PMEnv, 0, Ip[0] ), PV( &PMEnv, 1, Ip[1] )  };
  //PV *PMxk[NS+1], *PMyk[NS+1], *PMxpk[NS+1];
  //for( unsigned k=0; k<=NS; k++ ){
  //  PMxk[k] = new PV[NX];
  //  PMyk[k] = new PV[NF*NX];
  //  PMxpk[k] = new PV[NP*NX];
  //}
  //PV PMf[NF], PMfp[NF*NP];

  std::ofstream ofSTA, ofFSA[NP], ofASA[NF];
  char fname[50];

  /////////////////////////////////////////////////////////////////////////////
  //// DIFFERENTIAL INEQUALITIES

  mc::ODEBNDS_SUNDIALS<I,PM,PV> CO;

  CO.options.DISPLAY   = 1;
#if defined( SAVE_RESULTS )
  CO.options.RESRECORD = 500;
#endif
  CO.options.ATOL      = CO.options.ATOLB      = CO.options.ATOLS  = 1e-10;
  CO.options.RTOL      = CO.options.RTOLB      = CO.options.RTOLS  = 1e-8;
  CO.options.ETOL      = CO.options.ETOLB      = CO.options.ETOLS  = 1e-20;
  CO.options.NMAX      = CO.options.ASACHKPT   = 25000;
  CO.options.MAXFAIL   = 10;
  CO.options.FSAERR    = true;
  CO.options.QERRS     = true;
  CO.options.INTMETH   = mc::ODEBNDS_SUNDIALS<I,PM,PV>::Options::MSBDF;//MSADAMS;//MSBDF;
  CO.options.JACAPPROX = mc::ODEBNDS_SUNDIALS<I,PM,PV>::Options::CV_DIAG;//CV_DENSE;//CV_DIAG;
  CO.options.FSACORR   = mc::ODEBNDS_SUNDIALS<I,PM,PV>::Options::SIMULTANEOUS;
  CO.options.ORDMIT    = -2; //PMp->nord();
  CO.options.WRAPMIT   = mc::ODEBNDS_SUNDIALS<I,PM,PV>::Options::ELLIPS;//NONE;//DINEQ;
  CO.options.QERRB     = true;
  CO.options.QSCALE    = 1e-5;
  CO.options.HMIN      = 1e-12;
  CO.options.ODESLVS   = CO.options;

  //CO.options.ATOL        = CO.options.ODESCO.ATOLB      = CO.ODESCO.options.ATOLS  = 1e-10;
  //CO.options.RTOL        = CO.options.ODESCO.RTOLB      = CO.ODESCO.options.RTOLS  = 1e-8;
  //CO.options.ODESCO.NMAX = CO.options.ODESCO.ASACHKPT = 25000;
  //CO.options.ODESCO.INTMETH   = mc::ODESLVS_SUNDIALS::Options::MSBDF;//MSADAMS;//MSBDF;
  //CO.options.ODESCO.JACAPPROX = mc::ODESLVS_SUNDIALS::Options::CV_DIAG;//CV_DENSE;//CV_DIAG;

  CO.set_dag( &IVP );
  CO.set_time( NS, TS );
  CO.set_state( NX, X );
  CO.set_parameter( NP, P );
  CO.set_differential( NX, RHS );
  CO.set_initial( NX, IC );
  CO.set_quadrature( NQ, QUAD, Q );
  CO.set_function( NS, NF, FCT );


  std::cout << "\nNON_VALIDATED INTEGRATION - INNER-APPROXIMATION OF REACHABLE SET AND ADJOINT SENSITIVITY:\n\n";
  CO.bounds_ASA( NSAMP, Ip );//, Ixk, If, Iyk, Ifp );
#if defined( SAVE_RESULTS )
  ofSTA.open( "test2_APPROX_STA.dat", std::ios_base::out );
  for( unsigned i=0; i<NF; ++i ){
    sprintf( fname, "test2_APPROX_ASA%d.dat",i );  
    ofASA[i].open( fname, std::ios_base::out );
  }
  CO.record( ofSTA, ofASA );
  ofSTA.close();
  for( unsigned i=0; i<NF; ++i ) ofASA[i].close();
#endif

  std::cout << "\nCONTINUOUS SET-VALUED INTEGRATION - INTERVAL ENCLOSURE OF REACHABLE SET AND ADJOINT SENSITIVITY:\n\n";
  CO.bounds_ASA( Ip );//, Ixk, If, Iyk, Ifp );
#if defined( SAVE_RESULTS )
  ofSTA.open( "test2_DINEQI_STA.dat", std::ios_base::out );
  for( unsigned i=0; i<NF; ++i ){
    sprintf( fname, "test2_DINEQI_ASA%d.dat",i );  
    ofASA[i].open( fname, std::ios_base::out );
  }
  CO.record( ofSTA, ofASA );
  ofSTA.close();
  for( unsigned i=0; i<NF; ++i ) ofASA[i].close();
#endif

  std::cout << "\nCONTINUOUS SET-VALUED INTEGRATION - POLYNOMIAL MODEL ENCLOSURE OF REACHABLE SET AND ADJOINT SENSITIVITY:\n\n";
  CO.bounds_ASA( PMp );//, PMxk, PMf, PMyk, PMfp );
#if defined( SAVE_RESULTS )
  ofSTA.open( "test2_DINEQPM_STA.dat", std::ios_base::out );
  for( unsigned i=0; i<NF; ++i ){
    sprintf( fname, "test2_DINEQPM_ASA%d.dat",i );  
    ofASA[i].open( fname, std::ios_base::out );
  }
  CO.record( ofSTA, ofASA );
  ofSTA.close();
  for( unsigned i=0; i<NF; ++i ) ofASA[i].close();
#endif


  std::cout << "\nNON_VALIDATED INTEGRATION - INNER-APPROXIMATION OF REACHABLE SET AND FORWARD SENSITIVITY:\n\n";
  CO.bounds_FSA( NSAMP, Ip );//, Ixk, If, Ixpk, Ifp );
#if defined( SAVE_RESULTS )
  ofSTA.open( "test2_APPROX_STA.dat", std::ios_base::out );
  for( unsigned i=0; i<NP; ++i ){
    sprintf( fname, "test2_APPROX_FSA%d.dat",i );  
    ofFSA[i].open( fname, std::ios_base::out );
  }
  CO.record( ofSTA, ofFSA );
  ofSTA.close();
  for( unsigned i=0; i<NP; ++i ) ofFSA[i].close();
#endif

  std::cout << "\nCONTINUOUS SET-VALUED INTEGRATION - INTERVAL ENCLOSURE OF REACHABLE SET AND FORWARD SENSITIVITY:\n\n";
  CO.bounds_FSA( Ip );//, Ixk, If, Ixpk, Ifp );
#if defined( SAVE_RESULTS )
  ofSTA.open( "test2_DINEQI_STA.dat", std::ios_base::out );
  for( unsigned i=0; i<NP; ++i ){
    sprintf( fname, "test2_DINEQI_FSA%d.dat",i );  
    ofFSA[i].open( fname, std::ios_base::out );
  }
  CO.record( ofSTA, ofFSA );
  ofSTA.close();
  for( unsigned i=0; i<NP; ++i ) ofFSA[i].close();
#endif

  std::cout << "\nCONTINUOUS SET-VALUED INTEGRATION - POLYNOMIAL MODEL ENCLOSURE OF REACHABLE SET AND FORWARD SENSITIVITY:\n\n";
  CO.bounds_FSA( PMp );//, PMxk, PMf, PMxpk, PMfp );
#if defined( SAVE_RESULTS )
  ofSTA.open( "test2_DINEQPM_STA.dat", std::ios_base::out );
  for( unsigned i=0; i<NP; ++i ){
    sprintf( fname, "test2_DINEQPM_FSA%d.dat",i );  
    ofFSA[i].open( fname, std::ios_base::out );
  }
  CO.record( ofSTA, ofFSA );
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
