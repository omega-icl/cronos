/// EXOTHERMIC BATCH REACTOR ///
const unsigned int NPM   = 3;	// <- Order of poynomial expansion
const unsigned int NSAMP = 2;	// <- Number of sampling points for inner approx.
#define SAVE_RESULTS		    // <- Whether to save bounds to file
#define USE_CMODEL		        // <- whether to use Chebyshev models or Taylor models

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

  double T0 = 0., TF = 900.;  // Time span
  const unsigned int NS = 8;  // Time stages
  double TS[NS+1]; TS[0] = T0;
  for( unsigned k=0; k<NS; k++ ) TS[k+1] = TS[k] + (TF-T0)/(double)NS;

  const unsigned NP = NS+1;  // Number of parameters
  const unsigned NX = 2;  // Number of states
  const unsigned NF = 1;  // Number of state functions

  mc::FFVar P[NP];  // Parameters
  for( unsigned int i=0; i<NP; i++ ) P[i].set( &IVP );

  mc::FFVar X[NX];  // States
  for( unsigned int i=0; i<NX; i++ ) X[i].set( &IVP );
  mc::FFVar CA = X[0];
  mc::FFVar T  = X[1];

  //Constant Parameters
  const double k0 = 0.022;
  const double CA0 = 10.;
  const double V = 0.1;
  const double Cp = 60.;
  const double Ea = 6000.;
  const double R = 8.314;
  const double DHR = -140000.;
  const double UA = 3.;

  mc::FFVar RHS[NX*NS];  // Right-hand side functions
  for( unsigned i=0; i<NS; i++ ){
    RHS[i*NX+0] = k0*exp(-Ea/(R*T))*(1.-CA);
    RHS[i*NX+1] = (P[1+i] - T)*UA/(CA0*V*Cp) - k0*exp(-Ea/(R*T))*(1.-CA)*DHR/Cp;
  }

  mc::FFVar IC[NX];   // Initial value function
  IC[0] = 0.;
  IC[1] = P[0];

  mc::FFVar FCT[NF];  // State functions
  FCT[0] = CA;

  I Ip[NP] = { I(350.,370.) };
  for( unsigned i=0; i<NS; i++ ) Ip[1+i] = I(290.,310.);
  //I *Ixk[NS+1], *Iyk[NS+1], *Ixpk[NS+1];
  //for( unsigned k=0; k<=NS; k++ ){
  //  Ixk[k] = new I[NX];
  //  Iyk[k] = new I[NF*NX];
  //  Ixpk[k] = new I[NP*NX];
  //}
  //I If[NF], Ifp[NF*NP];

  PM PMEnv( NP, NPM );
  PV PMp[NP];
  for( unsigned i=0; i<NP; i++ ) PMp[i] = PV( &PMEnv, i, Ip[i] );
  //PV *PMxk[NS+1], *PMyk[NS+1], *PMxpk[NS+1];
  //for( unsigned k=0; k<=NS; k++ ){
  //  PMxk[k] = new PV[NX];
  //  PMyk[k] = new PV[NF*NX];
  //  PMxpk[k] = new PV[NP*NX];
  //}
  //PV PMf[NF], PMfp[NF*NP];

  /////////////////////////////////////////////////////////////////////////////
  //// DIFFERENTIAL INEQUALITIES
  mc::ODEBNDS_SUNDIALS<I,PM,PV> EBR;

  EBR.options.DISPLAY   = 1;
#if defined( SAVE_RESULTS )
  EBR.options.RESRECORD = 100;
#endif
  EBR.options.ATOL      = EBR.options.ATOLB      = EBR.options.ATOLS  = 1e-10;
  EBR.options.RTOL      = EBR.options.RTOLB      = EBR.options.RTOLS  = 1e-8;
  EBR.options.ETOL      = EBR.options.ETOLB      = EBR.options.ETOLS  = 1e-20;
  EBR.options.NMAX      = EBR.options.ASACHKPT   = 25000;
  EBR.options.HMIN      = 1e-12;
  EBR.options.MAXFAIL   = 10;
  EBR.options.FSAERR    = true;
  EBR.options.QERRS     = true;
  EBR.options.INTMETH   = mc::ODEBNDS_SUNDIALS<I,PM,PV>::Options::MSADAMS;//MSBDF;
  EBR.options.JACAPPROX = mc::ODEBNDS_SUNDIALS<I,PM,PV>::Options::CV_DIAG;//CV_DENSE;//CV_DIAG;
  EBR.options.FSACORR   = mc::ODEBNDS_SUNDIALS<I,PM,PV>::Options::SIMULTANEOUS;
  EBR.options.ORDMIT    = -2; //PMp->nord();
  EBR.options.WRAPMIT   = mc::ODEBNDS_SUNDIALS<I,PM,PV>::Options::ELLIPS;//NONE;//DINEQ;
  EBR.options.QERRB     = true;
  EBR.options.QSCALE    = 1e-8;
  EBR.options.ODESLVS   = EBR.options;

  EBR.set_dag( &IVP );
  EBR.set_time( NS, TS );
  EBR.set_state( NX, X );
  EBR.set_parameter( NP, P );
  EBR.set_differential( NS, NX, RHS );
  EBR.set_initial( NX, IC );
  EBR.set_function( NF, FCT );


  std::ofstream ofSTA, ofFSA[NP], ofASA[NF];
  char fname[50];

  std::cout << "\nNON_VALIDATED INTEGRATION - INNER-APPROXIMATION OF REACHABLE SET AND ADJOINT SENSITIVITY:\n\n";
  EBR.bounds_ASA( NSAMP, Ip );//, Ixk, If, Iyk, Ifp );
#if defined( SAVE_RESULTS )
  ofSTA.open( "test5_APPROX_STA.dat", std::ios_base::out );
  for( unsigned i=0; i<NF; ++i ){
    sprintf( fname, "test5_APPROX_ASA%d.dat",i );  
    ofASA[i].open( fname, std::ios_base::out );
  }
  EBR.record( ofSTA, ofASA );
  ofSTA.close();
  for( unsigned i=0; i<NF; ++i ) ofASA[i].close();
#endif

  std::cout << "\nCONTINUOUS SET-VALUED INTEGRATION - INTERVAL ENCLOSURE OF REACHABLE SET AND ADJOINT SENSITIVITY:\n\n";
  EBR.bounds_ASA( Ip );//, Ixk, If, Iyk, Ifp );
#if defined( SAVE_RESULTS )
  ofSTA.open( "test5_DINEQI_STA.dat", std::ios_base::out );
  for( unsigned i=0; i<NF; ++i ){
    sprintf( fname, "test5_DINEQI_ASA%d.dat",i );  
    ofASA[i].open( fname, std::ios_base::out );
  }
  EBR.record( ofSTA, ofASA );
  ofSTA.close();
  for( unsigned i=0; i<NF; ++i ) ofASA[i].close();
#endif

  std::cout << "\nCONTINUOUS SET-VALUED INTEGRATION - POLYNOMIAL MODEL ENCLOSURE OF REACHABLE SET AND ADJOINT SENSITIVITY:\n\n";
  EBR.bounds_ASA( PMp );//, PMxk, PMf, PMyk, PMfp );
#if defined( SAVE_RESULTS )
  ofSTA.open( "test5_DINEQPM_STA.dat", std::ios_base::out );
  for( unsigned i=0; i<NF; ++i ){
    sprintf( fname, "test5_DINEQPM_ASA%d.dat",i );  
    ofASA[i].open( fname, std::ios_base::out );
  }
  EBR.record( ofSTA, ofASA );
  ofSTA.close();
  for( unsigned i=0; i<NF; ++i ) ofASA[i].close();
#endif


  std::cout << "\nNON_VALIDATED INTEGRATION - INNER-APPROXIMATION OF REACHABLE SET AND FORWARD SENSITIVITY:\n\n";
  EBR.bounds_FSA( NSAMP, Ip );//, Ixk, If, Ixpk, Ifp );
#if defined( SAVE_RESULTS )
  ofSTA.open( "test5_APPROX_STA.dat", std::ios_base::out );
  for( unsigned i=0; i<NP; ++i ){
    sprintf( fname, "test5_APPROX_FSA%d.dat",i );  
    ofFSA[i].open( fname, std::ios_base::out );
  }
  EBR.record( ofSTA, ofFSA );
  ofSTA.close();
  for( unsigned i=0; i<NP; ++i ) ofFSA[i].close();
#endif

  std::cout << "\nCONTINUOUS SET-VALUED INTEGRATION - INTERVAL ENCLOSURE OF REACHABLE SET AND FORWARD SENSITIVITY:\n\n";
  EBR.bounds_FSA( Ip );//, Ixk, If, Ixpk, Ifp );
#if defined( SAVE_RESULTS )
  ofSTA.open( "test5_DINEQI_STA.dat", std::ios_base::out );
  for( unsigned i=0; i<NP; ++i ){
    sprintf( fname, "test5_DINEQI_FSA%d.dat",i );  
    ofFSA[i].open( fname, std::ios_base::out );
  }
  EBR.record( ofSTA, ofFSA );
  ofSTA.close();
  for( unsigned i=0; i<NP; ++i ) ofFSA[i].close();
#endif

  std::cout << "\nCONTINUOUS SET-VALUED INTEGRATION - POLYNOMIAL MODEL ENCLOSURE OF REACHABLE SET AND FORWARD SENSITIVITY:\n\n";
  EBR.bounds_FSA( PMp );//, PMxk, PMf, PMxpk, PMfp );
#if defined( SAVE_RESULTS )
  ofSTA.open( "test5_DINEQPM_STA.dat", std::ios_base::out );
  for( unsigned i=0; i<NP; ++i ){
    sprintf( fname, "test5_DINEQPM_FSA%d.dat",i );  
    ofFSA[i].open( fname, std::ios_base::out );
  }
  EBR.record( ofSTA, ofFSA );
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
