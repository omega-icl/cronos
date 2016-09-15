// ESCAPE26 Case Study 1 - Exothermic Batch Reactor

const unsigned int NPM   = 2;	// <- Order of poynomial expansion
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

  double t0 = 0., tf = 600.;  // Time span
  const unsigned int NS = 8;  // Time stages
  double tk[NS+1]; tk[0] = t0;
  for( unsigned k=0; k<NS; k++ ) tk[k+1] = tk[k] + (tf-t0)/(double)NS;

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
  I *Ixk[NS+1], *Iyk[NS+1], *Ixpk[NS+1];
  for( unsigned k=0; k<=NS; k++ ){
    Ixk[k] = new I[NX];
    Iyk[k] = new I[NF*NX];
    Ixpk[k] = new I[NP*NX];
  }
  I If[NF], Ifp[NF*NP];

  PM PMEnv( NP, NPM );
  PV PMp[NP];
  for( unsigned i=0; i<NP; i++ ) PMp[i] = PV( &PMEnv, i, Ip[i] );
  PV *PMxk[NS+1], *PMyk[NS+1], *PMxpk[NS+1];
  for( unsigned k=0; k<=NS; k++ ){
    PMxk[k] = new PV[NX];
    PMyk[k] = new PV[NF*NX];
    PMxpk[k] = new PV[NP*NX];
  }
  PV PMf[NF], PMfp[NF*NP];

  /////////////////////////////////////////////////////////////////////////////
  //// DIFFERENTIAL INEQUALITIES
  mc::ODEBNDS_SUNDIALS<I,PM,PV> EBR;

  EBR.options.DISPLAY   = 1;
#if defined( SAVE_RESULTS )
  EBR.options.RESRECORD = true;
#endif
  EBR.options.ATOL      = EBR.options.ATOLB   = EBR.options.ATOLS   = 1e-12;
  EBR.options.RTOL      = EBR.options.RTOLB   = EBR.options.RTOLS   = 1e-8;
  EBR.options.INTMETH   = mc::ODEBNDS_SUNDIALS<I,PM,PV>::Options::MSADAMS;
  EBR.options.JACAPPROX = mc::ODEBNDS_SUNDIALS<I,PM,PV>::Options::CV_DIAG;//CV_DENSE;
  EBR.options.ORDMIT    = -2; //PMp->nord();
  EBR.options.WRAPMIT   = mc::ODEBNDS_SUNDIALS<I,PM,PV>::Options::ELLIPS;//DINEQ;//NONE;
  //EBR.options.QSCALE    = 1e0;
  EBR.options.HMIN      = 1e-12;
  EBR.options.NMAX      = 100000;
  EBR.options.ASACHKPT  = 100000;
  EBR.options.QERRB     = true;
  EBR.options.FSAERR    = true;
  EBR.options.QERRS     = true;

  EBR.options.ODESLV.NMAX      = 100000;
  EBR.options.ODESLV.ASACHKPT  = 100000;
  EBR.options.ODESLV.INTMETH   = mc::ODESLVS_SUNDIALS::Options::MSBDF;//MSADAMS;//MSBDF;
  EBR.options.ODESLV.JACAPPROX = mc::ODESLVS_SUNDIALS::Options::CV_DIAG;//CV_DENSE;//CV_DIAG;

  EBR.set_dag( &IVP );
  EBR.set_state( NX, X );
  EBR.set_parameter( NP, P );
  EBR.set_differential( NS, NX, RHS );
  EBR.set_initial( NX, IC );
  EBR.set_function( NF, FCT );

  std::ofstream ofSTA, ofSEN;

  std::cout << "\nNON_VALIDATED INTEGRATION - INNER-APPROXIMATION OF REACHABLE SET AND ADJOINT SENSITIVITY:\n\n";
  EBR.bounds_ASA( NS, tk, Ip, Ixk, If, Iyk, Ifp, NSAMP );
#if defined( SAVE_RESULTS )
  ofSTA.open( "test5b_APPROX_STA.dat", std::ios_base::out );
  ofSEN.open( "test5b_APPROX_ASA.dat", std::ios_base::out );
  EBR.record( ofSTA, ofSEN );
  ofSTA.close();
  ofSEN.close();
#endif
  { int dum; std::cout << "PAUSED--"; std::cin >> dum; }
/*
  std::cout << "\nCONTINUOUS SET-VALUED INTEGRATION - INTERVAL ENCLOSURE OF REACHABLE SET AND ADJOINT SENSITIVITY:\n\n";
  EBR.bounds_ASA( NS, tk, Ip, Ixk, If, Iyk, Ifp );
#if defined( SAVE_RESULTS )
  ofSTA.open( "test5b_DINEQI_STA.dat", std::ios_base::out );
  ofSEN.open( "test5b_DINEQI_ASA.dat", std::ios_base::out );
  EBR.record( ofSTA, ofSEN );
  ofSTA.close();
  ofSEN.close();
#endif
*/
  std::cout << "\nCONTINUOUS SET-VALUED INTEGRATION - POLYNOMIAL MODEL ENCLOSURE OF REACHABLE SET AND ADJOINT SENSITIVITY:\n\n";
  //EBR.bounds( NS, tk, PMp, PMxk, PMf );
  EBR.bounds_ASA( NS, tk, PMp, PMxk, PMf, PMyk, PMfp );
#if defined( SAVE_RESULTS )
  ofSTA.open( "test5b_DINEQPM_STA.dat", std::ios_base::out );
  ofSEN.open( "test5b_DINEQPM_ASA.dat", std::ios_base::out );
  EBR.record( ofSTA, ofSEN );
  ofSTA.close();
  ofSEN.close();
#endif
  { int dum; std::cout << "--PAUSED "; std::cin >> dum; }

  std::cout << "\nNON_VALIDATED INTEGRATION - INNER-APPROXIMATION OF REACHABLE SET AND FORWARD SENSITIVITY:\n\n";
  EBR.bounds_FSA( NS, tk, Ip, Ixk, If, Ixpk, Ifp, NSAMP );
#if defined( SAVE_RESULTS )
  ofSTA.open( "test5b_APPROX_STA.dat", std::ios_base::out );
  ofSEN.open( "test5b_APPROX_FSA.dat", std::ios_base::out );
  EBR.record( ofSTA, ofSEN );
  ofSTA.close();
  ofSEN.close();
#endif
  { int dum; std::cout << "PAUSED--"; std::cin >> dum; }
/*
  std::cout << "\nCONTINUOUS SET-VALUED INTEGRATION - INTERVAL ENCLOSURE OF REACHABLE SET AND FORWARD SENSITIVITY:\n\n";
  //EBR.bounds( NS, tk, Ip, Ixk, If );
  EBR.bounds_FSA( NS, tk, Ip, Ixk, If, Ixpk, Ifp );
#if defined( SAVE_RESULTS )
  ofSTA.open( "test5b_DINEQI_STA.dat", std::ios_base::out );
  ofSEN.open( "test5b_DINEQI_FSA.dat", std::ios_base::out );
  EBR.record( ofSTA, ofSEN );
  ofSTA.close();
  ofSEN.close();
#endif
*/
  std::cout << "\nCONTINUOUS SET-VALUED INTEGRATION - POLYNOMIAL MODEL ENCLOSURE OF REACHABLE SET AND FORWARD SENSITIVITY:\n\n";
  EBR.bounds_FSA( NS, tk, PMp, PMxk, PMf, PMxpk, PMfp );
#if defined( SAVE_RESULTS )
  ofSTA.open( "test5b_DINEQPM_STA.dat", std::ios_base::out );
  ofSEN.open( "test5b_DINEQPM_FSA.dat", std::ios_base::out );
  EBR.record( ofSTA, ofSEN );
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
