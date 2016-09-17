/// CUBIC OSCILLATOR ///
const unsigned int NPM   = 3;	// <- Order of poynomial expansion
const unsigned int NSAMP = 2;	// <- Number of sampling points for inner approx.
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

  double t0 = 0., tf = 50.;  // Time span
  const unsigned int NS = 200;  // Time stages
  double tk[NS+1]; tk[0] = t0;
  for( unsigned k=0; k<NS; k++ ) tk[k+1] = tk[k] + (tf-t0)/(double)NS;

  const unsigned NP = 2;  // Number of parameters
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
  const double q = 0.1;
  RHS[0] = q*X[0]*(1.-X[0]*X[0]) + X[1]*(1.-q*X[0]*X[1]);
  RHS[1] = q*X[1]*(1.-X[1]*X[1]) - X[0]*(1.+q*X[0]*X[1]) - 0.2*X[1];

  mc::FFVar IC[NX];   // Initial value function
  IC[0] = P[0];
  IC[1] = P[1];

  mc::FFVar QUAD[NQ];  // Quadrature function
  QUAD[0] = X[1];

  mc::FFVar FCT[NF*NS];  // State functions
  for( unsigned k=0; k<NS-1; k++ ) FCT[k*NF] = 0.;
  FCT[(NS-1)*NF+0] = X[0] * X[1];
  for( unsigned k=0; k<NS; k++ ) FCT[k*NF+1] = Q[0];

  // Parameter and state bounds
  I Ip[NP] = { I(1.75,2.25), I(-0.25,0.25) };
  I *Ixk[NS+1], *Iyk[NS+1], *Ixpk[NS+1];
  for( unsigned k=0; k<=NS; k++ ){
    Ixk[k] = new I[NX];
    Iyk[k] = new I[NF*NX];
    Ixpk[k] = new I[NP*NX];
  }
  I If[NF], Ifp[NF*NP];

  PM PMEnv( NP, NPM );
  PV PMp[NP] = { PV( &PMEnv, 0, Ip[0] ), PV( &PMEnv, 1, Ip[1] )  };
  PV *PMxk[NS+1], *PMyk[NS+1], *PMxpk[NS+1];
  for( unsigned k=0; k<=NS; k++ ){
    PMxk[k] = new PV[NX];
    PMyk[k] = new PV[NF*NX];
    PMxpk[k] = new PV[NP*NX];
  }
  PV PMf[NF], PMfp[NF*NP];

  /////////////////////////////////////////////////////////////////////////////
  //// DIFFERENTIAL INEQUALITIES
  mc::ODEBNDS_SUNDIALS<I,PM,PV> CO;

  CO.options.DISPLAY   = 1;
#if defined( SAVE_RESULTS )
  CO.options.RESRECORD = true;
#endif
  CO.options.ATOL      = CO.options.ATOLB      = CO.options.ATOLS  = 1e-10;
  CO.options.RTOL      = CO.options.RTOLB      = CO.options.RTOLS  = 1e-8;
  CO.options.ETOL      = CO.options.ETOLB      = CO.options.ETOLS  = 1e-20;
  CO.options.NMAX      = CO.options.ASACHKPT   = 25000;
  CO.options.MAXFAIL   = 10;
  CO.options.FSAERR    = true;
  CO.options.QERRS     = true;
  //CO.options.AUTOTOLS  = true;
  CO.options.INTMETH   = mc::ODEBNDS_SUNDIALS<I,PM,PV>::Options::MSBDF;//MSADAMS;//MSBDF;
  CO.options.JACAPPROX = mc::ODEBNDS_SUNDIALS<I,PM,PV>::Options::CV_DIAG;//CV_DENSE;//CV_DIAG;
  CO.options.ORDMIT    = -2; //PMp->nord();
  CO.options.WRAPMIT   = mc::ODEBNDS_SUNDIALS<I,PM,PV>::Options::ELLIPS;//NONE;//DINEQ;
  CO.options.QERRB     = true;
  CO.options.QSCALE    = 1e-5;

  CO.options.ODESLV.NMAX = CO.options.ODESLV.ASACHKPT = 25000;
  CO.options.ODESLV.INTMETH   = mc::ODESLVS_SUNDIALS::Options::MSBDF;//MSADAMS;//MSBDF;
  CO.options.ODESLV.JACAPPROX = mc::ODESLVS_SUNDIALS::Options::CV_DIAG;//CV_DENSE;//CV_DIAG;

  CO.set_dag( &IVP );
  CO.set_state( NX, X );
  CO.set_parameter( NP, P );
  CO.set_differential( NX, RHS );
  CO.set_initial( NX, IC );
  CO.set_quadrature( NQ, QUAD, Q );
  CO.set_function( NS, NF, FCT );

  std::ofstream ofSTA, ofSEN;

  std::cout << "\nNON_VALIDATED INTEGRATION - INNER-APPROXIMATION OF REACHABLE SET AND ADJOINT SENSITIVITY:\n\n";
  CO.bounds_ASA( NS, tk, Ip, Ixk, If, Iyk, Ifp, NSAMP );
#if defined( SAVE_RESULTS )
  ofSTA.open( "test2_APPROX_STA.dat", std::ios_base::out );
  ofSEN.open( "test2_APPROX_ASA.dat", std::ios_base::out );
  CO.record( ofSTA, ofSEN );
  ofSTA.close();
  ofSEN.close();
#endif
  { int dum; std::cout << "PAUSED--"; std::cin >> dum; }
/*
  std::cout << "\nCONTINUOUS SET-VALUED INTEGRATION - INTERVAL ENCLOSURE OF REACHABLE SET AND ADJOINT SENSITIVITY:\n\n";
  CO.bounds_ASA( NS, tk, Ip, Ixk, If, Iyk, Ifp );
#if defined( SAVE_RESULTS )
  ofSTA.open( "test2_DINEQI_STA.dat", std::ios_base::out );
  ofSEN.open( "test2_DINEQI_ASA.dat", std::ios_base::out );
  CO.record( ofSTA, ofSEN );
  ofSTA.close();
  ofSEN.close();
#endif
*/
  std::cout << "\nCONTINUOUS SET-VALUED INTEGRATION - POLYNOMIAL MODEL ENCLOSURE OF REACHABLE SET AND ADJOINT SENSITIVITY:\n\n";
  //CO.bounds( NS, tk, PMp, PMxk, PMf );
  CO.bounds_ASA( NS, tk, PMp, PMxk, PMf, PMyk, PMfp );
#if defined( SAVE_RESULTS )
  ofSTA.open( "test2_DINEQPM_STA.dat", std::ios_base::out );
  ofSEN.open( "test2_DINEQPM_ASA.dat", std::ios_base::out );
  CO.record( ofSTA, ofSEN );
  ofSTA.close();
  ofSEN.close();
#endif
  { int dum; std::cout << "--PAUSED "; std::cin >> dum; }

  std::cout << "\nNON_VALIDATED INTEGRATION - INNER-APPROXIMATION OF REACHABLE SET AND FORWARD SENSITIVITY:\n\n";
  CO.bounds_FSA( NS, tk, Ip, Ixk, If, Ixpk, Ifp, NSAMP );
#if defined( SAVE_RESULTS )
  ofSTA.open( "test2_APPROX_STA.dat", std::ios_base::out );
  ofSEN.open( "test2_APPROX_FSA.dat", std::ios_base::out );
  CO.record( ofSTA, ofSEN );
  ofSTA.close();
  ofSEN.close();
#endif
  { int dum; std::cout << "PAUSED--"; std::cin >> dum; }
/*
  std::cout << "\nCONTINUOUS SET-VALUED INTEGRATION - INTERVAL ENCLOSURE OF REACHABLE SET AND FORWARD SENSITIVITY:\n\n";
  //CO.bounds( NS, tk, Ip, Ixk, If );
  CO.bounds_FSA( NS, tk, Ip, Ixk, If, Ixpk, Ifp );
#if defined( SAVE_RESULTS )
  ofSTA.open( "test2_DINEQI_STA.dat", std::ios_base::out );
  ofSEN.open( "test2_DINEQI_FSA.dat", std::ios_base::out );
  CO.record( ofSTA, ofSEN );
  ofSTA.close();
  ofSEN.close();
#endif
*/
  std::cout << "\nCONTINUOUS SET-VALUED INTEGRATION - POLYNOMIAL MODEL ENCLOSURE OF REACHABLE SET AND FORWARD SENSITIVITY:\n\n";
  CO.bounds_FSA( NS, tk, PMp, PMxk, PMf, PMxpk, PMfp );
#if defined( SAVE_RESULTS )
  ofSTA.open( "test2_DINEQPM_STA.dat", std::ios_base::out );
  ofSEN.open( "test2_DINEQPM_FSA.dat", std::ios_base::out );
  CO.record( ofSTA, ofSEN );
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
