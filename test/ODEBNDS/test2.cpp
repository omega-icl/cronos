/// CUBIC OSCILLATOR EXAMPLE ///
const unsigned int NPM   = 3;	// <- Order of poynomial expansion
const unsigned int NSAMP = 2;	// <- Number of sampling points for inner approx.
#define SAVE_RESULTS		// <- Whether to save bounds to file
#define USE_CMODEL		// <- whether to use Chebyshev models or Taylor models
#define USE_SUNDIALS		// <- whether to use SUNDIALS or GSL integrator

#include "odeslvs_gsl.hpp"
#ifndef USE_SUNDIALS
  #include "odebnds_gsl.hpp"
#else
  #include "odebnds_sundials.hpp"
#endif

#include "interval.hpp"
typedef mc::Interval I;
typedef mc::Ellipsoid E;

#ifdef USE_CMODEL
  #include "cmodel.hpp"
  typedef mc::CModel<I> PM;
  typedef mc::CVar<I> PV;
#else
  #include "tmodel.hpp"
  typedef mc::TModel<I> PM;
  typedef mc::TVar<I> PV;
#endif

int main()
{
  mc::FFGraph IVP;  // DAG describing the problem

  double t0 = 0., tf = 100.;  // Time span
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

  mc::FFVar FCT[NF];  // State functions
  FCT[0] = X[0] * X[1];
  FCT[1] = P[0] * pow( X[0], 2 );

//  mc::FFVar FCT[NF*NS];  // State functions
//  for( unsigned k=0; k<NF*NS; k++ ) FCT[k] = 0.;
//  //FCT[(NS/2)*NF+0] = X[0];
//  FCT[(NS-1)*NF+0] = X[0] * X[1];
//  FCT[(NS-1)*NF+1] = P[0] * pow( X[0], 2 );
//  for( unsigned k=0; k<NS; k++ ) FCT[k*NF+1] += Q[0];

  I Ip[NP] = { I(1.75,2.25), I(-0.25,0.25) };
  I *Ixk[NS+1], *Iyk[NS+1];
  for( unsigned k=0; k<=NS; k++ ){
    Ixk[k] = new I[NX];
    Iyk[k] = new I[NF*NX];
  }
  
  I Iq[NQ], If[NF], Idf[NF*NP];

  PM PMEnv( NP, NPM );
  PV PMp[NP] = { PV( &PMEnv, 0, Ip[0] ), PV( &PMEnv, 1, Ip[1] )  };
  PV *PMxk[NS+1], *PMyk[NS+1];
  for( unsigned k=0; k<=NS; k++ ){
    PMxk[k] = new PV[NX];
    PMyk[k] = new PV[NF*NX];
  }
  PV PMq[NQ], PMf[NF], PMdf[NF*NP];

  /////////////////////////////////////////////////////////////////////////////
  //// SAMPLING
  mc::ODESLVS_GSL<I> CO0;

  CO0.set_dag( &IVP );
  CO0.set_state( NX, X );
  CO0.set_parameter( NP, P );
  CO0.set_differential( NX, RHS );
  CO0.set_initial( NX, IC );
//  CO0.set_quadrature( NQ, QUAD, Q );
  CO0.set_function( NF, FCT );
//  CO0.set_function( NF, NS, FCT );

  CO0.options.DISPLAY   = 1;
#if defined( SAVE_RESULTS )
  CO0.options.RESRECORD = true;
#endif
  CO0.options.ATOL      = CO0.options.RTOL = 1e-10;
  CO0.options.INTMETH   = mc::ODESLV_GSL<I>::Options::MSBDF;

  // Approximate adjoint bounds
  std::cout << "\nNON_VALIDATED INTEGRATION - APPROXIMATE ENCLOSURE OF REACHABLE SET:\n\n";
  //CO0.bounds( NS, tk, Ip, Ixk, Iq, If, NSAMP );
  CO0.bounds_ASA( NS, tk, Ip, Ixk, 0, If, Iyk, Idf, NSAMP );
#if defined( SAVE_RESULTS )
  std::ofstream apprecSTA("test2_APPROX_STA.dat", std::ios_base::out );
  std::ofstream apprecADJ("test2_APPROX_ADJ.dat", std::ios_base::out );
  CO0.record( apprecSTA, apprecADJ ); 
#endif

  /////////////////////////////////////////////////////////////////////////////
  //// DIFFERENTIAL INEQUALITIES
#ifndef USE_SUNDIALS // GSL integrator
  mc::ODEBNDS_GSL<I,PM,PV> CO;

  CO.options.DISPLAY = 1;
#if defined( SAVE_RESULTS )
  CO.options.RESRECORD = true;
#endif
  CO.options.ATOL      = CO.options.RTOL = 1e-12;
  CO.options.ORDMIT    = 1; //PMp->nord();
  CO.options.WRAPMIT   = mc::ODEBND_GSL<I,PM,PV>::Options::ELLIPS;//ELLIPS//DINEQ
  CO.options.HMAX      = 1e-3;
  //CO.options.H0      = 1e-6;

#else // SUNDIALS integrator
  mc::ODEBNDS_SUNDIALS<I,PM,PV> CO;

  CO.options.DISPLAY   = 1;
#if defined( SAVE_RESULTS )
  CO.options.RESRECORD = true;
#endif
  CO.options.NMAX      = 100000;
  CO.options.ASACHKPT  = 100000;
  CO.options.ATOL      = CO.options.ATOLB  = 1e-10;
  CO.options.RTOL      = CO.options.RTOLB  = 1e-8;
  CO.options.INTMETH   = mc::ODEBNDS_SUNDIALS<I,PM,PV>::Options::MSADAMS;
  CO.options.JACAPPROX = mc::ODEBNDS_SUNDIALS<I,PM,PV>::Options::CV_DENSE;//CV_DIAG;
  CO.options.ORDMIT    = 1; //PMp->nord();
  CO.options.WRAPMIT   = mc::ODEBNDS_SUNDIALS<I,PM,PV>::Options::ELLIPS;//DINEQ;//NONE;
  CO.options.QSCALE    = 1e0;
  CO.options.HMIN      = 1e-12;
#endif

  CO.set_dag( &IVP );
  CO.set_state( NX, X );
  CO.set_parameter( NP, P );
  CO.set_differential( NX, RHS );
  CO.set_initial( NX, IC );
//  CO.set_quadrature( NQ, QUAD, Q );
  CO.set_function( NF, FCT );
//  CO.set_function( NF, NS, FCT );

  std::cout << "\nCONTINUOUS SET-VALUED INTEGRATION - INTERVAL ENCLOSURE OF REACHABLE SET:\n\n";
  CO.bounds_ASA( NS, tk, Ip, Ixk, Iq, If, Iyk, Idf );
#if defined( SAVE_RESULTS )
  std::ofstream direcISTA( "test2_DINEQI_STA.dat", std::ios_base::out );
  std::ofstream direcIADJ( "test2_DINEQI_ADJ.dat", std::ios_base::out );
  CO.record( direcISTA, direcIADJ );
#endif

  std::cout << "\nCONTINUOUS SET-VALUED INTEGRATION - POLYNOMIAL MODEL ENCLOSURE OF REACHABLE SET:\n\n";
  CO.bounds_ASA( NS, tk, PMp, PMxk, PMq, PMf, PMyk, PMdf );
#if defined( SAVE_RESULTS )
  std::ofstream direcPMSTA( "test2_DINEQPM_STA.dat", std::ios_base::out );
  std::ofstream direcPMADJ( "test2_DINEQPM_ADJ.dat", std::ios_base::out );
  CO.record( direcPMSTA, direcPMADJ );
#endif

  // Clean up
  for( unsigned k=0; k<=NS; k++ ){
    delete[] Ixk[k];
    delete[] Iyk[k];
    delete[] PMxk[k];
    delete[] PMyk[k];
  }

  return 0;
}
