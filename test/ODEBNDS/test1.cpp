const unsigned int NPM   = 4;	// <- Order of poynomial expansion
const unsigned int NSAMP = 50;	// <- Number of sampling points for inner approx.
#define SAVE_RESULTS		// <- Whether to save bounds to file
#define USE_CMODEL		// <- whether to use Chebyshev models or Taylor models
#undef  USE_SUNDIALS		// <- whether to use SUNDIALS or GSL integrator

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
  IC[1] = 1.1 + 0.01*P[0];

  mc::FFVar QUAD[NQ];  // Quadrature function
  QUAD[0] = X[1];

  //mc::FFVar FCT[NF];  // State functions
  //FCT[0] = X[0] * X[1];
  //FCT[1] = P[0] * pow( X[0], 2 );

  mc::FFVar FCT[NF*NS];  // State functions
  for( unsigned k=0; k<NF*NS; k++ ) FCT[k] = 0.;
  //if( NS > 1 ) FCT[(NS/2)*NF+0] = X[0] + 0.1*P[0];
  FCT[(NS-1)*NF+0] = X[0] * X[1];
  FCT[(NS-1)*NF+1] = P[0] * pow( X[0], 2 );
  for( unsigned k=0; k<NS; k++ ) FCT[k*NF+1] += Q[0];

  I Ip[NP] = { I(2.95,3.05) };
  I *Ixk[NS+1], *Iyk[NS+1];
  for( unsigned k=0; k<=NS; k++ ){
    Ixk[k] = new I[NX];
    Iyk[k] = new I[NF*NX];
  }
  I Iq[NQ], If[NF], Idf[NF*NP];

  PM PMEnv( NP, NPM );
  PV PMp[NP] = { PV( &PMEnv, 0, Ip[0] ) };
  PV *PMxk[NS+1], *PMyk[NS+1];
  for( unsigned k=0; k<=NS; k++ ){
    PMxk[k] = new PV[NX];
    PMyk[k] = new PV[NF*NX];
  }
  PV PMq[NQ], PMf[NF], PMdf[NF*NP];

  /////////////////////////////////////////////////////////////////////////////
  //// SAMPLING
  mc::ODESLVS_GSL<I> LV0;

  LV0.set_dag( &IVP );
  LV0.set_state( NX, X );
  LV0.set_parameter( NP, P );
  LV0.set_differential( NX, RHS );
  LV0.set_initial( NX, IC );
  LV0.set_quadrature( NQ, QUAD, Q );
  //LV0.set_function( NF, FCT );
  LV0.set_function( NF, NS, FCT );

  LV0.options.DISPLAY   = 1;
#if defined( SAVE_RESULTS )
  LV0.options.RESRECORD = true;
#endif
  LV0.options.ATOL      = LV0.options.RTOL = 1e-10;
  LV0.options.INTMETH   = mc::ODESLV_GSL<I>::Options::MSBDF;
  //LV0.options.HMAX      = 1e-3;

  // Approximate adjoint bounds
  std::cout << "\nNON_VALIDATED INTEGRATION - APPROXIMATE ENCLOSURE OF REACHABLE SET:\n\n";
  //LV0.bounds( NS, tk, Ip, Ixk, Iq, If, NSAMP );
  LV0.bounds_ASA( NS, tk, Ip, Ixk, 0, If, Iyk, Idf, NSAMP );
  //{ int dum; std::cin >> dum; }
#if defined( SAVE_RESULTS )
  std::ofstream apprecSTA("test1_APPROX_STA.dat", std::ios_base::out );
  std::ofstream apprecADJ("test1_APPROX_ADJ.dat", std::ios_base::out );
  LV0.record( apprecSTA, apprecADJ ); 
#endif

  /////////////////////////////////////////////////////////////////////////////
  //// DIFFERENTIAL INEQUALITIES
#ifndef USE_SUNDIALS // GSL integrator
  mc::ODEBNDS_GSL<I,PM,PV> LV;

  LV.options.DISPLAY = 1;
#if defined( SAVE_RESULTS )
  LV.options.RESRECORD = true;
#endif
  LV.options.ATOL      = LV.options.RTOL = 1e-13;
  LV.options.ORDMIT    = 1; //PMp->nord();
  LV.options.WRAPMIT   = mc::ODEBND_GSL<I,PM,PV>::Options::ELLIPS;//NONE;//DINEQ
  LV.options.HMAX      = 1e-3;
  //LV.options.H0      = 1e-6;
  LV.options.QSCALE    = false;

#else // SUNDIALS integrator
  mc::ODEBNDS_SUNDIALS<I,PM,PV> LV;

  LV.options.DISPLAY   = 1;
#if defined( SAVE_RESULTS )
  LV.options.RESRECORD = true;
#endif
  LV.options.ATOL      = LV.options.ATOLB  = 1e-10;
  LV.options.RTOL      = LV.options.RTOLB  = 1e-8;
  LV.options.NMAX      = 10000;
  LV.options.INTMETH   = mc::ODEBNDS_SUNDIALS<I,PM,PV>::Options::MSADAMS;
  LV.options.JACAPPROX = mc::ODEBNDS_SUNDIALS<I,PM,PV>::Options::CV_DENSE;//CV_DIAG;
  LV.options.ORDMIT    = 1;//PMp->nord();
  LV.options.WRAPMIT   = mc::ODEBNDS_SUNDIALS<I,PM,PV>::Options::ELLIPS;//NONE;//DINEQ;
  LV.options.QERRB     = true;
  LV.options.QSCALE    = false;
#endif

  LV.set_dag( &IVP );
  LV.set_state( NX, X );
  LV.set_parameter( NP, P );
  LV.set_differential( NX, RHS );
  LV.set_initial( NX, IC );
  LV.set_quadrature( NQ, QUAD, Q );
  //LV.set_function( NF, FCT );
  LV.set_function( NF, NS, FCT );

  std::cout << "\nCONTINUOUS SET-VALUED INTEGRATION - INTERVAL ENCLOSURE OF REACHABLE SET:\n\n";
  LV.bounds_ASA( NS, tk, Ip, Ixk, Iq, If, Iyk, Idf );
#if defined( SAVE_RESULTS )
  std::ofstream direcISTA( "test1_DINEQI_STA.dat", std::ios_base::out );
  std::ofstream direcIADJ( "test1_DINEQI_ADJ.dat", std::ios_base::out );
  LV.record( direcISTA, direcIADJ );
#endif

  std::cout << "\nCONTINUOUS SET-VALUED INTEGRATION - POLYNOMIAL MODEL ENCLOSURE OF REACHABLE SET:\n\n";
  LV.bounds_ASA( NS, tk, PMp, PMxk, PMq, PMf, PMyk, PMdf );
#if defined( SAVE_RESULTS )
  std::ofstream direcPMSTA( "test1_DINEQPM_STA.dat", std::ios_base::out );
  std::ofstream direcPMADJ( "test1_DINEQPM_ADJ.dat", std::ios_base::out );
  LV.record( direcPMSTA, direcPMADJ );
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
