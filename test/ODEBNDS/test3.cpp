const unsigned int NPM   = 3;	// <- Order of poynomial expansion
const unsigned int NSAMP = 20;	// <- Number of sampling points for inner approx.
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

  double t0 = 0., tf = 1.5;  // Time span
  const unsigned int NS = 100;  // Time stages
  double tk[NS+1]; tk[0] = t0;
  for( unsigned k=0; k<NS; k++ ) tk[k+1] = tk[k] + (tf-t0)/(double)NS;

  const unsigned NP = 1;  // Number of parameters
  const unsigned NX = 4;  // Number of states
  const unsigned NQ = 1;  // Numboer of state quadratures
  const unsigned NF = 2;  // Number of state functions

  mc::FFVar P[NP];  // Parameters
  for( unsigned int i=0; i<NP; i++ ) P[i].set( &IVP );

  mc::FFVar X[NX];  // States
  for( unsigned int i=0; i<NX; i++ ) X[i].set( &IVP );

  mc::FFVar Q[NQ];  // State quadratures
  for( unsigned i=0; i<NQ; i++ ) Q[i].set( &IVP );
  
  const double g       = 9.81e0;		// m/s^2
  const double m1      = 1.e0;		// kg
  const double m2      = 1.e0;		// kg
  const double l1      = 1.e0;		// m
  const double l2      = 1.e0;		// m
  const double psi1_0  = 3.*mc::PI/4.;	// rad
  const double psi2_0  = -11.*mc::PI/20.;	// rad
  const double psi3_0  = .43e0;		// rad/s
  const double psi4_0  = .67e0;		// rad/s

  mc::FFVar RHS[NX];  // Right-hand side function
  mc::FFVar a11 = m1*l1 + m2*(l1+l2*cos(X[1]));
  mc::FFVar a21 = m2*(l1*cos(X[1])+l2);
  mc::FFVar a12 = m2*l2*cos(X[1]);
  mc::FFVar a22 = m2*l2;
  mc::FFVar b1  = -g*(m1+m2)*sin(X[0]) + m2*l2*sin(X[1])*((X[2]+X[3])*(X[2]+X[3]));
  mc::FFVar b2  = -g*m2*sin(X[0]+X[1]) - m2*l1*sin(X[1])*(X[2]*X[2]);
  RHS[0] = X[2];
  RHS[1] = X[3];
  RHS[2] = ( a22*b1 - a12*b2 ) / ( a11*a22 - a12*a21 );
  RHS[3] = ( a11*b2 - a21*b1 ) / ( a11*a22 - a12*a21 );

  mc::FFVar IC[NX];   // Initial value function
  IC[0] = psi1_0 * P[0];
  IC[1] = psi2_0;
  IC[2] = psi3_0;
  IC[3] = psi4_0;

  mc::FFVar QUAD[NQ];  // Quadrature function
  QUAD[0] = X[1];

//  mc::FFVar FCT[NF];  // State functions
//  FCT[0] = X[0] * X[1];
//  FCT[1] = P[0] * X[2] + exp(X[3]);

  mc::FFVar FCT[NF*NS];  // State functions
  for( unsigned k=0; k<NF*NS; k++ ) FCT[k] = 0.;
  //FCT[(NS/2)*NF+0] = X[0];
  FCT[(NS-1)*NF+0] = X[0] * X[1];
  FCT[(NS-1)*NF+1] = P[0] * pow( X[0], 2 );
  for( unsigned k=0; k<NS; k++ ) FCT[k*NF+1] += Q[0];

  I Ip[NP] = { I(0.995,1.005) };
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
  mc::ODESLVS_GSL<I> DP0;

  DP0.set_dag( &IVP );
  DP0.set_state( NX, X );
  DP0.set_parameter( NP, P );
  DP0.set_differential( NX, RHS );
  DP0.set_initial( NX, IC );
  DP0.set_quadrature( NQ, QUAD, Q );
//  DP0.set_function( NF, FCT );
  DP0.set_function( NF, NS, FCT );

  DP0.options.DISPLAY   = 1;
#if defined( SAVE_RESULTS )
  DP0.options.RESRECORD = true;
#endif
  DP0.options.ATOL      = DP0.options.RTOL = 1e-10;
  DP0.options.INTMETH   = mc::ODESLV_GSL<I>::Options::MSBDF;

  // Approximate adjoint bounds
  std::cout << "\nNON_VALIDATED INTEGRATION - APPROXIMATE ENCLOSURE OF REACHABLE SET:\n\n";
  //DP0.bounds( NS, tk, Ip, Ixk, Iq, If, NSAMP );
  DP0.bounds_ASA( NS, tk, Ip, Ixk, 0, If, Iyk, Idf, NSAMP );
#if defined( SAVE_RESULTS )
  std::ofstream apprecSTA("test3_APPROX_STA.dat", std::ios_base::out );
  std::ofstream apprecADJ("test3_APPROX_ADJ.dat", std::ios_base::out );
  DP0.record( apprecSTA, apprecADJ ); 
#endif

  /////////////////////////////////////////////////////////////////////////////
  //// DIFFERENTIAL INEQUALITIES
#ifndef USE_SUNDIALS // GSL integrator
  mc::ODEBNDS_GSL<I,PM,PV> DP;

  DP.options.DISPLAY = 1;
#if defined( SAVE_RESULTS )
  DP.options.RESRECORD = true;
#endif
  DP.options.ATOL      = DP.options.RTOL = 1e-8;
  DP.options.ORDMIT    = 1; //PMp->nord();
  DP.options.WRAPMIT   = mc::ODEBND_GSL<I,PM,PV>::Options::ELLIPS;//DINEQ
  DP.options.HMAX      = 1e-2;
  //DP.options.H0      = 1e-6;

#else // SUNDIALS integrator
  mc::ODEBNDS_SUNDIALS<I,PM,PV> DP;

  DP.options.DISPLAY   = 1;
#if defined( SAVE_RESULTS )
  DP.options.RESRECORD = true;
#endif
  DP.options.NMAX      = 100000;
  DP.options.ASACHKPT  = 100000;
  DP.options.ATOL      = DP.options.ATOLB  = 1e-12;
  DP.options.RTOL      = DP.options.RTOLB  = 1e-8;
  DP.options.INTMETH   = mc::ODEBNDS_SUNDIALS<I,PM,PV>::Options::MSADAMS;
  //DP.options.JACAPPROX = mc::ODEBNDS_SUNDIALS<I,PM,PV>::Options::CV_DENSE;//CV_DIAG;
  DP.options.ORDMIT    = 1; //PMp->nord();
  DP.options.WRAPMIT   = mc::ODEBNDS_SUNDIALS<I,PM,PV>::Options::ELLIPS;//DINEQ;//NONE;
  //DP.options.QSCALE    = 1e0;
  DP.options.HMIN      = 1e-14;
#endif

  DP.set_dag( &IVP );
  DP.set_state( NX, X );
  DP.set_parameter( NP, P );
  DP.set_differential( NX, RHS );
  DP.set_initial( NX, IC );
  DP.set_quadrature( NQ, QUAD, Q );
//  DP.set_function( NF, FCT );
  DP.set_function( NF, NS, FCT );

  std::cout << "\nCONTINUOUS SET-VALUED INTEGRATION - INTERVAL ENCLOSURE OF REACHABLE SET:\n\n";
  DP.bounds_ASA( NS, tk, Ip, Ixk, Iq, If, Iyk, Idf );
#if defined( SAVE_RESULTS )
  std::ofstream direcISTA( "test3_DINEQI_STA.dat", std::ios_base::out );
  std::ofstream direcIADJ( "test3_DINEQI_ADJ.dat", std::ios_base::out );
  DP.record( direcISTA, direcIADJ );
#endif

  std::cout << "\nCONTINUOUS SET-VALUED INTEGRATION - POLYNOMIAL MODEL ENCLOSURE OF REACHABLE SET:\n\n";
  DP.bounds_ASA( NS, tk, PMp, PMxk, PMq, PMf, PMyk, PMdf );
#if defined( SAVE_RESULTS )
  std::ofstream direcPMSTA( "test3_DINEQPM_STA.dat", std::ios_base::out );
  std::ofstream direcPMADJ( "test3_DINEQPM_ADJ.dat", std::ios_base::out );
  DP.record( direcPMSTA, direcPMADJ );
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
