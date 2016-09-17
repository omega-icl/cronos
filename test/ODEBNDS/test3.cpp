/// DOUBLE PENDULUM ///
const unsigned int NPM   = 3;	// <- Order of poynomial expansion
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

  double t0 = 0., tf = 1.5;  // Time span
  const unsigned int NS = 150;  // Time stages
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
  
  const double g       = 9.81e0;	// m/s^2
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

  mc::FFVar FCT[NF*NS];  // State functions
  for( unsigned k=0; k<NF*NS; k++ ) FCT[k] = 0.;
  FCT[(NS-1)*NF+0] = X[0];
  for( unsigned k=0; k<NS; k++ ) FCT[k*NF+1] += Q[0];

  I Ip[NP] = { I(0.995,1.005) };
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
  mc::ODEBNDS_SUNDIALS<I,PM,PV> DP;

  DP.options.DISPLAY   = 1;
#if defined( SAVE_RESULTS )
  DP.options.RESRECORD = true;
#endif
  DP.options.ATOL      = DP.options.ATOLB      = DP.options.ATOLS  = 1e-10;
  DP.options.RTOL      = DP.options.RTOLB      = DP.options.RTOLS  = 1e-8;
  DP.options.ETOL      = DP.options.ETOLB      = DP.options.ETOLS  = 1e-20;
  DP.options.NMAX      = DP.options.ASACHKPT   = 25000;
  DP.options.MAXFAIL   = 20;
  DP.options.FSAERR    = true;
  DP.options.QERRS     = true;
  //DP.options.AUTOTOLS  = true;
  DP.options.INTMETH   = mc::ODEBNDS_SUNDIALS<I,PM,PV>::Options::MSBDF;//MSSADAMS;//MSBDF;
  DP.options.JACAPPROX = mc::ODEBNDS_SUNDIALS<I,PM,PV>::Options::CV_DIAG;//CV_DENSE;//CV_DIAG;
  DP.options.ORDMIT    = -2; //PMp->nord();
  DP.options.WRAPMIT   = mc::ODEBNDS_SUNDIALS<I,PM,PV>::Options::ELLIPS;//NONE;//DINEQ;
  DP.options.QERRB     = true;
  DP.options.QSCALE    = 1e-5;
  DP.options.HMIN      = 1e-14;

  DP.options.ODESLV.NMAX = DP.options.ODESLV.ASACHKPT = 25000;
  DP.options.ODESLV.INTMETH   = mc::ODESLVS_SUNDIALS::Options::MSBDF;//MSADAMS;//MSBDF;
  DP.options.ODESLV.JACAPPROX = mc::ODESLVS_SUNDIALS::Options::CV_DIAG;//CV_DENSE;//CV_DIAG;

  DP.set_dag( &IVP );
  DP.set_state( NX, X );
  DP.set_parameter( NP, P );
  DP.set_differential( NX, RHS );
  DP.set_initial( NX, IC );
  DP.set_quadrature( NQ, QUAD, Q );
  DP.set_function( NS, NF, FCT );

  std::ofstream ofSTA, ofSEN;

  std::cout << "\nNON_VALIDATED INTEGRATION - INNER-APPROXIMATION OF REACHABLE SET AND ADJOINT SENSITIVITY:\n\n";
  DP.bounds_ASA( NS, tk, Ip, Ixk, If, Iyk, Ifp, NSAMP );
#if defined( SAVE_RESULTS )
  ofSTA.open( "test3_APPROX_STA.dat", std::ios_base::out );
  ofSEN.open( "test3_APPROX_ASA.dat", std::ios_base::out );
  DP.record( ofSTA, ofSEN );
  ofSTA.close();
  ofSEN.close();
#endif
  { int dum; std::cout << "PAUSED--"; std::cin >> dum; }
/*
  std::cout << "\nCONTINUOUS SET-VALUED INTEGRATION - INTERVAL ENCLOSURE OF REACHABLE SET AND ADJOINT SENSITIVITY:\n\n";
  DP.bounds_ASA( NS, tk, Ip, Ixk, If, Iyk, Ifp );
#if defined( SAVE_RESULTS )
  ofSTA.open( "test3_DINEQI_STA.dat", std::ios_base::out );
  ofSEN.open( "test3_DINEQI_ASA.dat", std::ios_base::out );
  DP.record( ofSTA, ofSEN );
  ofSTA.close();
  ofSEN.close();
#endif
*/
  std::cout << "\nCONTINUOUS SET-VALUED INTEGRATION - POLYNOMIAL MODEL ENCLOSURE OF REACHABLE SET AND ADJOINT SENSITIVITY:\n\n";
  //DP.bounds( NS, tk, PMp, PMxk, PMf );
  DP.bounds_ASA( NS, tk, PMp, PMxk, PMf, PMyk, PMfp );
#if defined( SAVE_RESULTS )
  ofSTA.open( "test3_DINEQPM_STA.dat", std::ios_base::out );
  ofSEN.open( "test3_DINEQPM_ASA.dat", std::ios_base::out );
  DP.record( ofSTA, ofSEN );
  ofSTA.close();
  ofSEN.close();
#endif
  { int dum; std::cout << "--PAUSED "; std::cin >> dum; }

  std::cout << "\nNON_VALIDATED INTEGRATION - INNER-APPROXIMATION OF REACHABLE SET AND FORWARD SENSITIVITY:\n\n";
  DP.bounds_FSA( NS, tk, Ip, Ixk, If, Ixpk, Ifp, NSAMP );
#if defined( SAVE_RESULTS )
  ofSTA.open( "test3_APPROX_STA.dat", std::ios_base::out );
  ofSEN.open( "test3_APPROX_FSA.dat", std::ios_base::out );
  DP.record( ofSTA, ofSEN );
  ofSTA.close();
  ofSEN.close();
#endif
  { int dum; std::cout << "PAUSED--"; std::cin >> dum; }
/*
  std::cout << "\nCONTINUOUS SET-VALUED INTEGRATION - INTERVAL ENCLOSURE OF REACHABLE SET AND FORWARD SENSITIVITY:\n\n";
  //DP.bounds( NS, tk, Ip, Ixk, If );
  DP.bounds_FSA( NS, tk, Ip, Ixk, If, Ixpk, Ifp );
#if defined( SAVE_RESULTS )
  ofSTA.open( "test3_DINEQI_STA.dat", std::ios_base::out );
  ofSEN.open( "test3_DINEQI_FSA.dat", std::ios_base::out );
  DP.record( ofSTA, ofSEN );
  ofSTA.close();
  ofSEN.close();
#endif
*/
  std::cout << "\nCONTINUOUS SET-VALUED INTEGRATION - POLYNOMIAL MODEL ENCLOSURE OF REACHABLE SET AND FORWARD SENSITIVITY:\n\n";
  DP.bounds_FSA( NS, tk, PMp, PMxk, PMf, PMxpk, PMfp );
#if defined( SAVE_RESULTS )
  ofSTA.open( "test3_DINEQPM_STA.dat", std::ios_base::out );
  ofSEN.open( "test3_DINEQPM_FSA.dat", std::ios_base::out );
  DP.record( ofSTA, ofSEN );
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

