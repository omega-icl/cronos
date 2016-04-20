const unsigned int NPM   = 3;	// <- Order of polynomial expansion
const unsigned int NSAMP = 3;	// <- Number of sampling points for inner approx.
#define SAVE_RESULTS		// <- Saving the results to file
#define USE_CMODEL		// <- whether to use Chebyshev models or Taylor models
#define USE_SUNDIALS		// <- whether to use SUNDIALS or GSL integrator

#ifndef USE_SUNDIALS
  #include "odeslv_gsl.hpp"
  #include "odebnd_gsl.hpp"
#else
  #include "odebnd_sundials.hpp"
#endif
#include "odebnd_val.hpp"

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

int main()
{
  //// SETUP PROBLEM AND DAG
  mc::FFGraph IVP;  // DAG describing the problem

  double t0 = 0., tf = 50.;   // Time span
  const unsigned int NS = 100;  // Time stages
  double tk[NS+1]; tk[0] = t0;
  for( unsigned k=0; k<NS; k++ ) tk[k+1] = tk[k] + (tf-t0)/(double)NS;

  const unsigned NP = 6;  // Parameter dimension
  const unsigned NX = 3;  // State dimension
  const unsigned NQ = 1;  // Quadrature dimension
  const unsigned NF = 1;  // Function dimension

  mc::FFVar P[NP];  // Parameter array
  for( unsigned i=0; i<NP; i++ ) P[i].set( &IVP );

  mc::FFVar X[NX];  // State array
  for( unsigned i=0; i<NX; i++ ) X[i].set( &IVP );

  mc::FFVar RHS[NX];  // Right-hand side function
  mc::FFVar &X1 = X[0], &X2 = X[1], &X3 = X[2],
            &X10 = P[0], &X20 = P[1], &X30 = P[2],
            &MUM = P[3], &KS = P[4], &A = P[5];              
  double X3M = 50., D = 0.202, X2F = 20, Y = 0.4, B = 0.2;
  mc::FFVar MU = MUM * ( 1. - X3 / X3M ) * X2 / ( KS + X2 );
  RHS[0] = ( MU - D ) * X1;
  RHS[1] = D * ( X2F - X2 ) - MU * X1 / Y;
  RHS[2] = - D * X3 + ( A * MU + B ) * X1;

  mc::FFVar IC[NX*NS];  // Initial value function
  IC[0] = X10;
  IC[1] = X20;
  IC[2] = X30;

  mc::FFVar Q[NQ];  // State quadratures
  for( unsigned i=0; i<NQ; i++ ) Q[i].set( &IVP );

  mc::FFVar QUAD[NQ];  // Quadrature function
  double X2ref = 3.;
  QUAD[0] = pow( X2-X2ref, 2 );

  mc::FFVar FCT[NF];  // State functions
  FCT[0] = Q[0];

  //// PARAMETER AND STATE BOUNDS
  I Ip[NP] = { I(6.48,6.52), I(4.98,5.02), I(14.98,15.02),
               I(0.46,0.47), I(1.074,1.076), I(2.19,2.21) };
  I *Ixk[NS+1], If[NF];
  for( unsigned k=0; k<=NS; k++ ) Ixk[k] = new I[NX];

  PM PMEnv( NP, NPM );
  PV PMp[NP];
  for( unsigned i=0; i<NP; i++ ) PMp[i].set( &PMEnv, i, Ip[i] );
  PV *PMxk[NS+1], PMf[NF];
  for( unsigned k=0; k<=NS; k++ ) PMxk[k] = new PV[NX];
/*
  /////////////////////////////////////////////////////////////////////////
  // Bound ODE trajectories - sampling
  mc::ODESLV_GSL<I> LV0;

  LV0.set_dag( &IVP );
  LV0.set_state( NX, X );
  LV0.set_parameter( NP, P );
  LV0.set_differential( NX, RHS );
  LV0.set_initial( NX, IC );
  LV0.set_quadrature( NQ, QUAD, Q );
  LV0.set_function( NF, FCT );

  LV0.options.DISPLAY = 1;
  LV0.options.ATOL = LV0.options.RTOL = 1e-10;
  LV0.options.INTMETH = mc::ODESLV_GSL<I>::Options::MSBDF;
#if defined( SAVE_RESULTS )
  LV0.options.RESRECORD = true;
#endif

  LV0.bounds( NS, tk, Ip, Ixk, Iq, If, NSAMP );
#if defined( SAVE_RESULTS )
  std::ofstream apprec( "test2_APPROX_STA.dat", std::ios_base::out );
  LV0.record( apprec );
#endif

  /////////////////////////////////////////////////////////////////////////
  // Bound ODE trajectories - validated integrator
  mc::ODEBND_VAL<I,PM,PV> LV1;

  LV1.set_dag( &IVP );
  LV1.set_state( NX, X );
  LV1.set_parameter( NP, P );
  LV1.set_differential( NX, RHS );
  LV1.set_initial( NX, IC );

  LV1.options.DISPLAY   = 1;
#if defined( SAVE_RESULTS )
  LV1.options.RESRECORD = true;
#endif
  LV1.options.TSTOP     = false;
  LV1.options.TSORDER   = 5;
  LV1.options.PMVALID   = false;
  LV1.options.HSTAB     = false;
  LV1.options.ORDMIT    = 1; //TMp->nord();
  LV1.options.WRAPMIT   = mc::ODEBND_VAL<I,PM,PV>::Options::ELLIPS;//NONE;

  LV1.bounds( NS, tk, PMp, PMxk );
#if defined( SAVE_RESULTS )
  std::ofstream bnd1rec( "test2_outer_val_PM.out", std::ios_base::out );
  LV1.record( bnd1rec );
#endif
*/

  /////////////////////////////////////////////////////////////////////////
  // Bound ODE trajectories - differential inequalities
#ifndef USE_SUNDIALS             // GSL integrator
  mc::ODEBND_GSL<I,PM,PV> LV2;
#else                            // SUNDIALS integrator
  mc::ODEBND_SUNDIALS<I,PM,PV> LV2;
#endif
  LV2.set_dag( &IVP );
  LV2.set_state( NX, X );
  LV2.set_parameter( NP, P );
  LV2.set_differential( NX, RHS );
  LV2.set_initial( NX, IC );
  LV2.set_quadrature( NQ, QUAD, Q );
  LV2.set_function( NF, FCT );

#ifndef USE_SUNDIALS
  LV2.options.WRAPMIT   = mc::ODEBND_GSL<I,PM,PV>::Options::ELLIPS;//DINEQ;//NONE;
#else
  LV2.options.INTMETH   = mc::ODEBND_SUNDIALS<I,PM,PV>::Options::MSADAMS;
  LV2.options.JACAPPROX = mc::ODEBND_SUNDIALS<I,PM,PV>::Options::CV_DIAG;//CV_DENSE;//
  LV2.options.WRAPMIT   = mc::ODEBND_SUNDIALS<I,PM,PV>::Options::ELLIPS;//DINEQ;//NONE;
#endif
  LV2.options.DISPLAY   = 1;
#if defined( SAVE_RESULTS )
  LV2.options.RESRECORD = true;
#endif
  LV2.options.ORDMIT       = -1; //NPM;
  LV2.options.ATOL         = 1e-10;
  LV2.options.RTOL         = 1e-10;
  LV2.options.ODESLV.ATOL  = 1e-10;
  LV2.options.ODESLV.RTOL  = 1e-8;

  std::cout << "\nNON_VALIDATED INTEGRATION - INNER-APPROXIMATION OF REACHABLE SET:\n\n";
  LV2.bounds( NS, tk, Ip, Ixk, If, NSAMP );
#if defined( SAVE_RESULTS )
  std::ofstream apprec( "test2_APPROX_STA.dat", std::ios_base::out );
  LV2.record( apprec );
#endif
  { int dum; std::cout << "PAUSED--"; std::cin >> dum; }

  std::cout << "\nCONTINUOUS SET-VALUED INTEGRATION - INTERVAL ENCLOSURE OF REACHABLE SET:\n\n";
  LV2.bounds( NS, tk, Ip, Ixk, If );
#if defined( SAVE_RESULTS )
  std::ofstream bnd2recI( "test2_DINEQI_STA.dat", std::ios_base::out );
  LV2.record( bnd2recI );
#endif
  { int dum; std::cout << "PAUSED--"; std::cin >> dum; }

  std::cout << "\nCONTINUOUS SET-VALUED INTEGRATION - POLYNOMIAL MODEL ENCLOSURE OF REACHABLE SET:\n\n";
  LV2.options.PMNOREM = false;
  LV2.options.DMAX    = 5.;
  LV2.bounds( NS, tk, PMp, PMxk, PMf );
#if defined( SAVE_RESULTS )
  std::ofstream bnd2recPM( "test2_DINEQPM_STA.dat", std::ios_base::out );
  LV2.record( bnd2recPM );
#endif

  for( unsigned k=0; k<=NS; k++ ){
    delete[] Ixk[k];
    delete[] PMxk[k];
  }

  return 0;
}

