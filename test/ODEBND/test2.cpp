const unsigned int NPM   = 3;	// <- Order of polynomial expansion
const unsigned int NSAMP = 4;	// <- Number of sampling points for inner approx.
#define SAVE_RESULTS		// <- Saving the results to file
#define USE_CMODEL		// <- whether to use Chebyshev models or Taylor models

#include "odebnd_sundials.hpp"
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

  double T0 = 0., TF = 50.;   // Time span
  const unsigned int NS = 1;  // Time stages
  double TS[NS+1]; TS[0] = T0;
  for( unsigned k=0; k<NS; k++ ) TS[k+1] = TS[k] + (TF-T0)/(double)NS;

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

  PM PMEnv( NP, NPM );
  PV PMp[NP];
  for( unsigned i=0; i<NP; i++ ) PMp[i].set( &PMEnv, i, Ip[i] );
  PV *PMxk[NS+1], PMf[NF];
  for( unsigned k=0; k<=NS; k++ ) PMxk[k] = new PV[NX];

  /////////////////////////////////////////////////////////////////////////
  // Bound ODE trajectories - validated integrator
  mc::ODEBND_VAL<I,PM,PV> LV1;

  LV1.set_dag( &IVP );
  LV1.set_time( NS, TS );
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

  LV1.bounds( PMp, PMxk );
#if defined( SAVE_RESULTS )
  std::ofstream bnd1rec( "test2_VAL_STA.dat", std::ios_base::out );
  LV1.record( bnd1rec );
#endif

  /////////////////////////////////////////////////////////////////////////
  // Bound ODE trajectories - differential inequalities

  mc::ODEBND_SUNDIALS<I,PM,PV> LV2;

  LV2.set_dag( &IVP );
  LV2.set_time( NS, TS );
  LV2.set_state( NX, X );
  LV2.set_parameter( NP, P );
  LV2.set_differential( NX, RHS );
  LV2.set_initial( NX, IC );
  LV2.set_quadrature( NQ, QUAD, Q );
  LV2.set_function( NF, FCT );

#if defined( SAVE_RESULTS )
  LV2.options.RESRECORD = 100;
#endif
  LV2.options.INTMETH   = mc::ODEBND_SUNDIALS<I,PM,PV>::Options::MSADAMS;
  LV2.options.JACAPPROX = mc::ODEBND_SUNDIALS<I,PM,PV>::Options::CV_DIAG;//CV_DENSE;//
  LV2.options.WRAPMIT   = mc::ODEBND_SUNDIALS<I,PM,PV>::Options::ELLIPS;//DINEQ;//NONE;
  LV2.options.DISPLAY   = 1;
  LV2.options.ORDMIT    = -2; //NPM;
  LV2.options.ATOL      = 1e-10;
  LV2.options.RTOL      = 1e-10;
  LV2.options.ETOL      = 1e-20;
  LV2.options.ODESLVS   = LV2.options;

  std::cout << "\nNON_VALIDATED INTEGRATION - INNER-APPROXIMATION OF REACHABLE SET:\n\n";
  LV2.bounds( NSAMP, Ip );
#if defined( SAVE_RESULTS )
  std::ofstream apprec( "test2_APPROX_STA.dat", std::ios_base::out );
  LV2.record( apprec );
#endif

  std::cout << "\nCONTINUOUS SET-VALUED INTEGRATION - INTERVAL ENCLOSURE OF REACHABLE SET:\n\n";
  LV2.bounds( Ip );
#if defined( SAVE_RESULTS )
  std::ofstream bnd2recI( "test2_DINEQI_STA.dat", std::ios_base::out );
  LV2.record( bnd2recI );
#endif


  std::cout << "\nCONTINUOUS SET-VALUED INTEGRATION - POLYNOMIAL MODEL ENCLOSURE OF REACHABLE SET:\n\n";
  LV2.options.PMNOREM = false;
  LV2.options.DMAX    = 5.;
  LV2.bounds( PMp );
#if defined( SAVE_RESULTS )
  std::ofstream bnd2recPM( "test2_DINEQPM_STA.dat", std::ios_base::out );
  LV2.record( bnd2recPM );
#endif

  return 0;
}

