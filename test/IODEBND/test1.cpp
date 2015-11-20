///////////////////////////////////////////////////////////////////////////////
// FIRST EXAMPLE IN:
//   Joseph K. Scott · Paul I. Barton
//   Interval bounds on the solutions of semi-explicit index-one DAEs. Part 2: computation
//   Numer. Math. (2013) 125:27–60, DOI 10.1007/s00211-013-0532-x
///////////////////////////////////////////////////////////////////////////////
const unsigned int NPM   = 5;	// <- Order of polynomial model
const unsigned int NSAMP = 10;	// <- Number of sampling points
#define SAVE_RESULTS		// <- Whether to save bounds to file
#undef  TEST_CONVERGENCE	// <- Whether to test Hausdorff convergence of bounds
#define USE_CMODEL		// <- whether to use Chebyshev models or Taylor models
#define USE_SUNDIALS		// <- whether to use SUNDIALS or GSL integrator
///////////////////////////////////////////////////////////////////////////////

#include "odeslv_gsl.hpp"
#ifndef USE_SUNDIALS
  #include "iodebnd_gsl.hpp"
#else
  #include "iodebnd_sundials.hpp"
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

///////////////////////////////////////////////////////////////////////////////
int main()
{
  mc::FFGraph IVP;  // DAG describing the problem

  double t0 = 0., tf = 1.;   // Time span
  const unsigned int NS = 50;  // Time stages
  double tk[NS+1]; tk[0] = t0;
  for( unsigned k=0; k<NS; k++ ) tk[k+1] = tk[k] + (tf-t0)/(double)NS;

  const unsigned NP = 1;  // Parameter dimension
  const unsigned NX = 2;  // State dimension
  const unsigned ND = 1;  // Differential state dimension

  mc::FFVar P[NP];  // Parameter array
  for( unsigned i=0; i<NP; i++ ) P[i].set( &IVP );

  mc::FFVar X[NX], DX[NX];  // State and state derivative arrays
  for( unsigned i=0; i<NX; i++ ) X[i].set( &IVP );
  for( unsigned i=0; i<NX; i++ ) DX[i].set( &IVP );

  mc::FFVar IC[ND];  // Initial value function
  IC[0] = X[0] - 1.;

  mc::FFVar RHS[ND];  // RHS of differential equations
  RHS[0] = DX[0] + P[0] * X[0] + .1 * X[1];

  mc::FFVar AE[NX-ND];  // Algebraic equations
  AE[0] = X[1] - sin(P[0]) / sqrt(X[1]) - 25. * X[0];

  I Ip[NP] = { I(0.5,4.0) }, Ix0[NX] = { I(1.,1.), I(24.,26.) };
  I *Ixk[NS+1];
  for( unsigned k=0; k<=NS; k++ )
    Ixk[k] = new I[NX];

  PM PMEnv( NP, NPM );
  PV PMp[NP] = { PV( &PMEnv, 0, Ip[0] ) };
  PV *PMxk[NS+1];
  double *Hxk[NS+1];
  for( unsigned k=0; k<=NS; k++ ){
    PMxk[k] = new PV[NX];
    Hxk[k] = new double[NX];
  }
/*
  /////////////////////////////////////////////////////////////////////////
  // Bound ODE trajectories - sampling
  mc::ODESLV_GSL<I> LV0;

  LV0.set_dag( &IVP );
  LV0.set_state( NX, X );
  LV0.set_parameter( NP, P );
  LV0.set_dynamic( RHS );
  LV0.set_initial( IC );

  LV0.options.DISPLAY = 1;
  LV0.options.ATOL    = LV0.options.RTOL = 1e-10;
  LV0.options.H0      = 1e-6;
  LV0.options.INTMETH = mc::ODESLV_GSL<I>::Options::MSBDF;
#if defined( SAVE_RESULTS )
  LV0.options.RESRECORD = true;
#endif

  std::cout << "\nNON_VALIDATED INTEGRATION - APPROXIMATE ENCLOSURE OF REACHABLE SET:\n\n";
  LV0.bounds( NS, tk, Ip, Ixk, 0, NSAMP );
#if defined( SAVE_RESULTS )
  std::ofstream apprec( "approximate.out", std::ios_base::out );
  LV0.record( apprec );
#endif
*/
  /////////////////////////////////////////////////////////////////////////
  // Bound ODE trajectories - differential inequalities
#ifndef USE_SUNDIALS // GSL integrator
  mc::IODEBND_GSL<I,PM,PV> LV2;
#else // SUNDIALS integrator
  mc::IODEBND_SUNDIALS<I,PM,PV> LV2;
#endif
  LV2.set_dag( &IVP );
  LV2.set_state( NX, X, DX );
  LV2.set_parameter( NP, P );
  LV2.set_differential( ND, RHS );
  LV2.set_algebraic( NX-ND, AE );
  LV2.set_initial( ND, IC );
  LV2.set_apriori( Ix0 );

#ifndef USE_SUNDIALS
  LV2.options.WRAPMIT   = mc::IODEBND_GSL<I,PM,PV>::Options::ELLIPS;//DINEQ;//NONE;
#else
  LV2.options.INTMETH   = mc::IODEBND_SUNDIALS<I,PM,PV>::Options::MSADAMS;
  LV2.options.JACAPPROX = mc::IODEBND_SUNDIALS<I,PM,PV>::Options::CV_DIAG;//CV_DENSE;//
  LV2.options.WRAPMIT   = mc::IODEBND_SUNDIALS<I,PM,PV>::Options::ELLIPS;//DINEQ;//NONE;
#endif
  LV2.options.DISPLAY   = 1;
#if defined( SAVE_RESULTS )
  LV2.options.RESRECORD = true;
#endif
  LV2.options.HMIN      = 1e-5;
  LV2.options.ORDMIT    = 5;//NPM; //1;
  LV2.options.ATOL      = LV2.options.RTOL = 1e-8;
  LV2.options.RHSNUMER  = true;
  LV2.options.RHSBNDOPT.DISPLAY = 0;
  LV2.options.RHSBNDOPT.BLKDEC  = false;
  LV2.options.RHSBNDOPT.MAXIT   = 20;
  LV2.options.RHSBNDOPT.RTOL    = 1e-10;
  LV2.options.RHSBNDOPT.BOUNDER = mc::AEBND<I,PM,PV>::Options::GE;
  LV2.options.RHSBNDOPT.PRECOND = mc::AEBND<I,PM,PV>::Options::INVMD;//QRMID;//NONE;
  LV2.options.ICBNDOPT.DISPLAY = 0;
  LV2.options.ICBNDOPT.BLKDEC  = false;
  LV2.options.ICBNDOPT.MAXIT   = 20;
  LV2.options.ICBNDOPT.RTOL    = 1e-10;
  LV2.options.ICBNDOPT.BOUNDER = mc::AEBND<I,PM,PV>::Options::AUTO;
  LV2.options.ICBNDOPT.PRECOND = mc::AEBND<I,PM,PV>::Options::INVMD;//QRMID;//NONE;
/*
  std::cout << "\nCONTINUOUS SET-VALUED INTEGRATION - INTERVAL ENCLOSURE OF REACHABLE SET:\n\n";
  LV2.bounds( NS, tk, Ip, Ixk );
#if defined( SAVE_RESULTS )
  std::ofstream bnd2recI( "test1_DINEQI_STA.dat", std::ios_base::out );
  LV2.record( bnd2recI );
#endif
  //LV2.hausdorff( NS, tk, Ip, Hxk, LV0, NSAMP );
*/
  std::cout << "\nCONTINUOUS SET-VALUED INTEGRATION - POLYNOMIAL MODEL ENCLOSURE OF REACHABLE SET:\n\n";
  LV2.bounds( NS, tk, PMp, PMxk );
#if defined( SAVE_RESULTS )
  std::ofstream bnd2recPM( "test1_DINEQPM_STA.dat", std::ios_base::out );
  LV2.record( bnd2recPM );
#endif
  //LV2.hausdorff( NS, tk, PMp, Hxk, LV0, NSAMP );

#if defined( TEST_CONVERGENCE )
  LV2.options.DISPLAY = 0;
  std::cout << std::scientific << std::setprecision(5);
  for( unsigned int k=0; k<20; k++ ){
    PV PMp0[NP] = { PV( &PMEnv, 0, Ip[0]/pow(1.5,k) ) };
    //LV2.hausdorff( NS, tk, Ip, Hxk, NSAMP );
    LV2.hausdorff( NS, tk, PMp0, Hxk, NSAMP );
    std::cout << mc::Op<I>::diam( PMp0[0].B() ) << "  " << Hxk[NS][0] << std::endl;
  }
#endif

  for( unsigned k=0; k<=NS; k++ ){
    delete[] Ixk[k];
    delete[] PMxk[k];
    delete[] Hxk[k];
  }

  return 0;
}


