///////////////////////////////////////////////////////////////////////////////
// FIRST EXAMPLE IN:
//   Jai Rayjaguru, Mario E Villanueva, Boris Houska, Benoit Chachuat
//   Continuous-time enclosures for uncertain implicit differential equations
//   IFAC-PapersOnLine (2015) 48(8):94-99, DOI 10.1016/j.ifacol.2015.08.163
///////////////////////////////////////////////////////////////////////////////
const unsigned int NPM   = 7;	// <- Order of polynomial model
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

// Model parameters
const double g       = 9.81e0;	// m/s^2
const double m1      = 1.e0;		// kg
const double m2      = 1.e0;		// kg
const double l1      = 1.e0;		// m
const double l2      = 1.e0;		// m
const double psi1_0  = 3.*mc::PI/4.;	// rad
const double psi2_0  = -11.*mc::PI/20.;	// rad
const double psi3_0  = .43e0;		// rad/s
const double psi4_0  = .67e0;		// rad/s

///////////////////////////////////////////////////////////////////////////////
int main()
{
  mc::FFGraph IVP;  // DAG describing the problem

  double t0 = 0., tf = 5.;   // Time span
  const unsigned int NS = 500;  // Time stages
  double tk[NS+1]; tk[0] = t0;
  for( unsigned k=0; k<NS; k++ ) tk[k+1] = tk[k] + (tf-t0)/(double)NS;

  const unsigned NP = 1;  // Parameter dimension
  const unsigned NX = 4;  // State dimension
  const unsigned ND = 4;  // Differential state dimension

  mc::FFVar P[NP];  // Parameter array
  for( unsigned i=0; i<NP; i++ ) P[i].set( &IVP );

  mc::FFVar X[NX], DX[NX];  // State and state derivative arrays
  for( unsigned i=0; i<NX; i++ ) X[i].set( &IVP );
  for( unsigned i=0; i<NX; i++ ) DX[i].set( &IVP );

  mc::FFVar IC[ND];  // Initial value function
  IC[0] = X[0] - psi1_0 * P[0];
  IC[1] = X[1] - psi2_0;
  IC[2] = X[2] - psi3_0;
  IC[3] = X[3] - psi4_0;

  mc::FFVar RHS[ND];  // RHS of differential equations
  mc::FFVar a11 = m1*l1 + m2*(l1+l2*cos(X[1]));
  mc::FFVar a21 = m2*(l1*cos(X[1])+l2);
  mc::FFVar a12 = m2*l2*cos(X[1]);
  mc::FFVar a22 = m2*l2;
  mc::FFVar b1  = -g*(m1+m2)*sin(X[0]) + m2*l2*sin(X[1])*((X[2]+X[3])*(X[2]+X[3]));
  mc::FFVar b2  = -g*m2*sin(X[0]+X[1]) - m2*l1*sin(X[1])*(X[2]*X[2]);
  RHS[0] = DX[0] - X[2];
  RHS[1] = DX[1] - X[3];
  RHS[2] = DX[2]*a11 + DX[3]*a12 - b1;
  RHS[3] = DX[2]*a21 + DX[3]*a22 - b2;

  I Ip[NP]  = { I(0.99,1.01) },
    Ix0[NX] = { psi1_0*Ip[0], I(psi2_0), I(psi3_0), I(psi4_0) };
  I *Ixk[NS+1];
  for( unsigned k=0; k<=NS; k++ )
    Ixk[k] = new I[NX];

  PM PMEnv( NP, NPM );
  PV PMp[NP];
  for( unsigned i=0; i<NP; i++ ) PMp[i].set( &PMEnv, i, Ip[i] );
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
  LV2.options.HMIN      = 1e-5;//NPM; //1;
  LV2.options.ORDMIT    = 3;//NPM; //1;
  LV2.options.ATOL      = LV2.options.RTOL = 1e-8;
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

  std::cout << "\nCONTINUOUS SET-VALUED INTEGRATION - INTERVAL ENCLOSURE OF REACHABLE SET:\n\n";
  LV2.bounds( NS, tk, Ip, Ixk );
#if defined( SAVE_RESULTS )
  std::ofstream bnd2recI( "test3_DINEQI_STA.dat", std::ios_base::out );
  LV2.record( bnd2recI );
#endif
  //LV2.hausdorff( NS, tk, Ip, Hxk, LV0, NSAMP );

  std::cout << "\nCONTINUOUS SET-VALUED INTEGRATION - POLYNOMIAL MODEL ENCLOSURE OF REACHABLE SET:\n\n";
  LV2.bounds( NS, tk, PMp, PMxk );
#if defined( SAVE_RESULTS )
  std::ofstream bnd2recPM( "test3_DINEQPM_STA.dat", std::ios_base::out );
  LV2.record( bnd2recPM );
#endif
  //LV2.hausdorff( NS, tk, PMp, Hxk, LV0, NSAMP );

#if defined( TEST_CONVERGENCE )
  LV2.options.DISPLAY = 0;
  std::cout << std::scientific << std::setprecision(5);
  for( unsigned int k=0; k<20; k++ ){
    TV TMp0[NP] = { TV( &TMEnv, 0, Ip[0]/pow(1.5,k) ) };
    //LV2.hausdorff( NS, tk, Ip, Hxk, NSAMP );
    LV2.hausdorff( NS, tk, TMp0, Hxk, NSAMP );
    std::cout << mc::Op<I>::diam( TMp0[0].B() ) << "  " << Hxk[NS][0] << std::endl;
  }
#endif

  for( unsigned k=0; k<=NS; k++ ){
    delete[] Ixk[k];
    delete[] PMxk[k];
    delete[] Hxk[k];
  }

  return 0;
}


