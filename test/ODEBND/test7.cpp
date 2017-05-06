/// 2-D VAN DER POL EQUATIONS
const unsigned int NPM   = 3;	// <- Order of poynomial expansion
const unsigned int NSAMP = 20;	// <- Number of sampling points for inner approx.
#define SAVE_RESULTS		    // <- Whether to save bounds to file
#define USE_CMODEL		        // <- whether to use Chebyshev models or Taylor models

#include "odebnd_sundials.hpp"

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
  mc::FFGraph IVP;           // DAG describing the problem

  const unsigned NP = 2;     // Number of parameters
  const unsigned NX = 2;     // Number of states

  mc::FFVar P[NP];           // Parameter array
  for( unsigned i=0; i<NP; i++ ) P[i].set( &IVP );

  mc::FFVar X[NX];           // State array
  for( unsigned i=0; i<NX; i++ ) X[i].set( &IVP );

  mc::FFVar RHS[NX];         // Dynamics
  RHS[0] = X[1];
  RHS[1] = ( 1 - X[0]*X[0] ) * X[1] - X[0];

  mc::FFVar IC[NX];          // Initial value function
  IC[0] = P[0];
  IC[1] = P[1];

  I Ip[NP];                  // Parameter bounds
  Ip[0] = I( 1.3, 1.4 );
  Ip[1] = I( 2.2, 2.3 );

  PM PMEnv( NP, NPM );
  PV PMp[NP];
  for( unsigned i=0; i<NP; i++ ) PMp[i].set( &PMEnv, i, Ip[i] );

  /////////////////////////////////////////////////////////////////////////
  // ODE trajectories bounding

  mc::ODEBND_SUNDIALS<I,PM,PV> OC;
  OC.set_dag( &IVP );
  OC.set_time( 0., 50. );
  OC.set_state( NX, X );
  OC.set_parameter( NP, P );
  OC.set_differential( NX, RHS );
  OC.set_initial( NX, IC );

#if defined( SAVE_RESULTS )
  OC.options.RESRECORD = 5000;
#endif
  OC.options.INTMETH   = mc::ODEBND_SUNDIALS<I,PM,PV>::Options::MSBDF;
  OC.options.JACAPPROX = mc::ODEBND_SUNDIALS<I,PM,PV>::Options::CV_DIAG;//CV_DENSE;//
  OC.options.WRAPMIT   = mc::ODEBND_SUNDIALS<I,PM,PV>::Options::ELLIPS;//DINEQ;//NONE;
  OC.options.NMAX      = 5000;
  OC.options.DISPLAY   = 1;
  OC.options.ATOL      = 1e-10;
  OC.options.RTOL      = 1e-8;
  OC.options.ETOL      = 1e-20;
  OC.options.HMIN      = 1e-12;
  OC.options.ODESLVS   = OC.options;

  std::cout << "\nNON_VALIDATED INTEGRATION - INNER-APPROXIMATION OF REACHABLE SET:\n\n";
  OC.bounds( NSAMP, Ip );
#if defined( SAVE_RESULTS )
  std::ofstream apprec( "test7_APPROX_STA.dat", std::ios_base::out );
  OC.record( apprec );
#endif

  std::cout << "\nCONTINUOUS SET-VALUED INTEGRATION - INTERVAL ENCLOSURE OF REACHABLE SET:\n\n";
  OC.bounds( Ip );
#if defined( SAVE_RESULTS )
  std::ofstream bnd2recI( "test7_DINEQI_STA.dat", std::ios_base::out );
  OC.record( bnd2recI );
#endif

  std::cout << "\nCONTINUOUS SET-VALUED INTEGRATION - POLYNOMIAL MODEL ENCLOSURE OF REACHABLE SET:\n\n";
  OC.bounds( PMp );
#if defined( SAVE_RESULTS )
  std::ofstream bnd2recPM( "test7_DINEQPM_STA.dat", std::ios_base::out );
  OC.record( bnd2recPM );
#endif

  return 0;
}

