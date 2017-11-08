const unsigned int NPM   = 3;	// <- Order of poynomial expansion
const unsigned int NSAMP = 2;	// <- Number of sampling points for inner approx.
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
  mc::FFGraph IVP;  // DAG describing the problem

  const unsigned int NS = 5;  // Time stages
  double TS[NS+1]; TS[0] = 0.;
  for( unsigned k=0; k<NS; k++ ) TS[k+1] = TS[k] + 1./(double)NS;

  const unsigned NP = 2+NS;  // Number of parameters
  const unsigned NX = 2;     // Number of states

  mc::FFVar P[NP];  // Parameter array
  for( unsigned i=0; i<NP; i++ ) P[i].set( &IVP );
  mc::FFVar TF = P[0], X10 = P[1], *U = P+2;

  mc::FFVar X[NX];  // State array
  for( unsigned i=0; i<NX; i++ ) X[i].set( &IVP );

  mc::FFVar RHS[NX*NS];                                            // Dynamics
  for( unsigned k=0; k<NS; k++ ){
    RHS[NX*k+0] = TF * X[1];
    RHS[NX*k+1] = TF * ( U[k]*X[0] - 2.*X[1] );
  }

  mc::FFVar IC[NX];   // Initial value function
  IC[0] = 0.;
  IC[1] = X10;

  I Ip[NP];
  //Ip[0] = I( 2., 6. );
  //Ip[1] = I( -1., 1. );
  //for( unsigned int is=0; is<NS; is++ ) Ip[2+is]  = I( -0.5, 0.5 );
  Ip[0] = I( 2.3828175413204792, 2.4008477663885901 );
  Ip[1] = I( 0.9904726225727921, 1. );
  Ip[2] = I( 1.3653398239103813, 1.5 );
  Ip[3] = I( 1.4238268878303090, 1.5 );
  Ip[4] = I( 1.4157608376380206, 1.5 );
  Ip[5] = I( 0.8211900796950281, 0.9317881555536508 );
  Ip[6] = I( -0.5, -0.476635168577946 );

  PM PMEnv( NP, NPM );
  PV PMp[NP];
  for( unsigned i=0; i<NP; i++ ) PMp[i].set( &PMEnv, i, Ip[i] );
  std::vector<PV*> PMx( NS+1, 0 );

  /////////////////////////////////////////////////////////////////////////
  // ODE trajectories bounding

  mc::ODEBND_SUNDIALS<I,PM,PV> OC;
  OC.set_dag( &IVP );
  OC.set_time( NS, TS );
  OC.set_state( NX, X );
  OC.set_parameter( NP, P );
  OC.set_differential( NS, NX, RHS );
  OC.set_initial( NX, IC );

#if defined( SAVE_RESULTS )
  OC.options.RESRECORD = 100;
#endif
  OC.options.INTMETH   = mc::ODEBND_SUNDIALS<I,PM,PV>::Options::MSADAMS;
  OC.options.JACAPPROX = mc::ODEBND_SUNDIALS<I,PM,PV>::Options::CV_DIAG;//CV_DENSE;//
  OC.options.WRAPMIT   = mc::ODEBND_SUNDIALS<I,PM,PV>::Options::ELLIPS;//DINEQ;//NONE;
  OC.options.NMAX      = 3000;
  OC.options.DISPLAY   = 1;
  OC.options.ATOL      = 1e-10;
  OC.options.RTOL      = 1e-8;
  OC.options.ETOL      = 1e-20;
  OC.options.HMIN      = 1e-10;
  OC.options.ODESLVS   = OC.options;

  std::cout << "\nNON_VALIDATED INTEGRATION - INNER-APPROXIMATION OF REACHABLE SET:\n\n";
  OC.bounds( NSAMP, Ip );
#if defined( SAVE_RESULTS )
  std::ofstream apprec( "test5_APPROX_STA.dat", std::ios_base::out );
  OC.record( apprec );
#endif
/*
  std::cout << "\nCONTINUOUS SET-VALUED INTEGRATION - INTERVAL ENCLOSURE OF REACHABLE SET:\n\n";
  OC.bounds( Ip );
#if defined( SAVE_RESULTS )
  std::ofstream bnd2recI( "test5_DINEQI_STA.dat", std::ios_base::out );
  OC.record( bnd2recI );
#endif
*/
  std::cout << "\nCONTINUOUS SET-VALUED INTEGRATION - POLYNOMIAL MODEL ENCLOSURE OF REACHABLE SET:\n\n";
  OC.bounds( PMp, PMx.data() );

  double p0[NP] = { 2.400848e+00, 1.000000e+00, 1.500000e+00, 1.500000e+00, 1.500000e+00, 8.716292e-01, -5.000000e-01 };
  for( unsigned i=0; i<NX; i++ )
    std::cout << "PMx[" << i << "] @p =" << PMx.back()[i].P(p0) + PMx.back()[i].R() << std::endl;

#if defined( SAVE_RESULTS )
  std::ofstream bnd2recPM( "test5_DINEQPM_STA.dat", std::ios_base::out );
  OC.record( bnd2recPM );
#endif

  return 0;
}

