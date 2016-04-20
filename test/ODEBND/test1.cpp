const unsigned int NPM   = 4;	// <- Order of poynomial expansion
const unsigned int NSAMP = 50;	// <- Number of sampling points for inner approx.
#define SAVE_RESULTS		// <- Whether to save bounds to file
#undef  TEST_CONVERGENCE	// <- Whether to test Hausdorff convergence of bounds
#define USE_CMODEL		// <- whether to use Chebyshev models or Taylor models
#define USE_SUNDIALS		// <- whether to use SUNDIALS or GSL integrator

//#include "odeslv_gsl.hpp"
#ifndef USE_SUNDIALS
  #include "odebnd_gsl.hpp"
#else
  #include "odebnd_sundials.hpp"
#endif
//#include "odebnd_val.hpp"

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

  double t0 = 0., tf = 10.;     // Time span
  const unsigned int NS = 2;  // Time stages
  double tk[NS+1]; tk[0] = t0;
  for( unsigned k=0; k<NS; k++ ) tk[k+1] = tk[k] + (tf-t0)/(double)NS;

  const unsigned NP = 1;  // Number of parameters
  const unsigned NX = 2;  // Number of states
  const unsigned NQ = 1;  // Number of state quadratures
  const unsigned NF = 2;  // Number of state functions

  mc::FFVar P[NP];  // Parameter array
  for( unsigned i=0; i<NP; i++ ) P[i].set( &IVP );

  mc::FFVar X[NX];  // State array
  for( unsigned i=0; i<NX; i++ ) X[i].set( &IVP );

  mc::FFVar Q[NQ];  // State quadratures
  for( unsigned i=0; i<NQ; i++ ) Q[i].set( &IVP );

  mc::FFVar RHS[NX];  // Right-hand side function
  RHS[0] = P[0] * X[0] * ( 1. - X[1] );
  RHS[1] = P[0] * X[1] * ( X[0] - 1. );

  mc::FFVar IC[NX];   // Initial value function
  IC[0] = 1.2;
  IC[1] = 1.1 + 0.01*(P[0]-3.);

  //mc::FFVar IC[NX*NS];  // Initial value and transition functions
  //IC[0] = 1.2+(P[0]-3.);
  //IC[1] = 1.1;
  //for( unsigned k=1; k<NS; k++ )
  //  for( unsigned i=0; i<NX; i++ ) IC[k*NX+i] = X[i];
  //if( NS > 1 ) IC[(NS/2)*NX+0] += 0.5; // <- uncomment to simulate a discontinuity

  mc::FFVar QUAD[NQ];  // Quadrature function
  QUAD[0] = X[1];

  //mc::FFVar FCT[NF];  // State functions
  //FCT[0] = X[0] * X[1];
  //FCT[1] = P[0] * pow( X[0], 2 );

  mc::FFVar FCT[NF*NS];  // State functions
  for( unsigned k=0; k<NF*NS; k++ ) FCT[k] = 0.;
  //if( NS > 1 ) FCT[((NS-1)/NF)*NF+0] = X[0] + 0.1*P[0];
  FCT[(NS-1)*NF+0] = X[0] * X[1];
  FCT[(NS-1)*NF+1] = P[0] * pow( X[0], 2 );
  for( unsigned k=0; k<NS; k++ ) FCT[k*NF+1] += Q[0];

  I Ip[NP] = { I(2.9,3.1) };
  I *Ixk[NS+1];
  for( unsigned k=0; k<=NS; k++ )
    Ixk[k] = new I[NX];
  I If[NF];

  PM PMEnv( NP, NPM );
  PV PMp[NP] = { PV( &PMEnv, 0, Ip[0] ) };
  PV *PMxk[NS+1];
  double *Hxk[NS+1];
  for( unsigned k=0; k<=NS; k++ ){
    PMxk[k] = new PV[NX];
    Hxk[k] = new double[NX];
  }
  PV PMf[NF];
  double Hf[NF];

/*
  /////////////////////////////////////////////////////////////////////////
  // Bound ODE trajectories - sampling
  mc::ODESLV_GSL<I> LV0;

  LV0.set_dag( &IVP );
  LV0.set_state( NX, X );
  LV0.set_parameter( NP, P );
  LV0.set_differential( NX, RHS );
  LV0.set_initial( NS, NX, IC );
  //LV0.set_initial( NX, IC );

  LV0.options.DISPLAY = 1;
  LV0.options.ATOL    = LV0.options.RTOL = 1e-12;
  LV0.options.INTMETH = mc::ODESLV_GSL<I>::Options::MSBDF;
#if defined( SAVE_RESULTS )
  LV0.options.RESRECORD = true;
#endif

  std::cout << "\nNON_VALIDATED INTEGRATION - APPROXIMATE ENCLOSURE OF REACHABLE SET:\n\n";
  LV0.bounds( NS, tk, Ip, Ixk, 0, 0, NSAMP );
#if defined( SAVE_RESULTS )
  std::ofstream apprec( "test1_APPROX_STA.dat", std::ios_base::out );
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
  LV1.options.SCALING   = false;
  LV1.options.TSTOP     = false;
  LV1.options.TSORDER   = 5;
  LV1.options.PMVALID   = false;
  LV1.options.HSTAB     = false;
  LV1.options.TOL       = 1e-10;
  //LV1.options.QTOL      = 1e-12;
  LV1.options.ORDMIT    = 1; //PMp->nord();
  LV1.options.WRAPMIT   = mc::ODEBND_VAL<I,PM,PV>::Options::ELLIPS;//NONE;

  std::cout << "\nDISCRETIZED SET-VALUED INTEGRATION - POLYNOMIAL MODEL ENCLOSURE OF REACHABLE SET:\n\n";
  LV1.bounds( NS, tk, PMp, PMxk );
#if defined( SAVE_RESULTS )
  std::ofstream bnd1rec( "test1_VAL_STA.dat", std::ios_base::out );
      LV1.record( bnd1rec );
#endif

#if defined( TEST_CONVERGENCE )
  LV1.options.DISPLAY = 0;
  std::cout << std::scientific << std::setprecision(5);
  for( unsigned int k=0; k<20; k++ ){
    PV PMp0[NP] = { PV( &PMEnv, 0, Ip[0]/pow(1.5,k) ) };
    LV1.hausdorff( NS, tk, PMp0, Hxk, NSAMP );
    std::cout << mc::Op<I>::diam( PMp0[0].B() ) << "  " << Hxk[NS][0] << std::endl;
  }
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
  //LV2.set_initial( NS, NX, IC );
  LV2.set_quadrature( NQ, QUAD, Q );
  //LV2.set_function( NF, FCT );
  LV2.set_function( NS, NF, FCT );

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
  LV2.options.ATOL         = 1e-20;
  LV2.options.RTOL         = 1e-10;
  LV2.options.ODESLV.ATOL  = 1e-10;
  LV2.options.ODESLV.RTOL  = 1e-8;
/*
  std::cout << "\nNON_VALIDATED INTEGRATION - INNER-APPROXIMATION OF REACHABLE SET:\n\n";
  LV2.bounds( NS, tk, Ip, Ixk, If, NSAMP );
#if defined( SAVE_RESULTS )
  std::ofstream apprec( "test1_APPROX_STA.dat", std::ios_base::out );
  LV2.record( apprec );
#endif
  { int dum; std::cout << "PAUSED--"; std::cin >> dum; }

  std::cout << "\nCONTINUOUS SET-VALUED INTEGRATION - INTERVAL ENCLOSURE OF REACHABLE SET:\n\n";
  LV2.bounds( NS, tk, Ip, Ixk, If );
#if defined( SAVE_RESULTS )
  std::ofstream bnd2recI( "test1_DINEQI_STA.dat", std::ios_base::out );
  LV2.record( bnd2recI );
#endif
  //LV2.hausdorff( NS, tk, Ip, Hxk, LV0, NSAMP );
  { int dum; std::cout << "PAUSED--"; std::cin >> dum; }
*/
  std::cout << "\nCONTINUOUS SET-VALUED INTEGRATION - POLYNOMIAL MODEL ENCLOSURE OF REACHABLE SET:\n\n";
  LV2.options.PMNOREM = false;
  LV2.options.DMAX    = 5.;
  LV2.bounds( NS, tk, PMp, PMxk, PMf );
#if defined( SAVE_RESULTS )
  std::ofstream bnd2recPM( "test1_DINEQPM_STA.dat", std::ios_base::out );
  LV2.record( bnd2recPM );
#endif
  //LV2.hausdorff( NS, tk, PMp, Hxk, Hf, NSAMP );

#if defined( TEST_CONVERGENCE )
  LV2.options.DISPLAY = 0;
  std::cout << std::scientific << std::setprecision(5);
  for( unsigned int k=0; k<10; k++ ){
    I Ipred = mc::Op<I>::mid(Ip[0])+I(-0.5,0.5)*mc::Op<I>::diam(Ip[0])*pow(0.8,k);
    PV PMp0[NP] = { PV( &PMEnv, 0, Ipred ) };
    //LV2.hausdorff( NS, tk, Ip, Hxk, NSAMP );
    LV2.hausdorff( NS, tk, PMp0, Hxk, Hf, NSAMP );
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


