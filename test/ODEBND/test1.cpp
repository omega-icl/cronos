const unsigned int NPM   = 4;	// <- Order of poynomial expansion
const unsigned int NSAMP = 50;	// <- Number of sampling points for inner approx.
#define SAVE_RESULTS		    // <- Whether to save bounds to file
#undef  TEST_CONVERGENCE	    // <- Whether to test Hausdorff convergence of bounds
#define USE_CMODEL		        // <- whether to use Chebyshev models or Taylor models

#include "odebnd_sundials.hpp"
#include "odebnd_val.hpp"

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

  const unsigned int NS = 3;  // Time stages
  double T0 = 0., TF = 15.;   // Time span
  double TS[NS+1]; TS[0] = T0;
  for( unsigned k=0; k<NS; k++ ) TS[k+1] = TS[k] + (TF-T0)/(double)NS;

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
    Ixk[k] = new I[NX+NQ];
  I If[NF];

  PM PMEnv( NP, NPM );
  PV PMp[NP] = { PV( &PMEnv, 0, Ip[0] ) };
  PV *PMxk[NS+1];
  double *Hxk[NS+1];
  for( unsigned k=0; k<=NS; k++ ){
    PMxk[k] = new PV[NX+NQ];
    Hxk[k] = new double[NX+NQ];
  }
  PV PMf[NF];


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
  LV1.bounds( PMp, PMxk );
#if defined( SAVE_RESULTS )
  std::ofstream bnd1rec( "test1_VAL_STA.dat", std::ios_base::out );
      LV1.record( bnd1rec );
#endif

#if defined( TEST_CONVERGENCE )
  LV1.options.DISPLAY = 0;
  std::cout << std::scientific << std::setprecision(5);
  for( unsigned int k=0; k<20; k++ ){
    PV PMp0[NP] = { PV( &PMEnv, 0, Ip[0]/pow(1.5,k) ) };
    LV1.hausdorff( PMp0, Hxk, NSAMP );
    std::cout << mc::Op<I>::diam( PMp0[0].B() ) << "  " << Hxk[NS][0] << std::endl;
  }
#endif


  /////////////////////////////////////////////////////////////////////////
  // ODE trajectories bounding

  mc::ODEBND_SUNDIALS<I,PM,PV> LV2;
  LV2.set_dag( &IVP );
  LV2.set_time( NS, TS );
  LV2.set_state( NX, X );
  LV2.set_parameter( NP, P );
  LV2.set_differential( NX, RHS );
  LV2.set_initial( NX, IC );
  //LV2.set_initial( NS, NX, IC );
  LV2.set_quadrature( NQ, QUAD, Q );
  //LV2.set_function( NF, FCT );
  LV2.set_function( NS, NF, FCT );

  LV2.options.INTMETH   = mc::ODEBND_SUNDIALS<I,PM,PV>::Options::MSADAMS;
  LV2.options.JACAPPROX = mc::ODEBND_SUNDIALS<I,PM,PV>::Options::CV_DIAG;//CV_DENSE;//
  LV2.options.WRAPMIT   = mc::ODEBND_SUNDIALS<I,PM,PV>::Options::ELLIPS;//DINEQ;//NONE;
  LV2.options.NMAX      = 3000;
  LV2.options.DISPLAY   = 1;
#if defined( SAVE_RESULTS )
  LV2.options.RESRECORD = 100;
#endif

  std::cout << "\nNON_VALIDATED INTEGRATION - INNER-APPROXIMATION OF REACHABLE SET:\n\n";
  LV2.bounds( NSAMP, Ip, Ixk, If );
#if defined( SAVE_RESULTS )
  std::ofstream apprec( "test1_APPROX_STA.dat", std::ios_base::out );
  LV2.record( apprec );
#endif

  std::cout << "\nCONTINUOUS SET-VALUED INTEGRATION - INTERVAL ENCLOSURE OF REACHABLE SET:\n\n";
  LV2.bounds( Ip, Ixk, If );
#if defined( SAVE_RESULTS )
  std::ofstream bnd2recI( "test1_DINEQI_STA.dat", std::ios_base::out );
  LV2.record( bnd2recI );
#endif

  std::cout << "\nCONTINUOUS SET-VALUED INTEGRATION - POLYNOMIAL MODEL ENCLOSURE OF REACHABLE SET:\n\n";
  LV2.options.PMNOREM = false;
  LV2.options.DMAX    = 5.;
  LV2.bounds( PMp, PMxk, PMf );
#if defined( SAVE_RESULTS )
  std::ofstream bnd2recPM( "test1_DINEQPM_STA.dat", std::ios_base::out );
  LV2.record( bnd2recPM );
#endif

#if defined( TEST_CONVERGENCE )
  LV2.options.DISPLAY = 0;
  std::cout << std::scientific << std::setprecision(5);
  for( unsigned int k=0; k<10; k++ ){
    I Ipred = mc::Op<I>::mid(Ip[0])+I(-0.5,0.5)*mc::Op<I>::diam(Ip[0])*pow(0.8,k);
    PV PMp0[NP] = { PV( &PMEnv, 0, Ipred ) };
    LV2.hausdorff( NSAMP, PMp0, Hxk );
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


