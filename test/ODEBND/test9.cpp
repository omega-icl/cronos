/// DENBIGH TEST PROBLEM
const unsigned int NPM   = 4;	// <- Order of poynomial expansion
const unsigned int NSAMP = 10;	// <- Number of sampling points for inner approx.
#define SAVE_RESULTS		    // <- Whether to save bounds to file
#define USE_CMODEL		        // <- whether to use Chebyshev models or Taylor models
#define MC__AEBND_SHOW_INTER_FAIL
#define MC__AEBND_IGNORE_INTER_FAIL
#define  USE_PROFIL	    // specify to use PROFIL for interval arithmetic

#include "odebnd_sundials.hpp"
#include "odebnd_expand.hpp"

#ifdef USE_PROFIL
  #include "mcprofil.hpp"
  typedef INTERVAL I;
#else
  #ifdef USE_FILIB
    #include "mcfilib.hpp"
    typedef filib::interval<double> I;
  #else
    #include "interval.hpp"
    typedef mc::Interval I;
  #endif
#endif

#ifdef USE_CMODEL
  #include "cmodel.hpp"
  typedef mc::CModel<I> PM;
  typedef mc::CVar<I> PV;
#else
  #include "tmodel.hpp"
  typedef mc::TModel<I> PM;
  typedef mc::TVar<I> PV;
#endif

//Constant Parameters
const double b1dR = 5.0e3/1.9872;
const double b2dR = 1.0e4/1.9872;
const double a1 = 4.0e3;
const double a2 = 6.2e5;
const double TL = 298.;
const double TU = 398.;
const double tf = 1.;

////////////////////////////////////////////////////////////////////////
// DENBIGH TEST PROBLEM
////////////////////////////////////////////////////////////////////////
int main()
{
  mc::FFGraph DAG;                                                 // DAG

  const unsigned int NS = 1;                                       // Time stages
  double TS[NS+1]; TS[0] = 0.;
  for( unsigned k=0; k<NS; k++ ) TS[k+1] = TS[k] + tf/(double)NS;
  mc::FFVar T; T.set( &DAG );                                      // Time variable

  const unsigned NP = NS;                                          // Parameters
  mc::FFVar P[NP];
  for( unsigned int i=0; i<NP; i++ ) P[i].set( &DAG );

  const unsigned NX = 2;                                           // States
  mc::FFVar X[NX];
  for( unsigned int i=0; i<NX; i++ ) X[i].set( &DAG );

  mc::FFVar RHS[NX*NS];                                            // Dynamics
  for( unsigned k=0; k<NS; k++ ){
    //mc::FFVar K1 = a1 * exp( -b1dR * P[k] / TL );
    //mc::FFVar K2 = a2 * exp( -b2dR * P[k] / TL );
    mc::FFVar K1 = a1 * exp( -b1dR / P[k] );
    mc::FFVar K2 = a2 * exp( -b2dR / P[k] );
    RHS[NX*k+0] = - K1 * sqr( X[0] );
    RHS[NX*k+1] =   K1 * sqr( X[0] ) - K2 * X[1];
  }

  mc::FFVar IC[NX];                                                // Initial
  IC[0] = 1.;
  IC[1] = 0.;

  I Ip[NP];                  // Parameter bounds
  for( unsigned k=0; k<NP; k++ )
    //Ip[k] = TL / I( TL, TU );
    Ip[k] = I( TL, TU );

  PM PMEnv( NP, NPM );
  PV PMp[NP];
  for( unsigned i=0; i<NP; i++ ) PMp[i].set( &PMEnv, i, Ip[i] );

  /////////////////////////////////////////////////////////////////////////
  // ODE trajectories bounding

  mc::ODEBND_SUNDIALS<I,PM,PV> OC;
  OC.set_dag( &DAG );
  OC.set_time( NS, TS, &T );
  OC.set_parameter( NP, P );
  OC.set_state( NX, X );
  OC.set_differential( NS, NX, RHS );
  OC.set_initial( NX, IC );

#if defined( SAVE_RESULTS )
  OC.options.RESRECORD = 500;
#endif
  OC.options.INTMETH   = mc::ODEBND_SUNDIALS<I,PM,PV>::Options::MSBDF;
  OC.options.JACAPPROX = mc::ODEBND_SUNDIALS<I,PM,PV>::Options::CV_DIAG;//CV_DENSE;//
  OC.options.WRAPMIT   = mc::ODEBND_SUNDIALS<I,PM,PV>::Options::ELLIPS;//DINEQ;//NONE;//ELLIPS;
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
  std::ofstream apprec( "test9_APPROX_STA.dat", std::ios_base::out );
  OC.record( apprec );
#endif

  std::cout << "\nCONTINUOUS SET-VALUED INTEGRATION - INTERVAL ENCLOSURE OF REACHABLE SET:\n\n";
  OC.bounds( Ip );
#if defined( SAVE_RESULTS )
  std::ofstream bnd2recI( "test9_DINEQI_STA.dat", std::ios_base::out );
  OC.record( bnd2recI );
#endif

  std::cout << "\nCONTINUOUS SET-VALUED INTEGRATION - POLYNOMIAL MODEL ENCLOSURE OF REACHABLE SET:\n\n";
  OC.options.PMNOREM = false;
  OC.bounds( PMp );
#if defined( SAVE_RESULTS )
  std::ofstream bnd2recPM( "test9_DINEQPM_STA.dat", std::ios_base::out );
  OC.record( bnd2recPM );
#endif


  /////////////////////////////////////////////////////////////////////////
  // ODE discretized bounds

  mc::ODEBND_EXPAND<I,PM,PV> OC2;
  OC2.set_dag( &DAG );
  OC2.set_time( NS, TS, &T );
  OC2.set_parameter( NP, P );
  OC2.set_state( NX, X );
  OC2.set_differential( NS, NX, RHS );
  OC2.set_initial( NX, IC );

  OC2.options.INTMETH   = mc::ODEBND_EXPAND<I,PM,PV>::Options::METHOD::RK;//TS
  OC2.options.TORD      = 4;
  OC2.options.H0        = 0.02;
  OC2.options.LBLK      =
  OC2.options.DBLK      = 2/NS;
  OC2.options.DISPLAY   = 1;
  OC2.options.RESRECORD = true;
  OC2.options.ODESLVS.RESRECORD = 100;
  OC2.options.AEBND.MAXIT   = 100;
  OC2.options.AEBND.DISPLAY = 1;
  OC2.options.AEBND.RTOL    =
  OC2.options.AEBND.ATOL    = 1e-10;
  OC2.options.AEBND.INTERBND = true;
  OC2.options.AEBND.BOUNDER = mc::AEBND<I,PM,PV>::Options::ALGORITHM::GS;//KRAW;//AUTO;
  OC2.options.AEBND.PRECOND = mc::AEBND<I,PM,PV>::Options::PRECONDITIONING::INVMB;//INVBD;//NONE;
  OC2.options.AEBND.BLKDEC  = mc::AEBND<I,PM,PV>::Options::DECOMPOSITION::RECUR;//DIAG;//NONE;
  OC2.setup();

  std::cout << "\nDISCRETIZED SET-VALUED INTEGRATION - INTERVAL ENCLOSURE OF REACHABLE SET:\n\n";
  I Ix0[NX] = { I(0,1), I(0,1) }, *Ixk[NS+1];
  for( unsigned k=0; k<=NS; k++ ){
    Ixk[k] = new I[NX];
    for( unsigned i=0; i<NX; i++ )
      Ixk[k][i] = Ix0[i];
  }
  OC2.bounds( Ip, Ixk );
  std::ofstream ofileI( "test9_EXPANDI_STA.dat", std::ios_base::out );
  OC2.record( ofileI );

  std::cout << "\nDISCRETIZED SET-VALUED INTEGRATION - POLYNOMIAL MODEL ENCLOSURE OF REACHABLE SET:\n\n";
  PV PMx0[NX] = { I(0,1), I(0,1) }, *PMxk[NS+1];
  for( unsigned k=0; k<=NS; k++ ){
    PMxk[k] = new PV[NX];
    for( unsigned i=0; i<NX; i++ )
      PMxk[k][i] = PMx0[i];
  }
  OC2.bounds( PMp, PMxk );
  std::ofstream ofilePM( "test9_EXPANDPM_STA.dat", std::ios_base::out );
  OC2.record( ofilePM );

  return 0;
}

