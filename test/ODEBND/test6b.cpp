/// TWO-PHASE COUNTER-CURRENT MULTISTAGE LIQUID-LIQUID EXTRACTION WITH SINGLE SOLUTE
const unsigned int NPM   = 3;	// <- Order of poynomial expansion
const unsigned int NSAMP = 2;	// <- Number of sampling points for inner approx.
#define SAVE_RESULTS		    // <- Whether to save bounds to file
#undef  USE_SPARSE	   	        // <- whether to use sparse Chebyshev models

#include "odebnd_expand.hpp"

#include "interval.hpp"
typedef mc::Interval I;
typedef mc::Ellipsoid E;

#ifdef USE_SPARSE
  #include "scmodel.hpp"
  typedef mc::SCModel<I> PM;
  typedef mc::SCVar<I>   PV;
#else
  #include "cmodel.hpp"
  typedef mc::CModel<I> PM;
  typedef mc::CVar<I>   PV;
#endif

int main()
{
  mc::FFGraph IVP;           // DAG describing the problem

  const unsigned NS = 5;     // Number of extraction stages
  const unsigned NP = 6;     // Number of parameters
  const unsigned NX = 2*NS;  // Number of states

  //Constant Parameters
  const double VL = 2.;      // Liquid phase volume [m3]
  const double VG = 2.;      // Gas phase volume [m3]
  const double L  = 5.;      // Liquid phase flowrate [m3/h]
  const double G  = 5.;      // Gas phase flowrate [m3/h]
  const double LF = 10.;     // Solute concentration in liquid phase feed [kg/m3]
  const double GF = 1.;      // Solute concentration in gas phase feed [kg/m3]

  mc::FFVar P[NP];           // Parameter array
  for( unsigned i=0; i<NP; i++ ) P[i].set( &IVP );
  mc::FFVar KLA = P[0];      // Mass transfer constant [1/h]

  mc::FFVar X[NX];           // State array
  for( unsigned i=0; i<NX; i++ ) X[i].set( &IVP );

  mc::FFVar RHS[NX];         // Dynamics
  for( unsigned k=1; k<=NS; k++ ){
    mc::FFVar XEQ = P[5] + X[2*k-1] * ( P[4] + X[2*k-1] * ( P[3] + X[2*k-1] * ( P[2] + X[2*k-1] * P[1] ) ) );
    mc::FFVar Q   = KLA * ( X[2*k-2] - XEQ ) * ( VL + VG );
    RHS[2*k-2] = ( L * ( (k>1?  X[2*k-4]: LF) - X[2*k-2] ) - Q ) / VL;
    RHS[2*k-1] = ( G * ( (k<NS? X[2*k+1]: GF) - X[2*k-1] ) + Q ) / VG;
  }

  mc::FFVar IC[NX];          // Initial value function
  for( unsigned k=1; k<=NS; k++ )
    IC[2*k-2] = IC[2*k-1] = 0.;

  I Ip[NP], Ix0[NX], Ixf[NX]; // Parameter bounds
  Ip[0] = I( 8., 9. );
  Ip[1] = I( 1.48, 1.49 ) * 1e-5;
  Ip[2] = I( -1.11, -1.05 ) * 1e-3;
  Ip[3] = I( 3.28, 3.30 ) * 1e-3;
  Ip[4] = I( 7.56, 7.58 ) * 1e-1;
  Ip[5] = I( 4.93, 4.95 ) * 1e-2;
  for( unsigned i=0; i<NX; i++ ) Ix0[i] = Ixf[i] = I(-0.1,10.);
  I* Ixk[2] = {Ix0, Ixf};

#ifdef USE_SPARSE
  PM PMenv( NPM );
#else
  PM PMenv( NP, NPM );
#endif
  PV PMp[NP], PMx0[NX], PMxf[NX];
  for( unsigned i=0; i<NP; i++ ) PMp[i].set( &PMenv, i, Ip[i] );
  for( unsigned i=0; i<NX; i++ ) PMx0[i] = PMxf[i] = I(-0.1,10.);
  PV* PMxk[2] = {PMx0, PMxf};

  /////////////////////////////////////////////////////////////////////////
  // ODE trajectories bounding

  mc::ODEBND_EXPAND<I,PM,PV> OC;
  OC.set_dag( &IVP );
  OC.set_time( 0., 10. );
  OC.set_state( NX, X );
  OC.set_parameter( NP, P );
  OC.set_differential( NX, RHS );
  OC.set_initial( NX, IC );

#if defined( SAVE_RESULTS )
  OC.options.RESRECORD = true;
  OC.options.ODESLVS.RESRECORD = 1000;
#endif
  OC.options.H0        = 2e-2;
  OC.options.TORD      = 4;
  OC.options.LBLK      = OC.options.DBLK = 4;
  OC.options.DISPLAY   = 1;
  OC.options.AEBND.MAXIT   = 1;
  OC.options.AEBND.DISPLAY = 1;
  OC.options.AEBND.RTOL    = OC.options.AEBND.ATOL = 1e-10;
  OC.options.AEBND.BOUNDER = mc::AEBND<I,PM,PV>::Options::ALGORITHM::GS;//KRAW;//AUTO;
  OC.options.AEBND.PRECOND = mc::AEBND<I,PM,PV>::Options::PRECONDITIONING::INVMB;//INVMB;//INVBD;//NONE;
  OC.options.AEBND.BLKDEC  = mc::AEBND<I,PM,PV>::Options::DECOMPOSITION::RECUR;//DIAG;//NONE;
  OC.setup();

  std::cout << "\nNON_VALIDATED INTEGRATION - INNER-APPROXIMATION OF REACHABLE SET:\n\n";
  OC.bounds( NSAMP, Ip );
#if defined( SAVE_RESULTS )
  std::ofstream apprec( "test6_APPROX_STA.dat", std::ios_base::out );
  OC.record( apprec );
#endif

  std::cout << "\nCONTINUOUS SET-VALUED INTEGRATION - INTERVAL ENCLOSURE OF REACHABLE SET:\n\n";
  OC.bounds( Ip, Ixk );
#if defined( SAVE_RESULTS )
  std::ofstream bnd2recI( "test6_EXPANDI_STA.dat", std::ios_base::out );
  OC.record( bnd2recI );
#endif

  std::cout << "\nCONTINUOUS SET-VALUED INTEGRATION - POLYNOMIAL MODEL ENCLOSURE OF REACHABLE SET:\n\n";
  OC.bounds( PMp, PMxk );
#if defined( SAVE_RESULTS )
  std::ofstream bnd2recPM( "test6_EXPANDPM_STA.dat", std::ios_base::out );
  OC.record( bnd2recPM );
#endif

  return 0;
}

