/// Belousov-Zhabotinskii reaction system PDE ///
const unsigned int NPM   = 5;	// <- Order of poynomial expansion
const unsigned int NSAMP = 2;	// <- Number of sampling points for inner approx.
#define SAVE_RESULTS            // <- Whether to save bounds to file
#define USE_CMODEL              // <- whether to use Chebyshev models or Taylor models

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

#include "odebnds_sundials.hpp"

int main()
{
  mc::FFGraph IVP;  // DAG describing the problem

  double t0 = 0., tf = 20.;   // Time span
  const unsigned int NS = 1; // Time stages
  double tk[NS+1]; tk[0] = t0;
  for( unsigned k=0; k<NS; k++ ) tk[k+1] = tk[k] + (tf-t0)/(double)NS;
  
  double v0 = -30., vf = 20.; // Space range
  const unsigned int NM = 41; // Mesh stages -- Keep this as an odd number, please
  double h = (vf-v0)/(double)NM;

  const unsigned NP = 2;  // Number of parameters
  const unsigned nu = 2;  // Number of states
  const unsigned NF = 1;  // Number of state functions
  
  const unsigned NX = nu*NM; // Number of ODEs

  mc::FFVar P[NP];  // Parameters
  for( unsigned int i=0; i<NP; i++ ) P[i].set( &IVP );

  mc::FFVar X[NX];  // States
  for( unsigned int i=0; i<NX; i++ ) X[i].set( &IVP );
  
  // Parameters
  double    rc = 10.;     double    bc = 1.25; 
  mc::FFVar r  = rc*P[0]; mc::FFVar b  = bc*P[1];
  double pc = 2e-2;
  
  // Map u -> x
  mc::FFVar u[nu*(NM+2)];
  for(unsigned int i = 0; i < NM; i++ ) { u[i+1] = X[i]; u[NM+3+i] = X[NM+i]; }
  
  // Boundary Conditions
  u[0]    = u[1];    u[NM+1]        = u[NM];          // u_1
  u[NM+2] = u[NM+3]; u[nu*(NM+2)-1] = u[nu*(NM+2)-2]; // u_2
  
  mc::FFVar IC[NX];   // Initial value function
  IC[NM/2] = 0.5; IC[NM/2 + NM] = 0.5;
  for (unsigned int i = 0      ; i < NM/2; i++) { IC[i] = 0.; IC[i+NM] = 1.; }
  for (unsigned int i = NM/2 +1; i < NM  ; i++) { IC[i] = 1.; IC[i+NM] = 0.; }

  mc::FFVar RHS[NX];  // Right-hand side function
  unsigned int j;
  for( unsigned int i = 1; i < NM+1; i++ ){
    j = i+NM+2;
    RHS[i-1] = ( u[i+1] - 2*u[i] + u[i-1] )/( h*h ) + u[i]*( 1. - u[i] - r*u[j] );
    RHS[j-3] = ( u[j+1] - 2*u[j] + u[j-1] )/( h*h ) - b*u[i]*u[j];
  }

  mc::FFVar FCT[2];  // State functions
  FCT[0] = 0.; for(unsigned int i = 1; i < NM+1; i++) { FCT[0] += u[i]/double(NM);      }
  FCT[1] = 0.; for(unsigned int i = 1; i < NM+1; i++) { FCT[1] += u[i+NM+1]/double(NM); }

  I Ip[NP] = { I(1.-pc,1.+pc), I(1.-pc,1.+pc) };
  I *Ixk[NS+1], *Iyk[NS+1], *Ixpk[NS+1];
  for( unsigned k=0; k<=NS; k++ ){
    Ixk[k] = new I[NX];
    Iyk[k] = new I[NF*NX];
    Ixpk[k] = new I[NP*NX];
  }
  I If[NF], Ifp[NF*NP];

  PM PMEnv( NP, NPM );
  PV PMp[NP] = { PV( &PMEnv, 0, Ip[0] ),  PV( &PMEnv, 1, Ip[1] ) };
  PV *PMxk[NS+1], *PMyk[NS+1], *PMxpk[NS+1];
  for( unsigned k=0; k<=NS; k++ ){
    PMxk[k] = new PV[NX];
    PMyk[k] = new PV[NF*NX];
    PMxpk[k] = new PV[NP*NX];
  }
  PV PMf[NF], PMfp[NF*NP];

  /////////////////////////////////////////////////////////////////////////////
  //// DIFFERENTIAL INEQUALITIES
  mc::ODEBNDS_SUNDIALS<I,PM,PV> BZ;

  BZ.options.DISPLAY   = 1;
#if defined( SAVE_RESULTS )
  BZ.options.RESRECORD = true;
#endif
  BZ.options.ATOL      = BZ.options.ATOLB      = BZ.options.ATOLS  = 1e-9;
  BZ.options.RTOL      = BZ.options.RTOLB      = BZ.options.RTOLS  = 1e-7;
  BZ.options.ETOL      = BZ.options.ETOLB      = BZ.options.ETOLS  = 1e-20;
  BZ.options.NMAX      = BZ.options.ASACHKPT   = 100000;
  BZ.options.MAXFAIL   = 20;
  BZ.options.FSAERR    = true;
  BZ.options.QERRS     = true;
  BZ.options.INTMETH   = mc::ODEBNDS_SUNDIALS<I,PM,PV>::Options::MSBDF;//MSADAMS;//MSBDF;
  BZ.options.JACAPPROX = mc::ODEBNDS_SUNDIALS<I,PM,PV>::Options::CV_DIAG;//CV_DENSE;
  BZ.options.FSACORR   = mc::ODEBNDS_SUNDIALS<I,PM,PV>::Options::SIMULTANEOUS;
  BZ.options.ORDMIT    = -2; //PMp->nord();
  BZ.options.WRAPMIT   = mc::ODEBNDS_SUNDIALS<I,PM,PV>::Options::ELLIPS;//NONE;//DINEQ;
  BZ.options.QERRB     = true;
  BZ.options.QSCALE    = 1e-5;

  BZ.options.ODESLV.NMAX = 0;
  BZ.options.ODESLV.INTMETH   = mc::ODESLVS_SUNDIALS::Options::MSBDF;//MSADAMS;//MSBDF;
  BZ.options.ODESLV.JACAPPROX = mc::ODESLVS_SUNDIALS::Options::CV_DIAG;//CV_DENSE;

  BZ.set_dag( &IVP );
  BZ.set_state( NX, X );
  BZ.set_parameter( NP, P );
  BZ.set_differential( NX, RHS );
  BZ.set_initial( NX, IC );
  BZ.set_function( NF, FCT );

  std::ofstream ofSTA, ofSEN;

  std::cout << "\nCONTINUOUS SET-VALUED INTEGRATION - POLYNOMIAL MODEL ENCLOSURE OF REACHABLE SET AND ADJOINT SENSITIVITY:\n\n";
  //BZ.bounds( NS, tk, PMp, PMxk, PMf );
  BZ.bounds_ASA( NS, tk, PMp, PMxk, PMf, PMyk, PMfp );
#if defined( SAVE_RESULTS )
  ofSTA.open( "test7_STA_q5_41.dat", std::ios_base::out );
  ofSEN.open( "test7_ASA_q5_41.dat", std::ios_base::out );
  BZ.record( ofSTA, ofSEN );
  ofSTA.close();
  ofSEN.close();
#endif

  std::cout << "\nCONTINUOUS SET-VALUED INTEGRATION - POLYNOMIAL MODEL ENCLOSURE OF REACHABLE SET AND FORWARD SENSITIVITY:\n\n";
  //BZ.bounds( NS, tk, PMp, PMxk, PMf );
  BZ.bounds_FSA( NS, tk, PMp, PMxk, PMf, PMxpk, PMfp );
#if defined( SAVE_RESULTS )
  ofSTA.open( "test7_STA_q5_41.dat", std::ios_base::out );
  ofSEN.open( "test7_FSA_q5_41.dat", std::ios_base::out );
  BZ.record( ofSTA, ofSEN );
  ofSTA.close();
  ofSEN.close();
#endif

  // Clean up
  for( unsigned k=0; k<=NS; k++ ){
    delete[] Ixk[k];
    delete[] Iyk[k];
    delete[] Ixpk[k];
    delete[] PMxk[k];
    delete[] PMyk[k];
    delete[] PMxpk[k];
  }

  return 0;
}
