/// Belousov-Zhabotinskii reaction system PDE ///
const unsigned int NPM   = 4;	// <- Order of poynomial expansion
const unsigned int NSAMP = 20;	// <- Number of sampling points for inner approx.
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

  double T0 = 0., TF = 20.;   // Time span
  const unsigned NS = 1;      // Time stages
  double TS[NS+1]; TS[0] = T0;
  for( unsigned k=0; k<NS; k++ ) TS[k+1] = TS[k] + (TF-T0)/(double)NS;
  
  double V0 = -30., VF = 20.; // Space range
  const unsigned int NM = 41; // Mesh stages -- Keep this as an odd NUmber, please
  double H = (VF-V0)/(double)NM;

  const unsigned NP = 2;      // Number of parameters
  const unsigned NU = 2;      // Number of states
  const unsigned NF = 1;      // Number of state functions
  const unsigned NX = NU*NM;  // Number of discretized states (ODEs)

  mc::FFVar P[NP];  // Parameters
  for( unsigned int i=0; i<NP; i++ ) P[i].set( &IVP );

  mc::FFVar X[NX];  // States
  for( unsigned int i=0; i<NX; i++ ) X[i].set( &IVP );
  
  // Parameters
  double    rc = 10.,     bc = 1.25; 
  mc::FFVar r  = rc*P[0], b  = bc*P[1];
  double    pc = 2e-2;
  
  // Map u -> x
  mc::FFVar U[NU*(NM+2)];
  for(unsigned int i = 0; i < NM; i++ ) { U[i+1] = X[i]; U[NM+3+i] = X[NM+i]; }
  
  // Boundary Conditions
  U[0]    = U[1];    U[NM+1]        = U[NM];          // u_1
  U[NM+2] = U[NM+3]; U[NU*(NM+2)-1] = U[NU*(NM+2)-2]; // u_2
  
  mc::FFVar IC[NX];   // Initial value function
  IC[NM/2] = 0.5; IC[NM/2 + NM] = 0.5;
  for (unsigned int i = 0      ; i < NM/2; i++) { IC[i] = 0.; IC[i+NM] = 1.; }
  for (unsigned int i = NM/2 +1; i < NM  ; i++) { IC[i] = 1.; IC[i+NM] = 0.; }

  mc::FFVar RHS[NX];  // Right-hand side function
  unsigned int j;
  for( unsigned int i = 1; i < NM+1; i++ ){
    j = i+NM+2;
    RHS[i-1] = ( U[i+1] - 2*U[i] + U[i-1] )/( H*H ) + U[i]*( 1. - U[i] - r*U[j] );
    RHS[j-3] = ( U[j+1] - 2*U[j] + U[j-1] )/( H*H ) - b*U[i]*U[j];
  }

  mc::FFVar FCT[2];  // State functions
  FCT[0] = 0.; for(unsigned int i = 1; i < NM+1; i++) { FCT[0] += U[i]/double(NM);      }
  FCT[1] = 0.; for(unsigned int i = 1; i < NM+1; i++) { FCT[1] += U[i+NM+1]/double(NM); }

  I Ip[NP] = { I(1.-pc,1.+pc), I(1.-pc,1.+pc) };
  //I *Ixk[NS+1], *Iyk[NS+1], *Ixpk[NS+1];
  //for( unsigned k=0; k<=NS; k++ ){
  //  Ixk[k] = new I[NX];
  //  Iyk[k] = new I[NF*NX];
  //  Ixpk[k] = new I[NP*NX];
  //}
  //I If[NF], Ifp[NF*NP];

  PM PMEnv( NP, NPM );
  PV PMp[NP] = { PV( &PMEnv, 0, Ip[0] ),  PV( &PMEnv, 1, Ip[1] ) };
  //PV *PMxk[NS+1], *PMyk[NS+1], *PMxpk[NS+1];
  //for( unsigned k=0; k<=NS; k++ ){
  //  PMxk[k] = new PV[NX];
  //  PMyk[k] = new PV[NF*NX];
  //  PMxpk[k] = new PV[NP*NX];
  //}
  //PV PMf[NF], PMfp[NF*NP];

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
  BZ.options.ODESLVS   = BZ.options;

  BZ.set_dag( &IVP );
  BZ.set_time( NS, TS );
  BZ.set_state( NX, X );
  BZ.set_parameter( NP, P );
  BZ.set_differential( NX, RHS );
  BZ.set_initial( NX, IC );
  BZ.set_function( NF, FCT );

  std::ofstream ofSTA, ofFSA[NP], ofASA[NF];
  char fname[50];

  std::cout << "\nNON_VALIDATED INTEGRATION - INNER-APPROXIMATION OF REACHABLE SET AND ADJOINT SENSITIVITY:\n\n";
  BZ.bounds_ASA( NSAMP, Ip );//, Ixk, If, Iyk, Ifp );
#if defined( SAVE_RESULTS )
  ofSTA.open( "test5_APPROX_STA.dat", std::ios_base::out );
  for( unsigned i=0; i<NF; ++i ){
    sprintf( fname, "test5_APPROX_ASA%d.dat",i );  
    ofASA[i].open( fname, std::ios_base::out );
  }
  BZ.record( ofSTA, ofASA );
  ofSTA.close();
  for( unsigned i=0; i<NF; ++i ) ofASA[i].close();
#endif

  std::cout << "\nCONTINUOUS SET-VALUED INTEGRATION - POLYNOMIAL MODEL ENCLOSURE OF REACHABLE SET AND ADJOINT SENSITIVITY:\n\n";
  BZ.bounds_ASA( PMp );//, PMxk, PMf, PMyk, PMfp );
#if defined( SAVE_RESULTS )
  ofSTA.open( "test5_DINEQPM_STA.dat", std::ios_base::out );
  for( unsigned i=0; i<NF; ++i ){
    sprintf( fname, "test5_DINEQPM_ASA%d.dat",i );  
    ofASA[i].open( fname, std::ios_base::out );
  }
  BZ.record( ofSTA, ofASA );
  ofSTA.close();
  for( unsigned i=0; i<NF; ++i ) ofASA[i].close();
#endif


  std::cout << "\nNON_VALIDATED INTEGRATION - INNER-APPROXIMATION OF REACHABLE SET AND FORWARD SENSITIVITY:\n\n";
  BZ.bounds_FSA( NSAMP, Ip );//, Ixk, If, Ixpk, Ifp );
#if defined( SAVE_RESULTS )
  ofSTA.open( "test5_APPROX_STA.dat", std::ios_base::out );
  for( unsigned i=0; i<NP; ++i ){
    sprintf( fname, "test5_APPROX_FSA%d.dat",i );  
    ofFSA[i].open( fname, std::ios_base::out );
  }
  BZ.record( ofSTA, ofFSA );
  ofSTA.close();
  for( unsigned i=0; i<NP; ++i ) ofFSA[i].close();
#endif

  std::cout << "\nCONTINUOUS SET-VALUED INTEGRATION - POLYNOMIAL MODEL ENCLOSURE OF REACHABLE SET AND FORWARD SENSITIVITY:\n\n";
  BZ.bounds_FSA( PMp );//, PMxk, PMf, PMxpk, PMfp );
#if defined( SAVE_RESULTS )
  ofSTA.open( "test5_DINEQPM_STA.dat", std::ios_base::out );
  for( unsigned i=0; i<NP; ++i ){
    sprintf( fname, "test5_DINEQPM_FSA%d.dat",i );  
    ofFSA[i].open( fname, std::ios_base::out );
  }
  BZ.record( ofSTA, ofFSA );
  ofSTA.close();
  for( unsigned i=0; i<NP; ++i ) ofFSA[i].close();
#endif


  // Clean up
  //for( unsigned k=0; k<=NS; k++ ){
  //  delete[] Ixk[k];
  //  delete[] Iyk[k];
  //  delete[] Ixpk[k];
  //  delete[] PMxk[k];
  //  delete[] PMyk[k];
  //  delete[] PMxpk[k];
  //}

  return 0;
}
