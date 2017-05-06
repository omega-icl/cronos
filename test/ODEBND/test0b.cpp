const unsigned int NPM   = 5;	// <- Order of poynomial expansion
#define USE_SPARSE		    // <- whether to use sparse Chebyshev models

#include <fstream>
#include "odebnd_expand.hpp"

#include "interval.hpp"
typedef mc::Interval I;

#ifdef USE_SPARSE
  #include "scmodel.hpp"
  typedef mc::SCModel<I> PMI;
  typedef mc::SCVar<I>   PVI;
#else
  #include "cmodel.hpp"
  typedef mc::CModel<I> PMI;
  typedef mc::CVar<I>   PVI;
#endif

int main()
{
  mc::FFGraph IVP;  // DAG describing the problem

  const unsigned int NP = 1;  // Parameter dimension
  const unsigned int NX = 2;  // State dimension

  mc::FFVar P[NP];  // Parameter array
  for( unsigned int i=0; i<NP; i++ ) P[i].set( &IVP );

  mc::FFVar X[NX];  // State array
  for( unsigned int i=0; i<NX; i++ ) X[i].set( &IVP );

  mc::FFVar RHS[NX];  // Right-hand side function
  RHS[0] = ( 3. + P[0] ) * X[0] * ( 1. - X[1] );
  RHS[1] = ( 3. + P[0] ) * X[1] * ( X[0] - 1. );

  mc::FFVar IC[NX];   // Initial value function
  IC[0] = 1.2;// + 0.1*P[0];
  IC[1] = 1.1;// + P[0];

  mc::ODEBND_EXPAND<I,PMI,PVI> LV;
  LV.options.H0        = 0.1;
  LV.options.LBLK      =
  LV.options.DBLK      = 10;
  LV.options.DISPLAY   = 1;
  LV.options.RESRECORD = true;
  LV.options.ODESLVS.RESRECORD = 1000;
  LV.options.AEBND.MAXIT   = 100;
  LV.options.AEBND.DISPLAY = 1;
  LV.options.AEBND.RTOL    =
  LV.options.AEBND.ATOL    = 1e-10;
  LV.options.AEBND.BOUNDER = mc::AEBND<I,PMI,PVI>::Options::ALGORITHM::GS;//KRAW;//AUTO;
  LV.options.AEBND.PRECOND = mc::AEBND<I,PMI,PVI>::Options::PRECONDITIONING::INVMB;//INVMB;//INVBD;//NONE;
  LV.options.AEBND.BLKDEC  = mc::AEBND<I,PMI,PVI>::Options::DECOMPOSITION::RECUR;//DIAG;//NONE;

  LV.set_dag( &IVP );
  LV.set_time( 0., 25. );
  LV.set_state( NX, X );
  LV.set_parameter( NP, P );
  LV.set_differential( NX, RHS );
  LV.set_initial( NX, IC );
  LV.setup();


  // Compute Interval Bounds
  I Ip[NP];
  Ip[0] = 0.05 * I(-1.,1);
  I Ix0[NX] = { I(0.,5.), I(0.,5.) }, Ixf[NX] = { I(0.,5.), I(0.,5.) };
  I* Ixk[2] = { Ix0, Ixf };

  LV.bounds( Ip, Ixk );
  std::ofstream ofileI( "test0b_EXPANDI_STA.dat", std::ios_base::out );
  LV.record( ofileI );


  // Compute Polynomial Bounds
#ifdef USE_SPARSE
  PMI PMenv( NPM );   // Polynomial model environment
#else
  PMI PMenv( NP, NPM );   // Polynomial model environment
#endif
  //PMI PMenv( NPM );   // Polynomial model environment
  PVI PMp[NP];
  PMp[0] = PVI( &PMenv, 0, Ip[0] );
  //PVI PMx0[NX] = { I(0.5,1.5), I(0.5,1.5) }, PMxf[NX] = { I(0.5,1.5), I(0.5,1.5) };
  PVI PMx0[NX] = { I(0.,5.), I(0.,5.) }, PMxf[NX] = { I(0.,5.), I(0.,5.) };
  PVI* PMxk[2] = {PMx0, PMxf};

  LV.bounds( PMp, PMxk );
  std::ofstream ofilePM( "test0b_EXPANDPM_STA.dat", std::ios_base::out );
  LV.record( ofilePM );


  // Compute Approximate Bounds (Sampling)
  const unsigned int NSAMP = 50; // Number of sample points
  LV.bounds( NSAMP, Ip );
  std::ofstream ofile0( "test0b_APPROX_STA.dat", std::ios_base::out );
  LV.record( ofile0 );

  // Compute Hausdorff Distance for Enclosures
  //double Hx0[NX], Hxf[NX];
  //double* Hxk[2] = {Hx0, Hxf};

  //LV.hausdorff( 1, tk, Ip, Hxk, 0, NSAMP );
  //LV.hausdorff( 1, tk, TMp, Hxk, 0, NSAMP );

  return 0;
}
