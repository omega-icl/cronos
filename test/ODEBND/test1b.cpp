const unsigned int NPM   = 5;	// <- Order of poynomial expansion
const unsigned int NSAMP = 50;	// <- Number of sampling points for inner approx.
#define SAVE_RESULTS		    // <- Whether to save bounds to file
#undef  TEST_CONVERGENCE	    // <- Whether to test Hausdorff convergence of bounds
#define USE_SPARSE		    // <- whether to use sparse Chebyshev models

#include "odebnd_expand.hpp"

#include "interval.hpp"
typedef mc::Interval I;
typedef mc::Ellipsoid E;

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

  const unsigned int NS = 1;  // Time stages
  double T0 = 0., TF = 15.;   // Time span
  double TK[NS+1]; TK[0] = T0;
  for( unsigned k=0; k<NS; k++ ) TK[k+1] = TK[k] + (TF-T0)/(double)NS;

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
  for( unsigned k=0; k<=NS; k++ ){
    Ixk[k] = new I[NX];
    for( unsigned i=0; i<NX; i++ )
      Ixk[k][i] = I(0.,5.);
  }
  I If[NF];

#ifdef USE_SPARSE
  PMI PMenv( NPM );
#else
  PMI PMenv( NP, NPM );
#endif
  PVI PMp[NP] = { PVI( &PMenv, 0, Ip[0] ) };
  PVI *PMxk[NS+1];
  double *Hxk[NS+1];
  for( unsigned k=0; k<=NS; k++ ){
    PMxk[k] = new PVI[NX];
    for( unsigned i=0; i<NX; i++ )
      PMxk[k][i] = I(0.,5.);
    Hxk[k] = new double[NX];
  }
  PVI PMf[NF];
  double Hf[NF];

  /////////////////////////////////////////////////////////////////////////
  // ODE trajectories bounding

  mc::ODEBND_EXPAND<I,PMI,PVI> LV;

  LV.options.H0        = 0.1;
  LV.options.LBLK      =
  LV.options.DBLK      = 10;
  LV.options.DISPLAY   = 1;
#if defined( SAVE_RESULTS )
  LV.options.RESRECORD = true;
#endif

  LV.options.AEBND.MAXIT   = 100;
  LV.options.AEBND.DISPLAY = 1;
  LV.options.AEBND.RTOL    =
  LV.options.AEBND.ATOL    = 1e-10;
  LV.options.AEBND.BOUNDER = mc::AEBND<I,PMI,PVI>::Options::ALGORITHM::GS;//KRAW;//AUTO;
  LV.options.AEBND.PRECOND = mc::AEBND<I,PMI,PVI>::Options::PRECONDITIONING::INVMB;//INVMB;//INVBD;//NONE;
  LV.options.AEBND.BLKDEC  = mc::AEBND<I,PMI,PVI>::Options::DECOMPOSITION::RECUR;//DIAG;//NONE;

  LV.set_dag( &IVP );
  LV.set_time( NS, TK );
  LV.set_state( NX, X );
  LV.set_parameter( NP, P );
  LV.set_differential( NX, RHS );
  LV.set_initial( NX, IC );
  //LV.set_initial( NS, NX, IC );
  //LV.set_quadrature( NQ, QUAD, Q );
  //LV.set_function( NF, FCT );
  //LV.set_function( NS, NF, FCT );
  LV.setup();

  std::cout << "\nNON_VALIDATED INTEGRATION - INNER-APPROXIMATION OF REACHABLE SET:\n\n";
  LV.bounds( NSAMP, Ip, Ixk, If );
#if defined( SAVE_RESULTS )
  std::ofstream apprec( "test1_APPROX_STA.dat", std::ios_base::out );
  LV.record( apprec );
#endif
  { int dum; std::cout << "PAUSED--"; std::cin >> dum; }

  std::cout << "\nCONTINUOUS SET-VALUED INTEGRATION - INTERVAL ENCLOSURE OF REACHABLE SET:\n\n";
  LV.bounds( Ip, Ixk, If );
#if defined( SAVE_RESULTS )
  std::ofstream bnd2recI( "test1_DINEQI_STA.dat", std::ios_base::out );
  LV.record( bnd2recI );
#endif
  //LV.hausdorff( Ip, Hxk, LV0, NSAMP );
  { int dum; std::cout << "PAUSED--"; std::cin >> dum; }

  std::cout << "\nCONTINUOUS SET-VALUED INTEGRATION - POLYNOMIAL MODEL ENCLOSURE OF REACHABLE SET:\n\n";
  LV.options.DMAX    = 5.;
  LV.bounds( PMp, PMxk, PMf );
#if defined( SAVE_RESULTS )
  std::ofstream bnd2recPM( "test1_DINEQPM_STA.dat", std::ios_base::out );
  LV.record( bnd2recPM );
#endif
  //LV.hausdorff( PMp, Hxk, Hf, NSAMP );

#if defined( TEST_CONVERGENCE )
  LV.options.DISPLAY = 0;
  std::cout << std::scientific << std::setprecision(5);
  for( unsigned int k=0; k<10; k++ ){
    I Ipred = mc::Op<I>::mid(Ip[0])+I(-0.5,0.5)*mc::Op<I>::diam(Ip[0])*pow(0.8,k);
    PV PMp0[NP] = { PV( &PMEnv, 0, Ipred ) };
    //LV.hausdorff( Ip, Hxk, NSAMP );
    LV.hausdorff( PMp0, Hxk, Hf, NSAMP );
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


