///////////////////////////////////////////////////////////////////////////////
// SECOND EXAMPLE IN:
//   Jai Rayjaguru, Mario E Villanueva, Boris Houska, Benoit Chachuat
//   Continuous-time enclosures for uncertain implicit differential equations
//   IFAC-PapersOnLine (2015) 48(8):94-99, DOI 10.1016/j.ifacol.2015.08.163
///////////////////////////////////////////////////////////////////////////////
const unsigned int NPM   = 5;	// <- Order of polynomial model
const unsigned int NSAMP = 10;	// <- Number of sampling points
#define SAVE_RESULTS		// <- Whether to save bounds to file
#undef  TEST_CONVERGENCE	// <- Whether to test Hausdorff convergence of bounds
#define USE_CMODEL		// <- whether to use Chebyshev models or Taylor models
#define USE_SUNDIALS		// <- whether to use SUNDIALS or GSL integrator
///////////////////////////////////////////////////////////////////////////////

#include "odeslv_gsl.hpp"
#ifndef USE_SUNDIALS
  #include "iodebnd_gsl.hpp"
#else
  #include "iodebnd_sundials.hpp"
#endif

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

  // Model parameters and initial conditions
  double D = 0.1;       // /day
  double f1 =  0.3;     // gCOD/gCOD
  double f2 =  0.4;     // gCOD/gCOD
  double Sin = 10.0;    // gCOD/L
  double Cin = 19.0;    // mmol/L
  double Nin = 11.0;    // mmol/L
  double Zin = 17.0;    // mmol/L

  double k1 = 12.5;     // gCOD/gCOD
  double k2 =  6.2;     // mmol/gCOD
  double k3 = 11.5;     // gCOD/gCOD
  double k4 = 30.0;     // mmol/gCOD
  double k5 =  9.1;     // gCOD/gCOD
  double k6 =  8.1;     // gCOD/gCOD
  double k7 = 54.0;     // mmol/gCOD
  double k8 = 30.0;     // mmol/gCOD
  double k9 = 20.0;     // gCOD/gCOD
  double k10 = 6.2;     // mmol/gCOD
  double k11 = 300.0;   // mmol/gCOD
  double k12 = 200.0;   // mmol/gCOD
  double mu1m = 0.30;   // /day
  double KS1 = 2.11;    // gCOD/L
  double mu2m = 0.053;  // /day
  double KS2 = 0.056;   // gCOD/L
  double mu3m = 0.14;   // /day
  double KS3 = 0.02;    // gCOD/L
  double KI3 = 16.4;    // gCOD/L
  double KIN = 1.8;     // mmol/L
  double KC = 4.9e-1;   // µmol/L
  double KN = 1.1e-3;   // µmol/L
  double KVFA = 1.7e1;  // µmol/L
  double KH2O = 2.1e-2; // (µmol/L)^2
  double KHCO2 = 2.7e1; // mmol/Bar
  double MVFA = 0.064;  // gCOD/mmol
  double kla = 5.0;     // /day
  double PT = 1.0;      // Bar

  double S1_0 = Sin/2.; // gCOD/L
  double S2_0 = Sin/2.; // gCOD/L
  double S3_0 = 0.2;    // mmol/L
  double X1_0 = 1.0;    // gCOD/L
  double X2_0 = 1.0;    // gCOD/L
  double X3_0 = 1.0;    // gCOD/L
  double N_0 = Nin;     // mmol/L
  double C_0 = Cin;     // mmol/L
  double Z_0 = Zin;     // mmol/L

///////////////////////////////////////////////////////////////////////////////
int main()
{
  mc::FFGraph IVP;  // DAG describing the problem

  double t0 = 0., tf = 10.;   // Time span
  const unsigned int NS = 200;  // Time stages
  double tk[NS+1]; tk[0] = t0;
  for( unsigned k=0; k<NS; k++ ) tk[k+1] = tk[k] + (tf-t0)/(double)NS;

  const unsigned NP = 3;  // Parameter dimension
  const unsigned NX = 11; // State dimension [ S1 | S2 | S3 | X1 | X2 | X3 | N | C | Z | h | PCO2 ]
  const unsigned ND = 9;  // Differential state dimension [ S1 | S2 | S3 | X1 | X2 | X3 | N | C | Z ]

  mc::FFVar P[NP];  // Parameter array
  for( unsigned i=0; i<NP; i++ ) P[i].set( &IVP );

  mc::FFVar X[NX], DX[NX];  // State and state derivative arrays
  for( unsigned i=0; i<NX; i++ ) X[i].set( &IVP );
  for( unsigned i=0; i<NX; i++ ) DX[i].set( &IVP );

  mc::FFVar IC[NX];  // Initial value function
  IC[0] = X[0] - S1_0;
  IC[1] = X[1] - S2_0;
  IC[2] = X[2] - S3_0;
  IC[3] = X[3] - X1_0 * P[0];
  IC[4] = X[4] - X2_0 * P[1];
  IC[5] = X[5] - X3_0 * P[2];
  IC[6] = X[6] - N_0;
  IC[7] = X[7] - C_0;
  IC[8] = X[8] - Z_0;

  mc::FFVar RHS[NX], DE[NX], AE[NX-ND];  // Differential and algebraic equations
  mc::FFVar NH3 = KN / ( KN + X[9] ) * X[6];
  mc::FFVar HCO3 = KC / ( KC + X[9] ) * X[7];
  mc::FFVar VFA = KVFA / ( KVFA + X[9] ) * X[2] / MVFA;
  mc::FFVar H = X[9] * 1e3;
  mc::FFVar OH = KH2O / X[9] * 1e3;
  mc::FFVar mu1 = mu1m * X[0] / ( KS1*X[3] + X[0] );
  mc::FFVar mu2 = mu2m * X[1] / ( KS2*X[4] + X[1] );
  mc::FFVar mu3 = mu3m * X[2] / ( KS3 + X[2] + X[2]*X[2]/KI3 ) * KIN / ( KIN + NH3 );
  mc::FFVar qCO2 = kla * ( X[7] - HCO3 - KHCO2 * X[10] );
  mc::FFVar qCH4 = k11 * mu3 * X[5];
  RHS[0] = D * ( f1*Sin - X[0] ) - k1*mu1*X[3];
  RHS[1] = D * ( f2*Sin - X[1] ) - k5*mu2*X[4];
  RHS[2] = - D * X[2] + k3*mu1*X[3] + k6*mu2*X[4] - k9*mu3*X[5];
  RHS[3] = ( mu1 - D ) * X[3];
  RHS[4] = ( mu2 - D ) * X[4];
  RHS[5] = ( mu3 - D ) * X[5];
  RHS[6] = D * ( Nin - X[6] ) - k2*mu1*X[3] + k7*mu2*X[4] - k10*mu3*X[5];
  RHS[7] = D * ( Cin - X[7] ) + k4*mu1*X[3] + k8*mu2*X[4] + k12*mu3*X[5] - qCO2;
  RHS[8] = D * ( Zin - X[8] );
  AE[0]  = OH + HCO3 + VFA - ( X[6] - NH3 ) - X[8] - H; // Ion balance
  AE[1]  = PT * qCO2 - X[10] * ( qCO2 + qCH4 );         // Pressure equilibrium
  for( unsigned int id=0; id<ND; id++) DE[id] = DX[id] - RHS[id];

  I Ip[NP]  = { I(0.99,1.01), I(0.99,1.01), I(0.99,1.01) }, 
    Ix0[NX] = { I(S1_0), I(S2_0), I(S3_0), X1_0*Ip[0], X2_0*Ip[1], X3_0*Ip[2],
                I(N_0), I(C_0), I(Z_0), I(0.01,10.0), I(0.1*PT,0.3*PT) };
  I *Ixk[NS+1];
  for( unsigned k=0; k<=NS; k++ )
    Ixk[k] = new I[NX];

  PM PMEnv( NP, NPM );
  PV PMp[NP];
  for( unsigned i=0; i<NP; i++ ) PMp[i].set( &PMEnv, i, Ip[i] );
  PV *PMxk[NS+1];
  double *Hxk[NS+1];
  for( unsigned k=0; k<=NS; k++ ){
    PMxk[k] = new PV[NX];
    Hxk[k] = new double[NX];
  }

#if defined( GET_SOLUTION_SET )
  /////////////////////////////////////////////////////////////////////////
  // Bound ODE trajectories - sampling
  const mc::FFVar* AEJAC = IVP.FAD( NX-ND, AE, NX, X );//, DX );
  mc::FFVar UODE[NX-ND];
  for( unsigned ia=0; ia<NX-ND; ia++ ){
    UODE[ia] = 0.;
    for( unsigned jx=0; jx<ND; jx++ )
      UODE[ia] += AEJAC[ia*NX+jx] * RHS[jx];
    for( unsigned jx=0; jx<NX-ND; jx++ )
      UODE[ia] += AEJAC[ia*NX+ND+jx] * DX[ND+jx];
  }
  delete[] AEJAC;
  //std::ofstream o_sol1( "test4_UODE.dot", std::ios_base::out );
  //IVP.dot_script( (NX-ND)*NX, UODE, o_sol1 );
  //o_sol1.close();
  //std::list<const mc::FFOp*> op_UODE = IVP.subgraph( NX-ND, UODE );
  //IVP.output( op_UODE );
  //throw( 1 );
  mc::AEBND<I> DALG;
  DALG.set_dag( &IVP );
  mc::FFVar PX[NP+NX];
  for( unsigned i=0; i<NP; i++ ) PX[i] = P[i];
  for( unsigned i=0; i<NX; i++ ) PX[NP+i] = X[i];
  DALG.set_par( NP+NX, PX );
  DALG.set_dep( NX-ND, DX+ND, UODE );
  DALG.options.DISPLAY = 2;
  DALG.options.BLKDEC  = false;
  DALG.options.BOUNDER = mc::AEBND<I>::Options::GE;
  DALG.setup();
  std::cout << "\nSuccessful? " << (DALG.solve( RHS+ND )==mc::AEBND<I>::NORMAL?"Y\n":"N\n");
  //std::ofstream o_sol2( "test4_ODE.dot", std::ios_base::out );
  //IVP.dot_script( NX-ND, DE+ND, o_sol2 );
  //IVP.dot_script( NX, DE, o_sol2 );
  //o_sol2.close();
  //std::list<const mc::FFOp*> op_DE = IVP.subgraph( NX-ND, DE+ND );
  //IVP.output( op_DE );
  //throw( 1. );

  IC[9]  = 1.39982e-01;
  IC[10] = 1.19966e-01 - 2.70665e-02*(P[2]-1.);

  mc::ODESLV_GSL<I> LV0;
  LV0.set_dag( &IVP );
  LV0.set_state( NX, X );
  LV0.set_parameter( NP, P );
  LV0.set_differential( NX, RHS );
  LV0.set_initial( NX, IC );

  LV0.options.DISPLAY = 1;
  LV0.options.ATOL    = LV0.options.RTOL = 1e-8;
  LV0.options.H0      = 1e-6;
  LV0.options.INTMETH = mc::ODESLV_GSL<I>::Options::MSBDF;
#if defined( SAVE_RESULTS )
  LV0.options.RESRECORD = true;
#endif

  std::cout << "\nNON_VALIDATED INTEGRATION - APPROXIMATE ENCLOSURE OF REACHABLE SET:\n\n";
  LV0.bounds( NS, tk, Ip, Ixk, 0, NSAMP );
#if defined( SAVE_RESULTS )
  std::ofstream apprec( "approximate.out", std::ios_base::out );
  LV0.record( apprec );
#endif
  throw(1);
#endif

  /////////////////////////////////////////////////////////////////////////
  // Bound ODE trajectories - differential inequalities
#ifndef USE_SUNDIALS // GSL integrator
  mc::IODEBND_GSL<I,PM,PV> LV2;
#else // SUNDIALS integrator
  mc::IODEBND_SUNDIALS<I,PM,PV> LV2;
#endif
  LV2.set_dag( &IVP );
  LV2.set_state( NX, X, DX );
  LV2.set_parameter( NP, P );
  LV2.set_differential( ND, DE );
  LV2.set_algebraic( NX-ND, AE );
  LV2.set_initial( ND, IC );
  LV2.set_apriori( Ix0 );

#ifndef USE_SUNDIALS
  LV2.options.WRAPMIT   = mc::IODEBND_GSL<I,PM,PV>::Options::ELLIPS;//DINEQ;//NONE;
#else
  LV2.options.INTMETH   = mc::IODEBND_SUNDIALS<I,PM,PV>::Options::MSADAMS;
  LV2.options.JACAPPROX = mc::IODEBND_SUNDIALS<I,PM,PV>::Options::CV_DIAG;//CV_DENSE;//
  LV2.options.WRAPMIT   = mc::IODEBND_SUNDIALS<I,PM,PV>::Options::ELLIPS;//DINEQ;//NONE;
#endif
  LV2.options.DISPLAY   = 1;
#if defined( SAVE_RESULTS )
  LV2.options.RESRECORD = true;
#endif
  LV2.options.NMAX      = 10000;
  LV2.options.HMIN      = 1e-6;
  LV2.options.ORDMIT    = 4;//NPM; //1;
  LV2.options.ATOL      = 1e-8;
  LV2.options.RTOL      = 1e-7;
  LV2.options.RHSBNDOPT.DISPLAY = 0;
  LV2.options.RHSBNDOPT.BLKDEC  = false;
  LV2.options.RHSBNDOPT.MAXIT   = 20;
  LV2.options.RHSBNDOPT.RTOL    = 1e-10;
  LV2.options.RHSBNDOPT.BOUNDER = mc::AEBND<I,PM,PV>::Options::GE;
  LV2.options.RHSBNDOPT.PRECOND = mc::AEBND<I,PM,PV>::Options::INVMD;//QRMID;//NONE;
  LV2.options.ICBNDOPT.DISPLAY = 0;
  LV2.options.ICBNDOPT.BLKDEC  = false;
  LV2.options.ICBNDOPT.MAXIT   = 20;
  LV2.options.ICBNDOPT.RTOL    = 1e-10;
  LV2.options.ICBNDOPT.BOUNDER = mc::AEBND<I,PM,PV>::Options::AUTO;
  LV2.options.ICBNDOPT.PRECOND = mc::AEBND<I,PM,PV>::Options::INVMD;//QRMID;//NONE;
/*
  std::cout << "\nCONTINUOUS SET-VALUED INTEGRATION - INTERVAL ENCLOSURE OF REACHABLE SET:\n\n";
  LV2.bounds( NS, tk, Ip, Ixk );
#if defined( SAVE_RESULTS )
  std::ofstream bnd2recI( "test4_DINEQI_STA.dat", std::ios_base::out );
  LV2.record( bnd2recI );
#endif
  //LV2.hausdorff( NS, tk, Ip, Hxk, LV0, NSAMP );
*/
  std::cout << "\nCONTINUOUS SET-VALUED INTEGRATION - POLYNOMIAL MODEL ENCLOSURE OF REACHABLE SET:\n\n";
  LV2.bounds( NS, tk, PMp, PMxk );
#if defined( SAVE_RESULTS )
  std::ofstream bnd2recPM( "test4_DINEQPM_STA.dat", std::ios_base::out );
  LV2.record( bnd2recPM );
#endif
  //LV2.hausdorff( NS, tk, PMp, Hxk, LV0, NSAMP );

#if defined( TEST_CONVERGENCE )
  LV2.options.DISPLAY = 0;
  std::cout << std::scientific << std::setprecision(5);
  for( unsigned int k=0; k<20; k++ ){
    TV TMp0[NP] = { TV( &TMEnv, 0, Ip[0]/pow(1.5,k) ) };
    //LV2.hausdorff( NS, tk, Ip, Hxk, NSAMP );
    LV2.hausdorff( NS, tk, TMp0, Hxk, NSAMP );
    std::cout << mc::Op<I>::diam( TMp0[0].B() ) << "  " << Hxk[NS][0] << std::endl;
  }
#endif

  for( unsigned k=0; k<=NS; k++ ){
    delete[] Ixk[k];
    delete[] PMxk[k];
    delete[] Hxk[k];
  }

  return 0;
}


