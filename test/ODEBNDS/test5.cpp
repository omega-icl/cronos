// ESCAPE26 Case Study 1 - Exothermic Batch Reactor

const unsigned int NPM   = 3;	// <- Order of poynomial expansion
const unsigned int NSAMP = 20;	// <- Number of sampling points for inner approx.
#define SAVE_RESULTS		// <- Whether to save bounds to file
#define USE_CMODEL		// <- whether to use Chebyshev models or Taylor models
#define USE_SUNDIALS		// <- whether to use SUNDIALS or GSL integrator

#include "odeslvs_gsl.hpp"
#ifndef USE_SUNDIALS
  #include "odebnds_gsl.hpp"
#else
  #include "odebnds_sundials.hpp"
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

int main()
{
  mc::FFGraph IVP;  // DAG describing the problem

  double t0 = 0., tf = 600.;  // Time span
  const unsigned int NS = 100;  // Time stages
  double tk[NS+1]; tk[0] = t0;
  for( unsigned k=0; k<NS; k++ ) tk[k+1] = tk[k] + (tf-t0)/(double)NS;

  const unsigned NP = 2;  // Number of parameters
  const unsigned NX = 2;  // Number of states
  const unsigned NF = 1;  // Number of state functions

  mc::FFVar P[NP];  // Parameters
  for( unsigned int i=0; i<NP; i++ ) P[i].set( &IVP );

  mc::FFVar X[NX];  // States
  for( unsigned int i=0; i<NX; i++ ) X[i].set( &IVP );
  mc::FFVar X1 = X[0];
  mc::FFVar T =  X[1];

  //Constant Parameters
  const double k0 = 0.022;
  const double CA0 = 10.;
  const double V = 0.1;
  const double Cp = 60.;
  const double Ea = 6000.;
  const double R = 8.314;
  const double DHR = -140000.;
  const double UA = 3.;

  mc::FFVar RHS[NX];  // Right-hand side function
  RHS[0] = k0*exp(-Ea/(R*T))*(1. - X1);
  RHS[1] = (P[1] - T)*UA/(CA0*V*Cp) - (1 - X1)*exp(-Ea/(R*T))*(DHR*k0)/Cp;

  mc::FFVar IC[NX];   // Initial value function
  IC[0] = 0.;
  IC[1] = P[0];

  mc::FFVar FCT[NF];  // State functions
  FCT[0] = X1;

  I Ip[NP] = { I(350.,370.), I(290.,310.) };
  I *Ixk[NS+1], *Iyk[NS+1];
  for( unsigned k=0; k<=NS; k++ ){
    Ixk[k] = new I[NX];
    Iyk[k] = new I[NF*NX];
  }
  I If[NF], Idf[NF*NP];

  PM PMEnv( NP, NPM );
  PV PMp[NP] = { PV( &PMEnv, 0, Ip[0] ), PV( &PMEnv, 1, Ip[1] ) };
  PV *PMxk[NS+1], *PMyk[NS+1];
  for( unsigned k=0; k<=NS; k++ ){
    PMxk[k] = new PV[NX];
    PMyk[k] = new PV[NF*NX];
  }
  PV PMf[NF], PMdf[NF*NP];

  /////////////////////////////////////////////////////////////////////////////
  //// SAMPLING
  mc::ODESLVS_GSL<I> LV0;

  LV0.set_dag( &IVP );
  LV0.set_state( NX, X );
  LV0.set_parameter( NP, P );
  LV0.set_differential( NX, RHS );
  LV0.set_initial( NX, IC );
  LV0.set_function( NF, FCT );

  LV0.options.DISPLAY   = 1;
#if defined( SAVE_RESULTS )
  LV0.options.RESRECORD = true;
#endif
  LV0.options.ATOL      = LV0.options.RTOL = 1e-7;
  LV0.options.INTMETH   = mc::ODESLV_GSL<I>::Options::MSBDF;

  // Approximate adjoint bounds
  std::cout << "\nNON_VALIDATED INTEGRATION - APPROXIMATE ENCLOSURE OF REACHABLE SET:\n\n";
  //LV0.bounds( NS, tk, Ip, Ixk, 0, If, NSAMP );
  LV0.bounds_ASA( NS, tk, Ip, Ixk, 0, If, Iyk, Idf, NSAMP );
#if defined( SAVE_RESULTS )
  std::ofstream apprecSTA("test5_APPROX_STA.dat", std::ios_base::out );
  std::ofstream apprecADJ("test5_APPROX_ADJ.dat", std::ios_base::out );
  LV0.record( apprecSTA, apprecADJ ); 
#endif

  /////////////////////////////////////////////////////////////////////////////
  //// DIFFERENTIAL INEQUALITIES
#ifndef USE_SUNDIALS // GSL integrator
  mc::ODEBNDS_GSL<I,PM,PV> LV;

  LV.options.DISPLAY = 1;
#if defined( SAVE_RESULTS )
  LV.options.RESRECORD = true;
#endif
  LV.options.ATOL      = LV.options.RTOL = 1e-8;
  LV.options.ORDMIT    = 1; //PMp->nord();
  LV.options.WRAPMIT   = mc::ODEBND_GSL<I,PM,PV>::Options::ELLIPS;//DINEQ; //NONE;
  //LV.options.HMAX      = 1e-2;
  //LV.options.H0        = 1e-6;

#else // SUNDIALS integrator
  mc::ODEBNDS_SUNDIALS<I,PM,PV> LV;

  LV.options.DISPLAY   = 1;
#if defined( SAVE_RESULTS )
  LV.options.RESRECORD = true;
#endif
  LV.options.ATOL      = LV.options.ATOLB = 1e-10;
  LV.options.RTOL      = LV.options.RTOLB = 1e-8;
  LV.options.INTMETH   = mc::ODEBNDS_SUNDIALS<I,PM,PV>::Options::MSADAMS;
  LV.options.JACAPPROX = mc::ODEBNDS_SUNDIALS<I,PM,PV>::Options::CV_DIAG;//CV_DENSE;
  LV.options.ORDMIT    = 1;//PMp->nord();
  LV.options.WRAPMIT   = mc::ODEBNDS_SUNDIALS<I,PM,PV>::Options::ELLIPS;//DINEQ;//NONE;
  //LV.options.QSCALE    = 1e0;
  LV.options.HMIN      = 1e-12;
#endif

  LV.set_dag( &IVP );
  LV.set_state( NX, X );
  LV.set_parameter( NP, P );
  LV.set_differential( NX, RHS );
  LV.set_initial( NX, IC );
  LV.set_function( NF, FCT );

  std::cout << "\nCONTINUOUS SET-VALUED INTEGRATION - INTERVAL ENCLOSURE OF REACHABLE SET:\n\n";
  LV.bounds_ASA( NS, tk, Ip, Ixk, 0, If, Iyk, Idf );
#if defined( SAVE_RESULTS )
  std::ofstream direcISTA( "test5_DINEQI_STA.dat", std::ios_base::out );
  std::ofstream direcIADJ( "test5_DINEQI_ADJ.dat", std::ios_base::out );
  LV.record( direcISTA, direcIADJ );
#endif

  std::cout << "\nCONTINUOUS SET-VALUED INTEGRATION - POLYNOMIAL MODEL ENCLOSURE OF REACHABLE SET:\n\n";
  //LV.bounds( NS, tk, PMp, PMxk, 0, PMf );
  LV.bounds_ASA( NS, tk, PMp, PMxk, 0, PMf, PMyk, PMdf );
#if defined( SAVE_RESULTS )
  std::ofstream direcPMSTA( "test5_DINEQPM_STA.dat", std::ios_base::out );
  std::ofstream direcPMADJ( "test5_DINEQPM_ADJ.dat", std::ios_base::out );
  LV.record( direcPMSTA, direcPMADJ );
#endif

  // Clean up
  for( unsigned k=0; k<=NS; k++ ){
    delete[] Ixk[k];
    delete[] Iyk[k];
    delete[] PMxk[k];
    delete[] PMyk[k];
  }

  return 0;
}
