// ESCAPE Case Study - AD Parameter Estimation

const unsigned int NPM   = 2;	// <- Order of polynomial expansion
const unsigned int NSAMP = 3;	// <- Number of sampling points for inner approx.
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

  const unsigned int NM = 12;
  const double yS1[NM] ={0.74, 0.63, 0.61, 0.58, 0.59, 0.60, 0.61, 0.63, 0.64, 0.66, 0.67, 0.70}; 
  const double yS2[NM] ={2.3, 2.2, 2.2, 2.1, 2.2, 2.3, 2.3, 2.3, 2.4, 2.4, 2.4, 2.5}; 
  const double yC[NM]  ={51.1, 53.4, 54.0, 54.0, 54.0, 54.0, 53.9, 53.9, 53.8, 53.8, 53.8, 53.7};

int main()
{
  mc::FFGraph IVP;  // DAG describing the problem

  double t0 = 0., tf = 3.;  // Time span
  const unsigned NS = 480; //12; // Time stages
  double tk[NS+1]; tk[0] = t0;
  for( unsigned k=0; k<NS; k++ ) tk[k+1] = tk[k] + (tf-t0)/(double)NS;

  const unsigned NP = 3;  // Number of parameters
  const unsigned NX = 6;  // Number of states
  const unsigned NQ = 0;  // Number of state quadratures
  const unsigned NF = 1;  // Number of state functions

  mc::FFVar P[NP];  // Parameters
  for( unsigned int i=0; i<NP; i++ ) P[i].set( &IVP );

  mc::FFVar X[NX];  // States
  for( unsigned int i=0; i<NX; i++ ) X[i].set( &IVP );
  mc::FFVar S1 = X[0];
  mc::FFVar S2 = X[1];
  mc::FFVar C  = X[2];
  mc::FFVar X1 = X[3];
  mc::FFVar X2 = X[4];
  mc::FFVar Z  = X[5];

  //Constant Parameters
  const double mu1 = 1.2;
  const double KS1 = 7.1;
  const double mu2 = 0.74;
  const double KS2 = 9.28;
  const double KI2 = 256.;
  const double kLa = 19.8;;
  const double KH = 16.;
  const double Pt = 1.;
  const double a = 0.5;
  const double D = 0.4;
  const double k1 = 42.14;
  const double k2 = 116.5;
  const double k3 = 268.;
  const double k4 = 50.6;
  const double k5 = 343.6;
  const double k6 = 453.;
  const double S1in = 5.;
  const double S2in = 80.;
  const double Zin = 50.;
  const double Cin = 0.;

  //Pieces
  mc::FFVar mu1S1  = mu1*S1/(S1 + KS1);
  mc::FFVar mu2S2  = mu2*S2/(S2 + KS2 + S2*S2/KI2);
  mc::FFVar phiCO2 = C + S2 - Z + KH*Pt + (k6/kLa)*mu2S2*X2;
  mc::FFVar PCO2   = (phiCO2 - pow(phiCO2*phiCO2 - 4.*KH*Pt*(C + S2 - Z),0.5) )/(2.*KH);
  mc::FFVar qCO2   = kLa*(C + S2 - Z - KH*PCO2);
  
  mc::FFVar RHS[NX];  // Right-hand side function
  RHS[0] = D*(S1in - S1) - k1*mu1S1*X1; //S1
  RHS[1] = D*(S2in - S2) + k2*mu1S1*X1 - k3*mu2S2*X2; //S2
  RHS[2] = D*(Cin - C) - qCO2 + k4*mu1S1*X1 + k5*mu2S2*X2; //C
  RHS[3] = (mu1S1 - a*D)*X1; //X1
  RHS[4] = (mu2S2 - a*D)*X2; //X2
  RHS[5] = D*(Zin - Z); //Z

  mc::FFVar IC[NX];   // Initial value function
  IC[0] = 1.;   //S1
  IC[1] = 5.;   //S2
  IC[2] = P[2]; //C
  IC[3] = P[0]; //X1
  IC[4] = P[1]; //X2
  IC[5] = 50.;  //Z

  if( NS % NM ){std::cout << "MISMATCH, FAILED!" << std::endl; return 1;}
  mc::FFVar FCT[NF*NS];  // State functions
  for( unsigned k=0; k<NF*NS; k++ ) FCT[k] = 0.;
  for( unsigned k=0; k<NM; k++ ){
    FCT[(k+1)*NS/NM-1] = pow((yS1[k] - S1),2) + pow((yS2[k] - S2),2) + pow((yC[k] - C),2);
    //FCT[k*(NS/NM)] = pow((yS1[k] - S1),2) + pow((yS2[k] - S2),2) + pow((yC[k] - C),2);
  }

  I Ip[NP]; Ip[0]= 0.5*I(0.98,1.02); Ip[1] = I(0.98,1.02); Ip[2] = 40.*I(0.98,1.02);
  I *Ixk[NS+1], *Iyk[NS+1];
  for( unsigned k=0; k<=NS; k++ ){
    Ixk[k] = new I[NX];
    Iyk[k] = new I[NF*NX];
  }
  I Iq[NQ], If[NF], Idf[NF*NP];

  PM PMEnv( NP, NPM );
  PV PMp[NP] = { PV( &PMEnv, 0, Ip[0] ), PV( &PMEnv, 1, Ip[1] ), PV( &PMEnv, 2, Ip[2] ) };
  PV *PMxk[NS+1], *PMyk[NS+1];
  for( unsigned k=0; k<=NS; k++ ){
    PMxk[k] = new PV[NX];
    PMyk[k] = new PV[NF*NX];
  }
  PV PMq[NQ], PMf[NF], PMdf[NF*NP];

  /////////////////////////////////////////////////////////////////////////////
  //// SAMPLING
  mc::ODESLVS_GSL<I> LV0;

  LV0.set_dag( &IVP );
  LV0.set_state( NX, X );
  LV0.set_parameter( NP, P );
  LV0.set_differential( NX, RHS );
  LV0.set_initial( NX, IC );
  LV0.set_function( NF, NS, FCT );

  LV0.options.DISPLAY   = 1;
#if defined( SAVE_RESULTS )
  LV0.options.RESRECORD = true;
#endif
  LV0.options.ATOL      = LV0.options.RTOL = 1e-8;
  LV0.options.HMAX      = 1e-3;
  LV0.options.INTMETH   = mc::ODESLV_GSL<I>::Options::MSBDF;

  // Approximate adjoint bounds
  std::cout << "\nNON_VALIDATED INTEGRATION - APPROXIMATE ENCLOSURE OF REACHABLE SET:\n\n";
  //LV0.bounds( NS, tk, Ip, Ixk, Iq, If, NSAMP );
  LV0.bounds_ASA( NS, tk, Ip, Ixk, 0, If, Iyk, Idf, NSAMP );
#if defined( SAVE_RESULTS )
  std::ofstream apprecSTA("test6_APPROX_STA.dat", std::ios_base::out );
  std::ofstream apprecADJ("test6_APPROX_ADJ.dat", std::ios_base::out );
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
  LV.options.WRAPMIT   = mc::ODEBND_GSL<I,PM,PV>::Options::DINEQ;//ELLIPS;//
  LV.options.HMAX      = 1e-2;
  LV.options.H0      = 1e-6;

#else // SUNDIALS integrator
  mc::ODEBNDS_SUNDIALS<I,PM,PV> LV;

  LV.options.DISPLAY   = 1;
#if defined( SAVE_RESULTS )
  LV.options.RESRECORD = true;
#endif
  LV.options.ATOL      = LV.options.RTOL  = 1e-8;
  LV.options.ATOLB     = LV.options.RTOLB = 1e-8;
  LV.options.INTMETH   = mc::ODEBNDS_SUNDIALS<I,PM,PV>::Options::MSADAMS;
  LV.options.JACAPPROX = mc::ODEBNDS_SUNDIALS<I,PM,PV>::Options::CV_DIAG;//CV_DENSE;
  LV.options.ORDMIT    = 1; //PMp->nord();
  LV.options.WRAPMIT   = mc::ODEBNDS_SUNDIALS<I,PM,PV>::Options::DINEQ;//ELLIPS;//NONE;
#endif

  LV.set_dag( &IVP );
  LV.set_state( NX, X );
  LV.set_parameter( NP, P );
  LV.set_differential( NX, RHS );
  LV.set_initial( NX, IC );
  LV.set_function( NF, NS, FCT );

  std::cout << "\nCONTINUOUS SET-VALUED INTEGRATION - INTERVAL ENCLOSURE OF REACHABLE SET:\n\n";
  //LV.bounds( NS, tk, Ip, Ixk, Iq, If );
  LV.bounds_ASA( NS, tk, Ip, Ixk, Iq, If, Iyk, Idf );
#if defined( SAVE_RESULTS )
  std::ofstream direcISTA( "test6_DINEQI_STA.dat", std::ios_base::out );
  std::ofstream direcIADJ( "test6_DINEQI_ADJ.dat", std::ios_base::out );
  LV.record( direcISTA, direcIADJ );
#endif

  std::cout << "\nCONTINUOUS SET-VALUED INTEGRATION - POLYNOMIAL MODEL ENCLOSURE OF REACHABLE SET:\n\n";
  //LV.bounds( NS, tk, PMp, PMxk, PMq, PMf );
  LV.bounds_ASA( NS, tk, PMp, PMxk, PMq, PMf, PMyk, PMdf );
#if defined( SAVE_RESULTS )
  std::ofstream direcPMSTA( "test6_DINEQPM_STA.dat", std::ios_base::out );
  std::ofstream direcPMADJ( "test6_DINEQPM_ADJ.dat", std::ios_base::out );
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
