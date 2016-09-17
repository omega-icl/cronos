/// ANAEROBIC DIGESTER ///
const unsigned int NPM   = 3;	// <- Order of polynomial expansion
const unsigned int NSAMP = 3;	// <- Number of sampling points for inner approx.
#define SAVE_RESULTS    		// <- Whether to save bounds to file
#define USE_CMODEL		        // <- whether to use Chebyshev models or Taylor models

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

const unsigned int NM = 12;
const double yS1[NM] ={0.74, 0.63, 0.61, 0.58, 0.59, 0.60, 0.61, 0.63, 0.64, 0.66, 0.67, 0.70}; 
const double yS2[NM] ={2.3, 2.2, 2.2, 2.1, 2.2, 2.3, 2.3, 2.3, 2.4, 2.4, 2.4, 2.5}; 
const double yC[NM]  ={51.1, 53.4, 54.0, 54.0, 54.0, 54.0, 53.9, 53.9, 53.8, 53.8, 53.8, 53.7};

int main()
{
  mc::FFGraph IVP;  // DAG describing the problem

  double t0 = 0., tf = 3.;  // Time span
  const unsigned NS = 1*NM; //12; // Time stages
  assert( !(NS % NM) );
  double tk[NS+1]; tk[0] = t0;
  for( unsigned k=0; k<NS; k++ ) tk[k+1] = tk[k] + (tf-t0)/(double)NS;

  const unsigned NP = 3;  // Number of parameters
  const unsigned NX = 6;  // Number of states
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

  mc::FFVar FCT[NF*NS];  // State functions
  for( unsigned k=0; k<NF*NS; k++ ) FCT[k] = 0.;
  for( unsigned k=0; k<NM; k++ ){
    FCT[(k+1)*NS/NM-1] = pow((yS1[k] - S1),2) + pow((yS2[k] - S2),2) + pow((yC[k] - C),2);
    //FCT[k*(NS/NM)] = pow((yS1[k] - S1),2) + pow((yS2[k] - S2),2) + pow((yC[k] - C),2);
  }

  I Ip[NP]; Ip[0]= 0.5*I(0.98,1.02); Ip[1] = I(0.98,1.02); Ip[2] = 40.*I(0.98,1.02);
  I *Ixk[NS+1], *Iyk[NS+1], *Ixpk[NS+1];
  for( unsigned k=0; k<=NS; k++ ){
    Ixk[k] = new I[NX];
    Iyk[k] = new I[NF*NX];
    Ixpk[k] = new I[NP*NX];
  }
  I If[NF], Ifp[NF*NP];

  PM PMEnv( NP, NPM );
  PV PMp[NP] = { PV( &PMEnv, 0, Ip[0] ), PV( &PMEnv, 1, Ip[1] ), PV( &PMEnv, 2, Ip[2] ) };
  PV *PMxk[NS+1], *PMyk[NS+1], *PMxpk[NS+1];
  for( unsigned k=0; k<=NS; k++ ){
    PMxk[k] = new PV[NX];
    PMyk[k] = new PV[NF*NX];
    PMxpk[k] = new PV[NP*NX];
  }
  PV PMf[NF], PMfp[NF*NP];

  /////////////////////////////////////////////////////////////////////////////
  //// DIFFERENTIAL INEQUALITIES
  mc::ODEBNDS_SUNDIALS<I,PM,PV> AD;

  AD.options.DISPLAY   = 1;
#if defined( SAVE_RESULTS )
  AD.options.RESRECORD = true;
#endif
  AD.options.ATOL      = AD.options.ATOLB      = AD.options.ATOLS  = 1e-10;
  AD.options.RTOL      = AD.options.RTOLB      = AD.options.RTOLS  = 1e-8;
  AD.options.ETOL      = AD.options.ETOLB      = AD.options.ETOLS  = 1e-20;
  AD.options.NMAX      = AD.options.ASACHKPT   = 25000;
  AD.options.MAXFAIL   = 20;
  AD.options.FSAERR    = true;
  AD.options.QERRS     = true;
  //AD.options.AUTOTOLS  = true;
  AD.options.INTMETH   = mc::ODEBNDS_SUNDIALS<I,PM,PV>::Options::MSBDF;//MSSADAMS;//MSBDF;
  AD.options.JACAPPROX = mc::ODEBNDS_SUNDIALS<I,PM,PV>::Options::CV_DIAG;//CV_DENSE;//CV_DIAG;
  AD.options.ORDMIT    = -2; //PMp->nord();
  AD.options.WRAPMIT   = mc::ODEBNDS_SUNDIALS<I,PM,PV>::Options::ELLIPS;//NONE;//DINEQ;
  AD.options.QERRB     = true;
  AD.options.QSCALE    = 1e-5;
  AD.options.HMIN      = 1e-14;

  AD.options.ODESLV.NMAX = AD.options.ODESLV.ASACHKPT = 25000;
  AD.options.ODESLV.INTMETH   = mc::ODESLVS_SUNDIALS::Options::MSBDF;//MSADAMS;//MSBDF;
  AD.options.ODESLV.JACAPPROX = mc::ODESLVS_SUNDIALS::Options::CV_DIAG;//CV_DENSE;//CV_DIAG;

  AD.set_dag( &IVP );
  AD.set_state( NX, X );
  AD.set_parameter( NP, P );
  AD.set_differential( NX, RHS );
  AD.set_initial( NX, IC );
  AD.set_function( NS, NF, FCT );

  std::ofstream ofSTA, ofSEN;

  std::cout << "\nNON_VALIDATED INTEGRATION - INNER-APPROXIMATION OF REACHABLE SET AND ADJOINT SENSITIVITY:\n\n";
  AD.bounds_ASA( NS, tk, Ip, Ixk, If, Iyk, Ifp, NSAMP );
#if defined( SAVE_RESULTS )
  ofSTA.open( "test4_APPROX_STA.dat", std::ios_base::out );
  ofSEN.open( "test4_APPROX_ASA.dat", std::ios_base::out );
  AD.record( ofSTA, ofSEN );
  ofSTA.close();
  ofSEN.close();
#endif
  { int dum; std::cout << "PAUSED--"; std::cin >> dum; }
/*
  std::cout << "\nCONTINUOUS SET-VALUED INTEGRATION - INTERVAL ENCLOSURE OF REACHABLE SET AND ADJOINT SENSITIVITY:\n\n";
  AD.bounds_ASA( NS, tk, Ip, Ixk, If, Iyk, Ifp );
#if defined( SAVE_RESULTS )
  ofSTA.open( "test4_DINEQI_STA.dat", std::ios_base::out );
  ofSEN.open( "test4_DINEQI_ASA.dat", std::ios_base::out );
  AD.record( ofSTA, ofSEN );
  ofSTA.close();
  ofSEN.close();
#endif
*/
  std::cout << "\nCONTINUOUS SET-VALUED INTEGRATION - POLYNOMIAL MODEL ENCLOSURE OF REACHABLE SET AND ADJOINT SENSITIVITY:\n\n";
  //AD.bounds( NS, tk, PMp, PMxk, PMf );
  AD.bounds_ASA( NS, tk, PMp, PMxk, PMf, PMyk, PMfp );
#if defined( SAVE_RESULTS )
  ofSTA.open( "test4_DINEQPM_STA.dat", std::ios_base::out );
  ofSEN.open( "test4_DINEQPM_ASA.dat", std::ios_base::out );
  AD.record( ofSTA, ofSEN );
  ofSTA.close();
  ofSEN.close();
#endif
  { int dum; std::cout << "--PAUSED "; std::cin >> dum; }

  std::cout << "\nNON_VALIDATED INTEGRATION - INNER-APPROXIMATION OF REACHABLE SET AND FORWARD SENSITIVITY:\n\n";
  AD.bounds_FSA( NS, tk, Ip, Ixk, If, Ixpk, Ifp, NSAMP );
#if defined( SAVE_RESULTS )
  ofSTA.open( "test4_APPROX_STA.dat", std::ios_base::out );
  ofSEN.open( "test4_APPROX_FSA.dat", std::ios_base::out );
  AD.record( ofSTA, ofSEN );
  ofSTA.close();
  ofSEN.close();
#endif
  { int dum; std::cout << "PAUSED--"; std::cin >> dum; }
/*
  std::cout << "\nCONTINUOUS SET-VALUED INTEGRATION - INTERVAL ENCLOSURE OF REACHABLE SET AND FORWARD SENSITIVITY:\n\n";
  //AD.bounds( NS, tk, Ip, Ixk, If );
  AD.bounds_FSA( NS, tk, Ip, Ixk, If, Ixpk, Ifp );
#if defined( SAVE_RESULTS )
  ofSTA.open( "test4_DINEQI_STA.dat", std::ios_base::out );
  ofSEN.open( "test4_DINEQI_FSA.dat", std::ios_base::out );
  AD.record( ofSTA, ofSEN );
  ofSTA.close();
  ofSEN.close();
#endif
*/
  std::cout << "\nCONTINUOUS SET-VALUED INTEGRATION - POLYNOMIAL MODEL ENCLOSURE OF REACHABLE SET AND FORWARD SENSITIVITY:\n\n";
  AD.bounds_FSA( NS, tk, PMp, PMxk, PMf, PMxpk, PMfp );
#if defined( SAVE_RESULTS )
  ofSTA.open( "test4_DINEQPM_STA.dat", std::ios_base::out );
  ofSEN.open( "test4_DINEQPM_FSA.dat", std::ios_base::out );
  AD.record( ofSTA, ofSEN );
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

