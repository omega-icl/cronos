const unsigned int NPM   = 1;	// <- Order of Taylor model
#define USE_CMODEL		// <- Use Chebyshev models?

#include "aebnd.hpp"

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

int main()
{
  mc::FFGraph NLE;  // DAG describing the problem

  const unsigned NP = 1;  // Parameter dimension [phiB, AB, BB, CB, AT, BT, CT]
  const unsigned NX = 2;  // State dimension [ h | PCO2 ]

  mc::FFVar P[NP];  // Parameters p
  for( unsigned i=0; i<NP; i++ ) P[i].set( &NLE );

  mc::FFVar X[NX];  // Dependents x
  for( unsigned i=0; i<NX; i++ ) X[i].set( &NLE );

  double Cin   = 19.0;   // mmol/L
  double Nin   = 11.0;   // mmol/L
  double Zin   = 17.0;   // mmol/L
  double k11   = 300.0;  // mmol/gCOD
  double mu3m  = 0.14;   // /day
  double KS3   = 0.02;   // gCOD/L
  double KI3   = 16.4;   // gCOD/L
  double KIN   = 1.8;    // mmol/L
  double KC    = 4.9e-1; // µmol/L
  double KN    = 1.1e-3; // µmol/L
  double KVFA  = 1.7e1;  // µmol/L
  double KH2O  = 2.1e-2; // (µmol/L)^2
  double KHCO2 = 2.7e1;  // mmol/Bar
  double MVFA  = 0.064;  // gCOD/mmol
  double kla   = 5.0;    // /day
  double PT    = 1.0;    // Bar

  double S3 = 0.2;    // mmol/L
  double X3 = 1.0;    // gCOD/L
  double N  = Nin;     // mmol/L
  double C  = Cin;     // mmol/L
  double Z  = Zin;     // mmol/L

  mc::FFVar F[NX];  // Equations f(x,p)=0
  mc::FFVar NH3  = KN / ( KN + X[0] ) * N;
  mc::FFVar HCO3 = KC / ( KC + X[0] ) * C;
  mc::FFVar VFA  = KVFA / ( KVFA + X[0] ) * S3 / MVFA;
  mc::FFVar H    = X[0] * 1e3;
  mc::FFVar OH   = KH2O / X[0] * 1e3;
  mc::FFVar mu3  = mu3m * S3 / ( KS3 + S3 + S3*S3/KI3 ) * KIN / ( KIN + NH3 );
  mc::FFVar qCO2 = kla * ( C - HCO3 - KHCO2 * X[1] );
  mc::FFVar qCH4 = k11 * mu3 * X3 * P[0];

  F[0] = OH + HCO3 + VFA - ( N - NH3 ) - Z - H; // Ion balance
  F[1] = PT * qCO2 - X[1] * ( qCO2 + qCH4 );    // Pressure equilibrium

  I Ip[NP]  = { I(0.99,1.01) },
    Ix0[NX] = { I(0.01,10.0), I(0.1*PT,0.3*PT) },
    Ix[NX];

  PM PMEnv( NP, NPM );
  PV PMp[NP], PMx[NX], PMx0[NX];
  for( unsigned i=0; i<NP; i++ ) PMp[i].set( &PMEnv, i, Ip[i] );
  for( unsigned i=0; i<NX; i++ ) PMx0[i] = Ix0[i];

  /////////////////////////////////////////////////////////////////////////
  // Bound AE solution set
  mc::AEBND<I,PM,PV> EX1;

  EX1.set_dag( &NLE );
  EX1.set_var( NP, P );
  EX1.set_dep( NX, X, F );

  EX1.options.DISPLAY = 2;
  EX1.options.BLKDEC  = true;
  EX1.options.MAXIT   = 20;
  EX1.options.RTOL    =
  EX1.options.ATOL    = 0e0;
  EX1.options.BOUNDER = mc::AEBND<I,PM,PV>::Options::AUTO;
  EX1.options.PRECOND = mc::AEBND<I,PM,PV>::Options::INVMD;//QRM;//NONE;

  EX1.setup();
  std::cout << "\nSuccessful? " << (EX1.solve( Ip, Ix, Ix0 )==mc::AEBND<I,PM,PV>::NORMAL?"Y\n":"N\n");
  std::cout << "\nSuccessful? " << (EX1.solve( PMp, PMx, Ix )==mc::AEBND<I,PM,PV>::NORMAL?"Y\n":"N\n");

  return 0;
}


