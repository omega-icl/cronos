const unsigned int NPM   = 3;	// <- Order of polynomial model

#ifdef MC__USE_PROFIL
 #include "mcprofil.hpp"
 typedef INTERVAL I;
#else
 #ifdef MC__USE_FILIB
  #include "mcfilib.hpp"
  typedef filib::interval<double,filib::native_switched,filib::i_mode_extended> I;
 #else
  #ifdef MC__USE_BOOST
   #include "mcboost.hpp"
   typedef boost::numeric::interval_lib::save_state<boost::numeric::interval_lib::rounded_transc_opp<double>> T_boost_round;
   typedef boost::numeric::interval_lib::checking_base<double> T_boost_check;
   typedef boost::numeric::interval_lib::policies<T_boost_round,T_boost_check> T_boost_policy;
   typedef boost::numeric::interval<double,T_boost_policy> I;
  #else
   #include "interval.hpp"
   typedef mc::Interval I;
  #endif
 #endif
#endif

#include "scmodel.hpp"
typedef mc::SCModel<I> PM;
typedef mc::SCVar<I> PV;

#include "aebnd.hpp"
typedef mc::FFGraph<> t_DAG;
typedef mc::AEBND<I,PM,PV> t_AEBND;

int main()
{
 try{
  t_DAG NLE;  // DAG describing the problem

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

  I Ip[NP]  = { I(0.5,1.5) },
    Ix0[NX] = { I(0.01,10.0), I(0.1*PT,0.3*PT) },
    Ix[NX];

  PM PMEnv( NPM );
  PV PMp[NP], PMx[NX], PMx0[NX];
  for( unsigned i=0; i<NP; i++ ) PMp[i].set( &PMEnv, i, Ip[i] );
  for( unsigned i=0; i<NX; i++ ) PMx0[i] = Ix0[i];

  /////////////////////////////////////////////////////////////////////////
  // Bound AE solution set
  t_AEBND BND;

  BND.set_dag( &NLE );
  BND.set_var( NP, P );
  BND.set_dep( NX, X, F );

  BND.options.DISPLEVEL = 1;
  BND.options.MAXIT     = 20;
  BND.options.RTOL      =
  BND.options.ATOL      = 1e-10;
  BND.options.BOUNDER   = t_AEBND::Options::ALGORITHM::AUTO;//GE;//KRAW;//GS;
  BND.options.PRECOND   = t_AEBND::Options::PRECONDITIONING::INVMD;//QRM;//NONE;
  BND.options.BLKDEC    = t_AEBND::Options::DECOMPOSITION::RECUR;//NONE;//DIAG;

  BND.setup();
  std::cout << "\nSuccessful? " << (BND.solve( Ip, Ix, Ix0 )==t_AEBND::NORMAL?"Y\n":"N\n");
  std::cout << "\nSuccessful? " << (BND.solve( PMp, PMx, Ix )==t_AEBND::NORMAL?"Y\n":"N\n");

 }

 catch( t_AEBND::Exceptions &eObj ){
  std::cerr << "Error " << eObj.ierr()
            << eObj.what() << std::endl;
 }

 return 0;
}


