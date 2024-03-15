const unsigned int NPM   = 5;	// <- Order of Taylor/Chebyshev model
#define USE_CMODEL		// <- Use Chebyshev models?
#undef  MC__AEBND_SHOW_PRECONDITIONING

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

#ifdef USE_CMODEL
  #include "cmodel.hpp"
  typedef mc::CModel<I> PM;
  typedef mc::CVar<I> PV;
#else
  #include "tmodel.hpp"
  typedef mc::TModel<I> PM;
  typedef mc::TVar<I> PV;
#endif

#include "aebnd.hpp"
typedef mc::FFGraph<> t_DAG;
typedef mc::AEBND<I,PM,PV> t_AEBND;

  
int main()
{
 try{
  t_DAG NLE;  // DAG describing the problem

  const unsigned NP = 1;  // Parameter dimension
  const unsigned NX = 3;  // State dimension

  mc::FFVar P[NP];  // Parameters p
  for( unsigned i=0; i<NP; i++ ) P[i].set( &NLE );

  mc::FFVar X[NX];  // Dependents x
  for( unsigned i=0; i<NX; i++ ) X[i].set( &NLE );

  mc::FFVar F[NX];  // Equations f(x,p)=0
  //F[0] = (P[0]+3.)*X[0] + (P[0]+2.)*X[1] - 3.*P[0];
  F[0] = (P[0]+3.)*X[0] + P[0]*X[1] - 3.*P[0];
  F[1] = P[0]*X[0] + (P[0]+2.)*X[1];
  F[2] = (P[0]+2.)*X[0] + (P[0]+3.)*X[1] - X[2];
  //F[2] = (P[0]+3.)*X[1] - X[2];

  I Ip[NP]  = { I(0e0,1e0) },
    Ix0[NX] = { I(-1e1,1e1), I(-1e1,1e1), I(-2e1,2e1) },
    Ix[NX];

  PM PMEnv( NP, NPM );
  //PMEnv.options.CENTER_REMAINDER = true;
  //PMEnv.options.REF_MIDPOINT = true;
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
  BND.options.INTERBND  = true; //false;
  BND.options.MAXIT     = 20;
  BND.options.RTOL      =
  BND.options.ATOL      = 1e-7;
  BND.options.BOUNDER   = t_AEBND::Options::ALGORITHM::AUTO;//GE;//KRAW;//GS;
  BND.options.PRECOND   = t_AEBND::Options::PRECONDITIONING::INVMD;//QRMD;//NONE;
  BND.options.AUTODIFF  = t_AEBND::Options::DIFFERENTIATION::ASA;//FSA;
  BND.setup();

  BND.options.BLKDEC    = t_AEBND::Options::DECOMPOSITION::NONE;
  BND.setup();
  std::cout << "\nDECOMPOSITION: NONE\n";
  BND.solve( Ip, Ix, Ix0 );
  BND.solve( PMp, PMx, Ix );

  BND.options.BLKDEC    = t_AEBND::Options::DECOMPOSITION::DIAG;
  BND.setup();
  std::cout << "\nDECOMPOSITION: DIAGONAL BLOCKS\n";
  BND.solve( Ip, Ix, Ix0 );
  BND.solve( PMp, PMx, Ix );

  BND.options.BLKDEC    = t_AEBND::Options::DECOMPOSITION::RECUR;
  BND.setup();
  std::cout << "\nDECOMPOSITION: RECURSIVE BLOCKS\n";
  BND.solve( Ip, Ix, Ix0 );
  BND.solve( PMp, PMx, Ix );
 }

 catch( t_AEBND::Exceptions &eObj ){
  std::cerr << "Error " << eObj.ierr()
            << eObj.what() << std::endl;
 }

 return 0;
}

