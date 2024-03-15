const unsigned int NPM   = 5;	// <- Order of Taylor model
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

  const unsigned NP = 2;  // Parameter dimension
  const unsigned NX = 1;  // State dimension

  mc::FFVar P[NP];  // Parameters p
  for( unsigned i=0; i<NP; i++ ) P[i].set( &NLE );

  mc::FFVar X[NX];  // Dependents x
  for( unsigned i=0; i<NX; i++ ) X[i].set( &NLE );

  mc::FFVar F[NX];  // Equations f(x,p)=0
  F[0] = X[0] - sin(P[0]) / sqrt(X[0]) - 25. * P[1];

  I Ip[NP]  = { I(0.5,2.), I(1.,2.) },
    Ix0[NX] = { I(20.,60.) },
    Ix[NX];

  PM PMEnv( NP, NPM );
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
  BND.options.BLKDEC    = t_AEBND::Options::DECOMPOSITION::NONE;//DIAG;

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


