#define SAVE_RESULTS		// <- Whether to save bounds to file
#define USE_CMODEL		// <- whether to use Chebyshev models or Taylor models
const unsigned NPM   = 3;	// <- Order of polynomial expansion
const unsigned NSAMP = 10;	// <- Number of sampling points for inner approx.
const double   REDUC = 1.4;	// <- Reduction ratio for convergence analysis

#include "odebnd_sundials.hpp"
//#include "odebnd_val.hpp"

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

  double T0 = 0., TF = 2000.;   // Time span
  const unsigned NS = 1;  // Time stages
  double TS[NS+1]; TS[0] = T0;
  for( unsigned k=0; k<NS; k++ ) TS[k+1] = TS[k] + (TF-T0)/(double)NS;

  const unsigned NP = 2;  // Number of parameters
  const unsigned NX = 2;  // Number of states
  const unsigned NQ = 0;  // Number of state quadratures
  const unsigned NF = 1;  // Number of state functions

  const double k0 = 0.022;
  const double CA0 = 10.;
  const double V = 0.1;
  const double Cp = 60.;
  const double Ea = 6000.;
  const double R = 8.314;
  const double DHR = -140000.;
  const double UA = 3.;

  mc::FFVar P[NP];  // Parameter array
  for( unsigned i=0; i<NP; i++ ) P[i].set( &IVP );

  mc::FFVar X[NX];  // State array
  for( unsigned i=0; i<NX; i++ ) X[i].set( &IVP );

  mc::FFVar IC[NX];  // Initial value function
  IC[0] = 0.;
  IC[1] = P[0];

  mc::FFVar RHS[NX];  // RHS of differential equations
  RHS[0] = k0*exp(-Ea/(R*X[1]))*(1.-X[0]);
  RHS[1] = (P[1]-X[1])*UA/(CA0*V*Cp) - (1.-X[0])*exp(-Ea/(R*X[1]))*(DHR*k0/Cp);

  mc::FFVar FCT[NF];  // State functionals
  FCT[0] = X[0];

  I Ip[NP] = { I(350.,370.), I(290.,310.) };

  PM PMEnv( NP, NPM );
  PV PMp[NP];
  for( unsigned i=0; i<NP; i++ ) PMp[i].set( &PMEnv, i, Ip[i] );

  /////////////////////////////////////////////////////////////////////////
  // Bound ODE trajectories - differential inequalities

  mc::ODEBND_SUNDIALS<I,PM,PV> CSTR;

  CSTR.set_dag( &IVP );
  CSTR.set_time( NS, TS );
  CSTR.set_state( NX, X );
  CSTR.set_parameter( NP, P );
  CSTR.set_differential( NX, RHS );
  CSTR.set_initial( NX, IC );
  CSTR.set_function( NF, FCT );

#if defined( SAVE_RESULTS )
  CSTR.options.RESRECORD = 1000;
#endif
  CSTR.options.INTMETH   = mc::ODEBND_SUNDIALS<I,PM,PV>::Options::MSADAMS;
  CSTR.options.JACAPPROX = mc::ODEBND_SUNDIALS<I,PM,PV>::Options::CV_DIAG;//CV_DENSE;//
  CSTR.options.WRAPMIT   = mc::ODEBND_SUNDIALS<I,PM,PV>::Options::ELLIPS;//DINEQ;//NONE;
  CSTR.options.DISPLAY   = 1;
  CSTR.options.ORDMIT    = -2; //NPM;
  CSTR.options.ATOL      = 1e-10;
  CSTR.options.RTOL      = 1e-10;
  CSTR.options.ETOL      = 1e-20;
  CSTR.options.HMIN      = 1e-10;
  CSTR.options.ODESLVS   = CSTR.options;

  std::cout << "\nNON_VALIDATED INTEGRATION - INNER-APPROXIMATION OF REACHABLE SET:\n\n";
  CSTR.bounds( NSAMP, Ip );
#if defined( SAVE_RESULTS )
  std::ofstream apprec( "test4_APPROX_STA.dat", std::ios_base::out );
  CSTR.record( apprec );
#endif

  std::cout << "\nCONTINUOUS SET-VALUED INTEGRATION - INTERVAL ENCLOSURE OF REACHABLE SET:\n\n";
  CSTR.bounds( Ip );
#if defined( SAVE_RESULTS )
  std::ofstream bnd2recI( "test4_DINEQI_STA.dat", std::ios_base::out );
  CSTR.record( bnd2recI );
#endif

  std::cout << "\nCONTINUOUS SET-VALUED INTEGRATION - POLYNOMIAL MODEL ENCLOSURE OF REACHABLE SET:\n\n";
  //CSTR.options.PMNOREM = false;
  //CSTR.options.DMAX    = 5.;
  CSTR.bounds( PMp );
#if defined( SAVE_RESULTS )
  std::ofstream bnd2recPM( "test4_DINEQPM_STA.dat", std::ios_base::out );
  CSTR.record( bnd2recPM );
#endif

  return 0;
}


