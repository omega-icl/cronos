#define SAVE_RESULTS		// <- Whether to save bounds to file
#undef  TEST_CONVERGENCE	// <- Whether to test Hausdorff convergence of bounds
#define USE_CMODEL		// <- whether to use Chebyshev models or Taylor models
#define USE_SUNDIALS		// <- whether to use SUNDIALS or GSL integrator
const unsigned NPM   = 3;	// <- Order of polynomial expansion
const unsigned NSAMP = 10;	// <- Number of sampling points for inner approx.
const double   REDUC = 1.4;	// <- Reduction ratio for convergence analysis

#include "odeslv_gsl.hpp"
#ifndef USE_SUNDIALS
  #include "odebnd_gsl.hpp"
#else
  #include "odebnd_sundials.hpp"
#endif
#include "odebnd_val.hpp"

#include "interval.hpp"
typedef mc::Interval I;
typedef mc::Ellipsoid E;

#ifdef USE_CMODEL
  #include "cmodel.hpp"
  typedef mc::CModel<I> TM;
  typedef mc::CVar<I> TV;
#else
  #include "tmodel.hpp"
  typedef mc::TModel<I> TM;
  typedef mc::TVar<I> TV;
#endif

int main()
{
  mc::FFGraph IVP;  // DAG describing the problem

  double t0 = 0., tf = 4000.;   // Time span
  const unsigned NS = 100;  // Time stages
  double tk[NS+1]; tk[0] = t0;
  for( unsigned k=0; k<NS; k++ ) tk[k+1] = tk[k] + (tf-t0)/(double)NS;

  const unsigned NP = 2;  // Number of parameters
  const unsigned NX = 2;  // Number of states
  const unsigned NQ = 0;  // Numboer of state quadratures
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
  I *Ixk[NS+1], Iq[NQ], If[NF];
  for( unsigned k=0; k<=NS; k++ )
    Ixk[k] = new I[NX];

  TM TMEnv( NP, NPM );
  TV TMp[NP];
  for( unsigned i=0; i<NP; i++ ) TMp[i].set( &TMEnv, i, Ip[i] );
  TV *TMxk[NS+1], TMq[NQ], TMf[NF];
  double *Hxk[NS+1];
  for( unsigned k=0; k<=NS; k++ ){
    TMxk[k] = new TV[NX];
    Hxk[k] = new double[NX];
  }

  /////////////////////////////////////////////////////////////////////////
  // Bound ODE trajectories - sampling
  mc::ODESLV_GSL<I> LV0;

  LV0.set_dag( &IVP );
  LV0.set_state( NX, X );
  LV0.set_parameter( NP, P );
  LV0.set_differential( NX, RHS );
  LV0.set_initial( NX, IC );
  LV0.set_function( NF, FCT );

  LV0.options.DISPLAY = 1;
  LV0.options.ATOL    = LV0.options.RTOL = 1e-10;
  LV0.options.H0      = 1e-6;
  //LV0.options.INPMETH = mc::ODESLV_GSL<I>::Options::MSBDF;
#if defined( SAVE_RESULTS )
  LV0.options.RESRECORD = true;
#endif

  std::cout << "\nNON_VALIDATED INTEGRATION - APPROXIMATE ENCLOSURE OF REACHABLE SET:\n\n";
  //LV0.states( NS, tk, p, xk, q, f );
  LV0.bounds( NS, tk, Ip, Ixk, Iq, If, NSAMP );
#if defined( SAVE_RESULTS )
  std::ofstream apprec( "test4_APPROX_STA.dat", std::ios_base::out );
  LV0.record( apprec );
#endif

  /////////////////////////////////////////////////////////////////////////
  // Bound ODE trajectories - differential inequalities
#ifndef USE_SUNDIALS // GSL integrator
  mc::ODEBND_GSL<I,TM,TV> LV2;

  LV2.set_dag( &IVP );
  LV2.set_state( NX, X );
  LV2.set_parameter( NP, P );
  LV2.set_differential( NX, RHS );
  LV2.set_initial( NX, IC );
  LV2.set_function( NF, FCT );

  LV2.options.DISPLAY   = 1;
#if defined( SAVE_RESULTS )
  LV2.options.RESRECORD = true;
#endif
  LV2.options.ORDMIT    = 1; //TMp->nord()+1;
  LV2.options.WRAPMIT   = mc::ODEBND_GSL<I,TM,TV>::Options::ELLIPS;//NONE;

#else // SUNDIALS integrator
  mc::ODEBND_SUNDIALS<I,TM,TV> LV2;

  LV2.set_dag( &IVP );
  LV2.set_state( NX, X );
  LV2.set_parameter( NP, P );
  LV2.set_differential( NX, RHS );
  LV2.set_initial( NX, IC );
  LV2.set_function( NF, FCT );

  LV2.options.DISPLAY   = 1;
#if defined( SAVE_RESULTS )
  LV2.options.RESRECORD = true;
#endif
  LV2.options.ORDMIT    = 1; //TMp->nord();
  //LV2.options.HMIN      = 1e-10;
  //LV2.options.ATOL      = LV2.options.RTOL = 1e-10;
  //LV2.options.INTMETH   = mc::ODEBND_SUNDIALS<I,TM,TV>::Options::MSBDF;
  //LV2.options.JACAPPROX = mc::ODEBND_SUNDIALS<I,TM,TV>::Options::CV_DIAG;
  LV2.options.WRAPMIT   = mc::ODEBND_SUNDIALS<I,TM,TV>::Options::ELLIPS;//DINEQ;//NONE;
#endif

  std::cout << "\nCONTINUOUS SET-VALUED INTEGRATION - INTERVAL ENCLOSURE OF REACHABLE SET:\n\n";
  LV2.options.DISPLAY = 1;
  LV2.bounds( NS, tk, Ip, Ixk, Iq, If );
#if defined( SAVE_RESULTS )
  std::ofstream bndrecIA( "test4_DINEQI_STA.dat", std::ios_base::out );
  LV2.record( bndrecIA );
#endif
#if defined( TEST_CONVERGENCE )
  LV2.hausdorff( NS, tk, Ip, Hxk, LV0, NSAMP );
  LV2.options.DISPLAY = 0;
  std::cout << std::scientific << std::setprecision(5);
  for( unsigned k=0; k<20; k++ ){
    I Ip0[NP] = { Ip[0]/std::pow(REDUC,k) };
    double Hxk0 = (k? Hxk[NS][0]: 0. );
    LV2.hausdorff( NS, tk, Ip0, Hxk, LV0, NSAMP );
    std::cout << mc::Op<I>::diam( Ip0[0] ) << "  " << Hxk[NS][0]
              << "  " << std::log(Hxk0/Hxk[NS][0])/log(REDUC)
              << std::endl;
  }
#endif

  std::cout << "\nCONTINUOUS SET-VALUED INTEGRATION - POLYNOMIAL MODEL ENCLOSURE OF REACHABLE SET:\n\n";
  LV2.options.DISPLAY = 1;
  LV2.bounds( NS, tk, TMp, TMxk, TMq, TMf );
#if defined( SAVE_RESULTS )
  std::ofstream bndrecPM( "test4_DINEQPM_STA.dat", std::ios_base::out );
  LV2.record( bndrecPM );
#endif
#if defined( TEST_CONVERGENCE )
  LV2.hausdorff( NS, tk, TMp, Hxk, LV0, NSAMP );
  LV2.options.DISPLAY = 0;
  std::cout << std::scientific << std::setprecision(5);
  for( unsigned k=0; k<20; k++ ){
    TV TMp0[NP] = { TV( &TMEnv, 0, Ip[0]/std::pow(REDUC,k) ) };
    double Hxk0 = (k? Hxk[NS][0]: 0. );
    LV2.hausdorff( NS, tk, TMp0, Hxk, LV0, NSAMP );
    std::cout << mc::Op<I>::diam( TMp0[0].B() ) << "  " << Hxk[NS][0]
              << "  " << std::log(Hxk0/Hxk[NS][0])/log(REDUC)
              << std::endl;
  }
#endif

  for( unsigned k=0; k<=NS; k++ ){
    delete[] Ixk[k];
    delete[] TMxk[k];
    delete[] Hxk[k];
  }

  return 0;
}


