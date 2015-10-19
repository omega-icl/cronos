#define SAVE_RESULTS		// <- Whether to save bounds to file
#undef  TEST_CONVERGENCE	// <- Whether to test Hausdorff convergence of bounds
#define USE_CMODEL		// <- whether to use Chebyshev models or Taylor models
const unsigned NPM   = 3;	// <- Order of polynomial expansion
const unsigned NSAMP = 10;	// <- Number of sampling points for inner approx.
const double   REDUC = 1.4;	// <- Reduction ratio for convergence analysis

#include "odeslv_gsl.hpp"
#include "odebnd_gsl.hpp"

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

  double t0 = 0., tf = .3;   // Time span
  const unsigned NS = 10;  // Time stages
  double tk[NS+1]; tk[0] = t0;
  for( unsigned k=0; k<NS; k++ ) tk[k+1] = tk[k] + (tf-t0)/(double)NS;

  const unsigned NP = 1;  // Parameter dimension
  const unsigned NX = 4;  // State dimension
  const unsigned NQ = 1;  // State dimension
  const unsigned NF = 2;  // State dimension

  const double g       = 9.81e0;	// m/s^2
  const double m1      = 1.e0;		// kg
  const double m2      = 1.e0;		// kg
  const double l1      = 1.e0;		// m
  const double l2      = 1.e0;		// m
  const double psi1_0  = 3.*mc::PI/4.;	// rad
  const double psi2_0  = -11.*mc::PI/20.;	// rad
  const double psi3_0  = .43e0;		// rad/s
  const double psi4_0  = .67e0;		// rad/s

  mc::FFVar P[NP];  // Parameter array
  for( unsigned i=0; i<NP; i++ ) P[i].set( &IVP );

  mc::FFVar X[NX];  // State array
  for( unsigned i=0; i<NX; i++ ) X[i].set( &IVP );

  mc::FFVar Q[NQ];  // Quadrature array
  for( unsigned i=0; i<NQ; i++ ) Q[i].set( &IVP );

  mc::FFVar IC[NX];  // Initial value function
  IC[0] = psi1_0 * P[0];
  IC[1] = psi2_0;
  IC[2] = psi3_0;
  IC[3] = psi4_0;

  mc::FFVar RHS[NX];  // RHS of differential equations
  mc::FFVar a11 = m1*l1 + m2*(l1+l2*cos(X[1]));
  mc::FFVar a21 = m2*(l1*cos(X[1])+l2);
  mc::FFVar a12 = m2*l2*cos(X[1]);
  mc::FFVar a22 = m2*l2;
  mc::FFVar b1  = -g*(m1+m2)*sin(X[0]) + m2*l2*sin(X[1])*((X[2]+X[3])*(X[2]+X[3]));
  mc::FFVar b2  = -g*m2*sin(X[0]+X[1]) - m2*l1*sin(X[1])*(X[2]*X[2]);
  RHS[0] = X[2];
  RHS[1] = X[3];
  RHS[2] = ( a22*b1 - a12*b2 ) / ( a11*a22 - a12*a21 );
  RHS[3] = ( a11*b2 - a21*b1 ) / ( a11*a22 - a12*a21 );

  mc::FFVar QUAD[NQ];  // Quadrature equations
  QUAD[0] = X[3];

  //mc::FFVar FCT[NF];  // State functionals
  //FCT[0] = Q[0];
  //FCT[1] = sqr(X[2]);
  mc::FFVar FCT[NF*NS];  // State functionals
  for( unsigned k=0; k<NS; k++ ){
    FCT[k*NF+0] = 0.;
    FCT[k*NF+1] = X[3] * (tf-t0)/(double)NS;
  }
  FCT[(NS-1)*NF+0] = Q[0];

  I Ip[NP]  = { I(0.9,1.1) };
  I *Ixk[NS+1], Iq[NQ], If[NF];
  for( unsigned k=0; k<=NS; k++ )
    Ixk[k] = new I[NX];

  //double p[NP]  = { 1. };
  //double *xk[NS+1], q[NQ], f[NF];
  //for( unsigned k=0; k<=NS; k++ )
  //  xk[k] = new double[NX];

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
  LV0.set_quadrature( NQ, QUAD, Q );
  LV0.set_function( NF, NS, FCT );

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
  std::ofstream apprec( "test3_APPROX_STA.dat", std::ios_base::out );
  LV0.record( apprec );
#endif

  /////////////////////////////////////////////////////////////////////////
  // Bound ODE trajectories - differential inequalities
#ifdef USE_CMODEL
  mc::ODEBND_GSL<I,TM,TV> LV2;
#else
  mc::ODEBND_GSL<I> LV2;
#endif
  LV2.set_dag( &IVP );
  LV2.set_state( NX, X );
  LV2.set_parameter( NP, P );
  LV2.set_differential( NX, RHS );
  LV2.set_initial( NX, IC );
  LV2.set_quadrature( NQ, QUAD, Q );
  LV2.set_function( NF, NS, FCT );

  LV2.options.H0        = 1e-4;
  //LV2.options.DMAX      = 1e4;
  LV2.options.NMAX      = 10000;
  LV2.options.ORDMIT    = 1; //TMp->nord();
  LV2.options.ATOL      = LV2.options.RTOL = 1e-10;
  LV2.options.WRAPMIT   = mc::ODEBND_GSL<I,TM,TV>::Options::ELLIPS;//DINEQ;//NONE;//
#if defined( SAVE_RESULTS )
  LV2.options.RESRECORD = true;
#endif

  std::cout << "\nCONTINUOUS SET-VALUED INTEGRATION - INTERVAL ENCLOSURE OF REACHABLE SET:\n\n";
  LV2.options.DISPLAY = 1;
  LV2.bounds( NS, tk, Ip, Ixk, Iq, If );
#if defined( SAVE_RESULTS )
  std::ofstream bndrecIA( "test3_DINEQI_STA.dat", std::ios_base::out );
  LV2.record( bndrecIA );
#endif
  LV2.hausdorff( NS, tk, Ip, Hxk, LV0, NSAMP );
#if defined( TEST_CONVERGENCE )
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
  std::ofstream bndrecPM( "test3_DINEQPM_STA.dat", std::ios_base::out );
  LV2.record( bndrecPM );
#endif
  LV2.hausdorff( NS, tk, TMp, Hxk, LV0, NSAMP );
#if defined( TEST_CONVERGENCE )
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


