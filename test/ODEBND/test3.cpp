const unsigned NPM   = 8;	// <- Order of polynomial expansion
const unsigned NSAMP = 10;	// <- Number of sampling points for inner approx.
const double   REDUC = 1.4;	// <- Reduction ratio for convergence analysis
#define SAVE_RESULTS		// <- Whether to save bounds to file
#undef  TEST_CONVERGENCE	// <- Whether to test Hausdorff convergence of bounds
#define USE_CMODEL		// <- whether to use Chebyshev models or Taylor models

#include "odebnd_sundials.hpp"

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
  mc::FFGraph IVP;  // DAG describing the problem

  double T0 = 0., TF = 3.;   // Time span
  const unsigned NS = 120;    // Time stages
  double TS[NS+1]; TS[0] = T0;
  for( unsigned k=0; k<NS; k++ ) TS[k+1] = TS[k] + (TF-T0)/(double)NS;

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
    FCT[k*NF+1] = X[3] * (TF-T0)/(double)NS;
  }
  FCT[(NS-1)*NF+0] = Q[0];

  I Ip[NP]  = { I(0.95,1.05) };
  I *Ixk[NS+1], Iq[NQ], If[NF];
  for( unsigned k=0; k<=NS; k++ )
    Ixk[k] = new I[NX];

  //double p[NP]  = { 1. };
  //double *xk[NS+1], q[NQ], f[NF];
  //for( unsigned k=0; k<=NS; k++ )
  //  xk[k] = new double[NX];

  PM PMEnv( NP, NPM );
  PV PMp[NP];
  for( unsigned i=0; i<NP; i++ ) PMp[i].set( &PMEnv, i, Ip[i] );
  PV *PMxk[NS+1], PMq[NQ], PMf[NF];
  double *Hxk[NS+1];
  for( unsigned k=0; k<=NS; k++ ){
    PMxk[k] = new PV[NX];
    Hxk[k] = new double[NX];
  }

  /////////////////////////////////////////////////////////////////////////
  // Bound ODE trajectories - differential inequalities

  mc::ODEBND_SUNDIALS<I,PM,PV> PEND;

  PEND.set_dag( &IVP );
  PEND.set_time( NS, TS );
  PEND.set_state( NX, X );
  PEND.set_parameter( NP, P );
  PEND.set_differential( NX, RHS );
  PEND.set_initial( NX, IC );
  PEND.set_quadrature( NQ, QUAD, Q );
  PEND.set_function( NS, NF, FCT );

  PEND.options.INTMETH   = mc::ODEBND_SUNDIALS<I,PM,PV>::Options::MSADAMS;
  PEND.options.JACAPPROX = mc::ODEBND_SUNDIALS<I,PM,PV>::Options::CV_DIAG;//CV_DENSE;//
  PEND.options.WRAPMIT   = mc::ODEBND_SUNDIALS<I,PM,PV>::Options::ELLIPS;//DINEQ;//NONE;
  PEND.options.DISPLAY   = 1;
#if defined( SAVE_RESULTS )
  PEND.options.RESRECORD = true;
#endif
  PEND.options.ORDMIT       = -2; //NPM;
  PEND.options.ATOL         = 1e-10;
  PEND.options.RTOL         = 1e-10;
  PEND.options.ETOL         = 1e-20;
  PEND.options.ODESLV.ATOL  = 1e-10;
  PEND.options.ODESLV.RTOL  = 1e-8;
  PEND.options.NMAX         = 10000;
  PEND.options.HMIN         = 1e-12;

  std::cout << "\nNON_VALIDATED INTEGRATION - INNER-APPROXIMATION OF REACHABLE SET:\n\n";
  PEND.bounds( Ip, Ixk, If, NSAMP );
#if defined( SAVE_RESULTS )
  std::ofstream apprec( "test3_APPROX_STA.dat", std::ios_base::out );
  PEND.record( apprec );
#endif
  { int dum; std::cout << "PAUSED--"; std::cin >> dum; }

  std::cout << "\nCONTINUOUS SET-VALUED INTEGRATION - INTERVAL ENCLOSURE OF REACHABLE SET:\n\n";
  PEND.bounds( Ip, Ixk, If );
#if defined( SAVE_RESULTS )
  std::ofstream bnd2recI( "test3_DINEQI_STA.dat", std::ios_base::out );
  PEND.record( bnd2recI );
#endif
  { int dum; std::cout << "PAUSED--"; std::cin >> dum; }

  std::cout << "\nCONTINUOUS SET-VALUED INTEGRATION - POLYNOMIAL MODEL ENCLOSURE OF REACHABLE SET:\n\n";
  PEND.options.PMNOREM = false;
  PEND.options.DMAX    = 1e2;
  PEND.bounds( PMp, PMxk, PMf );
#if defined( SAVE_RESULTS )
  std::ofstream bnd2recPM( "test3_DINEQPM_STA.dat", std::ios_base::out );
  PEND.record( bnd2recPM );
#endif

#if defined( TEST_CONVERGENCE )
  PEND.hausdorff( PMp, Hxk, 0, NSAMP );
  PEND.options.DISPLAY = 0;
  std::cout << std::scientific << std::setprecision(5);
  for( unsigned k=0; k<20; k++ ){
    PV PMp0[NP] = { PV( &PMEnv, 0, Ip[0]/std::pow(REDUC,k) ) };
    double Hxk0 = (k? Hxk[NS][0]: 0. );
    PEND.hausdorff( PMp0, Hxk, 0, NSAMP );
    std::cout << mc::Op<I>::diam( PMp0[0].B() ) << "  " << Hxk[NS][0]
              << "  " << std::log(Hxk0/Hxk[NS][0])/log(REDUC)
              << std::endl;
  }
#endif

  for( unsigned k=0; k<=NS; k++ ){
    delete[] Ixk[k];
    delete[] PMxk[k];
    delete[] Hxk[k];
  }

  return 0;
}

