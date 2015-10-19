const unsigned int NPM   = 2;	// <- Order of Taylor model
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

using namespace CPPL;

int main()
{
/*dgematrix A(4,4);
  A(0,0) = 1.5; A(0,1) = 1.;  A(0,2) = 0.;  A(0,3) = 0.; 
  A(1,0) = 1.;  A(1,1) = 0.5; A(1,2) = 1.;  A(1,3) = 0.; 
  A(2,0) = 0.;  A(2,1) = 1.;  A(2,2) = 0.5; A(2,3) = 1.; 
  A(3,0) = 0.;  A(3,1) = 0.;  A(3,2) = 1.;  A(3,3) = 1.5; 

  dgematrix A(3,3);
  A(0,0) = 1.;     A(0,1) = 1./2.; A(0,2) = 1./3.;
  A(1,0) = 1./2.;  A(1,1) = 1./3.; A(1,2) = 1./4.;
  A(2,0) = 1./3.;  A(2,1) = 1./4.; A(2,2) = 1./5.; 

  dgematrix Q, R;
  return dgeqrf( A, Q, R ); */

  mc::FFGraph NLE;  // DAG describing the problem

  const unsigned NP = 1;  // Parameter dimension
  const unsigned NT = 101;  // State dimension

  mc::FFVar P[NP];  // Parameters p
  for( unsigned i=0; i<NP; i++ ) P[i].set( &NLE );

  mc::FFVar T[NT];  // Dependents x
  for( unsigned i=0; i<NT; i++ ) T[i].set( &NLE );

  mc::FFVar F[NT];  // Equations f(x,p)=0
  F[0] = T[0] - 500.;
  for( unsigned k=1; k<NT-1; k++){
    if( k < 0.5*(NT-1.) || k > 0.6*(NT-1.) )
      F[k] = mc::sqr(NT-1.) * ( T[k-1] - 2.*T[k] + T[k+1] ) - ( T[k] + 5e3 ) / P[0];
    else
      F[k] = mc::sqr(NT-1.) * ( T[k-1] - 2.*T[k] + T[k+1] ) - ( T[k] - 15e3 ) / P[0];
  }
  F[NT-1] = T[NT-1] - 600.;

  I Ip[NP]  = { I(9e-1,1e0) },
    It[NT], It0[NT];
  for( unsigned k=0; k<NT; k++) It0[k] = I(1e2,1e3);

  PM PMEnv( NP, NPM );
  PV PMp[NP], PMt[NT];
  for( unsigned i=0; i<NP; i++ ) PMp[i].set( &PMEnv, i, Ip[i] );

  /////////////////////////////////////////////////////////////////////////
  // Bound AE solution set
  mc::AEBND<I,PM,PV> EX1;

  EX1.set_dag( &NLE );
  EX1.set_var( NP, P );
  EX1.set_dep( NT, T, F );

  EX1.options.DISPLAY = 1;
  EX1.options.BLKDEC  = false;//true;
  EX1.options.MAXIT   = 20;
  EX1.options.RTOL    =
  EX1.options.ATOL    = 1e-10;
  EX1.options.BOUNDER = mc::AEBND<I,PM,PV>::Options::GS;//AUTO;
  EX1.options.PRECOND = mc::AEBND<I,PM,PV>::Options::INVMB;//INVMB;//INVBD;//NONE;

  EX1.setup();
  std::cout << "\nSuccessful? " << (EX1.solve( Ip, It, It0 )==mc::AEBND<I,PM,PV>::NORMAL?"Y\n":"N\n");
  std::cout << "\nSuccessful? " << (EX1.solve( PMp, PMt, It )==mc::AEBND<I,PM,PV>::NORMAL?"Y\n":"N\n");

  return 0;
}


