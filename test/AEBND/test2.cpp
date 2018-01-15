const unsigned int NPM   = 2;	// <- Order of Taylor model
#define USE_CMODEL		// <- Use Chebyshev models?
#undef  MC__AEBND_SHOW_PRECONDITIONING

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
  mc::FFGraph NLE;  // DAG describing the problem

  const unsigned NP = 1;  // Parameter dimension
  const unsigned NT = 11;  // State dimension

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
  mc::AEBND<I,PM,PV> BND;

  BND.set_dag( &NLE );
  BND.set_var( NP, P );
  BND.set_dep( NT, T, F );

  BND.options.DISPLAY = 1;
  BND.options.MAXIT   = 20;
  BND.options.RTOL    =
  BND.options.ATOL    = 1e-10;
  BND.options.BOUNDER  = mc::AEBND<I,PM,PV>::Options::ALGORITHM::GS;//GE;//KRAW;//GS;
  BND.options.PRECOND  = mc::AEBND<I,PM,PV>::Options::PRECONDITIONING::INVMD;//QRM;//NONE;
  BND.options.BLKDEC   = mc::AEBND<I,PM,PV>::Options::DECOMPOSITION::RECUR;//NONE;//DIAG;

  BND.setup();
  std::cout << "\nSuccessful? " << (BND.solve( Ip, It, It0 )==mc::AEBND<I,PM,PV>::NORMAL?"Y\n":"N\n");
  std::cout << "\nSuccessful? " << (BND.solve( PMp, PMt, It )==mc::AEBND<I,PM,PV>::NORMAL?"Y\n":"N\n");

  return 0;
}


