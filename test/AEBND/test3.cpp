const unsigned int NPM   = 5;	// <- Order of Taylor model
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

  const unsigned NP = 2;  // Parameter dimension
  const unsigned NX = 1;  // State dimension

  mc::FFVar P[NP];  // Parameters p
  for( unsigned i=0; i<NP; i++ ) P[i].set( &NLE );

  mc::FFVar X[NX];  // Dependents x
  for( unsigned i=0; i<NX; i++ ) X[i].set( &NLE );

  mc::FFVar F[NX];  // Equations f(x,p)=0
  F[0] = X[0] - sin(P[0]) / sqrt(X[0]) - 25. * P[1];

  //I Ip[NP]  = { I(0.5,4.0), I(1.,1.) },
  I Ip[NP]  = { I(0.5,0.6), I(1.,1.) },
    Ix0[NX] = { I(2e1,3e1) },
    Ix[NX];

  PM PMEnv( NP, NPM );
  PV PMp[NP], PMx[NX], PMx0[NX];
  for( unsigned i=0; i<NP; i++ ) PMp[i].set( &PMEnv, i, Ip[i] );
  for( unsigned i=0; i<NX; i++ ) PMx0[i] = Ix0[i];

  /////////////////////////////////////////////////////////////////////////
  // Bound AE solution set
  mc::AEBND<I,PM,PV> BND;

  BND.set_dag( &NLE );
  BND.set_var( NP, P );
  BND.set_dep( NX, X, F );

  BND.options.DISPLAY = 1;
  BND.options.MAXIT   = 20;
  BND.options.RTOL     =
  BND.options.ATOL     = 1e-10;
  BND.options.BOUNDER  = mc::AEBND<I,PM,PV>::Options::ALGORITHM::AUTO;//GE;//KRAW;//GS;
  BND.options.PRECOND  = mc::AEBND<I,PM,PV>::Options::PRECONDITIONING::INVMD;//QRM;//NONE;
  BND.options.BLKDEC   = mc::AEBND<I,PM,PV>::Options::DECOMPOSITION::RECUR;//NONE;//DIAG;

  BND.setup();
  std::cout << "\nSuccessful? " << (BND.solve( Ip, Ix, Ix0 )==mc::AEBND<I,PM,PV>::NORMAL?"Y\n":"N\n");
  std::cout << "\nSuccessful? " << (BND.solve( PMp, PMx, Ix )==mc::AEBND<I,PM,PV>::NORMAL?"Y\n":"N\n");

  return 0;
}


