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

int main()
{
  mc::FFGraph NLE;  // DAG describing the problem

  const unsigned NP = 4;  // Parameter dimension [phiB, AB, BB, CB, AT, BT, CT]
  const unsigned NX = 4;  // State dimension [T, phiT, psiB, psiT]

  mc::FFVar P[NP];  // Parameters p
  for( unsigned i=0; i<NP; i++ ) P[i].set( &NLE );

  mc::FFVar X[NX];  // Dependents x
  for( unsigned i=0; i<NX; i++ ) X[i].set( &NLE );

  mc::FFVar F[NX];  // Equations f(x,p)=0
  double AB=6.87987, BB=1196.76, CB=219.161, AT=6.95087, BT=1342.31, CT=219.187,
         P0=759.81, phiB0=0.5;
  F[0] = phiB0 + X[1] - 1.;
  F[1] = X[2] + X[3] - 1.;
  mc::FFVar PBsat = pow(10.,AB*P[0]-BB*P[2]/(X[0]+CB));
  F[2] = P0*X[2] - PBsat*phiB0;
  mc::FFVar PTsat = pow(10.,AT*P[1]-BT*P[3]/(X[0]+CT));
  F[3] = P0*X[3] - PTsat*X[1];

  I Ip[NP]  = { I(0.998,1.002), I(0.998,1.002), I(0.998,1.002), I(0.998,1.002) }, //, I(1.), I(1.) },
  //I Ip[NP]  = { I(1.), I(0.998,1.002), I(0.998,1.002), I(0.998,1.002), I(0.998,1.002), I(0.998,1.002), I(0.998,1.002) },
  //I Ip[NP]  = { I(1.), I(0.998,1.002), I(0.998,1.002), I(1.), I(0.998,1.002), I(0.998,1.002), I(1.) },
    Ix0[NX] = { I(80.,100.), I(0.45,0.55), I(0.5,.9), I(0.1,0.5) },
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

  BND.options.DISPLAY = 2;
  BND.options.MAXIT   = 20;
  BND.options.RTOL     =
  BND.options.ATOL     = 1e-10;
  BND.options.BOUNDER  = mc::AEBND<I,PM,PV>::Options::ALGORITHM::AUTO;//GE;//KRAW;//GS;
  BND.options.PRECOND  = mc::AEBND<I,PM,PV>::Options::PRECONDITIONING::INVMD;//QRM;//NONE;
  BND.options.BLKDEC   = mc::AEBND<I,PM,PV>::Options::DECOMPOSITION::RECUR;//NONE;//DIAG;

  BND.setup();
  //std::cout << "\nSuccessful? " << (BND.solve( Ip, Ix, Ix0 )==mc::AEBND<I,PM,PV>::NORMAL?"Y\n":"N\n");
  std::cout << "\nSuccessful? " << (BND.solve( PMp, PMx, Ix0 )==mc::AEBND<I,PM,PV>::NORMAL?"Y\n":"N\n");

  return 0;
}


