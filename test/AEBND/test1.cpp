const unsigned int NPM   = 5;	// <- Order of Taylor/Chebyshev model
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

int main()
{
  mc::FFGraph NLE;  // DAG describing the problem

  const unsigned NP = 1;  // Parameter dimension
  const unsigned NX = 3;  // State dimension

  mc::FFVar P[NP];  // Parameters p
  for( unsigned i=0; i<NP; i++ ) P[i].set( &NLE );

  mc::FFVar X[NX], SOL[NX];  // Dependents x
  for( unsigned i=0; i<NX; i++ ) X[i].set( &NLE );

  mc::FFVar F[NX];  // Equations f(x,p)=0
  //F[0] = (P[0]+3.)*X[1] - 3.*P[0];
  F[0] = (P[0]+3.)*X[0] + P[0]*X[1] - 3.*P[0];
  F[1] = P[0]*X[0] + (P[0]+3.)*X[1];
  F[2] = (P[0]+3.)*X[2] - X[1];// + (P[0]+3.)*X[1];

  I Ip[NP]  = { I(-1e0,1e0) },
    //Ix0[NX] = { I(-1e1,1e1), I(-1e1,1e1) },
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
  mc::AEBND<I,PM,PV> BND;

  BND.set_dag( &NLE );
  BND.set_var( NP, P );
  BND.set_dep( NX, X, F );

  BND.options.DISPLAY  = 1;
  BND.options.INTERBND = true; //false;
  BND.options.MAXIT    = 20;
  BND.options.RTOL     =
  BND.options.ATOL     = 1e-8;
  BND.options.BOUNDER  = mc::AEBND<I,PM,PV>::Options::ALGORITHM::GS;//KRAW;//GE;
  BND.options.PRECOND  = mc::AEBND<I,PM,PV>::Options::PRECONDITIONING::INVMD;//QRM;//NONE;
  BND.setup();

  //BND.setup();
  //std::cout << "\nSuccessful? " << (BND.solve( PMp, PMx, Ix0 )==mc::AEBND<I,PM,PV>::NORMAL?"Y\n":"N\n");
  //std::cout << "\nPMx2: " << (3+PMp[0])*PMx[1];
  //return 0;

  BND.options.BLKDEC = mc::AEBND<I,PM,PV>::Options::DECOMPOSITION::RECUR;
  std::cout << "\nSuccessful? " << (BND.solve( Ip, Ix, Ix0 )==mc::AEBND<I,PM,PV>::NORMAL?"Y\n":"N\n");
  std::cout << "\nSuccessful? " << (BND.solve( PMp, PMx, Ix )==mc::AEBND<I,PM,PV>::NORMAL?"Y\n":"N\n");

  BND.options.BLKDEC = mc::AEBND<I,PM,PV>::Options::DECOMPOSITION::NONE;
  std::cout << "\nSuccessful? " << (BND.solve( Ip, Ix, Ix0 )==mc::AEBND<I,PM,PV>::NORMAL?"Y\n":"N\n");
  std::cout << "\nSuccessful? " << (BND.solve( PMp, PMx, Ix )==mc::AEBND<I,PM,PV>::NORMAL?"Y\n":"N\n");

  BND.options.BLKDEC = mc::AEBND<I,PM,PV>::Options::DECOMPOSITION::DIAG;
  std::cout << "\nSuccessful? " << (BND.solve( Ip, Ix, Ix0 )==mc::AEBND<I,PM,PV>::NORMAL?"Y\n":"N\n");
  std::cout << "\nSuccessful? " << (BND.solve( PMp, PMx, Ix )==mc::AEBND<I,PM,PV>::NORMAL?"Y\n":"N\n");

  std::cout << "\nPMx2: " << PMx[1]/(3+PMp[0]);
  return 0;

  std::cout << "\nSuccessful? " << (BND.solve( SOL )==mc::AEBND<I,PM,PV>::NORMAL?"Y\n":"N\n");
  std::ofstream o_sol( "test1_DAG.dot", std::ios_base::out );
  NLE.dot_script( NX, SOL, o_sol );
  o_sol.close();
  NLE.eval( NX, SOL, PMx, NP, P, PMp );
  for( unsigned i=0; i<NX; i++ )
    std::cout << " X" << i << ": " << PMx[i] << std::endl;
  NLE.eval( NX, SOL, Ix, NP, P, Ip );
  for( unsigned i=0; i<NX; i++ )
    std::cout << " X" << i << ": " << Ix[i] << std::endl;

  return 0;
}


