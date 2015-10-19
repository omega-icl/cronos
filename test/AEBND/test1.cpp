const unsigned int NPM   = 15;	// <- Order of Taylor/Chebyshev model
//#define USE_CMODEL		// <- Use Chebyshev models?

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
  //F[0] = P[0]*X[1] - 3.*P[0];
  F[0] = (P[0]+3.)*X[0] + P[0]*X[1] - 3.*P[0];
  F[1] = P[0]*X[0] + (P[0]+3.)*X[1];
  F[2] = X[2] - (P[0]+3.)*X[1];

  I Ip[NP]  = { I(-1e0,1e0) },
    //Ix0[NX] = { I(-1e1,1e1), I(-1e1,1e1) },
    Ix0[NX] = { I(-1e1,1e1), I(-1e1,1e1), I(-1e1,1e1) },
    Ix[NX];

  PM PMEnv( NP, NPM );
  PV PMp[NP], PMx[NX], PMx0[NX];
  for( unsigned i=0; i<NP; i++ ) PMp[i].set( &PMEnv, i, Ip[i] );
  for( unsigned i=0; i<NX; i++ ) PMx0[i] = Ix0[i];

  /////////////////////////////////////////////////////////////////////////
  // Bound AE solution set
  mc::AEBND<I,PM,PV> EX1;

  EX1.set_dag( &NLE );
  EX1.set_var( NP, P );
  EX1.set_dep( NX, X, F );

  EX1.options.DISPLAY = 2;
  EX1.options.BLKDEC  = false;//true;
  EX1.options.MAXIT   = 20;
  EX1.options.RTOL    =
  EX1.options.ATOL    = 0e0;
  EX1.options.BOUNDER = mc::AEBND<I,PM,PV>::Options::GE;//KRAW;//GS;
  EX1.options.PRECOND = mc::AEBND<I,PM,PV>::Options::INVMD;//QRM;//NONE;

  EX1.setup();
  std::cout << "\nSuccessful? " << (EX1.solve( Ip, Ix, Ix0 )==mc::AEBND<I,PM,PV>::NORMAL?"Y\n":"N\n");
  std::cout << "\nSuccessful? " << (EX1.solve( PMp, PMx, Ix )==mc::AEBND<I,PM,PV>::NORMAL?"Y\n":"N\n");

  std::cout << "\nSuccessful? " << (EX1.solve( SOL )==mc::AEBND<I,PM,PV>::NORMAL?"Y\n":"N\n");
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


