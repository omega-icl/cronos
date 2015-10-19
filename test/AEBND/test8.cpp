const unsigned int NPM   = 10;	// <- Order of Taylor model
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

  const unsigned NP = 1;  // Parameter dimension
  const unsigned NX = 3;  // State dimension

  mc::FFVar P[NP];  // Parameters p
  for( unsigned i=0; i<NP; i++ ) P[i].set( &NLE );

  mc::FFVar X[NX];  // Dependents x
  for( unsigned i=0; i<NX; i++ ) X[i].set( &NLE );

  mc::FFVar F[NX];  // Equations f(x,p)=0
  F[0]=(0.5*(mc::exp(-X[1])-mc::exp(X[1]))*P[0]+0.5*(mc::exp(-X[1])+mc::exp(X[1]))*X[0])*(-0.5*(mc::exp(-X[1])-mc::exp(X[1]))*P[0]-0.5*(mc::exp(-X[1])+mc::exp(X[1]))*X[0])-(0.5*(mc::exp(-X[1])-mc::exp(X[1]))*P[0]+0.5*(mc::exp(-X[1])+mc::exp(X[1]))*X[0])*(0.5*(mc::exp(-X[1])+mc::exp(X[1]))*P[0]+0.5*(mc::exp(-X[1])-mc::exp(X[1]))*X[0]+1.2);

  F[1]=((mc::exp(-X[2]+X[1])-mc::exp(X[2]-X[1]))*(0.5*(mc::exp(-X[1])+mc::exp(X[1]))*P[0]+0.5*(mc::exp(-X[1])-mc::exp(X[1]))*X[0]) + mc::exp(-X[2]+X[1])*(0.5*(mc::exp(-X[1])-mc::exp(X[1]))*P[0]+0.5*(mc::exp(-X[1])+mc::exp(X[1]))*X[0]) - 1.2*mc::exp(X[2]-X[1]) +1.2)*(mc::exp(X[2]-X[1])*(0.5*(mc::exp(-X[1])+mc::exp(X[1]))*P[0]+0.5*(mc::exp(-X[1])-mc::exp(X[1]))*X[0]+1.2))-((mc::exp(-X[2]+X[1])-mc::exp(X[2]-X[1]))*(0.5*(mc::exp(-X[1])+mc::exp(X[1]))*P[0]+0.5*(mc::exp(-X[1])-mc::exp(X[1]))*X[0]) + mc::exp(-X[2]+X[1])*(0.5*(mc::exp(-X[1])-mc::exp(X[1]))*P[0]+0.5*(mc::exp(-X[1])+mc::exp(X[1]))*X[0]) - 1.2*mc::exp(X[1]-X[1]) +1.2)*(-mc::exp(X[2]-X[1])*(0.5*(mc::exp(-X[1])+mc::exp(X[1]))*P[0]+0.5*(mc::exp(-X[1])-mc::exp(X[1]))*X[0]+1.2)+1.2);

  F[2]= 0.5*(mc::exp(-2+X[2])+mc::exp(2-X[2]))*(mc::exp(X[2]-X[1])*(0.5*(mc::exp(-X[1])+mc::exp(X[1]))*P[0]+0.5*(mc::exp(-X[1])-mc::exp(X[1]))*X[0]+1.2)-1.2)+0.5*(mc::exp(-2+X[2])-mc::exp(2-X[2]))*((mc::exp(-X[2]+X[1])-mc::exp(X[2]-X[1]))*(0.5*(mc::exp(-X[1])+mc::exp(X[1]))*P[0]+0.5*(mc::exp(-X[1])-mc::exp(X[1]))*X[0]) + mc::exp(-X[2]+X[1])*(0.5*(mc::exp(-X[1])-mc::exp(X[1]))*P[0]+0.5*(mc::exp(-X[1])+mc::exp(X[1]))*X[0]) - 1.2*mc::exp(X[1]-X[1]) +1.2)-(0.5*(mc::exp(-2+X[2])-mc::exp(2-X[2]))*(mc::exp(X[2]-X[1])*(0.5*(mc::exp(-X[1])+mc::exp(X[1]))*P[0]+0.5*(mc::exp(-X[1])-mc::exp(X[1]))*X[0]+1.2)-1.2)+0.5*(mc::exp(-2+X[2])+mc::exp(2-X[2]))*((mc::exp(-X[2]+X[1])-mc::exp(X[2]-X[1]))*(0.5*(mc::exp(-X[1])+mc::exp(X[1]))*P[0]+0.5*(mc::exp(-X[1])-mc::exp(X[1]))*X[0]) + mc::exp(-X[2]+X[1])*(0.5*(mc::exp(-X[1])-mc::exp(X[1]))*P[0]+0.5*(mc::exp(-X[1])+mc::exp(X[1]))*X[0]) - 1.2*mc::exp(X[1]-X[1]) +1.2));

  I Ip[NP]  = { I(-0.94,-0.9399)},
    Ix0[NX] = { I(-0.5,0.5), I(0.,0.1), I(0.,0.1)},
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
  EX1.options.BLKDEC  = false;
  EX1.options.MAXIT   = 5;
  EX1.options.RTOL    =
  EX1.options.ATOL    = 1e-6;
  EX1.options.BOUNDER = mc::AEBND<I,PM,PV>::Options::GS;
  EX1.options.PRECOND = mc::AEBND<I,PM,PV>::Options::INVMD;//QRM;//NONE;
  try{
    EX1.setup();
    std::cout << "\nSuccessful? " << (EX1.solve( Ip, Ix, Ix0 )==mc::AEBND<I,PM,PV>::NORMAL?"Y\n":"N\n");
    std::cout << "\nSuccessful? " << (EX1.solve( PMp, PMx, Ix )==mc::AEBND<I,PM,PV>::NORMAL?"Y\n":"N\n");
  }
  catch(mc::AEBND<mc::Interval, mc::TModel<mc::Interval>, mc::TVar<mc::Interval> >::Exceptions e){
    std::cout << "Error number: " << e.ierr() << endl;
    std::cout << e.what() << endl;
  }

  return 0;
}




