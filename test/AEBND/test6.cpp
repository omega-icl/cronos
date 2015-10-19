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

const unsigned NP = 1;   // Parameter dimension
const unsigned NX = 2;   // State dimension
const unsigned NS = 250;  // Time stages
const double   DT = 0.002; // Time step

template <class T>
void Residual
( const T*xkp1, const T*xk, const T*p, T*res )
{
  //res[0] = xkp1[0] - xk[0] - DT*p[0]*xk[0]*(1.-xk[1]);
  //res[1] = xkp1[1] - xk[1] - DT*p[0]*xk[1]*(xk[0]-1.);
  res[0] = xkp1[0] - xk[0] - DT*p[0]*xkp1[0]*(1.-xkp1[1]);
  res[1] = xkp1[1] - xk[1] - DT*p[0]*xkp1[1]*(xkp1[0]-1.);
}

template <class T>
void Initial
( const T*p, const T*x0, T*res )
{
  res[0] = x0[0] - 1.2;
  res[1] = x0[1] - 1.1;
}

int main()
{
  mc::FFGraph IVP;  // DAG describing the problem

  mc::FFVar P[NP];  // Parameters p
  for( unsigned i=0; i<NP; i++ ) P[i].set( &IVP );

  mc::FFVar X[NX+NX*NS];  // States x
  for( unsigned i=0; i<NX+NX*NS; i++ ) X[i].set( &IVP );

  mc::FFVar F[NX+NX*NS];  // IC+IVP residuals
  Initial( P, X, F );
  for( unsigned k=0, ires=NX; k<NS; k++, ires+=NX )
    Residual( X+ires, X+ires-2, P, F+ires );
/*
  std::ofstream o_ivp( "test6_DAG.dot", std::ios_base::out );
  IVP.dot_script( NX+NX*NS, F, o_ivp );
  o_ivp.close();
*/

  /////////////////////////////////////////////////////////////////////////
  // Bound AE solution set
  mc::AEBND<I,PM,PV> EX6;

  EX6.options.DISPLAY = 1;
  EX6.options.MAXIT   = 30;
  EX6.options.RTOL    =
  EX6.options.ATOL    = 1e-10;
  EX6.options.BOUNDER = mc::AEBND<I,PM,PV>::Options::GS;//AUTO;
  EX6.options.PRECOND = mc::AEBND<I,PM,PV>::Options::INVMB;//INVMB;//INVBD;//NONE;

  EX6.set_dag( &IVP );
  EX6.set_var( NP, P );
  EX6.set_dep( NX+NX*NS, X, F );

  I Ip[NP]  = { I(2.95,3.05) },
    Ix[NX+NX*NS], Ix0[NX+NX*NS];
  for( unsigned k=0; k<NX+NX*NS; k++) Ix0[k] = I(0e0,2e0);

  EX6.options.BLKDEC  = true;
  EX6.setup();
  std::cout << "\nSuccessful? " << (EX6.solve( Ip, Ix, Ix0 )==mc::AEBND<I,PM,PV>::NORMAL?"Y\n":"N\n");

  std::ofstream o_I0sol( "test6_I0.out", std::ios_base::out );
  o_I0sol << std::scientific << std::setprecision(5);
  double t=0;
  for( unsigned k=0, px=0; k<=NS; k++, t+=DT, px+=NX ){
    o_I0sol << t;
    for( unsigned i=0; i<NX; i++ )
      o_I0sol << "  " << mc::Op<I>::l(Ix[px+i]) << "  " << mc::Op<I>::u(Ix[px+i]);
    o_I0sol << std::endl;
  }
  o_I0sol.close();

  EX6.options.BLKDEC  = false;
  EX6.setup();
  std::cout << "\nSuccessful? " << (EX6.solve( Ip, Ix, Ix )==mc::AEBND<I,PM,PV>::NORMAL?"Y\n":"N\n");

  std::ofstream o_Isol( "test6_I.out", std::ios_base::out );
  o_Isol << std::scientific << std::setprecision(5);
  t=0;
  for( unsigned k=0, px=0; k<=NS; k++, t+=DT, px+=NX ){
    o_Isol << t;
    for( unsigned i=0; i<NX; i++ )
      o_Isol << "  " << mc::Op<I>::l(Ix[px+i]) << "  " << mc::Op<I>::u(Ix[px+i]);
    o_Isol << std::endl;
  }
  o_Isol.close();

  PM PMEnv( NP, NPM );
  PV PMp[NP], PMx[NX+NX*NS];
  for( unsigned i=0; i<NP; i++ ) PMp[i].set( &PMEnv, i, Ip[i] );

  std::cout << "\nSuccessful? " << (EX6.solve( PMp, PMx, Ix )==mc::AEBND<I,PM,PV>::NORMAL?"Y\n":"N\n");

  std::ofstream o_PMsol( "test6_PM.out", std::ios_base::out );
  o_PMsol << std::scientific << std::setprecision(5);
  t=0;
  for( unsigned k=0, px=0; k<=NS; k++, t+=DT, px+=NX ){
    o_PMsol << t;
    for( unsigned i=0; i<NX; i++ )
      o_PMsol << "  " << mc::Op<PV>::l(PMx[px+i]) << "  " << mc::Op<PV>::u(PMx[px+i]);
    o_PMsol << std::endl;
  }
  o_PMsol.close();

  return 0;
}


