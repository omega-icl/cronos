const unsigned int NPM   = 4;	// <- Order of Taylor model
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

const unsigned NP = 1;    // Parameter dimension
const unsigned NX = 2;    // State dimension
const double   DT = 1e-1; // Time step
const unsigned NT = 7;    // Time expansion order
const unsigned NH = 40;   // Length of receeding horizon
const unsigned DH = 40;   // Shift of receeding horizon
const unsigned NS = 6;   // Number of shifts
//const unsigned NH = 100;   // Length of receeding horizon
//const unsigned DH = 10;   // Shift of receeding horizon
//const unsigned NS = 1;    // Number of shifts

mc::FFGraph IVP;  // DAG describing the problem

void RightHandSide
( const mc::FFVar*x, const mc::FFVar*p, mc::FFVar*rhs )
{
  rhs[0] = p[0]*x[0]*(1.-x[1]);
  rhs[1] = p[0]*x[1]*(x[0]-1.);
}

void RightHandSide
( const mc::FFVar*x, const mc::FFVar*p, const double h, mc::FFVar*phi )
{
  mc::FFVar rhs[NX]; RightHandSide( x, p, rhs );
  const mc::FFVar* TRHS = IVP.TAD( NT, NX, rhs, NX, x );
  unsigned pT = NT*NX;
  for( unsigned i=0; i<NX; i++ ) phi[i] = TRHS[pT+i];
  pT -= NX;
  for( unsigned k=0; k<NT; k++, pT-=NX )
    for( unsigned i=0; i<NX; i++ ){ phi[i] *= h; phi[i] += TRHS[pT+i]; }
  delete[] TRHS;
}

void Residual
( const mc::FFVar*xkp1, const mc::FFVar*xk, const mc::FFVar*p, const double h, mc::FFVar*res )
{
  mc::FFVar rhs[NX]; RightHandSide( xk, p, h, rhs );
  res[0] = xkp1[0] - rhs[0];
  res[1] = xkp1[1] - rhs[1];
  //mc::FFVar rhs[NX]; RightHandSide( xkp1, p, rhs );
  //mc::FFVar rhs[NX]; RightHandSide( xkp1, p, -h, rhs );
  //res[0] = xk[0] - rhs[0];
  //res[1] = xk[1] - rhs[1];
}

template <class T>
void Initial
( const T*xk, const T*x0, T*res )
{
  res[0] = xk[0] - x0[0];
  res[1] = xk[1] - x0[1];
}

int main()
{
  mc::FFVar PX0[NP+NX];  // Parameters p and initial states x0
  for( unsigned i=0; i<NP+NX; i++ ) PX0[i].set( &IVP );

  mc::FFVar X[NX+NX*NH];  // States x
  for( unsigned i=0; i<NX+NX*NH; i++ ) X[i].set( &IVP );

  mc::FFVar F[NX+NX*NH];  // ODE residuals
  Initial( X, PX0+NP, F );
  for( unsigned k=0, ires=NX; k<NH; k++, ires+=NX )
    Residual( X+ires, X+ires-NX, PX0, DT, F+ires );

  IVP.output( IVP.subgraph( NX, F+NX ) );

  std::ofstream o_ivp( "test7_DAG.dot", std::ios_base::out );
  IVP.dot_script( NX+NX*NH, F, o_ivp );
  o_ivp.close();

  /////////////////////////////////////////////////////////////////////////
  // Bound AE solution set
  mc::AEBND<I,PM,PV> EX7;

  EX7.options.DISPLAY = 1;
  EX7.options.MAXIT   = 100;
  EX7.options.RTOL    =
  EX7.options.ATOL    = 1e-10;
  EX7.options.BOUNDER = mc::AEBND<I,PM,PV>::Options::KRAW;//AUTO;
  EX7.options.PRECOND = mc::AEBND<I,PM,PV>::Options::INVMB;//INVMB;//INVBD;//NONE;

  EX7.set_dag( &IVP );
  EX7.set_var( NP+NX, PX0 );
  EX7.set_dep( NX+NX*NH, X, F );
  EX7.options.BLKDEC  = false;
  EX7.setup();

  //H.set( &IVP, DT );
  I Ipx0[NP+NX]  = { I(2.95,3.05), I(1.2), I(1.1) },
    Ix[NX*(1+NH+(NS-1)*DH)];
  for( unsigned i=0; i<NX; i++ ) Ix[i] = Ipx0[NP+i];
  for( unsigned k=NX; k<NX*(1+NH+(NS-1)*DH); k++) Ix[k] = I(.5e0,1.5e0);

  for( unsigned k=0, ires=0; k<NS; k++, ires+=NX*DH ){
    std::cout << "Stages " << k*DH << "-" << k*DH+NH-1 << ": "
              << (EX7.solve( Ipx0, Ix+ires, Ix+ires )==mc::AEBND<I,PM,PV>::NORMAL?"OK\n":"FAIL\n");
    for( unsigned i=0; i<NX; i++ ) Ipx0[NP+i] = Ix[ires+NX*DH+i];
  }

  std::ofstream o_Isol( "test7_I.out", std::ios_base::out );
  o_Isol << std::scientific << std::setprecision(5);
  double t=0;
  for( unsigned k=0, px=0; k<1+NH+(NS-1)*DH; k++, t+=DT, px+=NX ){
    o_Isol << t;
    for( unsigned i=0; i<NX; i++ )
      o_Isol << "  " << mc::Op<I>::l(Ix[px+i]) << "  " << mc::Op<I>::u(Ix[px+i]);
    o_Isol << std::endl;
  }
  o_Isol.close();

  PM PMEnv( NP, NPM );
  PV PMpx0[NP+NX] = { PV(&PMEnv, 0, Ipx0[0]), PV(1.2), PV(1.1) },
     PMx[NX*(1+NH+(NS-1)*DH)];

  for( unsigned k=0, ires=0; k<NS; k++, ires+=NX*DH ){
    std::cout << "Stages " << k*DH << "-" << k*DH+NH-1 << ": "
              << (EX7.solve( PMpx0, PMx+ires, Ix+ires )==mc::AEBND<I,PM,PV>::NORMAL?"OK\n":"FAIL\n");
    for( unsigned i=0; i<NX; i++ ) PMpx0[NP+i] = PMx[ires+NX*DH+i];
  }
/*
  EX7.options.BLKDEC  = false;
  EX7.options.DISPLAY = 2;
  EX7.setup();
  for( unsigned k=0, ires=0; k<NS; k++, ires+=NX*DH ){
    std::cout << "Stages " << k*DH << "-" << k*DH+NH-1 << ": "
              << (EX7.solve( PMpx0, PMx+ires, PMx+ires )==mc::AEBND<I,PM,PV>::NORMAL?"OK\n":"FAIL\n");
    for( unsigned i=0; i<NX; i++ ) PMpx0[NP+i] = PMx[ires+NX*DH+i];
  }
*/
  std::ofstream o_PMsol( "test7_PM.out", std::ios_base::out );
  o_PMsol << std::scientific << std::setprecision(5);
  t=0;
  for( unsigned k=0, px=0; k<1+NH+(NS-1)*DH; k++, t+=DT, px+=NX ){
    o_PMsol << t;
    for( unsigned i=0; i<NX; i++ )
      o_PMsol << "  " << mc::Op<PV>::l(PMx[px+i]) << "  " << mc::Op<PV>::u(PMx[px+i]);
    o_PMsol << std::endl;
  }
  o_PMsol.close();

  return 0;
}


