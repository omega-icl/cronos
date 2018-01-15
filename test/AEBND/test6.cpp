const unsigned int NPM   = 5;	// <- Order of Taylor/Chebyshev model
#define USE_CMODEL		// <- Use Chebyshev models?
#undef  USE_IMPLICIT    // <- Whether to use explicit or implicitEuler formula?
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

  const unsigned NP = 2;  // Parameter dimension
  const unsigned NX = 200;  // State dimension

  mc::FFVar P[NP];  // Parameters p
  for( unsigned i=0; i<NP; i++ ) P[i].set( &NLE );

  mc::FFVar X[NX], SOL[NX];  // Dependents x
  for( unsigned i=0; i<NX; i++ ) X[i].set( &NLE );

  double h = 0.05;
  mc::FFVar F[NX];  // Equations f(x,p)=0
  F[0] = X[0] - P[0];
  F[1] = X[1] - (1.+P[1]);
  //F[1] = X[1] - P[1];
  for( unsigned i=0; i<NX/2-1; i++ ){
#ifndef USE_IMPLICIT
    F[2+2*i] = X[2+2*i] - X[0+2*i] - h*X[1+2*i];
    F[3+2*i] = X[3+2*i] - X[1+2*i] + h*X[0+2*i];
#else
    F[2+2*i] = X[2+2*i] - X[0+2*i] - h*X[3+2*i];
    F[3+2*i] = X[3+2*i] - X[1+2*i] + h*X[2+2*i];
#endif
  }

  I Ip[NP] = { 5e-2*I(-1e0,1e0), 5e-2*I(-1e0,1e0) }, Ix0[NX], Ix1[NX], Ix2[NX], Ix3[NX];
  for( unsigned i=0; i<NX; i++ ) Ix0[i] = I(-1e1,1e1);
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

  BND.options.DISPLAY  = 1;
  BND.options.INTERBND = true; //false;
  BND.options.MAXIT    = 20;
  BND.options.RTOL     =
  BND.options.ATOL     = 1e-8;
  BND.options.BOUNDER  = mc::AEBND<I,PM,PV>::Options::ALGORITHM::GS;//KRAW;//GE;
  BND.options.PRECOND  = mc::AEBND<I,PM,PV>::Options::PRECONDITIONING::INVMB;//INVMD;//NONE;
  BND.setup();

  BND.options.BLKDEC = mc::AEBND<I,PM,PV>::Options::DECOMPOSITION::NONE;
  std::cout << "\nSuccessful? " << (BND.solve( Ip, Ix1, Ix0 )==mc::AEBND<I,PM,PV>::NORMAL?"Y\n":"N\n");

  BND.options.BLKDEC = mc::AEBND<I,PM,PV>::Options::DECOMPOSITION::DIAG;
  std::cout << "\nSuccessful? " << (BND.solve( Ip, Ix2, Ix0 )==mc::AEBND<I,PM,PV>::NORMAL?"Y\n":"N\n");

  BND.options.BLKDEC = mc::AEBND<I,PM,PV>::Options::DECOMPOSITION::RECUR;
  std::cout << "\nSuccessful? " << (BND.solve( Ip, Ix3, Ix0 )==mc::AEBND<I,PM,PV>::NORMAL?"Y\n":"N\n");

  std::ofstream ofile( "test6.out" );
  for( unsigned i=0; i<NX; i+=2 )
    ofile << i/2+1 << "  "
          << mc::Op<I>::l(Ix1[i]) << "  "  << mc::Op<I>::u(Ix1[i]) << "  "
          << mc::Op<I>::l(Ix1[i+1]) << "  "  << mc::Op<I>::u(Ix1[i+1]) << "  "
          << mc::Op<I>::l(Ix2[i]) << "  "  << mc::Op<I>::u(Ix2[i]) << "  " 
          << mc::Op<I>::l(Ix2[i+1]) << "  "  << mc::Op<I>::u(Ix2[i+1]) << "  " 
          << mc::Op<I>::l(Ix3[i]) << "  "  << mc::Op<I>::u(Ix3[i]) << "  " 
          << mc::Op<I>::l(Ix3[i+1]) << "  "  << mc::Op<I>::u(Ix3[i+1]) << std::endl;
  ofile.close();

  return 0;


/*
  std::cout << "\nSuccessful? " << (BND.solve( PMp, PMx, Ix )==mc::AEBND<I,PM,PV>::NORMAL?"Y\n":"N\n");
  std::cout << "\nPMx2: " << (3+PMp[0])*PMx[1];
  return 0;

  std::cout << "\nSuccessful? " << (BND.solve( SOL )==mc::AEBND<I,PM,PV>::NORMAL?"Y\n":"N\n");
  std::ofstream o_sol( "test6_DAG.dot", std::ios_base::out );
  NLE.dot_script( NX, SOL, o_sol );
  o_sol.close();
  NLE.eval( NX, SOL, PMx, NP, P, PMp );
  for( unsigned i=0; i<NX; i++ )
    std::cout << " X" << i << ": " << PMx[i] << std::endl;
  NLE.eval( NX, SOL, Ix, NP, P, Ip );
  for( unsigned i=0; i<NX; i++ )
    std::cout << " X" << i << ": " << Ix[i] << std::endl;

  return 0;
*/
}


