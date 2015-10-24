#define SAVE_RESULTS    // whether or not to save results to file
#undef USE_PROFIL	// specify to use PROFIL for interval arithmetic
#undef USE_FILIB	// specify to use FILIB++ for interval arithmetic
#undef DEBUG            // whether to output debug information
////////////////////////////////////////////////////////////////////////////////

#include <fstream>
#include "nlegpe.hpp"

#ifdef USE_PROFIL
  #include "mcprofil.hpp"
  typedef INTERVAL I;
#else
  #ifdef USE_FILIB
    #include "mcfilib.hpp"
    typedef filib::interval<double> I;
  #else
    #include "interval.hpp"
    typedef mc::Interval I;
  #endif
#endif

////////////////////////////////////////////////////////////////////////////////
// Find all solutions of the system of nonlinear equations:
//   (1.-R) * (D/10/(1+b1)-x) * exp(10*x/(1+10*x/g)) - x = 0
//   x - (1+b2)*y + (1.-R)*(D/10-b1*x-(1+b2)*y)*exp(10*y/(1+10*y/g)) = 0
// for (x,y) in [0,1]^2, R in [0.93,0.99], and g=1e3, D=22, b1=2, b2=2
////////////////////////////////////////////////////////////////////////////////

int main() 
{
  mc::FFGraph DAG;
  const unsigned NP = 3, NY = 2;
  mc::FFVar T, P[NP], Y[NY];
  T.set( &DAG );
  for( unsigned int i=0; i<NP; i++ ) P[i].set( &DAG );
  double g=1e3, D=22., b1=2., b2=2.;
  mc::FFVar R = P[2];
  Y[0] = (1.-R)*(D/10/(1+b1)-P[0])*mc::exp(10*P[0]/(1+10*P[0]/g))-P[0];
  Y[1] = P[0]-(1+b2)*P[1]+(1.-R)*(D/10-b1*P[0]-(1+b2)*P[1])*mc::exp(10*P[1]/(1+10*P[1]/g));

  mc::NLEGPE<I> problem;
  problem.set_dag( &DAG );
  problem.set_indep( &T );
  problem.set_par( NP, P );
  problem.set_dep( NY, Y );

  problem.options.SETINV.DISPLAY = 1;
  problem.options.SETINV.MAX_NODES = 2000;
  problem.options.SETINV.ABSOLUTE_TOLERANCE = 1e-6;
  problem.options.SETINV.RELATIVE_TOLERANCE = 1e-6;
  problem.options.SETINV.BRANCHING_VARIABLE_CRITERION = mc::SetInv<I>::Options::RGREL;
  problem.options.SETINV.MEASURE = mc::SetInv<I>::Options::LENGTH;

  problem.options.OUTPUT_BOUND = mc::NLEGPE<I>::Options::CM;
  problem.options.CM_ORDER     = 1;
  problem.options.OUTRED_MAX   = 10;
  problem.options.OUTRED_THRES = 2e-2;
  problem.options.OUTRED_TOL   = 1e-9;
  problem.options.INRED_MAX    = 0;
  problem.options.INRED_THRES  = 1e-2;

  const I Ip[NP] = { I(0.,1.), I(0.,1.), I(0.93,0.99) };

  const double TOL = 1e-6;
  typename std::list<mc::NLEGPE<I>::Data> Iym;
  for( unsigned int k=0; k<NY; k++ )
    Iym.push_back( typename mc::NLEGPE<I>::Data( -TOL, TOL, k ) );

  problem.solve( Ip, Iym, std::cout );

#if defined(SAVE_RESULTS )
  ofstream K_un( "undetermined.out", ios_base::out );
  problem.output_nodes( K_un ); //, true );
  K_un.close();
#endif

  return 0;
}
