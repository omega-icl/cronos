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
//   4*x^3 + 4*x*y + 2*y^2 - 42*x = 14
//   4*y^3 + 2*x^2 + 4*x*y - 26*y = 22
// for (x,y) in [-5,5] x [-5x5]
////////////////////////////////////////////////////////////////////////////////

int main() 
{
  mc::FFGraph DAG;
  const unsigned NP = 2, NY = 2;
  mc::FFVar T, P[NP], Y[NY];
  T.set( &DAG );
  for( unsigned int i=0; i<NP; i++ ) P[i].set( &DAG );
  Y[0] = 4*mc::pow(P[0],3)+4*P[0]*P[1]+2*mc::pow(P[1],2)-42*P[0]-14;
  Y[1] = 4*mc::pow(P[1],3)+2*mc::pow(P[0],2)+4*P[0]*P[1]-26*P[1]-22;

  mc::NLEGPE<I> problem;
  problem.set_dag( &DAG );
  problem.set_indep( &T );
  problem.set_par( NP, P );
  problem.set_dep( NY, Y );

  problem.options.SETINV.DISPLAY = 2;
  problem.options.SETINV.MAX_NODES = 100;
  problem.options.SETINV.ABSOLUTE_TOLERANCE = 1e-6;
  problem.options.SETINV.RELATIVE_TOLERANCE = 1e-6;
  problem.options.SETINV.BRANCHING_VARIABLE_CRITERION = mc::SetInv<I>::Options::RGABS;
  problem.options.SETINV.MEASURE = mc::SetInv<I>::Options::LENGTH;

  problem.options.OUTPUT_BOUND = mc::NLEGPE<I>::Options::CM;
  problem.options.CM_ORDER     = 3;
  problem.options.OUTRED_MAX   = 10;
  problem.options.OUTRED_THRES = 2e-2;
  problem.options.OUTRED_TOL   = 1e-9;
  problem.options.INRED_MAX    = 0;
  problem.options.INRED_THRES  = 1e-2;

  const I Ip[NP] = { I(-5.,5.), I(-5.,5.) };

  const double TOL = 1e-6;
  typename std::list<mc::NLEGPE<I>::Data> Iym;
  for( unsigned int k=0; k<NY; k++ )
    Iym.push_back( typename mc::NLEGPE<I>::Data( -TOL, TOL, k ) );

  problem.solve( Ip, Iym, std::cout );

  std::list<I*> L_clus = problem.clusters();
  std::cout << "No clusters: " << L_clus.size() << std::endl;
  for( auto it=L_clus.begin(); it!=L_clus.end(); ++it ) delete[] *it;

#if defined(SAVE_RESULTS )
  ofstream K_un( "undetermined.out", ios_base::out );
  problem.output_nodes( K_un ); //, true );
  K_un.close();

  ofstream K_cl( "clusters.out", ios_base::out );
  problem.output_clusters( K_cl ); //, true );
  K_cl.close();
#endif

  return 0;
}
