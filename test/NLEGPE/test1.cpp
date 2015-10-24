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
//   4.731e-3*x1*x3 - 0.3578*x2*x3 - 0.1238*x1 + x7 - 1.637e-3*x2 - 0.9338*x4 = 0.3571
//   0.2238*x1*x3   + 0.7623*x2*x3 + 0.2638*x1 - x7 - 0.07745*x2  - 0.6734*x4 = 0.6022;
//   x6*x8 + 0.3578*x1 + 4.731e-3*x2 = 0;
//   -0.7623*x1 + 0.2238*x2 + 0.3461 = 0;
//   x1^2 + x2^2 = 1;
//   x3^2 + x4^2 = 1;
//   x5^2 + x6^2 = 1;
//   x7^2 + x8^2 = 1;
// for (x1,...,x8) in [-1,1]^8
////////////////////////////////////////////////////////////////////////////////

int main() 
{
  mc::FFGraph DAG;
  const unsigned NP = 8, NY = 8;
  mc::FFVar T, P[NP], Y[NY];
  T.set( &DAG );
  for( unsigned int i=0; i<NP; i++ ) P[i].set( &DAG );
  Y[0] = 4.731e-3*P[0]*P[2]-0.3578*P[1]*P[2]-0.1238*P[0]+P[6]-1.637e-3*P[1]-0.9338*P[3]-0.3571;
  Y[1] = 0.2238*P[0]*P[2]+0.7623*P[1]*P[2]+0.2638*P[0]-P[6]-0.07745*P[1]-0.6734*P[3]-0.6022;
  Y[2] = P[5]*P[7]+0.3578*P[0]+4.731e-3*P[1];
  Y[3] = -0.7623*P[0]+0.2238*P[1]+0.3461;
  Y[4] = sqr(P[0])+sqr(P[1])-1.;
  Y[5] = sqr(P[2])+sqr(P[3])-1.;
  Y[6] = sqr(P[4])+sqr(P[5])-1.;
  Y[7] = sqr(P[6])+sqr(P[7])-1.;

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
  problem.options.CM_ORDER     = 2;
  problem.options.OUTRED_MAX   = 10;
  problem.options.OUTRED_THRES = 2e-2;
  problem.options.OUTRED_TOL   = 1e-9;
  problem.options.INRED_MAX    = 0;
  problem.options.INRED_THRES  = 1e-2;

  const I Ip[NP] = { I(-1.,1.), I(-1.,1.), I(-1.,1.), I(-1.,1.),
                     I(-1.,1.), I(-1.,1.), I(-1.,1.), I(-1.,1.) };

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
