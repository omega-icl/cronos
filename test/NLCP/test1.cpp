#define SAVE_RESULTS    // whether or not to save results to file
#define USE_PROFIL	    // specify to use PROFIL for interval arithmetic
#undef USE_FILIB	    // specify to use FILIB++ for interval arithmetic
#undef DEBUG            // whether to output debug information
#define MC__USE_CPLEX   // whether to use CPLEX or GUROBI
#undef MC__CSEARCH_SHOW_BOXES
#undef MC__CSEARCH_SHOW_REDUC
////////////////////////////////////////////////////////////////////////////////

#include <fstream>
#include "nlcp.hpp"

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
  const unsigned NP = 8, NC = 8;
  mc::FFVar P[NP], C[NC];
  for( unsigned int i=0; i<NP; i++ ) P[i].set( &DAG );
  C[0] = 4.731e-3*P[0]*P[2]-0.3578*P[1]*P[2]-0.1238*P[0]+P[6]-1.637e-3*P[1]-0.9338*P[3]-0.3571;
  C[1] = 0.2238*P[0]*P[2]+0.7623*P[1]*P[2]+0.2638*P[0]-P[6]-0.07745*P[1]-0.6734*P[3]-0.6022;
  C[2] = P[5]*P[7]+0.3578*P[0]+4.731e-3*P[1];
  C[3] = -0.7623*P[0]+0.2238*P[1]+0.3461;
  C[4] = sqr(P[0])+sqr(P[1])-1.;
  C[5] = sqr(P[2])+sqr(P[3])-1.;
  C[6] = sqr(P[4])+sqr(P[5])-1.;
  C[7] = sqr(P[6])+sqr(P[7])-1.;

  mc::NLCP<I> CP;
  CP.set_dag( &DAG );
  CP.set_var( NP, P );
  for( unsigned ic=0; ic<NC; ic++ )
    CP.add_ctr( mc::BASE_OPT::EQ, C[ic] );

  //CP.options.MIPFILE     = "test1.lp";
  CP.options.DISPLAY     = 2;
  CP.options.MAXITER     = 100;
  CP.options.CVTOL       = 1e-6;
  CP.options.BRANCHVAR   = mc::SBP<I>::Options::RGREL;
  CP.options.NODEMEAS    = mc::SBP<I>::Options::RELMEANLEN;
  CP.options.DOMREDMAX   = 10;
  CP.options.DOMREDTHRES = 1e-1;
  CP.options.DOMREDBKOFF = 1e-7;
  CP.options.RELMETH     = mc::NLCP<I>::Options::CHEB;
  CP.options.CMODPROP    = 2;
  CP.options.STGBCHDEPTH = 0;
  CP.options.STGBCHRTOL  = 1e-2;
  CP.options.BLKDECUSE   = true;
  std::cout << CP;

  const I Ip[NP] = { I(-1.,1.), I(-1.,1.), I(-1.,1.), I(-1.,1.),
                     I(-1.,1.), I(-1.,1.), I(-1.,1.), I(-1.,1.) };

  CP.setup();
  CP.solve( Ip );
  CP.stats.display();

  auto L_clus = CP.clusters();
  std::cout << "No clusters: " << L_clus.size() << std::endl;
  for( auto it=L_clus.begin(); it!=L_clus.end(); ++it ) delete[] *it;

#if defined(SAVE_RESULTS )
  std::ofstream K_un( "test1.out", std::ios_base::out );
  CP.output_nodes( K_un ); //, true );
  K_un.close();

  std::ofstream K_cl( "test1b.out", std::ios_base::out );
  CP.output_clusters( K_cl ); //, true );
  K_cl.close();
#endif

  return 0;
}
