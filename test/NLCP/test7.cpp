#define SAVE_RESULTS    // whether or not to save results to file
#define USE_PROFIL	// specify to use PROFIL for interval arithmetic
#undef  USE_FILIB	// specify to use FILIB++ for interval arithmetic

#define MC__USE_CPLEX   // whether to use CPLEX or GUROBI
#define MC__CSEARCH_SHOW_BOXES
#undef  MC__CSEARCH_SHOW_REDUC
#undef  MC__CSEARCH_SHOW_INCLUSION
#undef  MC__CSEARCH_PAUSE_INFEASIBLE

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
// Find all Fritz-John points of the nonlinear program:
//   min  x1 + x2
//   s.t. x1·x2 <= 4
//        0 <= x1 <= 6
//        1 <= x2 <= 4
////////////////////////////////////////////////////////////////////////////////

int main()
{  
  mc::FFGraph DAG;
  const unsigned NP = 2; mc::FFVar P[NP];
  for( unsigned i=0; i<NP; i++ ) P[i].set( &DAG );

  mc::NLCP<I> NLP;
  NLP.set_dag( &DAG );                       // DAG
  NLP.set_var( NP, P );                      // decision variables
  NLP.set_obj( mc::BASE_NLP::MAX, P[0]+P[1] );   // objective
  //NLP.set_obj( mc::BASE_NLP::MAX, sqr(P[0]+P[1]-1.) );   // objective
  NLP.add_ctr( mc::BASE_NLP::LE, P[0]*P[1]-4. ); // constraints

  NLP.options.MIPFILE     = ""; //"test7.lp";
  NLP.options.DISPLAY     = 2;
  NLP.options.MAXITER     = 1000;
  NLP.options.BLKDECUSE   = true; //false;
  NLP.options.NODEMEAS    = mc::SBP<I>::Options::RELMAXLEN;
  NLP.options.CVTOL       = 1e-2;
  NLP.options.FEASTOL     = 1e-9;
  NLP.options.BRANCHVAR   = mc::SBP<I>::Options::RGREL;
  NLP.options.VARMEAS     = mc::NLCP<I>::Options::ALL;//INDEP;//ALL;
  NLP.options.STGBCHDEPTH = 0;
  NLP.options.STGBCHRTOL  = 1e-2;
  NLP.options.DOMREDMAX   = 10;
  NLP.options.DOMREDTHRES = 5e-2;
  NLP.options.RELMETH     = mc::NLCP<I>::Options::CHEB;
  NLP.options.NCOCUTS     = false;
  NLP.options.CMODPROP    = 2;
  NLP.options.CMODDEPS    = 0;
  NLP.options.AEBND.DISPLAY = 0;
  std::cout << NLP;

  // Describe the feasible set
  I Ip[NP] = { I(0.,6.), I(1.,4.) };
  NLP.setup( Ip );
  NLP.solve( Ip );
  NLP.stats.display();

#if defined(SAVE_RESULTS )
  std::ofstream K_bnd( "test7_bnd.out", std::ios_base::out );
  std::ofstream K_inn( "test7_inn.out", std::ios_base::out );
  NLP.output_nodes( K_bnd, K_inn ); //, true );
  K_bnd.close();
  K_inn.close();
#endif

  // Enclose all the Fritz-John points
  NLP.options.MAXITER     = 10;
  NLP.options.NCOCUTS     = true;
  NLP.options.MIPDISPLAY  = 0;
  //NLP.options.POLIMG.AGGREG_LIN = false;
  NLP.setup( Ip );
  //unsigned nred;
  //I Ip1[NP] = { I(2.,3.), I(1.,2.) };
  //NLP.contract( Ip1, nred );
  NLP.solve( Ip );
  NLP.stats.display();

  auto L_clus = NLP.clusters();
  std::cout << "No clusters: " << L_clus.size() << std::endl;
  for( auto it=L_clus.begin(); it!=L_clus.end(); ++it ) delete[] *it;

#if defined(SAVE_RESULTS )
  std::ofstream K_clus( "test7_clus.out", std::ios_base::out );
  NLP.output_clusters( K_clus ); //, true );
  K_clus.close();
#endif

  return 0;
}

