#define SAVE_RESULTS    // whether or not to save results to file
#define USE_PROFIL	// specify to use PROFIL for interval arithmetic
#undef  USE_FILIB	// specify to use FILIB++ for interval arithmetic
#undef  DEBUG            // whether to output debug information
#define MC__USE_CPLEX   // whether to use CPLEX or GUROBI
#undef  MC__CSEARCH_SHOW_BOXES
#undef  MC__CSEARCH_SHOW_INCLUSION
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
// Find all solutions of the nonlinear equation:
//   x^2 + y^2 = 1
// for (x,y) in [-5,5] x [-5x5]
////////////////////////////////////////////////////////////////////////////////

int main()
{
  mc::FFGraph DAG;
  const unsigned NP = 2, NC = 1;
  mc::FFVar P[NP], C[NC];
  for( unsigned int i=0; i<NP; i++ ) P[i].set( &DAG );
  C[0] = mc::sqr(P[0])+mc::sqr(P[1])-1;

  mc::NLCP<I> CP;
  CP.set_dag( &DAG );
  CP.set_var( NP, P );
  CP.add_ctr( mc::BASE_OPT::EQ, C[0] );

  CP.options.MIPFILE     = "test0.lp";
  CP.options.DISPLAY     = 2;
  CP.options.MAXITER     = 100;
  CP.options.NODEMEAS    = mc::SBP<I>::Options::RELMAXLEN;
  CP.options.CVTOL       = 1e-2;
  CP.options.FEASTOL     = 1e-7;
  CP.options.BRANCHVAR   = mc::SBP<I>::Options::RGREL;
  CP.options.VARMEAS     = mc::NLCP<I>::Options::DEP;
  CP.options.STGBCHDEPTH = 0;
  CP.options.STGBCHRTOL  = 1e-2;
  CP.options.DOMREDMAX   = 10;
  CP.options.DOMREDTHRES = 1e-1;
  CP.options.DOMREDBKOFF = 1e-7;
  CP.options.RELMETH     = mc::NLCP<I>::Options::CHEB;
  CP.options.CMODPROP    = 2;
  CP.options.CMODDEPS    = 5;
  CP.options.CMODATOL    = 1e-5;
  CP.options.CMODRTOL    = 2e-2;
  //CP.options.CMODWARMS   = false;
  CP.options.CMODRED     = mc::NLCP<I>::Options::APPEND;
  //CP.options.AEBND.ATOL    = 
  //CP.options.AEBND.RTOL    = 1e-10;
  CP.options.AEBND.DISPLAY = 0;
  std::cout << CP;

  const I Ip[NP] = { I(-5.,5.), I(-5.,5.) };

  CP.setup( Ip );
  CP.solve( Ip );
  CP.stats.display();

  auto L_clus = CP.clusters();
  std::cout << "No clusters: " << L_clus.size() << std::endl;
  for( auto it=L_clus.begin(); it!=L_clus.end(); ++it ) delete[] *it;

#if defined(SAVE_RESULTS )
  std::ofstream K_un( "test0.out", std::ios_base::out );
  CP.output_nodes( K_un ); //, true );
  K_un.close();

  std::ofstream K_cl( "test0b.out", std::ios_base::out );
  CP.output_clusters( K_cl ); //, true );
  K_cl.close();
#endif

  return 0;
}
