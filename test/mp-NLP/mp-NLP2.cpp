#define SAVE_RESULTS    // whether or not to save results to file
#define USE_PROFIL	// specify to use PROFIL for interval arithmetic
#undef USE_FILIB	// specify to use FILIB++ for interval arithmetic
#undef DEBUG            // whether to output debug information
#define USE_DEPS	// whether to use dependents
#define MC__USE_CPLEX   // whether to use CPLEX or GUROBI
#define MC__NLCP_SHOW_BOXES
//#define MC__NLCP_SHOW_INCLUSIONS
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
typedef mc::CVar<I> CVI;

////////////////////////////////////////////////////////////////////////////////
// Find all solutions of the nonlinear constraints:
//   3x1^2 + 4x1 - 5  + 2y1 + 0.5y2  = 0
//                2x2 - 3  + y1 + y2 = 0
//                         2x1 + x2  = 2.5 + p1
//                       0.5x1 + x2  = 1.5 + p2
//                       0 <= p1,p2 <= 1
//                      x1,x2,y1,y2 >= 0
////////////////////////////////////////////////////////////////////////////////

int main() 
{
  mc::FFGraph DAG;
  const unsigned NP = 6, NC = 4;
  mc::FFVar P[NP], C[NC];
  for( unsigned int i=0; i<NP; i++ ) P[i].set( &DAG );

  mc::FFVar *CEQ = C, *X = P+2, *Y = P+4;
  CEQ[0] = 3*sqr(X[0]) + 4*X[0] - 5 + 2*Y[0] + 0.5*Y[1];
  CEQ[1] = 2*X[1] - 3 + Y[0] + Y[1];
  CEQ[2] = 2*X[0] + X[1] - P[0] - 2.5;
  CEQ[3] = .5*X[0] + X[1] - P[1] - 1.5;

  mc::NLCP<I> CP;
  CP.set_dag( &DAG );
#ifndef USE_DEPS
  CP.set_var( NP, P );
  CP.add_ctr( mc::BASE_OPT::EQ, CEQ[0] );
  CP.add_ctr( mc::BASE_OPT::EQ, CEQ[1] );
  CP.add_ctr( mc::BASE_OPT::EQ, CEQ[2] );
  CP.add_ctr( mc::BASE_OPT::EQ, CEQ[3] );
#else
  CP.set_var( 2, P );
  CP.set_dep( 4, X, CEQ );
#endif

  //CP.options.MIPFILE   = "test0.lp";
  CP.options.DISPLAY     = 2;
  CP.options.MAXITER     = 2;
  CP.options.CVATOL      = 1e-3;
  CP.options.CVRTOL      = 1e-3;
  CP.options.BRANCHVAR   = mc::SetInv<CVI>::Options::RGREL;
  CP.options.NODEMEAS    = mc::SetInv<CVI>::Options::MAXWIDTH;
  CP.options.STGBCHDEPTH = 0;
  CP.options.STGBCHDRMAX = 1;
  CP.options.STGBCHRTOL  = 1e-2;
  CP.options.DOMREDMAX   = 20;
  CP.options.DOMREDTHRES = 1e-4;
  CP.options.DOMREDBKOFF = 1e-5;
  CP.options.RELMETH     = mc::NLCP<I>::Options::CHEB;
  CP.options.CMODPROP    = 2;
  CP.options.CMODSPAR    = true;
  CP.options.CMREDORD    = 2;
  CP.options.CMREDTHRES  = 1e-3;
  CP.options.CMREDWARMS  = false;//true;
  std::cout << CP;

  const I IP[NP] = { I(0.,1.), I(0.,1.), I(0.,5.), I(0.,5.), I(0.,5.), I(0.,5.) };
  CP.setup();
  CP.solve( IP );

#if defined(SAVE_RESULTS )
  std::ofstream ores;
  ores.open( "mp-NLP2.out", std::ios_base::out );
  CP.output_nodes( ores, 7 );
  ores.close();
#endif

  return 0;
}
