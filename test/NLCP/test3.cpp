#define SAVE_RESULTS    // whether or not to save results to file
#undef  USE_PROFIL	// specify to use PROFIL for interval arithmetic
#undef  USE_FILIB	// specify to use FILIB++ for interval arithmetic
#undef  DEBUG		// whether to output debug information
#define USE_DEPS	// whether to use dependents
#undef MC__USE_CPLEX   // whether to use CPLEX or GUROBI
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
// Find all solutions of the system of nonlinear equations:
//   (1.-R) * (D/10/(1+b1)-x) * exp(10*x/(1+10*x/g)) - x = 0
//   x - (1+b2)*y + (1.-R)*(D/10-b1*x-(1+b2)*y)*exp(10*y/(1+10*y/g)) = 0
// for (x,y) in [0,1]^2, R in [0.93,0.99], and g=1e3, D=22, b1=2, b2=2
////////////////////////////////////////////////////////////////////////////////

int main() 
{
  mc::FFGraph DAG;
  const unsigned NP = 3, NC = 2;
  mc::FFVar P[NP], C[NC];
  for( unsigned int i=0; i<NP; i++ ) P[i].set( &DAG );
  double g=1e3, D=22., b1=2., b2=2.;
  mc::FFVar R = P[0];
  C[0] = (1.-R)*(D/10/(1+b1)-P[1])*mc::exp(10*P[1]/(1+10*P[1]/g))-P[1];
  C[1] = P[1]-(1+b2)*P[2]+(1.-R)*(D/10-b1*P[1]-(1+b2)*P[2])*mc::exp(10*P[2]/(1+10*P[2]/g));

  mc::NLCP<I> CP;
  CP.set_dag( &DAG );
#ifndef USE_DEPS
  CP.set_var( NP, P );
  CP.add_ctr( mc::BASE_OPT::EQ, C[0] );
  CP.add_ctr( mc::BASE_OPT::EQ, C[1] );
#else
  CP.set_var( NP-NC, P );
  CP.set_dep( NC, P+NP-NC, C );
#endif
  CP.setup();

  //CP.options.MIPFILE   = "test3.lp";
  CP.options.DISPLAY     = 2;
  CP.options.MAXITER     = 1000;
  CP.options.CVATOL      = 1e-6;
  CP.options.CVRTOL      = 1e-2;
  CP.options.BRANCHVAR   = mc::SetInv<CVI>::Options::RGREL;
  CP.options.NODEMEAS    = mc::SetInv<CVI>::Options::MAXWIDTH;
  CP.options.STGBCHDEPTH = 0;
  CP.options.STGBCHDRMAX = 0;
  CP.options.STGBCHRTOL  = 1e-2;
  CP.options.DOMREDMAX   = 10;
  CP.options.DOMREDTHRES = 1e-2;
  CP.options.DOMREDBKOFF = 1e-8;
  //CP.options.CTRBACKOFF  = 1e-7;
  CP.options.RELMETH     = mc::NLCP<I>::Options::CHEB;
  CP.options.CMODPROP    = 3;
//  CP.options.CMODSPAR    = true;
  CP.options.CMREDORD    = 5;
  CP.options.CMREDTHRES  = 1e-4;
  CP.options.CMREDWARMS  = false;//true;
  std::cout << CP;

  const I Ip[NP] = { I(0.93,0.99), I(0.,1.), I(0.,1.) };
  //const I Ip[NP] = { I(0.98,0.99), I(0.4,0.45), I(0.15,0.175) };
  //const I Ip[NP] = { I(0.97125,0.973125), I(2.219548e-02,2.648165e-02), I(6.895857e-01,6.908187e-01) };
  //const I Ip[NP] = { I(0.973125,0.975), I(2.219548e-02,2.648165e-02), I(6.895857e-01,6.908187e-01) };
  //const I Ip[NP] = { I(0.97125,0.975), I(2.219548e-02,2.648165e-02), I(6.895857e-01,6.908187e-01) };

  CP.solve( Ip );

#if defined(SAVE_RESULTS )
  std::ofstream ores;
  ores.open( "test3.out", std::ios_base::out );
  CP.output_nodes( ores, 7 );
  ores.close();
  ores.open( "test3b.out", std::ios_base::out );
  CP.output_nodes( ores, 20, 7 );
  ores.close();
#endif

  return 0;
}
