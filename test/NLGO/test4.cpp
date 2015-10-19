#define USE_PROFIL

#include <fstream>
#include <iomanip>
#include "nlgo_gurobi.hpp"
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

// PROBLEM 71 IN THE HOCK-SCHITTKOWSKY TEST SUITE
// (Lecture Notes in Economics and Mathematical Systems, 187, 1981
//  doi: 10.1007/978-3-642-48320-2)
////////////////////////////////////////////////////////////////////////
int main()
////////////////////////////////////////////////////////////////////////
{
  mc::FFGraph DAG;
  const unsigned NP = 4; mc::FFVar p[NP];
  for( unsigned i=0; i<NP; i++ ) p[i].set( &DAG );

  mc::NLGO_GUROBI<I> NLP;
  //NLP.options.CVATOL = NLP.options.CVRTOL = 1e-5;
  //NLP.options.MIPABSGAP = NLP.options.MIPRELGAP = 1e-9;
  //NLP.options.MAXITER = 0;
  NLP.options.DOMREDMAX = 10;
  NLP.options.DOMREDTHRES = 0.1;
  NLP.options.MAXITER = 0;
  NLP.options.MIPDISPLAY = 0;

/*
  NLP.options.POLIMG.SANDWICH_MAXCUT = 10;
  NLP.options.POLIMG.SANDWICH_ATOL   = NLP.options.POLIMG.SANDWICH_RTOL  = 1e-5;
  NLP.options.POLIMG.BREAKPOINT_TYPE = mc::PolImg<I>::Options::BIN;//SOS2;//NONE;
  NLP.options.POLIMG.DCDECOMP_SCALE  = false;//true;
  //NLP.options.MIPFILE = "test4.lp";
  NLP.options.MIPDISPLAY = 0;
  NLP.options.DISPLAY = 2;
  NLP.options.NLPSLV.DISPLAY = 0;
  NLP.options.MAXITER = 21;
  NLP.options.PREPROC = true;
  NLP.options.DOMREDMAX = 10;
  NLP.options.PRESOS2BIGM = -1;
  NLP.options.CVATOL    = NLP.options.CVRTOL    = 1e-5;
*/
  NLP.set_dag( &DAG );  // DAG
  NLP.set_var( NP, p ); // decision variables
  NLP.set_obj( mc::BASE_NLP::MIN,
               (p[0]*p[3])*(p[0]+p[1]+p[2])+p[2] );
  NLP.add_ctr( mc::BASE_NLP::GE,
               (p[0]*p[3])*p[1]*p[2]-25 );
  NLP.add_ctr( mc::BASE_NLP::EQ,
               sqr(p[0])+sqr(p[1])+sqr(p[2])+sqr(p[3])-40 );
  NLP.setup();
  //std::cout << NLP;

  I Ip[NP] = { I(1,5), I(1,5), I(1,5), I(1,5) };
  unsigned Tp[NP] = { 0, 1, 1, 0 };
  // double p0[NP] = { 1, 5, 5, 1 };

  // Global optimization using SBB
  NLP.options.CSALGO  = mc::NLGO_GUROBI<I>::Options::SBB;
  NLP.options.RELMETH = mc::NLGO_GUROBI<I>::Options::DRL;
  NLP.solve( Ip );
  NLP.solve( Ip, Tp );

  NLP.options.RELMETH = mc::NLGO_GUROBI<I>::Options::CHEB;
  NLP.options.CMODPROP = 3;
  NLP.solve( Ip );
  NLP.solve( Ip, Tp );

  NLP.options.RELMETH = mc::NLGO_GUROBI<I>::Options::HYBRID;
  NLP.options.CMODPROP = 3;
  NLP.solve( Ip );
  NLP.solve( Ip, Tp );

  // Global optimization using PWL
  NLP.options.CSALGO = mc::NLGO_GUROBI<I>::Options::PWL;
  NLP.options.RELMETH = mc::NLGO_GUROBI<I>::Options::CHEB;
  NLP.options.POLIMG.BREAKPOINT_TYPE = mc::PolImg<I>::Options::BIN;//SOS2;//NONE;
  //NLP.solve( Ip );
  NLP.solve( Ip, Tp );

  return 0;
}
