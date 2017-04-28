#include <fstream>
#include <iostream>

#include "interval.hpp"
typedef mc::Interval I;

#undef MC__NLGO_SHOW_DEPS
#undef MC__NLGO_SHOW_BOXES
#define MC__USE_CPLEX
#include "nlgo.hpp"

////////////////////////////////////////////////////////////////////////
int main()
////////////////////////////////////////////////////////////////////////
{  
  mc::FFGraph DAG;
  const unsigned NP = 2; mc::FFVar P[NP];
  for( unsigned i=0; i<NP; i++ ) P[i].set( &DAG );

  mc::NLGO<I> NLP;
  NLP.set_dag( &DAG );                       // DAG
  NLP.set_var( NP, P );                      // decision variables
  //NLP.set_obj( mc::BASE_NLP::MAX, P[0]+P[1] );   // objective
  NLP.set_obj( mc::BASE_NLP::MAX, sqr(P[0]+P[1]-1.) );   // objective
  NLP.add_ctr( mc::BASE_NLP::EQ, P[0]*P[1]-4. ); // constraints

  NLP.options.POLIMG.SANDWICH_MAXCUT = 4;
  //NLP.options.POLIMG.BREAKPOINT_TYPE = mc::PolImg<I>::Options::SOS2;//NONE;
  //NLP.options.POLIMG.DCDECOMP_SCALE = false;
  NLP.options.MIPFILE = "test1.lp";
  NLP.options.MAXITER = 10;
  NLP.options.BLKDECUSE = true; //false;
  NLP.options.CMODDEPS = 0;
  NLP.options.AEBND.DISPLAY = 0;
  NLP.options.PREPROC = true;
  NLP.options.DISPLAY = 2;
  NLP.options.MIPDISPLAY = 0;
  NLP.options.NLPSLV.DISPLAY = 0;
  NLP.setup();

  I Ip[NP] = { I(0.,6.), I(1.,4.) };
  int status = NLP.solve( Ip );
  NLP.stats.display();

  return status;
}
