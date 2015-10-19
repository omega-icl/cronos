#include <fstream>
#include <iomanip>
#include "nlgo_gurobi.hpp"
#include "interval.hpp"
typedef mc::Interval I;

////////////////////////////////////////////////////////////////////////
int main()
////////////////////////////////////////////////////////////////////////
{  
  mc::FFGraph DAG;
  const unsigned NP = 2; mc::FFVar P[NP];
  for( unsigned i=0; i<NP; i++ ) P[i].set( &DAG );

  mc::NLGO_GUROBI<I> NLP;
  NLP.options.POLIMG.SANDWICH_MAXCUT = 4;
  NLP.options.POLIMG.BREAKPOINT_TYPE = mc::PolImg<I>::Options::SOS2;//NONE;
  NLP.options.POLIMG.DCDECOMP_SCALE  = false;
  NLP.options.MIPFILE = "test1.lp";
  NLP.options.MAXITER = 10;

  NLP.set_dag( &DAG );                       // DAG
  NLP.set_var( NP, P );                      // decision variables
  NLP.set_obj( mc::BASE_NLP::MAX, P[0]+P[1] );   // objective
  //NLP.set_obj( mc::BASE_NLP::MAX, P[0]+P[1]+P[0]/P[1] );   // objective
  NLP.add_ctr( mc::BASE_NLP::LE, P[0]*P[1]-4. ); // constraints

  I Ip[NP] = { I(0.,6.), I(1.,4.) };

  NLP.setup();
  NLP.options.PREPROC = true;
  NLP.options.MIPDISPLAY = 0;
  NLP.options.NLPSLV.DISPLAY = 0;
  //NLP.tighten( Ip );
  //int status = NLP.relax( Ip );
  ////if( status == Ipopt::Solve_Succeeded ){
  //std::cout << "RELAXED NLP SOLUTION: " << std::endl;
  //std::cout << "  f* = " << NLP.get_objective() << std::endl;
  //  for( unsigned ip=0; ip<NP; ip++ )
  //    std::cout << "  p*(" << ip << ") = " << NLP.get_variable(P[ip])
  //              << std::endl;
  ////}
  int status = NLP.solve( Ip );

  return 0;
}
