#include <fstream>
#include <iomanip>

#include "interval.hpp"
typedef mc::Interval I;

#define MC__USE_CPLEX
#include "nlgo.hpp"

// EXAMPLE 7 IN RYOO & SAHINIDIS, C&CE, 1995
////////////////////////////////////////////////////////////////////////
int main()
////////////////////////////////////////////////////////////////////////
{  
  mc::FFGraph DAG;
  const unsigned NP = 10; mc::FFVar p[NP];
  for( unsigned i=0; i<NP; i++ ) p[i].set( &DAG );

  mc::NLGO<I> NLP;
  //NLP.options.POLIMG.SANDWICH_MAXCUT = 10;
  //NLP.options.POLIMG.SANDWICH_ATOL   = NLP.options.POLIMG.SANDWICH_RTOL  = 1e-5;
  //NLP.options.POLIMG.BREAKPOINT_TYPE = mc::PolImg<I>::Options::SOS2;//NONE;
  //NLP.options.POLIMG.DCDECOMP_SCALE  = true;
  //NLP.options.MIPFILE = "test2.lp";
  NLP.options.MIPDISPLAY = 0;
  NLP.options.NLPSLV.DISPLAY = 5;
  NLP.options.DISPLAY   = 2;
  NLP.options.MAXITER   = 10;
  NLP.options.CVATOL    = NLP.options.CVRTOL    = 1e-4;
  NLP.options.MIPABSGAP = NLP.options.MIPRELGAP = 0.;//1e-7;
  NLP.options.DOMREDMAX = 10;
  NLP.options.RELMETH   = mc::NLGO<I>::Options::CHEB;//DRL;//HYBRID;
  NLP.options.BLKDECUSE = true; //false;
  NLP.options.CMODDEPS  = 0;

  NLP.set_dag( &DAG );  // DAG
  NLP.set_var( NP, p ); // decision variables
  NLP.set_obj( mc::BASE_NLP::MIN, -9*p[4]-15*p[8]+6*p[0]+16*p[1]+10*p[5] );
  NLP.add_ctr( mc::BASE_NLP::EQ, p[0]+p[1]-p[2]-p[3] );
  NLP.add_ctr( mc::BASE_NLP::EQ, p[2]+p[6]-p[4] );
  NLP.add_ctr( mc::BASE_NLP::EQ, p[3]+p[7]-p[8] );
  NLP.add_ctr( mc::BASE_NLP::EQ, p[6]+p[7]-p[5] );
  NLP.add_ctr( mc::BASE_NLP::LE, p[9]*p[2]+2*p[6]-2.5*p[4] );
  NLP.add_ctr( mc::BASE_NLP::LE, p[9]*p[3]+2*p[7]-1.5*p[8] );
  NLP.add_ctr( mc::BASE_NLP::LE, 3*p[0]+p[1]-p[9]*p[2]-p[9]*p[3] );
  //NLP.add_ctr( mc::BASE_NLP::LE, 3*p[0]+p[1]-p[9]*(p[2]+p[3]) );

  I Ip[NP] = { I(0,300), I(0,300), I(0,100), I(0,200), I(0,100), I(0,300),
               I(0,100), I(0,200), I(0,200), I(1,3) };
  //I Ip[NP] = { I(0,66.67), I(75,133.34), I(0,72.73), I(75,166.67), I(0,100), I(66.66,175),
  //             I(0,75), I(33.33,100), I(150,200), I(1,2.25) };
  //double p0[NP] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 };

  NLP.setup();
  //int status = NLP.relax( Ip );
  int status = NLP.solve( Ip );
  NLP.stats.display();

  //if( status == Ipopt::Solve_Succeeded ){
  //  std::cout << "RELAXED NLP SOLUTION: " << std::endl;
  //  std::cout << "  f* = " << NLP.get_objective() << std::endl;
  //  for( unsigned ip=0; ip<NP; ip++ )
  //    std::cout << "  p*(" << ip << ") = " << NLP.get_variable(p[ip])
  //              << std::endl;
  //}

  return status;
}
