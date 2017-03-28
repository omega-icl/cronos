#include <fstream>
#include <iostream>

#include "interval.hpp"
typedef mc::Interval I;

#define MC__USE_CPLEX
#include "nlgo.hpp"

// PROBLEM 5 IN BEN-TAL ET AL., MATH PROG, 1994
////////////////////////////////////////////////////////////////////////
int main()
////////////////////////////////////////////////////////////////////////
{  
  mc::FFGraph DAG;
  const unsigned NP = 32; mc::FFVar p[NP];
  for( unsigned i=0; i<NP; i++ ) p[i].set( &DAG );

  mc::NLGO<I> NLP;
  //NLP.options.POLIMG.SANDWICH_MAXCUT = 5;
  //NLP.options.POLIMG.SANDWICH_ATOL   = NLP.options.POLIMG.SANDWICH_RTOL  = 1e-5;
  //NLP.options.POLIMG.BREAKPOINT_TYPE = mc::PolImg<I>::Options::NONE;//BIN;//SOS2;//NONE;
  NLP.options.NLPSLV.MAXITER = 100;//NONE;
  NLP.options.NLPSLV.DISPLAY = 0;
  NLP.options.MIPFILE = "";//"test3.lp";
  NLP.options.MIPDISPLAY = 0;
  NLP.options.CMODPROP = 2;
  NLP.options.DISPLAY = 2;
  NLP.options.MAXITER = 0;
  NLP.options.CVATOL = NLP.options.CVRTOL = 1e-3;
  NLP.options.RELMETH = mc::NLGO<I>::Options::CHEB;
  NLP.options.DOMREDMAX = 10;
  NLP.options.DOMREDTHRES = 0.1;
  NLP.options.BRANCHVAR   = mc::SBB<I>::Options::RGREL;
  NLP.options.STGBCHDEPTH = 0;//5;
  NLP.options.STGBCHDRMAX = 2;

  NLP.set_dag( &DAG );  // DAG
  NLP.set_var( NP, p ); // decision variables
  mc::FFVar *q = p, *y = q+12, *z = y+15;  
  NLP.set_obj( mc::NLPSLV_IPOPT::MAX,
                8*z[0]+5*z[1]+9*z[2]+6*z[3]+4*z[4]
                 + (18-6*q[0]-16*q[1]-15*q[2] -12*q[3]) *y[0]
                 + (18-6*q[4]-16*q[5]-15*q[6] -12*q[7]) *y[1]
                 + (18-6*q[8]-16*q[9]-15*q[10]-12*q[11])*y[2]
                 + (15-6*q[0]-16*q[1]-15*q[2] -12*q[3]) *y[3]
                 + (15-6*q[4]-16*q[5]-15*q[6] -12*q[7]) *y[4]
                 + (15-6*q[8]-16*q[9]-15*q[10]-12*q[11])*y[5]
                 + (19-6*q[0]-16*q[1]-15*q[2] -12*q[3]) *y[6]
                 + (19-6*q[4]-16*q[5]-15*q[6] -12*q[7]) *y[7]
                 + (19-6*q[8]-16*q[9]-15*q[10]-12*q[11])*y[8]
                 + (16-6*q[0]-16*q[1]-15*q[2] -12*q[3]) *y[9]
                 + (16-6*q[4]-16*q[5]-15*q[6] -12*q[7]) *y[10]
                 + (16-6*q[8]-16*q[9]-15*q[10]-12*q[11])*y[11]
                 + (14-6*q[0]-16*q[1]-15*q[2] -12*q[3]) *y[12]
                 + (14-6*q[4]-16*q[5]-15*q[6] -12*q[7]) *y[13]
                 + (14-6*q[8]-16*q[9]-15*q[10]-12*q[11])*y[14] );
  NLP.add_ctr( mc::NLPSLV_IPOPT::LE,
                q[2]*y[0] +q[2]*y[3] +q[2]*y[6] +q[2]*y[9]  +q[2]*y[12]
                 +q[6]*y[1] +q[6]*y[4] +q[6]*y[7] +q[6]*y[10] +q[6]*y[13]
                 +q[10]*y[2]+q[10]*y[5]+q[10]*y[8]+q[10]*y[11]+q[10]*y[14]-50 );
  NLP.add_ctr( mc::NLPSLV_IPOPT::LE,
                y[0]+y[1]+y[2]+z[0]-100 );
  NLP.add_ctr( mc::NLPSLV_IPOPT::LE,
                y[3]+y[4]+y[5]+z[1]-200 );
  NLP.add_ctr( mc::NLPSLV_IPOPT::LE,
                y[6]+y[7]+y[8]+z[2]-100 );
  NLP.add_ctr( mc::NLPSLV_IPOPT::LE,
                y[9]+y[10]+y[11]+z[3]-100 );
  NLP.add_ctr( mc::NLPSLV_IPOPT::LE,
                y[12]+y[13]+y[14]+z[4]-100 );
  NLP.add_ctr( mc::NLPSLV_IPOPT::LE, 
                (3*q[0]+q[1]+q[2]+1.5*q[3]-2.5)*y[0]
                 +(3*q[4]+q[5]+q[6]+1.5*q[7]-2.5)*y[1]
                 +(3*q[8]+q[9]+q[10]+1.5*q[11]-2.5)*y[2] - 0.5*z[0] );
  NLP.add_ctr( mc::NLPSLV_IPOPT::LE, 
                (q[0]+3*q[1]+2.5*q[2]+2.5*q[3]-2)*y[0]
                 +(q[4]+3*q[5]+2.5*q[6]+2.5*q[7]-2)*y[1]
                 +(q[8]+3*q[9]+2.5*q[10]+2.5*q[11]-2)*y[2] + 0.5*z[0] );
  NLP.add_ctr( mc::NLPSLV_IPOPT::LE,
                (3*q[0]+q[1]+q[2]+1.5*q[3]-1.5)*y[3]
                 +(3*q[4]+q[5]+q[6]+1.5*q[7]-1.5)*y[4]
                 +(3*q[8]+q[9]+q[10]+1.5*q[11]-1.5)*y[5] + 0.5*z[1] );
  NLP.add_ctr( mc::NLPSLV_IPOPT::LE,
                (q[0]+3*q[1]+2.5*q[2]+2.5*q[3]-2.5)*y[3]
                 +(q[4]+3*q[5]+2.5*q[6]+2.5*q[7]-2.5)*y[4]
                 +(q[8]+3*q[9]+2.5*q[10]+2.5*q[11]-2.5)*y[5] );
  NLP.add_ctr( mc::NLPSLV_IPOPT::LE,
                (3*q[0]+q[1]+q[2]+1.5*q[3]-2)*y[6]
                 +(3*q[4]+q[5]+q[6]+1.5*q[7]-2)*y[7]
                 +(3*q[8]+q[9]+q[10]+1.5*q[11]-2)*y[8] );
  NLP.add_ctr( mc::NLPSLV_IPOPT::LE,
                (q[0]+3*q[1]+2.5*q[2]+2.5*q[3]-2.6)*y[6]
                 +(q[4]+3*q[5]+2.5*q[6]+2.5*q[7]-2.6)*y[7]
                 +(q[8]+3*q[9]+2.5*q[10]+2.5*q[11]-2.6)*y[8] - 0.1*z[2] );
  NLP.add_ctr( mc::NLPSLV_IPOPT::LE,
                (3*q[0]+q[1]+q[2]+1.5*q[3]-2)*y[9]
                 +(3*q[4]+q[5]+q[6]+1.5*q[7]-2)*y[10]
                 +(3*q[8]+q[9]+q[10]+1.5*q[11]-2)*y[11]);
  NLP.add_ctr( mc::NLPSLV_IPOPT::LE,
                (q[0]+3*q[1]+2.5*q[2]+2.5*q[3]-2)*y[9]
                 +(q[4]+3*q[5]+2.5*q[6]+2.5*q[7]-2)*y[10]
                 +(q[8]+3*q[9]+2.5*q[10]+2.5*q[11]-2)*y[11] + 0.5*z[3] );
  NLP.add_ctr( mc::NLPSLV_IPOPT::LE,
                (3*q[0]+q[1]+q[2]+1.5*q[3]-2)*y[12]
                 +(3*q[4]+q[5]+q[6]+1.5*q[7]-2)*y[13]
                 +(3*q[8]+q[9]+q[10]+1.5*q[11]-2)*y[14]);
  NLP.add_ctr( mc::NLPSLV_IPOPT::LE,
                (q[0]+3*q[1]+2.5*q[2]+2.5*q[3]-2)*y[12]
                 +(q[4]+3*q[5]+2.5*q[6]+2.5*q[7]-2)*y[13]
                 +(q[8]+3*q[9]+2.5*q[10]+2.5*q[11]-2)*y[14] + 0.5*z[4] );
  NLP.add_ctr( mc::NLPSLV_IPOPT::EQ,
                q[0]+q[1]+q[2]+q[3]-1 );
  NLP.add_ctr( mc::NLPSLV_IPOPT::EQ,
                q[4]+q[5]+q[6]+q[7]-1 );
  NLP.add_ctr( mc::NLPSLV_IPOPT::EQ,
                q[8]+q[9]+q[10]+q[11]-1 );
  for( unsigned i=0; i<5; i++ ){
    NLP.add_ctr( mc::NLPSLV_IPOPT::EQ,
                  (q[0]+q[1]+q[2]+q[3]-1)*y[3*i] );
    NLP.add_ctr( mc::NLPSLV_IPOPT::EQ,
                  (q[4]+q[5]+q[6]+q[7]-1)*y[3*i+1] );
    NLP.add_ctr( mc::NLPSLV_IPOPT::EQ,
                  (q[8]+q[9]+q[10]+q[11]-1)*y[3*i+2] );
  }

  typedef mc::Interval I;
  I Ip[NP] = { I(0,1), I(0,1), I(0,1), I(0,1), I(0,1), I(0,1),
              I(0,1), I(0,1), I(0,1), I(0,1), I(0,1), I(0,1),
              I(0,100), I(0,100), I(0,100), I(0,200), I(0,200), I(0,200),
              I(0,100), I(0,100), I(0,100), I(0,100), I(0,100), I(0,100),
              I(0,100), I(0,100), I(0,100), I(0,100), I(0,200), I(0,100),
              I(0,100), I(0,100) };
  double p0[NP] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 0 };

  NLP.setup();
  std::cout << NLP;

  //int status = NLP.relax( Ip );
  //if( status == Ipopt::Solve_Succeeded ){
  //  std::cout << "RELAXED NLP SOLUTION: " << std::endl;
  //  std::cout << "  f* = " << NLP.get_objective() << std::endl;
  //  for( unsigned ip=0; ip<NP; ip++ )
  //    std::cout << "  p*(" << ip << ") = " << NLP.get_variable(p[ip])
  //              << std::endl;
  //}

  int status = NLP.solve( Ip, 0, p0 );
  NLP.stats.display();

  return status;
}
