#include <fstream>
#include <iomanip>
#include "nlpslv_ipopt.hpp"
#include "interval.hpp"

// EXAMPLE 2 IN BEN-TAL ET AL., MATH PROG, 1994
////////////////////////////////////////////////////////////////////////
int main()
////////////////////////////////////////////////////////////////////////
{  
  mc::FFGraph DAG;
  const unsigned NP = 32; mc::FFVar p[NP];
  for( unsigned i=0; i<NP; i++ ) p[i].set( &DAG );

  Ipopt::SmartPtr<mc::NLPSLV_IPOPT> NLP = new mc::NLPSLV_IPOPT;
  NLP->set_dag( &DAG );  // DAG
  NLP->set_var( NP, p ); // decision variables
  mc::FFVar *q = p, *y = q+12, *z = y+15;  
  NLP->set_obj( mc::NLPSLV_IPOPT::MAX,
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
  NLP->add_ctr( mc::NLPSLV_IPOPT::LE,
                q[2]*y[0] +q[2]*y[3] +q[2]*y[6] +q[2]*y[9]  +q[2]*y[12]
                 +q[6]*y[1] +q[6]*y[4] +q[6]*y[7] +q[6]*y[10] +q[6]*y[13]
                 +q[10]*y[2]+q[10]*y[5]+q[10]*y[8]+q[10]*y[11]+q[10]*y[14]-50 );
  NLP->add_ctr( mc::NLPSLV_IPOPT::LE,
                y[0]+y[1]+y[2]+z[0]-100 );
  NLP->add_ctr( mc::NLPSLV_IPOPT::LE,
                y[3]+y[4]+y[5]+z[1]-200 );
  NLP->add_ctr( mc::NLPSLV_IPOPT::LE,
                y[6]+y[7]+y[8]+z[2]-100 );
  NLP->add_ctr( mc::NLPSLV_IPOPT::LE,
                y[9]+y[10]+y[11]+z[3]-100 );
  NLP->add_ctr( mc::NLPSLV_IPOPT::LE,
                y[12]+y[13]+y[14]+z[4]-100 );
  NLP->add_ctr( mc::NLPSLV_IPOPT::LE, 
                (3*q[0]+q[1]+q[2]+1.5*q[3]-2.5)*y[0]
                 +(3*q[4]+q[5]+q[6]+1.5*q[7]-2.5)*y[1]
                 +(3*q[8]+q[9]+q[10]+1.5*q[11]-2.5)*y[2] - 0.5*z[0] );
  NLP->add_ctr( mc::NLPSLV_IPOPT::LE, 
                (q[0]+3*q[1]+2.5*q[2]+2.5*q[3]-2)*y[0]
                 +(q[4]+3*q[5]+2.5*q[6]+2.5*q[7]-2)*y[1]
                 +(q[8]+3*q[9]+2.5*q[10]+2.5*q[11]-2)*y[2] + 0.5*z[0] );
  NLP->add_ctr( mc::NLPSLV_IPOPT::LE,
                (3*q[0]+q[1]+q[2]+1.5*q[3]-1.5)*y[3]
                 +(3*q[4]+q[5]+q[6]+1.5*q[7]-1.5)*y[4]
                 +(3*q[8]+q[9]+q[10]+1.5*q[11]-1.5)*y[5] + 0.5*z[1] );
  NLP->add_ctr( mc::NLPSLV_IPOPT::LE,
                (q[0]+3*q[1]+2.5*q[2]+2.5*q[3]-2.5)*y[3]
                 +(q[4]+3*q[5]+2.5*q[6]+2.5*q[7]-2.5)*y[4]
                 +(q[8]+3*q[9]+2.5*q[10]+2.5*q[11]-2.5)*y[5] );
  NLP->add_ctr( mc::NLPSLV_IPOPT::LE,
                (3*q[0]+q[1]+q[2]+1.5*q[3]-2)*y[6]
                 +(3*q[4]+q[5]+q[6]+1.5*q[7]-2)*y[7]
                 +(3*q[8]+q[9]+q[10]+1.5*q[11]-2)*y[8] );
  NLP->add_ctr( mc::NLPSLV_IPOPT::LE,
                (q[0]+3*q[1]+2.5*q[2]+2.5*q[3]-2.6)*y[6]
                 +(q[4]+3*q[5]+2.5*q[6]+2.5*q[7]-2.6)*y[7]
                 +(q[8]+3*q[9]+2.5*q[10]+2.5*q[11]-2.6)*y[8] - 0.1*z[2] );
  NLP->add_ctr( mc::NLPSLV_IPOPT::LE,
                (3*q[0]+q[1]+q[2]+1.5*q[3]-2)*y[9]
                 +(3*q[4]+q[5]+q[6]+1.5*q[7]-2)*y[10]
                 +(3*q[8]+q[9]+q[10]+1.5*q[11]-2)*y[11]);
  NLP->add_ctr( mc::NLPSLV_IPOPT::LE,
                (q[0]+3*q[1]+2.5*q[2]+2.5*q[3]-2)*y[9]
                 +(q[4]+3*q[5]+2.5*q[6]+2.5*q[7]-2)*y[10]
                 +(q[8]+3*q[9]+2.5*q[10]+2.5*q[11]-2)*y[11] + 0.5*z[3] );
  NLP->add_ctr( mc::NLPSLV_IPOPT::LE,
                (3*q[0]+q[1]+q[2]+1.5*q[3]-2)*y[12]
                 +(3*q[4]+q[5]+q[6]+1.5*q[7]-2)*y[13]
                 +(3*q[8]+q[9]+q[10]+1.5*q[11]-2)*y[14]);
  NLP->add_ctr( mc::NLPSLV_IPOPT::LE,
                (q[0]+3*q[1]+2.5*q[2]+2.5*q[3]-2)*y[12]
                 +(q[4]+3*q[5]+2.5*q[6]+2.5*q[7]-2)*y[13]
                 +(q[8]+3*q[9]+2.5*q[10]+2.5*q[11]-2)*y[14] + 0.5*z[4] );
  NLP->add_ctr( mc::NLPSLV_IPOPT::EQ,
                q[0]+q[1]+q[2]+q[3]-1 );
  NLP->add_ctr( mc::NLPSLV_IPOPT::EQ,
                q[4]+q[5]+q[6]+q[7]-1 );
  NLP->add_ctr( mc::NLPSLV_IPOPT::EQ,
                q[8]+q[9]+q[10]+q[11]-1 );

  typedef mc::Interval I;
  I Ip[NP] = { I(0,1), I(0,1), I(0,1), I(0,1), I(0,1), I(0,1),
              I(0,1), I(0,1), I(0,1), I(0,1), I(0,1), I(0,1),
              I(0,100), I(0,100), I(0,100), I(0,200), I(0,200), I(0,200),
              I(0,100), I(0,100), I(0,100), I(0,100), I(0,100), I(0,100),
              I(0,100), I(0,100), I(0,100), I(0,100), I(0,200), I(0,100),
              I(0,100), I(0,100) };
//  double p0[NP] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
//                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
//                    0, 0, 0, 0, 0 };
  double p0[NP] = { 1, 0, 0, 0, 0, 0, 0, 1, 0.2757, 0, 0, 0.7243,
                    51.5527, 0, 17.9504, 0, 200, 0, 0, 0, 0,
                    7.9609, 0, 92.0391, 15.7866, 20.5622, 63.6512,
                    30.4969, 0, 100, 0, 0 };

  NLP->options.DISPLAY = 5;
  NLP->options.MAXITER = 100;
  NLP->options.CVTOL = 1e-9;
  NLP->options.GRADIENT = mc::NLPSLV_IPOPT::Options::FORWARD;//BACKWARD;
  NLP->options.HESSIAN  = mc::NLPSLV_IPOPT::Options::LBFGS;
  NLP->setup();
  Ipopt::ApplicationReturnStatus status = NLP->solve( Ip, p0 );

  if( status == Ipopt::Solve_Succeeded ){
    std::cout << "NLP (LOCAL) SOLUTION: " << std::endl;
    std::cout << "  f* = " << NLP->solution().f << std::endl;
    for( unsigned int ip=0; ip<NP; ip++ )
      std::cout << "  p*(" << ip << ") = " << NLP->solution().p[ip]
                << std::endl;
  }

  NLP->options.HESSIAN  = mc::NLPSLV_IPOPT::Options::EXACT;
  NLP->options.TESTDER  = true;
  NLP->setup();
  status = NLP->solve( Ip, p0 );

  if( status == Ipopt::Solve_Succeeded ){
    std::cout << "NLP (LOCAL) SOLUTION: " << std::endl;
    std::cout << "  f* = " << NLP->solution().f << std::endl;
    for( unsigned int ip=0; ip<NP; ip++ )
      std::cout << "  p*(" << ip << ") = " << NLP->solution().p[ip]
                << std::endl;
  }

  return 0;
}
