#include <fstream>
#include <iomanip>
#include "nlpslv_ipopt.hpp"
#include "interval.hpp"

// EXAMPLE 7 IN RYOO & SAHINIDIS, C&CE, 1995
////////////////////////////////////////////////////////////////////////
int main()
////////////////////////////////////////////////////////////////////////
{  
  mc::FFGraph DAG;
  const unsigned NP = 10; mc::FFVar p[NP];
  for( unsigned i=0; i<NP; i++ ) p[i].set( &DAG );

  Ipopt::SmartPtr<mc::NLPSLV_IPOPT> NLP = new mc::NLPSLV_IPOPT;
  NLP->set_dag( &DAG );  // DAG
  NLP->set_var( NP, p ); // decision variables
  NLP->set_obj( mc::NLPSLV_IPOPT::MIN, -9*p[4]-15*p[8]+6*p[0]+16*p[1]+10*p[5] );
  NLP->add_ctr( mc::NLPSLV_IPOPT::EQ, p[0]+p[1]-p[2]-p[3] );
  NLP->add_ctr( mc::NLPSLV_IPOPT::EQ, p[2]+p[6]-p[4] );
  NLP->add_ctr( mc::NLPSLV_IPOPT::EQ, p[3]+p[7]-p[8] );
  NLP->add_ctr( mc::NLPSLV_IPOPT::EQ, p[6]+p[7]-p[5] );
  NLP->add_ctr( mc::NLPSLV_IPOPT::LE, p[9]*p[2]+2*p[6]-2.5*p[4] );
  NLP->add_ctr( mc::NLPSLV_IPOPT::LE, p[9]*p[3]+2*p[7]-1.5*p[8] );
  NLP->add_ctr( mc::NLPSLV_IPOPT::LE, 3*p[0]+p[1]-p[9]*p[2]-p[9]*p[3] );

  typedef mc::Interval I;
  I Ip[NP] = { I(0,300), I(0,300), I(0,100), I(0,200), I(0,100), I(0,300),
               I(0,100), I(0,200), I(0,200), I(1,3) };
  double p0[NP] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 };

  NLP->options.DISPLAY = 5;
  NLP->options.MAXITER = 100;
  NLP->options.CVTOL = 1e-9;
  NLP->options.GRADIENT = mc::NLPSLV_IPOPT::Options::BACKWARD;
  NLP->options.HESSIAN  = mc::NLPSLV_IPOPT::Options::LBFGS;
  NLP->setup();
  int status = NLP->solve( Ip, p0 );

  if( status == Ipopt::Solve_Succeeded ){
    std::cout << "NLP (LOCAL) SOLUTION: " << std::endl;
    std::cout << "  f* = " << NLP->solution().f << std::endl;
    for( unsigned int ip=0; ip<NP; ip++ )
      std::cout << "  p*(" << ip << ") = " << NLP->solution().p[ip]
                << std::endl;
  }

  NLP->options.HESSIAN  = mc::NLPSLV_IPOPT::Options::EXACT;
  NLP->setup();
  status = NLP->solve( Ip, p0 );

  if( status == Ipopt::Solve_Succeeded ){
    std::cout << "NLP (LOCAL) SOLUTION: " << std::endl;
    std::cout << "  f* = " << NLP->solution().f << std::endl;
    for( unsigned int ip=0; ip<NP; ip++ )
      std::cout << "  p*(" << ip << ") = " << NLP->solution().p[ip]
                << std::endl;
    std::cout << "FEASIBLE:   " << NLP->is_feasible( 1e-5 ) << std::endl;
    std::cout << "STATIONARY: " << NLP->is_stationary( 1e-5 ) << std::endl;
  }

  return 0;
}
