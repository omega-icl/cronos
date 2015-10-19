#include <fstream>
#include <iomanip>
#include "nlpslv_ipopt.hpp"
#include "interval.hpp"

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

  Ipopt::SmartPtr<mc::NLPSLV_IPOPT> NLP = new mc::NLPSLV_IPOPT;
  NLP->set_dag( &DAG );  // DAG
  NLP->set_var( NP, p ); // decision variables
  NLP->set_obj( mc::NLPSLV_IPOPT::MIN,
                p[0]*p[3]*(p[0]+p[1]+p[2])+p[2] );
  NLP->add_ctr( mc::NLPSLV_IPOPT::GE,
                p[0]*p[1]*p[2]*p[3]-25 );
  NLP->add_ctr( mc::NLPSLV_IPOPT::EQ,
                sqr(p[0])+sqr(p[1])+sqr(p[2])+sqr(p[3])-40 );

  typedef mc::Interval I;
  I Ip[NP] = { I(1,5), I(1,5), I(1,5), I(1,5) };
  double p0[NP] = { 1, 5, 5, 1 };

  NLP->options.DISPLAY = 5;
  NLP->options.MAXITER = 100;
  NLP->options.CVTOL = 1e-9;
  NLP->options.GRADIENT = mc::NLPSLV_IPOPT::Options::BACKWARD;
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
