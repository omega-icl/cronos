#include <fstream>
#include <iomanip>
#include "nlpslv_ipopt.hpp"
#include "interval.hpp"

////////////////////////////////////////////////////////////////////////
int main()
////////////////////////////////////////////////////////////////////////
{  
  mc::FFGraph DAG;
  const unsigned NP = 2; mc::FFVar P[NP];
  for( unsigned i=0; i<NP; i++ ) P[i].set( &DAG );

  Ipopt::SmartPtr<mc::NLPSLV_IPOPT> NLP = new mc::NLPSLV_IPOPT;
  NLP->set_dag( &DAG );                       // DAG
  NLP->set_var( NP, P );                      // decision variables
  NLP->set_obj( mc::NLPSLV_IPOPT::MAX, P[0]+P[1] );   // objective
  NLP->add_ctr( mc::NLPSLV_IPOPT::LE, P[0]*P[1]-4. ); // constraints

  typedef mc::Interval I;
  I Ip[NP] = { I(0.,6.), I(0.,4.) };
  double p0[NP] = { 5., 1. };

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
