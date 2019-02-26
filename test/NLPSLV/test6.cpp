#include <fstream>
#include <iomanip>
#include "nlpslv_ipopt.hpp"
#include "interval.hpp"

const unsigned NPTS = 4;
const double p1[NPTS] = { 2, 4, 4, 6 };
const double p2[NPTS] = { 3, 4, 2, 2 };

////////////////////////////////////////////////////////////////////////
int main()
////////////////////////////////////////////////////////////////////////
{  
  mc::FFGraph DAG;
  const unsigned NP = 5; mc::FFVar P[NP];
  for( unsigned i=0; i<NP; i++ ) P[i].set( &DAG );
  mc::FFVar &alpha = P[0], &y1L = P[1], &y1U = P[2], &y2L = P[3], &y2U = P[4];

  Ipopt::SmartPtr<mc::NLPSLV_IPOPT> NLP = new mc::NLPSLV_IPOPT;
  NLP->set_dag( &DAG );                       // DAG
  NLP->set_var( NP, P );                      // decision variables
  NLP->set_obj( mc::NLPSLV_IPOPT::MIN, (y1U-y1L)*(y2U-y2L) );   // objective
  for( unsigned k=0; k<NPTS; k++ ){
    NLP->add_ctr( mc::NLPSLV_IPOPT::LE, p1[k]*cos(alpha)+p2[k]*sin(alpha)-y1U );
    NLP->add_ctr( mc::NLPSLV_IPOPT::GE, p1[k]*cos(alpha)+p2[k]*sin(alpha)-y1L );
    NLP->add_ctr( mc::NLPSLV_IPOPT::LE, p2[k]*cos(alpha)-p1[k]*sin(alpha)-y2U );
    NLP->add_ctr( mc::NLPSLV_IPOPT::GE, p2[k]*cos(alpha)-p1[k]*sin(alpha)-y2L );
  }
  typedef mc::Interval I;
  I Ip[NP] = { I(0.,mc::PI/2), I(-10.,10.), I(-10.,10.), I(-10.,10.), I(-10.,10.) };

  NLP->options.DISPLAY = 0;
  NLP->options.MAXITER = 100;
  NLP->options.CVTOL = 1e-9;
  NLP->options.GRADIENT = mc::NLPSLV_IPOPT::Options::BACKWARD;
  NLP->options.HESSIAN  = mc::NLPSLV_IPOPT::Options::LBFGS;
  NLP->setup();

  int status = NLP->solve( 10, Ip );

  if( status == Ipopt::Solve_Succeeded ){
    std::cout << "NLP (MULTISTART) SOLUTION: " << std::endl;
    std::cout << "  f* = " << NLP->solution().f << std::endl;
    for( unsigned int ip=0; ip<NP; ip++ )
      std::cout << "  p*(" << ip << ") = " << NLP->solution().p[ip]
                << std::endl;
    std::cout << "FEASIBLE:   " << NLP->is_feasible( 1e-7 ) << std::endl;
    std::cout << "STATIONARY: " << NLP->is_stationary( 1e-7 ) << std::endl;
  }

  return 0;
}
