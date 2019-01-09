#include <fstream>
#include <iomanip>
#include "nlpslv_ipopt.hpp"
#include "interval.hpp"
typedef mc::Interval I;

const double k1    = 5.0;  // [/min]
const double k2    = 2.0;  // [/min]
const double V     = 1.0;  // [L]
const double CAU   = 0.2;  // [mol/L]
const double CAinU = 1.0;  // [mol/L]
const double FU    = 20.0; // [L/min]

// EXAMPLE 7 IN RYOO & SAHINIDIS, C&CE, 1995
////////////////////////////////////////////////////////////////////////
int main()
////////////////////////////////////////////////////////////////////////
{  
  mc::FFGraph DAG;
  const unsigned NP = 2; mc::FFVar p[NP];
  for( unsigned i=0; i<NP; i++ ) p[i].set( &DAG );
  const unsigned NX = 2; mc::FFVar x[NX], f[NX];
  for( unsigned i=0; i<NX; i++ ) x[i].set( &DAG );
  f[0] = p[0]*(p[1]-x[0])-k1*x[0]*V;
  f[1] = -p[0]*x[1]+(k1*x[0]-k2*x[1])*V;

  Ipopt::SmartPtr<mc::NLPSLV_IPOPT> NLP = new mc::NLPSLV_IPOPT;
  NLP->set_dag( &DAG );  // DAG
  NLP->set_var( NP, p ); // decision variables
  NLP->set_dep( NX, x, f ); // decision variables

  NLP->set_obj( mc::NLPSLV_IPOPT::MAX, x[1] );
  NLP->add_ctr( mc::NLPSLV_IPOPT::LE, x[0]-CAU );

  I Ip[NP+NX] = { I(0,FU), I(0,CAinU), I(0,1e1), I(0,1e1) };
  double p0[NP+NX] = { 0., 0., 0., 0. };

  NLP->options.DISPLAY = 5;
  NLP->options.MAXITER = 100;
  NLP->options.CVTOL = 1e-9;
  NLP->options.GRADIENT = mc::NLPSLV_IPOPT::Options::BACKWARD;
  NLP->options.HESSIAN  = mc::NLPSLV_IPOPT::Options::LBFGS;
  NLP->setup();
  //int status = NLP->solve( Ip, p0 );
  int status = NLP->solve( Ip );

  if( status == Ipopt::Solve_Succeeded ){
    std::cout << "NLP (LOCAL) SOLUTION: " << std::endl;
    std::cout << "  f* = " << NLP->solution().f << std::endl;
    for( unsigned ip=0; ip<NP; ip++ )
      std::cout << "  p*(" << ip << ") = " << NLP->solution().p[ip]
                << std::endl;
    for( unsigned ix=0; ix<NX; ix++ )
      std::cout << "  x*(" << ix << ") = " << NLP->solution().p[NP+ix]
                << std::endl;
  }

  return 0;
}
