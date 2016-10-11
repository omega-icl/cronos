#define SAVE_RESULTS
////////////////////////////////////////////////////////////////////////////////
// Find all solutions of the system of nonlinear inequalities:
//   yk - 0.5 <= p1*exp(-p2*tk) + p3*exp(-p4*tk) <= yk + 0.5
//   for 10 measurement pairs (t1,y1),...,(t10,y10) as follows:
//      tk := [ 0.75, 1.5,  2.25, 3.,     6.,    9.,   13.,   17.,   21.,   25. ]
//      yk := [ 7.39, 4.09, 1.74, 0.097, -2.57, -2.71, -2.07, -1.44, -0.98, -0.66 ]
//   and (p1,...,p4) in [0,1]^4
////////////////////////////////////////////////////////////
#include <iostream>
#include <fstream>
#include <iomanip>
#include "nlgo.hpp"
#include "nlcp.hpp"
#include "interval.hpp"
typedef mc::Interval I;
typedef mc::CVar<I> CVI;
//data
const unsigned NK = 10;
const double Tm[NK] = { 0.75, 1.5,  2.25, 3.,     6.,    9.,   13.,   17.,   21.,   25. };
const double Ym[NK] = { 7.39, 4.09, 1.74, 0.097, -2.57, -2.71, -2.07, -1.44, -0.98, -0.66 };//0.66
const double dYm = 0.5;


//predictive model g(Input)=p1*exp(-p2*input) + p3*exp(-p4*input)
template <class T>
T Ymod
( double Input , const T*P)
{
//return  P[0]*exp(-P[1]*Input) + P[2]*exp(-P[3]*Input);
return  (58.*P[0]+2.)*exp(-P[1]*Input) + (29.*P[2]-30.)*exp(-0.5*P[3]*Input);

}

//number of parameters in model;
const unsigned NP=4;

////////////////////////////////////////////////////////////////////////////////
// min max problem to find upper bound U
//   OBJ U=max sum(E(k)) for k=1 - size(data) 
//   CONSTRAINTS for least squares error for 10 measurement pairs as follows
//      E(k) >= (y(k)-e(k)-g(k))^2 for k=1 - size(data)
//   with e(k) within error bounds [-dYm,+dYm] and (p1,...,p4) in [0,1]^4

////////////////////////////////////////////////////////////
int main(){
  //DAG construction
  mc::FFGraph DAG;
  const unsigned  NT=NP+NK; mc::FFVar k[NT];
  for( unsigned i=0; i<NT; i++ ) k[i].set( &DAG );

  mc::NLGO<I> NLP; //defining nlp optimizer class 
  NLP.set_dag( &DAG );  // DAG
  NLP.set_var( NT, k ); // decision variables

  mc::FFVar SUM = 0.;
  for( unsigned i=0; i<NK; i++ ){
    NLP.add_ctr( mc::BASE_NLP::GE, k[NP+i] - mc::sqr( Ym[i] + dYm - Ymod(Tm[i],k) ) ); //constraints
    NLP.add_ctr( mc::BASE_NLP::GE, k[NP+i] - mc::sqr( Ym[i] - dYm - Ymod(Tm[i],k) ) ); //constraints
    SUM += k[NP+i];
  }
  NLP.set_obj( mc::NLGO<I>::MIN, SUM ); //objective
  NLP.setup();

  I Ip[NT];
  Ip[0] = I(0.,1); Ip[1] = I(0.,1.); Ip[2] = I(0.,1.); Ip[3] = I(0.,1.);
  for( unsigned i=0; i<NK; i++ ) Ip[NP+i] = dYm * I(-1.,1.);

  double k0[NT]; //initial guesses
  k0[0] = 0.5; k0[1] = 0.5; k0[2] = 0.5; k0[3] = 0.5;
  for( unsigned i=0; i<NK; i++ ) k0[NP+i] = 0e0;

  // Local optimization
  NLP.options.NLPSLV.DISPLAY = 5;
  NLP.options.NLPSLV.MAXITER = 100;
  NLP.local( Ip, k0 );
  std::cout << std::setprecision(4);
  std::cout << "  f = " << NLP.get_local_solution().f;
  for( unsigned i=0; i<NT; i++ )
    std::cout << "  k(" << i << ") = " << NLP.get_local_solution().p[i];
  std::cout << std::endl;
  { int dum; std::cin >> dum; }
/*
  // Global optimization using SBB
  NLP.options.DISPLAY = 2;
  //NLP.options.MIPFILE = "test4b.lp";
  NLP.options.NLPSLV.DISPLAY = 0;
  NLP.options.CSALGO  = mc::NLGO<I>::Options::SBB;
  NLP.options.RELMETH = mc::NLGO<I>::Options::DRL;//HYBRID;//CHEB;
  NLP.options.CMODPROP = 1;
  NLP.solve( Ip, 0, NLP.get_local_solution().p );
  //NLP.solve( Ip, 0, k0 );
*/

  //Constraint projection problem to find all partitions of the parameter set Pe such that for any p in Pe 
  //sum((y(k)-e(k)-g(k))^2)<U for  k = 1 - size(data)

  mc::NLCP<I> CP; // defining constraint projection class
  CP.set_dag( &DAG ); //DAG
  CP.set_var( NP, k ); //decision variables

  SUM = 0.;
  for( unsigned i=0; i<NK; i++ )
    SUM += mc::sqr( Ym[i] - Ymod(Tm[i],k) ); //constraint 
  CP.add_ctr( mc::BASE_OPT::LE, SUM - NLP.get_local_solution().f );
  CP.setup();

  //options
  CP.options.DISPLAY     = 2;
  CP.options.MAXITER     = 5000;
  CP.options.CVATOL      = 1e-6;
  CP.options.CVRTOL      = 1e-6;
  CP.options.BRANCHVAR   = mc::SetInv<CVI>::Options::RGREL;
  CP.options.NODEMEAS    = mc::SetInv<CVI>::Options::MEANWIDTH;
  CP.options.DOMREDMAX   = 10;
  CP.options.DOMREDTHRES = 2e-2;
  CP.options.DOMREDBKOFF = 1e-8;
  CP.options.RELMETH     = mc::NLCP<I>::Options::DRL;
  CP.options.CMODPROP    = 1;
  std::cout << CP;

  CP.solve( Ip );

#if defined(SAVE_RESULTS )
  std::ofstream K_un( "test4b.out", std::ios_base::out );
  CP.output_nodes( K_un );
  K_un.close();
#endif

return 0;
}







