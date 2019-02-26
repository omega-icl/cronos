#define SAVE_RESULTS    // whether or not to save results to file
#define USE_PROFIL	    // specify to use PROFIL for interval arithmetic
#undef  USE_FILIB	    // specify to use FILIB++ for interval arithmetic

#define MC__USE_CPLEX   // whether to use CPLEX or GUROBI
#undef  MC__CSEARCH_SHOW_BOXES
#undef  MC__CSEARCH_SHOW_REDUC
#undef  MC__CSEARCH_SHOW_INCLUSION
#undef  MC__CSEARCH_PAUSE_INFEASIBLE

////////////////////////////////////////////////////////////////////////////////

#include <fstream>
#include "nlcp.hpp"

#ifdef USE_PROFIL
  #include "mcprofil.hpp"
  typedef INTERVAL I;
#else
  #ifdef USE_FILIB
    #include "mcfilib.hpp"
    typedef filib::interval<double> I;
  #else
    #include "interval.hpp"
    typedef mc::Interval I;
  #endif
#endif

const unsigned NT = 10;
const double Tm[NT] = { 0.75, 1.5,  2.25, 3.,     6.,    9.,   13.,   17.,   21.,   25. };
const double Ym[NT] = { 7.39, 4.09, 1.74, 0.097, -2.57, -2.71, -2.07, -1.44, -0.98, -0.66 };
const double dYm = 0.5;

////////////////////////////////////////////////////////////////////////////////
// Find all solutions of the system of nonlinear inequalities:
//   yk - 0.5 <= (58.*p1+2.)*exp(-p2*tk) + (29.*p3-30.)*exp(-0.5*p4*tk) <= yk + 0.5
//   for 10 measurement pairs (t1,y1),...,(t10,y10) as follows:
//      tk := [ 0.75, 1.5,  2.25, 3.,     6.,    9.,   13.,   17.,   21.,   25. ]
//      yk := [ 7.39, 4.09, 1.74, 0.097, -2.57, -2.71, -2.07, -1.44, -0.98, -0.66 ]
//   and (p1,...,p4) in [0,1]^4
////////////////////////////////////////////////////////////////////////////////

int main() 
{
  mc::FFGraph DAG;
  const unsigned NP = 4, NY = NT;
  mc::FFVar P[NP], Y[NY];
  for( unsigned i=0; i<NP; i++ ) P[i].set( &DAG );

  mc::NLCP<I> CP;
  CP.set_dag( &DAG );
  CP.set_var( NP, P );
  for( unsigned k=0; k<NY; k++ ){
    Y[k] = (58.*P[0]+2.)*exp(-P[1]*Tm[k]) + (29.*P[2]-30.)*exp(-0.5*P[3]*Tm[k]);
    CP.add_ctr( mc::BASE_OPT::LE, Y[k]-Ym[k]-dYm );
    CP.add_ctr( mc::BASE_OPT::GE, Y[k]-Ym[k]+dYm );
  }

  //CP.options.MIPFILE     = "test4.lp";
  CP.options.DISPLAY     = 2;
  CP.options.MAXITER     = 1000;
  CP.options.CVTOL       = 1e-3;
  CP.options.BRANCHVAR   = mc::SBP<I>::Options::RGREL;
  CP.options.NODEMEAS    = mc::SBP<I>::Options::RELMAXLEN;
  CP.options.DOMREDMAX   = 10;
  CP.options.DOMREDTHRES = 5e-2;
  CP.options.RELMETH     = mc::NLCP<I>::Options::CHEB;
  CP.options.CMODPROP    = 3;
  std::cout << CP;

  const I Ip[NP] = { I(0.,1.), I(0.,1.), I(0.,1.), I(0.,1.) };
  CP.setup( Ip );
  CP.solve( Ip );
  CP.stats.display();

#if defined(SAVE_RESULTS )
  std::ofstream K_bnd( "test4_bnd.out", std::ios_base::out );
  std::ofstream K_inn( "test4_inn.out", std::ios_base::out );
  CP.output_nodes( K_bnd, K_inn );
  K_bnd.close();
  K_inn.close();
#endif

  return 0;
}

