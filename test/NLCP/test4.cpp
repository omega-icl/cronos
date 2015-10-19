#define SAVE_RESULTS    // whether or not to save results to file
#undef USE_PROFIL	// specify to use PROFIL for interval arithmetic
#undef USE_FILIB	// specify to use FILIB++ for interval arithmetic
#undef DEBUG            // whether to output debug information
////////////////////////////////////////////////////////////////////////////////

#include <fstream>
#include "nlcp_gurobi.hpp"

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
typedef mc::CVar<I> CVI;

const unsigned NT = 10;
const double Tm[NT] = { 0.75, 1.5,  2.25, 3.,     6.,    9.,   13.,   17.,   21.,   25. };
const double Ym[NT] = { 7.39, 4.09, 1.74, 0.097, -2.57, -2.71, -2.07, -1.44, -0.98, -0.66 };
const double dYm = 0.5;

////////////////////////////////////////////////////////////////////////////////
// Find all solutions of the system of nonlinear inequalities:
//   yk - 0.5 <= p1*exp(-p2*tk) + p3*exp(-p4*tk) <= yk + 0.5
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

  mc::NLCP_GUROBI<I> CP;
  CP.set_dag( &DAG );
  CP.set_var( NP, P );
  for( unsigned k=0; k<NY; k++ ){
    Y[k] = (58.*P[0]+2.)*exp(-P[1]*Tm[k]) + (29.*P[2]-30.)*exp(-0.5*P[3]*Tm[k]);
    CP.add_ctr( mc::BASE_OPT::LE, Y[k]-Ym[k]-dYm );
    CP.add_ctr( mc::BASE_OPT::GE, Y[k]-Ym[k]+dYm );
  }
  CP.setup();

  //CP.options.MIPFILE     = "test1.lp";
  CP.options.DISPLAY     = 1;
  CP.options.MAXITER     = 5000;
  CP.options.CVATOL      = 1e-6;
  CP.options.CVRTOL      = 1e-6;
  CP.options.BRANCHVAR   = mc::SetInv<CVI>::Options::RGREL;
  CP.options.NODEMEAS    = mc::SetInv<CVI>::Options::LENGTH;
  CP.options.DOMREDMAX   = 10;
  CP.options.DOMREDTHRES = 2e-2;
  CP.options.DOMREDBKOFF = 1e-8;
  CP.options.RELMETH     = mc::NLCP_GUROBI<I>::Options::CHEB;
  CP.options.CMODPROP    = 1;
  std::cout << CP;

  const I Ip[NP] = { I(0.,1.), I(0.,1.), I(0.,1.), I(0.,1.) };

  CP.solve( Ip );

#if defined(SAVE_RESULTS )
  std::ofstream K_un( "test4.out", std::ios_base::out );
  CP.output_nodes( K_un );
  K_un.close();
#endif

  return 0;
}

