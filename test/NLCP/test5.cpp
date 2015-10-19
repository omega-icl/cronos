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

////////////////////////////////////////////////////////////////////////
// Parameter estimation problem (microalgae PI curve)

template <class T>
T Pmod
( const T*k, double q, double Ir )
{
  //const double tau = 5.5e-3, K = 3.57e-2;
  //return k[0]*Ir/(1.+tau*k[1]*pow(q,k[2])*Ir+K*tau*sqr(k[1]*pow(q,k[2])*Ir));
  const double tau = 5.5e-3;
  return k[0]*Ir/(1.+tau*k[1]*pow(q,k[2])*Ir+k[3]*tau*sqr(k[1]*pow(q,k[2])*Ir));
}

const unsigned NT = 14;
const double Im_50[NT] = {   20.,   40.,   60.,   80.,  100.,  150.,  200.,  300., 400.,
                            600.,  800., 1000., 1200., 1400. }; // [µE/m2·s]
const double Pm_50[NT] = { 0.31417, 0.61492, 0.90003, 1.16775, 1.41685, 1.95433, 2.37109,
                           2.88988, 3.10166, 3.05590, 2.79247, 2.50380, 2.24271, 2.01873 }; // [gchl/gC·h]
const double Im_1200[NT] = {   20.,   40.,   60.,   80.,  100.,  150.,  200.,  300., 400.,
                              600.,  800., 1000., 1200., 1400. }; // [µE/m2·s]
const double Pm_1200[NT] = { 0.31725, 0.62856, 0.93329, 1.23086, 1.52077, 2.20916, 2.84152,
                             3.92743, 4.77558, 5.84660, 6.30032, 6.37668, 6.24496, 6.00915 }; // [gchl/gC·h]
const double q_50 = 8.2e-2, q_1200 = 1.8e-2; // [gchl/gC]
const double dPm_rel = 0.01;

////////////////////////////////////////////////////////////////////////////////

int main() 
{
  mc::FFGraph DAG;
  //const unsigned NP = 3, NY = 2;
  const unsigned NP = 4, NY = 2;
  mc::FFVar P[NP], Y[NY];
  for( unsigned i=0; i<NP; i++ ) P[i].set( &DAG );

  mc::NLCP_GUROBI<I> CP;
  CP.set_dag( &DAG );
  CP.set_var( NP, P );
  for( unsigned k=0; k<NT; k++ ){
    Y[0] = Pmod( P, q_50,   Im_50[k]   );
    Y[1] = Pmod( P, q_1200, Im_1200[k] );
    CP.add_ctr( mc::BASE_OPT::LE, Y[0]-Pm_50[k]*(1+dPm_rel) );
    CP.add_ctr( mc::BASE_OPT::GE, Y[0]-Pm_50[k]*(1-dPm_rel) );
    CP.add_ctr( mc::BASE_OPT::LE, Y[1]-Pm_1200[k]*(1+dPm_rel) );
    CP.add_ctr( mc::BASE_OPT::GE, Y[1]-Pm_1200[k]*(1-dPm_rel) );
  }
  CP.setup();

  //CP.options.MIPFILE     = "test5.lp";
  CP.options.DISPLAY     = 2;
  CP.options.MAXITER     = 10000;
  CP.options.CVATOL      = 1e-6;
  CP.options.CVRTOL      = 1e-6;
  CP.options.BRANCHVAR   = mc::SetInv<CVI>::Options::RGREL;
  CP.options.NODEMEAS    = mc::SetInv<CVI>::Options::LENGTH;
  CP.options.DOMREDMAX   = 10;
  CP.options.DOMREDTHRES = 2e-2;
  CP.options.DOMREDBKOFF = 1e-8;
  CP.options.RELMETH     = mc::NLCP_GUROBI<I>::Options::CHEB;
  CP.options.CMODPROP    = 2;
  CP.options.CMODCUTS    = 1;
  std::cout << CP;

  //const I Ip[NP] = { I(0.015,0.017), I(0.4,0.6), I(0.4,0.6) };
  const I Ip[NP] = { I(0.015,0.017), I(0.4,0.6), I(0.4,0.6), I(1e-2,6e-2),  };

  CP.solve( Ip );

#if defined(SAVE_RESULTS )
  std::ofstream K_un( "test5.out", std::ios_base::out );
  CP.output_nodes( K_un );
  K_un.close();
#endif

  return 0;
}

