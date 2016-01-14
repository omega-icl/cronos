#include <fstream>
#include <iomanip>
#include "nlgo_gurobi.hpp"
#include "interval.hpp"
typedef mc::Interval I;

// Parameter estimation problem (microalgae PI curve)
////////////////////////////////////////////////////////////////////////

template <class T>
T Pmod
( const T*k, double q, double Ir )
{
  const double tau = 5.5e-3, K = 3.57e-2;
  return k[0]*Ir/(1.+tau*k[1]*pow(q,k[2])*Ir+K*tau*sqr(k[1]*pow(q,k[2])*Ir));
}

const unsigned NT = 22;
const double Im_50[NT] = { 10.9249, 18.16, 36.46, 73.12, 29.06, 69.3541, 40, 43.54, 102.196, 94.78, 112.995, 142.184, 207.9, 255.617, 310.528, 387.566, 416.857, 611.368, 640.827, 754.596, 971.193, 1059.39 }; // [µE/m2·s]
const double Pm_50[NT] = { 0.197, 0.412, 0.51955, 0.609, 0.6448, 0.806, 0.788139, 1.074, 1.2006, 1.3617, 1.6483, 2.00669, 2.7235, 2.706, 2.99324, 3.0478, 3.19136, 3.15773, 2.94322, 2.92659, 2.74999, 2.4824 }; // [gchl/gC·h]
const double Im_1200[NT] = { 21.98, 47.521, 36.41, 50.97, 98.6, 83.83, 160.676, 193.358, 167.544, 222.11, 313.827, 269.622, 419.289, 610.021, 602.581, 624.616, 668.694, 866.909, 998.867, 1186.03, 1310.83 ,1347.39 }; // [µE/m2·s]
const double Pm_1200[NT] = { 0.0897199, 0.41226, 0.627, 0.878, 1.0394, 1.2183, 1.703, 2.4369, 2.705, 3.726, 3.781, 4.1386, 5.82, 6.0222, 6.23695, 6.201, 6.11237, 6.00717, 6.33, 6.315, 6.2628, 6.5496 }; // [gchl/gC·h]
const double q_50 = 8.2e-2, q_1200 = 1.8e-2; // [gchl/gC]
const double ePm = 1e-2; // relative error

////////////////////////////////////////////////////////////////////////
int main()
////////////////////////////////////////////////////////////////////////
{
  mc::FFGraph DAG;
  const unsigned NK = 3, N = NK+NT; mc::FFVar k[N];
  for( unsigned i=0; i<N; i++ ) k[i].set( &DAG );

  mc::NLGO_GUROBI<I> NLP;
  NLP.options.POLIMG.SANDWICH_MAXCUT = 7;
  NLP.options.POLIMG.SANDWICH_ATOL   = NLP.options.POLIMG.SANDWICH_RTOL  = 1e-5;
  //NLP.options.CVATOL = NLP.options.CVRTOL = 1e-5;
  //NLP.options.MIPABSGAP = NLP.options.MIPRELGAP = 1e-9;
  //NLP.options.MAXITER = 0;

/*
  NLP.options.POLIMG.SANDWICH_MAXCUT = 10;
  NLP.options.POLIMG.SANDWICH_ATOL   = NLP.options.POLIMG.SANDWICH_RTOL  = 1e-5;
  NLP.options.POLIMG.BREAKPOINT_TYPE = mc::PolImg<I>::Options::BIN;//SOS2;//NONE;
  NLP.options.POLIMG.DCDECOMP_SCALE  = false;//true;
  //NLP.options.MIPFILE = "test4.lp";
  NLP.options.MIPDISPLAY = 0;
  NLP.options.DISPLAY = 2;
  NLP.options.NLPSLV.DISPLAY = 0;
  NLP.options.MAXITER = 21;
  NLP.options.PREPROC = true;
  NLP.options.DOMREDMAX = 10;
  NLP.options.PRESOS2BIGM = -1;
  NLP.options.CVATOL    = NLP.options.CVRTOL    = 1e-5;
*/
  NLP.set_dag( &DAG );  // DAG
  NLP.set_var( N, k ); // decision variables

  mc::FFVar OBJ=0.;
  for( unsigned i=0; i<NT; i++ ){
    NLP.add_ctr( mc::BASE_NLP::GE, k[NK+i] - mc::sqr( Pmod( k, q_1200, Im_1200[i] ) - Pm_1200[i]*(1.+ePm) ) );
    NLP.add_ctr( mc::BASE_NLP::GE, k[NK+i] - mc::sqr( Pmod( k, q_1200, Im_1200[i] ) - Pm_1200[i]*(1.-ePm) ) );
    OBJ += k[NK+i];
  }
  NLP.set_obj( mc::NLGO_GUROBI<I>::MIN, OBJ );

  //const mc::FFVar* dLSQdk = DAG.BAD( 1, &LSQ, N, k );
  //for( unsigned i=0; i<N; i++ )
  //  NLP.add_ctr( mc::BASE_NLP::EQ, dLSQdk[i] );
  NLP.setup();
  std::cout << NLP;

  typedef mc::Interval I;
  I Ik[N];
  Ik[0] = I(0.,0.02); Ik[1] = I(0.,1.); Ik[2] = I(0.,1.);
  for( unsigned i=0; i<NT; i++ ) Ik[NK+i] = I(0.,1.e2);

  double k0[N];
  k0[0] = 0.0160; k0[1] = 0.492; k0[2] = 0.469;
  for( unsigned i=0; i<NT; i++ ) k0[NK+i] = 0.;

  // Local optimization
  NLP.options.NLPSLV.DISPLAY = 5;
  NLP.options.NLPSLV.MAXITER = 100;
  NLP.local( Ik, k0 );
  std::cout << std::setprecision(4);
  for( unsigned i=0; i<N; i++ )
    std::cout << "  k(" << i << ") = " << NLP.get_local_solution().p[i];
  std::cout << std::endl;

  // Global optimization using SBB
  NLP.options.DISPLAY = 2;
  //NLP.options.MIPFILE = "test6.lp";
  NLP.options.NLPSLV.DISPLAY = 0;
  NLP.options.CSALGO  = mc::NLGO_GUROBI<I>::Options::SBB;
  NLP.options.RELMETH = mc::NLGO_GUROBI<I>::Options::HYBRID;//CHEB;
  NLP.options.CMODPROP = 1;
  NLP.solve( Ik, 0, k0 );

  return 0;
}
