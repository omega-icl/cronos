#define SAVE_RESULTS    // whether or not to save results to file
#undef USE_PROFIL	// specify to use PROFIL for interval arithmetic
#undef USE_FILIB	// specify to use FILIB++ for interval arithmetic
#undef DEBUG            // whether to output debug information
////////////////////////////////////////////////////////////////////////////////

#include <fstream>
#include "nlegpe.hpp"

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
  //const unsigned NP = 3, NY = 2, NT = 28;
  const unsigned NP = 4, NY = 2, NT = 28;
  mc::FFVar P[NP], Y[NY], Ir;
  Ir.set( &DAG );
  for( unsigned i=0; i<NP; i++ ) P[i].set( &DAG );
  //const double tau = 5.5e-3, K = 3.57e-2, thetaL = 8.2e-2, thetaH = 1.8e-2;
  const double tau = 5.5e-3, thetaL = 8.2e-2, thetaH = 1.8e-2;
  mc::FFVar K = P[3];
  Y[0] = P[0]*Ir/(1.+tau*P[1]*pow(thetaL,P[2])*Ir+K*tau*sqr(P[1]*pow(thetaL,P[2])*Ir));
  Y[1] = P[0]*Ir/(1.+tau*P[1]*pow(thetaH,P[2])*Ir+K*tau*sqr(P[1]*pow(thetaH,P[2])*Ir));

  mc::NLEGPE<I> problem;
  problem.set_dag( &DAG );
  problem.set_indep( &Ir );
  problem.set_par( NP, P );
  problem.set_dep( NY, Y );

  problem.options.SETINV.DISPLAY = 1;
  problem.options.SETINV.MAX_NODES = 20000;
  problem.options.SETINV.ABSOLUTE_TOLERANCE = 1e-6;
  problem.options.SETINV.RELATIVE_TOLERANCE = 1e-6;
  problem.options.SETINV.BRANCHING_VARIABLE_CRITERION = mc::SetInv<I>::Options::RGABS;
  problem.options.SETINV.MEASURE = mc::SetInv<I>::Options::LENGTH;

  problem.options.OUTPUT_BOUND = mc::NLEGPE<I>::Options::CM;
  problem.options.CM_ORDER     = 2;
  problem.options.OUTRED_MAX   = 10;
  problem.options.OUTRED_THRES = 2e-2;
  problem.options.OUTRED_TOL   = 1e-9;
  problem.options.INRED_MAX    = 0;
  problem.options.INRED_THRES  = 1e-2;

  const double Im[NT] = {  20.,  20.,  40.,  40.,  60.,  60.,  80.,  80., 100., 100.,
                          150., 150., 200., 200., 300., 300., 400., 400., 600., 600.,
                          800., 800., 1000., 1000., 1200., 1200., 1400., 1400. };
  const double Ym[NT] = { 0.31417,  .31725, 0.61492, 0.62856, 0.90003, 0.93329,
                          1.16775, 1.23086, 1.41685, 1.52077, 1.95433, 2.20916,
                          2.37109, 2.84152, 2.88988, 3.92743, 3.10166, 4.77558,
                          3.05590, 5.8466,  2.79247, 6.30032, 2.5038,  6.37668, 
                          2.24271, 6.24496, 2.01873, 6.00915 };
  const unsigned iy[NT] = { 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1,
                            0, 1, 0, 1, 0, 1, 0, 1, 0, 1 };
  const double dy_rel = 0.01;

  //const I Ip[NP] = { I(0.015,0.017), I(0.4,0.6), I(0.4,0.6) };
  const I Ip[NP] = { I(0.015,0.017), I(0.4,0.6), I(0.4,0.6), I(1e-2,6e-2),  };

  typename std::list<mc::NLEGPE<I>::Data> Iym;
  for( unsigned int k=0; k<NT; k++ )
    Iym.push_back( typename mc::NLEGPE<I>::Data( Ym[k]*(1-dy_rel), Ym[k]*(1+dy_rel), iy[k], Im[k] ) );

  problem.solve( Ip, Iym, std::cout );

#if defined(SAVE_RESULTS )
  ofstream K_un( "undetermined.out", ios_base::out );
  problem.output_nodes( K_un ); //, true );
  K_un.close();
#endif

  return 0;
}
