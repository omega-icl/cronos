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
  const unsigned NP = 4, NY = 1, NT = 10;
  mc::FFVar P[NP], Y[NY], T;
  T.set( &DAG );
  for( unsigned i=0; i<NP; i++ ) P[i].set( &DAG );
  //Y[0] = P[0]*exp(-P[1]*T) + P[2]*exp(-P[3]*T);
  Y[0] = (58.*P[0]+2.)*exp(-P[1]*T) + (29.*P[2]-30.)*exp(-0.5*P[3]*T);

  mc::NLEGPE<I> problem;
  problem.set_dag( &DAG );
  problem.set_indep( &T );
  problem.set_par( NP, P );
  problem.set_dep( NY, Y );

  problem.options.SETINV.DISPLAY = 1;
  problem.options.SETINV.MAX_NODES = 5000;
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

  const double t[NT] = { 0.75, 1.5,  2.25, 3.,     6.,    9.,   13.,   17.,   21.,   25. };
  const double y[NT] = { 7.39, 4.09, 1.74, 0.097, -2.57, -2.71, -2.07, -1.44, -0.98, -0.66 };
  const double dy_abs = 0.5;

  //const I Ip[NP] = { I(2.,60.), I(0.,1.), I(-30.,-1.), I(0.,0.5) };
  const I Ip[NP] = { I(0.,1.), I(0.,1.), I(0.,1.), I(0.,1.) };

  typename std::list<mc::NLEGPE<I>::Data> Iym;
  for( unsigned int k=0; k<NT; k++ )
    Iym.push_back( typename mc::NLEGPE<I>::Data( y[k]-dy_abs, y[k]+dy_abs, 0, t[k] ) );

  problem.solve( Ip, Iym, std::cout );

#if defined(SAVE_RESULTS )
  ofstream K_un( "undetermined.out", ios_base::out );
  problem.output_nodes( K_un ); //, true );
  K_un.close();
#endif

  return 0;
}
