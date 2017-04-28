#define USE_PROFIL	    // specify to use PROFIL for interval arithmetic
#undef  USE_FILIB	    // specify to use FILIB++ for interval arithmetic
#undef  USE_DEPS        // whether to use dependents
#define MC__USE_CPLEX   // whether to use CPLEX or GUROBI
#define  MC__CSEARCH_SHOW_BOXES
#undef  MC__CSEARCH_SHOW_DEPS
#define  MC__CSEARCH_SHOW_REDUC
#undef  MC__CSEARCH_SHOW_OUTER
#undef  MC__SBP_SHOW_SCOREBRANCHING
////////////////////////////////////////////////////////////////////////////////

#include <fstream>
#include <iomanip>
#include "doseqgo.hpp"

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

////////////////////////////////////////////////////////////////////////
// MINIMUM TIME CONTROL PROBLEM WITH TERMINAL EQUALITY CONSTRAINTS
// 
//    MIN  tf
//   tf, u, p
//     s.t.   x1(tf) = x1f
//            x2(tf) = x2f                 _
//             dx1dt = x2                   |  t in (0,tf]
//             dx2dt = (u-x1-2*x2)         _|
//             x1(0) = x10, x2(0) = x20
//            umin <= u  <= umax
//           tfmin <= tf <= tfmax
// 
//     where:  x1f   = 1.0
//             x2f   = 0.0
//             x10   = 0.0
//             x20   = p
//             umin  = -0.5
//             umax  =  1.5
//             tfmin = 0.1
//             tfmax = 10.0
////////////////////////////////////////////////////////////////////////
int main()
////////////////////////////////////////////////////////////////////////
{
  mc::DOSEQGO<I> DO;
  mc::FFGraph DAG;                                                 // DAG
  DO.set_dag( &DAG );

  const unsigned int NS = 8;                                       // Time stages
  double tk[NS+1]; tk[0] = 0.;
  for( unsigned k=0; k<NS; k++ ) tk[k+1] = tk[k] + 1./(double)NS;
  //mc::FFVar T; T.set( &DAG );                                    // Time variable
  DO.set_time( NS, tk );//, &T );

  const unsigned NP = 2+NS;                                        // Parameters
  mc::FFVar P[NP];
  for( unsigned int i=0; i<NP; i++ ) P[i].set( &DAG );
  mc::FFVar TF = P[0], X10 = P[1], *U = P+2;
  DO.set_parameter( NP, P );

  const unsigned NX = 2;                                           // States
  mc::FFVar X[NX];
  for( unsigned int i=0; i<NX; i++ ) X[i].set( &DAG );
  DO.set_state( NX, X );

  mc::FFVar RHS[NX*NS];                                            // Dynamics
  for( unsigned k=0; k<NS; k++ ){
    RHS[NX*k+0] = TF * X[1];
    RHS[NX*k+1] = TF * ( U[k]*X[0] - 2.*X[1] );
    //RHS[NX*k+1] = TF * ( U[k]*X[0] - 2.*X[1] ) * T;
  }
  DO.set_differential( NS, NX, RHS );
 
  mc::FFVar IC[NX];                                                // Initial
  IC[0] = 0.;
  IC[1] = X10;
  DO.set_initial( NX, IC );

  DO.set_obj( mc::BASE_DO::MIN, TF );                              // Cost
  DO.add_ctr( mc::BASE_DO::EQ, std::make_pair( NS-1, X[0]-1. ) );  // Terminal constraint
  DO.add_ctr( mc::BASE_DO::EQ, std::make_pair( NS-1, X[1] ) );     // Terminal constraint

  I Ip[NP];
  double p0[NP];
  p0[0] = 3.;   Ip[0] = I( 1., 10. );
  p0[1] = 0.5;  Ip[1] = I( -1., 1. );
  for( unsigned int is=0; is<NS; is++ ){
    p0[2+is] = 0.5;  Ip[2+is]  = I( -0.5, 1.5 );
  }

  DO.options.POLIMG.SANDWICH_MAXCUT = 4;
  DO.options.MIPFILE = "test1.lp";
  DO.options.MAXITER = 10;
  DO.options.CMODPROP =
  DO.options.CMODDEPS = 3;
  DO.options.ODEBNDS.DISPLAY = 0;
  DO.options.PREPROC = true;
  DO.options.DISPLAY = 2;
  DO.options.MIPDISPLAY = 0;
  DO.options.ODEBNDS.DISPLAY = 0;
  DO.options.ODEBNDS.ATOL = 1e-10;
  DO.options.ODEBNDS.RTOL = 1e-9;
  DO.options.DOSEQSLV.DISPLAY = 0;
  DO.options.DOSEQSLV.MAXITER = 100;
  DO.options.DOSEQSLV.CVTOL = 1e-7;
  DO.options.DOSEQSLV.GRADIENT = mc::DOSEQSLV_IPOPT::Options::BACKWARD; //FORWARD; //BACKWARD;
  DO.options.DOSEQSLV.ODESLVS.NMAX = DO.options.DOSEQSLV.ODESLVS.ASACHKPT = 5000;
  DO.options.DOSEQSLV.ODESLVS.ATOL = DO.options.DOSEQSLV.ODESLVS.ATOLB = DO.options.DOSEQSLV.ODESLVS.ATOLS  = 1e-10;
  DO.options.DOSEQSLV.ODESLVS.RTOL = DO.options.DOSEQSLV.ODESLVS.RTOLB = DO.options.DOSEQSLV.ODESLVS.RTOLS  = 1e-9;
  DO.setup();

/*
  int statloc = DO.local( Ip, p0 );
  std::cout << "DO LOCAL SOLUTION: " << statloc << std::endl;
  std::cout << "  f* = " << DO.get_local_solution().f << std::endl;
  for( unsigned int ip=0; ip<NP; ip++ )
    std::cout << "  p*(" << ip << ") = " << DO.get_local_solution().p[ip] << std::endl;

  int statrel = DO.relax( Ip );
  std::cout << "DO RELAXED SOLUTION: " << statrel << std::endl;
  std::cout << "  f* = " << DO.get_objective() << std::endl;
  for( unsigned int ip=0; ip<NP; ip++ )
    std::cout << "  p*(" << ip << ") = " << DO.get_variable(P[ip]) << std::endl;

  unsigned nred;
  int statred = DO.contract( Ip, nred );
  std::cout << "DO CONTRACTED BOUNDS: " << nred << std::endl;
  for( unsigned int ip=0; ip<NP; ip++ )
    std::cout << "  P(" << ip << ") = " << Ip[ip] << std::endl;
*/
  int status = DO.solve( Ip, 0, 0, p0 );


  return 0;
}
