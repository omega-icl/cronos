#define USE_PROFIL	    // specify to use PROFIL for interval arithmetic
#undef  USE_FILIB	    // specify to use FILIB++ for interval arithmetic
#undef  USE_DEPS        // whether to use dependents
#define MC__USE_CPLEX   // whether to use CPLEX or GUROBI
#undef  MC__CSEARCH_SHOW_BOXES
#undef  MC__CSEARCH_SHOW_DEPS
#undef  MC__CSEARCH_SHOW_REDUC
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
// TEST PROBLEM #2 IN FLOUDAS ET AL (SINGULAR CONTROL PROBLEM)
// 
//    MIN  z4(1)
//     u                                                     _
//     s.t.  dz1dt = z2,                                      |
//           dz2dt = -z3*u+16*t-8,                            |
//           dz3dt = u,                                       |  t in (0,1]
//           dz4dt = z1^2+z2^2+5e-4*(z2+16*t-8-0.1*z3*u^2)^2,_|
//           z(0) = [ 0, -1, -sqrt(5), 0 ]
//           -4 <= u  <= 10
////////////////////////////////////////////////////////////////////////
int main()
////////////////////////////////////////////////////////////////////////
{
  mc::DOSEQGO<I> OC;
  mc::FFGraph DAG;                                                 // DAG
  OC.set_dag( &DAG );

  const unsigned int NS = 2;                                       // Time stages
  double TS[NS+1]; TS[0] = 0.;
  for( unsigned k=0; k<NS; k++ ) TS[k+1] = TS[k] + 1./(double)NS;
  mc::FFVar T; T.set( &DAG );                                      // Time variable
  OC.set_time( NS, TS, &T );

  const unsigned NP = NS;                                          // Parameters
  mc::FFVar P[NP];
  for( unsigned int i=0; i<NP; i++ ) P[i].set( &DAG );
  OC.set_parameter( NP, P );

  const unsigned NX = 3;                                           // States
  mc::FFVar X[NX];
  for( unsigned int i=0; i<NX; i++ ) X[i].set( &DAG );
  OC.set_state( NX, X );

  mc::FFVar RHS[NX*NS];                                            // Dynamics
  for( unsigned k=0; k<NS; k++ ){
    RHS[NX*k+0] = X[1];
    RHS[NX*k+1] = -X[2]*P[k] + 16*T - 8;
    RHS[NX*k+2] = P[k];
  }
  OC.set_differential( NS, NX, RHS );
 
  const unsigned NQ = 1;                                           // States
  mc::FFVar Q[NQ];
  for( unsigned int i=0; i<NQ; i++ ) Q[i].set( &DAG );

  mc::FFVar QUAD[NQ*NS];                                           // Quadratures
  for( unsigned k=0; k<NS; k++ )
    QUAD[NQ*k+0] = sqr(X[0])+sqr(X[1])+5e-4*sqr(X[1]+16.*T-8.-0.1*X[2]*sqr(P[k]));
  OC.set_quadrature( NS, NQ, QUAD, Q );

  mc::FFVar IC[NX];                                                // Initial
  IC[0] = 0.;
  IC[1] = -1.;
  IC[2] = -sqrt(5);
  OC.set_initial( NX, IC );

  std::map<unsigned,mc::FFVar> COST;
  for( unsigned k=0; k<NS; k++ ) COST[k] = Q[0];
  OC.set_obj( mc::BASE_DO::MIN, COST );                            // Cost

  I Ip[NP];
  double p0[NP];
  for( unsigned int is=0; is<NS; is++ ){
    p0[is] = 0.;  Ip[is]  = I( -4., 10. );
  }

  OC.options.POLIMG.SANDWICH_MAXCUT = 4;
  OC.options.MIPFILE = "";//test4.lp";
  OC.options.MAXITER = 10000;
  OC.options.CMODPROP = OC.options.CMODDEPS = 3;
  OC.options.CMODMIG = 1e-7;
  OC.options.ODEBNDS.DISPLAY = 0;
  OC.options.PREPROC = true;
  OC.options.DISPLAY = 2;
  OC.options.MIPDISPLAY = 0;
  OC.options.ODEBNDS.DISPLAY = 0;
  OC.options.ODEBNDS.ATOL = 1e-10;
  OC.options.ODEBNDS.RTOL = 1e-8;
  OC.options.DOSEQSLV.DISPLAY = 0;
  OC.options.DOSEQSLV.MAXITER = 10;
  OC.options.DOSEQSLV.CVTOL = 1e-5;
  OC.options.DOSEQSLV.GRADIENT = mc::DOSEQSLV_IPOPT::Options::BACKWARD; //FORWARD; //BACKWARD;
  OC.options.DOSEQSLV.ODESLVS.NMAX = OC.options.DOSEQSLV.ODESLVS.ASACHKPT = 5000;
  //OC.options.DOSEQSLV.ODESLVS.ATOL = OC.options.DOSEQSLV.ODESLVS.ATOLB = OC.options.DOSEQSLV.ODESLVS.ATOLS  = 1e-;
  //OC.options.DOSEQSLV.ODESLVS.RTOL = OC.options.DOSEQSLV.ODESLVS.RTOLB = OC.options.DOSEQSLV.ODESLVS.RTOLS  = 1e-9;
  OC.setup();

/*
  int statloc = OC.local( Ip, p0 );
  std::cout << "DO LOCAL SOLUTION: " << statloc << std::endl;
  std::cout << "  f* = " << OC.get_local_solution().f << std::endl;
  for( unsigned int ip=0; ip<NP; ip++ )
    std::cout << "  p*(" << ip << ") = " << OC.get_local_solution().p[ip] << std::endl;

  int statrel = OC.relax( Ip );
  std::cout << "DO RELAXED SOLUTION: " << statrel << std::endl;
  std::cout << "  f* = " << OC.get_objective() << std::endl;
  for( unsigned int ip=0; ip<NP; ip++ )
    std::cout << "  p*(" << ip << ") = " << OC.get_variable(P[ip]) << std::endl;

  unsigned nred;
  int statred = OC.contract( Ip, nred );
  std::cout << "DO CONTRACTED BOUNDS: " << nred << std::endl;
  for( unsigned int ip=0; ip<NP; ip++ )
    std::cout << "  P(" << ip << ") = " << Ip[ip] << std::endl;
*/
  int status = OC.solve( Ip, 0, 0, p0 );
  OC.stats.display();


  return 0;
}
