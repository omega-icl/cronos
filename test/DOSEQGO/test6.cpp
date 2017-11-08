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

//Constant Parameters
const double b1dR = 5.0e3/1.9872;
const double b2dR = 1.0e4/1.9872;
const double a1 = 4.0e3;
const double a2 = 6.2e5;
const double TL = 298.;
const double TU = 398.;//323.;
const double tf = 1.;

////////////////////////////////////////////////////////////////////////
// DENBIGH TEST PROBLEM
////////////////////////////////////////////////////////////////////////
int main()
////////////////////////////////////////////////////////////////////////
{
  mc::DOSEQGO<I> OC;
  mc::FFGraph DAG;                                                 // DAG
  OC.set_dag( &DAG );

  const unsigned int NS = 6;                                       // Time stages
  double TS[NS+1]; TS[0] = 0.;
  for( unsigned k=0; k<NS; k++ ) TS[k+1] = TS[k] + tf/(double)NS;
  mc::FFVar T; T.set( &DAG );                                      // Time variable
  OC.set_time( NS, TS, &T );

  const unsigned NP = NS;                                          // Parameters
  mc::FFVar P[NP];
  for( unsigned int i=0; i<NP; i++ ) P[i].set( &DAG );
  OC.set_parameter( NP, P );

  const unsigned NX = 2;                                           // States
  mc::FFVar X[NX];
  for( unsigned int i=0; i<NX; i++ ) X[i].set( &DAG );
  OC.set_state( NX, X );

  mc::FFVar RHS[NX*NS];                                            // Dynamics
  for( unsigned k=0; k<NS; k++ ){
    //mc::FFVar K1 = a1 * exp( -b1dR * P[k] / TL );
    //mc::FFVar K2 = a2 * exp( -b2dR * P[k] / TL );
    mc::FFVar K1 = a1 * exp( -b1dR / P[k] );
    mc::FFVar K2 = a2 * exp( -b2dR / P[k] );
    RHS[NX*k+0] = - K1 * sqr( X[0] );
    RHS[NX*k+1] =   K1 * sqr( X[0] ) - K2 * X[1];
  }
  OC.set_differential( NS, NX, RHS );

  mc::FFVar IC[NX];                                                // Initial
  IC[0] = 1.;
  IC[1] = 0.;
  OC.set_initial( NX, IC );

  OC.set_obj( mc::BASE_DO::MAX, std::make_pair( NS-1, X[1] ) );    // cost

  I Ip[NP];                                                        // Parameter bounds
  double p0[NP];
  for( unsigned k=0; k<NP; k++ ){
    Ip[k] = I( TL, TU );
    p0[k] = TL;
    //Ip[k] = TL / I( TL, TU );
    //p0[k] = 1.;//mc::Op<I>::mid( Ip[k] ); //1.;
  }

  OC.options.POLIMG.SANDWICH_MAXCUT = 4;
  OC.options.MIPFILE = "";//test6.lp";
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
  OC.options.DOSEQSLV.GRADIENT = mc::DOSEQSLV_IPOPT::Options::FORWARD; //BACKWARD;//
  OC.options.DOSEQSLV.ODESLVS.NMAX = OC.options.DOSEQSLV.ODESLVS.ASACHKPT = 50000;
  OC.options.DOSEQSLV.ODESLVS.ATOL = OC.options.DOSEQSLV.ODESLVS.ATOLB = OC.options.DOSEQSLV.ODESLVS.ATOLS  = 1e-7;
  OC.options.DOSEQSLV.ODESLVS.RTOL = OC.options.DOSEQSLV.ODESLVS.RTOLB = OC.options.DOSEQSLV.ODESLVS.RTOLS  = 1e-9;
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
