#include <fstream>
#include <iomanip>
#include "doseqslv_ipopt.hpp"

using mc::sqr;

////////////////////////////////////////////////////////////////////////
// TEST PROBLEM #2 IN FLOUDAS ET AL (NEAR-SINGULAR CONTROL PROBLEM)
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
  Ipopt::SmartPtr<mc::DOSEQSLV_IPOPT> OC = new mc::DOSEQSLV_IPOPT;
  mc::FFGraph DAG;                                                 // DAG
  OC->set_dag( &DAG );

  const unsigned int NS = 4;                                       // Time stages
  double TS[NS+1]; TS[0] = 0.;
  for( unsigned k=0; k<NS; k++ ) TS[k+1] = TS[k] + 1./(double)NS;
  mc::FFVar T; T.set( &DAG );                                      // Time variable
  OC->set_time( NS, TS, &T );

  const unsigned NP = NS;                                          // Parameters
  mc::FFVar P[NP];
  for( unsigned int i=0; i<NP; i++ ) P[i].set( &DAG );
  OC->set_parameter( NP, P );

  const unsigned NX = 3;                                           // States
  mc::FFVar X[NX];
  for( unsigned int i=0; i<NX; i++ ) X[i].set( &DAG );
  OC->set_state( NX, X );

  mc::FFVar RHS[NX*NS];                                            // Dynamics
  for( unsigned k=0; k<NS; k++ ){
    RHS[NX*k+0] = X[1];
    RHS[NX*k+1] = -X[2]*P[k] + 16*T - 8;
    RHS[NX*k+2] = P[k];
  }
  OC->set_differential( NS, NX, RHS );
 
  const unsigned NQ = 1;                                           // States
  mc::FFVar Q[NQ];
  for( unsigned int i=0; i<NQ; i++ ) Q[i].set( &DAG );

  mc::FFVar QUAD[NQ*NS];                                           // Quadratures
  for( unsigned k=0; k<NS; k++ )
    QUAD[NQ*k+0] = sqr(X[0])+sqr(X[1])+5e-4*sqr(X[1]+16.*T-8.-0.1*X[2]*sqr(P[k]));
  OC->set_quadrature( NS, NQ, QUAD, Q );

  mc::FFVar IC[NX];                                                // Initial
  IC[0] = 0.;
  IC[1] = -1.;
  IC[2] = -sqrt(5);
  OC->set_initial( NX, IC );

  std::map<unsigned,mc::FFVar> COST;
  for( unsigned k=0; k<NS; k++ ) COST[k] = Q[0];
  OC->set_obj( mc::BASE_DO::MIN, COST );                            // Cost

  typedef mc::Interval I;
  I Ip[NP];
  double p0[NP];
  for( unsigned int is=0; is<NS; is++ ){
    p0[is] = 0.;  Ip[is]  = I( -4., 10. );
  }

  OC->options.DISPLAY  = 5;
  OC->options.MAXITER  = 100;
  OC->options.TESTDER  = false;
  OC->options.CVTOL    = 1e-7;
  OC->options.GRADIENT = mc::DOSEQSLV_IPOPT::Options::FORWARD; // BACKWARD;
  OC->options.ODESLVS.JACAPPROX = mc::ODESLVS_SUNDIALS::Options::CV_DENSE;//CV_DIAG;
  OC->options.ODESLVS.INTMETH   = mc::ODESLVS_SUNDIALS::Options::MSBDF;   // MSADAMS;
  OC->options.ODESLVS.NMAX      = OC->options.ODESLVS.ASACHKPT = 10000;
  OC->options.ODESLVS.ATOL      = OC->options.ODESLVS.ATOLB = OC->options.ODESLVS.ATOLS = 1e-10;
  OC->options.ODESLVS.RTOL      = OC->options.ODESLVS.RTOLB = OC->options.ODESLVS.RTOLS = 1e-9;

  OC->setup();
  int status = OC->solve( Ip, p0 );

  //if( status == Ipopt::Solve_Succeeded ){
    std::cout << "OC (LOCAL) SOLUTION: " << std::endl;
    std::cout << "  f* = " << OC->solution().f << std::endl;
    for( unsigned int ip=0; ip<NP; ip++ )
      std::cout << "  p*(" << ip << ") = " << OC->solution().p[ip] << std::endl;
  //}

  OC->options.ODESLVS.DISPLAY = 1;
  OC->states( OC->solution().p );

  return 0;
}
