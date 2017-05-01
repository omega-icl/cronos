#include <fstream>
#include <iomanip>
#include "doseqslv_ipopt.hpp"

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
  Ipopt::SmartPtr<mc::DOSEQSLV_IPOPT> OC = new mc::DOSEQSLV_IPOPT;
  mc::FFGraph DAG;             // DAG of the DAG
  OC->set_dag( &DAG );

  const unsigned int NS = 50;   // Time stages
  double tk[NS+1]; tk[0] = 0.;
  for( unsigned k=0; k<NS; k++ ) tk[k+1] = tk[k] + 1./(double)NS;
  //mc::FFVar T; T.set( &DAG );  // Time
  OC->set_time( NS, tk ); //,&T );

  const unsigned NP = 2+NS;    // Number of parameters
  mc::FFVar P[NP];             // Parameters
  for( unsigned int i=0; i<NP; i++ ) P[i].set( &DAG );
  mc::FFVar TF = P[0], X10 = P[1], *U = P+2;
  OC->set_parameter( NP, P );

  const unsigned NX = 2;       // Number of states
  mc::FFVar X[NX];             // States
  for( unsigned int i=0; i<NX; i++ ) X[i].set( &DAG );
  OC->set_state( NX, X );

  mc::FFVar RHS[NX*NS];        // Right-hand side function
  for( unsigned k=0; k<NS; k++ ){
    RHS[NX*k+0] = TF * X[1];
    RHS[NX*k+1] = TF * ( U[k]*X[0] - 2.*X[1] );
    //RHS[NX*k+1] = TF * ( U[k]*X[0] - 2.*X[1] ) * T;
  }
  OC->set_differential( NS, NX, RHS );
 
  mc::FFVar IC[NX];            // Initial value function
  IC[0] = 0.;
  IC[1] = X10;
  OC->set_initial( NX, IC );

  OC->set_obj( mc::DOSEQSLV_IPOPT::MIN, TF );                              // objective
  OC->add_ctr( mc::DOSEQSLV_IPOPT::EQ, std::make_pair( NS-1, X[0]-1. ) );  // constraint #1
  OC->add_ctr( mc::DOSEQSLV_IPOPT::EQ, std::make_pair( NS-1, X[1] ) );     // constraint #2

  typedef mc::Interval I;
  I Ip[NP];
  double p0[NP];
  p0[0] = 6.;   Ip[0] = I( 0.1, 10. );
  p0[1] = 0.5;  Ip[1] = I( -1., 1. );
  for( unsigned int is=0; is<NS; is++ ){
    p0[2+is] = 0.5;  Ip[2+is]  = I( -0.5, 1.5 );
  }

  OC->options.DISPLAY   = 5;
  OC->options.MAXITER   = 100;
  OC->options.TESTDER   = false;
  OC->options.CVTOL     = 1e-7;
  OC->options.GRADIENT  = mc::DOSEQSLV_IPOPT::Options::BACKWARD; //FORWARD; //BACKWARD;
  OC->options.ODESLVS.INTMETH   = mc::BASE_SUNDIALS::Options::MSBDF; // MSADAMS;
  OC->options.ODESLVS.JACAPPROX = mc::BASE_SUNDIALS::Options::CV_DENSE; //CV_DIAG;
  OC->options.ODESLVS.NMAX = 0;
  OC->options.ODESLVS.DISPLAY = 0;
  OC->options.ODESLVS.ATOL = OC->options.ODESLVS.ATOLB = OC->options.ODESLVS.ATOLS = 1e-9;
  OC->options.ODESLVS.RTOL = OC->options.ODESLVS.RTOLB = OC->options.ODESLVS.RTOLS = 1e-9;

  OC->setup();
  OC->solve( Ip, p0 );

  // Piecewise constant optimal control profile
  std::ofstream direcCON( "test1_CON.dat", std::ios_base::out );
  direcCON << std::scientific << std::setprecision(5);
  double t = 0.;
  direcCON << t << "  " << OC->solution().p[2] << std::endl;
  for( unsigned int is=0; is<NS; is++ ){
    t += OC->solution().p[0]/NS;
    direcCON << t << "  " << OC->solution().p[2+is] << std::endl;
  }
  direcCON.close();

  // Optimal response trajectories
  OC->options.ODESLVS.DISPLAY   = 1;
  OC->options.ODESLVS.RESRECORD = 20;
  OC->states( OC->solution().p );
  std::ofstream direcSTA( "test1_STA.dat", std::ios_base::out );
  OC->record( direcSTA );
  direcSTA.close();

  return 0;
}
