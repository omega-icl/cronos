#include <fstream>
#include <iomanip>
#include "doseqslv_ipopt.hpp"

char fname[50];
std::ofstream ofile;

////////////////////////////////////////////////////////////////////////
// MINIMUM TIME CONTROL PROBLEM WITH TERMINAL EQUALITY CONSTRAINTS
// 
//    MIN  1/T \int_0^T a·x(t) + b·u(t) - c·x(t)·u(t)
//     u
//     s.t.   dx/dt = x(t)·[xs-x(t)-u(t)],  t in (0,T]
//            x(0) = x0
//            0 <= u <= U
// 
//     where:  a = 1
//             b = c = 2
//             x0 = 1
//             xs = 5
//             U  = 5 
////////////////////////////////////////////////////////////////////////
int main()
////////////////////////////////////////////////////////////////////////
{
  const double X0 = 1., TF = 2., XS = 5., A = 1., B = 2., C = 2.;
  mc::FFGraph DAG;             // DAG of the DAG

  Ipopt::SmartPtr<mc::DOSEQSLV_IPOPT> OC = new mc::DOSEQSLV_IPOPT;
  OC->set_dag( &DAG );

  const unsigned int NS = 3;   // Time stages
  double tk[NS+1]; tk[0] = 0.;
  for( unsigned k=0; k<NS; k++ ) tk[k+1] = tk[k] + 1;
  OC->set_time( NS, tk );

  const unsigned int NP = 2*NS;   // Decision variables [ U1, U2, U3, DT1, DT2, DT3 ]
  mc::FFVar P[NP], *U = P, *DT = P+NS;  // Discretized controls
  for( unsigned int i=0; i<NP; i++ ) P[i].set( &DAG );
  OC->set_parameter( NP, P );

  const unsigned NX = 2;       // Number of states
  mc::FFVar X[NX];             // States
  for( unsigned int i=0; i<NX; i++ ) X[i].set( &DAG );
  OC->set_state( NX, X );

  mc::FFVar RHS[NX*NS];        // Right-hand side function
  for( unsigned k=0; k<NS; k++ ){
    RHS[NX*k+0] = DT[k] * X[0] * ( XS - X[0] - U[k] );
    RHS[NX*k+1] = DT[k] * ( A * X[0] + B * U[k] - C * X[0] * U[k] );
  }
  OC->set_differential( NS, NX, RHS );
 
  mc::FFVar IC[NX];            // Initial value function
  IC[0] = X0;
  IC[1] = 0.;
  OC->set_initial( NX, IC );

  OC->set_obj( mc::DOSEQSLV_IPOPT::MIN, std::make_pair( NS-1, X[1] ) ); // objective
  OC->add_ctr( mc::DOSEQSLV_IPOPT::EQ, DT[0]+DT[1]+DT[2]-TF ); // constraint

  typedef mc::Interval I;
  I Ip[NP];
  double p0[NP];
  for( unsigned int is=0; is<NS; is++ ){
    Ip[is] = I( 0., 5. ); p0[is] = mc::Op<I>::mid(Ip[is]);  
  }
  Ip[1] = XS/2.-(B-A)/(2.*C); p0[1] = mc::Op<I>::mid(Ip[1]);  
  for( unsigned int is=0; is<NS; is++ ){
    Ip[NS+is] = I( TF/20., TF ); p0[NS+is] = mc::Op<I>::mid(Ip[NS+is]);  
  }

  OC->options.DISPLAY   = 5;
  OC->options.MAXITER   = 100;
  OC->options.TESTDER   = false;
  OC->options.CVTOL     = 1e-7;
  OC->options.GRADIENT  = mc::DOSEQSLV_IPOPT::Options::FORWARD; // BACKWARD;
  OC->options.ODESLVS.INTMETH   = mc::BASE_SUNDIALS::Options::MSBDF;   // MSADAMS;
  OC->options.ODESLVS.JACAPPROX = mc::BASE_SUNDIALS::Options::CV_DENSE;// CV_DIAG;
  OC->options.ODESLVS.DISPLAY = 0;
  OC->options.ODESLVS.NMAX = OC->options.ODESLVS.ASACHKPT = 4000;
  OC->options.ODESLVS.HMIN = 1e-12;
  OC->options.ODESLVS.ATOL = OC->options.ODESLVS.ATOLB = OC->options.ODESLVS.ATOLS = 1e-8;
  OC->options.ODESLVS.RTOL = OC->options.ODESLVS.RTOLB = OC->options.ODESLVS.RTOLS = 1e-8;

  OC->setup();
  OC->solve( Ip, p0 );
 
  // Piecewise constant optimal control profile
  std::ofstream direcCON( "test3b_CON.dat", std::ios_base::out );
  direcCON << std::scientific << std::setprecision(5);
  double t = 0.;
  direcCON << t << "  " << OC->solution().p[0] << std::endl;
  for( unsigned int is=0; is<NS; is++ ){
    t += OC->solution().p[NS+is];
    direcCON << t << "  " << OC->solution().p[is] << std::endl;
  }
  direcCON.close();

  // Optimal response trajectories
  OC->options.ODESLVS.DISPLAY   = 1;
  OC->options.ODESLVS.RESRECORD = 50;
  OC->states( OC->solution().p );
  std::ofstream direcSTA( "test3b_STA.dat", std::ios_base::out );
  OC->record( direcSTA );
  direcSTA.close();

  return 0;
}
