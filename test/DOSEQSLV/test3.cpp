#include <fstream>
#include <iomanip>
#include "doseqslv_ipopt.hpp"

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
  Ipopt::SmartPtr<mc::DOSEQSLV_IPOPT> OC = new mc::DOSEQSLV_IPOPT;
  mc::FFGraph DAG;             // DAG of the DAG
  OC->set_dag( &DAG );

  //const double X0 = 1., TF = 1., XS = 5., A = 1., B = 2., C = 2.;
  const double X0 = 1., TF = 5., XS = 5., A = 1., B = 2., C = 2.;

  const unsigned int NS = 50;   // Time stages
  double tk[NS+1]; tk[0] = 0.;
  for( unsigned k=0; k<NS; k++ ) tk[k+1] = tk[k] + TF/(double)NS;
  OC->set_time( NS, tk );

  mc::FFVar U[NS];             // Discretized controls
  for( unsigned int i=0; i<NS; i++ ) U[i].set( &DAG );
  OC->set_parameter( NS, U );

  const unsigned NX = 2;       // Number of states
  mc::FFVar X[NX];             // States
  for( unsigned int i=0; i<NX; i++ ) X[i].set( &DAG );
  OC->set_state( NX, X );

  mc::FFVar RHS[NX*NS];        // Right-hand side function
  for( unsigned k=0; k<NS; k++ ){
    //RHS[NX*k+0] = X[0] * ( XS - X[0] - U[k] );
    RHS[NX*k+0] = X[0] * ( (k<NS/2? XS: XS-1.) - X[0] - U[k] );
    RHS[NX*k+1] = A * X[0] + B * U[k] - C * X[0] * U[k];
  }
  OC->set_differential( NS, NX, RHS );
 
  mc::FFVar IC[NX];            // Initial value function
  IC[0] = X0;
  IC[1] = 0.;
  OC->set_initial( NX, IC );

  OC->set_obj( mc::DOSEQSLV_IPOPT::MIN, std::make_pair( NS-1, X[1]/TF ) ); // cost

  typedef mc::Interval I;
  I Iu[NS];
  double u0[NS];
  for( unsigned int is=0; is<NS; is++ ){
    Iu[is] = I( 0., 5. ); u0[is] = mc::Op<I>::mid(Iu[is]);  
  }

  OC->options.DISPLAY   = 5;
  OC->options.MAXITER   = 200;
  OC->options.TESTDER   = false;
  OC->options.CVTOL     = 1e-7;
  OC->options.GRADIENT  = mc::DOSEQSLV_IPOPT::Options::BACKWARD;//FORWARD; // BACKWARD;
  OC->options.ODESLVS.INTMETH   = mc::BASE_SUNDIALS::Options::MSBDF;   // MSADAMS;
  OC->options.ODESLVS.JACAPPROX = mc::BASE_SUNDIALS::Options::CV_DENSE;// CV_DIAG;
  OC->options.ODESLVS.DISPLAY = 0;
  OC->options.ODESLVS.NMAX = OC->options.ODESLVS.ASACHKPT = 4000;
  OC->options.ODESLVS.HMIN = 1e-12;
  OC->options.ODESLVS.ATOL = OC->options.ODESLVS.ATOLB = OC->options.ODESLVS.ATOLS = 1e-8;
  OC->options.ODESLVS.RTOL = OC->options.ODESLVS.RTOLB = OC->options.ODESLVS.RTOLS = 1e-8;

  OC->setup();
  OC->solve( Iu, u0 );

  // Piecewise constant optimal control profile
  std::ofstream direcCON( "test3_CON.dat", std::ios_base::out );
  direcCON << std::scientific << std::setprecision(5);
  direcCON << tk[0] << "  " << OC->solution().p[0] << std::endl;
  for( unsigned int is=0; is<NS; is++ )
    direcCON << tk[is+1] << "  " << OC->solution().p[is] << std::endl;
  direcCON.close();

  // Optimal response trajectories
  OC->options.ODESLVS.DISPLAY   = 1;
  OC->options.ODESLVS.RESRECORD = 20;
  OC->states( OC->solution().p.data() );
  std::ofstream direcSTA( "test3_STA.dat", std::ios_base::out );
  OC->record( direcSTA );
  direcSTA.close();

  return 0;
}
