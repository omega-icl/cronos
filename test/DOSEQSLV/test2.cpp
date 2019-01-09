#include <fstream>
#include <iomanip>
#include "doseqslv_ipopt.hpp"

const double cX0   = 1.0;  // [g/L]
const double V0    = 2.0;  // [L]
const double cSin  = 20.0; // [g/L]
const double mumax = 0.53; // [/h]
const double KX    = 1.2;  // [g/L]
const double KI    = 22.0; // [g/L]
const double numax = 0.5;  // [/h]
const double KP    = 0.1;  // [g/L]
const double YX    = 0.4;  // [-]
const double YP    = 1.0;  // [-]
const double TF    = 10.0; // [h]

////////////////////////////////////////////////////////////////////////
// MINIMUM TIME CONTROL PROBLEM WITH TERMINAL EQUALITY CONSTRAINTS
// 
//     MIN    cP(T)·V(T)
// F(t), 0<t<T                                                                            _
//     s.t.    dcXdt = mu(cS(t))·cX(t) - F(t)/V(t)·cX(t)                                   |
//             dcSdt = -mu(cS(t))·cX(t)/YX - nu(cS(t))·cX(t)/YP + F(t)/V(t)·(cSin-cS(t))   | t in (0,T]
//             dcPdt = nu(cS(t))·cX(t) - F(t)/V(t)·cP(t)                                   |
//             dVdt  = F(t)                                                               _|
//             cX(0) = cX0, cS(0) = cP(0) = 0, V(0) = V0
//             0 <= F(t) <= 1
// 
//     with:   mu(S) := mumax·S/(S+KX+S^2/KI)
//             nu(S) := numax·S/(S+KP)
////////////////////////////////////////////////////////////////////////
int main()
////////////////////////////////////////////////////////////////////////
{
  Ipopt::SmartPtr<mc::DOSEQSLV_IPOPT> OC = new mc::DOSEQSLV_IPOPT;
  mc::FFGraph DAG;             // DAG of the DAG
  OC->set_dag( &DAG );

  const unsigned int NS = 100;   // Time stages
  double tk[NS+1]; tk[0] = 0.;
  for( unsigned k=0; k<NS; k++ ) tk[k+1] = tk[k] + TF/(double)NS;
  OC->set_time( NS, tk );

  const unsigned NP = NS;      // Number of parameters
  mc::FFVar P[NP];             // Parameters
  for( unsigned int i=0; i<NP; i++ ) P[i].set( &DAG );
  OC->set_parameter( NP, P );

  const unsigned NX = 4;       // Number of states
  mc::FFVar X[NX];             // States
  for( unsigned int i=0; i<NX; i++ ) X[i].set( &DAG );
  OC->set_state( NX, X );

  mc::FFVar RHS[NX*NS];        // Right-hand side function
  mc::FFVar MU = mumax*X[1]/(X[1]+KX+X[1]*X[1]/KI);
  mc::FFVar NU = numax*X[1]/(X[1]+KP);
  for( unsigned k=0; k<NS; k++ ){
    RHS[NX*k+0] = MU*X[0] - P[k]/X[3]*X[0];
    RHS[NX*k+1] = -MU*X[0]/YX - NU*X[0]/YP + P[k]/X[3]*(cSin-X[1]);
    RHS[NX*k+2] = NU*X[0] - P[k]/X[3]*X[2];
    RHS[NX*k+3] = P[k];
  }
  OC->set_differential( NS, NX, RHS );
 
  mc::FFVar IC[NX];            // Initial value function
  IC[0] = cX0;
  IC[1] = 0.;
  IC[2] = 0.;
  IC[3] = V0;
  OC->set_initial( NX, IC );

  OC->set_obj( mc::DOSEQSLV_IPOPT::MAX, std::make_pair( NS-1, X[2] ) );  // cost

  typedef mc::Interval I;
  I Ip[NP];
  double p0[NP];
  for( unsigned int is=0; is<NS; is++ ){
    p0[is] = 0.5;  Ip[is]  = I( 0., 1. );
  }

  OC->options.DISPLAY   = 5;
  OC->options.MAXITER   = 100;
  OC->options.TESTDER   = false;
  OC->options.CVTOL     = 1e-7;
  OC->options.GRADIENT  = mc::DOSEQSLV_IPOPT::Options::BACKWARD; //FORWARD;
  OC->options.ODESLVS.INTMETH   = mc::BASE_SUNDIALS::Options::MSBDF;   // MSADAMS;
  OC->options.ODESLVS.JACAPPROX = mc::BASE_SUNDIALS::Options::CV_DENSE;//CV_DIAG;
  OC->options.ODESLVS.NMAX = OC->options.ODESLVS.ASACHKPT = 10000;
  OC->options.ODESLVS.ATOL = OC->options.ODESLVS.ATOLB = OC->options.ODESLVS.ATOLS  = 1e-9;
  OC->options.ODESLVS.RTOL = OC->options.ODESLVS.RTOLB = OC->options.ODESLVS.RTOLS  = 1e-9;

  OC->setup();
  OC->solve( Ip, p0 );

  // Piecewise constant optimal control profile
  std::ofstream direcCON( "test2_CON.dat", std::ios_base::out );
  direcCON << std::scientific << std::setprecision(5);
  direcCON << tk[0] << "  " << OC->solution().p[0] << std::endl;
  for( unsigned int is=0; is<NS; is++ )
    direcCON << tk[is+1] << "  " << OC->solution().p[is] << std::endl;
  direcCON.close();

  // Optimal response trajectories
  OC->options.ODESLVS.DISPLAY   = 1;
  OC->options.ODESLVS.RESRECORD = 20;
  OC->states( OC->solution().p.data() );
  std::ofstream direcSTA( "test2_STA.dat", std::ios_base::out );
  OC->record( direcSTA );
  direcSTA.close();

  return 0;
}
