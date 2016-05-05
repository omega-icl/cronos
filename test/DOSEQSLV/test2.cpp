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
  mc::FFGraph IVP;             // DAG of the IVP
  OC->set_dag( &IVP );

  const unsigned int NS = 50;   // Time stages
  double tk[NS+1]; tk[0] = 0.;
  for( unsigned k=0; k<NS; k++ ) tk[k+1] = tk[k] + TF/(double)NS;

  const unsigned NP = NS;      // Number of parameters
  mc::FFVar P[NP];             // Parameters
  for( unsigned int i=0; i<NP; i++ ) P[i].set( &IVP );
  OC->set_parameter( NP, P );

  const unsigned NX = 4;       // Number of states
  mc::FFVar X[NX];             // States
  for( unsigned int i=0; i<NX; i++ ) X[i].set( &IVP );
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

  std::map<unsigned,mc::FFVar> OBJ; OBJ.insert( std::make_pair( NS-1, X[2] ) );
  OC->set_obj( mc::DOSEQSLV_IPOPT::MAX, OBJ );         // objective

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
  OC->options.INTMETH   = mc::DOSEQSLV_IPOPT::Options::MSBDF;   // MSADAMS;
  OC->options.JACAPPROX = mc::DOSEQSLV_IPOPT::Options::CV_DENSE;//CV_DIAG;
  OC->options.NMAX      = 10000;
  OC->options.ASACHKPT  = 10000;
  OC->options.ATOL      = OC->options.ATOLB      = OC->options.ATOLS  = 1e-9;
  OC->options.RTOL      = OC->options.RTOLB      = OC->options.RTOLS  = 1e-9;

  OC->setup();
  Ipopt::ApplicationReturnStatus status = OC->solve( NS, tk, Ip, p0 );

  //if( status == Ipopt::Solve_Succeeded ){
    std::cout << "OC (LOCAL) SOLUTION: " << std::endl;
    std::cout << "  f* = " << OC->solution().f << std::endl;
    for( unsigned int ip=0; ip<NP; ip++ )
      std::cout << "  p*(" << ip << ") = " << OC->solution().p[ip] << std::endl;
  //}

  double*xk[NS+1], f;
  OC->ODESLVS_SUNDIALS::options.DISPLAY = 1;
  OC->states( NS, tk, OC->solution().p, xk, &f );

  return 0;
}
