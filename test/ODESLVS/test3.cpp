#define SAVE_RESULTS		// <- Whether to save bounds to file

#include "odeslv_cvodes.hpp"

using std::cout;
using std::endl;

int main()
{

  /////////////////////////////////////////////////////////////////////////
  // Define NIFTE dynamics - NTP model

  mc::FFGraph NIFTE;  // DAG describing the problem

  // Constant parameters  
  const double P0     = 1.013e5;  // Pa
  const double U0     = 8.00e-4;  // m3/s
  const double tau    = 5.;       // s
  const double Lambda = 330.;     // -
  const double K      = 1.37;	  // -, control gain
  const double psi    = 1.39e5;   // Pa, for converting P_d to P_th
  const double chi    = 3.28e-3;  // /Pa, for converting P_d to P_th

  double t0 = 0., tf = 5. ;       // Time span

  const unsigned NP = 10;  // Number of parameters
  const unsigned NX = 5;   // Number of states
  const unsigned NQ = 3;   // Number of state quadratures

  mc::FFVar P[NP];  // Parameters
  for( unsigned int i=0; i<NP; i++ ) P[i].set( &NIFTE );
  
  // Design parameters: R_f, R_l, R_th, L_d, L_f, L_l, L_p, C_ad, C_d, C_p
  double p[NP] = { 2.13e6, 4.08e6, 5.02e8, 1.80e5, 4.74e6, 1.27e7, 3.77e5, 1.76e-9, 7.35e-8, 7.43e-8 };

  mc::FFVar X[NX];  // States
  for( unsigned int i=0; i<NX; i++ ) X[i].set( &NIFTE );
  
  mc::FFVar &P_ad = X[0], &P_d = X[1], &P_p = X[2], &U_f = X[3], &U_p = X[4]; 
  mc::FFVar &R_f = P[0], &R_l = P[1], &R_th = P[2], &L_d = P[3], &L_f = P[4], &L_l = P[5], &L_p = P[6], &C_ad = P[7], &C_d = P[8], &C_p = P[9];
  
  mc::FFVar RHS[NX];  // Right-hand side function
  RHS[0] = ( tau*(K*tanh(Lambda*P_d)-P_ad) )/(R_th*C_ad) + ( tau*U0*(U_f + U_p) )/(P0*C_ad);
  RHS[1] = ( tau*U0*U_f )/(P0*C_d);
  RHS[2] = ( tau*U0*U_p )/(P0*C_p);
  RHS[3] = ( (tau*P0/U0)*( L_l*P_p-L_p*P_ad-(L_p+L_l)*P_d ) - tau*(L_l*R_f+L_p*(R_f+R_l))*U_f - tau*L_p*R_l*U_p ) / ( L_l*(L_d+L_f)+L_p*(L_d+L_f+L_l) );
  RHS[4] = ( (tau*P0/U0)*( L_l*P_d - (L_d+L_f)*P_ad - (L_d+L_f+L_l)*P_p ) + tau*(L_l*R_f-(L_d+L_f)*R_l)*U_f - tau*(L_d+L_f)*R_l*U_p ) / ( L_l*(L_d+L_f)+L_p*(L_d+L_f+L_l) );
	  
  mc::FFVar Q[NQ];  // State quadratures
  for( unsigned i=0; i<NQ; i++ ) Q[i].set( &NIFTE );

  mc::FFVar QUAD[NQ];
  QUAD[0] = U_f * U0;	// to find V_f
  QUAD[1] = -(U_f + U_p)*U0;
  mc::FFVar dPaddt = ( tau*(K*((mc::exp(2*Lambda*P_d)-1.)/(mc::exp(2*Lambda*P_d)+1.))-P_ad) )/(R_th*C_ad) + ( tau*U0*(U_f + U_p) )/(P0*C_ad);
  mc::FFVar U_ad = C_ad *dPaddt *P0/tau;
  QUAD[2] = U_ad - U_f*U0 - U_p*U0;

  mc::FFVar &U_l = QUAD[1], &U_th = QUAD[2];

  mc::FFVar IC[NX];   // Initial value function
  IC[0] = 0.01;
  IC[1] = 0.;
  IC[2] = 0.01;
  IC[3] = 0.;
  IC[4] = 0.;

  /////////////////////////////////////////////////////////////////////////
  // Solve NIFTE dynamics

  mc::ODESLV_CVODES IVP;
  IVP.options.INTMETH   = mc::BASE_CVODES::Options::MSBDF;//MSADAMS;
  IVP.options.NLINSOL   = mc::BASE_CVODES::Options::NEWTON;//FIXEDPOINT;
  IVP.options.LINSOL    = mc::BASE_CVODES::Options::DENSE;
  IVP.options.DISPLAY   = 1;
  IVP.options.NMAX      = 10000;
  IVP.options.ATOL      = IVP.options.ATOLB      = IVP.options.ATOLS  = 1e-9;
  IVP.options.RTOL      = IVP.options.RTOLB      = IVP.options.RTOLS  = 1e-9;

  IVP.set_dag( &NIFTE );
  IVP.set_time( t0, tf );
  IVP.set_state( NX, X );
  IVP.set_parameter( NP, P );
  IVP.set_differential( NX, RHS );
  IVP.set_initial( NX, IC );
  IVP.set_quadrature( NQ, QUAD, Q );
  IVP.setup();
  
#if defined( SAVE_RESULTS )
  IVP.options.RESRECORD = 500;
#endif

  std::cout << "\nCONTINUOUS-TIME INTEGRATION:\n\n";
  IVP.states( p );//, xk );

#if defined( SAVE_RESULTS )
  std::ofstream direcSTA;
  direcSTA.open( "test3_STA.dat", std::ios_base::out );
  IVP.record( direcSTA );
#endif

  return 0;
}


