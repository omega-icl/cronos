const unsigned NSAMP = 100;    // <- Number of sampling points for inner approx.
#define SAVE_RESULTS	       // <- Whether to save bounds to file
#define USE_CMODEL		       // <- whether to use Chebyshev models or Taylor models
#define USE_PROFIL	           // specify to use PROFIL for interval arithmetic

#include "odebnd_expand.hpp"

#ifdef USE_PROFIL
  #include "mcprofil.hpp"
  typedef INTERVAL I;
#else
  #include "interval.hpp"
  typedef mc::Interval I;
#endif

#include "cmodel.hpp"
typedef mc::CModel<I> PM;
typedef mc::CVar<I> PV;

////////////////////////////////////////////////////////////////////////
// CYCLIC STEADY_STATE IN A NIFTE ENGINE
////////////////////////////////////////////////////////////////////////

// constant parameters  
const double P0     = 1.013e5;   // Pa
const double U0     = 8.00e-4;   // m3/s
const double tau    = 5.;        // s
const double Lambda = 330.;
const double K      = 1.37;      // control gain
const double R_l    = 4.0809e6;
const double R_f    = 2.1307e6;
const double R_th   = 5.0212e8;
const double L_d    = 1.8027e5;
const double L_l    = 1.2710e7;
const double L_p    = 3.7739e5;
const double L_f    = 4.7428e6;
const double C_ad   = 1.7618e-9;
const double C_d    = 7.3513e-8;
const double C_p    = 7.4280e-8;
const double R_TL1  = 8.8203e+06;
const double C_TL1  = 8.6087e-10;
const double R_TL2  = 1.7178e+07;
const double C_TL2  = 5.9528e-09;
const double R_TL3  = 3.8277e+07;
const double C_TL3  = 1.7939e-08;
const double R_TL4  = 2.8117e+08;
const double C_TL4  = 2.6578e-08;
const double R_TL5  = 3.1500e+07;

const double TF_REF    =  5.840515155e-01;
const double P_ad_REF  = -4.348489106e-02;
const double P_d_REF   = -4.690881398e-02;
const double P_p_REF   =  3.772299664e-02;
const double U_f_REF   =  0.000000000e+00;
const double U_p_REF   =  3.806980015e-01;
const double P_TL1_REF =  8.826363938e-04;
const double P_TL2_REF = -3.099606866e-03;
const double P_TL3_REF = -2.492678902e-02;
const double P_TL4_REF = -1.979265621e-02;

int main()
{
	
  mc::FFGraph NIFTE;

  // States
  const unsigned NX = 9;
  mc::FFVar X[NX];
  for( unsigned int i=0; i<NX; i++ ) X[i].set( &NIFTE );
  mc::FFVar &P_ad  = X[0],
            &P_d   = X[1], 
            &P_p   = X[2],
            &U_f   = X[3],
            &U_p   = X[4],
            &P_TL1 = X[5],
            &P_TL2 = X[6],
            &P_TL3 = X[7], 
            &P_TL4 = X[8];

  // Parameters
  const unsigned NP = NX;
  mc::FFVar P[NP];  for( unsigned int i=0; i<NP; i++ ) P[i].set( &NIFTE );
  mc::FFVar  TF     = TF_REF,
            &P_ad0  = P[0],
            &P_d0   = P[1],
            &P_p0   = P[2],
            &U_f0   = P[3],
            &U_p0   = P[4], 
            &P_TL10 = P[5],
            &P_TL20 = P[6],
            &P_TL30 = P[7],
            &P_TL40 = P[8];

  // Dynamics
  mc::FFVar SWITCH, RHS[NX];
  SWITCH = max( min( Lambda*P_d, 1. ), -1. );
  RHS[0] = TF * ( ( tau*(K*SWITCH-P_ad) )/(R_th*C_ad) 
                + ( tau*U0*(U_f + U_p) )/(P0*C_ad)
                - ( tau*(P_ad-P_TL1-P_TL2-P_TL3-P_TL4) )/(R_TL5*C_ad) );
  RHS[1] = TF * ( tau*U0*U_f ) / ( P0*C_d );
  RHS[2] = TF * ( tau*U0*U_p ) / ( P0*C_p );
  RHS[3] = TF * ( (tau*P0/U0)*( L_l*P_p-L_p*P_ad-(L_p+L_l)*P_d )
                 - tau*(L_l*R_f+L_p*(R_f+R_l))*U_f - tau*L_p*R_l*U_p )
                / ( L_l*(L_d+L_f)+L_p*(L_d+L_f+L_l) );
  RHS[4] = TF * ( (tau*P0/U0)*( L_l*P_d - (L_d+L_f)*P_ad - (L_d+L_f+L_l)*P_p )
                 + tau*(L_l*R_f-(L_d+L_f)*R_l)*U_f - tau*(L_d+L_f)*R_l*U_p )
                / ( L_l*(L_d+L_f)+L_p*(L_d+L_f+L_l) );
  RHS[5] = TF * ( (tau*(P_ad-P_TL1-P_TL2-P_TL3-P_TL4))/(R_TL5*C_TL1)
                - (tau*P_TL1)/(R_TL1*C_TL1) );
  RHS[6] = TF * ( (tau*(P_ad-P_TL1-P_TL2-P_TL3-P_TL4))/(R_TL5*C_TL2)
                - (tau*P_TL2)/(R_TL2*C_TL2) );
  RHS[7] = TF * ( (tau*(P_ad-P_TL1-P_TL2-P_TL3-P_TL4))/(R_TL5*C_TL3)
                - (tau*P_TL3)/(R_TL3*C_TL3) );
  RHS[8] = TF * ( (tau*(P_ad-P_TL1-P_TL2-P_TL3-P_TL4))/(R_TL5*C_TL4)
                - (tau*P_TL4)/(R_TL4*C_TL4) );
	
  // Quadratures
  const unsigned NQ = 3+3;
  mc::FFVar Q[NQ];
  for( unsigned int i=0; i<NQ; i++ ) Q[i].set( &NIFTE );

  mc::FFVar U_l  = -(U_f+U_p)*U0,
            P_l  = U_l*R_l,
            U_ad = C_ad*RHS[0]*P0/(tau*TF),
            U_TL = (P_ad-P_TL1-P_TL2-P_TL3-P_TL4)*P0/R_TL5,
            U_th = U_ad + U_l + U_TL,
            P_th = R_th*U_th + P_ad*P0,
            P_f  = R_f*U_f*U0;

  mc::FFVar QUAD[NQ];
  QUAD[0] = P_l * U_l * tau;      // Average power dissipated in the load * tau
  QUAD[1] = P_th * U_th * tau;    // Average power input
  QUAD[2] = P_f * U_f * U0 * tau; // Average power dissipated in the feedback tube
  QUAD[3] = TF * tau * U_f * U0;  // scaled back, V_f
  QUAD[4] = TF * tau * U_l;       // scaled back, V_l
  QUAD[5] = TF * tau * U_th;      // scaled back, V_th

  // Initial conditions
  mc::FFVar IC[NX];
  IC[0] = P_ad0;
  IC[1] = P_d0;
  IC[2] = P_p0;
  IC[3] = U_f0;
  IC[4] = U_p0;
  IC[5] = P_TL10;
  IC[6] = P_TL20;
  IC[7] = P_TL30;
  IC[8] = P_TL40;

  // Parameter bounds
  I Ip[NP];
  Ip[0] = P_ad_REF  * I( 0.9, 1.1 );
  Ip[1] = P_d_REF   * I( 0.9, 1.1 );
  Ip[2] = P_p_REF   * I( 0.9, 1.1 );
  Ip[3] = U_f_REF   * I( 0.9, 1.1 );
  Ip[4] = U_p_REF   * I( 0.9, 1.1 );
  Ip[5] = P_TL1_REF * I( 0.9, 1.1 );
  Ip[6] = P_TL2_REF * I( 0.9, 1.1 );
  Ip[7] = P_TL3_REF * I( 0.9, 1.1 );
  Ip[8] = P_TL4_REF * I( 0.9, 1.1 );


  /////////////////////////////////////////////////////////////////////////
  // ODE Discretization

  mc::ODEBND_EXPAND<I,PM,PV> bnd;
  bnd.set_dag( &NIFTE );
  bnd.set_time( 0., 1. );
  bnd.set_parameter( NP, P );
  bnd.set_state( NX, X );
  bnd.set_differential( NX, RHS );
  //bnd.set_quadrature( NQ, QUAD, Q );
  bnd.set_initial( NX, IC );

#if defined( SAVE_RESULTS )
  bnd.options.RESRECORD = 1000;
#endif
  bnd.options.INTMETH   = mc::ODEBND_EXPAND<I,PM,PV>::Options::METHOD::RK;//TS
  bnd.options.TORD      = 4;
  bnd.options.H0        = 0.005;
  bnd.options.LBLK      = 10;
  bnd.options.DBLK      = 10;//50/NS;
  bnd.options.DISPLAY   = 1;
  bnd.options.RESRECORD = true;
  bnd.options.AEBND.DISPLAY = 2;
  bnd.options.AEBND.INTERBND = true;
  bnd.options.AEBND.BOUNDER = mc::AEBND<I,PM,PV>::Options::ALGORITHM::GE;//KRAW;//AUTO;
  bnd.options.AEBND.PRECOND = mc::AEBND<I,PM,PV>::Options::PRECONDITIONING::INVMB;//INVBD;//NONE;
  bnd.options.AEBND.BLKDEC  = mc::AEBND<I,PM,PV>::Options::DECOMPOSITION::RECUR;//DIAG;//NONE;
  //bnd.options.ODESLVS   = bnd.options;
  bnd.setup();
 
  std::ofstream oRES( "NIFTE_RES.dot", std::ios_base::out );
  //NIFTE.dot_script( (NX+NQ)*(bnd.options.LBLK+1), bnd.vRES().data(), oRES );
  //NIFTE.dot_script( (NX+NQ)*bnd.options.LBLK, bnd.vRES().data()+NX+NQ, oRES );
  NIFTE.dot_script( (NX)*bnd.options.LBLK, bnd.vRES().data()+NX, oRES );
  oRES.close();
  //NIFTE.output( NIFTE.subgraph( (NX+NQ)*bnd.options.LBLK, bnd.vRES().data()+NX+NQ ) ); 
  NIFTE.output( NIFTE.subgraph( (NX)*bnd.options.LBLK, bnd.vRES().data()+NX ) ); 
 
//  std::cout << "\nNON_VALIDATED INTEGRATION - INNER-APPROXIMATION OF REACHABLE SET:\n\n";
//  bnd.bounds( -NSAMP, Ip );
//#if defined( SAVE_RESULTS )
//  std::ofstream apprec( "test8_APPROX_STA.dat", std::ios_base::out );
//  bnd.record( apprec );
//#endif

  return 0;
}
