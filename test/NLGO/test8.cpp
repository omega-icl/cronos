#define SAVE_RESULTS	       // <- Whether to save bounds to file
#define USE_CMODEL		       // <- whether to use Chebyshev models or Taylor models
#undef  USE_PROFIL	           // specify to use PROFIL for interval arithmetic
#define  MC__FFUNC_DEBUG_EVAL

#include "odebnd_expand.hpp"
#include "nlgo.hpp"

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
//const double L_p    = 3.7739e5 * 39.81071706;
const double L_p    = 3.7739e5;
const double L_f    = 4.7428e6;
const double C_ad   = 1.7618e-9;
const double C_d    = 7.3513e-8;
//const double C_p    = 7.4280e-8 * 0.1202264435;
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

//////// css1 - stable
//const double TF_REF    =  5.103390707e-01;
//const double P_ad_REF  =  9.365730873e-02;
//const double P_d_REF   =  4.776426344e-02;
//const double P_p_REF   = -3.088034797e-01;
//const double U_f_REF   =  0.000000000e+00;
//const double U_p_REF   = -2.962018115e-01;
//const double P_TL1_REF =  5.166030950e-03;
//const double P_TL2_REF =  1.728732917e-02;
//const double P_TL3_REF =  3.287392167e-02;
//const double P_TL4_REF =  2.067277687e-02;

//////// css2 - unstable
//const double TF_REF    =  2.622451394e-01;
//const double P_ad_REF  =  1.408898506e-01;
//const double P_d_REF   =  3.813808745e-03;
//const double P_p_REF   = 	1.339098843e-02;
//const double U_f_REF   =  0.000000000e+00;
//const double U_p_REF   =  7.662042110e-03;
//const double P_TL1_REF =  1.873265679e-02;
//const double P_TL2_REF =  3.525918952e-02;
//const double P_TL3_REF =  1.482245760e-02;
//const double P_TL4_REF =  6.166426847e-03;

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
  const unsigned NP = NX-1;
  mc::FFVar P[NP];
  for( unsigned int i=0; i<NP; i++ ) P[i].set( &NIFTE );
  mc::FFVar  TF     = TF_REF,
            &P_ad0  = P[0],
            &P_d0   = P[1],
            &P_p0   = P[2],
             U_f0   = U_f_REF,
            &U_p0   = P[3], 
            &P_TL10 = P[4],
            &P_TL20 = P[5],
            &P_TL30 = P[6],
            &P_TL40 = P[7];

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
  I Ip[NP], UNC( 0.5, 2.0 ); //UNC( 0.9, 1.1 ) //UNC( 0.99, 1.01 )
  Ip[0] = P_ad_REF  * UNC;
  Ip[1] = P_d_REF   * UNC;
  Ip[2] = P_p_REF   * UNC;
  Ip[3] = U_p_REF   * UNC;
  Ip[4] = P_TL1_REF * UNC;
  Ip[5] = P_TL2_REF * UNC;
  Ip[6] = P_TL3_REF * UNC;
  Ip[7] = P_TL4_REF * UNC;

  // Initial state bounds
  I Ix0[NX];
  NIFTE.eval( NX, IC, Ix0, NP, P, Ip );

  // A priori state bounds
  I Ix[NX];
  for( unsigned i=0; i<NX; i++ )
    Ix[i] = I( -2., 2. );
    
  /////////////////////////////////////////////////////////////////////////
  // ODE Discretization

  mc::ODEBND_EXPAND<I,PM,PV> ODEDISCR;
  ODEDISCR.set_dag( &NIFTE );
  ODEDISCR.set_time( 0., 1. );
  ODEDISCR.set_parameter( NP, P );
  ODEDISCR.set_state( NX, X );
  ODEDISCR.set_differential( NX, RHS );
  //ODEDISCR.set_quadrature( NQ, QUAD, Q );
  ODEDISCR.set_initial( NX, IC );

#if defined( SAVE_RESULTS )
  ODEDISCR.options.RESRECORD = 1000;
#endif
  ODEDISCR.options.INTMETH   = mc::ODEBND_EXPAND<I,PM,PV>::Options::METHOD::RK;//TS
  ODEDISCR.options.TORD      = 4; // 2;
  ODEDISCR.options.H0        = 0.005; //0.002;
  ODEDISCR.options.LBLK      = 200; //200;
  ODEDISCR.options.DBLK      = 200; //200;//50/NS;
  ODEDISCR.options.DISPLAY   = 1;
  ODEDISCR.options.RESRECORD = true;
  ODEDISCR.options.AEBND.DISPLAY = 2;
  ODEDISCR.options.AEBND.INTERBND = true;
  ODEDISCR.options.AEBND.BOUNDER = mc::AEBND<I,PM,PV>::Options::ALGORITHM::GE;//KRAW;//AUTO;
  ODEDISCR.options.AEBND.PRECOND = mc::AEBND<I,PM,PV>::Options::PRECONDITIONING::INVMB;//INVBD;//NONE;
  ODEDISCR.options.AEBND.BLKDEC  = mc::AEBND<I,PM,PV>::Options::DECOMPOSITION::RECUR;//DIAG;//NONE;

  ODEDISCR.setup();
  /*
  std::ofstream oRES( "NIFTE_RES.dot", std::ios_base::out );
  //NIFTE.dot_script( (NX+NQ)*(ODEDISCR.options.LBLK+1), ODEDISCR.vRES().data(), oRES );
  //NIFTE.dot_script( (NX+NQ)*ODEDISCR.options.LBLK, ODEDISCR.vRES().data()+NX+NQ, oRES );
  NIFTE.dot_script( (NX)*ODEDISCR.options.LBLK, ODEDISCR.vRES().data()+NX, oRES );
  oRES.close();
  //NIFTE.output( NIFTE.subgraph( (NX+NQ)*ODEDISCR.options.LBLK, ODEDISCR.vRES().data()+NX+NQ ) ); 
  NIFTE.output( NIFTE.subgraph( (NX)*ODEDISCR.options.LBLK, ODEDISCR.vRES().data()+NX ) ); 
 */

  /////////////////////////////////////////////////////////////////////////
  // MILP Representation
  /*
  mc::PolImg<I> PolEnv;
  PolEnv.options.AGGREG_LIN = true;
  PolEnv.options.BREAKPOINT_DISC = true;

  ODEDISCR.vVAR()[NP+NX+1].set( ODEDISCR.options.H0 );
  std::vector<mc::PolVar<I>> PVp, PVx, PVf(ODEDISCR.vRES().size());
  std::vector<I> Ipall, Ixall, If(ODEDISCR.vRES().size());
  for( unsigned i=0; i<NP; i++ ){
    PVp.push_back( mc::PolVar<I>( &PolEnv, ODEDISCR.vVAR()[i], Ip[i] ) );
    Ipall.push_back( Ip[i] );
  }
//  for( unsigned i=0; i<NX+1; i++ ){
//    PVp.push_back( mc::PolVar<I>(0.) );
//    Ipall.push_back( I(0.) );
//  }
//  PVp.push_back( mc::PolVar<I>(&PolEnv, ODEDISCR.vVAR()[NP+NX+1], ODEDISCR.options.H0) );
//  Ipall.push_back( I(ODEDISCR.options.H0) );
  
  for( unsigned k=0; k<=ODEDISCR.options.LBLK; k++ )
    for( unsigned i=0; i<NX; i++ ){
      PVx.push_back( mc::PolVar<I>( &PolEnv, ODEDISCR.vDEP()[NX*k+i], Ix[i] ) );
      Ixall.push_back( Ix[i] );
    }
  std::cout << PVx.size() << " =?= " << ODEDISCR.vDEP().size() << std::endl;

  NIFTE.output( NIFTE.subgraph( 1, ODEDISCR.vRES().data()+NX ) ); 
  for( unsigned k=0; k<PVp.size(); k++ )
    std::cout << ODEDISCR.vVAR()[k] << " in " << Ipall[k] << std::endl;
  for( unsigned k=0; k<PVx.size(); k++ )
    std::cout << ODEDISCR.vDEP()[k] << " in " << Ixall[k] << std::endl;
  
//  NIFTE.eval( 1, ODEDISCR.vRES().data()+NX, If.data()+NX, 
//              PVp.size(), ODEDISCR.vVAR().data(), Ipall.data(), 
//              PVx.size(), ODEDISCR.vDEP().data(), Ixall.data() );
//  std::cout << "If = " << If[NX] << std::endl;

//  NIFTE.eval( 1, ODEDISCR.vRES().data()+NX, PVf.data()+NX, 
  NIFTE.eval( ODEDISCR.vDEP().size(), ODEDISCR.vRES().data(), PVf.data(),   
              PVp.size(), ODEDISCR.vVAR().data(), PVp.data(), 
              PVx.size(), ODEDISCR.vDEP().data(), PVx.data() );
  std::cout << PolEnv;
//  PolEnv.generate_cuts( 1, PVf.data()+NX, true );
  PolEnv.generate_cuts( ODEDISCR.vDEP().size(), PVf.data(), true );
  std::cout << PolEnv;

  return 0;
  */

#if defined( SAVE_RESULTS )
  std::ofstream bndSTA;
  bndSTA.open( "test8_STA3.dat", std::ios_base::app );
#endif

  ODEDISCR.vVAR()[NP+NX+1].set( ODEDISCR.options.H0 );
  std::vector<I> IALL( Ip, Ip+NP );
  IALL.insert( IALL.end(), Ix0, Ix0+NX );
  for( unsigned k=0; k<=ODEDISCR.options.LBLK; k++ )
    IALL.insert( IALL.end(), Ix, Ix+NX );

  int status( 0 );    
  for( unsigned NT=1; NT<=ODEDISCR.options.LBLK; NT+=5 ){
    const double tf = TF_REF * NT * ODEDISCR.options.H0;
#if defined( SAVE_RESULTS )
    bndSTA << tf;
#endif

    mc::NLGO<I> NLP;  
    NLP.options.MIPFILE         = "";//test8.lp";
    NLP.options.RELMETH         = mc::NLGO<I>::Options::DRL;
    NLP.options.MAXITER         = 10;
    NLP.options.DEPSUSE         = false;//true; //false;
    NLP.options.CMODDEPS        = 0;
    NLP.options.NCOCUTS         = false;//true; //false;
    NLP.options.AEBND.DISPLAY   = 0;
    NLP.options.PREPROC         = true;
    NLP.options.DISPLAY         = 2;
    NLP.options.MIPDISPLAY      = 1;
    NLP.options.NLPSLV.DISPLAY  = 0;
    NLP.options.POLIMG.AGGREG_LIN      = true;
    NLP.options.POLIMG.BREAKPOINT_DISC = true;
    //NLP.options.POLIMG.SANDWICH_MAXCUT = 4;
    //NLP.options.POLIMG.BREAKPOINT_TYPE = mc::PolImg<I>::Options::SOS2;//NONE;
    //NLP.options.POLIMG.DCDECOMP_SCALE = false;

    NLP.set_dag( &NIFTE );
    NLP.set_var( NP+NX, ODEDISCR.vVAR().data() );
    NLP.set_dep( NX*(NT+1), ODEDISCR.vDEP().data(), ODEDISCR.vRES().data() );
 
    // Minimize the value of P_ad at final time
    std::cout << "setting-up" << std::endl;
    NLP.set_obj( mc::BASE_NLP::MIN, ODEDISCR.vDEP()[NX*NT] );   
    NLP.setup( IALL.data() ); 
    std::cout << " OK!" << std::endl;
    status = NLP.relax( IALL.data() );
    const double P_ad_L = NLP.get_objective();
    std::cout << "relaxed objective: " << P_ad_L << std::endl;
#if defined( SAVE_RESULTS )
    bndSTA << "  " << P_ad_L;
#endif

    // Maximize the value of P_ad at final time
    std::cout << "setting-up" << std::endl;
    NLP.set_obj( mc::BASE_NLP::MAX, ODEDISCR.vDEP()[NX*NT] );   
    NLP.setup( IALL.data() ); 
    std::cout << " OK!" << std::endl;
    status = NLP.relax( IALL.data() );
    const double P_ad_U = NLP.get_objective();
    std::cout << "relaxed objective: " << P_ad_U << std::endl;
#if defined( SAVE_RESULTS )
    bndSTA << "  " << P_ad_U << std::endl;
#endif
  }

#if defined( SAVE_RESULTS )
    bndSTA.close();
#endif
  return status;
}
