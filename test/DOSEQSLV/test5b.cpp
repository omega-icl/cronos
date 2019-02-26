#define SAVE_TRAJECTORIES	// <- Whether to save trajectories to file

#include <fstream>
#include <iomanip>
#include "doseqslv_ipopt.hpp"
#include "interval.hpp"
#include "mclapack.hpp"

////////////////////////////////////////////////////////////////////////
// TEST PROBLEM #5: CYCLIC STEADY_STATE DETECTION IN A NIFTE ENGINE
// http://www.sciencedirect.com/science/article/pii/0898122187900459
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

// optimization criterion
const unsigned crit = 1;

////////////////////////////////////////////////////////////////////////
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

  // Decision variables
  const unsigned NP = NX+1;
  mc::FFVar P[NP];
  for( unsigned int i=0; i<NP; i++ ) P[i].set( &NIFTE );
  mc::FFVar &TF     = P[0],
            &P_ad0  = P[1],
            &P_d0   = P[2],
            &P_p0   = P[3],
            &U_f0   = P[4],
            &U_p0   = P[5], 
            &P_TL10 = P[6],
            &P_TL20 = P[7],
            &P_TL30 = P[8],
            &P_TL40 = P[9];

  // Dynamics
  mc::FFVar RHS[NX];
  RHS[0] = TF * ( ( tau*(K*tanh(Lambda*P_d)-P_ad) )/(R_th*C_ad) 
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

  // Dynamic optimization problem
  Ipopt::SmartPtr<mc::DOSEQSLV_IPOPT> CSS = new mc::DOSEQSLV_IPOPT;
  CSS->set_dag( &NIFTE );
  CSS->set_time( 0., 1. );
  CSS->set_parameter( NP, P );
  CSS->set_state( NX, X );
  CSS->set_differential( NX, RHS );
  CSS->set_quadrature( NQ, QUAD, Q );
  CSS->set_initial( NX, IC );

  switch( crit ){
    case 1:  
      CSS->set_obj( mc::DOSEQSLV_IPOPT::MAX, Q[0]/Q[1] );
      CSS->add_ctr( mc::DOSEQSLV_IPOPT::GE, std::make_pair( 0, Q[0]/tau ) ); break;
    case 2:
      CSS->set_obj( mc::DOSEQSLV_IPOPT::MAX, Q[0]/tau );
      CSS->add_ctr( mc::DOSEQSLV_IPOPT::GE, std::make_pair( 0, Q[0]/Q[1] ) ); break;
    default:
      CSS->set_obj( mc::DOSEQSLV_IPOPT::MIN, TF );
      CSS->add_ctr( mc::DOSEQSLV_IPOPT::GE, std::make_pair( 0, Q[0]/tau ) );
      CSS->add_ctr( mc::DOSEQSLV_IPOPT::GE, std::make_pair( 0, Q[0]/Q[1] ) ); break;
  }

  CSS->add_ctr( mc::DOSEQSLV_IPOPT::EQ, std::make_pair( 0, P_ad  - P_ad0 ) );
  CSS->add_ctr( mc::DOSEQSLV_IPOPT::EQ, std::make_pair( 0, P_d   - P_d0 ) );
  CSS->add_ctr( mc::DOSEQSLV_IPOPT::EQ, std::make_pair( 0, P_p   - P_p0 ) );
  CSS->add_ctr( mc::DOSEQSLV_IPOPT::EQ, std::make_pair( 0, U_f   - U_f0 ) );
  CSS->add_ctr( mc::DOSEQSLV_IPOPT::EQ, std::make_pair( 0, U_p   - U_p0 ) );
  CSS->add_ctr( mc::DOSEQSLV_IPOPT::EQ, std::make_pair( 0, P_TL1 - P_TL10 ) );
  CSS->add_ctr( mc::DOSEQSLV_IPOPT::EQ, std::make_pair( 0, P_TL2 - P_TL20 ) );
  CSS->add_ctr( mc::DOSEQSLV_IPOPT::EQ, std::make_pair( 0, P_TL3 - P_TL30 ) );
  CSS->add_ctr( mc::DOSEQSLV_IPOPT::EQ, std::make_pair( 0, P_TL4 - P_TL40 ) );
  CSS->add_ctr( mc::DOSEQSLV_IPOPT::GE, std::make_pair( 0, Q[1]  - 1e0 ) );

  CSS->options.DISPLAY   = 0;
  CSS->options.MAXITER   = 50;
  CSS->options.TESTDER   = false;
  CSS->options.CVTOL     = 1e-9;
  CSS->options.GRADIENT  = mc::DOSEQSLV_IPOPT::Options::FORWARD; //BACKWARD;
  CSS->options.ODESLVS.INTMETH   = mc::BASE_SUNDIALS::Options::MSBDF; // MSADAMS;
  CSS->options.ODESLVS.JACAPPROX = mc::BASE_SUNDIALS::Options::CV_DENSE; //CV_DIAG;
  CSS->options.ODESLVS.NMAX = CSS->options.ODESLVS.ASACHKPT = 10000;
  CSS->options.ODESLVS.DISPLAY = 0;
  CSS->options.ODESLVS.RESRECORD = 0;
  CSS->options.ODESLVS.ATOL = CSS->options.ODESLVS.ATOLB = CSS->options.ODESLVS.ATOLS = 1e-8;
  CSS->options.ODESLVS.RTOL = CSS->options.ODESLVS.RTOLB = CSS->options.ODESLVS.RTOLS = 1e-7;
  CSS->setup();

  // Parameter bounds
  typedef mc::Interval I;
  I Ip[NP] = { I( 0.1, 2.0 ) };
  for( unsigned i=1; i<NP; i++ ) Ip[i] = ( i==4? I( 0. ): I( -1., 1. ) );
  //double p0[NP] = { 0.5 };
  //for( unsigned i=1; i<NX; i++ ) p0[i] = ( i==4? 0.: 0.05 );
	
  // Numerical CSS solution
  //int stat = CSS->solve( Ip, p0 );
  int stat = CSS->solve( 10, Ip );

  const unsigned IPREC = 9;
  std::cout << std::endl
            << "STATUS: " << std::setw(4) << stat << std::endl
            << std::scientific << std::setprecision(IPREC);
  for( unsigned i=0; i<CSS->solution().p.size(); i++ )
    std::cout << std::setw(IPREC+9) << CSS->solution().p[i];
  std::cout << std::setw(IPREC+9) << CSS->solution().f << std::endl;
  for( unsigned i=0; i<CSS->solution().g.size(); i++ )
    std::cout << std::setw(IPREC+9) << CSS->solution().g[i] << std::endl;
  std::cout << std::endl;
 
#ifdef SAVE_TRAJECTORIES
  std::ofstream otraj( "test5b.traj", std::ios_base::out );
  CSS->options.ODESLVS.DISPLAY   = 0;
  CSS->options.ODESLVS.RESRECORD = 1000;	// number of sample points in the .dat file
  CSS->set_time( 0., 1. );
  CSS->setup();
  CSS->states( CSS->solution().p.data() );
  CSS->record( otraj );
  otraj.close();
  //return 0;
#endif

  // CSS analysis
  double *xk[2] = { 0, 0 }, *xpk[2] = { 0, 0 };
  CSS->options.ODESLVS.DISPLAY = 0;
  CSS->states_FSA( CSS->solution().p.data(), xk, 0, xpk );
  CPPL::dgematrix JACF( NX, NX );
  for( unsigned i=0; i<NX; i++ )
    for( unsigned j=0; j<NX; j++ )
      JACF(i,j) = xpk[1][(j+1)*(NX+NQ)+i];
  std::cout << "\nMonodromy Matrix:\n" << JACF;
  std::vector<double> wr, wi;
  JACF.dgeev( wr, wi );
  std::cout << "\nEigenvalues of Monodromy Matrix:\n";
  for( unsigned i=0; i<NX; i++ )
    std::cout <<  "(" << wr[i] << ", " << wi[i] << ")\n";
  for( unsigned k=0; k<2; k++ ){
    delete[] xk[k];
    delete[] xpk[k];
  }

  return 0;
}

