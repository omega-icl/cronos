// with WARM-START
// Varying L_f only
// http://www.sciencedirect.com/science/article/pii/0898122187900459

#define SAVE_RESULTS		// <- Whether to save results to file
#define SAVE_TRAJECTORIES	// <- Whether to save trajectories to file

#include <fstream>
#include <iomanip>
#include "doseqslv_ipopt.hpp"
#include "interval.hpp"
#include "mclapack.hpp"

using std::cout;
using std::endl;
typedef mc::Interval I;

////////////////////////////////////////////////////////////////////////
int main()
{
#ifdef SAVE_RESULTS
  std::ofstream ores( "test5b_1.dat", std::ios_base::out );
#endif

  mc::FFGraph NIFTE;

  // constant parameters  
  const double P0 = 1.013e5;		// Pa
  const double U0 = 8.00e-4;		// m3/s
  const double tau = 5.;				// s
  const double Lambda = 330.;
  const double K = 1.37;				// control gain
  
  const unsigned NDOF = 10;    // Number of degrees of freedom
  const double R_l_0  = 4.0809e6;
  const double R_f_0  = 2.1307e6;
  const double R_th_0 = 5.0212e8;
  const double L_d_0  = 1.8027e5;
  const double L_l_0  = 1.2710e7;
  const double L_p_0  = 3.7739e5;
  const double L_f_0  = 4.7428e6;
  const double C_ad_0 = 1.7618e-9;
  const double C_d_0  = 7.3513e-8;
  const double C_p_0  = 7.4280e-8;
	
  const unsigned NX = 5;       // Number of states
  mc::FFVar X[NX];             // States
  for( unsigned int i=0; i<NX; i++ ) X[i].set( &NIFTE );

  const unsigned NP = 1+NX+NDOF;    // Number of parameters
  mc::FFVar P[NP];             // Parameters
  for( unsigned int i=0; i<NP; i++ ) P[i].set( &NIFTE );
  mc::FFVar &TF  = P[0],  &Pad_0 = P[1],  &Pd_0 = P[2],  &Pp_0 = P[3],  &Uf_0 = P[4],  &Up_0 = P[5],
            &R_l = P[6],  &R_f   = P[7],  &R_th = P[8],  &L_d  = P[9],  &L_l  = P[10], &L_p  = P[11],
            &L_f = P[12], &C_ad  = P[13], &C_d  = P[14], &C_p  = P[15];

  mc::FFVar &P_ad = X[0], &P_d = X[1], &P_p = X[2], &U_f = X[3], &U_p = X[4]; 
  mc::FFVar RHS[NX];        // Right-hand side function
  RHS[0] = TF * ( ( tau*(K*tanh(Lambda*P_d)-P_ad) )/(R_th*C_ad) + ( tau*U0*(U_f + U_p) )/(P0*C_ad) );
  RHS[1] = TF * ( tau*U0*U_f ) / ( P0*C_d );
  RHS[2] = TF * ( tau*U0*U_p ) / ( P0*C_p );
  RHS[3] = TF * ( (tau*P0/U0)*( L_l*P_p-L_p*P_ad-(L_p+L_l)*P_d ) - tau*(L_l*R_f+L_p*(R_f+R_l))*U_f - tau*L_p*R_l*U_p )
              / ( L_l*(L_d+L_f)+L_p*(L_d+L_f+L_l) );
  RHS[4] = TF * ( (tau*P0/U0)*( L_l*P_d - (L_d+L_f)*P_ad - (L_d+L_f+L_l)*P_p ) + tau*(L_l*R_f-(L_d+L_f)*R_l)*U_f - tau*(L_d+L_f)*R_l*U_p )
              / ( L_l*(L_d+L_f)+L_p*(L_d+L_f+L_l) );	

  const unsigned NQ = 3;
  mc::FFVar Q[NQ];
  for( unsigned int i=0; i<NQ; i++ ) Q[i].set( &NIFTE );
  mc::FFVar U_l = -(U_f + U_p)*U0;
  mc::FFVar P_l = U_l * R_l;
  mc::FFVar U_ad = C_ad * RHS[0]*P0/(tau*TF);
  mc::FFVar U_th = U_ad + U_l;
  mc::FFVar P_th = R_th * U_th + P_ad*P0;
  mc::FFVar P_f = R_f * U_f * U0;
  mc::FFVar QUAD[NQ];  
  QUAD[0] = P_l * U_l * tau;		// Average power dissipated in the load
  QUAD[1] = P_th * U_th * tau;	// Average power input
  QUAD[2] = P_f * U_f *U0 * tau;		// Average power dissipated in the feedback tube
	
  mc::FFVar IC[NX];
  IC[0] = Pad_0;
  IC[1] = Pd_0;
  IC[2] = Pp_0;
  IC[3] = Uf_0;
  IC[4] = Up_0;
		
  Ipopt::SmartPtr<mc::DOSEQSLV_IPOPT> OC = new mc::DOSEQSLV_IPOPT;
  OC->set_dag( &NIFTE );
  OC->set_time( 0., 1. );
  OC->set_parameter( NP, P );
  OC->set_state( NX, X );
  OC->set_differential( NX, RHS );
  OC->set_quadrature( NQ, QUAD, Q );
  OC->set_initial( NX, IC );

  OC->set_obj( mc::DOSEQSLV_IPOPT::MAX, Q[0]/Q[1] );  // max efficiency
  //OC->set_obj( mc::DOSEQSLV_IPOPT::MAX, std::make_pair( 0, Q[0]/Q[1] + mc::sqr(P_ad-Pad_0) + mc::sqr(P_d-Pd_0) + mc::sqr(P_p-Pp_0) + mc::sqr(U_f-Uf_0) + mc::sqr(U_p-Up_0) ) );  // max efficiency
  OC->add_ctr( mc::DOSEQSLV_IPOPT::EQ, std::make_pair( 0, P_ad - Pad_0 ) );
  OC->add_ctr( mc::DOSEQSLV_IPOPT::EQ, std::make_pair( 0, P_d - Pd_0 ) );
  OC->add_ctr( mc::DOSEQSLV_IPOPT::EQ, std::make_pair( 0, P_p - Pp_0 ) );
  OC->add_ctr( mc::DOSEQSLV_IPOPT::EQ, std::make_pair( 0, U_f - Uf_0 ) );
  OC->add_ctr( mc::DOSEQSLV_IPOPT::EQ, std::make_pair( 0, U_p - Up_0 ) );
  OC->add_ctr( mc::DOSEQSLV_IPOPT::GE, std::make_pair( 0, Q[1] - 1e0 ) );

  OC->options.DISPLAY   = 5;
  OC->options.MAXITER   = 50;
  OC->options.TESTDER   = false;
  OC->options.CVTOL     = 1e-9;
  OC->options.GRADIENT  = mc::DOSEQSLV_IPOPT::Options::FORWARD; //FORWARD; //BACKWARD;
  OC->options.ODESLVS.INTMETH   = mc::BASE_SUNDIALS::Options::MSBDF; // MSADAMS;
  OC->options.ODESLVS.JACAPPROX = mc::BASE_SUNDIALS::Options::CV_DENSE; //CV_DIAG;
  OC->options.ODESLVS.NMAX = 10000;
  OC->options.ODESLVS.DISPLAY = 0;
  OC->options.ODESLVS.RESRECORD = 0;
  OC->options.ODESLVS.ATOL = OC->options.ODESLVS.ATOLB = OC->options.ODESLVS.ATOLS = 1e-9;
  OC->options.ODESLVS.RTOL = OC->options.ODESLVS.RTOLB = OC->options.ODESLVS.RTOLS = 1e-9;
  OC->setup();

  // Parameter bounds and values
  I Ip[NP] = { I( 0.1, 5.0 ), I( -0.8, 0.8 ), I( -0.2, 0.2 ), I( -0.2, 0.2 ), 0., I( -2.5, 2.5 ),
                    4.080900000e+06, 2.130700000e+06, 5.021200000e+08, 1.802700000e+05, 1.271000000e+07, 3.773900000e+05,
                    4.742800000e+05, 1.761800000e-09, 7.351300000e-08, 7.428000000e-08 }; //, R_l_0, R_f_0, R_th_0, L_d_0, L_l_0, L_p_0, L_f_0, C_ad_0, C_d_0, C_p_0 };
  //double p0[NP] = { 0.5, -0.3, 0.05, -0.06, 0., -0.4, R_l_0, R_f_0, R_th_0, L_d_0, L_l_0, L_p_0, L_f_0, C_ad_0, C_d_0, C_p_0 };

  //double p0[NP] = { 0.5, 1e-4, 1e-4, 1e-4, 0., 1e-4,
  //                  4.080900000e+06, 2.130700000e+06, 5.021200000e+08, 1.802700000e+05, 1.271000000e+07, 3.773900000e+05,
  //                  4.742800000e+05, 1.761800000e-09, 7.351300000e-08, 7.428000000e-08 };      

  // CSS #1 for L_f = 4.7428e5
  double p0[NP] = { 2.617581565e-01, 1.449924893e-01, 7.698273696e-03, -2.036559682e-03, 0.000000000e+00, -7.559298732e-01,
                    4.080900000e+06, 2.130700000e+06, 5.021200000e+08, 1.802700000e+05, 1.271000000e+07, 3.773900000e+05,
                    4.742800000e+05, 1.761800000e-09, 7.351300000e-08, 7.428000000e-08 }; //  3.945433878e+00   1.285006632e+02   4.051040089e-01   3.070360713e-02  -4.904965323e-13  -3.988996633e-15  -1.628861976e-14  -2.797810798e-13  -5.095923683e-14

  // perturbed x[0] in CSS #1
  //double p0[NP] = { 2.617581565e-01, 0.449924893e-01, 7.698273696e-03, -2.036559682e-03, 0.000000000e+00, -7.559298732e-01,
  //                  4.080900000e+06, 2.130700000e+06, 5.021200000e+08, 1.802700000e+05, 1.271000000e+07, 3.773900000e+05,
  //                  4.742800000e+05, 1.761800000e-09, 7.351300000e-08, 7.428000000e-08 }; //  3.945433878e+00   1.285006632e+02   4.051040089e-01   3.070360713e-02  -4.904965323e-13  -3.988996633e-15  -1.628861976e-14  -2.797810798e-13  -5.095923683e-14

  // CSS #2 for L_f = 4.7428e5
  //double p0[NP] = { 1.761471327e-01, 6.061150098e-01, 2.941732991e-03, 9.965898403e-03, 0.000000000e+00, 5.180310276e-01,
  //                  4.080900000e+06, 2.130700000e+06, 5.021200000e+08, 1.802700000e+05, 1.271000000e+07, 3.773900000e+05,
  //                  4.742800000e+05, 1.761800000e-09, 7.351300000e-08, 7.428000000e-08 }; //  6.639838830e+00   2.716861146e+01   1.301496733e-01   2.443937498e-01   3.718136909e-13   1.105452535e-15   6.550315845e-15   1.338842077e-13   1.203481759e-12

  // perturbed x[0] in CSS #2
  //double p0[NP] = { 1.761471327e-01, 6.001150098e-01, 2.941732991e-03, 9.965898403e-03, 0.000000000e+00, 5.180310276e-01,
  //                  4.080900000e+06, 2.130700000e+06, 5.021200000e+08, 1.802700000e+05, 1.271000000e+07, 3.773900000e+05,
  //                  4.742800000e+05, 1.761800000e-09, 7.351300000e-08, 7.428000000e-08 }; //  6.639838830e+00   2.716861146e+01   1.301496733e-01   2.443937498e-01   3.718136909e-13   1.105452535e-15   6.550315845e-15   1.338842077e-13   1.203481759e-12

  // CSS #1 for L_f = 4.7428e6
  //double p0[NP] = { 5.621991177708545e-01, -3.689528171983159e-01, 5.389414281510629e-02, -6.109921182890796e-02, 0.000000000000000e+00, -4.389440961841696e-01,
  //                  4.080900000e+06, 2.130700000e+06, 5.021200000e+08, 1.802700000e+05, 1.271000000e+07, 3.773900000e+05,
  //                  4.742800000e+06, 1.761800000e-09, 7.351300000e-08, 7.428000000e-08 };      

//  std::vector<double> dRHS(NX);
//  for( unsigned i=0; i<NX; i++ ) P[i+1] = X[i];
//  NIFTE.eval( 1, RHS, dRHS.data(), NP, P, p0 ); 
//  for( unsigned i=0; i<1; i++ ) std::cout << "RHS[" << i << "] = " << dRHS[i] << std::endl;

//  auto dRHS0dP = NIFTE.FAD( 1, RHS, NX+1, P );
//  std::vector<double> ddRHS0dP(NX+1);
//  NIFTE.eval( NX+1, dRHS0dP, ddRHS0dP.data(), NP, P, p0 ); 
//  for( unsigned i=0; i<NX+1; i++ ) std::cout << "dRHS[0]dP[" << i << "] = " << ddRHS0dP[i] << std::endl;

//  NIFTE.output( NIFTE.subgraph( 1, RHS ) );
//  std::ofstream o_RHS( "RHS[0].dot", std::ios_base::out );
//  NIFTE.dot_script( 1, RHS, o_RHS );
//  o_RHS.close();

//  exit(0);  

  int stat = OC->solve( Ip, p0 );
  cout << "\n" << std::scientific << std::setprecision(15);
  for( unsigned int i=0; i<NP; i++ )
     cout << "P[" << i << "] = " << OC->solution().p[i] << "  " << p0[i] << endl;
  cout << "\n";

  double *xk[2] = { 0, 0 }, *xpk[2] = { 0, 0 };
  OC->options.ODESLVS.DISPLAY   = 1;
  OC->states_FSA( OC->solution().p, xk, 0, xpk );
  CPPL::dgematrix JACF( NX, NX );
  for( unsigned i=0; i<NX; i++ )
    for( unsigned j=0; j<NX; j++ )
      JACF(i,j) = xpk[1][(j+1)*(NX+NQ)+i];
  std::cout << "Monodromy Matrix:\n" << JACF;
  std::vector<double> wr, wi;
  JACF.dgeev( wr, wi );
  std::cout << "Eignevalues of Monodromy Matrix:\n";
  for( unsigned i=0; i<NX; i++ )
    std::cout <<  "(" << wr[i] << ", " << wi[i] << ")\n";
 
#ifdef SAVE_TRAJECTORIES
  std::ofstream otraj( "test5b_CSS2.traj", std::ios_base::out );
  OC->options.ODESLVS.DISPLAY   = 1;	// 1; // 0;
  OC->options.ODESLVS.RESRECORD = 200000;	// number of sample points in the .dat file
  OC->set_time( 0., 10. );
  OC->setup();
  OC->states( p0 );
  //OC->states( OC->solution().p );//p0 );
  OC->record( otraj );
  otraj.close();
  return 0;
#endif

  const int LORD = -2, UORD = 1, PTORD = 40;
  for( int i=0; i<=(UORD-LORD)*PTORD; i++ ){

    Ip[12] = p0[12] = std::pow( 10., std::log10( L_f_0 ) + LORD + i/double(PTORD) ); 		
    cout << "L_f = " << p0[12] << endl;	
					
    OC->options.ODESLVS.DISPLAY = 0;
    int stat = OC->solve( Ip, p0 );
    if( stat ){
      p0[0] = 0.3;
      p0[1] = -0.3;
      p0[2] = 0.05;
      p0[3] = -0.06;
      p0[4] = 0.;
      p0[5] = -0.4;
	  stat = OC->solve( Ip, p0 );
    }

    // Optimal decision variables
    cout << "\n" << std::scientific << std::setprecision(15);
    for( unsigned int i=0; i<NP; i++ )
       cout << "P[" << i << "] = " << OC->solution().p[i] << endl;
    cout << "\n";

    OC->options.ODESLVS.DISPLAY   = 1;	// 1; // 0;
    //OC->set_time( 0., 50. );
    //OC->setup();
    OC->states( OC->solution().p, xk );
/*
    OC->options.ODESLVS.DISPLAY = 0;
    OC->setup();
    stat = OC->solve( Ip, OC->solution().p );
	
    // Optimal decision variables
    cout << "\n" << std::scientific << std::setprecision(15);
    for( unsigned int i=0; i<NP; i++ )
       cout << "P[" << i << "] = " << OC->solution().p[i] << endl;
    cout << "\n";

    //double* xk[2] = { 0, 0 };
    OC->options.ODESLVS.DISPLAY   = 1;	// 1; // 0;
    OC->options.ODESLVS.RESRECORD = 100;	// number of sample points in the .dat file
    //OC->set_time( 0., 50. );
    //OC->setup();
    OC->states( OC->solution().p, xk );
	OC->record( otraj );
*/
#ifdef SAVE_RESULTS
    if( !stat ){
      const unsigned IPREC = 9;
      ores << std::setw(4) << stat  
            << std::scientific << std::setprecision(IPREC);
      for( unsigned i=0; i<NP; i++ )
        ores	<< std::setw(IPREC+9) << OC->solution().p[i];
      for( unsigned i=0; i<NX+NQ; i++ )
	    ores << std::setw(IPREC+9) << xk[1][i];
      ores << std::setw(IPREC+9) << OC->solution().f;
      for( unsigned i=0; i<NX; i++ )
	    ores << std::setw(IPREC+9) << OC->solution().g[i];
      ores << std::endl;
   }
#endif

	// Update the initial conditions		
    for( unsigned i=0; i<NP; i++ )
      p0[i] = OC->solution().p[i];
  }	

  for( unsigned k=0; k<2; k++ )
    delete[] xk[k];
	
#ifdef SAVE_RESULTS
  ores.close();
#endif
  
  return 0;
}

