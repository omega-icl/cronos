// with WARM-START
// Varying L_f only

#undef SAVE_RESULTS		// <- Whether to save bounds to file

#include <fstream>
#include <iomanip>
#include "doseqslv_ipopt.hpp"

using std::cout;
using std::endl;
typedef mc::Interval I;

////////////////////////////////////////////////////////////////////////
int main()
{
#ifdef SAVE_RESULTS
  std::ofstream ores( "test5_4.dat", std::ios_base::out );
#endif
  std::ofstream otraj( "test5_CSS.traj", std::ios_base::out );


  mc::FFGraph NIFTE;

  // constant parameters  
  const double P0 = 1.013e5;		// Pa
  const double U0 = 8.00e-4;		// m3/s
  const double tau = 5.;				// s
  const double Lambda = 330.;
  const double K = 1.37;				// control gain
  
  const double R_l = 4.0809e6;
  const double R_f = 2.1307e6;
  const double R_th = 5.0212e8;
  const double L_d = 1.8027e5;
  const double L_l = 1.2710e7;
  const double L_p = 3.7739e5;
  const double C_ad = 1.7618e-9;
  const double C_d = 7.3513e-8;
  const double C_p = 7.4280e-8;
	
  // Varying the value of R_f -- nominal value = 2.1307e6
  // const double power_low = 0.;	
  // const double power_up = 1.;
  // const double step = 0.01;	// step size
  // const int NSTEP = 0;// (power_up-power_low)/step;
  const double L_f_nom = 4.7428e6;
	
  const double L_f_low = 4.8728e6;
  const double L_f_up = 4.7428e6;
  const double step = -0.01e6;	// step size
  const int NSTEP = (L_f_up-L_f_low)/step;
  double L_f;	 

  const unsigned NP = 6;    // Number of parameters
  mc::FFVar P[NP];             // Parameters
  for( unsigned int i=0; i<NP; i++ ) P[i].set( &NIFTE );
  mc::FFVar &TF = P[0], &Pad_0 = P[1], &Pd_0 = P[2], &Pp_0 = P[3], &Uf_0 = P[4], &Up_0 = P[5];

  const unsigned NX = 5;       // Number of states
  mc::FFVar X[NX];             // States
 for( unsigned int i=0; i<NX; i++ ) X[i].set( &NIFTE );
  mc::FFVar &P_ad = X[0], &P_d = X[1], &P_p = X[2], &U_f = X[3], &U_p = X[4]; 
  mc::FFVar RHS[NX];        // Right-hand side function
		 	
  const unsigned NQ = 3;
  mc::FFVar Q[NQ];
  for( unsigned int i=0; i<NQ; i++ ) Q[i].set( &NIFTE );
  mc::FFVar QUAD[NQ];  
	
  mc::FFVar U_l, P_l, U_ad, U_th, P_th, P_f;
	
  mc::FFVar IC[NX], CTRF[NX];

  // Initial conditions for the first value of R_f
  //double p0[NP];
  //p0[0] = 0.5;   			
  //p0[1] = -0.3;//-0.36;		
  //p0[2] = 0.05;			
  //p0[3] = -0.06;  		
  //p0[4] = 0.;
  //p0[5] = -0.4;//-0.42;  
  //double p0[NP] = { 5.621991178e-01, -3.689528172e-01,  5.389414282e-02, -6.109921183e-02, 0.000000000e+00, -4.389440962e-01 };
  //double p0[NP] = { 5.643141691e-01, -3.693161457e-01,  5.428951305e-02, -6.149191448e-02, 0.000000000e+00, -4.388824295e-01 };
  double p0[NP] = { 5.658952689e-01, -3.695822734e-01,  5.458591504e-02, -6.178617235e-02, 0.000000000e+00, -4.388350828e-01 };

  I Ip[NP];
  Ip[0] = I( 0.1, 1.0 );
  Ip[1] = I( -0.8, 0.8 ); 
  Ip[2] = I( -0.2, 0.2 ); 
  Ip[3] = I( -0.2, 0.2 ); 
  Ip[4] = I( 0., 0. );
  Ip[5] = I( -2.5, 2.5 ); 

  for (int j=NSTEP/2; j<=NSTEP; j++){
  //for (int j=NSTEP/2; j<=NSTEP; j++){
  //for (int j=NSTEP/2; j>=0; j--){

    //L_f = L_f_nom * 3.63078;//pow( 10, power_low+j*step );
    L_f = L_f_low + j*step;
		
    cout << "L_f = " << L_f << endl;	
		
    RHS[0] = TF * ( ( tau*(K*((mc::exp(2*Lambda*P_d)-1.)/(mc::exp(2*Lambda*P_d)+1.))-P_ad) )/(R_th*C_ad) + ( tau*U0*(U_f + U_p) )/(P0*C_ad) );
    RHS[1] = TF * ( ( tau*U0*U_f )/(P0*C_d) );
    RHS[2] = TF * ( ( tau*U0*U_p )/(P0*C_p) );
    RHS[3] = TF * ( ( (tau*P0/U0)*( L_l*P_p-L_p*P_ad-(L_p+L_l)*P_d ) - tau*(L_l*R_f+L_p*(R_f+R_l))*U_f - tau*L_p*R_l*U_p ) / ( L_l*(L_d+L_f)+L_p*(L_d+L_f+L_l) ) );
    RHS[4] = TF * ( ( (tau*P0/U0)*( L_l*P_d - (L_d+L_f)*P_ad - (L_d+L_f+L_l)*P_p ) + tau*(L_l*R_f-(L_d+L_f)*R_l)*U_f - tau*(L_d+L_f)*R_l*U_p ) / ( L_l*(L_d+L_f)+L_p*(L_d+L_f+L_l) ) );
	 	
    U_l = -(U_f + U_p)*U0;
    P_l = U_l * R_l;
    U_ad = C_ad * RHS[0]*P0/(tau*TF);
    U_th = U_ad + U_l;
    P_th = R_th * U_th + P_ad*P0;
    P_f = R_f * U_f * U0;
	
    QUAD[0] = P_l * U_l * tau;		// Average power dissipated in the load
    QUAD[1] = P_th * U_th * tau;	// Average power input
    QUAD[2] = P_f * U_f *U0 * tau;		// Average power dissipated in the feedback tube
			
    IC[0] = Pad_0;
    IC[1] = Pd_0;
    IC[2] = Pp_0;
    IC[3] = Uf_0;
    IC[4] = Up_0;
	
    CTRF[0] = P_ad - Pad_0;
    CTRF[1] = P_d - Pd_0;
    CTRF[2] = P_p - Pp_0;
    CTRF[3] = U_f - Uf_0;
    CTRF[4] = U_p - Up_0;		
					
    Ipopt::SmartPtr<mc::DOSEQSLV_IPOPT> OC = new mc::DOSEQSLV_IPOPT;

    OC->set_dag( &NIFTE );
    OC->set_time( 0., 1. );
    OC->set_parameter( NP, P );
    OC->set_state( NX, X );
    OC->set_differential( NX, RHS );
    OC->set_quadrature( NQ, QUAD, Q );
    OC->set_initial( NX, IC );

    OC->set_obj( mc::DOSEQSLV_IPOPT::MAX, Q[0]/Q[1] );            		   // objective
    OC->add_ctr( mc::DOSEQSLV_IPOPT::EQ, std::make_pair( 0, CTRF[0] ) );   // constraint #1
    OC->add_ctr( mc::DOSEQSLV_IPOPT::EQ, std::make_pair( 0, CTRF[1] ) );   // constraint #2
    OC->add_ctr( mc::DOSEQSLV_IPOPT::EQ, std::make_pair( 0, CTRF[2] ) );   // constraint #3
    OC->add_ctr( mc::DOSEQSLV_IPOPT::EQ, std::make_pair( 0, CTRF[3] ) );   // constraint #4	
    OC->add_ctr( mc::DOSEQSLV_IPOPT::EQ, std::make_pair( 0, CTRF[4] ) );   // constraint #5

    OC->options.DISPLAY   = 5;
    OC->options.MAXITER   = 500;
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

    double* xk[2] = { 0, 0 };
    OC->options.ODESLVS.DISPLAY   = 1;	// 1; // 0;
    OC->options.ODESLVS.RESRECORD = 100;	// number of sample points in the .dat file
    OC->states( p0, xk );
	//OC->record( otraj );
    //return 0.;

    OC->options.ODESLVS.DISPLAY = 0;
    int stat = OC->solve( Ip, p0 );
	
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

#ifdef SAVE_RESULTS
	OC->record( otraj );
    const unsigned IPREC = 9;
    ores << std::setw(4) << stat  
          << std::scientific << std::setprecision(IPREC)
	      << std::setw(IPREC+9) << L_f;
    for( unsigned i=0; i<NP; i++ )
      ores	<< std::setw(IPREC+9) << OC->solution().p[i];
    for( unsigned i=0; i<NX+NQ; i++ )
	  ores << std::setw(IPREC+9) << xk[1][i];
    ores << std::setw(IPREC+9) << OC->solution().f;
    for( unsigned i=0; i<NX; i++ )
	  ores << std::setw(IPREC+9) << OC->solution().g[i];
    ores << std::endl;
#endif
    for( unsigned k=0; k<2; k++ )
      delete[] xk[k];
		
	// Update the initial conditions		
    //for( unsigned i=0; i<NP; i++ )
    //  p0[i] = OC->solution().p[i];
		
	//for(unsigned i=0; i<NP; i++)	cout << "p0[" << i << "] = " << p0[i] << endl;	
  }	
	
#ifdef SAVE_RESULTS
  otraj.close();
  ores.close();
#endif



    mc::ODESLV_SUNDIALS IVP;
    IVP.options = OC->options.ODESLVS;
    IVP.set_dag( &NIFTE );
    IVP.set_time( 0., 1. );
    IVP.set_parameter( NP, P );
    IVP.set_state( NX, X );
    IVP.set_differential( NX, RHS );
    IVP.set_quadrature( NQ, QUAD, Q );
    IVP.set_initial( NX, IC );
    IVP.states( OC->solution().p );

  
  return 0;
}

