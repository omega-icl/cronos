#define SAVE_RESULTS    // whether or not to save results to file
#undef USE_PROFIL	// specify to use PROFIL for interval arithmetic
#undef USE_FILIB	// specify to use FILIB++ for interval arithmetic
#undef DEBUG            // whether to output debug information

////////////////////////////////////////////////////////////////////////////////

#include <fstream>
#include "nlegpe.hpp"

#ifdef USE_PROFIL
  #include "mcprofil.hpp"
  typedef INTERVAL I;
#else
  #ifdef USE_FILIB
    #include "mcfilib.hpp"
    typedef filib::interval<double> I;
  #else
    #include "interval.hpp"
    typedef mc::Interval I;
  #endif
#endif


mc::FFVar Leverrier_det( const unsigned int N, const mc::FFVar* A , mc::FFVar* c_array )
{
  mc::FFVar Bkm1[N*N];
  mc::FFVar Bk[N*N];
  //Initialization 
  for( unsigned int i=0; i<(N*N); i++ ){
  	Bkm1[i] = A[i] ; // Bk = Jacobian
  	Bk[i]   = 0.0  ;
  }  
  c_array[0] = 0.0;
  for( unsigned int i=0; i<N; i++ ){
  	c_array[0] += Bkm1[ i*N + i ];	// pcoeff_stack[0] = tr(Bk) 
  }
  // Main Loop  
  for( unsigned int l = 1; l<N; ++l ) { // for l = 1 to NF - 1
	// Auxmat = (B_{k-1} - c_{k_1}*I)  
	for( unsigned int i = 0; i<N; ++i ){
		Bkm1[ i*N + i ] -= c_array[l-1]; 
	}
	// B_k = Jacobian*Auxmat
	for( unsigned int i=0; i<N; i++ ) {
		for( unsigned int j=0; j<N; j++ ) {
			for( unsigned int k=0; k<N; k++ ) Bk[ i*N+j ] += (A[i*N+k] * Bkm1[k*N+j]) ;  
		}
	}
	// pcoeff_stack[l] = tr(Bk) 
	c_array[l] = 0.0;
	for( unsigned int i=0; i<N; i++ ){
		c_array[l] += Bk[ i*N + i ];	// pcoeff_stack[l] = tr(Bk) 
	}
	c_array[l] *= 1./double(l+1);
	// Reset C and Bk
	for( unsigned int i=0; i<(N*N); i++ ){
		Bkm1[i] = Bk[i];
		Bk[i]   = 0.0;
	} 
  }
  //Compute determinant
  mc::FFVar Determinant;
  if( N % 2 ) Determinant = c_array[N-1]; // odd
  else         Determinant = -1.0*c_array[N-1];
  return Determinant;
} 

int main() 
{
  mc::FFGraph DAG;

  const unsigned NP = 7, NF = 6;
  mc::FFVar T, P[NP], F[NF+1];
  T.set( &DAG );
  for( unsigned int i=0; i<NP; i++ ) P[i].set( &DAG );
  // Model Parameters
  double mu_1b = 1.2 , K_S1 = 7.1  , mu_2b = 0.74 , K_S2  = 9.28 ; 
  double K_I2  = 256., kLa  = 19.8 , K_H   = 16.  , P_t   = 1.   ; 
  double alpha = 0.5 , k1   = 42.14, k2    = 116.5, k3    = 268.0;
  double k4    = 50.8, k5   = 343.6, k6    = 453.0, S1_in = 5.   ;
  double S2_in = 80.0, Z_in = 50.0 , C_in  = 0.   ;  
  // Equilibrium Parameters
  mc::FFVar X1 = P[1], X2 = P[2], S1 = P[3], S2 = P[4], Z = P[5], C = P[6], D = P[0];
  // Auxiliary Equations
  mc::FFVar mu_1    = mu_1b * ( S1 / (S1 + K_S1) );
  mc::FFVar mu_2    = mu_2b * ( S2 / (S2 + K_S2 + S2*S2/K_I2 ) );
  mc::FFVar phi_CO2 = C + S2 - Z + K_H*P_t + (k6/kLa)*mu_2*X2;
  mc::FFVar P_CO2   = ( phi_CO2 - mc::sqrt( phi_CO2*phi_CO2 - 4.*K_H*P_t*( C + S2 - Z ) ) )/( 2.*K_H );
  mc::FFVar q_CO2   = kLa * ( C + S2 - Z - K_H * P_CO2 );
  // Equilibrium Conditions
  F[0] = ( mu_1 - alpha*D )*X1;
  F[1] = ( mu_2 - alpha*D )*X2;
  F[2] = D*( S1_in - S1 ) - k1*mu_1*X1;
  F[3] = D*( S2_in - S2 ) + k2*mu_1*X1 - k3*mu_2*X2;
  F[4] = Z_in  - Z; //D*( Z_in  - Z  );
  F[5] = D*( C_in  - C  ) - q_CO2 + k4*mu_1*X1 + k5*mu_2*X2; 
  const mc::FFVar* dFdX = DAG.FAD( NF,F, NP-1, P , false ); // Row-wise
  mc::FFVar c_array[NF];
  F[6] = Leverrier_det( NF , dFdX , c_array );


  mc::NLEGPE<I> EQSet;
  EQSet.set_dag( &DAG );
  EQSet.set_indep( &T );
  EQSet.set_par( NP-1, P );
  EQSet.set_dep( NF-1, F );

  EQSet.options.SETINV.DISPLAY = 2;
  EQSet.options.SETINV.MAX_NODES = 20000;
  EQSet.options.SETINV.ABSOLUTE_TOLERANCE = 1e-7;
  EQSet.options.SETINV.RELATIVE_TOLERANCE = 1e-7;
  EQSet.options.SETINV.BRANCHING_VARIABLE_CRITERION = mc::SetInv<I>::Options::RGREL;
  EQSet.options.SETINV.MEASURE = mc::SetInv<I>::Options::LENGTH;

  EQSet.options.OUTPUT_BOUND = mc::NLEGPE<I>::Options::CM;
  EQSet.options.CM_ORDER     = 3;
  EQSet.options.OUTRED_MAX   = 10;
  EQSet.options.OUTRED_THRES = 2e-2;
  EQSet.options.OUTRED_TOL   = 1e-9;
  EQSet.options.INRED_MAX    = 0;
  EQSet.options.INRED_THRES  = 2e-2;
  EQSet.options.ROTATION_USE = false;

  //const I Ip[NP] = { I(0.9,1.), I(0.,0.6), I(0.4,1.1), I(0.,5.), I(0.,25.), I(40.,50.), I(30.,55.) };
  const I Ip[NP] = { I(1e-3,1.5), I(0.,0.6), I(0.,0.8), I(0.,5.), I(0.,80.), I(40.,50.), I(30.,55.) };
  const double TOL = 1e-8;
  typename std::list<mc::NLEGPE<I>::Data> Iym;
  for( unsigned int k=0; k<NF-1; k++ )
    Iym.push_back( typename mc::NLEGPE<I>::Data( -TOL, TOL, k ) );

/*
  mc::CModel<I> CM( NP, 7 );
  mc::CVar<I> CV_P[NP], CV_F[NF+1];
  for( unsigned int i=0; i<NP; i++ ) CV_P[i].set( &CM, i, Ip[i] );
  //for( unsigned int i=0; i<NP; i++ ) std::cout << "R: " << CV_P[i].R() << "  B: " << CV_P[i].B() << std::endl;
  DAG.eval( NF, F, CV_F, NP, P, CV_P );
  for( unsigned int i=0; i<NF-1; i++ ) std::cout << "R: " << CV_F[i].R() << "  B: " << CV_F[i].B() << std::endl;
*/
  EQSet.solve( Ip, Iym, std::cout );

#if defined(SAVE_RESULTS )
  ofstream K_un( "undetermined_500000.out", ios_base::out );
  EQSet.output_stacks( K_un ); //, true );
  K_un.close();
#endif

////// Steady-State Bifurcation Set

////  mc::NLEGPE<I> SSBifurcation;
////  SSBifurcation.set_dag( &DAG );
////  SSBifurcation.set_indep( &T );
////  SSBifurcation.set_par( NP, P );
////  SSBifurcation.set_dep( NF+1, F );

////  SSBifurcation.options.SETINV.DISPLAY = 1;
////  SSBifurcation.options.SETINV.MAX_NODES = 30000;
////  SSBifurcation.options.SETINV.ABSOLUTE_TOLERANCE = 1e-6;
////  SSBifurcation.options.SETINV.RELATIVE_TOLERANCE = 1e-6;
////  SSBifurcation.options.SETINV.BRANCHING_VARIABLE_CRITERION = mc::SetInv<I>::Options::RGREL;
////  SSBifurcation.options.SETINV.MEASURE = mc::SetInv<I>::Options::LENGTH;

////  SSBifurcation.options.OUTPUT_BOUND = mc::NLEGPE<I>::Options::CM;
////  SSBifurcation.options.CM_ORDER     = 3;
////  SSBifurcation.options.OUTRED_MAX   = 10;
////  SSBifurcation.options.OUTRED_THRES = 2e-2;
////  SSBifurcation.options.OUTRED_TOL   = 1e-9;
////  SSBifurcation.options.INRED_MAX    = 0;
////  SSBifurcation.options.INRED_THRES  = 1e-2;
////  SSBifurcation.options.ROTATION_USE = false;

////  typename std::list<mc::NLEGPE<I>::Data> IymB;
////  for( unsigned int k=0; k<NF+1; k++ )
////    IymB.push_back( typename mc::NLEGPE<I>::Data( -TOL, TOL, k ) );

////  SSBifurcation.solve( Ip, IymB, std::cout );

////#if defined(SAVE_RESULTS )
////  std::list<I*> L_clus = SSBifurcation.clusters();
////  std::cout << "No clusters: " << L_clus.size() << std::endl;
////  ofstream K_cl( "clusters10000.out", ios_base::out );
////  SSBifurcation.output_clusters( K_cl );
////  K_cl.close();
////#endif



  return 0;
}
