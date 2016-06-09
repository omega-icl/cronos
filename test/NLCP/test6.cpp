#define SAVE_RESULTS    // whether or not to save results to file
#undef  USE_PROFIL	// specify to use PROFIL for interval arithmetic
#undef  USE_FILIB	// specify to use FILIB++ for interval arithmetic
#undef  DEBUG		// whether to output debug information
#define USE_DEPS	// whether to use dependents
////////////////////////////////////////////////////////////////////////////////

#include <fstream>
#include "nlcp.hpp"

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
typedef mc::CVar<I> CVI;
//typedef I CVI;

////////////////////////////////////////////////////////////////////////////////
mc::FFVar Leverrier_det
( const unsigned int N, const mc::FFVar* A , mc::FFVar* c_array )
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
  if( N % 2 ) Determinant =  c_array[N-1]; // odd
  else        Determinant = -c_array[N-1];
  return Determinant;
} 

////////////////////////////////////////////////////////////////////////////////
// Find all solutions of the system of nonlinear equations:
//   (1.-R) * (D/10/(1+b1)-x) * exp(10*x/(1+10*x/g)) - x = 0
//   x - (1+b2)*y + (1.-R)*(D/10-b1*x-(1+b2)*y)*exp(10*y/(1+10*y/g)) = 0
// for (x,y) in [0,1]^2, R in [0.93,0.99], and g=1e3, D=22, b1=2, b2=2
////////////////////////////////////////////////////////////////////////////////

int main() 
{
  mc::FFGraph DAG;
  const unsigned NP = 7, NF = 6;
  mc::FFVar P[NP], F[NF];
  for( unsigned int i=0; i<NP; i++ ) P[i].set( &DAG );

  // Model Parameters
  double mu_1b = 1.2 , K_S1 = 7.1  , mu_2b = 0.74 , K_S2  = 9.28 ; 
  double K_I2  = 256., kLa  = 19.8 , K_H   = 16.  , P_t   = 1.   ; 
  double alpha = 0.5 , k1   = 42.14, k2    = 116.5, k3    = 268.0;
  double k4    = 50.8, k5   = 343.6, k6    = 453.0, S1_in = 5.   ;
  double S2_in = 80.0, Z_in = 50.0 , C_in  = 0.   ;  

  // Equilibrium Parameters
  mc::FFVar D = P[0];
  mc::FFVar X1 = P[1], X2 = P[2], S1 = P[3], S2 = P[4], Z = P[5], C = P[6];

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
  F[4] = D*( Z_in  - Z  ); //Z_in  - Z; //
  F[5] = D*( C_in  - C  ) - q_CO2 + k4*mu_1*X1 + k5*mu_2*X2; 

  // Constraint Projection
  mc::NLCP<I> CP;
  CP.set_dag( &DAG );
#ifndef USE_DEPS
  CP.set_var( NP, P );
  for( unsigned i=0; i<NF; i++ )
    CP.add_ctr( mc::BASE_OPT::EQ, F[i] );
#else
  CP.set_var( NP-NF, P );
  CP.set_dep( NF, P+NP-NF, F );
#endif

  //CP.options.MIPFILE   = "test6.lp";
  CP.options.DISPLAY     = 2;
  CP.options.MAXITER     = 1000;
  CP.options.CVATOL      = 1e-6;
  CP.options.CVRTOL      = 5e-3;
  CP.options.BRANCHVAR   = mc::SetInv<CVI>::Options::SCORES;//RGREL;
  CP.options.NODEMEAS    = mc::SetInv<CVI>::Options::MEANWIDTH;
  CP.options.STGBCHDEPTH = 0;
  CP.options.STGBCHDRMAX = 0;
  CP.options.STGBCHRTOL  = 1e-2;
  CP.options.DOMREDMAX   = 10;
  CP.options.DOMREDTHRES = 1e-2;
  CP.options.DOMREDBKOFF = 1e-7;
  CP.options.CTRBACKOFF  = 1e-5;
  CP.options.RELMETH     = mc::NLCP<I>::Options::CHEB;//HYBRID;//CHEB;
  CP.options.CMODSPAR    = true;
  CP.options.CMODPROP    = 2;
  CP.options.CMODCUTS    = 2;
  CP.options.CMREDORD    = 5;
  CP.options.CMREDTHRES  = 1e-4;
  CP.options.CMREDALL    = true;
  CP.options.CMODEL.MIXED_IA = true ;
  std::cout << CP;

  //const I Ip[NP] = { I(0.9,1.), I(0.,0.6), I(0.4,1.1), I(0.,5.), I(0.,25.), I(40.,50.), I(30.,55.) };
  //const I Ip[NP] = { I(.85,.9), I(0.,0.6), I(0.4,0.8), I(0.,5.), I(0.,40.), I(40.,50.), I(30.,55.) };
  const I Ip[NP] = { I(1e-2,1.5), I(0.,0.6), I(0.,0.8), I(0.,5.), I(0.,80.), I(40.,50.), I(0.,60.) };
  //const I Ip[NP] = { I(1.00000e-03,6.33054e-01), I(0.00000e+00,1.50000e-01), I(4.00000e-01,8.00000e-01),
  //                   I(1.25000e+00,2.50000e+00), I(0.00000e+00,2.00000e+01), I(5.00000e+01,5.00001e+01),
  //                   I(4.25000e+01,4.87500e+01) };
  //const I Ip[NP] = { I(1.00000e-03,6.33054e-01), I(0.00000e+00,1.50000e-01), I(4.00000e-01,8.00000e-01),
  //                   I(1.25000e+00,2.50000e+00), I(0.00000e+00,2.00000e+01), I(4.99999e+01,5.00001e+01),
  //                   I(3.62500e+01,4.875000e+01) };

  //std::set<unsigned> CMexcl;
  //CMexcl.insert( 5 );
  //CP.setup(CMexcl);
  CP.setup();
  CP.solve( Ip );

#if defined(SAVE_RESULTS )
  std::ofstream ores;
  ores.open( "test6.out", std::ios_base::out );
  CP.output_nodes( ores, 7 );
  ores.close();
#ifdef USE_DEPS
  ores.open( "test6b.out", std::ios_base::out );
  CP.output_nodes( ores, 20, 7 );
  ores.close();
#endif
#endif
/*
  // Steady-State Bifurcation condition
  const mc::FFVar* dFdX = DAG.FAD( NF, F, NP-1, P+1, false ); // Row-wise
  mc::FFVar c_array[NF];
  mc::FFVar SSEQ = Leverrier_det( NF , dFdX , c_array );
  CP.add_ctr( mc::BASE_OPT::EQ, SSEQ );

  const I Ip[NP] = { I(1e-2,1.5), I(0.,0.6), I(0.,0.8), I(0.,5.), I(0.,80.), I(40.,50.), I(0.,60.) };

  //CP.options.MIPFILE   = "test6.lp";
  CP.options.DISPLAY     = 2;
  CP.options.MAXITER     = 10000;
  CP.options.CVATOL      = 1e-6;
  CP.options.CVRTOL      = 1e-6;
  CP.options.BRANCHVAR   = mc::SetInv<CVI>::Options::RGREL;
  CP.options.NODEMEAS    = mc::SetInv<CVI>::Options::LENGTH;
  CP.options.STGBCHDEPTH = 0;
  CP.options.STGBCHDRMAX = 0;
  CP.options.STGBCHRTOL  = 1e-2;
  CP.options.DOMREDMAX   = 10;
  CP.options.DOMREDTHRES = 2e-2;
  CP.options.DOMREDBKOFF = 1e-8;
  CP.options.CTRBACKOFF  = 1e-4;
  CP.options.RELMETH     = mc::NLCP<I>::Options::HYBRID;//DRL;
  CP.options.CMODSPAR    = true;
  CP.options.CMODPROP    = 3;
  CP.options.CMODCUTS    = 2;
  CP.options.CMREDORD    = 0;
  CP.options.CMREDTHRES  = 1e-5;
  CP.options.CMREDALL    = true;
  CP.options.CMODEL.MIXED_IA = true ;
  std::cout << CP;

  CP.setup();
  CP.solve( Ip );

  auto L_clus = CP.clusters();
  std::cout << "No clusters: " << L_clus.size() << std::endl;
  for( auto it=L_clus.begin(); it!=L_clus.end(); ++it ) delete[] *it;

#if defined(SAVE_RESULTS )
  std::ofstream K_un( "test6c.out", std::ios_base::out );
  CP.output_nodes( K_un ); //, true );
  K_un.close();

  std::ofstream K_cl( "test6d.out", std::ios_base::out );
  CP.output_clusters( K_cl ); //, true );
  K_cl.close();
#endif
*/
  return 0;
}
