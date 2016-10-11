#define SAVE_RESULTS    // whether or not to save results to file
#undef  USE_PROFIL	// specify to use PROFIL for interval arithmetic
#undef  USE_FILIB	// specify to use FILIB++ for interval arithmetic
#undef  DEBUG		// whether to output debug information
#define USE_DEPS	// whether to use dependents
////////////////////////////////////////////////////////////////////////////////

#include <fstream>
#include "nlcp.hpp"//.backup"

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
  const unsigned NP = 4, NF = 3;
  mc::FFVar T, P[NP], F[NF];
  T.set( &DAG );
  for( unsigned int i=0; i<NP; i++ ) P[i].set( &DAG );

  // Model Parameters 
  double D = 0.02;
  double mu_NR = 0.5 , mu_RC = 0.2  , I_NR = 1.25 , I_RC = 0.333 ; 
  double k_NR  = 8.0 , k_RC  = 9.0  , m_R  = 0.025, m_C  = 0.01  ; 

  // Equilibrium Parameters
  mc::FFVar Nr = P[0];
  mc::FFVar N  = P[1], R = P[2], C = P[3];

  // Auxiliary Equations
  mc::FFVar r_N    = N / (N + k_NR);
  mc::FFVar r_R    = R / (R + k_RC);

  // Equilibrium Conditions
  F[0] = D * ( Nr - N ) - I_NR * r_N * R;
  F[1] = ( mu_NR * r_N - D - m_R ) * R - I_RC * r_R * C;
  F[2] = ( mu_RC * r_R - D - m_C ) * C;

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

  //CP.options.MIPFILE   = "test9.lp";
  CP.options.DISPLAY     = 1;
  CP.options.MAXITER     = 10000;
  CP.options.CVATOL      = 1e-5;
  CP.options.CVRTOL      = 1e-5;
  CP.options.BRANCHVAR   = mc::SetInv<CVI>::Options::RGREL;
  CP.options.NODEMEAS    = mc::SetInv<CVI>::Options::MEANWIDTH;
  CP.options.STGBCHDEPTH = 0;
  CP.options.STGBCHDRMAX = 0;
  CP.options.STGBCHRTOL  = 1e-2;
  CP.options.DOMREDMAX   = 10;
  CP.options.DOMREDTHRES = 2e-2;
  CP.options.DOMREDBKOFF = 1e-8;
  CP.options.CTRBACKOFF  = 1e-5;
  CP.options.RELMETH     = mc::NLCP<I>::Options::HYBRID;//CHEB;
//  CP.options.CMODSPAR    = true;//false;
  CP.options.CMODPROP    = 3;
  CP.options.CMODCUTS    = 2;
  CP.options.CMREDORD    = 3;
  CP.options.CMREDTHRES  = 1e-5;
  CP.options.CMREDALL    = true;
  CP.options.CMODEL.MIXED_IA = true ;
  std::cout << CP;

  const I Ip[NP] = { I(0,80), I(0,30), I(0,30), I(0,24) };

  CP.setup();
  CP.solve( Ip );

#if defined(SAVE_RESULTS )
  std::ofstream ores;
  ores.open( "test9.out", std::ios_base::out );
  CP.output_nodes( ores, 7 );
  ores.close();
#ifdef USE_DEPS
  ores.open( "test9b.out", std::ios_base::out );
  CP.output_nodes( ores, 20, 7 );
  ores.close();
#endif
#endif

  // Steady-State Bifurcation condition
  const mc::FFVar* dFdX = DAG.FAD( NF, F, NP-1, P+1, false ); // Row-wise
  mc::FFVar c_array[NF];
  mc::FFVar SSEQ = Leverrier_det( NF , dFdX , c_array );
  CP.add_ctr( mc::BASE_OPT::EQ, SSEQ );

  CP.options.DISPLAY     = 2;
  CP.options.CMREDORD    = 0;

  CP.setup();
  CP.solve( Ip );

  auto L_clus = CP.clusters();
  std::cout << "No clusters: " << L_clus.size() << std::endl;
  for( auto it=L_clus.begin(); it!=L_clus.end(); ++it ) delete[] *it;

#if defined(SAVE_RESULTS )
  std::ofstream K_un( "test9c.out", std::ios_base::out );
  CP.output_nodes( K_un ); //, true );
  K_un.close();

  std::ofstream K_cl( "test9d.out", std::ios_base::out );
  CP.output_clusters( K_cl ); //, true );
  K_cl.close();
#endif

  return 0;
}
