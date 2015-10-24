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
  const unsigned NP = 3, NF = 2;
  mc::FFVar T, P[NP], F[NF+1];
  T.set( &DAG );
  for( unsigned int i=0; i<NP; i++ ) P[i].set( &DAG );
  double alpha = 1.0, mu_0 = 1.0, k_s = 1.0, k_i = 1.0, k_1 = 1.0, s_in = 10.0;
  mc::FFVar X = P[0], S = P[1], D = P[2] ;
  mc::FFVar mu_H = mu_0*(S/(S+k_s+(mc::sqr(S)/k_i)));  
  F[0] = mu_H*X-alpha*D*X;
  F[1] = -k_1*mu_H*X+D*(s_in-S);
  const mc::FFVar* dFdX = DAG.FAD( NF,F, NP-1, P , false ); // Row-wise
  mc::FFVar c_array[NF];
  F[2] = Leverrier_det( NF , dFdX , c_array );

  mc::NLEGPE<I> EQSet;
  EQSet.set_dag( &DAG );
  EQSet.set_indep( &T );
  EQSet.set_par( NP, P );
  EQSet.set_dep( NF, F );

  EQSet.options.SETINV.DISPLAY = 1;
  EQSet.options.SETINV.MAX_NODES = 10000;
  EQSet.options.SETINV.ABSOLUTE_TOLERANCE = 1e-6;
  EQSet.options.SETINV.RELATIVE_TOLERANCE = 1e-6;
  EQSet.options.SETINV.BRANCHING_VARIABLE_CRITERION = mc::SetInv<I>::Options::RGREL;
  EQSet.options.SETINV.MEASURE = mc::SetInv<I>::Options::LENGTH;

  EQSet.options.OUTPUT_BOUND = mc::NLEGPE<I>::Options::CM;
  EQSet.options.CM_ORDER     = 3;
  EQSet.options.OUTRED_MAX   = 10;
  EQSet.options.OUTRED_THRES = 2e-2;
  EQSet.options.OUTRED_TOL   = 1e-9;
  EQSet.options.INRED_MAX    = 0;
  EQSet.options.INRED_THRES  = 1e-2;
  EQSet.options.ROTATION_USE = false;

  const I Ip[NP] = { I(0.,10.), I(0.,10.), I(1e-3,1.) };

  const double TOL = 1e-7;
  typename std::list<mc::NLEGPE<I>::Data> Iym;
  for( unsigned int k=0; k<NF; k++ )
    Iym.push_back( typename mc::NLEGPE<I>::Data( -TOL, TOL, k ) );

  EQSet.solve( Ip, Iym, std::cout );

#if defined(SAVE_RESULTS )
  ofstream K_un( "undetermined.out", ios_base::out );
  EQSet.output_nodes( K_un ); //, true );
  K_un.close();
#endif

// Steady-State Bifurcation Set

  mc::NLEGPE<I> SSBifurcation;
  SSBifurcation.set_dag( &DAG );
  SSBifurcation.set_indep( &T );
  SSBifurcation.set_par( NP, P );
  SSBifurcation.set_dep( NF+1, F );

  SSBifurcation.options.SETINV.DISPLAY = 1;
  SSBifurcation.options.SETINV.MAX_NODES = 10000;
  SSBifurcation.options.SETINV.ABSOLUTE_TOLERANCE = 1e-6;
  SSBifurcation.options.SETINV.RELATIVE_TOLERANCE = 1e-6;
  SSBifurcation.options.SETINV.BRANCHING_VARIABLE_CRITERION = mc::SetInv<I>::Options::RGREL;
  SSBifurcation.options.SETINV.MEASURE = mc::SetInv<I>::Options::LENGTH;

  SSBifurcation.options.OUTPUT_BOUND = mc::NLEGPE<I>::Options::CM;
  SSBifurcation.options.CM_ORDER     = 3;
  SSBifurcation.options.OUTRED_MAX   = 10;
  SSBifurcation.options.OUTRED_THRES = 2e-2;
  SSBifurcation.options.OUTRED_TOL   = 1e-9;
  SSBifurcation.options.INRED_MAX    = 0;
  SSBifurcation.options.INRED_THRES  = 1e-2;
  SSBifurcation.options.ROTATION_USE = false;

  typename std::list<mc::NLEGPE<I>::Data> IymB;
  for( unsigned int k=0; k<NF+1; k++ )
    IymB.push_back( typename mc::NLEGPE<I>::Data( -TOL, TOL, k ) );

  SSBifurcation.solve( Ip, IymB, std::cout );

#if defined(SAVE_RESULTS )
  std::list<I*> L_clus = SSBifurcation.clusters();
  std::cout << "No clusters: " << L_clus.size() << std::endl;
  ofstream K_cl( "clusters.out", ios_base::out );
  SSBifurcation.output_clusters( K_cl );
  K_cl.close();
#endif



  return 0;
}
