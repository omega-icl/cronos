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
typedef mc::CModel<I> CM;
typedef mc::CVar<I> CV;

#include "specbnd.hpp"
typedef mc::Specbnd<I> SB;

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
  else        Determinant = -1.0*c_array[N-1];
  return Determinant;
} 

int main() 
{
  mc::FFGraph DAG;

  ////// Equilibrium Set /// 

  const unsigned NP = 8, NF = 7;
  mc::FFVar T, P[NP], F[NF+1];
  T.set( &DAG );
  for( unsigned int i=0; i<NP; i++ ) P[i].set( &DAG );
  // Model Parameters
  double mu1b  = 1.2 , KS1 = 7.1  , mu2b = 0.74 , KS2  = 9.28 ; 
  double KI2   = 256., kLa = 19.8 , KH   = 16.  , Pt    = 1.   ; 
  double alpha = 0.5 , k1  = 42.14, k2   = 116.5, k3    = 268.0;
  double k4    = 50.8, k5  = 343.6, k6   = 453.0, S1in = 5.   ;
  //double S2in  = 80.0, Zin = 50.0 , Cin  = 20.0   ;  
  double S2in  = 80.0, Zin = 50.0 , Cin  = 0.0   ;  
  // Equilibrium Parameters
  mc::FFVar X1 = P[1], X2 = P[2], S1 = P[3], S2 = P[4], Z = P[5], C = P[6], PCO2 = P[7], D = P[0];
  // Auxiliary Equations
  mc::FFVar mu1    = mu1b * ( S1 / (S1 + KS1) );
  mc::FFVar mu2    = mu2b * ( S2 / (S2 + KS2 + S2*S2/KI2 ) );
  mc::FFVar CO2aq  = C + S2 - Z;
  mc::FFVar phiCO2 = CO2aq + KH*Pt + (k6/kLa)*mu2*X2;
  mc::FFVar qCO2   = kLa * ( CO2aq - KH * PCO2 );
  // Equilibrium Conditions
  F[0] = ( mu1 - alpha*D )*X1;				// dX1dt
  F[1] = ( mu2 - alpha*D )*X2;				// dX2dt
  F[2] = D*( S1in - S1 ) - k1*mu1*X1;			// dS1dt
  F[3] = D*( S2in - S2 ) + k2*mu1*X1 - k3*mu2*X2;	// dS2dt
  F[4] = Zin - Z; //D*( Zin - Z );			// dZdt
  F[5] = KH * PCO2*PCO2 - phiCO2 * PCO2 + CO2aq * Pt;   // definition of PCO2
  F[6] = D*( Cin - C ) - qCO2 + k4*mu1*X1 + k5*mu2*X2;	// dCdt 

  mc::NLEGPE<I> EQSet;
  EQSet.set_dag( &DAG );
  EQSet.set_indep( &T );
  EQSet.set_par( NP, P );
  EQSet.set_dep( NF, F );

  EQSet.options.SETINV.DISPLAY = 1;
  EQSet.options.SETINV.MAX_NODES = 20000;
  EQSet.options.SETINV.ABSOLUTE_TOLERANCE = 1e-7;
  EQSet.options.SETINV.RELATIVE_TOLERANCE = 1e-7;
  EQSet.options.SETINV.BRANCHING_VARIABLE_CRITERION = mc::SetInv<I>::Options::RGREL;
  EQSet.options.SETINV.MEASURE = mc::SetInv<I>::Options::LENGTH;

  EQSet.options.OUTPUT_BOUND = mc::NLEGPE<I>::Options::CM;
  EQSet.options.CM_ORDER     = 2;
  EQSet.options.OUTRED_MAX   = 10;
  EQSet.options.OUTRED_THRES = 2e-2;
  EQSet.options.OUTRED_TOL   = 1e-9;
  EQSet.options.INRED_MAX    = 0;
  EQSet.options.INRED_THRES  = 2e-2;
  EQSet.options.ROTATION_USE = false;

  const I Ip[NP] = { I(1e-3,1.5), I(0.,0.6), I(0.,0.8), I(0.,5.), I(0.,80.), I(Zin,Zin+1.), I(0.,80.), I(0.,Pt) };
  const double TOL = 1e-8;
  typename std::list<mc::NLEGPE<I>::Data> Iym;
  for( unsigned int k=0; k<NF; k++ )
    Iym.push_back( typename mc::NLEGPE<I>::Data( -TOL, TOL, k ) );

  EQSet.solve( Ip, Iym, std::cout );

#if defined(SAVE_RESULTS )
  ofstream K_eq( "equilibrium.out", ios_base::out );
  EQSet.output_nodes( K_eq ); //, true );
  K_eq.close();
  
  ////// Stability Analysis /// 
  
  const unsigned NP_ODE = 7, NF_ODE = 6 ;
  mc::FFVar F_ODE[NF_ODE];
  // Auxiliary Equations
  PCO2 = ( phiCO2 - mc::sqrt( phiCO2*phiCO2 - 4.*KH*Pt*CO2aq ) )/(2.*KH); 
  qCO2 = kLa * ( CO2aq - KH * PCO2 );
  // Equilibrium Conditions
  F_ODE[0] = F[0]; 						// dX1dt
  F_ODE[1] = F[1];						// dX2dt
  F_ODE[2] = F[2];						// dS1dt
  F_ODE[3] = F[3];						// dS2dt
  F_ODE[4] = D*(Zin - Z);					// dZdt
  F_ODE[5] = D*( Cin - C ) - qCO2 + k4*mu1*X1 + k5*mu2*X2;	// dCdt 
  const mc::FFVar* dFdX = DAG.FAD( NF_ODE, F_ODE, NF_ODE, P+1, false ); // Row-wise

  //I P_I[ NP_ODE ] ;
  CPPL::dgematrix dFdX_mid( NF_ODE, NF_ODE ); 
  I dFdX_I[ NF_ODE*NF_ODE ];
  CM cmod( NP_ODE, EQSet.options.CM_ORDER );
  CV P_CM[ NP_ODE ] ; 
  CV dFdX_CM[ NF_ODE*NF_ODE ];

  ofstream K_sta( "stable.out", ios_base::out );
  ofstream K_uns( "unstable.out", ios_base::out );  
  K_sta << std::scientific << std::setprecision(7);
  K_uns << std::scientific << std::setprecision(7); 

  SB::options.HESSBND = SB::Options::HERTZROHN;
  const mc::SetInv<I>::t_Nodes& EQ_Nodes = EQSet.open_nodes();
  mc::SetInv<I>::t_Nodes::const_iterator it = EQ_Nodes.cbegin() ;
  for( ; it != EQ_Nodes.cend(); ++it ){
    //for( unsigned i=0; i<NP_ODE; ++i ) P_I[i] = (*it)->P()[i];
    for( unsigned i=0; i<NP_ODE; ++i ) P_CM[i].set( &cmod, i, (*it)->P()[i] );
    try{
      //DAG.eval(NF_ODE*NF_ODE, dFdX, dFdX_I, NP_ODE, P, P_I);
      DAG.eval(NF_ODE*NF_ODE, dFdX, dFdX_CM, NP_ODE, P, P_CM);
      for( unsigned i=0; i<NF_ODE*NF_ODE; ++i ) dFdX_I[i] = dFdX_CM[i].B();
      std::pair<double,double> rebndHR = SB::spectral_bound_re( NF_ODE, dFdX_I );
      for( unsigned i=0; i<NF_ODE; ++i )
        for( unsigned j=0; j<NF_ODE; ++j ) dFdX_mid(i,j) = mc::Op<CV>::mid(dFdX_CM[i+NF_ODE*j]);
      std::cerr << "Node #" << (*it)->index() << ": eigL >= " << rebndHR.first << "  eigU <= " << rebndHR.second << std::endl;
/*      if( rebndHR.second <= 0. ){
        for( unsigned i=0; i<NP ; ++i )
          K_sta << mc::Op<I>::l((*it)->P()[i]) <<"\t"<< mc::Op<I>::u((*it)->P()[i]) << "  ";
        K_sta << "\n";
      }
      else if( rebndHR.first >= 0. ){
        for( unsigned i=0; i<NP ; ++i )
          K_uns <<  mc::Op<I>::l((*it)->P()[i]) <<"\t"<< mc::Op<I>::u((*it)->P()[i]) << "  ";
        K_uns << "\n";
      }
*/
      std::vector< double > wr, wi;
      //std::vector< CPPL::dcovector > vrr(), vri();
      dFdX_mid.dgeev( wr, wi );//, vrr, vri );      
      std::cerr << "Node #" << (*it)->index() << ": ";
      for( unsigned i=0; i<NF_ODE; ++i )
        std::cerr << "  " << wr[i];
      std::cerr << std::endl;
      { int dum; std::cin >> dum; }
      rebndHR = std::make_pair(wr[0],wr[0]);
      for( unsigned i=0; i<NF_ODE; ++i ){
        if( rebndHR.first  > wr[i] ) rebndHR.first  = wr[i];
        if( rebndHR.second < wr[i] ) rebndHR.second = wr[i];
      }
      if( rebndHR.second <= 0. ){
        for( unsigned i=0; i<NP ; ++i )
          K_sta << mc::Op<I>::l((*it)->P()[i]) <<"\t"<< mc::Op<I>::u((*it)->P()[i]) << "  ";
        K_sta << "\n";
      }
      else if( rebndHR.second >= 0. ){
        for( unsigned i=0; i<NP ; ++i )
          K_uns <<  mc::Op<I>::l((*it)->P()[i]) <<"\t"<< mc::Op<I>::u((*it)->P()[i]) << "  ";
        K_uns << "\n";
      }
    }
    catch(...){
      //std::cerr << "Node #" << (*it)->index() << ": Stability test failed (D = " << P_I[0] << ")\n";
      std::cerr << "Node #" << (*it)->index() << ": Stability test failed (D = " << P_CM[0].B() << ")\n";
    }	
  }

  K_sta.close(); 
  K_uns.close(); 
#endif

  return 0;
}
