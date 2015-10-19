const unsigned int NPM   = 7;	// <- Order of Taylor/Chebyshev model
#define USE_CMODEL		// <- Use Chebyshev models?
//#define MC__AEBND_DEBUG
#include "aebnd.hpp"

#include "interval.hpp"
typedef mc::Interval I;

#ifdef USE_CMODEL
  #include "cmodel.hpp" 
  typedef mc::CModel<I> PM;
  typedef mc::CVar<I> PV;
#else
  #include "tmodel.hpp"
  typedef mc::TModel<I> PM;
  typedef mc::TVar<I> PV;
#endif

int main()
{
  mc::FFGraph NLE;  // DAG describing the problem

  const unsigned NP = 22 , NF = 21 ;

  mc::FFVar T, P[NP], F[NF], SOL[NF];
  T.set( &NLE);
  for( unsigned int i=0; i<NP; i++ ) P[i].set( &NLE );
  // Model Parameters
  const unsigned nc = 3, ns = 3, pf = 1; 	//No. stages, feed location, No. components 
  double V = 1.38; 							//Vapour flow rate
  /*
  	let a[1] := 23.4832; 	let a[2] := 20.5110; 	let a[3] := 20.9064;
	let b[1] := -3634.01;	let b[2] := -2664.30;	let b[3] := -3096.52;
	let c[1] := -33.768;	let c[2] := -79.483;	let c[3] := -53.668;
	
	let r[1,1] := 0.0;		let r[1,2] :=  0.7411;	let r[1,3] :=  0.9645;	
	let r[2,1] := -1.0250;	let r[2,2] := 0.0;		let r[2,3] := -1.4350;
	let r[3,1] := -0.9645;	let r[3,2] :=  2.7470;	let r[3,3] := 0.0;

	let s[1,1] := 0.0;		let s[1,2] := -477.00;	let s[1,3] := -903.1024;
	let s[2,1] :=  72.78;	let s[2,2] := 0.0;		let s[2,3] :=  768.20;
	let s[3,1] := -140.9995;let s[3,2] := -1419.0;	let s[3,3] := 0.0;
  */
  double a[nc] = { 23.4832,  20.5110, 20.9064 }; 
  double b[nc] = { -3634.01, -2664.30, -3096.52}; 
  double c[nc] = { -33.768, -79.483, -53.668}; 
  double r[nc*nc] = { 0.0,  0.7411, 0.9645, -1.0250, 0.0, -1.4350, -0.9645, 2.7470, 0.0};
  double s[nc*nc] = { 0.0, -477.00, -903.1024, 72.78, 0.0, 768.20, -140.9995, -1419.0, 0.0};
  /*
    for {i in 1..C, j in 0..N}
	let f[i,j] := 0.0;

	let f[1, N_F] := 0.4098370;
	let f[2, N_F] := 0.01229769;
	let f[3, N_F] := 0.06090665;
  */
  double f[nc*(ns+1)];   
  for( unsigned i = 0; i < nc; ++i ){
  	for( unsigned j = 0; j < ns+1; ++j ) f[i*ns+j] = 0.0;
  }
  f[ 0*ns+pf ] = 0.4098370;
  f[ 1*ns+pf ] = 0.01229769;
  f[ 2*ns+pf ] = 0.06090665;

  double FF[ns+1], sumd;
  //for {j in 0..N} let F[j] := sum{i in 1..C} f[i,j];
  for( unsigned j = 0; j < ns+1; ++j ){
  	sumd = 0.0;
  	for( unsigned i = 0; i < nc; ++i ) sumd += f[ i*ns + j ];
  	FF[j] = sumd;
  } 
  // Equilibrium Parameters
  mc::FFVar D, X[nc*ns], K[nc*ns], Te[ns];
  D = P[0];
  for( unsigned i = 0; i < nc*ns ; ++i  ) X[i] = P[i+  1];
  for( unsigned i = 0; i < nc*ns ; ++i  ) K[i] = P[i+ 10];
  for( unsigned i = 0; i <    ns ; ++i  ) Te[i] = P[i+ 19];
////  for( unsigned i = 0; i <    NP ; ++i  ) std::cout << "P[" << i <<"]=" << P[i] << "\n";
////  for( unsigned i = 0; i < nc*ns ; ++i  ) std::cout << "X[" << i+1 <<"]=" << P[i+1] << "\n";
////  for( unsigned i = 0; i < nc*ns ; ++i  ) std::cout << "K[" << i+1 <<"]=" << P[i+10] << "\n";
////  for( unsigned i = 0; i <    ns ; ++i  ) std::cout << "Te[" << i+1 <<"]=" << P[i+19] << "\n";
  // Auxiliary Equations 
  mc::FFVar B, L[ns+1], p[nc*ns], rcp_T[ns], Lambda[nc*nc*ns], sum_xLambda[nc*ns];
  mc::FFVar rcp_sum_xLambda[nc*ns], gamma[nc*ns], sum;
  // param B := F[N_F] - D;
  B = FF[pf] - D;
  /*
  	param L{j in 0..N} = V - D + sum{k in 0..j} F[k];
  */
  for( unsigned j = 0; j < ns+1 ; ++j ){
  	sumd = 0.0;
  	for( unsigned k = 0; k <= j; ++k ) sumd += FF[k];
  	L[j] = V - D + sumd ;	
  }
  /*
  	var p{i in 1..C, j in 1..N} = exp(a[i]+b[i]/(T[j]+c[i]));
  */
  for( unsigned i = 0; i < nc; ++i  ){
  	for( unsigned j = 0; j < ns; ++j ) p[i*ns+j] = exp( a[i]+b[i]/(Te[j]+c[i]));  	
  }
  /*
  	var rcp_T{j in 1..N} = 1.0/T[j];
  */
  for( unsigned i = 0; i < ns; ++i ) rcp_T[i] = 1.0/Te[i]; 
  /*
  	var Lambda{i1 in 1..C, i2 in 1..C, j in 1..N} = exp(r[i1,i2]+s[i1,i2]*rcp_T[j]);
  */
  for( unsigned i1 = 0; i1 < nc; ++i1 ){
  	for( unsigned i2 = 0; i2 < nc; ++i2 ){
  		for( unsigned j = 0; j < ns; ++j ) Lambda[nc*ns*i1+ns*i2+j] = exp(r[nc*i1+i2]
  									+s[nc*i1+i2]*rcp_T[j]);
  	}
  }
  /*
  	var sum_xLambda{i in 1..C, j in 1..N} = sum{i1 in 1..C} (x[i1,j]*Lambda[i,i1,j]);
  */
  for( unsigned i = 0; i < nc; ++i  ){
  	for( unsigned j = 0; j < ns; ++j ){
  		sum = 0.0;
  		for( unsigned i1 = 0; i1 < nc; ++i1 )  sum += X[i1*ns+j]*Lambda[i*nc*ns+i1*ns+j];	
  		sum_xLambda[i*ns+j] = sum;
  	}
  }
  /*
  	var rcp_sum_xLambda{i in 1..C, j in 1..N} = 1.0/sum_xLambda[i,j];
  */
  for( unsigned i = 0; i < nc; ++i  ){
  	for( unsigned j = 0; j < ns; ++j ) rcp_sum_xLambda[i*ns+j] = 1.0/sum_xLambda[i*ns+j];
  }
  /*
  	var gamma{i in 1..C, j in 1..N} =
  	exp( -log(sum_xLambda[i,j]) + 1.0 - 
  	(sum{i2 in 1..C} (x[i2,j]*Lambda[i2,i,j]*rcp_sum_xLambda[i2,j])) );
  */
  for( unsigned i = 0; i < nc; ++i  ){
  	for( unsigned j = 0; j < ns; ++j ){
  	    sum = 0.0;
  		for( unsigned i2 = 0; i2 < nc; ++i2 ) sum += X[i2*ns+j]*Lambda[nc*ns*i2+ns*i+j]*rcp_sum_xLambda[i2*ns+j];
  		gamma[i*ns+j] = exp( -1.0*log( sum_xLambda[i*ns+j] ) + 1.0 - sum ); 
	}
  }
  // Equilibrium Conditions
  /*
  	AUXILIARY EQUATIONS
  	E_aux_K{j in 1..N, i in 1..C}: 	K[i,j] - gamma[i,j]*(p[i,j]/P) = 0.0;  
  */
  for( unsigned j = 0; j < ns; ++j  ){
  	for( unsigned i = 0; i < nc; ++i ){
  		F[ j*nc+i ] = K[i*ns+j]-gamma[i*ns+j]*( p[i*ns+j]/100000.0);
  		std::cout << "E_aux_K->F[" << j*nc+i <<"]=" << F[j*nc+i] << std::endl;
  	}
  }
  /*
  	MATERIAL BALANCES
  	M_tot{i in 1..C}: D*(K[i,1]*x[i,1]) + B*x[i,N] - f[i,N_F] = 0.0;
  */
  for( unsigned i = 0; i < nc; ++i  ){
  	F[ i+(ns*nc) ] = D*(K[i*ns+0]*X[i*ns+0]) + B*X[i*ns+(ns-1)] - f[i*ns+pf];
  	std::cout << "M_tot->F[" << i+(ns*nc) <<"]=" << F[i+(ns*nc)] << std::endl;
  } 
  /* 
  	NOTE THE UNUSUAL FORMULATION
	M_eq{j in 1..N-1, i in 1..C}:
	L[j]*x[i,j] + sum{i1 in j+1..N} f[i,i1] - B*x[i,N] - V*(K[i,j+1]*x[i,j+1]) = 0.0;
  */
  for( unsigned j = 0; j < ns-1; ++j  ){
  	for( unsigned i = 0; i < nc; ++i ){
  		sum = 0.0; 
  		for( unsigned i2 = j+1; i2 < ns; ++i2 ) sum += f[ i*ns+i2 ] ; 
  		F[ j*nc+i+(ns*nc)+nc ] = L[j]*X[i*ns+j] + sum - B*X[i*ns+(ns-1)] - V*( K[i*ns+j+1]*X[i*ns+j+1] ) ;
        std::cout << "M_eq->F[" << j*nc+i+(ns*nc)+nc <<"]=" << F[j*nc+i+(ns*nc)+nc] << std::endl;
  	} 
  } 
  /*
  	SUMMATION EQUATIONS
	S_x_eq{j in 1..N}: 	sum{i in 1..C} x[i,j] - 1.0 = 0.0;
  */
  for( unsigned j = 0; j < ns; ++j  ){
    F[ j+(ns*nc)+nc+(ns-1)*nc ] = 0.0;
  	for( unsigned i = 0; i < nc; ++i ) F[ j+(ns*nc)+nc+(ns-1)*nc ] += X[i*ns+j] ;
  	F[ j+(ns*nc)+nc+(ns-1)*nc ] -=  1.0;
	std::cout << "S_x_eq->F[" << j+(ns*nc)+nc+(ns-1)*nc <<"]=" << F[j+(ns*nc)+nc+(ns-1)*nc] << std::endl;
  } 

  //I Ip[1]  = { I(0.455 ,0.4550005 ) };
  //  Ix0[NF] = { I(0.850766,0.8507665), I(0.79568,0.795685	), I(0.421913,0.4219135), I(0.0303188,0.03031885), I(0.0528961,0.05289615), I(0.167538,0.1675385), I(0.118915,0.1189155), I(0.151424,0.1514245), I(0.41055,0.410555), I(1.02818,1.028185), I(1.07917,1.079175), I(1.90389,1.903895), I(0.5509,0.55095), I(0.488305,0.4883055), I(0.301822,0.3018225), I(0.912914,0.9129145), I(0.762759,0.7627595), I(0.356008,0.3560085), I(336.712,336.7125), I(336.983,336.9835), I(339.169,339.1695) },
  //  Ix[NF];

//  I Ip[1]   = { I(4.5781117e-01, 4.5875000e-01 ) },
//    Ix0[NF] = { I(8.4827231e-01, 8.4899345e-01), I(7.8822429e-01, 7.9040964e-01), I(3.6754895e-01, 3.8227192e-01),
//                I(3.1114743e-02, 3.1418465e-02), I(5.4773693e-02, 5.5512808e-02), I(1.7764506e-01, 1.8129180e-01),
//                I(1.1989499e-01, 1.2030255e-01), I(1.5483542e-01, 1.5624415e-01), I(4.4010512e-01, 4.5113715e-01),
//                I(1.0296346e+00, 1.0302134e+00), I(1.0847033e+00, 1.0869695e+00), I(2.0872690e+00, 2.1645833e+00),
//                I(5.4772557e-01, 5.4862903e-01), I(4.8087946e-01, 4.8299378e-01), I(2.9397662e-01, 2.9570445e-01),
//                I(9.0507351e-01, 9.0729411e-01), I(7.4567386e-01, 7.5051718e-01), I(3.3479993e-01, 3.3999198e-01),
//                I(3.3672058e+02, 3.3672392e+02), I(3.3700793e+02, 3.3701809e+02), I(3.3952994e+02, 3.3968132e+02) },

  I Ip[1]   = { I(4.57e-01, 4.6e-01 ) },
    Ix0[NF] = { I(8.4e-01, 8.5e-01), I(7.8e-01, 8.0e-01), I(3.6e-01, 3.9e-01),
                I(3.0e-02, 3.2e-02), I(5.4e-02, 5.6e-02), I(1.7e-01, 1.9e-01),
                I(1.1e-01, 1.3e-01), I(1.4e-01, 1.6e-01), I(4.4e-01, 4.6e-01),
                I(1.0e+00, 1.1e+00), I(1.0e+00, 1.1e+00), I(2.0e+00, 2.2e+00),
                I(5.4e-01, 5.5e-01), I(4.8e-01, 4.9e-01), I(2.9e-01, 3.0e-01),
                I(9.0e-01, 9.1e-01), I(7.4e-01, 7.6e-01), I(3.3e-01, 3.4e-01),
                I(3.2e+02, 3.4e+02), I(3.2e+02, 3.4e+02), I(3.2e+02, 3.4e+02) },


    Ix[NF];

  PM PMEnv( 1, NPM );
  PV PMp[1], PMx[NF], PMx0[NF];
  PMp[0].set( &PMEnv, 0, Ip[0] );
  for( unsigned i=0; i<NF; i++ ) PMx0[i] = Ix0[i];

  /////////////////////////////////////////////////////////////////////////
  // Bound AE solution set
  mc::AEBND<I,PM,PV> EX1;

  EX1.set_dag( &NLE );
  EX1.set_var( 1, P );
  EX1.set_dep( NF, &P[1], F );

  EX1.options.DISPLAY = 2;
  EX1.options.BLKDEC  = false;//true;
  EX1.options.MAXIT   = 20;
  EX1.options.RTOL    =
  EX1.options.ATOL    = 0e0;
  EX1.options.BOUNDER = mc::AEBND<I,PM,PV>::Options::GS;//KRAW;//GS;
  EX1.options.PRECOND = mc::AEBND<I,PM,PV>::Options::INVMD;//INVMD;//QRM;//NONE;
try{
  EX1.setup();
  //std::cout << "\nSuccessful? " << (EX1.solve( Ip, Ix, Ix0 )==mc::AEBND<I,PM,PV>::NORMAL?"Y\n":"N\n");
  std::cout << "\nSuccessful? " << (EX1.solve( PMp, PMx, Ix0 )==mc::AEBND<I,PM,PV>::NORMAL?"Y\n":"N\n");
}
catch( mc::AEBND<I,PM,PV>::Exceptions &eObj ){
    cerr << "Error " << eObj.ierr()
         << " in AEBND:" << endl
	 << eObj.what() << endl
         << "Aborts." << endl;
    return eObj.ierr();
  }

  return 0;
}


