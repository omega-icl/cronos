
#include "ffunc.hpp"
#include <iomanip>

#ifndef MC__STABH_HPP
#define MC__STABH_HPP

#undef NEVILLE_HURWITZ_DEBUG

namespace mc
{

/////////////////////////////////////////////////////////////////////////////////////////
void doPause()
/////////////////////////////////////////////////////////////////////////////////////////  
{
    std::cout << "Press [ Enter ] to continue...";
    std::cin.clear(); 
    std::cin.ignore();
}

/////////////////////////////////////////////////////////////////////////////////////////
enum STABILITY_TYPE
/////////////////////////////////////////////////////////////////////////////////////////  
{
  STABLE       = 0,
  UNSTABLE        ,
  UNDETERMINED
};

/////////////////////////////////////////////////////////////////////////////////////////
FFVar FL_Determinant( const unsigned int N, const FFVar* A , FFVar* c_array )
/////////////////////////////////////////////////////////////////////////////////////////
{
  FFVar Bkm1[N*N];
  FFVar Bk[N*N];
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
  for( unsigned int i=0; i<N; i++ ) c_array[i] *= (-1.0);
  //Compute determinant
  FFVar Determinant;
  if( N % 2 ) Determinant = -1.0*c_array[N-1]; // odd
  else        Determinant = c_array[N-1];
  return Determinant;
}

/////////////////////////////////////////////////////////////////////////////////////////
FFVar* Hurwitz_Matrix( const unsigned int N, FFVar* c_array, FFVar* H )
/////////////////////////////////////////////////////////////////////////////////////////
{
  for( unsigned int i=0; i<N; i++ ) {
  	for( unsigned int j=0; j<N; j++ ) {
		 if( int(2*j - i) == -1 )								H[ i*N+j ] = 1. ; // std::cout <<  1. << std::endl;
		 else if( int(2*j - i) < 0 || int(2*j - i) > int(N-1) ) H[ i*N+j ] = 0. ; //std::cout <<  0. << std::endl;
		 else 													H[ i*N+j ] = c_array[ 2*j - i ]; //std::cout << "a" << int( 2*j - i ) << std::endl;  
		}
	}
  return H;
}


/////////////////////////////////////////////////////////////////////////////////////////
FFVar* Row_Update( const unsigned int N, const unsigned int row, FFVar multiple, FFVar* H )
/////////////////////////////////////////////////////////////////////////////////////////
{
  for( unsigned int j=0; j<N; j++ ) {
  	H[ (row*N)+j ]	-= multiple*H[ (row-1)*N + j ]; 	
  }
  return H;
}

/////////////////////////////////////////////////////////////////////////////////////////
FFVar* NH_Elimination( const unsigned int N, FFVar* H )
/////////////////////////////////////////////////////////////////////////////////////////
{
#ifdef NEVILLE_HURWITZ_DEBUG
  for( unsigned int i=0; i<N*N; i++ ){
  	std::cout << "H[" << i << "]= " ;
  	if( H[ i ].id().first == FFVar::CREAL  ) std::cout << H[ i ].num();
	else std::cout << H[ i ] ;
	std::cout << std::endl;  
  }
#endif
  FFVar multiple_odd, multiple_even;
  for (unsigned int k=1; k<N; k++ ){ 
	  if( k % 2 ){
#ifdef NEVILLE_HURWITZ_DEBUG
	    std::cout << "multiple_odd = H[" << (N*k)+(k-1) << "]/H[" << N*(k-1)+(k-1) << "]= " ;
	    if( H[ (N*k)+(k-1) ].id().first == FFVar::CREAL  ) std::cout << H[ (N*k)+(k-1) ].num();
	    else std::cout << H[ (N*k)+(k-1) ] ;
	    std::cout << "/" ;
	    if( H[ (N*(k-1))+(k-1) ].id().first == FFVar::CREAL  )std::cout<< H[ (N*(k-1))+(k-1) ].num();
	    else std::cout<< H[ (N*(k-1))+(k-1) ];  
	    std::cout<< std::endl;
#endif
		multiple_odd = H[(N*k)+(k-1)]/H[N*(k-1)+(k-1)];
	  	for ( unsigned int i=k; i<N; i++ ) {
			// odd rows 
		  	if( i % 2 ) Row_Update( N, i, multiple_odd, H ); 	
		  }
	  }
	  else{
#ifdef NEVILLE_HURWITZ_DEBUG
	    std::cout << "multiple_odd = H[" << (N*k)+(k-1) << "]/H[" << N*(k-1)+(k-1) << "]= " ;
	    if( H[ (N*k)+(k-1) ].id().first == FFVar::CREAL  ) std::cout << H[ (N*k)+(k-1) ].num();
	    else std::cout << H[ (N*k)+(k-1) ] ;
	    std::cout << "/" ;
	    if( H[ (N*(k-1))+(k-1) ].id().first == FFVar::CREAL  )std::cout<< H[ (N*(k-1))+(k-1) ].num();
	    else std::cout<< H[ (N*(k-1))+(k-1) ];  
	    std::cout<< std::endl;
#endif
	    multiple_even = H[ (N*k)+(k-1) ] / H[ N*(k-1)+(k-1) ];
	  	for ( unsigned int i=k; i<N; i++ ){
	  		// even rows
	  		if ( !(i % 2)  ) Row_Update( N, i, multiple_even, H );
	  	}		
	  }
#ifdef NEVILLE_HURWITZ_DEBUG
	std::cout << N <<"-dimensional Hurwitz Matrix(" <<k<<"):\n" ;
	for( unsigned int i=0; i<N; i++ ){
		for( unsigned int j=0; j<N; j++ ){
			 std::cout<< "H[" << i*N+j << "] = " ;
			 if( H[ i*N+j ].id().first == FFVar::CREAL ) std::cout << H[ i*N+j ].num() << "\t" ;
			 else 	std::cout << H[ i*N+j ] << "\t" ;	
		}
		std::cout << std::endl;
	}
	std::cout << *(H->dag());
	doPause();
#endif	  
  }
  return H;
}

/////////////////////////////////////////////////////////////////////////////////////////
template <typename T>
unsigned int RHP_test( const unsigned int N, T* H_T, bool und_flag = false )
/////////////////////////////////////////////////////////////////////////////////////////
{
  bool exit_flag = false;
  int  steps     = 0;
  for( unsigned int i = 0; i < N ; i++){
  	for( unsigned int j = i; j < N; j++ ){
  		if( !(int(2*j - i) > int(N-1)) ){
	  		if( Op<T>::l(H_T[i*N+j]) <= 0. && Op<T>::u(H_T[i*N+j])<=0. ){ // unstable  
	  			exit_flag = true;
	  			break;
	  		}
	  		else if( Op<T>::l(H_T[i*N+j]) <= 0. && Op<T>::u(H_T[i*N+j]) >= 0. ){ // undetermined
	  			exit_flag = true;
	  			und_flag  = true;
	  			break;
	  		}
	  	}
  	}
  if( exit_flag == true || und_flag == true ) break;	
  steps++;		
  }
return steps;
}

/////////////////////////////////////////////////////////////////////////////////////////
template <typename T>
unsigned int RH_test( const unsigned int N, T* H_T, T* c_T, bool& und_flag )
/////////////////////////////////////////////////////////////////////////////////////////
{
  int  steps     = 0;
  bool coef_pos  = false;
  for( unsigned int i = 0; i < N ; i++){
  	if( Op<T>::l(c_T[i])     < 0. && 
	    Op<T>::u(c_T[i])     < 0.    ){
		und_flag = false;
		return steps;
	} 
    else if( Op<T>::l(c_T[i])     > 0. && 
		     Op<T>::u(c_T[i])     > 0.    ){
		coef_pos = true;
  	} 
	else{
  		und_flag  = true;
	}               
  }
  if( coef_pos || und_flag ){
  	for( unsigned int i = 0; i < N ; i++){
		if( Op<T>::l(H_T[i*N+i]) > 0. &&  
	        Op<T>::u(H_T[i*N+i]) > 0.    ){
	    	steps++;
 		} 
    	else if( Op<T>::l(H_T[i*N+i]) < 0. && 						
	  	         Op<T>::u(H_T[i*N+i]) < 0.    ){
	  		und_flag  = false;
	  	    return steps;
	  	}  	    
	  	else{ 			    											    
			und_flag  = true;
		}
  	}
  }
  return steps; 	
}

/////////////////////////////////////////////////////////////////////////////////////////
unsigned int RH_test( const unsigned int N, double* H_double, double* c_double, bool& und_flag )
/////////////////////////////////////////////////////////////////////////////////////////
{
  int  steps     = 0;
  bool coef_pos  = false;
  for( unsigned int i = 0; i < N ; i++){
	  if( c_double[i] < 0.    ){
//		  std::cout << "unstable\n";
		  return steps;
		  } 
		  
	  else if( c_double[i] > 0. ) coef_pos = true; 
	  else{
//	  		std::cout << "undetermined\n";
	  		und_flag  = true;
			return steps;	
	  }               
  }
  if( coef_pos ){
	  for( unsigned int i = 0; i < N ; i++){
	  		 if( H_double[i*N+i] > 0. ) steps++;
	  	else if( H_double[i*N+i] < 0. ) return steps;             // If the i-th diagonal element of the
	  	              	                                          // negative, the matrix is unstable.
	  	else{ 			    											    
				und_flag  = true;
				return steps;
		}
	  }
  }
  return steps; 	
}

}// end of namespace mc
#endif
