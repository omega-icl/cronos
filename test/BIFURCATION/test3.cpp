#define SAVE_RESULTS    // whether or not to save results to file
#define REPRESSILATOR
#undef USE_PROFIL	// specify to use PROFIL for interval arithmetic
#undef USE_FILIB	// specify to use FILIB++ for interval arithmetic
#undef DEBUG            // whether to output debug information

////////////////////////////////////////////////////////////////////////////////

#include <fstream>
#include "nlegpe.hpp"
//#include "stabh.hpp"

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



////////////////////////////////////////////////////////////////////////
#if defined( GOODWIN_REDUCED )
////////////////////////////////////////////////////////////////////////
  /*
  	Taken from: http://iopscience.iop.org/1478-3975/11/4/045002

	dx/dt = k1*(mc::pow(Ki,n)/(mc::pow(Ki,n)+mc::pow(z,n) )) - k2*x;
	dy/dt = k3*x - k4*y;
	dz/dt = k5*y - k6*z;
	
	reduced
	
	dx/dt = a/(1+mc::pow(z,n)) - x;
	dy/dt = x - y;
	dz/dt = y - z;
	
	a = k1*k3*k5/(mc::pow(k2,3)*Ki) \in I(1.5,2);
	n = 20 WTF!!!;

  */
  // Constant definitions
  const unsigned NP = 4 , NF = 3;  		// Parameter dimension
  // Ranges
  const I Ip[NP] = { I(1.5,2.) , I(0.,4.), I(0.,4.), I(0.,4.) };
  // Right Hand Side Function
  mc::FFVar RHSF
  ( unsigned int i, const mc::FFVar* P ){
  	// Parameters
  	int n  = 20;
   	// Variables
  	mc::FFVar a = P[0], x = P[1], y = P[2], z = P[3];
  	mc::FFVar result;    
  	// RHS
    switch( i ){
  		case 0: result =  a/(1.+mc::pow(z,n)) - x   ; break;
  		case 1: result =  x - y ; break;
  		case 2: result =  y - z ; break;
   		default: 
              std::cout << "Error in Function evaluation";
              break;
    }
    return result;
  };

////////////////////////////////////////////////////////////////////////
#elif defined( GOODWIN_REDUCED2 )
////////////////////////////////////////////////////////////////////////
  /*
  	Taken from: http://iopscience.iop.org/1478-3975/11/4/045002

	dx/dt = k1*(mc::pow(Ki,n)/(mc::pow(Ki,n)+mc::pow(z,n) )) - k2*x;
	dy/dt = k3*x - k4*y;
	dz/dt = k5*y - k6*z;
	
	reduced
	
	dx/dt = a/(1+mc::pow(z,n)) - x;
	dy/dt = x - y;
	dz/dt = y - z;
	
	a = k1*k3*k5/(mc::pow(k2,3)*Ki) \in I(1.5,2);
	n = 20 WTF!!!;

  */
  // Constant definitions
  const unsigned NP = 5 , NF = 4;  		// Parameter dimension
  // Ranges
  const I Ip[NP] = { I(1.5,2.) , I(1e-1,2.), I(1e-1,2.), I(1e-1,2.), I(1e-3,20.) };
  // Right Hand Side Function
  mc::FFVar RHSF
  ( unsigned int i, const mc::FFVar* P ){
  	// Parameters
  	int n  = 20;
   	// Variables
  	mc::FFVar a = P[0], x = P[1], y = P[2], z = P[3], w = P[4];
  	mc::FFVar result;    
  	// RHS
    switch( i ){
  		case 0: return  a/(1.+w) - x   ; break;
  		case 1: return  x - y ; break;
  		case 2: return  y - z ; break;
  		case 3: return  n*mc::log(z) - mc::log(w) ; break;
   		default: 
              std::cout << "Error in Function evaluation";
              break;
    }
    return 0.;
  };
////////////////////////////////////////////////////////////////////////
#elif defined( REPRESSILATOR )
////////////////////////////////////////////////////////////////////////
  /*
  	Taken from: http://www.nature.com/nature/journal/v403/n6767/pdf/403335a0.pdf

	dm_i/dt = -c2*m_i + ( c1/( Km + mc::pow(p_j,n) ) ) ; for i = lacl,tetR, cl; j = cl, lacl, tetR;
	dp_i/dt = -c4*p_i + c3*m_i;

  	mc::FFVar a      = P[0], b      = P[1], m_lacl = P[2], m_tetR = P[3], m_cl = P[4];
  	mc::FFVar p_lacl = P[5], p_tetr = P[6], p_cl   = P[7]; 

	c4 = 0.06 min^-1 -> protein half life 10-60 min.
	c3 = 0.16 min^-1 -> translation efficiency (might be increased to 20 proteins per transcript)
	c2 = 0.12 min^-1 -> mRNA degradation rate, (range 30 s to 50 min).
	c1 = 1-10 min^-1 -> calibrated
	Km = 1 
	
  */
  // Constant definitions
  const unsigned NP = 8 , NF = 6;  		// Parameter dimension
  // Ranges
  const I Ip[NP] = { I(1e-2,1e2), I(1e-2,1e2), I(0.,100.), I(0.,100.), I(0.,100.), I(0.,100.), I(0.,100.), I(0.,100.) };
  // Right Hand Side Function
  mc::FFVar RHSF
  ( unsigned int i, const mc::FFVar* P ){
  	// Parameters
  	int    n  = 2; 
  	double c1 = 1.;
//  	double c2 = 0.12; 
  	double c3 = 0.16; 
//  	double c4 = 0.06;
    double Km = 1.;
  	// Variables
  	mc::FFVar c2     = P[0], c4     = P[1], m_lacl = P[2], m_tetR = P[3], m_cl = P[4];
  	mc::FFVar p_lacl = P[5], p_tetR = P[6], p_cl   = P[7];  
  	mc::FFVar result;    
  	// RHS
    switch( i ){
  		case 0: result =  -c2*m_lacl + ( c1/( Km + mc::pow(p_cl,n) ) )   ; break;
  		case 1: result =  -c2*m_tetR + ( c1/( Km + mc::pow(p_lacl,n) ) ) ; break;
  		case 2: result =  -c2*m_cl   + ( c1/( Km + mc::pow(p_tetR,n) ) ) ; break;
   		case 3: result =  -c4*p_lacl + c3*m_lacl ; break;
      	case 4: result =  -c4*p_tetR + c3*m_tetR ; break;
      	case 5: result =  -c4*p_cl   + c3*m_cl   ; break;
      	default: 
              std::cout << "Error in Function evaluation";
              break;
    }
    return result;
  };

////////////////////////////////////////////////////////////////////////
#endif
////////////////////////////////////////////////////////////////////////

int main() 
{
  mc::FFGraph DAG;
  mc::FFVar T, P[NP], F[NF];
  T.set( &DAG );
  for( unsigned int i=0; i<NP; i++ ) P[i].set( &DAG );
  for( unsigned int i=0; i<NF; i++)  F[i] = RHSF( i , P );
  ////// Equilibrium Set /// 

  mc::NLEGPE<I> EQSet;
  EQSet.set_dag( &DAG );
  EQSet.set_indep( &T );
  EQSet.set_par( NP, P );
  EQSet.set_dep( NF, F );

  EQSet.options.SETINV.DISPLAY = 1;
  EQSet.options.SETINV.MAX_NODES = 30000;
  EQSet.options.SETINV.ABSOLUTE_TOLERANCE = 1e-7;
  EQSet.options.SETINV.RELATIVE_TOLERANCE = 1e-7;
  EQSet.options.SETINV.BRANCHING_VARIABLE_CRITERION = mc::SetInv<I>::Options::RGREL;
  EQSet.options.SETINV.MEASURE = mc::SetInv<I>::Options::LENGTH;

  EQSet.options.OUTPUT_BOUND = mc::NLEGPE<I>::Options::CM;
  EQSet.options.CM_ORDER     = 3;
  EQSet.options.OUTRED_MAX   = 10;
  EQSet.options.OUTRED_THRES = 2e-2; // 2e-2
  EQSet.options.OUTRED_TOL   = 1e-9; // 1e-9 
  EQSet.options.INRED_MAX    = 0;
  EQSet.options.INRED_THRES  = 2e-2;
  EQSet.options.ROTATION_USE = false;
  const double TOL = 1.e-6;
  typename std::list<mc::NLEGPE<I>::Data> Iym;
  for( unsigned int k=0; k<NF; k++ )
    Iym.push_back( typename mc::NLEGPE<I>::Data( -TOL, TOL, k ) );

  EQSet.solve( Ip, Iym, std::cout );

#if defined(SAVE_RESULTS )
  ofstream K_eq( "equilibrium.out", ios_base::out );
  EQSet.output_nodes( K_eq ); //, true );
  K_eq.close();
#endif  

  return 0;
}
