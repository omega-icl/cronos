#define SAVE_RESULTS    // whether or not to save results to file
#undef USE_PROFIL	// specify to use PROFIL for interval arithmetic
#undef USE_FILIB	// specify to use FILIB++ for interval arithmetic
#undef DEBUG            // whether to output debug information
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
using mc::cheb;

////////////////////////////////////////////////////////////////////////////////
// Plot an enclosure of the reachable set of the Lotka-Volterra system:
//
// @t = 1.0000e+01 :
// x[0] = 
//   a0    =  1.08115e+00     0
//   a1    =  1.68988e-01     1
//   a2    = -1.49733e-02     2
//   a3    = -8.83146e-03     3
//   a4    = -4.77404e-04     4
//   R     =  [ -3.60239e-02 :  3.60239e-02 ]
//   B     =  [  8.51853e-01 :  1.28050e+00 ]
//
// x[1] = 
//   a0    =  8.62378e-01     0
//   a1    =  7.79273e-02     1
//   a2    =  4.00543e-02     2
//   a3    = -3.75604e-04     3
//   a4    = -5.90753e-04     4
//   R     =  [ -3.41736e-02 :  3.41736e-02 ]
//   B     =  [  7.68232e-01 :  1.01550e+00 ]
////////////////////////////////////////////////////////////////////////////////

int main() 
{
  mc::FFGraph DAG;
  const unsigned NP = 3, NX = 2;
  mc::FFVar P[NP], *X=P+1, Pol[NX];
  for( unsigned int i=0; i<NP; i++ ) P[i].set( &DAG );
  Pol[0] = 1.08115e+00 * cheb(P[0],0) + 1.68988e-01 * cheb(P[0],1) 
         - 1.49733e-02 * cheb(P[0],2) - 8.83146e-03 * cheb(P[0],3)
         - 4.77404e-04 * cheb(P[0],4);
  Pol[1] = 8.62378e-01 * cheb(P[0],0) + 7.79273e-02 * cheb(P[0],1)
         + 4.00543e-02 * cheb(P[0],2) - 3.75604e-04 * cheb(P[0],3)
         - 5.90753e-04 * cheb(P[0],4);
  double Rem[2] = { 3.60239e-02, 3.41736e-02 };

  //const I Ip[NP] = { I(-1.,1.), I(8.51853e-01,1.28050e+00), I(7.68232e-01,1.01550e+00) };
  const I Ip[NP] = { I(-1.,1.), I(0.,10.), I(0.,10.) };

  // Image set with remainder

  mc::NLCP<I> CP;
  CP.set_dag( &DAG );
  CP.set_var( NP, P );
  for( unsigned i=0; i<NX; i++ ){
    CP.add_ctr( mc::BASE_OPT::GE, X[i] - Pol[i] + Rem[i] );
    CP.add_ctr( mc::BASE_OPT::LE, X[i] - Pol[i] - Rem[i] );
  }
  CP.setup();

  CP.options.MIPFILE     = "";//"test13.lp";
  CP.options.DISPLAY     = 1;
  CP.options.MAXITER     = 1000;
  CP.options.CVATOL      = 1e-4;
  CP.options.CVRTOL      = 1e-5;
  CP.options.BRANCHVAR   = mc::SetInv<CVI>::Options::RGREL;
  CP.options.STGBCHDEPTH = 10;
  CP.options.STGBCHRTOL  = 1e-2;
  CP.options.NODEMEAS    = mc::SetInv<CVI>::Options::MEANWIDTH;
  CP.options.DOMREDMAX   = 10;
  CP.options.DOMREDTHRES = 5e-2;
  CP.options.DOMREDBKOFF = 1e-8;
  CP.options.RELMETH     = mc::NLCP<I>::Options::CHEB;//DRL;
  CP.options.CMODPROP    = 4;
  CP.options.CMODCUTS    = 2;
  std::cout << CP;

  try{  
    CP.solve( Ip );
  }
  catch(GRBException e){
    std::cout << "Error code = " << e.getErrorCode() << std::endl;
    std::cout << e.getMessage() << std::endl;
  }
  catch(...){
    std::cout << "Exception during optimization" << std::endl;
  }

#if defined(SAVE_RESULTS )
  std::ofstream file( "test13.out", std::ios_base::out );
  CP.output_nodes( file );
  file.close();
#endif

  // Image set without remainder

  mc::NLCP<I> CP2;
  CP2.set_dag( &DAG );
  CP2.set_var( NP, P );
  for( unsigned i=0; i<NX; i++ ){
    CP2.add_ctr( mc::BASE_OPT::EQ, X[i] - Pol[i] );
  }
  CP2.setup();

  CP2.options.MIPFILE     = "";//"test13.lp";
  CP2.options.DISPLAY     = 1;
  CP2.options.MAXITER     = 1000;
  CP2.options.CVATOL      = 1e-4;
  CP2.options.CVRTOL      = 1e-5;
  CP2.options.BRANCHVAR   = mc::SetInv<CVI>::Options::RGREL;
  CP2.options.STGBCHDEPTH = 10;
  CP2.options.STGBCHRTOL  = 1e-2;
  CP2.options.NODEMEAS    = mc::SetInv<CVI>::Options::MEANWIDTH;
  CP2.options.DOMREDMAX   = 10;
  CP2.options.DOMREDTHRES = 5e-2;
  CP2.options.DOMREDBKOFF = 1e-8;
  CP2.options.RELMETH     = mc::NLCP<I>::Options::CHEB;//DRL;
  CP2.options.CMODPROP    = 4;
  CP2.options.CMODCUTS    = 2;
  std::cout << CP2;

  try{  
    CP2.solve( Ip );
  }
  catch(GRBException e){
    std::cout << "Error code = " << e.getErrorCode() << std::endl;
    std::cout << e.getMessage() << std::endl;
  }
  catch(...){
    std::cout << "Exception during optimization" << std::endl;
  }

#if defined(SAVE_RESULTS )
  std::ofstream file2( "test13b.out", std::ios_base::out );
  CP2.output_nodes( file2 );
  file2.close();
#endif

  return 0;
}
