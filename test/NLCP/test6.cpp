#define SAVE_RESULTS    // whether or not to save results to file
#define USE_PROFIL	    // specify to use PROFIL for interval arithmetic
#undef USE_FILIB	    // specify to use FILIB++ for interval arithmetic
#undef  USE_DEPs        // whether to use dependents
#define MC__USE_CPLEX   // whether to use CPLEX or GUROBI
#undef MC__CSEARCH_SHOW_BOXES
#undef MC__CSEARCH_SHOW_DEPS
#undef MC__CSEARCH_SHOW_REDUC
#undef MC__CSEARCH_SHOW_OUTER
#undef MC__SBP_SHOW_SCOREBRANCHING
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

////////////////////////////////////////////////////////////////////////
// Equilibrium problem: getting the equilibrium manifold of an AD model

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

  CP.options.MIPFILE     = "";//"test6.lp";
  CP.options.DISPLAY     = 2;
  CP.options.MIPDISPLAY  = 0;
  CP.options.MAXITER     = 0;
  CP.options.NODEMEAS    = mc::SBP<I>::Options::RELMAXLEN;
  CP.options.VARMEAS     = mc::NLCP<I>::Options::DEP;
  CP.options.CVTOL       = 1e-3;
  CP.options.FEASTOL     = 1e-6;
  CP.options.PREPROC     = true;
  CP.options.BLKDECUSE   = true;
  CP.options.BRANCHVAR   = mc::SBP<I>::Options::RGREL;
  CP.options.STGBCHDEPTH = 0;
  CP.options.STGBCHRTOL  = 1e-2;
  CP.options.DOMREDMAX   = 10;
  CP.options.DOMREDTHRES = 1e-1;
  CP.options.DOMREDBKOFF = 1e-6;
  CP.options.RELMETH     = mc::NLCP<I>::Options::CHEB;
  CP.options.CMODPROP    = 2;
  CP.options.CMODRED     = mc::NLCP<I>::Options::APPEND;
  CP.options.CMODDEPS    = 2;
  CP.options.CMODATOL    = 1e-8;
  CP.options.CMODRTOL    = 1e-2;
  CP.options.CMODWARMS   = false;//true;
  CP.options.AEBND.ATOL    = 
  CP.options.AEBND.RTOL    = 1e-8;
  CP.options.AEBND.DISPLAY = 0;
  std::cout << CP;

  //const I Ip[NP] = { I(0.8,1.1), I(0.,0.6), I(0.4,1.1), I(0.,5.), I(0.,25.), I(40.,50.), I(30.,55.) };
  //const I Ip[NP] = { I(.85,.9), I(0.,0.6), I(0.4,0.8), I(0.,5.), I(0.,40.), I(40.,50.), I(30.,55.) };
  //const I Ip[NP] = { I(0.5,1.2), I(0.,0.6), I(0.,0.8), I(0.,5.), I(0.,80.), I(40.,50.), I(0.,60.) };
  const I Ip[NP] = { I(1e-2,1.5), I(0.,0.6), I(0.,0.8), I(0.,5.), I(0.,80.), I(40.,50.), I(0.,60.) };
  //const I Ip[NP] = { I(0.4756250000000050,0.5221875000000074),
  //                   I(0.0000000000000000,0.0000000000000001),
  //                   I(0.5512859583286112,0.5750890733044568),
  //                   I(4.9999999999999955,5.0000000000000036),
  //                   I(4.4292381079482296,5.8778967258135291),
  //                   I(49.9999999999999786,50.000000000001214),
  //                   I(49.5210620393331097,51.9928120454771446) };

  try{  
    CP.setup();
    CP.solve( Ip );
    CP.stats.display();
  }
  catch(...){
    std::cout << "Exception during optimization" << std::endl;
  }

#if defined(SAVE_RESULTS )
  std::ofstream K_un( "test6.out", std::ios_base::out );
  CP.output_nodes( K_un );
  K_un.close();
#endif

  return 0;

  // Steady-State Bifurcation condition
  const mc::FFVar* jacF = DAG.FAD( NF, F, NP-1, P+1, false ); // Row-wise
  CP.add_ctr( mc::BASE_OPT::EQ, DAG.det( NF, jacF ) );

  try{  
    CP.options.CVTOL       = 1e-3;
    CP.setup();
    CP.solve( Ip );
    CP.stats.display();
  }
#ifdef MC__USE_CPLEX
  catch(IloException& e){
    std::cout << "Error code = " << e.getMessage() << std::endl;
#else
  catch(GRBException& e){
    std::cout << "Error code = " << e.getErrorCode() << std::endl;
    std::cout << e.getMessage() << std::endl;
#endif
  }
  catch(...){
    std::cout << "Exception during optimization" << std::endl;
  }

  auto L_clus = CP.clusters();
  std::cout << "No clusters: " << L_clus.size() << std::endl;
  for( auto it=L_clus.begin(); it!=L_clus.end(); ++it ) delete[] *it;

#if defined(SAVE_RESULTS )
  std::ofstream K_bi( "test6b.out", std::ios_base::out );
  CP.output_nodes( K_bi ); //, true );
  K_bi.close();

  std::ofstream K_cl( "test6c.out", std::ios_base::out );
  CP.output_clusters( K_cl ); //, true );
  K_cl.close();
#endif

  return 0;
}
