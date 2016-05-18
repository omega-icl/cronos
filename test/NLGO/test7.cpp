#define USE_PROFIL
#define MC__USE_CPLEX

#include <fstream>
#include <iomanip>

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
#include "nlgo.hpp"

////////////////////////////////////////////////////////////////////////
// Parameter estimation problem (OSN permeability)

const unsigned NT = 6;
const double xR_C16[NT] = { 8.616e-3, 5.051e-2, 9.114e-2, 1.593e-1, 2.108e-1, 2.578e-1 }; // [kg/kg]
const double xP_C16[NT] = { 3.678e-3, 2.170e-2, 4.051e-2, 7.577e-2, 1.084e-1, 1.703e-1 }; // [kg/kg]
const double J[NT] = { 3.015e-1, 2.867e-1, 1.953e-1, 1.687e-1, 1.149e-1, 9.071e-2 }; // [kg/m2·s·bar]
const double v_C7 = 1.475e-4, v_C16 = 2.928e-4; // [m3/mol]
const double T  = 298e0; // [K]
const double dP = 30e5;  // [Pa]
const double R  = 8.314; // [J/K·mol]

template <class U>
U Jmod
( const U&k, double xR, double xP, double nu, double act=1. )
{
  return k * ( xR - act * xP * exp( -nu * dP / R / T ) );
}

////////////////////////////////////////////////////////////////////////
int main()
////////////////////////////////////////////////////////////////////////
{
  mc::FFGraph DAG;
  enum Species { C7=0, C16 };
  const unsigned NP = 2; mc::FFVar P[NP];
  for( unsigned i=0; i<NP; i++ ) P[i].set( &DAG );

  mc::NLGO<I> NLP;
  NLP.options.POLIMG.SANDWICH_MAXCUT = 7;
  NLP.options.POLIMG.SANDWICH_ATOL   = NLP.options.POLIMG.SANDWICH_RTOL  = 1e-3;
  NLP.options.CVATOL = NLP.options.CVRTOL = 1e-7;
  //NLP.options.MIPABSGAP = NLP.options.MIPRELGAP = 1e-9;
  //NLP.options.MAXITER = 0;
/*
  NLP.options.POLIMG.SANDWICH_MAXCUT = 10;
  NLP.options.POLIMG.SANDWICH_ATOL   = NLP.options.POLIMG.SANDWICH_RTOL  = 1e-5;
  NLP.options.POLIMG.BREAKPOINT_TYPE = mc::PolImg<I>::Options::BIN;//SOS2;//NONE;
  NLP.options.POLIMG.DCDECOMP_SCALE  = false;//true;
  //NLP.options.MIPFILE = "test4.lp";
  NLP.options.MIPDISPLAY = 0;
  NLP.options.DISPLAY = 2;
  NLP.options.NLPSLV.DISPLAY = 0;
  NLP.options.MAXITER = 21;
  NLP.options.PREPROC = true;
  NLP.options.DOMREDMAX = 10;
  NLP.options.PRESOS2BIGM = -1;
  NLP.options.CVATOL    = NLP.options.CVRTOL    = 1e-5;
*/
  NLP.set_dag( &DAG );  // DAG
  NLP.set_var( NP, P ); // decision variables


  mc::FFVar LSQ=0.;
  for( unsigned k=0; k<NT; k++ ){
    LSQ += mc::sqr( Jmod( P[C16], xR_C16[k], xP_C16[k], v_C16 ) - J[k] * xP_C16[k] )
         + (mc::sqr( Jmod( P[C7], 1.-xR_C16[k], 1.-xP_C16[k], v_C7 ) - J[k] * (1.-xP_C16[k]) ));
  }
  NLP.set_obj( mc::BASE_OPT::MIN, LSQ );
  // KKT cuts
  const mc::FFVar* dLSQdP = DAG.BAD( 1, &LSQ, NP, P );
  //for( unsigned i=0; i<NK; i++ )
  //  NLP.add_ctr( mc::BASE_NLP::EQ, dLSQdP[i] );
  NLP.setup();
  std::cout << NLP;

  I Ip[NP];
  Ip[C7] = I(0.,10.); Ip[C16] = I(0.,10.);

  double p0[NP];
  p0[C7] = 1.; p0[C16] = 1.;

  // Local optimization
  NLP.options.NLPSLV.DISPLAY = 5;
  NLP.options.NLPSLV.MAXITER = 100;
  NLP.local( Ip, p0 );
  std::cout << std::setprecision(5);
  for( unsigned i=0; i<NP; i++ )
    std::cout << "  p(" << i << ") = " << NLP.get_local_solution().p[i];
  std::cout << std::endl;
  double dLSQdp0[NP];
  for( unsigned i=0; i<NP; i++ ){
    DAG.eval( NP, dLSQdP, dLSQdp0, NP, P, NLP.get_local_solution().p );
    std::cout << "  dLSQdP(" << i << ") = " << dLSQdp0[i];
  }
  std::cout << std::endl;
  { int dum; std::cin >> dum; }

  // Global optimization using SBB
  NLP.options.DISPLAY = 2;
  //NLP.options.MIPFILE = "test5.lp";
  NLP.options.NLPSLV.DISPLAY = 0;
  NLP.options.CSALGO  = mc::NLGO<I>::Options::SBB;
  NLP.options.RELMETH = mc::NLGO<I>::Options::HYBRID;//CHEB;//
  NLP.options.CMODPROP = 3;
  NLP.options.CMODSPAR = true;

  NLP.solve( Ip, 0, p0 );
  std::cout << std::fixed << std::setprecision(1);
  std::cout << "POLIMG: " << NLP.stats.tPOLIMG << " sec\n";
  std::cout << "LPSET:  " << NLP.stats.tLPSET << " sec\n";
  std::cout << "LPSOL:  " << NLP.stats.tLPSOL << " sec, " << NLP.stats.nLPSOL << " problems\n";

  return 0;
}
