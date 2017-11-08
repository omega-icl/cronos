#define SAVE_RESULTS    // whether or not to save results to file
#define USE_PROFIL	    // specify to use PROFIL for interval arithmetic
#undef  USE_FILIB	    // specify to use FILIB++ for interval arithmetic
#undef  USE_DEPS        // whether to use dependents
#define MC__USE_CPLEX   // whether to use CPLEX or GUROBI
#undef  MC__CSEARCH_SHOW_BOXES
#undef  MC__CSEARCH_SHOW_DEPS
#undef  MC__CSEARCH_SHOW_REDUC
#undef  MC__CSEARCH_SHOW_OUTER
#undef  MC__SBP_SHOW_SCOREBRANCHING
////////////////////////////////////////////////////////////////////////////////

#include <fstream>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include "nlcp.hpp"
#include "scmodel.hpp"

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

///////////////////////////////////////////////////////////////////////////////
// Set-Membership Estimation: BOD model wtih noisy simulated measurements

const double Tmax = 8.; // [day]
const unsigned NDAT = 16; // 2^n
const unsigned NP = 2;

const double varYm = 1.;
const double conf  = 0.90; // % confidence level
//const double Chi2  = gsl_cdf_chisq_Pinv( conf, (double)NDAT );
const double Chi2  = gsl_cdf_chisq_Pinv( conf, (double)NP );

const double preal[NP] = { 20., 0.5 };  // Assumed "true" parameter values
I Ip[NP+NDAT+1] = { I(0.,40.), I(0.,3.) }; // Parameter range
double p0[NP+NDAT+1];

const unsigned int NDATMAX = 128;//1024;
char fname[50];
std::ofstream ofile;

////////////////////////////////////////////////////////////////////////////////
template <class U>
U Ymod
( const U*P, double T )
{
  return P[0]*(1.-exp(-P[1]*T));
}

////////////////////////////////////////////////////////////////////////////////
std::set<unsigned> two_phase_branching
( const typename mc::NLCP<I>::NODE*pNode )
{
  std::set<unsigned> sset;
  bool thres = true;
  for( unsigned i=0; i<NP; i++ ){
    sset.insert(i);
    //if( mc::Op<I>::diam(pNode->P(i)) > 4e-3*mc::Op<I>::diam(pNode->P0(i)) )
      thres = false;
  }
  if( !thres ) return sset;
  
  for( unsigned i=0; i<NDAT; i++ ){
    sset.insert(NP+i);
  }
  return sset;
/*
  unsigned int NT = NP+NDAT;
  double thres[NT];
  thres[0] = 5e-2; thres[1] = 5e-3;
  for( unsigned i=NP; i<NT; i++ ) thres[i] = 5e-2;

  std::set<unsigned> sset;
  for( unsigned i=0; i<NT; i++ )
    if( mc::Op<I>::diam(pNode->P(i)) > thres[i] ) sset.insert(i);  
  return sset;
*/
}

////////////////////////////////////////////////////////////////////////////////
int main() 
{

  // Noisy data generation

  double GaussE[NDATMAX];
  gsl_rng *r;
  r = gsl_rng_alloc( gsl_rng_taus );
  gsl_rng_set( r, NDATMAX );
  for(unsigned k=0; k<NDATMAX; k++){ GaussE[k] = gsl_ran_gaussian_ziggurat( r, varYm ); }
  gsl_rng_free( r );
  
  double Tm[NDAT], Ym[NDAT];
  for( unsigned k=0; k<NDAT; k++){
    Tm[k] = Tmax*(double)(k+1)/(double)NDAT;
    Ym[k] = Ymod( preal, Tm[k] ) + GaussE[(k+1)*NDATMAX/NDAT-1];
  }

  std::cout << std::fixed << std::setprecision(3) << std::right;
  std::cout << std::setfill('_') << std::setw(8*NDAT+2) << " " << std::endl << std::endl << std::setfill(' ');
  for( unsigned k=0; k<NDAT; k++) std::cout << std::setw(8) << Tm[k];
  std::cout << std::endl;
  for( unsigned k=0; k<NDAT; k++) std::cout << std::setw(8) << Ym[k];
  std::cout << std::endl;
  std::cout << std::setfill('_') << std::setw(8*NDAT+2) << " " << std::endl << std::endl << std::setfill(' ');

  // DAG
  
  mc::FFGraph DAG;
  mc::FFVar P[NP+NDAT+1], *Em = P+NP, &lambda = P[NP+NDAT];
  for( unsigned i=0; i<NP+NDAT+1; i++ ) P[i].set( &DAG );

  mc::FFVar J = 0., Jerr = 0., Eerr = 0.;
  mc::FFVar MMk[NDAT], dgk[NP*NDAT];
  for( unsigned k=0; k<NDAT; k++ ){
    J += sqr( Ym[k] - Ymod( P, Tm[k] ) );
    Jerr += sqr( Ym[k] + Em[k] - Ymod( P, Tm[k] ) );
    Eerr += sqr( Em[k] );
  }
  const mc::FFVar* dJdP = DAG.BAD( 1, &J, NP, P );
  const mc::FFVar* dJerrdP = DAG.BAD( 1, &Jerr, NP, P );
  mc::FFVar Lerr = Jerr + lambda * ( Eerr - Chi2*varYm );
  const mc::FFVar* dLerrdE = DAG.BAD( 1, &Lerr, NDAT, Em );

  // Bayesian Estimation with a Flat Prior

#if defined(SAVE_RESULTS )
  sprintf( fname, "test7_%2.0f%s_Nm%d_CRDATA.out",100.*conf,"%",NDAT );
  ofile.open( fname, std::ios_base::out );
  ofile << std::scientific << std::setprecision(6);
#endif
  const unsigned NPTS = 500;
  std::multimap<double,std::vector<double>> Jsamp;
  std::vector<double> dP(NP);
  double dJ;
  for( unsigned i=0; i<=NPTS; i++ ){
    dP[0] = mc::Op<I>::l(Ip[0]) + i/double(NPTS)*mc::Op<I>::diam(Ip[0]);
    for( unsigned j=0; j<=NPTS; j++ ){
      dP[1] = mc::Op<I>::l(Ip[1]) + j/double(NPTS)*mc::Op<I>::diam(Ip[1]);
      DAG.eval( 1, &J, &dJ, NP, P, dP.data() );
#if defined(SAVE_RESULTS )
      ofile << dP[0] << "  " << dP[1] << "  " << dJ << std::endl;
#endif
      Jsamp.insert( std::make_pair( dJ, dP ) );
    }
#if defined(SAVE_RESULTS )
    ofile << std::endl;
#endif
  }
#if defined(SAVE_RESULTS )
  ofile.close();
#endif
  double Jsum = 0., Jconf = 0.;
  auto it=Jsamp.begin();
  for( ; it!=Jsamp.end(); ++it )
    Jsum += std::exp(-it->first);
  for( it=Jsamp.begin(); Jconf < conf*Jsum; ++it )
    Jconf += std::exp(-it->first);
  const double JBAYES = it->first;
  std::cout << "JBAYES = " << JBAYES << std::endl;

  { int dum; std::cout << "PAUSED--"; std::cin >> dum; }

  // Bayesian Credibility Region

  mc::NLCP<I> SMMLE4;
  SMMLE4.set_dag( &DAG );
  SMMLE4.set_var( NP, P );
  SMMLE4.add_ctr( mc::BASE_OPT::LE, J - JBAYES );

  SMMLE4.options.DISPLAY     = 1;
  SMMLE4.options.MAXITER     = 0;
  SMMLE4.options.NODEMEAS    = mc::SBP<I>::Options::RELMAXLEN;
  SMMLE4.options.CVTOL       = 1e-3;
  SMMLE4.options.FEASTOL     = 1e-6;
  SMMLE4.options.BLKDECUSE   = false;
  SMMLE4.options.BRANCHVAR   = mc::SBP<I>::Options::RGREL;
  SMMLE4.options.DOMREDMAX   = 10;
  SMMLE4.options.DOMREDTHRES = 1e-1;
  SMMLE4.options.DOMREDBKOFF = 1e-6;
  SMMLE4.options.RELMETH     = mc::NLCP<I>::Options::CHEB;
  SMMLE4.options.CMODPROP    = 3;
  SMMLE4.options.CMODCUTS    = 2;
  SMMLE4.options.CMODDEPS    = 0;

  SMMLE4.setup();
  SMMLE4.solve( Ip );
  SMMLE4.stats.display();

#if defined(SAVE_RESULTS )
  sprintf( fname, "test7_%2.0f%s_Nm%d_CR.out",100.*conf,"%",NDAT );
  ofile.open( fname, std::ios_base::out );
  SMMLE4.output_nodes( ofile );
  ofile.close();
#endif

  { int dum; std::cout << "PAUSED--"; std::cin >> dum; }
  //return 0;

  // Least-Squares Estimation

  mc::NLGO<I> LSE;
  LSE.set_dag( &DAG );
  LSE.set_var( NP, P );
  LSE.set_obj( mc::BASE_OPT::MIN, J );
  for( unsigned i=0; i<NP; i++ )
    LSE.add_ctr( mc::BASE_OPT::EQ, dJdP[i], true );
 
  LSE.options.NLPSLV.DISPLAY = 5;
  LSE.options.NLPSLV.MAXITER = 200;
  LSE.setup();
  LSE.local( Ip );
/*
  LSE.options.NLPSLV.DISPLAY = 0;
  LSE.options.DISPLAY        = 2;
  LSE.options.CVATOL         = 1e-6;
  LSE.options.CVRTOL         = 1e-3;
  LSE.options.RELMETH        = mc::NLGO<I>::Options::CHEB;//DRL;
  LSE.options.CMODPROP       = 3;
  LSE.options.CMODCUTS       = 2;
  LSE.setup();
  LSE.solve( Ip );
*/
  const double JLS = LSE.get_local_solution().f;
  std::cout << std::setprecision(8);
  for( unsigned i=0; i<NP; i++ )
    std::cout << "  p(" << i << ") = " << LSE.get_local_solution().p[i];
  std::cout << std::endl << std::endl;
  std::cout << "JLS = " << JLS << std::endl;

  { int dum; std::cout << "PAUSED--"; std::cin >> dum; }
  //return 0;

  // Log-Likelihood Ratio Confidence Region

  mc::NLCP<I> SMMLE3;
  SMMLE3.set_dag( &DAG );
  SMMLE3.set_var( NP, P );
  SMMLE3.add_ctr( mc::BASE_OPT::LE, J - JLS - Chi2*varYm );
  std::cout << "JLR = " << JLS+Chi2*varYm << std::endl;

  SMMLE3.options.DISPLAY     = 1;
  SMMLE3.options.MAXITER     = 0;
  SMMLE3.options.NODEMEAS    = mc::SBP<I>::Options::RELMAXLEN;
  SMMLE3.options.CVTOL       = 1e-3;
  SMMLE3.options.FEASTOL     = 1e-6;
  SMMLE3.options.BLKDECUSE   = false;
  SMMLE3.options.BRANCHVAR   = mc::SBP<I>::Options::RGREL;
  SMMLE3.options.DOMREDMAX   = 10;
  SMMLE3.options.DOMREDTHRES = 1e-1;
  SMMLE3.options.DOMREDBKOFF = 1e-6;
  SMMLE3.options.RELMETH     = mc::NLCP<I>::Options::CHEB;
  SMMLE3.options.CMODPROP    = 3;
  SMMLE3.options.CMODCUTS    = 2;
  SMMLE3.options.CMODDEPS    = 0;

  SMMLE3.setup();
  SMMLE3.solve( Ip );
  SMMLE3.stats.display();

#if defined(SAVE_RESULTS )
  sprintf( fname, "test7d_%2.0f%s_Nm%d.out",100.*conf,"%",NDAT );
  ofile.open( fname, std::ios_base::out );
  SMMLE3.output_nodes( ofile );
  ofile.close();
#endif

  { int dum; std::cout << "PAUSED--"; std::cin >> dum; }
  return 0;

  // Worst-case upper bound: Interval set

  LSE.set_var( NP+NDAT, P );
  LSE.reset_ctr();
  mc::FFVar Serr = 0.;
  for( unsigned k=0; k<NDAT; k++ ){
    Serr += P[2+k]; Ip[2+k] = I(0e0,1e5); // eta
    LSE.add_ctr(mc::BASE_OPT::GE, P[2+k] - mc::sqr( Ym[k] - std::sqrt(Chi2*varYm) - Ymod( P, Tm[k] ) ));
    LSE.add_ctr(mc::BASE_OPT::GE, P[2+k] - mc::sqr( Ym[k] + std::sqrt(Chi2*varYm) - Ymod( P, Tm[k] ) ));
  }
  LSE.set_obj( mc::BASE_OPT::MIN, Serr );

  for( unsigned i=0; i<NP; i++ ){
    //Ip[i] = LSE.get_local_solution().p[i]; // Parameters
    p0[i] = preal[i];//LSE.get_local_solution().p[i];
  }
  for( unsigned k=0; k<NDAT; k++ ){
    Ip[NP+k] = I(0.,1e2); // errors
    p0[NP+k] = 0.;
  }

  LSE.options.NLPSLV.DISPLAY = 5;
  LSE.options.NLPSLV.MAXITER = 200;
  LSE.setup();
  LSE.local( Ip, p0 );
/*
  LSE.options.NLPSLV.DISPLAY = 0;
  LSE.options.DISPLAY        = 2;
  LSE.options.CVATOL         = 1e-6;
  LSE.options.CVRTOL         = 1e-3;
  LSE.options.BLKDECUSE      = true;
  LSE.options.RELMETH        = mc::NLGO<I>::Options::CHEB;//DRL;
  LSE.options.CMODPROP       = 3;
  LSE.options.CMODCUTS       = 2;
  LSE.options.AEBND.ATOL     = 
  LSE.options.AEBND.RTOL     = 1e-8;
  LSE.options.AEBND.DISPLAY  = 0;
  //LSE.options.AEBND.BOUNDER  = mc::NLGO<I>::t_AEBND::Options::ALGORITHM::GE;
  //LSE.options.AEBND.BLKDEC   = mc::NLGO<I>::t_AEBND::Options::DECOMPOSITION::DIAG;
  LSE.setup();
  LSE.solve( Ip );//, 0, p0 );
*/
  std::cout << std::setprecision(8);
  for( unsigned i=0; i<NP; i++ )
    std::cout << "  P(" << i << ") = " << LSE.get_local_solution().p[i];
  for( unsigned i=0; i<NDAT; i++ )
    std::cout << "  Em(" << i << ") = " << LSE.get_local_solution().p[NP+i];
  std::cout << std::endl << std::endl;
  std::cout << "U = " << LSE.get_local_solution().f << std::endl;

  { int dum; std::cout << "PAUSED--"; std::cin >> dum; }

  // Worst-case upper bound: Ellipsoidal set

  LSE.set_var( NP+NDAT+1, P );
  LSE.set_obj( mc::BASE_OPT::MIN, Jerr );
  LSE.reset_ctr();
  LSE.add_ctr( mc::BASE_OPT::EQ, Eerr - Chi2*varYm );
  for( unsigned k=0; k<NDAT; k++ )
    LSE.add_ctr( mc::BASE_OPT::EQ, dLerrdE[k] );

  for( unsigned i=0; i<NP; i++ ){
    //Ip[i] = LSE.get_local_solution().p[i]; // Parameters
    p0[i] = preal[i];//LSE.get_local_solution().p[i];
  }
  for( unsigned k=0; k<NDAT; k++ ){
    Ip[NP+k] = std::sqrt(Chi2*varYm) * I(-1.,1.); // errors
    p0[NP+k] = 0.;
  }
  Ip[NP+NDAT] = I{-1e2,-varYm*varYm}; // lambda
  p0[NP+NDAT] = -varYm*varYm;

  LSE.options.NLPSLV.DISPLAY = 5;
  LSE.options.NLPSLV.MAXITER = 200;
  LSE.setup();
  LSE.local( Ip, p0 );
/*
  LSE.options.NLPSLV.DISPLAY = 0;
  LSE.options.DISPLAY        = 2;
  LSE.options.CVATOL         = 1e-6;
  LSE.options.CVRTOL         = 1e-3;
  LSE.options.BLKDECUSE      = true;
  LSE.options.RELMETH        = mc::NLGO<I>::Options::CHEB;//DRL;
  LSE.options.CMODPROP       = 3;
  LSE.options.CMODCUTS       = 2;
  LSE.options.AEBND.ATOL     = 
  LSE.options.AEBND.RTOL     = 1e-8;
  LSE.options.AEBND.DISPLAY  = 0;
  //LSE.options.AEBND.BOUNDER  = mc::NLGO<I>::t_AEBND::Options::ALGORITHM::GE;
  //LSE.options.AEBND.BLKDEC   = mc::NLGO<I>::t_AEBND::Options::DECOMPOSITION::DIAG;
  LSE.setup();
  LSE.solve( Ip );//, 0, p0 );
*/
  const double JU = LSE.get_local_solution().f;
  std::cout << std::setprecision(8);
  for( unsigned i=0; i<NP; i++ )
    std::cout << "  P(" << i << ") = " << LSE.get_local_solution().p[i];
  for( unsigned i=0; i<NDAT; i++ )
    std::cout << "  Em(" << i << ") = " << LSE.get_local_solution().p[NP+i];
  std::cout << "  lambda = " << LSE.get_local_solution().p[NP+NDAT];
  std::cout << std::endl << std::endl;
  std::cout << "JU = " << JU << std::endl;

  { int dum; std::cout << "PAUSED--"; std::cin >> dum; }

  // Set-Membership Regression Region - Interval bounds

  LSE.set_var( NP+NDAT, P );
  LSE.reset_ctr();
  LSE.add_ctr( mc::BASE_OPT::LE, Eerr - Chi2*varYm );
  LSE.add_ctr( mc::BASE_OPT::LE, J - JU );
  //LSE.add_ctr( mc::BASE_OPT::LE, Jerr - JU );
  for( unsigned k=0; k<NP; k++ )
    LSE.add_ctr( mc::BASE_OPT::EQ, dJerrdP[k] );
  std::vector<I> Ip_BND( Ip, Ip+NP+NDAT );

  for( unsigned ip=0; ip<2*NP; ip++ ){
    if( ip%2 ) LSE.set_obj( mc::BASE_OPT::MIN, P[ip/2] );
    else       LSE.set_obj( mc::BASE_OPT::MAX, P[ip/2] );

    LSE.options.DISPLAY        = 0;
    LSE.options.NLPSLV.DISPLAY = 0;
    LSE.options.NLPSLV.MAXITER = 200;
    LSE.setup();
    LSE.local( Ip_BND.data(), p0 );
/*
    LSE.options.NLPSLV.DISPLAY = 0;
    LSE.options.DISPLAY        = 2;
    LSE.options.CVATOL         = 1e-6;
    LSE.options.CVRTOL         = 1e-3;
    LSE.options.BLKDECUSE      = true;
    LSE.options.RELMETH        = mc::NLGO<I>::Options::CHEB;//DRL;
    LSE.options.CMODPROP       = 3;
    LSE.options.CMODCUTS       = 2;
    LSE.options.AEBND.ATOL     = 
    LSE.options.AEBND.RTOL     = 1e-8;
    LSE.options.AEBND.DISPLAY  = 0;
    //LSE.options.AEBND.BOUNDER  = mc::NLGO<I>::t_AEBND::Options::ALGORITHM::GE;
    //LSE.options.AEBND.BLKDEC   = mc::NLGO<I>::t_AEBND::Options::DECOMPOSITION::DIAG;
    LSE.setup();
    LSE.solve( Ip_BND.data() );//, 0, p0 );
*/
    std::cout << std::setprecision(8);
    if( ip%2 ){
      Ip_BND[ip/2] = I( LSE.get_local_solution().f, mc::Op<I>::u(Ip_BND[ip/2]) );
      std::cout << "  pL(" << ip/2 << ") = " << LSE.get_local_solution().f << std::endl;
    }
    else{
      Ip_BND[ip/2] = I( mc::Op<I>::l(Ip_BND[ip/2]),LSE.get_local_solution().f );
      std::cout << "  pU(" << ip/2 << ") = " << LSE.get_local_solution().f << std::endl;
    }
  }

  { int dum; std::cout << "PAUSED--"; std::cin >> dum; }

#if defined(SAVE_RESULTS )
  sprintf( fname, "test7c_%2.0f%s_Nm%d.out",100.*conf,"%",NDAT );
  ofile.open( fname, std::ios_base::out );
  ofile << std::scientific << std::setprecision(8);
  for( unsigned ip=0; ip<NP; ip++ )
    ofile << mc::Op<I>::l(Ip_BND[ip]) << "  " << mc::Op<I>::u(Ip_BND[ip]) << "  ";
  ofile.close();
#endif

  // Set-Membership Regression Region - Exact

  mc::NLCP<I> SMMLE;
  SMMLE.set_dag( &DAG );

#ifndef USE_DEPS
  SMMLE.set_var( NP+NDAT, P );
  SMMLE.reset_ctr();
  for( unsigned k=0; k<NP; k++ )
    SMMLE.add_ctr( mc::BASE_OPT::EQ, dJerrdP[k] );
#else
  SMMLE.set_var( NDAT, P );
  SMMLE.set_dep( NP, P+NDAT, dJerrdP );
  SMMLE.reset_ctr();
#endif
  SMMLE.add_ctr( mc::BASE_OPT::LE, Eerr - Chi2*varYm );
  SMMLE.add_ctr( mc::BASE_OPT::LE, J - JU );

  SMMLE.options.DISPLAY     = 2;
  SMMLE.options.MAXITER     = 0;
  SMMLE.options.NODEMEAS    = mc::SBP<I>::Options::RELMAXLEN;
  SMMLE.options.VARMEAS     = mc::NLCP<I>::Options::USER;
  SMMLE.options.BRANCHVAR   = mc::SBP<I>::Options::RGREL;
  SMMLE.options.BRANCHSEL   = two_phase_branching;
  SMMLE.options.CVTOL       = 8e-3;
  SMMLE.options.FEASTOL     = 1e-6;
  SMMLE.options.BLKDECUSE   = false;
  SMMLE.options.DOMREDMAX   = 10;
  SMMLE.options.DOMREDTHRES = 1e-1;
  SMMLE.options.DOMREDBKOFF = 1e-6;
  SMMLE.options.SCOBCHUSE   = false;//true;
  SMMLE.options.SCOBCHRTOL  = 1e-1;
  SMMLE.options.RELMETH     = mc::NLCP<I>::Options::CHEB;
  SMMLE.options.CMODPROP    = 3;
  SMMLE.options.CMODCUTS    = 2;
  SMMLE.options.CMODDEPS    = 0;

  SMMLE.setup();
  SMMLE.solve( Ip_BND.data() );
  SMMLE.stats.display();

#if defined(SAVE_RESULTS )
  sprintf( fname, "test7b_%2.0f%s_Nm%d.out",100.*conf,"%",NDAT );
  ofile.open( fname, std::ios_base::out );
  SMMLE.output_nodes( ofile );
  ofile.close();
#endif

  { int dum; std::cout << "PAUSED--"; std::cin >> dum; }

  // Set-Membership Regression Region - Enclosure

  mc::NLCP<I> SMMLE2;
  SMMLE2.set_dag( &DAG );
  SMMLE2.set_var( NP, P );
  SMMLE2.add_ctr( mc::BASE_OPT::LE, J - JU );

  SMMLE2.options.DISPLAY     = 2;
  SMMLE2.options.MAXITER     = 0;
  SMMLE2.options.NODEMEAS    = mc::SBP<I>::Options::RELMAXLEN;
  SMMLE2.options.CVTOL       = 5e-4;
  SMMLE2.options.FEASTOL     = 1e-6;
  SMMLE2.options.BLKDECUSE   = false;
  SMMLE2.options.BRANCHVAR   = mc::SBP<I>::Options::RGREL;
  SMMLE2.options.DOMREDMAX   = 10;
  SMMLE2.options.DOMREDTHRES = 1e-1;
  SMMLE2.options.DOMREDBKOFF = 1e-6;
  SMMLE2.options.RELMETH     = mc::NLCP<I>::Options::CHEB;
  SMMLE2.options.CMODPROP    = 3;
  SMMLE2.options.CMODCUTS    = 2;
  SMMLE2.options.CMODDEPS    = 0;

  SMMLE2.setup();
  SMMLE2.solve( Ip );
  SMMLE2.stats.display();

#if defined(SAVE_RESULTS )
  sprintf( fname, "test7a_%2.0f%s_Nm%d.out",100.*conf,"%",NDAT );
  ofile.open( fname, std::ios_base::out );
  SMMLE2.output_nodes( ofile );
  ofile.close();
#endif

  { int dum; std::cout << "PAUSED--"; std::cin >> dum; }

  return 0;
}

