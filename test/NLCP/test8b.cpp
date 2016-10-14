#define SAVE_RESULTS    // whether or not to save results to file
#undef  USE_PROFIL	// specify to use PROFIL for interval arithmetic
#undef  USE_FILIB	// specify to use FILIB++ for interval arithmetic
#undef  DEBUG		// whether to output debug information
#define USE_DEPS	// whether to use dependents
#define USE_RRLT	// whether to use reduced RLT
#define MC__USE_CPLEX
#undef  MC__NLCP_SHOW_CHEBSCORES
#undef  MC__SBB_SHOW__SCOREBRANCHING
////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include "nlcp.hpp"
#include "rltred.hpp"

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
typedef mc::CModel<I> CMI;

//No. stages, feed location, No. components 
//const unsigned nc = 3, ns = 3, pf = 1;
const unsigned nc = 3, ns = 5, pf = 3;
//const unsigned nc = 3, ns = 7, pf = 4;
//const unsigned nc = 3, ns = 10, pf = 5;

////////////////////////////////////////////////////////////////////////////////

std::set<unsigned> branching_strategy
( const typename mc::NLCP<I>::NODE*pNode )
{
  // Variable selection rule based on Baharev et al [AIChE Journal 2011, Vol. 57, No. 6]
  // The widest interval corresponding to the bulk composition is bisected exclusively
  // as long as its width is above 'thres'
  std::set<unsigned> sset;
  double thres = 2e-2, maxw = 0.;
  for( unsigned i=1; i<=nc*ns; i++ ){
    //std::cout << pNode->P(i) << std::endl;
    if( pNode->depend().find(i) == pNode->depend().end()
     && mc::Op<CVI>::diam(pNode->P(i)) > maxw )
      maxw = mc::Op<CVI>::diam(pNode->P(i));
  }
  //std::cout << std::setprecision(3) << std::scientific << "  " << maxw;
  if( maxw > thres ) //pNode->index() < 1000 )
    for( unsigned i=1; i<=nc*ns; i++ ) sset.insert( i );
  else
    for( unsigned i=0; i<=ns*(nc*2+1); i++ ) sset.insert( i );
  return sset;
}
/*
////////////////////////////////////////////////////////////////////////////////
std::set<unsigned> branching_strategy
( const typename mc::NLCP<I>::NODE*pNode )
{
  // Variable selection rule based on Baharev et al [AIChE Journal 2011, Vol. 57, No. 6]
  // The widest interval corresponding to the bulk composition is bisected exclusively
  // as long as its width is above 'thres'
  std::set<unsigned> sset;
  double thres = 2e-2, maxw = 0.;
  for( unsigned i=nc*ns+1; i<=2*nc*ns; i++ ){
    //std::cout << pNode->P(i) << std::endl;
    if( pNode->depend().find(i) == pNode->depend().end()
     && mc::Op<CVI>::diam(pNode->P(i)) > maxw )
      maxw = mc::Op<CVI>::diam(pNode->P(i));
  }
  std::cout << std::setprecision(3) << std::scientific << "  " << maxw;
  if( maxw > thres ) //pNode->index() < 1000 )
    for( unsigned i=nc*ns+1; i<=2*nc*ns; i++ ) sset.insert( i );
  else
    for( unsigned i=0; i<=ns*(nc*2+1); i++ ) sset.insert( i );
  return sset;
}
*/
////////////////////////////////////////////////////////////////////////////////

int main() 
{
  // Model Parameters
  double Patm = 1e5;				//Atmospheric pressure
  double V    = 1.38;				//Vapour flow rate

  /*
  	let a[1] := 23.4832; 	let a[2] := 20.5110; 	let a[3] := 20.9064;
	let b[1] := -3634.01;	let b[2] := -2664.30;	let b[3] := -3096.52;
	let c[1] := -33.768;	let c[2] := -79.483;	let c[3] := -53.668;
	
	let r[1,1] := 0.0;	let r[1,2] :=  0.7411;	let r[1,3] :=  0.9645;	
	let r[2,1] := -1.0250;	let r[2,2] := 0.0;	let r[2,3] := -1.4350;
	let r[3,1] := -0.9645;	let r[3,2] :=  2.7470;	let r[3,3] := 0.0;

	let s[1,1] := 0.0;	let s[1,2] := -477.00;	let s[1,3] := -903.1024;
	let s[2,1] :=  72.78;	let s[2,2] := 0.0;	let s[2,3] :=  768.20;
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

  // Equilibrium Model
  mc::FFGraph DAG;
  const unsigned NF = ns*(nc*3+1), NP = NF+1;
  mc::FFVar P[NP], F[NF];
  for( unsigned int i=0; i<NP; i++ ) P[i].set( &DAG );

  // Equilibrium Parameters
  mc::FFVar  D  = *P;
  mc::FFVar *X  = P+1;
  mc::FFVar *K  = X+nc*ns;
  mc::FFVar *G  = K+nc*ns;
  mc::FFVar *Te = G+nc*ns;
  //for( unsigned i = 0; i <    NP ; ++i  ) std::cout << "P[" << i <<"]=" << P[i] << "\n";
  //std::cout << "D=" << D << "\n";
  //for( unsigned i = 0; i < nc*ns ; ++i  ) std::cout << "X[" << i <<"]=" << X[i] << "\n";
  //for( unsigned i = 0; i < nc*ns ; ++i  ) std::cout << "K[" << i <<"]=" << K[i] << "\n";
  //for( unsigned i = 0; i <    ns ; ++i  ) std::cout << "Te[" << i <<"]=" << Te[i] << "\n";

  // Auxiliary Equations 
  mc::FFVar B, L[ns+1], p[nc*ns], Lambda[nc*nc*ns], sum_xLambda[nc*ns], sum;
  // param B := F[N_F] - D;
  B = FF[pf] - D;
  /*
  	param L{j in 0..N} = V - D + sum{k in 0..j} F[k];
  */
  for( unsigned j = 0; j < ns+1 ; ++j ){
    sumd = 0.0;
    for( unsigned k = 0; k <= j; ++k )
      sumd += FF[k];
    L[j] = V - D + sumd ;	
  }
  /*
  	var p{i in 1..C, j in 1..N} = exp(a[i]+b[i]/(T[j]+c[i]));
  */
  for( unsigned i = 0; i < nc; ++i  )
    for( unsigned j = 0; j < ns; ++j )
      p[i*ns+j] = exp( a[i] + b[i] / (Te[j]+c[i]) );  	
  /*
  	var Lambda{i1 in 1..C, i2 in 1..C, j in 1..N} = exp(r[i1,i2]+s[i1,i2]*rcp_T[j]);
  */
  for( unsigned i1 = 0; i1 < nc; ++i1 )
    for( unsigned i2 = 0; i2 < nc; ++i2 )
      for( unsigned j = 0; j < ns; ++j )
        Lambda[nc*ns*i1+ns*i2+j] = exp(r[nc*i1+i2] + s[nc*i1+i2]/Te[j]);
  /*
  	var sum_xLambda{i in 1..C, j in 1..N} = sum{i1 in 1..C} (x[i1,j]*Lambda[i,i1,j]);
  */
  for( unsigned i = 0; i < nc; ++i  )
    for( unsigned j = 0; j < ns; ++j ){
      sum = 0.;
      for( unsigned i1 = 0; i1 < nc; ++i1 )
        sum += X[i1*ns+j] * Lambda[i*nc*ns+i1*ns+j];	
      sum_xLambda[i*ns+j] = sum;
    }

  // Equilibrium Conditions
  unsigned ieq=0;
  /*
  	ACTIVITY EQUATIONS
  	E_aux_G{j in 1..N, i in 1..C}: 	G[i,j] - exp( -log(sum_xLambda[i,j]) + 1.0 - (sum{i2 in 1..C} (x[i2,j]*Lambda[i2,i,j]*rcp_sum_xLambda[i2,j])) );
  */
  for( unsigned i = 0; i < nc; ++i  )
    for( unsigned j = 0; j < ns; ++j ){
      sum = 0.;
      for( unsigned i2 = 0; i2 < nc; ++i2 )
        sum += X[i2*ns+j] * Lambda[nc*ns*i2+ns*i+j] / sum_xLambda[i2*ns+j];
      //F[ ieq++ ] = G[i*ns+j] - exp( 1.0 - sum - log( sum_xLambda[i*ns+j] ) );
      F[ ieq++ ] = G[i*ns+j] * sum_xLambda[i*ns+j] - exp( 1.0 - sum );
      //F[ ieq++ ] = 1. - log( G[i*ns+j] * sum_xLambda[i*ns+j] ) + sum;
    }
  /*
  	AUXILIARY EQUATIONS
  	E_aux_K{j in 1..N, i in 1..C}: 	K[i,j] - gamma[i,j]*(p[i,j]/P) = 0.0;  
  */
  for( unsigned j = 0; j < ns; ++j  ){
    for( unsigned i = 0; i < nc; ++i ){
      F[ ieq++ ] = K[i*ns+j] - G[i*ns+j] * p[i*ns+j] / Patm;
      //std::cout << "E_aux_K->F[" << j*nc+i <<"]=" << F[j*nc+i] << std::endl;
    }
  }
  /*
  	MATERIAL BALANCES
  	M_tot{i in 1..C}: D*(K[i,1]*x[i,1]) + B*x[i,N] - f[i,N_F] = 0.0;
  */
  for( unsigned i = 0; i < nc; ++i  ){
    F[ ieq++ ] = D*(K[i*ns+0]*X[i*ns+0]) + B*X[i*ns+(ns-1)] - f[i*ns+pf];
    //std::cout << "M_tot->F[" << i+(ns*nc) <<"]=" << F[i+(ns*nc)] << std::endl;
  } 
  /* 
  	NOTE THE UNUSUAL FORMULATION
	M_eq{j in 1..N-1, i in 1..C}:
	L[j]*x[i,j] + sum{i1 in j+1..N} f[i,i1] - B*x[i,N] - V*(K[i,j+1]*x[i,j+1]) = 0.0;
  */
  for( unsigned j = 0; j < ns-1; ++j  ){
    for( unsigned i = 0; i < nc; ++i ){
      sum = 0.; 
      for( unsigned i2 = j+1; i2 < ns; ++i2 )
        sum += f[ i*ns+i2 ];
      F[ ieq++ ] = L[j]*X[i*ns+j] + sum - B*X[i*ns+(ns-1)] - V*( K[i*ns+j+1]*X[i*ns+j+1] );
      //std::cout << "M_eq->F[" << j*nc+i+(ns*nc)+nc <<"]=" << F[j*nc+i+(ns*nc)+nc] << std::endl;
    } 
  } 
  /*
  	SUMMATION EQUATIONS
	S_x_eq{j in 1..N}: 	sum{i in 1..C} x[i,j] - 1.0 = 0.0;
  */
  for( unsigned j = 0; j < ns; ++j  ){
    sum = -1.;
    for( unsigned i = 0; i < nc; ++i )
      sum += X[i*ns+j];
    F[ ieq++ ] = sum;
    //std::cout << "S_x_eq->F[" << j+(ns*nc)+nc+(ns-1)*nc <<"]=" << F[j+(ns*nc)+nc+(ns-1)*nc] << std::endl;
  }
 
  // Variable bounds
  I Ip[NP], &ID = *Ip, *IX = Ip+1, *IK = IX+nc*ns, *IG = IK+nc*ns, *ITe = IG+nc*ns;
  //ID = I( 0.455, 0.455 );
  //ID = I( 0.45, 0.46 );
  ID = I( 0.44, 0.47 );
  for( unsigned j=0; j<ns; ++j ){
    IX[j]  = I( 0.   , 0.9998 ); IX[ns+j] = I( 0.0001, 0.9999 ); IX[2*ns+j] = I( 0.0001, 0.9999 );
    IK[j]  = I( 0.97 , 40.52  ); IK[ns+j] = I( 0.2445, 1.317  ); IK[2*ns+j] = I( 0.2745, 1.975  );
    IG[j]  = I( 0.01 , 100.   ); IG[ns+j] = I( 0.01  , 100.   ); IG[2*ns+j] = I( 0.01  , 100.   );
    ITe[j] = I( 336.3, 383.4  );
  }
  std::cout << DAG;

#ifdef USE_RRLT
  // Test for reduction constraints
  mc::RLTRED RRLT( &DAG );
  RRLT.options.LEVEL = mc::RLTRED::Options::FULLSEQ;
  RRLT.options.DISPLAY = 1;
  RRLT.search( NF, F );
  std::vector<const mc::FFVar*> RED = RRLT.constraints();  
#endif

/*
  // Function bounding
  //CMI CMenv( NP, 3, true );
  //CMenv.options.BOUNDER_TYPE = CMI::Options::NAIVE;
  //CVI CVF, CVp[NP];
  typedef mc::SCVar<I> SCVI;
  typedef mc::SCModel<I> SCMI;
  SCMI SCMenv( 5, NP );
  SCVI SCVF, SCVp[NP];
  SCMenv.options.MIXED_IA = true;
  //CMenv.options.MIXED_IA = true;
  for( unsigned i=0; i<NP; i++ ){
    //CVp[i].set( &CMenv, i, Ip[i] );
    SCVp[i].set( &SCMenv, i, Ip[i] );
  }
  for( unsigned i=0; i<NF; i++ ){
    //DAG.output( DAG.subgraph( 1, F+i ) );
    //{int dum; std::cout << "PAUSED";std::cin >> dum;}
    //try{
    //  DAG.eval( 1, F+i, &CVF, NP, P, CVp );
    //  std::cout << CVF;
    //}
    //catch(...){
    //  std::cout << "Error in CVF(" << i << ")" << std::endl;
    //}
    try{
      DAG.eval( 1, F+i, &SCVF, NP, P, SCVp );
      std::cout << SCVF;
    }
    catch(...){
      std::cout << "Error in SCVF(" << i << ")" << std::endl;
    }
    {int dum; std::cout << "PAUSED";std::cin >> dum;}
  }

  // Relaxation
  mc::NLGO<I> RELAX;
  RELAX.set_dag( &DAG );
  RELAX.set_var( NP, P );
  for( unsigned i=0; i<NF; i++ )
    RELAX.add_ctr( mc::BASE_OPT::EQ, F[i] );
  RELAX.set_obj( mc::BASE_NLP::MIN, D );//gamma[0] );
  RELAX.options.RELMETH     = mc::NLGO<I>::Options::CHEB;
  RELAX.options.CMODPROP    = 3; //2
  //RELAX.options.CMODSPAR    = true;
  //RELAX.options.LPPRESOLVE  = false;
  RELAX.options.MIPFILE     = "test8.lp";
  //std::cout << RELAX;
  RELAX.setup();
  try{ 
    int status = RELAX.relax( Ip );
    //std::cout << RELAX.contract( Ip, 0, false ) << std::endl;
    std::cout << "  f* = " << RELAX.get_objective() << std::endl;
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
  return 0;
*/
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
#ifdef USE_RRLT
  for( auto it=RED.begin(); it!=RED.end(); ++it )
    CP.add_ctr( mc::BASE_NLP::EQ, **it, true );
#endif
  CP.options.MIPFILE     = "";//"test8.lp"
  CP.options.DISPLAY     = 2;
  CP.options.MAXITER     = 0;
  CP.options.NODEMEAS    = mc::SetInv<CVI>::Options::MEANWIDTH;
  CP.options.CVATOL      = 5e-5; // 1e-5
  CP.options.CVRTOL      = 5e-5; // 1e-5
  CP.options.BRANCHVAR   = mc::SetInv<CVI>::Options::RGREL;
  //CP.options.BRANCHSEL   = branching_strategy;
  CP.options.SCOBCHUSE   = true;
  CP.options.SCOBCHRTOL  = 1e-1;
  CP.options.STGBCHDEPTH = 0;//10000;
  CP.options.STGBCHDRMAX = 2;
  CP.options.STGBCHRTOL  = 1e-2;
  CP.options.DOMREDMAX   = 20;
  CP.options.DOMREDTHRES = 5e-2; //2e-2
  CP.options.DOMREDBKOFF = 1e-5; //2e-2
  //CP.options.CTRBACKOFF  = 1e-4; //1e-6
  //CP.options.LPALGO      = -1;
  //CP.options.LPPRESOLVE  = -1;
  CP.options.RELMETH     = mc::NLCP<I>::Options::CHEB;//DRL//HYBRID;;
  CP.options.CMODDMAX    = 1e6; //2
  CP.options.CMODPROP    = 2; //2
  CP.options.CMODCUTS    = 2; //2
  CP.options.CMREDORD    = 5; //3
  CP.options.CMREDTHRES  = 1e-4;//1e-3
  CP.options.CMREDWARMS  = false;
  CP.options.MAXCPU      = 1e6 ;
  CP.options.CMODEL.MIXED_IA = true ;
  std::cout << CP;

try{  
  CP.setup();
  CP.solve( Ip );
  std::cout << std::fixed << std::setprecision(1);
  std::cout << "POLIMG: " << CP.stats.tPOLIMG << " sec\n";
  std::cout << "LPSET:  " << CP.stats.tLPSET << " sec\n";
  std::cout << "LPSOL:  " << CP.stats.tLPSOL << " sec, " << CP.stats.nLPSOL << " problems\n";
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
#if defined(SAVE_RESULTS )
  std::ofstream ores;
  ores.open( "test8.out", std::ios_base::out );
  CP.output_nodes( ores, 7 );
  ores.close();
#ifdef USE_DEPS
  ores.open( "test8b.out", std::ios_base::out );
  CP.output_nodes( ores, 20, 7 );
  ores.close();
#endif
#endif

  return 0;
}
