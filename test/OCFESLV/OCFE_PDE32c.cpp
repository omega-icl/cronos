// OCFE_PDE32feed_march.cpp  ---  PSA rung 7a with a TIME-VARYING INPUT FEED
// ===========================================================================
// Variant of OCFE_PDE32_march that exercises SOLVE_MARCHING item 2: PER-ELEMENT
// INPUT RE-SEEDING.  The two Danckwerts inlet feeds g1_in(t), g2_in(t) are here
// posed as distributed INPUTS over the evolution domain {t} (not symbolic DAG
// expressions), each given a time-varying reference via update_ref.
//
//   monolithic : the feed input is sampled once over the full [0,1] grid.
//   marching   : the evolution domain is collapsed to one element, so the feed
//                input must be RE-SAMPLED from its reference at each window's
//                absolute t-nodes -- exactly item 2.  Without it, every element
//                would reuse element-0's feed values (t in [0,1/3]) and the
//                solution on later windows would be wrong.
//
// The IC inputs (c1_ic..q2_ic) are auto-detected transfer targets (item 1.2) and
// are NOT re-seeded; the feed inputs have references and ARE re-seeded.  The
// marched solution is validated against the manufactured exact -- a pass means
// the time-varying feed advanced correctly window-to-window.
// ===========================================================================
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <string>

#ifndef OCFE_OCFESLV_HEADER
#define OCFE_OCFESLV_HEADER "ocfeslv.hpp"
#endif
#include OCFE_OCFESLV_HEADER

using namespace mc;

static double const U_VEL = 1.0, D_ax = 0.1, F_ph = 0.5, k_ldf = 2.0, Rg = 1.0;
static double const qs1 = 1.0, qs2 = 1.0, b1_L = 1.0, b2_L = 0.5;
static double const Ac1 = 0.3, Cc1 = 0.2;
static double const Q10 = 0.4, Qz1 = 0.1, Qt1 = 0.3;
static double const C20 = 0.5, Ac2 = 0.2, Cc2 = 0.15;
static double const Q20 = 0.25, Qz2 = 0.08, Qt2 = 0.12;

static inline double c1M_exact( double z, double t ){ return 1.0 + Ac1*(1.0-z)*(1.0-z) + Cc1*t; }
static inline double c2M_exact( double z, double t ){ return C20 + Ac2*(1.0-z)*(1.0-z) + Cc2*t; }
static inline double q1M_exact( double z, double t ){ return Q10 + Qz1*z + Qt1*t; }
static inline double q2M_exact( double z, double t ){ return Q20 + Qz2*z + Qt2*t; }
static inline double PM_exact ( double z, double t ){ return Rg*( c1M_exact(z,t) + c2M_exact(z,t) ); }
// time-varying inlet feeds (functions of t), consistent with the Danckwerts BC at z=0
static inline double g1_feed_exact( double t ){ return U_VEL*( 1.0 + Ac1 + Cc1*t ) + 2.0*D_ax*Ac1; }
static inline double g2_feed_exact( double t ){ return U_VEL*( C20 + Ac2 + Cc2*t ) + 2.0*D_ax*Ac2; }

static int g_pass = 0, g_fail = 0;
static void check_close( char const* name, double got, double want, double tol )
{
  bool ok = std::fabs( got - want ) <= tol;
  ( ok ? g_pass : g_fail )++;
  std::cout << "  " << std::left << std::setw(44) << name
            << std::scientific << std::setprecision(3)
            << " err=" << std::fabs( got - want ) << " tol=" << tol
            << "  " << ( ok ? "PASS" : "FAIL" ) << "\n";
}
static void check_true( char const* name, bool ok )
{
  ( ok ? g_pass : g_fail )++;
  std::cout << "  " << std::left << std::setw(44) << name << "  " << ( ok ? "PASS" : "FAIL" ) << "\n";
}

static bool run_psa7a_feed( FFDom::TYPE coltype, std::string const& cname, bool march )
{
  size_t const n_el = 3, n_nd = 6;
  std::cout << "\n------------------------------------------------------------\n"
            << "  PSA rung 7a + input feed  " << cname
            << "  (" << ( march ? "SOLVE_MARCHING" : "monolithic" ) << ")\n"
            << "------------------------------------------------------------\n";

  FFGraph DAG;
  FFVar t  = DAG.add_var( "t" );
  FFVar z  = DAG.add_var( "z" );
  FFVar c1 = DAG.add_var( "c1(t,z)" );
  FFVar c2 = DAG.add_var( "c2(t,z)" );
  FFVar q1 = DAG.add_var( "q1(t,z)" );
  FFVar q2 = DAG.add_var( "q2(t,z)" );
  FFVar P  = DAG.add_var( "P(t,z)"  );
  FFVar c1_ic = DAG.add_var( "c1_ic(z)" );
  FFVar c2_ic = DAG.add_var( "c2_ic(z)" );
  FFVar q1_ic = DAG.add_var( "q1_ic(z)" );
  FFVar q2_ic = DAG.add_var( "q2_ic(z)" );
  // time-varying inlet feeds as distributed inputs over the evolution domain {t}
  FFVar g1_in = DAG.add_var( "g1_in(t)" );
  FFVar g2_in = DAG.add_var( "g2_in(t)" );

  FFPartial  OpP;

  FFVar c1Man = 1.0 + Ac1*(1.0-z)*(1.0-z) + Cc1*t;
  FFVar c2Man = C20 + Ac2*(1.0-z)*(1.0-z) + Cc2*t;
  FFVar q1Man = Q10 + Qz1*z + Qt1*t;
  FFVar q2Man = Q20 + Qz2*z + Qt2*t;
  FFVar denMan = 1.0 + b1_L*c1Man + b2_L*c2Man;
  FFVar q1starMan = qs1*b1_L*c1Man/denMan;
  FFVar q2starMan = qs2*b2_L*c2Man/denMan;
  FFVar s_c1 = Cc1 - 2.0*U_VEL*Ac1*(1.0-z) - 2.0*D_ax*Ac1 + F_ph*Qt1;
  FFVar s_c2 = Cc2 - 2.0*U_VEL*Ac2*(1.0-z) - 2.0*D_ax*Ac2 + F_ph*Qt2;
  FFVar s_q1 = Qt1 - k_ldf*( q1starMan - q1Man );
  FFVar s_q2 = Qt2 - k_ldf*( q2starMan - q2Man );

  FFVar den    = 1.0 + b1_L*c1 + b2_L*c2;
  FFVar q1star = qs1*b1_L*c1/den;
  FFVar q2star = qs2*b2_L*c2/den;
  FFVar CONT1 = OpP( c1, t ) + U_VEL*OpP( c1, z ) - D_ax*OpP( OpP( c1, z ), z ) + F_ph*OpP( q1, t ) - s_c1;
  FFVar CONT2 = OpP( c2, t ) + U_VEL*OpP( c2, z ) - D_ax*OpP( OpP( c2, z ), z ) + F_ph*OpP( q2, t ) - s_c2;
  FFVar LDF1  = OpP( q1, t ) - k_ldf*( q1star - q1 ) - s_q1;
  FFVar LDF2  = OpP( q2, t ) - k_ldf*( q2star - q2 ) - s_q2;
  FFVar EOS   = P - Rg*( c1 + c2 );

  FFVar IC_c1 = march ? ( c1 - c1_ic ) : ( c1 - ( 1.0 + Ac1*(1.0-z)*(1.0-z) ) );
  FFVar IC_c2 = march ? ( c2 - c2_ic ) : ( c2 - ( C20 + Ac2*(1.0-z)*(1.0-z) ) );
  FFVar IC_q1 = march ? ( q1 - q1_ic ) : ( q1 - ( Q10 + Qz1*z ) );
  FFVar IC_q2 = march ? ( q2 - q2_ic ) : ( q2 - ( Q20 + Qz2*z ) );

  // Danckwerts inlet BCs now driven by the INPUT feeds g1_in(t), g2_in(t):
  FFVar BC_L1 = U_VEL*c1 - D_ax*OpP( c1, z ) - g1_in;
  FFVar BC_L2 = U_VEL*c2 - D_ax*OpP( c2, z ) - g2_in;
  FFVar BC_U1 = OpP( c1, z );
  FFVar BC_U2 = OpP( c2, z );

  OCFESLV oc( &DAG );
  oc.add_domain( t, FFDom( 0., 1., n_el, coltype, n_nd ) );
  oc.add_domain( z, FFDom( 0., 1., n_el, coltype, n_nd ) );
  oc.add_state ( c1, {t,z} );  oc.add_state ( c2, {t,z} );
  oc.add_state ( q1, {t,z} );  oc.add_state ( q2, {t,z} );
  oc.add_state ( P,  {t,z} );
  oc.update_ref( c1, [&]( OCFESLV::t_Coord const& cr ){ return c1M_exact( cr.at(z), cr.at(t) ); } );
  oc.update_ref( c2, [&]( OCFESLV::t_Coord const& cr ){ return c2M_exact( cr.at(z), cr.at(t) ); } );
  oc.update_ref( q1, [&]( OCFESLV::t_Coord const& cr ){ return q1M_exact( cr.at(z), cr.at(t) ); } );
  oc.update_ref( q2, [&]( OCFESLV::t_Coord const& cr ){ return q2M_exact( cr.at(z), cr.at(t) ); } );
  oc.update_ref( P,  [&]( OCFESLV::t_Coord const& cr ){ return PM_exact ( cr.at(z), cr.at(t) ); } );

  // time-varying feed inputs over {t}, with references (re-seeded per window under marching)
  oc.add_input ( g1_in, {t} );
  oc.add_input ( g2_in, {t} );
  oc.update_ref( g1_in, [&]( OCFESLV::t_Coord const& cr ){ return g1_feed_exact( cr.at(t) ); } );
  oc.update_ref( g2_in, [&]( OCFESLV::t_Coord const& cr ){ return g2_feed_exact( cr.at(t) ); } );

  if( march ){
    oc.add_input ( c1_ic, {z} );  oc.add_input ( c2_ic, {z} );
    oc.add_input ( q1_ic, {z} );  oc.add_input ( q2_ic, {z} );
    oc.update_ref( c1_ic, [&]( OCFESLV::t_Coord const& cr ){ return c1M_exact( cr.at(z), 0.0 ); } );
    oc.update_ref( c2_ic, [&]( OCFESLV::t_Coord const& cr ){ return c2M_exact( cr.at(z), 0.0 ); } );
    oc.update_ref( q1_ic, [&]( OCFESLV::t_Coord const& cr ){ return q1M_exact( cr.at(z), 0.0 ); } );
    oc.update_ref( q2_ic, [&]( OCFESLV::t_Coord const& cr ){ return q2M_exact( cr.at(z), 0.0 ); } );
  }

  oc.set_evolution_domain( t );

  int const T_NO_LB = FFDom::ALL - FFDom::LB;
  int const Z_INT   = FFDom::ALL - FFDom::LB - FFDom::UB;
  oc.add_equation( CONT1, {t,z}, {T_NO_LB,   Z_INT},      OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( CONT2, {t,z}, {T_NO_LB,   Z_INT},      OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( LDF1,  {t,z}, {T_NO_LB,   FFDom::ALL}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( LDF2,  {t,z}, {T_NO_LB,   FFDom::ALL}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( EOS,   {t,z}, {FFDom::ALL,FFDom::ALL}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( IC_c1, {t,z}, {FFDom::LB, FFDom::ALL}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INITIAL,  0 ) );
  oc.add_equation( IC_c2, {t,z}, {FFDom::LB, FFDom::ALL}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INITIAL,  0 ) );
  oc.add_equation( IC_q1, {t,z}, {FFDom::LB, FFDom::ALL}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INITIAL,  0 ) );
  oc.add_equation( IC_q2, {t,z}, {FFDom::LB, FFDom::ALL}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INITIAL,  0 ) );
  oc.add_equation( BC_L1, {t,z}, {T_NO_LB,   FFDom::LB},  OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );
  oc.add_equation( BC_L2, {t,z}, {T_NO_LB,   FFDom::LB},  OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );
  oc.add_equation( BC_U1, {t,z}, {T_NO_LB,   FFDom::UB},  OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );
  oc.add_equation( BC_U2, {t,z}, {T_NO_LB,   FFDom::UB},  OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );

  oc.options.REDUCE.ORDER     = OCFESLV::Options::RED_FULL;
  oc.options.CLASSIFY.MODE         = OCFESLV::Options::CLASS_AUTO;
  oc.options.INTERFACE.TYPE   = OCFESLV::Options::IC_AUTO;
  oc.options.INTERFACE.IMPOSITION  = OCFESLV::Options::IC_STRONG;
  oc.options.INTERFACE.SAT_SIGMA0       = 10.0;
  oc.options.DISPLAY_LEVEL    = march ? 1 : 0;
#if defined(CRONOS__WITH_SPQR)
  oc.options.SOLVE.FACTORIZATION = OCFESLV::Options::SOLVE_SPQR;
#else
  oc.options.SOLVE.FACTORIZATION = OCFESLV::Options::SOLVE_SUPERLU;
#endif
  if( march )
    oc.options.SOLVE.WARMSTART = OCFESLV::Options::BROADCAST_IC;   // transfers auto-detected
  else
    oc.options.SOLVE.MARCHING  = false;

  if( !oc.setup() ){ std::cerr << "  setup() FAILED\n"; ++g_fail; return false; }
  if( march ){
    check_true( "auto-detected 4 transfer pairs", oc.n_marching_transfers() == 4 );
    check_true( "collapsed to one evolution element", oc.n_evolution_elem() == 1 );
  }

  std::vector<double> xv, inp;
  if( !oc.init( xv, inp, nullptr ) ){ std::cerr << "  init() FAILED\n"; ++g_fail; return false; }
  double const* ip = inp.empty() ? nullptr : inp.data();
  for( size_t i = 0; i < xv.size(); ++i ) xv[i] += 0.05*std::sin( 0.7*double(i) + 0.2 );

  OCFESLV::SolveReport const rep = oc.solve( xv.data(), ip, nullptr );
  std::cout << "  solve: converged=" << ( rep.converged ? "yes" : "no" )
            << " iters=" << rep.iterations << "\n";
  check_true( "solve converged", rep.converged );
  if( !rep.converged ) return false;

  // terminal(t=1) [marching] or full field [monolithic] vs exact
  double c1e=0., c2e=0., q1e=0., q2e=0., Pe=0.;
  double const zs_arr[4] = { 0.15, 0.35, 0.65, 0.85 };
  size_t const nts = march ? 1 : 4;
  double const ts_mono[4] = { 0.15, 0.35, 0.65, 0.85 };
  for( double zs : zs_arr ) for( size_t it = 0; it < nts; ++it ){
    double const ts = march ? 1.0 : ts_mono[it];
    OCFESLV::t_Coord pt; pt[z]=zs; pt[t]=ts;
    c1e = std::max( c1e, std::fabs( oc.eval_colloc<double>( c1, pt, xv.data(), ip, nullptr ) - c1M_exact(zs,ts) ) );
    c2e = std::max( c2e, std::fabs( oc.eval_colloc<double>( c2, pt, xv.data(), ip, nullptr ) - c2M_exact(zs,ts) ) );
    q1e = std::max( q1e, std::fabs( oc.eval_colloc<double>( q1, pt, xv.data(), ip, nullptr ) - q1M_exact(zs,ts) ) );
    q2e = std::max( q2e, std::fabs( oc.eval_colloc<double>( q2, pt, xv.data(), ip, nullptr ) - q2M_exact(zs,ts) ) );
    Pe  = std::max( Pe,  std::fabs( oc.eval_colloc<double>( P,  pt, xv.data(), ip, nullptr ) - PM_exact (zs,ts) ) );
  }
  std::cout << std::scientific << std::setprecision(3)
            << "  " << ( march ? "terminal(t=1)" : "field" ) << " err  c1=" << c1e << " c2=" << c2e
            << " q1=" << q1e << " q2=" << q2e << " P=" << Pe << "\n";
  double const tol = 1e-8;
  check_close( "c1 vs exact", c1e, 0., tol );
  check_close( "c2 vs exact", c2e, 0., tol );
  check_close( "q1 vs exact (feed-independent)", q1e, 0., tol );
  check_close( "q2 vs exact", q2e, 0., tol );
  check_close( "P  vs exact", Pe, 0., tol );
  return g_fail == 0;
}

int main()
{
  std::cout << "================================================================\n"
            << "  PSA rung 7a + TIME-VARYING INPUT FEED  (item-2 re-seeding test)\n"
            << "  g1_in(t), g2_in(t) are distributed inputs over {t}; marching\n"
            << "  must re-sample them from update_ref at each window's t-nodes.\n"
            << "================================================================\n";
  run_psa7a_feed( FFDom::LGL, "LGL", /*march=*/false );   // monolithic (feed sampled once)
  run_psa7a_feed( FFDom::LGL, "LGL", /*march=*/true  );   // SOLVE_MARCHING (feed re-seeded/element)

  std::cout << "\n============================================================\n"
            << "  RESULT: " << g_pass << " passed, " << g_fail << " failed  -- "
            << ( g_fail == 0 ? "ALL PASS" : "SOME FAILED" ) << "\n"
            << "============================================================\n";
  return g_fail == 0 ? 0 : 1;
}
