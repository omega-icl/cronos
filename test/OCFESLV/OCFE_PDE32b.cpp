// OCFE_PDE32_march.cpp  ---  PSA rung 7a (MULTI-COMPONENT) via SOLVE_MARCHING
// ===========================================================================
// Adaptation of OCFE_PDE32 (binary gas c1,c2 on a packed bed; adsorbed q1,q2 via
// LDF; total-pressure EOS P=Rg(c1+c2)) to exercise the SOLVE_MARCHING item-1 core:
// MULTI-STATE terminal transfer.
//
//   differential states (carry an IC, transferred):  c1, c2, q1, q2
//   algebraic  state    (no IC, reinit'd per block):  P   (EOS-slaved)
//
// The four ICs are posed as distributed inputs (c1_ic..q2_ic over z) and
// designated with set_marching_transfer; SOLVE_MARCHING is on by default, so the
// four transfers make it available.  The manufactured sources/BC are symbolic in
// t, so window-sliding auto-evaluates them at each element's absolute window --
// no input re-seeding (item 2) is needed here.  BROADCAST_IC is used as the
// DAE-default warm-start so P (and the differential states) reinit near the
// consistent LB manifold at every element transition.
//
// Validation: the MARCHED terminal at t=1 is compared against the manufactured
// exact solution, and a MONOLITHIC solve of the same model is run as a
// cross-check.  Outputs (time-integral purity/recovery) are NOT evaluated here --
// they require the full-grid function machinery (item 3), which marching does not
// build; this driver validates the multi-state STATE march only.
// ===========================================================================
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <fstream>
#include <string>

#ifndef OCFE_OCFESLV_HEADER
#define OCFE_OCFESLV_HEADER "ocfeslv.hpp"
#endif
#include OCFE_OCFESLV_HEADER

using namespace mc;

// transport / kinetics (single shared bed)
static double const U_VEL = 1.0, D_ax = 0.1, F_ph = 0.5, k_ldf = 2.0, Rg = 1.0;
// competitive Langmuir
static double const qs1 = 1.0, qs2 = 1.0, b1_L = 1.0, b2_L = 0.5;
// manufactured coefficients -- component 1
static double const Ac1 = 0.3, Cc1 = 0.2;
static double const Q10 = 0.4, Qz1 = 0.1, Qt1 = 0.3;
// manufactured coefficients -- component 2
static double const C20 = 0.5, Ac2 = 0.2, Cc2 = 0.15;
static double const Q20 = 0.25, Qz2 = 0.08, Qt2 = 0.12;

static inline double c1M_exact( double z, double t ){ return 1.0 + Ac1*(1.0-z)*(1.0-z) + Cc1*t; }
static inline double c2M_exact( double z, double t ){ return C20 + Ac2*(1.0-z)*(1.0-z) + Cc2*t; }
static inline double q1M_exact( double z, double t ){ return Q10 + Qz1*z + Qt1*t; }
static inline double q2M_exact( double z, double t ){ return Q20 + Qz2*z + Qt2*t; }
static inline double PM_exact ( double z, double t ){ return Rg*( c1M_exact(z,t) + c2M_exact(z,t) ); }

static int g_pass = 0, g_fail = 0;
static void check_close( char const* name, double got, double want, double tol )
{
  bool ok = std::fabs( got - want ) <= tol;
  ( ok ? g_pass : g_fail )++;
  std::cout << "  " << std::left << std::setw(42) << name
            << std::scientific << std::setprecision(3)
            << " err=" << std::fabs( got - want ) << " tol=" << tol
            << "  " << ( ok ? "PASS" : "FAIL" ) << "\n";
}
static void check_true( char const* name, bool ok )
{
  ( ok ? g_pass : g_fail )++;
  std::cout << "  " << std::left << std::setw(42) << name << "  " << ( ok ? "PASS" : "FAIL" ) << "\n";
}

// march=false: symbolic ICs, monolithic solve (cross-check).
// march=true : input ICs + transfers, SOLVE_MARCHING block march.
static bool run_psa7a( FFDom::TYPE coltype, std::string const& cname, bool march )
{
  size_t const n_el = 3, n_nd = 6;
  std::cout << "\n------------------------------------------------------------\n"
            << "  PSA rung 7a  " << cname << "  (" << ( march ? "SOLVE_MARCHING" : "monolithic" ) << ")\n"
            << "------------------------------------------------------------\n";

  FFGraph DAG;
  FFVar t  = DAG.add_var( "t" );
  FFVar z  = DAG.add_var( "z" );
  FFVar c1 = DAG.add_var( "c1(t,z)" );
  FFVar c2 = DAG.add_var( "c2(t,z)" );
  FFVar q1 = DAG.add_var( "q1(t,z)" );
  FFVar q2 = DAG.add_var( "q2(t,z)" );
  FFVar P  = DAG.add_var( "P(t,z)"  );
  // IC inputs (distributed over the spatial grid) -- internal to the march
  FFVar c1_ic = DAG.add_var( "c1_ic(z)" );
  FFVar c2_ic = DAG.add_var( "c2_ic(z)" );
  FFVar q1_ic = DAG.add_var( "q1_ic(z)" );
  FFVar q2_ic = DAG.add_var( "q2_ic(z)" );

  FFPartial  OpP;
  FFIntegral OpI;

  // manufactured fields (for symbolic sources)
  FFVar c1Man = 1.0 + Ac1*(1.0-z)*(1.0-z) + Cc1*t;
  FFVar c2Man = C20 + Ac2*(1.0-z)*(1.0-z) + Cc2*t;
  FFVar q1Man = Q10 + Qz1*z + Qt1*t;
  FFVar q2Man = Q20 + Qz2*z + Qt2*t;
  FFVar denMan = 1.0 + b1_L*c1Man + b2_L*c2Man;
  FFVar q1starMan = qs1*b1_L*c1Man/denMan;
  FFVar q2starMan = qs2*b2_L*c2Man/denMan;

  // symbolic sources / inlet fluxes (functions of t,z -> window-slide handles them)
  FFVar s_c1 = Cc1 - 2.0*U_VEL*Ac1*(1.0-z) - 2.0*D_ax*Ac1 + F_ph*Qt1;
  FFVar s_c2 = Cc2 - 2.0*U_VEL*Ac2*(1.0-z) - 2.0*D_ax*Ac2 + F_ph*Qt2;
  FFVar s_q1 = Qt1 - k_ldf*( q1starMan - q1Man );
  FFVar s_q2 = Qt2 - k_ldf*( q2starMan - q2Man );
  FFVar g1_in = U_VEL*( 1.0 + Ac1 + Cc1*t ) + 2.0*D_ax*Ac1;
  FFVar g2_in = U_VEL*( C20 + Ac2 + Cc2*t ) + 2.0*D_ax*Ac2;

  // governing residuals
  FFVar den    = 1.0 + b1_L*c1 + b2_L*c2;
  FFVar q1star = qs1*b1_L*c1/den;
  FFVar q2star = qs2*b2_L*c2/den;
  FFVar CONT1 = OpP( c1, t ) + U_VEL*OpP( c1, z ) - D_ax*OpP( OpP( c1, z ), z ) + F_ph*OpP( q1, t ) - s_c1;
  FFVar CONT2 = OpP( c2, t ) + U_VEL*OpP( c2, z ) - D_ax*OpP( OpP( c2, z ), z ) + F_ph*OpP( q2, t ) - s_c2;
  FFVar LDF1  = OpP( q1, t ) - k_ldf*( q1star - q1 ) - s_q1;
  FFVar LDF2  = OpP( q2, t ) - k_ldf*( q2star - q2 ) - s_q2;
  FFVar EOS   = P - Rg*( c1 + c2 );

  // ICs: input-based for the march (transferable), symbolic for monolithic
  FFVar IC_c1 = march ? ( c1 - c1_ic ) : ( c1 - ( 1.0 + Ac1*(1.0-z)*(1.0-z) ) );
  FFVar IC_c2 = march ? ( c2 - c2_ic ) : ( c2 - ( C20 + Ac2*(1.0-z)*(1.0-z) ) );
  FFVar IC_q1 = march ? ( q1 - q1_ic ) : ( q1 - ( Q10 + Qz1*z ) );
  FFVar IC_q2 = march ? ( q2 - q2_ic ) : ( q2 - ( Q20 + Qz2*z ) );

  FFVar BC_L1 = U_VEL*c1 - D_ax*OpP( c1, z ) - g1_in;
  FFVar BC_L2 = U_VEL*c2 - D_ax*OpP( c2, z ) - g2_in;
  FFVar BC_U1 = OpP( c1, z );
  FFVar BC_U2 = OpP( c2, z );

  // output functionals: one evolution-direction integral (F1in, OpI over t) and several
  // terminal/final-time outputs.
  FFVar PUR1 = c1 / ( c1 + c2 );      // purity of comp 1
  FFVar F1in = OpI( U_VEL*c1, t );    // cumulative inlet feed of comp 1, integrated over t
  FFVar Q1T  = OpI( q1, z );          // bed loading of comp 1, integrated over z
  FFVar C1T  = OpI( c1, z );          // gas inventory of comp 1, integrated over z

  OCFESLV oc( &DAG );
  oc.add_domain( t, FFDom( 0., 1., n_el, coltype, n_nd ) );
  oc.add_domain( z, FFDom( 0., 1., n_el, coltype, n_nd ) );
  oc.add_state ( c1, {t,z} );
  oc.add_state ( c2, {t,z} );
  oc.add_state ( q1, {t,z} );
  oc.add_state ( q2, {t,z} );
  oc.add_state ( P,  {t,z} );
  oc.update_ref( c1, [&]( OCFESLV::t_Coord const& cr ){ return c1M_exact( cr.at(z), cr.at(t) ); } );
  oc.update_ref( c2, [&]( OCFESLV::t_Coord const& cr ){ return c2M_exact( cr.at(z), cr.at(t) ); } );
  oc.update_ref( q1, [&]( OCFESLV::t_Coord const& cr ){ return q1M_exact( cr.at(z), cr.at(t) ); } );
  oc.update_ref( q2, [&]( OCFESLV::t_Coord const& cr ){ return q2M_exact( cr.at(z), cr.at(t) ); } );
  oc.update_ref( P,  [&]( OCFESLV::t_Coord const& cr ){ return PM_exact ( cr.at(z), cr.at(t) ); } );

  if( march ){
    oc.add_input( c1_ic, {z} );
    oc.add_input( c2_ic, {z} );
    oc.add_input( q1_ic, {z} );
    oc.add_input( q2_ic, {z} );
    // element-0 IC = manufactured profile at t=0
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

  // outputs (same order for both modes): 4 terminal point outputs, then F1in (evolution
  // integral, accumulated), then Q1T/C1T (spatial integrals at t=1, terminal).
  oc.add_output( c1,   {z,t}, {1.0, 1.0} );   // c1_out  = c1(1,1)
  oc.add_output( c2,   {z,t}, {1.0, 1.0} );   // c2_out  = c2(1,1)
  oc.add_output( P,    {z,t}, {1.0, 1.0} );   // P_out   = P(1,1)
  oc.add_output( PUR1, {z,t}, {1.0, 1.0} );   // y1_out  = c1/(c1+c2)|out
  oc.add_output( F1in, {z},   {0.0} );        // F1_in   = int_0^1 U c1 dt |z=0   (accumulated)
  oc.add_output( Q1T,  {t},   {1.0} );        // Q1_T    = int q1 dz |t=1
  oc.add_output( C1T,  {t},   {1.0} );        // C1_T    = int c1 dz |t=1

  oc.options.REDUCE.ORDER     = OCFESLV::Options::RED_FULL;
  oc.options.CLASSIFY.MODE         = OCFESLV::Options::CLASS_AUTO;
  oc.options.INTERFACE.TYPE   = OCFESLV::Options::IC_AUTO;
  oc.options.INTERFACE.IMPOSITION  = OCFESLV::Options::IC_STRONG;   // PSA default
  oc.options.INTERFACE.SAT_SIGMA0       = 10.0;
  oc.options.DISPLAY_LEVEL    = march ? 1 : 0;
#if defined(CRONOS__WITH_SPQR)
  oc.options.SOLVE.FACTORIZATION = OCFESLV::Options::SOLVE_SPQR;
#else
  oc.options.SOLVE.FACTORIZATION = OCFESLV::Options::SOLVE_SUPERLU;
#endif

  if( march ){
    // NO explicit set_marching_transfer: the 4 (state, ic_input) pairs are AUTO-DETECTED
    // from the INITIAL equations (c1-c1_ic, c2-c2_ic, q1-q1_ic, q2-q2_ic).  P has no
    // INITIAL equation -> not detected -> reinit'd per block.  SOLVE_MARCHING stays true
    // (default); the auto-detected pairs make it AVAILABLE.
    oc.options.SOLVE.WARMSTART = OCFESLV::Options::BROADCAST_IC;  // DAE reinit warm-start
    oc.options.OUTPUT.MARCH_STORE = true;                  // capture blocks for plotting/inspection
  }
  else{
    oc.options.SOLVE.MARCHING = false;                         // monolithic cross-check
  }

  if( !oc.setup() ){ std::cerr << "  setup() FAILED\n"; ++g_fail; return false; }

  if( march ){
    check_true( "auto-detected 4 transfer pairs", oc.n_marching_transfers() == 4 );
    check_true( "collapsed to one evolution element", oc.n_evolution_elem() == 1 );
    check_true( "march grid recorded (N+1)", oc.march_grid().size() == n_el + 1 );
    check_true( "n_march_steps() == n_el", oc.n_march_steps() == n_el );
  }

  std::vector<double> xv, inp;
  if( !oc.init( xv, inp, nullptr ) ){ std::cerr << "  init() FAILED\n"; ++g_fail; return false; }
  double const* ip = inp.empty() ? nullptr : inp.data();

  // perturb the guess so the solve is non-trivial
  for( size_t i = 0; i < xv.size(); ++i ) xv[i] += 0.05*std::sin( 0.7*double(i) + 0.2 );

  OCFESLV::SolveReport const rep = oc.solve( xv.data(), ip, nullptr );
  std::cout << "  solve: converged=" << ( rep.converged ? "yes" : "no" )
            << " iters=" << rep.iterations
            << " final|r|=" << std::scientific << std::setprecision(3) << rep.final_residual << "\n";
  check_true( "solve converged", rep.converged );
  if( !rep.converged ) return false;

  // error check.  Monolithic: full field.  Marching: terminal at t=1 (the final
  // element's window is [ (n_el-1)/n_el, 1 ], so t=1 is its UB).
  double c1e=0., c2e=0., q1e=0., q2e=0., Pe=0.;
  double const zs_arr[4] = { 0.15, 0.35, 0.65, 0.85 };
  double const ts_arr_mono[4] = { 0.15, 0.35, 0.65, 0.85 };
  double const ts_arr_march[1] = { 1.0 };
  size_t const nts = march ? 1 : 4;
  for( double zs : zs_arr ) for( size_t it = 0; it < nts; ++it ){
    double const ts = march ? ts_arr_march[0] : ts_arr_mono[it];
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
  check_close( "q1 vs exact", q1e, 0., tol );
  check_close( "q2 vs exact", q2e, 0., tol );
  check_close( "P  vs exact (algebraic reinit)", Pe, 0., tol );

  // --- output-function validation: marched accumulation vs analytic (and vs monolithic) ---
  double const c1_out_exact = c1M_exact( 1.0, 1.0 );                       // 1.20
  double const c2_out_exact = c2M_exact( 1.0, 1.0 );                       // 0.65
  double const P_out_exact  = PM_exact ( 1.0, 1.0 );                       // 1.85
  double const y1_out_exact = c1_out_exact / ( c1_out_exact + c2_out_exact );
  double const F1_in_exact  = U_VEL * ( ( 1.0 + Ac1 ) + 0.5*Cc1 );         // 1.40 (int over t)
  double const Q1_T_exact   = ( Q10 + Qt1 ) + 0.5*Qz1;                     // 0.75
  double const C1_T_exact   = ( 1.0 + Cc1 ) + Ac1/3.0;                     // 1.30
  char const*  fname[7] = { "c1_out", "c2_out", "P_out", "y1_out",
                            "F1_in (evolution integral)", "Q1_T", "C1_T" };
  double const fexp [7] = { c1_out_exact, c2_out_exact, P_out_exact, y1_out_exact,
                            F1_in_exact, Q1_T_exact, C1_T_exact };
  std::vector<double> fct;
  if( march ){
    fct = oc.val_functions();         // per-element accumulation (sum F1_in, last-elem others)
  }
  else{
    fct.assign( oc.n_colloc_fct(), 0. );
    std::vector<double> eqn( oc.n_colloc_eqn(), 0. );
    oc.eval( eqn.data(), fct.data(), xv.data(), ip, nullptr );
  }
  if( fct.size() >= 7 ){
    std::cout << "  outputs:";
    for( int i = 0; i < 7; ++i )
      std::cout << " " << fname[i] << "=" << std::fixed << std::setprecision(4) << fct[i];
    std::cout << "\n" << std::scientific;
    for( int i = 0; i < 7; ++i ) check_close( fname[i], fct[i], fexp[i], 1e-6 );
  }
  else check_true( "7 output rows produced", false );

  // --- trajectory store validation (marching only): interior-point interpolation ---
  if( march ){
    check_true( "trajectory has n_el blocks", oc.march_trajectory().size() == n_el );
    // Dense full-horizon sweep of ALL states (t=i/30 hits the element boundaries 1/3,2/3,
    // so continuity across the block hand-offs is exercised), written as a gnuplot data
    // file: whitespace-separated, '#'-commented header, in Z-SCAN order (z outer, t inner)
    // with a blank line after each constant-z scan.  That way each fixed-z slice is one
    // contiguous datablock (index j -> z=j/NZs), so `with lines` works, and splot/pm3d
    // still sees the (t,z) grid.
    std::ofstream dat( "OCFE_PDE32_traj.dat" );
    dat << "# columns: 1=t 2=z  3=c1 4=c1_exact  5=c2 6=c2_exact  7=q1 8=q1_exact"
           "  9=q2 10=q2_exact  11=P 12=P_exact\n";
    dat << "# z-scan order: datablock index j is the slice z = j/10 (0..10)\n";
    dat << std::setprecision(10);
    int const NTs = 30, NZs = 10;
    double terr = 0.;
    for( int j = 0; j <= NZs; ++j ){
      double const zv = double(j)/NZs;
      for( int i = 0; i <= NTs; ++i ){
        double const tv = double(i)/NTs;
        OCFESLV::t_Coord pt; pt[z] = zv; pt[t] = tv;
        double const c1v = oc.eval_solution( c1, pt ), c2v = oc.eval_solution( c2, pt );
        double const q1v = oc.eval_solution( q1, pt ), q2v = oc.eval_solution( q2, pt );
        double const Pv  = oc.eval_solution( P,  pt );
        terr = std::max( { terr,
          std::fabs( c1v - c1M_exact(zv,tv) ), std::fabs( c2v - c2M_exact(zv,tv) ),
          std::fabs( q1v - q1M_exact(zv,tv) ), std::fabs( q2v - q2M_exact(zv,tv) ),
          std::fabs( Pv  - PM_exact (zv,tv) ) } );
        dat << tv <<" "<< zv
            <<"  "<< c1v <<" "<< c1M_exact(zv,tv) <<"  "<< c2v <<" "<< c2M_exact(zv,tv)
            <<"  "<< q1v <<" "<< q1M_exact(zv,tv) <<"  "<< q2v <<" "<< q2M_exact(zv,tv)
            <<"  "<< Pv  <<" "<< PM_exact (zv,tv) <<"\n";
      }
      dat << "\n";   // blank line between z-scans -> contiguous slices + gnuplot grid
    }
    dat.close();
    std::cout << "  wrote OCFE_PDE32_traj.dat (" << (NTs+1)*(NZs+1)
              << " points; gnuplot: load 'OCFE_PDE32_traj.gp')\n";
    check_close( "dense trajectory (all states, incl. block boundaries) vs exact", terr, 0., 1e-8 );
  }

  return g_fail == 0;
}

int main()
{
  std::cout << "================================================================\n"
            << "  PSA rung 7a marching test (item-1 multi-state core)\n"
            << "  binary (c1,c2)+adsorbed (q1,q2) differential; P algebraic (EOS)\n"
            << "  4 differential transfers, P reinit'd per block, BROADCAST_IC\n"
            << "================================================================\n";
  run_psa7a( FFDom::LGL, "LGL", /*march=*/false );   // monolithic cross-check
  run_psa7a( FFDom::LGL, "LGL", /*march=*/true  );   // SOLVE_MARCHING

  std::cout << "\n============================================================\n"
            << "  RESULT: " << g_pass << " passed, " << g_fail << " failed  -- "
            << ( g_fail == 0 ? "ALL PASS" : "SOME FAILED" ) << "\n"
            << "============================================================\n";
  return g_fail == 0 ? 0 : 1;
}
