// OCFE_PDE36_march.cpp  ---  FULL PHYSICAL PSA (rung 9) under SOLVE_MARCHING
// ===========================================================================
// Marching adaptation of OCFE_PDE36: the first PHYSICAL multi-component run
// (NO manufactured solution).  A regenerated (clean) bed is fed a binary mixture;
// component 1 (strongly adsorbed) breaks through late, component 2 rolls up, and
// the adsorption exotherm drives a thermal wave relaxed by wall cooling.  Five
// DIFFERENTIAL states c1,c2,q1,q2,T; competitive van't Hoff Langmuir; a smooth
// exponential-ramp inlet feed.
//
//   mass_i : dci/dt + U dci/dz - D d2ci/dz2 + F dqi/dt          = 0   (i=1,2)
//   LDF_i  : dqi/dt - k ( qi*(c1,c2,T) - qi )                   = 0
//   energy : Cp dT/dt + G dT/dz - lam d2T/dz2
//              - F ( dH1 dq1/dt + dH2 dq2/dt ) + hw (T - Tw)    = 0
//   IC (clean bed): c1=c2=q1=q2=0, T=T0.   feed: ci = c0_i (1-e^{-t/tau}).
//
// MARCHING: the clean-bed ICs are posed as distributed z-inputs (c1_ic..q2_ic=0,
// T_ic=T0); their (state,input) transfer pairs are AUTO-DETECTED; all five states
// are differential, so all five transfer (no algebraic reinit here).  The feed is
// symbolic in t, so window-sliding evaluates it at each element's absolute time --
// no input re-seeding needed.
//
// VALIDATION (no exact solution): the MARCHED field is compared point-by-point
// against the MONOLITHIC field on the same model over a (t,z) grid; agreement to
// the solve tolerance means marching reproduces the full physical solution.
// ===========================================================================
#include <iostream>
#include <iomanip>
#include <vector>
#include <array>
#include <cmath>
#include <string>

#ifndef OCFE_OCFESLV_HEADER
#define OCFE_OCFESLV_HEADER "ocfeslv.hpp"
#endif
#include OCFE_OCFESLV_HEADER

using namespace mc;

static double const U_VEL = 1.0, D_ax = 0.1, F_ph = 0.5, k_ldf = 2.0;
static double const qs1 = 1.0, qs2 = 1.0, b01_L = 3.0, b02_L = 1.0;
static double const beta1 = 2.0, beta2 = 1.0, T0_ref = 1.0;
static double const Cp_e = 1.0, G_cv = 1.0, lam = 0.1, dH1 = 1.0, dH2 = 0.5, hw = 0.5, Tw = 1.0;
static double const c0_1 = 0.5, c0_2 = 0.5, tau_in = 0.15, T_end = 5.0;

static inline double c1feed_d( double t ){ return c0_1*( 1.0 - std::exp( -t/tau_in ) ); }
static inline double c2feed_d( double t ){ return c0_2*( 1.0 - std::exp( -t/tau_in ) ); }
static inline double q1star_T0( double c1g, double c2g )
{ return qs1*b01_L*c1g/( 1.0 + b01_L*c1g + b02_L*c2g ); }
static inline double q2star_T0( double c1g, double c2g )
{ return qs2*b02_L*c2g/( 1.0 + b01_L*c1g + b02_L*c2g ); }
// amount of species i fed over [0,T_end] at the inlet (Danckwerts flux U*ci,feed)
static inline double feed_throughput( double c0 )
{ return U_VEL*c0*( T_end - tau_in*( 1.0 - std::exp( -T_end/tau_in ) ) ); }

static int g_pass = 0, g_fail = 0;
static void check_close( char const* name, double got, double want, double tol )
{
  bool ok = std::fabs( got - want ) <= tol;
  ( ok ? g_pass : g_fail )++;
  std::cout << "  " << std::left << std::setw(46) << name
            << std::scientific << std::setprecision(3)
            << " err=" << std::fabs( got - want ) << " tol=" << tol
            << "  " << ( ok ? "PASS" : "FAIL" ) << "\n";
}
static void check_true( char const* name, bool ok )
{
  ( ok ? g_pass : g_fail )++;
  std::cout << "  " << std::left << std::setw(46) << name << "  " << ( ok ? "PASS" : "FAIL" ) << "\n";
}

// Solve the physical model (monolithic or marching) and sample c1,c2,q1,q2,T at a
// fixed (t,z) grid: field[i*(NZs+1)+j] holds the five states at (t=i/NTs*T_end, z=j/NZs).
static bool run_pde36( FFDom::TYPE coltype, size_t n_el, size_t n_nd, bool march,
                       std::vector<std::array<double,5>>& field, std::array<double,4>& funcs,
                       int NTs, int NZs )
{
  std::cout << "\n------------------------------------------------------------\n"
            << "  PSA rung 9 (physical breakthrough)  n_el=" << n_el << " n_nd=" << n_nd
            << "  (" << ( march ? "SOLVE_MARCHING" : "monolithic" ) << ")\n"
            << "------------------------------------------------------------\n";

  FFGraph DAG;
  FFVar t  = DAG.add_var( "t" );
  FFVar z  = DAG.add_var( "z" );
  FFVar c1 = DAG.add_var( "c1(t,z)" );
  FFVar c2 = DAG.add_var( "c2(t,z)" );
  FFVar q1 = DAG.add_var( "q1(t,z)" );
  FFVar q2 = DAG.add_var( "q2(t,z)" );
  FFVar T  = DAG.add_var( "T(t,z)"  );
  FFVar c1_ic = DAG.add_var( "c1_ic(z)" );
  FFVar c2_ic = DAG.add_var( "c2_ic(z)" );
  FFVar q1_ic = DAG.add_var( "q1_ic(z)" );
  FFVar q2_ic = DAG.add_var( "q2_ic(z)" );
  FFVar T_ic  = DAG.add_var( "T_ic(z)"  );

  FFPartial  OpP;
  FFIntegral OpI;

  FFVar c1feed = c0_1*( 1.0 - exp( -t/tau_in ) );   // symbolic feed (window-slide handles it)
  FFVar c2feed = c0_2*( 1.0 - exp( -t/tau_in ) );
  FFVar b1 = b01_L*exp( beta1*( 1.0/T - 1.0/T0_ref ) );
  FFVar b2 = b02_L*exp( beta2*( 1.0/T - 1.0/T0_ref ) );
  FFVar den = 1.0 + b1*c1 + b2*c2;
  FFVar q1star = qs1*b1*c1/den;
  FFVar q2star = qs2*b2*c2/den;

  FFVar CONT1 = OpP( c1, t ) + U_VEL*OpP( c1, z ) - D_ax*OpP( OpP( c1, z ), z ) + F_ph*OpP( q1, t );
  FFVar CONT2 = OpP( c2, t ) + U_VEL*OpP( c2, z ) - D_ax*OpP( OpP( c2, z ), z ) + F_ph*OpP( q2, t );
  FFVar LDF1  = OpP( q1, t ) - k_ldf*( q1star - q1 );
  FFVar LDF2  = OpP( q2, t ) - k_ldf*( q2star - q2 );
  FFVar ENE_T = Cp_e*OpP( T, t ) + G_cv*OpP( T, z ) - lam*OpP( OpP( T, z ), z )
              - F_ph*( dH1*OpP( q1, t ) + dH2*OpP( q2, t ) ) + hw*( T - Tw );

  // clean-bed ICs: input-based (transferable) for marching, direct for monolithic
  FFVar IC_c1 = march ? ( c1 - c1_ic ) : ( c1 );
  FFVar IC_c2 = march ? ( c2 - c2_ic ) : ( c2 );
  FFVar IC_q1 = march ? ( q1 - q1_ic ) : ( q1 );
  FFVar IC_q2 = march ? ( q2 - q2_ic ) : ( q2 );
  FFVar IC_T  = march ? ( T  - T_ic  ) : ( T - T0_ref );

  FFVar BC_L1 = U_VEL*c1 - D_ax*OpP( c1, z ) - U_VEL*c1feed;   // Danckwerts inlet
  FFVar BC_L2 = U_VEL*c2 - D_ax*OpP( c2, z ) - U_VEL*c2feed;
  FFVar BC_LT = G_cv*T  - lam*OpP( T,  z ) - G_cv*T0_ref;      // thermal feed at T0
  FFVar BC_U1 = OpP( c1, z );                                  // zero-gradient outlets
  FFVar BC_U2 = OpP( c2, z );
  FFVar BC_UT = OpP( T,  z );

  // output functionals: species inventory (spatial integral at final time, TERMINAL) and
  // cumulative effluent (evolution/time integral at the outlet, ACCUMULATED).
  FFVar Inv1 = OpI( c1 + F_ph*q1, z );   // bed inventory of species 1
  FFVar Inv2 = OpI( c2 + F_ph*q2, z );
  FFVar Eff1 = OpI( U_VEL*c1, t );        // cumulative effluent of species 1
  FFVar Eff2 = OpI( U_VEL*c2, t );

  OCFESLV oc( &DAG );
  oc.add_domain( t, FFDom( 0., T_end, n_el, coltype, n_nd ) );
  oc.add_domain( z, FFDom( 0., 1.0,   n_el, coltype, n_nd ) );
  oc.add_state ( c1, {t,z} );  oc.add_state ( c2, {t,z} );
  oc.add_state ( q1, {t,z} );  oc.add_state ( q2, {t,z} );
  oc.add_state ( T,  {t,z} );
  // physical (feed-scaled) reference guesses, as in PDE36, to seed the Newton solve
  auto c1g_ref = [&]( OCFESLV::t_Coord const& cr ){ return c1feed_d( cr.at(t) )*( 1.0 - 0.5*cr.at(z) ); };
  auto c2g_ref = [&]( OCFESLV::t_Coord const& cr ){ return c2feed_d( cr.at(t) )*( 1.0 - 0.5*cr.at(z) ); };
  oc.update_ref( c1, c1g_ref );
  oc.update_ref( c2, c2g_ref );
  oc.update_ref( q1, [&]( OCFESLV::t_Coord const& cr ){ return q1star_T0( c1g_ref(cr), c2g_ref(cr) ); } );
  oc.update_ref( q2, [&]( OCFESLV::t_Coord const& cr ){ return q2star_T0( c1g_ref(cr), c2g_ref(cr) ); } );
  oc.update_ref( T,  [&]( OCFESLV::t_Coord const& cr ){
    return T0_ref + F_ph*( dH1*q1star_T0(c1g_ref(cr),c2g_ref(cr)) + dH2*q2star_T0(c1g_ref(cr),c2g_ref(cr)) )/Cp_e; } );

  if( march ){
    oc.add_input ( c1_ic, {z} );  oc.add_input ( c2_ic, {z} );
    oc.add_input ( q1_ic, {z} );  oc.add_input ( q2_ic, {z} );
    oc.add_input ( T_ic,  {z} );
    oc.update_ref( c1_ic, [&]( OCFESLV::t_Coord const& cr ){ return 0.0; } );     // clean bed
    oc.update_ref( c2_ic, [&]( OCFESLV::t_Coord const& cr ){ return 0.0; } );
    oc.update_ref( q1_ic, [&]( OCFESLV::t_Coord const& cr ){ return 0.0; } );
    oc.update_ref( q2_ic, [&]( OCFESLV::t_Coord const& cr ){ return 0.0; } );
    oc.update_ref( T_ic,  [&]( OCFESLV::t_Coord const& cr ){ return T0_ref; } );  // T0
  }

  oc.set_evolution_domain( t );

  int const T_NO_LB = FFDom::ALL - FFDom::LB;
  int const Z_INT   = FFDom::ALL - FFDom::LB - FFDom::UB;
  oc.add_equation( CONT1, {t,z}, {T_NO_LB,   Z_INT},      OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( CONT2, {t,z}, {T_NO_LB,   Z_INT},      OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( LDF1,  {t,z}, {T_NO_LB,   FFDom::ALL}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( LDF2,  {t,z}, {T_NO_LB,   FFDom::ALL}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( ENE_T, {t,z}, {T_NO_LB,   Z_INT},      OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( IC_c1, {t,z}, {FFDom::LB, FFDom::ALL}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INITIAL,  0 ) );
  oc.add_equation( IC_c2, {t,z}, {FFDom::LB, FFDom::ALL}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INITIAL,  0 ) );
  oc.add_equation( IC_q1, {t,z}, {FFDom::LB, FFDom::ALL}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INITIAL,  0 ) );
  oc.add_equation( IC_q2, {t,z}, {FFDom::LB, FFDom::ALL}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INITIAL,  0 ) );
  oc.add_equation( IC_T,  {t,z}, {FFDom::LB, FFDom::ALL}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INITIAL,  0 ) );
  oc.add_equation( BC_L1, {t,z}, {T_NO_LB,   FFDom::LB},  OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );
  oc.add_equation( BC_L2, {t,z}, {T_NO_LB,   FFDom::LB},  OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );
  oc.add_equation( BC_LT, {t,z}, {T_NO_LB,   FFDom::LB},  OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );
  oc.add_equation( BC_U1, {t,z}, {T_NO_LB,   FFDom::UB},  OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );
  oc.add_equation( BC_U2, {t,z}, {T_NO_LB,   FFDom::UB},  OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );
  oc.add_equation( BC_UT, {t,z}, {T_NO_LB,   FFDom::UB},  OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );

  oc.add_output( Inv1, {t}, {T_end} );   // fct 0: inventory of sp1 at t=T_end   (terminal)
  oc.add_output( Inv2, {t}, {T_end} );   // fct 1
  oc.add_output( Eff1, {z}, {1.0}   );   // fct 2: cumulative effluent of sp1|z=1 (accumulated)
  oc.add_output( Eff2, {z}, {1.0}   );   // fct 3

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
  if( march ){
    oc.options.SOLVE.WARMSTART        = OCFESLV::Options::BROADCAST_IC;   // transfers auto-detected
    oc.options.OUTPUT.MARCH_STORE = true;                          // for field sampling
  }
  else
    oc.options.SOLVE.MARCHING = false;

  if( !oc.setup() ){ std::cerr << "  setup() FAILED\n"; ++g_fail; return false; }
  if( march ){
    check_true( "auto-detected 5 transfer pairs", oc.n_marching_transfers() == 5 );
    check_true( "collapsed to one evolution element", oc.n_evolution_elem() == 1 );
  }

  std::vector<double> xv, inp;
  if( !oc.init( xv, inp, nullptr ) ){ std::cerr << "  init() FAILED\n"; ++g_fail; return false; }
  double const* ip = inp.empty() ? nullptr : inp.data();

  OCFESLV::SolveReport const rep = oc.solve( xv.data(), ip, nullptr );
  std::cout << "  solve: converged=" << ( rep.converged ? "yes" : "no" )
            << " iters=" << rep.iterations << "\n";
  check_true( "solve converged", rep.converged );
  if( !rep.converged ) return false;

  // sample the field
  field.assign( size_t(NTs+1)*size_t(NZs+1), { 0.,0.,0.,0.,0. } );
  FFVar const S[5] = { c1, c2, q1, q2, T };
  for( int i = 0; i <= NTs; ++i ){
    double const tv = double(i)/NTs * T_end;
    for( int j = 0; j <= NZs; ++j ){
      double const zv = double(j)/NZs;
      std::array<double,5>& out = field[ size_t(i)*size_t(NZs+1) + size_t(j) ];
      for( int s = 0; s < 5; ++s ){
        OCFESLV::t_Coord pt; pt[z] = zv; pt[t] = tv;
        if( march )
          out[s] = oc.eval_solution( S[s], pt );          // buffer-free; routes through the trajectory
        else
          out[s] = oc.eval_colloc<double>( S[s], pt, xv.data(), ip, nullptr );
      }
    }
  }

  // output functions: marched accumulation vs monolithic eval(fct); then the physical
  // mass balance  Feed_i = Inv_i + Eff_i  (fed = stored at T_end + effluent out).
  std::vector<double> fct;
  if( march ) fct = oc.val_functions();
  else{
    fct.assign( oc.n_colloc_fct(), 0. );
    std::vector<double> eqn( oc.n_colloc_eqn(), 0. );
    oc.eval( eqn.data(), fct.data(), xv.data(), ip, nullptr );
  }
  if( fct.size() >= 4 ){
    funcs = { fct[0], fct[1], fct[2], fct[3] };            // Inv1, Inv2, Eff1, Eff2
    double const feed1 = feed_throughput( c0_1 ), feed2 = feed_throughput( c0_2 );
    double const mb1 = std::fabs( feed1 - ( fct[0] + fct[2] ) )/feed1;
    double const mb2 = std::fabs( feed2 - ( fct[1] + fct[3] ) )/feed2;
    std::cout << std::scientific << std::setprecision(3)
              << "  Inv1=" << fct[0] << " Eff1=" << fct[2] << " Feed1=" << feed1 << "  (mass-bal rel=" << mb1 << ")\n"
              << "  Inv2=" << fct[1] << " Eff2=" << fct[3] << " Feed2=" << feed2 << "  (mass-bal rel=" << mb2 << ")\n";
    check_close( "species-1 mass balance closes (Feed=Inv+Eff)", mb1, 0., 5e-3 );
    check_close( "species-2 mass balance closes (Feed=Inv+Eff)", mb2, 0., 5e-3 );
  }
  else{ funcs = { 0.,0.,0.,0. }; check_true( "4 output functions produced", false ); }
  return true;
}

int main()
{
  std::cout << "================================================================\n"
            << "  FULL PHYSICAL PSA (rung 9) marching test -- marched vs monolithic\n"
            << "  5 differential states c1,c2,q1,q2,T; non-isothermal van't Hoff;\n"
            << "  clean-bed ICs (auto-detected transfers); symbolic ramp feed.\n"
            << "================================================================\n";

  size_t const n_el = 5;             // 5 time elements over [0,5] -> 5 march steps
  int const NTs = 20, NZs = 10;
  char const* sname[5] = { "c1", "c2", "q1", "q2", "T" };

  auto grid_diff = [&]( size_t n_nd, double& allmax ) -> bool {
    std::vector<std::array<double,5>> fmono, fmarch;
    std::array<double,4> gmono = {0,0,0,0}, gmarch = {0,0,0,0};
    bool const okM = run_pde36( FFDom::LGL, n_el, n_nd, /*march=*/false, fmono,  gmono,  NTs, NZs );
    bool const okm = run_pde36( FFDom::LGL, n_el, n_nd, /*march=*/true,  fmarch, gmarch, NTs, NZs );
    if( !( okM && okm && fmono.size() == fmarch.size() ) ){
      check_true( "both solves produced comparable fields", false ); allmax = 1e30; return false;
    }
    double sdiff[5] = { 0.,0.,0.,0.,0. };
    for( size_t p = 0; p < fmono.size(); ++p )
      for( int s = 0; s < 5; ++s )
        sdiff[s] = std::max( sdiff[s], std::fabs( fmono[p][s] - fmarch[p][s] ) );
    allmax = 0.;
    std::cout << "\n  marched-vs-monolithic max |diff| (n_nd=" << n_nd << "):\n"
              << std::scientific << std::setprecision(3);
    for( int s = 0; s < 5; ++s ){
      std::cout << "    " << sname[s] << ": " << sdiff[s] << "\n";
      allmax = std::max( allmax, sdiff[s] );
    }
    // marched-accumulated outputs (Inv terminal, Eff evolution-integral) vs monolithic
    char const* fname[4] = { "Inv1", "Inv2", "Eff1", "Eff2" };
    double fdiff = 0.;
    std::cout << "  outputs (marched):";
    for( int k = 0; k < 4; ++k ){
      fdiff = std::max( fdiff, std::fabs( gmono[k] - gmarch[k] ) );
      std::cout << " " << fname[k] << "=" << gmarch[k];
    }
    std::cout << "  (max|diff vs mono|=" << fdiff << ")\n";
    check_close( "marched outputs ~ monolithic (Inv terminal, Eff accumulated)", fdiff, 0., 3e-3 );
    return true;
  };

  // The marched (purely causal) and monolithic (IC_STRONG continuity) solutions of a
  // PHYSICAL problem agree only to DISCRETIZATION level, and that gap must SHRINK under
  // refinement -- both converge to the same physical trajectory.
  double diff_coarse = 0., diff_fine = 0.;
  bool const okc = grid_diff( 6, diff_coarse );
  bool const okf = grid_diff( 9, diff_fine   );
  if( okc && okf ){
    check_close( "marched ~ monolithic at n_nd=6 (discretization level)", diff_coarse, 0., 3e-3 );
    // When marching and monolithic solve the SAME discretisation (causal time interfaces: CRONOS_EVOLUTION_CAUSAL,
    // 2026-10-02) they agree to round-off at every refinement, and there is nothing left to improve.
    check_true ( "agreement improves under refinement (n_nd 6 -> 9), or is already at round-off",
                 diff_fine < diff_coarse || std::max( diff_coarse, diff_fine ) < 1e-9 );
    std::cout << std::scientific << std::setprecision(3)
              << "  refinement: n_nd=6 diff=" << diff_coarse
              << "  ->  n_nd=9 diff=" << diff_fine
              << "  (ratio " << std::fixed << std::setprecision(1)
              << ( diff_fine>0 ? diff_coarse/diff_fine : 0. ) << "x)\n";
  }

  std::cout << "\n============================================================\n"
            << "  RESULT: " << g_pass << " passed, " << g_fail << " failed  -- "
            << ( g_fail == 0 ? "ALL PASS" : "SOME FAILED" ) << "\n"
            << "============================================================\n";
  return g_fail == 0 ? 0 : 1;
}
