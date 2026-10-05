// OCFE_FFIPDAERed_PSA.cpp  ---  reduced-space FFOCFESLV on the PSA oracle (rung-9)
// ===========================================================================
// Validates the reduced-space external DAG operation FFOCFESLV (ffocfe.hpp),
// the input-only sibling of FFOCFERES, on the PSA marching model of OCFE_PDE36c.
//
// FFOCFESLV exposes ONLY the control inputs -> output functions map; the state
// collocation system is solved INTERNALLY (marching).  Controls: b01 (lumped
// isotherm affinity) + T_ic (distributed initial temperature).  Outputs: the two
// functions Inv1 (fct 0) and Eff1 (fct 1).
//
// Battery:
//   T1  value via direct FFOCFESLV::eval<double>          vs direct oc.solve + marched_functions
//   T2  value via DAG round-trip rdag.eval(...)            vs same reference (exercises the copy ctor path)
//   T3  forward Jacobian via FFOCFESLV::eval<fadbad::F>   seeded identity over the ncd controls
//        3v  primal values carried by the F-sweep          vs reference
//        3a  b01 column                                    vs analytic solve_marching_fsens(b01,{1})   [tight]
//        3b  sum of T_ic columns (uniform IC shift)        vs analytic solve_marching_fsens(T_ic,ones) [tight]
//        3c  b01 column                                    vs central-difference FD through the op      [loose]
// The tight checks (3a/3b) tie the FFOp forward mode to the already-validated
// forward-sensitivity march; the FD check (3c) is an independent sanity bound.
// SHALLOW policy: the DAG shares the caller's registered OCFESLV (registry intact).
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

#include "ffocfe.hpp"

using namespace mc;

static double const U_VEL = 1.0, D_ax = 0.1, F_ph = 0.5, k_ldf = 2.0;
static double const qs1 = 1.0, qs2 = 1.0, b02_L = 1.0;
static double const beta1 = 2.0, beta2 = 1.0, T0_ref = 1.0;
static double const Cp_e = 1.0, G_cv = 1.0, lam = 0.1, dH1 = 1.0, dH2 = 0.5, hw = 0.5, Tw = 1.0;
static double const c0_1 = 0.5, c0_2 = 0.5, tau_in = 0.15, T_end = 5.0;
// --- inlet composition STEP (manufactured discontinuity) ---
static double const t_step = 2.0, eps_step = 0.05;      // step time and (small) smoothing width
static double const cA1 = 0.2, cB1 = 0.5;               // c1 feed: cA1 -> cB1 at t_step
static double const cA2 = 0.5, cB2 = 0.2;               // c2 feed: cA2 -> cB2 (composition swap)
static double const b01_nom = 3.0;

static inline double c1feed_d( double t ){ return cA1 + ( cB1 - cA1 )*0.5*( 1.0 + std::tanh( ( t - t_step )/eps_step ) ); }
static inline double c2feed_d( double t ){ return cA2 + ( cB2 - cA2 )*0.5*( 1.0 + std::tanh( ( t - t_step )/eps_step ) ); }
static inline double q1star_g( double c1g, double c2g, double b01 )
{ return qs1*b01*c1g/( 1.0 + b01*c1g + b02_L*c2g ); }
static inline double q2star_g( double c1g, double c2g, double b01 )
{ return qs2*b02_L*c2g/( 1.0 + b01*c1g + b02_L*c2g ); }

static int g_pass = 0, g_fail = 0;
static void check_close( char const* name, double got, double want, double tol )
{
  bool ok = std::fabs( got - want ) <= tol;
  ( ok ? g_pass : g_fail )++;
  std::cout << "  " << std::left << std::setw(46) << name
            << std::scientific << std::setprecision(3)
            << " |got-want|=" << std::fabs( got - want ) << " tol=" << tol
            << "  " << ( ok ? "PASS" : "FAIL" ) << "\n";
}
static void check_true( char const* name, bool ok )
{
  ( ok ? g_pass : g_fail )++;
  std::cout << "  " << std::left << std::setw(46) << name << "  " << ( ok ? "PASS" : "FAIL" ) << "\n";
}

int main()
{
  std::cout << "================================================================\n"
            << "  STEP-INLET graded-march PSA: manufactured inlet-composition discontinuity\n"
            << "  controls: b01 (lumped) + T_ic (distributed);  outputs: Inv1, Eff1\n"
            << "================================================================\n";

  size_t const n_el = 5, n_nd = 6;
  FFDom::TYPE const coltype = FFDom::LGL;

  // ---------------- build the PSA marching model (== OCFE_PDE36c, march=true) ----------------
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
  FFVar b01   = DAG.add_var( "b01" );

  FFPartial  OpP;
  FFIntegral OpI;

  FFVar c1feed = cA1 + ( cB1 - cA1 )*0.5*( 1.0 + tanh( ( t - t_step )/eps_step ) );
  FFVar c2feed = cA2 + ( cB2 - cA2 )*0.5*( 1.0 + tanh( ( t - t_step )/eps_step ) );
  FFVar b1 = b01  *exp( beta1*( 1.0/T - 1.0/T0_ref ) );
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

  FFVar IC_c1 = c1 - c1_ic;
  FFVar IC_c2 = c2 - c2_ic;
  FFVar IC_q1 = q1 - q1_ic;
  FFVar IC_q2 = q2 - q2_ic;
  FFVar IC_T  = T  - T_ic;

  FFVar BC_L1 = U_VEL*c1 - D_ax*OpP( c1, z ) - U_VEL*c1feed;
  FFVar BC_L2 = U_VEL*c2 - D_ax*OpP( c2, z ) - U_VEL*c2feed;
  FFVar BC_LT = G_cv*T  - lam*OpP( T,  z ) - G_cv*T0_ref;
  FFVar BC_U1 = OpP( c1, z );
  FFVar BC_U2 = OpP( c2, z );
  FFVar BC_UT = OpP( T,  z );

  FFVar Inv1 = OpI( c1 + F_ph*q1, z );   // fct 0 (terminal)
  FFVar Eff1 = OpI( U_VEL*c1, t );        // fct 1 (accumulated)

  OCFESLV oc( &DAG );
  // graded t-mesh: coarse away from t_step, fine (0.1-wide) straddling the step at t_step=2.0
  std::vector<double> t_bnd = { 0.0, 1.0, 1.8, 1.9, 2.0, 2.1, 2.2, 3.0, 4.0, T_end };
  oc.add_domain( t, FFDom( t_bnd, coltype, n_nd ) );   // explicit-boundary (non-uniform) constructor
  oc.add_domain( z, FFDom( 0., 1.0,   n_el, coltype, n_nd ) );
  oc.add_state ( c1, {t,z} );  oc.add_state ( c2, {t,z} );
  oc.add_state ( q1, {t,z} );  oc.add_state ( q2, {t,z} );
  oc.add_state ( T,  {t,z} );
  oc.add_input ( b01, b01_nom, /*is_decision=*/true );

  auto c1g = [&]( OCFESLV::t_Coord const& cr ){ return c1feed_d( cr.at(t) )*( 1.0 - 0.5*cr.at(z) ); };
  auto c2g = [&]( OCFESLV::t_Coord const& cr ){ return c2feed_d( cr.at(t) )*( 1.0 - 0.5*cr.at(z) ); };
  oc.update_ref( c1, c1g );
  oc.update_ref( c2, c2g );
  oc.update_ref( q1, [&]( OCFESLV::t_Coord const& cr ){ return q1star_g( c1g(cr), c2g(cr), b01_nom ); } );
  oc.update_ref( q2, [&]( OCFESLV::t_Coord const& cr ){ return q2star_g( c1g(cr), c2g(cr), b01_nom ); } );
  oc.update_ref( T,  [&]( OCFESLV::t_Coord const& cr ){
    return T0_ref + F_ph*( dH1*q1star_g(c1g(cr),c2g(cr),b01_nom) + dH2*q2star_g(c1g(cr),c2g(cr),b01_nom) )/Cp_e; } );

  oc.add_input ( c1_ic, {z} );  oc.add_input ( c2_ic, {z} );
  oc.add_input ( q1_ic, {z} );  oc.add_input ( q2_ic, {z} );
  oc.add_input ( T_ic, {z}, [&]( OCFESLV::t_Coord const& cr ){ return T0_ref; }, /*is_decision=*/true );
  oc.update_ref( c1_ic, [&]( OCFESLV::t_Coord const& cr ){ return 0.0;    } );
  oc.update_ref( c2_ic, [&]( OCFESLV::t_Coord const& cr ){ return 0.0;    } );
  oc.update_ref( q1_ic, [&]( OCFESLV::t_Coord const& cr ){ return 0.0;    } );
  oc.update_ref( q2_ic, [&]( OCFESLV::t_Coord const& cr ){ return 0.0;    } );

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
  oc.add_output( Inv1, {t}, {T_end} );   // fct 0
  oc.add_output( Eff1, {z}, {1.0}   );   // fct 1

  oc.options.REDUCE.ORDER     = OCFESLV::Options::RED_FULL;
  oc.options.CLASSIFY.MODE         = OCFESLV::Options::CLASS_AUTO;
  oc.options.INTERFACE.TYPE   = OCFESLV::Options::IC_AUTO;
  oc.options.INTERFACE.IMPOSITION  = OCFESLV::Options::IC_STRONG;
  oc.options.INTERFACE.SAT_SIGMA0       = 10.0;
  oc.options.DISPLAY_LEVEL    = 0;
#if defined(CRONOS__WITH_SPQR)
  oc.options.SOLVE.FACTORIZATION = OCFESLV::Options::SOLVE_SPQR;
#else
  oc.options.SOLVE.FACTORIZATION = OCFESLV::Options::SOLVE_SUPERLU;
#endif
  oc.options.SOLVE.WARMSTART = OCFESLV::Options::BROADCAST_IC;   // marching

  if( !oc.setup() ){ std::cerr << "  setup() FAILED\n"; return 1; }
  std::vector<double> xv, inp;
  if( !oc.init( xv, inp, nullptr ) ){ std::cerr << "  init() FAILED\n"; return 1; }

  // controls b01, T_ic were flagged is_decision at add_input() and auto-registered by setup()
  size_t const ncd = oc.n_control_dof();
  size_t const ncf = oc.n_colloc_fct();
  auto const& C = oc.controls();
  size_t const ib01 = C.at( b01 ).offset;                // b01 column (keyed -> canonical FFVar order)
  size_t const itic0 = C.at( T_ic ).offset, iticN = C.at( T_ic ).ndof;   // T_ic block

  std::vector<double> p0; oc.encode_controls( inp.data(), p0 );
  std::cout << "  reduced control space: ncd=" << ncd << " (b01 @ " << ib01
            << ", T_ic block @ " << itic0 << " x " << iticN << "),  outputs ncf=" << ncf << "\n\n";

  // ---------------- reference: direct solve at the nominal controls ----------------
  std::vector<double> xvR( xv ), inpR( inp );
  OCFESLV::SolveReport const rep0 = oc.solve( xvR.data(), inpR.data(), nullptr );
  if( !rep0.converged ){ std::cerr << "  reference solve FAILED\n"; return 1; }
  std::vector<double> Fref = oc.val_functions();
  if( Fref.size() < ncf ){ std::cerr << "  reference functions missing\n"; return 1; }

  std::cout << "  graded t-mesh windows (n_march_steps) = " << oc.n_march_steps()
            << " ;  inlet step at t*=" << t_step << " (eps=" << eps_step << ")\n"
            << "  reduced control space ncd=" << ncd << ",  outputs ncf=" << ncf << "\n\n";

  // ---- S0: the graded marching solve converges THROUGH the step ----
  check_true( "S0 graded step-inlet march converged", rep0.converged );

  // ---- S1: forward == adjoint reduced Jacobian across the discontinuity (tight) ----
  std::vector<std::vector<double>> Jf( ncf, std::vector<double>( ncd, 0. ) ),
                                   Ja( ncf, std::vector<double>( ncd, 0. ) );
  {
    std::vector<double> xvf( xv ), inpf( inp );
    if( !oc.solve_fsens( xvf.data(), inpf.data(), nullptr ) ){ std::cerr << "  solve_fsens FAILED\n"; return 1; }
    std::vector<double> const& J = oc.sens_jacobian();
    for( size_t i=0;i<ncf;++i ) for( size_t j=0;j<ncd;++j ) Jf[i][j]=J[i*ncd+j];
  }
  {
    std::vector<double> xva( xv ), inpa( inp );
    if( !oc.solve_asens( xva.data(), inpa.data(), nullptr ) ){ std::cerr << "  solve_asens FAILED\n"; return 1; }
    std::vector<double> const& J = oc.sens_jacobian();
    for( size_t i=0;i<ncf;++i ) for( size_t j=0;j<ncd;++j ) Ja[i][j]=J[i*ncd+j];
  }
  { double e=0.; for( size_t i=0;i<ncf;++i ) for( size_t j=0;j<ncd;++j ) e=std::max(e,std::fabs(Jf[i][j]-Ja[i][j]));
    check_close( "S1 forward == adjoint Jacobian (through step)", e, 0., 1e-9 ); }

  // ---- S2: b01 column vs central-difference FD (independent, through the step) ----
  {
    double const h = 1e-4;
    std::vector<double> pp( p0 ), pm( p0 ); pp[ib01]+=h; pm[ib01]-=h;
    std::vector<double> Fp, Fm, inpP( inp ), inpM( inp ), xvP( xv ), xvM( xv );
    oc.decode_controls( pp, inpP.data() ); oc.solve( xvP.data(), inpP.data(), nullptr ); Fp = oc.val_functions();
    oc.decode_controls( pm, inpM.data() ); oc.solve( xvM.data(), inpM.data(), nullptr ); Fm = oc.val_functions();
    for( size_t i=0;i<ncf;++i ){
      double fd = ( Fp[i]-Fm[i] )/( 2.0*h );
      double tol = 1e-3*std::max( 1e-6, std::fabs(Jf[i][ib01]) ) + 1e-6;
      std::ostringstream nm; nm << "S2 dF" << i << "/db01 ~ FD (through step)";
      check_close( nm.str().c_str(), Jf[i][ib01], fd, tol );
    }
  }

  // ---- S3: warm-start robustness across the step -- REUSE vs BROADCAST_IC agree ----
  {
    std::vector<double> xvR( xv ), inpR( inp ), xvB( xv ), inpB( inp );
    oc.options.SOLVE.WARMSTART = OCFESLV::Options::REUSE;
    OCFESLV::SolveReport rR = oc.solve( xvR.data(), inpR.data(), nullptr );
    std::vector<double> Fr = oc.val_functions();
    oc.options.SOLVE.WARMSTART = OCFESLV::Options::BROADCAST_IC;
    OCFESLV::SolveReport rB = oc.solve( xvB.data(), inpB.data(), nullptr );
    std::vector<double> Fb = oc.val_functions();
    oc.options.SOLVE.WARMSTART = OCFESLV::Options::REUSE;
    check_true( "S3 both warm-starts converge through step", rR.converged && rB.converged );
    double e=0.; for( size_t i=0;i<ncf && i<Fr.size() && i<Fb.size();++i ) e=std::max(e,std::fabs(Fr[i]-Fb[i]));
    check_close( "S3 REUSE == BROADCAST_IC functions", e, 0., 1e-7 );
  }

  std::cout << "\n============================================================\n"
            << "  RESULT: " << g_pass << " passed, " << g_fail << " failed  -- "
            << ( g_fail == 0 ? "ALL PASS" : "SOME FAILED" ) << "\n"
            << "============================================================\n";
  return g_fail == 0 ? 0 : 1;
}
