// OCFE_mono.cpp  ---  MONOLITHIC reduced-space sensitivity validation on the PSA model (rung-9)
// ===========================================================================
// Validates the reduced-space FFOCFESLV / OCFESLV sensitivity API on a MONOLITHIC
// solve (SOLVE_MARCHING=false): the full space-time PSA system is solved at once
// (no marching), exercising OCFESLV::mono_fsens_setup / mono_fsens_apply /
// mono_asens_apply through the unified solve_fsens/solve_fsens_apply/solve_asens_apply dispatch.
// Controls: b01 (lumped) + T_ic (distributed).  Outputs: Inv1 (fct 0), Eff1 (fct 1).
//
//   M0  FFOp value (monolithic solve + eval) == direct reference
//   M1  mono adjoint Jacobian == mono forward Jacobian (full nf x ncd)     [tight]
//   M2  mono forward b01 column == central-difference FD                   [loose]
//   M3  FFOp FADType == mono forward, FFOp BADType == mono adjoint         [tight]
// M1/M3 confirm the monolithic forward and adjoint agree; M2 is the independent
// ground truth; M3 confirms the FFOp dispatches to the monolithic backend.
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

#ifndef PSA_REDUCE
#define PSA_REDUCE RED_FULL
#endif

using namespace mc;

// ---- FFOCFESLV two-map migration (2026-09-27) ----------------------------------------------------------
// The reduced operation takes EVERY declared input and constant: map 1 = the controls, map 2 = the rest.
// Map 2 here: every declared input that is NOT a control, held at its current value in @p inp as LITERAL
// constants -- what the old inp0 argument carried.  @p vals receives those values in op-input order.
static std::vector<FFModel::InputArg> psa_rest_map
( OCFESLV& oc, double const* inp, std::vector<double>& vals )
{
  std::vector<FFModel::InputArg> m;  vals.clear();
  for( auto const& [w,dom] : oc.var_declared_input() ){
    if( oc.controls().count( w ) ) continue;
    std::vector<double> const v = oc.get_input_values( w, inp );
    std::vector<FFVar> lit;  for( double x : v ){ lit.push_back( FFVar( x ) ); vals.push_back( x ); }
    m.emplace_back( w, lit );
  }
  return m;
}
// Op inputs for a direct eval(): a control-space vector followed by the map-2 values.
template <typename T> static std::vector<T> psa_opin
( std::vector<T> const& pc, std::vector<double> const& vals )
{ std::vector<T> v( pc ); for( double x : vals ) v.push_back( T( x ) ); return v; }


static double const U_VEL = 1.0, D_ax = 0.1, F_ph = 0.5, k_ldf = 2.0;
static double const qs1 = 1.0, qs2 = 1.0, b02_L = 1.0;
static double const beta1 = 2.0, beta2 = 1.0, T0_ref = 1.0;
static double const Cp_e = 1.0, G_cv = 1.0, lam = 0.1, dH1 = 1.0, dH2 = 0.5, hw = 0.5, Tw = 1.0;
static double const c0_1 = 0.5, c0_2 = 0.5, tau_in = 0.15, T_end = 5.0;
static double const b01_nom = 3.0;

static inline double c1feed_d( double t ){ return c0_1*( 1.0 - std::exp( -t/tau_in ) ); }
static inline double c2feed_d( double t ){ return c0_2*( 1.0 - std::exp( -t/tau_in ) ); }
static inline double q1star_g( double c1g, double c2g, double b01 )
{ return qs1*b01*c1g/( 1.0 + b01*c1g + b02_L*c2g ); }
static inline double q2star_g( double c1g, double c2g, double b01 )
{ return qs2*b02_L*c2g/( 1.0 + b01*c1g + b02_L*c2g ); }

static bool g_quiet = false;   // 2026-09-19: the matrix re-runs the body; only the reference cell counts
static int g_pass = 0, g_fail = 0;
static void check_close( char const* name, double got, double want, double tol )
{
  if( g_quiet ) return;
  bool ok = std::fabs( got - want ) <= tol;
  ( ok ? g_pass : g_fail )++;
  std::cout << "  " << std::left << std::setw(46) << name
            << std::scientific << std::setprecision(3)
            << " |got-want|=" << std::fabs( got - want ) << " tol=" << tol
            << "  " << ( ok ? "PASS" : "FAIL" ) << "\n";
}

static void check_true( char const* name, bool ok )
{
  if( g_quiet ) return;
  ( ok ? g_pass : g_fail )++;
  std::cout << "  " << std::left << std::setw(46) << name << "  " << ( ok ? "PASS" : "FAIL" ) << "\n";
}

// 2026-09-19: the body is now a CELL -- the model is rebuilt per (imposition, reduction) so the matrix can ask
// whether either changes the reduced-space outputs.  Defaults reproduce the original run exactly.
// EACH CELL OWNS ITS FFGraph.  The model is rebuilt per (imposition, reduction) and the FFIPDAE operator
// registers against the graph it was built in, so a cell cannot inherit another cell's registrations -- the
// same discipline OCFE_partialfold needs for FOLD_NESTED.  `DAG` below is local to this function already; the
// point is recorded here because it is the assumption the matrix rests on.
//
// 2026-09-19 RESOLVED: the reference cell passes every FFIPDAE check on Benoit's build (17 passed, 0 failed) and
// the matrix is 5/5.  A "mono solve_fsens FAILED" seen while developing this was an artefact of the slow
// sandbox, not this wrap and not the header.  Measured: both EXACT modes reproduce the reference reduced-space
// outputs to 2.6e-16 -- the imposition does not perturb F -- while IC_WEAK sits at ~3e-06, the penalty's
// O(1/sigma), consistent with OCFE_PDE27's 2.26e-05 on a coarser mesh.
static int run_cell( OCFESLV::Options::ImpositionType imp = OCFESLV::Options::IC_STRONG,
                     OCFESLV::Options::ReductionType red = OCFESLV::Options::PSA_REDUCE,
                     std::vector<double>* Fout = nullptr,
                     bool run_sensitivities = true )
{
  std::cout << "================================================================\n"
            << "  MONOLITHIC reduced-space sensitivity (mono_fsens/asens) on the PSA model (rung-9)\n"
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

  FFVar c1feed = c0_1*( 1.0 - exp( -t/tau_in ) );
  FFVar c2feed = c0_2*( 1.0 - exp( -t/tau_in ) );
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
  oc.add_domain( t, FFDom( 0., T_end, n_el, coltype, n_nd ) );
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

  oc.options.REDUCE.ORDER     = red;
  oc.options.CLASSIFY.MODE         = OCFESLV::Options::CLASS_AUTO;
  oc.options.INTERFACE.TYPE   = OCFESLV::Options::IC_AUTO;
  oc.options.INTERFACE.IMPOSITION  = imp;
  oc.options.INTERFACE.SAT_SIGMA0       = 10.0;
  oc.options.DISPLAY_LEVEL    = 0;
#if defined(CRONOS__WITH_SPQR)
  oc.options.SOLVE.FACTORIZATION = OCFESLV::Options::SOLVE_SPQR;
#else
  oc.options.SOLVE.FACTORIZATION = OCFESLV::Options::SOLVE_SUPERLU;
#endif
  oc.options.SOLVE.WARMSTART = OCFESLV::Options::BROADCAST_IC;   // marching
  oc.options.SOLVE.MARCHING     = false;   // MONOLITHIC: solve the full space-time system (no marching)

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

  check_true( "model is monolithic (is_marching()==false)", !oc.is_marching() );

  // ---------------- direct monolithic reference: solve + eval outputs ----------------
  std::vector<double> xvR( xv ), inpR( inp );
  OCFESLV::SolveReport const rep0 = oc.solve( xvR.data(), inpR.data(), nullptr );
  if( !rep0.converged ){ std::cerr << "  monolithic reference solve FAILED\n"; return 1; }
  std::vector<double> eqnR( oc.n_colloc_eqn(), 0. ), Fref( ncf, 0. );
  if( !oc.eval( eqnR.data(), Fref.data(), xvR.data(), inpR.data(), nullptr ) ){ std::cerr << "  monolithic eval FAILED\n"; return 1; }

  // ---------------- reduced-space operation on the monolithic OCFESLV ----------------
  FFGraph rdag;
  std::vector<FFVar> pctrl( ncd );
  for( size_t i=0; i<ncd; ++i ){ std::ostringstream os; os << "p" << i; pctrl[i] = rdag.add_var( os.str() ); }
  std::vector<double> vrest;                        // map-2 values, op-input order
  FFOCFESLV ffred;
  std::vector<FFVar> F = ffred( FFOCFESLV::controls_map( oc, pctrl ), psa_rest_map( oc, inp.data(), vrest ), &oc, FFOCFESLV::SHALLOW, "PSA_mono" );
  check_true( "FFOCFESLV returns ncf outputs", F.size() == ncf );

  // M0: FFOp value (monolithic solve + eval) == direct reference
  {
    std::vector<double> Fval( ncf, 0. );
    ffred.eval( (unsigned)ncf, Fval.data(), (unsigned)( ncd + vrest.size() ), psa_opin( p0, vrest ).data(), nullptr );
    double e=0.; for( size_t j=0;j<ncf;++j ) e=std::max(e,std::fabs(Fval[j]-Fref[j]));
    check_close( "M0 FFOp value == monolithic solve", e, 0., 1e-9 );
  }

  // ---------------- monolithic reduced Jacobian, forward and adjoint (unified API) ----------------
  std::cout << "  reduced Jacobian dF/dp : " << ncf << " x " << ncd
            << "  (mono forward = " << ncd << " applies, mono adjoint = " << ncf << " sweeps)\n\n";
  std::vector<std::vector<double>> Jfwd( ncf, std::vector<double>( ncd, 0. ) );
  {
    std::vector<double> xvf( xv ), inpf( inp );
    if( !run_sensitivities ){ if( Fout ) *Fout = Fref; return 0; }   // matrix cells compare F only
    if( !oc.solve_fsens( xvf.data(), inpf.data(), nullptr ) ){ std::cerr << "  mono solve_fsens FAILED\n"; return 1; }
    std::vector<double> const& Jf = oc.sens_jacobian();   // nf x ncd, row-major
    if( Jf.size() < ncf*ncd ){ std::cerr << "  mono forward Jacobian missing\n"; return 1; }
    for( size_t i=0; i<ncf; ++i ) for( size_t j=0; j<ncd; ++j ) Jfwd[i][j] = Jf[ i*ncd + j ];
  }
  std::vector<std::vector<double>> Jadj( ncf, std::vector<double>( ncd, 0. ) );
  {
    std::vector<double> xva( xv ), inpa( inp );
    if( !oc.solve_asens( xva.data(), inpa.data(), nullptr ) ){ std::cerr << "  mono solve_asens FAILED\n"; return 1; }
    std::vector<double> const& Ja = oc.sens_jacobian();   // nf x ncd, row-major
    if( Ja.size() < ncf*ncd ){ std::cerr << "  mono adjoint Jacobian missing\n"; return 1; }
    for( size_t i=0; i<ncf; ++i ) for( size_t j=0; j<ncd; ++j ) Jadj[i][j] = Ja[ i*ncd + j ];
  }

  // M1: mono adjoint == mono forward (full Jacobian) [tight]
  { double e=0.; for( size_t i=0;i<ncf;++i ) for( size_t j=0;j<ncd;++j ) e=std::max(e,std::fabs(Jadj[i][j]-Jfwd[i][j]));
    check_close( "M1 mono adjoint Jacobian == mono forward", e, 0., 1e-9 ); }

  // ---- SP: reduced dependency pattern -- must COVER the actual Jacobian nonzeros (no false
  // negatives, the dangerous bug); and for the fully-coupled PSA it must be dense. ----
  {
    std::vector<std::vector<bool>> const& P = oc.reduced_dependency_pattern();
    bool shape = ( P.size() == ncf );
    for( size_t i=0;i<P.size();++i ) shape = shape && ( P[i].size() == ncd );
    check_true( "SP0 reduced_dependency_pattern shape nf x ncd", shape );
    bool cover = true, dense = true;
    for( size_t i=0;i<ncf;++i ) for( size_t j=0;j<ncd;++j ){
      bool const nz = std::fabs( Jfwd[i][j] ) > 1e-9;
      bool const p  = ( i<P.size() && j<P[i].size() && P[i][j] );
      if( nz && !p ) cover = false;    // structural pattern missed an actual dependency
      if( !p )       dense = false;
    }
    check_true( "SP1 pattern covers actual Jacobian nonzeros", cover );
    check_true( "SP2 coupled PSA pattern is dense", dense );
  }

  // M2: mono forward b01 column == central-difference FD (independent) [loose]
  {
    double const eps = 1e-4;
    std::vector<double> pp( p0 ), pm( p0 ); pp[ib01]+=eps; pm[ib01]-=eps;
    std::vector<double> Fp( ncf, 0. ), Fm( ncf, 0. );
    std::vector<double> inpP( inp ), inpM( inp ), xvP( xv ), xvM( xv ), eqnb( oc.n_colloc_eqn(), 0. );
    oc.decode_controls( pp, inpP.data() ); oc.solve( xvP.data(), inpP.data(), nullptr ); oc.eval( eqnb.data(), Fp.data(), xvP.data(), inpP.data(), nullptr );
    oc.decode_controls( pm, inpM.data() ); oc.solve( xvM.data(), inpM.data(), nullptr ); oc.eval( eqnb.data(), Fm.data(), xvM.data(), inpM.data(), nullptr );
    for( size_t i=0; i<ncf; ++i ){
      double fd = ( Fp[i]-Fm[i] )/( 2.0*eps );
      double tol = 1e-3*std::max( 1e-6, std::fabs(Jfwd[i][ib01]) ) + 1e-6;
      std::ostringstream nm; nm << "M2 mono dF" << i << "/db01 ~ FD";
      check_close( nm.str().c_str(), Jfwd[i][ib01], fd, tol );
    }
  }

  // M3: FFOp forward (FADType) == Jfwd, and reverse (BADType) == Jadj -- confirms FFOp
  // dispatches to the monolithic backend for both AD modes.
  {
    std::vector<FADType<double>> xF( ncd ), yF( ncf );
    for( size_t i=0;i<ncd;++i ){ xF[i]=p0[i]; xF[i].diff((unsigned)i,(unsigned)ncd); }
    ffred.eval( (unsigned)ncf, yF.data(), (unsigned)( ncd + vrest.size() ), psa_opin( xF, vrest ).data(), nullptr );
    double ef=0.; for( size_t i=0;i<ncf;++i ) for( size_t j=0;j<ncd;++j ) ef=std::max(ef,std::fabs(yF[i].deriv((unsigned)j)-Jfwd[i][j]));
    check_close( "M3 FFOp FADType == mono forward Jacobian", ef, 0., 1e-9 );

    double eb=0.;
    for( size_t irow=0; irow<ncf; ++irow ){
      std::vector<BADType<double>> xB( ncd ), yB( ncf );
      for( size_t j=0;j<ncd;++j ) xB[j]=p0[j];
      ffred.eval( (unsigned)ncf, yB.data(), (unsigned)( ncd + vrest.size() ), psa_opin( xB, vrest ).data(), nullptr );
      yB[irow].diff( 0, 1 );
      for( size_t j=0;j<ncd;++j ) eb=std::max(eb,std::fabs(xB[j].d(0)-Jadj[irow][j]));
    }
    check_close( "M3 FFOp BADType == mono adjoint Jacobian", eb, 0., 1e-9 );
  }

  // ---- R: primal-solve reuse (SOLVE_REUSE) -- correctness + cache hit/miss ----
  {
    // second solve at the SAME nominal input must be served from the cache (no Newton), same functions
    std::vector<double> xa( xv ), ia( inp );
    OCFESLV::SolveReport ra = oc.solve( xa.data(), ia.data(), nullptr );   // cold or warm depending on prior
    std::vector<double> Fa = oc.val_functions();
    std::vector<double> xb( xv ), ib( inp );
    OCFESLV::SolveReport rb = oc.solve( xb.data(), ib.data(), nullptr );   // same input -> reuse
    std::vector<double> Fb = oc.val_functions();
    check_true( "R1 repeat solve at same input is reused", rb.reused );
    double e1=0.; for( size_t j=0;j<ncf && j<Fa.size() && j<Fb.size();++j ) e1=std::max(e1,std::fabs(Fa[j]-Fb[j]));
    check_close( "R1 reused functions == recomputed", e1, 0., 1e-12 );

    // solve at a DIFFERENT input must miss the cache (fresh Newton)
    std::vector<double> pp( p0 ); pp[ib01]+=1e-2;
    std::vector<double> ic( inp ), xc( xv );
    oc.decode_controls( pp, ic.data() );
    OCFESLV::SolveReport rc = oc.solve( xc.data(), ic.data(), nullptr );
    check_true( "R2 solve at a different input is not reused", !rc.reused && rc.converged );

    // with SOLVE_REUSE off, a repeat at the same input recomputes (still correct)
    bool const save = oc.options.SOLVE.REUSE;
    oc.options.SOLVE.REUSE = false;
    std::vector<double> xd( xv ), id( inp ); oc.decode_controls( p0, id.data() );
    OCFESLV::SolveReport rd = oc.solve( xd.data(), id.data(), nullptr );
    std::vector<double> Fd = oc.val_functions();
    oc.options.SOLVE.REUSE = save;
    check_true( "R3 SOLVE_REUSE=off recomputes (not reused)", !rd.reused );
    double e3=0.; for( size_t j=0;j<ncf && j<Fd.size() && j<Fref.size();++j ) e3=std::max(e3,std::fabs(Fd[j]-Fref[j]));
    check_close( "R3 recomputed functions == reference", e3, 0., 1e-9 );
  }

  if( Fout ) *Fout = Fref;   // the reduced-space outputs, for the matrix's comparison
  return 0;
}

int main()
{
  std::vector<double> Fref0;
  if( run_cell( OCFESLV::Options::IC_STRONG, OCFESLV::Options::PSA_REDUCE, &Fref0 ) != 0 ) return 1;

  // ---------------------------------------------------------------------------------------------------------
  // 2026-09-19 -- THE SYSTEMATIC MATRIX: 3 impositions x 2 reductions.  PSA4 is the corpus's reduced-space
  // (FFIPDAE) test -- forward/adjoint Jacobians and solve reuse -- and the whole PSA family ran IC_STRONG with
  // PSA_REDUCE only.  Each cell runs the model build + reference solve + output evaluation ONLY (the Jacobian
  // and reuse checks belong to the reference cell); what it compares is the reduced-space outputs F, which
  // every Jacobian and reuse check downstream depends on.  Compared against the reference cell above:
  // exact modes to 1e-6 relative, IC_WEAK to 1e-2,
  // because a penalty imposes continuity only to O(1/sigma) and must differ at that level (measured on PDE27:
  // 2.26e-05).  The cell's own check_* calls are silenced (g_quiet) so only the reference cell scores.
  // A cell expected to differ gets a documented XFAIL naming the mechanism, never a loosened bar.
  // ---------------------------------------------------------------------------------------------------------
  std::cout << "\n---- SYSTEMATIC MATRIX: imposition x reduction (reduced-space outputs F) ----\n";
  {
    struct MCell { char const* imp; OCFESLV::Options::ImpositionType it;
                   char const* red; OCFESLV::Options::ReductionType rt; };
    MCell const cells[5] = {
      { "WEAK",   OCFESLV::Options::IC_WEAK,   "PSA_REDUCE", OCFESLV::Options::PSA_REDUCE },
      { "TRACE",  OCFESLV::Options::IC_TRACE,  "PSA_REDUCE", OCFESLV::Options::PSA_REDUCE },
      { "WEAK",   OCFESLV::Options::IC_WEAK,   "RED_MAIN",   OCFESLV::Options::RED_MAIN   },
      { "TRACE",  OCFESLV::Options::IC_TRACE,  "RED_MAIN",   OCFESLV::Options::RED_MAIN   },
      { "STRONG", OCFESLV::Options::IC_STRONG, "RED_MAIN",   OCFESLV::Options::RED_MAIN   } };
    bool matrix_ok = true;
    g_quiet = true;
    std::cout << "  " << std::left << std::setw(10) << "IMPOSITION" << std::setw(13) << "REDUCTION"
              << std::setw(8) << "setup" << std::setw(13) << "max rel.dev of F" << "  VERDICT\n";
    for( auto const& mc : cells ){
      std::vector<double> Fc;
      int const rc = run_cell( mc.it, mc.rt, &Fc, false );   // F only: the Jacobian/reuse checks are the
                                                            // reference cell's business, not the matrix's
      double dev = 0.;
      if( rc == 0 && Fc.size() == Fref0.size() )
        for( size_t i = 0; i < Fc.size(); ++i ){
          double const m = std::max( std::fabs( Fc[i] ), std::fabs( Fref0[i] ) );
          if( m > 0. ) dev = std::max( dev, std::fabs( Fc[i] - Fref0[i] ) / m ); }
      double const tol = ( mc.it == OCFESLV::Options::IC_WEAK ) ? 1e-2 : 1e-6;
      bool const xfail = false;
      bool const cell  = ( rc == 0 ) && Fc.size() == Fref0.size() && dev <= tol;
      if( !cell && !xfail ) matrix_ok = false;
      std::cout << "  " << std::left << std::setw(10) << mc.imp << std::setw(13) << mc.red
                << std::setw(8) << ( rc == 0 ? "ok" : "FAILED" ) << std::scientific << std::setprecision(3)
                << std::setw(15) << dev
                << ( cell ? ( xfail ? "PASS (unexpected: the documented defect is gone?)" : "PASS" )
                          : ( xfail ? "XFAIL (documented)" : "FAIL" ) ) << "\n";
    }
    g_quiet = false;
    check_true( "matrix: every cell reproduces the reference reduced-space outputs", matrix_ok );
    std::cout << "  MATRIX (vs STRONG/PSA_REDUCE; exact 1e-6, WEAK 1e-2): "
              << ( matrix_ok ? "PASS" : "FAIL" ) << "\n";
  }

  std::cout << "\n============================================================\n"
            << "  RESULT: " << g_pass << " passed, " << g_fail << " failed  -- "
            << ( g_fail == 0 ? "ALL PASS" : "SOME FAILED" ) << "\n"
            << "============================================================\n";
  return g_fail == 0 ? 0 : 1;
}
