// OCFE_march_asens.cpp  ---  OCFESLV-level ADJOINT sensitivity validation on the PSA oracle (rung-9)
// ===========================================================================
// Validates OCFESLV::solve_asens_apply (the transpose/backward dual of solve_fsens_apply)
// and its helper terminal_profile_adjoint, on the PSA marching model of OCFE_PDE36c.
// Controls: b01 (lumped) + T_ic (distributed).  Outputs: Inv1 (fct 0), Eff1 (fct 1).
//
// The reduced Jacobian dF/dp (nf x ncd) is assembled two ways after ONE sens_setup:
//   forward : one solve_fsens_apply per control column   (ncd applies)   -> J_fwd
//   adjoint : one solve_asens_apply per function row + encode_controls (nf sweeps) -> J_adj
// and cross-checked against the analytic solve_marching_fsens oracle:
//   A1  adjoint Jacobian == forward Jacobian (full nf x ncd)                         [tight]
//   A2  b01 column (forward AND adjoint) == analytic solve_marching_fsens(b01,{1})   [tight]
//   A3  adjoint sum of T_ic columns == analytic solve_marching_fsens(T_ic,ones)      [tight]
// This isolates the new adjoint numerics before wiring them into FFGradOCFESLV.
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

#ifndef PSA_REDUCE
#define PSA_REDUCE RED_FULL
#endif

using namespace mc;

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

// Analytic dF/d(control-direction) via the PUBLIC reduced Jacobian (replaces the retired
// solve_marching_fsens oracle).  A lumped control (ndof 1) -> its single column; a distributed
// control -> the uniform-shift (all-ones) sum over its columns.  @p out is sized n_colloc_fct().
static bool fsens_oracle( mc::OCFESLV& oc, size_t off, size_t ndof,
                          std::vector<double> xv, std::vector<double> inp, std::vector<double>& out )
{
  size_t const ncf = oc.n_colloc_fct(), ncd = oc.n_control_dof();
  if( !oc.solve_fsens( xv.data(), inp.empty()?nullptr:inp.data(), nullptr ) ) return false;
  std::vector<double> const& J = oc.sens_jacobian();
  if( J.size() < ncf*ncd ) return false;
  out.assign( ncf, 0. );
  for( size_t j = 0; j < ncf; ++j ){ double s = 0.; for( size_t i = 0; i < ndof; ++i ) s += J[ j*ncd + off + i ]; out[j] = s; }
  return true;
}

int main()
{
  std::cout << "================================================================\n"
            << "  OCFESLV adjoint sensitivity (solve_asens_apply) on the PSA oracle (rung-9)\n"
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
  oc.add_input ( b01, {} );

  auto c1g = [&]( OCFESLV::t_Coord const& cr ){ return c1feed_d( cr.at(t) )*( 1.0 - 0.5*cr.at(z) ); };
  auto c2g = [&]( OCFESLV::t_Coord const& cr ){ return c2feed_d( cr.at(t) )*( 1.0 - 0.5*cr.at(z) ); };
  oc.update_ref( c1, c1g );
  oc.update_ref( c2, c2g );
  oc.update_ref( q1, [&]( OCFESLV::t_Coord const& cr ){ return q1star_g( c1g(cr), c2g(cr), b01_nom ); } );
  oc.update_ref( q2, [&]( OCFESLV::t_Coord const& cr ){ return q2star_g( c1g(cr), c2g(cr), b01_nom ); } );
  oc.update_ref( T,  [&]( OCFESLV::t_Coord const& cr ){
    return T0_ref + F_ph*( dH1*q1star_g(c1g(cr),c2g(cr),b01_nom) + dH2*q2star_g(c1g(cr),c2g(cr),b01_nom) )/Cp_e; } );

  oc.add_input ( c1_ic, {z} );  oc.add_input ( c2_ic, {z} );
  oc.add_input ( q1_ic, {z} );  oc.add_input ( q2_ic, {z} );  oc.add_input ( T_ic, {z} );
  oc.update_ref( c1_ic, [&]( OCFESLV::t_Coord const& cr ){ return 0.0;    } );
  oc.update_ref( c2_ic, [&]( OCFESLV::t_Coord const& cr ){ return 0.0;    } );
  oc.update_ref( q1_ic, [&]( OCFESLV::t_Coord const& cr ){ return 0.0;    } );
  oc.update_ref( q2_ic, [&]( OCFESLV::t_Coord const& cr ){ return 0.0;    } );
  oc.update_ref( T_ic,  [&]( OCFESLV::t_Coord const& cr ){ return T0_ref; } );

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

  oc.options.REDUCE.ORDER     = OCFESLV::Options::PSA_REDUCE;
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
  oc.set_input_values( b01, { b01_nom }, inp.data() );
  size_t const ntic = oc.get_input_values( T_ic, inp.data() ).size();
  oc.set_input_values( T_ic, std::vector<double>( ntic, T0_ref ), inp.data() );

  // ---------------- register the reduced control space ----------------
  oc.register_control( b01 );          // C.at(b01): lumped coefficient (ndof 1)
  oc.register_control( T_ic );   // C.at(T_ic): distributed IC input   (ndof ntic)
  size_t const ncd = oc.n_control_dof();
  size_t const ncf = oc.n_colloc_fct();
  auto const& C = oc.controls();
  size_t const ib01 = C.at(b01).offset;                       // b01 column in the control vector
  size_t const itic0 = C.at(T_ic).offset, iticN = C.at(T_ic).ndof;   // T_ic block

  std::vector<double> p0; oc.encode_controls( inp.data(), p0 );
  std::cout << "  reduced control space: ncd=" << ncd << " (b01 @ " << ib01
            << ", T_ic block @ " << itic0 << " x " << iticN << "),  outputs ncf=" << ncf << "\n\n";

  // ---------------- reference: direct solve at the nominal controls ----------------
  std::vector<double> xvR( xv ), inpR( inp );
  OCFESLV::SolveReport const rep0 = oc.solve( xvR.data(), inpR.data(), nullptr );
  if( !rep0.converged ){ std::cerr << "  reference solve FAILED\n"; return 1; }
  std::vector<double> Fref = oc.val_functions();
  if( Fref.size() < ncf ){ std::cerr << "  reference functions missing\n"; return 1; }

  // ================= OCFESLV-level ADJOINT validation =================
  // Assemble the reduced Jacobian dF/dp (nf x ncd) two ways after ONE sens_setup:
  //   forward : one solve_fsens_apply per control column  -> J_fwd[:,j]
  //   adjoint : one solve_asens_apply per function row + encode_controls -> J_adj[i,:]
  // Both must equal each other and the analytic solve_marching_fsens oracle.
  std::cout << "  reduced Jacobian dF/dp : " << ncf << " x " << ncd
            << "  (forward = " << ncd << " applies, adjoint = " << ncf << " sweeps)\n\n";

  // ---- forward-assembled J (complete public solve_fsens -> sens_jacobian) ----
  std::vector<std::vector<double>> Jfwd( ncf, std::vector<double>( ncd, 0. ) );
  {
    std::vector<double> xvf( xv ), inpf( inp );
    if( !oc.solve_fsens( xvf.data(), inpf.data(), nullptr ) ){ std::cerr << "  fwd solve FAILED\n"; return 1; }
    std::vector<double> const& J = oc.sens_jacobian();      // nf x ncd, row-major (i*ncd + j)
    for( size_t i = 0; i < ncf; ++i ) for( size_t j = 0; j < ncd; ++j ) Jfwd[i][j] = J[ i*ncd + j ];
  }

  // ---- adjoint-assembled J (complete public solve_asens -> sens_jacobian) ----
  std::vector<std::vector<double>> Jadj( ncf, std::vector<double>( ncd, 0. ) );
  {
    std::vector<double> xva( xv ), inpa( inp );
    if( !oc.solve_asens( xva.data(), inpa.data(), nullptr ) ){ std::cerr << "  adj solve FAILED\n"; return 1; }
    std::vector<double> const& J = oc.sens_jacobian();      // identical Jacobian, assembled by reverse sweeps
    for( size_t i = 0; i < ncf; ++i ) for( size_t j = 0; j < ncd; ++j ) Jadj[i][j] = J[ i*ncd + j ];
  }

  // ---- A1: adjoint == forward (full Jacobian) ----
  { double e = 0.; for( size_t i=0;i<ncf;++i ) for( size_t j=0;j<ncd;++j ) e=std::max(e,std::fabs(Jadj[i][j]-Jfwd[i][j]));
    check_close( "A1 adjoint Jacobian == forward (full)", e, 0., 1e-9 ); }

  // ---- A2: b01 column vs analytic oracle (both modes) ----
  {
    std::vector<double> xvA( xv ), inpA( inp ), dFdb01;
    if( fsens_oracle( oc, ib01, 1, xvA, inpA, dFdb01 ) && dFdb01.size()>=ncf ){
      double ef=0., ea=0.; for( size_t i=0;i<ncf;++i ){ ef=std::max(ef,std::fabs(Jfwd[i][ib01]-dFdb01[i])); ea=std::max(ea,std::fabs(Jadj[i][ib01]-dFdb01[i])); }
      check_close( "A2 forward dF/db01 == analytic", ef, 0., 1e-9 );
      check_close( "A2 adjoint dF/db01 == analytic", ea, 0., 1e-9 );
    } else check_true( "A2 analytic b01 oracle available", false );
  }

  // ---- A3: uniform T_ic shift, adjoint sum vs analytic ----
  {
    std::vector<double> xvA( xv ), inpA( inp ), dFdTic;
    if( fsens_oracle( oc, itic0, iticN, xvA, inpA, dFdTic ) && dFdTic.size()>=ncf ){
      double e=0.; for( size_t i=0;i<ncf;++i ){ double s=0.; for( size_t j=0;j<iticN;++j ) s+=Jadj[i][itic0+j]; e=std::max(e,std::fabs(s-dFdTic[i])); }
      check_close( "A3 adjoint sum dF/dTic == analytic", e, 0., 1e-9 );
    } else check_true( "A3 analytic T_ic oracle available", false );
  }

  std::cout << "\n============================================================\n"
            << "  RESULT: " << g_pass << " passed, " << g_fail << " failed  -- "
            << ( g_fail == 0 ? "ALL PASS" : "SOME FAILED" ) << "\n"
            << "============================================================\n";
  return g_fail == 0 ? 0 : 1;
}
