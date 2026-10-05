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
//   T3  forward Jacobian via FFOCFESLV::eval<FADType>     seeded identity over the ncd controls
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
            << "  FFOCFESLV reduced-space operation on the PSA oracle (rung-9)\n"
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
  size_t const ntic = oc.get_input_values( T_ic, inp.data() ).size();

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

  // ---------------- build the reduced-space operation ----------------
  FFGraph rdag;
  std::vector<FFVar> pctrl( ncd );
  for( size_t i=0; i<ncd; ++i ){ std::ostringstream os; os << "p" << i; pctrl[i] = rdag.add_var( os.str() ); }

  std::vector<double> vrest;                        // map-2 values, op-input order
  FFOCFESLV ffred;
  std::vector<FFVar> F = ffred( FFOCFESLV::controls_map( oc, pctrl ), psa_rest_map( oc, inp.data(), vrest ), &oc, FFOCFESLV::SHALLOW, "PSA_rung9" );
  check_true( "FFOCFESLV returns ncf outputs", F.size() == ncf );

  // ================= T1: value via direct FFOCFESLV::eval<double> =================
  {
    std::vector<double> Fval( ncf, 0. );
    ffred.eval( (unsigned)ncf, Fval.data(), (unsigned)( ncd + vrest.size() ), psa_opin( p0, vrest ).data(), nullptr );
    double e = 0.; for( size_t j=0; j<ncf; ++j ) e = std::max( e, std::fabs( Fval[j] - Fref[j] ) );
    check_close( "T1 value (direct eval) == solve", e, 0., 1e-9 );
  }

  // ================= T2: value via DAG round-trip (exercises copy ctor) =================
  {
    std::vector<double> FvalDag( ncf, 0. );
    rdag.eval( F, FvalDag, pctrl, p0 );
    double e = 0.; for( size_t j=0; j<ncf; ++j ) e = std::max( e, std::fabs( FvalDag[j] - Fref[j] ) );
    check_close( "T2 value (DAG round-trip) == solve", e, 0., 1e-9 );
  }

  // ================= T3: forward Jacobian via FFOCFESLV::eval<FADType<double>> =================
  std::vector<FADType<double>> xF( ncd ), yF( ncf );
  for( size_t i=0; i<ncd; ++i ){ xF[i] = p0[i]; xF[i].diff( (unsigned)i, (unsigned)ncd ); }
  ffred.eval( (unsigned)ncf, yF.data(), (unsigned)( ncd + vrest.size() ), psa_opin( xF, vrest ).data(), nullptr );

  // 3v: primal values carried by the forward sweep
  {
    double e = 0.; for( size_t j=0; j<ncf; ++j ) e = std::max( e, std::fabs( yF[j].val() - Fref[j] ) );
    check_close( "T3v F-sweep values == solve", e, 0., 1e-9 );
  }

  // 3a: b01 column vs analytic forward-sensitivity march
  {
    std::vector<double> xvA( xv ), inpA( inp ), dFdb01;
    bool ok = fsens_oracle( oc, ib01, 1, xvA, inpA, dFdb01 );
    if( ok && dFdb01.size() >= ncf ){
      double e = 0.; for( size_t j=0; j<ncf; ++j ) e = std::max( e, std::fabs( yF[j].deriv((unsigned)ib01) - dFdb01[j] ) );
      check_close( "T3a dF/db01 (FFOp) == analytic march", e, 0., 1e-9 );
    }
    else check_true( "T3a analytic b01 march available", false );
  }

  // 3b: uniform T_ic shift (sum of T_ic columns) vs analytic march
  {
    std::vector<double> xvA( xv ), inpA( inp ), dFdTic;
    bool ok = fsens_oracle( oc, itic0, iticN, xvA, inpA, dFdTic );
    if( ok && dFdTic.size() >= ncf ){
      double e = 0.;
      for( size_t j=0; j<ncf; ++j ){
        double s = 0.; for( size_t i=0; i<iticN; ++i ) s += yF[j].deriv( (unsigned)( itic0 + i ) );
        e = std::max( e, std::fabs( s - dFdTic[j] ) );
      }
      check_close( "T3b sum dF/dTic cols (FFOp) == analytic", e, 0., 1e-9 );
    }
    else check_true( "T3b analytic T_ic march available", false );
  }

  // 3c: b01 column vs central-difference FD through the operation
  {
    double const eps = 1e-4;
    std::vector<double> pp( p0 ), pm( p0 );
    pp[ib01] += eps; pm[ib01] -= eps;
    std::vector<double> Fp( ncf, 0. ), Fm( ncf, 0. );
    ffred.eval( (unsigned)ncf, Fp.data(), (unsigned)( ncd + vrest.size() ), psa_opin( pp, vrest ).data(), nullptr );
    ffred.eval( (unsigned)ncf, Fm.data(), (unsigned)( ncd + vrest.size() ), psa_opin( pm, vrest ).data(), nullptr );
    for( size_t j=0; j<ncf; ++j ){
      double fd = ( Fp[j] - Fm[j] )/( 2.0*eps );
      double an = yF[j].deriv( (unsigned)ib01 );
      double tol = 1e-3*std::max( 1e-6, std::fabs( an ) ) + 1e-6;
      std::ostringstream nm; nm << "T3c dF" << j << "/db01 (FFOp) ~ FD";
      check_close( nm.str().c_str(), an, fd, tol );
    }
  }

  // ================= T4: symbolic Jacobian via SFAD (deriv -> FFGradOCFESLV) =================
  // SFAD triggers FFOCFESLV::eval<FADType<FFVar>>, which creates FFGradOCFESLV
  // derivative nodes; evaluating that DAG in double arithmetic triggers
  // FFGradOCFESLV::eval<double>.  The result must equal the numeric forward
  // Jacobian yF from T3 (both are march_fsens on one-hot control directions).
  {
    auto sJac = rdag.SFAD( F, pctrl );
    std::vector<unsigned> const& si = std::get<0>( sJac );
    std::vector<unsigned> const& sj = std::get<1>( sJac );
    std::vector<FFVar>    const& sd = std::get<2>( sJac );

    check_true( "T4 SFAD nnz == ncf*ncd (dense Jacobian)", sd.size() == ncf*ncd );

    std::vector<double> sdv( sd.size(), 0. );
    rdag.eval( sd, sdv, pctrl, p0 );

    double e = 0.; size_t cnt = 0; bool idx_ok = true;
    for( size_t k=0; k<sd.size(); ++k ){
      if( si[k] >= ncf || sj[k] >= ncd ){ idx_ok = false; continue; }
      e = std::max( e, std::fabs( sdv[k] - yF[si[k]].deriv( (unsigned)sj[k] ) ) );
      ++cnt;
    }
    check_true ( "T4 SFAD indices in range", idx_ok && cnt == ncf*ncd );
    check_close( "T4 SFAD Jacobian == FADType<double>", e, 0., 1e-9 );
  }

  // ================= T5: COPY policy (deep-copied OCFESLV) alongside SHALLOW =================
  // Exercises OCFESLV deep-copy for BOTH the primal op and the FFGradOCFESLV nodes
  // created inside SFAD (gradPolicy inherits COPY when the primal owns its OCFESLV).
  // COPY correctness hinges on deep_copy_from carrying the control registry; if it
  // does not, decode_controls on the copy fails -- caught and reported as a clean FAIL.
  try {
    FFGraph cdag;
    std::vector<FFVar> pc( ncd );
    for( size_t i=0; i<ncd; ++i ){ std::ostringstream os; os << "pc" << i; pc[i] = cdag.add_var( os.str() ); }

    FFOCFESLV ffredC;
    std::vector<FFVar> FC = ffredC( FFOCFESLV::controls_map( oc, pc ), psa_rest_map( oc, inp.data(), vrest ), &oc, FFOCFESLV::COPY, "PSA_copy" );
    check_true( "T5 COPY constructs ncf outputs", FC.size() == ncf );

    // value via direct eval on the deep-copied OCFESLV
    std::vector<double> FvalC( ncf, 0. );
    ffredC.eval( (unsigned)ncf, FvalC.data(), (unsigned)( ncd + vrest.size() ), psa_opin( p0, vrest ).data(), nullptr );
    { double e=0.; for( size_t j=0; j<ncf; ++j ) e=std::max(e,std::fabs(FvalC[j]-Fref[j]));
      check_close( "T5 COPY value == solve", e, 0., 1e-9 ); }

    // value via DAG round-trip on the copy
    std::vector<double> FvalCdag( ncf, 0. );
    cdag.eval( FC, FvalCdag, pc, p0 );
    { double e=0.; for( size_t j=0; j<ncf; ++j ) e=std::max(e,std::fabs(FvalCdag[j]-Fref[j]));
      check_close( "T5 COPY DAG round-trip == solve", e, 0., 1e-9 ); }

    // forward Jacobian on the copy vs the SHALLOW Jacobian yF (full)
    std::vector<FADType<double>> xC( ncd ), yC( ncf );
    for( size_t i=0; i<ncd; ++i ){ xC[i]=p0[i]; xC[i].diff((unsigned)i,(unsigned)ncd); }
    ffredC.eval( (unsigned)ncf, yC.data(), (unsigned)( ncd + vrest.size() ), psa_opin( xC, vrest ).data(), nullptr );
    { double e=0.; for( size_t j=0; j<ncf; ++j ) for( size_t i=0;i<ncd;++i )
        e=std::max(e,std::fabs(yC[j].deriv((unsigned)i)-yF[j].deriv((unsigned)i)));
      check_close( "T5 COPY forward Jacobian == SHALLOW", e, 0., 1e-9 ); }

    // symbolic Jacobian on the copy (SFAD -> FFGradOCFESLV COPY nodes) vs yC
    auto sJC = cdag.SFAD( FC, pc );
    std::vector<unsigned> const& ci = std::get<0>( sJC );
    std::vector<unsigned> const& cj = std::get<1>( sJC );
    std::vector<FFVar>    const& cd = std::get<2>( sJC );
    std::vector<double> cdv( cd.size(), 0. );
    cdag.eval( cd, cdv, pc, p0 );
    { double e=0.; for( size_t k=0;k<cd.size();++k ) if( ci[k]<ncf && cj[k]<ncd )
        e=std::max(e,std::fabs(cdv[k]-yC[ci[k]].deriv((unsigned)cj[k])));
      check_close( "T5 COPY SFAD Jacobian == COPY FADType", e, 0., 1e-9 ); }
  }
  catch( std::exception const& ex ){
    std::string msg = std::string("T5 COPY policy threw (deep_copy carries registry?): ") + ex.what();
    check_true( msg.c_str(), false );
  }

  // ================= T6: reverse-mode (BADType) gradient == forward Jacobian =================
  // Exercises FFOCFESLV::eval<BADType<double>> (adjoint march + linearization).  Per output
  // row: fresh control B-vars, seed that output's adjoint, read the input adjoints = row of dF/dp.
  {
    double e = 0.;
    for( size_t irow=0; irow<ncf; ++irow ){
      std::vector<BADType<double>> xB( ncd ), yB( ncf );
      for( size_t j=0; j<ncd; ++j ) xB[j] = p0[j];
      ffred.eval( (unsigned)ncf, yB.data(), (unsigned)( ncd + vrest.size() ), psa_opin( xB, vrest ).data(), nullptr );
      yB[irow].diff( 0, 1 );
      for( size_t j=0; j<ncd; ++j )
        e = std::max( e, std::fabs( xB[j].d(0) - yF[irow].deriv((unsigned)j) ) );
    }
    check_close( "T6 BADType gradient == forward Jacobian", e, 0., 1e-9 );
  }

  // ================= T7: FFGradOCFESLV FORWARD vs ADJOINT mode (via SFAD) =================
  // The symbolic Jacobian DAG evaluates FFGradOCFESLV::eval<double>, which branches on
  // options.GRADIENT.  Both modes must give the same Jacobian, and match the FADType reference.
  {
    auto sJ = rdag.SFAD( F, pctrl );
    std::vector<unsigned> const& si = std::get<0>( sJ );
    std::vector<unsigned> const& sj = std::get<1>( sJ );
    std::vector<FFVar>    const& sd = std::get<2>( sJ );
    std::vector<double> vF( sd.size(), 0. ), vA( sd.size(), 0. );

    int const save = FFBaseOCFE::options.GRADIENT;
    FFBaseOCFE::options.GRADIENT = FFBaseOCFE::FORWARD; rdag.eval( sd, vF, pctrl, p0 );
    FFBaseOCFE::options.GRADIENT = FFBaseOCFE::ADJOINT; rdag.eval( sd, vA, pctrl, p0 );
    FFBaseOCFE::options.GRADIENT = save;

    double emode = 0., eref = 0.;
    for( size_t k=0; k<sd.size(); ++k ){
      emode = std::max( emode, std::fabs( vF[k] - vA[k] ) );
      if( si[k] < ncf && sj[k] < ncd )
        eref = std::max( eref, std::fabs( vF[k] - yF[si[k]].deriv((unsigned)sj[k]) ) );
    }
    check_close( "T7 SFAD FORWARD == ADJOINT mode", emode, 0., 1e-9 );
    check_close( "T7 SFAD FORWARD == FADType ref",  eref,  0., 1e-9 );
  }

  std::cout << "\n============================================================\n"
            << "  RESULT: " << g_pass << " passed, " << g_fail << " failed  -- "
            << ( g_fail == 0 ? "ALL PASS" : "SOME FAILED" ) << "\n"
            << "============================================================\n";
  return g_fail == 0 ? 0 : 1;
}
