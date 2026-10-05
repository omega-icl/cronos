// ============================================================================
//  OCFE_DAE2.cpp  --  DAE with a PIECEWISE-CONSTANT time-grid control, wrapped
//                     in FFOCFESLV; evaluation and derivative accuracy w.r.t.
//                     the per-element control values, vs a CLOSED-FORM oracle.
//
//  Adds a genuinely time-varying control to the DAE family.  The control u(t)
//  is a DISTRIBUTED input carrying ONE Legendre-Gauss node per element
//  (add_input(u,{t},FFDom::LG,1)) -- i.e. piecewise-constant on the time grid,
//  jumping at element interfaces.  register_control(u) then exposes ONE control
//  DOF per element, so the reduced control vector p = (u_0,...,u_{Nel-1}) is the
//  vector of per-element control levels, and FFOCFESLV differentiates the
//  outputs with respect to each of them.
//
//  MODEL (1 differential + 1 algebraic state, Nel PWC control DOFs)
//        dc/dt = -a c + u(t) ,   c(0) = c0          (differential)
//        r - c^2 = 0                                 (algebraic; index-1)
//    outputs   G_m = c(t_m)  at each element boundary t_m = m h  (m=1..Nel)
//              G_last = r(T) = c(T)^2
//
//  Because u is constant = u_k on element k, the exact solution is the linear
//  recurrence  c_{k+1} = b c_k + (1-b) u_k / a ,  b = e^{-a h}, so the reduced
//  Jacobian is CLOSED-FORM and CAUSAL-TRIANGULAR (a later control cannot move an
//  earlier output):
//        d c(t_m)/d u_j = (1-b)/a * b^{m-1-j}   for j < m,   0 for j >= m
//        d r(T)/d u_j    = 2 c(T) * d c(T)/d u_j
//  The strict lower-triangular zero pattern (j>=m) is itself a test that the PWC
//  control is applied on the correct element and that causality is respected.
//
//  TESTS (per (imposition,mode) cell; run under IC_WEAK/IC_STRONG x mono/march)
//    MODE is_marching()==requested ;  ncd == Nel (one control DOF per element)
//    EV   c(t_m) == analytic recurrence                       [evaluation accuracy]
//    T1   FFOCFESLV::eval<double> == solve/val_functions
//    T3a  forward AD (eval<FADType>) == CLOSED-FORM oracle     [PWC derivative accuracy]
//    T3z  causal zeros: d c(t_m)/d u_j == 0 for j >= m         [triangular structure]
//    T3b  forward AD == oc.solve_fsens/sens_jacobian
//    T4   symbolic SFAD == forward AD
//  CROSS-CELL
//    XM   J(WEAK[march]) == J(STRONG[march]) GATED exact ;  J WEAK~STRONG(mono) &
//         J mono~march INFORMATIONAL (gated property = each cell vs oracle, T3a)
//
//  Build:
//    g++ -std=c++17 <suite flags> -DOCFE_OCFESLV_HEADER='"ocfeslv.hpp"' \
//        OCFE_DAE2.cpp -o OCFE_DAE2 <libs>
// ============================================================================

#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#ifndef OCFE_OCFESLV_HEADER
#define OCFE_OCFESLV_HEADER "ocfeslv.hpp"
#endif
#include OCFE_OCFESLV_HEADER
#include "ffocfe.hpp"

using namespace mc;

static double const A_dec = 1.0;                 // decay rate
static double const C0    = 0.5;                 // initial condition (fixed)
static size_t const N_EL  = 5;                   // elements == PWC control DOFs
static size_t const N_ND  = 7;                   // spectral nodes per element
static double const H_EL  = 1.0;                 // element width
static double const T_end = H_EL * double(N_EL); // = 5.0
static double const U_NOM[5] = { 1.0, 2.0, 0.5, 1.5, 1.0 };   // nominal per-element control

// ---- closed-form oracle ----
static void analytic_c( std::vector<double> const& uk, std::vector<double>& cm )
{
  double const b = std::exp( -A_dec * H_EL );
  cm.assign( uk.size() + 1, 0. );
  cm[0] = C0;
  for( size_t k = 0; k < uk.size(); ++k ) cm[k+1] = b*cm[k] + ( 1.0 - b )*uk[k]/A_dec;
}
static double dcm_duj( size_t m, size_t j )   // d c(t_m)/d u_j
{
  if( j >= m ) return 0.0;
  double const b = std::exp( -A_dec * H_EL );
  return ( 1.0 - b )/A_dec * std::pow( b, double( m - 1 - j ) );
}

static int g_pass = 0, g_fail = 0;
static void check_close( char const* name, double got, double want, double tol )
{
  bool const ok = std::fabs( got - want ) <= tol;
  ( ok ? g_pass : g_fail )++;
  std::cout << "  " << std::left << std::setw(40) << name
            << std::right << std::scientific << std::setprecision(3)
            << " got=" << std::setw(11) << got << " want=" << std::setw(11) << want
            << " |d|=" << std::setw(9) << std::fabs( got - want )
            << "  " << ( ok ? "PASS" : "FAIL" ) << "\n";
}
static void check_true( char const* name, bool ok )
{
  ( ok ? g_pass : g_fail )++;
  std::cout << "  " << std::left << std::setw(40) << name << "  " << ( ok ? "PASS" : "FAIL" ) << "\n";
}

struct Cell {
  bool                converged = false, ok = false;
  std::vector<double> Fval;               // ncf
  std::vector<double> J;                  // ncf*ncd forward-AD Jacobian
};

static Cell run_mode( OCFESLV::Options::ImpositionType imp, char const* strimp, bool marching )
{
  Cell R;
  char const* mstr = marching ? "march" : "mono";
  std::cout << "\n================ IC_" << strimp << " [" << mstr << "] ================\n";

  FFGraph DAG;
  FFVar t = DAG.add_var( "t"    );
  FFVar c = DAG.add_var( "c(t)" );
  FFVar r = DAG.add_var( "r(t)" );
  FFVar u = DAG.add_var( "u(t)" );      // PWC control (distributed on t)
  FFPartial OpP;

  OCFESLV oc( &DAG );
  oc.add_domain( t, FFDom( 0., T_end, N_EL, FFDom::LGR, N_ND ) );
  oc.set_evolution_domain( t );
  oc.add_state( c, {t} );
  oc.add_state( r, {t} );
  // ONE Legendre-Gauss node per element => piecewise-constant control, Nel DOFs.
  oc.add_input( u, {t}, FFDom::LG, 1, std::optional<double>( 1.0 ) );

  oc.update_ref( c, []( OCFESLV::t_Coord const& ){ return C0; } );
  oc.update_ref( r, []( OCFESLV::t_Coord const& ){ return C0*C0; } );

  FFVar EVOL = OpP( c, t ) - ( -A_dec * c + u );
  FFVar ALG  = r - c*c;
  FFVar IC   = c - C0;

  int const T_NO_LB = FFDom::ALL - FFDom::LB;
  oc.add_equation( EVOL, {t}, {T_NO_LB},    OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( ALG,  {t}, {FFDom::ALL}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( IC,   {t}, {FFDom::LB},  OCFESLV::EqnOptions( OCFESLV::EqnRole::INITIAL,  0 ) );

  for( size_t m = 1; m <= N_EL; ++m ) oc.add_output( c, {t}, { double(m)*H_EL } );  // c(t_m)
  oc.add_output( r, {t}, { T_end } );                                                // r(T)

//  oc.options.REDUCE.ORDER    = OCFESLV::Options::RED_NONE;
  oc.options.CLASSIFY.MODE        = OCFESLV::Options::CLASS_AUTO;
  oc.options.INTERFACE.IMPOSITION = imp;
  oc.options.INTERFACE.TYPE  = OCFESLV::Options::IC_AUTO;
  oc.options.SOLVE.MARCHING  = marching;
  oc.options.SOLVE.MAX_ITER  = 80;
  oc.options.SOLVE.RES_TOL   = 1.0e-10;
  oc.options.DISPLAY_LEVEL   = 1;

  if( !oc.setup() ){
    std::cerr << "  setup() FAILED: " << OCFESLV::setup_status_str( oc.setup_status() ) << "\n";
    return R;
  }

  bool const marching_active = oc.is_marching();
  std::cout << "  is_marching()=" << ( marching_active ? "true" : "false" )
            << "  n_march_steps=" << oc.n_march_steps() << "\n";
  check_true( "MODE is_marching()==requested", marching_active == marching );

  std::vector<double> xv, inp;
  if( !oc.init( xv, inp, nullptr ) ){ std::cerr << "  init() FAILED\n"; return R; }

  // set the nominal per-element control levels, then register as a control
  std::vector<double> uvec( U_NOM, U_NOM + N_EL );
  oc.set_input_values( u, uvec, inp.data() );
  oc.register_control( u );

  size_t const ncd = oc.n_control_dof();
  size_t const ncf = oc.n_colloc_fct();
  std::cout << "  states=" << oc.n_colloc_sta() << " eqns=" << oc.n_colloc_eqn()
            << " control DOFs ncd=" << ncd << " (expect Nel=" << N_EL << ")"
            << " outputs ncf=" << ncf << "\n";
  check_true( "ncd == Nel (one control DOF per element)", ncd == N_EL );
  if( ncd != N_EL || ncf != N_EL + 1 ){ std::cerr << "  ERROR: expected ncd=Nel, ncf=Nel+1\n"; return R; }

  std::vector<double> p0; oc.encode_controls( inp.data(), p0 );

  // ---- reference solve ----
  std::vector<double> xvR( xv ), inpR( inp );
  OCFESLV::SolveReport const rep = oc.solve( xvR.data(), inpR.data(), nullptr );
  R.converged = rep.converged;
  std::cout << "  reference solve: converged=" << ( rep.converged ? "yes" : "no" )
            << " iters=" << rep.iterations
            << " final|r|=" << std::scientific << std::setprecision(3) << rep.final_residual << "\n";
  if( !rep.converged ){ std::cerr << "  ERROR: reference solve did not converge\n"; return R; }
  R.Fval = oc.val_functions();
  if( R.Fval.size() < ncf ){ std::cerr << "  ERROR: functions missing\n"; return R; }

  // analytic recurrence at the nominal control
  std::vector<double> cm; analytic_c( uvec, cm );   // cm[0..Nel]

  std::cout << "\n  [EV] evaluation accuracy: c(t_m) vs analytic recurrence:\n";
  for( size_t m = 1; m <= N_EL; ++m ){
    std::ostringstream nm; nm << "EV c(t" << m << ") == recurrence";
    check_close( nm.str().c_str(), R.Fval[m-1], cm[m], 5e-6 );  // IC_WEAK not exact (~disc scale)
  }
  check_close( "EV r(T) == c(T)^2", R.Fval[N_EL], cm[N_EL]*cm[N_EL], 5e-6 );

  // closed-form reduced Jacobian oracle Jan[f][j]
  std::vector<std::vector<double>> Jan( ncf, std::vector<double>( ncd, 0. ) );
  for( size_t m = 1; m <= N_EL; ++m )
    for( size_t j = 0; j < ncd; ++j )
      Jan[m-1][j] = dcm_duj( m, j );
  for( size_t j = 0; j < ncd; ++j )
    Jan[N_EL][j] = 2.0 * cm[N_EL] * dcm_duj( N_EL, j );

  // ---- FFOCFESLV ----
  FFGraph rdag;
  std::vector<FFVar> pctrl( ncd );
  for( size_t i = 0; i < ncd; ++i ){ std::ostringstream os; os << "p" << i; pctrl[i] = rdag.add_var( os.str() ); }
  FFOCFESLV ffred;
  std::ostringstream tagos; tagos << "DAE2_" << strimp << "_" << mstr;
  std::string const tag = tagos.str();
  std::vector<FFVar> F = ffred( FFOCFESLV::controls_map( oc, pctrl ), &oc, FFOCFESLV::COPY, tag.c_str() );
  check_true( "FFOCFESLV returns ncf outputs", F.size() == ncf );

  // T1: value
  {
    std::vector<double> Fv( ncf, 0. );
    ffred.eval( (unsigned)ncf, Fv.data(), (unsigned)ncd, p0.data(), nullptr );
    double e = 0.; for( size_t j = 0; j < ncf; ++j ) e = std::max( e, std::fabs( Fv[j] - R.Fval[j] ) );
    check_close( "T1 value (direct eval) == solve", e, 0., 1e-9 );
  }

  // T3: forward AD Jacobian
  std::vector<FADType<double>> xF( ncd ), yF( ncf );
  for( size_t i = 0; i < ncd; ++i ){ xF[i] = p0[i]; xF[i].diff( (unsigned)i, (unsigned)ncd ); }
  ffred.eval( (unsigned)ncf, yF.data(), (unsigned)ncd, xF.data(), nullptr );
  R.J.assign( ncf*ncd, 0. );
  for( size_t f = 0; f < ncf; ++f )
    for( size_t j = 0; j < ncd; ++j )
      R.J[ f*ncd + j ] = yF[f].deriv( (unsigned)j );

  std::cout << "\n  [T3a] forward AD  vs  CLOSED-FORM oracle (PWC control):\n";
  { double e = 0.;
    for( size_t f = 0; f < ncf; ++f )
      for( size_t j = 0; j < ncd; ++j )
        e = std::max( e, std::fabs( R.J[ f*ncd + j ] - Jan[f][j] ) );
    check_close( "T3a reduced Jacobian == analytic", e, 0., 1e-5 );
    // DIAGNOSTIC: dump the reduced Jacobian (got vs analytic) to localize the marched-sensitivity error.
    std::cout << "  [T3a-DIAG] reduced Jacobian rows f=0.." << (ncf-1) << " cols j=0.." << (ncd-1) << "\n";
    for( size_t f = 0; f < ncf; ++f ){
      std::cout << "    got[" << f << "]";
      for( size_t j = 0; j < ncd; ++j ) std::cout << " " << std::setw(11) << std::setprecision(4) << R.J[f*ncd+j];
      std::cout << "   | want";
      for( size_t j = 0; j < ncd; ++j ) std::cout << " " << std::setw(11) << std::setprecision(4) << Jan[f][j];
      std::cout << "\n";
    }
  }

  std::cout << "\n  [T3z] causal-triangular zero pattern (j>=m => 0):\n";
  { double zmax = 0.;
    for( size_t m = 1; m <= N_EL; ++m )
      for( size_t j = m; j < ncd; ++j )                 // j >= m: later control, earlier output
        zmax = std::max( zmax, std::fabs( R.J[ (m-1)*ncd + j ] ) );
    check_close( "T3z d c(t_m)/d u_j == 0 for j>=m", zmax, 0., 1e-5 );  // IC_WEAK SAT couples both ways; STRONG/march ~exact
  }

  std::cout << "\n  [T3b] forward AD  vs  oc.solve_fsens / sens_jacobian:\n";
  {
    std::vector<double> xvS( xv ), inpS( inp );
    if( oc.solve_fsens( xvS.data(), inpS.empty()?nullptr:inpS.data(), nullptr ) ){
      std::vector<double> const& J = oc.sens_jacobian();
      if( J.size() >= ncf*ncd ){
        double e = 0.;
        for( size_t f = 0; f < ncf; ++f )
          for( size_t j = 0; j < ncd; ++j )
            e = std::max( e, std::fabs( R.J[ f*ncd + j ] - J[ f*ncd + j ] ) );
        check_close( "T3b FAD Jacobian == reduced fwd-sens", e, 0., 1e-8 );
      }
      else check_true( "T3b sens_jacobian sized ncf*ncd", false );
    }
    else check_true( "T3b solve_fsens available", false );
  }

  std::cout << "\n  [T4] symbolic SFAD Jacobian  vs  forward AD:\n";
  {
    auto sJac = rdag.SFAD( F, pctrl );
    std::vector<unsigned> const& si = std::get<0>( sJac );
    std::vector<unsigned> const& sj = std::get<1>( sJac );
    std::vector<FFVar>    const& sd = std::get<2>( sJac );
    std::vector<double> sdv( sd.size(), 0. );
    rdag.eval( sd, sdv, pctrl, p0 );
    double e = 0.; bool idx_ok = true;
    for( size_t kk = 0; kk < sd.size(); ++kk ){
      if( si[kk] >= ncf || sj[kk] >= ncd ){ idx_ok = false; continue; }
      e = std::max( e, std::fabs( sdv[kk] - R.J[ si[kk]*ncd + sj[kk] ] ) );
    }
    check_true ( "T4 SFAD indices in range", idx_ok );
    check_close( "T4 SFAD Jacobian == forward AD", e, 0., 1e-8 );
  }

  std::cout << "\n  [T5] ADJOINT sensitivity  vs  forward AD + CLOSED-FORM oracle:\n";
  {
    // solve_asens() computes the SAME reduced Jacobian dF/dp (n_colloc_fct x n_control_dof, row-major)
    // as the forward path, but by one reverse sweep per output.  It must agree with BOTH the forward AD
    // Jacobian R.J (transpose-consistency: forward and adjoint compute the same dF/dp) and the analytic
    // oracle Jan.  Exercised in every cell -> {IC_WEAK,IC_STRONG} x {mono,march}.
    std::vector<double> xvA( xv ), inpA( inp );
    if( oc.solve_asens( xvA.data(), inpA.empty()?nullptr:inpA.data(), nullptr ) ){
      std::vector<double> const& Ja = oc.sens_jacobian();
      if( Ja.size() >= ncf*ncd ){
        double ef = 0., eo = 0.;
        for( size_t f = 0; f < ncf; ++f )
          for( size_t j = 0; j < ncd; ++j ){
            ef = std::max( ef, std::fabs( Ja[ f*ncd + j ] - R.J[ f*ncd + j ] ) );   // adjoint == forward AD
            eo = std::max( eo, std::fabs( Ja[ f*ncd + j ] - Jan[f][j] ) );          // adjoint == oracle
          }
        check_close( "T5 adjoint Jacobian == forward AD",     ef, 0., 1e-8 );
        check_close( "T5 adjoint Jacobian == analytic oracle", eo, 0., 1e-5 );
        // [T5-DIAG] full adjoint Jacobian vs forward AD -- localise any marching-adjoint mismatch
        std::cout << "  [T5-DIAG] adjoint Ja (left) vs forward R.J (right), rows f=0.." << ncf-1 << "\n";
        for( size_t f = 0; f < ncf; ++f ){
          std::cout << "    adj[" << f << "] ";
          for( size_t j = 0; j < ncd; ++j ) std::cout << std::setw(11) << std::setprecision(4) << Ja[f*ncd+j];
          std::cout << "   | fwd ";
          for( size_t j = 0; j < ncd; ++j ) std::cout << std::setw(11) << std::setprecision(4) << R.J[f*ncd+j];
          std::cout << "\n";
        }
      }
      else check_true( "T5 sens_jacobian sized ncf*ncd", false );
    }
    else check_true( "T5 solve_asens available", false );
  }

  R.ok = true;
  return R;
}

static void compare_J( char const* name, Cell const& a, Cell const& b, double tol )
{
  if( a.J.empty() || a.J.size() != b.J.size() ){ check_true( name, false ); return; }
  double e = 0.; for( size_t i = 0; i < a.J.size(); ++i ) e = std::max( e, std::fabs( a.J[i] - b.J[i] ) );
  check_close( name, e, 0., tol );
}

// informational cross-cell Jacobian delta: WEAK-vs-STRONG(mono) and mono-vs-march carry the
// same discretisation-scale difference as the primal, so they are printed, not gated; the
// gated property is each cell's Jacobian vs the CLOSED-FORM oracle (T3a) inside run_mode.
static void info_J( char const* name, Cell const& a, Cell const& b )
{
  double e = 0.;
  if( !a.J.empty() && a.J.size() == b.J.size() )
    for( size_t i = 0; i < a.J.size(); ++i ) e = std::max( e, std::fabs( a.J[i] - b.J[i] ) );
  std::cout << "  " << std::left << std::setw(44) << name
            << " delta=" << std::right << std::scientific << std::setprecision(3) << e
            << "  (discretisation-scale; informational)\n";
}

int main()
{
  std::cout << "================================================================\n"
            << "  OCFE_DAE2 : DAE with PIECEWISE-CONSTANT time-grid control\n"
            << "  dc/dt = -a c + u(t),  r = c^2 ;  u PWC (Nel=" << N_EL << " DOFs)\n"
            << "  evaluation + derivative accuracy vs closed-form triangular oracle\n"
            << "  swept over {IC_WEAK, IC_STRONG} x {monolithic, marching}\n"
            << "================================================================\n";

  Cell const wm = run_mode( OCFESLV::Options::IC_WEAK,   "WEAK",   false );
  Cell const sm = run_mode( OCFESLV::Options::IC_STRONG, "STRONG", false );
  Cell const wc = run_mode( OCFESLV::Options::IC_WEAK,   "WEAK",   true  );
  Cell const sc = run_mode( OCFESLV::Options::IC_STRONG, "STRONG", true  );

  std::cout << "\n---------------- cross-cell Jacobian comparisons ----------------\n";
  // Gated only where exact: marching windows have no interior interface, so the WEAK and
  // STRONG reduced Jacobians are the same computation.
  compare_J( "XM  J(WEAK[march]) == J(STRONG[march]) (exact)", wc, sc, 1e-9 );
  // WEAK-vs-STRONG(mono) and mono-vs-march differ at the discretisation scale; the gated
  // property is each cell's Jacobian vs the CLOSED-FORM oracle (T3a), inside run_mode.
  info_J( "XM  J(WEAK[mono])   ~ J(STRONG[mono])", wm, sm );
  info_J( "MM  J(WEAK):   mono ~ march",           wm, wc );
  info_J( "MM  J(STRONG): mono ~ march",           sm, sc );

  bool built = wm.ok && sm.ok && wc.ok && sc.ok;
  std::cout << "\n================================================================\n"
            << "  OCFE_DAE2: " << g_pass << " passed, " << g_fail << " failed -- "
            << ( ( g_fail == 0 && built ) ? "ALL PASS" : "SOME FAILED" ) << "\n"
            << "================================================================\n";
  return ( g_fail == 0 && built ) ? 0 : 1;
}
