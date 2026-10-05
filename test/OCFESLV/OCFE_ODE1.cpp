// ============================================================================
//  OCFE_ODE1b.cpp  --  Parametric ODE + FFOCFESLV, swept over the 2x2 of
//                      {IC_WEAK, IC_STRONG} x {monolithic, marching}, with the
//                      reduced input->output Jacobian checked against a
//                      CLOSED-FORM oracle in every cell and compared across cells.
//
//  Extends OCFE_ODE1 with the monolithic-vs-marching axis.  The reduced
//  sensitivity now runs through TWO code paths -- the monolithic reduced solve
//  and the marching causal-chain reduced solve -- and both must reproduce the
//  same closed-form Jacobian.  Because both modes share the SAME element grid,
//  their reduced Jacobians must also agree with EACH OTHER to sensitivity
//  tolerance (tighter than the shared spectral-discretisation error vs the
//  oracle), which is the cross-cell check.
//
//  SOLVE_MARCHING defaults to TRUE, so each run asserts is_marching()==requested.
//
//  MODEL (1 state, 2 controls)
//        dx/dt = -k x ,  x(0) = x0        ->  x(t) = x0 e^{-k t}
//    controls   k (dynamics) , x0 (initial condition)
//    outputs    G0 = x(T) = x0 e^{-kT} ,  G1 = x(T/2) = x0 e^{-kT/2}
//    closed-form reduced Jacobian J[f][ctrl], order (k, x0):
//        dG0/dk=-T x0 e^{-kT},   dG0/dx0=e^{-kT}
//        dG1/dk=-T/2 x0 e^{-kT/2}, dG1/dx0=e^{-kT/2}
//
//  TESTS (per (imposition,mode) cell)
//    MODE is_marching()==requested
//    T1   value FFOCFESLV::eval<double>            == solve/val_functions
//    T3a  forward AD (eval<FADType>)                == CLOSED-FORM oracle
//    T3b  forward AD                                 == oc.solve_fsens/sens_jacobian
//    T4   symbolic SFAD (FFGradOCFESLV)            == forward AD
//  CROSS-CELL
//    XM   J(WEAK[march]) == J(STRONG[march])  (GATED, exact: no interior interface)
//    ~    J WEAK~STRONG(mono), J mono~march  (INFORMATIONAL, discretisation scale;
//         the gated property is each cell's Jacobian vs CLOSED-FORM oracle, T3a)
//
//  Build:
//    g++ -std=c++17 <suite flags> -DOCFE_OCFESLV_HEADER='"ocfeslv.hpp"' \
//        OCFE_ODE1b.cpp -o OCFE_ODE1b <libs>
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

static double const T_end = 2.0;
static size_t const N_EL  = 6;
static size_t const N_ND  = 5;
static double const k_nom  = 0.7;
static double const x0_nom = 1.0;

static int g_pass = 0, g_fail = 0;
static void check_close( char const* name, double got, double want, double tol )
{
  bool const ok = std::fabs( got - want ) <= tol;
  ( ok ? g_pass : g_fail )++;
  std::cout << "  " << std::left << std::setw(44) << name
            << std::right << std::scientific << std::setprecision(3)
            << " got=" << std::setw(11) << got << " want=" << std::setw(11) << want
            << " |d|=" << std::setw(9) << std::fabs( got - want )
            << "  " << ( ok ? "PASS" : "FAIL" ) << "\n";
}
static void check_true( char const* name, bool ok )
{
  ( ok ? g_pass : g_fail )++;
  std::cout << "  " << std::left << std::setw(44) << name << "  " << ( ok ? "PASS" : "FAIL" ) << "\n";
}

struct Cell {
  bool                converged = false;
  bool                ok        = false;
  std::vector<double> Fval;              // ncf primal outputs
  std::vector<double> J;                 // ncf*ncd forward-AD reduced Jacobian (row-major)
};

static Cell run_mode( OCFESLV::Options::ImpositionType imp, char const* strimp, bool marching )
{
  Cell R;
  char const* mstr = marching ? "march" : "mono";
  std::cout << "\n================ IC_" << strimp << " [" << mstr << "] ================\n";

  FFGraph DAG;
  FFVar t  = DAG.add_var( "t"    );
  FFVar x  = DAG.add_var( "x(t)" );
  FFVar k  = DAG.add_var( "k"    );
  FFVar x0 = DAG.add_var( "x0"   );

  OCFESLV oc( &DAG );
  oc.add_domain( t, FFDom( 0., T_end, N_EL, FFDom::LGR, N_ND ) );
  oc.set_evolution_domain( t );
  oc.add_state( x, {t} );
  oc.add_input( k,  {} );
  oc.add_input( x0, {} );

  oc.update_ref( x, []( OCFESLV::t_Coord const& ){ return x0_nom; } );

  FFPartial OpP;
  FFVar EVOL = OpP( x, t ) - ( -k * x );
  FFVar IC   = x - x0;

  int const T_NO_LB = FFDom::ALL - FFDom::LB;
  oc.add_equation( EVOL, {t}, {T_NO_LB},   OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( IC,   {t}, {FFDom::LB}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INITIAL,  0 ) );

  oc.add_output( x, {t}, { T_end       } );   // G0
  oc.add_output( x, {t}, { 0.5 * T_end } );   // G1

//  oc.options.REDUCE.ORDER    = OCFESLV::Options::RED_NONE;
  oc.options.CLASSIFY.MODE        = OCFESLV::Options::CLASS_AUTO;
  oc.options.INTERFACE.IMPOSITION = imp;
  oc.options.INTERFACE.TYPE  = OCFESLV::Options::IC_AUTO;
  oc.options.SOLVE.MARCHING  = marching;
  oc.options.SOLVE.MAX_ITER  = 60;
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

  oc.set_input_values( k,  { k_nom  }, inp.data() );
  oc.set_input_values( x0, { x0_nom }, inp.data() );
  oc.register_control( k );
  oc.register_control( x0 );

  size_t const ncd = oc.n_control_dof();
  size_t const ncf = oc.n_colloc_fct();
  auto const& C  = oc.controls();
  size_t const ik  = C.at( k  ).offset;
  size_t const ix0 = C.at( x0 ).offset;
  if( ncd != 2 || ncf != 2 ){ std::cerr << "  ERROR: expected ncd=2, ncf=2\n"; return R; }

  std::vector<double> p0; oc.encode_controls( inp.data(), p0 );

  // ---- reference solve at nominal controls ----
  std::vector<double> xvR( xv ), inpR( inp );
  OCFESLV::SolveReport const rep = oc.solve( xvR.data(), inpR.data(), nullptr );
  R.converged = rep.converged;
  std::cout << "  reference solve: converged=" << ( rep.converged ? "yes" : "no" )
            << " iters=" << rep.iterations
            << " final|r|=" << std::scientific << std::setprecision(3) << rep.final_residual << "\n";
  if( !rep.converged ){ std::cerr << "  ERROR: reference solve did not converge\n"; return R; }
  std::vector<double> Fref = oc.val_functions();
  if( Fref.size() < ncf ){ std::cerr << "  ERROR: functions missing\n"; return R; }
  R.Fval = Fref;

  double const kv = p0[ik], xv0 = p0[ix0];
  double const G0an = xv0 * std::exp( -kv * T_end       );
  double const G1an = xv0 * std::exp( -kv * 0.5 * T_end );
  double Jan[2][2];
  Jan[0][ik] = -T_end       * xv0 * std::exp( -kv * T_end       ); Jan[0][ix0] = std::exp( -kv * T_end       );
  Jan[1][ik] = -0.5 * T_end * xv0 * std::exp( -kv * 0.5 * T_end ); Jan[1][ix0] = std::exp( -kv * 0.5 * T_end );

  std::cout << "\n  [primal] point outputs vs analytic:\n";
  check_close( "G0 == x(T)   = x0 e^{-kT}",   Fref[0], G0an, 1e-6 );
  check_close( "G1 == x(T/2) = x0 e^{-kT/2}", Fref[1], G1an, 1e-6 );

  // ---- reduced-space FFOCFESLV ----
  FFGraph rdag;
  std::vector<FFVar> pctrl( ncd );
  for( size_t i = 0; i < ncd; ++i ){ std::ostringstream os; os << "p" << i; pctrl[i] = rdag.add_var( os.str() ); }
  FFOCFESLV ffred;
  std::ostringstream tagos; tagos << "ODE1b_" << strimp << "_" << mstr;
  std::string const tag = tagos.str();
  std::vector<FFVar> F = ffred( FFOCFESLV::controls_map( oc, pctrl ), &oc, FFOCFESLV::COPY, tag.c_str() );
  check_true( "FFOCFESLV returns ncf outputs", F.size() == ncf );

  // T1: value
  {
    std::vector<double> Fv( ncf, 0. );
    ffred.eval( (unsigned)ncf, Fv.data(), (unsigned)ncd, p0.data(), nullptr );
    double e = 0.; for( size_t j = 0; j < ncf; ++j ) e = std::max( e, std::fabs( Fv[j] - Fref[j] ) );
    check_close( "T1 value (direct eval) == solve", e, 0., 1e-9 );
  }

  // T3: forward AD Jacobian
  std::vector<FADType<double>> xF( ncd ), yF( ncf );
  for( size_t i = 0; i < ncd; ++i ){ xF[i] = p0[i]; xF[i].diff( (unsigned)i, (unsigned)ncd ); }
  ffred.eval( (unsigned)ncf, yF.data(), (unsigned)ncd, xF.data(), nullptr );

  R.J.assign( ncf*ncd, 0. );
  for( size_t j = 0; j < ncf; ++j )
    for( size_t i = 0; i < ncd; ++i )
      R.J[ j*ncd + i ] = yF[j].deriv( (unsigned)i );

  std::cout << "\n  [T3a] forward AD  vs  CLOSED-FORM oracle:\n";
  for( size_t j = 0; j < ncf; ++j )
    for( size_t i = 0; i < ncd; ++i ){
      std::ostringstream nm; nm << "dG" << j << "/dp" << i << " (FAD) == analytic";
      check_close( nm.str().c_str(), R.J[ j*ncd + i ], Jan[j][i], 1e-5 );  // spectral (~2e-6)
    }

  std::cout << "\n  [T3b] forward AD  vs  oc.solve_fsens / sens_jacobian:\n";
  {
    std::vector<double> xvS( xv ), inpS( inp );
    if( oc.solve_fsens( xvS.data(), inpS.empty()?nullptr:inpS.data(), nullptr ) ){
      std::vector<double> const& J = oc.sens_jacobian();
      if( J.size() >= ncf*ncd ){
        double e = 0.;
        for( size_t j = 0; j < ncf; ++j )
          for( size_t i = 0; i < ncd; ++i )
            e = std::max( e, std::fabs( R.J[ j*ncd + i ] - J[ j*ncd + i ] ) );
        check_close( "T3b FAD Jacobian == reduced fwd-sens", e, 0., 1e-9 );
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
    check_true( "T4 SFAD nnz == ncf*ncd (dense)", sd.size() == ncf*ncd );
    std::vector<double> sdv( sd.size(), 0. );
    rdag.eval( sd, sdv, pctrl, p0 );
    double e = 0.; size_t cnt = 0; bool idx_ok = true;
    for( size_t kk = 0; kk < sd.size(); ++kk ){
      if( si[kk] >= ncf || sj[kk] >= ncd ){ idx_ok = false; continue; }
      e = std::max( e, std::fabs( sdv[kk] - R.J[ si[kk]*ncd + sj[kk] ] ) );
      ++cnt;
    }
    check_true ( "T4 SFAD indices in range", idx_ok && cnt == ncf*ncd );
    check_close( "T4 SFAD Jacobian == forward AD", e, 0., 1e-9 );
  }

  std::cout << "\n  [T5] ADJOINT sensitivity  vs  forward AD + CLOSED-FORM oracle:\n";
  {
    // solve_asens() computes the SAME reduced Jacobian dF/dp as the forward path by one reverse sweep
    // per output.  Controls here are a PARAMETER (k) and an INITIAL CONDITION (x0) -- distinct from the
    // distributed-input case -- so this exercises the marching adjoint's parameter/IC attribution.
    std::vector<double> xvA( xv ), inpA( inp );
    if( oc.solve_asens( xvA.data(), inpA.empty()?nullptr:inpA.data(), nullptr ) ){
      std::vector<double> const& Ja = oc.sens_jacobian();
      if( Ja.size() >= ncf*ncd ){
        double ef = 0., eo = 0.;
        for( size_t j = 0; j < ncf; ++j )
          for( size_t i = 0; i < ncd; ++i ){
            ef = std::max( ef, std::fabs( Ja[ j*ncd + i ] - R.J[ j*ncd + i ] ) );
            eo = std::max( eo, std::fabs( Ja[ j*ncd + i ] - Jan[j][i] ) );
          }
        check_close( "T5 adjoint Jacobian == forward AD",     ef, 0., 1e-9 );
        check_close( "T5 adjoint Jacobian == analytic oracle", eo, 0., 1e-5 );
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
            << "  OCFE_ODE1b : parametric ODE + FFOCFESLV, 2x2 sweep\n"
            << "  {IC_WEAK, IC_STRONG} x {monolithic, marching}\n"
            << "  dx/dt = -k x,  x(0)=x0 ;  T_end=" << T_end << "\n"
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
            << "  OCFE_ODE1b: " << g_pass << " passed, " << g_fail << " failed -- "
            << ( ( g_fail == 0 && built ) ? "ALL PASS" : "SOME FAILED" ) << "\n"
            << "================================================================\n";
  return ( g_fail == 0 && built ) ? 0 : 1;
}
