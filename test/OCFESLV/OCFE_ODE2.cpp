// ============================================================================
//  OCFE_ODE2.cpp  --  Practical NONLINEAR batch-reactor ODE + FFOCFESLV, with
//                     the reduced Jacobian checked against a CLOSED-FORM oracle,
//                     swept over {IC_WEAK, IC_STRONG} x {monolithic, marching}.
//
//  A step up from OCFE_ODE1 (scalar linear decay): a genuinely NONLINEAR,
//  multi-state, physically meaningful system that nonetheless retains a
//  closed-form solution and Jacobian, so the derivative test stays a HARD oracle
//  (not just FD).  Complementary to the CSTR tests (those are a continuous
//  index-1 DAE validated by forward/adjoint + FD; this is a transient BATCH ODE
//  validated by a closed-form reduced Jacobian and an exact invariant).
//
//  MODEL  --  isothermal 2nd-order batch reaction  A -> C
//        dCA/dt = -k CA^2 ,  CA(0) = CA0
//        dCC/dt = +k CA^2 ,  CC(0) = 0
//    controls   k (dynamics) ,  CA0 (initial condition)
//    outputs    G0 = CA(T) ,  G1 = CC(T) ,  G2 = CA(T)+CC(T)   (mass; invariant)
//
//  CLOSED FORM (equimolar single reactant, exact):
//        CA(t) = CA0 / D(t) ,   D(t) = 1 + k CA0 t ,   CC(t) = CA0 - CA(t)
//    exact mass invariant  CA(t)+CC(t) = CA0  for all t.
//
//  CLOSED-FORM reduced Jacobian J[f][ctrl], control order (k, CA0), at t=T,
//  D = 1 + k CA0 T:
//        dG0/dk  = -CA0^2 T / D^2      dG0/dCA0 = 1/D^2
//        dG1/dk  = +CA0^2 T / D^2      dG1/dCA0 = 1 - 1/D^2
//        dG2/dk  = 0                   dG2/dCA0 = 1        (conservation sensitivity)
//
//  TESTS (per (imposition,mode) cell)
//    MODE is_marching()==requested
//    EV   CA(T),CC(T) == closed form ;  mass G2 == CA0            [invariant]
//    T1   FFOCFESLV::eval<double> == solve/val_functions
//    T3a  forward AD (eval<FADType>) == CLOSED-FORM oracle
//    T3b  forward AD == oc.solve_fsens/sens_jacobian
//    T4   symbolic SFAD == forward AD
//  CROSS-CELL
//    XM   J(WEAK[march]) == J(STRONG[march]) GATED exact ;  J WEAK~STRONG(mono) &
//         J mono~march INFORMATIONAL (gated property = each cell vs oracle, T3a)
//
//  Build:
//    g++ -std=c++17 <suite flags> -DOCFE_OCFESLV_HEADER='"ocfeslv.hpp"' \
//        OCFE_ODE2.cpp -o OCFE_ODE2 <libs>
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

static double const T_end   = 2.0;
static size_t const N_EL    = 6;
static size_t const N_ND    = 6;
static double const k_nom   = 1.0;   // rate constant
static double const CA0_nom = 1.0;   // initial reactant

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
  bool                converged = false, ok = false;
  std::vector<double> Fval;
  std::vector<double> J;
};

static Cell run_mode( OCFESLV::Options::ImpositionType imp, char const* strimp, bool marching )
{
  Cell R;
  char const* mstr = marching ? "march" : "mono";
  std::cout << "\n================ IC_" << strimp << " [" << mstr << "] ================\n";

  FFGraph DAG;
  FFVar t   = DAG.add_var( "t"     );
  FFVar CA  = DAG.add_var( "CA(t)" );
  FFVar CC  = DAG.add_var( "CC(t)" );
  FFVar k   = DAG.add_var( "k"     );      // control: dynamics
  FFVar CA0 = DAG.add_var( "CA0"   );      // control: initial condition
  FFPartial OpP;

  OCFESLV oc( &DAG );
  oc.add_domain( t, FFDom( 0., T_end, N_EL, FFDom::LGR, N_ND ) );   // Radau: natural for evolution
  oc.set_evolution_domain( t );
  oc.add_state( CA, {t} );
  oc.add_state( CC, {t} );
  oc.add_input( k,   {} );
  oc.add_input( CA0, {} );

  oc.update_ref( CA, []( OCFESLV::t_Coord const& ){ return CA0_nom; } );
  oc.update_ref( CC, []( OCFESLV::t_Coord const& ){ return 0.0;     } );

  FFVar EVOL_CA = OpP( CA, t ) + k * CA * CA;      // dCA/dt = -k CA^2
  FFVar EVOL_CC = OpP( CC, t ) - k * CA * CA;      // dCC/dt = +k CA^2
  FFVar IC_CA   = CA - CA0;
  FFVar IC_CC   = CC - 0.0;

  int const T_NO_LB = FFDom::ALL - FFDom::LB;
  OCFESLV::EqnOptions const io( OCFESLV::EqnRole::INTERIOR, 0 );
  OCFESLV::EqnOptions const ii( OCFESLV::EqnRole::INITIAL,  0 );
  oc.add_equation( EVOL_CA, {t}, {T_NO_LB},   io );
  oc.add_equation( EVOL_CC, {t}, {T_NO_LB},   io );
  oc.add_equation( IC_CA,   {t}, {FFDom::LB}, ii );
  oc.add_equation( IC_CC,   {t}, {FFDom::LB}, ii );

  FFVar OUT_MASS = CA + CC;
  oc.add_output( CA,       {t}, { T_end } );   // G0
  oc.add_output( CC,       {t}, { T_end } );   // G1
  oc.add_output( OUT_MASS, {t}, { T_end } );   // G2 = mass (invariant)

  oc.options.REDUCE.ORDER    = OCFESLV::Options::RED_NONE;
  oc.options.CLASSIFY.MODE        = OCFESLV::Options::CLASS_AUTO;
  oc.options.INTERFACE.IMPOSITION = imp;
  oc.options.INTERFACE.TYPE  = OCFESLV::Options::IC_AUTO;
  oc.options.SOLVE.MARCHING  = marching;
  oc.options.SOLVE.MAX_ITER  = 80;
  oc.options.SOLVE.RES_TOL   = 1.0e-10;
  oc.options.DISPLAY_LEVEL   = 0;

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

  oc.set_input_values( k,   { k_nom   }, inp.data() );
  oc.set_input_values( CA0, { CA0_nom }, inp.data() );
  oc.register_control( k );
  oc.register_control( CA0 );

  size_t const ncd = oc.n_control_dof();
  size_t const ncf = oc.n_colloc_fct();
  auto const& C  = oc.controls();
  size_t const ik  = C.at( k   ).offset;
  size_t const iC0 = C.at( CA0 ).offset;
  std::cout << "  states=" << oc.n_colloc_sta() << " eqns=" << oc.n_colloc_eqn()
            << " controls ncd=" << ncd << " (k@" << ik << ", CA0@" << iC0 << ")"
            << " outputs ncf=" << ncf << "\n";
  if( ncd != 2 || ncf != 3 ){ std::cerr << "  ERROR: expected ncd=2, ncf=3\n"; return R; }

  std::vector<double> p0; oc.encode_controls( inp.data(), p0 );

  std::vector<double> xvR( xv ), inpR( inp );
  OCFESLV::SolveReport const rep = oc.solve( xvR.data(), inpR.data(), nullptr );
  R.converged = rep.converged;
  std::cout << "  reference solve: converged=" << ( rep.converged ? "yes" : "no" )
            << " iters=" << rep.iterations
            << " final|r|=" << std::scientific << std::setprecision(3) << rep.final_residual << "\n";
  if( !rep.converged ){ std::cerr << "  ERROR: reference solve did not converge\n"; return R; }
  R.Fval = oc.val_functions();
  if( R.Fval.size() < ncf ){ std::cerr << "  ERROR: functions missing\n"; return R; }

  // closed form at nominal
  double const kv = p0[ik], c0 = p0[iC0];
  double const D  = 1.0 + kv * c0 * T_end;
  double const CA_T = c0 / D, CC_T = c0 - CA_T;
  double Jan[3][2];
  Jan[0][ik] = -c0*c0*T_end/(D*D); Jan[0][iC0] = 1.0/(D*D);
  Jan[1][ik] =  c0*c0*T_end/(D*D); Jan[1][iC0] = 1.0 - 1.0/(D*D);
  Jan[2][ik] =  0.0;               Jan[2][iC0] = 1.0;

  std::cout << "\n  [EV] evaluation accuracy vs closed form:\n";
  check_close( "EV CA(T) == CA0/(1+k CA0 T)", R.Fval[0], CA_T, 1e-6 );
  check_close( "EV CC(T) == CA0 - CA(T)",     R.Fval[1], CC_T, 1e-6 );
  check_close( "EV mass CA(T)+CC(T) == CA0",  R.Fval[2], c0,   5e-6 );  // IC_WEAK: SAT, not exact (~3e-7)

  // ---- FFOCFESLV ----
  FFGraph rdag;
  std::vector<FFVar> pctrl( ncd );
  for( size_t i = 0; i < ncd; ++i ){ std::ostringstream os; os << "p" << i; pctrl[i] = rdag.add_var( os.str() ); }
  FFOCFESLV ffred;
  std::ostringstream tagos; tagos << "ODE2_" << strimp << "_" << mstr;
  std::string const tag = tagos.str();
  std::vector<FFVar> F = ffred( FFOCFESLV::controls_map( oc, pctrl ), &oc, FFOCFESLV::COPY, tag.c_str() );
  check_true( "FFOCFESLV returns ncf outputs", F.size() == ncf );

  {
    std::vector<double> Fv( ncf, 0. );
    ffred.eval( (unsigned)ncf, Fv.data(), (unsigned)ncd, p0.data(), nullptr );
    double e = 0.; for( size_t j = 0; j < ncf; ++j ) e = std::max( e, std::fabs( Fv[j] - R.Fval[j] ) );
    check_close( "T1 value (direct eval) == solve", e, 0., 1e-9 );
  }

  std::vector<FADType<double>> xF( ncd ), yF( ncf );
  for( size_t i = 0; i < ncd; ++i ){ xF[i] = p0[i]; xF[i].diff( (unsigned)i, (unsigned)ncd ); }
  ffred.eval( (unsigned)ncf, yF.data(), (unsigned)ncd, xF.data(), nullptr );
  R.J.assign( ncf*ncd, 0. );
  for( size_t f = 0; f < ncf; ++f )
    for( size_t j = 0; j < ncd; ++j )
      R.J[ f*ncd + j ] = yF[f].deriv( (unsigned)j );

  std::cout << "\n  [T3a] forward AD  vs  CLOSED-FORM oracle:\n";
  for( size_t f = 0; f < ncf; ++f )
    for( size_t j = 0; j < ncd; ++j ){
      std::ostringstream nm; nm << "dG" << f << "/dp" << j << " (FAD) == analytic";
      check_close( nm.str().c_str(), R.J[ f*ncd + j ], Jan[f][j], 1e-5 );
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
    // solve_asens() reduced dF/dp via reverse sweeps; controls are a rate parameter (k) and a feed
    // concentration (CA0).  Asserted equal to the forward AD Jacobian and the closed-form oracle, in
    // every cell -> {IC_WEAK,IC_STRONG} x {mono,march}.
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
        check_close( "T5 adjoint Jacobian == forward AD",     ef, 0., 1e-8 );
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
            << "  OCFE_ODE2 : nonlinear 2nd-order batch reactor A->C + FFOCFESLV\n"
            << "  dCA/dt = -k CA^2, dCC/dt = +k CA^2 ;  closed-form Jacobian oracle\n"
            << "  swept over {IC_WEAK, IC_STRONG} x {monolithic, marching}\n"
            << "  T_end=" << T_end << "  (LGR time discretisation)\n"
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
            << "  OCFE_ODE2: " << g_pass << " passed, " << g_fail << " failed -- "
            << ( ( g_fail == 0 && built ) ? "ALL PASS" : "SOME FAILED" ) << "\n"
            << "================================================================\n";
  return ( g_fail == 0 && built ) ? 0 : 1;
}
