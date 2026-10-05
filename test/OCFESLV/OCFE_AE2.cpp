// ============================================================================
//  OCFE_AE2.cpp  --  OCFESLV / FFOCFESLV on a two-CSTR STEADY-STATE model.
//
//  Reuses the reaction model of test7.cpp (two continuous reactors in series,
//  A + D -> B, B -> C, exothermic).  test7 integrates the ODE to steady state;
//  here the SAME model is posed as the STEADY-STATE ALGEBRAIC system RHS(X,U)=0
//  (at steady state dX/dt=0, so the ODE's steady state is the root of RHS).  This
//  keeps AE2 in the zero-domain algebraic family of OCFE_AE / OCFE_AE1 while being
//  a strongly-coupled, Arrhenius-nonlinear model -- a genuine test of FFOCFESLV
//  evaluation and derivatives on a realistic system.
//
//    states (10)  reactor 1: CA1 CB1 CC1 CD1 T1 ;  reactor 2: CA2 CB2 CC2 CD2 T2
//    inputs  (6)  U = [ CA0, CA0_CD0, q0, V1, V2, T0 ]     (registered controls)
//    params  (8)  fixed, baked in as constants (dP)
//    outputs (3)  Y0 = CA2/(CA2+CB2+CC2+CD2)
//                 Y1 = CC2/(CA2+CB2+CC2+CD2)
//                 Y2 = T2
//
//  Per reactor r with inlet (CAi,CBi,CCi,CDi,Ti), volume Vr:
//    k1 = k0_1 exp(-Ea_1/(T+273.15)),  k2 = k0_2 exp(-Ea_2/(T+273.15))
//    R1 = k1 CA CD,  R2 = k2 CB
//    0 = q0 (CAi-CA) + Vr(-R1)
//    0 = q0 (CBi-CB) + Vr( R1-R2)
//    0 = q0 (CCi-CC) + Vr( R2)
//    0 = q0 (CDi-CD) + Vr(-R1)
//    0 = q0 rho cp (Ti-T) + Vr(-R1 DH1 - R2 DH2)
//  with reactor-1 inlet = feed (CA0, 0, 0, CA0/CA0_CD0, T0) and reactor-2 inlet =
//  reactor-1 outlet.
//
//  TESTS
//    T1   value      FFOCFESLV::eval<double>          == direct oc.solve + val_functions
//    Tres residual   F0..F9 at the solution (DAG.eval)  == 0  (proves a genuine root of RHS=0)
//    Tphy physical    all concentrations >= 0  (the system is multi-rooted; a negative-
//                     concentration root also closes the residuals, so guard the physical one)
//    T3a  d/dU       forward AD (eval<FADType>)        == oc.solve_fsens/sens_jacobian
//    T3b  d/dU       forward AD                         == central differences (ground truth)
//    T4   d/dU       symbolic SFAD (FFGradOCFESLV)     == forward AD
//
//  Build:
//    g++ -std=c++17 <suite flags> -DOCFE_OCFESLV_HEADER='"ocfeslv.hpp"' \
//        OCFE_AE2.cpp -o OCFE_AE2 <libs>
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

static int g_pass = 0, g_fail = 0;
static void check_close( char const* name, double got, double want, double tol )
{
  bool const ok = std::fabs( got - want ) <= tol;
  ( ok ? g_pass : g_fail )++;
  std::cout << "  " << std::left << std::setw(40) << name
            << std::right << std::scientific << std::setprecision(4)
            << " got=" << std::setw(12) << got << " want=" << std::setw(12) << want
            << " |d|=" << std::setw(10) << std::fabs( got - want )
            << "  " << ( ok ? "PASS" : "FAIL" ) << "\n";
}
static void check_true( char const* name, bool ok )
{
  ( ok ? g_pass : g_fail )++;
  std::cout << "  " << std::left << std::setw(40) << name << "  " << ( ok ? "PASS" : "FAIL" ) << "\n";
}

int main()
{
  std::cout << "================================================================\n"
            << "  OCFE_AE2 : two-CSTR steady-state (algebraic) via OCFESLV/FFOCFESLV\n"
            << "  10 states, 6 input controls, 3 outputs; eval + derivatives\n"
            << "================================================================\n";

  // ---- fixed model parameters (test7 dP), baked in as constants ----
  double const k0_1 = 2800., k0_2 = 12., Ea_1 = 2995., Ea_2 = 4427.;
  double const DH1  = -80.,  DH2  = 0.,  cp   = 1.7,   rho  = 0.8;
  double const CB0  = 0.,    CC0  = 0.;

  FFGraph DAG;
  // states
  FFVar CA1 = DAG.add_var("CA1"), CB1 = DAG.add_var("CB1"), CC1 = DAG.add_var("CC1"),
        CD1 = DAG.add_var("CD1"), T1  = DAG.add_var("T1");
  FFVar CA2 = DAG.add_var("CA2"), CB2 = DAG.add_var("CB2"), CC2 = DAG.add_var("CC2"),
        CD2 = DAG.add_var("CD2"), T2  = DAG.add_var("T2");
  // inputs (controls)
  FFVar CA0     = DAG.add_var("CA0");
  FFVar CA0_CD0 = DAG.add_var("CA0_CD0");
  FFVar q0      = DAG.add_var("q0");
  FFVar V1      = DAG.add_var("V1");
  FFVar V2      = DAG.add_var("V2");
  FFVar T0      = DAG.add_var("T0");

  FFVar CD0 = CA0 / CA0_CD0;   // feed of D

  // ---- steady-state residuals RHS = 0 ----
  // reactor 1 (inlet = feed)
  FFVar k11 = k0_1 * exp( -Ea_1/( T1 + 273.15 ) );
  FFVar k21 = k0_2 * exp( -Ea_2/( T1 + 273.15 ) );
  FFVar R11 = k11 * CA1 * CD1;
  FFVar R21 = k21 * CB1;
  FFVar F0 = q0*( CA0 - CA1 ) + V1*( -R11 );
  FFVar F1 = q0*( CB0 - CB1 ) + V1*(  R11 - R21 );
  FFVar F2 = q0*( CC0 - CC1 ) + V1*(  R21 );
  FFVar F3 = q0*( CD0 - CD1 ) + V1*( -R11 );
  FFVar F4 = q0*rho*cp*( T0 - T1 ) + V1*( -R11*DH1 - R21*DH2 );
  // reactor 2 (inlet = reactor-1 outlet)
  FFVar k12 = k0_1 * exp( -Ea_1/( T2 + 273.15 ) );
  FFVar k22 = k0_2 * exp( -Ea_2/( T2 + 273.15 ) );
  FFVar R12 = k12 * CA2 * CD2;
  FFVar R22 = k22 * CB2;
  FFVar F5 = q0*( CA1 - CA2 ) + V2*( -R12 );
  FFVar F6 = q0*( CB1 - CB2 ) + V2*(  R12 - R22 );
  FFVar F7 = q0*( CC1 - CC2 ) + V2*(  R22 );
  FFVar F8 = q0*( CD1 - CD2 ) + V2*( -R12 );
  FFVar F9 = q0*rho*cp*( T1 - T2 ) + V2*( -R12*DH1 - R22*DH2 );

  // ---- outputs ----
  FFVar sum2 = CA2 + CB2 + CC2 + CD2;
  FFVar Y0 = CA2 / sum2;
  FFVar Y1 = CC2 / sum2;
  FFVar Y2 = T2;

  OCFESLV oc( &DAG );
  FFVar const st[10] = { CA1,CB1,CC1,CD1,T1, CA2,CB2,CC2,CD2,T2 };
  for( auto const& s : st ) oc.add_state( s, {} );
  FFVar const uu[6] = { CA0, CA0_CD0, q0, V1, V2, T0 };
  for( auto const& u : uu ) oc.add_input( u, {} );

  // ---- nominal controls (test7 dU #2) and state initial guesses ----
  double const dU[6] = { 0.8, 0.8, 0.036, 1.743, 1.826, 22.0 };  // CA0,CA0_CD0,q0,V1,V2,T0
  double const CD0n = dU[0]/dU[1];                                 // feed D at nominal
  // seed states at the feed (the ODE's initial condition, in the physical basin)
  oc.update_ref( CA1, dU[0] ); oc.update_ref( CB1, 0. ); oc.update_ref( CC1, 0. );
  oc.update_ref( CD1, CD0n  ); oc.update_ref( T1,  dU[5] );
  oc.update_ref( CA2, dU[0] ); oc.update_ref( CB2, 0. ); oc.update_ref( CC2, 0. );
  oc.update_ref( CD2, CD0n  ); oc.update_ref( T2,  dU[5] );

  OCFESLV::EqnOptions alg( OCFESLV::EqnRole::INTERIOR, 0 );
  FFVar const res[10] = { F0,F1,F2,F3,F4, F5,F6,F7,F8,F9 };
  for( auto const& f : res ) oc.add_equation( f, {}, {}, alg );

  oc.add_output( Y0, std::vector<FFVar>{}, std::vector<double>{} );   // fct 0
  oc.add_output( Y1, std::vector<FFVar>{}, std::vector<double>{} );   // fct 1
  oc.add_output( Y2, std::vector<FFVar>{}, std::vector<double>{} );   // fct 2

  oc.options.REDUCE.ORDER   = OCFESLV::Options::RED_NONE;
  oc.options.CLASSIFY.MODE       = OCFESLV::Options::CLASS_NONE;
  oc.options.SOLVE.MAX_ITER = 200;
  oc.options.SOLVE.RES_TOL  = 1.0e-9;
  oc.options.DISPLAY_LEVEL  = 0;

  if( !oc.setup() ){
    std::cerr << "ERROR: setup() failed: " << OCFESLV::setup_status_str( oc.setup_status() ) << "\n";
    return 2;
  }

  std::vector<double> xv, inp;
  if( !oc.init( xv, inp, nullptr ) ){ std::cerr << "ERROR: init() failed\n"; return 3; }

  for( size_t i = 0; i < 6; ++i ) oc.set_input_values( uu[i], { dU[i] }, inp.data() );
  for( size_t i = 0; i < 6; ++i ) oc.register_control( uu[i] );

  size_t const ncd = oc.n_control_dof();
  size_t const ncf = oc.n_colloc_fct();
  std::cout << "  states=" << oc.n_colloc_sta() << " eqns=" << oc.n_colloc_eqn()
            << " controls ncd=" << ncd << " outputs ncf=" << ncf << "\n";
  if( ncd != 6 || ncf != 3 ){ std::cerr << "ERROR: expected ncd=6, ncf=3\n"; return 4; }

  std::vector<double> p0; oc.encode_controls( inp.data(), p0 );

  // ---- reference: direct solve at nominal controls ----
  std::vector<double> xvR( xv ), inpR( inp );
  OCFESLV::SolveReport const rep = oc.solve( xvR.data(), inpR.data(), nullptr );
  std::cout << "  reference solve: converged=" << ( rep.converged ? "yes" : "no" )
            << " iters=" << rep.iterations
            << " final|r|=" << std::scientific << std::setprecision(3) << rep.final_residual << "\n\n";
  if( !rep.converged ){ std::cerr << "ERROR: reference solve did not converge\n"; return 5; }
  std::vector<double> Fref = oc.val_functions();
  if( Fref.size() < ncf ){ std::cerr << "ERROR: functions missing\n"; return 6; }

  std::cout << "  [primal] steady-state outputs (informational; test7 dU#2 row ~ 0.0018, 0.020, 67.9):\n";
  std::cout << std::scientific << std::setprecision(6)
            << "    Y0=CA2/S=" << Fref[0] << "  Y1=CC2/S=" << Fref[1] << "  Y2=T2=" << Fref[2] << "\n";

  // ---- residual closure: independently evaluate the ORIGINAL residuals F0..F9 at the
  //      OCFESLV solution and confirm they vanish.  This is the self-contained proof that
  //      the solve found a genuine root of RHS=0 -- independent of WHICH steady state it
  //      is (the exothermic system admits multiple roots; the ODE picks one by dynamics,
  //      Newton another by basin, and both are valid).  Uses a fresh DAG.eval of the
  //      user expressions, not OCFESLV's internal (possibly row-scaled) residual.
  {
    OCFESLV::t_Coord pt;
    std::vector<double> sx( 10 );
    for( size_t i = 0; i < 10; ++i )
      sx[i] = oc.eval_colloc<double>( st[i], pt, xvR.data(), inpR.data(), nullptr );

    // print the 10 solved steady-state values
    char const* xn[10] = { "CA1","CB1","CC1","CD1","T1", "CA2","CB2","CC2","CD2","T2" };
    std::cout << "  [states] steady-state solution:\n"
              << std::scientific << std::setprecision(6);
    for( size_t i = 0; i < 10; ++i ){
      std::cout << "    " << std::left << std::setw(4) << xn[i] << "= "
                << std::right << std::setw(14) << sx[i] << "\n";
      if( i == 4 ) std::cout << "    ---- (reactor 2 below) ----\n";
    }

    // residual closure: F0..F9 at the solution (fresh DAG.eval of the user expressions)
    std::vector<FFVar>  Xvec( st, st+10 ), Uvec( uu, uu+6 );
    std::vector<FFVar>  Fexpr = { F0,F1,F2,F3,F4,F5,F6,F7,F8,F9 };
    std::vector<double> dUvec( dU, dU+6 ), Fout( 10, 0. );
    DAG.eval( Fexpr, Fout, Xvec, sx, Uvec, dUvec );
    double maxres = 0.;
    for( size_t i = 0; i < 10; ++i ) maxres = std::max( maxres, std::fabs( Fout[i] ) );
    std::cout << "    max|F0..F9|(x*) = " << std::scientific << std::setprecision(3) << maxres << "\n";
    check_close( "residuals F0-F9 closed at solution", maxres, 0., 1e-8 );

    // physicality: this system is multi-rooted -- a spurious root with NEGATIVE
    // concentrations also closes the residuals, so residual closure alone does not pin
    // the physical steady state.  Guard that OCFESLV found the physical root by requiring
    // all 8 concentrations non-negative (temperatures T1,T2 are unconstrained here).
    size_t const conc[8] = { 0,1,2,3, 5,6,7,8 };   // CA1,CB1,CC1,CD1, CA2,CB2,CC2,CD2
    double cmin = sx[ conc[0] ];
    for( size_t k = 1; k < 8; ++k ) cmin = std::min( cmin, sx[ conc[k] ] );
    std::cout << "    min concentration = " << std::scientific << std::setprecision(3) << cmin << "\n";
    check_true( "physical root (all concentrations >= 0)", cmin >= -1e-9 );
  }

  // ---- build the reduced-space FFOCFESLV operation ----
  FFGraph rdag;
  std::vector<FFVar> pctrl( ncd );
  for( size_t i = 0; i < ncd; ++i ){ std::ostringstream os; os << "p" << i; pctrl[i] = rdag.add_var( os.str() ); }

  FFOCFESLV ffred;
  std::vector<FFVar> F = ffred( FFOCFESLV::controls_map( oc, pctrl ), &oc, FFOCFESLV::COPY, "AE2" );
  check_true( "FFOCFESLV returns ncf outputs", F.size() == ncf );

  // ===== T1: value via direct eval<double> =====
  {
    std::vector<double> Fval( ncf, 0. );
    ffred.eval( (unsigned)ncf, Fval.data(), (unsigned)ncd, p0.data(), nullptr );
    double e = 0.; for( size_t j = 0; j < ncf; ++j ) e = std::max( e, std::fabs( Fval[j] - Fref[j] ) );
    check_close( "T1 value (direct eval) == solve", e, 0., 1e-9 );
  }

  // ===== T3: forward AD Jacobian (18 = 3x6) =====
  std::vector<FADType<double>> xF( ncd ), yF( ncf );
  for( size_t i = 0; i < ncd; ++i ){ xF[i] = p0[i]; xF[i].diff( (unsigned)i, (unsigned)ncd ); }
  ffred.eval( (unsigned)ncf, yF.data(), (unsigned)ncd, xF.data(), nullptr );

  std::cout << "\n  [T3v] forward-sweep primal values:\n";
  { double e = 0.; for( size_t j = 0; j < ncf; ++j ) e = std::max( e, std::fabs( yF[j].val() - Fref[j] ) );
    check_close( "T3v F-sweep values == solve", e, 0., 1e-9 ); }

  std::cout << "\n  [T3a] forward AD  vs  oc.solve_fsens / sens_jacobian:\n";
  {
    std::vector<double> xvS( xv ), inpS( inp );
    if( oc.solve_fsens( xvS.data(), inpS.empty()?nullptr:inpS.data(), nullptr ) ){
      std::vector<double> const& J = oc.sens_jacobian();
      if( J.size() >= ncf*ncd ){
        double e = 0.;
        for( size_t j = 0; j < ncf; ++j )
          for( size_t i = 0; i < ncd; ++i )
            e = std::max( e, std::fabs( yF[j].deriv((unsigned)i) - J[ j*ncd + i ] ) );
        check_close( "T3a FAD Jacobian == reduced fwd-sens", e, 0., 1e-7 );
      }
      else check_true( "T3a sens_jacobian sized ncf*ncd", false );
    }
    else check_true( "T3a solve_fsens available", false );
  }

  std::cout << "\n  [T3b] forward AD  vs  central differences (per control):\n";
  {
    for( size_t i = 0; i < ncd; ++i ){
      double const eps = 1e-6*std::max( 1.0, std::fabs( p0[i] ) );
      std::vector<double> pp( p0 ), pm( p0 ); pp[i] += eps; pm[i] -= eps;
      std::vector<double> Fp( ncf, 0. ), Fm( ncf, 0. );
      ffred.eval( (unsigned)ncf, Fp.data(), (unsigned)ncd, pp.data(), nullptr );
      ffred.eval( (unsigned)ncf, Fm.data(), (unsigned)ncd, pm.data(), nullptr );
      double emax = 0.;
      for( size_t j = 0; j < ncf; ++j ){
        double const fd = ( Fp[j] - Fm[j] )/( 2.0*eps );
        double const an = yF[j].deriv( (unsigned)i );
        double const rel = std::fabs( an - fd )/( std::max( 1e-8, std::fabs( an ) ) );
        emax = std::max( emax, rel );
      }
      std::ostringstream nm; nm << "dY/dp" << i << " (FAD) ~ FD (max rel)";
      check_close( nm.str().c_str(), emax, 0., 5e-4 );
    }
  }

  // ===== T4: symbolic Jacobian via SFAD =====
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
    for( size_t k = 0; k < sd.size(); ++k ){
      if( si[k] >= ncf || sj[k] >= ncd ){ idx_ok = false; continue; }
      e = std::max( e, std::fabs( sdv[k] - yF[si[k]].deriv( (unsigned)sj[k] ) ) );
      ++cnt;
    }
    check_true ( "T4 SFAD indices in range", idx_ok && cnt == ncf*ncd );
    check_close( "T4 SFAD Jacobian == forward AD", e, 0., 1e-9 );
  }

  std::cout << "\n================================================================\n"
            << "  OCFE_AE2: " << g_pass << " passed, " << g_fail << " failed -- "
            << ( g_fail == 0 ? "ALL PASS" : "SOME FAILED" ) << "\n"
            << "================================================================\n";
  return g_fail == 0 ? 0 : 1;
}
