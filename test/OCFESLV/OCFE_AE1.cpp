// ============================================================================
//  OCFE_AE1.cpp  --  OCFESLV on a PARAMETRIC algebraic system + FFOCFESLV
//                    reduced-space outputs, with derivative-machinery tests.
//
//  Extends OCFE_AE.cpp (completely algebraic, zero-domain) with INPUTS and
//  OUTPUTS, wraps the OCFESLV solve in the reduced-space external DAG operation
//  FFOCFESLV (ffocfe.hpp), and validates the input->output derivative machinery
//  against a CLOSED-FORM oracle -- so this is a genuine sensitivity test, not just
//  an internal-consistency check.
//
//  MODEL (2 states, 2 scalar inputs/parameters, no domain)
//  -------------------------------------------------------
//    states  x1, x2      inputs (controls)  a, b
//      F1 :  x1 - a          = 0     ->  x1 = a
//      F2 :  x1*x2 - b       = 0     ->  x2 = b/a          (nonlinear, couples x2,b,a)
//    outputs
//      G0 :  x1 + x2  = a + b/a
//      G1 :  x1*x2    = b
//
//  CLOSED-FORM reduced Jacobian  dF/d(a,b)  (F=[G0,G1], p=[a,b]):
//      dG0/da = 1 - b/a^2      dG0/db = 1/a
//      dG1/da = 0              dG1/db = 1
//  At the nominal (a,b)=(2,6):  x=(2,3),  F=(5,6),  J = [[-0.5, 0.5],[0, 1]].
//
//  WHAT THIS EXERCISES
//  -------------------
//    * add_input(p,{}) / set_input_values / register_control / encode_controls
//      on a zero-domain system;
//    * add_output(expr,{},{}) scalar outputs and val_functions();
//    * FFOCFESLV( controls_map( oc, pctrl ), &oc, COPY ) reduced-space wrapping (one-map form);
//    * three independent derivative paths, each vs the CLOSED-FORM oracle:
//        T3a  forward AD           FFOCFESLV::eval<FADType<double>>
//        T3b  reduced fwd-sens     oc.solve_fsens + oc.sens_jacobian
//        T3c  central differences  through FFOCFESLV::eval<double>
//        T4   symbolic Jacobian    rdag.SFAD -> FFGradOCFESLV, evaluated
//
//  Build:
//    g++ -std=c++17 <suite flags> -DOCFE_OCFESLV_HEADER='"ocfeslv.hpp"' \
//        OCFE_AE1.cpp -o OCFE_AE1 <libs>
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
  std::cout << "  " << std::left << std::setw(44) << name
            << std::right << std::scientific << std::setprecision(3)
            << " got=" << std::setw(11) << got << " want=" << std::setw(11) << want
            << " |d|=" << std::setw(9) << std::fabs( got - want )
            << "  " << ( ok ? "PASS" : "FAIL" ) << "\n";
}
static void check_true( char const* name, bool ok )
{
  ( ok ? g_pass : g_fail )++;
  std::cout << "  " << std::left << std::setw(44) << name << "  "
            << ( ok ? "PASS" : "FAIL" ) << "\n";
}

int main()
{
  std::cout << "================================================================\n"
            << "  OCFE_AE1 : parametric algebraic system + FFOCFESLV outputs\n"
            << "  zero-domain; input->output derivative machinery vs closed form\n"
            << "================================================================\n";

  FFGraph DAG;
  FFVar x1 = DAG.add_var( "x1" );
  FFVar x2 = DAG.add_var( "x2" );
  FFVar a  = DAG.add_var( "a" );      // input / control
  FFVar b  = DAG.add_var( "b" );      // input / control

  OCFESLV oc( &DAG );
  // No add_domain(): zero collocation domains.
  oc.add_state( x1, {} );
  oc.add_state( x2, {} );
  oc.add_input( a, {} );
  oc.add_input( b, {} );

  // Seed the states off the origin (F2's Jacobian [[.,.],[x2,x1]] is singular there).
  oc.update_ref( x1, 1.0 );
  oc.update_ref( x2, 1.0 );

  FFVar F1 = x1 - a;
  FFVar F2 = x1 * x2 - b;
  FFVar G0 = x1 + x2;      // output fct 0
  FFVar G1 = x1 * x2;      // output fct 1

  OCFESLV::EqnOptions alg( OCFESLV::EqnRole::INTERIOR, 0 );
  oc.add_equation( F1, {}, {}, alg );
  oc.add_equation( F2, {}, {}, alg );
  oc.add_output( G0, std::vector<FFVar>{}, std::vector<double>{} );   // fct 0 (empty obs -> double overload)
  oc.add_output( G1, std::vector<FFVar>{}, std::vector<double>{} );   // fct 1

//  oc.options.REDUCE.ORDER   = OCFESLV::Options::RED_NONE;
//  oc.options.CLASSIFY.MODE       = OCFESLV::Options::CLASS_NONE;
  oc.options.SOLVE.MAX_ITER = 60;
  oc.options.SOLVE.RES_TOL  = 1.0e-9;
  oc.options.DISPLAY_LEVEL  = 1;

  if( !oc.setup() ){
    std::cerr << "ERROR: setup() failed: "
              << OCFESLV::setup_status_str( oc.setup_status() ) << "\n";
    return 2;
  }

  std::vector<double> xv, inp;
  if( !oc.init( xv, inp, nullptr ) ){ std::cerr << "ERROR: init() failed\n"; return 3; }

  double const a_nom = 2.0, b_nom = 6.0;
  oc.set_input_values( a, { a_nom }, inp.data() );
  oc.set_input_values( b, { b_nom }, inp.data() );

  oc.register_control( a );
  oc.register_control( b );
  size_t const ncd = oc.n_control_dof();
  size_t const ncf = oc.n_colloc_fct();
  auto const& C  = oc.controls();
  size_t const ia = C.at( a ).offset;
  size_t const ib = C.at( b ).offset;
  std::cout << "  states=" << oc.n_colloc_sta() << " eqns=" << oc.n_colloc_eqn()
            << " inputs->controls ncd=" << ncd << " (a@" << ia << ", b@" << ib << ")"
            << " outputs ncf=" << ncf << "\n\n";
  if( ncd != 2 || ncf != 2 ){ std::cerr << "ERROR: expected ncd=2, ncf=2\n"; return 4; }

  std::vector<double> p0; oc.encode_controls( inp.data(), p0 );

  // ---- reference: direct solve at nominal controls ----
  std::vector<double> xvR( xv ), inpR( inp );
  OCFESLV::SolveReport const rep = oc.solve( xvR.data(), inpR.data(), nullptr );
  std::cout << "  reference solve: converged=" << ( rep.converged ? "yes" : "no" )
            << " iters=" << rep.iterations
            << " final|r|=" << std::scientific << std::setprecision(3)
            << rep.final_residual << "\n";
  if( !rep.converged ){ std::cerr << "ERROR: reference solve did not converge\n"; return 5; }
  std::vector<double> Fref = oc.val_functions();
  if( Fref.size() < ncf ){ std::cerr << "ERROR: functions missing\n"; return 6; }

  OCFESLV::t_Coord pt;
  double const x1v = oc.eval_colloc<double>( x1, pt, xvR.data(), inpR.data(), nullptr );
  double const x2v = oc.eval_colloc<double>( x2, pt, xvR.data(), inpR.data(), nullptr );

  std::cout << "\n  [primal] states and outputs:\n";
  check_close( "x1 == a (2)",        x1v,     2.0, 1e-9 );
  check_close( "x2 == b/a (3)",      x2v,     3.0, 1e-9 );
  check_close( "G0 == x1+x2 (5)",    Fref[0], 5.0, 1e-9 );
  check_close( "G1 == x1*x2 (6)",    Fref[1], 6.0, 1e-9 );

  // ---- closed-form reduced Jacobian oracle  Jan[f][control] ----
  double const av = p0[ia], bv = p0[ib];
  double Jan[2][2];
  Jan[0][ia] = 1.0 - bv/(av*av);  Jan[0][ib] = 1.0/av;   // dG0/da, dG0/db
  Jan[1][ia] = 0.0;               Jan[1][ib] = 1.0;      // dG1/da, dG1/db

  // ---- build the reduced-space FFOCFESLV operation ----
  FFGraph rdag;
  std::vector<FFVar> pctrl( ncd );
  for( size_t i = 0; i < ncd; ++i ){ std::ostringstream os; os << "p" << i; pctrl[i] = rdag.add_var( os.str() ); }

  FFOCFESLV ffred;
  std::vector<FFVar> F = ffred( FFOCFESLV::controls_map( oc, pctrl ), &oc, FFOCFESLV::COPY, "AE1" );
  check_true( "FFOCFESLV returns ncf outputs", F.size() == ncf );

  // ===== T1: value via direct eval<double> =====
  {
    std::vector<double> Fval( ncf, 0. );
    ffred.eval( (unsigned)ncf, Fval.data(), (unsigned)ncd, p0.data(), nullptr );
    double e = 0.; for( size_t j = 0; j < ncf; ++j ) e = std::max( e, std::fabs( Fval[j] - Fref[j] ) );
    check_close( "T1 value (direct eval) == solve", e, 0., 1e-9 );
  }
  // ===== T2: value via DAG round-trip =====
  {
    std::vector<double> Fdag( ncf, 0. );
    rdag.eval( F, Fdag, pctrl, p0 );
    double e = 0.; for( size_t j = 0; j < ncf; ++j ) e = std::max( e, std::fabs( Fdag[j] - Fref[j] ) );
    check_close( "T2 value (DAG round-trip) == solve", e, 0., 1e-9 );
  }

  // ===== T3: forward AD Jacobian, checked three ways =====
  std::vector<FADType<double>> xF( ncd ), yF( ncf );
  for( size_t i = 0; i < ncd; ++i ){ xF[i] = p0[i]; xF[i].diff( (unsigned)i, (unsigned)ncd ); }
  ffred.eval( (unsigned)ncf, yF.data(), (unsigned)ncd, xF.data(), nullptr );

  std::cout << "\n  [T3v] forward-sweep primal values:\n";
  { double e = 0.; for( size_t j = 0; j < ncf; ++j ) e = std::max( e, std::fabs( yF[j].val() - Fref[j] ) );
    check_close( "T3v F-sweep values == solve", e, 0., 1e-9 ); }

  std::cout << "\n  [T3a] forward AD  vs  CLOSED-FORM oracle:\n";
  for( size_t j = 0; j < ncf; ++j )
    for( size_t i = 0; i < ncd; ++i ){
      std::ostringstream nm; nm << "dG" << j << "/dp" << i << " (FAD) == analytic";
      check_close( nm.str().c_str(), yF[j].deriv( (unsigned)i ), Jan[j][i], 1e-8 );
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
            e = std::max( e, std::fabs( yF[j].deriv( (unsigned)i ) - J[ j*ncd + i ] ) );
        check_close( "T3b FAD Jacobian == reduced fwd-sens", e, 0., 1e-9 );
      }
      else check_true( "T3b sens_jacobian sized ncf*ncd", false );
    }
    else check_true( "T3b solve_fsens available", false );
  }

  std::cout << "\n  [T3c] forward AD  vs  central differences:\n";
  {
    double const eps = 1e-5;
    for( size_t i = 0; i < ncd; ++i ){
      std::vector<double> pp( p0 ), pm( p0 ); pp[i] += eps; pm[i] -= eps;
      std::vector<double> Fp( ncf, 0. ), Fm( ncf, 0. );
      ffred.eval( (unsigned)ncf, Fp.data(), (unsigned)ncd, pp.data(), nullptr );
      ffred.eval( (unsigned)ncf, Fm.data(), (unsigned)ncd, pm.data(), nullptr );
      for( size_t j = 0; j < ncf; ++j ){
        double const fd  = ( Fp[j] - Fm[j] )/( 2.0*eps );
        double const an  = yF[j].deriv( (unsigned)i );
        double const tol = 1e-4*std::max( 1e-6, std::fabs( an ) ) + 1e-6;
        std::ostringstream nm; nm << "dG" << j << "/dp" << i << " (FAD) ~ FD";
        check_close( nm.str().c_str(), an, fd, tol );
      }
    }
  }

  // ===== T4: symbolic Jacobian via SFAD (FFGradOCFESLV) =====
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

  // ===== T5: ADJOINT sensitivity (reverse sweep) vs forward AD + oracle =====
  // AE1 is a pure algebraic (steady, no evolution) system, so solve_asens() routes through the
  // MONOLITHIC adjoint apply -- this validates the mono adjoint attribution for input->control DOFs.
  std::cout << "\n  [T5] ADJOINT sensitivity  vs  forward AD + CLOSED-FORM oracle:\n";
  {
    std::vector<double> xvA( xv ), inpA( inp );
    if( oc.solve_asens( xvA.data(), inpA.empty()?nullptr:inpA.data(), nullptr ) ){
      std::vector<double> const& Ja = oc.sens_jacobian();
      if( Ja.size() >= ncf*ncd ){
        double ef = 0., eo = 0.;
        for( size_t j = 0; j < ncf; ++j )
          for( size_t i = 0; i < ncd; ++i ){
            ef = std::max( ef, std::fabs( Ja[ j*ncd + i ] - yF[j].deriv( (unsigned)i ) ) );
            eo = std::max( eo, std::fabs( Ja[ j*ncd + i ] - Jan[j][i] ) );
          }
        check_close( "T5 adjoint Jacobian == forward AD",     ef, 0., 1e-9 );
        check_close( "T5 adjoint Jacobian == analytic oracle", eo, 0., 1e-8 );
      }
      else check_true( "T5 sens_jacobian sized ncf*ncd", false );
    }
    else check_true( "T5 solve_asens available", false );
  }

  std::cout << "\n================================================================\n"
            << "  OCFE_AE1: " << g_pass << " passed, " << g_fail << " failed -- "
            << ( g_fail == 0 ? "ALL PASS" : "SOME FAILED" ) << "\n"
            << "================================================================\n";
  return g_fail == 0 ? 0 : 1;
}
