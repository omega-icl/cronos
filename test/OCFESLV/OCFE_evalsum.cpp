// ============================================================================
//  OCFE_EVALSUM.cpp  --  outputs that COMBINE several evaluations and integrals,
//                        on OCFESLV, monolithic and marching.
//
//  MODEL (Lotka-Volterra, 2 states, 1 input p = control)
//        dx0/dt = p x0 (1 - x1),  dx1/dt = p x1 (x0 - 1),   x0(0) = 1.2,  x1(0) = 1.1 + 0.01 p
//        t in [0, tf], N_EL elements, LGR collocation; T1 = tf/2 is an element boundary
//
//  REFERENCE: a model with FOUR SEPARATE outputs (each one reduction)
//        e0 = x0^2(T1),  e1 = x0^2(tf),  e2 = x1(tf),  q = int_0^tf x1 dt
//  Each COMBINED output is declared alone in its own model and must reproduce, ON THE SAME MESH, the same
//  combination of the references -- value, forward reduced Jacobian (solve_fsens) and adjoint (solve_asens):
//        C1  e0 + e1                   evaluations at two different times
//        C2  e1 + e1                   the same evaluation twice (one record)
//        C3  2 e0 - 3 q + 5            evaluation, integral and a constant
//        C4  p q + e2                  a PARAMETER coefficient on an integral
//        C5  e0 * e2                   NONLINEAR (product rule for the derivative)
//  Tolerances are at solver level (same discretisation on both sides), not discretisation error.
//  ORACLES independent of OCFESLV's self-consistency:
//        O1  the 4 reference values and dF/dp == CVODES (ODESLVS on the same model, RTOL = ATOL = 1e-12), at
//            discretisation tolerance
//        O2  marching == monolithic on the 4 references (value and dF/dp), at discretisation tolerance
//  Every case in both modes: monolithic, then marching (a value captured at T1 in one window must reach an output
//  assembled in a later window).
//
//  Build:
//    g++ -std=c++17 <suite flags> -DOCFE_OCFESLV_HEADER='"ocfeslv.hpp"' OCFE_EVALSUM.cpp -o OCFE_EVALSUM <libs>
// ============================================================================

#include <cmath>
#include <functional>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#ifndef OCFE_OCFESLV_HEADER
#define OCFE_OCFESLV_HEADER "ocfeslv.hpp"
#endif
#include OCFE_OCFESLV_HEADER

using namespace mc;

static double const tf = 6.0, T1 = 3.0, P = 3.0;
static size_t const N_EL = 24, N_ND = 8;          // fine enough for the CVODES oracle below

// CVODES reference (ODESLVS, same model and p, RTOL = ATOL = 1e-12; computed 2026-09-28): e0, e1, e2, q and dF/dp
static double const CV_F[4] = { 6.319583e-01, 1.490090e+00, 9.112639e-01, 5.994301e+00 };
static double const CV_D[4] = { 1.650006e-01, 1.606430e+00, 1.209847e+00, -1.777797e-01 };
static std::vector<double> g_monoF, g_monoD;      // monolithic references, for O2

static int g_pass = 0, g_fail = 0;
static void check_close( std::string const& name, double got, double want, double tol )
{
  double const d = std::fabs( got - want ) / std::max( 1., std::fabs( want ) );
  bool const ok = d <= tol;  ( ok ? g_pass : g_fail )++;
  std::cout << "  " << std::left << std::setw(56) << name << std::right << std::scientific << std::setprecision(3)
            << " rel=" << d << " tol=" << tol << "  " << ( ok ? "PASS" : "FAIL" ) << "\n";
}
static void check_true( std::string const& name, bool ok )
{
  ( ok ? g_pass : g_fail )++;
  std::cout << "  " << std::left << std::setw(56) << name << "  " << ( ok ? "PASS" : "FAIL" ) << "\n";
}

struct M { FFGraph DAG; OCFESLV oc{ &DAG }; FFVar t, x0, x1, p; };
//! @brief Build the LV model; @p outputs declares the outputs.  Returns setup().
static bool build( M& m, bool marching, std::function<void( M& )> outputs )
{
  m.t = m.DAG.add_var( "t" );  m.x0 = m.DAG.add_var( "x0(t)" );  m.x1 = m.DAG.add_var( "x1(t)" );  m.p = m.DAG.add_var( "p" );
  OCFESLV& oc = m.oc;
  oc.add_domain( m.t, FFDom( 0., tf, N_EL, FFDom::LGR, N_ND ) );  oc.set_evolution_domain( m.t );
  oc.add_state( m.x0, {m.t} );  oc.add_state( m.x1, {m.t} );  oc.add_input( m.p, {} );
  oc.update_ref( m.x0, 1.2 );  oc.update_ref( m.x1, 1.1 );
  FFPartial OpP;
  int const T_INT = FFDom::ALL - FFDom::LB;
  oc.add_equation( OpP( m.x0, m.t ) - m.p*m.x0*( 1. - m.x1 ), {m.t}, {T_INT},       OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( OpP( m.x1, m.t ) - m.p*m.x1*( m.x0 - 1. ), {m.t}, {T_INT},       OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( m.x0 - 1.2,                                {m.t}, {FFDom::LB},   OCFESLV::EqnOptions( OCFESLV::EqnRole::INITIAL,  0 ) );
  oc.add_equation( m.x1 - ( 1.1 + 0.01*m.p ),                 {m.t}, {FFDom::LB},   OCFESLV::EqnOptions( OCFESLV::EqnRole::INITIAL,  0 ) );
  outputs( m );
  oc.options.SOLVE.MARCHING = marching;  oc.options.SOLVE.MAX_ITER = 60;  oc.options.SOLVE.RES_TOL = 1e-12;
  oc.options.DISPLAY_LEVEL = 0;
  return oc.setup();
}

struct R { bool ok = false; std::vector<double> F, Jf, Ja; };      // J: nf x 1 (the control p)
static R solve( M& m )
{
  R r;  OCFESLV& oc = m.oc;
  std::vector<double> xv, inp;
  if( !oc.init( xv, inp, nullptr ) ) return r;
  oc.set_input_values( m.p, { P }, inp.data() );
  oc.register_control( m.p );
  std::vector<double> x1( xv ), x2( xv ), x3( xv );
  auto const rep = oc.solve( x1.data(), inp.data(), nullptr );
  r.F = oc.val_functions();
  bool const fs = oc.solve_fsens( x2.data(), inp.data(), nullptr );  r.Jf = oc.sens_jacobian();
  bool const as = oc.solve_asens( x3.data(), inp.data(), nullptr );  r.Ja = oc.sens_jacobian();
  r.ok = rep.converged && fs && as && r.F.size() >= 1 && r.Jf.size() == r.F.size() && r.Ja.size() == r.F.size();
  return r;
}

static void run_mode( bool marching )
{
  char const* mstr = marching ? "MARCHING" : "MONOLITHIC";
  std::cout << "\n================ " << mstr << " ================\n";
  FFEval OpE;  FFIntegral OpI;
  auto Q = [&]( M& m ){ return OpI( std::vector<FFVar>{ m.x1 }, { m.t } )[0]; };

  M ref;
  bool const rs = build( ref, marching, [&]( M& m ){
    m.oc.add_output( OpE( sqr( m.x0 ), m.t, T1 ) );  m.oc.add_output( OpE( sqr( m.x0 ), m.t, tf ) );
    m.oc.add_output( OpE( m.x1, m.t, tf ) );          m.oc.add_output( Q( m ) ); } );
  check_true( std::string( "REF setup, 4 separate outputs [" ) + mstr + "]", rs && ref.oc.n_colloc_fct() == 4 );
  if( !rs ) return;
  if( ref.oc.is_marching() != marching ){ check_true( "REF is_marching()==requested", false ); return; }
  R const r0 = solve( ref );
  check_true( "REF solve, fsens and asens", r0.ok );
  if( !r0.ok ) return;
  { double w = 0.; for( size_t k = 0; k < 4; ++k ) w = std::max( w, std::fabs( r0.Ja[k] - r0.Jf[k] ) / std::max( 1., std::fabs( r0.Jf[k] ) ) );
    check_close( "REF adjoint == forward on the 4 references", w, 0., 1e-8 ); }
  { double wv = 0., wd = 0.;
    for( size_t k = 0; k < 4; ++k ){ wv = std::max( wv, std::fabs( r0.F[k]  - CV_F[k] ) / std::max( 1., std::fabs( CV_F[k] ) ) );
                                     wd = std::max( wd, std::fabs( r0.Jf[k] - CV_D[k] ) / std::max( 1., std::fabs( CV_D[k] ) ) ); }
    check_close( "O1 references: value == CVODES (worst)", wv, 0., 1e-5 );
    check_close( "O1 references: dF/dp == CVODES (worst)", wd, 0., 1e-4 ); }
  if( !marching ){ g_monoF = r0.F; g_monoD = r0.Jf; }
  else if( g_monoF.size() == 4 ){
    double wv = 0., wd = 0.;
    for( size_t k = 0; k < 4; ++k ){ wv = std::max( wv, std::fabs( r0.F[k]  - g_monoF[k] ) / std::max( 1., std::fabs( g_monoF[k] ) ) );
                                     wd = std::max( wd, std::fabs( r0.Jf[k] - g_monoD[k] ) / std::max( 1., std::fabs( g_monoD[k] ) ) ); }
    check_close( "O2 references: marching value == monolithic (worst)", wv, 0., 1e-5 );
    check_close( "O2 references: marching dF/dp == monolithic (worst)", wd, 0., 1e-4 );
  }
  double const e0 = r0.F[0], e1 = r0.F[1], e2 = r0.F[2], q = r0.F[3];
  double const d0 = r0.Jf[0], d1 = r0.Jf[1], d2 = r0.Jf[2], dq = r0.Jf[3];

  struct C { std::string name; std::function<FFVar( M& )> g; double val, der; };
  std::vector<C> cases{
    { "C1 e0 + e1",          [&]( M& m ){ return OpE( sqr( m.x0 ), m.t, T1 ) + OpE( sqr( m.x0 ), m.t, tf ); },           e0 + e1,          d0 + d1 },
    { "C2 e1 + e1",          [&]( M& m ){ return OpE( sqr( m.x0 ), m.t, tf ) + OpE( sqr( m.x0 ), m.t, tf ); },           2.*e1,            2.*d1 },
    { "C3 2 e0 - 3 q + 5",   [&]( M& m ){ return 2.*OpE( sqr( m.x0 ), m.t, T1 ) - 3.*Q( m ) + 5.; },                   2.*e0 - 3.*q + 5., 2.*d0 - 3.*dq },
    { "C4 p q + e2",         [&]( M& m ){ return m.p*Q( m ) + OpE( m.x1, m.t, tf ); },                                  P*q + e2,         q + P*dq + d2 },
    { "C5 e0 * e2",          [&]( M& m ){ return OpE( sqr( m.x0 ), m.t, T1 ) * OpE( m.x1, m.t, tf ); },                 e0*e2,            d0*e2 + e0*d2 } };
  for( auto const& c : cases ){
    M m;  bool const su = build( m, marching, [&]( M& mm ){ mm.oc.add_output( c.g( mm ) ); } );
    check_true( c.name + ": setup, ONE output", su && m.oc.n_colloc_fct() == 1 );
    if( !su ) continue;
    R const r = solve( m );
    check_true( c.name + ": solve, fsens and asens", r.ok );
    if( !r.ok ) continue;
    check_close( c.name + ": value   == reference combination", r.F[0],  c.val, 1e-9 );
    check_close( c.name + ": forward == reference combination", r.Jf[0], c.der, 1e-8 );
    check_close( c.name + ": adjoint == forward",               r.Ja[0], r.Jf[0], 1e-8 );
  }
}

int main()
{
  std::cout << "================================================================\n"
            << "  OCFE_EVALSUM : outputs combining evaluations and integrals (OCFESLV)\n"
            << "  Lotka-Volterra on [0," << tf << "], T1 = " << T1 << ", " << N_EL << " LGR elements x " << N_ND << " nodes\n"
            << "================================================================\n";
  run_mode( false );
  run_mode( true );
  std::cout << "\n============================================================\n"
            << "  RESULT: " << g_pass << " passed, " << g_fail << " failed  -- " << ( g_fail == 0 ? "ALL PASS" : "SOME FAILED" ) << "\n"
            << "============================================================\n";
  return g_fail == 0 ? 0 : 1;
}
