// TRANS_ode.cpp -- ODESLV transitions (add_transition), against a CLOSED FORM: an impulsive dose.
//   dx/dt = -k x, x(0) = x0;  dy/dt = x, y(0) = 0;   at tau = 0.4:  x(tau+) = x(tau-) + d   (y NOT mentioned: continuous)
//   x(tau-) = x0 e^{-k tau};  x(tau+) = x(tau-) + d;  x(tf) = x(tau+) e^{-k (tf - tau)};  y(tf) = int_0^tf x dt
// Checks: values and FORWARD gradients w.r.t. (k, d, x0) against the closed form; the IDENTITY transition reproduces
// the model without one; ADJOINT gradients (incl. a registered subset of controls) == closed form; implicit maps REFUSED.
#include <cmath>
#include <cstdio>
#include <sstream>
#include <iostream>
#include <string>
#include <vector>
#include "ffode.hpp"
using namespace mc;
static int npass = 0, nfail = 0;
static void check( bool c, std::string const& w, double v = -1. ){ std::printf( "  %s  %-66s", c? "PASS": "FAIL", w.c_str() ); if( v >= 0. ) std::printf( " (%.2e)", v ); std::printf( "\n" ); c? ++npass: ++nfail; }
double const TAU = 0.4, TF = 1.0;
static std::vector<double> closed( double k, double d, double x0 ){
  double const xm = x0*std::exp( -k*TAU ), xp = xm + d, xf = xp*std::exp( -k*( TF-TAU ) );
  double const I = x0*( 1.-std::exp( -k*TAU ) )/k + xp*( 1.-std::exp( -k*( TF-TAU ) ) )/k;
  return { xm, xp, xf, I, I }; }                                            // x(tau-), x(tau+), x(tf), y(tf), int x
struct Model { FFGraph G; ODESLVS_CVODES I{ &G }; FFVar t, x, y, k, d, x0; };
static bool build( Model& m, int kind, bool only_d = false ){   // kind 0: dose; 1: identity; 2: no transition; 3: implicit; only_d: register d alone
  m.t = m.G.add_var( "t" ); m.x = m.G.add_var( "x(t)" ); m.y = m.G.add_var( "y(t)" ); m.k = m.G.add_var( "k" ); m.d = m.G.add_var( "d" ); m.x0 = m.G.add_var( "x0" );
  FFPartial OpP; FFEval OpE; FFIntegral OpI;
  m.I.add_domain( m.t, FFDom( std::vector<double>{ 0., 0.5, 1. }, FFDom::LGR, 4 ) ); m.I.set_evolution_domain( m.t );
  m.I.add_state( m.x, {m.t} ); m.I.add_state( m.y, {m.t} ); m.I.add_input( m.k ); m.I.add_input( m.d ); m.I.add_input( m.x0 );
  m.I.update_ref( m.x, 1. ); m.I.update_ref( m.y, 0. );
  int const T_INT = FFDom::ALL - FFDom::LB;
  m.I.add_equation( OpP( m.x, m.t ) + m.k*m.x, {m.t}, {T_INT}, FFModel::EqnOptions( FFModel::EqnRole::INTERIOR, 0 ) );
  m.I.add_equation( OpP( m.y, m.t ) - m.x,     {m.t}, {T_INT}, FFModel::EqnOptions( FFModel::EqnRole::INTERIOR, 0 ) );
  m.I.add_equation( m.x - m.x0, {m.t}, {FFDom::LB}, FFModel::EqnOptions( FFModel::EqnRole::INITIAL, 0 ) );
  m.I.add_equation( m.y,        {m.t}, {FFDom::LB}, FFModel::EqnOptions( FFModel::EqnRole::INITIAL, 0 ) );
  if( kind == 0 ) m.I.add_transition( m.x + m.d, m.x, m.t, TAU );
  if( kind == 1 ) m.I.add_transition( m.x, m.x, m.t, TAU );
  if( kind == 3 ) m.I.add_transition( m.x + m.d, 2.*m.x, m.t, TAU );
  m.I.add_output( OpE( m.x, m.t, TAU, FFDom::MINUS ) ); m.I.add_output( OpE( m.x, m.t, TAU, FFDom::PLUS ) );
  m.I.add_output( OpE( m.x, m.t, TF ) ); m.I.add_output( OpE( m.y, m.t, TF ) ); m.I.add_output( OpI( std::vector<FFVar>{ m.x }, { m.t } )[0] );
  m.I.options.DISPLAY = 0; m.I.options.ATOL = m.I.options.ATOLS = m.I.options.ATOLB = 1e-12; m.I.options.RTOL = m.I.options.RTOLS = m.I.options.RTOLB = 1e-12;
  if( only_d ) m.I.register_control( m.d );                                 // registered BEFORE setup
  std::ostringstream os; auto* old = std::cerr.rdbuf( os.rdbuf() ); bool const ok = m.I.setup(); std::cerr.rdbuf( old );
  return ok; }
static std::vector<double> P_of( Model& m, double k, double d, double x0 ){
  std::vector<double> P( m.I.np() ); P[ m.I.parameter_index( m.k )[0] ] = k; P[ m.I.parameter_index( m.d )[0] ] = d; P[ m.I.parameter_index( m.x0 )[0] ] = x0; return P; }
int main(){
  double const K = 1.3, D = 0.7, X0 = 1.1;
  char const* nm[5] = { "x(tau-)", "x(tau+)", "x(tf)", "y(tf) (y continuous)", "int x dt" };
  { Model m; bool const su = build( m, 0 ); check( su, "dose: setup (transition accepted by ODESLV)" );
    if( su ){
      auto const P = P_of( m, K, D, X0 );
      bool const ok = m.I.solve_fsens( P ) == ODESLVS_CVODES::STATUS::NORMAL; auto const F = m.I.val_function(); auto const G = m.I.val_function_gradient();
      auto const F0 = closed( K, D, X0 );
      for( size_t j = 0; j < 5; ++j ) check( ok && std::fabs( F[j]-F0[j] ) < 1e-9, std::string( "dose: " ) + nm[j] + " == closed form", ok? std::fabs( F[j]-F0[j] ): 1. );
      double wg = 0.; FFVar const* pv[3] = { &m.k, &m.d, &m.x0 }; double const pr[3] = { K, D, X0 };
      for( int q = 0; q < 3; ++q ){ double h = 1e-6*std::max( 1., pr[q] ); double a[3] = { K, D, X0 }, b[3] = { K, D, X0 }; a[q] += h; b[q] -= h;
        auto const fp = closed( a[0], a[1], a[2] ), fm = closed( b[0], b[1], b[2] ); size_t const ip = m.I.parameter_index( *pv[q] )[0];
        for( size_t j = 0; j < 5 && ok; ++j ) wg = std::max( wg, std::fabs( G[ip][j] - ( fp[j]-fm[j] )/( 2.*h ) ) ); }
      check( ok && wg < 1e-7, "dose: FORWARD gradients w.r.t. (k, d, x0) == closed form", wg );
      bool const oa = m.I.solve_asens( P ) == ODESLVS_CVODES::STATUS::NORMAL; auto const Ga = m.I.val_function_gradient();
      double wa = 0.;
      for( int q = 0; q < 3; ++q ){ double h = 1e-6*std::max( 1., pr[q] ); double a[3] = { K, D, X0 }, b[3] = { K, D, X0 }; a[q] += h; b[q] -= h;
        auto const fp = closed( a[0], a[1], a[2] ), fm = closed( b[0], b[1], b[2] ); size_t const ip = m.I.parameter_index( *pv[q] )[0];
        for( size_t j = 0; j < 5 && oa; ++j ) wa = std::max( wa, std::fabs( Ga[ip][j] - ( fp[j]-fm[j] )/( 2.*h ) ) ); }
      check( oa && wa < 1e-7, "dose: ADJOINT gradients w.r.t. (k, d, x0) == closed form", oa? wa: 1. ); } }
  { Model m;                                                                  // SELECTIVE: d only -- the jump increment must land on it
    bool const su = build( m, 0, true ); auto const P = P_of( m, K, D, X0 );
    bool const o1 = su && m.I.solve_fsens( P ) == ODESLVS_CVODES::STATUS::NORMAL; auto const Gf = m.I.val_function_gradient();
    bool const o2 = su && m.I.solve_asens( P ) == ODESLVS_CVODES::STATUS::NORMAL;     auto const Ga = m.I.val_function_gradient();
    double const h = 1e-6; auto const fp = closed( K, D+h, X0 ), fm = closed( K, D-h, X0 ); double w = 0.;
    for( size_t j = 0; o1 && o2 && Gf.size() == 1 && Ga.size() == 1 && j < 5; ++j )
      w = std::max( { w, std::fabs( Gf[0][j] - ( fp[j]-fm[j] )/( 2.*h ) ), std::fabs( Ga[0][j] - ( fp[j]-fm[j] )/( 2.*h ) ) } );
    check( o1 && o2 && Gf.size() == 1 && Ga.size() == 1 && w < 1e-7, "dose, ONLY d registered: forward and adjoint dF/dd == closed form", ( o1 && o2 )? w: 1. ); }
  { Model mi, m0; bool const s1 = build( mi, 1 ), s2 = build( m0, 2 );
    bool ok = s1 && s2 && mi.I.solve_fsens( P_of( mi, K, D, X0 ) ) == ODESLVS_CVODES::STATUS::NORMAL && m0.I.solve_fsens( P_of( m0, K, D, X0 ) ) == ODESLVS_CVODES::STATUS::NORMAL;
    double w = 0.; if( ok ){ auto const a = mi.I.val_function(), b = m0.I.val_function(); for( size_t j = 0; j < 5; ++j ) w = std::max( w, std::fabs( a[j]-b[j] ) ); }
    check( ok && w < 1e-9, "IDENTITY transition == no transition (values)", w ); }
  { Model m; std::ostringstream os; bool const su = build( m, 3 ); check( !su && m.I.extract_error().find( "EXPLICIT" ) != std::string::npos, "an IMPLICIT map (right = 2x) is REFUSED: " + m.I.extract_error().substr( 0, 60 ) ); }
  std::printf( "\n  TRANS_ode: %d passed, %d failed -- %s\n", npass, nfail, nfail? "FAILURES": "ALL PASS" );
  return nfail? 1: 0;
}
