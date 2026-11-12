// TRANS_ode.cpp -- ODESLV transitions (add_transition), against a CLOSED FORM: an impulsive dose.
//   dx/dt = -k x, x(0) = x0;  dy/dt = x, y(0) = 0;   at tau = 0.4:  x(tau+) = x(tau-) + d   (y NOT mentioned: continuous)
//   x(tau-) = x0 e^{-k tau};  x(tau+) = x(tau-) + d;  x(tf) = x(tau+) e^{-k (tf - tau)};  y(tf) = int_0^tf x dt
// Checks: values and FORWARD gradients w.r.t. (k, d, x0) against the closed form; the IDENTITY transition reproduces
// the model without one; ADJOINT gradients (incl. a registered subset of controls) == closed form; implicit maps REFUSED;
// FFModel::fdiff carries the transition (s(tau+) = s(tau-) + d_d), in the returning form before and after setup, in the
// legacy layout, and with a LIFTED output; as a DAG operation, the SYMBOLIC derivative (SYMDIFF) == solve_fsens.
#include <cmath>
#include <cstdio>
#include <sstream>
#include <iostream>
#include <string>
#include <vector>
#include "ffde.hpp"
using namespace mc;
static int npass = 0, nfail = 0;
static void check( bool c, std::string const& w, double v = -1. ){ std::printf( "  %s  %-66s", c? "PASS": "FAIL", w.c_str() ); if( v >= 0. ) std::printf( " (%.2e)", v ); std::printf( "\n" ); c? ++npass: ++nfail; }
double const TAU = 0.4, TF = 1.0, K = 1.3, D = 0.7, X0 = 1.1;
static std::vector<double> closed( double k, double d, double x0 ){
  double const xm = x0*std::exp( -k*TAU ), xp = xm + d, xf = xp*std::exp( -k*( TF-TAU ) );
  double const I = x0*( 1.-std::exp( -k*TAU ) )/k + xp*( 1.-std::exp( -k*( TF-TAU ) ) )/k;
  return { xm, xp, xf, I, I }; }                                            // x(tau-), x(tau+), x(tf), y(tf), int x
struct Model { FFGraph G; ODESLVS_CVODES I{ &G }; FFVar t, x, y, k, d, x0; };
static bool build( Model& m, int kind, bool only_d = false,     // kind 0: dose; 1: identity; 2: no transition; 3: implicit; only_d: register d alone
                   bool lift = false, bool dosetup = true ){        // lift: a sixth output [x(0.25)]^2 + d, lifted at setup; dosetup: false leaves setup() to the caller
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
  if( lift ) m.I.add_output( sqr( OpE( m.x, m.t, 0.25 ) ) + m.d );
  m.I.options.SOLVE.VERBOSE = false; m.I.options.SOLVE.ATOL = m.I.options.SOLVE.ATOLS = m.I.options.SOLVE.ATOLB = 1e-12; m.I.options.SOLVE.RTOL = m.I.options.SOLVE.RTOLS = m.I.options.SOLVE.RTOLB = 1e-12;
  if( only_d ) m.I.register_control( m.d );                                 // registered BEFORE setup
  if( !dosetup ) return true;
  std::ostringstream os; auto* old = std::cerr.rdbuf( os.rdbuf() ); bool const ok = m.I.setup(); std::cerr.rdbuf( old );
  return ok; }
static std::vector<double> P_of( Model& m, double k, double d, double x0 ){
  std::vector<double> P( m.I.np() ); P[ m.I.parameter_index( m.k )[0] ] = k; P[ m.I.parameter_index( m.d )[0] ] = d; P[ m.I.parameter_index( m.x0 )[0] ] = x0; return P; }
// FFModel::fdiff of the dose model, every direction (k, d, x0) seeded with fdiff_seed_all: max |product - solve_fsens| over
// values and gradients (1 on any failure).  before: the returning form is called BEFORE setup (setup() is then the caller's).
static double gateaux( Model& m, bool before, size_t& ntr ){
  std::vector<std::vector<FFVar>> vDU, vS;  std::string err;  ODESLVS_CVODES* S = nullptr;  ntr = 0;
  std::vector<FFVar> const U{ m.k, m.d, m.x0 };
  if( before ) S = m.I.fdiff( U, 0, vDU, vS, err );
  std::ostringstream os;  auto* old = std::cerr.rdbuf( os.rdbuf() );
  bool const su = m.I.is_setup() || m.I.setup();
  if( !before && su ) S = m.I.fdiff( U, 0, vDU, vS, err );
  bool const ss = su && S && S->setup();  std::cerr.rdbuf( old );
  if( !ss ){ delete S;  return 1.; }
  auto const P = P_of( m, K, D, X0 );
  bool const ok = m.I.solve_fsens( P ) == ODESLVS_CVODES::STATUS::NORMAL;  auto const F = m.I.val_function();  auto const G = m.I.val_function_gradient();
  bool const os2 = ok && S->solve( FFModel::fdiff_seed_all( *S, vDU, m.I, P ) ) == ODESLVS_CVODES::STATUS::NORMAL;
  double w = 1.;  ntr = S->var_transition().size();
  if( os2 ){ auto const& FS = S->val_function();  size_t const nf = m.I.nf();  FFVar const* pv[3] = { &m.k, &m.d, &m.x0 };  w = 0.;
    for( size_t j = 0; j < nf; ++j ){ w = std::max( w, std::fabs( FS[j] - F[j] ) );
      for( size_t q = 0; q < 3; ++q ) w = std::max( w, std::fabs( FS[(q+1)*nf+j] - G[ m.I.parameter_index( *pv[q] )[0] ][j] ) ); } }
  delete S;  return w; }
// the LEGACY layout (after setup, every parameter): f[i*nf+j] = d f_j / d p_i, against solve_fsens
static double legacy( Model& m ){
  std::string err;  ODESLVS_CVODES* L = m.I.fdiff( m.I.var_parameter().size(), m.I.var_parameter().data(), err );
  std::ostringstream os;  auto* old = std::cerr.rdbuf( os.rdbuf() );  bool const sl = L && L->setup();  std::cerr.rdbuf( old );
  if( !sl ){ delete L;  return 1.; }
  auto const P = P_of( m, K, D, X0 );
  bool const ok = m.I.solve_fsens( P ) == ODESLVS_CVODES::STATUS::NORMAL;  auto const G = m.I.val_function_gradient();
  std::vector<double> PL( L->np(), 0. );
  for( auto const& [w, dm] : L->var_declared_input() ) for( size_t j : L->parameter_index( w ) )
    PL[j] = w.id() == m.k.id()? K: w.id() == m.d.id()? D: w.id() == m.x0.id()? X0: 0.;
  bool const ol = ok && L->solve( PL ) == ODESLVS_CVODES::STATUS::NORMAL;  double e = 1.;
  if( ol ){ auto const& FL = L->val_function();  size_t const nf = m.I.nf();  e = 0.;
    for( size_t i = 0; i < m.I.np(); ++i ) for( size_t j = 0; j < nf; ++j ) e = std::max( e, std::fabs( FL[i*nf+j] - G[i][j] ) ); }
  delete L;  return e; }
int main(){
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
  { Model m; bool const su = build( m, 0 );                                   // the dose model as a DAG operation (FFODESLV)
    double wn = 1., ws = 1.;
    if( su ){
      auto const P = P_of( m, K, D, X0 );
      bool const ok = m.I.solve_fsens( P ) == ODESLVS_CVODES::STATUS::NORMAL;  auto const G = m.I.val_function_gradient();
      FFGraph NLP;  std::vector<FFVar> V( m.I.np() );  for( auto& v : V ) v.set( &NLP );
      std::vector<FFModel::InputArg> arg;
      for( auto const& [w,dom] : m.I.var_declared_input() ){ std::vector<FFVar> dv; for( size_t i : m.I.parameter_index( w ) ) dv.push_back( V[i] ); arg.emplace_back( w, dv ); }
      FFODESLV OpODE;  std::vector<FFVar> const F = OpODE( arg, &m.I, FFODESLV::SHALLOW );
      FFODESLV::options.SYMDIFF.clear();                                      // NUMERIC: the forward sensitivities, across the jump
      try{ auto const dF = NLP.FAD( F, V );  std::vector<double> dv( dF.size() );  NLP.eval( dF, dv, V, P );
           wn = 0.;  for( size_t j = 0; ok && j < F.size(); ++j ) for( size_t i = 0; i < V.size(); ++i ) wn = std::max( wn, std::fabs( dv[j*V.size()+i] - G[i][j] ) ); }
      catch( std::exception const& ){ wn = 1.; }
      FFODESLV::options.SYMDIFF = V;                                          // SYMBOLIC: fdiff, the transition carried into the product
      try{ auto const dF = NLP.FAD( F, V );  std::vector<double> dv( dF.size() );  NLP.eval( dF, dv, V, P );
           ws = 0.;  for( size_t j = 0; ok && j < F.size(); ++j ) for( size_t i = 0; i < V.size(); ++i ) ws = std::max( ws, std::fabs( dv[j*V.size()+i] - G[i][j] ) ); }
      catch( std::exception const& e ){ ws = 1.;  std::printf( "  SYMDIFF threw: %s\n", e.what() ); }
      FFODESLV::options.SYMDIFF.clear(); }
    check( su && wn < 1e-9, "as a DAG operation: the NUMERIC derivative (SYMDIFF empty) == solve_fsens", su? wn: 1. );
    check( su && ws < 1e-9, "as a DAG operation: the SYMBOLIC derivative (SYMDIFF = all) == solve_fsens (was refused)", su? ws: 1. ); }
  for( bool lift : { false, true } ){                                         // FFModel::fdiff carries the transition: the product holds
    // it and its derivative per direction (1 + 3), and lifts its own nonlinear outputs (the original's and 3 derivatives) at its setup
    std::string const tag = lift? "dose + LIFTED output, ": "dose, ";
    { Model m;  build( m, 0, false, lift, false );  size_t nt;  double const w = gateaux( m, true, nt );
      check( w < 1e-9 && nt == ( lift? 8u: 4u ), tag + "fdiff (Gateaux form) BEFORE setup == solve_fsens", w ); }
    { Model m;  bool const su = build( m, 0, false, lift );  size_t nt;  double const w = su? gateaux( m, false, nt ): 1.;
      check( su && w < 1e-9 && nt == ( lift? 8u: 4u ), tag + "fdiff (Gateaux form) AFTER setup == solve_fsens", w ); }
    { Model m;  bool const su = build( m, 0, false, lift );  double const w = su? legacy( m ): 1.;
      check( su && w < 1e-9, tag + "fdiff, LEGACY layout (every parameter) == solve_fsens", w ); } }
  std::printf( "\n  TRANS_ode: %d passed, %d failed -- %s\n", npass, nfail, nfail? "FAILURES": "ALL PASS" );
  return nfail? 1: 0;
}
