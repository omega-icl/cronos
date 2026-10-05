// TRANS_eval.cpp -- ODESLV transitions declared through EVALUATIONS, with no domain or time point given:
//   add_transition( OpE( x1, t, 0.3 ) + d,         OpE( x1, t, 0.3, FFDom::PLUS ) )
//   add_transition( OpE( x2 + c*x1*x1, t, 0.5 ),   OpE( x2, t, 0.5, FFDom::PLUS ) )
// tau and the direction are read off the evaluations; left is at tau^- (MINUS, the default), right at tau^+ (PLUS).
// Model and closed form as TRANS_multi.  Checks: this form is BITWISE identical to the explicit form (values, forward
// and adjoint gradients); it matches the closed form; var_transition() reports the inferred tau and direction; and
// five malformed declarations are REFUSED, each with its reason.  The SAME model through OCFESLV (monolithic and
// marching; mesh {0, 0.3, 0.5, 1}): the evaluation form is BITWISE identical to OCFESLV's explicit form, and matches
// the closed form and ODESLV (three-way agreement).
#include <cmath>
#include <cstdio>
#include <sstream>
#include <iostream>
#include <string>
#include <vector>
#include "ffode.hpp"
#include "ocfeslv.hpp"
using namespace mc;
static int npass = 0, nfail = 0;
static void check( bool c, std::string const& w, double v = -1. ){ std::printf( "  %s  %-70s", c? "PASS": "FAIL", w.c_str() ); if( v >= 0. ) std::printf( " (%.2e)", v ); std::printf( "\n" ); c? ++npass: ++nfail; }
double const T1 = 0.3, T2 = 0.5;
static std::vector<double> closed( std::vector<double> const& P ){          // P = a, b, c, d, x10, x20
  double const a = P[0], b = P[1], c = P[2], d = P[3], x10 = P[4], x20 = P[5];
  double const x1m = x10*std::exp( -a*T1 ), x1p = x1m + d;
  auto x1 = [&]( double t ){ return t < T1? x10*std::exp( -a*t ): x1p*std::exp( -a*( t-T1 ) ); };
  double const x2m = x20*std::exp( -b*T2 ), x2p = x2m + c*x1( T2 )*x1( T2 );
  double const I1 = x10*( 1.-std::exp( -a*T1 ) )/a + x1p*( 1.-std::exp( -a*( 1.-T1 ) ) )/a;
  double const I2 = x20*( 1.-std::exp( -b*T2 ) )/b + x2p*( 1.-std::exp( -b*( 1.-T2 ) ) )/b;
  return { x1( 1. ), x2p*std::exp( -b*( 1.-T2 ) ), x2m, x2p, I1 + I2 }; }
struct Run { bool ok = false; std::string msg; std::vector<double> F; std::vector<std::vector<double>> Gf, Ga; std::vector<std::pair<double,bool>> trn; };
//! form 0: explicit (dom, tau given); 1: through evaluations; 10+k: malformed declaration k
static Run run( int form, std::vector<double> const& P0 ){
  Run R; FFGraph G; ODESLVS_CVODES I( &G ); FFPartial OpP; FFEval OpE; FFIntegral OpI;
  FFVar t = G.add_var( "t" ), x1 = G.add_var( "x1(t)" ), x2 = G.add_var( "x2(t)" );
  char const* pn[6] = { "a", "b", "c", "d", "x10", "x20" }; std::vector<FFVar> p; for( auto n : pn ) p.push_back( G.add_var( n ) );
  I.add_domain( t, FFDom( std::vector<double>{ 0., 0.5, 1. }, FFDom::LGR, 4 ) ); I.set_evolution_domain( t );
  I.add_state( x1, {t} ); I.add_state( x2, {t} ); for( auto const& v : p ) I.add_input( v );
  I.update_ref( x1, 1. ); I.update_ref( x2, 0.5 );
  int const T_INT = FFDom::ALL - FFDom::LB;
  I.add_equation( OpP( x1, t ) + p[0]*x1, {t}, {T_INT}, FFModel::EqnOptions( FFModel::EqnRole::INTERIOR, 0 ) );
  I.add_equation( OpP( x2, t ) + p[1]*x2, {t}, {T_INT}, FFModel::EqnOptions( FFModel::EqnRole::INTERIOR, 0 ) );
  I.add_equation( x1 - p[4], {t}, {FFDom::LB}, FFModel::EqnOptions( FFModel::EqnRole::INITIAL, 0 ) );
  I.add_equation( x2 - p[5], {t}, {FFDom::LB}, FFModel::EqnOptions( FFModel::EqnRole::INITIAL, 0 ) );
  switch( form ){
    case 0:  I.add_transition( x1 + p[3], x1, t, T1 );  I.add_transition( x2 + p[2]*x1*x1, x2, t, T2 );  break;
    case 1:  I.add_transition( OpE( x1, t, T1 ) + p[3], OpE( x1, t, T1, FFDom::PLUS ) );
             I.add_transition( OpE( x2 + p[2]*x1*x1, t, T2 ), OpE( x2, t, T2, FFDom::PLUS ) );  break;
    case 10: I.add_transition( OpE( x1, t, T1 ) + p[3], OpE( x1, t, T1 ) );  break;                        // right without PLUS
    case 11: I.add_transition( OpE( x1, t, T1, FFDom::PLUS ) + p[3], OpE( x1, t, T1, FFDom::PLUS ) );  break; // left with PLUS
    case 12: I.add_transition( OpE( x1, t, T1 ) + p[3], OpE( x1, t, T2, FFDom::PLUS ) );  break;            // two points
    case 13: I.add_transition( x1 + p[3], OpE( x1, t, T1, FFDom::PLUS ) );  break;                          // bare state
    case 14: I.add_transition( OpE( OpP( x1, t ), t, T1 ), OpE( x1, t, T1, FFDom::PLUS ) );  break;         // derivative inside
  }
  I.add_output( OpE( x1, t, 1. ) ); I.add_output( OpE( x2, t, 1. ) );
  I.add_output( OpE( x2, t, T2, FFDom::MINUS ) ); I.add_output( OpE( x2, t, T2, FFDom::PLUS ) );
  I.add_output( OpI( std::vector<FFVar>{ x1 + x2 }, { t } )[0] );
  I.options.DISPLAY = 0; I.options.ATOL = I.options.ATOLS = I.options.ATOLB = 1e-12; I.options.RTOL = I.options.RTOLS = I.options.RTOLB = 1e-12;
  std::ostringstream os; auto* old = std::cerr.rdbuf( os.rdbuf() ); bool const su = I.setup(); std::cerr.rdbuf( old ); R.msg = os.str();
  if( !su ) return R;
  for( auto const& tr : I.var_transition() ) R.trn.push_back( { tr.tau, tr.dom.id() == t.id() || tr.dom.name() == t.name() } );
  std::vector<double> Q( I.np() ); std::vector<size_t> ix; for( size_t q = 0; q < 6; ++q ){ ix.push_back( I.parameter_index( p[q] )[0] ); Q[ix[q]] = P0[q]; }
  bool const of = I.solve_fsens( Q ) == ODESLVS_CVODES::STATUS::NORMAL; R.F = I.val_function(); auto const Gf = I.val_function_gradient();
  bool const oa = I.solve_asens( Q ) == ODESLVS_CVODES::STATUS::NORMAL;     auto const Ga = I.val_function_gradient();
  if( !of || !oa ) return R;
  for( size_t q = 0; q < 6; ++q ){ R.Gf.push_back( Gf[ix[q]] ); R.Ga.push_back( Ga[ix[q]] ); }     // rows ordered a..x20
  R.ok = true; return R; }

//! The model through OCFESLV, declared in form 0 (explicit) or 1 (evaluations): values, control Jacobians (5 x 6)
struct OcRun { bool ok = false; std::vector<double> F, Jf, Ja; };
static OcRun run_oc( int form, bool march, std::vector<double> const& P0 ){
  OcRun R; FFGraph G; OCFESLV I( &G ); FFPartial OpP; FFEval OpE; FFIntegral OpI;
  FFVar t = G.add_var( "t" ), x1 = G.add_var( "x1(t)" ), x2 = G.add_var( "x2(t)" );
  char const* pn[6] = { "a", "b", "c", "d", "x10", "x20" }; std::vector<FFVar> p; for( auto n : pn ) p.push_back( G.add_var( n ) );
  I.add_domain( t, FFDom( std::vector<double>{ 0., T1, T2, 1. }, FFDom::LGR, 8 ) ); I.set_evolution_domain( t );
  I.add_state( x1, {t} ); I.add_state( x2, {t} ); for( auto const& v : p ) I.add_input( v, {} );
  I.update_ref( x1, 1. ); I.update_ref( x2, 0.5 );
  int const T_INT = FFDom::ALL - FFDom::LB;
  I.add_equation( OpP( x1, t ) + p[0]*x1, {t}, {T_INT}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  I.add_equation( OpP( x2, t ) + p[1]*x2, {t}, {T_INT}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  I.add_equation( x1 - p[4], {t}, {FFDom::LB}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INITIAL, 0 ) );
  I.add_equation( x2 - p[5], {t}, {FFDom::LB}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INITIAL, 0 ) );
  if( form == 0 ){ I.add_transition( x1 + p[3], x1, t, T1 );  I.add_transition( x2 + p[2]*x1*x1, x2, t, T2 ); }
  else{ I.add_transition( OpE( x1, t, T1 ) + p[3], OpE( x1, t, T1, FFDom::PLUS ) );
        I.add_transition( OpE( x2 + p[2]*x1*x1, t, T2 ), OpE( x2, t, T2, FFDom::PLUS ) ); }
  I.add_output( OpE( x1, t, 1. ) ); I.add_output( OpE( x2, t, 1. ) );
  I.add_output( OpE( x2, t, T2, FFDom::MINUS ) ); I.add_output( OpE( x2, t, T2, FFDom::PLUS ) );
  I.add_output( OpI( std::vector<FFVar>{ x1 + x2 }, { t } )[0] );
  I.options.SOLVE.MARCHING = march; I.options.SOLVE.RES_TOL = 1e-12; I.options.DISPLAY_LEVEL = 0;
  std::ostringstream os; auto* old = std::cerr.rdbuf( os.rdbuf() );
  std::vector<double> xv, inp; bool ok = I.setup() && I.init( xv, inp, nullptr );
  if( ok ) for( size_t q = 0; q < 6; ++q ){ I.set_input_values( p[q], { P0[q] }, inp.data() ); I.register_control( p[q] ); }
  std::vector<double> y1( xv ), y2( xv ), y3( xv );
  if( ok ) ok = I.solve( y1.data(), inp.data(), nullptr ).converged;  if( ok ) R.F = I.val_functions();
  if( ok ) ok = I.solve_fsens( y2.data(), inp.data(), nullptr );     if( ok ) R.Jf = I.sens_jacobian();
  if( ok ) ok = I.solve_asens( y3.data(), inp.data(), nullptr );     if( ok ) R.Ja = I.sens_jacobian();
  std::cerr.rdbuf( old );
  R.ok = ok && R.F.size() == 5 && R.Jf.size() == 30 && R.Ja.size() == 30;
  return R; }

int main(){
  std::vector<double> const P0{ 1.3, 0.8, 0.6, 0.7, 1.1, 0.4 };
  Run const A = run( 0, P0 ), B = run( 1, P0 );
  check( A.ok, "explicit form (dom, tau given): setup and solves" );
  check( B.ok, "EVALUATION form (no dom, no tau): setup and solves" + ( B.ok? std::string(): ": " + B.msg.substr( 0, 60 ) ) );
  if( A.ok && B.ok ){
    check( B.trn.size() == 2 && B.trn[0].first == T1 && B.trn[1].first == T2 && B.trn[0].second && B.trn[1].second,
           "var_transition(): tau (0.3, 0.5) and the direction INFERRED from the evaluations" );
    check( A.F == B.F,   "values BITWISE identical to the explicit form" );
    check( A.Gf == B.Gf, "FORWARD gradients BITWISE identical to the explicit form" );
    check( A.Ga == B.Ga, "ADJOINT gradients BITWISE identical to the explicit form" );
    auto const F0 = closed( P0 ); double ev = 0.; for( size_t j = 0; j < 5; ++j ) ev = std::max( ev, std::fabs( B.F[j]-F0[j] ) );
    check( ev < 1e-9, "values == closed form", ev );
    double eg = 0.;
    for( size_t q = 0; q < 6; ++q ){ auto Pp = P0, Pm = P0; double const h = 1e-6*std::max( 1., std::fabs( P0[q] ) ); Pp[q] += h; Pm[q] -= h;
      auto const fp = closed( Pp ), fm = closed( Pm );
      for( size_t j = 0; j < 5; ++j ){ double const D = ( fp[j]-fm[j] )/( 2.*h ); eg = std::max( { eg, std::fabs( B.Gf[q][j]-D ), std::fabs( B.Ga[q][j]-D ) } ); } }
    check( eg < 1e-7, "forward and adjoint gradients w.r.t. all 6 parameters == closed form", eg );
  }
  struct Bad { int form; char const* what; char const* expect; };
  for( Bad b : { Bad{ 10, "right-hand evaluation without PLUS", "give its evaluations FFDom::PLUS" },
                 Bad{ 11, "left-hand evaluation with PLUS",     "an evaluation there is FFDom::PLUS" },
                 Bad{ 12, "evaluations at two points",          "DIFFERENT points" },
                 Bad{ 13, "a state outside any evaluation",     "OUTSIDE an evaluation" },
                 Bad{ 14, "a derivative inside an evaluation",  "POINTWISE" } } ){
    Run const R = run( b.form, P0 );
    check( !R.ok && R.msg.find( b.expect ) != std::string::npos, std::string( "REFUSED, " ) + b.what );
  }
  for( bool march : { false, true } ){                   // ---- the SAME model through OCFESLV
    std::string const tag = march? "OCFESLV marching  : ": "OCFESLV monolithic: ";
    OcRun const X = run_oc( 0, march, P0 ), E = run_oc( 1, march, P0 );
    if( !X.ok || !E.ok ){ check( false, tag + "setup and solves, both forms" ); continue; }
    check( X.F == E.F && X.Jf == E.Jf && X.Ja == E.Ja, tag + "EVALUATION form BITWISE identical to the explicit form" );
    auto const F0 = closed( P0 ); double ev = 0., e3 = A.ok? 0.: 1.;
    for( size_t j = 0; j < 5; ++j ){ ev = std::max( ev, std::fabs( E.F[j]-F0[j] ) ); if( A.ok ) e3 = std::max( e3, std::fabs( E.F[j]-A.F[j] ) ); }
    check( ev < 1e-9, tag + "values == closed form", ev );
    check( e3 < 1e-9, tag + "values == ODESLV's (three-way agreement)", e3 );
    double eg = 0.;
    for( size_t q = 0; q < 6; ++q ){ auto Pp = P0, Pm = P0; double const h = 1e-6*std::max( 1., std::fabs( P0[q] ) ); Pp[q] += h; Pm[q] -= h;
      auto const fp = closed( Pp ), fm = closed( Pm );
      for( size_t j = 0; j < 5; ++j ){ double const Dq = ( fp[j]-fm[j] )/( 2.*h ); eg = std::max( { eg, std::fabs( E.Jf[j*6+q]-Dq ), std::fabs( E.Ja[j*6+q]-Dq ) } ); } }
    check( eg < 1e-7, tag + "forward and adjoint gradients (6 parameters) == closed form", eg );
  }
  std::printf( "\n  TRANS_eval: %d passed, %d failed -- %s\n", npass, nfail, nfail? "FAILURES": "ALL PASS" );
  return nfail? 1: 0;
}
