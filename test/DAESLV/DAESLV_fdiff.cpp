// DAESLV_fdiff -- the solver as a DIRECTIONAL DERIVATIVE of itself.
// =====================================================================================================
// fdiff() returns a NEW DAESLV_IDAS whose model is the derivative of this one along the given
// parameters.  Solving THAT solver gives the same numbers solve_fsens() gives, by a different route:
// fsens differentiates the integration, fdiff differentiates the MODEL and integrates the result.
// Agreement between the two is the gate, and it is what FFDAESLV will rest on.
//
//   x' = -k x , x(0) = a ;  f = int_0^1 x dt
//   d f/da = (1-e^-k)/k ,  d f/dk = a[ k e^-k - (1-e^-k) ]/k^2
//
//   A. fdiff over both parameters (the legacy layout): the derivatives appear among the product's outputs.
//   B. TRANSITIONS are carried (refused until 12 Nov 2026): the product holds each transition and, per direction, its
//      derivative -- the same rows on the sensitivity states and the direction inputs -- so the stop the restart solves
//      is the differentiated stop.  x' = -k x, x(0) = x0, with
//        explicit   x(tau+) = x(tau-) + d at tau = 0.4, OFF the grid;
//        implicit   x(tau+)^2 = x(tau-) + d at tau = 0.5;
//        lifted     no transition declared, an output [x(0.4)]^2 + d, whose lift declares one at setup;
//        algebraic  the explicit dose with z = x^2 algebraic, re-solved at the stop, observed at tau+, tf and in integral;
//        evaluated  the explicit dose declared as add_transition( E(x, 0.4) + d, E(x, 0.4, PLUS) ).
//      Each: the Gateaux form (fdiff called before, then after setup; fdiff_seed_all) and the legacy layout == solve_fsens.
//   C. METHOD OF LINES with a transition and a REGISTERED control (the registry was lost to the swap of declarations
//      until 12 Nov 2026): the direction input, fdiff_seed_all and the legacy layout == solve_fsens.
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>
#include "daeslv_idas.hpp"

using namespace mc;

static int n_pass = 0, n_fail = 0;
static void check( std::string const& w, double got, double want, double tol )
{
  double const d = std::fabs( got - want );
  ( d <= tol )? ++n_pass: ++n_fail;
  std::cout << "    " << std::left << std::setw( 52 ) << w << std::right
            << " got=" << std::scientific << std::setprecision(10) << got
            << " want=" << want << " |d|=" << std::setprecision(2) << d
            << ( d <= tol? "  PASS\n": "  FAIL\n" );
}
static void check_bool( std::string const& w, bool ok )
{
  ok? ++n_pass: ++n_fail;
  std::cout << "    " << std::left << std::setw( 52 ) << w << ( ok? "  PASS\n": "  FAIL\n" );
}

static bool normal( DAESLV_IDAS::STATUS s ){ return s == DAESLV_IDAS::NORMAL; }

// B: one model with a transition; kind 0 explicit off-grid, 1 implicit, 2 lifted output, 3 explicit with an algebraic state,
// 4 the explicit one declared through EVALUATIONS (its derivative is declared through the sensitivities' evaluations)
static void transitions( char const* name, int kind )
{
  std::cout << "\n--- B. fdiff with a transition: " << name << "\n";
  FFGraph G;  DAESLV_IDAS I( &G );
  FFVar t = G.add_var( "t" ), x = G.add_var( "x(t)" ), z = G.add_var( "z(t)" ), k = G.add_var( "k" ), d = G.add_var( "d" ), x0 = G.add_var( "x0" );
  FFPartial OpP;  FFEval OpE;  FFIntegral OpI;
  double const tau = kind == 1? 0.5: 0.4;
  I.add_domain( t, FFDom( 0., 1., 2, FFDom::LGR, 4 ) );  I.set_evolution_domain( t );
  I.add_state( x, {t} );  I.add_input( k );  I.add_input( d );  I.add_input( x0 );  I.update_ref( x, 1. );
  I.add_equation( OpP( x, t ) + k*x, {t}, {FFDom::ALL - FFDom::LB}, FFModel::EqnOptions( FFModel::EqnRole::INTERIOR, 0 ) );
  I.add_equation( x - x0,            {t}, {FFDom::LB},             FFModel::EqnOptions( FFModel::EqnRole::INITIAL, 0 ) );
  if( kind == 3 ){ I.add_state( z, {t} );  I.update_ref( z, 1. );
    I.add_equation( z - x*x, {t}, {FFDom::ALL}, FFModel::EqnOptions( FFModel::EqnRole::INTERIOR, 0 ) ); }
  if( kind == 0 || kind == 3 ) I.add_transition( x + d, x, t, tau );
  if( kind == 4 )              I.add_transition( OpE( x, t, tau ) + d, OpE( x, t, tau, FFDom::PLUS ) );
  if( kind == 1 )              I.add_transition( x + d, x*x, t, tau );
  if( kind == 2 )      I.add_output( sqr( OpE( x, t, tau ) ) + d );
  else if( kind == 3 ){ I.add_output( OpE( z, t, tau, FFDom::PLUS ) );  I.add_output( OpE( z, t, 1. ) );  I.add_output( OpI( std::vector<FFVar>{ z }, { t } )[0] ); }
  else{                 I.add_output( OpE( x, t, tau, FFDom::PLUS ) );  I.add_output( OpE( x, t, 1. ) ); }
  I.options.DISPLAY_LEVEL = 0;  I.options.SOLVE.VERBOSE = 0;  I.options.SOLVE.ATOL = I.options.SOLVE.RTOL = I.options.SOLVE.ATOLS = I.options.SOLVE.RTOLS = 1e-12;
  std::vector<FFVar> const U{ k, d, x0 };  double const pv[3] = { 1.3, 0.7, 1.1 };
  for( bool before : { true, false } ){
    std::string const tag = before? "fdiff BEFORE setup: ": "fdiff AFTER setup: ";
    std::vector<std::vector<FFVar>> vDU, vS;  std::string err;  DAESLV_IDAS* S = nullptr;
    if( before ) S = I.fdiff( U, 0, vDU, vS, err );
    std::ostringstream os;  auto* old = std::cerr.rdbuf( os.rdbuf() );
    bool const su = I.is_setup() || I.setup();
    if( !before && su ) S = I.fdiff( U, 0, vDU, vS, err );
    bool const ss = su && S && S->setup();  std::cerr.rdbuf( old );
    check_bool( tag + "the product sets up", ss );
    if( !ss ){ std::cout << "      " << ( !su? I.extract_error(): !S? err: S->extract_error() ) << "\n";  delete S;  continue; }
    std::vector<double> P( I.np(), 0. );  for( size_t q = 0; q < 3; ++q ) P[ I.parameter_index( U[q] )[0] ] = pv[q];
    bool const ok = normal( I.solve_fsens( P ) );  auto const F = I.val_function();  auto const Gr = I.val_function_gradient();
    bool const sv = ok && normal( S->solve( FFModel::fdiff_seed_all( *S, vDU, I, P ) ) );
    double wf = 1., wg = 1.;
    if( sv ){ auto const& FS = S->val_function();  size_t const nf = I.nf();  wf = wg = 0.;
      for( size_t j = 0; j < nf; ++j ){ wf = std::max( wf, std::fabs( FS[j] - F[j] ) );
        for( size_t q = 0; q < 3; ++q ) wg = std::max( wg, std::fabs( FS[(q+1)*nf+j] - Gr[ I.parameter_index( U[q] )[0] ][j] ) ); } }
    check( tag + "values", wf, 0., 1e-9 );
    check( tag + "d/d(k,d,x0) == fsens", wg, 0., 1e-8 );
    delete S; }
  { std::string err;  DAESLV_IDAS* L = I.fdiff( I.var_parameter().size(), I.var_parameter().data(), err );      // the LEGACY layout
    std::ostringstream os;  auto* old = std::cerr.rdbuf( os.rdbuf() );  bool const sl = L && L->setup();  std::cerr.rdbuf( old );
    check_bool( "legacy layout: the product sets up", sl );
    if( !sl ){ std::cout << "      " << ( L? L->extract_error(): err ) << "\n";  delete L;  return; }
    std::vector<double> P( I.np(), 0. );  for( size_t q = 0; q < 3; ++q ) P[ I.parameter_index( U[q] )[0] ] = pv[q];
    bool const ok = normal( I.solve_fsens( P ) );  auto const Gr = I.val_function_gradient();
    std::vector<double> PL( L->np(), 0. );
    for( auto const& [w, dm] : L->var_declared_input() ) for( size_t j : L->parameter_index( w ) )
      PL[j] = w.id() == k.id()? pv[0]: w.id() == d.id()? pv[1]: w.id() == x0.id()? pv[2]: 0.;
    bool const lv = ok && normal( L->solve( PL ) );  double wl = 1.;
    if( lv ){ auto const& FL = L->val_function();  size_t const nf = I.nf();  wl = 0.;
      for( size_t i = 0; i < I.np(); ++i ) for( size_t j = 0; j < nf; ++j ) wl = std::max( wl, std::fabs( FL[i*nf+j] - Gr[i][j] ) ); }
    check( "legacy layout: f[i*nf+j] == fsens", wl, 0., 1e-8 );
    delete L; }
}

// C: a heated rod by the method of lines, halved at t = 0.3, its source p a REGISTERED control
static void lines()
{
  std::cout << "\n--- C. fdiff under the method of lines: a transition and a registered control\n";
  FFGraph G;  DAESLV_IDAS S( &G );
  FFVar t = G.add_var("t"), z = G.add_var("z"), u = G.add_var("u"), p = G.add_var("p");
  S.add_domain( t, FFDom( 0., 1., 2, FFDom::LGR, 3 ) );  S.add_domain( z, FFDom( 0., 1., 3, FFDom::LGL, 5 ) );
  S.set_evolution_domain( t );  S.add_state( u, {t,z} );  S.update_ref( u, 1. );  S.add_input( p, {} );
  FFPartial OpP;  FFIntegral OpI;
  int const TI = FFDom::ALL - FFDom::LB, ZI = FFDom::ALL - FFDom::LB - FFDom::UB;
  FFModel::EqnOptions io( FFModel::EqnRole::INTERIOR, 0 ), ib( FFModel::EqnRole::BOUNDARY, 0 ), ii( FFModel::EqnRole::INITIAL, 0 );
  S.add_equation( OpP(u,t) - OpP( 0.1*OpP(u,z), z ) - p*exp( -t ), {t,z}, {TI, ZI},                io );
  S.add_equation( OpP(u,z),                                       {t,z}, {TI, FFDom::LB},         ib );
  S.add_equation( OpP(u,z),                                       {t,z}, {TI, FFDom::UB},         ib );
  S.add_equation( u - p,                                          {t,z}, {FFDom::LB, FFDom::ALL}, ii );
  S.add_transition( { 0.5*u }, { u }, t, 0.3 );
  S.add_output( OpI( u, z ), {t}, {1.} );
  S.options.INTERFACE.IMPOSITION = OCFEMOL::Options::IC_STRONG;  S.options.DISPLAY_LEVEL = 0;  S.options.SOLVE.VERBOSE = 0;
  S.options.SOLVE.RTOL = S.options.SOLVE.RTOLS = 1e-10;  S.options.SOLVE.ATOL = S.options.SOLVE.ATOLS = 1e-12;
  S.register_control( p );
  std::ostringstream os;  auto* old = std::cerr.rdbuf( os.rdbuf() );  bool const su = S.setup();  std::cerr.rdbuf( old );
  check_bool( "the rod sets up, p registered", su );
  if( !su ){ std::cout << "      " << S.extract_error() << "\n";  return; }
  check_bool( "the registry reaches the semi-discrete model", S.n_control_dof() == 1 && S.sensitivity_index().size() == 1 );
  std::vector<FFModel::InputVal> const in{ { p, { 0.7 } } };
  bool const ok = normal( S.solve_fsens( in ) );
  double const J = ok? S.val_function()[0]: 0., dJ = ok && S.val_function_gradient().size() == 1? S.val_function_gradient()[0][0]: 0.;
  check_bool( "solve_fsens: one gradient row, dJ/dp", ok && S.val_function_gradient().size() == 1 );
  std::vector<std::vector<FFVar>> vDU, vS;  std::string err;
  DAESLV_IDAS* D = S.fdiff( std::vector<FFVar>{ p }, 0, vDU, vS, err );
  old = std::cerr.rdbuf( os.rdbuf() );  bool const sd = D && D->setup();  std::cerr.rdbuf( old );
  check_bool( "fdiff( {p} ): the product sets up", sd );
  if( !sd ) std::cout << "      " << ( D? D->extract_error(): err ) << "\n";
  else{
    bool const s1 = normal( D->solve( std::vector<FFModel::InputVal>{ { p, { 0.7 } }, { vDU[0][0], { 1. } } } ) );
    check( "direction by name: J",              s1? D->val_function()[0]: 1., J,  1e-9 );
    check( "direction by name: dJ/dp == fsens", s1? D->val_function()[1]: 1., dJ, 1e-7 );
    bool const s2 = normal( D->solve( FFModel::fdiff_seed_all( *D, vDU, S, std::vector<double>{ 0.7 } ) ) );
    check( "fdiff_seed_all: dJ/dp == fsens",    s2? D->val_function()[1]: 1., dJ, 1e-7 ); }
  delete D;
  std::string e2;  DAESLV_IDAS* L = S.fdiff( std::vector<FFVar>{ p }, e2 );
  old = std::cerr.rdbuf( os.rdbuf() );  bool const sl = L && L->setup();  std::cerr.rdbuf( old );
  check_bool( "legacy fdiff( {p} ): the product sets up", sl );
  if( !sl ) std::cout << "      " << ( L? L->extract_error(): e2 ) << "\n";
  else{ bool const lv = normal( L->solve( std::vector<FFModel::InputVal>{ { p, { 0.7 } } } ) );
    check( "legacy fdiff( {p} ): dJ/dp == fsens", lv? L->val_function()[0]: 1., dJ, 1e-7 ); }
  delete L;
}

int main()
{
  FFPartial OpP;  FFIntegral OpI;
  int const TI = FFDom::ALL - FFDom::LB;
  auto io = FFModel::EqnOptions( FFModel::EqnRole::INTERIOR, 0 );
  auto ii = FFModel::EqnOptions( FFModel::EqnRole::INITIAL,  0 );

  FFGraph G;  DAESLV_IDAS M;  M.set_dag( &G );
  FFVar t = G.add_var("t"), x = G.add_var("x(t)"), k = G.add_var("k"), a = G.add_var("a");
  M.add_domain( t, FFDom( 0., 1., 1, FFDom::LGR, 5 ) );
  M.add_state( x, {t} );  M.add_input( k );  M.add_input( a );
  M.set_evolution_domain( t );  M.update_ref( x, 1. );
  M.add_equation( OpP( x, t ) + k*x, {t}, {TI}, io );
  M.add_equation( x - a,             {t}, {FFDom::LB}, ii );
  M.add_output( OpI( std::vector<FFVar>{ x }, { t } )[0] );
  M.options.SOLVE.RTOL = M.options.SOLVE.ATOL = 1e-11;
  M.options.DISPLAY_LEVEL = 0;
  if( !M.setup() ){ std::cout << "setup failed: " << M.extract_error() << "\n"; return 1; }

  double const K = 1.3, A = 2.1, e = std::exp(-K);
  size_t const jk = M.parameter_index( k )[0], ja = M.parameter_index( a )[0];
  std::vector<double> P( M.np(), 0. );  P[jk] = K;  P[ja] = A;

  // the reference: forward sensitivities of the ORIGINAL solver, already gated elsewhere
  M.solve_fsens( P );
  auto const Gref = M.val_function_gradient();
  size_t const b = M.blk_fct(0).first;
  std::cout << "\n--- reference: solve_fsens on the original model\n";
  check( "d f/da = (1-e^-k)/k", Gref[ja][b], ( 1.-e )/K, 1e-7 );
  check( "d f/dk",              Gref[jk][b], A*( K*e - (1.-e) )/(K*K), 1e-7 );

  // ------------------------------------------------------------------ A: fdiff over both parameters
  { std::cout << "\n--- A. fdiff( {k,a} ): a NEW solver whose model is the derivative\n";
    std::string err;
    std::vector<FFVar> vU{ k, a };
    DAESLV_IDAS* S = M.fdiff( vU, err );
    check_bool( "fdiff returned a solver", S != nullptr );
    if( !S ){ std::cout << "    err: " << err << "\n"; }
    else{
      check_bool( "it is a DAESLV_IDAS", dynamic_cast<DAESLV_IDAS*>( S ) != nullptr );
      if( !S->setup() ){ std::cout << "    setup of the derivative model failed: " << S->extract_error() << "\n"; ++n_fail; }
      else{
        std::cout << "    derivative model: np=" << S->np() << " nf=" << S->nf()
                  << " nx=" << S->nx() << "\n";
        check_bool( "it carries more outputs than the original", S->nf() > M.nf() );
        // The derivative model's parameters are the ORIGINALS plus one DIRECTION input per
        // parameter differentiated.  Every one needs a value -- a missing d1[k] is what the dag
        // reports, not a defaulted zero.  They are filled by NAME, as ODESLV_test8 does.
        // var_declared_input() is the generic accessor (ODESLV spells it var_parameter); each input
        // may carry several DOFs, so values go in by parameter_index().
        std::vector<double> Q( S->np(), 0. );
        size_t ndir = 0;
        for( auto const& [w, dom] : S->var_declared_input() ){
          std::string const n = w.name();
          double v = 0.;
          if      ( n == "k" ) v = K;
          else if ( n == "a" ) v = A;
          else { ++ndir; v = 0.; }           // a direction input, seeded by fdiff through fix_input
          for( size_t j : S->parameter_index( w ) ) if( j < Q.size() ) Q[j] = v;
        }
        std::cout << "    inputs: " << S->var_declared_input().size()
                  << " (" << ndir << " direction input(s)), np=" << S->np() << "\n";
        if( S->solve( Q ) != DAESLV_IDAS::NORMAL ){ std::cout << "    solve failed\n"; ++n_fail; }
        else{
          auto const& F = S->val_function();
          std::cout << "    derivative outputs:";
          for( size_t i = 0; i < F.size() && i < 6; ++i ) std::cout << " " << std::scientific
                    << std::setprecision(6) << F[i];
          std::cout << "\n";
          // the directional derivatives appear among the outputs: look for the two we know
          bool got_da = false, got_dk = false;
          for( double v : F ){
            if( std::fabs( v - ( 1.-e )/K ) < 1e-6 ) got_da = true;
            if( std::fabs( v - A*( K*e - (1.-e) )/(K*K) ) < 1e-6 ) got_dk = true;
          }
          check_bool( "d f/da appears among the derivative outputs", got_da );
          check_bool( "d f/dk appears among the derivative outputs", got_dk );
        }
      }
      delete S;
    } }

  transitions( "EXPLICIT, off the grid: x(0.4+) = x(0.4-) + d", 0 );
  transitions( "IMPLICIT: x(0.5+)^2 = x(0.5-) + d", 1 );
  transitions( "the lift of an output [x(0.4)]^2 + d", 2 );
  transitions( "EXPLICIT, an algebraic z = x^2 re-solved at the stop", 3 );
  transitions( "EXPLICIT, declared through evaluations", 4 );
  lines();

  std::cout << "\n==========  DAESLV_fdiff: " << n_pass << " passed, " << n_fail << " failed -- "
            << ( n_fail? "FAILURES": "ALL PASS" ) << "  ==========\n\n";
  return n_fail? 1: 0;
}
