// AUTO_role.cpp -- resolution of the AUTO equation role and ODESLV's refusal of rows it does not use (2026-10-01).
// An equation on LB alone of the EVOLUTION domain resolves to INITIAL (it used to resolve to BOUNDARY, which ODESLV
// silently skipped: x(0) stayed 0 and solve() returned NORMAL), also when the evolution domain is declared after the
// equation; spatial LB/UB rows stay BOUNDARY.  ODESLV refuses a row with a role it does not use (BOUNDARY,
// INTERFACE, LINK, SURFACE, DIAGNOSTIC) instead of ignoring it.
#include <cmath>
#include <cstdio>
#include <sstream>
#include <string>
#include "odeslvs_cvodes.hpp"
#include "ocfeslv.hpp"
using namespace mc;
static int npass = 0, nfail = 0;
static void check( bool c, std::string const& w ){ std::printf( "  %s  %s\n", c? "PASS": "FAIL", w.c_str() ); c? ++npass: ++nfail; }
typedef FFModel::EqnRole Role;
double const P = 0.7, X1 = std::exp( -0.7 );

// x' = -p x, x(0) = 1 with the initial condition given with role `ic` (AUTO by default); evolution domain set
// before (late = false) or after (late = true) the equations; returns x(1), or NaN if setup or solve fails
static double ode( bool late, Role ic, bool diag, std::string& err, Role& resolved ){
  FFGraph G; FFPartial OpP; FFEval OpE;
  FFVar t = G.add_var( "t" ), x = G.add_var( "x" ), p = G.add_var( "p" );
  ODESLVS_CVODES S( &G );
  S.add_domain( t, FFDom( 0., 1., 2 ) );
  if( !late ) S.set_evolution_domain( t );
  S.add_state( x, {t} );  S.add_input( p );
  S.add_equation( OpP( x, t ) + p * x, {t}, {FFDom::ALL - FFDom::LB}, FFModel::EqnOptions( Role::INTERIOR ) );
  S.add_equation( x - 1., {t}, {FFDom::LB}, FFModel::EqnOptions( ic ) );
  if( diag ) S.add_equation( x - 2., {t}, {FFDom::UB}, FFModel::EqnOptions( Role::DIAGNOSTIC ) );
  if( late ) S.set_evolution_domain( t );
  resolved = S.var_equation()[1].opt->role;              // the declared record, before setup
  S.add_output( OpE( x, t, 1. ) );
  S.options.RTOL = S.options.ATOL = 1e-10;  S.options.DISPLAY = 0;  S.options.DISPLAY_LEVEL = 0;
  std::ostringstream os; auto* old = std::cerr.rdbuf( os.rdbuf() );
  bool ok = S.setup();
  std::cerr.rdbuf( old );
  err = S.extract_error() + os.str();
  if( !ok ) return std::nan( "" );
  if( S.solve( std::vector<double>{ P } ) != ODESLVS_CVODES::NORMAL ) return std::nan( "" );
  return S.val_function()[0];
}

int main(){
  std::string err;  Role r;
  double v = ode( false, Role::AUTO, false, err, r );
  check( r == Role::INITIAL, "AUTO on LB alone of the evolution domain -> INITIAL (domain set before)" );
  check( std::fabs( v - X1 ) < 1e-8, "  ODESLV: x(1) == exp(-0.7)" );
  v = ode( true, Role::AUTO, false, err, r );
  check( r == Role::INITIAL, "AUTO re-resolved -> INITIAL when the evolution domain is set AFTER the equation" );
  check( std::fabs( v - X1 ) < 1e-8, "  ODESLV: x(1) == exp(-0.7)" );
  v = ode( false, Role::BOUNDARY, false, err, r );
  check( std::isnan( v ) && err.find( "BOUNDARY" ) != std::string::npos && err.find( "does not use" ) != std::string::npos,
         "ODESLV refuses an explicit BOUNDARY row on t = 0, naming it" );
  v = ode( false, Role::INITIAL, true, err, r );
  check( std::isnan( v ) && err.find( "DIAGNOSTIC" ) != std::string::npos, "ODESLV refuses a DIAGNOSTIC row (it does not evaluate it), naming it" );

  // heat equation in OCFESLV with NO role given: the initial condition INITIAL, the spatial conditions BOUNDARY
  {
    FFGraph G; FFPartial OpP; FFEval OpE; OCFESLV S( &G );
    FFVar t = G.add_var( "t" ), z = G.add_var( "z" ), u = G.add_var( "u" );
    S.add_domain( t, FFDom( std::vector<double>{ 0., .4, 1. }, FFDom::LGR, 8 ) );
    S.add_domain( z, FFDom( 0., 1., 4, FFDom::LGL, 7 ) );
    S.set_evolution_domain( t );  S.add_state( u, {t,z} );  S.update_ref( u, .5 );
    int const TI = FFDom::ALL - FFDom::LB, ZI = FFDom::ALL - FFDom::LB - FFDom::UB;
    S.add_equation( OpP( u, t ) - .25 * OpP( u, {{z,2}} ), {t,z}, {TI,ZI} );
    S.add_equation( u, {t,z}, {TI,FFDom::LB} );
    S.add_equation( u, {t,z}, {TI,FFDom::UB} );
    S.add_equation( u - sin( M_PI * z ), {t,z}, {FFDom::LB,FFDom::ALL} );
    S.add_output( OpE( u, {{t,1},{z,1}}, {{t,1.},{z,.5}} ) );
    auto const& E = S.var_equation();
    check( E[0].opt->role == Role::INTERIOR && E[1].opt->role == Role::BOUNDARY && E[2].opt->role == Role::BOUNDARY
           && E[3].opt->role == Role::INITIAL, "OCFESLV heat: AUTO -> INTERIOR, BOUNDARY, BOUNDARY, INITIAL" );
    S.options.SOLVE.MARCHING = false;  S.options.SOLVE.RES_TOL = 1e-12;  S.options.DISPLAY_LEVEL = 0;
    bool ok = S.setup();  std::vector<double> var, inp;  S.init( var, inp );
    ok = ok && S.solve( var.data(), inp.data() ).converged;
    check( ok && std::fabs( S.val_functions()[0] - std::exp( -.25 * M_PI * M_PI ) ) < 1e-6,
           "  OCFESLV heat: u(1,1/2) == exp(-pi^2/4) (1e-6)" );
  }
  std::printf( "\n  AUTO_role: %d passed, %d failed -- %s\n", npass, nfail, nfail? "FAILURES": "ALL PASS" );
  return nfail? 1: 0;
}
