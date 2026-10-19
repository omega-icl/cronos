// OCFE_TRANSderiv.cpp -- a transition whose law reads the state TIME DERIVATIVE at tau.
//
//   dx/dt = -k x on [0, TF], x(0) = x0, with a jump at tau:   x(tau+) = -a dx/dt(tau-)
//
// dx/dt(tau-) = -k x0 e^{-k tau}, so x(tau+) = a k x0 e^{-k tau} and x(TF) = x(tau+) e^{-k(TF-tau)}.
// Choosing a != 1/k makes the jump a genuine discontinuity, so a law that read the RIGHT-limit
// derivative, or no derivative at all, cannot reach the same answer.
//
// The same model is run through DAESLV_IDAS, which applies the jump as a ROW of its stop system solved
// AT tau -- so the one-sided derivative is an unknown already there.  OCFESLV instead LIFTS a transition
// into an auxiliary and appends a link row, which a derivative operator makes non-pointwise: its
// interface audit cannot square the system and, before _transition_derivatives_supported() gated it,
// it ABORTED rather than refused.  So the capability is DAESLV's, and what this driver pins down is
// that OCFESLV refuses it cleanly instead of returning a wrong answer.
//
// tau is an ELEMENT BOUNDARY of the mesh, which both solvers require.
#include <cmath>
#include <cstdio>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include "daeslv_idas.hpp"
#include "ocfeslv.hpp"
using namespace mc;

static int npass = 0, nfail = 0;
static void check( bool c, std::string const& w, double v = -1. )
{ std::printf( "  %s  %-62s", c? "PASS": "FAIL", w.c_str() );
  if( v >= 0. ) std::printf( " (%.2e)", v );
  std::printf( "\n" );  c? ++npass: ++nfail; }

double const TAU = 0.4, TF = 1.0, K = 1.3, A = 0.8, X0 = 1.1;

//! x(tau-), x(tau+) and x(TF) in closed form
static std::vector<double> closed( double k, double a, double x0 )
{
  double const xm = x0*std::exp( -k*TAU );
  double const dm = -k*xm;                       // dx/dt(tau-)
  double const xp = -a*dm;                       // the jump law
  return { xm, xp, xp*std::exp( -k*( TF-TAU ) ) };
}

//! The model through OCFESLV, which LIFTS the transition
static bool run_oc( bool march, std::vector<double>& F )
{
  FFGraph G; OCFESLV I( &G ); FFPartial OpP; FFEval OpE;
  FFVar t = G.add_var("t"), x = G.add_var("x(t)"),
        pk = G.add_var("k"), pa = G.add_var("a"), p0 = G.add_var("x0");
  I.add_domain( t, FFDom( std::vector<double>{ 0., TAU, 1. }, FFDom::LGR, 8 ) );
  I.set_evolution_domain( t );
  I.add_state( x, {t} );  for( auto const& v : { pk, pa, p0 } ) I.add_input( v, {} );
  I.update_ref( x, 1. );
  int const T_INT = FFDom::ALL - FFDom::LB;
  I.add_equation( OpP( x, t ) + pk*x, {t}, {T_INT}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  I.add_equation( x - p0, {t}, {FFDom::LB}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INITIAL, 0 ) );
  I.add_transition( -pa*OpP( x, t ), x, t, TAU );             // reads dx/dt at tau^-
  I.add_output( OpE( x, t, TAU, FFDom::MINUS ) );
  I.add_output( OpE( x, t, TAU, FFDom::PLUS ) );
  I.add_output( OpE( x, t, TF ) );
  I.options.SOLVE.MARCHING = march;  I.options.SOLVE.RES_TOL = 1e-12;  I.options.DISPLAY_LEVEL = 0;
  std::ostringstream os;  auto* old = std::cerr.rdbuf( os.rdbuf() );
  std::vector<double> xv, inp;
  bool ok = I.setup() && I.init( xv, inp, nullptr );
  if( ok ){ I.set_input_values( pk, { K }, inp.data() );
            I.set_input_values( pa, { A }, inp.data() );
            I.set_input_values( p0, { X0 }, inp.data() ); }
  std::vector<double> x1( xv );
  if( ok ) ok = I.solve( x1.data(), inp.data(), nullptr ).converged;
  if( ok ) F = I.val_functions();
  std::cerr.rdbuf( old );
  if( !ok ) std::cout << "    [ocfeslv] " << os.str().substr( 0, 200 ) << "\n";
  return ok && F.size() == 3;
}

//! The same model through DAESLV, which applies the jump as a ROW of the stop system
static bool run_dae( double& xm, double& xf )
{
  FFGraph G; DAESLV_IDAS M; M.set_dag( &G );  FFPartial OpP;
  FFVar t = G.add_var("t"), x = G.add_var("x(t)");
  M.add_domain( t, FFDom( std::vector<double>{ 0., TAU, 1. }, FFDom::LGR, 8 ) );
  M.add_state( x, {t} );  M.set_evolution_domain( t );  M.update_ref( x, X0 );
  int const T_INT = FFDom::ALL - FFDom::LB;
  M.add_equation( OpP( x, t ) + K*x, {t}, {T_INT}, FFModel::EqnOptions( FFModel::EqnRole::INTERIOR, 0 ) );
  M.add_equation( x - X0, {t}, {FFDom::LB}, FFModel::EqnOptions( FFModel::EqnRole::INITIAL, 0 ) );
  M.add_transition( -A*OpP( x, t ), x, t, TAU );
  M.options.SOLVE.RTOL = M.options.SOLVE.ATOL = 1e-12;  M.options.DISPLAY_LEVEL = 0;
  if( !M.setup() ){ std::cout << "    [daeslv] " << M.extract_error() << "\n"; return false; }
  if( M.solve( std::vector<double>{} ) != DAESLV_IDAS::NORMAL ) return false;
  auto const& X = M.val_state();
  if( X.size() < 3 ) return false;
  xm = X[1][0];  xf = X.back()[0];
  return true;
}

int main()
{
  std::printf( "\n  OCFE_TRANSderiv: a jump law reading dx/dt(tau-)\n\n" );
  auto const F0 = closed( K, A, X0 );
  std::printf( "  closed form: x(tau-)=%.10e  x(tau+)=%.10e  x(TF)=%.10e\n\n", F0[0], F0[1], F0[2] );

  { std::vector<double> F;  bool threw = false, ok = false;
    try{ ok = run_oc( false, F ); } catch( ... ){ threw = true; }
    check( !ok && !threw, "OCFESLV  : REFUSES a derivative-valued jump, cleanly" ); }

  double xm = 0., xf = 0.;
  if( !run_dae( xm, xf ) ) check( false, "DAESLV   : setup and solves" );
  else{
    check( true, "DAESLV   : setup and solves" );
    check( std::fabs( xm-F0[0] ) < 1e-8, "DAESLV   : x(tau-) == closed form", std::fabs( xm-F0[0] ) );
    check( std::fabs( xf-F0[2] ) < 1e-8, "DAESLV   : x(TF)   == closed form", std::fabs( xf-F0[2] ) );
    check( std::fabs( F0[1]-F0[0] ) > 1e-2, "DAESLV   : the jump is a real discontinuity",
           std::fabs( F0[1]-F0[0] ) );
  }

  std::printf( "\n  OCFE_TRANSderiv: %d passed, %d failed -- %s\n\n", npass, nfail,
               nfail? "FAILURES": "ALL PASS" );
  return nfail? 1: 0;
}
