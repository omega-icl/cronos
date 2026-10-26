// DAESLV_jypprobe -- the adjoint on models WITH an algebraic row.
// =====================================================================================================
// The terminal condition assembles Jyp = dF/dy' into a dense _nx x _nx matrix and factorises it.  An
// algebraic row contributes a ZERO ROW to Jyp, so that matrix should be structurally singular for ANY
// DAE with one -- index 1 included.  Yet every index-1 adjoint case in the corpus passes.  Either the
// factorisation is succeeding for a reason I have not understood, or those models have no algebraic
// row by the time index reduction is done with them.  This driver asks which, on three models:
//   A  pure ODE, no algebraic row at all          -- Jyp should be full rank
//   B  index 1, one algebraic state               -- Jyp should have one zero row
//   C  index 2, the tutorial's vented reactor     -- Jyp should have one zero row, and ASA fails here
//
// C is the DISCRIMINATING case, and A and B are not.  The multiplier of an algebraic equation is
// determined by the adjoint's own algebraic equation, so it vanishes whenever the algebraic variable
// does not feed back into the differential ones.  In B, z = k x is computed FROM x and never used, so
// lambda_z = 0 and consistent initialisation has nothing to do.  In C the vent rate q appears in both
// differential equations, and IDACalcICB moves its multiplier from 0 to 0.2495 (measured).  So C is
// the case that actually exercises the algebraic part of the terminal condition.
#include <cmath>
#include <iomanip>
#include <iostream>
#include "daeslv_idas.hpp"

using namespace mc;

static int n_pass = 0, n_fail = 0;
static void note( std::string const& w, std::string const& v )
{ std::cout << "    " << std::left << std::setw( 46 ) << w << v << "\n" << std::flush; }
static void check_bool( std::string const& w, bool ok )
{ ok? ++n_pass: ++n_fail;
  std::cout << "    " << std::left << std::setw( 46 ) << w << ( ok? "  PASS\n": "  FAIL\n" ) << std::flush; }

// how many rows of dF/dy' are entirely zero, and does solve_asens succeed?
static void probe( std::string const& tag, DAESLV_IDAS& M, std::vector<double> const& P, bool expect_zero_row )
{
  // the adjoint must not merely RUN: it must agree with the forward route, which is the thing the
  // original dense terminal solve could not do on any model with an algebraic row.
  std::cout << "\n--- " << tag << "\n" << std::flush;
  note( "nx (augmented, minted levels included)", std::to_string( M.nx() ) );
  M.solve( P );
  auto const st = M.solve_asens( P );
  note( "solve_asens", st == DAESLV_IDAS::NORMAL? "NORMAL": ( "FAILED: " + M.extract_error() ) );
  if( expect_zero_row )
    note( "expected", "an algebraic row -> a zero row in dF/dy'" );
  else
    note( "expected", "no algebraic row -> dF/dy' full rank" );
  check_bool( "adjoint " + std::string( expect_zero_row? "with": "without" ) + " an algebraic row runs",
              ( st == DAESLV_IDAS::NORMAL ) );
  if( st != DAESLV_IDAS::NORMAL ) return;
  auto const Ga = M.val_function_gradient();
  M.solve_fsens( P );
  auto const& Gf = M.val_function_gradient();
  double worst = 0.;
  for( size_t j = 0; j < Gf.size(); ++j )
    for( size_t r = 0; r < Gf[j].size() && r < Ga[j].size(); ++r )
      worst = std::max( worst, std::fabs( Gf[j][r] - Ga[j][r] ) );
  std::cout << "    " << std::left << std::setw( 46 ) << "max |adjoint - forward|"
            << std::scientific << std::setprecision(2) << worst << "\n" << std::flush;
  check_bool( "  and AGREES with solve_fsens", worst < 1e-6 );
}

int main()
{
  FFPartial OpP;  FFEval OpE;
  int const TI = FFDom::ALL - FFDom::LB;
  auto io = FFModel::EqnOptions( FFModel::EqnRole::INTERIOR, 0 );
  auto ii = FFModel::EqnOptions( FFModel::EqnRole::INITIAL,  0 );

  // ---------------------------------------------------------------- A: a pure ODE
  { FFGraph G;  DAESLV_IDAS M;  M.set_dag( &G );
    FFVar t = G.add_var("t"), x = G.add_var("x(t)"), k = G.add_var("k");
    M.add_domain( t, FFDom( 0., 1., 2, FFDom::LGR, 3 ) );
    M.add_state( x, {t} );  M.add_input( k );  M.set_evolution_domain( t );  M.update_ref( x, 1. );
    M.add_equation( OpP( x, t ) + k*x, {t}, {TI}, io );
    M.add_equation( x - 1.,            {t}, {FFDom::LB}, ii );
    M.add_output( OpE( x, t, 1. ) );
    M.options.SOLVE.RTOL = M.options.SOLVE.ATOL = 1e-10;  M.options.DISPLAY_LEVEL = 0;
    if( !M.setup() ){ std::cout << "    setup A failed: " << M.extract_error() << "\n"; ++n_fail; }
    else{ std::vector<double> P( M.np(), 0. ); P[ M.parameter_index(k)[0] ] = 0.8;
          probe( "A: pure ODE, no algebraic row", M, P, false ); } }

  // ---------------------------------------------------------------- B: index 1, one algebraic state
  // z is fixed by an algebraic row that MENTIONS it, so one differentiation is not needed: index 1.
  { FFGraph G;  DAESLV_IDAS M;  M.set_dag( &G );
    FFVar t = G.add_var("t"), x = G.add_var("x(t)"), z = G.add_var("z(t)"), k = G.add_var("k");
    M.add_domain( t, FFDom( 0., 1., 2, FFDom::LGR, 3 ) );
    M.add_state( x, {t} );  M.add_state( z, {t} );  M.add_input( k );
    M.set_evolution_domain( t );  M.update_ref( x, 1. );  M.update_ref( z, 1. );
    M.add_equation( OpP( x, t ) + k*x, {t}, {TI}, io );
    M.add_equation( z - k*x,           {t}, {FFDom::ALL}, io );     // algebraic, solvable for z
    M.add_equation( x - 1.,            {t}, {FFDom::LB}, ii );
    M.add_output( OpE( x, t, 1. ) );
    M.options.SOLVE.RTOL = M.options.SOLVE.ATOL = 1e-10;  M.options.DISPLAY_LEVEL = 0;
    if( !M.setup() ){ std::cout << "    setup B failed: " << M.extract_error() << "\n"; ++n_fail; }
    else{ std::vector<double> P( M.np(), 0. ); P[ M.parameter_index(k)[0] ] = 0.8;
          probe( "B: index 1, one algebraic state", M, P, true ); } }

  // ---------------------------------------------------------------- C: index 2, the vented reactor
  { FFGraph G;  DAESLV_IDAS M;  M.set_dag( &G );
    FFVar t = G.add_var("t"), nA = G.add_var("nA(t)"), nB = G.add_var("nB(t)"), q = G.add_var("q(t)"),
          k = G.add_var("k");
    double const NT = 1.;
    M.add_domain( t, FFDom( 0., 1., 2, FFDom::LGR, 3 ) );
    M.add_state( nA, {t} );  M.add_state( nB, {t} );  M.add_state( q, {t} );  M.add_input( k );
    M.set_evolution_domain( t );
    M.update_ref( nA, 0.5 );  M.update_ref( nB, 0.5 );  M.update_ref( q, 0.1 );
    M.add_equation( OpP( nA, t ) + k*nA + q*nA/NT,    {t}, {TI}, io );
    M.add_equation( OpP( nB, t ) - 2.*k*nA + q*nB/NT, {t}, {TI}, io );
    M.add_equation( nA + nB - NT,                     {t}, {FFDom::ALL}, io );   // no q: index 2
    M.add_equation( nA - 0.8,                         {t}, {FFDom::LB}, ii );
    M.add_output( OpE( nA, t, 1. ) );
    M.options.SOLVE.RTOL = M.options.SOLVE.ATOL = 1e-10;  M.options.DISPLAY_LEVEL = 0;
    if( !M.setup() ){ std::cout << "    setup C failed: " << M.extract_error() << "\n"; ++n_fail; }
    else{ std::vector<double> P( M.np(), 0. ); P[ M.parameter_index(k)[0] ] = 0.8;
          probe( "C: index 2, the vented reactor", M, P, true ); } }

  // ---------------------------------------------------------------- D: index 1 with a declared JUMP
  // The same algebraic state as B, but z FEEDS BACK ( x' = -k x - z ) so its multiplier is nonzero,
  // and a transition at t = 0.5 adds 0.2 to x.  Every adjoint-across-a-jump case in the corpus is an
  // ODE; this is the first on a DAE, and it exercises the lambda^- recovery after the stop system.
  { FFGraph G;  DAESLV_IDAS M;  M.set_dag( &G );
    FFVar t = G.add_var("t"), x = G.add_var("x(t)"), z = G.add_var("z(t)"), k = G.add_var("k");
    M.add_domain( t, FFDom( 0., 1., 2, FFDom::LGR, 3 ) );
    M.add_state( x, {t} );  M.add_state( z, {t} );  M.add_input( k );
    M.set_evolution_domain( t );  M.update_ref( x, 1. );  M.update_ref( z, 1. );
    M.add_equation( OpP( x, t ) + k*x + z, {t}, {TI}, io );
    M.add_equation( z - k*x,               {t}, {FFDom::ALL}, io );
    M.add_equation( x - 1.,                {t}, {FFDom::LB}, ii );
    M.add_transition( x + 0.2, x, t, 0.5 );              // left(tau^-) + 0.2 = right(tau^+)
    M.add_output( OpE( x, t, 1. ) );
    M.options.SOLVE.RTOL = M.options.SOLVE.ATOL = 1e-10;  M.options.DISPLAY_LEVEL = 0;
    if( !M.setup() ){ std::cout << "    setup D failed: " << M.extract_error() << "\n"; ++n_fail; }
    else{ std::vector<double> P( M.np(), 0. ); P[ M.parameter_index(k)[0] ] = 0.8;
          probe( "D: index 1, z feeds back, a declared jump at 0.5", M, P, true ); } }

  std::cout << "\n==========  DAESLV_jypprobe: " << n_pass << " passed, " << n_fail << " failed  ==========\n\n";
  return n_fail? 1: 0;
}
