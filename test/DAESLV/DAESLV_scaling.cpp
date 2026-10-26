// DAESLV_scaling -- a few hundred equations, and DENSE against SPARSE (KLU).
// =====================================================================================================
// Every ODESLV and DAESLV driver so far is small -- a handful of states -- where a dense factorisation
// is never the wrong choice and a banded Jacobian never pays.  This one is the opposite: a
// reaction-diffusion system discretised by the method of lines, so the state count is a parameter and
// the Jacobian is TRIDIAGONAL by construction.
//
// The model is the Belousov-Zhabotinskii system of Peric, Villanueva & Chachuat (SIAM J. Sci. Comput.
// 39(6), 2017, section 5.4), in dimensionless form on -30 <= z <= 20:
//
//     du1/dt = d2u1/dz2 + u1 ( 1 - u1 - r u2 )
//     du2/dt = d2u2/dz2 - b u1 u2
//
// with no-flux boundaries and a moving front as the initial condition.  N interior points give 2N
// states for ODESLV.
//
// DAESLV gets the same physics in FLUX form, which is how a chemical engineer would write it and
// which is a genuine PDAE rather than an ODE in disguise:
//
//     du_i/dt = -( J_i,k+1 - J_i,k ) / h + R_i ,      0 = J_i,k + ( u_i,k - u_i,k-1 ) / h
//
// The fluxes at the N-1 interior faces are algebraic -- no dJ/dt anywhere -- so the system is index 1
// with 2N differential and 2(N-1) algebraic equations, about 4N in total.  (A first cut carried the
// two zero boundary fluxes as states pinned by J = 0; CRONOS's index analysis read those trivial rows
// as reduction triggers and materialised eight hidden levels, overdetermining the IC system.)  Both solvers therefore integrate the same problem,
// and their trajectories can be compared directly.
//
// What this measures: wall time for DENSE against SPARSE at several N, for each solver.  The dense
// factorisation is O(n^3) per step on a matrix that is almost all zeros; KLU should overtake it
// somewhere, and where it does is the number worth knowing.
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
#include "daeslv_idas.hpp"
#include "odeslv_cvodes.hpp"

using namespace mc;

static int DISP = 0;      // --display turns the per-phase setup timings on
static int n_pass = 0, n_fail = 0;
static void check_bool( std::string const& w, bool ok )
{
  ok? ++n_pass: ++n_fail;
  std::cout << "    " << std::left << std::setw( 50 ) << w << ( ok? "  PASS\n": "  FAIL\n" ) << std::flush;
}
static void check( std::string const& w, double got, double want, double tol )
{
  double const d = std::fabs( got - want );
  ( d <= tol )? ++n_pass: ++n_fail;
  std::cout << "    " << std::left << std::setw( 50 ) << w << std::right
            << " got=" << std::scientific << std::setprecision(8) << got << " want=" << want
            << ( d <= tol? "  PASS\n": "  FAIL\n" ) << std::flush;
}

double const ZL = -30., ZR = 20., RK = 10., BK = 1.25, TEND = 2.;

// ---------------------------------------------------------------- the ODE form, for ODESLV
static double run_ode( size_t const N, int const linsol, double& secs, double& setup_s, bool& ok )
{
  FFPartial OpP;  FFEval OpE;
  int const TI = FFDom::ALL - FFDom::LB;
  auto const INTER = FFModel::EqnRole::INTERIOR;
  auto const INIT  = FFModel::EqnRole::INITIAL;
  double const h = ( ZR - ZL ) / (double)N;

  if( DISP ) std::cout << "Starting ODE model\n";
  FFGraph G;  ODESLV_CVODES M;  M.set_dag( &G );
  FFVar t = G.add_var("t");
  std::vector<FFVar> u1( N ), u2( N );
  for( size_t k = 0; k < N; ++k ){
    u1[k] = G.add_var( "u1_" + std::to_string(k) );
    u2[k] = G.add_var( "u2_" + std::to_string(k) );
  }
  M.add_domain( t, FFDom( 0., TEND, 2, FFDom::LGR, 3 ) );
  M.set_evolution_domain( t );
  for( size_t k = 0; k < N; ++k ){
    M.add_state( u1[k], {t} );  M.update_ref( u1[k], 0.5 );
    M.add_state( u2[k], {t} );  M.update_ref( u2[k], 0.5 );
  }
  auto lap = [&]( std::vector<FFVar> const& u, size_t k ){      // no-flux at both ends
    FFVar const& um = ( k == 0 )?     u[0]   : u[k-1];
    FFVar const& up = ( k == N-1 )?   u[N-1] : u[k+1];
    return ( um - 2.*u[k] + up ) / ( h*h );
  };
  for( size_t k = 0; k < N; ++k ){
    M.add_equation( OpP( u1[k], t ) - lap( u1, k ) - u1[k]*( 1. - u1[k] - RK*u2[k] ), {t}, {TI}, INTER );
    M.add_equation( OpP( u2[k], t ) - lap( u2, k ) + BK*u1[k]*u2[k],                  {t}, {TI}, INTER );
    double const f1 = ( k <  N/2 )? 0.: ( k == N/2? 0.5: 1. );
    double const f2 = ( k <  N/2 )? 1.: ( k == N/2? 0.5: 0. );
    M.add_equation( u1[k] - f1, {t}, {FFDom::LB}, INIT );
    M.add_equation( u2[k] - f2, {t}, {FFDom::LB}, INIT );
  }
  FFVar avg( 0. );
  for( size_t k = 0; k < N; ++k ) avg += u1[k];
  M.add_output( OpE( avg/(double)N, t, TEND ) );                // the spatial average at T

  if( DISP ) std::cout << "Starting ODE setup\n";
  M.options.SOLVE.RTOL = 1e-8;  M.options.SOLVE.ATOL = 1e-10;
  M.options.SOLVE.LINSOL = (BASE_CVODES::Options::LINEAR_SOLVER)linsol;
  M.options.SOLVE.VERBOSE = 0;  M.options.DISPLAY_LEVEL = DISP;
  auto const tsu = std::chrono::steady_clock::now();
  bool const su = M.setup();
  setup_s = std::chrono::duration<double>( std::chrono::steady_clock::now() - tsu ).count();
  if( !su ){ std::cout << "      setup failed: " << M.extract_error() << "\n"; ok = false; return 0.; }

  if( DISP ) std::cout << "Starting ODE solve\n";
  std::vector<double> P( M.np(), 0. );
  auto const t0 = std::chrono::steady_clock::now();
  ok = ( M.solve( P ) == ODESLV_CVODES::NORMAL );
  secs = std::chrono::duration<double>( std::chrono::steady_clock::now() - t0 ).count();
  return ok? M.val_function()[0]: 0.;
}

// ---------------------------------------------------------------- the PDAE form, for DAESLV
static double run_dae( size_t const N, int const linsol, double& secs, double& setup_s, bool& ok, size_t& nx )
{
  FFPartial OpP;  FFEval OpE;
  int const TI = FFDom::ALL - FFDom::LB;
  auto const INTER = FFModel::EqnRole::INTERIOR;
  auto const INIT  = FFModel::EqnRole::INITIAL;
  double const h = ( ZR - ZL ) / (double)N;

  if( DISP ) std::cout << "Starting DAE model\n";
  FFGraph G;  DAESLV_IDAS M;  M.set_dag( &G );
  FFVar t = G.add_var("t");
  std::vector<FFVar> u1( N ), u2( N ), J1( N+1 ), J2( N+1 );   // J[0], J[N] unused: no-flux boundaries
  for( size_t k = 0; k < N; ++k ){
    u1[k] = G.add_var( "u1_" + std::to_string(k) );
    u2[k] = G.add_var( "u2_" + std::to_string(k) );
  }
  for( size_t k = 1; k < N; ++k ){           // interior faces only
    J1[k] = G.add_var( "J1_" + std::to_string(k) );
    J2[k] = G.add_var( "J2_" + std::to_string(k) );
  }
  M.add_domain( t, FFDom( 0., TEND, 2, FFDom::LGR, 3 ) );
  M.set_evolution_domain( t );
  for( size_t k = 0; k < N; ++k ){
    M.add_state( u1[k], {t} );  M.update_ref( u1[k], 0.5 );
    M.add_state( u2[k], {t} );  M.update_ref( u2[k], 0.5 );
  }
  for( size_t k = 1; k < N; ++k ){           // the ALGEBRAIC states: no derivative of their own
    M.add_state( J1[k], {t} );  M.update_ref( J1[k], 0. );
    M.add_state( J2[k], {t} );  M.update_ref( J2[k], 0. );
  }
  auto div = [&]( std::vector<FFVar> const& J, size_t k ) -> FFVar {   // (J_k+1 - J_k)/h, zero flux at the ends
    if( N == 1 ) return FFVar( 0. );
    if( k == 0 )   return   J[1]/h;
    if( k == N-1 ) return - J[N-1]/h;
    return ( J[k+1] - J[k] )/h;
  };
  for( size_t k = 0; k < N; ++k ){
    M.add_equation( OpP( u1[k], t ) + div( J1, k ) - u1[k]*( 1. - u1[k] - RK*u2[k] ), {t}, {TI}, INTER );
    M.add_equation( OpP( u2[k], t ) + div( J2, k ) + BK*u1[k]*u2[k],                  {t}, {TI}, INTER );
    double const f1 = ( k <  N/2 )? 0.: ( k == N/2? 0.5: 1. );
    double const f2 = ( k <  N/2 )? 1.: ( k == N/2? 0.5: 0. );
    M.add_equation( u1[k] - f1, {t}, {FFDom::LB}, INIT );
    M.add_equation( u2[k] - f2, {t}, {FFDom::LB}, INIT );
  }
  // Fick's law at the interior faces -- algebraic, everywhere
  for( size_t k = 1; k < N; ++k ){
    M.add_equation( J1[k] + ( u1[k] - u1[k-1] )/h, {t}, {FFDom::ALL}, INTER );
    M.add_equation( J2[k] + ( u2[k] - u2[k-1] )/h, {t}, {FFDom::ALL}, INTER );
  }

  FFVar avg( 0. );
  for( size_t k = 0; k < N; ++k ) avg += u1[k];
  M.add_output( OpE( avg/(double)N, t, TEND ) );

  if( DISP ) std::cout << "Starting DAE setup\n";
  M.options.SOLVE.RTOL = 1e-8;  M.options.SOLVE.ATOL = 1e-10;
  M.options.SOLVE.LINSOL = (BASE_IDAS::Options::LINEAR_SOLVER)linsol;
  M.options.SOLVE.VERBOSE = 0;  M.options.DISPLAY_LEVEL = DISP;
  auto const tsu = std::chrono::steady_clock::now();
  bool const su = M.setup();
  setup_s = std::chrono::duration<double>( std::chrono::steady_clock::now() - tsu ).count();
  if( !su ){ std::cout << "      DAE setup failed at N=" << N << ": " << M.extract_error() << "\n" << std::flush; ok = false; return 0.; }
  nx = M.nx();

  if( DISP ) std::cout << "Starting DAE solve\n";
  std::vector<double> P( M.np(), 0. );
  auto const t0 = std::chrono::steady_clock::now();
  auto const st = M.solve( P );
  ok = ( st == DAESLV_IDAS::NORMAL );
  secs = std::chrono::duration<double>( std::chrono::steady_clock::now() - t0 ).count();
  if( !ok ) std::cout << "      DAE solve failed at N=" << N << ": " << M.extract_error() << "\n" << std::flush;
  return ok? M.val_function()[0]: 0.;
}

int main( int argc, char* argv[] )
{
  for( int i = 1; i < argc; ++i )
    if( std::string( argv[i] ) == "--display" ) DISP = 2;
  std::cout << "\n--- Belousov-Zhabotinskii by the method of lines: DENSE against SPARSE (KLU)\n"
            << "    ODESLV integrates the divided form (2N states); DAESLV the FLUX form, a genuine\n"
            << "    index-1 PDAE (2N differential + 2(N-1) algebraic).  Same physics, same answer.\n\n";
#if !defined( CRONOS__WITH_KLU )
  std::cout << "    *** built without KLU: SPARSE is unavailable, so only DENSE is timed ***\n";
  bool const klu = false;
#else
  bool const klu = true;
#endif

  // Two lines per mesh: the SOLVE times, then the SETUP times under them.  Setup used to dominate
  // -- at 2400 states it was 5.7 s against a 2 s solve -- so both belong in the table until that
  // stays true (2026-10-25).
  std::cout << "    " << std::left << std::setw(7) << "N" << std::setw(8) << "states" << std::setw(8) << ""
            << std::setw(13) << "ODE dense" << std::setw(13) << "ODE sparse"
            << std::setw(13) << "DAE dense" << std::setw(13) << "DAE sparse" << "\n" << std::flush;

  double ref_ode = 0., ref_dae = 0.;
  // Above NDENSE only SPARSE is timed: a dense Jacobian costs n^2 doubles per factorisation, which
  // is the thing KLU exists to avoid, and timing it at every mesh is not the point.
  //
  // MEMORY.  A dense Jacobian is n^2 doubles PER factorisation: 9600 states at N=1600 is ~740 MB,
  // so the dense columns may fail there on a modest machine -- which is the point KLU exists for and
  // a result worth recording rather than a problem.  NDENSE caps the mesh at which dense is timed;
  // lower it if the dense runs are what fail.
  //
  // Setup was the wall until 2026-10-25 -- model construction, not integration: 5.7 s against a 2 s
  // solve at 2400 states, and a kill at N=400 on a 4 GB box.  The setup row of the table is here to
  // keep that honest.  Run with --display for the per-phase breakdown.
  size_t const NDENSE = 1600;
  for( size_t N : { 20, 50, 100, 200, 400, 800 } ){//, 1600 } ){
    double td_o = 0., ts_o = 0., td_d = 0., ts_d = 0.;
    double su_td_o = 0., su_ts_o = 0., su_td_d = 0., su_ts_d = 0.;
    bool ok1 = false, ok2 = true, ok3 = false, ok4 = true;
    size_t nxd = 0;
    bool const dense = ( N <= NDENSE );
    double vo = 0., vd = 0.;
    if( dense ){
      vo = run_ode( N, (int)BASE_CVODES::Options::DENSE, td_o, su_td_o, ok1 );
      vd = run_dae( N, (int)BASE_IDAS::Options::DENSE,   td_d, su_td_d, ok3, nxd );
    }
    double vos = vo, vds = vd;
    if( klu ){
      vos = run_ode( N, (int)BASE_CVODES::Options::SPARSE, ts_o, su_ts_o, ok2 );
      vds = run_dae( N, (int)BASE_IDAS::Options::SPARSE,   ts_d, su_ts_d, ok4, nxd );
      if( !dense ){ ok1 = ok2;  ok3 = ok4;  vo = vos;  vd = vds; }
    }
    auto cell = [&]( bool shown, double v ){
      if( shown ) std::cout << std::setw(13) << v;
      else        std::cout << std::setw(13) << "-";
    };
    std::cout << "    " << std::left << std::setw(7) << N << std::setw(8) << ( 2*N + nxd )
              << std::setw(8) << "solve" << std::right << std::fixed << std::setprecision(3);
    cell( dense, td_o );  cell( klu, ts_o );  cell( dense, td_d );  cell( klu, ts_d );
    std::cout << "\n    " << std::left << std::setw(7) << "" << std::setw(8) << ""
              << std::setw(8) << "setup" << std::right;
    cell( dense, su_td_o );  cell( klu, su_ts_o );  cell( dense, su_td_d );  cell( klu, su_ts_d );
    std::cout << "\n" << std::flush;
    check_bool( "N=" + std::to_string(N) + ": every solve succeeded", ok1 && ok2 && ok3 && ok4 );
    if( ok1 && ok3 ){
      // the two formulations describe the same physics: same average at T
      check( "N=" + std::to_string(N) + ": ODE form == PDAE form", vd, vo, 1e-6 );
      if( klu && dense ){
        check( "N=" + std::to_string(N) + ": DENSE == SPARSE (ODE)", vos, vo, 1e-8 );
        check( "N=" + std::to_string(N) + ": DENSE == SPARSE (DAE)", vds, vd, 1e-8 );
      }
    }
    if( N == 20 )  ref_ode = vo;
    if( N == 800 ) ref_dae = vo;
  }
  // it IS a discretisation: the answer moves with the mesh, and by less as the mesh refines
  check_bool( "the answer moves with the mesh", std::fabs( ref_ode - ref_dae ) > 1e-4 );

  std::cout << "\n==========  DAESLV_scaling: " << n_pass << " passed, " << n_fail << " failed -- "
            << ( n_fail? "FAILURES": "ALL PASS" ) << "  ==========\n\n";
  return n_fail? 1: 0;
}
