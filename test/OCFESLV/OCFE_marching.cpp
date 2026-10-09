// Copyright (C) Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the EPL-2.0 with GPL-2.0-or-later as a Secondary License; see the LICENSE file.
// SPDX-License-Identifier: EPL-2.0 OR GPL-2.0-or-later
//
// OCFE_marching.cpp
// -----------------
// Regression test for the Stage-0 element-by-element marching primitives on an
// EVOL_HYPERBOLIC problem (a scalar linear advection, the same class as the
// blocks in OCFE_PDE8).
//
//   d_t u + a d_z u = 0     u_exact(t,z) = u0 sin( k ( z - a t ) )   (traveling wave, f=0)
//   IC (t=LB):  u(0,z) = u0 sin( k z )
//   BC (z=LB):  u(t,0) = u0 sin( -k a t )     (inflow, a>0)
//
// It exercises the marching read primitives added to OCFESLV:
//   * n_evolution_elem()           -- element count along the evolution domain
//   * evolution_slab_selfcheck()   -- the slab column-map partition check
//   * terminal_value()             -- interpolated state value at the evolution UB
//   * terminal_profile()           -- terminal profile over all spatial nodes
//
// Crucially, it runs the SAME physical problem with two DOF layouts:
//   (a) evolution domain declared FIRST  -> t has the lower id -> t is the innermost
//       stride, so evolution slabs are STRIDED;
//   (b) evolution domain declared LAST   -> t is the outermost stride, so evolution
//       slabs are CONTIGUOUS.
// The marching primitives must be layout-agnostic: the terminal profile must be
// IDENTICAL (to solver tolerance) across the two orderings.  This is the real
// check that nothing assumes the incidental "evolution declared first" convention.

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include <armadillo>
#ifndef OCFE_OCFESLV_HEADER
#define OCFE_OCFESLV_HEADER "ocfeslv.hpp"
#endif
#include OCFE_OCFESLV_HEADER

using namespace mc;

namespace {

struct Par {
  double u0 = 0.7;
  double k  = 6.283185307179586;   // 2*pi
  double a  = 0.9;                  // advection speed (>0 -> inflow at z=LB)
  double tf = 0.4;                  // evolution (t) upper bound
  double zf = 1.0;                  // spatial   (z) upper bound
};

constexpr size_t NEL_T = 3, NT = 6;   // evolution elements / nodes-per-element
constexpr size_t NEL_Z = 3, NZ = 8;   // spatial   elements / nodes-per-element

static double u_exact( double t, double z, Par const& p )
{ return p.u0 * std::sin( p.k * ( z - p.a * t ) ); }

int g_pass = 0, g_fail = 0;

static bool report( std::string const& label, bool ok, double err = -1., double tol = -1. )
{
  std::cout << "  " << std::left << std::setw(50) << label;
  if( err >= 0. )
    std::cout << " err=" << std::scientific << std::setprecision(3) << err
              << " tol=" << tol;
  std::cout << "  " << ( ok ? "PASS" : "FAIL" ) << "\n";
  ( ok ? g_pass : g_fail )++;
  return ok;
}
static bool check_true ( std::string const& l, bool ok ) { return report( l, ok ); }
static bool check_close( std::string const& l, double got, double exp, double tol )
{ double e = std::abs( got - exp ); return report( l, e <= tol, e, tol ); }

// Flat physical node coordinates of a domain, element-major / node-inner.
static std::vector<double> physical_nodes( FFDom const& dom )
{
  std::vector<double> xs;
  for( size_t ie = 0; ie < dom.n_elem; ++ie ){
    auto xe = dom.lgnodes( dom.elem_lo( ie ), dom.elem_up( ie ) );
    xs.insert( xs.end(), xe.begin(), xe.end() );
  }
  return xs;
}

// Build + solve the advection problem for one DOF layout, then exercise the
// marching read primitives.  Returns the terminal profile of u (empty on failure)
// so the caller can cross-check the two layouts against each other.
static std::vector<double> run( bool evolution_first, Par const& p )
{
  std::string const tag = evolution_first
    ? "(t,z): evolution declared FIRST  (t innermost -> strided slabs)"
    : "(z,t): evolution declared LAST   (t outermost -> contiguous slabs)";
  std::cout << "\n------------------------------------------------------------\n"
            << "  " << tag << "\n"
            << "------------------------------------------------------------\n";

  // Declaration order sets the id order, hence the stride order (add_state uses a
  // set sorted by id, so the initializer-list order is irrelevant -- only the
  // add_var order matters).
  FFGraph DAG;
  FFVar t, z;
  if( evolution_first ){ t = DAG.add_var( "t" ); z = DAG.add_var( "z" ); }
  else                 { z = DAG.add_var( "z" ); t = DAG.add_var( "t" ); }
  FFVar u    = DAG.add_var( "u(t,z)" );
  FFVar u_ic = DAG.add_var( "u_ic(z)" );   // IC posed as a distributed spatial input
  FFPartial OpP;

  FFDom zdom( 0., p.zf, NEL_Z, FFDom::LGL, NZ );

  OCFESLV oc( &DAG );
  oc.add_domain( t, FFDom( 0., p.tf, NEL_T, FFDom::LGR, NT ) );
  oc.add_domain( z, zdom );
  oc.add_state ( u, { t, z } );
  oc.update_ref( u, [&]( OCFESLV::t_Coord const& cd ){ return u_exact( cd.at(t), cd.at(z), p ); } );
  oc.add_input ( u_ic, { z }, FFDom::LGL, NZ );                           // explicitly the state z-grid
  oc.update_ref( u_ic, [&]( OCFESLV::t_Coord const& cd ){ return u_exact( 0., cd.at(z), p ); } );

  // d_t u + a d_z u = 0
  FFVar PDE = OpP( u, t ) + p.a * OpP( u, z );
  // IC at t=LB posed as an input:  u(0,z) = u_ic(z)   (u_ic is re-seeded per march step)
  FFVar INI = u - u_ic;
  // BC at z=LB (inflow):  u(t,0) = u0 sin(-k a t)
  FFVar BC  = u - p.u0 * sin( p.k * ( 0.0 - p.a * t ) );

  oc.add_equation( PDE, { t, z }, { FFDom::ALL - FFDom::LB, FFDom::ALL - FFDom::LB },
                   OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( INI, { t, z }, { FFDom::LB,              FFDom::ALL              },
                   OCFESLV::EqnOptions( OCFESLV::EqnRole::INITIAL, 0 ) );
  oc.add_equation( BC,  { t, z }, { FFDom::ALL - FFDom::LB, FFDom::LB               },
                   OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );

  oc.options.INTERFACE.TYPE = OCFESLV::Options::IC_AUTO;
  oc.options.CLASSIFY.MODE       = OCFESLV::Options::CLASS_AUTO;
  oc.options.DISPLAY_LEVEL  = 1;   // surface [probe] / [slabmap] / [evolution]
  oc.options.SOLVE.MARCHING = false;   // monolithic reference + read primitives on the FULL (uncollapsed) domain
  oc.set_evolution_domain( t );

  if( !oc.setup() ){
    std::cerr << "  setup() FAILED\n"; ++g_fail; return {};
  }

  check_true( "classified EVOL_HYPERBOLIC",
              oc.pde_type().type == OCFESLV::EVOL_HYPERBOLIC );

  // Monolithic reference solve.  Input values (u_ic) flow through inp[].
  std::vector<double> xv, inp;
  if( !oc.init( xv, inp ) ){ std::cerr << "  init() FAILED\n"; ++g_fail; return {}; }
  double const* ip = inp.empty() ? nullptr : inp.data();
  OCFESLV::SolveReport const srep = oc.solve( xv.data(), ip );
  check_true( "monolithic solve converged", srep.converged );

  // ---- marching read primitives ----
  check_true( "n_evolution_elem() == NEL_T", oc.n_evolution_elem() == NEL_T );
  check_true( "evolution_slab_selfcheck() partitions", oc.evolution_slab_selfcheck() );

  // terminal_value at a few interior spatial points, vs the exact terminal.
  for( double frac : { 0.15, 0.50, 0.85 } ){
    double const zt = frac * p.zf;
    OCFESLV::t_Coord pt; pt[z] = zt;
    double const got = oc.terminal_value( u, pt, xv.data(), ip );
    std::ostringstream lab; lab << "terminal_value(u, z=" << std::fixed
                                << std::setprecision(2) << zt << ") vs exact(UB)";
    check_close( lab.str(), got, u_exact( p.tf, zt, p ), 1e-3 );
  }

  // terminal_profile: size == spatial DOF count, values vs exact at UB.
  // Take the z nodes from the PROCESSED domain inside oc (a locally-constructed
  // FFDom is not node-populated until add_domain copies + processes it).
  auto const znodes = physical_nodes( oc.var_domain().at( z ) );
  auto prof = oc.terminal_profile( xv.data(), ip );
  auto itp = prof.find( u );
  bool const have = ( itp != prof.end() );
  check_true( "terminal_profile has state u", have );

  std::vector<double> pu = have ? itp->second : std::vector<double>{};
  std::cout << "    (terminal_profile[u].size=" << pu.size()
            << ", NEL_Z*NZ=" << ( NEL_Z * NZ )
            << ", physical_nodes(z)=" << znodes.size() << ")\n";
  check_true( "terminal_profile[u].size == NEL_Z*NZ",
              pu.size() == znodes.size() );

  if( pu.size() == znodes.size() ){
    double emax = 0.;
    for( size_t i = 0; i < pu.size(); ++i )
      emax = std::max( emax, std::abs( pu[i] - u_exact( p.tf, znodes[i], p ) ) );
    check_close( "terminal_profile[u] vs exact at UB", emax, 0., 1e-3 );
  }

  // ---- IC re-seed: the transfer that advances a march one element ----
  // Element k's terminal becomes element k+1's IC.  transfer_terminal() samples the
  // state's terminal at the IC input's OWN nodes and seeds it, so the re-seeded IC is
  // exact for how eval()/solve() read the input (the input and state grids need not
  // coincide).  We first sanity-check the raw write path (set/get round-trip), then
  // the exact transfer.
  if( !pu.empty() ){
    bool const wrote = oc.set_input_values( u_ic, pu, inp.data() );
    check_true( "set_input_values(u_ic, profile) succeeded", wrote );
    auto rt = oc.get_input_values( u_ic, inp.data() );
    bool rtok = ( rt.size() == pu.size() );
    if( rtok )
      for( size_t i = 0; i < rt.size(); ++i )
        rtok = rtok && ( std::abs( rt[i] - pu[i] ) <= 1e-14 );
    check_true( "set/get_input_values round-trip exact", rtok );

    // Exact transfer: seed u_ic from u's terminal at u_ic's own nodes.
    bool const xok = oc.transfer_terminal( u, u_ic, xv.data(), inp.data() );
    check_true( "transfer_terminal(u -> u_ic) succeeded", xok );
    if( xok ){
      // The re-seeded IC input now represents u(t=UB,.) : check it reproduces the
      // exact terminal at several interior points (interpolated, so discretisation tol).
      double rmax = 0.;
      for( double frac : { 0.15, 0.35, 0.55, 0.75, 0.95 } ){
        double const zt = frac * p.zf;
        OCFESLV::t_Coord pt; pt[z] = zt;
        double const uic = oc.eval_colloc<double>( u_ic, pt, xv.data(), inp.data() );
        rmax = std::max( rmax, std::abs( uic - u_exact( p.tf, zt, p ) ) );
      }
      check_close( "re-seeded IC input represents the terminal u(UB,.)", rmax, 0., 1e-3 );
    }
  }

  return pu;
}

} // namespace

// ---------------------------------------------------------------------------
// Drive an actual element-by-element march on a ONE-evolution-element template and
// validate the marched trajectory against the exact solution.  Each step: set the
// inflow BC over the element's absolute-time window, solve, compare the terminal to
// exact, then transfer the terminal onto the IC input for the next element.
namespace {
static bool run_march( Par const& p, std::vector<double> widths = {} )
{
  // Per-element evolution widths (must sum to tf); empty => uniform tf/NEL_T.
  std::vector<double> W = widths;
  if( W.empty() ) W.assign( NEL_T, p.tf / double( NEL_T ) );
  bool const nonuniform = ( widths.size() == NEL_T );
  std::cout << "\n============================================================\n"
            << "  march loop: one-element template, " << W.size() << " steps"
            << ( nonuniform ? "  (NON-UNIFORM widths)" : "  (uniform)" ) << "\n"
            << "============================================================\n";
  double const w0 = W.front();

  FFGraph DAG;
  FFVar t    = DAG.add_var( "t" );
  FFVar z    = DAG.add_var( "z" );
  FFVar u    = DAG.add_var( "u(t,z)" );
  FFVar u_ic = DAG.add_var( "u_ic(z)" );
  FFVar u_bc = DAG.add_var( "u_bc(t)" );
  FFPartial OpP;
  FFDom zdom( 0., p.zf, NEL_Z, FFDom::LGL, NZ );

  OCFESLV oc( &DAG );
  oc.add_domain( t, FFDom( 0., w0, 1, FFDom::LGR, NT ) );     // ONE evolution element [0,w0]
  oc.add_domain( z, zdom );
  oc.add_state ( u, { t, z } );
  oc.add_input ( u_ic, { z } );                              // IC   (transferred per step)
  oc.add_input ( u_bc, { t } );                              // inflow BC (absolute-time)
  oc.update_ref( u,    [&]( OCFESLV::t_Coord const& c ){ return u_exact( c.at(t), c.at(z), p ); } );
  oc.update_ref( u_ic, [&]( OCFESLV::t_Coord const& c ){ return u_exact( 0., c.at(z), p ); } );
  oc.update_ref( u_bc, [&]( OCFESLV::t_Coord const& c ){ return u_exact( c.at(t), 0., p ); } );

  FFVar PDE = OpP( u, t ) + p.a * OpP( u, z );
  FFVar INI = u - u_ic;
  FFVar BC  = u - u_bc;
  oc.add_equation( PDE, { t, z }, { FFDom::ALL - FFDom::LB, FFDom::ALL - FFDom::LB },
                   OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( INI, { t, z }, { FFDom::LB,              FFDom::ALL              },
                   OCFESLV::EqnOptions( OCFESLV::EqnRole::INITIAL, 0 ) );
  oc.add_equation( BC,  { t, z }, { FFDom::ALL - FFDom::LB, FFDom::LB               },
                   OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );
  oc.options.INTERFACE.TYPE = OCFESLV::Options::IC_AUTO;
  oc.options.CLASSIFY.MODE       = OCFESLV::Options::CLASS_AUTO;
  oc.options.DISPLAY_LEVEL  = 0;                              // quiet during the march
  oc.options.SOLVE.MARCHING = false;                         // MANUAL march: each solve() is a single-window solve
  oc.set_evolution_domain( t );

  if( !oc.setup() ){ std::cerr << "  template setup() FAILED\n"; ++g_fail; return false; }

  std::vector<double> xv, inp;
  if( !oc.init( xv, inp ) ){ std::cerr << "  init() FAILED\n"; ++g_fail; return false; }

  // Initial condition for element 0: u(0,.) sampled at the IC input's own nodes.
  oc.sample_input( u_ic, [&]( OCFESLV::t_Coord const& c ){ return u_exact( 0., c.at(z), p ); },
                   inp.data() );

  // ---- warm-start primitive: broadcast makes the state constant-in-evolution ----
  {
    auto const zn = physical_nodes( oc.var_domain().at( z ) );
    std::vector<double> icp( zn.size() );
    for( size_t i = 0; i < zn.size(); ++i ) icp[i] = u_exact( 0., zn[i], p );
    std::vector<double> vw( xv.size(), 0. );
    bool const wb = oc.warmstart_broadcast( u, icp, vw.data() );
    check_true( "warmstart_broadcast succeeded", wb );
    if( wb ){
      double dmax = 0., emax = 0.;
      for( double frac : { 0.2, 0.5, 0.8 } ){
        double const zt = frac * p.zf;
        OCFESLV::t_Coord pa; pa[t] = 0.25 * w0; pa[z] = zt;
        OCFESLV::t_Coord pb; pb[t] = 0.85 * w0; pb[z] = zt;
        double const a = oc.eval_colloc<double>( u, pa, vw.data(), nullptr );
        double const b = oc.eval_colloc<double>( u, pb, vw.data(), nullptr );
        dmax = std::max( dmax, std::abs( a - b ) );                     // constant in evolution
        emax = std::max( emax, std::abs( a - u_exact( 0., zt, p ) ) );  // equals the IC profile
      }
      check_close( "broadcast state is constant in evolution", dmax, 0., 1e-12 );
      check_close( "broadcast state equals the IC profile",    emax, 0., 1e-3  );
    }
  }

  double maxerr = 0.; bool all_conv = true; bool all_width = true;
  double tk = 0.;
  for( size_t k = 0; k < W.size(); ++k ){
    double const wk = W[k];
    // Rewrite the template's evolution width to this element's width (no re-setup).
    bool const wok = oc.reset_evolution_width( wk );
    all_width = all_width && wok;
    // inflow BC over the element's absolute-time window [tk, tk+wk].
    oc.sample_input( u_bc, [&]( OCFESLV::t_Coord const& c ){ return u_exact( tk + c.at(t), 0., p ); },
                     inp.data() );
    // (u_ic already holds u(tk,.): the initial profile for k=0, the transferred terminal for k>0.)

    OCFESLV::SolveReport const rep = oc.solve( xv.data(), inp.data() );
    all_conv = all_conv && rep.converged;

    // Terminal of this element == u(tk+wk,.) : compare to exact at a few points.
    double eterm = 0.;
    for( double frac : { 0.1, 0.3, 0.5, 0.7, 0.9 } ){
      double const zt = frac * p.zf;
      OCFESLV::t_Coord pt; pt[z] = zt;
      eterm = std::max( eterm, std::abs( oc.terminal_value( u, pt, xv.data(), inp.data() )
                                         - u_exact( tk + wk, zt, p ) ) );
    }
    maxerr = std::max( maxerr, eterm );
    std::cout << "  step " << k << "  window [" << std::fixed << std::setprecision(3)
              << tk << "," << ( tk + wk ) << "]  w=" << wk
              << "  converged=" << ( rep.converged ? "y" : "n" )
              << "  max|term-exact|=" << std::scientific << std::setprecision(3) << eterm << "\n";

    // Transfer this element's terminal onto the IC input for the next element.
    oc.transfer_terminal( u, u_ic, xv.data(), inp.data() );
    tk += wk;
  }
  if( nonuniform ) check_true( "all reset_evolution_width() succeeded", all_width );

  check_true ( "all march steps converged", all_conv );
  check_close( "marched trajectory matches exact (terminals)", maxerr, 0., 1e-3 );
  return g_fail == 0;
}

// solve()-driven march: SOLVE_MARCHING coordinates setup (collapse) + solve (block march).
// IC posed as an input (transferred per element); symbolic inflow BC auto-evaluates at each
// window.  Validates convergence and the marched final terminal against exact.
static bool run_solve_marching( Par const& p, OCFESLV::Options::SolveWarmstart ws, char const* wslabel )
{
  std::cout << "\n------------------------------------------------------------\n"
            << "  solve()-driven march  (warm-start = " << wslabel << ")\n"
            << "------------------------------------------------------------\n";
  FFGraph DAG;
  FFVar t = DAG.add_var( "t" ), z = DAG.add_var( "z" ), u = DAG.add_var( "u(t,z)" );
  FFVar u_ic = DAG.add_var( "u_ic(z)" );
  FFPartial OpP;
  OCFESLV oc( &DAG );
  oc.add_domain( t, FFDom( 0., p.tf, NEL_T, FFDom::LGR, NT ) );   // N-element decl = march grid
  oc.add_domain( z, FFDom( 0., p.zf, NEL_Z, FFDom::LGL, NZ ) );
  oc.add_state ( u, { t, z } );
  oc.add_input ( u_ic, { z } );                                  // IC input (transferred per element)
  oc.update_ref( u,    [&]( OCFESLV::t_Coord const& c ){ return u_exact( c.at(t), c.at(z), p ); } );
  oc.update_ref( u_ic, [&]( OCFESLV::t_Coord const& c ){ return u_exact( 0., c.at(z), p ); } );
  FFVar PDE = OpP( u, t ) + p.a * OpP( u, z );
  FFVar INI = u - u_ic;                                          // IC via the transferred input
  FFVar BC  = u - p.u0 * sin( p.k * ( 0.0 - p.a * t ) );         // symbolic inflow (auto per window)
  oc.add_equation( PDE, { t, z }, { FFDom::ALL - FFDom::LB, FFDom::ALL - FFDom::LB },
                   OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( INI, { t, z }, { FFDom::LB,              FFDom::ALL              },
                   OCFESLV::EqnOptions( OCFESLV::EqnRole::INITIAL, 0 ) );
  oc.add_equation( BC,  { t, z }, { FFDom::ALL - FFDom::LB, FFDom::LB               },
                   OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );
  oc.options.INTERFACE.TYPE  = OCFESLV::Options::IC_AUTO;
  oc.options.CLASSIFY.MODE        = OCFESLV::Options::CLASS_AUTO;
  oc.options.DISPLAY_LEVEL   = 0;
  oc.options.SOLVE.VERBOSE   = true;     // per-step [march] trace
  oc.options.SOLVE.WARMSTART = ws;        // SOLVE_MARCHING is on by default; a differential state makes it available
  oc.set_evolution_domain ( t );
  // NO set_marching_transfer: u is auto-detected as the differential state, and its canonical IC
  // u - u_ic auto-resolves the (u, u_ic) transfer.

  if( !oc.setup() ){ std::cerr << "  setup() FAILED\n"; ++g_fail; return false; }
  // Collapse happened by default (SOLVE_MARCHING on + differential state present -> available).
  check_true( "collapsed to one evolution element", oc.n_evolution_elem() == 1 );
  check_true( "march grid recorded (N+1 boundaries)", oc.march_grid().size() == NEL_T + 1 );
  check_true( "n_march_steps() == NEL_T", oc.n_march_steps() == NEL_T );

  std::vector<double> xv, inp;
  if( !oc.init( xv, inp ) ){ std::cerr << "  init() FAILED\n"; ++g_fail; return false; }
  double const* ip = inp.empty() ? nullptr : inp.data();

  OCFESLV::SolveReport const rep = oc.solve( xv.data(), ip );
  check_true( "solve()-march converged", rep.converged );

  // After the march xv holds the LAST element's solution on [t_{N-1}, tf]; its terminal
  // (at tf) must match exact -- errors from every step have propagated into it.
  double emax = 0.;
  for( double frac : { 0.1, 0.3, 0.5, 0.7, 0.9 } ){
    double const zt = frac * p.zf;
    OCFESLV::t_Coord pt; pt[z] = zt;
    emax = std::max( emax, std::abs( oc.terminal_value( u, pt, xv.data(), ip )
                                     - u_exact( p.tf, zt, p ) ) );
  }
  check_close( "solve()-march final terminal vs exact", emax, 0., 1e-3 );
  return g_fail == 0;
}
} // namespace

int main()
{
  std::cout << "\n============================================================\n"
            << "  OCFE_marching: Stage-0 terminal-read + slab-map primitives\n"
            << "  (scalar linear advection, EVOL_HYPERBOLIC)\n"
            << "============================================================\n";

  Par p;

  std::vector<double> prof_tz = run( /*evolution_first=*/true,  p );
  std::vector<double> prof_zt = run( /*evolution_first=*/false, p );

  std::cout << "\n------------------------------------------------------------\n"
            << "  cross-layout agreement (order-agnostic check)\n"
            << "------------------------------------------------------------\n";
  bool const same_size = ( !prof_tz.empty() && prof_tz.size() == prof_zt.size() );
  check_true( "both layouts produced a terminal profile of equal size", same_size );
  if( same_size ){
    double dmax = 0.;
    for( size_t i = 0; i < prof_tz.size(); ++i )
      dmax = std::max( dmax, std::abs( prof_tz[i] - prof_zt[i] ) );
    // Same physics, different DOF layout: the terminal profiles must match to
    // solver tolerance (this is the real test that nothing assumes evolution-first).
    check_close( "terminal_profile identical across (t,z) and (z,t)", dmax, 0., 1e-9 );
  }

  run_march( p );                                                    // uniform grid
  run_march( p, { 0.20 * p.tf, 0.35 * p.tf, 0.45 * p.tf } );         // non-uniform grid (sums to tf)

  run_solve_marching( p, OCFESLV::Options::REUSE,        "REUSE"        );
  run_solve_marching( p, OCFESLV::Options::BROADCAST_IC, "BROADCAST_IC" );

  std::cout << "\n============================================================\n"
            << "  RESULT: " << g_pass << " passed, " << g_fail << " failed  -- "
            << ( g_fail == 0 ? "ALL PASS" : "FAILURES" ) << "\n"
            << "============================================================\n";
  return g_fail == 0 ? 0 : 1;
}
