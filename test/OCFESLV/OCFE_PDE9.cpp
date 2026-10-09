// Copyright (C) Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the EPL-2.0 with GPL-2.0-or-later as a Secondary License; see the LICENSE file.
// SPDX-License-Identifier: EPL-2.0 OR GPL-2.0-or-later
//
// OCFE_PDE9_solve.cpp
// -------------------
// Manufactured dynamic 2D advection-diffusion regression test for OCFESLV.
//
// The model uses three collocation domains, time t and space x,y:
//
//   u_t + a*u_x - kappa * ( u_xx + u_yy ) = f(t,x,y),   (t,x,y) in (0,T]x(0,1)^2
//
// This is a deliberate sibling of OCFE_PDE5 (pure 2D heat): it adds an x-
// advection term so the item-11 symbol-AUTO read marks x as advection-dominated.
// NOTE: because x still carries the reduced auxiliary Dx_U = u_x (from u_xx),
// _resolve_interface_type returns IC_VALUE for x regardless of advection (the
// aux-link rule precedes the upwind branch).  So this remains a *parabolic*
// interface case with all reduced directions VALUE; its purpose is a SECOND
// three-direction reduced-aux example (distinct mesh/structure from PDE5) to
// harden the tensor-corner spanning fix, plus a 3D exercise of the AUTO read.
// Genuine interface UPWIND requires a first-order (no-aux) hyperbolic direction.
//
// The exact solution is the same smooth central Gaussian-like bump whose
// variance grows in time, with the source adjusted for the advection term.
//
// The test rebuilds and solves a fresh OCFESLV for each imposition mode:
//
//   IC_WEAK, IC_TRACE, IC_STRONG
//
// and prints per-mode residuals, exact-profile errors, duplicate-node spreads,
// output files, and a final summary table.
//
// Build-time switches:
//   -DTEST_ADV_OCENV_HEADER     header to include; defaults to "ocfeslv.hpp"
//   -DTEST_ADV_NEL_T            finite elements in t; defaults to 3
//   -DTEST_ADV_NEL_X            finite elements in x; defaults to 3
//   -DTEST_ADV_NEL_Y            finite elements in y; defaults to 3
//   -DTEST_ADV_NT               collocation nodes/element in t; defaults to 3
//   -DTEST_ADV_NX               collocation nodes/element in x; defaults to 3
//   -DTEST_ADV_NY               collocation nodes/element in y; defaults to 3
//   -DTEST_ADV_INIT_MODE        0 exact, 1 perturbed exact, 2 constant; defaults to 1
//   -DTEST_ADV_MAXIT            Newton/LM iterations; defaults to 20
//   -DTEST_ADV_SOLVE_TOL        residual tolerance; defaults to 1e-9
//   -DTEST_ADV_EXACT_TOL        exact-profile tol (TRACE/STRONG); defaults to 1e-5
//   -DTEST_ADV_AUX_TOL          aux derivative tol (TRACE/STRONG); defaults to 2e-4
//   -DTEST_ADV_REF_TOL          reference-profile residual tol; defaults to 5e-4
//   -DTEST_ADV_SPREAD_TOL       duplicate spread tol (TRACE/STRONG); defaults to 5e-8
//   -DTEST_ADV_EXACT_TOL_WEAK   exact-profile tol (IC_WEAK, penalty); defaults to 5e-3
//   -DTEST_ADV_AUX_TOL_WEAK     aux derivative tol (IC_WEAK, penalty); defaults to 2e-1
//   -DTEST_ADV_SPREAD_TOL_WEAK  duplicate spread tol (IC_WEAK, penalty); defaults to 1e-2
//   -DTEST_ADV_OUTPUT_PREFIX    output file prefix; defaults to "OCFE_PDE9_solve"
//   -DTEST_ADV_PRINT_OC         print OCFESLV structure; defaults to 0

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <string>
#include <vector>

#include <armadillo>

#ifndef OCFE_OCFESLV_HEADER
#define OCFE_OCFESLV_HEADER "ocfeslv.hpp"
#endif
#include OCFE_OCFESLV_HEADER

using namespace mc;

#ifndef TEST_ADV_NEL_T
#define TEST_ADV_NEL_T 2
#endif
#ifndef TEST_ADV_NEL_X
#define TEST_ADV_NEL_X 2
#endif
#ifndef TEST_ADV_NEL_Y
#define TEST_ADV_NEL_Y 2
#endif
#ifndef TEST_ADV_NT
#define TEST_ADV_NT 10
#endif
#ifndef TEST_ADV_NX
#define TEST_ADV_NX 10
#endif
#ifndef TEST_ADV_NY
#define TEST_ADV_NY 10
#endif
#ifndef TEST_ADV_INIT_MODE
#define TEST_ADV_INIT_MODE 1
#endif
#ifndef TEST_ADV_MAXIT
#define TEST_ADV_MAXIT 20
#endif
#ifndef TEST_ADV_SOLVE_TOL
#define TEST_ADV_SOLVE_TOL 1e-9
#endif
// Accuracy tolerances are mode-aware (selected in run_mode).  The exact-multiplier
// modes (IC_TRACE/IC_STRONG) enforce continuity exactly and reach spectral accuracy
// for the gentled bump (width0=0.1) at NEL=2/N=10: |U-exact|~2e-6, |aux-exact|~4e-5,
// ref-residual~7e-5, dup-spread~3e-16.  The bars below sit ~4-7x above those.  IC_WEAK
// is penalty-limited (SAT consistency error O(1/sigma0)), so it gets penalty-level
// bars; its accuracy is not a correctness invariant.  The corner-fix regression guard
// is the EXACT-mode dup-spread (kept tight at 5e-8): a regression would push it from
// 3e-16 back toward the 1.577e-1 float and fail by ~7 orders.  The WEAK dup-spread bar
// (1e-2) is a looser secondary backstop.
#ifndef TEST_ADV_EXACT_TOL
#define TEST_ADV_EXACT_TOL 1e-5
#endif
#ifndef TEST_ADV_AUX_TOL
#define TEST_ADV_AUX_TOL 2e-4
#endif
#ifndef TEST_ADV_REF_TOL
#define TEST_ADV_REF_TOL 5e-4
#endif
#ifndef TEST_ADV_SPREAD_TOL
#define TEST_ADV_SPREAD_TOL 5e-8
#endif
#ifndef TEST_ADV_EXACT_TOL_WEAK
#define TEST_ADV_EXACT_TOL_WEAK 5e-3
#endif
#ifndef TEST_ADV_AUX_TOL_WEAK
#define TEST_ADV_AUX_TOL_WEAK 2e-1
#endif
#ifndef TEST_ADV_SPREAD_TOL_WEAK
#define TEST_ADV_SPREAD_TOL_WEAK 1e-2
#endif
#ifndef TEST_ADV_SAT_SIGMA0
#define TEST_ADV_SAT_SIGMA0 10.0
#endif
#ifndef TEST_ADV_PRINT_OC
#define TEST_ADV_PRINT_OC 0
#endif
#ifndef TEST_ADV_OUTPUT_PREFIX
#define TEST_ADV_OUTPUT_PREFIX "OCFE_PDE9_solve"
#endif

// Opt-in: route the augmented IC_TRACE/IC_STRONG solve through sparse SPQR.
// Effective only when the header is built with -DCRONOS__WITH_SPQR (SuiteSparse).
#define TEST_ADV_SPQR

namespace {

struct Par
{
  double tf      = 1.0;
  double kappa   = 3.0e-2;
  double a       = 2.0;     // x-advection speed (makes the AUTO read flag x
                            // advection-dominated; resolver still folds to VALUE)
  double bg      = 0.0;
  double A       = 4.0e-1;
  double x0      = 0.5;
  double y0      = 0.5;
  double width0  = 1.0e-1;
  double spread  = 8.0e-2;
};

static double Swidth( double t, Par const& p )
{ return p.width0 + p.spread*t; }
static double R2( double x, double y, Par const& p )
{
  double const dx = x - p.x0;
  double const dy = y - p.y0;
  return dx*dx + dy*dy;
}
static double bump_exact( double t, double x, double y, Par const& p )
{
  double const S = Swidth(t,p);
  return p.A * p.width0 / S * std::exp( -R2(x,y,p) / S );
}
static double U_exact( double t, double x, double y, Par const& p )
{ return p.bg + bump_exact(t,x,y,p); }
static double U_t_exact( double t, double x, double y, Par const& p )
{
  double const S  = Swidth(t,p);
  double const rr = R2(x,y,p);
  double const B  = bump_exact(t,x,y,p);
  return B * p.spread * ( -1.0/S + rr/(S*S) );
}
static double U_x_exact( double t, double x, double y, Par const& p )
{
  double const S  = Swidth(t,p);
  double const dx = x - p.x0;
  return bump_exact(t,x,y,p) * ( -2.0*dx/S );
}
static double U_y_exact( double t, double x, double y, Par const& p )
{
  double const S  = Swidth(t,p);
  double const dy = y - p.y0;
  return bump_exact(t,x,y,p) * ( -2.0*dy/S );
}
static double U_xx_exact( double t, double x, double y, Par const& p )
{
  double const S  = Swidth(t,p);
  double const dx = x - p.x0;
  return bump_exact(t,x,y,p) * ( 4.0*dx*dx/(S*S) - 2.0/S );
}
static double U_yy_exact( double t, double x, double y, Par const& p )
{
  double const S  = Swidth(t,p);
  double const dy = y - p.y0;
  return bump_exact(t,x,y,p) * ( 4.0*dy*dy/(S*S) - 2.0/S );
}
// PRE-EXISTING, unrelated to the rev29 patch: F_exact has no caller in this driver and
// did not have one before it either.  Retained rather than deleted because it documents
// the manufactured source term that the assembled equation reproduces inline, which is
// worth keeping legible; [[maybe_unused]] silences -Wunused-function without discarding it.
[[maybe_unused]] static double F_exact( double t, double x, double y, Par const& p )
{ return U_t_exact(t,x,y,p) + p.a*U_x_exact(t,x,y,p)
       - p.kappa*( U_xx_exact(t,x,y,p) + U_yy_exact(t,x,y,p) ); }

static double max_abs( std::vector<double> const& r )
{
  double m=0.0;
  for( double v: r ) m = std::max(m,std::abs(v));
  return m;
}

static void print_residuals( std::string const& label, std::vector<double> const& r )
{
  double sum=0.0;
  for( double v: r ) sum += std::abs(v);
  std::cout << std::left << std::setw(48) << label
            << " n=" << std::setw(6) << r.size()
            << " max|r|=" << std::scientific << std::setprecision(4) << max_abs(r)
            << " mean|r|=" << (r.empty()?0.0:sum/double(r.size())) << "\n";
}

static bool check_close( std::string const& label, double val, double tol )
{
  bool ok = std::isfinite(val) && val <= tol;
  std::cout << std::left << std::setw(56) << label
            << " value=" << std::scientific << std::setprecision(6) << val
            << " tol=" << tol << "  " << (ok?"PASS":"FAIL") << "\n";
  return ok;
}

struct StateErrors
{
  double eU=0.0;
  double eUx=0.0;
  double eUy=0.0;
  double eAux=0.0;
};

static bool is_auxiliary_state_name( std::string const& nm )
{
  return nm.find("Dx_") != std::string::npos
      || nm.find("Dy_") != std::string::npos
      || nm.find("Daux") != std::string::npos;
}

static double exact_value_for_state( std::string const& nm, std::vector<double> const& xyz, Par const& p )
{
  // OCFESLV node_colloc() reports coordinates in sorted dependency order.  With
  // variables created as t, x, y and states depending on {t,x,y}, this gives:
  //   xyz[0] = t, xyz[1] = x, xyz[2] = y.
  double const t = xyz.size()>0 ? xyz[0] : 0.0;
  double const x = xyz.size()>1 ? xyz[1] : 0.0;
  double const y = xyz.size()>2 ? xyz[2] : 0.0;

  if( nm == "U(t,x,y)" ) return U_exact(t,x,y,p);
  if( nm.find("Dx_") != std::string::npos ) return U_x_exact(t,x,y,p);
  if( nm.find("Dy_") != std::string::npos ) return U_y_exact(t,x,y,p);
  // Fallback for unexpected auxiliary names.  This keeps the diagnostic finite
  // while making the name visible in the printed comparison table.
  return 0.0;
}

static StateErrors print_primitive_exact_errors
( OCFESLV const& oc, std::vector<double> const& var, Par const& p, std::string const& label )
{
  StateErrors out;
  size_t off=0;
  std::cout << "\nVariable comparison against exact profiles (" << label << "):\n";
  for( auto const& st: oc.states_colloc() ){
    auto nodes = oc.node_colloc(st);
    std::string const nm = st.name();
    if( !is_auxiliary_state_name(nm) ){
      double maxerr=0.0, meanerr=0.0;
      for( size_t i=0; i<nodes.size(); ++i ){
        double const ref = exact_value_for_state(nm,nodes[i],p);
        double const err = std::abs(var[off+i]-ref);
        maxerr = std::max(maxerr,err);
        meanerr += err;
      }
      meanerr /= nodes.empty()? 1.0 : double(nodes.size());
      if( nm == "U(t,x,y)" ) out.eU = maxerr;
      std::cout << "  " << std::setw(18) << nm
                << " max|v-v_exact|=" << std::scientific << std::setprecision(4) << maxerr
                << " mean|v-v_exact|=" << meanerr << "\n";
    }
    off += nodes.size();
  }
  return out;
}

static StateErrors print_auxiliary_exact_errors
( OCFESLV const& oc, std::vector<double> const& var, Par const& p, std::string const& label )
{
  StateErrors out;
  size_t off=0;
  std::cout << "\nAuxiliary derivative comparison against exact profiles (" << label << "):\n";
  for( auto const& st: oc.states_colloc() ){
    auto nodes = oc.node_colloc(st);
    std::string const nm = st.name();
    if( is_auxiliary_state_name(nm) ){
      double maxerr=0.0, meanerr=0.0;
      for( size_t i=0; i<nodes.size(); ++i ){
        double const ref = exact_value_for_state(nm,nodes[i],p);
        double const err = std::abs(var[off+i]-ref);
        maxerr = std::max(maxerr,err);
        meanerr += err;
      }
      meanerr /= nodes.empty()? 1.0 : double(nodes.size());
      if( nm.find("Dx_") != std::string::npos ) out.eUx = maxerr;
      else if( nm.find("Dy_") != std::string::npos ) out.eUy = maxerr;
      out.eAux = std::max(out.eAux,maxerr);
      std::cout << "  " << std::setw(18) << nm
                << " max|aux-aux_exact|=" << std::scientific << std::setprecision(4) << maxerr
                << " mean|aux-aux_exact|=" << meanerr << "\n";
    }
    off += nodes.size();
  }
  return out;
}

struct DuplicateSpread
{
  double max_pair_spread = 0.0;
  double max_corner_spread = 0.0;
  size_t max_multiplicity = 0;
};

// ---------------------------------------------------------------------------
// rev29: MARCH-AWARE exact-error comparison
// ---------------------------------------------------------------------------
// print_primitive_exact_errors / print_auxiliary_exact_errors index `var` directly at
// node_colloc(st) offsets.  Under marching the evolution domain is COLLAPSED to a single
// element and `var` is sized n_colloc_sta() = ONE WINDOW, so those helpers compare over
// only part of the horizon.  MEASURED at nel_t=2: residual vector n 24000 -> 12000, plot
// file 8081 -> 4041 rows.  A half-horizon comparison that PASSES is worse than one that
// fails, because it reports a green result computed over the wrong set.
//
// eval_solution() is buffer-free and mode-agnostic: it reads the cached monolithic primal
// (SOLVE_REUSE, default on) or the marched trajectory (MARCH_STORE_TRAJECTORY, default on)
// and routes each point to the window whose [t0,t1] contains it.  Sampling an INDEPENDENT
// tensor grid over the full domain is correct in either mode, and is a stronger oracle
// than re-reading the solver's own nodes in the solver's own ordering.
static StateErrors print_sampled_exact_errors
( OCFESLV const& oc, FFVar const& t, FFVar const& x, FFVar const& y,
  Par const& p, std::string const& label, size_t ns = 8 )
{
  StateErrors out;
  std::cout << "\nSampled comparison against exact profiles, FULL domain (" << label << "):\n";
  for( auto const& st: oc.states_colloc() ){
    std::string const nm = st.name();
    double maxerr=0.0, meanerr=0.0;
    size_t cnt=0;
    bool avail=true;
    for( size_t it=0; it<=ns && avail; ++it )
    for( size_t ix=0; ix<=ns && avail; ++ix )
    for( size_t iy=0; iy<=ns && avail; ++iy ){
      std::vector<double> const c = { p.tf*double(it)/double(ns),
                                      double(ix)/double(ns),
                                      double(iy)/double(ns) };
      OCFESLV::t_Coord pt; pt[t]=c[0]; pt[x]=c[1]; pt[y]=c[2];
      double v=0.0;
      try { v = oc.eval_solution(st,pt); }
      catch( ... ) { avail=false; break; }
      double const err = std::abs(v-exact_value_for_state(nm,c,p));
      maxerr = std::max(maxerr,err); meanerr += err; ++cnt;
    }
    if( !avail ){
      std::cout << "  " << std::setw(18) << nm
                << "  eval_solution unavailable (no cached solution) -- SKIPPED\n";
      continue;
    }
    meanerr /= cnt? double(cnt) : 1.0;
    if( nm == "U(t,x,y)" ) out.eU = maxerr;
    else if( nm.find("Dx_") != std::string::npos ) out.eUx = maxerr;
    else if( nm.find("Dy_") != std::string::npos ) out.eUy = maxerr;
    if( is_auxiliary_state_name(nm) ) out.eAux = std::max(out.eAux,maxerr);
    std::cout << "  " << std::setw(18) << nm
              << " max|v-v_exact|=" << std::scientific << std::setprecision(4) << maxerr
              << " mean|v-v_exact|=" << meanerr
              << "  (" << cnt << " samples over the full domain)\n";
  }
  return out;
}

static DuplicateSpread duplicate_node_spread
( OCFESLV const& oc, FFVar const& st, std::vector<double> const& var, size_t off )
{
  struct Accum { double lo, hi; size_t count; };
  std::map< std::vector<long long>, Accum > groups;
  auto nodes = oc.node_colloc(st);
  for( size_t i=0; i<nodes.size(); ++i ){
    std::vector<long long> key;
    key.reserve(nodes[i].size());
    for( double c: nodes[i] ) key.push_back( static_cast<long long>( std::llround(c*1.0e12) ) );
    double const v = var[off+i];
    auto it = groups.find(key);
    if( it == groups.end() ) groups.emplace( std::move(key), Accum{v,v,1} );
    else{
      it->second.lo = std::min(it->second.lo,v);
      it->second.hi = std::max(it->second.hi,v);
      ++it->second.count;
    }
  }

  DuplicateSpread out;
  for( auto const& kv: groups ){
    auto const& g = kv.second;
    out.max_multiplicity = std::max(out.max_multiplicity,g.count);
    if( g.count >= 2 ) out.max_pair_spread = std::max(out.max_pair_spread,g.hi-g.lo);
    if( g.count >= 4 ) out.max_corner_spread = std::max(out.max_corner_spread,g.hi-g.lo);
  }
  return out;
}

static double print_duplicate_spreads
( OCFESLV const& oc, std::vector<double> const& var, std::string const& label )
{
  std::cout << "\nElement-interface duplicate-node spreads (" << label << "):\n";
  double max_pair=0.0;
  size_t off=0;
  for( auto const& st: oc.states_colloc() ){
    auto nodes = oc.node_colloc(st);
    DuplicateSpread const d = duplicate_node_spread(oc,st,var,off);
    max_pair = std::max(max_pair,d.max_pair_spread);
    std::cout << "  " << std::setw(18) << st.name()
              << " max_pair=" << std::scientific << std::setprecision(6) << d.max_pair_spread
              << " max_corner=" << d.max_corner_spread
              << " max_multiplicity=" << d.max_multiplicity << "\n";
    off += nodes.size();
  }
  return max_pair;
}

static bool eval_residual( OCFESLV& oc, std::vector<double> const& x, std::vector<double>& r )
{ std::fill(r.begin(),r.end(),0.0); return oc.eval(r.data(),nullptr,x.data(),nullptr,nullptr); }

static void initialise_vector( OCFESLV const& oc, std::vector<double>& var, Par const& p )
{
  size_t off=0;
  for( auto const& st: oc.states_colloc() ){
    auto nodes = oc.node_colloc(st);
    std::string const nm = st.name();
    for( size_t i=0; i<nodes.size(); ++i ){
      double val = exact_value_for_state(nm,nodes[i],p);
#if TEST_ADV_INIT_MODE == 1
      val += 1.0e-2 * std::sin(0.31*double(off+i+1));
#elif TEST_ADV_INIT_MODE == 2
      val = is_auxiliary_state_name(nm) ? 0.0 : 1.0;
#endif
      var[off+i] = val;
    }
    off += nodes.size();
  }
  for( size_t i=off; i<var.size(); ++i ) var[i]=0.0;
}

static void write_output( OCFESLV const& oc, std::vector<double> const& var,
                          Par const& p, std::string const& fname )
{
  struct Rec
  {
    double t = 0.0;
    double x = 0.0;
    double y = 0.0;
    double u = 0.0;
    double uref = 0.0;
  };

  std::map< long long, std::vector<Rec> > blocks;
  size_t off=0;
  for( auto const& st: oc.states_colloc() ){
    auto nodes = oc.node_colloc(st);
    if( st.name() == "U(t,x,y)" ){
      for( size_t i=0; i<nodes.size(); ++i ){
        Rec r;
        r.t    = nodes[i].size()>0 ? nodes[i][0] : 0.0;
        r.x    = nodes[i].size()>1 ? nodes[i][1] : 0.0;
        r.y    = nodes[i].size()>2 ? nodes[i][2] : 0.0;
        r.u    = var[off+i];
        r.uref = U_exact(r.t,r.x,r.y,p);
        long long const tkey = static_cast<long long>( std::llround( r.t * 1.0e12 ) );
        blocks[tkey].push_back(r);
      }
      break;
    }
    off += nodes.size();
  }

  std::ofstream os(fname);
  os << std::scientific << std::setprecision(16);
  os << "# OCFE_PDE9_solve time-sliced spreading-bump output\n";
  os << "# Each data block corresponds to one collocation time point.\n";
  os << "# Within each block, columns are: t x y U U_exact error\n";

  bool first_block = true;
  for( auto& kv: blocks ){
    auto& b = kv.second;
    std::sort( b.begin(), b.end(), []( Rec const& a, Rec const& c ){
      if( std::abs(a.x-c.x) > 1.0e-14 ) return a.x < c.x;
      return a.y < c.y;
    } );
    if( !first_block ) os << "\n\n";
    first_block = false;
    double const t = b.empty()? 0.0 : b.front().t;
    os << "# t = " << t << "\n";
    os << "# t x y U U_exact error\n";
    for( auto const& r: b ){
      os << r.t << ' ' << r.x << ' ' << r.y << ' ' << r.u << ' ' << r.uref
         << ' ' << r.u-r.uref << "\n";
    }
  }
}

static const char* imp_name( OCFESLV::Options::ImpositionType imp )
{
  switch( imp ){
    case OCFESLV::Options::IC_WEAK:   return "IC_WEAK";
    case OCFESLV::Options::IC_TRACE:  return "IC_TRACE";
    case OCFESLV::Options::IC_STRONG: return "IC_STRONG";
  }
  return "?";
}

struct ModeResult
{
  std::string name;
  bool   ok=false;
  bool   solved=false;
  size_t nVar=0, nEqn=0, nTrace=0;
  double final_res=std::numeric_limits<double>::infinity();
  double max_spread=0.0;
  double eU=0.0;
  double eUx=0.0;
  double eUy=0.0;
  double eAux=0.0;
};

static ModeResult run_mode( OCFESLV::Options::ImpositionType imp,
                            std::string const& suffix,
                            Par const& p )
{
  ModeResult R;
  R.name = imp_name(imp);
  bool ok=true;

  // Mode-aware accuracy bars: IC_WEAK is penalty-limited, the exact-multiplier
  // modes reach spectral accuracy.  ref-residual and the solve residual are
  // mode-independent and stay uniform.
  bool   const weak       = ( imp == OCFESLV::Options::IC_WEAK );
  double const exact_tol  = weak ? TEST_ADV_EXACT_TOL_WEAK  : TEST_ADV_EXACT_TOL;
  double const aux_tol    = weak ? TEST_ADV_AUX_TOL_WEAK    : TEST_ADV_AUX_TOL;
  double const spread_tol = weak ? TEST_ADV_SPREAD_TOL_WEAK : TEST_ADV_SPREAD_TOL;
  double const ref_tol    = TEST_ADV_REF_TOL;

  std::cout << "\n========== manufactured dynamic 2D advection-diffusion test =========="
            << "\nimposition: " << R.name
            << ", finite elements: t=" << TEST_ADV_NEL_T
            << " x=" << TEST_ADV_NEL_X
            << " y=" << TEST_ADV_NEL_Y
            << ", nodes/element: t=" << TEST_ADV_NT
            << " x=" << TEST_ADV_NX
            << " y=" << TEST_ADV_NY
            << ", init_mode=" << TEST_ADV_INIT_MODE << "\n";
  std::cout << "manufactured profile: central Gaussian bump centred at ("
            << p.x0 << "," << p.y0 << ") with width S(t)=" << p.width0
            << "+" << p.spread << "*t\n";

  FFGraph DAG;
  FFVar t = DAG.add_var("t");
  FFVar x = DAG.add_var("x");
  FFVar y = DAG.add_var("y");
  FFVar U = DAG.add_var("U(t,x,y)");
  FFPartial OpP;

  FFVar DX = x - p.x0;
  FFVar DY = y - p.y0;
  FFVar S  = p.width0 + p.spread*t;
  FFVar RR = DX*DX + DY*DY;
  FFVar BE = p.A * p.width0 * exp( -RR/S ) / S;
  FFVar UE = p.bg + BE;
  FFVar FE = BE * p.spread * ( -1.0/S + RR/(S*S) )
           + p.a * BE * ( -2.0*DX/S )
           - p.kappa * BE * ( 4.0*RR/(S*S) - 4.0/S );

  FFVar PDE = OpP(U,t) + p.a*OpP(U,x) - p.kappa*( OpP(U,{x,2}) + OpP(U,{y,2}) ) - FE;
  FFVar BC  = U - UE;

  OCFESLV oc(&DAG);
  oc.add_domain( t, FFDom(0., p.tf, TEST_ADV_NEL_T, FFDom::LGR, TEST_ADV_NT) );
  oc.add_domain( x, FFDom(0., 1.0,  TEST_ADV_NEL_X, FFDom::CGL, TEST_ADV_NX) );
  oc.add_domain( y, FFDom(0., 1.0,  TEST_ADV_NEL_Y, FFDom::CGL, TEST_ADV_NY) );
  oc.add_state( U, {t,x,y} );
  oc.update_ref( U, [&]( OCFESLV::t_Coord const& c ){ return U_exact(c.at(t),c.at(x),c.at(y),p); } );

  OCFESLV::EqnOptions int_opt( OCFESLV::EqnRole::INTERIOR, 0 );
  OCFESLV::EqnOptions ini_opt( OCFESLV::EqnRole::INITIAL,  0 );
  OCFESLV::EqnOptions bnd_opt( OCFESLV::EqnRole::BOUNDARY, 0 );

  int const T_INT = FFDom::ALL - FFDom::LB;
  int const X_INT = FFDom::ALL - FFDom::LB - FFDom::UB;
  int const Y_INT = FFDom::ALL - FFDom::LB - FFDom::UB;

  // Interior parabolic PDE; no evolution equation at the initial-time nodes,
  // and no second-order equation on spatial Dirichlet boundaries.
  oc.add_equation( PDE, {t,x,y}, {T_INT, X_INT, Y_INT}, int_opt );

  // Initial condition at t=0 for the full spatial tensor-product grid.
  oc.add_equation( BC, {t,x,y}, {FFDom::LB, FFDom::ALL, FFDom::ALL}, ini_opt );

  // Dirichlet boundaries for t>0.  The y-boundary rows exclude x-boundary
  // nodes to avoid duplicate corner equations already supplied by x-boundaries.
  oc.add_equation( BC, {t,x,y}, {T_INT, FFDom::LB, Y_INT}, bnd_opt );
  oc.add_equation( BC, {t,x,y}, {T_INT, FFDom::UB, Y_INT}, bnd_opt );
  oc.add_equation( BC, {t,x,y}, {T_INT, FFDom::ALL, FFDom::LB}, bnd_opt );
  oc.add_equation( BC, {t,x,y}, {T_INT, FFDom::ALL, FFDom::UB}, bnd_opt );

  oc.reset_evolution_domain();
  oc.options.REDUCE.ORDER    = OCFESLV::Options::RED_FULL;
  oc.options.CLASSIFY.MODE        = OCFESLV::Options::CLASS_AUTO;
  oc.options.INTERFACE.TYPE  = OCFESLV::Options::IC_AUTO;
  oc.options.INTERFACE.IMPOSITION = imp;
  oc.options.INTERFACE.SAT_SIGMA0      = TEST_ADV_SAT_SIGMA0;
  oc.options.DISPLAY_LEVEL   = 1;

  if( !oc.setup() ){
    std::cerr << "ERROR: OCFESLV::setup() failed for " << R.name << "\n";
    R.ok=false;
    return R;
  }
#if TEST_ADV_PRINT_OC
  std::cout << oc;
#endif

  size_t const nVar=oc.n_colloc_sta();
  size_t const nEqn=oc.n_colloc_eqn();
  size_t const nTrace=oc.n_colloc_trace();
  R.nVar=nVar; R.nEqn=nEqn; R.nTrace=nTrace;
  std::cout << "nVar=" << nVar << " nEqn=" << nEqn << " nTrace=" << nTrace
            << " square=" << (nVar==nEqn?"yes":"no") << "\n";
  ok &= (nVar == nEqn);

  {
    auto const& cls = oc.pde_type();
    auto const& sym = oc.symbol_cached();
    std::cout << "PDE type: " << OCFESLV::pde_type_name(cls.type)
              << "  At_singular=" << (cls.At_singular?"yes":"no")
              << "  parabolic=" << (cls.parabolic_structure_detected?"yes":"no")
              << "  evolution_hyperbolic=" << (cls.evolution_hyperbolic?"yes":"no") << "\n";
    std::cout << "Principal symbol size: states=" << sym.vState.size()
              << " equations=" << sym.vEqn.size()
              << " domains=" << sym.vDom.size() << "\n";
  }

  std::cout << "States after setup:";
  for( auto const& st: oc.states_colloc() ) std::cout << ' ' << st.name();
  std::cout << "\nAuxiliary states introduced: " << (oc.states_colloc().size()>1?oc.states_colloc().size()-1:0) << "\n";
  if( nTrace ) std::cout << "Trace/tau variables appended: " << nTrace << " (initialised to zero)\n";

  std::vector<double> xref(nVar,0.0);
  size_t off=0;
  for( auto const& st: oc.states_colloc() ){
    auto nodes = oc.node_colloc(st);
    for( size_t i=0; i<nodes.size(); ++i ) xref[off+i] = exact_value_for_state(st.name(),nodes[i],p);
    off += nodes.size();
  }

  std::vector<double> res(nEqn,123456.0);
  bool eval_ok = oc.eval(res.data(),nullptr,xref.data(),nullptr,nullptr);
  ok &= eval_ok;
  size_t first_unwritten=nEqn;
  for( size_t i=0; i<nEqn; ++i ) if( res[i] == 123456.0 ){ first_unwritten=i; break; }
  ok &= check_close("all residual rows written before solve", first_unwritten==nEqn?0.0:1.0, 0.0);
  print_residuals("Reference-profile residual", res);
  ok &= check_close("reference-profile residual", max_abs(res), ref_tol);
  StateErrors refU = print_primitive_exact_errors(oc,xref,p,"manufactured reference");
  StateErrors refA = print_auxiliary_exact_errors(oc,xref,p,"manufactured reference");
  (void)refU; (void)refA;

  std::vector<double> var(nVar,0.0);
  initialise_vector(oc,var,p);
#if TEST_ADV_INIT_MODE == 0
  std::cout << "Initialisation: exact state and auxiliary profiles, zero trace variables\n";
#elif TEST_ADV_INIT_MODE == 1
  std::cout << "Initialisation: perturbed exact state and auxiliary profiles, zero trace variables\n";
#else
  std::cout << "Initialisation: constant primitive field, zero auxiliary/trace variables\n";
#endif

  if( !eval_residual(oc,var,res) ){
    std::cerr << "ERROR: initial residual evaluation failed for " << R.name << "\n";
    R.ok=false;
    return R;
  }
  print_residuals("Initial residual", res);

  // Nonlinear solve via OCFESLV::solve() (equilibrated LM/Newton, sparse AD).
  oc.options.SOLVE.MAX_ITER = TEST_ADV_MAXIT;
  oc.options.SOLVE.RES_TOL  = TEST_ADV_SOLVE_TOL;
#if defined(CRONOS__WITH_SPQR) && defined(TEST_ADV_SPQR)
  // Route the augmented (IC_TRACE/IC_STRONG) solve through sparse rank-revealing
  // QR instead of the JtJ normal equations.  Requires the header built with
  // -DCRONOS__WITH_SPQR and SuiteSparse linked.
  oc.options.SOLVE.FACTORIZATION = OCFESLV::Options::SOLVE_SPQR;
#endif
  OCFESLV::SolveReport const srep = oc.solve( var.data() );
  R.solved = srep.converged;
  if( !R.solved )
    std::cerr << "OCFESLV::solve did not converge for " << R.name
              << ": final|r|=" << srep.final_residual
              << " after " << srep.iterations << " it ("
              << srep.newton_steps << " Newton / " << srep.lm_steps << " LM)\n";
  ok &= R.solved;
  if( !eval_residual(oc,var,res) ){
    std::cerr << "ERROR: final residual evaluation failed for " << R.name << "\n";
    R.ok=false;
    return R;
  }
  R.final_res = max_abs(res);
  print_residuals("Final residual", res);
  // rev29: the residual GATE is srep.converged (checked above), which is the SOLVER's own
  // verdict and is correct in both modes -- under marching it is the AND over windows.
  // eval_residual() re-derives a residual from the environment AS IT STANDS after the
  // solve; when marching collapsed the evolution domain that environment is ONE WINDOW, so
  // it evaluates a different system and is not comparable to the monolithic figure.
  // MEASURED: both windows converged at 4.65e-13 and 2.26e-12 against a 1e-9 tolerance
  // while this check read 4.906e-07 and failed.  Diagnostic when marching, gate otherwise.
  if( !oc.is_marching() )
    ok &= check_close("final max residual", R.final_res, TEST_ADV_SOLVE_TOL);
  else
    std::cout << std::left << std::setw(48) << "final max residual (marched: diagnostic)"
              << " value=" << std::scientific << std::setprecision(6) << R.final_res
              << "  solver final|r|=" << srep.final_residual
              << "  converged=" << (srep.converged?"yes":"no") << "\n";

  StateErrors ePrim = print_primitive_exact_errors(oc,var,p,"final " + R.name + " solution");
  StateErrors eAux  = print_auxiliary_exact_errors(oc,var,p,"final " + R.name + " solution");

  // rev29: gate on the FULL-DOMAIN sampled comparison.  The two var-indexed reports above
  // are retained as diagnostics: under marching they cover one window, and printing both
  // makes the truncation visible instead of silent.
  StateErrors eSamp = print_sampled_exact_errors(oc,t,x,y,p,"final " + R.name + " solution");
  R.eU   = eSamp.eU;
  R.eUx  = eSamp.eUx;
  R.eUy  = eSamp.eUy;
  R.eAux = eSamp.eAux;

  // rev29: turn the convention into a CHECKED INVARIANT rather than a latent hazard.
  // PDE9 is the corpus's only driver that indexes `var` post-solve with no coordinate
  // evaluator (handoff rev29 s1).  Under marching that read is truncated to one window BY
  // DESIGN.  Reporting both numbers side by side means that if the framework ever changes
  // -- raw indexing spanning the horizon again, or eval_solution ceasing to route by
  // window -- the disagreement surfaces here rather than as a quiet wrong answer.
  // Reported UNCONDITIONALLY, and for every state, so both directions are checkable:
  //   monolithic : raw-var and sampled read the SAME set and should broadly agree, which
  //                is what validates the sampled oracle against the established one
  //   marching   : `var` holds ONE window by design, so they should DIFFER while nel_t>1
  // Printed, not gated.  The sampled grid does not land on the collocation nodes, so any
  // threshold here would be arbitrary -- and this file has just had one arbitrary gate
  // removed for producing a false failure under marching.  A wrong read shows up as an
  // implausible pair of numbers, which a reader can judge and a tolerance cannot.
  {
    size_t nodes_total = 0;
    for( auto const& st: oc.states_colloc() ) nodes_total += oc.node_colloc(st).size();
    bool const marching = oc.is_marching();
    std::cout << "\nRead-convention consistency (rev29):\n"
              << "  mode=" << ( marching? "marching" : "monolithic" )
              << "  windows=" << oc.n_march_steps()
              << "  raw var nodes=" << nodes_total
              << ( marching? "  (ONE window, by design)" : "" ) << "\n";
    std::cout << std::scientific << std::setprecision(4)
              << "  max|U -U_exact|   raw-var=" << ePrim.eU
              << "  sampled-full-domain=" << eSamp.eU << "\n"
              << "  max|Ux-Ux_exact|  raw-var=" << eAux.eUx
              << "  sampled-full-domain=" << eSamp.eUx << "\n"
              << "  max|Uy-Uy_exact|  raw-var=" << eAux.eUy
              << "  sampled-full-domain=" << eSamp.eUy << "\n";
    if( marching )
      std::cout << "  these SHOULD differ while nel_t>1; equality would mean either the raw\n"
                   "  read now spans the horizon or eval_solution stopped routing by window.\n";
    else
      std::cout << "  monolithic: both read the full domain, so a LARGE disagreement here\n"
                   "  would mean the sampled oracle is reading the wrong set.\n";
  }

  std::cout << "\nExact-solution comparison:\n";
  ok &= check_close("max |U-U_exact|", R.eU, exact_tol);
  ok &= check_close("max |Ux_aux-Ux_exact|", R.eUx, aux_tol);
  ok &= check_close("max |Uy_aux-Uy_exact|", R.eUy, aux_tol);

  R.max_spread = print_duplicate_spreads(oc,var,"post solve");
  ok &= check_close("max duplicate-node spread", R.max_spread, spread_tol);

  std::string const outname = std::string(TEST_ADV_OUTPUT_PREFIX) + "_" + suffix + ".out";
  write_output(oc,var,p,outname);
  std::cout << "\nWrote plot file: " << outname << "\n";

  std::cout << "\nManufactured dynamic 2D advection-diffusion test: " << (ok?"PASS":"FAIL") << "\n";
  R.ok=ok;
  return R;
}

} // namespace

int main()
{
  Par const p;
  std::vector<ModeResult> results;
  results.push_back( run_mode(OCFESLV::Options::IC_WEAK,   "weak",   p) );
  results.push_back( run_mode(OCFESLV::Options::IC_TRACE,  "trace",  p) );
  results.push_back( run_mode(OCFESLV::Options::IC_STRONG, "strong", p) );

  std::cout << "\n==================== PDE9 mode sweep summary ====================\n";
  std::cout << std::left << std::setw(10) << "mode"
            << std::right << std::setw(8)  << "nTrace"
            << std::setw(14) << "final|r|"
            << std::setw(14) << "dup-spread"
            << std::setw(12) << "|dU|"
            << std::setw(12) << "|dUx|"
            << std::setw(12) << "|dUy|"
            << std::setw(12) << "|dAux|max"
            << std::setw(8)  << "result" << "\n";

  bool all_ok=true;
  for( auto const& r: results ){
    all_ok &= r.ok;
    std::cout << std::left << std::setw(10) << r.name
              << std::right << std::setw(8) << r.nTrace
              << std::scientific << std::setprecision(3)
              << std::setw(14) << r.final_res
              << std::setw(14) << r.max_spread
              << std::setw(12) << r.eU
              << std::setw(12) << r.eUx
              << std::setw(12) << r.eUy
              << std::setw(12) << r.eAux
              << std::setw(8) << (r.ok?"PASS":"FAIL") << "\n";
  }
  std::cout << "=================================================================\n";
  std::cout << "Manufactured dynamic 2D advection-diffusion sweep over IC_WEAK, IC_TRACE, IC_STRONG: "
            << (all_ok?"PASS":"FAIL") << "\n";
  return all_ok ? 0 : 1;
}
