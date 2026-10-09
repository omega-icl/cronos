// Copyright (C) Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the EPL-2.0 with GPL-2.0-or-later as a Secondary License; see the LICENSE file.
// SPDX-License-Identifier: EPL-2.0 OR GPL-2.0-or-later
//
// OCFE_PDE5_solve.cpp
// -------------------
// Manufactured dynamic 2D Laplace/heat-equation regression test for OCFESLV.
//
// The model uses three collocation domains, time t and space x,y:
//
//   u_t - kappa * ( u_xx + u_yy ) = f(t,x,y),       (t,x,y) in (0,T]x(0,1)^2
//
// with manufactured initial and Dirichlet boundary data.  The exact solution is
// a smooth central Gaussian-like bump whose variance grows in time, so the
// output snapshots show a source-like profile spreading from the centre of the
// square.  The source and the reduced-order auxiliary derivatives remain
// available analytically.  This keeps the test compact while exercising
// finite-element interfaces in all three domains and the reduced first-order
// LINK equations introduced for the spatial Laplacian.
//
// The test rebuilds and solves a fresh OCFESLV for each imposition mode:
//
//   IC_WEAK, IC_TRACE, IC_STRONG
//
// and prints per-mode residuals, exact-profile errors, duplicate-node spreads,
// output files, and a final summary table.
//
// Build-time switches:
//   -OCFE_OCFESLV_HEADER          header to include; defaults to "ocfeslv.hpp"
//   -DTEST_LAP_NEL_T            finite elements in t; defaults to 3
//   -DTEST_LAP_NEL_X            finite elements in x; defaults to 3
//   -DTEST_LAP_NEL_Y            finite elements in y; defaults to 3
//   -DTEST_LAP_NT               collocation nodes/element in t; defaults to 3
//   -DTEST_LAP_NX               collocation nodes/element in x; defaults to 3
//   -DTEST_LAP_NY               collocation nodes/element in y; defaults to 3
//   -DTEST_LAP_INIT_MODE        0 exact, 1 perturbed exact, 2 constant; defaults to 1
//   -DTEST_LAP_MAXIT            Newton/LM iterations; defaults to 20
//   -DTEST_LAP_SOLVE_TOL        residual tolerance; defaults to 1e-9
//   -DTEST_LAP_EXACT_TOL        exact-profile tol (TRACE/STRONG); defaults to 1e-5
//   -DTEST_LAP_AUX_TOL          aux derivative tol (TRACE/STRONG); defaults to 2e-4
//   -DTEST_LAP_REF_TOL          reference-profile residual tol; defaults to 5e-4
//   -DTEST_LAP_SPREAD_TOL       duplicate spread tol (TRACE/STRONG); defaults to 5e-8
//   -DTEST_LAP_EXACT_TOL_WEAK   exact-profile tol (IC_WEAK, penalty); defaults to 5e-3
//   -DTEST_LAP_AUX_TOL_WEAK     aux derivative tol (IC_WEAK, penalty); defaults to 2e-1
//   -DTEST_LAP_SPREAD_TOL_WEAK  duplicate spread tol (IC_WEAK, penalty); defaults to 1e-2
//   -DTEST_LAP_OUTPUT_PREFIX    output file prefix; defaults to "OCFE_PDE5_solve"
//   -DTEST_LAP_PRINT_OC         print OCFESLV structure; defaults to 0

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

//#ifndef MC__OCFESLV_INTERFACE_DECISION_V2
//#define MC__OCFESLV_INTERFACE_DECISION_V2
//#endif

#ifndef OCFE_OCFESLV_HEADER
#define OCFE_OCFESLV_HEADER "ocfeslv.hpp"
#endif

//#define MC__OCFESLV_INTERFACE_DECISION_LEGACY
//#define MC__OCFESLV_STRONG_IMPLICIT_TAU

#include OCFE_OCFESLV_HEADER

using namespace mc;

#ifndef TEST_LAP_NEL_T
#define TEST_LAP_NEL_T 2
#endif
#ifndef TEST_LAP_NEL_X
#define TEST_LAP_NEL_X 2
#endif
#ifndef TEST_LAP_NEL_Y
#define TEST_LAP_NEL_Y 2
#endif
#ifndef TEST_LAP_NT
#define TEST_LAP_NT 10
#endif
#ifndef TEST_LAP_NX
#define TEST_LAP_NX 10
#endif
#ifndef TEST_LAP_NY
#define TEST_LAP_NY 10
#endif
#ifndef TEST_LAP_INIT_MODE
#define TEST_LAP_INIT_MODE 1
#endif
#ifndef TEST_LAP_MAXIT
#define TEST_LAP_MAXIT 20
#endif
#ifndef TEST_LAP_SOLVE_TOL
#define TEST_LAP_SOLVE_TOL 1e-9
#endif
// Accuracy tolerances are mode-aware (selected in run_mode).  The exact-multiplier
// modes (IC_TRACE/IC_STRONG) enforce continuity exactly and reach spectral accuracy
// for the gentled bump (width0=0.1) at NEL=2/N=10: |U-exact|~2e-6, |aux-exact|~4e-5,
// ref-residual~7e-5, dup-spread~3e-16.  The bars below sit ~4-7x above those.  IC_WEAK
// is penalty-limited (its SAT consistency error is O(1/sigma0), here ~6e-2 on the aux),
// so it gets penalty-level bars; its accuracy is not a correctness invariant.  The
// corner-fix regression guard is the EXACT-mode dup-spread (kept tight at 5e-8): a
// regression would push it from 3e-16 back toward the 1.577e-1 float and fail by ~7
// orders.  The WEAK dup-spread bar (1e-2) is a looser secondary backstop.
#ifndef TEST_LAP_EXACT_TOL
#define TEST_LAP_EXACT_TOL 1e-5
#endif
#ifndef TEST_LAP_AUX_TOL
#define TEST_LAP_AUX_TOL 2e-4
#endif
#ifndef TEST_LAP_DERIV_TOL
// eval_colloc_deriv spectrally differentiates the collocated U -- inherently ~one order coarser
// than U itself (and dU/dt on the wide window is the limiter), so it needs its own bar, distinct
// from the SOLVED-aux tol.  ~2.8e-4 observed at NT=10; 1e-3 leaves headroom yet a wrong evolution
// rescale (~50% error, O(0.1-1)) still fails loudly.  Mode-independent (interior spectral deriv).
#define TEST_LAP_DERIV_TOL 1e-3
#endif
#ifndef TEST_LAP_REF_TOL
#define TEST_LAP_REF_TOL 5e-4
#endif
#ifndef TEST_LAP_SPREAD_TOL
#define TEST_LAP_SPREAD_TOL 5e-8
#endif
#ifndef TEST_LAP_EXACT_TOL_WEAK
#define TEST_LAP_EXACT_TOL_WEAK 5e-3
#endif
#ifndef TEST_LAP_AUX_TOL_WEAK
#define TEST_LAP_AUX_TOL_WEAK 2e-1
#endif
#ifndef TEST_LAP_SPREAD_TOL_WEAK
#define TEST_LAP_SPREAD_TOL_WEAK 1e-2
#endif
#ifndef TEST_LAP_SAT_SIGMA0
#define TEST_LAP_SAT_SIGMA0 100.0
#endif
#ifndef TEST_LAP_PRINT_OC
#define TEST_LAP_PRINT_OC 0
#endif
#ifndef TEST_LAP_OUTPUT_PREFIX
#define TEST_LAP_OUTPUT_PREFIX "OCFE_PDE5_solve"
#endif

#define TEST_LAP_SPQR

namespace {

struct Par
{
  double tf      = 1.0;
  double kappa   = 3.0e-2;
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
static double F_exact( double t, double x, double y, Par const& p )
{ return U_t_exact(t,x,y,p) - p.kappa*( U_xx_exact(t,x,y,p) + U_yy_exact(t,x,y,p) ); }

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
#if TEST_LAP_INIT_MODE == 1
      val += 1.0e-2 * std::sin(0.31*double(off+i+1));
#elif TEST_LAP_INIT_MODE == 2
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
  os << "# OCFE_PDE5_solve time-sliced spreading-bump output\n";
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
  double parity=0.0;    // max |eval_colloc - eval_solution| over the interior field (both modes)
  double eU_field=0.0;  // full (t,x,y) field: max |eval_colloc - exact| via the SAME call path in both modes
};

// Max error of eval_colloc_deriv(U, wrt=t/x/y, ORD=1) vs the manufactured exact derivatives,
// sampled at interior (t,x,y).  Under marching this routes through the trajectory (and the dU/dt
// read exercises the evolution block-width rescale).
static double deriv_field_err( OCFESLV& oc, FFVar const& U, FFVar const& t, FFVar const& x,
                               FFVar const& y, std::vector<double> const& var, Par const& p,
                               char const* tag )
{
  double et=0., ex=0., ey=0.;
  for( int it=1; it<=5; ++it ){
    double const tv = double(it)/6.0*p.tf;   // it=1,2 -> 40% window; it=3,4,5 -> 60% window
    for( double xv : {0.25,0.5,0.75} )
      for( double yv : {0.25,0.5,0.75} ){
        OCFESLV::t_Coord pt; pt[t]=tv; pt[x]=xv; pt[y]=yv;
        double const ut = oc.eval_colloc_deriv<double>( U, pt, t, 1, var.data(), nullptr, nullptr );
        double const ux = oc.eval_colloc_deriv<double>( U, pt, x, 1, var.data(), nullptr, nullptr );
        double const uy = oc.eval_colloc_deriv<double>( U, pt, y, 1, var.data(), nullptr, nullptr );
        et = std::max( et, std::fabs( ut - U_t_exact(tv,xv,yv,p) ) );
        ex = std::max( ex, std::fabs( ux - U_x_exact(tv,xv,yv,p) ) );
        ey = std::max( ey, std::fabs( uy - U_y_exact(tv,xv,yv,p) ) );
      }
  }
  std::cout << "eval_colloc_deriv vs exact (" << tag << "): dU/dt=" << et
            << " dU/dx=" << ex << " dU/dy=" << ey << "  (max=" << std::max({et,ex,ey}) << ")\n";
  return std::max( { et, ex, ey } );
}

static ModeResult run_mode( OCFESLV::Options::ImpositionType imp,
                            std::string const& suffix,
                            Par const& p,
                            bool marching = false )
{
  ModeResult R;
  R.name = imp_name(imp) + std::string( marching ? " [march]" : " [mono]" );
  bool ok=true;

  // Mode-aware accuracy bars: IC_WEAK is penalty-limited, the exact-multiplier
  // modes reach spectral accuracy.  ref-residual and the solve residual are
  // mode-independent and stay uniform.
  bool   const weak       = ( imp == OCFESLV::Options::IC_WEAK );
  double const exact_tol  = weak ? TEST_LAP_EXACT_TOL_WEAK  : TEST_LAP_EXACT_TOL;
  double const aux_tol    = weak ? TEST_LAP_AUX_TOL_WEAK    : TEST_LAP_AUX_TOL;
  double const deriv_tol  = TEST_LAP_DERIV_TOL;   // spectral-derivative bar for eval_colloc_deriv
  double const spread_tol = weak ? TEST_LAP_SPREAD_TOL_WEAK : TEST_LAP_SPREAD_TOL;
  double const ref_tol    = TEST_LAP_REF_TOL;

  std::cout << "\n========== manufactured dynamic 2D Laplace/heat test =========="
            << "\nimposition: " << R.name
            << ", finite elements: t=" << TEST_LAP_NEL_T
            << " x=" << TEST_LAP_NEL_X
            << " y=" << TEST_LAP_NEL_Y
            << ", nodes/element: t=" << TEST_LAP_NT
            << " x=" << TEST_LAP_NX
            << " y=" << TEST_LAP_NY
            << ", init_mode=" << TEST_LAP_INIT_MODE << "\n";
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
           - p.kappa * BE * ( 4.0*RR/(S*S) - 4.0/S );

  FFVar PDE = OpP(U,t) - p.kappa*( OpP(U,{x,2}) + OpP(U,{y,2}) ) - FE;
  FFVar BC  = U - UE;

  OCFESLV oc(&DAG);
  // Evolution grid: for the default 2-element case, split the two marching windows 40%/60% (a
  // NON-UNIFORM grid) so the eval_colloc_deriv block-width rescale is exercised -- the dU/dt read
  // in the narrow first window must use that window's width, not the resident (last) window's.
  oc.add_domain( t, ( TEST_LAP_NEL_T == 2 )
    ? FFDom( std::vector<double>{ 0., 0.4*p.tf, p.tf }, FFDom::LGR, TEST_LAP_NT )
    : FFDom( 0., p.tf, TEST_LAP_NEL_T, FFDom::LGR, TEST_LAP_NT ) );
  oc.add_domain( x, FFDom(0., 1.0,  TEST_LAP_NEL_X, FFDom::CGL, TEST_LAP_NX) );
  oc.add_domain( y, FFDom(0., 1.0,  TEST_LAP_NEL_Y, FFDom::CGL, TEST_LAP_NY) );
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

  oc.set_evolution_domain( t );
  oc.options.SOLVE.MARCHING  = marching;   // both modes exercised: monolithic (full-domain solve) and
                                           // marching over the DISTRIBUTED state U(t,x,y).  Marching takes
                                           // the general-IC override path, now generalised to a distributed
                                           // initial condition (one continuity row per spatial LB node).
  //oc.reset_evolution_domain();
  oc.options.REDUCE.ORDER    = OCFESLV::Options::RED_FULL;
  oc.options.CLASSIFY.MODE        = OCFESLV::Options::CLASS_AUTO;
  oc.options.INTERFACE.TYPE  = OCFESLV::Options::IC_AUTO;
  oc.options.INTERFACE.IMPOSITION = imp;
  oc.options.INTERFACE.SAT_SIGMA0      = TEST_LAP_SAT_SIGMA0;
  oc.options.DISPLAY_LEVEL   = 2;
  oc.options.SOLVE.WARMSTART = OCFESLV::Options::REUSE;//BROADCAST_IC;
  //oc.options.MARCH_IC_C0     = true;

  if( !oc.setup() ){
    std::cerr << "ERROR: OCFESLV::setup() failed for " << R.name << "\n";
    R.ok=false;
    return R;
  }
#if TEST_LAP_PRINT_OC
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
#if TEST_LAP_INIT_MODE == 0
  std::cout << "Initialisation: exact state and auxiliary profiles, zero trace variables\n";
#elif TEST_LAP_INIT_MODE == 1
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

  // item 12: nonlinear solve via OCFESLV::solve() (equilibrated sparse damped LM).
  oc.options.SOLVE.MAX_ITER = TEST_LAP_MAXIT;
  oc.options.SOLVE.RES_TOL  = TEST_LAP_SOLVE_TOL;
#if defined(CRONOS__WITH_SPQR) && defined(TEST_LAP_SPQR)
  // item 13: route the augmented (IC_TRACE/IC_STRONG) solve through sparse
  // rank-revealing QR instead of the JtJ normal equations.  Requires the header
  // built with -DCRONOS__WITH_SPQR and SuiteSparse linked.
  oc.options.SOLVE.FACTORIZATION = OCFESLV::Options::SOLVE_SPQR;
#endif
  OCFESLV::SolveReport const srep = oc.solve( var.data() );
  R.solved = srep.converged;
  if( !R.solved )
    std::cerr << "OCFESLV::solve did not converge: final|r|=" << srep.final_residual
              << " after " << srep.iterations << " it\n";
  ok &= R.solved;
  if( !marching ){
    if( !eval_residual(oc,var,res) ){
      std::cerr << "ERROR: final residual evaluation failed for " << R.name << "\n";
      R.ok=false;
      return R;
    }
    R.final_res = max_abs(res);
    print_residuals("Final residual", res);
    ok &= check_close("final max residual", R.final_res, TEST_LAP_SOLVE_TOL);
  }
  else {
    // A full-domain residual re-eval is NOT meaningful under marching: the marched solution
    // satisfies inter-window continuity at each window's IC (var[col]-terminal), not the
    // as-written IC (U-UE), so the IC rows legitimately carry a nonzero as-written residual.
    // Per-window convergence (srep) is the residual criterion under marching.
    R.final_res = srep.final_residual;
    std::cout << "Marching per-window final residual (max over windows) = " << R.final_res << "\n";
  }

  if( !marching ){
    // Monolithic: the whole space-time solution lives in var, so compare every collocation
    // node against the manufactured exact and check duplicate-node continuity.
    StateErrors ePrim = print_primitive_exact_errors(oc,var,p,"final " + R.name + " solution");
    StateErrors eAux  = print_auxiliary_exact_errors(oc,var,p,"final " + R.name + " solution");
    R.eU = ePrim.eU;
    R.eUx = eAux.eUx;
    R.eUy = eAux.eUy;
    R.eAux = eAux.eAux;

    std::cout << "\nExact-solution comparison:\n";
    ok &= check_close("max |U-U_exact|", R.eU, exact_tol);
    ok &= check_close("max |Ux_aux-Ux_exact|", R.eUx, aux_tol);
    ok &= check_close("max |Uy_aux-Uy_exact|", R.eUy, aux_tol);

    R.max_spread = print_duplicate_spreads(oc,var,"post solve");
    ok &= check_close("max duplicate-node spread", R.max_spread, spread_tol);

    // Read-consolidation: the SAME plain eval_colloc(state,{t,x,y},var) reads the full field
    // regardless of solve mode.  Monolithic baseline -- sample interior (t,x,y) vs exact.
    {
      double emax=0.0, dmax=0.0;
      for( int it=0; it<=6; ++it ){
        double const tv = double(it)/6.0*p.tf;
        for( double xv : {0.25,0.5,0.75} )
          for( double yv : {0.25,0.5,0.75} ){
            OCFESLV::t_Coord pt; pt[t]=tv; pt[x]=xv; pt[y]=yv;
            double const got = oc.eval_colloc<double>(U, pt, var.data(), nullptr, nullptr);
            double const sol = oc.eval_solution(U, pt);   // buffer-free; reads the cached primal
            emax = std::max( emax, std::abs( got - U_exact(tv,xv,yv,p) ) );
            dmax = std::max( dmax, std::abs( got - sol ) );
          }
      }
      R.eU_field = emax; R.parity = dmax;
      std::cout << "Full-field (monolithic): eval_colloc vs exact=" << emax
                << "  eval_colloc vs eval_solution=" << dmax << "\n";
      ok &= check_close("monolithic full-field eval_colloc vs exact", R.eU_field, exact_tol);
      ok &= check_close("monolithic eval_colloc == eval_solution parity", R.parity, 1e-9);
      double const dderr = deriv_field_err( oc, U, t, x, y, var, p, "monolithic" );
      ok &= check_close("monolithic eval_colloc_deriv vs exact", dderr, deriv_tol);
    }

    std::string const outname = std::string(TEST_LAP_OUTPUT_PREFIX) + "_" + suffix + ".out";
    write_output(oc,var,p,outname);
    std::cout << "\nWrote plot file: " << outname << "\n";
  }
  else {
    // Marching: var holds the final window's collapsed template.  Validate the marched
    // trajectory through its terminal profile at t=tf -- terminal_profile()/terminal_value()
    // interpolate the state at the evolution UB, so a correct distributed-IC override (each
    // spatial LB node continued from the previous window's terminal) reproduces the exact
    // manufactured bump at the final time across the full spatial grid.
    auto const prof = oc.terminal_profile( var.data() );

    // Spatial-node positions in the same (x fastest, y slower) order terminal_profile uses.
    // Reconstruct the x/y domains from the same parameters passed to add_domain() above; the
    // reference nodes are computed lazily, so set_nodes() must be called before lgnodes().  Shared
    // by the U, auxiliary-derivative and duplicate-spread terminal checks below.
    FFDom dx(0., 1.0, TEST_LAP_NEL_X, FFDom::CGL, TEST_LAP_NX); dx.set_nodes();
    FFDom dy(0., 1.0, TEST_LAP_NEL_Y, FFDom::CGL, TEST_LAP_NY); dy.set_nodes();
    std::vector<double> xs, ys;
    for( size_t ie=0; ie<dx.n_elem; ++ie ){ auto e=dx.lgnodes(dx.elem_lo(ie),dx.elem_up(ie)); xs.insert(xs.end(),e.begin(),e.end()); }
    for( size_t ie=0; ie<dy.n_elem; ++ie ){ auto e=dy.lgnodes(dy.elem_lo(ie),dy.elem_up(ie)); ys.insert(ys.end(),e.begin(),e.end()); }

    auto const itU  = prof.find( U );
    double emax_prof = 0.0;
    if( itU != prof.end() && !itU->second.empty() ){
      auto const& pu = itU->second;
      double emax_prof_val = 0.0;   // cross-check: same nodes via terminal_value() interpolation
      size_t kmax = 0; double xmax=0, ymax=0, pumax=0, uemax=0;
      for( size_t iy=0, k=0; iy<ys.size(); ++iy )
        for( size_t ix=0; ix<xs.size(); ++ix, ++k )
          if( k < pu.size() ){
            double const ue = U_exact(p.tf, xs[ix], ys[iy], p);
            double const e  = std::abs( pu[k] - ue );
            if( e > emax_prof ){ emax_prof=e; kmax=k; xmax=xs[ix]; ymax=ys[iy]; pumax=pu[k]; uemax=ue; }
            OCFESLV::t_Coord pt; pt[x]=xs[ix]; pt[y]=ys[iy];
            emax_prof_val = std::max( emax_prof_val, std::abs( oc.terminal_value(U,pt,var.data()) - ue ) );
          }
      std::cout << "\nMarched terminal profile: nodes=" << pu.size()
                << " grid=" << xs.size() << "x" << ys.size()
                << " max|U(tf)-U_exact(tf)|=" << emax_prof
                << "  (at node k=" << kmax << " (x,y)=(" << xmax << "," << ymax
                << ") U=" << pumax << " exact=" << uemax << ")\n";
      std::cout << "  cross-check same nodes via terminal_value(): max err=" << emax_prof_val << "\n";
    }
    else {
      std::cout << "\nERROR: terminal_profile() returned no U profile under marching\n";
      ok = false;
    }
    // Independent spot-check via terminal_value() interpolation at interior points, printed per point.
    double emax_pt = 0.0;
    std::cout << "Marched terminal spot-check (interior points):\n";
    for( double xv : {0.25,0.5,0.75} )
      for( double yv : {0.25,0.5,0.75} ){
        OCFESLV::t_Coord pt; pt[x]=xv; pt[y]=yv;
        double const got = oc.terminal_value(U,pt,var.data());
        double const ue  = U_exact(p.tf,xv,yv,p);
        emax_pt = std::max( emax_pt, std::abs( got - ue ) );
        std::cout << "   (" << xv << "," << yv << ")  U=" << got << "  exact=" << ue
                  << "  |err|=" << std::abs(got-ue) << "\n";
      }
    std::cout << "spot-check max|U(tf)-U_exact(tf)|=" << emax_pt << "\n";

    R.eU = std::max( emax_prof, emax_pt );

    // Auxiliary-derivative accuracy at t=tf from the terminal profiles of the aux states
    // (Dx_U -> U_x_exact, Dy_U -> U_y_exact), mirroring the monolithic aux comparison so the
    // summary reports real terminal errors rather than placeholder zeros.
    double eUx_t=0.0, eUy_t=0.0;
    for( auto const& kv : prof ){
      std::string const nm = kv.first.name();
      if( !is_auxiliary_state_name(nm) ) continue;
      bool const isx = nm.find("Dx_") != std::string::npos;
      auto const& pv = kv.second;
      double emax = 0.0;
      for( size_t iy=0,k=0; iy<ys.size(); ++iy )
        for( size_t ix=0; ix<xs.size(); ++ix,++k )
          if( k < pv.size() ){
            double const ref = isx ? U_x_exact(p.tf,xs[ix],ys[iy],p)
                                   : U_y_exact(p.tf,xs[ix],ys[iy],p);
            emax = std::max( emax, std::abs( pv[k]-ref ) );
          }
      if( isx ) eUx_t = std::max(eUx_t,emax); else eUy_t = std::max(eUy_t,emax);
    }
    R.eUx = eUx_t; R.eUy = eUy_t; R.eAux = std::max(eUx_t,eUy_t);

    // Interface duplicate-node spread of the marched solution: computed from the raw collocation
    // DOFs of the final window (var), so it reflects genuine element-interface continuity -- unlike
    // terminal_profile(), which interpolates one value per physical node and collapses duplicates.
    R.max_spread = print_duplicate_spreads(oc,var,"marched final window");

    std::cout << "\nExact-solution comparison (marched terminal at t=tf):\n";
    ok &= check_close("marched terminal max |U-U_exact|",     R.eU,        exact_tol);
    ok &= check_close("marched terminal max |Ux_aux-Ux_exact|", R.eUx,      aux_tol);
    ok &= check_close("marched terminal max |Uy_aux-Uy_exact|", R.eUy,      aux_tol);
    ok &= check_close("marched final-window duplicate spread",  R.max_spread, spread_tol);

    // --- Read-consolidation parity (Stage B) ---------------------------------------------------
    // After a march the SAME plain eval_colloc(state,{t,x,y},var) over the WHOLE t in [0,tf] must:
    //   (i) succeed at INTERIOR t (previously it aborted -- var/domain were the last window only);
    //  (ii) match the buffer-free eval_solution() to ~machine precision (both read the trajectory);
    // (iii) hit the exact solution to the marched tolerance across the full field, not just t=tf.
    {
      double dmax=0.0, emax=0.0;
      for( int it=0; it<=6; ++it ){
        double const tv = double(it)/6.0*p.tf;   // spans both windows incl. the interface t=tf/2
        for( double xv : {0.25,0.5,0.75} )
          for( double yv : {0.25,0.5,0.75} ){
            OCFESLV::t_Coord pt; pt[t]=tv; pt[x]=xv; pt[y]=yv;
            double const via_colloc   = oc.eval_colloc<double>(U, pt, var.data(), nullptr, nullptr);
            double const via_solution = oc.eval_solution(U, pt);   // buffer-free; also routes the trajectory
            dmax = std::max( dmax, std::abs( via_colloc - via_solution ) );
            emax = std::max( emax, std::abs( via_colloc - U_exact(tv,xv,yv,p) ) );
          }
      }
      R.parity = dmax; R.eU_field = emax;
      std::cout << "\nRead-consolidation (marching): full-field eval_colloc over t in [0," << p.tf << "]\n";
      std::cout << "  eval_colloc vs eval_solution max diff = " << dmax << "\n";
      std::cout << "  eval_colloc vs exact         max err  = " << emax << "\n";
      ok &= check_close("eval_colloc == eval_solution parity",       R.parity,   1e-9);
      ok &= check_close("marching full-field eval_colloc vs exact",  R.eU_field, exact_tol);
      double const dderr = deriv_field_err( oc, U, t, x, y, var, p, "marching" );
      ok &= check_close("marching eval_colloc_deriv vs exact", dderr, deriv_tol);
    }
  }

  std::cout << "\nManufactured dynamic 2D Laplace/heat test: " << (ok?"PASS":"FAIL") << "\n";
  R.ok=ok;
  return R;
}

} // namespace

int main()
{
  Par const p;
  std::vector<ModeResult> results;
  results.push_back( run_mode(OCFESLV::Options::IC_WEAK,   "weak",   p, false) );
  results.push_back( run_mode(OCFESLV::Options::IC_WEAK,   "weak",   p, true ) );
  results.push_back( run_mode(OCFESLV::Options::IC_TRACE,  "trace",  p, false) );
  results.push_back( run_mode(OCFESLV::Options::IC_TRACE,  "trace",  p, true ) );
  results.push_back( run_mode(OCFESLV::Options::IC_STRONG, "strong", p, false) );
  results.push_back( run_mode(OCFESLV::Options::IC_STRONG, "strong", p, true ) );

  std::cout << "\n==================== PDE5 mode sweep summary ====================\n";
  std::cout << std::left << std::setw(10) << "mode"
            << std::right << std::setw(8)  << "nTrace"
            << std::setw(14) << "final|r|"
            << std::setw(14) << "dup-spread"
            << std::setw(12) << "|dU|"
            << std::setw(12) << "|dUx|"
            << std::setw(12) << "|dUy|"
            << std::setw(12) << "|dAux|max"
            << std::setw(12) << "rd-parity"
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
              << std::setw(12) << r.parity
              << std::setw(8) << (r.ok?"PASS":"FAIL") << "\n";
  }
  std::cout << "=================================================================\n";
  std::cout << "Manufactured dynamic 2D Laplace/heat sweep over IC_WEAK, IC_TRACE, IC_STRONG"
            << " x {monolithic, marching}: "
            << (all_ok?"PASS":"FAIL") << "\n";
  return all_ok ? 0 : 1;
}
