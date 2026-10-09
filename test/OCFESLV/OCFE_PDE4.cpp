// Copyright (C) Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the EPL-2.0 with GPL-2.0-or-later as a Secondary License; see the LICENSE file.
// SPDX-License-Identifier: EPL-2.0 OR GPL-2.0-or-later
//
// OCFE_PDE4_solve.cpp
// --------------------
// Manufactured-solution membrane-only regression test for OCFESLV.
//
// This is a reduced version of OCFE_PDE3_solve.cpp that keeps only the dry
// and wet membrane compartments:
//
//   dry membrane : Cd(rd,z)
//   wet membrane : Cw(rw,z), Sw(rw,z)
//
// The gas and liquid compartments are deliberately removed.  The dry-side
// gas/membrane boundary is replaced by a constant Dirichlet condition
//
//   Cd(rd=0,z) = Cdry0.
//
// The dry/wet interface keeps the same Henry value condition, CO2 flux
// continuity, and zero solvent flux used in PDE3.  The wet outer side is closed
// by manufactured Dirichlet values for Cw and Sw.  This gives a compact test
// aimed at isolating IC_TRACE behaviour at axial z element interfaces in the
// tensor-product dry/wet membrane blocks.
//
// Build-time switches:
//   -DTEST_MEM_IMPOSE        use IC_WEAK (0), IC_STRONG (1), or IC_TRACE (2)
//   -DTEST_MEM_NEL_Z         number of axial finite elements
//   -DTEST_MEM_NEL_R         number of radial finite elements per membrane
//   -DTEST_MEM_NZ            axial nodes per element
//   -DTEST_MEM_NR            radial nodes per element
//   -DTEST_MEM_INIT_MODE     0 exact, 1 perturbed exact, 2 constant primitives
//   -DTEST_MEM_OCENV_HEADER  header to include; defaults to "ocfeslv.hpp"

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

#ifndef MC__OCFESLV_INTERFACE_ASSIGN_DEBUG
#define MC__OCFESLV_INTERFACE_ASSIGN_DEBUG
#endif
#ifndef MC__OCFESLV_INTERFACE_DEBUG
#define MC__OCFESLV_INTERFACE_DEBUG
#endif

#ifndef OCFE_OCFESLV_HEADER
#define OCFE_OCFESLV_HEADER "ocfeslv.hpp"
#endif

//#define MC__OCFESLV_STRONG_EXPLICIT_TAU
//#define MC__OCFESLV_INTERFACE_DECISION_V2
//#define MC__OCFESLV_STRONG_DUMMY_DIAG
//#define MC__OCFESLV_SYMBOL_TAG_PROBE
//#define MC__OCFESLV_SYMBOL_DROP_GATE

//#define MC__OCFESLV_INTERFACE_DECISION_V2
//#define MC__OCFESLV_GEOM_ROLE_DIAG
//#define MC__OCFESLV_STRONG_DUMMY_DIAG
//#define MC__OCFESLV_STRONG_EXPLICIT_TAU
#include OCFE_OCFESLV_HEADER

using namespace mc;

#ifndef TEST_MEM_NEL_Z
#define TEST_MEM_NEL_Z 3
#endif
#ifndef TEST_MEM_NEL_R
#define TEST_MEM_NEL_R 3
#endif
#ifndef TEST_MEM_NZ
#define TEST_MEM_NZ 7
#endif
#ifndef TEST_MEM_NR
#define TEST_MEM_NR 7
#endif
#ifndef TEST_MEM_IMPOSE
#define TEST_MEM_IMPOSE 2
#endif
#ifndef TEST_MEM_REDUCE_ORDER
#define TEST_MEM_REDUCE_ORDER 1
#endif
#ifndef TEST_MEM_CLASSIFY
#define TEST_MEM_CLASSIFY 1
#endif
#ifndef TEST_MEM_OUTPUT_PREFIX
#define TEST_MEM_OUTPUT_PREFIX "OCFE_PDE4_solve"
#endif
#ifndef TEST_MEM_MAXIT
#define TEST_MEM_MAXIT 30
#endif
#ifndef TEST_MEM_SOLVE_TOL
#define TEST_MEM_SOLVE_TOL 1e-9
#endif
#ifndef TEST_MEM_SAT_SIGMA0
#define TEST_MEM_SAT_SIGMA0 1.0
#endif
#ifndef TEST_MEM_PRINT_OC
#define TEST_MEM_PRINT_OC 1
#endif
#ifndef TEST_MEM_INIT_MODE
#define TEST_MEM_INIT_MODE 2
#endif

#define TEST_LAP_SPQR

namespace {

struct Par
{
  double L      = 1.0;
  double H      = 0.55;
  double Dmd    = 2.5e-2;
  double DmwC   = 1.1e-2;
  double DmwS   = 6.0e-3;
  double krw    = 0.28;
  double nuS    = 1.0;

  // Constant dry-side value replacing the gas/dry coupling in PDE3.
  double Cdry0  = 1.0;

  // Manufactured dry/wet interface and wet-wall amplitudes.  The exact
  // solution retains nontrivial z-dependence so that axial element-interface
  // continuity is still exercised.
  double Acd0   = 0.16;
  double AcdZ   = 0.04;
  double qdw0   = -5.0e-3;
  double qdwZ   = -1.0e-3;
  double dCd0   = -0.10;
  double CwU0   = 0.18;
  double CwUZ   = 0.03;
  double CwD1Z  = 0.020;
  double Sw0    = 0.95;
  double SwZ    = -0.35;
  double SwD1Z  = 0.040;
};

static double sqr( double x ){ return x*x; }
static double Zs( double z ){ return 3.0*z*z - 2.0*z*z*z; }
static double Zs_z( double z ){ return 6.0*z - 6.0*z*z; }

static double h00( double r ){ return 2.0*r*r*r - 3.0*r*r + 1.0; }
static double h10( double r ){ return r*r*r - 2.0*r*r + r; }
static double h01( double r ){ return -2.0*r*r*r + 3.0*r*r; }
static double h11( double r ){ return r*r*r - r*r; }
static double h00r( double r ){ return 6.0*r*r - 6.0*r; }
static double h10r( double r ){ return 3.0*r*r - 4.0*r + 1.0; }
static double h01r( double r ){ return -6.0*r*r + 6.0*r; }
static double h11r( double r ){ return 3.0*r*r - 2.0*r; }
static double h00rr( double r ){ return 12.0*r - 6.0; }
static double h10rr( double r ){ return 6.0*r - 4.0; }
static double h01rr( double r ){ return -12.0*r + 6.0; }
static double h11rr( double r ){ return 6.0*r - 2.0; }

static double herm( double r, double f0, double d0, double f1, double d1 )
{ return h00(r)*f0 + h10(r)*d0 + h01(r)*f1 + h11(r)*d1; }
static double herm_r( double r, double f0, double d0, double f1, double d1 )
{ return h00r(r)*f0 + h10r(r)*d0 + h01r(r)*f1 + h11r(r)*d1; }
static double herm_rr( double r, double f0, double d0, double f1, double d1 )
{ return h00rr(r)*f0 + h10rr(r)*d0 + h01rr(r)*f1 + h11rr(r)*d1; }

template <typename T>
static T HERM( T const& r, T const& f0, T const& d0, T const& f1, T const& d1 )
{ return (2.0*r*r*r-3.0*r*r+1.0)*f0 + (r*r*r-2.0*r*r+r)*d0
       + (-2.0*r*r*r+3.0*r*r)*f1 + (r*r*r-r*r)*d1; }
template <typename T>
static T HERM_R( T const& r, T const& f0, T const& d0, T const& f1, T const& d1 )
{ return (6.0*r*r-6.0*r)*f0 + (3.0*r*r-4.0*r+1.0)*d0
       + (-6.0*r*r+6.0*r)*f1 + (3.0*r*r-2.0*r)*d1; }
template <typename T>
static T HERM_RR( T const& r, T const& f0, T const& d0, T const& f1, T const& d1 )
{ return (12.0*r-6.0)*f0 + (6.0*r-4.0)*d0
       + (-12.0*r+6.0)*f1 + (6.0*r-2.0)*d1; }

static double CdU_exact( double z, Par const& p ){ return p.Acd0 + p.AcdZ*Zs(z); }
static double qdw_exact( double z, Par const& p ){ return p.qdw0 + p.qdwZ*Zs(z); }
static double Cd_d0_exact( double z, Par const& p ){ return p.dCd0 * ( 0.5 + 0.5*Zs(z) ); }
static double Cd_exact( double r, double z, Par const& p )
{ return herm( r, p.Cdry0, Cd_d0_exact(z,p), CdU_exact(z,p), qdw_exact(z,p)/p.Dmd ); }
static double Cd_r_exact( double r, double z, Par const& p )
{ return herm_r( r, p.Cdry0, Cd_d0_exact(z,p), CdU_exact(z,p), qdw_exact(z,p)/p.Dmd ); }
static double Cd_rr_exact( double r, double z, Par const& p )
{ return herm_rr( r, p.Cdry0, Cd_d0_exact(z,p), CdU_exact(z,p), qdw_exact(z,p)/p.Dmd ); }

static double Cw0_exact( double z, Par const& p ){ return p.H*CdU_exact(z,p); }
static double Cw_d0_exact( double z, Par const& p ){ return qdw_exact(z,p)/p.DmwC; }
static double CwU_exact( double z, Par const& p ){ return p.CwU0 + p.CwUZ*Zs(z); }
static double Cw_d1_exact( double z, Par const& p ){ return p.CwD1Z*Zs(z); }
static double Cw_exact( double r, double z, Par const& p )
{ return herm( r, Cw0_exact(z,p), Cw_d0_exact(z,p), CwU_exact(z,p), Cw_d1_exact(z,p) ); }
static double Cw_r_exact( double r, double z, Par const& p )
{ return herm_r( r, Cw0_exact(z,p), Cw_d0_exact(z,p), CwU_exact(z,p), Cw_d1_exact(z,p) ); }
static double Cw_rr_exact( double r, double z, Par const& p )
{ return herm_rr( r, Cw0_exact(z,p), Cw_d0_exact(z,p), CwU_exact(z,p), Cw_d1_exact(z,p) ); }

static double Sw0_exact( double z, Par const& p ){ return p.Sw0 + 0.08*Zs(z); }
static double SwU_exact( double z, Par const& p ){ return p.Sw0 + p.SwZ*Zs(z); }
static double Sw_d1_exact( double z, Par const& p ){ return p.SwD1Z*Zs(z); }
static double Sw_exact( double r, double z, Par const& p )
{ return herm( r, Sw0_exact(z,p), 0.0, SwU_exact(z,p), Sw_d1_exact(z,p) ); }
static double Sw_r_exact( double r, double z, Par const& p )
{ return herm_r( r, Sw0_exact(z,p), 0.0, SwU_exact(z,p), Sw_d1_exact(z,p) ); }
static double Sw_rr_exact( double r, double z, Par const& p )
{ return herm_rr( r, Sw0_exact(z,p), 0.0, SwU_exact(z,p), Sw_d1_exact(z,p) ); }

static double dry_c_src_exact( double r, double z, Par const& p )
{ return -p.Dmd*Cd_rr_exact(r,z,p); }
static double wet_c_src_exact( double r, double z, Par const& p )
{ return -p.DmwC*Cw_rr_exact(r,z,p) + p.krw*Cw_exact(r,z,p)*Sw_exact(r,z,p); }
static double wet_s_src_exact( double r, double z, Par const& p )
{ return -p.DmwS*Sw_rr_exact(r,z,p) + p.nuS*p.krw*Cw_exact(r,z,p)*Sw_exact(r,z,p); }

static double max_abs( std::vector<double> const& r )
{ double m=0.; for( double v: r ) m = std::max(m,std::abs(v)); return m; }

static std::vector<double> physical_nodes( FFDom const& dom )
{
  std::vector<double> x;
  for( size_t ie=0; ie<dom.n_elem; ++ie ){
    auto xe = dom.lgnodes( dom.elem_lo(ie), dom.elem_up(ie) );
    x.insert( x.end(), xe.begin(), xe.end() );
  }
  return x;
}

static void print_residuals( std::string const& label, std::vector<double> const& r )
{
  double sum=0.; for( double v: r ) sum += std::abs(v);
  std::cout << std::left << std::setw(46) << label
            << " n=" << std::setw(5) << r.size()
            << " max|r|=" << std::scientific << std::setprecision(4) << max_abs(r)
            << " mean|r|=" << (r.empty()?0.:sum/r.size()) << "\n";
}

static void print_largest_residuals( std::string const& label, std::vector<double> const& r, size_t nprint=12 )
{
  if( r.empty() ) return;
  std::vector<size_t> idx(r.size());
  std::iota(idx.begin(),idx.end(),size_t(0));
  std::partial_sort(idx.begin(),idx.begin()+std::min(nprint,idx.size()),idx.end(),
    [&](size_t a,size_t b){return std::abs(r[a])>std::abs(r[b]);});
  std::cout << label << " largest residual rows:\n";
  for( size_t k=0; k<std::min(nprint,idx.size()); ++k ){
    size_t const i=idx[k];
    std::cout << "  row " << std::setw(6) << i << "  r=" << std::scientific << std::setprecision(8) << r[i] << "\n";
  }
}

static bool check_close( std::string const& label, double val, double tol )
{
  bool ok = std::isfinite(val) && val <= tol;
  std::cout << std::left << std::setw(54) << label
            << " value=" << std::scientific << std::setprecision(6) << val
            << " tol=" << tol << "  " << (ok?"PASS":"FAIL") << "\n";
  return ok;
}

static double interp( OCFESLV const& oc, FFVar const& v, std::map<FFVar,double,lt_FFVar> const& pt,
                      std::vector<double> const& var )
{ return oc.eval_colloc<double>( v, pt, var.data(), nullptr, nullptr ); }

static double deriv_interp( OCFESLV const& oc, FFVar const& v, FFVar const& dom,
                            std::map<FFVar,double,lt_FFVar> const& pt,
                            std::vector<double> const& var )
{
  double lo = oc.var_domain().at(dom).lo_dom;
  double up = oc.var_domain().at(dom).up_dom;
  double x  = pt.at(dom);
  double h  = 1e-6*std::max(1.0,up-lo);
  std::map<FFVar,double,lt_FFVar> pm=pt, pp=pt;
  if( x-h < lo ){ pp[dom]=x+h; pm[dom]=x; return (interp(oc,v,pp,var)-interp(oc,v,pm,var))/h; }
  if( x+h > up ){ pp[dom]=x; pm[dom]=x-h; return (interp(oc,v,pp,var)-interp(oc,v,pm,var))/h; }
  pp[dom]=x+h; pm[dom]=x-h; return (interp(oc,v,pp,var)-interp(oc,v,pm,var))/(2*h);
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
  auto nodes = oc.node_colloc( st );
  for( size_t i=0; i<nodes.size(); ++i ){
    std::vector<long long> key;
    key.reserve(nodes[i].size());
    for( double c: nodes[i] ) key.push_back( static_cast<long long>( std::llround( c * 1.0e12 ) ) );
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
  double max_pair = 0.0;
  size_t off = 0;
  for( auto const& st: oc.states_colloc() ){
    auto nodes = oc.node_colloc(st);
    DuplicateSpread const d = duplicate_node_spread( oc, st, var, off );
    max_pair = std::max(max_pair, d.max_pair_spread);
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


static double exact_value_for_state( std::string const& nm, std::vector<double> const& xy, Par const& p )
{
  // OCFESLV::node_colloc() orders tensor coordinates by the sorted dependency
  // set; with variables created as z,rd,rw this gives xy[0]=z and xy[1]=r
  // for the two-dimensional membrane states, matching OCFE_PDE3_solve.cpp.
  double const z = xy[0];
  double const r = xy.size() > 1 ? xy[1] : 0.0;
  if( nm == "Cd(rd,z)" ) return Cd_exact(r,z,p);
  if( nm == "Cw(rw,z)" ) return Cw_exact(r,z,p);
  if( nm == "Sw(rw,z)" ) return Sw_exact(r,z,p);
  if( nm.find("Cd") != std::string::npos && nm.find("Drd_") != std::string::npos ) return Cd_r_exact(r,z,p);
  if( nm.find("Cw") != std::string::npos && nm.find("Drw_") != std::string::npos ) return Cw_r_exact(r,z,p);
  if( nm.find("Sw") != std::string::npos && nm.find("Drw_") != std::string::npos ) return Sw_r_exact(r,z,p);
  return 0.0;
}

static double constant_initial_value_for_state( std::string const& nm )
{
  if( nm == "Cd(rd,z)" ) return 0.45;
  if( nm == "Cw(rw,z)" ) return 0.10;
  if( nm == "Sw(rw,z)" ) return 0.80;
  return 0.0;
}

} // namespace

// ============================================================================
// REPLACEMENT for OCFE_PDE4_solve.cpp main()  (paste over everything from the
// existing `int main()` to end of file; keep all helpers and #defines above).
//
// What changed vs the original
// -----------------------------
//  * The single compile-time mode (TEST_MEM_IMPOSE) is replaced by a runtime
//    sweep over IC_WEAK, IC_TRACE and IC_STRONG.  Each mode builds its OWN
//    fresh FFGraph + OCFESLV, so the three solves are fully independent.
//  * Per-mode pass/fail, final residual, post-solve duplicate-node spread and
//    the Cd/Cw/Sw exact-solution errors are collected into a summary table so
//    the z-seam behaviour across modes is visible at a glance.
//  * Plot files are written per mode (suffix _weak/_trace/_strong) so they do
//    not clobber each other.
//
// This is the validation vehicle for the header fix.  Without the header fix
// the membrane primitives Cd/Cw/Sw will still float in z and the per-mode
// "post-solve max duplicate spread" / "max |Cd-Cd_exact|" checks will FAIL.
// With the fix all three modes should drive the z-jump and the exact errors
// to tolerance.
// ============================================================================

struct ModeResult
{
  std::string name;
  bool        ok          = false;
  bool        solved      = false;
  size_t      nVar=0, nEqn=0, nTrace=0;
  double      final_res   = 0.0;
  double      max_spread  = 0.0;
  double      eCd=0.0, eCw=0.0, eSw=0.0;
};

static ModeResult run_mode( OCFESLV::Options::ImpositionType imp,
                            std::string const& name,
                            std::string const& suffix,
                            Par const& p )
{
  ModeResult R; R.name = name;
  bool ok = true;

  std::cout << "\n========== manufactured membrane-only dry/wet test =========="
            << "\nimposition: " << name
            << ", finite elements: z=" << TEST_MEM_NEL_Z
            << " r=" << TEST_MEM_NEL_R
            << ", nodes/element: z=" << TEST_MEM_NZ
            << " r=" << TEST_MEM_NR
            << ", init_mode=" << TEST_MEM_INIT_MODE << "\n";
  std::cout << "dry-side boundary: Cd(rd=0,z) = "
            << std::scientific << p.Cdry0 << "\n";

  // --- model (identical to the original, rebuilt per mode) ------------------
  FFGraph DAG;
  FFVar z  = DAG.add_var( "z" );
  FFVar rd = DAG.add_var( "rd" );
  FFVar rw = DAG.add_var( "rw" );
  FFVar Cd = DAG.add_var( "Cd(rd,z)" );
  FFVar Cw = DAG.add_var( "Cw(rw,z)" );
  FFVar Sw = DAG.add_var( "Sw(rw,z)" );

  FFPartial OpP;

  FFVar Z   = 3.0*z*z - 2.0*z*z*z;
  FFVar CdU = p.Acd0 + p.AcdZ*Z;
  FFVar qdw = p.qdw0 + p.qdwZ*Z;
  FFVar CdD0 = p.dCd0*(0.5+0.5*Z);
  FFVar CdE    = HERM<FFVar>( rd, FFVar(p.Cdry0), CdD0, CdU, qdw/p.Dmd );
  FFVar CdE_rr = HERM_RR<FFVar>( rd, FFVar(p.Cdry0), CdD0, CdU, qdw/p.Dmd );

  FFVar Cw0 = p.H*CdU;
  FFVar CwD0 = qdw/p.DmwC;
  FFVar CwU = p.CwU0 + p.CwUZ*Z;
  FFVar CwD1 = p.CwD1Z*Z;
  FFVar CwE    = HERM<FFVar>( rw, Cw0, CwD0, CwU, CwD1 );
  FFVar CwE_rr = HERM_RR<FFVar>( rw, Cw0, CwD0, CwU, CwD1 );

  FFVar Sw0 = p.Sw0 + 0.08*Z;
  FFVar SwU = p.Sw0 + p.SwZ*Z;
  FFVar SwD1 = p.SwD1Z*Z;
  FFVar SwE    = HERM<FFVar>( rw, Sw0, FFVar(0.0), SwU, SwD1 );
  FFVar SwE_rr = HERM_RR<FFVar>( rw, Sw0, FFVar(0.0), SwU, SwD1 );

  FFVar Rw  = p.krw*Cw*Sw;
  FFVar RwE = p.krw*CwE*SwE;

  FFVar FD_C = -p.Dmd*CdE_rr;
  FFVar FW_C = -p.DmwC*CwE_rr + RwE;
  FFVar FW_S = -p.DmwS*SwE_rr + p.nuS*RwE;

  FFVar DRY_C = -p.Dmd*OpP(Cd,{rd,2}) - FD_C;
  FFVar WET_C = -p.DmwC*OpP(Cw,{rw,2}) + Rw - FW_C;
  FFVar WET_S = -p.DmwS*OpP(Sw,{rw,2}) + p.nuS*Rw - FW_S;

  FFVar D_C_IN   = Cd - p.Cdry0;
  FFVar DW_VAL   = Cw - p.H*Cd;
  FFVar DW_FLUX  = p.Dmd*OpP(Cd,rd) - p.DmwC*OpP(Cw,rw);
  FFVar DW_SFLUX = OpP(Sw,rw);
  FFVar W_C_OUT  = Cw - CwU;
  FFVar W_S_OUT  = Sw - SwU;

  OCFESLV oc(&DAG);
  oc.add_domain( z,  FFDom(0., p.L, TEST_MEM_NEL_Z, FFDom::CGL, TEST_MEM_NZ) );
  oc.add_domain( rd, FFDom(0., 1.0, TEST_MEM_NEL_R, FFDom::CGL, TEST_MEM_NR) );
  oc.add_domain( rw, FFDom(0., 1.0, TEST_MEM_NEL_R, FFDom::CGL, TEST_MEM_NR) );
  oc.add_state( Cd, {rd,z} );
  oc.add_state( Cw, {rw,z} );
  oc.add_state( Sw, {rw,z} );

  oc.update_ref( Cd, [&]( OCFESLV::t_Coord const& c ){ return Cd_exact(c.at(rd),c.at(z),p); } );
  oc.update_ref( Cw, [&]( OCFESLV::t_Coord const& c ){ return Cw_exact(c.at(rw),c.at(z),p); } );
  oc.update_ref( Sw, [&]( OCFESLV::t_Coord const& c ){ return Sw_exact(c.at(rw),c.at(z),p); } );

  OCFESLV::EqnOptions dry_opt( OCFESLV::EqnRole::INTERIOR, 1 );
  OCFESLV::EqnOptions wet_opt( OCFESLV::EqnRole::INTERIOR, 2 );
  OCFESLV::EqnOptions bnd_opt( OCFESLV::EqnRole::BOUNDARY, 3 );
  OCFESLV::EqnOptions int_opt( OCFESLV::EqnRole::INTERFACE, 4 );

  int const R_INT = FFDom::ALL - FFDom::LB - FFDom::UB;

  oc.add_equation( DRY_C,    {rd,z},    {R_INT,FFDom::ALL},               dry_opt );
  oc.add_equation( WET_C,    {rw,z},    {R_INT,FFDom::ALL},               wet_opt );
  oc.add_equation( WET_S,    {rw,z},    {R_INT,FFDom::ALL},               wet_opt );
  oc.add_equation( D_C_IN,   {rd,z},    {FFDom::LB,FFDom::ALL},           bnd_opt );
  oc.add_equation( DW_VAL,   {rd,rw,z}, {FFDom::UB,FFDom::LB,FFDom::ALL}, int_opt );
  oc.add_equation( DW_FLUX,  {rd,rw,z}, {FFDom::UB,FFDom::LB,FFDom::ALL}, int_opt );
  oc.add_equation( DW_SFLUX, {rw,z},    {FFDom::LB,FFDom::ALL},           int_opt );
  oc.add_equation( W_C_OUT,  {rw,z},    {FFDom::UB,FFDom::ALL},           bnd_opt );
  oc.add_equation( W_S_OUT,  {rw,z},    {FFDom::UB,FFDom::ALL},           bnd_opt );

  oc.reset_evolution_domain();
  oc.options.REDUCE.ORDER   = TEST_MEM_REDUCE_ORDER ? OCFESLV::Options::RED_FULL : OCFESLV::Options::RED_NONE;
  oc.options.CLASSIFY.MODE       = TEST_MEM_CLASSIFY ? OCFESLV::Options::CLASS_AUTO : OCFESLV::Options::CLASS_NONE;
  oc.options.INTERFACE.TYPE = OCFESLV::Options::IC_AUTO;
  oc.options.INTERFACE.IMPOSITION = imp;                 // <-- swept, not compile-time
  oc.options.INTERFACE.SAT_SIGMA0 = TEST_MEM_SAT_SIGMA0;

  if( !oc.setup() ){
    std::cerr << "ERROR: OCFESLV::setup() failed for " << name << "\n";
    R.ok = false;
    return R;
  }
#if TEST_MEM_PRINT_OC
  std::cout << oc;
#endif

  size_t const nVar=oc.n_colloc_sta();
  size_t const nEqn=oc.n_colloc_eqn();
  size_t const nTrace=oc.n_colloc_trace();
  R.nVar=nVar; R.nEqn=nEqn; R.nTrace=nTrace;
  std::cout << "nVar=" << nVar << " nEqn=" << nEqn << " nTrace=" << nTrace
            << " square=" << (nVar==nEqn?"yes":"no") << "\n";
  ok &= (nVar == nEqn);

#if TEST_MEM_CLASSIFY
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
#endif

  std::cout << "States after setup:";
  for( auto const& st : oc.states_colloc() ) std::cout << " " << st.name();
  std::cout << "\n";
  std::cout << "Auxiliary states introduced: "
            << (oc.states_colloc().size() >= 3 ? oc.states_colloc().size()-3 : 0) << "\n";
  if( nTrace ) std::cout << "Trace/tau variables appended: " << nTrace << " (initialised to zero)\n";

  size_t const nState = nVar >= nTrace ? nVar - nTrace : nVar;
  std::vector<double> xExact;
  xExact.reserve(nVar);
  for( auto const& st : oc.states_colloc() ){
    std::string const nm = st.name();
    auto nodes = oc.node_colloc(st);
    for( size_t i=0; i<nodes.size(); ++i )
      xExact.push_back( exact_value_for_state(nm,nodes[i],p) );
  }
  if( xExact.size() != nState ){
    std::cerr << "ERROR(" << name << "): exact ordinary-state vector size mismatch: got="
              << xExact.size() << " expected=" << nState << "\n";
    R.ok = false; return R;
  }
  xExact.resize(nVar,0.0);

  std::vector<double> res(nEqn,123456.0);
  if( !oc.eval(res.data(),nullptr,xExact.data(),nullptr,nullptr) ){
    std::cerr << "ERROR(" << name << "): exact residual evaluation failed\n";
    R.ok = false; return R;
  }
  size_t first_unwritten=nEqn;
  for( size_t i=0; i<nEqn; ++i ) if( res[i] == 123456.0 ){ first_unwritten=i; break; }
  ok &= check_close("all residual rows written before solve", first_unwritten==nEqn?0.0:1.0, 0.0);
  print_residuals("Reference-profile residual", res);

  print_duplicate_spreads( oc, xExact, "manufactured reference" );

  // --- initial guess (same INIT_MODE convention as the original) -----------
  std::vector<double> x;
#if TEST_MEM_INIT_MODE == 0
  x = xExact;
  std::cout << "Initialisation: exact collocation profile\n";
#elif TEST_MEM_INIT_MODE == 1
  x = xExact;
  for( size_t i=0; i<x.size(); ++i )
    x[i] += 1.0e-3*std::sin(0.37*double(i+1))*std::max(1.0,std::abs(x[i]));
  std::cout << "Initialisation: exact profile plus small perturbation\n";
#else // TEST_MEM_INIT_MODE == 2
  x.reserve(nVar);
  for( auto const& st : oc.states_colloc() ){
    std::string const nm=st.name();
    auto nodes=oc.node_colloc(st);
    for( size_t i=0; i<nodes.size(); ++i ) x.push_back( constant_initial_value_for_state(nm) );
  }
  if( x.size() != nState ){
    std::cerr << "ERROR(" << name << "): constant initial vector size mismatch\n";
    R.ok = false; return R;
  }
  x.resize(nVar,0.0);
  std::cout << "Initialisation: constant primitive fields, zero auxiliary/trace variables\n";
#endif

  if( !oc.eval(res.data(),nullptr,x.data(),nullptr,nullptr) ){
    std::cerr << "ERROR(" << name << "): initial residual evaluation failed\n";
    R.ok = false; return R;
  }
  print_residuals("Initial residual", res);

  // --- solve ----------------------------------------------------------------
  oc.options.SOLVE.MAX_ITER = TEST_MEM_MAXIT;
  oc.options.SOLVE.RES_TOL  = TEST_MEM_SOLVE_TOL;
 #if defined(CRONOS__WITH_SPQR) && defined(TEST_LAP_SPQR)
  // item 13: route the augmented (IC_TRACE/IC_STRONG) solve through sparse
  // rank-revealing QR instead of the JtJ normal equations.  Requires the header
  // built with -DCRONOS__WITH_SPQR and SuiteSparse linked.
  oc.options.SOLVE.FACTORIZATION = OCFESLV::Options::SOLVE_SPQR;
#endif

 OCFESLV::SolveReport srep = oc.solve( x.data() );
  R.solved = srep.converged;
  ok &= R.solved;
  if( !eval_residual(oc,x,res) ){ R.ok=false; return R; }
  print_residuals("Final residual", res);
  print_largest_residuals("Final residual", res);
  R.final_res  = max_abs(res);
  R.max_spread = print_duplicate_spreads( oc, x, "post solve" );
  ok &= check_close("final max residual", R.final_res, TEST_MEM_SOLVE_TOL);
  ok &= check_close("post-solve max duplicate spread", R.max_spread, 3e-6);

  // --- exact-solution comparison + interface diagnostics --------------------
  auto zNodes = physical_nodes( oc.var_domain().at(z) );
  auto rNodes = physical_nodes( oc.var_domain().at(rd) );

  double dry_bnd=0., dwv=0., dwf=0., dws=0., woutC=0., woutS=0.;
  for( double zz: zNodes ){
    std::map<FFVar,double,lt_FFVar> pdL{{rd,0.0},{z,zz}}, pdU{{rd,1.0},{z,zz}};
    std::map<FFVar,double,lt_FFVar> pwL{{rw,0.0},{z,zz}}, pwU{{rw,1.0},{z,zz}};
    dry_bnd=std::max(dry_bnd,std::abs(interp(oc,Cd,pdL,x)-p.Cdry0));
    dwv=std::max(dwv,std::abs(interp(oc,Cw,pwL,x)-p.H*interp(oc,Cd,pdU,x)));
    dwf=std::max(dwf,std::abs(p.Dmd*deriv_interp(oc,Cd,rd,pdU,x)-p.DmwC*deriv_interp(oc,Cw,rw,pwL,x)));
    dws=std::max(dws,std::abs(deriv_interp(oc,Sw,rw,pwL,x)));
    woutC=std::max(woutC,std::abs(interp(oc,Cw,pwU,x)-CwU_exact(zz,p)));
    woutS=std::max(woutS,std::abs(interp(oc,Sw,pwU,x)-SwU_exact(zz,p)));
    for( double rr: rNodes ){
      std::map<FFVar,double,lt_FFVar> pd{{rd,rr},{z,zz}}, pw{{rw,rr},{z,zz}};
      R.eCd=std::max(R.eCd,std::abs(interp(oc,Cd,pd,x)-Cd_exact(rr,zz,p)));
      R.eCw=std::max(R.eCw,std::abs(interp(oc,Cw,pw,x)-Cw_exact(rr,zz,p)));
      R.eSw=std::max(R.eSw,std::abs(interp(oc,Sw,pw,x)-Sw_exact(rr,zz,p)));
    }
  }

  std::cout << "\nExact-solution comparison:\n";
  ok &= check_close("max |Cd-Cd_exact|", R.eCd, 3e-3);
  ok &= check_close("max |Cw-Cw_exact|", R.eCw, 2e-3);
  ok &= check_close("max |Sw-Sw_exact|", R.eSw, 2e-3);

  std::cout << "\nIndependent interface/boundary diagnostics:\n";
  ok &= check_close("dry-side constant Cd value", dry_bnd, 3e-6);
  ok &= check_close("dry/wet value", dwv, 3e-6);
  ok &= check_close("dry/wet CO2 flux", dwf, 3e-6);
  ok &= check_close("dry/wet solvent zero-flux", dws, 3e-6);
  ok &= check_close("wet outer Cw value", woutC, 3e-6);
  ok &= check_close("wet outer Sw value", woutS, 3e-6);

  // --- per-mode plot files --------------------------------------------------
  std::string const prefix = std::string(TEST_MEM_OUTPUT_PREFIX) + "_" + suffix;
  {
    std::ofstream out(prefix + std::string("_dry.out"));
    out << "# z rd Cd Cd_exact\n";
    for( double zz: zNodes ){ for( double rr: rNodes ){
      std::map<FFVar,double,lt_FFVar> pt{{rd,rr},{z,zz}};
      out << std::setprecision(16) << zz << " " << rr << " "
          << interp(oc,Cd,pt,x) << " " << Cd_exact(rr,zz,p) << "\n";
    } out << "\n"; }
  }
  {
    std::ofstream out(prefix + std::string("_wet.out"));
    out << "# z rw Cw Sw Cw_exact Sw_exact\n";
    for( double zz: zNodes ){ for( double rr: rNodes ){
      std::map<FFVar,double,lt_FFVar> pt{{rw,rr},{z,zz}};
      out << std::setprecision(16) << zz << " " << rr << " "
          << interp(oc,Cw,pt,x) << " " << interp(oc,Sw,pt,x)
          << " " << Cw_exact(rr,zz,p) << " " << Sw_exact(rr,zz,p) << "\n";
    } out << "\n"; }
  }
  std::cout << "\nWrote plot files: " << prefix << "_{dry,wet}.out\n";
  std::cout << "Mode " << name << ": " << (ok?"PASS":"FAIL") << "\n";

  R.ok = ok;
  return R;
}

int main()
{
  Par const p;

  std::vector<ModeResult> results;
  results.push_back( run_mode( OCFESLV::Options::IC_WEAK,   "IC_WEAK",   "weak",   p ) );
  results.push_back( run_mode( OCFESLV::Options::IC_TRACE,  "IC_TRACE",  "trace",  p ) );
  results.push_back( run_mode( OCFESLV::Options::IC_STRONG, "IC_STRONG", "strong", p ) );

  std::cout << "\n==================== PDE4 mode sweep summary ====================\n";
  std::cout << std::left << std::setw(10) << "mode"
            << std::right << std::setw(8)  << "nTrace"
            << std::setw(14) << "final|r|"
            << std::setw(14) << "dup-spread"
            << std::setw(12) << "|dCd|"
            << std::setw(12) << "|dCw|"
            << std::setw(12) << "|dSw|"
            << std::setw(8)  << "result" << "\n";
  bool all_ok = true;
  for( auto const& r : results ){
    all_ok &= r.ok;
    std::cout << std::left << std::setw(10) << r.name
              << std::right << std::setw(8) << r.nTrace
              << std::scientific << std::setprecision(3)
              << std::setw(14) << r.final_res
              << std::setw(14) << r.max_spread
              << std::setw(12) << r.eCd
              << std::setw(12) << r.eCw
              << std::setw(12) << r.eSw
              << std::setw(8)  << (r.ok?"PASS":"FAIL") << "\n";
  }
  std::cout << "=================================================================\n";
  std::cout << "Manufactured membrane-only dry/wet sweep: "
            << (all_ok?"PASS":"FAIL") << "\n";
  return all_ok?0:1;
}
