// Copyright (C) Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the EPL-2.0 with GPL-2.0-or-later as a Secondary License; see the LICENSE file.
// SPDX-License-Identifier: EPL-2.0 OR GPL-2.0-or-later
//
// OCFE_PDE3_solve.cpp
// --------------------
// Manufactured-solution multidomain membrane-contactor regression test for
// OCFESLV.  This file restores the dry and wet membrane compartments that were
// omitted in OCFE_PDE2_solve.cpp, while retaining an exact smooth reference solution
// for all primitive states:
//
//   gas tube      : Cg(z), Vg(z)
//   dry membrane  : Cd(rd,z)
//   wet membrane  : Cw(rw,z), Sw(rw,z)
//   liquid shell  : Cl(rl,z), Sl(rl,z)
//
// The exact radial profiles are constructed with Hermite polynomials so that
// gas/dry, dry/wet, wet/liquid and wall value/flux conditions are satisfied by
// construction.  Manufactured source terms force the interior equations, and
// the final numerical solution is compared with the exact profile field by
// field.  Default iterations start from a small perturbation of the exact profile;
// TEST_MDW_INIT_MODE=2 can be used to start from constant primitive fields and
// zero auxiliary derivative/trace variables.
//
// Build-time switches:
//   Runtime sweep: IC_WEAK, IC_TRACE, and IC_STRONG are run sequentially.
//   -DTEST_MDW_NEL_Z        number of axial finite elements
//   -DTEST_MDW_NEL_R        number of radial finite elements per membrane/liquid
//   -DTEST_MDW_NZ           axial nodes per element
//   -DTEST_MDW_NR           radial nodes per element
//   -DTEST_MDW_INIT_MODE    0 exact, 1 perturbed exact, 2 constant primitives

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

#ifndef TEST_MDW_NEL_Z
#define TEST_MDW_NEL_Z 3
#endif
#ifndef TEST_MDW_NEL_R
#define TEST_MDW_NEL_R 3
#endif
#ifndef TEST_MDW_NZ
#define TEST_MDW_NZ 7
#endif
#ifndef TEST_MDW_NR
#define TEST_MDW_NR 7
#endif
#ifndef TEST_MDW_REDUCE_ORDER
#define TEST_MDW_REDUCE_ORDER 1
#endif
#ifndef TEST_MDW_CLASSIFY
#define TEST_MDW_CLASSIFY 1
#endif
#ifndef TEST_MDW_OUTPUT_PREFIX
#define TEST_MDW_OUTPUT_PREFIX "OCFE_PDE3_solve"
#endif
#ifndef TEST_MDW_MAXIT
#define TEST_MDW_MAXIT 30
#endif
#ifndef TEST_MDW_SOLVE_TOL
#define TEST_MDW_SOLVE_TOL 1e-9
#endif
#ifndef TEST_MDW_SAT_SIGMA0
#define TEST_MDW_SAT_SIGMA0 1.0
#endif
#ifndef TEST_MDW_PRINT_OC
#define TEST_MDW_PRINT_OC 1
#endif
#ifndef TEST_MDW_INIT_MODE
#define TEST_MDW_INIT_MODE 2
#endif
// Gate 1 (rank / well-posedness): smallest singular value of the full Jacobian
// at the manufactured solution.  This costs a dense O(n^3) SVD per (epszL,mode),
// so it is OFF by default; build with -DTEST_MDW_RANK_GATE=1 to enable.  When on
// the check is mode-aware: sigma_min is a genuine conditioning number only for
// non-augmented systems (nTrace==0), so the pass/fail reduction is taken over
// those alone.  Augmented (trace/tau) modes are rank-deficient by construction
// (structural null space; QR min-norm makes that benign), so their sigma_min is
// reported for information only and excluded from the verdict.
#ifndef TEST_MDW_RANK_GATE
#define TEST_MDW_RANK_GATE 0
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
  double DlC    = 1.6e-2;
  double DlS    = 7.0e-3;
  double epszL  = 1.0e-2;
  double krw    = 0.28;
  double krl    = 0.42;
  double nuS    = 1.0;
  double Kvg    = 0.25;
  double Kcg    = 0.65;
  double betaG  = 0.10;
  double CgIn   = 1.0;
  double VgIn   = 1.0;
  double ClIn   = 0.015;
  double SlIn   = 1.0;

  // Exact-profile amplitudes.  These keep all concentrations positive while
  // giving visible solvent depletion near the wet/liquid interface.
  double Acg    = 0.12;
  double Avg    = 0.05;
  double Acd0   = 0.16;
  double AcdZ   = 0.04;
  double Acl0   = 0.08;
  double AclR   = -0.020;
  double As     = 0.45;
  double Asw    = 0.05;
  double qdw0   = -5.0e-3;
  double qdwZ   = -1.0e-3;
  double dCd0   = -0.10;
};

static double sqr( double x ){ return x*x; }
static double Zs( double z ){ return 3.0*z*z - 2.0*z*z*z; }
static double Zs_z( double z ){ return 6.0*z - 6.0*z*z; }
static double Zs_zz( double z ){ return 6.0 - 12.0*z; }
static double Bz( double z ){ return z*z*sqr(1.0-z); }
static double Bz_z( double z ){ return 2.0*z - 6.0*z*z + 4.0*z*z*z; }
static double Bz_zz( double z ){ return 2.0 - 12.0*z + 12.0*z*z; }
static double Rr( double r ){ return r*(2.0-r); }
static double Rr_r( double r ){ return 2.0 - 2.0*r; }
static double Rr_rr(){ return -2.0; }
static double Hs( double r ){ return 1.0 - 3.0*r*r + 2.0*r*r*r; }
static double Hs_r( double r ){ return -6.0*r + 6.0*r*r; }
static double Hs_rr( double r ){ return -6.0 + 12.0*r; }

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

static double Cg_exact( double z, Par const& p ){ return p.CgIn - p.Acg*sqr(1.0-z); }
static double Cg_z_exact( double z, Par const& p ){ return 2.0*p.Acg*(1.0-z); }
static double Vg_exact( double z, Par const& p ){ return p.VgIn - p.Avg*sqr(1.0-z); }
static double Vg_z_exact( double z, Par const& p ){ return 2.0*p.Avg*(1.0-z); }
static double CdU_exact( double z, Par const& p ){ return p.Acd0 + p.AcdZ*Zs(z); }
static double CdU_z_exact( double z, Par const& p ){ return p.AcdZ*Zs_z(z); }
static double qdw_exact( double z, Par const& p ){ return p.qdw0 + p.qdwZ*Zs(z); }
static double qdw_z_exact( double z, Par const& p ){ return p.qdwZ*Zs_z(z); }
static double Cd_d0_exact( double z, Par const& p ){ return p.dCd0 * ( 0.5 + 0.5*Zs(z) ); }
static double Cd_d0_z_exact( double z, Par const& p ){ return 0.5*p.dCd0*Zs_z(z); }
static double Cd_exact( double r, double z, Par const& p )
{ return herm( r, Cg_exact(z,p), Cd_d0_exact(z,p), CdU_exact(z,p), qdw_exact(z,p)/p.Dmd ); }
static double Cd_r_exact( double r, double z, Par const& p )
{ return herm_r( r, Cg_exact(z,p), Cd_d0_exact(z,p), CdU_exact(z,p), qdw_exact(z,p)/p.Dmd ); }
static double Cd_rr_exact( double r, double z, Par const& p )
{ return herm_rr( r, Cg_exact(z,p), Cd_d0_exact(z,p), CdU_exact(z,p), qdw_exact(z,p)/p.Dmd ); }
static double Cd_z_exact( double r, double z, Par const& p )
{ return herm( r, Cg_z_exact(z,p), Cd_d0_z_exact(z,p), CdU_z_exact(z,p), qdw_z_exact(z,p)/p.Dmd ); }

static double Cl0_exact( double z, Par const& p ){ return p.ClIn + p.Acl0*Zs(z); }
static double Cl0_z_exact( double z, Par const& p ){ return p.Acl0*Zs_z(z); }
static double Cl0_zz_exact( double z, Par const& p ){ return p.Acl0*Zs_zz(z); }
static double ClA_exact( double z, Par const& p ){ return p.AclR*Bz(z); }
static double ClA_z_exact( double z, Par const& p ){ return p.AclR*Bz_z(z); }
static double ClA_zz_exact( double z, Par const& p ){ return p.AclR*Bz_zz(z); }
static double Cl_exact( double r, double z, Par const& p )
{ return Cl0_exact(z,p) + ClA_exact(z,p)*Rr(r); }
static double Cl_r_exact( double r, double z, Par const& p )
{ return ClA_exact(z,p)*Rr_r(r); }
static double Cl_rr_exact( double, double z, Par const& p )
{ return ClA_exact(z,p)*Rr_rr(); }
static double Cl_z_exact( double r, double z, Par const& p )
{ return Cl0_z_exact(z,p) + ClA_z_exact(z,p)*Rr(r); }
static double Cl_zz_exact( double r, double z, Par const& p )
{ return Cl0_zz_exact(z,p) + ClA_zz_exact(z,p)*Rr(r); }

static double Cw0_exact( double z, Par const& p ){ return p.H*CdU_exact(z,p); }
static double Cw0_z_exact( double z, Par const& p ){ return p.H*CdU_z_exact(z,p); }
static double Cw_d0_exact( double z, Par const& p ){ return qdw_exact(z,p)/p.DmwC; }
static double Cw_d0_z_exact( double z, Par const& p ){ return qdw_z_exact(z,p)/p.DmwC; }
static double Cw_d1_exact( double z, Par const& p ){ return p.DlC*Cl_r_exact(0.0,z,p)/p.DmwC; }
static double Cw_d1_z_exact( double z, Par const& p ){ return p.DlC*(p.AclR*Bz_z(z)*2.0)/p.DmwC; }
static double Cw_exact( double r, double z, Par const& p )
{ return herm( r, Cw0_exact(z,p), Cw_d0_exact(z,p), Cl0_exact(z,p), Cw_d1_exact(z,p) ); }
static double Cw_r_exact( double r, double z, Par const& p )
{ return herm_r( r, Cw0_exact(z,p), Cw_d0_exact(z,p), Cl0_exact(z,p), Cw_d1_exact(z,p) ); }
static double Cw_rr_exact( double r, double z, Par const& p )
{ return herm_rr( r, Cw0_exact(z,p), Cw_d0_exact(z,p), Cl0_exact(z,p), Cw_d1_exact(z,p) ); }
static double Cw_z_exact( double r, double z, Par const& p )
{ return herm( r, Cw0_z_exact(z,p), Cw_d0_z_exact(z,p), Cl0_z_exact(z,p), Cw_d1_z_exact(z,p) ); }

static double Sl_exact( double r, double z, Par const& p )
{ return p.SlIn - p.As*Zs(z)*Hs(r); }
static double Sl_r_exact( double r, double z, Par const& p )
{ return -p.As*Zs(z)*Hs_r(r); }
static double Sl_rr_exact( double r, double z, Par const& p )
{ return -p.As*Zs(z)*Hs_rr(r); }
static double Sl_z_exact( double r, double z, Par const& p )
{ return -p.As*Zs_z(z)*Hs(r); }
static double Sl_zz_exact( double r, double z, Par const& p )
{ return -p.As*Zs_zz(z)*Hs(r); }
static double Sw0_exact( double z, Par const& p ){ return Sl_exact(0.0,z,p) + p.Asw*Zs(z); }
static double Sw0_z_exact( double z, Par const& p ){ return Sl_z_exact(0.0,z,p) + p.Asw*Zs_z(z); }
static double Sw_exact( double r, double z, Par const& p )
{ return HERM<double>( r, Sw0_exact(z,p), 0.0, Sl_exact(0.0,z,p), 0.0 ); }
static double Sw_r_exact( double r, double z, Par const& p )
{ return HERM_R<double>( r, Sw0_exact(z,p), 0.0, Sl_exact(0.0,z,p), 0.0 ); }
static double Sw_rr_exact( double r, double z, Par const& p )
{ return HERM_RR<double>( r, Sw0_exact(z,p), 0.0, Sl_exact(0.0,z,p), 0.0 ); }
static double Sw_z_exact( double r, double z, Par const& p )
{ return HERM<double>( r, Sw0_z_exact(z,p), 0.0, Sl_z_exact(0.0,z,p), 0.0 ); }

static double gas_v_src_exact( double z, Par const& p )
{ return Vg_z_exact(z,p) + p.Kvg*p.Dmd*Cd_r_exact(0.0,z,p); }
static double gas_c_src_exact( double z, Par const& p )
{ return Cg_z_exact(z,p) + p.Kcg*p.Dmd*(1.0-p.betaG*Cg_exact(z,p))/Vg_exact(z,p)*Cd_r_exact(0.0,z,p); }
static double dry_c_src_exact( double r, double z, Par const& p )
{ return -p.Dmd*Cd_rr_exact(r,z,p); }
static double wet_c_src_exact( double r, double z, Par const& p )
{ return -p.DmwC*Cw_rr_exact(r,z,p) + p.krw*Cw_exact(r,z,p)*Sw_exact(r,z,p); }
static double wet_s_src_exact( double r, double z, Par const& p )
{ return -p.DmwS*Sw_rr_exact(r,z,p) + p.nuS*p.krw*Cw_exact(r,z,p)*Sw_exact(r,z,p); }
static double liq_c_src_exact( double r, double z, Par const& p )
{ return (1.0+0.20*r)*Cl_z_exact(r,z,p) - p.DlC*(Cl_rr_exact(r,z,p)+p.epszL*Cl_zz_exact(r,z,p)) + p.krl*Cl_exact(r,z,p)*Sl_exact(r,z,p); }
static double liq_s_src_exact( double r, double z, Par const& p )
{ return (1.0+0.20*r)*Sl_z_exact(r,z,p) - p.DlS*(Sl_rr_exact(r,z,p)+p.epszL*Sl_zz_exact(r,z,p)) + p.nuS*p.krl*Cl_exact(r,z,p)*Sl_exact(r,z,p); }

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
    max_pair = std::max( max_pair, d.max_pair_spread );
    std::cout << "  " << std::setw(18) << st.name()
              << " max_pair=" << std::scientific << std::setprecision(6) << d.max_pair_spread
              << " max_corner=" << d.max_corner_spread
              << " max_multiplicity=" << d.max_multiplicity << "\n";
    off += nodes.size();
  }
  return max_pair;
}

#if TEST_MDW_RANK_GATE
// Dense Jacobian assembly (sparse AD pattern + values via deriv(), densified).
// Retained for the Gate-1 sigma_min(J) well-posedness check; the production solve
// uses the sparse path inside OCFESLV::solve().
static bool dense_jacobian( OCFESLV& oc, std::vector<double> const& x, arma::mat& J )
{
  size_t const nEqn=oc.n_colloc_eqn(), nVar=oc.n_colloc_sta();
  std::vector<size_t> nnz(nEqn,0);
  std::vector<std::vector<size_t>> col(nEqn);
  std::vector<size_t*> colptr(nEqn,nullptr);
  if( !oc.deriv(nnz.data(),nullptr) ) return false;
  for( size_t i=0; i<nEqn; ++i ){ col[i].assign(nnz[i],0); if(nnz[i]) colptr[i]=col[i].data(); }
  if( !oc.deriv(nnz.data(),colptr.data()) ) return false;
  size_t nnzsum=0; for( size_t v: nnz ) nnzsum += v;
  std::vector<double> grad(nnzsum,0.0);
  if( !oc.deriv(grad.data(),nullptr,x.data(),nullptr,nullptr) ) return false;
  J.zeros(nEqn,nVar);
  size_t off=0;
  for( size_t i=0; i<nEqn; ++i ){
    for( size_t k=0; k<nnz[i]; ++k ) if( col[i][k] < nVar ) J(i,col[i][k]) = grad[off+k];
    off += nnz[i];
  }
  return true;
}
#endif // TEST_MDW_RANK_GATE

static bool eval_residual( OCFESLV& oc, std::vector<double> const& x, std::vector<double>& r )
{ std::fill(r.begin(),r.end(),0.0); return oc.eval(r.data(),nullptr,x.data(),nullptr,nullptr); }

static double exact_initial_value_for_state( std::string const& nm, std::vector<double> const& xy, Par const& p )
{
  double const z=xy[0];
  double const r=xy.size()>1?xy[1]:0.0;
  if( nm == "Cg(z)" ) return Cg_exact(z,p);
  if( nm == "Vg(z)" ) return Vg_exact(z,p);
  if( nm == "Cd(rd,z)" ) return Cd_exact(r,z,p);
  if( nm == "Cw(rw,z)" ) return Cw_exact(r,z,p);
  if( nm == "Sw(rw,z)" ) return Sw_exact(r,z,p);
  if( nm == "Cl(rl,z)" ) return Cl_exact(r,z,p);
  if( nm == "Sl(rl,z)" ) return Sl_exact(r,z,p);
  if( nm.find("Cd") != std::string::npos ){
    if( nm.find("Drd_") != std::string::npos ) return Cd_r_exact(r,z,p);
    if( nm.find("Dz_")  != std::string::npos ) return Cd_z_exact(r,z,p);
  }
  if( nm.find("Cw") != std::string::npos ){
    if( nm.find("Drw_") != std::string::npos ) return Cw_r_exact(r,z,p);
    if( nm.find("Dz_")  != std::string::npos ) return Cw_z_exact(r,z,p);
  }
  if( nm.find("Sw") != std::string::npos ){
    if( nm.find("Drw_") != std::string::npos ) return Sw_r_exact(r,z,p);
    if( nm.find("Dz_")  != std::string::npos ) return Sw_z_exact(r,z,p);
  }
  if( nm.find("Cl") != std::string::npos ){
    if( nm.find("Drl_") != std::string::npos ) return Cl_r_exact(r,z,p);
    if( nm.find("Dz_")  != std::string::npos ) return Cl_z_exact(r,z,p);
  }
  if( nm.find("Sl") != std::string::npos ){
    if( nm.find("Drl_") != std::string::npos ) return Sl_r_exact(r,z,p);
    if( nm.find("Dz_")  != std::string::npos ) return Sl_z_exact(r,z,p);
  }
  return 0.0;
}
static double constant_initial_value_for_state( std::string const& nm )
{
  if( nm == "Cg(z)" ) return 0.94;
  if( nm == "Vg(z)" ) return 0.975;
  if( nm == "Cd(rd,z)" ) return 0.45;
  if( nm == "Cw(rw,z)" ) return 0.075;
  if( nm == "Sw(rw,z)" ) return 0.78;
  if( nm == "Cl(rl,z)" ) return 0.07;
  if( nm == "Sl(rl,z)" ) return 0.75;
  return 0.0;
}

struct StateExactErrors
{
  double eCg  = 0.0;
  double eVg  = 0.0;
  double eCd  = 0.0;
  double eCw  = 0.0;
  double eSw  = 0.0;
  double eCl  = 0.0;
  double eSl  = 0.0;
  double eAux = 0.0;
};

static bool is_auxiliary_state_name( std::string const& nm )
{
  return !nm.empty() && nm[0] == 'D';
}

static StateExactErrors print_variable_exact_errors
( OCFESLV const& oc, std::vector<double> const& var, Par const& p, std::string const& label )
{
  StateExactErrors out;
  size_t off = 0;
  std::cout << "\nVariable comparison against exact profiles (" << label << "):\n";
  for( auto const& st : oc.states_colloc() ){
    auto nodes = oc.node_colloc( st );
    std::string const nm = st.name();
    if( !is_auxiliary_state_name( nm ) ){
      double maxerr = 0.0, meanerr = 0.0;
      for( size_t i=0; i<nodes.size(); ++i ){
        double const ref = exact_initial_value_for_state( nm, nodes[i], p );
        double const err = std::abs( var[off+i] - ref );
        maxerr = std::max( maxerr, err );
        meanerr += err;
      }
      meanerr /= nodes.empty()? 1.0 : double(nodes.size());
      if( nm == "Cg(z)" ) out.eCg = maxerr;
      else if( nm == "Vg(z)" ) out.eVg = maxerr;
      else if( nm == "Cd(rd,z)" ) out.eCd = maxerr;
      else if( nm == "Cw(rw,z)" ) out.eCw = maxerr;
      else if( nm == "Sw(rw,z)" ) out.eSw = maxerr;
      else if( nm == "Cl(rl,z)" ) out.eCl = maxerr;
      else if( nm == "Sl(rl,z)" ) out.eSl = maxerr;
      std::cout << "  " << std::setw(18) << nm
                << " max|v-v_exact|=" << std::scientific << std::setprecision(4) << maxerr
                << " mean|v-v_exact|=" << meanerr << "\n";
    }
    off += nodes.size();
  }
  return out;
}

static StateExactErrors print_auxiliary_flux_exact_errors
( OCFESLV const& oc, std::vector<double> const& var, Par const& p, std::string const& label )
{
  StateExactErrors out;
  size_t off = 0;
  std::cout << "\nAuxiliary flux/derivative comparison against exact profiles (" << label << "):\n";
  for( auto const& st : oc.states_colloc() ){
    auto nodes = oc.node_colloc( st );
    std::string const nm = st.name();
    if( is_auxiliary_state_name( nm ) ){
      double maxerr = 0.0, meanerr = 0.0;
      for( size_t i=0; i<nodes.size(); ++i ){
        double const ref = exact_initial_value_for_state( nm, nodes[i], p );
        double const err = std::abs( var[off+i] - ref );
        maxerr = std::max( maxerr, err );
        meanerr += err;
      }
      meanerr /= nodes.empty()? 1.0 : double(nodes.size());
      out.eAux = std::max( out.eAux, maxerr );
      std::cout << "  " << std::setw(18) << nm
                << " max|aux-aux_exact|=" << std::scientific << std::setprecision(4) << maxerr
                << " mean|aux-aux_exact|=" << meanerr << "\n";
    }
    off += nodes.size();
  }
  return out;
}

} // namespace

struct ModeResult
{
  std::string name;
  bool        ok      = false;
  bool        solved  = false;
  size_t      nVar    = 0;
  size_t      nEqn    = 0;
  size_t      nTrace  = 0;
  double      final_res  = std::numeric_limits<double>::infinity();
  double      max_spread = 0.0;
  double      eCg = 0.0;
  double      eVg = 0.0;
  double      eCd = 0.0;
  double      eCw = 0.0;
  double      eSw = 0.0;
  double      eCl = 0.0;
  double      eSl = 0.0;
  double      eAux = 0.0;
  double      epszL = 0.0;
  double      sigma_min = 0.0;
};

static ModeResult run_mode( OCFESLV::Options::ImpositionType imp, std::string const& strimp, std::string const& suffix, double epszL_arg )
{
  ModeResult MR;
  MR.name = strimp;
  MR.epszL = epszL_arg;
  bool ok=true;
  Par p; p.epszL = epszL_arg;


  std::cout << "\n========== manufactured dry/wet/liquid MBC test ==========" << "\n";
  std::cout << "imposition: " << strimp
            << ", finite elements: z=" << TEST_MDW_NEL_Z
            << " r=" << TEST_MDW_NEL_R
            << ", nodes/element: z=" << TEST_MDW_NZ
            << " r=" << TEST_MDW_NR
            << ", init_mode=" << TEST_MDW_INIT_MODE << "\n";

  FFGraph DAG;
  FFVar z  = DAG.add_var( "z" );
  FFVar rd = DAG.add_var( "rd" );
  FFVar rw = DAG.add_var( "rw" );
  FFVar rl = DAG.add_var( "rl" );
  FFVar Cg = DAG.add_var( "Cg(z)" );
  FFVar Vg = DAG.add_var( "Vg(z)" );
  FFVar Cd = DAG.add_var( "Cd(rd,z)" );
  FFVar Cw = DAG.add_var( "Cw(rw,z)" );
  FFVar Sw = DAG.add_var( "Sw(rw,z)" );
  FFVar Cl = DAG.add_var( "Cl(rl,z)" );
  FFVar Sl = DAG.add_var( "Sl(rl,z)" );

  FFPartial OpP;

  FFVar omz = 1.0-z;
  FFVar Z  = 3.0*z*z - 2.0*z*z*z;
  FFVar Z_z = 6.0*z - 6.0*z*z;
  FFVar Z_zz = 6.0 - 12.0*z;
  FFVar B = z*z*omz*omz;
  FFVar B_z = 2.0*z - 6.0*z*z + 4.0*z*z*z;
  FFVar B_zz = 2.0 - 12.0*z + 12.0*z*z;

  FFVar CgE = p.CgIn - p.Acg*omz*omz;
  FFVar CgE_z = 2.0*p.Acg*omz;
  FFVar VgE = p.VgIn - p.Avg*omz*omz;
  FFVar VgE_z = 2.0*p.Avg*omz;

  FFVar CdU = p.Acd0 + p.AcdZ*Z;
  FFVar CdU_z = p.AcdZ*Z_z;
  FFVar qdw = p.qdw0 + p.qdwZ*Z;
  FFVar qdw_z = p.qdwZ*Z_z;
  FFVar CdD0 = p.dCd0*(0.5+0.5*Z);
  FFVar CdD0_z = 0.5*p.dCd0*Z_z;
  FFVar CdE = HERM<FFVar>( rd, CgE, CdD0, CdU, qdw/p.Dmd );
  FFVar CdE_r = HERM_R<FFVar>( rd, CgE, CdD0, CdU, qdw/p.Dmd );
  FFVar CdE_rr = HERM_RR<FFVar>( rd, CgE, CdD0, CdU, qdw/p.Dmd );

  FFVar R = rl*(2.0-rl);
  FFVar Cl0 = p.ClIn + p.Acl0*Z;
  FFVar Cl0_z = p.Acl0*Z_z;
  FFVar Cl0_zz = p.Acl0*Z_zz;
  FFVar ClA = p.AclR*B;
  FFVar ClA_z = p.AclR*B_z;
  FFVar ClA_zz = p.AclR*B_zz;
  FFVar ClE = Cl0 + ClA*R;
  FFVar ClE_r = ClA*(2.0-2.0*rl);
  FFVar ClE_rr = -2.0*ClA;
  FFVar ClE_z = Cl0_z + ClA_z*R;
  FFVar ClE_zz = Cl0_zz + ClA_zz*R;

  FFVar Cw0 = p.H*CdU;
  FFVar Cw0_z = p.H*CdU_z;
  FFVar CwD0 = qdw/p.DmwC;
  FFVar CwD0_z = qdw_z/p.DmwC;
  FFVar CwD1 = p.DlC*(p.AclR*B*2.0)/p.DmwC;
  FFVar CwE = HERM<FFVar>( rw, Cw0, CwD0, Cl0, CwD1 );
  FFVar CwE_r = HERM_R<FFVar>( rw, Cw0, CwD0, Cl0, CwD1 );
  FFVar CwE_rr = HERM_RR<FFVar>( rw, Cw0, CwD0, Cl0, CwD1 );

  FFVar H = 1.0 - 3.0*rl*rl + 2.0*rl*rl*rl;
  FFVar H_r = -6.0*rl + 6.0*rl*rl;
  FFVar H_rr = -6.0 + 12.0*rl;
  FFVar SlE = p.SlIn - p.As*Z*H;
  FFVar SlE_r = -p.As*Z*H_r;
  FFVar SlE_rr = -p.As*Z*H_rr;
  FFVar SlE_z = -p.As*Z_z*H;
  FFVar SlE_zz = -p.As*Z_zz*H;

  FFVar Sw0 = (p.SlIn - p.As*Z) + p.Asw*Z;
  FFVar Sw1 = p.SlIn - p.As*Z;
  FFVar SwE = HERM<FFVar>( rw, Sw0, FFVar(0.0), Sw1, FFVar(0.0) );
  FFVar SwE_r = HERM_R<FFVar>( rw, Sw0, FFVar(0.0), Sw1, FFVar(0.0) );
  FFVar SwE_rr = HERM_RR<FFVar>( rw, Sw0, FFVar(0.0), Sw1, FFVar(0.0) );

  FFVar vl = 1.0 + 0.20*rl;
  FFVar Rw = p.krw*Cw*Sw;
  FFVar Rl = p.krl*Cl*Sl;
  FFVar RwE = p.krw*CwE*SwE;
  FFVar RlE = p.krl*ClE*SlE;

  FFVar FG_V = VgE_z + p.Kvg*p.Dmd*HERM_R<FFVar>( FFVar(0.0), CgE, CdD0, CdU, qdw/p.Dmd );
  FFVar FG_C = CgE_z + p.Kcg*p.Dmd*(1.0-p.betaG*CgE)/VgE*HERM_R<FFVar>( FFVar(0.0), CgE, CdD0, CdU, qdw/p.Dmd );
  FFVar FD_C = -p.Dmd*CdE_rr;
  FFVar FW_C = -p.DmwC*CwE_rr + RwE;
  FFVar FW_S = -p.DmwS*SwE_rr + p.nuS*RwE;
  FFVar FL_C = vl*ClE_z - p.DlC*(ClE_rr+p.epszL*ClE_zz) + RlE;
  FFVar FL_S = vl*SlE_z - p.DlS*(SlE_rr+p.epszL*SlE_zz) + p.nuS*RlE;

  FFVar GAS_V = OpP(Vg,z) + p.Kvg*p.Dmd*OpP(Cd,rd) - FG_V;
  FFVar GAS_C = OpP(Cg,z) + p.Kcg*p.Dmd*(1.0-p.betaG*Cg)/Vg*OpP(Cd,rd) - FG_C;
  FFVar DRY_C = -p.Dmd*OpP(Cd,{rd,2}) - FD_C;
  FFVar WET_C = -p.DmwC*OpP(Cw,{rw,2}) + Rw - FW_C;
  FFVar WET_S = -p.DmwS*OpP(Sw,{rw,2}) + p.nuS*Rw - FW_S;
  FFVar LIQ_C = vl*OpP(Cl,z) - p.DlC*(OpP(Cl,{rl,2})+p.epszL*OpP(Cl,{z,2})) + Rl - FL_C;
  FFVar LIQ_S = vl*OpP(Sl,z) - p.DlS*(OpP(Sl,{rl,2})+p.epszL*OpP(Sl,{z,2})) + p.nuS*Rl - FL_S;

  FFVar GAS_C_IN = Cg - CgE;
  FFVar GAS_V_IN = Vg - VgE;
  FFVar GD_VAL   = Cd - Cg;
  FFVar DW_VAL   = Cw - p.H*Cd;
  FFVar DW_FLUX  = p.Dmd*OpP(Cd,rd) - p.DmwC*OpP(Cw,rw);
  FFVar DW_SFLUX = OpP(Sw,rw);
  FFVar WL_CVAL  = Cw - Cl;
  FFVar WL_CFLUX = p.DmwC*OpP(Cw,rw) - p.DlC*OpP(Cl,rl);
  FFVar WL_SVAL  = Sw - Sl;
  FFVar WL_SFLUX = p.DmwS*OpP(Sw,rw) - p.DlS*OpP(Sl,rl);
  FFVar L_C_IN   = Cl - p.ClIn;
  FFVar L_S_IN   = Sl - p.SlIn;
  FFVar L_C_OUT  = OpP(Cl,z);
  FFVar L_S_OUT  = OpP(Sl,z);
  FFVar L_C_WALL = OpP(Cl,rl);
  FFVar L_S_WALL = OpP(Sl,rl);

  OCFESLV oc(&DAG);
  oc.add_domain( z,  FFDom(0., p.L, TEST_MDW_NEL_Z, FFDom::CGL, TEST_MDW_NZ) );
  oc.add_domain( rd, FFDom(0., 1.0, TEST_MDW_NEL_R, FFDom::CGL, TEST_MDW_NR) );
  oc.add_domain( rw, FFDom(0., 1.0, TEST_MDW_NEL_R, FFDom::CGL, TEST_MDW_NR) );
  oc.add_domain( rl, FFDom(0., 1.0, TEST_MDW_NEL_R, FFDom::CGL, TEST_MDW_NR) );
  oc.add_state( Cg, {z} );
  oc.add_state( Vg, {z} );
  oc.add_state( Cd, {rd,z} );
  oc.add_state( Cw, {rw,z} );
  oc.add_state( Sw, {rw,z} );
  oc.add_state( Cl, {rl,z} );
  oc.add_state( Sl, {rl,z} );

  oc.update_ref( Cg, [&]( OCFESLV::t_Coord const& c ){ return Cg_exact(c.at(z),p); } );
  oc.update_ref( Vg, [&]( OCFESLV::t_Coord const& c ){ return Vg_exact(c.at(z),p); } );
  oc.update_ref( Cd, [&]( OCFESLV::t_Coord const& c ){ return Cd_exact(c.at(rd),c.at(z),p); } );
  oc.update_ref( Cw, [&]( OCFESLV::t_Coord const& c ){ return Cw_exact(c.at(rw),c.at(z),p); } );
  oc.update_ref( Sw, [&]( OCFESLV::t_Coord const& c ){ return Sw_exact(c.at(rw),c.at(z),p); } );
  oc.update_ref( Cl, [&]( OCFESLV::t_Coord const& c ){ return Cl_exact(c.at(rl),c.at(z),p); } );
  oc.update_ref( Sl, [&]( OCFESLV::t_Coord const& c ){ return Sl_exact(c.at(rl),c.at(z),p); } );

  OCFESLV::EqnOptions bulk_opt( OCFESLV::EqnRole::INTERIOR  );//,  0, OCFESLV::Options::IC_AUTO, TEST_MDW_CLASSIFY!=0, true );
  OCFESLV::EqnOptions bnd_opt(  OCFESLV::EqnRole::BOUNDARY  );//,  0, OCFESLV::Options::IC_AUTO, false, false );
  OCFESLV::EqnOptions int_opt(  OCFESLV::EqnRole::INTERFACE );//, 0, OCFESLV::Options::IC_AUTO );
  int const Z_NO_UB = FFDom::ALL - FFDom::UB;
  int const Z_INT   = FFDom::ALL - FFDom::LB - FFDom::UB;
  int const R_INT   = FFDom::ALL - FFDom::LB - FFDom::UB;

  oc.add_equation( GAS_V,    {z,rd},    {Z_NO_UB,FFDom::LB},              bulk_opt );
  oc.add_equation( GAS_C,    {z,rd},    {Z_NO_UB,FFDom::LB},              bulk_opt );
  oc.add_equation( GAS_V_IN, {z},       {FFDom::UB},                      bnd_opt );
  oc.add_equation( GAS_C_IN, {z},       {FFDom::UB},                      bnd_opt );
  oc.add_equation( DRY_C,    {rd,z},    {R_INT,FFDom::ALL},               bulk_opt );
  oc.add_equation( WET_C,    {rw,z},    {R_INT,FFDom::ALL},               bulk_opt );
  oc.add_equation( WET_S,    {rw,z},    {R_INT,FFDom::ALL},               bulk_opt );
  oc.add_equation( LIQ_C,    {rl,z},    {R_INT,Z_INT},                    bulk_opt );
  oc.add_equation( LIQ_S,    {rl,z},    {R_INT,Z_INT},                    bulk_opt );
  oc.add_equation( GD_VAL,   {rd,z},    {FFDom::LB,FFDom::ALL},           int_opt );
  oc.add_equation( DW_VAL,   {rd,rw,z}, {FFDom::UB,FFDom::LB,FFDom::ALL}, int_opt );
  oc.add_equation( DW_FLUX,  {rd,rw,z}, {FFDom::UB,FFDom::LB,FFDom::ALL}, int_opt );
  oc.add_equation( DW_SFLUX, {rw,z},    {FFDom::LB,FFDom::ALL},           int_opt );
  oc.add_equation( WL_CVAL,  {rw,rl,z}, {FFDom::UB,FFDom::LB,FFDom::ALL}, int_opt );
  oc.add_equation( WL_CFLUX, {rw,rl,z}, {FFDom::UB,FFDom::LB,FFDom::ALL}, int_opt );
  oc.add_equation( WL_SVAL,  {rw,rl,z}, {FFDom::UB,FFDom::LB,FFDom::ALL}, int_opt );
  oc.add_equation( WL_SFLUX, {rw,rl,z}, {FFDom::UB,FFDom::LB,FFDom::ALL}, int_opt );
  oc.add_equation( L_C_WALL, {rl,z},    {FFDom::UB,FFDom::ALL},           bnd_opt );
  oc.add_equation( L_S_WALL, {rl,z},    {FFDom::UB,FFDom::ALL},           bnd_opt );
  oc.add_equation( L_C_IN,   {rl,z},    {R_INT,FFDom::LB},                bnd_opt );
  oc.add_equation( L_S_IN,   {rl,z},    {R_INT,FFDom::LB},                bnd_opt );
  oc.add_equation( L_C_OUT,  {rl,z},    {R_INT,FFDom::UB},                bnd_opt );
  oc.add_equation( L_S_OUT,  {rl,z},    {R_INT,FFDom::UB},                bnd_opt );

  //oc.set_evolution_domain(z);
  oc.reset_evolution_domain();
  oc.options.REDUCE.ORDER = TEST_MDW_REDUCE_ORDER ? OCFESLV::Options::RED_FULL : OCFESLV::Options::RED_NONE;
  oc.options.CLASSIFY.MODE     = TEST_MDW_CLASSIFY ? OCFESLV::Options::CLASS_AUTO : OCFESLV::Options::CLASS_NONE;
  oc.options.INTERFACE.TYPE = OCFESLV::Options::IC_AUTO;//VALUE;
  oc.options.INTERFACE.IMPOSITION = imp;
  oc.options.INTERFACE.SAT_SIGMA0 = TEST_MDW_SAT_SIGMA0;

  if( !oc.setup() ){
    std::cerr << "ERROR: OCFESLV::setup() failed\n";
    MR.ok = false;
    return MR;
  }
#if TEST_MDW_PRINT_OC
  std::cout << oc;
  //oc.display_interface_plan_diagnostics(std::cout);
#endif

  size_t const nVar=oc.n_colloc_sta();
  size_t const nEqn=oc.n_colloc_eqn();
  size_t const nTrace=oc.n_colloc_trace();
  MR.nVar = nVar;
  MR.nEqn = nEqn;
  MR.nTrace = nTrace;
  std::cout << "nVar=" << nVar << " nEqn=" << nEqn << " nTrace=" << nTrace
            << " square=" << (nVar==nEqn?"yes":"NO") << "\n";
  ok &= (nVar==nEqn);

#if TEST_MDW_CLASSIFY
  auto const& cls = oc.pde_type();
  auto const& sym = oc.symbol_cached();
  std::cout << "PDE type: " << OCFESLV::pde_type_name(cls.type)
            << "  At_singular=" << (cls.At_singular?"yes":"no")
            << "  parabolic=" << (cls.parabolic_structure_detected?"yes":"no")
            << "  evolution_hyperbolic=" << (cls.evolution_hyperbolic?"yes":"no") << "\n";
  std::cout << "Principal symbol size: states=" << sym.vState.size()
            << " equations=" << sym.vEqn.size()
            << " domains=" << sym.vDom.size() << "\n";
#endif

  std::vector<double> xExact;
  xExact.reserve(nVar);
  size_t n_aux=0;
  std::cout << "States after setup:";
  for( auto const& st : oc.states_colloc() ){
    std::cout << " " << st;
    std::string const nm=st.name();
    if( nm.find("D") == 0 ) ++n_aux;
    auto nodes=oc.node_colloc(st);
    for( auto const& xy : nodes ) xExact.push_back( exact_initial_value_for_state(nm,xy,p) );
  }
  std::cout << "\nAuxiliary states introduced: " << n_aux << "\n";
  size_t const nState = nVar >= nTrace ? nVar-nTrace : nVar;
  if( xExact.size() != nState ){
    std::cerr << "ERROR: exact vector size mismatch: ordinary states=" << xExact.size()
              << " expected=" << nState << " nVar=" << nVar << " nTrace=" << nTrace << "\n";
    MR.ok = false;
    return MR;
  }
  if( nTrace ){
    std::cout << "Trace/tau variables appended: " << nTrace << " (initialised to zero)\n";
    xExact.resize(nVar,0.0);
  }

  std::vector<double> res(nEqn,123456.0);
  if( !oc.eval(res.data(),nullptr,xExact.data(),nullptr,nullptr) ){
    std::cerr << "ERROR: exact residual evaluation failed\n";
    MR.ok = false;
    return MR;
  }
  size_t first_unwritten=nEqn;
  for( size_t i=0; i<nEqn; ++i ) if( res[i] == 123456.0 ){ first_unwritten=i; break; }
  ok &= check_close("all residual rows written before solve", first_unwritten==nEqn?0.0:1.0, 0.0);
  print_residuals("Reference-profile residual", res);
  print_variable_exact_errors( oc, xExact, p, "manufactured reference" );
  print_auxiliary_flux_exact_errors( oc, xExact, p, "manufactured reference" );
  print_duplicate_spreads( oc, xExact, "manufactured reference" );

  // Gate 1 (rank / well-posedness): smallest singular value of the full
  // Jacobian at the manufactured solution.  O(1) => well-posed; ~0 => structural
  // null (expected for the augmented modes, benign under QR min-norm).  Captured
  // per (epszL, mode) for the mode-aware verdict in main().  Costs a dense SVD,
  // so it is compiled only under -DTEST_MDW_RANK_GATE=1; otherwise sigma_min
  // stays at its sentinel and is reported as "--".
#if TEST_MDW_RANK_GATE
  {
    arma::mat Jx;
    if( dense_jacobian( oc, xExact, Jx ) ){
      arma::vec const sv = arma::svd( Jx );
      MR.sigma_min = sv.is_empty() ? 0.0 : sv( sv.n_elem-1 );
      std::cout << "Gate-1 sigma_min(J @ exact) = "
                << std::scientific << std::setprecision(3) << MR.sigma_min
                << "  (epszL=" << epszL_arg << ", " << strimp << ")\n";
    }
  }
#endif

  std::vector<double> x;
#if TEST_MDW_INIT_MODE == 0
  x = xExact;
  std::cout << "Initialisation: exact collocation profile\n";
#elif TEST_MDW_INIT_MODE == 1
  x = xExact;
  for( size_t i=0; i<x.size(); ++i ) x[i] += 1.0e-3*std::sin(0.37*double(i+1))*std::max(1.0,std::abs(x[i]));
  std::cout << "Initialisation: exact profile plus small perturbation\n";
#elif TEST_MDW_INIT_MODE == 2
  x.reserve(nVar);
  for( auto const& st : oc.states_colloc() ){
    std::string const nm=st.name();
    auto nodes=oc.node_colloc(st);
    for( size_t i=0; i<nodes.size(); ++i ) x.push_back( constant_initial_value_for_state(nm) );
  }
  if( x.size() != nState ){
    std::cerr << "ERROR: constant initial vector size mismatch: ordinary states=" << x.size()
              << " expected=" << nState << " nVar=" << nVar << " nTrace=" << nTrace << "\n";
    MR.ok = false;
    return MR;
  }
  x.resize(nVar,0.0);
  std::cout << "Initialisation: constant primitive fields, zero auxiliary/trace variables\n";
#else
#error "TEST_MDW_INIT_MODE must be 0, 1, or 2"
#endif

  if( !oc.eval(res.data(),nullptr,x.data(),nullptr,nullptr) ){
    std::cerr << "ERROR: initial residual evaluation failed\n";
    MR.ok = false;
    return MR;
  }
  print_residuals("Initial residual", res);

  // item 12: nonlinear solve via OCFESLV::solve() (equilibrated LM/Newton, sparse AD).
  oc.options.SOLVE.MAX_ITER = TEST_MDW_MAXIT;
  oc.options.SOLVE.RES_TOL  = TEST_MDW_SOLVE_TOL;
#if defined(CRONOS__WITH_SPQR) && defined(TEST_LAP_SPQR)
  // item 13: route the augmented (IC_TRACE/IC_STRONG) solve through sparse
  // rank-revealing QR instead of the JtJ normal equations.  Requires the header
  // built with -DCRONOS__WITH_SPQR and SuiteSparse linked.
  oc.options.SOLVE.FACTORIZATION = OCFESLV::Options::SOLVE_SPQR;
#endif

  OCFESLV::SolveReport const srep = oc.solve( x.data() );
  bool solved = srep.converged;
  if( !solved )
    std::cerr << "OCFESLV::solve did not converge: final|r|=" << srep.final_residual
              << " after " << srep.iterations << " it ("
              << srep.newton_steps << " Newton / " << srep.lm_steps << " LM)\n";
  MR.solved = solved;
  ok &= solved;
  if( !eval_residual(oc,x,res) ){
    MR.ok = false;
    return MR;
  }
  print_residuals("Final residual", res);
  print_largest_residuals("Final residual", res);
  MR.final_res = max_abs( res );
  StateExactErrors const varErr = print_variable_exact_errors( oc, x, p, "final "+strimp+" solution" );
  StateExactErrors const auxErr = print_auxiliary_flux_exact_errors( oc, x, p, "final "+strimp+" solution" );
  MR.eCg = varErr.eCg;
  MR.eVg = varErr.eVg;
  MR.eCd = varErr.eCd;
  MR.eCw = varErr.eCw;
  MR.eSw = varErr.eSw;
  MR.eCl = varErr.eCl;
  MR.eSl = varErr.eSl;
  MR.eAux = auxErr.eAux;
  MR.max_spread = print_duplicate_spreads( oc, x, "post solve" );
  ok &= check_close("final max residual", MR.final_res, TEST_MDW_SOLVE_TOL);

  std::vector<size_t> nnz(nEqn,0);
  ok &= oc.deriv(nnz.data(),nullptr);
  size_t nnzsum=0,maxrownnz=0;
  for( size_t v: nnz ){ nnzsum += v; maxrownnz=std::max(maxrownnz,v); }
  std::cout << "Jacobian sparsity: nnz=" << nnzsum << " max_row_nnz=" << maxrownnz << "\n";
  ok &= (nnzsum > 0);

  auto zNodes = physical_nodes( oc.var_domain().at(z) );
  auto rNodes = physical_nodes( oc.var_domain().at(rl) );

  double max_gv_src=0., max_gc_src=0., max_d_src=0., max_wc_src=0., max_ws_src=0., max_lc_src=0., max_ls_src=0.;
  for( double zz: zNodes ){
    max_gv_src=std::max(max_gv_src,std::abs(gas_v_src_exact(zz,p)));
    max_gc_src=std::max(max_gc_src,std::abs(gas_c_src_exact(zz,p)));
    for( double rr: rNodes ){
      max_d_src=std::max(max_d_src,std::abs(dry_c_src_exact(rr,zz,p)));
      max_wc_src=std::max(max_wc_src,std::abs(wet_c_src_exact(rr,zz,p)));
      max_ws_src=std::max(max_ws_src,std::abs(wet_s_src_exact(rr,zz,p)));
      max_lc_src=std::max(max_lc_src,std::abs(liq_c_src_exact(rr,zz,p)));
      max_ls_src=std::max(max_ls_src,std::abs(liq_s_src_exact(rr,zz,p)));
    }
  }
  std::cout << "Exact manufactured source/helper magnitudes:\n"
            << "  max|gas_v_src_exact|=" << std::scientific << std::setprecision(6) << max_gv_src << "\n"
            << "  max|gas_c_src_exact|=" << max_gc_src << "\n"
            << "  max|dry_c_src_exact|=" << max_d_src << "\n"
            << "  max|wet_c_src_exact|=" << max_wc_src << "\n"
            << "  max|wet_s_src_exact|=" << max_ws_src << "\n"
            << "  max|liq_c_src_exact|=" << max_lc_src << "\n"
            << "  max|liq_s_src_exact|=" << max_ls_src << "\n";

  double eCg=0., eVg=0., eCd=0., eCw=0., eSw=0., eCl=0., eSl=0.;
  double minSl=std::numeric_limits<double>::infinity();
  double maxSl=-std::numeric_limits<double>::infinity();
  double max_sl_interface_minus_wall=0.0;
  double gd=0., dwv=0., dwf=0., dws=0., wlv=0., wlf=0., wlsv=0., wlsf=0., wallc=0., walls=0., linC=0., linS=0., loutC=0., loutS=0.;

  for( double zz: zNodes ){
    std::map<FFVar,double,lt_FFVar> pg{{z,zz}};
    double Cgv=interp(oc,Cg,pg,x), Vgv=interp(oc,Vg,pg,x);
    eCg=std::max(eCg,std::abs(Cgv-Cg_exact(zz,p)));
    eVg=std::max(eVg,std::abs(Vgv-Vg_exact(zz,p)));
    std::map<FFVar,double,lt_FFVar> pdL{{rd,0.0},{z,zz}}, pdU{{rd,1.0},{z,zz}};
    std::map<FFVar,double,lt_FFVar> pwL{{rw,0.0},{z,zz}}, pwU{{rw,1.0},{z,zz}};
    std::map<FFVar,double,lt_FFVar> plL{{rl,0.0},{z,zz}}, plU{{rl,1.0},{z,zz}};
    gd=std::max(gd,std::abs(interp(oc,Cd,pdL,x)-Cgv));
    dwv=std::max(dwv,std::abs(interp(oc,Cw,pwL,x)-p.H*interp(oc,Cd,pdU,x)));
    dwf=std::max(dwf,std::abs(p.Dmd*deriv_interp(oc,Cd,rd,pdU,x)-p.DmwC*deriv_interp(oc,Cw,rw,pwL,x)));
    dws=std::max(dws,std::abs(deriv_interp(oc,Sw,rw,pwL,x)));
    wlv=std::max(wlv,std::abs(interp(oc,Cw,pwU,x)-interp(oc,Cl,plL,x)));
    wlf=std::max(wlf,std::abs(p.DmwC*deriv_interp(oc,Cw,rw,pwU,x)-p.DlC*deriv_interp(oc,Cl,rl,plL,x)));
    wlsv=std::max(wlsv,std::abs(interp(oc,Sw,pwU,x)-interp(oc,Sl,plL,x)));
    wlsf=std::max(wlsf,std::abs(p.DmwS*deriv_interp(oc,Sw,rw,pwU,x)-p.DlS*deriv_interp(oc,Sl,rl,plL,x)));
    wallc=std::max(wallc,std::abs(deriv_interp(oc,Cl,rl,plU,x)));
    walls=std::max(walls,std::abs(deriv_interp(oc,Sl,rl,plU,x)));
    max_sl_interface_minus_wall=std::max(max_sl_interface_minus_wall, interp(oc,Sl,plL,x)-interp(oc,Sl,plU,x));
    for( double rr: rNodes ){
      std::map<FFVar,double,lt_FFVar> pd{{rd,rr},{z,zz}}, pw{{rw,rr},{z,zz}}, pl{{rl,rr},{z,zz}};
      eCd=std::max(eCd,std::abs(interp(oc,Cd,pd,x)-Cd_exact(rr,zz,p)));
      eCw=std::max(eCw,std::abs(interp(oc,Cw,pw,x)-Cw_exact(rr,zz,p)));
      eSw=std::max(eSw,std::abs(interp(oc,Sw,pw,x)-Sw_exact(rr,zz,p)));
      double const Clv=interp(oc,Cl,pl,x), Slv=interp(oc,Sl,pl,x);
      eCl=std::max(eCl,std::abs(Clv-Cl_exact(rr,zz,p)));
      eSl=std::max(eSl,std::abs(Slv-Sl_exact(rr,zz,p)));
      minSl=std::min(minSl,Slv); maxSl=std::max(maxSl,Slv);
      if( rr > 1e-12 && rr < 1.0-1e-12 ){
        if( zz < 1e-12 ){
          linC=std::max(linC,std::abs(Clv-p.ClIn));
          linS=std::max(linS,std::abs(Slv-p.SlIn));
        }
        if( std::abs(zz-p.L) < 1e-12 ){
          loutC=std::max(loutC,std::abs(deriv_interp(oc,Cl,z,pl,x)));
          loutS=std::max(loutS,std::abs(deriv_interp(oc,Sl,z,pl,x)));
        }
      }
    }
  }

  std::cout << "\nExact-solution comparison:\n";
  ok &= check_close("max |Cg-Cg_exact|", eCg, 3e-3);
  ok &= check_close("max |Vg-Vg_exact|", eVg, 2e-3);
  ok &= check_close("max |Cd-Cd_exact|", eCd, 3e-3);
  ok &= check_close("max |Cw-Cw_exact|", eCw, 2e-3);
  ok &= check_close("max |Sw-Sw_exact|", eSw, 2e-3);
  ok &= check_close("max |Cl-Cl_exact|", eCl, 2e-3);
  ok &= check_close("max |Sl-Sl_exact|", eSl, 2e-3);
  ok &= check_close("Sl remains positive", std::max(0.0,-minSl), 0.0);
  ok &= check_close("Sl lower at wet/liquid interface than wall", std::max(0.0,max_sl_interface_minus_wall), 2e-6);
  std::cout << "Sl computed range: [" << std::scientific << std::setprecision(6) << minSl << ", " << maxSl << "]\n";

  std::cout << "\nIndependent interface/boundary diagnostics:\n";
  ok &= check_close("gas/dry value", gd, 3e-6);
  ok &= check_close("dry/wet value", dwv, 3e-6);
  ok &= check_close("dry/wet CO2 flux", dwf, 3e-6);
  ok &= check_close("dry/wet solvent zero-flux", dws, 3e-6);
  ok &= check_close("wet/liquid CO2 value", wlv, 3e-6);
  ok &= check_close("wet/liquid CO2 flux", wlf, 3e-6);
  ok &= check_close("wet/liquid solvent value", wlsv, 3e-6);
  ok &= check_close("wet/liquid solvent flux", wlsf, 3e-6);
  ok &= check_close("liquid wall CO2 no-flux", wallc, 3e-6);
  ok &= check_close("liquid wall solvent no-flux", walls, 3e-6);
  ok &= check_close("liquid inlet CO2 interior", linC, 3e-6);
  ok &= check_close("liquid inlet solvent interior", linS, 3e-6);
  ok &= check_close("liquid outlet CO2 no-dispersion", loutC, 1e-5);
  ok &= check_close("liquid outlet solvent no-dispersion", loutS, 1e-5);

  std::string const prefix = std::string(TEST_MDW_OUTPUT_PREFIX) + "_" + suffix;
  {
    std::ofstream out(prefix + std::string("_gas.out"));
    out << "# z Cg Vg Cg_exact Vg_exact\n";
    for( double zz: zNodes ){
      std::map<FFVar,double,lt_FFVar> pt{{z,zz}};
      out << std::setprecision(16) << zz << " " << interp(oc,Cg,pt,x) << " " << interp(oc,Vg,pt,x)
          << " " << Cg_exact(zz,p) << " " << Vg_exact(zz,p) << "\n";
    }
  }
  {
    std::ofstream out(prefix + std::string("_dry.out"));
    out << "# z rd Cd Cd_exact\n";
    for( double zz: zNodes ){ for( double rr: rNodes ){
      std::map<FFVar,double,lt_FFVar> pt{{rd,rr},{z,zz}};
      out << std::setprecision(16) << zz << " " << rr << " " << interp(oc,Cd,pt,x) << " " << Cd_exact(rr,zz,p) << "\n";
    } out << "\n"; }
  }
  {
    std::ofstream out(prefix + std::string("_wet.out"));
    out << "# z rw Cw Sw Cw_exact Sw_exact\n";
    for( double zz: zNodes ){ for( double rr: rNodes ){
      std::map<FFVar,double,lt_FFVar> pt{{rw,rr},{z,zz}};
      out << std::setprecision(16) << zz << " " << rr << " " << interp(oc,Cw,pt,x) << " " << interp(oc,Sw,pt,x)
          << " " << Cw_exact(rr,zz,p) << " " << Sw_exact(rr,zz,p) << "\n";
    } out << "\n"; }
  }
  {
    std::ofstream out(prefix + std::string("_liquid.out"));
    out << "# z rl Cl Sl Cl_exact Sl_exact\n";
    for( double zz: zNodes ){ for( double rr: rNodes ){
      std::map<FFVar,double,lt_FFVar> pt{{rl,rr},{z,zz}};
      out << std::setprecision(16) << zz << " " << rr << " " << interp(oc,Cl,pt,x) << " " << interp(oc,Sl,pt,x)
          << " " << Cl_exact(rr,zz,p) << " " << Sl_exact(rr,zz,p) << "\n";
    } out << "\n"; }
  }
  std::cout << "\nWrote plot files: " << prefix << "_{gas,dry,wet,liquid}.out\n";
  std::cout << "\nManufactured dry/wet/liquid MBC test: " << (ok?"PASS":"FAIL") << "\n";
  MR.ok = ok;
  return MR;
}


int main()
{
  // Acceptance-gate sweep (plan §2 / §9 item 3): run all three imposition modes
  // across the epszL range so the gates are checked for loop-invariance (§2.5),
  // not just at the single default point.  epszL scales the z-second-derivative
  // coupling in the liquid equations, taking the interface from weakly to
  // strongly axially coupled.
  std::vector<double> const epszL_sweep = { 1.0e-2, 1.0, 1.5e2 };
  std::vector<ModeResult> results;
  for( double e : epszL_sweep ){
    std::cout << "\n################ epszL = " << std::scientific << std::setprecision(2)
              << e << " ################\n";
    results.push_back( run_mode( OCFESLV::Options::IC_WEAK,   "IC_WEAK",   "weak",   e ) );
    results.push_back( run_mode( OCFESLV::Options::IC_TRACE,  "IC_TRACE",  "trace",  e ) );
    results.push_back( run_mode( OCFESLV::Options::IC_STRONG, "IC_STRONG", "strong", e ) );
  }

  std::cout << "\n============== PDE3 acceptance-gate sweep (epszL x mode) ==============\n";
  std::cout << std::left  << std::setw(11) << "epszL"
            << std::setw(10) << "mode"
            << std::right << std::setw(8)  << "nTrace"
            << std::setw(13) << "final|r|"
            << std::setw(13) << "sigma_min"
            << std::setw(13) << "|dAux|max"
            << std::setw(9)  << "result" << "\n";

  bool   all_ok        = true;
  double sig_min_unaug = std::numeric_limits<double>::infinity();
  bool   any_unaug     = false;
  for( auto const& r : results ){
    all_ok &= r.ok;
    if( r.nTrace==0 ){ sig_min_unaug = std::min( sig_min_unaug, r.sigma_min ); any_unaug = true; }
    std::cout << std::left  << std::setw(11) << std::scientific << std::setprecision(1) << r.epszL
              << std::setw(10) << r.name
              << std::right << std::setprecision(3)
              << std::setw(8)  << r.nTrace
              << std::setw(13) << r.final_res
#if TEST_MDW_RANK_GATE
              << std::setw(13) << r.sigma_min
#else
              << std::setw(13) << "--"
#endif
              << std::setw(13) << r.eAux
              << std::setw(9)  << ( r.ok ? "PASS" : "FAIL" ) << "\n";
  }
  std::cout << "======================================================================\n";

#if TEST_MDW_RANK_GATE
  // Gate 1 (mode-aware): sigma_min is a true conditioning number only for the
  // non-augmented (nTrace==0) systems; reduce the pass/fail over those alone.
  // Augmented modes are rank-deficient by construction and listed above for
  // information only.  Reported, not folded into the return code; gates 2-5 are
  // the hard pass/fail.
  bool const rank_ok = any_unaug && sig_min_unaug > 1.0e-8;
  std::cout << "Gate 1 (rank, non-augmented modes): min sigma_min = "
            << std::scientific << std::setprecision(3) << ( any_unaug ? sig_min_unaug : 0.0 )
            << "  => " << ( rank_ok ? "OK (full rank, no soft mode)" : "WEAK (<1e-8 -- inspect)" ) << "\n";
#else
  (void)sig_min_unaug; (void)any_unaug;
  std::cout << "Gate 1 (rank): disabled -- build -DTEST_MDW_RANK_GATE=1 to enable"
               " (dense O(n^3) SVD per run)\n";
#endif
  std::cout << "Gates 2-5 (consistency / oracle-match / no-regression / loop-invariance): "
            << ( all_ok ? "PASS" : "FAIL" ) << "  over epszL in {1e-2, 1, 1.5e2}\n";
  std::cout << "PDE3 acceptance-gate sweep: " << ( all_ok ? "PASS" : "FAIL" ) << "\n";
  return all_ok ? 0 : 1;
}
