// Copyright (C) Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the EPL-2.0 with GPL-2.0-or-later as a Secondary License; see the LICENSE file.
// SPDX-License-Identifier: EPL-2.0 OR GPL-2.0-or-later
//
// OCFE_PDE2_solve.cpp
// -------------
// Simplified direct gas/liquid contactor. The dry and wet
// membrane compartments are omitted. The remaining compartments are coupled
// directly through the liquid-side boundary rl=LB:
//
//   gas tube      : Cg(z), Vg(z)
//   liquid shell  : Cl(rl,z), Sl(rl,z)
//
// The test uses a smooth manufactured reference profile to force the
// equations and validate the final solution.  By default, the nonlinear
// iterations start from constant primitive fields and zero derivative/trace
// auxiliary variables.  It then validates the square OCFESLV residual system,
// direct gas/liquid boundary equations, and the extra trace/tau variables
// introduced by IC_TRACE.
//
// Build-time switches:
//   Runtime sweep: IC_WEAK, IC_TRACE, and IC_STRONG are run sequentially.
//   -DTEST_GL_NEL_Z        number of axial finite elements
//   -DTEST_GL_NEL_R        number of liquid radial finite elements
//   -DTEST_GL_NZ           axial nodes per element
//   -DTEST_GL_NR           radial nodes per element
//   -DTEST_GL_INIT_MODE    0 exact, 1 perturbed exact, 2 constant primitives

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
#include OCFE_OCFESLV_HEADER

using namespace mc;

#ifndef TEST_GL_NEL_Z
#define TEST_GL_NEL_Z 3
#endif
#ifndef TEST_GL_NEL_R
#define TEST_GL_NEL_R 3
#endif
#ifndef TEST_GL_NZ
#define TEST_GL_NZ 7
#endif
#ifndef TEST_GL_NR
#define TEST_GL_NR 7
#endif
#ifndef TEST_GL_REDUCE_ORDER
#define TEST_GL_REDUCE_ORDER 1
#endif
#ifndef TEST_GL_CLASSIFY
#define TEST_GL_CLASSIFY 1
#endif
// 0 = forced IC_VALUE (default, current behaviour);
// 1 = IC_AUTO (let the resolver derive the interface type per direction).
#ifndef TEST_GL_INTERFACE_AUTO
#define TEST_GL_INTERFACE_AUTO 1
#endif
#ifndef TEST_GL_OUTPUT_PREFIX
#define TEST_GL_OUTPUT_PREFIX "OCFE_PDE2_solve"
#endif
#ifndef TEST_GL_MAXIT
#define TEST_GL_MAXIT 10
#endif
#ifndef TEST_GL_SAT_SIGMA0
#define TEST_GL_SAT_SIGMA0 1.0
#endif
#ifndef TEST_GL_PRINT_OC
#define TEST_GL_PRINT_OC 1
#endif
#ifndef TEST_GL_SOLVE_TOL
#define TEST_GL_SOLVE_TOL 1e-9
#endif
#ifndef TEST_GL_INIT_MODE
#define TEST_GL_INIT_MODE 2
#endif

#define TEST_GL_SPQR

namespace {

struct Par
{
  double L      = 1.0;
  double Hgl    = 0.12;    // direct gas/liquid equilibrium Cl = Hgl*Cg at rl=0
  double DlC    = 1.6e-2;
  double DlS    = 7.0e-3;
  double epszL  = 1.0e-2;
  double krl    = 0.42;
  double nuS    = 1.0;
  double Kvg    = 0.25;
  double Kcg    = 0.65;
  double betaG  = 0.10;
  double CgIn   = 1.0;     // gas inlet at z=L
  double VgIn   = 1.0;
  double SlIn   = 1.0;
  double Acg    = 0.12;    // exact gas concentration drop from z=L to z=0
  double Avg    = 0.05;    // exact gas velocity drop from z=L to z=0
  double Ar     = 0.025;   // radial liquid concentration curvature amplitude

  double ClIn() const { return Hgl * ( CgIn - Acg ); }
};

static double sqr( double x ){ return x*x; }
static double Rr( double r ){ return r * ( 2.0 - r ); }
static double Rr_r( double r ){ return 2.0 - 2.0*r; }
static double Rr_rr(){ return -2.0; }
static double Bz( double z, Par const& p ){ return p.Ar * z*z * sqr(1.0-z); }
static double Bz_z( double z, Par const& p ){ return p.Ar * ( 2.0*z - 6.0*z*z + 4.0*z*z*z ); }
static double Bz_zz( double z, Par const& p ){ return p.Ar * ( 2.0 - 12.0*z + 12.0*z*z ); }

static double Cg_exact( double z, Par const& p ){ return p.CgIn - p.Acg * sqr(1.0-z); }
static double Cg_z_exact( double z, Par const& p ){ return 2.0 * p.Acg * (1.0-z); }
static double Cg_zz_exact( Par const& p ){ return -2.0 * p.Acg; }
static double Vg_exact( double z, Par const& p ){ return p.VgIn - p.Avg * sqr(1.0-z); }
static double Vg_z_exact( double z, Par const& p ){ return 2.0 * p.Avg * (1.0-z); }
static double Cl_exact( double r, double z, Par const& p )
{ return p.Hgl * Cg_exact(z,p) + Bz(z,p) * Rr(r); }
static double Cl_r_exact( double r, double z, Par const& p )
{ return Bz(z,p) * Rr_r(r); }
static double Cl_rr_exact( double, double z, Par const& p )
{ return Bz(z,p) * Rr_rr(); }
static double Cl_z_exact( double r, double z, Par const& p )
{ return p.Hgl * Cg_z_exact(z,p) + Bz_z(z,p) * Rr(r); }
static double Cl_zz_exact( double r, double z, Par const& p )
{ return p.Hgl * Cg_zz_exact(p) + Bz_zz(z,p) * Rr(r); }
// Non-constant manufactured solvent profile.  The previous constant Sl field
// made Dz_Sl and Drl_Sl trivially zero, so donor/acceptor errors in the solvent
// auxiliary rows could be hidden while the CO2 auxiliary rows failed.
static double Hs( double r ){ return 1.0 - 3.0*r*r + 2.0*r*r*r; }
static double Hs_r( double r ){ return -6.0*r + 6.0*r*r; }
static double Hs_rr( double r ){ return -6.0 + 12.0*r; }
static double Zs( double z ){ return 3.0*z*z - 2.0*z*z*z; }
static double Zs_z( double z ){ return 6.0*z - 6.0*z*z; }
static double Zs_zz( double z ){ return 6.0 - 12.0*z; }
static double Sl_exact( double r, double z, Par const& p )
{ return p.SlIn - 0.08 * Zs(z) * Hs(r); }
static double Sl_r_exact( double r, double z, Par const& )
{ return -0.08 * Zs(z) * Hs_r(r); }
static double Sl_z_exact( double r, double z, Par const& )
{ return -0.08 * Zs_z(z) * Hs(r); }
static double Sl_rr_exact( double r, double z, Par const& )
{ return -0.08 * Zs(z) * Hs_rr(r); }
static double Sl_zz_exact( double r, double z, Par const& )
{ return -0.08 * Zs_zz(z) * Hs(r); }

[[maybe_unused]] static double gas_v_res_exact( double z, Par const& p )
{ return Vg_z_exact(z,p) + p.Kvg*p.DlC*Cl_r_exact(0.0,z,p); }
[[maybe_unused]] static double gas_c_res_exact( double z, Par const& p )
{
  double const C = Cg_exact(z,p);
  double const V = Vg_exact(z,p);
  return Cg_z_exact(z,p) + p.Kcg*p.DlC*(1.0-p.betaG*C)/V*Cl_r_exact(0.0,z,p);
}
[[maybe_unused]] static double liq_c_res_exact( double r, double z, Par const& p )
{
  double const vl = 1.0 + 0.20*r;
  double const C  = Cl_exact(r,z,p);
  double const S  = Sl_exact(r,z,p);
  return vl*Cl_z_exact(r,z,p)
       - p.DlC*( Cl_rr_exact(r,z,p) + p.epszL*Cl_zz_exact(r,z,p) )
       + p.krl*C*S;
}
[[maybe_unused]] static double liq_s_res_exact( double r, double z, Par const& p )
{
  double const C  = Cl_exact(r,z,p);
  double const S  = Sl_exact(r,z,p);
  double const vl = 1.0 + 0.20*r;
  return vl*Sl_z_exact(r,z,p)
       - p.DlS*( Sl_rr_exact(r,z,p) + p.epszL*Sl_zz_exact(r,z,p) )
       + p.nuS*p.krl*C*S;
}

static double max_abs( std::vector<double> const& r )
{
  double m = 0.;
  for( double v: r ) m = std::max( m, std::abs(v) );
  return m;
}

static std::vector<double> physical_nodes( FFDom const& dom )
{
  std::vector<double> x;
  for( size_t ie=0; ie<dom.n_elem; ++ie ){
    auto xe = dom.lgnodes( dom.elem_lo(ie), dom.elem_up(ie) );
    x.insert( x.end(), xe.begin(), xe.end() );
  }
  return x;
}

[[maybe_unused]] static bool is_element_boundary( FFDom const& dom, double x, double tol=1e-10 )
{
  for( double b : dom.elem_bnd )
    if( std::abs( x - b ) <= tol * std::max( 1.0, std::abs(dom.up_dom-dom.lo_dom) ) )
      return true;
  return false;
}

static void print_residuals( std::string const& label, std::vector<double> const& r )
{
  double sum = 0.;
  for( double v: r ) sum += std::abs(v);
  std::cout << std::left << std::setw(44) << label
            << " n=" << std::setw(5) << r.size()
            << " max|r|=" << std::scientific << std::setprecision(4) << max_abs(r)
            << " mean|r|=" << ( r.empty()? 0.:sum/r.size() ) << "\n";
}

[[maybe_unused]] static void print_largest_residuals
( std::string const& label, std::vector<double> const& r, size_t nprint=12 )
{
  std::vector<size_t> idx( r.size() );
  std::iota( idx.begin(), idx.end(), size_t(0) );
  std::partial_sort( idx.begin(), idx.begin()+std::min(nprint,idx.size()), idx.end(),
    [&]( size_t a, size_t b ){ return std::abs(r[a]) > std::abs(r[b]); } );
  std::cout << label << " largest residual rows:\n";
  for( size_t k=0; k<std::min(nprint,idx.size()); ++k ){
    size_t const i = idx[k];
    std::cout << "  row " << std::setw(6) << i
              << "  r=" << std::scientific << std::setprecision(8) << r[i] << "\n";
  }
}

static bool check_close( std::string const& label, double val, double tol )
{
  bool ok = std::isfinite(val) && val <= tol;
  std::cout << std::left << std::setw(52) << label
            << " value=" << std::scientific << std::setprecision(6) << val
            << " tol=" << tol << "  " << ( ok?"PASS":"FAIL" ) << "\n";
  return ok;
}

static double interp
( OCFESLV const& oc, FFVar const& v, std::map<FFVar,double,lt_FFVar> const& pt,
  std::vector<double> const& var )
{
  return oc.eval_colloc<double>( v, pt, var.data(), nullptr, nullptr );
}

static double deriv_interp
( OCFESLV const& oc, FFVar const& v, FFVar const& dom,
  std::map<FFVar,double,lt_FFVar> const& pt,
  std::vector<double> const& var )
{
  double lo = oc.var_domain().at( dom ).lo_dom;
  double up = oc.var_domain().at( dom ).up_dom;
  double x  = pt.at(dom);
  double h  = 1e-6 * std::max(1.0, up-lo);
  std::map<FFVar,double,lt_FFVar> pm = pt, pp = pt;
  if( x - h < lo ){
    pp[dom] = x + h; pm[dom] = x;
    return ( interp(oc,v,pp,var) - interp(oc,v,pm,var) ) / h;
  }
  if( x + h > up ){
    pp[dom] = x; pm[dom] = x - h;
    return ( interp(oc,v,pp,var) - interp(oc,v,pm,var) ) / h;
  }
  pp[dom] = x + h; pm[dom] = x - h;
  return ( interp(oc,v,pp,var) - interp(oc,v,pm,var) ) / (2*h);
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
    key.reserve( nodes[i].size() );
    for( double c : nodes[i] )
      key.push_back( static_cast<long long>( std::llround( c * 1.0e12 ) ) );
    double const v = var[off+i];
    auto it = groups.find( key );
    if( it == groups.end() ) groups.emplace( std::move(key), Accum{v,v,1} );
    else{
      it->second.lo = std::min( it->second.lo, v );
      it->second.hi = std::max( it->second.hi, v );
      ++it->second.count;
    }
  }

  DuplicateSpread out;
  for( auto const& kv : groups ){
    auto const& g = kv.second;
    out.max_multiplicity = std::max( out.max_multiplicity, g.count );
    if( g.count >= 2 ) out.max_pair_spread = std::max( out.max_pair_spread, g.hi-g.lo );
    if( g.count >= 4 ) out.max_corner_spread = std::max( out.max_corner_spread, g.hi-g.lo );
  }
  return out;
}

static double print_duplicate_spreads
( OCFESLV const& oc, std::vector<double> const& var, std::string const& label )
{
  std::cout << "\nElement-interface duplicate-node spreads (" << label << "):\n";
  double max_pair = 0.0;
  size_t off = 0;
  for( auto const& st : oc.states_colloc() ){
    auto nodes = oc.node_colloc( st );
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

static bool eval_residual( OCFESLV& oc, std::vector<double> const& x, std::vector<double>& r )
{
  std::fill( r.begin(), r.end(), 0.0 );
  return oc.eval( r.data(), nullptr, x.data(), nullptr, nullptr );
}


static double initial_value_for_state( std::string const& nm, std::vector<double> const& xy, Par const& p )
{
  double const z = xy[0];
  double const r = xy.size() > 1 ? xy[1] : 0.0;
  if( nm == "Cg(z)" ) return Cg_exact(z,p);
  if( nm == "Vg(z)" ) return Vg_exact(z,p);
  if( nm == "Cl(rl,z)" ) return Cl_exact(r,z,p);
  if( nm == "Sl(rl,z)" ) return Sl_exact(r,z,p);
  if( nm.find("Cl") != std::string::npos ){
    if( nm.find("Dz_")  != std::string::npos ) return Cl_z_exact(r,z,p);
    if( nm.find("Drl_") != std::string::npos ) return Cl_r_exact(r,z,p);
  }
  if( nm.find("Sl") != std::string::npos ){
    if( nm.find("Dz_")  != std::string::npos ) return Sl_z_exact(r,z,p);
    if( nm.find("Drl_") != std::string::npos ) return Sl_r_exact(r,z,p);
  }
  return 0.0;
}

static double constant_initial_value_for_state( std::string const& nm )
{
  // Deliberately not an exact profile: primitive states are spatially uniform,
  // while reduced derivative states and generated trace/tau variables start at
  // zero.  This exercises the nonlinear iterations from a physically scaled but
  // non-solution initial point.
  if( nm == "Cg(z)" )      return 1.00;
  if( nm == "Vg(z)" )      return 1.00;
  if( nm == "Cl(rl,z)" )   return 0.12;
  if( nm == "Sl(rl,z)" )   return 1.00;
  return 0.0;
}

struct StateExactErrors
{
  double eCg  = 0.0;
  double eVg  = 0.0;
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
        double const ref = initial_value_for_state( nm, nodes[i], p );
        double const err = std::abs( var[off+i] - ref );
        maxerr = std::max( maxerr, err );
        meanerr += err;
      }
      meanerr /= nodes.empty()? 1.0 : double(nodes.size());
      if( nm == "Cg(z)" ) out.eCg = maxerr;
      else if( nm == "Vg(z)" ) out.eVg = maxerr;
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
        double const ref = initial_value_for_state( nm, nodes[i], p );
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
  double      eCl = 0.0;
  double      eSl = 0.0;
  double      eAux = 0.0;
};

static ModeResult run_mode( OCFESLV::Options::ImpositionType imp, std::string const& strimp, std::string const& suffix )
{
  ModeResult R;
  R.name = strimp;
  bool ok = true;
  Par const p;


  std::cout << "\n========== simplified direct gas/liquid MBC test ==========" << "\n";
  std::cout << "imposition: " << strimp
            << ", finite elements: z=" << TEST_GL_NEL_Z
            << " r=" << TEST_GL_NEL_R
            << ", nodes/element: z=" << TEST_GL_NZ
            << " r=" << TEST_GL_NR << "\n";

  FFGraph DAG;
  FFVar z  = DAG.add_var( "z" );
  FFVar rl = DAG.add_var( "rl" );
  FFVar Cg = DAG.add_var( "Cg(z)" );
  FFVar Vg = DAG.add_var( "Vg(z)" );
  FFVar Cl = DAG.add_var( "Cl(rl,z)" );
  FFVar Sl = DAG.add_var( "Sl(rl,z)" );

  FFPartial OpP;

  FFVar omz  = 1.0 - z;
  FFVar CgE  = p.CgIn - p.Acg * omz * omz;
  FFVar VgE  = p.VgIn - p.Avg * omz * omz;
  FFVar RE   = rl * ( 2.0 - rl );
  FFVar BE   = p.Ar * z * z * omz * omz;
  FFVar ClE  = p.Hgl * CgE + BE * RE;
  FFVar SlE  = p.SlIn;

  FFVar CgE_z  = 2.0 * p.Acg * omz;
  FFVar VgE_z  = 2.0 * p.Avg * omz;
  FFVar ClE_r  = BE * ( 2.0 - 2.0*rl );
  FFVar ClE_rr = -2.0 * BE;
  FFVar BE_z   = p.Ar * ( 2.0*z - 6.0*z*z + 4.0*z*z*z );
  FFVar BE_zz  = p.Ar * ( 2.0 - 12.0*z + 12.0*z*z );
  FFVar ClE_z  = p.Hgl * CgE_z + BE_z * RE;
  FFVar ClE_zz = p.Hgl * ( -2.0 * p.Acg ) + BE_zz * RE;

  // Keep the symbolic manufactured solvent profile consistent with the
  // numerical helper functions Sl_exact(), Sl_r_exact(), etc.
  FFVar HsE    = 1.0 - 3.0*rl*rl + 2.0*rl*rl*rl;
  FFVar HsE_r  = -6.0*rl + 6.0*rl*rl;
  FFVar HsE_rr = -6.0 + 12.0*rl;
  FFVar ZsE    = 3.0*z*z - 2.0*z*z*z;
  FFVar ZsE_z  = 6.0*z - 6.0*z*z;
  FFVar ZsE_zz = 6.0 - 12.0*z;
  SlE          = p.SlIn - 0.08 * ZsE * HsE;
  FFVar SlE_r  = -0.08 * ZsE * HsE_r;
  FFVar SlE_rr = -0.08 * ZsE * HsE_rr;
  FFVar SlE_z  = -0.08 * ZsE_z * HsE;
  FFVar SlE_zz = -0.08 * ZsE_zz * HsE;

  FFVar vl = 1.0 + 0.20*rl;
  FFVar Rl = p.krl * Cl * Sl;
  FFVar RlE = p.krl * ClE * SlE;

  FFVar FG_V = VgE_z + p.Kvg*p.DlC*ClE_r;
  FFVar FG_C = CgE_z + p.Kcg*p.DlC*(1.0-p.betaG*CgE)/VgE*ClE_r;
  FFVar FL_C = vl*ClE_z - p.DlC*( ClE_rr + p.epszL*ClE_zz ) + RlE;
  FFVar FL_S = vl*SlE_z - p.DlS*( SlE_rr + p.epszL*SlE_zz ) + p.nuS*RlE;

  FFVar GAS_V = OpP( Vg, z ) + p.Kvg * p.DlC * OpP( Cl, rl ) - FG_V;
  FFVar GAS_C = OpP( Cg, z ) + p.Kcg * p.DlC * ( 1.0 - p.betaG*Cg ) / Vg * OpP( Cl, rl ) - FG_C;
  FFVar LIQ_C = vl*OpP( Cl, z ) - p.DlC*( OpP( Cl, {rl,2} ) + p.epszL*OpP( Cl, {z,2} ) ) + Rl - FL_C;
  FFVar LIQ_S = vl*OpP( Sl, z ) - p.DlS*( OpP( Sl, {rl,2} ) + p.epszL*OpP( Sl, {z,2} ) ) + p.nuS*Rl - FL_S;

  FFVar GAS_C_IN = Cg - CgE;
  FFVar GAS_V_IN = Vg - VgE;
  FFVar GL_CVAL  = Cl - p.Hgl*Cg;
  FFVar GL_SFLUX = OpP(Sl,rl) - SlE_r;
  FFVar L_C_IN   = Cl - p.ClIn();
  FFVar L_S_IN   = Sl - p.SlIn;
  FFVar L_C_OUT  = OpP(Cl,z) - ClE_z;
  FFVar L_S_OUT  = OpP(Sl,z) - SlE_z;
  FFVar L_C_WALL = OpP(Cl,rl) - ClE_r;
  FFVar L_S_WALL = OpP(Sl,rl) - SlE_r;

  OCFESLV oc( &DAG );
  oc.add_domain( z,  FFDom( 0., p.L, TEST_GL_NEL_Z, FFDom::CGL, TEST_GL_NZ ) );
  oc.add_domain( rl, FFDom( 0., 1.0, TEST_GL_NEL_R, FFDom::CGL, TEST_GL_NR ) );
  oc.add_state( Cg, {z} );
  oc.add_state( Vg, {z} );
  oc.add_state( Cl, {rl,z} );
  oc.add_state( Sl, {rl,z} );

  oc.update_ref( Cg, [&]( OCFESLV::t_Coord const& c ){ return Cg_exact(c.at(z),p); } );
  oc.update_ref( Vg, [&]( OCFESLV::t_Coord const& c ){ return Vg_exact(c.at(z),p); } );
  oc.update_ref( Cl, [&]( OCFESLV::t_Coord const& c ){ return Cl_exact(c.at(rl),c.at(z),p); } );
  oc.update_ref( Sl, [&]( OCFESLV::t_Coord const& c ){ return Sl_exact(c.at(rl),c.at(z),p); } );

  OCFESLV::EqnOptions bulk_opt( OCFESLV::EqnRole::INTERIOR );//, 0, OCFESLV::Options::IC_AUTO, TEST_GL_CLASSIFY!=0, true );
  OCFESLV::EqnOptions bnd_opt ( OCFESLV::EqnRole::BOUNDARY );//, 0, OCFESLV::Options::IC_AUTO, false, false );

  int const Z_NO_UB = FFDom::ALL - FFDom::UB;
  int const Z_INT   = FFDom::ALL - FFDom::LB - FFDom::UB;
  int const R_INT   = FFDom::ALL - FFDom::LB - FFDom::UB;

  oc.add_equation( GAS_V,    {z,rl}, {Z_NO_UB, FFDom::LB}, bulk_opt );
  oc.add_equation( GAS_C,    {z,rl}, {Z_NO_UB, FFDom::LB}, bulk_opt );
  oc.add_equation( GAS_V_IN, {z},    {FFDom::UB}, bnd_opt );
  oc.add_equation( GAS_C_IN, {z},    {FFDom::UB}, bnd_opt );

  oc.add_equation( LIQ_C,    {rl,z}, {R_INT, Z_INT}, bulk_opt );
  oc.add_equation( LIQ_S,    {rl,z}, {R_INT, Z_INT}, bulk_opt );

  oc.add_equation( GL_CVAL,  {rl,z}, {FFDom::LB, FFDom::ALL}, bnd_opt );
  oc.add_equation( GL_SFLUX, {rl,z}, {FFDom::LB, FFDom::ALL}, bnd_opt );
  oc.add_equation( L_C_WALL, {rl,z}, {FFDom::UB, FFDom::ALL}, bnd_opt );
  oc.add_equation( L_S_WALL, {rl,z}, {FFDom::UB, FFDom::ALL}, bnd_opt );
  oc.add_equation( L_C_IN,   {rl,z}, {R_INT, FFDom::LB}, bnd_opt );
  oc.add_equation( L_S_IN,   {rl,z}, {R_INT, FFDom::LB}, bnd_opt );
  oc.add_equation( L_C_OUT,  {rl,z}, {R_INT, FFDom::UB}, bnd_opt );
  oc.add_equation( L_S_OUT,  {rl,z}, {R_INT, FFDom::UB}, bnd_opt );

  oc.reset_evolution_domain();
  //oc.set_evolution_domain( z );
  oc.options.REDUCE.ORDER    = TEST_GL_REDUCE_ORDER
                             ? OCFESLV::Options::RED_FULL : OCFESLV::Options::RED_NONE;
  oc.options.CLASSIFY.MODE        = TEST_GL_CLASSIFY
                             ? OCFESLV::Options::CLASS_AUTO : OCFESLV::Options::CLASS_NONE;
  oc.options.INTERFACE.TYPE  = TEST_GL_INTERFACE_AUTO
                             ? OCFESLV::Options::IC_AUTO : OCFESLV::Options::IC_VALUE;
  oc.options.INTERFACE.IMPOSITION = imp;
  oc.options.INTERFACE.SAT_SIGMA0      = TEST_GL_SAT_SIGMA0;
  oc.options.DISPLAY_LEVEL   = 2;

  if( !oc.setup() ){
    std::cerr << "ERROR: OCFESLV::setup() failed\n";
    R.ok = false;
    return R;
  }

#if TEST_GL_PRINT_OC
  std::cout << oc;
  // Interface-plan diagnostics are printed by OCFESLV at DISPLAY_LEVEL>0 in current headers.
#endif
  size_t const nVar = oc.n_colloc_sta();
  size_t const nEqn = oc.n_colloc_eqn();
  size_t const nTrace = oc.n_colloc_trace();
  R.nVar = nVar;
  R.nEqn = nEqn;
  R.nTrace = nTrace;
  std::cout << "nVar=" << nVar << " nEqn=" << nEqn << " nTrace=" << nTrace
            << " square=" << (nVar==nEqn?"yes":"NO") << "\n";
  ok &= ( nVar == nEqn );

#if TEST_GL_CLASSIFY
  auto const& cls = oc.pde_type();
  auto const& sym = oc.symbol_cached();
  std::cout << "PDE type: " << OCFESLV::pde_type_name( cls.type )
            << "  At_singular=" << (cls.At_singular?"yes":"no")
            << "  parabolic=" << (cls.parabolic_structure_detected?"yes":"no")
            << "  evolution_hyperbolic=" << (cls.evolution_hyperbolic?"yes":"no") << "\n";
  std::cout << "Principal symbol size: states=" << sym.vState.size()
            << " equations=" << sym.vEqn.size()
            << " domains=" << sym.vDom.size() << "\n";
#endif

  std::vector<double> xExact;
  xExact.reserve( nVar );
  size_t n_aux = 0;
  std::cout << "States after setup:";
  for( auto const& st : oc.states_colloc() ){
    std::cout << " " << st;
    std::string const nm = st.name();
    if( nm.find("D") == 0 ) ++n_aux;
    auto nodes = oc.node_colloc( st );
    for( auto const& xy : nodes ) xExact.push_back( initial_value_for_state( nm, xy, p ) );
  }
  std::cout << "\nAuxiliary states introduced: " << n_aux << "\n";

  size_t const nState = nVar >= nTrace ? nVar - nTrace : nVar;
  if( xExact.size() != nState ){
    std::cerr << "ERROR: exact vector size mismatch: ordinary states=" << xExact.size()
              << " expected=" << nState << " nVar=" << nVar << " nTrace=" << nTrace << "\n";
    R.ok = false;
    return R;
  }
  if( nTrace ){
    std::cout << "Trace/tau variables appended: " << nTrace
              << " (initialised to zero)\n";
    xExact.resize( nVar, 0.0 );
  }

  std::vector<double> res( nEqn, 123456.0 );
  if( !oc.eval( res.data(), nullptr, xExact.data(), nullptr, nullptr ) ){
    std::cerr << "ERROR: exact residual evaluation failed\n";
    R.ok = false;
    return R;
  }
  size_t first_unwritten = nEqn;
  for( size_t i=0; i<nEqn; ++i ) if( res[i] == 123456.0 ){ first_unwritten=i; break; }
  ok &= check_close( "all residual rows written before solve", first_unwritten==nEqn?0.0:1.0, 0.0 );
  print_residuals( "Reference-profile residual", res );
  print_variable_exact_errors( oc, xExact, p, "manufactured reference" );
  print_auxiliary_flux_exact_errors( oc, xExact, p, "manufactured reference" );

  std::vector<double> x = xExact;
#if TEST_GL_INIT_MODE == 0
  std::cout << "Initialisation: exact manufactured profile\n";
#elif TEST_GL_INIT_MODE == 1
  std::cout << "Initialisation: perturbed manufactured profile\n";
  for( size_t i=0; i<x.size(); ++i ){
    double const s = std::sin( 0.37 * double(i+1) );
    x[i] += 1.0e-3 * s * std::max( 1.0, std::abs(x[i]) );
  }
#elif TEST_GL_INIT_MODE == 2
  std::cout << "Initialisation: constant primitive fields, zero auxiliary/trace variables\n";
  x.clear();
  x.reserve( nVar );
  for( auto const& st : oc.states_colloc() ){
    std::string const nm = st.name();
    auto nodes = oc.node_colloc( st );
    for( size_t i=0; i<nodes.size(); ++i ) x.push_back( constant_initial_value_for_state( nm ) );
  }
  if( x.size() != nState ){
    std::cerr << "ERROR: constant initial vector size mismatch: ordinary states=" << x.size()
              << " expected=" << nState << " nVar=" << nVar << " nTrace=" << nTrace << "\n";
    R.ok = false;
    return R;
  }
  x.resize( nVar, 0.0 );
#else
#error "TEST_GL_INIT_MODE must be 0, 1, or 2"
#endif

  if( !oc.eval( res.data(), nullptr, x.data(), nullptr, nullptr ) ){
    std::cerr << "ERROR: initial residual evaluation failed\n";
    R.ok = false;
    return R;
  }
  print_residuals( "Initial residual", res );

  // item 12: nonlinear solve via OCFESLV::solve() (equilibrated sparse damped LM).
  oc.options.SOLVE.MAX_ITER = TEST_GL_MAXIT;
  oc.options.SOLVE.RES_TOL  = TEST_GL_SOLVE_TOL;
#if defined(CRONOS__WITH_SPQR) && defined(TEST_GL_SPQR)
  // item 13: route the augmented (IC_TRACE/IC_STRONG) solve through sparse
  // rank-revealing QR instead of the JtJ normal equations.  Requires the header
  // built with -DCRONOS__WITH_SPQR and SuiteSparse linked.
  oc.options.SOLVE.FACTORIZATION = OCFESLV::Options::SOLVE_SPQR;
#endif

  OCFESLV::SolveReport const srep = oc.solve( x.data() );
  bool solved = srep.converged;
  if( !solved )
    std::cerr << "OCFESLV::solve did not converge: final|r|=" << srep.final_residual
              << " after " << srep.iterations << " it\n";
  R.solved = solved;
  ok &= solved;
  if( !eval_residual( oc, x, res ) ){
    R.ok = false;
    return R;
  }
  print_residuals( "Final residual", res );
  R.final_res = max_abs( res );
  ok &= check_close( "final max residual", R.final_res, TEST_GL_SOLVE_TOL );
  StateExactErrors const varErr = print_variable_exact_errors( oc, x, p, "final "+strimp+" solution" );
  StateExactErrors const auxErr = print_auxiliary_flux_exact_errors( oc, x, p, "final "+strimp+" solution" );
  R.eCg = varErr.eCg;
  R.eVg = varErr.eVg;
  R.eCl = varErr.eCl;
  R.eSl = varErr.eSl;
  R.eAux = auxErr.eAux;

  std::vector<size_t> nnz( nEqn, 0 );
  ok &= oc.deriv( nnz.data(), nullptr );
  size_t nnzsum = 0, maxrownnz = 0;
  for( size_t v: nnz ){ nnzsum += v; maxrownnz = std::max(maxrownnz,v); }
  std::cout << "Jacobian sparsity: nnz=" << nnzsum << " max_row_nnz=" << maxrownnz << "\n";
  ok &= ( nnzsum > 0 );

  double max_gl_cval=0., max_gl_sflux=0., max_wall_c=0., max_wall_s=0.;
  double max_liq_in_c=0., max_liq_in_s=0., max_liq_out_c=0., max_liq_out_s=0.;

  auto zNodes = physical_nodes( oc.var_domain().at(z) );
  auto rNodes = physical_nodes( oc.var_domain().at(rl) );
  for( double zz : zNodes ){
    std::map<FFVar,double,lt_FFVar> pg {{z,zz}};
    double const Cgv = interp(oc,Cg,pg,x);
    std::map<FFVar,double,lt_FFVar> plL {{rl,0.0},{z,zz}}, plU {{rl,1.0},{z,zz}};
    max_gl_cval  = std::max(max_gl_cval,  std::abs(interp(oc,Cl,plL,x)-p.Hgl*Cgv));
    max_gl_sflux = std::max(max_gl_sflux, std::abs(deriv_interp(oc,Sl,rl,plL,x)));
    max_wall_c   = std::max(max_wall_c,   std::abs(deriv_interp(oc,Cl,rl,plU,x)-Cl_r_exact(1.0,zz,p)));
    max_wall_s   = std::max(max_wall_s,   std::abs(deriv_interp(oc,Sl,rl,plU,x)));


    for( double rr : rNodes ){
      std::map<FFVar,double,lt_FFVar> pl {{rl,rr},{z,zz}};
      double const Clv = interp(oc,Cl,pl,x);
      double const Slv = interp(oc,Sl,pl,x);
      if( rr > 1e-12 && rr < 1.0-1e-12 ){
        if( zz < 1e-12 ){
          max_liq_in_c = std::max(max_liq_in_c, std::abs(Clv-p.ClIn()));
          max_liq_in_s = std::max(max_liq_in_s, std::abs(Slv-p.SlIn));
        }
        if( std::abs(zz-p.L) < 1e-12 ){
          max_liq_out_c = std::max(max_liq_out_c, std::abs(deriv_interp(oc,Cl,z,pl,x)-Cl_z_exact(rr,zz,p)));
          max_liq_out_s = std::max(max_liq_out_s, std::abs(deriv_interp(oc,Sl,z,pl,x)));
        }
      }
    }
  }

  std::cout << "\nIndependent direct gas/liquid diagnostics:\n";
  ok &= check_close( "gas/liquid CO2 equilibrium", max_gl_cval, 2e-7 );
  ok &= check_close( "gas/liquid solvent zero-flux", max_gl_sflux, 2e-6 );
  ok &= check_close( "liquid wall CO2 no-flux", max_wall_c, 5e-7 );
  ok &= check_close( "liquid wall solvent no-flux", max_wall_s, 2e-6 );
  ok &= check_close( "liquid inlet CO2 interior", max_liq_in_c, 5e-7 );
  ok &= check_close( "liquid inlet solvent interior", max_liq_in_s, 5e-7 );
  ok &= check_close( "liquid outlet CO2 no-dispersion", max_liq_out_c, 5e-6 );
  ok &= check_close( "liquid outlet solvent no-dispersion", max_liq_out_s, 5e-6 );

  double max_Cg_err=0., max_Vg_err=0., max_Cl_err=0., max_Sl_err=0.;
  for( double zz : zNodes ){
    std::map<FFVar,double,lt_FFVar> pg {{z,zz}};
    max_Cg_err = std::max( max_Cg_err, std::abs( interp(oc,Cg,pg,x) - Cg_exact(zz,p) ) );
    max_Vg_err = std::max( max_Vg_err, std::abs( interp(oc,Vg,pg,x) - Vg_exact(zz,p) ) );
    for( double rr : rNodes ){
      std::map<FFVar,double,lt_FFVar> pl {{rl,rr},{z,zz}};
      max_Cl_err = std::max( max_Cl_err, std::abs( interp(oc,Cl,pl,x) - Cl_exact(rr,zz,p) ) );
      max_Sl_err = std::max( max_Sl_err, std::abs( interp(oc,Sl,pl,x) - Sl_exact(rr,zz,p) ) );
    }
  }

  std::cout << "\nExact-solution comparison:\n";
  ok &= check_close( "max |Cg-Cg_exact|", max_Cg_err, 2e-8 );
  ok &= check_close( "max |Vg-Vg_exact|", max_Vg_err, 2e-5 );
  ok &= check_close( "max |Cl-Cl_exact|", max_Cl_err, 2e-8 );
  ok &= check_close( "max |Sl-Sl_exact|", max_Sl_err, 2e-8 );

  R.max_spread = print_duplicate_spreads( oc, x, "post solve" );
  ok &= check_close( "max duplicate-node spread", R.max_spread, 2e-8 );

  std::string const prefix = std::string(TEST_GL_OUTPUT_PREFIX) + "_" + suffix;
  {
    std::ofstream out( prefix + std::string("_gas.out") );
    out << "# z Cg Vg Cg_exact Vg_exact\n";
    for( double zz : zNodes ){
      std::map<FFVar,double,lt_FFVar> pt {{z,zz}};
      out << std::setprecision(16) << zz << " " << interp(oc,Cg,pt,x)
          << " " << interp(oc,Vg,pt,x)
          << " " << Cg_exact(zz,p) << " " << Vg_exact(zz,p) << "\n";
    }
  }
  {
    std::ofstream out( prefix + std::string("_liquid.out") );
    out << "# z rl Cl Sl Cl_exact Sl_exact\n";
    for( double zz : zNodes ){
      for( double rr : rNodes ){
        std::map<FFVar,double,lt_FFVar> pt {{rl,rr},{z,zz}};
        out << std::setprecision(16) << zz << " " << rr << " " << interp(oc,Cl,pt,x)
            << " " << interp(oc,Sl,pt,x)
            << " " << Cl_exact(rr,zz,p) << " " << Sl_exact(rr,zz,p) << "\n";
      }
      out << "\n";
    }
  }
  std::cout << "\nWrote plot files: " << prefix << "_gas.out, " << prefix << "_liquid.out\n";

  std::cout << "\nSimplified direct gas/liquid MBC test: " << (ok?"PASS":"FAIL") << "\n";
  R.ok = ok;
  return R;
}


int main()
{
  std::vector<ModeResult> results;
  results.push_back( run_mode( OCFESLV::Options::IC_WEAK,   "IC_WEAK",   "weak"   ) );
  results.push_back( run_mode( OCFESLV::Options::IC_TRACE,  "IC_TRACE",  "trace"  ) );
  results.push_back( run_mode( OCFESLV::Options::IC_STRONG, "IC_STRONG", "strong" ) );

  std::cout << "\n==================== PDE2 mode sweep summary ====================\n";
  std::cout << std::left << std::setw(10) << "mode"
            << std::right << std::setw(8)  << "nTrace"
            << std::setw(14) << "final|r|"
            << std::setw(14) << "dup-spread"
            << std::setw(12) << "|dCg|"
            << std::setw(12) << "|dVg|"
            << std::setw(12) << "|dCl|"
            << std::setw(12) << "|dSl|"
            << std::setw(12) << "|dAux|max"
            << std::setw(8)  << "result" << "\n";

  bool all_ok = true;
  for( auto const& r : results ){
    all_ok &= r.ok;
    std::cout << std::left << std::setw(10) << r.name
              << std::right << std::setw(8) << r.nTrace
              << std::scientific << std::setprecision(3)
              << std::setw(14) << r.final_res
              << std::setw(14) << r.max_spread
              << std::setw(12) << r.eCg
              << std::setw(12) << r.eVg
              << std::setw(12) << r.eCl
              << std::setw(12) << r.eSl
              << std::setw(12) << r.eAux
              << std::setw(8)  << ( r.ok ? "PASS" : "FAIL" ) << "\n";
  }

  std::cout << "=================================================================\n";
  std::cout << "Simplified direct gas/liquid MBC sweep over IC_WEAK, IC_TRACE, IC_STRONG: "
            << ( all_ok ? "PASS" : "FAIL" ) << "\n";
  return all_ok ? 0 : 1;
}
