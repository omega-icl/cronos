// Copyright (C) Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the EPL-2.0 with GPL-2.0-or-later as a Secondary License; see the LICENSE file.
// SPDX-License-Identifier: EPL-2.0 OR GPL-2.0-or-later
//
// OCFE_PDE1_solve.cpp
// ---------
// Gas-phase tubular-reactor case study for OCFESLV with homogeneous reaction
// creating moles, hence increasing the axial superficial velocity.
//
// The model is a transient axially dispersed ideal-gas reactor for the
// exothermic reaction
//
//     A  ->  2 B
//
// represented through a gas-expansion velocity law
//
//     U(C,T) = U0 * (T/Tfeed) * (1 + epsV*(Cfeed-C)/Cfeed),   epsV > 0.
//
// This captures both thermal expansion and net mole creation as reactant A is
// consumed.  The balances use conservative axial convection:
//
//     C_t + d/dz( U(C,T) C ) - d/dz( Dax   C_z ) + r(C,T)      = F_C(t,z)
//     T_t + d/dz( U(C,T) T ) - d/dz( alpha T_z ) - beta r(C,T) = F_T(t,z)
//
// with r(C,T) = Da*C*exp(gamma*T/(1+epsT*T)).  F_C and F_T are manufactured
// source terms from a smooth exact profile.  The test therefore remains a
// well-defined nonlinear square residual system while exercising:
//   * two primitive states C,T;
//   * two reduced auxiliary diffusive fluxes;
//   * nonlinear reaction coupling;
//   * nonlinear velocity coupling from gas expansion/mole creation;
//   * conservative first-order convection involving dU/dz;
//   * IC_WEAK reduced-interface exact replacement.
//
// Build-time switches:
//   -DGAS_HARD=1                  use a stronger expansion/reaction case
//   Runtime sweep: IC_WEAK, IC_TRACE, and IC_STRONG are run sequentially.
//   -DSKIP_DERIV=1                skip derivative regression checks
//   -DFAIL_ON_INTERFACE_SPREAD=1  fail if reduced z-interface trace
//                                       spreads exceed the tolerance
//   -DOUTPUT_FILE=\"name.out\"     write t,z,C,T quadruplets to name.out
//
// This file is self-contained and does not require NLPSLV/IPOPT/SNOPT headers.
// It converges the square residual system with a damped Newton method using
// OCFESLV::deriv().  In a full MC++ build, solve_problem() can be replaced by
// the project IPOPT/SNOPT harness used in other case studies.

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <string>
#include <utility>
#include <vector>

#ifndef OCFE_OCFESLV_HEADER
#define OCFE_OCFESLV_HEADER "ocfeslv.hpp"
#endif

#define MC__OCFESLV_INTERFACE_DECISION_V2
#define MC__OCFESLV_STRONG_EXPLICIT_TAU
//#define MC__OCFESLV_STRONG_DUMMY_DIAG
#define MC__OCFESLV_AUGMENTED_DROP_PROBE

#include OCFE_OCFESLV_HEADER
#include "test_deriv_utils.hpp"

using namespace mc;

#ifndef GAS_HARD
#define GAS_HARD 0
#endif


#ifndef SKIP_DERIV
#define SKIP_DERIV 0
#endif

#ifndef FAIL_ON_INTERFACE_SPREAD
#define FAIL_ON_INTERFACE_SPREAD 0
#endif

#ifndef OUTPUT_FILE
#define OUTPUT_FILE "OCFE_PDE1_solve.out"
#endif

namespace {

struct Params
{
  double L      = 1.0;
  double tf     = 2.0;
  double U0     = 1.0;
  double Dax    = 2.0e-2;
  double alpha  = 1.0e-2;
  double Da     = 0.35;
  double gamma  = 3.0;
  double epsT   = 0.20;
  double beta   = 1.5;
  double epsV   = 0.80;   // gas expansion factor: A -> 2B gives epsV ~= y_A0
  double Cbase  = 1.00;
  double Tbase  = 0.25;
  double Ac     = 0.16;
  double At     = 0.08;
  double lamC   = 0.55;
  double lamT   = 0.35;

  double Cfeed() const { return Cbase + Ac; }
  double Tfeed() const { return Tbase; }
};

struct Disc
{
  size_t n_el_t = 3;
  size_t n_nd_t = 7;
  size_t n_el_z = 3;
  size_t n_nd_z = 7;
  FFDom::TYPE typ_t = FFDom::LGR;
  FFDom::TYPE typ_z = FFDom::CGL;
};

Params params()
{
  Params p;
#if GAS_HARD
  p.U0     = 1.2;
  p.Dax    = 7.5e-3;
  p.alpha  = 6.0e-3;
  p.Da     = 1.25;
  p.gamma  = 6.0;
  p.epsT   = 0.22;
  p.beta   = 3.5;
  p.epsV   = 1.20;
  p.Cbase  = 1.00;
  p.Tbase  = 0.28;
  p.Ac     = 0.18;
  p.At     = 0.10;
  p.lamC   = 0.75;
  p.lamT   = 0.50;
#endif
  return p;
}

Disc disc()
{
  Disc d;
#if GAS_HARD
  d.n_el_t = 3;
  d.n_nd_t = 4;
  d.n_el_z = 3;
  d.n_nd_z = 4;
#endif
  return d;
}

static double phi( double z )   { return ( 1.0 - z ) * ( 1.0 - z ); }
static double phiz( double z )  { return -2.0 * ( 1.0 - z ); }
static double phizz()           { return 2.0; }
static double theta( double z ) { return 1.0 - ( 1.0 - z ) * ( 1.0 - z ); }
static double thetaz( double z ){ return 2.0 * ( 1.0 - z ); }
static double thetazz()         { return -2.0; }

static double C_exact( double t, double z, Params const& p )
{ return p.Cbase + p.Ac * std::exp( -p.lamC * t ) * phi( z ); }

static double T_exact( double t, double z, Params const& p )
{ return p.Tbase + p.At * std::exp( -p.lamT * t ) * theta( z ); }

static double qC_exact( double t, double z, Params const& p )
{ return p.Dax * p.Ac * std::exp( -p.lamC * t ) * phiz( z ); }

static double qT_exact( double t, double z, Params const& p )
{ return p.alpha * p.At * std::exp( -p.lamT * t ) * thetaz( z ); }

static double U_exact( double t, double z, Params const& p )
{
  double const C = C_exact( t, z, p );
  double const T = T_exact( t, z, p );
  double const X = ( p.Cfeed() - C ) / p.Cfeed();
  return p.U0 * ( T / p.Tfeed() ) * ( 1.0 + p.epsV * X );
}

static double max_abs( std::vector<double> const& r )
{
  double m = 0.;
  for( double v : r ) m = std::max( m, std::abs( v ) );
  return m;
}

static void print_residuals( std::string const& label, std::vector<double> const& res )
{
  double sum = 0.;
  for( double v : res ) sum += std::abs( v );
  std::cout << std::left << std::setw(50) << label
            << " n=" << std::setw(5) << res.size()
            << " max|r|=" << std::scientific << std::setprecision(4)
            << std::setw(12) << max_abs( res )
            << " mean|r|=" << ( res.empty()? 0. : sum / res.size() ) << "\n";
}

static bool dense_solve( std::vector<double> A, std::vector<double> b,
                         std::vector<double>& x )
{
  size_t const n = b.size();
  x.assign( n, 0.0 );
  if( A.size() != n*n ) return false;

  double const eps = 1e-14;
  for( size_t k = 0; k < n; ++k ){
    size_t piv = k;
    double apiv = std::abs( A[k*n+k] );
    for( size_t i = k+1; i < n; ++i ){
      double const a = std::abs( A[i*n+k] );
      if( a > apiv ){ apiv = a; piv = i; }
    }
    if( apiv < eps ){
      std::cerr << "dense_solve: near-singular pivot at column " << k
                << " pivot=" << std::scientific << apiv << "\n";
      return false;
    }
    if( piv != k ){
      for( size_t j = k; j < n; ++j ) std::swap( A[k*n+j], A[piv*n+j] );
      std::swap( b[k], b[piv] );
    }

    double const akk = A[k*n+k];
    for( size_t i = k+1; i < n; ++i ){
      double const mult = A[i*n+k] / akk;
      A[i*n+k] = 0.0;
      for( size_t j = k+1; j < n; ++j ) A[i*n+j] -= mult * A[k*n+j];
      b[i] -= mult * b[k];
    }
  }

  for( size_t ii = 0; ii < n; ++ii ){
    size_t const i = n - 1 - ii;
    double s = b[i];
    for( size_t j = i+1; j < n; ++j ) s -= A[i*n+j] * x[j];
    x[i] = s / A[i*n+i];
  }
  return true;
}

struct SolveStats
{
  bool ok = false;
  size_t iter = 0;
  double initial_res = std::numeric_limits<double>::infinity();
  double final_res   = std::numeric_limits<double>::infinity();
};

static SolveStats solve_problem( OCFESLV& oc, std::vector<double>& var,
                                 std::vector<double> const& cval,
                                 std::string const& label )
{
  SolveStats st;
  size_t const n = oc.n_colloc_eqn();
  if( n != oc.n_colloc_sta() ){
    std::cerr << label << ": residual system is not square: nEqn=" << n
              << " nVar=" << oc.n_colloc_sta() << "\n";
    return st;
  }

  std::vector<size_t> nnz( n, 0 );
  if( !oc.deriv( nullptr, nullptr, nnz.data(), nullptr,
                 nullptr, nullptr, nullptr, nullptr,
                 nullptr, nullptr, cval.data() ) ){
    std::cerr << label << ": derivative sparsity query failed\n";
    return st;
  }

  std::vector<std::vector<size_t>> col( n );
  std::vector<std::vector<double>> grad( n );
  std::vector<size_t*> colptr( n, nullptr );
  std::vector<double*> gradptr( n, nullptr );
  for( size_t i = 0; i < n; ++i ){
    col[i].resize( nnz[i] );
    grad[i].resize( nnz[i] );
    colptr[i] = col[i].data();
    gradptr[i] = grad[i].data();
  }

  std::vector<double> res( n, 0.0 ), trial_res( n, 0.0 );
  if( !oc.eval( res.data(), nullptr, var.data(), nullptr, cval.data() ) ){
    std::cerr << label << ": initial residual evaluation failed\n";
    return st;
  }
  st.initial_res = max_abs( res );
  std::cout << "\n" << label << " damped Newton solve\n";
  print_residuals( "initial residual", res );

  double norm = st.initial_res;
  double const tol = 1e-10;
  size_t const maxit = 14;
  if( norm <= tol ){
    st.ok = true;
    st.final_res = norm;
    return st;
  }

  for( size_t it = 0; it < maxit; ++it ){
    if( !oc.deriv( res.data(), gradptr.data(), nnz.data(), colptr.data(),
                   nullptr, nullptr, nullptr, nullptr,
                   var.data(), nullptr, cval.data() ) ){
      std::cerr << label << ": derivative evaluation failed at iteration " << it << "\n";
      return st;
    }
    norm = max_abs( res );
    if( norm <= tol ){
      st.ok = true;
      st.iter = it;
      st.final_res = norm;
      return st;
    }

    std::vector<double> A( n*n, 0.0 );
    for( size_t i = 0; i < n; ++i ){
      for( size_t j = 0; j < col[i].size(); ++j ){
        if( col[i][j] >= n ){
          std::cerr << label << ": out-of-range Jacobian column " << col[i][j]
                    << " at row " << i << "\n";
          return st;
        }
        A[i*n + col[i][j]] += grad[i][j];
      }
    }

    std::vector<double> rhs( n );
    for( size_t i = 0; i < n; ++i ) rhs[i] = -res[i];
    std::vector<double> step;
    if( !dense_solve( std::move( A ), std::move( rhs ), step ) ){
      std::cerr << label << ": dense Newton solve failed at iteration " << it << "\n";
      return st;
    }

    double alpha = 1.0;
    std::vector<double> trial( var.size(), 0.0 );
    bool accepted = false;
    for( size_t ls = 0; ls < 16; ++ls ){
      for( size_t i = 0; i < var.size(); ++i ) trial[i] = var[i] + alpha * step[i];
      if( !oc.eval( trial_res.data(), nullptr, trial.data(), nullptr, cval.data() ) ){
        alpha *= 0.5;
        continue;
      }
      double const trial_norm = max_abs( trial_res );
      if( trial_norm < ( 1.0 - 1e-4 * alpha ) * norm || trial_norm < tol ){
        var.swap( trial );
        res.swap( trial_res );
        norm = trial_norm;
        accepted = true;
        break;
      }
      alpha *= 0.5;
    }

    std::cout << "  iter " << std::setw(2) << it+1
              << "  alpha=" << std::scientific << std::setprecision(3) << alpha
              << "  max|r|=" << norm << "\n";

    if( !accepted ){
      std::cerr << label << ": line search failed at iteration " << it << "\n";
      st.final_res = norm;
      return st;
    }
  }

  st.iter = maxit;
  st.final_res = norm;
  st.ok = ( norm <= tol );
  return st;
}

static std::vector<double> perturbed_initial_guess( std::vector<double> const& exact )
{
  std::vector<double> out( exact );
  for( size_t i = 0; i < out.size(); ++i ){
    double const s = std::sin( 0.31 * double( i + 1 ) );
    out[i] += 1.5e-2 * s * std::max( 1.0, std::abs( out[i] ) );
  }
  return out;
}

static size_t state_offset( OCFESLV const& oc, FFVar const& v )
{
  size_t off = 0;
  for( auto const& st : oc.states_colloc() ){
    auto nodes = oc.node_colloc( st );
    if( st.id() == v.id() ) return off;
    off += nodes.size();
  }
  return std::numeric_limits<size_t>::max();
}

static double z_interface_spread( OCFESLV const& oc, FFVar const& st,
                                  std::vector<double> const& var,
                                  Disc const& d )
{
  size_t const off = state_offset( oc, st );
  if( off == std::numeric_limits<size_t>::max() ) return std::numeric_limits<double>::infinity();

  auto nodes = oc.node_colloc( st );
  size_t const expected = d.n_el_z * d.n_el_t * d.n_nd_z * d.n_nd_t;
  if( nodes.size() != expected ){
    std::cerr << "z_interface_spread: unexpected node count for " << st
              << ": got " << nodes.size() << " expected " << expected << "\n";
    return std::numeric_limits<double>::infinity();
  }

  auto idx = [&]( size_t ez, size_t et, size_t iz, size_t it ){
    return off + ( ( ( ez * d.n_el_t + et ) * d.n_nd_z + iz ) * d.n_nd_t + it );
  };

  double maxspread = 0.0;
  for( size_t ez = 0; ez + 1 < d.n_el_z; ++ez ){
    for( size_t et = 0; et < d.n_el_t; ++et ){
      for( size_t it = 0; it < d.n_nd_t; ++it ){
        double const vl = var[ idx( ez,   et, d.n_nd_z-1, it ) ];
        double const vr = var[ idx( ez+1, et, 0,          it ) ];
        maxspread = std::max( maxspread, std::abs( vl - vr ) );
      }
    }
  }
  return maxspread;
}

static bool fill_exact_vector_with_fluxes( OCFESLV const& oc, FFVar const& C, FFVar const& T,
                                           std::vector<double>& var,
                                           Params const& p )
{
  size_t const nVar = oc.n_colloc_sta();
  size_t const nTrace = oc.n_colloc_trace();
  size_t const nState = nVar >= nTrace ? nVar - nTrace : nVar;
  var.assign( nVar, 0.0 );
  size_t off = 0;
  for( auto const& st : oc.states_colloc() ){
    auto nodes = oc.node_colloc( st );
    std::string const nm = st.name();
    for( size_t i=0; i<nodes.size(); ++i ){
      double const tt = nodes[i][0];
      double const zz = nodes[i][1];
      if( st.id() == C.id() )
        var[off+i] = C_exact( tt, zz, p );
      else if( st.id() == T.id() )
        var[off+i] = T_exact( tt, zz, p );
      else if( nm.find("_C") != std::string::npos || nm.find("C(") != std::string::npos )
        var[off+i] = qC_exact( tt, zz, p );
      else if( nm.find("_T") != std::string::npos || nm.find("T(") != std::string::npos )
        var[off+i] = qT_exact( tt, zz, p );
      else
        return false;
    }
    off += nodes.size();
  }
  if( off != nState ){
    std::cerr << "fill_exact_vector_with_fluxes: state count " << off
              << " != non-trace variable count " << nState << "\n";
    return false;
  }
  // Tau/trace variables, if any, are generated unknowns and are zero at the manufactured reference.
  return true;
}

struct PrimitiveExactErrors
{
  double eC = 0.0;
  double eT = 0.0;
};

struct FluxExactErrors
{
  double eQC = 0.0;
  double eQT = 0.0;
};

static PrimitiveExactErrors print_variable_exact_errors
( OCFESLV const& oc, FFVar const& C, FFVar const& T,
  std::vector<double> const& var, Params const& p, std::string const& label )
{
  PrimitiveExactErrors out;
  size_t off = 0;
  std::cout << "\nVariable comparison against C_exact/T_exact (" << label << "):\n";
  for( auto const& st : oc.states_colloc() ){
    auto nodes = oc.node_colloc( st );
    bool const is_C = st.id() == C.id();
    bool const is_T = st.id() == T.id();
    if( is_C || is_T ){
      double maxerr = 0.0, meanerr = 0.0;
      for( size_t i=0; i<nodes.size(); ++i ){
        double const ref = is_C ? C_exact( nodes[i][0], nodes[i][1], p )
                                : T_exact( nodes[i][0], nodes[i][1], p );
        double const err = std::abs( var[off+i] - ref );
        maxerr = std::max( maxerr, err );
        meanerr += err;
      }
      meanerr /= nodes.empty()? 1.0 : double(nodes.size());
      if( is_C ) out.eC = maxerr;
      else       out.eT = maxerr;
      std::cout << "  " << std::setw(18) << st.name()
                << " max|v-v_exact|=" << std::scientific << std::setprecision(4) << maxerr
                << " mean|v-v_exact|=" << meanerr << "\n";
    }
    off += nodes.size();
  }
  return out;
}

static FluxExactErrors print_flux_exact_errors( OCFESLV const& oc, FFVar const& C, FFVar const& T,
                                                std::vector<double> const& var,
                                                Params const& p,
                                                std::string const& label )
{
  FluxExactErrors out;
  size_t off = 0;
  std::cout << "\nAuxiliary flux comparison against qC_exact/qT_exact (" << label << "):\n";
  for( auto const& st : oc.states_colloc() ){
    auto nodes = oc.node_colloc( st );
    std::string const nm = st.name();
    bool const is_qC = !(st.id() == C.id()) && !(st.id() == T.id())
                   && ( nm.find("_C") != std::string::npos || nm.find("C(") != std::string::npos );
    bool const is_qT = !(st.id() == C.id()) && !(st.id() == T.id())
                   && ( nm.find("_T") != std::string::npos || nm.find("T(") != std::string::npos );
    if( is_qC || is_qT ){
      double maxerr = 0.0, meanerr = 0.0;
      for( size_t i=0; i<nodes.size(); ++i ){
        double const ref = is_qC ? qC_exact( nodes[i][0], nodes[i][1], p )
                                 : qT_exact( nodes[i][0], nodes[i][1], p );
        double const err = std::abs( var[off+i] - ref );
        maxerr = std::max( maxerr, err );
        meanerr += err;
      }
      meanerr /= nodes.empty()? 1.0 : double(nodes.size());
      if( is_qC ) out.eQC = std::max( out.eQC, maxerr );
      else        out.eQT = std::max( out.eQT, maxerr );
      std::cout << "  " << std::setw(18) << st.name()
                << " max|q-q_exact|=" << std::scientific << std::setprecision(4) << maxerr
                << " mean|q-q_exact|=" << meanerr << "\n";
    }
    off += nodes.size();
  }
  return out;
}

static bool write_solution_grid( OCFESLV const& oc, FFVar const& C, FFVar const& T,
                                 std::vector<double> const& var,
                                 std::vector<double> const& cval,
                                 std::string const& fname )
{
  auto nodes = oc.node_colloc( C );
  std::vector<double> tv, zv;
  tv.reserve( nodes.size() );
  zv.reserve( nodes.size() );
  for( auto const& p : nodes ){
    tv.push_back( p[0] );
    zv.push_back( p[1] );
  }
  auto unique_sort = []( std::vector<double>& v ){
    std::sort( v.begin(), v.end() );
    std::vector<double> u;
    for( double x : v ){
      if( u.empty() || std::abs( x - u.back() ) > 1e-11 ) u.push_back( x );
    }
    v.swap( u );
  };
  unique_sort( tv );
  unique_sort( zv );

  FFVar vt, vz;
  bool found_t = false, found_z = false;
  for( auto const& kv : oc.var_domain() ){
    std::string const nm = kv.first.name();
    if( nm == "t" ){ vt = kv.first; found_t = true; }
    if( nm == "z" ){ vz = kv.first; found_z = true; }
  }
  if( !found_t || !found_z ) return false;

  std::ofstream os( fname.c_str() );
  if( !os ) return false;
  os << std::scientific << std::setprecision(16);
  for( double tval : tv ){
    for( double zval : zv ){
      OCFESLV::t_Coord pt;
      pt[vt] = tval;
      pt[vz] = zval;
      double const Cv = oc.eval_colloc( C, pt, var, std::vector<double>(), cval );
      double const Tv = oc.eval_colloc( T, pt, var, std::vector<double>(), cval );
      os << tval << " " << zval << " " << Cv << " " << Tv << "\n";
    }
    os << "\n";
  }
  std::cout << "Wrote collocated C/T solution to " << fname
            << " (" << tv.size() << " t-blocks, " << zv.size()
            << " z-points per block).\n";
  return true;
}


static std::string suffixed_output_file( char const* fname, std::string const& suffix )
{
  std::string out( fname );
  std::string const tag = std::string("_") + suffix;
  std::string::size_type const slash = out.find_last_of( "/\\" );
  std::string::size_type const dot = out.find_last_of( '.' );
  if( dot != std::string::npos && ( slash == std::string::npos || dot > slash ) )
    out.insert( dot, tag );
  else
    out += tag;
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
  double      eC = 0.0;
  double      eT = 0.0;
  double      eQC = 0.0;
  double      eQT = 0.0;
};

static ModeResult run_mode( OCFESLV::Options::ImpositionType imp, std::string const& strimp, std::string const& suffix )
{
  ModeResult R;
  R.name = std::string("IC_") + strimp;
  bool ok = true;
  Params const p = params();
  Disc const d = disc();
  std::string const output_file = suffixed_output_file( OUTPUT_FILE, suffix );

  std::cout << "\n========== gas-phase expanding tubular reactor ==========" << "\n";
  std::cout << "imposition: IC_" << strimp << "\n";
  std::cout << "Case difficulty: " << ( GAS_HARD ? "HARD" : "SIMPLE" ) << "\n";
  std::cout << "Discretisation: t ELE=" << d.n_el_t << " NOD=" << d.n_nd_t
            << " LGR, z ELE=" << d.n_el_z << " NOD=" << d.n_nd_z << " CGL\n";
  std::cout << "Velocity law: U=U0*(T/Tfeed)*(1+epsV*(Cfeed-C)/Cfeed), epsV="
            << p.epsV << "\n";
  std::cout << "Representative exact velocity: U(t=0,z=0)=" << U_exact( 0., 0., p )
            << "  U(t=0,z=L)=" << U_exact( 0., p.L, p ) << "\n";

  FFGraph DAG;
  FFVar t      = DAG.add_var( "t" );
  FFVar z      = DAG.add_var( "z" );
  FFVar C      = DAG.add_var( "C(t,z)" );
  FFVar T      = DAG.add_var( "T(t,z)" );
  FFVar U0     = DAG.add_var( "U0" );
  FFVar Dax    = DAG.add_var( "Dax" );
  FFVar alpha  = DAG.add_var( "alpha" );
  FFVar Da     = DAG.add_var( "Da" );
  FFVar gamma  = DAG.add_var( "gamma" );
  FFVar epsT   = DAG.add_var( "epsT" );
  FFVar beta   = DAG.add_var( "beta" );
  FFVar epsV   = DAG.add_var( "epsV" );
  FFVar Cbase  = DAG.add_var( "Cbase" );
  FFVar Tbase  = DAG.add_var( "Tbase" );
  FFVar Ac     = DAG.add_var( "Ac" );
  FFVar At     = DAG.add_var( "At" );
  FFVar lamC   = DAG.add_var( "lamC" );
  FFVar lamT   = DAG.add_var( "lamT" );

  FFPartial OpP;

  FFVar ph    = ( 1. - z ) * ( 1. - z );
  FFVar ph_z  = -2. * ( 1. - z );
  FFVar ph_zz = 2.;
  FFVar th    = 1. - ph;
  FFVar th_z  = 2. * ( 1. - z );
  FFVar th_zz = -2.;
  FFVar Ec    = exp( -lamC * t );
  FFVar Et    = exp( -lamT * t );

  FFVar Cfeed = Cbase + Ac;
  FFVar Tfeed = Tbase;

  FFVar Cex  = Cbase + Ac * Ec * ph;
  FFVar Tex  = Tbase + At * Et * th;
  FFVar Cex_t  = -lamC * Ac * Ec * ph;
  FFVar Cex_z  = Ac * Ec * ph_z;
  FFVar Cex_zz = Ac * Ec * ph_zz;
  FFVar Tex_t  = -lamT * At * Et * th;
  FFVar Tex_z  = At * Et * th_z;
  FFVar Tex_zz = At * Et * th_zz;

  FFVar X     = ( Cfeed - C ) / Cfeed;
  FFVar Ugas  = U0 * ( T / Tfeed ) * ( 1. + epsV * X );
  FFVar Xex   = ( Cfeed - Cex ) / Cfeed;
  FFVar Uex   = U0 * ( Tex / Tfeed ) * ( 1. + epsV * Xex );
  FFVar Xex_z = -Cex_z / Cfeed;
  FFVar Uex_z = U0 / Tfeed * ( Tex_z * ( 1. + epsV * Xex ) + Tex * epsV * Xex_z );

  FFVar rate    = Da * C   * exp( gamma * T   / ( 1. + epsT * T   ) );
  FFVar rate_ex = Da * Cex * exp( gamma * Tex / ( 1. + epsT * Tex ) );

  // Manufactured conservative-convection derivatives.
  FFVar dUexCex_z = Uex_z * Cex + Uex * Cex_z;
  FFVar dUexTex_z = Uex_z * Tex + Uex * Tex_z;

  FFVar FC = Cex_t + dUexCex_z - Dax   * Cex_zz + rate_ex;
  FFVar FT = Tex_t + dUexTex_z - alpha * Tex_zz - beta * rate_ex;

  // Conservative convection, expanded as d(UY)/dz = U*Y_z + Y*U_z.
  // This avoids a nested Partial[z](nonlinear state expression), while still
  // retaining the dU/dz terms from mole creation and thermal expansion.
  FFVar Cz = OpP( C, z );
  FFVar Tz = OpP( T, z );
  FFVar Ugas_z = U0 / Tfeed * ( Tz * ( 1. + epsV * X ) - T * epsV * Cz / Cfeed );
  FFVar convC = Ugas * Cz + C * Ugas_z;
  FFVar convT = Ugas * Tz + T * Ugas_z;

  FFVar PDE_C = OpP( C, t ) + convC
              - OpP( Dax   * Cz, z ) + rate - FC;
  FFVar PDE_T = OpP( T, t ) + convT
              - OpP( alpha * Tz, z ) - beta * rate - FT;

  FFVar IC_C = C - ( Cbase + Ac * ph );
  FFVar IC_T = T - ( Tbase + At * th );

  // Exact Danckwerts-type inlet fluxes at z=0.
  FFVar Cex0  = Cbase + Ac * Ec;
  FFVar Tex0  = Tbase;
  FFVar Cez0  = -2. * Ac * Ec;
  FFVar Tez0  = 2. * At * Et;
  FFVar Xex0  = ( Cfeed - Cex0 ) / Cfeed;
  FFVar Uex0  = U0 * ( Tex0 / Tfeed ) * ( 1. + epsV * Xex0 );
  FFVar FinC  = Uex0 * Cex0 - Dax   * Cez0;
  FFVar FinT  = Uex0 * Tex0 - alpha * Tez0;

  FFVar BC_C_IN  = Ugas * C - Dax   * OpP( C, z ) - FinC;
  FFVar BC_T_IN  = Ugas * T - alpha * OpP( T, z ) - FinT;
  FFVar BC_C_OUT = OpP( C, z );
  FFVar BC_T_OUT = OpP( T, z );

  std::cout << "Original equations:\n";
  auto sg = DAG.subgraph( {PDE_C, PDE_T, IC_C, IC_T, BC_C_IN, BC_T_IN, BC_C_OUT, BC_T_OUT} );
  auto ss = FFExpr::subgraph( &DAG, sg );
  for( auto const& s : ss ) std::cout << "  " << s << "\n";

  std::vector<double> cval{
    p.U0, p.Dax, p.alpha, p.Da, p.gamma, p.epsT, p.beta, p.epsV,
    p.Cbase, p.Tbase, p.Ac, p.At, p.lamC, p.lamT
  };

  OCFESLV oc( &DAG );
  oc.add_domain( t, FFDom( 0., p.tf, d.n_el_t, d.typ_t, d.n_nd_t ) );
  oc.add_domain( z, FFDom( 0., p.L,  d.n_el_z, d.typ_z, d.n_nd_z ) );
  oc.add_state( C, {t,z} );
  oc.add_state( T, {t,z} );
  oc.update_ref( C, [&]( OCFESLV::t_Coord const& coord ){
    return C_exact( coord.at( t ), coord.at( z ), p );
  } );
  oc.update_ref( T, [&]( OCFESLV::t_Coord const& coord ){
    return T_exact( coord.at( t ), coord.at( z ), p );
  } );
  oc.set_constant( {U0,Dax,alpha,Da,gamma,epsT,beta,epsV,Cbase,Tbase,Ac,At,lamC,lamT}, cval );
  oc.set_evolution_domain( t );

  oc.add_equation( PDE_C, {t,z}, {FFDom::ALL-FFDom::LB, FFDom::ALL-FFDom::LB-FFDom::UB},
                   OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( PDE_T, {t,z}, {FFDom::ALL-FFDom::LB, FFDom::ALL-FFDom::LB-FFDom::UB},
                   OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( IC_C, {t,z}, {FFDom::LB, FFDom::ALL},
                   OCFESLV::EqnOptions( OCFESLV::EqnRole::INITIAL, 0 ) );
  oc.add_equation( IC_T, {t,z}, {FFDom::LB, FFDom::ALL},
                   OCFESLV::EqnOptions( OCFESLV::EqnRole::INITIAL, 0 ) );
  oc.add_equation( BC_C_IN, {t,z}, {FFDom::ALL-FFDom::LB, FFDom::LB},
                   OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );
  oc.add_equation( BC_T_IN, {t,z}, {FFDom::ALL-FFDom::LB, FFDom::LB},
                   OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );
  oc.add_equation( BC_C_OUT, {t,z}, {FFDom::ALL-FFDom::LB, FFDom::UB},
                   OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );
  oc.add_equation( BC_T_OUT, {t,z}, {FFDom::ALL-FFDom::LB, FFDom::UB},
                   OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );

  oc.options.REDUCE.ORDER     = OCFESLV::Options::RED_FULL;
  oc.options.CLASSIFY.MODE         = OCFESLV::Options::CLASS_AUTO;
  oc.options.REDUCE.ORDER     = OCFESLV::Options::RED_FULL;
  oc.options.CLASSIFY.MODE         = OCFESLV::Options::CLASS_AUTO;
  oc.options.INTERFACE.IMPOSITION  = imp;
  oc.options.INTERFACE.TYPE   = OCFESLV::Options::IC_AUTO;
  oc.options.INTERFACE.SAT_SIGMA0       = 1.0;
  oc.options.DISPLAY_LEVEL    = 0;

  if( !oc.setup() ){
    std::cerr << "ERROR: OCFESLV::setup() failed\n";
    R.ok = false;
    return R;
  }
  std::cout << oc;
  //oc.display_interface_plan_diagnostics(std::cout);

  size_t const nVar = oc.n_colloc_sta();
  size_t const nEqn = oc.n_colloc_eqn();
  size_t const nTrace = oc.n_colloc_trace();
  R.nVar = nVar;
  R.nEqn = nEqn;
  R.nTrace = nTrace;
  std::cout << "nVar=" << nVar << " nEqn=" << nEqn
            << " nTrace=" << nTrace
            << " square=" << ( nVar == nEqn ? "yes" : "NO" ) << "\n";
  ok &= ( nVar == nEqn );

  auto const& cls = oc.pde_type();
  std::cout << "PDE type: " << OCFESLV::pde_type_name( cls.type )
            << "  At singular: " << ( cls.At_singular ? "yes" : "no" ) << "\n";

  std::cout << "States after reduction:";
  for( auto const& st : oc.states_colloc() ) std::cout << " " << st;
  std::cout << "\n";

  std::vector<double> varInit, inpInit;
  if( !oc.init( varInit, inpInit, cval.data() ) ){
    std::cerr << "ERROR: OCFESLV::init() failed\n";
    R.ok = false;
    return R;
  }
  if( varInit.size() != nVar ){
    std::cerr << "ERROR: init state size mismatch\n";
    R.ok = false;
    return R;
  }

  std::vector<double> res( nEqn, 0.0 );
  oc.eval( res.data(), nullptr, varInit.data(), nullptr, cval.data() );
  print_residuals( "OCFESLV::init reference residual", res );

  std::vector<double> varExact;
  if( !fill_exact_vector_with_fluxes( oc, C, T, varExact, p ) ){
    std::cerr << "ERROR: could not build exact C/T/qC/qT vector\n";
    R.ok = false;
    return R;
  }
  std::fill( res.begin(), res.end(), 0.0 );
  oc.eval( res.data(), nullptr, varExact.data(), nullptr, cval.data() );
  print_residuals( "manual C/T/qC/qT exact residual", res );
  print_variable_exact_errors( oc, C, T, varExact, p, "manual exact vector" );
  print_flux_exact_errors( oc, C, T, varExact, p, "manual exact vector" );

#if !SKIP_DERIV
  mc_test::DerivCheckOptions deropt;
  deropt.max_columns = 24;
  deropt.abs_tol = 5e-4;
  deropt.rel_tol = 2e-3;
  deropt.missing_tol = 1e-5;
  deropt.res_tol = 1e-8;
  ok &= mc_test::check_oc_derivatives( oc, "reactor IC_"+strimp, varInit, nullptr, cval.data(), deropt );
#endif

  std::vector<double> varGuess = perturbed_initial_guess( varInit );
  SolveStats stat = solve_problem( oc, varGuess, cval, "IC_"+strimp+" IC_AUTO" );
  R.solved = stat.ok;
  ok &= stat.ok;
  std::cout << "IC_"+strimp+" solve: " << ( stat.ok ? "PASS" : "FAIL" )
            << "  initial=" << std::scientific << stat.initial_res
            << "  final=" << stat.final_res
            << "  it=" << stat.iter << "\n";

  std::fill( res.begin(), res.end(), 0.0 );
  oc.eval( res.data(), nullptr, varGuess.data(), nullptr, cval.data() );
  print_residuals( "final IC_"+strimp+" residual", res );
  R.final_res = max_abs( res );
  PrimitiveExactErrors const varErr = print_variable_exact_errors( oc, C, T, varGuess, p, "final IC_"+strimp+" solution" );
  R.eC = varErr.eC;
  R.eT = varErr.eT;
  FluxExactErrors const fluxErr = print_flux_exact_errors( oc, C, T, varGuess, p, "final IC_"+strimp+" solution" );
  R.eQC = fluxErr.eQC;
  R.eQT = fluxErr.eQT;

  std::cout << "\nReduced-direction z-interface trace spreads in final IC_"+strimp+" solution:\n";
  for( auto const& st : oc.states_colloc() ){
    double const sp = z_interface_spread( oc, st, varGuess, d );
    R.max_spread = std::max( R.max_spread, sp );
    std::cout << "  " << std::setw(18) << st.name()
              << " max trace spread=" << std::scientific << std::setprecision(4)
              << sp << "\n";
#if FAIL_ON_INTERFACE_SPREAD
    ok &= std::isfinite( sp ) && sp < 1e-7;
#endif
  }

  if( !write_solution_grid( oc, C, T, varGuess, cval, output_file.c_str() ) ){
    std::cerr << "ERROR: could not write " << output_file << "\n";
    ok = false;
  }

  std::cout << "\nGnuplot example:\n"
            << "  splot '" << output_file << "' using 1:2:3 with lines title 'C_A'\n"
            << "  splot '" << output_file << "' using 1:2:4 with lines title 'T'\n";

  R.ok = ok;
  return R;
}



int main()
{
  std::vector<ModeResult> results;
  results.push_back( run_mode( OCFESLV::Options::IC_WEAK,   "WEAK",   "weak"   ) );
  results.push_back( run_mode( OCFESLV::Options::IC_TRACE,  "TRACE",  "trace"  ) );
  results.push_back( run_mode( OCFESLV::Options::IC_STRONG, "STRONG", "strong" ) );

  std::cout << "\n==================== PDE1 mode sweep summary ====================\n";
  std::cout << std::left << std::setw(10) << "mode"
            << std::right << std::setw(8)  << "nTrace"
            << std::setw(14) << "final|r|"
            << std::setw(14) << "z-spread"
            << std::setw(12) << "|dC|"
            << std::setw(12) << "|dT|"
            << std::setw(12) << "|dqC|"
            << std::setw(12) << "|dqT|"
            << std::setw(8)  << "result" << "\n";

  bool all_ok = true;
  for( auto const& r : results ){
    all_ok &= r.ok;
    std::cout << std::left << std::setw(10) << r.name
              << std::right << std::setw(8) << r.nTrace
              << std::scientific << std::setprecision(3)
              << std::setw(14) << r.final_res
              << std::setw(14) << r.max_spread
              << std::setw(12) << r.eC
              << std::setw(12) << r.eT
              << std::setw(12) << r.eQC
              << std::setw(12) << r.eQT
              << std::setw(8)  << ( r.ok ? "PASS" : "FAIL" ) << "\n";
  }

  std::cout << "=================================================================\n";
  std::cout << "Gas-expanding tubular reactor sweep over IC_WEAK, IC_TRACE, IC_STRONG: "
            << ( all_ok ? "PASS" : "FAIL" ) << "\n";
  return all_ok ? 0 : 1;
}
