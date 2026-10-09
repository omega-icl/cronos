// Copyright (C) Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the EPL-2.0 with GPL-2.0-or-later as a Secondary License; see the LICENSE file.
// SPDX-License-Identifier: EPL-2.0 OR GPL-2.0-or-later
//
// OCFE_deepcopy.cpp
// -----------------
// Regression checks for OCFESLV::deep_copy_from(), copy construction and
// assignment after the OCBase/OCFESLV split and the setup-time interface-plan
// refactor.  The test deliberately uses multiple finite elements and
// reduce_order(), then repeats the checks for IC_WEAK, IC_STRONG and IC_TRACE:
//
//   * IC_WEAK exercises SAT assembly and the lazy _satCoupling rebuild in copies;
//   * IC_STRONG exercises frozen donor/acceptor row replacement;
//   * IC_TRACE exercises the appended trace/tau variables and residual rows.
//
// The copy checks are performed after the user FFGraph has gone out of scope,
// so the copied environments must be self-contained and use their private
// local DAGs, setup metadata, frozen interface plans, and rebuildable caches.

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#ifndef OCFE_OCFESLV_HEADER
#define OCFE_OCFESLV_HEADER "ocfeslv.hpp"
#endif
#include OCFE_OCFESLV_HEADER

using namespace mc;

namespace {

static double T_exact( double t, double x )
{
  return 1.0 + 0.5*t + 3.0*x + x*x;
}

static double Tx_exact( double /*t*/, double x )
{
  return 3.0 + 2.0*x;
}

static double a_exact( double /*t*/, double /*x*/ )
{
  // Chosen so that T_t - T_xx + a = 0 exactly for T_exact.
  return 1.5;
}

static std::pair<double,double> infer_tx
( std::vector<double> const& coord )
{
  if( coord.size() != 2 )
    throw std::runtime_error( "Expected two coordinates" );

  // t is in [2,3] and x is in [0,1], so this does not depend on an assumed
  // coordinate ordering returned by node_colloc().
  if( coord[0] > 1.5 && coord[1] <= 1.5 ) return { coord[0], coord[1] };
  if( coord[1] > 1.5 && coord[0] <= 1.5 ) return { coord[1], coord[0] };
  throw std::runtime_error( "Could not infer (t,x) from node coordinates" );
}

static void accumulate_error
( double& maxerr, double got, double expect )
{
  maxerr = std::max( maxerr, std::fabs( got - expect ) );
}

static bool check_bool
( std::string const& label, bool pass )
{
  std::cout << std::left << std::setw(66) << label
            << ( pass? "PASS": "FAIL" ) << "\n";
  return pass;
}

static bool check_error
( std::string const& label, double err, double tol )
{
  bool const pass = std::isfinite( err ) && err <= tol;
  std::cout << std::left << std::setw(66) << label
            << "maxerr=" << std::scientific << std::setprecision(3) << err
            << " tol=" << tol << "  " << ( pass? "PASS": "FAIL" ) << "\n";
  return pass;
}

static double max_abs
( std::vector<double> const& v )
{
  double m = 0.0;
  for( double x : v ) m = std::max( m, std::fabs( x ) );
  return m;
}

static double max_diff
( std::vector<double> const& a, std::vector<double> const& b )
{
  if( a.size() != b.size() ) return std::numeric_limits<double>::infinity();
  double m = 0.0;
  for( size_t i=0; i<a.size(); ++i )
    m = std::max( m, std::fabs( a[i] - b[i] ) );
  return m;
}

static double max_diff_vector
( std::vector<double> const& a, std::vector<double> const& b )
{
  return max_diff( a, b );
}

static double max_diff_nodes
( OCFESLV const& a, OCFESLV const& b )
{
  if( a.states_colloc().size() != b.states_colloc().size() )
    return std::numeric_limits<double>::infinity();

  double err = 0.0;
  for( size_t k=0; k<a.states_colloc().size(); ++k ){
    auto const na = a.node_colloc( a.states_colloc()[k] );
    auto const nb = b.node_colloc( b.states_colloc()[k] );
    if( na.size() != nb.size() )
      return std::numeric_limits<double>::infinity();
    for( size_t i=0; i<na.size(); ++i ){
      if( na[i].size() != nb[i].size() )
        return std::numeric_limits<double>::infinity();
      for( size_t j=0; j<na[i].size(); ++j )
        accumulate_error( err, na[i][j], nb[i][j] );
    }
  }
  return err;
}

static double max_diff_matrix
( arma::mat const& A, arma::mat const& B )
{
  if( A.n_rows != B.n_rows || A.n_cols != B.n_cols )
    return std::numeric_limits<double>::infinity();
  double err = 0.0;
  for( arma::uword i=0; i<A.n_rows; ++i )
    for( arma::uword j=0; j<A.n_cols; ++j )
      accumulate_error( err, A(i,j), B(i,j) );
  return err;
}

static double max_diff_cxvec
( arma::cx_vec const& a, arma::cx_vec const& b )
{
  if( a.n_elem != b.n_elem )
    return std::numeric_limits<double>::infinity();
  double err = 0.0;
  for( arma::uword i=0; i<a.n_elem; ++i ){
    err = std::max( err, std::abs( a(i).real() - b(i).real() ) );
    err = std::max( err, std::abs( a(i).imag() - b(i).imag() ) );
  }
  return err;
}

static double max_classification_diff
( OCFESLV const& a, OCFESLV const& b )
{
  auto const& ca = a.pde_type();
  auto const& cb = b.pde_type();
  double err = 0.0;
  auto flag = [&]( bool lhs, bool rhs ){
    if( lhs != rhs ) err = std::numeric_limits<double>::infinity();
  };

  if( ca.type != cb.type || ca.evolution_dom_idx != cb.evolution_dom_idx
   || ca.rank_evolution != cb.rank_evolution )
    return std::numeric_limits<double>::infinity();

  flag( ca.At_singular, cb.At_singular );
  flag( ca.descriptor, cb.descriptor );
  flag( ca.parabolic_structure_detected, cb.parabolic_structure_detected );
  flag( ca.evolution_hyperbolic, cb.evolution_hyperbolic );
  flag( ca.spatially_characteristic, cb.spatially_characteristic );
  flag( ca.weak_hyperbolic, cb.weak_hyperbolic );
  flag( ca.degenerate, cb.degenerate );
  flag( ca.evolution_user_supplied, cb.evolution_user_supplied );
  flag( ca.evolution_auto_detected, cb.evolution_auto_detected );
  if( !std::isfinite( err ) ) return err;

  accumulate_error( err, ca.sigma_min_evolution, cb.sigma_min_evolution );
  accumulate_error( err, ca.sigma_max_evolution, cb.sigma_max_evolution );
  if( std::isfinite( ca.cond_evolution ) && std::isfinite( cb.cond_evolution ) )
    accumulate_error( err, ca.cond_evolution, cb.cond_evolution );
  else if( std::isfinite( ca.cond_evolution ) != std::isfinite( cb.cond_evolution ) )
    return std::numeric_limits<double>::infinity();
  accumulate_error( err, ca.max_imag_eig, cb.max_imag_eig );
  accumulate_error( err, ca.max_eigvec_cond, cb.max_eigvec_cond );
  accumulate_error( err, ca.imag_tol, cb.imag_tol );

  if( ca.Ai.size() != cb.Ai.size() )
    return std::numeric_limits<double>::infinity();
  for( size_t i=0; i<ca.Ai.size(); ++i )
    err = std::max( err, max_diff_matrix( ca.Ai[i], cb.Ai[i] ) );
  err = std::max( err, max_diff_cxvec( ca.eigenvalues, cb.eigenvalues ) );
  if( ca.eigendata.size() != cb.eigendata.size() )
    return std::numeric_limits<double>::infinity();
  for( size_t i=0; i<ca.eigendata.size(); ++i ){
    err = std::max( err, max_diff_vector( ca.eigendata[i].first, cb.eigendata[i].first ) );
    err = std::max( err, max_diff_cxvec( ca.eigendata[i].second, cb.eigendata[i].second ) );
  }
  return err;
}

static double max_symbol_signature_diff
( OCFESLV const& a, OCFESLV const& b )
{
  auto const& sa = a.symbol_cached();
  auto const& sb = b.symbol_cached();
  double err = 0.0;

  // State/domain ids are user-visible variables copied by FFGraph::insert()
  // with stable ids.  Equation/coefficient entries are expression nodes in the
  // private local DAG; their numerical effect is covered by classification and
  // derivative checks, so only their dimensions are compared here.
  auto cmp_user_vars = [&]( std::vector<FFVar> const& va, std::vector<FFVar> const& vb ){
    if( va.size() != vb.size() ){
      err = std::numeric_limits<double>::infinity();
      return;
    }
    for( size_t i=0; i<va.size(); ++i ){
      if( va[i].id() != vb[i].id() )
        err = std::numeric_limits<double>::infinity();
    }
  };

  if( sa.vEqn.size() != sb.vEqn.size() || sa.vCoeff.size() != sb.vCoeff.size() )
    return std::numeric_limits<double>::infinity();
  cmp_user_vars( sa.vState, sb.vState );
  cmp_user_vars( sa.vDom, sb.vDom );
  if( !std::isfinite( err ) ) return err;
  for( size_t i=0; i<sa.vCoeff.size(); ++i ){
    if( sa.vCoeff[i].size() != sb.vCoeff[i].size() )
      return std::numeric_limits<double>::infinity();
  }
  return err;
}

static FFDom const* find_domain_by_name
( OCFESLV const& oc, std::string const& name )
{
  for( auto const& [var, dom] : oc.var_domain() )
    if( std::string( var.name() ) == name ) return &dom;
  return nullptr;
}

static double max_domain_diff
( OCFESLV const& a, OCFESLV const& b, std::string const& name )
{
  FFDom const* da = find_domain_by_name( a, name );
  FFDom const* db = find_domain_by_name( b, name );
  if( !da || !db ) return std::numeric_limits<double>::infinity();
  if( da->type != db->type || da->n_elem != db->n_elem || da->n_node != db->n_node )
    return std::numeric_limits<double>::infinity();
  double err = 0.0;
  accumulate_error( err, da->lo_dom, db->lo_dom );
  accumulate_error( err, da->up_dom, db->up_dom );
  accumulate_error( err, da->w_elem, db->w_elem );
  err = std::max( err, max_diff_vector( da->elem_bnd, db->elem_bnd ) );
  err = std::max( err, max_diff_vector( da->elem_len, db->elem_len ) );
  return err;
}

static double max_domain_diff_expected
( OCFESLV const& oc, std::string const& name,
  std::vector<double> const& elem_bnd, std::vector<double> const& elem_len )
{
  FFDom const* dom = find_domain_by_name( oc, name );
  if( !dom ) return std::numeric_limits<double>::infinity();
  double err = 0.0;
  err = std::max( err, max_diff_vector( dom->elem_bnd, elem_bnd ) );
  err = std::max( err, max_diff_vector( dom->elem_len, elem_len ) );
  return err;
}

struct RefCheck
{
  double err_state = 0.0;
  double err_aux   = 0.0;
  double err_input = 0.0;
  bool   found_state = false;
  bool   found_aux   = false;
  bool   found_input = false;
};

static bool is_primitive_state_name
( std::string const& nm )
{
  return nm == "T(t,x)";
}

static RefCheck check_initialized_references
( OCFESLV const& oc, std::vector<double> const& var, std::vector<double> const& inp )
{
  RefCheck chk;

  if( var.size() != oc.n_colloc_sta() )
    throw std::runtime_error( "state vector has wrong size" );
  if( inp.size() != oc.n_colloc_inp() )
    throw std::runtime_error( "input vector has wrong size" );

  size_t offset = 0;
  for( auto const& st : oc.states_colloc() ){
    auto const nodes = oc.node_colloc( st );
    std::string const nm = st.name();
    bool const primitive = is_primitive_state_name( nm );
    for( size_t i=0; i<nodes.size(); ++i ){
      auto const [tt,xx] = infer_tx( nodes[i] );
      if( primitive ){
        chk.found_state = true;
        accumulate_error( chk.err_state, var[offset+i], T_exact( tt, xx ) );
      }
      else{
        chk.found_aux = true;
        accumulate_error( chk.err_aux, var[offset+i], Tx_exact( tt, xx ) );
      }
    }
    offset += nodes.size();
  }

  size_t inp_offset = 0;
  for( auto const& [invar, domset] : oc.var_input() ){
    (void)domset;
    auto const input_nodes = oc.node_colloc( invar );
    for( size_t i=0; i<input_nodes.size(); ++i ){
      auto const [tt,xx] = infer_tx( input_nodes[i] );
      chk.found_input = true;
      accumulate_error( chk.err_input, inp[inp_offset+i], a_exact( tt, xx ) );
    }
    inp_offset += input_nodes.size();
  }

  if( inp_offset != inp.size() )
    throw std::runtime_error( "input node count does not match input vector size" );

  return chk;
}

static bool init_eval
( OCFESLV& oc, std::vector<double>& var, std::vector<double>& inp,
  std::vector<double>& eqn, std::vector<double>& fct )
{
  var.assign( oc.n_colloc_sta(), 0.0 );
  inp.assign( oc.n_colloc_inp(), 0.0 );
  eqn.assign( oc.n_colloc_eqn(), 0.0 );
  fct.assign( oc.n_colloc_fct(), 0.0 );

  if( !oc.init( var.empty()? nullptr: var.data(), inp.empty()? nullptr: inp.data(), nullptr ) )
    return false;
  if( !oc.eval( eqn.empty()? nullptr: eqn.data(), fct.empty()? nullptr: fct.data(),
                var.empty()? nullptr: var.data(), inp.empty()? nullptr: inp.data(), nullptr ) )
    return false;
  return true;
}

static bool eval_only
( OCFESLV& oc, std::vector<double> const& var, std::vector<double> const& inp,
  std::vector<double>& eqn, std::vector<double>& fct )
{
  eqn.assign( oc.n_colloc_eqn(), 0.0 );
  fct.assign( oc.n_colloc_fct(), 0.0 );
  return oc.eval( eqn.empty()? nullptr: eqn.data(), fct.empty()? nullptr: fct.data(),
                  var.empty()? nullptr: var.data(), inp.empty()? nullptr: inp.data(), nullptr );
}

static std::vector<double> perturbed_vector
( std::vector<double> const& x, double scale, double phase )
{
  std::vector<double> y( x );
  for( size_t i=0; i<y.size(); ++i ){
    double const s = std::sin( phase + 0.417 * double( i + 1 ) );
    y[i] += scale * s * std::max( 1.0, std::fabs( y[i] ) );
  }
  return y;
}

static bool perturb_one_trace_variable
( OCFESLV const& oc, std::vector<double>& var, double delta )
{
  if( !oc.n_colloc_trace() || var.size() != oc.n_colloc_sta() ) return false;
  size_t const first_trace = oc.n_colloc_sta() - oc.n_colloc_trace();
  if( first_trace >= var.size() ) return false;
  var[first_trace] += delta;
  return true;
}

struct DerivData
{
  std::vector<std::vector<size_t>> cols;
  std::vector<std::vector<double>> vals;
};

static DerivData deriv_eqn_data
( OCFESLV& oc, std::vector<double> const& var, std::vector<double> const& inp )
{
  DerivData out;
  std::vector<size_t> nnz( oc.n_colloc_eqn(), 0 );
  if( !oc.deriv( nullptr, nullptr, nnz.data(), nullptr,
                 nullptr, nullptr, nullptr, nullptr,
                 nullptr, nullptr, nullptr ) )
    throw std::runtime_error( "deriv() failed while querying equation nnz" );

  out.cols.resize( oc.n_colloc_eqn() );
  out.vals.resize( oc.n_colloc_eqn() );
  std::vector<size_t*> colptr( oc.n_colloc_eqn(), nullptr );
  std::vector<double*> valptr( oc.n_colloc_eqn(), nullptr );
  for( size_t k=0; k<out.cols.size(); ++k ){
    out.cols[k].resize( nnz[k] );
    out.vals[k].resize( nnz[k] );
    colptr[k] = out.cols[k].empty()? nullptr: out.cols[k].data();
    valptr[k] = out.vals[k].empty()? nullptr: out.vals[k].data();
  }

  if( !oc.deriv( nullptr, valptr.empty()? nullptr: valptr.data(), nnz.data(), colptr.empty()? nullptr: colptr.data(),
                 nullptr, nullptr, nullptr, nullptr,
                 var.empty()? nullptr: var.data(), inp.empty()? nullptr: inp.data(), nullptr ) )
    throw std::runtime_error( "deriv() failed while retrieving equation derivatives" );

  return out;
}

static DerivData deriv_fct_data
( OCFESLV& oc, std::vector<double> const& var, std::vector<double> const& inp )
{
  DerivData out;
  std::vector<size_t> nnz( oc.n_colloc_fct(), 0 );
  if( !oc.deriv( nullptr, nullptr, nullptr, nullptr,
                 nullptr, nullptr, nnz.data(), nullptr,
                 nullptr, nullptr, nullptr ) )
    throw std::runtime_error( "deriv() failed while querying output nnz" );

  out.cols.resize( oc.n_colloc_fct() );
  out.vals.resize( oc.n_colloc_fct() );
  std::vector<size_t*> colptr( oc.n_colloc_fct(), nullptr );
  std::vector<double*> valptr( oc.n_colloc_fct(), nullptr );
  for( size_t k=0; k<out.cols.size(); ++k ){
    out.cols[k].resize( nnz[k] );
    out.vals[k].resize( nnz[k] );
    colptr[k] = out.cols[k].empty()? nullptr: out.cols[k].data();
    valptr[k] = out.vals[k].empty()? nullptr: out.vals[k].data();
  }

  if( !oc.deriv( nullptr, nullptr, nullptr, nullptr,
                 nullptr, valptr.empty()? nullptr: valptr.data(), nnz.data(), colptr.empty()? nullptr: colptr.data(),
                 var.empty()? nullptr: var.data(), inp.empty()? nullptr: inp.data(), nullptr ) )
    throw std::runtime_error( "deriv() failed while retrieving output derivatives" );

  return out;
}

static bool has_columns_at_or_above
( DerivData const& d, size_t col0 )
{
  for( auto const& row : d.cols )
    for( size_t c : row )
      if( c >= col0 ) return true;
  return false;
}

static bool same_pattern
( DerivData const& a, DerivData const& b )
{
  return a.cols == b.cols;
}

static double max_diff_deriv_values
( DerivData const& a, DerivData const& b )
{
  if( a.vals.size() != b.vals.size() ) return std::numeric_limits<double>::infinity();
  double err = 0.0;
  for( size_t i=0; i<a.vals.size(); ++i )
    err = std::max( err, max_diff( a.vals[i], b.vals[i] ) );
  return err;
}

static char const* imposition_name
( OCFESLV::Options::ImpositionType imp )
{
  switch( imp ){
    case OCFESLV::Options::IC_WEAK:   return "IC_WEAK";
    case OCFESLV::Options::IC_STRONG: return "IC_STRONG";
    case OCFESLV::Options::IC_TRACE:  return "IC_TRACE";
  }
  return "UNKNOWN";
}

static std::unique_ptr<OCFESLV> make_source_after_user_dag_scope
( OCFESLV::Options::ImpositionType imp )
{
  std::unique_ptr<OCFESLV> oc;

  {
    FFGraph DAG;
    FFVar t = DAG.add_var( "t" );
    FFVar x = DAG.add_var( "x" );
    FFVar T = DAG.add_var( "T(t,x)" );
    FFVar a = DAG.add_var( "a(t,x)" );

    OCFESLV::t_Fun Tref = [t,x]( OCFESLV::t_Coord const& coord ) -> double {
      return T_exact( coord.at( t ), coord.at( x ) );
    };
    OCFESLV::t_Fun Aref = [t,x]( OCFESLV::t_Coord const& coord ) -> double {
      return a_exact( coord.at( t ), coord.at( x ) );
    };

    FFPartial  OpP;
    FFIntegral OpI;
    FFVar Tex = 1.0 + 0.5*t + 3.0*x + x*x;
    FFVar PDE = OpP( T, t ) - OpP( OpP( T, x ), x ) + a;
    FFVar Tx  = OpP( T, x );
    FFVar Qx  = OpI( T, x );

    oc.reset( new OCFESLV( &DAG ) );

    // Non-uniform domains exercise both new FFDom constructors:
    //   - t from explicit element boundaries,
    //   - x from a lower bound and element lengths.
    std::vector<double> const t_bnd{ 2.0, 2.40, 3.0 };
    std::vector<double> const x_len{ 0.30, 0.70 };

    oc->add_domain( t, FFDom( t_bnd, FFDom::LGL, 4 ), 2.25 );
    oc->add_domain( x, FFDom( 0.0, x_len, FFDom::CGL, 5 ), 0.40 );
    oc->add_state( T, {t,x}, Tref );
    oc->add_input( a, {t,x}, Aref );
    oc->add_equation( PDE, {t,x}, {FFDom::ALL-FFDom::LB, FFDom::ALL-FFDom::LB-FFDom::UB},
                      OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
    oc->add_equation( T - Tex, {t,x}, {FFDom::LB, FFDom::ALL},
                      OCFESLV::EqnOptions( OCFESLV::EqnRole::INITIAL, 0 ) );
    oc->add_equation( T - Tex, {t,x}, {FFDom::ALL-FFDom::LB, FFDom::LB},
                      OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );
    oc->add_equation( T - Tex, {t,x}, {FFDom::ALL-FFDom::LB, FFDom::UB},
                      OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );

    // Outputs cover point evaluation, distributed output grids, mixed fixed/grid
    // output evaluation, and integral output rows.  The copy test compares both
    // their values and their sparse derivative patterns/values.
    oc->add_output( T, {t,x}, {2.75, 0.35} );
    oc->add_output( T + a, {t,x}, {FFDom::UB, FFDom::ALL-FFDom::LB-FFDom::UB} );
    oc->add_output( Tx, {t}, {FFDom::ALL}, {x}, {0.35} );
    oc->add_output( Qx, {t}, {FFDom::ALL} );

    oc->options.REDUCE.ORDER    = OCFESLV::Options::RED_FULL;
    oc->options.CLASSIFY.MODE        = OCFESLV::Options::CLASS_AUTO;
    oc->options.INTERFACE.IMPOSITION = imp;
    oc->options.INTERFACE.SAT_SIGMA0      = 1.0;
    // rev105: this driver tests STRUCTURE, not solving.  It builds an OCFESLV with a specific
    // boundary-constructor t domain and a non-uniform x domain, then asserts that a deep
    // copy preserves the grid, the node coordinates, the derivative cache and the residual.
    //
    // Path B (CRONOS_HOIST_EVODETECT=1) infers the evolution domain before the localisation
    // collapse, so `t` gets COLLAPSED to a single element -- and the driver then compares
    // that collapsed grid against the multi-element one it constructed.  MEASURED:
    // max_domain_diff_expected( oc, "t", ... ) returns inf, n_colloc_nnz() and
    // n_deriv_colors() read 0, and init_eval fails outright.  Exit status went 0 -> 1 in the
    // full-corpus A/B, and CRONOS_HOIST_EVODETECT=0 passes while =1 fails.
    //
    // SOLVE_MARCHING=false is the exact per-driver opt-out: it gates the collapse itself
    // (the header states it at rev105:4588), so the grid stays as built regardless of how
    // the evolution domain was found or when.  A structural test should not be silently
    // re-gridded.
    oc->options.SOLVE.MARCHING  = false;

    if( !oc->setup() )
      throw std::runtime_error( std::string("setup() failed for ") + imposition_name( imp ) );
  }

  // The user FFGraph is now destroyed.  Any subsequent copy/eval/deriv must
  // rely on OCFESLV's private local DAG and copied setup data only.
  return oc;
}

static bool run_copy_suite
( OCFESLV::Options::ImpositionType imp )
{
  std::cout << "\n------------------------------------------------------------\n";
  std::cout << "  Deep-copy case: " << imposition_name( imp ) << "\n";
  std::cout << "------------------------------------------------------------\n";

  bool ok = true;
  auto oc = make_source_after_user_dag_scope( imp );

  std::vector<double> const t_bnd{ 2.0, 2.40, 3.0 };
  std::vector<double> const t_len{ 0.40, 0.60 };
  std::vector<double> const x_len{ 0.30, 0.70 };
  std::vector<double> const x_bnd{ 0.0, 0.30, 1.0 };

  std::cout << "source dimensions: nEqn=" << oc->n_colloc_eqn()
            << " nVar=" << oc->n_colloc_sta()
            << " nInp=" << oc->n_colloc_inp()
            << " nFct=" << oc->n_colloc_fct()
            << " nTrace=" << oc->n_colloc_trace() << "\n";

  ok = check_bool ( "source setup creates a private local DAG",
                    oc->dag() != nullptr ) && ok;
  ok = check_bool ( "source residual system is square",
                    oc->n_colloc_eqn() == oc->n_colloc_sta() ) && ok;
  ok = check_bool ( "trace variables appear only under IC_TRACE",
                    imp == OCFESLV::Options::IC_TRACE ? oc->n_colloc_trace() > 0
                                                    : oc->n_colloc_trace() == 0 ) && ok;
  ok = check_bool ( "source has point, distributed, mixed and integral outputs",
                    oc->n_colloc_fct() > 3
                 && oc->blk_fct(0).second == 1
                 && oc->blk_fct(1).second > 0
                 && oc->blk_fct(2).second > 0
                 && oc->blk_fct(3).second > 0
                 && oc->row_fct(0) == oc->n_colloc_eqn() ) && ok;
  ok = check_bool ( "source derivative cache reports nonzero count/colors",
                    oc->n_colloc_nnz() > 0 && oc->n_deriv_colors() > 0 ) && ok;
  ok = check_error( "source boundary-constructor t boundaries/lengths",
                    max_domain_diff_expected( *oc, "t", t_bnd, t_len ), 1e-14 ) && ok;
  ok = check_error( "source non-uniform x boundaries/lengths",
                    max_domain_diff_expected( *oc, "x", x_bnd, x_len ), 1e-14 ) && ok;

  std::vector<double> var, inp, eqn, fct;
  if( !init_eval( *oc, var, inp, eqn, fct ) ){
    std::cerr << "ERROR: init/eval failed for source " << imposition_name( imp ) << ".\n";
    return false;
  }

  RefCheck const src_ref = check_initialized_references( *oc, var, inp );
  double const expected_output = T_exact( 2.75, 0.35 );
  double err_out = fct.empty()? std::numeric_limits<double>::infinity()
                            : std::fabs( fct[0] - expected_output );

  ok = check_bool ( "source init found primitive state", src_ref.found_state ) && ok;
  ok = check_bool ( "source init found reduce_order auxiliary", src_ref.found_aux ) && ok;
  ok = check_bool ( "source init found distributed input", src_ref.found_input ) && ok;
  ok = check_error( "source init primitive state values", src_ref.err_state, 1e-10 ) && ok;
  ok = check_error( "source init distributed input values", src_ref.err_input, 1e-12 ) && ok;
  ok = check_error( "source init reduce_order auxiliary values", src_ref.err_aux, 1e-9 ) && ok;
  ok = check_error( "source residual at initialized reference", max_abs( eqn ), 1e-9 ) && ok;
  ok = check_error( "source output at initialized reference", err_out, 1e-10 ) && ok;

  DerivData const src_der  = deriv_eqn_data( *oc, var, inp );
  DerivData const src_fder = deriv_fct_data( *oc, var, inp );

  size_t const first_trace_col = oc->n_colloc_sta() - oc->n_colloc_trace();
  ok = check_bool( "source IC_TRACE equation Jacobian contains tau columns",
                   imp == OCFESLV::Options::IC_TRACE
                 ? has_columns_at_or_above( src_der, first_trace_col )
                 : true ) && ok;

  std::vector<double> pvar = perturbed_vector( var, 3.0e-4, 0.1 );
  std::vector<double> pinp = perturbed_vector( inp, 7.0e-4, 0.7 );
  std::vector<double> peqn, pfct;
  bool const pert_eval_ok = eval_only( *oc, pvar, pinp, peqn, pfct );
  ok = check_bool( "source perturbed residual/output evaluation succeeds", pert_eval_ok ) && ok;
  DerivData const src_pder  = deriv_eqn_data( *oc, pvar, pinp );
  DerivData const src_pfder = deriv_fct_data( *oc, pvar, pinp );

  std::vector<double> tau_var = var, tau_eqn, tau_fct;
  bool const trace_perturbed = perturb_one_trace_variable( *oc, tau_var, 1.0e-3 );
  bool const trace_eval_ok = trace_perturbed && eval_only( *oc, tau_var, inp, tau_eqn, tau_fct );
  double const trace_sensitivity = trace_eval_ok ? max_diff( eqn, tau_eqn ) : 0.0;
  ok = check_bool( "source IC_TRACE tau perturbation changes residual",
                   imp == OCFESLV::Options::IC_TRACE
                 ? trace_eval_ok && trace_sensitivity > 1.0e-12
                 : !trace_perturbed ) && ok;

  OCFESLV cc( *oc );
  OCFESLV asg;
  asg = *oc;
  OCFESLV expl;
  bool explicit_ok = false;
  try{
    explicit_ok = expl.deep_copy_from( *oc );
  }
  catch(...){
    explicit_ok = false;
  }
  ok = check_bool( "explicit deep_copy_from succeeds", explicit_ok ) && ok;

  struct CopyCase { char const* name; OCFESLV* env; };
  std::vector<CopyCase> cases = {
    { "copy constructor", &cc },
    { "assignment",       &asg },
    { "explicit copy",    &expl }
  };

  for( auto const& c : cases ){
    std::vector<double> v2, i2, e2, f2;
    bool const init_eval_ok = init_eval( *c.env, v2, i2, e2, f2 );
    ok = check_bool( std::string( c.name ) + ": init/eval succeeds", init_eval_ok ) && ok;
    if( !init_eval_ok ) continue;

    RefCheck const cr = check_initialized_references( *c.env, v2, i2 );
    ok = check_bool( std::string( c.name ) + ": dimensions/container sizes",
                     c.env->n_colloc_sta()   == oc->n_colloc_sta()
                  && c.env->n_colloc_inp()   == oc->n_colloc_inp()
                  && c.env->n_colloc_eqn()   == oc->n_colloc_eqn()
                  && c.env->n_colloc_fct()   == oc->n_colloc_fct()
                  && c.env->n_colloc_trace() == oc->n_colloc_trace()
                  && c.env->states_colloc().size() == oc->states_colloc().size() ) && ok;
    ok = check_bool( std::string( c.name ) + ": options and derivative cache dimensions",
                     c.env->options.INTERFACE.IMPOSITION == oc->options.INTERFACE.IMPOSITION
                  && c.env->options.INTERFACE.TYPE  == oc->options.INTERFACE.TYPE
                  && c.env->options.REDUCE.ORDER    == oc->options.REDUCE.ORDER
                  && c.env->options.CLASSIFY.MODE        == oc->options.CLASSIFY.MODE
                  && c.env->n_colloc_nnz()          == oc->n_colloc_nnz()
                  && c.env->n_colloc_eqn_nnz()      == oc->n_colloc_eqn_nnz()
                  && c.env->n_colloc_fct_nnz()      == oc->n_colloc_fct_nnz()
                  && c.env->n_deriv_colors()        == oc->n_deriv_colors() ) && ok;
    ok = check_bool( std::string( c.name ) + ": output row blocks preserved",
                     c.env->blk_fct(0) == oc->blk_fct(0)
                  && c.env->blk_fct(1) == oc->blk_fct(1)
                  && c.env->blk_fct(2) == oc->blk_fct(2)
                  && c.env->blk_fct(3) == oc->blk_fct(3)
                  && c.env->row_fct(0) == oc->row_fct(0) ) && ok;
    ok = check_error( std::string( c.name ) + ": copied PDE classification",
                      max_classification_diff( *oc, *c.env ), 1e-10 ) && ok;
    ok = check_error( std::string( c.name ) + ": copied principal-symbol signature",
                      max_symbol_signature_diff( *oc, *c.env ), 0.0 ) && ok;
    ok = check_error( std::string( c.name ) + ": copied collocation nodes",
                      max_diff_nodes( *oc, *c.env ), 1e-14 ) && ok;
    ok = check_error( std::string( c.name ) + ": copied boundary-constructor t domain",
                      max_domain_diff( *oc, *c.env, "t" ), 1e-14 ) && ok;
    ok = check_error( std::string( c.name ) + ": copied non-uniform x domain",
                      max_domain_diff( *oc, *c.env, "x" ), 1e-14 ) && ok;
    ok = check_error( std::string( c.name ) + ": init primitive state values",
                      cr.err_state, 1e-10 ) && ok;
    ok = check_error( std::string( c.name ) + ": init distributed input values",
                      cr.err_input, 1e-12 ) && ok;
    ok = check_error( std::string( c.name ) + ": init auxiliary values",
                      cr.err_aux, 1e-9 ) && ok;
    ok = check_error( std::string( c.name ) + ": initialized state vector matches source",
                      max_diff( var, v2 ), 1e-10 ) && ok;
    ok = check_error( std::string( c.name ) + ": initialized input vector matches source",
                      max_diff( inp, i2 ), 1e-12 ) && ok;
    ok = check_error( std::string( c.name ) + ": equation residual matches source",
                      max_diff( eqn, e2 ), 1e-10 ) && ok;
    ok = check_error( std::string( c.name ) + ": output value matches source",
                      max_diff( fct, f2 ), 1e-10 ) && ok;
    ok = check_error( std::string( c.name ) + ": residual at initialized reference",
                      max_abs( e2 ), 1e-9 ) && ok;

    DerivData const der2 = deriv_eqn_data( *c.env, v2, i2 );
    DerivData const fder2 = deriv_fct_data( *c.env, v2, i2 );
    ok = check_bool( std::string( c.name ) + ": equation derivative sparsity matches source",
                     same_pattern( src_der, der2 ) ) && ok;
    ok = check_error( std::string( c.name ) + ": equation derivative values match source",
                      max_diff_deriv_values( src_der, der2 ), 1e-10 ) && ok;
    ok = check_bool( std::string( c.name ) + ": output derivative sparsity matches source",
                     same_pattern( src_fder, fder2 ) ) && ok;
    ok = check_error( std::string( c.name ) + ": output derivative values match source",
                      max_diff_deriv_values( src_fder, fder2 ), 1e-10 ) && ok;

    std::vector<double> ceqn, cfct;
    bool const copy_pert_eval_ok = eval_only( *c.env, pvar, pinp, ceqn, cfct );
    ok = check_bool( std::string( c.name ) + ": perturbed eval succeeds",
                     copy_pert_eval_ok ) && ok;
    if( copy_pert_eval_ok ){
      ok = check_error( std::string( c.name ) + ": perturbed residual matches source",
                        max_diff( peqn, ceqn ), 1e-10 ) && ok;
      ok = check_error( std::string( c.name ) + ": perturbed outputs match source",
                        max_diff( pfct, cfct ), 1e-10 ) && ok;
    }
    DerivData const pder2 = deriv_eqn_data( *c.env, pvar, pinp );
    DerivData const pfder2 = deriv_fct_data( *c.env, pvar, pinp );
    ok = check_bool( std::string( c.name ) + ": perturbed equation derivative sparsity matches source",
                     same_pattern( src_pder, pder2 ) ) && ok;
    ok = check_error( std::string( c.name ) + ": perturbed equation derivative values match source",
                      max_diff_deriv_values( src_pder, pder2 ), 1e-10 ) && ok;
    ok = check_bool( std::string( c.name ) + ": perturbed output derivative sparsity matches source",
                     same_pattern( src_pfder, pfder2 ) ) && ok;
    ok = check_error( std::string( c.name ) + ": perturbed output derivative values match source",
                      max_diff_deriv_values( src_pfder, pfder2 ), 1e-10 ) && ok;

    if( imp == OCFESLV::Options::IC_TRACE ){
      std::vector<double> cteqn, ctfct;
      bool const copy_trace_eval_ok = eval_only( *c.env, tau_var, inp, cteqn, ctfct );
      ok = check_bool( std::string( c.name ) + ": tau-perturbed eval succeeds",
                       copy_trace_eval_ok ) && ok;
      if( copy_trace_eval_ok ){
        ok = check_error( std::string( c.name ) + ": tau-perturbed residual matches source",
                          max_diff( tau_eqn, cteqn ), 1e-10 ) && ok;
        ok = check_error( std::string( c.name ) + ": tau-perturbed outputs match source",
                          max_diff( tau_fct, ctfct ), 1e-10 ) && ok;
      }
      ok = check_bool( std::string( c.name ) + ": equation Jacobian still contains tau columns",
                       has_columns_at_or_above( der2, first_trace_col ) ) && ok;
    }
  }

  return ok;
}

// ---------------------------------------------------------------------------
// Index-2 reduction deep-copy suite (corpus M2).
//
// The heat-equation suite above is index-1, so its reduction_plan() is empty
// and reduced_dof_audit().ran is false -- it verifies the *empty* records copy
// faithfully, but does not exercise a populated reduction plan.  This suite
// builds the minimal index-2 DAE
//
//     x' - y   = 0          (ODE,   INTERIOR at ALL-LB)
//     x - a(t) = 0          (ALG, no y -> index 2, INTERIOR at ALL)
//     x - a(0) = 0          (consistency IC, INITIAL at LB)
//     a(t) = 1 + 0.5 t - 0.3 t^2       (exact: x = a, y = a')
//
// so the high-index reducer differentiates the ALG constraint once, pins y, and
// records one t_ReductionPlan::t_Assign plus a nonzero-work t_DofAudit.  It then
// checks that deep copies reproduce both records exactly -- including that the
// copied plan's captured constraint/witness FFVars are rebound to the *copy's*
// working DAG (the point of registering the reduction roots in _register).
// ---------------------------------------------------------------------------

static std::unique_ptr<OCFESLV> make_index2_source()
{
  std::unique_ptr<OCFESLV> oc;

  {
    FFGraph DAG;
    FFVar t = DAG.add_var( "t" );
    FFVar x = DAG.add_var( "x(t)" );
    FFVar y = DAG.add_var( "y(t)" );          // ALGEBRAIC, HIDDEN (absent from ALG)
    FFPartial OpP;

    FFVar Ae   = 1.0 + 0.5*t - 0.3*t*t;       // a(t)
    FFVar DIFF = OpP( x, t ) - y;             // x' - y = 0
    FFVar ALG  = x - Ae;                      // 0 = x - a  (no y -> index 2)
    FFVar ICx  = x - Ae;                      // consistency at t=0

    oc.reset( new OCFESLV( &DAG ) );
    oc->add_domain( t, FFDom( 0., 0.5, 3, FFDom::LGR, 6 ) );
    oc->add_state( x, {t} );
    oc->add_state( y, {t} );
    oc->update_ref( x, []( OCFESLV::t_Coord const& ){ return 1.0; } );  // adversarial (non-exact) seed
    oc->update_ref( y, []( OCFESLV::t_Coord const& ){ return 0.0; } );

    OCFESLV::EqnOptions io( OCFESLV::EqnRole::INTERIOR, 0 ), ii( OCFESLV::EqnRole::INITIAL, 0 );
    int const T_INT = FFDom::ALL - FFDom::LB;
    oc->add_equation( DIFF, {t}, {T_INT},      io );
    oc->add_equation( ALG,  {t}, {FFDom::ALL}, io );
    // 2026-09-23: index 2 with one differential state (x), so no free initial data -- ICx was a consistency
    // condition, now materialised at t=0 by REDUCE.HIDDEN_IC.
    (void)ICx;
    oc->set_evolution_domain( t );

    oc->options.REDUCE.ORDER     = OCFESLV::Options::RED_FULL;
    oc->options.REDUCE.HIDDEN_IC = true;
    oc->options.CLASSIFY.MODE         = OCFESLV::Options::CLASS_AUTO;
    oc->options.INTERFACE.TYPE   = OCFESLV::Options::IC_AUTO;
    oc->options.INTERFACE.IMPOSITION  = OCFESLV::Options::IC_WEAK;
    oc->options.SOLVE.MARCHING   = false;   // rev105: structural test -- see the note above
    oc->options.DISPLAY_LEVEL    = 0;

    if( !oc->setup() )
      throw std::runtime_error( "index-2 (corpus M2) setup() failed" );
  }

  // User FFGraph destroyed; the source now relies on its private local DAG.
  return oc;
}

static bool reduction_records_match
( std::string const& label, OCFESLV const& src, OCFESLV const& cp )
{
  bool ok = true;

  auto const& sa = src.reduced_dof_audit();
  auto const& ca = cp.reduced_dof_audit();
  bool const audit_ok =
       sa.ran==ca.ran && sa.trace==ca.trace && sa.rows==ca.rows && sa.cols==ca.cols
    && sa.rank==ca.rank && sa.deficiency==ca.deficiency && sa.square==ca.square
    && sa.ok()==ca.ok();
  ok = check_bool( label + ": reduced_dof_audit matches source", audit_ok ) && ok;

  auto const& sp = src.reduction_plan();
  auto const& cprp = cp.reduction_plan();
  bool const shape_ok =
       sp.assigns.size()==cprp.assigns.size()
    && sp.max_index==cprp.max_index
    && sp.resolved==cprp.resolved
    && sp.empty()==cprp.empty();
  ok = check_bool( label + ": reduction_plan shape matches source", shape_ok ) && ok;

  bool fields_ok = shape_ok;
  bool rebind_ok = shape_ok;
  for( size_t i=0; shape_ok && i<sp.assigns.size(); ++i ){
    auto const& s = sp.assigns[i];
    auto const& c = cprp.assigns[i];
    if( s.block_id != c.block_id || s.n_diff != c.n_diff ) fields_ok = false;

    // Every captured FFVar that is DAG-bound in the source must, in the copy, be
    // DAG-bound to the COPY's working DAG -- never null (dropped) and never still
    // pointing into the source DAG (dangling cross-DAG reference).  FFGraph::insert
    // preserves indices for leaf VARIABLES only; an operation/auxiliary node -- and
    // the differentiated constraint residual is one -- is re-created in the
    // destination DAG and legitimately receives a fresh index, so assert
    // id-preservation only for genuine variables (e.g. the pinned witness).
    auto rebound = [&]( FFVar const& s_v, FFVar const& c_v ){
      if( !s_v.dag() ) return true;                 // numeric constant: no rebind needed
      bool bound = c_v.dag() && c_v.dag() == cp.dag() && c_v.dag() != src.dag();
      if( s_v.id().first == FFVar::VAR )
        bound = bound && ( c_v.id() == s_v.id() );
      return bound;
    };
    if( !rebound( s.constraint, c.constraint ) ) rebind_ok = false;
    if( !rebound( s.pinned_var, c.pinned_var ) ) rebind_ok = false;
  }
  ok = check_bool( label + ": reduction_plan block_id/n_diff match source", fields_ok ) && ok;
  ok = check_bool( label + ": reduction_plan FFVars rebound to copy DAG (var id preserved)", rebind_ok ) && ok;

  return ok;
}

static bool run_index2_copy_suite()
{
  std::cout << "\n------------------------------------------------------------\n";
  std::cout << "  Deep-copy case: index-2 reduction (corpus M2)\n";
  std::cout << "------------------------------------------------------------\n";

  bool ok = true;
  auto oc = make_index2_source();

  std::cout << "source dimensions: nEqn=" << oc->n_colloc_eqn()
            << " nVar=" << oc->n_colloc_sta()
            << " reduction assigns=" << oc->reduction_plan().assigns.size()
            << " max_index=" << oc->reduction_plan().max_index
            << " audit(ran/def/sq)=" << oc->reduced_dof_audit().ran << "/"
            << oc->reduced_dof_audit().deficiency << "/"
            << oc->reduced_dof_audit().square << "\n";

  // Guards: the source must genuinely exercise the reducer, else the copy checks
  // below would be vacuous.
  ok = check_bool( "source triggered index reduction (non-empty plan)",
                   !oc->reduction_plan().empty() ) && ok;
  ok = check_bool( "source reduction resolved",
                   oc->reduction_plan().resolved ) && ok;
  ok = check_bool( "source DOF audit ran, square and full structural rank",
                   oc->reduced_dof_audit().ran && oc->reduced_dof_audit().ok() ) && ok;
  ok = check_bool( "source plan FFVars are bound to the source local DAG",
                   oc->reduction_plan().assigns.empty()
                 ? true
                 : ( oc->reduction_plan().assigns.front().constraint.dag() == oc->dag()
                  && oc->reduction_plan().assigns.front().pinned_var.dag() == oc->dag() ) ) && ok;

  std::vector<double> var, inp, eqn, fct;
  ok = check_bool( "source init/eval succeeds", init_eval( *oc, var, inp, eqn, fct ) ) && ok;
  DerivData const src_der = deriv_eqn_data( *oc, var, inp );

  OCFESLV cc( *oc );
  OCFESLV asg;
  asg = *oc;
  OCFESLV expl;
  bool explicit_ok = false;
  try{ explicit_ok = expl.deep_copy_from( *oc ); }
  catch(...){ explicit_ok = false; }
  ok = check_bool( "explicit deep_copy_from succeeds", explicit_ok ) && ok;

  struct CopyCase { char const* name; OCFESLV* env; };
  std::vector<CopyCase> cases = {
    { "index-2 copy constructor", &cc },
    { "index-2 assignment",       &asg },
    { "index-2 explicit copy",    &expl }
  };

  for( auto const& c : cases ){
    std::vector<double> v2, i2, e2, f2;
    bool const ie = init_eval( *c.env, v2, i2, e2, f2 );
    ok = check_bool( std::string( c.name ) + ": init/eval succeeds", ie ) && ok;
    if( !ie ) continue;

    ok = check_bool( std::string( c.name ) + ": reduced dimensions match source",
                     c.env->n_colloc_sta()   == oc->n_colloc_sta()
                  && c.env->n_colloc_eqn()   == oc->n_colloc_eqn()
                  && c.env->n_colloc_trace() == oc->n_colloc_trace() ) && ok;
    ok = check_error( std::string( c.name ) + ": reduced residual matches source",
                      max_diff( eqn, e2 ), 1e-10 ) && ok;

    DerivData const d2 = deriv_eqn_data( *c.env, v2, i2 );
    ok = check_bool( std::string( c.name ) + ": reduced equation derivative sparsity matches source",
                     same_pattern( src_der, d2 ) ) && ok;
    ok = check_error( std::string( c.name ) + ": reduced equation derivative values match source",
                      max_diff_deriv_values( src_der, d2 ), 1e-10 ) && ok;

    ok = reduction_records_match( c.name, *oc, *c.env ) && ok;
  }

  return ok;
}

} // namespace

int main()
{
  std::cout << "\n============================================================\n";
  std::cout << "  OCFE_deepcopy: revised OCFESLV setup/deep-copy checks\n";
  std::cout << "============================================================\n";

  bool ok = true;

  // rev194: a not-set-up source is REPORTED (false), not thrown -- the copy constructor and assignment cannot
  // report a throw, so an unguarded copy of such an environment used to terminate the process.
  OCFESLV unsetup;
  OCFESLV reject;
  bool rejected_unsetup = false;
  bool no_throw = true;
  try{
    rejected_unsetup = !reject.deep_copy_from( unsetup );
  }
  catch(...){
    no_throw = false;
  }
  ok = check_bool( "deep_copy_from rejects unsetup source (returns false)", rejected_unsetup ) && ok;
  ok = check_bool( "deep_copy_from does not throw on unsetup source",       no_throw ) && ok;

  // and the copy constructor of a not-set-up environment is usable: an empty, not-set-up environment.
  bool copy_ctor_ok = true;
  try{
    OCFESLV copy_of_unsetup( unsetup );
    copy_ctor_ok = !copy_of_unsetup.is_setup();
  }
  catch(...){ copy_ctor_ok = false; }
  ok = check_bool( "copy ctor of unsetup source does not terminate", copy_ctor_ok ) && ok;

  ok = run_copy_suite( OCFESLV::Options::IC_WEAK   ) && ok;
  ok = run_copy_suite( OCFESLV::Options::IC_STRONG ) && ok;
  ok = run_copy_suite( OCFESLV::Options::IC_TRACE  ) && ok;
  ok = run_index2_copy_suite() && ok;

  std::cout << "\nDeep-copy/reference initialization checks: "
            << ( ok? "PASS": "FAIL" ) << "\n";
  return ok? 0: 1;
}
