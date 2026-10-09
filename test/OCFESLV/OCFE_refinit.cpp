// Copyright (C) Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the EPL-2.0 with GPL-2.0-or-later as a Secondary License; see the LICENSE file.
// SPDX-License-Identifier: EPL-2.0 OR GPL-2.0-or-later
//
// OCFE_refinit.cpp
// ----------------
// Regression test for OCFESLV domain-dependent reference functions and
// OCFESLV::init() pointer-style initialization.  The test checks that
//
//   * reference functions use std::map<FFVar,double,lt_FFVar>, keyed by the
//     original domain variables rather than by an assumed coordinate order;
//   * OCFESLV::init(double*,double*,double const*) initializes collocated states
//     and distributed inputs from those reference functions; and
//   * auxiliary states introduced by setup-time reduce_order() receive
//     collocation-consistent reference values from their defining expression.

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <map>
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

static double state_ref_value( double t, double x )
{
  return 1.0 + 0.5*t + 3.0*x + x*x;
}

static double input_ref_value( double t, double x )
{
  return 2.0 - 0.25*t + 0.1*x;
}

static double aux_ref_value( double /*t*/, double x )
{
  // d/dx state_ref_value(t,x), which is represented exactly by the local
  // collocation derivative for the quadratic state reference used here.
  return 3.0 + 2.0*x;
}

static std::pair<double,double> infer_tx
( std::vector<double> const& coord )
{
  if( coord.size() != 2 )
    throw std::runtime_error( "Expected two coordinates" );

  // The test chooses t in [2,3] and x in [0,1], so the domain variable can be
  // identified without assuming node_colloc()'s coordinate ordering.
  if( coord[0] > 1.5 && coord[1] <= 1.5 ) return {coord[0], coord[1]};
  if( coord[1] > 1.5 && coord[0] <= 1.5 ) return {coord[1], coord[0]};
  throw std::runtime_error( "Could not infer (t,x) from node coordinates" );
}

static void accumulate_error
( double& maxerr, double got, double expect )
{
  maxerr = std::max( maxerr, std::fabs( got - expect ) );
}

static bool pass_fail
( std::string const& label, double err, double tol )
{
  bool const pass = ( err <= tol );
  std::cout << std::left << std::setw(54) << label
            << " maxerr=" << std::scientific << std::setprecision(3) << err
            << " tol=" << tol << "  " << ( pass? "PASS": "FAIL" ) << "\n";
  return pass;
}

} // namespace

int main()
{
  std::cout << "\n========== OCFE_refinit: FUNCTION REFERENCES AND init() ==========" << "\n";

  FFGraph DAG;
  FFVar t = DAG.add_var( "t" );
  FFVar x = DAG.add_var( "x" );
  FFVar T = DAG.add_var( "T(t,x)" );
  FFVar a = DAG.add_var( "a(t,x)" );

  OCFESLV::t_Fun Tref = [t,x]( OCFESLV::t_Coord const& coord ) -> double {
    return state_ref_value( coord.at(t), coord.at(x) );
  };

  OCFESLV::t_Fun Aref = [t,x]( OCFESLV::t_Coord const& coord ) -> double {
    return input_ref_value( coord.at(t), coord.at(x) );
  };

  FFPartial OpPartial;
  FFVar PDE = OpPartial( OpPartial( T, x ), x );

  OCFESLV oc( &DAG );
  oc.add_domain( t, FFDom( 2.0, 3.0, 1, FFDom::LGL, 4 ), 2.25 );
  oc.add_domain( x, FFDom( 0.0, 1.0, 2, FFDom::LGL, 5 ), 0.40 );
  oc.add_state( T, {t,x}, Tref );
  oc.add_input( a, {t,x}, Aref );
  oc.add_equation( PDE, {t,x}, {FFDom::ALL, FFDom::ALL-FFDom::LB-FFDom::UB},
                   OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );

  oc.options.REDUCE.ORDER = OCFESLV::Options::RED_FULL;
  oc.options.CLASSIFY.MODE     = OCFESLV::Options::CLASS_NONE;

  double const tref_scalar = 2.25;
  double const xref_scalar = 0.40;
  double err_ref_before = 0.0;
  accumulate_error( err_ref_before, oc.ref( T ), state_ref_value( tref_scalar, xref_scalar ) );
  accumulate_error( err_ref_before, oc.ref( a ), input_ref_value( tref_scalar, xref_scalar ) );

  if( !oc.setup() ){
    std::cerr << "ERROR: OCFESLV::setup() failed.\n";
    return 1;
  }

  std::cout << oc;

  std::vector<double> var( oc.n_colloc_sta(), 0.0 );
  std::vector<double> inp( oc.n_colloc_inp(), 0.0 );
  if( !oc.init( var.data(), inp.data(), nullptr ) ){
    std::cerr << "ERROR: OCFESLV::init(double*,double*,double const*) failed.\n";
    return 1;
  }

  double err_ref_after = 0.0;
  accumulate_error( err_ref_after, oc.ref( T ), state_ref_value( tref_scalar, xref_scalar ) );
  accumulate_error( err_ref_after, oc.ref( a ), input_ref_value( tref_scalar, xref_scalar ) );

  double err_state = 0.0;
  double err_aux   = 0.0;
  bool found_state = false;
  bool found_aux   = false;

  size_t offset = 0;
  for( auto const& st : oc.states_colloc() ){
    auto const nodes = oc.node_colloc( st );
    bool const is_user_state = ( std::string( st.name() ) == std::string( T.name() ) );
    for( size_t i=0; i<nodes.size(); ++i ){
      auto const [tt,xx] = infer_tx( nodes[i] );
      if( is_user_state ){
        found_state = true;
        accumulate_error( err_state, var[offset+i], state_ref_value( tt, xx ) );
      }
      else{
        found_aux = true;
        accumulate_error( err_aux, var[offset+i], aux_ref_value( tt, xx ) );
      }
    }
    offset += nodes.size();
  }

  // Exercise deep_copy_from()/copy-constructor handling of function references
  // and auxiliary reference definitions.
  OCFESLV oc_copy( oc );
  std::vector<double> var_copy( oc_copy.n_colloc_sta(), 0.0 );
  std::vector<double> inp_copy( oc_copy.n_colloc_inp(), 0.0 );
  if( !oc_copy.init( var_copy.data(), inp_copy.data(), nullptr ) ){
    std::cerr << "ERROR: copied OCFESLV::init(double*,double*,double const*) failed.\n";
    return 1;
  }

  double err_copy = 0.0;
  if( var_copy.size() != var.size() || inp_copy.size() != inp.size() ){
    std::cerr << "ERROR: copied OCFESLV dimensions do not match source.\n";
    return 1;
  }
  for( size_t i=0; i<var.size(); ++i ) accumulate_error( err_copy, var_copy[i], var[i] );
  for( size_t i=0; i<inp.size(); ++i ) accumulate_error( err_copy, inp_copy[i], inp[i] );

  double err_input = 0.0;
  auto const input_nodes = oc.node_colloc( a );
  if( input_nodes.size() != inp.size() ){
    std::cerr << "ERROR: input node count does not match initialized input vector size.\n";
    return 1;
  }
  for( size_t i=0; i<input_nodes.size(); ++i ){
    auto const [tt,xx] = infer_tx( input_nodes[i] );
    accumulate_error( err_input, inp[i], input_ref_value( tt, xx ) );
  }

  bool ok = true;
  ok = pass_fail( "function ref() before setup", err_ref_before, 1e-12 ) && ok;
  ok = pass_fail( "function ref() after setup/local DAG copy", err_ref_after, 1e-12 ) && ok;
  ok = pass_fail( "init() primitive state values", err_state, 1e-10 ) && ok;
  ok = pass_fail( "init() distributed input values", err_input, 1e-10 ) && ok;
  ok = pass_fail( "init() reduce_order auxiliary values", err_aux, 1e-9 ) && ok;
  ok = pass_fail( "copy constructor preserves initialized references", err_copy, 1e-10 ) && ok;

  if( !found_state ){
    std::cerr << "ERROR: original state not found after setup.\n";
    ok = false;
  }
  if( !found_aux ){
    std::cerr << "ERROR: reduce_order did not introduce an auxiliary state.\n";
    ok = false;
  }

  std::cout << "\nReference initialization checks: " << ( ok? "PASS": "FAIL" ) << "\n";
  return ok? 0: 1;
}
