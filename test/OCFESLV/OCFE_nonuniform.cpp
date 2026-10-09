// Copyright (C) Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the EPL-2.0 with GPL-2.0-or-later as a Secondary License; see the LICENSE file.
// SPDX-License-Identifier: EPL-2.0 OR GPL-2.0-or-later
//
// OCFE_nonuniform.cpp
// ----------------------------
// Regression test for nonuniform finite-element partitions in FFDom.
// It checks that physical node placement, partial-derivative scaling and
// integral quadrature use element-specific widths rather than a uniform width.

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#ifndef OCFE_OCFESLV_HEADER
#define OCFE_OCFESLV_HEADER "ocfeslv.hpp"
#endif
#include OCFE_OCFESLV_HEADER
#include "test_deriv_utils.hpp"

using namespace mc;

static double u_exact( double t )  { return t*t; }
static double q_exact()            { return 1.0/3.0; }

static double max_abs( std::vector<double> const& x )
{
  double m = 0.;
  for( double v : x ) m = std::max( m, std::fabs(v) );
  return m;
}

static bool run_case( std::string const& label, FFDom const& dom )
{
  std::cout << "\n========== " << label << " ==========" << std::endl;

  FFGraph DAG;
  FFVar t = DAG.add_var( "t" );
  FFVar u = DAG.add_var( "u(t)" );

  FFPartial  OpP;
  FFIntegral OpI;

  FFVar ODE = OpP( u, t ) - 2.0*t;
  FFVar IC  = u;
  FFVar Q   = OpI( u, t );
  FFVar PATH = OpP( u, t ) - 2.0*t;

  OCFESLV oc( &DAG );
  oc.add_domain( t, dom );
  oc.add_state( u, {t} );
  oc.update_ref( u, [&]( OCFESLV::t_Coord const& coord ){
    return u_exact( coord.at( t ) );
  } );
  oc.set_evolution_domain( t );
  oc.options.SOLVE.MARCHING = false;   // direct node_colloc/eval reference test over a nonuniform
                                       // partition: keep every element -- do not auto-march/collapse.
  oc.add_equation( ODE, {t}, {FFDom::ALL-FFDom::LB},
                   OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( IC,  {t}, {FFDom::LB},
                   OCFESLV::EqnOptions( OCFESLV::EqnRole::INITIAL, 0 ) );
  oc.add_output( Q );
  oc.add_output( PATH, {t}, {FFDom::ALL} );

  if( !oc.setup() ){
    std::cerr << "ERROR: setup failed\n";
    return false;
  }

  size_t const nVar = oc.n_colloc_sta();
  size_t const nEqn = oc.n_colloc_eqn();
  size_t const nOut = oc.n_colloc_fct();

  auto nodes = oc.node_colloc( u );
  std::vector<double> var;
  var.reserve( nodes.size() );
  for( auto const& node : nodes ) var.push_back( u_exact( node[0] ) );

  std::vector<double> varInit, inpInit;
  if( !oc.init( varInit, inpInit, nullptr ) ){
    std::cerr << "ERROR: init failed\n";
    return false;
  }

  double init_err = 0.;
  for( size_t i=0; i<var.size(); ++i )
    init_err = std::max( init_err, std::fabs( varInit[i] - var[i] ) );
  var = varInit;

  std::vector<double> res( nEqn, 0.0 ), out( nOut, 0.0 );
  if( !oc.eval( res.data(), out.data(), var.data(), nullptr, nullptr ) ){
    std::cerr << "ERROR: eval failed\n";
    return false;
  }

  auto blk_Q    = oc.blk_fct( 0 );
  auto blk_PATH = oc.blk_fct( 1 );

  double const res_max = max_abs( res );
  double const q_err   = std::fabs( out[blk_Q.first] - q_exact() );
  double path_max = 0.;
  for( size_t i=0; i<blk_PATH.second; ++i )
    path_max = std::max( path_max, std::fabs( out[blk_PATH.first+i] ) );

  // Check representative physical nodes.  With LGR the global lower bound is
  // present; each element's first local node is its lower boundary.
  bool node_ok = !nodes.empty();
  if( node_ok ){
    size_t pos = 0;
    for( size_t ie=0; ie<dom.n_elem; ++ie ){
      node_ok &= std::fabs( nodes[pos][0] - dom.elem_lo( ie ) ) <= 1e-13;
      pos += dom.n_node;
    }
  }

  bool const value_ok = init_err <= 1e-12
                     && res_max  <= 5e-11
                     && q_err    <= 5e-13
                     && path_max <= 5e-11
                     && node_ok;

  std::cout << "nVar=" << nVar << " nEqn=" << nEqn << " nOut=" << nOut << "\n";
  std::cout << "node placement: " << ( node_ok ? "PASS" : "FAIL" ) << "\n";
  std::cout << "max|init-exact|=" << std::scientific << std::setprecision(4) << init_err
            << " max|res|=" << res_max
            << " |Q-1/3|=" << q_err
            << " max|path|=" << path_max
            << "  " << ( value_ok ? "PASS" : "FAIL" ) << "\n";

  bool const deriv_ok = mc_test::check_oc_derivatives
    ( oc, label, var, nullptr, nullptr );

  return value_ok && deriv_ok;
}

int main()
{
  bool ok = true;
  ok &= run_case( "nonuniform boundaries", FFDom( {0.0,0.1,0.4,1.0}, FFDom::LGR, 4 ) );
  ok &= run_case( "nonuniform lengths", FFDom( 0.0, {0.1,0.3,0.6}, FFDom::LGR, 4 ) );

  std::cout << "\nOCFE_nonuniform: " << ( ok ? "PASS" : "FAIL" ) << "\n";
  return ok ? 0 : 1;
}
