// Copyright (C) Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the EPL-2.0 with GPL-2.0-or-later as a Secondary License; see the LICENSE file.
// SPDX-License-Identifier: EPL-2.0 OR GPL-2.0-or-later
//
// ============================================================================
//  OCFE_capture_expr.cpp -- a deferred value whose SOURCE IS AN EXPRESSION
//
//  Every capture in the corpus has a source that is already a state, so the path
//  that materialises an expression source as a state (with a defining LINK row)
//  is never exercised.  This driver exercises it: dx/dt = 1, x(0) = 0, with two
//  outputs -- INT x^2 dt (an expression) and INT x dt (a bare state) -- whose exact
//  values are 1/3 and 1/2.
//
//  It is the gate for moving that materialisation out of the model and into
//  OCFESLV: what OCFESLV computes, and the system it solves, must not change.
// ============================================================================
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include "ffunc.hpp"
#include OCFE_OCFESLV_HEADER
using namespace mc;

static int g_pass = 0, g_fail = 0;
static void check( std::string const& nm, bool ok )
{
  std::cout << "  " << std::left << std::setw(56) << nm << std::right << " " << ( ok? "PASS": "FAIL" ) << std::endl;
  ( ok? g_pass: g_fail )++;
}

int main()
{
  std::cout << "OCFE_capture_expr ** a capture whose source is an expression" << std::endl;
  FFGraph DAG; OCFESLV oc( &DAG );
  FFVar t = DAG.add_var("t"), x = DAG.add_var("x(t)");
  FFPartial OpP; FFIntegral OpI;
  oc.add_domain( t, FFDom( 0., 1., 4, FFDom::LGR, 4 ) );
  oc.add_state ( x, { t } );
  oc.set_evolution_domain( t );
  oc.add_equation( OpP(x,t) - 1., { t }, { FFDom::ALL - FFDom::LB },
                   OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( x - 0.,        { t }, { FFDom::LB },
                   OCFESLV::EqnOptions( OCFESLV::EqnRole::INITIAL, 0 ) );
  oc.add_output( OpI( x * x, t ) );          // an EXPRESSION source: materialised today
  oc.add_output( OpI( x, t ) );              // a bare state source: never materialised
  oc.update_ref( x, 0.5 );
  oc.options.DISPLAY_LEVEL = 1;

  check( "setup() succeeds", oc.setup() );
  std::cout << "  states=" << oc.var_state().size() << " inputs=" << oc.var_input().size()
            << " eqns=" << oc.var_equation().size() << " class="
            << FFModel::pde_type_name( oc.pde_type().type ) << std::endl;
  size_t nlink = 0;
  for( auto const& e : oc.var_equation() ) if( e.opt && e.opt->role == FFModel::EqnRole::LINK ) ++nlink;
  for( auto const& C : oc.var_deferred() )
    std::cout << "  deferred: input=" << C.input.name() << ( C.accumulate? "  INTEGRAL": "  VALUE" )
              << "  of " << C.source.name()
              << ( oc.var_state().count( C.source )? "  (a state)": "  (an expression)" ) << std::endl;
  check( "two deferred values", oc.var_deferred().size() == 2 );
  check( "OCFESLV solves a system with the source as a state", oc.var_state().size() == 2 && nlink == 1 );

  std::vector<double> xv( 4096, 0. );
  OCFESLV::SolveReport rep = oc.solve( xv.data(), nullptr, nullptr );
  check( "solve() converges", rep.converged );
  std::vector<double> const& F = oc.val_functions();
  std::cout << std::scientific << std::setprecision(12);
  for( size_t i = 0; i < F.size(); ++i ) std::cout << "  output[" << i << "] = " << F[i] << std::endl;
  check( "INT x^2 dt = 1/3 (the expression capture)", F.size() > 0 && std::fabs( F[0] - 1./3. ) < 1e-8 );
  check( "INT x dt = 1/2 (the state capture)",        F.size() > 1 && std::fabs( F[1] - 0.5   ) < 1e-8 );

  std::cout << "\n  OCFE_capture_expr: " << g_pass << " passed, " << g_fail << " failed" << std::endl;
  return g_fail ? 1 : 0;
}
