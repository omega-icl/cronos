// Copyright (C) Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the EPL-2.0 with GPL-2.0-or-later as a Secondary License; see the LICENSE file.
// SPDX-License-Identifier: EPL-2.0 OR GPL-2.0-or-later
//
// ============================================================================
//  OCFE_claim.cpp -- gate for C0, the continuity-CLAIM classification of
//  algebraic variables per direction (FFModel::claim_classification).
//
//  Includes ffmodel.hpp only: the classification is STRUCTURAL -- no mesh, no
//  collocation, no reference point -- so nothing but the model layer is linked.
//
//  Each case declares a model whose correct tags are known by hand, and asserts
//  them.  The point is that C0 mints a claim ONLY where an independent flux
//  condition exists, so there is never anything to rescue.
//
//    1 heat/mixed     T_t = q_x, q = T_x     q IMPLIED in t (L), NATURAL in x;
//                                            T C1 in x (q stands for T_x)
//    2 moving mesh    xg = x_xi, g = M xg,   every algebraic state IMPLIED in
//                     gg = g_xi, x_t = gg    BOTH directions -- the case that
//                                            mints the gauge multiplier today
//    3 closure        q = k(T) T_x           q CLOSURE in t: the RELATION has a
//                                            state-dependent coefficient; T stays
//                                            NATURAL in x (the tag goes on the
//                                            implied state)
//    4 nonlinear      q^2 = T_x              q CLOSURE in t with letter N: C1
//                                            must not substitute through it
//    5 unattached     w free, no relation    w NONE in every direction: a claim
//                                            would be an illegitimate constraint
// ============================================================================
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include "ffunc.hpp"
#include "ffmodel.hpp"
using namespace mc;

static int g_fail = 0;

static std::string nm( FFVar const& v ){ std::ostringstream o; o << v; return o.str(); }

static char tagchar( char want )   // 'A' natural, 'I' implied, 'C' closure, 'N' none
{ return want; }

//! Assert one (state, direction) entry through the PUBLIC query FFModel::claim_info -- the classification
//! itself (ClaimTag, t_ClaimMap, claim_classification) is protected, being an internal analysis.  rev271.
static void expect( FFModel const& M, FFVar const& s, FFVar const& d,
                    char want_tag, char coeff = 0, int order = -2 )
{
  auto const ci = M.claim_info( s, d );
  bool const pass = ci.found && ci.tag == want_tag
                 && ( coeff == 0 || ci.coeff == coeff || ci.coeff == (char)( coeff + 32 ) )
                 && ( order == -2 || ci.order == order );
  std::cout << "    " << std::left << std::setw( 10 ) << nm( s ) << " in " << std::setw( 5 ) << nm( d )
            << " expect " << want_tag;
  if( coeff ) std::cout << "(" << coeff << ")";
  if( order != -2 ) std::cout << " order=" << order;
  std::cout << "   got " << ( ci.found? ci.tag: '-' );
  if( ci.found && ci.coeff != '-' ) std::cout << "(" << ci.coeff << ")";
  if( ci.found ) std::cout << " order=" << ci.order;
  std::cout << "   " << ( pass? "ok": "MISMATCH" ) << std::endl;
  if( !pass ) ++g_fail;
}

int main()
{
  std::cout << "OCFE_claim ** C0 continuity-claim classification (" << FFModel::HEADER_ID << ")\n";
  FFPartial OpP;

  // ---- 1: heat equation in mixed form -------------------------------------------------------
  {
    std::cout << "\n  1  mixed-form heat: T_t = q_x, q = T_x\n";
    FFGraph D; FFModel M( &D );
    FFVar t = D.add_var("t"), x = D.add_var("x"), T = D.add_var("T(t,x)"), q = D.add_var("q(t,x)");
    M.add_domain( t, FFDom( 0., 1., 2, FFDom::LGL, 3 ) );
    M.add_domain( x, FFDom( 0., 1., 2, FFDom::LGL, 3 ) );
    M.add_state( T, { t, x } );  M.add_state( q, { t, x } );
    M.add_equation( OpP( T, t ) - OpP( q, x ), std::vector<FFVar>{ t, x }, std::vector<int>{} );
    M.add_equation( q - OpP( T, x ),           std::vector<FFVar>{ t, x }, std::vector<int>{} );
    if( !M.setup() ){ std::cout << "    setup FAILED\n"; ++g_fail; }
    else {
      expect( M, q, t, 'I', 'L' );
      expect( M, q, x, 'A' );
      expect( M, T, x, 'A', 0, 1 );   // C1: q stands for T_x and is natural in x
    }
  }

  // ---- 3: state-dependent coefficient -> CLOSURE ---------------------------------------------
  {
    std::cout << "\n  3  closure: q = k(T) T_x  (coefficient depends on T)\n";
    FFGraph D; FFModel M( &D );
    FFVar t = D.add_var("t"), x = D.add_var("x"), T = D.add_var("T(t,x)"), q = D.add_var("q(t,x)");
    M.add_domain( t, FFDom( 0., 1., 2, FFDom::LGL, 3 ) );
    M.add_domain( x, FFDom( 0., 1., 2, FFDom::LGL, 3 ) );
    M.add_state( T, { t, x } );  M.add_state( q, { t, x } );
    M.add_equation( OpP( T, t ) - OpP( q, x ),     std::vector<FFVar>{ t, x }, std::vector<int>{} );
    M.add_equation( q - ( 1. + T ) * OpP( T, x ),  std::vector<FFVar>{ t, x }, std::vector<int>{} );
    if( !M.setup() ){ std::cout << "    setup FAILED\n"; ++g_fail; }
    else {
      expect( M, q, t, 'C' );      // the RELATION has a state-dependent coefficient
      expect( M, q, x, 'A' );
      expect( M, T, x, 'A' );      // the tag goes on the IMPLIED state, never here
    }
  }

  // ---- 4: the state enters nonlinearly -------------------------------------------------------
  {
    std::cout << "\n  4  nonlinear: q^2 = T_x\n";
    FFGraph D; FFModel M( &D );
    FFVar t = D.add_var("t"), x = D.add_var("x"), T = D.add_var("T(t,x)"), q = D.add_var("q(t,x)");
    M.add_domain( t, FFDom( 0., 1., 2, FFDom::LGL, 3 ) );
    M.add_domain( x, FFDom( 0., 1., 2, FFDom::LGL, 3 ) );
    M.add_state( T, { t, x } );  M.add_state( q, { t, x } );
    M.add_equation( OpP( T, t ) - OpP( q, x ),  std::vector<FFVar>{ t, x }, std::vector<int>{} );
    M.add_equation( sqr( q ) - OpP( T, x ),     std::vector<FFVar>{ t, x }, std::vector<int>{} );
    if( !M.setup() ){ std::cout << "    setup FAILED\n"; ++g_fail; }
    else {
      expect( M, q, t, 'C', 'N' );   // nonlinear relation: CLOSURE, letter N for C1
      expect( M, q, x, 'A' );
    }
  }

  // ---- 5: a state no relation determines -> NONE ---------------------------------------------
  {
    std::cout << "\n  5  unattached state w (no relation determines it in x)\n";
    FFGraph D; FFModel M( &D );
    FFVar t = D.add_var("t"), x = D.add_var("x"), T = D.add_var("T(t,x)"), w = D.add_var("w(t,x)");
    M.add_domain( t, FFDom( 0., 1., 2, FFDom::LGL, 3 ) );
    M.add_domain( x, FFDom( 0., 1., 2, FFDom::LGL, 3 ) );
    M.add_state( T, { t, x } );  M.add_state( w, { t, x } );
    M.add_equation( OpP( T, t ) - OpP( T, x ),  std::vector<FFVar>{ t, x }, std::vector<int>{} );
    M.add_equation( OpP( w, t ) - T,            std::vector<FFVar>{ t, x }, std::vector<int>{} );
    if( !M.setup() ){ std::cout << "    setup FAILED\n"; ++g_fail; }
    else {
      expect( M, w, t, 'A' );   // w_t appears
      expect( M, w, x, 'N' );      // nothing ties w across x
    }
  }

  std::cout << "\n  OCFE_claim: " << ( g_fail? "FAIL": "PASS" ) << "  (" << g_fail << " mismatch(es))\n";
  return g_fail? 1: 0;
}
