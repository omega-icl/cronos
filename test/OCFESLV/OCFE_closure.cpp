// Copyright (C) Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the EPL-2.0 with GPL-2.0-or-later as a Secondary License; see the LICENSE file.
// SPDX-License-Identifier: EPL-2.0 OR GPL-2.0-or-later
//
// ============================================================================
//  OCFE_closure.cpp -- manufacture a CLOSURE-tagged state that RECEIVES
//  forced-rescue edges, to exercise the C0-S bucket of the rev226 probe.
//
//  The corpus never produces one: on MBC6 every recoverable forced edge is
//  C0-L (a constant coefficient), and the 16 CLOSURE tags sit on states that
//  receive no forced edges.  So "a state-dependent coupling never arises" is
//  unproven, not established -- this driver tries to construct the case.
//
//  Model: nonlinear-conduction heat equation in mixed form, on 3 elements,
//         T_t = q_x,  q = (1 + T) * T_x
//  q is determined pointwise by a relation whose coefficient on T_x depends on
//  the state T, so C0 must tag q CLOSURE (letter S).  With order reduction the
//  relation becomes q = (1 + T) * Dx_T, and Dx_T is an aux-algebraic receiver --
//  exactly the population that takes the forced rescue.
//
//  Run with CRONOS_RESCUE_C0_PROBE=1 CRONOS_RESCUE_ROWCOEF=1.
//  PASS here means only "the case was constructed and classified"; the probe
//  census is the actual output.
// ============================================================================
#include <iostream>
#include <iomanip>
#include <vector>
#include "ffunc.hpp"
#include OCFE_OCFESLV_HEADER
using namespace mc;

int main()
{
  std::cout << "OCFESLV ** header: " << OCFESLV::HEADER_ID << std::endl;
  FFGraph DAG; OCFESLV oc( &DAG );

  FFVar t = DAG.add_var("t"), x = DAG.add_var("x");
  FFVar T = DAG.add_var("T(t,x)"), q = DAG.add_var("q(t,x)");
  oc.add_domain( t, FFDom( 0., 0.1, 1, FFDom::LGR, 4 ) );
  oc.add_domain( x, FFDom( 0., 1., 3, FFDom::LGL, 4 ) );
  oc.add_state( T, { t, x } );
  oc.add_state( q, { t, x } );

  FFPartial OpP;
  oc.add_equation( OpP( T, t ) - OpP( q, x ), std::vector<FFVar>{ t, x }, std::vector<int>{} );
  oc.add_equation( q - ( 1. + T ) * OpP( T, x ), std::vector<FFVar>{ t, x }, std::vector<int>{} );
  oc.add_equation( T - 1., std::vector<FFVar>{ t, x }, std::vector<int>{ FFDom::LB, 0 },
                   OCFESLV::EqnOptions( OCFESLV::EqnRole::INITIAL, 0 ) );
  oc.add_equation( T - 1., std::vector<FFVar>{ t, x }, std::vector<int>{ 0, FFDom::LB },
                   OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );
  oc.add_equation( q - 0., std::vector<FFVar>{ t, x }, std::vector<int>{ 0, FFDom::UB },
                   OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );
  oc.update_ref( T, 1.0 );
  oc.update_ref( q, 0.1 );

  oc.options.INTERFACE.IMPOSITION = OCFESLV::Options::IC_STRONG;
  oc.options.REDUCE.ORDER    = OCFESLV::Options::RED_FULL;
  oc.options.DISPLAY_LEVEL   = 1;

  bool const ok = oc.setup();
  std::cout << "  setup: " << ok << std::endl;
  if( ok ){
    oc.claim_report( std::cout );
    std::vector<double> xv, inp;
    if( oc.init( xv, inp, nullptr ) ){
      OCFESLV::SolveReport sr = oc.solve( xv.data(), inp.data(), nullptr );
      double m = 0.; for( double v : xv ) m = std::max( m, std::abs( v ) );
      std::cout << "  solve: converged=" << sr.converged << " max|x|=" << std::setprecision(8) << m << std::endl;
    }
  }
  // ---- CASE 2: an ALIAS of an auxiliary whose relation has a STATE-DEPENDENT coefficient ----
  // The forced rescue reaches a state either as a reduction auxiliary (whose LINK row is
  // aux - d(parent)/d(dir), coefficient 1 by construction -- always C0-L) or as an ALIAS:
  // a row holding the state plus EXACTLY ONE other state.  The alias route dropped its
  // FFDep-linearity gate (fix5), so the coefficient there need NOT be constant.  That is the
  // only way a CLOSURE state can take a forced edge, and this case constructs it:
  //     s - (1 + Dx_T) * Dx_T = 0
  // one other state (Dx_T, appearing in both the coefficient and the factor), non-constant
  // coefficient.  If the probe reports C0-S here, the S branch is REQUIRED; if the rescue
  // never fires, a documented refusal suffices.
  std::cout << "\n  --- case 2: alias with a state-dependent coefficient ---" << std::endl;
  {
    FFGraph D2; OCFESLV oc2( &D2 );
    FFVar t2 = D2.add_var("t"), x2 = D2.add_var("x");
    FFVar T2 = D2.add_var("T(t,x)"), s2 = D2.add_var("s(t,x)");
    oc2.add_domain( t2, FFDom( 0., 0.1, 1, FFDom::LGR, 4 ) );
    oc2.add_domain( x2, FFDom( 0., 1., 3, FFDom::LGL, 4 ) );
    oc2.add_state( T2, { t2, x2 } );
    oc2.add_state( s2, { t2, x2 } );
    FFPartial OpQ;
    // T_t = T_xx  (reduction mints Dx_T);  s aliases Dx_T through a non-constant coefficient
    oc2.add_equation( OpQ( T2, t2 ) - OpQ( OpQ( T2, x2 ), x2 ), std::vector<FFVar>{ t2, x2 }, std::vector<int>{} );
    oc2.add_equation( s2 - ( 1. + OpQ( T2, x2 ) ) * OpQ( T2, x2 ), std::vector<FFVar>{ t2, x2 }, std::vector<int>{} );
    oc2.add_equation( T2 - 1., std::vector<FFVar>{ t2, x2 }, std::vector<int>{ FFDom::LB, 0 },
                      OCFESLV::EqnOptions( OCFESLV::EqnRole::INITIAL, 0 ) );
    oc2.add_equation( T2 - 1., std::vector<FFVar>{ t2, x2 }, std::vector<int>{ 0, FFDom::LB },
                      OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );
    oc2.add_equation( OpQ( T2, x2 ) - 0., std::vector<FFVar>{ t2, x2 }, std::vector<int>{ 0, FFDom::UB },
                      OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );
    oc2.update_ref( T2, 1.0 );  oc2.update_ref( s2, 0.5 );
    oc2.options.INTERFACE.IMPOSITION = OCFESLV::Options::IC_STRONG;
    oc2.options.REDUCE.ORDER    = OCFESLV::Options::RED_FULL;
    oc2.options.DISPLAY_LEVEL   = 1;
    bool const ok2 = oc2.setup();
    std::cout << "  setup: " << ok2 << std::endl;
    if( ok2 ) oc2.claim_report( std::cout );
  }

  std::cout << "\n  OCFE_closure: " << ( ok? "constructed": "SETUP FAILED" ) << std::endl;
  return ok? 0: 1;
}
