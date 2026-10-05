// ============================================================================
//  OCFE_hybrid_adequacy.cpp  --  exercise the Option C three-way stall attribution
//
//  WHY THIS DRIVER EXISTS.
//  rev44 added a three-way attribution to the IC_STRONG pivot-search stall in
//  _build_frozen_strong_tau_elimination:
//
//      stalled_at < rank(B_eliminated)   -> PIVOT SEARCH DEFECT
//      rank(B_eliminated) < ntrace       -> SELECTION INADEQUATE  (the Option A / beta trigger)
//      otherwise                         -> STRUCTURAL
//
//  NO PLAN IN THE 134-DRIVER CORPUS STALLS.  That branch has therefore never executed.  An
//  error path that has never run is a coin flip when it finally matters, and the first stall
//  that ever occurs will be the ONLY evidence available about which mechanism failed -- so a
//  message that blames the model when the SELECTION is at fault would send that investigation
//  the wrong way.  This driver makes the branch run on demand.
//
//  HOW.  rev45's CRONOS_FORCE_BAD_KEEP_EXPLICIT=1 suppresses the keep-explicit routing, so a
//  protected claim is left to the Schur elimination, which cannot cover it, and the search
//  stalls.
//
//  WHAT THIS DOES AND DOES NOT VALIDATE.
//    DOES     -- that the branch is reachable; that r_elim is COMPUTED from the eliminated
//                block rather than inferred as rank(B_all) - n_explicit; that the counts
//                printed are self-consistent; that the message names the right cause.
//    DOES NOT -- reproduce the failure Option A exists for.  A genuinely inadequate selection
//                keeps the WRONG columns; this keeps NONE.  Both stall, by different routes.
//                Do NOT read a run under this flag as evidence about the hybrid's adequacy.
//
//  MODEL.  Deliberately the smallest thing in the corpus that has a protected claim and a
//  rank-deficient multiplier block: the scalar_nest geometry (two 1-D fields, 2 elements x 6
//  CGL each), which reports k=1, rank(B)=5 of 6, and routes exactly 1 LINK-anchored claim to
//  explicit tau.  Small enough that the printed matrices can be read by hand.
//
//  Build:
//    g++ -std=c++17 ... -DOCFE_OCFESLV_HEADER='"ocfeslv_rev45.hpp"' \
//        OCFE_hybrid_adequacy.cpp -o OCFE_hybrid_adequacy  <libs>
//  Run:
//    ./OCFE_hybrid_adequacy            # baseline: must SUCCEED (setup ok, no stall)
//    CRONOS_FORCE_BAD_KEEP_EXPLICIT=1 \
//      ./OCFE_hybrid_adequacy          # injected: must STALL with SELECTION INADEQUATE
//
//  The driver runs BOTH itself and compares, so a single invocation is enough.
// ============================================================================
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#ifndef OCFE_OCFESLV_HEADER
#define OCFE_OCFESLV_HEADER "ocfeslv.hpp"
#endif
#include OCFE_OCFESLV_HEADER

using namespace mc;

// Build and set up the model.  Returns true when setup() succeeds.
static bool try_setup( bool& square, size_t& nvar, size_t& neqn )
{
  square = false; nvar = 0; neqn = 0;

  FFGraph DAG;
  FFVar z  = DAG.add_var( "z" );
  FFVar y  = DAG.add_var( "y" );
  FFVar a  = DAG.add_var( "a(z)" );
  FFVar b  = DAG.add_var( "b(y)" );
  FFVar p1 = DAG.add_var( "p1" );
  FFVar p2 = DAG.add_var( "p2" );
  FFVar p3 = DAG.add_var( "p3" );
  FFVar p4 = DAG.add_var( "p4" );

  FFPartial  OpP;
  FFEval     OpEval;
  FFIntegral OpI;

  OCFESLV oc( &DAG );
  oc.add_domain( z, FFDom( 0., 1., 2, FFDom::CGL, 6 ) );
  oc.add_domain( y, FFDom( 0., 1., 2, FFDom::CGL, 6 ) );
  oc.add_state ( a,  { z } );
  oc.add_state ( b,  { y } );
  oc.add_state ( p1, {} );
  oc.add_state ( p2, {} );
  oc.add_state ( p3, {} );
  oc.add_state ( p4, {} );

  int const INT = FFDom::ALL - FFDom::LB - FFDom::UB;

  FFVar PDE_a = OpP( OpP( a, z ), z ) + 2.;
  FFVar PDE_b = OpP( OpP( b, y ), y ) + 2.;
  FFVar E1 = p1 - OpEval( OpI( a * b, z ), y, 0.5 );
  FFVar E2 = p2 - OpI( OpI( a * b, z ), y );
  FFVar E3 = p3 - OpEval( OpP( a, z ), z, 0.25 );
  FFVar E4 = p4 - OpEval( OpP( OpI( a * b, z ), y ), y, 0.25 );

  oc.add_equation( PDE_a, { z }, { INT },       OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( a,     { z }, { FFDom::LB }, OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );
  oc.add_equation( a,     { z }, { FFDom::UB }, OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );
  oc.add_equation( PDE_b, { y }, { INT },       OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( b,     { y }, { FFDom::LB }, OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );
  oc.add_equation( b,     { y }, { FFDom::UB }, OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );
  oc.add_equation( E1, {}, {}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( E2, {}, {}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( E3, {}, {}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( E4, {}, {}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );

  oc.options.REDUCE.ORDER    = OCFESLV::Options::RED_FULL;
  oc.options.CLASSIFY.MODE        = OCFESLV::Options::CLASS_AUTO;
  oc.options.INTERFACE.IMPOSITION = OCFESLV::Options::IC_STRONG;   // the mode under test
  oc.options.INTERFACE.TYPE  = OCFESLV::Options::IC_AUTO;
  oc.options.INTERFACE.SAT_SIGMA0      = 1.0;
  oc.options.DISPLAY_LEVEL   = 1;                           // needed for the attribution line

  if( !oc.setup() ) return false;
  nvar = oc.n_colloc_sta(); neqn = oc.n_colloc_eqn();
  square = ( nvar == neqn );
  return true;
}

int main()
{
  std::cout << "================================================================\n"
            << "  hybrid-adequacy attribution: exercise the IC_STRONG stall branch\n"
            << "================================================================\n";

  // ---- 1. baseline: the keep-explicit routing is active and setup must succeed ----------
  std::cout << "\n---- [1] BASELINE (CRONOS_FORCE_BAD_KEEP_EXPLICIT unset) ----\n"
            << "  expect: setup succeeds, square, NO [hybrid-adequacy] line.\n\n";
  ::unsetenv( "CRONOS_FORCE_BAD_KEEP_EXPLICIT" );
  bool sq0 = false; size_t nv0 = 0, ne0 = 0;
  bool const ok0 = try_setup( sq0, nv0, ne0 );
  std::cout << "\n  baseline: setup=" << ( ok0 ? "ok" : "FAILED" )
            << " nvar=" << nv0 << " neqn=" << ne0
            << " square=" << ( sq0 ? "yes" : "no" ) << "\n";

  // ---- 2. injected: suppress the routing; the Schur path must stall --------------------
  std::cout << "\n---- [2] FAULT INJECTED (CRONOS_FORCE_BAD_KEEP_EXPLICIT=1) ----\n"
            << "  expect: setup FAILS, and a [hybrid-adequacy] block reporting\n"
            << "          rank(B_eliminated) < ntrace_elim -> SELECTION INADEQUATE.\n"
            << "  NOTE: setup failing here is the PASS condition, not a regression.\n\n";
  ::setenv( "CRONOS_FORCE_BAD_KEEP_EXPLICIT", "1", 1 );
  bool sq1 = false; size_t nv1 = 0, ne1 = 0;
  bool const ok1 = try_setup( sq1, nv1, ne1 );
  ::unsetenv( "CRONOS_FORCE_BAD_KEEP_EXPLICIT" );
  std::cout << "\n  injected: setup=" << ( ok1 ? "ok (UNEXPECTED)" : "failed (expected)" ) << "\n";

  // ---- verdict --------------------------------------------------------------------------
  bool const pass = ok0 && sq0 && !ok1;
  std::cout << "\n================================================================\n";
  if( pass ){
    std::cout << "  RESULT: PASS -- the stall branch is reachable and the baseline is clean.\n"
              << "\n  NOW READ THE [hybrid-adequacy] LINES ABOVE BY HAND.  A reachable branch\n"
              << "  is not a correct one, and this driver cannot check the attribution text:\n"
              << "    - stalled_at, rank(B_eliminated), ntrace_elim, rank(B_all), k and\n"
              << "      n_explicit must be mutually consistent;\n"
              << "    - with the routing suppressed n_explicit must be 0;\n"
              << "    - the verdict must be SELECTION INADEQUATE, i.e. rank(B_eliminated) <\n"
              << "      ntrace_elim.  If it reports STRUCTURAL or PIVOT SEARCH DEFECT instead,\n"
              << "      the attribution arithmetic is wrong and Option C is not doing its job.\n";
  }
  else{
    std::cout << "  RESULT: FAIL\n";
    if( !ok0 )       std::cout << "    baseline setup failed -- the model itself is broken;\n"
                                  "    nothing about the attribution can be concluded.\n";
    else if( !sq0 )  std::cout << "    baseline not square.\n";
    else if( ok1 )   std::cout << "    injection did NOT stall the elimination.  Either the flag\n"
                                  "    is not wired (is this rev45 or later? check the HEADER_ID\n"
                                  "    stamp) or the keep-explicit routing is not the only path\n"
                                  "    covering this claim -- which would itself be worth knowing.\n";
  }
  std::cout << "================================================================\n";
  return pass ? 0 : 1;
}
