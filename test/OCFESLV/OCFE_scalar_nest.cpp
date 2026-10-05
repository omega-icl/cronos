// ============================================================================
//  OCFE_scalar_nest2.cpp  --  reduction/derivative nesting coverage
//
//  Two elliptic 1-D fields (spatial, no evolution direction):
//      a(z) = z(1-z)  ( INT a = 1/6 , a'(z) = 1-2z )
//      b(y) = y(1-y)  ( INT b = 1/6 )
//  Nested reductions/derivatives on the product a*b (distributed over {z,y}):
//
//    p1 = OpEval( OpI(a*b,z), y, 0.5 )   = b(0.5)/6 = 1/24   [ integral-inside-point ]
//    p2 = OpI( OpI(a*b,z), y )           = (1/6)^2  = 1/36   [ double integral ]
//    p3 = OpEval( OpP(a,z), z, 0.25 )    = a'(0.25) = 0.5    [ OpEval OF a derivative ]
//    p4 = OpEval( OpP(OpI(a*b,z),y), y, 0.25 ) = b'(0.25)/6 = 1/12  [ OpP over an integral ]
//
//  p1/p2 exercise the distributed integral LINK (free_dom-for-distributed);
//  p3/p4 exercise OpP-inside-a-reduction, which requires materialising the OpP
//  operand as a state (FFIntegral/FFEval/FFPartial have no deriv()).
//
//  Build:
//    g++ -std=c++17 ... -DOCFE_OCFESLV_HEADER='"ocfeslv_nest5.hpp"'
//        OCFE_scalar_nest2.cpp -o OCFE_scalar_nest2  <libs>
// ============================================================================
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
#include <limits>
#include <utility>
#include <map>
#include <cstdlib>

#ifndef OCFE_OCFESLV_HEADER
#define OCFE_OCFESLV_HEADER "ocfeslv.hpp"
#endif
#include OCFE_OCFESLV_HEADER

using namespace mc;

static bool g_ok = true;
static void check( char const* nm, double got, double want, double tol )
{
  double const err = std::fabs( got - want );
  bool const pass = ( err < tol );
  g_ok &= pass;
  std::cout << "  " << std::left << std::setw(34) << nm
            << std::right << std::scientific << std::setprecision(8)
            << " got=" << std::setw(15) << got << " want=" << std::setw(15) << want
            << std::setprecision(2) << "  err=" << std::setw(9) << err
            << "  " << ( pass ? "PASS" : "FAIL" ) << "\n";
}

// ---------------------------------------------------------------------------------------
//  Imposition-mode comparison
// ---------------------------------------------------------------------------------------
// The rev43 corpus sweep flagged this driver: under IC_STRONG, converged,
//     Dz_a(z) 2.500e-01   Dy_b(y) 2.500e-01   Intz_a(z) 1.624e-02   Dpy_Intz_a(z) 1.389e-02
// while a(z)=0 and b(y)=1.943e-16.  The [claims] readout shows ALL SIX of those states carry
// a retained continuity claim, so this is not a case of measuring something never claimed:
//     [claims] states WITH a retained claim: a(z) b(y) Dz_a(z) Dy_b(y) Intz_a(z) Dpy_Intz_a(z)
//     [claims] states WITHOUT any claim:     p1 p2 p3 p4 Evz_* Evy_* Inty_*
// and every state WITHOUT a claim sits at exactly 0.
//
// The violated four are exactly the states materialised to give FFPartial/FFIntegral an
// operand (those operators have no deriv()).  The clean two are the declared fields.  So the
// question this run answers is whether continuity of an OPERATOR-MATERIALISED state is
// enforced under any imposition mode, or under none.
//
// This driver is small (2x6 CGL on each of two 1-D domains, 6 claims total, k=1) which makes
// it the cheapest place in the corpus to settle that.
// ---------------------------------------------------------------------------------------
//  Per-node dump at coincidence groups  --  the aliasing discriminator
// ---------------------------------------------------------------------------------------
// The spread alone cannot distinguish two very different situations:
//
//   (A) the continuity claim is INERT -- the two copies hold two plausible but different
//       values of the quantity, and nothing is forcing them together;
//   (B) the reported "duplicate group" is an INDEXING artifact -- node_colloc() for a
//       materialised operand state does not line up with the var[] walk
//       (off += node_colloc(st).size()), so the instrument is differencing unrelated entries.
//
// (B) is a live possibility because the numbers are suspiciously structured: a(z)=z(1-z) is
// QUADRATIC and exactly representable in degree-5 CGL, so a'(0.5)=0 should be exact from both
// sides and the spread should be ZERO -- yet it reads exactly 0.25, which is a(0.5).  And
// 1.388889e-02 is exactly 1/72.
//
// So print the values themselves, next to the analytic a(z) and a'(z) at the same coordinate.
// If a node of Dz_a holds a(z) instead of a'(z), that is (B) and there is no continuity
// failure at all.  If both copies hold plausible derivative values that simply differ, it is
// (A) and the claim is not being enforced.
static void dump_coincidence( OCFESLV const& oc, FFVar const& V, std::vector<double> const& xv,
                              char const* mode_name )
{
  // v3: compute the offset BOTH ways and print both.
  //   walk_off  -- the old manual walk `off += node_colloc(st).size()`.  Adds 0 for a SCALAR
  //                state (node_colloc returns an empty vector) even though the state owns a
  //                var[] slot, so it under-counts by the number of preceding scalar states.
  //   true_off  -- pos_state(V, all-domains-zero), which is block 0 for scalars and
  //                distributed states alike, and is correct for every state.
  // If they differ, the shift IS the bug, quantified, on the same line as the values it
  // corrupted.  Keeping the wrong one visible is deliberate: this is the diagnostic that
  // identified the fault, and a rerun that shows walk_off == true_off on a driver with no
  // scalar states is itself a useful negative control.
  size_t walk_off = 0; bool found = false;
  for( auto const& st : oc.states_colloc() ){
    if( st.id() == V.id() ){ found = true; break; }
    walk_off += oc.node_colloc( st ).size();
  }
  if( !found ){ std::cout << "  [nodes] " << V.name() << ": NOT FOUND in states_colloc()\n"; return; }
  // v4: pos_state() does ndx_el.at(d) for EVERY domain of the state, so an EMPTY index map
  // throws std::out_of_range for any DISTRIBUTED state -- it is only valid for scalars.  v3
  // passed an empty map here and aborted after the (correct) header report had printed, which
  // made a working fix look like it was still broken.  Pin every domain of the state to
  // element 0 instead: that is block 0 for scalars and distributed states alike.
  //
  // The header has a private _state_base_offset() doing exactly this; drivers must rebuild it
  // from the public var_state() accessor until it is exposed.
  std::map<FFVar,size_t,lt_FFVar> ndx0;
  {
    auto const& vs = oc.var_state();
    auto const it = vs.find( V );
    if( it != vs.end() )
      for( auto const& d : it->second ) ndx0[d] = 0;
  }
  size_t const off = oc.pos_state( V, ndx0 );

  std::vector<std::vector<double>> const nodes = oc.node_colloc( V );
  std::map<std::vector<long long>, std::vector<size_t>> groups;
  for( size_t i=0; i<nodes.size(); ++i ){
    std::vector<long long> key;
    for( double c : nodes[i] ) key.push_back( (long long)std::llround( c * 1.0e12 ) );
    groups[key].push_back( i );
  }

  std::cout << "  [nodes] " << mode_name << "  " << V.name()
            << "  base_off=" << off << " (pos_state)"
            << "  walk_off=" << walk_off
            << ( walk_off == off ? "  [agree]" : "  [** DISAGREE: manual walk is short **]" )
            << "  n_node=" << nodes.size() << "\n";
  for( auto const& kv : groups ){
    if( kv.second.size() < 2 ) continue;
    for( size_t idx : kv.second ){
      double const val = ( off+idx < xv.size() ) ? xv[off+idx] : std::numeric_limits<double>::quiet_NaN();
      std::cout << "  [nodes]     i=" << std::setw(4) << idx
                << "  var[" << std::setw(4) << (off+idx) << "]"
                << "  coord=(";
      for( size_t c=0; c<nodes[idx].size(); ++c )
        std::cout << ( c? "," : "" ) << std::fixed << std::setprecision(6) << nodes[idx][c];
      std::cout << ")  value=" << std::scientific << std::setprecision(10) << val;
      // Analytic reference at this coordinate: a(z)=z(1-z), a'(z)=1-2z, same for b(y).
      if( nodes[idx].size() == 1 ){
        double const w = nodes[idx][0];
        std::cout << "   [ w(1-w)=" << std::setprecision(6) << w*(1.0-w)
                  << "  1-2w=" << (1.0-2.0*w) << " ]";
      }
      std::cout << "\n";
    }
  }
}

struct ImpMode { OCFESLV::Options::ImpositionType t; char const* name; };
static ImpMode const IMP_MODES[3] = {
  { OCFESLV::Options::IC_WEAK,   "IC_WEAK"   },
  { OCFESLV::Options::IC_TRACE,  "IC_TRACE"  },
  { OCFESLV::Options::IC_STRONG, "IC_STRONG" }
};

struct ModeResult {
  char const* name = "";
  bool   setup_ok = false, converged = false;
  double final_r = 0.;
  double worst = -1.;
  std::vector<std::pair<std::string,double>> per_state;
  // CRASH FIX (v2).  This was `std::vector<double> pval{ 4, NaN };`, which selects the
  // INITIALIZER-LIST constructor and builds a TWO-element vector {4.0, NaN} -- not four NaNs.
  // Writing pval[2] and pval[3] then corrupted the heap and the run died with
  // "munmap_chunk(): invalid pointer" after IC_TRACE, losing IC_STRONG and the whole
  // comparison table.  A plain array removes the ambiguity entirely.
  double pval[4] = { std::numeric_limits<double>::quiet_NaN(),
                     std::numeric_limits<double>::quiet_NaN(),
                     std::numeric_limits<double>::quiet_NaN(),
                     std::numeric_limits<double>::quiet_NaN() };
};

static ModeResult run_mode( ImpMode const& mode )
{
  ModeResult MR; MR.name = mode.name;
  std::cout << "\n---------------- imposition = " << mode.name << " ----------------\n";

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

  oc.add_equation( PDE_a, { z }, { INT }, OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( a, { z }, { FFDom::LB }, OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );
  oc.add_equation( a, { z }, { FFDom::UB }, OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );
  oc.add_equation( PDE_b, { y }, { INT }, OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( b, { y }, { FFDom::LB }, OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );
  oc.add_equation( b, { y }, { FFDom::UB }, OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );
  oc.add_equation( E1, {}, {}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( E2, {}, {}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( E3, {}, {}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( E4, {}, {}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );

  oc.options.REDUCE.ORDER    = OCFESLV::Options::RED_FULL;
  oc.options.CLASSIFY.MODE        = OCFESLV::Options::CLASS_AUTO;
  oc.options.INTERFACE.IMPOSITION = mode.t;
  oc.options.INTERFACE.TYPE  = OCFESLV::Options::IC_AUTO;
  oc.options.INTERFACE.SAT_SIGMA0      = 1.0;
  oc.options.DISPLAY_LEVEL   = 1;

  if( !oc.setup() ){ std::cerr << "  ERROR: setup() failed under " << mode.name << "\n"; return MR; }
  MR.setup_ok = true;
  std::cout << "  setup: states=" << oc.n_colloc_sta()
            << " equations=" << oc.n_colloc_eqn()
            << " square=" << ( oc.n_colloc_sta() == oc.n_colloc_eqn() ? "yes" : "NO" ) << "\n";
  if( oc.n_colloc_sta() != oc.n_colloc_eqn() ){ std::cerr << "  ERROR: not square\n"; return MR; }

  std::vector<double> varInit, inpInit;
  if( !oc.init( varInit, inpInit, nullptr ) ){ std::cerr << "  ERROR: init() failed\n"; return MR; }

  std::vector<double> xv = varInit;
  OCFESLV::SolveReport rep = oc.solve( xv.data(), nullptr, nullptr );
  MR.converged = rep.converged; MR.final_r = rep.final_residual;
  std::cout << "  solve: converged=" << ( rep.converged ? "yes" : "no" )
            << " iters=" << rep.iterations
            << " final|r|=" << std::scientific << std::setprecision(3) << rep.final_residual << "\n";

  // ---- per-state duplicate-node spreads ------------------------------------------------
  // Captured even when the solve fails: a non-converged spread is structural information,
  // and losing the row for the one mode that does not converge would be the worst possible
  // gap in a three-way comparison.  Labelled with conv so it cannot be read as accuracy.
  MR.worst = oc.report_duplicate_node_spreads(
               xv, ( std::string("scalar_nest ") + mode.name
                     + ( rep.converged ? "  conv=yes" : "  conv=NO -- STRUCTURAL ONLY" ) ).c_str(),
               &std::cout );
  {
    // v3: let the HEADER derive each state's base offset (rev44 defaults `off` to
    // pos_state).  v2 passed a manually-walked offset here -- `off += node_colloc().size()`
    // -- which adds 0 for a SCALAR state even though the state owns a var[] slot, so every
    // state after p1..p4 was read 4 slots short.  That is the bug rev44 fixes in the header;
    // passing our own offset would have kept the driver on the broken path and made a
    // correctly-fixed header look unchanged.
    for( auto const& st : oc.states_colloc() )
      MR.per_state.emplace_back( st.name(),
                                 oc.duplicate_node_spread( st, xv ).max_pair_spread );
  }

  // ---- the aliasing discriminator ------------------------------------------------------
  // a(z) is the CONTROL: it is a declared field, its spread is at round-off, and its nodal
  // values must equal z(1-z).  If the walk is sound for a(z) but Dz_a's nodes hold z(1-z)
  // rather than 1-2z, the fault is indexing, not continuity.  Dz_a is the suspect: exactly
  // representable, so BOTH copies should read 1-2z = 0 at the interface z=0.5.
  std::cout << "\n  ---- per-node values at coincidence groups (" << mode.name << ") ----\n";
  dump_coincidence( oc, a,    xv, mode.name );   // control: declared field, spread ~1e-16
  dump_coincidence( oc, b,    xv, mode.name );   // control
  {
    // The materialised operands are not driver-side FFVars, so reach them by name.
    for( auto const& st : oc.states_colloc() ){
      std::string const nm = st.name();
      if( nm.rfind("Dz_a",0)==0 || nm.rfind("Dy_b",0)==0
       || nm.rfind("Intz_a",0)==0 || nm.rfind("Dpy_Intz_a",0)==0 )
        dump_coincidence( oc, st, xv, mode.name );
    }
  }

  if( !rep.converged ) return MR;

  OCFESLV::t_Coord sp;
  MR.pval[0] = oc.eval_colloc<double>( p1, sp, xv.data(), nullptr, nullptr );
  MR.pval[1] = oc.eval_colloc<double>( p2, sp, xv.data(), nullptr, nullptr );
  MR.pval[2] = oc.eval_colloc<double>( p3, sp, xv.data(), nullptr, nullptr );
  MR.pval[3] = oc.eval_colloc<double>( p4, sp, xv.data(), nullptr, nullptr );

  std::cout << "\n  checks (" << mode.name << "):\n";
  check( "p1 OpEval(OpI(a*b,z),y,0.5)==1/24", MR.pval[0], 1.0/24.0, 1e-8 );
  check( "p2 OpI(OpI(a*b,z),y)      ==1/36", MR.pval[1], 1.0/36.0, 1e-8 );
  check( "p3 OpEval(OpP(a,z),z,0.25)==0.5",  MR.pval[2], 0.5,      1e-8 );
  check( "p4 OpEval(OpP(OpI),y,0.25)==1/12", MR.pval[3], 1.0/12.0, 1e-8 );
  return MR;
}

int main()
{
  std::cout << "================================================================\n"
            << "  reduction/derivative nesting coverage\n"
            << "  + imposition-mode comparison of duplicate-node spreads\n"
            << "================================================================\n";

  ModeResult MR[3];
  for( int m=0; m<3; ++m ) MR[m] = run_mode( IMP_MODES[m] );

  // ---- side-by-side, per state ---------------------------------------------------------
  // Per STATE, not just the worst: the whole point is that a(z)/b(y) are clean while the
  // operator-materialised states are not, and a single worst-case number hides that.
  std::cout << "\n================================================================\n"
            << "  duplicate-node spread by state and imposition mode\n"
            << "================================================================\n";
  std::cout << "  " << std::left << std::setw(24) << "state";
  for( int m=0; m<3; ++m ) std::cout << std::right << std::setw(13) << MR[m].name;
  std::cout << "\n";
  size_t const nst = MR[2].per_state.size();
  for( size_t i=0; i<nst; ++i ){
    std::cout << "  " << std::left << std::setw(24) << MR[2].per_state[i].first;
    for( int m=0; m<3; ++m ){
      if( i < MR[m].per_state.size() )
        std::cout << std::right << std::scientific << std::setprecision(3)
                  << std::setw(13) << MR[m].per_state[i].second;
      else
        std::cout << std::right << std::setw(13) << "-";
    }
    std::cout << "\n";
  }
  std::cout << "  " << std::left << std::setw(24) << "converged";
  for( int m=0; m<3; ++m ) std::cout << std::right << std::setw(13) << ( MR[m].converged ? "yes" : "NO" );
  std::cout << "\n";

  std::cout << "\n  READ THIS AS:\n"
            << "    a(z), b(y) are the DECLARED fields; the rest are materialised operands for\n"
            << "    FFPartial/FFIntegral.  All six of a, b, Dz_a, Dy_b, Intz_a, Dpy_Intz_a carry a\n"
            << "    retained continuity claim -- check the [claims] block above to confirm it for\n"
            << "    each mode rather than assuming it is mode-independent.\n"
            << "      operator states dirty under IC_STRONG only -> Schur/keep-explicit path;\n"
            << "                                                    Stage 2 work.\n"
            << "      dirty under all three modes                -> claim generation for\n"
            << "                                                    materialised operands.\n"
            << "      clean under IC_TRACE only                  -> the projection is load-bearing\n"
            << "                                                    and Stage 2 inherits it.\n"
            << "    p1..p4 PASS in every mode is NOT evidence the continuity is fine: those are\n"
            << "    reductions, which average over the jump.  This is the k/SYMMETRIC/residual/PASS\n"
            << "    trap that Stage A passed while destroying MBC5 by 24 orders.\n";

  bool any_conv = false;
  for( int m=0; m<3; ++m ) any_conv |= MR[m].converged;
  std::cout << "\n  RESULT: " << ( g_ok && any_conv ? "ALL PASS" : "FAIL" ) << "\n";
  return ( g_ok && any_conv ) ? 0 : 1;
}
