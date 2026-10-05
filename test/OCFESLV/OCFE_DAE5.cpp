// OCFE_DAE5.cpp
// ===========================================================================
//  A MINIMAL, WELL-POSED INDEX-1 DAE WHOSE DIFFERENTIAL CORE IS RECTANGULAR
// ===========================================================================
//
//  WHY THIS DRIVER EXISTS
//  ----------------------
//  OCFE_MMPDE27 classifies as UNDETERMINED with sym=3x4: three retained differential
//  rows differentiating four states, so A_e is 3x4 and neither branch of the
//  causality test (A_e nonsingular -> ODE, singular -> DAE) can be taken.  Rank and
//  singularity are only defined on a square matrix.
//
//  That was easy to dismiss as an ALE/moving-mesh peculiarity.  It is not.  Benoit's
//  observation: the same shape arises from an ORDINARY DAE written on a SUM --
//
//      dx1/dt + dx2/dt = ...        one row, TWO differentiated states
//      x1 + x2         = ...        no derivative -> routed to vAlgEqn
//
//  giving vState = {x1,x2} (2 columns) and diff_eqn = {Eq1} (1 row): sym=1x2.
//
//  THE CRITERION IS NOT "moving mesh" OR "ALE".  It is: does any differential row
//  differentiate MORE THAN ONE state, without another row to balance it?  One row
//  carrying two derivatives is enough.  OCFE_DAE0 escapes only because each of its
//  four evolution rows differentiates exactly one state (x'=u, y'=v, u'=..., v'=...),
//  and lam is never differentiated so holds no column -- 4x4, square, classifiable.
//
//  This shape is COMMON: total mass or energy balances written on a sum, Kirchhoff
//  current laws, any conservation law imposed on an aggregate while a constraint fixes
//  the split.  So DIFFERENTIAL_RECTANGULAR should be expected well beyond ALE models.
//
//  MODEL  (deliberately as small as a DAE can be)
//  ---------------------------------------------
//      SUM  :  dx1/dt + dx2/dt = -(x1 + x2)          -> s' = -s,  s = x1+x2
//      SPLIT:  x1 - x2         = d0 * exp(-2 t)      -> d      = d0 e^{-2t}
//
//  s(0) = s0 given, so  s(t) = s0 e^{-t}  and  x1 = (s+d)/2, x2 = (s-d)/2, both in
//  closed form.  WELL-POSED and index-1: SPLIT differentiated once gives d', which with
//  SUM determines both x1' and x2'.  Note the NAIVE version x1+x2 = g(t) is structurally
//  SINGULAR -- differentiating it collides with SUM unless f = g' -- so the difference
//  form is used instead.
//
//  Only x1 needs an initial condition: SPLIT fixes the other combination algebraically.
//  Registering ICs for both would over-determine the system.
//
//  WHAT TO LOOK FOR
//    [classify] type=DIFFERENTIAL_RECTANGULAR  sym=1x2   -> rev70 names the shape, and
//        the rule fires on a model with no mesh, no ALE, no front and no auxiliaries.
//    [classify] type=UNDETERMINED              sym=1x2   -> pre-rev70 behaviour.
//    [classify] type=DIFFERENTIAL_ALGEBRAIC    sym=2x2   -> the symbol builder pairs
//        rows and states differently than expected, and the whole rectangular-core
//        account (handoff rev20 section 3) needs revisiting.
//
//  The solve must succeed in ALL cases: this is a well-posed index-1 DAE, and
//  classification does not gate the solve.  If |x-x*| is large the MODEL is wrong,
//  not the classifier.
//
//  USAGE
//    ./OCFE_DAE5 [--evolve] [--nomarch] [--icx2] [--weak|--trace] [--tol R] [-v] [nel_t ...]
// ===========================================================================

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <limits>

#ifndef OCFE_OCFESLV_HEADER
  #define OCFE_OCFESLV_HEADER "ocfeslv.hpp"
#endif
#include OCFE_OCFESLV_HEADER

using namespace mc;

static double const TF = 1.0;
static double const S0 = 2.0;      // s(0) = x1(0) + x2(0)
static double const D0 = 0.5;      // d(0) = x1(0) - x2(0)

static size_t NND_T  = 5;
static double RESTOL = 1.0e-9;
static bool   g_trace = false, g_weak = false, EVOLVE = false, NOMARCH = false;
static bool   IC_ON_X2 = false;   // --icx2: close x2 instead of x1 (see the note at IC_X1)
static int    g_disp  = 0;

static double const XTOL = 1.0e-8;
static double const ORDTOL = 3.0;   // rev106: minimum observed convergence order

// closed-form solution
static double s_exact( double t ){ return S0 * std::exp( -t ); }
static double d_exact( double t ){ return D0 * std::exp( -2.0 * t ); }
static double x1_exact( double t ){ return 0.5 * ( s_exact(t) + d_exact(t) ); }
static double x2_exact( double t ){ return 0.5 * ( s_exact(t) - d_exact(t) ); }

struct Result
{
  bool   setup_ok = false, square = false, conv = false, threw = false;
  double err = std::numeric_limits<double>::infinity();
  size_t nVar = 0;
};

// ---------------------------------------------------------------------------
static Result run( size_t nel_t, int display )
{
  Result R;

  FFGraph DAG;
  FFVar t  = DAG.add_var( "t" );
  FFVar x1 = DAG.add_var( "x1(t)" );
  FFVar x2 = DAG.add_var( "x2(t)" );

  FFPartial OpP;   // OpP is an OBJECT, declared per driver (OCFE_DAE0.cpp:190), not a
                   // free function -- omitting this gave "'OpP' was not declared".

  // SUM differentiates BOTH x1 and x2 in ONE row -- the whole point of this driver.
  FFVar SUM   = OpP( x1, t ) + OpP( x2, t ) + ( x1 + x2 );
  // SPLIT is algebraic: no OpP, so it is routed to vAlgEqn and holds no symbol row.
  FFVar SPLIT = ( x1 - x2 ) - D0 * exp( -2.0 * t );
  // ONE initial condition only: SPLIT fixes the other combination algebraically, so
  // registering ICs for both x1 and x2 would over-determine the system.
  //
  // --icx2 DISCRIMINATOR (handoff rev29).  MEASURED at auto-detect: the interface plan
  // mints 6/14/30 continuity claims at nel_t=4/8/16, of which rank(B) supports only half,
  // and the orphaned tau columns are ALWAYS x2(t) -- one per interior t-seam, giving
  // k = nel_t - 1 exactly.  Two hypotheses fit that equally well:
  //
  //   (a) RECTANGULAR CORE.  SUM is the only row carrying a t-derivative and it
  //       differentiates BOTH states, so it can support one claim per seam; which state
  //       gets it is decided by the Schur elimination's pivoting, not by the model.
  //   (b) UNCLOSED STATE.  x1 has its own INITIAL equation and x2 does not; the state
  //       without its own closure is the one left orphaned.
  //
  // DAE5 as written cannot separate them -- x2 is simultaneously the second column of a
  // 1x2 core AND the state without an IC.  OCFE_DAE0 does not separate them either: it
  // has four differentiating rows, four ICs, and a fifth state (lam) that is both
  // underdetermined and unclosed, the same confound with a different letter.
  //
  // Swapping the IC to x2 breaks it while changing nothing else.  SPLIT fixes
  // x1 - x2 = d0*exp(-2t), so closing x2 determines x1(0) exactly as the reverse does;
  // the model stays well-posed and the core stays 1x2.
  //
  //   orphans STAY on x2  -> (a): the core shape decides, the IC is irrelevant
  //   orphans MOVE to x1  -> (b): the closure decides, and "k = nel-1 from a rectangular
  //                          core" is the wrong account
  //
  // Read [claims] '*** N claim(s) kept EXPLICIT' and the ORPHAN lines, not just k --
  // k is expected to stay at nel_t-1 under BOTH hypotheses.
  FFVar IC_X1 = IC_ON_X2 ? ( x2 - x2_exact( 0. ) )
                         : ( x1 - x1_exact( 0. ) );

  OCFESLV oc( &DAG );
  oc.options.DISPLAY_LEVEL = display;
  oc.add_domain( t, FFDom( 0., TF, nel_t, FFDom::LGR, NND_T ) );
  oc.add_state ( x1, { t } );
  oc.add_state ( x2, { t } );

  if( EVOLVE ) oc.set_evolution_domain( t );
  else         oc.reset_evolution_domain();

  OCFESLV::EqnOptions interior( OCFESLV::EqnRole::INTERIOR, 0 );
  OCFESLV::EqnOptions initial ( OCFESLV::EqnRole::INITIAL,  0 );

  int const T_NO_LB = FFDom::ALL - FFDom::LB;
  oc.add_equation( SUM,   { t }, { T_NO_LB    }, interior );
  oc.add_equation( SPLIT, { t }, { FFDom::ALL }, interior );
  oc.add_equation( IC_X1, { t }, { FFDom::LB  }, initial  );

  oc.options.REDUCE.ORDER    = OCFESLV::Options::RED_NONE;
  oc.options.INTERFACE.TYPE  = OCFESLV::Options::IC_AUTO;
  oc.options.INTERFACE.IMPOSITION = g_weak  ? OCFESLV::Options::IC_WEAK
                             : g_trace ? OCFESLV::Options::IC_TRACE
                                       : OCFESLV::Options::IC_STRONG;
  oc.options.SOLVE.RES_TOL   = RESTOL;
  oc.options.CLASSIFY.MODE        = OCFESLV::Options::CLASS_AUTO;
  if( NOMARCH ) oc.options.SOLVE.MARCHING = false;   // rev106: suppress the COLLAPSE only

  try { R.setup_ok = oc.setup(); }
  catch( ... ) { R.threw = true; return R; }
  if( !R.setup_ok ) return R;

  size_t nEqn = 0;
  try { R.nVar = oc.n_colloc_sta(); nEqn = oc.n_colloc_eqn(); } catch( ... ) {}
  R.square = ( R.nVar && R.nVar == nEqn );

  std::vector<double> varInit, inpInit;
  if( !oc.init( varInit, inpInit, nullptr ) ) return R;
  std::vector<double> xv = varInit;

  OCFESLV::SolveReport rep = oc.solve( xv.data(), nullptr, nullptr );
  R.conv = rep.converged;

  {
    double emax = 0.;
    size_t const NS = 41;
    for( size_t k = 0; k <= NS; ++k ){
      double const tt = TF * double( k ) / double( NS );
      OCFESLV::t_Coord pt; pt[t] = tt;
      double a = 0., b = 0.;
      try { a = oc.eval_colloc<double>( x1, pt, xv.data(), nullptr, nullptr );
            b = oc.eval_colloc<double>( x2, pt, xv.data(), nullptr, nullptr ); }
      catch( ... ) { continue; }
      emax = std::max( emax, std::max( std::fabs( a - x1_exact(tt) ),
                                       std::fabs( b - x2_exact(tt) ) ) );
    }
    R.err = emax;
  }

  return R;
}

// ---------------------------------------------------------------------------
int main( int argc, char* argv[] )
{
  std::vector<size_t> nelList;

  for( int i = 1; i < argc; ++i ){
    if( !std::strcmp( argv[i], "--trace"  ) ){ g_trace = true; continue; }
    if( !std::strcmp( argv[i], "--weak"   ) ){ g_weak  = true; continue; }
    if( !std::strcmp( argv[i], "--evolve" ) ){ EVOLVE  = true; continue; }
  // rev106: --nomarch makes the THIRD configuration measurable.
  //
  // CRONOS_HOIST_EVODETECT=0 and options.SOLVE.MARCHING=false both give a monolithic
  // solve, but they are NOT the same configuration:
  //
  //   HOIST=0        the evolution domain is found LATE, at classify -- so the interface
  //                  plan is built with _evolution_dom_set still FALSE, over the full
  //                  space-time tensor grid and with no evolution direction labelled.
  //   --nomarch      the hoist still runs, so the domain IS known when the plan is built;
  //                  only the COLLAPSE is suppressed.  Claim directions, the delta-(A)
  //                  gate and k are all computed with the evolution direction known.
  //
  // The second is what a driver author writes to opt out of re-gridding while keeping
  // everything else -- it is the fix applied to OCFE_deepcopy -- so its behaviour needs to
  // be a measurement rather than an assumption.
  if( !std::strcmp( argv[i], "--nomarch" ) ){ NOMARCH = true; continue; }
    if( !std::strcmp( argv[i], "--icx2"  ) ){ IC_ON_X2 = true; continue; }
    if( !std::strcmp( argv[i], "-v"       ) ){ g_disp  = 1;    continue; }
    if( !std::strcmp( argv[i], "--verbose") ){ g_disp  = 1;    continue; }
    if( !std::strcmp( argv[i], "--tol" ) && i+1 < argc ){ RESTOL = std::atof( argv[++i] ); continue; }
    if( argv[i][0] != '-' ){
      double const v = std::atof( argv[i] );
      if( v > 0. ) nelList.push_back( (size_t)v );
      continue;
    }
  }
  if( nelList.empty() ) nelList = { 4, 8, 16 };

  std::cout
    << "================================================================\n"
    << "  DAE5 -- minimal index-1 DAE with a RECTANGULAR differential core\n"
    << "    SUM  :  dx1/dt + dx2/dt = -(x1 + x2)     ONE row, TWO derivatives\n"
    << "    SPLIT:  x1 - x2 = " << D0 << " exp(-2t)               algebraic\n"
    << "  exact: s = " << S0 << " e^-t,  d = " << D0 << " e^-2t,"
       "  x1 = (s+d)/2,  x2 = (s-d)/2\n"
    << "  No mesh, no ALE, no front, no auxiliaries.  vState = {x1,x2} (2 columns)\n"
    << "  while diff_eqn = {SUM} (1 row), so A_e is 1x2 and neither branch of the\n"
    << "  causality test applies -- rank and singularity need a SQUARE matrix.\n"
    << "  IC on " << ( IC_ON_X2 ? "x2 (--icx2: rev29 discriminator)" : "x1 (default)" ) << "\n"
    << "  mode=" << ( g_weak ? "IC_WEAK" : g_trace ? "IC_TRACE" : "IC_STRONG" )
    << "  cfg=" << ( EVOLVE ? "evolve (user-supplied)" : "auto-detect" )
    << ( NOMARCH ? "  --nomarch (SOLVE_MARCHING=false: no collapse)" : "" )
    << "  tol=" << RESTOL << "\n"
    << "  The solve MUST succeed regardless: classification does not gate it.\n"
    << "================================================================\n"
    << "  nel_t  nVar   square  conv   max|x-x*|   order   status\n";

  // ===================================================================
  // rev106: the gate is CONVERGENCE ORDER + finest-mesh accuracy.
  // ===================================================================
  // The old gate required R.err < XTOL at EVERY nel_t.  That asks a 4-element
  // discretisation for mesh-independent accuracy, which no consistent scheme delivers.
  //
  // MEASURED at rev106 (both promotions live): 6.73e-07, 4.03e-08, 2.44e-09 at
  // nel_t 4, 8, 16 -- every mesh CONVERGED, nVar flat at 10 because Path B collapses the
  // evolution grid, and the observed order is 4.06 then 4.04.  The old gate read FAIL.
  // Before Path B the same driver was NOT CONVERGED at nel_t 4 and 8 with nVar 43 and 87,
  // and also read FAIL.  A gate that cannot tell those two apart is not measuring the
  // property it exists to protect.
  //
  // What actually matters for a manufactured-solution test:
  //   (a) every mesh SOLVES -- setup square, Newton converged, error finite;
  //   (b) the finest mesh meets XTOL;
  //   (c) the error falls at the expected RATE.
  // (c) is the real check: an accidentally-accurate coarse answer or a scheme that has
  // silently dropped an order both fail it, and neither is caught by an absolute
  // threshold applied uniformly.
  //
  // ORDTOL is deliberately loose (3.0 against an observed 4.05): the point is to catch a
  // LOST order, not to certify the constant.
  bool all_ok   = true;   // (a) structural: every mesh must solve
  bool fine_ok  = false;  // (b) finest mesh meets XTOL
  double prev_h = 0., prev_e = 0., worst_ord = 1e30;
  size_t nord   = 0;

  for( size_t n : nelList ){
    Result const R = run( n, g_disp );
    bool const solved = R.setup_ok && R.square && R.conv && std::isfinite( R.err );
    all_ok  = all_ok && solved;
    fine_ok = solved && R.err < XTOL;          // overwritten each pass; ends on the last

    double ord = 0.;
    if( solved && prev_e > 0. && R.err > 0. ){
      ord = std::log( prev_e / R.err ) / std::log( double(n) / prev_h );
      worst_ord = std::min( worst_ord, ord );
      ++nord;
    }
    if( solved ){ prev_h = double(n); prev_e = R.err; }

    std::cout << "  " << std::setw(5) << n
              << std::setw(7) << R.nVar
              << std::setw(8) << ( R.square ? "yes" : "NO" )
              << std::setw(7) << ( R.conv ? "yes" : "no" )
              << std::setw(12) << std::scientific << std::setprecision(2) << R.err;
    if( ord > 0. ) std::cout << std::setw(8) << std::fixed << std::setprecision(2) << ord;
    else           std::cout << std::setw(8) << "-";
    std::cout << std::setw(9) << ( R.threw ? "THREW" : !R.setup_ok ? "SETUP"
                                           : solved ? "solved" : "[--]" )
              << std::defaultfloat << "\n";
  }

  // nord == 0 means fewer than two meshes solved, so no rate is observable.  That must
  // NOT pass vacuously -- it is precisely the pre-Path-B situation (nel_t 4 and 8 did not
  // converge), where a vacuous (c) would leave the verdict resting on the single finest
  // mesh.  Require at least one observed refinement.
  bool const ord_ok = ( nord >= 1 ) && ( worst_ord >= ORDTOL );
  bool const gate_ok= all_ok && fine_ok && ord_ok;

  std::cout << "  ----------------------------------------------------------------\n"
            << "  (a) every mesh solved            : " << ( all_ok  ? "yes" : "NO" ) << "\n"
            << "  (b) finest mesh |x-x*| < " << std::scientific << std::setprecision(0)
            << XTOL << " : " << ( fine_ok ? "yes" : "NO" ) << "\n"
            << "  (c) worst observed order >= " << std::fixed << std::setprecision(1)
            << ORDTOL << "  : " << ( ord_ok ? "yes" : "NO" );
  if( nord ) std::cout << "  (worst " << std::setprecision(2) << worst_ord
                       << " over " << nord << " refinement(s))";
  std::cout << std::defaultfloat << "\n";

  std::cout << "================================================================\n"
            << "  DAE5: " << ( gate_ok ? "PASS" : "FAIL" ) << "\n"
            << "  Read the [classify] line, not just this verdict:\n"
            << "    DIFFERENTIAL_RECTANGULAR sym=1x2 -> rev70 names the shape, on a model\n"
            << "      with none of MMPDE27's complications.  The rule generalises.\n"
            << "    UNDETERMINED sym=1x2             -> pre-rev70 behaviour.\n"
            << "    DIFFERENTIAL_ALGEBRAIC sym=2x2   -> the symbol builder pairs rows and\n"
            << "      states differently than expected; handoff rev20 section 3 needs revisiting.\n"
            << "  This shape is COMMON -- total balances on a sum, Kirchhoff current laws,\n"
            << "  any conservation law on an aggregate with a constraint fixing the split.\n"
            << "================================================================\n";
  return gate_ok ? 0 : 1;
}
