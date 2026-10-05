// OCFE_MESHONLY10.cpp
// ===========================================================================
//  MMPDE5 MESH ALONE -- the r-adaptive machinery with the physics removed.
// ===========================================================================
//
//  THE MODEL.  Four states, no physics:
//      GRAD :  xg - x_xi   = 0        (the one place x_xi is taken)
//      DEF  :  g  - M*xg   = 0        (M = 1/x*_xi, PRESCRIBED -- no u feedback)
//             or, with --gdirect (DEFAULT), g - M*x_xi = 0, which removes xg from the chain
//      GGDEF:  gg - g_xi   = 0
//      RELAX:  tau*x_t - gg = 0       (MMPDE5; x is genuinely differential)
//      X_IC :  x(xi,0) = x*(xi,0)     X_LB: x(0,t)=0     X_UB: x(1,t)=L
//  The oracle gate is |x - x*| against the closed-form sinh mesh.
//
//  WHY IT EXISTS.  MMPDE26 confounds "moving mesh" with "sharp interior front".  This driver
//  keeps the mesh machinery and deletes everything else, so a failure can be attributed.
//
//  WHAT IT ESTABLISHED (2026-09-09, sandbox + corpus):
//    - The four-link chain x -> xg -> g -> gg makes the interface plan mint ONE REDUNDANT
//      MULTIPLIER PER INTERIOR SEAM: deficient clusters = nel_xi - 1 exactly (3, 7, 15, 31 at
//      nel_xi 4, 8, 16, 32).  The cause is that g = M*xg is POINTWISE IN STATES, so g's
//      continuity is implied by xg's, and the plan claims it anyway.
//    - --gdirect removes that link and the deficiency goes to 0; the determinacy gap rises
//      from 6-9 (one tolerance-nudge from rank-deficient) to 8.4e+04.  Hence the DEFAULT.
//      --chain restores the four-link form for bisection.
//    - IC_WEAK is also clean, because it builds no tau block at all.  So the deficiency lives
//      in the tau block, not in the claims themselves.
//    - |x - x*| = 0.169 * TAU exactly, mesh-independent: the MMPDE5 relaxation lag, not a
//      discretisation error.  The absolute oracle gate is therefore only meaningful at small
//      TAU, hence TAU = 1e-6 by default.
//
//  USAGE
//    ./OCFE_MESHONLY10 [--weak|--trace] [--chain] [--evolve] [--showinit] [--seedoracle]
//                      [--tau T] [--tol R] [--nelt N] [--dump] [--verbose] [Pe] [nel_xi ...]
//    --chain       restore the four-link chain (g = M*xg); default is --gdirect
//    --showinit    report the range init() produced per state block; flags an xg range
//                  containing zero, the singular point of the mesh chain
//    --seedoracle  start x at the sinh oracle and xg at x*_xi, both closed form
//    The first bare numeric >= 10 with the mesh list still empty is Pe; every other bare
//    numeric is a nel_xi.  Pe is a MESH-CLUSTERING knob here (it enters only via the monitor).
// ===========================================================================

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>            // MESHONLY7: pos_state takes a std::map<FFVar,size_t,lt_FFVar>
#include <utility>        // std::pair for the block lookup
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

// ---------------------------------------------------------------------------
//  Model constants -- IDENTICAL to MMPDE26 so the oracle mesh is the same field.
// ---------------------------------------------------------------------------
static double const L_DOM = 1.0;
static double const C_ADV = 1.0;
static double const S0    = 0.25;
static double const TF    = 0.5;
static double const DM_OVER_DELTA = 1.0;

static size_t NND_XI = 8, NEL_T = 8, NND_T = 5;   // non-const: --small shrinks these
// TAU: the mesh lags the oracle by a fixed multiple of TAU -- MEASURED |x-x*| = 0.169*TAU
// exactly (1.69e-03 / 1.69e-04 / 1.69e-05 / 1.69e-06 / 1.70e-07 at TAU 1e-2 ... 1e-6), and
// mesh-independent at fixed TAU.  The oracle gate is absolute, so it is only meaningful at a
// TAU whose lag falls below it; hence 1e-6 by default.
static double TAU    = 1.0e-6;                    // --tau  : mesh relaxation time
static double RESTOL = 1.0e-9;                    // --tol  : SOLVE_RES_TOL override
static double SIGMA0 = 1.0;                       // --sigma0
static bool   g_trace = false;                    // --trace : IC_TRACE
static bool   g_weak  = false;                    // --weak  : IC_WEAK reference
static bool   EVOLVE  = false;                    // --evolve: set_evolution_domain(t)
static bool   g_dump  = false;                    // --dump  : gnuplot files
static char const* SOLDUMP = nullptr;             // --savesol <f>
static char const* SOLLOAD = nullptr;             // --loadsol <f>
static bool   RESIDPROBE = false;                 // --residprobe : eval/deriv at the seed
static bool   SHOWINIT   = false;                 // --showinit : what init() actually gave
static bool   SEEDORACLE = false;                 // --seedoracle : start from x*(xi,t)
// 2026-09-09: DEFAULT CHANGED to true.  MEASURED: the four-link chain x -> xg -> g -> gg
// makes the plan mint one redundant multiplier per interior seam (deficient clusters = nel_xi-1
// exactly: 3, 7, 15, 31 at nel_xi 4, 8, 16, 32), because g = M*xg is POINTWISE in states and
// g's continuity is therefore implied by xg's.  --gdirect defines g from x_xi directly, which
// removes that link and takes the deficiency to 0 and the determinacy gap from 6-9 to 8.4e+04.
// --chain restores the four-link form for bisection.
static bool   GDIRECT    = true;                 // --gdirect : g = M*x_xi, bypassing xg
static bool   XGCONSUME  = false;                 // --xgconsume : keep g=M*xg AND consume xg
static int    MAXIT      = -1;                    // --maxit N : SOLVE_MAX_ITER override
static size_t MARCHWIN   = 0;                     // --marchwin N : which marched window to seed
static int    g_disp  = 0;                        // --verbose

static double const XTOL = 1.0e-6;   // |x-x*| against the sinh oracle   // oracle gate on |x - x*|, same as MMPDE26

struct Params
{
  double L = L_DOM, c = C_ADV, s0 = S0, tf = TF, Pe = 1.0e2;
  double D()     const { return c * L / Pe; }
  double delta() const { return D() / c; }
  double dm()    const { return DM_OVER_DELTA * delta(); }
  double s( double t ) const { return s0 + c * t; }
};

// --- oracle map x*(xi,t): M1a's sinh map, copied verbatim from MMPDE26 -----
struct Oracle
{
  Params p;
  double sm( double t ) const { return p.s( t ); }
  double A1( double t ) const { return std::asinh( sm( t ) / p.dm() ); }
  double a ( double t ) const { return std::asinh( sm( t ) / p.dm() )
                                     + std::asinh( ( p.L - sm( t ) ) / p.dm() ); }
  double b ( double t ) const { return A1( t ) / a( t ); }
  double x   ( double xi, double t ) const
    { return sm( t ) + p.dm() * std::sinh( a( t ) * ( xi - b( t ) ) ); }
  double x_xi( double xi, double t ) const
    { return a( t ) * p.dm() * std::cosh( a( t ) * ( xi - b( t ) ) ); }
};

struct Result
{
  bool   setup_ok = false, square = false, conv = false, threw = false;
  double errMesh  = std::numeric_limits<double>::infinity();
  double xgmin    = std::numeric_limits<double>::infinity();
  size_t nVar = 0;
};

// ---------------------------------------------------------------------------
static Result run( Params const& p, size_t nel_xi, int display )
{
  Result R;
  Oracle O{ p };

  FFGraph DAG;
  FFVar t  = DAG.add_var( "t" );
  FFVar xi = DAG.add_var( "xi" );
  FFVar x  = DAG.add_var( "x(t,xi)" );      // SOLVED mesh state (DIFFERENTIAL in t)
  FFVar xg = DAG.add_var( "xg(t,xi)" );     // mesh gradient  xg = x_xi
  FFVar g  = DAG.add_var( "g(t,xi)" );      // equidistribution flux  g = M xg
  FFVar gg = DAG.add_var( "gg(t,xi)" );     // flux divergence  gg = g_xi
  FFPartial OpP;

  double const dm = p.dm();

  // --- prescribed monitor M = 1/x*_xi, as a DAG expression -------------------
  // A known function of (xi,t) only -- it does NOT reference the solved x.  This is
  // what keeps the model free of u<->x feedback, and it is why removing the physics
  // leaves the mesh block intact rather than under-determined.
  auto ash = []( FFVar const& y ){ return log( y + sqrt( y * y + 1.0 ) ); };
  FFVar const sm  = p.s0 + p.c * t;
  FFVar const a1  = ash( sm / dm );
  FFVar const a2  = ash( ( p.L - sm ) / dm );
  FFVar const aa  = a1 + a2;                          // a(t)
  FFVar const bb  = a1 / aa;                          // b(t)
  FFVar const yv  = aa * ( xi - bb );
  FFVar const xstar_xi = aa * dm * cosh( yv );        // x*_xi(xi,t)
  FFVar const W   = 1.0 / xstar_xi;                   // monitor M

  // --- the mesh block, verbatim from MMPDE26 --------------------------------
  FFVar GRAD  = xg - OpP( x, xi );            // xg = x_xi

  // --------------------------------------------------------------------------
  // --gdirect:  g = M * x_xi   instead of   g = M * xg
  // --------------------------------------------------------------------------
  // ANALYTIC RESULT this tests.  At steady state (x_t = 0):
  //     RELAX -> gg = 0 ;  GGDEF -> g_xi = 0  =>  g is constant ON EACH ELEMENT,
  //     and GLOBALLY constant only if g is CONTINUOUS at the seams.
  //     DEF + GRAD -> x_xi = g_e * x*_xi, so integrating over element e and using
  //     x's continuity (x IS claimed) to telescope:
  //           sum_e g_e * D_e = x(1) - x(0) = L ,   D_e = x*(xi_e+1) - x*(xi_e) > 0
  //     ONE equation, nel unknowns  ->  a solution family of dimension nel-1, whose
  //     members with g_e < 0 have xg < 0: a TANGLED mesh.  The uniform member g == 1
  //     recovers x = x* exactly.
  //
  // nel-1 is 3, 7, 15 at nel_xi = 4, 8, 16 -- EXACTLY the measured intra-domain kernel
  // count, and [null-anatomy] dir 0 is the g/gg pair at weights +/-1/sqrt(2).
  //
  // WHY THE MODES DIFFER.  The claim census reads
  //     states WITH a retained claim: x, gg      states WITHOUT any claim: xg, g
  // so under EXACT imposition g's continuity is never imposed and the family is real.
  // IC_WEAK penalises every interface edge indiscriminately (the sat-survival table
  // shows g with 30 edges), so a jump in g costs residual, g is driven to a single
  // global constant, and the tangled members are not roots of the weak system at all.
  // That is exactly the warm-start asymmetry: the two systems have different ROOT SETS.
  //
  // The header names this mechanism itself, at auxiliary_algebraic_receiver_state:
  // "A's continuity follows from the primitive P's continuity plus the defining LINK...
  //  v2: the 'follows from primitive + LINK' assumption HOLDS ONLY FOR THE PENALIZED
  //  (IC_WEAK) SYSTEM ... The exact modes remove that coupling, leaving the transverse
  //  aux genuinely null."  Neither CRONOS_NO_PIN_CLOSURE nor CRONOS_LINK_PAIR_SUPPRESS
  //  changes the census, so the suppression is unconditional on this path.
  //
  // THE TEST, driver-side and reversible.  Defining g from x directly removes one
  // unclaimed link from the chain: g then depends on a state that IS claimed, rather
  // than on the unclaimed alias xg.  xg survives (GRAD still defines it) but no longer
  // sits between x and g.
  //     tangling GOES AWAY  -> the unclaimed alias in the chain is the cause, and the
  //         remedy is to claim it (a header change, narrowly scoped).
  //     tangling PERSISTS   -> the nel-1 family analysis is incomplete, and a header
  //         flag would have been premature.
  FFVar DEF   = GDIRECT ? ( g - W * OpP( x, xi ) )   // g = M x_xi   (x is claimed)
                        : ( g - W * xg );            // g = M xg     (xg is not)
  FFVar GGDEF = gg - OpP( g, xi );            // gg = g_xi   (ALL-t defining eqn)

  // --------------------------------------------------------------------------
  // --xgconsume:  SEPARATE THE TWO THINGS --gdirect CHANGED AT ONCE
  // --------------------------------------------------------------------------
  // --gdirect worked -- both exact modes went from min(xg) = -9.65e+02 [TANGLED] to
  // 1.69e-03 with min(xg) = +8.92e-02, matching IC_WEAK -- but it changed TWO things:
  //
  //   (a) the CLAIM SET:  base "WITH: x, gg / WITHOUT: xg, g"
  //                    gdirect "WITH: x, xg, g, gg / WITHOUT: (none)"
  //   (b) the DISCRETISATION: g is built from the polynomial derivative x_xi rather
  //       than from the state xg, and xg becomes CONSUMED BY NOTHING -- GRAD still
  //       defines it, but no equation reads it.  The chain shortens from
  //       x -> xg -> g -> gg  to  x -> g -> gg with xg dangling.
  //
  // So --gdirect is a different NUMERICAL SCHEME, not merely a different interface
  // plan, and the improvement may owe as much to the shorter chain as to the claims.
  //
  // The mechanism is NOT a predicate suppressing these claims, as first supposed.  The
  // emission loop reads `if( receiver_states.empty() ) continue;` -- a claim needs a
  // RECEIVER ROW, and in the base chain xg and g have none.  --gdirect gives them one
  // as a side effect of rewiring the dependencies.
  //
  // This flag keeps the ORIGINAL scheme (g = M*xg, so xg is genuinely consumed) and
  // adds one redundant consumer of xg, changing the receiver graph WITHOUT changing
  // what is discretised.  XGID is the identity xg - x_xi = 0 again, registered on the
  // interior only so it cannot fight GRAD at the faces:
  //     claim set changes AND tangling goes  -> the CLAIMS are what matter
  //     claim set changes AND tangling stays -> the claims are not sufficient; the
  //         shortened chain in --gdirect did the work, and the header change would
  //         have been aimed at the wrong thing
  FFVar XGID = xg - OpP( x, xi );            // same identity as GRAD, extra consumer
  FFVar RELAX = TAU * OpP( x, t ) - gg;       // tau x_t - gg = 0   (MMPDE5)

  double const _a0 = O.a( 0.0 ), _b0 = O.b( 0.0 );
  FFVar X_IC = x - ( p.s0 + dm * sinh( _a0 * ( xi - _b0 ) ) );   // x(xi,0) = x*(xi,0)
  FFVar X_LB = x - 0.0;                                          // x(0,t) = 0
  FFVar X_UB = x - p.L;                                          // x(1,t) = L

  // --- environment ----------------------------------------------------------
  // The DAG is passed to the CONSTRUCTOR -- there is no set_dag().  MMPDE26 line 571.
  OCFESLV oc( &DAG );
  oc.options.DISPLAY_LEVEL = display;
  // t uses LGR (Radau), xi uses LGL.  This is NOT interchangeable: Radau in the
  // evolution direction is what MMPDE26 uses for its differential state, and x is
  // differential here too (RELAX: tau x_t = gg).  Copied exactly from MMPDE26 line 572.
  oc.add_domain( t,  FFDom( 0., p.tf, NEL_T,  FFDom::LGR, NND_T  ) );
  oc.add_domain( xi, FFDom( 0., 1.0,  nel_xi, FFDom::LGL, NND_XI ) );
  // add_state( Var, vDom ) -- the domain list is REQUIRED; only the classification
  // reference value defaults.  Copied from MMPDE26 lines 577-587.
  oc.add_state ( x,  { t, xi } );
  oc.add_state ( xg, { t, xi } );
  oc.add_state ( g,  { t, xi } );
  oc.add_state ( gg, { t, xi } );

  OCFESLV::EqnOptions int_opt( OCFESLV::EqnRole::INTERIOR, 0 );
  OCFESLV::EqnOptions ini_opt( OCFESLV::EqnRole::INITIAL,  0 );
  OCFESLV::EqnOptions bnd_opt( OCFESLV::EqnRole::BOUNDARY, 0 );

  // Interior masks are built by ARITHMETIC on int, exactly as MMPDE26 does (lines
  // 650-651).  There is no FFDom::INT enumerator.
  int const T_INT  = FFDom::ALL - FFDom::LB;              // interior + UB in t
  int const XI_INT = FFDom::ALL - FFDom::LB - FFDom::UB;  // strict interior in xi

  // --xgconsume (corrected).  The first attempt ADDED XGID alongside GRAD, which
  // duplicates the same identity and made the system non-square (square=NO) -- the test
  // was void.  As Benoit pointed out, adding xg - x_xi = 0 means another equation has to
  // go, and that equation is GRAD itself.
  //
  // So: SPLIT the registration rather than duplicate it.  GRAD covers the xi FACES and
  // XGID the xi INTERIOR -- the same identity, the same discretisation, the same row and
  // column counts -- but xg is now consumed by TWO equation objects instead of one, which
  // is the receiver-graph change without any change to what is discretised.
  //
  // This isolates the confound in --gdirect, which altered BOTH the claim set (xg and g
  // gained claims) AND the scheme (g built from the polynomial derivative x_xi, with xg
  // consumed by nothing).  A claim needs a receiver row -- the emission loop reads
  // `if( receiver_states.empty() ) continue;` -- so the question is whether a receiver
  // alone is enough:
  //     square=yes, claims gain xg/g, tangling GOES  -> the CLAIMS are what matter, and
  //         the header fix is about giving these states a receiver
  //     square=yes, claims unchanged                 -> a second consumer is not what
  //         mints a receiver edge; the mechanism is elsewhere
  //     square=yes, claims gain but tangling STAYS   -> the claims are not sufficient and
  //         --gdirect's shortened chain did the work
  if( XGCONSUME ){
    // ATTEMPT 4.  Attempt 3 used XI_FACE = FFDom::ALL - XI_INT, i.e. subtracting a DERIVED
    // mask from ALL, and setup() rejected it: "MISSPECIFIED BOUNDARY/DOMAIN xi IN EQUATION
    // Z28".  Masks compose by subtracting the NAMED positions LB/UB from ALL, not by
    // complementing another mask.  A single add_equation takes ONE mask per domain, so
    // "both faces" needs TWO registrations.
    //
    // The three registrations below tile xi exactly once -- LB, UB, strict interior -- so
    // the row count equals the single ALL registration and the system stays square, while
    // xg gains a second consuming equation object.
    oc.add_equation( GRAD, { t, xi }, { FFDom::ALL, FFDom::LB }, int_opt );
    oc.add_equation( GRAD, { t, xi }, { FFDom::ALL, FFDom::UB }, int_opt );
    oc.add_equation( XGID, { t, xi }, { FFDom::ALL, XI_INT    }, int_opt );
  }
  else
    oc.add_equation( GRAD, { t, xi }, { FFDom::ALL, FFDom::ALL }, int_opt );
  oc.add_equation( DEF,   { t, xi }, { FFDom::ALL, FFDom::ALL }, int_opt );
  oc.add_equation( GGDEF, { t, xi }, { FFDom::ALL, FFDom::ALL }, int_opt );
  oc.add_equation( RELAX, { t, xi }, { T_INT,      XI_INT     }, int_opt );
  oc.add_equation( X_IC,  { t, xi }, { FFDom::LB,  FFDom::ALL }, ini_opt );
  oc.add_equation( X_LB,  { t, xi }, { T_INT,      FFDom::LB  }, bnd_opt );
  oc.add_equation( X_UB,  { t, xi }, { T_INT,      FFDom::UB  }, bnd_opt );

  // --evolve: see handoff rev14 section 2.  auto-detection and user-supply are NOT
  // equivalent and both print the same legend line; on MMPDE26 supplying it removed
  // all 42 cross-domain kernel directions and took errNode from 7.54e-01 to 9.32e-02.
  // Whether that carries over WITHOUT the physics is one of the questions here.
  if( EVOLVE ) oc.set_evolution_domain( t );
  else         oc.reset_evolution_domain();

  oc.options.REDUCE.ORDER    = OCFESLV::Options::RED_MAIN;
  oc.options.INTERFACE.TYPE  = OCFESLV::Options::IC_AUTO;
  oc.options.INTERFACE.IMPOSITION = g_weak  ? OCFESLV::Options::IC_WEAK
                             : g_trace ? OCFESLV::Options::IC_TRACE
                                       : OCFESLV::Options::IC_STRONG;
  oc.options.INTERFACE.SAT_SIGMA0      = SIGMA0;
  oc.options.SOLVE.RES_TOL   = RESTOL;
  // --maxit N: cap the LM/Newton iterations.  Measured on IC_STRONG at nel_xi=8, seeded
  // from the weak solution: ONE step takes max|r| from 1.4708e+01 to 5.97e-07 -- a factor
  // 2.5e+07 -- and the solve then STALLS, wandering between 9e-08 and 1.8e-07 for eight
  // more iterations with alpha collapsing to 1.2e-04.  So the exact system is very nearly
  // satisfied a short distance from the weak solution, and the interesting rows are the
  // ones the first step CANNOT clear.  --maxit 1 stops there so the post-step probe can
  // name them.
  if( MAXIT >= 0 ) oc.options.SOLVE.MAX_ITER = MAXIT;

  try { R.setup_ok = oc.setup(); }
  catch( ... ) { R.threw = true; return R; }
  if( !R.setup_ok ) return R;

  // n_colloc_sta(), not n_colloc_var() -- the latter does not exist.  MMPDE26 line 847.
  size_t nEqn = 0;
  try { R.nVar = oc.n_colloc_sta(); nEqn = oc.n_colloc_eqn(); } catch( ... ) {}
  R.square = ( R.nVar && R.nVar == nEqn );

  std::vector<double> varInit, inpInit;
  if( !oc.init( varInit, inpInit, nullptr ) ) return R;
  std::vector<double> xv = varInit;

  // ---------------------------------------------------------------------------
  // WHAT DOES init() ACTUALLY GIVE?   (--showinit, --seedoracle)
  // ---------------------------------------------------------------------------
  // MEASURED on OCFE_ALGFIELD2 (2026-08-16): init() returns ALL ZEROS.  This driver has
  // never seeded anything -- it takes varInit verbatim -- so if the same holds here then
  // EVERY cold-start MESHONLY run began from x == 0 at every node: not merely a poor
  // guess but a COLLAPSED mesh, with xg == 0 identically.
  //
  // That matters because xg = 0 is the singular point of the mesh chain.  ALGFIELD2 made
  // the mechanism explicit on a model with no derivatives at all: for a closure of the
  // form p*a = w, d/dp is a, so a == 0 leaves EVERY p column free rather than
  // one-pivot-from-free, and tier-0 Newton is singular at the first iterate by
  // construction (measured: --a0 0 -> "tier-0 Newton singular ... using LM steps";
  // --a0 1 -> ordinary NT it=1).
  //
  // WHAT THIS DOES AND DOES NOT EXPLAIN.  It does NOT explain the warm-start asymmetry:
  // IC_TRACE seeded with the CONVERGED IC_WEAK field (min(xg) = +8.92e-02, a good mesh)
  // still returned to the tangled root, and IC_WEAK seeded from the tangled field
  // returned to 1.69e-03 bit-identically.  Those two systems have different ROOT SETS.
  // But it may well explain the COLD-START tangling, and until it is excluded every cold
  // MESHONLY result was taken from a degenerate starting mesh -- a defect in the
  // experiment, not in the framework.
  //
  // --seedoracle starts x at the sinh oracle x*(xi,t) and xg at x*_xi, both known in
  // closed form, leaving g and gg at init()'s value.  This is the analogue of MMPDE26's
  // own "perturbed exact state and auxiliary profiles" initialisation -- that driver does
  // NOT start from zero, which is one more reason the two are not comparable cold.
  {
    auto blockof = [&]( FFVar const& V )->std::pair<size_t,size_t>{
      std::map<FFVar,size_t,lt_FFVar> ndx0;
      auto const& vs = oc.var_state();
      auto const it = vs.find( V );
      if( it != vs.end() ) for( auto const& d : it->second ) ndx0[d] = 0;
      size_t off = 0;
      try { off = oc.pos_state( V, ndx0 ); } catch( ... ) { return { 0, 0 }; }
      return { off, oc.node_colloc( V ).size() };
    };

    if( SHOWINIT || SEEDORACLE ){
      FFVar const* SV[4] = { &x, &xg, &g, &gg };
      char const*  SN[4] = { "x ", "xg", "g ", "gg" };
      std::pair<size_t,size_t> blk[4];
      for( int i = 0; i < 4; ++i ) blk[i] = blockof( *SV[i] );

      if( SHOWINIT ){
        std::cerr << "  [showinit] init() ranges per state block:\n";
        for( int i = 0; i < 4; ++i ){
          if( !blk[i].second ){ std::cerr << "  [showinit]    " << SN[i]
                                          << " : block not located\n"; continue; }
          double lo = xv[blk[i].first], hi = lo;
          for( size_t k = 0; k < blk[i].second && blk[i].first+k < xv.size(); ++k ){
            lo = std::min( lo, xv[blk[i].first+k] );
            hi = std::max( hi, xv[blk[i].first+k] ); }
          std::cerr << "  [showinit]    " << SN[i] << " : off=" << blk[i].first
                    << " n=" << blk[i].second << " range=[" << std::scientific
                    << std::setprecision(3) << lo << "," << hi << "]"
                    << ( ( i == 1 && lo <= 0. && hi >= 0. )
                           ? "   *** xg range CONTAINS ZERO: the mesh chain is singular"
                             " at this iterate ***" : "" )
                    << std::defaultfloat << "\n";
        }
      }

      if( SEEDORACLE ){
        // node_colloc gives each node's (t,xi) coordinates in add_domain order, so the
        // oracle can be evaluated per node without assuming any index layout.
        auto nx  = oc.node_colloc( x ), nxg = oc.node_colloc( xg );
        size_t seeded = 0;
        for( size_t k = 0; k < blk[0].second && k < nx.size(); ++k ){
          if( nx[k].size() < 2 ) continue;
          double const tt = nx[k][0], qq = nx[k][1];
          xv[blk[0].first+k] = O.x( qq, tt );  ++seeded;
        }
        for( size_t k = 0; k < blk[1].second && k < nxg.size(); ++k ){
          if( nxg[k].size() < 2 ) continue;
          double const tt = nxg[k][0], qq = nxg[k][1];
          xv[blk[1].first+k] = O.x_xi( qq, tt );
        }
        std::cerr << "  [seedoracle] x <- x*(xi,t) and xg <- x*_xi over " << seeded
                  << " node(s); g and gg left at init() value.\n";
      }
    }
  }

  // ---------------------------------------------------------------------------
  // WARM START  (--savesol / --loadsol)
  // ---------------------------------------------------------------------------
  // MESHONLY1 measured (2026-08-16): IC_WEAK reaches |x-x*| = 1.69e-07 and meets the
  // oracle gate, while IC_TRACE and IC_STRONG converge -- cleanly, NT it=1 alpha=1.0,
  // max|r| ~ 1e-11 -- to a field wrong by three orders (5.57e+02 / 8.54e+01 / 2.18e+01),
  // AGREEING WITH EACH OTHER TO EVERY PRINTED DIGIT, with k=0, cond(C)=1.00,
  // cond(B)=1.41 and a SYMMETRIC verdict.  Nothing the interface diagnostics can
  // measure is deficient.  The exact-mode field also has xg < 0 in places -- a TANGLED
  // mesh -- where IC_WEAK's stays positive.
  //
  // THE TEST.  Seed the exact modes with the CONVERGED IC_WEAK state vector, which here
  // is the TRUE solution rather than merely a plausible field (unlike the MMPDE13 warm
  // start, where the target was itself only IC_WEAK-accurate).
  //
  //   stays near it and converges  -> the correct field IS a root of the exact system
  //                                   and the solver lands on a second, tangled root
  //                                   from the standard guess => ROOT SELECTION, and
  //                                   the remedy is initialisation/continuation plus
  //                                   the x_xi > 0 enforcement of roadmap C5.
  //   migrates back to 5.57e+02    -> the correct field is NOT a root of the exact
  //                                   system => FORMULATION DEFECT, and no
  //                                   initialisation strategy can help.
  //
  // Only the STATE block transfers: the collocation mesh is identical across modes,
  // while the tau block exists only under IC_TRACE/IC_STRONG and differs in size.
  // n_colloc_sta() returns the TOTAL var count INCLUDING taus, so the state block is
  // the leading n_colloc_sta() - n_colloc_trace() entries.  The file records its own
  // length and the load is REFUSED on mismatch -- a warm start that half-applies would
  // produce a plausible wrong answer, which is this arc's characteristic failure.
  size_t const n_sta = oc.n_colloc_sta() - oc.n_colloc_trace();

  if( SOLLOAD ){
    std::ifstream f( SOLLOAD );
    // A MARCHED file begins with the literal "MARCH <nblocks>"; a monolithic one with a
    // bare coefficient count.  Peeking the first token lets a mismatch be REFUSED rather
    // than silently half-applied.
    std::string tok;
    if( f >> tok && tok == "MARCH" ){
      size_t nblk = 0; f >> nblk;
      if( !EVOLVE ){
        std::cerr << "  [warmstart] REFUSED: '" << SOLLOAD << "' holds a MARCHED trajectory ("
                  << nblk << " window(s)) but this run is MONOLITHIC.  A single window is"
                  " not a monolithic solution.  Not seeding.\n";
      }
      else{
        size_t const want = ( MARCHWIN < nblk ? MARCHWIN : nblk - 1 );
        bool seeded = false;
        for( size_t b = 0; b < nblk; ++b ){
          double t0 = 0., t1 = 0.; size_t nb = 0;
          if( !( f >> t0 >> t1 >> nb ) ) break;
          if( b == want ){
            if( nb != n_sta ){
              std::cerr << "  [warmstart] REFUSED: window " << b << " holds " << nb
                        << " coefficient(s) but this model has " << n_sta << ".\n";
              break;
            }
            for( size_t i = 0; i < n_sta; ++i ) f >> xv[i];
            std::cerr << "  [warmstart] seeded from MARCHED window " << b
                      << " [" << t0 << "," << t1 << "], " << n_sta
                      << " coefficient(s); tau block left at init() value.\n"
                      << "  [warmstart] NOTE the probe below assembles THIS WINDOW's system"
                         " only -- that is the point of --marchwin, but every number is"
                         " window-local and must not be compared with a monolithic run.\n";
            seeded = true;
            break;
          }
          for( size_t i = 0; i < nb; ++i ){ double d; f >> d; }
        }
        if( !seeded )
          std::cerr << "  [warmstart] window " << want << " not found in '"
                    << SOLLOAD << "'\n";
      }
    }
    else{
    f.clear(); f.seekg( 0 );
    size_t n = 0;
    if( !f || !( f >> n ) )
      std::cerr << "  [warmstart] LOAD FAILED: cannot read '" << SOLLOAD << "'\n";
    else if( n != n_sta )
      std::cerr << "  [warmstart] REFUSED: file holds " << n << " state coefficient(s) but"
                   " this model has " << n_sta << ".  Same nel_xi/NEL_T/order/tau as the"
                   " saving run?  Not seeding.\n";
    else{
      for( size_t i = 0; i < n_sta; ++i ) f >> xv[i];
      std::cerr << "  [warmstart] seeded " << n_sta << " state coefficient(s) from '"
                << SOLLOAD << "'; tau block left at init() value.\n";
    }
    }
  }

  // ---------------------------------------------------------------------------
  // RESIDUAL PROBE  (--residprobe)   -- run BEFORE solve, on the seeded vector
  // ---------------------------------------------------------------------------
  // THE QUESTION (Benoit, 2026-08-16): use the IC_WEAK solution to analyse the
  // residual/Jacobian under IC_TRACE and find what is missing or incorrect.
  //
  // MEASURED so far: IC_WEAK reaches |x-x*| = 1.69e-07 (the true solution); the exact
  // modes converge to a TANGLED field (min(xg) = -4.77e+02) with k=0, cond(C)=1.00,
  // cond(B)=1.41, SYMMETRIC and max|r| ~ 1e-11.  The two warm starts are ASYMMETRIC:
  // IC_TRACE seeded from the correct field returns to the tangled root, while IC_WEAK
  // seeded from the TANGLED field returns to 1.69e-03 bit-identically.  So the two
  // systems have DIFFERENT ROOT SETS -- not different basins -- and only the exact ones
  // admit the tangled solution.
  //
  // This probe evaluates the residual AT the loaded solution, without solving, using the
  // PUBLIC OCFESLV::eval and OCFESLV::deriv.  If the exact formulation were a faithful
  // restatement of the weak one, the PHYSICAL rows would be at truncation there.
  //
  // *** READING GUARD -- THE TAU BLOCK IS NOT A DEFECT SIGNAL. ***
  // Under IC_TRACE/IC_STRONG, var[] carries n_colloc_trace() appended multiplier columns
  // that the weak solution does not supply; --loadsol leaves them at init()'s zero.  The
  // continuity/trace rows therefore carry whatever lambda SHOULD have been, and a large
  // residual there means only "lambda != 0 at the solution", which is EXPECTED.  It is
  // the PHYSICAL rows -- those whose Jacobian support lies entirely in the state block --
  // that must be small at a correct solution regardless of lambda.  The split below is
  // computed from the sparsity pattern, not assumed.
  auto residprobe = [&]( char const* when ){
    size_t const nEq = oc.n_colloc_eqn(), nFc = oc.n_colloc_fct();
    size_t const nAll = oc.n_colloc_sta(), nTr = oc.n_colloc_trace();
    size_t const nSta = nAll - nTr;

    std::vector<double> req( nEq, 0. ), rfc( nFc > 0 ? nFc : 1, 0. );
    if( !oc.eval( req.data(), nFc > 0 ? rfc.data() : nullptr,
                  xv.data(), nullptr, nullptr ) ){
      std::cerr << "  [residprobe] eval() FAILED\n";
    }
    else{
      // Sparsity pattern.  deriv(nnz,colnz) is a TWO-CALL protocol, exactly as PDE7 uses
      // it (OCFE_PDE7.cpp:226-228): the first call fills the per-row nonzero COUNTS, the
      // caller then allocates one array per row, and the second call fills the column
      // INDICES through an array of row pointers.  MESHONLY3 made a single call with
      // scalar arguments, which returned false -- every row was counted as physical and
      // the split never happened.
      //
      // A row is TAU-TOUCHING if any of its columns is >= nSta.  PDE7's own note confirms
      // the split point: "the physical/trace split is at trace_var_offset =
      // n_colloc_sta() - n_colloc_trace()".  Under IC_STRONG the taus are ELIMINATED so
      // n_colloc_trace()==0 and every row is physical by construction -- which is why the
      // IC_STRONG reading needs no split to be trustworthy, and the IC_TRACE one does.
      std::vector<size_t> nnz( nEq, 0 );
      std::vector<std::vector<size_t>> col( nEq );
      std::vector<size_t*> colptr( nEq, nullptr );
      std::vector<char> touches_tau( nEq, 0 );
      bool have_pattern = false;
      try {
        if( oc.deriv( nnz.data(), nullptr ) ){
          for( size_t i = 0; i < nEq; ++i ){
            col[i].assign( nnz[i], 0 );
            if( nnz[i] ) colptr[i] = col[i].data();
          }
          have_pattern = oc.deriv( nnz.data(), colptr.data() );
        }
      } catch( ... ) { have_pattern = false; }
      if( have_pattern )
        for( size_t i = 0; i < nEq; ++i )
          for( size_t k = 0; k < nnz[i]; ++k )
            if( col[i][k] >= nSta ){ touches_tau[i] = 1; break; }

      // ---- ROW EQUILIBRATION -------------------------------------------------
      // The solver scales the Jacobian by rw[i] = 1/max(1, row inf-norm) INSIDE the linear
      // solve (rev66:22099) but tests convergence on the RAW residual (rev66:22085).  A row
      // whose coefficients are O(1/tau) is therefore down-weighted for the STEP and counted
      // at full magnitude for CONVERGENCE.
      //
      // RELAX is tau*x_t - gg, so at tau=1e-6 that row is six orders out of balance against
      // gg.  This is not hypothetical: the tau sweep showed a seed FOUR orders more accurate
      // (|x-x*| 1.69e-03 -> 1.70e-07, reaching the oracle gate [OK]) producing a post-step
      // max|r| 57x LARGER (6.41e-07 -> 3.67e-05).  A correct solution with a large raw
      // residual is what an unequilibrated norm does to a badly scaled row.
      //
      // Both norms are reported below.  The RAW one is what the solver's convergence test
      // sees; the EQUILIBRATED one is what the step sees, and is the fair comparison across
      // different tau.  If the seam concentration survives equilibration it is structural;
      // if it disappears, it was a scaling artefact and the rev16 headline must be withdrawn.
      std::vector<double> rw( nEq, 1.0 );
      bool have_values = false;
      if( have_pattern ){
        size_t nnzsum = 0;
        for( size_t i = 0; i < nEq; ++i ) nnzsum += nnz[i];
        std::vector<double> grad( nnzsum ? nnzsum : 1, 0. );
        try {
          if( oc.deriv( grad.data(), nullptr, xv.data(), nullptr, nullptr ) ){
            have_values = true;
            size_t off = 0;
            for( size_t i = 0; i < nEq; ++i ){
              double rmax = 0.;
              for( size_t k = 0; k < nnz[i]; ++k )
                if( col[i][k] < nSta ) rmax = std::max( rmax, std::fabs( grad[off+k] ) );
              rw[i] = 1.0 / std::max( 1.0, rmax );   // the solver's own weight
              off += nnz[i];
            }
          }
        } catch( ... ) { have_values = false; }
      }

      double phys_max = 0., tau_max = 0., phys_max_eq = 0.;
      size_t phys_arg = 0, n_tau_rows = 0, phys_arg_eq = 0;
      for( size_t i = 0; i < nEq; ++i ){
        double const a = std::fabs( req[i] );
        if( touches_tau[i] ){ ++n_tau_rows; tau_max = std::max( tau_max, a ); }
        else{
          if( a > phys_max ){ phys_max = a; phys_arg = i; }
          double const ae = a * rw[i];
          if( ae > phys_max_eq ){ phys_max_eq = ae; phys_arg_eq = i; }
        }
      }

      std::cerr << "  [residprobe] " << when << "  mode="
                << ( g_weak ? "IC_WEAK" : g_trace ? "IC_TRACE" : "IC_STRONG" )
                << "  nEqn=" << nEq << " nSta=" << nSta << " nTau=" << nTr << "\n"
                << "  [residprobe] pattern=" << ( have_pattern ? "ok" : "UNAVAILABLE -- every"
                   " row counted as physical, so the split below is NOT trustworthy" ) << "\n";
      std::cerr << std::scientific << std::setprecision(4)
                << "  [residprobe] PHYSICAL rows (no tau column): count="
                << ( nEq - n_tau_rows ) << "  max|r|=" << phys_max
                << "  at row " << phys_arg << "\n"
                << "  [residprobe] TAU-TOUCHING rows:            count="
                << n_tau_rows << "  max|r|=" << tau_max
                << "   [EXPECTED to be large: lambda is seeded at 0]\n";

      std::cerr << "  [residprobe] EQUILIBRATED (row weight 1/max(1,|row|_inf), the solver's"
                   " own):  max|r|=" << ( have_values ? phys_max_eq : -1.0 )
                << "  at row " << phys_arg_eq
                << ( have_values ? "" : "   *** deriv VALUES unavailable: weights are 1,"
                                        " equilibrated == raw ***" ) << "\n";

      // Histogram of the physical rows -- one large row among thousands is a very
      // different finding from a broad floor, and a max alone cannot tell them apart.
      // Bucketed on the EQUILIBRATED magnitude when values are available.
      int bucket[7] = {0,0,0,0,0,0,0};
      for( size_t i = 0; i < nEq; ++i ){
        if( touches_tau[i] ) continue;
        double const a = std::fabs( req[i] ) * rw[i];
        int b = a < 1e-12 ? 0 : a < 1e-9 ? 1 : a < 1e-6 ? 2
              : a < 1e-3  ? 3 : a < 1e-1 ? 4 : a < 1e+1 ? 5 : 6;
        ++bucket[b];
      }
      char const* lbl[7] = { "<1e-12","<1e-9 ","<1e-6 ","<1e-3 ","<1e-1 ","<1e+1 ",">=1e+1" };
      std::cerr << "  [residprobe] physical-row EQUILIBRATED magnitude histogram:\n";
      for( int b = 0; b < 7; ++b )
        if( bucket[b] ) std::cerr << "  [residprobe]    " << lbl[b] << " : " << bucket[b] << "\n";

      // WHICH STATES do the bad rows touch?  The model has four states with claims on x
      // and gg only (xg and g have none), and the coverage audit flags gg as
      // "singly-covered ... AT_RISK".  If the inconsistency concentrates on the gg
      // columns, that matches; if it is spread evenly, it does not.  Column ownership is
      // read from the pattern, so this section is printed ONLY when the pattern is live.
      if( have_pattern ){
        size_t const nStates = 4;                       // x, xg, g, gg -- add order
        char const* snm[4] = { "x ", "xg", "g ", "gg" };
        size_t const per = nSta / nStates;              // contiguous per-state blocks
        std::vector<size_t> hit( nStates, 0 );
        size_t bad = 0;
        for( size_t i = 0; i < nEq; ++i ){
          if( touches_tau[i] ) continue;
          if( std::fabs( req[i] ) * rw[i] < 1e-9 ) continue;
          ++bad;
          std::vector<char> seen( nStates, 0 );
          for( size_t k = 0; k < nnz[i]; ++k ){
            size_t const c = col[i][k];
            if( c < nSta && per ){ size_t const sidx = c / per;
              if( sidx < nStates && !seen[sidx] ){ seen[sidx] = 1; ++hit[sidx]; } }
          }
        }
        std::cerr << "  [residprobe] states touched by the " << bad
                  << " physical row(s) with EQUILIBRATED |r| >= 1e-9"
                     "  [block guess: nSta/4 -- verify against the claim list]:\n";
        for( size_t sIx = 0; sIx < nStates; ++sIx )
          std::cerr << "  [residprobe]    " << snm[sIx] << " : " << hit[sIx] << "\n";

        // ---- WHERE are the bad rows?  --------------------------------------
        // MESHONLY6 factorised the column index by hand and got it wrong: it assumed
        // SHARED seam nodes (nt*nxi = 33*57 = 1881) when the true block is 2560 =
        // nel_xi*NND_XI*NEL_T*NND_T -- every element carries its OWN full node set with
        // DUPLICATES at the seams.  Assuming duplicates away is precisely the error this
        // whole arc is about.
        //
        // node_colloc(V) returns the physical coordinate of every collocation node for a
        // state, indexed the same way as that state's coefficient block, so the position
        // is READ from the framework rather than inferred.  The header records THREE
        // wrong answers from hand-rolled variants of this (rev41-43 accumulated
        // node_colloc().size(), which adds zero for scalar states; rev44-46 used
        // pos_state(st,{}), which throws for distributed ones) and says "do not replace
        // it without reading them".
        //
        // pos_state() needs EVERY domain of the state pinned -- an empty index map throws
        // for a distributed state.  Pattern copied from OCFE_PDE38.cpp:130-139.
        {
          auto blockof = [&]( FFVar const& V )->std::pair<size_t,size_t>{
            std::map<FFVar,size_t,lt_FFVar> ndx0;
            auto const& vs = oc.var_state();
            auto const it = vs.find( V );
            if( it != vs.end() ) for( auto const& d : it->second ) ndx0[d] = 0;
            size_t off = 0;
            try { off = oc.pos_state( V, ndx0 ); } catch( ... ) { return { 0, 0 }; }
            return { off, oc.node_colloc( V ).size() };
          };

          FFVar const* SV[4] = { &x, &xg, &g, &gg };
          char const*  SN[4] = { "x ", "xg", "g ", "gg" };
          std::cerr << "  [residprobe] state blocks from pos_state/node_colloc"
                       "  (MEASURED, not assumed):\n";
          bool blocks_ok = true;
          size_t boff[4] = {0,0,0,0}, blen[4] = {0,0,0,0};
          for( int sIx = 0; sIx < 4; ++sIx ){
            auto const bo = blockof( *SV[sIx] );
            boff[sIx] = bo.first; blen[sIx] = bo.second;
            std::cerr << "  [residprobe]    " << SN[sIx] << " off=" << boff[sIx]
                      << " nodes=" << blen[sIx] << "\n";
            if( !blen[sIx] ) blocks_ok = false;
          }

          if( !blocks_ok )
            std::cerr << "  [residprobe] node_colloc returned EMPTY for at least one state"
                         " -- positions suppressed.\n";
          else{
            // xi-coordinate histogram of the bad rows.  Each row is attributed to its
            // lowest-indexed state column as a representative location (a row couples
            // several nodes, so this is a proxy, not an exact position).  Seams are the
            // interior element boundaries k/nel_xi.
            std::vector<std::vector<double>> nd[4];
            for( int sIx = 0; sIx < 4; ++sIx ) nd[sIx] = oc.node_colloc( *SV[sIx] );
            size_t n_seam = 0, n_interior = 0, n_unmapped = 0;
            double const seam_tol = 1e-9;
            for( size_t i = 0; i < nEq; ++i ){
              if( touches_tau[i] ) continue;
              if( std::fabs( req[i] ) * rw[i] < 1e-9 ) continue;
              size_t cmin = nSta;
              for( size_t k = 0; k < nnz[i]; ++k )
                if( col[i][k] < nSta && col[i][k] < cmin ) cmin = col[i][k];
              if( cmin >= nSta ){ ++n_unmapped; continue; }
              int owner = -1; size_t loc = 0;
              for( int sIx = 0; sIx < 4; ++sIx )
                if( cmin >= boff[sIx] && cmin < boff[sIx] + blen[sIx] )
                  { owner = sIx; loc = cmin - boff[sIx]; break; }
              if( owner < 0 || loc >= nd[owner].size() ){ ++n_unmapped; continue; }
              // node_colloc returns vector<vector<double>>: OUTER indexed by node, INNER
              // holding that node's coordinates in add_domain order (t, then xi here).
              // MESHONLY7 wrote nd[owner].back(), which takes the last NODE rather than
              // the last COORDINATE of node loc.
              std::vector<double> const& ncoord = nd[owner][loc];
              if( ncoord.empty() ){ ++n_unmapped; continue; }
              double const qq = ncoord.back();
              double const scaled = qq * double( nel_xi );
              double const dist = std::fabs( scaled - std::round( scaled ) );
              if( dist < seam_tol ) ++n_seam; else ++n_interior;
            }
            std::cerr << "  [residprobe] bad-row xi positions:  AT ELEMENT SEAM=" << n_seam
                      << "   element INTERIOR=" << n_interior
                      << "   unmapped=" << n_unmapped << "\n";
          }
        }
      }
      std::cerr << std::defaultfloat;
    }
  };

  // The pre-solve probe is valid under --evolve ONLY when the seed came from a marched
  // window of the SAME model (--marchwin), because eval() then assembles that same window's
  // system.  Seeding from a monolithic file, or not seeding at all, leaves xv at init()
  // and the number means nothing.
  if( RESIDPROBE ) residprobe( EVOLVE ? "PRE-SOLVE (marched window)"
                                      : "PRE-SOLVE (at the loaded vector)" );

  OCFESLV::SolveReport rep = oc.solve( xv.data(), nullptr, nullptr );
  R.conv = rep.converged;

  // POST-SOLVE probe.  With --maxit 1 this reports the residual after exactly one Newton
  // step, so the rows still large are the ones the step CANNOT clear -- few enough to
  // identify individually, unlike the ~1000-4000 rows the pre-solve probe reports.
  // Guarded on the march flag: after a marched solve xv holds the LAST WINDOW only, and
  // eval() assembles the monolithic system, so the comparison would be against a stale
  // vector -- that artefact produced a spurious max|r|=1.79e+00 on the IC_WEAK control
  // before it was caught.
  if( RESIDPROBE ){
    if( EVOLVE )
      std::cerr << "  [residprobe] POST-SOLVE SKIPPED under --evolve: this solve marched, so"
                   " xv now holds the LAST window of THIS run, which is not the window the"
                   " pre-solve probe assembled.  The PRE-solve probe above IS valid when"
                   " seeded via --marchwin.  For a post-step reading, run monolithically.\n";
    else
      residprobe( "POST-SOLVE" );
  }

  // --- MARCHED SAVE (--savesol under --evolve) -------------------------------
  // MESHONLY8's --savesol wrote xv, which after a MARCHED solve holds the LAST WINDOW
  // only.  --loadsol then accepted it because n_colloc_sta() is also one window's length,
  // so the length guard PASSED on a window-local vector and the probe reported a spurious
  // max|r| = 1.7853e+00 -- the same artefact value that contaminated the IC_WEAK control
  // earlier.  A length guard that passes is not a validity guard.
  //
  // march_trajectory() returns one MarchBlock per evolution element, each with its window
  // [t0,t1] and its own collocated var[].  Writing the whole trajectory with a "MARCH"
  // tag and a block count lets --loadsol select a window explicitly, and lets a monolithic
  // load REFUSE a marched file rather than silently accept one block of it.
  if( SOLDUMP && rep.converged && !oc.march_trajectory().empty() ){
    auto const& traj = oc.march_trajectory();
    std::ofstream f( SOLDUMP );
    f << "MARCH " << traj.size() << "\n" << std::setprecision(17) << std::scientific;
    for( auto const& B : traj ){
      f << B.t0 << " " << B.t1 << " " << B.var.size() << "\n";
      for( double v : B.var ) f << v << "\n";
    }
    std::cerr << "  [warmstart] saved MARCHED trajectory: " << traj.size()
              << " window(s), " << ( traj.empty() ? 0 : traj[0].var.size() )
              << " coefficient(s) each, to '" << SOLDUMP << "'\n";
  }
  else if( SOLDUMP ){
    if( !rep.converged )
      std::cerr << "  [warmstart] NOT SAVING: solve did not converge, so this vector is"
                   " not a solution and seeding from it would test nothing.\n";
    else{
      std::ofstream f( SOLDUMP );
      f << n_sta << "\n" << std::setprecision(17) << std::scientific;
      for( size_t i = 0; i < n_sta; ++i ) f << xv[i] << "\n";
      std::cerr << "  [warmstart] saved " << n_sta << " state coefficient(s) to '"
                << SOLDUMP << "'\n";
    }
  }

  // --- oracle error |x - x*| over a sampling grid ---------------------------
  // Same construction as MMPDE26's errMesh: max over an (t,xi) lattice of the solved
  // mesh position against the sinh map.
  {
    // min(xg) is a FIRST-CLASS result.  MESHONLY1's exact modes converged to a field
    // with xg = -37.9 at a xi-seam -- a TANGLED mesh -- while every interface
    // diagnostic read clean (k=0, cond 1.0, SYMMETRIC).  x_xi > 0 is enforced NOWHERE
    // (roadmap C5, still [DESIGN]).  A negative value invalidates |x-x*| whatever it says.
    double xmax = 0., xgmin = std::numeric_limits<double>::infinity();
    size_t const NT_S = 11, PPE = 6;
    for( size_t k = 0; k <= NT_S; ++k ){
      double const tt = p.tf * double( k ) / double( NT_S );
      for( size_t e = 0; e < nel_xi; ++e )
        for( size_t sI = 0; sI <= PPE; ++sI ){
          double const qq = ( double( e ) + double( sI ) / double( PPE ) ) / double( nel_xi );
          OCFESLV::t_Coord pt; pt[t] = tt; pt[xi] = qq;
          double xn = 0., xgn = 0.;
          try { xn  = oc.eval_colloc<double>( x,  pt, xv.data(), nullptr, nullptr );
                xgn = oc.eval_colloc<double>( xg, pt, xv.data(), nullptr, nullptr ); }
          catch( ... ) { continue; }
          xmax  = std::max( xmax, std::fabs( xn - O.x( qq, tt ) ) );
          xgmin = std::min( xgmin, xgn );
        }
    }
    R.errMesh = xmax; R.xgmin = xgmin;
  }

  // --- dump -----------------------------------------------------------------
  if( g_dump ){
    std::string const mode_tag = g_weak ? "weak" : ( g_trace ? "trace" : "strong" );
    std::string const cfg_tag  = EVOLVE ? "evolve" : "mono";
    std::ostringstream tag; tag << mode_tag << "_" << cfg_tag << "_" << nel_xi;
    std::string const fs = "meshonly_sol_"  + tag.str() + ".dat";
    std::string const fg = "meshonly_plot_" + tag.str() + ".gp";

    std::ofstream os( fs );
    // Provenance IN the file: a Pe-parsing slip contaminated a results table in the
    // previous session and was only caught by this line.
    os << "# Pe=" << p.Pe << " nel=" << nel_xi << " delta=" << p.delta()
       << " NEL_T=" << NEL_T << " tau=" << TAU
       << " mode=" << mode_tag << " cfg=" << cfg_tag
       << " conv=" << ( R.conv ? "yes" : "NO" ) << " tol=" << RESTOL << "\n"
       << "# t  xi  x  x_star  xg  xg_star\n";
    size_t const NT_S = 11, PPE = 6;
    for( size_t k = 0; k <= NT_S; ++k ){
      double const tt = p.tf * double( k ) / double( NT_S );
      for( size_t e = 0; e < nel_xi; ++e ){
        for( size_t sI = 0; sI <= PPE; ++sI ){
          double const qq = ( double( e ) + double( sI ) / double( PPE ) ) / double( nel_xi );
          OCFESLV::t_Coord pt; pt[t] = tt; pt[xi] = qq;
          double xn = 0., xgn = 0.;
          try { xn  = oc.eval_colloc<double>( x,  pt, xv.data(), nullptr, nullptr );
                xgn = oc.eval_colloc<double>( xg, pt, xv.data(), nullptr, nullptr ); }
          catch( ... ) { continue; }
          os << std::fixed << std::setprecision(6) << tt << " " << qq << " "
             << std::setprecision(10) << xn << " " << O.x( qq, tt ) << " "
             << xgn << " " << O.x_xi( qq, tt ) << "\n";
        }
        os << "\n";
      }
      os << "\n";
    }

    std::ofstream og( fg );
    og << "# gnuplot -p " << fg << "\n"
       << "set term qt size 1400,500\n"
       << "set multiplot layout 1,3 title 'MESHONLY " << mode_tag << "/" << cfg_tag
       << "  Pe=" << p.Pe << " nel=" << nel_xi << " NEL_T=" << NEL_T
       << " tau=" << TAU << "  conv=" << ( R.conv ? "yes" : "NO" ) << "'\n"
       << "tmid=" << std::fixed << std::setprecision(6) << p.tf * 6.0 / 11.0 << "\n"
       << "set title 'mesh trajectories z(t)'\n"
       << "set xlabel 't'; set ylabel 'x'\n"
       << "plot '" << fs << "' using 1:3 with lines lc 'grey' notitle, "
          "" << p.s0 << "+" << p.c << "*x with lines lw 2 lc 'red' title 'front s(t)'\n"
       << "set title 'x vs xi at t=tmid'\n"
       << "set xlabel 'xi'; set ylabel 'x'\n"
       << "plot '" << fs << "' using ($1==tmid?$2:1/0):3 with linespoints pt 7 ps 0.6 "
          "title 'x (solved)', '" << fs
       << "' using ($1==tmid?$2:1/0):4 with lines lw 2 title 'x* oracle'\n"
       << "set title 'xg = x_xi vs xi at t=tmid'\n"
       << "set xlabel 'xi'; set ylabel 'xg'\n"
       << "plot '" << fs << "' using ($1==tmid?$2:1/0):5 with linespoints pt 7 ps 0.6 "
          "title 'xg (solved)', '" << fs
       << "' using ($1==tmid?$2:1/0):6 with lines lw 2 title 'xg* oracle'\n"
       << "unset multiplot\n";
    std::cerr << "    [dump] wrote " << fs << ", " << fg
              << "  (run:  gnuplot -p " << fg << ")\n";
  }

  return R;
}

// ---------------------------------------------------------------------------
int main( int argc, char* argv[] )
{
  Params p;
  std::vector<size_t> nelList;
  bool pe_set = false;

  for( int i = 1; i < argc; ++i ){
    if( !std::strcmp( argv[i], "--small" ) ){ NND_XI = 5; NEL_T = 2; NND_T = 3; continue; }
    if( !std::strcmp( argv[i], "--trace"  ) ){ g_trace = true;  continue; }
    if( !std::strcmp( argv[i], "--weak"   ) ){ g_weak  = true;  continue; }
    if( !std::strcmp( argv[i], "--evolve" ) ){ EVOLVE  = true;  continue; }
    if( !std::strcmp( argv[i], "--dump"   ) ){ g_dump  = true;  continue; }
    if( !std::strcmp( argv[i], "--savesol" ) && i+1 < argc ){ SOLDUMP = argv[++i]; continue; }
    if( !std::strcmp( argv[i], "--loadsol" ) && i+1 < argc ){ SOLLOAD = argv[++i]; continue; }
    if( !std::strcmp( argv[i], "--residprobe" ) ){ RESIDPROBE = true; continue; }
    if( !std::strcmp( argv[i], "--showinit"   ) ){ SHOWINIT   = true; continue; }
    if( !std::strcmp( argv[i], "--seedoracle" ) ){ SEEDORACLE = true; continue; }
    if( !std::strcmp( argv[i], "--gdirect"    ) ){ GDIRECT    = true; continue; }
    if( !std::strcmp( argv[i], "--chain"      ) ){ GDIRECT    = false; continue; }   // 2026-09-09: restore the four-link chain
    if( !std::strcmp( argv[i], "--xgconsume"  ) ){ XGCONSUME  = true; continue; }
    if( !std::strcmp( argv[i], "--maxit" ) && i+1 < argc ){ MAXIT = std::atoi( argv[++i] ); continue; }
    if( !std::strcmp( argv[i], "--marchwin" ) && i+1 < argc )
      { MARCHWIN = (size_t)std::atol( argv[++i] ); continue; }
    if( !std::strcmp( argv[i], "--verbose") ){ g_disp  = 1;     continue; }
    if( !std::strcmp( argv[i], "--nelt"   ) && i+1 < argc ){ NEL_T  = (size_t)std::atol( argv[++i] ); continue; }
    if( !std::strcmp( argv[i], "--tau"    ) && i+1 < argc ){ TAU    = std::atof( argv[++i] ); continue; }
    if( !std::strcmp( argv[i], "--tol"    ) && i+1 < argc ){ RESTOL = std::atof( argv[++i] ); continue; }
    if( !std::strcmp( argv[i], "--sigma0" ) && i+1 < argc ){ SIGMA0 = std::atof( argv[++i] ); continue; }
    // Bare numerics: first one >= 10 with the mesh list still empty is Pe.
    // SAME CONVENTION AS MMPDE26, and the same trap -- "--trace 16" sets Pe, not nel_xi.
    if( argv[i][0] != '-' ){
      double const v = std::atof( argv[i] );
      if( !pe_set && nelList.empty() && v >= 10.0 ){ p.Pe = v; pe_set = true; continue; }
      if( v > 0. ) nelList.push_back( (size_t)v );
      continue;
    }
  }
  if( nelList.empty() ) nelList = { 4, 8, 16, 32 };

  std::cout
    << "================================================================\n"
    << "  MESHONLY1 -- MMPDE5 mesh with the PHYSICS REMOVED\n"
    << "  states: x (differential), xg, g, gg.  No u, no q, no PDE.\n"
    << ( GDIRECT ? "  --gdirect: g = M*x_xi (from the CLAIMED state x), not g = M*xg\n"
                 : "" )
    << "  mode=" << ( g_weak ? "IC_WEAK" : g_trace ? "IC_TRACE" : "IC_STRONG" )
    << "   cfg=" << ( EVOLVE ? "evolve" : "monolithic" ) << "\n"
    << "  Pe=" << p.Pe << " delta=" << p.delta() << " tau=" << TAU
    << " sigma0=" << SIGMA0 << " tol=" << RESTOL
    << "  t-mesh " << NEL_T << "x" << NND_T << "  xi order " << NND_XI << "\n"
    << "  gate: setup square + solve converges + |x-x*| < " << XTOL << "\n"
    << "  NOTE Pe here is a MESH-CLUSTERING knob (it enters only via the monitor),\n"
    << "  not a physics knob -- there is no PDE in this model.\n"
    << "================================================================\n"
    << "  nel   nVar    square  conv   |x-x*|      min(xg)     status\n";

  bool all_ok = true;
  for( size_t n : nelList ){
    Result const R = run( p, n, g_disp );
    bool const ok = R.setup_ok && R.square && R.conv
                 && std::isfinite( R.errMesh ) && R.errMesh < XTOL
                 && R.xgmin > 0.;
    all_ok = all_ok && ok;
    std::cout << "  " << std::setw(4) << n
              << std::setw(8) << R.nVar
              << std::setw(9) << ( R.square ? "yes" : "NO" )
              << std::setw(7) << ( R.conv ? "yes" : "no" )
              << std::setw(12) << std::scientific << std::setprecision(2) << R.errMesh
              << std::setw(12) << R.xgmin
              << std::setw(11) << ( R.threw ? "THREW"
                                  : R.xgmin <= 0. ? "[TANGLED]"
                                  : ok ? "[OK]" : "[--]" )
              << std::defaultfloat << "\n";
  }

  std::cout << "================================================================\n"
            << "  MESHONLY1: " << ( all_ok ? "PASS" : "FAIL" ) << "\n"
            << "  ANSWERED 2026-09-09: the nel_xi threshold was neither \"mesh machinery\" nor\n"
            << "  \"physics coupling\" -- it is a CONDITIONING boundary that nel_xi moves the\n"
            << "  model across.  With the four-link chain (--chain) the plan mints one\n"
            << "  redundant multiplier per interior seam (deficient clusters = nel_xi-1) and\n"
            << "  the determinacy gap sits at 6-9, a tolerance-nudge from rank-deficient; the\n"
            << "  mesh then tangles.  The default --gdirect removes the implied link, the\n"
            << "  deficiency goes to 0 and the gap to 8.4e+04.  IC_WEAK is clean either way,\n"
            << "  because it builds no tau block at all.\n"
            << "================================================================\n";
  return all_ok ? 0 : 1;
}
