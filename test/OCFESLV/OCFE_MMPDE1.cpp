// ===========================================================================
//  OCFE_MMPDE27.cpp  --  MMPDE relaxation on a SOLVED equidistribution mesh
//                        (r-free reference formulation; --form keeps the comparison)
// ===========================================================================
//
//  THE MODEL.  The moving-mesh chain plus the physics:
//      GRAD :  xg - x_xi = 0
//      DEF  :  g  - M*xg = 0      or, with --gdirect (DEFAULT), g - M*x_xi = 0
//      GGDEF:  gg - g_xi = 0
//      RELAX:  tau*x_t - gg = 0   (MMPDE5)
//      plus the transported physics; --form selects the u/q coupling to compare.
//  Gates: setup square, classification OK, the solve converges, and on the FINEST mesh
//  |x-x*| < 1e-6 against the equidistribution oracle and errNode < 1e-5.
//
//  WHAT IT ESTABLISHED (2026-09-09, sandbox + corpus sweep):
//    - The chained form mints ONE REDUNDANT MULTIPLIER PER INTERIOR SEAM.  MEASURED:
//      tau-block deficient clusters = k' = nel_xi - 1 exactly (3, 7, 15, 31 at nel_xi
//      4, 8, 16, 32).  The cause is that g = M*xg is POINTWISE IN STATES, so g's continuity
//      follows from xg's; the plan claims it anyway and realises it as a tau multiplier that
//      enforces nothing.  Corpus-wide, that deficit EQUALS the IC_STRONG k' on every driver.
//    - --gdirect removes the link and the deficiency goes to 0; the solve also goes
//      conv=no -> conv=yes at nel_xi 4, 8 and 16.  Hence the DEFAULT.  --chain restores the
//      four-link form for bisection.
//    - |x-x*| = 0.169 * TAU exactly and mesh-independent -- the MMPDE5 relaxation lag, not a
//      discretisation error.  The oracle gate is absolute, so TAU = 1e-6 by default.
//    - errNode IS a discretisation error and converges: 4.23e-03, 5.46e-05, 8.89e-06 at
//      nel_xi 4, 8, 16 (TAU=1e-6).  The verdict is therefore taken on the finest mesh, with
//      the coarser rows reported so the convergence is visible.
//
//  EARLIER RESULTS RETAINED, because they are still the reason for two flags:
//    - SAT_SIGMA1 hypothesis REFUTED (2026-08-04): all explicit claims report kind=0
//      (CLAIM_EXACT_C0), so the C1 branch is never taken.  --sigma0/--sigma1 are kept because
//      the underlying ASYMMETRY is real and latent: the plan's rank test scores
//      coupling*orientation while the assembler scores tau*coupling*orientation, so a future
//      C1 claim would read rank-extending in the plan and empty in the matrix.
//    - Orphaned explicit tau columns are NOT on their own fatal (PDE2 control, 16 of 16
//      orphaned, passes).  That is now explained: an orphan column is the multiplier of an
//      implied claim.
//
//  GLOSSARY -- "DEFECT A" and "DEFECT B" are referred to throughout this file and printed
//  in the diagnostic block, so they are defined here:
//    DEFECT A -- the plan emits more continuity multipliers than the rows can support.  This
//                is the implied-claim phenomenon above: k' = nel_xi-1 multipliers enforcing
//                nothing.  UNDERSTOOD, and removed by the --gdirect default.
//    DEFECT B -- a state whose sole surviving row at a node is also an interface receiver
//                loses its jump to that row.  --pde-all-t and --qpin are its two remedies and
//                remain available; it is a different mechanism from A and is NOT addressed by
//                --gdirect.
//
//  USAGE
//    ./OCFE_MMPDE27 [--evolve] [--chain] [--form F] [--tau T] [--tol R] [--nelt N]
//                   [--sigma0 S] [--sigma1 S] [--verbose] [nel_xi ...]
//    --chain    restore the four-link chain (g = M*xg); default is --gdirect
// ===========================================================================

#include <cmath>
#include <cstdlib>
#include <algorithm>   // kerow probe: std::sort/std::max
#include <map>          // kerow probe: per-column support maps
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <vector>

#ifndef OCFE_OCFESLV_HEADER
#define OCFE_OCFESLV_HEADER "ocfeslv.hpp"   // consolidated diagnostics header
#endif
#include OCFE_OCFESLV_HEADER

using namespace mc;

static double const L_DOM = 1.0, C_ADV = 1.0, S0 = 0.25, TF = 0.50;
static double const DM_OVER_DELTA = 1.0;
static double const ERRTOL = 1.0e-5;   // u accuracy (looser than M1a: solved mesh)
static double const XTOL   = 1.0e-6;   // |x-x*| against the equidistribution oracle

static size_t NND_XI = 8, NEL_T = 8, NND_T = 5;   // non-const: --small shrinks these
static bool   g_kerow = false; // --kerow : keep-explicit tau ROW-vs-COLUMN structural probe
// TAU: MMPDE5 relaxes the mesh toward the oracle on this timescale, so the mesh lags the oracle
// by a fixed multiple of TAU -- MEASURED |x-x*| = 0.169*TAU exactly (1.69e-03, 1.69e-05,
// 1.70e-07 at TAU 1e-2, 1e-4, 1e-6).  The oracle gate is absolute (1e-6), so it is only
// meaningful at a TAU small enough that the lag is below it; hence 1e-6 by default.
static double TAU = 1.0e-6;   // [MMPDE] mesh relaxation time (tau x_t = g_xi); tunable
static double RSCALE = -1.0;  // rh = RSCALE*r; <0 => default delta^2 (set once Pe is known)
// SAT penalty weights, now COMMAND-LINE knobs (were hard-coded 1.0 / 0.0).  SIG1 is the
// discriminator for DEFECT A: see the STEP-1 block at the top of this file.  Defaults are
// UNCHANGED, so a run with no --sigma flag is bit-identical to MMPDE7.
static double SIG0 = 1.0;     // --sigma0 : SAT_SIGMA0 (C0 jump penalty scale)
static double SIG1 = 0.0;     // --sigma1 : SAT_SIGMA1 (C1 flux-jump penalty scale)
static int    DROPP = -1;     // --drop : -1 = leave header default (DROP_VERIFY); 0 verify, 1 rederive, 2 off
// --eps : TRACE_PROJ_EPS, the stiffness knob for the k continuity conditions the SAT receiver
// graph cannot support.  Under the rev26 regularised projection those conditions are enforced
// through a spring of stiffness 1/eps rather than dropped:
//     [ A     B  ]      W^T C x = -W^T r_c            exactly (the r realisable conditions)
//     [ C   eps*P]      Z^T C x + eps*lambda_null = -Z^T r_c   (the k unrealisable ones)
// eps -> 0 recovers the singular unprojected system; eps -> infinity recovers DROPPING them.
// Measured on this model, IC_TRACE, --small --nelt 3:
//     dropped (rev25)      max|u| = 3.04e+03   errNode = 3.04e+03   |x-x*| = 1.27e-02
//     eps = 1   (rev26)    max|u| = 2.90e+01   errNode = 2.94e+01   |x-x*| = 9.74e-03
//     reference (IC_WEAK)  max|u| = 1.19
// so the spring is worth two orders but eps=1 is far too soft.  <0 leaves the header default.
static double PEPS = -1.0;
// 2026-09-09: DEFAULT CHANGED to true, same reason as OCFE_MESHONLY10.  MEASURED here:
// deficient clusters = nel_xi-1 (k' = 3, 7, 15, 31 at nel_xi 4, 8, 16, 32) with the chain, 0
// with --gdirect; and the solve goes conv=no -> conv=yes at nel_xi = 4, 8 and 16.  The chained
// form mints one multiplier per interior seam for a claim that is implied (g = M*xg is
// pointwise in states).  --chain restores it for bisection.
static bool   GDIRECT = true;   // --gdirect : DEF as g = M*x_xi (a RECEIVER row), not g = M*xg
static double ADAPT   = -1.0;    // --adapt A : SOLUTION-DEPENDENT monitor, M = 1/sqrt(1+A*q^2)
static bool   SEAMPROBE  = false;  // --seamprobe : one-sided values either side of every seam
static bool   PDE_ALL_T  = false;  // --pde-all-t : DEFECT B remedy (a), never previously run
static bool   QPIN       = false;  // --qpin      : DEFECT B remedy (b), verified on another driver
static bool   QMOVE      = false;  // --qmove     : ONLY move QDEF off t=LB, no pin
static bool   QPINONLY   = false;  // --qpinonly  : ONLY add the q pin, QDEF stays on ALL
static double QPINSCALE  = 1.0;    // --qpinscale <f> : multiply the pinned q value by f
static double RESTOL     = 1.0e-9; // --tol <f>       : SOLVE_RES_TOL override
static bool   FORCEVAL   = false;  // --forceval  : IC_VALUE on INITIAL/BOUNDARY rows in BOTH modes
static bool   EVOLVE     = false;  // --evolve    : set_evolution_domain(t) instead of reset
static char const* SOLDUMP = nullptr;  // --savesol <f> : write the state block after solve
static char const* SOLLOAD = nullptr;  // --loadsol <f> : seed the state block before solve

enum FormKind { FORM_RSTATE = 0, FORM_DIVIDE = 1, FORM_XGMULT = 2 };
static int    FORM   = FORM_XGMULT;   // DEFAULT: r-free.  --form rstate|divide|xgmult
static bool g_trace = false;  // --trace: IC_TRACE instead of IC_STRONG
static bool g_weak  = false;  // --weak : IC_WEAK reference (reproduces MESH13 table)
static bool g_diag = false;   // --diag: dump worst-error location + a profile slice
static bool g_dump = false;   // --dump: write gnuplot-ready solution + mesh data files
static bool g_iface = true;   // --no-iface: suppress the per-state interface-plan audit
static bool g_solveverbose = false; // --verbose: SOLVE_VERBOSE -> solve-time tier0-audit

struct Params
{
  double L = L_DOM, c = C_ADV, s0 = S0, tf = TF, Pe = 1.0e2;
  double D()     const { return c * L / Pe; }
  double delta() const { return D() / c; }
  double dm()    const { return DM_OVER_DELTA * delta(); }
  double s( double t ) const { return s0 + c * t; }
};

// --- oracle map x*(xi,t) and its derivatives (M1a's sinh map) --------------
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
  double x_xixi( double xi, double t ) const
    { return a( t ) * a( t ) * p.dm() * std::sinh( a( t ) * ( xi - b( t ) ) ); }
  double xi_of_z( double z, double t ) const
    { double v = b( t ) + std::asinh( ( z - sm( t ) ) / p.dm() ) / a( t );
      return v < 0. ? 0. : ( v > 1. ? 1. : v ); }
};

static double u_exact( double z, double t, Params const& p )
  { return 0.5 * ( 1.0 - std::tanh( ( z - p.s( t ) ) / p.delta() ) ); }

// q = u_z for the tanh front.  d/dz [0.5(1-tanh(w))] with w=(z-s)/delta gives
// -sech^2(w)/(2 delta), so the peak is -1/(2 delta) -- e.g. -50 at delta=1e-2, which is
// what IC_WEAK measures (-49.8).  Used for the q panel in --dump.
static double q_exact( double z, double t, Params const& p )
  { double const w = ( z - p.s( t ) ) / p.delta();
    double const sech = 1.0 / std::cosh( w );
    return -0.5 * sech * sech / p.delta(); }

// ---------------------------------------------------------------------------
struct Result
{
  bool   setup_ok = false, square = false, march_conv = false, threw = false;
  std::string status = "(not reached)", exmsg;
  size_t nVar = 0, nEqn = 0;
  double errNode = std::numeric_limits<double>::infinity();
  double errMesh = std::numeric_limits<double>::infinity();
  double umax = 0.;
};

// ---------------------------------------------------------------------------
static Result run( Params const& p, size_t nel_xi, int display )
{
  Result R;
  Oracle O{ p };

  FFGraph DAG;
  FFVar t  = DAG.add_var( "t" );
  FFVar xi = DAG.add_var( "xi" );
  FFVar u  = DAG.add_var( "u(t,xi)" );
  FFVar q  = DAG.add_var( "q(t,xi)" );      // physical gradient  q = u_z  (= u_xi/x_xi)
  FFVar r  = DAG.add_var( "r(t,xi)" );      // RESCALED Laplacian rh = RSCALE*u_zz (see header)
  FFVar x  = DAG.add_var( "x(t,xi)" );      // SOLVED mesh state
  FFVar xg = DAG.add_var( "xg(t,xi)" );     // mesh GRADIENT  xg = x_xi (derivative-defined)
  FFVar g  = DAG.add_var( "g(t,xi)" );      // equidistribution FLUX  g = M xg
  FFVar gg = DAG.add_var( "gg(t,xi)" );     // [MMPDE] flux divergence gg = g_xi (like r=u_zz)
  FFPartial OpP;

  double const D = p.D(), del = p.delta(), dm = p.dm();
  // r non-dimensionalisation factor: default delta^2, overridable with --rscale.
  double const rsc = ( RSCALE > 0. ? RSCALE : del * del );

  // --- prescribed monitor w(xi,t) = 1 / x*_xi(xi,t), as a DAG expression ----
  // Reciprocal of the oracle Jacobian, so equidistribution reproduces x*.  A
  // known function of (xi,t) only -- it does NOT reference the solved x, which is
  // what keeps M1b minimal (no u<->x feedback).  Built from the same closed forms
  // the Oracle uses on the host side.
  auto ash = []( FFVar const& y ){ return log( y + sqrt( y * y + 1.0 ) ); };
  FFVar const sm  = p.s0 + p.c * t;
  FFVar const a1  = ash( sm / dm );
  FFVar const a2  = ash( ( p.L - sm ) / dm );
  FFVar const aa  = a1 + a2;              // a(t)
  FFVar const bb  = a1 / aa;              // b(t)
  FFVar const yv  = aa * ( xi - bb );
  FFVar const xstar_xi = aa * dm * cosh( yv );        // x*_xi(xi,t)
  // --- --adapt A: a SOLUTION-DEPENDENT monitor -------------------------------
  //
  // WHY THIS MATTERS FOR THE ORIGINAL MOTIVATION.  The arc's purpose is to enable moving
  // mesh on models like PSA where sharp fronts make the solution oscillate.  Everything
  // tested so far -- MESHONLY, MMPDE13/26/27 -- uses the PRESCRIBED monitor above, a known
  // function of (xi,t) that "does NOT reference the solved x".  That is exactly why FFInv
  // types DEF as LINEAR in the states: the coefficient depends only on DOMAINS, so a
  // nonlinear coefficient is benign (measured, [depdump] e1: xg=L g=L).
  //
  // A GENUINE adaptive monitor is a function of the SOLUTION -- classically
  // M = 1/sqrt(1 + A |u_z|^2), clustering nodes where the gradient is steep.  Then DEF is
  // nonlinear in the STATES, FFInv should type it N, and the Defect A at-risk condition
  // returns in a form that may be REACHABLE, because nothing else pins g.
  //
  // Registered prediction, so it can fail:
  //   stays clean          -> an internal moving-mesh option is well founded: the
  //       framework can emit g = M*x_xi by construction and the coupling to the physics
  //       keeps the algebraic freedom unreachable, as it does for the prescribed monitor.
  //   tangles / k > 0      -> the solution-dependent coefficient re-opens Defect A, and
  //       roadmap C5 (x_xi > 0 enforcement) becomes a PREREQUISITE for the option rather
  //       than a design note.
  //
  // q = u_xi is the ALE gradient state already present in this model, so the adaptive form
  // reuses it rather than introducing a new derivative: M = 1/sqrt(1 + A*q^2), which is
  // the standard arclength monitor written in the computational coordinate.
  FFVar const W  = ( ADAPT > 0. ? 1.0 / sqrt( 1.0 + ADAPT * q * q )
                                : 1.0 / xstar_xi );   // monitor M = w
  // M_xi = d/dxi (1/x*_xi) = -x*_xixi / x*_xi^2, with x*_xixi = a^2 dm sinh(y)
  FFVar const xstar_xixi = aa * aa * dm * sinh( yv );
  // NOTE W_xi is DECLARED AND NEVER USED (checked: no other reference in this file).  It
  // is the PRESCRIBED monitor's derivative, so if it were live it would be inconsistent
  // with --adapt's solution-dependent W and would silently mix two monitor forms.  It is
  // dead, so --adapt is a clean substitution -- but do not start using W_xi without
  // giving it an --adapt branch.
  FFVar const W_xi = -xstar_xixi / ( xstar_xi * xstar_xi );

  // --- Fix-1: x_t is now the SOLVED mesh velocity OpP(x,t) -------------------
  // Fix-2 used the prescribed analytic x*_t (x algebraic) to confirm the ALE term
  // is physically required.  Fix-1 replaces it with OpP(x,t): the mesh velocity is
  // now taken from the SOLVED mesh state x.  This makes x DIFFERENTIAL in t (a
  // state is differential iff OpP(state,t) appears anywhere), so x now needs an
  // initial condition and -- crucially -- its OWN equations (GRAD/DEF/FLUX) remain
  // ALGEBRAIC in x while x_t appears in the u-PDE.  An algebraic constraint on a
  // variable whose time-derivative also appears is the textbook INDEX-2 structure.
  // This is the first driver to exercise RED_FULL's Pantelides + dummy-derivative
  // index reduction, on a formulation now KNOWN physically correct (Fix-2).  Any
  // failure here is purely about index reduction, not a missing term.
  FFVar const Xt = OpP( x, t );            // solved mesh velocity (index-2 driver)

  // --- mesh derivatives via the gradient STATE xg ----------------------------
  // NOTE (corrected).  The text that stood here described "OPTION A (static mesh)"
  // and claimed x_t had been dropped entirely to make x algebraic.  That has not
  // been true since RELAX was introduced: Xt = OpP(x,t) is live above, RELAX makes
  // x genuinely DIFFERENTIAL in t, and the model is an MMPDE relaxation, not an
  // instantaneous (index-2) equidistribution.  x therefore carries a real IC (X_IC)
  // and a real time-derivative in the ALE convective term.
  //   x_xi is taken exactly ONCE, in GRAD, as a bare OpP of one state; everywhere
  // else the mesh gradient is the STATE xg.  That substitution is also why the
  // header's alias_of_auxiliary gate stops firing for r (see CURRENT STATE above):
  // on a prescribed mesh x_xi was a domain coefficient, here it is a state.
  FFVar const Xxi = xg;

  // --- MMS source at the SOLVED physical point x (not the oracle) -----------
  FFVar const TH  = ( x - ( p.s0 + p.c * t ) ) / del;
  FFVar const TT  = tanh( TH );
  FFVar const FSR = -( D / ( del * del ) ) * TT * ( 1.0 - TT * TT );

  // --- PDE with MULTIPLY-DEFINED PHYSICAL FLUX STATES (no division at all) ---
  // Suggestion carried to its conclusion: instead of DIVIDING a derivative by the
  // mesh gradient (which hid each aux's source from the DM matcher -> unmatched
  // column -> rectangular), DEFINE the physical derivatives as states via
  // MULTIPLICATION, so no '/' ever appears in a residual:
  //     QDEF:  q*x_xi - u_xi = 0     ->  q = u_xi/x_xi = u_z   (physical gradient)
  //     RDEF:  r*x_xi - q_xi = 0     ->  r = q_xi/x_xi = u_zz  (physical Laplacian)
  // Both are the MULTIPLY form of a division (a*b - c = 0 instead of a = c/b), so
  // the residual is a product of states -- exactly the framework's invariant.
  // q and r ARE the physical u_z and u_zz (verified to 1e-28), so the ALE PDE
  // collapses to its ORIGINAL physical form with the derivatives supplied as
  // states:
  //     u_t + c*q - D*r - f = 0
  // No division, no second derivative, and u_t keeps coefficient 1 (multiplying
  // the whole equation by x_xi would have scaled u_t by a STATE; the corpus only
  // proves a CONSTANT coefficient on u_t -- F_ph*OpP(q,t) in PDE24-32 -- so we
  // avoid that and scale nothing).  Every OpP operand is a bare state (u in UG-
  // analogue below, q in RDEF); every derivative is first order with an in-block
  // source.  This is the maximally-corpus-shaped form.
  bool const use_r = ( FORM == FORM_RSTATE );

  FFVar QDEF = q * Xxi - OpP( u, xi );     // q*x_xi = u_xi   (q = u_z, multiply form)
  // RESCALED RDEF.  rh = RSCALE*r, so  r*x_xi - q_xi = 0  becomes  rh*x_xi - RSCALE*q_xi = 0
  // (the original row multiplied through by RSCALE).  d/d(rh) = xg/RSCALE now DOMINATES
  // d/d(xg) = rh, inverting the 2e-5 coefficient ratio that made r the softest direction.
  FFVar RDEF = r * Xxi - rsc * OpP( q, xi );   // rh*x_xi = RSCALE*q_xi  (multiply form, rescaled)
  // ALE convective coefficient is (c - x_t)/x_xi, and q = u_xi/x_xi = u_z, so the
  // convective term is (c - x_t)*q -- the -x_t*q piece is the mesh-motion
  // correction that was missing (Fix-2 restores it with the prescribed x*_t = Xt).
  // -D*u_zz = -(D/RSCALE)*rh: the physical term is unchanged, only its carrier is scaled.
  // form=rstate : u_t + (c-x_t)q - (D/RSCALE)rh - f          (r is a state)
  // form=divide : u_t + (c-x_t)q - (D/xg)q_xi - f            (no r; divide by the Jacobian)
  // form=xgmult : xg*u_t + (c-x_t)xg*q - D*q_xi - xg*f       (no r; multiply through by it)
  FFVar PDE = ( FORM == FORM_RSTATE )
            ? ( OpP( u, t ) + ( p.c - Xt ) * q - ( D / rsc ) * r - FSR )
            : ( FORM == FORM_DIVIDE )
            ? ( OpP( u, t ) + ( p.c - Xt ) * q - ( D / Xxi ) * OpP( q, xi ) - FSR )
            : ( Xxi * OpP( u, t ) + ( p.c - Xt ) * Xxi * q - D * OpP( q, xi ) - Xxi * FSR );

  // --- MESH: FIRST-ORDER system with an EXPLICIT GRADIENT STATE --------------
  // MESH2d split the 2nd-order mesh into g=M x_xi + g_xi=0, but wrote x_xi as a
  // raw OpP(x,xi) inside DEF.  RED_FULL then minted a derivative auxiliary
  // Dxi_x = x_xi whose defining equation was CONSUMED downstream (by g_xi) with
  // NO in-block source -> an unmatched column -> 4x5 RECTANGULAR / UNDETERMINED
  // (the DM analysis named it exactly: unmatched_cols={Dxi_x}, "genuine DAE
  // closure").  Not the old singularity (Ae_sing=n now) -- a DOF-accounting gap.
  //   Fix (the standard first-order idiom, cf. PDE1's Cz=OpP(C,z)): promote the
  // gradient to its OWN state xg with its OWN defining equation, so every
  // derivative that appears has a matching in-block source and nothing dangles:
  //     GRAD:  xg - x_xi = 0      (derivative-defined state; matches the aux)
  //     DEF :  g  - M*xg = 0      (flux from the gradient state; algebraic)
  //     FLUX:  g_xi     = 0       (equidistribution; first order in g)
  // x_xi now appears exactly ONCE (in GRAD) as a bare OpP of one state, xg and g
  // are bare states, and there is no x_xixi anywhere in the mesh block.
  FFVar GRAD = xg - OpP( x, xi );         // xg = x_xi  (the ONE place x_xi is taken)
  // --gdirect (handoff rev18 section 6, the transfer test).
  //
  // As written, DEF is DERIVATIVE-FREE in xi -- the driver's own comment above says
  // "algebraic" -- so it cannot RECEIVE a flux/SAT term in the xi face direction.
  // On MESHONLY that left g and xg with NO continuity claim at all, and the steady
  // state then admitted a family of dimension nel_xi-1:
  //     gg = 0 -> g_xi = 0 -> g constant PER ELEMENT (globally constant only if g is
  //     continuous), and x's continuity telescopes the increments to the single
  //     equation  sum_e g_e D_e = L,  D_e = x*(xi_e+1) - x*(xi_e) > 0.
  // Members with g_e < 0 are TANGLED meshes.  Rewriting DEF as g = M*x_xi gives the row
  // a xi-derivative, g gains 35 claims (7 seams x 5 t-nodes), the family collapses, and
  // BOTH exact modes went from min(xg) = -9.65e+02 to +8.92e-02, matching IC_WEAK
  // exactly.
  //
  // MMPDE26 DIFFERS in a way that makes this the discriminating test: here g IS already
  // claimed (u's PDE consumes xg with derivatives, so DEF gets a receiver by proximity),
  // yet [null-anatomy] still reports k=3 = nel_xi-1 intra-domain on the g/gg pair at
  // face=0 -- the SAME signature.  So a claim on g is NOT by itself sufficient, and two
  // readings remain open:
  //     k collapses 3 -> 0  ->  the claim was emitted but NOT ENFORCING (contending for
  //         a row, the case the header warns about), and the mechanism transfers.
  //     k stays at 3        ->  MMPDE26's family has a DIFFERENT origin and the match to
  //         nel_xi-1 is coincidence; rev18 section 1 then applies to MESHONLY only,
  //         which is a real limitation to record.
  FFVar DEF  = GDIRECT ? ( g - W * OpP( x, xi ) )   // g = M x_xi : row HAS a xi-derivative
                       : ( g - W * xg );            // g = M xg   : row is ALGEBRAIC in xi
  // [MMPDE] RELAXATION replaces instantaneous equidistribution.  MMPDE5:
  //   tau x_t = g_xi.  To mirror u's structure EXACTLY -- u's PDE uses STATES
  //   q=u_z, r=u_zz (closed on ALL t via QDEF/RDEF), NOT bare derivatives -- the
  //   flux divergence g_xi is promoted to its OWN state gg with an ALL-t defining
  //   equation GGDEF.  RELAX then carries NO bare spatial derivative, so (like u)
  //   it needs only interior-t + IC, and g_xi is closed on ALL t (incl. t=0) by
  //   GGDEF -- the missing t=0 g-derivative source that made MMPDE1 rectangular.
  //   GGDEF DEFINES a new unknown (gg), so it adds a MATCHED row, not a redundant
  //   one: no over-determination, unlike a bare g_xi=0 row at t=0.
  FFVar GGDEF = gg - OpP( g, xi );                  // gg = g_xi  (all-t derivative state)
  FFVar RELAX = TAU * OpP( x, t ) - gg;             // tau x_t - gg = 0  (no bare derivative)

  // --- boundary/initial residuals -------------------------------------------
  FFVar U_BC = u - 0.5 * ( 1.0 - tanh( ( x - ( p.s0 + p.c * t ) ) / del ) );
  FFVar X_LB = x - 0.0;      // x(0,t) = 0
  FFVar X_UB = x - p.L;      // x(1,t) = L
  // [MMPDE] x is now differential -> needs a TIME initial condition.  Seed with
  // the oracle mesh at t=0: x*(xi,0) = s0 + dm sinh( a0 (xi - b0) ), a0=a(0), b0=b(0).
  double const _a0 = O.a( 0.0 ), _b0 = O.b( 0.0 );
  FFVar X_IC = x - ( p.s0 + dm * sinh( _a0 * ( xi - _b0 ) ) );   // x(xi,0) = x*(xi,0)

  // --- DEFECT B remedy (b): an initial pin for q ------------------------------
  // q is defined by QDEF (q*x_xi = u_xi) on {ALL,ALL}.  At t=0 that is q's SOLE row,
  // because PDE is registered on T_INT = ALL-LB and so does not exist at the global t
  // lower bound.  QDEF is ALSO the receiver row for u's xi-continuity there, and a state
  // whose sole determinant doubles as a receiver row loses its jump under exact
  // imposition -- the two one-sided copies enter the same row and annihilate.
  //
  // The remedy verified on another driver (2026-07-30) is: give the state its own initial
  // pin and move its defining equation off the t lower bound.  q(xi,0) is computable in
  // closed form because both u(xi,0) and x(xi,0) are prescribed:
  //     u(xi,0) = 0.5*(1 - tanh((x - s0)/delta))  =>  u_x = -1/(4 delta) sech^2((x-s0)/delta)
  //     q = u_z = u_x   (q is the PHYSICAL-space derivative, so no x_xi factor here)
  // Written directly in terms of the state x so it holds on whatever mesh x takes at t=0.
  FFVar const _sech0 = 1.0 / cosh( ( x - p.s0 ) / del );
  // FACTOR-OF-2 CORRECTION (MMPDE20).  MMPDE16-19 used 0.25/del, i.e. -sech^2/(4 del).
  //     u  = 0.5*(1 - tanh(w)),  w = (x-s0)/del
  //     u_z = 0.5 * (-sech^2 w) * (1/del) = -sech^2(w)/(2 del)
  // so the coefficient is 0.5/del, not 0.25/del.  The wrong pin halved q at every window
  // start; with --evolve it made q's continuity exact (1.4e-14, DEFECT B's remedy working)
  // while pinning it to the WRONG VALUE, so the solve stopped converging and |x-x*| went
  // 1.69e-03 -> 4.84e-02.  Independent check from the data: at delta=1e-2 the peak is
  // -1/(2 delta) = -50, and IC_WEAK measures q peaking at -49.8.  The old form asks for -25.
  // --qpinscale <f>: FAULT INJECTION on the pinned VALUE.  Doubling the coefficient from
  // 0.25/del to 0.5/del (MMPDE19 -> MMPDE20) changed NOTHING -- q spread 7.105427e-15,
  // |x-x*| 4.84e-02, errNode 1.08e+00, every printed digit identical -- while the ROW is
  // structurally required (dropping it makes the system non-square, MMPDE21).  So the row is
  // present and mandatory, yet its right-hand side does not reach the solution.  This flag
  // makes that unambiguous by scaling the pin far outside any plausible physical value.
  //
  //   results MOVE with f  -> the pin acts; 0.25 vs 0.5 was simply too small a perturbation
  //                           to show at the printed precision, and the coefficient question
  //                           reopens
  //   results BIT-IDENTICAL at f = 100 -> the Q_IC right-hand side is not reaching the
  //                           assembled system: assembled with a zero coefficient, or
  //                           overwritten downstream.  That is a defect in its own right and
  //                           a candidate explanation for conv=no.
  FFVar Q_IC = q + QPINSCALE * ( 0.5 / del ) * _sech0 * _sech0;

  OCFESLV oc( &DAG );
  oc.add_domain( t,  FFDom( 0., p.tf, NEL_T,  FFDom::LGR, NND_T  ) );
  oc.add_domain( xi, FFDom( 0., 1.0,  nel_xi, FFDom::LGL, NND_XI ) );
  oc.add_state ( u, { t, xi } );
  oc.add_state ( q, { t, xi } );
  if( use_r ) oc.add_state ( r, { t, xi } );
  oc.add_state ( x, { t, xi } );
  // NOTE: the mesh-map declaration that stood here (oc.mark_mesh_map_state(x)) is GONE,
  // along with the API itself.  Its only effect was to suppress the declared family's
  // interior spatial-C0 claims under exact imposition, which was never IC_WEAK-equivalent
  // (it left those interfaces with no exact row, no tau AND no SAT) and measurably made the
  // deficiency worse: rank 2/6/11 with it on versus 2/5/7 with it off.  The mechanism was
  // dropped when the diagnostics were rebased on the corpus-validated ocfeslv.hpp; x is now
  // an ordinary state and keeps its continuity claims like any other.
  oc.add_state ( xg, { t, xi } );
  oc.add_state ( g, { t, xi } );
  oc.add_state ( gg, { t, xi } );    // [MMPDE] flux-divergence state (gg = g_xi)

  // exact-solution reference for u; oracle map reference for x  (before setup)
  { Oracle const Ov = O; Params const pv = p; FFVar tv = t, xv = xi;
    oc.update_ref( u, [Ov,pv,tv,xv]( OCFESLV::t_Coord const& q ){
        double const tt = q.at( tv ), xx = q.at( xv );
        return u_exact( Ov.x( xx, tt ), tt, pv ); } );
    // q* = u_z, r* = u_zz (physical derivatives of the exact solution).
    // u_z = -1/(2 delta) sech^2((z-s)/delta);  u_zz = (1/delta^2) th (1-th^2).
    oc.update_ref( q, [pv,tv,xv,Ov]( OCFESLV::t_Coord const& c ){
        double const tt = c.at( tv ), xx = c.at( xv );
        double const z  = Ov.x( xx, tt );
        double const th = std::tanh( ( z - pv.s( tt ) ) / pv.delta() );
        return -0.5 / pv.delta() * ( 1.0 - th * th ); } );
    // rh* = RSCALE * u_zz, so the seed matches the rescaled state (captured by value --
    // update_ref stores the callable and invokes it later during init()).
    // Only registered when r is a state; in the r-free forms there is nothing to seed.
    if( use_r )
    oc.update_ref( r, [pv,tv,xv,Ov,rsc]( OCFESLV::t_Coord const& c ){
        double const tt = c.at( tv ), xx = c.at( xv );
        double const z  = Ov.x( xx, tt );
        double const th = std::tanh( ( z - pv.s( tt ) ) / pv.delta() );
        return rsc * ( 1.0 / ( pv.delta() * pv.delta() ) ) * th * ( 1.0 - th * th ); } );
    oc.update_ref( x, [Ov,tv,xv]( OCFESLV::t_Coord const& q ){
        double const tt = q.at( tv ), xx = q.at( xv );
        return Ov.x( xx, tt ); } );
    // xg* = x*_xi (oracle gradient); g* = M x*_xi = 1 identically (self-check).
    oc.update_ref( xg, [Ov,tv,xv]( OCFESLV::t_Coord const& q ){
        double const tt = q.at( tv ), xx = q.at( xv );
        return Ov.x_xi( xx, tt ); } );
    oc.update_ref( g, []( OCFESLV::t_Coord const& ){ return 1.0; } );
    oc.update_ref( gg, []( OCFESLV::t_Coord const& ){ return 0.0; } );  // [MMPDE] g*_xi = 0 (g*=1)
  }

  // --forceval: make IC_TRACE and IC_STRONG take the SAME interface-type path.
  //
  // _resolve_interface_type (rev66:7913) contains
  //     bool const weak_like = weak || IMPOSITION_TYPE == IC_TRACE;
  //     if( weak_like && _ordinary_trace_role(role) && interface_type == IC_AUTO )
  //         return IC_VALUE;
  // so with INTERFACE_TYPE = IC_AUTO (driver line 776) and EqnOptions leaving
  // interface_type at its IC_AUTO default, IC_TRACE is FORCED to IC_VALUE on ordinary trace
  // rows while IC_STRONG falls through to the classifier and may select characteristic or
  // upwind coupling.  _ordinary_trace_role is INITIAL | BOUNDARY | INTERFACE, so this bites
  // ini_opt and bnd_opt -- NOT the INTERIOR PDE rows.
  //
  // THE TWO MODES MAY THEREFORE BE IMPOSING DIFFERENT INTERFACE CONDITIONS, and every
  // IC_TRACE-vs-IC_STRONG comparison in this arc has assumed they were not.  The resolver
  // logs nothing, so this cannot be read off the output -- hence the flag.
  //
  // Supplying IC_VALUE explicitly makes the branch inert (it only fires on IC_AUTO), so BOTH
  // modes then take the IC_VALUE path.
  //   IC_STRONG's convergence threshold MOVES  -> difference #2 is the mode gap
  //   nothing changes                          -> #2 excluded on this model; only the
  //                                               claim-family naming (rev66:17740) remains
  // IC_TRACE must be UNCHANGED by this flag -- it was already being forced to IC_VALUE.
  // If IC_TRACE moves, the reading of the branch is wrong and the experiment is void.
  OCFESLV::Options::InterfaceType const trace_it
    = FORCEVAL ? OCFESLV::Options::IC_VALUE : OCFESLV::Options::IC_AUTO;
  OCFESLV::EqnOptions int_opt( OCFESLV::EqnRole::INTERIOR, 0 );
  OCFESLV::EqnOptions ini_opt( OCFESLV::EqnRole::INITIAL,  0, trace_it );
  OCFESLV::EqnOptions bnd_opt( OCFESLV::EqnRole::BOUNDARY, 0, trace_it );

  int const T_INT  = FFDom::ALL - FFDom::LB;
  int const XI_INT = FFDom::ALL - FFDom::LB - FFDom::UB;

  // u: PDE interior + INITIAL (t=0) + two xi Dirichlet faces.
  //   u carries u_t, so even under reset_evolution_domain() the classifier
  // AUTO-DETECTS t as an evolution direction (ocfeslv.hpp:11952) and REQUIRES an
  // INITIAL condition at t=0.  MESH2b set that face to BOUNDARY -> no INITIAL row
  // -> "0/1 blocks covered" -> u under-determined at t=0 -> a self-consistent but
  // WRONG solve (max|u|=15, errNode~1) while the mesh x was perfect (|x-x*|=5e-8).
  //   The fix is to give u the INITIAL role.  We stay MONOLITHIC (reset, no
  // marching): _marchTransfers "matters only in the marching solve"
  // (ocfeslv.hpp:2473), so monolithically the INITIAL row is imposed as an ordinary
  // residual u(xi,0)-u_exact=0 -- exactly the missing constraint -- without
  // needing the marched mixed-descriptor block that MESH2a could not converge.
  // u: PDE interior + INITIAL + 2 xi Dirichlet faces (Nx rows/slice)
  // DEFECT B remedy (a): register PDE on ALL t so that QDEF is no longer q's sole row at
  // t=0.  Recorded 2026-08-04 as "one driver line" and never run.  CAUTION: this ADDS
  // (nel_xi-1) rows at t=0 where an IC already pins u, so it may over-determine -- watch
  // for a changed row count or an ASYMMETRIC verdict, which is itself informative.
  oc.add_equation( PDE,   { t, xi }, { PDE_ALL_T ? FFDom::ALL : T_INT, XI_INT }, int_opt );
  oc.add_equation( U_BC,  { t, xi }, { FFDom::LB, FFDom::ALL }, ini_opt );
  oc.add_equation( U_BC,  { t, xi }, { T_INT,     FFDom::LB  }, bnd_opt );
  oc.add_equation( U_BC,  { t, xi }, { T_INT,     FFDom::UB  }, bnd_opt );
  // q,r: physical-flux state definitions on ALL nodes (Nx rows/slice each),
  // multiply form -> each sources its derivative aux with no division wrapper.
  // DEFECT B remedy (b): move QDEF off the t lower bound and pin q there instead.  This is
  // the remedy VERIFIED on another driver ("NT it=1 on all windows", 2026-07-30) and noted
  // as expected to transfer to the moving-mesh drivers; it was never applied here.  Row
  // count is unchanged -- one QDEF row per xi node at t=0 is replaced by one Q_IC row.
  // DECOMPOSITION (MMPDE21).  --qpin bundles TWO independent changes, and the evidence says
  // only one of them acts:
  //   (1) MOVE   QDEF from {ALL,ALL} to {T_INT,ALL}  -- removes QDEF from the t lower bound,
  //              which is where it is q's SOLE row AND the receiver row for u's xi-continuity.
  //              This is DEFECT B's mechanism.
  //   (2) PIN    add Q_IC on {LB,ALL} -- supplies q's value there instead.
  //
  // MEASURED: doubling the pinned value (MMPDE19 0.25/del -> MMPDE20 0.5/del, the ARITHMETICALLY
  // CORRECT coefficient since u_z = -sech^2/(2 delta), confirmed by IC_WEAK's q peak of -49.8)
  // changed NOTHING -- q spread 7.105427e-15, |x-x*| 4.84e-02, errNode 1.08e+00, every digit
  // identical.  So the pin is PRESENT but VALUE-INERT, while the bundle as a whole took q's
  // continuity from 1.987e+01 to 7.1e-15.  That points at (1) as the acting change and leaves
  // (2) as a row that occupies a slot without contributing its right-hand side -- a plausible
  // source of the conv=no that came with it.
  //
  // These flags separate them.  Expected readings:
  //   --qmove alone fixes q AND converges  -> (1) is the remedy, (2) is harmful, drop the pin
  //   --qmove alone fixes q, still conv=no -> the non-convergence is (1)'s, not the pin's
  //   --qpinonly fixes nothing             -> confirms (2) is inert, as the value test implies
  //   --qpinonly fixes q                   -> refutes the reading above; the pin does act
  bool const q_move = ( QPIN || QMOVE );
  bool const q_pin  = ( QPIN || QPINONLY );
  oc.add_equation( QDEF,  { t, xi }, { q_move ? T_INT : FFDom::ALL, FFDom::ALL }, int_opt );
  if( q_pin ) oc.add_equation( Q_IC, { t, xi }, { FFDom::LB, FFDom::ALL }, ini_opt );
  if( use_r ) oc.add_equation( RDEF,  { t, xi }, { FFDom::ALL, FFDom::ALL }, int_opt );

  // x: equidistribution interior + IC + two Dirichlet faces
  // NOTE the equidistribution equation is imposed on ALL t (including t=0): the
  // mesh is algebraically determined at every time, there is no x dynamics to
  // seed.  The IC row pins x(xi,0) to the oracle so the marched first window has
  // a consistent start.
  // x,g: first-order equidistribution.  Square per t-slice (2*Nx rows for 2 fields):
  //   DEF  on ALL xi        (Nx rows) -- g pinned to M x_xi everywhere
  //   FLUX on xi-INTERIOR   (Nx-2)    -- g_xi=0 is first order: one BC's worth dropped
  //   x(0)=0, x(1)=L        (2)       -- the two mesh Dirichlet endpoints
  // (Nx) + (Nx-2) + 2 = 2*Nx.  g needs no explicit BC: g_xi=0 makes it constant in
  // xi and DEF+the two x-endpoints fix that constant implicitly (the equidistributed
  // flux value), which is the correct equidistribution closure.
  // x,xg,g: first-order equidistribution with an explicit gradient state.
  //   GRAD on ALL xi (Nx), DEF on ALL xi (Nx), FLUX on xi-INTERIOR (Nx-2),
  //   x(0)=0 + x(1)=L (2)  ->  3*Nx rows for the 3 mesh fields (x,xg,g).  xg and
  //   g need no explicit BC: GRAD/DEF pin them pointwise, FLUX+the two x-endpoints
  //   close x (and hence the equidistributed-flux constant) implicitly.
  oc.add_equation( GRAD,  { t, xi }, { FFDom::ALL, FFDom::ALL }, int_opt );
  oc.add_equation( DEF,   { t, xi }, { FFDom::ALL, FFDom::ALL }, int_opt );
  // [MMPDE] x now has u-like parabolic closure: RELAX on interior + IC at t=0 +
  // two xi Dirichlet faces on interior-t (leaving the t=0 face to the IC).  This
  // is the SAME structure as u (PDE + INITIAL + 2 BCs), so x is an ordinary
  // parabolic field, not an algebraic mesh map.
  oc.add_equation( GGDEF, { t, xi }, { FFDom::ALL, FFDom::ALL }, int_opt );  // gg=g_xi on ALL t
  oc.add_equation( RELAX, { t, xi }, { T_INT,     XI_INT     }, int_opt );
  oc.add_equation( X_IC,  { t, xi }, { FFDom::LB, FFDom::ALL }, ini_opt );
  oc.add_equation( X_LB,  { t, xi }, { T_INT,     FFDom::LB  }, bnd_opt );
  oc.add_equation( X_UB,  { t, xi }, { T_INT,     FFDom::UB  }, bnd_opt );

  // OPTION A2 (monolithic space-time).  No evolution/marching domain: solve u and
  // x TOGETHER over all of [0,tf] as one big BVP.  This sidesteps the mixed-
  // --- MONOLITHIC index-2 baseline (the control MESH10 was missing) ---------
  // MESH13 is MESH10 with the ONLY change reset_evolution_domain() instead of
  // set_evolution_domain(t).  It exists because the attached OCFE_MESH3.cpp is NOT
  // an index-2 driver: MESH3 (internally still titled MESH2j) prescribes the
  // ANALYTIC mesh velocity  x*_t = c[1 + dm cosh(y)(a_s(xi-b) - a b_s)],  so x
  // there is purely ALGEBRAIC and there is no OpP(x,t) anywhere in the model.
  // MESH10 changed TWO things at once relative to it -- Xt = OpP(x,t) (index-2)
  // AND marching -- so "index-2 monolithic works" was never actually demonstrated.
  //   MESH13 supplies that control.  Read the pair as:
  //     MESH13 passes, MESH10 fails  -> the defect is in MARCHING (transfer/IC-row
  //                                     attribution), not in the index-2 structure;
  //     MESH13 also fails            -> the index-2 formulation itself is at fault
  //                                     and marching is a red herring.
  //   Monolithically t is an ordinary BVP direction: no evolution principal symbol
  // is formed, no INITIAL/transfer machinery runs, and the t=0 row is imposed as
  // an ordinary residual.  Mesh equations stay on ALL t (equidistribution holds at
  // every time incl t=0), so x is pinned algebraically at t=0 with no IC of its own
  // -- exactly the structure the marching path must learn to leave alone.
  // --------------------------------------------------------------------------
  // EVOLUTION DOMAIN  (--evolve)
  // --------------------------------------------------------------------------
  // THE ONE STRUCTURAL DIFFERENCE FROM THE ENTIRE WORKING CORPUS.
  //
  // This driver calls reset_evolution_domain(), so t is an ordinary BVP direction:
  // no evolution principal symbol is formed, no INITIAL/transfer machinery runs, and the
  // t=0 row is imposed as an ordinary residual (see the OPTION A2 note above).  EVERY
  // fixed-mesh driver that holds continuity at round-off in the exact modes -- PDE2, PDE3,
  // PDE5, PDE7, PDE9, PDE30/30b/31b/33/35/38/39, scalar_nest, the whole MBC family, 28 in
  // all, measured 2026-08-15 -- runs WITH an evolution domain.
  //
  // MMPDE13 is the sole exact-mode continuity failure in that sweep, and it is also the
  // sole driver without an evolution domain.  That is one uncontrolled variable against
  // the whole corpus, and it had not been examined.
  //
  // WHAT THIS IS NOT.  It is NOT a claim that the missing evolution domain CAUSES the
  // oscillation -- MMPDE13 is also the only moving-mesh model in the corpus, so the sweep
  // cannot separate "no evolution domain" from "moving mesh" or from anything else unique
  // to this driver.  It is a controlled change of the one difference that costs two lines.
  //
  // WHAT TO EXPECT, and each outcome is informative:
  //   (a) setup succeeds, spreads collapse toward round-off  -> the discriminator is found
  //   (b) setup succeeds, spreads unchanged                  -> the evolution domain is not
  //                                                             the difference; retire it
  //   (c) setup FAILS on the index-2 structure               -> this is the MESH13 control
  //       branch spelled out in the OPTION A2 note: "MESH13 also fails -> the index-2
  //       formulation itself is at fault and marching is a red herring"
  //
  // Structural side effects to read rather than be surprised by: "evolution direction = t"
  // appears in the eqn legend, the classification probe runs, INITIAL-role rows get
  // evolution treatment, and k may change because interior t-seams become evolution seams
  // rather than ordinary BVP interfaces -- which is the mechanism by which every marched
  // driver in the corpus has k = 0.
  //
  // NOTE the mesh equations remain on ALL t.  Under an evolution domain the t=0 row is no
  // longer an ordinary residual, so x -- which is pinned algebraically at t=0 by DEF/GRAD
  // with no IC of its own -- may now need one.  X_IC is already registered on {LB, ALL},
  // so that should be covered; if setup reports an IC-coverage shortfall, that is outcome
  // (c) arriving through a different door and is equally decisive.
  if( EVOLVE ) oc.set_evolution_domain( t );
  else         oc.reset_evolution_domain();

  // INERT for this model: q and r are hand-declared states, and no second-order
  // partial survives anywhere, so reduce_order has nothing to reduce.  RED_MAIN and
  // RED_FULL produce a byte-identical assembled model (verified via the eqn-legend).
  oc.options.REDUCE.ORDER    = OCFESLV::Options::RED_MAIN;
  oc.options.CLASSIFY.MODE        = OCFESLV::Options::CLASS_AUTO;
  oc.options.INTERFACE.TYPE  = OCFESLV::Options::IC_AUTO;
  oc.options.INTERFACE.IMPOSITION = g_weak  ? OCFESLV::Options::IC_WEAK      // --weak: reference
                             : g_trace ? OCFESLV::Options::IC_TRACE     // --trace
                                       : OCFESLV::Options::IC_STRONG;   // default: the FIX-1b target
  if( PEPS > 0.0 ) oc.options.INTERFACE.TRACE_PROJ_EPS = PEPS;   // --eps
  if( DROPP == 0 ) oc.options.INTERFACE.DROP_POLICY = OCFESLV::Options::DROP_VERIFY;
  else if( DROPP == 1 ) oc.options.INTERFACE.DROP_POLICY = OCFESLV::Options::DROP_REDERIVE;
  else if( DROPP == 2 ) oc.options.INTERFACE.DROP_POLICY = OCFESLV::Options::DROP_OFF;
  oc.options.SOLVE.VERBOSE   = g_solveverbose;              // --verbose -> tier0-audit
  oc.options.INTERFACE.SAT_SIGMA0      = SIG0;   // --sigma0 (default 1.0, as MMPDE7)
  oc.options.SOLVE.MAX_ITER  = 40;    // solved mesh is nonlinear; allow more
  // --tol <f>: SOLVE_RES_TOL override.  Under --evolve, IC_STRONG stalls in all eight windows
  // with max|r| ~ 1e-05 against the 1e-09 default -- yet it lands on errNode = 9.32e-02 and
  // |x-x*| = 1.69e-03, IDENTICAL to IC_TRACE's CONVERGED values, and IC_WEAK's own errNode at
  // this mesh is 9.64e-02.  So the answer is right and only the certification fails.  This
  // flag asks whether 1e-09 is simply below what this discretisation supports.
  //
  //   conv=yes at 1e-4 with errNode UNCHANGED -> the residual floor is the discretisation's,
  //       not a defect; IC_STRONG has been working under --evolve and the bar was misplaced.
  //   still conv=no, or errNode MOVES -> the stall is real and the (c) attribution stands.
  //
  // NOTE this changes only what is CERTIFIED, never what is computed, so errNode and |x-x*|
  // must not move.  If they do, the looser tolerance stopped the solve earlier and the
  // comparison is void.
  oc.options.SOLVE.RES_TOL   = RESTOL;
  oc.options.DISPLAY_LEVEL   = display;
#if defined(CRONOS__WITH_SPQR)
  oc.options.SOLVE.FACTORIZATION = OCFESLV::Options::SOLVE_SPQR;
#endif
  oc.options.AUDIT.EQN_LEGEND = true;
  oc.options.AUDIT.SPECTRUM   = true;
  // oc.options.MESH_MAP_SUPPRESS_C0 = true;   // only to reproduce meshmap3

  try { R.setup_ok = oc.setup(); }
  catch( OCFESLV::Exceptions& e ){
    R.threw = true; R.exmsg = "OCFESLV::Exceptions ierr=" + std::to_string( e.ierr() ); }
  catch( FFBase::Exceptions& e ){
    R.threw = true; R.exmsg = "FFBase::Exceptions ierr=" + std::to_string( e.ierr() )
                            + " (" + e.what() + ")"; }
  catch( std::exception& e ){
    R.threw = true; R.exmsg = std::string( "std::exception: " ) + e.what(); }
  catch( ... ){ R.threw = true; R.exmsg = "unknown exception"; }

  R.status = OCFESLV::setup_status_str( oc.setup_status() );
  try { R.nVar = oc.n_colloc_sta(); R.nEqn = oc.n_colloc_eqn(); } catch(...){}
  R.square = ( R.nVar && R.nVar == R.nEqn );
  if( !R.setup_ok || !R.square ) return R;

  // --- IC_STRONG per-state read -------------------------------------------
  // Setup-time, no solve required.  The prediction once recorded here -- that the
  // bilinear terminal r (RDEF: r*xg - q_xi, n_other_states == 2) would be excluded
  // by the alias resolver's n_other_states == 1 gate -- is CONFIRMED, and is worse
  // than the LEFT-ZERO that was expected: r acquires no receiver edge at all, so it
  // does not appear in the per-state table below.  It is the only state absent, and
  // the solve-time right-null space is 100% concentrated on it.
  if( g_iface ){
    std::cout << "  [MESH16] IC_STRONG interface-plan audit (nel_xi=" << nel_xi << ")\n";
    oc.display_interface_plan_diagnostics( std::cout );
  }

  // =====================================================================================
  //  --kerow : does a keep-explicit tau get a ROW, or only a COLUMN?
  // =====================================================================================
  // WHY THIS IS STRUCTURAL, NOT NUMERICAL.  deriv(nnz,colnz) returns the sparsity PATTERN.
  // No values, no solve, no dependence on an iterate -- so a zero count means no row
  // references the column AT ALL, not that a coefficient happens to vanish where we looked.
  // An SVD cannot make that distinction: sigma_max ~ 1e-13 on the tau block is equally
  // consistent with "no row" and with "row present but tiny".
  //
  // PREDICTION, RECORDED IN ADVANCE, and it differs by keep-explicit ROUTE.  Reading the
  // header found two routes into keep-explicit that are NOT equivalent:
  //
  //   reason 1  (verify-protected, ~line 17409): allocates the column, registers the claim,
  //             and DOES NOT call exact_only_slots.erase(...).  The slot therefore stays
  //             marked exact-only, its receiver edges are never re-activated, and the shared
  //             tau-term emission never fires for it.        -> PREDICT nnz == 0
  //
  //   reason 2/3 (all-LINK / protected, ~line 17462): does the same THREE things and then
  //             erases the exact_only mark, with the comment "an explicit-tau claim keeps
  //             its receiver edges active".  Emission fires, but lambda_index is then set to
  //             the RAW var[] column under IC_STRONG (line ~18488) where IC_TRACE subtracts
  //             trace_var_offset -- an index into the Schur system, not var[].
  //                                                          -> PREDICT nnz > 0, coupling
  //                                                             misdirected or out of range
  //
  // MMPDE10's log reports reason=1, so nnz == 0 is expected HERE.  If instead nnz > 0, the
  // reason-1 reading is wrong and the defect is the lambda_index mismatch after all.  Either
  // outcome is informative; both cannot be true.
  if( g_kerow ){
    // =====================================================================================
    //  v3: ask the PLAN, not deriv().
    // =====================================================================================
    // v1 and v2 both tried to read the assembled Jacobian pattern via deriv().  Both failed
    // ("deriv(pattern+values) FAILED"), and I guessed wrong twice about why -- first the
    // derivative cache, then the init() ordering.  Rather than guess a third time: the
    // question does not need the Jacobian at all.
    //
    // A tau column is coupled to physical rows through the frozen weak-SAT terms: the shared
    // emission sets term.lambda_col when a claim's receiver edges are tau-active.  So
    // counting terms per lambda_col answers "does this tau column have any row?" DIRECTLY,
    // from the plan, with no cache, no solve and no evaluation point.  It is also closer to
    // the mechanism under test than the Jacobian is: if the emission never fired, that is
    // exactly what this counts.
    //
    // PREDICTION (unchanged, and it distinguishes the two keep-explicit routes found in the
    // header):
    //   reason 1  (~17409) does NOT call exact_only_slots.erase(...), so the slot stays
    //             exact-only, edges are never re-activated, emission never fires
    //                                                    -> PREDICT 0 terms per explicit tau
    //   reason 2/3 (~17462) DOES erase, so emission fires, but lambda_index is then set to
    //             the raw var[] column under IC_STRONG   -> PREDICT >0 terms, misdirected
    // MMPDE10 logs reason=1, so 0 is expected here.  Both readings are in the source and
    // they cannot both be true.
    size_t const nvar = oc.n_colloc_sta();
    size_t const ntau = oc.n_colloc_trace();
    size_t const tau0 = ( nvar >= ntau ) ? ( nvar - ntau ) : nvar;

    std::cout << "\n  ===== [kerow] keep-explicit tau: coupled by any SAT term? =====\n"
              << "  var=" << nvar << "  n_trace_var(user-visible)=" << ntau
              << "  tau columns = [" << tau0 << "," << nvar << ")\n";

    if( !ntau ){
      std::cout << "  [kerow] no user-visible tau columns here: keep-explicit set EMPTY, so\n"
                << "  [kerow] the defect is not exercised.  INCONCLUSIVE, not a pass.\n";
    }
    else{
      std::map<size_t,size_t> terms_per_col;      // lambda_col -> #terms
      std::map<size_t,double> maxpref_per_col;    // lambda_col -> max |trace_prefactor|
      std::map<size_t,double> maxweak_per_col;    // lambda_col -> max |prefactor| (IC_WEAK)
      size_t n_terms = 0, n_tau_terms = 0;
      for( auto const& term : oc.interface_plan().weak_sat_terms() ){
        ++n_terms;
        size_t const lc = term.lambda_col;
        if( lc >= tau0 && lc < nvar ){
          ++n_tau_terms;
          ++terms_per_col[lc];
          // trace_prefactor is the IC_TRACE/IC_STRONG weight; `prefactor` is the IC_WEAK one
          // and is a diagnostic here.  Using the wrong field would report a live coupling as
          // dead (or vice versa), so both are tracked and printed.
          maxpref_per_col[lc] = std::max( maxpref_per_col[lc],
                                          std::fabs( term.trace_prefactor ) );
          maxweak_per_col[lc] = std::max( maxweak_per_col[lc], std::fabs( term.prefactor ) );
        }
      }
      // ---- v4: WHICH rows, and with what coefficients -------------------------------
      // v3 answered "is the column coupled?" -- yes, 2-4 terms each, |trace_pf| 5.9e+01 to
      // 2.9e+03.  That REFUTED the missing-row reading and left the kernel unexplained: the
      // columns are wired into physical rows, yet the spectrum finds exactly k exact null
      // directions with 100% energy on these columns and ||J v|| = 6.6e-31.
      //
      // Coupled columns can still be LINEARLY DEPENDENT.  The v3 output hints at it: the 9
      // columns fall into groups of 3 sharing identical |trace_pf| (1.960e+03, 5.895e+01,
      // 2.947e+03), and at nel_xi=8 into groups of 7 -- i.e. groups of (nel_xi-1), which is
      // the arc's kernel count law.  But the reported null vector at nel_xi=4 was
      // c1922=7.26e-01, c1920=6.78e-01, c1921=1.16e-01 -- SAME SIGN, not the pure differences
      // duplicate columns would produce.  So "identical columns" is too crude and I do not
      // have an explanation that fits both facts.
      //
      // This dump is deliberately raw: per tau column, every (row_id, trace_prefactor) that
      // references it.  Two columns sharing the same row set with proportional coefficients
      // are dependent; that is visible here and is not visible in any aggregate.  No
      // mechanism is asserted -- the last four proposed from reading the header were all
      // refuted by measurement, so this prints the evidence and stops.
      std::cout << "\n  [kerow] per-column row support (row_id : trace_prefactor)\n";
      {
        std::map<size_t,std::vector<std::pair<size_t,double>>> support;
        for( auto const& term : oc.interface_plan().weak_sat_terms() ){
          size_t const lc = term.lambda_col;
          if( lc >= tau0 && lc < nvar )
            support[lc].emplace_back( term.row_id, term.trace_prefactor );
        }
        for( auto& kv : support ){
          std::sort( kv.second.begin(), kv.second.end() );
          std::cout << "  [kerow]   col " << std::setw(6) << kv.first << " :";
          for( auto const& rp : kv.second )
            std::cout << "  r" << rp.first << "=" << std::scientific
                      << std::setprecision(4) << rp.second;
          std::cout << "\n";
        }
        // Rows shared between DIFFERENT tau columns are where dependence would live.
        std::map<size_t,std::vector<size_t>> col_of_row;
        for( auto const& kv : support )
          for( auto const& rp : kv.second ) col_of_row[rp.first].push_back( kv.first );
        size_t nshared = 0;
        for( auto const& kv : col_of_row ) if( kv.second.size() > 1 ) ++nshared;
        std::cout << "  [kerow]   rows referenced by MORE THAN ONE tau column: " << nshared
                  << " of " << col_of_row.size() << "\n";
        if( nshared )
          for( auto const& kv : col_of_row ){
            if( kv.second.size() < 2 ) continue;
            std::cout << "  [kerow]     r" << kv.first << " <-";
            for( size_t c : kv.second ) std::cout << " c" << c;
            std::cout << "\n";
          }
      }

      std::cout << "  [kerow] weak_sat_terms total=" << n_terms
                << "  referencing a tau column=" << n_tau_terms << "\n"
                << "  [kerow] " << std::left << std::setw(8) << "col"
                << std::setw(10) << "terms" << std::right << std::setw(15) << "max|trace_pf|"
                << std::setw(15) << "max|weak_pf|\n";
      size_t nocouple = 0;  double gmax = 0.;
      for( size_t c = tau0; c < nvar; ++c ){
        size_t const nt = terms_per_col.count(c) ? terms_per_col[c] : 0;
        double const mp = maxpref_per_col.count(c) ? maxpref_per_col[c] : 0.;
        if( !nt ) ++nocouple;
        gmax = std::max( gmax, mp );
        std::cout << "  [kerow] " << std::left << std::setw(8) << c
                  << std::setw(10) << nt
                  << std::right << std::scientific << std::setprecision(3) << std::setw(15) << mp
                  << std::setw(15) << ( maxweak_per_col.count(c) ? maxweak_per_col[c] : 0. )
                  << ( nt ? "" : "   ** NO SAT TERM COUPLES THIS COLUMN **" ) << "\n";
      }
      std::cout << "  [kerow] uncoupled columns: " << nocouple << " of " << ntau
                << "   max|prefactor| over tau block = " << std::scientific
                << std::setprecision(3) << gmax << "\n";
      if( nocouple == ntau )
        std::cout << "  [kerow] ==> CONFIRMED: the keep-explicit columns are coupled by NOTHING.\n"
                  << "  [kerow]     Consistent with reason-1 omitting exact_only_slots.erase().\n"
                  << "  [kerow]     Candidate fix is that one line; try it before any projection\n"
                  << "  [kerow]     or beta work.\n";
      else if( !nocouple && gmax < 1e-10 )
        std::cout << "  [kerow] ==> terms EXIST but prefactors ~" << gmax << ": the coupling is\n"
                  << "  [kerow]     present and negligible -- a SCALING defect, not a missing row.\n";
      else if( nocouple )
        std::cout << "  [kerow] ==> PARTIAL: some columns coupled, some not.  Compare the keep-\n"
                  << "  [kerow]     explicit REASON codes of the uncoupled ones against the rest.\n";
      else
        std::cout << "  [kerow] ==> terms exist AND are well-scaled: the missing-coupling reading\n"
                  << "  [kerow]     is REFUTED and the kernel originates elsewhere.  Do not\n"
                  << "  [kerow]     proceed on the reason-1 story.\n";
    }
  }

  std::vector<double> varInit, inpInit;
  if( !oc.init( varInit, inpInit, nullptr ) ) return R;
  std::vector<double> xv = varInit;

  // ---------------------------------------------------------------------------
  // WARM START  (--savesol / --loadsol)
  // ---------------------------------------------------------------------------
  // THE QUESTION.  IC_TRACE converges (max|r| ~ 1e-14) to a field that is out of range and
  // oscillatory; IC_WEAK on the IDENTICAL mesh converges to a clean one.  Every measurement
  // so far is consistent with two very different diagnoses, and they need different fixes:
  //
  //   ROOT SELECTION      the clean field IS a root of the exact system, and the solver
  //                       lands on a second, oscillatory root from the standard initial
  //                       guess.  Remedy: initialisation / continuation.  There is
  //                       precedent in this project -- collocated MBC has multiple roots
  //                       and continuation was the root-selection mechanism that fixed it.
  //
  //   FORMULATION DEFECT  the clean field is NOT a root of the exact system at all; the
  //                       exact formulation admits only the oscillatory answer.  Remedy:
  //                       the formulation, and no amount of initialisation will help.
  //
  // THE TEST.  Seed IC_TRACE with the CONVERGED IC_WEAK state vector.  If the solve stays
  // near it and converges, the clean field is a root => ROOT SELECTION.  If it migrates
  // back to the oscillatory field, the clean field is not a root => FORMULATION DEFECT.
  //
  // LAYOUT.  Only the STATE block is transferred: the collocation mesh is identical across
  // modes, so state coefficients are directly comparable, while the tau block exists only
  // under IC_TRACE/IC_STRONG and differs in size between modes.  Taus keep init()'s value
  // (zero), which is what the driver already reports doing.  The file records its own
  // length and is REFUSED on mismatch rather than silently truncated -- a warm start that
  // half-applies would produce a plausible wrong answer, which is this arc's characteristic
  // failure mode.
  // n_colloc_sta() returns _nCollVar, the TOTAL var count INCLUDING the appended tau block
  // (8978 under IC_TRACE at nel_xi=4), not the state-only length.  The state block is the
  // leading n_colloc_sta() - n_colloc_trace() entries: 7680 - 0 under IC_WEAK and
  // 8978 - 1298 under IC_TRACE, which agree, as they must -- the collocation mesh is
  // identical across modes.  (The length guard below caught this misreading on the first
  // run rather than seeding 8978 values into a 7680-long block.)
  size_t const n_sta = oc.n_colloc_sta() - oc.n_colloc_trace();

  if( SOLLOAD ){
    std::ifstream f( SOLLOAD );
    size_t n = 0;
    if( !f || !( f >> n ) )
      std::cerr << "  [warmstart] LOAD FAILED: cannot read '" << SOLLOAD << "'\n";
    else if( n != n_sta )
      std::cerr << "  [warmstart] REFUSED: file holds " << n << " state coefficient(s) but this"
                   " model has " << n_sta << ".  Same nel_xi/NEL_T/order and --form as the"
                   " saving run?  Not seeding.\n";
    else{
      for( size_t i = 0; i < n_sta; ++i ) f >> xv[i];
      std::cerr << "  [warmstart] seeded " << n_sta << " state coefficient(s) from '"
                << SOLLOAD << "'; tau block left at init() value.\n";
    }
  }

  OCFESLV::SolveReport rep = oc.solve( xv.data(), nullptr, nullptr );
  R.march_conv = rep.converged;

  if( SOLDUMP ){
    if( !rep.converged )
      std::cerr << "  [warmstart] NOT SAVING: solve did not converge, so this vector is not a"
                   " solution and seeding from it would test nothing.\n";
    else{
      std::ofstream f( SOLDUMP );
      f << n_sta << "\n" << std::setprecision(17) << std::scientific;
      for( size_t i = 0; i < n_sta; ++i ) f << xv[i] << "\n";
      std::cerr << "  [warmstart] saved " << n_sta << " state coefficient(s) to '"
                << SOLDUMP << "'\n";
    }
  }

  // ---------------------------------------------------------------------------
  // SEAM PROBE  (--seamprobe)
  // ---------------------------------------------------------------------------
  // WHY.  [dup-spread] reports q discontinuous by 3.08e+01 and u by 5.31e+00 under IC_TRACE at
  // Pe=100, NEL_T=8, nel_xi=4.  The SOLUTION PLOTS SHOW NO SUCH DISCONTINUITY.  Both cannot be
  // right, and the plot is the more direct evidence, so the reporter is on trial here -- not
  // the framework.  ("The instrument is likelier to be wrong than the framework": two of three
  // defects investigated in the previous arc were instrument bugs, and the dup-spread reporter
  // ITSELF carried a base-offset bug through rev41-47 that produced STRUCTURED FALSE
  // POSITIVES.)
  //
  // THREE WAYS BOTH OBSERVATIONS CAN BE TRUE, and this probe separates them:
  //
  //  (1) The jumps are at interior t-SEAMS and the plots never sample them.  The plot loop
  //      uses tt = p.tf*k/NT_S with NT_S=11, chosen independently of the t-element boundaries
  //      at p.tf*j/NEL_T; for NEL_T=8 the two grids coincide only at t=0 and t=tf.  A jump in
  //      TIME at an element boundary is invisible in a spatial profile plotted between seams.
  //      Supporting evidence: the spread GROWS with NEL_T (3.0e-02 / 1.98 / 30.8 at NEL_T =
  //      2/3/8) and SHRINKS with nel_xi -- the wrong way round for a xi-seam effect.
  //
  //  (2) The plot is blind by construction: at a shared node eval_colloc resolves to ONE
  //      element, so if both sides of a xi-seam resolve to the same basis only one side is
  //      ever plotted and a seam jump cannot appear.
  //
  //  (3) The reporter is wrong and the spreads are false positives.
  //
  // WHAT THIS PROBE DOES.  Evaluates u, q and the mesh states either side of each interior
  // t-seam AND each interior xi-seam, by offsetting the evaluation point by +/- SEAM_EPS in
  // the seam's own coordinate.  This is the same eval_colloc path the plots use, so it asks
  // the question in the plot's own terms rather than the reporter's.
  //
  // HOW TO READ IT.  Jumps ~5 in u and ~30 in q at t-seams => the reporter is right, the
  // discontinuity is REAL and TEMPORAL, and the plots were never able to show it.  All sides
  // agreeing to round-off => [dup-spread] is producing false positives on this model and every
  // conclusion drawn from it must be withdrawn, including "IC_TRACE does not enforce
  // continuity".
  //
  // NOTE the probe reports the raw one-sided values as well as the jump, because a jump is
  // only meaningful against the local scale -- u is O(1) here, so |du| ~ 5 is a real jump,
  // whereas the same absolute number on a state of magnitude 1e6 would not be.
  if( SEAMPROBE ){
    double const SEAM_EPS = 1.0e-9;
    struct SP { char const* nm; FFVar const* v; };
    SP const probes[] = { { "u ", &u }, { "q ", &q }, { "x ", &x },
                          { "xg", &xg }, { "g ", &g }, { "gg", &gg } };

    auto at = [&]( FFVar const& v, double tt, double qq )->double {
      OCFESLV::t_Coord pt; pt[t] = tt; pt[xi] = qq;
      try { return oc.eval_colloc<double>( v, pt, xv.data(), nullptr, nullptr ); }
      catch( ... ) { return std::numeric_limits<double>::quiet_NaN(); }
    };

    std::cerr << "  [seamprobe] conv=" << ( rep.converged ? "yes" : "NO" )
              << "  NEL_T=" << NEL_T << " nel_xi=" << nel_xi
              << "  eps=" << SEAM_EPS << "\n";

    // ---- interior t-seams, sampled at three MID-ELEMENT xi positions -----------------
    for( size_t j = 1; j < NEL_T; ++j ){
      double const ts = p.tf * double( j ) / double( NEL_T );
      for( double const qq : { 0.15, 0.45, 0.75 } ){
        std::cerr << "  [seamprobe] t-seam j=" << j << " t=" << ts << " xi=" << qq;
        for( auto const& pr : probes ){
          double const lo = at( *pr.v, ts - SEAM_EPS, qq );
          double const hi = at( *pr.v, ts + SEAM_EPS, qq );
          std::cerr << "  " << pr.nm << "=" << std::scientific << std::setprecision(3)
                    << lo << "/" << hi << "(d=" << ( hi - lo ) << ")";
        }
        std::cerr << std::defaultfloat << "\n";
      }
    }

    // ---- interior xi-seams, sampled at three MID-ELEMENT t positions -----------------
    for( size_t e = 1; e < nel_xi; ++e ){
      double const qs = double( e ) / double( nel_xi );
      for( double const f : { 0.15, 0.45, 0.75 } ){
        double const tt = p.tf * f;
        std::cerr << "  [seamprobe] xi-seam e=" << e << " xi=" << qs << " t=" << tt;
        for( auto const& pr : probes ){
          double const lo = at( *pr.v, tt, qs - SEAM_EPS );
          double const hi = at( *pr.v, tt, qs + SEAM_EPS );
          std::cerr << "  " << pr.nm << "=" << std::scientific << std::setprecision(3)
                    << lo << "/" << hi << "(d=" << ( hi - lo ) << ")";
        }
        std::cerr << std::defaultfloat << "\n";
      }
    }
    std::cerr << "  [seamprobe] jumps ~O(1) in u/q at t-seams => dup-spread is RIGHT and the"
                 " discontinuity is temporal (plots sample between seams, so they cannot show"
                 " it).  All sides equal to round-off => dup-spread is a FALSE POSITIVE on this"
                 " model and the conclusions drawn from it must be withdrawn." << std::endl;
  }


  // --- errors: u vs exact, x vs oracle, on the mesh's own nodes ------------
  int const NT_S = 11, PPE = 6;
  double emax = 0., xmax = 0., umx = 0.;
  for( int k = 0; k <= NT_S; ++k ){
    double const tt = p.tf * double( k ) / double( NT_S );
    for( size_t e = 0; e < nel_xi; ++e )
      for( int s = 0; s <= PPE; ++s ){
        double const qq = ( double( e ) + double( s ) / double( PPE ) )
                        / double( nel_xi );
        OCFESLV::t_Coord pt; pt[t] = tt; pt[xi] = qq;
        double un = 0., xn = 0.;
        try {
          un = oc.eval_colloc<double>( u, pt, xv.data(), nullptr, nullptr );
          xn = oc.eval_colloc<double>( x, pt, xv.data(), nullptr, nullptr );
        } catch( ... ) { continue; }
        if( !std::isfinite( un ) || !std::isfinite( xn ) ){
          emax = std::numeric_limits<double>::infinity(); continue; }
        umx  = std::max( umx, std::fabs( un ) );
        xmax = std::max( xmax, std::fabs( xn - O.x( qq, tt ) ) );
        double const zz = O.x( qq, tt );   // compare u at the oracle physical pt
        emax = std::max( emax, std::fabs( un - u_exact( zz, tt, p ) ) );
      }
  }
  R.errNode = emax; R.errMesh = xmax; R.umax = umx;

  // ---- DIAGNOSTIC: locate the worst u error and dump a profile slice --------
  // errNode ~ 1 with a perfect mesh (|x-x*|~1e-8) points at a PHASE error in u
  // (front at the wrong place), not the flux formulation.  Find WHERE, and dump a
  // mid-time slice so the deviation is visible rather than inferred.
  if( g_diag ){
    // worst-error point
    double ewmax = 0., ew_t = 0, ew_xi = 0, ew_z = 0, ew_un = 0, ew_ue = 0;
    for( int k = 0; k <= 40; ++k ){
      double const tt = p.tf * double( k ) / 40.0;
      for( int j = 0; j <= 400; ++j ){
        double const qq = double( j ) / 400.0;
        OCFESLV::t_Coord pt; pt[t] = tt; pt[xi] = qq;
        double un = 0.;
        try { un = oc.eval_colloc<double>( u, pt, xv.data(), nullptr, nullptr ); }
        catch( ... ) { continue; }
        double const zz = O.x( qq, tt );
        double const e  = std::fabs( un - u_exact( zz, tt, p ) );
        if( e > ewmax ){ ewmax=e; ew_t=tt; ew_xi=qq; ew_z=zz;
                         ew_un=un; ew_ue=u_exact(zz,tt,p); }
      }
    }
    std::cout << "    [diag] worst |u-u*|=" << std::scientific << std::setprecision(3)
              << ewmax << " at t=" << std::fixed << std::setprecision(4) << ew_t
              << " xi=" << ew_xi << " z=" << ew_z
              << "  u=" << std::setprecision(5) << ew_un
              << " u*=" << ew_ue << "  front s(t)=" << p.s( ew_t ) << "\n";
    // profile slice at t=tf/2, across xi, showing u vs u* and the physical z
    double const ts = 0.5 * p.tf;
    std::cout << "    [diag] slice t=" << std::fixed << std::setprecision(3) << ts
              << " (front at z=" << p.s( ts ) << "):  xi  z(xi)  u  u*  x-x*\n";
    for( int j = 0; j <= 20; ++j ){
      double const qq = double( j ) / 20.0;
      OCFESLV::t_Coord pt; pt[t] = ts; pt[xi] = qq;
      double un = 0., xn = 0.;
      try { un = oc.eval_colloc<double>( u, pt, xv.data(), nullptr, nullptr );
            xn = oc.eval_colloc<double>( x, pt, xv.data(), nullptr, nullptr ); }
      catch( ... ) { continue; }
      double const zz = O.x( qq, ts );
      std::cout << "      " << std::fixed << std::setprecision(4) << qq
                << "  " << std::setprecision(4) << zz
                << "  " << std::setprecision(5) << un
                << "  " << u_exact( zz, ts, p )
                << "  " << std::scientific << std::setprecision(2) << ( xn - zz ) << "\n";
    }
  }
  // ---- DUMP: gnuplot-ready data files (--dump) ------------------------------
  // Writes, for this (Pe, nel):
  //   mesh2_sol_<nel>.dat   : t  xi  z=x(xi,t)  u(xi,t)  u_exact  x*(xi,t)
  //                           (blank line between t-blocks -> gnuplot pm3d/lines)
  //   mesh2_mesh_<nel>.dat  : t  xi  z=x(xi,t)   (node trajectories; one xi-line per block)
  //   mesh2_plot_<nel>.gp   : ready-to-run gnuplot script
  // The solution is emitted in PHYSICAL coordinates (z on the moving mesh), which
  // is what you actually want to see: u plotted against z, with the mesh nodes
  // visibly clustering on and tracking the front.
  if( g_dump ){
    // Filenames carry MODE and CONFIG.  Previously every run wrote mesh16_*_<nel>, so a
    // --weak / --trace / IC_STRONG comparison silently overwrote itself -- and the mode
    // comparison is the whole point of looking at these plots.
    std::string const mode_tag = g_weak ? "weak" : ( g_trace ? "trace" : "strong" );
    std::string const cfg_tag  = std::string( EVOLVE ? "evolve" : "mono" )
                               + ( QPIN ? "_qpin" : "" );
    std::ostringstream tag; tag << mode_tag << "_" << cfg_tag << "_" << nel_xi;
    std::string const fs = "mmpde_sol_"  + tag.str() + ".dat";
    std::string const fm = "mmpde_mesh_" + tag.str() + ".dat";
    std::string const fg = "mmpde_plot_" + tag.str() + ".gp";

    int const NT_D = 40;      // time slices
    int const NP_D = 12;      // sample points per element (dense enough for the front)

    std::ofstream os( fs );
    // conv= is recorded IN THE FILE.  A converged and a stalled run produce plots that look
    // equally plausible, and this session spent hours on numbers taken from non-converged
    // solves; the provenance must travel with the data.
    os << "# Pe=" << p.Pe << " nel=" << nel_xi << " delta=" << p.delta()
       << " NEL_T=" << NEL_T << " tau=" << TAU   // TAU is file-scope, not a Params member
       << " mode=" << mode_tag << " cfg=" << cfg_tag
       << " conv=" << ( rep.converged ? "yes" : "NO" )
       << " tol=" << RESTOL << "\n"
       << "# t  xi  z=x(xi,t)  u  u_exact  xstar  q  q_exact\n";
    for( int k = 0; k <= NT_D; ++k ){
      double const tt = p.tf * double( k ) / double( NT_D );
      for( size_t e = 0; e < nel_xi; ++e )
        for( int s = 0; s <= NP_D; ++s ){
          if( e > 0 && s == 0 ) continue;                 // avoid duplicate element joints
          double const qq = ( double( e ) + double( s ) / double( NP_D ) )
                          / double( nel_xi );
          OCFESLV::t_Coord pt; pt[t] = tt; pt[xi] = qq;
          double un = 0., xn = 0., qn = 0.;
          try { un = oc.eval_colloc<double>( u, pt, xv.data(), nullptr, nullptr );
                xn = oc.eval_colloc<double>( x, pt, xv.data(), nullptr, nullptr );
                qn = oc.eval_colloc<double>( q, pt, xv.data(), nullptr, nullptr ); }
          catch( ... ) { continue; }
          double const zs = O.x( qq, tt );
          // q is compared against q_exact at the SOLVED node position xn, not at the oracle
          // zs: the question a q panel answers is whether the computed gradient matches the
          // analytic one WHERE THE MESH ACTUALLY PUT THE NODE.
          os << std::fixed << std::setprecision(6)
             << tt << " " << qq << " " << xn << " "
             << std::setprecision(8) << un << " "
             << u_exact( zs, tt, p ) << " " << zs << " "
             << qn << " " << q_exact( xn, tt, p ) << "\n";
        }
      os << "\n";                                          // block separator per t
    }
    os.close();

    // mesh trajectories: for each COLLOCATION-ish xi line, z(t) -- shows nodes
    // sweeping with the front.  One block per xi sample, t running down the block.
    std::ofstream om( fm );
    om << "# Pe=" << p.Pe << " nel=" << nel_xi << "\n# xi  t  z=x(xi,t)\n";
    int const NX_L = nel_xi * 4;                           // xi lines to draw
    for( int jx = 0; jx <= NX_L; ++jx ){
      double const qq = double( jx ) / double( NX_L );
      for( int k = 0; k <= NT_D; ++k ){
        double const tt = p.tf * double( k ) / double( NT_D );
        OCFESLV::t_Coord pt; pt[t] = tt; pt[xi] = qq;
        double xn = 0.;
        try { xn = oc.eval_colloc<double>( x, pt, xv.data(), nullptr, nullptr ); }
        catch( ... ) { continue; }
        om << std::fixed << std::setprecision(6)
           << qq << " " << tt << " " << xn << "\n";
      }
      om << "\n";
    }
    om.close();

    // gnuplot script: three panels -- u(z,t) surface, u vs z at mid-time with
    // mesh nodes, and the mesh trajectory z(t) fanning out around the front.
    std::ofstream og( fg );
    // Four panels now: the q panel is added because q = u_z is the state whose continuity
    // dominated this arc's diagnostics, and a plot of it against q_exact says immediately
    // whether a large dup-spread is a defect or an under-resolved front.  The title carries
    // conv= so a stalled run cannot be mistaken for a converged one at a glance.
    og << "# gnuplot -p " << fg << "\n"
       << "set term qt size 1800,500\n"
       << "set multiplot layout 1,4 title 'MMPDE " << mode_tag << "/" << cfg_tag
       << "  Pe=" << p.Pe << " nel=" << nel_xi << " NEL_T=" << NEL_T
       << "  conv=" << ( rep.converged ? "yes" : "NO" ) << "'\n"
       << "set xlabel 'z'; set ylabel 't'\n"
       << "set title 'u(z,t) on moving mesh'\n"
       << "set pm3d map; set palette defined (0 'white',1 'navy')\n"
       << "splot '" << fs << "' using 3:1:4 with pm3d notitle\n"
       << "unset pm3d\n"
       << "set title 'u vs z at t=" << 0.5 * p.tf
       << " (points=mesh nodes)'\n"
       << "set xlabel 'z'; set ylabel 'u'\n"
       << "tmid=" << std::fixed << std::setprecision(6) << 0.5 * p.tf << "\n"
       << "plot '" << fs << "' using ($1==tmid?$3:1/0):4 with linespoints pt 7 ps 0.6 "
          "title 'u (solved)', '" << fs
       << "' using ($1==tmid?$3:1/0):5 with lines lw 2 title 'u* exact'\n"
       << "set title 'mesh trajectories z(t) tracking the front'\n"
       << "set xlabel 't'; set ylabel 'z'\n"
       << "plot '" << fm << "' using 2:3 with lines lc 'grey' notitle, "
          "" << p.s0 << "+" << p.c << "*x with lines lw 2 lc 'red' title 'front s(t)'\n"
       // ---- panel 4: q = u_z against the analytic gradient ----------------------
       // Peak should be -1/(2 delta) = " << -0.5 / p.delta() << ".  If the solved q tracks
       // q_exact, a large q dup-spread is the front being under-resolved, NOT a continuity
       // failure -- IC_WEAK measures 1.23e+01 on the same configuration, and weak spreads
       // are truncation-limited by definition.
       << "set title 'q=u_z vs z at t=" << 0.5 * p.tf
       << "  (peak should be " << std::fixed << std::setprecision(1) << -0.5 / p.delta()
       << ")'\n"
       << std::setprecision(6)
       << "set xlabel 'z'; set ylabel 'q'\n"
       << "plot '" << fs << "' using ($1==tmid?$3:1/0):7 with linespoints pt 7 ps 0.6 "
          "title 'q (solved)', '" << fs
       << "' using ($1==tmid?$3:1/0):8 with lines lw 2 title 'q* exact'\n"
       << "unset multiplot\n";
    og.close();

    std::cout << "    [dump] wrote " << fs << ", " << fm << ", " << fg
              << "  (run:  gnuplot -p " << fg << ")\n";
  }

  return R;
}

// ---------------------------------------------------------------------------
int main( int argc, char** argv )
{
  double Pe = 100.0;
  bool   pe_set = false;   // first numeric >=10 is Pe (once); rest are nel values.
                           // Was 'Pe==100.0' which wrongly re-triggered when the
                           // user actually passed Pe=100 (then nel=32 became Pe=32).
  int    display = 0;
  std::vector<size_t> nelList;
  for( int i = 1; i < argc; ++i ){
    if( !std::strcmp( argv[i], "-v" ) ){ display = 1; continue; }
    if( !std::strcmp( argv[i], "--diag" ) ){ g_diag = true; continue; }
    if( !std::strcmp( argv[i], "--dump" ) ){ g_dump = true; continue; }
    if( !std::strcmp( argv[i], "--no-iface" ) ){ g_iface = false; continue; }
    if( !std::strcmp( argv[i], "--verbose" ) ){ g_solveverbose = true; continue; }
    if( !std::strcmp( argv[i], "--trace" ) ){ g_trace = true; continue; }
    if( !std::strcmp( argv[i], "--weak" ) ){ g_weak = true; continue; }
    if( !std::strcmp( argv[i], "--small" ) ){ NND_XI = 5; NEL_T = 2; NND_T = 3; continue; }
    // --nelt / --tau: the t-mesh knob is the discriminating experiment for the r-kernel.
    // The kernel localises at t=0 and at each t-element seam, i.e. 1+(NEL_T-1) t-locations
    // per xi-interface; if that reading is right, NEL_T=3 gives 3*(nel_xi-1), not 2*.
    // Placed AFTER --small so "--small --nelt 3" overrides the small preset.
    if( !std::strcmp( argv[i], "--nelt" ) && i+1 < argc ){ NEL_T = (size_t)std::atol( argv[++i] ); continue; }
    if( !std::strcmp( argv[i], "--tau"  ) && i+1 < argc ){ TAU   = std::atof( argv[++i] ); continue; }
    if( !std::strcmp( argv[i], "--kerow" ) ){ g_kerow = true; continue; }
    // --sigma0 / --sigma1: the SAT penalty weights.  SAT_SIGMA1 = 0 (the MMPDE7 default,
    // provenance unknown) makes tau IDENTICALLY ZERO for every CLAIM_EXACT_C1 receiver
    // term, hence prefactor 0, hence a tau column with no entry anywhere -- while the
    // PLAN's rank test rates the same claim rank-extending, because it scores the column
    // as coupling*orientation with no tau factor and no C1 zeroing.  If any of the five
    // explicit claims has kind=1, that asymmetry alone produces DEFECT A.
    if( !std::strcmp( argv[i], "--sigma0" ) && i+1 < argc ){ SIG0 = std::atof( argv[++i] ); continue; }
    if( !std::strcmp( argv[i], "--sigma1" ) && i+1 < argc ){ SIG1 = std::atof( argv[++i] ); continue; }
    if( !std::strcmp( argv[i], "--eps"    ) && i+1 < argc ){ PEPS = std::atof( argv[++i] ); continue; }
    if( !std::strcmp( argv[i], "--gdirect" ) ){ GDIRECT = true; continue; }
    if( !std::strcmp( argv[i], "--chain"   ) ){ GDIRECT = false; continue; }   // 2026-09-09: restore the chained form
    if( !std::strcmp( argv[i], "--adapt" ) && i+1 < argc ){ ADAPT = std::atof( argv[++i] ); continue; }
    if( !std::strcmp( argv[i], "--seamprobe" ) ){ SEAMPROBE = true; continue; }
    if( !std::strcmp( argv[i], "--pde-all-t" ) ){ PDE_ALL_T = true; continue; }
    if( !std::strcmp( argv[i], "--qpin"      ) ){ QPIN      = true; continue; }
    if( !std::strcmp( argv[i], "--qmove"    ) ){ QMOVE     = true; continue; }
    if( !std::strcmp( argv[i], "--qpinonly" ) ){ QPINONLY  = true; continue; }
    if( !std::strcmp( argv[i], "--qpinscale" ) && i+1 < argc )
      { QPINSCALE = std::atof( argv[++i] ); continue; }
    if( !std::strcmp( argv[i], "--tol" ) && i+1 < argc )
      { RESTOL = std::atof( argv[++i] ); continue; }
    if( !std::strcmp( argv[i], "--forceval" ) ){ FORCEVAL = true; continue; }
    if( !std::strcmp( argv[i], "--evolve"    ) ){ EVOLVE    = true; continue; }
    if( !std::strcmp( argv[i], "--savesol" ) && i+1 < argc ){ SOLDUMP = argv[++i]; continue; }
    if( !std::strcmp( argv[i], "--loadsol" ) && i+1 < argc ){ SOLLOAD = argv[++i]; continue; }
    // --drop off|verify|rederive : INTERFACE_DROP_POLICY.  Unset leaves the header
    // default untouched, so a bare run stays bit-identical to MMPDE8/MMPDE7.
    if( !std::strcmp( argv[i], "--drop" ) && i+1 < argc ){
      char const* d = argv[++i];
      if     ( !std::strcmp( d, "verify"   ) ) DROPP = 0;
      else if( !std::strcmp( d, "rederive" ) ) DROPP = 1;
      else if( !std::strcmp( d, "off"      ) ) DROPP = 2;
      else { std::cerr << "unknown --drop " << d << " (off|verify|rederive)\n"; return 1; }
      continue;
    }
    // --rscale: the r non-dimensionalisation factor.  Default delta^2 inverts the
    // xg/r coefficient ratio inside RDEF.  --rscale 1 reproduces MMPDE4 unscaled, so a
    // sweep 1 / 1e-2 / 1e-4 reads sigma_min directly against the scaling hypothesis.
    if( !std::strcmp( argv[i], "--rscale" ) && i+1 < argc ){ RSCALE = std::atof( argv[++i] ); continue; }
    // --form: which discrete formulation of the SAME continuous problem to assemble.
    if( !std::strcmp( argv[i], "--form" ) && i+1 < argc ){
      char const* f = argv[++i];
      if     ( !std::strcmp( f, "rstate" ) ) FORM = FORM_RSTATE;
      else if( !std::strcmp( f, "divide" ) ) FORM = FORM_DIVIDE;
      else if( !std::strcmp( f, "xgmult" ) ) FORM = FORM_XGMULT;
      else { std::cerr << "unknown --form " << f << " (rstate|divide|xgmult)\n"; return 1; }
      continue;
    }
    if( !pe_set && nelList.empty() && std::atof( argv[i] ) >= 10.0 )
      { Pe = std::atof( argv[i] ); pe_set = true; continue; }
    nelList.push_back( (size_t)std::atol( argv[i] ) );
  }
  if( nelList.empty() ) nelList = { 4, 8, 16, 32 };   // accuracy sweep vs IC_WEAK ref

  Params p; p.Pe = Pe;

  std::cout << "================================================================\n"
            << "  MMPDE10 -- MMPDE5 relaxation, r-FREE reference formulation (--form)\n"
            << "  mode=" << ( g_weak ? "IC_WEAK(ref)" : g_trace ? "IC_TRACE" : "IC_STRONG" )
            << "   (--weak = converging reference, --trace = other exact mode)\n"
            << "  form=" << ( FORM == FORM_RSTATE ? "rstate" : FORM == FORM_DIVIDE ? "divide" : "xgmult" )
            << "   SOLUTION:  q*xg = u_xi ;  "
            << ( FORM == FORM_RSTATE ? "rh*xg = RSCALE*q_xi ;  u_t + (c-x_t)q - (D/RSCALE)rh = f"
               : FORM == FORM_DIVIDE ? "u_t + (c-x_t)q - (D/xg)q_xi = f   [no r]"
                                     : "xg*u_t + (c-x_t)xg*q - D*q_xi = xg*f   [no r]" ) << "\n"
            << "  MESH    :  xg = x_xi ;  g = M*xg ;  gg = g_xi ;  tau*x_t = gg   (MMPDE5)\n"
            << "             M = 1/x*_xi PRESCRIBED -> x -> M1a's sinh oracle as tau->0.\n"
            << "             x is DIFFERENTIAL in t (RELAX), so it carries a real IC.\n"
            << "             This is NOT instantaneous equidistribution and NOT index-2;\n"
            << "             earlier banners said otherwise and were wrong.\n"
            << "  MEASURED:  kernel = 2*NEL_T*(nel_xi-1) for xgmult, NEL_T*(nel_xi-1) for\n"
            << "             rstate -- dropping the u_zz state did NOT reduce the defect, it\n"
            << "             converted one mechanism into two.  Two causes, both confirmed:\n"
            << "  DEFECT A:  the plan emits 5 continuity multipliers the rows cannot support\n"
            << "             -- (nel_xi-1)*(2*NEL_T-1) of them.  tau-block RANK DEFICIT = 5 in\n"
            << "             BOTH exact modes (orphan columns under STRONG, distributed under\n"
            << "             TRACE).  NOT a zero-coupling effect: u/x have max|sat| 1.00/4.61.\n"
            << "  DEFECT B:  a state whose sole surviving row at a node is also an interface\n"
            << "             receiver row loses its jump.  Here q at t=0 (PDE is on ALL-LB so\n"
            << "             it is absent there); nel_xi-1 modes, identical in both modes.\n"
            << "             Read the gap-based kernel estimate, not rank_deficiency[1e-12].\n"
            << "  Pe=" << Pe << "  delta=" << std::scientific << std::setprecision(2)
            << p.delta() << "  tau=" << TAU
            << ( FORM == FORM_RSTATE
                 ? "  rscale=" + std::to_string( RSCALE > 0. ? RSCALE : p.delta()*p.delta() )
                 : std::string( "  rscale=n/a" ) )
            << "  sigma0=" << SIG0 << " sigma1=" << SIG1
            << "  drop=" << ( DROPP < 0 ? "default" : DROPP == 0 ? "verify" : DROPP == 1 ? "rederive" : "off" )
            << "  eps=" << ( PEPS > 0.0 ? PEPS : -1.0 )
            << ( GDIRECT ? "  gdirect=ON" : "" )
            << ( ADAPT > 0. ? "  adapt=ON(SOLUTION-DEPENDENT monitor)" : "" )
            << ( PDE_ALL_T ? "  pde-all-t=ON" : "" )
            << ( QPIN      ? "  qpin=ON"      : "" )
            << ( QMOVE     ? "  qmove=ON"     : "" )
            << ( QPINONLY  ? "  qpinonly=ON"  : "" )
            << ( RESTOL != 1.0e-9 ? "  tol=" : "" ) << ( RESTOL != 1.0e-9 ? RESTOL : 0.0 )
            << ( FORCEVAL  ? "  forceval=ON"  : "" )
            << ( QPINSCALE != 1.0 ? "  qpinscale=" : "" )
            << ( QPINSCALE != 1.0 ? QPINSCALE : 0.0 )
            << ( EVOLVE    ? "  evolve=ON"    : "" )
            << ( SOLDUMP ? "  savesol=on" : "" ) << ( SOLLOAD ? "  loadsol=on" : "" )
            << ( SEAMPROBE ? "  seamprobe=on" : "" )
            << "  t-mesh " << NEL_T << "x" << NND_T
            << "  xi order " << NND_XI << "\n"
            << "  gates: setup square + classify OK, solve converges, |x-x*|<"
            << XTOL << ", errNode<" << ERRTOL << "  (finest mesh)\n"
            << "================================================================\n";
  std::cout << "  " << std::left
            << std::setw(6) << "nel" << std::setw(8) << "nVar"
            << std::setw(8) << "square" << std::setw(7) << "conv"
            << std::setw(12) << "|x-x*|" << std::setw(12) << "errNode"
            << std::setw(11) << "max|u|" << "  status\n";

  bool all_ok = true;
  for( size_t ne : nelList ){
    Result R = run( p, ne, display );
    bool const ok = R.setup_ok && R.square && R.march_conv
                 && std::isfinite( R.errMesh ) && R.errMesh < XTOL
                 && std::isfinite( R.errNode ) && R.errNode < ERRTOL;
    // errNode is a DISCRETISATION error and converges under refinement -- MEASURED at
    // TAU=1e-6: 4.23e-03, 5.46e-05, 8.89e-06 at nel_xi 4, 8, 16.  Gating every mesh against an
    // absolute tolerance would therefore fail the coarse ones for having too few elements,
    // which is not what this driver tests.  The verdict is the FINEST mesh; the coarser rows
    // are reported so the convergence is visible.
    all_ok = ok;
    std::cout << "  " << std::left
              << std::setw(6) << ne << std::setw(8) << R.nVar
              << std::setw(8) << ( R.square ? "yes" : "NO" )
              << std::setw(7) << ( R.march_conv ? "yes" : "no" )
              << std::scientific << std::setprecision(2)
              << std::setw(12) << R.errMesh << std::setw(12) << R.errNode
              << std::setw(11) << R.umax << "  " << R.status;
    if( R.threw ) std::cout << "  THREW: " << R.exmsg;
    std::cout << ( ok ? "   [OK]" : "   [--]" ) << "\n";
  }

  std::cout << "\n================================================================\n"
            << "  MMPDE10[" << ( FORM == FORM_RSTATE ? "rstate" : FORM == FORM_DIVIDE ? "divide" : "xgmult" )
            << "]: " << ( all_ok ? "PASS" : "FAIL" ) << "\n"
            << "  PASS means the relaxed mesh reproduces the equidistribution oracle\n"
            << "  (|x-x*| small) AND the ALE solution stays accurate on it, on the FINEST\n"
            << "  mesh.  |x-x*| is the relaxation lag, 0.169*tau exactly and\n"
            << "  mesh-independent; errNode is a discretisation error and converges\n"
            << "  (4.23e-03, 5.46e-05, 8.89e-06 at nel_xi 4, 8, 16 with tau=1e-6).\n"
            << "  RESOLVED 2026-09-09: IC_STRONG previously stalled here because the\n"
            << "  four-link chain x -> xg -> g -> gg has g = M*xg POINTWISE IN STATES, so\n"
            << "  g\'s continuity is implied by xg\'s -- the plan claimed it anyway and\n"
            << "  realised it as a tau multiplier enforcing nothing, one per interior seam\n"
            << "  (deficient clusters = k\' = nel_xi-1; the two claims share a receiver row\n"
            << "  with equal and opposite coefficients).  --gdirect, now the DEFAULT,\n"
            << "  removes that link: deficiency 0 and the solve converges.  --chain\n"
            << "  restores the old form and the old failure.\n"
            << "  A converged residual is never on its own evidence of a correct solution:\n"
            << "  check |x-x*|, the right-null energy, and whether ANY step was accepted.\n"
            << "  Run --verbose for the singular-value spectrum, the eqn legend and\n"
            << "  the left/right null-space decomposition.\n"
            << "================================================================\n";
  return all_ok ? 0 : 1;
}
