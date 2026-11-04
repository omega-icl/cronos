// Copyright (C) Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the EPL-2.0 with GPL-2.0-or-later as a Secondary License; see the LICENSE file.
// SPDX-License-Identifier: EPL-2.0 OR GPL-2.0-or-later

#ifndef CRONOS__OCENV_HPP
#define CRONOS__OCENV_HPP

// <windows.h>, in its FULL form (not WIN32_LEAN_AND_MEAN), leaves INTERFACE defined as a macro; this header
// declares OCFESLV::Options::INTERFACE.  Set aside for the header and restored after it (2026-10-09).
#if defined(_WIN32)
# pragma push_macro("INTERFACE")
# undef INTERFACE
#endif

// Hyperbolic boundary conditions.  The incoming-BC direction guard and its flux/derivative-BC extension are always
// on (unconditional since rev319; validated on PDE14, PDE15, PDE16 and PDE18): an inflow condition must pin an
// incoming characteristic through its value OR its face-normal derivative, and is rejected only when it touches none.
// The boundary closures that complete a face to square are environment-only knobs, both default on:
// CRONOS_AUTO_HYP_CLOSURE (outgoing characteristics, hyperbolic blocks; appends only the per-face deficit since
// rev317) and CRONOS_AUTO_ALG_CLOSURE (algebraic constraints, parabolic/DAE blocks).  The symbol-balance probe
// (MC__OCFESLV_SYMBOL_BALANCE_PROBE) stays an opt-in, print-only compile flag.

#include "ocbase.hpp"
#include "ffmodel.hpp"   // rev195: the model layer (class FFModel)
#include "ocplan.hpp"
#include "ocfebase.hpp"  // the direction-agnostic base: OCVar environment, interface plan, kernels
// rev265: the ocplan family is one header now   // rev152a: plan decision-layer contract (phase 0)
// rev264: ocplan_build_legacy.hpp DELETED (954 lines) -- its builder was never applied.
   // rev153a: GLOBAL builder (phase 1a: shadow only)
         // rev153a: PlanReport printer
#include "ffinv.hpp"    // rev82: invertible-structure detection (L/S/N/U) for the claim-
                        // inheritance test.  FFDep cannot serve: it types the WHOLE
                        // dependence, so g in DEF (g - W*xg) reports xg=N merely because the
                        // monitor W = 1/x*_xi inlines the sinh oracle (33 DAG ops) -- MEASURED
                        // via [depdump].  FFInv's S ("separably linear") is the property that
                        // matters: xg appears once, multiplied by a coefficient, so the row is
                        // invertible for xg whatever that coefficient contains.

#include <chrono>       // setup() phase timing (gated behind DISPLAY_LEVEL); homotopy stage timing
#include <functional>   // std::function -- homotopy s -> value maps (add_homotopy)
#include <tuple>
#include <set>          // rev12: projection fill-in -- support unions for C_tilde
#include <algorithm>    // std::find/sort/copy/max -- homotopy stage ordering and stepping
#include <sstream>      // std::ostringstream -- decode_state_column / equation legend formatting
#include <cstdlib>      // std::getenv -- environment defaults for the AUDIT_* options
#include <string>       // std::string -- environment flag parsing
#include <map>
#include <mutex>
#include <memory>          // rev41: geometric coincidence grouping in duplicate_node_spread
#include <iomanip>      // rev41: std::setw/setprecision in report_duplicate_node_spreads
#include <cmath>        // rev41: std::llround -- node-coordinate quantisation (1e-12)
#include <iostream>
#include <fstream>     // rev41: std::ostream default argument on the spread reporter
#include <cstdint>      // rev142: int64_t in the _staircase_live_columns DECLARATION, which
                        // sits in the class body and so is compiled even when neither the
                        // SPQR nor the UMFPACK block (the only other includers of <cstdint>)
                        // is active.  Relying on a transitive include here would make a
                        // no-backend build fail on a line that has nothing to do with the
                        // missing backend.

// POSIX file-descriptor headers used by OCFESLV::solve() to scope-silence the
// benign reference-BLAS chatter emitted by a direct sparse Newton solve on a
// rank-deficient interface Jacobian.  Guarded so non-POSIX builds still compile
// (the guard degrades to a no-op there).
#if defined(__unix__) || defined(__APPLE__)
# include <cstdio>
# include <fcntl.h>
# include <unistd.h>
# include <dlfcn.h>     // rev142b: MAXTHREAD resolves the OpenMP/BLAS thread setters by
                        // NAME at runtime rather than declaring them.  Declaring them would
                        // add a link dependency and could collide with <omp.h> or a vendor
                        // BLAS header; looking them up cannot.  glibc >= 2.34 carries dlsym
                        // in libc, so no -ldl is required.
#endif

// SuiteSparseQR (sparse rank-revealing QR) backend for OCFESLV::solve()'s
// SOLVE_SPQR factorization option.  Entirely opt-in behind CRONOS__WITH_SPQR:
// compile with -DCRONOS__WITH_SPQR and link -lspqr -lcholmod -lsuitesparseconfig
// -lamd -lcolamd -lccolamd -lcamd (plus BLAS/LAPACK).  When the macro is absent
// the SOLVE_SPQR enum value, the QR helper, and all QR branches are excised and
// the solver is exactly the current SuperLU/LAPACK Newton/LM path (no SuiteSparse
// dependency whatsoever).
#if defined(CRONOS__WITH_SPQR)
# include <SuiteSparseQR.hpp>   // pulls in cholmod.h
# include <cstring>             // std::memset for the cholmod struct wrappers
# include <cstdint>            // int64_t index arrays for the long cholmod API
#endif

// UMFPACK (sparse unsymmetric multifrontal LU, AMD fill-reducing ordering) backend
// for the IC_STRONG tau-elimination's B_S factorization.  Opt-in behind
// CRONOS__WITH_UMFPACK: compile with -DCRONOS__WITH_UMFPACK and link -lumfpack -lamd
// (-lcholmod -lsuitesparseconfig are already pulled in by the SPQR link above).
// When the macro is absent the tau elimination uses the dense Armadillo inv path
// (no UMFPACK dependency).  Used only at setup to build the frozen projection.
#if defined(CRONOS__WITH_UMFPACK)
# include <umfpack.h>
# include <cstdint>            // int64_t index arrays for the long UMFPACK API
#endif

// rev133: optional Eigen SparseQR backend (MPL2, header-only) for the setup-time
// determinacy audit's numerical rank.  SPQR is preferred when built; this is the
// permissive fallback.  Neither is required: a dense SVD under DETERMINACY_MAX_DENSE
// is the last resort.
#if defined(CRONOS__WITH_EIGEN)
# include <Eigen/Sparse>
# include <Eigen/SparseQR>
# include <Eigen/OrderingMethods>
#endif

namespace mc
{
//! @brief Orthogonal collocation on finite elements (OCFE) solver for IPDAE models.
////////////////////////////////////////////////////////////////////////
//! mc::OCFESLV discretises and solves integro-partial-differential-algebraic equation (IPDAE) models by
//! orthogonal collocation on finite elements.  Expressions are MC++ DAGs.
//!
//! WHAT LIVES WHERE.  The split with the base class is now settled, and it is the map to this file:
//!   mc::FFModel   THE MODEL.  Declarations (constants, domains, states, inputs, equations, outputs), the
//!                 principal symbol, the PDE classification, order reduction, the continuity-claim
//!                 classification, and the setup() skeleton.  It is mesh-free: it sets up and classifies
//!                 with no discretisation, which is what lets a model be analysed without a solver.
//!   mc::OCFESLV   THIS CLASS.  Collocation and the node cache, the interface plan (claims -> realisations
//!                 -> rows and columns), and the Newton, marching, continuation and sensitivity solvers.
//!                 It contributes to setup() through the hooks FFModel calls, never by overriding it.
//!
//! So a question about WHAT THE MODEL MEANS is answered in ffmodel.hpp; a question about HOW IT IS
//! DISCRETISED OR SOLVED is answered here.
//!
//! ON THE rev### TAGS in this file: they are PROVENANCE, not contract -- each marks a change and, where it
//! matters, the measurement that justified it.  Read the surrounding sentence for what the code does; the
//! tag tells you when and why it came to be that way.  A tag that only names a revision carries nothing and
//! is being removed as the code around it is revisited.


////////////////////////////////////////////////////////////////////////
class OCFESLV
////////////////////////////////////////////////////////////////////////
: public OCFEBASE
{

  friend  std::ostream& operator<< ( std::ostream&, OCFESLV const& );
  template <typename U> friend class OCVar;

public:
  //! @brief The solver's options: OCFEBASE::Options -- INTERFACE, OUTPUT, AUDIT, DETERMINACY, FATAL, the
  //! groups the collocation base reads -- plus SOLVE, SENSITIVITY and HOMOTOPY, which only the monolithic /
  //! marching solver reads.  OCFEMOL carries none of these three (MoL phase 3, 31 Oct 2026).
  struct Options
  : public OCFEBASE::Options
  {
    enum DerivCacheType
    {
      DERIV_LAZY  = 0,
      DERIV_EAGER = 1
    };
    enum LinearEvalCacheType
    {
      LINEAR_NONE  = 0,
      LINEAR_LAZY  = 1,
      LINEAR_EAGER = 2
    };
    //! @brief Sparse linear factorization backend for the internal LM/Newton solve.
    //! SOLVE_SUPERLU uses arma::spsolve with SuperLU (requires the Armadillo build to
    //! define ARMA_USE_SUPERLU); when it is not linked the solver transparently falls
    //! back to a densifying LAPACK solve (one-time notice).  SOLVE_LAPACK forces the
    //! densifying dense-LU path (small/debug systems).  The enum is left extensible for
    //! a future UMFPACK/other direct backend without touching call sites.
    enum SolveFactorization
    {
      SOLVE_SUPERLU = 0,  //!< sparse SuperLU (arma::spsolve), default. The solver only ever
                          //!< factorizes the SPD/nonsingular damped operator (Js^T Js + lambda*D),
                          //!< so there are no zero-pivot issues even on rank-deficient interface
                          //!< Jacobians; sparse + fast on the large multi-block systems. Falls
                          //!< back to dense LAPACK (with a one-time notice) if not linked.
      SOLVE_LAPACK  = 1,  //!< densify the damped operator and use dense LU (debug / small systems)
      SOLVE_SPQR    = 2   //!< sparse rank-revealing QR (SuiteSparseQR) min-2-norm least squares on
                          //!< Js directly -- no Js^T Js, so no fill/condition squaring, and the
                          //!< rank-revealing factor resolves the singular interface-augmented saddle
                          //!< natively (the sparse analogue of the dense min-norm tier).  Requires a build
                          //!< with SuiteSparseQR (-DCRONOS__WITH_SPQR, the GPL build); otherwise setup()
                          //!< REFUSES it (BACKEND_UNAVAILABLE).  Always defined (2026-10-09), so a model
                          //!< compiles -- and Python lists it -- in every build.
    };
    //! @brief Warm-start strategy for the per-element step of a SOLVE.MARCHING solve.
    enum SolveWarmstart
    {
      REUSE       = 0,  //!< reuse the previous element's converged state as the next guess (free; default)
      BROADCAST_IC = 1, //!< seed every evolution node of the next element with the transferred IC profile
      EXTRAPOLATE = 2   //!< linear/quadratic extrapolation from the last one/two elements' terminals
    };
    enum NT_SKIP_MODE { NT_SKIP_ALWAYS=0, NT_SKIP_AUTO=1, NT_SKIP_NEVER=2 };

    //! @brief The nonlinear solve: Newton/Levenberg-Marquardt control, linear backend, marching, reuse, evaluation caches, and the solver's own verbosity.
    struct t_Solve
    {
      SolveFactorization FACTORIZATION;  //!< linear backend (default SOLVE_SUPERLU; see note)
      bool   MARCHING;       //!< auto-march when AVAILABLE (evolution domain + >=1 differential state); else monolithic (default true)
      double RES_TOL;        //!< accept on max|r| (1e-9)
      //! @name Internal nonlinear-solve (LM/Newton) controls — see OCFESLV::solve().
      //! Defaults reproduce the legacy driver harness numerics.
      //!@{
      int    MAX_ITER;       //!< max LM/Newton iterations (30)
      SolveWarmstart WARMSTART; //!< per-element warm-start strategy for a SOLVE.MARCHING solve (default REUSE)
      bool   VERBOSE;        //!< per-iteration LM/Newton trace to stdout (true)
      bool   REUSE;          //!< reuse last converged monolithic solution at the same input (default true)
      double LAMBDA0;        //!< initial LM damping (1e-8)
      double LAMBDA_DROP;    //!< lambda *= this on an accepted LM step (0.25)
      double LAMBDA_RAISE;   //!< lambda *= this on a failed factorization/step (10.0)
      double LAMBDA_MIN;     //!< lambda floor (1e-14)
      int    LS_MAX;         //!< backtracking line-search steps (24)
      double ARMIJO_C;       //!< Armijo sufficient-decrease constant (1e-4)
      bool   EQUILIBRATE;    //!< row/col inf-norm equilibration (true)
      //! @brief Derivative sparsity/coloring cache policy.
      //! DERIV_LAZY postpones FFDep sparsity propagation and coloring until a
      //! derivative/sparsity query needs them; DERIV_EAGER preserves the older
      //! setup-time behaviour.
      DerivCacheType DERIV_CACHE;
      //! @brief Linear residual cache policy.
      //! When enabled, purely linear rows are evaluated by a cached sparse
      //! affine map r = A z + b whenever the supplied constants match the
      //! setup constants.  Nonlinear rows and changing constants fall back to
      //! normal OCVar evaluation.
      LinearEvalCacheType LINEAR_EVAL_CACHE;
      //! @brief Solver-side verbosity, independent of the inherited DISPLAY_LEVEL: the model reports its plan,
      //! claims and classification, the solver its iterations, factorisation and marching.  The default, -1,
      //! follows the model's DISPLAY_LEVEL, resolved each time it is read (so a driver may set DISPLAY_LEVEL after
      //! construction); set it in code to make the two differ.  Not an environment knob.  MAXTHREAD is deliberately not
      //! duplicated -- it is a process-global RAII cap.
      int DISPLAY_LEVEL;
    } SOLVE;

    //! @brief The linear solve inside the forward/adjoint sensitivity computation.
    struct t_Sensitivity
    {
      bool   SPARSE;          //!< assemble the sensitivity Jacobian SPARSELY and solve with
                                   //!< arma::spsolve (SuperLU) instead of a dense arma::mat + arma::lu.
                                   //!<
                                   //!< WHY.  The dense Jacobian is ~2.6 GB at 18202^2 and O(n^3) to
                                   //!< factorise -- a scalability wall INDEPENDENT of the singularity
                                   //!< question, and one that affects IC_WEAK fixed-mesh models too.
                                   //!< It has simply not bitten yet because every fsens/asens driver in
                                   //!< the corpus is small.  The primal path has been sparse throughout.
                                   //!<
                                   //!< WHAT IT DOES NOT DO.  Sparse LU is still LU.  It does NOT make
                                   //!< the interface-augmented systems that fail today succeed -- it
                                   //!< reports the same singularity, though it fails cleanly (SuperLU
                                   //!< throws) instead of falling into a dense approximate solve that
                                   //!< does not terminate usefully.  That is better DIAGNOSABILITY, not
                                   //!< better capability.  Reliability is tier 1 (SPQR).
                                   //!<
                                   //!< This is TIER 0 of the eventual design, mirroring the primal:
                                   //!< SuperLU normally, rank-revealing QR when nTrace > 0.
                                   //!< Default ON; set false for the dense path.
    } SENSITIVITY;

    //! @brief Continuation in the parameters registered with add_homotopy().
    struct t_Homotopy
    {
      bool   ENABLE;           //!< honour homotopy parameters registered with add_homotopy() (true).
                                   //!< When false the schedule is IGNORED and solve() runs directly at the
                                   //!< current input values -- useful for with/without contrast runs
                                   //!< without unregistering the model's parameters.
      bool   SOLVE_S0;         //!< solve the s = 0 problem BEFORE the ramp, so every trial restarts from an actual
                               //!< solution rather than from the caller's initial vector (false: the ramp takes
                               //!< that vector as its start, as always).  A failure at s = 0 is reported as "the
                               //!< continuation has no starting point".  Costs one solve; on a model whose initial
                               //!< vector is already near the s = 0 solution it can shift the path (MBC5fig_ebal:
                               //!< one extra backtrack).  (false)
      double STEP_CAP;         //!< max homotopy step in s, and the step the ramp STARTS at (0.1).
                                   //!< PATH FIDELITY KNOB, and the one that matters: a converged step is
                                   //!< NOT evidence of correctness -- too large a cap can move the
                                   //!< corrector onto a different solution branch with no warning.
      double STEP_GROW;        //!< step growth after an accepted step (1.4).  Only acts as the
                                   //!< RECOVERY rate after a backtrack, since the ramp starts at the cap.
      double STEP_SHRINK;      //!< step shrink after a failed step (0.5).
                                   //!< *** REQUIRED: HOMOTOPY.STEP_GROW * HOMOTOPY.STEP_SHRINK < 1. ***  Otherwise a
                                   //!< step near the critical size alternates fail/succeed forever (shrink
                                   //!< then grow returns it above the failing value), wasting one solve per
                                   //!< accepted step.  The asymmetry -- shrink fast, grow slowly -- is what
                                   //!< makes the step ratchet down and settle.
      double STEP_MIN;         //!< abandon stage below this step (1e-5)
      bool   AUTO_DESTAGE;     //!< when a multi-parameter stage fails, split it into sequential
                                   //!< single-parameter sub-stages from the point it reached, rather than
                                   //!< giving up (true).  Announced LOUDLY on stderr and recorded in
                                   //!< StageReport::destaged: needing it usually says something about the
                                   //!< MODEL -- those parameters interact too strongly to move together --
                                   //!< so it should be noticed, not silently absorbed.
      bool   PREDICT;          //!< guarded extrapolation predictor (true)
      int    PREDICT_ORDER;    //!< max order of the extrapolation from accepted points (2).
                                   //!< 0 = plain reuse, 1 = secant (2 points), 2 = quadratic (3 points).
                                   //!< Every candidate up to this order is TRIED and the smallest-residual
                                   //!< one is kept, so raising it can only help: a worse extrapolation is
                                   //!< discarded for the price of one residual evaluation.
                                   //!< NOTE a TANGENT predictor is deliberately ABSENT.  It is redundant
                                   //!< with the corrector's first Newton step: at a converged x(s) with the
                                   //!< parameters moved to s+ds, r(x,p(s+ds)) ~ ds*dr/ds, so the first
                                   //!< Newton step is -J^-1 r ~ ds*dx/ds -- the tangent itself.  Computing
                                   //!< it explicitly would cost a Jacobian and a factorization to save the
                                   //!< very iteration that produces it.  Extrapolation from already-accepted
                                   //!< points needs neither.
    } HOMOTOPY;


    Options
      ()
      { reset(); }

    Options& operator=
      ( Options const& opt )
      {
        OCFEBASE::Options::operator=( opt );
        SOLVE.DISPLAY_LEVEL = opt.SOLVE.DISPLAY_LEVEL;
        SOLVE.DERIV_CACHE         = opt.SOLVE.DERIV_CACHE;
        SOLVE.LINEAR_EVAL_CACHE   = opt.SOLVE.LINEAR_EVAL_CACHE;
        SOLVE.MAX_ITER      = opt.SOLVE.MAX_ITER;
        SOLVE.RES_TOL       = opt.SOLVE.RES_TOL;
        SOLVE.LAMBDA0       = opt.SOLVE.LAMBDA0;
        SOLVE.LAMBDA_DROP   = opt.SOLVE.LAMBDA_DROP;
        SOLVE.LAMBDA_RAISE  = opt.SOLVE.LAMBDA_RAISE;
        SOLVE.LAMBDA_MIN    = opt.SOLVE.LAMBDA_MIN;
        SOLVE.LS_MAX        = opt.SOLVE.LS_MAX;
        SOLVE.ARMIJO_C      = opt.SOLVE.ARMIJO_C;
        SOLVE.EQUILIBRATE   = opt.SOLVE.EQUILIBRATE;
        SOLVE.FACTORIZATION = opt.SOLVE.FACTORIZATION;
        SOLVE.VERBOSE       = opt.SOLVE.VERBOSE;
        _force_verbose_from_env( SOLVE.VERBOSE );   // rev115: AFTER the driver, not only before
        SENSITIVITY.SPARSE      = opt.SENSITIVITY.SPARSE;
        HOMOTOPY.ENABLE          = opt.HOMOTOPY.ENABLE;
        HOMOTOPY.STEP_CAP        = opt.HOMOTOPY.STEP_CAP;
        HOMOTOPY.STEP_GROW       = opt.HOMOTOPY.STEP_GROW;
        HOMOTOPY.STEP_SHRINK     = opt.HOMOTOPY.STEP_SHRINK;
        HOMOTOPY.STEP_MIN        = opt.HOMOTOPY.STEP_MIN;
        HOMOTOPY.AUTO_DESTAGE    = opt.HOMOTOPY.AUTO_DESTAGE;
        HOMOTOPY.SOLVE_S0        = opt.HOMOTOPY.SOLVE_S0;
        HOMOTOPY.PREDICT         = opt.HOMOTOPY.PREDICT;
        HOMOTOPY.PREDICT_ORDER   = opt.HOMOTOPY.PREDICT_ORDER;
        SOLVE.MARCHING      = opt.SOLVE.MARCHING;
        SOLVE.REUSE         = opt.SOLVE.REUSE;
        SOLVE.WARMSTART     = opt.SOLVE.WARMSTART;
        return *this;
      }

    void apply_environment_defaults
      ()
      {
        OCFEBASE::Options::apply_environment_defaults();
        // 2026-10-07 (WORKPLAN 3.A): the environment no longer sets an option that changes what is computed or decided --
        // CRONOS_DETERMINACY_FATAL, _TOL, _BACKEND, CRONOS_SOLVE_FACTORIZATION, CRONOS_TRACE_PROJ_EPS, CRONOS_SAT_SIGMA0
        // (and FFModel's CRONOS_CLASSIFY_ROBUST) are retired: the option is the interface.  The print-only switches above
        // (audit reports, forced verbosity) remain, as diagnostics that instrument a whole sweep without editing drivers.
        SOLVE.VERBOSE       = true;
      }

    void reset
      ()
      {
        OCFEBASE::Options::reset();
        SOLVE.DERIV_CACHE         = DERIV_LAZY;
        SOLVE.LINEAR_EVAL_CACHE   = LINEAR_NONE;
        SOLVE.MAX_ITER      = 30;
        SOLVE.RES_TOL       = 1e-9;
        SOLVE.LAMBDA0       = 1e-8;
        SOLVE.LAMBDA_DROP   = 0.25;
        SOLVE.LAMBDA_RAISE  = 10.0;
        SOLVE.LAMBDA_MIN    = 1e-14;
        SOLVE.LS_MAX        = 24;
        SOLVE.ARMIJO_C      = 1e-4;
        SOLVE.EQUILIBRATE   = true;
        // rev113: run-time selection of the linear backend.
        //
        // SOLVE_SPQR is the sparse rank-revealing QR min-norm solve on Js DIRECTLY -- no
        // Js^T Js, so no condition squaring -- and item 13 already skips the doomed tier-0
        // SuperLU factorization a priori when _nTraceVar > 0.  That is precisely the
        // configuration a gauge-singular interface-augmented system wants, and until now
        // it could only be selected by editing a driver.
        //
        // WHY IT MATTERS (measured, OCFE_DAE5 --nomarch, IC_STRONG):
        //   the kernel is 100% multiplier-block, so the primal is DETERMINED and the answer
        //   is correct (2.56e-07 / 2.07e-08 / 1.51e-09, the most accurate of every
        //   configuration tried).  What fails is the CONVERGENCE TEST: Js is numerically
        //   singular, solve_linear refuses, nt_dead latches PERMANENTLY, tier 1 forms the
        //   damped normal equations and squares cond 2.5e16 into ~6e32, and LM stalls.
        //   Perturbing the singular directions above machine precision fixes convergence
        //   (IC_TRACE eps=1e-14 converges; eps=1e-20 and eps=0 do not) but costs 17%
        //   accuracy: 2.56e-07 -> 2.99e-07.  A rank-tolerant solve should need no
        //   perturbation and keep the better answer.  UNVERIFIED -- this hook exists to
        //   test it.
        SOLVE.FACTORIZATION = SOLVE_SUPERLU;

        // rev114: let the environment FORCE solver verbosity on.
        //
        // The [step-fail] attribution -- step energy by state, dominant components, merit
        // vs alpha -- is the instrument that says WHICH STATES carry a failing step.  It is
        // gated on SOLVE.VERBOSE, which several drivers set false (MMPDE27, MBC5, ...), so
        // on exactly the models where a step failure is hard to explain the attribution is
        // unavailable without editing the driver.
        //
        // Measured need: OCFE_MBC5 at CONSTANT dP reports "LM/Newton step failed at
        // iteration 28" and nothing else, while the z-varying case does not fail at all.
        // The difference is that constant dP makes omega(z) constant -- omega[min,max] =
        // [0.0405,0.0405] against [0.0405,0.1596] -- so the hypothesis is a degeneracy from
        // the lost z-dependence.  That is a statement about WHICH STATES go singular, and
        // [step-fail] is the only thing that reports it.
        //
        // Force-on only.  A driver that sets SOLVE.VERBOSE=false does so to keep its own
        // output readable; the variable overrides that for a diagnostic run and never the
        // other way, so no driver loses output it expects.
        _force_verbose_from_env( SOLVE.VERBOSE );
        // rev42: warn threshold for the automatic duplicate-node spread report.  Exact
        // imposition is round-off limited at any mesh; weak is truncation limited.  1e-10 sits
        // far above the former and far below the latter, so it separates them without firing
        // on IC_WEAK.  Report-only -- see the declaration.
        // rev49: REFUTED, kept only so the dead end is documented at the point of temptation.
        // Measured on MMPDE10 nel_xi=8 IC_TRACE, delta-sweep at fixed h_rel:
        //     delta_rel 1e-5 -> analytic 4.14230536 (ratio 1.2434 vs FD)
        //     delta_rel 1e-6 -> analytic 3.49817529 (ratio 1.0501 vs FD)
        //     delta_rel 1e-7 -> SINGULAR (rcond 3.06e-17), no result
        // The answer drifts strongly with delta and would need roughly delta_rel 1e-7 to reach
        // the FD value 3.3314 -- which is exactly where the factorization fails.  THE WINDOW IS
        // EMPTY.  Worse, rcond is NON-MONOTONE in delta (3.06e-17 at 2.165e-03 but 1.43e-16 at
        // 2.165e-04), which a shift that genuinely lifted sigma_min could not produce: the null
        // directions are COMBINATIONS spanning tau columns, not aligned with individual tau
        // coordinates, so a diagonal shift perturbs the spectrum without removing the kernel.
        // Do not re-enable this expecting a usable delta.  Use SENS_MINNORM_PINV / SPQR.
        SENSITIVITY.SPARSE     = true;     // rev50: ON by default -- the dense path is a memory wall
        // Stage 3: extend the pointwise-pin suppressor to states pinned through OTHER states
        // pinned at the same node.  Default ON: unlike Stage A this rests on a row implication
        // by physical rows, the same footing as FIX 1.  CRONOS_NO_PIN_CLOSURE=1 disables it
        // for before/after comparison without a rebuild.
        HOMOTOPY.ENABLE          = true;   //!< honour registered homotopy parameters inside solve()
        HOMOTOPY.STEP_CAP        = 0.1;    //!< max homotopy step in s; also the STARTING step (path fidelity: larger steps can SKIP BRANCHES)
        HOMOTOPY.STEP_GROW       = 1.4;    //!< step growth after an accepted step (recovery rate after a backtrack)
        HOMOTOPY.STEP_SHRINK     = 0.5;    //!< step shrink after a failed step; GROW*SHRINK MUST BE < 1
        HOMOTOPY.STEP_MIN        = 1e-5;   //!< abandon the stage below this step
        HOMOTOPY.AUTO_DESTAGE    = true;   //!< on stage failure, retry its parameters one at a time (LOUD)
        HOMOTOPY.SOLVE_S0        = false;  //!< the ramp starts from the caller's vector (2026-10-10: option, 4R.18)
        HOMOTOPY.PREDICT         = true;   //!< guarded extrapolation predictor (kept only when it beats plain reuse)
        HOMOTOPY.PREDICT_ORDER   = 2;  //!< max order of the extrapolation from accepted points
        SOLVE.MARCHING      = true;   //!< auto-march when AVAILABLE (evolution domain + >=1 differential state); else monolithic
        SOLVE.REUSE         = true;   //!< reuse the last converged monolithic solution when solved again at the same input (skips the redundant primal Newton across value+gradient at one iterate)
        SOLVE.WARMSTART     = REUSE;  //!< per-element warm-start for a SOLVE.MARCHING solve
        SOLVE.DISPLAY_LEVEL = -1;   // rev314: follow DISPLAY_LEVEL -- see _solve_display_level()
        apply_environment_defaults();
      }
  };

  Options options;

  //! @brief OCFEBASE::_solve_display_level(): SOLVE.DISPLAY_LEVEL when set, else the model's.
  int _solve_display_level() const override
    { return options.SOLVE.DISPLAY_LEVEL < 0 ? options.DISPLAY_LEVEL : options.SOLVE.DISPLAY_LEVEL; }
  //! @brief OCFEBASE::_march_requested(): SOLVE.MARCHING.
  bool _march_requested() const override { return options.SOLVE.MARCHING; }

  //! @brief The collocation grid of direction @p dir (OCBase::colloc_dom): the grid the equations are collocated on,
  //! onto which a nonlinear operation lifts a coarser multi-node input first (2026-09-28)
  FFDom const* colloc_dom
    ( FFVar const& dir )
    const override
    { auto const it = _mDom.find( dir );  return it == _mDom.end()? nullptr: &it->second; }


  //! @brief OCFESLV imposes transitions through lifted auxiliaries (FFModel::_transitions_lifted): the continuity claim
  //! of a post-jump state across tau joins x_i(tau^+) to its auxiliary w_i instead of to x_i(tau^-) -- see
  //! _trace_claim_coeffs -- in every imposition mode (penalty under IC_WEAK, exact under IC_STRONG / IC_TRACE).
  bool _transitions_supported() const override { return true; }
  bool _transitions_lifted() const override { return !options.SOLVE.MARCHING; }   // marching: a transfer map (_apply_march_jumps)
  //! @brief v1: tau on an element boundary of the evolution direction -- in monolithic AND marching solves (marching
  //! validates against the user's DECLARED grid, as its working domain is collapsed to one window: a tau between two
  //! seams was silently ignored before 2026-10-06, WORKPLAN 1.6)
  bool _validate_transition_tau_unlifted() const override { return true; }
  bool _accept_transition_tau
    ( double const tau, std::string& why )
    const override
    {
      FFDom const* pd = nullptr;
      if( _transitions_lifted() ){
        auto const it = _mDom.find( _evolution_dom_var );  if( it != _mDom.end() ) pd = &it->second;
      }
      else{
        auto const it = _usr._mDomUsr.find( _usr._evolution_dom_varUsr );  if( it != _usr._mDomUsr.end() ) pd = &it->second;
      }
      if( !pd ){ why = "no evolution domain"; return false; }
      for( size_t k = 1; k < pd->n_elem; ++k )
        if( std::fabs( pd->elem_bnd[k] - tau ) <= 64.*DBL_EPSILON*std::max( 1., std::fabs( tau ) ) ) return true;
      why = "tau must be an ELEMENT BOUNDARY of the evolution direction (add one there)";  return false;
    }


  //! @brief Build the marching transfer maps (see t_MarchJump); false with the reason printed if v1 cannot march them
  bool _build_march_jumps
    ()
    {
      _marchJumps.clear();      // frees the Jacobians (unique_ptr)
      _marchJumpsReady = true;
      if( _mTrn.empty() ) return true;
      FFModel::t_TransitionJumps G;  std::string why;
      auto is_state = [&]( FFVar const& v ){ return _mVar.count( v ) > 0; };
      if( !_transition_jumps( _mTrn, is_state, G, why ) ){ _disp( 2 ) << "OCFESLV::solve ** TRANSITION: " << why << std::endl; return false; }
      FFVar const t = _evolution_dom_var;
      for( auto const& g : G ){
        t_MarchJump M;  M.tau = g.tau;  M.xs = g.jump;
        for( auto const& [x, dm] : _mVar ) if( dm.count( t ) ) M.evo.push_back( x );   // v1.5: profiles, node by node
        std::vector<FFVar> tgt, rep;
        for( auto const& x : M.xs ){ std::ostringstream nm; nm << "Xp_" << x.name() << "@" << g.tau; M.xp.push_back( _dag->add_var( nm.str() ) ); tgt.push_back( x ); rep.push_back( M.xp.back() ); }
        for( size_t const k : g.comp ){
          auto const& tr = _mTrn[ k / 65536 ];  size_t const c = k % 65536;
          M.R.push_back( tr.left[c] - _dag->compose( std::vector<FFVar>{ tr.right[c] }, tgt, rep )[0] );
        }
        FFSubgraph sg = _dag->subgraph( M.R.size(), M.R.data() );
        for( auto const* op : sg.l_op ){
          if( !op || op->type != FFOp::VAR || !op->varout[0] ) continue;
          FFVar const& v = *op->varout[0];
          if( _mDom.count( v ) ){ _disp( 2 ) << "OCFESLV::solve ** TRANSITION at tau = " << g.tau << ": a marching map may not read the"
                                              " domain variable " << v.name() << " (v1.5: its node coordinates are not available at the seam)" << std::endl; return false; }
          auto const ii = _mInp.find( v );
          if( ii == _mInp.end() ) continue;
          if( !ii->second.empty() ){ _disp( 2 ) << "OCFESLV::solve ** TRANSITION at tau = " << g.tau << ": marching v1 transfers maps whose inputs"
                                                  " are time-invariant (" << v.name() << " is distributed)" << std::endl; return false; }
          M.inp.push_back( v );
        }
        M.J .reset( _dag->FAD( M.R.size(), M.R.data(), M.xp.size(),  M.xp.data() ) );
        M.Jx.reset( _dag->FAD( M.R.size(), M.R.data(), M.evo.size(), M.evo.data() ) );
        if( !M.inp.empty() ) M.Jp.reset( _dag->FAD( M.R.size(), M.R.data(), M.inp.size(), M.inp.data() ) );
        _marchJumps.push_back( std::move( M ) );
      }
      return true;
    }

  //! @brief At window k's seam tau = _marchGrid[k]: if transitions sit there, solve R( x-, x+, p ) = 0 for x+ (Newton;
  //! one step for an explicit map) from the previous window's terminals x-, and hand x+ on instead of the terminal --
  //! into the IC input (data path) or _marchICTerminal (general-IC override path).  False if Newton fails.
  bool _apply_march_jumps
    ( size_t const k, double const* var, double* pinp, double const* cst )
    {
      if( !_marchJumpsReady && !_build_march_jumps() ) return false;
      double const tau = _marchGrid[k];
      for( auto& M : _marchJumps ){
        if( std::fabs( M.tau - tau ) > 64.*DBL_EPSILON*std::max( 1., std::fabs( tau ) ) ) continue;
        size_t const m = M.xs.size(), ne = M.evo.size();
        // the previous window's terminals, node by node (a scalar state is broadcast)
        auto const prof = terminal_profile( var, pinp, cst );
        std::vector< std::vector<double> > xe( ne );  size_t nn = 1;
        for( size_t e = 0; e < ne; ++e ){
          auto const ip = prof.find( M.evo[e] );
          xe[e] = ( ip != prof.end() && !ip->second.empty() )? ip->second: std::vector<double>{ terminal_value( M.evo[e], t_Coord{}, var, pinp, cst ) };
          if( xe[e].size() > 1 ){
            if( nn > 1 && xe[e].size() != nn ){ _disp( 2 ) << "OCFESLV::solve ** TRANSITION at tau = " << tau << ": profiles with different"
                                                           " node counts in one marching map (v1.5)" << std::endl; return false; }
            nn = xe[e].size(); }
        }
        std::vector<FFVar> vv;  for( auto const& x : M.evo ) vv.push_back( x );  for( auto const& u : M.xp ) vv.push_back( u );
        std::vector<double> tail;
        for( auto const& u : M.inp ){ vv.push_back( u ); tail.push_back( pinp[ pos_input( u, std::map<FFVar,size_t,lt_FFVar>() ) ] ); }
        for( size_t c = 0; c < _vCst.size(); ++c ){ vv.push_back( _vCst[c] ); tail.push_back( cst? cst[c]: 0. ); }
        M.vars = vv;  M.at_n.assign( nn, std::vector<double>() );
        std::vector< std::vector<double> > wprof( m, std::vector<double>( nn ) );
        for( size_t j = 0; j < nn; ++j ){
          std::vector<double> vl( ne + m );
          for( size_t e = 0; e < ne; ++e ) vl[e] = xe[e].size() > 1? xe[e][j]: xe[e][0];
          std::vector<double> xp( m );
          for( size_t i = 0; i < m; ++i ) for( size_t e = 0; e < ne; ++e ) if( M.evo[e].id() == M.xs[i].id() ) xp[i] = vl[e];   // guess: x-
          vl.insert( vl.end(), tail.begin(), tail.end() );
          bool conv = false;
          for( int it = 0; it < 30 && !conv; ++it ){
            for( size_t i = 0; i < m; ++i ) vl[ne+i] = xp[i];
            std::vector<double> r( m ), a( m*m );
            _dag->eval( m, M.R.data(), r.data(), vv.size(), vv.data(), vl.data() );
            _dag->eval( m*m, M.J.get(), a.data(), vv.size(), vv.data(), vl.data() );
            double rn = 0.;  for( double const q : r ) rn = std::max( rn, std::fabs( q ) );
            double xn = 1.;  for( double const q : xp ) xn = std::max( xn, std::fabs( q ) );
            if( rn <= 1e-14 * xn ){ conv = true; break; }
            std::vector<double> b( m );  for( size_t i = 0; i < m; ++i ) b[i] = -r[i];
            if( !_small_solve( a, b, m ) ){ _disp( 2 ) << "OCFESLV::solve ** TRANSITION at tau = " << tau << ": singular dR/dx+" << std::endl; return false; }
            for( size_t i = 0; i < m; ++i ) xp[i] += b[i];
          }
          if( !conv ){ _disp( 2 ) << "OCFESLV::solve ** TRANSITION at tau = " << tau << ": the jump map did not converge" << std::endl; return false; }
          for( size_t i = 0; i < m; ++i ){ vl[ne+i] = xp[i];  wprof[i][j] = xp[i]; }
          M.at_n[j] = vl;
        }
        M.at = M.at_n[0];
        for( size_t i = 0; i < m; ++i ){                 // hand x+ on through the state's own transfer path
          FFVar const& st = M.xs[i];
          for( auto const& [s2, ic] : _marchTransfers ){
            if( s2.id() != st.id() ) continue;
            if( ic.dag() ){
              // the canonical IC-input path: terminal_profile() and set_input_values() share the element-major /
              // node-inner layout, so a post-jump profile on the SAME spatial grid as the input is written directly
              if( nn > 1 ){
                size_t ndof = 1;  auto const ii = _mInp.find( ic );
                if( ii != _mInp.end() ) for( auto const& d : ii->second ){ auto const& od = _mDom.at( d ); ndof *= od.n_elem * od.n_node; }
                if( ii == _mInp.end() || ndof != nn ){
                  _disp( 2 ) << "OCFESLV::solve ** TRANSITION at tau = " << tau << ": the IC input of " << st.name() << " is not on the"
                               " state's spatial grid (" << ndof << " vs " << nn << " nodes) -- v1.5 does not resample" << std::endl; return false; }
                set_input_values( ic, wprof[i], pinp );
              }
              else set_input_values( ic, std::vector<double>{ wprof[i][0] }, pinp );
            }
            else{
              auto itt = _marchICTerminal.find( st );
              if( itt != _marchICTerminal.end() ){
                if( nn > 1 && itt->second.size() == nn ) itt->second = wprof[i];
                else std::fill( itt->second.begin(), itt->second.end(), wprof[i][0] );
              }
            }
          }
        }
      }
      return true;
    }

  //! @brief Solve a x = b in place for a small dense m x m system (partial pivoting); false if singular
  static bool _small_solve
    ( std::vector<double> a, std::vector<double>& b, size_t const m, bool const transpose = false )
    {
      if( transpose ) for( size_t i = 0; i < m; ++i ) for( size_t j = i+1; j < m; ++j ) std::swap( a[i*m+j], a[j*m+i] );
      for( size_t c = 0; c < m; ++c ){
        size_t pr = c;  for( size_t i = c+1; i < m; ++i ) if( std::fabs( a[i*m+c] ) > std::fabs( a[pr*m+c] ) ) pr = i;
        if( std::fabs( a[pr*m+c] ) < 1e-300 ) return false;
        if( pr != c ){ for( size_t j = 0; j < m; ++j ) std::swap( a[c*m+j], a[pr*m+j] ); std::swap( b[c], b[pr] ); }
        for( size_t i = c+1; i < m; ++i ){ double const f = a[i*m+c]/a[c*m+c]; for( size_t j = c; j < m; ++j ) a[i*m+j] -= f*a[c*m+j]; b[i] -= f*b[c]; }
      }
      for( size_t c = m; c-- > 0; ){ double q = b[c]; for( size_t j = c+1; j < m; ++j ) q -= a[c*m+j]*b[j]; b[c] = q/a[c*m+c]; }
      return true;
    }

  //! @brief ADJOINT through the jump at window k's seam -- the transpose of _march_jump_fsens.  With y = J^{-T} mu+ (mu+
  //! the adjoint on the post-jump initial values, from @p lam_term), the costate handed to window k-1 is
  //!   mu-_e = [e continuous] mu+_e - ( dR/dx-^T y )_e,   and the input-space gradient gains  - dR/dinp^T y.
  bool _march_jump_asens
    ( size_t const k, std::map<FFVar,std::vector<double>,lt_FFVar>& lam_term, std::vector<double>& ldinp )
    {
      double const tau = _marchGrid[k];
      for( auto const& M : _marchJumps ){
        if( std::fabs( M.tau - tau ) > 64.*DBL_EPSILON*std::max( 1., std::fabs( tau ) ) || M.at_n.empty() ) continue;
        size_t const m = M.xs.size(), ne = M.evo.size(), ni = M.inp.size(), nn = M.at_n.size();
        std::map<FFVar,std::vector<double>,lt_FFVar> incoming = lam_term;   // mu+ as handed back from window k
        for( size_t i = 0; i < m; ++i ){ auto& v = lam_term[ M.xs[i] ]; if( v.empty() ) v.assign( nn, 0. ); else std::fill( v.begin(), v.end(), 0. ); }
        for( size_t e = 0; e < ne; ++e ){ auto& v = lam_term[ M.evo[e] ]; if( v.empty() ) v.assign( nn, 0. ); }
        for( size_t j = 0; j < nn; ++j ){                // node by node (v1.5 profiles; nn = 1 for scalars)
          std::vector<double> a( m*m ), ax( m*ne ), ap( m*ni );
          _dag->eval( m*m,  M.J.get(),  a.data(),  M.vars.size(), M.vars.data(), M.at_n[j].data() );
          _dag->eval( m*ne, M.Jx.get(), ax.data(), M.vars.size(), M.vars.data(), M.at_n[j].data() );
          if( ni ) _dag->eval( m*ni, M.Jp.get(), ap.data(), M.vars.size(), M.vars.data(), M.at_n[j].data() );
          std::vector<double> y( m, 0. );
          for( size_t i = 0; i < m; ++i ){ auto it = incoming.find( M.xs[i] );
            if( it != incoming.end() && !it->second.empty() ) y[i] = it->second.size() > j? it->second[j]: it->second[0]; }
          if( !_small_solve( a, y, m, /*transpose=*/true ) ) return false;
          for( size_t e = 0; e < ne; ++e ){
            double q = 0.;  for( size_t i = 0; i < m; ++i ) q -= ax[i*ne+e]*y[i];
            auto& v = lam_term[ M.evo[e] ];
            if( v.size() > j ) v[j] += q; else for( auto& z : v ) z += q;
          }
          for( size_t qi = 0; qi < ni; ++qi ){
            double q = 0.;  for( size_t i = 0; i < m; ++i ) q -= ap[i*ni+qi]*y[i];
            ldinp[ pos_input( M.inp[qi], std::map<FFVar,size_t,lt_FFVar>() ) ] += q;
          }
        }
      }
      return true;
    }

  //! @brief FORWARD sensitivity through the jump at window k's seam: for the post-jump states,
  //!   s+ = -J^{-1} ( dR/dx- . s- + dR/dinp . dinp )      (implicit function theorem; J = dR/dx+)
  //! replaces the terminal sensitivity carried in @p sens_term.  @p dinp: the direction, input layout.
  bool _march_jump_fsens
    ( size_t const k, std::map<FFVar,std::vector<double>,lt_FFVar>& sens_term, std::vector<double> const& dinp )
    {
      double const tau = _marchGrid[k];
      for( auto const& M : _marchJumps ){
        if( std::fabs( M.tau - tau ) > 64.*DBL_EPSILON*std::max( 1., std::fabs( tau ) ) || M.at_n.empty() ) continue;
        size_t const m = M.xs.size(), ne = M.evo.size(), ni = M.inp.size(), nn = M.at_n.size();
        std::vector< std::vector<double> > out( m, std::vector<double>( nn ) );
        for( size_t j = 0; j < nn; ++j ){                // node by node (v1.5 profiles; nn = 1 for scalars)
          std::vector<double> a( m*m ), ax( m*ne ), ap( m*ni );
          _dag->eval( m*m,  M.J.get(),  a.data(),  M.vars.size(), M.vars.data(), M.at_n[j].data() );
          _dag->eval( m*ne, M.Jx.get(), ax.data(), M.vars.size(), M.vars.data(), M.at_n[j].data() );
          if( ni ) _dag->eval( m*ni, M.Jp.get(), ap.data(), M.vars.size(), M.vars.data(), M.at_n[j].data() );
          std::vector<double> b( m, 0. );
          for( size_t i = 0; i < m; ++i ){
            for( size_t e = 0; e < ne; ++e ){ auto it = sens_term.find( M.evo[e] ); double se = 0.;
              if( it != sens_term.end() && !it->second.empty() ) se = it->second.size() > j? it->second[j]: it->second[0];
              b[i] -= ax[i*ne+e]*se; }
            for( size_t q = 0; q < ni; ++q ) b[i] -= ap[i*ni+q]*dinp[ pos_input( M.inp[q], std::map<FFVar,size_t,lt_FFVar>() ) ];
          }
          if( !_small_solve( a, b, m ) ) return false;
          for( size_t i = 0; i < m; ++i ) out[i][j] = b[i];
        }
        for( size_t i = 0; i < m; ++i ){ auto& v = sens_term[ M.xs[i] ];
          if( v.size() == nn ) v = out[i]; else if( v.empty() ) v = out[i]; else std::fill( v.begin(), v.end(), out[i][0] ); }
      }
      return true;
    }


  //! @brief W-test the claims dropped by the redundancy pass on the built system, and restore as explicit the subset
  //! that is needed (one re-derive, revert-guarded).  Default on; (The CRONOS_DROP_WDECIDE variable is retired, 2026-10-07.)
  static bool _knob_DROP_WDECIDE
    ()
    { return true; }


  //! @brief Tier-0 skip policy on interface-augmented systems.  ALWAYS (default): skip plain Newton whenever the
  //! system carries trace variables under SPQR.  AUTO consults the setup audit's right-deficiency for the window and
  //! attempts tier-0 only on an affirmative full-rank verdict -- no verdict, or a deficiency, keeps the skip.  NEVER attempts
  //! tier-0 regardless; a dead tier-0 falls through to the QR tier, so a wrong NEVER costs one failed factorization, not
  //! correctness.  Environment: CRONOS_SKIP_NT_ON_TRACE = always|auto|never (or 0|1|2), applied at its single use site.
  static int _knob_SKIP_NT_ON_TRACE
    ()
    { return Options::NT_SKIP_ALWAYS; }

  //! @brief Diagonal regularisation added to the TAU BLOCK ONLY of the
  //! sensitivity Jacobian before factorisation, as a MULTIPLE of
  //! ||J||_inf.  0 (default) = off.
  //! WHY.  An interface-augmented Jacobian is STRUCTURALLY SINGULAR
  //! (see the note at the SPQR tier-0 skip): the primal never
  //! factorises it, it uses the QR min-norm tier.  The sensitivity
  //! path uses a dense LU, which reports rcond ~1e-17 and falls back
  //! to an approximate solve that does not terminate usefully --
  //! measured on MMPDE10 at nel_xi=8, IC_TRACE (18202^2) while
  //! IC_WEAK at the same mesh (15360^2) is clean.  The 2842 tau
  //! columns are the whole difference.
  //! WHAT IT ASSUMES.  Orphan multiplier columns are a GAUGE
  //! freedom: null(J) lies in the tau block and dF/dvar is zero
  //! there, so the tau components of the sensitivity solution do not
  //! enter dF/dp.  Pinning them to zero IS the min-norm choice, and
  //! dF/dp is unchanged to O(delta).
  //! THE ASSUMPTION IS NOT PROVEN.  It has two parts and they are
  //! different claims:
  //! (1) dF/dvar _|_ null(J)      -> the gradient is well defined;
  //! (2) null(J) subset tau span  -> pinning s_tau removes ALL of
  //! it, so J_reg is actually nonsingular.
  //! (1) follows from the physical solution being unique under
  //! IC_TRACE (its tau-table matches IC_WEAK at every tau).  (2)
  //! does NOT follow from (1), and the fixed ~2.4e-09 near-null
  //! direction recorded in the roadmap is unexplained and may have a
  //! physical component.  _knob_AUDIT_SENS_NULL (CRONOS_AUDIT_SENS)
  //! measures (1); a sweep of this knob's delta measures (2).
  static double _knob_SENS_TAU_REG
    ()
    { return 0.0; }

  //! @brief Solve the sensitivity systems with a dense pseudo-inverse
  //! (arma::pinv) instead of LU.  DIAGNOSTIC ONLY -- O(n^3) dense SVD,
  //! usable only at nel small enough to fit (a few thousand columns).
  //! WHY IT EXISTS.  The interface-augmented Jacobian is structurally
  //! singular, so J s = b determines s only up to null(J).  The
  //! reduced gradient is nonetheless WELL DEFINED -- measured: with
  //! AUDIT_SENS_NULL, dF/dvar has ZERO nonzero entries on the tau
  //! columns, so dF/dvar _|_ null(J) and every solution s gives the
  //! same dF/dvar . s.  What is needed is a solver that returns ANY
  //! valid solution rather than failing; the min-norm one is the
  //! natural choice and is what the primal's QR tier already returns.
  //! WHAT IT IS FOR.  pinv is the cheap STAND-IN for SPQR: if the
  //! min-norm solution reproduces the finite-difference gradient at a
  //! small mesh, then adopting SPQR in the sensitivity path will
  //! reproduce it at scale -- and that is worth knowing BEFORE paying
  //! for the factorization-ownership work SPQR needs (a handle in
  //! FSensBlock, destructor/move semantics, and a SECOND
  //! factorization for J^T since QR does not transpose the way LU
  //! does).
  //! NOT A PRODUCTION PATH.  Dense SVD of an 18202^2 system is ~2.6 GB
  //! and O(n^3); use nel_xi 2-3.
  static bool _knob_SENS_MINNORM_PINV
    ()
    { return false; }


  //! return false from setup() immediately after the audit reports,
  //! without
  //! building the derivative/linear caches or running any solve.
  //! For structural corpus sweeps: every audited quantity is fixed
  //! at setup time.  Environment: CRONOS_AUDIT_SETUP_ONLY.
  static bool _knob_AUDIT_SETUP_ONLY
    ()
    { return Options::_env_flag( "CRONOS_AUDIT_SETUP_ONLY", false ); }


  //! @brief Report ||dF/dvar|| restricted to the tau columns, once per sensitivity setup -- claim (1) of
  //! _knob_SENS_TAU_REG, and cheap: dF/dvar is already assembled as B.gfct/B.cfct, so no factorisation is involved.
  //! Off by default; CRONOS_AUDIT_SENS.
  static bool _knob_AUDIT_SENS_NULL
    ()
    { return Options::_env_flag( "CRONOS_AUDIT_SENS", false ); }

  //! after every solve(), report per-state duplicate-node
  //! spreads (report_duplicate_node_spreads).  Nodes that coincide
  //! geometrically but belong to different elements carry separate
  //! unknowns, so the spread over each coincidence group measures
  //! the interface conditions ACTUALLY satisfied by the solution --
  //! independently of rank, residual, verdict or convergence, all
  //! four of which were clean when Stage A destroyed MBC5's radial
  //! fluxes by 24 orders.  Off by default; CRONOS_DUP_SPREAD=1.
  //! NOT under the CRONOS_AUDIT umbrella: it changes driver output.
  static bool _knob_AUDIT_DUP_SPREAD
    ()
    { return Options::_env_flag( "CRONOS_DUP_SPREAD", false ); }

  //! warn (non-fatally) when the worst max_pair spread
  //! exceeds this under an EXACT imposition mode, where continuity
  //! is round-off limited (~1e-17 at any mesh) rather than
  //! truncation limited.  Default 1e-10.  REPORT ONLY -- promoting
  //! it to a failure is a correctness gate and needs an over-fire
  //! scan across the corpus first; see the detect-then-promote
  //! discipline that governed every other gate in this arc.
  static double _knob_DUP_SPREAD_WARN
    ()
    { return Options::_env_dbl( "CRONOS_DUP_SPREAD_WARN", 1e-10 ); }


  //! @brief the model's rules, then the solver's on the solver's own fields.
  void _normalise_equation_options
    ( FFModel::EqnOptions& opt, std::map<FFVar,int,lt_FFVar> const& eqndom )
    const override
    {
      FFModel::_normalise_equation_options( opt, eqndom );
      auto* p = dynamic_cast<EqnOptions*>( &opt );
      if( !p ) throw std::logic_error( "OCFESLV::_normalise_equation_options ** not the solver's options type" );
      if( p->role == EqnRole::DIAGNOSTIC ){ p->receive_sat = false; p->donate_for_state_continuity = false; }
      if( p->role == EqnRole::LINK ) p->donate_for_state_continuity = false;
    }

  //! @brief slice the requested configuration into the applied one, and host the imposition
  //! instrument that used to live inside the model.
  void _apply_options
    ()
    override
    {
      // 2026-10-12: the CRONOS_DISPLAY_LEVEL floor is applied to the SOLVER's object first, so the slice below
      // carries it into FFModel::options too and both layers report at the same level.
      FFModel::Options::_force_display_from_env( options.DISPLAY_LEVEL );
      // slice the requested configuration into the APPLIED one (FFModel::options is the member the
      // model layer reads; `options` here is the solver's richer, shadowing object).
      FFModel::options = static_cast<FFModel::Options const&>( options );

      // IC_TRACE system instead -- the instrument for the MBC-family TRACE defect
      // analysis (pre-existing upstream: MBC3/MBC5 fail under IC_TRACE while
      // legacy IC_STRONG passes).  Off by default: rev144 behaviour untouched.
      // (CRONOS_STRONG_VIA_TRACE (retired 2026-10-07, WORKPLAN 3.B batch 2b): a driver wanting the trace system sets IC_TRACE)

    }


  //! @brief Physical equation-node metadata after interface filtering.
  //!
  //! `retained` is true when the physical equation row is present in the
  //! residual system.  Under IC_STRONG, `interface_donor` marks rows removed
  //! because exact interface continuity rows replace them.  This metadata is
  //! intended for diagnostics so they do not check equations at rows that OCFESLV
  //! intentionally replaced by interface constraints.
  struct t_EqnNode
  {
    size_t eqn_index = 0;
    size_t flat_index = 0;
    t_Coord coord;
    std::map<FFVar,size_t,lt_FFVar> elem;
    std::map<FFVar,size_t,lt_FFVar> inode;
    bool base_keep = false;
    bool interface_donor = false;
    bool retained = false;
  };

  //! @brief Return physical nodes for equation `eqn_index`; optionally only rows retained by OCFESLV.
  std::vector<t_EqnNode> equation_nodes
    ( size_t eqn_index, bool retained_only=false )
    const;

  //! @brief Return physical nodes for equation expression `eqn`; optionally only rows retained by OCFESLV.
  std::vector<t_EqnNode> equation_nodes
    ( FFVar const& eqn, bool retained_only=false )
    const;


  //! @brief Reason for the most recent setup() call (OK if it succeeded).

  // setup_status() and setup_status_str() live in FFModel (shared with ODESLV).

protected:


protected:


  // ---------------------------------------------------------------------------
  // rev107 [taucol]: observation state for the direct named-column instrument.
  // ---------------------------------------------------------------------------
  // All MUTABLE for the same reason _stage2PreflightReported is: _eval_eqn() is
  // const and these record only what a diagnostic has seen.  NO BEHAVIOUR READS
  // THEM -- grep for each name before changing that.  Deliberately NOT propagated
  // by deep_copy_from(): a copied environment observes its own assembly.
  //
  // _tauHostRows is written ONCE, by the plan half, and only READ by the assembly
  // half.  That direction matters: the assembly half must not re-derive which rows
  // host a tau or the two halves stop being independent observations.








  // ------------------------------------------------------------------
  // Interface plan — owned by OCPlan, driven by OCFESLV::setup()
  // ------------------------------------------------------------------


  //! @brief Internal exact-interface claim occurrence key.
  struct t_InterfaceClaimOccurrenceKey {
    int    block_id;
    size_t eqn_id;
    size_t dom_id;
    size_t iel_lo;
    size_t face;
    int    kind;
    bool operator<( t_InterfaceClaimOccurrenceKey const& o ) const
    {
      if( block_id != o.block_id ) return block_id < o.block_id;
      if( eqn_id   != o.eqn_id   ) return eqn_id   < o.eqn_id;
      if( dom_id   != o.dom_id   ) return dom_id   < o.dom_id;
      if( iel_lo   != o.iel_lo   ) return iel_lo   < o.iel_lo;
      if( face     != o.face     ) return face     < o.face;
      return kind < o.kind;
    }
  };


public:


protected:




  //! @brief Forward-AD output function values used by deriv()
  mutable std::vector<t_FAD>         _FFct;


  //! @brief Cached sparse affine residual/output map for purely linear rows.
  struct t_LinearRow
  {
    double                           offset = 0.;
    std::vector<double>              val;
  };
  mutable std::vector<t_LinearRow>   _linEqn;
  mutable std::vector<t_LinearRow>   _linFct;
  mutable std::vector<double>        _linCstVal;
  mutable bool                       _linEqnReady = false;
  mutable bool                       _linFctReady = false;



  //! @brief Clear what a re-setup invalidates on the solver side (derivative caches, the interface plan, the
  //! per-block symbol/classification/face caches).  Overrides FFModel::_on_setup_begin().
  bool _on_setup_begin
    ()
    override;

  //! @brief Preflight the declared model before it is imported (overrides FFModel::_validate_model): the
  //! solver's check, which also validates against the marching grid.
  //! @brief Under SOLVE.MARCHING the working evolution domain is collapsed to its first element, but the march covers
  //! the full recorded grid: an output point on the evolution coordinate is checked against [t_0, t_N].
  void _output_point_bounds
    ( FFVar const& var, double& lo, double& up )
    const override
    {
      if( options.SOLVE.MARCHING && !_marchGrid.empty() && _evolution_dom_var.dag()
          && var.id().second == _evolution_dom_var.id().second ){
        lo = _marchGrid.front();
        up = _marchGrid.back();
      }
    }

  bool _validate_model
    ()
    override;

  //! @brief Before order reduction: hoist the evolution structure and, when marching, collapse the evolution
  //! domain to the first window -- reduction then runs against the horizon actually being solved.
  bool _on_before_reduction
    ()
    override;

  //! @brief OCFEBASE::_validate_input_grid_marching: under SOLVE.MARCHING the STATE evolution domain was collapsed
  //! to a single window by _collapse_evolution_for_marching(), but a distributed input keeps its FULL per-window
  //! grid (== the recorded march grid).  The residual eval locates the input element by absolute coordinate
  //! during each window solve, and the DOF count uses the input's OWN full grid, so register_control still sees
  //! all N control DOFs.  So the input is validated against the full march grid rather than the collapsed state
  //! domain.  Inert off the marching path (_marchGrid non-empty only after a genuine collapse; evolution
  //! direction of a distributed input only).
  bool _validate_input_grid_marching( FFVar const& v, FFDom const& d ) const override
  {
    if( !( options.SOLVE.MARCHING && !_marchGrid.empty() && _evolution_dom_var.dag()
        && v.id().second == _evolution_dom_var.id().second ) ) return false;
    if( d.n_elem + 1 != _marchGrid.size()
     || std::fabs( d.lo_dom - _marchGrid.front() ) > 64.*DBL_EPSILON
     || std::fabs( d.up_dom - _marchGrid.back()  ) > 64.*DBL_EPSILON )
      throw Exceptions( Exceptions::INDEX );
    return true;
  }

  //! @brief OCFEBASE::_init_input_nodes_marching: the twin of the above for the input's OWN node cache.  By this
  //! point _collapse_evolution_for_marching() has collapsed the STATE evolution domain in _mDom to a single window,
  //! but a distributed input must KEEP its full per-window partition (see _interp: a per-window solve leaves
  //! _lastSolveMarched==false, so _point_in_domain selects the element containing the window coordinate; window k
  //! reads the input's element k directly, no per-window value injection).  The partition is validated against the
  //! RECORDED march grid, kept, and only this input's node cache is (re)built.  Provably inert off the marching path.
  int _init_input_nodes_marching( FFVar const& dvar, FFDom& idom ) override
  {
    if( !( options.SOLVE.MARCHING && !_marchGrid.empty() && _evolution_dom_var.dag()
        && dvar.id().second == _evolution_dom_var.id().second ) ) return 0;
    if( idom.elem_bnd.size() != _marchGrid.size()
     || idom.n_elem + 1     != _marchGrid.size() )
      throw Exceptions( Exceptions::INDEX );      // input grid must match the full march grid
    for( size_t ie = 0; ie < _marchGrid.size(); ++ie )
      if( std::fabs( idom.elem_bnd[ie] - _marchGrid[ie] ) > 64.*DBL_EPSILON )
        throw Exceptions( Exceptions::INDEX );
    idom.lo_dom = _marchGrid.front();
    idom.up_dom = _marchGrid.back();
    if( !idom.set_nodes() ){ _setupStatus = SetupStatus::NODE_SETUP_FAILED; return -1; }
    return 1;                                     // keep the full input partition (no collapse)
  }

  bool _on_model_discretise() override
    { return _discretise_colloc_model(); }


  //! @brief Opt-in validation of user-supplied incoming BCs at hyperbolic boundary faces (count and
  //! characteristic direction), so a mis-placed inflow condition fails setup cleanly instead of producing a
  //! rank-deficient solve.  Always on (unconditional).
  bool _on_before_closures
    ()
    override;

  //! @brief After classification, the reduction plan and the closures: counts and subgraph caches, dependency
  //! map and evaluation plans, SAT coupling, the interface plan, the collocated caches and the audits.
  bool _on_model_ready
    ()
    override;

  // _on_model_imported -- the working input collocation (_mInpDisc) -- is OCFEBASE's: OCFEMOL needs it as much.

  //! @brief The declared model was replaced: drop the moving-mesh specifications, which name variables of the
  //! model that has just been discarded.
  void _on_model_replaced
    ()
    override
    { _adaptDomUsr.clear(); }


  //! @brief The classification is determined but not yet published: derive the per-(primitive,direction) symbol
  //! interface reads (cached in _blockInterfaceRead, consumed by the symbol-AUTO resolver).  Only needed when the
  //! interface type is actually derived (IC_AUTO); the balance probe forces it on as a diagnostic.
  void _on_classify_ready
    ( std::vector<double> const& state_ref, std::vector<double> const& input_ref,
      std::vector<double> const& cst_ref, std::vector<double> const& dom_ref )
    override
    {
      bool want_reads = ( options.INTERFACE.TYPE == Options::IC_AUTO );
#ifdef MC__OCFESLV_SYMBOL_BALANCE_PROBE
      want_reads = true;
#endif
      if( want_reads )
        _compute_interface_reads( state_ref, cst_ref, dom_ref, input_ref );
    }

  //! @brief The classification is published: it can change IC_AUTO dispatch and therefore the residual
  //! ordering/pattern for strong interface rows, so invalidate the derivative and linear-evaluation caches;
  //! rebuild now only in eager mode.
  //! @return false if an eager rebuild fails
  bool _on_classify_published
    ()
    override
    {
    // Classification can change IC_AUTO dispatch and therefore the residual
    // ordering/pattern for strong interface rows.  Invalidate derivative and
    // linear-evaluation caches; rebuild now only in eager mode.
    if( _issetup ){
      _jacReady = false;
      _linEqnReady = _linFctReady = false;
      _reset_colloc( _OCFVar, _OCFInp, _OCFDom, _OCFCst );
      _OCFEqn = OCVar<t_FAD>();
      _wrkOCF.clear();
      if( options.SOLVE.DERIV_CACHE == Options::DERIV_EAGER
       || options.SOLVE.LINEAR_EVAL_CACHE == Options::LINEAR_EAGER ){
        if( !_setup_deriv_cache() ) return false;
        _reset_autodiff( _FVar, _FInp, _FEqn, _FFct );
        if( options.SOLVE.LINEAR_EVAL_CACHE == Options::LINEAR_EAGER
         && !_setup_linear_cache( true, true, _vCstVal.empty()? nullptr: _vCstVal.data() ) )
          return false;
      }
    }
      return true;
    }

  //! @brief Invalidate the derivative/Jacobian caches after a model change (the FFModel mutators call this
  //! hook).
  void _on_model_changed( ModelChange const chg ) override
  {
    switch( chg ){
     case ModelChange::DERIVATIVES:  _jacReady = false; break;
    }
  }

  //! @brief The solver's Jacobian cache: the linear-residual cache reset, then OCFEBASE::_build_jacobian_pattern.
  bool _setup_deriv_cache
    ();

  //! @brief Ensure cached Jacobian sparsity/coloring match current setup/options.
  bool _ensure_deriv_cache
    ();



  //! @brief Build cached sparse affine maps for purely linear rows and fixed setup constants.
  bool _setup_linear_cache
    ( bool need_eqn, bool need_fct, double const* cst )
    const;

  //! @brief Try to evaluate rows through cached sparse affine maps.
  bool _eval_linear_cache
    ( double* eqn, double* fct, double const* var, double const* inp, double const* cst )
    const;

  //! @brief Check whether supplied constants match setup constants.
  bool _linear_cache_constants_match
    ( double const* cst )
    const;



  //! @brief Reset local working environment: the solver's state, then FFModel's (overrides FFModel::_reset).
  void _reset
    ()
    override;



  //! @brief Templated initialization of AD variable vectors
  template <typename U>
  void _reset_autodiff
    ( std::vector<FADType<U>>& Fvar, std::vector<FADType<U>>& Finp,
      std::vector<FADType<U>>& Feqn, std::vector<FADType<U>>& Ffct )
    const;






  //! @brief Templated sparse derivative evaluation of collocated output functions in compressed FADBAD forward mode
  template <typename U>
  void _deriv_fct
    ( U* fct, U** gradfct, size_t* nnzfct, size_t** colnzfct,
      OCVar<FADType<U>>& OCFfct, std::vector<OCVar<FADType<U>>>& OCFvar,
      std::vector<OCVar<FADType<U>>>& OCFinp, std::vector<OCVar<FADType<U>>>& OCFdom,
      std::vector<OCVar<FADType<U>>>& OCFcst, std::vector<OCVar<FADType<U>>>& OCFwrk,
      FADType<U>* Ffct, FADType<U>* Fvar, FADType<U>* Finp )
    const;


public:

  //! @brief Default Constructor
  OCFESLV
    ( FFGraph* dag = nullptr )
    : OCFEBASE( options )
    {
      set( dag );
    }

  //! @brief Copy constructor: rebuild the working DAG from the source OC model.  A copy of an environment that
  //! is not set up is itself an empty, not-set-up environment (is_setup() reports which).
  OCFESLV
    ( OCFESLV const& src )
    : OCFEBASE( options )
    {
      deep_copy_from( src );
    }

  //! @brief Assignment: rebuild the working DAG from the source OC model.  Assigning a not-set-up source leaves
  //! this environment empty and not set up (see the copy constructor).
  OCFESLV& operator=
    ( OCFESLV const& src )
    {
      if( this != &src ) deep_copy_from( src );
      return *this;
    }

  //! @brief Destructor: clears the solver state; FFModel's destructor then releases the declared model and
  //! the working DAG.
  virtual ~OCFESLV
    ()
    {
      _reset();
    }

  //! @brief The working domain map, as OCBase's collocation layer reads it: FFModel's (owner) map.  Always the
  //! working one -- the collocation layer runs during setup(), before _issetup is set.
  t_Dom const& _dom_colloc
    ()
    const
    override
    { return _mDom; }

  //! @brief Deep copy from another set-up environmentred and working model, options, analyses and the homotopy settings (parameters,
  //! remapped to this DAG, and per-stage step caps); the last continuation report is not copied.
  //! @return false if @p src is not set up, or if the DAG import fails; this environment is then left empty
  //! (not set up), never half-built.
  bool deep_copy_from
    ( OCFESLV const& src );


  //! @brief Human-readable name for an interface receiver-edge source.

//==========================================================================================================
// INTERFACE PLAN -- claims -> realisations -> rows and columns.  The DECISIONS live in ocplan.hpp; what is here
// is the assembly: which claims exist, which rows receive them, and how the chosen realisation becomes
// matrix structure.
// (banner placed from the definitions themselves: this region is lines 1459-1889 of the pre-rev268 file)
//==========================================================================================================


  //! @brief Attach the user DAG @p dag, clearing the solver state and the model (FFModel::set()).
  void set
    ( FFGraph* dag )
    override
    {
      _reset();
      _adaptDomUsr.clear();
      FFModel::set( dag );
    }


  //! @brief Flag a domain direction as MOVING-MESH (r-adaptive).
  //! Mirrors set_evolution_domain: records intent in _usr; validated at setup().
  //! The direction keeps its registered FFDom (finite-element count, collocation
  //! type, node style, and PHYSICAL extent lo_dom/up_dom all unchanged); at setup
  //! the header injects a solved mesh state x(xi,t) mapping this direction's
  //! computational xi onto [lo_dom, up_dom] via equidistribution with the given
  //! monitor.  Multiple directions may be flagged (tensor-product multi-D).
  //! No pre-check that add_domain(dom,...) preceded this -- deferred to setup().
  void add_adaptive_domain
    ( FFVar const& dom, MeshAdaptSpec const& spec )
    {
      _adaptDomUsr[ dom ] = spec;
      _issetup    = false;
      _classified = false;
      _jacReady   = false;
    }

  //! @brief Convenience: flag an adaptive direction with an explicit prescribed
  //! monitor expression M(xi,t) (the only monitor mode supported in M2a).
  void add_adaptive_domain
    ( FFVar const& dom, FFVar const& monitorExpr )
    {
      MeshAdaptSpec spec;
      spec.monitorExpr    = monitorExpr;
      spec.hasMonitorExpr = true;
      add_adaptive_domain( dom, spec );
    }

  //! @brief Clear the moving-mesh flag on one direction.
  void reset_adaptive_domain
    ( FFVar const& dom )
    {
      _adaptDomUsr.erase( dom );
      _issetup    = false;
      _classified = false;
      _jacReady   = false;
    }

  //! @brief Clear all moving-mesh flags.
  void reset_adaptive_domain
    ()
    {
      _adaptDomUsr.clear();
      _issetup    = false;
      _classified = false;
      _jacReady   = false;
    }

  //! @brief Read the registered moving-mesh specifications (per direction).
  std::map<FFVar,MeshAdaptSpec,lt_FFVar> const& adaptive_domains
    ()
    const
    { return _adaptDomUsr; }


  //! @brief Number of finite elements along the evolution/marching domain
  //! (0 when no evolution domain is set).  Basis for element-by-element
  //! marching (solve_marching(), forthcoming).
  size_t n_evolution_elem
    () const;

  //! @brief Number of designated/auto-detected marching transfer pairs (populated at
  //! setup() under SOLVE.MARCHING).  0 for a monolithic setup.
  size_t n_marching_transfers
    () const
    { return _marchTransfers.size(); }

  //! The control registry and its LAYOUT live in FFModel: add_input() fixes each control's function
  //! space (the domain's element partition and the input's own InpColloc{type,n_node}), so the DOF count
  //! is a property of the DECLARATION, not of this solver.  register_control(), clear_controls(),
  //! controls(), n_control_dof(), control_block() and ControlSpec are all inherited; OCFESLV used to keep
  //! a second copy and rediscover each count via get_input_values( input, dummy ).size(), which is the
  //! duplication this removes.  Only the cache invalidation stays here, as the _on_controls_changed()
  //! override below.

  //! @brief The reduced-space dependency pattern is built against the control vector, so any change to
  //! the registry or its layout invalidates it.
  void _on_controls_changed
    ()
    override
    { _reducedDepValid = false; }


  //! @brief Structural reduced-space dependency pattern: P[f][j] = true iff output function f can
  //! depend on control DOF j (j indexes the reduced control vector p, in canonical control order).
  //! Built from the residual/function structural Jacobians (_jacPatternEqn/_jacPatternFct, columns
  //! [states|inputs]).  Each control DOF maps to one input column; f depends on it (a) DIRECTLY if
  //! that column is in dF's pattern, or (b) THROUGH THE SOLVE if the column enters some residual
  //! whose states -- via J^{-1}, structurally DENSE within a state-coupling component -- reach a
  //! state f depends on.  For a fully-coupled monolithic block this is dense (correct for e.g. the
  //! MBC); for a structurally decomposable system it is block-sparse.  For a MARCHING solve the
  //! per-window pattern misses cross-window causality, so this returns a conservative dense pattern
  //! (the per-control first-window / per-function support-window causality refinement is TODO).
  std::vector<std::vector<bool>> const& reduced_dependency_pattern
    ()
  {
    if( _reducedDepValid ) return _reducedDepPattern;   // structural: computed once, reused by reference
    size_t const nf = n_colloc_fct(), ncd = n_control_dof();
    auto& P = _reducedDepPattern;
    P.assign( nf, std::vector<bool>( ncd, false ) );
    _reducedDepValid = true;
    if( !nf || !ncd ) return P;

    // Marching: cross-window causality not yet analysed -> conservative dense.
    if( is_marching() ){
      for( size_t f = 0; f < nf; ++f ) P[f].assign( ncd, true );
      return P;
    }

    if( !_ensure_deriv_cache() ){ for( size_t f=0; f<nf; ++f ) P[f].assign( ncd, true ); return P; }
    size_t const nv = n_colloc_sta(), ni = n_colloc_inp();

    // Map each control DOF (p-index) to its input column (nv + input position).
    std::vector<double> marker( ni );
    for( size_t i = 0; i < ni; ++i ) marker[i] = double(i);
    std::vector<size_t> p_incol( ncd, size_t(-1) );
    for( auto const& [input, spec] : _controls ){
      std::vector<double> const pos = get_input_values( input, marker.data() );
      for( size_t k = 0; k < spec.ndof && k < pos.size(); ++k )
        p_incol[ spec.offset + k ] = nv + size_t( pos[k] + 0.5 );
    }

    // State-coupling components: union states that co-occur in a residual (J^{-1} dense within).
    std::vector<size_t> comp( nv ); for( size_t s = 0; s < nv; ++s ) comp[s] = s;
    auto find = [&comp]( size_t x ){ while( comp[x] != x ){ comp[x] = comp[comp[x]]; x = comp[x]; } return x; };
    for( auto const& row : _jacPatternEqn ){
      size_t first = size_t(-1);
      for( size_t c : row ) if( c < nv ){ if( first == size_t(-1) ) first = c; else comp[find(c)] = find(first); }
    }

    // Input column -> set of state-components it perturbs (components of states in residuals it enters).
    std::map<size_t,std::set<size_t>> incol_comp;
    for( auto const& row : _jacPatternEqn ){
      std::set<size_t> rc; std::vector<size_t> ic;
      for( size_t c : row ){ if( c < nv ) rc.insert( find(c) ); else ic.push_back( c ); }
      for( size_t c : ic ) incol_comp[c].insert( rc.begin(), rc.end() );
    }

    // Captured reductions: a captured input cap = INT src carries the state-coupling of its SOURCE.
    // Map cap's input column(s) -> the state-components of src, so an output reading cap inherits src's
    // components (and thus every control that perturbs src through the solve).  Without this the pattern
    // sees f=cap depending only on the cap input and wrongly predicts dF/dp = 0.
    std::map<size_t,std::set<size_t>> capcol_comp;
    for( auto const& C : _deferredValue ){
      if( !C.input.dag() || !C.source.dag() ) continue;
      auto itb = _blkVar.find( C.source ); auto itv = _mVar.find( C.source );
      if( itb == _blkVar.end() || itv == _mVar.end() ) continue;
      size_t const off = std::get<0>( itb->second );
      size_t cnt = 1; for( auto const& v : itv->second ){ auto id = _mDom.find( v ); if( id != _mDom.end() ) cnt *= id->second.n_elem * id->second.n_node; }
      std::set<size_t> srccomp;
      for( size_t s = off; s < off + cnt && s < nv; ++s ) srccomp.insert( find(s) );
      std::vector<double> const cappos = get_input_values( C.input, marker.data() );
      for( double d : cappos ) capcol_comp[ nv + size_t( d + 0.5 ) ] = srccomp;
    }

    // Per function: its state-components and its direct input columns.
    for( size_t f = 0; f < nf && f < _jacPatternFct.size(); ++f ){
      std::set<size_t> fcomp, finp;
      for( size_t c : _jacPatternFct[f] ){
        if( c < nv ) fcomp.insert( find(c) );
        else { finp.insert( c ); auto ic = capcol_comp.find( c ); if( ic != capcol_comp.end() ) fcomp.insert( ic->second.begin(), ic->second.end() ); }
      }
      for( size_t j = 0; j < ncd; ++j ){
        size_t const ic = p_incol[j];
        if( ic == size_t(-1) ) continue;
        bool dep = finp.count( ic ) > 0;                    // direct dF/dinp
        if( !dep ){                                          // indirect through the solve
          auto it = incol_comp.find( ic );
          if( it != incol_comp.end() )
            for( size_t cc : it->second ) if( fcomp.count( cc ) ){ dep = true; break; }
        }
        P[f][j] = dep;
      }
    }
    return P;
  }

  //! @brief Decode a control vector @p p into the input array @p target (of size
  //! n_colloc_inp()): scatter each control's DOF block to its input via set_input_values.
  //! Works for VALUES (target=inp) and DIRECTIONS (target=dinp_dp, since the map is linear)
  //! -- for an IC-input control this writes the element-0 seed the march then propagates.
  void decode_controls
    ( std::vector<double> const& p, double* target ) const
    {
      for( auto const& [input, spec] : _controls ){
        if( spec.offset + spec.ndof > p.size() ) continue;
        std::vector<double> block( p.begin() + spec.offset, p.begin() + spec.offset + spec.ndof );
        const_cast<OCFESLV*>( this )->set_input_values( input, block, target );
      }
    }

  //! @brief Encode the current control-input values from @p inp into the control vector @p p.
  void encode_controls
    ( double const* inp, std::vector<double>& p ) const
    {
      p.assign( _nControlDof, 0. );
      for( auto const& [input, spec] : _controls ){
        auto block = const_cast<OCFESLV*>( this )->get_input_values( input, inp );
        for( size_t i = 0; i < spec.ndof && i < block.size(); ++i ) p[spec.offset + i] = block[i];
      }
    }

  //! @brief Output-function values accumulated across the elements of the most recent
  //! SOLVE.MARCHING solve (see _valFct).  Evolution-direction integrals are summed over
  //! elements to give the full-horizon value; terminal / final-time outputs take the last
  //! element.  Empty for a monolithic solve or a model without outputs.
  std::vector<double> const& val_functions
    () const
    { return _valFct; }

  //! @brief The recorded march grid (evolution-element boundaries, size N+1) when
  //! SOLVE.MARCHING collapsed a multi-element evolution declaration at setup(); empty
  //! for a monolithic setup.
  std::vector<double> const& march_grid
    () const
    { return _marchGrid; }

  //! @brief Number of march steps (evolution elements the user declared) under
  //! SOLVE.MARCHING; 0 for a monolithic setup.
  size_t n_march_steps
    () const
    { return _marchGrid.empty() ? 0 : _marchGrid.size() - 1; }

  //! @brief One evolution element's captured block solution from a SOLVE.MARCHING solve
  //! with OUTPUT.MARCH_STORE set: its window [t0,t1] and the collocated state/input
  //! values on that element.
  struct MarchBlock
  {
    double              t0 = 0., t1 = 0.;   //!< evolution-element window
    std::vector<double> var;                //!< collocated state solution on [t0,t1]
    std::vector<double> inp;                //!< input values on [t0,t1] (transferred ICs / feeds)
  };

  //! @brief Per-element block solutions captured during the most recent SOLVE.MARCHING
  //! solve when OUTPUT.MARCH_STORE was set (size n_march_steps()); empty otherwise.
  std::vector<MarchBlock> const& march_trajectory
    () const
    { return _marchTrajectory; }

  //! @brief Buffer-free read of the converged solution at point @p pt, for EITHER solve mode.
  //! Reads the framework's stored result: the cached primal (monolithic; needs SOLVE.REUSE) or
  //! the trajectory (marching; needs OUTPUT.MARCH_STORE) -- both are defaults.  const and
  //! side-effect-free.  Throws Exceptions::NOSTORE if no solution is cached.
  //! @param state  the state OR input to evaluate
  //! @param pt  the physical point, keyed by domain variable; under a marching solve the point itself selects the
  //! window, so the whole stored field is readable and not merely the last window
  //! @param side  2026-10-14: optional per-direction side (FFDom::MINUS, the default where absent, or FFDom::PLUS),
  //! for a coordinate where the value differs on either side of an element boundary or a transition.  Under WEAK
  //! imposition the two element traces are genuinely different numbers; under the exact modes they agree to the
  //! solver tolerance and the side only picks which element's polynomial is evaluated.  Forwarded to eval_colloc().
  double eval_solution
    ( FFVar const& state, t_Coord const& pt, t_Side const* side=nullptr ) const;

  protected:

  protected:
  //! @brief Forward-sensitivity SETUP: run the primal (nonlinear) march ONCE and, per
  //! evolution element, store the block solution, the sparse dR/dinp & dF Jacobians, and
  //! the dense LU factorization of the square block Jacobian J_k = dR_k/dvar_k.  These
  //! depend only on the primal solution, so they are computed once and reused across every
  //! seed by _march_fsens_apply().  Returns false if the primal march fails to converge.
  bool _march_fsens_setup
    ( double* var, double const* inp, double const* cst );

  //! @brief Forward-sensitivity APPLY: given an input-space seed @p dinp_dp (the element-0
  //! control direction d(inp)/d(scalar)), propagate the LINEAR state sensitivity through the
  //! stored blocks -- per element, RHS = -(dR/dinp . dinp_dp) with the IC-input entries
  //! overwritten by terminal(s_{k-1}); s solved by the stored LU (triangular solves only);
  //! dFdp += dF/dinp . dinp_dp + dF/dvar . s (evolution-integral summed / else last).  Cheap
  //! and reusable: call once per FADType seed after a single _march_fsens_setup().
  bool _march_fsens_apply
    ( std::vector<double> const& dinp_dp, std::vector<double>& dFdp );

  //! @brief Adjoint-sensitivity APPLY: transpose (reverse) of _march_fsens_apply.  Given a
  //! function-space weight @p w (size n_colloc_fct()), propagate BACKWARD through the stored
  //! blocks to accumulate the gradient of w.F in INPUT space @p g (size n_colloc_inp()):
  //! per element (reverse order) lb = J_k^{-T} (dF/dvar^T wk) using the stored LU transposed,
  //! g += dF/dinp^T wk - dR/dinp^T lb, with the IC-input entries gathered back into the
  //! previous block's terminal adjoint via terminal_profile_adjoint() (the transpose of the
  //! forward terminal transfer).  wk applies the evolution-integral-summed / else-last rule
  //! of the forward.  Cheap and reusable: call once per adjoint seed after a single
  //! _march_fsens_setup().  The reduced-space gradient in control space is then encode_controls(g).
  bool _march_asens_apply
    ( std::vector<double> const& w, std::vector<double>& g );

  //! @brief Adjoint (transpose) of terminal_profile(): scatter per-state terminal-node
  //! adjoints @p lam (ordered exactly as terminal_profile()'s output) into the state block
  //! adjoint @p gvar (size n_colloc_var(), ACCUMULATED into) using the evolution-direction
  //! Lagrange weights at the UB.  Since terminal_value == eval_colloc pinned to the UB with
  //! spatial nodes at collocation points, M is L_m(UB) along evolution x selection in space;
  //! this applies M^T.  Window-invariant (reference-element weights), so no reset needed.
  void terminal_profile_adjoint
    ( std::map<FFVar,std::vector<double>,lt_FFVar> const& lam, double* gvar ) const;


  public:
  //! @brief Whether the most-recent setup recorded a marchable evolution grid AND marching
  //! is enabled: dispatches the reduced-space sensitivity API between the marching and the
  //! monolithic single-block backends.
  bool is_marching
    () const
    { return options.SOLVE.MARCHING && n_march_steps() > 0; }

  protected:
  //! @brief Monolithic forward-sensitivity SETUP: solve the FULL (non-marched) system, then
  //! capture the single-block linearization (LU of J=dR/dvar, sparse dR/dinp & dF, and the
  //! primal function values in _monoFct).  The N=1, no-coupling counterpart of
  //! _march_fsens_setup().  Returns false if the monolithic solve fails.
  bool _mono_fsens_setup
    ( double* var, double const* inp, double const* cst );

  //! @brief Monolithic forward-sensitivity APPLY: dFdp = dF/dinp.dinp + dF/dvar . (J^{-1}(-dR/dinp.dinp)).
  //! One triangular solve per seed after a single _mono_fsens_setup().  No fctsum accumulation
  //! (single block: every output row is direct).
  bool _mono_fsens_apply
    ( std::vector<double> const& dinp, std::vector<double>& dFdp );

  //! @brief Monolithic adjoint-sensitivity APPLY: transpose of _mono_fsens_apply.  Given a
  //! function-space weight @p w (size n_colloc_fct()), g = dF/dinp^T w - dR/dinp^T (J^{-T}(dF/dvar^T w))
  //! in INPUT space (size n_colloc_inp()).  No terminal coupling (single block).
  bool _mono_asens_apply
    ( std::vector<double> const& w, std::vector<double>& g );


  //! @brief Reduced-space sensitivity SETUP dispatch (marching vs monolithic per is_marching()):
  //! the primal solve + factorization shared by forward and adjoint.  Internal building block
  //! of solve_fsens()/solve_asens(); also leaves sens_functions() populated.
  bool _sens_setup
    ( double* var, double const* inp, double const* cst )
    { return is_marching() ? _march_fsens_setup( var, inp, cst ) : _mono_fsens_setup( var, inp, cst ); }

  //! @brief Forward-sensitivity directional APPLY (auto): dFdp for one input-space direction.
  bool _solve_fsens_apply
    ( std::vector<double> const& dinp, std::vector<double>& dFdp )
    { return is_marching() ? _march_fsens_apply( dinp, dFdp ) : _mono_fsens_apply( dinp, dFdp ); }

  //! @brief Adjoint-sensitivity weighted APPLY (auto): input-space gradient for one function weight.
  bool _solve_asens_apply
    ( std::vector<double> const& w, std::vector<double>& g )
    { return is_marching() ? _march_asens_apply( w, g ) : _mono_asens_apply( w, g ); }

  public:
  //! @brief COMPLETE reduced-space forward-sensitivity solve.  Like solve() evaluates the primal
  //! and its functions, solve_fsens() evaluates the functions AND their derivatives w.r.t. the
  //! registered controls: it runs _sens_setup() then one forward apply per control DOF, populating
  //! sens_functions() and the reduced Jacobian sens_jacobian() (n_colloc_fct x n_control_dof,
  //! row-major).  Auto-dispatches marching vs monolithic.  Requires registered controls.
  bool solve_fsens
    ( double* var, double const* inp, double const* cst );

  //! @brief COMPLETE reduced-space adjoint-sensitivity solve.  As solve_fsens() but assembling the
  //! reduced Jacobian by ONE reverse sweep per output function (cheaper when n_colloc_fct <
  //! n_control_dof).  Identical functions and Jacobian to solve_fsens(); populates the same caches.
  bool solve_asens
    ( double* var, double const* inp, double const* cst );

  //! @brief Output-function values captured by the last solve_fsens()/solve_asens()/_sens_setup().
  std::vector<double> const& sens_functions
    () const
    { return _sensFct; }

  //! @brief Reduced Jacobian dF/dp (row-major, index i*n_control_dof + j = dF_i/dp_j) captured by
  //! the last solve_fsens()/solve_asens().  Empty until one has been called.
  std::vector<double> const& sens_jacobian
    () const
    { return _sensJac; }


  //! @brief Self-check that the evolution slabs partition the evolution-bearing
  //! state DOFs exactly (pairwise disjoint, union covers every such DOF).
  //! Read-only diagnostic to validate the slab column map before the windowed
  //! marching solve is built; DISPLAY_LEVEL>=1 prints a per-slab summary.
  bool evolution_slab_selfcheck
    () const;

  //! @brief Column -> evolution-element index for every collocated state DOF (via the
  //! slab column map); DOFs not in any evolution slab (tau/interface/lumped) are -1.
  std::vector<int> evolution_col_to_elem
    () const;

  //! @brief Structure audit of the evolution-slab decomposition of the assembled square
  //! system.  Assigns each residual row to the HIGHEST evolution element among the state
  //! columns it touches (Jacobian sparsity), then reports per-element row/column counts,
  //! rows/columns not attributable to one element (interface/tau coupling), whether the
  //! per-element blocks are square, and whether the coupling is block-lower-BIDIAGONAL
  //! (rows of element k touch only elements {k-1,k}) -- the prerequisites for a
  //! forward-substitution marching solve.  DISPLAY_LEVEL>=1 prints the per-element map.
  bool evolution_slab_map_selfcheck
    () const;

  //! @brief Value of an evolution-bearing state at the evolution UPPER bound (the
  //! marching terminal), interpolated from the solved coefficients @p var at the
  //! given spatial point.  Family-agnostic: it interpolates to the UB rather than
  //! assuming a UB collocation node, so it is correct for LG/LGL/CGL/LGR.  The
  //! atomic transfer read for element-by-element marching (Stage-0).
  double terminal_value
    ( FFVar const& state, t_Coord spatial_pt,
      double const* var, double const* inp=nullptr, double const* cst=nullptr ) const;

  //! @brief Terminal profile -- state values at the evolution UB for every
  //! evolution-bearing state, over all spatial collocation nodes.  Ordered by the
  //! spatial-node tensor product in var_dom order (first spatial domain innermost),
  //! element-major/node-inner per domain, so the re-seed step (forthcoming) can
  //! write it back onto the next element's LB-face IC.  Built on terminal_value().
  std::map<FFVar,std::vector<double>,lt_FFVar> terminal_profile
    ( double const* var, double const* inp=nullptr, double const* cst=nullptr ) const;


  //! @brief The forward-differentiated model along the directions @p dirs, each a (declared input, DOF index) pair
  //! or a (constant, 0) pair: a NEW OCFESLV (not set up) on the same user DAG, whose nDir*nf outputs are
  //! f[i*nf+j] = dF_j/d(dir_i) -- the legacy layout of ODESLVS::fdiff(nPar,pPar), which this mirrors.  The product
  //! takes exactly the original's inputs and constants: every direction input is FIXED to its unit seed (which
  //! OCFESLV honours), and the original's own fixed inputs stay fixed.  Options are copied.  nullptr and a reason
  //! in @p err on an unknown entry or a DOF out of range.
  OCFESLV* fdiff
    ( std::vector<std::pair<FFVar,size_t>> const& dirs, std::string& err )
    const
    {
      err.clear();
      size_t const n = dirs.size();
      if( !n ){ err = "fdiff: no direction given"; return nullptr; }
      std::vector<FFVar> vU, vC;  std::vector<int> iu( n, -1 ), ic( n, -1 );
      for( size_t i=0; i<n; ++i ){
        FFVar const& v = dirs[i].first;
        bool isin = false; for( auto const& [w,dom] : var_declared_input() ) if( w.id().second == v.id().second ){ isin = true; break; }
        if( isin ){
          if( dirs[i].second >= input_ndof( v ) ){ err = "fdiff: DOF " + std::to_string( dirs[i].second ) + " out of range for " + v.name(); return nullptr; }
          for( size_t c=0; c<vU.size(); ++c ) if( vU[c].id().second == v.id().second ){ iu[i] = (int)c; break; }
          if( iu[i] < 0 ){ iu[i] = (int)vU.size(); vU.push_back( v ); }
          continue;
        }
        bool iscst = false; for( auto const& k : _usr._vCstUsr ) if( k.id().second == v.id().second ){ iscst = true; break; }
        if( !iscst ){ err = "fdiff: " + v.name() + " is neither a declared input nor a constant"; return nullptr; }
        for( size_t c=0; c<vC.size(); ++c ) if( vC[c].id().second == v.id().second ){ ic[i] = (int)c; break; }
        if( ic[i] < 0 ){ ic[i] = (int)vC.size(); vC.push_back( v ); }
      }
      std::vector<std::vector<double>> seedC( n, std::vector<double>( vC.size(), 0. ) );
      for( size_t i=0; i<n; ++i ) if( ic[i] >= 0 ) seedC[i][ ic[i] ] = 1.;
      std::vector<std::vector<FFVar>> vDU, vS;
      OCFESLV* sens = new OCFESLV( _usr._dagUsr );
      if( !FFModel::fdiff( *sens, vU, n, vDU, vS, err, /*keep_originals=*/false, vC, seedC ) ){ delete sens; return nullptr; }
      for( size_t i=0; i<n; ++i )                            // direction i: unit seed at its DOF, zero elsewhere
        for( size_t c=0; c<vU.size(); ++c ){
          std::vector<double> e( sens->input_ndof( vDU[i][c] ), 0. );
          if( iu[i] == (int)c && dirs[i].second < e.size() ) e[ dirs[i].second ] = 1.;
          if( !sens->fix_input( vDU[i][c], e ) ){ err = "fdiff: cannot fix the direction input " + vDU[i][c].name(); delete sens; return nullptr; }
        }
      for( auto const& [w,vals] : _usr._mInpFixUsr ) sens->fix_input( w, vals );
      sens->options = options;
      return sens;
    }

  //! @brief Call depth of the public entry points that apply fixed-input values.  The values are applied to the
  //! CALLER's inp at the OUTERMOST call only.  Calls nested inside a solve work on the solver's own inp, which
  //! already carries them -- and which the solver may legitimately have written since: marching's canonical
  //! transfer writes each window's terminal into the IC input's slot for the next window, and for an fdiff product
  //! that IC input is a FIXED direction input.  Re-applying the fixed values there reset every sensitivity state to
  //! its seed at each window (OCFE_SYMDIFF1, 2026-09-27: derivatives of exactly one window).
  mutable int _fixedInputDepth = 0;
  struct FixedInputScope {
    int& d;
    explicit FixedInputScope( int& d_ ) : d( d_ ) { ++d; }
    ~FixedInputScope() { --d; }
  };

  //! @brief The input vector with the values of every input held by fix_input written over @p inp: returns @p inp
  //! itself when no input is fixed, else a copy in @p buf.  Every public entry point taking a caller's inp goes
  //! through this, so a fixed input's values always come from the model.
  double const* _with_fixed_inputs
    ( double const* inp, std::vector<double>& buf ) const
    {
      if( _fixedInputDepth > 1 ) return inp;           // nested inside a solve: the solver's own inp, see below
      if( _usr._mInpFixUsr.empty() || !_nCollInp ) return inp;
      buf.assign( _nCollInp, 0. );
      if( inp ) std::copy( inp, inp + _nCollInp, buf.begin() );
      for( auto const& [w,vals] : _usr._mInpFixUsr ) set_input_values( w, vals, buf.data() );
      return buf.data();
    }


  //! @brief Exact marching transfer: seed distributed input @p input (in inp[]) with
  //! @p state evaluated at the evolution UB, sampled at the INPUT's own collocation
  //! nodes.  Because an input is read by eval()/solve() on its own grid, sampling the
  //! state's terminal there makes the re-seeded IC exact even if the input and state
  //! grids differ.  Combines terminal_value() + set_input_values().  This is the
  //! per-element transfer that advances a Stage-0 march.
  bool transfer_terminal
    ( FFVar const& state, FFVar const& input,
      double const* var, double* inp, double const* cst=nullptr ) const;

  //! @brief Broadcast a spatial IC profile across every evolution node of @p state in
  //! var[] -- a constant-in-evolution initial guess for the next marching element's
  //! solve (the BROADCAST_IC warm-start).  @p ic_profile is indexed like
  //! terminal_profile() (spatial nodes, element-major / node-inner).  Returns false on
  //! unknown state, no evolution domain, or a size mismatch.
  bool warmstart_broadcast
    ( FFVar const& state, std::vector<double> const& ic_profile, double* var ) const;

  //! @brief Slide a ONE-element evolution template to the window [t0, t0+w] in place,
  //! with NO re-setup().  Moving both origin and width lets a symbolic forcing/BC (a DAG
  //! function of the evolution coordinate) auto-evaluate at the correct ABSOLUTE window,
  //! which is what lets solve()-marching advance a non-autonomous problem without input
  //! re-seeding.  A no-op when the window is unchanged.  The differentiation reads the
  //! geometry live via the OCVar eval path, so the new window takes effect on the next
  //! solve().  Requires options.SOLVE.LINEAR_EVAL_CACHE == LINEAR_NONE (the default).
  bool reset_evolution_window
    ( double t0, double w );

  //! @brief Rewrite the width of a ONE-element evolution template to @p w in place,
  //! keeping the current origin (convenience wrapper over reset_evolution_window).
  bool reset_evolution_width
    ( double w );


  //! @brief Evaluate collocated equation residuals and post-evaluation outputs in double arithmetic
  bool eval
    ( double* eqn, double* fct, double const* var, double const* inp, double const* cst );

  //! @brief Evaluate collocated equation residuals and post-evaluation outputs in FFDep arithmetic
  template <typename CST>
  bool eval
    ( FFDep* eqn, FFDep* fct, FFDep const* var, FFDep const* inp, CST const* cst );

  //! @brief Evaluate collocated residuals and outputs and their sparse derivatives wrt state and input values. Nonzero derivative columns are stacked as [var, inp], so an input column i is returned as _nCollVar+i in colnzeqn/colnzfct.
  bool deriv
    ( double* eqn, double** gradeqn, size_t* nnzeqn, size_t** colnzeqn,
      double* fct, double** gradfct, size_t* nnzfct, size_t** colnzfct,
      double const* var, double const* inp, double const* cst );

  //! @brief Return the sparse derivative pattern for all rows. Rows are stacked as [equations, outputs], and columns as [collocated states, collocated inputs].
  bool deriv
    ( size_t* nnz, size_t** colnz );

  //! @brief Evaluate the sparse derivative values for all rows, stacked row-by-row according to deriv(nnz,colnz).
  bool deriv
    ( double* gradeqn, double* gradfct, double const* var, double const* inp, double const* cst );

#ifdef MC__OCFESLV_JAC_FD_PROBE
  //! @brief [PROBE, compile with -DMC__OCFESLV_JAC_FD_PROBE] Compare the analytic (FADBAD) equation
  //! Jacobian against central finite differences for every column of one STATE, at the given point.
  //! Localizes cross-domain derivative errors (e.g. a {z}-state whose derivative is not carried
  //! through the domain-union broadcast into 2-D equations).  Prints per-column max|fd-analytic|.
  void probe_jac_fd
    ( FFVar const& state, double* var, double const* inp, double const* cst, double h );
#endif

  //! @brief Outcome and diagnostics of the internal LM/Newton solve (OCFESLV::solve()).
  struct SolveReport
  {
    bool   converged       = false;  //!< max|r| <= options.SOLVE.RES_TOL at return
    bool   reused          = false;  //!< returned from the SOLVE.REUSE primal cache (no Newton run)
    int    iterations      = 0;      //!< accepted outer iterations
    int    newton_steps    = 0;      //!< accepted undamped-Newton steps
    int    lm_steps        = 0;      //!< accepted damped-LM steps
    double initial_residual= 0.0;    //!< max|r| at the supplied guess
    double final_residual  = 0.0;    //!< max|r| at return
    double final_merit     = 0.0;    //!< 0.5 * ||row-scaled r||^2 at return
  };

  //! @brief Solve the square collocated residual system r(var)=0 on the setup() layout
  //!        with an equilibrated Levenberg-Marquardt / Newton hybrid (sparse AD Jacobian).
  //!
  //! Operates in place on @p var (size n_colloc_sta(), tau columns already zero-initialised
  //! by setup(); the caller supplies the initial guess).  Mode-agnostic: IC_TRACE trace taus
  //! are ordinary var[] columns and IC_STRONG Schur taus are eliminated inside the residual /
  //! Jacobian evaluation, so the solver never branches on the interface mode.  Behaviour is
  //! governed by the Options SOLVE_* fields.  Does NOT call verify_interface_drop(): the
  //! post-solve continuity assertion remains a separate, explicit step.
  SolveReport solve
    ( double* var, double const* inp = nullptr, double const* cst = nullptr );

  //! @brief Status of a by-name solve, as BASE_DE's: NORMAL on convergence, FAILURE when the solve does not
  //! converge, FATAL when the model is not set up or init() fails.
  enum STATUS{ NORMAL=0, FAILURE, FATAL };
  //! @brief Solve with the input and constant values given BY NAME, as ODESLV_CVODES and DAESLV_IDAS do:
  //! { {u, {v1,...,vn}}, ..., {K, {k}} } -- every declared input (except a fixed one) and every constant exactly
  //! once, an input's values a flat list in input_dofs() order (OCFESLV's own element-major / node-inner order)
  //! or a generator, a constant one value.  The initial guess comes from the references (init()); the solution is
  //! left in val_state(), the outputs in val_function() -- in declaration order, as eval() returns them -- and the
  //! solver's report in last_report().  Throws std::invalid_argument, naming the entry, on an incomplete or
  //! inconsistent list.
  STATUS solve
    ( std::vector<FFModel::InputVal> const& vIn );
  //! @brief The collocated solution of the last by-name solve.
  std::vector<double> const& val_state() const { return _valSta; }
  //! @brief The outputs of the last solve, in declaration order -- the name ODESLV_CVODES and DAESLV_IDAS use for
  //! val_functions(): filled by the solve itself, with captures (monolithic) or accumulated across elements
  //! (SOLVE.MARCHING).
  std::vector<double> const& val_function() const { return _valFct; }
  //! @brief The solver's report of the last by-name solve.
  SolveReport const& last_report() const { return _lastReport; }


  // the by-name solve's results (its outputs are the solver's own _valFct)
  std::vector<double> _valSta;
  SolveReport         _lastReport;

  // ============================ HOMOTOPY / CONTINUATION ============================
  //
  //  A homotopy parameter is an input driven from `from` to `to` by a scalar s in [0,1].  When
  //  options.HOMOTOPY.ENABLE is set and parameters are registered, solve() runs the continuation
  //  (_solve_homotopy()), stepping s adaptively; each step re-enters solve().  Stages define the
  //  schedule: parameters sharing a stage move together and stages run in ascending order
  //  (all stage 0: simultaneous; 0,1,2,...: staircase; mixed: partial grouping).
  //
  //  Continuation also SELECTS THE ROOT on multi-rooted models: a different schedule or step cap
  //  can converge to a different solution with no warning, so validate such a change against a
  //  known-good final state.  (Measured on MBC: two paths converged to |r| ~ 1e-11 yet differed
  //  by 6.4 in the infinity norm.)
  //
  //  Solver configuration, not model declaration: parameters are looked up among the inputs of
  //  the set-up model.  History: in FFModel rev187-rev188, returned to OCFESLV in rev189
  //  (NOTES_20260911h_rev189).
  bool add_homotopy
    ( FFVar const& param, double from, double to, int stage = 0 );

  //! @brief Register input @p param as a homotopy parameter with a map s -> value, e.g. for a rate
  //! constant spanning decades: [](double s){ return 2e4*std::%pow(50.,s); }.  Must be called
  //! after setup().
  //! To step through a LIST OF VALUES, interpolate them in the map -- values[i] at s = i/(n-1), joined
  //! linearly, is the waypoint idiom, and `{from, to}` reproduces the ramp overload above:
  //! \code
  //!   std::vector<double> const v{ 0.1, 0.4, 1.0 };          // CAPTURE BY VALUE: the map outlives the call
  //!   oc.add_homotopy( p, [v]( double s ){
  //!     double const x = std::min( std::max( s, 0. ), 1. ) * ( v.size() - 1 );
  //!     size_t const i = std::min( (size_t)x, v.size() - 2 );
  //!     return v[i] + ( x - (double)i ) * ( v[i+1] - v[i] );
  //!   } );
  //! \endcode
  //! A bent path like that is useful where a straight ramp leaves the branch.  Note the stepper chooses its
  //! own s, so the PATH passes through the values but a solve need not LAND on one: to force a stop at a
  //! value, give it its own stage or cap that stage with set_homotopy_cap().
  //! @return false if @p param is not an input of the set-up model or @p map is empty
  bool add_homotopy
    ( FFVar const& param, std::function<double(double)> map, int stage = 0 );

  //! @brief Remove all homotopy parameters and the per-stage step caps set by set_homotopy_cap().
  void reset_homotopy
    ();

  //! @brief Override options.HOMOTOPY.STEP_CAP for a single stage.
  //!
  //! The cap is a PATH-FIDELITY setting, and different stages of the same model can tolerate very
  //! different steps: a stage that merely deforms a profile may run happily at 0.1 while one that
  //! switches on a stiff coupling needs 0.05 or finer to stay on the same branch.  With one global
  //! cap the whole schedule pays for the worst stage.  Pass @p cap <= 0 to clear the override.
  void set_homotopy_cap
    ( int stage, double cap );

  //! @brief Per-stage cost breakdown.
  struct StageReport
  {
    int    stage       = 0;      //!< stage index
    int    nparam      = 0;      //!< parameters driven together in this stage
    bool   converged   = false;  //!< stage reached s = 1
    int    solves      = 0;      //!< solve() calls
    int    accepts     = 0;      //!< accepted steps
    int    backtracks  = 0;      //!< rejected steps (step halved)
    int    iterations  = 0;      //!< summed Newton/LM iterations
    int    pred_used   = 0;      //!< predictor accepted (some extrapolation beat plain reuse)
    int    pred_rejected = 0;    //!< predictor rejected (plain reuse was best)
    int    pred_ord1   = 0;      //!< accepted steps where the SECANT (2-point) won
    int    pred_ord2   = 0;      //!< accepted steps where the QUADRATIC (3-point) won.  Watch this:
                                 //!< pred_used alone cannot distinguish a better prediction from the
                                 //!< same one, so a rising ord2 (with falling `iterations`) is what
                                 //!< shows the higher order actually paying.
    double seconds     = 0.0;    //!< wall time in solve()
    double s_reached   = 0.0;    //!< s attained (1.0 on success)
    double last_residual = 1e30; //!< max|r| at the last solve of this stage
    int    destaged    = 0;      //!< 0 normally; else the number of sequential sub-stages the stage
                                 //!< was split into after failing to move its parameters together
  };

  //! @brief Detail of a continuation run; retrieve with continuation_report() after solve().
  struct ContinuationReport
  {
    bool   converged      = false;  //!< every stage reached s = 1
    int    solves         = 0;      //!< total solve() calls
    int    accepts        = 0;
    int    backtracks     = 0;
    int    iterations     = 0;      //!< total Newton/LM iterations
    int    pred_used      = 0;
    int    pred_rejected  = 0;
    int    pred_ord1      = 0;      //!< accepted steps won by the secant
    int    pred_ord2      = 0;      //!< accepted steps won by the quadratic
    double seconds        = 0.0;    //!< total wall time in solve()
    double final_residual = 1e30;   //!< max|r| at return
    int    failed_stage   = -1;     //!< index of the stage that failed (-1 if none)
    double failed_at_s    = 0.0;    //!< s reached when it failed
    std::vector<StageReport> stage; //!< per-stage breakdown
  };

  //! @brief Per-stage detail of the continuation run inside the last solve().
  //!
  //! solve() returns the usual aggregate SolveReport; this returns the stage-by-stage breakdown
  //! (solves, backtracks, predictor accept/reject, timings, and which stage failed).  Empty when
  //! the last solve() ran no continuation.
  ContinuationReport const& continuation_report
    ()
    const
    { return _homReport; }

private:

  //! @brief Homotopy parameter record.
  struct t_Homotopy
  {
    FFVar  param;                        //!< the driven input
    double from = 0.0, to = 1.0;         //!< linear endpoints (ignored when map is set)
    std::function<double(double)> map;   //!< optional s -> value override
    int    stage = 0;                    //!< schedule stage
  };

  //! @brief Homotopy parameters in registration order (carried by deep copies).
  std::vector<t_Homotopy> _vHom;

  //! @brief Value of homotopy entry @p h at continuation coordinate @p s.
  double _homotopy_value
    ( t_Homotopy const& h, double s )
    const
    { return h.map ? h.map( s ) : h.from + s * ( h.to - h.from ); }

  //! @brief Continuation driver: ramp the registered homotopy schedule, calling the ordinary
  //! solve() at each step (guarded by _inHomotopy so the inner call does not re-enter).  Mirrors
  //! the SOLVE.MARCHING pattern: dispatched from solve(), returns an aggregate SolveReport
  //! (converged only if every stage reached s=1), and stashes the detail in _homReport.
  //!
  //! The caller's @p inp is NOT modified: the parameter values are driven on an internal copy, so
  //! a continuation solve has no side effect on the caller's input array.
  SolveReport _solve_homotopy
    ( double* var, double const* inp, double const* cst );

  ContinuationReport _homReport;        //!< detail of the last continuation run
  bool               _inHomotopy = false; //!< guard: inner solve() must not re-enter the ramp
  std::map<int,double> _mHomCap;         //!< per-stage HOMOTOPY.STEP_CAP overrides (set_homotopy_cap; carried by deep copies)

  //! @brief Step cap in force for @p stage: the override if one was set, else options.HOMOTOPY.STEP_CAP.
  double _homotopy_cap
    ( int stage )
    const
    { auto it = _mHomCap.find( stage ); return it != _mHomCap.end() ? it->second : options.HOMOTOPY.STEP_CAP; }

  //! @brief Ramp @p act (the parameters active in this leg) from @p s0 to @p s1 with the adaptive
  //! step and the guarded predictor, accumulating cost into @p st.  Returns true on reaching @p s1;
  //! on failure @p s_out holds the coordinate attained and the state/inputs are restored to it.
  //! Used both for a whole stage and, after an auto-destage, for each single parameter in turn.
  bool _homotopy_ramp
    ( std::vector<t_Homotopy const*> const& act, double s0, double s1, double cap,
      double* var, double* pinp, double const* cst,
      std::vector<double>& saved, StageReport& st, double& s_out );

  //! @brief max|r| at a candidate state (one residual evaluation, no Jacobian).  Used to guard
  //! the predictor so it can never hand the corrector a guess worse than plain reuse.
  double _homotopy_resid
    ( double const* var, double const* inp, double const* cst )
    const;
public:

private:
  //! @brief SOLVE.MARCHING driver: march the recorded grid, sliding the single-element
  //! template per step (reset_evolution_window), transferring the IC (transfer_terminal
  //! on the auto-detected transfer state), applying the SOLVE.WARMSTART strategy, and
  //! calling the monolithic solve() on each block (guarded by _inMarch).  Returns an
  //! aggregate report (converged only if every step converged).
  SolveReport _solve_marching
    ( double* var, double const* inp, double const* cst );
public:

  //! @brief Outcome of the post-solve interface-drop verification gate.
  //!
  //!   OK       : no dropped continuity claim is violated on the solution.
  //!   OVERDROP : a dropped claim's continuity residual exceeds rel_tol, i.e.
  //!              the redundancy DROP removed a NEEDED constraint and the solve
  //!              converged (cleanly, no residual signal) to a continuity-
  //!              violating solution.  This is PDE7's silent-wrong failure mode.
  //!   UNDERDROP: reserved for the conditioning probe (a retained-but-redundant
  //!              claim leaving the saddle near-singular); NOT computed here.
  enum class InterfaceDropStatus { OK, OVERDROP, UNDERDROP };

  //! @brief Post-solve over-drop gate for the redundancy DROP (refactor item 10,
  //!        eval-stage verification).
  //!
  //! Each claim in the committed drop set was removed at setup on the assumption
  //! that its continuity is *implied* by the retained constraints plus the PDE
  //! rows it couples through.  That assumption is exact for symmetric (value-C0 /
  //! parabolic) coupling but can be wrong for directional (upwind) coupling,
  //! where the solver then converges with a clean residual to a continuity-
  //! violating solution (no signal).  This gate re-evaluates each dropped claim's
  //! continuity residual on the converged @p var and flags OVERDROP if the worst
  //! relative violation exceeds @p rel_tol, turning that silent failure loud.
  //!
  //! @param var      converged collocated state vector (same layout as eval()'s var)
  //! @param rel_tol  relative continuity-violation threshold.  Must sit ABOVE the
  //!                 discretization truncation floor (a soundly-dropped claim's
  //!                 implied continuity holds only up to truncation, e.g. a flux
  //!                 dup-spread ~1e-3 at the current resolution) and BELOW genuine
  //!                 physical jumps (a wrong drop violates by O(1)); default 1e-2.
  //! Always prints a one-line summary (worst relative violation, claims checked)
  //! to std::cerr so the threshold can be calibrated from run logs.
  //!
  //! @param overdropped  if non-null, receives the frozen keys of the claims
  //!                 whose continuity is violated beyond @p rel_tol -- exactly
  //!                 the set the DROP_REDERIVE driver loop protects and re-derives
  //!                 (item 10 step 3).
  InterfaceDropStatus verify_interface_drop
    ( double const* var, double rel_tol = 1.0e-2,
      std::set<t_FrozenInterfaceClaimKey>* overdropped = nullptr ) const;

  //! @brief Pin continuity claims so the redundancy detector keeps them explicit
  //!        on the next setup() (item 10 step 3, DROP_REDERIVE).  The driver
  //!        feeds it the over-drop set from verify_interface_drop, then re-runs
  //!        setup()/solve.  Additive across calls (protected claims accumulate
  //!        over re-derive iterations); persists across setup() until cleared.
  void protect_interface_claims( std::set<t_FrozenInterfaceClaimKey> const& claims )
    { _strong_protected_claims.insert( claims.begin(), claims.end() ); }


  //! @brief The continuity claims currently pinned against dropping.
  std::set<t_FrozenInterfaceClaimKey> const& protected_interface_claims() const
    { return _strong_protected_claims; }

  //! @brief Driver-facing re-derive step (item 10 step 3).  Runs the post-solve
  //!        over-drop gate (verify_interface_drop) and, if it fires while the
  //!        policy is DROP_REDERIVE, pins the offending claims so the next
  //!        setup() keeps them explicit, then returns OVERDROP to signal
  //!        "re-derive and re-solve".  Returns OK when the drop is sound (or
  //!        nothing was dropped).  Lets a driver loop setup()/solve without ever
  //!        naming the internal claim-key type.  Typical use:
  //!            oc.options.INTERFACE.DROP_POLICY = Options::DROP_REDERIVE;
  //!            oc.clear_protected_interface_claims();
  //!            for( int k=0;; ++k ){
  //!              oc.setup(); /* re-size+init xv at the new layout */ solve(xv);
  //!              auto st = oc.verify_and_protect( xv.data() );
  //!              if( st != InterfaceDropStatus::OVERDROP || k>=CAP ) break;
  //!            }
  InterfaceDropStatus verify_and_protect( double const* var, double rel_tol = 1.0e-2 )
    {
      std::set<t_FrozenInterfaceClaimKey> over;
      InterfaceDropStatus const st = verify_interface_drop( var, rel_tol, &over );
      if( st == InterfaceDropStatus::OVERDROP
       && options.INTERFACE.DROP_POLICY == Options::DROP_REDERIVE
       && !over.empty() )
        protect_interface_claims( over );
      return st;
    }

  //! @brief Evaluate a collocated state/input variable at a physical point. The point map is keyed by physical domain variables.  For distributed variables, the method locates the finite-element block containing pt and evaluates the tensor-product Lagrange interpolant from the flat state/input coefficient arrays.  Lumped variables are returned directly.
  template <typename U>
  U eval_colloc
    ( FFVar const& Var, t_Coord const& pt,
      U const* val_state, U const* val_input=nullptr, U const* val_cst=nullptr,
      t_Side const* side=nullptr )
    const;

  //! @brief Vector wrapper for eval_colloc().
  template <typename U>
  U eval_colloc
    ( FFVar const& Var, t_Coord const& pt,
      std::vector<U> const& val_state, std::vector<U> const& val_input=std::vector<U>(),
      std::vector<U> const& val_cst=std::vector<U>() )
    const
    {
      return eval_colloc( Var, pt,
        val_state.empty()? nullptr: val_state.data(),
        val_input.empty()? nullptr: val_input.data(),
        val_cst.empty()? nullptr: val_cst.data() );
    }

  //! @brief Evaluate the ORD-th physical derivative of a collocated state/input along a
  //! single domain direction `wrt` at a physical point.  Mirrors eval_colloc() but applies
  //! the order-ORD Lagrange differential along `wrt` (with the exact per-element physical
  //! rescale (2/elem_len[iel])^ORD) and order 0 along every other direction.  ORD<=0 is
  //! identical to eval_colloc(); the derivative of a variable that does not depend on `wrt`
  //! (or of a lumped variable) is 0.
  template <typename U>
  U eval_colloc_deriv
    ( FFVar const& Var, t_Coord const& pt, FFVar const& wrt, int ORD,
      U const* val_state, U const* val_input=nullptr, U const* val_cst=nullptr,
      t_Side const* side=nullptr )
    const;

  //! @brief Select the equation-row node donated to a strong exact-interface
  //! continuity row on a one-dimensional normal collocation stencil.
  //!
  //! The returned node may differ from the physical interface node.  In
  //! particular, endpoint-based Lobatto/Chebyshev-Lobatto stencils can become
  //! locally rank-deficient when endpoint derivative rows on both sides of an
  //! element are replaced by exact value rows; this helper shifts the donor to
  //! a nearest interior node when the local rank check requires it.
  size_t select_strong_interface_donor_inode
    ( FFDom const& dom, size_t iel, size_t interface_inode, int lim,
      bool has_left_exact_value, bool has_right_exact_value )
    const;

  //! @brief Initialize collocated state and input coefficient arrays from stored reference data.
  //! When supplied, cst provides the current constant values and refreshes the
  //! values stored in this OCFESLV before auxiliary references are evaluated.
  bool init
    ( double* var, double* inp=nullptr, double const* cst=nullptr );

  //! @brief Convenience wrapper for vector-based initialization.
  bool init
    ( std::vector<double>& var, std::vector<double>& inp, double const* cst=nullptr )
    {
      var.assign( _nCollVar, 0. );
      inp.assign( _nCollInp, 0. );
      return init( var.empty()? nullptr: var.data(), inp.empty()? nullptr: inp.data(), cst );
    }

  //! @brief Convenience wrapper for vector-based initialization with constant values.
  bool init
    ( std::vector<double>& var, std::vector<double>& inp, std::vector<double> const& cst )
    {
      var.assign( _nCollVar, 0. );
      inp.assign( _nCollInp, 0. );
      if( !cst.empty() && cst.size() != _vCst.size() ) return false;
      return init( var.empty()? nullptr: var.data(),
                   inp.empty()? nullptr: inp.data(),
                   cst.empty()? nullptr: cst.data() );
    }


  //! @brief Compute and cache the PDE classification using stored/default references. Uses the user-specified evolution domain, or infers one from INITIAL equations or domain names. State references default to zero, domain references to midpoints, and constant references to set_constant values when available.
  bool classify_pde
    ()
    {
      if( !_classify_pde() ) return false;
      // Classification changes IC_AUTO dispatch and therefore the frozen
      // interface/SAT plan.  Do not rebuild or audit the plan outside setup();
      // require the caller to run setup() again before evaluation.
      if( _issetup ){
        _clear_frozen_interface_consumers();
        _issetup = false;
      }
      return true;
    }
    
  // ======================================================================
  // Types for principal symbol and PDE classification -
  // declared here (before protected:) so they are visible to member
  // variable declarations that follow.
  // ======================================================================


protected:


  //! @brief Traces of all states in one classified block on an internal face.
  //!
  //! The helper stores the left-element trace at +1 and the right-element
  //! trace at -1 in the selected face direction.  The same object is used by
  //! strong IC_UPWIND rows and weak IC_UPWIND SAT penalties so that both paths
  //! use the same block, state ordering and characteristic face matrices.
  template <typename U>
  struct t_InterfaceTraceData
  {
    t_Symbol const*             sym   = nullptr;
    t_FaceData const*           face  = nullptr;
    size_t                      nface = 0;
    std::map<long,size_t>       state_index;
    std::map<long,size_t>       eqn_index;
    std::vector<size_t>         sym_to_vvar;
    std::vector<bool>           active;
    std::vector< std::vector<U> > u_left;
    std::vector< std::vector<U> > u_right;
  };


  //! @brief Recorded evolution-element boundaries (the march grid) when SOLVE.MARCHING
  //! collapsed a multi-element evolution declaration to a single-element template at
  //! setup().  Size N+1 for an N-step march; empty for a monolithic setup.
  std::vector<double>                       _marchGrid;

  //! @brief SOLVE.MARCHING transfer pairs (state, IC-input): between elements each
  //! state's terminal is written into its IC input.  One per evolution-bearing state (a
  //! system marches all of them).  User vars are stored directly (VAR ids are preserved
  //! into the local DAG, so terminal_profile/transfer_terminal look them up).
  std::vector<std::pair<FFVar,FFVar>>       _marchTransfers;
  //! @brief MARCHING transitions (add_transition): at a window seam tau the post-jump states x+ solve
  //! R( x-, x+, p ) = left(x-) - right(x+) = 0, with x- the previous window's terminals, and x+ is handed on instead
  //! of the terminal (see _apply_march_jumps).  Built once by _build_march_jumps.
  struct t_MarchJump
  {
    double             tau = 0.;
    std::vector<FFVar> evo;        //!< every evolving state (its terminal is x-)
    std::vector<FFVar> xs;         //!< the post-jump states S+, in matching order
    std::vector<FFVar> xp;         //!< placeholders standing for x+ of xs
    std::vector<FFVar> R;          //!< the residuals, in terms of evo (at tau^-), xp and the inputs
    // the FAD Jacobians (FFGraph::FAD allocates them with new[]): owned here, freed with the struct (2026-10-03:
    // a solver that marched across a transition leaked them at destruction -- valgrind on OCFE_TRANSrep)
    std::unique_ptr<FFVar[]> J;    //!< dR/dxp, m x m
    std::unique_ptr<FFVar[]> Jx;   //!< dR/devo, m x |evo|
    std::unique_ptr<FFVar[]> Jp;   //!< dR/dinp, m x |inp| (null if no inputs)
    std::vector<FFVar> inp;        //!< time-invariant inputs read by R
    std::vector<FFVar> vars;       //!< evaluation list: evo, xp, inp, constants
    std::vector<double> at;        //!< its values at the last solved jump (x-, x+, inputs, constants) -- node 0
    std::vector< std::vector<double> > at_n;   //!< ... at every spatial node of the terminal profile (v1.5)
  };
  std::vector<t_MarchJump> _marchJumps;
  bool                     _marchJumpsReady = false;


  //! @brief Marching general-IC continuity override (for transfer states designated with the
  //! override path, i.e. null ic_input).  _marchICRows maps each such state
  //! to its INITIAL-closure residual rows and the state's LB-node column in each row (resolved once
  //! from _rowEqnNdx + _structural_dae_decomposition, guarded by _marchICRowsReady).  During the
  //! march, _marchICTerminal holds the current window's terminal per state and _marchICActive gates
  //! the in-place override (residual = var[lb_col] - terminal, Jacobian row = unit on lb_col) that
  //! solve() applies at k>0.
  struct t_MarchICRow { size_t row = 0; size_t lb_col = 0; };
  std::map<FFVar,std::vector<t_MarchICRow>,lt_FFVar> _marchICRows;
  std::map<FFVar,std::vector<double>,lt_FFVar>       _marchICTerminal;
  //! For each override state: LB-face DOF column -> spatial-node index k (the terminal_profile order).
  //! Lets the solve-time override read each INITIAL row's continuity column from the row's own
  //! Jacobian nonzeros (robust to row-emission vs var_dom ordering) and index _marchICTerminal by node.
  std::map<FFVar,std::map<size_t,size_t>,lt_FFVar>   _marchICColToNode;
  bool                                      _marchICActive    = false;
  bool                                      _marchICRowsReady = false;

  //! @brief Whether a non-integral output (point/terminal/spatial) should be captured (overwritten)
  //! in the marching window [t0,t1]: true iff its evolution-t coordinate lies in the window -- so an
  //! initial-time output (t=0) is captured at window 0 and a final-time output at the last window --
  //! or the output has no evolution-t coordinate (spatial-only -> keep the last window's value).
  //! Evolution-direction integrals are summed and do not consult this.
  //! @brief Is output f DISTRIBUTED over the evolution direction?  (2026-10-06, WORKPLAN 1.1)  Such an output spans
  //! the whole horizon: under marching every window emits the rows of ITS time element (the emission applies the
  //! evolution mask with the global element index), and _march_field_rows places them in the monolithic layout.
  bool _march_field_over_evo( t_Fct const& f ) const override
  {
    if( f.kind != FctKind::DISTRIBUTED || !_evolution_dom_var.dag() ) return false;
    for( auto const& [dv, lim] : f.grid ){ (void)lim; if( dv.id() == _evolution_dom_var.id() ) return true; }
    return false;
  }
  //! @brief Nodes of one element kept by a grid mask (the emission's rule: the domain's ends by mask, the unique-node
  //! policy for the left trace of every non-first LGL/CGL element; LB/UB are one interpolated row per block).
  size_t _grid_sel_count( FFDom const& dom, int lim, size_t iel, size_t nel ) const
  {
    if( lim == FFDom::LB || lim == FFDom::UB ) return 1;
    size_t n = 0;
    for( size_t k = 0; k < dom.n_node; ++k ){
      if( ( lim == FFDom::ALL-FFDom::LB || lim == FFDom::ALL-FFDom::LB-FFDom::UB ) && iel == 0 && k == 0 ) continue;
      if( ( lim == FFDom::ALL-FFDom::UB || lim == FFDom::ALL-FFDom::LB-FFDom::UB ) && iel + 1 == nel && k + 1 == dom.n_node ) continue;
      if( options.OUTPUT.GRID_POLICY == Options::OUT_UNIQUE_PHYSICAL && ( dom.type == FFDom::LGL || dom.type == FFDom::CGL )
       && iel > 0 && k == 0 ) continue;
      ++n;
    }
    return n;
  }
  //! @brief For output f distributed over the evolution direction and marching window k: for each row the WINDOW's
  //! emission writes (in its order: the full grid's blocks restricted to time element k, in _next_block order --
  //! first grid variable fastest, LB/UB pinned -- each block laid out dimension by dimension, first dimension
  //! fastest, with ALL the window's evolution nodes), the GLOBAL row it occupies in the monolithic layout (relative
  //! to f.row0), or npos when the global evolution mask drops that node (t = 0 under NLB, the final node under NUB,
  //! the left trace of a non-first LGL/CGL element under OUT_UNIQUE_PHYSICAL) or when the window does not own the
  //! row (an LB/UB evolution mask: the first/last window does).
  static constexpr size_t npos = static_cast<size_t>( -1 );
  std::vector<size_t> _march_field_rows( t_Fct const& f, size_t k ) const override
  {
    std::vector<size_t> rows;
    size_t const nw = _marchGrid.size() > 1? _marchGrid.size() - 1: 1;
    std::vector<FFVar> dims;  std::vector<int> lims;  std::vector<size_t> nel, cur;
    for( auto const& [dv, lim] : f.grid ){
      dims.push_back( dv );  lims.push_back( lim );
      nel.push_back( dv.id() == _evolution_dom_var.id()? nw: _mDom.at( dv ).n_elem );
      cur.push_back( lim == FFDom::UB? nel.back() - 1: 0 );
    }
    size_t ievo = dims.size();
    for( size_t i = 0; i < dims.size(); ++i ) if( dims[i].id() == _evolution_dom_var.id() ) ievo = i;
    if( ievo == dims.size() ) return rows;
    FFDom const& tdom = _mDom.at( dims[ievo] );
    bool const evo_pinned = ( lims[ievo] == FFDom::LB || lims[ievo] == FFDom::UB );
    // the evolution node list of global element e: kept node -> compressed index, dropped -> npos
    auto evo_kept = [&]( size_t e ){
      std::vector<size_t> m( tdom.n_node, npos );  size_t c = 0;
      for( size_t q = 0; q < tdom.n_node; ++q ){
        int const lim = lims[ievo];
        if( ( lim == FFDom::ALL-FFDom::LB || lim == FFDom::ALL-FFDom::LB-FFDom::UB ) && e == 0 && q == 0 ) continue;
        if( ( lim == FFDom::ALL-FFDom::UB || lim == FFDom::ALL-FFDom::LB-FFDom::UB ) && e + 1 == nw && q + 1 == tdom.n_node ) continue;
        if( options.OUTPUT.GRID_POLICY == Options::OUT_UNIQUE_PHYSICAL && ( tdom.type == FFDom::LGL || tdom.type == FFDom::CGL ) && e > 0 && q == 0 ) continue;
        m[q] = c++;
      }
      return m; };
    size_t off = 0;
    while( true ){
      // the global block's per-dimension selection counts, and the window block's (all evolution nodes)
      std::vector<size_t> gsel( dims.size() ), wsel( dims.size() );
      for( size_t i = 0; i < dims.size(); ++i ){
        gsel[i] = _grid_sel_count( _mDom.at( dims[i] ), lims[i], cur[i], nel[i] );
        wsel[i] = ( i == ievo && !evo_pinned )? tdom.n_node: gsel[i];
      }
      size_t gsize = 1, wsize = 1;
      for( size_t i = 0; i < dims.size(); ++i ){ gsize *= gsel[i]; wsize *= wsel[i]; }
      bool const owns = ( cur[ievo] == k );
      if( owns || evo_pinned ){                       // the window visits this block
        std::vector<size_t> const m = evo_pinned? std::vector<size_t>{ 0 }: evo_kept( cur[ievo] );
        for( size_t w = 0; w < wsize; ++w ){
          // decompose w (first dimension fastest) and recompose the global in-block index
          size_t rem = w, g = 0, gstride = 1;  bool drop = !owns;
          for( size_t i = 0; i < dims.size(); ++i ){
            size_t const idx = rem % wsel[i];  rem /= wsel[i];
            size_t gidx = idx;
            if( i == ievo && !evo_pinned ){ gidx = m[idx]; if( gidx == npos ) drop = true; }
            g += ( drop? 0: gidx ) * gstride;  gstride *= gsel[i];
          }
          rows.push_back( drop? npos: off + g );
        }
      }
      off += gsize;
      size_t i = 0;
      for( ; i < dims.size(); ++i ){
        if( lims[i] == FFDom::LB || lims[i] == FFDom::UB ) continue;
        if( ++cur[i] < nel[i] ) break;
        cur[i] = 0;
      }
      if( i == dims.size() ) break;
    }
    return rows;
  }
  bool _march_point_capture( t_Fct const& f, double t0, double t1 ) const
  {
    if( !_evolution_dom_var.dag() ) return true;
    auto it = f.point.find( _evolution_dom_var );
    if( it == f.point.end() ) return true;
    double const eps = 1e-9 * std::max( 1.0, std::fabs( t1 - t0 ) );
    // Own each point output to the SINGLE window where its t is the UB -- the window that integrated up
    // to it, whose block Jacobian for that output is correctly formed.  LB-exclusive, so an interior seam
    // point t_m is claimed only by the window ending at t_m and NOT also by the successor starting at t_m,
    // whose fresh block contributes zero for that row and (under last-wins / the adjoint's highest-k latch)
    // would clobber the correct value -- the interior-seam reduced-Jacobian-row zeroing.  The FIRST window
    // keeps its LB so a point at the global start (t == _marchGrid.front()) still lands.
    bool const first = _marchGrid.empty() || std::fabs( t0 - _marchGrid.front() ) <= eps;
    // 4.10 (2026-10-09): a PLUS-sided point output (add_output(f,{t},{tau},{FFDom::PLUS})) is the RIGHT limit at tau:
    // at an interior seam only the window STARTING there holds it -- owned with the LB inclusive and the UB exclusive
    // (the last window keeps a point at the global end).  The rule above gave it the LEFT window: u(tau+) read u(tau-).
    auto const isd = f.side.find( _evolution_dom_var );
    if( isd != f.side.end() && isd->second == FFDom::PLUS ){
      bool const last = _marchGrid.empty() || std::fabs( t1 - _marchGrid.back() ) <= eps;
      return ( it->second >= t0 - eps ) && ( last ? ( it->second <= t1 + eps ) : ( it->second < t1 - eps ) );
    }
    bool const lo_ok = first ? ( it->second >= t0 - eps ) : ( it->second > t0 + eps );
    return lo_ok && ( it->second <= t1 + eps );
  }
  //! @brief Per-element-accumulated output-function values from the most recent
  //! SOLVE.MARCHING solve: evolution-direction integrals (OpI over the evolution var) are
  //! summed across elements; all other outputs take the last element's value.  Size
  //! n_colloc_fct(); empty if the model has no outputs or the last solve was monolithic.
  std::vector<double>                       _valFct;

  //! @brief Marching captured-reduction cache (Approach A): at each converged SOLVE.MARCHING solve the
  //! accumulated input array (with cap[] filled per window) and the last-window var[] are stored here, so
  //! a subsequent DIRECT eval() of a captured output on a marching model can serve the completed caps
  //! (which a last-window var[] alone cannot reconstruct) instead of the (monolithic-only) self-fill.
  //! Served only when the eval's var[] matches _marchCapVar (the solved point); a copy the size of the
  //! input / state vectors, refreshed once per marching solve -- negligible vs the solve.
  std::vector<double>                       _marchCapInp;   //!< accumulated inp[] (caps filled) at the solve
  std::vector<double>                       _marchCapVar;   //!< last-window var[] at the solve (scope key)
  bool                                      _marchCapValid = false;

  //! @brief Primal-solve reuse cache (SOLVE.REUSE): the last converged MONOLITHIC solution and
  //! the input it was solved at.  A subsequent monolithic solve() at the same input returns this
  //! directly, so the gradient path (_sens_setup -> _mono_fsens_setup -> solve) skips the redundant
  //! primal Newton and factorizes at the value path's converged point.  Not used for marching or
  //! nested (_inMarch) solves.
  std::vector<double>                       _primalInp;
  std::vector<double>                       _primalVar;
  bool                                      _primalValid = false;


  //! @brief Per-element block solutions captured when OUTPUT.MARCH_STORE is set.
  std::vector<MarchBlock>                   _marchTrajectory;

  //! @brief Per-element linearization captured by _march_fsens_setup() and reused by
  //! _march_fsens_apply(): primal block solution/inputs, window, the LU of J_k=dR_k/dvar_k
  //! (P*J=L*U), and the sparse dR/dinp & dF Jacobians (per-row columns+values).
  struct FSensBlock {
    double                            t0 = 0., t1 = 0.;
    std::vector<double>               var, inp;
    arma::mat                         L, U, P;
    //! sparse Jacobian, stored INSTEAD of the dense LU factors when SENSITIVITY.SPARSE is on.
    //!
    //! DESIGN POINT, recorded because it is the one thing here that can go the wrong way.
    //! The dense path factorises ONCE in setup and every forward/adjoint apply reuses L/U/P.
    //! arma::spsolve exposes no reusable factorisation object, so this path REFACTORISES PER
    //! APPLY -- one solve per control (forward) and per output (adjoint).  With few controls
    //! that is far cheaper than a dense O(n^3) factorisation of an 18202^2 matrix; with many
    //! it could lose.  Option (c) in the scope note -- SuperLU's own factor/solve API with a
    //! stored handle -- is the answer if it does, and SPQR (tier 1) needs that machinery
    //! anyway, so building it now would duplicate work tier 1 subsumes.
    //! MEASURE _ts_sens BEFORE concluding either way.
    arma::sp_mat                      Jsp;
    bool                              sparse = false;
    //! when true, L holds the dense pseudo-inverse
    //! J^+ and U/P are unused.  The forward apply is then s = J^+ b and the adjoint is
    //! lambda = (J^+)^T ls, so both apply sites branch on this ONE flag rather than
    //! acquiring a parallel code path.
    bool                              pinv = false;
    std::vector<std::vector<size_t>>  ceqn, cfct;
    std::vector<std::vector<double>>  geqn, gfct;
    //! General-IC override continuity rows in THIS block (empty for k=0 and canonical windows):
    //! (residual row, transfer state, LB-node index into that state's terminal profile).  The
    //! sensitivity RHS (forward) / adjoint seed for these rows carries the previous window's
    //! terminal sensitivity, mirroring the canonical ic_input coupling.
    std::vector<std::tuple<size_t,FFVar,size_t>> icRows;
  };
  std::vector<FSensBlock>                   _fsBlocks;

  //! @brief Monolithic (single-block) sensitivity linearization captured by
  //! _mono_fsens_setup() and reused by _mono_fsens_apply()/_mono_asens_apply(): the LU of the
  //! full J=dR/dvar, the sparse dR/dinp & dF Jacobians, and the primal function values.
  //! This is the N=1, no-terminal-coupling counterpart of _fsBlocks for systems that are
  //! not marched (e.g. boundary-value / genuinely coupled space-time models).
  FSensBlock                                _monoBlock;
  bool                                      _monoValid = false;
  std::vector<double>                       _monoFct;

  //! @brief Function values and reduced Jacobian dF/dp (row-major, n_colloc_fct x n_control_dof)
  //! captured by the last complete solve_fsens()/solve_asens(); sens_functions() is also refreshed
  //! by _sens_setup() (marching accumulates per-window, monolithic caches the single-block eval).
  std::vector<double>                       _sensFct;
  std::vector<double>                       _sensJac;

  // The control registry (_controls, _nControlDof) is FFModel's -- see the note by
  // _on_controls_changed() above.  The uses below resolve to the inherited members unchanged.

  //! @brief OCFEBASE::_march_context(), answered from the marching members below.
  MarchContext _march_context
    () const override
    {
      MarchContext m;
      m.in_march = _inMarch;  m.window = _marchWindow;
      m.n_grid   = _marchGrid.size();  m.n_steps = n_march_steps();
      m.active   = options.SOLVE.MARCHING && m.n_steps >= 1;
      return m;
    }

  //! @brief Reentrancy guard: true while _solve_marching() drives per-element solve()s.
  bool                                      _inMarch = false;
  size_t                                    _marchWindow = 0;      //!< resident window index during a marching solve; a distributed input matched to the evolution direction remaps its evolution element to this window (see the INPUT cases in the assembly)
  bool                                      _lastSolveMarched = false; //!< last solve() marched with a stored trajectory -> route eval_colloc reads through it

  //! @brief SOLVE.MARCHING setup hook: record the evolution grid and collapse the
  //! evolution domain in _mDom to a single element (decoupled marching template).  A
  //! no-op that clears _marchGrid when SOLVE.MARCHING is off; fails only when marching
  //! model, since no transfer pair is designated).
  void _collapse_evolution_for_marching();

  //! @brief Resolve the marching transfer set from the DIFFERENTIAL states (single source of truth,
  //! via _states_with_evolution_derivative -- the same rule the DAE decomposition uses for dyn),
  //! classifying each by its INITIAL equation (canonical `state - ic_input` -> data path; any other
  //! closure -> override path).  The set is structurally determined; there is no user override.
  void _autodetect_marching_transfers();


  // rev92 PATH B probe -- set in setup() BEFORE the collapse, so they must live outside
  // _evoProbe, which _prepare_interface_plan_tables resets much later.
  int         _evoHoistLevel       = 0;      //!< CRONOS_HOIST_EVODETECT as parsed
  bool        _evoHoistApplied     = false;  //!< detection was actually hoisted
  bool        _evoHoistWouldChange = false;  //!< a domain WAS inferable before the collapse
  std::string _evoHoistDom;                  //!< which one


  //! @brief Resolve, from the residual-row provenance (_rowEqnNdx) and the DAE decomposition, the
  //! INITIAL-closure rows + LB column for each override-path transfer state (null ic_input).  Builds
  //! _marchICRows; requires _rowProvenanceReady (one eval must have run).  Idempotent via _marchICRowsReady.
  void _build_march_ic_rows();

  //! @brief LB-face DOF columns of override-path state @p st in the collapsed marching template,
  //! one per spatial node, enumerated in the SAME spatial-node order as terminal_profile()
  //! (spatial domains in var_dom order, element-major/node-inner, first spatial domain fastest).
  //! For a time-only state (no spatial domains) this returns exactly one column == the old
  //! pos_state(st,{all->0}) corner, so the time-only path is byte-for-byte unchanged.  Used by
  //! _build_march_ic_rows() to give each distributed INITIAL row its own continuity column.
  std::vector<size_t> _march_ic_lb_columns( FFVar const& st ) const;

  //! @brief evolution-direction TRACE of the window lower bound, for the marching
  //! IC pin.  Returns the DOF stride between consecutive evolution nodes of the LB element
  //! and the Lagrange weights of that element's nodes AT the lower bound.
  //! _march_ic_lb_columns() returns the column with evolution node index 0, on the assumption
  //! that node 0 SITS at the lower bound.  That holds for LGL and LGR (both include the left
  //! endpoint) and FAILS for LG, whose first node is strictly interior -- so the pin lands on
  //! an interior coefficient and every window after the first is handed a wrong IC.  MEASURED
  //! (OCFE_PDE25, 2026-09-08): marched LG converges to machine-zero residual in every window
  //! and is 6.1e-3 wrong, deterministically; monolithic LG is exact to 2.2e-16.
  //! With (stride, w) the pin column c generalises to the trace {c + j*stride, w[j]}, which
  //! reduces to {c, 1} on LGL/LGR by construction.
  bool _march_ic_lb_trace( FFVar const& st, size_t& stride_evo, std::vector<double>& w ) const;


  //! @brief Fill captured-reduction INPUTS from a converged solution, writing into @p inp.
  //! ACCUM: cap += INT R over the CURRENT evolution domain (an explicit LGR element quadrature
  //! summed over its elements -- one collapsed window under marching, the full horizon under a
  //! monolithic solve).  LATCH: cap = R(tau).  @p k is the march-window index, or -1 for a
  //! monolithic (single, overwriting) fill; @p cap_window schedules LATCH writes under marching.
  void _fill_captures( double const* var, double* inp, double const* cst,
                       int k, std::vector<int> const& cap_window ) const;

  //! @brief Capture-SENSITIVITY fill (Stage 1): the exact LGR quadrature of _fill_captures run on the
  //! STATE-SENSITIVITY buffer @p s_var instead of the primal var[], writing d(cap)/dp into the cap
  //! slots of @p dcap.  Because each capture source is a bare (materialised) state whose forward
  //! sensitivity s[src] IS d(src)/dp, the window quadrature of s[src] is d(cap)/dp.  ACCUM accumulates
  //! per window (k>=0 read-modify-write of dcap; k<0 single overwrite); LATCH writes d(R(tau))/dp at
  //! tau's window.  @p dcap is indexed like the input array; only cap slots are touched.
  void _fill_capture_sens( double const* s_var, int k,
                           std::vector<int> const& cap_window, double* dcap ) const;


  //! @brief Evaluate output functions into _valFct, first filling captured-reduction inputs over
  //! the (full-domain) solution on a mutable copy of @p inp -- the monolithic analogue of the
  //! per-window capture+accumulate performed inside _solve_marching.
  void _eval_functions_captured( double const* var, double const* inp, double const* cst );


  //! @brief W-test every dropped claim on the built system; protect the
  //! needed subset and re-derive once (revert-guarded).  Returns true iff a
  //! restore re-derive ran (the inner setup already audited and printed).
  //! most recent audit record for (marching, window), or nullptr.  Reverse
  //! scan because _determinacy appends per setup(); the newest record for the window
  //! belongs to the plan being solved.  A stale hit degrades safely: a wrong "attempt
  //! tier-0" fails one factorization and falls through via nt_dead.
  t_Determinacy const* _nt_skip_verdict( bool marching, size_t win ) const
    { for( auto it = _determinacy.rbegin(); it != _determinacy.rend(); ++it )
        if( it->marching == marching && it->window == win ) return &*it;
      return nullptr; }

  bool _wdecide_restore();


  //! @brief Helper for copying optional input discretisation maps
  static t_Coll _copy_inp_disc
    ( t_Coll const& in, std::vector<FFVar> const& all,
      std::map<FFVar::pt_idVar,size_t> const& ndx );


  //! @brief Assert that supplied input data honour the per-direction interface-continuity
  //! declarations in _mInpCont.  For each declared (input, direction, level>=C0) it samples
  //! every interior interface from both sides (eval_colloc +/- eps over a transverse
  //! element-midpoint grid) and returns false on a C0 (or, for SMOOTH_C1, a C1) violation
  //! -- over-declaring smoother-than-the-data is the dangerous direction.  Inert (returns
  //! true immediately) when _mInpCont is empty.
  bool _check_input_continuity
    ( double const* var, double const* inp, double const* cst,
      double tol = 1e-6 ) const;


  //! @brief Count nonzero entries in a sparse row pattern
  static size_t _n_jac_nnz
    ( std::vector<std::vector<size_t>> const& pattern )
    {
      size_t nnz = 0;
      for( auto const& row : pattern ) nnz += row.size();
      return nnz;
    }

public:

  //! @brief Supply an extra STATE reference for the reference-robustness probe
  //! (MC__OCFESLV_REFERENCE_ROBUSTNESS_PROBE), as a partial map {state -> value}.
  //! Listed states override the base reference; unlisted states keep it.  Supply
  //! states that span any anticipated flow reversal so the probe can catch a
  //! characteristic sign flip (e.g. a PSA velocity reversing across cycle steps).
  //! Without any sample the probe falls back to a LOCAL magnitude sweep of the
  //! single reference (sensitivity only, no sign-flip coverage).  Inert unless the
  //! probe macro is defined.
  void add_classification_reference_sample
    ( std::map<FFVar,double,lt_FFVar> const& state_vals )
    { _classRefSamples.push_back( state_vals ); }

  //! @brief Verdict of the reference-robustness probe after setup(): true if the
  //! PDE type and per-face characteristic split were invariant across all sampled
  //! states.  Always true when the probe is not compiled in.
  bool reference_robustness_ok() const { return _refRobustOk; }

  //! @brief Identity of THIS header revision, stamped into the log at first setup().
  //!
  //! A stale object file that was not relinked against a new header produces a log which
  //! is indistinguishable from a correct run of the OLD header -- the new diagnostics
  //! simply do not appear, which reads as "the instrument found nothing" rather than "the
  //! instrument is not there".  That has cost this project a session before.  Bump the
  //! string with every header revision; a log without the expected id is a stale binary,
  //! full stop, and no conclusion should be drawn from it.
  static constexpr char const* HEADER_ID
    // BUMPED WITH THE FILE, ALWAYS.  This string is the header's answer to "which
    // revision produced this log", and it read rev106 while the makefile was building
    // rev110 -- so a solve log carried a revision three steps stale.  Nothing else in
    // the build chain would have caught it: the makefile names the file, the build
    // oracle checks the instrument's format marker, and neither reads this.
    = "ocfeslv  rev373  2026-10-11";

  //! @brief The revision of this ocfeslv.hpp (HEADER_ID); FFModel::revision() gives the model layer's, which a
  //! binary may mix with another solver revision.
  static char const* revision() { return HEADER_ID; }

  //! @brief One-line summary of what this revision changed, for the process stamp.
  //!
  //! Why a summary and not a changelog: a stamp carrying the cumulative changelog (4.5 kB, printed
  //! once per process and on every [cov1] coverage report) was MEASURED on the corpus logs at 70%
  //! of OCFE_AE0's total output and 42% of OCFE_DAE0's.
  //!
  //! The stamp exists as a BUILD ORACLE: "a log without the expected id is a stale
  //! binary, full stop".  A revision tag serves that purpose exactly; the changelog does
  //! not, and belongs in the header comments where it can be read rather than in every
  //! log where it must be skipped.  Sweep scripts derive EXPECT_HDR from this same
  //! constant, so a short id keeps working -- rank_sweep7.sh greps the tag, not the prose.
  //!
  //! Set CRONOS_HEADER_FULL=1 to append the summary line to the stamp.
  static constexpr char const* HEADER_SUMMARY
    = "rev99: PATH B PROMOTED -- CRONOS_HOIST_EVODETECT now DEFAULTS TO 1, so the evolution"
      " domain is inferred BEFORE the localisation collapse and an auto-detect model gets the"
      " same collapsed plan a user-supplied one gets.  Measured on the 139-target corpus: 10"
      " drivers applied, 7 march, 3 decline (empty transfer set); k 45;105;225;465 -> 3;7;15;31"
      " on MMPDE27, 3;7;15 -> 0 on DAE5, 36 -> 18 on PDE9; two non-convergences fixed; systems"
      " 4-17x smaller; adjoint-vs-forward sensitivity intact at round-off on both drivers that"
      " change solver path; 14 drivers index var post-solve without a coordinate evaluator and"
      " NONE of them march.  ALSO records the INFERENCE PROVENANCE: rule 1 finds the domain an"
      " INITIAL row pins (structural), rule 2 matches a domain by NAME (t/T/time/Time, a"
      " convention).  Both set the same member and rev98 reported them identically; at the"
      " promoted default the inference decides whether the solve MARCHES, and a wrong causal"
      " direction is silently wrong where monolithic is merely slow, so the hoist probe now"
      " prints infer=INITIAL-pinned or infer=NAME-BASED (weak).  CRONOS_HOIST_EVODETECT=0"
      " remains the exact opt-out.  "
      "rev96: DELTA (A) PROMOTED -- CRONOS_CAUSAL_EVOFACE now DEFAULTS TO 1, so the causal"
      " evolution-direction block is live under its original gate (evolution_input_may_jump)."
      "  Measured on the 139-target corpus: 7 affected drivers, 0 verdict regressions, 1"
      " improvement (DAE3 58/2 -> 60/0), and 30 numeric checks better with NONE worse across"
      " the three that report quantities (DAE2 13/0, DAE3 12/0, DAE4 5/0) -- including"
      " adjoint-vs-forward Jacobian agreement improving to round-off on all three."
      "  `acted` is 0 on all 135 usable drivers at level 0, so CRONOS_CAUSAL_EVOFACE=0 is a"
      " verified exact opt-out to rev86 behaviour.  rev95: fixes _differential_core_shape being called before _evolution_dom_var is"
      " assigned, which made it return 0x0 on every auto-detect model -- core=0x0 on all"
      " 139 corpus drivers, and level 2's square-core gate had no data.  The candidate"
      " domain is now set before the measurement and restored if the hoist does not apply,"
      " so levels 0 and refused-2 stay inert.  rev94: CRONOS_HOIST_EVODETECT=2 hoists only when the differential core is SQUARE,"
      " replicating rev94:9096's marchability rule at the collapse site because the"
      " classifier runs after it.  That rule's comment argues against Path B on the grounds"
      " that auto-detected causality is unverified; its EVIDENCE (MMPDE27 auto failing) is"
      " refuted by rev92 as the ordering accident, its CONCERN is not.  Reports the"
      " user-side core shape, which may DISAGREE with the [classify] sym=RxC because it"
      " runs pre-localisation and pre-reduction -- both are printed."
      "  rev93: unconditional-under-EVOFACE marching stamp.  The only prior evidence that a"
      " model marched was gated on options.SOLVE.VERBOSE, so `grep -c \"march] step\"`"
      " measured SOLVE.VERBOSE and not marching -- and handoff rev25 s0.3 used that grep on"
      " OCFE_MMPDE27, which disables SOLVE.VERBOSE by default, to conclude marching was not"
      " involved.  That conclusion is RETRACTED.  rev92: PATH B probe.  CRONOS_HOIST_EVODETECT=1 (rev99 default, apply) / 0 (opt-out, report only)"
      " moves evolution-domain detection ahead of _collapse_evolution_for_marching, which"
      " is the ordering accident that keeps auto-detect models from marching."
      "  MEASURED basis: at NEL_T=3 on OCFE_PDE0, IC_TRACE[march] |dC|=3.493e-06 against"
      " IC_TRACE[mono] 4.593e-05 -- marching is fully causal and 13x BETTER exactly where"
      " delta (B) is 26.8x worse, so causality is not what costs accuracy; the seam is."
      "  NOTE: rev91's comment attributing PDE0's 26x loss to the aux branch is WRONG --"
      " measured per-branch, `diff` costs 26.8x and `aux` costs 9%.  (rev91 CRONOS_CAUSAL_BRANCH removed in rev100 -- delta (B) is superseded by Path B, so the branch selector had no remaining consumer.)"
      "  MEASURED on OCFE_PDE0 level 2, per-state z-interface spread: C and T UNCHANGED at"
      " round-off while Daux12_C and Daux56_T degrade ~1e13, so the differential branch is"
      " free on that diagnostic and the aux-flux branch is not.  The auxiliaries' z-claims"
      " all survive the drop, so the dropped t-face RECEIVERS were carrying z-face"
      " continuity rows: receivers do double duty across directions."
      "  rev90: removes a stray second newline in the rev89 [evoface] gate-probe chain."
      "  Cosmetic in itself, but it emitted one BLANK LINE per plan build, which the L0"
      " inertness gate correctly flagged as a difference from rev86 -- the gate caught the"
      " instrument, not the plan.  rev89: predicate (ii) no longer tests options.SOLVE.MARCHING -- OCFE_PDE0's"
      " monolithic configs set it false BY DESIGN, so the predicate refused on the one"
      " model with a published before-number.  What remains is structural causality."
      "  Every report field is now COMPUTED before any return (rev88 printed three"
      " uncomputed zeros on DAE3 mono that read as measurements).  Adds per-direction"
      " RECEIVER EDGE counts and the committed-plan n_trace_var / rank(B) / k, because"
      " MEASURED on DAE5: dropping 6 receivers left the claim census unchanged at 6 --"
      " both receivers of a seam canonicalise to one trace slot, so claims are blind to"
      " the differential branch.  rev88: DANGLING-IF FIX.  The causal evolution-direction block was the swallowed body"
      " of a diagnostic `if` whose braces lived inside #ifdef MC__OCFESLV_PIN_LEAK_DIAG -- a"
      " macro defined in NO build -- so it has been UNREACHABLE corpus-wide, on every model,"
      " regardless of evolution_input_may_jump.  The #ifdef now encloses the whole statement."
      "  Reviving it is the same act as fixing it, so revival is behind CRONOS_CAUSAL_EVOFACE"
      " (0=rev86-identical default, 1=original gate, 2=predicate (ii)).  Adds an UNGATED"
      " visits-by-direction liveness tally, because rev87's counters were themselves the"
      " swallowed statement and their zeros proved nothing.  rev87: PROBE ONLY -- predicate"
      " (ii) computed alongside the live gate, plus the per-direction claim census that"
      " MEASURED 606 = 276(t) + 330(xi) on MMPDE27 against 81 all-xi collapsed."
      "  rev87 note: the causal evolution-direction interface"
      " block already implements handoff rev25 section 4.1 ('make the monolithic solve"
      " reproduce the marching solve in the evolution direction') but is gated on"
      " evolution_input_may_jump, a predicate about DISTRIBUTED INPUTS -- false by"
      " construction on DAE5/MMPDE27/MESHONLY10/PDE0 (add_input=0), live only on DAE3."
      "  rev87 computes predicate (ii) 'marching would have been available' alongside it,"
      " counts the receiver visits it WOULD drop split by branch, prints a per-direction"
      " claim census, and then takes the rev86 branch regardless.  CRONOS_AUDIT_EVOFACE=1."
      "  rev68: classification stamp carries its EVIDENCE (sym=RxC, nDom, nCoeff, Ae_singular,"
      " sigma_min/max_evol, imag_tol, evolution user-supplied vs auto-detected) -- a verdict"
      " without evidence cannot be checked, and OCFE_DAE0 reported DIFFERENTIAL_ORDINARY for"
      " an unambiguous index-1 DAE with nothing in the log to say why.  ALGEBRAIC_* demoted"
      " from evolution ISSUE to NOTE; classify stamp moved out of the evolution validator"
      " (which early-returns with no evolution domain, so OCFE_AE0 printed nothing).  rev67:"
      " CLASSIFICATION: ALGEBRAIC_FIELD + ALGEBRAIC_LUMPED added; ODE->DIFFERENTIAL_ORDINARY,"
      " DAE->DIFFERENTIAL_ALGEBRAIC.  CRASH FIX: _eval_symbol indexed sym.vCoeff[i] for"
      " i<sym.vDom.size() (DIFFERENT vectors), so a distributed purely-algebraic field"
      " aborted inside setup().  An empty principal symbol is now a CONCLUSIVE verdict."
      "  REPORTING: classification stamped unconditionally once per block; [cov1] no longer"
      " repeats the full header string.";


  //! @brief Print the eN -> equation legend consumed by every audit that reports "by eqn eN".
  //!
  //! The audit's eN is an index into _mEqn, and nothing in the log has ever said which
  //! equation that is -- so "the null couples e6 and e10" was never checkable against the
  //! model.  For each equation this prints its role, block, classification participation,
  //! the states it uses ALGEBRAICALLY (bare, i.e. not under a partial/integral), and every
  //! partial derivative it applies as state'direction^order.  That last column is what
  //! distinguishes an INSTANTANEOUS equidistribution ( d/dxi(M x_xi)=0, no x't anywhere )
  //! from an MMPDE RELAXATION ( tau*x_t = g_xi, x't present ) -- i.e. it answers, from the
  //! assembled DAG rather than from a driver banner, which equation set the binary solves.
  //! eqn -1 in an audit means a plan-owned interface/trace row, not an _mEqn row.
  void print_equation_legend
    ( std::ostream& os = std::cerr )
    const
    {
      os << "  [eqn-legend] " << _mEqn.size() << " equation(s) in _mEqn"
         << " (audit 'eN' indexes this list; e-1 = plan-owned interface/trace row)\n";
      if( !_dag ){ os << "  [eqn-legend] no DAG bound\n"; return; }

      auto role_name = []( EqnRole r ) -> char const* {
        switch( r ){
          case EqnRole::AUTO:       return "AUTO";
          case EqnRole::INTERIOR:   return "INTERIOR";
          case EqnRole::INITIAL:    return "INITIAL";
          case EqnRole::BOUNDARY:   return "BOUNDARY";
          case EqnRole::INTERFACE:  return "INTERFACE";
          case EqnRole::LINK:       return "LINK";
          case EqnRole::SURFACE:    return "SURFACE";
          case EqnRole::DIAGNOSTIC: return "DIAGNOSTIC";
        }
        return "?";
      };

      for( size_t i = 0; i < _mEqn.size(); ++i ){
        auto const& E = _mEqn[i];
        std::set<std::string> bare;
        std::set<std::string> parts;
        try{
          auto sg = _dag->subgraph( 1, &E.var );
          for( auto const& op : sg.l_op ){
            bool const is_par = op->sameid( typeid(FFPartial) );
            bool const is_int = op->sameid( typeid(FFIntegral) );
            if( is_par ){
              auto const* pop = mc::type_cast<FFPartial const>( op );
              for( auto const* od : op->varin ){
                if( !od || _mVar.find( *od ) == _mVar.end() ) continue;
                if( !pop ) continue;
                for( auto const& [iv,ord] : pop->Indep().expr ){
                  std::ostringstream ss;
                  ss << od->name() << "'" << iv.name();
                  if( ord != 1 ) ss << "^" << ord;
                  parts.insert( ss.str() );
                }
              }
              continue;
            }
            if( is_int ) continue;
            for( auto const* in : op->varin )
              if( in && _mVar.find( *in ) != _mVar.end() ) bare.insert( in->name() );
          }
        }
        catch( ... ){ /* diagnostic only: never let the legend abort a solve */ }

        os << "    e" << std::setw(3) << i
           << "  role=" << std::left << std::setw(10) << role_name( E.opt->role ) << std::right
           << " blk=" << E.opt->block_id
           << " classify=" << ( E.opt->participate_in_classification ? "y" : "n" )
           << " sat=" << ( _sopt( *E.opt ).receive_sat ? "y" : "n" );
        // rev125: print the DOMAIN SPEC.  The legend gave role, block, bare states and
        // partials, and never the bound mask -- so on OCFE_PDE6_blk0 the MANUAL and AUTO
        // LINK rows read as structurally IDENTICAL (both role=LINK, bare={aux},
        // partials={T'x}) while their claim counts differ by 100.  If the two are
        // registered over different bound sets, this line says so; nothing else does.
        //
        // Same notation as the plan renderer: [d]=ALL  |d=LB  d|=UB  |d]=ALL-LB
        // [d|=ALL-UB  |d|=interior only.
        os << "  dom=";
        {
          bool dfirst = true;
          for( auto const& kv : E.dom ){
            if( !dfirst ) os << " x ";
            dfirst = false;
            switch( kv.second ){
              case FFDom::ALL:                     os << "[" << kv.first << "]"; break;
              case FFDom::LB:                      os << "|" << kv.first;        break;
              case FFDom::UB:                      os <<        kv.first << "|"; break;
              case FFDom::ALL-FFDom::LB:           os << "|" << kv.first << "]"; break;
              case FFDom::ALL-FFDom::UB:           os << "[" << kv.first << "|"; break;
              case FFDom::ALL-FFDom::LB-FFDom::UB: os << "|" << kv.first << "|"; break;
              default:                             os << "<" << kv.second << ">"; break;
            }
          }
          if( dfirst ) os << "(none)";
        }
        os << "\n           bare={";
        for( auto const& s : bare ) os << " " << s;
        os << " }  partials={";
        for( auto const& s : parts ) os << " " << s;
        os << " }\n";
      }

      if( _evolution_dom_set && _evolution_dom_var.dag() )
        os << "  [eqn-legend] evolution direction = " << _evolution_dom_var.name() << "\n";
    }

  //! @brief Structural coverage of every state column in the FINAL residual system.
  //!
  //! The direct test the pivot report only approximates.  _jacPatternEqn is the structural
  //! pattern of the system as actually assembled -- after Schur pivot rows are dropped and
  //! replaced by trace rows (IC_STRONG), after exact replacements (IC_TRACE), with SAT
  //! contributions included.  A state column that appears in NO row of that pattern is not
  //! determined by anything: it is a free degree of freedom, and it will show up as a right
  //! null vector concentrated on exactly that column.
  //!
  //! This localises the defect at SETUP, with no solve, no SVD and no reasoning about which
  //! pivot took which row -- and it is imposition-mode agnostic, so IC_STRONG, IC_TRACE and
  //! IC_WEAK are directly comparable.
  //!
  //! Necessary, not sufficient: a column covered by rows that are mutually dependent can
  //! still be free.  So uncovered_count <= kernel dimension.  Equality means the kernel is
  //! fully explained by structural coverage; a shortfall means part of it is a rank defect
  //! among covered rows and needs the SVD.  Columns covered exactly ONCE are reported too --
  //! they are one pivot away from becoming free and are the fragility margin.
  void print_column_coverage_report
    ( std::ostream& os = std::cerr )
    const
    {
      // Retained for sweep scripts that key on the [cov1] block.  rev35 stamps HEADER_ID
      // once per process anyway, so this repeat was pure duplication -- and while
      // HEADER_ID carried the full changelog it was 70% of OCFE_AE0's log bytes and 42%
      // of OCFE_DAE0's.  rev67 reduced HEADER_ID to the revision tag, which is what the
      // stamp is FOR: a log without the expected id is a stale binary.
      // With HEADER_ID reduced to the revision tag (rev67), there is nothing left to
      // shorten here: the tag IS the build oracle and sweep scripts grep it.  The
      // CRONOS_COV1_FULL_HDR escape hatch that an earlier rev67 draft added is gone
      // with the long string it was there to suppress.
      os << "  [cov1] hdr=" << HEADER_ID << "\n";
      if( !const_cast<OCFESLV*>( this )->_ensure_deriv_cache() ){
        os << "  [coverage] derivative cache unavailable; skipped\n";
        return;
      }
      size_t const ncol = _nCollVar;
      std::vector<size_t> cover( ncol, 0 );
      for( auto const& row : _jacPatternEqn )
        for( size_t c : row )
          if( c < ncol ) ++cover[c];

      std::map<std::string,size_t> free_by_state, once_by_state;
      std::vector<size_t> free_cols;
      size_t n_free = 0, n_once = 0;
      for( size_t c = 0; c < ncol; ++c ){
        if( cover[c] > 1 ) continue;
        std::string const d = decode_state_column( c );
        size_t const br = d.find( '[' );
        std::string const nm = ( br == std::string::npos ) ? d : d.substr( 0, br );
        if( cover[c] == 0 ){ ++n_free; ++free_by_state[nm]; free_cols.push_back( c ); }
        else               { ++n_once; ++once_by_state[nm]; }
      }

      os << "  [coverage] " << ncol << " state column(s), "
         << _jacPatternEqn.size() << " row(s) in the assembled pattern\n";
      os << "  [coverage] UNCOVERED (free) columns = " << n_free;
      if( n_free ){
        os << "  by state:";
        for( auto const& kv : free_by_state ) os << " " << kv.first << "x" << kv.second;
      }
      os << "\n";
      for( size_t k = 0; k < free_cols.size() && k < 24; ++k )
        os << "      FREE c" << std::setw(5) << free_cols[k] << "   "
           << decode_state_column( free_cols[k] ) << "\n";
      if( free_cols.size() > 24 )
        os << "      ... " << ( free_cols.size() - 24 ) << " more\n";

      os << "  [coverage] singly-covered columns = " << n_once;
      if( n_once ){
        os << "  by state:";
        for( auto const& kv : once_by_state ) os << " " << kv.first << "x" << kv.second;
        os << "   [one pivot from free]";
      }
      os << "\n";
      if( n_free )
        os << "  [coverage] *** compare UNCOVERED against the right-null dimension:"
              " equal => kernel fully structural; fewer => residual rank defect too. ***\n";
      os << "  [coverage] NOTE: structural only.  A column present in the pattern but"
            " NUMERICALLY zero everywhere still reads as covered here; the spectrum audit's"
            " empty-column scan catches those.\n";

      // ---- one-line machine-readable summary, for corpus sweeps ---------------------
      // Fixed prefix, fixed field order, no wrapping: `grep '\[cov1\]' sweep.log` gives one
      // row per setup() across the whole corpus.  CRONOS_AUDIT_TAG labels the run (a sweep
      // script sets it to the driver name), since the header cannot see argv.
      char const* tag = std::getenv( "CRONOS_AUDIT_TAG" );
      os << "  [cov1] tag=" << ( tag && *tag ? tag : "-" )
         << " imp=" << ( options.INTERFACE.IMPOSITION == Options::IC_WEAK   ? "WEAK"
                       : options.INTERFACE.IMPOSITION == Options::IC_TRACE  ? "TRACE"
                       : options.INTERFACE.IMPOSITION == Options::IC_STRONG ? "STRONG" : "?" )
         << " ncol=" << ncol
         << " nrow=" << _jacPatternEqn.size()
         << " free=" << n_free
         << " single=" << n_once
         << " states=";
      if( once_by_state.empty() && free_by_state.empty() ) os << "-";
      else{
        bool first = true;
        for( auto const& kv : once_by_state ){
          if( !first ) os << ",";
          first = false;
          std::string nm = kv.first;
          size_t const p = nm.find( '(' );            // strip "(t,xi)" for a compact field
          if( p != std::string::npos ) nm.erase( p );
          os << nm << ":" << kv.second;
        }
      }
      os << ( n_free ? " EXPOSED" : ( n_once ? " AT_RISK" : " CLEAR" ) ) << "\n";
    }

  //! @brief Claim census: which states carry a retained continuity claim, and which do not.
  //!
  //! Under IC_STRONG every retained claim is imposed by eliminating its tau against a PIVOT
  //! ROW, and the plan then "drops pivot rows entirely, replacing them with the exact trace
  //! rows" (ocplan.hpp).  That consumption is real but invisible in the existing accounting:
  //! it is counted neither by dropped_physical_rows (explicit drops only) nor by
  //! exact_replacement_rows (the IC_TRACE route), both of which read 0 on a plan that
  //! consumes dozens of rows this way.
  //!
  //! An earlier version of this report flagged every pivot landing on a row that defines a
  //! claim-less state.  It over-counted badly -- 39 flagged where the true free-DOF count was
  //! 3 -- because a state appearing bare in several equations absorbs the loss.
  //! print_column_coverage_report() is the exact test.  What survives here is the census:
  //! "which states have no continuity claim at all" is a one-line fact worth reading next to
  //! the coverage numbers, because a state with no claim AND a single defining row is the
  //! shape that has produced every exact-imposition kernel found so far.
  void print_tau_pivot_report
    ( std::ostream& os = std::cerr )
    const
    {
      std::set<size_t> claimed;
      for( auto const& kv : _plan.claims_frozen() )
        claimed.insert( std::get<1>( kv.first ) );

      if( options.INTERFACE.IMPOSITION == Options::IC_STRONG ){
        auto const& elim = _plan.strong_projection();
        size_t const n_exp = _plan.n_trace_var_explicit();
        os << "  [claims] " << _plan.n_trace_constraint() << " claim(s): "
           << elim.pivot_rows.size() << " Schur-eliminated (rank=" << elim.rank
           << " n_trace=" << elim.n_trace << "), " << n_exp << " kept EXPLICIT\n";

        // A claim kept explicit is one for which NO pivot row could be found -- which
        // happens when its coupling is zero, so the tau injects into nothing.  Such a tau
        // is a column of the augmented system that no equation touches: an exact kernel
        // direction, and a constraint row that cannot be imposed, putting a FLOOR under the
        // least-squares residual.  Measured on the ALE mesh model: 99 claims = 94 eliminated
        // + 5 explicit, and the spectrum audit reported exactly those 5 columns as
        // `numerically EMPTY`, each a unit right-null vector with ||J*v|| = 0.
        if( n_exp ){
          // CORRECTED.  An earlier version said these claims "had no eliminable pivot".
          // They did not: both keep-explicit routes are DELIBERATE (see the
          // _keepExplicitReason doc).  The defect is that the route promises an explicit
          // tau column AND a real continuity row, and the row does not materialise --
          // leaving the column unsupported.  Reasons are printed per claim.
          os << "  [claims] *** " << n_exp << " claim(s) kept EXPLICIT (reason: 1=protected"
                "+rank-extending, 2=all-LINK replacements, 3=protected): ***\n";
          size_t shown = 0;
          for( auto const& ck : _plan.explicit_tau_claims() ){
            if( shown++ >= 16 ){ os << "      ... " << ( n_exp - 16 ) << " more\n"; break; }
            size_t const sid = std::get<1>( ck ), did = std::get<3>( ck );
            std::string dn = "dom?";
            for( auto const& kv : _mDom )
              if( static_cast<size_t>( kv.first.id().second ) == did ){ dn = kv.first.name(); break; }
            auto itr = _keepExplicitReason.find( ck );
            os << "      state=" << state_name_of_id( sid )
               << " dir=" << dn
               << " iel_lo=" << std::get<4>( ck )
               << " face="   << std::get<5>( ck )
               << " kind="   << std::get<6>( ck )
               << " blk="    << std::get<0>( ck )
               << "  reason=" << ( itr == _keepExplicitReason.end() ? 0 : itr->second )
               << "\n";
          }
        }
      }
      else
        os << "  [claims] imposition is not IC_STRONG; no Schur pivots\n";

      os << "  [claims] states WITH a retained claim:";
      for( size_t sid : claimed ) os << " " << state_name_of_id( sid );
      os << "\n  [claims] states WITHOUT any claim:";
      size_t n_unclaimed = 0;
      for( auto const& kv : _mVar ){
        size_t const sid = static_cast<size_t>( kv.first.id().second );
        if( !claimed.count( sid ) ){ os << " " << kv.first.name(); ++n_unclaimed; }
      }
      if( !n_unclaimed ) os << " (none)";
      os << "\n";
    }


  //! @brief Canonical trace-slot form of a claim key, identical to the plan builder's.

//==========================================================================================================
// INTERFACE PLAN -- claims -> realisations -> rows and columns.  The DECISIONS live in ocplan.hpp; what is here
// is the assembly: which claims exist, which rows receive them, and how the chosen realisation becomes
// matrix structure.
// (banner placed from the definitions themselves: this region is lines 4167-5948 of the pre-rev268 file)
//==========================================================================================================


  // =====================================================================
  // rev107 [taucol] -- DIRECT NAMED-COLUMN INSTRUMENT.  Instrument only.
  // =====================================================================
  //
  // WHY THIS EXISTS.  handoff rev32 section 4.2 records three certificates for the
  // orphan-tau kernel, all three PROPOSED FROM READING THE HEADER and all three
  // refuted by instruments that already existed and already disagreed:
  //
  //   structural coverage (_jacPatternEqn)   claimed "column in no row"
  //                                          measured "singly-covered = 3"  -- present
  //   blk.bcol coupling                      claimed "zero coefficient in B"
  //                                          measured "with_row=3, ORPHAN=0" -- backwards
  //   trace_constraints == trace_variables   claimed "short by k"
  //                                          never evaluated (twice gated)
  //
  // Ground truth never moved: [tier0-audit] numerically EMPTY columns = 3/7/15, k = 3/7/15.
  //
  // THIS IS NOT A FOURTH CERTIFICATE.  It asserts nothing about WHY the column is
  // empty.  It follows ONE NAMED COLUMN from the plan to the assembled system and
  // reports what it finds, distinguishing the three outcomes handoff rev32 section 5
  // named as the ones that need separating:
  //
  //   (1) row absent       -- the claim never produced a receiver term at all;
  //   (2) row present,     -- a term exists but its coefficient is zero, or the
  //       constrains nothing   eval-time dispatch never fires for it;
  //   (3) row dropped      -- a term exists, fires, and lands on a row LOOP 2
  //       downstream           removes; the [stage2-preflight] emits/drops accounting
  //                            is then where the coupling went.
  //
  // The three are NOT distinguishable from any existing report, because every
  // existing report is computed on ONE SIDE of the plan/assembly boundary.  B is
  // the PROMISE -- what a claim was minted WITH.  tier0-audit is the DELIVERY --
  // what reaches the assembled Jacobian.  This instrument spans the boundary and
  // prints both halves side by side, so the two can be compared without inference.
  //
  // TWO HALVES, DELIBERATELY SEPARATE, DELIBERATELY LABELLED:
  //   [taucol/plan]      emitted from setup(), off the COMMITTED plan.
  //   [taucol/assembly]  emitted from the first _eval_eqn(), off the rows actually
  //                      built, after the Schur substitution.
  // Reading one and attributing it to the other is the failure mode that printed
  // `infer=none` for five revisions (rev99 set _evoInferProvenance inside
  // _infer_evolution_domain and read it from the hoist, which is a REPLAY of those
  // rules, not a call).  Hence: two tags, two sites, no shared derivation.
  //
  // Gated on CRONOS_AUDIT_TAUCOL.  Off => not one branch taken outside the report
  // itself, so the corpus is byte-identical.
  // =====================================================================


  //! The returned string must be CONTIGUOUS in .rodata -- the build oracle greps the
  //! linked binary for it.  Streaming the instrument name and its tag separately would put
  //! the two halves in different literals and the grep would silently never match.


  //! @brief Total number of collocation state coefficients (valid after setup())
  size_t n_colloc_sta
    ()
    const
    { return _nCollVar; }

  //! @brief Number of generated trace/tau unknowns appended after state coefficients.
  //!
  //! For IC_TRACE, this is the number of tau unknowns appended to var[].
  //! For IC_STRONG, this returns 0 because tau variables are internal
  //! Schur variables not exposed in var[] (set during commit).
  size_t n_colloc_trace
    ()
    const
    { return _nTraceVar; }

  //! @brief Total number of collocation input coefficients (valid after setup())
  size_t n_colloc_inp
    ()
    const
    { return _nCollInp; }

  //! @brief Number of columns (state followed by input coefficients) in sparse derivative (returned by deriv())
  size_t n_colloc_cols
    ()
    const
    { return _nCollVar + _nCollInp; }

  //! @brief Total number of collocation equation residuals, excluding scalar outputs (valid after setup())
  size_t n_colloc_eqn
    ()
    const
    { return _nCollEqn; }

  //! @brief Total number of scalar output functions
  size_t n_colloc_fct
    ()
    const
    { return _nCollFct; }

  //! @brief Total number of equation residual and output rows
  size_t n_colloc_rows
    ()
    const
    { return _nCollEqn + _nCollFct; }

  //! @brief Number of nonzeros in the sparse derivative of equation residuals
  size_t n_colloc_eqn_nnz
    ()
    const
    { return const_cast<OCFESLV*>(this)->_ensure_deriv_cache()? _n_jac_nnz( _jacPatternEqn ): 0; }

  //! @brief Number of nonzeros in the sparse derivative of output functions
  size_t n_colloc_fct_nnz
    ()
    const
    { return const_cast<OCFESLV*>(this)->_ensure_deriv_cache()? _n_jac_nnz( _jacPatternFct ): 0; }

  //! @brief Number of nonzeros in the full sparse derivative returned by deriv()
  size_t n_colloc_nnz
    ()
    const
    { return n_colloc_eqn_nnz() + n_colloc_fct_nnz(); }

  //! @brief Number of colors in the current compressed-forward derivative seed.
  size_t n_deriv_colors
    ()
    const
    { return const_cast<OCFESLV*>(this)->_ensure_deriv_cache()? _jacNColor: 0; }
  
  //! @brief First row index of user output function ndx in the stacked evaluation array.
  //! Output indices follow add_output() insertion order.
  size_t row_fct
    ( size_t const ndx )
    const
    {
      return _issetup && ndx < _mFct.size()? _nCollEqn + _mFct[ndx].row0
                                           : std::numeric_limits<size_t>::max();
    }

  //! @brief Row index of local row loc in user output function ndx.
  size_t row_fct
    ( size_t const ndx, size_t const loc )
    const
    {
      return _issetup && ndx < _mFct.size() && loc < _mFct[ndx].nrow
           ? _nCollEqn + _mFct[ndx].row0 + loc
           : std::numeric_limits<size_t>::max();
    }

  //! @brief Output-only row block [first,nrow] for user output function ndx.
  std::pair<size_t,size_t> blk_fct
    ( size_t const ndx )
    const
    {
      return _issetup && ndx < _mFct.size()? std::make_pair( _mFct[ndx].row0, _mFct[ndx].nrow )
                                           : std::make_pair( std::numeric_limits<size_t>::max(), size_t(0) );
    }
  


  //! @brief factorise (or store) the sensitivity Jacobian for a block.
  //!
  //! Called from BOTH _mono_fsens_setup and _march_fsens_setup so the two paths cannot drift.
  //! Assembles sparsely when SENSITIVITY.SPARSE, using the SAME triplet idiom the primal uses in
  //! solve(); otherwise the dense path.
  //!
  //! @return false on factorisation failure (dense/pinv paths only; the sparse path defers
  //!         its solve to apply time and so cannot fail here).
  bool _sens_factorise
    ( FSensBlock& B, size_t neq, size_t nv,
      std::vector<size_t> const& nnzeqn, char const* tag ) const
    {
      B.sparse = false; B.pinv = false;

      // ---- sparse (tier 0) -------------------------------------------------------------
      if( options.SENSITIVITY.SPARSE && !_knob_SENS_MINNORM_PINV() ){
        size_t nnz = 0;
        for( size_t r = 0; r < neq; ++r )
          for( size_t idx = 0; idx < nnzeqn[r] && idx < B.ceqn[r].size(); ++idx )
            if( B.ceqn[r][idx] < nv ) ++nnz;
        arma::umat loc( 2, nnz ? nnz : 1 );
        arma::vec  val( nnz ? nnz : 1 );
        size_t e = 0;
        for( size_t r = 0; r < neq; ++r )
          for( size_t idx = 0; idx < nnzeqn[r] && idx < B.ceqn[r].size(); ++idx ){
            size_t const c = B.ceqn[r][idx];
            if( c >= nv ) continue;
            loc(0,e) = (arma::uword)r; loc(1,e) = (arma::uword)c;
            val(e)   = B.geqn[r][idx]; ++e;
          }
        B.Jsp = arma::sp_mat( (arma::uword)neq, (arma::uword)nv );
        if( e > 0 )
          B.Jsp = arma::sp_mat( loc.cols(0,e-1), val.subvec(0,e-1),
                                (arma::uword)neq, (arma::uword)nv );
        B.sparse = true;
        // Free the dense factors: on a marched run these are held per window, so leaving
        // them allocated would defeat the point of the change.
        B.L.reset(); B.U.reset(); B.P.reset();
        if( _knob_AUDIT_SENS_NULL() )
          _disp( 3 ) << "  [sens-sparse] " << ( tag ? tag : "" ) << " " << neq << "x" << nv
                    << "  nnz=" << e
                    << "  density=" << std::scientific << std::setprecision(2)
                    << ( double(e) / ( double(neq)*double(nv) ) )
                    << "  (dense would be "
                    << ( double(neq)*double(nv)*8.0/1024.0/1024.0 ) << " MB)"
                    << std::defaultfloat << std::endl;
        return true;
      }

      // ---- dense: pinv (diagnostic) or LU (rev49 path) ---------------------------------
      arma::mat J( neq, nv, arma::fill::zeros );
      for( size_t r = 0; r < neq; ++r )
        for( size_t idx = 0; idx < nnzeqn[r] && idx < B.ceqn[r].size(); ++idx ){
          size_t const c = B.ceqn[r][idx];
          if( c < nv ) J( r, c ) = B.geqn[r][idx];
        }
      _sens_regularise_tau_block( J, B.cfct, B.gfct, tag );
      if( _knob_SENS_MINNORM_PINV() ){
        B.pinv = true;
        if( !arma::pinv( B.L, J ) ) return false;
        if( _knob_AUDIT_SENS_NULL() )
          _disp( 3 ) << "  [sens-pinv] " << ( tag ? tag : "" ) << " " << J.n_rows << "x"
                    << J.n_cols << " pseudo-inverse formed (DIAGNOSTIC; dense SVD)" << std::endl;
        return true;
      }
      return arma::lu( B.L, B.U, B.P, J );
    }

  //! @brief solve J s = b for a block, honouring whichever representation is stored.
  //! Sparse refactorises per call -- see the FSensBlock::Jsp note.
  static bool _sens_solve
    ( FSensBlock const& B, arma::vec const& b, arma::vec& s, bool transpose )
    {
      if( B.sparse ){
        try {
          arma::sp_mat const A = transpose ? B.Jsp.t() : B.Jsp;
          return arma::spsolve( s, A, b, "superlu" );
        }
        catch( ... ) { return false; }
      }
      if( B.pinv ){ s = transpose ? arma::vec( B.L.t() * b ) : arma::vec( B.L * b ); return true; }
      s = transpose
        ? arma::vec( B.P.t() * arma::solve( arma::trimatu( B.L.t() ),
                       arma::solve( arma::trimatl( B.U.t() ), b ) ) )
        : arma::vec( arma::solve( arma::trimatu( B.U ),
                       arma::solve( arma::trimatl( B.L ), B.P * b ) ) );
      return true;
    }

  //! @brief tau-block diagonal regularisation + null-space probe for the sensitivity
  //! Jacobian.  Called from BOTH _mono_fsens_setup and _march_fsens_setup so the two paths
  //! cannot drift; the marching path calls it once per window.
  //!
  //! Returns ||J||_inf (before regularisation), for reporting.
  double _sens_regularise_tau_block
    ( arma::mat& J,
      std::vector<std::vector<size_t>> const& cfct,
      std::vector<std::vector<double>> const& gfct,
      char const* tag ) const
    {
      size_t const nv = _nCollVar;
      size_t const ntau = _nTraceVar;
      size_t const tau0 = ( nv >= ntau ) ? ( nv - ntau ) : nv;   // tau columns are the tail

      double jinf = 0.;
      for( arma::uword r = 0; r < J.n_rows; ++r ){
        double rs = 0.;
        for( arma::uword c = 0; c < J.n_cols; ++c ) rs += std::fabs( J(r,c) );
        jinf = std::max( jinf, rs );
      }

      // ---- CLAIM (1): is dF/dvar zero on the tau columns? ------------------------------
      // If it is not, dF/dp is NOT unique -- s is determined only up to null(J), and a
      // component of dF/dvar along a null direction makes dF/dvar . s depend on which
      // solution the factorisation happens to return.  No factorisation fixes that; it would
      // be a formulation question about what F is allowed to depend on.
      //
      // Cheap by construction: dF/dvar is already assembled in gfct/cfct, so this costs one
      // pass over the output-function gradients and no linear algebra at all.
      if( _knob_AUDIT_SENS_NULL() && ntau ){
        double worst = 0.; size_t worst_f = 0, nhit = 0;
        for( size_t f = 0; f < gfct.size(); ++f ){
          double nrm = 0.;
          for( size_t idx = 0; idx < gfct[f].size() && idx < cfct[f].size(); ++idx )
            if( cfct[f][idx] >= tau0 && cfct[f][idx] < nv ){
              nrm += gfct[f][idx]*gfct[f][idx]; ++nhit;
            }
          nrm = std::sqrt( nrm );
          if( nrm > worst ){ worst = nrm; worst_f = f; }
        }
        _disp( 3 ) << "  [sens-null] " << ( tag ? tag : "" )
                  << " n_tau=" << ntau << " tau0=" << tau0
                  << " nonzero_dFdtau_entries=" << nhit
                  << " worst||dF/dvar|tau||=" << std::scientific << std::setprecision(6) << worst
                  << " (fct " << worst_f << ")"
                  << ( worst > 1e-10 ? "  ** NONZERO: dF/dp is NOT well defined; the tau block is"
                                       " not pure gauge for this objective **"
                                     : "  [gauge -- dF/dp well defined]" )
                  << std::defaultfloat << std::endl;
      }

      // ---- CLAIM (2): pin the gauge directions ----------------------------------------
      // delta on the tau diagonal ONLY.  s_tau -> 0, which IS the min-norm choice on gauge
      // directions, and dF/dp is unchanged to O(delta) BECAUSE dF/dvar is zero there.
      //
      // This is valid ONLY IF null(J) is contained in the tau span.  That does NOT follow
      // from claim (1): a null direction with a physical component would survive this and
      // J_reg would stay singular.  The delta-sweep is the test -- dF/dp must be FLAT in
      // delta and must agree with the IC_WEAK value for the same model.  If it drifts with
      // delta, this regularisation is invalid and the SPQR route is required.
      if( _knob_SENS_TAU_REG() > 0.0 && ntau ){
        double const delta = _knob_SENS_TAU_REG() * ( jinf > 0. ? jinf : 1.0 );
        for( size_t c = tau0; c < nv && c < J.n_cols; ++c )
          if( c < J.n_rows ) J( c, c ) += delta;
        if( _knob_AUDIT_SENS_NULL() )
          _disp( 3 ) << "  [sens-reg] " << ( tag ? tag : "" )
                    << " tau-block diagonal += " << std::scientific << std::setprecision(3)
                    << delta << " (" << _knob_SENS_TAU_REG() << " x ||J||_inf=" << jinf << ")"
                    << std::defaultfloat << std::endl;
      }
      return jinf;
    }


  //! @brief internal wrapper used by solve() when AUDIT_DUP_SPREAD is set.
  //!
  //! Adapts the raw pointer solve() carries to the vector the public reporter takes, stamps
  //! the imposition mode onto the label (the interpretation of the number depends entirely
  //! on it -- IC_WEAK is truncation limited and legitimately ~1e-11, the exact modes are
  //! round-off limited and anything above ~1e-17 is a finding), and emits the non-fatal
  //! DUP_SPREAD_WARN line under the exact modes only.
  void _report_dup_spread_gated
    ( double const* var, char const* label )
    const;

  //! @brief post-solve readout of the violation IC_TRACE absorbs into beta.
  //! No-op unless AUDIT_BETA is set and the projection was applied.  Read-only.
  //!
  //! `converged` is a REQUIRED argument, not a label.  The identity
  //! `eps*beta = -Z^T r_c^raw` holds AT A SOLUTION.  At a general iterate the continuity row
  //! is `R_i = (Cx)_i + eps*(P lambda)_i + r_c,i` and projecting gives
  //!
  //!     eps*beta = -Z^T r_c^raw + Z^T R
  //!
  //! so a non-converged reading carries the unconverged residual IN THE SAME SUBSPACE the
  //! measurement is about.  That contamination is not a small correction and is
  //! indistinguishable from the signal, so the interpretation is withheld rather than
  //! qualified.  MMPDE13 at `--small` does not converge under IC_TRACE (and sometimes not
  //! under IC_WEAK), which is exactly the run this guard exists for.
  void _report_beta_gated
    ( double const* var, char const* label, bool converged )
    const;


  //! @brief The reduction auxiliary representing d(@p parent)/d(@p dom), or nullptr when
  //! order reduction created no such auxiliary for that (state, direction) pair.
  //!
  //! DIAGNOSTIC / VERIFICATION USE.  RED_FULL already builds the mixed (LDG-shaped) system
  //! internally -- the auxiliary IS the first-order gradient unknown -- but it is otherwise
  //! invisible to a driver: eval_colloc() throws on an FFPartial node, and eval_solution()
  //! needs a state FFVar the driver never created.  Without this accessor the total interface
  //! flux
  //!        N = a*u - d*(du/dz)
  //! can only be formed by spectrally re-differentiating the primitive, which is a DIFFERENT
  //! object from the gradient the scheme actually solves for.  Verifying any flux-form
  //! interface mode, and measuring interface-jump diagnostics honestly, both require reading
  //! the auxiliary itself.
  //!
  //! The returned pointer aliases the environment's auxiliary table and is invalidated by the
  //! next setup().  Returns nullptr when REDUCE.ORDER made no auxiliary, when @p parent is not
  //! a differentiated state, or when @p dom is not one of the directions the auxiliary carries.
  FFVar const* auxiliary_state
    ( FFVar const& parent, FFVar const& dom )
    const
    {
      for( auto const& a : _auxDef ){
        if( a.parent.id().second < 0 )        continue;   // no parent recorded
        if( a.parent.id() != parent.id() )    continue;
        if( a.diff_dom.find( dom ) != a.diff_dom.end() ) return &a.aux;
      }
      return nullptr;
    }

  //! @brief Every reduction auxiliary in creation order, as (auxiliary, parent) pairs.
  //! Empty unless REDUCE.ORDER created some.  Companion to auxiliary_state() for drivers that
  //! need to discover the auxiliaries rather than ask for a specific one.  Diagnostic use.
  std::vector< std::pair<FFVar,FFVar> > auxiliary_states
    ()
    const
    {
      std::vector< std::pair<FFVar,FFVar> > v;
      v.reserve( _auxDef.size() );
      for( auto const& a : _auxDef ) v.emplace_back( a.aux, a.parent );
      return v;
    }


  //! @brief The post-index-reduction structural DOF audit result.
  t_DofAudit const& reduced_dof_audit() const { return _dofAudit; }


  //! @brief re-run the determinacy audit at a caller-supplied point (a
  //! physical initial guess, e.g. a monotone mesh) instead of the reference point;
  //! null -> reference point.  Prints the same [determinacy] lines, ref=user.
  bool audit_determinacy( double const* var=nullptr, double const* inp=nullptr );

  //! @brief residual and sparse Jacobian at a supplied var vector -- the same
  //! assembly solve() performs per iteration, with no Newton state attached.  Uses the
  //! environment's current marching state (override rows live inside a march window).
  //! Pattern is cached in A across calls.  Returns false if setup() has not run.
  bool assemble
    ( double const* var, double const* inp, double const* cst, t_JacAssembly& A,
      bool equilibrate=false );

  //! @brief same, for marching window w at setup time.  w == 0: the problem's
  //! IC rows; w > 0: the _marchICRows continuity-pin override rows (pinned against var
  //! itself when no terminal profile is stored).  Requires a march grid.
  bool assemble_window
    ( size_t w, double const* var, double const* inp, double const* cst, t_JacAssembly& A,
      bool equilibrate=false );


protected:

  // ======================================================================
  // Principal symbol and PDE classification - method declarations
  // ======================================================================

/*
  //! @brief Compute and cache the PDE classification at a given reference point.
  //! Must be called after reduce_order() and setup().
  //! @param state_ref  Reference values for state variables (length = _vVar.size())
  //! @param cst_ref    Reference values for constants (length = _vCst.size())
  //! @param dom_ref    Reference values for domain variables (length = _vDom.size())
  //! @param time_dom   Pointer to the time-like domain variable, or nullptr
  //! @param n_sample   Directions sampled for no-time classification
  //! @param imag_tol   Imaginary-part tolerance
  //! @return true on success
  bool classify_pde
    ( std::vector<double> const& state_ref,
      std::vector<double> const& cst_ref,
      std::vector<double> const& dom_ref,
      FFVar const*               time_dom  = nullptr,
      unsigned const             n_sample  = 32,
      double   const             imag_tol  = 1e-8 );
*/


  //! @brief Stage-4: post-index-reduction structural DOF audit (max-bipartite
  //! matching on _jacPatternEqn); fills _dofAudit.  See t_DofAudit.
  void _audit_reduced_dof();

  //! @brief assembly core shared by assemble(), assemble_window() and solve().
  //! r_known: the caller's residual at var (solve() passes its accepted line-search
  //! residual, which is eval(var), to avoid a second evaluation); null -> evaluated here.
  bool _assemble_core
    ( double const* var, double const* inp, double const* cst, t_JacAssembly& A,
      bool equilibrate, bool ic_override, double const* r_known );
  //! @brief query the fixed equation-row sparsity pattern into A (once per pattern).
  bool _assemble_pattern( t_JacAssembly& A );
  //! @brief residual at var with the marching override rows (the former apply_ic_res).
  bool _residual_core
    ( double const* var, double const* inp, double const* cst, double* r, t_JacAssembly& A, bool ic_override );
  //! @brief cache the IC_WEAK SAT re-injection stencil for the override rows.
  void _ensure_ic_sat( t_JacAssembly& A, bool ic_override );
  //! @brief override row -> pin column, from the row's own sparsity.
  long _ic_col_for_row( t_JacAssembly const& A, FFVar const& st, t_MarchICRow const& rr ) const;

  //! @brief deterministic nonzero reference point (init() profile, zeros -> 1).
  bool _reference_point
    ( std::vector<double>& x, std::vector<double>& u, std::string& tag );

  //! @brief setup-time numerical determinacy audit; fills _determinacy and
  //! prints one [determinacy] line per audited window.  Null pointers -> the
  //! deterministic reference point.
  void _audit_determinacy( double const* var_user=nullptr, double const* inp_user=nullptr );


};

//! @brief mc::OCFESLV is the former name of mc::OCFESLV, kept so that existing code compiles unchanged.
// rev270: the OCEnv spelling is GONE.  The alias that carried it was removed once the drivers were migrated
// by f0_tools/migrate_drivers.py (146 of 154 files; 4,917 OCEnv, 192 MC__OCENV_, 422 OCFE_OCENV_HEADER).
// The MACROS are why a shim was not enough: drivers #define MC__OCENV_* BEFORE including this header, so after
// the rev266 rename those defines reached nothing and silently disabled features -- PDE17's robustness probe
// (expected FAIL became PASS) and PDE18's guard suppression (setup then failed) -- while every gate driver,
// defining none of them, stayed byte-identical.

////////////////////////////////////////////////////////////////////////

inline std::ostream&
operator<<
( std::ostream& os, OCFESLV const& env )
{
  os << std::endl;

  if( !env._dag ){
    os << "NO DAG" << std::endl;
    return os;
  }
  
  if( env._mDom.empty() )
    os << "NO DOMAIN" << std::endl;
  else if( env._mDom.size() == 1 )
    os << "DOMAIN:" << std::endl;
  else
    os << "DOMAINS:" << std::endl;
  for( auto const& [var,dom] : env._mDom )
    os << "  " << var << " IN " << dom << std::endl; 
  os << std::endl;

  if( env._mVar.empty() )
    os << "NO VARIABLE";
  else if( env._mVar.size() == 1 )
    os << "VARIABLE:";
  else
    os << "VARIABLES:";
  for( auto const& [var,sdom] : env._mVar )
    os << " " << var;
  os << std::endl << std::endl;
 
  if( env._mInp.empty() )
    os << "NO INPUT";
  else if( env._mInp.size() == 1 )
    os << "INPUT:";
  else
    os << "INPUTS:";
  for( auto const& [var, sdom] : env._mInp )
    os << " " << var;
  os << std::endl << std::endl;

  if( env._vCst.empty() )
    os << "NO CONSTANT";
  else if( env._vCst.size() == 1 )
    os << "CONSTANT:";
  else
    os << "CONSTANTS:";
  for( auto const& var : env._vCst )
    os << " " << var;
  os << std::endl << std::endl;

  if( env._mEqn.empty() )
    os << "NO EXPRESSION" << std::endl;
  else if( env._mEqn.size() == 1 )
    os << "EXPRESSION:" << std::endl;
  else
    os << "EXPRESSIONS:" << std::endl;
  for( auto const& eqnrec : env._mEqn ){
    auto const& eqn = eqnrec.var;
    auto const& dom = eqnrec.dom;
    auto const& opt = *eqnrec.opt;
    try{
      auto sgeqn  = env._dag->subgraph( {eqn} );
      auto streqn = FFExpr::subgraph( env._dag, sgeqn );
      os << "  0 = " << streqn[0];
      if( !dom.empty() ) os << "  <== COLLOCATED ON ";
      bool first = true;
      for( auto const& [d, l] : dom ){
        if( !first ) os << " U ";
        first = false;
        switch( l ){
          case FFDom::ALL:                     os << "[" << d << "]"; break;
          case FFDom::LB:                      os << "|" << d;        break;
          case FFDom::UB:                      os <<        d << "|"; break;
          case FFDom::ALL-FFDom::LB:           os << "|" << d << "]"; break;
          case FFDom::ALL-FFDom::UB:           os << "[" << d << "|"; break;
          case FFDom::ALL-FFDom::LB-FFDom::UB: os << "|" << d << "|"; break;
          default:                             os << "**ERROR";       break;
        }
      }
      auto role_name = []( OCFESLV::EqnRole r ){
        switch( r ){
          case OCFESLV::EqnRole::INTERIOR:   return "INTERIOR";
          case OCFESLV::EqnRole::INITIAL:    return "INITIAL";
          case OCFESLV::EqnRole::BOUNDARY:   return "BOUNDARY";
          case OCFESLV::EqnRole::INTERFACE:  return "INTERFACE";
          case OCFESLV::EqnRole::LINK:       return "LINK";
          case OCFESLV::EqnRole::SURFACE:    return "SURFACE";
          case OCFESLV::EqnRole::DIAGNOSTIC: return "DIAGNOSTIC";
          default:                         return "AUTO";
        }
      };
      os << "  {role=" << role_name(opt.role)
         << ", block=" << opt.block_id
         << ", classify=" << (opt.participate_in_classification?"yes":"no")
         << ", sat=" << (OCFESLV::_sopt( opt ).receive_sat?"yes":"no")
         << ", donate=" << (OCFESLV::_sopt( opt ).donate_for_state_continuity?"yes":"no") << "}";
      os << std::endl;
    }
    catch(...){
      os << "  **ERROR" << std::endl;
    }
  }

  if( env._mFct.empty() )
    os << std::endl << "NO OUTPUT" << std::endl;
  else if( env._mFct.size() == 1 )
    os << std::endl << "OUTPUT:" << std::endl;
  else
    os << std::endl << "OUTPUTS:" << std::endl;
  for( auto const& fctrec : env._mFct ){
    auto const& fct = fctrec.var;
    try{
      auto sgfct  = env._dag->subgraph( {fct} );
      auto strfct = FFExpr::subgraph( env._dag, sgfct );
      os << "  " << strfct[0];
      if( fctrec.kind == OCFESLV::FctKind::DISTRIBUTED ){
        os << "  <== DISTRIBUTED ON ";
        bool first = true;
        for( auto const& [d, lim] : fctrec.grid ){
          if( !first ) os << ", ";
          first = false;
          os << d << ":" << lim;
        }
        if( !fctrec.point.empty() ){
          os << " AT ";
          first = true;
          for( auto const& [d, v] : fctrec.point ){
            if( !first ) os << ", ";
            first = false;
            os << d << "=" << v;
          }
        }
        os << "  [rows=" << fctrec.row0 << ":" << fctrec.row0+fctrec.nrow << "]";
      }
      else{
        auto const& dom = fctrec.point;
        if( !dom.empty() ) os << "  <== EVALUATED AT ";
        bool first = true;
        for( auto const& [d, v] : dom ){
          if( !first ) os << ", ";
          first = false;
          os << d << "=" << v;
        }
      }
      os << std::endl;
    }
    catch(...){
      os << "  **ERROR" << std::endl;
    }
  }

  if( env._issetup ){
    env._disp( 3 ) << std::endl;
    env._disp( 2 ) << "TOTAL COLLOCATED INPUTS:    " << env._nCollInp << std::endl;
    env._disp( 2 ) << "TOTAL COLLOCATED VARIABLES: " << env._nCollVar << std::endl;
    if( env._nTraceVar ){
      env._disp( 2 ) << "TOTAL TRACE/TAU VARIABLES: " << env._nTraceVar << std::endl;
      // IC_STRONG hybrid: the explicit subset is reported above as var[] trace
      // variables; the remaining accepted taus are eliminated Schur variables.
      if( env._plan.n_trace_var() > env._nTraceVar )
        env._disp( 2 ) << "TOTAL SCHUR TAU VARIABLES: "
                  << env._plan.n_trace_var() - env._nTraceVar
                  << " (internal, not in var[])" << std::endl;
    }
    else if( env._plan.n_trace_var() ){
      env._disp( 2 ) << "TOTAL SCHUR TAU VARIABLES: " << env._plan.n_trace_var()
                << " (internal, not in var[])" << std::endl;
    }
    env._disp( 2 ) << "TOTAL COLLOCATED EQUATIONS: " << env._nCollEqn << std::endl;
    env._disp( 2 ) << "TOTAL COLLOCATED OUTPUTS:   " << env._nCollFct << std::endl;
    //std::cout << "TOTAL COLLOCATED ROWS:      " << env._nCollEqn + env._nCollFct << std::endl;
  }

  return os;
}


inline void
OCFESLV::_reset
()
{
  _reset_model();
  _reset_colloc_state();   // the base's state: domains, classification, SAT coupling, auxiliaries, plans, counts, the plan





  _jacPatternEqn.clear();
  _jacPatternFct.clear();
  _jacLinearEqn.clear();
  _jacLinearFct.clear();
  _jacAllEqnLinear = false;
  _jacAllFctLinear = false;
  _linEqn.clear();
  _linFct.clear();
  _linCstVal.clear();
  _linEqnReady = false;
  _linFctReady = false;
  _jacColor.clear();
  _jacNColor = 0;
  
  _FEqn.clear();
  _FFct.clear();
  _FVar.clear();
  _FInp.clear();

  _jacReady      = false;
  _jacClassified = false;
}


inline bool
OCFESLV::init
( double* var, double* inp, double const* cst )
{
  if( !_issetup ) return false;
  try{
    if( cst && !_vCst.empty() ){
      _vCstVal.assign( cst, cst + _vCst.size() );

      // Keep the user-side storage in sync as well, so that a later setup()
      // or another copy from the user model sees the same constant values.
      if( _usr._vCstUsr.size() == _vCst.size() )
        _usr._vCstValUsr = _vCstVal;
    }
    else if( !_vCstVal.empty() && _vCstVal.size() != _vCst.size() ){
      throw Exceptions( Exceptions::CSTVAL );
    }

    _init_reference( var, inp, cst );
    if( inp ) for( auto const& [w,vals] : _usr._mInpFixUsr ) set_input_values( w, vals, inp );   // fix_input
    return true;
  }
  catch(...){
    return false;
  }
}


inline size_t
OCFESLV::select_strong_interface_donor_inode
( FFDom const& dom_in, size_t iel, size_t interface_inode, int lim,
  bool has_left_exact_value, bool has_right_exact_value )
const
{
  if( !dom_in.n_node ) throw Exceptions( Exceptions::DOMAIN );
  if( iel >= dom_in.n_elem ) throw Exceptions( Exceptions::DOMAIN );

  // Work with a node-initialised copy when the helper is called outside
  // setup(); setup-time calls normally use the already initialised _mDom entry.
  FFDom dom_work;
  FFDom const* pdom = &dom_in;
  if( dom_in.lgnodes().size() != dom_in.n_node ){
    dom_work = dom_in;
    if( !dom_work.set_nodes() ) throw Exceptions( Exceptions::DOMAIN );
    pdom = &dom_work;
  }
  FFDom const& dom = *pdom;

  size_t const n = dom.n_node;
  size_t donor = std::min( interface_inode, n-1 );
  bool const at_left_face  = ( donor == 0 );
  bool const at_right_face = ( donor + 1 == n );

  auto physically_allowed = [&]( size_t k ) -> bool
  {
    // The selected donor must be an equation row that survives the user domain
    // mask before interface replacement.  Physical endpoint exclusions remain
    // hard exclusions; interface replacement is handled separately by claims.
    switch( lim ){
      case FFDom::ALL:
        return true;
      case FFDom::ALL-FFDom::LB:
        return !( iel == 0 && k == 0 );
      case FFDom::ALL-FFDom::UB:
        return !( iel + 1 == dom.n_elem && k + 1 == n );
      case FFDom::ALL-FFDom::LB-FFDom::UB:
        return !( iel == 0 && k == 0 )
            && !( iel + 1 == dom.n_elem && k + 1 == n );
      case FFDom::LB:
      case FFDom::UB:
        return false;
      default:
        throw Exceptions( Exceptions::DOMAIN );
    }
  };

  if( !physically_allowed( donor ) ){
    // Fall back to the closest physically retained row.  This should only be
    // needed for unusual masks or future collocation variants.
    for( size_t r = 1; r < n; ++r ){
      if( donor >= r && physically_allowed( donor-r ) ){ donor -= r; break; }
      if( donor + r < n && physically_allowed( donor+r ) ){ donor += r; break; }
    }
  }

  bool const two_exact_faces = has_left_exact_value && has_right_exact_value;

  // Endpoint-first policy: first-order/direct normal derivative rows donate
  // the duplicated endpoint row. Interior shifting is permitted only if the
  // state-specific local rank check below rejects this endpoint choice.

  auto companion_donor = [&]() -> size_t
  {
    // Approximate the donor that the opposite exact face of the same element
    // will consume. This is used only for the tiny local normal-stencil rank
    // test.
    if( !two_exact_faces ) return n;
    if( at_left_face ){
      return physically_allowed( n-1 ) ? n-1 : n;
    }
    if( at_right_face ){
      if( ( dom.type == FFDom::LGL || dom.type == FFDom::CGL )
       && n > 2 && physically_allowed( 1 ) ) return 1;
      return physically_allowed( 0 ) ? 0 : n;
    }
    return n;
  };

  auto full_rank_with = [&]( size_t cand ) -> bool
  {
    if( cand >= n || !physically_allowed( cand ) ) return false;

    std::vector< std::vector<double> > rows;
    rows.reserve( n + 2 );

    std::vector<double> w;
    if( has_left_exact_value ){
      if( !dom.w_eval_lagrange( w, -1.0, 0 ) ) return false;
      rows.push_back( w );
    }
    if( has_right_exact_value ){
      if( !dom.w_eval_lagrange( w,  1.0, 0 ) ) return false;
      rows.push_back( w );
    }

    std::set<size_t> drop;
    drop.insert( cand );
    size_t const comp = companion_donor();
    if( comp < n ) drop.insert( comp );

    for( size_t k = 0; k < n; ++k ){
      if( !physically_allowed( k ) ) continue;
      if( drop.find( k ) != drop.end() ) continue;
      if( !dom.w_eval_lagrange( w, dom.lgnodes()[k], 1 ) ) return false;
      rows.push_back( w );
    }

    if( rows.empty() ) return false;
    arma::mat A( (arma::uword)rows.size(), (arma::uword)n, arma::fill::zeros );
    for( arma::uword i = 0; i < A.n_rows; ++i )
      for( arma::uword j = 0; j < A.n_cols; ++j )
        A(i,j) = rows[(size_t)i][(size_t)j];

    arma::vec sv;
    if( !arma::svd( sv, A ) ) return false;
    double const smax = sv.n_elem ? sv.max() : 0.;
    double const tol  = 256. * DBL_EPSILON * std::max( 1.0, smax )
                      * std::max( rows.size(), n );
    arma::uword rank = 0;
    for( arma::uword i = 0; i < sv.n_elem; ++i )
      if( sv(i) > tol ) ++rank;

    return rank >= std::min<arma::uword>( (arma::uword)n, A.n_rows );
  };

  if( full_rank_with( donor ) ) return donor;

  // Robust fallback: choose the nearest row that restores the local stencil
  // rank.  Prefer interior rows before the opposite endpoint, since keeping at
  // least one endpoint derivative row is the intended LGL/CGL safeguard.
  for( size_t r = 1; r < n; ++r ){
    if( donor >= r ){
      size_t const k = donor - r;
      if( k != 0 && k + 1 != n && full_rank_with( k ) ) return k;
    }
    if( donor + r < n ){
      size_t const k = donor + r;
      if( k != 0 && k + 1 != n && full_rank_with( k ) ) return k;
    }
  }
  for( size_t r = 1; r < n; ++r ){
    if( donor >= r && full_rank_with( donor-r ) ) return donor-r;
    if( donor + r < n && full_rank_with( donor+r ) ) return donor+r;
  }

  return donor;
}


inline std::vector<OCFESLV::t_EqnNode>
OCFESLV::equation_nodes
( FFVar const& eqn, bool retained_only )
const
{
  for( size_t i = 0; i < _mEqn.size(); ++i )
    if( _mEqn[i].var.id() == eqn.id() )
      return equation_nodes( i, retained_only );
  throw Exceptions( Exceptions::INDEX );
}

inline std::vector<OCFESLV::t_EqnNode>
OCFESLV::equation_nodes
( size_t eqn_index, bool retained_only )
const
{
  if( eqn_index >= _mEqn.size() ) throw Exceptions( Exceptions::INDEX );
  auto const& rec = _mEqn[eqn_index];
  bool const weak = ( options.INTERFACE.IMPOSITION == Options::IC_WEAK );

  struct DimNode {
    FFVar var;
    FFDom const* dom = nullptr;
    int lim = FFDom::ALL;
    size_t iel = 0;
    size_t ncur = 0;
  };

  std::vector<DimNode> dims;
  dims.reserve( rec.dom.size() );
  for( auto const& [dvar, lim] : rec.dom ){
    auto itd = _mDom.find( dvar );
    if( itd == _mDom.end() ) throw Exceptions( Exceptions::DOMAIN );
    DimNode dn;
    dn.var = dvar;
    dn.dom = &itd->second;
    dn.lim = lim;
    dims.push_back( dn );
  }

  std::vector<t_EqnNode> out;
  if( dims.empty() ){
    t_EqnNode nd;
    nd.eqn_index = eqn_index;
    nd.flat_index = 0;
    nd.base_keep = true;
    nd.retained = true;
    if( !retained_only || nd.retained ) out.push_back( std::move(nd) );
    return out;
  }

  auto physical_inode = []( FFDom const& dom, int lim, size_t local_inode ) -> size_t
  {
    if( lim == FFDom::LB ) return 0;
    if( lim == FFDom::UB ) return dom.n_node - 1;
    return local_inode;
  };

  std::function<void(size_t,size_t&)> visit_elem = [&]( size_t d, size_t& flat )
  {
    if( d == dims.size() ){
      for( auto& dn : dims )
        dn.ncur = ( dn.lim == FFDom::LB || dn.lim == FFDom::UB ) ? 1 : dn.dom->n_node;

      size_t ntot = 1;
      for( auto const& dn : dims ) ntot *= dn.ncur;

      std::vector<size_t> midx( dims.size(), 0 );
      for( size_t loc = 0; loc < ntot; ++loc, ++flat ){
        size_t rem = loc;
        for( size_t j = 0; j < dims.size(); ++j ){
          midx[j] = rem % dims[j].ncur;
          rem /= dims[j].ncur;
        }

        t_EqnNode nd;
        nd.eqn_index = eqn_index;
        nd.flat_index = flat;
        nd.base_keep = true;

        for( size_t j = 0; j < dims.size(); ++j ){
          auto const& dn = dims[j];
          size_t iel = dn.iel;
          size_t inode = physical_inode( *dn.dom, dn.lim, midx[j] );
          if( dn.lim == FFDom::LB ) iel = 0;
          if( dn.lim == FFDom::UB ) iel = dn.dom->n_elem - 1;

          switch( dn.lim ){
            case FFDom::LB:
            case FFDom::UB:
            case FFDom::ALL:
              break;
            case FFDom::ALL-FFDom::LB:
              if( iel == 0 && inode == 0 ) nd.base_keep = false;
              break;
            case FFDom::ALL-FFDom::UB:
              if( iel + 1 == dn.dom->n_elem && inode + 1 == dn.dom->n_node ) nd.base_keep = false;
              break;
            case FFDom::ALL-FFDom::LB-FFDom::UB:
              if( iel == 0 && inode == 0 ) nd.base_keep = false;
              if( iel + 1 == dn.dom->n_elem && inode + 1 == dn.dom->n_node ) nd.base_keep = false;
              break;
            default:
              throw Exceptions( Exceptions::DOMAIN );
          }

          nd.elem[dn.var] = iel;
          nd.inode[dn.var] = inode;
          double const xi = dn.dom->lgnodes()[inode];
          nd.coord[dn.var] = dn.dom->elem_lo( iel )
                           + 0.5 * dn.dom->elem_width( iel ) * ( xi + 1.0 );
        }

        out.push_back( std::move(nd) );
      }
      return;
    }

    if( dims[d].lim == FFDom::LB || dims[d].lim == FFDom::UB ){
      dims[d].iel = ( dims[d].lim == FFDom::UB ? dims[d].dom->n_elem - 1 : 0 );
      visit_elem( d+1, flat );
    }
    else{
      for( size_t iel = 0; iel < dims[d].dom->n_elem; ++iel ){
        dims[d].iel = iel;
        visit_elem( d+1, flat );
      }
    }
  };

  size_t flat = 0;
  visit_elem( 0, flat );

  auto row_code = [&]( t_EqnNode const& nd, size_t override_dim,
                       size_t override_elem, size_t override_inode ) -> size_t
  {
    size_t code = 0;
    size_t mult = 1;
    for( size_t j = 0; j < dims.size(); ++j ){
      auto const& dn = dims[j];
      size_t ed = nd.elem.at( dn.var );
      size_t kd = nd.inode.at( dn.var );
      if( j == override_dim ){
        ed = override_elem;
        kd = override_inode;
      }
      size_t const radix = dn.dom->n_elem * dn.dom->n_node;
      code += ( ed * dn.dom->n_node + kd ) * mult;
      mult *= radix;
    }
    return code;
  };

  std::map<size_t,size_t> row_to_node;
  for( size_t i = 0; i < out.size(); ++i )
    row_to_node[row_code( out[i], dims.size(), 0, 0 )] = i;

  std::set<size_t> donor_indices;
  for( auto const& nd : out ){
    if( !nd.base_keep ) continue;
    for( size_t j = 0; j < dims.size(); ++j ){
      auto const& dn = dims[j];
      if( dn.lim == FFDom::LB || dn.lim == FFDom::UB ) continue;
      size_t const inode = nd.inode.at( dn.var );
      size_t const iel   = nd.elem.at( dn.var );
      if( _keep_eqn_node( rec.var, *rec.opt, dn.var, *dn.dom,
                          iel, inode, dn.lim, weak, eqn_index ) )
        continue;

      // Public diagnostics do not have a unique donor state for this row.
      // Treat only neighbouring element interfaces as exact value constraints;
      // the actual state-specific physical-face logic is frozen in the setup
      // donor/acceptor plan consumed by _eval_eqn().
      bool const has_left_exact_value  = ( iel > 0 );
      bool const has_right_exact_value = ( iel + 1 < dn.dom->n_elem );
      size_t const donor_inode = select_strong_interface_donor_inode
        ( *dn.dom, iel, inode, dn.lim, has_left_exact_value, has_right_exact_value );
      size_t const dk = row_code( nd, j, iel, donor_inode );
      auto it = row_to_node.find( dk );
      if( it != row_to_node.end() ) donor_indices.insert( it->second );
    }
  }

  for( size_t i = 0; i < out.size(); ++i ){
    out[i].interface_donor = donor_indices.find( i ) != donor_indices.end();
    out[i].retained = out[i].base_keep && !out[i].interface_donor;
  }

  if( retained_only ){
    std::vector<t_EqnNode> kept;
    kept.reserve( out.size() );
    for( auto const& nd : out ) if( nd.retained ) kept.push_back( nd );
    return kept;
  }
  return out;
}


inline OCFESLV::t_Coll
OCFESLV::_copy_inp_disc
( t_Coll const& in, std::vector<FFVar> const& all,
  std::map<FFVar::pt_idVar,size_t> const& ndx )
{
  t_Coll out;
  for( auto const& [inp, disc] : in ){
    t_Dom disc_local;
    //for( auto const& [dom, opt] : disc )
    //  disc_local.insert( { _local( dom, all, ndx ), opt } );
    out.insert( { _local( inp, all, ndx ), std::move( disc_local ) } );
    for( auto const& [dom, opt] : disc )
      // Preserve input-only n_node=1 discretisations without invoking
      // FFDom's validating copy constructor.
      disc_local[ _local( dom, all, ndx ) ] = opt;
    out[ _local( inp, all, ndx ) ] = std::move( disc_local );
  }
  return out;
}


inline bool
OCFESLV::deep_copy_from
( OCFESLV const& src )
{
  if( this == &src ) return true;
  // rev309b: the SOLVER's options are part of the environment being copied.  The frozen interface plan below
  // was built from src.options, and the guard in _eval_eqn compares the frozen record against options -- so a
  // copy that kept its own default-constructed options would fail that guard on its first evaluation.
  options = src.options;
  // A source that is not set up has no working model to copy.  Reported, not thrown: the copy constructor and
  // assignment cannot report a throw, and an unguarded copy of such an environment would terminate the process.
  if( !src._issetup || !src._dag ){
    _reset();
    _release_working_dag();
    return false;
  }

  _reset();

  // The primary cached symbol (_symbol) is usually aliased by one entry of
  // _blockSymbol, but this is not guaranteed after failed/partial
  // classification or after callers have explicitly queried/pruned symbol
  // data.  Register both roots so every FFVar that may be copied below has a
  // corresponding clone in the destination DAG.
  std::map<int,t_Symbol> symRoots( src._blockSymbol );
  if( src._classified )
    symRoots.emplace( std::numeric_limits<int>::min(), src._symbol );

  // Solver roots, in the order the former single registration used: symbol coefficients per block, then
  // auxiliary definitions.  FFModel imports them with the model roots in ONE insert.
  std::vector<FFVar> extraRoots;
  for( auto const& [blk, sym] : symRoots )
    for( auto const& row : sym.vCoeff )
      for( auto const& var : row ) extraRoots.push_back( var );
  for( auto const& aux : src._auxDef ){
    extraRoots.push_back( aux.aux );
    extraRoots.push_back( aux.expr );
    for( auto const& d : aux.dom )      extraRoots.push_back( d );
    for( auto const& d : aux.diff_dom ) extraRoots.push_back( d );
    extraRoots.push_back( aux.parent );
  }

  std::map<FFVar::pt_idVar,size_t> ndxAllVars;
  std::vector<FFVar> vAllVarsLoc;
  if( !_deep_copy_model( src, extraRoots, vAllVarsLoc, ndxAllVars ) ){
    _reset();
    _release_working_dag();
    return false;
  }

  // Declared solver data (user-DAG keys, like _usr: copied as-is, no remapping).
  _adaptDomUsr = src._adaptDomUsr;

  // Deep-copy local DAG containers.
  _mDom = _copy( src._mDom, vAllVarsLoc, ndxAllVars );
  for( auto& [var, dom] : _mDom )
    if( !dom.set_nodes() ) return false;

  _mInpDisc = _copy_inp_disc( src._mInpDisc, vAllVarsLoc, ndxAllVars );
  _mInpCont = _copy_inp_cont( src._mInpCont, vAllVarsLoc, ndxAllVars );
  for( auto& [inp, disc] : _mInpDisc )
    for( auto& [dom, opt] : disc )
      if( !opt.set_nodes() ) return false;


  
  _evolution_dom_user = src._evolution_dom_user;
  _evolution_dom_set  = src._evolution_dom_set;
  if( src._evolution_dom_set )
    _evolution_dom_var = _local( src._evolution_dom_var, vAllVarsLoc, ndxAllVars );


  // Solver configuration: homotopy parameters (remapped to this DAG) and per-stage step caps.
  // Assigned, not appended: a destination that already held settings is overwritten.
  _vHom.clear();
  _vHom.reserve( src._vHom.size() );
  for( auto const& h : src._vHom ){
    t_Homotopy hc = h;
    hc.param = _local( h.param, vAllVarsLoc, ndxAllVars );
    _vHom.push_back( std::move( hc ) );
  }
  _mHomCap = src._mHomCap;

  _auxDef.clear();
  _auxVarID.clear();
  _auxDef.reserve( src._auxDef.size() );
  for( auto const& a : src._auxDef ){
    t_AuxDef b;
    b.aux         = _local( a.aux,  vAllVarsLoc, ndxAllVars );
    b.expr        = _local( a.expr, vAllVarsLoc, ndxAllVars );
    b.dom         = _copy( a.dom, vAllVarsLoc, ndxAllVars );
    b.diff_dom    = _copy( a.diff_dom, vAllVarsLoc, ndxAllVars );
    b.parent      = _local( a.parent, vAllVarsLoc, ndxAllVars );
    b.block_id    = a.block_id;
    _auxVarID.insert( b.aux.id() );
    _auxDef.push_back( std::move( b ) );
  }

  // Copy setup metadata and id-based caches.  These should not be
  // recomputed in a different traversal order during deep copy.
  _nCollEqn = src._nCollEqn;
  _evoDerivRows = src._evoDerivRows;   // the rows the row-level causal design excludes at element starts (2 Oct)
  _nCollFct = src._nCollFct;
  _nCollVar = src._nCollVar;
  _nTraceVar = src._nTraceVar;
  _traceVarOffset = src._traceVarOffset;
  _nTraceConstraint = src._nTraceConstraint;
  _nCollInp = src._nCollInp;

  // Copy the interface plan using OCPlan::copy_from, which handles the
  // equation-id remapping needed after DAG::insert clones expression nodes.
  {
    std::map<size_t,size_t> eqn_id_map;
    for( size_t i = 0; i < src._mEqn.size() && i < _mEqn.size(); ++i )
      eqn_id_map[ static_cast<size_t>( src._mEqn[i].var.id().second ) ]
        = static_cast<size_t>( _mEqn[i].var.id().second );
    _plan.copy_from( src._plan, eqn_id_map );
  }
  _interfacePlanBuildAllowed = false;

  // Copy the setup-time high-index / DOF-analysis records so a deep copy reports
  // the same reduction_plan() / reduced_dof_audit() as its source WITHOUT a
  // re-setup.  The reduced *executable* state (differentiated equations, pinned
  // states, minted dummy-derivative auxiliaries) already rides in the _mEqn /
  // _mVar / _auxDef containers copied above; these two members are the
  // diagnostic records on top of that.
  //
  //  * _dofAudit is plain data (bool/size_t) -> a direct copy is safe.
  //  * _reductionPlan holds FFVars (the ORIGINAL differentiated constraint
  //    residual + the pinned witness).  _reduce_high_index can rewrite the
  //    constraint out of _mEqn during setup, so it is not guaranteed reachable
  //    from the containers copied above.  Rather than special-case the
  //    registration, copy these vars straight into the local DAG with
  //    FFGraph::insert: it clones any nodes not already present (reusing the
  //    rest), and copied variables keep the same .id() in the new DAG.
  //  * _idxCache (structural-index memoization) is deliberately NOT copied: it
  //    is a lazy performance cache rebuilt on demand, and a copy that does not
  //    re-run setup never re-queries it, so the reset()-cleared empty state is
  //    correct.
  _dofAudit = src._dofAudit;
  _determinacy = src._determinacy;

  _reductionPlan.clear();
  _reductionPlan.max_index = src._reductionPlan.max_index;
  _reductionPlan.resolved  = src._reductionPlan.resolved;
  _reductionPlan.assigns.reserve( src._reductionPlan.assigns.size() );
  {
    std::vector<FFVar> planSrc, planLoc;
    planSrc.reserve( 2 * src._reductionPlan.assigns.size() );
    for( auto const& a : src._reductionPlan.assigns ){
      planSrc.push_back( a.constraint );
      planSrc.push_back( a.pinned_var );
    }
    if( !planSrc.empty() )
      _dag->insert( src._dag, planSrc, planLoc );
    for( size_t i = 0; i < src._reductionPlan.assigns.size(); ++i ){
      t_ReductionPlan::t_Assign b;
      b.block_id   = src._reductionPlan.assigns[i].block_id;
      b.n_diff     = src._reductionPlan.assigns[i].n_diff;
      b.constraint = planLoc[2*i];
      b.pinned_var = planLoc[2*i+1];
      _reductionPlan.assigns.push_back( std::move( b ) );
    }
  }

  _mDep   = src._mDep;
  _depEqn = src._depEqn;
  _depFct = src._depFct;

  _vVar = _copy( src._vVar, vAllVarsLoc, ndxAllVars );
  _vInp = _copy( src._vInp, vAllVarsLoc, ndxAllVars );
  _vDom = _copy( src._vDom, vAllVarsLoc, ndxAllVars );

  _depVar = src._depVar;
  _depInp = src._depInp;
  _depCst = src._depCst;
  _depDom = src._depDom;

  _blkVar.clear();
  for( auto const& [var, dat] : src._blkVar )
    _blkVar.insert( {_local( var, vAllVarsLoc, ndxAllVars ), dat} );

  _blkInp.clear();
  for( auto const& [var, dat] : src._blkInp )
    _blkInp.insert( {_local( var, vAllVarsLoc, ndxAllVars ), dat} );

  // Subgraphs contain pointers into the source DAG.  Rebuild only those
  // pointer lists on the cloned DAG; this is not a setup pass.
  _sgEqn.resize( _mEqn.size() );
  size_t ndx = 0;
  for( auto const& eqn : _mEqn )
    _sgEqn[ndx++] = _dag->subgraph( 1, &eqn.var );

  _sgFct.resize( _mFct.size() );
  ndx = 0;
  for( auto const& fct : _mFct )
    _sgFct[ndx++] = _dag->subgraph( 1, &fct.var );

  // Deep-copy classification and rebuild face-data pointers into the
  // destination symbol vectors.  When classification is not active, do not
  // propagate stale symbol objects left by an unsuccessful/aborted classifier
  // attempt.
  _classified          = src._classified;
  _classificationSerial= src._classificationSerial;
  _classification      = src._classified ? src._classification : t_Classify();
  _blockClassification = src._classified ? src._blockClassification
                                         : std::map<int,t_Classify>();
  _symbol              = t_Symbol();
  _face_data.clear();
  _blockSymbol.clear();
  _blockFaceData.clear();

  if( _classified ){
    _symbol    = _copy( src._symbol, vAllVarsLoc, ndxAllVars );
    _face_data = _build_face_data( _symbol, _classification );

    for( auto const& [bid,sym] : src._blockSymbol )
      _blockSymbol[bid] = _copy( sym, vAllVarsLoc, ndxAllVars );

    for( auto const& [bid,cls] : _blockClassification ){
      auto itsym = _blockSymbol.find( bid );
      if( itsym != _blockSymbol.end() )
        _blockFaceData[bid] = _build_face_data( itsym->second, cls );
    }
  }

  // Copy the setup-time SAT coupling table when available.  It is indexed
  // by equation/state positions and stores domain ids; FFVar ids are
  // preserved by the cloned local DAG.
  _satCoupling           = src._satCoupling;
  _satCouplingBuilt      = src._satCouplingBuilt;
  _satCouplingClassified = src._satCouplingClassified;
  _satCouplingSerial     = src._satCouplingSerial;

  // Copy sparse derivative/coloring cache.
  _jacPatternEqn     = src._jacPatternEqn;
  _jacPatternFct     = src._jacPatternFct;
  _jacLinearEqn      = src._jacLinearEqn;
  _jacLinearFct      = src._jacLinearFct;
  _jacAllEqnLinear   = src._jacAllEqnLinear;
  _jacAllFctLinear   = src._jacAllFctLinear;
  _jacColor          = src._jacColor;
  _jacNColor         = src._jacNColor;
  _jacReady          = src._jacReady;
  _jacInterfaceType  = src._jacInterfaceType;
  _jacImpositionType = src._jacImpositionType;
  _jacClassified     = src._jacClassified;
  _linEqnReady = _linFctReady = false;
  _linEqn.clear();
  _linFct.clear();
  _linCstVal.clear();

  _setup_node_cache();
  _setup_eval_plans();

  // Size and initialise pointer-bearing collocation and AD workspaces.
  _reset_colloc( _OCVar, _OCInp, _OCDom, _OCCst );
  _reset_colloc( _OCDepVar, _OCDepInp, _OCDepDom, _OCDepCst );
  _reset_colloc( _OCFVar, _OCFInp, _OCFDom, _OCFCst );

  _OCEqn    = OCVar<double>();
  _OCDepEqn = OCVar<FFDep>();
  _OCFEqn   = OCVar<t_FAD>();

  _wrkOC.clear();
  _wrkOCDep.clear();
  _wrkOCF.clear();

  _reset_autodiff( _FVar, _FInp, _FEqn, _FFct );

  // rev107 [taucol/plan]: the copied environment gets its OWN table.  The observation
  // state is deliberately not propagated by the copy (a copy is a different
  // process-visible object), so without this call _tauColPlanReady stays false and
  // OCFE_deepcopy -- an AT-RISK driver with its own history here -- would report
  // nothing at all under the flag.  Must precede the audit below, which evals.
  if( _taucol_enabled() ) _report_tau_column_plan( _disp( 3 ) );

  // Structural self-check on the copied state: cheap, and catches any
  // theoretical misalignment between the source's setup and what was
  // actually copied across the DAG boundary.
  if( !_audit_strong_interface_rows() ){
    _clear_frozen_interface_consumers();
    _issetup = false;
    return false;
  }

  // ------------------------------------------------------------------
  // Reduced-space control registry (register_control) and SOLVE.MARCHING
  // collapse state.  deep_copy_from does NOT re-run setup(), so these members --
  // populated at/after setup() -- must be carried explicitly.  FFGraph::insert
  // preserves FFVar::id, and every input/state lookup here (set_input_values,
  // _mInp / _blkInp, the marching helpers) resolves by id through lt_FFVar, so a
  // verbatim copy is correct -- no _local() remapping is needed.  Without these,
  // a copy reports no controls (n_control_dof()==0, decode/encode no-op) and
  // cannot march (_marchGrid empty => _solve_marching early-returns =>
  // val_functions() empty), which is what breaks a COPY-policy FFOCFESLV.
  // _fsBlocks is transient (rebuilt by _march_fsens_setup) and left cleared.
  // ------------------------------------------------------------------
  _controls       = src._controls;
  _nControlDof    = src._nControlDof;
  _decisionFlag       = src._decisionFlag;
  _primalValid    = false; _primalInp.clear(); _primalVar.clear();   // copy starts with a cold primal cache
  _reducedDepValid = false;   // and recomputes the dependency pattern on first access
  _rowProvenanceReady = false; _rowEqnNdx.clear();   // and rebuilds row provenance on first eval
  _eqnDiffInReady = false; _eqnDiffIn.clear();       // and the row x direction Partial table
  _marchICRowsReady = false; _marchICRows.clear(); _marchICTerminal.clear(); _marchICColToNode.clear(); _marchICActive = false;
  _marchGrid      = src._marchGrid;
  _marchTransfers = src._marchTransfers;
  _deferredValue   = src._deferredValue;

  _issetup = true;
  return true;
}


// ======================================================================
// OCFESLV::setup
// ======================================================================
inline bool
OCFESLV::_validate_model
()
{
  // validates through FFModel::_check_model, then sets the declared outputs' rows (as _on_model_ready does the working ones)
  return _check_consistency( _usr._dagUsr, _usr._vCstUsr, &_usr._vCstValUsr,
                             _usr._mDomUsr, _usr._mVarUsr, _usr._mInpUsr,
                             _usr._mEqnUsr, _usr._mFctUsr );
}

inline bool
OCFESLV::_on_setup_begin
()
{
  _setupStatus = SetupStatus::OK;
  _jacReady = false;
  _jacPatternEqn.clear();
  _jacPatternFct.clear();
  _jacColor.clear();
  _jacNColor = 0;
  _jacClassified = false;
  _plan.reset();
  _blockSymbol.clear();
  _blockClassification.clear();
  _blockFaceData.clear();
  _idxCache.clear();
  _face_data.clear();
  // 2026-10-09: an option that EXPLICITLY requests a backend this build lacks is refused here, before any work.
  // (Default / AUTO choices fall back: SOLVE_SUPERLU to dense LAPACK, DET_AUTO to Eigen or dense.)
  [[maybe_unused]] auto refuse = [&]( char const* what, char const* need ){   // unused when every backend is built
    std::cerr << "OCFESLV::setup ** " << what << " requires " << need << ", which this build does not include"
              << " (build_info()['backends']); choose another value or use a build with it" << std::endl;
    _setupStatus = SetupStatus::BACKEND_UNAVAILABLE;
    return false;
  };
#if !defined(CRONOS__WITH_SPQR)
  if( options.SOLVE.FACTORIZATION == Options::SOLVE_SPQR )
    return refuse( "SOLVE.FACTORIZATION = SOLVE_SPQR", "SuiteSparseQR (the GPL build)" );
  if( options.DETERMINACY.BACKEND == Options::DET_SPQR )
    return refuse( "DETERMINACY.BACKEND = DET_SPQR", "SuiteSparseQR (the GPL build)" );
#endif
#if !defined(CRONOS__WITH_EIGEN)
  if( options.DETERMINACY.BACKEND == Options::DET_EIGEN )
    return refuse( "DETERMINACY.BACKEND = DET_EIGEN", "Eigen" );
#endif
  return true;
}

inline bool
OCFESLV::_on_before_reduction
()
{
  // SOLVE.MARCHING: record the user's N-element evolution grid and collapse the
  // evolution domain to a single-element template BEFORE any structural stage (node
  // cache, classification, interface plan), so the assembled block has no inter-element
  // evolution interface -- decoupled by construction.  solve() then marches the recorded
  // grid by sliding this one element (reset_evolution_window).  Requires that the user
  // set the evolution domain via set_evolution_domain().
  // =====================================================================
  // rev92 PATH B PROBE -- hoist evolution-domain detection AHEAD of the collapse
  // =====================================================================
  // WHY.  handoff rev27 section 5 blamed a "resolution floor" for delta (B)'s coarse-mesh
  // damage.  That is contradicted by the marching data at the SAME NEL_T=3 on the SAME
  // model: IC_TRACE[march] |dC| = 3.493e-06 against IC_TRACE[mono] 4.593e-05 -- marching is
  // FULLY causal and 13x BETTER exactly where (B) is 26.8x worse.  So causality does not
  // cost accuracy below some resolution; what distinguishes them is that MARCHING HAS NO
  // SEAM (the grid collapses, the IC arrives as substituted data on a domain boundary,
  // exactly imposed) while (B) keeps the seam and weakens it to a one-sided variational
  // constraint carried by a tau.
  //
  // And since rev89, predicate (ii) IS structural marching availability -- so (B) can only
  // fire where marching could have run but did not.  Two such cases exist:
  //   * SOLVE.MARCHING=false, a deliberate user choice (OCFE_PDE0's mono configs are
  //     diagnostic, there to compare mono against march in one process);
  //   * the evolution domain was AUTO-DETECTED, where marching never runs because
  //     _collapse_evolution_for_marching() executes HERE, and _infer_evolution_domain()
  //     does not run until _classify_pde().  That is rev25's Defect B, and those ARE
  //     production configurations: k=3 on DAE5-auto and k=12 on MMPDE27-auto against k=0
  //     under --evolve, with MMPDE27-auto failing to converge where --evolve converges.
  //
  // The second case is an ORDERING ACCIDENT, not a modelling choice, and it does not need
  // (B): it needs marching. Hence Path B -- move the detection, do not approximate marching
  // inside the monolithic plan.
  //
  // THE STATED OBSTACLE DOES NOT HOLD.  rev86:11354 records that
  // "_infer_evolution_domain walks _mEqn for INITIAL equations pinned at a face, so it
  // depends on localisation having run".  True, and irrelevant: localisation HAS run --
  // _copy_usr_to_local() populates _mDom, _mVar and _mEqn and returns immediately above
  // this line.  The function consumes only eqn.opt->role, the lim mask, and domain NAMES,
  // all copied verbatim.  MEASURED, not argued: this probe replays its two rules here and
  // reports what they find.
  //
  // rev99: CRONOS_HOIST_EVODETECT DEFAULTS TO 1 (Path B promoted).  =0 is the opt-out.
  //                        =1          apply: set the evolution domain before the collapse,
  //                                    so marching becomes available on auto-detect models.
  //
  // PREDICTIONS, registered before running:
  //   H1 DAE5-auto / MMPDE27-auto / MESHONLY10-auto detect `t` here, and at level 1 collapse
  //      ( n_elem@plan=1 ), march, reach k=0, and match their --evolve counterparts.
  //   H2 the --evolve controls are UNCHANGED between levels: _evolution_dom_set is already
  //      true, so the probe does nothing on them.  If they move, the hoist is not inert
  //      where it should be and nothing else is interpretable.
  //   H3 a driver that detects nothing here detects nothing later either -- the two call
  //      sites read the same fields.  A driver that classifies with an evolution domain but
  //      reports would_change=n here would refute the whole ordering account.
  {
    static int const kHoist = []() -> int {
        // ===============================================================
        // rev106: DEFAULT 1.  PATH B RE-PROMOTED, on 138 rows.
        // ===============================================================
        // rev99 promoted this and rev105 reverted it, because the promotion rested on two
        // statements computed over 135 rows while presented as 139: "14 exposed, none
        // march" and "PDE9 is the only intersection".  OCFE_deepcopy, OCFE_fadbad and
        // OCFE_parallel abort under CRONOS_AUDIT_SETUP_ONLY, the sweep discarded their
        // rows, and OCFE_deepcopy -- exposed, hoist-applied -- regressed 0 -> 1 unseen.
        //
        // Both defects are fixed and the evidence is now complete:
        //
        //  * The sweep salvages rows from drivers that abort AFTER reporting (rev115), so
        //    the corpus reports crash=0, 138 usable, and the blast radius is 11 applied /
        //    92 already-set, not 10 / 90.
        //  * The exposure criterion is EXPOSED AND APPLIED, not "and marching" (rev116).
        //    OCFE_deepcopy regressed with marched 0->0: collapsing the evolution grid is
        //    sufficient to break a driver that inspects it; marching is not required.
        //    rev29 s1.1's criterion was reasoned, not measured, and it was too narrow.
        //  * OCFE_deepcopy sets SOLVE.MARCHING=false -- it is a STRUCTURAL test, and that
        //    option gates the collapse itself.  Verified: fail 0->0, agree 408/408.
        //  * OCFE_refinit, the other AT RISK row, was tested directly at both settings and
        //    passes.  It uses node_colloc to LOCATE values, not to verify a grid it built,
        //    and its own source says so.  Its grid collapses and it does not care.
        //
        // MEASURED on the complete set: 11 applied, 7 march, 0 regressions.  k
        // 45;105;225;465 -> 3;7;15;31 on MMPDE27, 3;7;15 -> 0 on DAE5 with nVar flat at 10
        // instead of 43/87/175, two non-convergences fixed, adjoint-vs-forward sensitivity
        // intact at round-off on both drivers that change solver path.
        //
        // CRONOS_HOIST_EVODETECT=0 remains the exact opt-out.  For a single driver that
        // must not be re-gridded, options.SOLVE.MARCHING=false is the finer-grained one and
        // is preferable: it needs no environment variable and travels with the source.
        return 1;   // CRONOS_HOIST_EVODETECT (retired 2026-10-07, WORKPLAN 3.C) (SOLVE.MARCHING=false is the per-driver alternative)
  }();
    _evoHoistLevel       = kHoist;
    _evoHoistApplied     = false;
    _evoHoistWouldChange = false;
    _evoHoistDom.clear();

    bool const already = ( _evolution_dom_set && _evolution_dom_var.dag() );
    if( !already ){
      // Non-mutating replay of _infer_evolution_domain()'s two rules, in its order and
      // with its exact comparisons ( lim == LB or lim == UB, not a bitmask test ).
      //
      // rev104: the replay records the PROVENANCE here.  rev99 added it to
      // _infer_evolution_domain and read it from this probe -- but this is a REPLAY, not a
      // call, so the member was still 0 when the probe printed and `infer=none` appeared on
      // all ten applied drivers.  Instrumenting one code path and reading it from another.
      FFVar const* found = nullptr;
      int prov = 0;                       // 1 = INITIAL-pinned (structural), 2 = name-based
      for( auto const& eqn : _mEqn ){
        if( eqn.opt->role != EqnRole::INITIAL ) continue;
        for( auto it = eqn.dom.cbegin(); it != eqn.dom.cend(); ++it )
          if( it->second == FFDom::LB || it->second == FFDom::UB ){ found = &it->first; prov = 1; break; }
        if( found ) break;
      }
      for( auto it = _mDom.cbegin(); it != _mDom.cend() && !found; ++it ){
        std::string const nm = it->first.name();
        if( nm == "t" || nm == "T" || nm == "time" || nm == "Time" ){ found = &it->first; prov = 2; }
      }
      _evoInferProvenance = prov;
      if( found ){
        _evoHoistWouldChange = true;
        _evoHoistDom         = found->name();

        // rev94 level 2: refuse unless the differential core is SQUARE.
        //
        // WHY THIS OPTION EXISTS.  rev94:9096 makes DIFFERENTIAL_RECTANGULAR marchable
        // only when the evolution domain is USER-SUPPLIED, and its comment argues the
        // case directly: "Admitting it unconditionally would let marching proceed on
        // unverified causality, and a wrong causal direction under marching gives a
        // SILENTLY wrong answer."  Hoisting removes exactly that guard.
        //
        // The comment's EVIDENCE is now refuted -- it cites MMPDE27 auto failing while
        // --evolve converges, which rev92 showed is the ordering accident, not unverified
        // causality: hoisted, auto converges to the same answer as --evolve on every
        // metric.  The comment's CONCERN is not refuted.  _infer_evolution_domain picks
        // the first INITIAL-role equation pinned at a face, else a domain literally NAMED
        // t/T/time/Time.  That is a heuristic, and a heuristic that is wrong turns a slow
        // monolithic solve into a silently wrong marched one.
        //
        // The classifier cannot be consulted here -- _classify_pde runs AFTER the
        // collapse, which is the same ordering that caused the defect.  So the shape test
        // is replicated on user-side structures instead.
        //
        // NOTE WHERE THIS LANDS, measured before shipping it: OCFE_MMPDE27 is sym=3x4 and
        // OCFE_DAE5 is sym=1x2 -- BOTH rectangular.  Level 2 therefore refuses the hoist
        // on both models it fixed at level 1.  This option is here to be MEASURED against
        // level 1, not because it is the recommendation.
        // ---------------------------------------------------------------------
        // rev95 FIX -- the shape could not be measured before the domain was set.
        // ---------------------------------------------------------------------
        // _differential_core_shape reaches _time_derivative_state, whose first line is
        //   if( !op || !op->sameid(typeid(FFPartial)) || !_evolution_dom_set
        //    || !_evolution_dom_var.dag() ) return FFVar();
        // rev94 called it BEFORE assigning those two, five lines below.  On the
        // auto-detect path -- the only path that reaches here -- they are unset, so it
        // returned rows=0, cols=0 UNCONDITIONALLY.
        //
        // MEASURED consequence: core=0x0 on all 139 drivers of the corpus sweep.  An
        // invariant that reads identically on every model in a corpus is not a
        // measurement, and that uniformity was visible in the output before anything was
        // concluded from it.  Level 2's square-core gate had no data at all and would
        // have refused every model for the same reason, which would have looked like a
        // finding about the corpus rather than about the probe.
        //
        // Set the candidate, measure, restore if not applying -- so level 0 stays inert.
        FFVar const evo_save     = _evolution_dom_var;
        bool  const evo_set_save = _evolution_dom_set;
        _evolution_dom_var = *found;
        _evolution_dom_set = true;

        // rev103: the level-2 SQUARE-CORE GATE is removed; the SHAPE is kept.
        //
        // rev94 added a gate that hoisted only when the differential core was square,
        // replicating the collapse site's marchability rule ahead of the classifier.  It
        // was measured and it protected the wrong set: rectangular is exactly the shape
        // Path B helps (DAE5 1x2, k 3;7;15 -> 0; MMPDE27 3x4), while OCFE_PDE9 -- the one
        // driver in 139 that needed care, being the sole post-solve eval() gate the hoist
        // reaches -- is sym=3x3 SQUARE and the gate would have admitted it.
        //
        // _differential_core_shape stays.  core=RxC is what made DAE5's rectangular core
        // legible and is how the orphaned tau columns were localised (handoff rev29 s2).
        // A reported quantity, not a decision.
        _differential_core_shape( _evoCoreRows, _evoCoreCols );

        if( kHoist == 1 ){
          _evoHoistApplied = true;   // domain assigned above; keep it
          // _evolution_dom_user stays FALSE on purpose: the provenance is still
          // auto-detection, only the TIMING moved.  Reporting it as user-supplied would
          // hide exactly the distinction this probe exists to measure.
        }
        else{
          _evolution_dom_var = evo_save;      // restore -- level 0 must be
          _evolution_dom_set = evo_set_save;  // byte-for-byte inert
        }
      }
    }

    static bool const kEvoFaceH = Options::_env_flag( "CRONOS_AUDIT_EVOFACE", false );
    if( kEvoFaceH )
      _disp( 3 ) << "  [evoface] rev92 hoist probe: hoist level=" << kHoist
                << "  domain_already_set=" << ( already ? "y (user-supplied; probe inert)" : "n" )
                << "  inferable_here=" << ( _evoHoistWouldChange ? _evoHoistDom : std::string("<none>") )
                << "  core=" << _evoCoreRows << "x" << _evoCoreCols
                << "  infer=" << ( _evoInferProvenance == 1 ? "INITIAL-pinned"
                                 : _evoInferProvenance == 2 ? "NAME-BASED (weak)"
                                                            : "none" )
                << ( _evoCoreRows && _evoCoreRows == _evoCoreCols ? " (square)" : " (RECTANGULAR)" )
                << "  APPLIED=" << ( _evoHoistApplied ? "YES -> marching may now collapse"
                                                      : "no" )
                << "\n";
  }

  _collapse_evolution_for_marching();   // no-op unless marching is available


  // A local working DAG and local model are now available. Mark the model as
  // set up before order reduction, since reduce_order() is a protected local
  // transformation and requires a valid local setup model.
  return true;
}


inline bool
OCFESLV::_on_before_closures
()
{
  // Validate user-supplied incoming BCs at hyperbolic boundary faces (count +
  // characteristic direction) before generating the outgoing closure, so a
  // mis-placed inflow condition fails setup cleanly rather than producing a
  // rank-deficient solve.  Always on (unconditional since rev319).
  if( _classified && !_validate_hyperbolic_incoming_bcs() ){
    std::cerr << "OCFESLV::setup ** hyperbolic incoming-BC validation failed "
                 "(see messages above); not square / well-posed for this "
                 "classification\n";
    _issetup = false;
    // The validator sets HYP_BC_MISDIRECTED on a criterion #3 direction reject;
    // fall back to HYP_INCOMING_BC (count mismatch) only if it left status OK.
    if( _setupStatus == SetupStatus::OK )
      _setupStatus = SetupStatus::HYP_INCOMING_BC;
    return false;
  }
  return true;
}

inline bool
OCFESLV::_on_model_ready
()
{
  // the collocated layout's consistency count, dependency maps, eval plans, SAT coupling, the interface
  // plan (draft, validate, commit), the collocation buffers, the tau report and the FFDep audit of the
  // strong interface rows: OCFEBASE::_build_colloc_plan, shared with OCFEMOL (MoL phase 3, P3d)
  if( !_build_colloc_plan() ) return false;

#ifdef CRONOS_DEBUG_OUT
  _disp( 3 ) << "setup: before deriv cache rows=" << _nCollEqn << " fct=" << _nCollFct << "\n";
#endif
  // Derivative sparsity/coloring can be built eagerly for the legacy
  // behaviour, or lazily on the first derivative/sparsity query.
  _jacReady = false;
  _linEqnReady = _linFctReady = false;
  if( options.SOLVE.DERIV_CACHE == Options::DERIV_EAGER
   || options.SOLVE.LINEAR_EVAL_CACHE == Options::LINEAR_EAGER ){
    if( !_setup_deriv_cache() )
      { _setupStatus = SetupStatus::DERIV_CACHE_FAILED; return false; }
    _reset_autodiff( _FVar, _FInp, _FEqn, _FFct );
    if( options.SOLVE.LINEAR_EVAL_CACHE == Options::LINEAR_EAGER
     && !_setup_linear_cache( true, true, _vCstVal.empty()? nullptr: _vCstVal.data() ) )
      { _setupStatus = SetupStatus::LINEAR_CACHE_FAILED; return false; }
  }
  _phase( "deriv_cache (if eager)" );

  // Mark setup complete BEFORE the audit: _audit_reduced_dof() -> _ensure_deriv_cache()
  // requires _issetup (it builds the reduced jac-pattern on demand).
  _issetup = true;

  // Stage-4: post-index-reduction structural DOF audit.  A structural rank
  // deficiency in the reduced system is a genuine singularity (structural rank
  // bounds numerical rank), so by default setup() rejects the model; the
  // FATAL.REDUCED_DOF=false opt-out downgrades it to a warning.
  _audit_reduced_dof();
  _phase( "audit_reduced_dof" );

  // 2026-10-10: invalidate the PREVIOUS setup's caches HERE, before the W-test (which recurses into setup()) and
  // the determinacy audit can read them.  They were invalidated at the END of this hook, after both: a first setup
  // never noticed (nothing cached yet), but on a re-setup the window audit ran with the previous plan's marching IC
  // override rows (_marchICRows, indexed into the OLD plan) and found a deficiency the new plan does not have --
  // the W-test then REVERTED a correct restore (OCFE_MMPDE1, adaptive monitor: "worst excess 9 >= 3").
  _primalValid = false;   // model (re)built -> any cached primal solution is stale
  _reducedDepValid = false;   // ... and the structural dependency pattern
  _rowProvenanceReady = false; _rowEqnNdx.clear();   // ... and the residual-row provenance
  _eqnDiffInReady = false; _eqnDiffIn.clear();       // ... and the row x direction Partial table
  _marchICRowsReady = false; _marchICRows.clear(); _marchICTerminal.clear(); _marchICColToNode.clear(); _marchICActive = false;

  // rev133: numerical determinacy audit on the assembled Jacobian at a deterministic
  // reference point.  Unconditional on imposition mode and on index reduction (the
  // structural audit above runs only on reduced IC_WEAK systems).  Report-only unless
  // FATAL.DETERMINACY.
  // rev139 (uniform-rule step 1): test the dropped claims on the built system and
  // restore the needed subset via one revert-guarded re-derive.  When the restore
  // recursion runs, the INNER setup has already audited and printed the FINAL
  // system, so the audit here is skipped (exactly one [determinacy] line set per
  // instance either way).
  if( _wdecide_restore() ){
    _phase( "wdecide_restore" );
  }
  else{
    _audit_determinacy( nullptr, nullptr );
    _phase( "audit_determinacy" );
    // C4 (rev222): the STRUCTURAL determinacy.  With the numeric audit OFF this line stands in for it (same
    // shape, tagged STRUCTURAL); with it ON and CRONOS_KPRIME_STRUCTURAL=1 the numbers are printed alongside,
    // which is the corpus-wide equality check.
    // rev225: printed by DEFAULT.  With the audit off this is the only determinacy record in the log, and the
    // sweep tooling parses these lines; the knob now only forces it ALONGSIDE the audited line.
    static bool const kStruct = Options::_env_flag( "CRONOS_KPRIME_STRUCTURAL", false );
    if( !options.DETERMINACY.AUDIT || kStruct ){
      size_t n_implied = 0, n_states = 0;
      auto const cm = claim_classification();
      for( auto const& [s, byd] : cm ){ ++n_states;
        for( auto const& [d, c] : byd )
          if( c.tag == ClaimTag::IMPLIED || c.tag == ClaimTag::CLOSURE ) ++n_implied; }
      // k_struct = m - rank(B) from the trace projection: MEASURED equal to the audited k' on 373/373 fresh
      // IC_STRONG pairings (the shadow deficit: 368/373, over-reporting because it precedes the explicit-tau
      // hybrid).  Computed in both exact modes; printed here for both, which is the TRACE evidence the earlier
      // analysis lacked.  The shadow deficit stays as the CLAIM-level diagnostic C0 predicts.
      _disp( 1 ) << "[determinacy] STRUCTURAL mode=" << ( options.INTERFACE.IMPOSITION == Options::IC_WEAK? "IC_WEAK":
                                                        options.INTERFACE.IMPOSITION == Options::IC_STRONG? "IC_STRONG": "IC_TRACE" )
                << " k_struct=" << _traceProjK
                << " proj_applied=" << ( _traceProjApplied? 1: 0 )
                << " k_pred=" << ( _traceProjApplied ? size_t(0) : _traceProjK )   // rev281: applied in either exact mode
                << " deficit=" << _structDeficit << " deficient_clusters=" << _structDeficientClusters
                << " implied_entries=" << n_implied << " over " << n_states << " state(s)"
                << ( options.DETERMINACY.AUDIT? "  (audit ON: compare with the audited defR above)":
                                                "  (audit OFF: this line replaces the SPQR audit)" ) << std::endl;
      // Fill the record too, when the numeric audit did not: an EMPTY _determinacy reads as "no deficiency"
      // to any determinacy_audit() caller, which is the "ran and found nothing" vs "never ran" trap again.
      // gap / sigma_cut / tol_abs are left ZERO: they come from a factorisation and are not fabricated here.
      if( !options.DETERMINACY.AUDIT ){
        t_Determinacy d;
        d.ran      = true;
        d.window   = 0;
        d.mode     = ( options.INTERFACE.IMPOSITION == Options::IC_WEAK? "IC_WEAK":
                       options.INTERFACE.IMPOSITION == Options::IC_STRONG? "IC_STRONG": "IC_TRACE" );
        d.backend  = "STRUCTURAL";
        d.ref_point= "none (structural)";
        d.k_prime  = ( options.INTERFACE.IMPOSITION == Options::IC_TRACE && _traceProjApplied )? size_t(0): _traceProjK;
        d.def_right= d.k_prime;
        d.n_tau    = _traceProjNTau;
        _determinacy.push_back( d );
      }
    }
  }


  // rev142 (FAIL-CLOSED).  FATAL.DETERMINACY asks setup() to REJECT a model the audit
  // did not certify.  Through rev141 the test was `D.ran && verdict == STATE CONTENT`,
  // so a window that never ran passed the gate in silence -- and rev141 is exactly the
  // revision that introduced a way for windows not to run (CRONOS_AUDIT_MAX_N).  A
  // large model under the budget would sail through the fatal option with no rank ever
  // computed, which is the one case the option exists to stop.  A MISSING verdict is
  // now a failure, with the reason named on stderr.  The same rule catches "assembly
  // failed" and "no rank backend": each of those is an absent certificate, not a pass.
  // 2026-10-10: a MEASURED state null space is said at display level 1 whatever FATAL.DETERMINACY is.  Until now a
  // model whose audit reported "** STATE CONTENT **" (a solution family, not a solution) set up with status OK and
  // nothing at the default display: OCFE_MMPDE1 with an adaptive monitor, where a wrongly reverted restore had
  // left rank 957 of 960 in the states.  FATAL.DETERMINACY (below) turns this into a refusal.
  for( auto const& D : _determinacy )
    if( D.ran && D.verdict == t_Determinacy::STATE_CONTENT )
      _disp( 1 ) << "OCFESLV::setup ** determinacy: window " << D.window << " has a NULL SPACE IN THE STATES"
                 << " (rank " << ( D.n > D.def_right ? D.n - D.def_right : 0 ) << " of " << D.n << "): the discrete"
                 << " problem has a family of solutions, not one -- setup continues (set FATAL.DETERMINACY to refuse)\n";
  // Corpus-inert -- no driver sets the option and the default is off.
  if( options.FATAL.DETERMINACY ){
    std::string why;
    if( !options.DETERMINACY.AUDIT )
      why = "FATAL.DETERMINACY set but DETERMINACY.AUDIT is off -- no verdict to act on";
    else if( _determinacy.empty() )
      why = "no determinacy record: the audit did not run";
    else{
      for( auto const& D : _determinacy ){
        if( !D.ran ){
          why = "window " + std::to_string( D.window ) + " NOT AUDITED ("
              + ( D.skipped_reason.empty() ? std::string( "no reason recorded" ) : D.skipped_reason ) + ")";
          break;
        }
        if( D.verdict == t_Determinacy::STATE_CONTENT ){
          why = "window " + std::to_string( D.window ) + " ** STATE CONTENT **";
          break;
        }
      }
    }
    if( !why.empty() ){
      std::cerr << "[determinacy] FATAL: " << why << std::endl;
      _setupStatus = SetupStatus::DETERMINACY_UNDETERMINED;
      _issetup = false;
      return false;
    }
  }
  if( _dofAudit.ran && !_dofAudit.ok() && options.FATAL.REDUCED_DOF ){
    _setupStatus = SetupStatus::REDUCED_DOF_INCONSISTENT;
    _issetup = false;
    return false;
  }

  // Stage-5: evolution-domain suitability + consistent-IC validation.  Warn-by-
  // default (returns true); fatal only under FATAL.EVOLUTION_DOMAIN.
  if( !_validate_evolution_domain() ){
    _setupStatus = SetupStatus::EVOLUTION_DOMAIN_UNSUITABLE;
    _issetup = false;
    return false;
  }
  _phase( "validate_evolution_domain" );

  // Auto-register inputs flagged is_decision at add_input() (declaration order).  Rebuilt on each
  // setup() so it is idempotent under re-setup; an explicit register_control() after setup() still
  // appends (deduplicated).  Skipped implicitly on the early-return setup failures above.
  _reindex_controls();   // layout from _decisionFlag; FFModel::input_ndof reads the declaration
  // (the previous setup's caches were invalidated ABOVE, before the W-test and the determinacy audit -- 2026-10-10)

  // Header identity, once per process.  A stale object file not relinked against a new
  // header produces a log indistinguishable from a correct run of the OLD header -- the
  // new diagnostics simply do not appear, which reads as "the instrument found nothing"
  // rather than "the instrument is not there".  That has cost this project a session.
  //
  // rev35: UNCONDITIONAL, once per process.  The stamp is a BUILD ORACLE, not a diagnostic,
  // and it has now caught five stale binaries in this arc -- including one that produced two
  // byte-identical logs which looked like a valid before/after comparison.  Gating it behind
  // DISPLAY_LEVEL>=1 meant that marched drivers (PSA9 sets DISPLAY_LEVEL=0) and any run
  // without CRONOS_AUDIT had NO build identification at all, so `make` skipping a relink --
  // the -DOCFE_OCFESLV_HEADER flag is not tracked as a dependency -- was undetectable.
  //
  // One line per process, no flags, no measurable cost.  Do not re-gate it: the whole point
  // is that `grep -c "<expected rev>"` is a valid guard on a log from ANY run.
  {
    static bool id_stamped = false;
    if( !id_stamped ){
      id_stamped = true;
      // rev67: the tag alone is the build oracle.  CRONOS_HEADER_FULL=1 appends the
      // revision summary for anyone who wants it in the log.
      // 2026-10-04: at DISPLAY_LEVEL >= 2 -- out of normal output and the tutorials.  The revisions are queried
      // with revision() (FFModel::revision(), OCFESLV::revision()), and the corpus driver OCFE_revision prints
      // them once per sweep log.
      _disp( 2 ) << "OCFESLV ** header: " << HEADER_ID;
      // The model layer carries its own stamp and a binary can mix the two (a sweep once ran a rev208
      // ffmodel.hpp under a rev207 solver stamp).  Printed only when the revisions DIFFER, so that matched
      // builds keep their output byte-identical and a mismatch cannot pass unnoticed.
      {
        auto rev = []( char const* s ){
          std::string t( s ); auto p = t.find( "rev" );
          if( p == std::string::npos ) return std::string();
          auto q = t.find( ' ', p );
          return t.substr( p, q == std::string::npos? q: q-p );
        };
        if( rev( HEADER_ID ) != rev( FFModel::HEADER_ID ) )
          _disp( 2 ) << "  [ffmodel: " << FFModel::HEADER_ID << "]";
      }
      if( Options::_env_flag( "CRONOS_HEADER_FULL", false ) )
        _disp( 2 ) << "  " << HEADER_SUMMARY;
      _disp( 2 ) << std::endl;
    }
  }

  // rev67: CLASSIFICATION STAMP.  Provenance, not diagnostics -- stamped unconditionally,
  // once per block, on the same footing as HEADER_ID.
  //
  // An earlier rev67 draft put this inside the evolution validator, which EARLY-RETURNS
  // when no evolution domain is in play (`if( !_classified || !_evolution_dom_set ... )
  // return true;`).  OCFE_AE0 is zero-domain, so it printed no verdict at all -- the exact
  // invisibility this stamp exists to remove.  Emitting it here covers every model.
  //
  // Why unconditional: drivers advertise "gates: setup square + classify OK", but under
  // CLASS_AUTO a failed classification proceeds silently and setup() still returns true,
  // so that gate was VACUOUS and unverifiable from a log.  The only prior emission site
  // was gated on the compile-time macro MC__OCFESLV_SYMBOL_BALANCE_PROBE, which no build in
  // the corpus defines.
  //
  // The EMPTY case is reported explicitly rather than passed over in silence.  The map is
  // filled by walking _blockSymbol, and a model with no collocation domains may never
  // build a principal symbol at all -- printing nothing there would reproduce the exact
  // failure this stamp removes, and "no line" is indistinguishable from "stale binary".
  // ---------------------------------------------------------------------------
  // rev86 DEFECT B WARNING -- the collapse-point divergence, made diagnosable
  // ---------------------------------------------------------------------------
  // MEASURED (handoff rev22 section 0).  The interface plan is built over a DIFFERENT GRID
  // depending on HOW the evolution domain became known, and nothing said so:
  //
  //     NEL_T          2      4      8
  //     --evolve      81     81     81      <- independent of the t-mesh
  //     auto-detect  260    606   1298      <- 173*NEL_T - 86, linear in it
  //
  // set_evolution_domain() records the domain at localisation (rev86:10651), BEFORE the
  // plan is built, so the t-direction is collapsed and the plan is constructed for ONE
  // t-window and reused: 81 claims, all dir=xi.  _infer_evolution_domain() assigns it at
  // 15315/15328, called from _classify_pde at 10943 -- AFTER that collapse -- so the plan
  // is built over the FULL SPACE-TIME TENSOR GRID: ~87 claims per t-element plus ~86 per
  // interior t-seam.
  //
  // Consequences, all measured: k = nel_xi*(NEL_T-1) cross-domain kernel directions (12 at
  // nel_xi=4, NEL_T=4), n_trace_var 606 vs 81, and conv=no at nel_xi=8 where the
  // user-supplied path converges.
  //
  // FIVE mechanisms were proposed from reading the source and refuted by measurement before
  // the NEL_T sweep found this one; see rev22 section 2.  In particular _evolution_dom_set
  // is TRUE on both paths, _evolution_dom_user has no consumer in the plan, and the causal
  // interface block is skipped in both (evolution_input_may_jump=FALSE, n_inputs=0).
  //
  // THE FIX IS NOT "set the flag earlier".  The collapse is a STRUCTURAL operation on the
  // domain grid; making auto-detection trigger it means moving DETECTION ahead of the
  // collapse, and _infer_evolution_domain walks _mEqn for INITIAL equations pinned at a
  // face, so it depends on localisation having run.  Options (a) detect earlier on the
  // user-side structures and (b) re-run the collapse after detection are both assessed in
  // handoff rev22 section 3 and both carry corpus-wide behavioural risk.
  //
  // This is option (c): NO behavioural change, just say so.  A silent 7.5x claim inflation
  // becomes a diagnosable one, and the caller is told the remedy.
  if( _evolution_dom_set && _evolution_dom_var.dag() && !_evolution_dom_user )
    _disp( 1 ) << "OCFESLV::setup [evolution] NOTE: the evolution domain '"
              << _evolution_dom_var.name() << "' was AUTO-DETECTED, which happens AFTER the"
                 " localisation collapse, so the interface plan is built over the FULL"
                 " space-time tensor grid rather than one collapsed evolution window."
                 "  Measured on OCFE_MMPDE27: 173*NEL_T-86 continuity claims instead of a"
                 " constant 81, a cross-domain kernel of nel_xi*(NEL_T-1), and loss of"
                 " convergence at mesh sizes where the collapsed plan converges."
                 "  Call set_evolution_domain( " << _evolution_dom_var.name() << " ) to get"
                 " the collapsed plan." << std::endl;

  {
    if( _blockClassification.empty() )
      _disp( 3 ) << "OCFESLV::setup [classify] no classified block: "
                << ( _vDom.empty()
                       ? "zero collocation domains (ALGEBRAIC_LUMPED by construction)"
                       : "no principal symbol was built" )
                << std::endl;
    else
      for( auto const& [bid, bc] : _blockClassification ){
        // The VERDICT alone is not enough to check.  OCFE_DAE0 -- five states, four
        // evolution equations and one algebraic constraint with no OpP, i.e. an
        // unambiguous index-1 DAE -- reported DIFFERENTIAL_ORDINARY, and nothing in the
        // log said why: the only line carrying sym=RxC and Ae_rank is the
        // _symbol_balance_probe output, gated on a compile-time macro no build defines.
        //
        // The ODE/DAE split is decided by ONE test, `res.type = ri.singular ? DAE : ODE`,
        // where ri is the rank of the evolution matrix A_e.  So A_e's singularity, its
        // smallest singular value and the tolerance that judged it are exactly the
        // evidence needed to tell "the rank test saw the singularity and the verdict was
        // overwritten downstream" from "the algebraic row never reached A_e".  They cost
        // one line and they travel with the verdict.
        // Symbol DIMENSIONS travel with the verdict too.  Measured on OCFE_DAE0 -- five
        // states, four evolution equations, one algebraic constraint with no OpP, so A_e
        // must be diag(1,1,1,1,0) and sigma_min must be 0 -- the stamp reported
        // sigma_min = sigma_max = 1, i.e. A_e is an IDENTITY and the algebraic row is not
        // in it.  Without RxC there is no way to tell whether the row was dropped from
        // the symbol, or the symbol is smaller than the model, or a square-ing pass
        // removed it.  ODE0 and DAE0 print identical numbers, so the dimensions are the
        // only thing that can separate them.
        auto itsym = _blockSymbol.find( bid );
        size_t const symR = ( itsym != _blockSymbol.end() ) ? itsym->second.vEqn.size()   : 0;
        size_t const symC = ( itsym != _blockSymbol.end() ) ? itsym->second.vState.size() : 0;
        size_t const symD = ( itsym != _blockSymbol.end() ) ? itsym->second.vDom.size()   : 0;
        _disp( 1 ) << "OCFESLV::setup [classify] blk" << bid << " type="
                  << pde_type_name( bc.type )
                  << "  sym=" << symR << "x" << symC << " nDom=" << symD
                  << " nCoeff=" << ( itsym != _blockSymbol.end()
                                       ? itsym->second.vCoeff.size() : 0 )
                  << "  Ae_singular=" << ( bc.At_singular ? "y" : "n" )
                  << " sigma_min_evol=" << bc.sigma_min_evolution
                  << " sigma_max_evol=" << bc.sigma_max_evolution
                  << " imag_tol=" << bc.imag_tol
                  << " evol_dom_idx="
                  << ( bc.evolution_dom_idx == std::numeric_limits<size_t>::max()
                         ? std::string( "none" ) : std::to_string( bc.evolution_dom_idx ) )
                  << ( bc.evolution_user_supplied ? " [user-supplied]" : "" )
                  << ( bc.evolution_auto_detected ? " [auto-detected]" : "" )
                  << std::endl;
      }
  }

  if( options.AUDIT.EQN_LEGEND ) print_equation_legend( std::cout );
  if( options.AUDIT.COVERAGE ){
    print_column_coverage_report( std::cout );
    print_tau_pivot_report( std::cout );
  }

  // ---------------------------------------------------------------------------
  // rev89 [evoface] PLAN METRIC -- the numbers handoff rev25 section 3.3 is actually about
  // ---------------------------------------------------------------------------
  // rev87 and rev88 reported claims and receiver visits; NEITHER reported n_trace_var or
  // k, which is what "0 taus against 6, k=0 against nel_t-1" refers to.  On OCFE_DAE5 the
  // claim census was IDENTICAL at switch levels 0 and 2 while six receivers were dropped,
  // so the whole effect was invisible to every number then being captured.  Emitted from
  // the COMMITTED plan, so it is the value the solver will use, not a draft field.
  //
  // interface_rank_symmetry() is two sparse Gram-Schmidt passes (rev7) and returns early
  // under IC_WEAK, where there are no multipliers.  Gated on CRONOS_AUDIT_EVOFACE so it
  // costs nothing in a normal run.
  {
    static bool const kEvoFace = Options::_env_flag( "CRONOS_AUDIT_EVOFACE", false );
    if( kEvoFace ){
      size_t _rB_e = 0, _rC_e = 0;
      bool const _sym_e = interface_rank_symmetry( nullptr, &_rB_e, &_rC_e );
      size_t const _ntv = _plan.n_trace_var();
      _disp( 3 ) << "  [evoface] PLAN METRIC (committed): n_trace_var=" << _ntv
                << "  rank(B)=" << _rB_e << "  rank(C)=" << _rC_e
                << "  k=" << ( _ntv > _rB_e ? _ntv - _rB_e : 0 )
                << "  symmetric=" << ( _sym_e ? "y" : "n" )
                << "  imposition="
                << ( options.INTERFACE.IMPOSITION == Options::IC_WEAK   ? "IC_WEAK"
                   : options.INTERFACE.IMPOSITION == Options::IC_TRACE  ? "IC_TRACE"
                   : options.INTERFACE.IMPOSITION == Options::IC_STRONG ? "IC_STRONG" : "?" )
                << std::endl;
    }
  }

  // rev7: saddle-point rank symmetry (B vs C).  Report-only by default; fatal only
  // under FATAL.INTERFACE_RANK, because a column-rank-deficient B makes the augmented
  // system singular WITHOUT necessarily making it unsolvable -- the residual may still
  // satisfy the compatibility conditions.  Reached whenever the plan built, which is
  // exactly the silent case: the loud one already fails inside the tau elimination.
  if( options.AUDIT.RANK_SYMMETRY || options.FATAL.INTERFACE_RANK ){
    size_t _rB = 0, _rC = 0;
    bool const _sym_ok = interface_rank_symmetry(
      options.AUDIT.RANK_SYMMETRY ? &std::cerr : nullptr, &_rB, &_rC );
    if( !_sym_ok && options.FATAL.INTERFACE_RANK ){
      std::cerr << "OCFESLV::setup ** interface rank symmetry: the multiplier column "
                   "block B has rank " << _rB << " < n_trace_var "
                << _plan.n_trace_var() << " (continuity row block C has rank " << _rC
                << ").  The augmented system is singular for any C; the plan commits "
                   "to multipliers the receiver graph cannot represent.  Set "
                   "options.FATAL.INTERFACE_RANK = false to proceed anyway."
                << std::endl;
      _setupStatus = SetupStatus::INTERFACE_PLAN_INVALID;
      _issetup = false;
      return false;
    }
  }

  // Structural-sweep abort.  Everything the audits report is frozen by this point, so a
  // corpus sweep has no reason to continue into the caches and the solve.  Returns false;
  // the audit lines above are the deliverable.
  //
  // rev38 tried std::exit(0) here to stop three drivers (OCFE_deepcopy, OCFE_fadbad,
  // OCFE_parallel) THROWING on the false return and dying on SIGABRT.  It fixed the crashes
  // and BROKE THE SWEEP: many drivers construct several OCFESLV instances -- sweeping
  // IC_WEAK/IC_TRACE/IC_STRONG, or several meshes -- and exiting at the first setup() call
  // killed the rest.  Corpus rows fell 597 -> 135 and SYMMETRIC 257 -> 15, because the first
  // configuration in most drivers is IC_WEAK and reports nothing structural.
  //
  // Reverted deliberately.  The three aborts are a KNOWN, HARMLESS artefact of this flag:
  // they occur AFTER the audit output, so those drivers still contribute their rows, and
  // rc=134 with a core dump is noise rather than data loss.  Do not "fix" it by exiting the
  // process.  If it is worth fixing at all, the shape is to keep `return false` and force
  // the exit status separately, which needs care around drivers that throw.
  //
  // Consequence to remember when reading a sweep: under this flag rc is NOT a failure
  // signal -- 101 drivers exit rc=1 by design and three abort.  Judge stage 1 by the rows
  // produced, and use the full corpus sweep (no flag) for pass/fail.
  // rev144: the w-decide nested setups (restore attempt and revert) are EXEMPT from
  // this stop.  MEASURED (rev143 [ref-audit] on PDE6): the stop's `_issetup=false;
  // return false;` made the restore re-derive fail BY CONSTRUCTION under audit sweeps
  // (ok=false for every model, reference innocent), and because a reverting
  // _wdecide_restore() returns the revert setup's ok, the outer setup then took the
  // re-audit branch with _issetup clobbered -> init() refused -> zeros reference ->
  // "assembly failed" SKIPPED lines.  Letting the nested setups complete restores the
  // solve-mode w-decide semantics under CRONOS_AUDIT_SETUP_ONLY; the OUTER setup
  // (retry flag cleared by then) still stops here, so sweep behaviour is otherwise
  // unchanged.  det_expect_setuponly-changing for PDE6/PDE6_blk0/MESHONLY10 by
  // design: their entries move to the solve-baseline verdicts.  Re-emit after gate.
  if( _knob_AUDIT_SETUP_ONLY() && !_wdecideRetry ){
    std::cerr << "OCFESLV::setup ** AUDIT_SETUP_ONLY: stopping after the audit reports "
                 "(CRONOS_AUDIT_SETUP_ONLY is set).  This is not an error." << std::endl;
    _issetup = false;
    return false;
  }

  if( _solve_display_level() >= 2 ){
    double const tot = std::chrono::duration<double,std::milli>(
                         std::chrono::steady_clock::now() - _t_setup0 ).count();
    _disp( 2 ) << "OCFESLV::setup [timing] " << std::left << std::setw(30) << "TOTAL"
              << std::right << std::fixed << std::setprecision(2)
              << std::setw(11) << tot << " ms" << std::defaultfloat << "\n";
  }

  return _issetup;
  return _issetup;
  return _issetup;
}


inline bool
OCFESLV::_ensure_deriv_cache
()
{
  if( !_issetup ) return false;
  if( _jacReady
   && _jacInterfaceType  == options.INTERFACE.TYPE
   && _jacImpositionType == options.INTERFACE.IMPOSITION
   && _jacClassified     == _classified )
    return true;

  if( !_setup_deriv_cache() )
    return false;

  // Rebuild the compressed FAD seed arrays and the pointer-bearing OCVar
  // wrappers because the color count may have changed.
  _reset_colloc( _OCFVar, _OCFInp, _OCFDom, _OCFCst );
  _OCFEqn = OCVar<t_FAD>();
  _wrkOCF.clear();
  _reset_autodiff( _FVar, _FInp, _FEqn, _FFct );

  // rev122: make display_interface_plan_diagnostics() reachable WITHOUT a driver edit.
  //
  // It carries the per-state SAT receiver-edge audit
  //     state : edges | sat==0 | ->AUX | ->GEO | LEFT-ZERO | max|sat|
  // which is the only thing that says whether a claim's coupling came from the PRINCIPAL
  // SYMBOL or from one of the two rescues that FORCE coupling = 1.  That distinction is
  // load-bearing: on OCFE_PDE30b the [taucol/plan] trace shows coupling=1.000000e+00
  // exactly on one receiver row, which is either a genuine unit coefficient or a rescue,
  // and nothing else in the output separates them.
  //
  // Only 8 of 139 drivers call it, and none of the ones carrying k>0 that we have been
  // reading.  A diagnostic that needs a driver edit is not usable at corpus scale -- the
  // same reason CRONOS_SOLVE_VERBOSE and CRONOS_TRACE_PROJ_EPS exist.
  //
  // To std::cerr, so it does not interleave with a driver's own std::cout table.
  if( Options::_env_flag( "CRONOS_AUDIT_IFACEPLAN", false ) )
    display_interface_plan_diagnostics( std::cout );

  // -------------------------------------------------------------------------------
  // rev127: PER-STATE ROLE -- differential (transferred) vs algebraic (recovered).
  // -------------------------------------------------------------------------------
  // The one quantity the MANUAL/AUTO divergence implies and NOTHING currently prints.
  //
  // OCFE_PDE6_blk0 writes the same parabolic problem two ways.  Their claim censuses
  // differ by 100 in the EVOLUTION direction:
  //     MANUAL  q        t-claims = 0     x-claims = 100
  //     AUTO    Daux3_T  t-claims = 100   x-claims = 84
  // and MANUAL's monolithic solve is PRIMAL UNDETERMINED while AUTO's is not.
  //
  // Five candidate causes have been eliminated BY MEASUREMENT, not by argument:
  //     domain specs      -- identical, [t] x [x] on both LINK rows (rev125 dom= field)
  //     registration order -- --linklast: census, rank and answer all unchanged
  //     _auxVarID routing  -- inert here; [symrow] shows value-slaved states = 0, so the
  //                           (b2) branch that consults it is never reached
  //     boundary-cond kind -- --bcuderiv: rank_deficiency = 4 either way
  //     the RECTANGULAR warning -- an instrument fault, fixed in rev126
  //
  // _states_with_evolution_derivative IS the differential-state identification -- the
  // header's own words, "evaluated wherever needed (globally for marching transfers,
  // per-block in the decomposition)".
  //
  // MEASURED, and it REFUTES the expectation this was built to test:
  //     T(t,x)         role=DIFFERENTIAL  in_auxVarID=n   t-claims =  92
  //     q(t,x)         role=ALGEBRAIC     in_auxVarID=n   t-claims =   0
  //     Daux3_T(t,x)   role=ALGEBRAIC     in_auxVarID=y   t-claims = 100
  //
  // Daux3_T is ALGEBRAIC -- the SAME role as q -- and still carries 100 evolution-face
  // claims where q carries none.  So the differential/algebraic identification is not
  // what drives the census difference.  in_auxVarID still separates them, so aux-ness
  // remains the only property known to differ, but it reaches the claim count by some
  // path other than this role and other than the (b2) routing branch already shown
  // inert by [symrow].
  //
  // It also breaks the marching-reference reading (see the [evoface] MARCHING REFERENCE
  // line): under that reading an ALGEBRAIC state should carry NO evolution-face claim,
  // yet AUTO -- the DETERMINATE and more accurate formulation -- gives one 100 of them,
  // while MANUAL, which obeys the reference, is the undetermined one.  The reference is
  // a correct description of what MARCHING does; it is NOT a correct predictor of what
  // the monolithic plan should mint.  Do not treat it as an invariant.
  //
  // Also prints _auxVarID membership, because that is the one property MANUAL's q
  // provably lacks (add_state does not call register_aux_def) -- so if the differential
  // sets differ, the two columns side by side say whether aux-ness is what drives it.

  return true;
}



inline bool
OCFESLV::_linear_cache_constants_match
( double const* cst )
const
{
  if( _vCst.empty() ) return true;
  if( !cst ) return false;
  if( _vCstVal.size() != _vCst.size() ) return false;
  for( size_t i = 0; i < _vCst.size(); ++i )
    if( cst[i] != _vCstVal[i] ) return false;
  return true;
}

inline bool
OCFESLV::_setup_linear_cache
( bool need_eqn, bool need_fct, double const* cst )
const
{
  if( options.SOLVE.LINEAR_EVAL_CACHE == Options::LINEAR_NONE ) return false;
  if( !const_cast<OCFESLV*>(this)->_ensure_deriv_cache() ) return false;
  if( !_linear_cache_constants_match( cst ) ) return false;

  bool const build_eqn = need_eqn && !_linEqnReady && _jacAllEqnLinear;
  bool const build_fct = need_fct && !_linFctReady && _jacAllFctLinear;
  if( !build_eqn && !build_fct ) return ( !need_eqn || _linEqnReady ) && ( !need_fct || _linFctReady );

  // Re-seed compressed FAD arrays at the zero state/input vector.  Constants
  // are kept at the setup values; the affine offset is then r(0).
  const_cast<OCFESLV*>(this)->_reset_autodiff( _FVar, _FInp, _FEqn, _FFct );
  for( size_t i = 0; i < _nCollVar; ++i ) _FVar[i].x() = 0.;
  for( size_t i = 0; i < _nCollInp; ++i ) _FInp[i].x() = 0.;
  if( !_vCst.empty() && !cst ) return false;
  for( size_t i = 0; i < _vCst.size(); ++i ) _OCFCst[i].update( &cst[i] );

  size_t const ncol = _nCollVar + _nCollInp;

  try{
    if( build_eqn ){
      _linEqn.assign( _nCollEqn, t_LinearRow() );
      std::vector< std::vector<double> > grow( _nCollEqn );
      std::vector<double*> gptr( _nCollEqn, nullptr );
      for( size_t k = 0; k < _nCollEqn; ++k ){
        grow[k].assign( _jacPatternEqn[k].size(), 0. );
        gptr[k] = grow[k].empty()? nullptr: grow[k].data();
      }
      _deriv_eqn( static_cast<double*>(nullptr), _nCollEqn? gptr.data(): nullptr, static_cast<size_t*>(nullptr), static_cast<size_t**>(nullptr),
                  _OCFEqn, _OCFVar, _OCFInp, _OCFDom, _OCFCst, _wrkOCF,
                  _FEqn.empty()? nullptr: _FEqn.data(),
                  _FVar.empty()? nullptr: _FVar.data(),
                  _FInp.empty()? nullptr: _FInp.data() );
      for( size_t k = 0; k < _nCollEqn; ++k ){
        _linEqn[k].offset = _FEqn[k].val();
        _linEqn[k].val    = std::move( grow[k] );
        if( _linEqn[k].val.size() != _jacPatternEqn[k].size() ) throw Exceptions( Exceptions::INTERNAL );
      }
      _linEqnReady = true;
    }

    if( build_fct ){
      _linFct.assign( _nCollFct, t_LinearRow() );
      std::vector< std::vector<double> > grow( _nCollFct );
      std::vector<double*> gptr( _nCollFct, nullptr );
      for( size_t k = 0; k < _nCollFct; ++k ){
        grow[k].assign( _jacPatternFct[k].size(), 0. );
        gptr[k] = grow[k].empty()? nullptr: grow[k].data();
      }
      _deriv_fct( static_cast<double*>(nullptr), _nCollFct? gptr.data(): nullptr, static_cast<size_t*>(nullptr), static_cast<size_t**>(nullptr),
                  _OCFEqn, _OCFVar, _OCFInp, _OCFDom, _OCFCst, _wrkOCF,
                  _FFct.empty()? nullptr: _FFct.data(),
                  _FVar.empty()? nullptr: _FVar.data(),
                  _FInp.empty()? nullptr: _FInp.data() );
      for( size_t k = 0; k < _nCollFct; ++k ){
        _linFct[k].offset = _FFct[k].val();
        _linFct[k].val    = std::move( grow[k] );
        if( _linFct[k].val.size() != _jacPatternFct[k].size() ) throw Exceptions( Exceptions::INTERNAL );
      }
      _linFctReady = true;
    }

    _linCstVal.clear();
    if( !_vCst.empty() ) _linCstVal.assign( cst, cst + _vCst.size() );
    (void)ncol;
    return true;
  }
  catch(...){
    if( build_eqn ){ _linEqn.clear(); _linEqnReady = false; }
    if( build_fct ){ _linFct.clear(); _linFctReady = false; }
    return false;
  }
}

inline bool
OCFESLV::_eval_linear_cache
( double* eqn, double* fct, double const* var, double const* inp, double const* cst )
const
{
  if( options.SOLVE.LINEAR_EVAL_CACHE == Options::LINEAR_NONE ) return false;
  if( !_linear_cache_constants_match( cst ) ) return false;
  bool const need_eqn = eqn != nullptr;
  bool const need_fct = fct != nullptr;
  if( _nCollVar && !var ) return false;
  if( _nCollInp && !inp ) return false;
  if( need_eqn && !_jacAllEqnLinear ) return false;
  if( need_fct && !_jacAllFctLinear ) return false;
  if( !_setup_linear_cache( need_eqn, need_fct, cst ) ) return false;

  auto eval_rows = [&]( double* out, std::vector<t_LinearRow> const& rows,
                        std::vector< std::vector<size_t> > const& pat )
  {
    for( size_t k = 0; k < rows.size(); ++k ){
      double r = rows[k].offset;
      for( size_t j = 0; j < pat[k].size(); ++j ){
        size_t const col = pat[k][j];
        double const z = col < _nCollVar ? var[col] : inp[col - _nCollVar];
        r += rows[k].val[j] * z;
      }
      out[k] = r;
    }
  };

  if( need_eqn ) eval_rows( eqn, _linEqn, _jacPatternEqn );
  if( need_fct ) eval_rows( fct, _linFct, _jacPatternFct );
  return true;
}



template <typename U>
inline void
OCFESLV::_reset_autodiff
( std::vector<FADType<U>>& Fvar, std::vector<FADType<U>>& Finp,
  std::vector<FADType<U>>& Feqn, std::vector<FADType<U>>& Ffct )
const
{
  // Seed compressed forward AD.  State columns occupy 0.._nCollVar-1 and
  // input columns occupy _nCollVar.._nCollVar+_nCollInp-1.  Columns sharing a
  // color never occur in the same equation/output row by construction of
  // _jacColor over the union of both row sets.
  // Re-create all FADBAD objects instead of resizing in place.  A cached
  // object may already own a derivative array with the old color count;
  // FADType::diff() assumes the size is unchanged and otherwise corrupts
  // memory when assertions are disabled.
  std::vector<FADType<U>>().swap( Fvar );
  Fvar.resize( _nCollVar );
  for( size_t i = 0; i < _nCollVar; ++i )
    if( _jacNColor ) Fvar[i].diff( _jacColor[i], _jacNColor );

  std::vector<FADType<U>>().swap( Finp );
  Finp.resize( _nCollInp );
  for( size_t i = 0; i < _nCollInp; ++i )
    if( _jacNColor ) Finp[i].diff( _jacColor[_nCollVar+i], _jacNColor );

  std::vector<FADType<U>>().swap( Feqn );
  Feqn.resize( _nCollEqn );
  std::vector<FADType<U>>().swap( Ffct );
  Ffct.resize( _nCollFct );
}


// ======================================================================
// OCFESLV::_audit_reduced_dof   (Stage 4, post-reduction structural DOF audit)
// Recovered 2026-07-02 from the 2026-06-30/07-01 transcripts (lost to a fork).
// Max-bipartite matching on the reduced _jacPatternEqn; skips interface-trace
// systems (nTrace>0) conservatively.  Fills _dofAudit.
// ======================================================================
inline void
OCFESLV::_audit_reduced_dof()
{
  _dofAudit.clear();
  if( _reductionPlan.empty() ) return;            // only audit reduced systems
  if( n_colloc_trace() != 0 ){                    // interface-tau: extension
    if( _solve_display_level() >= 1 )
      _disp( 2 ) << "OCFESLV::setup ** DOF audit: skipped: interface-trace system"
                   " (nTrace=" << n_colloc_trace() << ") is outside the"
                   " structural DOF audit (trace rows/cols not in the pattern).\n";
    return;
  }
  if( !_ensure_deriv_cache() ) return;            // need the sparsity pattern

  size_t const nR = _jacPatternEqn.size();        // equation rows
  size_t nC = 0;                                  // variable columns (max idx+1)
  for( auto const& row : _jacPatternEqn )
    for( size_t c : row ) if( c + 1 > nC ) nC = c + 1;

  // Maximum bipartite matching: equation rows -> variable columns (augmenting
  // path; same primitive used by _build_reduction_plan / the iplan matcher).
  std::vector<long> colMate( nC, -1 );            // column -> matched row, or -1
  std::function<bool(size_t,std::vector<char>&)> aug =
    [&]( size_t r, std::vector<char>& seen )->bool {
      for( size_t c : _jacPatternEqn[r] ){
        if( seen[c] ) continue;
        seen[c] = 1;
        if( colMate[c] < 0 || aug( (size_t)colMate[c], seen ) ){
          colMate[c] = (long)r;
          return true;
        }
      }
      return false;
    };
  size_t matched = 0;
  for( size_t r = 0; r < nR; ++r ){
    std::vector<char> seen( nC, 0 );
    if( aug( r, seen ) ) ++matched;
  }

  size_t const nfull = std::min( nR, nC );
  _dofAudit.ran        = true;
  _dofAudit.rows       = nR;
  _dofAudit.cols       = nC;
  _dofAudit.rank       = matched;
  _dofAudit.square     = ( nR == nC );
  _dofAudit.deficiency = ( matched < nfull ? nfull - matched : 0 );

  if( _dofAudit.deficiency ){
    _disp( 2 ) << "OCFESLV::setup ** DOF audit: post-reduction structural rank "
              << matched << " < " << nfull << " (rows=" << nR << " cols=" << nC
              << "): " << _dofAudit.deficiency << " redundant or missing"
                 " initial/boundary condition(s) for the reduced degrees of"
                 " freedom -- index reduction lowered the consistent IC/BC"
                 " count.  Check the IC/BC specification.\n";
  }
  else if( !_dofAudit.square ){
    std::cerr << "OCFESLV::setup ** DOF audit: post-reduction system not square"
                 " (rows=" << nR << " cols=" << nC << "): "
              << ( nR > nC ? "over" : "under" )
              << "-specified initial/boundary conditions for the reduced DOF.\n";
  }
  else if( _solve_display_level() >= 1 ){
    _disp( 2 ) << "OCFESLV::setup ** DOF audit: OK: reduced system square and"
                 " structurally full rank (" << nR << "), IC/BC count consistent"
                 " with the reduced degrees of freedom.\n";
  }
}


// ======================================================================
// rev133: solver-independent Jacobian assembly, numerical rank, and the
// setup-time determinacy audit.
//
// WHY.  Stage 0 of the determinacy workplan (2026-08-27) measured, on the
// rev132 corpus, that the existing structural DOF audit runs only on
// index-reduced IC_WEAK models, that the SVD spectrum audit is solve-time and
// silent on full-rank and rectangular systems, and that five drivers carry a
// primal null space the plan-level rank checks cannot see (PDE6_blk0 TRACE
// MANUAL, PDE16 IC_WEAK x3, PDE20e DEFICIENT; ALGFIELD2 and PDE40b are
// reference-point / tolerance artefacts).  The audit below runs at setup(),
// unconditionally on imposition mode, on the SAME assembled Jacobian solve()
// would see, with ONE rank tolerance for both ranks it compares.
//
// DISCRIMINATOR.  With M = [[A,B],[C,D]] (columns: states | tau; rows:
// physical | claims), the x=0 null vectors are exactly null([B;D]), so
//   k' = n_tau - rank([B;D]),   def(M) >= k' always,
//   def(M) == k'  <=>  the null space is pure tau (gauge only, primal determined),
//   def(M) >  k'  <=>  some null vector has state content (primal undetermined).
// Under IC_WEAK n_tau = 0 and k' = 0.  def(M) < k' is impossible and is
// reported as a tolerance inconsistency between the two rank calls.
//
// REV-A DELIVERY RULE.  assemble() is a SECOND path here; solve() keeps its
// inline assembly untouched.  CRONOS_AUDIT_ASSEMBLE_XCHECK=1 makes solve()
// compare the two at iterate 0 and print the max difference (expected 0
// exactly: same deriv() calls, same scaling code).  solve()/eval() switch to
// assemble() in a later revision, after a corpus sweep shows 0 everywhere.
// ======================================================================

inline bool
OCFESLV::_assemble_pattern
( t_JacAssembly& A )
{
  if( !_issetup ) return false;
  size_t const n = n_colloc_eqn();
  size_t const m = n_colloc_sta();
  if( A.pattern_ready && A.n == n && A.m == m ) return true;
  // Fixed sparsity pattern (unchanged across iterations): per-row nnz then column
  // indices.  EQUATION ROWS ONLY -- the all-row deriv(nnz,colnz) overload stacks
  // _nCollEqn + _nCollFct rows and would overflow these n=_nCollEqn buffers by
  // _nCollFct whenever the model declares scalar outputs (heap corruption).
  A.pattern_ready = false;
  A.rownnz.assign( n, 0 );
  if( !deriv( nullptr, nullptr, A.rownnz.data(), nullptr,
              nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr ) ) return false;
  A.rowcol.assign( n, std::vector<size_t>() );
  std::vector<size_t*> colptr( n, nullptr );
  for( size_t i=0; i<n; ++i ){
    A.rowcol[i].assign( A.rownnz[i], 0 );
    if( A.rownnz[i] ) colptr[i] = A.rowcol[i].data();
  }
  if( !deriv( nullptr, nullptr, A.rownnz.data(), colptr.data(),
              nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr ) ) return false;
  A.nnz = 0; for( size_t v : A.rownnz ) A.nnz += v;
  A.grad.assign( A.nnz, 0. );
  A.gradptr.assign( n, nullptr );
  { size_t off=0; for( size_t i=0; i<n; ++i ){ if( A.rownnz[i] ){ A.gradptr[i] = A.grad.data()+off; off += A.rownnz[i]; } } }
  A.n = n; A.m = m; A.pattern_ready = true;
  A.ic_sat_ready = false;
  return true;
}

inline void
OCFESLV::_ensure_ic_sat
( t_JacAssembly& A, bool ic_override )
{
  // IC_WEAK re-injects the interface SAT penalty onto the pinned override rows: the
  // continuity pin is agnostic to the IC form, but under IC_WEAK the INITIAL interface
  // rows also carry a SAT penalty that couples the element-face duplicates, and a bare
  // pin drops it -> a spurious converged-but-wrong interface mode (blow-up at
  // multi-element corners).  The solver's own frozen weak-SAT stencil (compiled per
  // receiver row) is added back on top of the pin.  Cached per (window, override).
  if( A.ic_sat_ready && A.ic_sat_window == _marchWindow && A.ic_sat_override == ic_override ) return;
  A.ic_sat.clear();
  size_t const n = A.n;
  if( ic_override && options.INTERFACE.IMPOSITION == Options::IC_WEAK ){
    std::set<size_t> ovr_rows;
    for( auto const& [st, rows] : _marchICRows ) for( auto const& rr : rows ) if( rr.row < n ) ovr_rows.insert( rr.row );
    for( auto const& term : _plan.weak_sat_terms() ){
      if( term.sat_coeffs.empty() ) continue;                        // prefactor == 0
      size_t const row = term.receiver_emit_flat;                    // global emitted row (matches rr.row)
      if( row >= n || !ovr_rows.count( row ) ) continue;
      auto& acc = A.ic_sat[row];
      for( auto const& [col, c] : term.sat_coeffs ) if( col < n ) acc.emplace_back( col, c );
    }
  }
  A.ic_sat_ready = true; A.ic_sat_window = _marchWindow; A.ic_sat_override = ic_override;
}

inline long
OCFESLV::_ic_col_for_row
( t_JacAssembly const& A, FFVar const& st, t_MarchICRow const& rr )
const
{
  // Row -> pin column, read from the row's OWN Jacobian nonzeros (the unique nonzero
  // that is an LB-face DOF of st, via _marchICColToNode) rather than a precomputed
  // index: robust to any mismatch between the equation's row-emission order and st's
  // var_dom/DOF ordering, and the column is always in the row's sparsity.
  auto itc = _marchICColToNode.find( st );
  if( itc != _marchICColToNode.end() && rr.row < A.n )
    for( size_t k = 0; k < A.rownnz[rr.row]; ++k )
      if( itc->second.count( A.rowcol[rr.row][k] ) ) return (long)A.rowcol[rr.row][k];
  return ( rr.lb_col < A.n ) ? (long)rr.lb_col : -1L;   // structural fallback (same-order / time-only)
}

inline bool
OCFESLV::_residual_core
( double const* var, double const* inp, double const* cst, double* r, t_JacAssembly& A, bool ic_override )
{
  // Residual at var, with the marching general-IC continuity override (window k>0):
  // each transfer state's INITIAL residual rows become  var[col] - terminal  (+ the
  // re-injected SAT under IC_WEAK).  The single residual path solve() and the audit
  // share -- the former apply_ic_res lambda.
  if( !_assemble_pattern( A ) ) return false;
  std::vector<double> inpauto;
  if( !inp && n_colloc_inp() ){ inpauto.assign( n_colloc_inp(), 0. ); inp = inpauto.data(); }
  if( !eval( r, nullptr, var, inp, cst ) ) return false;
  bool const ovr = ic_override && !_marchICRows.empty();
  if( !ovr ) return true;
  _ensure_ic_sat( A, ic_override );
  for( auto const& [st, rows] : _marchICRows ){
    auto itt = _marchICTerminal.find( st );
    auto itc = _marchICColToNode.find( st );
    for( auto const& rr : rows ){
      if( rr.row >= A.n ) continue;
      long const col = _ic_col_for_row( A, st, rr );
      if( col < 0 ) continue;
      double term = 0.0;                                   // ic_term_for_col
      if( itt != _marchICTerminal.end() && !itt->second.empty() ){
        term = itt->second.front();
        if( itc != _marchICColToNode.end() ){
          auto itk = itc->second.find( (size_t)col );
          if( itk != itc->second.end() && itk->second < itt->second.size() ) term = itt->second[ itk->second ];
        }
      }
      // rev166: the pin is a TRACE at the window lower bound, not a nodal identity.
      // CRONOS_MARCH_TRACE_PIN=1 enables it; default 0 keeps rev165 exactly.  On LGL/LGR the
      // weights are delta_{j0}, so the two forms coincide there.
      double rv = var[ (size_t)col ] - term;               // form-agnostic continuity pin
      {
        // rev169 D4: DEFAULT ON.  The nodal pin is wrong on any collocation family without a
        // node at the window lower bound (LG), and silently so -- machine-zero per-window
        // residuals with a 6.1e-3 solution error (OCFE_PDE25, 2026-09-08).  MEASURED inert
        // elsewhere: corpus sweep 20260908_211746 with the pin on is 1199/1199 determinacy
        // lines identical, PASS/FAIL unchanged, including the two LGR-evolution drivers
        // (PDE5b, PDE20d) whose override rows are multi-column.  CRONOS_MARCH_TRACE_PIN=0
        // restores the rev165 nodal pin for bisection.
        static int const kTracePin = 1;   // CRONOS_MARCH_TRACE_PIN (retired 2026-10-07, WORKPLAN 3.B batch 2b)
        size_t stride_evo = 0; std::vector<double> wtr;
        if( kTracePin && _march_ic_lb_trace( st, stride_evo, wtr ) && wtr.size() > 1 ){
          double sw = 0., vv = 0.;
          for( size_t j = 0; j < wtr.size(); ++j ){
            size_t const cj = (size_t)col + j * stride_evo;
            if( cj >= A.m ) { sw = 0.; break; }
            sw += wtr[j];
            vv += wtr[j] * var[ cj ];
          }
          if( sw != 0. ) rv = vv - sw * term;
        }
      }
      auto is = A.ic_sat.find( rr.row );                   // + re-injected interface SAT (IC_WEAK)
      if( is != A.ic_sat.end() ) for( auto const& [sc, c] : is->second ) rv += c * var[ sc ];
      r[ rr.row ] = rv;
    }
  }
  return true;
}

inline bool
OCFESLV::_assemble_core
( double const* var, double const* inp, double const* cst, t_JacAssembly& A,
  bool equilibrate, bool ic_override, double const* r_known )
{
  if( !_issetup ) return false;
  size_t const n = n_colloc_eqn();
  size_t const m = n_colloc_sta();
  if( n && !var ) return false;
  auto const tA = std::chrono::steady_clock::now();
  auto _ms = []( std::chrono::steady_clock::time_point a, std::chrono::steady_clock::time_point b )
             { return std::chrono::duration<double,std::milli>( b - a ).count(); };

  if( !_assemble_pattern( A ) ) return false;
  A.n_tau  = _nTraceVar;
  A.tau0   = ( _nTraceVar ? _traceVarOffset : m );
  A.window = _marchWindow;

  std::vector<double> inpauto;
  if( !inp && n_colloc_inp() ){ inpauto.assign( n_colloc_inp(), 0. ); inp = inpauto.data(); }

  // Residual: the caller's (solve() carries the accepted line-search residual, which
  // IS eval(var) -- no second evaluation) or a fresh one.
  std::vector<double> rbuf;
  double const* r = r_known;
  if( !r ){
    rbuf.assign( n, 0. );
    if( !_residual_core( var, inp, cst, rbuf.data(), A, ic_override ) ) return false;
    r = rbuf.data();
  }

  // Jacobian values at var, then the override rows (unit row on the pin column, plus
  // the SAT stencil under IC_WEAK) -- the former apply_ic_jac lambda.
  std::fill( A.grad.begin(), A.grad.end(), 0. );
  if( !deriv( nullptr, A.gradptr.data(), nullptr, nullptr,
              nullptr, nullptr, nullptr, nullptr, var, inp, cst ) ) return false;
  A.ic_override = ( ic_override && !_marchICRows.empty() );
  if( A.ic_override ){
    _ensure_ic_sat( A, ic_override );
    for( auto const& [st, rows] : _marchICRows )
      for( auto const& rr : rows ){
        if( rr.row >= n || !A.gradptr[rr.row] ) continue;
        long const col = _ic_col_for_row( A, st, rr );
        static int const kTracePinJ = 1;   // CRONOS_MARCH_TRACE_PIN (retired 2026-10-07, WORKPLAN 3.B batch 2b);   // rev169 D4: default ON, see the residual site
        size_t strideJ = 0; std::vector<double> wtrJ;
        if( kTracePinJ ){
          if( !_march_ic_lb_trace( st, strideJ, wtrJ ) || wtrJ.size() <= 1 ){ wtrJ.clear(); strideJ = 0; }
        }
        auto is = A.ic_sat.find( rr.row );
        for( size_t k = 0; k < A.rownnz[rr.row]; ++k ){
          size_t const c = A.rowcol[rr.row][k];
          double val = ( (long)c == col ) ? 1.0 : 0.0;                  // pin
          if( kTracePinJ && !wtrJ.empty() && strideJ ){                  // rev166: trace pin
            val = 0.;
            for( size_t j = 0; j < wtrJ.size(); ++j )
              if( (long)c == col + (long)( j * strideJ ) ){ val = wtrJ[j]; break; }
          }
          if( is != A.ic_sat.end() )                                    // + SAT stencil (accumulate on col)
            for( auto const& [sc, cf] : is->second ) if( sc == c ) val += cf;
          A.gradptr[rr.row][k] = val;
        }
      }
  }
  A.ms_deriv = _ms( tA, std::chrono::steady_clock::now() );
  auto const tB = std::chrono::steady_clock::now();

  // Equilibration: row inf-norms over state columns -> row weights; column inf-norms
  // of the row-scaled matrix -> column weights.  All ones when not requested.
  A.rw.assign( n, 1.0 );
  A.cw.assign( m, 1.0 );
  A.equilibrated = equilibrate;
  if( equilibrate ){
    { size_t off=0;
      for( size_t i=0; i<n; ++i ){
        double rmax=0.0;
        for( size_t k=0; k<A.rownnz[i]; ++k )
          if( A.rowcol[i][k] < m ) rmax = std::max( rmax, std::fabs( A.grad[off+k] ) );
        A.rw[i] = 1.0/std::max(1.0,rmax);
        off += A.rownnz[i];
      } }
    std::vector<double> cmax( m, 0.0 );
    { size_t off=0;
      for( size_t i=0; i<n; ++i ){
        for( size_t k=0; k<A.rownnz[i]; ++k ){
          size_t const c = A.rowcol[i][k];
          if( c < m ) cmax[c] = std::max( cmax[c], std::fabs( A.grad[off+k] )*A.rw[i] );
        }
        off += A.rownnz[i];
      } }
    for( size_t j=0; j<m; ++j ) A.cw[j] = 1.0/std::max(1.0,cmax[j]);
  }

  // Sparse assembly J = diag(rw) * J * diag(cw)  (n x m, state+tau columns only).  The
  // mask (rev137: empty; rev138: the drop) removes rows/columns at this point.
  bool const masked = !A.mask_rows.empty() || !A.mask_cols.empty();
  std::vector<char> mrow, mcol;
  if( masked ){
    mrow.assign( n, 0 ); mcol.assign( m, 0 );
    for( size_t i : A.mask_rows ) if( i < n ) mrow[i] = 1;
    for( size_t j : A.mask_cols ) if( j < m ) mcol[j] = 1;
  }
  arma::umat loc( 2, (arma::uword)std::max<size_t>( A.nnz, 1 ) );
  arma::vec  val( (arma::uword)std::max<size_t>( A.nnz, 1 ) );
  arma::uword e = 0;
  { size_t off=0;
    for( size_t i=0; i<n; ++i ){
      for( size_t k=0; k<A.rownnz[i]; ++k ){
        size_t const c = A.rowcol[i][k];
        if( c < m && !( masked && ( mrow[i] || mcol[c] ) ) ){
          loc(0,e)=(arma::uword)i; loc(1,e)=(arma::uword)c; val(e)=A.grad[off+k]*A.rw[i]*A.cw[c]; ++e; }
      }
      off += A.rownnz[i];
    } }
  A.J = arma::sp_mat( (arma::uword)n, (arma::uword)m );
  if( e > 0 ){
    arma::umat L = loc.cols( 0, e-1 );
    arma::vec  V = val.subvec( 0, e-1 );
    A.J = arma::sp_mat( L, V, (arma::uword)n, (arma::uword)m );
  }
  A.r.set_size( (arma::uword)n );
  for( size_t i=0; i<n; ++i ) A.r(i) = ( masked && mrow[i] ) ? 0. : r[i]*A.rw[i];
  A.ms_asm = _ms( tB, std::chrono::steady_clock::now() );
  return true;
}

inline bool
OCFESLV::assemble
( double const* var, double const* inp, double const* cst, t_JacAssembly& A, bool equilibrate )
{
  FixedInputScope fixed_scope_( _fixedInputDepth );  std::vector<double> fixed_buf_;  inp = _with_fixed_inputs( inp, fixed_buf_ );   // fix_input: model-supplied values, outermost call only
  // Uses the environment's CURRENT marching state: inside a march window the
  // override rows are live exactly as solve() would see them; outside, not.
  bool const ovr = _inMarch && _marchICActive && !_marchICRows.empty();
  return _assemble_core( var, inp, cst, A, equilibrate, ovr, nullptr );
}

inline bool
OCFESLV::assemble_window
( size_t w, double const* var, double const* inp, double const* cst, t_JacAssembly& A, bool equilibrate )
{
  FixedInputScope fixed_scope_( _fixedInputDepth );  std::vector<double> fixed_buf_;  inp = _with_fixed_inputs( inp, fixed_buf_ );   // fix_input: model-supplied values, outermost call only
  if( !_issetup ) return false;
  if( _marchGrid.size() < 2 || w >= n_march_steps() ) return false;

  // Save the per-window state solve_marching() mutates, set it for window w,
  // restore on every exit path.
  bool   const sIn  = _inMarch;
  size_t const sWin = _marchWindow;
  bool   const sAct = _marchICActive;
  auto   const sTerm = _marchICTerminal;
  auto restore = [&](){ _inMarch = sIn; _marchWindow = sWin; _marchICActive = sAct; _marchICTerminal = sTerm; };

  _inMarch = true; _marchWindow = w; _marchICActive = false;
  if( w > 0 ){
    if( !_marchICRowsReady ){
      // _build_march_ic_rows() needs row provenance, which one eval builds.
      std::vector<double> inpauto;
      double const* ip = inp;
      if( !ip && n_colloc_inp() ){ inpauto.assign( n_colloc_inp(), 0. ); ip = inpauto.data(); }
      std::vector<double> rtmp( n_colloc_eqn(), 0. );
      if( !eval( rtmp.data(), nullptr, var, ip, cst ) ){ restore(); return false; }
      _build_march_ic_rows();
    }
    _marchICActive = !_marchICRows.empty();
    _marchICTerminal.clear();      // no terminal profile at setup: pin against var itself
  }
  bool const ok = _assemble_core( var, inp, cst, A, equilibrate, _marchICActive, nullptr );
  A.window = w;
  restore();
  return ok;
}


inline bool
OCFESLV::_wdecide_restore
()
{
  if( !_knob_DROP_WDECIDE() ) return false;
  if( _wdecideRetry ) return false;                       // inside the restore recursion
  if( options.INTERFACE.DROP_POLICY == Options::DROP_OFF ) return false;
  if( _strong_dropped_claims.empty() ) return false;
  if( options.INTERFACE.IMPOSITION != Options::IC_TRACE
   && options.INTERFACE.IMPOSITION != Options::IC_STRONG ) return false;
  // (2026-10-07, WORKPLAN 3.D: CRONOS_AUDIT_MAX_N no longer skips this W-test -- over its budget the drop used to
  // "stand unverified", a DECISION changed by an environment variable named as an audit budget.  The W-test always runs
  // (the default all along: the variable defaults to unlimited; above the dense cap the QR basis keeps it affordable,
  // rev140).  The variable remains the determinacy AUDIT's budget, its other use.)

  // rev142a: n_colloc_eqn() was read into an unused `n` here, residue of the rank pass
  // that moved out in rev140.  m stays -- the W build tests c0 < m.
  size_t const m = n_colloc_sta();
  std::ostringstream rep;
  rep << "[w-decide] dropped=" << _strong_dropped_claims.size();
  if( _knob_CLAIM_FROM_C0() ){
    std::map<size_t,size_t> by_state;
    for( auto const& k : _strong_dropped_claims ) ++by_state[ std::get<1>( k ) ];
    rep << "  by state:";
    for( auto const& [sid, n] : by_state ){
      FFVar const* sv = _state_by_id( sid );
      rep << " "; if( sv ) rep << *sv; else rep << "#" << sid; rep << "=" << n;
    }
  }

  // ---- per-window BEFORE deficiencies (win 0, and win 1 when marching) -------
  auto const mctx = _march_context();
  bool const marching = mctx.active;
  size_t const nwin = marching ? std::min<size_t>( 2, mctx.n_steps ) : 1;
  std::vector<double> x, u; std::string reftag;
  _reference_point( x, u, reftag );
  double const* xv = x.empty() ? nullptr : x.data();
  double const* uv = u.empty() ? nullptr : u.data();
  double const* cv = _vCstVal.empty() ? nullptr : _vCstVal.data();
  std::vector<size_t> def_before( nwin, 0 );
  t_JacAssembly A0;
  t_RankResult  RR0;
  bool have0 = false;
  for( size_t w = 0; w < nwin; ++w ){
    t_JacAssembly A;
    bool ok = ( w == 0 && !marching ) ? assemble( xv, uv, cv, A, true )
                                      : assemble_window( w, xv, uv, cv, A, true );
    if( !ok ){ rep << " win" << w << "=ASSEMBLE-FAILED (no decision)"; _disp( 3 ) << rep.str() << std::endl; return false; }
    t_RankResult const RR = numeric_rank( A.J, options.DETERMINACY.TOL, 0.0 );
    def_before[w] = ( (size_t)A.J.n_cols > RR.rank ? (size_t)A.J.n_cols - RR.rank : 0 );
    if( w == 0 ){ A0 = std::move( A ); RR0 = RR; have0 = true; }
  }
  rep << " def0=" << def_before[0];
  if( nwin > 1 ) rep << " def1=" << def_before[1];
  if( !def_before[0] ){
    rep << "  ALL-IMPLIED (full rank with the drop)"; _disp( 3 ) << rep.str() << std::endl; return false;
  }

  // ---- W over the dropped claims on the win-0 null space ---------------------
  // rev140: dense SVD stays authoritative BELOW the cap (bit-compatible with every
  // rev139 decision); ABOVE the cap the QR dead-column basis takes over.
  static bool const kForceQR = false;   // CRONOS_WDECIDE_FORCE_QR (retired 2026-10-07, WORKPLAN 3.B batch 2b)
  bool const over_cap = std::max<size_t>( A0.J.n_rows, A0.J.n_cols ) > _knob_DETERMINACY_MAX_DENSE();
  arma::mat NB;              // n_cols x k null basis, original column coordinates
  if( over_cap || kForceQR ){
    size_t rk2 = 0;
    if( !have0 || !_null_basis( A0.J, RR0.tol_abs, rk2, NB ) ){
      rep << "  QR-NULL-FAILED (no decision)"; _disp( 3 ) << rep.str() << std::endl;
      // 4R.1: never silent -- a dropped claim may be needed and cannot be verified
      _disp( 1 ) << "OCFESLV::setup ** WARNING: the W-test of dropped interface claims found no null basis above its dense"
                 << " cap (" << A0.J.n_rows << "x" << A0.J.n_cols << "): dropped claims stand UNVERIFIED -- a build with"
                 << " SPQR or Eigen decides them" << std::endl;
      return false;
    }
    rep << ( kForceQR && !over_cap ? " [qr-null forced]" : " [qr-null]" );
  }
  else{
    // rev140: right-vectors-only economy SVD -- U was computed and discarded, at
    // ~3x the cost on the 1500-2300-square drop systems (the sweep-slowdown fix).
    arma::mat U2, V2; arma::vec sv;
    if( !have0 || !arma::svd_econ( U2, sv, V2, arma::mat( A0.J ), "right" ) ){
      rep << "  SVD-FAILED (no decision)"; _disp( 3 ) << rep.str() << std::endl; return false;
    }
    std::vector<arma::uword> nullk;
    for( arma::uword k2 = 0; k2 < sv.n_elem; ++k2 ) if( sv(k2) <= RR0.tol_abs ) nullk.push_back( k2 );
    for( arma::uword k2 = (arma::uword)sv.n_elem; k2 < (arma::uword)A0.J.n_cols; ++k2 ) nullk.push_back( k2 );
    NB.zeros( V2.n_rows, (arma::uword)std::max<size_t>( nullk.size(), 1 ) );
    for( size_t kk = 0; kk < nullk.size(); ++kk ) NB.col( (arma::uword)kk ) = V2.col( nullk[kk] );
    if( nullk.empty() ) NB.reset(), NB.zeros( V2.n_rows, 0 );
  }
  std::vector<t_FrozenInterfaceClaimKey> tcand;
  std::vector<std::vector<std::pair<size_t,double>>> tcoef;
  { std::vector<std::pair<size_t,double>> buf;
    for( auto const& k : _strong_dropped_claims ){
      t_InterfaceClaimKey const ck = _unfreeze_interface_claim_key( k );
      buf.clear();
      if( _trace_claim_coeffs( ck, buf ) ){ tcand.push_back( k ); tcoef.push_back( buf ); }
    } }
  size_t const knull = (size_t)NB.n_cols;
  arma::mat W( (arma::uword)tcand.size(), (arma::uword)std::max<size_t>( knull, 1 ),
               arma::fill::zeros );
  for( size_t i = 0; i < tcand.size(); ++i ){
    double scale = 0.0;
    for( auto const& [c0, cf] : tcoef[i] ) if( c0 < m ) scale += std::fabs( cf );
    for( size_t kk = 0; kk < knull; ++kk ){
      double r = 0.0;
      for( auto const& [c0, cf] : tcoef[i] ){
        if( c0 >= m ) continue;
        r += cf * NB( (arma::uword)c0, (arma::uword)kk ) * A0.cw[c0];
      }
      W( i, kk ) = r / std::max( scale, 64. * DBL_EPSILON );
    }
  }
  std::set<t_FrozenInterfaceClaimKey> needed;
  { arma::mat Q3, R3; arma::uvec P3;
    if( W.n_cols && arma::qr( Q3, R3, P3, W.t(), "vector" ) ){
      double const r0 = ( R3.n_rows && R3.n_cols ) ? std::fabs( R3(0,0) ) : 0.0;
      for( arma::uword i = 0; i < std::min( R3.n_rows, R3.n_cols ); ++i )
        if( std::fabs( R3(i,i) ) > 1e-8 * std::max( r0, 1e-300 ) )
          needed.insert( tcand[ (size_t)P3(i) ] );
    } }
  rep << " null=" << knull << " tested=" << tcand.size() << " needed=" << needed.size();
  if( needed.empty() ){
    rep << "  ALL-IMPLIED (deficiency is not the dropped rows')";
    _disp( 3 ) << rep.str() << std::endl; return false;
  }

  // ---- restore: protect the needed subset, re-derive once, verify, else revert
  rep << " -> KEEP-EXPLICIT (restore pass)";
  _disp( 3 ) << rep.str() << std::endl;
  std::set<t_FrozenInterfaceClaimKey> added;
  for( auto const& ck : needed )
    if( _strong_protected_claims.insert( ck ).second ) added.insert( ck );
  _wdecideRetry = true;
  _auditQuiet   = true;             // rev140: the restore attempt audits quietly;
  bool ok = setup();                // the kept system's lines are replayed below,
  _auditQuiet   = false;            // a revert prints via the second recursion.
  _wdecideRetry = false;
  // 2026-10-10: the RESTORED system is measured exactly as the dropped one was above -- its windows assembled at the
  // reference and numeric_rank'ed by this W-test itself.  It used to be read from _determinacy, which the numerical
  // audit fills only when DETERMINACY.AUDIT is on (the structural count otherwise): the verdict, hence the plan the
  // solver got, depended on a DIAGNOSTIC flag (OCFE_MMPDE1 with an adaptive monitor, found 9 Oct).
  size_t worst_after = 0, worst_before = 0;
  for( size_t w = 0; w < nwin; ++w ) worst_before = std::max( worst_before, def_before[w] );
  if( ok ){
    std::vector<double> xa, ua; std::string reftag_a;
    _reference_point( xa, ua, reftag_a );
    double const* xav = xa.empty() ? nullptr : xa.data();
    double const* uav = ua.empty() ? nullptr : ua.data();
    double const* cav = _vCstVal.empty() ? nullptr : _vCstVal.data();
    bool const marching_a = mctx.active;
    size_t const nwin_a = marching_a ? std::min<size_t>( 2, mctx.n_steps ) : 1;
    for( size_t w = 0; w < nwin_a; ++w ){
      t_JacAssembly A;
      bool const oka = ( w == 0 && !marching_a ) ? assemble( xav, uav, cav, A, true )
                                                 : assemble_window( w, xav, uav, cav, A, true );
      if( !oka ){ ok = false; break; }
      t_RankResult const RRa = numeric_rank( A.J, options.DETERMINACY.TOL, 0.0 );
      size_t const def_after = ( (size_t)A.J.n_cols > RRa.rank ? (size_t)A.J.n_cols - RRa.rank : 0 );
      worst_after = std::max( worst_after, def_after );
    }
  }
  if( !ok || worst_after >= worst_before ){
    _disp( 3 ) << "[w-decide] restore did not strictly improve (worst excess "
              << worst_after << " >= " << worst_before
              << ( ok ? "" : "; re-derive FAILED" )
              << "): REVERTING to the dropped plan." << std::endl;
    for( auto const& ck : added ) _strong_protected_claims.erase( ck );
    _wdecideRetry = true;
    ok = setup();
    _wdecideRetry = false;
    _disp( 3 ) << "[w-decide] reverted; the pre-restore verdict stands." << std::endl;
    return ok;
  }
  _disp( 3 ) << "[w-decide] restore VERIFIED: worst excess " << worst_before
            << " -> " << worst_after << "; " << added.size()
            << " claim(s) keep-explicit." << std::endl;
  for( size_t w2 = 0; w2 < _determinacy.size(); ++w2 )
    _print_determinacy_window( _determinacy[w2], w2 );
  return true;
}

inline bool
OCFESLV::_reference_point
( std::vector<double>& x, std::vector<double>& u, std::string& tag )
{
  // Deterministic nonzero reference point (agreed 2026-08-26): the stored
  // reference profile (update_ref / init()) where one exists, with every exact
  // zero replaced -- STATES by the deterministic NON-UNIFORM value
  // 1 + 0.5 sin(0.7 (i+1)), INPUTS by the uniform 1.0.  No randomness anywhere.
  //
  // Three measured reasons for the rule (Stage 0-2, 2026-08-27):
  //  * ALGFIELD2: init() seeds a=p=0, dLINK/dp=a, dLINK/da=p -> 1280 empty rows and
  //    columns at the seed that the formulation does not have.  Any nonzero fill
  //    removes it.
  //  * MESHONLY10 (rev133 sweep): a UNIFORM fill of 1.0 on mesh-coordinate states put
  //    every node at the same position -- a degenerate mesh -- and the audit reported
  //    def=8 at gap 9 that the solve-time oracle never saw.  Hence non-uniform.
  //  * AE2 (rev133 sweep): the six registered inputs (q0, V1, V2, ...) are zero at
  //    setup, every material balance collapses to 0*(...), and the 10-state system
  //    read rank 1.  Hence the inputs are filled too.
  size_t const m = n_colloc_sta(), ni = n_colloc_inp();
  x.assign( m, 0. ); u.assign( ni, 0. );
  bool ok = false;
  std::string ref_err;                       // rev143: report-only failure capture
  try{
    ok = init( m ? x.data() : nullptr, ni ? u.data() : nullptr,
               _vCstVal.empty() ? nullptr : _vCstVal.data() );
    if( !ok ) ref_err = _issetup ? "init() returned false"
                                 : "init() returned false (setup incomplete: _issetup=false)"; // rev144
  }
  catch( std::exception const& e ){ ok = false; ref_err = std::string( "init() threw: " ) + e.what(); }
  catch( ... ){ ok = false; ref_err = "init() threw a non-std exception"; }
  if( !ok ){ x.assign( m, 0. ); u.assign( ni, 0. ); }
  tag = ok ? "init" : "zeros";
  // rev143: locate WHERE the restore re-derive loses the init reference.  MEASURED:
  // primary audits run at ref=init while _wdecide_restore's nested setup() re-derives
  // at ref=zeros (PDE6/blk0/MESHONLY10 under AUDIT_SETUP_ONLY), and the zeros assembly
  // then fails on EOS terms.  This print answers WHY init() fails there -- returned
  // false, or threw what.  Report-only, default off: CRONOS_REF_AUDIT=1 (also under
  // CRONOS_DEBUG_RANK).
  { static bool const kRefAudit = Options::_env_flag( "CRONOS_REF_AUDIT", false )
                               || Options::_env_flag( "CRONOS_DEBUG_RANK", false );
    if( kRefAudit )
      _disp( 3 ) << "  [ref-audit] base=" << ( ok ? "init" : "zeros" )
                << " m=" << m << " ni=" << ni
                << " wdecide_retry=" << ( _wdecideRetry ? 1 : 0 )
                << " audit_quiet="   << ( _auditQuiet   ? 1 : 0 )
                << ( ok ? std::string() : "  ** " + ref_err + " **" ) << std::endl; }
  //  * PDE23 / PDE23b (rev134 sweep): inputs declared with an InputContinuity level
  //    are asserted by eval() against the supplied data; a non-uniform input fill is
  //    discontinuous by construction and eval() refuses it (14 SKIPPED lines).  Inputs
  //    only ever needed to be NONZERO (AE2), so they get the uniform 1.0; the
  //    non-uniform fill is for coordinate-like STATES only (MESHONLY10).
  size_t nzx = 0, nzu = 0;
  for( size_t i = 0; i < x.size(); ++i )
    if( x[i] == 0. || !std::isfinite( x[i] ) ){ x[i] = 1. + 0.5 * std::sin( 0.7 * (double)( i + 1 ) ); ++nzx; }
  for( size_t i = 0; i < u.size(); ++i )
    if( u[i] == 0. || !std::isfinite( u[i] ) ){ u[i] = 1.; ++nzu; }
  if( nzx ) tag += "+fillsin(x:" + std::to_string( nzx ) + ")";
  if( nzu ) tag += "+fill1(u:" + std::to_string( nzu ) + ")";
  return true;
}


inline bool
OCFESLV::audit_determinacy
( double const* var, double const* inp )
{
  if( !_issetup ) return false;
  _audit_determinacy( var, inp );
  return !_determinacy.empty() && _determinacy.front().ran;
}


inline void
OCFESLV::_audit_determinacy
( double const* var_user, double const* inp_user )
{
  _determinacy.clear();
  if( !var_user && !options.DETERMINACY.AUDIT ) return;

  // Point: the caller's own (a driver's physical initial guess -- MESHONLY-class
  // models want a monotone mesh, which no generic fill can supply) or the
  // deterministic reference point.
  std::vector<double> x, u; std::string reftag;
  if( var_user ){
    x.assign( var_user, var_user + n_colloc_sta() );
    u.assign( n_colloc_inp(), 0. );
    if( inp_user && n_colloc_inp() ) u.assign( inp_user, inp_user + n_colloc_inp() );
    reftag = "user";
  }
  else _reference_point( x, u, reftag );
  double const* cst = _vCstVal.empty() ? nullptr : _vCstVal.data();
  auto const mctx = _march_context();
  bool const marching = mctx.active;
  size_t const nwin = marching ? std::min<size_t>( 2, mctx.n_steps ) : 1;

  auto imp_name = [&]() -> char const* {
    return options.INTERFACE.IMPOSITION == Options::IC_WEAK   ? "IC_WEAK"
         : options.INTERFACE.IMPOSITION == Options::IC_TRACE  ? "IC_TRACE"
         : options.INTERFACE.IMPOSITION == Options::IC_STRONG ? "IC_STRONG" : "?"; };
  // rev142a: A_n_ne_m and verdict_name were duplicated here and in
  // _print_determinacy_window.  rev140 moved the printing out and left these behind,
  // where they compiled as dead code and emitted -Wunused-but-set-variable.  Deleted
  // rather than silenced: two copies of a verdict-naming switch is exactly how a
  // report and the record it reports on drift apart.  The live copies are the
  // printer's.


  for( size_t w = 0; w < nwin; ++w ){
    t_Determinacy D;
    D.window = w; D.marching = marching; D.ref_point = reftag; D.mode = imp_name();
    // Equilibrated exactly as solve() equilibrates (row inf-norms, then column
    // inf-norms): the rank threshold is a fraction of the largest column norm, so
    // on the RAW Jacobian a legitimately small-scaled column (a tau column against
    // a 1/h^2 stencil column) would read as dead.  On the scaled matrix column
    // norms are O(1) and the threshold means what it says; it is also the matrix
    // the Stage 0 SVD oracle audited.
    bool const eq = options.SOLVE.EQUILIBRATE;
    t_JacAssembly A;
    // rev141: routine-sweep budget -- above CRONOS_AUDIT_MAX_N the window is reported
    // SKIPPED reason=budget(...) instead of factorized (default 0 = unlimited).
    static size_t const kMaxAuditN = (size_t)Options::_env_dbl( "CRONOS_AUDIT_MAX_N", 0. );
    bool const over_budget = kMaxAuditN && std::max( n_colloc_eqn(), n_colloc_sta() ) > kMaxAuditN;
    bool const ok = !over_budget && ( marching ? assemble_window( w, x.data(), u.empty() ? nullptr : u.data(), cst, A, eq )
                             : assemble( x.data(), u.empty() ? nullptr : u.data(), cst, A, eq ) );
    if( !ok ){
      if( over_budget ){
        D.n = n_colloc_eqn(); D.m = n_colloc_sta();
        D.skipped_reason = "budget(n>" + std::to_string( kMaxAuditN ) + ")";
      }
      else D.skipped_reason = _issetup ? "assembly failed"
                                       : "assembly failed (setup incomplete)"; // rev144
    }
    else{
      D.ran = true; D.n = A.n; D.m = A.m; D.n_tau = A.n_tau; D.ms_assemble = A.ms_deriv + A.ms_asm;
      if( _determinacy_verdict( A, D ) ){
        // Marginal gap -- on ANY verdict: re-evaluate at the same point scaled by
        // 1.3 (deterministic).  On a deficient verdict, agreement keeps it and a
        // change is reported as STATE CONTENT with the flag set.  On a DETERMINED
        // verdict a small gap means near-singular at the reference point (PDE40b:
        // sigma_min = 20 x threshold at n = 1008); the second point tells whether
        // that is the formulation or the point.  Both are printed either way.
        if( D.gap < _knob_DETERMINACY_GAP_MIN() ){
          std::vector<double> x2( x ); for( double& v : x2 ) v *= 1.3;
          t_JacAssembly A2;
          bool const ok2 = marching ? assemble_window( w, x2.data(), u.empty() ? nullptr : u.data(), cst, A2, eq )
                                    : assemble( x2.data(), u.empty() ? nullptr : u.data(), cst, A2, eq );
          t_Determinacy D2 = D;
          if( ok2 && _determinacy_verdict( A2, D2 ) ){
            D.marginal_reeval = true;
            D.rank2 = D2.rank; D.gap2 = D2.gap; D.verdict2 = D2.verdict;
            if( D2.verdict != D.verdict
             && ( D.def_right > D.free_expected || D2.def_right > D2.free_expected ) )
              D.verdict = t_Determinacy::STATE_CONTENT;
          }
        }
      }
      else D.ran = false;
    }
    // rev137 instrument (CRONOS_AUDIT_DROPNULL=1): WHICH dropped claims does the null
    // space need back?  Continuity is linear in var, so the residuals of the dropped
    // claim rows on the right-null vectors form W (n_dropped x n_null).  The number of
    // claims that must be restored is rank(W); a minimal set is the pivot rows of a
    // column-pivoted QR of W^T.  Identified at setup, solver-free, without a re-derive.
    // Dense SVD (cap) here; the QR backends can supply the null vectors from their dead
    // columns.  MEASURED: PDE6_blk0 W is 4x4 full rank -> the 4 dropped claims are the
    // 4 null directions, exactly.
    static bool const kDropNull = Options::_env_flag( "CRONOS_AUDIT_DROPNULL", false );
    if( kDropNull && D.ran && D.verdict == t_Determinacy::STATE_CONTENT
     && !_strong_dropped_claims.empty() && std::max( A.n, A.m ) <= _knob_DETERMINACY_MAX_DENSE() ){
      arma::mat U, Vm; arma::vec sv;
      if( arma::svd( U, sv, Vm, arma::mat( A.J ) ) ){
        std::vector<arma::uword> nullk;
        for( arma::uword k = 0; k < sv.n_elem; ++k ) if( sv(k) <= D.tol_abs ) nullk.push_back( k );
        std::vector<t_FrozenInterfaceClaimKey> keys( _strong_dropped_claims.begin(), _strong_dropped_claims.end() );
        arma::mat W( (arma::uword)keys.size(), (arma::uword)nullk.size(), arma::fill::zeros );
        std::vector<std::pair<size_t,double>> coeffs;
        for( size_t i = 0; i < keys.size(); ++i ){
          t_InterfaceClaimKey const claim = _unfreeze_interface_claim_key( keys[i] );
          coeffs.clear();
          if( !_trace_claim_coeffs( claim, coeffs ) ) continue;
          double scale = 0.;
          for( auto const& [col, c] : coeffs ) if( col < A.m ) scale += std::fabs( c );
          for( size_t kk = 0; kk < nullk.size(); ++kk ){
            double r = 0.;
            for( auto const& [col, c] : coeffs ) if( col < A.m ) r += c * Vm( col, nullk[kk] ) * A.cw[col];
            W( i, kk ) = r / std::max( scale, 64. * DBL_EPSILON );
          }
        }
        arma::mat Q, R; arma::uvec P;
        size_t needed = 0; std::vector<size_t> sel;
        if( W.n_rows && W.n_cols && arma::qr( Q, R, P, W.t(), "vector" ) ){
          double const r0 = ( R.n_rows && R.n_cols ) ? std::fabs( R(0,0) ) : 0.;
          for( arma::uword i = 0; i < std::min( R.n_rows, R.n_cols ); ++i )
            if( std::fabs( R(i,i) ) > 1e-8 * r0 ){ ++needed; sel.push_back( (size_t)P(i) ); }
        }
        _disp( 3 ) << "  [drop-null] win=" << w << " null vectors=" << nullk.size() << " dropped claims=" << keys.size()
                  << " rank(W)=" << needed << " -> claims that must be restored=" << needed << ":";
        for( size_t i : sel ){
          auto const& ck = keys[i];
          _disp( 3 ) << "  [st=" << std::get<1>(ck) << " dom=" << std::get<3>(ck) << " iel=" << std::get<4>(ck)
                    << " face=" << std::get<5>(ck) << " kind=" << std::get<6>(ck) << "]";
        }
        _disp( 3 ) << std::endl;
      }
    }
    _determinacy.push_back( D );

    // ONE LINE PER WINDOW, unconditional on DISPLAY_LEVEL.  Silence is never a
    // verdict: a full-rank system prints DETERMINED explicitly (the solve-time
    // spectrum audit was found to print nothing on full-rank instances, which
    // made "no line" indistinguishable from "did not run").
    // Stream state is saved/restored: the rev133 sweep's only numerics diff against
    // rev132 was a driver table's column justification inherited from this print.
    if( !_auditQuiet ) _print_determinacy_window( D, w );
  }
}


// ======================================================================
// OCFESLV::classify_pde
// ======================================================================


#ifdef MC__OCFESLV_REFERENCE_ROBUSTNESS_PROBE
#endif // MC__OCFESLV_REFERENCE_ROBUSTNESS_PROBE


//-----------------------------------------------------------------------------
// OCFESLV::_validate_hyperbolic_incoming_bcs
//
// Auto-closure supplies the OUTGOING characteristic rows but assumes the user
// placed the INCOMING boundary conditions correctly.  Nothing else checks that,
// so a mis-placed inflow condition (right count, wrong characteristic / wrong
// end) silently yields a rank-deficient system that diverges instead of failing
// cleanly.  This guard, per EVOL_HYPERBOLIC block and spatial face, verifies:
//   (1) COUNT: #user BOUNDARY rows at the face == rank of the incoming projector
//       (A+ at LB / A- at UB);
//   (2) DIRECTION: each BC's state coupling d(BC)/d(state) lies in the incoming
//       characteristic subspace rowspace(incoming projector) = left-eigenvector
//       span (V columns of its SVD) -- i.e. the BC prescribes an incoming, not
//       outgoing, characteristic.  (Rowspace, not colspace: the coupling is the
//       LEFT eigenvector; the two coincide only for a symmetric symbol.)
// Section C of PDE14 (a<0 with forward-placed BCs) is the oracle: it must turn a
// divergent solve into a clean setup rejection here.  PDE15 (non-symmetric) is
// the oracle for the rowspace/colspace distinction: with colspace the guard
// wrongly rejects the correct left-eigenvector BCs.
//
// FLUX / DERIVATIVE BCs: the direction test extracts the VALUE coupling
// a = d(BC)/d(u) via DAG FAD, which is undefined through OpP (FFPartial).  Under
// the flux-BC guard (always on since rev319) the surviving first-order
// OpP(u_j,{x_i,1}) nodes are proxied with the principal-symbol derivative proxies
// vDerivProxy[i*nState+j] (so FAD is well-defined) and the proxies seeded to a
// nonzero nominal (1.0).  The acceptance test is CRITERION #3 (value OR
// derivative): a BC pins an incoming mode either through its value coupling
// a = d(BC)/d(u) OR through its face-normal derivative coupling b = d(BC)/d(d_z u)
// (FAD w.r.t. the face-normal proxies).  ACCEPT iff a in Vin OR b in Vin; REJECT
// only when the BC constrains something yet touches NO incoming characteristic in
// either coupling.  Rationale: the spectral boundary derivative is a full nodal
// combination that includes the boundary value, so a derivative coupling on an
// incoming characteristic legitimately pins it (e.g. the bilinear c*d_z P, where
// a=(0,1) is outgoing but b=(c_ref,0) is incoming -- PDE18 B, well-posed, recovers).
// This b-branch RESCUES (OR); it is the MIRROR of the removed b-test, which wrongly
// made b a REQUIREMENT (AND) and false-rejected PDE16 section C.  Genuinely
// ill-posed forms touch no incoming characteristic in either coupling: PDE16 D
// (w_out + b d_z w_out, a and b both outgoing) and PDE18 B2 (c*d_z c, pure
// outgoing) are rejected.  Seeding to 1.0 is exact in DIRECTION for a single-
// derivative flux because the test is scale-invariant.  Without the macro, OpP BCs
// skip the test and resid_b stays trivial so the logic reduces EXACTLY to the
// value-only test for plain value BCs (count still enforced).  FFIntegral
// (nonlocal) BCs always skip.
//-----------------------------------------------------------------------------


////////////////////////////////////////////////////////////////////////////////
// OCFESLV::duplicate_node_spread / report_duplicate_node_spreads
//
// Geometric coincidence is detected by rounding each node coordinate to 1e-12 and using
// the rounded tuple as a map key -- the same rule the driver-side implementations have
// used since the arc began, kept identical deliberately so that moving the readout into
// the header cannot silently change a number a driver has been gating on.
//
// The `off` argument is the state's base offset into the flat collocation vector.  The
// reporter walks states in states_colloc() order accumulating node_colloc().size(), which
// is the same layout pos_state() produces; a caller that only wants one state can pass its
// own offset.
////////////////////////////////////////////////////////////////////////////////


inline void
OCFESLV::_report_beta_gated
( double const* var, char const* label, bool converged )
const
{
  if( !_issetup || !var || !_betaReady ) return;
  if( _betaZ.n_rows == 0 || _betaZ.n_cols == 0 ) return;
  // rev63 CORRECTION.  rev62 guarded this on `_inMarch`, which is DEAD CODE at the site that
  // needed it: solve_marching clears _inMarch before emitting its assembled-trajectory
  // report, so the guard never fired and PDE5's march row was reported anyway (it happened to
  // agree with the monolithic value, which is exactly how a silent wrong reading survives).
  // The refusal now lives at the call site, where the caller knows what it is holding.

  size_t const m  = _betaColOf.size();
  size_t const kk = (size_t)_betaZ.n_cols;
  if( (size_t)_betaZ.n_rows != m ){
    _disp( 3 ) << "  [beta] Z rows=" << _betaZ.n_rows << " but colof size=" << m
              << " -- not reporting" << std::endl;
    return;
  }

  // lambda from the solved var[]: under IC_TRACE every tau IS an ordinary column of var[].
  // var[] under IC_TRACE is [ state coefficients | tau columns ], so the valid upper bound
  // for a tau column is n_colloc_sta() + n_colloc_trace().  n_colloc_trace() is 0 under
  // IC_STRONG, which is correct here: _betaReady is only ever set on the IC_TRACE path.
  size_t const ncol = n_colloc_sta() + n_colloc_trace();
  arma::vec lam( m, arma::fill::zeros );
  size_t n_oob = 0;
  for( size_t j = 0; j < m; ++j ){
    if( _betaColOf[j] >= ncol ){ ++n_oob; continue; }
    lam[j] = var[ _betaColOf[j] ];
  }
  if( n_oob ){
    _disp( 3 ) << "  [beta] " << n_oob << " of " << m
              << " tau columns out of range -- not reporting" << std::endl;
    return;
  }

  arma::vec const beta = _betaZ.t() * lam;          // beta = Z^T lambda
  double const nbeta = arma::norm( beta, 2 );
  double const nlam  = arma::norm( lam,  2 );

  // rev63: split the absorbed violation between the two kernel families.  The hypothesis
  // this tests is that ALL of it lives in the cross-domain directions -- the family unique to
  // the moving mesh -- and none in the lifting-alias family that PDE2/PDE5/PDE9 also carry
  // and that measures zero violation on them (beta ~ 1e-16 with beta/lambda ~ 1e-8..1e-11).
  // If confirmed, the promotion criterion narrows from "any Z-support tau" to "taus in a
  // cross-domain null component", which is a strictly smaller and better-motivated set.
  double b2cross = 0., b2intra = 0.;
  size_t ncross = 0, nintra = 0;
  bool const have_class = ( _betaDirCross.size() == (size_t)beta.n_elem );
  if( have_class )
    for( arma::uword j = 0; j < beta.n_elem; ++j ){
      double const b2 = beta[j]*beta[j];
      if( _betaDirCross[(size_t)j] ){ b2cross += b2; ++ncross; }
      else                          { b2intra += b2; ++nintra; }
    }
  double const nbcross = std::sqrt( b2cross ), nbintra = std::sqrt( b2intra );

  // rev64 TAG FIX.  rev63 tagged "[violation is ENTIRELY cross-domain]" whenever the intra
  // energy was negligible -- which is trivially true when nintra == 0, REGARDLESS of whether
  // there is any violation at all.  PDE7 fired it at ||beta|| = 7.1e-16 (beta/lambda = 6e-12),
  // i.e. a share statement about a quantity that is zero, reading as a finding.  A share is
  // only meaningful once the violation clears round-off, measured RELATIVE to ||lambda||
  // because beta is a component of it: gauge models sit at 1e-8..1e-12, MMPDE13 at 0.07..0.9.
  bool const gauge = !( nbeta > 1.0e-8 * std::max( nlam, 1.0e-300 ) );
  // W-component by difference: ||lambda||^2 = ||W^T lambda||^2 + ||beta||^2 (Z orthonormal,
  // W its orthogonal complement).  Reported for scale, so a large ||beta|| cannot be read as
  // large when lambda itself is large.
  double const nW = std::sqrt( std::max( 0.0, nlam*nlam - nbeta*nbeta ) );

  _disp( 3 ) << "  [beta] " << ( label ? label : "" )
            << "  conv=" << ( converged ? "yes" : "NO" )
            << "  m=" << m << " k=" << kk
            << "  eps=" << std::scientific << std::setprecision(3) << _betaEps
            << "  ||lambda||=" << nlam
            << "  ||W^T lambda||=" << nW
            << "  ||beta||=" << nbeta
            << ( converged ? "  VIOLATION eps*||beta||=" : "  eps*||beta||=" )
            << ( _betaEps * nbeta )
            << std::defaultfloat << std::endl;
  if( have_class )
    _disp( 3 ) << "  [beta-split] directions: cross=" << ncross << " intra=" << nintra
              << "  ||beta_cross||=" << std::scientific << std::setprecision(3) << nbcross
              << "  ||beta_intra||=" << nbintra
              << "  cross share=" << ( nbeta > 0. ? ( b2cross / ( nbeta*nbeta ) ) : 0. )
              << std::defaultfloat
              << ( gauge
                     ? "  [GAUGE: violation at round-off -- the share is a ratio of noise and"
                       " means nothing]"
                   : ( b2intra <= 1e-20 * nbeta*nbeta )
                     ? "  [violation is ENTIRELY cross-domain]" : "" )
              << std::endl;
  else if( converged )
    _disp( 3 ) << "  [beta-split] direction classes unavailable (size mismatch "
              << _betaDirCross.size() << " vs " << beta.n_elem << ")" << std::endl;

  if( converged )
    _disp( 3 ) << "  [beta] eps*||beta|| = ||Z^T r_c^raw||: the continuity violation this mode"
                 " ABSORBS rather than enforces.  Round-off => the claim set is consistent and"
                 " both exact modes impose it; O(1) => IC_TRACE converges by absorbing a"
                 " conflict IC_STRONG has nowhere to put.  Sweep INTERFACE.TRACE_PROJ_EPS to confirm the"
                 " product is eps-independent before reading it." << std::endl;
  else
    _disp( 3 ) << "  [beta] *** NOT A VIOLATION READING -- THE SOLVE DID NOT CONVERGE. ***  The"
                 " identity eps*beta = -Z^T r_c^raw holds AT A SOLUTION; at this iterate"
                 " eps*beta = -Z^T r_c^raw + Z^T R, and Z^T R is the unconverged residual in"
                 " the SAME subspace being measured.  The number above is the sum of the two"
                 " and cannot be separated.  Re-run on a mesh where this mode converges."
              << std::endl;
}

inline void
OCFESLV::_report_dup_spread_gated
( double const* var, char const* label )
const
{
  if( !_issetup || !var ) return;

  char const* mode = "IC_WEAK";
  bool exact = false;
  switch( options.INTERFACE.IMPOSITION ){
    case Options::IC_TRACE:  mode = "IC_TRACE";  exact = true;  break;
    case Options::IC_STRONG: mode = "IC_STRONG"; exact = true;  break;
    default:                 mode = "IC_WEAK";   exact = false; break;
  }

  std::vector<double> const v( var, var + n_colloc_sta() );

  std::ostringstream lab;
  lab << ( label ? label : "" ) << "  mode=" << mode;
  double const worst = report_duplicate_node_spreads( v, lab.str().c_str(), &std::cerr );

  // Non-fatal warning, EXACT modes only.  Under IC_WEAK continuity is imposed by penalty and
  // the spread is truncation limited -- PDE2 IC_WEAK sits at 9.458e-11, which is correct
  // behaviour and would fire this every time.  Under an exact mode the same number would mean
  // a continuity condition has quietly stopped being enforced, which is the Stage A failure
  // mode and is invisible to k, SYMMETRIC, the residual and the PASS flag alike.
  if( exact && worst > _knob_DUP_SPREAD_WARN() )
    _disp( 3 ) << "  [dup-spread] *** WARNING: worst max_pair=" << std::scientific
              << std::setprecision(6) << worst << " exceeds " << _knob_DUP_SPREAD_WARN()
              << " under " << mode << ".  Exact imposition is round-off limited at any mesh"
                 " (~1e-17), so this is a continuity condition that is NOT being enforced --"
                 " not a truncation effect and not a mesh-resolution question.  Report only:"
                 " this does not fail the solve. ***" << std::defaultfloat << std::endl;
}

template <typename U>
inline U
OCFESLV::eval_colloc
( FFVar const& Var, t_Coord const& pt,
  U const* val_state, U const* val_input, U const* val_cst, t_Side const* side )
const
{
  auto side_of = [&]( FFVar const& d ){ if( !side ) return -1; auto const i = side->find( d );   // -1: legacy
                                        return i == side->end()? -1: i->second; };
  if( !_issetup ) throw Exceptions( Exceptions::SETUP );

  auto itdep = _mDep.find( Var.id().second );
  if( itdep == _mDep.end() ) throw Exceptions( Exceptions::INDEX );

  int const    typ = itdep->second.first;
  size_t const pos = itdep->second.second;

  auto coord_value = [&]( FFVar const& dvar ) -> double
  {
    auto it = pt.find( dvar );
    if( it != pt.end() ) return it->second;

    // Be tolerant of user-DAG variables passed after setup(), when the
    // internal collocation model may be held on a private copied DAG.
    for( auto const& kv : pt )
      if( kv.first.id().second == dvar.id().second ) return kv.second;

    throw Exceptions( Exceptions::INDEX );
  };

  switch( typ ){
    case DOMAIN:
      if( pos >= _vDom.size() ) throw Exceptions( Exceptions::INDEX );
      return U( coord_value( _vDom[pos] ) );

    case CONSTANT:
      if( pos >= _vCst.size() ) throw Exceptions( Exceptions::INDEX );
      if( val_cst ) return val_cst[pos];
      if( pos < _vCstVal.size() ) return U( _vCstVal[pos] );
      throw Exceptions( Exceptions::CSTVAL );

    case STATE:
    case INPUT:
      break;

    default:
      throw Exceptions( Exceptions::INDEX );
  }

  bool const is_state = ( typ == STATE );
  FFVar const& v      = is_state ? _vVar.at( pos ) : _vInp.at( pos );
  U const* data       = is_state ? val_state : val_input;
  if( !data ) throw Exceptions( Exceptions::INTERNAL );

  auto const& dep = is_state ? _mVar.at( v ) : _mInp.at( v );
  auto const& blk = is_state ? _blkVar.at( v ) : _blkInp.at( v );
  size_t const blk_size = std::get<1>( blk );

  std::map<FFVar,size_t,lt_FFVar> ndx_el;

  // --- Marching read routing: after a SOLVE.MARCHING solve the resident evolution domain is one
  //     sliding window and the passed buffer holds only the last window.  If this state/input
  //     depends on the evolution direction, read from the trajectory block whose [t0,t1] contains
  //     pt[evolution].  Side-effect-free: the block window only forms the reference coordinate; the
  //     resident single-element layout matches every stored block, so pos_state()/pos_input() index
  //     the block buffer directly.  (Trajectory values are double; a non-double U is value-only,
  //     which the derivative paths do not use here.)
  MarchBlock const* mblk = nullptr;
  if( _lastSolveMarched && _evolution_dom_var.dag() ){
    bool dep_on_evo = false;
    for( auto const& d : dep )
      if( d.id().second == _evolution_dom_var.id().second ){ dep_on_evo = true; break; }
    if( dep_on_evo ){
      double const tval = coord_value( _evolution_dom_var );
      if( !_marchTrajectory.empty() ){
        // Trajectory stored (default): route to the block whose window [t0,t1] contains tval.
        double const tol = 64.*DBL_EPSILON*std::max( 1., std::fabs( tval ) );
        for( auto const& b : _marchTrajectory )
          if( tval >= b.t0 - tol && tval <= b.t1 + tol ){ mblk = &b; break; }
        if( mblk && side_of( _evolution_dom_var ) == FFDom::PLUS && std::fabs( tval - mblk->t1 ) <= tol )
          for( auto const& b : _marchTrajectory )           // tau^+ at a window seam: the window STARTING there
            if( std::fabs( tval - b.t0 ) <= tol ){ mblk = &b; break; }
        if( !mblk ) throw Exceptions( Exceptions::INDEX );   // evolution coordinate outside [0,T]
      }
      else{
        // Marched with OUTPUT.MARCH_STORE disabled: only the last window is resident in the
        // passed buffer.  An in-window read falls through and reads it; an out-of-window read cannot
        // be served, so fail with an explicit diagnostic instead of a bare INDEX abort downstream.
        FFDom const& edom = _mDom.at( _evolution_dom_var );
        double const tol  = 64.*DBL_EPSILON*std::max( 1., std::max( std::fabs(edom.lo_dom), std::fabs(edom.up_dom) ) );
        if( tval < edom.lo_dom - tol || tval > edom.up_dom + tol ){
          _disp( 3 ) << "OCFESLV::eval_colloc: evolution coordinate " << tval << " lies outside the resident "
                       "window [" << edom.lo_dom << "," << edom.up_dom << "] after a marching solve with "
                       "OUTPUT.MARCH_STORE disabled -- only the last window is available.  Enable "
                       "options.OUTPUT.MARCH_STORE (the default) to read the full marched field, or "
                       "read during the march / use terminal_value()." << std::endl;
          throw Exceptions( Exceptions::INDEX );   // opt-out: full-field marched read unavailable
        }
      }
    }
  }

  // Lumped state/input: one coefficient, no interpolation.
  if( dep.empty() ){
    size_t const inc = is_state ? pos_state( v, ndx_el )
                                : pos_input( v, ndx_el );
    return data[inc];
  }

  struct DimEval
  {
    FFVar const* pvar;
    FFDom const* pdom;
    double       xi;
    size_t       cur_n;
  };

  std::vector<DimEval> dims;
  dims.reserve( dep.size() );

  for( auto const& dvar : dep ){
    FFDom const& dom = is_state ? _mDom.at( dvar ) : _input_dom( v, dvar );
    size_t iel; double xi;
    if( mblk && dvar.id().second == _evolution_dom_var.id().second ){
      // Evolution dim: map onto this block's window [t0,t1] (single element -> iel=0); the
      // reference nodes come from the resident template (window-position-independent in xi).
      double const w = mblk->t1 - mblk->t0;
      iel = 0;
      xi  = ( w > 0. ? 2.*( coord_value( dvar ) - mblk->t0 )/w - 1. : 0. );
    }
    else{
      auto const pr = _point_in_domain( dom, coord_value( dvar ), side_of( dvar ) );
      iel = pr.first; xi = pr.second;
    }
    ndx_el[dvar] = iel;
    dims.push_back( { &dvar, &dom, xi, dom.n_node } );
  }

  size_t const inc = is_state ? pos_state( v, ndx_el )
                              : pos_input( v, ndx_el );
  std::vector<U> buf;
  if( mblk ){
    auto const& src = is_state ? mblk->var : mblk->inp;   // double block from the trajectory
    if( inc + blk_size > src.size() ) throw Exceptions( Exceptions::INTERNAL );
    buf.assign( src.begin() + inc, src.begin() + inc + blk_size );  // double -> U (value)
  }
  else
    buf.assign( data + inc, data + inc + blk_size );

  // Collapse one tensor-product dimension at a time.  The dimension ordering
  // matches dep, pos_state()/pos_input(), node_colloc(), and OCVar flattening:
  // the first domain dimension varies fastest.
  for( size_t idim = 0; idim < dims.size(); ++idim ){
    auto& dim = dims[idim];

    size_t stride = 1;
    for( size_t k = 0; k < idim; ++k ) stride *= dims[k].cur_n;

    size_t outer = 1;
    for( size_t k = idim+1; k < dims.size(); ++k ) outer *= dims[k].cur_n;

    std::vector<U> newbuf;
    if( !dim.pdom->eval_lagrange( newbuf, dim.xi, 0, buf, stride, outer ) )
      throw Exceptions( Exceptions::INTERNAL );

    buf = std::move( newbuf );
    dim.cur_n = 1;
  }

  if( buf.size() != 1 ) throw Exceptions( Exceptions::INTERNAL );
  return buf[0];
}


template <typename U>
inline U
OCFESLV::eval_colloc_deriv
( FFVar const& Var, t_Coord const& pt, FFVar const& wrt, int ORD,
  U const* val_state, U const* val_input, U const* val_cst, t_Side const* side )
const
{
  auto side_of = [&]( FFVar const& d ){ if( !side ) return -1; auto const i = side->find( d );   // -1: legacy
                                        return i == side->end()? -1: i->second; };
  if( !_issetup ) throw Exceptions( Exceptions::SETUP );
  if( ORD <= 0 )
    return eval_colloc<U>( Var, pt, val_state, val_input, val_cst, side );

  auto itdep = _mDep.find( Var.id().second );
  if( itdep == _mDep.end() ) throw Exceptions( Exceptions::INDEX );

  int const    typ = itdep->second.first;
  size_t const pos = itdep->second.second;

  auto coord_value = [&]( FFVar const& dvar ) -> double
  {
    auto it = pt.find( dvar );
    if( it != pt.end() ) return it->second;
    for( auto const& kv : pt )
      if( kv.first.id().second == dvar.id().second ) return kv.second;
    throw Exceptions( Exceptions::INDEX );
  };

  switch( typ ){
    case DOMAIN:
      // d^ORD of a coordinate w.r.t. wrt: 1 if wrt is this coordinate and ORD==1, else 0.
      if( pos >= _vDom.size() ) throw Exceptions( Exceptions::INDEX );
      return U( ( _vDom[pos].id().second == wrt.id().second && ORD == 1 )? 1.: 0. );

    case CONSTANT:
      return U( 0. );

    case STATE:
    case INPUT:
      break;

    default:
      throw Exceptions( Exceptions::INDEX );
  }

  bool const is_state = ( typ == STATE );
  FFVar const& v      = is_state ? _vVar.at( pos ) : _vInp.at( pos );
  U const* data       = is_state ? val_state : val_input;
  if( !data ) throw Exceptions( Exceptions::INTERNAL );

  auto const& dep = is_state ? _mVar.at( v ) : _mInp.at( v );
  auto const& blk = is_state ? _blkVar.at( v ) : _blkInp.at( v );
  size_t const blk_size = std::get<1>( blk );

  // Lumped, or not distributed along wrt: derivative is identically zero.
  bool wrt_in_dep = false;
  for( auto const& dvar : dep )
    if( dvar.id().second == wrt.id().second ){ wrt_in_dep = true; break; }
  if( dep.empty() || !wrt_in_dep ) return U( 0. );

  // Marching read routing (mirrors eval_colloc): read from the trajectory block containing
  // pt[evolution].  The derivative rescale along the evolution direction must use the BLOCK's
  // window width (t1-t0), not the resident window's -- carried in dim.w below.
  MarchBlock const* mblk = nullptr;
  if( _lastSolveMarched && _evolution_dom_var.dag() ){
    bool dep_on_evo = false;
    for( auto const& d : dep )
      if( d.id().second == _evolution_dom_var.id().second ){ dep_on_evo = true; break; }
    if( dep_on_evo ){
      double const tval = coord_value( _evolution_dom_var );
      if( !_marchTrajectory.empty() ){
        double const tol = 64.*DBL_EPSILON*std::max( 1., std::fabs( tval ) );
        for( auto const& b : _marchTrajectory )
          if( tval >= b.t0 - tol && tval <= b.t1 + tol ){ mblk = &b; break; }
        if( mblk && side_of( _evolution_dom_var ) == FFDom::PLUS && std::fabs( tval - mblk->t1 ) <= tol )
          for( auto const& b : _marchTrajectory )           // tau^+ at a window seam: the window STARTING there
            if( std::fabs( tval - b.t0 ) <= tol ){ mblk = &b; break; }
        if( !mblk ) throw Exceptions( Exceptions::INDEX );   // evolution coordinate outside [0,T]
      }
      else{
        FFDom const& edom = _mDom.at( _evolution_dom_var );
        double const tol  = 64.*DBL_EPSILON*std::max( 1., std::max( std::fabs(edom.lo_dom), std::fabs(edom.up_dom) ) );
        if( tval < edom.lo_dom - tol || tval > edom.up_dom + tol ){
          _disp( 3 ) << "OCFESLV::eval_colloc_deriv: evolution coordinate " << tval << " lies outside the "
                       "resident window after a marching solve with OUTPUT.MARCH_STORE disabled -- "
                       "enable options.OUTPUT.MARCH_STORE (the default)." << std::endl;
          throw Exceptions( Exceptions::INDEX );
        }
      }
    }
  }

  struct DimEval
  {
    FFVar const* pvar;
    FFDom const* pdom;
    double       xi;
    size_t       iel;
    double       w;      // physical element width for the derivative rescale (block window if marched)
    size_t       cur_n;
  };

  std::vector<DimEval> dims;
  dims.reserve( dep.size() );
  std::map<FFVar,size_t,lt_FFVar> ndx_el;

  for( auto const& dvar : dep ){
    FFDom const& dom = is_state ? _mDom.at( dvar ) : _input_dom( v, dvar );
    size_t iel; double xi, w;
    if( mblk && dvar.id().second == _evolution_dom_var.id().second ){
      double const wblk = mblk->t1 - mblk->t0;
      iel = 0;
      xi  = ( wblk > 0. ? 2.*( coord_value( dvar ) - mblk->t0 )/wblk - 1. : 0. );
      w   = wblk;   // rescale along the evolution dim uses the block window width
    }
    else{
      auto const pr = _point_in_domain( dom, coord_value( dvar ), side_of( dvar ) );
      iel = pr.first; xi = pr.second;
      w   = dom.elem_len.at( iel );
    }
    ndx_el[dvar] = iel;
    dims.push_back( { &dvar, &dom, xi, iel, w, dom.n_node } );
  }

  size_t const inc = is_state ? pos_state( v, ndx_el )
                              : pos_input( v, ndx_el );
  std::vector<U> buf;
  if( mblk ){
    auto const& src = is_state ? mblk->var : mblk->inp;
    if( inc + blk_size > src.size() ) throw Exceptions( Exceptions::INTERNAL );
    buf.assign( src.begin() + inc, src.begin() + inc + blk_size );
  }
  else
    buf.assign( data + inc, data + inc + blk_size );

  // Collapse one tensor-product dimension at a time (same ordering as eval_colloc):
  // apply the order-ORD differential along wrt with the per-element physical rescale,
  // and order 0 (plain interpolation) along all other directions.
  for( size_t idim = 0; idim < dims.size(); ++idim ){
    auto& dim = dims[idim];

    size_t stride = 1;
    for( size_t k = 0; k < idim; ++k ) stride *= dims[k].cur_n;

    size_t outer = 1;
    for( size_t k = idim+1; k < dims.size(); ++k ) outer *= dims[k].cur_n;

    int const ord_dim = ( dim.pvar->id().second == wrt.id().second )? ORD : 0;
    double const rescale = ord_dim
      ? std::pow( 2. / dim.w, ord_dim )
      : 1.;

    std::vector<U> newbuf;
    if( !dim.pdom->eval_lagrange( newbuf, dim.xi, ord_dim, buf, stride, outer, DBL_EPSILON, rescale ) )
      throw Exceptions( Exceptions::INTERNAL );

    buf = std::move( newbuf );
    dim.cur_n = 1;
  }

  if( buf.size() != 1 ) throw Exceptions( Exceptions::INTERNAL );
  return buf[0];
}


// rev152a: EqnRole <-> PlanRowRole bridge.  The contract header cannot see EqnRole
// (it lives in OCFESLV), so it mirrors the enum; these asserts pin the mirror.
static_assert( static_cast<int>( OCFESLV::EqnRole::AUTO )       == static_cast<int>( PlanRowRole::AUTO ),       "PlanRowRole mirror" );
static_assert( static_cast<int>( OCFESLV::EqnRole::INTERIOR )   == static_cast<int>( PlanRowRole::INTERIOR ),   "PlanRowRole mirror" );
static_assert( static_cast<int>( OCFESLV::EqnRole::INITIAL )    == static_cast<int>( PlanRowRole::INITIAL ),    "PlanRowRole mirror" );
static_assert( static_cast<int>( OCFESLV::EqnRole::BOUNDARY )   == static_cast<int>( PlanRowRole::BOUNDARY ),   "PlanRowRole mirror" );
static_assert( static_cast<int>( OCFESLV::EqnRole::INTERFACE )  == static_cast<int>( PlanRowRole::INTERFACE ),  "PlanRowRole mirror" );
static_assert( static_cast<int>( OCFESLV::EqnRole::LINK )       == static_cast<int>( PlanRowRole::LINK ),       "PlanRowRole mirror" );
static_assert( static_cast<int>( OCFESLV::EqnRole::SURFACE )    == static_cast<int>( PlanRowRole::SURFACE ),    "PlanRowRole mirror" );
static_assert( static_cast<int>( OCFESLV::EqnRole::DIAGNOSTIC ) == static_cast<int>( PlanRowRole::DIAGNOSTIC ), "PlanRowRole mirror" );


// ======================================================================
// OCFESLV::_sat_claim_coeffs
// ======================================================================
// Delegates to _trace_claim_coeffs for the unscaled Lagrange weights,
// then multiplies each by prefactor.  Keeping a single implementation of
// the face-decode / stride / Lagrange logic in _trace_claim_coeffs means
// any future change propagates automatically to both paths.
// ======================================================================


// _clear_frozen_interface_consumers is defined inline in the class body.


// _validate_interface_plan_draft is defined inline in the class body.

// _commit_interface_plan_draft is defined inline in the class body.


// ---------------------------------------------------------------------------
// Post-solve over-drop verification gate (refactor item 10, eval stage).
//
// _strong_dropped_claims is populated by the IC_STRONG explicit-tau setup
// detector; the gate is unguarded and trivially returns OK when nothing was
// dropped.  It needs no committed-plan-only state: the continuity residual of
// a dropped claim is a function of the physical DOFs alone, recovered via the
// same _trace_claim_coeffs() that built the (now-absent) trace row.
// ---------------------------------------------------------------------------
inline OCFESLV::InterfaceDropStatus
OCFESLV::verify_interface_drop( double const* var, double rel_tol,
                              std::set<t_FrozenInterfaceClaimKey>* overdropped ) const
{
  if( overdropped ) overdropped->clear();
  if( !var || _strong_dropped_claims.empty() )
    return InterfaceDropStatus::OK;

  double              worst_rel = 0.0;
  t_InterfaceClaimKey worst_claim{};
  size_t              n_checked = 0, n_skipped = 0, n_over = 0;

  std::vector<std::pair<size_t,double>> coeffs;
  // rev124: NAME every dropped claim, not just the worst one.
  //
  // The summary below reports a count and the single worst violation.  That is enough
  // to answer "was the drop safe" and NOT enough to answer "WHICH claims were dropped",
  // which is the open question on OCFE_PDE6_blk0: the redundancy detector fires on the
  // hand-written formulation, removes 4 claims (census 276 -> tau 272), and leaves a
  // rank deficiency of exactly 4.  Whether those are the SAME 4 decides whether the
  // detector is under-dropping or dropping the wrong ones, and no existing output says.
  //
  // The overdropped-set overload cannot serve a driver here: t_FrozenInterfaceClaimKey
  // is PROTECTED, so OCFE_PDE6_blk0 has to pass nullptr.  Naming belongs on this side.
  static bool const kDropSet = Options::_env_flag( "CRONOS_AUDIT_DROPSET", false );
  if( kDropSet )
    _disp( 3 ) << "  [dropset] " << _strong_dropped_claims.size()
              << " claim(s) removed by the redundancy detector"
                 " (rel = continuity violation on the converged solution;"
                 " small means the drop was sound):" << std::endl;

  for( auto const& fck : _strong_dropped_claims ){
    t_InterfaceClaimKey const claim = _unfreeze_interface_claim_key( fck );
    coeffs.clear();
    // Rebuild the dropped claim's continuity row in DOF space.  Failure means
    // the claim geometry is no longer resolvable -> cannot assess; count it.
    if( !_trace_claim_coeffs( claim, coeffs ) || coeffs.empty() ){
      ++n_skipped;
      continue;
    }

    // Continuity residual r = sum_k coeff_k * var[col_k].  Relative scale is
    // the sum of absolute term magnitudes (= |value_lo| + |value_hi| for a C0
    // claim), which self-normalises value- and flux-form claims alike.
    double r = 0.0, scale = 0.0;
    bool   col_ok = true;
    for( auto const& [col, c] : coeffs ){
      if( col >= _nCollVar ){ col_ok = false; break; }  // physical DOFs only
      double const term = c * var[col];
      r     += term;
      scale += std::fabs( term );
    }
    if( !col_ok ){ ++n_skipped; continue; }

    ++n_checked;
    double const rel = std::fabs( r ) / std::max( scale, 64. * DBL_EPSILON );

    if( kDropSet ){
      std::string dn = "dom?";
      for( auto const& kv : _mDom )
        if( static_cast<size_t>( kv.first.id().second ) == claim.dom_id ){
          dn = kv.first.name(); break; }
      _disp( 3 ) << "      state=" << state_name_of_id( claim.state_id )
                << " dir=" << dn
                << " iel_lo=" << claim.iel_lo
                << " face="   << claim.face
                << " kind="   << claim.kind
                << "   rel="  << std::scientific << std::setprecision(3) << rel
                << std::defaultfloat << std::endl;
    }
    if( rel > worst_rel ){ worst_rel = rel; worst_claim = claim; }
    if( rel > rel_tol ){ ++n_over; if( overdropped ) overdropped->insert( fck ); }
  }

  // Always emit a one-line summary so rel_tol can be calibrated from the log:
  // a sound drop floors worst_rel at the discretization truncation scale, a
  // wrong drop (PDE7) pushes it to O(1).
  _disp( 3 ) << "OCFESLV::verify_interface_drop: checked " << n_checked
            << " dropped continuity claim(s)"
            << ( n_skipped ? ( " (" + std::to_string( n_skipped )
                               + " unresolvable, skipped)" ) : std::string() )
            << ", worst relative continuity violation=" << worst_rel
            << " (rel_tol=" << rel_tol << ")." << std::endl;

  if( n_over ){
    _disp( 2 ) << "OCFESLV::verify_interface_drop ** OVERDROP: " << n_over
              << " dropped claim(s) violated on the converged solution; the "
                 "redundancy DROP removed a NEEDED continuity constraint and the "
                 "solve converged to a continuity-violating solution. Worst: "
                 "state_id=" << worst_claim.state_id
              << " dom_id="  << worst_claim.dom_id
              << " iel_lo="  << worst_claim.iel_lo
              << " face="    << worst_claim.face
              << " kind="    << worst_claim.kind
              << " (relative violation " << worst_rel << ")." << std::endl;
    return InterfaceDropStatus::OVERDROP;
  }
  return InterfaceDropStatus::OK;
}

//! @brief RAII guard that silences the C-level stdout/stderr file descriptors
//! for its lifetime.  Used only to swallow the benign reference-BLAS xerbla
//! chatter ("** On entry to DTRSV/DGEMV ... had an illegal value") that SuperLU
//! emits when it factorizes a structurally rank-deficient *direct Newton*
//! system.  The resulting step is checked for finiteness and accepted only if
//! the Armijo line search reduces the merit, so suppressing the cosmetic message
//! is safe.  Scoped to the Newton solve alone -- the surrounding solver output
//! and line-search evaluation are unaffected.  No-op on non-POSIX builds.
//!
//! The file descriptors are PROCESS-wide, so the guard is reference-counted under a mutex (2026-10-03): the first
//! guard to enter silences, the last to leave restores.  Concurrent solves (FFGraph::veval threads) used to
//! interleave dup/dup2/close -- one thread "restoring" another's /dev/null -- leaving the whole process's stdout
//! and stderr silenced (the OCFESLV tutorial's veval with MAXTHREAD = 0).  While any thread is inside, every thread's
//! C-level output is silenced; C++ streams are flushed on entry and exit.
struct OCEnvScopedFdSilence
{
#if defined(__unix__) || defined(__APPLE__)
  struct State { std::mutex mtx; size_t depth = 0; int out_fd = -1, err_fd = -1, null_fd = -1; bool active = false; };
  static State& _state(){ static State st; return st; }
  OCEnvScopedFdSilence()
  {
    State& st = _state();
    std::lock_guard<std::mutex> lock( st.mtx );
    if( st.depth++ > 0 ) return;
    std::fflush( stdout ); std::fflush( stderr );
    st.null_fd = ::open( "/dev/null", O_WRONLY );
    if( st.null_fd >= 0 ){
      st.out_fd = ::dup( 1 ); st.err_fd = ::dup( 2 );
      if( st.out_fd >= 0 && st.err_fd >= 0 ){
        ::dup2( st.null_fd, 1 ); ::dup2( st.null_fd, 2 ); st.active = true;
      }
    }
  }
  ~OCEnvScopedFdSilence()
  {
    State& st = _state();
    std::lock_guard<std::mutex> lock( st.mtx );
    if( st.depth == 0 || --st.depth > 0 ) return;
    if( st.active ){
      std::fflush( stdout ); std::fflush( stderr );
      ::dup2( st.out_fd, 1 ); ::dup2( st.err_fd, 2 );
    }
    if( st.out_fd  >= 0 ) ::close( st.out_fd );
    if( st.err_fd  >= 0 ) ::close( st.err_fd );
    if( st.null_fd >= 0 ) ::close( st.null_fd );
    st.out_fd = st.err_fd = st.null_fd = -1;  st.active = false;
  }
#else
  OCEnvScopedFdSilence() {}
#endif
  OCEnvScopedFdSilence( OCEnvScopedFdSilence const& ) = delete;
  OCEnvScopedFdSilence& operator=( OCEnvScopedFdSilence const& ) = delete;
};

//! Internal nonlinear solve (item 12): equilibrated Newton-first descent with Armijo
//! backtracking and sparse damped Levenberg-Marquardt / dense minimum-norm fallbacks.
//! Each iteration assembles the row/column-scaled sparse AD Jacobian Js via deriv() and
//! tries step sources in priority order until the line search accepts one:
//!   (1) a *direct sparse Newton* solve Js dy = -rs through SOLVE.FACTORIZATION (sparse
//!       SuperLU by default).  This forms NO normal equations, so it incurs no JtJ fill-in
//!       and no kappa^2 -- it is what keeps genuinely 3D problems fast.  A direct solve of a
//!       rank-deficient interface Jacobian makes the reference BLAS emit benign chatter,
//!       which is scope-silenced (see OCEnvScopedFdSilence); the step is checked and used
//!       only if it reduces the merit.
//!   (2) a damped normal-equations LM step (Js^T Js + lambda*D) dy = -Js^T rs, formed
//!       lazily only when Newton is unavailable/rejected, so well-conditioned systems never
//!       pay the JtJ cost.  Resolves steps a direct solve cannot take (rank-deficient
//!       interface-augmented Jacobians).
//!   (3) a dense minimum-norm (pseudo-inverse) step via arma::solve's least-squares path,
//!       the last resort for genuinely singular Schur systems (e.g. the degenerate PDE3
//!       IC_STRONG operator).  Once an accepted LM step stalls (barely reduces the merit
//!       while still above tolerance) the solve sticks to this tier for the remainder.
// ======================================================================
//  Captured-reduction fill (shared by marching per-window and monolithic)
// ======================================================================
inline void
OCFESLV::_fill_captures( double const* var, double* inp, double const* cst,
                       int k, std::vector<int> const& cap_window ) const
{
  if( _deferredValue.empty() || !inp ) return;
  auto const edit = ( _evolution_dom_var.dag() ? _mDom.find( _evolution_dom_var ) : _mDom.end() );
  for( size_t li = 0; li < _deferredValue.size(); ++li ){
    t_DeferredValue const& C = _deferredValue[li];
    if( !C.input.dag() || !C.source.dag() ) continue;
    FFVar const cap = C.input, src = C.source;
    if( C.accumulate ){
      if( edit == _mDom.end() ) continue;
      FFDom const&              ed   = edit->second;
      std::vector<double> const wref = ed.lgweights();       // reference LGR weights (sum ~ 2)
      size_t const              nn   = ed.n_node;
      FFVar const               tvar = _evolution_dom_var;
      sample_input( cap,
        [&]( t_Coord const& pt ){
          double win = 0.;
          for( size_t ie = 0; ie < ed.n_elem; ++ie ){
            double const lo = ed.elem_bnd[ie], hi = ed.elem_bnd[ie+1];
            double const jac = 0.5 * ( hi - lo );
            std::vector<double> const tn = ed.lgnodes( lo, hi );
            for( size_t j = 0; j < nn && j < tn.size() && j < wref.size(); ++j ){
              t_Coord q = pt; q[tvar] = tn[j];
              // the node belongs to element [lo,hi]: read THAT element, whatever the rounding of tn[j] -- a lower-half
              // node from above (PLUS), an upper-half node from below (MINUS).  A physical lookup put a left-end node
              // that rounds a hair below lo into the PREVIOUS element: exact for continuous states only to round-off,
              // an O(1) error for a DISCONTINUOUS input (2026-09-28, INPUTS_gate: int x u)
              t_Side const sd{ { tvar, tn[j] < 0.5*( lo + hi )? (int)FFDom::PLUS: (int)FFDom::MINUS } };
              win += jac * wref[j] * eval_colloc<double>( src, q, var, inp, cst, &sd );
            }
          }
          // k<0 (monolithic): single OVERWRITING full-horizon fill (prior=0).
          // k>=0 (marching):  read-modify-write ACCUMULATE of this window's sub-integral.
          double const prior = ( k < 0 ) ? 0. : eval_colloc<double>( cap, pt, var, inp, cst );
          return prior + win;
        }, inp );
    }
    else{   // LATCH
      bool const due = ( k < 0 ) || ( li < cap_window.size() && cap_window[li] == k );
      if( !due || !_evolution_dom_var.dag() ) continue;
      FFVar  const tvar = _evolution_dom_var;
      double const tau  = C.tau;
      sample_input( cap,
        [&]( t_Coord const& pt ){ t_Coord q = pt; q[tvar] = tau;  t_Side const sd{ { tvar, C.side } };
                                  return eval_colloc<double>( src, q, var, inp, cst, &sd ); },
        inp );
    }
  }
}

// Capture-SENSITIVITY fill (Stage 1): same quadrature as _fill_captures, on the state-sensitivity
// buffer s_var, writing d(cap)/dp into dcap.  src is a bare (materialised) state, so eval_colloc(src,
// s_var, dcap) = s_var[src] = d(src)/dp; eval_colloc(cap, s_var, dcap) = dcap[cap] (cap is an input,
// so it reads the dcap buffer -- the running d(cap)/dp for read-modify-write).
inline void
OCFESLV::_fill_capture_sens( double const* s_var, int k,
                           std::vector<int> const& cap_window, double* dcap ) const
{
  if( _deferredValue.empty() || !dcap ) return;
  auto const edit = ( _evolution_dom_var.dag() ? _mDom.find( _evolution_dom_var ) : _mDom.end() );
  for( size_t li = 0; li < _deferredValue.size(); ++li ){
    t_DeferredValue const& C = _deferredValue[li];
    if( !C.input.dag() || !C.source.dag() ) continue;
    FFVar const cap = C.input, src = C.source;
    if( C.accumulate ){
      if( edit == _mDom.end() ) continue;
      FFDom const&              ed   = edit->second;
      std::vector<double> const wref = ed.lgweights();
      size_t const              nn   = ed.n_node;
      FFVar const               tvar = _evolution_dom_var;
      sample_input( cap,
        [&]( t_Coord const& pt ){
          double win = 0.;
          for( size_t ie = 0; ie < ed.n_elem; ++ie ){
            double const lo = ed.elem_bnd[ie], hi = ed.elem_bnd[ie+1];
            double const jac = 0.5 * ( hi - lo );
            std::vector<double> const tn = ed.lgnodes( lo, hi );
            for( size_t j = 0; j < nn && j < tn.size() && j < wref.size(); ++j ){
              t_Coord q = pt; q[tvar] = tn[j];
              t_Side const sd{ { tvar, tn[j] < 0.5*( lo + hi )? (int)FFDom::PLUS: (int)FFDom::MINUS } };   // this node's own element
              win += jac * wref[j] * eval_colloc<double>( src, q, s_var, dcap, nullptr, &sd );   // d(src)/dp
            }
          }
          double const prior = ( k < 0 ) ? 0. : eval_colloc<double>( cap, pt, s_var, dcap, nullptr );
          return prior + win;
        }, dcap );
    }
    else{   // LATCH
      bool const due = ( k < 0 ) || ( li < cap_window.size() && cap_window[li] == k );
      if( !due || !_evolution_dom_var.dag() ) continue;
      FFVar  const tvar = _evolution_dom_var;
      double const tau  = C.tau;
      sample_input( cap,
        [&]( t_Coord const& pt ){ t_Coord q = pt; q[tvar] = tau;  t_Side const sd{ { tvar, C.side } };
                                  return eval_colloc<double>( src, q, s_var, dcap, nullptr, &sd ); },
        dcap );
    }
  }
}


inline void
OCFESLV::_eval_functions_captured( double const* var, double const* inp, double const* cst )
{
  size_t const nf = n_colloc_fct();
  if( !nf ){ _valFct.clear(); return; }
  size_t const ni = n_colloc_inp();
  std::vector<double> eqndum( n_colloc_eqn(), 0. ), fctv( nf, 0. );
  double const* useinp = inp;
  std::vector<double> capbuf;
  if( ni && !_deferredValue.empty() ){
    capbuf.assign( ni, 0. );
    if( inp ) std::copy( inp, inp + ni, capbuf.begin() );
    _fill_captures( var, capbuf.data(), cst, -1, std::vector<int>() );   // monolithic: single full-horizon fill
    useinp = capbuf.data();
  }
  if( eval( eqndum.data(), fctv.data(), var, useinp, cst ) ) _valFct = std::move( fctv );
  else _valFct.clear();
}

//! Mode-agnostic: operates on the setup() var[] layout; never branches on the imposition mode.
inline OCFESLV::SolveReport
OCFESLV::_solve_marching( double* var, double const* inp, double const* cst )
{
  SolveReport rep;
  size_t const N = n_march_steps();

  // ---------------------------------------------------------------------
  // rev93: UNCONDITIONAL-under-EVOFACE marching stamp.
  // ---------------------------------------------------------------------
  // The only previous evidence that a model marched was the "[march] step" line ~150
  // lines below, which is gated on options.SOLVE.VERBOSE.  That made `grep -c "march]
  // step"` a measurement of SOLVE.VERBOSE, not of marching -- and handoff rev25 section
  // 0.3 used exactly that grep to conclude "Marching is not involved ... so the failure is
  // not a solve-path difference", on OCFE_MMPDE27, which sets SOLVE.VERBOSE=false unless
  // --verbose is passed (OCFE_MMPDE27.cpp:868).  The count was 0 because the print was
  // off.  That conclusion is RETRACTED; the arc's own rule -- a diagnostic gated on
  // verbosity is not a diagnostic -- applies to it.
  //
  // This stamp is gated only on the audit flag, and reports BOTH outcomes, so a zero is
  // distinguishable from a silence.
  {
    static bool const kEvoFaceM = Options::_env_flag( "CRONOS_AUDIT_EVOFACE", false );
    if( kEvoFaceM ){
      if( !N || _marchGrid.size() < 2 )
        _disp( 3 ) << "  [evoface] march: NOT MARCHED (n_march_steps=" << N
                  << ", grid=" << _marchGrid.size()
                  << ") -- monolithic block solve\n";
      else
        _disp( 3 ) << "  [evoface] march: MARCHED " << N << " window(s) over ["
                  << _marchGrid.front() << "," << _marchGrid.back() << "]\n";
    }
  }

  if( !N || _marchGrid.size() < 2 ) return rep;

  // Captured evolution-direction reductions.  LATCH (accumulate=false): static schedule
  // window(tau) = first window whose UB >= tau (tau in (t_k,t_{k+1}] -> window k; a boundary tau at an
  // interior grid node tie-broken to the LEFT window).  ACCUM (accumulate=true): running sum over
  // windows, no schedule.  Both written POST-SOLVE (below).
  std::vector<int> cap_window( _deferredValue.size(), -1 );
  {
    double const span = std::fabs( _marchGrid.back() - _marchGrid.front() );
    double const eps  = 1e-9 * std::max( 1.0, span );
    for( size_t li = 0; li < _deferredValue.size(); ++li )
      if( !_deferredValue[li].accumulate )
        for( size_t kk = 0; kk < N; ++kk )
          if( _deferredValue[li].side == FFDom::PLUS
              ? ( _deferredValue[li].tau < _marchGrid[kk+1] - eps || kk+2 == _marchGrid.size() )   // tau^+: the window STARTING at tau (at tf: the last)
              : _deferredValue[li].tau <= _marchGrid[kk+1] + eps ){  // tau^-: the window ENDING at tau
            cap_window[li] = (int)kk; break; }
  }

  // The IC input is (re)written per element, so march on a mutable copy of inp[].
  size_t const ni = n_colloc_inp();
  std::vector<double> inpwork;
  if( ni ){
    inpwork.assign( ni, 0. );
    if( inp ) std::copy( inp, inp + ni, inpwork.begin() );
  }
  double* pinp = ni ? inpwork.data() : nullptr;

  _inMarch = true;
  bool all_conv = true;
  bool const broadcast = ( options.SOLVE.WARMSTART == Options::BROADCAST_IC );
  std::vector<std::vector<double>> prevterm( _marchTransfers.size() );  // per-pair terminal (broadcast)

  // Output-function accumulation.  eval(fct) on each converged block gives that element's
  // contribution -- the collapsed window is exactly [t_k, t_{k+1}], so an evolution-
  // direction integral (OpI over the evolution var) yields the sub-interval integral,
  // which we SUM; every other output (terminal / final-time / spatial-only) is OVERWRITTEN
  // so it retains the LAST element's value.
  size_t const nfct = n_colloc_fct();
  size_t const neqn = n_colloc_eqn();
  std::vector<double> fctacc( nfct, 0. ), fctblk( nfct, 0. ), eqndum( neqn, 0. );
  std::vector<char>   fctsum( nfct, 0 );   // 1 = accumulate (evolution-direction integral)
  for( auto const& f : _mFct ){
    bool const ei = _fct_is_evolution_integral( f.var );
    for( size_t r = f.row0; r < f.row0 + f.nrow && r < nfct; ++r ) fctsum[r] = ei ? 1 : 0;
  }
  if( options.OUTPUT.MARCH_STORE ) _marchTrajectory.clear();
  size_t const nsta = n_colloc_sta();
  for( size_t k = 0; k < N; ++k ){
    double const t0 = _marchGrid[k], t1 = _marchGrid[k+1], w = t1 - t0;
    _marchWindow = k;         // resident window: distributed-input evolution element remaps to k

    _marchICActive = false;   // general-IC override is live only at k>0 for override-path states

    // Between elements: while the evolution UB is still the previous element's terminal
    // (before sliding the window), transfer every state's terminal into its IC input, and
    // (for BROADCAST_IC) capture each terminal profile for the warm-start.
    if( k > 0 ){
      if( !_marchICRowsReady ) _build_march_ic_rows();   // provenance is ready after window 0
      // terminal_profile() is the per-spatial-node UB profile.  It is needed for the BROADCAST_IC
      // warm-start AND as the continuity data for every general-IC override state (one terminal per
      // INITIAL row); compute it once when either consumer is present.
      bool const need_prof = broadcast || !_marchICRows.empty();
      auto prof = ( need_prof ? terminal_profile( var, pinp, cst )
                              : std::map<FFVar,std::vector<double>,lt_FFVar>{} );
      for( size_t j = 0; j < _marchTransfers.size(); ++j ){
        FFVar const& st  = _marchTransfers[j].first;
        FFVar const& inp = _marchTransfers[j].second;
        if( !st.dag() ) continue;
        if( inp.dag() ){                                 // canonical: data-driven ic_input write
          if( broadcast ){ auto it = prof.find( st ); prevterm[j] = ( it != prof.end() ) ? it->second : std::vector<double>{}; }
          transfer_terminal( st, inp, var, pinp, cst );
        }
        else {                                           // general-IC override: capture terminal for the row swap
          auto itr = _marchICRows.find( st );
          if( itr != _marchICRows.end() && !itr->second.empty() ){
            auto itp = prof.find( st );
            auto itc = _marchICColToNode.find( st );
            size_t const nnode = ( itc != _marchICColToNode.end() ) ? itc->second.size() : 0;
            // The terminal profile must cover every LB-face NODE, not every override ROW.  An
            // INITIAL equation may legitimately be imposed on a SUBSET of the face -- a moving-mesh
            // IC on the xi-INTERIOR, with the two Dirichlet ends carrying the remaining rows -- so
            // rows.size() < nodes.  That is fine: ic_term_for_col indexes the profile by each row's
            // OWN column -> node, so a subset of rows still reads the right terminals.  Requiring
            // rows == profile sent those cases to the scalar-broadcast fallback, and
            // terminal_value(st,{}) on a DISTRIBUTED state has no coordinate to evaluate at, so it
            // threw OCBase::Exceptions -- an abort several layers from the actual mismatch.
            if( itp != prof.end() && nnode && itp->second.size() == nnode )
              _marchICTerminal[st] = itp->second;      // node-indexed profile (rows may be a subset)
            else if( itp != prof.end() && itp->second.size() == itr->second.size() )
              _marchICTerminal[st] = itp->second;      // exact row/profile match (legacy path)
            else {                                     // last resort: broadcast a single UB value
              try {
                double const term = terminal_value( st, t_Coord{}, var, pinp, cst );
                _marchICTerminal[st].assign( itr->second.size(), term );
              }
              catch( ... ){
                std::cerr << "OCFESLV::solve ** marching IC transfer: cannot form a terminal for "
                          << st << " -- override rows=" << itr->second.size()
                          << " LB-face nodes=" << nnode << " profile="
                          << ( itp != prof.end() ? itp->second.size() : 0 )
                          << "; a distributed state has no scalar terminal.\n";
                throw;
              }
            }
            _marchICActive = true;
          }
        }
      }
    }

    if( k > 0 && !_apply_march_jumps( k, var, pinp, cst ) ) throw Exceptions( Exceptions::INTERNAL );   // transitions at this seam

    reset_evolution_window( t0, w );

    // Item 2 -- per-element input re-seeding: set every input that has an update_ref
    // reference from that reference at the CURRENT window's absolute nodes (time-varying
    // feed/BC).  IC-transfer targets are set by the transfer above and skipped; inputs
    // without a reference keep their values.  For a spatial-only input the reference is
    // t-independent, so this is a harmless no-op; for a t-varying feed it advances it.
    for( auto const& kv : _classInpRefFun ){
      FFVar const& inpvar = kv.first;
      if( !kv.second ) continue;
      bool is_ic = false;
      for( auto const& pr : _marchTransfers )
        if( pr.second.dag() && pr.second.id().second == inpvar.id().second ){ is_ic = true; break; }
      // A REGISTERED CONTROL's DOFs are decision variables carried in the input buffer (init seeds them
      // for every window via the same fun + interior-side correction sample_input would apply, so the
      // base solve is byte-identical).  Re-sampling would clobber a perturbed/optimiser-set control with
      // its own reference, hiding it from the window solve -- so skip controls here; genuine time-varying
      // feeds/BCs (not registered) still advance to this window's absolute nodes.
      if( !is_ic && !_is_registered_control( inpvar ) ) sample_input( inpvar, kv.second, pinp );
    }

    if( k > 0 && broadcast )
      for( size_t j = 0; j < _marchTransfers.size(); ++j )
        if( _marchTransfers[j].first.dag() && !prevterm[j].empty() )
          warmstart_broadcast( _marchTransfers[j].first, prevterm[j], var );
    // REUSE (default) and EXTRAPOLATE (-> reuse for now) keep the previous window's var[].

    SolveReport const r = solve( var, pinp, cst );   // monolithic block solve (guarded by _inMarch)
    all_conv = all_conv && r.converged;
    if( k == 0 ) rep.initial_residual = r.initial_residual;
    rep.iterations   += r.iterations;
    rep.lm_steps     += r.lm_steps;
    rep.newton_steps += r.newton_steps;
    rep.final_residual = r.final_residual;
    rep.final_merit    = r.final_merit;
    if( options.SOLVE.VERBOSE )
      _disp( 1 ) << "OCFESLV::solve [march] step " << k << " window ["
                << t0 << "," << t1 << "] converged=" << ( r.converged ? "y" : "n" )
                << " max|r|=" << r.final_residual << "\n";

    // Captured reductions (post-solve, from the just-solved window): ACCUM sums this window's
    // sub-integral into cap; LATCH writes cap=R(tau) at tau's window.  Shared with the monolithic
    // finalize (which calls the same helper once over the full horizon, k=-1).
    _fill_captures( var, pinp, cst, (int)k, cap_window );

    // Accumulate this element's output-function contribution (domain is now [t_k,t_{k+1}]).
    // Evolution-direction integrals SUM.  A point/terminal output is OVERWRITTEN only in the window
    // whose evolution range [t0,t1] contains its evolution-t point (so an initial-time output like
    // c(0) is captured at window 0, a final-time output at the last window); an output with no
    // evolution-t coordinate (spatial-only, t-independent) keeps the last-window value as before.
    if( nfct ){
      eval( eqndum.data(), fctblk.data(), var, pinp, cst );
      for( auto const& f : _mFct ){
        bool const summed = ( f.row0 < nfct && fctsum[f.row0] );
        // An output DISTRIBUTED over the evolution direction (2026-10-06, WORKPLAN 1.1): this window emitted the rows
        // of ITS time element, in its emission order; place them in the monolithic layout
        if( !summed && _march_field_over_evo( f ) ){
          std::vector<size_t> const rows = _march_field_rows( f, k );
          for( size_t i = 0; i < rows.size() && f.row0 + i < nfct; ++i )
            if( rows[i] != npos && f.row0 + rows[i] < nfct ) fctacc[ f.row0 + rows[i] ] = fctblk[ f.row0 + i ];
          continue;
        }
        // Value capture routed through the SAME predicate as the sensitivity capture (UB-owned), so a
        // point output lands in exactly one window on both paths.  For continuous states LB and UB agree
        // by IC continuity (unchanged), so this only re-sides a DISCONTINUOUS seam output to its left-
        // limit -- matching where its gradient is now formed.  _march_point_capture returns true every
        // window for a no-evolution-point output, which last-wins to the final window as before.
        bool const capture = summed || _march_point_capture( f, t0, t1 );
        for( size_t r = f.row0; r < f.row0 + f.nrow && r < nfct; ++r ){
          if( summed )       fctacc[r] += fctblk[r];
          else if( capture ) fctacc[r]  = fctblk[r];
        }
      }
    }

    // Optionally capture this element's block solution for plotting/inspection.
    if( options.OUTPUT.MARCH_STORE ){
      MarchBlock blk;
      blk.t0 = t0; blk.t1 = t1;
      blk.var.assign( var, var + nsta );
      if( pinp ) blk.inp.assign( pinp, pinp + ni );
      _marchTrajectory.push_back( std::move( blk ) );
    }
  }
  _inMarch = false;
  _marchICActive = false;   // override is a marching-only, per-window state
  _valFct = std::move( fctacc );
  _lastSolveMarched = true;   // last solve marched; eval_colloc routing gates additionally on a stored trajectory
  // Approach A: cache the accumulated caps + last-window var so a direct eval() of a captured output can
  // serve them (valid only at this solved point; scope-keyed on var[] in eval()).
  _marchCapValid = ( all_conv && !_deferredValue.empty() && ni > 0 );
  if( _marchCapValid ){
    _marchCapInp.assign( pinp, pinp + ni );
    _marchCapVar.assign( var, var + n_colloc_sta() );
  }
  else { _marchCapInp.clear(); _marchCapVar.clear(); }
  rep.converged = all_conv;

  // =====================================================================
  // rev97 BLOCKER-1 INSTRUMENT -- residual on the resident window, in place
  // =====================================================================
  // OPEN QUESTION (handoff rev28 s4.1, rev30 s1).  OCFE_PDE9's post-solve
  // eval_residual() read max|r| = 4.906440e-07 after a marched solve, while BOTH windows
  // had converged (4.651e-13 and 2.258e-12) and the collapsed system it evaluates is only
  // n=12000 rather than 24000.  So it is neither a stall nor conditioning -- but why the
  // COLLAPSED system reads 4.9e-07 instead of ~1e-12 is unexplained, and max/mean = 192
  // says it sits on a handful of rows.
  //
  // WHY IT MATTERS.  If a marched solution genuinely leaves a large residual on the
  // assembled system, then EVERY driver that checks one after solving sees
  // discretisation-error magnitude rather than 1e-14.  OCFE_PDE9 was patched to gate on
  // srep.converged instead; if the effect is real rather than an artefact of evaluating a
  // collapsed environment AFTER the fact, other drivers need the same treatment and we do
  // not know how many.
  //
  // THE DISCRIMINATOR.  Evaluate the residual HERE -- inside _solve_marching, on the state
  // the final window just converged, before anything is reset or returned -- and compare it
  // with what the driver measures afterwards.
  //
  //   the two AGREE   -> the collapsed system really does carry that residual.  A marched
  //                      solution is a solution of a DIFFERENT system, the caveat is real,
  //                      and every residual-gating driver needs PDE9's treatment.
  //   the two DIFFER  -> something changes between this point and the driver's call --
  //                      window position, a stale transferred IC input, a reset -- and the
  //                      caveat is an artefact of the post-hoc read, not of marching.
  //
  // Two numbers settle a question that four revisions of reasoning did not.
  {
    static bool const kMarchRes = Options::_env_flag( "CRONOS_AUDIT_EVOFACE", false );
    if( kMarchRes ){
      size_t const nr = n_colloc_eqn();
      std::vector<double> rr( nr, 0. );
      if( nr && eval( rr.data(), nullptr, var, inp, cst ) ){
        double mx = 0., sm = 0.;
        size_t arg = 0;
        for( size_t i = 0; i < nr; ++i ){
          double const a = std::fabs( rr[i] );
          sm += a;
          if( a > mx ){ mx = a; arg = i; }
        }
        double const mean = sm / double( nr );
        _disp( 3 ) << "  [evoface] march: residual on the RESIDENT window, evaluated IN PLACE"
                  << " immediately after the final step:  n=" << nr
                  << "  max|r|=" << mx << " at row " << arg
                  << "  mean|r|=" << mean
                  << "  max/mean=" << ( mean > 0. ? mx / mean : 0. ) << "\n"
                  << "  [evoface]   compare against the driver's own post-solve residual."
                     "  EQUAL -> the collapsed system carries it and the caveat is real;"
                     "  DIFFERENT -> something changes after this point and it is an"
                     " artefact of the post-hoc read.\n";
      }
      else
        _disp( 3 ) << "  [evoface] march: in-place residual unavailable"
                     " (n_colloc_eqn=" << nr << ", eval failed)\n";
    }
  }

  // rev43: duplicate-node spreads for MARCHED solves.
  //
  // rev42 skipped the report whenever _inMarch was set and left a comment saying the
  // meaningful frame is "the assembled trajectory the marching loop reports" -- without
  // adding a report here.  The result was that marching drivers were SILENTLY uninstrumented:
  // six of the 28 flagged drivers (PDE30, PDE30b, PDE31b, PDE33, PDE35, PDE38) produced no
  // spread row at all, and the sweep summary gave no indication they were missing.
  //
  // The per-window report really is the wrong frame -- it would emit one block per window and
  // measure continuity inside a window rather than across the trajectory -- but "wrong frame"
  // is not a reason to emit nothing.  This fires once, after the last window, on the assembled
  // var[], which is the same object a monolithic solve reports on.
  //
  // This matters more than the six rows: roadmap C3 makes moving mesh UNDER MARCHING the
  // highest-value item in Part C, precisely because collapsing the evolution domain removes
  // interior t-seams and hence the corner claims this whole arc was about.  A gate that goes
  // dark on exactly the configuration the framework is heading toward is worse than no gate.
  if( _knob_AUDIT_DUP_SPREAD() )
    _report_dup_spread_gated( var, all_conv
        ? "post solve [march, assembled]  conv=yes"
        : "post solve [march, assembled]  conv=NO -- STRUCTURAL ONLY, NOT AN ACCURACY NUMBER" );
  // rev63: the assembled trajectory is NOT the monolithic [states | tau] layout the beta
  // column map was built against, so the tau columns do not mean there what they mean in a
  // monolithic solve.  rev62 tried to guard this inside the reporter on `_inMarch`, which is
  // already false by the time this line runs -- so the guard was dead and the number was
  // printed regardless.  Refuse here, where the caller knows what it is holding.
  if( _knob_AUDIT_BETA() )
    _disp( 3 ) << "  [beta] not reported for the assembled march trajectory: its layout is not"
                 " the monolithic [states|tau] var[] this measurement indexes." << std::endl;

  return rep;
}

inline double
OCFESLV::eval_solution
( FFVar const& state, t_Coord const& pt, t_Side const* side ) const
{
  if( !_issetup ) throw Exceptions( Exceptions::SETUP );

  // Marching solve: read from the stored trajectory.  eval_colloc() routes an evolution-dependent
  // state to the block whose window contains pt[evolution]; the passed buffer is a non-null
  // fallback used only for a non-evolution state (which reads the last window).
  if( _lastSolveMarched ){
    if( _marchTrajectory.empty() )
      throw Exceptions( Exceptions::NOSTORE );   // marched with OUTPUT.MARCH_STORE disabled
    auto const& b = _marchTrajectory.back();
    return eval_colloc<double>( state, pt, b.var.data(),
                                b.inp.empty() ? nullptr : b.inp.data(), nullptr, side );
  }

  // Monolithic solve: read from the cached primal solution (SOLVE.REUSE).
  if( _primalValid && _primalVar.size() == n_colloc_sta() )
    return eval_colloc<double>( state, pt, _primalVar.data(),
                                _primalInp.empty() ? nullptr : _primalInp.data(), nullptr, side );

  throw Exceptions( Exceptions::NOSTORE );   // no cached solution: enable SOLVE.REUSE, or use eval_colloc() with a buffer
}


inline bool
OCFESLV::_march_fsens_setup
( double* var, double const* inp, double const* cst )
{
  _fsBlocks.clear();
  size_t const N = n_march_steps();
  if( !N || _marchGrid.size() < 2 ) return false;
  size_t const ni  = n_colloc_inp();
  size_t const nv  = _nCollVar;
  size_t const neq = n_colloc_eqn();
  size_t const nf  = n_colloc_fct();

  std::vector<double> inpwork( ni, 0. );
  if( inp ) std::copy( inp, inp + ni, inpwork.begin() );
  double* pinp = ni ? inpwork.data() : nullptr;

  bool const broadcast = ( options.SOLVE.WARMSTART == Options::BROADCAST_IC );
  std::vector<std::vector<double>> prevterm( _marchTransfers.size() );

  // Per-window function accumulation (fixes _sensFct so a complete solve_fsens()/solve_asens()
  // reports fresh functions without a separate primal march): evolution-integral rows sum over
  // windows, others take the last window -- mirrors _solve_marching's _valFct accumulation.
  std::vector<char>   fctsum( nf, 0 );
  for( auto const& f : _mFct ){ bool const ei = _fct_is_evolution_integral( f.var ); for( size_t r = f.row0; r < f.row0 + f.nrow && r < nf; ++r ) fctsum[r] = ei ? 1 : 0; }
  std::vector<double> fctacc( nf, 0. );

  // LATCH schedule (as in _solve_marching) so captures fill correctly during the setup march; without
  // this cap[] stays 0 and deriv() would compute df/dcap at cap=0 (wrong for a nonlinear f(cap)).
  std::vector<int> cap_window( _deferredValue.size(), -1 );
  if( !_deferredValue.empty() && _marchGrid.size() >= 2 ){
    double const span = std::fabs( _marchGrid.back() - _marchGrid.front() );
    double const eps  = 1e-9 * std::max( 1.0, span );
    for( size_t li = 0; li < _deferredValue.size(); ++li )
      if( !_deferredValue[li].accumulate )
        for( size_t kk = 0; kk + 1 < _marchGrid.size(); ++kk )
          if( _deferredValue[li].side == FFDom::PLUS
              ? ( _deferredValue[li].tau < _marchGrid[kk+1] - eps || kk+2 == _marchGrid.size() )   // tau^+: the window STARTING at tau (at tf: the last)
              : _deferredValue[li].tau <= _marchGrid[kk+1] + eps ){  // tau^-: the window ENDING at tau
            cap_window[li] = (int)kk; break; }
  }
  _inMarch = true;
  bool all_conv = true;
  for( size_t k = 0; k < N; ++k ){
    double const t0 = _marchGrid[k], t1 = _marchGrid[k+1], w = t1 - t0;
    _marchWindow = k;         // resident window: distributed-input evolution element remaps to k
    _marchICActive = false;   // general-IC override live only at k>0 override states (mirrors _solve_marching)
    if( k > 0 ){
      if( !_marchICRowsReady ) _build_march_ic_rows();
      // Per-node UB profile is needed for BROADCAST_IC warm-start AND for the general-IC override
      // continuity swap -- capture it once, exactly as _solve_marching does (terminal_value() with an
      // empty spatial coord would throw on a DISTRIBUTED override state).
      bool const need_prof = broadcast || !_marchICRows.empty();
      auto prof = ( need_prof ? terminal_profile( var, pinp, cst )
                              : std::map<FFVar,std::vector<double>,lt_FFVar>{} );
      for( size_t j = 0; j < _marchTransfers.size(); ++j ){
        FFVar const& st = _marchTransfers[j].first;
        FFVar const& in = _marchTransfers[j].second;
        if( !st.dag() ) continue;
        if( in.dag() ){
          if( broadcast ){ auto it = prof.find( st ); prevterm[j] = ( it != prof.end() ) ? it->second : std::vector<double>{}; }
          transfer_terminal( st, in, var, pinp, cst );
        }
        else {   // general-IC override: capture this window's per-node terminal profile for the continuity swap
          auto itr = _marchICRows.find( st );
          if( itr != _marchICRows.end() && !itr->second.empty() ){
            auto itp = prof.find( st );
            auto itc = _marchICColToNode.find( st );
            size_t const nnode = ( itc != _marchICColToNode.end() ) ? itc->second.size() : 0;
            // The terminal profile must cover every LB-face NODE, not every override ROW.  An
            // INITIAL equation may legitimately be imposed on a SUBSET of the face -- a moving-mesh
            // IC on the xi-INTERIOR, with the two Dirichlet ends carrying the remaining rows -- so
            // rows.size() < nodes.  That is fine: ic_term_for_col indexes the profile by each row's
            // OWN column -> node, so a subset of rows still reads the right terminals.  Requiring
            // rows == profile sent those cases to the scalar-broadcast fallback, and
            // terminal_value(st,{}) on a DISTRIBUTED state has no coordinate to evaluate at, so it
            // threw OCBase::Exceptions -- an abort several layers from the actual mismatch.
            if( itp != prof.end() && nnode && itp->second.size() == nnode )
              _marchICTerminal[st] = itp->second;      // node-indexed profile (rows may be a subset)
            else if( itp != prof.end() && itp->second.size() == itr->second.size() )
              _marchICTerminal[st] = itp->second;      // exact row/profile match (legacy path)
            else {                                     // last resort: broadcast a single UB value
              try {
                double const term = terminal_value( st, t_Coord{}, var, pinp, cst );
                _marchICTerminal[st].assign( itr->second.size(), term );
              }
              catch( ... ){
                std::cerr << "OCFESLV::solve ** marching IC transfer: cannot form a terminal for "
                          << st << " -- override rows=" << itr->second.size()
                          << " LB-face nodes=" << nnode << " profile="
                          << ( itp != prof.end() ? itp->second.size() : 0 )
                          << "; a distributed state has no scalar terminal.\n";
                throw;
              }
            }
            _marchICActive = true;
          }
        }
      }
    }
    if( k > 0 && !_apply_march_jumps( k, var, pinp, cst ) ) throw Exceptions( Exceptions::INTERNAL );   // transitions at this seam
    reset_evolution_window( t0, w );
    for( auto const& kv : _classInpRefFun ){
      FFVar const& iv = kv.first;
      if( !kv.second ) continue;
      bool is_ic = false;
      for( auto const& pr : _marchTransfers ) if( pr.second.dag() && pr.second.id().second == iv.id().second ){ is_ic = true; break; }
      if( !is_ic && !_is_registered_control( iv ) ) sample_input( iv, kv.second, pinp );   // controls live in the buffer -- see _solve_marching
    }
    if( k > 0 && broadcast )
      for( size_t j = 0; j < _marchTransfers.size(); ++j )
        if( _marchTransfers[j].first.dag() && !prevterm[j].empty() )
          warmstart_broadcast( _marchTransfers[j].first, prevterm[j], var );

    SolveReport const r = solve( var, pinp, cst );
    all_conv = all_conv && r.converged;
    if( !r.converged ){ _inMarch = false; _fsBlocks.clear(); return false; }

    // Fill captured-reduction inputs from THIS window's solution (accumulate into cap[]), so the
    // stored block linearization below differentiates f(cap) at the correct cap value.
    _fill_captures( var, pinp, cst, (int)k, cap_window );

    // linearize at the converged block, factorize J_k, and store for reuse across seeds.
    FSensBlock B;
    B.t0 = t0; B.t1 = t1;
    B.var.assign( var, var + nv );
    if( pinp ) B.inp.assign( pinp, pinp + ni );
    std::vector<size_t> nnzeqn( neq, 0 ), nnzfct( nf, 0 );
    if( !deriv( nullptr, nullptr, nnzeqn.data(), nullptr, nullptr, nullptr, nnzfct.data(), nullptr, var, pinp, cst ) ){ _inMarch = false; _fsBlocks.clear(); return false; }
    B.geqn.resize( neq ); B.ceqn.resize( neq ); B.gfct.resize( nf ); B.cfct.resize( nf );
    std::vector<double*> geqnp( neq, nullptr ), gfctp( nf, nullptr );
    std::vector<size_t*> ceqnp( neq, nullptr ), cfctp( nf, nullptr );
    for( size_t r = 0; r < neq; ++r ){ B.geqn[r].resize( nnzeqn[r] ); B.ceqn[r].resize( nnzeqn[r] ); geqnp[r] = B.geqn[r].data(); ceqnp[r] = B.ceqn[r].data(); }
    for( size_t r = 0; r < nf;  ++r ){ B.gfct[r].resize( nnzfct[r] ); B.cfct[r].resize( nnzfct[r] ); gfctp[r] = B.gfct[r].data(); cfctp[r] = B.cfct[r].data(); }
    std::vector<double> eqnbuf( neq, 0. ), fctbuf( nf, 0. );
    if( !deriv( eqnbuf.data(), geqnp.data(), nnzeqn.data(), ceqnp.data(), fctbuf.data(), gfctp.data(), nnzfct.data(), cfctp.data(), var, pinp, cst ) ){ _inMarch = false; _fsBlocks.clear(); return false; }
    // Match the primal's general-IC override in the linearised block: at k>0 the INITIAL rows are the
    // continuity rows var[lb_col] - terminal, so their Jacobian is the unit row on lb_col -- not the
    // general-IC row deriv() just produced.  (The RHS terminal coupling is applied in *_apply.)
    if( _marchICActive && !_marchICRows.empty() )
      for( auto const& [st, rows] : _marchICRows )
        for( size_t i = 0; i < rows.size(); ++i ){
          auto const& rr = rows[i];
          if( rr.row >= neq ) continue;
          for( size_t idx = 0; idx < B.geqn[rr.row].size(); ++idx )
            B.geqn[rr.row][idx] = ( B.ceqn[rr.row][idx] == rr.lb_col ) ? 1.0 : 0.0;
          B.icRows.emplace_back( rr.row, st, i );
        }
    if( !_sens_factorise( B, neq, nv, nnzeqn, "march" ) ){
      _inMarch = false; _fsBlocks.clear(); return false; }
    for( auto const& f : _mFct ){
      bool const summed = ( f.row0 < nf && fctsum[f.row0] );
      bool const cap = summed || _march_point_capture( f, t0, t1 );
      for( size_t r = f.row0; r < f.row0 + f.nrow && r < nf; ++r ){
        if( summed )   fctacc[r] += fctbuf[r];
        else if( cap ) fctacc[r]  = fctbuf[r];
      }
    }
    _fsBlocks.push_back( std::move( B ) );
  }
  _inMarch = false;
  _marchICActive = false;
  _sensFct = std::move( fctacc );
  return all_conv;
}

inline bool
OCFESLV::_march_fsens_apply
( std::vector<double> const& dinp_dp0, std::vector<double>& dFdp )
{
  size_t const ni = n_colloc_inp(), nv = _nCollVar, nf = n_colloc_fct();
  dFdp.assign( nf, 0. );
  if( _fsBlocks.empty() ) return false;
  std::vector<char> fctsum( nf, 0 );
  for( auto const& f : _mFct ){ bool const ei = _fct_is_evolution_integral( f.var ); for( size_t r = f.row0; r < f.row0 + f.nrow && r < nf; ++r ) fctsum[r] = ei ? 1 : 0; }
  std::map<FFVar,std::vector<double>,lt_FFVar> sens_term;

  // Stage 1: persistent d(cap)/dp accumulator (input-array layout; only cap slots used) and the LATCH
  // schedule window(tau), recomputed exactly as _solve_marching does.
  std::vector<double> dcap_acc( ni, 0. );
  std::vector<int>    cap_window( _deferredValue.size(), -1 );
  if( !_deferredValue.empty() && _marchGrid.size() >= 2 ){
    double const span = std::fabs( _marchGrid.back() - _marchGrid.front() );
    double const eps  = 1e-9 * std::max( 1.0, span );
    for( size_t li = 0; li < _deferredValue.size(); ++li )
      if( !_deferredValue[li].accumulate )
        for( size_t kk = 0; kk + 1 < _marchGrid.size(); ++kk )
          if( _deferredValue[li].side == FFDom::PLUS
              ? ( _deferredValue[li].tau < _marchGrid[kk+1] - eps || kk+2 == _marchGrid.size() )   // tau^+: the window STARTING at tau (at tf: the last)
              : _deferredValue[li].tau <= _marchGrid[kk+1] + eps ){  // tau^-: the window ENDING at tau
            cap_window[li] = (int)kk; break; }
  }
  for( size_t k = 0; k < _fsBlocks.size(); ++k ){
    FSensBlock const& B = _fsBlocks[k];
    // block seed: element-0 control direction, with canonical IC-inputs overwritten by terminal(s_{k-1}).
    std::vector<double> dinp_dp = dinp_dp0;
    if( k > 0 && !_march_jump_fsens( k, sens_term, dinp_dp0 ) ) return false;   // through a transition at this seam
    if( k > 0 )
      for( auto const& pr : _marchTransfers ){ if( !pr.second.dag() ) continue; auto it = sens_term.find( pr.first ); if( it != sens_term.end() && !it->second.empty() ) set_input_values( pr.second, it->second, dinp_dp.data() ); }
    // RHS b = -(dR/dinp . dinp_dp)
    arma::vec b( B.geqn.size(), arma::fill::zeros );
    for( size_t r = 0; r < B.geqn.size(); ++r )
      for( size_t idx = 0; idx < B.geqn[r].size(); ++idx ){ size_t c = B.ceqn[r][idx]; if( c >= nv && c - nv < ni ) b( r ) -= B.geqn[r][idx] * dinp_dp[c-nv]; }
    // General-IC override rows: the continuity row var[lb_col] - terminal_{k-1} has unit Jacobian on
    // lb_col, so s(lb_col) = d(terminal_{k-1})/dp.  There is no ic_input to carry it (null second),
    // so set the RHS directly from the previous window's terminal sensitivity.
    for( auto const& [row, st, nidx] : B.icRows ){
      auto it = sens_term.find( st );
      b( row ) = ( it != sens_term.end() && nidx < it->second.size() ) ? it->second[nidx] : 0.0;
    }
    // s = J^{-1} b via the stored LU (P*J = L*U): s = U \ ( L \ ( P*b ) )
    // rev50: representation-agnostic solve.  The SPARSE path can fail HERE rather than at
  // setup (spsolve throws on a singular factorisation), which is the tier-0 failure mode --
  // clean and attributable, unlike the dense approximate fallback it replaces.
  arma::vec s;
  if( !_sens_solve( B, b, s, /*transpose=*/false ) ) return false;
    std::vector<double> svec( nv, 0. );
    for( size_t i = 0; i < nv; ++i ) svec[i] = s( i );
    // Collapse the evolution domain to THIS window before the capture quadrature (and the terminal
    // profile below), exactly as the primal march does before _fill_captures.
    reset_evolution_window( B.t0, B.t1 - B.t0 );
    // Stage 1: accumulate this window's d(cap)/dp = INT_win s[R] and inject it into the cap input
    // slots so the contract's dF/dinp . dinp_dp term yields df/dcap . d(cap)/dp.  Captures feed
    // outputs only (dR/dinp[cap]=0), so this does not perturb the state sensitivity s.
    if( !_deferredValue.empty() ){
      _fill_capture_sens( svec.data(), (int)k, cap_window, dcap_acc.data() );
      for( auto const& C : _deferredValue ){
        if( !C.input.dag() ) continue;
        std::vector<double> const cv = get_input_values( C.input, dcap_acc.data() );
        set_input_values( C.input, cv, dinp_dp.data() );
      }
    }
    // contract dF/dp += dF/dinp . dinp_dp + dF/dvar . s  (point outputs captured in-window)
    for( auto const& f : _mFct ){
      bool const summed = ( f.row0 < nf && fctsum[f.row0] );
      // a field output over the evolution direction (2026-10-06, WORKPLAN 1.1): this window's rows, emitted in
      // window order at f.row0.., go to their global rows
      if( !summed && _march_field_over_evo( f ) ){
        std::vector<size_t> const rows = _march_field_rows( f, k );
        for( size_t i = 0; i < rows.size(); ++i ){
          if( rows[i] == npos ) continue;
          size_t const r = f.row0 + i, g = f.row0 + rows[i];
          if( r >= nf || g >= nf || r >= B.gfct.size() ) break;
          double contr = 0.;
          for( size_t idx = 0; idx < B.gfct[r].size(); ++idx ){ size_t c = B.cfct[r][idx]; double v = B.gfct[r][idx]; if( c < nv ) contr += v * s( c ); else if( c - nv < ni ) contr += v * dinp_dp[c-nv]; }
          dFdp[g] = contr;
        }
        continue;
      }
      bool const cap = summed || _march_point_capture( f, B.t0, B.t1 );
      for( size_t r = f.row0; r < f.row0 + f.nrow && r < nf && r < B.gfct.size(); ++r ){
        double contr = 0.;
        for( size_t idx = 0; idx < B.gfct[r].size(); ++idx ){ size_t c = B.cfct[r][idx]; double v = B.gfct[r][idx]; if( c < nv ) contr += v * s( c ); else if( c - nv < ni ) contr += v * dinp_dp[c-nv]; }
        if( summed )   dFdp[r] += contr;
        else if( cap ) dFdp[r]  = contr;
      }
    }
    // terminal(s_k) for the next block's transfer coupling (evolution window already collapsed above)
    auto sprof = terminal_profile( svec.data(), B.inp.empty() ? nullptr : B.inp.data(), nullptr );
    sens_term.clear();
    for( auto const& pr : _marchTransfers ){ auto it = sprof.find( pr.first ); sens_term[pr.first] = ( it != sprof.end() ) ? it->second : std::vector<double>{}; }
  }
  return true;
}

inline void
OCFESLV::terminal_profile_adjoint
( std::map<FFVar,std::vector<double>,lt_FFVar> const& lam, double* gvar )
const
{
  if( !gvar || !_evolution_dom_set || !_evolution_dom_var.dag() ) return;
  auto const evo_id = _evolution_dom_var.id();

  // Evolution-direction Lagrange weights at the UB (reference element; window-invariant).
  auto ievo = _mDom.find( _evolution_dom_var );
  if( ievo == _mDom.end() ) return;
  std::vector<double> w_evo;
  if( !ievo->second.w_eval_lagrange( w_evo, 1.0, 0 ) ) return;

  for( auto const& [state, lamvec] : lam ){
    if( lamvec.empty() ) continue;
    auto itv = _mVar.find( state );
    if( itv == _mVar.end() ) continue;
    auto const& var_dom = itv->second;                 // dep order (== eval_colloc dims order)

    // Per-domain node counts and first-innermost node strides (matches eval_colloc's block).
    size_t const D = var_dom.size();
    std::vector<size_t> nnode( D, 0 ), nStride( D, 0 );
    size_t evoDim = D;
    {
      size_t s = 1, d = 0;
      for( auto const& v : var_dom ){
        size_t const nn = _mDom.at( v ).n_node;
        nnode[d] = nn; nStride[d] = s; s *= nn;
        if( v.id() == evo_id ) evoDim = d;
        ++d;
      }
    }
    if( evoDim == D ) continue;                        // state not on the evolution domain

    // Spatial domains in var_dom order (skip evo): global-node count = n_elem*n_node each,
    // ordered first-spatial-domain innermost -- exactly terminal_profile()'s vals ordering.
    std::vector<FFVar> sVar; std::vector<size_t> sPos, sElem, sNode;
    {
      size_t d = 0;
      for( auto const& v : var_dom ){
        if( v.id() != evo_id ){ FFDom const& dd = _mDom.at( v ); sVar.push_back( v ); sPos.push_back( d ); sElem.push_back( dd.n_elem ); sNode.push_back( dd.n_node ); }
        ++d;
      }
    }
    size_t total = 1;
    for( size_t k = 0; k < sVar.size(); ++k ) total *= sElem[k] * sNode[k];
    if( lamvec.size() < total ) continue;

    for( size_t lin = 0; lin < total; ++lin ){
      std::map<FFVar,size_t,lt_FFVar> ndx_el;
      ndx_el[ _evolution_dom_var ] = 0;                // collapsed template: single evo element
      size_t node_base = 0, rem = lin;
      for( size_t k = 0; k < sVar.size(); ++k ){
        size_t const gtot = sElem[k] * sNode[k];
        size_t const g  = rem % gtot; rem /= gtot;
        size_t const ie = g / sNode[k];                // spatial element
        size_t const in = g % sNode[k];                // within-element node
        ndx_el[ sVar[k] ] = ie;
        node_base += in * nStride[ sPos[k] ];
      }
      size_t const inc = pos_state( state, ndx_el );
      double const lj  = lamvec[lin];
      size_t const nme = std::min( w_evo.size(), nnode[evoDim] );
      for( size_t m = 0; m < nme; ++m )
        gvar[ inc + node_base + m * nStride[evoDim] ] += w_evo[m] * lj;
    }
  }
}

inline bool
OCFESLV::_march_asens_apply
( std::vector<double> const& w, std::vector<double>& g )
{
  size_t const ni = n_colloc_inp(), nv = _nCollVar, nf = n_colloc_fct();
  g.assign( ni, 0. );
  // the capture -> window schedule, computed EXACTLY as _solve_marching / _march_fsens_apply do: the first window
  // whose upper bound reaches tau (UB-inclusive, LB-exclusive)
  std::vector<int> cap_window( _deferredValue.size(), -1 );
  if( !_deferredValue.empty() && _marchGrid.size() >= 2 ){
    double const span = std::fabs( _marchGrid.back() - _marchGrid.front() );
    double const eps  = 1e-9 * std::max( 1.0, span );
    for( size_t li = 0; li < _deferredValue.size(); ++li )
      if( !_deferredValue[li].accumulate )
        for( size_t kk = 0; kk + 1 < _marchGrid.size(); ++kk )
          if( _deferredValue[li].side == FFDom::PLUS
              ? ( _deferredValue[li].tau < _marchGrid[kk+1] - eps || kk+2 == _marchGrid.size() )   // tau^+: the window STARTING at tau (at tf: the last)
              : _deferredValue[li].tau <= _marchGrid[kk+1] + eps ){  // tau^-: the window ENDING at tau
            cap_window[li] = (int)kk; break; }
  }
  if( _fsBlocks.empty() || w.size() < nf ) return false;

  // Evolution-integral rows accumulate over all windows; else only the last window writes
  // (transpose of the forward's dFdp[r] = fctsum[r] ? dFdp[r]+contr : contr rule).
  std::vector<char> fctsum( nf, 0 );
  for( auto const& f : _mFct ){ bool const ei = _fct_is_evolution_integral( f.var ); for( size_t r = f.row0; r < f.row0 + f.nrow && r < nf; ++r ) fctsum[r] = ei ? 1 : 0; }

  size_t const N = _fsBlocks.size();
  std::map<FFVar,std::vector<double>,lt_FFVar> lam_term;   // terminal adjoint carried from window k+1

  // Stage 1b: adjoint of the capture quadrature.  lam_cap[R] = df/dcap^T w, read from the LAST window
  // (cap complete there).  ACCUM caps aggregate by source and scatter with jac*w_ref (quadrature)
  // every window.  LATCH caps must NOT aggregate by source (that loses tau); each is kept per-capture
  // and scattered with the Lagrange weights L_m(tau) at tau's window (the exact forward transpose).
  std::map<FFVar,std::vector<double>,lt_FFVar> lam_cap;                 // ACCUM only
  std::vector< std::pair<size_t,std::vector<double>> > lamLatch;        // (capture idx, adjoint on cap)
  if( !_deferredValue.empty() && N ){
    FSensBlock const& BL = _fsBlocks[N-1];
    std::vector<double> dummy( ni, 0. );
    for( size_t ci = 0; ci < _deferredValue.size(); ++ci ){
      t_DeferredValue const& C = _deferredValue[ci];
      if( !C.input.dag() || !C.source.dag() ) continue;
      auto itblk = _blkInp.find( C.input ); if( itblk == _blkInp.end() ) continue;
      size_t const cap_off = std::get<0>( itblk->second );
      size_t const cap_tot = get_input_values( C.input, dummy.data() ).size();
      if( !cap_tot ) continue;
      std::vector<double> lam( cap_tot, 0. );
      for( size_t r = 0; r < BL.gfct.size() && r < nf; ++r ){
        if( w[r] == 0. ) continue;
        for( size_t idx = 0; idx < BL.gfct[r].size(); ++idx ){
          size_t const c = BL.cfct[r][idx];
          if( c >= nv && (c-nv) >= cap_off && (c-nv) < cap_off + cap_tot )
            lam[ (c-nv) - cap_off ] += BL.gfct[r][idx] * w[r];
        }
      }
      if( C.accumulate ){
        auto& dst = lam_cap[ C.source ];
        if( dst.empty() ) dst = std::move( lam );
        else for( size_t i = 0; i < lam.size() && i < dst.size(); ++i ) dst[i] += lam[i];
      }
      else lamLatch.emplace_back( ci, std::move( lam ) );
    }
  }

  // _march_point_capture now owns each point output to a SINGLE window -- the one where its t is the UB
  // (LB-exclusive; the first window keeps its LB for a global-start point).  So an interior seam point t_m
  // is captured only in the window ending at t_m, matching the forward's in-window dFdp assignment, and
  // fwd == adjoint holds.  The capDone latch below is thus belt-and-braces: it guarantees a point row is
  // written in at most one window even if the predicate ever admitted two, and it lands a no-evolution-
  // point output (captured every window) on the highest-k == last window, exactly as the forward last-wins.
  // Summed (evolution-integral) rows still accumulate every window.  Inert for mono, integral-only, and
  // single-final-point models.
  std::vector<char> capDone( nf, 0 );
  for( size_t kk = 0; kk < N; ++kk ){
    size_t const k = N - 1 - kk;                          // reverse sweep
    bool const last = ( k + 1 == N );
    FSensBlock const& B = _fsBlocks[k];
    reset_evolution_window( B.t0, B.t1 - B.t0 );          // window domain for the capture scatter below

    // Per-row window gate: a point output is written only in the window whose evolution range
    // contains its t (so c(0) is picked up at window 0), evolution integrals every window.  This is
    // the transpose of the forward's in-window dFdp capture and of _solve_marching's value capture.
    std::vector<char> capW( nf, 0 );
    std::vector<double> weff( nf, 0. );    // this window's weight per WINDOW row (w[global row] for a captured row)
    for( auto const& f : _mFct ){
      bool const summed = ( f.row0 < nf && fctsum[f.row0] );
      // a field output over the evolution direction (2026-10-06, WORKPLAN 1.1): window row f.row0+i is global row
      // f.row0+rows[i] -- it takes that row's weight; the window's other rows of f carry none
      if( !summed && _march_field_over_evo( f ) ){
        std::vector<size_t> const rows = _march_field_rows( f, k );
        for( size_t i = 0; i < rows.size(); ++i ){
          if( rows[i] == npos ) continue;
          size_t const r = f.row0 + i, g = f.row0 + rows[i];
          if( r < nf && g < nf ){ capW[r] = 1; weff[r] = w[g]; }
        }
        continue;
      }
      bool const pt     = _march_point_capture( f, B.t0, B.t1 );
      for( size_t r = f.row0; r < f.row0 + f.nrow && r < nf; ++r ){
        if( summed )                  capW[r] = 1;                       // evolution integral: every window
        else if( pt && !capDone[r] ){ capW[r] = 1; capDone[r] = 1; }     // point: highest-k window only
        else                          capW[r] = 0;
        weff[r] = capW[r]? w[r]: 0.;
      }
    }

    // state adjoint  ls = dF/dvar^T wk
    arma::vec ls( nv, arma::fill::zeros );
    for( size_t r = 0; r < B.gfct.size() && r < nf; ++r ){
      double const wk = weff[r];
      if( wk == 0. ) continue;
      for( size_t idx = 0; idx < B.gfct[r].size(); ++idx ){ size_t const c = B.cfct[r][idx]; if( c < nv ) ls( c ) += B.gfct[r][idx] * wk; }
    }
    // Stage 1b: add the capture adjoint -- lam_cap[R] scattered onto ls[R] with jac*w_ref over this
    // window (transpose of d(cap)/dp = INT jac*w_ref*s[R]).
    if( !lam_cap.empty() ){
      std::vector<double> lsc( nv, 0. );
      _fill_capture_sens_adjoint( lam_cap, lsc.data() );
      for( size_t i = 0; i < nv; ++i ) ls( i ) += lsc[i];
    }
    // LATCH captures: exact transpose of the forward eval_colloc(src,tau) -- scatter each latch
    // cap with the Lagrange weights L_m(tau) at tau's element.  The window domain is set above by
    // reset_evolution_window, so the method no-ops unless tau lies in THIS window (transpose of the
    // forward cap_window==k gate).  Per-capture, so tau is not lost to source-aggregation.
    for( auto const& [ci, lamC] : lamLatch ){
      // ONE window per capture, the forward's: _fill_capture_sens_adjoint_latch alone accepts tau in any element of
      // the window with BOTH ends inclusive, so a tau on a window seam (the UB of window k, the LB of k+1) was
      // scattered in both and its adjoint DOUBLED (OCFE_EVALSUM, x0^2 at T1, 2026-09-28).
      if( ci < cap_window.size() && cap_window[ci] != (int)k ) continue;
      _fill_capture_sens_adjoint_latch( _deferredValue[ci].source, _deferredValue[ci].tau, lamC, ls.memptr(), _deferredValue[ci].side );
    }

    // coupling from window k+1: += M^T ( gathered ic adjoint )  into the state adjoint
    if( !last && !lam_term.empty() ){
      std::vector<double> lsv( nv, 0. );
      terminal_profile_adjoint( lam_term, lsv.data() );
      for( size_t i = 0; i < nv; ++i ) ls( i ) += lsv[i];
    }

    // lb = J_k^{-T} ls = P^T ( L^{-T} ( U^{-T} ls ) ) via the stored LU (P J = L U)
    arma::vec lb;
    if( !_sens_solve( B, ls, lb, /*transpose=*/true ) ) return false;

    // input-space adjoint  ldinp = dF/dinp^T wk - dR/dinp^T lb
    std::vector<double> ldinp( ni, 0. );
    for( size_t r = 0; r < B.gfct.size() && r < nf; ++r ){
      double const wk = weff[r];
      if( wk == 0. ) continue;
      for( size_t idx = 0; idx < B.gfct[r].size(); ++idx ){ size_t const c = B.cfct[r][idx]; if( c >= nv && c - nv < ni ) ldinp[c-nv] += B.gfct[r][idx] * wk; }
    }
    for( size_t r = 0; r < B.geqn.size(); ++r )
      for( size_t idx = 0; idx < B.geqn[r].size(); ++idx ){ size_t const c = B.ceqn[r][idx]; if( c >= nv && c - nv < ni ) ldinp[c-nv] -= B.geqn[r][idx] * lb( r ); }

    // IC split: for k>0 the ic-input entries came from s_{k-1}'s terminal, not the control
    // seed -- gather them into lam_term (for window k-1) and zero them out of ldinp so they
    // do not leak into the control-seed gradient g.  For k==0 the whole seed is the control.
    if( k > 0 ){
      lam_term.clear();
      for( auto const& pr : _marchTransfers ){
        if( !pr.second.dag() ) continue;               // general-IC override handled below
        std::vector<double> blk = get_input_values( pr.second, ldinp.data() );
        lam_term[ pr.first ] = blk;
        set_input_values( pr.second, std::vector<double>( blk.size(), 0. ), ldinp.data() );
      }
      // General-IC override rows: the continuity row's adjoint lb(row) IS the terminal adjoint for
      // window k-1 (transpose of the forward b(row)=d(terminal)/dp).  There is no ic_input to carry
      // it, so gather directly; nothing leaks into ldinp since the continuity row has no input column.
      for( auto const& [row, st, nidx] : B.icRows ){
        auto& lt = lam_term[st];
        if( lt.size() <= nidx ) lt.resize( nidx + 1, 0. );
        lt[nidx] = lb( row );
      }
    }

    if( k > 0 && !_march_jump_asens( k, lam_term, ldinp ) ) return false;   // through a transition at this seam
    for( size_t i = 0; i < ni; ++i ) g[i] += ldinp[i];
  }
  return true;
}


inline bool
OCFESLV::_mono_fsens_setup
( double* var, double const* inp, double const* cst )
{
  _monoValid = false;
  _monoFct.clear();
  size_t const ni  = n_colloc_inp();
  size_t const nv  = _nCollVar;
  size_t const neq = n_colloc_eqn();
  size_t const nf  = n_colloc_fct();

  std::vector<double> inpwork( ni, 0. );
  if( inp ) std::copy( inp, inp + ni, inpwork.begin() );
  double* pinp = ni ? inpwork.data() : nullptr;

  // Monolithic primal solve (marchGrid empty / SOLVE.MARCHING off => solve() is monolithic).
  SolveReport const r = solve( var, pinp, cst );
  if( !r.converged ) return false;

  // Fill captured-reduction inputs over the full horizon so the block linearization differentiates
  // f(cap) at the correct cap value (single overwriting fill, k=-1).
  _fill_captures( var, pinp, cst, -1, std::vector<int>() );

  // Single-block linearization at the converged solution.
  FSensBlock B;
  B.var.assign( var, var + nv );
  if( pinp ) B.inp.assign( pinp, pinp + ni );
  std::vector<size_t> nnzeqn( neq, 0 ), nnzfct( nf, 0 );
  if( !deriv( nullptr, nullptr, nnzeqn.data(), nullptr, nullptr, nullptr, nnzfct.data(), nullptr, var, pinp, cst ) ) return false;
  B.geqn.resize( neq ); B.ceqn.resize( neq ); B.gfct.resize( nf ); B.cfct.resize( nf );
  std::vector<double*> geqnp( neq, nullptr ), gfctp( nf, nullptr );
  std::vector<size_t*> ceqnp( neq, nullptr ), cfctp( nf, nullptr );
  for( size_t r = 0; r < neq; ++r ){ B.geqn[r].resize( nnzeqn[r] ); B.ceqn[r].resize( nnzeqn[r] ); geqnp[r] = B.geqn[r].data(); ceqnp[r] = B.ceqn[r].data(); }
  for( size_t r = 0; r < nf;  ++r ){ B.gfct[r].resize( nnzfct[r] ); B.cfct[r].resize( nnzfct[r] ); gfctp[r] = B.gfct[r].data(); cfctp[r] = B.cfct[r].data(); }
  std::vector<double> eqnbuf( neq, 0. ), fctbuf( nf, 0. );
  if( !deriv( eqnbuf.data(), geqnp.data(), nnzeqn.data(), ceqnp.data(), fctbuf.data(), gfctp.data(), nnzfct.data(), cfctp.data(), var, pinp, cst ) ) return false;

  {
    auto const _t0 = std::chrono::steady_clock::now();
    bool const ok = _sens_factorise( B, neq, nv, nnzeqn, "mono" );
    if( _knob_AUDIT_SENS_NULL() )
      _disp( 3 ) << "  [sens-time] mono factorise/store "
                << std::chrono::duration<double,std::milli>(
                     std::chrono::steady_clock::now() - _t0 ).count() << " ms"
                << ( options.SENSITIVITY.SPARSE && !_knob_SENS_MINNORM_PINV() ? " (sparse)" : " (dense)" )
                << std::endl;
    if( !ok ) return false;
  }

  _monoBlock = std::move( B );
  _monoFct.assign( fctbuf.begin(), fctbuf.end() );
  _sensFct = _monoFct;                       // unified sens function-value cache
  _monoValid = true;
  return true;
}

inline bool
OCFESLV::_mono_fsens_apply
( std::vector<double> const& dinp, std::vector<double>& dFdp )
{
  size_t const ni = n_colloc_inp(), nv = _nCollVar, nf = n_colloc_fct();
  dFdp.assign( nf, 0. );
  if( !_monoValid ) return false;
  FSensBlock const& B = _monoBlock;

  // RHS b = -(dR/dinp . dinp)
  arma::vec b( B.geqn.size(), arma::fill::zeros );
  for( size_t r = 0; r < B.geqn.size(); ++r )
    for( size_t idx = 0; idx < B.geqn[r].size(); ++idx ){ size_t c = B.ceqn[r][idx]; if( c >= nv && c - nv < ni ) b( r ) -= B.geqn[r][idx] * dinp[c-nv]; }
  // s = J^{-1} b via stored LU (P J = L U): s = U \ ( L \ ( P b ) )
  // rev50: representation-agnostic solve.  The SPARSE path can fail HERE rather than at
  // setup (spsolve throws on a singular factorisation), which is the tier-0 failure mode --
  // clean and attributable, unlike the dense approximate fallback it replaces.
  arma::vec s;
  if( !_sens_solve( B, b, s, /*transpose=*/false ) ) return false;
  // Stage 1: inject d(cap)/dp = INT s[R] (single full-horizon fill) into the cap input slots on a
  // mutable copy, so the contract's dF/dinp term yields df/dcap . d(cap)/dp.
  std::vector<double> dinp_w = dinp;
  if( !_deferredValue.empty() ){
    std::vector<double> svec( nv, 0. );
    for( size_t i = 0; i < nv; ++i ) svec[i] = s( i );
    _fill_capture_sens( svec.data(), -1, std::vector<int>(), dinp_w.data() );
  }
  // dFdp = dF/dinp . dinp + dF/dvar . s   (single block: direct, no fctsum accumulation)
  for( size_t r = 0; r < B.gfct.size() && r < nf; ++r ){
    double contr = 0.;
    for( size_t idx = 0; idx < B.gfct[r].size(); ++idx ){ size_t c = B.cfct[r][idx]; double v = B.gfct[r][idx]; if( c < nv ) contr += v * s( c ); else if( c - nv < ni ) contr += v * dinp_w[c-nv]; }
    dFdp[r] = contr;
  }
  return true;
}

inline bool
OCFESLV::_mono_asens_apply
( std::vector<double> const& w, std::vector<double>& g )
{
  size_t const ni = n_colloc_inp(), nv = _nCollVar, nf = n_colloc_fct();
  g.assign( ni, 0. );
  if( !_monoValid || w.size() < nf ) return false;
  FSensBlock const& B = _monoBlock;

  // state adjoint ls = dF/dvar^T w
  arma::vec ls( nv, arma::fill::zeros );
  for( size_t r = 0; r < B.gfct.size() && r < nf; ++r ){
    if( w[r] == 0. ) continue;
    for( size_t idx = 0; idx < B.gfct[r].size(); ++idx ){ size_t c = B.cfct[r][idx]; if( c < nv ) ls( c ) += B.gfct[r][idx] * w[r]; }
  }
  // Stage 1b: capture adjoint -- lam_cap[R] = df/dcap^T w scattered onto ls[R] with jac*w_ref over the
  // full horizon (transpose of the monolithic d(cap)/dp fill).
  if( !_deferredValue.empty() ){
    std::map<FFVar,std::vector<double>,lt_FFVar> lam_cap;   // ACCUM only
    std::vector<double> dummy( ni, 0. );
    std::vector<double> ls_latch( nv, 0. );                // LATCH contribution (L_m(tau) transpose)
    for( auto const& C : _deferredValue ){
      if( !C.input.dag() || !C.source.dag() ) continue;
      auto itblk = _blkInp.find( C.input ); if( itblk == _blkInp.end() ) continue;
      size_t const cap_off = std::get<0>( itblk->second );
      size_t const cap_tot = get_input_values( C.input, dummy.data() ).size();
      if( !cap_tot ) continue;
      std::vector<double> lam( cap_tot, 0. );
      for( size_t r = 0; r < B.gfct.size() && r < nf; ++r ){
        if( w[r] == 0. ) continue;
        for( size_t idx = 0; idx < B.gfct[r].size(); ++idx ){
          size_t const c = B.cfct[r][idx];
          if( c >= nv && (c-nv) >= cap_off && (c-nv) < cap_off + cap_tot )
            lam[ (c-nv) - cap_off ] += B.gfct[r][idx] * w[r];
        }
      }
      if( C.accumulate ){
        auto& dst = lam_cap[ C.source ];
        if( dst.empty() ) dst = std::move( lam );
        else for( size_t i = 0; i < lam.size() && i < dst.size(); ++i ) dst[i] += lam[i];
      }
      else                                                  // LATCH: interp transpose at tau (full horizon)
        _fill_capture_sens_adjoint_latch( C.source, C.tau, lam, ls_latch.data(), C.side );
    }
    if( !lam_cap.empty() ){
      std::vector<double> lsc( nv, 0. );
      _fill_capture_sens_adjoint( lam_cap, lsc.data() );
      for( size_t i = 0; i < nv; ++i ) ls( i ) += lsc[i];
    }
    for( size_t i = 0; i < nv; ++i ) ls( i ) += ls_latch[i];   // add LATCH interp-transpose contribution
  }
  // lb = J^{-T} ls = P^T ( L^{-T} ( U^{-T} ls ) )
  // rev50: adjoint solves J^T lambda = ls.  NOTE the transpose is free for LU (same
  // triangular pair), free for pinv, and cheap-but-not-free for sparse (an explicit .t()).
  // It is NOT free for QR -- tier 1 will need a SECOND factorisation of J^T.
  arma::vec lb;
  if( !_sens_solve( B, ls, lb, /*transpose=*/true ) ) return false;
  // g = dF/dinp^T w - dR/dinp^T lb
  for( size_t r = 0; r < B.gfct.size() && r < nf; ++r ){
    if( w[r] == 0. ) continue;
    for( size_t idx = 0; idx < B.gfct[r].size(); ++idx ){ size_t c = B.cfct[r][idx]; if( c >= nv && c - nv < ni ) g[c-nv] += B.gfct[r][idx] * w[r]; }
  }
  for( size_t r = 0; r < B.geqn.size(); ++r )
    for( size_t idx = 0; idx < B.geqn[r].size(); ++idx ){ size_t c = B.ceqn[r][idx]; if( c >= nv && c - nv < ni ) g[c-nv] -= B.geqn[r][idx] * lb( r ); }
  return true;
}

inline bool

//==========================================================================================================
// SOLVERS -- Newton / Levenberg-Marquardt, the marching loop, continuation and sensitivities.
// (banner placed from the definitions themselves: this region is lines 18301-20162 of the pre-rev268 file)
//==========================================================================================================

OCFESLV::solve_fsens
( double* var, double const* inp, double const* cst )
{
  FixedInputScope fixed_scope_( _fixedInputDepth );  std::vector<double> fixed_buf_;  inp = _with_fixed_inputs( inp, fixed_buf_ );   // fix_input: model-supplied values, outermost call only
  t_ThreadCap const _tcap( options.MAXTHREAD, options.MAXTHREAD != 1 || options.DISPLAY_LEVEL >= 1 );   // rev142b: covers the sensitivity
                                                  // factorizations, not just the primal
  // Primal solve + factorization (auto march/mono); leaves sens_functions() populated.
  _sensJac.clear();
  if( !_sens_setup( var, inp, cst ) ) return false;

  size_t const ncd = n_control_dof(), nf = n_colloc_fct(), nci = n_colloc_inp();
  _sensJac.assign( nf * ncd, 0. );
  if( !ncd ) return true;                       // functions still available; no controls registered

  // One forward directional apply per control DOF -> a COLUMN of dF/dp (control space).
  std::vector<double> dp( ncd, 0. ), dinp( nci ), g;
  for( size_t j = 0; j < ncd; ++j ){
    std::fill( dp.begin(), dp.end(), 0. ); dp[j] = 1.0;
    std::fill( dinp.begin(), dinp.end(), 0. );
    decode_controls( dp, dinp.data() );
    if( !_solve_fsens_apply( dinp, g ) ){ _sensJac.clear(); return false; }
    for( size_t i = 0; i < nf && i < g.size(); ++i ) _sensJac[ i*ncd + j ] = g[i];
  }
  return true;
}

inline bool
OCFESLV::solve_asens
( double* var, double const* inp, double const* cst )
{
  FixedInputScope fixed_scope_( _fixedInputDepth );  std::vector<double> fixed_buf_;  inp = _with_fixed_inputs( inp, fixed_buf_ );   // fix_input: model-supplied values, outermost call only
  t_ThreadCap const _tcap( options.MAXTHREAD, options.MAXTHREAD != 1 || options.DISPLAY_LEVEL >= 1 );   // rev142b: covers the sensitivity
                                                  // factorizations, not just the primal
  // Primal solve + factorization (auto march/mono); leaves sens_functions() populated.
  _sensJac.clear();
  if( !_sens_setup( var, inp, cst ) ) return false;

  size_t const ncd = n_control_dof(), nf = n_colloc_fct();
  _sensJac.assign( nf * ncd, 0. );
  if( !ncd ) return true;

  // One reverse sweep per output function -> a ROW of dF/dp (control space, via encode_controls,
  // the transpose of decode_controls).  Cheaper than forward when nf < ncd.
  std::vector<double> w( nf, 0. ), g, gp;
  for( size_t i = 0; i < nf; ++i ){
    std::fill( w.begin(), w.end(), 0. ); w[i] = 1.0;
    if( !_solve_asens_apply( w, g ) ){ _sensJac.clear(); return false; }
    encode_controls( g.data(), gp );
    for( size_t j = 0; j < ncd && j < gp.size(); ++j ) _sensJac[ i*ncd + j ] = gp[j];
  }
  return true;
}


////////////////////////////////////////////////////////////////////////////////////////////////
//                          HOMOTOPY / CONTINUATION IMPLEMENTATION
////////////////////////////////////////////////////////////////////////////////////////////////

inline bool
OCFESLV::add_homotopy( FFVar const& param, double from, double to, int stage )
{
  if( _mInp.find( param ) == _mInp.end() ) return false;   // not a registered input
  t_Homotopy h; h.param = param; h.from = from; h.to = to; h.stage = stage;
  _vHom.push_back( h );
  return true;
}

inline bool
OCFESLV::add_homotopy( FFVar const& param, std::function<double(double)> map, int stage )
{
  if( _mInp.find( param ) == _mInp.end() || !map ) return false;
  t_Homotopy h; h.param = param; h.map = std::move( map ); h.stage = stage;
  _vHom.push_back( h );
  return true;
}

inline void
OCFESLV::reset_homotopy()
{
  _vHom.clear();
  _mHomCap.clear();
}

inline double
OCFESLV::_homotopy_resid( double const* var, double const* inp, double const* cst )
const
{
  size_t const nr = n_colloc_eqn();
  if( !nr ) return 0.;
  std::vector<double> r( nr, 0. );
  if( !const_cast<OCFESLV*>(this)->eval( r.data(), nullptr, var, inp, cst ) ) return 1e300;
  double m = 0.;
  for( double v : r ) m = std::max( m, std::fabs( v ) );
  return m;
}

inline void
OCFESLV::set_homotopy_cap( int stage, double cap )
{
  if( cap > 0. ) _mHomCap[ stage ] = cap;
  else           _mHomCap.erase( stage );
}

inline bool
OCFESLV::_homotopy_ramp( std::vector<t_Homotopy const*> const& act, double s0, double s1, double cap,
                       double* var, double* pinp, double const* cst,
                       std::vector<double>& saved, StageReport& st, double& s_out )
{
  // 2026-10-10: the FULL solve vector -- the state coefficients AND the trace/tau unknowns appended to var[] (the
  // claims kept explicit under IC_STRONG, every tau under IC_TRACE).  With n_colloc_sta() alone the traces were
  // neither saved nor restored: after a failed trial they kept the diverged trial's values, and every later trial
  // started from them -- the ramp never recovered from its first failure (OCFE_MMPDE2: 3840 states + 30 traces).
  size_t const nv = n_colloc_sta() + n_colloc_trace();
  double cur = s0, d = cap;
  s_out = cur;

  // History of accepted points, most recent LAST, for the guarded extrapolation: two points give the
  // secant, three a quadratic.  Local to this leg -- the path bends where the active set changes.
  std::vector<double> sHist{ cur };
  std::vector< std::vector<double> > xHist{ saved };
  std::vector<double> cand( nv, 0. );

  while( std::abs( cur - s1 ) > 1e-9 ){
    double const trial = cur + std::min( d, s1 - cur );
    for( auto const* h : act )
      set_input_values( h->param, { _homotopy_value( *h, trial ) }, pinp );

    // DO-NO-HARM PREDICTOR: score plain reuse and every available extrapolation order at `trial`,
    // keep the best.  A worse candidate costs one residual evaluation and is discarded.
    std::copy( saved.begin(), saved.end(), var );
    if( options.HOMOTOPY.PREDICT && options.HOMOTOPY.PREDICT_ORDER > 0 && xHist.size() >= 2 ){
      double best = _homotopy_resid( saved.data(), pinp, cst );
      bool   used = false;
      int    bestord = 0;
      size_t const H = xHist.size();
      int const ordmax = std::min( options.HOMOTOPY.PREDICT_ORDER, (int)H - 1 );
      for( int ord = 1; ord <= ordmax && ord <= 2; ++ord ){
        size_t const n = (size_t)ord + 1;
        size_t const b = H - n;
        bool ok = true;
        for( size_t i = 0; i < nv; ++i ) cand[i] = 0.;
        for( size_t k = 0; k < n && ok; ++k ){
          double L = 1.;
          for( size_t j = 0; j < n; ++j ){
            if( j == k ) continue;
            double const den = sHist[b+k] - sHist[b+j];
            if( std::abs( den ) < 1e-12 ){ ok = false; break; }
            L *= ( trial - sHist[b+j] ) / den;
          }
          if( !ok ) break;
          for( size_t i = 0; i < nv; ++i ) cand[i] += L * xHist[b+k][i];
        }
        if( !ok ) continue;
        double const rC = _homotopy_resid( cand.data(), pinp, cst );
        if( rC < best ){ best = rC; std::copy( cand.begin(), cand.end(), var ); used = true; bestord = ord; }
      }
      if( used ){
        ++st.pred_used;
        if     ( bestord == 1 ) ++st.pred_ord1;
        else if( bestord == 2 ) ++st.pred_ord2;
      }
      else ++st.pred_rejected;
    }

    auto tic = std::chrono::steady_clock::now();
    SolveReport sr = solve( var, pinp, cst );
    st.seconds += std::chrono::duration<double>( std::chrono::steady_clock::now() - tic ).count();
    ++st.solves; st.iterations += sr.iterations; st.last_residual = sr.final_residual;

    if( sr.converged ){
      cur = trial; s_out = cur;
      saved.assign( var, var + nv );
      sHist.push_back( cur ); xHist.push_back( saved );
      if( sHist.size() > 3 ){ sHist.erase( sHist.begin() ); xHist.erase( xHist.begin() ); }
      ++st.accepts;
      _disp( 2 ) << "OCFESLV::solve_homotopy s=" << cur << " accepted (" << sr.iterations << " iterations)\n";
      d = std::min( cap, d * options.HOMOTOPY.STEP_GROW );
    }
    else{
      std::copy( saved.begin(), saved.end(), var );      // restore the last converged state
      for( auto const* h : act )
        set_input_values( h->param, { _homotopy_value( *h, cur ) }, pinp );
      ++st.backtracks;
      d *= options.HOMOTOPY.STEP_SHRINK;
      if( d < options.HOMOTOPY.STEP_MIN ){ s_out = cur; return false; }
    }
  }
  return true;
}

inline OCFESLV::SolveReport
OCFESLV::_solve_homotopy( double* var, double const* inp, double const* cst )
{
  ContinuationReport rep;
  SolveReport agg;

  // Drive the parameters on a PRIVATE copy of the caller's input array, so a continuation solve
  // leaves the caller's inp untouched.  The parameters finish at their registered `to` values,
  // which the caller already knows, so nothing is lost by not writing back.
  std::vector<double> inpwrk;
  size_t const ni = n_colloc_inp();
  if( ni ){
    inpwrk.assign( ni, 0. );
    if( inp ) std::copy( inp, inp + ni, inpwrk.begin() );
  }
  double* pinp = ni ? inpwrk.data() : nullptr;

  _inHomotopy = true;                                   // inner solve() must not re-enter here
  struct Guard { bool& f; ~Guard(){ f = false; } } guard{ _inHomotopy };

  // Stage order: ascending, preserving declaration order within a stage.
  std::vector<int> stages;
  for( auto const& h : _vHom )
    if( std::find( stages.begin(), stages.end(), h.stage ) == stages.end() ) stages.push_back( h.stage );
  std::sort( stages.begin(), stages.end() );

  // All parameters start at their s = 0 value, so the model is consistent before stage 0 runs.
  for( auto const& h : _vHom )
    set_input_values( h.param, { _homotopy_value( h, 0. ) }, pinp );

  // 2026-10-10 (4R.18): SOLVE the s = 0 problem before the ramp.  The ramp takes `saved` as the solution at its
  // start and restarts every trial from it; the caller's vector is only an initial guess, so without this every
  // trial was a cold start, and a failed one shrank the step from a point that was never a solution (OCFE_MMPDE2).
  // A failure here is said as such: the continuation has no starting point, which is not a failure to move.
  if( options.HOMOTOPY.SOLVE_S0 ){
    SolveReport const sr0 = solve( var, pinp, cst );
    rep.solves += 1; rep.iterations += sr0.iterations; rep.final_residual = sr0.final_residual;
    if( !sr0.converged ){
      _disp( 1 ) << "OCFESLV::solve_homotopy ** the s = 0 problem did not converge (final residual "
                 << sr0.final_residual << "): the continuation has no starting point\n";
      rep.converged = false; rep.failed_at_s = 0.;
      _homReport = rep;
      agg.converged = false; agg.iterations = rep.iterations; agg.final_residual = rep.final_residual;
      return agg;
    }
    _disp( 2 ) << "OCFESLV::solve_homotopy s=0 solved (" << sr0.iterations << " iterations)\n";
  }
  std::vector<double> saved( var, var + n_colloc_sta() + n_colloc_trace() );

  auto accumulate = [&]( StageReport const& st ){
    rep.stage.push_back( st );
    rep.solves += st.solves; rep.accepts += st.accepts; rep.backtracks += st.backtracks;
    rep.iterations += st.iterations; rep.pred_used += st.pred_used;
    rep.pred_rejected += st.pred_rejected; rep.seconds += st.seconds;
    rep.pred_ord1 += st.pred_ord1; rep.pred_ord2 += st.pred_ord2;
    rep.final_residual = st.last_residual;
  };

  for( int sg : stages ){
    StageReport st; st.stage = sg;
    std::vector<t_Homotopy const*> act;
    for( auto const& h : _vHom ) if( h.stage == sg ) act.push_back( &h );
    st.nparam = (int)act.size();
    double const cap = _homotopy_cap( sg );

    double s_out = 0.;
    bool ok = _homotopy_ramp( act, 0., 1., cap, var, pinp, cst, saved, st, s_out );

    // AUTO-DESTAGE: a multi-parameter stage that cannot move its parameters together is retried with
    // them moved one at a time, from wherever the joint ramp got to.  Announced loudly, because the
    // usual cause is the MODEL -- those parameters interact too strongly to share a stage.
    if( !ok && options.HOMOTOPY.AUTO_DESTAGE && act.size() > 1 ){
      _disp( 2 ) << "OCFESLV::solve_homotopy ** stage " << sg << " FAILED to move its "
                << act.size() << " parameters together (stalled at s=" << s_out
                << ") -- AUTO-DESTAGING into " << act.size() << " sequential sub-stage(s).\n"
                << "   This is a RECOVERY, not a fix: those parameters interact too strongly to share\n"
                << "   a stage.  Consider declaring them in separate stages, or lowering the step cap\n"
                << "   for this stage with set_homotopy_cap(" << sg << ", ...).\n";
      double const s_join = s_out;
      ok = true;
      for( size_t k = 0; k < act.size() && ok; ++k ){
        std::vector<t_Homotopy const*> one{ act[k] };
        double sub = s_join;
        ok = _homotopy_ramp( one, s_join, 1., cap, var, pinp, cst, saved, st, sub );
        if( !ok ) s_out = sub;
      }
      if( ok ) st.destaged = (int)act.size();
    }

    if( !ok ){
      st.s_reached = s_out;
      accumulate( st );
      rep.failed_stage = sg; rep.failed_at_s = s_out; rep.converged = false;
      _homReport = rep;
      agg.converged = false; agg.iterations = rep.iterations;
      agg.final_residual = rep.final_residual;
      if( _solve_display_level() >= 1 )
        _disp( 2 ) << "OCFESLV::solve_homotopy ** stage " << sg << " FAILED at s=" << s_out
                  << " after " << st.solves << " solves (" << st.iterations << " its)\n";
      return agg;
    }

    st.converged = true; st.s_reached = 1.;
    accumulate( st );

    if( _solve_display_level() >= 1 )
      _disp( 3 ) << "OCFESLV::solve_homotopy [stage " << sg << "] " << st.nparam << " param(s)"
                << ( st.destaged ? " DESTAGED" : "" ) << "  cap=" << cap << "  "
                << st.solves << " solves  " << st.accepts << " acc  " << st.backtracks << " bt  "
                << st.pred_used << " pred+ (" << st.pred_ord1 << " lin/" << st.pred_ord2 << " quad)  "
                << st.pred_rejected << " pred-  "
                << st.iterations << " its  " << st.seconds << " s\n";
  }

  rep.converged = true;
  _homReport = rep;
  agg.converged = true; agg.iterations = rep.iterations;
  agg.final_residual = rep.final_residual;
  if( _solve_display_level() >= 1 )
    _disp( 2 ) << "OCFESLV::solve_homotopy ** " << rep.stage.size() << " stage(s), "
              << rep.solves << " solves, " << rep.iterations << " its, "
              << rep.seconds << " s\n";
  return agg;
}

inline OCFESLV::SolveReport
OCFESLV::solve( double* var, double const* inp, double const* cst )
{
  FixedInputScope fixed_scope_( _fixedInputDepth );  std::vector<double> fixed_buf_;  inp = _with_fixed_inputs( inp, fixed_buf_ );   // fix_input: model-supplied values, outermost call only
  t_ThreadCap const _tcap( options.MAXTHREAD, options.MAXTHREAD != 1 || options.DISPLAY_LEVEL >= 1 );   // rev142b: see setup()
  // rev115: force solver verbosity from the environment HERE, at the entry to solve(),
  // because this is the last point before the flag is read.
  //
  // TWO EARLIER SITES FAILED, both silently:
  //   rev114, Options constructor  -- OCFE_MBC5.cpp:204 does `oc.options.SOLVE.VERBOSE
  //                                   = false;` afterwards and overwrote it.
  //   rev115a, Options adopt op    -- the driver writes the FIELD directly, so no copy
  //                                   operator runs and the hook never executed.
  // Both produced a run that looked exactly like the run without the variable, which is
  // the failure mode a diagnostic switch must not have.  Force-ON only: a driver that
  // quietens itself keeps that unless a diagnostic run asks otherwise.
  Options::_force_verbose_from_env( options.SOLVE.VERBOSE );

  // HOMOTOPY: when continuation parameters were registered with add_homotopy() and we are not
  // already inside the ramp, hand off to the continuation driver.  This sits OUTSIDE the marching
  // dispatch on purpose -- the ramp is the outer loop, and each of its steps re-enters solve()
  // here, where a marching model then takes the marching path as usual.  Sensitivity solves that
  // re-enter solve() (e.g. via _mono_fsens_setup) therefore inherit the continuation for free.
  if( options.HOMOTOPY.ENABLE && !_vHom.empty() && !_inHomotopy )
    return _solve_homotopy( var, inp, cst );

  // SOLVE.MARCHING: when a march grid was recorded at setup() and we are not already
  // inside the march driver, hand off to the element-by-element march.  Otherwise (and
  // for each per-element block solve) run the monolithic solver below unchanged.
  if( options.SOLVE.MARCHING && !_marchGrid.empty() && !_inMarch )
    return _solve_marching( var, inp, cst );
  _lastSolveMarched = false;   // reached only for a monolithic (or per-window) solve

  // Captured reductions (and any other reduce_order-generated inputs) make n_colloc_inp()>0 without
  // the driver knowing.  A monolithic solve needs a non-null input array (deriv/eval require one);
  // marching threads its own mutable pinp, so this only bites the direct path.  Self-allocate a zero
  // buffer when the caller passed none -- the captures are filled post-solve (_eval_functions_captured)
  // and the equations do not reference them, so zeros are correct during the Newton.
  std::vector<double> inpauto;
  if( !inp && n_colloc_inp() ){ inpauto.assign( n_colloc_inp(), 0. ); inp = inpauto.data(); }

  // Primal-solve reuse: a genuine (non-marching, non-nested) monolithic solve at the SAME input
  // as the last converged one returns the stored solution directly, skipping the Newton -- so the
  // gradient path (which re-enters solve() through _mono_fsens_setup) reuses the value path's
  // converged point and only pays the factorization.  Functions are still (cheaply) refreshed.
  if( options.SOLVE.REUSE && !_inMarch && _primalValid
   && _primalVar.size() == n_colloc_sta() && _primalInp.size() == n_colloc_inp()
   && ( n_colloc_inp() == 0 || ( inp && std::equal( _primalInp.begin(), _primalInp.end(), inp ) ) ) ){
    std::copy( _primalVar.begin(), _primalVar.end(), var );
    SolveReport rep; rep.converged = true; rep.reused = true;
    _eval_functions_captured( var, inp, cst );
    return rep;
  }

  SolveReport rep;
  size_t const n = n_colloc_eqn();
  size_t const m = n_colloc_sta();

  // rev137 (rev B): the assembly lives in _assemble_core(); solve() owns one
  // t_JacAssembly for the whole Newton -- pattern queried once, buffers reused, the
  // marching override and the equilibration in exactly one place (the code the
  // setup-time determinacy audit assembles with; cross-checked IDENTICAL against the
  // former inline block on 5337 solves before the switch).
  t_JacAssembly A;
  if( !_assemble_pattern( A ) ) return rep;
  std::vector<double> r( n, 0.0 ), rtrial( n, 0.0 ), xt( m, 0.0 );
  bool const ic_override = _inMarch && _marchICActive && !_marchICRows.empty();
  _ensure_ic_sat( A, ic_override );

  // ---- Marching override diagnostic (DISPLAY_LEVEL>=1): localise the IC_WEAK interface failure ----
  if( ic_override && _solve_display_level() >= 1 ){
    size_t n_terms = _plan.weak_sat_terms().size(), n_valid = 0;
    for( auto const& term : _plan.weak_sat_terms() ){
      if( term.sat_coeffs.empty() ) continue;
      if( term.receiver_emit_flat < n ) ++n_valid;
    }
    size_t ovr_n = 0, multi = 0, ambig = 0; double tmax = 0.0;
    for( auto const& [st, rows] : _marchICRows ){
      auto itt = _marchICTerminal.find( st );
      if( itt != _marchICTerminal.end() ) for( double v : itt->second ) tmax = std::max( tmax, std::fabs(v) );
      auto itc = _marchICColToNode.find( st );
      for( auto const& rr : rows ){
        if( rr.row >= n ) continue;
        ++ovr_n;
        if( A.rownnz[rr.row] > 1 ) ++multi;
        size_t hits = 0;
        if( itc != _marchICColToNode.end() )
          for( size_t k = 0; k < A.rownnz[rr.row]; ++k ) if( itc->second.count( A.rowcol[rr.row][k] ) ) ++hits;
        if( hits != 1 ) ++ambig;
      }
    }
    _disp( 2 ) << "[march-override] imp=" << (int)options.INTERFACE.IMPOSITION
              << " weak_sat_terms=" << n_terms << " (valid_rows=" << n_valid << ")"
              << " sat_on_override_rows=" << A.ic_sat.size()
              << " | override_rows=" << ovr_n << " multi_col=" << multi
              << " ambiguous_pin(hits!=1)=" << ambig
              << " max|terminal|=" << tmax << "\n";
  }

  if( !_residual_core( var, inp, cst, r.data(), A, ic_override ) ) return rep;
  double norm = 0.0; for( double v : r ) norm = std::max( norm, std::fabs(v) );
  rep.initial_residual = norm;

  // Newton linear solve A dy = b through the configured sparse factorization backend.
  auto solve_linear = [&]( arma::sp_mat const& A, arma::vec const& b, arma::vec& xsol ) -> bool {
#if defined(ARMA_USE_SUPERLU)
    // Tier-0 is the *direct square* Newton solve and is sparse regardless of the
    // chosen rank-tolerant fallback: SOLVE_SPQR only changes the tier-1 step (QR
    // min-norm) reached when this solve is rejected, so keep tier-0 on sparse
    // SuperLU rather than densifying a ~1e4^2 Jacobian under SPQR.
    bool tier0_sparse = ( options.SOLVE.FACTORIZATION == Options::SOLVE_SUPERLU );
#if defined(CRONOS__WITH_SPQR)
    if( options.SOLVE.FACTORIZATION == Options::SOLVE_SPQR ) tier0_sparse = true;
#endif
    if( tier0_sparse ){
      try { return arma::spsolve( xsol, A, b, "superlu" ); }
      catch( ... ) { return false; }
    }
#else
    if( options.SOLVE.FACTORIZATION == Options::SOLVE_SUPERLU ){
      static bool warned = false;
      if( !warned ){
        _disp( 2 ) << "OCFESLV::solve ** SuperLU backend not linked (ARMA_USE_SUPERLU "
                     "undefined); using dense LAPACK factorization." << std::endl;
        warned = true;
      }
    }
#endif
    arma::mat Ad( A );   // densify the (square) scaled Jacobian
    return arma::solve( xsol, Ad, b, arma::solve_opts::fast + arma::solve_opts::no_approx );
  };

#if defined(CRONOS__WITH_SPQR)
  // Sparse rank-revealing QR min-2-norm least squares:  x = argmin ||x||  over
  // argmin ||A x - b||.  Solves the singular interface-augmented saddle directly
  // (no Js^T Js), returning the same min-norm step the dense tier does, sparsely.
  static bool spqr_rank_noted = false;   // rev118: companion one-shot, see below
  auto solve_lsq_spqr = [&]( arma::sp_mat const& A, arma::vec const& b, arma::vec& x ) -> bool {
    arma::sp_mat As = A;
    As.sync();                         // materialize CSC (values/row_indices/col_ptrs)

    cholmod_common cc;
    cholmod_l_start( &cc );
    // rev142b: CHOLMOD's own OpenMP regions honour this; guarded on the CHOLMOD major
    // version because nthreads_max does not exist before CHOLMOD 4 (SuiteSparse 6).
#if defined(CHOLMOD_MAIN_VERSION) && (CHOLMOD_MAIN_VERSION >= 4)
    if( int const _cap = t_ThreadCap::active() ) cc.nthreads_max = _cap;
#endif

    // --- wrap A as an unsymmetric real/double cholmod_sparse (CSC, long indices).
    // arma stores col_ptrs/row_indices as uword; SuiteSparse's C++ (template) API
    // is long-only, so copy/cast the index arrays to int64 (values shared by ptr).
    std::vector<int64_t> cp( As.n_cols + 1 ), ri( As.n_nonzero ? As.n_nonzero : 1 );
    for( size_t c = 0; c <= As.n_cols;  ++c ) cp[c] = (int64_t) As.col_ptrs[c];
    for( size_t k = 0; k < As.n_nonzero; ++k ) ri[k] = (int64_t) As.row_indices[k];

    cholmod_sparse Ac;  std::memset( &Ac, 0, sizeof(Ac) );
    Ac.nrow = As.n_rows; Ac.ncol = As.n_cols; Ac.nzmax = As.n_nonzero;
    Ac.p = cp.data(); Ac.i = ri.data(); Ac.x = const_cast<double*>( As.values );
    Ac.stype = 0;               // unsymmetric: full pattern stored
    Ac.itype = CHOLMOD_LONG;    // matches the int64 index arrays + cholmod_l_* routines
    Ac.xtype = CHOLMOD_REAL; Ac.dtype = CHOLMOD_DOUBLE;
    Ac.sorted = 1;              // arma CSC row indices are ascending within a column
    Ac.packed = 1;

    cholmod_dense Bc;  std::memset( &Bc, 0, sizeof(Bc) );
    Bc.nrow = b.n_elem; Bc.ncol = 1; Bc.nzmax = b.n_elem; Bc.d = b.n_elem;
    Bc.x = const_cast<double*>( b.memptr() );
    Bc.xtype = CHOLMOD_REAL; Bc.dtype = CHOLMOD_DOUBLE;

    // rev117: the rank tolerance is now settable, because it is a SUSPECT.
    //
    // SPQR_DEFAULT_TOL (-2) lets SPQR pick: 20*(m+n)*eps*max||A(:,j)||.  On a 12042-unknown
    // system that is ~1e-10.  A singular value ABOVE it is INVERTED (amplification 1/sigma),
    // not truncated -- so a direction at 1e-9 contributes 1e9 to the min-norm step, which is
    // the opposite of what "min-norm" suggests to a reader.
    //
    // MEASURED, OCFE_MBC5 constant-dP: QR (not a fall-through -- rev116's note does not
    // fire) produces a step 77x the solution vector, 77% of its energy on Drwc_Pw at one
    // inlet element, merit 1.1e+09 at full step against phi0 1.0e-05, and no alpha improves
    // on phi0.  Iterations 26-27 already crawl at alpha ~1e-7 with |r| pinned at 1.7756.
    // That is the signature of a direction just above the tolerance being inverted.
    //
    // SPQR semantics: tol >= 0 is an absolute column-norm threshold; -1 disables rank
    // detection (treat as full rank); -2 is the automatic default.  RAISING the tolerance
    // truncates more, which for a gauge-ish direction is what min-norm should do anyway.
    // Default unchanged (-2), so this is inert unless set.
    double spqr_tol = SPQR_DEFAULT_TOL;
    // (CRONOS_SPQR_TOL (retired 2026-10-07, WORKPLAN 3.B batch 2b): the automatic default always)
    // rev118: ONE-SHOT UNCONDITIONAL TRACE OF THE ACTUAL CALL.
    //
    // Not gated on any verbosity flag, and to std::cerr, deliberately.  rev117 made the
    // rank tolerance settable and CRONOS_SPQR_TOL produced output that was IDENTICAL
    // across -2, 1e-10, 1e-8 and 1e-6 -- down to best_alpha=1.1921e-07.  That is equally
    // consistent with "the tolerance does not matter here" and with "the hook never
    // reached the call", and this session has already had TWO hooks that silently failed
    // to fire (CRONOS_SOLVE_VERBOSE at the Options constructor, then at the adopt site).
    // A value read back at the call site is the only thing that separates them.
    //
    // Also reports the rank SPQR actually chose: if the tolerance is biting, the rank
    // moves with it; if the rank is identical across tolerances, it is not.
    static bool spqr_call_noted = false;
    if( !spqr_call_noted ){
      spqr_call_noted = true;
      _disp( 3 ) << "OCFESLV ** [spqr-call] A=" << As.n_rows << "x" << As.n_cols
                << "  tol_passed=" << std::scientific << std::setprecision(6) << spqr_tol
                << "  (SPQR_DEFAULT_TOL=" << (double)SPQR_DEFAULT_TOL << ")" << std::endl;
    }

    cholmod_dense* X = SuiteSparseQR_min2norm<double>(
                         SPQR_ORDERING_DEFAULT, spqr_tol, &Ac, &Bc, &cc );

    if( spqr_call_noted && !spqr_rank_noted ){
      spqr_rank_noted = true;
      _disp( 3 ) << "OCFESLV ** [spqr-call] estimated rank=" << cc.SPQR_istat[4]
                << " of " << As.n_cols
                << "   (rank < n means the tolerance TRUNCATED; rank == n means it did"
                   " not, and no tolerance in range will change this step)" << std::endl;
    }

    bool ok = false;
    if( X != nullptr && X->nrow == As.n_cols ){
      x.set_size( As.n_cols );
      double const* xv = static_cast<double const*>( X->x );
      for( size_t j = 0; j < As.n_cols; ++j ) x(j) = xv[j];
      ok = x.is_finite();
    }
    if( X ) cholmod_l_free_dense( &X, &cc );
    cholmod_l_finish( &cc );
    return ok;
  };
#endif // CRONOS__WITH_SPQR

  if( options.SOLVE.VERBOSE )
    _disp( 2 ) << std::left << std::setw(46) << "LM/Newton initial residual"
              << " n=" << std::setw(5) << n
              << " max|r|=" << std::scientific << std::setprecision(4) << norm << "\n";

  double lambda = options.SOLVE.LAMBDA0;
  bool   dense_minnorm = false;   // sticky: engage once the sparse damped step stalls on a
                                  // (near-)singular system (see stall detection below).
  bool   nt_dead = false;         // item 13: sticky -- once tier-0 Newton's factorization
                                  // fails (singular augmented saddle), skip re-attempting it
                                  // every iteration (the factorization is structurally doomed).
  bool   nt_disabled_noted = false;  // emit the "tier-0 disabled" note exactly once
  [[maybe_unused]] bool qr_fell_through_noted = false;  // rev116: emit the "QR unavailable" note exactly once
  bool   spectrum_audited  = false;  // AUDIT.SPECTRUM: run the audit at most once per solve()
  bool   step_fail_noted   = false;  // emit the line-search failure report exactly once per solve

  // ---- singular-value / null-space audit of the equilibrated Jacobian ---------------
  //
  // Supersedes the old "rank_deficiency=N" one-liner, which counted sigma_i <= 1e-12*sigma_max
  // on a matrix whose condition number is ~1e16 -- i.e. the cut sat INSIDE the noise, so the
  // number mixed a structural kernel with a conditioning tail and moved around for reasons that
  // are not structural.  Three changes:
  //
  //  (1) PRINT THE SPECTRUM, NOT A COUNT.  The smallest ~20 singular values plus the largest
  //      RELATIVE GAP among them.  A genuine kernel of dimension d shows a large gap after the
  //      d-th smallest value; a conditioning tail shows a smooth geometric decay and no gap.
  //      The reported deficiency is now the gap-based estimate, with the legacy fixed-threshold
  //      count printed alongside for continuity with earlier logs.
  //
  //  (2) RIGHT NULL SPACE.  The audit only ever reported LEFT null vectors (dependent rows).
  //      For a square deficient system the RIGHT null space (V) is what distinguishes an
  //      OVER-determined system (dependent rows, solution still pinned) from an UNDER-determined
  //      one (a direction the residual cannot see, which a min-norm step will fill arbitrarily --
  //      exactly the failure mode where max|r| reaches 1e-11 on a solution that is 80% wrong).
  //      Reported as a per-STATE energy share, which localises the kernel to a state family in
  //      one line, plus the dominant columns decoded to state/element/node/coordinate.
  //
  //  (3) EQUATION LEGEND.  Printed first, so "by eqn e6" is resolvable against the model.
  //
  // Debug path: dense SVD, small cases only.
  auto run_spectrum_audit = [&]( arma::sp_mat const& Jm,
                                std::vector<double> const& colw, char const* tag,
                                arma::vec const* pres = nullptr )
  {
    size_t const smax_dim = ( options.AUDIT.SPECTRUM_MAX > 0 ? options.AUDIT.SPECTRUM_MAX : 3000 );
    if( Jm.n_rows > smax_dim || Jm.n_cols > smax_dim ){
      _disp( 3 ) << "  [" << tag << "-audit] skipped: " << Jm.n_rows << "x" << Jm.n_cols
                << " exceeds the dense-SVD limit " << smax_dim
                << " (raise CRONOS_SPECTRUM_MAX, or reduce nel for a structural run)\n";
      return;
    }
    arma::mat const JD( Jm );
    arma::mat U, V; arma::vec sv;
    if( !arma::svd( U, sv, V, JD ) || !sv.n_elem ){
      _disp( 3 ) << "  [" << tag << "-audit] SVD failed\n";
      return;
    }

    print_equation_legend( std::cout );
    print_column_coverage_report( std::cout );
    print_tau_pivot_report( std::cout );

    double const smax = sv(0);
    arma::uword const ns = sv.n_elem;

    // Legacy fixed-threshold count, kept so old and new logs can be compared directly.
    size_t nz_fixed = 0;
    for( arma::uword i = 0; i < ns; ++i ) if( sv(i) <= 1.0e-12 * smax ) ++nz_fixed;

    _disp( 3 ) << "  [" << tag << "-audit] Js " << JD.n_rows << "x" << JD.n_cols
              << "  sigma_max=" << std::scientific << std::setprecision(4) << smax
              << "  sigma_min=" << sv(ns-1)
              << "  cond=" << ( sv(ns-1) > 0. ? smax/sv(ns-1) : INFINITY )
              << "  rank_deficiency[1e-12*smax]=" << nz_fixed << "\n";

    // Smallest values, ascending, with the ratio to the next-larger one.
    //
    // rev65 WINDOW FIX.  This was `min(20, ns)`, a hard literal, and the gap scan below never
    // looked past it.  Once the kernel dimension reached 20 the ENTIRE window sat INSIDE the
    // kernel: every value ~2.5e-16, every consecutive ratio 1.0, so the scan reported
    // "largest gap = 1.000e+00x ... kernel dim estimate 1  [WEAK GAP]" and the compatibility
    // verdict was suppressed -- while the real gap sat between s[k-1] and s[k], just outside.
    //
    // MEASURED (MMPDE13, --small, nel_xi=4, k = 3*(2*NEL_T-1)):
    //   NEL_T=2  k=9   verdict issued, rel=4.14e-04
    //   NEL_T=3  k=15  verdict issued, rel=6.29e-03
    //   NEL_T=4  k=21  SUPPRESSED  (rank_deficiency=21, correct; only the gap scan failed)
    //   NEL_T=5  k=27  SUPPRESSED  (rank_deficiency=27, correct)
    //
    // The consequence is worse than a missing line: rel becomes unmeasurable exactly where k
    // grows, which is where the question lives, and the suppression message advised "refine
    // the mesh and re-measure" -- the opposite of what helps, since refining raises k.  Any
    // historical rel reading was necessarily taken at k < 20.
    //
    // The window now always extends past the threshold-based deficiency count, which is
    // independent of the gap scan and was correct in every suppressed case above.
    arma::uword const nwant = ( nz_fixed > 0 )
                            ? (arma::uword)nz_fixed + 8 : (arma::uword)20;
    arma::uword const nshow = std::min<arma::uword>( std::max<arma::uword>( 20, nwant ), ns );
    bool const window_exhausted = ( nshow < (arma::uword)nz_fixed + 2 && ns >= nshow );
    _disp( 3 ) << "  [" << tag << "-audit] smallest " << nshow
              << " singular values (ascending, with ratio to the next):\n";
    for( arma::uword j = 0; j < nshow; ++j ){
      arma::uword const i = ns - 1 - j;                       // ascending
      double const ratio = ( i > 0 && sv(i) > 0. ) ? sv(i-1)/sv(i) : 0.;
      _disp( 3 ) << "      s[" << std::setw(3) << j << "] = " << std::scientific
                << std::setprecision(4) << sv(i)
                << "   next/this = " << ratio << "\n";
    }

    // Gap-based kernel-dimension estimate: the largest jump among the small end.
    arma::uword gap_at = 0; double gap_max = 0.;
    for( arma::uword j = 0; j + 1 < nshow; ++j ){
      arma::uword const i = ns - 1 - j;
      if( sv(i) <= 0. ) continue;
      double const ratio = sv(i-1)/sv(i);
      if( ratio > gap_max ){ gap_max = ratio; gap_at = j + 1; }
    }
    _disp( 3 ) << "  [" << tag << "-audit] largest gap = " << std::scientific
              << std::setprecision(3) << gap_max << "x after the " << gap_at
              << " smallest value(s)  =>  kernel dim estimate " << gap_at;
    if( gap_max < 1.0e3 && window_exhausted )
      _disp( 3 ) << "   [WINDOW EXHAUSTED -- the scan window (" << nshow << ") does not extend"
                   " past the deficiency count (" << nz_fixed << "), so no gap CAN be seen."
                   "  This is an instrument limit, not a property of the spectrum;"
                   " do NOT read it as a weak gap]";
    else if( gap_max < 1.0e3 )
      _disp( 3 ) << "   [WEAK GAP -- no clean separation; the deficiency count is a"
                   " conditioning artefact as much as a structural one]";
    _disp( 3 ) << "\n";

    size_t const nk = static_cast<size_t>( std::max<arma::uword>( gap_at, (arma::uword)nz_fixed ) );

    // ---- NUMERICALLY EMPTY COLUMNS --------------------------------------------------
    // The setup-time coverage report is STRUCTURAL: it reads _jacPatternEqn, so a column in
    // the pattern but zero in every row passes it.  Not hypothetical -- a continuity claim
    // whose coupling is zero mints a trace/tau column with no support at all, and such
    // columns have been observed `touched by 0 row(s)` while coverage read free=0.  Each is
    // an exact kernel direction (a unit basis vector), so scan for them numerically.
    {
      std::vector<arma::uword> empty_cols;
      for( arma::uword c = 0; c < JD.n_cols; ++c ){
        double cmax = 0.;
        for( arma::uword rr = 0; rr < JD.n_rows; ++rr )
          cmax = std::max( cmax, std::fabs( JD(rr,c) ) );
        if( cmax <= 0. ) empty_cols.push_back( c );
      }
      _disp( 3 ) << "  [" << tag << "-audit] numerically EMPTY columns = " << empty_cols.size();
      for( size_t k = 0; k < empty_cols.size() && k < 12; ++k )
        _disp( 3 ) << ( k ? "," : ":  " ) << "c" << empty_cols[k]
                  << "(" << decode_state_column( (size_t)empty_cols[k] ) << ")";
      if( empty_cols.size() > 12 ) _disp( 3 ) << " ...";
      if( !empty_cols.empty() )
        _disp( 3 ) << "\n      *** each is an exact kernel direction; a zero-coupling claim"
                     " mints a multiplier column with no support ***";
      _disp( 3 ) << "\n";
    }

    // ---- NUMERICALLY EMPTY ROWS ------------------------------------------------------
    // The mirror of the empty-column scan, and the one that tests the keep-explicit
    // contract.  The IC_STRONG hybrid promises, for each claim it keeps explicit, "an
    // explicit tau unknown (var[] column) PLUS a real continuity row", and
    // n_coll_eqn is incremented to make room for that row.  If the row is emitted with
    // content, the tau column cannot be empty -- the row would reference it.  The columns
    // ARE empty, so either the row is not emitted or it does not reference its tau.  An
    // empty ROW is a free residual slot that can never be driven to zero by any choice of
    // unknowns; paired with an empty column it is one unit of rank deficiency and one
    // unsatisfiable compatibility condition.  eqn=-1 marks a plan-owned interface/trace row.
    {
      std::vector<arma::uword> empty_rows;
      for( arma::uword rr = 0; rr < JD.n_rows; ++rr ){
        double rmx = 0.;
        for( arma::uword c = 0; c < JD.n_cols; ++c )
          rmx = std::max( rmx, std::fabs( JD(rr,c) ) );
        if( rmx <= 0. ) empty_rows.push_back( rr );
      }
      _disp( 3 ) << "  [" << tag << "-audit] numerically EMPTY rows = " << empty_rows.size();
      for( size_t k = 0; k < empty_rows.size() && k < 12; ++k ){
        int const e = ( (size_t)empty_rows[k] < _rowEqnNdx.size() )
                        ? _rowEqnNdx[(size_t)empty_rows[k]] : -99;
        _disp( 3 ) << ( k ? "," : ":  " ) << "r" << empty_rows[k] << "(eqn" << e << ")";
      }
      if( empty_rows.size() > 12 ) _disp( 3 ) << " ...";
      if( !empty_rows.empty() )
        _disp( 3 ) << "\n      *** an empty row is an unsatisfiable residual slot; if the count"
                     " matches the empty-column count, the keep-explicit contract"
                     " (tau column + continuity row) delivered the column but not the row ***";
      _disp( 3 ) << "\n";
    }

    // ---- TAU / MULTIPLIER BLOCK RANK ------------------------------------------------
    // Every trace/tau column (index >= trace_var_offset) is a Lagrange multiplier whose job
    // is to enforce one continuity claim.  If the tau BLOCK J[:,tau] is rank deficient, that
    // many multipliers are not independently supported by the rows, and each shortfall is one
    // kernel direction and one unsatisfiable compatibility condition.
    //
    // This replaces a pairwise "identical support, constant ratio" heuristic that was too
    // narrow.  It caught the IC_TRACE pair c603/c632 (same two rows, ratio -1) but MISSED a
    // second, distributed degeneracy in which {c543,c544,c551,c552,c564,c565,c567,c568} each
    // touch a DIFFERENT row pair yet are collectively dependent -- the null vectors combine
    // them in antisymmetric couples.  Measured shortfall: the heuristic reported 1 where the
    // kernel needed 5.  A rank test subsumes every signature, including orphans (a zero
    // column is rank-deficient by construction), and needs no pattern matching.
    //
    // The per-column support listing is retained for the columns the rank test implicates,
    // because "which rows carry this multiplier" is what localises the cause.
    if( _plan.trace_var_offset() < (size_t)JD.n_cols ){
      size_t const toff = _plan.trace_var_offset();
      arma::uword const ntau = JD.n_cols - (arma::uword)toff;
      arma::mat const T = JD.cols( (arma::uword)toff, JD.n_cols - 1 );

      arma::vec tsv;
      size_t n_orphan = 0;
      std::vector<arma::uword> orphan_cols;
      for( arma::uword c = 0; c < ntau; ++c ){
        double cmax = 0.;
        for( arma::uword rr = 0; rr < T.n_rows; ++rr )
          cmax = std::max( cmax, std::fabs( T(rr,c) ) );
        if( cmax <= 0. ){ ++n_orphan; orphan_cols.push_back( (arma::uword)toff + c ); }
      }

      _disp( 3 ) << "  [" << tag << "-audit] tau/multiplier block: " << ntau
                << " column(s) from offset " << toff << "\n";
      if( arma::svd( tsv, T ) && tsv.n_elem ){
        double const tmax = tsv(0);

        // Deficit = COUNT of null singular values; the gap only PLACES the tolerance.
        //
        // An earlier version returned the gap position itself, which under-counted by one
        // whenever the block is entirely null: with every sigma == 0 there is no ratio to
        // form at the smallest value, so the loop reported gap_at = ntau-1.  Measured on
        // IC_STRONG (5 orphan columns, sigma_max = 0): "RANK DEFICIT = 4 of 5" against an
        // orphan count of 5 and four tau-null vectors for five independent zero columns.
        // Handle the all-zero block explicitly and otherwise count against a gap-placed
        // tolerance, so the number is a rank deficit rather than an index.
        size_t deficit = 0;
        double gap_max = 0.; arma::uword gap_at = 0;
        if( tmax <= 0. ){
          deficit = (size_t)tsv.n_elem;          // entire block is null
          gap_max = INFINITY; gap_at = (arma::uword)deficit;
        }
        else{
          for( arma::uword k = 0; k + 1 < tsv.n_elem; ++k ){
            arma::uword const idx = tsv.n_elem - 1 - k;
            double const ratio = ( tsv(idx) > 0. ) ? tsv(idx-1)/tsv(idx) : INFINITY;
            if( ratio > gap_max ){ gap_max = ratio; gap_at = k + 1; }
          }
          if( gap_max > 1.0e3 ){
            // tolerance sits inside the gap; count everything below it
            arma::uword const above = tsv.n_elem - gap_at;
            double const tol = ( above < tsv.n_elem )
                             ? 0.5 * ( tsv(above) + tsv(above-1) ) : 0.;
            for( arma::uword k = 0; k < tsv.n_elem; ++k ) if( tsv(k) <= tol ) ++deficit;
          }
        }
        _disp( 3 ) << "      sigma_max=" << std::scientific << std::setprecision(4) << tmax
                  << "  sigma_min=" << tsv(tsv.n_elem-1)
                  << "  largest gap=" << std::setprecision(3) << gap_max
                  << "x after " << gap_at << "\n";
        _disp( 3 ) << "      RANK DEFICIT = " << deficit << " of " << ntau
                  << "   (orphan columns: " << n_orphan << ")\n";
        if( deficit && deficit != n_orphan && n_orphan )
          _disp( 3 ) << "      NOTE: deficit " << deficit << " vs orphan count " << n_orphan
                    << " -- orphans are a subset; the remainder are distributed dependencies\n";
        if( deficit )
          _disp( 3 ) << "      *** " << deficit << " multiplier(s) are not independently"
                       " supported by the rows; each is one kernel direction and one"
                       " compatibility condition on the residual.  Whether that condition is"
                       " VIOLATED is model-dependent and is MEASURED by the compatibility report"
                       " below (u_k^T r), not guessed: when the continuity rows are IMPLIED by"
                       " the physical rows the conditions hold identically, the kernel is a gauge"
                       " freedom and the primal answer is correct (measured: PDE2, 16 dependent"
                       " claims, all asserting continuity of a derivative auxiliary at a"
                       " transverse boundary face already pinned on both sides -- max|r|=3e-10)."
                       "  Compare deficit against the kernel dim: the remainder is not a tau"
                       " defect. ***\n";

        // Which tau columns carry the dependency: dominant entries of the tau-block's own
        // right-null vectors, decoded back to absolute column indices.
        arma::mat TU, TV; arma::vec tsv2;
        if( deficit && arma::svd( TU, tsv2, TV, T ) )
          for( size_t k = 0; k < deficit && k < 4; ++k ){
            arma::vec const vn = TV.col( TV.n_cols - 1 - k );
            arma::uvec const idx = arma::sort_index( arma::abs( vn ), "descend" );
            _disp( 3 ) << "      tau-null " << k << ":";
            for( arma::uword j2 = 0; j2 < 8 && j2 < idx.n_elem; ++j2 ){
              if( std::fabs( vn(idx(j2)) ) < 1.0e-6 ) break;
              _disp( 3 ) << " c" << ( (arma::uword)toff + idx(j2) ) << "="
                        << std::setprecision(3) << vn(idx(j2));
            }
            _disp( 3 ) << "\n";
          }
      }
      for( size_t k = 0; k < orphan_cols.size() && k < 12; ++k )
        _disp( 3 ) << "      ORPHAN c" << orphan_cols[k]
                  << ": touched by 0 rows -- multiplier enforces nothing\n";
    }

    // ---- LEFT null space: which ROWS are dependent (unchanged in substance) ----------
    for( size_t k = 0; k < nk && k < ns; ++k ){
      arma::vec const un = U.col( U.n_cols - 1 - k );
      arma::uvec const idx = arma::sort_index( arma::abs( un ), "descend" );
      size_t n_claim = 0, n_phys = 0; std::map<int,size_t> byeqn;
      for( arma::uword j = 0; j < idx.n_elem; ++j ){
        if( std::fabs( un(idx(j)) ) < 1.0e-6 ) break;
        int const e = ( (size_t)idx(j) < _rowEqnNdx.size() ) ? _rowEqnNdx[(size_t)idx(j)] : -99;
        if( e < 0 ) ++n_claim; else ++n_phys;
        ++byeqn[e];
      }
      _disp( 3 ) << "  [" << tag << "-audit] left-null " << k << " composition: claims=" << n_claim
                << " physical=" << n_phys << "  by eqn:";
      for( auto const& kv : byeqn ) _disp( 3 ) << " e" << kv.first << "x" << kv.second;
      _disp( 3 ) << "\n";
    }
    for( size_t k = 0; k < nk && k < 3; ++k ){
      arma::vec const un = U.col( U.n_cols - 1 - k );
      arma::uvec const idx = arma::sort_index( arma::abs( un ), "descend" );
      _disp( 3 ) << "  [" << tag << "-audit] left-null " << k << " dominant rows:";
      for( arma::uword j = 0; j < 6 && j < idx.n_elem; ++j ){
        arma::uword const rr = idx(j);
        _disp( 3 ) << " r" << rr << "(" << un(rr);
        if( (size_t)rr < _rowEqnNdx.size() ) _disp( 3 ) << ",eqn" << _rowEqnNdx[(size_t)rr];
        _disp( 3 ) << ")";
      }
      _disp( 3 ) << "\n";
      for( arma::uword j = 0; j < 3 && j < idx.n_elem; ++j ){
        arma::uword const rr = idx(j);
        if( std::fabs( un(rr) ) < 1.0e-6 ) break;
        double rmax = 0.;
        for( arma::uword c = 0; c < JD.n_cols; ++c )
          rmax = std::max( rmax, std::fabs( JD(rr,c) ) );
        double const tolr = 1.0e-10 * ( rmax > 0. ? rmax : 1.0 );
        _disp( 3 ) << "    r" << rr << " rmax=" << rmax << " nnz:";
        size_t cnt = 0, tot = 0;
        for( arma::uword c = 0; c < JD.n_cols; ++c )
          if( std::fabs( JD(rr,c) ) > tolr ){
            ++tot;
            if( cnt < 8 ){ _disp( 3 ) << " c" << c << "=" << JD(rr,c); ++cnt; }
          }
        _disp( 3 ) << "   (total nnz_rel=" << tot << ")\n";
      }
    }

    // ---- RIGHT null space: which UNKNOWNS the residual cannot see --------------------
    // Report the WHOLE kernel (bounded at 8), not the first three.  A kernel can mix
    // causes -- measured on the ALE mesh model: 6 directions, of which 5 were unsupported
    // tau columns and the 6th was something else that the old cap of 3 hid entirely.
    // Energy share per state answers "is the kernel in the mesh family or the solution
    // family" in one line; the decoded dominant columns then say which elements/nodes.
    // The COLUMN SCALING is undone first (dx = diag(cw)*dy) so the reported direction is
    // the one in the model's own variables, not the equilibrated ones.
    // Columns the kernel actually sits on, accumulated across the reported null
    // vectors.  The SVD returns an arbitrary basis of the null space, so any single
    // vector mixes the modes; the UNION of dominant columns is basis-independent and
    // is what the support dump below is keyed on.  Kept in a set so the same column
    // list appears in IC_WEAK, IC_TRACE and IC_STRONG runs and the three logs diff
    // line-for-line.
    // ---- CONSISTENCY: pair each left-null vector with the RESIDUAL -----------------
    //
    // rev8.  A rank-deficient square system has one compatibility condition per left-null
    // direction:  u_k^T r = 0.  Whether those conditions HOLD is what separates a benign
    // deficiency from a fatal one, and nothing in this arc has ever measured it -- the
    // separation was inferred from whether the solve happened to converge.
    //
    //   sum_k (u_k^T r)^2  IS the least-squares floor of the merit at this iterate:
    //       min_dx 0.5||r + J dx||^2 = 0.5 * sum_k (u_k^T r)^2
    // so this predicts the number the line search will stall at, and the prediction is
    // falsifiable against the [step-fail] merit-vs-alpha table.
    //
    //   ~0  => the continuity rows the dead multipliers belong to are IMPLIED by the
    //          physical rows.  The kernel is a gauge freedom: the multiplier is arbitrary,
    //          the primal is untouched, min-norm picks lambda=0 and the answer is correct.
    //          (PDE2's 16: every dependent claim asserts continuity of a derivative
    //          auxiliary at a transverse boundary face where a Dirichlet BC already pins
    //          both sides -- the same predicate the emission-side "pinned pointwise on both
    //          sides" suppressor applies, but reaching one side through a LINK row.)
    //
    //   >0  => the claims genuinely extend rank and the dead column makes the system
    //          INCONSISTENT.  No exact solution exists and QR stops at the LS optimum.
    //
    // Cost: one dot product per null direction against a vector we already hold.
    if( pres && (arma::uword)pres->n_elem == JD.n_rows && nk ){
      // rev9 GATE.  The compatibility test is only meaningful on directions that are
      // NUMERICALLY NULL.  rev8 gated it on nk -- the gap-based kernel estimate -- which
      // on a weak gap includes ill-conditioned but perfectly solvable directions.  Those
      // have u_k^T r != 0 legitimately, and the resulting "floor" is fictitious.
      //
      // Measured failure of rev8: PDE2 under IC_WEAK, sigma_min = 2.13e-06 against
      // sigma_max = 3.51 (the audit's own "WEAK GAP" annotation fires), reported
      // INCONSISTENT with a floor of 1.43e-09 -- and then converged to max|r| = 1.7e-10,
      // merit 5.8e-19.  The verdict was wrong; the gate was the bug.
      //
      // A direction now enters the sum only if sigma_k <= 1e-10 * sigma_max, and the
      // verdict is suppressed entirely when the gap is weak (gap_max < 1e3), because on
      // such a spectrum "null" is not a well-posed classification at all.
      double const s_null_tol = 1.0e-10 * smax;
      double s2 = 0.;
      size_t n_null = 0;
      _disp( 3 ) << "  [" << tag << "-audit] compatibility u_k^T r (per NULL direction; "
                   "sigma_k <= " << std::scientific << std::setprecision(2) << s_null_tol
                << " required):";
      for( size_t k = 0; k < nk && k < ns; ++k ){
        double const sk = sv( (arma::uword)( ns - 1 - k ) );
        if( sk > s_null_tol ){
          if( k < 12 ) _disp( 3 ) << "  u" << k << "[sigma=" << std::scientific
                                 << std::setprecision(2) << sk << " NOT NULL, skipped]";
          continue;
        }
        arma::vec const un = U.col( U.n_cols - 1 - k );
        double const ur = arma::dot( un, *pres );
        s2 += ur * ur;
        ++n_null;
        if( k < 12 ) _disp( 3 ) << "  u" << k << "^Tr=" << std::scientific
                               << std::setprecision(3) << ur
                               << "(sigma=" << std::setprecision(2) << sk << ")";
      }
      if( nk > 12 ) _disp( 3 ) << "  ...";
      _disp( 3 ) << "\n";

      double const rn = arma::norm( *pres );
      if( !n_null ){
        _disp( 3 ) << "  [" << tag << "-audit] no direction among the " << nk
                  << " smallest is numerically null: this is a CONDITIONING tail, not a"
                     " kernel, and there are no compatibility conditions to test.\n";
      }
      else if( gap_max < 1.0e3 && window_exhausted ){
        _disp( 3 ) << "  [" << tag << "-audit] verdict suppressed because the gap scan WINDOW"
                     " was exhausted (" << nshow << " values vs deficiency " << nz_fixed
                  << "), not because the spectrum lacks a gap.  rev65 widens the window; if"
                     " this still appears, the deficiency count itself is suspect.\n";
      }
      else if( gap_max < 1.0e3 ){
        _disp( 3 ) << "  [" << tag << "-audit] " << n_null << " of " << nk
                  << " direction(s) below the null tolerance, but the spectral gap is WEAK ("
                  << std::scientific << std::setprecision(3) << gap_max
                  << "x): 'null' is not a well-posed classification on this spectrum, so no"
                     " compatibility verdict is issued.  Improve the conditioning and"
                     " re-measure.  NOTE: refining the mesh RAISES k on interface-kernel"
                     " models, so it is not the remedy here.\n";
      }
      else{
        // ---------------------------------------------------------------------------
        // rev111: DETERMINACY, measured BEFORE the consistency verdict is worded.
        // ---------------------------------------------------------------------------
        // CONSISTENCY and DETERMINACY are different properties and this block used to
        // report only the first while asserting the second.  u_k^T r = 0 says the system
        // is SOLVABLE.  It says nothing about whether the solution is UNIQUE -- that
        // depends on whether the RIGHT-null space has state content.
        //
        // MEASURED, OCFE_scalar_nest4 under Z-support promotion: the CONSISTENT branch
        // fired at rel = 3.0e-16 and asserted "the primal solution is unaffected", three
        // lines above its own report that 64.8% of the right-null energy sat on Dz_a(z).
        // The answer was wrong by 4-20%.  The audit had the number and did not use it.
        //
        // The premise of the old wording -- "the DEAD MULTIPLIER COLUMN(S) are a gauge
        // freedom" -- holds only when the null space lies IN the multiplier block.  Where
        // it lies in the states there are no dead multiplier columns (that run reported
        // orphan columns = 0) and the sentence describes a configuration that is not the
        // one in front of it.
        //
        // So: measure the state share of every null direction first, and let it choose
        // the wording.  Same vectors, same weighting as the per-direction report below.
        double state_share_max = 0.;
        size_t state_dir_worst = 0;
        for( size_t kk = 0; kk < nk; ++kk ){
          arma::vec vk = V.col( V.n_cols - 1 - kk );
          for( arma::uword j = 0; j < vk.n_elem && (size_t)j < colw.size(); ++j )
            vk(j) *= colw[(size_t)j];
          double const nk2 = arma::norm( vk );
          if( nk2 <= 0. ) continue;
          vk /= nk2;
          double e_state = 0.;
          for( arma::uword j = 0; j < vk.n_elem; ++j ){
            double const e = vk(j)*vk(j);
            if( e <= 0. ) continue;
            // a STATE column decodes with a "[" index suffix; multiplier/trace columns
            // do not.  Same discriminator the per-direction report uses.
            if( decode_state_column( (size_t)j ).find( '[' ) != std::string::npos )
              e_state += e;
          }
          if( e_state > state_share_max ){ state_share_max = e_state; state_dir_worst = kk; }
        }
        bool const primal_determined = ( state_share_max <= 1.0e-2 );   // 1% of the energy

        _disp( 3 ) << "  [" << tag << "-audit] least-squares FLOOR at this iterate = "
                  << std::scientific << std::setprecision(4) << ( 0.5 * s2 )
                  << "   over " << n_null << " null direction(s)"
                  << "   (||r||=" << rn << ", phi=" << ( 0.5 * rn * rn ) << ")\n";
        double const rel = ( rn > 0. ? std::sqrt( s2 ) / rn : 0. );
        if( rel <= 1e-10 && primal_determined )
          _disp( 3 ) << "      *** CONSISTENT and PRIMAL-DETERMINED (rel=" << rel
                    << ", max state share of any null direction = "
                    << std::fixed << std::setprecision(1) << ( 100.*state_share_max ) << "%"
                    << std::scientific << "): every compatibility condition holds AND the"
                       " null space lies in the multiplier block, so the dead multiplier"
                       " column(s) are a GAUGE freedom -- the continuity rows they belong to"
                       " are IMPLIED by the physical rows, the primal solution is unaffected,"
                       " and min-norm returns the correct answer.  The right fix is"
                       " emission-side (do not emit an implied claim), not the multiplier"
                       " block.  (Measured: PDE2 IC_STRONG, 16 dead columns, rel = 2.2e-16 at"
                       " iterate 0 against ||r|| = 1.38.) ***\n";
        else if( rel <= 1e-10 )
          _disp( 3 ) << "      *** CONSISTENT but PRIMAL UNDETERMINED (rel=" << rel
                    << "): the compatibility conditions hold, so the system is SOLVABLE -- but"
                       " " << std::fixed << std::setprecision(1) << ( 100.*state_share_max )
                    << "% of right-null direction " << state_dir_worst << "'s energy is in the"
                       " STATES" << std::scientific << ", so the solution is NOT UNIQUE.  The"
                       " solver returns a point on a solution manifold and WHICH point depends"
                       " on the step sequence; the answer can be wrong while the residual is at"
                       " round-off.  This is NOT a gauge freedom and min-norm does NOT return"
                       " the correct answer.  See the per-direction energy below for the state"
                       " carrying it.  (Measured: OCFE_scalar_nest4 under Z-support promotion,"
                       " rel = 3.0e-16, 91.4% state energy, converged |r| = 8.4e-15 to an answer"
                       " wrong by 4-20%.) ***\n";
        else if( rel >= 0.5 )
          std::cout << "      *** INCONSISTENT and FLOOR-BOUND (rel=" << rel << "): the"
                       " compatibility conditions are VIOLATED and the residual is essentially"
                       " ALL compatibility defect, so no exact solution exists, the solver has"
                       " already extracted everything the equations permit, and the stall is"
                       " the floor.  Fixing the interface plan is what removes it.  (Measured:"
                       " MMPDE9 IC_STRONG floor 1.5630e-06 vs stalled merit 1.5630e-06;"
                       " IC_TRACE 1.0596e-06 vs 1.0596e-06; MMPDE10 IC_STRONG at Pe=100"
                       " rel=1.0000.) ***\n";
        else
          // rev64 VERDICT SPLIT.  rev63 said "the line search will stall at merit >= the floor
          // above" for EVERY rel > 1e-10.  True, and useless when rel << 1: MMPDE13 at Pe=10
          // measured rel = 6.3e-03 with floor 1.58e-09 against a stalled merit of 2.51e-07 --
          // the solver stopped 158x ABOVE the floor, so the kernel accounted for 0.6% of the
          // stall and the rest was the (c) linearisation breakdown.  Both quoted examples in
          // the old text were rel = 1.0000 cases, which reinforced reading every INCONSISTENT
          // verdict as floor-bound.  The distinction decides whether fixing the interface plan
          // is worth anything on this model, so the report now makes it instead of the reader.
          std::cout << "      *** INCONSISTENT but NOT FLOOR-BOUND (rel=" << rel << "): the"
                       " compatibility conditions are VIOLATED, so no exact solution exists --"
                       " but only rel^2 = " << ( rel * rel ) << " of the residual ENERGY lies in"
                       " the left-null space.  If the solve stalls well above the floor printed"
                       " above, the kernel is NOT what stopped it and fixing the interface plan"
                       " buys at most this fraction; look to the step-failure attribution"
                       " instead.  (Measured: MMPDE13 IC_STRONG at Pe=10, rel=6.3e-03, floor"
                       " 1.58e-09 vs stalled merit 2.51e-07 -- a factor 158.) ***\n";

        // AFTER the chain, never inside it.  An earlier revision of this edit put a bare
        // `if` between the brace-less `else if( rel >= 0.5 )` and its `else`, which
        // re-bound that `else` to the new `if` -- so two mutually exclusive verdicts
        // printed for the same rel.  Braces balanced and the preflight passed; only
        // running it showed the double print.  Hence the explicit braces below.
        if( rel > 1e-10 && !primal_determined ){
          _disp( 3 ) << "      *** ALSO PRIMAL UNDETERMINED: "
                    << std::fixed << std::setprecision(1) << ( 100.*state_share_max )
                    << "% of right-null direction " << state_dir_worst << "'s energy is in the"
                       " STATES" << std::scientific << ".  Whatever the solver reaches is one"
                       " point of many; do not read the primal until this is removed.\n";
        }

        // ---------------------------------------------------------------------------
        // rev119: TOLERANCE REACHABILITY -- a SEPARATE question from the chain above.
        // ---------------------------------------------------------------------------
        // The FLOOR-BOUND / NOT FLOOR-BOUND split (rev64) answers: does the kernel account
        // for MOST OF THE STALL?  `rel` is the fraction of residual ENERGY in the left-null
        // space, and on MMPDE13 at Pe=10 it correctly said no -- rel = 6.3e-03 with the
        // solver stopping 158x above the floor.  That verdict is right for what it measures.
        //
        // It does NOT answer: can SOLVE.RES_TOL be reached at all?  Those come apart when the
        // floor is small in absolute terms but larger than the tolerance.
        //
        // MEASURED, OCFE_DAE5 --nomarch IC_STRONG.  rel = 1.79e-07 lands in the middle band,
        // so the chain prints "NOT FLOOR-BOUND" -- and the kernel really is only rel^2 =
        // 3.2e-14 of the energy.  But SPQR reaches the least-squares optimum in ONE step at
        // max|r| = 1.94e-07, which equals rel to 8%, against SOLVE.RES_TOL = 1e-09.  The solve
        // is AT its floor and the tolerance is unreachable, at nel_t 4 and 8.  At nel_t 16 the
        // floor has fallen to 2.32e-10, below the tolerance, and the same driver converges --
        // the verdict flips on refinement, not on anything the solver did.
        //
        // I read "NOT FLOOR-BOUND" as "the floor is not what stops it" and used it to rule out
        // a floor-aware convergence test.  The chain never made that claim; the reachability
        // statement below is the one that does, and it is decidable AT ITERATE 0 because
        // SOLVE.RES_TOL is known.
        //
        // BOUND: max|r| >= ||r||_2 / sqrt(n) and ||r||_2 >= sqrt(s2) at the optimum, so
        // sqrt(s2)/sqrt(n) is a RIGOROUS lower bound on the max-norm the driver tests.
        // Conservative on purpose -- it fires only when the tolerance is provably out of reach.
        {
          double const floor_2norm  = std::sqrt( s2 );
          double const floor_maxabs = floor_2norm / std::sqrt( (double) JD.n_rows );
          if( floor_maxabs > options.SOLVE.RES_TOL )
            _disp( 3 ) << "      *** TOLERANCE UNREACHABLE: least-squares floor is "
                      << std::scientific << std::setprecision(4) << floor_2norm
                      << " in ||r||_2, so max|r| >= " << floor_maxabs
                      << " against SOLVE.RES_TOL = " << options.SOLVE.RES_TOL
                      << ".  NO step can satisfy the convergence test -- a stall at this level"
                         " is the TOLERANCE, not the solver and not the step quality.  Either"
                         " the claim set is over-determined (the floor is real and scales with"
                         " the discretisation, so it converges away under refinement) or the"
                         " tolerance should be relative to it.  (Measured: OCFE_DAE5 --nomarch,"
                         " floor 1.94e-07 / 6.94e-09 / 2.32e-10 at nel_t 4/8/16 against a fixed"
                         " 1e-09 -- conv=no, no, YES.) ***\n";
        }
      }
    }

    std::set<arma::uword> kernel_cols;

    size_t const nk_show = std::min<size_t>( nk, 8 );
    if( nk > nk_show )
      _disp( 3 ) << "  [" << tag << "-audit] kernel dim " << nk
                << "; showing the first " << nk_show << " right-null vector(s)\n";
    for( size_t k = 0; k < nk_show; ++k ){
      arma::vec vn = V.col( V.n_cols - 1 - k );
      for( arma::uword j = 0; j < vn.n_elem && (size_t)j < colw.size(); ++j )
        vn(j) *= colw[(size_t)j];
      double const nrm = arma::norm( vn );
      if( nrm > 0. ) vn /= nrm;

      std::map<std::string,double> energy;                   // state name -> sum of squares
      double e_other = 0.;
      for( arma::uword j = 0; j < vn.n_elem; ++j ){
        double const e = vn(j)*vn(j);
        if( e <= 0. ) continue;
        std::string const d = decode_state_column( (size_t)j );
        size_t const br = d.find( '[' );
        if( br == std::string::npos ){ e_other += e; continue; }
        energy[ d.substr( 0, br ) ] += e;
      }
      _disp( 3 ) << "  [" << tag << "-audit] right-null " << k << " energy by state:";
      std::vector<std::pair<double,std::string>> ord;
      for( auto const& kv : energy ) ord.emplace_back( kv.second, kv.first );
      std::sort( ord.begin(), ord.end(), []( auto const& a, auto const& b ){ return a.first > b.first; } );
      for( auto const& p : ord )
        if( p.first > 1.0e-6 )
          _disp( 3 ) << " " << p.second << "=" << std::fixed << std::setprecision(1)
                    << 100.*p.first << "%";
      if( e_other > 1.0e-6 )
        _disp( 3 ) << " trace/tau=" << std::fixed << std::setprecision(1) << 100.*e_other << "%";
      _disp( 3 ) << std::scientific << "\n";

      arma::uvec const vidx = arma::sort_index( arma::abs( vn ), "descend" );
      _disp( 3 ) << "  [" << tag << "-audit] right-null " << k << " dominant columns:\n";
      for( arma::uword j = 0; j < 8 && j < vidx.n_elem; ++j ){
        arma::uword const cc2 = vidx(j);
        if( std::fabs( vn(cc2) ) < 1.0e-6 ) break;
        if( j < 4 ) kernel_cols.insert( cc2 );        // top few per mode -> support dump
        _disp( 3 ) << "      c" << std::setw(5) << cc2 << " = " << std::scientific
                  << std::setprecision(4) << vn(cc2)
                  << "   " << decode_state_column( (size_t)cc2 ) << "\n";
      }
      // A kernel direction that leaves the residual untouched: report how much of the
      // residual it actually moves, as the sanity check on the whole diagnosis.
      arma::vec const Jv = JD * V.col( V.n_cols - 1 - k );
      _disp( 3 ) << "      ||J*v|| = " << std::scientific << std::setprecision(3)
                << arma::norm( Jv ) << "  (||v||=1, scaled system)\n";
    }

    // ---- ROW SUPPORT of the kernel-implicated columns -------------------------------
    //
    // Read the matrix instead of modelling it.  Three mechanisms for this defect were
    // proposed from the assembly's structure and all three were refuted by measurement
    // (corner dedup: suppressed=0; missing claim: absent under IC_WEAK too; consumed
    // pivot row: IC_TRACE has no pivots and the same kernel; structural coverage: 0 free
    // columns in every mode).  What every one of those had in common was reasoning
    // FORWARD from what the plan should do.  This block does the opposite: for each
    // column the kernel sits on, it lists every row that touches it, with the coefficient
    // and that row's own nonzero pattern.
    //
    // Run in IC_WEAK, IC_TRACE and IC_STRONG and diff.  Coverage is identical across the
    // three, so the rows are the same rows; the dependency has to be visible in the
    // coefficients.  If the two rows covering a column differ only by a SAT term present
    // under IC_WEAK and absent under exact imposition, that is the mechanism, read off
    // directly.  If they look independent, the model of the assembly is wrong somewhere
    // else and this says so too.
    if( !kernel_cols.empty() ){
      _disp( 3 ) << "  [" << tag << "-audit] row support of " << kernel_cols.size()
                << " kernel column(s)  (imposition="
                << ( options.INTERFACE.IMPOSITION == Options::IC_WEAK   ? "IC_WEAK"
                   : options.INTERFACE.IMPOSITION == Options::IC_TRACE  ? "IC_TRACE"
                   : options.INTERFACE.IMPOSITION == Options::IC_STRONG ? "IC_STRONG" : "?" )
                << ")\n";

      for( arma::uword cc2 : kernel_cols ){
        // Relative threshold per column so a genuinely small-but-real entry is not cut
        // while numerical dust is.
        double cmax = 0.;
        for( arma::uword rr = 0; rr < JD.n_rows; ++rr )
          cmax = std::max( cmax, std::fabs( JD(rr,cc2) ) );
        double const tolc = 1.0e-10 * ( cmax > 0. ? cmax : 1.0 );

        std::vector<arma::uword> rows;
        for( arma::uword rr = 0; rr < JD.n_rows; ++rr )
          if( std::fabs( JD(rr,cc2) ) > tolc ) rows.push_back( rr );

        _disp( 3 ) << "    c" << std::setw(5) << cc2 << "  " << decode_state_column( (size_t)cc2 )
                  << "   touched by " << rows.size() << " row(s), colmax=" << std::scientific
                  << std::setprecision(3) << cmax << "\n";

        for( size_t m = 0; m < rows.size() && m < 6; ++m ){
          arma::uword const rr = rows[m];
          int const e = ( (size_t)rr < _rowEqnNdx.size() ) ? _rowEqnNdx[(size_t)rr] : -99;
          double rmax = 0.;
          for( arma::uword c3 = 0; c3 < JD.n_cols; ++c3 )
            rmax = std::max( rmax, std::fabs( JD(rr,c3) ) );
          double const tolr = 1.0e-10 * ( rmax > 0. ? rmax : 1.0 );
          size_t tot = 0;
          for( arma::uword c3 = 0; c3 < JD.n_cols; ++c3 )
            if( std::fabs( JD(rr,c3) ) > tolr ) ++tot;

          _disp( 3 ) << "        r" << std::setw(5) << rr << " eqn" << std::setw(3) << e
                    << "  coef=" << std::scientific << std::setprecision(4) << JD(rr,cc2)
                    << "  rmax=" << rmax << "  nnz=" << tot << "  cols:";
          size_t shown = 0;
          for( arma::uword c3 = 0; c3 < JD.n_cols && shown < 8; ++c3 )
            if( std::fabs( JD(rr,c3) ) > tolr ){
              _disp( 3 ) << " c" << c3 << "=" << std::setprecision(3) << JD(rr,c3);
              ++shown;
            }
          if( tot > shown ) _disp( 3 ) << " ...";
          _disp( 3 ) << "\n";
        }
        if( rows.size() > 6 )
          _disp( 3 ) << "        ... " << ( rows.size() - 6 ) << " more row(s)\n";
      }
    }
  };

  // Per-iteration cost breakdown (printed at the end when SOLVE.VERBOSE): which
  // part of each Newton/LM iteration dominates -- the AD Jacobian (deriv), the
  // equilibration, the sparse assembly, the linear step solve, or the line
  // search residual evaluations.  step_solve is split by tier (item 13) to
  // separate the wasted tier-0 factorization from the LM JtJ solve.
  double _ts_deriv=0.0, _ts_equil=0.0, _ts_asm=0.0, _ts_solve=0.0, _ts_ls=0.0;
  double _ts_nt=0.0, _ts_lm=0.0, _ts_mn=0.0;
  double _ts_lmform=0.0, _ts_lmsolve=0.0;  // item 13: LM internal split -- JtJ product+assembly vs SuperLU factor/solve
  auto _swnow = []{ return std::chrono::steady_clock::now(); };
  auto _swms  = []( std::chrono::steady_clock::time_point a,
                    std::chrono::steady_clock::time_point b )
                { return std::chrono::duration<double,std::milli>( b - a ).count(); };

#if defined(CRONOS__WITH_SPQR)
  // item 13: an interface-augmented system (appended trace/tau columns) has a
  // structurally singular Jacobian, so the tier-0 SuperLU Newton factorization
  // would fail every iteration.  Under SPQR the QR min-norm tier resolves it
  // directly, so skip the doomed tier-0 attempt a priori (the sticky nt_dead
  // skip never helps here -- the QR step converges in one iteration, leaving no
  // later iteration for the skip to benefit).
  // rev143: the skip fired on _nTraceVar>0 alone, which assumes augmented => singular.
  // IC_TRACE is the measured counterexample: full rank (PDE5 mono rank 27272/27272),
  // and the accidental no-SPQR build showed plain Newton succeeding there (15.0 s vs
  // 20.2 s QR, cross-binary and suggestive only -- this switch exists to measure it
  // pinned, single-binary).  Policy per _knob_SKIP_NT_ON_TRACE() / env override.
  static int const kSkipNTEnv = -1;   // CRONOS_SKIP_NT_ON_TRACE retired 2026-10-07: _knob_SKIP_NT_ON_TRACE() decides
  int const skipnt_mode = ( kSkipNTEnv >= 0 ? kSkipNTEnv : _knob_SKIP_NT_ON_TRACE() );
  bool spqr_skip_nt =
        ( options.SOLVE.FACTORIZATION == Options::SOLVE_SPQR && _nTraceVar > 0 );
  if( spqr_skip_nt && skipnt_mode == Options::NT_SKIP_NEVER )
    spqr_skip_nt = false;
  else if( spqr_skip_nt && skipnt_mode == Options::NT_SKIP_AUTO ){
    t_Determinacy const* D = _nt_skip_verdict( _inMarch, _inMarch ? _marchWindow : 0 );
    bool const affirm_full = D && D->ran && D->verdict != t_Determinacy::SKIPPED
                          && D->def_right == 0;
    if( affirm_full ){
      spqr_skip_nt = false;
      _disp( 2 ) << "  note: skip-nt=auto: audit "
                << ( _inMarch ? "win=" + std::to_string( _marchWindow ) + " " : "" )
                << "defR=0 -> attempting tier-0 Newton on the augmented system\n";
    }
    else if( options.SOLVE.VERBOSE )
      _disp( 2 ) << "  note: skip-nt=auto: "
                << ( D && D->ran && D->verdict != t_Determinacy::SKIPPED
                     ? "audit defR=" + std::to_string( D->def_right )
                     : std::string( "no usable audit verdict" ) )
                << " -> keeping the a-priori skip\n";
  }
  if( options.SOLVE.VERBOSE && spqr_skip_nt )
    _disp( 2 ) << "  note: interface-augmented system (nTrace=" << _nTraceVar
              << "); tier-0 Newton skipped, using QR min-norm steps\n";
#else
  bool const spqr_skip_nt = false;
#endif

  for( int it=0; it < options.SOLVE.MAX_ITER && norm > options.SOLVE.RES_TOL; ++it ){
    if( !_assemble_core( var, inp, cst, A, options.SOLVE.EQUILIBRATE, ic_override, r.data() ) ) return rep;
    _ts_deriv += A.ms_deriv;
    _ts_asm   += A.ms_asm;
    arma::sp_mat&        Js = A.J;     // diag(rw) * J * diag(cw), n x m (state cols)
    std::vector<double>& rw = A.rw;
    std::vector<double>& cw = A.cw;
    arma::vec&           rs = A.r;     // row-scaled residual
    double const phi = 0.5*arma::dot( rs, rs );

    // AUDIT.SPECTRUM: audit the FIRST assembled Jacobian regardless of whether tier-0
    // Newton later fails.  The historical trigger fired only on a singular factorization,
    // so a run that converged to a wrong answer through an unconstrained kernel was never
    // audited -- which is precisely the comparison that has to be made.
    if( options.AUDIT.SPECTRUM && !spectrum_audited ){
      spectrum_audited = true;
      run_spectrum_audit( Js, cw, "spectrum", &rs );
    }

#ifdef MC__OCFESLV_CONDITIONING_PROBE
    // [march-cond] sparse smallest-singular-value + null-vector localiser for a marching window.
    // OPT-IN behind MC__OCFESLV_CONDITIONING_PROBE (compile flag): each call runs ~12 sparse solves
    // (inverse iteration on Js^T Js), so it is NOT compiled into normal builds.  Runs once per window
    // (first iterate; the linear heat-window Jacobian is iteration-independent).  Reports the smallest
    // singular pair of the SAME equilibrated Js the solver factorizes: sig_min quantifies the
    // conditioning and the per-state max|null-vector| names the DOFs of the offending mode -- the
    // diagnostic that localised the IC_WEAK marching rank-deficiency to a too-soft interface SAT.
    if( _inMarch && it == 0 && n == m && m > 0 ){
      arma::sp_mat const JsT = Js.t();
      arma::arma_rng::set_seed( 1u );
      arma::vec xv( (arma::uword)m, arma::fill::randn );
      double nx = arma::norm( xv ); if( nx > 0. ) xv /= nx;
      bool ok_probe = ( nx > 0. );
      for( int p = 0; p < 6 && ok_probe; ++p ){          // (Js^T Js)^{-1} x = Js^{-1} (Js^{-T} x)
        arma::vec w, y;
        { OCEnvScopedFdSilence const silence;
          ok_probe = solve_linear( JsT, xv, w ) && w.is_finite()
                  && solve_linear( Js,  w,  y ) && y.is_finite(); }
        if( !ok_probe ) break;
        double const ny = arma::norm( y ); if( !( ny > 0. ) ){ ok_probe = false; break; }
        xv = y / ny;
      }
      if( ok_probe ){
        double const sigmin = arma::norm( Js * xv );     // ||xv||=1 -> ||Js xv|| = sig_min
        _disp( 2 ) << "  [march-cond] sig_min(Js)=" << std::scientific << std::setprecision(3) << sigmin;
        for( auto const& [st, dom] : _mVar ){
          auto itb = _blkVar.find( st ); if( itb == _blkVar.end() ) continue;
          size_t const bpos = std::get<0>( itb->second ), btot = std::get<1>( itb->second );
          size_t nblk = 1;
          for( auto const& d : dom ){ auto itd = _mDom.find( d ); if( itd != _mDom.end() ) nblk *= itd->second.n_elem; }
          size_t const hi = std::min<size_t>( bpos + btot*nblk, (size_t)m );
          double mx = 0.0; size_t imx = bpos;
          for( size_t i = bpos; i < hi; ++i ){ double const a = std::fabs( xv((arma::uword)i) ); if( a > mx ){ mx = a; imx = i; } }
          _disp( 2 ) << "  " << st << "=" << mx << "@" << imx;
        }
        _disp( 2 ) << "\n";
      }
      else _disp( 2 ) << "  [march-cond] inverse-iteration probe failed (Js singular/solve failed)\n";
    }
#endif // MC__OCFESLV_CONDITIONING_PROBE

#ifdef MC__OCFESLV_CONDITIONING_PROBE
    // ----- conditioning probe (OPT-IN, PRINT-ONLY) ----------------------------
    // Dense SVD of the row/col-equilibrated Jacobian Js the solver factorizes.
    // Recovered 2026-07-02 from the 2026-07-01 transcript (lost to a fork).
    // A failing run's trajectory shows whether kappa blows up (conditioning-
    // limited basin) or stays modest while merit stalls (guess-distance limited).
    // Densify only at test sizes; skip on large systems.
    if( n == m && m <= 3000 ){
      arma::vec sv;
      bool const ok = arma::svd( sv, arma::mat( Js ) );
      double const smax = ( ok && sv.n_elem ) ? sv.max() : 0.0;
      double const smin = ( ok && sv.n_elem ) ? sv.min() : 0.0;
      double const kap  = ( smin > 0.0 ) ? smax/smin : arma::datum::inf;
      _disp( 3 ) << "  [cond] it=" << it << std::scientific << std::setprecision(3)
                << " sig_max=" << smax << " sig_min=" << smin
                << " kappa=" << kap << " (JtJ kappa^2=" << ( smin>0.0 ? kap*kap : arma::datum::inf ) << ")"
                << " merit=" << phi << " |r|=" << norm << "\n";
    }
#endif

#ifdef MC__OCFESLV_TIER0_PROBE
    // ----- tier-0 / saddle diagnostic (item 13, OPT-IN, PRINT-ONLY) -----------
    // For an interface-augmented system (nTrace>0, IC_TRACE/IC_STRONG) the tier-0
    // direct Newton factorization is normally skipped (SPQR) or fails (SuperLU),
    // forcing the expensive QR / JtJ-LM tier every iteration.  This SHADOW probe
    // attempts the tier-0 solve WITHOUT using its step, to decide WHY tier-0 is
    // unavailable and hence which remedy applies:
    //   (a) solve_linear FLAG -- does the direct factorization of the augmented
    //       saddle succeed?  Consistently FAIL => structurally singular saddle =>
    //       KKT/saddle-point direct solve (Schur on the tau block).
    //   (b) line-search MERIT of the Newton step (alpha=1) + Armijo verdict -- if
    //       solve_linear OKs a step that REDUCES merit yet the solver still uses
    //       LM/QR, a Newton step is being thrown away => relaxed Newton acceptance.
    //   (c) tau-column vs state-column equilibration SPREAD -- if the saddle only
    //       goes singular because the tau columns are badly scaled => better
    //       tau-column equilibration (numerical, not structural).
    // Cost: one extra factorization + one residual eval per iteration; not added to
    // the _ts_* solver-cost totals, so it does not pollute the real cost breakdown.
    if( _nTraceVar > 0 && n == m ){
      size_t const tau0 = m - (size_t)_nTraceVar;
      double sc_min=1e300, sc_max=0.0, tc_min=1e300, tc_max=0.0;
      for( size_t j=0; j<m; ++j ){
        double const w = cw[j];
        if( j < tau0 ){ sc_min=std::min(sc_min,w); sc_max=std::max(sc_max,w); }
        else          { tc_min=std::min(tc_min,w); tc_max=std::max(tc_max,w); }
      }
      arma::vec dyp;
      bool ok_nt;
      { OCEnvScopedFdSilence const silence;
        ok_nt = solve_linear( Js, -rs, dyp ) && dyp.is_finite(); }
      double phit_nt = std::numeric_limits<double>::infinity();
      bool   armijo_nt = false, reduces = false;
      if( ok_nt ){
        arma::vec dxp = dyp;
        for( arma::uword j=0; j<dxp.n_elem; ++j ) dxp(j) *= cw[(size_t)j];
        for( size_t i=0; i<m; ++i ) xt[i] = var[i] + dxp(i);
        if( _residual_core( xt.data(), inp, cst, rtrial.data(), A, ic_override ) ){
          double rt2=0.0; for( size_t i=0; i<n; ++i ){ double const s=rtrial[i]*rw[i]; rt2+=s*s; }
          phit_nt   = 0.5*rt2;
          reduces   = ( phit_nt < phi );
          armijo_nt = ( phit_nt < (1.0 - options.SOLVE.ARMIJO_C)*phi );
        }
      }
      _disp( 3 ) << "  [tier0-probe] it=" << std::setw(2) << it+1
                << " imp=" << imp_type_name( options.INTERFACE.IMPOSITION )
                << " nTrace=" << _nTraceVar
                << " solve_linear=" << ( ok_nt ? "OK  " : "FAIL" )
                << std::scientific << std::setprecision(3)
                << " phi=" << phi << " phi_NT=" << phit_nt
                << " NT_reduces=" << ( reduces ? "yes" : "no " )
                << " armijo=" << ( armijo_nt ? "pass" : "fail" )
                << " cw_state=[" << sc_min << "," << sc_max << "]"
                << " cw_tau=[" << tc_min << "," << tc_max << "]\n";
    }
#endif

    bool   accepted = false;
    double best = std::numeric_limits<double>::infinity();
    double best_alpha = 0.0;
    double accepted_phit = std::numeric_limits<double>::infinity();
    int    accepted_tier = -1;   // 0 = direct Newton, 1 = damped LM/QR, 2 = dense min-norm

    // Try step sources in priority order, line-searching each candidate, until one
    // is accepted.  Tier 0 (direct sparse Newton) is primary and forms NO normal
    // equations, so well-conditioned (including genuinely 3D) systems stay fast.
    // Tiers 1/2 are formed lazily, only if the previous tier's step is rejected,
    // so a converging Newton path never pays the JtJ fill-in or dense-solve cost.
    // Once an LM step stalls on a singular system we stick to the dense min-norm
    // tier for the remainder (see stall detection below).
    int const first_tier = ( dense_minnorm ? 2 : ( (nt_dead || spqr_skip_nt) ? 1 : 0 ) );
#if defined(CRONOS__WITH_SPQR)
    // Under SPQR the QR min-norm solve (tier 1) is the terminal rank-tolerant
    // solver; the dense min-norm tier (2) would densify a 1e4^2 system, so cap there.
    int const last_tier  = ( options.SOLVE.FACTORIZATION == Options::SOLVE_SPQR ) ? 1 : 2;
#else
    int const last_tier  = 2;
#endif
    for( int tier = first_tier; tier <= last_tier && !accepted; ++tier ){
      arma::vec dy;
      bool have_step = false;
      char const* label = "";
      auto _tstep = _swnow();

      if( tier == 0 ){
        if( n != m ) continue;        // direct Newton is for square systems only
        // Direct sparse Newton: solve Js dy = -rs.  No JtJ, no fill-in, no kappa^2.
        // Scope-silence the benign reference-BLAS chatter a rank-deficient Js
        // provokes; the step is validated and used only if the line search accepts.
        OCEnvScopedFdSilence const silence;
        have_step = solve_linear( Js, -rs, dy ) && dy.is_finite();
        if( !have_step ) nt_dead = true;  // singular augmented saddle: skip tier-0 hereafter
                                          // (note emitted below, outside the fd silence)
        label = "NT";
      }
      else if( tier == 1 ){
        bool spqr_used = false;
#if defined(CRONOS__WITH_SPQR)
        if( options.SOLVE.FACTORIZATION == Options::SOLVE_SPQR ){
          // Sparse rank-revealing QR min-norm least squares of Js directly:
          // no Js^T Js, so no fill/condition squaring.  lambda is unused here
          // (the QR is the min-norm Gauss-Newton step); the line search globalizes.
          auto _tlsol = _swnow();
          bool const qr_factored = solve_lsq_spqr( Js, -rs, dy );
          have_step = qr_factored && dy.is_finite();
          _ts_lmsolve += _swms( _tlsol, _swnow() );
          label = "QR";
          spqr_used = have_step;   // false -> fall through to the SuperLU path below

          // rev116: SAY SO WHEN THE FALL-THROUGH HAPPENS.
          //
          // On failure this drops into the damped normal equations below, `label` is
          // overwritten with the LM label, and nothing is emitted -- so the log reads as
          // though QR was never selected.  The two situations "SPQR was not chosen for
          // this solve" and "QR was chosen, failed, and fell through" then look identical,
          // and they are different problems.
          //
          // MEASURED, OCFE_MBC5 constant-dP: the failing solve reports tier=1 with
          // lambda=1e-13 -- the LM path -- while every SUBSEQUENT solve on the same
          // 12042-unknown system reports "using QR min-norm steps" and converges in 7
          // iterations to 1.5e-11.  Which of the two situations produced the failure could
          // not be read from the log.  The fall-through matters because the normal
          // equations SQUARE the condition number, and the failing step was 77x the
          // solution vector.
          if( !have_step && options.SOLVE.VERBOSE && !qr_fell_through_noted ){
            qr_fell_through_noted = true;
            _disp( 2 ) << "  note: QR min-norm step UNAVAILABLE at iteration " << it+1
                      << " ("
                      << ( qr_factored ? "SuiteSparseQR returned a step that is not finite"
                                       : "SuiteSparseQR factorization failed" )
                      << "); falling through to damped normal equations (Js^T Js + lambda*D),"
                         " which squares the condition number.  Steps from here are LM, not"
                         " QR, whatever the tier label says.\n";
          }
        }
#endif
        if( !spqr_used ){
          // Damped normal equations (Js^T Js + lambda*D) dy = -Js^T rs, formed
          // lazily.  SPD/nonsingular for lambda>0, so the sparse factorization is
          // quiet even on rank-deficient interface-augmented Jacobians.
          auto _tform = _swnow();
          arma::sp_mat const JtJ = Js.t()*Js;
          arma::vec    const g   = Js.t()*rs;
          arma::umat dloc( 2, (arma::uword)m );
          arma::vec  dval( (arma::uword)m );
          for( arma::uword j=0; j<(arma::uword)m; ++j ){
            dloc(0,j)=j; dloc(1,j)=j;
            dval(j) = lambda*std::max( 1.0, std::fabs( JtJ(j,j) ) );
          }
          arma::sp_mat const D( dloc, dval, (arma::uword)m, (arma::uword)m );
          arma::sp_mat const H = JtJ + D;
          _ts_lmform += _swms( _tform, _swnow() );
          auto _tlsol = _swnow();
          have_step = solve_linear( H, -g, dy ) && dy.is_finite();
          _ts_lmsolve += _swms( _tlsol, _swnow() );
          label = "LM";
        }
      }
      else{
        // Dense minimum-norm (pseudo-inverse) step: arma::solve's least-squares
        // path truncates the ~0 singular values and resolves the consistent part
        // of a genuinely singular Schur system that a damped step cannot.
        arma::mat const JsD( Js );
        arma::vec dyd;
        have_step = arma::solve( dyd, JsD, -rs ) && dyd.is_finite();
        if( have_step ) dy = dyd;
        label = "MN";
      }
      double const _dstep = _swms( _tstep, _swnow() );
      _ts_solve += _dstep;
      if( tier == 0 ) _ts_nt += _dstep; else if( tier == 1 ) _ts_lm += _dstep; else _ts_mn += _dstep;

      // One-time note (chronological: printed when tier-0 is first found singular,
      // before the step that the iteration actually accepts is reported).
      if( options.SOLVE.VERBOSE && tier == 0 && nt_dead && !nt_disabled_noted ){
        nt_disabled_noted = true;
        char const* via = "LM";
#if defined(CRONOS__WITH_SPQR)
        if( options.SOLVE.FACTORIZATION == Options::SOLVE_SPQR ) via = "QR";
#endif
        _disp( 2 ) << "  note: tier-0 Newton singular on the augmented interface system; "
                  << "using " << via << " steps\n";

        run_spectrum_audit( Js, cw, "tier0", &rs );
      }

      if( !have_step ) continue;

      arma::vec dx = dy;                                   // un-scale columns: dx = diag(cw)*dy
      for( arma::uword j=0; j<dx.n_elem; ++j ) dx(j) *= cw[(size_t)j];

      auto _tls = _swnow();
      double alpha = 1.0;
      std::vector<std::pair<double,double>> ls_trace;   // (alpha, merit) for the failure report
      for( int ls=0; ls<options.SOLVE.LS_MAX; ++ls, alpha*=0.5 ){
        for( size_t i=0; i<m; ++i ) xt[i] = var[i] + alpha*dx(i);
        if( !_residual_core( xt.data(), inp, cst, rtrial.data(), A, ic_override ) ) continue;
        double rt2=0.0, tnorm=0.0;
        for( size_t i=0; i<n; ++i ){ double const s = rtrial[i]*rw[i]; rt2 += s*s; tnorm = std::max(tnorm,std::fabs(rtrial[i])); }
        double const phit = 0.5*rt2;
        if( ls_trace.size() < 24 ) ls_trace.emplace_back( alpha, phit );
        best = std::min( best, phit ); best_alpha = alpha;
        if( std::isfinite(phit) && ( phit < (1.0 - options.SOLVE.ARMIJO_C*alpha)*phi || tnorm <= options.SOLVE.RES_TOL ) ){
          std::copy( xt.begin(), xt.end(), var );
          r.swap( rtrial );
          norm = tnorm;
          if( tier != 0 )   // a damped/min-norm step drives lambda; Newton does not
            lambda = std::max( lambda*options.SOLVE.LAMBDA_DROP, options.SOLVE.LAMBDA_MIN );
          accepted = true; ++rep.iterations;
          accepted_phit = phit; accepted_tier = tier;
          if( tier == 1 ) ++rep.lm_steps; else ++rep.newton_steps;
          rep.final_merit = phit;
          if( options.SOLVE.VERBOSE ){
            _disp( 3 ) << label << " it=" << std::setw(2) << it+1;
            if( label[0] == 'L' )   // only the damped LM step is driven by lambda
              _disp( 3 ) << " lambda=" << std::scientific << std::setprecision(2) << lambda;
            _disp( 3 ) << std::scientific << std::setprecision(4)
                      << " alpha=" << alpha
                      << " max|r|=" << norm
                      << " merit=" << phit << "\n";
          }
          break;
        }
      }
      _ts_ls += _swms( _tls, _swnow() );

      // ---- LINE-SEARCH FAILURE REPORT ------------------------------------------------
      //
      // A rejected line search is currently reported only as "step failed at iteration N",
      // which does not distinguish the two causes that need different remedies:
      //
      //   (a) BAD DIRECTION -- the linear solve returned a step dominated by a near-null
      //       direction of an ill-conditioned Jacobian.  Signature: dx concentrated on one
      //       or two states, ||dx||/||var|| enormous, merit rising monotonically as alpha
      //       shrinks and only flattening at alpha ~ 0.  Remedy: a step limiter, a trust
      //       region, or better conditioning -- NOT more interface work.
      //
      //   (b) GENUINE NONLINEARITY -- the direction is fine but the merit is not decreasing
      //       along it.  Signature: dx spread across states, ||dx||/||var|| O(1), merit
      //       decreasing then turning.  Remedy: continuation (tau ramp), a better seed.
      //
      // Reported once per solve, at the first failing iteration only, so a stalling solve
      // does not flood the log.  SOLVE.VERBOSE only; zero cost otherwise.
      if( !accepted && options.SOLVE.VERBOSE && !step_fail_noted ){
        step_fail_noted = true;
        double dxn = 0., varn = 0., dxmax = 0.;
        for( size_t i = 0; i < m; ++i ){
          dxn  += dx(i)*dx(i);
          varn += var[i]*var[i];
          dxmax = std::max( dxmax, std::fabs( dx(i) ) );
        }
        dxn = std::sqrt( dxn ); varn = std::sqrt( varn );

        _disp( 2 ) << "  [step-fail] it=" << it+1 << " tier=" << tier
                  << " |dx|=" << std::scientific << std::setprecision(3) << dxn
                  << " |var|=" << varn
                  << " |dx|/|var|=" << ( varn > 0. ? dxn/varn : INFINITY )
                  << " max|dx|=" << dxmax
                  << "  phi0=" << phi << " best=" << best << "\n";

        // Step composition by state: the discriminator between (a) and (b).
        std::map<std::string,double> en;
        double e_other = 0.;
        for( size_t i = 0; i < m; ++i ){
          double const e = dx(i)*dx(i);
          if( e <= 0. ) continue;
          std::string const d = decode_state_column( i );
          size_t const br = d.find( '[' );
          if( br == std::string::npos ){ e_other += e; continue; }
          en[ d.substr( 0, br ) ] += e;
        }
        double const tot = ( dxn > 0. ? dxn*dxn : 1. );
        std::vector<std::pair<double,std::string>> ord;
        for( auto const& kv : en ) ord.emplace_back( kv.second, kv.first );
        std::sort( ord.begin(), ord.end(), []( auto const& a, auto const& b ){ return a.first > b.first; } );
        _disp( 2 ) << "  [step-fail] step energy by state:";
        for( auto const& p : ord )
          if( p.first/tot > 1.0e-4 )
            _disp( 2 ) << " " << p.second << "=" << std::fixed << std::setprecision(1)
                      << 100.*p.first/tot << "%";
        if( e_other/tot > 1.0e-4 )
          _disp( 2 ) << " trace/tau=" << std::fixed << std::setprecision(1) << 100.*e_other/tot << "%";
        _disp( 2 ) << std::scientific << "\n";

        // Dominant components, decoded -- says WHERE the excursion is, not just in what.
        {
          std::vector<size_t> idx( m );
          for( size_t i = 0; i < m; ++i ) idx[i] = i;
          std::partial_sort( idx.begin(), idx.begin() + std::min<size_t>( 6, m ), idx.end(),
                             [&]( size_t a, size_t b ){ return std::fabs(dx(a)) > std::fabs(dx(b)); } );
          _disp( 2 ) << "  [step-fail] dominant step components:\n";
          for( size_t k = 0; k < 6 && k < m; ++k )
            _disp( 2 ) << "      c" << std::setw(5) << idx[k]
                      << " dx=" << std::scientific << std::setprecision(3) << dx(idx[k])
                      << "  var=" << var[idx[k]]
                      << "  " << decode_state_column( idx[k] ) << "\n";
        }

        // Merit along the ray.  Monotone-rising all the way down is (a); a turning point is (b).
        _disp( 2 ) << "  [step-fail] merit vs alpha:";
        for( auto const& [a,ph] : ls_trace )
          _disp( 2 ) << "  " << std::scientific << std::setprecision(1) << a << ":" << ph;
        _disp( 2 ) << "\n";

        // rev8: re-run the spectrum audit AT THE FAILING ITERATE.  The compatibility
        // measurement (u_k^T r) is what separates a gauge kernel from an inconsistent one,
        // and the iterate that matters is the one the line search stalled on -- at iterate 0
        // the residual is dominated by the initial guess, so the floor there says little about
        // where the merit settles.  Costs one extra dense SVD, only on a failed step and only
        // under AUDIT.SPECTRUM.  The predicted floor should match phi0/best above.
        if( options.AUDIT.SPECTRUM )
          run_spectrum_audit( Js, cw, "stepfail", &rs );
      }
    }

    // Stall detection: an accepted *LM* step that barely reduces the merit while
    // still above tolerance is the signature of a genuinely singular system the
    // damped step cannot resolve -> stick to the dense minimum-norm tier.
    // Skipped under SPQR: tier-1 QR is already the rank-tolerant min-norm solver,
    // and the dense tier it would switch to is disabled at scale.
    if(
#if defined(CRONOS__WITH_SPQR)
        options.SOLVE.FACTORIZATION != Options::SOLVE_SPQR &&
#endif
        !dense_minnorm && accepted && accepted_tier == 1
        && norm > options.SOLVE.RES_TOL && accepted_phit > 0.5*phi ){
      dense_minnorm = true;
      if( options.SOLVE.VERBOSE )
        _disp( 2 ) << "OCFESLV::solve ** damped step stalled (merit ratio "
                  << std::scientific << std::setprecision(2) << (accepted_phit/phi)
                  << "); switching to dense minimum-norm steps\n";
    }

    if( !accepted ){
      lambda *= options.SOLVE.LAMBDA_RAISE;
      // rev52 (B2): the old message reported only the SCALED merit and the best alpha, which
      // is not enough to tell three different failures apart.  All three print identically
      // without the fields added here.
      //
      //  (a) NOTHING LEFT TO DO.  max|r| already at or below SOLVE.RES_TOL: the step failed
      //      because there is no descent left, not because the solve is stuck.  Would be an
      //      ACCEPTANCE-TEST defect, not a numerical one.  NOTE this is a priori unlikely --
      //      133 corpus drivers converge under both exact modes and would misfire the same
      //      way -- but max|r| is the one quantity that settles it and it was not printed.
      //  (b) NEVER STARTED.  accepted=0: the initial guess is outside the region of
      //      attraction and no step was ever taken.  A genuine basin failure.
      //  (c) STARTED THEN STALLED.  accepted>0: progress was made and then the linearisation
      //      stopped being usable.  On a MOVING MESH this is the mechanism to expect and it
      //      is unique to this class of model: an accepted step moves x, which moves the
      //      COLLOCATION POINTS, which changes the discretisation the next Jacobian is built
      //      on -- so the linearisation is valid over a smaller neighbourhood than on any
      //      fixed-mesh model.  No other corpus model exercises this at all.
      //
      // lambda is printed because (c) and a merely ill-conditioned step are distinguished by
      // whether damping was already large when the line search gave up.
      std::cerr << "OCFESLV::solve ** LM/Newton step failed at iteration " << it+1
                << " best_scaled_merit=" << std::scientific << std::setprecision(4) << best
                << " best_alpha=" << best_alpha << std::endl;
      _disp( 2 ) << "OCFESLV::solve **   max|r|=" << std::scientific << std::setprecision(4) << norm
                << "  SOLVE.RES_TOL=" << options.SOLVE.RES_TOL
                << "  ratio=" << ( options.SOLVE.RES_TOL > 0. ? norm/options.SOLVE.RES_TOL : 0. )
                << "  lambda=" << lambda
                << "  accepted_steps=" << rep.iterations
                << ( norm <= options.SOLVE.RES_TOL
                       ? "   [(a) max|r| ALREADY AT TOLERANCE -- acceptance test, not a basin]"
                   : rep.iterations == 0
                       ? "   [(b) NEVER ACCEPTED A STEP -- initial guess outside the basin]"
                       : "   [(c) PROGRESSED THEN STALLED -- linearisation ceased to be usable;"
                         " on a moving mesh, expect the step-moves-the-nodes mechanism]" )
                << std::endl;
      rep.final_residual = norm; rep.converged = false;
      // rev42: report spreads even on a FAILED solve.  Omitting them would silently exclude
      // exactly the models that most need the readout -- MMPDE10 IC_STRONG among them, which
      // is conv=no at production mesh.  The label carries the convergence state so the number
      // cannot be read as an accuracy result by accident: this arc drew wrong conclusions from
      // non-converged iterates at least twice, and "check conv before reading any accuracy
      // number" is a standing method note.
      if( _knob_AUDIT_DUP_SPREAD() )
        _report_dup_spread_gated( var, "post solve  conv=NO -- STRUCTURAL ONLY,"
                                       " NOT AN ACCURACY NUMBER" );
      if( _knob_AUDIT_BETA() )
        _report_beta_gated( var, "post solve", false );
      return rep;
    }
  }

  rep.final_residual = norm;
  rep.converged = ( norm <= options.SOLVE.RES_TOL );

  // [march-window] mechanism-2 localiser (DISPLAY_LEVEL>=1): after a marching window solve, report
  // max|var| per state.  A near-singular window whose null mode a residual-small solve populates
  // shows up as one state (or its interface DOFs) carrying an O(10..100) blow-up while the others
  // stay O(1) -- this says whether the unconstrained/mis-coupled mode lives in U or in an auxiliary
  // derivative state, which is what the (so-far mis-guessed) override fix must actually close.
  if( _inMarch && _solve_display_level() >= 1 ){
    _disp( 2 ) << "  [march-window] max|var| by state:";
    for( auto const& [st, dom] : _mVar ){
      auto itb = _blkVar.find( st ); if( itb == _blkVar.end() ) continue;
      size_t const bpos = std::get<0>( itb->second );
      size_t const btot = std::get<1>( itb->second );
      size_t nblk = 1;
      for( auto const& d : dom ){ auto itd = _mDom.find( d ); if( itd != _mDom.end() ) nblk *= itd->second.n_elem; }
      size_t const hi = std::min<size_t>( bpos + btot*nblk, (size_t)n );
      double mx = 0.0; size_t imx = bpos;
      for( size_t i = bpos; i < hi; ++i ){ double const a = std::fabs( var[i] ); if( a > mx ){ mx = a; imx = i; } }
      _disp( 2 ) << "  " << st << "=" << std::scientific << std::setprecision(3) << mx << "@" << imx;
    }
    _disp( 2 ) << "\n";
  }

  if( options.SOLVE.VERBOSE ){
#if defined(CRONOS__WITH_SPQR)
    bool const qr = ( options.SOLVE.FACTORIZATION == Options::SOLVE_SPQR );
#else
    bool const qr = false;
#endif
    _disp( 2 ) << "  [solve-cost] deriv=" << _ts_deriv << "ms equil=" << _ts_equil
              << "ms assembly=" << _ts_asm << "ms step_solve=" << _ts_solve
              << " (NT=" << _ts_nt << ( qr ? " QR=" : " LM=" ) << _ts_lm
              << " MN=" << _ts_mn << ")";
    if( !qr )   // the JtJ form-vs-factor split is meaningful only for the normal-equations LM
      _disp( 2 ) << " [LM: form=" << _ts_lmform << " solve=" << _ts_lmsolve << "]";
    _disp( 2 ) << "ms line_search=" << _ts_ls << "ms (over " << rep.iterations
              << " accepted it)\n";
  }
  // Unified value-function cache: on a genuine (non-marching) monolithic solve, evaluate the
  // output functions so val_functions() works for monolithic just as _solve_marching populates
  // it for marching.  Skipped inside the march driver (_inMarch=true): the marching loop
  // accumulates the per-window contributions itself.
  if( rep.converged && !_inMarch ){
    _eval_functions_captured( var, inp, cst );
    // Cache the converged primal for reuse by a subsequent solve at the same input.
    if( options.SOLVE.REUSE ){
      _primalVar.assign( var, var + n_colloc_sta() );
      if( inp ) _primalInp.assign( inp, inp + n_colloc_inp() ); else _primalInp.clear();
      _primalValid = true;
    }
  }
  // rev42: automatic per-state duplicate-node spreads, the Stage A gate.  Placed AFTER the
  // value-function cache so a driver's own reporting still sees an untouched var[]; this is
  // read-only and cannot perturb the solve.  Skipped inside the march driver, where the
  // per-window solves would emit one block each and the meaningful frame is the assembled
  // trajectory the marching loop reports.
  if( _knob_AUDIT_DUP_SPREAD() && !_inMarch )
    _report_dup_spread_gated( var, rep.converged ? "post solve  conv=yes"
                                                 : "post solve  conv=NO -- STRUCTURAL ONLY,"
                                                   " NOT AN ACCURACY NUMBER" );
  if( _knob_AUDIT_BETA() && !_inMarch )
    _report_beta_gated( var, "post solve", rep.converged );
  return rep;
}



template <typename CST>
inline bool
OCFESLV::eval
( FFDep* eqn, FFDep* fct, FFDep const* var, FFDep const* inp, CST const* cst )
{
  if( !_issetup ) return false;

  try{
    // Initialise constants
    if( !_vCst.empty() && !cst )
      throw Exceptions( Exceptions::CSTVAL );
    for( size_t i = 0; i < _vCst.size(); ++i )
      _OCDepCst[i].set( cst[i] );

    // Evaluate equation residuals and output functions through the exact
    // OCVar<FFDep> path.  This preserves the minimal FFDep pattern used by
    // external DAG operations such as FFOCFERES.
    if( eqn )
      _eval_eqn( eqn, var, inp, _OCDepEqn, _OCDepVar, _OCDepInp, _OCDepDom, _OCDepCst, _wrkOCDep );

    if( fct )
      _eval_fct( fct, var, inp, _OCDepEqn, _OCDepVar, _OCDepInp, _OCDepDom, _OCDepCst, _wrkOCDep );

    // If this FFDep evaluation used the canonical stacked independent-vector
    // seeding, reuse the exact rows to initialise the derivative cache.  This
    // avoids paying a second OCVar<FFDep> propagation when a test or solver asks
    // for eval<FFDep>() immediately followed by deriv() sparsity/values.
    auto canonical = [&]() -> bool {
      if( _nCollVar && !var ) return false;
      if( _nCollInp && !inp ) return false;
      for( size_t i = 0; i < _nCollVar; ++i ){
        auto const& d = var[i].dep();
        if( d.size() != 1 ) return false;
        auto it = d.begin();
        if( it->first != static_cast<int>( i ) || it->second != FFDep::L ) return false;
      }
      for( size_t i = 0; i < _nCollInp; ++i ){
        auto const& d = inp[i].dep();
        if( d.size() != 1 ) return false;
        auto it = d.begin();
        if( it->first != static_cast<int>( _nCollVar + i ) || it->second != FFDep::L ) return false;
      }
      return true;
    };

    if( canonical()
     && ( !_jacReady
       || _jacInterfaceType  != options.INTERFACE.TYPE
       || _jacImpositionType != options.INTERFACE.IMPOSITION
       || _jacClassified     != _classified ) ){
      size_t const ncol = _nCollVar + _nCollInp;
      _jacPatternEqn.assign( _nCollEqn, std::vector<size_t>() );
      _jacPatternFct.assign( _nCollFct, std::vector<size_t>() );
      _jacLinearEqn.assign( _nCollEqn, true );
      _jacLinearFct.assign( _nCollFct, true );

      auto fill = [&]( FFDep const* dep_rows, size_t nrow,
                       std::vector<std::vector<size_t>>& pat,
                       std::vector<bool>& lin ){
        if( !dep_rows ) return false;
        for( size_t k = 0; k < nrow; ++k ){
          auto const& dep = dep_rows[k].dep();
          pat[k].reserve( dep.size() );
          bool islin = true;
          for( auto const& kv : dep ){
            if( kv.first >= 0 && static_cast<size_t>( kv.first ) < ncol )
              pat[k].push_back( static_cast<size_t>( kv.first ) );
            if( kv.second != FFDep::L ) islin = false;
          }
          std::sort( pat[k].begin(), pat[k].end() );
          pat[k].erase( std::unique( pat[k].begin(), pat[k].end() ), pat[k].end() );
          lin[k] = islin;
        }
        return true;
      };

      bool const have_eqn = !_nCollEqn || fill( eqn, _nCollEqn, _jacPatternEqn, _jacLinearEqn );
      bool const have_fct = !_nCollFct || fill( fct, _nCollFct, _jacPatternFct, _jacLinearFct );
      if( have_eqn && have_fct ){
        _jacAllEqnLinear = true;
        for( bool b : _jacLinearEqn ) if( !b ){ _jacAllEqnLinear = false; break; }
        _jacAllFctLinear = true;
        for( bool b : _jacLinearFct ) if( !b ){ _jacAllFctLinear = false; break; }
        std::vector< std::vector<size_t> > all_pattern( _jacPatternEqn );
        all_pattern.insert( all_pattern.end(), _jacPatternFct.begin(), _jacPatternFct.end() );
        _jacColor = _color_jacobian_columns( all_pattern, ncol, _jacNColor );
        _jacInterfaceType  = options.INTERFACE.TYPE;
        _jacImpositionType = options.INTERFACE.IMPOSITION;
        _jacClassified     = _classified;
        _jacReady          = true;
        // Match _ensure_deriv_cache(): if a canonical eval<FFDep>() built the
        // cache, the compressed FAD work arrays must be resized to the new
        // color count before deriv() reuses this cache.
        _reset_colloc( _OCFVar, _OCFInp, _OCFDom, _OCFCst );
        _OCFEqn = OCVar<t_FAD>();
        _wrkOCF.clear();
        _reset_autodiff( _FVar, _FInp, _FEqn, _FFct );
        _linEqnReady = _linFctReady = false;
        _linEqn.clear();
        _linFct.clear();
        _linCstVal.clear();
      }
    }

    return true;
  }

  catch(...){
    return false;
  }
}

inline bool
OCFESLV::_check_input_continuity
( double const* var, double const* inp, double const* cst, double tol ) const
{
  if( _mInpCont.empty() ) return true;

  for( auto const& [pinp, dirs] : _mInpCont ){
    auto itinp = _mInp.find( pinp );
    if( itinp == _mInp.end() ) continue;
    std::set<FFVar,lt_FFVar> const& idoms = itinp->second;

    for( auto const& [pdir, level] : dirs ){
      if( level < static_cast<int>( CONTINUOUS_C0 ) ) continue;

      FFDom const& dd = _input_dom( pinp, pdir );
      if( dd.n_elem < 2 ) continue;                         // no interior interface
      double const span = dd.up_dom - dd.lo_dom;
      double const eps  = ( span > 0.? span : 1. ) * 1e-5 / double( dd.n_elem );

      // Transverse sample coordinates: element midpoints of every other input domain.
      std::vector<FFVar>               tdir;
      std::vector<std::vector<double>> tsamp;
      for( auto const& dv : idoms ){
        if( dv.id() == pdir.id() ) continue;
        FFDom const& od = _input_dom( pinp, dv );
        std::vector<double> s;
        for( size_t e=0; e<od.n_elem; ++e )
          s.push_back( 0.5*( od.elem_bnd[e] + od.elem_bnd[e+1] ) );
        if( s.empty() ) s.push_back( 0.5*( od.lo_dom + od.up_dom ) );
        tdir.push_back( dv );
        tsamp.push_back( std::move( s ) );
      }

      std::vector<size_t> idx( tdir.size(), 0 );
      while( true ){
        for( size_t e=1; e<dd.n_elem; ++e ){
          double const c = dd.elem_bnd[e];
          t_Coord pt;
          for( size_t k=0; k<tdir.size(); ++k )
            pt[ tdir[k] ] = tsamp[k][ idx[k] ];

          // Sample each side just inside its element, then use the ANALYTIC first
          // derivative (eval_colloc_deriv) to Taylor-extrapolate the value to the exact
          // face coordinate c.  This removes the O(eps) slope term, so a continuous but
          // sloped input is not mistaken for a jump.
          pt[pdir] = c - eps;
          double const vL = eval_colloc<double>      ( pinp, pt, var, inp, cst );
          double const dL = eval_colloc_deriv<double>( pinp, pt, pdir, 1, var, inp, cst );
          pt[pdir] = c + eps;
          double const vR = eval_colloc<double>      ( pinp, pt, var, inp, cst );
          double const dR = eval_colloc_deriv<double>( pinp, pt, pdir, 1, var, inp, cst );

          double const vLf = vL + eps*dL;   // left interpolant extrapolated to the face
          double const vRf = vR - eps*dR;   // right interpolant extrapolated to the face
          double const vscale = std::max( 1., std::max( std::fabs(vLf), std::fabs(vRf) ) );
          if( std::fabs( vLf - vRf ) > tol*vscale ){
            _disp( 3 ) << "OCFESLV: input continuity violated (declared C" << (level-1)
                      << ") across an element interface at coord=" << c
                      << " (value jump~" << std::fabs(vLf - vRf) << ")\n";
            return false;
          }

          if( level >= static_cast<int>( SMOOTH_C1 ) ){
            // Analytic one-sided slopes at the face: a C1 input matches them to O(eps),
            // a kink leaves an O(1) jump.  Relative tolerance separates the two cleanly.
            double const ss = std::max( 1., std::max( std::fabs(dL), std::fabs(dR) ) );
            if( std::fabs( dL - dR ) > 1.e-3*ss ){
              _disp( 3 ) << "OCFESLV: input smoothness violated (declared C1) across an "
                        << "element interface at coord=" << c
                        << " (slope jump~" << std::fabs(dL - dR) << ")\n";
              return false;
            }
          }
        }
        // advance the transverse odometer (no-op when there are no transverse directions)
        if( tdir.empty() ) break;
        size_t k = 0;
        for( ; k<idx.size(); ++k ){ if( ++idx[k] < tsamp[k].size() ) break; idx[k] = 0; }
        if( k == idx.size() ) break;
      }
    }
  }
  return true;
}

inline bool

//==========================================================================================================
// EVALUATION AND DERIVATIVES -- collocated residual and Jacobian evaluation over the working DAG.
// (banner placed from the definitions themselves: this region is lines 25692-26135 of the pre-rev268 file)
//==========================================================================================================

OCFESLV::eval
( double* eqn, double* fct, double const* var, double const* inp, double const* cst )
{
  FixedInputScope fixed_scope_( _fixedInputDepth );  std::vector<double> fixed_buf_;  inp = _with_fixed_inputs( inp, fixed_buf_ );   // fix_input: model-supplied values, outermost call only
  if( !_issetup ) return false;

  try{ 
    // Assert any declared input interface-continuity against the supplied data.  Inert
    // unless the model used the InputContinuity overload of add_input (_mInpCont empty).
    if( inp && !_check_input_continuity( var, inp, cst ) )
      return false;

    // Initialise constants
    if( !_vCst.empty() && !cst )
      throw Exceptions( Exceptions::CSTVAL );

    // Direct eval with reduce_order-generated inputs (captures) and a null buffer: self-allocate zeros
    // -- equations never reference captures, so zeros give correct residuals.  If OUTPUT functions are
    // requested, fill the captured reductions so f(cap) reads completed integrals: MONOLITHIC via a
    // single full-horizon self-fill keyed by var[]; MARCHING via the cap cache from the last converged
    // marching solve (a last-window var[] cannot reconstruct the full-horizon integral), served only
    // when var[] matches the solved point (_marchCapVar) -- else fall through unfilled.
    double const*       useinp = inp;
    std::vector<double> capbuf;
    if( _nCollInp && !inp ){ capbuf.assign( _nCollInp, 0. ); useinp = capbuf.data(); }
    if( fct && _nCollInp && !_deferredValue.empty() && !_inMarch ){
      if( !is_marching() ){                                // monolithic: self-fill from var[]
        if( capbuf.empty() ){ capbuf.assign( _nCollInp, 0. ); if( inp ) std::copy( inp, inp + _nCollInp, capbuf.begin() ); }
        _fill_captures( var, capbuf.data(), cst, -1, std::vector<int>() );
        useinp = capbuf.data();
      }
      else if( _marchCapValid && var && _marchCapVar.size() == _nCollVar
            && std::equal( var, var + _nCollVar, _marchCapVar.begin() ) ){   // marching: serve cap cache
        // CAVEAT: this marching cap-serving lives in eval() only; deriv() stays monolithic-only (marching
        // sensitivity flows through _march_fsens_setup, not direct deriv()).  FFOCFESLV reads outputs via
        // val_functions() and FFOCFERES is the MONOLITHIC full-collocation map (is_marching()==false),
        // so neither wrapper triggers this branch today.  If a SOLVE.MARCHING OCFESLV is ever wrapped in
        // FFOCFERES, eval() would serve caps at the solved point while deriv() would not -- an
        // eval/deriv inconsistency; deriv() would then need the same cache-serving for consistency.
        if( capbuf.empty() ){ capbuf.assign( _nCollInp, 0. ); if( inp ) std::copy( inp, inp + _nCollInp, capbuf.begin() ); }
        for( auto const& C : _deferredValue ){
          if( !C.input.dag() ) continue;
          std::vector<double> const cv = get_input_values( C.input, _marchCapInp.data() );
          set_input_values( C.input, cv, capbuf.data() );
        }
        useinp = capbuf.data();
      }
    }

    // Fast path for linear residual/output maps with setup constants.
    if( _eval_linear_cache( eqn, fct, var, useinp, cst ) )
      return true;

    for( size_t i = 0; i < _vCst.size(); ++i )
      _OCCst[i].set( cst[i] );

    // Evaluate equation residuals
    if( eqn )
      _eval_eqn( eqn, var, useinp, _OCEqn, _OCVar, _OCInp, _OCDom, _OCCst, _wrkOC );

    // Evaluate output functions
    if( fct )
      _eval_fct( fct, var, useinp, _OCEqn, _OCVar, _OCInp, _OCDom, _OCCst, _wrkOC );

    return true;
  }

  catch(...){
    return false;
  }
}




template <typename U>
inline void
OCFESLV::_deriv_fct
( U* fct, U** gradfct, size_t* nnzfct, size_t** colnzfct,
  OCVar<FADType<U>>& OCFfct, std::vector<OCVar<FADType<U>>>& OCFvar,
  std::vector<OCVar<FADType<U>>>& OCFinp, std::vector<OCVar<FADType<U>>>& OCFdom,
  std::vector<OCVar<FADType<U>>>& OCFcst, std::vector<OCVar<FADType<U>>>& OCFwrk,
  FADType<U>* Ffct, FADType<U>* Fvar, FADType<U>* Finp )
const
{
  if( !_jacReady
   || _jacInterfaceType  != options.INTERFACE.TYPE
   || _jacImpositionType != options.INTERFACE.IMPOSITION
   || _jacClassified     != _classified )
    throw Exceptions( Exceptions::SETUP );

  size_t const ncol = _nCollVar + _nCollInp;
  if( gradfct && _jacColor.size() != ncol )
    throw Exceptions( Exceptions::INTERNAL );

  // Populate the sparsity pattern. FFDep pass and coloring were completed by setup().
  _deriv_copy_pattern( _jacPatternFct, nnzfct, colnzfct );

  if( !fct && !gradfct )
    return;

  // This is the only output propagation performed in _deriv_fct().
  _eval_fct( Ffct, Fvar, Finp, OCFfct, OCFvar, OCFinp, OCFdom, OCFcst, OCFwrk );

  // Populate the value and gradient arrays.
  _deriv_extract_rows( Ffct, _jacPatternFct, _jacColor, ncol, fct, gradfct );
}

#ifdef MC__OCFESLV_JAC_FD_PROBE
inline void
OCFESLV::probe_jac_fd
( FFVar const& state, double* var, double const* inp, double const* cst, double h )
{
  FixedInputScope fixed_scope_( _fixedInputDepth );  std::vector<double> fixed_buf_;  inp = _with_fixed_inputs( inp, fixed_buf_ );   // fix_input: model-supplied values, outermost call only
  size_t const neq = _nCollEqn;
  auto itblk = _blkVar.find( state );
  if( itblk == _blkVar.end() ){ _disp( 3 ) << "[jacfd] state not found\n"; return; }
  size_t const col0  = std::get<0>( itblk->second );
  size_t const ncols = std::get<1>( itblk->second );

  // Analytic sparse equation Jacobian at the given point.
  std::vector<size_t> nnz( neq, 0 );
  if( !deriv( nullptr, nullptr, nnz.data(), nullptr, nullptr, nullptr, nullptr, nullptr, var, inp, cst ) ){
    _disp( 3 ) << "[jacfd] deriv(nnz) failed\n"; return; }
  std::vector<std::vector<double>> gvec( neq );
  std::vector<std::vector<size_t>> cvec( neq );
  std::vector<double*> gp( neq ); std::vector<size_t*> cp( neq );
  for( size_t r = 0; r < neq; ++r ){
    gvec[r].resize( nnz[r] ); cvec[r].resize( nnz[r] );
    gp[r] = gvec[r].data(); cp[r] = cvec[r].data(); }
  std::vector<double> eqnbuf( neq, 0. );
  if( !deriv( eqnbuf.data(), gp.data(), nnz.data(), cp.data(), nullptr, nullptr, nullptr, nullptr, var, inp, cst ) ){
    _disp( 3 ) << "[jacfd] deriv(val) failed\n"; return; }

  // Central-difference each column of `state`; compare to the analytic entry row by row.  Also
  // measure STRUCTURAL completeness: entries where FD is nonzero but the analytic pattern has no
  // slot for column d (an "absent-but-nonzero" miss) -- these are dropped columns, not wrong values.
  std::vector<double> fctbuf( _nCollFct, 0. );
  double* fctp = _nCollFct ? fctbuf.data() : nullptr;
  std::vector<double> rp( neq, 0. ), rm( neq, 0. );
  double const struct_tol = 1e-6;   // |fd| above this but absent from pattern = structural miss
  double worst = 0.; size_t worst_col = 0, worst_row = 0; double worst_fd = 0., worst_an = 0.;
  size_t n_missing = 0; double wm_fd = 0.; size_t wm_col = 0, wm_row = 0;
  _disp( 3 ) << "[jacfd] state=" << state << " columns [" << col0 << "," << col0+ncols << ") h=" << h << "\n";
  for( size_t d = col0; d < col0 + ncols; ++d ){
    double const save = var[d];
    var[d] = save + h; if( !eval( rp.data(), fctp, var, inp, cst ) ){ var[d] = save; continue; }
    var[d] = save - h; if( !eval( rm.data(), fctp, var, inp, cst ) ){ var[d] = save; continue; }
    var[d] = save;
    double colworst = 0.; size_t cwrow = 0; double cwfd = 0., cwan = 0.;
    for( size_t r = 0; r < neq; ++r ){
      double const fd = ( rp[r] - rm[r] ) / ( 2.0 * h );
      double an = 0.; bool in_pattern = false;
      for( size_t i = 0; i < nnz[r]; ++i ) if( cvec[r][i] == d ){ an = gvec[r][i]; in_pattern = true; break; }
      if( !in_pattern && std::abs(fd) > struct_tol ){                 // structural miss: dropped column
        ++n_missing;
        if( std::abs(fd) > std::abs(wm_fd) ){ wm_fd = fd; wm_col = d; wm_row = r; }
      }
      double const diff = std::abs( fd - an );
      if( diff > colworst ){ colworst = diff; cwrow = r; cwfd = fd; cwan = an; }
    }
    if( colworst > 1e-4 )
      _disp( 3 ) << "[jacfd]   col " << d << ": max|fd-an|=" << colworst
                << " at row " << cwrow << " (fd=" << cwfd << " an=" << cwan << ")\n";
    if( colworst > worst ){ worst = colworst; worst_col = d; worst_row = cwrow; worst_fd = cwfd; worst_an = cwan; }
  }
  // Structural coverage: how many (row, column) entries the analytic pattern DROPPED despite a
  // nonzero true derivative.  Nonzero here means the analytic Jacobian is structurally incomplete.
  _disp( 3 ) << "[jacfd] STRUCTURAL: " << n_missing << " nonzero-FD entries absent from pattern (|fd|>" << struct_tol << ")";
  if( n_missing ) _disp( 3 ) << " worst fd=" << wm_fd << " at col " << wm_col << " row " << wm_row;
  _disp( 3 ) << "\n";
  _disp( 3 ) << "[jacfd] VALUE:     WORST |fd-an|=" << worst << " at col " << worst_col
            << " row " << worst_row << " (fd=" << worst_fd << " an=" << worst_an << ")\n";
}
#endif // MC__OCFESLV_JAC_FD_PROBE

inline bool
OCFESLV::deriv
( double* eqn, double** gradeqn, size_t* nnzeqn, size_t** colnzeqn,
  double* fct, double** gradfct, size_t* nnzfct, size_t** colnzfct,
  double const* var, double const* inp, double const* cst )
{
  FixedInputScope fixed_scope_( _fixedInputDepth );  std::vector<double> fixed_buf_;  inp = _with_fixed_inputs( inp, fixed_buf_ );   // fix_input: model-supplied values, outermost call only
  if( !_issetup ) return false;

  try{
    if( !_ensure_deriv_cache() )
      throw Exceptions( Exceptions::SETUP );

    size_t const ncol = _nCollVar + _nCollInp;
    if( _jacColor.size() != ncol )
      throw Exceptions( Exceptions::INTERNAL );

    bool const need_eqn = ( eqn || gradeqn );
    bool const need_fct = ( fct || gradfct );
    bool const need_val = ( need_eqn || need_fct );

    // MARCHING + distributed evolution input: colnz here is copied straight from the setup-time pattern
    // (collapsed ELEMENT-0 input layout), but _deriv_eqn/_deriv_fct extract the VALUE at this window's
    // element k=_marchWindow (color _jacColor[nv+k]).  colnz must be remapped to element k to match --
    // otherwise the (now-nonzero) window-k input derivative is labelled column nv+0 and every window
    // collapses onto control DOF 0.  Remap the distributed-evolution-input columns of both patterns; the
    // element-k DOFs are contiguous (offset basek-base0), so this covers n_node>=1 (PWC and higher-order
    // time controls).  Inert off the marching path.
    std::vector<std::vector<size_t>> patEqnRemap, patFctRemap;
    std::vector<std::vector<size_t>> const* pPatEqn = &_jacPatternEqn;
    std::vector<std::vector<size_t>> const* pPatFct = &_jacPatternFct;
    if( _inMarch && _marchWindow > 0 && _evolution_dom_var.dag() && !_mInpDisc.empty() ){
      patEqnRemap = _jacPatternEqn; patFctRemap = _jacPatternFct;
      for( auto const& [uvar, disc] : _mInpDisc ){
        if( !disc.count( _evolution_dom_var ) ) continue;
        std::map<FFVar,size_t,lt_FFVar> e0, ek; for( auto const& kv : disc ){ e0[kv.first]=0; ek[kv.first]=0; }
        ek[ _evolution_dom_var ] = _marchWindow;
        size_t const base0 = pos_input( uvar, e0 ), basek = pos_input( uvar, ek );
        if( basek <= base0 ) continue;
        size_t const shift = basek - base0, span = shift / _marchWindow;
        auto remap = [&]( std::vector<std::vector<size_t>>& P ){
          for( auto& row : P ) for( auto& c : row )
            if( c >= _nCollVar && ( c - _nCollVar ) >= base0 && ( c - _nCollVar ) < base0 + span ) c += shift;
        };
        remap( patEqnRemap ); remap( patFctRemap );
      }
      pPatEqn = &patEqnRemap; pPatFct = &patFctRemap;
    }
    // Pattern-only calls are answered entirely from the setup-time cache.
    if( nnzeqn || colnzeqn )
      _deriv_copy_pattern( *pPatEqn, nnzeqn, colnzeqn );
    if( nnzfct || colnzfct )
      _deriv_copy_pattern( *pPatFct, nnzfct, colnzfct );

    if( !need_val )
      return true;

    // The AD objects and OCVar<t_FAD> workspaces are sized and seeded by
    // setup()/classify_pde().  deriv() only updates primal values via x(),
    // preserving the compressed derivative buffers and avoiding _reset_colloc().
    if( _FVar.size() != _nCollVar || _FInp.size() != _nCollInp
     || _FEqn.size() != _nCollEqn || _FFct.size() != _nCollFct
     || _OCFVar.size() != _mVar.size() || _OCFInp.size() != _mInp.size()
     || _OCFDom.size() != _mDom.size() || _OCFCst.size() != _vCst.size() )
      throw Exceptions( Exceptions::INTERNAL );

    if( _nCollVar && !var )
      throw Exceptions( Exceptions::INTERNAL );
    for( size_t i = 0; i < _nCollVar; ++i )
      _FVar[i].x() = var[i];

    // Direct deriv() with reduce_order-generated inputs (captures) and a null buffer: self-allocate
    // zeros -- equations never reference captures, so zeros are correct for the equation Jacobian.
    // If OUTPUT derivatives are requested, self-fill the captured reductions from the primal so
    // d f/dcap is differentiated at the completed integral (deriv() analogue of eval()'s self-fill
    // and _mono_fsens_setup's pre-deriv fill).  Monolithic only (var[] is a full solution).
    std::vector<double> inpauto;
    if( _nCollInp && !inp ){ inpauto.assign( _nCollInp, 0. ); inp = inpauto.data(); }
    if( need_fct && _nCollInp && !_deferredValue.empty() && !_inMarch && !is_marching() ){
      if( inpauto.empty() ) inpauto.assign( inp, inp + _nCollInp );
      _fill_captures( var, inpauto.data(), cst, -1, std::vector<int>() );
      inp = inpauto.data();
    }

    if( _nCollInp && !inp )
      throw Exceptions( Exceptions::INTERNAL );
    for( size_t i = 0; i < _nCollInp; ++i )
      _FInp[i].x() = inp[i];

    // MARCHING + distributed evolution input: site B loaded element k=_marchWindow into the FADType input
    // slots, but that slot carries color _jacColor[nv+basek] -- which the setup coloring made orthogonal
    // only to element k's MONO neighbours.  Under the collapsed window it co-occurs with state columns
    // sharing that color, so the compressed extract at nv+basek is polluted (|b_input| was 6.87, not 1).
    // Element 0's color _jacColor[nv+base0] IS orthogonal to the collapsed row's states by construction,
    // so re-seed the resident element's DOFs with element-0's colors: the control derivative then lands on
    // the clean element-0 channel, _deriv_extract_rows (original element-0 pattern) returns element k's
    // UNPOLLUTED derivative, and colnz (remapped to element k above) attributes it to control DOF k.
    // Restored after the eval so later windows start clean.  Element 0 itself is never loaded in window
    // k>0, so nothing else contributes to that channel.
    std::vector<std::pair<size_t,unsigned>> _reseedSave;
    if( _inMarch && _marchWindow > 0 && _evolution_dom_var.dag() && !_mInpDisc.empty() && _jacNColor && need_val ){
      for( auto const& [uvar, disc] : _mInpDisc ){
        if( !disc.count( _evolution_dom_var ) ) continue;
        std::map<FFVar,size_t,lt_FFVar> e0, ek; for( auto const& kv : disc ){ e0[kv.first]=0; ek[kv.first]=0; }
        ek[ _evolution_dom_var ] = _marchWindow;
        size_t const base0 = pos_input( uvar, e0 ), basek = pos_input( uvar, ek );
        if( basek <= base0 ) continue;
        size_t const span = ( basek - base0 ) / _marchWindow;      // DOFs per evolution element
        for( size_t j = 0; j < span && basek + j < _nCollInp; ++j ){
          _reseedSave.emplace_back( basek + j, _jacColor[ _nCollVar + basek + j ] );
          _FInp[ basek + j ].diff( _jacColor[ _nCollVar + base0 + j ], _jacNColor );   // element-0's clean color
        }
      }
    }

    if( !_vCst.empty() && !cst )
      throw Exceptions( Exceptions::CSTVAL );
    for( size_t i = 0; i < _vCst.size(); ++i )
      _OCFCst[i].update( &cst[i] );

    if( need_eqn ){
      _deriv_eqn( eqn, gradeqn, nullptr, nullptr,
                  _OCFEqn, _OCFVar, _OCFInp, _OCFDom, _OCFCst, _wrkOCF,
                  _FEqn.empty()? nullptr: _FEqn.data(),
                  _FVar.empty()? nullptr: _FVar.data(),
                  _FInp.empty()? nullptr: _FInp.data() );
    }

    if( need_fct ){
      _deriv_fct( fct, gradfct, nullptr, nullptr,
                  _OCFEqn, _OCFVar, _OCFInp, _OCFDom, _OCFCst, _wrkOCF,
                  _FFct.empty()? nullptr: _FFct.data(),
                  _FVar.empty()? nullptr: _FVar.data(),
                  _FInp.empty()? nullptr: _FInp.data() );
    }

    for( auto const& [dof, col] : _reseedSave ) _FInp[dof].diff( col, _jacNColor );   // restore element k's own color

    return true;
  }
  catch(...){
    return false;
  }
}

inline bool
OCFESLV::deriv
( size_t* nnz, size_t** colnz )
{
  if( !_issetup ) return false;

  // All-row pattern wrapper around the split derivative API.  Rows are stacked
  // as [equation residuals, scalar outputs].
  return deriv( nullptr, nullptr, nnz, colnz,
                nullptr, nullptr,
                nnz   ? nnz   + _nCollEqn : nullptr,
                colnz ? colnz + _nCollEqn : nullptr,
                nullptr, nullptr, nullptr );
}

inline bool
OCFESLV::deriv
( double* gradeqn, double* gradfct, double const* var, double const* inp, double const* cst )
{
  FixedInputScope fixed_scope_( _fixedInputDepth );  std::vector<double> fixed_buf_;  inp = _with_fixed_inputs( inp, fixed_buf_ );   // fix_input: model-supplied values, outermost call only
  if( !_issetup ) return false;
  if( n_colloc_eqn_nnz() && !gradeqn ) return false;
  if( n_colloc_fct_nnz() && !gradfct ) return false;

  try{
    if( !_ensure_deriv_cache() )
      throw Exceptions( Exceptions::SETUP );

    // Build temporary row pointers into the user-supplied contiguous gradient
    // array.  The layout matches the all-row sparsity query exactly:
    //   gradeqn[ row offset + j ] = d row / d colnz[row][j].
    std::vector<double*> vgradeqn( _nCollEqn, nullptr );
    double* row_gradeqn = gradeqn;
    for( size_t k = 0; k < _nCollEqn; ++k ){
      if( !_jacPatternEqn[k].empty() ){
        vgradeqn[k] = row_gradeqn;
        row_gradeqn += _jacPatternEqn[k].size();
      }
    }

    std::vector<double*> vgradfct( _nCollFct, nullptr );
    double* row_gradfct = gradfct;
    for( size_t k = 0; k < _nCollFct; ++k ){
      if( !_jacPatternFct[k].empty() ){
        vgradfct[k] = row_gradfct;
        row_gradfct += _jacPatternFct[k].size();
      }
    }

    return deriv( nullptr, _nCollEqn ? vgradeqn.data() : nullptr, nullptr, nullptr,
                  nullptr, _nCollFct ? vgradfct.data() : nullptr, nullptr, nullptr,
                  var, inp, cst );
  }
  catch(...){
    return false;
  }
}


inline size_t
OCFESLV::n_evolution_elem
()
const
{
  if( !_evolution_dom_set || !_evolution_dom_var.dag() ) return 0;
  auto it = _mDom.find( _evolution_dom_var );
  if( it == _mDom.end() ) return 0;
  return it->second.n_elem;
}

inline void
OCFESLV::_autodetect_marching_transfers
()
{
  // Single source of truth: the marching-transfer states ARE the differential states (states
  // carrying d/d(evolution-var)), identified by the SAME rule the DAE decomposition uses for its
  // dyn set (_states_with_evolution_derivative).  For each such state the mechanism follows the
  // INITIAL equation that CLOSES it -- canonical `state - ic_input` -> data path (write terminal
  // into ic_input); any other closure (e.g. a general steady-state D(cf0-c)-r) -> override path
  // (null ic_input, continuity swap at k>0).  There is no user selection: the set is structurally
  // determined.
  //   The closure is resolved by _match_initial_closures, NOT by "the first INITIAL equation whose
  // subgraph mentions st".  The latter over-attributes: an IC that evaluates another differential
  // state as an ARGUMENT (a moving mesh's u(xi,0) = f(x(xi,0))) was handed to that argument state
  // too, so two states claimed one row set and the k>0 override became inconsistent.  A
  // differential state with no closure of its OWN is pinned at the evolution-LB by the model's
  // algebraic equations (the index-2 equidistribution case): it needs neither an IC row nor a
  // transfer, and is correctly omitted here.
  _marchTransfers.clear();
  if( !_evolution_dom_var.dag() ) return;

  std::set<FFVar,lt_FFVar> const dyn = _states_with_evolution_derivative( _usr._mEqnUsr, _usr._mVarUsr );
  std::map<FFVar,size_t,lt_FFVar> const closure =
    _match_initial_closures( _usr._mEqnUsr, _usr._mVarUsr, dyn, "autodetect_marching_transfers" );

  size_t ncanon = 0, novr = 0, nomit = 0;
  for( auto const& st : dyn ){
    auto const itc = closure.find( st );
    if( itc == closure.end() ){                        // no INITIAL closure of its own
      ++nomit;
      if( _solve_display_level() >= 1 )
        _disp( 2 ) << "OCFESLV::setup ** marching transfers: differential state " << st
                  << " is pinned by the model's own rows, not by initial data -- an all-node algebraic relation, or a hidden constraint level materialised at the evolution boundary (REDUCE.HIDDEN_IC), which is an INVARIANT and is re-imposed at every window start.  No transfer is needed and none is created"
                     " (evolution-LB value fixed by the algebraic closure).\n";
      continue;
    }
    // Classify the CLOSING equation.  Matching against the USER model (stable ids, preserved into
    // the local DAG) is what terminal_profile/transfer_terminal expect.
    FFVar expr = _usr._mEqnUsr[ itc->second ].var, ic_input;
    FFGraph* dag = static_cast<FFGraph*>( expr.dag() );
    int nstate = 0, ninput = 0;
    FFVar st_var;
    if( dag ){
      FFSubgraph sg = dag->subgraph( 1, &expr );
      for( auto const* op : sg.l_op ){
        if( !op || op->type != FFOp::VAR ) continue;
        FFVar* pv = op->varout[0];
        if( !pv ) continue;
        if     ( _usr._mVarUsr.find( *pv ) != _usr._mVarUsr.end() ){ st_var = *pv; ++nstate; }
        else if( _usr._mInpUsr.find( *pv ) != _usr._mInpUsr.end() ){ ic_input = *pv; ++ninput; }
      }
    }
    // The DATA path writes each window's terminal value into ic_input's slot, so it is only valid when
    //  (1) the closure IS state - ic_input (up to sign): unit coefficients of opposite sign, zero remainder --
    //      counting one state and one input is not enough: x1 - (1.1 + 0.01 p) has one of each, and taking it
    //      for canonical overwrote the PARAMETER p with x1's terminal at every window (OCFE_EVALSUM, 2026-09-28);
    //  (2) ic_input occurs in NO other equation and no output -- otherwise overwriting its slot corrupts them.
    // Anything else takes the general-IC override path, which is built for arbitrary closures.
    bool canonical = false;
    if( nstate == 1 && ninput == 1 && ic_input.dag() && st_var.dag() ){
      std::vector<FFVar> const vv{ st_var, ic_input };
      std::vector<FFVar> A, B;
      try{
        A = dag->FAD( std::vector<FFVar>{ expr }, vv );
        B = dag->compose( std::vector<FFVar>{ expr }, vv, std::vector<FFVar>{ FFVar( 0. ), FFVar( 0. ) } );
      }
      catch( ... ){ A.clear(); B.clear(); }
      bool const unit = A.size() == 2 && A[0].cst() && A[1].cst() && std::fabs( A[0].num().val() ) == 1.
                        && A[0].num().val() == -A[1].num().val();
      bool const zero = B.size() == 1 && B[0].cst() && B[0].num().val() == 0.;
      auto uses = [&]( FFVar const& e ){
        if( !e.dag() ) return false;
        FFSubgraph g = static_cast<FFGraph*>( e.dag() )->subgraph( 1, &e );
        for( auto const* op : g.l_op )
          if( op && op->type == FFOp::VAR && op->varout[0] && op->varout[0]->id() == ic_input.id() ) return true;
        return false; };
      bool alone = true;
      for( size_t e = 0; alone && e < _usr._mEqnUsr.size(); ++e )
        if( e != (size_t)itc->second && uses( _usr._mEqnUsr[e].var ) ) alone = false;
      for( auto const& f : var_output() ) if( alone && uses( f.var ) ) alone = false;
      canonical = unit && zero && alone;
    }
    if( canonical ){
      _marchTransfers.emplace_back( st, ic_input );    // data path
      ++ncanon;
    }
    else {
      _marchTransfers.emplace_back( st, FFVar() );     // general-IC override path
      ++novr;
    }
  }

  if( _solve_display_level() >= 1 )
    _disp( 2 ) << "[march-transfers] differential=" << dyn.size()
              << " canonical=" << ncanon << " override=" << novr
              << " omitted(no own IC)=" << nomit << "\n";
  // rev98: the same line, unconditional under the audit flag.
  //
  // The DISPLAY_LEVEL>=1 emission above is invisible on any driver that runs quiet, which
  // is most of them -- so `grep march-transfers` across the corpus returned "(none)" for
  // OCFE_PDE0, DAE5, MESHONLY10 and CSTR1 and looked like "all canonical" when it only
  // meant "not verbose".  That is the THIRD verbosity-gated diagnostic to mislead in this
  // arc, after `[march] step` (SOLVE.VERBOSE, rev93) which invalidated handoff rev25 s0.3.
  //
  // It matters here because the canonical/override split is the leading candidate for why
  // OCFE_PDE9 leaves max|r|=4.9e-07 on the assembled system after marching while OCFE_PDE0
  // -- which also gates on a post-solve oc.eval() (OCFE_PDE0.cpp:830) -- leaves 2.772e-12.
  // PDE9 is canonical=0 override=1.  If PDE0 is canonical=1 override=0, the residual is a
  // property of OVERRIDE transfers, not of marching, which is both narrower and fixable.
  {
    static bool const kAuditMT = Options::_env_flag( "CRONOS_AUDIT_EVOFACE", false );
    if( kAuditMT )
      _disp( 3 ) << "  [evoface] march-transfers: differential=" << dyn.size()
                << "  canonical=" << ncanon << "  override=" << novr
                << "  omitted(no own IC)=" << nomit << "\n";
  }
}


// ======================================================================
// OCFESLV::_build_march_ic_rows
//   For each override-path transfer state (null ic_input), locate the INITIAL-role equation that
//   closes it and read off its residual rows from the frozen provenance map _rowEqnNdx, together
//   with the state's LB-node column.  Anchors on the EQUATION (INITIAL role + the differential
//   state that equation CLOSES, per _match_initial_closures -- not merely one it references),
//   never on a variable->row search, so the mapping is unambiguous and one row serves one state.
// ======================================================================
inline void
OCFESLV::_build_march_ic_rows
()
{
  _marchICRows.clear();
  _marchICColToNode.clear();
  _marchICRowsReady = false;
  if( !_rowProvenanceReady ) return;   // provenance not built yet (one eval must run first)

  // Resolve each override-path state's OWN INITIAL closure.  The map is INJECTIVE, so no two states
  // can be handed the same equation and therefore the same residual rows -- the defect that made the
  // k>0 override write two value-pins into one row set (dropping one state's continuity and leaving
  // the other unclosed at the window LB).  `claimed` below is the belt-and-braces row-level guard.
  std::set<FFVar,lt_FFVar> ovr_states;
  for( auto const& pr : _marchTransfers )
    if( !pr.second.dag() && pr.first.dag() ) ovr_states.insert( pr.first );
  std::map<FFVar,size_t,lt_FFVar> const closure =
    _match_initial_closures( _mEqn, _mVar, ovr_states, "build_march_ic_rows" );

  std::set<size_t> claimed;            // residual rows already assigned to some override state

  for( auto const& pr : _marchTransfers ){
    if( pr.second.dag() ) continue;    // canonical (has ic_input): keep the data-driven path
    FFVar const& st = pr.first;
    if( !st.dag() ) continue;

    // The INITIAL-role equation that CLOSES st (never one that merely references it).
    auto const itc = closure.find( st );
    if( itc == closure.end() ) continue;   // no INITIAL closure of its own -> nothing to override
    int const eqn_ndx = (int)itc->second;

    // LB-face columns of st in the collapsed marching template, one per spatial node, in
    // terminal_profile() order.  These define the continuity column <-> spatial-node bijection:
    // colToNode[lbcols[k]] = k, so the solve-time override can (a) read each INITIAL row's own
    // continuity column from that row's Jacobian nonzeros -- robust to any row-emission vs
    // var_dom-set ordering -- and (b) index _marchICTerminal (== terminal_profile order) by node.
    // The stored per-row lb_col below is only a structural fallback for the same-order case.
    std::vector<size_t> const lbcols = _march_ic_lb_columns( st );
    { auto& c2n = _marchICColToNode[st];
      for( size_t k = 0; k < lbcols.size(); ++k ) c2n[ lbcols[k] ] = k; }

    std::vector<t_MarchICRow> rows;
    { size_t i = 0;
      for( size_t row = 0; row < _rowEqnNdx.size(); ++row ){
        if( _rowEqnNdx[row] != eqn_ndx ) continue;
        size_t const col = ( i < lbcols.size() ) ? lbcols[i]
                         : ( lbcols.empty() ? 0 : lbcols.back() );     // structural fallback only
        ++i;                                                           // keep lbcols alignment
        if( !claimed.insert( row ).second ){                           // unreachable while closure is injective
          if( _solve_display_level() >= 1 )
            std::cerr << "OCFESLV::solve ** marching IC override: residual row " << row
                      << " already claimed by another transfer state; skipped for " << st
                      << " (the override would have been inconsistent).\n";
          continue;
        }
        rows.push_back({ row, col });
      }
    }
    if( !rows.empty() ) _marchICRows[st] = std::move( rows );
  }
  _marchICRowsReady = true;

  if( _solve_display_level() >= 1 ){
    _disp( 2 ) << "[march-ic-rows] override_states=" << _marchICRows.size();
    for( auto const& [st, rows] : _marchICRows ) _disp( 2 ) << " " << st << ":" << rows.size();
    _disp( 2 ) << " total_rows=" << claimed.size() << "\n";
  }
}


inline bool
OCFESLV::_march_ic_lb_trace
( FFVar const& st, size_t& stride_evo, std::vector<double>& w )
const
{
  stride_evo = 0; w.clear();
  if( !_evolution_dom_set || !_evolution_dom_var.dag() ) return false;
  auto itv = _mVar.find( st );
  if( itv == _mVar.end() ) return false;
  auto const& var_dom = itv->second;
  auto const  evo_id  = _evolution_dom_var.id();
  size_t stride = 1;
  bool   found  = false;
  for( auto const& v : var_dom ){
    auto itd = _mDom.find( v );
    size_t const nn = ( itd != _mDom.end() ) ? itd->second.n_node : 1;
    if( v.id() == evo_id ){
      if( itd == _mDom.end() ) return false;
      // reference-element left endpoint: p(t_lo) = dot( w, y ) over that element's nodes
      if( !itd->second.w_eval_lagrange( w, -1.0, 0 ) ) return false;
      if( w.size() != nn ) { w.clear(); return false; }
      stride_evo = stride;
      found = true;
    }
    stride *= nn;
  }
  return found && !w.empty();
}

inline std::vector<size_t>
OCFESLV::_march_ic_lb_columns
( FFVar const& st )
const
{
  std::vector<size_t> cols;
  if( !_evolution_dom_set || !_evolution_dom_var.dag() ) return cols;
  auto itv = _mVar.find( st );
  if( itv == _mVar.end() ) return cols;
  auto const& var_dom = itv->second;
  auto const  evo_id  = _evolution_dom_var.id();

  // Per-domain nodes/element in var_dom order (for the intra-block node offset, first domain
  // fastest -- see eval_colloc); spatial domains (var_dom minus the evolution var) also feed the
  // spatial-node enumeration, element-major/node-inner, first spatial domain fastest, matching
  // terminal_profile() exactly.
  std::vector<size_t> dom_nn;                        // nodes/element per domain (var_dom order)
  std::vector<FFVar>  sdom;                           // spatial domains (var_dom order)
  std::vector<size_t> snode, snelem;                  // nodes/element and #elements per spatial domain
  for( auto const& v : var_dom ){
    auto itd = _mDom.find( v );
    size_t const nn = ( itd != _mDom.end() ) ? itd->second.n_node : 1;
    size_t const ne = ( itd != _mDom.end() ) ? itd->second.n_elem : 1;
    dom_nn.push_back( nn );
    if( v.id() == evo_id ) continue;                  // evolution is collapsed; its LB node index is 0
    sdom.push_back( v ); snode.push_back( nn ); snelem.push_back( ne );
  }

  size_t total = 1;
  for( size_t s = 0; s < sdom.size(); ++s ) total *= snode[s] * snelem[s];

  cols.reserve( total );
  for( size_t k = 0; k < total; ++k ){
    std::map<FFVar,size_t,lt_FFVar> ndx_el;           // element index per domain -> pos_state
    std::map<FFVar,size_t,lt_FFVar> inode;            // node-within-element per spatial domain
    ndx_el[ _evolution_dom_var ] = 0;                 // LB element of the collapsed evolution direction
    size_t rem = k;
    for( size_t s = 0; s < sdom.size(); ++s ){        // first spatial domain fastest (terminal_profile order)
      size_t const Ns = snode[s] * snelem[s];
      size_t const g  = rem % Ns; rem /= Ns;
      ndx_el[ sdom[s] ] = g / snode[s];               // element index (element-major)
      inode[ sdom[s] ]  = g % snode[s];               // node within element (node-inner)
    }
    size_t const blk = pos_state( st, ndx_el );       // block start (node 0 in every direction)
    size_t off = 0, stride = 1;                       // intra-block offset, first domain fastest, evo node = 0
    { size_t di = 0;                                  // var_dom is a std::set -> iterate by range, index dom_nn in step
      for( auto const& v : var_dom ){
        size_t const nd = ( v.id() == evo_id ) ? 0 : inode.at( v );
        off    += nd * stride;
        stride *= dom_nn[di];
        ++di;
      }
    }
    cols.push_back( blk + off );
  }
  return cols;
}


inline void
OCFESLV::_collapse_evolution_for_marching
()
{
  _marchGrid.clear();
  // Marching is AVAILABLE only when enabled AND a transfer pair exists: an evolution
  // domain is set and at least one (state, IC-input) pair was designated OR auto-detected
  // from the INITIAL equations.  When not available -- e.g. every ordinary monolithic
  // model -- fall through to the monolithic path with no collapse and no failure.
  if( !options.SOLVE.MARCHING || !_evolution_dom_set || !_evolution_dom_var.dag() ) return;
  _autodetect_marching_transfers();      // populate _marchTransfers from INITIAL eqns if unset
  if( _marchTransfers.empty() ) return;  // not available -> monolithic
  auto itd = _mDom.find( _evolution_dom_var );
  if( itd == _mDom.end() ) return;
  FFDom& d = itd->second;
  if( d.elem_bnd.size() < 2 || d.elem_len.empty() ) return;   // cannot collapse -> monolithic

  // section 4 causality guard: a domain-less (scalar/lumped) equation whose residual reduces a
  // distributed field over the EVOLUTION direction couples ALL windows together -- its value
  // depends on the whole trajectory, which is exactly the structure marching assumes is absent.
  // A marched free-accumulate (see _eval_eqn LOOP 2) would sum only the resident window, silently
  // returning a per-window integral instead of the horizon integral.  This coupling is invisible
  // to the classifier (the algebraic scalar row is excluded from the principal symbol), so detect
  // it here and stay monolithic -- where the whole evolution domain is present and the
  // free-accumulate is correct.  (A scalar OUTPUT integrated over time is fine: that is the
  // functional path, window-summed via val_functions(); this guard is specific to state closures.)
  {
    auto const evo_id = _evolution_dom_var.id();
    for( size_t ne = 0; ne < _sgEqn.size() && ne < _mEqn.size(); ++ne ){
      if( !_mEqn[ne].dom.empty() ) continue;   // only domain-less (lumped/scalar) equations
      bool couples_evo = false;
      for( auto itop = _sgEqn[ne].l_op.cbegin();
           itop != _sgEqn[ne].l_op.cend() && !couples_evo; ++itop ){
        auto const& op = *itop;
        if( op->type != FFOp::VAR || _constant( *op->varout[0] ) ) continue;
        std::set<FFVar,lt_FFVar> const* pdom = nullptr;
        auto it_state = _mVar.find( *op->varout[0] );
        if( it_state != _mVar.end() ) pdom = &it_state->second;
        else{
          auto it_input = _mInp.find( *op->varout[0] );
          if( it_input != _mInp.end() ) pdom = &it_input->second;
        }
        if( !pdom ) continue;
        for( auto const& dv : *pdom ) if( dv.id() == evo_id ){ couples_evo = true; break; }
      }
      if( couples_evo ){
        _disp( 2 ) << "OCFESLV::_collapse_evolution_for_marching ** domain-less (scalar/lumped) "
                     "equation #" << ne << " reduces a distributed field over the evolution "
                     "direction '" << _evolution_dom_var << "'; this all-windows coupling is "
                     "non-causal for marching -- falling back to a monolithic solve." << std::endl;
        return;   // stay monolithic: _marchGrid empty -> n_march_steps()==0 -> is_marching()==false
      }
    }
  }

  // Record the full evolution grid (N+1 element boundaries) as the march grid.
  _marchGrid = d.elem_bnd;

  // A distributed input over the EVOLUTION direction that uses the PHYSICAL domain
  // discretisation -- the base add_input(Var,{t}) form, which records NO _mInpDisc entry
  // -- would otherwise resolve _input_dom(Var,t) to the _mDom[t] we are about to collapse
  // to n_elem=1.  That (a) fails the full-grid validation in setup() ( d.n_elem+1 !=
  // _marchGrid.size() -> throw INDEX at the input-block loop ) and (b) mis-locates the input
  // element in the per-window residual eval.  Pin its FULL (pre-collapse) grid into _mInpDisc
  // now, so _input_dom() returns the complete per-window grid exactly as it already does for
  // explicit-discretisation inputs ( add_input(Var,{t},type,n_node), which DAE2/DAE3 used and
  // which therefore are byte-for-byte unaffected -- they already carry an _mInpDisc entry ).
  // Non-evolution input directions are untouched; runs only on the genuine marching-collapse
  // path (this function early-returns off it).
  {
    FFDom const full_evo = d;   // copy the FULL evolution domain BEFORE it is collapsed below
    for( auto const& kv : _mInp ){
      FFVar const& inp = kv.first;
      bool over_evo = false;
      for( auto const& dv : kv.second )
        if( dv.id() == _evolution_dom_var.id() ){ over_evo = true; break; }
      if( !over_evo ) continue;
      auto& disc = _mInpDisc[ inp ];                    // inserts an empty direction map if absent
      if( disc.find( _evolution_dom_var ) == disc.end() )
        disc[ _evolution_dom_var ] = full_evo;          // base-form input: pin the full grid
    }
  }

  // Collapse the evolution domain to its FIRST element [t0,t1] in place.  The reference
  // collocation-node cache (built by add_domain via set_nodes) is element-count
  // independent and is preserved -- only the element partition shrinks.  This runs before
  // any OCVar/node cache is built, so the vector resize is safe.
  double const t0 = d.elem_bnd.front();
  double const t1 = d.elem_bnd[1];
  d.elem_bnd = { t0, t1 };
  d.elem_len = { t1 - t0 };
  d.n_elem   = 1;
  d.w_elem   = t1 - t0;
  d.lo_dom   = t0;
  d.up_dom   = t1;
}

  // rev101 CONSOLIDATION: predicate (ii) / delta (B) probe removed.
  //
  // rev87 computed `_march_would_be_available` alongside the live gate to answer whether a
  // WIDER causal gate would help.  It answered: delta (B) fires only where marching could
  // have run but did not, and on auto-detect models that is an ORDERING ACCIDENT, which
  // Path B fixes at the source (handoff rev28 s0).  With Path B promoted in rev99 the
  // predicate has no remaining non-trivial domain.
  //
  // Nothing measurable is lost.  The hoist probe reports `inferable_here`, and rev98's
  // [evoface] march-transfers line reports differential/canonical/override -- together they
  // say everything gate_proposed said, at the site where it now matters.


inline bool
OCFESLV::evolution_slab_selfcheck
()
const
{
  size_t const nel = n_evolution_elem();
  if( !nel ){
    if( _solve_display_level() >= 1 )
      _disp( 3 ) << "OCFESLV::marching [slabmap] no evolution domain -- nothing to partition\n";
    return true;
  }

  // Total evolution-bearing state DOFs (slab-independent reference count).
  size_t total_evo_dof = 0;
  auto const evo_id = _evolution_dom_var.id();
  for( auto const& [var, var_dom] : _mVar ){
    bool on_evo = false;
    for( auto const& v : var_dom ) if( v.id() == evo_id ){ on_evo = true; break; }
    if( !on_evo ) continue;
    auto itblk = _blkVar.find( var );
    if( itblk == _blkVar.end() ) continue;
    size_t sz = std::get<1>( itblk->second );
    for( auto const& v : var_dom ) sz *= _mDom.at( v ).n_elem;
    total_evo_dof += sz;
  }

  // Union the slabs; require pairwise-disjoint and exact cover.
  std::set<size_t> seen;
  size_t sum = 0;
  bool disjoint = true;
  for( size_t k = 0; k < nel; ++k ){
    auto const c = _evolution_slab_cols( k );
    sum += c.size();
    for( size_t idx : c ) if( !seen.insert( idx ).second ) disjoint = false;
    if( _solve_display_level() >= 1 )
      _disp( 3 ) << "OCFESLV::marching [slabmap] slab " << k << ": " << c.size() << " cols"
                << ( c.empty() ? std::string()
                     : " [" + std::to_string(c.front()) + ".." + std::to_string(c.back()) + "]" )
                << "\n";
  }
  bool const cover = ( sum == total_evo_dof ) && ( seen.size() == total_evo_dof );
  bool const ok    = disjoint && cover;
  if( _solve_display_level() >= 1 || !ok )
    _disp( 3 ) << "OCFESLV::marching [slabmap] n_evo_elem=" << nel
              << " union=" << seen.size() << " sum=" << sum
              << " evo_dof=" << total_evo_dof
              << " disjoint=" << ( disjoint ? "y" : "n" )
              << " cover=" << ( cover ? "y" : "n" )
              << " => " << ( ok ? "OK" : "MISMATCH" ) << "\n";
  return ok;
}

inline std::vector<int>
OCFESLV::evolution_col_to_elem
() const
{
  std::vector<int> c2e( n_colloc_sta(), -1 );
  size_t const nel = n_evolution_elem();
  for( size_t k = 0; k < nel; ++k )
    for( size_t c : _evolution_slab_cols( k ) )
      if( c < c2e.size() ) c2e[c] = static_cast<int>( k );
  return c2e;
}

inline bool
OCFESLV::evolution_slab_map_selfcheck
() const
{
  if( !_evolution_dom_set || !_evolution_dom_var.dag() ) return false;
  if( !const_cast<OCFESLV*>( this )->_ensure_deriv_cache() ) return false;
  auto const   c2e  = evolution_col_to_elem();
  size_t const nel  = n_evolution_elem();
  size_t const nrow = _jacPatternEqn.size();
  if( !nel ) return false;

  std::vector<size_t> rows_per( nel, 0 ), cols_per( nel, 0 );
  size_t cols_unclassified = 0, rows_unclassified = 0;
  bool   bidiagonal = true;

  for( int e : c2e ){
    if( e < 0 )                ++cols_unclassified;
    else if( (size_t)e < nel ) ++cols_per[ (size_t)e ];
  }

  for( size_t r = 0; r < nrow; ++r ){
    int row_elem = -1;
    for( size_t c : _jacPatternEqn[r] )
      if( c < c2e.size() && c2e[c] > row_elem ) row_elem = c2e[c];
    if( row_elem < 0 ){ ++rows_unclassified; continue; }
    ++rows_per[ (size_t)row_elem ];
    // block-lower-bidiagonal: every classified column lies in {row_elem-1, row_elem}
    for( size_t c : _jacPatternEqn[r] ){
      if( c >= c2e.size() ) continue;
      int const ce = c2e[c];
      if( ce >= 0 && ce < row_elem - 1 ) bidiagonal = false;
    }
  }

  bool square = true;
  for( size_t k = 0; k < nel; ++k ) if( rows_per[k] != cols_per[k] ) square = false;

  if( _solve_display_level() >= 1 ){
    _disp( 3 ) << "OCFESLV::marching [slabmap-rows] n_evo_elem=" << nel
              << " n_row=" << nrow << " n_col=" << n_colloc_sta() << "\n";
    for( size_t k = 0; k < nel; ++k )
      _disp( 3 ) << "  elem " << k << ": rows=" << rows_per[k] << " cols=" << cols_per[k]
                << ( rows_per[k] == cols_per[k] ? " square" : " NON-SQUARE" ) << "\n";
    _disp( 3 ) << "  unclassified rows=" << rows_unclassified
              << " cols=" << cols_unclassified
              << "  square=" << ( square ? "y" : "n" )
              << "  bidiagonal=" << ( bidiagonal ? "y" : "n" ) << "\n";
  }
  return square && bidiagonal;
}

inline double
OCFESLV::terminal_value
( FFVar const& state, t_Coord spatial_pt,
  double const* var, double const* inp, double const* cst )
const
{
  FixedInputScope fixed_scope_( _fixedInputDepth );  std::vector<double> fixed_buf_;  inp = _with_fixed_inputs( inp, fixed_buf_ );   // fix_input: model-supplied values, outermost call only
  // Pin the evolution coordinate to the domain upper bound, then interpolate the
  // solved coefficients there (eval_colloc locates the containing element and
  // evaluates the tensor-product Lagrange interpolant -- so this is exact at the
  // UB regardless of whether the collocation family carries a UB node).
  if( _evolution_dom_set && _evolution_dom_var.dag() ){
    auto it = _mDom.find( _evolution_dom_var );
    if( it != _mDom.end() && !it->second.elem_bnd.empty() )
      spatial_pt[ _evolution_dom_var ] = it->second.elem_bnd.back();
  }
  return eval_colloc<double>( state, spatial_pt, var, inp, cst );
}

inline std::map<FFVar,std::vector<double>,lt_FFVar>
OCFESLV::terminal_profile
( double const* var, double const* inp, double const* cst )
const
{
  std::map<FFVar,std::vector<double>,lt_FFVar> out;
  if( !_evolution_dom_set || !_evolution_dom_var.dag() ) return out;
  auto const evo_id = _evolution_dom_var.id();

  for( auto const& [state, var_dom] : _mVar ){
    bool on_evo = false;
    for( auto const& v : var_dom ) if( v.id() == evo_id ){ on_evo = true; break; }
    if( !on_evo ) continue;

    // Each spatial NODE of each spatial ELEMENT reads ITS OWN element's terminal (2026-10-06, WORKPLAN 1.2/1.3): the
    // Lagrange extrapolation to the window's UB of that element's DOFs, exactly the map terminal_profile_adjoint()
    // transposes.  Formerly each node was evaluated at its COORDINATE, so at a spatial element boundary both copies
    // read one element's value: harmless under IC_STRONG/IC_TRACE (the copies are equal), but under IC_WEAK the
    // other element's own terminal was ignored -- marching then differed from monolithic's causal continuity
    // (each copy pinned to its own element) and its adjoint, the transpose of a per-element map, was wrong.
    std::vector<double> w_evo;
    auto ievo = _mDom.find( _evolution_dom_var );
    if( ievo == _mDom.end() || !ievo->second.w_eval_lagrange( w_evo, 1.0, 0 ) ) continue;
    size_t const D = var_dom.size();
    std::vector<size_t> nnode( D, 0 ), nStride( D, 0 );
    size_t evoDim = D;
    { size_t st = 1, d = 0;
      for( auto const& v : var_dom ){
        size_t const nn = _mDom.at( v ).n_node;  nnode[d] = nn;  nStride[d] = st;  st *= nn;
        if( v.id() == evo_id ) evoDim = d;
        ++d; } }
    if( evoDim == D ) continue;
    std::vector<FFVar> sVar; std::vector<size_t> sPos, sElem, sNode;
    { size_t d = 0;
      for( auto const& v : var_dom ){
        if( v.id() != evo_id ){ FFDom const& dd = _mDom.at( v ); sVar.push_back( v ); sPos.push_back( d ); sElem.push_back( dd.n_elem ); sNode.push_back( dd.n_node ); }
        ++d; } }
    size_t total = 1;
    for( size_t k = 0; k < sVar.size(); ++k ) total *= sElem[k] * sNode[k];
    std::vector<double> vals;  vals.reserve( total );
    size_t const nme = std::min( w_evo.size(), nnode[evoDim] );
    for( size_t lin = 0; lin < total; ++lin ){
      std::map<FFVar,size_t,lt_FFVar> ndx_el;
      // the LAST evolution element: the resident window when marching (collapsed: element 0), the final element on
      // the full domain (a monolithic caller, e.g. OCFE_marching's read of the terminal)
      ndx_el[ _evolution_dom_var ] = ievo->second.n_elem? ievo->second.n_elem - 1: 0;
      size_t rem = lin, node_base = 0;
      for( size_t k = 0; k < sVar.size(); ++k ){
        size_t const ng = sElem[k] * sNode[k];
        size_t const g  = rem % ng;  rem /= ng;
        ndx_el[ sVar[k] ] = g / sNode[k];
        node_base += ( g % sNode[k] ) * nStride[ sPos[k] ];
      }
      size_t const inc = pos_state( state, ndx_el );
      double v = 0.;
      for( size_t m = 0; m < nme; ++m ) v += w_evo[m] * var[ inc + node_base + m * nStride[evoDim] ];
      vals.push_back( v );
    }
    out[ state ] = std::move( vals );
  }
  return out;
}


inline bool
OCFESLV::transfer_terminal
( FFVar const& state, FFVar const& input,
  double const* var, double* inp, double const* cst )
const
{
  // Exact marching transfer: seed `input` with `state` at the evolution UB, sampled
  // at the input's own nodes.  terminal_value() pins the evolution coordinate to the
  // UB, so the sampled `pt` need only carry spatial coordinates.
  _sampleOwnElement = true;     // each element-boundary copy reads its own element (see sample_input)
  bool const ok = sample_input( input,
    [&]( t_Coord const& pt ){ return terminal_value( state, pt, var, inp, cst ); },
    inp );
  _sampleOwnElement = false;
  return ok;
}

inline bool
OCFESLV::warmstart_broadcast
( FFVar const& state, std::vector<double> const& ic_profile, double* var )
const
{
  // Broadcast a spatial IC profile across EVERY evolution node of @p state in var[]:
  // a constant-in-evolution initial guess for the next marching element's solve.  The
  // profile is indexed like terminal_profile() (spatial domains in var_dom order,
  // element-major / node-inner, first spatial domain innermost).
  if( var == nullptr ) return false;
  if( !_evolution_dom_set || !_evolution_dom_var.dag() ) return false;
  auto itv = _mVar.find( state );
  if( itv == _mVar.end() ) return false;
  std::vector<FFVar> const doms( itv->second.begin(), itv->second.end() );
  auto const evo_id = _evolution_dom_var.id();

  size_t const nd = doms.size();
  std::vector<size_t> nel( nd ), nnode( nd );
  std::vector<bool>   is_sp( nd );
  size_t sp_total = 1;
  for( size_t j = 0; j < nd; ++j ){
    FFDom const& d = _mDom.at( doms[j] );
    nel[j]   = d.n_elem;
    nnode[j] = d.n_node;
    is_sp[j] = !( doms[j].id() == evo_id );
    if( is_sp[j] ) sp_total *= nel[j] * nnode[j];
  }
  if( ic_profile.size() != sp_total ) return false;   // shape guard

  size_t nblk = 1; for( size_t j = 0; j < nd; ++j ) nblk *= nel[j];
  size_t blk  = 1; for( size_t j = 0; j < nd; ++j ) blk  *= nnode[j];

  for( size_t lin = 0; lin < nblk; ++lin ){
    std::map<FFVar,size_t,lt_FFVar> ndx_el;
    std::vector<size_t> eidx( nd );
    size_t rem = lin;
    for( size_t j = 0; j < nd; ++j ){
      eidx[j] = ( nel[j] ? rem % nel[j] : 0 );
      ndx_el[ doms[j] ] = eidx[j];
      if( nel[j] ) rem /= nel[j];
    }
    size_t const pos = pos_state( state, ndx_el );
    for( size_t nk = 0; nk < blk; ++nk ){
      // within-block node multi-index (first domain innermost == pos_state layout)
      std::vector<size_t> nidx( nd );
      size_t r = nk;
      for( size_t j = 0; j < nd; ++j ){
        nidx[j] = ( nnode[j] ? r % nnode[j] : 0 );
        if( nnode[j] ) r /= nnode[j];
      }
      // spatial-flat index (first spatial domain innermost, element-major/node-inner)
      size_t sflat = 0, stride = 1;
      for( size_t j = 0; j < nd; ++j ) if( is_sp[j] ){
        size_t const k = eidx[j] * nnode[j] + nidx[j];
        sflat  += k * stride;
        stride *= nel[j] * nnode[j];
      }
      var[ pos + nk ] = ic_profile[ sflat ];
    }
  }
  return true;
}

inline bool
OCFESLV::reset_evolution_window
( double t0, double w )
{
  // Slide a ONE-element evolution template to the window [t0, t0+w] in place, with NO
  // re-setup().  Both the origin and the width move, so a symbolic forcing/BC (a DAG
  // function of the evolution coordinate) auto-evaluates at the correct ABSOLUTE window
  // -- this is what lets solve()-marching advance a non-autonomous problem without
  // re-seeding.  The residual/Jacobian differentiate live through the OCVar path
  // (_eval_eqn reads _mDom via _node_block), so the new window takes effect on the next
  // solve().  Classification / index reduction / interface plan are geometry-independent
  // and untouched.  Requires options.SOLVE.LINEAR_EVAL_CACHE == LINEAR_NONE (the default).
  if( !_issetup || !_evolution_dom_set || !_evolution_dom_var.dag() ) return false;
  if( !( w > 0. ) ) return false;
  auto itd = _mDom.find( _evolution_dom_var );
  if( itd == _mDom.end() ) return false;
  FFDom& d = itd->second;
  if( d.n_elem != 1 || d.elem_bnd.size() < 2 || d.elem_len.empty() )
    return false;                                  // single-element template only

  // No-op when the window is unchanged (uniform/repeated step): preserves node cache
  // and any linear-eval cache.
  double const cur0 = d.elem_bnd.front();
  double const curw = d.elem_len.front();
  double const tol  = 64. * DBL_EPSILON
                    * std::max( 1., std::max( std::fabs( t0 ), std::fabs( w ) ) );
  if( std::fabs( cur0 - t0 ) <= tol && std::fabs( curw - w ) <= tol ) return true;

  // In-place mutation only (OCFESLV is friend FFDom): never reconstruct the FFDom -- its
  // constructor does not build the collocation-node cache, so a rebuilt domain would
  // have empty lgnodes() and crash downstream.  In-place keeps the node-complete domain.
  d.elem_bnd[0] = t0;
  d.elem_bnd[1] = t0 + w;
  d.elem_len[0] = w;
  d.w_elem      = w;
  d.lo_dom      = t0;
  d.up_dom      = t0 + w;

  // _node_block() reads _nodeBlockCache (returns copies), so refresh this domain's cached
  // node block in place -- cannot dangle any caller.
  auto itnc = _nodeBlockCache.find( _evolution_dom_var );
  if( itnc != _nodeBlockCache.end() && !itnc->second.empty() )
    itnc->second[0] = d.lgnodes( d.elem_lo( 0 ), d.elem_up( 0 ) );

  _linEqnReady = _linFctReady = false;              // force live residual eval
  return true;
}

inline bool
OCFESLV::reset_evolution_width
( double w )
{
  // Change the width, keeping the current origin (a convenience over
  // reset_evolution_window()).
  if( !_issetup || !_evolution_dom_set || !_evolution_dom_var.dag() ) return false;
  auto itd = _mDom.find( _evolution_dom_var );
  if( itd == _mDom.end() || itd->second.elem_bnd.empty() ) return false;
  return reset_evolution_window( itd->second.elem_bnd.front(), w );
}



} // namespace mc

// ---------------------------------------------------------------------------
// v2 interface-plan decision layer (Phase 1), consolidated inline (formerly
// ocfeslv_refactor.hpp).  Placed after the OCFESLV / OCPlan definitions so the
// types are complete.  Defines OCFESLV::_apply_interface_edge_flags_v2: every
// admissible duplicated DOF that is NOT a genuine boundary-condition row
// becomes an exact, tau-able value-continuity constraint; coupling/source are
// left exactly as the caller computed them, so IC_WEAK stays byte-identical.
// The BC-role classification mirrors the trace_hard_physical_row() lambda so
// genuine boundary rows are still protected.
// ---------------------------------------------------------------------------
#include <cfloat>
#include <cmath>

namespace mc {


inline bool
OCFESLV::_setup_deriv_cache
()
{
  _linEqnReady = _linFctReady = false;
  _linEqn.clear();
  _linFct.clear();
  _linCstVal.clear();
  return _build_jacobian_pattern();
}


inline OCFESLV::STATUS
OCFESLV::solve
( std::vector<FFModel::InputVal> const& vIn )
{
  if( !_issetup ) return FATAL;
  std::vector< std::pair< FFVar, std::vector<double> > > given;
  std::vector<double> cst( var_constant().size(), 0. );
  std::string err;
  bool const ok = resolve_values( vIn,
    [&]( FFVar const& u, std::vector<double> const& v, std::string& e ){
      size_t const n = input_dofs( u ).size();
      if( !n ){ e = u.name() + " is neither an input nor a constant of the model"; return false; }
      if( fixed_input( u ) ){ e = u.name() + " is neither an input nor a constant of the model (it is fixed)"; return false; }
      if( v.size() != n ){ e = u.name() + " has " + std::to_string( n ) + " DOFs, but " + std::to_string( v.size() ) + " values were given"; return false; }
      given.emplace_back( u, v );  return true; },
    [&]( size_t k, double v ){ cst[k] = v; },
    err );
  if( !ok ) throw std::invalid_argument( "OCFESLV::solve ** " + err + "\n" );
  std::vector<double> inp;
  if( !init( _valSta, inp, cst.empty()? nullptr: cst.data() ) ) return FATAL;
  for( auto const& [u, v] : given )
    if( !set_input_values( u, v, inp.empty()? nullptr: inp.data() ) )
      throw std::invalid_argument( "OCFESLV::solve ** the values of " + u.name() + " could not be set\n" );
  // the solve fills the outputs itself (_valFct), with captures or accumulated over a march: never re-evaluated here
  _lastReport = solve( _valSta.data(), inp.empty()? nullptr: inp.data(), cst.empty()? nullptr: cst.data() );
  return _lastReport.converged? NORMAL: FAILURE;
}

} // namespace mc

#if defined(_WIN32)
# pragma pop_macro("INTERFACE")
#endif

#endif
