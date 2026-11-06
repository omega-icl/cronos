// Copyright (C) Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the EPL-2.0 with GPL-2.0-or-later as a Secondary License; see the LICENSE file.
// SPDX-License-Identifier: EPL-2.0 OR GPL-2.0-or-later

#ifndef CRONOS__OCFEBASE_HPP
#define CRONOS__OCFEBASE_HPP

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
//! @brief Direction-agnostic collocation machinery shared by the OCFE solvers.
////////////////////////////////////////////////////////////////////////
//! mc::OCFEBASE holds what sits between MC++'s collocation primitives (OCBase, BASE_OC, FFDom) and a
//! solver's time assembly: the OCVar environment and the evaluation kernels, the interface plan
//! (construction, frozen realisations, audits), imposition, the determinacy verdict core, and the
//! options.  mc::OCFESLV (orthogonal collocation in every direction) and, in the method-of-lines work,
//! mc::OCFEMOL (collocation in space, an integrator in time) derive from it.  It still reads the
//! evolution direction through FFModel; parameterising the directions it plans over is the first step
//! of the MoL phase that needs it.  Extracted from ocfeslv.hpp in eight relocation commits, each gated
//! on the OCFESLV corpus being byte-identical.
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
class OCFEBASE
////////////////////////////////////////////////////////////////////////
: public FFModel, public OCBase
{
  template <typename U> friend class OCVar;

public:
  //! @brief Constructor: the derived solver hands its own `options` member, to which OCFEBASE::options is bound.
  //! (Forming a reference to a not-yet-constructed member is well-formed; nothing here reads through it.)
  struct Options;                       // defined below; without this, `Options` here would name FFModel::Options
  explicit OCFEBASE( Options& opts ) : options( opts ) {}
  virtual ~OCFEBASE() = default;

public:
  using OCBase::VarType;
  using OCBase::Exceptions;
  using OCBase::t_Dom;

  typedef std::map< FFVar, std::tuple< size_t, size_t, std::vector<size_t> >, lt_FFVar >               t_Blk;
  typedef std::map< long, std::pair< int, size_t > >                                                   t_Typ;
  typedef FADType<double>                                                                              t_FAD;

protected:
  //! @brief Map of dependencies against variable type and index
  t_Typ                              _mDep;

  //! @brief vector of constant dependencies
  std::vector<FFDep>                 _depCst;

  //! @brief vector of constant collocations
  std::vector<OCVar<double>>         _OCCst;

  //! @brief vector of constant collocations
  std::vector<OCVar<FFDep>>          _OCDepCst;

  //! @brief vector of domain dependencies
  std::vector<FFDep>                 _depDom;

  //! @brief vector of domain collocations
  std::vector<OCVar<double>>         _OCDom;

  //! @brief vector of domain collocations
  std::vector<OCVar<FFDep>>          _OCDepDom;

  //! @brief Cache for row-specific auxiliary-direction expression lookups.
  mutable std::map<std::tuple<size_t,int,size_t>,bool> _eqnAuxDirCache;

  //! @brief Map from input to the FFDom of each of its domains: the domain's grid with the input's declared
  //! collocation applied.  Derived at setup from FFModel's declaration; also pinned per window when marching.
  typedef std::map< FFVar, t_Dom, lt_FFVar >  t_Coll;

  //! @brief optional per-input collocation discretisations. Keys are input variables, values are maps from physical domain variables to the collocation rule used to parameterise that input along that domain. Missing entries mean: use the state/domain discretisation in _mDom.
  t_Coll                             _mInpDisc;

  //! @brief Moving-mesh (r-adaptive) specification for one domain direction.
  //! A direction registered through OCFESLV::add_adaptive_domain() is discretised on a fixed uniform
  //! computational mesh xi in [0,1] (its registered FFDom); a solved mesh state x(xi,t) maps xi onto
  //! the physical extent [lo_dom, up_dom] of that same FFDom, and follows the equidistribution
  //! d/dxi(M x_xi) = 0 for the monitor M.  Only an explicit monitor expression is supported.
  //! Solver-side type, held in the model's UsrData because the declaration carries it; the discretisation
  //! it describes is applied here.
  struct MeshAdaptSpec {
    FFVar monitorExpr;  //!< explicit monitor M(xi,t); used iff hasMonitorExpr
    bool  hasMonitorExpr = false;  //!< must be true: there is no built-in monitor
  };

  //! @brief Declared moving-mesh (r-adaptive) specifications per domain direction, as given to
  //! add_adaptive_domain(); the physical extent is read from the direction's FFDom.
  std::map<FFVar,MeshAdaptSpec,lt_FFVar>  _adaptDomUsr;

  //! @brief vector of input dependencies
  std::vector<FFDep>                 _depInp;

  //! @brief vector of input collocations
  std::vector<OCVar<double>>         _OCInp;

  //! @brief vector of input collocations
  std::vector<OCVar<FFDep>>          _OCDepInp;

  //! @brief vector of input block data
  t_Blk                              _blkInp;

  //! @brief Stage 1: k = n_trace_var - rank(B) as computed from the DRAFT by
  //! _build_frozen_trace_projection.  Cross-checked against the value
  //! interface_rank_symmetry computes from the COMMITTED plan; a disagreement means the
  //! two are seeing different blocks and the projection would act on the wrong one.
  //! @brief (state, direction, element, node) triples whose C0 claim was suppressed by
  //! FIX 1 or FIX 3.  Keyed exactly as a claim is keyed, with no equation-row component,
  //! so a claim emitted for a suppressed triple is detectable as a leak.
  //! MUTABLE: the claim-emission scan that populates this runs in a const context, and
  //! this is a diagnostic-only record that no behaviour depends on.
  mutable std::set<std::tuple<size_t,size_t,size_t,size_t>> _pinSuppressedTriples;

  //! @brief Per-state counts of C0 claims suppressed by FIX 1 (pointwise pin) and FIX 3
  //! (closure over jointly-pinned states), summarised once per plan build.  NOTE these are
  //! per equation-row VISIT, not per claim: a claim may be emitted from several rows, and
  //! FIX 3 RELOCATES a claim rather than removing it, so the claim count does not move.
  mutable std::map<size_t,size_t> _pinSuppressCount1, _pinSuppressCount3;

  size_t                             _traceProjK    = 0;
  size_t                             _traceProjNTau = 0;
  //! @brief Stage 1: true once C_tilde=(I-P)C and the P block have been written into the
  //! plan's continuity rows.  The rank-symmetry audit compares against rank(C_tilde)=r
  //! rather than n_trace_var once this is set, since the projected system enforces r
  //! conditions by design rather than by deficiency.
  bool                               _traceProjApplied = false;
  size_t                             _traceProjRank    = 0;

  //! @brief one-shot guard for the `[stage2-preflight]` line
  //! emitted from the _eval_eqn row-count audit.
  //!
  //! MUTABLE: _eval_eqn() is const and the flag records only that a diagnostic line has
  //! already been printed for this OCFESLV instance; no behaviour depends on it.  Deliberately
  //! NOT propagated by deep_copy_from(): a copied environment is a different process-visible
  //! object and should emit its own line, and a duplicate diagnostic line is harmless where
  //! a missing one is not.
  mutable bool                       _stage2PreflightReported = false;

  //! @brief receiver row -> the explicit tau column(s) injected into it (plan half).
  //!
  //! FULLY QUALIFIED, deliberately.  This block sits ~115 lines ABOVE the
  //! `using FrozenPhysicalResidualRowKey = OCPlan::FrozenPhysicalResidualRowKey;`
  //! alias, so the short name is not in scope here.  OCPlan is complete at this
  //! point (ocplan.hpp is included before the class), and OCFESLV is a friend of
  //! OCPlan, so the nested type is reachable.  Writing it out costs nothing and
  //! removes any dependence on where the alias block happens to sit -- moving
  //! this declaration and moving the alias are then independent edits.
  mutable std::map<OCPlan::FrozenPhysicalResidualRowKey,
                   std::vector<std::pair<size_t,double>>> _tauHostRows;
  //! @brief pivot index -> (tau column -> its coefficient in that pivot row).
  //! Needed to evaluate the DELIVERED coefficient: a surviving row emits
  //! R_k - sum_p y_p * R_p, so a tau's net coefficient there is its coefficient
  //! in row k MINUS y_p times its coefficient in each pivot row p.  Structural
  //! survival does not imply a non-zero net (measured).
  mutable std::map<size_t,std::map<size_t,double>>  _tauPivotCoef;
  //! @brief tau column -> largest |net delivered coefficient| over surviving rows.
  mutable std::map<size_t,double>                   _tauNetMax;
  //! @brief tau column -> number of times the eval-time injection actually executed.
  mutable std::map<size_t,size_t>            _tauInjectCount;
  //! @brief tau column -> host rows that survived LOOP 2 / were dropped as pivot / as exact.
  mutable std::map<size_t,size_t>            _tauFateEmit;
  mutable std::map<size_t,size_t>            _tauFatePivot;
  mutable std::map<size_t,size_t>            _tauFateExact;
  //! @brief tau column -> pivot indices its dropped host rows occupy.
  mutable std::map<size_t,std::set<size_t>>  _tauPivotOfCol;
  //! @brief pivot index -> number of SURVIVING rows whose Schur correction references it.
  //! This is the "after Schur substitution" half: a dropped pivot row is not necessarily
  //! lost, because its value is subtracted into every corrected row.  Zero here means the
  //! coupling really did leave the system.
  mutable std::map<size_t,size_t>            _tauCorrRefs;
  //! @brief one-shot guard for the [taucol/assembly] report.
  mutable bool                               _tauColAssemblyReported = false;
  //! @brief set by the PLAN half once _tauHostRows is populated.
  //!
  //! ORDER GUARD, and it is load-bearing.  _audit_strong_interface_rows() runs a dry
  //! FFDep pass through _eval_eqn -- twice from setup() and once from deep_copy_from()
  //! -- and at least one of those calls precedes the plan half in source order.  Without
  //! this flag the assembly recorders fire on a pass where _tauHostRows is still empty,
  //! print a report that observed nothing, and LATCH _tauColAssemblyReported, suppressing
  //! the real one.  That is a diagnostic producing plausible output from code that was not
  //! measuring anything -- the failure mode this instrument exists to avoid.
  //!
  //! Requiring it makes the call ORDER NON-LOAD-BEARING: move the plan half anywhere
  //! after the commit and the assembly half still attaches to the first eval that follows.
  mutable bool                               _tauColPlanReady = false;

  //! @brief vector of state dependencies
  std::vector<FFDep>                 _depVar;

  //! @brief vector of state collocations
  std::vector<OCVar<double>>         _OCVar;

  //! @brief vector of state collocations
  std::vector<OCVar<FFDep>>          _OCDepVar;

  //! @brief vector of state block data
  t_Blk                              _blkVar;

  //! @brief vector of equation dependencies
  std::vector<FFDep>                 _depEqn;

  //! @brief equation/output collocation workspace
  OCVar<double>                      _OCEqn;

  //! @brief equation/output dependency collocation workspace
  OCVar<FFDep>                       _OCDepEqn;

  //! @brief equation/output forward-AD collocation workspace
  mutable OCVar<t_FAD>               _OCFEqn;

  //! @brief Number of generated trace/tau unknowns appended after state coefficients.
  size_t                             _nTraceVar = 0;

  //! @brief Offset of generated trace/tau unknowns inside the var[] vector.
  size_t                             _traceVarOffset = 0;

  //! @brief Number of setup-retained exact trace constraints.
  size_t                             _nTraceConstraint = 0;

  //! @brief Sealed frozen interface plan consumed by _eval_eqn / _eval_fct.
  //!
  //! All plan data (receiver edges, weak-SAT terms, exact trace rows,
  //! IC_STRONG projection, plan dimensions) live in this object.  OCFESLV
  //! builder methods populate a transactional OCPlan::InterfacePlanDraft,
  //! validate it via _plan.validate_draft(), and commit it atomically via
  //! _plan.commit_draft().  _eval_eqn() reads the plan exclusively through
  //! OCPlan's const query interface.
  mutable OCPlan _plan;

  //! @brief True only during OCFESLV::_build_interface_plan_draft().
  //!
  //! Kept as a separate OCFESLV flag (not inside OCPlan) so that plan
  //! builders — which are OCFESLV private methods — can assert the build
  //! is active without exposing a mutable setter through OCPlan's public API.
  mutable bool _interfacePlanBuildAllowed = false;

  // Bring OCPlan nested types into OCFESLV's protected scope so that all
  // private builder methods and inline implementations can use the short
  // names without qualification, preserving maximum source compatibility.
  using InterfacePlanState           = OCPlan::State;
  using t_FrozenInterfaceClaimKey    = OCPlan::t_FrozenInterfaceClaimKey;
  using t_InterfaceClaimKey          = OCPlan::t_InterfaceClaimKey;
  using InterfaceReceiverEdge        = OCPlan::InterfaceReceiverEdge;
  using FrozenWeakSatTerm            = OCPlan::FrozenWeakSatTerm;
  using FrozenTraceConstraintTerm    = OCPlan::FrozenTraceConstraintTerm;
  using FrozenPhysicalResidualRowKey = OCPlan::FrozenPhysicalResidualRowKey;
  using FrozenStrongTauElimination   = OCPlan::FrozenStrongTauElimination;
  using InterfacePlanDraft           = OCPlan::InterfacePlanDraft;
  using InterfaceEdgeSource          = OCPlan::InterfaceEdgeSource;

public:
  //! @brief One recorded coupling decision: the value, whether it was RESCUED (source is not the
  //! principal symbol, i.e. the number was substituted rather than read from the model), and how many edges.
  struct t_CouplingRecord { double value = 0.; bool refused = false; size_t n = 0; };

protected:
  mutable std::map<std::pair<size_t,size_t>,t_CouplingRecord> _couplingLog;

  //! @brief Two-pass IC_STRONG DROP (explicit-tau hybrid): redundant
  //! (column-dependent / implied) continuity claim keys found in pass 1 and
  //! excluded from allocation in pass 2.  Empty for every column-rank-full plan
  //! (PDE2/3/4), so the two-pass is a no-op there.
  mutable std::set<t_FrozenInterfaceClaimKey> _strong_dropped_claims;
  //! @brief C2b: claims whose rescued receivers must NOT be silenced under EXACT_NATURAL, because the
  //! redundancy pass found their multiplier column dependent with only the natural receiver -- restored on the
  //! two-pass re-run ("silence where rank allows").  Cleared with _strong_dropped_claims.
  mutable std::set<t_FrozenInterfaceClaimKey> _exact_natural_restored;
  //! @brief 2026-10-16: claims whose natural receiver was PROMOTED from a defining row by the per-node matching.
  //! Natural for the pivot, NOT for the exact-natural silencing of their other edges: on OCFE_zchain silencing
  //! those edges removed the projector's support (9 PASS -> 3), so a promoted claim keeps its duplicates live.
  mutable std::set<t_FrozenInterfaceClaimKey> _promotedNatural;

  //! @brief Why each IC_STRONG claim was kept as an EXPLICIT tau rather than eliminated.
  //!
  //! 1 = verify-protected AND rank-extending  (item 10 step 3, DROP_REDERIVE): a needed
  //!     continuity the receiver-B proxy mis-flagged as redundant.  Eliminating it through
  //!     the Schur would demand a corner copy the joint tau/geometric-C0 assignment cannot
  //!     supply, so it is realised the IC_TRACE way instead.
  //! 2 = every replacement candidate is an auxiliary LINK row: exact-replacing would consume
  //!     that auxiliary's defining equation and turn the auxiliary into a null direction.
  //! 3 = verify-protected but NOT rank-extending (same reason as 1, different branch).
  //!
  //! NOTE: none of these is "no pivot was available" -- an earlier version of the census
  //! said that and it was wrong.  These are DELIBERATE keep-explicit routes chosen to avoid
  //! a worse failure elsewhere.  What is NOT deliberate is that the promised "explicit tau
  //! unknown PLUS a real continuity row" appears to deliver the column but not the row.
  mutable std::map<t_FrozenInterfaceClaimKey,int> _keepExplicitReason;
  //! @brief Set when pass 1 identifies droppable redundant claims; signals
  //! _build_interface_plan_draft to re-run the pipeline excluding them.
  mutable bool                       _strong_redundancy_pending = false;
  //! @brief A RESTORE rebuild is pending (2026-10-06, WORKPLAN 1.5): EXACT_NATURAL returned silenced claims to their
  //! natural receivers and the plan must be re-derived with them -- no claim dropped.  Booked on its own retry, so
  //! that a DROP (a genuinely redundant claim excluded) can still follow it: a model whose pass 0 restores and whose
  //! pass 1 flags redundant claims (OCFE_zindex's chain, high-index in z, under IC_STRONG) used to exhaust the single
  //! DROP retry on the restore and fail at the exclusion it needed.
  mutable bool                       _strong_restore_pending    = false;

  //! @brief one-shot guard so the promotion requester cannot re-fire on its own
  //! rebuild pass.  Cleared once per `_build_interface_plan_draft` call, alongside
  //! `_strong_dropped_claims`.
  mutable bool                       _promotionRequested = false;

  //! @brief state retained so the post-solve audit can evaluate the
  //! violation IC_TRACE absorbs into beta.  Populated ONLY under IC_TRACE with the projection
  //! applied and CRONOS_AUDIT_BETA set, so production pays nothing.
  //!
  //! THE IDENTITY BEING MEASURED.  The applied form appends eps*P coefficients to the
  //! continuity rows, so row i reads  (C x)_i + eps*(P lambda)_i + r_c,i = 0.  With P = Z Z^T
  //! and Z^T Z = I, projecting on Z^T gives
  //!
  //!     eps * beta  =  - Z^T ( C x + r_c )  =  - Z^T r_c^raw ,      beta = Z^T lambda
  //!
  //! i.e. `eps*||beta||` IS the norm of the UNPROJECTED continuity violation lying in
  //! null(B) -- the part of the claim set IC_TRACE does not enforce but absorbs.  Roadmap B5
  //! proves that block imposes nothing on x; this measures how much it is carrying.
  //!
  //! FALSIFIABLE SELF-TEST, and it should be run before the number is believed: sweep
  //! INTERFACE.TRACE_PROJ_EPS.  `eps*||beta||` must be eps-INDEPENDENT (it equals a fixed violation)
  //! while `||beta||` scales as 1/eps.  If it does not, the instrument is wrong, not the model.
  mutable arma::mat                  _betaZ;         //!< m x k null basis of B
  mutable std::vector<size_t>        _betaColOf;     //!< dense tau index -> var[] column
  mutable double                     _betaEps = 0.;  //!< eps actually applied
  mutable bool                       _betaReady = false;

  //! @brief per-null-direction class, CROSS-domain (true) or intra-domain (false),
  //! so the absorbed violation can be split between the two kernel families.
  //!
  //! WHY THE SPLIT IS WELL DEFINED.  A basis of null(B) is not unique, so "which direction"
  //! is generally basis-dependent.  Here it is not: Z is built COMPONENT-WISE over the
  //! connected components of the claim/receiver-row bipartite graph, so every column has
  //! support inside ONE component, and the cross/intra label is a property of the component
  //! rather than of the basis chosen inside it.
  //!
  //! MEASURED DECOMPOSITION (MMPDE13, confirmed by forward prediction at two t-meshes):
  //!     intra = (nel_xi - 1)                     -- NEL_T-independent; the lifting-alias
  //!                                                 family that PDE2/PDE5/PDE9 also carry
  //!     cross = 2*(nel_xi - 1)*(NEL_T - 1)       -- two per interior t-seam node; no
  //!                                                 counterpart in any fixed-mesh control
  mutable std::vector<char>          _betaDirCross;

  //! @brief Continuity claim keys the redundancy detector must NOT drop, kept
  //! explicit instead (item 10 step 3, DROP_REDERIVE).  Populated by the driver
  //! from verify_interface_drop's over-drop set between re-derive iterations and
  //! consumed by _detect_redundant_continuity_claims.  Unlike _strong_dropped_-
  //! claims it is NOT cleared at build start (it must persist across the
  //! re-setup of a re-derive sequence); it is cleared only on a full reset /
  //! deep_copy / via clear_protected_interface_claims().
  mutable std::set<t_FrozenInterfaceClaimKey> _strong_protected_claims;
  //! @brief IC_STRONG taus the trace projector reaches; kept explicit on the promotion rebuild.
  mutable std::set<t_FrozenInterfaceClaimKey> _strong_promoted_claims;
  mutable bool                                _strong_promotion_pending = false;

  //! @brief Structural determinacy numbers of the last applied GLOBAL plan build (C4): the tau-block rank
  //! deficit the builder predicts, and the cluster count.  Measured corpus-wide to equal the
  //! SPQR defR on 96 of 101 programs and to over-report on the other 5 -- never the reverse.
  mutable long _structDeficit = -1;
  mutable long _structDeficientClusters = -1;

  //! @brief per state, [0] edges whose coupling was DERIVED from the receiving row's constant c0,
  //! [1] edges REFUSED because that coefficient is not constant (forced value kept, edge labelled undefined).
  mutable std::map<size_t,std::array<size_t,2>> _rescueC2;
  //! @brief per state, WHY C2 refused -- 1 the receiving row is in no symbol table (boundary or
  //! demoted), 2 that table was not computed, 3 the state is not one of its columns, 4 the coefficient is not
  //! constant, 5 it is the CONSTANT ZERO (the state enters the row only differentiated: no c0 exists).
  mutable std::map<size_t,std::array<size_t,8>> _rescueC2Why;
  //! @brief per state, for each zero-absent refusal, what the receiving row contains -- whether the
  //! claimed state occurs in it at all, in which directions it is differentiated there, and the claim's own
  //! direction.  Report-only; it decides whether those edges are a coupling question or a receiver-graph one.
  mutable std::map<size_t,std::map<std::string,size_t>> _rescueC2Zero;

  //! @brief Per-state receiver-edge audit for the SAT-coupling path (diagnostic only).
  //!
  //! A receiver edge acquires a nonzero coupling in one of three ways: from the block
  //! principal symbol (states present in sym.vState), or -- when _sat_coupling returns
  //! numerically ZERO -- by one of two rescues, auxiliary_algebraic_receiver_state or
  //! geometric_c0_receiver_state, each of which forces coupling = 1.  If the SAT coupling
  //! is zero and NEITHER rescue fires, the edge keeps coupling = 0 and its DEFAULT label
  //! PRINCIPAL_SYMBOL: a null row in the augmented interface system, i.e. exactly the
  //! "tier-0 Newton singular on the augmented interface system" failure.
  //!   Counters per state id: [0] edges on the SAT path, [1] SAT coupling was zero,
  //! [2] rescued as AUXILIARY_ALGEBRAIC, [3] rescued as GEOMETRIC_C0, [4] LEFT AT ZERO.
  //! a[4] > 0 is the defect.
  mutable std::map<size_t,std::array<size_t,5>> _ifaceEdgeAudit;
  mutable std::map<size_t,double>               _ifaceEdgeMaxSat;

  //! @brief Work vector for collocation propagataion through DAG
  std::vector<OCVar<double>>         _wrkOC;

  //! @brief Work vector for collocation propagataion through DAG
  std::vector<OCVar<FFDep>>          _wrkOCDep;

  //! @brief Forward-AD state collocation cache used by deriv()
  mutable std::vector<OCVar<t_FAD>>  _OCFVar;

  //! @brief Forward-AD input collocation cache used by deriv()
  mutable std::vector<OCVar<t_FAD>>  _OCFInp;

  //! @brief Forward-AD domain collocation cache used by deriv()
  mutable std::vector<OCVar<t_FAD>>  _OCFDom;

  //! @brief Forward-AD constant collocation cache used by deriv()
  mutable std::vector<OCVar<t_FAD>>  _OCFCst;

  //! @brief Forward-AD equation residual values used by deriv()
  mutable std::vector<t_FAD>         _FEqn;

  //! @brief Forward-AD state seeds used by deriv()
  mutable std::vector<t_FAD>         _FVar;

  //! @brief Forward-AD input seeds used by deriv()
  mutable std::vector<t_FAD>         _FInp;

  //! @brief Lightweight dependency-update plan used by eval() to avoid map lookups in every element block.
  struct t_DepUpdate
  {
    int                              typ = 0;
    size_t                           pos = 0;
  };
  std::vector< std::vector<t_DepUpdate> > _depPlanEqn;
  std::vector< std::vector<t_DepUpdate> > _depPlanFct;

  //! @brief Per-equation evaluation domain: the distributed directions each equation's residual
  //! touches (state/input domains referenced by its subgraph).  Equation analogue of _fctEvalDom.
  //! For a domain-less equation these are the FFIntegral-reduced ("free") directions accumulated
  //! over in _eval_eqn LOOP 2 (section 4).
  std::vector< std::set<FFVar,lt_FFVar> > _eqnEvalDom;
  //! @brief Per-equation OpEval (FFEval) point reductions: direction -> target physical coordinate.
  std::vector< std::map<FFVar,double,lt_FFVar> > _eqnEvalPoint;
  std::vector< std::map<FFVar,int,lt_FFVar> > _eqnEvalSide;   //!< per equation, the SIDE of each _eqnEvalPoint direction

  //! @brief Cached physical node coordinates by domain variable and element.
  std::map<FFVar,std::vector<std::vector<double>>,lt_FFVar> _nodeBlockCache;

public:
  //! @brief Post-index-reduction structural DOF audit (max-bipartite matching on
  //! the reduced collocation sparsity).  Filled by _audit_reduced_dof() at end of
  //! setup(); a nonzero deficiency is a genuine structural rank deficiency.
  struct t_DofAudit
  {
    bool   ran        = false;  //!< audit executed (reduction happened)
    bool   trace      = false;  //!< ran on an interface-trace augmented pattern (advisory verdict)
    size_t rows       = 0;      //!< collocation equation rows
    size_t cols       = 0;      //!< variable columns (max sparsity index + 1)
    size_t rank       = 0;      //!< structural (max-matching) rank
    size_t deficiency = 0;      //!< min(rows,cols) - rank  (0 => full structural rank)
    bool   square     = false;  //!< rows == cols
    bool   ok() const { return ran && square && deficiency == 0; }
    void   clear() { *this = t_DofAudit{}; }
  };

protected:
  //! @brief Cached structural reduced dependency pattern (reduced_dependency_pattern()).  Computed
  //! once on first access and reused by reference; invalidated by setup(), register_control(),
  //! clear_controls() and deep_copy_from() -- anything that changes the model structure or the
  //! control registry.
  std::vector<std::vector<bool>>            _reducedDepPattern;
  bool                                      _reducedDepValid = false;

  //! @brief Cached residual-row provenance: _rowEqnNdx[out_row] = source equation index into _mEqn
  //! (or -1 for plan-owned interface/trace rows).  Built once by _eval_eqn (the frozen plan makes it
  //! eval-invariant) behind _rowProvenanceReady; invalidated by setup()/deep_copy_from().  Used to
  //! resolve each marching-transfer state's INITIAL-closure rows for the k>0 continuity override.
  mutable std::vector<int>                  _rowEqnNdx;
  //! @brief Row provenance, the node half (MoL phase 3, P3f): for each emitted row, the collocation node it stands at --
  //! per direction of the equation's domain, (direction id, element, node within the element, the mask limit the
  //! direction was collocated under).  Empty for a plan-owned row (_rowEqnNdx == -1).  Built with _rowEqnNdx, on the
  //! first evaluation.  OCFEMOL places Mc's entries with it: the row's equation has a mass entry for a field, the field's
  //! value at THIS node is the column.
  struct t_RowNode
  {
    struct Dim { long dir_id; size_t iel; size_t node; int lim; };
    std::vector<Dim> dims;
  };
  mutable std::vector<t_RowNode>            _rowNode;
  //! @brief Row provenance, the fold half (MoL phase 5): an emitted IC_STRONG row is its own row MINUS a combination of
  //! pivot rows the exact elimination dropped -- value - sum_p coeff_p * pivot_p.  For each emitted row, the pivots folded
  //! in: (their equation index, their node, coeff_p).  Empty for an unfolded row.  OCFEMOL folds Mc the same way: the
  //! plan's row operations act on Rc, and Mc must receive the same operations.
  struct t_RowFold { int eqn; t_RowNode node; double coeff; };
  mutable std::vector< std::vector<t_RowFold> > _rowFold;
  mutable bool                              _rowProvenanceReady = false;
  //! @brief claim keys the legacy corner rule suppressed in the last receiver
  //! synthesis (recorded at the suppression point; no behavioural effect).  Packed into
  //! PlanInput::corner_suppressed_edges for the GLOBAL builder's shadow check.
  mutable std::set<t_InterfaceClaimKey>    _cornerSuppressedClaims;
  //! @brief claim keys suppressed as "pinned pointwise by an ALL-mask IC/BC"
  //! (never realised by a claim); packed as PlanInput::unrealised_claims.
  mutable std::set<t_InterfaceClaimKey>    _pinSuppressedClaims;
  //! @brief claims whose exact trace row the tensor-corner MGS rank filter REJECTED
  //! as linearly dependent on the accepted rows.  Homogeneous continuity rows, so dependence
  //! means the condition is IMPLIED by rows the plan imposes -- enforced, not missing.  Packed
  //! as PlanInput::rank_implied_claims and consumed by design §3.5.  NOT the same category as
  //! _pinSuppressedClaims (those are genuinely unenforced).
  mutable std::set<t_InterfaceClaimKey>    _rankImpliedClaims;
  //! @brief true iff the plan being materialised carries GLOBAL decisions
  //! (set when they are applied, cleared at each plan build).  The STRONG-side matcher
  //! keys on THIS, never on the environment (keying on the environment breaks legacy
  //! STRONG plans whenever CRONOS_PLAN_APPLY_GLOBAL is set).
  mutable bool                             _planIsGlobal = false;

  mutable bool                              _sampleOwnElement = false;   //!< sample_input: element-boundary nodes from their own element (transfer_terminal)

  t_DofAudit                                _dofAudit;

  // rev139 (uniform-rule step 1): W-verdict restore
  bool                                      _wdecideRetry = false;   //!< inside the one restore/revert re-derive
  mutable bool                              _auditQuiet = false;     //!< audit stores, defers printing (restore-attempt recursion)

  //! @brief Item-11 Part-2 per-(primitive,direction) symbol interface read.
  //! order/magnitude balance used by the symbol-AUTO resolver.  `order` is the
  //! detected operator order in the direction (1 or 2); `order_gt2` flags a
  //! nested aux chain (order>2, not handled — resolver must guard).  `via_aux`
  //! is true when order>=2 was detected through a reduce_order auxiliary link
  //! (i.e. a reduced direction).  `advection_dominated` is true when the
  //! top-order (diffusive) term is sub-grid (k_h < k* = adv/diff), i.e. the
  //! singular-perturbation corrector fires.
  struct t_DirRead {
    int    prim_col   = -1;
    int    dir_idx    = -1;
    int    order      = 0;
    bool   order_gt2  = false;
    bool   via_aux    = false;
    double adv   = 0.;
    double diff  = 0.;
    double kstar = std::numeric_limits<double>::quiet_NaN();
    double eps   = std::numeric_limits<double>::quiet_NaN();
    double kh    = std::numeric_limits<double>::quiet_NaN();
    bool   advection_dominated = false;
  };
  //! @brief Symbol interface reads per block id.  Populated read-only by
  //! _compute_interface_reads() at classify time; consumed by the symbol-AUTO
  //! resolver and the optional balance probe.  Always compiled.
  std::map<int, std::vector<t_DirRead>>     _blockInterfaceRead;

  //! @brief One cached principal-symbol SAT coupling entry.
  struct t_SatCoupling
  {
    size_t state_pos;   //!< Position of the coupled state in _vVar.
    long   dom_id;      //!< Identity (.id().second) of the coupled domain variable.
    double value;       //!< Classified coupling A_dom[eqn_row, state_col].
  };

  //! @brief Per-equation SAT coupling table, indexed by position in _mEqn.
  //! Built during setup after PDE classification has fixed the principal
  //! symbols.  It is then consumed by the frozen interface-plan builders and
  //! residual evaluation.  Empty when the model is not classified, in which
  //! case _sat_coupling() falls back to structural derivative incidence.
  mutable std::vector< std::vector<t_SatCoupling> > _satCoupling;
  mutable bool                                      _satCouplingBuilt      = false;
  mutable bool                                      _satCouplingClassified = false;
  mutable size_t                                    _satCouplingSerial     = 0;

protected:

  //! @brief WHAT OF THE EVOLUTION DIRECTION the current assembly collocates, and of what horizon -- as seen by the
  //! direction-agnostic code (the OCVar environment, the interface plan, the audits, the output evaluation), which must
  //! not read a solver's strategy: it asks this, and each derived class answers.  Three answers:
  //!   ALL  -- every element of the evolution domain (a monolithic solve: the default);
  //!   ONE  -- one element of a horizon solved element by element (OCFESLV's marching window);
  //!   NONE -- no element: t is a scalar input of the collocated system (OCFEMOL, the method of lines).
  //! (Was MarchContext / _march_context until 4 Nov 2026; NONE behaves as its "not marching" default did.)
  struct EvolutionWindow {
    enum Kind { ALL = 0, ONE, NONE } kind = ALL;
    size_t index      = 0;      //!< ONE: the collocated element's GLOBAL index in the horizon
    size_t n_elements = 0;      //!< the horizon's elements when it is solved element by element (0: it is not)
    bool   stepwise   = false;  //!< the solve proceeds element by element: requested, with n_elements >= 1
  };
  virtual EvolutionWindow _evolution_window
    () const
    { return EvolutionWindow(); }
  //! @brief Whether a solve element by element along the evolution direction is REQUESTED -- known before any
  //! horizon is recorded, unlike EvolutionWindow::stepwise.  Default: no.
  virtual bool _stepwise_requested() const { return false; }
  //! @brief Whether output @p f is distributed over the evolution direction: it spans the whole horizon.
  bool _output_over_evolution
    ( t_Fct const& f )
    const
    {
      if( f.kind != FctKind::DISTRIBUTED || !_evolution_dom_var.dag() ) return false;
      for( auto const& [dv, lim] : f.grid ){ (void)lim; if( dv.id() == _evolution_dom_var.id() ) return true; }
      return false;
    }
  //! @brief For an output distributed over the evolution direction: the rows horizon element @p k produces, placed in
  //! the full-horizon layout (the derived class that solves element by element answers).  Default: none.
  virtual std::vector<size_t> _window_output_rows( t_Fct const& f, size_t k ) const { (void)f; (void)k; return {}; }

public:
  //! @brief the SOLVER's options.  DERIVES from FFModel::Options, which brings the model's 12
  //! fields, its 14 enums and 35 enumerators, and its static environment helpers along unchanged -- so
  //! Options::IC_WEAK here, OCFESLV::Options::RED_FULL in a driver, and oc.options.REDUCE.ORDER in all 138
  //! driver files keep resolving exactly as before, with no migration.
  //!
  //! The inherited part is the REQUESTED model configuration; FFModel::options is the APPLIED one, and
  //! setup() is the single moment of transfer (a slicing assignment).  A model option changed after setup()
  //! visibly does nothing until the next setup(), rather than taking effect mid-solve as a pointer would --
  //! measured safe: no corpus driver mutates a model option after setup().
  //!
  //! What the MODEL layer can name is now only its own 12 fields: a solver option referenced inside
  //! ffmodel.hpp fails to compile instead of passing silently, which is the point of the split.
  struct Options
  : public FFModel::Options
  {
    //! @brief the interface-condition type.  ONE definition, OCPlan's: ocplan.hpp sits below this
    //! header in the include order, and it had carried an identical copy whose comment said it "must stay in
    //! sync" with FFModel's.  Aliased here, with its enumerators re-exported, so Options::IC_VALUE and
    //! OCFESLV::Options::IC_AUTO resolve exactly as before.
    using InterfaceType = OCPlan::InterfaceType;
    static constexpr InterfaceType IC_VALUE  = OCPlan::IC_VALUE;
    static constexpr InterfaceType IC_FLUX   = OCPlan::IC_FLUX;
    static constexpr InterfaceType IC_UPWIND = OCPlan::IC_UPWIND;
    static constexpr InterfaceType IC_AUTO   = OCPlan::IC_AUTO;

    // rev313: the solver's enums, moved verbatim from FFModel::Options -- they describe how this solver
    // discretises, imposes, factorises and caches, which a model has no business knowing.
    //! the imposition type.  ONE definition, OCPlan's (ocplan.hpp sits below this header), aliased here with
    //! its values re-exported -- as InterfaceType is -- so no static_assert has to keep two copies in step.
    using ImpositionType = OCPlan::ImpositionSentinel;
    static constexpr ImpositionType IC_WEAK   = OCPlan::IMPOSITION_IC_WEAK;
    static constexpr ImpositionType IC_STRONG = OCPlan::IMPOSITION_IC_STRONG;
    static constexpr ImpositionType IC_TRACE  = OCPlan::IMPOSITION_IC_TRACE;
    //! @brief What setup() does with interface continuity claims the structural analysis finds redundant -- the
    //! values of INTERFACE.DROP_POLICY: drop and verify after the solve, re-derive what verification flags, or keep all.
    enum InterfaceDropPolicy
    {
      DROP_VERIFY   = 0,  //!< structural drop + post-solve verify gate (the default until 2026-10-14)
      DROP_REDERIVE = 1,  //!< verify-driven re-derive: protect over-dropped claims, re-solve
      DROP_OFF      = 2   //!< disable the drop; keep every continuity claim explicit (DEFAULT since 2026-10-14)
    };
    enum OutputGridPolicy
    {
      OUT_ELEMENT_LOCAL   = 0,
      OUT_UNIQUE_PHYSICAL = 1
    };
    //! @brief Rank backend of the determinacy audit (2026-10-09: one value per line, so that the Python binder registers
    //! them -- the one-line form left the enum empty in Python).
    enum DETERMINACY_BACKEND_T
    {
      DET_AUTO  = 0,  //!< SuiteSparseQR if built, else Eigen if built, else dense under the cap (default)
      DET_SPQR  = 1,  //!< SuiteSparseQR only; refused by setup() in a build without it (BACKEND_UNAVAILABLE)
      DET_EIGEN = 2,  //!< Eigen's sparse QR only; refused by setup() in a build without Eigen (BACKEND_UNAVAILABLE)
      DET_DENSE = 3   //!< dense SVD only, under the size cap
    };

    //! @brief How continuity between elements is expressed and imposed: the condition type, the imposition (penalty, trace multiplier, strong), SAT penalties and trace-multiplier scaling.
    struct t_Interface
    {
      //! @brief Imposition strategy for interface conditions: IC_STRONG (default since 2026-10-06; exact continuity rows),
      //! IC_WEAK (penalty terms of strength SAT_SIGMA0) or IC_TRACE (trace multipliers)
      ImpositionType IMPOSITION;
      //! @brief Interface condition type for multi-element continuity
      InterfaceType TYPE;
      //! @brief SAT penalty coefficient for C0 (value) continuity.
      //! Effective penalty = INTERFACE.SAT_SIGMA0 * tau_d  per direction d.  Default 10.0 (since 2026-10-03;
      //! 1.0 was below threshold at moderate refinement -- a ~10x worse answer, silently -- and 100 can stall the
      //! monolithic Newton solve).
      //! Default penalty parameter for direction d:
      //!   tau_d = (n_node - 1)^2 / w_elem    (spectral scaling)
      //! scaled by INTERFACE.SAT_SIGMA0 (C0 continuity) or CRONOS_SAT_SIGMA1 (C1).
      //! When classify_pde() has been called, the penalty is coupled
      //! through the principal symbol: sum_j A_d[k,j] * sigma * tau_d * Delta u_j.
      double SAT_SIGMA0;
      //! @brief Policy for the IC_TRACE/IC_STRONG redundant-continuity drop
      //! (item 10).  DROP_VERIFY keeps the structural receiver-B drop and treats a
      //! post-solve over-drop (verify_interface_drop) as a hard failure.
      //! DROP_REDERIVE additionally lets the driver protect the over-dropped
      //! claims and re-derive/re-solve until the verify gate is clean (correct on
      //! the linear suite, and the safe path for the nonlinear case).  DROP_OFF
      //! disables the drop entirely (every continuity claim kept explicit) -- the
      //! most conservative, always-correct baseline, and the DEFAULT since 2026-10-14 (4R.28).
      //! MEASURED, full OCFE corpus with CRONOS_DROP_POLICY=2: 182/182 programs, PASS and FAIL
      //! identical (3646 / 17), classification identical, 23.0 min against 23.2 -- the drop changes
      //! nothing the corpus asserts and costs nothing to disable.  Against that, where it DOES fire
      //! it removed claims that were not implied: OCFE_Burger_move --form q (square, PARABOLIC, hence
      //! drop-eligible) left the mesh flux g discontinuous by 14 across seams and failed after 60
      //! iterations, where keeping the claims converges in 7 with every state continuous to round-off.
      //! The eligibility rule's parabolic branch argues that flux continuity follows from value
      //! continuity plus a constitutive LINK; a mesh relation g = M(q) x_xi is not constitutive in
      //! that sense and the implication fails.  CRONOS_DROP_POLICY=0 restores the drop.
      InterfaceDropPolicy DROP_POLICY;
      //! @brief Column scaling for generated IC_TRACE tau variables.
      //! Claim-centred multipliers enter receiver rows as INTERFACE.TRACE_TAU_SCALE*lambda.
      //! The projected u-solution is unchanged, but a non-unit value can improve
      //! conditioning when tau variables are much larger/smaller than state rows.
      double TRACE_TAU_SCALE;
      size_t TRACE_PROJ_MAX;       //!< Largest n_trace_var for which the trace projector is computed;
                                   //!< 0 disables the guard.  Default 6000.
      double TRACE_PROJ_EPS;       //!< Multiplier-block regularisation for the trace projection.
                                   //!< The k continuity conditions the receiver graph cannot support
                                   //!< are enforced through a spring of stiffness 1/eps rather than
                                   //!< dropped.  Default 1.0.
    } INTERFACE;

    //! @brief What the solver retains and emits after a solve.
    struct t_Output
    {
      bool   MARCH_STORE; //!< SOLVE.MARCHING: keep each element's block solution (O(N)) for plotting/inspection via march_trajectory()/eval_solution() (default true)
      //! @brief Grid policy for distributed post-evaluation outputs only.
      //! OUT_ELEMENT_LOCAL keeps one output row per element-local collocation
      //! coefficient, including duplicate physical interface points for LGL/CGL.
      //! OUT_UNIQUE_PHYSICAL skips the duplicate left trace of every non-first
      //! LGL/CGL element in distributed output dimensions. Equations, states,
      //! node_colloc(), interface/SAT handling and scalar point outputs remain
      //! element-local.
      OutputGridPolicy GRID_POLICY;
    } OUTPUT;

    //! @brief Reports that diagnose the user's MODEL.  They change what is printed, never the result.
    struct t_Audit
    {
      bool   EQN_LEGEND;     //!< print the eN -> equation legend (roles, bare states, partial orders per
                                   //!< direction) at the end of setup() (false).  Printed unconditionally at the
                                   //!< head of a spectrum audit regardless of this flag.  The legend is what makes
                                   //!< an audit's "by eqn eN" resolvable against the model, and its `partials`
                                   //!< column identifies the assembled formulation independently of any banner.
                                   //!< Environment default: CRONOS_AUDIT_LEGEND, else CRONOS_AUDIT.
      bool   COVERAGE;       //!< print the structural column-coverage report and the claim census at the end
                                   //!< of setup() (false).  Setup-time only: no solve, no SVD.  A state column
                                   //!< covered by NO row is a free degree of freedom; singly-covered columns are
                                   //!< the fragility margin.  Necessary, not sufficient -- see the method doc.
                                   //!< Environment default: CRONOS_AUDIT_COVERAGE, else CRONOS_AUDIT.
      bool   SPECTRUM;       //!< run the singular-value / left+right null-space audit at iteration 0 of every
                                   //!< solve(), not only when tier-0 Newton is found singular (false).  Needed to
                                   //!< compare a CONVERGING run against a failing one: the historical trigger fired
                                   //!< only on failure, so a wrong-but-converged solve was never audited.
                                   //!< Environment default: CRONOS_AUDIT_SPECTRUM (NOT implied by CRONOS_AUDIT).
      size_t SPECTRUM_MAX;         //!< largest dense-SVD the spectrum audit will attempt,
                                   //!< in either dimension.  Was a hard literal 3000, which put
                                   //!< `rel` out of reach on exactly the configurations where the
                                   //!< question lives (MMPDE13 at production order reaches the
                                   //!< limit by NEL_T=4).  Cost is O(n^3): 3000 is seconds, 6000
                                   //!< a few minutes.  Environment: CRONOS_SPECTRUM_MAX.
      bool   RANK_SYMMETRY;  //!< print the saddle-point rank-symmetry report at the end of setup():
                                   //!< rank of the multiplier column block B against rank of the continuity
                                   //!< row block C.  rank(B) < n_trace_var makes the augmented system singular
                                   //!< for ANY C and cannot be repaired by dropping or pinning rows.
                                   //!< Environment default: CRONOS_AUDIT_RANK, else CRONOS_AUDIT.
    } AUDIT;

    //! @brief The numerical-rank determinacy audit run at the end of setup().
    struct t_Determinacy
    {
      // rev133: setup-time determinacy audit (numerical rank of the assembled Jacobian).
      bool   AUDIT;    //!< run the numerical-rank determinacy audit at the end of setup() (false).
                                   //!< One report line per audited window, unconditional on DISPLAY_LEVEL.
                                   //!< Environment: CRONOS_AUDIT_DETERMINACY.
      double TOL;      //!< the SINGLE relative rank tolerance, applied identically to rank(M) and
                                   //!< rank([B;D]) (1e-12; the value the spectrum audit's rank_deficiency[]
                                   //!< uses).
      DETERMINACY_BACKEND_T BACKEND; //!< rank backend: DET_AUTO = SPQR if built, else Eigen if
                                   //!< built, else dense under the cap.  Values:
                                   //!< AUTO|SPQR|EIGEN|DENSE.  An explicit DET_SPQR / DET_EIGEN in a
                                   //!< build without that backend is REFUSED by setup()
                                   //!< (BACKEND_UNAVAILABLE) -- it used to skip the check silently.
    } DETERMINACY;

    //! @brief Which conditions make setup() FAIL rather than warn.
    struct t_Fatal
    {
      bool   DETERMINACY;    //!< setup() returns false (SetupStatus::DETERMINACY_UNDETERMINED) on a
                                   //!< STATE CONTENT verdict (false).
      bool   INTERFACE_RANK; //!< setup() returns false (INTERFACE_PLAN_INVALID) when the multiplier block
                                   //!< B is column rank-deficient.  Default FALSE: the resulting kernel is
                                   //!< singular but not always fatal (a model whose residual satisfies the
                                   //!< compatibility conditions still converges), so this stays opt-in until
                                   //!< the corpus has been re-baselined.
      bool   EVOLUTION_DOMAIN;  //!< setup() returns false on an unsuitable evolution domain or an IC that cannot be set consistently (true); false (default) downgrades to a DISPLAY_LEVEL warning and proceeds
      bool   REDUCED_DOF; //!< setup() returns false on a structurally IC/BC-inconsistent reduced system (true); false downgrades to a warning and proceeds (rank-revealing solve)
    } FATAL;


    Options
      ()
      { reset(); }

    Options& operator=
      ( Options const& opt )
      {
        FFModel::Options::operator=( opt );

        INTERFACE.TYPE      = opt.INTERFACE.TYPE;
        INTERFACE.IMPOSITION     = opt.INTERFACE.IMPOSITION;
        INTERFACE.DROP_POLICY = opt.INTERFACE.DROP_POLICY;
        INTERFACE.SAT_SIGMA0          = opt.INTERFACE.SAT_SIGMA0;
        INTERFACE.TRACE_TAU_SCALE     = opt.INTERFACE.TRACE_TAU_SCALE;
        OUTPUT.GRID_POLICY  = opt.OUTPUT.GRID_POLICY;
        AUDIT.EQN_LEGEND    = opt.AUDIT.EQN_LEGEND;
        AUDIT.COVERAGE      = opt.AUDIT.COVERAGE;
        AUDIT.SPECTRUM      = opt.AUDIT.SPECTRUM;
        AUDIT.RANK_SYMMETRY = opt.AUDIT.RANK_SYMMETRY;
        DETERMINACY.AUDIT   = opt.DETERMINACY.AUDIT;
        FATAL.DETERMINACY   = opt.FATAL.DETERMINACY;
        DETERMINACY.TOL     = opt.DETERMINACY.TOL;
        DETERMINACY.BACKEND = opt.DETERMINACY.BACKEND;
        FATAL.INTERFACE_RANK = opt.FATAL.INTERFACE_RANK;
        INTERFACE.TRACE_PROJ_EPS = opt.INTERFACE.TRACE_PROJ_EPS;
        INTERFACE.TRACE_PROJ_MAX = opt.INTERFACE.TRACE_PROJ_MAX;
        AUDIT.SPECTRUM_MAX     = opt.AUDIT.SPECTRUM_MAX;
        FATAL.REDUCED_DOF = opt.FATAL.REDUCED_DOF;
        OUTPUT.MARCH_STORE = opt.OUTPUT.MARCH_STORE;
        FATAL.EVOLUTION_DOMAIN  = opt.FATAL.EVOLUTION_DOMAIN;
        return *this;
      }

    void apply_environment_defaults
      ()
      {
        FFModel::Options::apply_environment_defaults();

          AUDIT.EQN_LEGEND = _env_flag( "CRONOS_AUDIT_LEGEND", AUDIT.EQN_LEGEND );
          AUDIT.COVERAGE   = _env_flag( "CRONOS_AUDIT_COVERAGE", AUDIT.COVERAGE );
          AUDIT.SPECTRUM   = _env_flag( "CRONOS_AUDIT_SPECTRUM", AUDIT.SPECTRUM );
          AUDIT.RANK_SYMMETRY = _env_flag( "CRONOS_AUDIT_RANK", AUDIT.RANK_SYMMETRY );
          DETERMINACY.AUDIT     = _env_flag( "CRONOS_AUDIT_DETERMINACY", DETERMINACY.AUDIT );
      }

    void reset
      ()
      {
        FFModel::Options::reset();

        INTERFACE.TYPE      = IC_AUTO;
        INTERFACE.IMPOSITION     = IC_STRONG;  // since 2026-10-06 (formerly IC_WEAK): exact continuity rows
        INTERFACE.DROP_POLICY = DROP_OFF;      // since 2026-10-14 (was DROP_VERIFY); see the field doc
        INTERFACE.SAT_SIGMA0          = 10.0;
        INTERFACE.TRACE_TAU_SCALE     = 1.0;
        OUTPUT.GRID_POLICY  = OUT_ELEMENT_LOCAL;
        {
          // CRONOS_AUDIT enables the two CHEAP setup-time reports (legend + coverage).
          // CRONOS_AUDIT_SPECTRUM is deliberately excluded: it runs a dense SVD.
          bool const _aud  = _env_flag( "CRONOS_AUDIT", false );
          AUDIT.EQN_LEGEND =   _aud;
          AUDIT.COVERAGE   = _aud;
          AUDIT.SPECTRUM   = false;
          // rev7: saddle-point rank symmetry.  Two sparse Gram-Schmidt passes over
          // blocks with a handful of nonzeros per column/row -- cheap enough to ride
          // along with CRONOS_AUDIT, unlike the dense SVD.
          AUDIT.RANK_SYMMETRY = _aud;
          // rev11: structural corpus sweeps need setup() only.  Aborting after the audits
          // avoids running 118 solves (and any IPOPT work) to obtain numbers that are all
          // fixed at setup time.  Off unless explicitly requested.
          // rev133: determinacy audit defaults (report-only, ON).
          // rev272 (C4 CLOSED): DEFAULT FALSE.  The structural rule reproduces this audit exactly --
          //   k_pred == the audited k' on 455 of 455 exact-mode pairs (rev224 sweep 20260915_115625),
          //   IC_STRONG 373/373 and IC_TRACE 82/82, zero under-reports, zero over-reports;
          //   an audit-off corpus sweep (20260915_153646) reported the SAME 48 deficiencies structurally,
          //   151/151 completed, and wall time fell 2,517 -> 1,815 s (28%).
          // The record is filled structurally when the audit does not run (rev225), so determinacy_audit()
          // callers still get a truthful answer with backend="STRUCTURAL".
          // The one blocker ever raised against this -- that OCFE_PDE6_blk0 consumes the audit -- was MY
          // MISREADING: its det_na counts a FAILED RUN, not a skipped audit; cmp() compares solution errors.
          // Verified directly: blk0 PASSES with CRONOS_AUDIT_DETERMINACY=0 (24 structural lines, 0 audited).
          // CRONOS_AUDIT_DETERMINACY=1 restores the SPQR audit as the cross-check that must agree.
          DETERMINACY.AUDIT     = false;
          FATAL.DETERMINACY     = false;
          DETERMINACY.TOL       = 1e-12;
          DETERMINACY.BACKEND   = DET_AUTO;
        }
        // rev42: automatic per-state duplicate-node spreads after every solve().
        //
        // DELIBERATELY OUTSIDE the CRONOS_AUDIT umbrella and defaulted OFF: it changes the
        // output of every driver that enables it, and a corpus of 135 drivers with golden
        // logs should not have its output move because an unrelated audit was switched on.
        //
        // WHY THE HEADER CALLS THIS ITSELF.  Per-state spreads are the only instrument that
        // caught the withdrawn Stage A rule -- which passed k / SYMMETRIC / residual / PASS
        // while destroying MBC5's radial fluxes by 24 orders.  Implemented driver-side, it
        // reached 3 of the 29 drivers the rank sweep flags, and 28 separate driver edits
        // would leave three independent reimplementations of one measurement plus a fresh
        // gap every time a driver is added.  One env flag instruments the whole corpus and
        // cannot be forgotten.
        AUDIT.SPECTRUM_MAX     = 3000;
        if( char const* v = std::getenv( "CRONOS_SPECTRUM_MAX" ) ){
          long const n = std::atol( v );
          if( n > 0 ) AUDIT.SPECTRUM_MAX = (size_t)n;
        }
        FATAL.INTERFACE_RANK = false;
        // Stage A -- WITHDRAWN, DEFAULT OFF.  Retained only as a reproducible negative
        // result; set CRONOS_LINK_PAIR_SUPPRESS=1 to re-enable it.
        //
        // The rule withdrew a derivative auxiliary's C0 claim where its defining LINK row is
        // active and the primitive is itself claimed.  It reached k=0 -- PDE2 482->466
        // multipliers, MBC5 175->172, both SYMMETRIC -- and destroyed the property exact
        // imposition exists to provide.  MBC5fig_ebal at PRODUCTION mesh, IC_STRONG, per-state
        // duplicate-node spread, baseline vs rule (vert; horiz identical in character):
        //     Drl_Cl   2.69e-26 -> 3.02e-08
        //     Drl_Ml   1.21e-27 -> 2.09e-03
        //     Drl_Pl   1.18e-26 -> 9.39e-03
        //     15 other states       unchanged at round-off
        // Perfectly localised to the three withdrawn claims.  Drl_* are radial fluxes in the
        // graded reaction film, so this is a physical error in the absorption model, not a
        // bookkeeping one.  PDE2 IC_STRONG shows the same signature (3.469e-18 -> 2.133e-01,
        // FAIL).  PDE2 IC_TRACE moved 1.388e-17 -> 1.216e-10, which is NOT itself harmful --
        // it sits at the level of everything else in that run and two orders under the gate --
        // but it does mean IC_TRACE stopped beating its own IC_WEAK baseline (4.31e-11) while
        // still paying for 482 extra unknowns and a QR solve.
        //
        // In all cases the primitive-variable errors, the residual, the rank symmetry and the
        // convergence were CLEAN.  Only the auxiliary spreads moved.  A k/SYMMETRIC/PASS
        // triple would have read as success.
        //
        // Three further defects, all measured:
        //   * OVER-BROAD.  Fires 1704x on PSA9, which has k=0 and no deficiency at all, and on
        //     7 state ids in MBC5 where only 3 claims are deficient.  The (primitive,
        //     derivative-auxiliary, LINK row) structure is NECESSARY but NOT SUFFICIENT for the
        //     multiplier deficiency, so no purely structural emission-time rule can separate
        //     the deficient pairs from the healthy ones.
        //   * UNRELIABLE.  Suppression is per node-visit, not per claim key: PSA9 withdrew
        //     1704 times with ZERO net change to n_trace_var, the claims being re-emitted from
        //     another visit.
        //   * INAPPLICABLE where it was wanted.  MMPDE9 defines gg = g'xi as an INTERIOR
        //     equation, not a reduce_order LINK row, so the predicate never fires on the model
        //     that motivated the arc.
        //
        // ROOT CAUSE, and why Option D is NOT affected.  The compatibility test u_k^T r ~ 0
        // says a COMBINATION of continuity rows lies in the row space of the physical block --
        // the left-null direction.  It does NOT say any individual participating claim is
        // redundant.  Participation was measured at 0.5/0.5, so the implied object is roughly
        // [a] - alpha*[Da], not [Da] alone.  Deleting a claim ROW removes a genuine constraint;
        // the projector C_tilde=(I-P)C keeps every row and removes only the component along the
        // implied combination.  Design rev4 2.2 conflated the two and needs rewriting.
        //
        // Direct evidence for the distinction: MBC5's SHIPPING configuration is ASYMMETRIC,
        // k=3, three orphan multiplier columns -- and its Drl_* continuity is exact to 1e-26.
        // The orphan COLUMNS are a harmless gauge freedom; the continuity ROWS are load-bearing.
        // Stiffness knob for the k unrealisable continuity conditions: they are enforced
        // through a spring of stiffness 1/INTERFACE.TRACE_PROJ_EPS.  eps -> 0 recovers the singular
        // unprojected system; eps -> infinity recovers dropping them, which diverged on
        // MMPDE9.  1.0 is a starting point, not a tuned value.
        INTERFACE.TRACE_PROJ_EPS = 1.0;
        // Largest n_trace_var the projector will be computed for.  0 = no limit.  6000 covers
        // PDE5/PDE9 (m=3272) and MBC at production mesh (m=2116) with headroom.
        INTERFACE.TRACE_PROJ_MAX = 6000;
        FATAL.REDUCED_DOF = true; //!< reject setup() on a structurally IC/BC-inconsistent reduced system
        FATAL.EVOLUTION_DOMAIN  = false; //!< warn-by-default (v1): true rejects setup() on an unsuitable evolution domain / inconsistent IC
        OUTPUT.MARCH_STORE = true;   //!< default ON: keep each window block so eval_solution()/eval_colloc() can read the full marched field (Stage A of read consolidation)
        apply_environment_defaults();
      }
  };

  //! @brief The solver's options.  Shadows FFModel::options deliberately: this is what the caller sets,
  //! FFModel::options is what setup() applies.
  //! @brief The options the base reads, as a REFERENCE to the derived solver's own `Options options` (each
  //! solver's Options derives from this one): one object, two views.  Bound in the constructor; never read there.
  Options& options;

protected:
  //! @brief vector of equation subgraphs
  mutable std::vector<FFSubgraph>    _sgEqn;

  //! @brief Count of collocated equation residuals during setup
  size_t                             _nCollEqn = 0;

  //! @brief Count of collocated states and generated trace/tau unknowns during setup (constants in _vCst excluded)
  size_t                             _nCollVar = 0;

  //! @brief Count of collocated inputs during setup (constants in _vCst excluded)
  size_t                             _nCollInp = 0;

public:
  //! @brief The equations carrying an evolution derivative (filled when the rows are counted): the EVOLUTION rows,
  //! whose slot at an element start the continuity row takes.  An `ALL - LB` row WITHOUT one is algebraic (its
  //! mask only reflects the t = 0 convention, e.g. x = a(t) in OCFE_PDE20e) and stays collocated at element starts.
  //! Also excluded: BOUNDARY rows -- a differential state's closure at a boundary node, standing in for its
  //! evolution row there (the heated end of OCFE_CAUSAL: D(u) u_z = q with the continuity row for u at the same
  //! node would over-determine it).  An INTERIOR algebraic row (x = a(t)) holds at every node and stays.
  //! Keyed by the equation's INDEX in _mEqn (a deep copy preserves it; the equations' DAG ids it does not).
  mutable std::set<size_t> _evoDerivRows;

  //! @brief The stream for an informational message of display level @p n (DISPLAY_REVIEW_20261001): std::cout when
  //! options.DISPLAY_LEVEL >= n, a null stream otherwise.  Refusals and failures go to std::cerr at every level.
  std::ostream& _disp
    ( int const n )
    const
    { static std::ostream null( nullptr ); return options.DISPLAY_LEVEL >= n? std::cout: null; }

protected:
  //! @brief Verify that public interface options still match the frozen plan.
  void _assert_frozen_interface_options_current() const;

  //! @brief Select first/next block given current block ndx_el corresponding to the variable domain var_dom
  bool _next_block
    ( std::map<FFVar, size_t, lt_FFVar>& ndx_el, std::set< FFVar, lt_FFVar > const& var_dom )
    const;

  //! @brief Select first/next block given current block ndx_el corresponding to the equation domain eqn_dom
  bool _next_block
    ( std::map<FFVar, size_t, lt_FFVar>& ndx_el, std::map< FFVar, int, lt_FFVar > const& eqn_dom )
    const;

  //! @brief Locate the finite-element block and reference coordinate for a physical point
  std::pair<size_t,double> _point_in_domain
    ( FFDom const& dom, double const& val, int const side = -1 )   // -1: legacy (internal callers); else FFDom::SIDE
    const;

  //! @brief Generate collocation nodes for domain dom in selected element
  std::vector<double> _node_block
    ( FFDom const& dom, size_t const iel )
    const;

  //! @brief Generate collocation nodes for variable var in selected block
  std::vector<double> _node_block
    ( FFVar const& var, std::map<FFVar, size_t, lt_FFVar> const& ndx_el )
    const;

  //! @brief Effective discretisation of input inp along physical domain dom
  //! @brief The declared model has been imported onto a new working DAG (for OCFESLV and OCFEMOL alike): derive the working input collocation
  //! from FFModel's declaration (the input shares its domains' grid and differs only in local order and node
  //! type).  @p vAllVarsLoc / @p ndxAllVars are the positional index the import produced.
  void _on_model_imported
    ( std::vector<FFVar> const& vAllVarsLoc,
      std::map<FFVar::pt_idVar,size_t> const& ndxAllVars )
    override
    {
    // The input shares its domains' grid: take each domain's FFDom and apply the input's declared
    // collocation.  (Through rev192 the FFDom was captured when add_input() ran; deriving it here means a
    // domain REDECLARED after its inputs now governs them too -- NOTES_20260911m Sec.3.)
    _mInpDisc.clear();
    for( auto const& [inp, decl] : _usr._mInpDiscUsr ){
      t_Dom disc;
      for( auto const& [dvar, ic] : decl ){
        auto itdom = _usr._mDomUsr.find( dvar );
        if( itdom == _usr._mDomUsr.end() )
          throw Exceptions( Exceptions::INDEX );
        FFDom opt = itdom->second;
        opt.type = ic.type;
        if( ic.n_node ) opt.n_node = ic.n_node;
        if( opt.n_node == 0 )
          throw FFDom::Exceptions( FFDom::Exceptions::NODES );
        disc[ _local( dvar, vAllVarsLoc, ndxAllVars ) ] = opt;
      }
      _mInpDisc[ _local( inp, vAllVarsLoc, ndxAllVars ) ] = std::move( disc );
    }
    }

  FFDom const& _input_dom
    ( FFVar const& inp, FFVar const& dom )
    const;

  //! @brief Build the OCVar domain map for an input variable
  template <typename U>
  typename OCVar<U>::t_Dom _input_ocdom
    ( FFVar const& inp, std::set<FFVar,lt_FFVar> const& dep )
    const;

  //! @brief Templated initialization of collocation variable vectors
  template <typename U>
  void _reset_colloc
    ( std::vector<OCVar<U>>& OCvar, std::vector<OCVar<U>>& OCinp,
      std::vector<OCVar<U>>& OCdom, std::vector<OCVar<U>>& OCcst )
    const;

  //! @brief Templated evaluation of collocated equation residuals
  template <typename U>
  void _eval_eqn
    ( U* eqn, U const* var, U const* inp, 
      OCVar<U>& OCres, std::vector<OCVar<U>>& OCvar,  std::vector<OCVar<U>>& OCinp,
      std::vector<OCVar<U>>& OCdom, std::vector<OCVar<U>>& OCcst, std::vector<OCVar<U>>& OCwrk )
    const;

  //! @brief Templated evaluation of the collocated output functions (the base's, for OCFESLV and OCFEMOL alike; the
  //! evolution direction through _evolution_window and _window_output_rows)
  template <typename U>
  void _eval_fct
    ( U* fct, U const* var, U const* inp,
      OCVar<U>& OCres, std::vector<OCVar<U>>& OCvar,  std::vector<OCVar<U>>& OCinp,
      std::vector<OCVar<U>>& OCdom, std::vector<OCVar<U>>& OCcst, std::vector<OCVar<U>>& OCwrk )
    const;

public:
  //! @brief Get flat index to start of block ndx_el for input var in data array
  size_t pos_input
    ( FFVar const& var, std::map<FFVar, size_t, lt_FFVar> const& ndx_el )
    const;

  //! @brief Get flat index to start of block ndx_el for state var in data array
  size_t pos_state
    ( FFVar const& var, std::map<FFVar, size_t, lt_FFVar> const& ndx_el )
    const;

  //! @brief True if @p input is a registered control.  Matched by VARIABLE INDEX (id().second), not
  //! by full-key map lookup: _controls is keyed by the caller's original FFVar (register_control is
  //! called with it), whereas marching-internal structures (_classInpRefFun, _marchTransfers) hold
  //! post-setup DAG copies, so _controls.count() would compare mismatched DAG pointers and miss.  This
  //! mirrors the id().second convention already used for marching IC-transfer matching.
  bool _is_registered_control
    ( FFVar const& input ) const
    { for( auto const& cc : _controls ) if( cc.first.id().second == input.id().second ) return true;
      return false; }

  //! @brief Write @p values into the collocated DOFs of distributed input @p input
  //! within the input-value array @p inp, in the SAME element-major / node-inner
  //! ordering that terminal_profile() produces -- so a transferred terminal profile
  //! can be re-seeded directly onto a marching IC (or forcing) input between element
  //! solves, with no re-setup().  Returns false on unknown input or a size mismatch.
  bool set_input_values
    ( FFVar const& input, std::vector<double> const& values, double* inp ) const;

  //! @brief Set distributed input @p input (in inp[]) from a callable @p fun
  //! (a t_Coord -> double) evaluated at the input's own collocation nodes.  General
  //! per-element seeding of marching IC / boundary / forcing inputs from a physical
  //! (space and/or time) function.  Templated on the callable to avoid type erasure.
  template <typename FUN>
  bool sample_input
    ( FFVar const& input, FUN&& fun, double* inp ) const;

protected:
  //! @brief Initialize collocated state and input coefficient arrays from references.
  void _init_reference
    ( double* var, double* inp, double const* cst ) const;

  //! @brief Evaluate an auxiliary definition into the state reference array.
  void _init_auxiliary_reference
    ( t_AuxDef const& aux, double* var, double const* inp, double const* cst ) const;

  //! @brief Sample a collocated state/input reference at the scalar classification point.
  double _sample_reference
    ( FFVar const& var, std::set<FFVar,lt_FFVar> const& dom,
      std::tuple<size_t,size_t,std::vector<size_t>> const& blk, double const* data ) const;

  //! @brief Reference values for the classification, sampled from the COLLOCATED field at each domain's
  //! reference coordinate (overrides FFModel::_classification_reference, which evaluates the model directly).
  //! CRONOS_MODEL_REF=1 defers to the model-level reference instead -- identical for scalar references, and
  //! differing by the interpolation error where a reference FUNCTION is declared.
  void _classification_reference
    ( std::vector<double>& state_ref, std::vector<double>& input_ref,
      std::vector<double>& cst_ref,   std::vector<double>& dom_ref )
    override;

  //! @brief Sample a distributed-input reference function on the element's INTERIOR side at
  //! interior-interface boundary nodes, so a DISCONTINUOUS (jumping) input gives each element its
  //! own one-sided value rather than inheriting its neighbour's across the shared face coordinate.
  //! At a boundary node (first node of a non-first element, or last node of a non-last element in a
  //! direction) the face coordinate is nudged a tiny fraction of the element width toward the
  //! interior and the ref function re-evaluated; the interior value is used ONLY if it differs from
  //! the exact-face value beyond a small relative tolerance (a genuine jump).  A continuous input
  //! (both sides equal) keeps its exact-face value, so it -- and every input left at the
  //! DISCONTINUOUS default -- is byte-for-byte unchanged.  @p dirs/@p elem_idx/@p node_idx/@p
  //! node_cnt are aligned per-direction (dirs in the input's domain-set order).
  template <typename FUN>
  double _input_ref_interior_side
    ( FUN const& fun, t_Coord const& exact, double exact_val, FFVar const& inpvar,
      std::vector<FFVar> const& dirs, std::vector<size_t> const& elem_idx,
      std::vector<size_t> const& node_idx, std::vector<size_t> const& node_cnt ) const
    {
      t_Coord icoord = exact;
      bool perturbed = false;
      for( size_t j = 0; j < dirs.size(); ++j ){
        FFDom const& dd = _input_dom( inpvar, dirs[j] );
        size_t const e = elem_idx[j], ni = node_idx[j], nn = node_cnt[j];
        double const w = dd.elem_width( e ), eps = 1e-7 * ( w > 0. ? w : 1. );
        // Only re-sample a boundary node that ACTUALLY sits on the interior interface.
        // For LGL the first/last nodes coincide with elem_lo/elem_up, so the perturbation
        // recovers the element's one-sided value at a jump.  For LGR (left-Radau) the LB is
        // at elem_lo but the LAST node is an INTERIOR right-Radau point (e.g. tau~0.845), NOT
        // at elem_up -- perturbing it to elem_up-eps would move the sample off its own node
        // onto the interface, corrupting the per-element input polynomial.  Gate on the node
        // coordinate so LGR's interior last node is left at its exact position.
        double const pos = icoord[ dirs[j] ];
        bool const at_lo = std::fabs( pos - dd.elem_lo( e ) ) < eps;
        bool const at_up = std::fabs( pos - dd.elem_up( e ) ) < eps;
        if( ni == 0 && e > 0 && at_lo )
          { icoord[ dirs[j] ] = dd.elem_lo( e ) + eps; perturbed = true; }        // LB of interior element (at the interface)
        else if( nn > 0 && ni + 1 == nn && e + 1 < dd.n_elem && at_up )
          { icoord[ dirs[j] ] = dd.elem_up( e ) - eps; perturbed = true; }        // UB of interior element (at the interface; LGL only)
      }
      if( !perturbed ) return exact_val;
      double const vi = fun( icoord );
      return ( std::fabs( vi - exact_val ) > 1e-4 * ( 1.0 + std::fabs( exact_val ) ) )
               ? vi : exact_val;                                                  // interior side only on a genuine jump
    }

public:
  static bool _taucol_enabled()
    {
      static bool const k = Options::_env_flag( "CRONOS_AUDIT_TAUCOL", false );
      return k;
    }

  //! @brief effective INTERFACE.SAT_SIGMA0 / CRONOS_SAT_SIGMA1, with an environment override.
  //!
  //! WHY AN OVERRIDE AT ALL.  INTERFACE.SAT_SIGMA0 = 1.0 is BELOW THRESHOLD at moderate refinement:
  //! MMPDE10 at --nelt 6 gives errNode 4.19e-01 at sigma0 = 1 and 8.44e-03 at sigma0 >= 10,
  //! saturating (10/100/1000 identical).  A driver below threshold returns a ~10x worse
  //! answer SILENTLY -- there is no convergence failure and no diagnostic.  The corpus has
  //! 29 drivers still at 1.0 against 43 already raised to 10.0, so the population is already
  //! split and nobody recorded why.  The audit needs to sweep sigma0 across the whole corpus
  //! WITHOUT editing 29 drivers, each edit being a chance to introduce a difference.
  //!
  //! WHY IT IS READ HERE AND NOT IN Options::reset().  Drivers assign INTERFACE.SAT_SIGMA0 in their
  //! body, AFTER construction.  An override applied at reset() would be overwritten by every
  //! driver that sets it -- i.e. by exactly the 72 drivers the audit is about.  Reading it at
  //! the single consumption site is the only placement that actually overrides.
  //!
  //! NOTE ON WHAT sigma0 IS.  It is a DIMENSIONLESS MULTIPLIER on tau_d = (n_node-1)^2/w_elem,
  //! which is already the p^2/h spectral scaling.  So the roadmap's "a default that scales
  //! with p^2/h rather than a fixed constant" is ALREADY IMPLEMENTED -- what is unknown is the
  //! theory CONSTANT in front of it, and the evidence says the current 1.0 is roughly 10x too
  //! small.  That is a much narrower question than the roadmap states, and it is what the
  //! calibration run is for.
  //!
  //! Since 2026-10-07 the option alone decides: the CRONOS_SAT_SIGMA0 override, which replaced even a value the driver
  //! set, is retired (WORKPLAN 3.A).
  double _eff_sat_sigma0() const
    {
      return options.INTERFACE.SAT_SIGMA0;   // the CRONOS_SAT_SIGMA0 override (B1 audit) retired 2026-10-07
    }

  //! The time discretisation of the monolithic solve is CAUSAL whenever an evolution domain is set (2026-10-03: the
  //! option INTERFACE.EVOLUTION_CAUSAL and CRONOS_EVOLUTION_CAUSAL are retired): every element's first evolution node
  //! is an initial node, as a marching window's start -- the evolution and boundary rows are not collocated there,
  //! and each state with an evolution derivative gets an exact continuity row to the previous element (x(tau+) = w
  //! at a transition).  The evolution faces carry no claims.
  //! @brief Whether the evolution direction is COLLOCATED: true for OCFESLV, where t is one collocated direction
  //! among the others and the causal time discretisation below applies; false for OCFEMOL (MoL phase 3), where
  //! t is the integrator's and every branch that treats it as a collocated domain -- the element-start
  //! continuity rows, the evolution-derivative row census of _check_consistency, the evolution-input jump
  //! probe -- does not exist.  A flag rather than a test on _mDom: _check_consistency also runs on the user
  //! model before _mDom is filled, and the answer must not depend on when it is asked.
  bool _evoCollocated = true;

  bool _causal_rows_active() const
    { return _evolution_dom_set && _evolution_dom_var.dag() && _evoCollocated; }
  bool _causal_start_node( FFVar const& dom, size_t iel, size_t inode ) const
    { return _causal_rows_active() && dom.id() == _evolution_dom_var.id() && iel > 0 && inode == 0; }
  //! Several equations may share one DAG expression (an initial and a boundary condition T = T_ex), so rows are
  //! identified by their index, never by their variable.
  bool _causal_start_excludes( size_t eqn_ndx, FFVar const& dom, size_t iel, size_t inode ) const
    { return _causal_start_node( dom, iel, inode ) && _evoDerivRows.count( eqn_ndx ) > 0; }

  //! @brief as above for the C1 flux-jump coefficient.
  //! sigma1 is ENVIRONMENT-ONLY (CRONOS_SAT_SIGMA1, read once per process as).  It scales only
  //! C1 claims, which are minted only for IC_FLUX in a direction the block was NOT order-reduced; IC_AUTO never
  //! yields IC_FLUX, and an explicit IC_FLUX is resolved to IC_VALUE wherever reduction left an auxiliary.  In the
  //! full corpus sweep, all 1,010 plan censuses showed c1=0, so the 93 driver settings -- including the 41 that
  //! "disabled C1 penalties" with 0.0 -- never changed a number.
  //! The roadmap flags MMPDE10's as having "no recorded provenance"; it is in fact nearly half
  //! the corpus, so it is a convention rather than an oversight, and the audit should
  //! establish what that convention is rather than treat one driver as anomalous.
  double _eff_sat_sigma1() const
    {
      return 1.0;   // CRONOS_SAT_SIGMA1 '+R+'; 1.0 is reset()'s former default
    }

  //! @brief Compact summary of frozen interface-plan diagnostics.
  struct InterfacePlanDiagnosticSummary
  {
    size_t n_receiver_edges = 0;
    size_t n_weak_sat_terms = 0;
    size_t n_trace_constraints = 0;
    size_t n_tau_claim_keys = 0;

    // IC_TRACE-specific accounting.  These fields distinguish the number of
    // exact trace rows from the number of multiplier columns and make any
    // exact-only/hybrid trace reduction visible in diagnostics.
    size_t n_trace_variables = 0;
    size_t n_trace_constraint_declared = 0;
    size_t n_tau_claim_columns = 0;
    size_t n_exact_replacement_rows = 0;
    size_t n_physical_rows_dropped_by_trace = 0;
    size_t n_constraints_without_tau = 0;
    size_t n_tau_keys_without_receiver = 0;

    size_t n_receiver_principal = 0;
    size_t n_receiver_auxiliary = 0;
    size_t n_receiver_geometric = 0;

    size_t n_frozen_principal = 0;
    size_t n_frozen_auxiliary = 0;
    size_t n_frozen_geometric = 0;

    size_t n_duplicate_retained_trace_claims = 0;
    size_t n_duplicate_tau_claim_keys = 0;
    size_t n_duplicate_receiver_edge_identities = 0;

    size_t n_geometric_candidates = 0;
    size_t n_geometric_discarded_covered = 0;
    size_t n_geometric_emitted = 0;
  };

protected:
  //! @brief what SPQR reports about the factorization it just did.  Plain
  //! types only, so the declaration needs no CRONOS__WITH_SPQR guard.
  struct t_SpqrStat
  {
    size_t  nrow = 0, ncol = 0, nnzA = 0, rank = 0;
    int64_t nnzR = 0;          //!< cc.SPQR_istat[0] -- nonzeros in R, THE fill measure
    int64_t ncolsing = 0;      //!< cc.SPQR_istat[5]
    int64_t nrowsing = 0;      //!< cc.SPQR_istat[6]
    int     ordering = -1;     //!< cc.SPQR_istat[7] -- the ordering ACTUALLY used, which
                               //!< under SPQR_ORDERING_BEST names the concrete winner
    double  flops = 0., t_analyze = 0., t_factorize = 0.;
  };

public:
  //! @brief Column of the lifted auxiliary's value at the SAME spatial node as a claim's face (v1.5 profiles): the
  //! auxiliary has the state's layout with the evolution direction removed, so the face offset sface + oface*stride
  //! carries over, within the block of the face's spatial elements (@p ndx without the evolution direction).
  size_t _lifted_aux_col
    ( FFVar const& waux, std::map<FFVar,size_t,lt_FFVar> ndx, FFVar const& evo,
      size_t const sface, size_t const oface, size_t const stride )
    const
    { ndx.erase( evo );  return pos_state( waux, ndx ) + sface + oface * stride; }

  //! @brief The auxiliary standing for x(tau^+) of @p claim's state across its face, if that state is lifted there
  FFVar const* _lifted_claim_aux
    ( OCPlan::t_InterfaceClaimKey const& claim )
    const
    {
      if( _trnLift.empty() || claim.kind != OCPlan::CLAIM_EXACT_C0 ) return nullptr;
      if( static_cast<size_t>( _evolution_dom_var.id().second ) != claim.dom_id ) return nullptr;
      auto const it = _mDom.find( _evolution_dom_var );
      if( it == _mDom.end() || claim.iel_lo + 1 >= it->second.n_elem ) return nullptr;
      double const tb = it->second.elem_bnd[ claim.iel_lo + 1 ];
      for( auto const& L : _trnLift )
        if( static_cast<size_t>( L.state.id().second ) == claim.state_id
         && std::fabs( L.tau - tb ) <= 64.*DBL_EPSILON*std::max( 1., std::fabs( tb ) ) ) return &L.aux;
      return nullptr;
    }

  // ---------------------------------------------------------------------------------------------------------
  // rev309 -- INTERIM KNOBS, environment-only.  Each has ONE measured setting that is the desired behaviour;
  // no corpus driver ever set them, and a user should not.  They are NOT Options fields: reachable from the
  // environment (the C2/C2b arc was driven through them and every sweep note names them by these spellings),
  // but not part of the API.  The measured default lives here and nowhere else.
  // ---------------------------------------------------------------------------------------------------------
  //! @brief Derive a rescued interface coupling from the receiving row's constant coefficient c0 instead of forcing +1.
  //! Where that coefficient is not constant the forced value is kept and the edge is labelled coupling-undefined.
  //! Default on; CRONOS_RESCUE_C2=0 forces +1 everywhere.
  static bool _knob_RESCUE_C2
    ()
    { return static_cast<bool>( FFModel::Options::_env_flag( "CRONOS_RESCUE_C2", true ) ); }
  //! @brief IC_WEAK with multipliers -- an experiment, not a production setting.  0 (default) off; 1 imposes exactly
  //! every claim with a non-principal receiver edge; 2 only the claims whose coupling the rescue refused.  An integer 0..2,
  //! not a scale.  (The CRONOS_WEAK_TAU_RESCUED variable is retired, 2026-10-07: a measured constant.)
  static int _knob_WEAK_TAU_RESCUED
    ()
    { return 0; }
  //! @brief With WEAK_TAU_RESCUED set: an exact claim with a natural receiver feeds its multiplier only there.
  //! Default off.  (The CRONOS_WEAK_TAU_NATURAL variable is retired, 2026-10-07: a measured constant.)
  static bool _knob_WEAK_TAU_NATURAL
    ()
    { return false; }
  //! @brief Receiver policy for the multiplier under IC_TRACE (under IC_STRONG only via EXACT_NATURAL_STRONG).
  //! 0: every receiver carries the multiplier, rescued duplicates included, so the exact modes read the rescue's fabricated
  //! weight.  1: a claim with a natural (principal-symbol) tau-allowed receiver feeds its multiplier only there; a claim
  //! with none keeps its rescued receivers (silencing them would leave an empty multiplier column).  2 (default): as 1, but
  //! a natural receiver that is a LINK row does not qualify.  A silenced duplicate is restored wherever the claim's
  //! remaining column is dependent -- by the redundancy pass, or by the projector where that pass is vetoed.
  //! (The CRONOS_EXACT_NATURAL variable is retired, 2026-10-07: a measured constant.)
  static int _knob_EXACT_NATURAL
    ()
    { return 2; }
  //! @brief Apply EXACT_NATURAL under IC_STRONG as well (requires STRONG_PROJECT).  A silenced duplicate that turns
  //! out to be an eliminated tau's pivot row makes the Schur elimination fail; the plan is then rebuilt with every silenced
  //! claim restored.  Default on; (The CRONOS_EXACT_NATURAL_STRONG variable is retired, 2026-10-07.)
  static bool _knob_EXACT_NATURAL_STRONG
    ()
    { return true; }
  //! @brief Align IC_STRONG with IC_TRACE: when the trace projector (which IC_TRACE applies to reach a full-rank
  //! multiplier block) has support on Schur-eliminated taus, those taus are promoted to explicit on a promotion rebuild and
  //! the projector is applied to them as under IC_TRACE; the rest are eliminated as before.  Default on;
  //! (The CRONOS_STRONG_PROJECT variable is retired, 2026-10-07.)
  static bool _knob_STRONG_PROJECT
    ()
    { return true; }
  //! @brief IC_WEAK receiver policy for continuity claims.  0: every receiver edge carries a penalty (the rescued
  //! edges with a fabricated weight).  1: a claim that has a natural (principal-symbol) receiver is penalised only there;
  //! its rescued receiver terms are dropped; a claim with no natural receiver keeps them.  2 (default): as 1, but a natural
  //! receiver that is a LINK row (a reduction auxiliary's definition, which differentiates the claimed state) does not
  //! qualify -- a claim whose natural receivers are all LINK rows keeps every term.  Environment: CRONOS_WEAK_NATURAL_PENALTY.
  static int _knob_WEAK_NATURAL_PENALTY
    ()
    { return static_cast<int>( FFModel::Options::_env_int( "CRONOS_WEAK_NATURAL_PENALTY", 2, 0, 2 ) ); }
  //! @brief Corpus-sweep override for INTERFACE.DROP_POLICY (2026-10-12, WORKPLAN 4R.28).  CRONOS_DROP_POLICY=0|1|2
  //! forces DROP_VERIFY / DROP_REDERIVE / DROP_OFF for EVERY model, whatever the driver set; unset (the default, -1)
  //! leaves each driver's own value, so this changes nothing unless the variable is present.  Added because the drop
  //! is a per-driver option and a sweep needs one switch: on OCFE_Burger_move --form q (a SQUARE symbol, so the
  //! redundancy pass is not vetoed) the drop costs 6 rebuilds, 9 plan passes and a FAILED solve (60 iterations,
  //! 37.5 s), where DROP_OFF keeps all 216 claims, builds one plan and converges in 7 iterations (0.2 s).
  static int _knob_DROP_POLICY
    ()
    { return static_cast<int>( FFModel::Options::_env_int( "CRONOS_DROP_POLICY", -1, -1, 2 ) ); }

  // rev314 -- FOURTEEN more environment-only knobs (see the rev309 block for the rule).
  //! @brief Stage A -- WITHDRAWN, default off.  Reached k=0 on PDE2 and MBC but collapsed interface continuity
  //! (PDE2 IC_STRONG duplicate spread 3.469e-18 -> 2.133e-01).  Over-broad (per visit rather than per claim) and
  //! inapplicable to MMPDE9.  Kept as a reproducible negative result: CRONOS_LINK_PAIR_SUPPRESS=1 re-enables it.
  static bool _knob_INTERFACE_SUPPRESS_LINK_PAIR
    ()
    { return false; }   // CRONOS_LINK_PAIR_SUPPRESS (retired 2026-10-07, WORKPLAN 3.C) -- the withdrawn Stage A

  //! @brief Stage 3: suppress a C0 claim on a state that an ALL-mask IC/BC pins at this node (fixpoint closure).
  //! Row implication by physical rows, on the same footing as FIX 1 -- not Stage A's multiplier-column argument.
  //! Default on; CRONOS_NO_PIN_CLOSURE=1 disables it.
  static bool _knob_INTERFACE_SUPPRESS_PIN_CLOSURE
    ()
    { return true; }   // CRONOS_NO_PIN_CLOSURE (retired 2026-10-07, WORKPLAN 3.C)

  //! @brief Mint continuity claims from the C0 claim classification: a claim whose (state, direction) is IMPLIED,
  //! CLOSURE or NONE is not created, so nothing has to rescue it.  Exact modes (IC_STRONG, IC_TRACE) only -- under IC_WEAK
  //! a suppressed claim changes a penalty rather than a constraint.  Default off (the plan is then unchanged).
  //! (The CRONOS_CLAIM_C0 variable is retired, 2026-10-07: a measured constant.)
  static bool _knob_CLAIM_FROM_C0
    ()
    { return false; }

  //! @brief After every solve() under IC_TRACE with the projection applied, report beta = Z^T lambda and the
  //! violation eps*||beta|| = ||Z^T r_c^raw|| it absorbs.  Environment: CRONOS_AUDIT_BETA.
  static bool _knob_AUDIT_BETA
    ()
    { return Options::_env_flag( "CRONOS_AUDIT_BETA", false ); }

  //! @brief the solver's effective verbosity.  SOLVE.DISPLAY_LEVEL < 0 (the default) FOLLOWS the model's
  //! DISPLAY_LEVEL, resolved HERE, at read time -- not copied at reset(), because drivers set DISPLAY_LEVEL
  //! after construction and a copy would be stale.  Set SOLVE.DISPLAY_LEVEL >= 0 to make the solver's output
  //! independent of the model's.  Rule: each layer's output follows its own level.
  virtual int _solve_display_level
    ()
    const
    { return options.DISPLAY_LEVEL; }   // OCFESLV overrides: SOLVE.DISPLAY_LEVEL when set

  //! @brief every stored options object is OCFESLV's own type.
  std::shared_ptr<FFModel::EqnOptions> _make_equation_options
    ( FFModel::EqnOptions const& seed )
    const override
    {
      if( dynamic_cast<EqnOptions const*>( &seed ) ) return seed.clone();
      return std::make_shared<EqnOptions>( seed );
    }

  //! @brief Return interface plan.
  mc::OCPlan const& interface_plan
    ()
    const
    noexcept
    { return _plan; }

protected:
  //! @brief Claim key freeze/unfreeze — delegate to OCPlan static helpers.
  static t_FrozenInterfaceClaimKey _freeze_interface_claim_key
    ( t_InterfaceClaimKey const& k )
    { return OCPlan::freeze_claim_key( k ); }
  static t_InterfaceClaimKey _unfreeze_interface_claim_key
    ( t_FrozenInterfaceClaimKey const& k )
    { return OCPlan::unfreeze_claim_key( k ); }

  //! @brief The working state / domain variable carrying a given id, or nullptr.  The plan keys states and
  //! domains by id; C0 keys them by FFVar, so the C1 filter needs the map back.
  FFVar const* _state_by_id
    ( size_t const id )
    const
    {
      for( auto const& [v, dom] : _mVar ){ (void)dom; if( v.id().second == static_cast<long>( id ) ) return &v; }
      return nullptr;
    }

  //! @brief The working domain variable carrying a given id, or nullptr.
  FFVar const* _domain_by_id
    ( size_t const id )
    const
    {
      for( auto const& [v, dom] : _mDom ){ (void)dom; if( v.id().second == static_cast<long>( id ) ) return &v; }
      return nullptr;
    }

public:
  //! @brief The coupling decisions of the last setup, for a tailored driver to assert on.
  std::map<std::pair<size_t,size_t>,t_CouplingRecord> const& coupling_log() const { return _couplingLog; }

protected:
  //! @brief copy PlanDecisions into the draft (and _keepExplicitReason);
  //! the only writer of decision fields.
  void _apply_plan_decisions( PlanDecisions const& dec, InterfacePlanDraft& draft ) const;

  //! @brief Reset the sealed plan to the empty state.
  //!
  //! Called whenever the symbolic model or classification changes so
  //! that the stale plan cannot be consumed by eval()/deriv().
  void _clear_frozen_interface_consumers() const
    { _plan.reset(); _interfacePlanBuildAllowed = false; }

  //! @brief Validate a completed draft before commit (delegates to OCPlan).
  bool _validate_interface_plan_draft
    ( InterfacePlanDraft const& draft ) const
    { return _plan.validate_draft( draft, options.INTERFACE.IMPOSITION ); }

  //! @brief Atomically publish a validated draft (delegates to OCPlan).
  void _commit_interface_plan_draft
    ( InterfacePlanDraft&& draft )
    {
      _plan.commit_draft( std::move( draft ) );
      // Propagate plan dimensions back into OCFESLV scalar members so that
      // the rest of the class can continue to use _nCollEqn etc. directly.
      _nCollEqn        = _plan.n_coll_eqn();
      _nCollVar        = _plan.n_coll_var();
      _traceVarOffset  = _plan.trace_var_offset();
      _nTraceConstraint = _plan.n_trace_constraint();
      // For IC_STRONG, tau variables are internal Schur variables not
      // appended to var[], so the user-visible trace count is 0.  The
      // plan-internal accepted tau count remains available via
      // _plan.n_trace_var() for the Schur builder and diagnostics.
      // IC_STRONG-hybrid exception: the explicit (kept) taus ARE appended to
      // var[] (like IC_TRACE taus), so the user-visible trace count is the
      // explicit subset.  The eliminated taus stay internal and are excluded.
      // n_trace_var_explicit() is 0 for every non-hybrid plan, so pure
      // IC_STRONG still reports 0.
      _nTraceVar       = ( options.INTERFACE.IMPOSITION == Options::IC_STRONG )
                       ? _plan.n_trace_var_explicit() : _plan.n_trace_var();
    }

  //! @brief Principal-symbol block-eligibility predicate for dropping a
  //! continuity claim (item 11).  True iff the block's principal symbol makes a
  //! value/flux continuity claim redundant: parabolic/parabolic-closure (flux
  //! over-specified by value + LINK) or scalar hyperbolic (single real
  //! characteristic).  False for system-hyperbolic (>=2 real characteristics),
  //! degenerate/unresolved, merely-spatially-characteristic, or unknown blocks
  //! (keep all).  Single source of truth shared by the redundancy gate and the
  //! value-slaved drop, so neither can drop a claim in an ineligible block.
  bool _block_drop_eligible( int block_id ) const;

  //! @brief v2 interface decision layer (Phase 1; defined inline at end of this header).
  //!
  //! Sets the per-edge exactness flags (weak_sat_allowed, hard_physical_row,
  //! tangential_physical_boundary, trace_tau_allowed) per the rewritten rule:
  //! every non-BC continuity receiver becomes an exact tau-able value-C0
  //! constraint.  coupling/source are left untouched so IC_WEAK is preserved.
  void _apply_interface_edge_flags_v2
    ( InterfaceReceiverEdge& edge,
      FFModel::EqnOptions const&      eqnopt,
      bool                   aux_alg_receiver,
      bool                   is_geo_c0 ) const;

  //! @brief The classification is about to be recomputed: the SAT coupling table must be rebuilt before any
  //! interface plan is built or evaluated.
  void _on_classify_begin
    ()
    override
    { _satCouplingBuilt = false; }

  //! @brief Build per-element physical node cache for domain variables.
  void _setup_node_cache
    ();

  //! @brief Check consistency of a collocation model container set.
  //!
  //! When optional output pointers are supplied, this also rebuilds the
  //! equation/output subgraph caches and row/type counts used by setup().
  bool _check_consistency
    ( FFGraph* dag, std::vector<FFVar> const& vCst, std::vector<double> const* vCstVal,
      t_Dom const& mDom, t_Var const& mVar, t_Var const& mInp,
      t_Eqns const& mEqn, t_Fcts const& mFct,
      size_t* nCollEqn = nullptr, std::vector<FFSubgraph>* sgEqn = nullptr,
      size_t* nCollFct = nullptr, std::vector<FFSubgraph>* sgFct = nullptr )
    const;

  //! @brief True when order reduction introduced a LINK equation in a block/direction.
  bool _block_has_aux_link_in_direction
    ( int block_id, FFVar const& dom )
    const;

  //! @brief True when an auxiliary state represents a derivative in a given direction.
  bool _auxiliary_state_in_direction
    ( FFVar const& state, FFVar const& dom )
    const;

  //! @brief True when an equation expression contains an auxiliary generated in a given direction.
  bool _equation_has_aux_link_in_direction
    ( FFVar const& eqn, int block_id, FFVar const& dom )
    const;

  //! @brief True when IC_WEAK should use exact replacement in a reduced direction.
  bool _weak_exact_reduced_interface_in_direction
    ( int block_id, FFVar const& dom )
    const;

  //! @brief True when this block has an actual INITIAL condition on a domain face.
  bool _block_has_initial_face
    ( int block_id, FFVar const& dom, int lim )
    const;

  //! @brief True when this block has an INITIAL or BOUNDARY trace on a domain face.
  bool _block_has_trace_face
    ( int block_id, FFVar const& dom, int lim )
    const;

  //! @brief True when a physical INITIAL/BOUNDARY trace directly contains a
  //! specific state as a value variable on a given face.  Auxiliary/flux rows
  //! do not count as value traces for their primitive state.
  bool _state_has_value_trace_face
    ( int block_id, size_t state_id, FFVar const& dom, int lim )
    const;

public:
  static char const* interface_edge_source_name
    ( OCPlan::InterfaceEdgeSource const source ) noexcept
  {
    switch( source ){
      case OCPlan::InterfaceEdgeSource::PRINCIPAL_SYMBOL:    return "PRINCIPAL_SYMBOL";
      case OCPlan::InterfaceEdgeSource::AUXILIARY_ALGEBRAIC: return "AUXILIARY_ALGEBRAIC";
      case OCPlan::InterfaceEdgeSource::GEOMETRIC_C0:        return "GEOMETRIC_C0";
    }
    return "UNKNOWN";
  }

  //! @brief Human-readable name for an interface imposition mode.
  static char const* interface_imposition_name
    ( int const imposition ) noexcept
  {
    switch( imposition ){
      case Options::IC_WEAK:   return "IC_WEAK";
      case Options::IC_STRONG: return "IC_STRONG";
      case Options::IC_TRACE:  return "IC_TRACE";
    }
    return "UNKNOWN";
  }

  //! @brief Return a compact summary of the frozen interface plan.
  InterfacePlanDiagnosticSummary interface_plan_diagnostic_summary
    () const
  {
    InterfacePlanDiagnosticSummary out;

    auto const& edges       = _plan.receiver_edges();
    auto const& weak_terms  = _plan.weak_sat_terms();
    auto const& trace_terms = _plan.trace_constraints();
    auto const& tau_claims  = _plan.tau_claims_frozen();

    out.n_receiver_edges     = edges.size();
    out.n_weak_sat_terms     = weak_terms.size();
    out.n_trace_constraints  = trace_terms.size();
    out.n_tau_claim_keys     = tau_claims.size();
    out.n_trace_variables    = _plan.n_trace_var();
    out.n_trace_constraint_declared = _plan.n_trace_constraint();
    out.n_exact_replacement_rows = _plan.exact_replacement_rows().size();
    if( _plan.frozen_options().imposition_type == Options::IC_TRACE )
      out.n_physical_rows_dropped_by_trace = out.n_exact_replacement_rows;
    for( auto const& kv : tau_claims ) out.n_tau_claim_columns += kv.second.size();
    out.n_geometric_candidates = _plan.diagnostics().geometric_candidates;
    out.n_geometric_discarded_covered = _plan.diagnostics().geometric_discarded_covered;
    out.n_geometric_emitted = _plan.diagnostics().geometric_emitted;

    for( auto const& e : edges ){
      switch( e.source ){
        case OCPlan::InterfaceEdgeSource::PRINCIPAL_SYMBOL:    ++out.n_receiver_principal; break;
        case OCPlan::InterfaceEdgeSource::AUXILIARY_ALGEBRAIC: ++out.n_receiver_auxiliary; break;
        case OCPlan::InterfaceEdgeSource::GEOMETRIC_C0:        ++out.n_receiver_geometric; break;
      }
    }
    for( auto const& t : weak_terms ){
      switch( t.source ){
        case OCPlan::InterfaceEdgeSource::PRINCIPAL_SYMBOL:    ++out.n_frozen_principal; break;
        case OCPlan::InterfaceEdgeSource::AUXILIARY_ALGEBRAIC: ++out.n_frozen_auxiliary; break;
        case OCPlan::InterfaceEdgeSource::GEOMETRIC_C0:        ++out.n_frozen_geometric; break;
      }
    }

    std::map<OCPlan::t_InterfaceClaimKey,size_t> retained_trace_count;
    for( auto const& tr : trace_terms ) ++retained_trace_count[ tr.edge ];
    for( auto const& kv : retained_trace_count )
      if( kv.second > 1 ) ++out.n_duplicate_retained_trace_claims;

    std::set<OCPlan::t_FrozenInterfaceClaimKey> tau_key_seen;
    for( auto const& kv : tau_claims )
      if( !tau_key_seen.insert( kv.first ).second ) ++out.n_duplicate_tau_claim_keys;

    for( auto const& tr : trace_terms ){
      auto const fk = OCPlan::freeze_claim_key( tr.edge );
      auto const it = tau_claims.find( fk );
      if( it == tau_claims.end() || it->second.empty() )
        ++out.n_constraints_without_tau;
    }

    for( auto const& kv : tau_claims ){
      bool has_receiver = false;
      for( auto const& t : weak_terms ){
        if( !t.trace_tau_allowed ) continue;
        if( OCPlan::freeze_claim_key( t.claim_key ) == kv.first ){
          has_receiver = true;
          break;
        }
      }
      if( !has_receiver ) ++out.n_tau_keys_without_receiver;
    }

    std::map<
      std::tuple< size_t, std::vector<std::pair<size_t,size_t>>, size_t,
                  OCPlan::t_FrozenInterfaceClaimKey, int, int >,
      size_t > receiver_edge_identity_count;
    for( auto const& e : edges ){
      auto const id = std::make_tuple(
        e.row_id, e.receiver_block_el, e.receiver_emit_flat,
        OCPlan::freeze_claim_key( e.claim_key ),
        static_cast<int>( e.source ), static_cast<int>( e.sat_type ) );
      ++receiver_edge_identity_count[ id ];
    }
    for( auto const& kv : receiver_edge_identity_count )
      if( kv.second > 1 ) ++out.n_duplicate_receiver_edge_identities;

    return out;
  }

protected:
  //! @brief Flat state-DOF (column) indices belonging to evolution element
  //! @p iel_evo -- every evolution-domain-bearing state, over all spatial
  //! elements and nodes.  pos_state flattens the element multi-index in
  //! lt_FFVar (DAG-creation) order, so which domain is innermost depends on
  //! variable creation order, NOT on being the evolution direction (drivers
  //! create t before z, so t is currently innermost => strided slabs).  This
  //! accessor does NOT depend on that ordering: it fixes the evolution element,
  //! enumerates the spatial multi-index, and delegates flattening to pos_state,
  //! so the returned set is correct whether a slab is contiguous or strided.
  //! States NOT on the evolution domain are excluded (they are not slab-local).
  std::vector<size_t> _evolution_slab_cols
    ( size_t iel_evo ) const;

public:
  //! @brief Read the collocated DOF values of distributed input @p input from @p inp
  //! in the same element-major / node-inner ordering set_input_values() writes.  The
  //! inverse of set_input_values(); returns an empty vector on unknown input.
  std::vector<double> get_input_values
    ( FFVar const& input, double const* inp ) const;

  //! @brief Clear the protected-claims set (call before a fresh problem/run).
  void clear_protected_interface_claims()
    { _strong_protected_claims.clear(); }

protected:
  //! @brief differential-core SHAPE from the USER-side structures.
  //!
  //! Rows   = equations carrying at least one d/d(evolution) of a state.
  //! Cols   = distinct states so differentiated ( _states_with_evolution_derivative ).
  //! Square <=> rows == cols, which is the classifier's own marchability criterion
  //! (the retired DIFFERENTIAL_RECTANGULAR was marchable only when the domain was
  //! user-supplied), replicated where it can actually be consulted.
  //!
  //! IT WILL NOT ALWAYS AGREE WITH THE [classify] LINE'S sym=RxC, and that is expected:
  //! this runs on _usr._mEqnUsr BEFORE localisation and BEFORE _reduce_order, while the
  //! classifier counts RETAINED rows after both.  Both numbers are reported so a
  //! disagreement is visible rather than assumed away.  A model where they differ is a
  //! finding about the replica, not a licence to trust it.
  void _differential_core_shape
    ( size_t& rows, size_t& cols ) const;

  //! @brief True if output expression @p fctvar contains an OpI (FFIntegral) quadrature
  //! over the evolution variable -- i.e. an evolution-direction integral whose per-element
  //! contributions must be SUMMED across the march.  All other outputs (terminal /
  //! final-time / spatial-only) take the last element's value.
  bool _fct_is_evolution_integral( FFVar const& fctvar ) const;

  //! @brief Adjoint of the capture quadrature (Stage 1b): given lam_cap[R] (the adjoint on each cap's
  //! value, indexed like R's spatial dofs), scatter it onto the state adjoint @p ls at R's collocation
  //! dofs with the SAME LGR weights jac*w_ref[m] used by _fill_capture_sens, summed over evolution
  //! elements.  Transpose of d(cap)/dp = INT jac*w_ref*s[R].  ACCUM only (LATCH handled inline).
  void _fill_capture_sens_adjoint( std::map<FFVar,std::vector<double>,lt_FFVar> const& lam_cap,
                                   double* ls ) const;

  //! @brief LATCH transpose of the forward latch d(cap)/dp = eval_colloc(src,tau): scatter ONE latch
  //! capture's cap-adjoint @p lamC onto @p ls at TAU's evolution element with the value Lagrange
  //! weights L_m(tau) (NOT the ACCUM quadrature jac*w_ref).  Per-capture (tau not lost to
  //! source-aggregation).  No-ops unless tau lies in the current evolution domain (window under
  //! marching, full horizon monolithic).
  void _fill_capture_sens_adjoint_latch( FFVar const& src, double tau,
                                         std::vector<double> const& lamC, double* ls,
                                         int const side = FFDom::MINUS ) const;

  //! @brief Orientation of IC_TRACE multiplier receiver terms.
  //!
  //! Exact trace rows use the canonical orientation g = u_left(+1)-u_right(-1).
  //! Receiver rows on the left/low element get +lambda and rows on the
  //! right/high element get -lambda, matching the transpose saddle-point pair.
  static double _trace_receiver_orient
    ( bool current_element_is_right_of_interface )
    { return current_element_is_right_of_interface ? -1.0 : 1.0; }

  //! @brief Canonical block id for a physical IC_TRACE claim slot.
  int _trace_claim_block
    ( size_t state_id, size_t dom_id ) const;

  //! @brief Default SAT spectral penalty parameter for domain d.
  //! Returns tau = (n_node - 1)^2 / element_width.
  //! The effective C0 penalty coefficient is  INTERFACE.SAT_SIGMA0 * tau,
  //! and the effective C1 penalty coefficient is  CRONOS_SAT_SIGMA1 * tau / element_width.
  double _sat_tau
    ( FFDom const& dom, size_t const iel )
    const
    { double const nm1 = static_cast<double>( dom.n_node - 1 );
      return nm1 * nm1 / dom.elem_width( iel ); }

  //! @brief Default SAT spectral penalty parameter using representative width; retained for setup paths without a block index.
  double _sat_tau
    ( FFDom const& dom )
    const
    { double const nm1 = static_cast<double>( dom.n_node - 1 );
      return nm1 * nm1 / dom.w_elem; }

  //! @brief Interface-symmetric finite-element width shared by elements
  //! iel_lo and iel_lo+1.  Returns the harmonic mean of the two adjoining
  //! element widths, so that both elements bordering an internal interface
  //! evaluate an identical SAT penalty parameter.  This is what makes the
  //! two-sided IC_WEAK / IC_VALUE penalty discretely conservative on
  //! nonuniform finite-element grids.
  double _sat_iface_width
    ( FFDom const& dom, size_t const iel_lo )
    const
    { double const wL = dom.elem_width( iel_lo );
      double const wR = dom.elem_width( iel_lo + 1 );
      double const s  = wL + wR;
      return ( s > 0. ) ? ( 2. * wL * wR / s ) : 0.; }

  //! @brief Build _satCoupling for the current classification signature.
  //! Called by setup immediately after classification; evaluation must not
  //! rely on this as a hidden topology-building fallback.
  void _ensure_sat_coupling
    () const;

  //! @brief Dry structural pass of the residual assembler that verifies the
  //! strong-imposition interface/corner logic keeps the collocation system
  //! square: LOOP 1 interface rows + LOOP 2 equation nodes == _nCollEqn,
  //! and LOOP 1 emissions == LOOP 2 interface-replacement drops.  Index-based and value-
  //! independent.  Surfaces a squareness violation at setup()/classify_pde()
  //! time instead of at the first eval()/deriv() call.
  bool _audit_strong_interface_rows
    ();

public:
  //! @brief Name of a state by its variable id (empty when not a state).
  std::string state_name_of_id
    ( size_t sid )
    const
    { for( auto const& kv : _mVar )
        if( static_cast<size_t>( kv.first.id().second ) == sid ) return kv.first.name();
      return std::string(); }

  //! @brief Sparse modified Gram-Schmidt column rank, 1e-9 relative cut.
  //!
  //! Same algorithm and tolerance as the rank_of() lambda in interface_rank_symmetry, which
  //! additionally reports the dependent indices.  Used by the projection phase so that the
  //! k it computes is the k the audit reports.
  static size_t _sparse_column_rank
    ( std::vector<std::vector<std::pair<size_t,double>>> V )
    {
      typedef std::vector<std::pair<size_t,double>> SpV;
      std::vector<SpV> basis;
      size_t rank = 0;
      for( size_t j = 0; j < V.size(); ++j ){
        SpV c = V[j];
        std::sort( c.begin(), c.end(),
                   []( std::pair<size_t,double> const& a,
                       std::pair<size_t,double> const& b ){ return a.first < b.first; } );
        double n0 = 0.;
        for( size_t k = 0; k < c.size(); ++k ) n0 += c[k].second * c[k].second;
        n0 = std::sqrt( n0 );
        if( n0 <= 64. * DBL_EPSILON ) continue;
        for( size_t bi = 0; bi < basis.size(); ++bi ){
          SpV const& q = basis[bi];
          double d = 0.;
          size_t a = 0, b = 0;
          while( a < c.size() && b < q.size() ){
            if     ( c[a].first < q[b].first ) ++a;
            else if( c[a].first > q[b].first ) ++b;
            else { d += c[a].second * q[b].second; ++a; ++b; }
          }
          if( d == 0. ) continue;
          SpV out;
          out.reserve( c.size() + q.size() );
          a = 0; b = 0;
          while( a < c.size() && b < q.size() ){
            if     ( c[a].first < q[b].first ){ out.push_back( c[a] ); ++a; }
            else if( c[a].first > q[b].first ){ out.push_back(
                       std::make_pair( q[b].first, -d * q[b].second ) ); ++b; }
            else { out.push_back( std::make_pair(
                     c[a].first, c[a].second - d * q[b].second ) ); ++a; ++b; }
          }
          for( ; a < c.size(); ++a ) out.push_back( c[a] );
          for( ; b < q.size(); ++b )
            out.push_back( std::make_pair( q[b].first, -d * q[b].second ) );
          c.swap( out );
        }
        double nr = 0.;
        for( size_t k = 0; k < c.size(); ++k ) nr += c[k].second * c[k].second;
        nr = std::sqrt( nr );
        if( nr <= 1e-9 * std::max( 1.0, n0 ) ) continue;
        double const inv = 1.0 / nr;
        for( size_t k = 0; k < c.size(); ++k ) c[k].second *= inv;
        basis.push_back( c );
        ++rank;
      }
      return rank;
    }

  //! @brief 2026-10-15 (WORKPLAN 4R.21): the DEFINING ROW of each state, state id -> equation index.
  //! A declaration (EqnOptions::defines) always wins.  Otherwise INFERRED: among the participating rows that are
  //! not evolution rows, those whose FFInv letter in the state is L or S -- linear, or separably linear (the state
  //! appears once, multiplied by a coefficient that may involve other states).  The FFPartial overload types its
  //! OPERAND U, so a row that differentiates the state can never be L/S in it: "pointwise in the state" comes for
  //! free, and no proxy substitution is needed.  Tie-break: a unique L beats any S; a unique S is taken; otherwise
  //! the state is AMBIGUOUS and gets no inferred row (declare it).  Measured on OCFE_Burger_move --form q:
  //! q -> QDEF (S; the PDE is excluded as an evolution row), g -> DEF (L), xg -> XGDEF (L over DEF's S).
  std::map<size_t,size_t> _defining_rows() const;

  //! @brief true once CRONOS_AUDIT_TAUCOL has been read.  Static so the
  //! getenv cost is paid once per process, matching kEvoFace / kEvoFaceH.
  //! @brief report-format revision, embedded in the [taucol/plan] banner so a
  //! binary can be checked for the FORMAT it emits, not merely for having the instrument.
  static char const* taucol_instrument_rev() { return "instrument rev108"; }

  //! @brief resolve a state id to its printed name (same lookup the
  //! IC_TRACE detail report uses; duplicated rather than shared because that one
  //! is a lambda local to a print_details branch).
  std::string _taucol_state_name( size_t const sid ) const
    {
      for( auto const& kv : _mVar )
        if( static_cast<size_t>( kv.first.id().second ) == sid ){
          std::ostringstream osn; osn << kv.first; return osn.str();
        }
      std::ostringstream osn; osn << "state#" << sid; return osn.str();
    }

  //! @brief resolve a domain id to its printed name.
  std::string _taucol_dom_name( size_t const did ) const
    {
      for( auto const& kv : _mDom )
        if( static_cast<size_t>( kv.first.id().second ) == did ){
          std::ostringstream osn; osn << kv.first; return osn.str();
        }
      std::ostringstream osn; osn << "dom#" << did; return osn.str();
    }

  //! @brief printable form of a physical residual row key.
  std::string _taucol_row_str( FrozenPhysicalResidualRowKey const& rk ) const
    {
      std::ostringstream o;
      o << "row_id=" << rk.row_id << " blk=[";
      for( size_t i = 0; i < rk.block_el.size(); ++i )
        o << ( i ? "," : "" ) << rk.block_el[i].first << ":" << rk.block_el[i].second;
      o << "] flat=" << rk.flat;
      return o.str();
    }

  //! @brief Decode a state-coefficient column index into "state[el=..](nd=..)@coord".
  //!
  //! Layout mirrors pos_state(): each state owns a contiguous range starting at its
  //! _blkVar offset; within it, element blocks of blk_size coefficients are strided with
  //! the FIRST domain direction fastest, and within a block the node multi-index is
  //! flattened the same way (see the receiver_flat construction).  Returns "?" for a
  //! column outside every state range (a trace/tau column under IC_TRACE).
  std::string decode_state_column
    ( size_t col )
    const
    {
      for( auto const& [v, dom] : _mVar ){
        auto itb = _blkVar.find( v );
        if( itb == _blkVar.end() ) continue;
        size_t const base     = std::get<0>( itb->second );
        size_t const blk_size = std::get<1>( itb->second );
        if( !blk_size ) continue;
        size_t span = blk_size;
        for( auto const& d : dom ){
          auto itd = _mDom.find( d );
          if( itd == _mDom.end() ){ span = 0; break; }
          span *= itd->second.n_elem;
        }
        if( !span || col < base || col >= base + span ) continue;

        size_t local = col - base;
        size_t blkid = local / blk_size, nodeoff = local % blk_size;

        std::ostringstream ss;
        ss << v.name() << "[";
        std::ostringstream nn, cc;
        bool first = true;
        for( auto const& d : dom ){
          auto const& dd = _mDom.at( d );
          size_t const iel = blkid % ( dd.n_elem ? dd.n_elem : 1 );
          blkid /= ( dd.n_elem ? dd.n_elem : 1 );
          size_t const ind = nodeoff % ( dd.n_node ? dd.n_node : 1 );
          nodeoff /= ( dd.n_node ? dd.n_node : 1 );
          if( !first ){ ss << ","; nn << ","; cc << ","; }
          first = false;
          ss << d.name() << "=" << iel;
          nn << ind;
          // Physical coordinate of the node: reference nodes live in [-1,1].
          auto const& xr = dd.lgnodes();
          if( ind < xr.size() && iel < dd.n_elem )
            cc << d.name() << "="
               << std::setprecision(4) << ( dd.elem_lo( iel )
                  + 0.5 * ( xr[ind] + 1.0 ) * dd.elem_width( iel ) );
          else
            cc << d.name() << "=?";
        }
        ss << "](nd=" << nn.str() << ")@" << cc.str();
        return ss.str();
      }
      return std::string( "?col" ) + std::to_string( col );
    }

  //! @brief Ordered vector of state variables as registered after setup()
  //! Includes auxiliary variables introduced by reduce_order().
  std::vector<FFVar> const& states_colloc
    ()
    const
    { return _vVar; }

  //! @brief Return vector of node coordinates for the domain of variable var
  std::vector<std::vector<double>> node_colloc
    ( FFVar const& var )
    const;

  //! @brief 2026-10-02: the ROW-level causal design is active -- every element's first time node is an initial node:
  //! the `ALL - LB` rows are not collocated there and each state with an evolution derivative gets an exact
  //! continuity row against the previous element (as a marching window's start).  For now under the option /
  //! environment only (the jumping-input gate keeps the claims-based causal block).
  //! @brief The gate that admits the causal evolution-direction treatment: an INPUT distributed in the evolution
  //! direction with an interior evolution interface whose declared continuity is below C0 (undeclared: may jump).
  bool _evolution_input_may_jump() const
  {
    if( !_evolution_dom_set || !_evolution_dom_var.dag() || !_evoCollocated ) return false;   // a jump in t is the integrator's business under MoL
    for( auto const& kv : _mInp ){
      FFVar const& inp = kv.first;
      bool here = false;
      for( auto const& dv : kv.second )
        if( dv.id() == _evolution_dom_var.id() ){ here = true; break; }
      if( !here ) continue;                                          // not distributed in t
      if( _input_dom( inp, _evolution_dom_var ).n_elem < 2 ) continue; // no interior t-interface
      int level = static_cast<int>( DISCONTINUOUS );                 // default when undeclared
      auto itc = _mInpCont.find( inp );
      if( itc != _mInpCont.end() ){
        auto itd = itc->second.find( _evolution_dom_var );
        if( itd != itc->second.end() ) level = static_cast<int>( itd->second );
      }
      if( level < static_cast<int>( CONTINUOUS_C0 ) ) return true;   // may jump -> kinked
    }
    return false;
  }

  //! @brief base offset of a state's data in the flat collocation vector.
  //!
  //! THREE WRONG ANSWERS PRECEDED THIS ONE; do not replace it without reading them.
  //!
  //!  1. Walking the state list, accumulating `node_colloc( st ).size()`.  That adds
  //!     ZERO for a SCALAR state (node_colloc returns an empty vector) even though the state
  //!     owns a var[] slot, so every state after the first scalar one was read short.  It
  //!     produced plausible numbers, not errors: OCFE_scalar_nest reported a 2.5e-01
  //!     "continuity violation" that was really a'(z) sampled at the wrong two nodes.
  //!
  //!  2. Using `pos_state( st, {} )`.  pos_state does `ndx_el.at(v)` for EVERY domain
  //!     of the state, so an empty index map throws std::out_of_range for every DISTRIBUTED
  //!     state.  The empty-map idiom is only valid for scalars -- which is the branch it was
  //!     copied from (_init_var_block's `if( dom.empty() )` path).  This aborted PDE38 and
  //!     scalar_nest after a successful solve.
  //!
  //!  3. Taking _blkVar's blk_pos directly would be correct today but duplicates pos_state's
  //!     contract, so a change to the layout would silently desynchronise the two.
  //!
  //! Correct: pos_state with every domain of the state pinned to element 0, which is block 0
  //! for scalars and distributed states alike, and delegates the layout entirely.
  size_t _state_base_offset
    ( FFVar const& st )
    const
    {
      std::map<FFVar,size_t,lt_FFVar> ndx0;
      auto const itv = _mVar.find( st );
      if( itv != _mVar.end() )
        for( auto const& d : itv->second ) ndx0[d] = 0;
      return pos_state( st, ndx0 );
    }

protected:
  //! @brief condition estimate on the retained R11 block of a weak
  //! rank-revealing QR; deflates the rank by the number of sigma(R11) below tol_abs.
  static void _refine_rank_from_R
    ( arma::sp_mat const& R11, arma::sp_mat const& R12, double tol_abs,
      size_t& rank, double& gap, double& sigma_cut );

  //! @brief live (staircase-pivot) columns of a SQUEEZED SPQR R, in column
  //! order.  ONE definition for both callers (numeric_rank and _null_basis once ran two
  //! textual copies of the same loop); `sorted` selects the O(1)
  //! last-entry test, everything else takes the column scan.  Pointers are the
  //! raw CSC triple of the cholmod R (no cholmod types here, so the declaration needs
  //! no SPQR guard).
  static void _staircase_live_columns
    ( int64_t const* rp, int64_t const* rI, double const* rx,
      size_t ncol, size_t r, bool sorted, std::vector<size_t>& live );

  //! @brief dense-SVD rank oracle (CRONOS_RANK_ORACLE=1, default OFF).  Below
  //! CRONOS_RANK_ORACLE_MAXN it factorizes the SAME matrix densely at the SAME
  //! tol_abs and prints AGREE/DISAGREE against the refined rank.  Report-only: it
  //! never changes a rank, a verdict, or a plan.
  static void _rank_oracle_check
    ( arma::sp_mat const& M, double tol_abs, size_t rank_backend, size_t rank_refined,
      char const* backend );

  //! @brief one line of SPQR's own statistics (CRONOS_SPQR_STAT=1, or under
  //! CRONOS_DEBUG_RANK).  Report-only; default off.
  static void _spqr_stat_report( t_SpqrStat const& S, char const* tag );

public:
  //! @brief the SOLVER's per-equation options.  Derives from FFModel::EqnOptions (role,
  //! participate_in_classification, block_id) and adds how OCFESLV imposes each equation at interfaces.
  //! OCFESLV's factory guarantees that EVERY stored equation -- user-added or minted by reduce_order() or
  //! setup() -- carries this type, which is what makes _sopt() below safe.
  struct EqnOptions
  : public FFModel::EqnOptions
  {
    Options::InterfaceType interface_type;              //!< local interface-condition type (IC_AUTO: resolved)
    bool                   receive_sat;                 //!< may receive weak (SAT) interface penalties
    bool                   donate_for_state_continuity; //!< may be donated as a duplicate trace row

    //! @brief Role-derived defaults: INITIAL, BOUNDARY and
    //! INTERFACE are donated, DIAGNOSTIC receives no SAT, LINK is never donated.
    EqnOptions
      ( EqnRole r = EqnRole::AUTO, int b = 0, Options::InterfaceType it = Options::IC_AUTO )
      : FFModel::EqnOptions( r, b ), interface_type( it ), receive_sat( true ), donate_for_state_continuity( false )
      { _solver_role_defaults(); }

    //! @brief Explicit flags; the role invariants still apply: INTERFACE never enters classification, LINK and
    //! DIAGNOSTIC are never donated.
    EqnOptions
      ( EqnRole r, int b, Options::InterfaceType it, bool classify, bool sat, bool donate = false )
      : FFModel::EqnOptions( r, b ), interface_type( it ), receive_sat( sat ), donate_for_state_continuity( donate )
      {
        participate_in_classification = classify;
        if( role == EqnRole::INTERFACE ) participate_in_classification = false;
        if( role == EqnRole::LINK || role == EqnRole::DIAGNOSTIC ) donate_for_state_continuity = false;
      }

    //! @brief Promote a model-side options object: the solver fields take their role-derived defaults.  Used
    //! by OCFESLV's factory for every equation the model mints itself.
    explicit EqnOptions
      ( FFModel::EqnOptions const& base )
      : FFModel::EqnOptions( base ), interface_type( Options::IC_AUTO ), receive_sat( true ),
        donate_for_state_continuity( false )
      { _solver_role_defaults(); }

    EqnOptions( EqnOptions const& other ) = default;
    EqnOptions& operator=( EqnOptions const& other ) = default;

    //! @brief Contract: a derived options type MUST override clone(), or copies drop these fields.
    std::shared_ptr<FFModel::EqnOptions> clone
      ()
      const override
      { return std::make_shared<EqnOptions>( *this ); }

    void _solver_role_defaults
      ()
      {
        if( role == EqnRole::INITIAL || role == EqnRole::BOUNDARY || role == EqnRole::INTERFACE )
          donate_for_state_continuity = true;
        if( role == EqnRole::DIAGNOSTIC ) receive_sat = false;
        if( role == EqnRole::LINK ) donate_for_state_continuity = false;
      }
  };

protected:
  struct t_InterfaceTensorCornerDim {
    size_t       face_stride = 0;
    FFDom const* pdom        = nullptr;
    FFVar const* pvar        = nullptr;
    size_t       n_node      = 0;
    size_t       n_elem      = 0;
    size_t       iel         = 0;
    size_t       inode       = 0;
  };

public:
  //! @brief Expected behaviour for built-in interface-plan diagnostics.
  enum class InterfacePlanGeometricExpectation
  {
    ANY = 0,    //!< Do not check the number of geometric-C0 edges.
    ZERO,       //!< Require no geometric-C0 receiver/frozen terms.
    NONZERO     //!< Require at least one geometric-C0 receiver/frozen term.
  };

  //! @brief sparse residual Jacobian assembled at a supplied var vector by
  //! assemble() / assemble_window().  Solver-independent; the sparsity pattern is
  //! cached inside so repeated calls at new points pay only the value evaluation.
  struct t_JacAssembly
  {
    size_t n = 0, m = 0;          //!< equation rows, variable columns (states + trace/tau)
    size_t n_tau = 0, tau0 = 0;   //!< tau column count and offset (n_tau = 0 under IC_WEAK)
    size_t nnz = 0;               //!< pattern nonzeros (state + input columns, as deriv() reports)
    arma::sp_mat J;               //!< n x m; = diag(rw)*J*diag(cw) iff equilibrated
    arma::vec    r;               //!< residual at the point (row-scaled iff equilibrated)
    std::vector<double> rw, cw;   //!< equilibration weights (all 1 unless equilibrated)
    bool   equilibrated = false;
    bool   ic_override  = false;  //!< marching window k>0 continuity-pin rows applied
    size_t window       = 0;      //!< marching window index the assembly corresponds to
    double ms_deriv = 0., ms_asm = 0.;
    // cached sparsity pattern (filled on first use)
    bool pattern_ready = false;
    std::vector<size_t> rownnz;
    std::vector<std::vector<size_t>> rowcol;
    std::vector<double> grad;
    std::vector<double*> gradptr;
    // rev137: the drop mask (rows = claim rows, cols = tau columns) applied at assembly;
    // empty in rev137, set by the plan in rev138.
    std::vector<size_t> mask_rows, mask_cols;
    // rev137: cached IC_WEAK SAT re-injection stencil for the marching override rows
    std::map<size_t,std::vector<std::pair<size_t,double>>> ic_sat;
    bool   ic_sat_ready = false; size_t ic_sat_window = 0; bool ic_sat_override = false;
    void clear() { *this = t_JacAssembly{}; }
  };

  //! @brief assembled interface blocks B and C.
  //!
  //! B = multiplier columns over receiver rows; C = continuity rows over state columns.
  //! Shared by interface_rank_symmetry() and by the trace-projection phase, which must see
  //! the SAME blocks: two copies of this assembly that can drift is the defect class this
  //! whole arc diagnosed (a claim made in one pass about an object built in another, with
  //! no verification between).
  struct InterfaceBlocks {
    std::map<t_FrozenInterfaceClaimKey,size_t>    tau_ndx;   //!< claim -> dense tau index
    std::vector<t_FrozenInterfaceClaimKey>        tau_key;   //!< dense tau index -> claim
    std::map<FrozenPhysicalResidualRowKey,size_t> brow_ndx;  //!< receiver row -> dense index
    std::vector<std::map<size_t,double>>          bcol;      //!< B, by tau: row index -> coeff
    std::map<size_t,std::map<size_t,double>>      crow;      //!< C, by row: state col -> coeff
    std::map<size_t,t_FrozenInterfaceClaimKey>    ckey;      //!< C row -> claim
    size_t n_tau_col_in_C = 0;   //!< continuity coefficients that referenced a multiplier column
    size_t ntau() const noexcept { return tau_ndx.size(); }
  };

  //! @brief Per-state duplicate-node spread at an element interface.
  //!
  //! Nodes that coincide geometrically but belong to different elements carry SEPARATE
  //! collocation unknowns.  Continuity means those unknowns agree; the spread hi-lo over
  //! each geometric coincidence group is therefore a DIRECT measurement of the interface
  //! conditions actually being satisfied by the solution, independent of any rank, residual
  //! or convergence report.
  struct DuplicateNodeSpread
  {
    double max_pair_spread   = 0.0;  //!< over groups of multiplicity >= 2 (faces)
    double max_corner_spread = 0.0;  //!< over groups of multiplicity >= 4 (corners)
    size_t max_multiplicity  = 0;    //!< largest coincidence group seen
    size_t n_group           = 0;    //!< number of coincidence groups with multiplicity >= 2
  };

  //! @brief CAUSALITY of a lifted claim (a transition): its receivers must lie AFTER the jump.  An ordinary continuity
  //! claim x(tau+) = x(tau-) is symmetric, so receivers on both sides are right; a lifted claim x(tau+) = w, with w defined
  //! by pre-jump data, is not -- a receiver in the element ENDING at tau puts post-jump information into a pre-jump row
  //! (a penalty under IC_WEAK, a multiplier term under IC_STRONG / IC_TRACE), and the solution before tau depended on the
  //! jump (PSA1: 2.9e-3; 2026-09-30).  True for an edge that belongs to a lifted claim and receives in element iel_lo.
  bool _lifted_prejump_receiver
    ( OCPlan::InterfaceReceiverEdge const& e )
    const
    {
      if( !_lifted_claim_aux( e.claim_key ) ) return false;
      size_t const evo = _evolution_dom_var.id().second;
      for( auto const& [d, el] : e.receiver_block_el ) if( d == evo ) return el == e.claim_key.iel_lo;
      return e.receiver_dom_id == evo && e.receiver_iel == e.claim_key.iel_lo;
    }

  //! @brief view stored options as the solver's type.  THROWS if the object is not one -- an equation
  //! that bypassed _make_equation_options() is exactly the silent failure class this refactor keeps meeting,
  //! so it must fail loudly rather than read a base object as if it had solver fields.
  static EqnOptions const& _sopt
    ( FFModel::EqnOptions const& o )
    {
      auto const* p = dynamic_cast<EqnOptions const*>( &o );
      if( !p ) throw std::logic_error( "OCFESLV::_sopt ** equation options lack the solver's fields: "
                                       "an equation was stored without going through _make_equation_options()" );
      return *p;
    }

protected:
  bool _suppress_tensor_interface_corner
    ( Options::InterfaceType type, bool exact_reduced_weak, int block_id,
      FFVar const& face_dom,
      std::vector<t_InterfaceTensorCornerDim> const& od,
      FFVar const* recv_state = nullptr ) const;
  //! @brief C1: drop every claim whose (state, direction) the C0 classification does not tag NATURAL.
  //! Such a claim is not an independent flux condition -- its continuity follows from a pointwise relation, or
  //! is not owed at all -- so the plan never sees it and nothing has to rescue it.  Exact modes only; a census
  //! line per state names the tag.  @return the number of claims suppressed.
  size_t _apply_c0_claim_filter
    ( std::map<t_InterfaceClaimKey,size_t>& claims )
    const;

  bool _load_frozen_interface_plan
    ( std::map<t_InterfaceClaimKey,size_t>& claims ) const;

  //! @brief Build sparse coefficients for one exact trace claim.
  bool _trace_claim_coeffs
    ( t_InterfaceClaimKey const& claim,
      std::vector<std::pair<size_t,double>>& coeffs ) const;

  //! @brief the two VALUE functionals of a claim, one per side (lo = element
  //! iel_lo evaluated at +1, hi = iel_lo+1 at -1; C1 scaled by 2/width).  Same decode and
  //! threshold as _trace_claim_coeffs (which stays byte-identical: coeffs = lo - hi,
  //! interleaved); packed into PlanInput::DofEdge for the GLOBAL builder's trace DOFs.
  bool _trace_claim_sides
    ( t_InterfaceClaimKey const& claim,
      std::vector<std::pair<size_t,double>>& lo,
      std::vector<std::pair<size_t,double>>& hi ) const;

  //! @brief COLUMN-LEVEL pinned pairs.  For every state, seam direction, seam and
  //! transverse node, FIX1's predicate ("an INITIAL/BOUNDARY equation with mask ALL in the
  //! seam direction, whose transverse masks contain this node, pinning only this state")
  //! evaluated independently of any row visit.  Each pinned coincident pair is appended as a
  //! nodal DofEdge (identification/forest joins only, never realised).  MEASURED need: blk0's
  //! u at t=0 (no row visits the node) and PDE5's two BC-pinned y-claims at x=0/1 that the
  //! corner dedup claimed before the pin test could record them.
  void _pinned_column_pairs
    ( std::vector<PlanInput::DofEdge>& out ) const;

  //! @brief Build pre-compiled sparse SAT coefficients for one weak-SAT term.
  //!
  //! Each (col, coeff) entry satisfies coeff = prefactor * w_lagrange[k],
  //! so eval time reduces to buf[flat] += sum_k coeff_k * var[col_k].
  //! Returns true with empty coeffs when prefactor == 0 (valid, zero contribution).
  //! Returns false when the claim geometry is invalid.
  bool _sat_claim_coeffs
    ( t_InterfaceClaimKey const& claim,
      double prefactor,
      std::vector<std::pair<size_t,double>>& coeffs ) const;

  //! @brief Build executable frozen exact trace rows from retained claims.
  bool _build_frozen_trace_constraint_terms
    ( std::map<t_InterfaceClaimKey,size_t> const& claims,
      InterfacePlanDraft& draft ) const;

  //! @brief Shared redundancy pre-pass consumed by BOTH IC_TRACE and IC_STRONG.
  //!
  //! Builds per-column receiver coupling from the frozen weak-SAT terms (keyed
  //! by lambda_index, identical to the in-Schur column scheme but mode-agnostic)
  //! and flags every column-dependent continuity claim -- one whose continuity
  //! row is a linear combination of the others, hence Schur-uneliminable
  //! (IC_STRONG) and, kept explicit (IC_TRACE), a benign-but-singular null that
  //! stalls the solver.  Flagged claims are recorded in _strong_dropped_claims;
  //! returns false (with _strong_redundancy_pending set) to request the pass-2
  //! re-run that excludes them from allocation.  No-op for IC_WEAK and for
  //! column-rank-full plans (PDE2/3/4), so those builds are unchanged.
  bool _detect_redundant_continuity_claims
    ( InterfacePlanDraft& draft ) const;

  //! @brief Build setup-frozen tau-block audit/pivot/projection for IC_STRONG.
  bool _build_frozen_strong_tau_elimination
    ( InterfacePlanDraft& draft ) const;

  //! @brief True when a specific LINK equation defines an auxiliary in a given direction.
  bool _link_equation_in_direction
    ( FFVar const& eqn, int block_id, FFVar const& dom )
    const;

  //! @brief True when automatic same-domain interface handling is enabled in a block/direction.
  bool _block_receive_interface_in_direction
    ( int block_id, FFVar const& dom )
    const;

public:
  //! @brief Display frozen interface-plan diagnostics.
  //!
  //! This is intended for setup/debug regression checks of IC_WEAK/IC_TRACE/
  //! IC_STRONG interface-plan construction.  It does not change the plan and
  //! may be called after setup() on any imposition mode.
  bool display_interface_plan_diagnostics
    ( std::ostream& os,
      InterfacePlanGeometricExpectation const expectation = InterfacePlanGeometricExpectation::ANY,
      bool const print_details = true,
      size_t const max_duplicate_details = 20,
      size_t const max_geometric_details = 40 ) const
  {
    auto const& edges       = _plan.receiver_edges();
    auto const& weak_terms  = _plan.weak_sat_terms();
    auto const& trace_terms = _plan.trace_constraints();
    auto const& tau_claims  = _plan.tau_claims_frozen();

    InterfacePlanDiagnosticSummary const s = interface_plan_diagnostic_summary();

    bool const geo_zero_ok    = ( s.n_receiver_geometric == 0 && s.n_frozen_geometric == 0 );
    bool const geo_nonzero_ok = ( s.n_receiver_geometric >  0 && s.n_frozen_geometric >  0 );
    bool geo_expect_ok = true;
    if( expectation == InterfacePlanGeometricExpectation::ZERO )
      geo_expect_ok = geo_zero_ok;
    else if( expectation == InterfacePlanGeometricExpectation::NONZERO )
      geo_expect_ok = geo_nonzero_ok;

    os << "\n========== interface-plan diagnostics ==========\n";
    int const requested_imposition = _plan.frozen_options().imposition_type;
    int const actual_imposition = requested_imposition; // No automatic fallback is currently applied.
    int const current_imposition = static_cast<int>( options.INTERFACE.IMPOSITION );
    os << "imposition: requested=" << interface_imposition_name( requested_imposition )
       << " actual=" << interface_imposition_name( actual_imposition );
    if( current_imposition != requested_imposition )
      os << " current_option=" << interface_imposition_name( current_imposition );
    os << "\n";
    os << "receiver_edges=" << s.n_receiver_edges
       << " weak_sat_terms=" << s.n_weak_sat_terms
       << " trace_constraints=" << s.n_trace_constraints
       << " tau_claim_keys=" << s.n_tau_claim_keys << "\n";
    os << "trace/tau accounting:"
       << " trace_variables=" << s.n_trace_variables
       << " declared_constraints=" << s.n_trace_constraint_declared
       << " tau_columns=" << s.n_tau_claim_columns
       << " exact_replacements=" << s.n_exact_replacement_rows
       << " dropped_physical_rows=" << s.n_physical_rows_dropped_by_trace
       << " constraints_without_tau=" << s.n_constraints_without_tau
       << " tau_keys_without_receiver=" << s.n_tau_keys_without_receiver
       << "\n";
    os << "receiver edge sources:"
       << " principal=" << s.n_receiver_principal
       << " auxiliary=" << s.n_receiver_auxiliary
       << " geometric=" << s.n_receiver_geometric << "\n";
    os << "frozen weak/tau term sources:"
       << " principal=" << s.n_frozen_principal
       << " auxiliary=" << s.n_frozen_auxiliary
       << " geometric=" << s.n_frozen_geometric << "\n";
    os << "geometric-C0 candidates:"
       << " collected=" << s.n_geometric_candidates
       << " discarded_covered=" << s.n_geometric_discarded_covered
       << " emitted=" << s.n_geometric_emitted << "\n";

    if( !_ifaceEdgeAudit.empty() ){
      os << "per-state SAT receiver-edge audit  (edges | sat==0 | ->AUX | ->GEO | LEFT-ZERO"
            " | max|sat|)   rescue fallback coupling = "
         << OCPlan::kRescueFallbackCoupling << "  [constant since rev292]\n";
      for( auto const& kv : _ifaceEdgeAudit ){
        FFVar nm;
        for( auto const& vv : _mVar )
          if( static_cast<size_t>( vv.first.id().second ) == kv.first ){ nm = vv.first; break; }
        auto const& a = kv.second;
        double mx = 0.; auto itm = _ifaceEdgeMaxSat.find( kv.first );
        if( itm != _ifaceEdgeMaxSat.end() ) mx = itm->second;
        os << "  ";
        if( nm.dag() ) os << nm; else os << "id" << kv.first;
        os << " : " << a[0] << " | " << a[1] << " | " << a[2] << " | " << a[3]
           << " | " << a[4] << " | " << std::scientific << mx << std::defaultfloat;
        if( a[4] ) os << "   <== NULL ROWS: coupling 0 and no rescue";
        os << "\n";
      }
    }

    if( expectation == InterfacePlanGeometricExpectation::ZERO ){
      os << "CHECK n_geometric_receiver_edges == 0: "
         << ( s.n_receiver_geometric == 0 ? "PASS" : "FAIL" )
         << " value=" << s.n_receiver_geometric << "\n";
      os << "CHECK n_geometric_frozen_terms == 0: "
         << ( s.n_frozen_geometric == 0 ? "PASS" : "FAIL" )
         << " value=" << s.n_frozen_geometric << "\n";
    }
    else if( expectation == InterfacePlanGeometricExpectation::NONZERO ){
      os << "CHECK n_geometric_receiver_edges > 0: "
         << ( s.n_receiver_geometric > 0 ? "PASS" : "FAIL" )
         << " value=" << s.n_receiver_geometric << "\n";
      os << "CHECK n_geometric_frozen_terms > 0: "
         << ( s.n_frozen_geometric > 0 ? "PASS" : "FAIL" )
         << " value=" << s.n_frozen_geometric << "\n";
    }
    else{
      os << "n_geometric_receiver_edges=" << s.n_receiver_geometric << "\n";
      os << "n_geometric_frozen_terms=" << s.n_frozen_geometric << "\n";
    }

    os << "CHECK duplicate retained trace claims == 0: "
       << ( s.n_duplicate_retained_trace_claims == 0 ? "PASS" : "FAIL" )
       << " value=" << s.n_duplicate_retained_trace_claims << "\n";
    os << "CHECK duplicate tau claim map keys == 0: "
       << ( s.n_duplicate_tau_claim_keys == 0 ? "PASS" : "FAIL" )
       << " value=" << s.n_duplicate_tau_claim_keys << "\n";
    os << "CHECK duplicate receiver-edge identities == 0: "
       << ( s.n_duplicate_receiver_edge_identities == 0 ? "PASS" : "FAIL" )
       << " value=" << s.n_duplicate_receiver_edge_identities << "\n";

    bool trace_purity_ok = true;
    if( requested_imposition == Options::IC_TRACE ){
      bool const no_constraints_without_tau = ( s.n_constraints_without_tau == 0 );
      bool const no_exact_replacements      = ( s.n_exact_replacement_rows == 0 );
      bool const no_dropped_physical_rows   = ( s.n_physical_rows_dropped_by_trace == 0 );
      bool const no_tau_without_receiver    = ( s.n_tau_keys_without_receiver == 0 );
      bool const tau_columns_match_trace_variables = ( s.n_tau_claim_columns == s.n_trace_variables );
      bool const trace_constraints_match_tau = ( s.n_trace_constraints == s.n_trace_variables );

      os << "IC_TRACE audit:\n";
      os << "CHECK IC_TRACE constraints_without_tau == 0: "
         << ( no_constraints_without_tau ? "PASS" : "FAIL" )
         << " value=" << s.n_constraints_without_tau << "\n";
      os << "CHECK IC_TRACE exact_replacement_rows == 0: "
         << ( no_exact_replacements ? "PASS" : "FAIL" )
         << " value=" << s.n_exact_replacement_rows << "\n";
      os << "CHECK IC_TRACE physical_rows_dropped_by_trace == 0: "
         << ( no_dropped_physical_rows ? "PASS" : "FAIL" )
         << " value=" << s.n_physical_rows_dropped_by_trace << "\n";
      os << "CHECK IC_TRACE tau_keys_without_receiver == 0: "
         << ( no_tau_without_receiver ? "PASS" : "FAIL" )
         << " value=" << s.n_tau_keys_without_receiver << "\n";
      os << "CHECK IC_TRACE tau_columns == trace_variables: "
         << ( tau_columns_match_trace_variables ? "PASS" : "FAIL" )
         << " value=" << s.n_tau_claim_columns
         << " expected=" << s.n_trace_variables << "\n";
      os << "CHECK IC_TRACE trace_constraints == trace_variables: "
         << ( trace_constraints_match_tau ? "PASS" : "FAIL" )
         << " value=" << s.n_trace_constraints
         << " expected=" << s.n_trace_variables << "\n";

      trace_purity_ok = no_constraints_without_tau
                     && no_exact_replacements
                     && no_dropped_physical_rows
                     && no_tau_without_receiver
                     && tau_columns_match_trace_variables
                     && trace_constraints_match_tau;
    }

    if( print_details ){
      auto state_name_by_id = [&]( size_t const sid ) -> std::string
      {
        for( auto const& kv : _mVar )
          if( static_cast<size_t>( kv.first.id().second ) == sid ){
            std::ostringstream osn; osn << kv.first; return osn.str();
          }
        std::ostringstream osn; osn << "state#" << sid; return osn.str();
      };
      auto domain_name_by_id = [&]( size_t const did ) -> std::string
      {
        for( auto const& kv : _mDom )
          if( static_cast<size_t>( kv.first.id().second ) == did ){
            std::ostringstream osn; osn << kv.first; return osn.str();
          }
        std::ostringstream osn; osn << "dom#" << did; return osn.str();
      };
      auto claim_string = [&]( OCPlan::t_InterfaceClaimKey const& k ) -> std::string
      {
        std::ostringstream cs;
        cs << "{block=" << k.block_id
           << ",state=" << k.state_id
           << ",eqn=";
        if( k.eqn_id == std::numeric_limits<size_t>::max() ) cs << "*";
        else cs << k.eqn_id;
        cs << ",dom=" << k.dom_id
           << ",iel_lo=" << k.iel_lo
           << ",face=" << k.face
           << ",kind=" << k.kind << "}";
        return cs.str();
      };

      auto block_el_string = [&]( std::vector<std::pair<size_t,size_t>> const& bel ) -> std::string
      {
        std::ostringstream bs;
        bs << "[";
        for( size_t i = 0; i < bel.size(); ++i ){
          if( i ) bs << ",";
          bs << "(" << bel[i].first << ":" << bel[i].second << ")";
        }
        bs << "]";
        return bs.str();
      };

      std::map<std::pair<size_t,size_t>,size_t> geo_edges_by_state_dom;
      std::map<std::pair<size_t,size_t>,size_t> geo_terms_by_state_dom;
      std::map<OCPlan::t_InterfaceClaimKey,size_t> retained_trace_count;
      for( auto const& tr : trace_terms ) ++retained_trace_count[ tr.edge ];
      std::map<
        std::tuple< size_t, std::vector<std::pair<size_t,size_t>>, size_t,
                    OCPlan::t_FrozenInterfaceClaimKey, int, int >,
        size_t > receiver_edge_identity_count;
      for( auto const& e : edges ){
        auto const id = std::make_tuple(
          e.row_id, e.receiver_block_el, e.receiver_emit_flat,
          OCPlan::freeze_claim_key( e.claim_key ),
          static_cast<int>( e.source ), static_cast<int>( e.sat_type ) );
        ++receiver_edge_identity_count[ id ];
      }

      for( auto const& e : edges )
        if( e.source == OCPlan::InterfaceEdgeSource::GEOMETRIC_C0 )
          ++geo_edges_by_state_dom[ { e.claim_key.state_id, e.receiver_dom_id } ];
      for( auto const& t : weak_terms )
        if( t.source == OCPlan::InterfaceEdgeSource::GEOMETRIC_C0 )
          ++geo_terms_by_state_dom[ { t.claim_key.state_id, t.receiver_dom_id } ];

      if( !geo_edges_by_state_dom.empty() ){
        os << "geometric-C0 receiver edges by state/domain:\n";
        for( auto const& kv : geo_edges_by_state_dom )
          os << "  state=" << state_name_by_id( kv.first.first )
             << " direction=" << domain_name_by_id( kv.first.second )
             << " edges=" << kv.second << "\n";
      }
      if( !geo_terms_by_state_dom.empty() ){
        os << "geometric-C0 frozen terms by state/domain:\n";
        for( auto const& kv : geo_terms_by_state_dom )
          os << "  state=" << state_name_by_id( kv.first.first )
             << " direction=" << domain_name_by_id( kv.first.second )
             << " terms=" << kv.second << "\n";
      }

      size_t shown_geo = 0;
      for( auto const& e : edges ){
        if( e.source != OCPlan::InterfaceEdgeSource::GEOMETRIC_C0 ) continue;
        if( shown_geo++ == 0 ) os << "sample geometric-C0 receiver edges:\n";
        if( shown_geo > max_geometric_details ) break;
        os << "  row=" << e.row_id
           << " block_el=" << block_el_string( e.receiver_block_el )
           << " emit_flat=" << e.receiver_emit_flat
           << " state=" << state_name_by_id( e.claim_key.state_id )
           << " direction=" << domain_name_by_id( e.receiver_dom_id )
           << " claim=" << claim_string( e.claim_key )
           << " source=" << interface_edge_source_name( e.source )
           << " sat=" << static_cast<int>( e.sat_type )
           << " coupling=" << e.coupling
           << " trace=" << ( e.trace_tau_allowed ? 1 : 0 )
           << " weak=" << ( e.weak_sat_allowed ? 1 : 0 ) << "\n";
      }

      size_t shown_dup = 0;
      for( auto const& kv : retained_trace_count ){
        if( kv.second <= 1 ) continue;
        if( shown_dup++ == 0 ) os << "duplicate retained trace claims:\n";
        if( shown_dup > max_duplicate_details ) break;
        os << "  count=" << kv.second << " claim=" << claim_string( kv.first ) << "\n";
      }

      size_t shown_edge_dup = 0;
      for( auto const& kv : receiver_edge_identity_count ){
        if( kv.second <= 1 ) continue;
        if( shown_edge_dup++ == 0 ) os << "duplicate receiver-edge identities:\n";
        if( shown_edge_dup > max_duplicate_details ) break;
        auto const& id = kv.first;
        OCPlan::t_InterfaceClaimKey const claim = OCPlan::unfreeze_claim_key( std::get<3>( id ) );
        os << "  count=" << kv.second
           << " row=" << std::get<0>( id )
           << " block_el=" << block_el_string( std::get<1>( id ) )
           << " emit_flat=" << std::get<2>( id )
           << " claim=" << claim_string( claim )
           << " source=" << interface_edge_source_name( static_cast<OCPlan::InterfaceEdgeSource>( std::get<4>( id ) ) )
           << " sat=" << std::get<5>( id ) << "\n";
      }

      if( requested_imposition == Options::IC_TRACE && s.n_constraints_without_tau ){
        size_t shown_missing_tau = 0;
        for( auto const& tr : trace_terms ){
          auto const fk = OCPlan::freeze_claim_key( tr.edge );
          auto const it = tau_claims.find( fk );
          if( it != tau_claims.end() && !it->second.empty() ) continue;
          if( shown_missing_tau++ == 0 ) os << "trace constraints without tau:\n";
          if( shown_missing_tau > max_duplicate_details ) break;
          os << "  claim=" << claim_string( tr.edge ) << "\n";
        }
      }

      if( requested_imposition == Options::IC_TRACE && s.n_tau_keys_without_receiver ){
        size_t shown_tau_no_rx = 0;
        for( auto const& kv : tau_claims ){
          bool has_receiver = false;
          for( auto const& t : weak_terms ){
            if( !t.trace_tau_allowed ) continue;
            if( OCPlan::freeze_claim_key( t.claim_key ) == kv.first ){
              has_receiver = true;
              break;
            }
          }
          if( has_receiver ) continue;
          if( shown_tau_no_rx++ == 0 ) os << "tau claim keys without receiver:\n";
          if( shown_tau_no_rx > max_duplicate_details ) break;
          OCPlan::t_InterfaceClaimKey const claim = OCPlan::unfreeze_claim_key( kv.first );
          os << "  claim=" << claim_string( claim )
             << " tau_columns=" << kv.second.size() << "\n";
        }
      }
    }

    os << "==============================================\n";
    return geo_expect_ok
        && s.n_duplicate_retained_trace_claims == 0
        && s.n_duplicate_tau_claim_keys == 0
        && s.n_duplicate_receiver_edge_identities == 0
        && trace_purity_ok;
  }

protected:
  //! @brief null basis of sparse M at tolerance tol_abs via SPQR dead-column
  //! back-substitution + deflated inverse iteration for refinement demotions.
  static bool _null_basis( arma::sp_mat const& M, double tol_abs, size_t& rank_out, arma::mat& N );

  //! @brief Resolve the effective interface type for one equation/block/domain.
  Options::InterfaceType _resolve_interface_type
    ( FFModel::EqnOptions const* eqnopt, FFVar const* pdom, bool weak, int block_id )
    const;

  //! @brief Validate user-supplied incoming BCs at hyperbolic boundary faces:
  //! the count must equal rank of the incoming projector (A+ at LB / A- at UB)
  //! and each BC's state coupling must lie in the incoming characteristic space,
  //! rejecting outgoing-direction / mis-placed BCs (e.g. an inflow condition put
  //! on the wrong end).  Returns false on violation so setup() fails cleanly
  //! instead of producing a rank-deficient solve.  Always on (unconditional).
  bool _validate_hyperbolic_incoming_bcs();

  //! @brief Interface-symmetric SAT penalty parameter for the interface
  //! between elements iel_lo and iel_lo+1.  Returns (n_node-1)^2 divided by
  //! the interface-symmetric width, so the value is identical whether the
  //! interface is reached from its low- or high-side element.  Equal-and-
  //! opposite SAT contributions from the two elements then telescope exactly.
  double _sat_tau_face
    ( FFDom const& dom, size_t const iel_lo )
    const
    { double const nm1 = static_cast<double>( dom.n_node - 1 );
      double const w   = _sat_iface_width( dom, iel_lo );
      return ( w > 0. ) ? ( nm1 * nm1 / w ) : 0.; }

  //! @brief Principal-symbol SAT coupling for (symbolic equation ndx, state
  //! _vVar[state_pos], domain dom).  Returns the cached classified coupling
  //! when available, or 1.0 (identity) for unclassified models / uncoupled
  //! rows.  O(small) linear scan; cheaper than the prior inline triple loop.
  double _sat_coupling
    ( size_t const ndx, size_t const state_pos, FFVar const& dom )
    const;

  //! @brief Post-classification suitability + consistent-IC validator for the
  //! evolution/marching domain.  Runs at end of setup() (needs classification +
  //! reduction + DOF-audit results).  Standalone setup-hardening and the future
  //! marching/MOL applicability gate.  Warn-by-default; fatal only under
  //! Options::FATAL.EVOLUTION_DOMAIN.  Returns false only on a fatal violation.
  bool _validate_evolution_domain
    () const;

public:
  t_FrozenInterfaceClaimKey _trace_slot_key( t_InterfaceClaimKey key ) const
    {
      key.block_id = _trace_claim_block( key.state_id, key.dom_id );
      key.eqn_id   = std::numeric_limits<size_t>::max();
      if( key.kind == CLAIM_EXACT_ROW ) key.kind = CLAIM_EXACT_C0;
      return _freeze_interface_claim_key( key );
    }

  //! @brief Assemble B and C.  Takes its inputs explicitly so it serves the committed
  //!        _plan (audit, end of setup) and the draft (projection phase, mid-build) alike.
  void build_interface_blocks
    ( std::map<t_FrozenInterfaceClaimKey,std::vector<size_t>> const& tau_claims,
      std::vector<FrozenWeakSatTerm>         const& weak_sat_terms,
      std::vector<FrozenTraceConstraintTerm> const& trace_constraints,
      size_t trace_var_offset,
      InterfaceBlocks& out ) const
    {
      out = InterfaceBlocks{};

      // Dense tau index.  Covers BOTH the Schur-eliminated and the IC_STRONG keep-explicit
      // taus, so the count matches n_trace_var in either mode.
      for( auto const& kv : tau_claims ){
        if( kv.second.empty() ) continue;
        size_t const n = out.tau_ndx.size();
        out.tau_ndx.emplace( kv.first, n );
      }
      out.tau_key.assign( out.tau_ndx.size(), t_FrozenInterfaceClaimKey{} );
      for( auto const& kv : out.tau_ndx ) out.tau_key[ kv.second ] = kv.first;

      // B: multiplier columns over receiver rows.  Row indices are assigned in
      // weak_sat_terms order -- do not reorder, several reports quote these indices.
      out.bcol.assign( out.tau_ndx.size(), std::map<size_t,double>{} );
      for( auto const& term : weak_sat_terms ){
        if( !term.trace_tau_allowed ) continue;
        if( std::fabs( term.trace_prefactor ) <= 64. * DBL_EPSILON ) continue;
        auto const it = out.tau_ndx.find( _trace_slot_key( term.claim_key ) );
        if( it == out.tau_ndx.end() ) continue;
        FrozenPhysicalResidualRowKey rk;
        rk.row_id   = term.row_id;
        rk.block_el = term.receiver_block_el;
        rk.flat     = term.receiver_emit_flat;
        auto ir = out.brow_ndx.find( rk );
        if( ir == out.brow_ndx.end() ){
          size_t const n = out.brow_ndx.size();
          ir = out.brow_ndx.emplace( rk, n ).first;
        }
        out.bcol[ it->second ][ ir->second ] += term.trace_prefactor;
      }

      // C: continuity rows over state columns only.
      for( auto const& t : trace_constraints ){
        auto& r = out.crow[ t.row ];
        for( auto const& cv : t.coeff ){
          if( cv.first >= trace_var_offset ){ ++out.n_tau_col_in_C; continue; }
          r[ cv.first ] += cv.second;
        }
        out.ckey[ t.row ] = _trace_slot_key( t.edge );
      }
    }

  // =====================================================================
  // rev150c [statenull]: census of the STATE-SIDE null directions of the audited
  // Jacobian, i.e. the def_right - k' directions that are not pure-tau.
  //
  // WHY.  MBC3 under IC_WEAK is STATE_CONTENT (defR=8, growing 8->11->17 with
  // refinement) with no taus in the system; under IC_TRACE the same model reports
  // defR=k' -- every null direction pure-tau, the 8 state directions CLOSED.  TRACE
  // only adds tau columns, C rows and exact-replacement rows to the WEAK plan, so
  // one of those supplies rank the weak realisation lacks.  The record attributed
  // the 8 to mesh/boundary DOFs (Defect-A pattern) from the growth law alone;
  // nobody had looked at the vectors.  This census does.
  //
  // WHAT IT PRINTS.  The rev140 QR dead-column basis of the full J (orthonormal,
  // original column coordinates, un-equilibrated through cw).  The PER-DIRECTION
  // lines are basis-dependent (any rotation of the null space is another basis) and
  // exist for localisation only.  The MASS AGGREGATES are not: for an orthonormal
  // basis, sum_d v[c,d]^2 is the diagonal of the null-space projector, a property of
  // the space, and the state-block share of it is the number of state-side null
  // directions.  Read the aggregates first.
  //
  // NODE CLASS per column, per dimension of the state's domain:
  //   seam      first node of a non-first element, or last node of a non-last one
  //   boundary  first node of the first element, or last node of the last one
  //   interior  anything else
  // A column is 'seam' if any dimension is a seam, else 'boundary' if any is a
  // boundary, else 'interior'.  Mass on seam nodes of claimed states says
  // "interface realisation"; mass on boundary/interior says "not the interface plan".
  //
  // Gated on CRONOS_STATENULL_CENSUS=1; default off.  Runs one extra SPQR on the
  // audited J when def_right > k'.  det-inert: prints only.
  // =====================================================================
  void _statenull_census( t_JacAssembly const& A, double tol_abs, size_t k_prime ) const
    {
      size_t rk = 0; arma::mat NB;
      if( !_null_basis( A.J, tol_abs, rk, NB ) || !NB.n_cols ){
        _disp( 3 ) << "  [statenull] null basis unavailable (no SPQR, or empty)" << std::endl;
        return;
      }
      size_t const m    = (size_t)A.J.n_cols;
      size_t const tau0 = A.n_tau ? A.tau0 : m;          // state columns are [0,tau0)
      size_t const K    = (size_t)NB.n_cols;

      // ---- decode every state column once --------------------------------------
      struct ColInfo { int state = -1; int cls = -1; std::string seam; };   // cls 0 int,1 bnd,2 seam
      std::vector<ColInfo> info( tau0 );
      std::vector<std::string> sname;
      { int sid = 0;
        for( auto const& [v, dom] : _mVar ){
          { std::ostringstream o; o << v; sname.push_back( o.str() ); }
          auto itb = _blkVar.find( v );
          if( itb == _blkVar.end() ){ ++sid; continue; }
          size_t const base = std::get<0>( itb->second );
          size_t const bsz  = std::get<1>( itb->second );
          size_t span = bsz;
          for( auto const& d : dom ){ auto itd = _mDom.find( d ); if( itd == _mDom.end() ){ span = 0; break; } span *= itd->second.n_elem; }
          for( size_t col = base; span && bsz && col < base + span && col < tau0; ++col ){
            size_t local = col - base, blkid = local / bsz, nodeoff = local % bsz;
            int cls = 0; std::string seam;
            for( auto const& d : dom ){
              auto const& dd = _mDom.at( d );
              size_t const ne = dd.n_elem ? dd.n_elem : 1, nn = dd.n_node ? dd.n_node : 1;
              size_t const iel = blkid % ne;  blkid /= ne;
              size_t const ind = nodeoff % nn; nodeoff /= nn;
              bool const lo = ( ind == 0 ), hi = ( ind + 1 == nn );
              if( ( lo && iel > 0 ) || ( hi && iel + 1 < ne ) ){
                cls = 2; if( !seam.empty() ) seam += "+"; seam += d.name() + ( lo ? "-lo" : "-hi" );
              }
              else if( ( lo && iel == 0 ) || ( hi && iel + 1 == ne ) ){ if( cls < 1 ) cls = 1; }
            }
            info[col].state = sid; info[col].cls = cls; info[col].seam = seam;
          }
          ++sid;
        } }

      // ---- per-direction state/tau split; un-equilibrate through cw --------------
      std::vector<double> smass( K, 0. );
      std::vector<arma::vec> vv( K );
      size_t n_state_side = 0;
      for( size_t d = 0; d < K; ++d ){
        arma::vec v = NB.col( (arma::uword)d );
        if( A.equilibrated && A.cw.size() == m )
          for( size_t c = 0; c < m; ++c ) v( (arma::uword)c ) *= A.cw[c];
        double const nrm = arma::norm( v ); if( nrm > 0. ) v /= nrm;
        double s = 0.; for( size_t c = 0; c < tau0; ++c ) s += v( (arma::uword)c ) * v( (arma::uword)c );
        smass[d] = s; vv[d] = std::move( v );
        if( s > 0.5 ) ++n_state_side;
      }

      std::ios::fmtflags const cf = std::cerr.flags(); auto const cp = std::cerr.precision();
      _disp( 3 ) << "  [statenull] J " << A.J.n_rows << "x" << m << " rank=" << rk
                << " null_dim=" << K << " state-side=" << n_state_side
                << " (expected def_right-k'=" << ( K > k_prime ? K - k_prime : 0 ) << ")"
                << " tau0=" << tau0 << " equilibrated=" << ( A.equilibrated ? 1 : 0 )
                << "  [aggregates are basis-independent; dir# lines are not]" << std::endl;

      // ---- aggregates over the state-side directions ---------------------------
      std::map<int,double> by_state; std::map<int,double> by_cls; std::map<std::string,double> by_seam;
      double tot = 0.;
      for( size_t d = 0; d < K; ++d ){
        if( smass[d] <= 0.5 ) continue;
        for( size_t c = 0; c < tau0; ++c ){
          double const w = vv[d]( (arma::uword)c ) * vv[d]( (arma::uword)c );
          if( w == 0. ) continue;
          tot += w;
          by_state[ info[c].state ] += w; by_cls[ info[c].cls ] += w;
          if( info[c].cls == 2 ) by_seam[ info[c].seam ] += w;
        }
      }
      _disp( 3 ) << std::fixed << std::setprecision(3);
      _disp( 3 ) << "  [statenull] state-side mass total=" << tot << "  by state:";
      { std::vector<std::pair<double,int>> srt; for( auto const& kv : by_state ) srt.push_back( { kv.second, kv.first } );
        std::sort( srt.rbegin(), srt.rend() );
        for( auto const& [w, s] : srt ) _disp( 3 ) << "  " << ( s >= 0 && (size_t)s < sname.size() ? sname[s] : "?" ) << "=" << w; }
      _disp( 3 ) << std::endl;
      _disp( 3 ) << "  [statenull] by node class:  seam=" << by_cls[2] << "  boundary=" << by_cls[1]
                << "  interior=" << by_cls[0] << std::endl;
      if( !by_seam.empty() ){
        _disp( 3 ) << "  [statenull] seam mass by direction:";
        for( auto const& kv : by_seam ) _disp( 3 ) << "  " << kv.first << "=" << kv.second;
        _disp( 3 ) << std::endl;
      }

      // ---- per-direction detail (localisation only) ----------------------------
      size_t shown = 0;
      for( size_t d = 0; d < K && shown < 24; ++d ){
        if( smass[d] <= 0.5 ) continue;
        ++shown;
        arma::vec const av = arma::abs( vv[d] );
        arma::uvec const idx = arma::sort_index( av, "descend" );
        _disp( 3 ) << "  [statenull] dir#" << d << " state_mass=" << smass[d] << ":";
        for( size_t j = 0; j < idx.n_elem && j < 6; ++j ){
          size_t const c = (size_t)idx(j);
          if( av( idx(j) ) < 1e-6 ) break;
          _disp( 3 ) << "  " << std::setprecision(3) << vv[d]( idx(j) ) << "@"
                    << ( c < tau0 ? decode_state_column( c ) : ( "tau" + std::to_string( c - tau0 ) ) )
                    << ( c < tau0 && info[c].cls == 2 ? "[" + info[c].seam + "]" : "" );
        }
        _disp( 3 ) << std::endl;
      }
      if( n_state_side > shown )
        _disp( 3 ) << "  [statenull] ... (" << n_state_side << " state-side directions total)" << std::endl;

      // ---- rev150d [statenull/rows]: who sees the direction, and who cancels whom ---
      // The net J*v is zero by definition.  The GROSS row mass sum_c |J_rc v_c| is
      // nonzero exactly on the rows that see the direction; within such a row the
      // signed contributions J_rc*v_c, grouped by column state, name the cancellation.
      // Computed on the EQUILIBRATED J and its null vector (columns folded through cw,
      // rows through rw), so row ranking is by scaled row -- the right scale for "does
      // this row see it" since that is the scale the rank decision was taken at.
      // Row provenance: _rowEqnNdx[r] = source equation index (row_id in taucol keys),
      // -1 for plan-owned rows; only used when its length matches J.
      // CRONOS_STATENULL_ROWS=1; default off.
      { static bool const kRows = Options::_env_flag( "CRONOS_STATENULL_ROWS", false );
        if( kRows ){
          bool const have_prov = ( _rowProvenanceReady && _rowEqnNdx.size() == (size_t)A.J.n_rows );
          size_t const nr = (size_t)A.J.n_rows;
          size_t rshown = 0;
          for( size_t d = 0; d < K && rshown < 8; ++d ){
            if( smass[d] <= 0.5 ) continue;
            ++rshown;
            arma::vec const vs = NB.col( (arma::uword)d );      // equilibrated coordinates
            std::vector<double> gross( nr, 0. ), net( nr, 0. );
            for( arma::sp_mat::const_iterator it = A.J.begin(); it != A.J.end(); ++it ){
              double const t = (*it) * vs( it.col() );
              if( t == 0. ) continue;
              gross[ it.row() ] += std::fabs( t ); net[ it.row() ] += t;
            }
            std::map<int,double> by_rowid; double gtot = 0.;
            std::vector<size_t> order; order.reserve( nr );
            for( size_t r = 0; r < nr; ++r ){
              if( gross[r] <= 0. ) continue;
              gtot += gross[r]; order.push_back( r );
              by_rowid[ have_prov ? _rowEqnNdx[r] : -2 ] += gross[r];
            }
            std::sort( order.begin(), order.end(), [&]( size_t a, size_t b ){ return gross[a] > gross[b]; } );
            _disp( 3 ) << "  [statenull/rows] dir#" << d << " rows-seeing=" << order.size()
                      << " gross-total=" << gtot << "  gross by row_id:";
            { std::vector<std::pair<double,int>> srt; for( auto const& kv : by_rowid ) srt.push_back( { kv.second, kv.first } );
              std::sort( srt.rbegin(), srt.rend() ); size_t n = 0;
              for( auto const& [w, id] : srt ){ if( ++n > 6 ) break; _disp( 3 ) << "  row_id=" << id << ":" << w; } }
            _disp( 3 ) << std::endl;
            for( size_t q = 0; q < order.size() && q < 6; ++q ){
              size_t const r = order[q];
              // signed contributions of this row, by column
              std::vector<std::pair<double,size_t>> contrib;
              std::map<int,double> by_state;
              for( arma::sp_mat::const_row_iterator it = A.J.begin_row( r ); it != A.J.end_row( r ); ++it ){
                double const t = (*it) * vs( it.col() );
                if( t == 0. ) continue;
                contrib.push_back( { t, (size_t)it.col() } );
                by_state[ it.col() < tau0 ? info[ it.col() ].state : -3 ] += t;
              }
              std::sort( contrib.begin(), contrib.end(),
                         []( auto const& a, auto const& b ){ return std::fabs( a.first ) > std::fabs( b.first ); } );
              _disp( 3 ) << "      r=" << r << " row_id=" << ( have_prov ? _rowEqnNdx[r] : -2 )
                        << " gross=" << gross[r] << " net=" << std::scientific << std::setprecision(1) << net[r]
                        << std::fixed << std::setprecision(3) << "  by state:";
              for( auto const& kv : by_state )
                _disp( 3 ) << "  " << ( kv.first >= 0 && (size_t)kv.first < sname.size() ? sname[kv.first] : ( kv.first == -3 ? "tau" : "?" ) )
                          << "=" << std::showpos << kv.second << std::noshowpos;
              _disp( 3 ) << std::endl << "          ";
              for( size_t j = 0; j < contrib.size() && j < 6; ++j ){
                size_t const c = contrib[j].second;
                _disp( 3 ) << "  " << std::showpos << contrib[j].first << std::noshowpos << "@"
                          << ( c < tau0 ? decode_state_column( c ) : ( "tau" + std::to_string( c - tau0 ) ) )
                          << ( c < tau0 && info[c].cls == 2 ? "[" + info[c].seam + "]" : "" );
              }
              _disp( 3 ) << std::endl;
            }
          }
        } }
      std::cerr.flags( cf ); std::cerr.precision( cp );
    }

  //! @brief Duplicate-node spread for one state, over the flat collocation vector `var`.
  //!
  //! `off` is the state's base offset into `var`.  OMIT IT: the default derives it from
  //! pos_state(), which is the only correct source.  Computing it by walking states and
  //! accumulating node_colloc().size() is WRONG in the presence of scalar states -- see the
  //! note in report_duplicate_node_spreads().  The parameter is retained only for
  //! callers that already hold a verified offset.
  DuplicateNodeSpread duplicate_node_spread
    ( FFVar const& st, std::vector<double> const& var,
      size_t off = std::numeric_limits<size_t>::max() )
    const;

  //! @brief Report per-state duplicate-node spreads and return the worst max_pair.
  //!
  //! WHY THIS LIVES IN THE HEADER.  The withdrawn Stage A interface rule passed k=0,
  //! SYMMETRIC, a clean residual and clean primitive variables while destroying MBC5's
  //! radial fluxes by 24 orders of magnitude.  The ONLY instrument that caught it was the
  //! per-state duplicate-node spread -- and it was implemented driver-side, so only 3 of
  //! the 29 drivers flagged by the rank sweep could report it.  Moving it here gives every
  //! driver the gate for free and makes "no readout" impossible rather than merely unlikely.
  //!
  //! Returns the worst max_pair_spread across all states, so a caller can gate on one
  //! number:  ok &= ( oc.report_duplicate_node_spreads( x, "post solve" ) < tol );
  double report_duplicate_node_spreads
    ( std::vector<double> const& var, char const* label = "",
      std::ostream* os = &std::cerr )
    const;

protected:
  enum { CLAIM_EXACT_C0 = 0, CLAIM_EXACT_C1 = 1, CLAIM_EXACT_ROW = CLAIM_EXACT_C0 };

public:
  //! @brief result of numeric_rank().
  struct t_RankResult
  {
    enum Backend { NONE=0, SPQR, EIGEN_QR, DENSE_SVD };
    size_t  rank = 0;
    double  gap = 0.;         //!< SVD: sigma[rank-1]/sigma[rank] (or /tol when nothing cut); QR: min|R_jj|/tol_abs
    double  tol = 0.;         //!< the relative tolerance requested
    double  tol_abs = 0.;     //!< the absolute threshold actually applied: max(tol, 20(nr+nc)eps) * max column norm
    double  sigma_cut = 0.;   //!< smallest retained singular value (SVD) or |R_jj| (QR)
    double  ms = 0.;
    Backend backend = NONE;
  };

  //! @brief setup-time determinacy verdict for one audited window.
  struct t_Determinacy
  {
    enum Verdict { SKIPPED=0, DETERMINED, GAUGE_ONLY, STATE_CONTENT, INCONSISTENT_TOL };
    bool   ran = false;
    bool   marching = false;
    size_t window = 0;
    std::string skipped_reason, ref_point, mode, backend;
    size_t n = 0, m = 0, n_tau = 0, rank = 0;
    size_t def_right = 0;     //!< m - rank : free columns (determinacy)
    size_t def_left  = 0;     //!< n - rank : dependent rows (consistency)
    size_t k_prime   = 0;     //!< n_tau - rank([B;D]) : pure-tau null dimension
    size_t free_expected = 0; //!< m-n when m>n (decision variables, free by design)
    size_t dep_expected  = 0; //!< n-m when n>m (dependent rows by design)
    size_t dep_excess    = 0; //!< dependent rows beyond dep_expected (consistency risk)
    Verdict verdict = SKIPPED;
    double gap = 0., sigma_cut = 0., tol_abs = 0., ms_assemble = 0., ms_rank = 0.;
    bool   marginal_reeval = false;   //!< second evaluation at 1.3x the reference point was run
    size_t rank2 = 0; double gap2 = 0.; Verdict verdict2 = SKIPPED;
  };

protected:
  std::vector<t_Determinacy>                _determinacy;   //!< one entry per audited window

public:
  //! @brief setup-time determinacy audit results (one per audited window).
  std::vector<t_Determinacy> const& determinacy_audit() const { return _determinacy; }

  //! a gap at the cut below this triggers the re-evaluation at 1.3x the
  //! reference point (1e2).  A cut of 1e6 fires on nearly every
  //! DETERMINED instance in the corpus (sigma_min/tol ~ 1e3-1e5 is ordinary
  //! conditioning for a discretised PDE) and doubles the audit cost at n > 1e4
  //! for nothing; the genuinely marginal cases sit at 6-20.  Environment:
  //! CRONOS_DETERMINACY_GAP_MIN (an AUDIT parameter, kept as a diagnostic).
  static double _knob_DETERMINACY_GAP_MIN
    ()
    { return Options::_env_dbl( "CRONOS_DETERMINACY_GAP_MIN", 1e2 ); }

  //! largest dimension the dense-SVD fallback backend accepts (3000).
  //! (The CRONOS_DETERMINACY_MAX_DENSE variable is retired, 2026-10-07: a measured constant.)
  static size_t _knob_DETERMINACY_MAX_DENSE
    ()
    { return 3000; }

protected:
  //! @brief emit one stored [determinacy] line (verbatim the audit print;
  //! used by the audit itself and by the w-decide replay after a kept restore).
  void _print_determinacy_window( t_Determinacy const& D, size_t w ) const;

public:
  //! @brief numerical rank of a sparse matrix at ONE relative tolerance, through
  //! the configured backend (SPQR / Eigen SparseQR / dense SVD under the cap).
  //! tol_abs_override > 0 replaces the computed absolute threshold (the audit uses it
  //! so rank(M) and rank([B;D]) are cut at the SAME value).
  //! @brief One rank pass on an assembled Jacobian: rank, deficiencies, k' (the tau-block
  //! null dimension) and the DETERMINED / GAUGE_ONLY / STATE_CONTENT / INCONSISTENT_TOL
  //! verdict.  Knows nothing of what was assembled: OCFESLV hands it the space-time window
  //! Jacobian; OCFEMOL will hand it the pencil of the semi-discrete DAE at a point.
  bool _determinacy_verdict
    ( t_JacAssembly const& A, t_Determinacy& D )
    const;

  t_RankResult numeric_rank
    ( arma::sp_mat const& M, double reltol, double tol_abs_override=-1. )
    const;

protected:
  //! @brief does each kept-explicit tau column actually HAVE a row?
  //! The header has carried this suspicion -- "the promised explicit tau
  //! unknown PLUS a real continuity row appears to deliver the column but not the row"
  //! -- with no instrument to settle it.  MEASURED on PDE5 IC_STRONG mono, tau block
  //! 24036x36: nnzA=0, rank=0, 36 column singletons.  Empty, confirmed.  This routine
  //! reports the per-column occupancy alongside the _keepExplicitReason census, so the
  //! two halves of the question -- WHICH keep-explicit route, and WHETHER the row was
  //! delivered -- appear on one line.  CRONOS_KEEP_EXPLICIT_AUDIT=1, default off.
  void _keep_explicit_audit( arma::sp_mat const& T, size_t tau0, size_t rank_T ) const;

  size_t      _evoCoreRows         = 0;      //!< differential-core rows (user-side)
  size_t      _evoCoreCols         = 0;      //!< ... and differentiated states

  void _store_frozen_interface_plan
    ( std::map<t_InterfaceClaimKey,size_t> const& claims,
      InterfacePlanDraft& draft,
      PlanInput const& in ) const;   // rev152b: decision layer reads through the contract

  //! @brief pack everything the plan decision layer reads into a
  //! PlanInput.  Called before _store_frozen_interface_plan; the plan decision layer reads
  //! the result.
  PlanInput _pack_plan_input
    ( std::map<t_InterfaceClaimKey,size_t> const& claims,
      InterfacePlanDraft const& draft ) const;

  //! @brief Build the frozen weak-SAT consumer list from the receiver graph.
  bool _build_frozen_weak_sat_terms
    ( InterfacePlanDraft& draft ) const;

  //! @brief Whether an equation-node trace is retained after exact interface replacement.
  bool _keep_eqn_node
    ( FFVar const& eqn, FFModel::EqnOptions const& opt, FFVar const& dom, FFDom const& pdom,
      size_t iel, size_t inode, int lim, bool weak, size_t eqn_ndx = std::numeric_limits<size_t>::max() )
    const;

public:
  //! @brief [taucol/seam]: LOCALISATION LADDER -- claim count vs rank(B) at four
  //! grouping scales, each compared against k.
  //!
  //! WHY A LADDER.  Grouping by (direction, element, face) and comparing the summed
  //! deficit against k is not enough: on OCFE_PDE7 it reports deficit=0 while k=4 -- every group
  //! is full column rank and the dependency runs BETWEEN groups.  The arithmetic makes this
  //! inevitable rather than surprising --
  //!
  //!     sum_i ( |G_i| - rank(B|G_i) )  =  n - sum_i rank(B|G_i)  <=  n - rank(B)  =  k
  //!
  //! because rank(B) <= sum_i rank(B|G_i) for any partition.  So a partitioned deficit is
  //! always a LOWER BOUND on k, tight only when the groups' column spaces are independent.
  //! A single grouping can therefore under-report and never over-report, and reporting one
  //! number invites reading a 0 as "no deficiency" when it means "not at this scale".
  //!
  //! The ladder coarsens until the deficit reaches k.  The first level that accounts for k
  //! is the SCALE of the coupling, which is the localisation the census exists to give:
  //!   face      the claims at one transverse node of one interface
  //!   interface one (direction, element) seam, all its transverse nodes
  //!   direction all seams in one domain direction
  //!   global    every tau; equals k by construction, so it is a check on the arithmetic,
  //!             not evidence -- if it disagrees with PLAN METRIC something is wrong here.
  void _report_tau_seam_capacity( std::ostream& os ) const
    {
      InterfaceBlocks blk;
      build_interface_blocks( _plan.tau_claims_frozen(), _plan.weak_sat_terms(),
                              _plan.trace_constraints(), _plan.trace_var_offset(), blk );
      size_t const ntau = blk.ntau();
      if( !ntau ){ os << "  [taucol/seam] no tau columns.\n"; return; }

      // Per-tau coordinates, resolved once.
      struct TauAt { size_t dom, iel, face, state; };
      std::vector<TauAt> at( ntau );
      for( size_t j = 0; j < ntau; ++j ){
        t_InterfaceClaimKey const ck = _unfreeze_interface_claim_key( blk.tau_key[j] );
        at[j] = { ck.dom_id, ck.iel_lo, ck.face, ck.state_id };
      }

      // rank of B restricted to a set of tau columns.
      auto rank_of = [&]( std::vector<size_t> const& cols ) -> size_t {
        if( cols.empty() ) return 0;
        std::set<size_t> rowset;
        for( size_t j : cols )
          for( auto const& rv : blk.bcol[j] ) rowset.insert( rv.first );
        if( rowset.empty() ) return 0;
        std::vector<size_t> rows( rowset.begin(), rowset.end() );
        arma::mat Bs( (arma::uword)rows.size(), (arma::uword)cols.size(), arma::fill::zeros );
        for( arma::uword jj = 0; jj < (arma::uword)cols.size(); ++jj )
          for( auto const& rv : blk.bcol[ cols[jj] ] ){
            auto it = std::lower_bound( rows.begin(), rows.end(), rv.first );
            Bs( (arma::uword)( it - rows.begin() ), jj ) = rv.second;
          }
        return (size_t)arma::rank( Bs, 1.0e-10 * std::max( 1.0, arma::norm( Bs, 2 ) ) );
      };

      // level: 0 = face, 1 = interface, 2 = direction, 3 = global
      auto deficit_at = [&]( int level, std::map<std::string,std::pair<size_t,size_t>>* detail )
                        -> size_t {
        std::map<std::string,std::vector<size_t>> grp;
        for( size_t j = 0; j < ntau; ++j ){
          std::ostringstream key;
          if( level == 3 ) key << "all";
          else{
            key << _taucol_dom_name( at[j].dom );
            if( level <= 1 ) key << " iel_lo=" << at[j].iel;
            if( level == 0 ) key << " face=" << at[j].face;
          }
          grp[ key.str() ].push_back( j );
        }
        size_t tot = 0;
        for( auto const& kv : grp ){
          size_t const r = rank_of( kv.second );
          tot += ( kv.second.size() > r ) ? ( kv.second.size() - r ) : 0;
          if( detail ) (*detail)[ kv.first ] = { kv.second.size(), r };
        }
        return tot;
      };

      size_t const k_plan = ( _plan.n_trace_var() > _traceProjRank )
                          ? ( _plan.n_trace_var() - _traceProjRank ) : 0;

      std::map<std::string,std::pair<size_t,size_t>> det_face, det_iface;
      size_t const d_face  = deficit_at( 0, &det_face );
      size_t const d_iface = deficit_at( 1, &det_iface );
      size_t const d_dir   = deficit_at( 2, nullptr );
      size_t const d_all   = deficit_at( 3, nullptr );

      char const* scale = ( d_face  >= k_plan && k_plan ) ? "face"
                        : ( d_iface >= k_plan && k_plan ) ? "interface"
                        : ( d_dir   >= k_plan && k_plan ) ? "direction"
                        : ( d_all   >= k_plan && k_plan ) ? "global (spans directions)"
                        :                                   "none (k=0)";

      os << "  [taucol/seam] localisation ladder (summed deficit at each grouping):\n"
         << "      by (dir,elem,face) : " << d_face  << "\n"
         << "      by (dir,elem)      : " << d_iface << "\n"
         << "      by (dir)           : " << d_dir   << "\n"
         << "      global             : " << d_all   << "   [k=" << k_plan << "]\n"
         << "      => coupling scale  : " << scale << "\n";
      if( d_all != k_plan )
        os << "      ** global deficit != k.  These are computed by different paths"
              " (this census from B directly, k from the projection); a disagreement"
              " here is an arithmetic fault in one of them, not a model property.\n";

      // Detail at the finest level that shows anything, for localisation.
      auto const& det = ( d_face > 0 ) ? det_face : det_iface;
      char const* dlab = ( d_face > 0 ) ? "(dir,elem,face)" : "(dir,elem)";
      size_t shown = 0;
      for( auto const& kv : det ){
        if( kv.second.first <= kv.second.second ) continue;
        if( shown == 0 )
          os << "      groups with deficit>0, by " << dlab << ":\n";
        if( ++shown > 12 ){ os << "        ... (" << shown - 12 << " more)\n"; break; }
        os << "        " << kv.first << ":  claims=" << kv.second.first
           << "  rank=" << kv.second.second
           << "  deficit=" << ( kv.second.first - kv.second.second ) << "\n";
      }
      size_t rows_e = 0, cols_e = 0;
      (void)rows_e; (void)cols_e;
      os << "      TOTAL: claims=" << ntau << "  seam-deficit=" << d_face
         << "  interface-deficit=" << d_iface
         << "  global-deficit=" << d_all
         << "   [reference: differential core = ";
      if( _evoCoreRows || _evoCoreCols ) os << _evoCoreRows << "x" << _evoCoreCols;
      else                               os << "NOT MEASURED";
      os << ", n_trace_var=" << _plan.n_trace_var() << "]\n";
    }

  //! @brief Stage 1 (step 2): trace projection onto range(B^T) -- k=0 path only.
  //!
  //! Option D keeps every multiplier and every continuity row, and removes only the
  //!
  //! component along the implied combination:
  //!
  //!     [ A          B ]        P  = Z Z^T          projector onto null(B)
  //!     [ C_tilde    P ]        C_tilde = (I - P) C
  //!
  //! P and C_tilde are basis-independent, so there is no selection, no pivot order and no
  //! non-unique demoted set -- the property that distinguishes this from claim-dropping,
  //! which was implemented as Stage A and falsified.
  //!
  //! THIS REVISION IMPLEMENTS THE k = 0 PATH ONLY.  When null(B) is trivial the projector is
  //! the zero matrix and C_tilde == C, so the correct behaviour is to touch nothing at all.
  //! That covers 148 of 234 exact-mode corpus configurations including ALL of PSA, which is
  //! why it is worth landing on its own: it exercises the phase ordering and the shared
  //! assembly against the DRAFT (rather than the committed plan) with nothing at stake, and
  //! bit-identity on PSA is a gate that cannot be argued with.
  //!
  //! k > 0 is detected, reported and then left alone -- the projector itself is the next
  //! commit.
  bool _build_frozen_trace_projection
    ( OCPlan::InterfacePlanDraft& draft )
    {
      _traceProjK       = 0;
      _traceProjNTau    = 0;
      _traceProjApplied = false;
      _traceProjRank    = 0;
      if( options.INTERFACE.IMPOSITION != Options::IC_TRACE
       && options.INTERFACE.IMPOSITION != Options::IC_STRONG ) return true;
      if( !draft.weak_sat_terms_ready || !draft.trace_constraints_ready ) return true;

      // Same assembly the audit uses, run against the DRAFT.  If these two ever disagree,
      // the projection would be computed on a different B from the one measured.
      InterfaceBlocks blk;
      build_interface_blocks( draft.tau_claims, draft.weak_sat_terms,
                              draft.trace_constraints, draft.trace_var_offset, blk );
      size_t const m = blk.ntau();
      _traceProjNTau = m;
      if( !m ) return true;

      // rank(B) by the same sparse modified Gram-Schmidt, at the same 1e-9 relative cut, as
      // interface_rank_symmetry.  Measured headroom across the corpus: smallest RETAINED
      // singular value 1.6e-5*sigma_max, largest NULL value ~7e-17*sigma_max.
      std::vector<std::vector<std::pair<size_t,double>>> Bv( m );
      for( size_t j = 0; j < m; ++j )
        for( auto const& e : blk.bcol[j] )
          Bv[j].push_back( std::make_pair( e.first, e.second ) );
      size_t const rank = _sparse_column_rank( Bv );
      _traceProjK = m - rank;

      // rev40: record rank(B) UNCONDITIONALLY.  Until now _traceProjRank was written only on
      // the applied path (IC_TRACE, k>0), so it was 0 for every IC_STRONG plan and for every
      // k=0 plan -- and it is read nowhere in rev39, so this is a write-only field being
      // given its first consumer.  The Stage 2 pre-flight probe in
      // _build_frozen_strong_tau_elimination compares the pivot-search stall count against
      // it, which requires it to be live on exactly the path that currently fails.
      _traceProjRank = rank;

      if( !_traceProjK ){
        // k = 0: P = 0 and C_tilde = C.  Return without touching the draft.
        if( options.AUDIT.RANK_SYMMETRY )
          _disp( 3 ) << "  [trace-proj] k=0 (rank(B)=" << rank << " of " << m
                    << "): projector is the zero matrix, C_tilde == C, draft unmodified\n";
        return true;
      }

      // ---------------------------------------------------------------------------
      // k > 0: COMPUTE AND VERIFY the projector.  Still NOT APPLIED -- see 3b.
      // ---------------------------------------------------------------------------
      // Splitting compute-and-verify from apply means that if the applied version breaks
      // something, the projector arithmetic has already been checked in isolation.  This
      // revision must therefore still be bit-identical on EVERY model, k=0 or not.
      // null(B) = null(B^T B), and B^T B is m x m symmetric PSD -- so the null space can be
      // had from a SYMMETRIC eigendecomposition of an m x m matrix, never touching the
      // receiver-row dimension at all.  rev19 used the full SVD and MBC5 (B is 4568x2116)
      // had to allocate a 4568x4568 U it never used; rev20's svd_econ then failed to return
      // the full V.  rev37 removes the dependence on nr entirely by decomposing instead.
      //
      // Squaring costs conditioning, but only for DECIDING which directions are null, and
      // the gap there is enormous: across the corpus the smallest retained singular value is
      // 1.6e-5*sigma_max and the largest null one ~7e-17, so squared that is 1e-10 against
      // 1e-33.  Nowhere near the cut.
      arma::uword const nc = (arma::uword)m;
      // The guard was set at 3000 when the null space came from a FULL SVD whose discarded U
      // dominated the cost.  With svd_econ(B) the cost is O(nr*nc^2) and the guard is now
      // conservative: the corpus sweep left PDE5 and PDE9 unprojected at m=3272, k=36, purely
      // because they sat 9% over the limit.  INTERFACE.TRACE_PROJ_MAX makes it a decision rather than a
      // literal; 0 disables the guard entirely.
      if( options.INTERFACE.TRACE_PROJ_MAX && nc > options.INTERFACE.TRACE_PROJ_MAX ){
        if( options.AUDIT.RANK_SYMMETRY )
          _disp( 3 ) << "  [trace-proj] k=" << _traceProjK << " but m=" << nc
                    << " exceeds INTERFACE.TRACE_PROJ_MAX=" << options.INTERFACE.TRACE_PROJ_MAX
                    << ": projector not computed\n";
        return true;
      }

      // Z comes from an SVD of B ITSELF, not of the Gram matrix.
      //
      // rev23 took Z from eig_sym(B^T B) and P densified: PDE2 nnz 64 -> 5840 (max_row 2 ->
      // 142), MBC5 nnz -> 39687 (max_row 694), while MMPDE9 was unchanged at 204.  P is
      // basis-independent only if Z spans the same subspace, so it did not.
      //
      // Eigenvector accuracy for a near-null cluster is eps*lmax/gap, NOT eps -- and squaring
      // halves the digits in the vectors as well as squaring the gap.  Measured against the
      // eigenvalue lists rev23 prints:
      //     MBC5    gap 6.21e-13  ->  predicted 3.5e-4,  observed |B Z| 8.5e-4
      //     PDE2    gap 1.64e-10  ->  predicted 1.3e-6,  observed        1.0e-8
      //     MMPDE9  gap 3.63e-06  ->  predicted 6e-11,   observed        3.1e-8  (harmless)
      // rev22 defended the ORDERING of the eigenvalues, which does survive; the ACCURACY of
      // the eigenvectors is what builds Z, and it does not.
      // rev37: COMPONENT-WISE.  null(B) = (+) null(B_comp) over the connected components of
      // the claim <-> receiver-row bipartite graph, so the null space is a direct sum of
      // independent small problems and Z is block-diagonal by construction.
      //
      // rev36 densified B and ran one svd_econ on it.  For PDE5/PDE9 that is nr x 3272 --
      // O(nr*nc^2) ~ 6e10 flops and ~10^8 bytes -- and it showed: PDE5 went 58 s -> 75 s wall
      // and 7m15 -> 11m11 user when INTERFACE.TRACE_PROJ_MAX was raised to 6000 and the projection ran
      // on that model for the first time.  The answer it computes is 72 disjoint 2x2 blocks
      // (nnz(P)=144, max_row_nnz=2, density 1.35e-05).
      //
      // Measured block structure, four models: PDE2 nnz 64 / max_row 2, MBC5 72 / 2,
      // PDE9 144 / 2, MMPDE9 204 / 10.  The connectivity that matters is local.
      //
      // CAVEAT: P's block structure is the NULL-SPACE structure, not the connectivity
      // structure -- a large component may have a trivial null space.  So the guard is now
      // per COMPONENT: an oversized component is skipped with a report, and the rest of the
      // model is still projected, rather than the whole model being skipped for one block.
      arma::uword const kk = (arma::uword)_traceProjK;
      double smax = 0.;
      arma::mat Z( nc, kk, arma::fill::zeros );

      {
        // Union-find over tau columns: two claims are connected if they share a receiver row.
        std::vector<size_t> par( m );
        for( size_t i = 0; i < m; ++i ) par[i] = i;
        std::function<size_t(size_t)> findp = [&]( size_t a ){
          while( par[a] != a ){ par[a] = par[par[a]]; a = par[a]; } return a; };
        std::map<size_t,size_t> firstcol;                 // receiver row -> a tau touching it
        for( size_t j = 0; j < m; ++j )
          for( auto const& e : blk.bcol[j] ){
            auto it = firstcol.find( e.first );
            if( it == firstcol.end() ) firstcol.emplace( e.first, j );
            else { size_t const a = findp( it->second ), b = findp( j );
                   if( a != b ) par[a] = b; }
          }
        std::map<size_t,std::vector<size_t>> comp;
        for( size_t j = 0; j < m; ++j ) comp[ findp( j ) ].push_back( j );

        size_t zc = 0, ncomp_done = 0, ncomp_skip = 0, biggest = 0;
        for( auto const& kv : comp ){
          auto const& cols = kv.second;
          biggest = std::max( biggest, cols.size() );
          std::map<size_t,arma::uword> rmap;             // receiver rows of this component
          for( size_t j : cols )
            for( auto const& e : blk.bcol[j] )
              if( !rmap.count( e.first ) ){
                arma::uword const nn = (arma::uword)rmap.size(); rmap[ e.first ] = nn; }
          arma::uword const cr = (arma::uword)rmap.size(), cc = (arma::uword)cols.size();
          if( !cr || !cc ) continue;
          if( options.INTERFACE.TRACE_PROJ_MAX && cc > options.INTERFACE.TRACE_PROJ_MAX ){
            ++ncomp_skip;
            _disp( 2 ) << "OCFESLV::setup ** trace projection: component of " << cc
                      << " claims exceeds INTERFACE.TRACE_PROJ_MAX=" << options.INTERFACE.TRACE_PROJ_MAX
                      << "; its null directions are NOT projected" << std::endl;
            continue;
          }

          arma::mat Bc( cr, cc, arma::fill::zeros );
          for( arma::uword jj = 0; jj < cc; ++jj )
            for( auto const& e : blk.bcol[ cols[jj] ] ) Bc( rmap[ e.first ], jj ) += e.second;

          arma::mat Uc, Vc; arma::vec sc;
          if( !arma::svd_econ( Uc, sc, Vc, Bc ) || Vc.n_cols < cc ){
            std::cerr << "OCFESLV::setup ** trace projection: svd_econ failed on a component ("
                      << cr << "x" << cc << "); projector not computed, plan unaffected"
                      << std::endl;
            return true;
          }
          smax = std::max( smax, sc.n_elem ? sc(0) : 0. );
          // Rank of THIS component at the global relative cut, so k sums to the MGS value.
          double const cmax = sc.n_elem ? sc(0) : 0.;
          arma::uword rc_ = 0;
          for( arma::uword i = 0; i < sc.n_elem; ++i ) if( sc(i) > 1.0e-9 * cmax ) ++rc_;
          for( arma::uword t2 = rc_; t2 < cc; ++t2 ){
            if( zc >= kk ) break;
            for( arma::uword jj = 0; jj < cc; ++jj )
              Z( (arma::uword)cols[jj], zc ) = Vc( jj, t2 );
            ++zc;
          }
          ++ncomp_done;
        }

        if( options.AUDIT.RANK_SYMMETRY )
          _disp( 3 ) << "  [trace-proj] components=" << comp.size() << " (largest " << biggest
                    << " claims), decomposed=" << ncomp_done << " skipped=" << ncomp_skip
                    << ", null vectors collected=" << zc << " of k=" << kk << "\n";

        if( zc != kk ){
          _disp( 2 ) << "OCFESLV::setup ** trace projection: component decomposition produced "
                    << zc << " null direction(s), expected " << kk
                    << "; projector not computed, plan unaffected" << std::endl;
          return true;
        }
      }

      // Bd is still needed for the |B Z| verification below -- and only for that.
      arma::uword const nr = (arma::uword)blk.brow_ndx.size();
      arma::mat Bd( nr, nc, arma::fill::zeros );
      for( size_t j = 0; j < m; ++j )
        for( auto const& e : blk.bcol[j] ) Bd( (arma::uword)e.first, (arma::uword)j ) += e.second;

      arma::mat const P = Z * Z.t();                        // BASIS-INDEPENDENT

      // ---- verification -------------------------------------------------------------
      // Z really spans null(B); P really is an orthogonal projector; trace(P) == k.
      //
      // RELATIVE residual, against B itself.  rev23 used tolv = 1e-8*max(1,smax), which inflates with the
      // matrix norm: MBC's SAT coefficients reach 5.9e+05, so |B Z| = 8.5e-4 PASSED a check
      // that existed to catch exactly that.  A verification that accepts bad input is worse
      // than no verification.
      double const nBZ  = ( smax > 0. ? arma::abs( Bd * Z ).max() / smax
                                      : arma::abs( Bd * Z ).max() );

      // rev36: check ORTHONORMALITY OF Z, not idempotency of P.
      //
      // Z comes from svd_econ(B), so its columns are orthonormal by construction; given
      // Z^T Z = I, P = Z Z^T satisfies P^2 = P and P = P^T IDENTICALLY and trace(P) = k.
      // So one m x k check subsumes three m x m ones:
      //     |Z^T Z - I|   O(m k^2) = 2116 * 18^2  ~ 7e5 flops
      //     |P^2 - P|     O(m^3)   = 2116^3       ~ 9.5e9 flops   <- was run on EVERY setup
      // The full-corpus sweep spent 13589 s with MBC6fig_ebal at 4276 s and MBC5fig_ebal at
      // 1416 s, WITHOUT CRONOS_AUDIT -- so that was the unconditional path, and this was the
      // dominant term in it.  |P - P^T| measured exactly 0.000e+00 on every model: zero by
      // construction, carrying no information.
      //
      // trace(P) = k is retained: it is O(m), it is an independent arithmetic path, and it
      // caught nothing only because nothing was wrong.
      arma::mat const ZtZ = Z.t() * Z;
      double nZO = 0.;
      for( arma::uword i = 0; i < ZtZ.n_rows; ++i )
        for( arma::uword j = 0; j < ZtZ.n_cols; ++j )
          nZO = std::max( nZO, std::fabs( ZtZ(i,j) - ( i == j ? 1.0 : 0.0 ) ) );
      double const trP  = arma::trace( P );
      bool const ok_BZ   = ( nBZ <= 1.0e-10 );
      bool const ok_idem = ( nZO <= 1.0e-10 );   // Z^T Z = I  =>  P^2 = P and P = P^T
      bool const ok_sym  = true;                 // identical by construction, see above
      bool const ok_tr   = ( std::fabs( trP - double(_traceProjK) ) <= 1.0e-8 );
      double const nPP = nZO, nPS = 0.;          // reported below under their old names

      // ---- sparsity of P and the predicted C_tilde fill-in --------------------------
      double pmax = arma::abs( P ).max();
      double const ptol = 1.0e-10 * ( pmax > 0. ? pmax : 1.0 );
      std::vector<std::vector<size_t>> Prow( m );
      size_t nnzP = 0, maxrowP = 0, npart = 0;
      for( arma::uword i = 0; i < nc; ++i ){
        for( arma::uword j = 0; j < nc; ++j )
          if( std::fabs( P(i,j) ) > ptol ) Prow[i].push_back( (size_t)j );
        nnzP += Prow[i].size();
        maxrowP = std::max( maxrowP, Prow[i].size() );
        if( P(i,i) > 1.0e-10 ) ++npart;
      }

      if( options.AUDIT.RANK_SYMMETRY ){
        _disp( 3 ) << "  [trace-proj] k=" << _traceProjK << " (rank(B)=" << rank << " of "
                  << m << ")  participating claims=" << npart << "\n";

      // ---------------------------------------------------------------------------
      // rev120: NAME the participating claims, do not merely count them.
      // ---------------------------------------------------------------------------
      // supp(P) is the whole deficiency, and until now nothing said WHICH claims are in
      // it -- only how many.  That gap made two wrong answers possible in one session:
      //
      //   (a) a structural rule "a claim dir/Dx_S is implied by dir/S when x != dir"
      //       predicted 286 implied claims on PDE3 against a measured k=16 -- 18x over.
      //       Acting on it would have removed 270 genuine continuity conditions.
      //   (b) the refuting number, rank(C)=952 = the full claim count, was already in
      //       PLAN METRIC.  C has FULL ROW RANK: no claim is implied by another claim.
      //       The redundancy the audit reports is against the PHYSICAL rows.
      //
      // On PDE3 the participating set is 32 of 952 claims -- 3.4%.  Whether those 32 fall
      // into a nameable pattern decides whether a STRUCTURAL rule is possible at all, or
      // whether the deficiency is numerical and only a numerical test will find it.  That
      // question cannot be answered from a count.
      if( options.AUDIT.RANK_SYMMETRY && npart ){
        _disp( 3 ) << "  [zsupport] the " << npart << " claim(s) in supp(P), i.e. the ones"
                     " whose symbol couplings are collinear:\n";
        size_t shown = 0;
        for( arma::uword i = 0; i < nc; ++i ){
          if( !( P(i,i) > 1.0e-10 ) ) continue;
          if( shown++ >= 64 ){
            _disp( 3 ) << "      ... " << ( npart - 64 ) << " more (capped at 64)\n";
            break;
          }
          _disp( 3 ) << "      tau[" << std::setw(5) << i << "]  "
                    << ( (size_t)i < blk.tau_key.size()
                           ? [&]() -> std::string {
                               auto const& kk = blk.tau_key[(size_t)i];
                               size_t const sid = std::get<1>( kk ), did = std::get<3>( kk );
                               std::string dn = "dom?";
                               for( auto const& kv : _mDom )
                                 if( static_cast<size_t>( kv.first.id().second ) == did ){
                                   dn = kv.first.name(); break; }
                               std::ostringstream o;
                               o << "state=" << state_name_of_id( sid ) << " dir=" << dn
                                 << " iel_lo=" << std::get<4>( kk )
                                 << " face="   << std::get<5>( kk )
                                 << " kind="   << std::get<6>( kk );
                               return o.str();
                             }()
                           : std::string( "<index beyond tau_key>" ) )
                    << "   P(i,i)=" << std::scientific << std::setprecision(3) << P(i,i)
                    << "  partners:";
          for( size_t j : Prow[(size_t)i] )
            if( j != (size_t)i ) _disp( 3 ) << " tau[" << j << "]";
          _disp( 3 ) << "\n";
        }
      }
        _disp( 3 ) << "  [trace-proj] verify:  |B Z|max/smax=" << std::scientific
                  << std::setprecision(3) << nBZ << ( ok_BZ ? " OK" : " ** FAIL **" )
                  << "   |Z^TZ-I|max=" << nPP << ( ok_idem ? " OK" : " ** FAIL **" )
                  << "   |P-P^T|=0 by construction"
                  << "   trace(P)=" << std::fixed << std::setprecision(6) << trP
                  << ( ok_tr ? " OK" : " ** FAIL **" ) << std::defaultfloat << "\n";
        _disp( 3 ) << "  [trace-proj] P: nnz=" << nnzP << "  max_row_nnz=" << maxrowP
                  << "  density=" << std::scientific << std::setprecision(2)
                  << ( double(nnzP) / ( double(m)*double(m) ) ) << std::defaultfloat << "\n";
        if( options.INTERFACE.IMPOSITION != Options::IC_TRACE )
          _disp( 3 ) << "  [trace-proj] projector COMPUTED AND VERIFIED; application is"
                       " IC_TRACE-only at this stage, so the plan is unmodified here\n";
      }

      // A correct P is block-diagonal over the participating claims, so nnz(P) cannot exceed
      // npart^2.  Measured with a good basis: PDE2 64 against npart=32 (1024), MMPDE9 204
      // against npart=22 (484), MBC5 12 at probe against npart=6 (36).  rev23's contaminated
      // basis gave PDE2 5840 against 1024 and MBC5 39687 against 9409 -- so this is the
      // cheapest available tell that Z does not span null(B), and it fires before anything
      // downstream can consume the projector.
      bool const ok_nnz = ( nnzP <= npart * npart );

      if( !ok_nnz )
        _disp( 2 ) << "OCFESLV::setup ** trace projection: nnz(P)=" << nnzP
                  << " exceeds npart^2=" << ( npart * npart )
                  << " -- P is not block-diagonal over the participating claims, so Z does"
                     " not span null(B)" << std::endl;

      if( !ok_BZ || !ok_idem || !ok_sym || !ok_tr || !ok_nnz ){
        std::cerr << "OCFESLV::setup ** trace projection: verification FAILED"
                     " (|B Z|max/smax=" << nBZ << " |Z^TZ-I|max=" << nPP
                  << " |P-P^T|max=" << nPS << " trace(P)=" << trP
                  << " k=" << _traceProjK << ").  Projector not applied, plan unaffected."
                  << std::endl;
        return true;
      }

      // ---------------------------------------------------------------------------
      // rev60 PROMOTION-SET SCAN -- MEASUREMENT ONLY, the plan is not touched here.
      // ---------------------------------------------------------------------------
      // Handoff rev7 §2 names this as the measurement that must precede any promotion work:
      // `k > 0` under IC_STRONG occurs on 28 corpus drivers, not one, and whether supp(P)
      // spills onto SCHUR-ELIMINATED taus on those is UNMEASURED.  That spill is the entire
      // blast radius -- a plan whose participating claims are all already keep-explicit would
      // be untouched by promotion, whereas PDE2 (16 keep-explicit claims, kernel harmless as
      // gauge, max|r| = 3e-10) is a working plan that promotion could disturb.
      //
      // The partition is exactly: participating claim (P(i,i) > tol) that IS in
      // draft.explicit_tau_claims -> already explicit, promotion is a no-op for it;
      // participating claim that is NOT -> Schur-eliminated with Z-support, i.e. precisely the
      // tau the rule says cannot be eliminated.
      //
      // KEY FORM.  draft.tau_claims (whose keys index blk.tau_key) and draft.explicit_tau_claims
      // are both written with _freeze_interface_claim_key(cand.slot) at the split, so they are
      // directly comparable.  An earlier promotion attempt compared a _trace_slot_key against a
      // frozen key and silently matched nothing, so the agreement is CHECKED rather than
      // assumed: n_explicit_seen must equal draft.explicit_tau_claims.size().
      {
        size_t n_part_explicit = 0, n_part_elim = 0, n_explicit_seen = 0;
        for( auto const& ck : draft.explicit_tau_claims )
          if( blk.tau_ndx.count( ck ) ) ++n_explicit_seen;
        for( arma::uword i = 0; i < nc; ++i ){
          if( !( P(i,i) > 1.0e-10 ) ) continue;   // same predicate as npart above
          if( draft.explicit_tau_claims.count( blk.tau_key[ (size_t)i ] ) ) ++n_part_explicit;
          else ++n_part_elim;
        }
        bool const key_ok   = ( n_explicit_seen == draft.explicit_tau_claims.size() );
        bool const count_ok = ( n_part_explicit + n_part_elim == npart );
        // rev287 (R2c) -- EXACT_NATURAL: a participating claim with SILENCED rescued terms is restored, not sprung.
        // The redundancy pass may be vetoed for this block (rectangular symbol, PDE38); the projector's support is
        // the dependency measurement that always runs.  Same rebuild as rev279/rev286.
        if( _knob_EXACT_NATURAL() && !_wdecideRetry ){
          std::set<t_FrozenInterfaceClaimKey> silenced;
          for( auto const& term : draft.weak_sat_terms )
            if( term.multiplier_silenced ) silenced.insert( _trace_slot_key( term.claim_key ) );
          size_t n_new = 0;
          for( arma::uword i = 0; i < nc; ++i ){
            if( !( P(i,i) > 1.0e-10 ) ) continue;
            auto const& kk = blk.tau_key[ (size_t)i ];
            if( silenced.count( kk ) && !_strong_dropped_claims.count( kk ) && _exact_natural_restored.insert( kk ).second ) ++n_new;
          }
          if( n_new ){
            _strong_restore_pending = true;      // a RESTORE (projector support), not a drop: its own retry
            _disp( 2 ) << "OCFESLV::setup ** trace projection: EXACT_NATURAL -- projector support (k=" << _traceProjK
                      << ") reaches " << n_new << " claim(s) whose rescued receivers were silenced; RESTORING them and rebuilding"
                      << std::endl;
            return false;
          }
        }
        // rev281 -- STRONG_PROJECT: the participating eliminated taus become explicit on a promotion rebuild.
        // rev282: the keep-explicit fault injector (OCFE_hybrid_adequacy) exists to make the elimination stall;
        // a promotion would rescue exactly the claim it hides.  The injector wins, and says so.
        bool const fault_injected = Options::_env_flag( "CRONOS_FORCE_BAD_KEEP_EXPLICIT", false );
        if( options.INTERFACE.IMPOSITION == Options::IC_STRONG && _knob_STRONG_PROJECT() && fault_injected && n_part_elim > 0 )
          _disp( 2 ) << "OCFESLV::setup ** trace projection: IC_STRONG -- promotion SUPPRESSED (CRONOS_FORCE_BAD_KEEP_EXPLICIT is set)" << std::endl;
        if( options.INTERFACE.IMPOSITION == Options::IC_STRONG && _knob_STRONG_PROJECT() && !fault_injected
         && key_ok && count_ok && n_part_elim > 0 ){
          size_t n_new = 0;
          for( arma::uword i = 0; i < nc; ++i ){
            if( !( P(i,i) > 1.0e-10 ) ) continue;
            if( draft.explicit_tau_claims.count( blk.tau_key[ (size_t)i ] ) ) continue;
            if( _strong_promoted_claims.insert( blk.tau_key[ (size_t)i ] ).second ) ++n_new;
          }
          if( n_new ){
            _strong_promotion_pending = true;
            _disp( 2 ) << "OCFESLV::setup ** trace projection: IC_STRONG -- projector support reaches " << n_new
                      << " Schur-eliminated tau(s) (k=" << _traceProjK << "); PROMOTING them to explicit and rebuilding"
                      << std::endl;
            return false;
          }
        }
        // rev61 SCOPE FIX.  rev60 printed this partition in BOTH exact modes, and under
        // IC_TRACE it can only ever say one thing: the two keep-explicit insert sites are
        // guarded on IC_STRONG, so `explicit_tau_claims` is empty under IC_TRACE, hence
        // part_explicit = 0 and part_eliminated = participating BY CONSTRUCTION -- and rev60
        // then tagged every TRACE row "SPILL".  There are no Schur-eliminated taus under
        // IC_TRACE for the scan to be measuring in the first place.  Same defect class as the
        // row-count identity that derived both sides from one source: a partition with only
        // one reachable answer.  The corpus grep must therefore count STRONG rows only, and
        // the 29 TRACE rows from the rev60 sweep carry no information.
        bool const strong = ( options.INTERFACE.IMPOSITION == Options::IC_STRONG );
        _disp( 3 ) << "  [promote-scan] mode=" << ( strong ? "STRONG" : "TRACE" )
                  << " m=" << m << " k=" << _traceProjK
                  << " n_explicit=" << draft.explicit_tau_claims.size()
                  << " participating=" << npart
                  << " max_row_nnz=" << maxrowP;
        if( !strong )
          _disp( 3 ) << "  part_explicit=n/a part_eliminated=n/a"
                       "  [no Schur elimination in this mode -- promotion does not apply]";
        else
          _disp( 3 ) << "  part_explicit=" << n_part_explicit
                    << " part_eliminated=" << n_part_elim
                    << ( n_part_elim ? "  [SPILL: promotion would change this plan]"
                                     : "  [no spill: promotion inert here]" );
        _disp( 3 ) << std::endl;
        if( !key_ok )
          _disp( 2 ) << "OCFESLV::setup ** [promote-scan] KEY-FORM MISMATCH: only "
                    << n_explicit_seen << " of " << draft.explicit_tau_claims.size()
                    << " keep-explicit claim(s) resolve in the dense tau index -- the scan's"
                       " partition is not trustworthy on this model." << std::endl;
        if( !count_ok )
          _disp( 2 ) << "OCFESLV::setup ** [promote-scan] PARTITION MISMATCH: "
                    << n_part_explicit << "+" << n_part_elim << " != participating="
                    << npart << std::endl;
      }

      // ---------------------------------------------------------------------------
      // rev61 MEASUREMENT 1 -- ANATOMY OF THE NULL DIRECTIONS.  Read-only.
      // ---------------------------------------------------------------------------
      // THE QUESTION.  Is the claim set genuinely over-determined, or is the corner-claim
      // GENERATION wrong for a moving mesh?  Those need different fixes and the compatibility
      // test cannot tell them apart: `rel = 1.0000` says the posed conditions conflict, not
      // why they were posed.
      //
      // THE DISCRIMINATOR.  The kernel count law k = (nel_xi-1)*(2*NEL_T-1) is a CORNER count,
      // and a corner is where seams in two different directions meet.  So: does each null
      // direction couple claims living on DIFFERENT domains (a xi-face claim tied to a t-face
      // claim), or does it stay inside one domain?
      //   CROSS-DOMAIN  -> the direction is a corner artefact; the claims were generated
      //                    independently per seam and are inconsistent where seams intersect.
      //                    Promotion would be treating a symptom.
      //   INTRA-DOMAIN  -> a genuinely over-determined set on a single face; dropping or
      //                    promoting is the honest fix.
      // This is a PREDICTION, not a readout to be reasoned about afterwards: the corner
      // account requires CROSS-DOMAIN to dominate.  If the directions come back intra-domain,
      // the corner-count law is a coincidence and the count-law reasoning has to be retracted.
      //
      // Aggregates are always printed; per-direction detail is capped (k reaches 105 on
      // PDE38).  Gated on AUDIT.RANK_SYMMETRY -- this is a targeted probe, not a corpus line.
      if( options.AUDIT.RANK_SYMMETRY ){
        size_t const NDUMP = 6;
        size_t n_cross = 0, n_intra = 0, n_multistate = 0;
        std::map<size_t,size_t> dom_hist, state_hist, kind_hist;
        for( arma::uword j = 0; j < kk; ++j ){
          double zmax = 0.;
          for( arma::uword i = 0; i < nc; ++i ) zmax = std::max( zmax, std::fabs( Z(i,j) ) );
          double const ztol = 1.0e-10 * ( zmax > 0. ? zmax : 1.0 );
          std::vector<size_t> supp;
          for( arma::uword i = 0; i < nc; ++i )
            if( std::fabs( Z(i,j) ) > ztol ) supp.push_back( (size_t)i );
          std::set<size_t> doms, states;
          for( size_t i : supp ){
            t_InterfaceClaimKey const ck = _unfreeze_interface_claim_key( blk.tau_key[i] );
            doms.insert( ck.dom_id );
            states.insert( ck.state_id );
            ++dom_hist[ ck.dom_id ];
            ++state_hist[ ck.state_id ];
            ++kind_hist[ (size_t)ck.kind ];
          }
          if( doms.size() > 1 ) ++n_cross; else ++n_intra;
          if( states.size() > 1 ) ++n_multistate;
          if( j < NDUMP ){
            _disp( 3 ) << "  [null-anatomy] dir " << j << "  supp=" << supp.size()
                      << ( doms.size() > 1 ? "  CROSS-DOMAIN" : "  intra-domain" ) << "\n";
            for( size_t i : supp ){
              t_InterfaceClaimKey const ck = _unfreeze_interface_claim_key( blk.tau_key[i] );
              _disp( 3 ) << "      tau " << std::setw(5) << i
                        << ( draft.explicit_tau_claims.count( blk.tau_key[i] )
                               ? "  KEEP-EXPLICIT" : "  eliminated   " )
                        << "  state=" << state_name_of_id( ck.state_id )
                        << " dom=" << ck.dom_id << " iel=" << ck.iel_lo
                        << " face=" << ck.face << " kind=" << ck.kind
                        << "  Z=" << std::scientific << std::setprecision(3)
                        << Z( (arma::uword)i, j ) << std::defaultfloat << "\n";
            }
          }
        }
        _disp( 3 ) << "  [null-anatomy] k=" << _traceProjK
                  << "  cross-domain=" << n_cross << "  intra-domain=" << n_intra
                  << "  multi-state=" << n_multistate
                  << ( n_cross > 0
                         ? "  [CROSS-DOMAIN FAMILY PRESENT -- candidate for promotion]"
                         : "  [intra-only: lifting-alias family, no promotion candidate]" )
                  << "\n";
        // rev64 WORDING FIX.  rev61-63 printed "CORNER ACCOUNT NOT SUPPORTED -- retract the
        // count-law reasoning" whenever cross == 0.  That verdict was written while the corner
        // account was the live hypothesis for MMPDE13; on the 27 fixed-mesh drivers with
        // cross == 0 no count law was ever claimed, so it instructed a retraction of reasoning
        // nobody had advanced -- 92 of 97 corpus lines carried it.  The question is settled
        // (the count law decomposes as intra = nel_xi-1, cross = 2(nel_xi-1)(NEL_T-1),
        // confirmed by forward prediction at two t-meshes), so the line now reports what the
        // reader actually needs: whether a promotion candidate exists on this model.
        _disp( 3 ) << "  [null-anatomy] claims by domain:";
        for( auto const& kv : dom_hist ) _disp( 3 ) << "  dom" << kv.first << "x" << kv.second;
        _disp( 3 ) << "\n  [null-anatomy] claims by state:";
        for( auto const& kv : state_hist )
          _disp( 3 ) << "  " << state_name_of_id( kv.first ) << "x" << kv.second;
        _disp( 3 ) << "\n  [null-anatomy] claims by kind:";
        for( auto const& kv : kind_hist ) _disp( 3 ) << "  kind" << kv.first << "x" << kv.second;
        _disp( 3 ) << std::endl;
      }


      // ---------------------------------------------------------------------------
      // APPLY (Stage 1, IC_TRACE only).  C_tilde = (I-P)C, and the P block.
      // ---------------------------------------------------------------------------
      // The whole change is a BUILD-TIME rewrite of trace_constraints[].coeff.  _eval_eqn
      // already computes  val = sum_c coeff_c * var[c]  over arbitrary columns and validates
      // c < _nCollVar, and under IC_TRACE the tau columns ARE ordinary columns of var[].  So
      // the P block is expressed simply by APPENDING (tau_column, P_ij) pairs to the row's
      // coefficient list.  No solver edit, no eval edit, and the projection travels with the
      // plan automatically because it IS the plan -- no new plan fields, nothing to freeze,
      // nothing for commit_draft/copy_from/deepcopy to lose.  Plan rev5 assumed this stage
      // needed OCPlan changes and an _eval_eqn change; it needs neither.
      //
      // IC_STRONG is deliberately NOT applied here.  Its taus are not all real columns of
      // var[] -- the Schur-eliminated ones carry internal indices -- so the column mapping
      // below does not hold, and the elimination additionally assumes a zero (2,2) block.
      // That is Stage 2.
      bool strong_apply = false;
      if( options.INTERFACE.IMPOSITION == Options::IC_STRONG && _knob_STRONG_PROJECT() ){
        // rev281: every participating tau is explicit (the rebuild above), so the column mapping below holds
        // for the whole support of P and the eliminated block is untouched by the append.
        strong_apply = true;
        for( arma::uword i = 0; i < nc; ++i )
          if( P(i,i) > 1.0e-10 && !draft.explicit_tau_claims.count( blk.tau_key[ (size_t)i ] ) ) strong_apply = false;
        if( !strong_apply )
          _disp( 2 ) << "OCFESLV::setup ** trace projection: IC_STRONG -- support still on eliminated taus after "
                       "the promotion rebuild; not applied" << std::endl;
      }
      if( options.INTERFACE.IMPOSITION != Options::IC_TRACE && !strong_apply ){
        if( options.AUDIT.RANK_SYMMETRY )
          _disp( 3 ) << "  [trace-proj] IC_STRONG: projector verified but NOT applied;"
                       " the plan is unmodified\n";
        return true;
      }

      // Row and column for each tau index.  The dense tau index is claim-key order, which is
      // NOT the var[] column order -- the actual column comes from tau_claims.
      std::vector<size_t> rowof( m, std::numeric_limits<size_t>::max() );
      std::vector<size_t> colof( m, std::numeric_limits<size_t>::max() );
      for( auto const& kv : blk.ckey ){
        auto const it = blk.tau_ndx.find( kv.second );
        if( it != blk.tau_ndx.end() ) rowof[ it->second ] = kv.first;
      }
      for( size_t j = 0; j < m; ++j ){
        auto const it = draft.tau_claims.find( blk.tau_key[j] );
        if( it != draft.tau_claims.end() && !it->second.empty() ) colof[j] = it->second.front();
      }
      size_t n_unmapped = 0;
      for( size_t j = 0; j < m; ++j )
        if( rowof[j] == std::numeric_limits<size_t>::max()
         || colof[j] == std::numeric_limits<size_t>::max() ) ++n_unmapped;

      // Only TAU-BACKED rows are projected.  Rows enforced by exact replacement are not in B,
      // carry no multiplier, and must be left alone.  Measured at production mesh: MBC5 has
      // n_trace_constraint=2430 against n_trace_var=2116, i.e. 314 such rows -- two orders
      // more than the 24 seen at probe resolution, which is why this is an assertion and not
      // a comment.
      size_t n_backed = 0;
      for( auto const& kv : blk.ckey ) if( blk.tau_ndx.count( kv.second ) ) ++n_backed;
      if( draft.n_trace_constraint < n_backed
       || draft.n_trace_constraint - n_backed != draft.exact_replacement_rows.size() ){
        _disp( 2 ) << "OCFESLV::setup ** trace projection: row accounting mismatch,"
                     " n_trace_constraint=" << draft.n_trace_constraint
                  << " tau-backed=" << n_backed
                  << " exact_replacement_rows=" << draft.exact_replacement_rows.size()
                  << "; not applied" << std::endl;
        return true;
      }

      // ---------------------------------------------------------------------------
      // What is actually applied: keep ALL of C, regularise the MULTIPLIER block.
      // ---------------------------------------------------------------------------
      //     [ A     B  ]              instead of      [ A          B ]
      //     [ C   eps*P]                              [ (I-P)C     P ]
      //
      // rev25 applied the second form -- enforce r conditions exactly, drop the other k --
      // and it did everything it was designed to do: MMPDE9 tau-block deficit 5 -> 0, kernel
      // 6 -> 1, least-squares floor 1.06e-06 -> 4.75e-26, verdict INCONSISTENT -> CONSISTENT,
      // and PDE2 unharmed (dup-spread 1.388e-17 -> 2.776e-17).  But MMPDE9 then DIVERGED:
      // max|u| 1.19 -> 3.04e+03, errNode 3.04e+03.  The five conditions it stopped enforcing
      // were load-bearing, and with the multiplier pinned and the condition absent nothing
      // held those directions.  So for an INCONSISTENT model the loss is catastrophic, not
      // graceful -- which is the question plan rev5 0.2 deferred to "after Stage 1".
      //
      // The regularised form enforces EVERY condition instead.  Decomposing the constraint
      // rows:  W^T C x = -W^T r_c exactly (P W = 0, the r realisable conditions), and
      //        Z^T C x + eps*lambda_null = -Z^T r_c
      // so the k unrealisable conditions are imposed through a spring of stiffness 1/eps
      // rather than discarded.  lambda_null is now DETERMINED by the violation instead of
      // being pinned to zero, and B Z = 0 still keeps it out of the physical rows.
      //
      // eps -> 0 recovers the singular original; eps -> infinity recovers rev25's drop.
      // INTERFACE.TRACE_PROJ_EPS is the knob and its meaning is exactly that stiffness.
      //
      // Implementation is again a build-time append: C is left alone and (tau_column,
      // eps*P_ij) pairs are added to each row.  _eval_eqn already sums coeff*var[col] over
      // arbitrary columns, and under IC_TRACE the tau columns are ordinary columns of var[].
      double const eps = options.INTERFACE.TRACE_PROJ_EPS;

      // rev61 MEASUREMENT 2: retain exactly what the post-solve beta readout needs, and only
      // when asked.  Z is m x k (3272 x 36 at the largest corpus mesh, ~0.9 MB) and colof is
      // m size_t; both are dead weight in production, hence the gate.
      _betaReady = false; _betaZ.reset(); _betaColOf.clear(); _betaEps = 0.;
      _betaDirCross.clear();
      if( _knob_AUDIT_BETA() ){
        _betaZ     = Z;
        _betaColOf = colof;
        _betaEps   = eps;
        // rev63: classify each direction cross/intra here, independently of
        // AUDIT.RANK_SYMMETRY, so the beta split does not depend on a second audit flag
        // being set.  Same predicate as [null-anatomy], deliberately duplicated rather than
        // shared: if the two ever disagree the reports will say so, whereas one shared
        // helper reporting to both could be wrong in both at once.
        _betaDirCross.assign( (size_t)kk, 0 );
        for( arma::uword j = 0; j < kk; ++j ){
          double zmax = 0.;
          for( arma::uword i = 0; i < nc; ++i ) zmax = std::max( zmax, std::fabs( Z(i,j) ) );
          double const ztol = 1.0e-10 * ( zmax > 0. ? zmax : 1.0 );
          std::set<size_t> doms;
          for( arma::uword i = 0; i < nc; ++i )
            if( std::fabs( Z(i,j) ) > ztol )
              doms.insert( _unfreeze_interface_claim_key( blk.tau_key[(size_t)i] ).dom_id );
          _betaDirCross[(size_t)j] = ( doms.size() > 1 ) ? 1 : 0;
        }
        _betaReady = true;
      }

      std::map<size_t,std::vector<std::pair<size_t,double>>> addcoeff;
      size_t n_added = 0;
      for( size_t i = 0; i < m; ++i ){
        if( Prow[i].empty() ) continue;
        std::vector<std::pair<size_t,double>> out;
        out.reserve( Prow[i].size() );
        for( size_t j : Prow[i] )
          out.push_back( std::make_pair( colof[j],
                                         eps * P( (arma::uword)i, (arma::uword)j ) ) );
        n_added += out.size();
        addcoeff[ rowof[i] ] = out;
      }

      size_t n_rewritten = 0;
      for( auto& t : draft.trace_constraints ){
        auto const it = addcoeff.find( t.row );
        if( it == addcoeff.end() ) continue;      // untouched: not in the demoted subspace,
                                                  // or an exact-replacement row
        t.coeff.insert( t.coeff.end(), it->second.begin(), it->second.end() );
        ++n_rewritten;
      }
      if( n_rewritten != addcoeff.size() ){
        _disp( 2 ) << "OCFESLV::setup ** trace projection: augmented " << n_rewritten
                  << " rows, expected " << addcoeff.size()
                  << " (duplicate row ids in trace_constraints?)" << std::endl;
        return true;
      }

      _traceProjApplied = true;
      _traceProjRank    = rank;
      if( options.AUDIT.RANK_SYMMETRY )
        _disp( 3 ) << "  [trace-proj] APPLIED (" << ( strong_apply ? "IC_STRONG, promoted" : "IC_TRACE" ) << "): C kept in full; eps*P appended to "
                  << n_rewritten << " of " << m << " continuity row(s), " << n_added
                  << " coefficient(s), eps=" << std::scientific << std::setprecision(3)
                  << eps << std::defaultfloat << ".  The " << _traceProjK
                  << " unrealisable condition(s) are enforced through a spring of stiffness"
                     " 1/eps, not dropped; " << ( draft.n_trace_constraint - n_backed )
                  << " exact-replacement row(s) untouched.\n";
      return true;
    }

  //! @brief saddle-point RANK SYMMETRY audit -- rank(B) vs rank(C).
  //!
  //! The exact-imposition system is the saddle point
  //!     [ A  B ] [ dx     ]   [ -r_phys ]
  //!     [ C  0 ] [ dlambda ] = [ -r_cont ]
  //! where C is the continuity-row block over the STATE columns and B is the
  //! multiplier column block.  For a square augmented system with A nonsingular,
  //! nonsingularity requires the Schur complement C A^-1 B to be nonsingular, hence
  //!     rank(B) == rank(C) == n_trace_var.
  //!
  //! B is NOT built as C^T.  It is assembled from the IC_WEAK SAT receiver graph --
  //! principal-symbol couplings distributed over whatever receiver rows the graph
  //! offers -- while C is assembled from the claims themselves.  Nothing has ever
  //! checked that the two have equal rank.  Where a state's claims have only one host
  //! equation, or where two states share a single defining row, the receiver graph
  //! loses rank while C does not, and the augmented system silently stops being a
  //! saddle point: rank(B) < n_trace_var gives right-null vectors that are pure
  //! multiplier directions (measured: ||J v|| ~ 1e-16, 100% trace/tau energy).
  //!
  //! Whether that kernel is FATAL is model-dependent -- it makes the system singular
  //! either way, but the residual may still happen to satisfy the resulting
  //! compatibility conditions (PDE2 converges to 3e-10 carrying a 16-dimensional
  //! multiplier kernel).  So this reports by default and is fatal only under
  //! options.FATAL.INTERFACE_RANK.
  //!
  //! Cost: two sparse modified Gram-Schmidt passes over blocks whose columns/rows
  //! carry a handful of nonzeros each.  Setup-time, no solve, no dense allocation.
  //!
  //! @param os      stream to report on, or nullptr for a silent check
  //! @param p_rankB optional out: rank of the multiplier column block
  //! @param p_rankC optional out: rank of the continuity row block
  //! @return NOT "rank(B) == n_trace_var".  That was the rev-era wording and it is stale:
  //!         the code returns ONE OF TWO comparisons, chosen by whether the projector was
  //!         applied --
  //!             _traceProjApplied :  rank(B) == rank(C restricted to tau-bearing rows)
  //!             otherwise         :  rank(B) == ntau_elim  (= ntau - n_explicit_tau)
  //!         Neither is rank(B) == n_trace_var, and on OCFE_PDE3 the difference is what
  //!         the printed flag actually reports:
  //!             IC_STRONG  n_trace_var=838  rank(B)=822  n_explicit=16
  //!                        -> ntau_elim = 822, rank(B) == 822  -> symmetric=y
  //!             IC_TRACE   same ranks, projector applied
  //!                        -> compares against rank(C|tau) instead  -> symmetric=n
  //!         Same n_trace_var and same rank(B) on both lines with OPPOSITE flags, which
  //!         is impossible under the old wording and is exactly how the staleness showed.
  bool interface_rank_symmetry
    ( std::ostream* os, size_t* p_rankB = nullptr, size_t* p_rankC = nullptr )
    const
    {
      if( p_rankB ) *p_rankB = 0;
      if( p_rankC ) *p_rankC = 0;
      if( options.INTERFACE.IMPOSITION != Options::IC_TRACE
       && options.INTERFACE.IMPOSITION != Options::IC_STRONG ){
        if( os ) *os << "  [rank-sym] imposition is IC_WEAK: no multipliers, "
                        "nothing to compare\n";
        return true;
      }

      // ---- assembled blocks (shared helper; see build_interface_blocks) ---------
      InterfaceBlocks blk;
      build_interface_blocks( _plan.tau_claims_frozen(), _plan.weak_sat_terms(),
                              _plan.trace_constraints(), _plan.trace_var_offset(), blk );
      auto const& tau_ndx  = blk.tau_ndx;
      auto const& tau_key  = blk.tau_key;
      auto const& brow_ndx = blk.brow_ndx;
      auto const& bcol     = blk.bcol;
      auto const& crow     = blk.crow;
      auto const& ckey     = blk.ckey;
      size_t const ntau           = blk.ntau();
      size_t const n_tau_col_in_C = blk.n_tau_col_in_C;

      auto name_of = [&]( t_FrozenInterfaceClaimKey const& k ) -> std::string {
        size_t const sid = std::get<1>( k ), did = std::get<3>( k );
        std::string dn = "dom?";
        for( auto const& kv : _mDom )
          if( static_cast<size_t>( kv.first.id().second ) == did ){
            dn = kv.first.name(); break; }
        std::ostringstream o;
        o << "state=" << state_name_of_id( sid ) << " dir=" << dn
          << " iel_lo=" << std::get<4>( k ) << " face=" << std::get<5>( k )
          << " kind=" << std::get<6>( k );
        return o.str();
      };

      // ---- sparse modified Gram-Schmidt ---------------------------------------
      typedef std::vector<std::pair<size_t,double>> SpV;
      // Same algorithm and tolerance as the shared _sparse_column_rank(), which the
      // projection phase uses; this variant additionally collects the dependent indices for
      // reporting.  They MUST stay in step: if they disagree, the projection acts on a
      // different k from the one the audit reports.
      auto rank_of = []( std::vector<SpV> V, std::vector<size_t>* dep ) -> size_t {
        std::vector<SpV> basis;
        size_t rank = 0;
        for( size_t j = 0; j < V.size(); ++j ){
          SpV c = V[j];
          std::sort( c.begin(), c.end(),
                     []( std::pair<size_t,double> const& a,
                         std::pair<size_t,double> const& b ){ return a.first < b.first; } );
          double n0 = 0.;
          for( size_t k = 0; k < c.size(); ++k ) n0 += c[k].second * c[k].second;
          n0 = std::sqrt( n0 );
          if( n0 <= 64. * DBL_EPSILON ){ if( dep ) dep->push_back( j ); continue; }
          for( size_t bi = 0; bi < basis.size(); ++bi ){
            SpV const& q = basis[bi];
            double d = 0.;
            size_t a = 0, b = 0;
            while( a < c.size() && b < q.size() ){
              if     ( c[a].first < q[b].first ) ++a;
              else if( c[a].first > q[b].first ) ++b;
              else { d += c[a].second * q[b].second; ++a; ++b; }
            }
            if( d == 0. ) continue;
            SpV out;
            out.reserve( c.size() + q.size() );
            a = 0; b = 0;
            while( a < c.size() && b < q.size() ){
              if     ( c[a].first < q[b].first ){ out.push_back( c[a] ); ++a; }
              else if( c[a].first > q[b].first ){ out.push_back(
                         std::make_pair( q[b].first, -d * q[b].second ) ); ++b; }
              else { out.push_back( std::make_pair(
                       c[a].first, c[a].second - d * q[b].second ) ); ++a; ++b; }
            }
            for( ; a < c.size(); ++a ) out.push_back( c[a] );
            for( ; b < q.size(); ++b )
              out.push_back( std::make_pair( q[b].first, -d * q[b].second ) );
            c.swap( out );
          }
          double nr = 0.;
          for( size_t k = 0; k < c.size(); ++k ) nr += c[k].second * c[k].second;
          nr = std::sqrt( nr );
          if( nr <= 1e-9 * std::max( 1.0, n0 ) ){ if( dep ) dep->push_back( j ); continue; }
          double const inv = 1.0 / nr;
          for( size_t k = 0; k < c.size(); ++k ) c[k].second *= inv;
          basis.push_back( c );
          ++rank;
        }
        return rank;
      };

      std::vector<SpV> Bv( ntau );
      for( size_t j = 0; j < ntau; ++j )
        for( auto const& e : bcol[j] )
          Bv[j].push_back( std::make_pair( e.first, e.second ) );

      std::vector<SpV>    Cv;
      std::vector<size_t> Crow_id;
      std::vector<SpV>    Cv_tau;
      std::vector<size_t> Crow_id_tau;
      for( auto const& kv : crow ){
        SpV r;
        for( auto const& e : kv.second )
          r.push_back( std::make_pair( e.first, e.second ) );
        Cv.push_back( r );
        Crow_id.push_back( kv.first );
        auto const ik = ckey.find( kv.first );
        if( ik != ckey.end() && tau_ndx.count( ik->second ) ){
          Cv_tau.push_back( r );
          Crow_id_tau.push_back( kv.first );
        }
      }

      std::vector<size_t> depB, depC;
      size_t const rankB     = rank_of( Bv,     &depB );
      size_t const rankC     = rank_of( Cv,     &depC );
      size_t const rankC_tau = rank_of( Cv_tau, nullptr );

      if( p_rankB ) *p_rankB = rankB;
      if( p_rankC ) *p_rankC = rankC;

      // After the projection is applied the system enforces r = rank(B) conditions BY
      // DESIGN, so the right comparison is rank(B) == rank(C_tilde), not rank(B) == m.
      //
      // rev41 -- KEEP-EXPLICIT HYBRID CORRECTION.  Under IC_STRONG the plan may retain
      // n_explicit LINK-anchored taus as ordinary var[] columns; those are DELIBERATELY
      // excluded from the Schur-eliminated block, get their own rows, and drop no physical
      // row.  Only the remaining (ntau - n_explicit) taus are eliminated through B_S, so
      // the rank B must reach is that count, not ntau.
      //
      // Comparing against ntau was a defect, not a conservative choice.  The rev40 corpus
      // sweep settled it: on ALL 78 rows reported ASYMMETRIC,
      //
      //     n_explicit == k       and       rank(B) + n_explicit == n_trace_var
      //
      // held EXACTLY, across 28 drivers and five model families -- i.e. every one of them
      // was a hybrid plan whose elimination had already succeeded at exactly rank(B).  They
      // could not have been anything else: _build_frozen_strong_tau_elimination runs inside
      // _build_interface_plan_draft, and a stall there returns INTERFACE_PLAN_INVALID before
      // this report is ever reached, so a genuinely failing plan produces NO row here.  Not
      // one driver in the corpus stalls.
      //
      // Consequence to keep in view: with this corrected the corpus is symmetric in BOTH
      // modes today.  That does NOT mean the hybrid is validated -- it is a LINK-anchored
      // coordinate SELECTION where Stage 1 deliberately chose a basis-independent
      // projection, and selection is exactly what the withdrawn Stage A rule got wrong while
      // passing k / SYMMETRIC / residual / PASS.  The gate is report_duplicate_node_spreads,
      // not this verdict.
      size_t const n_explicit_tau =
        ( options.INTERFACE.IMPOSITION == Options::IC_STRONG ) ? _plan.n_trace_var_explicit()
                                                          : size_t(0);
      size_t const ntau_elim = ( ntau >= n_explicit_tau ) ? ( ntau - n_explicit_tau ) : ntau;
      bool const ok = _traceProjApplied ? ( rankB == rankC_tau ) : ( rankB == ntau_elim );


      // ---- rev10: CONDITIONING of C and B on their compressed nonzero support -------
      //
      // rank(C) = m is necessary for any exact-imposition scheme, and it is what the
      // equality-constrained least-squares formulation needs: with the KKT block
      // [A^T A, C^T; C, 0], the Schur complement C (A^T A)^-1 C^T is SPD and nonsingular
      // EXACTLY when rank(C) = m.  But "full rank" at a 1e-9 Gram-Schmidt cut says nothing
      // about how close C is to deficient, and that margin is what decides whether the
      // formulation is well posed at production meshes rather than at --small.
      //
      // Two numbers per block: the RAW spectrum, and the ROW-EQUILIBRATED one.  The raw
      // condition number of C is dominated by the spread of Lagrange trace coefficients
      // across claims and is not a meaningful invariant on its own; the equilibrated number
      // is the one to track against mesh refinement.
      //
      // Compressed support only: C's rows carry a handful of nonzeros each, so the union of
      // touched columns is far smaller than n.  Dense SVD on the compressed block, guarded
      // at 3000 in either dimension.
      auto block_spectrum = [&]( std::vector<SpV> const& V, char const* what,
                                 size_t rank_mgs )
      {
        if( !os || V.empty() ) return;
        std::map<size_t,arma::uword> cmap;
        for( auto const& r : V )
          for( auto const& e : r )
            if( !cmap.count( e.first ) ){
              arma::uword const nn = (arma::uword)cmap.size();
              cmap[ e.first ] = nn;
            }
        arma::uword const nr = (arma::uword)V.size(), nc = (arma::uword)cmap.size();
        if( !nc ) return;
        // Gram fallback: when only the SMALLER dimension fits, the singular values are the
        // square roots of the eigenvalues of the Gram matrix over that dimension.  Without
        // this, MBC at production mesh (C 2430x4580, B 2116x4568) reported nothing at all,
        // which is why its cond(C) has never been measured outside the reduced probe.
        bool const use_gram = ( nr > 3000 || nc > 3000 );
        if( use_gram && std::min( nr, nc ) > 3000 ){
          *os << "  [rank-sym] " << what << " spectrum skipped: " << nr << "x" << nc
              << " -- both dimensions exceed 3000\n";
          return;
        }
        arma::mat M( nr, nc, arma::fill::zeros );
        for( arma::uword i = 0; i < nr; ++i )
          for( auto const& e : V[i] ) M( i, cmap[ e.first ] ) += e.second;

        auto report = [&]( arma::mat const& X, char const* tagx ){
          arma::vec sv;
          bool ok = false;
          if( !use_gram ) ok = arma::svd( sv, X );
          else{
            // Gram over the smaller dimension, then sigma = sqrt(lambda), descending.
            arma::mat const Gm = ( X.n_rows <= X.n_cols ) ? arma::mat( X * X.t() )
                                                          : arma::mat( X.t() * X );
            arma::vec ge;
            ok = arma::eig_sym( ge, Gm );
            if( ok ){
              sv.set_size( ge.n_elem );
              for( arma::uword i = 0; i < ge.n_elem; ++i )
                sv(i) = std::sqrt( std::max( 0.0, ge( ge.n_elem - 1 - i ) ) );
            }
          }
          if( !ok || !sv.n_elem ){
            *os << "  [rank-sym] " << what << " " << tagx
                << ( use_gram ? " Gram eigendecomposition failed\n" : " SVD failed\n" );
            return;
          }
          double const smax = sv(0), smin = sv( sv.n_elem - 1 );
          size_t rk = 0;
          double smin_range = smax;              // smallest sigma ABOVE the null cut
          for( arma::uword i = 0; i < sv.n_elem; ++i )
            if( sv(i) > 1.0e-12 * smax ){ ++rk; smin_range = sv(i); }
          size_t const defic = (size_t)sv.n_elem - rk;
          *os << "  [rank-sym] " << what << " " << tagx << " " << nr << "x" << nc
              << "  sigma_max=" << std::scientific << std::setprecision(4) << smax
              << "  sigma_min=" << smin
              << "  cond=" << ( smin > 0. ? smax / smin : INFINITY )
              << "  rank[1e-12]=" << rk << " (MGS said " << rank_mgs << ")\n";
          // rev11: on a rank-DEFICIENT block the raw cond is just the null space and
          // carries no information (measured: cond(B) ~ 2e16 in every corpus row).  The
          // conditioning ON THE RANGE -- smallest sigma above the null cut -- is what is
          // comparable against a full-rank block such as C.
          if( defic )
            *os << "  [rank-sym] " << what << " " << tagx
                << "  RANGE (excluding " << defic << " null direction(s)):"
                << "  sigma_min=" << std::scientific << std::setprecision(4) << smin_range
                << "  cond=" << ( smin_range > 0. ? smax / smin_range : INFINITY ) << "\n";
          // Show enough values to see PAST the null block, otherwise every printed number
          // is a null direction and the first meaningful one is invisible.
          size_t const nshow = std::min<size_t>( (size_t)sv.n_elem, defic + 5 );
          *os << "      smallest: ";
          for( size_t k = 0; k < nshow; ++k )
            *os << " " << std::scientific << std::setprecision(3)
                << sv( (arma::uword)( sv.n_elem - 1 - k ) );
          *os << "\n";
        };

        report( M, "raw       " );
        arma::mat Me = M;
        for( arma::uword i = 0; i < nr; ++i ){
          double const rn2 = arma::norm( Me.row( i ) );
          if( rn2 > 0. ) Me.row( i ) /= rn2;
        }
        report( Me, "row-equil " );
      };

      block_spectrum( Cv, "C", rankC );
      block_spectrum( Bv, "B", rankB );

      // ---- rev16: does a WEAK SAT penalty survive the IC_STRONG Schur pass? -----------
      //
      // Option D demotes a component of the claim set to a weak SAT penalty weighted by
      // diag(P)_j.  Under IC_TRACE that penalty simply sits in its receiver row.  Under
      // IC_STRONG the Schur pass rewrites the receiver rows -- it DROPS pivot rows entirely
      // (they are replaced by the exact trace rows) and emits R_k - y.R_S for every non-pivot
      // row -- so a penalty's fate depends on which side of that split its row falls.
      //
      // This is the Stage 2 prerequisite.  It could NOT be answered by the withdrawn Stage A
      // experiment: Stage A filters receiver_states during claim EMISSION, upstream of the
      // receiver graph, so a suppressed claim leaves no receiver edge and no weak SAT term
      // either (measured: PDE2 receiver edges 1556 -> 1236).  There was no fallback to look
      // for.  This report measures the quantity directly instead.
      //
      // Note the penalty in a dropped pivot row is not annihilated: the projection re-injects
      // R_pivot into every non-pivot row with coefficient -y.  It is REDISTRIBUTED, not
      // preserved, and whether that retains its effect is the open question.
      if( os && options.INTERFACE.IMPOSITION == Options::IC_STRONG
       && _plan.strong_projection().ready ){
        auto const& sp = _plan.strong_projection();
        struct SatFate { size_t n=0, n_pivot=0, n_proj=0, n_plain=0;
                         double mx_pivot=0., mx_kept=0.; };
        std::map<size_t,SatFate> fate;
        size_t tot_pivot = 0, tot_proj = 0, tot_plain = 0;
        for( auto const& term : _plan.weak_sat_terms() ){
          if( std::fabs( term.prefactor ) <= 64. * DBL_EPSILON ) continue;
          FrozenPhysicalResidualRowKey rk;
          rk.row_id   = term.row_id;
          rk.block_el = term.receiver_block_el;
          rk.flat     = term.receiver_emit_flat;
          SatFate& f = fate[ term.claim_key.state_id ];
          ++f.n;
          double const a = std::fabs( term.prefactor );
          if( sp.is_pivot( rk ) ){
            ++f.n_pivot; ++tot_pivot; f.mx_pivot = std::max( f.mx_pivot, a );
          }
          else if( sp.projection.find( rk ) != sp.projection.end() ){
            ++f.n_proj; ++tot_proj;  f.mx_kept = std::max( f.mx_kept, a );
          }
          else{
            ++f.n_plain; ++tot_plain; f.mx_kept = std::max( f.mx_kept, a );
          }
        }
        *os << "  [sat-survival] IC_STRONG Schur pass: " << sp.pivot_rows.size()
            << " pivot row(s) DROPPED, " << sp.projection.size()
            << " non-pivot row(s) projected  (rank=" << sp.rank
            << " n_trace=" << sp.n_trace << ")\n";
        *os << "  [sat-survival] weak SAT edges by fate of their receiver row"
               "  (state: edges | in-DROPPED-pivot | in-projected | untouched"
               " | max|sat| dropped | max|sat| kept)\n";
        size_t sum_n = 0;
        for( auto const& kv : fate ){
          auto const& f = kv.second;
          sum_n += f.n;
          *os << "  [sat-survival]   " << state_name_of_id( kv.first )
              << " : " << f.n << " | " << f.n_pivot << " | " << f.n_proj
              << " | " << f.n_plain
              << " | " << std::scientific << std::setprecision(3) << f.mx_pivot
              << " | " << f.mx_kept << "\n";
        }

        // rev66 SUM CHECK.  The per-state table is built by a single pass over
        // _plan.weak_sat_terms(), incrementing exactly one f.n per term, so Sum(f.n) can never
        // exceed the container size.  On MMPDE13 (--small --nelt 3) the printed table summed
        // to 1756 with a 2-pass plan and 2634 with a 3-pass one, against an accounting line
        // reporting weak_sat_terms=878 in BOTH -- exactly passes x 878.  Those two facts do
        // not fit together and the mechanism is NOT yet known, so nothing is "fixed" here:
        // this line measures the discrepancy instead of guessing at it.  The two sides come
        // from genuinely different places (a running total over the report loop vs the plan's
        // own container size), so the check is not an invariant asserting itself.
        //
        // Sum > size  => the report loop runs more than once per emitted table, or the
        //                container holds more entries than the accounting line counts.
        // Sum == size => the earlier reading was a mis-parse of the log, and the plan is fine
        //                (which the a/b diff already suggested: receiver_edges,
        //                weak_sat_terms, trace_constraints and tau_claim_keys were all
        //                IDENTICAL across the extra pass).
        *os << "  [sat-survival] SUM CHECK: sum(per-state edges)=" << sum_n
            << "  weak_sat_terms().size()=" << _plan.weak_sat_terms().size()
            << ( sum_n == _plan.weak_sat_terms().size()
                   ? "  [consistent]"
                   : "  *** MISMATCH -- the per-state table does not partition the term list;"
                     " do not read absolute edge counts from it ***" ) << "\n";
        if( tot_pivot )
          *os << "      *** " << tot_pivot << " weak SAT edge(s) sit in rows the Schur pass"
                 " DROPS.  The penalty is redistributed into the non-pivot rows with"
                 " coefficient -y, not preserved in place.  A demoted Option D component"
                 " landing here would be weighted by the projection, not by diag(P). ***\n";
        else
          *os << "      *** no weak SAT edge sits in a dropped pivot row: a demoted component"
                 " keeps its penalty in place under IC_STRONG. ***\n";
      }

      // ---- rev12: FILL-IN of the sparse projection formulation --------------------
      //
      // Option D is stated as a rotation onto range(B^T), which is basis-independent but
      // whose explicit basis W (m x r, r ~ m) is DENSE: B*W and W^T*C would densify the
      // Jacobian catastrophically (PDE2: nnz(B)=838 -> ~415000).
      //
      // The sparse realisation keeps all m multipliers and uses the CANONICAL projector
      // onto null(B) instead:
      //
      //     [ A          B ] [ dx ]          P  = Z Z^T           (projector, m x m)
      //     [ C_tilde    P ] [ dl ]          C_tilde = (I - P) C
      //
      // whose two row projections are exactly  W^T C x = 0  (the r realisable conditions)
      // and  Z^T l = 0  (the k unrealisable multipliers pinned).  B is untouched.
      //
      // P and C_tilde are BASIS-INDEPENDENT -- they do not depend on which basis of null(B)
      // is computed -- so the property that justified Option D over claim-dropping survives,
      // and there is no selection, no pivot order and no non-unique demoted set.
      //
      // The design stands or falls on two numbers, measured here rather than predicted:
      //     nnz(P)                  the new (2,2) block
      //     nnz(C_tilde) - nnz(C)   the fill-in in the constraint rows
      // If P is block-diagonal with small blocks these are negligible; if null(B) is
      // delocalised over all m claims they are not, and the fallback is column-pivoted QR
      // selection with its non-uniqueness cost.
      //
      // diag(P)_j in [0,1] is claim j's participation in the lost subspace: canonical,
      // reported for EVERY participating claim, and the weight for its weak SAT penalty.
      if( os && !depB.empty() && ntau && ntau <= 3000 ){
        // Dense B over the compressed receiver-row support, then the null basis from its SVD.
        std::map<size_t,arma::uword> rmapz;
        for( auto const& col : Bv )
          for( auto const& e : col )
            if( !rmapz.count( e.first ) ){
              arma::uword const nn = (arma::uword)rmapz.size();
              rmapz[ e.first ] = nn;
            }
        arma::uword const nrz = (arma::uword)rmapz.size(), ncz = (arma::uword)ntau;
        if( !nrz ){
          *os << "  [rank-sym] projection fill-in: B has no numerically nonzero entry;"
                 " P is the identity and every claim is unrealisable\n";
        }
        else{
          arma::mat Bz( nrz, ncz, arma::fill::zeros );
          for( size_t j = 0; j < ntau; ++j )
            for( auto const& e : Bv[j] ) Bz( rmapz[ e.first ], (arma::uword)j ) += e.second;

          arma::mat Uz, Vz; arma::vec sz;
          if( !arma::svd( Uz, sz, Vz, Bz ) ){
            *os << "  [rank-sym] projection fill-in: SVD failed\n";
          }
          else{
            double const smaxz = sz.n_elem ? sz(0) : 0.;
            arma::uword rk = 0;
            for( arma::uword i = 0; i < sz.n_elem; ++i )
              if( sz(i) > 1.0e-9 * smaxz ) ++rk;
            arma::uword const kz = ncz - rk;         // = dim null(B)
            if( !kz ){
              *os << "  [rank-sym] projection fill-in: k=0, P=0, C_tilde=C (no-op)\n";
            }
            else{
              arma::mat const Z = Vz.cols( rk, ncz - 1 );      // m x k, any basis
              arma::mat const P = Z * Z.t();                   // BASIS-INDEPENDENT
              double pmax = 0.;
              for( arma::uword i = 0; i < ncz; ++i )
                for( arma::uword j = 0; j < ncz; ++j )
                  pmax = std::max( pmax, std::fabs( P(i,j) ) );
              double const ptol = 1.0e-10 * ( pmax > 0. ? pmax : 1.0 );

              size_t nnzP = 0, maxrowP = 0, npart = 0;
              std::vector<std::vector<size_t>> Prow( ncz );
              for( arma::uword i = 0; i < ncz; ++i ){
                for( arma::uword j = 0; j < ncz; ++j )
                  if( std::fabs( P(i,j) ) > ptol ) Prow[i].push_back( (size_t)j );
                nnzP += Prow[i].size();
                maxrowP = std::max( maxrowP, Prow[i].size() );
                if( P(i,i) > 1.0e-10 ) ++npart;
              }

              // Block structure of P: union-find over its nonzeros.
              std::vector<size_t> par( ncz );
              for( size_t i = 0; i < ncz; ++i ) par[i] = i;
              std::function<size_t(size_t)> find = [&]( size_t a ){
                while( par[a] != a ){ par[a] = par[par[a]]; a = par[a]; } return a; };
              for( arma::uword i = 0; i < ncz; ++i )
                for( size_t j : Prow[i] ){
                  size_t const ra = find( (size_t)i ), rb = find( j );
                  if( ra != rb ) par[ra] = rb;
                }
              std::map<size_t,size_t> blk;
              for( arma::uword i = 0; i < ncz; ++i )
                if( !Prow[i].empty() ) ++blk[ find( (size_t)i ) ];
              size_t maxblk = 0;
              for( auto const& kv : blk ) maxblk = std::max( maxblk, kv.second );

              // Fill-in of C_tilde = (I-P) C.  Row i gains the union of the supports of
              // every C row j with P(i,j) != 0.
              std::vector<SpV> Ctau( ncz );                    // C row per tau index
              for( size_t ri = 0; ri < Crow_id.size(); ++ri ){
                auto const ik = ckey.find( Crow_id[ri] );
                if( ik == ckey.end() ) continue;
                auto const it = tau_ndx.find( ik->second );
                if( it != tau_ndx.end() ) Ctau[ it->second ] = Cv[ri];
              }
              size_t nnzC = 0, nnzCt = 0;
              for( arma::uword i = 0; i < ncz; ++i ){
                nnzC += Ctau[i].size();
                if( Prow[i].empty() ){ nnzCt += Ctau[i].size(); continue; }
                std::set<size_t> u;
                for( auto const& e : Ctau[i] ) u.insert( e.first );
                for( size_t j : Prow[i] )
                  for( auto const& e : Ctau[j] ) u.insert( e.first );
                nnzCt += u.size();
              }

              *os << "  [rank-sym] projection fill-in (sparse Option D):  m=" << ncz
                  << " k=" << kz << "  participating claims=" << npart << "\n";
              *os << "      P = ZZ^T : nnz=" << nnzP
                  << "  density=" << std::scientific << std::setprecision(2)
                  << ( double(nnzP) / (double(ncz)*double(ncz)) )
                  << "  max_row_nnz=" << maxrowP
                  << "  blocks=" << blk.size() << " (largest " << maxblk << ")\n";
              *os << "      C_tilde  : nnz " << nnzC << " -> " << nnzCt
                  << "   fill-in +" << ( nnzCt >= nnzC ? nnzCt - nnzC : 0 ) << "\n";
              *os << "      participation diag(P) for the " << npart
                  << " participating claim(s):\n";
              size_t shown = 0;
              for( arma::uword i = 0; i < ncz && shown < 24; ++i )
                if( P(i,i) > 1.0e-10 ){
                  *os << "        " << std::fixed << std::setprecision(4) << P(i,i)
                      << "  " << name_of( tau_key[i] ) << "\n";
                  ++shown;
                }
              if( npart > shown ) *os << "        ... " << ( npart - shown ) << " more\n";
              *os << std::scientific;
            }
          }
        }
      }

      if( os && options.AUDIT.RANK_SYMMETRY ){
        // One line per claim, sorted by key, so two runs can be diffed directly.  Added
        // because FIX 3 was measured to remove 2 receiver edges and 0 claims (340->338,
        // 99->99) even though the leak detector showed those were the ONLY emissions at that
        // face -- so the claim is keyed differently from how the suppression is keyed, and
        // the list is the only way to see which claim survives.
        *os << "  [claim-list] " << blk.ntau() << " tau-backed claim(s):\n";
        for( size_t j = 0; j < blk.ntau(); ++j )
          *os << "      [claim] " << name_of( blk.tau_key[j] ) << "\n";
      }

      if( os ){
        *os << "  [rank-sym] imposition="
            << ( options.INTERFACE.IMPOSITION == Options::IC_TRACE ? "TRACE" : "STRONG" )
            << "  n_trace_var=" << _plan.n_trace_var()
            << " (explicit " << _plan.n_trace_var_explicit() << ")"
            << "  n_trace_constraint=" << _plan.n_trace_constraint() << "\n";
        *os << "  [rank-sym] B (multiplier columns, from the SAT receiver graph): "
            << ntau << " column(s) over " << brow_ndx.size()
            << " receiver row(s), rank = " << rankB << "\n";
        *os << "  [rank-sym] C (continuity rows over state columns): "
            << Cv.size() << " row(s), rank = " << rankC
            << "   [tau-backed subset: " << Cv_tau.size()
            << " row(s), rank = " << rankC_tau << "]\n";
        if( n_tau_col_in_C )
          *os << "  [rank-sym] note: " << n_tau_col_in_C
              << " continuity coefficient(s) referenced a multiplier column and were"
                 " excluded from C\n";

        if( !depB.empty() ){
          *os << "  [rank-sym] DEPENDENT multiplier column(s) = " << depB.size()
              << "   [NOTE: greedy Gram-Schmidt flags whichever member of a dependent GROUP"
                 " arrives last in claim-key order, so the COUNT and the SPAN are invariant"
                 " but this particular SET is not -- three passes over the same deficiency"
                 " have named three different sets.  Do not read it as 'the culprits'.]\n";
          for( size_t k = 0; k < depB.size() && k < 24; ++k )
            *os << "      B[" << depB[k] << "]  " << name_of( tau_key[ depB[k] ] )
                << "   support=" << Bv[ depB[k] ].size() << " row(s)\n";
          if( depB.size() > 24 ) *os << "      ... " << ( depB.size() - 24 ) << " more\n";
        }
        if( !depC.empty() ){
          *os << "  [rank-sym] DEPENDENT continuity row(s) = " << depC.size() << ":\n";
          for( size_t k = 0; k < depC.size() && k < 24; ++k ){
            size_t const rid = Crow_id[ depC[k] ];
            auto const ik = ckey.find( rid );
            *os << "      C row " << rid << "  "
                << ( ik == ckey.end() ? std::string( "(no claim)" ) : name_of( ik->second ) )
                << "\n";
          }
          if( depC.size() > 24 ) *os << "      ... " << ( depC.size() - 24 ) << " more\n";
        }

        if( _traceProjApplied )
          *os << "  [rank-sym] PROJECTED: rank(B)=" << rankB << " = rank(C)=" << rankC_tau
              << " of " << ntau << " multiplier(s).  The " << ( ntau - rankB )
              << " unrealisable condition(s) are enforced through the eps*P block, not"
                 " dropped; the augmented system is nonsingular by construction.\n";
        else if( ok && rankC == Cv.size() ){
          *os << "  [rank-sym] SYMMETRIC: both blocks have full rank; the augmented"
                 " system is a well-posed saddle point.\n";
          if( n_explicit_tau )
            *os << "  [rank-sym] NOTE: this plan is a keep-explicit HYBRID -- " << n_explicit_tau
                << " LINK-anchored tau(s) are retained as var[] columns and " << ntau_elim
                << " are eliminated, and rank(B)=" << rankB << " matches the ELIMINATED count."
                   "  Structural symmetry here rests on a coordinate SELECTION, not on the"
                   " basis-independent projection used under IC_TRACE, and this verdict"
                   " cannot tell a correct selection from a merely non-failing one."
                   "  Gate on report_duplicate_node_spreads().\n";
        }
        else if( !ok && rankC_tau > rankB )
          *os << "  [rank-sym] *** ASYMMETRIC: " << ( rankC_tau - rankB )
              << " claim(s) are independent as ROWS of C but dependent as COLUMNS of B:"
                 " the multiplier block cannot represent multipliers the constraint block"
                 " does emit.  rank(B) < " << ntau_elim
              << " (the eliminated tau count) makes the augmented system singular"
                 " for ANY C -- B is built from the SAT receiver graph, not from C^T."
                 "  Whether that singularity is FATAL is a separate question this report"
                 " cannot answer: if the dependent claims' rows are IMPLIED by the physical"
                 " rows the kernel is a harmless gauge freedom (PDE2), and if they are not"
                 " the system is inconsistent (MMPDE).  Read the compatibility report"
                 " (u_k^T r) in the spectrum audit at the FAILING iterate to tell them"
                 " apart. ***\n";
        else if( !ok )
          *os << "  [rank-sym] *** B is column rank-deficient by " << ( ntau_elim - rankB )
              << " and C is rank-deficient by " << ( Cv.size() - rankC )
              << ": the continuity conditions themselves are redundant, so the"
                 " deterministic fix is to drop the dependent claim(s) at EMISSION"
                 " rather than to re-route their multipliers. ***\n";
        else
          *os << "  [rank-sym] *** B has full column rank but C is rank-deficient by "
              << ( Cv.size() - rankC )
              << ": redundant continuity rows (over-determined constraint set). ***\n";
      }

      return ok;
    }

  //! @brief Effective interface dispatch for one block along one domain
  //! direction after IC_AUTO resolution (IC_VALUE / IC_FLUX / IC_UPWIND).
  //! `weak` selects the SAT (IC_WEAK) path vs the exact (IC_TRACE/IC_STRONG)
  //! path; the result is therefore imposition-mode aware.  Diagnostic use.
  Options::InterfaceType resolved_interface_type
    ( int block_id, FFVar const& dom, EqnRole role = EqnRole::INTERIOR,
      bool weak = false )
    const
    {
      // Seed from the STORED per-equation interface_type for this block, not the global
      // options.INTERFACE.TYPE -- otherwise a per-equation IC_FLUX/IC_VALUE override (mixed
      // models) is invisible to this diagnostic and it wrongly reports the global type.  Fall
      // back to the global type only if no INTERIOR equation for the block is found.
      Options::InterfaceType seed = options.INTERFACE.TYPE;
      for( auto const& eqn : _mEqn )
        if( eqn.opt->block_id == block_id && eqn.opt->role == role ){
          seed = _sopt( *eqn.opt ).interface_type; break; }
      EqnOptions opt( role, block_id, seed );
      return _resolve_interface_type( &opt, &dom, weak, block_id );
    }

protected:
  // =====================================================================
  // rev87 PROBE -- evolution-face causalisation, MEASUREMENT ONLY
  // =====================================================================
  //! @brief Collected evidence for the marching-gate probe.  NOTHING here changes behaviour.
  //!
  //! WHAT IS BEING PROBED.  The causal evolution-direction interface block (search
  //! "CAUSAL EVOLUTION-DIRECTION INTERFACE") exists to "make the monolithic solve
  //! reproduce the marching solve in the evolution direction".  It is gated on
  //! evolution_input_may_jump, a predicate about DISTRIBUTED INPUTS in t, not about the
  //! evolution direction.  MEASURED: OCFE_DAE5, OCFE_MMPDE27, OCFE_MESHONLY10 and
  //! OCFE_PDE0 all register ZERO inputs (`grep -c add_input` == 0), so on every one of
  //! them the gate is false BY CONSTRUCTION and the block is unreachable -- including on
  //! OCFE_PDE0, which the block's own comment cites as where its treatment matters
  //! ("the 26x PDE0 STRONG-mono loss").  OCFE_DAE3 has add_input=2 and IS reachable, so
  //! it is the control that must not move.
  //!
  //! THE PROPOSED PREDICATE, (ii): "marching would have been available here" --
  //! SOLVE.MARCHING on, an evolution domain known, a non-empty transfer set, a
  //! collapsible evolution grid, and no domain-less row reducing over the evolution
  //! direction.  This makes the change definitionally "build the plan marching would
  //! have built", self-limits on models where marching is refused as non-causal, and
  //! leaves SOLVE.MARCHING=false as an exact opt-out to the ungated behaviour.
  //!
  //! WHAT THIS REVISION DOES: computes the proposed predicate, counts the receiver
  //! visits it WOULD drop (split by branch: reduction-auxiliary flux vs upstream
  //! differential), and prints a per-direction claim census -- then takes the existing
  //! branch regardless.  Promoting its decision would be a separate change, gated on a corpus sweep.
  struct t_EvoCausalProbe {
    bool   ran                 = false;
    bool   gate_current        = false;  //!< evolution_input_may_jump (the marching gate)
    bool   march_enabled       = false;  //!< options.SOLVE.MARCHING
    bool   evo_dom_set         = false;  //!< _evolution_dom_set && dag()
    bool   evo_dom_user        = false;  //!< user-supplied vs auto-detected
    size_t n_dyn               = 0;      //!< differential states (d/d(evolution var))
    size_t n_transfers         = 0;      //!< ... of those, with an INITIAL closure of their own
    size_t n_lumped_evo_couple = 0;      //!< domain-less rows reducing over the evolution direction
    size_t evo_n_elem          = 0;      //!< elements in the evolution direction AT PLAN TIME
    size_t sgeqn_size          = 0;      //!< _sgEqn.size() at plan time (see the note below)
    //! UNGATED liveness tally -- every visit that reaches the causal site, keyed by
    //! face direction id.  Counters inside the gated body would themselves be the dangling code, so a
    //! zero there would prove nothing.  A zero HERE proves the site is not reached.
    std::map<size_t,size_t> visits_by_dir;
    bool   computed      = false; //!< the fields below were MEASURED, not skipped
    bool   grid_ok       = false; //!< evolution grid well-formed
    // (2026-10-06, WORKPLAN 1.4: causal_level -- CRONOS_CAUSAL_EVOFACE as parsed -- and the n_acted* counters of the
    // claims-based causal block are gone with that block, retired at rev349; the switch had controlled nothing since)
  };
  mutable t_EvoCausalProbe _evoProbe;

  //! @brief Build or load the setup-frozen IC_TRACE/IC_STRONG claim basis.
  //!
  //! This method derives retained trace claims and receiver edges from the
  //! IC_WEAK SAT receiver graph.  _eval_eqn() only consumes executable frozen
  //! vectors built from this setup-owned plan.
  void _prepare_interface_plan_tables
    ( std::map<t_InterfaceClaimKey,size_t>& claims,
      InterfacePlanDraft& draft,
      bool store_plan = true ) const;

  //! @brief Build the complete interface plan into a transactional draft.
  bool _build_interface_plan_draft
    ( InterfacePlanDraft& draft );

  //! @brief Build the frozen receiver graph used by exact-row and tau consumers.
  bool _build_interface_receiver_graph
    ( InterfacePlanDraft& draft ) const;

  //! @brief Count and freeze IC_TRACE exact constraints and tau columns.
  bool _build_frozen_trace_constraints_and_tau_terms
    ( InterfacePlanDraft& draft ) const;

public:
  //! @brief [taucol/plan]: per-column receiver census off the COMMITTED plan,
  //! plus the per-seam multiplier-capacity census.
  //!
  //! ALSO builds _tauHostRows, the lookup table the assembly half needs.  Building it
  //! here rather than at eval time is deliberate: the assembly half must not re-derive
  //! which rows host a tau, only observe what happens to the rows this half named.  An
  //! invariant that derives both sides from one source cannot detect a discrepancy
  //! between them -- but a LOOKUP KEY is not a derivation, it is the subject of the
  //! observation, and it is exactly what makes the two halves comparable.
  void _report_tau_column_plan( std::ostream& os ) const
    {
      _tauColPlanReady = false;   // cleared FIRST: a rebuild must not leave the
                                  // assembly half attached to the previous table.
      _tauHostRows.clear();
      _tauInjectCount.clear();
      _tauFateEmit.clear();  _tauFatePivot.clear();  _tauFateExact.clear();
      _tauPivotOfCol.clear();  _tauCorrRefs.clear();
      _tauPivotCoef.clear();   _tauNetMax.clear();
      _tauColAssemblyReported = false;

      bool const strong = ( options.INTERFACE.IMPOSITION == Options::IC_STRONG );
      bool const trace  = ( options.INTERFACE.IMPOSITION == Options::IC_TRACE  );
      if( !strong && !trace ){
        os << "  [taucol/plan] IC_WEAK: no multiplier columns; nothing to trace.\n";
        return;   // _tauColPlanReady stays false: nothing to attach to.
      }
      if( !_plan.ready() || !_plan.weak_sat_terms_ready() ){
        os << "  [taucol/plan] plan not ready -- NOT REACHED.\n";
        return;
      }

      size_t const toff  = _plan.trace_var_offset();
      size_t const nexp  = _plan.n_trace_var_explicit();
      auto const& proj   = _plan.strong_projection();
      auto const& repl   = _plan.exact_replacement_rows();

      // The revision is part of the marker ON PURPOSE.  rev107 and rev108 differ only in
      // what the seam census prints, so a build oracle keyed on the un-versioned string
      // passed happily on a rev107 binary and the extractor then reported every row as
      // "incomplete" -- a stale header showing up as missing data rather than as a stale
      // header.  Bump this string whenever the REPORT FORMAT changes, so the oracle and
      // the parser can both refuse rather than half-succeed.
      os << "\n#### [taucol/plan] committed-plan column trace"
            " (CRONOS_AUDIT_TAUCOL, " << taucol_instrument_rev()
         << ") ####\n";
      os << "  mode=" << ( strong ? "IC_STRONG" : "IC_TRACE" )
         << "  trace_var_offset=" << toff
         << "  n_trace_var=" << _plan.n_trace_var()
         << "  n_trace_var_explicit=" << nexp
         << "  n_coll_var=" << _plan.n_coll_var()
         << "  proj.ready=" << ( proj.ready ? 1 : 0 )
         << "  proj.rank=" << proj.rank
         << "  pivot_rows=" << proj.pivot_rows.size()
         << "  exact_replacement_rows=" << repl.size() << "\n";

      // ---- which var[] columns are we tracing? ---------------------------------
      // IC_STRONG: the kept-explicit taus, [toff, toff+nexp).  These are the ones
      // tier0-audit reports as ORPHAN cN.
      // IC_TRACE:  every tau is an explicit column, [toff, toff+n_trace_var).
      size_t const c_lo = toff;
      size_t const c_hi = strong ? ( toff + nexp ) : ( toff + _plan.n_trace_var() );
      if( c_hi <= c_lo ){
        os << "  [taucol/plan] no explicit multiplier columns in this plan"
              " (n_explicit=0) -- nothing to trace.  This is the EXPECTED reading for a"
              " marching solve: the collapse leaves no interior seam in the evolution"
              " direction, so no claim is minted there.\n";
        _tauColPlanReady = true;   // table is legitimately empty, not unbuilt
        os << "########################################################\n";
        return;
      }

      // ---- claim key for each traced column, from tau_claims -------------------
      std::map<size_t,t_FrozenInterfaceClaimKey> col2claim;
      for( auto const& kv : _plan.tau_claims_frozen() )
        for( size_t c : kv.second )
          if( c >= c_lo && c < c_hi ) col2claim[c] = kv.first;

      size_t n_orphan = 0, n_supported = 0, n_nopromise = 0;

      for( size_t c = c_lo; c < c_hi; ++c ){
        auto itk = col2claim.find( c );
        os << "  --- col=c" << c << "  ";
        if( itk == col2claim.end() ){
          os << "<< NO CLAIM MAPS TO THIS COLUMN >>  (tau_claims has no entry;"
                " the column was allocated but never registered)\n";
          ++n_nopromise;
          continue;
        }
        t_InterfaceClaimKey const ck = _unfreeze_interface_claim_key( itk->second );
        os << "state=" << _taucol_state_name( ck.state_id )
           << " dom=" << _taucol_dom_name( ck.dom_id )
           << " iel_lo=" << ck.iel_lo << " face=" << ck.face
           << " kind=" << ck.kind;
        auto itr = _keepExplicitReason.find( itk->second );
        if( itr != _keepExplicitReason.end() )
          os << " keep_explicit_reason=" << itr->second;
        os << "\n";

        // ---- PROMISE: every weak-SAT term naming this column ------------------
        size_t n_term = 0, n_term_zero = 0, n_term_disallowed = 0;
        size_t n_pivot = 0, n_exact = 0, n_ordinary = 0, n_ordinary_corr = 0;
        for( auto const& term : _plan.weak_sat_terms() ){
          if( term.lambda_col != c ) continue;
          ++n_term;
          FrozenPhysicalResidualRowKey rk;
          rk.row_id   = term.row_id;
          rk.block_el = term.receiver_block_el;
          rk.flat     = term.receiver_emit_flat;

          bool const is_piv = strong && proj.is_pivot( rk );
          bool const is_exc = ( repl.find( rk ) != repl.end() );
          bool const has_cor = ( proj.projection.find( rk ) != proj.projection.end() );
          char const* fate = is_piv ? "PIVOT (dropped by LOOP 2)"
                           : is_exc ? "EXACT-REPLACEMENT (dropped by LOOP 2)"
                           : has_cor ? "ordinary, Schur-corrected"
                                     : "ordinary";
          if( is_piv ){
            ++n_pivot;
            _tauPivotOfCol[c].insert( proj.pivot_index.at( rk ) );
          }
          else if( is_exc ) ++n_exact;
          else { ++n_ordinary; if( has_cor ) ++n_ordinary_corr; }

          if( !term.trace_tau_allowed ) ++n_term_disallowed;
          if( std::fabs( term.trace_prefactor ) <= 64. * DBL_EPSILON ) ++n_term_zero;

          os << "        " << _taucol_row_str( rk )
             << " recv_flat=" << term.receiver_flat
             << " emit_flat=" << term.receiver_emit_flat
             << " coupling=" << std::scientific << std::setprecision(6) << term.coupling
             << " orient=" << term.orientation
             << " trace_pref=" << term.trace_prefactor
             << " tau_allowed=" << ( term.trace_tau_allowed ? 1 : 0 )
             << "\n             -> fate: " << fate;
          if( is_piv ) os << " [pivot_index=" << proj.pivot_index.at( rk ) << "]";
          os << "\n";

          // The SAT injection indexes buf[receiver_flat]; LOOP 2 keys the row by the
          // POST-boundary-collapse flat.  They coincide unless a dimension of this
          // equation carries an LB/UB domain limit.  Say so rather than assume it.
          if( term.receiver_flat != term.receiver_emit_flat )
            os << "             ** recv_flat != emit_flat: the injection index and the"
                  " row-ownership index are in DIFFERENT tensor spaces here.\n";

          _tauHostRows[ rk ].push_back( { c, term.trace_prefactor } );
          if( is_piv )
            _tauPivotCoef[ proj.pivot_index.at( rk ) ][ c ] += term.trace_prefactor;
        }

        os << "      PROMISE: " << n_term << " receiver term(s)"
           << "  zero-coefficient=" << n_term_zero
           << "  tau_disallowed=" << n_term_disallowed
           << "  |  host-row fates: pivot=" << n_pivot
           << " exact=" << n_exact
           << " ordinary=" << n_ordinary << " (of which corrected=" << n_ordinary_corr << ")\n";

        if( n_term == 0 ){
          os << "      VERDICT(plan): NO RECEIVER TERM -- outcome (1), the claim never"
                " produced a coupling.\n";
          ++n_orphan; ++n_nopromise;
        }
        else if( n_term == n_term_zero ){
          os << "      VERDICT(plan): every term has a zero coefficient -- outcome (2).\n";
          ++n_orphan;
        }
        else if( n_ordinary == 0 ){
          os << "      VERDICT(plan): every host row is DROPPED by LOOP 2 -- outcome (3)."
                "  Whether the coupling survives depends on whether a surviving row's"
                " Schur correction references those pivots; [taucol/assembly] measures it.\n";
        }
        else{
          os << "      VERDICT(plan): " << n_ordinary
             << " host row(s) survive LOOP 2 -- the column should be non-empty.\n";
          ++n_supported;
        }
      }

      os << "  [taucol/plan] summary: traced=" << ( c_hi - c_lo )
         << "  no-promise=" << n_nopromise
         << "  plan-side-orphan=" << n_orphan
         << "  plan-side-supported=" << n_supported << "\n";

      // ---- per-seam multiplier capacity ----------------------------------------
      // Direct measurement, no model: group the tau columns by the seam they sit on
      // and take the rank of B restricted to that seam.  Reported next to the number
      // of claims minted there.  This is a CENSUS -- it decides nothing.
      _report_tau_seam_capacity( os );

      // Table complete.  Only now may the assembly half record against it.
      _tauColPlanReady = true;
      os << "  [taucol/plan] host-row table ready: " << _tauHostRows.size()
         << " row(s) named.  The assembly half attaches to the next _eval_eqn.\n";
      os << "########################################################\n";
    }

protected:
  //! @brief Item 11 Part-2: compute per-(primitive,direction) symbol interface
  //! reads for every classified block and store them in _blockInterfaceRead.
  //! Solve-free, header-only, always compiled; consumed by the symbol-AUTO
  //! resolver.  Under MC__OCFESLV_SYMBOL_BALANCE_PROBE it additionally prints the
  //! measure-before-implement PRIMARY/CORRECTOR/now table.  Order>2 is flagged.
  void _compute_interface_reads
    ( std::vector<double> const& state_ref,
      std::vector<double> const& cst_ref,
      std::vector<double> const& dom_ref,
      std::vector<double> const& input_ref );

protected:
  // -------------------------------------------------------------- MoL phase 3 (P3d): the function-side layout
  // and buffers, the counterparts of _depEqn / _nCollEqn / _sgEqn, and the two setup routines both solvers share
  //! @brief vector of output-function dependencies
  std::vector<FFDep>                 _depFct;
  //! @brief vector of output-function subgraphs
  mutable std::vector<FFSubgraph>    _sgFct;
  //! @brief Count of scalar output functions during setup
  size_t                             _nCollFct = 0;
  //! @brief Forward-AD DAG work vector used by deriv()
  mutable std::vector<OCVar<t_FAD>>  _wrkOCF;
  //! @brief Cached output evaluation domains used by eval_fct().
  std::vector< std::set<FFVar,lt_FFVar> > _fctEvalDom;
  //! @brief Build dependency update plans for residual and output evaluation.
  void _setup_eval_plans
    ();
  //! @brief Test whether a variable is part of the vector of constants
  bool _constant
    ( FFVar const& var )
    const;

  //! @brief The collocated layout: states, inputs, constants and domains counted, blocked and given their
  //! dependencies; the domains' nodes set; the inputs' own discretisations validated against the domains; the
  //! node cache built.  Everything FFModel::setup needs from _on_model_discretise, shared by the solvers.  A
  //! marching solver validates a distributed input on the evolution direction against ITS grid through the two
  //! virtuals below.
  bool _discretise_colloc_model();
  //! @brief Marching hook: validate input @p v's discretisation @p d on a direction the solver has collapsed.
  //! Returns true when it handled (validated) the input; the default never does.
  virtual bool _validate_input_grid_marching( FFVar const& v, FFDom const& d ) const { (void)v; (void)d; return false; }
  //! @brief Marching hook: (re)build input @p dvar's own node cache @p idom against the solver's grid.
  //! Returns 0 when not handled, 1 when handled, -1 when handled and node setup failed (status set).
  virtual int _init_input_nodes_marching( FFVar const& dvar, FFDom& idom ) { (void)dvar; (void)idom; return 0; }
  //! @brief The plan over the collocated layout: consistency count, dependency maps, eval plans, SAT coupling, the
  //! interface plan (draft, validate, commit), the collocation buffers, the tau-column report and the FFDep audit
  //! of the strong interface rows.  The generic head of _on_model_ready, shared by the solvers.
  bool _build_colloc_plan();
  //! @brief The collocation base's state reset -- domains, classification, SAT coupling, auxiliaries, eval plans,
  //! buffers, counts, the plan.  Each solver's _reset calls it before its own state.
  void _reset_colloc_state();

  // -------------------------------------------------------------- MoL phase 4 (P4b): the sparse Jacobian
  //! @brief Cached row-wise sparsity pattern of the equation Jacobian. Columns are stacked as [collocated states, collocated inputs]. Thus columns 0.._nCollVar-1 correspond to state coefficients and columns _nCollVar.._nCollVar+_nCollInp-1 correspond to input coefficients.
  std::vector< std::vector<size_t> > _jacPatternEqn;

  //! @brief Cached row-wise sparsity pattern of the output-function Jacobian. Columns are stacked as [collocated states, collocated inputs]. Thus columns 0.._nCollVar-1 correspond to state coefficients and columns _nCollVar.._nCollVar+_nCollInp-1 correspond to input coefficients.
  std::vector< std::vector<size_t> > _jacPatternFct;

  //! @brief Per-row linearity flags derived from FFDep.
  std::vector<bool>                  _jacLinearEqn;
  std::vector<bool>                  _jacLinearFct;
  bool                               _jacAllEqnLinear = false;
  bool                               _jacAllFctLinear = false;




  //! @brief Cached color assigned to each stacked derivative column
  std::vector<unsigned>              _jacColor;

  //! @brief Number of colors in the cached compressed-forward seed matrix
  unsigned                           _jacNColor = 0;

  //! @brief True when _jacPattern/_jacColor match the current setup/options
  bool                               _jacReady = false;

  //! @brief Option/classification signature used to build the derivative cache
  Options::InterfaceType             _jacInterfaceType  = Options::IC_AUTO;
  Options::ImpositionType            _jacImpositionType = Options::IC_STRONG;
  bool                               _jacClassified     = false;
  //! @brief Greedy coloring of the column-intersection graph induced by pattern
  std::vector<unsigned> _color_jacobian_columns
    ( std::vector< std::vector<size_t> > const& pattern, size_t const ncol, unsigned& ncolor )
    const;
  //! @brief Helper for copying derivative sparsity pattern
  static void _deriv_copy_pattern
    ( std::vector<std::vector<size_t>> const& pattern, size_t* nnz, size_t** colnz );
  //! @brief Helper for extracting results from AD arrays
  template <typename U>
  static void _deriv_extract_rows
    ( FADType<U> const* val, std::vector< std::vector<size_t> > const& pattern,
      std::vector<unsigned> const& color, size_t const ncol, U* res, U** grad );
  //! @brief Templated sparse derivative evaluation of collocated equation residuals in compressed FADBAD forward mode
  template <typename U>
  void _deriv_eqn
    ( U* eqn, U** gradeqn, size_t* nnzeqn, size_t** colnzeqn,
      OCVar<FADType<U>>& OCFeqn, std::vector<OCVar<FADType<U>>>& OCFvar,
      std::vector<OCVar<FADType<U>>>& OCFinp, std::vector<OCVar<FADType<U>>>& OCFdom,
      std::vector<OCVar<FADType<U>>>& OCFcst, std::vector<OCVar<FADType<U>>>& OCFwrk,
      FADType<U>* Feqn, FADType<U>* Fvar, FADType<U>* Finp )
    const;
  //! @brief The Jacobian's sparsity pattern per row and its column colouring, from one dependence pass through the
  //! OCVar environment (FFDep), with each row's linearity.  OCFESLV's _setup_deriv_cache wraps it; OCFEMOL's
  //! deriv_residual reads it.  The output-function pass goes through _eval_fct_dep, which the solver overrides.
  bool _build_jacobian_pattern
    ();
  //! @brief The output functions' dependence pass for _build_jacobian_pattern: _eval_fct in the dependence environment
  //! (OCFESLV's, moved here with _eval_fct -- OCFEMOL's outputs need it as much).
  virtual void _eval_fct_dep
    ( std::vector<FFDep>& fct_dep, std::vector<FFDep> const& var_dep, std::vector<FFDep> const& inp_dep )
    {
      _eval_fct( fct_dep.data(), var_dep.empty()? nullptr: var_dep.data(),
                 inp_dep.empty()? nullptr: inp_dep.data(),
                 _OCDepEqn, _OCDepVar, _OCDepInp, _OCDepDom, _OCDepCst, _wrkOCDep );
    }

};

inline OCFEBASE::Options::InterfaceType
OCFEBASE::_resolve_interface_type
( FFModel::EqnOptions const* eqnopt, FFVar const* pdom, bool weak, int block_id )
const
{
  // Ordinary IC/BC/interface trace residuals are coupled only by trace-value
  // continuity along their free dimensions.  IC_TRACE is weak-like here:
  // it distributes claim-centred multipliers through the same receiver graph
  // as IC_WEAK, so IC_AUTO must not select characteristic/upwind coupling
  // for physical trace rows unless the user requests it explicitly.
  bool const weak_like = weak || options.INTERFACE.IMPOSITION == Options::IC_TRACE;
  if( weak_like && eqnopt && _ordinary_trace_role( eqnopt->role )
   && _sopt( *eqnopt ).interface_type == Options::IC_AUTO )
    return Options::IC_VALUE;

  auto itcls = _blockClassification.find( block_id );
  t_Classify const* cls = nullptr;
  if( itcls != _blockClassification.end() ) cls = &itcls->second;
  else if( _classified ) cls = &_classification;

  auto itsym = _blockSymbol.find( block_id );
  t_Symbol const* sym = nullptr;
  if( itsym != _blockSymbol.end() ) sym = &itsym->second;
  else if( _classified ) sym = &_symbol;

  auto is_evolution_face = [&]() -> bool
  {
    if( !pdom ) return false;
    if( cls && sym && cls->evolution_dom_idx < sym->vDom.size() )
      return sym->vDom[ cls->evolution_dom_idx ].id() == pdom->id();
    return _evolution_dom_set && pdom->id() == _evolution_dom_var.id();
  };

  // Return the first DISTRIBUTED input that may jump across element interfaces
  // in direction dom -- i.e. one that is collocated in dom with an interior
  // interface (n_elem >= 2) and is NOT declared C0-continuous there (the
  // DISCONTINUOUS default, or an explicit level below CONTINUOUS_C0).  A jumping
  // input forces the matched state's first derivative to jump, so an imposed
  // C1/flux match would contradict the true (kinked) solution.  nullptr if none.
  //
  // Domain-aware by construction: iterating _mInp (distributed inputs only) and
  // requiring dom in the input's own domain set means constants (_vCst), inputs
  // with no domain, and inputs distributed only in OTHER directions are ignored
  // -- exactly as they must be, since none of them can jump across an interface
  // in dom.  Reuses the same _input_dom + n_elem<2 test _check_input_continuity
  // uses, so the notion of "distributed across an interface here" is identical.
  auto discontinuous_input_in = [&]( FFVar const& dom ) -> FFVar const*
  {
    for( auto const& [inp, idoms] : _mInp ){
      bool here = false;
      for( auto const& dv : idoms ) if( dv.id() == dom.id() ){ here = true; break; }
      if( !here ) continue;                              // not distributed in dom
      if( _input_dom( inp, dom ).n_elem < 2 ) continue;  // no interior interface
      int level = static_cast<int>( DISCONTINUOUS );     // default when undeclared
      auto itc = _mInpCont.find( inp );
      if( itc != _mInpCont.end() ){
        auto itd = itc->second.find( dom );
        if( itd != itc->second.end() ) level = itd->second;
      }
      if( level < static_cast<int>( CONTINUOUS_C0 ) ) return &inp;  // may jump
    }
    return nullptr;
  };

  Options::InterfaceType type = options.INTERFACE.TYPE;
  if( eqnopt && _sopt( *eqnopt ).interface_type != Options::IC_AUTO )
    type = _sopt( *eqnopt ).interface_type;

  if( type != Options::IC_AUTO ){
    // Evolution-slab interfaces are causal trace-continuity constraints.
    // Do not impose legacy C1/flux or upwind jumps in the evolution direction.
    if( ( type == Options::IC_FLUX || type == Options::IC_UPWIND )
     && is_evolution_face() )
      return Options::IC_VALUE;

    // In a reduced-order direction the auxiliary state already represents the
    // physical flux-like quantity.  Legacy IC_FLUX would add a separate C1
    // derivative jump on top of C0 auxiliary continuity and over-count rows.
    if( type == Options::IC_FLUX && pdom
     && _block_has_aux_link_in_direction( block_id, *pdom ) )
      return Options::IC_VALUE;

    // Legacy IC_FLUX matches the state's first derivative (C1) across the
    // interface.  That is valid only where the derivative is genuinely
    // continuous; a distributed input that may jump in this direction breaks it
    // (the physical flux stays continuous, so the state derivative kinks).  If
    // such an input exists, downgrade to the always-safe C0 value match rather
    // than impose a derivative match the true solution violates, and tell the
    // user how to recover C1.
    //
    // Conservative on purpose: this fires on any distributed, not-C0-declared
    // input in the direction, without proving the input sits specifically in
    // the flux coefficient (vs a source term) or is incident to this block.
    // IC_FLUX is a deprecated opt-in path; a safe C0 fallback plus an actionable
    // warning is preferable to a blind flux-term dependency walk here.  A future
    // tightening could restrict this to flux-coefficient / block-incident inputs
    // if the conservative form is found to over-fire.
    if( type == Options::IC_FLUX && pdom ){
      if( FFVar const* culprit = discontinuous_input_in( *pdom ) ){
        std::cerr << "OCFESLV: WARNING: IC_FLUX (C1/flux) requested in direction "
                  << pdom->name() << ", but distributed input '" << culprit->name()
                  << "' is not declared continuous there (DISCONTINUOUS default), "
                  << "so its jump would break first-derivative continuity. "
                  << "Downgrading IC_FLUX -> IC_VALUE (C0) for this direction. "
                  << "Declare the input CONTINUOUS_C0/SMOOTH_C1 in this direction "
                  << "(add_input overload) if it is actually continuous to keep C1."
                  << std::endl;
        return Options::IC_VALUE;
      }
    }

    return type;
  }

  // Weak AUTO remains deliberately conservative.  Upwind SAT requires a
  // separate stable SAT form and should be requested explicitly.
  bool const surface = eqnopt && eqnopt->role == EqnRole::SURFACE;
  if( weak_like && !surface ) return Options::IC_VALUE;

  if( is_evolution_face() ) return Options::IC_VALUE;
  if( !cls ) return Options::IC_VALUE;

  // In a reduced-order direction the auxiliary state already carries the
  // flux/C1 quantity through its own C0 continuity, so value continuity is all
  // the interface needs (legacy IC_FLUX would over-count rows).  This mirrors
  // the explicit-type path above and makes reduced order-2 directions resolve
  // to VALUE regardless of the (sometimes coarse) block label.
  if( pdom && _block_has_aux_link_in_direction( block_id, *pdom ) )
    return Options::IC_VALUE;

  // Block-level legacy verdict, used as the fallback whenever the per-direction
  // symbol read is unavailable or cannot govern: upwind only genuinely
  // evolution-hyperbolic transverse faces; everything else value-continuity.
  auto block_level = [&]() -> Options::InterfaceType {
    if( cls->evolution_hyperbolic && !cls->weak_hyperbolic )
      return Options::IC_UPWIND;
    return Options::IC_VALUE;
  };

  // Item 11 Part-2: per-direction principal-symbol read.  Symbol-AUTO presupposes
  // order reduction (the read is blind to non-reduced high-order operators), so
  // under RED_NONE defer to the block-level verdict and let an explicit
  // interface_type govern the direction instead.
  if( options.REDUCE.ORDER == Options::RED_NONE || !pdom )
    return block_level();

  auto itrd = _blockInterfaceRead.find( block_id );
  if( itrd == _blockInterfaceRead.end() || !sym )
    return block_level();                       // no read cached for this block

  int dir = -1;
  for( size_t i = 0; i < sym->vDom.size(); ++i )
    if( sym->vDom[i].id() == pdom->id() ){ dir = (int)i; break; }
  if( dir < 0 )
    return block_level();

  bool any = false, all_order1 = true, any_order_gt2 = false, any_adv_dom = false;
  for( auto const& rd : itrd->second ){
    if( rd.dir_idx != dir ) continue;
    any = true;
    if( rd.order_gt2 )           any_order_gt2 = true;
    if( rd.order >= 2 )          all_order1    = false;
    if( rd.advection_dominated ) any_adv_dom   = true;
  }
  // No primitive couples through this direction (algebraic), or an unhandled
  // order>2 chain is present: defer to the block-level verdict.
  if( !any || any_order_gt2 )
    return block_level();

  // Order-1 (characteristic) directions upwind only in a genuinely
  // evolution-hyperbolic block; order-2 (parabolic) directions take value
  // continuity.  A singularly advection-dominated order-2 direction behaves
  // order-1 at this mesh and is folded into the order-1 arm (a no-op under
  // reduction, where the aux-link rule above already returned VALUE).
  bool const hyp = cls->evolution_hyperbolic && !cls->weak_hyperbolic;
  bool const effectively_order1 = all_order1 || any_adv_dom;
  if( hyp && effectively_order1 )
    return Options::IC_UPWIND;
  return Options::IC_VALUE;
}

inline void
OCFEBASE::_init_reference
( double* var, double* inp, double const* cst )
const
{
  if( _nCollVar && !var ) throw Exceptions( Exceptions::INTERNAL );
  if( _nCollInp && !inp ) throw Exceptions( Exceptions::INTERNAL );

  auto lookup_scalar = []( FFVar const& var,
                           std::map<FFVar,double,lt_FFVar> const& ref,
                           double& val ) -> bool
  {
    for( auto const& [v,x] : ref )
      if( v.id() == var.id() ){ val = x; return true; }
    return false;
  };

  auto lookup_fun = []( FFVar const& var,
                        std::map<FFVar,t_Fun,lt_FFVar> const& ref,
                        t_Fun const*& fun ) -> bool
  {
    for( auto const& [v,f] : ref )
      if( v.id() == var.id() ){ fun = &f; return true; }
    fun = nullptr;
    return false;
  };

  t_Coord default_coord;
  for( auto const& [dvar,dom] : _mDom ){
    (void)dom;
    default_coord[dvar] = _default_dom_ref( dvar );
  }

  auto fill_block_reference = [&]( FFVar const& v,
                                   std::set<FFVar,lt_FFVar> const& dom,
                                   std::tuple<size_t,size_t,std::vector<size_t>> const& blk,
                                   std::map<FFVar,double,lt_FFVar> const& scalar_ref,
                                   std::map<FFVar,t_Fun,lt_FFVar> const& fun_ref,
                                   double* data,
                                   bool is_state )
  {
    double scalar = 0.;
    bool const have_scalar = lookup_scalar( v, scalar_ref, scalar );
    t_Fun const* fun = nullptr;
    bool const have_fun = lookup_fun( v, fun_ref, fun );

    auto const& [blk_pos, blk_size, blk_dim] = blk;
    (void)blk_pos; (void)blk_dim;
    if( dom.empty() ){
      size_t const inc = is_state ? pos_state( v, std::map<FFVar,size_t,lt_FFVar>() )
                                  : pos_input( v, std::map<FFVar,size_t,lt_FFVar>() );
      data[inc] = have_fun ? (*fun)( default_coord )
                           : ( have_scalar ? scalar : 0. );
      return;
    }

    std::vector<std::vector<double>> anodes( dom.size() );
    std::map<FFVar,size_t,lt_FFVar> ndx_el;
    for( ; _next_block( ndx_el, dom ) && !ndx_el.empty(); ){
      size_t nnodes = 1, idep = 0;
      for( auto const& d : dom ){
        auto const& cdom = is_state ? _mDom.at( d ) : _input_dom( v, d );
        auto lgn = _node_block( cdom, ndx_el.at( d ) );
        nnodes *= lgn.size();
        anodes[idep++] = std::move( lgn );
      }
      if( nnodes != blk_size ) throw Exceptions( Exceptions::INTERNAL );

      size_t const inc = is_state ? pos_state( v, ndx_el )
                                  : pos_input( v, ndx_el );
      // Per-element direction/element/count vectors (constant over this element's nodes), used by
      // the interior-side re-sampling of a jumping distributed input below.
      std::vector<FFVar> ifdirs; std::vector<size_t> ifelem, ifcnt;
      if( have_fun && !is_state ){
        ifdirs.reserve( dom.size() ); ifelem.reserve( dom.size() ); ifcnt.reserve( dom.size() );
        size_t jj = 0; for( auto const& d : dom ){ ifdirs.push_back( d );
          ifelem.push_back( ndx_el.count(d)? ndx_el.at(d) : 0 ); ifcnt.push_back( anodes[jj].size() ); ++jj; }
      }
      std::vector<size_t> idx( dom.size(), 0 );
      for( size_t inode = 0; inode < nnodes; ++inode ){
        t_Coord coord = default_coord;
        size_t idom = 0;
        for( auto const& d : dom )
          coord[d] = anodes[idom][idx[idom]], ++idom;

        double node_val = have_fun ? (*fun)( coord )
                                   : ( have_scalar ? scalar : 0. );
        // A distributed INPUT may be DISCONTINUOUS across element interfaces (the default).  The
        // two abutting boundary nodes at an interior interface are DISTINCT DOFs, but a scalar ref
        // function returns one value for the shared face coordinate -- so a jumping input would give
        // the downstream element its upstream neighbour's value.  Re-sample the ref function on the
        // element's INTERIOR side at interior-interface boundary nodes, but only if it actually jumps
        // there (else keep the exact-face value, so continuous inputs are byte-for-byte unchanged).
        if( have_fun && !is_state )
          node_val = _input_ref_interior_side( *fun, coord, node_val, v, ifdirs, ifelem, idx, ifcnt );
        data[inc+inode] = node_val;

        for( size_t i=0; i<dom.size(); ++i ){
          if( ++idx[i] < anodes[i].size() ) break;
          idx[i] = 0;
        }
      }
    }
  };

  if( var && _nCollVar ) std::fill( var, var + _nCollVar, 0. );
  if( inp && _nCollInp ) std::fill( inp, inp + _nCollInp, 0. );

  // Primitive states first.  Auxiliary states are filled afterwards from their
  // tracked defining expressions so they remain discretely consistent with the
  // primitive reference values.
  for( auto const& [v,dom] : _mVar ){
    if( _auxVarID.find( v.id() ) != _auxVarID.end() ) continue;
    auto itblk = _blkVar.find( v );
    if( itblk == _blkVar.end() ) throw Exceptions( Exceptions::INDEX );
    fill_block_reference( v, dom, itblk->second, _classVarRef, _classVarRefFun, var, true );
  }

  for( auto const& [v,dom] : _mInp ){
    auto itblk = _blkInp.find( v );
    if( itblk == _blkInp.end() ) throw Exceptions( Exceptions::INDEX );
    fill_block_reference( v, dom, itblk->second, _classInpRef, _classInpRefFun, inp, false );
  }

  for( auto const& aux : _auxDef )
    _init_auxiliary_reference( aux, var, inp, cst );

  // A distributed input discretised on the evolution direction is filled by fill_block_reference only
  // on the RESIDENT element: its block iteration follows _mDom, which the marching setup collapses to a
  // single-element template, so el>0 come out zero.  For a REGISTERED CONTROL that under-seed matters --
  // encode_controls would read a truncated profile and the base marching solve would depend on the
  // per-window re-sample to populate el>0.  Re-fill each registered control across its OWN full element
  // grid via sample_input (the exact fill the solve uses), so init seeds every window's control DOFs and
  // decode/encode round-trip the complete profile.  No-op for mono (fill_block_reference already spans
  // all elements) and for non-control inputs (left to the per-window feed/BC re-sample as before).
  for( auto const& kv : _classInpRefFun )
    if( kv.second && _is_registered_control( kv.first ) )
      sample_input( kv.first, kv.second, inp );

  // 2026-10-01: the same under-seed for an input that is NOT a registered control, or whose reference is a
  // CONSTANT: under marching it stayed 0 after the first window (q = 1 0 0 0), and the march then used those
  // zeros -- nothing re-samples it per window.  Re-fill every input distributed along the evolution direction
  // across its full element grid (monolithic: already filled, left untouched).
  if( inp && _evolution_window().stepwise && _evolution_dom_set && _evolution_dom_var.dag() ){   // a horizon solved element by element
    for( auto const& [v, doms] : _mInp ){
      bool on_evol = false;
      for( auto const& d : doms ) if( d.id() == _evolution_dom_var.id() ){ on_evol = true; break; }
      if( !on_evol || _usr._mInpFixUsr.count( v ) ) continue;
      auto const itf = _classInpRefFun.find( v );
      if( itf != _classInpRefFun.end() && itf->second ){
        if( !_is_registered_control( v ) ) sample_input( v, itf->second, inp );   // controls: done above
        continue;
      }
      auto const itc = _classInpRef.find( v );
      if( itc != _classInpRef.end() ){
        double const c = itc->second;
        sample_input( v, [c]( t_Coord const& ){ return c; }, inp );
      }
    }
  }
}

inline void
OCFEBASE::_init_auxiliary_reference
( t_AuxDef const& aux, double* var, double const* inp, double const* cst )
const
{
  // A TRANSITION auxiliary (w = x(tau+), lifted by FFModel::_lower_transitions / _derive_transition_lifts) is defined by
  // its transition row, not by the parent it is registered with -- whose domain also spans the evolution direction, so
  // evaluating it on w's spatial elements has no evolution index (v1.5 profiles threw map::at here).  Leave it at its
  // default: the solve determines it.
  for( auto const& L : _trnLift ) if( L.aux.id() == aux.aux.id() ) return;
  if( !var ) throw Exceptions( Exceptions::INTERNAL );

  auto itblk = _blkVar.find( aux.aux );
  if( itblk == _blkVar.end() ) throw Exceptions( Exceptions::INDEX );
  size_t const blk_size = std::get<1>( itblk->second );

  std::vector<OCVar<double>> OCvar, OCinp, OCdom, OCcst, OCwrk;
  _reset_colloc( OCvar, OCinp, OCdom, OCcst );

  if( !_vCst.empty() ){
    if( cst ){
      for( size_t i=0; i<_vCst.size(); ++i ) OCcst[i].set( cst[i] );
    }
    else if( !_vCstVal.empty() ){
      if( _vCstVal.size() != _vCst.size() ) throw Exceptions( Exceptions::CSTVAL );
      for( size_t i=0; i<_vCst.size(); ++i ) OCcst[i].set( _vCstVal[i] );
    }
    else{
      for( size_t i=0; i<_vCst.size(); ++i ) OCcst[i].set( 0. );
    }
  }

  FFSubgraph sg = _dag->subgraph( 1, &aux.expr );

  // A REDUCTION auxiliary (its expr contains an OpI/OpEval that consumes a
  // direction NOT present in aux.dom) cannot be initialised by this per-element
  // pass: the consumed direction is absent from ndx_el, so evaluating the expr
  // would look up a state on a direction ndx_el does not carry.  The init is only
  // a reference guess -- the solve computes the auxiliary exactly from its
  // defining LINK -- so leave it at its default and return.  (Scalar reduction
  // auxes were already safe because aux.dom is empty and the block loop below
  // never runs; this covers reduction auxes with a NON-empty residual domain,
  // e.g. w = OpEval(c,z,z0) on {y}.)
  for( auto const* sop : sg.l_op )
    if( sop && ( sop->sameid( typeid( FFIntegral ) ) || sop->sameid( typeid( FFEval ) ) ) )
      return;

  FFDep dep;
  _dag->eval( sg, 1, &aux.expr, &dep,
              _vVar.size(), _vVar.data(), _depVar.data(),
              _vInp.size(), _vInp.data(), _depInp.data(),
              _vCst.size(), _vCst.data(), _depCst.data(),
              _vDom.size(), _vDom.data(), _depDom.data() );

  OCVar<double> OCres;
  std::map<FFVar,size_t,lt_FFVar> ndx_el;
  for( ; _next_block( ndx_el, aux.dom ) && !ndx_el.empty(); ){
    for( auto const& [id,d] : dep.dep() ){
      (void)d;
      auto itdep = _mDep.find( id );
      if( itdep == _mDep.end() ) continue;
      auto const& [typ,pos] = itdep->second;
      size_t inc = 0;
      switch( typ ){
        case DOMAIN:
          OCdom[pos].update( _node_block( _vDom[pos], ndx_el ), ndx_el );
          break;
        case STATE:
          inc = pos_state( _vVar[pos], ndx_el );
          OCvar[pos].update( &var[inc], ndx_el );
          break;
        case INPUT: {
          // Marching: a distributed input matched to the evolution direction keeps its FULL
          // per-window grid, but the assembly runs on the collapsed state window (evolution
          // element 0).  Remap the evolution element to the resident window so pos_input()
          // reads THIS window's slot (control DOF k), not element 0.  Inert off marching.
          auto ndx_in = ndx_el;
          auto const ew = _evolution_window();
          if( ew.kind == EvolutionWindow::ONE && _evolution_dom_var.dag() && _mInpDisc.count( _vInp[pos] ) ){
            auto ite = ndx_in.find( _evolution_dom_var );
            if( ite != ndx_in.end() ) ite->second = ew.index;
          }
          inc = pos_input( _vInp[pos], ndx_in );
          OCinp[pos].update( &inp[inc], ndx_in );
          break;
        }
        case CONSTANT:
          break;
        default:
          throw Exceptions( Exceptions::INTERNAL );
      }
    }

    _dag->eval( sg, OCwrk,
                1, &aux.expr, &OCres,
                _vVar.size(), _vVar.data(), OCvar.data(),
                _vInp.size(), _vInp.data(), OCinp.data(),
                _vCst.size(), _vCst.data(), OCcst.data(),
                _vDom.size(), _vDom.data(), OCdom.data() );

    auto const& coeff = OCres.coef();
    size_t const inc_aux = pos_state( aux.aux, ndx_el );
    if( coeff.size() == blk_size )
      std::copy( coeff.begin(), coeff.end(), &var[inc_aux] );
    else if( !coeff.empty() ){
      // The definition carries FEWER coefficients than the auxiliary's nodes -- e.g. w = u(t) with u piecewise
      // constant (1 per element) against a 3-node state block (2026-09-28: this threw INTERNAL).  This is only the
      // reference guess -- the solve computes w from its defining LINK -- so fill the element with their mean
      // (exact for a piecewise-constant input).
      double mean = 0.;  for( double c : coeff ) mean += c;  mean /= double( coeff.size() );
      std::fill( &var[inc_aux], &var[inc_aux] + blk_size, mean );
    }
    else throw Exceptions( Exceptions::INTERNAL );
  }
}

inline double
OCFEBASE::_sample_reference
( FFVar const& var, std::set<FFVar,lt_FFVar> const& dom,
  std::tuple<size_t,size_t,std::vector<size_t>> const& blk, double const* data )
const
{
  if( !data ) return 0.;
  auto const& [blk_pos, blk_size, blk_dim] = blk;
  (void)blk_pos; (void)blk_dim;
  bool const is_state = _mVar.find( var ) != _mVar.end();
  if( dom.empty() ){
    std::map<FFVar,size_t,lt_FFVar> ndx_empty;
    size_t const inc0 = is_state ? pos_state( var, ndx_empty )
                                 : pos_input( var, ndx_empty );
    return data[inc0];
  }

  std::map<FFVar,size_t,lt_FFVar> ndx_el;
  std::map<FFVar,double,lt_FFVar> ref_el;
  for( auto const& d : dom ){
    auto const& od = is_state ? _mDom.at( d ) : _input_dom( var, d );
    auto const [iel,xi] = _point_in_domain( od, _default_dom_ref( d ) );
    ndx_el[d] = iel;
    ref_el[d] = xi;
  }

  size_t const inc = is_state ? pos_state( var, ndx_el )
                              : pos_input( var, ndx_el );
  std::vector<double> buf( data + inc, data + inc + blk_size );

  struct RDim { FFDom const* pdom; FFVar const* pvar; size_t cur_n; };
  std::vector<RDim> dims;
  dims.reserve( dom.size() );
  for( auto const& d : dom ){
    auto const& od = is_state ? _mDom.at( d ) : _input_dom( var, d );
    dims.push_back( { &od, &d, od.n_node } );
  }

  for( size_t idim = 0; idim < dims.size(); ++idim ){
    auto& rd = dims[idim];
    size_t stride = 1;
    for( size_t k=0; k<idim; ++k ) stride *= dims[k].cur_n;
    size_t outer = 1;
    for( size_t k=idim+1; k<dims.size(); ++k ) outer *= dims[k].cur_n;
    std::vector<double> newbuf;
    if( !rd.pdom->eval_lagrange( newbuf, ref_el.at( *rd.pvar ), 0, buf, stride, outer ) )
      throw Exceptions( Exceptions::INTERNAL );
    buf = std::move( newbuf );
    rd.cur_n = 1;
  }
  if( buf.size() != 1 ) throw Exceptions( Exceptions::INTERNAL );
  return buf[0];
}

inline bool
OCFEBASE::_block_has_aux_link_in_direction
( int block_id, FFVar const& dom )
const
{
  for( auto const& aux : _auxDef ){
    if( aux.block_id != block_id ) continue;
    for( auto const& d : aux.diff_dom )
      if( d.id() == dom.id() ) return true;
  }
  return false;
}

inline bool
OCFEBASE::_auxiliary_state_in_direction
( FFVar const& state, FFVar const& dom )
const
{
  for( auto const& aux : _auxDef ){
    if( aux.aux.id() != state.id() ) continue;
    for( auto const& d : aux.diff_dom )
      if( d.id() == dom.id() ) return true;
    return false;
  }
  return false;
}

inline bool
OCFEBASE::_link_equation_in_direction
( FFVar const& eqn, int block_id, FFVar const& dom )
const
{
  return _equation_has_aux_link_in_direction( eqn, block_id, dom );
}

inline bool
OCFEBASE::_equation_has_aux_link_in_direction
( FFVar const& eqn, int block_id, FFVar const& dom )
const
{
  if( !_dag || !eqn.dag() ) return false;

  auto const cache_key = std::make_tuple
    ( static_cast<size_t>( eqn.id().second ), block_id,
      static_cast<size_t>( dom.id().second ) );
  auto itcache = _eqnAuxDirCache.find( cache_key );
  if( itcache != _eqnAuxDirCache.end() ) return itcache->second;

  std::vector<t_AuxDef const*> candidates;
  candidates.reserve( _auxDef.size() );
  for( auto const& aux : _auxDef ){
    if( aux.block_id != block_id ) continue;
    bool normal = false;
    for( auto const& d : aux.diff_dom )
      if( d.id() == dom.id() ){ normal = true; break; }
    if( normal ) candidates.push_back( &aux );
  }
  if( candidates.empty() ){
    _eqnAuxDirCache[cache_key] = false;
    return false;
  }

  auto contains_candidate_aux = [&]( FFVar const& v ) -> bool
  {
    for( auto const* aux : candidates )
      if( v.id() == aux->aux.id() ) return true;
    return false;
  };

  if( contains_candidate_aux( eqn ) ){
    _eqnAuxDirCache[cache_key] = true;
    return true;
  }

  FFSubgraph sg = _dag->subgraph( 1, &eqn );
  for( auto const& op : sg.l_op ){
    for( auto const* vin : op->varin )
      if( vin && contains_candidate_aux( *vin ) ){
        _eqnAuxDirCache[cache_key] = true;
        return true;
      }
    for( auto const* vout : op->varout )
      if( vout && contains_candidate_aux( *vout ) ){
        _eqnAuxDirCache[cache_key] = true;
        return true;
      }
  }

  _eqnAuxDirCache[cache_key] = false;
  return false;
}

inline bool
OCFEBASE::_block_receive_interface_in_direction
( int block_id, FFVar const& dom )
const
{
  // EqnOptions::receive_sat controls automatic same-domain interface
  // handling in both SAT and exact-replacement paths.  This lets models with
  // explicit physical interface/boundary equations opt out cleanly.
  for( auto const& eqn : _mEqn ){
    if( eqn.opt->block_id != block_id || !_sopt( *eqn.opt ).receive_sat ) continue;
    auto it = eqn.dom.find( dom );
    if( it == eqn.dom.end() ) continue;
    if( it->second == FFDom::LB || it->second == FFDom::UB ) continue;
    return true;
  }
  return false;
}

inline bool
OCFEBASE::_weak_exact_reduced_interface_in_direction
( int block_id, FFVar const& dom )
const
{
  (void)block_id;
  (void)dom;
  // True weak formulation: reduced auxiliary interfaces are handled by SAT
  // penalties only.  No exact row replacement is performed under IC_WEAK.
  return false;
}

inline bool
OCFEBASE::_block_has_initial_face
( int block_id, FFVar const& dom, int lim )
const
{
  for( auto const& eqn : _mEqn ){
    if( eqn.opt->block_id != block_id || eqn.opt->role != EqnRole::INITIAL )
      continue;
    auto it = eqn.dom.find( dom );
    if( it != eqn.dom.end() && it->second == lim ) return true;
  }
  return false;
}

inline bool
OCFEBASE::_validate_evolution_domain
()
const
{
  // Runs at end of setup(), after _classify_pde() + _build_reduction_plan() +
  // _audit_reduced_dof().  Consolidated, human-readable checks that (1) a declared
  // or auto-detected evolution domain is a genuine one-way causal direction, and
  // (2) initial conditions can be set consistently.  Warn-by-default:
  // FATAL.EVOLUTION_DOMAIN=false downgrades every violation to a DISPLAY_LEVEL
  // diagnostic and returns true.  (Marching/MOL-only causality checks -- no
  // evolution-UB / periodic coupling -- are deferred until SOLVE_MODE exists.)
  //
  // _classify_pde() sets _evolution_dom_set / _evolution_dom_var from the
  // effective (user-set or auto-inferred) evolution domain, so _evolution_dom_set
  // is the single "an evolution direction is in play" signal; a pure spatial BVP
  // with none leaves it false and this validator no-ops.
  if( !_classified || !_evolution_dom_set || !_evolution_dom_var.dag() )
    return true;

  t_Classify const& cls = _classification;
  std::vector<std::string> issues;   // fatal-eligible (under FATAL.EVOLUTION_DOMAIN)
  std::vector<std::string> notes;    // never fatal (caveats)

  // ---- Check 1: the evolution direction is causal / suitable ----
  // Suitable (one-way causal): ODE/DAE/EVOL_HYPERBOLIC/PARABOLIC/DESCRIPTOR.
  // weak hyperbolicity (now a FLAG, rev167 D1) and UNDETERMINED are NOTES (never fatal): the former is
  // marchable-but-ill-conditioned; the latter is "non-square or inconclusive",
  // i.e. causality could not be decided (a marching mode would fall back to
  // monolithic) -- which is NOT the same as being definitely unsuitable.  Only
  // the genuinely two-way / no-evolution types (ELLIPTIC / SPATIALLY_CHARACTERISTIC
  // / COMPLEX_CHARACTERISTIC / DEGENERATE) are violations.
  // rev70: DIFFERENTIAL_RECTANGULAR is suitable ONLY on a USER-SUPPLIED evolution domain.
  //
  // The classifier could not confirm causality -- the differential core is rectangular so
  // the square-symbol test does not apply -- but the SHAPE IS KNOWN rather than unknown,
  // which is what separates it from UNDETERMINED.  Admitting it unconditionally would let
  // marching proceed on unverified causality, and a wrong causal direction under marching
  // gives a SILENTLY wrong answer -- this arc's characteristic failure.
  //
  // The measurements give the split.  On OCFE_MMPDE27 (--gdirect, nel_xi=8): the
  // user-supplied domain converges (set_evolution_domain(t) -> conv=yes) while
  // auto-detection does not (conv=no, and k=12 cross-domain kernel directions appear).
  // The caller ASSERTING the causality the classifier cannot verify is exactly the case
  // that works; auto-detection on the same model is exactly the case that does not.
  // rev163 STAGE B (knob CRONOS_RECT_PROPERTY=1, default 0 = rev162): read the SHAPE from the
  // shape field rather than from the type enum.  With STAGE A recording in force (implied by
  // this knob) the rev70 clause becomes SATISFIABLE for the first time -- measured corpus-wide
  // at rev160, 136 of 136 rectangular classify lines reported evol_dom_idx=none and ZERO
  // reported [user-supplied], so this clause has never fired and the note below, commented
  // "reached only when AUTO-detected", has always been reached.
  // The type is still overwritten (STAGE C has not run), so symbol_rectangular is true exactly
  // when the type is DIFFERENTIAL_RECTANGULAR: the ONLY behavioural difference here is that the
  // clause can now be true, on blocks whose driver asserted the evolution domain.
  // rev164: STAGE C implies STAGE B.  Once the type is no longer overwritten, the three
  // consumers MUST read symbol_rectangular or they lose these blocks entirely -- in
  // particular the At-non-singular class (PDE39, PDE39b, DAE5), whose type becomes
  // UNDETERMINED and which would then fail the suitability gate that today admits them.
  // rev168 STAGE D5: the shape lives in symbol_rectangular and nowhere else.
  bool const rect_here = cls.symbol_rectangular;
  bool const type_ok =
       cls.type == DIFFERENTIAL_ORDINARY || cls.type == DIFFERENTIAL_ALGEBRAIC
    || cls.type == DIFFERENTIAL_IMPLICIT || cls.type == EVOL_HYPERBOLIC
    || cls.type == PARABOLIC || cls.type == DESCRIPTOR
    || ( rect_here && cls.evolution_user_supplied );
  // rev167 STAGE D1: weak hyperbolicity is now a FLAG on an EVOL_HYPERBOLIC block, so the
  // block is type_ok and the warning must be emitted on the flag, outside the else-chain --
  // otherwise the note silently disappears with the enum value.  It is a NOTE, never fatal,
  // exactly as before.
  if( cls.weak_hyperbolic )
    notes.push_back( "evolution direction is weakly hyperbolic (defective/"
      "ill-conditioned eigenvectors); marching/initialization may be poorly conditioned" );

  if( type_ok ){
    /* suitable -- no diagnostic */
  }
  else if( cls.type == UNDETERMINED )
    notes.push_back( "classification is UNDETERMINED (non-square or inconclusive) "
      "-- causality could not be confirmed; a marching mode would fall back to monolithic" );
  // rev67: ALGEBRAIC_FIELD / ALGEBRAIC_LUMPED are NOTES, not violations.
  //
  // The list below treats "not one of the marchable types" as a violation, which is right
  // for ELLIPTIC / SPATIALLY_CHARACTERISTIC / COMPLEX_CHARACTERISTIC / DEGENERATE -- those are genuinely
  // two-way problems being asked to march.  An ALGEBRAIC system is different: there is no
  // evolution to be causal ABOUT.  Every node is determined pointwise, so an auto-detected
  // "evolution direction" is merely the domain the detector happened to pick, and saying
  // so is informative rather than a fault.
  //
  // ALGFIELD1 auto-detects t and then classifies ALGEBRAIC_FIELD, so the old wording
  // ("evolution domain declared/detected but the classification is ...") fired a WARNING
  // on every correct run of a model that solves to |y-y*| = 5.8e-12.
  else if( rect_here )   // reached only when AUTO-detected (rev163: true in practice at last)
    notes.push_back( "the principal symbol is RECTANGULAR (symbol_rectangular) -- the retained differential"
      " rows differentiate more states than there are rows, so the square-symbol causality"
      " test does not apply.  The evolution domain was AUTO-DETECTED and causality is"
      " unverified, so a marching mode falls back to monolithic; call"
      " set_evolution_domain() to assert it." );
  else if( cls.type == ALGEBRAIC_FIELD || cls.type == ALGEBRAIC_LUMPED )
    notes.push_back( std::string( "classification is " ) + pde_type_name( cls.type )
      + " -- there is no evolution to be causal about: every node is determined pointwise,"
        " so the detected evolution direction carries no marching semantics.  Not a fault." );
  else
    issues.push_back( std::string("evolution domain declared/detected but the "
      "classification is ") + pde_type_name( cls.type )
      + " -- no one-way causal evolution direction (marchable: "
        "DIFFERENTIAL_ORDINARY/DIFFERENTIAL_ALGEBRAIC/EVOL_HYPERBOLIC/PARABOLIC/DESCRIPTOR,"
        " plus a RECTANGULAR symbol when the evolution domain is user-supplied)" );

  // ---- Check 2: initial conditions can be set consistently ----
  // (a) reduced-system structural consistency (authoritative; reuse the audit).
  if( _dofAudit.ran && !_dofAudit.ok() )
    issues.push_back( "reduced DOF audit inconsistent (rows="
      + std::to_string(_dofAudit.rows) + ", cols=" + std::to_string(_dofAudit.cols)
      + ", deficiency=" + std::to_string(_dofAudit.deficiency)
      + ") -- no consistent IC for the reduced system" );
  // (b) high-index resolvability.
  if( !_reductionPlan.empty() && !_reductionPlan.resolved )
    issues.push_back( "high-index reduction unresolved -- consistent ICs for the "
      "hidden (differentiated) constraints cannot be constructed" );
  // (c) each evolutionary block carries an INITIAL condition at the evolution LB.
  size_t evo_blocks = 0, evo_blocks_ic = 0;
  for( auto const& [bid, bcls] : _blockClassification ){
    bool const evolutionary =
         bcls.type == DIFFERENTIAL_ORDINARY || bcls.type == DIFFERENTIAL_ALGEBRAIC
      || bcls.type == DIFFERENTIAL_IMPLICIT || bcls.type == EVOL_HYPERBOLIC
      || bcls.type == PARABOLIC || bcls.type == DESCRIPTOR
      // rev70: same user-supplied gate as the block above -- a rectangular differential
      // core is evolutionary on the caller's assertion, never on auto-detection.
      || ( bcls.symbol_rectangular && bcls.evolution_user_supplied );
    if( !evolutionary ) continue;
    ++evo_blocks;
    if( _block_has_initial_face( bid, _evolution_dom_var, FFDom::LB ) )
      ++evo_blocks_ic;
  }
  if( evo_blocks_ic < evo_blocks )
    issues.push_back( "an evolutionary block has no INITIAL condition at the "
      "evolution lower bound (" + std::to_string(evo_blocks_ic) + "/"
      + std::to_string(evo_blocks) + " blocks covered)" );

  // ---- report ----
  bool const fatal = !issues.empty() && options.FATAL.EVOLUTION_DOMAIN;
  char const* sev  = options.FATAL.EVOLUTION_DOMAIN ? "ERROR" : "WARNING";
  // rev67: the CLASSIFICATION VERDICT is provenance, not diagnostics, and is stamped
  // UNCONDITIONALLY once per block -- the same footing as the HEADER_ID stamp.
  //
  // Before this it was invisible at every runtime verbosity: the only emission site was
  // gated on the compile-time macro MC__OCFESLV_SYMBOL_BALANCE_PROBE, which no driver and
  // no build in the corpus defines.  MEASURED: OCFE_AE0 and OCFE_DAE0 printed no verdict
  // at all, while the [cov1] HEADER_ID dump ran at DEFAULT level and accounted for 42-70%
  // of their log bytes.  The signal needed a rebuild; the noise needed nothing.
  //
  // That invisibility has real cost.  Drivers advertise "gates: setup square + classify
  // OK", but under CLASS_AUTO a failed classification proceeds silently and setup()
  // still returns true, so the gate was VACUOUS and unverifiable from a log.  It also
  // hid a hard abort: a distributed purely-algebraic field crashed inside classification
  // and the only visible symptom was a SIGABRT with no verdict line to explain it.
  for( auto const& s : issues )                 // real violations: always surfaced
    _disp( 1 ) << "OCFESLV::setup [evolution] " << sev << ": " << s << "\n";
  if( _solve_display_level() >= 1 ){
    for( auto const& n : notes )
      _disp( 1 ) << "OCFESLV::setup [evolution] NOTE: " << n << "\n";
    if( issues.empty() && notes.empty() ){
      std::ostringstream ev; ev << _evolution_dom_var;
      _disp( 1 ) << "OCFESLV::setup [evolution] OK -- var=" << ev.str()
                << " type=" << pde_type_name( cls.type )
                << " sigma_min_evol=" << cls.sigma_min_evolution
                << " IC coverage=" << evo_blocks_ic << "/" << evo_blocks << "\n";
    }
  }

  return !fatal;
}

inline bool
OCFEBASE::_block_has_trace_face
( int block_id, FFVar const& dom, int lim )
const
{
  for( auto const& eqn : _mEqn ){
    if( eqn.opt->block_id != block_id ) continue;
    if( eqn.opt->role != EqnRole::INITIAL && eqn.opt->role != EqnRole::BOUNDARY
     && eqn.opt->role != EqnRole::INTERFACE )
      continue;
    auto it = eqn.dom.find( dom );
    if( it != eqn.dom.end() && it->second == lim ) return true;
  }
  return false;
}

inline bool
OCFEBASE::_state_has_value_trace_face
( int block_id, size_t state_id, FFVar const& dom, int lim )
const
{
  FFVar const* pst = nullptr;
  for( auto const& [st, dep] : _mVar ){
    (void)dep;
    if( static_cast<size_t>( st.id().second ) == state_id ){ pst = &st; break; }
  }
  if( !pst ) return false;

  // Only a direct state occurrence in the boundary/initial expression counts
  // as an exact value trace for that state.  A row such as Dz_Cl = 0 contains
  // the auxiliary state Dz_Cl directly but must not be interpreted as an
  // exact value condition for primitive Cl.
  auto expr_directly_contains_state = [&]( FFVar const& expr ) -> bool
  {
    if( expr.id() == pst->id() ) return true;
    if( !_dag || !expr.dag() ) return false;
    FFSubgraph sg = _dag->subgraph( 1, &expr );
    for( auto const& op : sg.l_op ){
      for( auto const* vin : op->varin )
        if( vin && vin->id() == pst->id() ) return true;
      for( auto const* vout : op->varout )
        if( vout && vout->id() == pst->id() ) return true;
    }
    return false;
  };

  for( auto const& eqn : _mEqn ){
    if( eqn.opt->block_id != block_id ) continue;
    if( eqn.opt->role != EqnRole::INITIAL && eqn.opt->role != EqnRole::BOUNDARY
     && eqn.opt->role != EqnRole::INTERFACE )
      continue;
    auto it = eqn.dom.find( dom );
    if( it == eqn.dom.end() || it->second != lim ) continue;
    if( expr_directly_contains_state( eqn.var ) ) return true;
  }
  return false;
}

inline bool
OCFEBASE::_keep_eqn_node
( FFVar const& eqn, FFModel::EqnOptions const& opt, FFVar const& dom, FFDom const& pdom,
  size_t iel, size_t inode, int lim, bool weak, size_t eqn_ndx )
const
{
  size_t const n_elem = pdom.n_elem;
  size_t const n_node = pdom.n_node;
  bool ok = true;

  switch( lim ){
    case FFDom::ALL:
      break;
    case FFDom::ALL-FFDom::LB:
      if( iel == 0 && inode == 0 ) ok = false;
      if( _causal_start_excludes( eqn_ndx, dom, iel, inode ) ) ok = false;   // row-level causal design: an initial node
      break;
    case FFDom::ALL-FFDom::UB:
      if( iel + 1 == n_elem && inode + 1 == n_node ) ok = false;
      break;
    case FFDom::ALL-FFDom::LB-FFDom::UB:
      if( iel == 0 && inode == 0 ) ok = false;
      if( _causal_start_excludes( eqn_ndx, dom, iel, inode ) ) ok = false;
      if( iel + 1 == n_elem && inode + 1 == n_node ) ok = false;
      break;
    case FFDom::LB:
    case FFDom::UB:
      // LB/UB are collapsed before node filtering in _eval_eqn().
      break;
    default:
      throw Exceptions( Exceptions::DOMAIN );
  }

  bool const is_link = ( opt.role == EqnRole::LINK );
  bool const auto_interface_dir = is_link
    ? _block_receive_interface_in_direction( opt.block_id, dom )
    : _sopt( opt ).receive_sat;
  bool const link_reduced_dir =
    is_link && _link_equation_in_direction( eqn, opt.block_id, dom );
  bool const eqn_reduced_dir =
    !is_link && _equation_has_aux_link_in_direction( eqn, opt.block_id, dom );
  bool const row_reduced_dir = is_link ? link_reduced_dir : eqn_reduced_dir;
  bool const exact_reduced_weak =
    weak && auto_interface_dir && row_reduced_dir;

  if( !ok ) return ok;
  // True weak formulation: retain every physical equation/LINK row that
  // survives the user domain mask; interface jumps are added by SAT below.
  // This removes the previous hybrid weak path that performed exact donor
  // replacement for reduced auxiliary directions.
  if( weak ) return ok;

  // Physical trace residuals.  Under IC_STRONG and IC_TRACE they must remain
  // available as endpoint donor candidates in the global setup-time assignment;
  // they are not locally removed by the generic exact-replacement filter merely
  // because they may receive tangential SAT under IC_WEAK.  The global donor
  // graph will use them only as admissible duplicate trace rows and will retain
  // at least one representative physical condition.
  if( opt.role == EqnRole::INTERFACE
   || opt.role == EqnRole::BOUNDARY
   || opt.role == EqnRole::INITIAL ) return ok;

  if( !is_link && !_sopt( opt ).receive_sat ) return ok;
  if( !auto_interface_dir ) return ok;

  Options::InterfaceType const dir_type = exact_reduced_weak
    ? Options::IC_VALUE
    : _resolve_interface_type( &opt, &dom, false, opt.block_id );

  bool const reduced_link_dir = exact_reduced_weak || row_reduced_dir;

  bool const is_interior_like =
       opt.role == EqnRole::INTERIOR
    || opt.role == EqnRole::SURFACE;

  if( reduced_link_dir ){
    if( is_link ){
      // Keep both one-sided traces of q - D(u) = 0 at reduced-order
      // interfaces.  This ties each element-local derivative trace to the
      // auxiliary state and prevents derivative jumps such as q_R != a*T_x,R.
      return ok;
    }
    if( opt.role == EqnRole::INITIAL ){
      // At initial corners of a reduced spatial interface, keep both
      // duplicate primitive traces.  A companion corner rule in _eval_eqn()
      // drops one reduced LINK trace, avoiding derivative compatibility
      // constraints on the elementwise polynomial initial profile while
      // preserving the square row count.
      return ok;
    }
    if( is_interior_like ){
      // The two LINK traces, together with state continuity rows, replace the
      // evolution/interior PDE traces while maintaining the same row count.
      if( dir_type == Options::IC_VALUE || dir_type == Options::IC_UPWIND ){
        if( iel > 0 && inode == 0 ) ok = false;
        if( iel + 1 < n_elem && inode + 1 == n_node ) ok = false;
      }
      else if( dir_type == Options::IC_FLUX ){
        if( iel > 0 && inode == 0 ) ok = false;
        if( iel + 1 < n_elem && inode + 1 == n_node ) ok = false;
      }
      return ok;
    }
  }

  // LINK equations that are not normal to the current interface direction
  // remain pointwise auxiliary-definition identities.  Under IC_STRONG/IC_TRACE
  // they should not be claimed by ordinary first-order replacement; under
  // IC_WEAK, SAT handling has already returned above.
  if( is_link ) return ok;

  // Legacy/default strong replacement for ordinary first-order systems.
  // For initial-value/marching closures, interface continuity replaces the
  // right element's lower trace.  For terminal closures such as
  //   g_z + f = 0 on ALL-UB,  g(UB)=gB,
  // replacing the left element's upper trace preserves the one-dimensional
  // collocation rank for LG/LGR/LGL/CGL multi-element BVPs.
  bool const has_lb_trace = _block_has_trace_face( opt.block_id, dom, FFDom::LB );
  bool const has_ub_trace = _block_has_trace_face( opt.block_id, dom, FFDom::UB );
  bool const terminal_oriented =
       ( lim == FFDom::ALL-FFDom::UB || lim == FFDom::ALL-FFDom::LB-FFDom::UB )
    && has_ub_trace && !has_lb_trace;

  if( dir_type == Options::IC_VALUE || dir_type == Options::IC_UPWIND ){
    if( terminal_oriented ){
      if( iel + 1 < n_elem && inode + 1 == n_node ) ok = false;
    }
    else{
      if( iel > 0 && inode == 0 ) ok = false;
    }
  }
  else if( dir_type == Options::IC_FLUX ){
    if( iel > 0 && inode == 0 ) ok = false;
    if( iel + 1 < n_elem && inode + 1 == n_node ) ok = false;
  }

  return ok;
}

inline bool
OCFEBASE::_check_consistency
( FFGraph* dag, std::vector<FFVar> const& vCst, std::vector<double> const* vCstVal,
  t_Dom const& mDom, t_Var const& mVar, t_Var const& mInp,
  t_Eqns const& mEqn, t_Fcts const& mFct,
  size_t* nCollEqn, std::vector<FFSubgraph>* sgEqn,
  size_t* nCollFct, std::vector<FFSubgraph>* sgFct )
const
{
  // Validation is the model's (FFModel::_check_model), which also stores the subgraphs; what follows is the
  // solver's: the collocated row counts and each output's rows.  The masks are valid by now.
  if( !_check_model( dag, vCst, vCstVal, mDom, mVar, mInp, mEqn, mFct, sgEqn, sgFct ) ) return false;

  // 2026-10-02: under the row-level causal design every element's first evolution node is an initial node: the
  // `ALL - LB` rows are not collocated there, and each state with an evolution derivative gets one exact continuity
  // row per (element start, transverse node) -- counted here, emitted at the end of _eval_eqn().
  bool const causal_rows = _causal_rows_active();
  // The set belongs to the WORKING model: filled only by the counting call (_on_model_ready), never by the user-model
  // validation (_validate_model, nCollEqn == nullptr), which FFOCFESLV embeddings also run after setup.
  if( nCollEqn && _evoCollocated ) _evoDerivRows.clear();
  // When t is not collocated (OCFEMOL) nothing below fills the set, and the clear is skipped: OCFEMOL filled it
  // BEFORE the plan is built, from peel_evolution()'s rows with a mass entry -- the same set by another route,
  // and _defining_rows reads it during the plan build.
  if( causal_rows && nCollEqn )
    for( size_t ie = 0; ie < mEqn.size(); ++ie ){
      auto const& eqn = mEqn[ie];
      if( eqn.opt && eqn.opt->role == EqnRole::BOUNDARY ){ _evoDerivRows.insert( ie ); continue; }
      FFSubgraph const sgq = dag->subgraph( 1, &eqn.var );
      for( auto const& op : sgq.l_op ){
        if( !op || !op->sameid( typeid( FFPartial ) ) ) continue;
        auto const* pop = mc::type_cast<FFPartial const>( op );
        if( !pop ) continue;
        for( auto const& [dv, o_] : pop->Indep().expr ){ (void)o_; if( dv.id() == _evolution_dom_var.id() ){ _evoDerivRows.insert( ie ); break; } }
      }
    }
  // causal_evo( dv ) is evaluated per EQUATION (see the loop): the evolution direction, and this row carries an
  // evolution derivative
  size_t cur_eqn_ndx = 0;
  auto causal_evo = [&]( FFVar const& dv ){ return causal_rows && dv.id() == _evolution_dom_var.id() && _evoDerivRows.count( cur_eqn_ndx ) > 0; };
  if( nCollEqn ){
    *nCollEqn = 0;
    if( causal_rows ){
      FFDom const& edom = mDom.at( _evolution_dom_var );
      for( auto const& st : _states_with_evolution_derivative( mEqn, mVar ) ){
        auto itv = mVar.find( st );
        if( itv == mVar.end() ) continue;
        size_t transverse = 1;
        for( auto const& dv : itv->second ){
          if( dv.id() == _evolution_dom_var.id() ) continue;
          FFDom const& d2 = mDom.at( dv );  transverse *= d2.n_elem * d2.n_node;
        }
        *nCollEqn += ( edom.n_elem - 1 ) * transverse;
      }
    }
    for( auto const& eqn : mEqn ){
      cur_eqn_ndx = static_cast<size_t>( &eqn - &mEqn[0] );
      size_t cnt = 1;
      for( auto const& [var, lim] : eqn.dom ){
        auto const& dom = mDom.at( var );
        switch( lim ){
          case FFDom::ALL:                     cnt *= dom.n_elem * dom.n_node;     break;
          case FFDom::ALL-FFDom::LB:
            cnt *= causal_evo( var )? dom.n_elem * ( dom.n_node - 1 ): dom.n_elem * dom.n_node - 1; break;
          case FFDom::ALL-FFDom::UB:           cnt *= dom.n_elem * dom.n_node - 1; break;
          case FFDom::ALL-FFDom::LB-FFDom::UB:
            cnt *= causal_evo( var )? dom.n_elem * ( dom.n_node - 1 ) - 1: dom.n_elem * dom.n_node - 2; break;
          default:                                                                 break;   // LB, UB: no change
        }
      }
      *nCollEqn += cnt;
    }
  }

#ifdef CRONOS_DEBUG_OUT
  _disp( 3 ) << "setup: output checks begin\n";
#endif
  auto keep_output_node = [&]( FFDom const& dom, int const lim,
                                size_t const iel, size_t const inode, size_t const nel ) -> bool
  {
    switch( lim ){
      case FFDom::ALL:
        break;
      case FFDom::ALL-FFDom::LB:
        if( iel == 0 && inode == 0 ) return false;
        break;
      case FFDom::ALL-FFDom::UB:
        if( iel + 1 == nel && inode + 1 == dom.n_node ) return false;
        break;
      case FFDom::ALL-FFDom::LB-FFDom::UB:
        if( iel == 0 && inode == 0 ) return false;
        if( iel + 1 == nel && inode + 1 == dom.n_node ) return false;
        break;
      case FFDom::LB:
        return iel == 0 && inode == 0;
      case FFDom::UB:
        return iel + 1 == nel && inode + 1 == dom.n_node;
      default:
        throw Exceptions( Exceptions::DOMAIN );
    }

    if( options.OUTPUT.GRID_POLICY == Options::OUT_UNIQUE_PHYSICAL
     && ( dom.type == FFDom::LGL || dom.type == FFDom::CGL )
     && iel > 0 && inode == 0 )
      return false;

    return true;
  };

  auto count_grid_rows = [&]( t_EqnDom const& grid, FFVar const& fctvar ) -> size_t
  {
    size_t cnt = 1;
    for( auto const& [var, lim] : grid ){
      auto itdom = mDom.find( var );
      if( itdom == mDom.end() )
        throw Exceptions( Exceptions::DOMAIN );
      auto const& dom = itdom->second;

      size_t nsel = 0;
      // marching (2026-10-06): the evolution direction is collapsed to the resident window (one element), but the
      // output's rows span the WHOLE horizon -- count them over the global windows, with the global element index,
      // as the per-window emission applies the masks (see _window_output_rows)
      auto const ew = _evolution_window();
      bool const evo_march = _evolution_dom_var.dag() && var.id() == _evolution_dom_var.id() && ew.n_elements >= 1;
      size_t const nel_eff = evo_march? ew.n_elements: dom.n_elem;
      switch( lim ){
        case FFDom::ALL:
        case FFDom::ALL-FFDom::LB:
        case FFDom::ALL-FFDom::UB:
        case FFDom::ALL-FFDom::LB-FFDom::UB:
          for( size_t iel = 0; iel < nel_eff; ++iel )
            for( size_t inode = 0; inode < dom.n_node; ++inode )
              if( keep_output_node( dom, lim, iel, inode, nel_eff ) ) ++nsel;
          break;
        case FFDom::LB:
        case FFDom::UB:
          nsel = 1;
          break;
        default:
          std::cerr << "  **ERROR: MISSPECIFIED BOUNDARY/DOMAIN " << var
                    << " IN OUTPUT " << fctvar << std::endl;
          throw Exceptions( Exceptions::DOMAIN );
      }
      cnt *= nsel;
    }
    return cnt;
  };

  if( nCollFct ) *nCollFct = 0;
  size_t rowf = 0;
  for( auto const& fct : mFct ){
    try{
      fct.row0 = rowf;
      fct.nrow = ( fct.kind == FctKind::DISTRIBUTED )? count_grid_rows( fct.grid, fct.var ) : size_t(1);
    }
    catch(...){
      return false;
    }
    rowf += fct.nrow;
    if( nCollFct ) *nCollFct += fct.nrow;
  }

  return true;
}

// ===========================================================================
// rev142c: SPQR's own accounting of the factorization, which we were already paying for
// and then discarding.
//
// WHY.  PDE5 makes IC_STRONG cost 1.87x IC_WEAK at the SAME n (50.2 s against 26.9 s at
// n~24000), while IC_TRACE -- 13.5% LARGER than STRONG -- costs 0.95x.  So the driver of
// audit cost is not the column count.  The leading explanation is fill: STRONG is the one
// mode that ELIMINATES trace variables (3236 of them on PDE5 mono, 36 kept explicit)
// rather than carrying them, and substituting a tau's defining relation couples the two
// sides of an interface directly, which lengthens paths in the graph of A'A.  TRACE
// keeping all 3272 taus explicit is the same argument in reverse: trivially sparse border
// rows that the ordering handles well, hence a LARGER system that factorizes faster.
//
// That is a hypothesis read off the logs, and this project has been burned repeatedly by
// exactly that -- most recently by a 250.9 s baseline that stood unchallenged for three
// revisions and did not reproduce.  SPQR already computes the number that settles it.
// nnz(R) at equal n is the fill; if STRONG's is ~2x WEAK's, the explanation holds, and if
// it is not, the cost is somewhere else entirely and the IC_STRONG plan is not the lever.
//
// istat indices verified against SPQR/Source/spqr_1factor.cpp (SuiteSparse 7.6.1):
// [0] nnz(R) (accumulated within one call; our cholmod_common is started fresh per call,
// so it cannot accumulate across calls), [4] rank estimate, [5] column singletons,
// [6] singleton rows, [7] ordering used.  The analyze/factorize timers are assigned
// unconditionally there, not behind a TIMING build flag.
inline void
OCFEBASE::_spqr_stat_report
( t_SpqrStat const& S, char const* tag )
{
  static bool const kStat = Options::_env_flag( "CRONOS_SPQR_STAT", false )
                         || Options::_env_flag( "CRONOS_DEBUG_RANK", false );
  if( !kStat ) return;
  auto ordname = []( int o ) -> char const* {
    switch( o ){
      case 0: return "FIXED";   case 1: return "NATURAL"; case 2: return "COLAMD";
      case 3: return "GIVEN";   case 4: return "CHOLMOD"; case 5: return "AMD";
      case 6: return "METIS";   case 7: return "DEFAULT"; case 8: return "BEST";
      case 9: return "BESTAMD"; default: return "?"; } };
  std::ostream& os = std::cerr;
  std::ios_base::fmtflags const fl = os.flags();
  std::streamsize const pr = os.precision();
  os << "  [spqr-stat] " << ( tag ? tag : "" )
     << " n=" << S.nrow << "x" << S.ncol << " rank=" << S.rank
     << " ord=" << ordname( S.ordering ) << "(" << S.ordering << ")"
     << " nnzA=" << S.nnzA << " nnzR=" << S.nnzR;
  os << std::fixed << std::setprecision(1);
  if( S.ncol ) os << " nnzR/n=" << (double)S.nnzR / (double)S.ncol;
  if( S.nnzA ) os << " fill=" << (double)S.nnzR / (double)S.nnzA << "x";
  os << std::scientific << std::setprecision(3) << " flops=" << S.flops;
  os << std::fixed << std::setprecision(2)
     << " analyze=" << S.t_analyze << "s factorize=" << S.t_factorize << "s";
  if( S.ncolsing || S.nrowsing )
    os << " singletons(col/row)=" << S.ncolsing << "/" << S.nrowsing;
  os << std::endl;
  os.flags( fl ); os.precision( pr );
}

// rev142d: the keep-explicit census, joined to the tau-column occupancy.
//
// _keepExplicitReason records WHY each IC_STRONG claim was kept as an explicit tau
// rather than eliminated.  None of the three routes is "no pivot was available" -- an
// earlier census said that and it was wrong.  Route 2 is the one that makes the column
// load-bearing: every replacement candidate is an auxiliary LINK row, so exact-replacing
// would consume that auxiliary's defining equation and turn the AUXILIARY into a null
// direction.  Eliminating the block would not remove a null direction, it would move one
// from the tau block (benign gauge, GAUGE_ONLY) into a state (STATE_CONTENT).  That is
// why the column stays.
//
// What is NOT deliberate is the column arriving without its row.  rev41's corpus sweep
// pinned n_explicit == k and rank(B) + n_explicit == n_trace_var across 78 rows, 28
// drivers and five model families -- arithmetic that assumes the rows exist.  This
// routine is the first thing that can tell the two apart per column.
inline void
OCFEBASE::_keep_explicit_audit
( arma::sp_mat const& T, size_t tau0, size_t rank_T ) const
{
  static bool const kAudit = Options::_env_flag( "CRONOS_KEEP_EXPLICIT_AUDIT", false )
                          || Options::_env_flag( "CRONOS_DEBUG_RANK", false );
  if( !kAudit ) return;

  T.sync();
  size_t nempty = 0, nzmin = std::numeric_limits<size_t>::max(), nzmax = 0;
  size_t first_empty = std::numeric_limits<size_t>::max();
  for( arma::uword j = 0; j < T.n_cols; ++j ){
    size_t const nz = (size_t)( T.col_ptrs[j+1] - T.col_ptrs[j] );
    if( !nz ){ ++nempty; if( first_empty == std::numeric_limits<size_t>::max() ) first_empty = tau0 + (size_t)j; }
    nzmin = std::min( nzmin, nz );
    nzmax = std::max( nzmax, nz );
  }
  if( !T.n_cols ) nzmin = 0;

  size_t r1 = 0, r2 = 0, r3 = 0, rx = 0;
  for( auto const& kv : _keepExplicitReason )
    switch( kv.second ){ case 1: ++r1; break; case 2: ++r2; break; case 3: ++r3; break; default: ++rx; }

  std::ostream& os = std::cerr;
  os << "  [keep-explicit] tau block " << T.n_rows << "x" << T.n_cols
     << " tau0=" << tau0 << " rank(T)=" << rank_T
     << " | cols empty=" << nempty << "/" << T.n_cols
     << " nnz_min=" << nzmin << " nnz_max=" << nzmax;
  if( nempty ) os << " first_empty=" << first_empty;
  os << " | reasons{1:verify+rank-ext=" << r1 << " 2:aux-LINK=" << r2
     << " 3:verify=" << r3 << " other=" << rx
     << " total=" << _keepExplicitReason.size() << "}";
  // The two failure signatures, named rather than left to the reader.
  if( nempty )
    os << "  ** " << nempty << " COLUMN(S) WITHOUT A ROW -- each is one gauge direction **";
  if( _keepExplicitReason.size() != T.n_cols && !_keepExplicitReason.empty() )
    os << "  ** census " << _keepExplicitReason.size() << " != n_tau " << T.n_cols
       << ": the tau block holds more than the kept-explicit claims **";
  os << std::endl;
}

inline void
OCFEBASE::_setup_node_cache
()
{
  _nodeBlockCache.clear();
  for( auto const& [var, dom] : _mDom ){
    auto& blocks = _nodeBlockCache[var];
    blocks.resize( dom.n_elem );
    for( size_t iel = 0; iel < dom.n_elem; ++iel )
      blocks[iel] = dom.lgnodes( dom.elem_lo( iel ), dom.elem_up( iel ) );
  }
}

template <typename U>
inline void
OCFEBASE::_reset_colloc
( std::vector<OCVar<U>>& OCvar,  std::vector<OCVar<U>>& OCinp,
  std::vector<OCVar<U>>& OCdom, std::vector<OCVar<U>>& OCcst )
const
{
  OCvar.resize( _mVar.size() );
  size_t ndx = 0;
  for( auto const& [var, dom] : _mVar )
    OCvar[ndx++].set( *this, dom, nullptr ); // no data passed yet

  OCinp.resize( _mInp.size() );
  ndx = 0;
  for( auto const& [var, dom] : _mInp ){
    auto const idom = _input_ocdom<U>( var, dom );
    OCinp[ndx++].set( *this, idom, nullptr ); // no data passed yet
  }

  OCcst.resize( _vCst.size() );
  for( size_t ndx = 0; ndx < _vCst.size(); )
    OCcst[ndx++].set( *this, std::set<FFVar,lt_FFVar>{}, nullptr ); // no data passed yet

  OCdom.resize( _mDom.size() );
  ndx = 0;
  for( auto& [var, dom] : _mDom )
    OCdom[ndx++].set( *this, std::set<FFVar,lt_FFVar>{var}, nullptr ); // no data passed yet
}

// ======================================================================
// OCFESLV::_differential_core_shape        (rev94, USER-side, pre-reduction)
// ======================================================================
inline void
OCFEBASE::_differential_core_shape
( size_t& rows, size_t& cols )
const
{
  rows = 0;
  cols = 0;
  if( !_evolution_dom_var.dag() ) return;

  std::set<FFVar,lt_FFVar> states;
  for( auto const& eqn : _usr._mEqnUsr ){
    FFGraph* dag = static_cast<FFGraph*>( eqn.var.dag() );
    if( !dag ) continue;
    FFVar expr = eqn.var;
    FFSubgraph sg = dag->subgraph( 1, &expr );
    bool row_has_dt = false;
    for( auto const* op : sg.l_op ){
      FFVar const st = _time_derivative_state( op, _usr._mVarUsr );
      if( st.dag() ){ row_has_dt = true; states.insert( st ); }
    }
    if( row_has_dt ) ++rows;      // one ROW per equation, however many derivatives it holds
  }
  cols = states.size();
}

// ===========================================================================
//
// Question (agreed 2026-08-28): if the drop decision were ONE computation --
// maximal independent subset of the tau columns (the full pass-0 receiver
// block B, before any eligibility gate), dependent complement as candidates,
// each candidate's continuity row tested against the null space of the reduced
// system -- what changes against the shipped composite (gated MGS + STRONG
// keep-explicit protocol)?
//
// Mechanics (no shadow instance, no re-setup: a deep-copied or re-set-up env
// segfaults at classify on blk0-class models -- filed as an open item):
//  * _detect_redundant_continuity_claims stashes its pass-0 (col, col2claim)
//    maps -- the full symbolic B over receiver rows for EVERY claim, identical
//    in both modes -- before the gates run.
//  * At the end of the shipped setup, the oracle re-runs the SAME sparse MGS on
//    the stash with NO gates -> the uniform candidate set.
//  * W-test: on IC_TRACE, the shipped assembly minus the still-kept candidates'
//    rows and tau columns is ranked; its null vectors score every candidate's
//    reconstructed continuity row; pivoted QR of W names the needed subset.
//    On IC_STRONG the kept candidates are inside the Schur and cannot be
//    excised by row/col surgery, so the W-test covers the shipped-DROPPED
//    subset only (the [drop-null] test) and the kept candidates are reported
//    UNTESTED -- resolving them under STRONG needs the TRACE twin.
inline bool
OCFEBASE::_null_basis
( arma::sp_mat const& M, double tol_abs, size_t& rank_out, arma::mat& N )
{
  rank_out = 0; N.reset();
  size_t const nr = M.n_rows, nc = M.n_cols;
  if( !nr || !nc ) return false;
  // The factorization M P = Q R by SPQR or, without it, by Eigen's sparse QR (2026-10-08, WORKPLAN 4R.1: the EPL build
  // has no SPQR, so above the dense cap the W-test had no null basis and a needed dropped claim stayed dropped).  Either
  // way R (CSC) and the column permutation go into the arrays below, and ONE construction builds the basis from them.
  std::vector<int64_t> Rp, Ri, Ep;  std::vector<double> Rx;  size_t ncolR = 0;  bool sortedR = false;  int64_t rk = -1;
#if defined(CRONOS__WITH_SPQR)
  { // ---- SPQR ----
    M.sync();
    std::vector<int64_t> cp( nc + 1 ), ri( M.n_nonzero ? M.n_nonzero : 1 );
    for( size_t j = 0; j <= nc; ++j ) cp[j] = (int64_t)M.col_ptrs[j];
    for( size_t k = 0; k < M.n_nonzero; ++k ) ri[k] = (int64_t)M.row_indices[k];
    cholmod_common cc;  cholmod_l_start( &cc );
    // rev142b: CHOLMOD's own OpenMP regions honour this; guarded on the CHOLMOD major
    // version because nthreads_max does not exist before CHOLMOD 4 (SuiteSparse 6).
  #if defined(CHOLMOD_MAIN_VERSION) && (CHOLMOD_MAIN_VERSION >= 4)
    if( int const _cap = t_ThreadCap::active() ) cc.nthreads_max = _cap;
  #endif
    cholmod_sparse Ac;  std::memset( &Ac, 0, sizeof(Ac) );
    Ac.nrow = nr; Ac.ncol = nc; Ac.nzmax = M.n_nonzero;
    Ac.p = cp.data(); Ac.i = ri.data(); Ac.x = const_cast<double*>( M.values );
    Ac.stype = 0; Ac.itype = CHOLMOD_LONG; Ac.xtype = CHOLMOD_REAL; Ac.dtype = CHOLMOD_DOUBLE;
    Ac.sorted = 1; Ac.packed = 1;
    cholmod_sparse* Rc = nullptr; int64_t* E = nullptr;
    // rev142: ORDERING_BEST is the default for the AUDIT factorizations (this call and
    // numeric_rank's).  It lets SPQR try both its orderings and keep the sparser R.
    // MEASURED 11-12.5% off the factorization at n ~ 3k with the margin growing in n;
    // VERDICT-INVARIANT over the corpus -- 1159 instances, det_check 0 deviations and an
    // empty defR/defL/k' diff, which is also datum (3) for the refinement (see
    // _refine_rank_from_R).  CRONOS_SPQR_ORDERING restores any SPQR_ORDERING_* value.
    // The SOLVE path (SuiteSparseQR_min2norm) is deliberately NOT changed here: it
    // returns a solution, not a rank, so a different pivot order moves the numbers in
    // the behavioural signature and would need its own sweep.
    static int const kSpqrOrd = SPQR_ORDERING_BEST;   // CRONOS_SPQR_ORDERING (retired 2026-10-07, WORKPLAN 3.B batch 2b)
    rk = SuiteSparseQR<double>( kSpqrOrd, tol_abs,
                                              (int64_t)nr, &Ac, &Rc, &E, &cc );
    // rev142c: same statistics on the w-decide null-basis factorization.
    { t_SpqrStat S;
      S.nrow = nr; S.ncol = nc; S.nnzA = M.n_nonzero; S.rank = ( rk >= 0 ? (size_t)rk : 0 );
      S.nnzR = cc.SPQR_istat[0]; S.ncolsing = cc.SPQR_istat[5]; S.nrowsing = cc.SPQR_istat[6];
      S.ordering = (int)cc.SPQR_istat[7]; S.flops = cc.SPQR_flopcount;
      S.t_analyze = cc.SPQR_analyze_time; S.t_factorize = cc.SPQR_factorize_time;
      _spqr_stat_report( S, "null" ); }
    if( rk >= 0 && Rc && Rc->x && Rc->p && Rc->i ){
      ncolR = (size_t)Rc->ncol;  sortedR = ( Rc->sorted != 0 );
      int64_t const* p = static_cast<int64_t const*>( Rc->p );
      Rp.assign( p, p + ncolR + 1 );
      Ri.assign( static_cast<int64_t const*>( Rc->i ), static_cast<int64_t const*>( Rc->i ) + p[ncolR] );
      Rx.assign( static_cast<double const*>( Rc->x ), static_cast<double const*>( Rc->x ) + p[ncolR] );
      if( E ) Ep.assign( E, E + nc );
    }
    if( Rc ) cholmod_l_free_sparse( &Rc, &cc );
    if( E )  cholmod_l_free( nc, sizeof(int64_t), E, &cc );
    cholmod_l_finish( &cc );
  }
#elif defined(CRONOS__WITH_EIGEN)
  { // ---- Eigen's sparse QR (no SPQR: the EPL build) ----
    typedef Eigen::SparseMatrix<double, Eigen::ColMajor, int> t_Sp;
    std::vector<Eigen::Triplet<double>> trip;  trip.reserve( M.n_nonzero );
    for( auto it = M.begin(); it != M.end(); ++it ) trip.emplace_back( (int)it.row(), (int)it.col(), *it );
    t_Sp A( (int)nr, (int)nc );  A.setFromTriplets( trip.begin(), trip.end() );  A.makeCompressed();
    Eigen::SparseQR<t_Sp, Eigen::COLAMDOrdering<int>> qr;
    qr.setPivotThreshold( tol_abs );     // absolute, as SPQR's tol
    qr.compute( A );
    if( qr.info() == Eigen::Success ){
      rk = (int64_t)qr.rank();
      t_Sp R = qr.matrixR();  R.makeCompressed();
      ncolR = (size_t)R.cols();
      Rp.assign( R.outerIndexPtr(), R.outerIndexPtr() + ncolR + 1 );
      Ri.assign( R.innerIndexPtr(), R.innerIndexPtr() + R.nonZeros() );
      Rx.assign( R.valuePtr(), R.valuePtr() + R.nonZeros() );
      // Eigen does not sort a column's entries by row: sort them here (linear overall), so that the staircase uses its
      // row-sorted walk -- an unsorted R took the old column scan and printed a diagnostic on every EPL run.
      { std::vector<std::pair<int64_t,double>> col;
        for( size_t j = 0; j < ncolR; ++j ){
          col.clear();
          for( int64_t k = Rp[j]; k < Rp[j+1]; ++k ) col.emplace_back( Ri[k], Rx[k] );
          std::sort( col.begin(), col.end(), []( auto const& x, auto const& y ){ return x.first < y.first; } );
          for( int64_t k = Rp[j], q = 0; k < Rp[j+1]; ++k, ++q ){ Ri[k] = col[q].first; Rx[k] = col[q].second; }
        }
        sortedR = true; }
      // E[k] = the original column of R's k-th column: the permutation applied to 0..nc-1 (convention-proof)
      Eigen::RowVectorXd v( (Eigen::Index)nc );  for( size_t j = 0; j < nc; ++j ) v( (Eigen::Index)j ) = (double)j;
      Eigen::RowVectorXd const vp = v * qr.colsPermutation();
      Ep.resize( nc );  for( size_t j = 0; j < nc; ++j ) Ep[j] = (int64_t)std::llround( vp( (Eigen::Index)j ) );
    }
  }
#else
  return false;
#endif
  int64_t const* E = Ep.empty()? nullptr: Ep.data();
  bool ok = false;
  do{
    if( rk < 0 || Rp.empty() ) break;
    size_t const r = (size_t)rk;
    int64_t const* rp = Rp.data();
    int64_t const* rI = Ri.data();
    double  const* rx = Rx.data();
    // staircase: live columns of the squeezed R (rev142: one shared walk, O(1)/column
    // when Rc is row-sorted; identical set by construction -- see the routine).
    std::vector<size_t> live;
    _staircase_live_columns( rp, rI, rx, ncolR, r, sortedR, live );
    if( live.size() != r ) break;
    std::vector<char> is_live( ncolR, 0 ); for( size_t j : live ) is_live[j] = 1;
    std::vector<size_t> dead; for( size_t j = 0; j < ncolR; ++j ) if( !is_live[j] ) dead.push_back( j );
    // R11 (r x r, CSC, upper triangular by staircase) and per-dead-column rhs
    std::vector<std::vector<std::pair<size_t,double>>> colR( r );
    for( size_t jj = 0; jj < r; ++jj ){ size_t const j = live[jj];
      for( int64_t k = rp[j]; k < rp[j+1]; ++k ) if( (size_t)rI[k] < r )
        colR[jj].push_back( { (size_t)rI[k], rx[k] } ); }
    std::vector<double> dg( r, 0. );
    for( size_t jj = 0; jj < r; ++jj ) for( auto const& [i2,v] : colR[jj] ) if( i2 == jj ) dg[jj] = v;
    for( size_t jj = 0; jj < r; ++jj ) if( dg[jj] == 0. ) { dg[jj] = 64.*DBL_EPSILON; }
    auto usolve = [&]( std::vector<double>& b ){            // R11 x = b (in-place)
      for( size_t jj = r; jj-- > 0; ){
        double const xj = b[jj] / dg[jj]; b[jj] = xj;
        for( auto const& [i2,v] : colR[jj] ) if( i2 < jj ) b[i2] -= v * xj; } };
    auto lsolve = [&]( std::vector<double>& c ){            // R11^T y = c (in-place)
      for( size_t jj = 0; jj < r; ++jj ){
        double acc = c[jj];
        for( auto const& [i2,v] : colR[jj] ) if( i2 < jj ) acc -= v * c[i2];
        c[jj] = acc / dg[jj]; } };
    // refined rank (deflation) via the shared routine, on freshly built sp_mats
    size_t rank_ref = r; double gap2 = 0., scut2 = 0.;
    { arma::sp_mat R11s, R12s;
      { size_t cnt=0; for( size_t jj=0; jj<r; ++jj ) cnt += colR[jj].size();
        arma::umat loc(2,(arma::uword)std::max<size_t>(cnt,1)); arma::vec val((arma::uword)std::max<size_t>(cnt,1));
        arma::uword e2=0;
        for( size_t jj=0; jj<r; ++jj ) for( auto const& [i2,v] : colR[jj] ){ loc(0,e2)=(arma::uword)i2; loc(1,e2)=(arma::uword)jj; val(e2)=v; ++e2; }
        R11s = e2 ? arma::sp_mat( loc.cols(0,e2-1), val.subvec(0,e2-1), (arma::uword)r, (arma::uword)r )
                  : arma::sp_mat( (arma::uword)r, (arma::uword)r ); }
      { size_t cnt=0; for( size_t j : dead ) for( int64_t k=rp[j]; k<rp[j+1]; ++k ) if( (size_t)rI[k]<r ) ++cnt;
        arma::umat loc(2,(arma::uword)std::max<size_t>(cnt,1)); arma::vec val((arma::uword)std::max<size_t>(cnt,1));
        arma::uword e2=0;
        for( size_t jj=0; jj<dead.size(); ++jj ){ size_t const j=dead[jj];
          for( int64_t k=rp[j]; k<rp[j+1]; ++k ) if( (size_t)rI[k]<r ){ loc(0,e2)=(arma::uword)rI[k]; loc(1,e2)=(arma::uword)jj; val(e2)=rx[k]; ++e2; } }
        R12s = e2 ? arma::sp_mat( loc.cols(0,e2-1), val.subvec(0,e2-1), (arma::uword)r, (arma::uword)dead.size() )
                  : arma::sp_mat( (arma::uword)r, (arma::uword)dead.size() ); }
      _refine_rank_from_R( R11s, R12s, tol_abs, rank_ref, gap2, scut2 );
    }
    size_t const ndemote = ( rank_ref < r ? r - rank_ref : 0 );
    size_t const kdim = dead.size() + ndemote;
    if( !kdim ){ rank_out = rank_ref; ok = true; break; }     // full column rank
    N.zeros( (arma::uword)nc, (arma::uword)kdim );
    // (a) exact dead-column vectors
    std::vector<double> b( r );
    for( size_t dj = 0; dj < dead.size(); ++dj ){
      std::fill( b.begin(), b.end(), 0. );
      size_t const j = dead[dj];
      for( int64_t k = rp[j]; k < rp[j+1]; ++k ) if( (size_t)rI[k] < r ) b[ (size_t)rI[k] ] = -rx[k];
      usolve( b );
      double nrm = 1.0; for( size_t jj = 0; jj < r; ++jj ) nrm += b[jj]*b[jj];
      nrm = std::sqrt( nrm );
      for( size_t jj = 0; jj < r; ++jj ) N( (arma::uword)( E ? (size_t)E[ live[jj] ] : live[jj] ), (arma::uword)dj ) = b[jj] / nrm;
      N( (arma::uword)( E ? (size_t)E[j] : j ), (arma::uword)dj ) = 1.0 / nrm;
    }
    // (b) demotion vectors: deflated inverse iteration on R11^T R11
    for( size_t q2 = 0; q2 < ndemote; ++q2 ){
      std::vector<double> w( r );
      for( size_t jj = 0; jj < r; ++jj ) w[jj] = std::sin( 0.7*(double)(jj+1)*(double)(q2+1) ) + 1e-3;
      for( int it = 0; it < 8; ++it ){
        lsolve( w ); usolve( w );
        for( size_t p2 = 0; p2 < q2; ++p2 ){                 // deflate vs previous demotions
          double dp = 0.;
          for( size_t jj = 0; jj < r; ++jj ) dp += w[jj] * N( (arma::uword)( E ? (size_t)E[ live[jj] ] : live[jj] ), (arma::uword)( dead.size()+p2 ) );
          for( size_t jj = 0; jj < r; ++jj ) w[jj] -= dp * N( (arma::uword)( E ? (size_t)E[ live[jj] ] : live[jj] ), (arma::uword)( dead.size()+p2 ) );
        }
        double nrm = 0.; for( double v : w ) nrm += v*v; nrm = std::sqrt( std::max( nrm, 1e-300 ) );
        for( double& v : w ) v /= nrm;
      }
      for( size_t jj = 0; jj < r; ++jj ) N( (arma::uword)( E ? (size_t)E[ live[jj] ] : live[jj] ), (arma::uword)( dead.size()+q2 ) ) = w[jj];
    }
    rank_out = rank_ref; ok = true;
  } while( false );
  return ok;
}

// rev133: post-QR rank refinement.  A weak rank-revealing QR (SPQR, Eigen) declares
// a column dead only when its REMAINING norm falls below the threshold; an exact
// dependency read through accumulated round-off can sit just ABOVE it (MEASURED:
// PDE6_blk0 IC_TRACE MANUAL, equilibrated, the 4th dead column retained at 1.1x the
// threshold while the dense SVD shows a 9.7e11 gap).  The remedy is a condition
// estimate on the retained-row block R = [R11 R12] -- NOT on R11 alone: a dead
// column keeps its entries above its pivot (its projection onto the earlier
// Householder directions), so R12 != 0 and sigma_min(R11) can sit far below
// sigma_r(A) (MEASURED: PDE6_blk0 IC_STRONG, sigma_min(R11) = 4.0e-07 against
// sigma_1246(A) = 3.6e-02).  sigma_r(A) = sigma_min(R) exactly (Q is orthogonal).
//
// METHOD.  (1) Block inverse iteration on R11 R11^T (two sparse triangular solves
// per vector, backward-stable) generates the p smallest LEFT singular directions of
// R11.  (2) Rayleigh-Ritz of the FULL operator R R^T on that block: Y = R^T V,
// G = Y^T Y, sigma = sqrt(eig(G)).  Ritz values are upper bounds on the smallest
// sigma(R)^2, so every deflation is safe (a Ritz value below the threshold proves a
// singular value below it); a near-null direction of R is a near-null direction of
// R11^T up to tol/sigma_{p+1}(R11), so with p > #missed dependencies the bound is
// tight.  A Woodbury application of (R R^T)^{-1} was tried first and rejected:
// Z = R11^{-1} R12 is O(1e11) exactly when R11 carries a missed dependency, and the
// cancellation produced Ritz values BELOW the true sigma_2 (1.8e-06 vs 9.1e-04).
// Cost: 2p sparse triangular solves per iteration, 30 iterations; skipped (QR rank
// kept) when d > 256.
//
// rev142 -- WHY THIS PATH IS LOAD-BEARING, AND WHY IT MUST NOT BE TIERED AROUND.
// The name says "refinement", which reads cosmetic.  It is not: it is what makes the
// rank CORRECT.  Recorded here because two optimisation schemes were proposed against
// this routine and both were refuted on exactly this ground, and because nothing in
// the code says so.  Three independent lines, each measured:
//
//  (1) CORRECT AGAINST GROUND TRUTH.  62 instances cross-checked against a dense SVD
//      of the same matrix at the same tol_abs: zero disagreements.  PDE16 n = 360 --
//      SPQR 358, refined 355, dense SVD 355, with sigma[r-1] = 9.4e-04 against
//      sigma[r] = 5.1e-16 and tolabs = 7.4e-12.  A TWELVE-order gap, and SPQR still
//      counted three columns too many; the QR rank is not a near-miss on these
//      systems.  PDE6_blk0 n = 1522 -- SPQR 1519, refined 1518, SVD 1518.  The oracle
//      that produced these numbers ships as _rank_oracle_check (default off).
//  (2) MACHINE-INDEPENDENT.  PDE6_blk0 n = 1250 gives SPQR 1246 in the sandbox and
//      1247 on the workstation; the refinement returns 1246 on BOTH.  The corpus is
//      cross-machine reproducible BECAUSE of this routine, not despite it -- the raw
//      QR rank is not reproducible across BLAS/SuiteSparse builds.
//  (3) PIVOTING-PATH-INDEPENDENT.  Under SPQR_ORDERING_BEST (the rev142 default) the
//      whole corpus returns identical deficiencies to the ORDERING_DEFAULT run: 1159
//      instances, empty defR/defL/k' diff.  A different elimination order reaches the
//      same refined rank.
//
// CONSEQUENCE: the refinement must never be weakened, made conditional on the QR
// reporting a deficiency, tiered behind a cheap pre-check, or skipped on a size
// heuristic.  In particular "SPQR says full rank, so there is nothing to refine" is
// FALSE -- MEASURED: n = 1297, spqr = 1297 (full), refined = 1291.  That single result
// also killed coarse-to-fine plan transfer and bordered-LU confirmation, both of which
// needed a cheap trustworthy rank on the fine system.  If a cheaper path is ever
// proposed, it must be validated against (1)-(3) on the corpus BEFORE it ships, not
// argued from the structure of the code.
inline void
OCFEBASE::_refine_rank_from_R
( arma::sp_mat const& R11, arma::sp_mat const& R12, double tol_abs,
  size_t& rank, double& gap, double& sigma_cut )
{
  size_t const r = R11.n_rows;
  size_t const d = R12.n_cols;
  if( !r || R11.n_cols != r || ( d && R12.n_rows != r ) || d > 256 ){ return; }
  R11.sync(); R12.sync();
  arma::uword const* cp = R11.col_ptrs;
  arma::uword const* ri = R11.row_indices;
  double      const* rx = R11.values;

  std::vector<double> diag( r, 0. );
  for( size_t j = 0; j < r; ++j )
    for( arma::uword k = cp[j]; k < cp[j+1]; ++k ) if( ri[k] == j ){ diag[j] = rx[k]; break; }
  for( size_t j = 0; j < r; ++j ) if( diag[j] == 0. ){ return; }   // not a live triangular block: keep the QR rank

  // upper-triangular solve  R11 x = b   (CSC, column-oriented backward substitution)
  auto solve_R = [&]( std::vector<double>& b ){
    for( size_t jj = r; jj-- > 0; ){
      double const xj = b[jj] / diag[jj];
      b[jj] = xj;
      for( arma::uword k = cp[jj]; k < cp[jj+1]; ++k ){ size_t const i = ri[k]; if( i < jj ) b[i] -= rx[k] * xj; }
    }
  };
  // lower-triangular solve  R11^T y = c   (row j of R11^T is column j of R11)
  auto solve_Rt = [&]( std::vector<double>& c ){
    for( size_t j = 0; j < r; ++j ){
      double acc = c[j];
      for( arma::uword k = cp[j]; k < cp[j+1]; ++k ){ size_t const i = ri[k]; if( i < j ) acc -= rx[k] * c[i]; }
      c[j] = acc / diag[j];
    }
  };

  size_t const p = std::min<size_t>( 8, r );
  arma::mat V( (arma::uword)r, (arma::uword)p );
  for( size_t q = 0; q < p; ++q )                       // deterministic, non-degenerate start
    for( size_t i = 0; i < r; ++i ) V( i, q ) = std::sin( 0.37 * (double)( i + 1 ) * (double)( q + 1 ) + 0.11 * (double)q );
  auto orthonormalise = [&](){ arma::mat Q, Rq; if( arma::qr_econ( Q, Rq, V ) ) V = Q; };
  orthonormalise();
  std::vector<double> w( r );
  for( int iter = 0; iter < 30; ++iter ){
    for( size_t q = 0; q < p; ++q ){                    // v <- (R11 R11^T)^{-1} v = R11^{-T} R11^{-1} v
      for( size_t i = 0; i < r; ++i ) w[i] = V( i, q );
      solve_R( w ); solve_Rt( w );
      for( size_t i = 0; i < r; ++i ) V( i, q ) = w[i];
    }
    orthonormalise();
  }
  // Rayleigh-Ritz of R R^T on span(V):  Y = R^T V ((r+d) x p),  G = Y^T Y.
  arma::mat Y( (arma::uword)( r + d ), (arma::uword)p, arma::fill::zeros );
  for( size_t q = 0; q < p; ++q ){
    for( size_t j = 0; j < r; ++j ){                 // (R11^T v)_j = sum_i R11(i,j) v_i
      double acc = 0.;
      for( arma::uword k = cp[j]; k < cp[j+1]; ++k ) acc += rx[k] * V( ri[k], q );
      Y( j, q ) = acc;
    }
    for( size_t c2 = 0; c2 < d; ++c2 ){               // (R12^T v)_c
      double acc = 0.;
      for( arma::sp_mat::const_col_iterator it = R12.begin_col( c2 ); it != R12.end_col( c2 ); ++it ) acc += (*it) * V( it.row(), q );
      Y( r + c2, q ) = acc;
    }
  }
  arma::mat const G = Y.t() * Y;
  arma::vec ev;
  if( !arma::eig_sym( ev, G ) ) return;
  std::vector<double> sig( p, 0. );
  for( size_t q = 0; q < p; ++q ) sig[q] = std::sqrt( std::max( ev(q), 0. ) );
  std::sort( sig.begin(), sig.end() );
  size_t ndeflate = 0;
  for( size_t q = 0; q < p; ++q ) if( sig[q] < tol_abs ) ++ndeflate;
  rank = ( rank > ndeflate ? rank - ndeflate : 0 );
  double const kept = ( ndeflate < p ) ? sig[ndeflate] : tol_abs;
  double const cut  = ( ndeflate > 0 ) ? std::max( sig[ndeflate-1], 0. ) : 0.;
  sigma_cut = kept;
  gap = ( ndeflate > 0 ) ? ( cut > 0. ? kept / cut : std::numeric_limits<double>::infinity() )
                         : ( tol_abs > 0. ? kept / tol_abs : std::numeric_limits<double>::infinity() );
}

// rev142: live columns of a squeezed SPQR R, O(1) per column instead of O(nnz).
//
// SPQR returns R as rank-by-n with the dead columns left IN PLACE (not permuted to
// the end), so the live columns are found by walking the staircase: column j is live
// iff it has a nonzero in row `next_row`, the next unfilled pivot row.  Through rev141
// that test scanned the WHOLE column to find the entry.
//
// The scan is unnecessary.  R is upper triangular with respect to the staircase: when
// column j is examined, every pivot row < next_row is filled and every row >= next_row
// is empty, so a live column's entry at next_row is its LAST entry and a dead column's
// entries all sit strictly above it.  With row-sorted CSC columns (cholmod's `sorted`
// flag) the last stored entry of the column decides it, and the two paths agree by
// construction -- including the degenerate cases: an explicitly stored zero at
// next_row reads dead under both (the scan applies the same rx != 0 test), and an
// empty column reads dead under both.
//
// MEASURED (sandbox, PDE5 n = 2688, one audit rank pass): the scan was 239 ms of a
// 613 ms pass -- the same order as the SPQR factorization itself (239 ms) and 39% of
// the pass.  It is quadratic in the wrong quantity: every dead column is scanned to
// its end, and the dense-ish trailing columns are exactly the expensive ones.
//
// The rev141 scan is KEPT as the unsorted fallback rather than deleted: nothing in the
// SuiteSparseQR contract promises sorted output, and a silently wrong live set would
// not announce itself -- it would show up as a moved rank three gates later.
// CRONOS_STAIRCASE_VERIFY=1 runs BOTH on every call and prints one line each way
// (AGREE or MISMATCH, never silence), returning the SCAN result so a verify run cannot
// itself change a verdict.
inline void
OCFEBASE::_staircase_live_columns
( int64_t const* rp, int64_t const* rI, double const* rx,
  size_t ncol, size_t r, bool sorted, std::vector<size_t>& live )
{
  live.clear();
  if( !rp || !rI || !rx || !r || !ncol ) return;
  live.reserve( r );
  static bool const kVerify = Options::_env_flag( "CRONOS_STAIRCASE_VERIFY", false );

  // O(1) per column: last stored entry is the staircase entry of a live column.
  auto fast_scan = [&]( std::vector<size_t>& out ){
    out.clear(); out.reserve( r );
    size_t next_row = 0;
    for( size_t j = 0; j < ncol && next_row < r; ++j ){
      int64_t const k = rp[j+1] - 1;
      if( k >= rp[j] && (size_t)rI[k] == next_row && rx[k] != 0. ){ out.push_back( j ); ++next_row; }
    }
  };
  // rev141 path: scan the column for an entry in the next unfilled row.
  auto full_scan = [&]( std::vector<size_t>& out ){
    out.clear(); out.reserve( r );
    size_t next_row = 0;
    for( size_t j = 0; j < ncol && next_row < r; ++j ){
      bool hit = false;
      for( int64_t k = rp[j]; k < rp[j+1]; ++k )
        if( (size_t)rI[k] == next_row && rx[k] != 0. ){ hit = true; break; }
      if( hit ){ out.push_back( j ); ++next_row; }
    }
  };

  if( sorted && !kVerify ){ fast_scan( live ); return; }
  full_scan( live );

  // An unsorted R silently costs the whole rev142 saving, so say so ONCE per process
  // even without the verify flag: a runtime gate that fails for this reason must not
  // look like a mismeasurement.  If this line appears in a sweep the signature diff is
  // non-empty BY DESIGN -- it is the fast path reporting that it did not engage.
  if( !sorted ){
    static bool warned = false;
    if( !warned ){
      warned = true;
      std::cerr << "  [staircase] UNSORTED R from the QR backend: rev141 column scan in use"
                << " (ncol=" << ncol << " r=" << r << ")" << std::endl;
    }
    return;
  }
  if( !kVerify ) return;

  std::vector<size_t> fast; fast_scan( fast );
  std::ostream& os = std::cerr;
  std::ios_base::fmtflags const fl = os.flags();
  if( fast == live )
    os << "  [staircase] AGREE ncol=" << ncol << " r=" << r << " live=" << live.size() << std::endl;
  else{
    size_t d = 0, nmin = std::min( fast.size(), live.size() );
    while( d < nmin && fast[d] == live[d] ) ++d;
    os << "  [staircase] ** MISMATCH ** ncol=" << ncol << " r=" << r
       << " scan=" << live.size() << " fast=" << fast.size() << " first-diff@" << d;
    if( d < nmin ) os << " scan[" << d << "]=" << live[d] << " fast[" << d << "]=" << fast[d];
    os << std::endl;
  }
  os.flags( fl );
}

// rev142: dense-SVD rank oracle.  Report-only, default OFF (CRONOS_RANK_ORACLE=1),
// promoted from the sandbox scratch build that produced the 62-instance cross-check
// recorded above _refine_rank_from_R.  Factorizes the SAME matrix the backend saw, at
// the SAME tol_abs, and prints the two ranks side by side with the singular values
// straddling the cut -- the numbers needed to tell "the backend is right" from "the
// tolerance is in the wrong place".  Capped by CRONOS_RANK_ORACLE_MAXN (default 2000)
// because the dense factorization is O(n^3): at the cap it is seconds per call.
inline void
OCFEBASE::_rank_oracle_check
( arma::sp_mat const& M, double tol_abs, size_t rank_backend, size_t rank_refined,
  char const* backend )
{
  static bool   const kOracle = Options::_env_flag( "CRONOS_RANK_ORACLE", false );
  if( !kOracle ) return;
  static size_t const kMaxN = (size_t)Options::_env_dbl( "CRONOS_RANK_ORACLE_MAXN", 2000. );
  size_t const nr = M.n_rows, nc = M.n_cols;
  if( !nr || !nc ) return;
  std::ostream& os = std::cerr;
  std::ios_base::fmtflags const fl = os.flags();
  std::streamsize const pr = os.precision();
  if( std::max( nr, nc ) > kMaxN ){
    os << "  [rank-oracle] SKIPPED n=" << nr << "x" << nc
       << " over cap (CRONOS_RANK_ORACLE_MAXN=" << kMaxN << ")" << std::endl;
    os.flags( fl ); os.precision( pr );
    return;
  }
  arma::vec sv;
  if( !arma::svd( sv, arma::mat( M ) ) || sv.n_elem == 0 ){
    os << "  [rank-oracle] FAILED n=" << nr << "x" << nc << " (dense svd did not converge)" << std::endl;
    os.flags( fl ); os.precision( pr );
    return;
  }
  size_t rk = 0;
  for( arma::uword k = 0; k < sv.n_elem; ++k ) if( sv(k) > tol_abs ) ++rk;
  bool const agree = ( rk == rank_refined );
  os << "  [rank-oracle] " << ( agree ? "AGREE" : "** DISAGREE **" )
     << " backend=" << ( backend ? backend : "?" )
     << " n=" << nr << "x" << nc
     << " qr=" << rank_backend << " refined=" << rank_refined << " dense=" << rk
     << std::scientific << std::setprecision(3);
  // rev142a: "none" rather than 0 when the index does not exist.  Printing
  // sigma[r]=0.000e+00 on a FULL-RANK square system reads as "there is an exactly zero
  // singular value" -- the precise opposite of what the line is reporting.  Every
  // full-rank line in the first PDE16 oracle log said that.
  os << " sigma[r-1]=";  if( rk )            os << sv(rk-1); else os << "none";
  os << " sigma[r]=";    if( rk < sv.n_elem ) os << sv(rk);   else os << "none";
  os << " tolabs=" << tol_abs << std::endl;
  os.flags( fl ); os.precision( pr );
}

inline OCFEBASE::t_RankResult
OCFEBASE::numeric_rank
( arma::sp_mat const& M, double reltol, double tol_abs_override )
const
{
  t_RankResult R;
  R.tol = reltol;
  size_t const nr = M.n_rows, nc = M.n_cols;
  if( !nr || !nc ){ R.rank = 0; R.backend = t_RankResult::DENSE_SVD; R.gap = std::numeric_limits<double>::infinity(); return R; }
  auto const t0 = std::chrono::steady_clock::now();

  // THE tolerance.  All backends compare against ONE absolute threshold:
  //   tol_abs = max( reltol, 20 (nr+nc) eps ) * max_j ||M(:,j)||_2
  // The second term is SPQR's own default (Davis, Algorithm 915): the round-off
  // floor an EXACTLY dependent column reaches after nr+nc Householder updates.
  // MEASURED, PDE6_blk0 IC_STRONG MANUAL, n = 1250: with tol = 1e-12 * maxcol the
  // dense SVD finds 4 null directions with a 5.4e12 gap, while SPQR declared 3
  // dead columns and RETAINED the fourth at |R_jj| = 8e-12 * maxcol -- an exact
  // dependency read through accumulated round-off, not a near-singular
  // direction.  The size floor (1.1e-11 * maxcol there) catches it; PDE40b's
  // genuine near-singular direction at 4.6e-10 relative (n = 1008, floor 9e-12)
  // stays retained, which is the other half of the requirement.
  // The audit computes tol_abs ONCE from the full Jacobian and passes it to the
  // rank([B;D]) call as the override, so both ranks see the same threshold.
  M.sync();
  double maxcol = 0.;
  for( arma::uword j = 0; j < M.n_cols; ++j ){
    double s2 = 0.;
    for( arma::sp_mat::const_col_iterator it = M.begin_col(j); it != M.end_col(j); ++it ) s2 += (*it)*(*it);
    maxcol = std::max( maxcol, std::sqrt( s2 ) );
  }
  double const eps  = std::numeric_limits<double>::epsilon();
  double const rel  = std::max( reltol, 20. * (double)( nr + nc ) * eps );
  double const tol_abs = ( tol_abs_override > 0. ) ? tol_abs_override : rel * ( maxcol > 0. ? maxcol : 1. );
  R.tol_abs = tol_abs;

  auto want = [&]( Options::DETERMINACY_BACKEND_T b ){
    return options.DETERMINACY.BACKEND == Options::DET_AUTO || options.DETERMINACY.BACKEND == b; };

#if defined(CRONOS__WITH_SPQR)
  if( want( Options::DET_SPQR ) ){
    // arma CSC -> cholmod (long) sparse, no copy of the values.
    std::vector<int64_t> cp( nc + 1 ), ri( M.n_nonzero ? M.n_nonzero : 1 );
    for( size_t j = 0; j <= nc; ++j ) cp[j] = (int64_t)M.col_ptrs[j];
    for( size_t k = 0; k < M.n_nonzero; ++k ) ri[k] = (int64_t)M.row_indices[k];
    cholmod_common cc;  cholmod_l_start( &cc );
    // rev142b: CHOLMOD's own OpenMP regions honour this; guarded on the CHOLMOD major
    // version because nthreads_max does not exist before CHOLMOD 4 (SuiteSparse 6).
#if defined(CHOLMOD_MAIN_VERSION) && (CHOLMOD_MAIN_VERSION >= 4)
    if( int const _cap = t_ThreadCap::active() ) cc.nthreads_max = _cap;
#endif
    cholmod_sparse Ac;  std::memset( &Ac, 0, sizeof(Ac) );
    Ac.nrow = nr; Ac.ncol = nc; Ac.nzmax = M.n_nonzero;
    Ac.p = cp.data(); Ac.i = ri.data(); Ac.x = const_cast<double*>( M.values );
    Ac.stype = 0; Ac.itype = CHOLMOD_LONG; Ac.xtype = CHOLMOD_REAL; Ac.dtype = CHOLMOD_DOUBLE;
    Ac.sorted = 1; Ac.packed = 1;
    cholmod_sparse* Rc = nullptr; int64_t* E = nullptr;
    // rev142: ORDERING_BEST default (corpus-validated verdict-invariant; see the same
    // note in _null_basis).  CRONOS_SPQR_ORDERING overrides.
    static int const kSpqrOrd = SPQR_ORDERING_BEST;   // CRONOS_SPQR_ORDERING (retired 2026-10-07, WORKPLAN 3.B batch 2b)
    int64_t const rk = SuiteSparseQR<double>( kSpqrOrd, tol_abs,
                                                       (int64_t)nr, &Ac, &Rc, &E, &cc );
    R.backend = t_RankResult::SPQR;
    R.rank = ( rk >= 0 ? (size_t)rk : 0 );
    size_t const rk_qr_raw = R.rank;    // rev142: as SPQR reported it, for the oracle
    // rev142c: SPQR's own statistics for the factorization just performed (default off).
    { t_SpqrStat S;
      S.nrow = nr; S.ncol = nc; S.nnzA = M.n_nonzero; S.rank = R.rank;
      S.nnzR = cc.SPQR_istat[0]; S.ncolsing = cc.SPQR_istat[5]; S.nrowsing = cc.SPQR_istat[6];
      S.ordering = (int)cc.SPQR_istat[7]; S.flops = cc.SPQR_flopcount;
      S.t_analyze = cc.SPQR_analyze_time; S.t_factorize = cc.SPQR_factorize_time;
      _spqr_stat_report( S, "rank" ); }

    // Leading r x r block of R (the retained, permuted columns) -> condition
    // estimate + deflation (see _refine_rank_from_R).
    R.gap = std::numeric_limits<double>::infinity(); R.sigma_cut = 0.;
    if( Rc && Rc->x && Rc->p && Rc->i && R.rank ){
      int64_t const* rp = static_cast<int64_t const*>( Rc->p );
      int64_t const* rI = static_cast<int64_t const*>( Rc->i );
      double const* rx = static_cast<double const*>( Rc->x );
      size_t const r = R.rank;
      // SPQR returns a SQUEEZED R: rank-by-n, dead columns left IN PLACE as zero
      // columns (not permuted to the end).  The live columns form a staircase: a
      // column is live iff its lowest-indexed entry is the next unfilled row.
      // R11 = R restricted to the live columns is r x r upper triangular.
      // rev142: the walk moved into _staircase_live_columns (O(1) per column on sorted
      // R, rev141 scan retained as the fallback); it was 39% of an audit rank pass.
      std::vector<size_t> live;
      _staircase_live_columns( rp, rI, rx, (size_t)Rc->ncol, r, Rc->sorted != 0, live );
      if( live.size() == r ){
        std::vector<char> is_live( Rc->ncol, 0 ); for( size_t j : live ) is_live[j] = 1;
        std::vector<size_t> dead; for( size_t j = 0; j < Rc->ncol; ++j ) if( !is_live[j] ) dead.push_back( j );
        auto build = [&]( std::vector<size_t> const& cols, arma::sp_mat& out ){
          size_t cnt = 0;
          for( size_t j : cols ) for( int64_t k = rp[j]; k < rp[j+1]; ++k ) if( (size_t)rI[k] < r ) ++cnt;
          arma::umat loc( 2, (arma::uword)std::max<size_t>( cnt, 1 ) ); arma::vec val( (arma::uword)std::max<size_t>( cnt, 1 ) );
          arma::uword e2 = 0;
          for( size_t jj = 0; jj < cols.size(); ++jj ){ size_t const j = cols[jj];
            for( int64_t k = rp[j]; k < rp[j+1]; ++k ) if( (size_t)rI[k] < r ){
              loc(0,e2) = (arma::uword)rI[k]; loc(1,e2) = (arma::uword)jj; val(e2) = rx[k]; ++e2; } }
          out = arma::sp_mat( (arma::uword)r, (arma::uword)cols.size() );
          if( e2 ) out = arma::sp_mat( loc.cols(0,e2-1), val.subvec(0,e2-1), (arma::uword)r, (arma::uword)cols.size() );
        };
        arma::sp_mat R11, R12;
        build( live, R11 ); build( dead, R12 );
        size_t const rk_qr = R.rank;
        _refine_rank_from_R( R11, R12, tol_abs, R.rank, R.gap, R.sigma_cut );
        static bool const kDbg = Options::_env_flag( "CRONOS_DEBUG_RANK", false );
        if( kDbg ){
          _disp( 3 ) << "  [rank-dbg] SPQR: n=" << nr << "x" << nc << " rk_qr=" << rk_qr << " live=" << live.size()
                    << " dead=" << dead.size() << " R12.nnz=" << R12.n_nonzero
                    << std::scientific << std::setprecision(3) << " tol_abs=" << tol_abs
                    << " -> rank=" << R.rank << " sigma_cut=" << R.sigma_cut << " gap=" << R.gap;
          if( r <= 2000 ){
            arma::vec sv; arma::mat Rd( arma::join_rows( arma::mat( R11 ), arma::mat( R12 ) ) );
            if( arma::svd( sv, Rd ) ){
              _disp( 3 ) << "  dense-svd([R11 R12]): s_small=";
              for( arma::uword q = 0; q < 6 && q < sv.n_elem; ++q ) _disp( 3 ) << " " << sv( sv.n_elem - 1 - q );
            }
          }
          _disp( 3 ) << std::defaultfloat << std::endl;
        }
      }
    }
    // rev142: dense cross-check, default off.  OUTSIDE the staircase block on purpose --
    // the paths that skip the refinement (rank 0, a staircase that did not close) are
    // exactly the ones an oracle should be able to see.
    _rank_oracle_check( M, tol_abs, rk_qr_raw, R.rank, "SPQR" );
    if( Rc ) cholmod_l_free_sparse( &Rc, &cc );
    if( E )  cholmod_l_free( nc, sizeof(int64_t), E, &cc );
    cholmod_l_finish( &cc );
    R.ms = std::chrono::duration<double,std::milli>( std::chrono::steady_clock::now() - t0 ).count();
    return R;
  }
#endif

#if defined(CRONOS__WITH_EIGEN)
  if( want( Options::DET_EIGEN ) ){
    typedef Eigen::SparseMatrix<double,Eigen::ColMajor,long> ESp;
    ESp Es( (long)nr, (long)nc );
    std::vector<Eigen::Triplet<double,long>> trip; trip.reserve( M.n_nonzero );
    for( arma::sp_mat::const_iterator it = M.begin(); it != M.end(); ++it )
      trip.emplace_back( (long)it.row(), (long)it.col(), *it );
    Es.setFromTriplets( trip.begin(), trip.end() );
    Es.makeCompressed();
    Eigen::SparseQR<ESp, Eigen::COLAMDOrdering<long>> qr;
    qr.setPivotThreshold( tol_abs );
    qr.compute( Es );
    R.backend = t_RankResult::EIGEN_QR;
    size_t rk_qr_raw = 0;               // rev142: as Eigen reported it, for the oracle
    if( qr.info() != Eigen::Success ){ R.rank = 0; R.gap = 0.; R.sigma_cut = 0.; }
    else{
      R.rank = (size_t)qr.rank();
      rk_qr_raw = R.rank;
      R.gap = std::numeric_limits<double>::infinity(); R.sigma_cut = 0.;
      if( R.rank ){
        size_t const r = R.rank;
        ESp const Rm = qr.matrixR();
        std::vector<Eigen::Triplet<double,long>> t11, t12;
        for( long j = 0; j < Rm.outerSize(); ++j )
          for( ESp::InnerIterator it( Rm, j ); it; ++it ){
            if( (size_t)it.row() >= r ) continue;
            if( (size_t)j < r ) t11.emplace_back( it.row(), j, it.value() );
            else                t12.emplace_back( it.row(), j - (long)r, it.value() );
          }
        auto build = [&]( std::vector<Eigen::Triplet<double,long>> const& tr, size_t ncols, arma::sp_mat& out ){
          arma::umat loc( 2, (arma::uword)std::max<size_t>( tr.size(), 1 ) ); arma::vec val( (arma::uword)std::max<size_t>( tr.size(), 1 ) );
          for( size_t k = 0; k < tr.size(); ++k ){ loc(0,k) = (arma::uword)tr[k].row(); loc(1,k) = (arma::uword)tr[k].col(); val(k) = tr[k].value(); }
          out = arma::sp_mat( (arma::uword)r, (arma::uword)ncols );
          if( !tr.empty() ) out = arma::sp_mat( loc.cols(0,tr.size()-1), val.subvec(0,tr.size()-1), (arma::uword)r, (arma::uword)ncols );
        };
        arma::sp_mat R11, R12;
        build( t11, r, R11 ); build( t12, ( (size_t)Rm.cols() > r ? (size_t)Rm.cols() - r : 0 ), R12 );
        _refine_rank_from_R( R11, R12, tol_abs, R.rank, R.gap, R.sigma_cut );
      }
    }
    _rank_oracle_check( M, tol_abs, rk_qr_raw, R.rank, "EIGEN" );   // rev142, default off
    R.ms = std::chrono::duration<double,std::milli>( std::chrono::steady_clock::now() - t0 ).count();
    return R;
  }
#endif

  // Dense SVD fallback under the size cap: exact singular values, exact gap.
  if( want( Options::DET_DENSE ) && std::max( nr, nc ) <= _knob_DETERMINACY_MAX_DENSE() ){
    arma::mat Dm( M );
    arma::vec s;
    R.backend = t_RankResult::DENSE_SVD;
    if( !arma::svd( s, Dm ) || s.n_elem == 0 ){ R.rank = 0; R.gap = 0.; }
    else{
      size_t rk = 0;
      for( arma::uword k = 0; k < s.n_elem; ++k ) if( s(k) > tol_abs ) ++rk;
      R.rank = rk;
      R.sigma_cut = ( rk ? s(rk-1) : 0. );
      // gap: retained/first-cut when something was cut, else retained/tolerance.
      if( rk == 0 )               R.gap = 0.;
      else if( rk < s.n_elem )    R.gap = ( s(rk) > 0. ? s(rk-1) / s(rk) : std::numeric_limits<double>::infinity() );
      else                        R.gap = ( tol_abs > 0. ? s(rk-1) / tol_abs : std::numeric_limits<double>::infinity() );
    }
    R.ms = std::chrono::duration<double,std::milli>( std::chrono::steady_clock::now() - t0 ).count();
    return R;
  }

  R.backend = t_RankResult::NONE;
  R.ms = std::chrono::duration<double,std::milli>( std::chrono::steady_clock::now() - t0 ).count();
  return R;
}

// rev140: the audit's per-window print, verbatim (byte-identical output), so the
// restore-attempt recursion can audit quietly and the kept system's lines can be
// replayed once after the keep/revert decision.
inline void
OCFEBASE::_print_determinacy_window
( t_Determinacy const& D, size_t w ) const
{
  auto A_n_ne_m = []( t_Determinacy const& D ){ return D.n != D.m; };   // square: defL == defR trivially
  auto verdict_name = []( t_Determinacy::Verdict v ) -> char const* {
    switch( v ){ case t_Determinacy::DETERMINED: return "DETERMINED"; case t_Determinacy::GAUGE_ONLY: return "DETERMINED+IMPLIED";
                 case t_Determinacy::STATE_CONTENT: return "** STATE CONTENT **";
                 case t_Determinacy::INCONSISTENT_TOL: return "INCONSISTENT_TOL"; default: return "SKIPPED"; } };
    std::ios::fmtflags const cf = std::cerr.flags();
    std::streamsize const cp = std::cerr.precision();
    _disp( 2 ) << "[determinacy] win=" << w << " mode=" << D.mode;
    if( !D.ran ){
      _disp( 2 ) << "  SKIPPED reason=" << D.skipped_reason << " ref=" << D.ref_point << std::endl;
      std::cerr.flags( cf ); std::cerr.precision( cp );
      return;
    }
    _disp( 2 ) << " n=" << D.n << " m=" << D.m << " rank=" << D.rank
              << " defR=" << D.def_right << " defL=" << D.def_left << " k'=" << D.k_prime;
    if( D.free_expected ) _disp( 2 ) << " free=" << D.free_expected;
    if( D.dep_expected )  _disp( 2 ) << " dep=" << D.dep_expected;
    _disp( 2 ) << "  " << verdict_name( D.verdict );
    // rev170: report the gauge as what it IS.  MEASURED corpus-wide (rev168 sweep): the
    // IC_STRONG k' EQUALS the interface plan's tau-block deficit on every driver that carries
    // one -- MBC1 12, PDE2 12, PDE3 12, PDE30/33/35 10, PDE30b 10-35, PDE38 21-105, PDE39 30,
    // scalar_nest 1, hybrid_adequacy 1 -- and equals 0 exactly where the claim drop removed the
    // claims (PDE20f, PDE31b, STRONGCORNER1).  So k' counts CONTINUITY CLAIMS THAT THE PLAN
    // IMPLIES: an algebraic state whose defining relation is pointwise in states inherits its
    // continuity from its arguments, and the claim on it is redundant.  Realised as tau it
    // mints a multiplier that enforces nothing -- which is exactly k'.
    // The STATE block is fully determined in this verdict; the deficiency is entirely in the
    // multiplier block.  "GAUGE_ONLY" read as a failure; it is not one, and 15 corpus drivers
    // carry it while passing.  The count stays visible and countable because a CHANGE in it is
    // a real signal.
    if( D.verdict == t_Determinacy::GAUGE_ONLY ){
      _disp( 2 ) << "(" << D.k_prime << ")";
      static bool once = false;
      if( !once ){
        once = true;
        _disp( 2 ) << "  -- " << D.k_prime << " continuity claim(s) implied by the plan, each"
                     " realised as a tau multiplier that enforces nothing; the STATE block is"
                     " fully determined";
      }
    }
    if( D.free_expected && D.verdict == t_Determinacy::DETERMINED ) _disp( 2 ) << " (FULL_ROW_RANK)";
    if( D.dep_excess && A_n_ne_m( D ) ) _disp( 2 ) << " [+" << D.dep_excess << " dependent row(s)]";
    _disp( 2 ) << std::scientific << std::setprecision(1)
              << "  gap=" << D.gap << " tol=" << options.DETERMINACY.TOL << " tolabs=" << D.tol_abs
              << " be=" << D.backend << " ref=" << D.ref_point
              << std::fixed << std::setprecision(1)
              << " ms=" << ( D.ms_assemble + D.ms_rank ) << "(A=" << D.ms_assemble << ",R=" << D.ms_rank << ")"
              << std::defaultfloat;
    if( D.marginal_reeval )
      _disp( 2 ) << "  [" << ( D.verdict != t_Determinacy::DETERMINED ? "marginal" : "near-singular" )
                << ": re-eval at 1.3x -> rank=" << D.rank2
                << " gap=" << std::scientific << std::setprecision(1) << D.gap2 << std::defaultfloat
                << " " << verdict_name( D.verdict2 ) << "]";
    _disp( 2 ) << std::endl;
    std::cerr.flags( cf ); std::cerr.precision( cp );
}

inline bool
OCFEBASE::_determinacy_verdict
( t_JacAssembly const& A, t_Determinacy& D )
const
{
auto be_name = []( t_RankResult::Backend b ) -> char const* {
  switch( b ){ case t_RankResult::SPQR: return "SPQR"; case t_RankResult::EIGEN_QR: return "EIGEN";
               case t_RankResult::DENSE_SVD: return "DENSE"; default: return "NONE"; } };
  t_RankResult const RM = numeric_rank( A.J, options.DETERMINACY.TOL );
  D.tol_abs = RM.tol_abs;
  if( RM.backend == t_RankResult::NONE ){
    D.skipped_reason = "no rank backend for n=" + std::to_string( std::max( A.n, A.m ) )
                     + " (dense cap " + std::to_string( _knob_DETERMINACY_MAX_DENSE() ) + ", no SPQR/Eigen)";
    return false;
  }
  D.rank = RM.rank; D.gap = RM.gap; D.sigma_cut = RM.sigma_cut; D.backend = be_name( RM.backend ); D.ms_rank = RM.ms;
  D.def_right = ( A.m > RM.rank ? A.m - RM.rank : 0 );
  D.def_left  = ( A.n > RM.rank ? A.n - RM.rank : 0 );
  D.k_prime = 0;
  if( A.n_tau && A.tau0 + A.n_tau <= A.m ){
    arma::sp_mat const T = A.J.cols( (arma::uword)A.tau0, (arma::uword)( A.tau0 + A.n_tau - 1 ) );
    t_RankResult const RT = numeric_rank( T, options.DETERMINACY.TOL, RM.tol_abs );   // SAME threshold
    D.ms_rank += RT.ms;
    D.k_prime = ( A.n_tau > RT.rank ? A.n_tau - RT.rank : 0 );
    _keep_explicit_audit( T, A.tau0, RT.rank );   // rev142d, default off
    // rev150b [taunull]: name the COLLECTIVE tau-block null directions.
    // MEASURED (MBC1/MBC3 IC_TRACE): every column individually coupled
    // (taucol lost=0) yet k'>0 -- B is collectively rank-deficient, and the
    // step-fail census shows equal-and-opposite column pairs.  This census
    // prints each null direction's top support so the claim pairs can be
    // identified.  CRONOS_TAUNULL_CENSUS=1; default off.
    { static bool const kCensus =
        Options::_env_flag( "CRONOS_TAUNULL_CENSUS", false );
      if( kCensus && D.k_prime ){
        size_t rkT = 0; arma::mat NB;
        if( _null_basis( T, RM.tol_abs, rkT, NB ) && NB.n_cols ){
          _disp( 3 ) << "  [taunull] tau-block " << A.m << "x" << A.n_tau
                    << " rank=" << RT.rank << " null_dim=" << NB.n_cols
                    << " tau0=" << A.tau0 << std::endl;
          for( size_t d = 0; d < NB.n_cols && d < 16; ++d ){
            arma::vec const av = arma::abs( NB.col( d ) );
            arma::uvec const idx = arma::sort_index( av, "descend" );
            _disp( 3 ) << "  [taunull] dir#" << d << ":";
            for( size_t j = 0; j < idx.n_elem && j < 6; ++j ){
              if( av( idx(j) ) < 1e-8 ) break;
              _disp( 3 ) << "  tau" << idx(j)
                        << "(col" << A.tau0 + idx(j) << ")="
                        << NB( idx(j), d );
            }
            _disp( 3 ) << std::endl;
          }
          if( NB.n_cols > 16 )
            _disp( 3 ) << "  [taunull] ... (" << NB.n_cols
                      << " directions total)" << std::endl;
        }
      }
    }
  }
  // rev150h: dump the nonzeros of selected columns of the audited (equilibrated) J.
  // CRONOS_DUMP_COLS=c1,c2,...  (var[] indices).  Prints row, row_id, value.
  if( char const* dc = std::getenv( "CRONOS_DUMP_COLS" ) ){
    std::string const spec( dc ); size_t pos = 0;
    bool const have_prov = ( _rowProvenanceReady && _rowEqnNdx.size() == (size_t)A.J.n_rows );
    while( pos < spec.size() ){
      size_t nxt = spec.find( ',', pos ); if( nxt == std::string::npos ) nxt = spec.size();
      size_t const c = (size_t)std::atol( spec.substr( pos, nxt - pos ).c_str() ); pos = nxt + 1;
      if( c >= (size_t)A.J.n_cols ) continue;
      _disp( 3 ) << "  [dumpcol] c" << c << ":";
      for( arma::sp_mat::const_col_iterator it = A.J.begin_col( (arma::uword)c ); it != A.J.end_col( (arma::uword)c ); ++it )
        _disp( 3 ) << "  r=" << it.row() << "(row_id=" << ( have_prov ? _rowEqnNdx[ it.row() ] : -2 ) << ")=" << std::scientific << std::setprecision(3) << (*it);
      _disp( 3 ) << std::endl;
    }
  }
  // rev150i/153c: dump the audited J (equilibrated) in coordinate form: CRONOS_DUMP_J=<file>
  if( char const* dj = std::getenv( "CRONOS_DUMP_J" ) ){
    std::ofstream fo( dj );
    bool const have_prov = ( _rowProvenanceReady && _rowEqnNdx.size() == (size_t)A.J.n_rows );
    fo << "# n_rows " << A.J.n_rows << " n_cols " << A.J.n_cols << " tau0 " << ( A.n_tau ? A.tau0 : (size_t)A.J.n_cols ) << "\n";
    for( arma::sp_mat::const_iterator it = A.J.begin(); it != A.J.end(); ++it )
      fo << it.row() << " " << it.col() << " " << std::setprecision(17) << (*it) << " " << ( have_prov ? _rowEqnNdx[ it.row() ] : -2 ) << "\n";
    _disp( 3 ) << "  [dumpJ] wrote " << A.J.n_nonzero << " nonzeros to " << dj << std::endl;
  }
  // rev150c [statenull]: name the STATE-SIDE null directions (def_right > k').
  // Runs under every imposition; IC_WEAK has n_tau=0 so every direction is
  // state-side.  CRONOS_STATENULL_CENSUS=1; default off; prints only.
  { static bool const kSCensus = Options::_env_flag( "CRONOS_STATENULL_CENSUS", false );
    if( kSCensus && D.def_right > D.k_prime ) _statenull_census( A, RM.tol_abs, D.k_prime );
  }
  // Rectangular systems (rev134): m > n leaves m-n free columns BY DESIGN (decision
  // variables: the q field of PDE6, refinit's reduced-space inputs -- MEASURED
  // rev133 sweep: PDE6 1875x2500 rank 1875, refinit 72x80 rank 72, both full row
  // rank and both reported STATE CONTENT by the square-system rule).  The verdict
  // is on the free directions BEYOND m-n; n > m expects n-m dependent rows.
  D.free_expected = ( A.m > A.n ? A.m - A.n : 0 );
  D.dep_expected  = ( A.n > A.m ? A.n - A.m : 0 );
  size_t const excess = ( D.def_right > D.free_expected ? D.def_right - D.free_expected : 0 );
  D.dep_excess = ( D.def_left > D.dep_expected ? D.def_left - D.dep_expected : 0 );
  if( D.def_right < D.free_expected )  D.verdict = t_Determinacy::INCONSISTENT_TOL;   // impossible: rank > n
  else if( excess == 0 )               D.verdict = t_Determinacy::DETERMINED;
  else if( excess == D.k_prime )       D.verdict = t_Determinacy::GAUGE_ONLY;
  else if( excess >  D.k_prime )       D.verdict = t_Determinacy::STATE_CONTENT;
  else                                 D.verdict = t_Determinacy::INCONSISTENT_TOL;
  return true;
}

// ======================================================================
// OCFESLV::_compute_interface_reads   (item 11 Part-2: per-direction symbol read)
// ======================================================================
// Solve-free, always compiled.  Populates _blockInterfaceRead with the
// per-(primitive,direction) order + singular-perturbation magnitude balance
// consumed by the symbol-AUTO resolver.  Under MC__OCFESLV_SYMBOL_BALANCE_PROBE
// it additionally prints the measure-before-implement PRIMARY/CORRECTOR/now
// table (the gated diagnostic).  Order>2 is flagged, not handled.
inline void
OCFEBASE::_compute_interface_reads
( std::vector<double> const& state_ref,
  std::vector<double> const& cst_ref,
  std::vector<double> const& dom_ref,
  std::vector<double> const& input_ref )
{
  using IT = Options::InterfaceType;
  // rev67 CORRECTION-OF-A-CORRECTION.  An earlier rev67 draft widened kPrint to
  // DISPLAY_LEVEL >= 1, intending to surface the CLASSIFICATION VERDICT.  That was the
  // wrong gate: kPrint controls the whole _symbol_balance_probe developer output
  // ("item-11 Part-2", the per-block corrector read, the corrector summary), NOT the
  // classify verdict, which is emitted elsewhere.  The change therefore turned on a
  // developer probe and left the actual problem unfixed.  Reverted to the macro.
  //
  // The verdict is surfaced instead where it belongs -- see the classification stamp
  // below the [evolution] report, which prints once per block unconditionally, on the
  // same footing as the HEADER_ID stamp: it is PROVENANCE, not diagnostics.
#ifdef MC__OCFESLV_SYMBOL_BALANCE_PROBE
  constexpr bool kPrint = true;
#else
  constexpr bool kPrint = false;
#endif
  double const PI = 3.14159265358979323846;
  auto it_name = []( IT t ) -> const char* {
    switch( t ){ case Options::IC_VALUE:  return "VALUE";
                 case Options::IC_FLUX:   return "FLUX";
                 case Options::IC_UPWIND: return "UPWIND";
                 default:                 return "AUTO"; } };
  (void)it_name; (void)PI;

  _blockInterfaceRead.clear();

  if( kPrint )
    _disp( 3 ) << "OCFESLV::_symbol_balance_probe ** item-11 Part-2 (order<=2): "
                 "primary = principal-symbol characteristic read; "
                 "corrector = scoped singular-perturbation guard (k_h vs k*=adv/diff)\n"
                 "  (an order-2 read prints IC_FLUX, i.e. flux continuity; under order reduction the resolver realises it "
                 "as IC_VALUE on the auxiliary, and IC_AUTO never yields IC_FLUX -- so 'now' shows VALUE there)\n";

  size_t fired_total = 0, ord2_total = 0;

  for( auto const& [bid, sym] : _blockSymbol ){
    auto itc = _blockClassification.find( bid );
    if( itc == _blockClassification.end() ) continue;
    t_Classify const& cls = itc->second;

    size_t const nDom   = sym.vDom.size();
    size_t const nState = sym.vState.size();
    size_t const nEqn   = sym.vEqn.size();
    if( !nDom || !nState || !nEqn ) continue;

    try {
    // Principal first-derivative coefficients at the reference point.
    std::vector<arma::mat> Ai = cls.Ai;
    if( Ai.size() != nDom )
      Ai = _eval_symbol( sym, state_ref, cst_ref, dom_ref, input_ref );

    // rev67: _eval_symbol returns EMPTY when the block has no principal coefficients at
    // all -- a distributed purely algebraic field, where vDom is populated but vCoeff is
    // not.  This routine only READS the classification (cls is a const reference); the
    // VERDICT is set in classify(), which distinguishes nDom==0 -> ALGEBRAIC_LUMPED from
    // Ai.empty() with domains -> ALGEBRAIC_FIELD.  Here there is simply no principal
    // symbol to derive interface reads from, so the block is skipped.
    if( Ai.empty() ) continue;

    // Zeroth-order coupling A0 (requires vCoeff0 from _principal_symbol).
    if( sym.vCoeff0.size() != nEqn * nState ){
      if( kPrint )
        _disp( 3 ) << "  blk" << bid << ": vCoeff0 unavailable; skipped\n";
      continue;
    }
    arma::mat A0( (arma::uword)nEqn, (arma::uword)nState, arma::fill::zeros );
    {
      auto by_id = []( FFVar const& v, std::vector<FFVar> const& vars,
                       std::vector<double> const& vals, size_t fb ) -> double {
        for( size_t k = 0; k < vars.size() && k < vals.size(); ++k )
          if( vars[k].id() == v.id() ) return vals[k];
        return fb < vals.size() ? vals[fb] : 0.; };
      std::vector<FFVar>  ev;  std::vector<double> evv;
      for( size_t j = 0; j < nState; ++j ){
        ev.push_back( sym.vState[j] );
        evv.push_back( by_id( sym.vState[j], _vVar, state_ref, j ) ); }
      // provision non-principal coefficient states (e.g. wetting omega) absent from sym.vState
      for( size_t p = 0; p < _vVar.size(); ++p ){
        bool principal = false;
        for( size_t j = 0; j < nState; ++j )
          if( sym.vState[j].id() == _vVar[p].id() ){ principal = true; break; }
        if( principal ) continue;
        ev.push_back( _vVar[p] );
        evv.push_back( by_id( _vVar[p], _vVar, state_ref, p ) ); }
      for( size_t i = 0; i < nDom; ++i ){
        ev.push_back( sym.vDom[i] );
        evv.push_back( by_id( sym.vDom[i], _vDom, dom_ref, i ) ); }
      for( size_t p = 0; p < _vInp.size(); ++p ){
        ev.push_back( _vInp[p] );
        evv.push_back( by_id( _vInp[p], _vInp, input_ref, p ) ); }
      for( size_t c = 0; c < _vCst.size(); ++c ){
        ev.push_back( _vCst[c] );
        evv.push_back( c < cst_ref.size() ? cst_ref[c] : 0. ); }
      // vCoeff0 may reference _principal_symbol's derivative proxies when a
      // derivative coefficient is state-dependent (e.g. gas-expanding velocity).
      // Their value does not affect the read (only aux columns of A0 are
      // consumed, and an aux appears solely inside Partial[...], so its algebraic
      // Jacobian entry is proxy-free); supply 0 so the eval does not throw.
      for( size_t q = 0; q < sym.vDerivProxy.size(); ++q ){
        ev.push_back( sym.vDerivProxy[q] );
        evv.push_back( 0. ); }
      std::vector<double> c0( nEqn * nState, 0. );
      _dag->eval( sym.vCoeff0, c0, ev, evv );
      for( size_t k = 0; k < nEqn; ++k )
        for( size_t j = 0; j < nState; ++j )
          A0( (arma::uword)k, (arma::uword)j ) = c0[ k * nState + j ];
    }

    // Per-row role (INTERIOR vs LINK) by matching sym.vEqn against _mEqn.
    std::vector<EqnRole> role( nEqn, EqnRole::INTERIOR );
    for( size_t k = 0; k < nEqn; ++k )
      for( auto const& e : _mEqn )
        if( e.var.id() == sym.vEqn[k].id() ){ role[k] = e.opt->role; break; }

    // Structural (primitive,direction) -> aux map from LINK rows:
    //   LINK row  d(prim)/dx_i - aux = 0  ->  Ai[i][k,prim]~+-1, A0[k,aux]~+-1.
    std::map<std::pair<size_t,size_t>,size_t> aux_of;
    std::set<size_t> aux_cols;
    for( size_t k = 0; k < nEqn; ++k ){
      if( role[k] != EqnRole::LINK ) continue;
      size_t bi = nDom, bp = nState; double bv = 0.;
      for( size_t i = 0; i < nDom; ++i )
        for( size_t j = 0; j < nState; ++j ){
          double a = std::fabs( Ai[i]( (arma::uword)k, (arma::uword)j ) );
          if( a > bv ){ bv = a; bi = i; bp = j; } }
      size_t ba = nState; double bav = 0.;
      for( size_t j = 0; j < nState; ++j ){
        double a = std::fabs( A0( (arma::uword)k, (arma::uword)j ) );
        if( a > bav ){ bav = a; ba = j; } }
      if( bi < nDom && bp < nState && ba < nState && bv > 1e-12 && bav > 1e-12 ){
        aux_of[ { bp, bi } ] = ba;
        aux_cols.insert( ba );
      }
    }

    // Mesh wavenumber k_h = pi / dx_min per direction.  Read the node cache
    // directly: node_colloc() throws (Exceptions::SETUP) when called during
    // setup, but _setup_node_cache() has already populated _nodeBlockCache by
    // the classify phase.  Missing entry -> k_h = NaN (corrector reports no k_h).
    std::vector<double> kh( nDom, std::numeric_limits<double>::quiet_NaN() );
    for( size_t i = 0; i < nDom; ++i ){
      auto itnc = _nodeBlockCache.find( sym.vDom[i] );
      if( itnc == _nodeBlockCache.end() ) continue;
      double dxmin = std::numeric_limits<double>::infinity();
      for( auto const& el : itnc->second )
        for( size_t t = 1; t < el.size(); ++t ){
          double d = std::fabs( el[t] - el[t-1] );
          if( d > 0. && d < dxmin ) dxmin = d; }
      if( std::isfinite( dxmin ) && dxmin > 0. ) kh[i] = PI / dxmin;
    }

    // Block-level hyperbolicity proxy for the primary order-1 verdict.
    bool const hyp = cls.evolution_hyperbolic && !cls.weak_hyperbolic;

    std::vector<t_DirRead>& reads = _blockInterfaceRead[bid];

    if( kPrint ){
      _disp( 3 ) << "  blk" << bid << " type=" << pde_type_name( cls.type )
                << "  (primary = char read; corrector fires only on two-scale structure)\n";
      _disp( 3 ) << "    " << std::left << std::setw(16) << "state"
                << std::setw(4) << "dir" << std::right
                << " ord  " << std::left << std::setw(7) << "primary" << std::right
                << "  eps_eff    k*         k_h        act "
                << std::left << std::setw(9) << "corrector" << std::right
                << " " << std::left << std::setw(7) << "final" << std::right
                << " now\n";
    }

    double const tol = 1e-30;
    for( size_t jp = 0; jp < nState; ++jp ){
      if( aux_cols.count( jp ) ) continue;            // report primitives only
      for( size_t i = 0; i < nDom; ++i ){
        auto it = aux_of.find( { jp, i } );
        int ord = 0; double adv = 0., diff = 0.;
        double kstar = std::numeric_limits<double>::quiet_NaN();
        double eps   = std::numeric_limits<double>::quiet_NaN();
        bool order_gt2 = false;
        bool via_aux   = false;

        if( it == aux_of.end() ){
          double d1 = 0.;
          for( size_t k = 0; k < nEqn; ++k )
            if( role[k] != EqnRole::LINK )
              d1 = std::max( d1, std::fabs( Ai[i]( (arma::uword)k, (arma::uword)jp ) ) );
          ord = ( d1 > tol ) ? 1 : 0;
        } else {
          via_aux = true;
          size_t ja = it->second;
          if( aux_of.find( { ja, i } ) != aux_of.end() ) order_gt2 = true;
          // adv = first-order transport of the primitive in dir i, captured in
          // EITHER representation: the aux used algebraically (A0[.,aux], e.g.
          // PDE2 vl*Dz_Cl) OR the primitive's own first derivative
          // (Ai[i][.,prim], e.g. PDE0 u*dC/dz).  reduce_order may route the
          // advection through the aux or leave it on the primitive, so both
          // slots must be checked or advection-dominated cases are missed.
          for( size_t k = 0; k < nEqn; ++k ){
            if( role[k] == EqnRole::LINK ) continue;
            diff = std::max( diff, std::fabs( Ai[i]( (arma::uword)k, (arma::uword)ja ) ) );
            adv  = std::max( adv, std::max(
                     std::fabs( A0    ( (arma::uword)k, (arma::uword)ja ) ),
                     std::fabs( Ai[i] ( (arma::uword)k, (arma::uword)jp ) ) ) );
          }
          ord = ( diff > tol ) ? 2 : 1;
          if( diff > tol ){ kstar = adv / diff; eps = diff / adv; }
        }

        if( ord == 0 ) continue;                       // algebraic in i; no face coupling

        // Singular-perturbation corrector: top (diffusive) term sub-grid when
        // k_h < k* = adv/diff.  Counted for non-order_gt2 order-2 only.
        bool adv_dom = false;
        if( ord >= 2 && !order_gt2 ){
          ++ord2_total;
          if( std::isfinite( kh[i] ) && std::isfinite( kstar ) ){
            adv_dom = ( kh[i] < kstar );
            if( adv_dom ) ++fired_total;
          }
        }

        // ---- cache the read (always) ----
        t_DirRead rd;
        rd.prim_col = (int)jp;  rd.dir_idx = (int)i;
        rd.order    = ord;      rd.order_gt2 = order_gt2;  rd.via_aux = via_aux;
        rd.adv  = adv;   rd.diff = diff;
        rd.kstar = kstar; rd.eps = eps; rd.kh = kh[i];
        rd.advection_dominated = adv_dom;
        reads.push_back( rd );

        // ---- diagnostic print (gated; output format unchanged) ----
        if( kPrint ){
          if( order_gt2 ){
            _disp( 3 ) << "    " << std::left << std::setw(16) << sym.vState[jp].name()
                      << std::setw(4) << sym.vDom[i].name() << std::right
                      << "  order>2: not handled (probe is order<=2)\n";
            continue;
          }
          IT primary = ( ord >= 2 ) ? Options::IC_FLUX
                                    : ( hyp ? Options::IC_UPWIND : Options::IC_VALUE );
          IT final_t = primary;
          const char* active = "-";
          const char* corr   = "(n/a)";
          if( ord >= 2 ){
            if( std::isfinite( kh[i] ) && std::isfinite( kstar ) ){
              bool const top_active = ( kh[i] >= kstar );
              active = top_active ? "yes" : "no";
              if( !top_active ){
                final_t = hyp ? Options::IC_UPWIND : Options::IC_VALUE;
                corr = "->VALUE";
              } else corr = "(inert)";
            } else { active = "?"; corr = "(no k_h)"; }
          }

          IT now = resolved_interface_type( bid, sym.vDom[i], EqnRole::INTERIOR, false );

          _disp( 3 ) << "    " << std::left << std::setw(16) << sym.vState[jp].name()
                    << std::setw(4) << sym.vDom[i].name() << std::right
                    << "  " << ord << "  "
                    << std::left << std::setw(7) << it_name( primary ) << std::right
                    << std::scientific << std::setprecision(2)
                    << "  " << std::setw(9) << eps
                    << "  " << std::setw(9) << kstar
                    << "  " << std::setw(9) << kh[i]
                    << std::defaultfloat
                    << "  " << std::setw(3) << active << " "
                    << std::left << std::setw(9) << corr << std::right
                    << " " << std::left << std::setw(7) << it_name( final_t ) << std::right
                    << " " << it_name( now )
                    << ( ( primary != final_t ) ? "  <FIRED>" : "" )
                    << "\n";
        }
      }
    }
    } catch( ... ) {
      _blockInterfaceRead.erase( bid );
      if( kPrint )
        _disp( 3 ) << "  blk" << bid << ": interface read failed; skipped\n";
    }
  }

  if( kPrint )
    _disp( 3 ) << "OCFESLV::_symbol_balance_probe ** summary: corrector FIRED on "
              << fired_total << " of " << ord2_total
              << " order-2 (state,dir); inert elsewhere\n";
}

inline void
OCFEBASE::_classification_reference
( std::vector<double>& state_ref,
  std::vector<double>& input_ref,
  std::vector<double>& cst_ref,
  std::vector<double>& dom_ref )
{
  // rev207: the model-level reference needs no mesh; the knob lets the corpus measure whether any
  // classification verdict actually moves (they agree exactly wherever references are scalar).
  static bool const kModelRef = false;   // CRONOS_MODEL_REF (retired 2026-10-07, WORKPLAN 3.B batch 2b)
  if( kModelRef ){
    FFModel::_classification_reference( state_ref, input_ref, cst_ref, dom_ref );
    return;
  }

  // Build full collocated reference arrays first.  This makes scalar
  // classification references consistent with domain-dependent references and
  // with auxiliary states introduced by reduce_order().  Auxiliaries are
  // evaluated from their stored defining expressions by _init_reference().
  std::vector<double> var_full( _nCollVar, 0. );
  std::vector<double> inp_full( _nCollInp, 0. );
  _init_reference( var_full.empty()? nullptr: var_full.data(),
                   inp_full.empty()? nullptr: inp_full.data(),
                   _vCstVal.empty()? nullptr: _vCstVal.data() );

  state_ref.assign( _vVar.size(), 0. );
  for( size_t i = 0; i < _vVar.size(); ++i ){
    auto itm = _mVar.find( _vVar[i] );
    auto itb = _blkVar.find( _vVar[i] );
    if( itm != _mVar.end() && itb != _blkVar.end() )
      state_ref[i] = _sample_reference( _vVar[i], itm->second, itb->second,
                                        var_full.empty()? nullptr: var_full.data() );
  }

  input_ref.assign( _vInp.size(), 0. );
  for( size_t i = 0; i < _vInp.size(); ++i ){
    auto itm = _mInp.find( _vInp[i] );
    auto itb = _blkInp.find( _vInp[i] );
    if( itm != _mInp.end() && itb != _blkInp.end() )
      input_ref[i] = _sample_reference( _vInp[i], itm->second, itb->second,
                                        inp_full.empty()? nullptr: inp_full.data() );
  }

  cst_ref.assign( _vCst.size(), 0. );
  if( !_vCstVal.empty() ){
    if( _vCstVal.size() != _vCst.size() )
      throw Exceptions( Exceptions::CSTVAL );
    cst_ref = _vCstVal;
  }

  dom_ref.assign( _vDom.size(), 0. );
  for( size_t i = 0; i < _vDom.size(); ++i )
    dom_ref[i] = _default_dom_ref( _vDom[i] );
}

inline bool
OCFEBASE::_validate_hyperbolic_incoming_bcs()
{
  bool any_hyp = false;
  for( auto const& bc : _blockClassification )
    if( bc.second.evolution_hyperbolic ){ any_hyp = true; break; }
  if( !any_hyp ) return true;

  std::vector<double> state_ref, input_ref, cst_ref, dom_ref;
  _classification_reference( state_ref, input_ref, cst_ref, dom_ref );

  auto value_by_id = []( FFVar const& v, auto const& vars,
                         std::vector<double> const& vals, size_t fb ) -> double {
    for( size_t k = 0; k < vars.size() && k < vals.size(); ++k )
      if( vars[k].id() == v.id() ) return vals[k];
    return fb < vals.size() ? vals[fb] : 0.;
  };

  bool ok = true;
  for( auto const& bc : _blockClassification ){
    int             const  bid = bc.first;
    t_Classify      const& cls = bc.second;
    if( !cls.evolution_hyperbolic ) continue;
    auto itsym = _blockSymbol.find( bid );
    auto itfd  = _blockFaceData.find( bid );
    if( itsym == _blockSymbol.end() || itfd == _blockFaceData.end() ) continue;
    t_Symbol const& sym    = itsym->second;
    auto const&     fdata  = itfd->second;
    size_t const    nState = sym.vState.size();
    if( nState == 0 ) continue;
    size_t const    evo_idx = cls.evolution_dom_idx;

    // Reference evaluation context (states + dom + inputs + consts), keyed by id.
    std::vector<FFVar>  eval_vars;
    std::vector<double> eval_vals;
    for( size_t j = 0; j < nState; ++j ){
      eval_vars.push_back( sym.vState[j] );
      eval_vals.push_back( value_by_id( sym.vState[j], _vVar, state_ref, j ) );
    }
    // provision non-principal coefficient states (e.g. wetting omega) absent from sym.vState
    for( size_t p = 0; p < _vVar.size(); ++p ){
      bool principal = false;
      for( size_t j = 0; j < nState; ++j )
        if( sym.vState[j].id() == _vVar[p].id() ){ principal = true; break; }
      if( principal ) continue;
      eval_vars.push_back( _vVar[p] );
      eval_vals.push_back( value_by_id( _vVar[p], _vVar, state_ref, p ) );
    }
    for( size_t i = 0; i < sym.vDom.size(); ++i ){
      eval_vars.push_back( sym.vDom[i] );
      eval_vals.push_back( value_by_id( sym.vDom[i], _vDom, dom_ref, i ) );
    }
    for( size_t p = 0; p < _vInp.size(); ++p ){
      eval_vars.push_back( _vInp[p] );
      eval_vals.push_back( value_by_id( _vInp[p], _vInp, input_ref, p ) );
    }
    for( size_t c = 0; c < _vCst.size(); ++c ){
      eval_vars.push_back( _vCst[c] );
      eval_vals.push_back( c < cst_ref.size() ? cst_ref[c] : 0. );
    }

    for( auto const& fd : fdata ){
      if( !fd.pdom_var || fd.dom_idx == evo_idx ) continue;
      if( fd.Aplus.n_rows  < (arma::uword)nState ||
          fd.Aminus.n_rows < (arma::uword)nState ) continue;

      struct EndSpec { int face; arma::mat const* A; char const* name; };
      EndSpec const ends[2] = {
        { FFDom::LB, &fd.Aplus,  "z=LB" },   // incoming @LB = right-going (A+)
        { FFDom::UB, &fd.Aminus, "z=UB" } }; // incoming @UB = left-going  (A-)

      for( auto const& end : ends ){
        // Incoming characteristic coupling space = ROW space of the incoming
        // projector = leading RIGHT-singular vectors (V columns).  A boundary
        // condition prescribes the characteristic w = l.(states), whose state
        // coupling d(BC)/d(states) is the LEFT eigenvector l = a row of A+/A-;
        // it therefore lives in rowspace(A), spanned by V -- NOT colspace
        // (U = right eigenvectors).  For a SYMMETRIC principal symbol the two
        // coincide (PDE14), so U.cols formerly sufficed; for a NON-SYMMETRIC
        // block (PDE15 / PSA pressure-velocity) left != right and only V is
        // correct -- this matches the auto-closure, which is also built from
        // rowspace(A) = V columns.
        arma::mat U, V; arma::vec s;
        if( !arma::svd( U, s, V, *end.A ) || s.is_empty() ) continue;
        double const tol = 1e-8 * std::max( 1.0, s(0) );
        arma::uword r_in = 0;
        for( arma::uword i = 0; i < s.n_elem; ++i ) if( s(i) > tol ) ++r_in;
        arma::mat const Vin = ( r_in > 0 ) ? V.cols( 0, r_in - 1 )
                                           : arma::mat( (arma::uword)nState, 0 );

        // User BOUNDARY rows collocated at this block / face direction / end.
        std::vector<FFVar> bcs;
        for( auto const& eqn : _mEqn ){
          if( eqn.opt->role != EqnRole::BOUNDARY || eqn.opt->block_id != bid ) continue;
          auto it = eqn.dom.find( *fd.pdom_var );
          if( it == eqn.dom.end() || it->second != end.face ) continue;
          bcs.push_back( eqn.var );
        }

        if( bcs.size() != (size_t)r_in ){
          _disp( 2 ) << "OCFESLV::setup ** hyperbolic BC guard: block " << bid
                    << " face " << end.name << " needs " << r_in
                    << " incoming condition(s) but " << bcs.size()
                    << " BOUNDARY row(s) are collocated there\n";
          ok = false;
          continue;
        }

        // 4.0b-i (2026-10-09): the test is the KREISS condition, not membership of Vin.  Writing the states in
        // characteristic variables, u = R w, a condition a.u = g reads a.R_in w_in + a.R_out w_out = g: it pins the
        // incoming w_in iff a.R_in != 0 -- a NONZERO component on the incoming RIGHT eigenvectors, the COLUMN space
        // of A+/A- (the U columns of this SVD).  Membership of Vin (the LEFT eigenvectors: "the condition is a
        // function of w_in only") is sufficient but stricter: it REJECTED u = g at z = 0 for u_t + v_z = 0,
        // v_t + u_z = 0 (incoming mode u + v; a = (1,0) is outside span{(1,1)} yet pins it, u - v being known from
        // the interior).  The two agree for scalars and decoupled systems -- all the corpus had.  PDE16's
        // w_out + beta d_z w_out (a = (2,-1), r_in = (1,2)) stays rejected: a.r_in = 0.
        arma::mat Qin;  if( r_in > 0 ) Qin = arma::orth( U.cols( 0, r_in - 1 ) );
        // fraction of a coupling vector that lies in the incoming space: 1 inside Vin, 0 orthogonal to it;
        // <0 for a trivial (zero) coupling = "no constraint"
        auto in_frac = [&]( arma::vec const& w )->double {
          double const wn = arma::norm( w );
          if( wn < 1e-30 ) return -1.0;          // trivial: imposes nothing
          if( r_in == 0 )  return 0.0;           // nothing incoming -> touches nothing
          return arma::norm( Qin.t() * w ) / wn;
        };
        std::vector<arma::rowvec> face_rows;     // the accepted couplings restricted to Vin (the Kreiss matrix)

        for( auto const& bcvar : bcs ){
          // Classify the BC's structure: FFPartial = LOCAL derivative (flux BC,
          // refineable via the principal-symbol proxies); FFIntegral = NONLOCAL
          // (no local characteristic meaning -- count stays enforced, skip).
          bool has_partial = false, has_integral = false;
          auto sg = _dag->subgraph( 1, &bcvar );
          for( auto const& op : sg.l_op ){
            if( op->sameid( typeid(FFPartial)  ) ) has_partial  = true;
            if( op->sameid( typeid(FFIntegral) ) ) has_integral = true;
          }
          if( has_integral ){
            if( _solve_display_level() >= 1 )
              _disp( 2 ) << "OCFESLV::setup ** hyperbolic BC guard: block " << bid
                        << " face " << end.name
                        << " BC is nonlocal (FFIntegral); direction test skipped\n";
            continue;
          }

          // VALUE coupling a = d(BC)/d(state-values).  For a flux BC we first
          // proxy the surviving face/transverse OpP(u_j,{x_i,1}) nodes with the
          // block symbol's derivative proxies v_{u_j_x_i}=vDerivProxy[i*nState+j]
          // so FAD becomes well-defined; the proxies are seeded to a nonzero
          // nominal (1.0), which is exact in DIRECTION for the common single-
          // derivative flux (Darcy v=d_z P) because the test is scale-invariant.
          FFVar bc_use = bcvar;
          std::vector<FFVar>  ev_vars = eval_vars;
          std::vector<double> ev_vals = eval_vals;
          bool const proxies_ok =
            ( sym.vDerivProxy.size() == sym.vDom.size() * nState );
          if( has_partial && proxies_ok ){
            std::vector<FFVar> sub_targ, sub_repl;
            for( auto const& op : sg.l_op ){
              if( !op->sameid( typeid(FFPartial) ) ) continue;
              auto const* pop = mc::type_cast<FFPartial const>( op );
              for( size_t jj = 0; jj < op->varin.size(); ++jj ){
                FFVar const* operand = op->varin[jj];
                FFVar const* dag_out = op->varout[jj];
                size_t j = nState;
                for( size_t s = 0; s < nState; ++s )
                  if( sym.vState[s].id() == operand->id() ){ j = s; break; }
                if( j == nState ) continue;          // derivative of a non-state
                for( auto const& [indep_var, ord] : pop->Indep().expr ){
                  if( ord != 1 ) continue;           // only clean first-order
                  size_t i = sym.vDom.size();
                  for( size_t d = 0; d < sym.vDom.size(); ++d )
                    if( sym.vDom[d].id() == indep_var.id() ){ i = d; break; }
                  if( i == sym.vDom.size() ) continue;
                  sub_targ.push_back( *dag_out );
                  sub_repl.push_back( sym.vDerivProxy[ i*nState + j ] );
                }
              }
            }
            if( !sub_targ.empty() ){
              auto subd = _dag->substitute( std::vector<FFVar>{ bcvar },
                                            sub_targ, sub_repl );
              if( !subd.empty() ) bc_use = subd.front();
            }
            for( auto const& p : sym.vDerivProxy ){
              ev_vars.push_back( p );
              ev_vals.push_back( 1.0 );            // nonzero nominal seed
            }
          }

          auto jacA = _dag->FAD( std::vector<FFVar>{ bc_use }, sym.vState );
          std::vector<double> aval( nState, 0. );
          _dag->eval( jacA, aval, ev_vars, ev_vals );
          arma::vec a( (arma::uword)nState );
          for( size_t k = 0; k < nState; ++k ) a((arma::uword)k) = aval[k];
          double const frac_a = in_frac( a );
          double frac_b = -1.0;           // face-normal DERIVATIVE coupling: its fraction in Vin
          arma::vec b_used;
                                          // (criterion #3); stays trivial (-1) unless
                                          // a flux BC's b is computed below.

          // Criterion #3 (value OR derivative).  A BC pins an incoming mode through
          // its VALUE coupling a=d(BC)/d(u) OR its face-normal DERIVATIVE coupling
          // b=d(BC)/d(d_z u); accept if EITHER lies in Vin.  A pure value test
          // (a in Vin only) FALSE-REJECTS a derivative-pinned inflow such as the
          // bilinear c*d_z P (PDE18 B): there a=(0,1) touches the OUTGOING value but
          // b=(c_ref,0) pins the INCOMING characteristic through d_z P -- the
          // spectral boundary derivative is a full nodal combination that includes
          // the boundary value with nonzero weight, so a derivative coupling on an
          // incoming characteristic is a legitimate pin.  This is the MIRROR of the
          // removed b-test: b RESCUES (OR), it is NOT a REQUIREMENT (AND).  Reject
          // only when NEITHER coupling touches Vin -- e.g. PDE16 D (w_out+b d_z w_out:
          // a and b both outgoing) or PDE18 B2 (c*d_z c: pure outgoing).  The proxies
          // (seeded 1.0) are already in ev_vars, so b is FAD-able w.r.t. the
          // FACE-NORMAL derivative proxies vDerivProxy[i_face*nState + j].
          if( has_partial && proxies_ok ){
            size_t i_face = sym.vDom.size();
            for( size_t d = 0; d < sym.vDom.size(); ++d )
              if( sym.vDom[d].id() == fd.pdom_var->id() ){ i_face = d; break; }
            if( i_face < sym.vDom.size() ){
              std::vector<FFVar> face_proxies;
              for( size_t j = 0; j < nState; ++j )
                face_proxies.push_back( sym.vDerivProxy[ i_face*nState + j ] );
              auto jacB = _dag->FAD( std::vector<FFVar>{ bc_use }, face_proxies );
              std::vector<double> bval( nState, 0. );
              _dag->eval( jacB, bval, ev_vars, ev_vals );
              b_used.set_size( (arma::uword)nState );
              for( size_t k = 0; k < nState; ++k ) b_used((arma::uword)k) = bval[k];
              frac_b = in_frac( b_used );
            }
          }
          if( _solve_display_level() >= 2 && has_partial )
            _disp( 2 ) << "OCFESLV::setup ** hyperbolic flux-BC guard: block " << bid
                      << " face " << end.name << " (flux BC) value-in="
                      << ( frac_a < 0 ? 0.0 : frac_a ) << " deriv-in="
                      << ( frac_b < 0 ? 0.0 : frac_b ) << "\n";

          // Criterion #3: imposes nothing if BOTH couplings are trivial; accept if EITHER the value or the
          // derivative coupling TOUCHES Vin (a nonzero component on it); reject only when the BC constrains
          // something yet touches no incoming characteristic in either.
          bool const a_triv = ( frac_a < 0 ), b_triv = ( frac_b < 0 );
          if( a_triv && b_triv ) continue;             // trivial: imposes nothing
          bool const a_in = ( !a_triv && frac_a >= 1e-6 );
          bool const b_in = ( !b_triv && frac_b >= 1e-6 );
          if( a_in || b_in ){
            arma::vec const& w_used = a_in ? a : b_used;
            face_rows.push_back( ( Qin.t() * w_used ).t() / arma::norm( w_used ) );
          }
          if( !a_in && !b_in ){
            _disp( 2 ) << "OCFESLV::setup ** hyperbolic BC guard: block " << bid
                      << " face " << end.name
                      << " incoming BC touches no incoming characteristic in value "
                      << "OR derivative (value-in " << ( frac_a < 0 ? 0.0 : frac_a )
                      << ", deriv-in " << ( frac_b < 0 ? 0.0 : frac_b )
                      << ") -- it prescribes an OUTGOING characteristic; check the "
                      << "flow-direction / which end this condition belongs on\n";
            ok = false;
            _setupStatus = SetupStatus::HYP_BC_MISDIRECTED;   // distinct from the count mismatch
          }
        }
        // Face level: the Kreiss matrix (the accepted couplings restricted to Vin, one row each) must be
        // nonsingular -- two conditions each touching Vin may still pin the SAME incoming characteristic.
        if( ok && r_in > 0 && face_rows.size() == (size_t)r_in ){
          arma::mat K( (arma::uword)r_in, (arma::uword)r_in );
          for( size_t k = 0; k < face_rows.size(); ++k ) K.row( (arma::uword)k ) = face_rows[k];
          arma::vec const sv = arma::svd( K );
          if( sv.n_elem && sv.min() <= 1e-8 * sv.max() ){
            _disp( 2 ) << "OCFESLV::setup ** hyperbolic BC guard: block " << bid << " face " << end.name
                      << ": its " << r_in << " incoming condition(s) pin the same incoming characteristic(s)"
                      << " (Kreiss matrix singular: sigma_min/sigma_max = " << sv.min() / sv.max() << ")\n";
            ok = false;
            _setupStatus = SetupStatus::HYP_BC_MISDIRECTED;
          }
        }
      }
    }
  }
  return ok;
}

// ======================================================================
// OCFESLV::_ensure_sat_coupling / _sat_coupling
// ======================================================================
//
// Per-equation SAT coupling cache.  Faithfully reproduces the inline triple-
// nested lookup that the legacy scalar SAT path performed once per (equation,
// state, direction) on every eval, but precomputes it into a flat per-
// equation vector at most once per classification change.  When the model is
// not classified, the cache stays empty and _sat_coupling() returns the
// identity coupling 1.0 (the same fallback the legacy inline code used).
//
inline void
OCFEBASE::_ensure_sat_coupling
()
const
{
  if( _satCouplingBuilt
   && _satCouplingClassified == _classified
   && _satCouplingSerial == _classificationSerial )
    return;

  _satCoupling.clear();
  _satCouplingClassified = _classified;
  _satCouplingSerial     = _classificationSerial;
  _satCouplingBuilt      = true;

  if( !_classified ) return;

  _satCoupling.resize( _mEqn.size() );

  for( size_t ndx = 0; ndx < _mEqn.size(); ++ndx ){
    auto const& eqnrec = _mEqn[ndx];
    int const   bid    = eqnrec.opt->block_id;

    auto itsym = _blockSymbol.find( bid );
    auto itcls = _blockClassification.find( bid );
    t_Symbol   const& sym_c = ( itsym != _blockSymbol.end() )
                            ? itsym->second : _symbol;
    t_Classify const& cls_c = ( itcls != _blockClassification.end() )
                            ? itcls->second : _classification;

    // Row index of this equation in the block symbol.  An equation that did
    // not participate in classification (INITIAL/BOUNDARY/DIAGNOSTIC) is not
    // a row of vEqn -> no entries -> 1.0 fallback in the lookup.
    size_t k_eq = sym_c.vEqn.size();
    for( size_t k = 0; k < sym_c.vEqn.size(); ++k )
      if( sym_c.vEqn[k].id() == eqnrec.var.id() ){ k_eq = k; break; }
    if( k_eq == sym_c.vEqn.size() ) continue;

    // Pre-resolve symbol state columns to _vVar positions once per equation.
    size_t const nState = sym_c.vState.size();
    std::vector<size_t> state_pos( nState, std::numeric_limits<size_t>::max() );
    for( size_t j = 0; j < nState; ++j ){
      auto it = _mDep.find( sym_c.vState[j].id().second );
      if( it != _mDep.end() && it->second.first == STATE )
        state_pos[j] = it->second.second;
    }

    auto& rec = _satCoupling[ndx];
    rec.reserve( nState * sym_c.vDom.size() );
    for( size_t i_d = 0; i_d < sym_c.vDom.size() && i_d < cls_c.Ai.size(); ++i_d ){
      arma::mat const& A = cls_c.Ai[i_d];
      if( (arma::uword)k_eq >= A.n_rows ) continue;
      long const dom_id = sym_c.vDom[i_d].id().second;
      for( size_t j = 0; j < nState; ++j ){
        if( state_pos[j] == std::numeric_limits<size_t>::max() ) continue;
        if( (arma::uword)j >= A.n_cols ) continue;
        rec.push_back( { state_pos[j], dom_id,
                         A( (arma::uword)k_eq, (arma::uword)j ) } );
      }
    }
  }
}

inline double
OCFEBASE::_sat_coupling
( size_t const ndx, size_t const state_pos, FFVar const& dom )
const
{
  FFVar const* pstate_sat = state_pos < _vVar.size() ? &_vVar[state_pos] : nullptr;
  auto structural_derivative_coupling = [&]() -> double
  {
    if( !_dag || ndx >= _mEqn.size() || !pstate_sat ) return 0.0;
    FFVar const& state = *pstate_sat;
    FFVar const& eqn   = _mEqn[ndx].var;

    auto expr_contains_state = [&]( FFVar const& expr ) -> bool
    {
      if( expr.id() == state.id() ) return true;
      FFSubgraph sg = _dag->subgraph( 1, &expr );
      for( auto const& op : sg.l_op ){
        for( auto const* vin : op->varin )
          if( pstate_sat && vin && vin->id() == pstate_sat->id() ) return true;
        for( auto const* vout : op->varout )
          if( pstate_sat && vout && vout->id() == pstate_sat->id() ) return true;
      }
      return false;
    };

    FFSubgraph sg = _dag->subgraph( 1, &eqn );
    for( auto const& op : sg.l_op ){
      if( !op->sameid( typeid(FFPartial) ) ) continue;
      auto const* pop = mc::type_cast<FFPartial const>( op );
      if( !pop ) continue;

      bool dom_match = false;
      for( auto const& kv : pop->Indep().expr ){
        if( kv.first.id() == dom.id() ){ dom_match = true; break; }
      }
      if( !dom_match ) continue;

      for( auto const& operand : pop->Var() )
        if( expr_contains_state( operand ) ) return 1.0;

      // Legacy fallback: some old external-operation nodes may not carry the
      // FFPartial::_Var metadata after graph copy; compare by FFVar::id().
      for( auto const* operand : op->varin )
        if( operand && operand->id() == state.id() ) return 1.0;
    }
    return 0.0;
  };

  if( !_classified || ndx >= _satCoupling.size() )
    return structural_derivative_coupling();

  long const dom_id = dom.id().second;
  for( auto const& r : _satCoupling[ndx] )
    if( r.state_pos == state_pos && r.dom_id == dom_id ) return r.value;

  // Reduced-order LINK rows are auxiliary definition identities such as
  // q - D*u_x = 0.  Classification may give the block a non-empty principal
  // symbol and thereby suppress the ordinary structural fallback, but the
  // reduced interface plan still needs LINK rows to receive primitive C0
  // trace multipliers in the reduced normal direction.  Use structural
  // derivative incidence for LINK rows in any direction with a reduced aux
  // link in the block; this covers both RED_MAIN and RED_FULL.
  if( ndx < _mEqn.size()
   && _mEqn[ndx].opt->role == EqnRole::LINK
   && _block_has_aux_link_in_direction( _mEqn[ndx].opt->block_id, dom ) ){
    // Prefer the defining auxiliary expression over the copied LINK residual:
    // the latter may have been simplified/substituted during RED_FULL, while
    // _auxDef records the primitive derivative that the LINK represents.
    for( auto const& aux : _auxDef ){
      if( aux.block_id != _mEqn[ndx].opt->block_id ) continue;
      bool normal = false;
      for( auto const& d : aux.diff_dom )
        if( d.id() == dom.id() ){ normal = true; break; }
      if( !normal ) continue;
      if( pstate_sat && aux.expr.id() == pstate_sat->id() ) return 1.0;
      if( _dag && aux.expr.dag() ){
        FFSubgraph sg_aux = _dag->subgraph( 1, &aux.expr );
        for( auto const& op : sg_aux.l_op ){
          for( auto const* vin : op->varin )
            if( pstate_sat && vin && vin->id() == pstate_sat->id() ) return 1.0;
          for( auto const* vout : op->varout )
            if( pstate_sat && vout && vout->id() == pstate_sat->id() ) return 1.0;
        }
      }
    }
    return structural_derivative_coupling();
  }

  // If this row has classified principal-symbol entries, absence of this
  // state/domain pair means zero coupling.  For rows outside the classified
  // symbol, fall back to the structural derivative-incidence test above.
  if( ndx < _satCoupling.size() && !_satCoupling[ndx].empty() ) return 0.0;
  return structural_derivative_coupling();
}

inline std::vector<std::vector<double>>
OCFEBASE::node_colloc
( FFVar const& var )
const
{
  if( !_issetup ) throw Exceptions( Exceptions::SETUP );
  std::vector<std::vector<double>> vnodes;

  // Search variable var dependencies
  auto const& [typ, ndx] = _mDep.at(var.id().second);
  std::set< FFVar, lt_FFVar > dep;
  switch( typ ){
    case STATE:
      dep = _mVar.at(var);
      break;
    case INPUT:
      dep = _mInp.at(var);
      break;
    case DOMAIN:
      dep = {var};
      break;
    default:
      return vnodes;
  }
  size_t const ndep = dep.size();

  // Loop over all relevant subdomain blocks of current variable
  std::vector<std::vector<double>> anodes( ndep );
  std::map<FFVar, size_t, lt_FFVar> ndx_el;
  size_t pos = 0;
  for( ; _next_block( ndx_el, dep ) && !ndx_el.empty(); ){

#ifdef MC__OCFESLV_EVAL_DEBUG
    _disp( 3 ) << "BLOCK";
    for( auto const& var : dep )
      _disp( 3 ) << " " << var << "=" << ndx_el[var];
    _disp( 3 ) << std::endl;
#endif

    // Fill-in nodes for all variables in current block 
    size_t nnodes = 1, idep = 0;
    for( auto const& v : dep ){
      std::vector<double> lgnodes;
      if( typ == INPUT )
        lgnodes = _node_block( _input_dom( var, v ), ndx_el.at( v ) );
      else
        lgnodes = _node_block( v, ndx_el );
      nnodes *= lgnodes.size();
      anodes[idep++] = std::move( lgnodes );
    }

    // Flatten anodes as a vector of coordinate vectors
    // First dimension varies fastest
    vnodes.resize( pos+nnodes, std::vector<double>( ndep ) );
    std::vector<size_t> idx( ndep, 0 );
    for( size_t inode = 0; inode < nnodes; ++inode ){

      // Current coordinate vector
      for( size_t i=0; i<ndep; ++i )
        vnodes[pos+inode][i] = anodes[i][idx[i]];

      // Increment mixed-radix counter
      for( size_t i=0; i<ndep; ++i ){
        if( ++idx[i] < anodes[i].size() ) break;
        idx[i] = 0;
      }
    }
    
    // Increment position
    pos += nnodes;
  }

  return vnodes;
}

inline OCFEBASE::DuplicateNodeSpread
OCFEBASE::duplicate_node_spread
( FFVar const& st, std::vector<double> const& var, size_t off )
const
{
  // rev44: derive the base offset unless the caller supplied a verified one.
  if( off == std::numeric_limits<size_t>::max() )
    off = _state_base_offset( st );

  struct Accum { double lo; double hi; size_t count; };
  std::map<std::vector<long long>,Accum> groups;

  std::vector<std::vector<double>> const nodes = node_colloc( st );
  for( size_t i = 0; i < nodes.size(); ++i ){
    if( off + i >= var.size() ) break;               // short vector: report what exists
    std::vector<long long> key;
    key.reserve( nodes[i].size() );
    for( double c : nodes[i] )
      key.push_back( static_cast<long long>( std::llround( c * 1.0e12 ) ) );
    double const v = var[off+i];
    auto it = groups.find( key );
    if( it == groups.end() ) groups.emplace( std::move(key), Accum{ v, v, 1 } );
    else{
      it->second.lo = std::min( it->second.lo, v );
      it->second.hi = std::max( it->second.hi, v );
      ++it->second.count;
    }
  }

  DuplicateNodeSpread out;
  for( auto const& kv : groups ){
    Accum const& g = kv.second;
    out.max_multiplicity = std::max( out.max_multiplicity, g.count );
    if( g.count >= 2 ){
      ++out.n_group;
      out.max_pair_spread = std::max( out.max_pair_spread, g.hi - g.lo );
    }
    if( g.count >= 4 )
      out.max_corner_spread = std::max( out.max_corner_spread, g.hi - g.lo );
  }
  return out;
}

inline double
OCFEBASE::report_duplicate_node_spreads
( std::vector<double> const& var, char const* label, std::ostream* os )
const
{
  if( !_issetup ) throw Exceptions( Exceptions::SETUP );

  double max_pair = 0.0, max_corner = 0.0;
  size_t n_group_tot = 0, max_mult = 0;
  std::string worst_state;

  std::vector<std::pair<std::string,DuplicateNodeSpread>> per_state;
  for( auto const& st : states_colloc() ){
    // rev44 -- OFFSET BUG FIX.  This walked the states accumulating
    //     off += node_colloc( st ).size();
    // which is WRONG whenever any state is SCALAR (declared with an empty domain list).
    // node_colloc() returns an EMPTY vector for such a state, so the walk adds 0 -- but the
    // state still owns one var[] slot.  Every state after the first scalar one was therefore
    // read at an offset short by the number of preceding scalar states, and the "spread"
    // differenced two unrelated entries.
    //
    // HOW IT WAS CAUGHT, because the failure was not obvious.  OCFE_scalar_nest reported
    // Dz_a(z) spread = 2.500000e-01 at the interface z=0.5, bit-identical under IC_WEAK,
    // IC_TRACE and IC_STRONG.  a(z)=z(1-z) is quadratic and exactly representable in
    // degree-5 CGL, so a'(0.5)=0 must be exact from both sides and the spread had to be
    // zero.  Dumping the individual node values showed 9.0450849719e-01 and
    // 6.5450849719e-01 -- which are a'(z)=1-2z at z=0.04774575 and z=0.17274575, the first
    // two nodes of element 0.  Four scalar states (p1..p4) preceded Dz_a, so the read was
    // short by exactly 4.  Independently confirmed by Dpy_Intz_a, whose two "copies" read
    // 1/24 and 1/36 -- the values of the scalar states Evy_Intz_a and Inty_Intz_a.
    //
    // pos_state( st, {} ) is the authoritative base offset: it is what _init_var_block()
    // writes through for a scalar state (data[pos_state(v,{})] = ...), so it accounts for
    // every state regardless of dimensionality.  Do NOT reintroduce a manual walk.
    //
    // The three driver-side implementations of this measurement predate the header version
    // and use the SAME manual walk, so their numbers agree with the old header path by
    // sharing its error, not by confirming it.  They must be migrated to this call.
    size_t const off = _state_base_offset( st );
    DuplicateNodeSpread const d = duplicate_node_spread( st, var, off );
    if( d.max_pair_spread > max_pair ){
      max_pair = d.max_pair_spread; worst_state = st.name();
    }
    max_corner  = std::max( max_corner, d.max_corner_spread );
    max_mult    = std::max( max_mult,   d.max_multiplicity );
    n_group_tot += d.n_group;
    per_state.emplace_back( st.name(), d );
  }

  if( os ){
    *os << "  [dup-spread] " << ( label && *label ? label : "(unlabelled)" )
        << ": states=" << per_state.size()
        << " groups=" << n_group_tot
        << " max_multiplicity=" << max_mult << "\n";
    std::ios_base::fmtflags const _f = os->flags();
    std::streamsize const _p = os->precision();
    for( auto const& ps : per_state )
      *os << "  [dup-spread]   " << std::setw(18) << ps.first
          << " max_pair="   << std::scientific << std::setprecision(6)
          << ps.second.max_pair_spread
          << " max_corner=" << ps.second.max_corner_spread
          << " max_mult="   << ps.second.max_multiplicity << "\n";
    os->flags( _f ); os->precision( _p );
    *os << "  [dup-spread] WORST max_pair=" << std::scientific << std::setprecision(6)
        << max_pair << " (" << ( worst_state.empty() ? "n/a" : worst_state ) << ")"
        << "  max_corner=" << max_corner << std::defaultfloat << std::endl;
  }

  return max_pair;
}

inline bool

//==========================================================================================================
// INTERFACE PLAN -- claims -> realisations -> rows and columns.  The DECISIONS live in ocplan.hpp; what is here
// is the assembly: which claims exist, which rows receive them, and how the chosen realisation becomes
// matrix structure.
// (banner placed from the definitions themselves: this region is lines 14475-15202 of the pre-rev268 file)
//==========================================================================================================

OCFEBASE::_suppress_tensor_interface_corner
( Options::InterfaceType type, bool exact_reduced_weak, int block_id,
  FFVar const& face_dom,
  std::vector<t_InterfaceTensorCornerDim> const& od,
  FFVar const* recv_state )
const
{
  // One shared predicate for both setup-time donor-candidate filtering and
  // LOOP-1 exact interface-row emission.  Keeping these decisions identical is
  // essential: a donor row may be reserved only for an exact row that can later
  // be emitted, and every emitted exact row must have a matching claim in
  // the frozen plan.

  // Initial-face arm: at the t=LB initial face the PRIMITIVE state is pinned by
  // its initial condition, so its cross-direction C0 continuity there is
  // redundant and must be suppressed to avoid an over-determined (singular)
  // augmented system.  A reduced auxiliary, however, has NO initial condition:
  // its t=LB copies are pinned only by their own element's LINK to dU/dx, which
  // differ across the interface, so its continuity is genuine (not redundant)
  // and must be retained.  The connectivity probe confirmed every reduced-aux
  // float lives exactly on this initial face.  Hence skip this arm for auxiliary
  // receiver states.
  bool const recv_is_aux = ( recv_state && _is_auxiliary_state( *recv_state ) );
  if( _block_has_aux_link_in_direction( block_id, face_dom ) && !recv_is_aux ){
    for( auto const& o : od ){
      if( o.iel == 0 && o.inode == 0
       && o.pvar && _block_has_initial_face( block_id, *o.pvar, FFDom::LB ) )
        return true;
    }
  }

  if( exact_reduced_weak ){
    // In IC_WEAK, exact rows are emitted only for reduced-link directions.
    // Suppress a tensor corner only when the preceding direction also emits an
    // exact reduced row.  Ordinary weak/SAT directions do not own exact rows.
    for( auto const& o : od ){
      if( !o.pvar || !_weak_exact_reduced_interface_in_direction( block_id, *o.pvar ) )
        continue;
      if( o.iel > 0 && o.inode == 0 ) return true;
    }
    return false;
  }

  for( auto const& o : od ){
    if( type == Options::IC_VALUE || type == Options::IC_UPWIND ){
      if( o.iel > 0 && o.inode == 0 ) return true;
    }
    else if( type == Options::IC_FLUX ){
      if( o.iel > 0 && o.inode == 0 ) return true;
      if( o.iel + 1 < o.n_elem && o.inode + 1 == o.n_node ) return true;
    }
  }
  return false;
}

inline int
OCFEBASE::_trace_claim_block
( size_t state_id, size_t dom_id ) const
{
  FFVar const* pst = nullptr;
  for( auto const& kv : _mVar ){
    if( static_cast<size_t>( kv.first.id().second ) == state_id ){
      pst = &kv.first;
      break;
    }
  }
  FFVar const* pdom = nullptr;
  for( auto const& kv : _mDom ){
    if( static_cast<size_t>( kv.first.id().second ) == dom_id ){
      pdom = &kv.first;
      break;
    }
  }
  if( !pst || !pdom ) return 0;

  // Prefer a classified block whose symbol contains both the constrained
  // state and the face direction.  If no such block exists, canonicalise to
  // block 0.
  for( auto const& kv : _blockSymbol ){
    bool has_state = false, has_dom = false;
    for( auto const& v : kv.second.vState )
      if( v.id() == pst->id() ){ has_state = true; break; }
    for( auto const& v : kv.second.vDom )
      if( v.id() == pdom->id() ){ has_dom = true; break; }
    if( has_state && has_dom ) return kv.first;
  }
  return 0;
}

inline bool
OCFEBASE::_load_frozen_interface_plan
( std::map<t_InterfaceClaimKey,size_t>& claims ) const
{
  if( !_plan.frozen() || _plan.state() != OCPlan::State::FROZEN ) return false;
  claims.clear();
  for( auto const& kv : _plan.claims_frozen() )
    claims[ _unfreeze_interface_claim_key( kv.first ) ] = kv.second;
  return true;
}

inline PlanInput
OCFEBASE::_pack_plan_input
( std::map<t_InterfaceClaimKey,size_t> const& claims,
  InterfacePlanDraft const& draft ) const
{
  PlanInput in;

  switch( options.INTERFACE.IMPOSITION ){
    case Options::IC_WEAK:   in.options.imposition = PlanOptions::Imposition::WEAK;   break;
    case Options::IC_TRACE:  in.options.imposition = PlanOptions::Imposition::TRACE;  break;
    case Options::IC_STRONG: in.options.imposition = PlanOptions::Imposition::STRONG; break;
  }
  // The rev153c experiment knobs CRONOS_PLAN_REALISATION, _NO_AUX_IMPLIED, _ROOT_HI and _FACE_DESC are retired
  // (2026-10-07, WORKPLAN 3.B batch 1), and with them (4.0, 2026-10-09) the PlanOptions fields they set: the builder
  // hard-codes those defaults -- the plan every sweep since rev156 validated.
  in.options.force_bad_keep_explicit = Options::_env_flag( "CRONOS_FORCE_BAD_KEEP_EXPLICIT", false );

  in.claims         = claims;
  in.receiver_edges = &draft.receiver_edges;

  in.row_role.reserve( _mEqn.size() );
  for( auto const& e : _mEqn )
    in.row_role.push_back( static_cast<PlanRowRole>( static_cast<int>( e.opt->role ) ) );

  for( auto const& kv : _mVar )
    in.state_name[ static_cast<size_t>( kv.first.id().second ) ] = kv.first.name();
  for( auto const& kv : _mDom )
    in.dom_name[ static_cast<size_t>( kv.first.id().second ) ] = kv.first.name();

  auto pack_tcb = [&]( size_t sid, size_t did ){
    auto const key = std::make_pair( sid, did );
    if( !in.trace_claim_block.count( key ) )
      in.trace_claim_block[ key ] = _trace_claim_block( sid, did );
  };
  for( auto const& kv : claims ) pack_tcb( kv.first.state_id, kv.first.dom_id );
  for( auto const& e : draft.receiver_edges ) pack_tcb( e.claim_key.state_id, e.claim_key.dom_id );

  in.eqn_differentiates_in = _eqn_differentiates_in();

  in.protected_claims = _strong_protected_claims;
  in.promote_explicit = _strong_promoted_claims;   // rev281
  in.dropped_claims   = _strong_dropped_claims;

  in.n_coll_eqn       = draft.n_coll_eqn;
  in.n_coll_var       = draft.n_coll_var;
  in.trace_var_offset = draft.trace_var_offset;

  // ---- rev153a: geometry for the GLOBAL builder ------------------------------------
  // Direction rank: evolution domain (if any) first, then registration order of _mDom.
  { std::map<size_t,int> dir_rank; int r = 1;
    size_t const evo_id = _evolution_dom_var.id().second >= 0
      ? static_cast<size_t>( _evolution_dom_var.id().second ) : std::numeric_limits<size_t>::max();
    for( auto const& kd : _mDom ){
      size_t const did = static_cast<size_t>( kd.first.id().second );
      dir_rank[ did ] = ( did == evo_id ) ? 0 : r++;
    }
    // One DofEdge per distinct claim key (occurrences share geometry).
    std::map<t_InterfaceClaimKey,size_t> edge_index_of_claim;
    std::vector<std::pair<size_t,double>> coeffs;
    for( auto const& kv : claims ){
      PlanInput::DofEdge e; e.key = kv.first;
      auto itr = dir_rank.find( kv.first.dom_id );
      e.direction_rank = ( itr != dir_rank.end() ) ? itr->second : 99;
      e.causal = ( kv.first.dom_id == evo_id );
      if( kv.first.kind == CLAIM_EXACT_C0 && _trace_claim_coeffs( kv.first, coeffs ) && coeffs.size() == 2 ){
        auto const& a = coeffs[0]; auto const& b = coeffs[1];
        if( std::fabs( std::fabs( a.second ) - 1. ) < 1e-12 && std::fabs( std::fabs( b.second ) - 1. ) < 1e-12
         && a.second * b.second < 0. ){
          e.dof_lo = ( a.second > 0. ) ? a.first : b.first;
          e.dof_hi = ( a.second > 0. ) ? b.first : a.first;
          e.nodal  = true;
        }
      }
      _trace_claim_sides( kv.first, e.lo_support, e.hi_support );   // rev155a (both kinds)
      edge_index_of_claim[ kv.first ] = in.dof_edges.size();
      in.dof_edges.push_back( e );
    }
    // Corner-rule-suppressed claims, geometry only (shadow check of the forest).
    for( auto const& key : _cornerSuppressedClaims ){
      if( claims.count( key ) ) continue;
      PlanInput::DofEdge e; e.key = key;
      auto itr = dir_rank.find( key.dom_id );
      e.direction_rank = ( itr != dir_rank.end() ) ? itr->second : 99;
      e.causal = ( key.dom_id == evo_id );
      if( key.kind == CLAIM_EXACT_C0 && _trace_claim_coeffs( key, coeffs ) && coeffs.size() == 2 ){
        auto const& a = coeffs[0]; auto const& b = coeffs[1];
        if( std::fabs( std::fabs( a.second ) - 1. ) < 1e-12 && std::fabs( std::fabs( b.second ) - 1. ) < 1e-12
         && a.second * b.second < 0. ){
          e.dof_lo = ( a.second > 0. ) ? a.first : b.first;
          e.dof_hi = ( a.second > 0. ) ? b.first : a.first;
          e.nodal  = true;
        }
      }
      _trace_claim_sides( key, e.lo_support, e.hi_support );   // rev155a
      in.corner_suppressed_edges.push_back( e );
    }
    // rev155a: pinned-suppressed pairs, packed like corner edges (identification joins only).
    for( auto const& key : _pinSuppressedClaims ){
      if( claims.count( key ) ) continue;
      PlanInput::DofEdge e; e.key = key;
      auto itr = dir_rank.find( key.dom_id );
      e.direction_rank = ( itr != dir_rank.end() ) ? itr->second : 99;
      e.causal = ( key.dom_id == evo_id );
      if( key.kind == CLAIM_EXACT_C0 && _trace_claim_coeffs( key, coeffs ) && coeffs.size() == 2 ){
        auto const& a = coeffs[0]; auto const& b = coeffs[1];
        if( std::fabs( std::fabs( a.second ) - 1. ) < 1e-12 && std::fabs( std::fabs( b.second ) - 1. ) < 1e-12
         && a.second * b.second < 0. ){
          e.dof_lo = ( a.second > 0. ) ? a.first : b.first;
          e.dof_hi = ( a.second > 0. ) ? b.first : a.first;
          e.nodal  = true;
        }
      }
      _trace_claim_sides( key, e.lo_support, e.hi_support );
      in.pinned_edges.push_back( e );
    }
    { // rev155d: column-level pinned pairs (visit-independent), deduplicated on the geometric
      // 5-tuple (recorded keys carry block/eqn ids; synthesized ones do not).
      typedef std::tuple<size_t,size_t,size_t,size_t,int> Geo;
      auto geo = []( t_InterfaceClaimKey const& k ){ return Geo{ k.state_id, k.dom_id, k.iel_lo, k.face, k.kind }; };
      std::set<Geo> have;
      for( auto const& e : in.pinned_edges ) have.insert( geo( e.key ) );
      std::vector<PlanInput::DofEdge> col;
      _pinned_column_pairs( col );
      for( auto& e : col ){
        if( have.count( geo( e.key ) ) ) continue;
        auto itr = dir_rank.find( e.key.dom_id );
        e.direction_rank = ( itr != dir_rank.end() ) ? itr->second : 99;
        e.causal = ( e.key.dom_id == evo_id );
        in.pinned_edges.push_back( e ); have.insert( geo( e.key ) );
      }
    }
    // C2 EXPERIMENT -- IC_WEAK WITH MULTIPLIERS.  CRONOS_WEAK_TAU_RESCUED is an INTEGER, unset or 0 = OFF
    // (rev258 read it as "is set", so =0 turned it ON), and its state is printed whenever it is set:
    //   1  impose exactly every claim with a non-principal receiver edge (the rescue fired; C2 derived or refused)
    //   2  only claims with an edge whose coupling C2 REFUSED (InterfaceReceiverEdge::coupling_undefined)
    // Only IC_WEAK acts on it; IC_TRACE and IC_STRONG already impose every claim exactly.
    { int const kWeakTau = _knob_WEAK_TAU_RESCUED();   // rev276: Options fields, read per setup
      if( kWeakTau ){
        bool const weak = ( options.INTERFACE.IMPOSITION == Options::IC_WEAK );
        bool const kWeakNat = _knob_WEAK_TAU_NATURAL();
        size_t n_edge = 0;
        if( weak ){
          in.options.weak_multipliers = true;
          in.options.weak_natural_only = kWeakNat;
          for( auto const& e : draft.receiver_edges ){
            if( e.source == OCPlan::InterfaceEdgeSource::PRINCIPAL_SYMBOL ) continue;
            if( kWeakTau == 2 && !e.coupling_undefined ) continue;
            in.rescued_claims.insert( e.claim_key ); ++n_edge;
          }
        }
        _disp( 3 ) << "OCFESLV::setup ** PLAN: CRONOS_WEAK_TAU_RESCUED=" << kWeakTau
                  << ( kWeakTau == 2 ? " (undefined-coupling claims)" : " (all rescued claims)" )
                  << " CRONOS_WEAK_TAU_NATURAL=" << ( kWeakNat ? 1 : 0 )
                  << ( kWeakNat ? " (multiplier on natural receivers where a claim has one)" : " (multiplier on every tau-allowed receiver)" );
        if( weak ) _disp( 3 ) << " -- selected " << in.rescued_claims.size() << " claim(s) from " << n_edge << " edge(s)";
        else       _disp( 3 ) << " -- inactive: imposition is not IC_WEAK";
        _disp( 3 ) << std::endl;
      } }

    // Claims the plan never realises (pinned-suppressed), and nodes per element per domain.
    for( auto const& key : _pinSuppressedClaims ) if( !claims.count( key ) ) in.unrealised_claims.insert( key );
    for( auto const& key : _rankImpliedClaims ) in.rank_implied_claims.insert( key );   // rev160
    // rev160a: make the knob's STATE visible in the log.  A sweep that silently ran with the
    // knob unset is indistinguishable from a sweep where the feature is a no-op -- the exact
    // failure this arc already hit once (the census run whose trigger was never set).  Printed
    // only when the knob is ON, so the default path stays byte-identical to rev159.
    {
      static int const kRankImpliedN = 0;   // CRONOS_AUX_IMPLIED_RANK (retired 2026-10-07, WORKPLAN 3.B batch 2b)
      if( kRankImpliedN )
        _disp( 3 ) << "OCFESLV::setup ** PLAN: CRONOS_AUX_IMPLIED_RANK=" << kRankImpliedN
                  << " rank-implied claims collected=" << in.rank_implied_claims.size()
                  << std::endl;
    }
    for( auto const& kd : _mDom ){
      in.dom_n_node[ static_cast<size_t>( kd.first.id().second ) ] = kd.second.n_node;
      in.dom_n_elem[ static_cast<size_t>( kd.first.id().second ) ] = kd.second.n_elem;
    }
    for( auto const& [v, doms] : _mVar ){
      auto& vd = in.state_dims[ static_cast<size_t>( v.id().second ) ];
      for( auto const& d : doms ) vd.push_back( static_cast<size_t>( d.id().second ) );
      std::sort( vd.begin(), vd.end() );
    }
    // Auxiliary definitions (design §3.5).
    for( auto const& a : _auxDef ){
      if( a.parent.id().second < 0 ) continue;
      PlanInput::AuxInfo& ai = in.aux_info[ static_cast<size_t>( a.aux.id().second ) ];
      ai.parent = static_cast<size_t>( a.parent.id().second );
      for( auto const& d : a.diff_dom ) ai.diff_dirs.insert( static_cast<size_t>( d.id().second ) );
    }
    // Each auxiliary's definition row: the LINK row that depends on the auxiliary AND on its parent.  Unique, chains
    // included: the row defining d_y(d_y U) depends on it and on d_y U, not on U.  Every leaf of the row is supplied
    // -- the auxiliaries and parents as independent variables, the rest (parameters, other states, t) as constants --
    // so the dependence evaluation never meets a missing variable.
    if( _dag && !_auxDef.empty() ){
      std::map<long,int> slot;  int ns = 0;
      for( auto const& a : _auxDef ) for( FFVar const* v : { &a.aux, &a.parent } )
        if( v->id().second >= 0 && !slot.count( v->id().second ) ) slot[ v->id().second ] = ns++;
      for( size_t k = 0; k < _mEqn.size(); ++k ){
        if( !_mEqn[k].opt || _mEqn[k].opt->role != EqnRole::LINK ) continue;
        FFSubgraph const sg = _dag->subgraph( 1, &_mEqn[k].var );
        std::vector<FFVar> leaves;  std::vector<FFDep> vals;
        for( FFVar const* v : sg.v_indep ){
          leaves.push_back( *v );  FFDep d;
          auto its = slot.find( v->id().second );  if( its != slot.end() ) d.indep( its->second );
          vals.push_back( d );
        }
        std::vector<FFDep> dep( 1 );
        try{ _dag->eval( std::vector<FFVar>{ _mEqn[k].var }, dep, leaves, vals ); } catch( ... ){ continue; }
        for( auto const& a : _auxDef ){
          if( a.parent.id().second < 0 || !slot.count( a.aux.id().second ) || !slot.count( a.parent.id().second ) ) continue;
          if( dep[0].dep().count( slot[ a.aux.id().second ] ) && dep[0].dep().count( slot[ a.parent.id().second ] ) ){
            in.aux_info[ static_cast<size_t>( a.aux.id().second ) ].def_eqn = k;  break; }
        }
      }
    }
    // Owner rows per state (INTERIOR equations depending on the state or its auxiliaries)
    // and the collocation-set fact per (row, direction).
    if( _dag ){
      std::map<size_t,size_t> parent_of;   // aux id -> primitive id (through chains)
      for( auto const& a : _auxDef ){
        if( a.parent.id().second < 0 ) continue;
        size_t p = static_cast<size_t>( a.parent.id().second );
        auto itp = parent_of.find( p ); if( itp != parent_of.end() ) p = itp->second;
        parent_of[ static_cast<size_t>( a.aux.id().second ) ] = p;
      }
      for( size_t row_id = 0; row_id < _mEqn.size(); ++row_id ){
        if( _mEqn[ row_id ].opt->role != EqnRole::INTERIOR ) continue;
        FFVar const& e = _mEqn[ row_id ].var;
        FFSubgraph sg = _dag->subgraph( 1, &e );
        for( auto const& op : sg.l_op ){
          if( !op || op->type != FFOp::VAR || op->varout.empty() ) continue;
          size_t const vid = static_cast<size_t>( op->varout[0]->id().second );
          auto itp = parent_of.find( vid );
          in.owner_rows[ itp != parent_of.end() ? itp->second : vid ].insert( row_id );
        }
        // Direct differential order of the residual in each direction: the highest
        // Partial order applied to ANY field.  Chain order through reduce_order()
        // auxiliaries deliberately does NOT count (rev154c): in mixed form the seam
        // copies of a chained row differ through one-sided stencils of the auxiliary,
        // and their redundancy would need the auxiliary's own derivative continuity,
        // which no plan imposes (MEASURED: PDE31b defR=3 = the 3 consumed HEAT copies;
        // PDE30b/38/PDE0 same class).
        std::map<size_t,int> order_in;
        for( auto const& op : sg.l_op ){
          if( !op || !op->sameid( typeid(FFPartial) ) ) continue;
          auto const* pop = mc::type_cast<FFPartial const>( op );
          if( !pop ) continue;
          for( auto const& [dv,ord] : pop->Indep().expr ){
            if( !ord ) continue;
            size_t const did = static_cast<size_t>( dv.id().second );
            order_in[ did ] = std::max( order_in[ did ], (int)ord );
          }
        }
        for( auto const& kv : _mEqn[ row_id ].dom ){
          size_t const did = static_cast<size_t>( kv.first.id().second );
          // A seam copy of this equation may be consumed for a d-claim only if the
          // equation is SECOND order in d (a copy is then redundant at the seam) and its
          // collocation set is symmetric (ALL or both ends excluded).  First order in d:
          // continuity must ADD a condition (TAU).  Measured: GAS_V/GAS_C (MBC1),
          // MMPDE27's PDE in t.
          bool const sym = ( kv.second == FFDom::ALL || kv.second == ( FFDom::ALL - FFDom::LB - FFDom::UB ) );
          // REFUTED (2026-09-05, PDE20f cell C): an order-0 row is NOT a seam duplicate
          // when it references other states that are themselves seam-discontinuous
          // (w - Dz_u: the two copies bind different Dz_u dofs).  Consumption stays
          // gated on direct second order with a symmetric set; order 0 and 1 are TAU.
          in.row_full_in_dir[ { row_id, did } ] = sym && ( order_in[ did ] >= 2 );
        }
      }
    }
    // Drop eligibility per block and W-protected claims (STRONG scan policy inputs).
    { std::set<int> blks;
      for( auto const& kv : claims ) blks.insert( kv.first.block_id );
      for( int b : blks ) in.block_drop_eligible[ b ] = _block_drop_eligible( b );
    }
    // True tau-column coefficient per receiver edge (mirrors the weak_sat_terms scalar:
    // sigma0/1 * face factor [/ iface width for C1] * coupling * orientation; the global
    // INTERFACE.TRACE_TAU_SCALE cancels in any dependency test and is omitted).
    in.edge_trace_coeff.reserve( draft.receiver_edges.size() );
    for( auto const& edge : draft.receiver_edges ){
      FFDom const* face_dom = nullptr;
      for( auto const& kd : _mDom )
        if( static_cast<size_t>( kd.first.id().second ) == edge.claim_key.dom_id ){ face_dom = &kd.second; break; }
      double tf = 1.;
      if( face_dom ){
        tf = _eff_sat_sigma0() * _sat_tau_face( *face_dom, edge.claim_key.iel_lo );
        if( edge.claim_key.kind == CLAIM_EXACT_C1 ){
          double const iface_w = _sat_iface_width( *face_dom, edge.claim_key.iel_lo );
          tf = ( iface_w > 0. ) ? ( _eff_sat_sigma1() * _sat_tau_face( *face_dom, edge.claim_key.iel_lo ) / iface_w ) : 0.;
        }
      }
      in.edge_trace_coeff.push_back( tf * edge.coupling * edge.orientation );
    }
    // Receiver row copies tagged with the DOF on their side.
    in.rows_at_dof.reserve( draft.receiver_edges.size() );
    for( size_t i = 0; i < draft.receiver_edges.size(); ++i ){
      auto const& edge = draft.receiver_edges[i];
      auto ite = edge_index_of_claim.find( edge.claim_key );
      if( ite == edge_index_of_claim.end() ) continue;
      PlanInput::DofEdge const& de = in.dof_edges[ ite->second ];
      PlanInput::RowAtDof r;
      r.row = FrozenPhysicalResidualRowKey{ edge.row_id, edge.receiver_block_el, edge.receiver_emit_flat };
      r.side_hi = ( edge.receiver_iel != edge.claim_key.iel_lo );
      // rev155a: non-nodal claims' rows are packed with dof=SIZE_MAX; the builder resolves
      // the side DOF under trace_dofs and skips them otherwise (pre-155a behaviour).
      r.dof = !de.nodal ? std::numeric_limits<size_t>::max()
            : ( r.side_hi ? de.dof_hi : de.dof_lo );
      r.edge_index = i;
      PlanRowRole const role = ( edge.row_id < in.row_role.size() ) ? in.row_role[ edge.row_id ] : PlanRowRole::AUTO;
      r.replaceable = ( role == PlanRowRole::INTERIOR ) && !edge.hard_physical_row && !edge.tangential_physical_boundary;
      r.tau_allowed = edge.trace_tau_allowed;
      r.sat_allowed = edge.weak_sat_allowed;
      r.coeff       = edge.coupling * edge.orientation;
      in.rows_at_dof.push_back( r );
    }
  }
  return in;
}

inline void
OCFEBASE::_store_frozen_interface_plan
( std::map<t_InterfaceClaimKey,size_t> const& claims,
  InterfacePlanDraft& draft,
  PlanInput const& in ) const
{
  // rev152d: pack -> build (legacy decision layer, ocplan_build_legacy.hpp) -> apply ->
  // materialise.  No decision is taken in this function.
  (void)claims;
  PlanDecisions dec; PlanReport rep;
  _planIsGlobal = false;
  // rev264 (C3 step 2, DELETION): the legacy builder is GONE.  Evidence, in order:
  //   rev262 sweep  legacy SKIPPED on all 822 plan decisions; 0 verdict differences over 151 programs;
  //                 IC_WEAK n/m/rank and k' identical.  GLOBAL decided every plan.
  //   rev263 sweep  the one path that still called legacy -- the deferred re-setup when the redundancy-drop
  //                 cycle runs under a GLOBAL plan -- was DISABLED corpus-wide: 11 drop-cycle events across
  //                 5 programs (hybrid_adequacy, STRONGCORNER1, PDE20f, PDE31b, closure), 0 deferrals, and
  //                 identical drop sets AND error figures on every one (STRONGCORNER1 6, PDE31b 6,
  //                 closure 8, PDE20f 20).  The rev154f defect it guarded (40 vs 20 drops, e_s2 16x) no
  //                 longer reproduces: the spanning-forest realisation, trace DOFs and the drop-policy work
  //                 fixed it elsewhere.
  // So OCPlanLegacyBuilder::build is not called, _globalDeferSetup is never set, and
  // ocplan_build_legacy.hpp (954 lines) is deletable.

  // rev153c EXPERIMENT (CRONOS_PLAN_APPLY_GLOBAL=1, default off): apply the GLOBAL
  // builder's decisions instead of the legacy ones (IC_TRACE/IC_STRONG only; the legacy
  // corner rule stays active upstream in this experiment).  Feel for the solve before
  // the rev154 switch.
  // rev156: the GLOBAL spanning-forest plan with trace DOFs is the DEFAULT (corpus-accepted
  // 2026-09-06, sweep 173913: 59 field-wise improvements vs rev150d, no regressions).  The
  // knobs remain as opt-outs (=0) for attribution.
  // 2026-10-07 (WORKPLAN 3.B batch 1): the opt-outs CRONOS_PLAN_APPLY_GLOBAL, _STRONG and _WEAK are retired.  Since
  // rev264 deleted the legacy builder, an opt-out applied the EMPTY default plan -- no valid alternative was left.
  // CRONOS_PLAN_APPLY_TRACE (which printed those knobs) goes with them.
  { bool const kApply = true, kApplyStrong = true;
    // rev259 (Benoit, 2026-09-15): WEAK was never covered, so under IC_WEAK the plan still comes from the
    // LEGACY builder and GLOBAL runs only as a shadow.  That is why a change to the GLOBAL builder's weak
    // branch did nothing (nVar unchanged at 1250, realisation mix identical) -- and patching the legacy
    // builder instead would mean changing code the workplan intends to retire.  So: extend the apply to WEAK
    // behind its own knob, verify the two builders AGREE there (the shadow already reports "vs legacy"), and
    // only then land the weak-TAU change in the builder that survives.
    bool const kApplyWeak = true;   // rev262: GLOBAL decides in all modes
    if( kApply
     && ( in.options.imposition == PlanOptions::Imposition::TRACE
       || ( kApplyStrong && in.options.imposition == PlanOptions::Imposition::STRONG )
       || ( kApplyWeak   && in.options.imposition == PlanOptions::Imposition::WEAK ) ) ){
      PlanDecisions gdec; PlanReport grep;
      // rev155d EXPERIMENT: CRONOS_PLAN_TRACE_DOFS_APPLY=1 applies the trace-DOF plan (item 2,
      // first applied step).  Gated per model by the determinacy line, never on by default.
      // (CRONOS_PLAN_TRACE_DOFS_APPLY retired 2026-10-07: the trace-DOF plan always applies)
      { PlanInput in_t = in; in_t.options.trace_dofs = true; OCPlanGlobalBuilder::build( in_t, gdec, grep ); }
      _structDeficit = static_cast<long>( grep.shadow.n_rank_deficit );
      _structDeficientClusters = static_cast<long>( grep.shadow.n_clusters_rank_deficient );
      grep.shadow.legacy_exact = dec.exact_replacement_rows.size();
      grep.shadow.legacy_tau   = dec.n_trace_var;
      OCPlanReport::print_shadow( _disp( 3 ), grep );
      if( char const* ds = std::getenv( "CRONOS_PLAN_DUMP_STATE" ) ){
        for( auto const& kv : in.state_name ) if( kv.second == ds ) OCPlanReport::dump_state( std::cout, grep, in, kv.first );
      }
      _disp( 1 ) << "OCFESLV::setup ** PLAN: GLOBAL decisions APPLIED: exact="
                << gdec.exact_replacement_rows.size() << " tau=" << gdec.n_trace_var
                << " constraints=" << gdec.n_trace_constraint << std::endl;
      // A weak plan WITHOUT multipliers can still number tau ids (standalone claims) that no column or row
      // realises; they reach _nTraceVar and switch the solver path.  Reported, not changed (NOTES_20260917a).
      if( in.options.imposition == PlanOptions::Imposition::WEAK && !in.options.weak_multipliers && gdec.n_trace_var )
        _disp( 2 ) << "OCFESLV::setup ** PLAN: NOTICE IC_WEAK plan numbers " << gdec.n_trace_var
                  << " tau id(s) that no column or row realises (standalone claims)" << std::endl;
      dec = std::move( gdec );
      _keepExplicitReason.clear();
      _planIsGlobal = true;
    }
  }

  // (CRONOS_PLAN_SHADOW retired 2026-10-07, WORKPLAN 3.D: it compared the GLOBAL plan with the LEGACY builder's
  // decisions, which have been the empty default since rev264 deleted that builder; the GLOBAL report it also printed
  // is printed on the main path at the same display level.)

  // (CRONOS_PLAN_TRACE_DOFS retired 2026-10-07, WORKPLAN 3.D: it compared a build WITHOUT trace DOFs against one with
  // them, "against the default" -- since batch 1 the trace-DOF build IS the applied plan.)  Its two censuses stay, now
  // standalone (they used to run only when CRONOS_PLAN_TRACE_DOFS=1 was set too), on the applied plan:
  // CRONOS_PLAN_RECV_CENSUS=1 and CRONOS_PLAN_DUMP_CLAIMS=<state>.
  { bool const kRecvCensus = std::getenv( "CRONOS_PLAN_RECV_CENSUS" ) != nullptr;
    bool const kDumpClaims = std::getenv( "CRONOS_PLAN_DUMP_CLAIMS" ) != nullptr;
    if( ( kRecvCensus || kDumpClaims ) && in.options.imposition != PlanOptions::Imposition::WEAK ){
      PlanInput in2 = in; in2.options.trace_dofs = true;
      PlanDecisions d1; PlanReport r1; OCPlanGlobalBuilder::build( in2, d1, r1 );
      // rev155d: receiver census by row ROLE per (state, dom) of the kept claims, and, for
      // auxiliaries, how many receivers are the auxiliary's OWN LINK row.  Instrument for the
      // tau-free-row criterion (Sec.3.5): an implication through a LINK row is exact only if that
      // row receives no tau.  CRONOS_PLAN_RECV_CENSUS=1.
      if( std::getenv( "CRONOS_PLAN_RECV_CENSUS" ) ){
        auto role_name = []( EqnRole r ){ switch( r ){ case EqnRole::INTERIOR: return "INTERIOR"; case EqnRole::LINK: return "LINK";
          case EqnRole::INTERFACE: return "INTERFACE"; case EqnRole::BOUNDARY: return "BOUNDARY"; case EqnRole::INITIAL: return "INITIAL"; default: return "OTHER"; } };
        std::map<size_t,FFVar const*> link_of_aux;   // aux state id -> its LINK equation var (expr)
        for( auto const& a : _auxDef ) link_of_aux[ static_cast<size_t>( a.aux.id().second ) ] = &a.expr;
        std::map<std::pair<size_t,size_t>, std::map<std::string,size_t>> cen;
        std::map<std::pair<size_t,size_t>, size_t> ownlink, nclaims;
        for( auto const& kv : in.claims ) ++nclaims[ { kv.first.state_id, kv.first.dom_id } ];
        for( auto const& re : *in.receiver_edges ){
          if( !in.claims.count( re.claim_key ) || re.row_id >= _mEqn.size() ) continue;
          auto key = std::make_pair( re.claim_key.state_id, re.claim_key.dom_id );
          ++cen[ key ][ role_name( _mEqn[ re.row_id ].opt->role ) ];
          auto itl = link_of_aux.find( re.claim_key.state_id );
          if( itl != link_of_aux.end() && _mEqn[ re.row_id ].var.id() == itl->second->id() ) ++ownlink[ key ];
        }
        for( auto const& [key, roles] : cen ){
          auto itn = in.state_name.find( key.first );
          _disp( 3 ) << "  [plan/recv] state=" << ( itn != in.state_name.end() ? itn->second : std::to_string( key.first ) )
                    << " dom=" << key.second << " claims=" << nclaims[ key ] << " receivers:";
          for( auto const& [r, n] : roles ) _disp( 3 ) << " " << r << "=" << n;
          if( link_of_aux.count( key.first ) ) _disp( 3 ) << " | own-LINK=" << ownlink[ key ];
          _disp( 3 ) << std::endl;
        }
      }
      // rev155a: claim-key census for one state (which seam nodes carry a claim).
      if( char const* dc = std::getenv( "CRONOS_PLAN_DUMP_CLAIMS" ) ){
        for( auto const& kv : in.claims ){
          auto itn = in.state_name.find( kv.first.state_id );
          if( itn == in.state_name.end() || itn->second != dc ) continue;
          _disp( 3 ) << "  [plan/claims] " << dc << " dom=" << kv.first.dom_id << " iel_lo=" << kv.first.iel_lo
                    << " face=" << kv.first.face << " kind=" << kv.first.kind << " occ=" << kv.second
                    << ( in.dropped_claims.count( OCPlan::freeze_claim_key( kv.first ) ) ? " DROPPED" : "" )
                    << ( in.unrealised_claims.count( kv.first ) ? " UNREALISED" : "" )
                    << ( !d1.claims_kept.count( kv.first ) ? " FOREST-REDUNDANT" : "" ) << std::endl;
        }
        for( auto const& e : in.corner_suppressed_edges ){
          auto itn = in.state_name.find( e.key.state_id );
          if( itn == in.state_name.end() || itn->second != dc ) continue;
          std::set<t_InterfaceClaimKey> const fdset( r1.shadow.forest_dropped.begin(), r1.shadow.forest_dropped.end() );
          bool const fd = fdset.count( e.key ) > 0;
          _disp( 3 ) << "  [plan/claims] " << dc << " dom=" << e.key.dom_id << " iel_lo=" << e.key.iel_lo
                    << " face=" << e.key.face << " kind=" << e.key.kind << " CORNER-SUPPRESSED" << ( fd ? " FOREST-DROPPED" : " KEPT-BY-FOREST" ) << std::endl;
        }
      }
    }
  }

  draft.claims.clear();
  for( auto const& kv : dec.claims_kept )
    draft.claims[ _freeze_interface_claim_key( kv.first ) ] = kv.second;
  draft.n_trace_constraint = dec.n_trace_constraint;

  _apply_plan_decisions( dec, draft );   // rev152c: decisions -> draft, once

  if( in.options.imposition == PlanOptions::Imposition::TRACE
   || in.options.imposition == PlanOptions::Imposition::STRONG ){
    // Materialise from the claims the builder KEPT (post two-pass DROP), not from the
    // full input map: legacy erased dropped claims before this call.
    if( !_build_frozen_trace_constraint_terms( dec.claims_kept, draft ) )
      throw Exceptions( Exceptions::INTERNAL );
  }
  else if( in.options.weak_multipliers && dec.n_trace_var ){
    // IC_WEAK WITH MULTIPLIERS: the builder kept exactly the claims it realised TAU, sized n_coll_var and
    // n_coll_eqn for them and numbered their taus as var[] columns.  One exact trace row per kept claim, as
    // IC_TRACE.  POSITIVE SIGNAL that the experiment ran: the line below, nvar and neqn both rising.
    if( !_build_frozen_trace_constraint_terms( dec.claims_kept, draft ) )
      throw Exceptions( Exceptions::INTERNAL );
    draft.weak_multipliers  = true;
    draft.weak_natural_only = in.options.weak_natural_only;
    _disp( 2 ) << "OCFESLV::setup ** PLAN: IC_WEAK+TAU -- " << dec.n_trace_var << " claim(s) imposed EXACTLY:"
              << " nvar " << in.n_coll_var << " -> " << dec.n_coll_var
              << ", neqn " << in.n_coll_eqn << " -> " << dec.n_coll_eqn
              << ", trace rows " << draft.trace_constraints.size() << std::endl;
  }
  else{
    draft.trace_constraints.clear();
    draft.trace_constraints_ready = true;
  }

  draft.frozen = true;
  draft.state  = OCPlan::State::FROZEN;
}

// rev152c: the single point where plan decisions enter the draft.
inline void
OCFEBASE::_apply_plan_decisions
( PlanDecisions const& dec, InterfacePlanDraft& draft ) const
{
  draft.tau_claims             = dec.tau_claims;
  draft.explicit_tau_claims    = dec.explicit_tau_claims;
  draft.exact_replacement_rows = dec.exact_replacement_rows;
  draft.reserved_definition_rows = dec.reserved_definition_rows;
  draft.n_trace_var            = dec.n_trace_var;
  draft.n_trace_var_explicit   = dec.n_trace_var_explicit;
  draft.n_coll_eqn             = dec.n_coll_eqn;
  draft.n_coll_var             = dec.n_coll_var;
  for( auto const& kv : dec.keep_explicit_reason ) _keepExplicitReason[ kv.first ] = kv.second;
  for( auto const& [i, allowed] : dec.trace_tau_allowed_updates )
    if( i < draft.receiver_edges.size() ) draft.receiver_edges[i].trace_tau_allowed = allowed;
}

inline bool
OCFEBASE::_trace_claim_coeffs
( t_InterfaceClaimKey const& claim,
  std::vector<std::pair<size_t,double>>& coeffs ) const
{
  coeffs.clear();

  auto state_by_id = [&]( size_t sid ) -> FFVar const*
  {
    for( auto const& kv : _mVar )
      if( static_cast<size_t>( kv.first.id().second ) == sid ) return &kv.first;
    return nullptr;
  };

  auto dom_by_id = [&]( size_t did ) -> FFVar const*
  {
    for( auto const& kv : _mDom )
      if( static_cast<size_t>( kv.first.id().second ) == did ) return &kv.first;
    return nullptr;
  };

  FFVar const* pstate = state_by_id( claim.state_id );
  FFVar const* pvar   = dom_by_id  ( claim.dom_id );
  if( !pstate || !pvar ) return false;

  auto itdep = _mVar.find( *pstate );
  auto itdom = _mDom.find( *pvar );
  if( itdep == _mVar.end() || itdom == _mDom.end() ) return false;

  FFDom const& face_dom = itdom->second;
  if( claim.iel_lo + 1 >= face_dom.n_elem ) return false;

  std::map<FFVar const*, FFDom const*, lt_FFVar> mpDom;
  for( auto const& dv : itdep->second ) mpDom[&dv] = &_mDom.at( dv );
  if( mpDom.find( pvar ) == mpDom.end() ) return false;

  std::map<FFVar, size_t, lt_FFVar> ndx_el1, ndx_el2;
  ndx_el1[*pvar] = claim.iel_lo;
  ndx_el2[*pvar] = claim.iel_lo + 1;

  size_t rem = claim.face;
  size_t f_local = 0;
  size_t f_mult  = 1;
  for( auto const& kv : mpDom ){
    FFVar const* opvar = kv.first;
    FFDom const* opdom = kv.second;
    if( opvar->id() == pvar->id() ) continue;
    size_t const radix = opdom->n_elem * opdom->n_node;
    size_t const enc = rem % radix;
    rem /= radix;
    size_t const ed = enc / opdom->n_node;
    size_t const kd = enc % opdom->n_node;
    if( ed >= opdom->n_elem || kd >= opdom->n_node ) return false;
    ndx_el1[*opvar] = ed;
    ndx_el2[*opvar] = ed;
    f_local += kd * f_mult;
    f_mult  *= opdom->n_node;
  }
  if( rem != 0 ) return false;

  auto const so = this->stride( mpDom, pvar, true );
  size_t const stride = so.first, outer = so.second;
  if( !stride || !outer || f_local >= stride * outer ) return false;

  size_t const inc1 = pos_state( *pstate, ndx_el1 );
  size_t const inc2 = pos_state( *pstate, ndx_el2 );
  size_t const sface = f_local % stride;
  size_t const oface = f_local / stride;

  int const ord = ( claim.kind == CLAIM_EXACT_C1 ) ? 1 : 0;
  std::vector<double> wl, wr;
  if( !face_dom.w_eval_lagrange( wl,  1.0, ord )
   || !face_dom.w_eval_lagrange( wr, -1.0, ord ) )
    return false;
  if( wl.size() != face_dom.n_node || wr.size() != face_dom.n_node ) return false;

  double scale_l = 1.0, scale_r = 1.0;
  if( ord == 1 ){
    scale_l = 2.0 / face_dom.elem_width( claim.iel_lo );
    scale_r = 2.0 / face_dom.elem_width( claim.iel_lo + 1 );
  }

  coeffs.reserve( 2 * face_dom.n_node );
  FFVar const* const waux = _lifted_claim_aux( claim );    // a transition at this face: x(tau^-) -> w
  double sum_l = 0.;
  for( size_t k = 0; k < face_dom.n_node; ++k ){
    size_t const il = inc1 + sface + k * stride + oface * face_dom.n_node * stride;
    size_t const ir = inc2 + sface + k * stride + oface * face_dom.n_node * stride;
    double const cl =  wl[k] * scale_l;
    double const cr = -wr[k] * scale_r;
    if( waux ) sum_l += cl;
    else if( std::fabs( cl ) > 64. * DBL_EPSILON ) coeffs.emplace_back( il, cl );
    if( std::fabs( cr ) > 64. * DBL_EPSILON ) coeffs.emplace_back( ir, cr );
  }
  if( waux ) coeffs.emplace_back( _lifted_aux_col( *waux, ndx_el1, *pvar, sface, oface, stride ), sum_l );
  return !coeffs.empty();
}

inline bool
OCFEBASE::_trace_claim_sides
( t_InterfaceClaimKey const& claim,
  std::vector<std::pair<size_t,double>>& lo,
  std::vector<std::pair<size_t,double>>& hi ) const
{
  lo.clear(); hi.clear();
  FFVar const* pstate = nullptr; FFVar const* pvar = nullptr;
  for( auto const& kv : _mVar ) if( static_cast<size_t>( kv.first.id().second ) == claim.state_id ){ pstate = &kv.first; break; }
  for( auto const& kv : _mDom ) if( static_cast<size_t>( kv.first.id().second ) == claim.dom_id ){ pvar = &kv.first; break; }
  if( !pstate || !pvar ) return false;
  auto itdep = _mVar.find( *pstate );
  auto itdom = _mDom.find( *pvar );
  if( itdep == _mVar.end() || itdom == _mDom.end() ) return false;
  FFDom const& face_dom = itdom->second;
  if( claim.iel_lo + 1 >= face_dom.n_elem ) return false;
  std::map<FFVar const*, FFDom const*, lt_FFVar> mpDom;
  for( auto const& dv : itdep->second ) mpDom[&dv] = &_mDom.at( dv );
  if( mpDom.find( pvar ) == mpDom.end() ) return false;
  std::map<FFVar, size_t, lt_FFVar> ndx_el1, ndx_el2;
  ndx_el1[*pvar] = claim.iel_lo;
  ndx_el2[*pvar] = claim.iel_lo + 1;
  size_t rem = claim.face, f_local = 0, f_mult = 1;
  for( auto const& kv : mpDom ){
    FFVar const* opvar = kv.first; FFDom const* opdom = kv.second;
    if( opvar->id() == pvar->id() ) continue;
    size_t const radix = opdom->n_elem * opdom->n_node;
    size_t const enc = rem % radix; rem /= radix;
    size_t const ed = enc / opdom->n_node, kd = enc % opdom->n_node;
    if( ed >= opdom->n_elem || kd >= opdom->n_node ) return false;
    ndx_el1[*opvar] = ed; ndx_el2[*opvar] = ed;
    f_local += kd * f_mult; f_mult *= opdom->n_node;
  }
  if( rem != 0 ) return false;
  auto const so = this->stride( mpDom, pvar, true );
  size_t const stride = so.first, outer = so.second;
  if( !stride || !outer || f_local >= stride * outer ) return false;
  size_t const inc1 = pos_state( *pstate, ndx_el1 ), inc2 = pos_state( *pstate, ndx_el2 );
  size_t const sface = f_local % stride, oface = f_local / stride;
  int const ord = ( claim.kind == CLAIM_EXACT_C1 ) ? 1 : 0;
  std::vector<double> wl, wr;
  if( !face_dom.w_eval_lagrange( wl,  1.0, ord ) || !face_dom.w_eval_lagrange( wr, -1.0, ord ) ) return false;
  if( wl.size() != face_dom.n_node || wr.size() != face_dom.n_node ) return false;
  double scale_l = 1.0, scale_r = 1.0;
  if( ord == 1 ){ scale_l = 2.0 / face_dom.elem_width( claim.iel_lo ); scale_r = 2.0 / face_dom.elem_width( claim.iel_lo + 1 ); }
  FFVar const* const waux = _lifted_claim_aux( claim );  // a transition at this face: the low side is w
  double sum_l = 0.;
  for( size_t k = 0; k < face_dom.n_node; ++k ){
    size_t const il = inc1 + sface + k * stride + oface * face_dom.n_node * stride;
    size_t const ir = inc2 + sface + k * stride + oface * face_dom.n_node * stride;
    double const cl = wl[k] * scale_l, cr = wr[k] * scale_r;
    if( waux ) sum_l += cl;
    else if( std::fabs( cl ) > 64. * DBL_EPSILON ) lo.emplace_back( il, cl );
    if( std::fabs( cr ) > 64. * DBL_EPSILON ) hi.emplace_back( ir, cr );
  }
  if( waux ) lo.emplace_back( _lifted_aux_col( *waux, ndx_el1, *pvar, sface, oface, stride ), sum_l );
  return !lo.empty() && !hi.empty();
}

inline void
OCFEBASE::_pinned_column_pairs
( std::vector<PlanInput::DofEdge>& out ) const
{
  // Bare state-id sets of the INITIAL/BOUNDARY equations, computed once.
  struct PinEq { t_Eqn const* er; std::set<size_t> bare; bool used = false; };
  std::vector<PinEq> pins;
  for( auto const& er : _mEqn ){
    if( er.opt->role != EqnRole::INITIAL && er.opt->role != EqnRole::BOUNDARY ) continue;
    FFVar expr = er.var;
    FFGraph* dg = static_cast<FFGraph*>( expr.dag() );
    if( !dg ) continue;
    FFSubgraph sg = dg->subgraph( 1, &expr );
    PinEq pe; pe.er = &er;
    for( auto const* op : sg.l_op ){
      if( !op || op->type != FFOp::VAR || !op->varout[0] ) continue;
      if( _mVar.find( *op->varout[0] ) == _mVar.end() ) continue;
      pe.bare.insert( static_cast<size_t>( op->varout[0]->id().second ) );
    }
    if( pe.bare.size() == 1 ) pins.push_back( pe );   // FIX1 class only (single-state pins)
  }
  if( pins.empty() ) return;
  auto in_mask = [&]( int m, bool isLB, bool isUB ) -> bool {
    if( m == FFDom::ALL ) return true;
    if( m == FFDom::LB  ) return isLB;
    if( m == FFDom::UB  ) return isUB;
    if( m == FFDom::ALL - FFDom::LB )             return !isLB;
    if( m == FFDom::ALL - FFDom::UB )             return !isUB;
    if( m == FFDom::ALL - FFDom::LB - FFDom::UB ) return !isLB && !isUB;
    return true;
  };
  std::vector<std::pair<size_t,double>> lo, hi;
  for( auto const& kv : _mVar ){
    FFVar const& st = kv.first;
    size_t const sid = static_cast<size_t>( st.id().second );
    std::map<FFVar const*, FFDom const*, lt_FFVar> mpDom;
    for( auto const& dv : kv.second ) mpDom[&dv] = &_mDom.at( dv );
    for( auto const& fkv : mpDom ){
      FFVar const* pvar = fkv.first; FFDom const& fdom = *fkv.second;
      if( fdom.n_elem < 2 ) continue;
      // transverse dims in mpDom order (the face radix order of _trace_claim_coeffs)
      std::vector<std::pair<FFVar const*, FFDom const*>> tr;
      for( auto const& okv : mpDom ) if( okv.first->id() != pvar->id() ) tr.push_back( { okv.first, okv.second } );
      size_t nface = 1; for( auto const& t : tr ) nface *= t.second->n_elem * t.second->n_node;
      for( size_t iel_lo = 0; iel_lo + 1 < fdom.n_elem; ++iel_lo ){
        for( size_t face = 0; face < nface; ++face ){
          // decode transverse (element, node) per dim
          size_t rem = face; bool pinned = false;
          std::vector<std::tuple<FFVar const*, bool, bool>> pos;   // (var, isLB, isUB)
          for( auto const& t : tr ){
            size_t const radix = t.second->n_elem * t.second->n_node;
            size_t const enc = rem % radix; rem /= radix;
            size_t const ed = enc / t.second->n_node, kd = enc % t.second->n_node;
            pos.push_back( { t.first, ed == 0 && kd == 0, ed + 1 == t.second->n_elem && kd + 1 == t.second->n_node } );
          }
          for( auto const& pe : pins ){
            if( !pe.bare.count( sid ) ) continue;
            // The two seam-side nodes are INTERIOR in the seam direction, so any mask covering an
            // interior node pins both (ALL, ALL-LB, ALL-UB, ALL-LB-UB); LB/UB masks do not.  FIX1's
            // "== ALL" is stricter than needed: PDE5's x=0/x=1 BCs are {T_INT, LB, Y_INT} and pin
            // both sides of every y-seam node (MEASURED: the 2 corner claims the forest kept).
            auto itm = pe.er->dom.find( *pvar );
            if( itm == pe.er->dom.end() || !in_mask( itm->second, false, false ) ) continue;
            bool here = true;
            for( auto const& [tv, isLB, isUB] : pos ){
              auto it2 = pe.er->dom.find( *tv );
              if( it2 == pe.er->dom.end() || !in_mask( it2->second, isLB, isUB ) ){ here = false; break; }
            }
            if( here ){ pinned = true; break; }
          }
          if( !pinned ) continue;
          t_InterfaceClaimKey key{};
          key.state_id = sid; key.dom_id = static_cast<size_t>( pvar->id().second );
          key.iel_lo = iel_lo; key.face = face; key.kind = CLAIM_EXACT_C0;
          if( !_trace_claim_sides( key, lo, hi ) ) continue;
          if( lo.size() != 1 || hi.size() != 1 ) continue;                 // nodal pairs only
          if( std::fabs( std::fabs( lo[0].second ) - 1. ) > 1e-12 || std::fabs( std::fabs( hi[0].second ) - 1. ) > 1e-12 ) continue;
          PlanInput::DofEdge e; e.key = key; e.nodal = true;
          e.dof_lo = lo[0].first; e.dof_hi = hi[0].first;
          e.lo_support = lo; e.hi_support = hi;
          out.push_back( e );
        }
      }
    }
  }
}

inline bool

//==========================================================================================================
// INTERFACE PLAN -- claims -> realisations -> rows and columns.  The DECISIONS live in ocplan.hpp; what is here
// is the assembly: which claims exist, which rows receive them, and how the chosen realisation becomes
// matrix structure.
// (banner placed from the definitions themselves: this region is lines 15300-16168 of the pre-rev268 file)
//==========================================================================================================

OCFEBASE::_sat_claim_coeffs
( t_InterfaceClaimKey const& claim,
  double const prefactor,
  std::vector<std::pair<size_t,double>>& coeffs ) const
{
  coeffs.clear();
  if( prefactor == 0. ) return true;  // valid but zero — nothing to add
  if( !_trace_claim_coeffs( claim, coeffs ) ) return false;
  for( auto& [col, c] : coeffs ) c *= prefactor;
  return !coeffs.empty();
}

inline bool
OCFEBASE::_build_frozen_trace_constraint_terms
( std::map<t_InterfaceClaimKey,size_t> const& claims,
  InterfacePlanDraft& draft ) const
{
  draft.trace_constraints.clear();
  size_t nrow = 0;
  for( auto const& kv : claims ) nrow += kv.second;
  draft.trace_constraints.reserve( nrow );

  size_t row = 0;
  for( auto const& kv : claims ){
    for( size_t i = 0; i < kv.second; ++i ){
      FrozenTraceConstraintTerm term;
      term.row  = row++;
      term.edge = kv.first;
      if( !_trace_claim_coeffs( kv.first, term.coeff ) ){
        std::cerr << "OCFESLV::_build_frozen_trace_constraint_terms ** failed "
                     "to build coefficients for exact trace claim: block="
                  << kv.first.block_id
                  << " state=" << kv.first.state_id
                  << " dom=" << kv.first.dom_id
                  << " iel=" << kv.first.iel_lo
                  << " face=" << kv.first.face
                  << " kind=" << kv.first.kind << std::endl;
        draft.trace_constraints.clear();
        draft.trace_constraints_ready = false;
        return false;
      }
      draft.trace_constraints.push_back( std::move( term ) );
    }
  }
  draft.trace_constraints_ready = true;
  return true;
}

inline void
OCFEBASE::_assert_frozen_interface_options_current() const
{
  OCPlan::FrozenOptions const& frozen = _plan.frozen_options();

  bool const changed =
       frozen.imposition_type != static_cast<int>( options.INTERFACE.IMPOSITION )
    || frozen.interface_type  != static_cast<int>( options.INTERFACE.TYPE )
    // rev51: compare against the EFFECTIVE values.  The plan is built from
    // _eff_sat_sigma0/1, so comparing against options.INTERFACE.SAT_SIGMA0/1 would record and check a
    // value the plan does not use -- and a driver that sweeps sigma0 internally would then
    // compare against the wrong number and either miss a rebuild or force a spurious one.
    || frozen.sat_sigma0      != _eff_sat_sigma0()
    || frozen.sat_sigma1      != _eff_sat_sigma1()
    || frozen.trace_tau_scale != options.INTERFACE.TRACE_TAU_SCALE;

  if( changed ){
    _disp( 2 ) << "OCFESLV::_eval_eqn ** interface-plan options changed after "
                 "setup(); SAT/TRACE coefficients are frozen in the "
                 "executable interface plan.  Rerun setup() after changing "
                 "INTERFACE.IMPOSITION, INTERFACE.TYPE, INTERFACE.SAT_SIGMA0 "
                 "or INTERFACE.TRACE_TAU_SCALE.  Frozen/current values: "
              << "imposition=" << frozen.imposition_type << "/"
              << static_cast<int>( options.INTERFACE.IMPOSITION )
              << " interface=" << frozen.interface_type << "/"
              << static_cast<int>( options.INTERFACE.TYPE )
              << " sigma0=" << frozen.sat_sigma0 << "/"
              << options.INTERFACE.SAT_SIGMA0
              << " sigma1=" << frozen.sat_sigma1 << "/"
              << _eff_sat_sigma1()
              << " trace_scale=" << frozen.trace_tau_scale << "/"
              << options.INTERFACE.TRACE_TAU_SCALE << std::endl;
    throw Exceptions( Exceptions::INTERNAL );
  }
}

inline bool
OCFEBASE::_build_interface_plan_draft
( InterfacePlanDraft& draft )
{
  // rev72 INSTRUMENT [evoface].  The rev71 placement -- inside the per-face causal block
  // -- printed NOTHING, not only on OCFE_MMPDE27 but also on OCFE_DAE3 and OCFE_PDE0,
  // the two models that block's own comment cites as where its treatment matters
  // ("mono==march, DAE3: 1e-14", "the 26x PDE0 STRONG-mono loss").  Zero output on the
  // cited models means the INSTRUMENT was wrong, not the path.  Placed here, at the top
  // of the draft builder, execution is guaranteed.
  //
  // WHAT THIS IS FOR.  Measured on MMPDE27 (--gdirect, nel_xi=4, NEL_T=4): with
  // set_evolution_domain(t) the plan emits 81 taus and NO claims on the t domain; with
  // auto-detection it emits 606, of which 48 are on dom0 (t), carrying k=12 cross-domain
  // kernel directions -- one per interior t-seam per xi-element.  So t-faces get ORDINARY
  // TWO-WAY continuity claims in the auto case and causal one-way treatment in the user
  // case, and the question is which predicate decides that.
  //
  // FIVE mechanisms have now been proposed from reading this header and refuted by
  // measurement.  _evolution_dom_set is NOT the discriminator: _infer_evolution_domain()
  // sets it true on the auto path too, and runs BEFORE this function (_classify_pde at
  // 10943, this at 11086).  The fields below are printed rather than reasoned about.
  // ---------------------------------------------------------------------------
  // rev81 RAW DUMP [depdump] -- print the DATA, not a verdict computed from it
  // ---------------------------------------------------------------------------
  // Four attempts at an automated "can this state inherit a claim" test have failed, each
  // one a verdict computed from dependence data that was never inspected:
  //   rev74  "bare rows define the state"      -> right for g, false positive on PSA's c1
  //   rev77  "own derivative defines it"       -> right for c1, LOST g (GGDEF consumes it)
  //   rev78  FFDep linearity                   -> wrong dep key; uniform "NOT linear"
  //   rev80  + partial-free requirement        -> now rejects EVERYTHING, including g/xg
  //
  // g and xg in DEF (g - M*xg) are manifestly linear in each with no partial on either, so
  // rev80's verdict is wrong on the merits and the next step is NOT a fifth patch.  This
  // dump prints, per equation: the dep map as it actually is (key, resolved name, type),
  // which states sit under an FFPartial, and the subgraph size.  Printed for EVERY equation
  // rather than only those involving at-risk states, so the calibration row is visible
  // whether or not the filters are right, and indexed to match [eqn-legend].
  {
    static bool const kDepDump = Options::_env_flag( "CRONOS_AUDIT_DEPDUMP", false );
    if( kDepDump ){
      _disp( 3 ) << "  [depdump] " << _mEqn.size() << " equation(s); dep types L/Q/P/R/N\n";
      for( size_t en = 0; en < _mEqn.size(); ++en ){
        FFSubgraph sgd = _dag->subgraph( 1, &_mEqn[en].var );
        std::set<std::string> under_p, bare_p;
        for( auto const* op : sgd.l_op ){
          if( !op ) continue;
          if( op->sameid( typeid(FFPartial) ) ){
            for( auto const* od : op->varin )
              if( od && _mVar.find( *od ) != _mVar.end() ) under_p.insert( od->name() );
            continue;
          }
          for( auto const* in : op->varin )
            if( in && _mVar.find( *in ) != _mVar.end() ) bare_p.insert( in->name() );
        }
        _disp( 3 ) << "  [depdump] e" << en << "  ops=" << sgd.l_op.size() << "  bare={";
        for( auto const& b : bare_p ) _disp( 3 ) << " " << b;
        _disp( 3 ) << " }  under-partial={";
        for( auto const& u : under_p ) _disp( 3 ) << " " << u;
        _disp( 3 ) << " }\n";
        if( en >= _depEqn.size() ){
          _disp( 3 ) << "  [depdump]      *** no _depEqn entry for this index ***\n";
          continue;
        }
        _disp( 3 ) << "  [depdump]      dep():";
        char const* tn[5] = { "L", "Q", "P", "R", "N" };
        for( auto const& [ky, ty] : _depEqn[en].dep() ){
          auto it = _mDep.find( ky );
          std::string nm = "<key " + std::to_string( ky ) + " unresolved>";
          if( it != _mDep.end() ){
            if( it->second.first == STATE && it->second.second < _vVar.size() )
              nm = _vVar[ it->second.second ].name();
            else nm = "<non-state kind " + std::to_string( it->second.first ) + ">";
          }
          int const ti = (int)ty;
          _disp( 3 ) << "  " << nm << "=" << ( ti >= 0 && ti < 5 ? tn[ti] : "?" );
        }
        _disp( 3 ) << "\n";

        // rev83: the same row through FFInv.  Five successive verdict-computing tests have
        // been wrong on this sub-problem, and the ONE time the raw data was printed it
        // explained all of them at once ([depdump] showed xg=N in DEF because the monitor
        // inlines the sinh oracle).  So print what FFInv RETURNS rather than what a test
        // concludes from it.
        //
        // States are the independent variables; domains, inputs and constants are left
        // default-constructed, i.e. CONSTANTS -- which is what should make a nonlinear
        // COEFFICIENT benign while a nonlinear relation between STATES is not.
        //
        // With the ocbase FFInv overloads in place, FFPartial/FFIntegral/FFEval return U,
        // so rows containing them are evaluable and report U for the affected variables
        // instead of throwing.  If a row still throws, the message says so rather than
        // silently reading as "not invertible" -- which is how the rev82 test failed.
        {
          char const* it_n[4] = { "L", "S", "N", "U" };
          _disp( 3 ) << "  [depdump]      inv():";
          try {
            // rev84: pass ALL FOUR variable groups, exactly as the FFDep call at
            // rev84:11038 does.  The rev83 version passed STATES ONLY, on the reasoning
            // that unpassed domains would "act as constants" -- they do not, they are an
            // ERROR, and the dump said so plainly:
            //     e1  inv(): Subgraph evaluation failed -- missing variable t
            //     e4  inv(): Subgraph evaluation failed -- missing variable xi
            // e1 is DEF (g - W*xg), whose monitor W = 1/x*_xi references t and xi, so the
            // calibration row was exactly the one that could not be evaluated.
            //
            // To make a group behave as CONSTANTS, pass it with DEFAULT-CONSTRUCTED FFInv
            // values: those carry an empty _inv map and so contribute no independent index.
            // Only the STATES are seeded via indep(), which is what makes a nonlinear
            // COEFFICIENT benign (it depends on domains, hence on nothing) while a
            // nonlinear relation between STATES is not.
            std::vector<FFInv> invVar( _vVar.size() );
            for( size_t k = 0; k < _vVar.size(); ++k ) invVar[k].indep( (int)k );
            std::vector<FFInv> invInp( _vInp.size() ), invCst( _vCst.size() ),
                               invDom( _vDom.size() );          // all default = constants
            FFInv invEqn;
            _dag->eval( sgd, 1, &_mEqn[en].var, &invEqn,
                        _vVar.size(), _vVar.data(), invVar.data(),
                        _vInp.size(), _vInp.data(), invInp.data(),
                        _vCst.size(), _vCst.data(), invCst.data(),
                        _vDom.size(), _vDom.data(), invDom.data() );
            if( invEqn.inv().empty() ) _disp( 3 ) << "  <empty>";
            for( auto const& [ind, ty] : invEqn.inv() ){
              int const ti = (int)ty;
              _disp( 3 ) << "  "
                        << ( ind >= 0 && (size_t)ind < _vVar.size()
                               ? _vVar[ind].name() : std::string("<idx?>") )
                        << "=" << ( ti >= 0 && ti < 4 ? it_n[ti] : "?" );
            }
          }
          catch( std::exception const& e ){ _disp( 3 ) << "  *** THREW: " << e.what() << " ***"; }
          catch( ... )                    { _disp( 3 ) << "  *** THREW (unknown) ***"; }
          _disp( 3 ) << "\n";
        }
      }
    }
  }

  // ---------------------------------------------------------------------------
  // rev73 SCAN [recvscan] -- Defect A, sized WITHOUT solving anything
  // ---------------------------------------------------------------------------
  // Defect A (handoff rev18/rev19): a continuity claim needs a RECEIVER ROW -- a row
  // carrying a derivative in the face direction, into which a flux/SAT term can be placed.
  // A state whose every defining and consuming row is DERIVATIVE-FREE in a face direction
  // gets no ENFORCING claim there, and if it carries a chain's integration constant the
  // exact modes admit a per-element freedom of dimension nel-1 that the weak system
  // excludes (IC_WEAK penalises every edge regardless of whether the row has a derivative).
  //
  // MEASURED on MESHONLY: DEF is g = M*xg with partials={} -- purely algebraic in xi -- and
  // both g and xg go unclaimed; the exact modes converge to a TANGLED mesh, min(xg) =
  // -9.65e+02.  Rewriting it as g = M*x_xi gives the row a xi-derivative, g gains 35
  // claims, the family collapses, and both exact modes reach IC_WEAK's answer exactly.  On
  // MMPDE26 the same rewrite closed rev14 section 3.1's nel_xi=32 threshold: conv=no at
  // nel_xi = 4, 8 AND 16 became conv=yes at all of them.
  //
  // WHY A SCAN BEFORE A FIX.  The framework remedy -- minting an enforcing claim for such a
  // state -- needs somewhere to put the flux term, and PDE2 carries 16 claims already
  // recorded as "continuity of a derivative auxiliary at a transverse boundary face already
  // pinned on both sides".  Over-determining PDE2 is the failure mode.  Sizing the problem
  // first tells us whether this is an MMPDE curiosity or corpus-wide, and it costs NOTHING:
  // "does any consuming row carry a face-direction derivative" is a STATIC property of the
  // equation legend, so every driver can be scanned at setup with no solve.
  //
  // Reports, per state and per domain, whether ANY equation mentioning that state carries a
  // partial of ANY state in that direction.  AT-RISK means no such row exists.
  {
    static bool const kRecvScan = Options::_env_flag( "CRONOS_AUDIT_RECVSCAN", false );
    if( kRecvScan ){
      // rev74 CRITERION FIX.  The rev73 scan asked "does ANY equation mentioning this state
      // carry a derivative in this direction".  That is too weak: it marked g as covered in
      // xi because GGDEF (gg = g_xi) mentions g AND carries 'xi -- but GGDEF CONSUMES g; it
      // is not g's defining row.  The scan therefore reported the SAME result for MESHONLY
      // base and --gdirect (2 at-risk, both dir=t) even though --gdirect is exactly the
      // change that fixes the model, and it missed the xi direction entirely.
      //
      // The distinction rev18 section 0 turns on is the DEFINING row: a state is defined by
      // the row where it appears BARE as the subject (DEF: g - M*xg, with g bare).  If that
      // row is derivative-free in a face direction it cannot receive a flux term there,
      // whatever other rows do with the state.  Being consumed by a derivative-carrying row
      // elsewhere supplies nothing.
      //
      // BARE vs UNDER-A-PARTIAL is read the way print_equation_legend reads it: partial ops
      // are recognised by op->sameid( typeid(FFPartial) ) and CONTINUE after recording, so
      // any state reaching the plain varin loop is a bare mention.
      //
      // A state may be bare in several rows; ALL of them are candidates to receive, so the
      // state is at risk only when EVERY row it appears bare in is derivative-free in the
      // direction.  That is the loose direction: it under-reports rather than over-reports
      // here, which is the safer bias for a criterion that must not miss the known case.
      // rev75 SECOND CONDITION -- "no receiver" alone is FAR too weak.
      //
      // MEASURED, rev74 criterion over the whole corpus: ~90 of 156 drivers flagged,
      // including the ENTIRE PSA and MBC families and PDE2 -- production models that
      // converge correctly, and PDE9 where exact imposition is 20-28x MORE accurate than
      // IC_WEAK.  So a missing receiver is NECESSARY but nowhere near SUFFICIENT, and any
      // blanket fix minting claims for every such state would touch most of the corpus.
      //
      // The discriminator is visible in that data.  Among the flagged states, only g in
      // MESHONLY/MMPDE27 is ALSO CONSUMED UNDER A DERIVATIVE IN THE SAME DIRECTION
      // (GGDEF: gg = g_xi).  That is the lifting-alias chain: g carries an integration
      // constant which g_xi propagates, so its continuity is LOAD-BEARING and a jump in it
      // changes the solution -- the nel-1 tangled family of handoff rev18 section 1.
      //
      // The states that flag harmlessly are POINTWISE in the flagged direction:
      //   PSA/PDE36-40  q1,q2 in z : bare in the mass balance, never differentiated in z
      //   PDE2          Drl_Sl in z: differentiated in rl, NOT in z
      //   ALGFIELD1/2   y,a,p      : no derivatives at all; all three modes agree exactly
      // A jump in a state nothing differentiates in that direction costs nothing.
      //
      // So: at risk iff (no defining row can receive in dir) AND (some row differentiates
      // it in dir).  The second condition is what makes the continuity load-bearing.
      size_t n_at_risk = 0, n_weak_only = 0;
      for( auto const& [sv, sdoms] : _mVar ){
        for( auto const& dv : sdoms ){
          bool defined_somewhere = false, some_definer_can_receive = false;
          bool consumed_under_partial_in_dir = false;
          for( auto const& eqn : _mEqn ){
            if( !eqn.opt->participate_in_classification ) continue;
            FFSubgraph sg = _dag->subgraph( 1, &eqn.var );
            bool bare_here = false, partial_in_dir_here = false;
            for( auto const* op : sg.l_op ){
              if( !op ) continue;
              if( op->sameid( typeid(FFPartial) ) ){
                auto const* pop = mc::type_cast<FFPartial const>( op );
                bool const is_self = [&]{ for( auto const* od : op->varin )
                                            if( od && od->id() == sv.id() ) return true;
                                          return false; }();
                if( pop )
                  for( auto const& [iv,ord] : pop->Indep().expr ){
                    (void)ord;
                    if( iv.id() != dv.id() ) continue;
                    partial_in_dir_here = true;
                    if( is_self ) consumed_under_partial_in_dir = true;
                  }
                continue;                     // recorded as a partial, NOT a bare mention
              }
              for( auto const* in : op->varin )
                if( in && in->id() == sv.id() ) bare_here = true;
            }
            // rev77 DEFINING-ROW FIX.  rev74-76 treated "appears BARE here" as "this row
            // defines the state".  That holds for an ALGEBRAIC state and FAILS for a
            // DIFFERENTIAL one, whose defining row contains its DERIVATIVE, not a bare
            // occurrence -- and whose bare occurrences are typically arguments of some
            // OTHER row's closure.
            //
            // MEASURED on OCFE_PSA1, which the rev75 criterion flagged as load-bearing:
            //   CONT1:  bare={ Dz_c1 }           partials={ Dz_c1'z, c1't, q1't }
            //   LDF1:   bare={ T, c1, c2, q1 }   partials={ q1't }
            // c1 appears bare ONLY in LDF1 -- which is dq1/dt = k(q1* - q1), the row that
            // defines q1, with c1 present merely as an argument of the isotherm
            // q1*(T,c1,c2).  c1's actual defining row is CONT1, its own mass balance, where
            // it appears under c1't -- and CONT1 CARRIES Dz_c1'z, so c1 HAS a receiver in z
            // and was never at risk.  The PSA flag was a FALSE POSITIVE, and with it the
            // rev23 section 3 conclusion that "Option B does not transfer to balance
            // couplings" -- that was reasoning about an artefact.
            //
            // Contrast g in MESHONLY: no g't anywhere, so g IS algebraic and its bare row
            // DEF (g - M*xg) genuinely defines it.  g appears under 'xi in GGDEF
            // (gg = g_xi), but that row's subject is gg -- consuming, not defining.
            //
            // So: a state under its OWN derivative in a row is DIFFERENTIAL and that row
            // defines it; otherwise it is ALGEBRAIC and its bare rows define it.
            // rev78: rev77 made "appears under its own derivative" imply "this row defines
            // it".  That fixed PSA's c1 (defined by CONT1, which carries c1't) and BROKE
            // MESHONLY's g: GGDEF (gg = g_xi) differentiates g in xi, so it counted as g's
            // definer and supplied a receiver -- but GGDEF's SUBJECT is gg; it CONSUMES g.
            // MESHONLY base then read 0 at-risk, losing the one separation that calibrates
            // this scan.  Reverted.
            //
            // The underlying problem is that "defining row" is NOT DERIVABLE for a
            // user-declared state: _auxDef records it (aux, expr, diff_dom, parent) but is
            // populated ONLY by _reduce_order, so g -- declared with add_state and defined
            // by a user-written DEF -- has no entry.  Structurally DEF is one equation among
            // six and nothing marks g as its subject.  Three heuristics were tried and two
            // failed on the calibration pair in OPPOSITE directions.
            //
            // The [inherit] probe below sidesteps subject attribution entirely.
            if( !bare_here ) continue;
            defined_somewhere = true;
            if( partial_in_dir_here ) some_definer_can_receive = true;
          }
          if( !defined_somewhere || some_definer_can_receive ) continue;
          if( !consumed_under_partial_in_dir ){
            ++n_weak_only;                    // pointwise: a jump here costs nothing
            continue;
          }
          ++n_at_risk;
          _disp( 3 ) << "  [recvscan] AT-RISK state=" << sv.name()
                    << " dir=" << dv.name()
                    << " -- no defining row can receive a continuity claim in this"
                       " direction, AND the state is differentiated in it elsewhere,"
                       " so its continuity is LOAD-BEARING\n";

          // rev76 UPSTREAM PROBE -- viability of Option B (handoff rev23).
          //
          // Two generic remedies were considered.  Option A, a model-level INLINE_ALIAS
          // transform substituting an alias by its definition (which is what --gdirect did
          // by hand: g = M*xg with xg = x_xi becomes g = M*x_xi, exact, and the row inherits
          // the derivative).  Option B, leave the model alone and PROPAGATE THE CLAIM: walk
          // the algebraic definition upstream to a state that CAN receive, and place the
          // claim there.
          //
          // Option A has a tension: REDUCE.ORDER exists to SPLIT high-order terms into
          // first-order alias chains, and that splitting is what creates derivative-free
          // LINK rows in the first place -- so INLINE_ALIAS would partly undo RED_MAIN and
          // the two would need a defined interaction.
          //
          // Option B is sound exactly when the algebraic link's coefficients are themselves
          // continuous: g's continuity is IMPLIED by xg's whenever M is continuous, which is
          // the header's own v2 note at auxiliary_algebraic_receiver_state -- and that note
          // records the implication holding under IC_WEAK while the exact modes drop the
          // coupling.  M = 1/x*_xi is a known smooth function here; a link whose coefficient
          // were a discontinuous STATE would not qualify.
          //
          // VIABILITY TEST, printed below: in the rows where this state appears bare, which
          // OTHER states appear bare, and does any of them have a receiver in this direction?
          // On MESHONLY the answer should be xg -- flagged NOT at risk, because GRAD
          // (xg = x_xi) carries 'xi and IS a receiver, so the claim has somewhere to go.
          // If a load-bearing state has NO such upstream anywhere in the corpus, Option B
          // does not cover it and the remedy differs by model class.
          for( auto const& eqn : _mEqn ){
            if( !eqn.opt->participate_in_classification ) continue;
            FFSubgraph sg2 = _dag->subgraph( 1, &eqn.var );
            bool subject_here = false;
            std::set<std::string> co_bare;
            std::vector<FFVar> co_vars;
            for( auto const* op : sg2.l_op ){
              if( !op || op->sameid( typeid(FFPartial) ) ) continue;
              for( auto const* in : op->varin ){
                if( !in || _mVar.find( *in ) == _mVar.end() ) continue;
                if( in->id() == sv.id() ){ subject_here = true; continue; }
                if( co_bare.insert( in->name() ).second ) co_vars.push_back( *in );
              }
            }
            if( !subject_here || co_vars.empty() ) continue;
            for( auto const& cv : co_vars ){
              // does THIS co-bare state have a defining row that receives in dir?
              bool up_receives = false;
              for( auto const& e2 : _mEqn ){
                if( !e2.opt->participate_in_classification ) continue;
                FFSubgraph sg3 = _dag->subgraph( 1, &e2.var );
                bool cv_bare = false, part_dir = false;
                for( auto const* op : sg3.l_op ){
                  if( !op ) continue;
                  if( op->sameid( typeid(FFPartial) ) ){
                    auto const* pop = mc::type_cast<FFPartial const>( op );
                    if( pop )
                      for( auto const& [iv,ord] : pop->Indep().expr ){
                        (void)ord; if( iv.id() == dv.id() ) part_dir = true; }
                    continue;
                  }
                  for( auto const* in : op->varin )
                    if( in && in->id() == cv.id() ) cv_bare = true;
                }
                if( cv_bare && part_dir ){ up_receives = true; break; }
              }
              // rev78 LINEARITY TEST [inherit] -- the soundness condition, from FFDep.
              //
              // A claim may be INHERITED from an upstream state only when the link is an
              // INVERTIBLE algebraic definition with continuous coefficients.  "Has a
              // receiver" is necessary and NOT sufficient: MESHONLY's g = M*xg qualifies
              // (M a known smooth function, so g's continuity follows from xg's), while
              // PSA's c1 and T merely CO-APPEAR in a balance -- neither determines the
              // other, and placing c1's claim on T would enforce continuity of the WRONG
              // QUANTITY.
              //
              // _depEqn is already computed for every equation (rev78:11028) and types each
              // dependence L/Q/P/R/N.  LINEAR IN BOTH is exactly invertibility of the link:
              //   DEF  : g - M*xg      g is L, xg is L        -> g = M*xg, invertible
              //   LDF1 : q1' - k(q1*(T,c1,c2) - q1)   c1 enters through the isotherm, N
              // So the test rests on a quantity the framework ALREADY COMPUTES, rather than
              // on any attribution of which row "defines" the state -- which is not
              // derivable for user-declared states (see the rev78 note above).
              // rev82 FFInv INVERTIBILITY TEST -- replaces the FFDep linearity test.
              //
              // WHY FFDep FAILED, measured by [depdump] on MESHONLY e1 (DEF: g - W*xg):
              //     ops=33  bare={ g, xg }  under-partial={ }
              //     dep():  <non-state>=N  <non-state>=N  xg=N  g=L
              // xg types N, so rev80's "not linear in both" was CORRECT ON THE DATA and the
              // claim that DEF is "manifestly linear" was wrong.  FFDep types the whole
              // dependence, and W = 1/x*_xi inlines the sinh oracle -- 33 ops of arcsinh/sinh
              // -- so anything multiplied by it comes out nonlinear.  A nonlinear but KNOWN
              // coefficient is perfectly fine for inheritance; FFDep cannot distinguish that
              // from a nonlinear relation with another STATE.
              //
              // FFInv can: its lattice is L (linear), S (separably linear), N (separably
              // nonlinear), U (undetermined).  S is exactly "appears once, multiplied by a
              // coefficient" -- true of g = M*xg AND of the multiply-form q*x_xi = u_xi, and
              // false of an isotherm q1*(T,c1,c2) where c1 sits inside a nonlinear function
              // alongside other states.
              //
              // SCOPE LIMIT, corrected 2026-10-15 (Benoit).  The earlier note here said FFInv is
              // NOT overloaded for FFPartial / FFIntegral / FFEval in ocbase.hpp.  That is OUT OF
              // DATE: MC++ 5.0.4's ocbase.hpp carries eval(... FFInv* ...) for all three (29 FFInv
              // references).  The limit is different, and its own include comment states it: "these
              // operations are not invertible, so their FFInv pass returns U" -- which is correct,
              // a derivative is not invertible in its operand.  So a whole-row FFInv pass over a
              // subgraph containing OpP still yields U, and the multiply-form q*x_xi = u_xi cannot
              // be typed S that way.  Typing it S needs a STATE-SCOPED query -- "is this row
              // separably linear in q, treating OpP(...) nodes as opaque coefficients?" -- which is
              // a different question from evaluating the row in FFInv arithmetic, and is what
              // EqnOptions::defines answers by declaration until such a query exists.
              // That costs nothing here: Defect A is BY DEFINITION about states whose defining
              // row is derivative-free, so the rows that matter are exactly the partial-free
              // ones.  Rows containing a partial are treated as "not an algebraic link"
              // WITHOUT evaluating, which sidesteps the missing overloads rather than
              // requiring them.
              char const* lin = "no shared algebraic row";
              {
                bool any_shared = false, any_invertible = false;
                for( size_t en = 0; en < _mEqn.size(); ++en ){
                  FFSubgraph sgq = _dag->subgraph( 1, &_mEqn[en].var );
                  bool has_partial = false, mentions_both = false;
                  bool sees_s = false, sees_u = false;
                  // rev85 BREAK REMOVED.  The `break` on the first FFPartial was left over
                  // from rev82, when rows containing partials were skipped entirely.  With
                  // rev84 evaluating them (ocbase returns U), the break instead ABORTS THE
                  // SCAN of the subgraph, so any state appearing after the first partial is
                  // never seen and mentions_both stays false.
                  //
                  // That silently produced PSA1's "no shared algebraic row" for c1/T: LDF1
                  // is q1' - k(q1*(T,c1,c2) - q1), whose FIRST operation is the partial
                  // q1't, so the loop exited before reaching c1 or T.  The verdict matched
                  // the expectation -- PSA should not qualify -- but for the wrong reason,
                  // and a result that is right by accident is worse than one that is wrong,
                  // because nothing flags it.
                  //
                  // has_partial is still recorded (it is informative) but no longer
                  // terminates the walk.
                  for( auto const* op : sgq.l_op ){
                    if( !op ) continue;
                    if( op->sameid( typeid(FFPartial) ) ){
                      has_partial = true;
                      // states UNDER the partial still count as mentioned: FFInv will type
                      // them U and the L/S test will reject them, which is the correct
                      // outcome reached by the right route.
                      for( auto const* od : op->varin ){
                        if( !od ) continue;
                        if( od->id() == sv.id() ) sees_s = true;
                        if( od->id() == cv.id() ) sees_u = true;
                      }
                      continue;
                    }
                    for( auto const* in : op->varin ){
                      if( !in ) continue;
                      if( in->id() == sv.id() ) sees_s = true;
                      if( in->id() == cv.id() ) sees_u = true;
                    }
                  }
                  mentions_both = sees_s && sees_u;
                  // rev84: rows containing a partial are NO LONGER SKIPPED.  With the
                  // ocbase FFInv overloads in place they evaluate, returning U for the
                  // differentiated variable and a real type for the others -- MEASURED:
                  //   e0 (GRAD: xg - x_xi)   inv():  x=U  xg=L
                  //   e2 (GGDEF: gg - g_xi)  inv():  g=U  gg=L
                  //   e3 (RELAX: tau x_t-gg) inv():  x=U  gg=L
                  // So a differential row can still certify a SEPARABLE state, and U (which
                  // the L/S test rejects) correctly blocks the differentiated one.  This is
                  // what makes the multiply-form q*x_xi = u_xi testable at all.
                  (void)has_partial;
                  if( !mentions_both ) continue;
                  any_shared = true;

                  // Evaluate this partial-free subgraph in FFInv arithmetic: states are the
                  // independent variables; domains, inputs and constants are left as
                  // default-constructed (i.e. constants), which is what makes a nonlinear
                  // COEFFICIENT benign while a nonlinear relation between STATES is not.
                  try {
                    // rev84: all four variable groups -- see the [depdump] note.  States
                    // seeded via indep(); inputs, constants and domains passed
                    // default-constructed so they act as CONSTANTS rather than being absent.
                    std::vector<FFInv> invVar( _vVar.size() );
                    for( size_t k = 0; k < _vVar.size(); ++k ) invVar[k].indep( (int)k );
                    std::vector<FFInv> invInp( _vInp.size() ), invCst( _vCst.size() ),
                                       invDom( _vDom.size() );
                    FFInv invEqn;
                    _dag->eval( sgq, 1, &_mEqn[en].var, &invEqn,
                                _vVar.size(), _vVar.data(), invVar.data(),
                                _vInp.size(), _vInp.data(), invInp.data(),
                                _vCst.size(), _vCst.data(), invCst.data(),
                                _vDom.size(), _vDom.data(), invDom.data() );
                    size_t is = _vVar.size(), iu = _vVar.size();
                    for( size_t k = 0; k < _vVar.size(); ++k ){
                      if( _vVar[k].id() == sv.id() ) is = k;
                      if( _vVar[k].id() == cv.id() ) iu = k;
                    }
                    if( is < _vVar.size() && iu < _vVar.size() ){
                      auto const ri = invEqn.inv( (int)is );
                      auto const ru = invEqn.inv( (int)iu );
                      bool const ok = ri.first && ru.first
                                   && ( ri.second == FFInv::TYPE::L || ri.second == FFInv::TYPE::S )
                                   && ( ru.second == FFInv::TYPE::L || ru.second == FFInv::TYPE::S );
                      if( ok ) any_invertible = true;
                    }
                  }
                  catch( ... ) { /* FFInv unavailable on this row: leave as not-invertible */ }
                }
                if( any_shared )
                  lin = any_invertible
                          ? "INVERTIBLE (FFInv L/S) in both -> algebraic link"
                          : "shared algebraic row(s), but NOT invertible in both";
              }
              _disp( 3 ) << "  [recvscan]     upstream " << cv.name()
                        << ( up_receives ? "  has a receiver" : "  has NO receiver" )
                        << ";  " << lin
                        << ( up_receives && lin[0] == 'L'
                               ? "   -> INHERIT viable: place the claim here" : "" )
                        << "\n";
            }
          }
        }
      }
      _disp( 3 ) << "  [recvscan] LOAD-BEARING at-risk = " << n_at_risk
                << "   (pointwise, harmless = " << n_weak_only << ")"
                << ( n_at_risk ? "   [Defect A candidate]"
                               : "   [clean: no load-bearing state lacks a receiver]" )
                << "\n";
    }
  }

  {
    static bool const kEvoFace = Options::_env_flag( "CRONOS_AUDIT_EVOFACE", false );
    if( kEvoFace )
      _disp( 3 ) << "  [evoface] plan build: _evolution_dom_set="
                << ( _evolution_dom_set ? "y" : "n" )
                << "  _evolution_dom_user=" << ( _evolution_dom_user ? "y" : "n" )
                << "  evolution_dom="
                << ( _evolution_dom_var.dag() ? _evolution_dom_var.name()
                                              : std::string( "<none>" ) )
                << "  n_inputs=" << _mInp.size()
                << "\n";
  }
  // Two-pass DROP: pass 1 may identify redundant (column-dependent / implied)
  // continuity claims (Schur-uncoverable, and dropping them is generic-correct
  // because the constraint is implied by the covered claims).  If so we re-run
  // the whole pipeline ONCE with those claims excluded from allocation, so
  // every count (n_trace_constraint, n_trace_var, n_coll_eqn) derives from the
  // reduced, full-rank set and the eval squareness audit holds with no special
  // casing.  Cleared once here; preserved across the pass-2 re-run.  PDE2/3/4
  // are column-rank-full -> pass 1 succeeds, the loop runs once, byte-identical.
  //
  // rev60 -- PASS BUDGET 3, ACCOUNTED PER REASON.
  //
  // The budget was 2 with the retry guarded on `_pass == 0`, i.e. "one rebuild, whatever asked
  // for it".  Handoff rev7 §3 records why that is the binding obstacle for the Z-support
  // promotion rule: the split that decides eliminate-vs-keep-explicit runs BEFORE the data
  // that would inform it exists (split -> tau_claims -> weak_sat_terms -> B -> Z), so a
  // rebuild is the only mechanism available, and MMPDE-type models have already spent theirs
  // on the DROP pass by the time the projector is computed.
  //
  // The budget is therefore raised to 3 and the two reasons are accounted SEPARATELY, one
  // retry each.  This matters for inertness: a bare `_pass < 3` would also grant DROP a second
  // retry, which today is a hard failure, so it would change behaviour on any model that asks
  // twice.  With per-reason accounting, DROP still gets exactly one -- `!drop_retry_used` on
  // pass 1 is the same predicate as `_pass == 0` was, because only pass 0 can have consumed it.
  //
  // CONSEQUENCE, and it is the whole inertness argument: pass 2 is reachable ONLY if a
  // more than twice and the header is behaviourally identical to rev59.
  _strong_dropped_claims.clear();
  _exact_natural_restored.clear();
  _promotedNatural.clear();
  _keepExplicitReason.clear();
  _strong_redundancy_pending = false;
  _strong_restore_pending    = false;
  _promotionRequested = false;
  _strong_promotion_pending  = false;
  _strong_promoted_claims.clear();
  bool drop_retry_used    = false;
  bool promo_retry_used   = false;
  bool restore_retry_used = false;   // the RESTORE rebuild (see _strong_restore_pending): one, independent of DROP
  int const pass_budget = 3;
  for( int _pass = 0; _pass < pass_budget; ++_pass ){
  // rev154g: the protected set is shared with protect_interface_claims() (driver-owned,
  // persistent by contract) and with the W-decide restore; a deferral may undo only what
  // THIS pass added to it.
  std::set<t_FrozenInterfaceClaimKey> const protected_at_pass_start = _strong_protected_claims;
  draft = InterfacePlanDraft{};
  draft.state = OCPlan::State::BUILDING;
  draft.frozen = false;
  draft.n_coll_eqn = _nCollEqn;
  draft.n_coll_var = _nCollVar;
  draft.n_trace_var = 0;
  draft.trace_var_offset = _nCollVar;
  draft.n_trace_constraint = 0;
  draft.options.imposition_type = static_cast<int>( options.INTERFACE.IMPOSITION );
  draft.options.interface_type  = static_cast<int>( options.INTERFACE.TYPE );
  draft.options.sat_sigma0      = _eff_sat_sigma0();   // rev51: record what the plan USES
  draft.options.sat_sigma1      = _eff_sat_sigma1();   // rev51: ditto
  draft.options.trace_tau_scale = options.INTERFACE.TRACE_TAU_SCALE;

  bool const prev_build_allowed = _interfacePlanBuildAllowed;
  _interfacePlanBuildAllowed = true;
  try{
    auto _ib = std::chrono::steady_clock::now();
    auto _ibphase = [&]( char const* nm ){
      if( _solve_display_level() >= 2 ){
        auto const now = std::chrono::steady_clock::now();
        _disp( 2 ) << "    [iplan] " << std::left << std::setw(38) << nm << std::right
                  << std::fixed << std::setprecision(2) << std::setw(11)
                  << std::chrono::duration<double,std::milli>( now - _ib ).count()
                  << " ms" << std::defaultfloat << "\n";
        _ib = now;
      }
    };
    _pinSuppressedTriples.clear();
    _pinSuppressCount1.clear();
    _pinSuppressCount3.clear();
    // rev58 (B2 instrumentation, on rev52): report WHICH phase returned false, and on which
    // pass.  DIAGNOSTIC ONLY -- no behavioural change anywhere.
    //
    // WHY.  `ok = ok && phase(draft)` short-circuits, so a failure anywhere leaves the caller
    // with a bare false and the driver printing "interface plan build/validate failed".  Two
    // diagnostic cycles were spent inferring which phase that was, and both inferences were
    // wrong.  The tracer already names every phase for timing; this reuses the names for
    // failure at the cost of one pointer.
    //
    // Deliberately UNCONDITIONAL, not gated on DISPLAY_LEVEL: a gated promotion message in an
    // earlier revision was invisible in exactly the run that needed it, which made a retry
    // that never happened look like a retry that failed.
    char const* _failed_phase = nullptr;
    auto _phase = [&]( bool okk, char const* nm ) -> bool {
      if( !okk && !_failed_phase ) _failed_phase = nm;
      _ibphase( nm );
      return okk;
    };

    bool ok = _phase( _build_frozen_trace_constraints_and_tau_terms( draft ),
                      "trace_constraints_and_tau_terms" );
    // 2026-10-15 (WORKPLAN 4R.21): the PER-NODE RECEIVER RULE.  A claim's receiver is the row differentiating the
  // state across the face where one is collocated at the node; otherwise the state's DECLARED defining row
  // (EqnOptions::defines), taken as a NATURAL receiver -- one edge, never silenced, never rescued.  Exact modes only
  // (under IC_WEAK every rescue edge is a live SAT term by design).
  //   MEASURED: at a causal element start the evolution rows are excluded, so a state whose only differentiating
  // row is one of them has no natural receiver there.  OCFE_Burger_move --form q: q's 108 claim edges were all
  // hosted by the PDE at time nodes 1-3 of each element and NONE at node 0; q's seam value then jumped by 14 at a
  // causal start while u and x stayed continuous to 1e-9, scaling with mesh motion (14 -> 2.8 -> 3.3e-3 as tau
  // froze it).  --form r: xg's rescue sprayed edges into every row it appears in, one of them RDEF, which is q's
  // natural receiver there -- "system is singular" after tau=450.
  if( options.INTERFACE.IMPOSITION == Options::IC_TRACE || options.INTERFACE.IMPOSITION == Options::IC_STRONG ){
    // the declared defining row of each state (state id -> row id)
    std::map<size_t,size_t> const defrow = _defining_rows();   // declared, else inferred by FFInv L/S
    // HELD (2026-10-15): the same rule could cover reduce_order() auxiliaries without any declaration, by taking
    // the LINK row that mentions the auxiliary as its defining row.  Measured: it fires (36 edges promoted on
    // OCFE_Burger_move_red --red main) but does NOT make that model converge -- RED_MAIN has a further problem --
    // while it DOES change the plan of corpus models that use order reduction (OCFE_PDE20f: same 19 PASS, residuals
    // moving at round-off, 5.27e-13 -> 4.01e-13).  A behaviour change on the corpus that buys nothing yet; it ships
    // when it is shown to fix something.  Until then only a DECLARED definition (EqnOptions::defines) is used, so a
    // model that declares nothing is untouched by construction.
    if( !defrow.empty() ){
      // (claim, node) -> has a natural edge?  and the index of the declared row's edge there
      std::map< std::pair<t_FrozenInterfaceClaimKey,size_t>, bool >   has_natural;
      std::map< std::pair<t_FrozenInterfaceClaimKey,size_t>, size_t > decl_edge;
      for( size_t i = 0; i < draft.receiver_edges.size(); ++i ){
        auto const& e = draft.receiver_edges[i];
        auto const k = std::make_pair( _freeze_interface_claim_key( e.claim_key ), e.receiver_emit_flat );
        if( e.source == OCPlan::InterfaceEdgeSource::PRINCIPAL_SYMBOL ) has_natural[k] = true;
        auto dr = defrow.find( e.claim_key.state_id );
        if( dr != defrow.end() && e.row_id == dr->second ) decl_edge[k] = i;
      }
      // PER-NODE MATCHING (2026-10-16, Benoit's suggestion).  Promoting each claim's defining row on its own is not
      // injective: on OCFE_zchain (p = g, q = g', r = g'') it put natural edges of DIFFERENT claims on the same
      // (row, node) at 8 of 30 pairs, the projection lost its pivots, and the driver went from 9 PASS to 3.  So the
      // assignment is a bipartite matching per node, claims x rows, by PRIORITY: class 2 = an existing natural edge,
      // class 1 = the claim's defining row (declared or inferred).  Kuhn's augmenting path, class by class, and an
      // augmentation may re-route an already-matched claim only along an edge of at least its current class, so a
      // higher class is never traded for a lower one.  Matched class-1 edges are promoted to natural; everything
      // else is left exactly as it was (unmatched rescue edges stay rescue, to be silenced and restorable).
      size_t n_promoted = 0, n_pruned = 0;
      std::map<size_t, std::vector<size_t>> by_node;                                  // node -> edge indices
      for( size_t i = 0; i < draft.receiver_edges.size(); ++i ) by_node[ draft.receiver_edges[i].receiver_emit_flat ].push_back( i );
      std::vector<OCPlan::InterfaceReceiverEdge> kept = draft.receiver_edges;
      for( auto const& [node, idx] : by_node ){
        // class of each edge at this node: 2 natural, 1 defining row, 0 other
        auto cls = [&]( size_t i ){
          auto const& e = kept[i];
          if( e.source == OCPlan::InterfaceEdgeSource::PRINCIPAL_SYMBOL ) return 2;
          auto dr = defrow.find( e.claim_key.state_id );
          return ( dr != defrow.end() && e.row_id == dr->second )? 1: 0; };
        std::map<t_FrozenInterfaceClaimKey,size_t> cid; std::map<size_t,size_t> rid;      // claim/row -> local index
        std::vector<t_FrozenInterfaceClaimKey> claims; std::vector<size_t> rows;
        for( size_t i : idx ){
          auto k = _freeze_interface_claim_key( kept[i].claim_key );
          if( !cid.count( k ) ){ cid[k] = claims.size(); claims.push_back( k ); }
          if( !rid.count( kept[i].row_id ) ){ rid[ kept[i].row_id ] = rows.size(); rows.push_back( kept[i].row_id ); }
        }
        std::vector<int> match_row( rows.size(), -1 ), match_cls( rows.size(), -1 ), match_edge( rows.size(), -1 );
        std::vector<int> claim_row( claims.size(), -1 );
        for( int C = 2; C >= 1; --C ){
          for( size_t c = 0; c < claims.size(); ++c ){
            if( claim_row[c] >= 0 ) continue;
            std::vector<char> seen( rows.size(), 0 );
            std::function<bool(size_t)> aug = [&]( size_t cc ) -> bool {
              for( size_t i : idx ){
                if( cid[ _freeze_interface_claim_key( kept[i].claim_key ) ] != cc ) continue;
                int const k = cls( i ); if( k < C ) continue;
                size_t const r = rid[ kept[i].row_id ];
                if( seen[r] ) continue; seen[r] = 1;
                if( match_row[r] < 0 || ( match_cls[r] <= k && aug( (size_t)match_row[r] ) ) ){
                  match_row[r] = (int)cc; match_cls[r] = k; match_edge[r] = (int)i; claim_row[cc] = (int)r; return true; }
              }
              return false; };
            aug( c );
          }
        }
        for( size_t r = 0; r < rows.size(); ++r )
          if( match_edge[r] >= 0 && match_cls[r] == 1 ){
            kept[ match_edge[r] ].source = OCPlan::InterfaceEdgeSource::PRINCIPAL_SYMBOL; ++n_promoted;
            _promotedNatural.insert( _freeze_interface_claim_key( kept[ match_edge[r] ].claim_key ) ); }
      }
      if( n_promoted || n_pruned ){
        draft.receiver_edges.swap( kept );
        _disp( 2 ) << "OCFESLV::setup ** interface: per-node receiver rule -- " << n_promoted
                   << " declared-definition edge(s) promoted to natural where no differentiating row is collocated, "
                   << n_pruned << " rescue edge(s) pruned there" << std::endl;
      }
    }
  }
  ok = ok && _phase( _build_interface_receiver_graph( draft ), "interface_receiver_graph" );
    ok = ok && _phase( _build_frozen_weak_sat_terms( draft ), "frozen_weak_sat_terms" );
    ok = ok && _phase( _detect_redundant_continuity_claims( draft ),
                       "detect_redundant_continuity_claims" );
    // The projection MUTATES trace_constraints, so it must run after the drop pass (which
    // would otherwise inspect rewritten rows) and before the IC_STRONG Schur elimination
    // (which must see a full-rank multiplier block).  This slot satisfies both.  A two-pass
    // DROP retry rebuilds the draft from scratch, so the projection cannot be applied twice.
    ok = ok && _phase( _build_frozen_trace_projection( draft ), "frozen_trace_projection" );
    if( ( _solve_display_level() >= 1 || options.AUDIT.COVERAGE )
     && ( !_pinSuppressCount1.empty() || !_pinSuppressCount3.empty() ) ){
      _disp( 2 ) << "OCFESLV::setup ** interface: C0 claims suppressed --";
      for( auto const& kv : _pinSuppressCount1 )
        _disp( 2 ) << "  FIX1 " << state_name_of_id( kv.first ) << "x" << kv.second;
      for( auto const& kv : _pinSuppressCount3 )
        _disp( 2 ) << "  FIX3 " << state_name_of_id( kv.first ) << "x" << kv.second;
      _disp( 2 ) << "  (per equation-row visit, not per claim)" << std::endl;
    }
    // rev303 (S3): if an earlier phase has already marked this draft for a rebuild (the redundancy pass
    // restoring silenced duplicates, rev286/287), the elimination would run on a draft that is about to be
    // discarded.  Measured on OCFE_PDE25/PDE26 with the receiver rule on under IC_STRONG: it stalled every
    // time, printed "SELECTION INADEQUATE" -- blaming the keep-explicit heuristic for a block the restore was
    // about to change -- and the rebuild then produced the right plan.  Skip it; the rebuild is the plan.
    bool const elim_ok = _strong_redundancy_pending
                       ? ( _disp( 2 ) << "OCFESLV::setup ** tau elimination: skipped -- a rebuild is already pending"
                                        " (restored duplicates); the rebuilt draft is eliminated instead" << std::endl, false )
                       : _phase( _build_frozen_strong_tau_elimination( draft ), "frozen_strong_tau_elimination" );
    // rev289 (R3): a failed elimination with silenced duplicates present -- a silenced rescued row was an
    // eliminated tau's pivot.  Restore every silenced claim and rebuild: the plan is then today's STRONG plan.
    if( ok && !elim_ok && _knob_EXACT_NATURAL() && _knob_EXACT_NATURAL_STRONG() && !_wdecideRetry ){
      size_t n_new = 0;
      for( auto const& term : draft.weak_sat_terms )
        if( term.multiplier_silenced && _exact_natural_restored.insert( _trace_slot_key( term.claim_key ) ).second ) ++n_new;
      if( n_new ){
        _strong_restore_pending = true;    // a RESTORE, not a drop: its own retry
        std::cerr << "OCFESLV::setup ** tau elimination: EXACT_NATURAL_STRONG -- elimination failed with " << n_new
                  << " silenced claim(s); RESTORING all of them and rebuilding" << std::endl;
      }
    }
    ok = ok && elim_ok;
    _interfacePlanBuildAllowed = prev_build_allowed;
    if( ok ) return true;
    // Pass 1 flagged cleanly-redundant claims (drop) or symbol-ineligible claims
    // (keep-explicit) -> retry once with the reduced/protected set.
    // rev58: name the failing phase and the pass, before deciding whether a retry is
    // available.  This is where "interface plan build/validate failed" becomes attributable.
    // rev59 WORDING FIX.  rev58 called every false return a FAILURE.  But
    // _detect_redundant_continuity_claims returns false BY DESIGN to request the pass-2
    // rebuild -- so the normal two-pass DROP flow printed "FAILED at phase ...", which is
    // exactly the crying-wolf the instrument was added to prevent.  Distinguish the two: a
    // pending flag means "wants a retry", anything else means the build genuinely failed.
    // rev60: the retry decision is now per-reason.  `can_drop` is the same predicate rev59
    // spelled `_strong_redundancy_pending && _pass == 0` -- only pass 0 can consume the DROP
    // retry, so `!drop_retry_used` on pass 1 is false exactly where `_pass == 0` was false.
    bool const can_drop    = _strong_redundancy_pending && !drop_retry_used;
    bool const can_promo   = _strong_promotion_pending && !promo_retry_used;   // rev281: revived for STRONG_PROJECT
    bool const can_restore = _strong_restore_pending && !restore_retry_used;

    if( _failed_phase ){
      bool const wants_retry = can_drop || can_promo || can_restore;
      std::cerr << "OCFESLV::setup ** interface plan: pass " << _pass
                << ( wants_retry ? " REQUESTED A REBUILD at phase '"
                                 : " FAILED at phase '" )
                << _failed_phase << "'"
                << "  (redundancy_pending=" << ( _strong_redundancy_pending ? 1 : 0 )
                << " restore_pending=" << ( _strong_restore_pending ? 1 : 0 )
                                << " pass_budget=" << pass_budget
                << " retries_used=" << ( drop_retry_used ? "drop" : "-" )
                << "/" << ( promo_retry_used ? "promote" : "-" )
                << "/" << ( restore_retry_used ? "restore" : "-" ) << ")"
                << ( wants_retry ? "  [normal: the retry is available]"
                   : ( _strong_redundancy_pending || _strong_restore_pending )
                       ? "  [a retry is wanted but that reason's budget is SPENT]"
                       : "  [no retry requested -- a real failure]" )
                << std::endl;
    }

    // rev263: COUNT how often the escape hatch fires.  The rev262 sweep showed "legacy builder SKIPPED" on
    // all 822 plan decisions, which says the corpus did not exercise this path IN THOSE RUNS -- not that it is
    // dead.  It is the one place legacy is still load-bearing: under a GLOBAL plan the redundancy-drop cycle
    // lands on a physics-worse drop set (PDE20f cell C: 40 vs legacy's 20, e_s2 16x), so the whole setup is
    // re-run under pure legacy.  Deleting legacy before this is fixed would reinstate that 16x.
    // rev288: the deferral text is retired with the path (rev264 deleted the legacy re-setup); the rebuild
    // proceeds with the GLOBAL plan's own drop / restore set, which is what always happened since rev264.
    if( can_drop && _planIsGlobal )
      _disp( 2 ) << "OCFESLV::setup ** PLAN: drop-cycle reached with a GLOBAL plan -- rebuilding with its own drop/restore set" << std::endl;
    // rev264: the deferred re-setup to pure legacy is DELETED, with its whole body -- it re-ran the setup
    // under a builder that no longer exists.  The phase-2 trigger it carried is recorded in the workplan:
    // "make the family-member drop choice an explicit builder decision (drop the derivable/slaved member,
    // keep the primitive)", which is C0's NATURAL-vs-IMPLIED distinction.
    // NOTE: the `else if` below belonged to the deleted `if`.  Removing the block alone left it DANGLING,
    // where it bound to the diagnostic `if( can_drop && _planIsGlobal )` above and silently changed the retry
    // logic -- four PDE20f RED_MAIN cells failed setup.  It is now a plain `if`, which is what it must be once
    // the branch before it is gone.
    if( can_drop || can_promo || can_restore ){
      if( can_restore ){
        _strong_restore_pending = false;
        restore_retry_used = true;
      }
      if( can_drop ){
        _strong_redundancy_pending = false;
        drop_retry_used = true;
        if( _solve_display_level() >= 1 ){
        if( !_strong_dropped_claims.empty() )
          _disp( 2 ) << "OCFESLV::setup ** interface plan: "
                    << ( options.INTERFACE.IMPOSITION == Options::IC_TRACE ? "IC_TRACE" : "IC_STRONG" )
                    << " two-pass DROP: re-running pipeline excluding "
                    << _strong_dropped_claims.size()
                    << " redundant (implied) continuity claim(s)." << std::endl;
        else
          _disp( 2 ) << "OCFESLV::setup ** interface plan: IC_STRONG"
                    << " two-pass KEEP-EXPLICIT: re-running pipeline so step 1 routes the"
                    << " symbol-ineligible continuity claim(s) to explicit tau." << std::endl;
        }
      }
      if( can_promo ){
        promo_retry_used   = true;
        _strong_promotion_pending = false;
        // UNCONDITIONAL by design.  A promotion message gated on DISPLAY_LEVEL was invisible
        // in exactly the run that needed it once already (handoff rev7 §7), which made a
        // rebuild that never happened look like one that failed.
        _disp( 2 ) << "OCFESLV::setup ** interface plan: PROMOTION rebuild granted (pass "
                  << _pass << " -> " << ( _pass + 1 ) << " of " << pass_budget << ")."
                  << std::endl;
      }
      continue;
    }
    return false;
  }
  catch( Exceptions& e ){
    _interfacePlanBuildAllowed = prev_build_allowed;
    std::cerr << "OCFESLV::setup ** interface plan: build "
                 "failed (" << e.what() << ")" << std::endl;
    return false;
  }
  catch( ... ){
    _interfacePlanBuildAllowed = prev_build_allowed;
    std::cerr << "OCFESLV::setup ** interface plan: build "
                 "failed with an unknown exception" << std::endl;
    return false;
  }
  }   // end two-pass loop
  return false;   // pass loop exhausted (or single-pass fall-through)
}

inline bool
OCFEBASE::_build_interface_receiver_graph
( InterfacePlanDraft& draft ) const
{
  // IC_WEAK, IC_TRACE, and IC_STRONG all consume the same setup-frozen
  // interface topology.  An empty receiver graph is a valid frozen plan
  // (for example, a single-element problem), so readiness is represented by
  // the explicit plan seal rather than by container non-emptiness.
  if( draft.frozen && draft.state == OCPlan::State::FROZEN )
    return true;

  if( !_interfacePlanBuildAllowed
   || draft.state != OCPlan::State::BUILDING ){
    _disp( 2 ) << "OCFESLV::_build_interface_receiver_graph ** interface plan "
                 "construction is not active; setup() is the only owner of "
                 "the frozen receiver graph." << std::endl;
    return false;
  }

  try{
    std::map<t_InterfaceClaimKey,size_t> claims;
    _prepare_interface_plan_tables( claims, draft );

    // IC_TRACE/IC_STRONG are sealed by _store_frozen_interface_plan().
    // IC_WEAK has no exact trace/tau plan to store, but the receiver graph
    // produced above is still a complete frozen plan, including the valid
    // no-interface case where the receiver-edge vector is empty.
    if( options.INTERFACE.IMPOSITION == Options::IC_WEAK ){
      draft.trace_constraints.clear();
      draft.trace_constraints_ready = true;
      draft.frozen = true;
      draft.state  = OCPlan::State::FROZEN;
    }
  }
  catch( Exceptions& e ){
    std::cerr << "OCFESLV::_build_interface_receiver_graph ** setup-owned "
                 "interface receiver graph failed (" << e.what() << ")"
              << std::endl;
    return false;
  }
  catch( ... ){
    std::cerr << "OCFESLV::_build_interface_receiver_graph ** setup-owned "
                 "interface receiver graph failed" << std::endl;
    return false;
  }
  return true;
}

inline bool
OCFEBASE::_build_frozen_weak_sat_terms
( InterfacePlanDraft& draft ) const
{
  draft.weak_sat_terms.clear();

  size_t const invalid_claim_id = std::numeric_limits<size_t>::max();
  auto trace_slot_key = [&]( t_InterfaceClaimKey key ) -> t_InterfaceClaimKey
  {
    key.block_id = _trace_claim_block( key.state_id, key.dom_id );
    key.eqn_id   = invalid_claim_id;
    if( key.kind == CLAIM_EXACT_ROW ) key.kind = CLAIM_EXACT_C0;
    return key;
  };

  auto dom_by_id = [&]( size_t did ) -> FFDom const*
  {
    for( auto const& [dv, dom] : _mDom )
      if( static_cast<size_t>( dv.id().second ) == did ) return &dom;
    return nullptr;
  };

  // The receiver graph is the authoritative frozen list of weak-SAT receiver
  // rows.  Store a fully executable consumer projection so numeric SAT/TRACE
  // evaluation applies the setup-owned plan directly: no lambda maps, receiver
  // slot sets, or topology lookups are rebuilt in _eval_eqn().
  draft.weak_sat_terms.reserve( draft.receiver_edges.size() );
  // IC_WEAK WITH MULTIPLIERS, natural-only: which claims have a natural (principal-symbol) receiver that can
  // carry the multiplier.  Their rescued receivers are duplicates with a fabricated weight (NOTES_20260917a M2).
  std::set<t_FrozenInterfaceClaimKey> natural_tau;
  size_t n_silenced = 0;
  // C2b (rev278): the same balance-receiver rule for the MULTIPLIER under IC_TRACE / IC_STRONG.
  int const kExactNat = _knob_EXACT_NATURAL();
  // rev280: IC_TRACE ONLY.  IC_STRONG eliminates each tau through a receiving row that contains it (the Schur
  // block of _build_frozen_strong_tau_elimination); a rescued row is that tau's elimination row, not a
  // duplicate, and silencing it fails the elimination phase (sweep G 2026-09-17: every k'>0 line and every
  // failure was IC_STRONG at 'frozen_strong_tau_elimination'; IC_TRACE clean).
  bool const exact_nat = kExactNat && ( options.INTERFACE.IMPOSITION == Options::IC_TRACE
        || ( options.INTERFACE.IMPOSITION == Options::IC_STRONG && _knob_STRONG_PROJECT() && _knob_EXACT_NATURAL_STRONG() ) );   // rev289
  std::set<t_FrozenInterfaceClaimKey> natural_exact;
  size_t n_exact_silenced = 0, n_exact_sole = 0;
  if( exact_nat ){
    std::set<t_FrozenInterfaceClaimKey> all_exact;
    for( auto const& edge : draft.receiver_edges ){
      if( !edge.weak_sat_allowed || !edge.trace_tau_allowed ) continue;
      auto const k = _freeze_interface_claim_key( trace_slot_key( edge.claim_key ) );
      all_exact.insert( k );
      if( edge.source != OCPlan::InterfaceEdgeSource::PRINCIPAL_SYMBOL ) continue;
      if( kExactNat >= 2 && edge.row_id < _mEqn.size() && _mEqn[ edge.row_id ].opt->role == EqnRole::LINK ) continue;
      if( _promotedNatural.count( k ) ) continue;   // promoted: natural for the pivot, not for silencing siblings
      natural_exact.insert( k );
    }
    for( auto const& k : all_exact ) if( !natural_exact.count( k ) ) ++n_exact_sole;
    _disp( 3 ) << "OCFESLV::setup ** PLAN: CRONOS_EXACT_NATURAL=" << kExactNat
              << ( kExactNat >= 2 ? " (multiplier on natural BALANCE receivers only)"
                                  : " (multiplier on natural receivers only)" ) << std::endl;
  }
  else if( kExactNat )
    _disp( 3 ) << "OCFESLV::setup ** PLAN: CRONOS_EXACT_NATURAL=" << kExactNat
              << ( options.INTERFACE.IMPOSITION == Options::IC_STRONG
                   ? " -- inactive under IC_STRONG (set CRONOS_EXACT_NATURAL_STRONG=1 with STRONG_PROJECT to enable, rev289)"
                   : " -- inactive: imposition is not IC_TRACE" ) << std::endl;
  if( draft.weak_multipliers && draft.weak_natural_only )
    for( auto const& edge : draft.receiver_edges )
      if( edge.weak_sat_allowed && edge.trace_tau_allowed
       && edge.source == OCPlan::InterfaceEdgeSource::PRINCIPAL_SYMBOL )
        natural_tau.insert( _freeze_interface_claim_key( trace_slot_key( edge.claim_key ) ) );
  // C2 EXPERIMENT (rev274) -- NATURAL-ONLY PENALTY.  CRONOS_WEAK_NATURAL_PENALTY=1, IC_WEAK without multipliers:
  // a claim with a natural (principal-symbol) weak receiver is penalised ONLY there, the standard interior
  // penalty with the symbol as weight; its rescued receiver terms -- a second injection into an algebraic row,
  // whose weight the model does not define (NOTES_20260917a M2) -- are DROPPED.  A SOLE claim, with no natural
  // receiver, keeps its rescued terms unchanged (dropping them drops the claim: rev255, PDE20f C 44x).
  // No columns, no rows, nVar unchanged.  The knob is read as a value (=0 restores) and printed whenever set.
  int const kNatPen = _knob_WEAK_NATURAL_PENALTY();   // rev276: an Options field, read per setup
  bool const nat_pen = kNatPen && options.INTERFACE.IMPOSITION == Options::IC_WEAK && !draft.weak_multipliers;
  std::set<t_FrozenInterfaceClaimKey> natural_pen, all_claims_pen, natural_link_only;
  size_t n_dropped = 0, n_sole_terms = 0;
  if( kNatPen ){
    _disp( 3 ) << "OCFESLV::setup ** PLAN: CRONOS_WEAK_NATURAL_PENALTY=" << kNatPen
              << ( nat_pen && kNatPen >= 2 ? " (penalty on natural BALANCE receivers only; a claim whose natural receivers are all LINK rows keeps every term)"
                 : nat_pen ? " (penalty on natural receivers only; sole claims keep rescued rows)"
                 : options.INTERFACE.IMPOSITION != Options::IC_WEAK ? " -- inactive: imposition is not IC_WEAK"
                 : " -- inactive: IC_WEAK with multipliers (CRONOS_WEAK_TAU_RESCUED) takes precedence" ) << std::endl;
    if( nat_pen )
      for( auto const& edge : draft.receiver_edges ){
        if( !edge.weak_sat_allowed ) continue;
        auto const k = _freeze_interface_claim_key( trace_slot_key( edge.claim_key ) );
        all_claims_pen.insert( k );
        if( edge.source != OCPlan::InterfaceEdgeSource::PRINCIPAL_SYMBOL ) continue;
        // rev275 (=2): a natural receiver that is a LINK row is the NEXT auxiliary's definition, not a balance --
        // the penalty there pollutes that auxiliary by sigma*tau*[jump] and barely controls the jump (PDE12,
        // NOTES_20260917f).  Such a receiver does not qualify the claim for natural-only treatment.
        bool const is_link = ( edge.row_id < _mEqn.size() && _mEqn[ edge.row_id ].opt->role == EqnRole::LINK );
        if( kNatPen >= 2 && is_link ){ natural_link_only.insert( k ); continue; }
        natural_pen.insert( k );
      }
    for( auto const& k : natural_pen ) natural_link_only.erase( k );   // has a balance receiver too
  }
  for( auto const& edge : draft.receiver_edges ){
    if( !edge.weak_sat_allowed ) continue;
    if( nat_pen && edge.source != OCPlan::InterfaceEdgeSource::PRINCIPAL_SYMBOL ){
      if( natural_pen.count( _freeze_interface_claim_key( trace_slot_key( edge.claim_key ) ) ) ){ ++n_dropped; continue; }
      ++n_sole_terms;
    }

    FrozenWeakSatTerm term;
    term.row_id            = edge.row_id;
    term.claim_key         = edge.claim_key;
    term.sat_type          = edge.sat_type;
    term.source            = edge.source;
    term.receiver_dom_id   = edge.receiver_dom_id;
    term.receiver_iel      = edge.receiver_iel;
    term.receiver_inode    = edge.receiver_inode;
    term.receiver_flat     = edge.receiver_flat;
    term.receiver_emit_flat = edge.receiver_emit_flat;
    term.receiver_block_el = edge.receiver_block_el;
    term.coupling          = edge.coupling;
    term.orientation       = edge.orientation;
    term.trace_tau_allowed = edge.trace_tau_allowed;
    term.uses_c1           = ( edge.claim_key.kind == CLAIM_EXACT_C1 );

    FFDom const* face_dom = dom_by_id( edge.claim_key.dom_id );
    if( !face_dom ){
      std::cerr << "OCFESLV::_build_frozen_weak_sat_terms ** unknown face "
                   "domain in frozen receiver edge: dom_id="
                << edge.claim_key.dom_id << std::endl;
      return false;
    }
    if( edge.claim_key.iel_lo + 1 >= face_dom->n_elem ){
      std::cerr << "OCFESLV::_build_frozen_weak_sat_terms ** invalid interface "
                   "element in frozen receiver edge: dom_id="
                << edge.claim_key.dom_id
                << " iel_lo=" << edge.claim_key.iel_lo
                << " n_elem=" << face_dom->n_elem << std::endl;
      return false;
    }

    // rev51: _eff_sat_sigma0/1 apply the B1 audit override (CRONOS_SAT_SIGMA0/1) if set,
    // otherwise return options.INTERFACE.SAT_SIGMA0/1 unchanged.  This is the ONLY consumption site.
    double tau = _eff_sat_sigma0() * _sat_tau_face( *face_dom, edge.claim_key.iel_lo );
    if( term.uses_c1 ){
      double const iface_w = _sat_iface_width( *face_dom, edge.claim_key.iel_lo );
      tau = ( iface_w > 0. )
          ? ( _eff_sat_sigma1() * _sat_tau_face( *face_dom, edge.claim_key.iel_lo )
            / iface_w )
          : 0.;
    }

    term.prefactor = tau * term.coupling * term.orientation;   // (rev320: the IC_CONS flux-form branch retired)
    term.trace_prefactor = term.prefactor * options.INTERFACE.TRACE_TAU_SCALE;

    // Pre-compile the sparse SAT coefficient list for IC_WEAK evaluation.
    // eval time: buf[flat] += sum_k sat_coeffs[k].second * var[sat_coeffs[k].first]
    if( !_sat_claim_coeffs( edge.claim_key, term.prefactor, term.sat_coeffs ) ){
      std::cerr << "OCFESLV::_build_frozen_weak_sat_terms ** failed to build "
                   "SAT coefficients for receiver term: row=" << term.row_id
                << " state=" << edge.claim_key.state_id
                << " dom=" << edge.claim_key.dom_id
                << " iel_lo=" << edge.claim_key.iel_lo
                << " face=" << edge.claim_key.face
                << " kind=" << edge.claim_key.kind << std::endl;
      return false;
    }

    if( exact_nat && term.trace_tau_allowed
     && edge.source != OCPlan::InterfaceEdgeSource::PRINCIPAL_SYMBOL
     && natural_exact.count( _freeze_interface_claim_key( trace_slot_key( term.claim_key ) ) )
     && !_exact_natural_restored.count( _freeze_interface_claim_key( trace_slot_key( term.claim_key ) ) ) ){
      // C2b: a duplicate of a natural receiver -- the multiplier does not enter here, so nothing reads the
      // rescue's fabricated weight in this row.  Sole claims (no natural receiver) fall through and keep theirs.
      term.trace_tau_allowed = false;
      term.multiplier_silenced = true;
      ++n_exact_silenced;
    }
    else if( ( options.INTERFACE.IMPOSITION == Options::IC_TRACE
       || options.INTERFACE.IMPOSITION == Options::IC_STRONG ) && term.trace_tau_allowed ){
      t_InterfaceClaimKey const slot = trace_slot_key( term.claim_key );
      auto it = draft.tau_claims.find( _freeze_interface_claim_key( slot ) );
      if( it == draft.tau_claims.end() || it->second.empty() ){
        _disp( 2 ) << "OCFESLV::_build_frozen_weak_sat_terms ** retained "
                     "IC_TRACE/IC_STRONG receiver term has no lambda slot: row="
                  << term.row_id
                  << " block=" << slot.block_id
                  << " state=" << slot.state_id
                  << " dom=" << slot.dom_id
                  << " iel=" << slot.iel_lo
                  << " face=" << slot.face
                  << " kind=" << slot.kind << std::endl;
        return false;
      }
      term.claim_key  = slot;
      term.lambda_col = it->second.front();
      term.lambda_index = ( options.INTERFACE.IMPOSITION == Options::IC_TRACE )
        ? term.lambda_col - draft.trace_var_offset
        : term.lambda_col;
    }
    else if( draft.weak_multipliers ){
      // IC_WEAK WITH MULTIPLIERS: a term whose claim has a tau is imposed EXACTLY -- marked, and when it may
      // carry the multiplier it takes the IC_TRACE column and slot convention.  A tau-allowed term of a
      // penalised claim cannot occur (the builder switched those edges off), so it is refused loudly.
      t_InterfaceClaimKey const slot = trace_slot_key( term.claim_key );
      auto it = draft.tau_claims.find( _freeze_interface_claim_key( slot ) );
      if( it != draft.tau_claims.end() && !it->second.empty() ){
        term.exact_claim = true;
        if( draft.weak_natural_only && term.source != OCPlan::InterfaceEdgeSource::PRINCIPAL_SYMBOL
         && natural_tau.count( it->first ) ){
          term.trace_tau_allowed = false;   // a duplicate of a natural receiver: silent
          ++n_silenced;
        }
        else if( term.trace_tau_allowed ){
          term.claim_key    = slot;
          term.lambda_col   = it->second.front();
          term.lambda_index = term.lambda_col - draft.trace_var_offset;
        }
      }
      else if( term.trace_tau_allowed ){
        _disp( 2 ) << "OCFESLV::_build_frozen_weak_sat_terms ** IC_WEAK+TAU: tau-allowed receiver term of a "
                     "penalised claim: row=" << term.row_id << " state=" << slot.state_id
                  << " dom=" << slot.dom_id << " iel=" << slot.iel_lo << " face=" << slot.face
                  << " kind=" << slot.kind << std::endl;
        return false;
      }
    }

    draft.weak_sat_terms.push_back( term );
  }

  if( exact_nat )
    _disp( 2 ) << "OCFESLV::setup ** PLAN: exact natural-only -- " << n_exact_silenced
              << " rescued receiver term(s) silenced; " << natural_exact.size()
              << " claim(s) fed by natural receivers, " << n_exact_sole << " by rescued receivers only"
              << ( _exact_natural_restored.empty() ? std::string()
                 : "; " + std::to_string( _exact_natural_restored.size() ) + " restored (column-dependent with the natural receiver alone)" )
              << std::endl;
  if( nat_pen )
    _disp( 2 ) << "OCFESLV::setup ** PLAN: IC_WEAK natural-only penalty -- " << n_dropped
              << " rescued receiver term(s) dropped; " << natural_pen.size() << " claim(s) penalised on natural receivers, "
              << ( all_claims_pen.size() - natural_pen.size() ) << " sole claim(s) keep " << n_sole_terms
              << " rescued term(s)" << ( kNatPen >= 2 ? " (of which " + std::to_string( natural_link_only.size() )
                            + " natural-LINK-only: a chain auxiliary whose only natural receiver is the NEXT link"
                              " row -- see OCFE_PDE12's documented RED_MAIN cells for what an unshared chain costs)"
                          : std::string() )
              << "; weak_sat_terms=" << draft.weak_sat_terms.size() << std::endl;
  if( draft.weak_multipliers ){
    size_t n_sole = 0;
    for( auto const& kv : draft.tau_claims ) if( !natural_tau.count( kv.first ) ) ++n_sole;
    _disp( 2 ) << "OCFESLV::setup ** PLAN: IC_WEAK+TAU receivers -- "
              << ( draft.weak_natural_only ? "natural-only: " : "all tau-allowed: " )
              << n_silenced << " rescued receiver term(s) silenced";
    if( draft.weak_natural_only )
      _disp( 2 ) << "; " << ( draft.tau_claims.size() - n_sole ) << " exact claim(s) fed by natural receivers, "
                << n_sole << " by rescued receivers only";
    _disp( 2 ) << std::endl;
  }
  draft.weak_sat_terms_ready = true;
  return true;
}

inline bool
OCFEBASE::_build_frozen_trace_constraints_and_tau_terms
( InterfacePlanDraft& draft ) const
{
  draft.n_trace_var        = 0;
  draft.trace_var_offset   = draft.n_coll_var;
  draft.n_trace_constraint = 0;
  draft.claims.clear();
  draft.tau_claims.clear();
  draft.receiver_edges.clear();
  draft.weak_sat_terms.clear();
  draft.weak_sat_terms_ready = false;
  draft.trace_constraints.clear();
  draft.trace_constraints_ready = false;
  draft.exact_replacement_rows.clear();
  draft.strong_tau_elim.clear();
  draft.frozen = false;

  if( draft.state != OCPlan::State::BUILDING ){
    _disp( 2 ) << "OCFESLV::setup ** trace/tau plan: "
                 "interface draft construction is not active; setup() is the "
                 "only owner of the frozen trace/tau plan." << std::endl;
    return false;
  }

  // rev260 (C3 step 1): THIS is what makes IC_WEAK decision-free.  Weak returns here, before
  // _store_frozen_interface_plan and therefore before any builder runs -- so under weak the plan is the raw
  // draft (receiver edges + weak-SAT terms) committed directly, and NEITHER builder decides anything.  That is
  // why the rev258 weak-TAU edit and the rev259 apply knob both did nothing: there was no decision layer to
  // change.  Letting weak fall through gives it one; the builder's weak branch already assigns SAT to every
  // claim, so the plan must be UNCHANGED -- that is the acceptance criterion, and the POSITIVE signal that it
  // ran is the builder-invoked count rising.
  // (CRONOS_WEAK_DECISION_LAYER, the opt-out restoring that early return, is retired 2026-10-07: weak always falls
  // through to the decision layer, as by default since rev262)

  // Build the retained exact trace basis directly during setup.  This is the
  // sole source of trace-row count, tau-column assignment, and executable
  // FrozenTraceConstraintTerm rows.  No residual-evaluation dry run or
  // trace-count side channel is used.
  try{
    std::map<t_InterfaceClaimKey,size_t> claims;
    _prepare_interface_plan_tables( claims, draft, false );

    // Two-pass DROP (pass 2): remove the redundant (column-dependent / implied)
    // continuity claims identified in pass 1 BEFORE any count is taken, so
    // n_trace_constraint, n_trace_var, exact_replacement_rows, and n_coll_eqn
    // all derive from the reduced, full-rank claim set and the squareness
    // invariants hold automatically.  Match by the same trace-slot canonical
    // form used when the dropped keys were recorded (block via _trace_claim_
    // block, eqn_id cleared, EXACT_ROW->EXACT_C0).  No-op when the set is empty
    // (every column-rank-full plan), so PDE2/3/4 allocate exactly as before.
    // rev152d: two-pass DROP (pass 2), n_trace_constraint and the IC_TRACE n_trace_var
    // reset now happen inside OCPlanLegacyBuilder::build; the builder gets the full
    // claims map and PlanInput::dropped_claims.
    _apply_c0_claim_filter( claims );                                 // rev218 (C1), no-op unless CRONOS_CLAIM_C0
    PlanInput const plan_in = _pack_plan_input( claims, draft );      // rev152a/b: the contract
    _store_frozen_interface_plan( claims, draft, plan_in );
  }
  catch( Exceptions& e ){
    std::cerr << "OCFESLV::setup ** trace/tau plan: "
                 "trace-plan construction failed (" << e.what() << ")"
              << std::endl;
    return false;
  }
  catch( ... ){
    std::cerr << "OCFESLV::setup ** trace/tau plan: "
                 "trace-plan construction failed (unknown exception)"
              << std::endl;
    return false;
  }

  return true;
}

inline bool
OCFEBASE::_block_drop_eligible
( int block_id ) const
{
  // Item 11 (principal-symbol block-eligibility gate).  A continuity claim is
  // drop-eligible ONLY when the block's principal symbol over-specifies it:
  //   * parabolic / parabolic-closure -> a flux/derivative continuity is implied
  //     by the value continuity + the constitutive (LINK) relation -> eligible;
  //   * scalar hyperbolic (<=1 real characteristic per sampled direction) -> the
  //     trace is fixed by the one upwind characteristic -> eligible;
  //   * system hyperbolic (>=2 real characteristics) -> each characteristic
  //     carries independent continuity -> NOT eligible (keep all);
  //   * degenerate / unresolved / merely-spatially-characteristic / unknown
  //     -> NOT eligible (keep all).
  // Shared by the redundancy gate AND the value-slaved drop (both in
  // _detect_redundant_continuity_claims): a single source of truth so neither
  // path can drop a claim in an ineligible block.  Conservative -- when in doubt
  // it returns false (keep the claim).
  auto it = _blockClassification.find( block_id );
  if( it == _blockClassification.end() ) return false;   // unknown -> keep all
  t_Classify const& cls = it->second;
  if( cls.degenerate || cls.dn_elliptic ) return false;  // unresolved symbol -> keep all (dn_elliptic: report-only verdict, 2026-10-04)
  // rev165: DECOUPLE the classification fix from the drop unlock.
  // STAGE C gives rectangular blocks a character type, which sets
  // parabolic_structure_detected and therefore makes them drop-eligible for the first time.
  // MEASURED (OCFE_PDE25, LGL exact modes, 6 phases x 4 magnitudes x 2 modes = 48 cells per
  // configuration): with the drop unlocked the convergence basin collapses --
  //     pert   0.01   0.05   0.10   0.20
  //     STAGE B (no unlock)  12/12  12/12  12/12  12/12
  //     STAGE C (unlocked)   12/12  10/12   3/12   0/12
  // The drop is rank-justified (def0=0, ALL-IMPLIED) and the manufactured-exact residual still
  // passes, so this is a solver-robustness cost, not a correctness one -- but it is a factor
  // >=4 in basin radius and it is not acceptable by default.
  // CRONOS_RECT_DROP=1 restores rev164 behaviour (rectangular blocks drop-eligible).
  if( cls.parabolic_structure_detected ){
    static int const kRectDrop = 0;   // CRONOS_RECT_DROP (retired 2026-10-07, WORKPLAN 3.B batch 2b)
    if( cls.symbol_rectangular && !kRectDrop ) return false;   // type fixed, drop NOT unlocked
    return true;                                               // flux over-specified -> eligible
  }
  if( cls.evolution_hyperbolic ){
    size_t maxreal = 0;
    for( auto const& ed : cls.eigendata ){
      size_t nr = 0;
      for( arma::uword i = 0; i < ed.second.n_elem; ++i ){
        double const re = ed.second(i).real(), im = ed.second(i).imag();
        if( std::fabs( im ) <= 1e-9 * std::max( 1.0, std::fabs( re ) ) ) ++nr;
      }
      maxreal = std::max( maxreal, nr );
    }
    return maxreal <= 1;                                 // scalar hyperbolic -> eligible
  }
  return false;   // spatially-characteristic-only / weak-hyperbolic / undetermined -> keep all
}

inline bool

//==========================================================================================================
// INTERFACE PLAN -- claims -> realisations -> rows and columns.  The DECISIONS live in ocplan.hpp; what is here
// is the assembly: which claims exist, which rows receive them, and how the chosen realisation becomes
// matrix structure.
// (banner placed from the definitions themselves: this region is lines 16481-17045 of the pre-rev268 file)
//==========================================================================================================

OCFEBASE::_detect_redundant_continuity_claims
( InterfacePlanDraft& draft ) const
{
  // Only the exact-imposition modes allocate continuity taus; IC_WEAK uses SAT
  // penalties (no tau columns) and has nothing to drop.
  if( options.INTERFACE.IMPOSITION != Options::IC_TRACE
   && options.INTERFACE.IMPOSITION != Options::IC_STRONG )
    return true;
  // DROP_OFF (item 10 step 3): disable the redundant-continuity drop entirely;
  // every continuity claim is kept explicit.  The most conservative baseline.
  int const dp_env = _knob_DROP_POLICY();
  auto const dp = dp_env < 0 ? options.INTERFACE.DROP_POLICY
                             : static_cast<Options::InterfaceDropPolicy>( dp_env );
  if( dp == Options::DROP_OFF )
    return true;
  if( draft.weak_sat_terms.empty() ) return true;

  auto make_row_key = []( FrozenWeakSatTerm const& t ){
    FrozenPhysicalResidualRowKey r;
    r.row_id   = t.row_id;
    r.block_el = t.receiver_block_el;
    r.flat     = t.receiver_emit_flat;
    return r;
  };

  // Per-tau-column receiver coupling, keyed by lambda_index.  lambda_index is a
  // 0-based tau index in BOTH modes (IC_TRACE stores lambda_col-trace_var_offset,
  // IC_STRONG stores the eliminated index), so this reproduces the in-Schur
  // column scheme exactly while remaining mode-agnostic.  Each column also
  // records its claim key (already the canonical trace-slot form on the term).
  std::map<size_t,std::map<FrozenPhysicalResidualRowKey,double>> col;
  std::map<size_t,t_FrozenInterfaceClaimKey>                     col2claim;
  for( auto const& term : draft.weak_sat_terms ){
    if( !term.trace_tau_allowed ) continue;
    if( term.lambda_index == std::numeric_limits<size_t>::max() ) continue;
    // Skip the IC_STRONG explicit-hybrid (all-LINK kept) taus, which live at
    // var[] columns >= trace_var_offset: they are handled by the keep-explicit
    // hybrid (benign null, e.g. PDE2's 16), not by the drop, and must be left
    // untouched.  IC_TRACE lambda_index values are all below trace_var_offset,
    // so this skips nothing there.
    if( term.lambda_index >= draft.trace_var_offset ) continue;
    if( std::fabs( term.trace_prefactor ) <= 64. * DBL_EPSILON ) continue;
    col[ term.lambda_index ][ make_row_key( term ) ] += term.trace_prefactor;
    col2claim[ term.lambda_index ] = _freeze_interface_claim_key( term.claim_key );
  }
  if( col.size() < 2 ) return true;

  // Stable row index over the union of receiver rows.
  std::map<FrozenPhysicalResidualRowKey,size_t> row_id;
  for( auto const& kv : col )
    for( auto const& rv : kv.second ) row_id.emplace( rv.first, row_id.size() );


  // Greedy column Gram-Schmidt in lambda_index (allocation) order: the std::map
  // iterates ascending, so the independent (lower-index) claims enter the basis
  // and the dependent (higher-index) ones are flagged -- the same selection the
  // in-Schur detector made, hence the same drop set on the validated cases.
  // Sparse modified Gram-Schmidt (mirrors the trace rank-filter in
  // _prepare_interface_plan_tables).  Each column touches only a few receiver
  // rows, so dense arma::vec(m) dot/axpy iterate mostly zeros -- O(cols*rank*m).
  // Keep columns sparse (compact-row-index sorted); disjoint-support columns
  // give an exact-zero dot and are skipped, so 'c -= 0*q' (a dense no-op) is
  // matched bit-for-bit while the cost scales with real support overlap.
  typedef std::vector<std::pair<size_t,double>> SparseVec;
  auto sparse_dot = []( SparseVec const& a, SparseVec const& b ) -> double {
    double d = 0.0; size_t i = 0, j = 0;
    while( i < a.size() && j < b.size() ){
      if( a[i].first < b[j].first ) ++i;
      else if( a[i].first > b[j].first ) ++j;
      else { d += a[i].second * b[j].second; ++i; ++j; }
    }
    return d;
  };
  auto sparse_axpy = []( SparseVec& r, double s, SparseVec const& q ){
    SparseVec out; out.reserve( r.size() + q.size() );   // r := r - s*q
    size_t i = 0, j = 0;
    while( i < r.size() && j < q.size() ){
      if( r[i].first < q[j].first ){ out.push_back( r[i] ); ++i; }
      else if( r[i].first > q[j].first ){ out.emplace_back( q[j].first, -s*q[j].second ); ++j; }
      else { out.emplace_back( r[i].first, r[i].second - s*q[j].second ); ++i; ++j; }
    }
    for( ; i < r.size(); ++i ) out.push_back( r[i] );
    for( ; j < q.size(); ++j ) out.emplace_back( q[j].first, -s*q[j].second );
    r.swap( out );
  };
  if( std::getenv( "CRONOS_PLAN_DUMP_DETCOLS" ) ){
    static size_t det_invocation = 0; ++det_invocation;
    _disp( 3 ) << "  [plan/detinv] inv=" << det_invocation
              << " cols=" << col.size()
              << " dropped=" << _strong_dropped_claims.size()
              << " protected=" << _strong_protected_claims.size() << std::endl;
    for( auto const& kv : col ){
      auto const& fk = col2claim[ kv.first ];
      _disp( 3 ) << "  [plan/detcol] lam=" << kv.first
                << " blk=" << std::get<0>( fk ) << " state=" << std::get<1>( fk )
                << " eqn=" << (long long)std::get<2>( fk )
                << " dom=" << std::get<3>( fk )
                << " iel=" << std::get<4>( fk ) << " face=" << std::get<5>( fk )
                << " kind=" << std::get<6>( fk ) << " :";
      for( auto const& rv : kv.second )
        _disp( 3 ) << "  r" << rv.first.row_id << "." << rv.first.flat
                  << "=" << std::scientific << std::setprecision( 3 ) << rv.second;
      _disp( 3 ) << std::endl;
    }
  }
  std::vector<SparseVec> basis;
  std::vector<t_FrozenInterfaceClaimKey> basis_claim;                       // rev286: claim of each basis member
  std::map<t_FrozenInterfaceClaimKey, std::vector<size_t>> dep_group;       // rev286: flagged claim -> basis members it reduced against
  std::set<t_FrozenInterfaceClaimKey> dep;
  for( auto const& kv : col ){
    std::vector<size_t> touched;                                            // rev286
    SparseVec c;
    c.reserve( kv.second.size() );
    for( auto const& rv : kv.second )
      c.emplace_back( row_id[ rv.first ], rv.second );
    std::sort( c.begin(), c.end(),
               []( std::pair<size_t,double> const& x,
                   std::pair<size_t,double> const& y ){ return x.first < y.first; } );
    double n0 = 0.0;
    for( auto const& e : c ) n0 += e.second*e.second;
    n0 = std::sqrt( n0 );
    if( n0 <= 64. * DBL_EPSILON ){ dep.insert( col2claim[ kv.first ] ); continue; }
    for( size_t qi = 0; qi < basis.size(); ++qi ){
      auto const& q = basis[qi];
      double const dot = sparse_dot( c, q );
      if( dot == 0.0 ) continue;     // disjoint support -> dense 'c -= 0*q' no-op
      sparse_axpy( c, dot, q );
      touched.push_back( qi );                                              // rev286
    }
    double nr = 0.0;
    for( auto const& e : c ) nr += e.second*e.second;
    nr = std::sqrt( nr );
    if( nr <= 1e-9 * std::max( 1.0, n0 ) ){ dep.insert( col2claim[ kv.first ] ); dep_group[ col2claim[ kv.first ] ] = touched; continue; }
    double const inv = 1.0 / nr;
    for( auto& e : c ) e.second *= inv;
    basis.push_back( std::move( c ) );
    basis_claim.push_back( col2claim[ kv.first ] );                         // rev286
  }

  // Set by the IC_STRONG arm of the block-eligibility gate below when it NEWLY
  // protects a claim.  _store_frozen_interface_plan (pipeline step 1, which hosts
  // the keep-explicit hybrid) has ALREADY run this pass without that protection,
  // so the claim was routed to the Schur corner; the plan must be re-derived once
  // (reusing the two-pass mechanism) so step 1 sees the protection and keeps it
  // explicit instead.  Cleared between setup() calls with _strong_redundancy_pending.
  [[maybe_unused]] bool gate_newly_protected = false;

  // Default-ON (promoted after validation: no-op on PDE0-4/6, fixes PDE7 IC_TRACE).
  // Compile with -DMC__OCFESLV_SYMBOL_DROP_GATE_OFF to disable for A/B comparison.
#ifndef MC__OCFESLV_SYMBOL_DROP_GATE_OFF
  // -------------------------------------------------------------------------
  // Item 11 (principal-symbol block-eligibility gate, setup-time).  The
  // receiver-B column-dependence proxy cannot tell a verify-sound drop from an
  // over-drop within a block (the within-(state,direction) selection is a global
  // rank fact through the singular physical operator, not a per-claim symbol
  // property).  But the BLOCK principal symbol soundly bounds drop-eligibility:
  //   * parabolic / parabolic-closure block -> a flux/derivative continuity is
  //     over-specified by the value continuity + the constitutive (LINK) relation
  //     -> drop-eligible;
  //   * scalar hyperbolic (a single real characteristic in every sampled spatial
  //     direction) -> the trace is fixed by the one upwind characteristic, so the
  //     explicit value continuity is redundant -> drop-eligible;
  //   * system hyperbolic (>=2 real characteristics, e.g. a 2x2 wave with an
  //     incoming AND an outgoing field) -> each characteristic carries independent
  //     continuity -> NOT eligible (keep all);
  //   * degenerate / unresolved / merely-spatially-characteristic -> NOT eligible.
  // Vetoing the receiver-B flags in non-eligible blocks removes the over-drop
  // (e.g. PDE7's 2x2 wave) at SETUP without any solve, so the verify-driven
  // re-derive never has to re-size the plan (plan invariant 2.5).  The gate is
  // conservative: it only ever shrinks the flagged set, so it can under-drop
  // (safe: a kept redundant claim) but never introduces an over-drop.
  //
  // SCOPE: both IC_TRACE and IC_STRONG.  The gate decides the claim set entirely
  // at SETUP, so neither path needs the post-solve verify-driven re-derive loop
  // (verify becomes a passive confirmation).  Routing of the symbol-ineligible
  // flags is mode-specific because the taus live in different places:
  //   * IC_TRACE: taus are in var[].  A flagged claim is kept as an ORDINARY
  //     trace constraint -- just un-flag it (erase from dep) so it is not dropped.
  //     The saddle absorbs the (benign) redundancy; nVar is the un-dropped plan's
  //     size, fixed once at setup -- no protect, no re-size (plan invariant 2.5).
  //   * IC_STRONG: taus are internal (Schur-eliminated, not in var[]).  A system-
  //     hyperbolic flag is a tau the receiver block cannot cover by an ordinary
  //     row, so it MUST be kept EXPLICIT (a bare veto leaves it uncovered and
  //     setup fails the tau-rank check).  Protect it HERE -- the same claim set
  //     and the same canonical key form verify_and_protect would pin post-solve --
  //     and leave it in dep so the protected-loop below skips the drop.  The
  //     keep-explicit hybrid in _build_frozen_strong_tau_elimination (next in the
  //     pipeline) then routes it to an explicit tau.  No drop -> no over-drop ->
  //     no re-derive; nVar is invariant under the drop either way.
  // The gate is conservative: it only ever moves a flag from "drop" to "keep"
  // (ordinary, IC_TRACE; explicit, IC_STRONG), so it can under-drop (safe: a kept
  // redundant claim) but never introduces an over-drop.
  {
    // Eligibility predicate hoisted to OCFESLV::_block_drop_eligible (single source
    // of truth shared with the value-slaved drop below).
    size_t const n_raw = dep.size();
    std::vector<t_FrozenInterfaceClaimKey> ineligible;
    for( auto const& ck : dep ){
      t_InterfaceClaimKey const k = _unfreeze_interface_claim_key( ck );
      if( !_block_drop_eligible( k.block_id ) ) ineligible.push_back( ck );
    }
    if( !ineligible.empty() ){
      bool const trace = ( options.INTERFACE.IMPOSITION == Options::IC_TRACE );
      for( auto const& ck : ineligible ){
        if( trace ) dep.erase( ck );                       // ordinary trace constraint (no drop)
        else if( _strong_protected_claims.insert( ck ).second )
          gate_newly_protected = true;                     // newly kept-explicit -> needs a re-derive pass
      }
      if( _solve_display_level() >= 1 )
      _disp( 2 ) << "OCFESLV::setup ** continuity redundancy: "
                << ( trace ? "IC_TRACE" : "IC_STRONG" )
                << " block-eligibility gate: "
                << ( trace ? "vetoed " : "protected (keep-explicit) " )
                << ineligible.size() << " of " << n_raw
                << " column-dependent flag(s) in block(s) where "
                   "parabolic_structure_detected=0 (the sole condition in "
                   "_block_drop_eligible) at setup; "
                << ( trace ? "kept as ordinary continuity (no drop)."
                           : "routed to explicit tau (setup-internal two-pass; no driver re-derive)." )
                << std::endl;
    }
  }
#endif

  // -------------------------------------------------------------------------
  // IC_TRACE all_link keep (mirror of the IC_STRONG explicit-tau hybrid, 9210-
  // 9227).  IC_STRONG keeps any flagged continuity claim whose EVERY replacement-
  // candidate row is a LINK row ("all_link"): exact-replacing it would consume
  // the auxiliary's own defining equation, so STRONG retains it as an explicit
  // tau.  That branch is IC_STRONG-only, so IC_TRACE was DROPPING the same claims
  // -> weak-SAT fallback -> a fixed ~2e-2 spread at the t=LB initial face (where
  // the corner fix correctly suppresses U's value continuity, leaving the aux's
  // flux continuity the only tie; the aux has no IC).  STRONG reaching machine on
  // these proves they are load-bearing, not redundant.  Give IC_TRACE the same
  // keep: un-flag all_link claims so they stay ORDINARY trace constraints -- the
  // saddle absorbs the benign redundancy and nVar is fixed at setup (plan
  // invariant 2.5), so no explicit-tau machinery is needed.
  //
  // The candidate-row set is draft.receiver_edges (weak_sat_allowed) -- the SAME
  // source STRONG's replacement_candidates is built from (8962) -- and slots use
  // the same trace_slot_key canonicalization, so the kept set matches STRONG by
  // construction.  all_link holds ONLY where every replacement is a LINK row (the
  // t=LB aux case); interior aux claims keep a non-LINK interior-PDE replacement
  // (not all_link -> inert, e.g. PDE2) and LINK-free blocks have none (e.g. PDE7).
  if( options.INTERFACE.IMPOSITION == Options::IC_TRACE && !dep.empty() ){
    size_t const invalid_claim_id = std::numeric_limits<size_t>::max();
    auto slot_of = [&]( t_InterfaceClaimKey key ) -> t_InterfaceClaimKey {
      key.block_id = _trace_claim_block( key.state_id, key.dom_id );
      key.eqn_id   = invalid_claim_id;
      if( key.kind == CLAIM_EXACT_ROW ) key.kind = CLAIM_EXACT_C0;
      return key;
    };
    std::map<t_InterfaceClaimKey,bool> slot_all_link;   // slot -> every cand row is LINK
    for( auto const& edge : draft.receiver_edges ){
      if( !edge.weak_sat_allowed ) continue;
      bool const is_link = ( edge.row_id < _mEqn.size()
                          && _mEqn[ edge.row_id ].opt->role == EqnRole::LINK );
      t_InterfaceClaimKey const slot = slot_of( edge.claim_key );
      auto it = slot_all_link.find( slot );
      if( it == slot_all_link.end() ) slot_all_link.emplace( slot, is_link );
      else it->second = it->second && is_link;
    }
    std::vector<t_FrozenInterfaceClaimKey> keep;
    for( auto const& ck : dep ){
      t_InterfaceClaimKey const slot = slot_of( _unfreeze_interface_claim_key( ck ) );
      auto it = slot_all_link.find( slot );
      if( it != slot_all_link.end() && it->second ) keep.push_back( ck );
    }
    for( auto const& ck : keep ) dep.erase( ck );
    if( !keep.empty() && _solve_display_level() >= 1 )
      _disp( 2 ) << "OCFESLV::setup ** continuity redundancy: IC_TRACE "
                   "all_link keep: retained " << keep.size()
                << " LINK-anchored continuity claim(s) as ordinary trace "
                   "constraints (mirror of IC_STRONG hybrid; not dropped)."
                << std::endl;
  }

#ifdef MC__OCFESLV_SYMBOL_TAG_PROBE
  // -------------------------------------------------------------------------
  // Item 11 measure-before-implement (symbol-tag probe).  For every tau-coupled
  // continuity claim, print its principal-symbol / structural fingerprint and
  // whether the column-dependence proxy flagged it (`dep`).  Goal: confirm a
  // GENERATION-time symbol rule reproduces the verify-sound drop set (suppress on
  // PDE6, keep on PDE7) without the prune+verify+re-derive loop.  Diagnostic only;
  // the live flagging below is unchanged.  Header-only, gated, default-off.
  {
    auto state_ffvar = [&]( size_t sid ) -> FFVar const* {
      for( auto const& kv : _mVar )
        if( (size_t)kv.first.id().second == sid ) return &kv.first;
      return nullptr;
    };
    auto dom_ffvar = [&]( size_t did ) -> FFVar const* {
      for( auto const& kv : _mDom )
        if( (size_t)kv.first.id().second == did ) return &kv.first;
      return nullptr;
    };
    auto real_chars = [&]( arma::cx_vec const& ev, size_t& pos, size_t& neg ) -> size_t {
      size_t nr = 0; pos = 0; neg = 0;
      for( arma::uword i = 0; i < ev.n_elem; ++i ){
        double const re = ev(i).real(), im = ev(i).imag();
        if( std::fabs( im ) <= 1e-9 * std::max( 1.0, std::fabs( re ) ) ){
          ++nr; if( re > 1e-12 ) ++pos; else if( re < -1e-12 ) ++neg;
        }
      }
      return nr;
    };

    _disp( 2 ) << "OCFESLV::_symbol_tag_probe ** "
              << ( options.INTERFACE.IMPOSITION == Options::IC_TRACE ? "IC_TRACE" : "IC_STRONG" )
              << " : " << col2claim.size() << " tau-coupled claim(s), "
              << dep.size() << " column-dependent (flagged):\n";

    // Per-block principal-symbol fingerprint (the candidate item-11 rule input).
    for( auto const& bc : _blockClassification ){
      t_Classify const& cls = bc.second;
      _disp( 2 ) << "  block " << bc.first << ": type=" << pde_type_name( cls.type )
                << " analysed=" << cls.symbol_analysed
                << " rect=" << cls.symbol_rectangular
                << " parab_struct=" << cls.parabolic_structure_detected
                << " evol_hyp=" << cls.evolution_hyperbolic
                << " spat_char=" << cls.spatially_characteristic
                << " weak_hyp=" << cls.weak_hyperbolic
                << " At_sing=" << cls.At_singular
                << " degen=" << cls.degenerate << "\n";
      for( auto const& ed : cls.eigendata ){
        size_t pos = 0, neg = 0; size_t const nr = real_chars( ed.second, pos, neg );
        _disp( 2 ) << "      dir(";
        for( size_t k = 0; k < ed.first.size(); ++k )
          _disp( 2 ) << ( k ? "," : "" ) << ed.first[k];
        _disp( 2 ) << "): eig=" << ed.second.n_elem << " real=" << nr
                  << " in/out=" << pos << "/" << neg << "\n";
      }
    }

    // Per-claim structural fingerprint: LINK-auxiliary (flux) vs primitive, and
    // whether the block has an aux LINK in this claim's direction (the parabolic
    // over-specification signal).  Cross-reference btype + the dispatch echo for
    // the resolved interface type (IC_VALUE / IC_UPWIND / IC_FLUX).
    for( auto const& kv : col2claim ){
      t_InterfaceClaimKey const ck = _unfreeze_interface_claim_key( kv.second );
      bool const flagged = dep.count( kv.second ) > 0;
      FFVar const* sv = state_ffvar( ck.state_id );
      FFVar const* dv = dom_ffvar( ck.dom_id );
      bool const is_aux  = ( sv ? _is_auxiliary_state( *sv ) : false );
      bool const aux_dir = ( sv && dv ? _auxiliary_state_in_direction( *sv, *dv ) : false );
      bool const blk_aux = ( dv ? _block_has_aux_link_in_direction( ck.block_id, *dv ) : false );
      auto itc = _blockClassification.find( ck.block_id );
      char const* bt = ( itc != _blockClassification.end() )
                       ? pde_type_name( itc->second.type ) : "UNCLASSIFIED";
      _disp( 2 ) << "    claim(state=" << ck.state_id << " dom=" << ck.dom_id
                << " face=" << ck.face << " kind=" << ck.kind << "): block=" << ck.block_id
                << " btype=" << bt << " aux_state=" << is_aux << " aux_in_dir=" << aux_dir
                << " block_link_dir=" << blk_aux
                << " -> " << ( flagged ? "FLAGGED(drop-cand)" : "kept" ) << "\n";
    }
  }
#endif

#ifndef MC__OCFESLV_VALUE_SLAVED_KEEP
  // Value-slaved drop.  An algebraic state defined by a derivative-free algebraic
  // constraint (P - Rg c = 0) is pinned pointwise at every node including
  // interfaces, so it never needs a C0-continuity row: the claim is redundant
  // when the source states are continuous and harmless otherwise (the state is
  // closed by its constraint either way).  The receiver-coupling proxy above is
  // structurally blind to this dependence (P's claim and c's claim sit on disjoint
  // DOFs; only the EOS links them, and the EOS is not part of the tau coupling the
  // proxy inspects), so flag those claims here.  GATED on _block_drop_eligible:
  // the algebraic implication is real, but in a symbol-INELIGIBLE/degenerate block
  // the corner/tau assignment is contended and dropping the (redundant) value
  // claim disturbs a sibling state's C0 (measured: PDE2/3 Vg interface jump under
  // IC_TRACE/IC_STRONG -- confirmed by an -DMC__OCFESLV_VALUE_SLAVED_KEEP A/B).  In
  // such blocks the claim is LEFT as ordinary continuity (not inserted into the
  // drop set), which the standard IC_TRACE/IC_STRONG machinery resolves to machine
  // precision; no explicit protect is needed because these claims were never
  // proxy-flagged into dep.  In ELIGIBLE (parabolic-PSA) blocks the drop fires
  // exactly as before -- PDE30 rung 6 (EOS P=Rg c) is unaffected.  Mirror of the
  // boundary-closure rule.  Compile -DMC__OCFESLV_VALUE_SLAVED_KEEP to disable (A/B).
  {
    std::set<size_t> const vs = _value_slaved_algebraic_state_ids();
    size_t vs_added = 0;
    if( !vs.empty() )
      for( auto const& kv : col2claim ){
        t_InterfaceClaimKey const k = _unfreeze_interface_claim_key( kv.second );
        if( vs.count( k.state_id ) && _block_drop_eligible( k.block_id )
            && dep.insert( kv.second ).second ) ++vs_added;
      }
    if( vs_added && _solve_display_level() >= 1 )
      _disp( 2 ) << "OCFESLV::setup ** continuity redundancy: "
                << ( options.INTERFACE.IMPOSITION == Options::IC_TRACE ? "IC_TRACE" : "IC_STRONG" )
                << " value-slaved drop: flagged " << vs_added << " implied continuity "
                   "claim(s) on " << vs.size() << " value-slaved algebraic state(s) "
                   "(derivative-free defining constraint)." << std::endl;
  }
#endif

  if( dep.empty() ) return true;

  // Respect the protected set (item 10 step 3, DROP_REDERIVE; and the IC_STRONG
  // block-eligibility gate above): protected claims are kept explicit, never
  // dropped.  Track whether any NEW claim is dropped THIS pass -- on a pass-2
  // re-run the already-dropped claims are excluded upstream and the already-
  // protected ones are kept explicit, so "no new drop" means the plan is stable
  // (this is what lets a mixed drop+keep-explicit block set settle in two passes).
  // C1 (rev220): a claim on a state that C0 tags NATURAL in the claim's direction OWNS its condition and has a
  // natural receiver; it is never a drop candidate.  Without this, removing the implied claims first makes the
  // natural ones the "redundant" set by rank, and the rank test drops them (PDE20f: Dz_u, accuracy 10x worse).
  size_t n_c0_protected = 0;
  if( _knob_CLAIM_FROM_C0() ){
    auto const cm = claim_classification();
    for( auto const& ck : dep ){
      FFVar const* sv = _state_by_id( std::get<1>( ck ) );
      FFVar const* dv = _domain_by_id( std::get<3>( ck ) );
      if( !sv || !dv ) continue;
      auto is = cm.find( *sv ); if( is == cm.end() ) continue;
      auto id = is->second.find( *dv ); if( id == is->second.end() ) continue;
      if( id->second.tag == ClaimTag::NATURAL && _strong_protected_claims.insert( ck ).second ) ++n_c0_protected;
    }
    if( n_c0_protected )
      _disp( 2 ) << "OCFESLV::setup ** claim[C0]: " << n_c0_protected
                << " NATURAL claim(s) protected from the W-decide drop\n";
  }
  // C2b (rev279): a flagged claim whose rescued receivers were SILENCED under EXACT_NATURAL is not redundant --
  // its column is dependent only because the duplicate that made it independent is silent (blk0 MANUAL: the 16
  // corner claims sharing a PDE row).  Restore its receivers on the re-run instead of dropping the claim.
  size_t n_restored = 0;
  if( _knob_EXACT_NATURAL() ){
    // rev285 (R2): a dependency is between TWO columns and the pass flags the later one; the silenced partner
    // may be the one that stays (bisect2 chain 3).  So restore (a) a flagged claim whose own terms were silenced
    // and (b) every silenced claim that shares a receiver ROW with a flagged claim -- the row the dependency
    // lives in.  Both are "column-dependent with the natural receiver alone".
    // rev285 (R2): keys in the TRACE-SLOT canonical form, as `dep` and the silencing itself use (rev279 froze
    // the raw term key, which matched blk0's corner claims by luck and none of bisect2's).
    std::set<t_FrozenInterfaceClaimKey> silenced;
    for( auto const& term : draft.weak_sat_terms )
      if( term.multiplier_silenced ) silenced.insert( _trace_slot_key( term.claim_key ) );
    // (b) a flagged claim whose dependent PARTNER is the silenced one: they meet in a COLLOCATED receiver row
    // (equation id AND flat node index -- the equation id alone is shared corpus-wide and restores everything).
    typedef std::pair<size_t,size_t> RowAt;
    std::map<RowAt, std::set<t_FrozenInterfaceClaimKey>> silenced_at;
    std::map<t_FrozenInterfaceClaimKey, std::set<RowAt>> rows_of;
    for( auto const& term : draft.weak_sat_terms ){
      auto const ck = _trace_slot_key( term.claim_key );
      RowAt const at( term.row_id, term.receiver_flat );
      rows_of[ ck ].insert( at );
      if( term.multiplier_silenced ) silenced_at[ at ].insert( ck );
    }
    std::set<t_FrozenInterfaceClaimKey> to_restore;
    size_t n_group = 0;
    for( auto const& ck : dep ){
      if( silenced.count( ck ) ) to_restore.insert( ck );
      // rev286 (R2b): the silenced members of the column's dependency group.
      auto itg = dep_group.find( ck );
      if( itg != dep_group.end() )
        for( size_t qi : itg->second )
          if( qi < basis_claim.size() && silenced.count( basis_claim[qi] ) && to_restore.insert( basis_claim[qi] ).second ) ++n_group;
      auto itr = rows_of.find( ck );
      if( itr == rows_of.end() ) continue;
      for( auto const& at : itr->second ){
        auto its = silenced_at.find( at );
        if( its != silenced_at.end() ) to_restore.insert( its->second.begin(), its->second.end() );
      }
    }
    for( auto const& ck : to_restore )
      if( !_strong_dropped_claims.count( ck ) && _exact_natural_restored.insert( ck ).second ) ++n_restored;
    if( n_restored ){
      _strong_restore_pending = true;      // a RESTORE, not a drop: its own retry
      _disp( 2 ) << "OCFESLV::setup ** continuity redundancy: EXACT_NATURAL -- " << n_restored
                << " claim(s) restored to their rescued receivers (" << n_group << " as dependency-group partners; re-deriving plan; no drop)" << std::endl;
      return false;   // pass-2 re-run with the restored set
    }
  }
  size_t n_protected = 0;
  bool   any_new_drop = false;
  for( auto const& ck : dep ){
    if( _strong_protected_claims.count( ck ) ){ ++n_protected; continue; }
    if( _exact_natural_restored.count( ck ) ){ ++n_protected; continue; }   // rev279: restored, not droppable
    if( _strong_dropped_claims.insert( ck ).second ){ any_new_drop = true;
      if( std::getenv( "CRONOS_PLAN_DUMP_DROPPED" ) )
        _disp( 3 ) << "  [plan/dropped] block=" << std::get<0>( ck ) << " state=" << std::get<1>( ck )
                  << " dom=" << std::get<3>( ck ) << " iel_lo=" << std::get<4>( ck )
                  << " face=" << std::get<5>( ck ) << " kind=" << std::get<6>( ck ) << std::endl;
    }
  }
  if( !any_new_drop ){
    if( gate_newly_protected ){
      // The IC_STRONG gate just protected claim(s) that _store_frozen_interface_plan
      // (step 1) already routed to the Schur corner this pass.  Re-run the pipeline
      // once (two-pass) so step 1 keeps them explicit; nothing is dropped.  On the
      // re-run the gate finds them already protected (no new protection) and returns
      // true, so the loop terminates after the second pass.
      _strong_redundancy_pending = true;
      if( _solve_display_level() >= 1 )
      _disp( 2 ) << "OCFESLV::setup ** continuity redundancy: IC_STRONG"
                << " block-eligibility gate: re-deriving plan to keep " << n_protected
                << " symbol-ineligible continuity claim(s) explicit (no drop)." << std::endl;
      return false;   // request the pass-2 re-run that keeps these claims explicit.
    }
    if( n_protected && _solve_display_level() >= 1 )
      _disp( 2 ) << "OCFESLV::setup ** continuity redundancy: "
                << ( options.INTERFACE.IMPOSITION == Options::IC_TRACE ? "IC_TRACE" : "IC_STRONG" )
                << " redundancy pre-pass: all " << n_protected
                << " flagged claim(s) are protected/kept-explicit; no new drop." << std::endl;
    return true;   // nothing new to exclude -> proceed normally
  }
  _strong_redundancy_pending = true;
  if( _solve_display_level() >= 1 ){
  _disp( 2 ) << "OCFESLV::setup ** continuity redundancy: "
            << ( options.INTERFACE.IMPOSITION == Options::IC_TRACE ? "IC_TRACE" : "IC_STRONG" )
            << " redundancy pre-pass: flagged " << _strong_dropped_claims.size()
            << " column-dependent (redundant/implied) continuity claim(s) for "
               "pass-2 exclusion (of " << col.size() << " tau-coupled columns";
  if( n_protected ) _disp( 2 ) << ", " << n_protected << " protected/kept";
  _disp( 2 ) << ")." << std::endl;
  }
  return false;   // request the pass-2 re-run that excludes these claims.
}

// Adjoint of the capture quadrature (Stage 1b): scatter lam_cap[R] onto ls at R's collocation dofs with
// weights jac*w_ref[m], summed over evolution elements -- transpose of _fill_capture_sens's ACCUM sum.
// Modelled on terminal_profile_adjoint's coordinate->dof stride map, but with the LGR quadrature weights.
inline void
OCFEBASE::_fill_capture_sens_adjoint( std::map<FFVar,std::vector<double>,lt_FFVar> const& lam_cap,
                                   double* ls ) const
{
  if( !ls || !_evolution_dom_set || !_evolution_dom_var.dag() ) return;
  auto const evo_id = _evolution_dom_var.id();
  auto ievo = _mDom.find( _evolution_dom_var );
  if( ievo == _mDom.end() ) return;
  FFDom const&              ed   = ievo->second;
  std::vector<double> const wref = ed.lgweights();

  for( auto const& [state, lamvec] : lam_cap ){        // state = R (capture source)
    if( lamvec.empty() ) continue;
    auto itv = _mVar.find( state );
    if( itv == _mVar.end() ) continue;
    auto const& var_dom = itv->second;                 // dep order (== eval_colloc dims order)

    size_t const D = var_dom.size();
    std::vector<size_t> nnode( D, 0 ), nStride( D, 0 );
    size_t evoDim = D;
    { size_t s = 1, d = 0;
      for( auto const& v : var_dom ){ size_t const nn = _mDom.at(v).n_node; nnode[d]=nn; nStride[d]=s; s*=nn; if(v.id()==evo_id) evoDim=d; ++d; } }
    if( evoDim == D ) continue;                        // R not on the evolution domain

    std::vector<FFVar> sVar; std::vector<size_t> sPos, sElem, sNode;
    { size_t d=0; for( auto const& v : var_dom ){ if(v.id()!=evo_id){ FFDom const& dd=_mDom.at(v); sVar.push_back(v); sPos.push_back(d); sElem.push_back(dd.n_elem); sNode.push_back(dd.n_node);} ++d; } }
    size_t total = 1; for( size_t k=0;k<sVar.size();++k) total *= sElem[k]*sNode[k];
    if( lamvec.size() < total ) continue;

    size_t const nme = std::min( wref.size(), nnode[evoDim] );
    for( size_t lin = 0; lin < total; ++lin ){
      double const lj = lamvec[lin];
      if( lj == 0. ) continue;
      std::map<FFVar,size_t,lt_FFVar> ndx_el;
      size_t node_base = 0, rem = lin;
      for( size_t k=0;k<sVar.size();++k){ size_t const gtot=sElem[k]*sNode[k]; size_t const g=rem%gtot; rem/=gtot; size_t const ie=g/sNode[k], in=g%sNode[k]; ndx_el[sVar[k]]=ie; node_base += in*nStride[sPos[k]]; }
      for( size_t ie = 0; ie < ed.n_elem; ++ie ){       // sum over evolution elements (1 under marching)
        double const jac = 0.5 * ( ed.elem_bnd[ie+1] - ed.elem_bnd[ie] );
        ndx_el[ _evolution_dom_var ] = ie;
        size_t const inc = pos_state( state, ndx_el );
        for( size_t m = 0; m < nme; ++m )
          ls[ inc + node_base + m * nStride[evoDim] ] += jac * wref[m] * lj;
      }
    }
  }
}

// LATCH variant of _fill_capture_sens_adjoint: the exact transpose of the forward latch
// d(cap)/dp = eval_colloc(src,tau).  Scatters ONE latch capture's cap-adjoint lamC onto ls[src] at
// TAU's evolution element with the value Lagrange weights L_m(tau) -- NOT the ACCUM quadrature
// jac*w_ref summed over all elements.  Called per-capture (so tau is not lost to source-aggregation).
// Under marching the evolution domain is the current window, so this no-ops unless tau lies in that
// window (the transpose of the forward cap_window==k gate); under a monolithic solve the full-horizon
// domain locates tau's element directly.
inline void
OCFEBASE::_fill_capture_sens_adjoint_latch( FFVar const& src, double tau,
                                         std::vector<double> const& lamC, double* ls, int const side ) const
{
  if( !ls || lamC.empty() || !_evolution_dom_set || !_evolution_dom_var.dag() ) return;
  auto const evo_id = _evolution_dom_var.id();
  auto ievo = _mDom.find( _evolution_dom_var );
  if( ievo == _mDom.end() ) return;
  FFDom const& ed = ievo->second;

  // Locate tau's evolution element and its reference coordinate xi in [-1,1] with the FORWARD's own rule
  // (_point_in_domain, the one eval_colloc uses): on an element boundary it takes the element STARTING there.  This
  // loop used to take the FIRST element containing tau -- the one ENDING there -- so for a tau on a boundary the
  // adjoint read the other side than the forward, off by the (solver-precision) continuity gap (OCFE_EVALSUM, ~1e-6).
  int ie_tau = -1; double xi = 0.;
  {
    double const eps = 1e-9 * std::max( 1.0, std::fabs( ed.up_dom - ed.lo_dom ) );
    if( tau >= ed.lo_dom - eps && tau <= ed.up_dom + eps ){
      auto const pr = _point_in_domain( ed, tau, side );
      ie_tau = (int)pr.first;  xi = pr.second;
    }
  }
  if( ie_tau < 0 ) return;                              // tau not in this (windowed) evolution domain
  std::vector<double> Lw;
  if( !ed.w_eval_lagrange( Lw, xi, 0 ) ) return;        // L_m(tau), value order 0

  auto itv = _mVar.find( src );
  if( itv == _mVar.end() ) return;
  auto const& var_dom = itv->second;
  size_t const D = var_dom.size();
  std::vector<size_t> nnode( D, 0 ), nStride( D, 0 );
  size_t evoDim = D;
  { size_t s = 1, d = 0;
    for( auto const& v : var_dom ){ size_t const nn = _mDom.at(v).n_node; nnode[d]=nn; nStride[d]=s; s*=nn; if(v.id()==evo_id) evoDim=d; ++d; } }
  if( evoDim == D ) return;

  std::vector<FFVar> sVar; std::vector<size_t> sPos, sElem, sNode;
  { size_t d=0; for( auto const& v : var_dom ){ if(v.id()!=evo_id){ FFDom const& dd=_mDom.at(v); sVar.push_back(v); sPos.push_back(d); sElem.push_back(dd.n_elem); sNode.push_back(dd.n_node);} ++d; } }
  size_t total = 1; for( size_t k=0;k<sVar.size();++k) total *= sElem[k]*sNode[k];
  if( lamC.size() < total ) return;

  size_t const nme = std::min( Lw.size(), nnode[evoDim] );
  for( size_t lin = 0; lin < total; ++lin ){
    double const lj = lamC[lin];
    if( lj == 0. ) continue;
    std::map<FFVar,size_t,lt_FFVar> ndx_el;
    size_t node_base = 0, rem = lin;
    for( size_t k=0;k<sVar.size();++k){ size_t const gtot=sElem[k]*sNode[k]; size_t const g=rem%gtot; rem/=gtot; size_t const ie=g/sNode[k], in=g%sNode[k]; ndx_el[sVar[k]]=ie; node_base += in*nStride[sPos[k]]; }
    ndx_el[ _evolution_dom_var ] = (size_t)ie_tau;
    size_t const inc = pos_state( src, ndx_el );
    for( size_t m = 0; m < nme; ++m )
      ls[ inc + node_base + m * nStride[evoDim] ] += Lw[m] * lj;
  }
}

inline bool

//==========================================================================================================
// INTERFACE PLAN -- claims -> realisations -> rows and columns.  The DECISIONS live in ocplan.hpp; what is here
// is the assembly: which claims exist, which rows receive them, and how the chosen realisation becomes
// matrix structure.
// (banner placed from the definitions themselves: this region is lines 20165-24276 of the pre-rev268 file)
//==========================================================================================================

OCFEBASE::_build_frozen_strong_tau_elimination
( InterfacePlanDraft& draft ) const
{
  draft.strong_tau_elim.clear();
  if( options.INTERFACE.IMPOSITION != Options::IC_STRONG )
    return true;

  // The Schur projection operates on the accepted tau basis only
  // (n_trace_var), not the full raw claim set (n_trace_constraint).
  // ...and excludes the IC_STRONG-hybrid explicit (kept) taus, which are real
  // var[] columns assembled IC_TRACE-style rather than eliminated.  Their
  // lambda_index lives at >= trace_var_offset, well above the eliminated block.
  size_t const ntrace = draft.n_trace_var - draft.n_trace_var_explicit;
  draft.strong_tau_elim.n_trace = ntrace;
  if( !ntrace ){
    draft.strong_tau_elim.ready = true;
    return true;
  }
  if( !draft.weak_sat_terms_ready ){
    _disp( 2 ) << "OCFESLV::setup ** tau elimination: frozen "
                 "receiver terms are not ready" << std::endl;
    return false;
  }

  auto make_row_key = []( FrozenWeakSatTerm const& t )
  {
    FrozenPhysicalResidualRowKey r;
    r.row_id   = t.row_id;
    r.block_el = t.receiver_block_el;
    r.flat     = t.receiver_emit_flat;
    return r;
  };

  std::map<FrozenPhysicalResidualRowKey,std::vector<double>> Brow;
  for( auto const& term : draft.weak_sat_terms ){
    if( !term.trace_tau_allowed ) continue;
    // Explicit (kept) taus sit at var[] columns >= trace_var_offset; they are
    // not part of the eliminated Schur block, so skip their receiver terms.
    if( draft.n_trace_var_explicit
     && term.lambda_index >= draft.trace_var_offset ) continue;
    if( term.lambda_index == std::numeric_limits<size_t>::max()
     || term.lambda_index >= ntrace ){
      std::cerr << "OCFESLV::setup ** tau elimination: invalid "
                   "lambda index in frozen receiver term: row=" << term.row_id
                << " lambda_index=" << term.lambda_index
                << " ntrace=" << ntrace << std::endl;
      return false;
    }
    if( std::fabs( term.trace_prefactor ) <= 64. * DBL_EPSILON ) continue;
    auto& row = Brow[ make_row_key( term ) ];
    if( row.empty() ) row.assign( ntrace, 0.0 );
    row[ term.lambda_index ] += term.trace_prefactor;
  }

  if( Brow.empty() ){
    _disp( 2 ) << "OCFESLV::setup ** tau elimination: tau "
                 "receiver block is empty for " << ntrace
              << " exact trace claim(s)" << std::endl;
    return false;
  }

#ifdef MC__OCFESLV_STRONG_DUMMY_DIAG
  // Print-only diagnostic (no behaviour change; the function still aborts below
  // when the block is rank-deficient).  Reports the eliminable / kept-explicit
  // partition the dummy-derivative hybrid will need:
  //   - rank(B over ALL rows)        = the true column rank of the tau block;
  //       ntrace - rank_all taus are GENUINELY column-dependent (parallel
  //       aux/primitive clusters) and MUST stay explicit (dummy derivatives).
  //   - rank(B over VALID rows only) = what the Schur pivot selection can reach
  //       (valid = not an exact-replacement / dropped row).  The extra shortfall
  //       rank_all - rank_valid is "eliminable in principle, but its only
  //       independent pivot lives on a dropped row" -- a pivot-availability
  //       conflict, a different fix from a true dummy.
  // Column order is lambda_index, so the greedy column elimination is
  // deterministic and mirrors the row pivot selection used below.
  {
    struct TauDesc { size_t state_id=0, dom_id=0, iel_lo=0, face=0; int kind=0;
                     OCPlan::InterfaceEdgeSource source =
                       OCPlan::InterfaceEdgeSource::PRINCIPAL_SYMBOL;
                     bool have=false; };
    std::vector<TauDesc> desc( ntrace );
    for( auto const& term : draft.weak_sat_terms ){
      if( !term.trace_tau_allowed ) continue;
      size_t const li = term.lambda_index;
      if( li >= ntrace || desc[li].have ) continue;
      desc[li] = { term.claim_key.state_id, term.claim_key.dom_id,
                   term.claim_key.iel_lo, term.claim_key.face,
                   term.claim_key.kind, term.source, true };
    }
    auto src_str = []( OCPlan::InterfaceEdgeSource s ) -> char const* {
      switch( s ){
        case OCPlan::InterfaceEdgeSource::PRINCIPAL_SYMBOL: return "PRINCIPAL";
        case OCPlan::InterfaceEdgeSource::AUXILIARY_ALGEBRAIC: return "AUX";
        case OCPlan::InterfaceEdgeSource::GEOMETRIC_C0:     return "GEOM_C0";
      }
      return "?";
    };
    // Build B as (rows x ntrace) for a chosen row set, then return the dependent
    // column indices via greedy column Gram-Schmidt (same scheme as the row
    // selection below, transposed).
    auto col_dependent = [&]( bool valid_only, size_t& rank_out )
                         -> std::vector<size_t> {
      std::vector<std::vector<double> const*> rows;
      for( auto const& kv : Brow ){
        if( valid_only
         && draft.exact_replacement_rows.find( kv.first )
            != draft.exact_replacement_rows.end() ) continue;
        rows.push_back( &kv.second );
      }
      size_t const m = rows.size();
      std::vector<arma::vec> qbasis_c;
      std::vector<size_t> dep;
      rank_out = 0;
      for( size_t j = 0; j < ntrace; ++j ){
        arma::vec c( (arma::uword)m );
        for( size_t i = 0; i < m; ++i ) c( (arma::uword)i ) = (*rows[i])[j];
        double const n0 = arma::norm( c );
        if( n0 <= 64. * DBL_EPSILON ){ dep.push_back( j ); continue; }
        for( auto const& q : qbasis_c ) c -= arma::dot( c, q ) * q;
        double const nr = arma::norm( c );
        if( nr <= 1e-9 * std::max( 1.0, n0 ) ){ dep.push_back( j ); continue; }
        c /= nr; qbasis_c.push_back( std::move( c ) ); ++rank_out;
      }
      return dep;
    };
    size_t rank_all = 0, rank_valid = 0;
    std::vector<size_t> dep_all   = col_dependent( false, rank_all );
    std::vector<size_t> dep_valid = col_dependent( true,  rank_valid );
    std::set<size_t> dep_all_set( dep_all.begin(), dep_all.end() );
    std::vector<size_t> dep_valid_only;
    for( size_t j : dep_valid )
      if( !dep_all_set.count( j ) ) dep_valid_only.push_back( j );

    auto dump = [&]( char const* label, std::vector<size_t> const& cols ){
      _disp( 3 ) << "  " << label << " (" << cols.size() << "):\n";
      for( size_t j : cols ){
        _disp( 3 ) << "    tau[" << j << "] ";
        if( desc[j].have )
          _disp( 3 ) << "state=" << desc[j].state_id
                    << " dom=" << desc[j].dom_id
                    << " iel_lo=" << desc[j].iel_lo
                    << " face=" << desc[j].face
                    << " kind=" << desc[j].kind
                    << " src=" << src_str( desc[j].source );
        else _disp( 3 ) << "(no descriptor)";
        _disp( 3 ) << "\n";
      }
    };
    _disp( 3 ) << "\n#### IC_STRONG dummy-derivative partition ####\n";
    _disp( 3 ) << "  ntrace=" << ntrace
              << "  rank(all rows)=" << rank_all
              << "  rank(valid rows)=" << rank_valid
              << "  exact_replacement_rows=" << draft.exact_replacement_rows.size()
              << "\n";
    _disp( 3 ) << "  GENUINE column-dependent (parallel clusters; must be dummy)"
              << " = " << dep_all.size() << "\n";
    dump( "column-dependent taus", dep_all );
    _disp( 3 ) << "  PIVOT-CONFLICTED (eliminable over all rows, pivot only on a"
                 " dropped row) = " << dep_valid_only.size() << "\n";
    dump( "pivot-conflicted taus", dep_valid_only );
    _disp( 3 ) << "  => hybrid keep-explicit lower bound = column-dependent"
                 " (" << dep_all.size() << ");"
                 " pivot-conflicted are recoverable by freeing the dropped row"
                 " or also kept explicit.\n";
    _disp( 3 ) << "##############################################\n";

    // ---- Stage-0 dummy-derivative safety probe -------------------------------
    // For each GENUINE column-dependent tau, characterise the redundancy so a
    // clean, locally-resolvable alias (where a dummy-derivative is safe) can be
    // told apart from an entangled one.  Per dependent column we report: the tau
    // it is most collinear with (its parallel partner), and every receiver row
    // its coupling touches together with that row's sharing degree (how many tau
    // columns are nonzero there).  A maximum sharing degree of 2 means an
    // isolated redundant pair -- one member can be made a dummy / kept explicit
    // with no collateral.  A higher degree means the redundancy is entangled
    // with other claims, so a receiver-consuming fix would silently break those
    // other claims (the in-place hazard we are avoiding), and the corner must be
    // resolved via the corner-DOF's own interior row or a keep-explicit
    // fallback.  Pure diagnostic; no plan state is changed.
    if( !dep_all.empty() ){
      std::vector<std::vector<double> const*> rp;
      std::vector<FrozenPhysicalResidualRowKey> rk;
      rp.reserve( Brow.size() ); rk.reserve( Brow.size() );
      for( auto const& kv : Brow ){ rp.push_back( &kv.second ); rk.push_back( kv.first ); }
      size_t const M = rp.size();
      double const tol = 1e-12;
      std::vector<size_t> row_degree( M, 0 );
      for( size_t i = 0; i < M; ++i ){
        size_t d = 0;
        for( size_t jj = 0; jj < ntrace; ++jj )
          if( std::abs( (*rp[i])[jj] ) > tol ) ++d;
        row_degree[i] = d;
      }
      _disp( 3 ) << "\n#### IC_STRONG dummy-derivative stage-0 safety probe ####\n";
      for( size_t j : dep_all ){
        double nj = 0.0;
        for( size_t i = 0; i < M; ++i ){ double const v = (*rp[i])[j]; nj += v*v; }
        nj = std::sqrt( nj );
        double best = 0.0; size_t bestk = ntrace;
        if( nj > 0.0 )
          for( size_t k = 0; k < ntrace; ++k ){
            if( k == j ) continue;
            double dot = 0.0, nk = 0.0;
            for( size_t i = 0; i < M; ++i ){
              double const a = (*rp[i])[j], b = (*rp[i])[k];
              dot += a*b; nk += b*b;
            }
            nk = std::sqrt( nk );
            if( nk > 0.0 ){
              double const c = std::abs( dot ) / ( nj*nk );
              if( c > best ){ best = c; bestk = k; }
            }
          }
        _disp( 3 ) << "  dependent tau[" << j << "] state=" << desc[j].state_id
                  << " dom=" << desc[j].dom_id << " iel_lo=" << desc[j].iel_lo
                  << " face=" << desc[j].face << " kind=" << desc[j].kind
                  << " src=" << src_str( desc[j].source ) << "\n";
        if( bestk < ntrace )
          _disp( 3 ) << "      parallel partner tau[" << bestk << "] |cos|="
                    << std::scientific << std::setprecision(4) << best
                    << "  state=" << desc[bestk].state_id
                    << " dom=" << desc[bestk].dom_id
                    << " iel_lo=" << desc[bestk].iel_lo
                    << " face=" << desc[bestk].face
                    << " kind=" << desc[bestk].kind
                    << " src=" << src_str( desc[bestk].source ) << "\n";
        size_t shared = 0, support = 0, maxdeg = 0;
        _disp( 3 ) << "      coupling support (per row: the other tau(s) sharing it):";
        for( size_t i = 0; i < M; ++i ){
          if( std::abs( (*rp[i])[j] ) <= tol ) continue;
          ++support; maxdeg = std::max( maxdeg, row_degree[i] );
          if( row_degree[i] > 1 ) ++shared;
          _disp( 3 ) << "\n        row[row_id=" << rk[i].row_id
                    << " flat=" << rk[i].flat << " deg=" << row_degree[i]
                    << "] shared with:";
          bool any = false;
          for( size_t k = 0; k < ntrace; ++k ){
            if( k == j ) continue;
            if( std::abs( (*rp[i])[k] ) <= tol ) continue;
            any = true;
            _disp( 3 ) << " tau[" << k << " state=" << desc[k].state_id
                      << " dom=" << desc[k].dom_id << " iel_lo=" << desc[k].iel_lo
                      << " face=" << desc[k].face << " kind=" << desc[k].kind << "]";
          }
          if( !any ) _disp( 3 ) << " (none -- private)";
        }
        _disp( 3 ) << "\n      => support=" << support << " shared(deg>1)=" << shared
                  << " max_deg=" << maxdeg << "  "
                  << ( maxdeg <= 2
                       ? "isolated redundant pair -- dummy/keep-explicit is local and safe"
                       : "entangled with other claims -- receiver-consuming fix unsafe; "
                         "use corner-DOF interior-row dummy or keep-explicit fallback" )
                  << "\n";
      }
      _disp( 3 ) << "########################################################\n";
    }
  }
#endif

  auto row_norm = []( std::vector<double> const& r ){
    double n = 0.0;
    for( double v : r ) n += v*v;
    return std::sqrt( n );
  };

  // Sparse Gram-Schmidt support for the tau-pivot basis.  B_S is ~0.1% dense
  // (each receiver row couples a tau to only a few others), so the dense GS in
  // try_accept_pivot/gs_extends over length-ntrace rows iterates almost all
  // zeros -- the IC_STRONG setup bottleneck (pivot 'select').  Keep the
  // orthonormal basis sparse (index-sorted); disjoint-support rows give an
  // exact-zero dot and are skipped, so 'r -= 0*q' (a dense no-op) is matched
  // bit-for-bit while the cost scales with real support overlap.
  typedef std::vector<std::pair<size_t,double>> SparseVec;
  auto to_sparse = []( std::vector<double> const& d ) -> SparseVec {
    SparseVec s;
    for( size_t j = 0; j < d.size(); ++j ) if( d[j] != 0.0 ) s.emplace_back( j, d[j] );
    return s;   // ascending index order by construction
  };
  auto sp_norm = []( SparseVec const& r ){
    double n = 0.0; for( auto const& e : r ) n += e.second*e.second; return std::sqrt( n );
  };
  auto sp_dot = []( SparseVec const& a, SparseVec const& b ) -> double {
    double d = 0.0; size_t i = 0, j = 0;
    while( i < a.size() && j < b.size() ){
      if( a[i].first < b[j].first ) ++i;
      else if( a[i].first > b[j].first ) ++j;
      else { d += a[i].second * b[j].second; ++i; ++j; }
    }
    return d;
  };
  auto sp_axpy = []( SparseVec& r, double s, SparseVec const& q ){
    SparseVec out; out.reserve( r.size() + q.size() );   // r := r - s*q
    size_t i = 0, j = 0;
    while( i < r.size() && j < q.size() ){
      if( r[i].first < q[j].first ){ out.push_back( r[i] ); ++i; }
      else if( r[i].first > q[j].first ){ out.emplace_back( q[j].first, -s*q[j].second ); ++j; }
      else { out.emplace_back( r[i].first, r[i].second - s*q[j].second ); ++i; ++j; }
    }
    for( ; i < r.size(); ++i ) out.push_back( r[i] );
    for( ; j < q.size(); ++j ) out.emplace_back( q[j].first, -s*q[j].second );
    r.swap( out );
  };

  // ---- Joint tau-pivot / geometric-C0 drop assignment (matroid-aware) -------
  // Pivots (tau eliminations) and drops (geometric-C0 exact continuity rows)
  // compete for the same redundant copies at corner nodes.  Choosing pivots
  // greedily first over-constrains that shared pool -- the rank-extending corner
  // rows get borrowed as pivots, leaving the conflicted claims with no free
  // sibling (the matched-102/114 failure).  Instead the two are assigned
  // jointly.  A unified consumer matching is provably feasible (matching probe,
  // 952/952), so a valid joint assignment exists; the construction below finds
  // one while guaranteeing the pivot block is non-singular.
  //
  // Maintain a geometric-C0 system of distinct representatives (claim -> drop
  // row) and a Gram-Schmidt pivot basis over the NON-drop rows.  While the pivot
  // rank is short, a rank-extending DROP row must exist (the missing directions
  // live only on dropped corner rows); free it by re-routing its claim along an
  // augmenting path that avoids the pivot rows, then take it as a pivot.
  // Re-routing only moves claims among non-pivot redundant copies, so the basis
  // is never disturbed and recovered pivots are never re-dropped; rank strictly
  // increases, terminating at ntrace.  Swapping which redundant copy a claim
  // drops is exact-imposition-neutral (continuity rows are keyed by claim, not
  // by physical row), so IC_TRACE is untouched and |exact_replacement_rows| is
  // preserved.  Gram-Schmidt keeps B_S non-singular without a numerical fallback.
  // A residual short-fall would be GENUINE column-dependence (parallel
  // aux/primitive clusters) needing the keep-explicit dummy-derivative hybrid.
  std::vector<SparseVec> qbasis;   // sparse orthonormal tau-pivot basis
  std::vector<std::vector<double>> Bsel;
  std::set<FrozenPhysicalResidualRowKey> pivot_set;

  // Accept row rk (with tau-receiver row brow) as a pivot iff it is numerically
  // independent of the rows accepted so far.  Returns true on acceptance.
  auto try_accept_pivot =
    [&]( FrozenPhysicalResidualRowKey const& rk,
         std::vector<double> const& brow ) -> bool {
      // a definition row a REDUNDANT decision relies on is never a pivot (PlanDecisions::reserved_definition_rows)
      if( draft.reserved_definition_rows.count( rk ) ) return false;
      SparseVec r = to_sparse( brow );
      double const n0 = sp_norm( r );
      if( n0 <= 64. * DBL_EPSILON ) return false;
      for( auto const& q : qbasis ){
        double const dot = sp_dot( r, q );
        if( dot == 0.0 ) continue;     // disjoint support -> dense 'r -= 0*q' no-op
        sp_axpy( r, dot, q );
      }
      double const nr = sp_norm( r );
      if( nr <= 1e-10 * std::max( 1.0, n0 ) ) return false;
      double const inv = 1.0 / nr;
      for( auto& e : r ) e.second *= inv;
      qbasis.push_back( std::move( r ) );
      draft.strong_tau_elim.pivot_index[ rk ] =
        draft.strong_tau_elim.pivot_rows.size();
      draft.strong_tau_elim.pivot_rows.push_back( rk );
      Bsel.push_back( brow );
      pivot_set.insert( rk );
      return true;
    };

  // Geometric-C0 claim candidates (redundant copies) and the candidate-row set.
  std::map<t_InterfaceClaimKey,std::vector<FrozenPhysicalResidualRowKey>> geom_cand;
  std::set<FrozenPhysicalResidualRowKey> geom_cand_rows;
  for( auto const& term : draft.weak_sat_terms ){
    if( term.source != OCPlan::InterfaceEdgeSource::GEOMETRIC_C0 ) continue;
    FrozenPhysicalResidualRowKey const rk = make_row_key( term );
    geom_cand[ term.claim_key ].push_back( rk );
    geom_cand_rows.insert( rk );
  }
  // Non-geometric exact replacements (e.g. parallel-cluster taus demoted to
  // exact rows): consumed, never pivots, never re-routed -- preserved as-is.
  std::set<FrozenPhysicalResidualRowKey> nongeom_exact;
  // rev154a: under an APPLIED GLOBAL plan every exact replacement is a decision of the
  // builder (design §3.2), not a candidate for re-routing: preserve all of them and let
  // the matcher see no geometric candidates.
  bool const kGlobalPlan = _planIsGlobal;
  for( auto const& rk : draft.exact_replacement_rows )
    if( kGlobalPlan || geom_cand_rows.find( rk ) == geom_cand_rows.end() )
      nongeom_exact.insert( rk );
  if( kGlobalPlan ){ geom_cand.clear(); geom_cand_rows.clear(); }

  // Index claims and candidate rows for the augmenting-path matcher.
  std::vector<std::vector<int>> claim_cands( geom_cand.size() );
  std::map<FrozenPhysicalResidualRowKey,int> rid_of;
  std::vector<FrozenPhysicalResidualRowKey> row_of;
  auto intern_row = [&]( FrozenPhysicalResidualRowKey const& rk ) -> int {
    auto it = rid_of.find( rk );
    if( it != rid_of.end() ) return it->second;
    int const id = (int)row_of.size();
    rid_of.emplace( rk, id ); row_of.push_back( rk );
    return id;
  };
  {
    size_t c = 0;
    for( auto const& [ck, cands] : geom_cand ){
      (void)ck;
      for( auto const& rk : cands ) claim_cands[c].push_back( intern_row( rk ) );
      ++c;
    }
  }
  size_t const nrid = row_of.size();
  std::vector<int> match_row( nrid, -1 );    // candidate row -> owning claim (or -1)

  // Augmenting-path search over claims.  A row is usable iff it is not a pivot
  // (pivot_set) and not the row currently being vacated (forbid_rid).
  std::function<bool(int,int,std::vector<char>&)> augment =
    [&]( int c, int forbid_rid, std::vector<char>& seen ) -> bool {
      for( int rid : claim_cands[c] ){
        if( rid == forbid_rid ) continue;
        if( seen[rid] ) continue;
        if( pivot_set.find( row_of[rid] ) != pivot_set.end() ) continue;
        seen[rid] = 1;
        if( match_row[rid] < 0 || augment( match_row[rid], forbid_rid, seen ) ){
          match_row[rid] = c; return true;
        }
      }
      return false;
    };

  // Initial system of distinct representatives (no pivots reserved yet).
  size_t sdr = 0;
  for( size_t c = 0; c < claim_cands.size(); ++c ){
    std::vector<char> seen( nrid, 0 );
    if( augment( (int)c, -1, seen ) ) ++sdr;
  }
  if( sdr != claim_cands.size() ){
    _disp( 2 ) << "OCFESLV::setup ** tau elimination: IC_STRONG "
                 "geometric-C0 claims have no system of distinct representatives: "
              << sdr << " / " << claim_cands.size() << std::endl;
    return false;
  }
  auto current_drops = [&]() -> std::set<FrozenPhysicalResidualRowKey> {
    std::set<FrozenPhysicalResidualRowKey> d;
    for( size_t rid = 0; rid < nrid; ++rid )
      if( match_row[rid] >= 0 ) d.insert( row_of[rid] );
    return d;
  };

  // Non-mutating Gram-Schmidt rank-extension test against the current basis.
  auto gs_extends = [&]( std::vector<double> const& brow ) -> bool {
    SparseVec r = to_sparse( brow );
    double const n0 = sp_norm( r );
    if( n0 <= 64. * DBL_EPSILON ) return false;
    for( auto const& q : qbasis ){
      double const dot = sp_dot( r, q );
      if( dot == 0.0 ) continue;
      sp_axpy( r, dot, q );
    }
    return sp_norm( r ) > 1e-10 * std::max( 1.0, n0 );
  };

  // Pivots over the NON-geometric, non-exact rows first.  Deferring ALL
  // geometric candidate rows keeps every redundant sibling copy free, so the
  // recovery below can take each rank-extending corner copy as a pivot while
  // its claim still drops a distinct sibling.  Taking siblings greedily here
  // instead would lock them as pivots and dead-lock the re-route (the V5 stall:
  // every conflicted claim's sibling came back PIV with no FREE alternative).
  auto _bs_now = []{ return std::chrono::steady_clock::now(); };
  auto _bs_ms  = []( std::chrono::steady_clock::time_point a,
                     std::chrono::steady_clock::time_point b )
                 { return std::chrono::duration<double,std::milli>( b - a ).count(); };
  auto _t_bs = _bs_now();
  for( auto const& kv : Brow ){
    if( Bsel.size() == ntrace ) break;
    if( geom_cand_rows.find( kv.first ) != geom_cand_rows.end() ) continue;
    if( nongeom_exact.find( kv.first ) != nongeom_exact.end() ) continue;
    try_accept_pivot( kv.first, kv.second );
  }

  // Augment toward full rank with the rank-extending geometric rows.  A FREE
  // such row is taken as a pivot directly; a DROPPED one is freed first by
  // re-routing its claim onto a sibling (avoiding the pivots), then taken.
  // Re-routing only moves claims among non-pivot copies, so the basis is never
  // disturbed and recovered pivots are never re-dropped; rank strictly
  // increases, terminating at ntrace.
  size_t recovered = 0;
  while( Bsel.size() != ntrace ){
    bool progress = false;
    for( auto const& kv : Brow ){
      auto it_r = rid_of.find( kv.first );
      if( it_r == rid_of.end() ) continue;            // not a geometric candidate
      int const rid = it_r->second;
      if( !gs_extends( kv.second ) ) continue;        // must extend the pivot rank
      if( match_row[rid] < 0 ){                        // FREE: take it directly
        try_accept_pivot( kv.first, kv.second );
        ++recovered; progress = true; break;
      }
      int const owner = match_row[rid];                // DROP: free it, then take
      match_row[rid] = -1;                             // vacate the row
      std::vector<char> seen( nrid, 0 );
      if( augment( owner, rid, seen ) ){              // re-route claim off it
        try_accept_pivot( kv.first, kv.second );      // take freed row as a pivot
        ++recovered; progress = true; break;
      }
      match_row[rid] = owner;                         // restore; try another row
    }
    if( !progress ) break;
  }

  if( Bsel.size() == ntrace ){
    std::set<FrozenPhysicalResidualRowKey> new_exact = current_drops();
    for( auto const& rk : nongeom_exact ) new_exact.insert( rk );
    draft.exact_replacement_rows = std::move( new_exact );
    if( recovered && _solve_display_level() >= 1 )
      _disp( 2 ) << "OCFESLV::setup ** tau elimination: recovered "
                << recovered << " pivot row(s) by joint tau/geometric-C0 "
                   "assignment (IC_STRONG)" << std::endl;
  }

  draft.strong_tau_elim.rank = Bsel.size();
  if( Bsel.size() != ntrace ){
    {
      // Decisive stall analysis: for each rank-extending DROP row, show its
      // claim's candidate copies and their roles, so we can tell a code bug
      // (a FREE sibling exists) from genuine contention (all siblings are
      // pivots), and whether a pivot-relocation could break the tie.
      _disp( 3 ) << "\n#### IC_STRONG joint-assignment stall analysis ####\n";
      size_t n_ext = 0, n_free_cand = 0;
      for( auto const& kv : Brow ){
        auto it_r = rid_of.find( kv.first );
        if( it_r == rid_of.end() ) continue;
        int const rid = it_r->second;
        if( match_row[rid] < 0 ) continue;             // only current drops
        if( !gs_extends( kv.second ) ) continue;       // only rank-extending
        ++n_ext;
        int const owner = match_row[rid];
        _disp( 3 ) << "  rank-ext DROP row_id=" << kv.first.row_id
                  << " flat=" << kv.first.flat << "  claim candidates:";
        for( int crid : claim_cands[owner] ){
          FrozenPhysicalResidualRowKey const& cr = row_of[crid];
          char const* role;
          if( pivot_set.find( cr ) != pivot_set.end() ) role = "PIV";
          else if( match_row[crid] == owner )           role = "DROP(self)";
          else if( match_row[crid] >= 0 )               role = "DROP(other)";
          else { role = "FREE"; ++n_free_cand; }
          _disp( 3 ) << " [row_id=" << cr.row_id << " flat=" << cr.flat
                    << " " << role << "]";
        }
        _disp( 3 ) << "\n";
      }
      // Tau rank reachable from rows that are NOT geometric candidates (i.e.
      // without consuming any corner copy a claim might need).
      std::vector<std::vector<double>> qb;
      size_t r_nc = 0;
      for( auto const& kv : Brow ){
        if( geom_cand_rows.find( kv.first ) != geom_cand_rows.end() ) continue;
        if( nongeom_exact.find( kv.first ) != nongeom_exact.end() ) continue;
        std::vector<double> r = kv.second;
        double const n0 = row_norm( r );
        if( n0 <= 64. * DBL_EPSILON ) continue;
        for( auto const& q : qb ){
          double dot = 0.0;
          for( size_t j = 0; j < ntrace; ++j ) dot += r[j] * q[j];
          for( size_t j = 0; j < ntrace; ++j ) r[j] -= dot * q[j];
        }
        double const nr = row_norm( r );
        if( nr <= 1e-10 * std::max( 1.0, n0 ) ) continue;
        for( double& v : r ) v /= nr;
        qb.push_back( std::move( r ) ); ++r_nc;
      }
      _disp( 3 ) << "  rank-extending drops still needed = " << n_ext
                << "   FREE candidate copies found = " << n_free_cand << "\n"
                << "  tau-rank from non-geometric rows alone = " << r_nc
                << " / " << ntrace << "  (taus needing a corner copy = "
                << ( ntrace - r_nc ) << ")\n"
                << "  geometric candidate rows = " << geom_cand_rows.size()
                << "   claims = " << claim_cands.size() << "\n";
      _disp( 3 ) << "  => FREE>0 means a re-route bug; FREE=0 means corner copies"
                   " are contended between pivots and drops (pivot relocation or"
                   " keep-explicit hybrid needed).\n";
      _disp( 3 ) << "##################################################\n";
    }
    _disp( 2 ) << "OCFESLV::setup ** tau elimination: IC_STRONG tau "
                 "receiver block rank-deficient after joint tau/geometric-C0 "
                 "assignment: n_trace=" << ntrace
              << " rank=" << Bsel.size()
              << " receiver_rows=" << Brow.size()
              << ".  " << ( ntrace - Bsel.size() ) << " tau(s) uncovered.";
    if( !geom_cand_rows.empty() )
      _disp( 2 ) << "  Corner copies are contended by " << claim_cands.size()
                << " geometric-C0 claim(s) (pivot relocation or keep-explicit "
                   "dummy-derivative hybrid).";
    else
      _disp( 2 ) << "  There are NO geometric-C0 claims here, so this is not corner "
                   "contention: the multiplier column block is intrinsically column "
                   "rank-deficient.  No assignment, relocation or coupling wiring can "
                   "close it -- see the rank-symmetry report (rank(B) vs rank(C)).";
    _disp( 2 ) << std::endl;

    // ---- rev40: STAGE 2 PRE-FLIGHT PROBE (report only; no behaviour change) ------------
    //
    // Stage 2 replaces the termination target of the pivot search with rank(B), on the
    // argument that B = B W W^T (since B Z = 0), so a set of rows of B is independent iff
    // the corresponding rows of B W are -- i.e. the existing sparse Gram-Schmidt already
    // selects exactly the r rows Stage 2 needs and merely has to stop at r instead of
    // failing at m.
    //
    // That argument predicts something that can fail: on EVERY rank-deficient plan the
    // search should stall at EXACTLY rank(B), never below it.  If it stalls below, the
    // pivot search has a second, independent defect and Stage 2's premise is wrong -- so
    // this is measured on the corpus BEFORE the termination target is touched.
    //
    // ntrace is the ELIMINATED tau count (n_trace_var - n_trace_var_explicit) while
    // _traceProjNTau is m over the FULL tau set, so the prediction is exact only when the
    // keep-explicit hybrid is empty.  All four counts are printed rather than differenced,
    // and the hybrid count is reported alongside so a mismatch can be attributed rather
    // than guessed at.
    // ---- rev44: HYBRID-ADEQUACY ATTRIBUTION (Option C) --------------------------------
    //
    // A stall here has THREE possible causes and the message above names only one of them.
    // Getting the attribution right is the whole point of Option C: the corpus currently has
    // no stalling plan, so the first stall that ever occurs will be the ONLY evidence
    // available about which mechanism failed, and a message that blames the model when the
    // SELECTION is at fault would send the investigation the wrong way.
    //
    // Background, so this is readable cold.  Under IC_STRONG the plan may keep n_explicit
    // taus as ordinary var[] columns (the "keep-explicit hybrid") and Schur-eliminate the
    // rest.  The pivot search must find ntrace = n_trace_var - n_explicit independent
    // receiver rows.  It can only do so if the RETAINED columns are themselves independent,
    // i.e. if the kept-explicit set spans a complement of null(B).  The selection rule is
    // LINK-anchored -- a heuristic.  Across 134 drivers it has always sufficed
    // (n_explicit == k and rank(B) + n_explicit == n_trace_var held on 78 of 78 deficient
    // rows), but nothing PROVES it must.  Option D was chosen for Stage 1 precisely to avoid
    // depending on a selection; the hybrid reintroduces one here.
    //
    // If this fires with SELECTION INADEQUATE, that is the trigger condition for Stage 2 /
    // Option A (beta as explicit unknowns spanning null(B), basis-independent).  See the
    // Stage 2 spec.  Until then Stage 2 has no failing case to be validated against, which
    // is why it is deferred rather than abandoned.
    {
      size_t const r_B   = _traceProjRank;      // rank(B) over the FULL tau set
      size_t const m_all = _traceProjNTau;      // full tau count
      size_t const n_exp = draft.n_trace_var_explicit;

      // rank(B) restricted to the ELIMINATED columns is the quantity the pivot search can
      // actually reach.  Recomputed here rather than inferred from r_B - n_exp: that
      // subtraction is only valid if the kept-explicit columns are independent of each other
      // AND of the retained ones, which is exactly the assumption under test.
      // _sparse_column_rank takes COLUMNS: V[j] is the sparse column for tau j, as
      // (receiver-row index, coefficient) pairs.  Brow is the transpose -- one DENSE row of
      // length ntrace per receiver row, keyed by FrozenPhysicalResidualRowKey -- so build the
      // columns explicitly rather than reusing the projection phase's Bv, which is assembled
      // over the FULL tau set and would answer a different question.
      std::vector<std::vector<std::pair<size_t,double>>> Bcol( ntrace );
      {
        size_t ir = 0;
        for( auto const& kv : Brow ){
          std::vector<double> const& row = kv.second;
          for( size_t j = 0; j < ntrace && j < row.size(); ++j )
            if( row[j] != 0.0 ) Bcol[j].emplace_back( ir, row[j] );
          ++ir;
        }
      }
      size_t const r_elim = _sparse_column_rank( Bcol );

      char const* verdict;
      if( Bsel.size() < r_elim )
        verdict = "** PIVOT SEARCH DEFECT: stalled BELOW rank(B_eliminated) **";
      else if( r_elim < ntrace )
        verdict = "** SELECTION INADEQUATE: the keep-explicit set does not span a complement"
                  " of null(B).  This is the Option A / beta trigger -- see the Stage 2 spec."
                  "  NOT a structural defect of the model. **";
      else
        verdict = "** STRUCTURAL: B_eliminated has full rank yet the search could not realise"
                  " it -- neither selection nor model rank explains this. **";

      _disp( 2 ) << "  [hybrid-adequacy] pivot stall: stalled_at=" << Bsel.size()
                << " rank(B_eliminated)=" << r_elim
                << " ntrace_elim=" << ntrace
                << " rank(B_all)=" << r_B
                << " m_all_tau=" << m_all
                << " k=" << _traceProjK
                << " n_explicit=" << n_exp << std::endl;
      _disp( 2 ) << "  [hybrid-adequacy] " << verdict << std::endl;
    }
    return false;
  }

  double const _t_select_ms = _bs_ms( _t_bs, _bs_now() );
  _t_bs = _bs_now();

#ifdef MC__OCFESLV_BS_STRUCTURE_PROBE
  // -------------------------------------------------------------------------
  // B_S structure probe (default-off, measure-before-fixing).  The dense
  // inversion below is O(ntrace^3); whether that is necessary depends on the
  // sparsity pattern of B_S = [Bsel].  Characterize it cheaply so the fix can
  // be chosen with data instead of a guess:
  //   - connected components of the symmetrized coupling graph (i~j if
  //     B_S(i,j) or B_S(j,i) != 0) == the irreducible diagonal blocks.  If the
  //     blocks are small (taus cluster per interface face), block-wise dense
  //     inversion is sum(b^3) << ntrace^3 and the inverse stays block-local
  //     (consistent with the projection surviving the 128*eps threshold).
  //   - if there is one giant component, report its half-bandwidth in the
  //     current ordering: small band -> band solve; otherwise sparse LU.
  // Deliberately does NOT trigger the O(ntrace^3) SVD (separate gate), so the
  // structure run is cheap.  Cost here is one O(ntrace^2) scan of the dense Bsel.
  {
    std::vector<size_t> uf( ntrace );
    for( size_t i = 0; i < ntrace; ++i ) uf[i] = i;
    auto find = [&]( size_t x ) -> size_t {
      while( uf[x] != x ){ uf[x] = uf[uf[x]]; x = uf[x]; } return x; };
    auto unite = [&]( size_t a, size_t b ){ a=find(a); b=find(b); if(a!=b) uf[a]=b; };
    size_t nnz = 0, band = 0;
    std::vector<std::vector<size_t>> adj( ntrace );
    for( size_t i = 0; i < ntrace; ++i )
      for( size_t j = 0; j < ntrace; ++j )
        if( Bsel[i][j] != 0.0 ){
          ++nnz; unite( i, j );
          size_t const d = ( i > j ) ? ( i - j ) : ( j - i );
          if( d > band ) band = d;
          if( i != j ){ adj[i].push_back( j ); adj[j].push_back( i ); }
        }
    std::map<size_t,size_t> comp_size;
    for( size_t i = 0; i < ntrace; ++i ) ++comp_size[ find(i) ];
    std::vector<size_t> sizes;
    sizes.reserve( comp_size.size() );
    for( auto const& kv : comp_size ) sizes.push_back( kv.second );
    std::sort( sizes.begin(), sizes.end() );
    size_t const maxb = sizes.empty() ? 0 : sizes.back();
    size_t const medb = sizes.empty() ? 0 : sizes[ sizes.size()/2 ];
    long double sum_b3 = 0.0L;
    for( size_t b : sizes ) sum_b3 += (long double)b*(long double)b*(long double)b;
    long double const dense3 = (long double)ntrace*(long double)ntrace*(long double)ntrace;
    _disp( 3 ) << "    [iplan-bs-struct] ntrace=" << ntrace
              << " nnz=" << nnz << " ("
              << ( 100.0 * (double)nnz / ( (double)ntrace*(double)ntrace ) ) << "% dense)"
              << " components=" << sizes.size()
              << " max_block=" << maxb << " median_block=" << medb
              << " half_bandwidth=" << band << "\n"
              << "    [iplan-bs-struct] block-invert cost ~sum(b^3)=" << (double)sum_b3
              << " vs dense ntrace^3=" << (double)dense3
              << "  speedup~" << ( sum_b3 > 0 ? (double)(dense3/sum_b3) : 0.0 ) << "x"
              << "  (max_block^3=" << (double)( (long double)maxb*maxb*maxb ) << ")\n";

    // Cuthill-McKee reordered bandwidth: is the large natural-order bandwidth
    // intrinsic, or just scrambled tau numbering?  Small CM bandwidth => a
    // band/sparse-LU factor-solve is O(ntrace*w^2) << dense O(ntrace^3) and the
    // dense inversion is replaceable; large CM bandwidth => the coupling is
    // genuinely long-range, dense is hard to beat, and the coupling itself is
    // the thing to question.  Also report max_degree (hubs inflate fill).
    {
      std::vector<size_t> degree( ntrace );
      size_t maxdeg = 0;
      for( size_t i = 0; i < ntrace; ++i ){
        std::sort( adj[i].begin(), adj[i].end() );
        adj[i].erase( std::unique( adj[i].begin(), adj[i].end() ), adj[i].end() );
        degree[i] = adj[i].size();
        if( degree[i] > maxdeg ) maxdeg = degree[i];
      }
      std::vector<long long> pos( ntrace, -1 );
      std::vector<size_t> roots( ntrace );
      for( size_t i = 0; i < ntrace; ++i ) roots[i] = i;
      std::sort( roots.begin(), roots.end(),
                 [&]( size_t a, size_t b ){ return degree[a] < degree[b]; } );
      size_t placed = 0;
      for( size_t r : roots ){
        if( pos[r] >= 0 ) continue;
        std::vector<size_t> q; size_t qh = 0;
        pos[r] = (long long)placed++; q.push_back( r );
        while( qh < q.size() ){
          size_t const u = q[qh++];
          std::vector<size_t> nb( adj[u].begin(), adj[u].end() );
          std::sort( nb.begin(), nb.end(),
                     [&]( size_t a, size_t b ){ return degree[a] < degree[b]; } );
          for( size_t v : nb )
            if( pos[v] < 0 ){ pos[v] = (long long)placed++; q.push_back( v ); }
        }
      }
      size_t cm_band = 0;
      for( size_t i = 0; i < ntrace; ++i )
        for( size_t j : adj[i] ){
          long long d = pos[i] - pos[j]; if( d < 0 ) d = -d;
          if( (size_t)d > cm_band ) cm_band = (size_t)d;
        }
      long double const band_lu = (long double)ntrace*(long double)cm_band*(long double)cm_band;
      _disp( 3 ) << "    [iplan-bs-struct] max_degree=" << maxdeg
                << " CM_bandwidth=" << cm_band << " (natural=" << band << ")\n"
                << "    [iplan-bs-struct] band-LU cost ~ntrace*CM_band^2=" << (double)band_lu
                << " vs dense ntrace^3=" << (double)dense3
                << "  speedup~" << ( band_lu > 0 ? (double)(dense3/band_lu) : 0.0 ) << "x\n";
    }
  }
#endif

#if defined(CRONOS__WITH_UMFPACK)
  // ===== Sparse AMD-ordered LU path (UMFPACK) =====
  // B_S is ~0.05% dense (CM bandwidth ~30% of ntrace, max-degree ~9): the dense
  // O(ntrace^3) inverse is wasteful.  Factor B_S once with UMFPACK's AMD fill-
  // reducing ordering and solve B_S^T y_k = b_k per non-pivot receiver row, since
  // y_k = b_k^T B_S^{-1} = B_S^{-T} b_k.  UMFPACK_At applies the transpose solve off
  // the *same* factorization -- no explicit transpose or inverse is built -- and the
  // result (draft.strong_tau_elim.projection) is identical to the dense path.
  //
  // Pack B_S into compressed-sparse-column (CSC): B_S(i,j)=Bsel[i][j].  Iterating i
  // ascending within each column j yields ascending row indices, as UMFPACK requires.
  std::vector<std::int64_t> Ap( ntrace + 1, 0 );
  std::vector<std::int64_t> Ai;
  std::vector<double>       Ax;
  Ai.reserve( 8 * ntrace );
  Ax.reserve( 8 * ntrace );
  for( size_t j = 0; j < ntrace; ++j ){
    for( size_t i = 0; i < ntrace; ++i ){
      double const v = Bsel[i][j];
      if( v != 0.0 ){ Ai.push_back( (std::int64_t)i ); Ax.push_back( v ); }
    }
    Ap[j+1] = (std::int64_t)Ai.size();
  }

  double Control[UMFPACK_CONTROL];
  double Info[UMFPACK_INFO];
  umfpack_dl_defaults( Control );          // AMD/COLAMD fill-reducing ordering (default)
  void* Symbolic = nullptr;
  void* Numeric  = nullptr;
  std::int64_t ust = umfpack_dl_symbolic( (std::int64_t)ntrace, (std::int64_t)ntrace,
                                          Ap.data(), Ai.data(), Ax.data(),
                                          &Symbolic, Control, Info );
  if( ust != UMFPACK_OK ){
    std::cerr << "OCFESLV::setup ** tau elimination: UMFPACK symbolic "
                 "factorization failed (status=" << ust << ")" << std::endl;
    if( Symbolic ) umfpack_dl_free_symbolic( &Symbolic );
    return false;
  }
  ust = umfpack_dl_numeric( Ap.data(), Ai.data(), Ax.data(),
                            Symbolic, &Numeric, Control, Info );
  umfpack_dl_free_symbolic( &Symbolic );
  if( ust != UMFPACK_OK ){
    _disp( 2 ) << "OCFESLV::setup ** tau elimination: selected tau pivot "
                 "block is numerically singular (UMFPACK numeric status=" << ust
              << ")" << std::endl;
    if( Numeric ) umfpack_dl_free_numeric( &Numeric );
    return false;
  }
  double const _t_invert_ms = _bs_ms( _t_bs, _bs_now() );

  // Conditioning: UMFPACK returns a reciprocal condition estimate for free (no SVD).
  if( _solve_display_level() >= 2 )
    _disp( 2 ) << "OCFESLV::setup ** trace conditioning: B_S UMFPACK rcond="
              << std::scientific << std::setprecision(3) << Info[UMFPACK_RCOND]
              << " (~1/cond)  (ntrace=" << ntrace << ")" << std::endl;

  // Solve B_S^T y = b_k per non-pivot receiver row; threshold to the sparse projection
  // (128*eps, identical to the dense path).  One factorization, many cheap solves.
  _t_bs = _bs_now();
  std::vector<double> bk( ntrace ), yk( ntrace );
  for( auto const& kv : Brow ){
    if( draft.strong_tau_elim.is_pivot( kv.first ) ) continue;
    for( size_t j = 0; j < ntrace; ++j ) bk[j] = kv.second[j];
    std::int64_t const sst = umfpack_dl_solve( UMFPACK_At, Ap.data(), Ai.data(), Ax.data(),
                                               yk.data(), bk.data(), Numeric, Control, Info );
    if( sst != UMFPACK_OK ){
      std::cerr << "OCFESLV::setup ** tau elimination: UMFPACK solve "
                   "failed (status=" << sst << ")" << std::endl;
      umfpack_dl_free_numeric( &Numeric );
      return false;
    }
    std::vector<std::pair<size_t,double>> corr;
    for( size_t p = 0; p < ntrace; ++p )
      if( std::fabs( yk[p] ) > 128. * DBL_EPSILON )
        corr.emplace_back( p, yk[p] );
    if( !corr.empty() )
      draft.strong_tau_elim.projection[ kv.first ] = std::move( corr );
  }
  umfpack_dl_free_numeric( &Numeric );
  double const _t_project_ms = _bs_ms( _t_bs, _bs_now() );

#else
  // ===== Dense Armadillo inverse path (fallback; no UMFPACK) =====
  // Invert B_S using Armadillo.  The selected pivot rows are already stored
  // in Bsel; pack them into an arma::mat and compute B_S^{-1}.
  arma::mat B_S( (arma::uword)ntrace, (arma::uword)ntrace );
  for( size_t i = 0; i < ntrace; ++i )
    for( size_t j = 0; j < ntrace; ++j )
      B_S( (arma::uword)i, (arma::uword)j ) = Bsel[i][j];

  arma::mat B_S_inv;
  if( !arma::inv( B_S_inv, B_S ) ){
    std::cerr << "OCFESLV::setup ** tau elimination: selected "
                 "tau pivot block is numerically singular (Armadillo inv failed)"
              << std::endl;
    return false;
  }
  double const _t_invert_ms = _bs_ms( _t_bs, _bs_now() );

  // Record pivot condition information.  This is diagnostic only -- pivot_max/min
  // are never consumed by the elimination -- and a values-only SVD of B_S is
  // O(ntrace^3).  Measured on PDE5 IC_STRONG it is ~6 s at ntrace=3236 (~6.6x the
  // inversion, ~85% of this function's time) and dominates setup on big systems.
  // Gate it behind DISPLAY_LEVEL *and* a size cap so the conditioning stays for
  // small systems (where it is instant and useful) but never inflates large-system
  // setup.  Override the cap with -DMC__OCFESLV_BS_COND_MAX_N=<n>.
#ifndef MC__OCFESLV_BS_COND_MAX_N
#define MC__OCFESLV_BS_COND_MAX_N 1500
#endif
  if( _solve_display_level() >= 2 ){
    if( ntrace > (size_t)MC__OCFESLV_BS_COND_MAX_N ){
      _disp( 2 ) << "OCFESLV::setup ** trace conditioning: B_S"
                   " SVD skipped (ntrace=" << ntrace << " > " << MC__OCFESLV_BS_COND_MAX_N
                << "; O(ntrace^3) diagnostic).  Raise -DMC__OCFESLV_BS_COND_MAX_N to force."
                << std::endl;
    } else {
      arma::vec sv = arma::svd( B_S );
      if( sv.n_elem > 0 ){
        draft.strong_tau_elim.pivot_max = sv(0);
        draft.strong_tau_elim.pivot_min = sv( sv.n_elem - 1 );
        _disp( 2 ) << "OCFESLV::setup ** trace conditioning: B_S sigma:"
                  << " max=" << std::scientific << std::setprecision(3) << sv(0)
                  << " min=" << sv( sv.n_elem - 1 )
                  << " cond=" << ( sv( sv.n_elem-1 ) > 0
                                   ? sv(0)/sv( sv.n_elem-1 )
                                   : std::numeric_limits<double>::infinity() )
                  << "  (ntrace=" << ntrace << ")" << std::endl;
      }
    }
  }

  // For every non-pivot receiver row k, freeze y = B_k * B_S^{-1} as a
  // sparse vector of (pivot_index, coefficient) pairs.  Evaluation emits
  // R_k - sum_p y_p * R_S[p].  Rows outside B have no correction.
  _t_bs = _bs_now();
  arma::mat const B_S_invT = B_S_inv.t();   // contiguous columns for the sparse accumulation below
  for( auto const& kv : Brow ){
    if( draft.strong_tau_elim.is_pivot( kv.first ) ) continue;
    // y = b_k * B_S^{-1}.  b_k is sparse (each receiver row couples to only a few
    // taus), and y(p) = sum_j b_k(j) * B_S_inv(j,p) = sum_j b_k(j) * B_S_invT(p,j),
    // so accumulate the few nonzero-weighted (contiguous) columns of B_S_invT
    // rather than forming the dense rowvec and a full O(ntrace^2) rowvec*matrix
    // product.  Skipping a zero b_k(j) is exact -- it adds 0 to every y(p).
    arma::vec y( (arma::uword)ntrace, arma::fill::zeros );
    for( size_t j = 0; j < ntrace; ++j ){
      double const bj = kv.second[j];
      if( bj != 0.0 ) y += bj * B_S_invT.col( (arma::uword)j );
    }
    std::vector<std::pair<size_t,double>> corr;
    for( size_t p = 0; p < ntrace; ++p ){
      if( std::fabs( y( (arma::uword)p ) ) > 128. * DBL_EPSILON )
        corr.emplace_back( p, y( (arma::uword)p ) );
    }
    if( !corr.empty() )
      draft.strong_tau_elim.projection[ kv.first ] = std::move( corr );
  }
  double const _t_project_ms = _bs_ms( _t_bs, _bs_now() );
#endif // CRONOS__WITH_UMFPACK

  if( _solve_display_level() >= 2 ){
    size_t bs_nnz = 0;
    for( size_t i = 0; i < ntrace; ++i )
      for( size_t j = 0; j < ntrace; ++j )
        if( Bsel[i][j] != 0.0 ) ++bs_nnz;
    _disp( 3 ) << "    [iplan-bs] select=" << _t_select_ms << "ms invert=" << _t_invert_ms
              << "ms project=" << _t_project_ms << "ms  B_S nnz=" << bs_nnz
              << "/" << ( (double)ntrace*(double)ntrace ) << " ("
              << ( 100.0 * (double)bs_nnz / ( (double)ntrace*(double)ntrace ) )
              << "% dense)  recv_rows=" << Brow.size() << " ntrace=" << ntrace << "\n";
  }

  draft.strong_tau_elim.ready = true;
#ifdef MC__OCFESLV_INTERFACE_ASSIGN_DEBUG
  if( _solve_display_level() >= 1 )
    _disp( 2 ) << "OCFESLV::setup ** IC_STRONG tau-eliminated reduction: n_trace="
              << ntrace << " n_trace_constraint=" << draft.n_trace_constraint
              << " exact_replacements=" << draft.exact_replacement_rows.size()
              << " receiver_rows=" << Brow.size()
              << " pivot_min=" << draft.strong_tau_elim.pivot_min
              << " pivot_max=" << draft.strong_tau_elim.pivot_max
              << " projected_rows=" << draft.strong_tau_elim.projection.size()
              << std::endl;
#endif
  return true;
}

inline std::map<size_t,size_t>
OCFEBASE::_defining_rows() const
{
  std::map<size_t,size_t> out;
  // declared
  for( size_t ei = 0; ei < _mEqn.size(); ++ei )
    if( _mEqn[ei].opt && _mEqn[ei].opt->defines.dag() )
      out[ static_cast<size_t>( _mEqn[ei].opt->defines.id().second ) ] = ei;
  // inferred
  std::map<size_t, std::vector<std::pair<size_t,int>>> cand;   // state -> [(row, letter)]
  for( size_t ei = 0; ei < _mEqn.size(); ++ei ){
    if( !_mEqn[ei].opt || !_mEqn[ei].opt->participate_in_classification ) continue;
    if( _evoDerivRows.count( ei ) ) continue;                     // a definition is not an evolution row
    try {
      std::vector<FFInv> invVar( _vVar.size() );
      for( size_t k = 0; k < _vVar.size(); ++k ) invVar[k].indep( (int)k );
      std::vector<FFInv> invInp( _vInp.size() ), invCst( _vCst.size() ), invDom( _vDom.size() );
      FFInv invEqn; FFVar ev = _mEqn[ei].var;
      auto sg = _dag->subgraph( 1, &ev );
      _dag->eval( sg, 1, &ev, &invEqn, _vVar.size(), _vVar.data(), invVar.data(),
                  _vInp.size(), _vInp.data(), invInp.data(), _vCst.size(), _vCst.data(), invCst.data(),
                  _vDom.size(), _vDom.data(), invDom.data() );
      for( auto const& [ind, ty] : invEqn.inv() ){
        if( ind < 0 || (size_t)ind >= _vVar.size() ) continue;
        if( ty != FFInv::TYPE::L && ty != FFInv::TYPE::S ) continue;
        cand[ static_cast<size_t>( _vVar[ind].id().second ) ].push_back( { ei, (int)ty } );
      }
    } catch( ... ){ continue; }
  }
  std::map<size_t,size_t> inferred;
  for( auto const& [sid, rows] : cand ){
    if( out.count( sid ) ) continue;                               // declared wins
    size_t nL = 0, nS = 0, rL = 0, rS = 0;
    for( auto const& [ei, ty] : rows ){
      if( ty == (int)FFInv::TYPE::L ){ ++nL; rL = ei; } else { ++nS; rS = ei; } }
    if( nL == 1 ) inferred[sid] = rL;                              // a unique L beats any S
    else if( nL == 0 && nS == 1 ) inferred[sid] = rS;              // a unique S is taken
    // otherwise ambiguous: no inference
  }
  // A ROW DEFINES AT MOST ONE STATE.  `p - g = 0` is L in both p and g; the chain means p := g, but nothing in the
  // row says which, and inferring it for BOTH admits g through a row that defines p -- OCFE_zchain went from 9 PASS
  // to 3 that way.  Where two states would take the same row, neither gets it: declare the one that is defined.
  std::map<size_t,size_t> row_users;
  for( auto const& [sid, ei] : inferred ) ++row_users[ei];
  for( auto const& [sid, ei] : inferred ) if( row_users[ei] == 1 ) out[sid] = ei;
  return out;
}

inline void
OCFEBASE::_prepare_interface_plan_tables
( std::map<t_InterfaceClaimKey,size_t>& claims_out,
  InterfacePlanDraft& draft,
  bool store_plan ) const
{
  _cornerSuppressedClaims.clear();   // rev153b
  _pinSuppressedClaims.clear();      // rev153c
  // The defining row of each state (declared, else inferred by FFInv L/S) is a MODEL-level fact,
  // invariant for the whole plan build.  It was computed per NODE inside the per-dim loop below
  // (10,338 calls on OCFE_PDE3), and each call evaluates every equation's DAG through FFInv:
  // setup went 1.2 s -> 22 s between v5.0.0 and v5.1.0 on that driver, every solve-side number
  // unchanged.  Computed once here; the loop reads it.
  std::map<size_t,size_t> const defrow_once = _defining_rows();
  _rankImpliedClaims.clear();        // rev160
  size_t const invalid_claim_id = std::numeric_limits<size_t>::max();

  std::map<t_InterfaceClaimKey,size_t> interface_claims_master;

  // Per-equation subgraph cache.  The receiver-predicate lambdas below test
  // state/direction incidence by building _dag->subgraph(1,&eqnvar), and they are
  // invoked once per interface node in the main enumeration loop.  That subgraph
  // depends only on the equation, never on the node, so rebuilding it per node is
  // an O(interface_nodes) redundancy that dominates the interface-plan build on
  // large high-order problems (the 18 s seen at n=12288).  Cache by equation id
  // and hand out stable references (std::map nodes are not invalidated by later
  // insertions).
  std::map<long, FFSubgraph> _eqn_subgraph_cache;
  auto eqn_subgraph = [&]( FFVar const& ev ) -> FFSubgraph const& {
    long const key = ev.id().second;
    auto it = _eqn_subgraph_cache.find( key );
    if( it == _eqn_subgraph_cache.end() )
      it = _eqn_subgraph_cache.emplace( key, _dag->subgraph( 1, &ev ) ).first;
    return it->second;
  };

  // --- TEMPORARY per-predicate timing accumulators (DISPLAY_LEVEL>=2) ---------
  // Localizes which per-node receiver predicate dominates _prepare's cost.
  double _acc_geom=0.0, _acc_aux=0.0, _acc_satc=0.0, _acc_weak=0.0, _acc_perdim=0.0, _acc_emit=0.0;
  size_t _cnt_geom=0,   _cnt_aux=0,   _cnt_satc=0,   _cnt_weak=0,   _cnt_perdim=0,   _cnt_emit=0;
  double _acc_keep=0.0, _acc_rtype=0.0;
  size_t _cnt_keep=0,   _cnt_rtype=0;
  struct _AccT {
    double& a; size_t& c; std::chrono::steady_clock::time_point t0;
    _AccT( double& aa, size_t& cc ) : a(aa), c(cc), t0(std::chrono::steady_clock::now()) {}
    ~_AccT(){ a += std::chrono::duration<double,std::milli>(
                     std::chrono::steady_clock::now() - t0 ).count(); ++c; }
  };
  struct _AccReport {
    bool v; double&g; size_t&gc; double&a; size_t&ac; double&s; size_t&sc; double&w; size_t&wc;
    double&pd; size_t&pdc; double&em; size_t&emc; double&kp; size_t&kpc; double&rt; size_t&rtc;
    ~_AccReport(){
      if( v )
        std::cout << "    [iplan-acc] per_dim=" << pd << "ms/" << pdc
                  << "  emit_states=" << em << "ms/" << emc
                  << "  keep_eqn=" << kp << "ms/" << kpc
                  << "  resolve_iface=" << rt << "ms/" << rtc
                  << "\n                geom_c0=" << g << "ms/" << gc
                  << "  aux_alg=" << a << "ms/" << ac
                  << "  sat_coup=" << s << "ms/" << sc
                  << "  weak_ids=" << w << "ms/" << wc << "\n";
    }
  } _accrep{ _solve_display_level()>=2, _acc_geom,_cnt_geom, _acc_aux,_cnt_aux,
             _acc_satc,_cnt_satc, _acc_weak,_cnt_weak,
             _acc_perdim,_cnt_perdim, _acc_emit,_cnt_emit,
             _acc_keep,_cnt_keep, _acc_rtype,_cnt_rtype };

  // The exact IC_STRONG/IC_TRACE plan is intentionally derived from the same
  // row/state/direction incidence as IC_WEAK SAT.  The coupling cache is a
  // setup-time prerequisite: setup() builds it immediately after
  // classification and before any interface receiver/trace plan is frozen.
  // Reaching this point with a stale cache means some caller tried to build
  // or use an interface plan outside the setup-time sequence.
  if( !_satCouplingBuilt
   || _satCouplingClassified != _classified
   || _satCouplingSerial != _classificationSerial ){
    _disp( 2 ) << "OCFESLV::_prepare_interface_plan_tables ** SAT coupling "
                 "cache is not current; setup() must call "
                 "_ensure_sat_coupling() before building interface plans."
              << std::endl;
    throw Exceptions( Exceptions::INTERNAL );
  }

  auto encode_tensor_face_from_eqn_dims = []
  ( size_t idim, auto const& dims, std::vector<size_t> const& midx ) -> size_t
  {
    // Globally identify a face degree of freedom by the complete
    // tensor-product coordinate of every tangential dimension.  The previous
    // key used only the tangential local node number; for multiple elements
    // this aliases, for example, local z-node 1 in all axial elements.  Here
    // each tangential coordinate is encoded as (element,node) with radix
    // n_elem*n_node, while collapsed LB/UB dimensions are pinned to their
    // physical endpoint element and endpoint node.
    size_t face = 0;
    size_t mult = 1;
    for( size_t d = 0; d < dims.size(); ++d ){
      if( d == idim ) continue;
      if( !dims[d].pdom ) continue;

      size_t ed = dims[d].iel;
      size_t kd = 0;
      switch( dims[d].lim ){
        case FFDom::LB:
          ed = 0;
          kd = 0;
          break;
        case FFDom::UB:
          ed = dims[d].pdom->n_elem - 1;
          kd = dims[d].pdom->n_node - 1;
          break;
        default:
          kd = midx[d];
          break;
      }

      size_t const radix = dims[d].pdom->n_elem * dims[d].pdom->n_node;
      face += ( ed * dims[d].pdom->n_node + kd ) * mult;
      mult *= radix;
    }
    return face;
  };

  auto make_claim_key = [&]
  ( int block_id, size_t state_id, size_t eqn_id, FFVar const& dom,
    size_t iel, size_t inode, size_t n_node, size_t n_elem,
    size_t idim, auto const& dims, std::vector<size_t> const& midx,
    int claim_kind, t_InterfaceClaimKey& key ) -> bool
  {
    if( n_node <= 1 ) return false;

    size_t iel_lo = 0;
    if( iel > 0 && inode == 0 )
      iel_lo = iel - 1;
    else if( iel + 1 < n_elem && inode + 1 == n_node )
      iel_lo = iel;
    else
      return false;

    size_t const face = encode_tensor_face_from_eqn_dims( idim, dims, midx );

    key = { block_id, state_id, eqn_id, static_cast<size_t>( dom.id().second ),
            iel_lo, face, claim_kind };
    return true;
  };


  auto push_unique_state_id = []
  ( std::vector<size_t>& ids, size_t id )
  {
    if( std::find( ids.begin(), ids.end(), id ) == ids.end() )
      ids.push_back( id );
  };


  // IC_STRONG row ownership and IC_TRACE exact-claim discovery are derived
  // from the IC_WEAK SAT receiver graph.  The structural donor taxonomy below
  // is deliberately secondary: it ranks admissible weak-SAT receiver rows, but
  // it is not allowed to create strong donors or trace claims that the weak
  // formulation would not have penalised in the same (row,state,direction) slot.


  auto state_by_id = [&]( size_t sid ) -> FFVar const*
  {
    for( auto const& [st, dep] : _mVar ){
      (void)dep;
      if( static_cast<size_t>( st.id().second ) == sid ) return &st;
    }
    return nullptr;
  };

  auto state_depends_on_dom = [&]( size_t sid, FFVar const& dom ) -> bool
  {
    FFVar const* pst = state_by_id( sid );
    if( !pst ) return false;
    auto it = _mVar.find( *pst );
    if( it == _mVar.end() ) return false;
    for( auto const& d : it->second )
      if( d.id() == dom.id() ) return true;
    return false;
  };

  // FIX 3.  Algebraic-alias resolver.  A hand-written reduction chain ends in a
  // state that is structurally an auxiliary but was not minted by reduce_order(),
  // so _is_auxiliary_state() (provenance) misses it and its interior-face copies
  // are never tied -- the chain-3 defect (OCFE_bisect2 --chain 3).  A NON-aux
  // collocated state s is admitted as an algebraic alias of a registered aux iff
  // some participating equation is AFFINE in s (FFDep type L), its remaining
  // collocated-state content is EXACTLY ONE state, that state is a reduce_order()
  // auxiliary that differentiates in face_dom, AND that aux is NOT the derivative
  // of s itself (else the equation is s's own LINK, and s is a primitive whose
  // continuity the principal-symbol path already carries -- the u LINK false
  // positive measured on PDE20f).  Chain 3 survives: RDEF is L*r - Dxi_q, r linear,
  // remainder {Dxi_q}, root(Dxi_q)=q != r.  Returns the aliased aux or nullptr.
  // Is state sid differentiated in direction dom (an FFPartial with sid as operand and dom among its directions) in
  // any participating equation?
  auto state_differentiated_in = [&]( std::pair<FFVar::TYPE,long> const& sid, FFVar const& dom ) -> bool
  {
    for( size_t ei = 0; ei < _mEqn.size(); ++ei ){
      if( !_mEqn[ei].opt->participate_in_classification ) continue;
      FFSubgraph const& sgq = eqn_subgraph( _mEqn[ei].var );
      for( auto const& op : sgq.l_op ){
        if( !op || !op->sameid( typeid( FFPartial ) ) || op->varin.empty() || !op->varin[0] ) continue;
        if( op->varin[0]->id() != sid ) continue;
        auto const* pop = mc::type_cast<FFPartial const>( op );
        if( pop ) for( auto const& [dv, ordr] : pop->Indep().expr ){ (void)ordr; if( dv.id() == dom.id() ) return true; }
      }
    }
    return false;
  };

  auto alias_of_auxiliary = [&]
  ( size_t sid, FFVar const& face_dom ) -> FFVar const*
  {
    FFVar const* ps = state_by_id( sid );
    if( !ps || _is_auxiliary_state( *ps ) ) return nullptr;
    for( size_t ei = 0; ei < _mEqn.size() && ei < _depEqn.size(); ++ei ){
      if( !_mEqn[ei].opt->participate_in_classification ) continue;
      bool s_present = false, s_linear = false;
      size_t n_other_states = 0;
      FFVar const* only_other = nullptr;
      for( auto const& kv : _depEqn[ei].dep() ){
        long const fdid = kv.first;
        mc::FFDep::TYPE const tp = kv.second;
        auto it = _mDep.find( fdid );
        if( it == _mDep.end() || it->second.first != STATE ) continue;
        FFVar const& sv = _vVar[ it->second.second ];
        if( static_cast<size_t>( sv.id().second ) == sid ){
          s_present = true;
          s_linear  = ( tp == mc::FFDep::L );
        }
        else {
          ++n_other_states;
          only_other = &_vVar[ it->second.second ];
        }
      }
      // NOTE (fix5): the strict FFDep-linearity gate (s_linear) is NO LONGER applied at
      // entry.  On a graded / ALE mesh the multiply-form defining row is  s*C(xi) - P_xi
      // with C a state-INDEPENDENT domain coefficient (e.g. x_xi = aa*dm*cosh(...)); FFDep
      // tags s as nonlinear (up to type N) purely because of that coefficient, even though
      // the row is affine in s.  The AUX case keeps requiring s_linear (unchanged from
      // fix4); only the PRIMITIVE case is made coefficient-agnostic, gated on the
      // structural certificate below instead.  Detect-only probe3 confirmed this admits
      // exactly the graded chain terminals (r in chain 3, q in chain 2) and NOTHING new on
      // the uniform corpus (PDE20f/PDE12m/PDE2/3/4/PSA1/MBC5 all clean).
      if( !( s_present && n_other_states == 1 && only_other ) ) continue;

      // Structural roles of s and the sole other state P, from the equation subgraph:
      //   s_under_deriv : s is the differentiated operand (varin[0]) of an FFPartial
      //   P_under_face  : P is differentiated in the FACE direction
      //   P_bare        : P appears bare (in a non-derivative op, or not as the
      //                   differentiated operand) -- an s*P bilinear coupling, which the
      //                   old s_linear==L gate used to exclude implicitly and which the
      //                   coefficient-agnostic primitive case must now exclude explicitly.
      bool s_under_deriv = false, P_under_face = false, P_bare = false;
      {
        FFSubgraph const& sgq = eqn_subgraph( _mEqn[ei].var );
        for( auto const& op : sgq.l_op ){
          if( !op ) continue;
          bool const isPart = op->sameid( typeid( FFPartial ) );
          for( size_t k = 0; k < op->varin.size(); ++k ){
            FFVar const* vin = op->varin[k];
            if( !vin ) continue;
            if( vin->id() == ps->id() && isPart && k == 0 ) s_under_deriv = true;
            if( vin->id() == only_other->id() ){
              if( isPart && k == 0 ){
                auto const* pop = mc::type_cast<FFPartial const>( op );
                if( pop ) for( auto const& [dv, ordr] : pop->Indep().expr ){
                  (void)ordr; if( dv.id() == face_dom.id() ) P_under_face = true; }
              }
              else P_bare = true;
            }
          }
        }
      }

      if( _is_auxiliary_state( *only_other ) ){
        // AUX CASE (reduce_order chain) -- unchanged from fix4, INCLUDING the strict
        // FFDep-linearity requirement (a reduce_order aux always enters its defining row
        // linearly; the graded-coefficient relaxation is specific to the primitive case).
        if( !s_linear ) continue;
        {
          FFVar const* root = _auxiliary_primitive_root( *only_other );
          if( root && root->id() == ps->id() ) continue;   // s's own LINK
        }
        if( !_auxiliary_state_in_direction( *only_other, face_dom ) ) continue;
        return only_other;
      }
      else {
        // PRIMITIVE CASE (hand-written / varying-coefficient first-order chain, e.g.
        // RDEF = r*x_xi - OpP(q,xi)): coefficient-agnostic structural certificate.  Admit
        // s iff it appears ONLY bare (never differentiated), its sole state-neighbour P is
        // differentiated in the FACE direction, and P never appears bare (no s*P bilinear).
        // s-bare discriminates against flux equations (candidate under a derivative);
        // P-bare exclusion replaces the linearity guard the entry s_linear used to provide.
        if( s_under_deriv ) continue;
        if( !P_under_face ) continue;
        if( P_bare ) continue;
        // 2026-10-02: in the EVOLUTION direction, s must ALSO not carry an evolution derivative in ANY equation.  A
        // DIFFERENTIAL state can appear bare in another state's equation (dCC/dt = k CA^2 in OCFE_ODE2; y' = x; x' = u
        // in the pendulum DAEs) and was taken for an alias of that state's time derivative -- under the causal
        // evolution interface its continuity in t was then dropped as a "flux" (aux_flux), decoupling its time
        // elements.  Spatial directions keep the previous rule: a hand-declared flux q = u_z IS differentiated in x
        // in u's balance and must stay an alias there.
        if( _evolution_dom_set && _evolution_dom_var.dag() && face_dom.id() == _evolution_dom_var.id()
         && state_differentiated_in( ps->id(), face_dom ) ) continue;
        return only_other;
      }
    }
    return nullptr;
  };

  // Does a distributed input JUMP across interior interfaces in the EVOLUTION
  // direction?  InputContinuity ( DISCONTINUOUS default / CONTINUOUS_C0 / SMOOTH_C1 )
  // is the user's assertion of solution smoothness at the time interface: it exists
  // precisely to enforce higher-order solution continuity when a time-varying input
  // is present but known smooth there.  A jumping (DISCONTINUOUS) input kinks the
  // solution -> the evolution interface must be causal (primitives) / LINK-pinned
  // (flux).  When every evolution-direction input is declared continuous/smooth, OR
  // there is no distributed input at all (a pure-forcing PDE has _mInp empty here),
  // the solution is smooth and the full SYMMETRIC higher-order continuity is the
  // correct, more accurate treatment -- so the causal/flux-drop guard below is made
  // inert.  Mirrors discontinuous_input_in() in _resolve_interface_type (same
  // _mInp / _input_dom / n_elem<2 / _mInpCont test), evaluated once for the
  // evolution domain.
  // rev72: the gate that admits the causal evolution-direction interface.  Its name and
  // body concern INPUTS distributed in t with an interior t-interface -- a condition about
  // _mInp, not about states -- so on a model with no such input the causal block is skipped
  // and t-faces receive ordinary two-way claims REGARDLESS of whether the evolution domain
  // is known.  Whether that is what separates MMPDE27's two configurations is exactly what
  // the [evoface] line below reports.
  bool const evolution_input_may_jump = _evolution_input_may_jump();
  {
    static bool const kEvoFace = Options::_env_flag( "CRONOS_AUDIT_EVOFACE", false );
    if( kEvoFace )
      _disp( 3 ) << "  [evoface] evolution_input_may_jump="
                << ( evolution_input_may_jump ? "TRUE  -> causal one-way t-interface ADMITTED"
                                              : "FALSE -> causal block SKIPPED; t-faces get"
                                                " ordinary two-way claims" )
                << "\n";
  }

  {
    // rev104: PLAN SHAPE restored.
    //
    // `n_elem@plan` and `evo_dom` were emitted by the rev89 gate-probe line, which rev101
    // removed wholesale with the rest of the predicate-(ii) reporting.  They are not
    // predicate-(ii) diagnostics: n_elem@plan=1 is the direct evidence that the evolution
    // grid COLLAPSED, and it is what verified Path B removes the seams from the plan
    // (OCFE_DAE5 4 -> 1 with claim census 6 -> 0; OCFE_MMPDE27 1 with 81 xi claims
    // surviving).  Losing it emptied the sweep's nelem column and made its H3 filter --
    // nelem_H0 == nelem_H1 -- match every applied driver vacuously.
    //
    // Removing a feature should not remove the measurement that judges it.
    static bool const kEvoFacePS = Options::_env_flag( "CRONOS_AUDIT_EVOFACE", false );
    if( kEvoFacePS ){
      size_t nel_evo = 0;
      if( _evolution_dom_set && _evolution_dom_var.dag() ){
        auto const it = _mDom.find( _evolution_dom_var );
        if( it != _mDom.cend() ) nel_evo = it->second.n_elem;
      }
      _disp( 3 ) << "  [evoface] plan shape: evo_dom="
                << ( _evolution_dom_set && _evolution_dom_var.dag()
                       ? _evolution_dom_var.name() : std::string("<none>") )
                << ( _evolution_dom_user ? " (user-supplied)" : " (auto-detected)" )
                << "  n_elem@plan=" << nel_evo
                << ( nel_evo == 1 ? "  <- COLLAPSED: no interior evolution seam exists" : "" )
                << "\n";
    }
  }
  {
    _evoProbe.ran           = true;
    _evoProbe.gate_current  = evolution_input_may_jump;
    _evoProbe.visits_by_dir.clear();

    static bool const kEvoFace = Options::_env_flag( "CRONOS_AUDIT_EVOFACE", false );
    if( kEvoFace ){
      if( !_evoProbe.computed )
        _disp( 3 ) << "  [evoface]   differential_states / with_own_IC_closure /"
                     " lumped_evo_coupling / n_elem@plan: NOT COMPUTED -- no evolution"
                     " direction is set, so these are undefined rather than zero."
                     "  (rev88 printed them as 0 here and the zeros were read as"
                     " measurements.)  _sgEqn@plan=" << _evoProbe.sgeqn_size << "\n";
      else
        _disp( 3 ) << "  [evoface]   differential_states=" << _evoProbe.n_dyn
                  << "  with_own_IC_closure=" << _evoProbe.n_transfers
                  << "  lumped_evo_coupling=" << _evoProbe.n_lumped_evo_couple
                  << "  grid_ok=" << ( _evoProbe.grid_ok ? "y" : "n" )
                  << "  _sgEqn@plan=" << _evoProbe.sgeqn_size
                  << "   [MEASURED]\n";
      // The guard divergence, stated rather than left to be inferred.  _sgEqn is populated
      // AFTER _collapse_evolution_for_marching runs and _reset() does not clear it, so the
      // lumped-coupling guard is vacuous at collapse time on a first setup.
      if( _evoProbe.sgeqn_size && _evoProbe.n_lumped_evo_couple )
        _disp( 3 ) << "  [evoface]   NOTE: " << _evoProbe.n_lumped_evo_couple
                  << " domain-less row(s) reduce over '" << _evolution_dom_var.name()
                  << "' -- refused here.  The SAME guard inside"
                     " _collapse_evolution_for_marching runs BEFORE _sgEqn is populated"
                     " (_check_consistency), and _reset() does not clear _sgEqn, so at"
                     " collapse time it is vacuous on a first setup and reads a released"
                     " DAG on a re-setup.  Reported, not fixed, in rev87.\n";
    }
  }

  auto auxiliary_algebraic_receiver_state = [&]
  ( FFVar const& eqnvar, FFModel::EqnOptions const& eqnopt,
    FFVar const& face_dom, size_t state_id ) -> bool
  {
    _AccT _tacc( _acc_aux, _cnt_aux );
    // Reduced-order auxiliary variables are ordinary collocated states.
    // Their defining/boundary rows often contain the auxiliary value
    // algebraically, not through a normal derivative, so the principal-symbol
    // coupling test can miss the C0 receiver needed to keep these generated
    // states continuous across finite-element interfaces.
    if( !_sopt( eqnopt ).receive_sat || !_dag ) return false;
    FFVar const* pst = state_by_id( state_id );
    if( !pst ) return false;
    // FIX 3: admit by provenance (reduce_order aux) OR by algebraic alias.
    // rev132: THE GATE, instrumented.  Traced as the origin of the 100-claim gap on
    // OCFE_PDE6_blk0 -- but traced, and seven traces have been refuted on this question,
    // so it reports rather than being believed.
    //
    // weak_sat_receiver_state calls this as an EARLY-TRUE, so a state that fails here can
    // still enter receiver_states by the principal-symbol or geometric-C0 routes; failing
    // this test is NOT by itself a skip.  The counters below therefore record the VERDICT
    // per state and direction, and the caller's outcome is visible in [sitestates].
    {
      static bool const kAux = Options::_env_flag( "CRONOS_AUDIT_AUXRECV", false );
      if( kAux ){
        bool const isaux = ( pst && _is_auxiliary_state( *pst ) );
        bool const isali = ( !isaux && alias_of_auxiliary( state_id, face_dom ) != nullptr );
        std::string dn = "dom?";
        for( auto const& kv : _mDom )
          if( kv.first.id().second == face_dom.id().second ){ dn = kv.first.name(); break; }
        static std::set<std::tuple<size_t,size_t,int>> seen;   // once per state/dir/verdict
        int const verdict = isaux ? 2 : ( isali ? 1 : 0 );
        auto const key = std::make_tuple( state_id,
                           static_cast<size_t>( face_dom.id().second ), verdict );
        if( seen.insert( key ).second )
          _disp( 3 ) << "  [auxrecv] state=" << std::left << std::setw(14)
                    << state_name_of_id( state_id ) << std::right
                    << " dir=" << std::setw(3) << dn
                    << "  _is_auxiliary_state=" << ( isaux ? "Y" : "n" )
                    << "  alias_of_auxiliary=" << ( isali ? "Y" : "n" )
                    << "  ->  " << ( ( isaux || isali )
                                       ? "eligible for the AUX-ALGEBRAIC receiver route"
                                       : "REFUSED here (may still enter via symbol/geo-C0)" )
                    << "\n";
      }
    }
    // 2026-10-15 (WORKPLAN 4R.21): a third admission route -- a participating row DECLARES this state as the one it
    // defines (EqnOptions::defines).  The alias resolver above anchors on a reduce_order() auxiliary and requires the
    // row to be affine in the state with exactly one other state in it; a hand-reduced model has no auxiliary to
    // anchor to, and a metric-weighted definition (q*x_xi = u_xi, g = M*x_xi) is not affine, so neither route admits
    // the moving-mesh states.  Admitting a declared state here builds its edges with the site's own geometry; the
    // post-pass before the receiver graph then keeps only the declared row's edge, as NATURAL, at a node where no
    // row differentiating the state across the face is collocated.
    // ONLY THROUGH ITS DEFINING ROW: admitting the state into every row it appears in sprays rescue edges and
    // changes the silence-and-restore balance of models the rescue already handles -- OCFE_zchain (p = g, q = g',
    // r = g'') went from 9 PASS to 3 with no promotion at all, purely from the extra edges.  One edge, from the row
    // that defines the state, is what the per-node matching needs and nothing else.
    bool declared_here = false;
    { auto const& dr = defrow_once; auto it = dr.find( static_cast<size_t>( pst->id().second ) );
      if( it != dr.end() && it->second < _mEqn.size() && _mEqn[ it->second ].var.id() == eqnvar.id() ){
        // A DECLARATION is an instruction and admits unconditionally.  An INFERENCE admits only the population
        // that actually lacks a receiver somewhere: the state IS differentiated in this face direction, and EVERY
        // row differentiating it is an evolution row (excluded at a causal element start).  q through the PDE and
        // g through RELAX qualify; a chain state differentiated by a LINK does not, and a state not differentiated
        // in this direction at all (OCFE_zchain's transition evaluations Evtz_*) does not -- admitting those added
        // 8 rescue terms the chain's silence-and-restore balance could not absorb (9 PASS -> 3).
        bool const declared = _mEqn[ it->second ].opt && _mEqn[ it->second ].opt->defines.dag()
                           && _mEqn[ it->second ].opt->defines.id() == pst->id();
        if( declared ) declared_here = true;
        else {
          size_t n_diff = 0, n_evo = 0;
          for( size_t ei = 0; ei < _mEqn.size(); ++ei ){
            if( !_mEqn[ei].opt->participate_in_classification ) continue;
            bool here = false;
            FFSubgraph const& sgq = eqn_subgraph( _mEqn[ei].var );
            for( auto const& op : sgq.l_op ){
              if( !op || !op->sameid( typeid( FFPartial ) ) || op->varin.empty() || !op->varin[0] ) continue;
              if( op->varin[0]->id() != pst->id() ) continue;
              auto const* pop = mc::type_cast<FFPartial const>( op );
              if( pop ) for( auto const& [dv, ordr] : pop->Indep().expr ){ (void)ordr; if( dv.id() == face_dom.id() ){ here = true; break; } }
              if( here ) break;
            }
            if( here ){ ++n_diff; if( _evoDerivRows.count( ei ) ) ++n_evo; }
          }
          declared_here = ( n_diff > 0 && n_evo == n_diff );
        }
      } }
    if( !declared_here && !_is_auxiliary_state( *pst )
     && alias_of_auxiliary( state_id, face_dom ) == nullptr ) return false;
    if( !state_depends_on_dom( state_id, face_dom ) ) return false;
    // rev155b (knob, default = pre-155b): CRONOS_AUX_EVO_CLAIMS=0 refuses this route in the
    // EVOLUTION direction for an auxiliary that is not differentiated in it.  Such an aux is
    // algebraic in t exactly like a hand-declared flux state (q in OCFE_PDE6_blk0), which the
    // principal-symbol test never admits as a t-receiver.  MEASURED (PDE6_blk0, IC_WEAK mono):
    // via this route Daux3_T receives 392 t-direction weak-SAT terms (prefactor 800 = T's) that
    // q never gets, and the corner dedup then removes 96 x-terms in their favour -> MANUAL and
    // AUTO are different weak discretisations (3.9x apart; 1.6x even at a=1 where q == Daux).
    // Under exact modes the aux t-claim is implied (T's t-trace at every transverse node + the
    // LINK collocated on both sides), so refusing it is a pure redundancy removal there.  The
    // v2 rationale above (PDE2 soft modes) concerns SPATIAL transverse directions only.
    {
      static int const kAuxEvo = 1;   // CRONOS_AUX_EVO_CLAIMS (retired 2026-10-07, WORKPLAN 3.B batch 2b)
      if( !kAuxEvo && _evolution_dom_set && _evolution_dom_var.dag()
       && face_dom.id() == _evolution_dom_var.id() && _is_auxiliary_state( *pst ) ){
        bool diff_in_evo = false;
        for( auto const& a : _auxDef )
          if( a.aux.id() == pst->id() ){ diff_in_evo = a.diff_dom.count( face_dom ) > 0; break; }
        if( !diff_in_evo ) return false;
      }
    }
    // rev155e (knob, default = pre-155e): CRONOS_AUX_OWNDIR_LINK_RECV=0 refuses the aux-algebraic
    // route for an auxiliary A = d(P)/d(e) in its OWN direction e when the row is a LINK row.
    // MEASURED (PDE3, CRONOS_PLAN_DUMP_DETCOLS): P's C0 z-claim tau and A's z-claim tau land in
    // P's z-LINK row with IDENTICAL entries (+-108 at the same flats); A's differs only by tiny
    // PDE-row entries.  Two multipliers in one row -> collinear by construction: the residual
    // k'=12 (PDE2/PDE3/MBC1 STRONG) and, with a cross-aux tau absorbing the PDE-row entries, the
    // defR=4 that makes Sec.3.5 load-bearing on PDE3.  A's own-direction claim keeps its principal
    // receivers (rows where d(A)/d(e) appears).  The alias route (declared q) is NOT touched.
    {
      static int const kOwn = 1;   // CRONOS_AUX_OWNDIR_LINK_RECV (retired 2026-10-07, WORKPLAN 3.B batch 2b)
      if( !kOwn && eqnopt.role == EqnRole::LINK && _is_auxiliary_state( *pst ) ){
        for( auto const& a : _auxDef )
          if( a.aux.id() == pst->id() ){ if( a.diff_dom.count( face_dom ) ) return false; break; }
      }
    }
    // A reduced auxiliary A = d(P)/d(e) only needs a continuity claim in its
    // own differentiation direction e (the genuine C1/flux continuity).  In
    // any other direction A's continuity follows from the primitive P's
    // continuity plus the defining LINK, so a claim here is redundant and
    // contends with P's geometric/principal claim for the same row.
    //
    // v2: the "follows from primitive + LINK" assumption holds only for the
    // penalized (IC_WEAK) system, whose weak coupling incidentally pins the
    // transverse aux.  The exact modes remove that coupling, leaving the
    // transverse aux genuinely null (PDE2 soft modes Dz_*/Drl_* at the exact
    // solution).  Under v2 we therefore GENERATE the transverse aux value-C0
    // edge (coupling=1 via the fallback at edge construction) and defer
    // redundancy handling to the structural rank-prune (Phase 2), instead of
    // suppressing it here by direction.

    if( eqnvar.id() == pst->id() ) return true;
    FFSubgraph const& sg = eqn_subgraph( eqnvar );
    for( auto const& op : sg.l_op ){
      for( auto const* vin : op->varin )
        if( vin && vin->id() == pst->id() ) return true;
      for( auto const* vout : op->varout )
        if( vout && vout->id() == pst->id() ) return true;
    }
    return false;
  };

  auto weak_sat_receiver_state = [&]
  ( size_t row_id, FFVar const& eqnvar, FFModel::EqnOptions const& eqnopt,
    std::map<FFVar,int,lt_FFVar> const& eqndom,
    FFVar const& face_dom, FFDom const& face_pdom,
    size_t iel, size_t inode, size_t state_id ) -> bool
  {
    // Exact strong/trace rows may only consume or target rows on which the
    // corresponding IC_WEAK discretisation would have placed the continuity
    // SAT term for this state and this finite-element face.  This is the
    // primary admissibility test; donor classes are only tie-breakers later.
    if( !_sopt( eqnopt ).receive_sat ) return false;
    if( inode >= face_pdom.n_node ) return false;
    if( !state_depends_on_dom( state_id, face_dom ) ) return false;

    size_t const ndx = row_id;

    bool const on_left_receiver  = ( iel > 0 && inode == 0 );
    bool const on_right_receiver = ( iel + 1 < face_pdom.n_elem
                                  && inode + 1 == face_pdom.n_node );
    if( !on_left_receiver && !on_right_receiver ) return false;

    auto itdom = eqndom.find( face_dom );
    if( itdom != eqndom.end()
     && ( itdom->second == FFDom::LB || itdom->second == FFDom::UB ) )
      return false;

    Options::InterfaceType const sat_type =
      _resolve_interface_type( &eqnopt, &face_dom, true, eqnopt.block_id );
    if( sat_type != Options::IC_VALUE
     && sat_type != Options::IC_FLUX
     && sat_type != Options::IC_UPWIND )
      return false;


    if( sat_type != Options::IC_UPWIND
     && auxiliary_algebraic_receiver_state( eqnvar, eqnopt, face_dom, state_id ) )
      return true;

    if( sat_type == Options::IC_UPWIND ){
      if( !_classified ) return false;
      t_Symbol const& sym = _symbol_for_block( eqnopt.block_id );
      auto const& fdata = _face_data_for_block( eqnopt.block_id );
      t_FaceData const* fd = _face_data_for_dom( fdata, &face_dom );
      if( !fd ) return false;

      size_t k_eq = sym.vEqn.size();
      for( size_t k = 0; k < sym.vEqn.size(); ++k )
        if( sym.vEqn[k].id() == eqnvar.id() ){ k_eq = k; break; }
      if( k_eq == sym.vEqn.size() ) return false;

      size_t j_st = sym.vState.size();
      for( size_t j = 0; j < sym.vState.size(); ++j )
        if( static_cast<size_t>( sym.vState[j].id().second ) == state_id ){ j_st = j; break; }
      if( j_st == sym.vState.size() ) return false;

      if( on_left_receiver
       && (arma::uword)k_eq < fd->Aplus.n_rows
       && (arma::uword)j_st < fd->Aplus.n_cols
       && std::fabs( fd->Aplus( (arma::uword)k_eq, (arma::uword)j_st ) )
          > 64. * DBL_EPSILON )
        return true;
      if( on_right_receiver
       && (arma::uword)k_eq < fd->Aminus.n_rows
       && (arma::uword)j_st < fd->Aminus.n_cols
       && std::fabs( fd->Aminus( (arma::uword)k_eq, (arma::uword)j_st ) )
          > 64. * DBL_EPSILON )
        return true;
      return false;
    }

    auto itstate = _mDep.find( static_cast<long>( state_id ) );
    if( itstate == _mDep.end() || itstate->second.first != STATE ) return false;
    size_t const pos_s = itstate->second.second;

    double const coupling = _sat_coupling( ndx, pos_s, face_dom );
    return std::fabs( coupling ) > 64. * DBL_EPSILON;
  };

  auto weak_sat_receiver_state_ids = [&]
  ( size_t row_id, FFVar const& eqnvar, FFModel::EqnOptions const& eqnopt,
    std::map<FFVar,int,lt_FFVar> const& eqndom,
    FFVar const& face_dom, FFDom const& face_pdom,
    size_t iel, size_t inode )
  {
    _AccT _tacc( _acc_weak, _cnt_weak );
    std::vector<size_t> ids_out;
    for( auto const& [st, dep] : _mVar ){
      (void)dep;
      size_t const sid = static_cast<size_t>( st.id().second );
      if( weak_sat_receiver_state( row_id, eqnvar, eqnopt, eqndom, face_dom,
                                   face_pdom, iel, inode, sid ) )
        push_unique_state_id( ids_out, sid );
    }
    return ids_out;
  };

  //! Geometric-C0 receiver predicate.
  //!
  //! Fires only for primitive collocated states that are geometrically
  //! duplicated across an internal element interface but have zero
  //! principal-symbol coupling AND are not auxiliary-algebraic receivers in
  //! that direction.  The calling code (Step 6) must invoke this ONLY for
  //! states not already found by weak_sat_receiver_state_ids, so it is
  //! structurally a third-priority fallback.
  //!
  //! Guards:
  //!   1. Internal element face (not a physical domain boundary)
  //!   2. True collocated primitive state (in _mVar, not auxiliary/input/const)
  //!   3. State depends on face_dom (geometrically duplicated)
  //!   4. Equation has receive_sat enabled
  //!   5. Equation owns the state through a tangential auxiliary: it
  //!      references a reduced auxiliary whose primitive root is the state,
  //!      and that auxiliary's derivative direction differs from face_dom.
  //!      This captures the degenerate membrane case (e.g. equation for
  //!      Cd(rd,z) references Drd_Cd, derivative in rd, face in z) while
  //!      excluding ordinary PDE rows where the state merely appears in
  //!      reaction/source/convection terms.
  auto geometric_c0_receiver_state = [&]
  ( size_t /*row_id*/, FFVar const& eqnvar, FFModel::EqnOptions const& eqnopt,
    std::map<FFVar,int,lt_FFVar> const& eqndom,
    FFVar const& face_dom, FFDom const& face_pdom,
    size_t iel, size_t inode, size_t state_id ) -> bool
  {
    _AccT _tacc( _acc_geom, _cnt_geom );
    // Guard 4: equation must accept SAT penalties.
    if( !_sopt( eqnopt ).receive_sat ) return false;

    // Guard 4b: a geometric C0 claim is realised by an exact replacement that
    // DROPS the row it owns.  Only a redundant interior duplicate may be
    // dropped.  Dropping a LINK row (Partial[e](P) - A = 0) would unpin the
    // auxiliary A it defines, so A's value would float in the degenerate
    // direction (the spurious flux drift seen when geometric edges owned LINK
    // rows); dropping an INTERFACE/BOUNDARY row would discard unique physics
    // (e.g. the dry/wet flux match).  With LINK rows kept intact on both
    // copies, primitive C0 continuity propagates to the auxiliaries through
    // the LINK exactly as the auxiliary-direction gate assumes.
    if( eqnopt.role != EqnRole::INTERIOR ) return false;

    // Guard 1a: node must be at an internal element interface.
    if( inode >= face_pdom.n_node ) return false;
    bool const on_left  = ( iel > 0 && inode == 0 );
    bool const on_right = ( iel + 1 < face_pdom.n_elem
                         && inode + 1 == face_pdom.n_node );
    if( !on_left && !on_right ) return false;

    // Guard 1b: equation domain mask must not restrict this face to a
    // physical boundary (LB or UB).  Internal faces have ALL or ALL-LB/UB.
    auto itdom = eqndom.find( face_dom );
    if( itdom != eqndom.end()
     && ( itdom->second == FFDom::LB || itdom->second == FFDom::UB ) )
      return false;

    // Guard 2+3: state must be a primitive collocated field depending on face_dom.
    FFVar const* pst = state_by_id( state_id );
    if( !pst ) return false;
    if( _mVar.find( *pst ) == _mVar.end() ) return false;
    if( _is_auxiliary_state( *pst ) ) return false;
    if( !state_depends_on_dom( state_id, face_dom ) ) return false;

    // Guard 5: tangential auxiliary-mediated ownership.
    //
    // The equation must reference a reduced auxiliary whose primitive root
    // is pst, and whose derivative direction is NOT face_dom.  This captures
    // degenerate tensor-product membrane equations (Cd in z-direction when
    // the equation only differentiates in rd) while excluding ordinary PDE
    // rows where the primitive state merely appears in reaction/source terms.
    if( !_dag ) return false;
    FFSubgraph const& sg = eqn_subgraph( eqnvar );

    // Guard 5a: the host equation must be TANGENTIAL to the interface, i.e. it
    // must not differentiate anything in the face direction.  A face-direction
    // (normal) derivative makes the two element copies of the row at a shared
    // interface node genuinely different one-sided stencils, so dropping one
    // removes real normal-flux information rather than a redundant duplicate.
    // Guard 5 below only inspects the continuity state's AUXILIARY direction,
    // so it cannot see a normal derivative of a PRIMITIVE state -- e.g. the gas
    // Vg ODE  Partial[z](Vg) + Drd_Cd*c  references Drd_Cd undifferentiated and
    // would pass Guard 5, letting a Cd z-continuity claim drop the Vg evolution
    // row (the PDE3 soft-mode steal: Vg loses an equation at the rd-boundary
    // corner where no tangential host exists).  Requiring a tangential host
    // excludes such rows; a claim with no tangential host emits no geometric-C0
    // receiver edge and is left to the redundant physics already present
    // (e.g. the hard Cd=Cg interface row plus Cg's own continuity), keeping the
    // collocation system square because the would-be-dropped row is retained.
    for( auto const& op : sg.l_op ){
      if( !op || !op->sameid( typeid(FFPartial) ) ) continue;
      auto const* pop = mc::type_cast<FFPartial const>( op );
      if( !pop ) continue;
      for( auto const& [dv,ord] : pop->Indep().expr )
        if( ord && dv.id() == face_dom.id() ) return false;
    }

    for( auto const& [aux, auxdoms] : _mVar ){
      if( !_is_auxiliary_state( aux ) ) continue;
      FFVar const* root = _auxiliary_primitive_root( aux );
      if( !root || root->id() != pst->id() ) continue;

      // Does the equation reference this auxiliary?
      bool eqn_has_aux = false;
      for( auto const& op : sg.l_op ){
        for( auto const* vin : op->varin )
          if( vin && vin->id() == aux.id() ){ eqn_has_aux = true; break; }
        if( eqn_has_aux ) break;
      }
      if( !eqn_has_aux ) continue;

      // The auxiliary must NOT be differentiated in the face direction.
      // If it is, this is an ordinary reduced flux/auxiliary continuity
      // case, not a degenerate tangential primitive-C0 fallback.
      bool aux_diff_in_face_dom = false;
      for( auto const& a : _auxDef ){
        if( a.aux.id() != aux.id() ) continue;
        for( auto const& dd : a.diff_dom )
          if( dd.id() == face_dom.id() ){ aux_diff_in_face_dom = true; break; }
        break;
      }
      if( aux_diff_in_face_dom ) continue;

      // Equation owns pst through a tangential auxiliary.
      return true;
    }

    return false;
  };

  auto make_state_claim_key_from_eqn_context = [&]
  ( int block_id, size_t state_id, size_t eqn_id, FFVar const& face_dom,
    size_t iel, size_t inode, size_t n_node, size_t n_elem,
    auto const& dims, std::vector<size_t> const& midx,
    int claim_kind, t_InterfaceClaimKey& key ) -> bool
  {
    // Physical boundary/interface donor equations may involve states defined
    // on different tensor-product domains, e.g. DW_VAL: Cw(rw,z)-H*Cd(rd,z).
    // The donor row itself is indexed by the boundary equation dimensions
    // {rd,rw,z}, but LOOP 1 emits the exact C0 row in the candidate state's
    // own coefficient space.  Therefore the emitted claim key must encode the
    // candidate state's tangential coordinates, not the full boundary-equation
    // tensor face.  Otherwise Cw@rw=0 can be mis-keyed as if it lived at
    // Cd@rd=1, and the acceptor row is silently skipped.
    if( n_node <= 1 ) return false;

    size_t iel_lo = 0;
    if( iel > 0 && inode == 0 )
      iel_lo = iel - 1;
    else if( iel + 1 < n_elem && inode + 1 == n_node )
      iel_lo = iel;
    else
      return false;

    FFVar const* pst = state_by_id( state_id );
    if( !pst ) return false;
    auto itdep = _mVar.find( *pst );
    if( itdep == _mVar.end() ) return false;

    size_t face = 0;
    size_t mult = 1;
    for( auto const& sdom : itdep->second ){
      if( sdom.id() == face_dom.id() ) continue;

      size_t q = std::numeric_limits<size_t>::max();
      for( size_t d = 0; d < dims.size(); ++d ){
        if( !dims[d].pvar || !dims[d].pdom ) continue;
        if( dims[d].pvar->id() == sdom.id() ){ q = d; break; }
      }
      if( q == std::numeric_limits<size_t>::max() ) return false;

      size_t ed = dims[q].iel;
      size_t kd = 0;
      switch( dims[q].lim ){
        case FFDom::LB:
          ed = 0;
          kd = 0;
          break;
        case FFDom::UB:
          ed = dims[q].pdom->n_elem - 1;
          kd = dims[q].pdom->n_node - 1;
          break;
        default:
          kd = midx[q];
          break;
      }

      size_t const radix = dims[q].pdom->n_elem * dims[q].pdom->n_node;
      face += ( ed * dims[q].pdom->n_node + kd ) * mult;
      mult *= radix;
    }

    key = { block_id, state_id, eqn_id, static_cast<size_t>( face_dom.id().second ),
            iel_lo, face, claim_kind };
    return true;
  };

  auto build_trace_claims_from_weak_sat_graph =
  [&]( std::map<t_InterfaceClaimKey,size_t>& claims )
  {
    claims.clear();
    draft.receiver_edges.clear();
    draft.weak_sat_terms.clear();
    draft.weak_sat_terms_ready = false;

    // IC_TRACE does not impose every raw weak-SAT receiver edge.  At
    // tensor-product FE corners, the raw graph contains cycles between the
    // same duplicate physical trace DOFs.  We therefore collect all raw edges,
    // group their endpoint coefficient copies by physical state/coordinate, and
    // retain only a spanning forest per duplicate group.  The retained edge key
    // is still the ordinary exact-row claim key, so LOOP 1 and the receiver
    // injection code can consume the frozen plan without rediscovering graph
    // topology during eval().
    struct t_TraceVertexKey {
      size_t state_id = 0;
      std::vector<std::tuple<size_t,size_t,size_t>> dof; // (dom id, element, local node)
      bool operator<( t_TraceVertexKey const& o ) const
      {
        if( state_id != o.state_id ) return state_id < o.state_id;
        return dof < o.dof;
      }
    };

    struct t_TraceGroupKey {
      size_t state_id = 0;
      int    kind     = CLAIM_EXACT_C0;
      std::vector<std::pair<size_t,long long>> coord; // (dom id, rounded physical coordinate)
      bool operator<( t_TraceGroupKey const& o ) const
      {
        if( state_id != o.state_id ) return state_id < o.state_id;
        if( kind     != o.kind     ) return kind     < o.kind;
        return coord < o.coord;
      }
    };

    struct t_TraceRawEdge {
      t_InterfaceClaimKey claim;
      t_TraceVertexKey    vlo;
      t_TraceVertexKey    vhi;
    };

    std::map<t_TraceGroupKey,std::vector<t_TraceRawEdge>> raw_edges;
    std::set<t_InterfaceClaimKey> raw_claim_seen;

    auto rounded_coord = [] ( FFDom const& dom, size_t iel, size_t inode ) -> long long
    {
      if( iel >= dom.n_elem || inode >= dom.n_node ) return 0;
      auto const nodes = dom.lgnodes( dom.elem_lo( iel ), dom.elem_up( iel ) );
      double const scale = std::max( 1.0, std::fabs( dom.up_dom - dom.lo_dom ) );
      return static_cast<long long>( std::llround( nodes.at( inode ) / scale * 1.0e12 ) );
    };

    auto build_trace_edge_keys = [&]
    ( t_InterfaceClaimKey const& key, FFVar const& face_dom,
      auto const& dims, std::vector<size_t> const& midx,
      t_TraceGroupKey& group, t_TraceVertexKey& vlo,
      t_TraceVertexKey& vhi ) -> bool
    {
      FFVar const* pst = state_by_id( key.state_id );
      if( !pst ) return false;
      auto itdep = _mVar.find( *pst );
      if( itdep == _mVar.end() ) return false;

      group.state_id = key.state_id;
      group.kind     = key.kind;
      group.coord.clear();
      vlo.state_id = key.state_id;
      vhi.state_id = key.state_id;
      vlo.dof.clear();
      vhi.dof.clear();

      for( auto const& sdom : itdep->second ){
        auto itdom = _mDom.find( sdom );
        if( itdom == _mDom.end() ) return false;
        FFDom const& pdom = itdom->second;
        size_t const dom_id = static_cast<size_t>( sdom.id().second );

        if( sdom.id() == face_dom.id() ){
          if( key.iel_lo + 1 >= pdom.n_elem || pdom.n_node == 0 ) return false;
          size_t const klo = pdom.n_node - 1;
          size_t const khi = 0;
          vlo.dof.emplace_back( dom_id, key.iel_lo,     klo );
          vhi.dof.emplace_back( dom_id, key.iel_lo + 1, khi );
          group.coord.emplace_back( dom_id, rounded_coord( pdom, key.iel_lo, klo ) );
          continue;
        }

        size_t q = std::numeric_limits<size_t>::max();
        for( size_t d = 0; d < dims.size(); ++d ){
          if( !dims[d].pvar || !dims[d].pdom ) continue;
          if( dims[d].pvar->id() == sdom.id() ){ q = d; break; }
        }
        if( q == std::numeric_limits<size_t>::max() ) return false;

        size_t ed = dims[q].iel;
        size_t kd = 0;
        switch( dims[q].lim ){
          case FFDom::LB:
            ed = 0;
            kd = 0;
            break;
          case FFDom::UB:
            ed = dims[q].pdom->n_elem - 1;
            kd = dims[q].pdom->n_node - 1;
            break;
          default:
            kd = midx[q];
            break;
        }
        if( ed >= pdom.n_elem || kd >= pdom.n_node ) return false;

        vlo.dof.emplace_back( dom_id, ed, kd );
        vhi.dof.emplace_back( dom_id, ed, kd );
        group.coord.emplace_back( dom_id, rounded_coord( pdom, ed, kd ) );
      }

      return !vlo.dof.empty() && !vhi.dof.empty();
    };

    auto add_raw_trace_edge = [&]
    ( t_InterfaceClaimKey key, FFVar const& face_dom,
      auto const& dims, std::vector<size_t> const& midx )
    {
      // Canonicalise to the physical trace slot used by LOOP 1 and receiver
      // lookups.  Multiple receiver equations may discover the same raw edge;
      // keep only one graph edge/unknown for that physical slot.
      key.block_id = _trace_claim_block( key.state_id, key.dom_id );
      key.eqn_id   = invalid_claim_id;
      if( !raw_claim_seen.insert( key ).second ) return;

      t_TraceGroupKey group;
      t_TraceVertexKey vlo, vhi;
      if( build_trace_edge_keys( key, face_dom, dims, midx, group, vlo, vhi ) ){
        if( !( vlo < vhi ) && !( vhi < vlo ) ) return; // ignore degenerate self loops
        raw_edges[group].push_back( { key, vlo, vhi } );
      }
      else{
        // Conservative fallback for unusual legacy equations whose receiver
        // state/equation tensor domains cannot be reconciled.  Keep the edge as
        // an independent singleton rather than dropping a continuity condition.
        group.state_id = key.state_id;
        group.kind     = key.kind;
        group.coord    = { { key.dom_id, static_cast<long long>( key.iel_lo ) },
                           { static_cast<size_t>(-1), static_cast<long long>( key.face ) } };
        vlo.state_id = vhi.state_id = key.state_id;
        vlo.dof = { { key.dom_id, key.iel_lo,     key.face } };
        vhi.dof = { { key.dom_id, key.iel_lo + 1, key.face } };
        raw_edges[group].push_back( { key, vlo, vhi } );
      }
    };

    [[maybe_unused]] auto trace_hard_physical_row = []( FFModel::EqnOptions const& opt ) -> bool
    {
      switch( opt.role ){
        case EqnRole::INTERIOR:
        case EqnRole::SURFACE:
        case EqnRole::LINK:
        case EqnRole::AUTO:
          return false;
        case EqnRole::INITIAL:
        case EqnRole::BOUNDARY:
        case EqnRole::INTERFACE:  // user transmission rows are unique physical
                                  // conditions, exactly like INITIAL/BOUNDARY:
                                  // weak-only, never a free-multiplier receiver.
        case EqnRole::DIAGNOSTIC:
        default:
          return true;
      }
    };

    [[maybe_unused]] auto trace_tangential_physical_boundary = [&]
    ( size_t state_id, size_t face_dim, auto const& dims,
      std::vector<size_t> const& ) -> bool
    {
      for( size_t q = 0; q < dims.size(); ++q ){
        if( q == face_dim ) continue;
        if( !dims[q].pdom || !dims[q].pvar ) continue;
        // A tangential hard physical boundary is present only when the
        // equation domain is explicitly collapsed to a face of a tangential
        // domain that is actually part of the candidate receiver state.  Do
        // not classify auxiliary/coupling coordinates as hard tangential
        // boundaries for lower-dimensional states such as Cg(z) or Vg(z)
        // embedded in a gas/liquid equation written on {z,rl} with rl=LB.
        if( !state_depends_on_dom( state_id, *dims[q].pvar ) ) continue;
        if( dims[q].lim == FFDom::LB || dims[q].lim == FFDom::UB )
          return true;
      }
      return false;
    };

    // Deferred geometric-C0 candidates.
    //
    // The geometric-C0 fallback must be a GLOBAL decision: a (state, dom, iel_lo)
    // combination needs geometric-C0 only if NO equation row provides a
    // principal-symbol or auxiliary-algebraic receiver edge for it.  Because the
    // emission loop processes rows sequentially, a per-row check cannot know
    // whether a later row will provide coverage.  Therefore geometric edges are
    // collected during the loop but deferred; after the loop completes, they are
    // filtered against a global coverage set and only the truly uncovered
    // survivors are emitted.
    struct DeferredGeoCandidate {
      InterfaceReceiverEdge edge;
      t_InterfaceClaimKey   claim_key;
      t_TraceGroupKey       group;
      t_TraceVertexKey      vlo;
      t_TraceVertexKey      vhi;
      bool                  has_trace_keys = false;
    };
    std::vector<DeferredGeoCandidate> deferred_geo;

    // Coverage set: (state_id, dom_id, iel_lo) triples that have at least one
    // principal-symbol or auxiliary-algebraic receiver edge from any equation row.
    std::set<std::tuple<size_t,size_t,size_t>> covered_state_face;

    size_t n_receiver_edges = 0;
    size_t n_trace_edges_allowed = 0;
    size_t n_weak_only_edges = 0;
    size_t n_hard_row_edges = 0;
    size_t n_hard_tangent_edges = 0;
    size_t n_geometric_c0_edges = 0;
    size_t n_geometric_c0_discarded_covered = 0;

#ifdef MC__OCFESLV_CORNER_SUPPRESS_PROBE
    // item-PDE5 probe (default-off): measure reduced-aux tensor-corner
    // suppression per (state,direction).  Tests the >2-direction hypothesis:
    // the ordered single-owner corner dedup in _suppress_tensor_interface_corner
    // defers to dimensions EARLIER in the ordering, which spans a 2D corner but
    // is expected to over-fire at 3D edges/corners (>=2 earlier dims firing).
    struct _CornerProbeRec {
      size_t emitted = 0, suppressed = 0;
      std::map<size_t,size_t> fired_hist;   // #earlier-dims that triggered
      std::map<size_t,size_t> mult_hist;    // duplicate node multiplicity (2^k)
    };
    std::map<std::pair<size_t,size_t>,_CornerProbeRec> _corner_probe;

    // Connectivity (option-a): prove + size the float.  For each state we record
    // the SURVIVING value-C0 ties (recv copy <-> donor copy across an interior
    // interface) and the per-dim mesh; at the end we union-find the copies of
    // every duplicate node and count connected components.  components>1 == a
    // genuine float (the dedup removed too many ties).  copy is encoded per-dim
    // as iel*n_node + inode.
    struct _ConnDim { size_t dom_id = 0, n_node = 0, n_elem = 0; };
    struct _ConnState {
      std::vector<_ConnDim> dims;
      bool have_dims = false;
      std::vector<std::pair<std::vector<size_t>,std::vector<size_t>>> ties;
    };
    std::map<size_t,_ConnState> _conn;
#endif


    // The states with an evolution derivative: pinned at an element start by their continuity rows (the pin rule).
    std::set<size_t> evo_diff_ids;
    if( _causal_rows_active() )
      for( auto const& st : _states_with_evolution_derivative( _mEqn, _mVar ) )
        evo_diff_ids.insert( static_cast<size_t>( st.id().second ) );

    size_t receiver_row_id = 0;
    for( auto const& eqnrec : _mEqn ){
      size_t const row_id = receiver_row_id++;
      auto const& eqnvar = eqnrec.var;
      auto const& eqndom = eqnrec.dom;
      auto const& eqnopt = *eqnrec.opt;
      if( !_sopt( eqnopt ).receive_sat ) continue;

      std::map<FFVar, size_t, lt_FFVar> ndx_el;
      for( ; _next_block( ndx_el, eqndom ) && !ndx_el.empty(); ){
        struct TraceClaimDim {
          FFDom const* pdom = nullptr;
          FFVar const* pvar = nullptr;
          size_t       cur_n = 0;
          int          lim = FFDom::ALL;
          size_t       iel = 0;
        };

        std::vector<TraceClaimDim> cdims;
        cdims.reserve( eqndom.size() );
        for( auto const& [pvar, lim] : eqndom ){
          auto itdom = _mDom.find( pvar );
          if( itdom == _mDom.end() ) continue;
          size_t const iel = ndx_el.count(pvar) ? ndx_el.at(pvar) : 0;
          size_t const cur_n = ( lim == FFDom::LB || lim == FFDom::UB )
                             ? 1 : itdom->second.n_node;
          cdims.push_back({ &itdom->second, &pvar, cur_n, lim, iel });
        }

        size_t ntot = 1;
        for( auto const& cd : cdims ) ntot *= cd.cur_n;

        std::vector<size_t> midx( cdims.size(), 0 );
        for( size_t flat = 0; flat < ntot; ++flat ){
          size_t rem = flat;
          for( size_t d = 0; d < cdims.size(); ++d ){
            midx[d] = rem % cdims[d].cur_n;
            rem    /= cdims[d].cur_n;
          }

          bool base_keep = true;
          for( size_t d = 0; d < cdims.size() && base_keep; ++d ){
            auto const& cd = cdims[d];
            if( cd.lim == FFDom::LB || cd.lim == FFDom::UB ) continue;
            switch( cd.lim ){
              case FFDom::ALL:
                break;
              case FFDom::ALL-FFDom::LB:
                if( cd.iel == 0 && midx[d] == 0 ) base_keep = false;
                if( _causal_start_excludes( static_cast<size_t>( &eqnrec - &_mEqn[0] ), *cd.pvar, cd.iel, midx[d] ) ) base_keep = false;   // row-level causal design
                break;
              case FFDom::ALL-FFDom::UB:
                if( cd.iel + 1 == cd.pdom->n_elem
                 && midx[d] + 1 == cd.pdom->n_node ) base_keep = false;
                break;
              case FFDom::ALL-FFDom::LB-FFDom::UB:
                if( cd.iel == 0 && midx[d] == 0 ) base_keep = false;
                if( _causal_start_excludes( static_cast<size_t>( &eqnrec - &_mEqn[0] ), *cd.pvar, cd.iel, midx[d] ) ) base_keep = false;
                if( cd.iel + 1 == cd.pdom->n_elem
                 && midx[d] + 1 == cd.pdom->n_node ) base_keep = false;
                break;
              default:
                throw Exceptions( Exceptions::DOMAIN );
            }
          }
          if( !base_keep ) continue;

          for( size_t d = 0; d < cdims.size(); ++d ){
            auto const& cd = cdims[d];
            if( cd.lim == FFDom::LB || cd.lim == FFDom::UB ) continue;
            // 2026-10-02: under the row-level causal design the evolution faces carry no claims at all -- each
            // element start is an initial node with its own continuity rows (as a marching window's start)
            if( _causal_rows_active() && cd.pvar->id() == _evolution_dom_var.id() ) continue;
            _AccT _tpd( _acc_perdim, _cnt_perdim );

            bool keep;
            { _AccT _tk( _acc_keep, _cnt_keep );
              keep = _keep_eqn_node
                ( eqnvar, eqnopt, *cd.pvar, *cd.pdom, cd.iel, midx[d], cd.lim, true ); }
            if( !keep ) continue;

            bool const normal_internal_face =
                 ( cd.iel > 0 && midx[d] == 0 )
              || ( cd.iel + 1 < cd.pdom->n_elem && midx[d] + 1 == cd.pdom->n_node );
            if( !normal_internal_face ) continue;

            Options::InterfaceType dir_type;
            { _AccT _tr( _acc_rtype, _cnt_rtype );
              dir_type = _resolve_interface_type( &eqnopt, cd.pvar, true, eqnopt.block_id ); }
            if( dir_type != Options::IC_VALUE
             && dir_type != Options::IC_FLUX
             && dir_type != Options::IC_UPWIND )
              continue;


            auto receiver_states = weak_sat_receiver_state_ids
              ( row_id, eqnvar, eqnopt, eqndom, *cd.pvar, *cd.pdom, cd.iel, midx[d] );


            // Geometric-C0 augmentation: find collocated states that are
            // geometrically duplicated across this face but have no
            // principal-symbol or auxiliary-algebraic coupling.  These are
            // third-priority fallback receiver states.  Each eligible equation
            // row still needs its own receiver edge so that IC_TRACE has a
            // multiplier injection slot in the B-matrix; do not deduplicate
            // receiver edges across equation rows for the same trace claim.
            std::set<size_t> geo_c0_states;
            for( auto const& [st, dep] : _mVar ){
              (void)dep;
              size_t const sid = static_cast<size_t>( st.id().second );
              if( std::find( receiver_states.begin(), receiver_states.end(), sid )
                  != receiver_states.end() )
                continue;
              if( geometric_c0_receiver_state( row_id, eqnvar, eqnopt, eqndom,
                                               *cd.pvar, *cd.pdom, cd.iel, midx[d], sid ) ){
                push_unique_state_id( receiver_states, sid );
                geo_c0_states.insert( sid );
              }
            }
            // ---- FIX 1: do not ISSUE a C0 claim that duplicates a pointwise pin ------
            // Measured defect (tier0-audit, chain 2, nel=2): rank_deficiency = 1, whose
            // left-null vector is exactly {1 claim row, 2 INITIAL rows} with weights
            // +/-1/sqrt(3).  An INITIAL/BOUNDARY equation whose mask is ALL in THIS
            // direction pins the state at BOTH copies of an interior face node, so the
            // C0 claim across that face is the difference of the two pin rows.
            //   It must be suppressed at issue, not dropped afterwards:
            // _detect_redundant_continuity_claims is a COLUMN-dependency test and is
            // structurally blind to a ROW dependency of this shape.
            {
              auto in_mask = [&]( size_t dd, int m ) -> bool {
                auto const& c2 = cdims[dd];
                bool const isLB = ( c2.iel == 0 && midx[dd] == 0 );
                bool const isUB = ( c2.iel + 1 == c2.pdom->n_elem
                                 && midx[dd] + 1 == c2.pdom->n_node );
                if( m == FFDom::ALL ) return true;
                if( m == FFDom::LB  ) return isLB;
                if( m == FFDom::UB  ) return isUB;
                if( m == FFDom::ALL - FFDom::LB )             return !isLB;
                if( m == FFDom::ALL - FFDom::UB )             return !isUB;
                if( m == FFDom::ALL - FFDom::LB - FFDom::UB ) return !isLB && !isUB;
                return true;                       // unrecognised mask: keep the claim
              };
              auto pinned_pointwise = [&]( size_t sid ) -> bool {
                // 2026-10-02: at a causal element start (row-level design) a state with an evolution derivative
                // is pinned pointwise by its continuity row, as by an INITIAL row at t = 0
                if( _causal_rows_active() && evo_diff_ids.count( sid ) )
                  for( size_t e = 0; e < cdims.size(); ++e )
                    if( e != d && _causal_start_node( *cdims[e].pvar, cdims[e].iel, midx[e] ) ) return true;
                for( auto const& er : _mEqn ){
                  if( er.opt->role != EqnRole::INITIAL && er.opt->role != EqnRole::BOUNDARY )
                    continue;
                  auto itm = er.dom.find( *cd.pvar );
                  if( itm == er.dom.end() || itm->second != FFDom::ALL ) continue;
                  bool here = true;                // ... and this node is inside its domain
                  for( size_t d2 = 0; d2 < cdims.size() && here; ++d2 ){
                    if( d2 == d ) continue;
                    auto it2 = er.dom.find( *cdims[d2].pvar );
                    if( it2 == er.dom.end() || !in_mask( d2, it2->second ) ) here = false;
                  }
                  if( !here ) continue;
                  FFVar expr = er.var;                        // ... and it pins ONLY sid
                  FFGraph* dg = static_cast<FFGraph*>( expr.dag() );
                  if( !dg ) continue;
                  FFSubgraph sg2 = dg->subgraph( 1, &expr );
                  size_t nst = 0; bool has = false;
                  for( auto const* op : sg2.l_op ){
                    if( !op || op->type != FFOp::VAR || !op->varout[0] ) continue;
                    if( _mVar.find( *op->varout[0] ) == _mVar.end() ) continue;
                    ++nst;
                    if( static_cast<size_t>( op->varout[0]->id().second ) == sid ) has = true;
                  }
                  if( nst == 1 && has ) return true;
                }
                return false;
              };

              // ---- FIX 3 (Stage 3): pinning that reaches a state through OTHER states
              //      pinned pointwise at the same node.
              //
              // FIX 1 requires the pin equation to involve exactly one state.  MMPDE9's
              // initial conditions are
              //     e1  role=INITIAL  bare={ u, x }      u(0,xi) = u0( x(0,xi) )
              //     e9  role=INITIAL  bare={ x }
              // so x is pinned pointwise by e9 and [x]=0 holds at t=0; e1 on BOTH copies then
              // forces [u]=0.  FIX 1 cannot see this because e1 mentions two states, and the
              // consequence was measured: u's xi-continuity claim at t=0 was emitted, its tau
              // landed in QDEF (q's only row there, PDE being on T_INT = ALL-LB), and the
              // resulting kernel direction was
              //     c102 = +0.7071 q[t=0,xi=0]   c135 = -0.7071 q[t=0,xi=1]   c556 = -2.28e-03 tau
              // with r201: 2.827e-03*0.7071 + 0.8779*(-2.2766e-03) = 0 exactly.  q jumps at
              // t=0, funded by the multiplier.  That is Defect B.
              //
              // WHY THIS IS NOT STAGE A.  Stage A suppressed a claim because its MULTIPLIER
              // COLUMN was dependent -- a numerical rank property, invisible structurally,
              // while the claim ROW was independent.  It deleted load-bearing rows and cost
              // 24 orders of interface continuity on MBC5.  Here the claim ROW is implied by
              // PHYSICAL rows: linearising, dG/du*[u] + dG/dx*[x] = 0 from e1 on both copies,
              // with [x]=0 from e9 on both copies.  Same justification class as FIX 1.  u
              // appears bare in e1, so dG/du is structurally nonzero.
              //
              // Fixpoint closure so chains (a pinned via b, b via c) are handled; one level
              // suffices for MMPDE9 and the iteration costs nothing.
              // Seed over ALL model states, not receiver_states.  "Is this state pinned
              // pointwise at this node" is a property of the NODE, not of the equation row
              // currently being visited: receiver_states is per-row, and the state doing the
              // pinning need not be a receiver of this row.  Measured -- rev27 seeded from
              // receiver_states and fired ZERO times on MMPDE9, because u's claim is emitted
              // on the QDEF row (e4, bare={q,xg}, partials={u'xi}) where x is not a receiver,
              // so x never entered the closure and e1={u,x} showed two missing states.
              std::set<size_t> pinned_closure;
              for( auto const& [st, dep] : _mVar ){
                (void)dep;
                size_t const sid = static_cast<size_t>( st.id().second );
                if( pinned_pointwise( sid ) ) pinned_closure.insert( sid );
              }
              {
                bool grew = true;
                while( grew ){
                  grew = false;
                  for( auto const& er : _mEqn ){
                    if( er.opt->role != EqnRole::INITIAL && er.opt->role != EqnRole::BOUNDARY )
                      continue;
                    auto itm = er.dom.find( *cd.pvar );
                    if( itm == er.dom.end() || itm->second != FFDom::ALL ) continue;
                    bool here = true;
                    for( size_t d2 = 0; d2 < cdims.size() && here; ++d2 ){
                      if( d2 == d ) continue;
                      auto it2 = er.dom.find( *cdims[d2].pvar );
                      if( it2 == er.dom.end() || !in_mask( d2, it2->second ) ) here = false;
                    }
                    if( !here ) continue;
                    FFVar expr = er.var;
                    FFGraph* dg = static_cast<FFGraph*>( expr.dag() );
                    if( !dg ) continue;
                    FFSubgraph sg3 = dg->subgraph( 1, &expr );
                    std::set<size_t> bare;
                    for( auto const* op : sg3.l_op ){
                      if( !op || op->type != FFOp::VAR || !op->varout[0] ) continue;
                      if( _mVar.find( *op->varout[0] ) == _mVar.end() ) continue;
                      bare.insert( static_cast<size_t>( op->varout[0]->id().second ) );
                    }
                    if( bare.size() < 2 ) continue;          // FIX 1 already covers nst==1
                    // Exactly one state not yet in the closure -> this equation pins it.
                    size_t missing = std::numeric_limits<size_t>::max(); size_t n_missing = 0;
                    for( size_t b : bare )
                      if( !pinned_closure.count( b ) ){ missing = b; ++n_missing; }
                    if( n_missing != 1 ) continue;
                    if( pinned_closure.insert( missing ).second ) grew = true;
                  }
                }
              }
              // ---- FIX 2 (Stage A): do not ISSUE a C0 claim on a DERIVATIVE AUXILIARY
              //      whose defining LINK row is active at this node, when the PRIMITIVE it
              //      is defined from already carries a C0 claim here in the same direction.
              //
              // Measured (ocfeslv rev12, projection fill-in report): the projector onto
              // null(B) is a direct sum of RANK-1 projectors on exactly such pairs, each at
              // participation 0.5/0.5 --
              //     PDE2 : (Cl,Dz_Cl)z (Cl,Drl_Cl)rl (Sl,Dz_Sl)z (Sl,Drl_Sl)rl  = 16 = k
              //     MBC5 : (Cl,Drl_Cl) (Ml,Drl_Ml) (Pl,Drl_Pl) at rl           =  3 = k
              //     MMPDE9: (g,gg) at face 0, plus two larger blocks that are NOT pairs
              // so on both industrial targets every unit of the multiplier-block rank
              // deficiency is one of these pairs.
              //
              // HONESTY NOTE -- this rule is NOT justified the way FIX 1 is.  FIX 1 rests on
              // a genuine ROW dependency: the claim is literally the difference of two
              // ALL-mask pin rows.  Here rank(C) is FULL everywhere measured, so [a]=0 and
              // [Da]=0 are independent AS CONSTRAINTS; what is dependent is their MULTIPLIER
              // COLUMNS, which route through the same receiver rows with proportional
              // couplings.  The pair's lost direction is implied by the physical rows on
              // every CONSISTENT model (compatibility test u_k^T r, rel ~ 1e-16), but that
              // cannot be established structurally at emission time.  On an INCONSISTENT
              // model (MMPDE9) suppressing here DOES give up a real condition -- the trade
              // accepted in interface_multiplier_design_rev4.md 2.3.  Hence: reported per
              // claim, and deliberately NOT sharing FIX 1's wording, which would claim a
              // stronger justification than this rule has.
              //
              // The auxiliary is always the member suppressed: Da is DEFINED by the LINK row
              // from a, so its continuity is derived rather than primary.  Participation is
              // symmetric (0.5/0.5), so the lost direction is the same either way.
              // Survivors of FIX 1.  The pair rule must be tested against THESE, not the raw
              // list: a primitive that was itself suppressed pointwise is not claimed here,
              // so its auxiliary must not be withdrawn on the strength of it.
              std::vector<size_t> kept1; kept1.reserve( receiver_states.size() );
              auto pinned_by_link = [&]( size_t sid ) -> bool {
                for( auto const& er : _mEqn ){
                  if( er.opt->role != EqnRole::LINK ) continue;
                  auto itm = er.dom.find( *cd.pvar );
                  if( itm == er.dom.end() || itm->second != FFDom::ALL ) continue;
                  bool here = true;                 // ... and this node is inside its domain
                  for( size_t d2 = 0; d2 < cdims.size() && here; ++d2 ){
                    if( d2 == d ) continue;
                    auto it2 = er.dom.find( *cdims[d2].pvar );
                    if( it2 == er.dom.end() || !in_mask( d2, it2->second ) ) here = false;
                  }
                  if( !here ) continue;

                  FFVar expr = er.var;
                  FFGraph* dg = static_cast<FFGraph*>( expr.dag() );
                  if( !dg ) continue;
                  FFSubgraph sg2 = dg->subgraph( 1, &expr );

                  // Require exactly  d(a)/d(dir) - Da  in shape: ONE model state under a
                  // first-order FFPartial w.r.t. THIS direction, ONE bare model state, and
                  // nothing else of either kind.
                  std::set<size_t> bare_ids, part_ids;
                  bool shape_ok = true;
                  for( auto const* op : sg2.l_op ){
                    if( !op || !shape_ok ) continue;
                    if( op->sameid( typeid(FFIntegral) ) ){ shape_ok = false; continue; }
                    if( op->sameid( typeid(FFPartial) ) ){
                      auto const* pop = mc::type_cast<FFPartial const>( op );
                      if( !pop ){ shape_ok = false; continue; }
                      size_t nind = 0; bool dir_ok = false;
                      for( auto const& iv_ord : pop->Indep().expr ){
                        ++nind;
                        if( iv_ord.first.id() == cd.pvar->id() && iv_ord.second == 1 )
                          dir_ok = true;
                      }
                      if( nind != 1 || !dir_ok ){ shape_ok = false; continue; }
                      for( auto const* od : op->varin ){
                        if( !od || _mVar.find( *od ) == _mVar.end() ) continue;
                        part_ids.insert( static_cast<size_t>( od->id().second ) );
                      }
                      continue;
                    }
                    for( auto const* in : op->varin ){
                      if( !in || _mVar.find( *in ) == _mVar.end() ) continue;
                      bare_ids.insert( static_cast<size_t>( in->id().second ) );
                    }
                  }
                  if( !shape_ok ) continue;
                  if( part_ids.size() != 1 || bare_ids.size() != 1 ) continue;

                  size_t const prim = *part_ids.begin();       // a   (under the partial)
                  size_t const aux  = *bare_ids.begin();       // Da  (bare)
                  if( aux != sid || prim == sid ) continue;    // sid must BE the auxiliary

                  // ... and the primitive must itself be claimed here, in this direction.
                  // Without this the rule would be a blanket ban on auxiliary claims.
                  if( std::find( kept1.begin(), kept1.end(), prim ) == kept1.end() ) continue;
                  return true;
                }
                return false;
              };

              // Pass 1 -- FIX 1 (pointwise pins) and FIX 3 (pinned through other pinned
              // states at the same node).  Both rest on a ROW implication by physical rows.
              for( size_t sid : receiver_states ){
                bool const p1 = pinned_pointwise( sid );
                bool const p3 = !p1 && _knob_INTERFACE_SUPPRESS_PIN_CLOSURE()
                                    && pinned_closure.count( sid );
                if( !p1 && !p3 ){ kept1.push_back( sid ); continue; }
                // Record the DECISION keyed the way a CLAIM is keyed -- state, direction,
                // element, node -- with no equation-row component, because the decision has
                // none.  The leak detector below then reports any claim emitted for a triple
                // that was suppressed, which is the only way to see a visit the filter never
                // reached.  Measured need: rev28 fired twice, removed 2 receiver edges
                // (340 -> 338) and NO claim (trace_constraints 99 -> 99).
#ifdef MC__OCFESLV_PIN_LEAK_DIAG
                _pinSuppressedTriples.insert( std::make_tuple(
                  sid, static_cast<size_t>( cd.pvar->id().second ), cd.iel, midx[d] ) );
                if( options.AUDIT.RANK_SYMMETRY ){
                  // midx is the node index WITHIN an element, so it does not identify a
                  // physical point on its own: with NEL_T=3, midx=[0,4] occurs at t-node 0
                  // of every t-element.  The element multi-index disambiguates, and row_id
                  // says which visit made the decision.
                  _disp( 3 ) << "  [pin-supp] state id " << sid << " dir id "
                            << static_cast<size_t>( cd.pvar->id().second )
                            << " iel=" << cd.iel << " node=" << midx[d] << " midx=[";
                  for( size_t dd = 0; dd < midx.size(); ++dd )
                    _disp( 3 ) << ( dd ? "," : "" ) << midx[dd];
                  _disp( 3 ) << "] elem=[";
                  for( size_t dd = 0; dd < cdims.size(); ++dd )
                    _disp( 3 ) << ( dd ? "," : "" ) << cdims[dd].iel;
                  _disp( 3 ) << "]  from row_id=" << row_id << "\n";
                }
#endif
                // Per-visit this fires once per equation row at every node -- 312 lines on
                // PSA9 under CRONOS_AUDIT.  Count instead, and summarise once per plan build.
                ++( p1 ? _pinSuppressCount1[ sid ] : _pinSuppressCount3[ sid ] );
                { // rev153c: record the suppressed claim key (C0), geometry-only
                  t_InterfaceClaimKey pk{};
                  if( make_state_claim_key_from_eqn_context
                        ( _trace_claim_block( sid, static_cast<size_t>( cd.pvar->id().second ) ),
                          sid, invalid_claim_id, *cd.pvar, cd.iel, midx[d], cd.pdom->n_node, cd.pdom->n_elem,
                          cdims, midx, CLAIM_EXACT_C0, pk ) )
                    _pinSuppressedClaims.insert( pk );
                }
                if( _solve_display_level() >= 1
                 && ( p1 ? _pinSuppressCount1[ sid ] : _pinSuppressCount3[ sid ] ) == 1 )
                  _disp( 2 ) << "OCFESLV::setup ** interface: C0 claim suppressed for state id "
                            << sid << ( p1 ? " (pinned pointwise on both sides by an ALL-mask"
                                             " IC/BC at this node)"
                                           : " (pinned at this node by an ALL-mask IC/BC whose"
                                             " other states are themselves pinned here --"
                                             " closure over jointly-pinned states)" )
                            << "\n";
              }
              // Pass 2 -- FIX 2, LINK pairs, against the pass-1 survivors.
              std::vector<size_t> kept2; kept2.reserve( kept1.size() );
              for( size_t sid : kept1 ){
                if( _knob_INTERFACE_SUPPRESS_LINK_PAIR() && pinned_by_link( sid ) ){
                  if( _solve_display_level() >= 1 || options.AUDIT.COVERAGE )
                    _disp( 2 ) << "OCFESLV::setup ** interface: C0 claim WITHDRAWN for state id "
                              << sid << " (derivative auxiliary of a state that is itself "
                                 "claimed here; its multiplier cannot be independently "
                                 "supported -- see design rev4 2.3: where the pair is not "
                                 "implied this gives up a real condition)\n";
                  continue;
                }
                kept2.push_back( sid );
              }
              receiver_states.swap( kept2 );
            }
            // -------------------------------------------------------------------------
            // rev131: WHICH STATES REACH THE SITE, and how many.
            // -------------------------------------------------------------------------
            // The last unmeasured link.  Measured chain on OCFE_PDE6_blk0, IC_STRONG mono:
            //
            //     site liveness   t: MANUAL 184  AUTO 576      <-- HERE, 3.1x apart
            //     receiver edges  t: MANUAL 184  AUTO 576
            //                     x: MANUAL 584  AUTO 488      (retention differs too)
            //     candidates         MANUAL 288  AUTO 356
            //     rank test          MANUAL 276 (+12 rejected)  AUTO 356 (0 rejected)
            //
            // The mint site and the Gram-Schmidt test are correct on what they are handed;
            // the divergence originates upstream, at the traversal.  MANUAL's t-visits are
            // ~92 faces x 2 sides for ONE state, AUTO's cover BOTH -- so the question is
            // which states land in receiver_states on an evolution face, and that set is
            // built just above.
            //
            // Prints the SET, not a count, so an empty-vs-missing distinction is visible:
            //   receiver_states EMPTY  -> the site is skipped entirely, and the `continue`
            //   below is why the direction never appears in visits_by_dir at all.
            //
            // Also prints reduced_flux_state_continuity, which is the _auxDef-driven
            // predicate (_block_has_aux_link_in_direction) evaluated right here -- so if
            // aux-ness reaches the traversal, this line and the state set will disagree
            // between the two formulations in the same place.
            if( Options::_env_flag( "CRONOS_AUDIT_SITESTATES", false ) ){
              std::string dn = "dom?";
              for( auto const& kv : _mDom )
                if( kv.first.id().second == cd.pvar->id().second ){ dn = kv.first.name(); break; }
              _disp( 3 ) << "  [sitestates] dir=" << std::setw(3) << dn
                        << " blk=" << eqnopt.block_id
                        << "  reduced_flux_cont="
                        << ( _block_has_aux_link_in_direction( eqnopt.block_id, *cd.pvar )
                               ? "Y" : "n" )
                        << "  receiver_states(" << receiver_states.size() << ")={";
              for( size_t sid : receiver_states )
                _disp( 3 ) << " " << state_name_of_id( sid );
              _disp( 3 ) << " }" << ( receiver_states.empty() ? "   <- SITE SKIPPED" : "" )
                        << "\n";
            }
            if( receiver_states.empty() ) continue;

            bool const reduced_flux_state_continuity =
              _block_has_aux_link_in_direction( eqnopt.block_id, *cd.pvar );
            size_t const n_claim_kinds =
              ( dir_type == Options::IC_FLUX && !reduced_flux_state_continuity ) ? 2 : 1;
            // C1 (rev219): with the classifier on, a state that C0 says is NATURAL at second order in this
            // direction gets its C1 (flux) claim HERE, on itself, where it has a natural receiver -- instead of
            // being re-expressed as a C0 claim on a reduced auxiliary that has none (the TAU/gauge/rescue path).
            // Decided per state below; the per-site default stays as it was.
            bool const c0_reattribute = _knob_CLAIM_FROM_C0()
              && ( options.INTERFACE.IMPOSITION == Options::IC_STRONG || options.INTERFACE.IMPOSITION == Options::IC_TRACE );
            FFModel::t_ClaimMap const c0_map = c0_reattribute ? claim_classification( eqnopt.block_id ) : FFModel::t_ClaimMap();
            for( size_t state_id : receiver_states ){
              _AccT _tem( _acc_emit, _cnt_emit );

              // SUPPRESSION LEAK DETECTOR.  If this (state, direction, element, node) was
              // suppressed at some other visit, the filter did not reach this one and the
              // claim survives with a thinner receiver graph -- the worst of both outcomes,
              // since the multiplier ends up NEARLY rather than exactly unsupported.
              //
              // ===================================================================
              // rev88 FIX -- DANGLING-IF CAPTURE.  The #ifdef used to open AFTER this
              // condition, so with MC__OCFESLV_PIN_LEAK_DIAG undefined -- and it is defined
              // NOWHERE, not in this header and not in any of the 135 drivers -- the `if`
              // was left with no body of its own and SWALLOWED THE NEXT STATEMENT, which
              // is the causal evolution-direction block below.
              //
              // MEASURED, rev87: the probe placed between the two counted 0 visits on every
              // model while the claim census showed 6 claims on t for DAE5, 276 for MMPDE27,
              // 186 for MESHONLY10 and 112 for PDE0.  A code point that never fires cannot
              // have counted those.  Reproduced standalone: the guarded block executes 0
              // times, and -Wall -Wextra -Wmisleading-indentation are ALL SILENT, so no
              // build log could ever have shown it.
              //
              // And the condition could not have been true in any case: _pinSuppressedTriples
              // is only ever written inside the SAME #ifdef (rev88:25925), so count() is
              // identically 0.  The whole construct is a diagnostic; the #ifdef now encloses
              // the ENTIRE statement, which is what it should always have done.
              //
              // CONSEQUENCE, and it is the real finding of this arc's last measurement:
              // the CAUSAL EVOLUTION-DIRECTION INTERFACE has been unreachable in every build
              // in the corpus, on every model, REGARDLESS of evolution_input_may_jump.  The
              // gate was never the binding constraint, and rev72's [evoface] instrument was
              // faithfully reporting a value that had no consumer.
              // ===================================================================
#ifdef MC__OCFESLV_PIN_LEAK_DIAG
              if( options.AUDIT.RANK_SYMMETRY
               && _pinSuppressedTriples.count( std::make_tuple(
                    state_id, static_cast<size_t>( cd.pvar->id().second ),
                    cd.iel, midx[d] ) ) )
              {
                // The triple deliberately OMITS the transverse index, so it over-reports:
                // a claim is keyed by (state, direction, element, node, FACE), and the face
                // is the transverse position.  Print the full multi-index at both sites so
                // the suppressed visit and the emitting visits can be compared like with
                // like -- if they differ in a transverse slot they are different claims and
                // there is no leak at all.
                _disp( 3 ) << "  [pin-leak] claim EMITTED for a suppressed triple: state id "
                          << state_id << " dir id "
                          << static_cast<size_t>( cd.pvar->id().second )
                          << " iel=" << cd.iel << " node=" << midx[d]
                          << " midx=[";
                for( size_t dd = 0; dd < midx.size(); ++dd )
                  _disp( 3 ) << ( dd ? "," : "" ) << midx[dd];
                _disp( 3 ) << "] elem=[";
                for( size_t dd = 0; dd < cdims.size(); ++dd )
                  _disp( 3 ) << ( dd ? "," : "" ) << cdims[dd].iel;
                _disp( 3 ) << "]  from row_id=" << row_id << "\n";
              }
#endif

              // The evolution direction carries no claims (causal time discretisation, rows-level: see
              // _causal_rows_active); the former claims-based causal block (rev86-rev101) is retired.


              // Reduced-order auxiliary flux variables are ordinary collocated
              // states after RED_MAIN/RED_FULL.  In a reduced parabolic
              // direction, flux continuity is represented as C0 continuity of
              // those auxiliary states, not as a separate C1 derivative trace
              // of the primitive state.  Therefore keep auxiliary receiver
              // states and emit only C0 claims when the block has a reduced
              // auxiliary link in the face direction.
              int const claim_block = _trace_claim_block
                ( state_id, static_cast<size_t>( cd.pvar->id().second ) );

              // Geometric-C0 states: value-C0 claims only (never C1 flux).
              bool const is_geo_c0 = geo_c0_states.count( state_id ) > 0;
              size_t state_claim_kinds = is_geo_c0 ? 1 : n_claim_kinds;
              if( false && c0_reattribute && !is_geo_c0 ){   // DISABLED: no realisation for a C1 claim on a
                                                              // reduced block (frozen_weak_sat_terms fails) --
                                                              // kept for the record of what was tried
                FFVar const* sv = _state_by_id( state_id );
                if( sv ){
                  auto is = c0_map.find( *sv );
                  if( is != c0_map.end() ){
                    auto id = is->second.find( *cd.pvar );
                    if( id != is->second.end() && id->second.tag == ClaimTag::NATURAL && id->second.order == 1 )
                      state_claim_kinds = 2;   // C0 + C1 on the primitive
                  }
                }
              }

              for( size_t ikind = 0; ikind < state_claim_kinds; ++ikind ){
                int const claim_kind = ikind == 0 ? CLAIM_EXACT_C0 : CLAIM_EXACT_C1;
                t_InterfaceClaimKey key{};
                bool have_key = make_state_claim_key_from_eqn_context
                  ( claim_block, state_id, invalid_claim_id, *cd.pvar,
                    cd.iel, midx[d], cd.pdom->n_node, cd.pdom->n_elem,
                    cdims, midx, claim_kind, key );
                if( !have_key ){
                  // Fallback to the equation tensor indexing for legacy models
                  // whose receiver state/equation domains are congruent.  IC_TRACE
                  // is the constrained analogue of the weak-SAT receiver graph, so
                  // an active weak-SAT receiver slot that cannot be keyed is a
                  // malformed trace plan and must fail during setup rather than
                  // disappearing silently.
                  have_key = make_claim_key
                    ( claim_block, state_id, invalid_claim_id, *cd.pvar,
                      cd.iel, midx[d], cd.pdom->n_node, cd.pdom->n_elem,
                      d, cdims, midx, claim_kind, key );
                }
                if( !have_key ){
                  std::cerr << "OCFESLV::setup ** ERROR: interface plan (IC_TRACE): malformed "
                               "trace plan -- an active weak-SAT receiver slot cannot be keyed "
                               "(state_id=" << state_id << " dom=" << cd.pvar->name()
                            << " iel=" << cd.iel << " inode=" << midx[d]
                            << " kind=" << claim_kind << ")" << std::endl;
                  throw Exceptions( Exceptions::INTERNAL );
                }

                // Suppress tensor-product corner duplicates before adding
                // IC_TRACE claims.  Exact-row rank filtering alone can leave
                // a projected multiplier block with dependent tau columns on
                // tensor-product interfaces because the weak SAT receiver is
                // endpoint-local.  Use the same corner policy as the exact row
                // emitter/donor graph so LOOP 1 and LOOP 2 retain the same
                // independent interface basis.
                std::vector<t_InterfaceTensorCornerDim> corner_dims;
                corner_dims.reserve( d );
                for( size_t q = 0; q < d; ++q ){
                  auto const& cp = cdims[q];
                  // 2026-10-02 (row-level causal design): the evolution direction carries no interface -- an
                  // element start is an initial node -- so it is not a corner dimension (the x-terms stay)
                  if( _causal_rows_active() && cp.pvar->id() == _evolution_dom_var.id() ) continue;
                  size_t const kq = ( cp.lim == FFDom::UB )
                                  ? cp.pdom->n_node - 1 : midx[q];
                  corner_dims.push_back
                    ({ 0, cp.pdom, cp.pvar, cp.pdom->n_node, cp.pdom->n_elem, cp.iel, kq });
                }
                bool const corner_suppressed = reduced_flux_state_continuity
                 && _suppress_tensor_interface_corner
                    ( dir_type, false, eqnopt.block_id, *cd.pvar, corner_dims,
                      state_by_id( state_id ) );
#ifdef MC__OCFESLV_CORNER_SUPPRESS_PROBE
                {
                  // #earlier-dims that fire the corner predicate (the 2D-vs-3D
                  // over-fire metric): mirrors _suppress_tensor_interface_corner's
                  // IC_VALUE/IC_UPWIND/IC_FLUX arms over corner_dims (earlier dims).
                  size_t fired = 0;
                  for( auto const& o : corner_dims ){
                    bool hit = false;
                    if( dir_type == Options::IC_VALUE
                     || dir_type == Options::IC_UPWIND ){
                      if( o.iel > 0 && o.inode == 0 ) hit = true;
                    } else if( dir_type == Options::IC_FLUX ){
                      if( o.iel > 0 && o.inode == 0 ) hit = true;
                      if( o.iel + 1 < o.n_elem && o.inode + 1 == o.n_node ) hit = true;
                    }
                    if( hit ) ++fired;
                  }
                  // True duplicate multiplicity of this node across ALL dims:
                  // 2^k where k = #dims at an interior element interface.
                  size_t kdim = 0;
                  for( size_t q = 0; q < cdims.size(); ++q ){
                    auto const& cp = cdims[q];
                    size_t const iq = ( q == d )            ? midx[q]
                                    : ( cp.lim == FFDom::UB ) ? cp.pdom->n_node - 1
                                    : ( cp.lim == FFDom::LB ) ? 0
                                                              : midx[q];
                    bool const lo = ( cp.iel > 0 && iq == 0 );
                    bool const hi = ( cp.iel + 1 < cp.pdom->n_elem
                                   && iq + 1 == cp.pdom->n_node );
                    if( lo || hi ) ++kdim;
                  }
                  size_t const mult = static_cast<size_t>( 1 ) << kdim;
                  auto& rec = _corner_probe[ { state_id,
                              static_cast<size_t>( cd.pvar->id().second ) } ];
                  if( corner_suppressed ){
                    ++rec.suppressed;
                    ++rec.fired_hist[ fired ];
                    ++rec.mult_hist[ mult ];
                  } else {
                    ++rec.emitted;
                  }
                }
#endif
                if( corner_suppressed ){
                  _cornerSuppressedClaims.insert( key );   // rev153b: recorded either way (forest shadow check)
                  // rev157a (PHASE-2 STEP 1): the legacy corner dedup is RETIRED at claim
                  // emission under IC_WEAK.  There the retained claims become weak SAT
                  // penalty terms only (no exact rows, no tau, no determinacy contention);
                  // the dedup was an AUTO-only weak-discretisation asymmetry (OCFE_PDE6_blk0
                  // asymmetry B: 96 x-terms removed at t-seam nodes from the reduce_order
                  // form only).  MEASURED 2026-09-07 (blk0, same binary, 2x2 with the aux
                  // evolution route): mono MANUAL/AUTO 3.92x -> 3.21x, march untouched at
                  // 3.09x (= asymmetry C alone); with the dedup off, the aux route's weak
                  // marginal is exactly zero (3.21x with the route on or off), so STEP 2
                  // cannot degrade IC_WEAK once this is in.
                  //
                  // The rule STAYS ACTIVE for IC_TRACE/IC_STRONG.  The intended exact-mode
                  // retirement ("the forest reproduces the drops structurally") is REFUTED
                  // as a corpus invariant: it holds on blk0 (applied plans byte-equal,
                  // 272/356/36) and on PDE5/PDE6 (174/174, 8/8, rev155d), but on OCFE_PDE3
                  // live corner claims change the applied plan (exact 114->102, tau
                  // 704->766) and grow IC_TRACE n=6038 -> 6100 (+62, DETERMINED, k'=0):
                  // 62 previously-suppressed claims are realised as independent DOFs, an
                  // n-larger direction never accepted in this arc.  Exact-mode retirement
                  // is deferred behind the attribution of those 62 (which claims, and
                  // whether legacy's drops were conditions or redundancies there).
                  // Under exact modes the forest continues to consume these claims via
                  // PlanInput::corner_suppressed_edges, which is what already makes blk0's
                  // exact cells MANUAL==AUTO EQUIV.
                  //
                  // Knobs: CRONOS_PLAN_CORNER_RULE=1 restores rev156 (rule active under
                  // WEAK too).  CRONOS_PLAN_NO_CORNER_RULE=1 (rev153c instrument, kept)
                  // forces retirement everywhere INCLUDING legacy exact application; the
                  // SuperLU dgssv segfault class (legacy STRONG keep-explicit GAUGE_ONLY
                  // system, dpivotL on a singular Jacobian, root-caused 2026-09-07) is
                  // reachable only through that instrument combined with
                  // CRONOS_PLAN_TRACE_DOFS_APPLY=0.
                  // (2026-10-07: the knobs CRONOS_PLAN_CORNER_RULE and _NO_CORNER_RULE are retired -- the rule is
                  // retired under IC_WEAK and stays under the exact modes, as by default)
                  bool const retire = ( options.INTERFACE.IMPOSITION == Options::IC_WEAK );
                  if( !retire ) continue;          // exact modes: rule stays (see above)
                }

#ifdef MC__OCFESLV_CORNER_SUPPRESS_PROBE
                if( claim_kind == CLAIM_EXACT_C0 ){
                  auto& cs = _conn[ state_id ];
                  if( !cs.have_dims ){
                    cs.dims.resize( cdims.size() );
                    for( size_t q = 0; q < cdims.size(); ++q )
                      cs.dims[q] = { static_cast<size_t>( cdims[q].pvar->id().second ),
                                     cdims[q].pdom->n_node, cdims[q].pdom->n_elem };
                    cs.have_dims = true;
                  }
                  std::vector<size_t> recv( cdims.size() ), don( cdims.size() );
                  bool real_tie = false;
                  for( size_t q = 0; q < cdims.size(); ++q ){
                    size_t const nn   = cdims[q].pdom->n_node;
                    size_t const ne   = cdims[q].pdom->n_elem;
                    size_t const ielq = cdims[q].iel;
                    size_t inq;
                    if( q == d )                       inq = midx[q];
                    else if( cdims[q].lim == FFDom::UB ) inq = nn - 1;
                    else if( cdims[q].lim == FFDom::LB ) inq = 0;
                    else                                inq = midx[q];
                    recv[q] = ielq * nn + inq;
                    if( q == d ){
                      if( inq == 0 && ielq > 0 ){
                        don[q] = ( ielq - 1 ) * nn + ( nn - 1 ); real_tie = true;
                      } else if( inq == nn - 1 && ielq + 1 < ne ){
                        don[q] = ( ielq + 1 ) * nn + 0;          real_tie = true;
                      } else {
                        don[q] = recv[q];   // domain boundary, no interior donor
                      }
                    } else {
                      don[q] = recv[q];
                    }
                  }
                  if( real_tie )
                    cs.ties.emplace_back( std::move( recv ), std::move( don ) );
                }
#endif

                InterfaceReceiverEdge edge;
                edge.row_id = row_id;
                edge.claim_key       = key;
                edge.sat_type        = static_cast<OCPlan::InterfaceType>( dir_type );
                // Geometric-C0 edges always use value-C0 coupling regardless
                // of the equation's resolved interface type in this direction.
                if( is_geo_c0 ) edge.sat_type = OCPlan::IC_VALUE;
                edge.receiver_dom_id = static_cast<size_t>( cd.pvar->id().second );
                edge.receiver_iel    = cd.iel;
                edge.receiver_inode  = midx[d];
                edge.orientation     = _trace_receiver_orient( cd.iel > 0 && midx[d] == 0 );

                size_t receiver_flat = 0;
                size_t receiver_mult = 1;
                for( size_t q = 0; q < cdims.size(); ++q ){
                  auto const& cp = cdims[q];
                  size_t kq = 0;
                  if( q == d ) kq = midx[q];
                  else if( cp.lim == FFDom::UB ) kq = cp.pdom->n_node - 1;
                  else if( cp.lim == FFDom::LB ) kq = 0;
                  else kq = midx[q];
                  receiver_flat += kq * receiver_mult;
                  receiver_mult *= cp.pdom->n_node;
                  edge.receiver_block_el.push_back
                    ( std::make_pair( static_cast<size_t>( cp.pvar->id().second ), cp.iel ) );
                }
                edge.receiver_flat = receiver_flat;
                edge.receiver_emit_flat = flat;

                if( dir_type == Options::IC_UPWIND && !is_geo_c0 ){
                  t_Symbol const& sym = _symbol_for_block( eqnopt.block_id );
                  auto const& fdata = _face_data_for_block( eqnopt.block_id );
                  t_FaceData const* fd = _face_data_for_dom( fdata, cd.pvar );
                  if( !fd ) throw Exceptions( Exceptions::INTERNAL );
                  size_t k_eq = sym.vEqn.size();
                  for( size_t k = 0; k < sym.vEqn.size(); ++k )
                    if( sym.vEqn[k].id() == eqnvar.id() ){ k_eq = k; break; }
                  size_t j_st = sym.vState.size();
                  for( size_t j = 0; j < sym.vState.size(); ++j )
                    if( static_cast<size_t>( sym.vState[j].id().second ) == state_id ){ j_st = j; break; }
                  if( k_eq == sym.vEqn.size() || j_st == sym.vState.size() )
                    throw Exceptions( Exceptions::INTERNAL );
                  bool const on_left_receiver  = ( cd.iel > 0 && midx[d] == 0 );
                  bool const on_right_receiver = ( cd.iel + 1 < cd.pdom->n_elem
                                                && midx[d] + 1 == cd.pdom->n_node );
                  if( on_left_receiver )
                    edge.coupling = fd->Aplus( (arma::uword)k_eq, (arma::uword)j_st );
                  else if( on_right_receiver )
                    edge.coupling = fd->Aminus( (arma::uword)k_eq, (arma::uword)j_st );
                  else
                    edge.coupling = 0.0;
                }
                else{
                  auto itstate = _mDep.find( static_cast<long>( state_id ) );
                  if( itstate == _mDep.end() || itstate->second.first != STATE )
                    throw Exceptions( Exceptions::INTERNAL );
                  { _AccT _tsc( _acc_satc, _cnt_satc );
                  edge.coupling = _sat_coupling( edge.row_id, itstate->second.second, *cd.pvar ); }
                  double const audit_sat0 = edge.coupling;
                  bool const audit_zero = ( std::fabs( edge.coupling ) <= 64. * DBL_EPSILON );
                  bool const audit_aux  = audit_zero
                    && auxiliary_algebraic_receiver_state( eqnvar, eqnopt, *cd.pvar, state_id );
                  if( audit_aux ){
                    bool   c0_ok  = false;  int c0_why = 1;   // why refused: 1 row in no symbol table   // C2: the receiving row's coefficient of the claimed state
                    if( _classified ){   // rev291 (R4): C2's c0 read, unconditional
                      t_Symbol const& sya = _symbol_for_block( eqnopt.block_id );
                      auto read = [&]( std::vector<FFVar> const& rows, std::vector<FFVar> const& cols,
                                       std::vector<FFVar> const& coef, int kind ) -> bool {
                        size_t gr = rows.size(), jc = cols.size();
                        for( size_t g = 0; g < rows.size(); ++g ) if( rows[g].id() == eqnvar.id() ){ gr = g; break; }
                        if( !( gr < rows.size() ) ) return false;                    // not this kind of row
                        if( coef.empty() ) return true;                 // this kind, not computed
                        for( size_t j = 0; j < cols.size(); ++j ) if( static_cast<size_t>( cols[j].id().second ) == state_id ){ jc = j; break; }
                        if( !( jc < cols.size() ) ) return true;        // state not found
                        FFVar const& c = coef[ gr * cols.size() + jc ];
                        bool const is_cst = ( c.id().first == FFVar::CINT || c.id().first == FFVar::CREAL || c.cst() );
                        double const v = is_cst ? c.num().val() : 0.;
                        if( is_cst && std::fabs( v ) > 64. * DBL_EPSILON ) c0_ok = true;
                        else if( is_cst ){
                          // rev294: is this row differentiating the claimed state in the claim's direction?
                          // If so the row HAS a flux coefficient for it (of lower order than the block's
                          // principal order, which is why Ai[dir] has no entry and the lookup said zero) and the
                          // edge is a MISROUTED NATURAL RECEIVER, not a rescue case.
                          c0_why = 5;   // rev304: refined to 7 (cross-direction) below where that applies
                          if( _dag && _state_by_id( state_id ) && eqnvar.dag() ){
                            FFVar const* pst = _state_by_id( state_id );
                            FFSubgraph sgr = _dag->subgraph( 1, &eqnvar );
                            for( auto const& opr : sgr.l_op ){
                              if( !opr->sameid( typeid(FFPartial) ) ) continue;
                              auto const* popr = mc::type_cast<FFPartial const>( opr );
                              if( !popr ) continue;
                              bool dom_match = false;
                              for( auto const& kv : popr->Indep().expr )
                                if( kv.first.id() == cd.pvar->id() ){ dom_match = true; break; }

                              if( !dom_match ) continue;
                              bool hits = false;
                              for( auto const& operand : popr->Var() ){
                                if( operand.id() == pst->id() ){ hits = true; break; }
                                if( !operand.dag() ) continue;
                                FFSubgraph sgo = _dag->subgraph( 1, &operand );
                                for( auto const& opo : sgo.l_op ){
                                  for( auto const* vin : opo->varin ) if( vin && vin->id() == pst->id() ) hits = true;
                                  for( auto const* vout : opo->varout ) if( vout && vout->id() == pst->id() ) hits = true;
                                  if( hits ) break; }
                                if( hits ) break; }
                              if( hits ){ c0_why = 6; break; }   // 6 = zero-diff (misrouted natural receiver)
                            }
                            // rev295: zero-absent -- say what the row DOES contain for this state.
                            if( c0_why == 5 ){
                              bool occurs = false; std::set<long> dirs;
                              for( auto const& opr : sgr.l_op ){
                                for( auto const* vin : opr->varin ) if( vin && vin->id() == pst->id() ) occurs = true;
                                for( auto const* vout : opr->varout ) if( vout && vout->id() == pst->id() ) occurs = true;
                                if( !opr->sameid( typeid(FFPartial) ) ) continue;
                                auto const* popr = mc::type_cast<FFPartial const>( opr );
                                if( !popr ) continue;
                                bool hits2 = false;
                                for( auto const& operand : popr->Var() ){
                                  if( operand.id() == pst->id() ){ hits2 = true; break; }
                                  if( !operand.dag() ) continue;
                                  FFSubgraph sgo2 = _dag->subgraph( 1, &operand );
                                  for( auto const& opo2 : sgo2.l_op ){
                                    for( auto const* vin : opo2->varin ) if( vin && vin->id() == pst->id() ) hits2 = true;
                                    for( auto const* vout : opo2->varout ) if( vout && vout->id() == pst->id() ) hits2 = true;
                                    if( hits2 ) break; }
                                  if( hits2 ) break; }
                                if( !hits2 ) continue;
                                occurs = true;
                                for( auto const& kv : popr->Indep().expr ) dirs.insert( kv.first.id().second );
                              }
                              // rev304 (D4): a row that differentiates the claimed state in its OWN direction
                              // has no coefficient of any order in the CLAIM'S direction -- C0's IMPLIED class,
                              // undefined by construction and measured inert (M7b; MBC3's matrix).  Name it so:
                              // `cross-direction`, not `zero-absent`, which reads like a defect.
                              if( occurs && !dirs.empty() ) c0_why = 7;
                              std::ostringstream rz;
                              rz << "occurs=" << ( occurs ? "y" : "n" ) << " diff-dirs=";
                              if( dirs.empty() ) rz << "-";
                              else { bool first = true; for( long d : dirs ){ FFVar const* dv = _domain_by_id( d );
                                       rz << ( first ? "" : "," ); first = false;
                                       if( dv ) rz << *dv; else rz << "#" << d; } }
                              rz << " claim-dir=" << *cd.pvar;
                              rz << " row=" << eqnvar << " role=" << (int)eqnopt.role;
                              // rev295c: WHICH table matched, whether the claimed state is one of its columns,
                              // and a neighbouring nonzero entry of the same row -- the discriminator between
                              // "column present but zero" (the proxy consumed the bare occurrence) and
                              // "column absent" (the no-column path misreported as zero).
                              { t_Symbol const& sy = _symbol_for_block( eqnopt.block_id );
                                auto probe = [&]( std::vector<FFVar> const& rows, std::vector<FFVar> const& cols,
                                                  std::vector<FFVar> const& coef, char const* nm ){
                                  size_t gr = rows.size();
                                  for( size_t g = 0; g < rows.size(); ++g ) if( rows[g].id() == eqnvar.id() ){ gr = g; break; }
                                  if( !( gr < rows.size() ) || coef.empty() ) return;
                                  size_t jc = cols.size(), nnz = 0; std::string first_nz;
                                  for( size_t j = 0; j < cols.size(); ++j ){
                                    if( static_cast<size_t>( cols[j].id().second ) == state_id ) jc = j;
                                    FFVar const& cc = coef[ gr*cols.size() + j ];
                                    bool const cst = ( cc.id().first == FFVar::CINT || cc.id().first == FFVar::CREAL || cc.cst() );
                                    double const vv = cst ? cc.num().val() : 1.;
                                    if( !cst || std::fabs( vv ) > 64.*DBL_EPSILON ){
                                      ++nnz; if( first_nz.empty() ){ std::ostringstream o; o << cols[j]; first_nz = o.str(); } } }
                                  rz << " [" << nm << ": col=" << ( jc < cols.size() ? "present" : "ABSENT" )
                                     << " row-nnz=" << nnz << ( first_nz.empty() ? "" : " e.g." + first_nz ) << "]"; };
                                probe( sy.vAlgEqn, sy.vAlgState, sy.vAlgCoeff0, "alg" );
                                probe( sy.vEqn,    sy.vAllState, sy.vCoeff0All, "all" ); }
                              // rev295b: which row, which table, and the row's role -- the discriminator for
                              // "occurs bare yet c0 == 0".  (The stripper only zeroes FFIntegral outputs, so a
                              // model with no integral cannot lose its coefficient that way.)
                              rz << " row=" << eqnvar << " role=" << (int)eqnopt.role;
                              ++_rescueC2Zero[ state_id ][ rz.str() ];
                            }
                          }
                        }
                        else c0_why = 4;                // non-constant
                        
                        
                        return true;
                      };
                      // [0] algebraic receiver, constant c0   [5] differential receiver, constant c0
                      // rev296: C2's own table first -- the rows as written.  The two symbol tables remain the
                      // fallback so a model whose C2 table failed to build behaves exactly as before.
                      if( !read( sya.vC2Eqn, sya.vC2State, sya.vC2Coef0, 0 ) )
                        if( !read( sya.vAlgEqn, sya.vAlgState, sya.vAlgCoeff0, 0 ) )
                          read( sya.vEqn, sya.vAllState, sya.vCoeff0All, 5 );
                    }
                    // ---- rev227 (C2): is the receiving row's coefficient of the claimed state a nonzero CONSTANT? ----
                    // rev298: c2_derived went with the S1/S1d assignment it fed (rev291).
                    // 2026-10-13 (4R.25): the derived VALUE is no longer assigned to the coupling.  It was, eight
                    // lines above the unconditional fallback below, so it was overwritten on every edge and had been
                    // dead since that fallback became unconditional -- the corpus has always run on the convention.
                    // MEASURED before removing it: CRONOS_RESCUE_C2=0 (skip the whole block) gave a run identical to
                    // the default, to the iteration count.  What c0 still decides is the LABEL and the counts: an
                    // edge whose c0 is not a nonzero constant is coupling_undefined, which the plan reads.
                    if( _knob_RESCUE_C2() ){
                      if( c0_ok ){ ++_rescueC2[ state_id ][0]; }
                      else       { edge.coupling_undefined = true;            ++_rescueC2[ state_id ][1];
                                   ++_rescueC2Why[ state_id ][ c0_why ]; }   // rev293: WHY, reported with the counts
                    }
                    // rev292 (R5a): where C2 could not derive c0 (the receiving row's coefficient of the
                    // claimed state is not a nonzero constant) the edge keeps the fallback and is LABELLED
                    // coupling_undefined.  The value is a CONVENTION, not a model quantity; the corpus reads
                    // it on one cell (sweep M7b: PDE25's stiff LGL/exp WEAK cell) and nowhere else.
                    edge.coupling = OCPlan::kRescueFallbackCoupling;
                    edge.source = OCPlan::InterfaceEdgeSource::AUXILIARY_ALGEBRAIC;
                    // (the S1/S1a/S1c/S1d commentary went with the code at rev291; C2 is the derivation.)
                    // rev178 S1: the row's own coefficient replaces the forced constant.
                    // rev179 S1a: PENALTY PATH ONLY.  Under IC_WEAK the coupling is a SAT weight
                    // and c0 is measured right (model C: 1.69e+00 -> 3.84e-06).  Under
                    // IC_TRACE/IC_STRONG it is an input to the exact-realisation rank test, and
                    // MEASURED (2026-09-10j): switching +1 -> -1 there turns 20 "tau-block 2/3"
                    // clusters full-rank, i.e. enforces 20 more exact constraints, moving model
                    // C's second state 1.29e-04 -> 4.49e-04.  Whether those constraints are owed
                    // is a design question (S1b); S1a leaves the exact path untouched.
                    // rev180 S1c: ALGEBRAIC RECEIVERS ONLY.  MEASURED (2026-09-10k): the same
                    // +1-vs-(-1) configuration is catastrophic on model C, where the receiver
                    // is a two-term algebraic identity (w - Dz_u = 0) and the +1 nulls it, and
                    // CORRECT on OCFE_PDE6_blk0, where the receiver is a differential row with
                    // other content and the MANUAL reference agrees with +1 (7.04e-06 EQUIV;
                    // switching to c0 there gives 2.26e-05, 3.2x, ** DIFFER **).  So the rule
                    // is confined to the population where the failure was measured and the
                    // mechanism is understood.  Differential receivers keep the forced rescue
                    // and are counted; they are S1b's question, with the matrix as its gate.
                    // rev298: weak_here guarded the composite sign rule, deleted at rev291 (R4).
                    // rev185 S1d (CRONOS_RESCUE_FROM_ROW=2, experiment): ALSO apply c0 to SIGN-DEFAULTED
                    // edges of ANY receiver kind.  T4 (2026-09-11a): on blk0 every c0<0 differential
                    // edge is sign-DECIDABLE and the composite rule already outputs c0 (400 == 400); the
                    // S1a regression was a DOUBLE NEGATION (override before the rule), not a row-kind
                    // effect.  Decidable edges are never touched here, so that cannot recur.  MBC6
                    // carries 1,204 {diff, DEFAULTED, c0<0} edges under IC_STRONG and passes at +1;
                    // this knob measures what -1 (= c0) does to them.
                  }
                  // Geometric-C0 fallback: identity coupling for states with
                  // no principal-symbol or auxiliary-algebraic receiver in this
                  // direction.  The state was already validated by
                  // geometric_c0_receiver_state() during augmentation above.
                  if( std::fabs( edge.coupling ) <= 64. * DBL_EPSILON && is_geo_c0 ){
                    // rev292 (R5a): the geometric-C0 fallback is an IDENTITY coupling by construction -- the
                    // claimed state enters its receiver undifferentiated with coefficient 1 -- so the constant
                    // is the model's value here, not a convention.
                    edge.coupling = OCPlan::kRescueFallbackCoupling;
                    edge.source = OCPlan::InterfaceEdgeSource::GEOMETRIC_C0;
                  }
                  { auto& rc2 = _couplingLog[ { state_id, static_cast<size_t>( cd.pvar->id().second ) } ];
                    rc2.value = edge.coupling; ++rc2.n;
                    rc2.refused = ( edge.source != OCPlan::InterfaceEdgeSource::PRINCIPAL_SYMBOL ); }
                  { auto& a = _ifaceEdgeAudit[ state_id ];
                    ++a[0];
                    if( audit_zero ) ++a[1];
                    if( audit_aux  ) ++a[2];
                    if( edge.source == OCPlan::InterfaceEdgeSource::GEOMETRIC_C0 ) ++a[3];
                    if( std::fabs( edge.coupling ) <= 64. * DBL_EPSILON ) ++a[4];
                    double& mx = _ifaceEdgeMaxSat[ state_id ];
                    mx = std::max( mx, std::fabs( audit_sat0 ) );
                  }
                }

                // ---- v2 interface decision layer (inlined at end of header) ---- //
                // Route the exactness flags through the rewritten rule: every
                // non-BC continuity receiver becomes an exact tau-able value-C0
                // constraint.  coupling/source above are untouched so IC_WEAK
                // stays byte-identical.
                bool const aux_alg_receiver = auxiliary_algebraic_receiver_state
                  ( eqnvar, eqnopt, *cd.pvar, state_id );
                _apply_interface_edge_flags_v2( edge, eqnopt, aux_alg_receiver, is_geo_c0 );


                if( is_geo_c0 ){
                  // Defer geometric-C0 edges for post-loop global filtering.
                  // Pre-compute trace edge keys now while cdims/midx are in scope.
                  DeferredGeoCandidate cand;
                  cand.edge      = edge;
                  cand.claim_key = key;
                  cand.claim_key.block_id = _trace_claim_block( key.state_id, key.dom_id );
                  cand.claim_key.eqn_id   = invalid_claim_id;
                  cand.has_trace_keys = build_trace_edge_keys
                    ( cand.claim_key, *cd.pvar, cdims, midx, cand.group, cand.vlo, cand.vhi );
                  deferred_geo.push_back( std::move( cand ) );
                }
                else{
                  // Emit non-geometric edge immediately and record coverage.
                  covered_state_face.insert(
                    { key.state_id, edge.receiver_dom_id, key.iel_lo } );

                  ++n_receiver_edges;
                  if( !_lifted_prejump_receiver( edge ) ){                // a transition's claim receives AFTER the jump only
                    if( edge.hard_physical_row ) ++n_hard_row_edges;
                    if( edge.tangential_physical_boundary ) ++n_hard_tangent_edges;
                    if( edge.trace_tau_allowed ){
                      ++n_trace_edges_allowed;
                      add_raw_trace_edge( key, *cd.pvar, cdims, midx );
                    }
                    else{
                      ++n_weak_only_edges;
                    }
                    draft.receiver_edges.push_back( edge );
                  }
                }
              }
            }
          }
        }
      }
    }

    // Post-loop filtering of deferred geometric-C0 candidates.
    // A candidate is emitted only if its (state_id, dom_id, iel_lo) is NOT
    // covered by any principal-symbol or auxiliary-algebraic receiver edge
    // from ANY equation row.  This ensures the geometric fallback fires only
    // for states that are genuinely uncovered in a direction, regardless of
    // which equation rows were processed first.
    for( auto& cand : deferred_geo ){
      auto const coverage_key = std::make_tuple(
        cand.claim_key.state_id, cand.edge.receiver_dom_id, cand.claim_key.iel_lo );
      if( covered_state_face.count( coverage_key ) > 0 ){
        ++n_geometric_c0_discarded_covered;
        continue;
      }

      // Emit the surviving geometric-C0 edge.
      ++n_receiver_edges;
      ++n_geometric_c0_edges;
      if( cand.edge.hard_physical_row ) ++n_hard_row_edges;
      if( cand.edge.tangential_physical_boundary ) ++n_hard_tangent_edges;
      if( cand.edge.trace_tau_allowed ){
        ++n_trace_edges_allowed;
        // Add raw trace edge from pre-computed keys.
        if( raw_claim_seen.insert( cand.claim_key ).second ){
          if( cand.has_trace_keys ){
            if( cand.vlo < cand.vhi || cand.vhi < cand.vlo )
              raw_edges[cand.group].push_back( { cand.claim_key, cand.vlo, cand.vhi } );
          }
          else{
            // Fallback: build keys directly from claim key fields.
            t_TraceGroupKey group;
            group.state_id = cand.claim_key.state_id;
            group.kind     = cand.claim_key.kind;
            group.coord    = { { cand.claim_key.dom_id,
                                 static_cast<long long>( cand.claim_key.iel_lo ) },
                               { static_cast<size_t>(-1),
                                 static_cast<long long>( cand.claim_key.face ) } };
            t_TraceVertexKey vlo, vhi;
            vlo.state_id = vhi.state_id = cand.claim_key.state_id;
            vlo.dof = { { cand.claim_key.dom_id, cand.claim_key.iel_lo,     cand.claim_key.face } };
            vhi.dof = { { cand.claim_key.dom_id, cand.claim_key.iel_lo + 1, cand.claim_key.face } };
            raw_edges[group].push_back( { cand.claim_key, vlo, vhi } );
          }
        }
      }
      else{
        ++n_weak_only_edges;
      }
      // Geometric C0 enforces value-continuity in a degenerate direction.  Its
      // trace claim has just been registered above (so the primitive is still
      // pinned by an exact C0 row), but the stored edge must be realised as an
      // exact replacement, never a free multiplier: a tau here injects into the
      // primitive/LINK rows and lets the derivative auxiliaries (Drd_Cd, ...)
      // drift, the spurious converged solution seen when geometric edges were
      // tau-able.  Marking it exact-only routes it through the pivot-free
      // exact-replacement selection in _store_frozen_interface_plan.
      cand.edge.trace_tau_allowed = false;
      if( !_lifted_prejump_receiver( cand.edge ) ) draft.receiver_edges.push_back( std::move( cand.edge ) );
    }

    draft.diagnostics.geometric_candidates = deferred_geo.size();
    draft.diagnostics.geometric_discarded_covered = n_geometric_c0_discarded_covered;
    draft.diagnostics.geometric_emitted = n_geometric_c0_edges;
    for( auto const& kv : _rescueC2 ){
      FFVar const* sv = _state_by_id( kv.first );
      _disp( 2 ) << "OCFESLV::setup ** interface: rescue C2 (rev227) -- state ";
      if( sv ) _disp( 2 ) << *sv; else _disp( 2 ) << "#" << kv.first;
      _disp( 2 ) << ": " << kv.second[0] << " edge(s) with a nonzero constant c0 in the receiving row"
                    " (the coupling itself is the fallback convention, not this value), "
                << kv.second[1] << " REFUSED";
      // rev293: WHY.  The old text said "coefficient not constant"; PDE25's refusals are all "zero" -- the state
      // enters the receiving row only differentiated, so there is no c0 to derive or evaluate.
      if( kv.second[1] ){
        auto const& y = _rescueC2Why[ kv.first ];
        // rev304 (D4): `cross-direction` is the IMPLIED class (the row differentiates the state in its own
        // direction, not the claim's): undefined by construction, measured inert, NOT a defect.  `zero-absent`
        // is now only the genuinely coefficient-free case.  `zero-diff` has never been observed (corpus: 0).
        char const* nm[8] = { "", "no-table", "table-empty", "no-column", "non-constant", "zero-absent",
                              "zero-diff", "cross-direction" };
        _disp( 2 ) << " (why:"; for( int k = 1; k < 8; ++k ) if( y[k] ) _disp( 2 ) << " " << nm[k] << "=" << y[k];
        _disp( 2 ) << "; fallback kept, labelled)";
      }
      _disp( 2 ) << std::endl;
      // rev295/rev304: the breakdown for this state, one line per distinct row content.  A `cross-direction`
      // line is informational -- the coupling is undefined by construction; a `zero-absent` line is the case
      // worth reading, and rev296's C2 table removed the last of those (PDE25/PDE26/nest4).
      auto itz = _rescueC2Zero.find( kv.first );
      if( itz != _rescueC2Zero.end() )
        for( auto const& z : itz->second )
          _disp( 2 ) << "OCFESLV::setup ** interface: rescue C2 zero-absent -- state "
                    << ( sv ? [&]{ std::ostringstream o; o << *sv; return o.str(); }() : "#" + std::to_string( kv.first ) )
                    << ": " << z.second << " edge(s)  " << z.first << std::endl;
    }
    _rescueC2Why.clear(); _rescueC2Zero.clear();

    size_t nraw = 0;
    for( auto const& kv : raw_edges ) nraw += kv.second.size();

    auto dom_by_id = [&]( size_t did ) -> FFVar const*
    {
      for( auto const& [dv, dom] : _mDom ){
        (void)dom;
        if( static_cast<size_t>( dv.id().second ) == did ) return &dv;
      }
      return nullptr;
    };

    auto trace_claim_coeffs = [&]( t_InterfaceClaimKey const& claim,
                                   std::vector<std::pair<size_t,double>>& coeffs ) -> bool
    {
      coeffs.clear();
      FFVar const* pstate = state_by_id( claim.state_id );
      FFVar const* pvar   = dom_by_id  ( claim.dom_id );
      if( !pstate || !pvar ) return false;

      auto itdep = _mVar.find( *pstate );
      auto itdom = _mDom.find( *pvar );
      if( itdep == _mVar.end() || itdom == _mDom.end() ) return false;

      FFDom const& face_dom = itdom->second;
      if( claim.iel_lo + 1 >= face_dom.n_elem ) return false;

      std::map<FFVar const*, FFDom const*, lt_FFVar> mpDom;
      for( auto const& dv : itdep->second ) mpDom[&dv] = &_mDom.at( dv );
      if( mpDom.find( pvar ) == mpDom.end() ) return false;

      std::map<FFVar, size_t, lt_FFVar> ndx_el1, ndx_el2;
      ndx_el1[*pvar] = claim.iel_lo;
      ndx_el2[*pvar] = claim.iel_lo + 1;

      size_t rem = claim.face;
      size_t f_local = 0;
      size_t f_mult  = 1;
      for( auto const& [opvar, opdom] : mpDom ){
        if( opvar->id() == pvar->id() ) continue;
        size_t const radix = opdom->n_elem * opdom->n_node;
        size_t const enc = rem % radix;
        rem /= radix;
        size_t const ed = enc / opdom->n_node;
        size_t const kd = enc % opdom->n_node;
        if( ed >= opdom->n_elem || kd >= opdom->n_node ) return false;
        ndx_el1[*opvar] = ed;
        ndx_el2[*opvar] = ed;
        f_local += kd * f_mult;
        f_mult  *= opdom->n_node;
      }
      if( rem != 0 ) return false;

      auto const so = this->stride( mpDom, pvar, true );
      size_t const stride = so.first, outer = so.second;
      if( !stride || !outer || f_local >= stride * outer ) return false;

      size_t const inc1 = pos_state( *pstate, ndx_el1 );
      size_t const inc2 = pos_state( *pstate, ndx_el2 );
      size_t const sface = f_local % stride;
      size_t const oface = f_local / stride;

      int const ord = ( claim.kind == CLAIM_EXACT_C1 ) ? 1 : 0;
      std::vector<double> wl, wr;
      if( !face_dom.w_eval_lagrange( wl,  1.0, ord )
       || !face_dom.w_eval_lagrange( wr, -1.0, ord ) )
        return false;
      if( wl.size() != face_dom.n_node || wr.size() != face_dom.n_node ) return false;

      double scale_l = 1.0, scale_r = 1.0;
      if( ord == 1 ){
        scale_l = 2.0 / face_dom.elem_width( claim.iel_lo );
        scale_r = 2.0 / face_dom.elem_width( claim.iel_lo + 1 );
      }

      coeffs.reserve( 2 * face_dom.n_node );
      for( size_t k = 0; k < face_dom.n_node; ++k ){
        size_t const il = inc1 + sface + k * stride + oface * face_dom.n_node * stride;
        size_t const ir = inc2 + sface + k * stride + oface * face_dom.n_node * stride;
        double const cl =  wl[k] * scale_l;
        double const cr = -wr[k] * scale_r;
        if( std::fabs( cl ) > 64. * DBL_EPSILON ) coeffs.emplace_back( il, cl );
        if( std::fabs( cr ) > 64. * DBL_EPSILON ) coeffs.emplace_back( ir, cr );
      }
      return !coeffs.empty();
    };

    std::vector<t_TraceRawEdge> candidates;
    candidates.reserve( nraw );
    for( auto const& kv : raw_edges )
      for( auto const& edge : kv.second ) candidates.push_back( edge );

    auto trace_edge_priority = [&]( t_TraceRawEdge const& e ) -> std::tuple<int,int,size_t,size_t,size_t,size_t,int>
    {
      FFVar const* pdom = dom_by_id( e.claim.dom_id );
      bool const lb_value_anchor = pdom
        && _state_has_value_trace_face( e.claim.block_id, e.claim.state_id, *pdom, FFDom::LB );
      bool const ub_value_anchor = pdom
        && _state_has_value_trace_face( e.claim.block_id, e.claim.state_id, *pdom, FFDom::UB );
      int const anchor_rank = lb_value_anchor ? 0 : ( ub_value_anchor ? 1 : 2 );
      int const kind_rank = ( e.claim.kind == CLAIM_EXACT_C0 || e.claim.kind == CLAIM_EXACT_ROW ) ? 0 : 1;
      return std::make_tuple( kind_rank, anchor_rank, e.claim.state_id, e.claim.dom_id,
                              e.claim.iel_lo, e.claim.face, e.claim.kind );
    };

    std::stable_sort( candidates.begin(), candidates.end(),
      [&]( t_TraceRawEdge const& a, t_TraceRawEdge const& b )
      { return trace_edge_priority( a ) < trace_edge_priority( b ); } );

    // Retain an algebraically independent set of exact trace rows.  The older
    // graph-only spanning forest treated each one-sided trace as if it were a
    // single endpoint coefficient.  That is not true for Radau/high-side and
    // other interpolated traces, where the exact row is a full Lagrange
    // functional.  Rank filtering on the actual row coefficients removes
    // tensor-corner cycles without deleting the rows that anchor neighbouring
    // slabs/strips.
    size_t const n_state_cols = draft.trace_var_offset ? draft.trace_var_offset : draft.n_coll_var;
    // Modified Gram-Schmidt to retain an algebraically independent set of exact
    // trace rows.  Each trace functional's coefficient vector touches only the
    // few collocation columns under its interface node, while n_state_cols is the
    // full state column count (~1e4).  A dense MGS does
    // O(candidates * rank * n_state_cols) work iterating mostly-zero entries --
    // the interface-plan bottleneck on large problems.  Keep r and the qbasis
    // vectors sparse (index-sorted) so dot/axpy cost scales with actual support
    // overlap; rows with disjoint support (different states/faces) cost ~nothing.
    // The arithmetic, thresholds, and accept/reject decisions are identical to
    // the dense version -- only the zeros are skipped.
    typedef std::vector<std::pair<size_t,double>> SparseVec;
    auto sparse_dot = []( SparseVec const& a, SparseVec const& b ) -> double {
      double d = 0.0; size_t i = 0, j = 0;
      while( i < a.size() && j < b.size() ){
        if( a[i].first < b[j].first ) ++i;
        else if( a[i].first > b[j].first ) ++j;
        else { d += a[i].second * b[j].second; ++i; ++j; }
      }
      return d;
    };
    auto sparse_axpy = []( SparseVec& r, double s, SparseVec const& q ){
      SparseVec out; out.reserve( r.size() + q.size() );   // r := r - s*q
      size_t i = 0, j = 0;
      while( i < r.size() && j < q.size() ){
        if( r[i].first < q[j].first ){ out.push_back( r[i] ); ++i; }
        else if( r[i].first > q[j].first ){ out.emplace_back( q[j].first, -s*q[j].second ); ++j; }
        else { out.emplace_back( r[i].first, r[i].second - s*q[j].second ); ++i; ++j; }
      }
      for( ; i < r.size(); ++i ) out.push_back( r[i] );
      for( ; j < q.size(); ++j ) out.emplace_back( q[j].first, -s*q[j].second );
      r.swap( out );
    };
    std::vector<SparseVec> qbasis;
    std::vector<std::pair<size_t,double>> sp;
    size_t n_rank_rejected = 0;
    // -----------------------------------------------------------------------------
    // rev130: THE CLAIM MINTING SITE, logged per candidate.
    // -----------------------------------------------------------------------------
    // This loop IS where claims are decided, and the decision is a GRAM-SCHMIDT RANK
    // TEST: a candidate's continuity row r is orthogonalised against the rows already
    // accepted, and kept only if what remains has norm > 1e-10*max(1,|r|).  So the
    // census counts are not "how many claims the model wants" -- they are "how many
    // candidate rows were linearly INDEPENDENT, in the order presented".
    //
    // That reframes the OCFE_PDE6_blk0 question.  MANUAL mints 276 and AUTO 356 with
    // every other examined property identical (roles, domain specs, interface kinds,
    // evolution domain, differential/algebraic split -- all measured equal).  The gap
    // therefore has to be visible HERE, as candidates that AUTO accepts and MANUAL
    // rejects, or that MANUAL never presents at all.
    //
    // Three outcomes are distinguished, because collapsing them hides the answer:
    //     ACCEPT   independent -> minted
    //     reject   rank-dependent on an already-accepted row
    //     NO-ROW   trace_claim_coeffs failed: no continuity row could be built, so the
    //              candidate never reached the rank test.  This early-out was
    //              previously invisible and is NOT counted in n_rank_rejected.
    //
    // Diff the two logs; the first differing line is the divergence.
    static bool const kMint = Options::_env_flag( "CRONOS_AUDIT_MINT", false );
    auto mint_log = [&]( t_InterfaceClaimKey const& ck, char const* verdict,
                         double nrm0, double nrm )
    {
      std::string dn = "dom?";
      for( auto const& kv : _mDom )
        if( static_cast<size_t>( kv.first.id().second ) == ck.dom_id ){
          dn = kv.first.name(); break; }
      _disp( 3 ) << "  [mint] " << verdict
                << " state=" << std::left << std::setw(14) << state_name_of_id( ck.state_id )
                << std::right
                << " dir=" << std::setw(3) << dn
                << " iel_lo=" << std::setw(3) << ck.iel_lo
                << " face=" << std::setw(4) << ck.face
                << " kind=" << ck.kind
                << " eqn=" << std::setw(3) << ck.eqn_id
                << "  |r|=" << std::scientific << std::setprecision(3) << nrm0
                << " resid=" << nrm
                << " ratio=" << ( nrm0 > 0. ? nrm/nrm0 : 0. )
                << std::defaultfloat << "\n";
    };

    // A BANNER PER PASS, so each block of [mint] lines is self-identifying.  Without it
    // the only way to tell one formulation's pass from another's is which state names
    // appear -- and lines mentioning only the shared state T() then land in whichever
    // block was last seen, which is a guess, not a split.
    if( kMint ){
      _disp( 3 ) << "  [mint] ===== pass: " << candidates.size() << " candidate(s), states:";
      for( auto const& kv : _mVar ) _disp( 3 ) << " " << kv.first.name();
      _disp( 3 ) << "  imposition="
                << ( options.INTERFACE.IMPOSITION == Options::IC_WEAK   ? "IC_WEAK"
                   : options.INTERFACE.IMPOSITION == Options::IC_TRACE  ? "IC_TRACE"
                                                                   : "IC_STRONG" )
                << ( _stepwise_requested() ? "  march" : "  mono" ) << " =====" << std::endl;
    }

    auto _t_mgs0 = std::chrono::steady_clock::now();
    for( auto const& edge : candidates ){
      // Retain primitive evolution-domain trace claims whenever they arise
      // from the same weak-SAT receiver graph.  IC_TRACE is intended to be the
      // constrained analogue of IC_WEAK, so evolution/time-interface penalties
      // must not be removed here unless IC_WEAK would also suppress them.
      if( !trace_claim_coeffs( edge.claim, sp ) ){
        if( kMint ) mint_log( edge.claim, "NO-ROW ", 0.0, 0.0 );
        continue;
      }
      SparseVec r;
      r.reserve( sp.size() );
      for( auto const& ic : sp )
        if( ic.first < n_state_cols ) r.emplace_back( ic.first, ic.second );
      std::sort( r.begin(), r.end(),
                 []( std::pair<size_t,double> const& x,
                     std::pair<size_t,double> const& y ){ return x.first < y.first; } );
      { SparseVec c; c.reserve( r.size() );   // sum duplicate columns (dense did r[col]+=)
        for( auto const& e : r ){
          if( !c.empty() && c.back().first == e.first ) c.back().second += e.second;
          else c.push_back( e ); }
        r.swap( c ); }
      double nrm0 = 0.0;
      for( auto const& e : r ) nrm0 += e.second*e.second;
      nrm0 = std::sqrt( nrm0 );
      if( nrm0 <= 64. * DBL_EPSILON ) continue;
      for( auto const& q : qbasis ){
        double const dot = sparse_dot( r, q );
        if( std::fabs( dot ) <= 64. * DBL_EPSILON ) continue;
        sparse_axpy( r, dot, q );
      }
      double nrm = 0.0;
      for( auto const& e : r ) nrm += e.second*e.second;
      nrm = std::sqrt( nrm );
      if( nrm <= 1e-10 * std::max( 1.0, nrm0 ) ){
        ++n_rank_rejected;
        // rev160: record the REJECT branch only.  The NO-ROW and zero-norm early-outs above are
        // NOT recorded: a claim with no exact row, or with an identically zero functional, is
        // not shown to be enforced by anything -- only rank dependence on the ACCEPTED rows is.
        _rankImpliedClaims.insert( edge.claim );
        if( kMint ) mint_log( edge.claim, "reject ", nrm0, nrm );
        continue;
      }
      double const inv = 1.0 / nrm;
      for( auto& e : r ) e.second *= inv;
      qbasis.push_back( std::move( r ) );
      claims[ edge.claim ] = 1;
      if( kMint ) mint_log( edge.claim, "ACCEPT ", nrm0, nrm );
    }
    if( _solve_display_level() >= 2 )
      _disp( 3 ) << "    [iplan-mgs] trace rank-filter "
                << std::chrono::duration<double,std::milli>(
                     std::chrono::steady_clock::now() - _t_mgs0 ).count()
                << " ms (cand=" << candidates.size() << " rank=" << qbasis.size()
                << " rej=" << n_rank_rejected << " ncol=" << n_state_cols << ")\n";

    // Freeze the receiver-edge classification against the retained exact-row
    // basis.  Earlier code marked every raw weak-SAT receiver candidate as
    // trace_tau_allowed before the rank filter.  That allowed LOOP 2 to
    // rediscover raw/rank-rejected receiver slots and either inject no tau or
    // silently fall back to ordinary weak jump SAT.  In the projected TRACE
    // formulation, a receiver edge may carry a free multiplier if and only if
    // its exact continuity claim survived into the retained trace basis.
    size_t n_trace_edges_retained = 0;
    for( auto& edge : draft.receiver_edges ){
      if( !edge.trace_tau_allowed ) continue;
      t_InterfaceClaimKey slot = edge.claim_key;
      slot.block_id = _trace_claim_block( slot.state_id, slot.dom_id );
      slot.eqn_id   = invalid_claim_id;
      auto it = claims.find( slot );
      if( it == claims.end() || !it->second ){
        edge.trace_tau_allowed = false;
        continue;
      }
      edge.claim_key = slot;
      ++n_trace_edges_retained;
    }

#ifdef MC__OCFESLV_CORNER_SUPPRESS_PROBE
    {
      _disp( 2 ) << "OCFESLV::_corner_suppress_probe ** reduced-aux tensor-corner "
                   "suppression per (state,direction).\n"
                   "  emitted/suppressed = edges that survived/were dropped by the "
                   "corner dedup; fired{n:c} = #earlier-dims triggering (<=1 in 2D, "
                   ">=2 only at 3D edges/corners); mult{m:c} = duplicate node "
                   "multiplicity 2^k.\n";
      for( auto const& [k, rec] : _corner_probe ){
        FFVar const* pst  = state_by_id( k.first );
        FFVar const* pdom = dom_by_id ( k.second );
        _disp( 2 ) << "    state=" << ( pst  ? pst->name()  : std::to_string( k.first ) )
                  << " dir="     << ( pdom ? pdom->name() : std::to_string( k.second ) )
                  << " emitted=" << rec.emitted
                  << " suppressed=" << rec.suppressed;
        if( !rec.fired_hist.empty() ){
          _disp( 2 ) << " fired{";
          for( auto const& [f,c] : rec.fired_hist ) _disp( 2 ) << f << ":" << c << " ";
          _disp( 2 ) << "}";
        }
        if( !rec.mult_hist.empty() ){
          _disp( 2 ) << " mult{";
          for( auto const& [m,c] : rec.mult_hist ) _disp( 2 ) << m << ":" << c << " ";
          _disp( 2 ) << "}";
        }
        _disp( 2 ) << "\n";
      }
    }
#endif

#ifdef MC__OCFESLV_CORNER_SUPPRESS_PROBE
    {
      _disp( 2 ) << "OCFESLV::_corner_connectivity_probe ** value-C0 continuity "
                   "connectivity per state.  Each duplicate node's copies are "
                   "union-found over the SURVIVING C0 ties; components>1 == a "
                   "genuine float (under-determined node).  Per multiplicity m: "
                   "nodes / float-nodes / worst component count.\n";
      // Per-state set of floating physical nodes (gcoord keys), for the
      // cross-state identity check printed after the loop: confirm the auxes'
      // floats are the SAME physical nodes as U's before scoping any fix.
      std::map<std::string,std::set<std::vector<size_t>>> float_nodes_by_state;
      for( auto const& [sid, cs] : _conn ){
        if( !cs.have_dims ) continue;
        size_t const D = cs.dims.size();

        // Evolution (time) continuity is the causal coupling, not a C0 receiver
        // tie, so the probe never captures it and every evolution-interface node
        // looks "floating".  Model the evolution interface as always-connected:
        // copies that differ ONLY across the evolution dim are unioned below, so
        // the residual float is the genuine SPATIAL under-determination.
        int evo_dim = -1;
        if( _evolution_dom_set ){
          size_t const evo_id = static_cast<size_t>( _evolution_dom_var.id().second );
          for( size_t q = 0; q < D; ++q )
            if( cs.dims[q].dom_id == evo_id ){ evo_dim = static_cast<int>( q ); break; }
        }

        std::map<std::vector<size_t>,size_t> id_of;
        std::vector<size_t> parent;
        auto getid = [&]( std::vector<size_t> const& c ) -> size_t {
          auto it = id_of.find( c );
          if( it != id_of.end() ) return it->second;
          size_t const id = parent.size();
          id_of.emplace( c, id ); parent.push_back( id );
          return id;
        };
        auto find = [&parent]( size_t x ) -> size_t {
          while( parent[x] != x ){ parent[x] = parent[parent[x]]; x = parent[x]; }
          return x;
        };

        // Pass 1: enumerate every interface node + its copies over the gcoord grid.
        std::map<std::vector<size_t>,std::vector<size_t>> node_copies;
        std::map<std::vector<size_t>,size_t> node_mult;
        std::vector<size_t> g( D, 0 ), gmax( D );
        for( size_t q = 0; q < D; ++q )
          gmax[q] = cs.dims[q].n_elem * ( cs.dims[q].n_node - 1 );
        bool done = ( D == 0 );
        while( !done ){
          std::vector<std::vector<size_t>> reps( D );
          size_t j = 0;
          for( size_t q = 0; q < D; ++q ){
            size_t const nn = cs.dims[q].n_node, ne = cs.dims[q].n_elem, gq = g[q];
            bool const iface = ( gq > 0 && gq < ne * ( nn - 1 ) && ( gq % ( nn - 1 ) ) == 0 );
            if( iface ){
              size_t const k = gq / ( nn - 1 );
              reps[q] = { k * nn + 0, ( k - 1 ) * nn + ( nn - 1 ) };
              ++j;
            } else {
              size_t iel = gq / ( nn - 1 ), inode;
              if( iel >= ne ){ iel = ne - 1; inode = nn - 1; }
              else inode = gq - iel * ( nn - 1 );
              reps[q] = { iel * nn + inode };
            }
          }
          if( j > 0 ){
            std::vector<size_t> idx( D, 0 ), copies_here;
            // non-evo key -> a representative copy id, for evolution-collapse
            std::map<std::vector<size_t>,size_t> evo_first;
            bool cdone = false;
            while( !cdone ){
              std::vector<size_t> copy( D );
              for( size_t q = 0; q < D; ++q ) copy[q] = reps[q][ idx[q] ];
              size_t const cid = getid( copy );
              copies_here.push_back( cid );
              if( evo_dim >= 0 ){
                std::vector<size_t> key = copy;
                key[ static_cast<size_t>( evo_dim ) ] = 0;   // canonicalise evo dim
                auto it = evo_first.find( key );
                if( it == evo_first.end() ) evo_first.emplace( key, cid );
                else { size_t const a = find( cid ), b = find( it->second );
                       if( a != b ) parent[a] = b; }
              }
              size_t q = 0;
              for( ; q < D; ++q ){ if( ++idx[q] < reps[q].size() ) break; idx[q] = 0; }
              if( q == D ) cdone = true;
            }
            node_copies[g] = copies_here;
            node_mult[g]   = ( static_cast<size_t>( 1 ) << j );
          }
          size_t q = 0;
          for( ; q < D; ++q ){ if( ++g[q] <= gmax[q] ) break; g[q] = 0; }
          if( q == D ) done = true;
        }

        // Pass 2: apply the surviving ties.
        for( auto const& t : cs.ties ){
          size_t const a = find( getid( t.first ) ), b = find( getid( t.second ) );
          if( a != b ) parent[a] = b;
        }

        // Pass 3: components per node, bucketed by multiplicity.
        struct Bucket { size_t nodes = 0, floats = 0, worst = 0; };
        std::map<size_t,Bucket> bym;
        std::set<std::vector<size_t>> float_nodes;
        size_t float_at_evo_lb = 0, float_interior = 0;   // localise the arm
        for( auto const& [gnode, copies] : node_copies ){
          std::set<size_t> roots;
          for( size_t id : copies ) roots.insert( find( id ) );
          size_t const comp = roots.size();
          Bucket& bk = bym[ node_mult.at( gnode ) ];
          ++bk.nodes;
          if( comp > 1 ){
            ++bk.floats; float_nodes.insert( gnode );
            // gcoord 0 in the evolution dim == the t=LB initial face.  Splitting
            // here tells initial-face-arm floats from interior ordered-dedup ones.
            if( evo_dim >= 0 && gnode[ static_cast<size_t>( evo_dim ) ] == 0 )
              ++float_at_evo_lb;
            else
              ++float_interior;
          }
          if( comp > bk.worst ) bk.worst = comp;
        }

        FFVar const* pst = state_by_id( sid );
        std::string const sname = pst ? pst->name() : std::to_string( sid );
        float_nodes_by_state[ sname ] = float_nodes;
        _disp( 3 ) << "    state=" << sname
                  << " C0_ties=" << cs.ties.size();
        for( auto const& [m, bk] : bym )
          _disp( 3 ) << "  mult" << m << "{nodes=" << bk.nodes
                    << " float=" << bk.floats << " worstComp=" << bk.worst << "}";
        _disp( 3 ) << "  [float@t=LB=" << float_at_evo_lb
                  << " interior=" << float_interior << "]\n";
      }

      // Cross-state identity check: is each state's floating-node set the SAME
      // physical nodes?  If U's floats == aux's floats, the corner suppression
      // is geometric/state-independent and the scoped reduced-aux fix targets
      // exactly U's float nodes (which U survives via an external pin).
      {
        auto setdiff_size = []( std::set<std::vector<size_t>> const& a,
                                std::set<std::vector<size_t>> const& b ){
          size_t n = 0;
          for( auto const& e : a ) if( !b.count( e ) ) ++n;
          return n;
        };
        std::string ref;
        for( auto const& [nm, st] : float_nodes_by_state )
          if( ref.empty() && !st.empty() ){ ref = nm; }
        if( !ref.empty() ){
          auto const& R = float_nodes_by_state[ ref ];
          _disp( 3 ) << "    float-node identity (reference=" << ref
                    << ", |float|=" << R.size() << "):\n";
          for( auto const& [nm, st] : float_nodes_by_state ){
            if( nm == ref ) continue;
            size_t const only_ref = setdiff_size( R, st );
            size_t const only_st  = setdiff_size( st, R );
            _disp( 3 ) << "      vs " << nm << ": |float|=" << st.size()
                      << "  identical=" << ( (only_ref==0 && only_st==0) ? "YES" : "NO" )
                      << "  (only_in_" << ref << "=" << only_ref
                      << ", only_in_" << nm << "=" << only_st << ")\n";
          }
        }
      }
    }
#endif

#ifdef MC__OCFESLV_INTERFACE_ASSIGN_DEBUG
    if( _solve_display_level() >= 2 )
      _disp( 2 ) << "OCFESLV::setup ** IC_TRACE receiver edges=" << n_receiver_edges
                << " trace-replaceable=" << n_trace_edges_allowed
                << " weak-only=" << n_weak_only_edges
                << " hard-row=" << n_hard_row_edges
                << " hard-tangent=" << n_hard_tangent_edges
                << " geometric-C0=" << n_geometric_c0_edges
                << " raw trace edges=" << nraw
                << " duplicate groups=" << raw_edges.size()
                << " rank-rejected=" << n_rank_rejected
                << " retained=" << claims.size()
                << " retained receiver edges=" << n_trace_edges_retained
                << std::endl;
    // Per-state/direction breakdown for geometric-C0 edges.
    if( _solve_display_level() >= 2 && n_geometric_c0_edges > 0 ){
      std::map<std::pair<size_t,size_t>, size_t> geo_detail;
      for( auto const& edge : draft.receiver_edges )
        if( edge.source == OCPlan::InterfaceEdgeSource::GEOMETRIC_C0 )
          ++geo_detail[{ edge.claim_key.state_id, edge.receiver_dom_id }];
      for( auto const& [key, count] : geo_detail ){
        FFVar const* pst = state_by_id( key.first );
        FFVar const* pdom = dom_by_id( key.second );
        _disp( 2 ) << "OCFESLV::setup **   geometric-C0: state="
                  << ( pst  ? pst->name()  : "?" )
                  << " direction=" << ( pdom ? pdom->name() : "?" )
                  << " edges=" << count << std::endl;
      }
    }
#endif

    // -----------------------------------------------------------------------
    // rev87 PROBE: per-DIRECTION claim census + would-drop summary.  READ-ONLY.
    // -----------------------------------------------------------------------
    // Handoff rev25 section 3.3 reports "6 taus, k=3 = nel_t-1" on OCFE_DAE5 without
    // saying which DIRECTION those claims sit on -- on a one-domain model that is
    // inferable, on MMPDE27 it is not.  Read the artefact, not the total: the same
    // discipline that showed dom0x48 to be a kernel-participation tally rather than a
    // claim count.  Emitted from inside the builder, where `claims` is the master map on
    // every path that builds one, rather than at one of the three return points.
    {
      static bool const kEvoFace = Options::_env_flag( "CRONOS_AUDIT_EVOFACE", false );
      if( kEvoFace ){
        // Self-contained lookups: no dependence on the scope of the debug-only helpers above.
        auto nm_state = [&]( size_t sid ) -> std::string {
          for( auto const& kv : _mVar )
            if( static_cast<size_t>( kv.first.id().second ) == sid ) return kv.first.name();
          return "state#" + std::to_string( sid );
        };
        auto nm_dom = [&]( size_t did ) -> std::string {
          for( auto const& kv : _mDom )
            if( static_cast<size_t>( kv.first.id().second ) == did ) return kv.first.name();
          return "dom#" + std::to_string( did );
        };
        size_t const evo_id = ( _evolution_dom_set && _evolution_dom_var.dag() )
                            ? static_cast<size_t>( _evolution_dom_var.id().second )
                            : std::numeric_limits<size_t>::max();

        std::map<size_t,size_t> by_dom;
        std::map<std::pair<size_t,size_t>,size_t> by_dom_state;
        // rev128: also break down by the resolved interface KIND.  The census reported
        // counts per (direction, state) and never which KIND of continuity condition
        // those claims are -- IC_VALUE, IC_FLUX or IC_UPWIND.  On
        // OCFE_PDE6_blk0 the two formulations differ by 100 evolution-face claims with
        // every other examined property identical, and "is one of them minting a
        // non-VALUE interface?" was not answerable from any existing output.
        std::map<std::tuple<size_t,size_t,int>,size_t> by_dom_state_kind;
        for( auto const& kv : claims ){
          ++by_dom[ kv.first.dom_id ];
          ++by_dom_state[ { kv.first.dom_id, kv.first.state_id } ];
          ++by_dom_state_kind[ { kv.first.dom_id, kv.first.state_id, kv.first.kind } ];
        }
        auto kind_name = []( int k ) -> char const* {
          switch( k ){
            case 0:  return "IC_VALUE";
            case 1:  return "IC_FLUX";
            case 2:  return "IC_UPWIND";
            case 3:  return "IC_AUTO(unresolved!)";
            default: return "IC_?";
          }
        };
        _disp( 3 ) << "  [evoface] claim census: total=" << claims.size() << "\n";
        for( auto const& kv : by_dom )
          _disp( 3 ) << "  [evoface]   direction=" << nm_dom( kv.first )
                    << ( kv.first == evo_id ? "  [EVOLUTION]" : "" )
                    << "  claims=" << kv.second << "\n";
        for( auto const& kv : by_dom_state ){
          _disp( 3 ) << "  [evoface]     " << nm_dom( kv.first.first )
                    << " / " << nm_state( kv.first.second )
                    << "  claims=" << kv.second;
          // rev128: the kind split, inline.  A state whose claims are all one kind
          // prints one entry; a mixture prints each, which is itself worth seeing.
          for( auto const& kk : by_dom_state_kind )
            if( std::get<0>( kk.first ) == kv.first.first
             && std::get<1>( kk.first ) == kv.first.second )
              _disp( 3 ) << "  [" << kind_name( std::get<2>( kk.first ) )
                        << "=" << kk.second << "]";
          _disp( 3 ) << "\n";
        }

        // -------------------------------------------------------------------------
        // rev125: THE MARCHING REFERENCE SET, stated rather than left to be inferred.
        // -------------------------------------------------------------------------
        // Marching defines the correct evolution-face constraint set, and the code says
        // so: _marchICRows maps each TRANSFER STATE to its INITIAL-closure rows and the
        // state's LB-node column, and the per-window override is
        //     residual = var[lb_col] - terminal,  Jacobian row = unit on lb_col
        // with the transfer states taken from _states_with_evolution_derivative -- "the
        // differential-state (dyn) identification".  So each window pins every
        // DIFFERENTIAL state at its lower bound and transfers NO algebraic state; the
        // algebraic ones are recovered by their own rows at every node, LB included.
        // That is the consistent-reinitialisation problem.
        //
        // The monolithic plan should therefore mint, at each interior evolution-face:
        //     every DIFFERENTIAL state, at every transverse node      <- reference
        //     no ALGEBRAIC state                                      <- recovered
        //
        // MEASURED, OCFE_PDE6_blk0 MANUAL, NEL_T=NEL_X=5 (4 faces x 25 x-nodes = 100):
        //     t / T = 92   -- EIGHT SHORT of the reference
        //     t / q =  0   -- correct: q is algebraic, not a transfer state
        // and that driver's monolithic solve is PRIMAL UNDETERMINED while its marching
        // solve is not.  The under-claim is on the DIFFERENTIAL state, not the algebraic
        // one -- which is the opposite of where this arc spent its time looking.
        //
        // This block does NOT assert the identity; it PRINTS both sides so the gap is
        // visible without a hand comparison.  Turning it into an invariant needs the
        // per-face transverse-node count, which is not available here.
        if( evo_id != (size_t)-1 ){
          std::set<size_t> diff_sta, alg_sta;
          for( auto const& kv : by_dom_state )
            if( kv.first.first == evo_id ) diff_sta.insert( kv.first.second );
          for( auto const& kv : by_dom_state )
            if( kv.first.first != evo_id && !diff_sta.count( kv.first.second ) )
              alg_sta.insert( kv.first.second );

          _disp( 3 ) << "  [evoface] MARCHING REFERENCE (evolution dir = " << nm_dom( evo_id )
                    << "):  states claimed on the evolution face =";
          if( diff_sta.empty() ) _disp( 3 ) << " <none>";
          for( size_t sid : diff_sta ) _disp( 3 ) << " " << nm_state( sid );
          _disp( 3 ) << "  |  claimed ONLY transversally (recovered, not transferred) =";
          if( alg_sta.empty() ) _disp( 3 ) << " <none>";
          for( size_t sid : alg_sta ) _disp( 3 ) << " " << nm_state( sid );
          _disp( 3 ) << "\n            marching would pin each of the first group at EVERY"
                       " transverse node of EVERY interior face; compare that against the"
                       " per-state counts above -- a shortfall is an UNDER-CLAIMED"
                       " differential state, and OCFE_PDE6_blk0 shows one (92 of 100).\n";
        }

        if( _evoProbe.ran ){
          // THE LIVENESS LINE.  Read this BEFORE any would-drop number.  It is ungated:
          // if the evolution direction is absent here, the causal site never sees it and
          // every gate question below is moot.
          _disp( 3 ) << "  [evoface] site liveness (UNGATED) -- visits by face direction:";
          if( _evoProbe.visits_by_dir.empty() ) _disp( 3 ) << "  <SITE NEVER REACHED>";
          for( auto const& kv : _evoProbe.visits_by_dir )
            _disp( 3 ) << "  " << nm_dom( kv.first )
                      << ( kv.first == evo_id ? "[EVOLUTION]" : "" ) << "=" << kv.second;
          _disp( 3 ) << "\n";
          // rev89: RECEIVER EDGES per direction.  THE metric for the differential branch.
          // MEASURED rev88 (DAE5 at L2): acted=24, drop_upstream=6, and the claim census
          // did NOT move -- still 6.  Both receivers of a seam canonicalise to ONE physical
          // trace slot ( add_raw_trace_edge nulls eqn_id and dedupes via raw_claim_seen ),
          // so dropping the upstream receiver removes an EDGE, not a CLAIM.  Counting claims
          // is therefore blind to the entire differential-branch effect; counting edges is
          // not.  The aux-flux branch does move claims, which is why rev88 saw MMPDE27 go
          // 606 -> 519 while DAE5 stayed at 6.
          {
            std::map<size_t,size_t> edges_by_dom;
            for( auto const& e : draft.receiver_edges ) ++edges_by_dom[ e.receiver_dom_id ];
            size_t etot = 0;
            for( auto const& kv : edges_by_dom ) etot += kv.second;
            _disp( 3 ) << "  [evoface] receiver edges: total=" << etot;
            for( auto const& kv : edges_by_dom )
              _disp( 3 ) << "  " << nm_dom( kv.first )
                        << ( kv.first == evo_id ? "[EVOLUTION]" : "" ) << "=" << kv.second;
            _disp( 3 ) << "\n";
          }
        }
      }
    }
  };

  if( options.INTERFACE.IMPOSITION == Options::IC_TRACE
   || options.INTERFACE.IMPOSITION == Options::IC_STRONG ){
    if( false && store_plan && _load_frozen_interface_plan( interface_claims_master ) ){
      // loaded frozen setup-time receiver/trace/tau plan
    }
    else{
      if( !_interfacePlanBuildAllowed
       || draft.state != OCPlan::State::BUILDING ){
        std::cerr << "OCFESLV::setup ** IC_TRACE/IC_STRONG interface receiver graph "
                     "is not being built by setup(); evaluation cannot rebuild "
                     "the frozen trace/tau plan." << std::endl;
        throw Exceptions( Exceptions::INTERNAL );
      }
      build_trace_claims_from_weak_sat_graph( interface_claims_master );
      if( store_plan ){
        _apply_c0_claim_filter( interface_claims_master );                             // rev218 (C1)
        PlanInput const plan_in = _pack_plan_input( interface_claims_master, draft );   // rev152b
        _store_frozen_interface_plan( interface_claims_master, draft, plan_in );
      }
    }

#ifdef MC__OCFESLV_INTERFACE_ASSIGN_DEBUG
    if( _solve_display_level() >= 2 )
      _disp( 2 ) << "OCFESLV::setup ** IC_TRACE/IC_STRONG claims from weak SAT graph: "
                << interface_claims_master.size()
                << " donor rows=0 min-cost assignment=skipped"
                << std::endl;
#endif
    claims_out = interface_claims_master;
    return;
  }

  if( options.INTERFACE.IMPOSITION == Options::IC_WEAK ){
    if( draft.frozen && draft.state == OCPlan::State::FROZEN ){
      claims_out = interface_claims_master;
      return;
    }
    if( !_interfacePlanBuildAllowed
     || draft.state != OCPlan::State::BUILDING ){
      std::cerr << "OCFESLV::setup ** IC_WEAK weak-SAT receiver graph "
                   "is not being built by setup(); evaluation cannot rebuild "
                   "the frozen receiver plan." << std::endl;
      throw Exceptions( Exceptions::INTERNAL );
    }
    build_trace_claims_from_weak_sat_graph( interface_claims_master );
    claims_out = interface_claims_master;
    return;
  }

  // No legacy IC_STRONG donor-row path remains here.  All supported
  // imposition modes have returned above: IC_WEAK keeps the weak receiver
  // graph, IC_TRACE uses the retained trace/tau graph, and IC_STRONG uses
  // the same trace/tau graph followed by setup-time Schur projection.
  claims_out = interface_claims_master;
}

template <typename U>
inline void

//==========================================================================================================
// EVALUATION AND DERIVATIVES -- collocated residual and Jacobian evaluation over the working DAG.
// (banner placed from the definitions themselves: this region is lines 24280-25176 of the pre-rev268 file)
//==========================================================================================================

OCFEBASE::_eval_eqn
( U* eqn, U const* var, U const* inp, 
  OCVar<U>& OCres, std::vector<OCVar<U>>& OCvar,  std::vector<OCVar<U>>& OCinp,
  std::vector<OCVar<U>>& OCdom, std::vector<OCVar<U>>& OCcst, std::vector<OCVar<U>>& OCwrk )
const
{
  // Validate external coefficient arrays before any residual assembly.
  if( _nCollVar && !var )
    throw Exceptions( Exceptions::INTERNAL );
  if( _nCollInp && !inp )
    throw Exceptions( Exceptions::INTERNAL );
  if( _nCollEqn && !eqn )
    throw Exceptions( Exceptions::INTERNAL );

  // ---------------------------------------------------------------------
  // Resolve effective interface condition type and imposition strategy.
  // Both affect all loops, so they are resolved once here.
  // weak must be resolved first since IC_AUTO dispatch depends on it.
  // ---------------------------------------------------------------------
  bool const weak   = ( options.INTERFACE.IMPOSITION == Options::IC_WEAK );
  bool const trace  = ( options.INTERFACE.IMPOSITION == Options::IC_TRACE );
  bool const strong = ( options.INTERFACE.IMPOSITION == Options::IC_STRONG );
  bool const weak_mult = weak && _plan.weak_multipliers();   // IC_WEAK with multipliers (C2 experiment)

  if( !_plan.frozen() || _plan.state() != OCPlan::State::FROZEN
   || _interfacePlanBuildAllowed ){
    _disp( 2 ) << "OCFESLV::_eval_eqn ** frozen interface plan is not sealed; "
                 "residual evaluation is not allowed to build or modify the "
                 "interface plan." << std::endl;
    throw Exceptions( Exceptions::INTERNAL );
  }

  _assert_frozen_interface_options_current();

  // The per-equation principal-symbol SAT coupling table is built during
  // setup.  Evaluation consumes it only; stale/missing coupling data would
  // indicate that setup was bypassed or classification was changed after
  // setup.
  if( ( weak || trace || strong )
   && ( !_satCouplingBuilt
     || _satCouplingClassified != _classified
     || _satCouplingSerial != _classificationSerial ) ){
    _disp( 2 ) << "OCFESLV::_eval_eqn ** SAT coupling cache is not current; "
                 "rerun setup() before evaluating interface SAT/TRACE terms."
              << std::endl;
    throw Exceptions( Exceptions::INTERNAL );
  }

  // ---------------------------------------------------------------------
  // Row-count self-check (Part 2): capture the residual cursor and a drop
  // counter so the tail of this function can verify LOOP 1 + LOOP 2 emit
  // exactly _nCollEqn rows, with LOOP 1 emissions == LOOP 2 drops.  The
  // identities are index-based and hold for every U; the runtime cost is
  // one pointer capture, one increment per dropped node, and one compare
  // at the end.
  // ---------------------------------------------------------------------
  U* const eqn0          = eqn;
  size_t   n_audit_emit  = 0;   // rows produced by LOOP 1 (interface conditions)
  size_t   n_audit_drop  = 0;   // IC_TRACE exact-replacement physical-row drops
  size_t   n_strong_pivot_drop = 0; // IC_STRONG Schur-projection pivot-row drops
  size_t   n_strong_exact_drop = 0; // IC_STRONG exact-replacement row drops (rank-rejected claims)

  // Row-provenance build: on the first eval after setup (the frozen plan makes emit order
  // eval-invariant), record for each emitted output row its source equation index into _mEqn
  // (or -1 for plan-owned interface/trace rows).  emit_eqn is the single sink for every row, so
  // tagging it captures the true final row order across IC_WEAK/TRACE/STRONG.
  bool const build_prov = !_rowProvenanceReady;
  if( build_prov ){ _rowEqnNdx.clear(); _rowNode.clear(); _rowFold.clear(); }

  auto emit_eqn = [&]( U const& val, int src_ndx = -1 )
  {
    if( eqn0 && static_cast<size_t>( eqn - eqn0 ) >= _nCollEqn ){
      std::cerr << "OCFESLV::_eval_eqn ** attempted to emit more residual rows "
                << "than _nCollEqn=" << _nCollEqn
                << "; interface row ownership is inconsistent."
                << " cursor=" << static_cast<size_t>( eqn - eqn0 )
                << " interface_emitted=" << n_audit_emit
                << " interface_dropped=" << n_audit_drop
                << std::endl;
      throw Exceptions( Exceptions::INTERNAL );
    }
    if( build_prov ){ _rowEqnNdx.push_back( src_ndx ); _rowNode.emplace_back(); _rowFold.emplace_back(); }
    *eqn++ = val;
  };

  if( trace || strong || weak_mult ){
    if( !_plan.trace_constraints_ready() ){
      _disp( 2 ) << "OCFESLV::_eval_eqn ** frozen exact trace constraints are "
                   "not ready; setup() must build executable trace rows."
                << std::endl;
      throw Exceptions( Exceptions::INTERNAL );
    }
    for( auto const& term : _plan.trace_constraints() ){
      if( term.row != n_audit_emit ){
        std::cerr << "OCFESLV::_eval_eqn ** frozen exact trace row order is "
                     "inconsistent: expected=" << n_audit_emit
                  << " stored=" << term.row << std::endl;
        throw Exceptions( Exceptions::INTERNAL );
      }
      U val( 0. );
      for( auto const& ic : term.coeff ){
        if( ic.first >= _nCollVar ){
          _disp( 2 ) << "OCFESLV::_eval_eqn ** frozen exact trace coefficient "
                       "column is out of range: col=" << ic.first
                    << " nvar=" << _nCollVar << std::endl;
          throw Exceptions( Exceptions::INTERNAL );
        }
        val += U( ic.second ) * var[ ic.first ];
      }
      emit_eqn( val );
      ++n_audit_emit;
    }
  }


  // Project the current ndx_el map to a BlockKey for two-level index lookup.
  // Both plan-build and eval iterate domains in lt_FFVar order, producing
  // identical (domain_var_id, element_index) sequences for the same block.
  // This projection is computed once per block (O(ndim)), not per term.
  auto project_block_key = [&]
  ( std::map<FFVar, size_t, lt_FFVar> const& ndx_block ) -> OCPlan::BlockKey
  {
    OCPlan::BlockKey bk;
    bk.reserve( ndx_block.size() );
    for( auto const& [var, iel] : ndx_block )
      bk.emplace_back( static_cast<size_t>( var.id().second ), iel );
    return bk;
  };

  // IC_STRONG: tagged flat row buffer.
  // pivot_idx: >= 0 = Schur pivot, -1 = ordinary, -2 = exact-replacement (dropped).
  struct StrongRow {
    U   value;
    int pivot_idx;   //!< >= 0 pivot, -1 ordinary, -2 exact-replacement
    std::vector<std::pair<size_t,double>> const* correction;  //!< nullptr or &frozen projection
    int ndx;         //!< source equation index into _mEqn (row-provenance tag), -1 if none
    t_RowNode node;  //!< the node half of the provenance, recorded at buffering time (the emit comes after the loop)
  };
  std::vector<StrongRow> strong_rows;
  auto make_strong_row_key = [&]( size_t eqn_id,
                                  std::map<FFVar,size_t,lt_FFVar> const& ndx_block,
                                  size_t flat )
  {
    FrozenPhysicalResidualRowKey k;
    k.row_id = eqn_id;
    k.flat   = flat;
    for( auto const& kv : ndx_block )
      k.block_el.push_back( std::make_pair( static_cast<size_t>( kv.first.id().second ), kv.second ) );
    return k;
  };

  // =========================================================================
  // LOOP 2: EQUATION RESIDUALS
  //
  // IC_STRONG: IC-aware node exclusion removes nodes claimed by loop 1.
  // IC_WEAK:   Ordinary first-order interfaces and reduced auxiliary-link
  //         directions keep all equation nodes and receive SAT penalties when
  //         enabled by equation metadata.
  //
  // SAT penalty for IC_WEAK is applied two-sidedly: at every internal
  // interface in direction d, the same penalty is added at one face of
  // EACH adjoining element, using an interface-symmetric tau
  // (_sat_tau_face) and an element-independent coupling.  The two
  // contributions are exact negatives and therefore telescope, so the
  // scheme is discretely conservative on uniform AND nonuniform grids.
  //
  //   For each state u_j active in this equation and dependent on d, and
  //   each face of the current element that lies on an internal interface:
  //     C0 jump  = u_curr(t=±1) - u_nbr(t=∓1)
  //
  //   Without classification: \delta F += \sigma_0 * \tau_face * C0_jump   (identity coupling)
  //   With classification:    \delta F += \sigma_0 * \tau_face * A_d[k,j] * C0_jump
  //
  //   For IC_FLUX, a C1 (flux) penalty is additionally added using each
  //   element's own physical-derivative rescale 2/elem_width(iel):
  //     \delta F += \sigma_1 * \tau_face / iface_width * C1_jump
  //
  //   For IC_UPWIND, the scalar state-by-state penalty is bypassed: the full
  //   classified block jump vector is formed first and multiplied by A+ on
  //   the left face of the right element and A- on the right face of the
  //   left element.
  // =========================================================================

  size_t ndx = 0;
  for( auto const& eqnrec : _mEqn ){
    auto const& eqnvar = eqnrec.var;
    auto const& eqndom = eqnrec.dom;
#ifdef MC__OCFESLV_EVAL_DEBUG
    _disp( 3 ) << "\nEQUATION " << ndx << ":\n";
#endif

    std::map<FFVar, size_t, lt_FFVar> ndx_el;
    // Iterate element-blocks of this equation's domain.  A domain-less (scalar / lumped algebraic)
    // equation has exactly one global block: emit it once with an empty ndx_el, which the body reads
    // as "global" (pos_state(v,{}), OCres.dom() empty -> ntot=1).  The non-empty branch is the exact
    // former loop (behaviour-preserving for all distributed equations; mirrors the _eval_fct idiom).
    bool scalar_done = false;
    while( true ){
      if( eqndom.empty() ){ if( scalar_done ) break; scalar_done = true; }
      else{ if( !_next_block( ndx_el, eqndom ) || ndx_el.empty() ) break; }
#ifdef MC__OCFESLV_EVAL_DEBUG
      _disp( 3 ) << "BLOCK";
      for( auto const& [v, b] : eqndom )
        _disp( 3 ) << " " << v << "=" << ndx_el[v];
      _disp( 3 ) << std::endl;
#endif

      // Update participating variables using the compact setup-time plan.
      auto const& plan_eqn = ( ndx < _depPlanEqn.size() ? _depPlanEqn[ndx]
                                                        : std::vector<t_DepUpdate>() );

      // section 4: FREE (integrated) directions = touched distributed directions NOT spent as a
      // fixed point.  A domain-less equation reduces each touched direction either by INTEGRATION
      // (free -> _next_block + accumulate, section 4a) or by POINT evaluation (fixed coordinate ->
      // eval_lagrange, section 4b).  Distributed equations leave free_dom empty (they use the
      // per-block path below).
      std::set<FFVar,lt_FFVar> free_dom;
      if( ndx < _eqnEvalDom.size() ){
        for( auto const& dv : _eqnEvalDom[ndx] ){
          if( eqndom.count( dv ) ) continue;                                            // this equation's own distributed direction (handled per-block by the LOOP-2)
          if( ndx < _eqnEvalPoint.size() && _eqnEvalPoint[ndx].count( dv ) ) continue; // fixed by OpEval node
          free_dom.insert( dv );
        }
      }

      // Zero buffer for the "g" (non-integrated) evaluation; populated lazily (see below).
      std::vector<U> s4_zeros;

      // Assemble the residual tensor for one element selection.  Factored so a domain-less
      // (scalar/lumped) equation that reduces distributed fields via FFIntegral can accumulate
      // over the integrated ("free") elements -- the equation analogue of _eval_fct's free loop.
      // When zero_integrated is set, every distributed state/input whose domain intersects free_dom
      // (i.e. sits under an FFIntegral) is loaded as ZEROS, so all integral terms vanish and OCres
      // holds the non-integrated part g(scalars) of the residual.
      auto assemble = [&]( std::map<FFVar,size_t,lt_FFVar> const& nel, bool zero_integrated )
      {
        for( auto const& upd : plan_eqn ){
          size_t inc = 0;
          switch( upd.typ ){
            case DOMAIN:  OCdom[upd.pos].update( _node_block( _vDom[upd.pos], nel ), nel ); break;
            case STATE: {
              inc = pos_state( _vVar[upd.pos], nel );
              bool zeroit = false;
              if( zero_integrated ){
                auto itd = _mVar.find( _vVar[upd.pos] );
                if( itd != _mVar.end() )
                  for( auto const& dv : itd->second ) if( free_dom.count( dv ) ){ zeroit = true; break; }
              }
              OCvar[upd.pos].update( zeroit ? s4_zeros.data() : &var[inc], nel );
              break;
            }
            case INPUT: {
              auto nel_in = nel;   // marching: remap evolution element of a distributed input to the resident window
              auto const ew = _evolution_window();
              if( ew.kind == EvolutionWindow::ONE && _evolution_dom_var.dag() && _mInpDisc.count( _vInp[upd.pos] ) ){
                auto ite = nel_in.find( _evolution_dom_var );
                if( ite != nel_in.end() ) ite->second = ew.index;
              }
              inc = pos_input( _vInp[upd.pos], nel_in );
              bool zeroit = false;
              if( zero_integrated ){
                auto iti = _mInp.find( _vInp[upd.pos] );
                if( iti != _mInp.end() )
                  for( auto const& dv : iti->second ) if( free_dom.count( dv ) ){ zeroit = true; break; }
              }
              OCinp[upd.pos].update( zeroit ? s4_zeros.data() : &inp[inc], nel_in );
              break;
            }
            case CONSTANT: break;
            default: throw Exceptions( Exceptions::INTERNAL );
          }
        }
        _dag->eval( _sgEqn[ndx], OCwrk,
                    1, &eqnvar, &OCres,
                    _vVar.size(), _vVar.data(), OCvar.data(),
                    _vInp.size(), _vInp.data(), OCinp.data(),
                    _vCst.size(), _vCst.data(), OCcst.data(),
                    _vDom.size(), _vDom.data(), OCdom.data() );
      };

      // section 4 (OpEval): point-reduced directions of this equation.  pt_el pins the resident
      // element to contain each target coordinate; the FFEval operator collapses the direction
      // during the DAG eval, so OCres.dom() comes back without it and the row-write stays trivial.
      std::map<FFVar,size_t,lt_FFVar> pt_el;
      if( ndx < _eqnEvalPoint.size() ){
        for( auto const& [dv, val] : _eqnEvalPoint[ndx] ){
          auto itd = _mDom.find( dv );
          if( itd == _mDom.end() ) throw Exceptions( Exceptions::INDEX );
          int const esd = ( ndx < _eqnEvalSide.size() && _eqnEvalSide[ndx].count( dv ) )? _eqnEvalSide[ndx].at( dv ): (int)FFDom::MINUS;
          auto const [iel, xi] = _point_in_domain( itd->second, val, esd );
          (void)xi;
          pt_el[dv] = iel;
        }
      }

      std::vector<U> buf;
      if( !free_dom.empty() ){
        // The residual of a domain-less equation is  R = g(scalars) - sum_k I_k, where I_k is the
        // FFIntegral contribution of free element k.  Accumulating the whole per-element residual
        // sums the NON-integrated part g once per element:  buf = n_el*g - sum_k I_k.  We correct
        // by subtracting (n_el-1)*g, recovering the true residual g - INT.  EXACT for residuals
        // linear in the integrated terms (e.g. Tlo - INT c dz, scalar = a*INT + b); a nonlinear-
        // in-integral residual is handled later by reduce_order extraction of each FFIntegral into
        // a lumped auxiliary state (section 4d), which also retires this g evaluation.
        // free_el iterates the CONSUMED (free) directions ONLY, starting empty so _next_block sees
        // exactly free_dom.  Each assemble combines the current block ndx_el (distributed directions
        // of THIS equation, e.g. {y}) with the free element, so a DISTRIBUTED integral LINK
        // (w1 - OpI(c,z) on {y}) free-accumulates z within each y block.  (A domain-less equation has
        // ndx_el empty, so asm_el == free_el -- unchanged from the scalar case.)
        std::map<FFVar,size_t,lt_FFVar> free_el;
        bool have_acc = false;
        size_t n_iter = 0;
        while( _next_block( free_el, free_dom ) && !free_el.empty() ){
          std::map<FFVar,size_t,lt_FFVar> asm_el = ndx_el;
          for( auto const& [fdv, fiel] : free_el ) asm_el[fdv] = fiel;
          assemble( asm_el, false );
          if( !have_acc ){
            buf.assign( OCres.coef().begin(), OCres.coef().end() );
            have_acc = true;
          }
          else{
            if( buf.size() != OCres.coef().size() ) throw Exceptions( Exceptions::INTERNAL );
            for( size_t i = 0; i < buf.size(); ++i ) buf[i] += OCres.coef()[i];
          }
          ++n_iter;
        }
        if( !have_acc ) throw Exceptions( Exceptions::INTERNAL );

        if( n_iter > 1 ){
          // g = non-integrated part: re-evaluate with the integrated fields zeroed (FFIntegral -> 0).
          // Any valid free element works (g is independent of the free direction).
          s4_zeros.assign( _nCollVar, U( 0. ) );
          std::map<FFVar,size_t,lt_FFVar> free_el0;
          if( !_next_block( free_el0, free_dom ) || free_el0.empty() )
            throw Exceptions( Exceptions::INTERNAL );
          std::map<FFVar,size_t,lt_FFVar> g_el = ndx_el;   // block dirs + one free element
          for( auto const& [fdv, fiel] : free_el0 ) g_el[fdv] = fiel;
          assemble( g_el, true );
          if( OCres.coef().size() != buf.size() ) throw Exceptions( Exceptions::INTERNAL );
          U const corr = U( static_cast<double>( n_iter - 1 ) );
          for( size_t i = 0; i < buf.size(); ++i ) buf[i] -= corr * OCres.coef()[i];
        }
      }
      else if( !pt_el.empty() ){
        // section 4 (OpEval): POINT reduction.  Pin the target point element(s) into a COPY of the
        // current block ndx_el (leaving ndx_el clean for the next _next_block), then assemble.  The
        // FFEval operator collapses each point direction during the DAG eval, so OCres.dom() comes
        // back without it.  This covers a domain-less scalar (ndx_el empty -> nel = pt_el) AND a
        // DISTRIBUTED equation carrying an OpEval over another direction (e.g. a nested reduction's
        // inner LINK w1 - OpEval(a*b,z,z0) on {y}): nel = {y:iel} + {z:iel_z0}, so z is pinned while
        // the block's distributed directions are preserved and the row-write emits their rows.
        std::map<FFVar,size_t,lt_FFVar> nel = ndx_el;
        for( auto const& [dv, iel] : pt_el ) nel[dv] = iel;
        assemble( nel, false );
        buf.assign( OCres.coef().begin(), OCres.coef().end() );
      }
      else{
        assemble( ndx_el, false );
        buf.assign( OCres.coef().begin(), OCres.coef().end() );
      }
#ifdef MC__OCFESLV_EVAL_DEBUG
      _disp( 3 ) << "COLLOCATED RESIDUAL: " << OCres << std::endl;
#endif

      // -- Apply SAT penalties to buf[] (IC_WEAK / IC_TRACE / IC_STRONG) - //

      if( weak ){
        if( !_plan.weak_sat_terms_ready() ){
          _disp( 2 ) << "OCFESLV::_eval_eqn ** IC_WEAK frozen SAT terms are "
                       "not ready; setup() must build the weak receiver plan "
                       "before residual evaluation." << std::endl;
          throw Exceptions( Exceptions::INTERNAL );
        }

        // IC_WEAK: two-level dispatch via (row_id, BlockKey).
        // Outer lookup by equation index, inner lookup by element-block.
        // Both are O(1) amortised.  No frozen_term_matches_block scan.
        auto const it_outer = _plan.weak_sat_index().find( ndx );
        if( it_outer != _plan.weak_sat_index().end() ){
          OCPlan::BlockKey const bk = project_block_key( ndx_el );
          auto const it_inner = it_outer->second.find( bk );
          if( it_inner != it_outer->second.end() ){
            auto const& terms = _plan.weak_sat_terms();
            for( size_t ti : it_inner->second ){
              FrozenWeakSatTerm const& term = terms[ti];
              if( weak_mult && term.exact_claim ){
                // The claim is imposed EXACTLY: inject its multiplier as IC_TRACE does, and no penalty.
                if( !term.trace_tau_allowed || term.prefactor == 0. ) continue;
                if( term.receiver_flat >= buf.size() || term.lambda_col >= _nCollVar ){
                  _disp( 2 ) << "OCFESLV::_eval_eqn ** IC_WEAK+TAU receiver term out of range: row="
                            << term.row_id << " flat=" << term.receiver_flat << " buf.size=" << buf.size()
                            << " lambda_col=" << term.lambda_col << " nvar=" << _nCollVar << std::endl;
                  throw Exceptions( Exceptions::INTERNAL );
                }
                buf[ term.receiver_flat ] += U( term.trace_prefactor ) * var[ term.lambda_col ];
                continue;
              }
              if( term.sat_coeffs.empty() ) continue;  // prefactor == 0
              if( term.receiver_flat >= buf.size() ){
                _disp( 2 ) << "OCFESLV::_eval_eqn ** frozen weak-SAT receiver index "
                             "is out of range: row=" << term.row_id
                          << " flat=" << term.receiver_flat
                          << " buf.size=" << buf.size() << std::endl;
                throw Exceptions( Exceptions::INTERNAL );
              }
              for( auto const& [col, c] : term.sat_coeffs )
                buf[ term.receiver_flat ] += U( c ) * var[ col ];
            }
          }
        }
      }
      else if( trace ){
        if( !_plan.weak_sat_terms_ready() ){
          _disp( 2 ) << "OCFESLV::_eval_eqn ** IC_TRACE frozen SAT terms are "
                       "not ready; setup() must build the shared receiver plan "
                       "before residual evaluation." << std::endl;
          throw Exceptions( Exceptions::INTERNAL );
        }

        // IC_TRACE: two-level dispatch via (row_id, BlockKey).
        // trace_sat_index contains only trace_tau_allowed terms.
        auto const it_outer = _plan.trace_sat_index().find( ndx );
        if( it_outer != _plan.trace_sat_index().end() ){
          OCPlan::BlockKey const bk = project_block_key( ndx_el );
          auto const it_inner = it_outer->second.find( bk );
          if( it_inner != it_outer->second.end() ){
            auto const& terms = _plan.weak_sat_terms();
            for( size_t ti : it_inner->second ){
              FrozenWeakSatTerm const& term = terms[ti];
              if( term.receiver_flat >= buf.size() ){
                _disp( 2 ) << "OCFESLV::_eval_eqn ** frozen IC_TRACE receiver index "
                             "is out of range: row=" << term.row_id
                          << " flat=" << term.receiver_flat
                          << " buf.size=" << buf.size() << std::endl;
                throw Exceptions( Exceptions::INTERNAL );
              }
              if( term.lambda_col >= _nCollVar ){
                std::cerr << "OCFESLV::_eval_eqn ** frozen IC_TRACE receiver term has "
                             "invalid lambda column: row=" << term.row_id
                          << " lambda_col=" << term.lambda_col
                          << " nvar=" << _nCollVar << std::endl;
                throw Exceptions( Exceptions::INTERNAL );
              }
              if( term.prefactor == 0. ) continue;
              // rev107 [taucol/assembly]: same counter, IC_TRACE side.  IC_TRACE keeps
              // EVERY tau explicit and drops NO physical row, so a k>0 here cannot be
              // an artefact of the keep-explicit hybrid.  That makes IC_TRACE the
              // discriminating configuration and it must be instrumented too.
              if( _taucol_enabled() && _tauColPlanReady && !_tauColAssemblyReported )
                ++_tauInjectCount[ term.lambda_col ];
              buf[ term.receiver_flat ] += U( term.trace_prefactor )
                                         * var[ term.lambda_col ];
            }
          }
        }
      } // end frozen trace multiplier block
      else if( strong && _plan.n_trace_var_explicit() ){
        // IC_STRONG hybrid: the eliminated taus are handled by the Schur
        // projection below, but the explicit (kept) taus are real var[]
        // unknowns and must couple into their receiver rows exactly as in
        // IC_TRACE.  Explicit tau columns live at >= trace_var_offset; all
        // lower lambda columns belong to the eliminated Schur block and are
        // skipped here.
        if( !_plan.weak_sat_terms_ready() ){
          _disp( 2 ) << "OCFESLV::_eval_eqn ** IC_STRONG explicit-tau SAT terms are "
                       "not ready; setup() must build the receiver plan before "
                       "residual evaluation." << std::endl;
          throw Exceptions( Exceptions::INTERNAL );
        }
        auto const it_outer = _plan.trace_sat_index().find( ndx );
        if( it_outer != _plan.trace_sat_index().end() ){
          OCPlan::BlockKey const bk = project_block_key( ndx_el );
          auto const it_inner = it_outer->second.find( bk );
          if( it_inner != it_outer->second.end() ){
            auto const& terms = _plan.weak_sat_terms();
            for( size_t ti : it_inner->second ){
              FrozenWeakSatTerm const& term = terms[ti];
              if( term.lambda_col < _traceVarOffset ) continue; // eliminated tau
              if( term.receiver_flat >= buf.size() ){
                _disp( 2 ) << "OCFESLV::_eval_eqn ** IC_STRONG explicit-tau receiver "
                             "index is out of range: row=" << term.row_id
                          << " flat=" << term.receiver_flat
                          << " buf.size=" << buf.size() << std::endl;
                throw Exceptions( Exceptions::INTERNAL );
              }
              if( term.lambda_col >= _nCollVar ){
                std::cerr << "OCFESLV::_eval_eqn ** IC_STRONG explicit-tau receiver "
                             "term has invalid lambda column: row=" << term.row_id
                          << " lambda_col=" << term.lambda_col
                          << " nvar=" << _nCollVar << std::endl;
                throw Exceptions( Exceptions::INTERNAL );
              }
              if( term.prefactor == 0. ) continue;
              // rev107 [taucol/assembly]: count the injection that ACTUALLY executed.
              // Placed AFTER every guard above, so a column whose term exists in the
              // plan but never reaches this line reads 0 -- which is precisely how
              // outcome (2) "a dispatch that never fires" becomes visible.  Nothing
              // below reads the counter.
              if( _taucol_enabled() && _tauColPlanReady && !_tauColAssemblyReported )
                ++_tauInjectCount[ term.lambda_col ];
              buf[ term.receiver_flat ] += U( term.trace_prefactor )
                                         * var[ term.lambda_col ];
            }
          }
        }
      } // end IC_STRONG explicit-tau multiplier block

      // -- Domain-filter and write buf[] to res[] ---------------------- //
      {
        struct CDim {
          FFDom const* pdom;
          FFVar const* pvar;
          size_t       cur_n;
          int          lim;
          size_t       iel;
        };
        std::vector<CDim> cdims;
        cdims.reserve( OCres.dom().size() );
        for( auto const& [pvar, pdom] : OCres.dom() ){
          int const lim = eqndom.count(*pvar) ? eqndom.at(*pvar) : FFDom::ALL;
          size_t const iel = ndx_el.count(*pvar) ? ndx_el.at(*pvar) : 0;
          cdims.push_back({ pdom, pvar, pdom->n_node, lim, iel });
        }

        // First collapse true physical boundary dimensions by interpolation.
        // Remaining ALL/... dimensions are filtered simultaneously below so
        // that multi-dimensional interface corners can be treated without
        // changing the tensor-product ordering.
        for( size_t idim = 0; idim < cdims.size(); ++idim ){
          auto& cd = cdims[idim];
          if( cd.lim != FFDom::LB && cd.lim != FFDom::UB ) continue;

          size_t stride = 1;
          for( size_t d = 0; d < idim; ++d ) stride *= cdims[d].cur_n;
          size_t outer = 1;
          for( size_t d = idim+1; d < cdims.size(); ++d ) outer *= cdims[d].cur_n;

          double const t_eval = ( cd.lim == FFDom::LB ) ? -1.0 : 1.0;
          std::vector<U> newbuf;
          if( !cd.pdom->eval_lagrange( newbuf, t_eval, 0, buf, stride, outer ) )
            throw Exceptions( Exceptions::INTERNAL );
          buf      = std::move( newbuf );
          cd.cur_n = 1;
        }

        size_t ntot = 1;
        for( auto const& cd : cdims ) ntot *= cd.cur_n;

        std::vector<U> newbuf;
        newbuf.reserve( ntot );
        std::vector<size_t> midx( cdims.size(), 0 );

        // Record which flat tensor indices survive the boundary filter.
        // IC_STRONG and IC_TRACE need the original flat index to build
        // FrozenPhysicalResidualRowKey; storing kept_flat avoids re-running
        // the boundary filter switch in a second pass.
        // rev107 [taucol]: under IC_TRACE with no exact-replacement rows this was false,
        // so LOOP 2 took the bulk emit path, never built a row key, and the instrument
        // recorded nothing -- then reported LOST for six columns that tier0-audit measures
        // as fully supported (empty=0).  A wrong verdict beside a warning is still a wrong
        // verdict.  Ask for the row keys when the instrument is live.  Flag-off the
        // expression is unchanged, and the trace branch emits exactly the values the bulk
        // path emits, in the same order.
        bool const need_kept_flat = strong
          || ( trace && ( !_plan.exact_replacement_rows().empty()
                       || ( _taucol_enabled() && _tauColPlanReady ) ) )
          || build_prov;   // the node half of the row provenance needs the surviving flat indices
        // the node a surviving flat index stands at, into the last provenance record
        auto node_rec = [&]( size_t flat, t_RowNode& rec ){
          rec.dims.clear();  rec.dims.reserve( cdims.size() );
          size_t rem = flat;
          for( size_t d = 0; d < cdims.size(); ++d ){
            size_t const m = rem % cdims[d].cur_n;  rem /= cdims[d].cur_n;
            rec.dims.push_back({ cdims[d].pvar->id().second, cdims[d].iel, m, cdims[d].lim });
          }
        };
        auto prov_node = [&]( size_t flat ){ if( build_prov && !_rowNode.empty() ) node_rec( flat, _rowNode.back() ); };
        std::vector<size_t> kept_flat;
        if( need_kept_flat ) kept_flat.reserve( ntot );

        for( size_t flat = 0; flat < ntot; ++flat ){
          size_t rem = flat;
          for( size_t d = 0; d < cdims.size(); ++d ){
            midx[d] = rem % cdims[d].cur_n;
            rem    /= cdims[d].cur_n;
          }

          bool base_keep = true;
          for( size_t d = 0; d < cdims.size() && base_keep; ++d ){
            auto const& cd = cdims[d];
            if( cd.lim == FFDom::LB || cd.lim == FFDom::UB ) continue;
            switch( cd.lim ){
              case FFDom::ALL:
                break;
              case FFDom::ALL-FFDom::LB:
                if( cd.iel == 0 && midx[d] == 0 ) base_keep = false;
                if( _causal_start_excludes( static_cast<size_t>( &eqnrec - &_mEqn[0] ), *cd.pvar, cd.iel, midx[d] ) ) base_keep = false;   // causal start: initial node
                break;
              case FFDom::ALL-FFDom::UB:
                if( cd.iel + 1 == cd.pdom->n_elem
                 && midx[d] + 1 == cd.pdom->n_node ) base_keep = false;
                break;
              case FFDom::ALL-FFDom::LB-FFDom::UB:
                if( cd.iel == 0 && midx[d] == 0 ) base_keep = false;
                if( _causal_start_excludes( static_cast<size_t>( &eqnrec - &_mEqn[0] ), *cd.pvar, cd.iel, midx[d] ) ) base_keep = false;
                if( cd.iel + 1 == cd.pdom->n_elem
                 && midx[d] + 1 == cd.pdom->n_node ) base_keep = false;
                break;
              default:
                throw Exceptions( Exceptions::DOMAIN );
            }
          }

          if( base_keep ){
            newbuf.push_back( buf[flat] );
            if( need_kept_flat ) kept_flat.push_back( flat );
          }
        }

        if( strong ){
          // IC_STRONG: tag pivot, exact-replacement, and ordinary rows.
          auto const& proj = _plan.strong_projection();
          auto const& repl = _plan.exact_replacement_rows();
          for( size_t i = 0; i < kept_flat.size(); ++i ){
            FrozenPhysicalResidualRowKey const rk
              = make_strong_row_key( ndx, ndx_el, kept_flat[i] );
            int piv_idx = -1;
            std::vector<std::pair<size_t,double>> const* p_corr = nullptr;
            auto it_piv = proj.pivot_index.find( rk );
            if( it_piv != proj.pivot_index.end() ){
              piv_idx = static_cast<int>( it_piv->second );
            }
            else if( !repl.empty() && repl.find( rk ) != repl.end() ){
              piv_idx = -2;  // exact-replacement — will be dropped
            }
            else{
              auto it_corr = proj.projection.find( rk );
              if( it_corr != proj.projection.end() ) p_corr = &it_corr->second;
            }
            // ---- rev107 [taucol/assembly] ------------------------------------
            // The observation the handoff asked for, at the one place where a row's
            // fate is actually decided.  Two things are recorded, both read-only:
            //
            //  (i) if this row hosts an explicit tau (per _tauHostRows, written by the
            //      PLAN half and only read here), which of the three fates it takes;
            //  (ii) for every SURVIVING row, which pivot indices its Schur correction
            //      references.  A pivot row is DROPPED but its value is subtracted into
            //      each corrected row, so "dropped" and "lost" are NOT the same thing.
            //      Without (ii) the instrument could not tell them apart, and the whole
            //      distinction between outcome (3) and a genuine orphan would be an
            //      inference rather than a measurement.
            if( _taucol_enabled() && _tauColPlanReady && !_tauColAssemblyReported ){
              auto ith = _tauHostRows.find( rk );
              if( ith != _tauHostRows.end() )
                for( auto const& cc : ith->second ){
                  size_t const lc = cc.first;
                  if     ( piv_idx >= 0  ) ++_tauFatePivot[ lc ];
                  else if( piv_idx == -2 ) ++_tauFateExact[ lc ];
                  else                     ++_tauFateEmit [ lc ];
                }
              if( piv_idx == -1 && p_corr )
                for( auto const& pc : *p_corr )
                  if( std::fabs( pc.second ) > 64. * DBL_EPSILON )
                    ++_tauCorrRefs[ pc.first ];
              // NET DELIVERED COEFFICIENT.  The row this loop is tagging will be
              // emitted as  R_k - sum_p y_p * R_p.  Every tau appearing in R_k or in
              // any referenced pivot row therefore arrives with coefficient
              //     c_k(tau) - sum_p y_p * c_p(tau)
              // and that sum can be ZERO even though the row survives and the pivot is
              // redistributed into it.  Measuring survival alone reports SUPPORTED for a
              // column tier0-audit measures as numerically EMPTY -- the disagreement that
              // refuted the first reading of this instrument.  Measure the number.
              if( piv_idx == -1 ){
                std::map<size_t,double> net;
                if( ith != _tauHostRows.end() )
                  for( auto const& cc : ith->second ) net[ cc.first ] += cc.second;
                if( p_corr )
                  for( auto const& pc : *p_corr ){
                    auto ip = _tauPivotCoef.find( pc.first );
                    if( ip == _tauPivotCoef.end() ) continue;
                    for( auto const& cv : ip->second )
                      net[ cv.first ] -= pc.second * cv.second;
                  }
                for( auto const& nv : net ){
                  double& m = _tauNetMax[ nv.first ];
                  m = std::max( m, std::fabs( nv.second ) );
                }
              }
            }
            strong_rows.push_back({ newbuf[i], piv_idx, p_corr, (int)ndx, {} });
            if( build_prov ) node_rec( kept_flat[i], strong_rows.back().node );
          }
        }
        else if( trace && !kept_flat.empty() ){
          // IC_TRACE: iterate kept_flat to drop exact-replacement rows.
          for( size_t i = 0; i < kept_flat.size(); ++i ){
            FrozenPhysicalResidualRowKey const rk
              = make_strong_row_key( ndx, ndx_el, kept_flat[i] );
            bool const dropped = ( _plan.exact_replacement_rows().find( rk )
                                != _plan.exact_replacement_rows().end() );
            // rev107 [taucol/assembly]: IC_TRACE has no Schur pivots, but it DOES drop
            // exact-replacement rows.  Recording both fates here keeps the two modes
            // reported on the same footing, so an IC_TRACE run can be read against an
            // IC_STRONG run without a second convention.
            if( _taucol_enabled() && _tauColPlanReady && !_tauColAssemblyReported ){
              auto ith = _tauHostRows.find( rk );
              if( ith != _tauHostRows.end() )
                for( auto const& cc : ith->second ){
                  size_t const lc = cc.first;
                  if( dropped ){ ++_tauFateExact[ lc ]; }
                  else{
                    ++_tauFateEmit[ lc ];
                    double& m = _tauNetMax[ lc ];       // IC_TRACE: no Schur, so the
                    m = std::max( m, std::fabs( cc.second ) );  // coefficient IS the net
                  }
                }
            }
            if( dropped ){
              ++n_audit_drop;
              continue;
            }
            emit_eqn( newbuf[i], (int)ndx );
            prov_node( kept_flat[i] );
          }
        }
        else if( build_prov ){   // the bulk path, on the first run: one record per row
          for( size_t i = 0; i < newbuf.size(); ++i ){ emit_eqn( newbuf[i], (int)ndx ); prov_node( kept_flat[i] ); }
        }
        else{
          for( U const& val : newbuf ) emit_eqn( val, (int)ndx );
        }
      }
    }

    ++ndx;
  }

  if( strong ){
    if( !_plan.strong_projection().ready ){
      _disp( 2 ) << "OCFESLV::_eval_eqn ** IC_STRONG tau-eliminated projection "
                   "is not ready; setup() must build the frozen Schur plan."
                << std::endl;
      throw Exceptions( Exceptions::INTERNAL );
    }

    size_t const ntrace = _plan.strong_projection().n_trace;

    // Pass 1: extract pivot values into a dense array.
    std::vector<U> pivot_vals( ntrace, U(0.) );
    std::vector<StrongRow const*> pivot_row( build_prov? ntrace: 0, nullptr );   // provenance: the pivot rows themselves
    for( auto const& sr : strong_rows ){
      if( sr.pivot_idx >= 0 ){
        if( static_cast<size_t>( sr.pivot_idx ) >= ntrace )
          throw Exceptions( Exceptions::INTERNAL );
        pivot_vals[ static_cast<size_t>( sr.pivot_idx ) ] = sr.value;
        if( build_prov ) pivot_row[ static_cast<size_t>( sr.pivot_idx ) ] = &sr;
      }
    }

    // Pass 2: emit non-pivot, non-exact-replacement rows with Schur correction.
    for( auto const& sr : strong_rows ){
      if( sr.pivot_idx >= 0 ){ ++n_strong_pivot_drop; continue; }
      if( sr.pivot_idx == -2 ){ ++n_strong_exact_drop; continue; }
      U val = sr.value;
      if( sr.correction ){
        for( auto const& [pidx, coeff] : *sr.correction )
          val -= U( coeff ) * pivot_vals[ pidx ];
      }
      emit_eqn( val, sr.ndx );
      if( build_prov && !_rowNode.empty() ){
        _rowNode.back() = sr.node;
        if( sr.correction )
          for( auto const& [pidx, coeff] : *sr.correction )
            if( pidx < pivot_row.size() && pivot_row[pidx] )
              _rowFold.back().push_back({ pivot_row[pidx]->ndx, pivot_row[pidx]->node, coeff });
      }
    }
  }

  // ---- Row-count self-check (Part 2) -----------------------------------
  // Frozen exact interface rows plus surviving physical equation nodes must
  // jointly emit exactly _nCollEqn residual rows.  For IC_STRONG, the number
  // of exact trace rows must match the number of pivot physical rows removed
  // by the frozen Schur projection.  Ordinary LB/UB domain-limit drops are
  // not counted in n_audit_drop.  IC_STRONG keeps its projection pivot count
  // separate from the IC_TRACE exact-replacement audit counter to avoid
  // reusing one variable for two structurally different meanings.
  // Both identities are purely index-based;
  // a violation means the interface/corner logic has made the collocation
  // system non-square.  The diagnostic is emitted on cerr first because
  // eval()/deriv() turn the exception into a bare false return.
  // 2026-10-02 (row-level causal design): the element-start continuity rows -- one per state with an evolution
  // derivative and per (element start, transverse node): the state at the element's first time node minus the
  // previous element's trace there (eval_colloc at t_k^-, side MINUS), as a marching window's start.  Counted in
  // _check_consistency(); they take the slots of the `ALL - LB` rows not collocated at those nodes.
  if( _causal_rows_active() ){
    FFDom const& edom = _mDom.at( _evolution_dom_var );
    // Lagrange weights of the time nodes at an element's RIGHT end (xi = +1): the previous element's trace at t_k
    std::vector<double> w( edom.n_node, 0. );
    for( size_t j = 0; j < edom.n_node; ++j ){
      std::vector<double> e( edom.n_node, 0. );  e[j] = 1.;  double wj = 0.;
      if( !edom.eval_lagrange( wj, 1., 0, e ) ) throw Exceptions( Exceptions::INTERNAL );
      w[j] = wj;
    }
    for( auto const& st : _states_with_evolution_derivative( _mEqn, _mVar ) ){
      auto itv = _mVar.find( st );
      if( itv == _mVar.end() ) continue;
      FFVar const& sv = itv->first;
      std::vector<FFVar> const dep( itv->second.begin(), itv->second.end() );   // set order == block layout order
      size_t evo_pos = std::numeric_limits<size_t>::max();
      std::vector<FFDom const*> tdd;  std::vector<size_t> tpos;
      for( size_t k = 0; k < dep.size(); ++k ){
        if( dep[k].id() == _evolution_dom_var.id() ){ evo_pos = k; continue; }
        tdd.push_back( &_mDom.at( dep[k] ) );  tpos.push_back( k );
      }
      if( evo_pos == std::numeric_limits<size_t>::max() ) continue;
      size_t ncombo = 1;
      for( auto const* d2 : tdd ) ncombo *= d2->n_elem * d2->n_node;
      for( size_t iel = 1; iel < edom.n_elem; ++iel ){
        // A TRANSITION at this element start that makes the state JUMP: its post-jump value is the lifted auxiliary
        // w (defined by the transition's LINK row from x(tau^-)), so the row is x(t_k^+) - w, not continuity -- as
        // the claims path does by replacing the claim's low side with w (_lifted_claim_aux).
        FFVar const* waux = nullptr;
        { double const tb = edom.elem_lo( iel );
          for( auto const& L : _trnLift )
            if( L.state.id() == sv.id() && std::fabs( L.tau - tb ) <= 64.*DBL_EPSILON*std::max( 1., std::fabs( tb ) ) )
              waux = &L.aux; }
        std::vector<FFVar> wdep;
        if( waux ){
          auto itw = _mVar.find( *waux );
          if( itw == _mVar.end() ) throw Exceptions( Exceptions::INTERNAL );
          wdep.assign( itw->second.begin(), itw->second.end() );
        }
        for( size_t c = 0; c < ncombo; ++c ){
          std::map<FFVar,size_t,lt_FFVar> ndx_hi, ndx_lo;
          ndx_hi[ _evolution_dom_var ] = iel;  ndx_lo[ _evolution_dom_var ] = iel - 1;
          std::vector<size_t> tnode( tdd.size() );
          size_t rem = c;
          for( size_t k = 0; k < tdd.size(); ++k ){
            size_t const nn = tdd[k]->n_node, ne = tdd[k]->n_elem, idx = rem % ( ne * nn );
            rem /= ( ne * nn );
            ndx_hi[ dep[ tpos[k] ] ] = ndx_lo[ dep[ tpos[k] ] ] = idx / nn;  tnode[k] = idx % nn;
          }
          // local flat index of (time node jt, transverse nodes): dims in block order, stride = product of earlier n_node
          auto flat_of = [&]( size_t jt ){
            size_t flat = 0, stride = 1, kt = 0;
            for( size_t k = 0; k < dep.size(); ++k ){
              size_t const nn = _mDom.at( dep[k] ).n_node;
              size_t const in = ( k == evo_pos )? jt: tnode[kt++];
              flat += in * stride;  stride *= nn;
            }
            return flat;
          };
          size_t const inc_hi = pos_state( sv, ndx_hi );
          U r = var[ inc_hi + flat_of( 0 ) ];
          if( waux ){
            // w's column: its own domains (the state's transverse ones), at the same transverse nodes
            std::map<FFVar,size_t,lt_FFVar> ndx_w;  size_t flat = 0, stride = 1;
            for( auto const& dv : wdep ){
              size_t kt = 0;  while( kt < tpos.size() && dep[ tpos[kt] ].id() != dv.id() ) ++kt;
              if( kt == tpos.size() ) throw Exceptions( Exceptions::INTERNAL );
              ndx_w[ dv ] = ndx_hi.at( dep[ tpos[kt] ] );
              flat += tnode[kt] * stride;  stride *= _mDom.at( dv ).n_node;
            }
            r -= var[ pos_state( *waux, ndx_w ) + flat ];
          }
          else{
            size_t const inc_lo = pos_state( sv, ndx_lo );
            for( size_t j = 0; j < edom.n_node; ++j ) r -= w[j] * var[ inc_lo + flat_of( j ) ];
          }
          emit_eqn( r, -1 );
        }
      }
    }
  }

  {
    size_t const n_total = eqn0 ? static_cast<size_t>( eqn - eqn0 ) : 0;
    bool const trace_rows = ( trace || weak_mult );   // IC_WEAK with multipliers follows the IC_TRACE identity
    size_t const expected_trace_drop = ( trace_rows && _nTraceConstraint >= _nTraceVar )
                                           ? ( _nTraceConstraint - _nTraceVar )
                                           : size_t(0);
    size_t const n_effective_audit_drop = strong
                                        ? ( n_strong_pivot_drop + n_strong_exact_drop )
                                        : n_audit_drop;
    // IC_STRONG hybrid: explicit (kept) tau continuity rows are ADDED to the
    // system (LOOP 1) without dropping a physical row (LOOP 2), since their
    // LINK row is retained.  So LOOP 1 emissions exceed LOOP 2 drops by exactly
    // the explicit-tau count.  Zero for every non-hybrid plan.
    size_t const n_explicit_added = _plan.n_trace_var_explicit();

    // rev40 (Stage 2 pre-flight): shape the identity for the beta term NOW, hard-zero, so
    // that the arithmetic being validated on the corpus is the arithmetic Stage 2 will use.
    //
    // Under Stage 2 the pivot search stops at r = rank(B) rather than m, so LOOP 2 drops
    // fall by k while LOOP 1 emissions are unchanged; squareness is restored by k explicit
    // beta unknowns and their k rows.  The identity therefore acquires a "+ k" term.
    // In rev40 this is 0 by construction and the condition is arithmetically identical to
    // rev39 -- only the SHAPE is under test.  rev41 replaces the literal with
    // _plan.n_trace_var_beta().
    //
    // Run this on the corpus before changing the pivot loop: if it fires anywhere today,
    // the identity is already violated somewhere that has not been found, and that must be
    // understood before a term is added to it.
    size_t const n_beta_added = 0;   // rev41: _plan.n_trace_var_beta()

    bool const identity_ok =
         ( n_total == _nCollEqn )
      && ( trace_rows || n_audit_emit == n_effective_audit_drop + n_explicit_added + n_beta_added )
      && ( !trace_rows || n_effective_audit_drop == expected_trace_drop );

    // rev107 [taucol/assembly]: the delivery half.  Emitted here because this is the
    // one point that is reached on EVERY eval, unconditionally, for every imposition
    // mode -- the same property that makes [stage2-preflight] a reliable build oracle.
    // One-shot per OCFESLV instance.
    if( _taucol_enabled() && _tauColPlanReady && !_tauColAssemblyReported ){
      _tauColAssemblyReported = true;   // set FIRST: the recorders above test it, and a
                                        // throw below must not leave them accumulating.
      auto const& proj = _plan.strong_projection();
      _disp( 3 ) << "\n#### [taucol/assembly] delivered coupling, after Schur substitution"
                   " (CRONOS_AUDIT_TAUCOL) ####\n";
      _disp( 3 ) << "  mode=" << ( strong ? "IC_STRONG" : trace ? "IC_TRACE" : "IC_WEAK" )
                << "  rows_emitted=" << n_total
                << "  loop1_emits=" << n_audit_emit
                << "  loop2_drops=" << n_effective_audit_drop
                << " (pivot=" << n_strong_pivot_drop
                << " exact=" << n_strong_exact_drop << ")"
                << "  host_rows_named_by_plan=" << _tauHostRows.size() << "\n";
      if( _tauHostRows.empty() )
        _disp( 3 ) << "  [taucol/assembly] the plan half named no host rows"
                     " -- either it did not run (is CRONOS_AUDIT_TAUCOL set before setup?)"
                     " or this plan has no explicit multiplier columns.\n";
      if( trace )
        _disp( 3 ) << "  [taucol/assembly] NOTE: IC_TRACE drops no pivot row, so every host"
                     " row survives and the net delivered coefficient is simply the term's"
                     " own coefficient.  A k>0 here therefore CANNOT be an artefact of the"
                     " IC_STRONG keep-explicit hybrid: it is a rank property of B.\n";
      size_t n_lost = 0, n_kept = 0;
      std::set<size_t> cols;
      for( auto const& kv : _tauHostRows )
        for( auto const& cc : kv.second ) cols.insert( cc.first );
      for( size_t c : cols ){
        size_t const inj = _tauInjectCount.count(c) ? _tauInjectCount.at(c) : 0;
        size_t const fe  = _tauFateEmit .count(c) ? _tauFateEmit .at(c) : 0;
        size_t const fp  = _tauFatePivot.count(c) ? _tauFatePivot.at(c) : 0;
        size_t const fx  = _tauFateExact.count(c) ? _tauFateExact.at(c) : 0;
        size_t redistributed = 0;
        std::ostringstream pv;
        auto itp = _tauPivotOfCol.find( c );
        if( itp != _tauPivotOfCol.end() )
          for( size_t p : itp->second ){
            size_t const nref = _tauCorrRefs.count(p) ? _tauCorrRefs.at(p) : 0;
            redistributed += nref;
            pv << " p" << p << "->" << nref << " corrected row(s)";
          }
        double const netmax = _tauNetMax.count(c) ? _tauNetMax.at(c) : 0.0;
        bool const cancels = ( netmax <= 1.0e-12 );
        char const* verdict =
            ( inj == 0 )                 ? "NOT INJECTED -- outcome (2): the dispatch never fired"
          : ( ( fe > 0 || redistributed > 0 ) && cancels )
                                         ? "CANCELLED -- the row survives and the pivot is"
                                           " redistributed into it, but the two contributions"
                                           " sum to zero: the column reaches the assembled"
                                           " system with NO coefficient anywhere"
          : ( fe > 0 )                   ? "SUPPORTED -- a host row survives with net coefficient"
          : ( redistributed > 0 )        ? "REDISTRIBUTED -- host row dropped, but its value"
                                           " reaches surviving rows through the Schur correction"
          :                                "LOST -- outcome (3): every host row dropped and no"
                                           " surviving row references its pivot";
        _disp( 3 ) << "  c" << c
                  << ": injected=" << inj
                  << "  host rows: emitted=" << fe << " pivot-dropped=" << fp
                  << " exact-dropped=" << fx
                  << "  schur:" << ( pv.str().empty() ? std::string(" (none)") : pv.str() )
                  << "  max|net delivered coeff|=" << std::scientific
                  << std::setprecision(3) << netmax
                  << "\n      => " << verdict << "\n";
        if( ( fe > 0 || redistributed > 0 ) && !cancels ) ++n_kept; else ++n_lost;
      }
      _disp( 3 ) << "  [taucol/assembly] summary: columns traced=" << cols.size()
                << "  supported=" << n_kept << "  lost=" << n_lost
                << "  proj.rank=" << proj.rank << "  n_trace_var=" << _plan.n_trace_var()
                << "\n  COMPARE: 'lost' against [tier0-audit] numerically EMPTY columns and"
                   " against k from [evoface] PLAN METRIC.  They are computed by three"
                   " independent paths; if they disagree, believe the one that reads the"
                   " assembled matrix and re-open this instrument.\n";
      _disp( 3 ) << "########################################################\n";
    }

    // One line per OCFESLV instance, unconditional: a corpus sweep produces the (r, m, k)
    // table with no new gating, and a log without this line is a stale binary.  Emitted on
    // IC_STRONG only -- IC_TRACE is Stage 1 and complete, IC_WEAK has no multipliers.
    if( strong && !_stage2PreflightReported ){
      _stage2PreflightReported = true;
      _disp( 3 ) << "  [stage2-preflight] rows: r=" << _traceProjRank
                << " m_all_tau=" << _traceProjNTau
                << " k=" << _traceProjK
                << " n_explicit=" << n_explicit_added
                << " emits=" << n_audit_emit
                << " drops=" << n_effective_audit_drop
                << " (pivot=" << n_strong_pivot_drop
                << " exact=" << n_strong_exact_drop << ")"
                << " beta=" << n_beta_added
                << " identity=" << ( identity_ok ? "OK" : "** FAIL **" )
                << std::endl;
    }

    if( !identity_ok ){
      size_t const n_keep = ( n_total >= n_audit_emit )
                          ? ( n_total - n_audit_emit ) : 0;
      _disp( 2 ) << "OCFESLV::_eval_eqn ** row-count imbalance: "
                << n_audit_emit << " interface row(s) + "
                << n_keep << " equation row(s) = " << n_total
                << " (expected _nCollEqn = " << _nCollEqn << "); "
                << "LOOP 2 interface-dropped " << n_effective_audit_drop << " node(s).  "
                << "The interface/corner logic has made "
                << "the collocation system non-square for the current "
                << "INTERFACE.TYPE/INTERFACE.IMPOSITION/classification."
                << std::endl;
      throw Exceptions( Exceptions::INTERNAL );
    }
  }

  // Provenance built (and row-count verified): freeze it for reuse across evals.
  if( build_prov ) _rowProvenanceReady = true;
}

// ======================================================================
// OCFESLV::_audit_strong_interface_rows
// ======================================================================
//
// Dry structural pass of the residual assembler on dummy FFDep coefficients.
// The row-count invariants enforced at the tail of _eval_eqn() are purely
// index-based: the corner predicates and node filters depend on element/
// node indices, not on coefficient values.  So the dummy values are
// immaterial -- the pass exists only to make a squareness violation surface
// here, at setup()/classify_pde() time, rather than at the first eval()/
// deriv() call.
//
// FFDep (vs. double) is used to avoid any chance that floating-point domain
// errors in the user residual (e.g. log(0)) cause a spurious setup failure:
// FFDep arithmetic propagates dependency structure, never nan/inf.  The
// pass reuses the FFDep evaluation workspaces, which are scratch and are
// refreshed by every subsequent real evaluation.
//
inline bool
OCFEBASE::_audit_strong_interface_rows
()
{
  if( !_nCollEqn ) return true;

  if( !_plan.frozen() || _plan.state() != OCPlan::State::FROZEN ){
    _disp( 2 ) << "OCFESLV::_audit_strong_interface_rows ** frozen interface plan "
                 "is not sealed; setup() must complete plan construction before "
                 "the structural residual audit." << std::endl;
    return false;
  }

  try{
    std::vector<FFDep> dvar( _nCollVar );
    std::vector<FFDep> dinp( _nCollInp );
    std::vector<FFDep> deqn( _nCollEqn );

    _eval_eqn( deqn.data(),
               dvar.empty() ? nullptr : dvar.data(),
               dinp.empty() ? nullptr : dinp.data(),
               _OCDepEqn, _OCDepVar, _OCDepInp, _OCDepDom, _OCDepCst,
               _wrkOCDep );
  }
  catch( Exceptions& e ){
    std::cerr << "OCFESLV::_audit_strong_interface_rows ** structural residual "
                 "pass failed (" << e.what() << "); the collocation system is "
                 "not square for the current INTERFACE.TYPE/INTERFACE.IMPOSITION/"
                 "classification." << std::endl;
    return false;
  }
  catch( ... ){
    std::cerr << "OCFESLV::_audit_strong_interface_rows ** structural residual "
                 "pass failed (unknown exception); the collocation system is "
                 "not square for the current INTERFACE.TYPE/INTERFACE.IMPOSITION/"
                 "classification." << std::endl;
    return false;
  }

  return true;
}

inline bool
OCFEBASE::_next_block
( std::map<FFVar, size_t, lt_FFVar>& ndx_el, std::set< FFVar, lt_FFVar > const& var_dom )
const
{
  // First block
  if( ndx_el.empty() ){
    for( auto& var : var_dom )
      ndx_el[var] = 0;
    return true;
  }
  
  // Next block
  if( ndx_el.size() != var_dom.size() ) return false;

  // Check limits to decide increment
  for( auto& [var, blk] : ndx_el ){
    auto const& dom = _mDom.at(var);
    ++blk;
    if( blk < dom.n_elem ) 
      return true;  // can still increment in current dimension
    blk = 0;        // upper limit for current dimension
  }

  // All blocks visited at this point
  ndx_el.clear();
  return true;
}

inline bool
OCFEBASE::_next_block
( std::map<FFVar, size_t, lt_FFVar>& ndx_el, std::map< FFVar, int, lt_FFVar > const& eqn_dom )
const
{
  // First block
  if( ndx_el.empty() ){
    // Check limits to decide starting point
    for( auto const& [var, bnd] : eqn_dom ){
      auto const& dom = _mDom.at(var);
      switch( bnd ){
        case FFDom::ALL:
        case FFDom::ALL-FFDom::LB:
        case FFDom::ALL-FFDom::UB:
        case FFDom::ALL-FFDom::LB-FFDom::UB:
        case FFDom::LB:                      ndx_el[var] = 0;            break;
        case FFDom::UB:                      ndx_el[var] = dom.n_elem-1; break;
        default: throw Exceptions( Exceptions::INTERNAL );
      }
    }
    return true;
  }
  
  // Next block
  if( ndx_el.size() != eqn_dom.size() ) return false;
  
  // Check limits to decide increment
  for( auto& [var, blk] : ndx_el ){
    auto const& dom = _mDom.at(var);
    switch( eqn_dom.at(var) ){
      case FFDom::ALL:
      case FFDom::ALL-FFDom::LB:
      case FFDom::ALL-FFDom::UB:
      case FFDom::ALL-FFDom::LB-FFDom::UB:
        ++blk;
        if( blk < dom.n_elem ) return true;
        blk = 0;               break;
      case FFDom::LB:
      case FFDom::UB:          break;
      default: throw Exceptions( Exceptions::INTERNAL );
    }
  }

  // All blocks visited at this point
  ndx_el.clear();
  return true;
}

inline std::pair<size_t,double>
OCFEBASE::_point_in_domain
( FFDom const& dom, double const& val, int const side )
const
{
  double const tol = 64. * DBL_EPSILON
                   * std::max( 1., std::max( std::fabs(dom.lo_dom), std::fabs(dom.up_dom) ) );
  if( val < dom.lo_dom - tol || val > dom.up_dom + tol ){
    // During a SOLVE.MARCHING block solve the evolution domain is one sliding element, so an
    // output point outside the current window (e.g. a terminal point on an early block) is
    // expected: clamp here and let the block that actually contains the point produce the
    // value (the marching accumulator overwrites terminal outputs, discarding this clamp).
    if( _evolution_window().kind != EvolutionWindow::ONE )   // not one element of a longer horizon
      throw Exceptions( Exceptions::INDEX );
  }

  double x = std::min( std::max( val, dom.lo_dom ), dom.up_dom );
  size_t iel = 0;
  if( x >= dom.up_dom )
    iel = dom.n_elem - 1;
  else{
    auto it = std::upper_bound( dom.elem_bnd.begin(), dom.elem_bnd.end(), x );
    iel = ( it == dom.elem_bnd.begin() ) ? 0
        : std::min( static_cast<size_t>( std::distance( dom.elem_bnd.begin(), it ) - 1 ),
                    dom.n_elem - 1 );
  }
  // SIDE (2026-09-28): upper_bound gives the element STARTING at an interior boundary (tau^+).  An explicit
  // FFDom::MINUS -- the default of every USER evaluation (FFEval, point outputs), which always passes its side --
  // takes the element ENDING there (tau^-, xi = +1); FFDom::PLUS the one starting there (xi = -1), also when x sits
  // a rounding error below that boundary.  side = -1 (internal callers: quadrature nodes, reference points, which
  // mean THE ELEMENT BEING PROCESSED) keeps the historical upper_bound location unchanged.
  if( side == FFDom::MINUS ){
    if( iel > 0 && std::fabs( x - dom.elem_bnd[iel] ) <= tol ) --iel;
  }
  else if( side == FFDom::PLUS && iel+1 < dom.n_elem && std::fabs( x - dom.elem_bnd[iel+1] ) <= tol ) ++iel;
  // (side == -1, legacy: NO adjustment -- this branch used to run for it too, moving a point a rounding error below
  //  an interior boundary into the next element: a round-off change in drivers that never ask for a side)

  double const lo_el = dom.elem_lo( iel );
  double const w_el  = dom.elem_width( iel );
  double xi = 2. * ( x - lo_el ) / w_el - 1.;
  if( xi < -1. && xi > -1. - 128.*DBL_EPSILON ) xi = -1.;
  if( xi >  1. && xi <  1. + 128.*DBL_EPSILON ) xi =  1.;
  return { iel, xi };
}

inline std::vector<double>
OCFEBASE::_node_block
( FFDom const& dom, size_t const iel )
const
{
  if( iel >= dom.n_elem ) throw Exceptions( Exceptions::INDEX );
  return dom.lgnodes( dom.elem_lo( iel ), dom.elem_up( iel ) );
}

inline std::vector<double>
OCFEBASE::_node_block
( FFVar const& var, std::map<FFVar, size_t, lt_FFVar> const& ndx_el )
const
{
  auto itdom = _mDom.find( var );
  if( itdom == _mDom.end() ) throw Exceptions( Exceptions::INDEX );
  size_t const iel = ndx_el.at(var);
  auto itcache = _nodeBlockCache.find( var );
  if( itcache != _nodeBlockCache.end() && iel < itcache->second.size() )
    return itcache->second[iel];
  return _node_block( itdom->second, iel );
}

inline FFDom const&
OCFEBASE::_input_dom
( FFVar const& inp, FFVar const& dom )
const
{
  auto itdisc = _mInpDisc.find( inp );
  if( itdisc != _mInpDisc.end() ){
    auto itdom = itdisc->second.find( dom );
    if( itdom != itdisc->second.end() )
      return itdom->second;
  }

  auto itdom = _mDom.find( dom );
  if( itdom == _mDom.end() ) throw Exceptions( Exceptions::INDEX );
  return itdom->second;
}

template <typename U>
inline typename OCVar<U>::t_Dom
OCFEBASE::_input_ocdom
( FFVar const& inp, std::set<FFVar,lt_FFVar> const& dep )
const
{
  typename OCVar<U>::t_Dom out;
  for( auto const& dvar : dep ){
    auto itdom = _mDom.find( dvar );
    if( itdom == _mDom.end() ) throw Exceptions( Exceptions::INDEX );
    out[&itdom->first] = &_input_dom( inp, itdom->first );
  }
  return out;
}

inline size_t
OCFEBASE::pos_state
( FFVar const& var, std::map<FFVar, size_t, lt_FFVar> const& ndx_el )
const
{
  // Locate var in _mVar and _blkVar
  auto itvar = _mVar.find( var );
  if( itvar == _mVar.end() ) throw Exceptions( Exceptions::INDEX );
  auto const& var_dom = itvar->second;

  auto itblk = _blkVar.find( var );
  if( itblk == _blkVar.end() ) throw Exceptions( Exceptions::INDEX );
  auto const& [blk_pos, blk_tot, blk_dim] = itblk->second;

  if( blk_dim.size() != var_dom.size() ) throw Exceptions( Exceptions::INDEX );

  // Iterate to target block
  std::size_t pos = blk_pos, stride = blk_tot;
  for( auto const& v : var_dom ){
    size_t const iel = ndx_el.at(v);
    auto const& d = _mDom.at(v);
    if( iel >= d.n_elem ) throw Exceptions( Exceptions::INDEX );
    pos += iel * stride;
    stride *= d.n_elem;
  }

  return pos;
}

inline size_t
OCFEBASE::pos_input
( FFVar const& var, std::map<FFVar, size_t, lt_FFVar> const& ndx_el )
const
{
  // Locate var in _mInp and _blkInp
  auto itvar = _mInp.find( var );
  if( itvar == _mInp.end() ) throw Exceptions( Exceptions::INDEX );
  auto const& var_dom = itvar->second;

  auto itblk = _blkInp.find( var );
  if( itblk == _blkInp.end() ) throw Exceptions( Exceptions::INDEX );
  auto const& [blk_pos, blk_tot, blk_dim] = itblk->second;

  if( blk_dim.size() != var_dom.size() ) throw Exceptions( Exceptions::INDEX );

  // Iterate to target block
  std::size_t pos = blk_pos, stride = blk_tot;
  for( auto const& v : var_dom ){
    size_t const iel = ndx_el.at(v);
    auto const& d = _input_dom( var, v );
    if( iel >= d.n_elem ) throw Exceptions( Exceptions::INDEX );
    pos += iel * stride;
    stride *= d.n_elem;
  }

  return pos;
}

inline bool
OCFEBASE::_fct_is_evolution_integral
( FFVar const& fctvar )
const
{
  if( !_evolution_dom_set || !_evolution_dom_var.dag() ) return false;
  FFVar expr = fctvar;
  FFSubgraph sg = _dag->subgraph( 1, &expr );
  long const evoid = _evolution_dom_var.id().second;
  for( auto const* op : sg.l_op ){
    if( !op || !op->sameid( typeid( FFIntegral ) ) ) continue;
    for( auto const& [v,e] : mc::type_cast<FFIntegral const>( op )->Indep().expr ){
      (void)e;
      if( v.id().second == evoid ) return true;   // integrates over the evolution variable
    }
  }
  return false;
}

inline std::vector<size_t>
OCFEBASE::_evolution_slab_cols
( size_t iel_evo )
const
{
  std::vector<size_t> cols;
  if( !_evolution_dom_set || !_evolution_dom_var.dag() ) return cols;
  auto const evo_id = _evolution_dom_var.id();

  for( auto const& [var, var_dom] : _mVar ){
    // Only states that live on the evolution domain are slab-local.
    bool on_evo = false;
    for( auto const& v : var_dom ) if( v.id() == evo_id ){ on_evo = true; break; }
    if( !on_evo ) continue;

    auto itblk = _blkVar.find( var );
    if( itblk == _blkVar.end() ) continue;
    size_t const blk_tot = std::get<1>( itblk->second );

    // Enumerate the non-evolution element multi-index; fix the evolution
    // element to iel_evo.  pos_state gives the flat start of each block.
    std::vector<FFVar> other;
    std::vector<size_t> nel;
    for( auto const& v : var_dom ){
      if( v.id() == evo_id ) continue;
      other.push_back( v );
      nel.push_back( _mDom.at( v ).n_elem );
    }
    size_t total = 1;
    for( size_t e : nel ) total *= e;

    for( size_t lin = 0; lin < total; ++lin ){
      std::map<FFVar,size_t,lt_FFVar> ndx_el;
      ndx_el[ _evolution_dom_var ] = iel_evo;
      size_t rem = lin;
      for( size_t j = 0; j < other.size(); ++j ){
        ndx_el[ other[j] ] = ( nel[j] ? rem % nel[j] : 0 );
        if( nel[j] ) rem /= nel[j];
      }
      size_t const pos = pos_state( var, ndx_el );
      for( size_t k = 0; k < blk_tot; ++k ) cols.push_back( pos + k );
    }
  }
  std::sort( cols.begin(), cols.end() );
  return cols;
}

inline bool
OCFEBASE::set_input_values
( FFVar const& input, std::vector<double> const& values, double* inp )
const
{
  if( inp == nullptr ) return false;
  auto itv = _mInp.find( input );
  if( itv == _mInp.end() ) return false;
  auto const& var_dom = itv->second;
  auto itblk = _blkInp.find( input );
  if( itblk == _blkInp.end() ) return false;
  size_t const blk_tot = std::get<1>( itblk->second );

  // Element multi-index over the input's domains (var_dom order, first innermost),
  // blk_tot scalars per block -- the same element-major / node-inner ordering that
  // terminal_profile() and pos_input use, so a terminal profile re-seeds directly.
  std::vector<FFVar> doms;
  std::vector<size_t> nel;
  for( auto const& v : var_dom ){
    doms.push_back( v );
    nel.push_back( _input_dom( input, v ).n_elem );
  }
  size_t nblk = 1;
  for( size_t e : nel ) nblk *= e;
  if( nblk * blk_tot != values.size() ) return false;   // shape guard

  size_t vi = 0;
  for( size_t lin = 0; lin < nblk; ++lin ){
    std::map<FFVar,size_t,lt_FFVar> ndx_el;
    size_t rem = lin;
    for( size_t j = 0; j < doms.size(); ++j ){
      ndx_el[ doms[j] ] = ( nel[j] ? rem % nel[j] : 0 );
      if( nel[j] ) rem /= nel[j];
    }
    size_t const pos = pos_input( input, ndx_el );
    for( size_t k = 0; k < blk_tot; ++k ) inp[ pos + k ] = values[ vi++ ];
  }
  return true;
}

inline std::vector<double>
OCFEBASE::get_input_values
( FFVar const& input, double const* inp )
const
{
  std::vector<double> values;
  if( inp == nullptr ) return values;
  auto itv = _mInp.find( input );
  if( itv == _mInp.end() ) return values;
  auto const& var_dom = itv->second;
  auto itblk = _blkInp.find( input );
  if( itblk == _blkInp.end() ) return values;
  size_t const blk_tot = std::get<1>( itblk->second );

  std::vector<FFVar> doms;
  std::vector<size_t> nel;
  for( auto const& v : var_dom ){
    doms.push_back( v );
    nel.push_back( _input_dom( input, v ).n_elem );
  }
  size_t nblk = 1;
  for( size_t e : nel ) nblk *= e;
  values.reserve( nblk * blk_tot );

  for( size_t lin = 0; lin < nblk; ++lin ){
    std::map<FFVar,size_t,lt_FFVar> ndx_el;
    size_t rem = lin;
    for( size_t j = 0; j < doms.size(); ++j ){
      ndx_el[ doms[j] ] = ( nel[j] ? rem % nel[j] : 0 );
      if( nel[j] ) rem /= nel[j];
    }
    size_t const pos = pos_input( input, ndx_el );
    for( size_t k = 0; k < blk_tot; ++k ) values.push_back( inp[ pos + k ] );
  }
  return values;
}

template <typename FUN>
inline bool
OCFEBASE::sample_input
( FFVar const& input, FUN&& fun, double* inp )
const
{
  // Set distributed input @p input (in inp[]) by evaluating callable @p fun
  // (a t_Coord -> double) at the input's OWN collocation nodes, element-major /
  // node-inner in pos_input ordering.  General per-element seeding of marching
  // IC / boundary / forcing inputs from a physical-space (and/or time) function.
  if( inp == nullptr ) return false;
  auto itv = _mInp.find( input );
  if( itv == _mInp.end() ) return false;
  std::vector<FFVar> const doms( itv->second.begin(), itv->second.end() );

  std::vector<size_t> nel, nnode;
  std::vector<std::vector<std::vector<double>>> np;   // np[dom][elem][node]
  for( auto const& v : doms ){
    FFDom const& d = _input_dom( input, v );
    nel.push_back( d.n_elem );
    nnode.push_back( d.n_node );
    std::vector<std::vector<double>> per;
    for( size_t ie = 0; ie < d.n_elem; ++ie )
      per.push_back( d.lgnodes( d.elem_lo( ie ), d.elem_up( ie ) ) );
    np.push_back( std::move( per ) );
  }
  size_t nblk = 1; for( size_t e : nel )   nblk *= e;
  size_t blk  = 1; for( size_t n : nnode ) blk  *= n;

  std::vector<double> values;
  values.reserve( nblk * blk );
  for( size_t lin = 0; lin < nblk; ++lin ){
    std::vector<size_t> eidx( doms.size() );
    size_t rem = lin;
    for( size_t j = 0; j < doms.size(); ++j ){
      eidx[j] = ( nel[j] ? rem % nel[j] : 0 );
      if( nel[j] ) rem /= nel[j];
    }
    for( size_t nk = 0; nk < blk; ++nk ){
      t_Coord pt;
      size_t r = nk;
      std::vector<size_t> nidx( doms.size() );
      for( size_t j = 0; j < doms.size(); ++j ){
        size_t const ni = ( nnode[j] ? r % nnode[j] : 0 );
        nidx[j] = ni;
        pt[ doms[j] ] = np[j][ eidx[j] ][ ni ];
        if( nnode[j] ) r /= nnode[j];
      }
      // A marching transfer (transfer_terminal, 2026-10-06, WORKPLAN 1.2/1.3): a node on an INTERIOR element boundary of
      // a spatial direction is sampled from just inside ITS OWN element, so each copy of the node reads its own
      // element's terminal -- the per-element map terminal_profile() and the sensitivities use.  Sampled at the
      // coordinate itself, both copies read one element: invisible under IC_STRONG/IC_TRACE (equal copies), but under
      // IC_WEAK marching differed from monolithic and its gradient from its own finite differences.
      if( _sampleOwnElement ){
        for( size_t j = 0; j < doms.size(); ++j ){
          if( _evolution_dom_var.dag() && doms[j].id() == _evolution_dom_var.id() ) continue;
          if( nnode[j] < 2 || nel[j] < 2 ) continue;
          FFDom const& d = _input_dom( input, doms[j] );
          double const eps = 1e-12 * std::fabs( d.elem_up( eidx[j] ) - d.elem_lo( eidx[j] ) );
          if( nidx[j] == 0 && eidx[j] > 0 )                       pt[ doms[j] ] += eps;
          if( nidx[j] + 1 == nnode[j] && eidx[j] + 1 < nel[j] )   pt[ doms[j] ] -= eps;
        }
      }
      double v = fun( pt );
      // Interior-side re-sampling so a DISCONTINUOUS input gets each element's own one-sided value
      // (see _input_ref_interior_side).  No-op for continuous inputs.
      v = _input_ref_interior_side( fun, pt, v, input, doms, eidx, nidx, nnode );
      values.push_back( v );
    }
  }
  return set_input_values( input, values, inp );
}

inline size_t
OCFEBASE::_apply_c0_claim_filter
( std::map<t_InterfaceClaimKey,size_t>& claims )
const
{
  if( !_knob_CLAIM_FROM_C0() ) return 0;
  // EXACT MODES are INTERFACE.IMPOSITION (weak/strong/trace).  INTERFACE.TYPE is the value/flux/upwind/AUTO
  // resolver and has nothing to do with it -- gating on that one made the filter a silent no-op.
  if( options.INTERFACE.IMPOSITION != Options::IC_STRONG && options.INTERFACE.IMPOSITION != Options::IC_TRACE ) return 0;

  auto const cm = claim_classification();
  std::map<std::pair<size_t,size_t>,std::pair<size_t,ClaimTag>> dropped;
  size_t kept = 0;
  for( auto it = claims.begin(); it != claims.end(); ){
    FFVar const* sv = _state_by_id( it->first.state_id );
    FFVar const* dv = _domain_by_id( it->first.dom_id );
    ClaimTag tag = ClaimTag::NATURAL;
    bool known = false;
    if( sv && dv ){
      auto is = cm.find( *sv );
      if( is != cm.end() ){
        auto id = is->second.find( *dv );
        if( id != is->second.end() ){ tag = id->second.tag; known = true; }
      }
    }
    if( known && tag != ClaimTag::NATURAL ){
      auto& d = dropped[ { it->first.state_id, it->first.dom_id } ];
      ++d.first; d.second = tag;
      // rev221: keep the ROW (it is load-bearing); drop only the multiplier COLUMN, structurally.
      _strong_dropped_claims.insert( _freeze_interface_claim_key( it->first ) );
      ++kept; ++it;
    }
    else { ++kept; ++it; }
  }
  size_t ndrop = 0;
  for( auto const& kv : dropped ) ndrop += kv.second.first;
  // Printed even when NOTHING is suppressed: "ran and found nothing" must be distinguishable from "never ran".
  // This project has already been bitten once by a knob whose absence looked exactly like a no-op feature.
  _disp( 2 ) << "OCFESLV::setup ** claim[C0]: " << ndrop << " claim(s) column-dropped (rows kept), " << kept
            << " retained; classifier entries=" << cm.size() << " (exact mode, CRONOS_CLAIM_C0=1)\n";
  for( auto const& kv : dropped ){
    FFVar const* sv = _state_by_id( kv.first.first );
    FFVar const* dv = _domain_by_id( kv.first.second );
    _disp( 2 ) << "OCFESLV::setup ** claim[C0]:   suppressed " << kv.second.first << " on state ";
    if( sv ) _disp( 2 ) << *sv; else _disp( 2 ) << "#" << kv.first.first;
    _disp( 2 ) << " in ";
    if( dv ) _disp( 2 ) << *dv; else _disp( 2 ) << "#" << kv.first.second;
    _disp( 2 ) << "  (" << ( kv.second.second == ClaimTag::IMPLIED? "IMPLIED":
                            kv.second.second == ClaimTag::CLOSURE? "CLOSURE": "NONE" ) << ")\n";
  }
  return ndrop;
}

inline void
OCFEBASE::_apply_interface_edge_flags_v2
( InterfaceReceiverEdge& edge,
  FFModel::EqnOptions const&      eqnopt,
  bool                   aux_alg_receiver,
  bool                   is_geo_c0 ) const
{
  bool const non_physical_receiver = aux_alg_receiver || is_geo_c0;

  // Only INITIAL/BOUNDARY/DIAGNOSTIC are genuine boundary-condition rows.
  // INTERIOR/SURFACE/INTERFACE/LINK/AUTO are continuity-bearing and may carry
  // an exact value-C0 tau.
  bool role_is_bc = false;
  switch( eqnopt.role ){
    case EqnRole::INTERIOR:
    case EqnRole::SURFACE:
    case EqnRole::INTERFACE:
    case EqnRole::LINK:
    case EqnRole::AUTO:
      role_is_bc = false;
      break;
    case EqnRole::INITIAL:
    case EqnRole::BOUNDARY:
    case EqnRole::DIAGNOSTIC:
    default:
      role_is_bc = true;
      break;
  }
  bool const is_bc_row = role_is_bc && !non_physical_receiver;

  bool const has_coupling = ( std::fabs( edge.coupling ) > 64. * DBL_EPSILON );
  edge.hard_physical_row            = is_bc_row;
  edge.tangential_physical_boundary = false;
  {
    edge.weak_sat_allowed  = has_coupling;
    edge.trace_tau_allowed = has_coupling && !is_bc_row;
  }
}

inline void
OCFEBASE::_setup_eval_plans
()
{
  auto build_plan = [&]( std::vector< std::vector<t_DepUpdate> >& plan,
                         std::vector<FFDep> const& dep_rows )
  {
    plan.clear();
    plan.resize( dep_rows.size() );
    for( size_t k = 0; k < dep_rows.size(); ++k ){
      plan[k].reserve( dep_rows[k].dep().size() );
      for( auto const& [id, dep] : dep_rows[k].dep() ){
        (void)dep;
        auto it = _mDep.find( id );
        if( it == _mDep.end() ) continue;
        plan[k].push_back( { it->second.first, it->second.second } );
      }
    }
  };

  build_plan( _depPlanEqn, _depEqn );
  build_plan( _depPlanFct, _depFct );

  // section 4: per-equation evaluation domain -- the distributed directions each equation's
  // residual touches.  Built exactly like _fctEvalDom below, but over the equation subgraphs
  // _sgEqn (so it is inherently consistent with _depPlanEqn, same DAG source).  A domain-less
  // equation reduces every touched direction via FFIntegral; the residual assembler loops over
  // these "free" elements and accumulates (see _eval_eqn LOOP 2).
  _eqnEvalDom.clear();
  _eqnEvalDom.resize( _sgEqn.size() );
  for( size_t ne = 0; ne < _sgEqn.size(); ++ne ){
    auto& evaldom = _eqnEvalDom[ne];
    for( auto itop = _sgEqn[ne].l_op.cbegin(); itop != _sgEqn[ne].l_op.cend(); ++itop ){
      auto const& op = *itop;
      if( op->type != FFOp::VAR || _constant( *op->varout[0] ) ) continue;
      auto it_state = _mVar.find( *op->varout[0] );
      if( it_state != _mVar.end() ){
        evaldom.insert( it_state->second.cbegin(), it_state->second.cend() );
        continue;
      }
      auto it_input = _mInp.find( *op->varout[0] );
      if( it_input != _mInp.end() ){
        evaldom.insert( it_input->second.cbegin(), it_input->second.cend() );
        continue;
      }
    }
  }

  // section 4 (OpEval): collect the point reductions carried by FFEval nodes in each equation's
  // residual -- direction -> target coordinate.  An FFEval consumes its direction during the DAG
  // eval (the operator does the eval_lagrange collapse internally), so the direction is a POINT,
  // not a free integrated direction: it is excluded from free_dom below and its resident element
  // is pinned to contain the target coordinate.
  _eqnEvalPoint.assign( _sgEqn.size(), {} );
  _eqnEvalSide .assign( _sgEqn.size(), {} );
  for( size_t ne = 0; ne < _sgEqn.size(); ++ne ){
    for( auto const& op : _sgEqn[ne].l_op ){
      if( !op || !op->sameid( typeid( FFEval ) ) ) continue;
      auto const* fe = mc::type_cast<FFEval const>( op );
      if( !fe ) continue;
      for( auto const& [dv, z0] : fe->Coord() ){
        _eqnEvalPoint[ne][dv] = z0;
        _eqnEvalSide [ne][dv] = fe->side( dv );
      }
    }
  }

  // section 4: augment each domain-less equation's plan with DOMAIN updates for its free
  // (FFIntegral-reduced) directions.  build_plan derives updates from the residual's FFDep, but
  // OpI CONSUMES the integrated direction -- it drops out of the dependency set -- so no DOMAIN
  // entry is emitted and OCdom[dir] would never advance per free element in the LOOP-2 free
  // accumulate.  _mDep already carries the {DOMAIN,pos} entry for each direction; add any that
  // the plan is missing.  (Domain-less equations only; distributed equations already carry the
  // DOMAIN updates for their declared domain via the FFDep.)
  for( size_t ne = 0; ne < _eqnEvalDom.size(); ++ne ){
    if( ne >= _mEqn.size() || !_mEqn[ne].dom.empty() ) continue;
    for( auto const& dv : _eqnEvalDom[ne] ){
      auto itdep = _mDep.find( dv.id().second );
      if( itdep == _mDep.end() || itdep->second.first != DOMAIN ) continue;  // _mDep value: {typ,pos}
      bool present = false;
      for( auto const& upd : _depPlanEqn[ne] )
        if( upd.typ == DOMAIN && upd.pos == itdep->second.second ){ present = true; break; }
      if( !present )
        _depPlanEqn[ne].push_back( t_DepUpdate{ DOMAIN, itdep->second.second } );
    }
  }

  _fctEvalDom.clear();
  _fctEvalDom.resize( _mFct.size() );
  for( size_t nfct = 0; nfct < _mFct.size(); ++nfct ){
    auto& evaldom = _fctEvalDom[nfct];
    for( auto itop = _sgFct[nfct].l_op.cbegin(); itop != _sgFct[nfct].l_op.cend(); ++itop ){
      auto const& op = *itop;
      if( op->type != FFOp::VAR || _constant( *op->varout[0] ) ) continue;
      auto it_state = _mVar.find( *op->varout[0] );
      if( it_state != _mVar.end() ){
        evaldom.insert( it_state->second.cbegin(), it_state->second.cend() );
        continue;
      }
      auto it_input = _mInp.find( *op->varout[0] );
      if( it_input != _mInp.end() ){
        evaldom.insert( it_input->second.cbegin(), it_input->second.cend() );
        continue;
      }
    }
    for( auto const& [dvar, val] : _mFct[nfct].point ){
      (void)val;
      evaldom.insert( dvar );
    }
    for( auto const& [dvar, lim] : _mFct[nfct].grid ){
      (void)lim;
      evaldom.insert( dvar );
    }
  }
}

inline bool
OCFEBASE::_constant
( FFVar const& var )
const
{
  if( !var.dag() ) return false;
  for( auto const& cst : _vCst )
    if( var.id() == cst.id() ) return true;
  return false;
}

inline bool
OCFEBASE::_discretise_colloc_model
()
{
  // A deferred value's source is an EXPRESSION (FFModel::t_DeferredValue).  THIS solver evaluates it from the
  // collocated solution, so a source that is not already a state is materialised as one here -- with a defining row
  // and an auxiliary record, before the classification sees the system.  An integrator with a quadrature
  // right-hand side needs none of this, which is why the model no longer does it.
  for( auto& C : _deferredValue ){
    if( !C.source.dag() || _mVar.count( C.source ) ) continue;
    std::vector<FFVar> dom( C.source_dom.cbegin(), C.source_dom.cend() );
    std::string base = C.input.name();
    if( base.size() > 2 && base[base.size()-2] == '_' ) base.resize( base.size()-2 );
    FFVar w = _dag->add_var( base + "_R" );
    _add_working_state( w, dom );
    _register_aux_def( w, C.source, dom, std::set<FFVar,lt_FFVar>(), C.block_id );
    std::map<FFVar,int,lt_FFVar> ld;
    for( auto const& d : dom ) ld[d] = FFDom::ALL;
    _mEqn.push_back( { C.source - w, ld, _make_equation_options( _link_eqn_options( C.block_id, true ) ) } );
    C.source = w;
  }
  // A materialised source that reads a state jumping at a transition jumps with it (see _derive_transition_lifts)
  if( !_derive_transition_lifts() ){ _setupStatus = SetupStatus::INCONSISTENT_MODEL;  return false; }


  _issetup = false;
  _mDep.clear(); // dependency map for participating variables

  // Collocated state count, blocks and dependencies
  size_t ndx = 0;
  _nCollVar = 0; // total number of states after discretization
  _vVar.clear(); _vVar.reserve( _mVar.size() );
  _depVar.resize( _mVar.size() );
  _blkVar.clear();
  ndx = 0;
  for( auto const& [var,dom] : _mVar ){
    _vVar.push_back( var );

    // Set dependency
    _depVar[ndx].indep( var.id().second );
    _mDep[var.id().second] = {STATE, ndx};

    // Track collocated coefficient count and blocks
    size_t blk_size = 1;
    std::vector<size_t> blk_dim( dom.size() );
    size_t blk = 0;
    size_t cnt = 1;
    for( auto const& v : dom ){
      auto const& d = _mDom[v];
      blk_size *= d.n_node;
      blk_dim[blk++] = d.n_node;
      cnt *= d.n_elem * d.n_node;
    }
    _blkVar[var] = {_nCollVar, blk_size, std::move(blk_dim)};
    _nCollVar += cnt;
    
    ndx++;
  }

  // Collocated parameter count, blocks and dependencies
  _nCollInp = 0; // total number of inputs after discretization
  _vInp.clear(); _vInp.reserve( _mInp.size() );
  _depInp.resize( _mInp.size() );
  _blkInp.clear();
  ndx = 0;
  for( auto const& [var, dom] : _mInp ){
    _vInp.push_back( var );

    // Set dependency
    _depInp[ndx].indep( var.id().second );
    _mDep[var.id().second] = {INPUT, ndx};

    // Track collocated coefficient count and blocks
    size_t blk_size = 1;
    std::vector<size_t> blk_dim( dom.size() );
    size_t blk = 0;
    size_t cnt = 1;
    for( auto const& v : dom ){
      auto const& d = _input_dom( var, v );
      auto const& dphys = _mDom.at( v );
      if( _validate_input_grid_marching( v, d ) ){}   // a marching solver validated it against its march grid
      else if( d.n_elem != dphys.n_elem
       || std::fabs( d.lo_dom - dphys.lo_dom ) > 64.*DBL_EPSILON
       || std::fabs( d.up_dom - dphys.up_dom ) > 64.*DBL_EPSILON )
        throw Exceptions( Exceptions::INDEX );
      blk_size *= d.n_node;
      blk_dim[blk++] = d.n_node;
      cnt *= d.n_elem * d.n_node;
    }
    _blkInp[var] = {_nCollInp, blk_size, std::move(blk_dim)};
    _nCollInp += cnt;
    
    ndx++;
  }

  // Constants
  _depCst.resize( _vCst.size() );
  ndx = 0;
  for( auto const& var : _vCst ){
    // Set dependency
    _depCst[ndx].indep( var.id().second );
    _mDep[var.id().second] = {CONSTANT, ndx};
    
    ndx++;
  }

  // Setup domain collocation
  _vDom.clear(); _vDom.reserve( _mDom.size() );
  _depDom.resize( _mDom.size() );
  ndx = 0;
  for( auto& [var, dom] : _mDom ){
    _vDom.push_back( var );

    // Set dependency
    _depDom[ndx].indep( var.id().second );
    _mDep[var.id().second] = {DOMAIN, ndx};

    // Initialize collocation nodes
    dom.set_nodes();
#ifdef CRONOS__FFCOLLOC_DEBUG
    _disp( 3 ) << var << ": ";
    _disp( 3 ) << std::scientific << std::setprecision(6);
    for( auto const& x : dom.lgnodes() )
      _disp( 3 ) << x << "  ";
    _disp( 3 ) << std::endl;
#endif
    ndx++;
  }

  // Initialize optional input collocation nodes and verify matching element partitions
  for( auto& [inp, disc] : _mInpDisc ){
    auto itin = _mInp.find( inp );
    if( itin == _mInp.end() ) throw Exceptions( Exceptions::INDEX );
    for( auto& [dvar, idom] : disc ){
      if( itin->second.find( dvar ) == itin->second.end() )
        throw Exceptions( Exceptions::INDEX );

      { int const h = _init_input_nodes_marching( dvar, idom );   // a marching solver: against its march grid
        if( h < 0 ) return false;                                   // handled, node setup failed (status set)
        if( h > 0 ) continue; }                                     // handled: keep the full input partition

      auto const& pdom = _mDom.at( dvar );
      if( idom.n_elem != pdom.n_elem
       || idom.elem_bnd.size() != pdom.elem_bnd.size()
       || idom.elem_len.size() != pdom.elem_len.size()
       || std::fabs( idom.lo_dom - pdom.lo_dom ) > 64.*DBL_EPSILON
       || std::fabs( idom.up_dom - pdom.up_dom ) > 64.*DBL_EPSILON )
        throw Exceptions( Exceptions::INDEX );
      for( size_t ie=0; ie<pdom.elem_bnd.size(); ++ie )
        if( std::fabs( idom.elem_bnd[ie] - pdom.elem_bnd[ie] ) > 64.*DBL_EPSILON )
          throw Exceptions( Exceptions::INDEX );
      idom.w_elem   = pdom.w_elem;
      idom.elem_bnd = pdom.elem_bnd;
      idom.elem_len = pdom.elem_len;
      if( !idom.set_nodes() ){ _setupStatus = SetupStatus::NODE_SETUP_FAILED; return false; }
    }
  }

  // Cache physical nodes for every domain element.  Must precede _classify_pde():
  // the classify-phase mesh-wavenumber probe (k_h = pi/dx_min) reads
  // _nodeBlockCache, so it has to be populated before classification runs.
  _setup_node_cache();

  return true;
}

inline bool
OCFEBASE::_build_colloc_plan
()
{
  size_t ndx = 0;

  // Reuse the same consistency helper on the local, possibly order-reduced
  // model, and populate setup's equation/output row counts and subgraph caches.
  if( !_check_consistency( _dag, _vCst, &_vCstVal, _mDom, _mVar, _mInp,
                           _mEqn, _mFct, &_nCollEqn, &_sgEqn,
                           &_nCollFct, &_sgFct ) )
    { _setupStatus = SetupStatus::INCONSISTENT_MODEL; return false; }

  _phase( "consistency + counts" );

  // Map dependencies
  _depEqn.resize( _mEqn.size() );
  ndx = 0;
  for( auto const& eqn : _mEqn ){
    _dag->eval( _sgEqn[ndx], 1, &eqn.var, &_depEqn[ndx], _vVar.size(), _vVar.data(), _depVar.data(),
                _vInp.size(), _vInp.data(), _depInp.data(), _vCst.size(), _vCst.data(), _depCst.data(),
                _vDom.size(), _vDom.data(), _depDom.data() );
    ndx++;
  }

#ifdef MC__OCFESLV_DEPEQN_PROBE
  // [depeqn] per-equation state-dependency dump (OPT-IN, compile flag).  For each equation, list its
  // domain variables and the STATES that _depEqn recorded it depends on.  Localizes cross-domain
  // coefficient-state registration (e.g. a {z}-state used in the coefficients of {rwc,z} equations):
  // if such a state is absent here, the drop is in the dependency eval; if present, look downstream.
  {
    size_t pndx = 0;
    for( auto const& eqnp : _mEqn ){
      _disp( 3 ) << "  [depeqn] eqn#" << pndx << " dom={";
      for( auto const& [dv,b] : eqnp.dom ){ (void)b; _disp( 3 ) << " " << dv; }
      _disp( 3 ) << " } states:";
      for( auto const& [id,d] : _depEqn[pndx].dep() ){
        (void)d;
        auto it = _mDep.find( id );
        if( it != _mDep.end() && it->second.first == STATE )
          _disp( 3 ) << " " << _vVar[ it->second.second ];
      }
      _disp( 3 ) << "\n";
      ++pndx;
    }
  }
#endif // MC__OCFESLV_DEPEQN_PROBE

  // Map output dependencies
  _depFct.resize( _mFct.size() );
  ndx = 0;
  for( auto const& fct : _mFct ){
    _dag->eval( _sgFct[ndx], 1, &fct.var, &_depFct[ndx], _vVar.size(), _vVar.data(), _depVar.data(),
                _vInp.size(), _vInp.data(), _depInp.data(), _vCst.size(), _vCst.data(), _depCst.data(),
                _vDom.size(), _vDom.data(), _depDom.data() );
    ndx++;
  }

  // Build compact dependency-update plans used during repeated eval() calls.
  _setup_eval_plans();
  _phase( "dependency_map + eval_plans" );

  // Build the principal-symbol SAT coupling table as a setup-time input to
  // both IC_WEAK and IC_TRACE.  From this point onward, evaluation must only
  // consume the resulting coupling/plan data, not rebuild interface topology.
  _ensure_sat_coupling();
  _phase( "sat_coupling" );

  // IC_TRACE appends one multiplier unknown per exact trace/interface row and
  // keeps all physical residual rows.  The complete interface plan is built
  // transactionally into a local draft, validated, and committed only once.
  // Thus a setup failure cannot leave stale or half-built receiver/claim/tau
  // structures behind, and eval()/deriv() can only consume a sealed plan.
  InterfacePlanDraft interface_plan_draft;
  if( !_build_interface_plan_draft( interface_plan_draft )
   || !_validate_interface_plan_draft( interface_plan_draft ) ){
    _issetup = false;
    { _setupStatus = SetupStatus::INTERFACE_PLAN_INVALID; return false; }
  }
  _commit_interface_plan_draft( std::move( interface_plan_draft ) );
  _phase( "interface_plan (build+validate+commit)" );

  // Reset collocation variables
  _reset_colloc( _OCVar, _OCInp, _OCDom, _OCCst );
  _reset_colloc( _OCDepVar, _OCDepInp, _OCDepDom, _OCDepCst );
  _reset_colloc( _OCFVar, _OCFInp, _OCFDom, _OCFCst );
  _OCFEqn = OCVar<t_FAD>();
  _wrkOCF.clear();
  _phase( "reset_colloc" );

  // ---------------------------------------------------------------------------
  // rev107 [taucol/plan] -- direct named-column instrument, plan half.
  // ---------------------------------------------------------------------------
  // Placed HERE, between the plan commit and the FFDep structural pass below, and the
  // position is a decision rather than a convenience:
  //
  //   * it must be AFTER _commit_interface_plan_draft -- it reads the committed plan;
  //   * putting it BEFORE the FFDep pass makes THAT pass the assembly observation.
  //     The pass runs _eval_eqn end to end and every quantity the assembly half records
  //     is INDEX-BASED and value-independent (row keys, pivot tags, correction indices),
  //     so FFDep gives the same answer a double run gives.  [stage2-preflight] already
  //     reports from this pass for exactly that reason.
  //   * and therefore the instrument produces its rows under CRONOS_AUDIT_SETUP_ONLY,
  //     which returns false further down without ever reaching a real eval().  Placed
  //     after PLAN METRIC instead, the assembly half would print nothing in a corpus
  //     sweep -- the configuration it most needs to cover.
  //
  // CONSEQUENCE FOR READING A LOG: [taucol/*] now precedes [evoface] PLAN METRIC, so k
  // arrives AFTER the trace that is supposed to account for it.  Scroll down for k; do
  // not read the trace as self-contained.
  //
  // Deliberately NOT folded into the kEvoFace block: that flag already gates several
  // reports and a reader cannot tell which of them a line came from.  A diagnostic
  // sharing a gate with another diagnostic is one deletion away from disappearing
  // silently -- n_elem@plan shared a std::cerr with the predicate-(ii) report and went
  // with it in rev101.
  if( _taucol_enabled() ) _report_tau_column_plan( _disp( 3 ) );

  // Structural self-check (Part 2): a dry FFDep residual pass verifies that
  // Frozen exact interface rows and surviving physical equation nodes emitted
  // by _eval_eqn() jointly produce exactly _nCollEqn rows for the active
  // interface and imposition options.  Index-based and value-independent; surfaces a
  // squareness violation here rather than at the first eval()/deriv() call.
  if( !_audit_strong_interface_rows() )
    { _setupStatus = SetupStatus::AUDIT_NONSQUARE; return false; }
  _phase( "audit_strong_interface_rows (FFDep pass)" );
  return true;
}

inline void
OCFEBASE::_reset_colloc_state
()
{
  _mDom.clear();
  _mInpDisc.clear();
  _blockSymbol.clear();
  _blockClassification.clear();
  _blockFaceData.clear();
  _satCoupling.clear();
  _satCouplingBuilt      = false;
  _satCouplingClassified = false;
  _satCouplingSerial     = 0;
  _auxDef.clear();
  _auxVarID.clear();
  _eqnAuxDirCache.clear();
  _evolution_dom_user = false;
  _evolution_dom_set  = false;
  _depPlanEqn.clear();
  _depPlanFct.clear();
  _fctEvalDom.clear();
  _nodeBlockCache.clear();
  _OCFVar.clear();
  _OCFInp.clear();
  _OCFDom.clear();
  _OCFCst.clear();
  _wrkOCF.clear();
  _classificationSerial = 0;
  _nCollEqn = _nCollVar = _nCollInp = _nCollFct = 0;
  _nTraceVar = _traceVarOffset = _nTraceConstraint = 0;
  _plan.reset();
  _interfacePlanBuildAllowed = false;
  _evoDerivRows.clear();
}

inline std::vector<unsigned>
OCFEBASE::_color_jacobian_columns
( std::vector< std::vector<size_t> > const& pattern,
  size_t const ncol, unsigned& ncolor )
const
{
  // Column-intersection coloring: columns that appear in a common residual row
  // must receive distinct colors.  We do not need to materialize the full
  // column-intersection graph.  Instead, keep the row incidence of each column
  // and, during greedy coloring, scan only the rows containing the current
  // column to find already colored conflicts.
  std::vector< std::vector<size_t> > rows_of_col( ncol );
  std::vector<size_t> degree_proxy( ncol, 0 );

  for( size_t k = 0; k < pattern.size(); ++k ){
    auto const& row = pattern[k];
    size_t const rowdeg = row.empty()? 0: row.size() - 1;
    for( size_t col : row ){
      if( col >= ncol ) throw Exceptions( Exceptions::INTERNAL );
      rows_of_col[col].push_back( k );
      degree_proxy[col] += rowdeg;
    }
  }

  std::vector<size_t> order( ncol );
  std::iota( order.begin(), order.end(), 0 );
  std::sort( order.begin(), order.end(),
             [&]( size_t a, size_t b ){
               if( degree_proxy[a] != degree_proxy[b] )
                 return degree_proxy[a] > degree_proxy[b];
               return a < b;
             } );

  std::vector<unsigned> color( ncol, std::numeric_limits<unsigned>::max() );
  std::vector<unsigned char> used;
  ncolor = 0;

  for( size_t col : order ){
    if( used.size() < ncolor ) used.resize( ncolor, 0 );
    std::vector<unsigned> touched;

    for( size_t k : rows_of_col[col] ){
      for( size_t nb : pattern[k] ){
        if( nb == col ) continue;
        unsigned const cnb = color[nb];
        if( cnb == std::numeric_limits<unsigned>::max() ) continue;
        if( cnb >= used.size() ) used.resize( cnb+1, 0 );
        if( !used[cnb] ){
          used[cnb] = 1;
          touched.push_back( cnb );
        }
      }
    }

    unsigned c = 0;
    while( c < used.size() && used[c] ) ++c;
    color[col] = c;
    if( c + 1 > ncolor ){
      ncolor = c + 1;
      if( used.size() < ncolor ) used.resize( ncolor, 0 );
    }

    for( unsigned ct : touched ) used[ct] = 0;
  }

#ifdef CRONOS__FFCOLLOC_CHECK
  std::vector<unsigned char> seen( ncolor, 0 );
  for( auto const& row : pattern ){
    std::vector<unsigned> touched;
    touched.reserve( row.size() );
    for( size_t col : row ){
      unsigned const c = color[col];
      if( c >= ncolor ) throw Exceptions( Exceptions::INTERNAL );
      if( seen[c] ) throw Exceptions( Exceptions::INTERNAL );
      seen[c] = 1;
      touched.push_back( c );
    }
    for( unsigned c : touched ) seen[c] = 0;
  }
#endif

  return color;
}

inline void
OCFEBASE::_deriv_copy_pattern
( std::vector<std::vector<size_t>> const& pattern, size_t* nnz, size_t** colnz )
{
  for( size_t k = 0; k < pattern.size(); ++k ){
    auto const& row = pattern[k];
    if( nnz ) nnz[k] = row.size();
    if( colnz ){
      if( !row.empty() && !colnz[k] )
        throw Exceptions( Exceptions::INTERNAL );
      for( size_t j = 0; j < row.size(); ++j )
        colnz[k][j] = row[j];
    }
  }
}

template <typename U>
inline void
OCFEBASE::_deriv_extract_rows
( FADType<U> const* val, std::vector< std::vector<size_t> > const& pattern,
  std::vector<unsigned> const& color, size_t const ncol, U* res, U** grad )
{
  for( size_t k = 0; k < pattern.size(); ++k ){
    if( res ) res[k] = val[k].val();

    if( !grad ) continue;
    auto const& row = pattern[k];
    if( !row.empty() && !grad[k] )
      throw Exceptions( Exceptions::INTERNAL );

    for( size_t j = 0; j < row.size(); ++j ){
      size_t const col = row[j];
      if( col >= ncol )
        throw Exceptions( Exceptions::INTERNAL );
      unsigned const c = color[col];
      grad[k][j] = val[k].deriv( c );
    }
  }
}

template <typename U>
inline void
OCFEBASE::_deriv_eqn
( U* eqn, U** gradeqn, size_t* nnzeqn, size_t** colnzeqn,
  OCVar<FADType<U>>& OCFeqn, std::vector<OCVar<FADType<U>>>& OCFvar,
  std::vector<OCVar<FADType<U>>>& OCFinp, std::vector<OCVar<FADType<U>>>& OCFdom,
  std::vector<OCVar<FADType<U>>>& OCFcst, std::vector<OCVar<FADType<U>>>& OCFwrk,
  FADType<U>* Feqn, FADType<U>* Fvar, FADType<U>* Finp )
const
{
  if( !_jacReady
   || _jacInterfaceType  != options.INTERFACE.TYPE
   || _jacImpositionType != options.INTERFACE.IMPOSITION
   || _jacClassified     != _classified )
    throw Exceptions( Exceptions::SETUP );

  size_t const ncol = _nCollVar + _nCollInp;
  if( gradeqn && _jacColor.size() != ncol )
    throw Exceptions( Exceptions::INTERNAL );

  // Populate the sparsity pattern. FFDep pass and coloring were completed by setup().
  _deriv_copy_pattern( _jacPatternEqn, nnzeqn, colnzeqn );

  if( !eqn && !gradeqn )
    return;

  // This is the only residual propagation performed in _deriv_eqn().
  _eval_eqn( Feqn, Fvar, Finp, OCFeqn, OCFvar, OCFinp, OCFdom, OCFcst, OCFwrk );

  // Populate the value and gradient arrays.
  _deriv_extract_rows( Feqn, _jacPatternEqn, _jacColor, ncol, eqn, gradeqn );
}

inline bool
OCFEBASE::_build_jacobian_pattern
()
{
  _jacReady = false;
  _jacPatternEqn.clear();
  _jacPatternFct.clear();
  _jacLinearEqn.clear();
  _jacLinearFct.clear();
  _jacColor.clear();
  _jacNColor = 0;
  _jacAllEqnLinear = false;
  _jacAllFctLinear = false;

  try{
    size_t const ncol = _nCollVar + _nCollInp;
    if( ncol > static_cast<size_t>( std::numeric_limits<int>::max() ) )
      throw Exceptions( Exceptions::INTERNAL );

    std::vector<FFDep> var_dep( _nCollVar );
    std::vector<FFDep> inp_dep( _nCollInp );
    for( size_t i = 0; i < _nCollVar; ++i )
      var_dep[i].indep( static_cast<int>( i ) );
    for( size_t i = 0; i < _nCollInp; ++i )
      inp_dep[i].indep( static_cast<int>( _nCollVar + i ) );

    std::vector<FFDep> eqn_dep( _nCollEqn );
    std::vector<FFDep> fct_dep( _nCollFct );

    _reset_colloc( _OCDepVar, _OCDepInp, _OCDepDom, _OCDepCst );
    if( _nCollEqn )
      _eval_eqn( eqn_dep.data(), var_dep.empty()? nullptr: var_dep.data(),
                 inp_dep.empty()? nullptr: inp_dep.data(),
                 _OCDepEqn, _OCDepVar, _OCDepInp, _OCDepDom, _OCDepCst, _wrkOCDep );
    if( _nCollFct )
      _eval_fct_dep( fct_dep, var_dep, inp_dep );   // the solver's output functions, if it has any

    _jacPatternEqn.assign( _nCollEqn, std::vector<size_t>() );
    _jacPatternFct.assign( _nCollFct, std::vector<size_t>() );
    _jacLinearEqn.assign( _nCollEqn, true );
    _jacLinearFct.assign( _nCollFct, true );

    auto fill = [&]( std::vector<FFDep> const& dep_rows,
                     std::vector<std::vector<size_t>>& pat,
                     std::vector<bool>& lin )
    {
      for( size_t k = 0; k < dep_rows.size(); ++k ){
        auto const& dep = dep_rows[k].dep();
        pat[k].reserve( dep.size() );
        bool islin = true;
        for( auto const& kv : dep ){
          if( kv.first >= 0 && static_cast<size_t>( kv.first ) < ncol )
            pat[k].push_back( static_cast<size_t>( kv.first ) );
          if( kv.second != FFDep::L ) islin = false;
        }
        // FFDep::dep() is ordered by column, but keep this robust if FFDep's
        // internal container changes in the future.
        std::sort( pat[k].begin(), pat[k].end() );
        pat[k].erase( std::unique( pat[k].begin(), pat[k].end() ), pat[k].end() );
        lin[k] = islin;
      }
    };
    fill( eqn_dep, _jacPatternEqn, _jacLinearEqn );
    fill( fct_dep, _jacPatternFct, _jacLinearFct );

    if( _jacPatternEqn.size() != _nCollEqn || _jacPatternFct.size() != _nCollFct )
      throw Exceptions( Exceptions::INTERNAL );

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
    return true;
  }
  catch(...){
    _jacPatternEqn.clear();
    _jacPatternFct.clear();
    _jacLinearEqn.clear();
    _jacLinearFct.clear();
    _jacAllEqnLinear = false;
    _jacAllFctLinear = false;
    _jacColor.clear();
    _jacNColor = 0;
    _jacReady = false;
    return false;
  }
}

template <typename U>
inline void
OCFEBASE::_eval_fct
( U* fct, U const* var, U const* inp,
  OCVar<U>& OCres, std::vector<OCVar<U>>& OCvar,  std::vector<OCVar<U>>& OCinp,
  std::vector<OCVar<U>>& OCdom, std::vector<OCVar<U>>& OCcst, std::vector<OCVar<U>>& OCwrk )
const
{
  auto const ew = _evolution_window();               // what of the evolution direction is collocated (the solver answers)
  // Validate external coefficient arrays before any output assembly.
  if( _nCollVar && !var )
    throw Exceptions( Exceptions::INTERNAL );
  if( _nCollInp && !inp )
    throw Exceptions( Exceptions::INTERNAL );
  if( _nCollFct && !fct )
    throw Exceptions( Exceptions::INTERNAL );

  // =========================================================================
  // LOOP 3: POST-EVALUATION OUTPUT FUNCTIONS
  //
  // Point outputs contribute one scalar row.  Distributed outputs contribute
  // one row per collocation node selected by their grid mask.  These rows are
  // appended to the fct[] array only; they never participate in equation row
  // counts, classification, or interface residual/SAT assembly.
  // =========================================================================
  size_t nfct = 0;
  U* fct_base = fct;

  for( auto const& fctrec : _mFct ){
    auto const& fctvar   = fctrec.var;
    auto const& fctpoint = fctrec.point;
    auto const& fctgrid  = fctrec.grid;

#ifdef CRONOS_DEBUG_OUT
    _disp( 3 ) << "eval output " << nfct << " begin\n";
#endif

    std::set<FFVar,lt_FFVar> evaldom =
      nfct < _fctEvalDom.size()? _fctEvalDom[nfct] : std::set<FFVar,lt_FFVar>();

    std::map<FFVar, size_t, lt_FFVar> fixed_el;
    std::map<FFVar, double, lt_FFVar> ref_el;
    for( auto const& [dvar, val] : fctpoint ){
      auto const& dom = _mDom.at( dvar );
      auto const isd = fctrec.side.find( dvar );
      auto const [iel, xi] = _point_in_domain( dom, val, isd == fctrec.side.end()? (int)FFDom::MINUS: isd->second );
      fixed_el[dvar] = iel;
      ref_el[dvar]   = xi;
    }

    // Element dimensions not fixed by a point and not part of the distributed
    // grid are looped over and accumulated.  This preserves the existing scalar
    // output behaviour for FFIntegral terms, where integrated dimensions no
    // longer appear in OCres.dom() after evaluation.
    // A direction consumed by an FFEval INSIDE the output expression is PINNED to the element holding its point (with
    // the evaluation's side), exactly as the equation path pins _eqnEvalPoint -- not accumulated.  Without this, under
    // REDUCE.ORDER = RED_NONE (no capture), OpE(x,t,tau) was summed over every element of t: n_elem times its value.
    if( fctrec.var.dag() ){
      FFSubgraph sgf = _dag->subgraph( 1, &fctrec.var );
      for( auto const* op : sgf.l_op ){
        auto const* fe = op? mc::type_cast<FFEval const>( op ): nullptr;
        if( !fe ) continue;
        for( auto const& [dv, z0] : fe->Coord() ){
          if( fixed_el.count( dv ) || !_mDom.count( dv ) ) continue;
          auto const [iel, xi] = _point_in_domain( _mDom.at( dv ), z0, fe->side( dv ) );  (void)xi;
          fixed_el[dv] = iel;
        }
      }
    }
    std::set<FFVar,lt_FFVar> free_dom;
    for( auto const& dvar : evaldom ){
      if( fctpoint.find( dvar ) != fctpoint.end() ) continue;
      if( fixed_el.count( dvar ) ) continue;                  // pinned by an FFEval in the expression (above)
      if( fctgrid.find( dvar )  != fctgrid.end()  ) continue;
      free_dom.insert( dvar );
    }

    auto transform_buffer = [&]( std::vector<U> buf, std::map<FFVar,size_t,lt_FFVar> const& ndx_el ) -> std::vector<U>
    {
      struct FDim {
        FFDom const* pdom;
        FFVar const* pvar;
        size_t       cur_n;
      };

      std::vector<FDim> fdims;
      fdims.reserve( OCres.dom().size() );
      for( auto const& [pvar, pdom] : OCres.dom() )
        fdims.push_back({ pdom, pvar, pdom->n_node });

      for( size_t idim = 0; idim < fdims.size(); ++idim ){
        auto& fd = fdims[idim];

        size_t stride = 1;
        for( size_t d = 0; d < idim; ++d ) stride *= fdims[d].cur_n;
        size_t outer = 1;
        for( size_t d = idim+1; d < fdims.size(); ++d ) outer *= fdims[d].cur_n;

        auto itref = ref_el.find( *fd.pvar );
        if( itref != ref_el.end() ){
          std::vector<U> newbuf;
          if( !fd.pdom->eval_lagrange( newbuf, itref->second, 0, buf, stride, outer ) )
            throw Exceptions( Exceptions::INTERNAL );
          buf      = std::move( newbuf );
          fd.cur_n = 1;
          continue;
        }

        auto itgrid = fctgrid.find( *fd.pvar );
        if( itgrid == fctgrid.end() )
          throw Exceptions( Exceptions::DOMAIN );

        int lim = itgrid->second;
        size_t const iel = ndx_el.count(*fd.pvar)? ndx_el.at(*fd.pvar) : size_t(0);
        size_t const n   = fd.cur_n;
        size_t const nel = fd.pdom->n_elem;
        // marching (2026-10-06, WORKPLAN 1.1): the evolution direction is collapsed to the resident window, ONE element
        // that would count as both the first and the last.  The window emits ALL its evolution nodes (no mask drop,
        // no unique-node drop) -- the same rows, hence the same Jacobian pattern, in every window -- and
        // _window_output_rows applies the GLOBAL mask when it places them in the monolithic layout.
        bool const evo_all = _evolution_dom_var.dag() && fd.pvar->id() == _evolution_dom_var.id() && ew.n_elements >= 1
                          && lim != FFDom::LB && lim != FFDom::UB;
        if( evo_all ) lim = FFDom::ALL;

        if( lim == FFDom::LB || lim == FFDom::UB ){
          double const t_eval = ( lim == FFDom::LB ) ? -1.0 : 1.0;
          std::vector<U> newbuf;
          if( !fd.pdom->eval_lagrange( newbuf, t_eval, 0, buf, stride, outer ) )
            throw Exceptions( Exceptions::INTERNAL );
          buf      = std::move( newbuf );
          fd.cur_n = 1;
        }
        else{
          std::vector<size_t> vi;
          vi.reserve( n );
          for( size_t k = 0; k < n; ++k ){
            bool ok = true;
            switch( lim ){
              case FFDom::ALL:              break;
              case FFDom::ALL-FFDom::LB:
                if( iel == 0      && k == 0   ) ok = false;
                break;
              case FFDom::ALL-FFDom::UB:
                if( iel == nel-1  && k == n-1 ) ok = false;
                break;
              case FFDom::ALL-FFDom::LB-FFDom::UB:
                if( iel == 0      && k == 0   ) ok = false;
                if( iel == nel-1  && k == n-1 ) ok = false;
                break;
              default: throw Exceptions( Exceptions::INTERNAL );
            }

            if( ok && !evo_all
             && options.OUTPUT.GRID_POLICY == Options::OUT_UNIQUE_PHYSICAL
             && ( fd.pdom->type == FFDom::LGL || fd.pdom->type == FFDom::CGL )
             && iel > 0 && k == 0 )
              ok = false;

            if( ok ) vi.push_back( k );
          }

          if( vi.size() < n ){
            std::vector<U> newbuf;
            newbuf.reserve( stride * vi.size() * outer );
            for( size_t o = 0; o < outer; ++o )
              for( size_t vk : vi )
                for( size_t s = 0; s < stride; ++s )
                  newbuf.push_back( buf[ s + vk * stride + o * n * stride ] );
            buf      = std::move( newbuf );
            fd.cur_n = vi.size();
          }
        }
      }

      return buf;
    };

    auto eval_block = [&]( std::map<FFVar,size_t,lt_FFVar> const& ndx_el ) -> std::vector<U>
    {
      auto const& plan_fct = ( nfct < _depPlanFct.size() ? _depPlanFct[nfct]
                                                         : std::vector<t_DepUpdate>() );
      for( auto const& upd : plan_fct ){
        size_t inc = 0;
        switch( upd.typ ){
          case DOMAIN:  OCdom[upd.pos].update( _node_block( _vDom[upd.pos], ndx_el ), ndx_el ); break;
          case STATE:   inc = pos_state( _vVar[upd.pos], ndx_el );
                        OCvar[upd.pos].update( &var[inc], ndx_el ); break;
          case INPUT: { auto ndx_in = ndx_el;   // marching: remap evolution element to the resident window
                        if( ew.kind == EvolutionWindow::ONE && _evolution_dom_var.dag() && _mInpDisc.count( _vInp[upd.pos] ) ){
                          auto ite = ndx_in.find( _evolution_dom_var );
                          if( ite != ndx_in.end() ) ite->second = ew.index; }
                        inc = pos_input( _vInp[upd.pos], ndx_in );
                        OCinp[upd.pos].update( &inp[inc], ndx_in ); break; }
          case CONSTANT: break;
          default: throw Exceptions( Exceptions::INTERNAL );
        }
      }

      _dag->eval( _sgFct[nfct], OCwrk,
                  1, &fctvar, &OCres,
                  _vVar.size(), _vVar.data(), OCvar.data(),
                  _vInp.size(), _vInp.data(), OCinp.data(),
                  _vCst.size(), _vCst.data(), OCcst.data(),
                  _vDom.size(), _vDom.data(), OCdom.data() );

      return transform_buffer( std::vector<U>( OCres.coef().begin(), OCres.coef().end() ), ndx_el );
    };

    size_t out_row = fctrec.row0;
    std::map<FFVar, size_t, lt_FFVar> grid_el;
    bool grid_done = false;
    bool grid_once = fctgrid.empty();

    while( true ){
      if( fctgrid.empty() ){
        if( grid_done ) break;
        grid_done = true;
      }
      else{
        if( !_next_block( grid_el, fctgrid ) || grid_el.empty() ) break;
      }

      std::vector<U> acc;
      bool have_acc = false;

      std::map<FFVar, size_t, lt_FFVar> free_el;
      bool free_done = false;
      while( true ){
        if( free_dom.empty() ){
          if( free_done ) break;
          free_done = true;
        }
        else{
          if( !_next_block( free_el, free_dom ) || free_el.empty() ) break;
        }

        std::map<FFVar, size_t, lt_FFVar> ndx_el( fixed_el );
        for( auto const& [dvar, iel] : grid_el ) ndx_el[dvar] = iel;
        for( auto const& [dvar, iel] : free_el ) ndx_el[dvar] = iel;

        std::vector<U> vals = eval_block( ndx_el );
        if( !have_acc ){
          acc = std::move( vals );
          have_acc = true;
        }
        else{
          if( acc.size() != vals.size() )
            throw Exceptions( Exceptions::INTERNAL );
          for( size_t i = 0; i < acc.size(); ++i ) acc[i] += vals[i];
        }
      }

      if( !have_acc )
        throw Exceptions( Exceptions::INTERNAL );

      if( fctrec.kind == FctKind::POINT ){
        if( acc.size() != 1 )
          throw Exceptions( Exceptions::INTERNAL );
        fct[out_row++] = acc[0];
      }
      else{
        for( U const& val : acc )
          fct[out_row++] = val;
      }

      if( grid_once ) break;
    }

    if( ew.n_elements >= 1 && _output_over_evolution( fctrec ) ){
      if( out_row != fctrec.row0 + _window_output_rows( fctrec, ew.kind == EvolutionWindow::ONE? ew.index: 0 ).size() ) throw Exceptions( Exceptions::INTERNAL );
    }
    else if( out_row != fctrec.row0 + fctrec.nrow )
      throw Exceptions( Exceptions::INTERNAL );

#ifdef CRONOS_DEBUG_OUT
    _disp( 3 ) << "eval output " << nfct << " assigned rows "
              << fctrec.row0 << ":" << out_row << "\n";
#endif
    ++nfct;
  }

  (void)fct_base;
}

} // namespace mc

#if defined(_WIN32)
# pragma pop_macro("INTERFACE")
#endif
#endif
