// Copyright (C) Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the EPL-2.0 with GPL-2.0-or-later as a Secondary License; see the LICENSE file.
// SPDX-License-Identifier: EPL-2.0 OR GPL-2.0-or-later

#ifndef CRONOS__OCPLAN_HPP
#define CRONOS__OCPLAN_HPP

// <windows.h>, in its FULL form (not WIN32_LEAN_AND_MEAN), leaves INTERFACE defined as a macro; this header
// declares OCFESLV::Options::INTERFACE.  Set aside for the header and restored after it (2026-10-09).
#if defined(_WIN32)
# pragma push_macro("INTERFACE")
# undef INTERFACE
#endif

#include <cstddef>
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <tuple>
#include <unordered_map>
#include <vector>
#include <iostream>
#include <cstdint>
#include <string>
#include <utility>
#include <algorithm>
#include <functional>
#include <armadillo>


namespace mc
{

// Forward-declare the interface-condition type enum so OCPlan structs can
// reference it without pulling in the full OCFESLV/OCBase headers.  The
// full definition lives in OCBase; this forward reference is resolved
// when ocplan.hpp is included from ocfeslv.hpp after ocbase.hpp.
class OCBase;


// ======================================================================
// OCPlan — frozen interface-condition plan for OCFESLV
// ======================================================================
//
// OCPlan owns all data that describes a validated, sealed interface plan:
//
//   * InterfaceReceiverEdge   — setup-time receiver-edge topology (shared
//                               by IC_WEAK, IC_TRACE and IC_STRONG).
//   * FrozenWeakSatTerm       — executable SAT/tau-injection consumer rows.
//   * FrozenTraceConstraintTerm — executable exact-trace rows (IC_TRACE,
//                               IC_STRONG).
//   * FrozenStrongTauElimination — Schur projection plan (IC_STRONG).
//   * Plan dimensions induced by trace variables (n_coll_eqn, n_coll_var,
//                               n_trace_var, n_trace_constraint, ...).
//   * Plan-state protocol (EMPTY / BUILDING / FROZEN).
//   * Validation logic.
//   * Key-freeze / key-unfreeze helpers.
//
// OCPlan does NOT own:
//   * The symbolic model (DAG, equations, domains, states).
//   * Classification results.
//   * OCFESLV options.
//
// The builder methods that populate the plan (topology discovery,
// _prepare_interface_plan_tables, etc.) remain in OCFESLV because they
// depend heavily on OCFESLV's symbolic model and classification internals.
// OCFESLV::setup() builds a transactional InterfacePlanDraft, validates it
// via OCPlan::validate(), and commits it into the live OCPlan member via
// OCPlan::commit_draft().
//
// _eval_eqn() accesses the sealed plan exclusively through OCPlan's
// const query interface, so it never reaches into raw OCFESLV members for
// plan data.
// ======================================================================


// ======================================================================
// OCPlan
// ======================================================================
//! @brief Frozen interface-condition plan produced by OCFESLV::setup().
//!
//! OCPlan holds all data consumed by OCFESLV::_eval_eqn() and
//! OCFESLV::_eval_fct() for interface/continuity enforcement.  It is
//! populated transactionally during setup via a local InterfacePlanDraft,
//! validated, and committed as a sealed unit.  After commit the plan is
//! immutable until the next call to OCFESLV::setup() or
//! OCFESLV::classify_pde().
class OCPlan
{
  // OCFESLV is the sole builder: it populates the draft, calls validate(),
  // and commits via commit_draft().  No other class may modify an OCPlan.
  friend class OCFESLV;

public:

  // ------------------------------------------------------------------
  // Public types
  // ------------------------------------------------------------------

  //! @brief Three-state lifecycle of the plan.
  enum class State { EMPTY, BUILDING, FROZEN };

  //! @brief The imposition type -- THE definition: OCFESLV::Options::ImpositionType is an alias of it,
  //! and OCFESLV::Options::IC_WEAK/IC_STRONG/IC_TRACE re-export these values.
  enum ImpositionSentinel {
    IMPOSITION_IC_WEAK   = 0,  //!< Equation residuals are enforced at ALL nodes (including interface nodes), and SAT penalty terms are added to interface-node residuals
    IMPOSITION_IC_STRONG = 1,  //!< Tau-eliminated reduction of IC_TRACE: exact continuity rows are retained and physical residual rows are projected through a setup-frozen Schur complement of the trace/tau receiver block.
    IMPOSITION_IC_TRACE  = 2   //!< Trace/tau formulation: keep all physical rows, add exact continuity rows, and add one multiplier per exact row distributed through the IC_WEAK SAT receiver graph
  };

  //! @brief The coupling a receiver-edge rescue uses when the model does not define one: C2's refused case
  //! (the receiving row's coefficient of the claimed state is not a nonzero constant -- a CONVENTION, and the
  //! edge is LABELLED coupling_undefined) and the geometric-C0 fallback (an identity coupling by construction).
  //! Constant (R5a); it was Options::RESCUE_COUPLING / CRONOS_RESCUE_COUPLING, a falsification
  //! knob whose corpus reach sweep M7b measured.
  static constexpr double kRescueFallbackCoupling = 1.0;

  //! @brief Serialisable exact-interface claim key (tuple form).
  //!
  //! Stored in maps/sets that must survive FFGraph::insert() remapping.
  //! All fields are plain integer ids so the key is trivially copyable
  //! and comparable without touching FFVar or DAG pointers.
  typedef std::tuple<int,size_t,size_t,size_t,size_t,size_t,int>
    t_FrozenInterfaceClaimKey;

  //! @brief Internal claim key used during plan construction.
  //!
  //! Converted to/from t_FrozenInterfaceClaimKey by the freeze/unfreeze
  //! helpers.  The eqn_id field carries numeric_limits<size_t>::max()
  //! for claim-centred (IC_TRACE/IC_STRONG) slots.
  struct t_InterfaceClaimKey {
    int    block_id = 0;
    size_t state_id = 0;
    size_t eqn_id   = 0;
    size_t dom_id   = 0;
    size_t iel_lo   = 0;
    size_t face     = 0;
    int    kind     = 0;
    bool operator<( t_InterfaceClaimKey const& o ) const noexcept
    {
      if( block_id != o.block_id ) return block_id < o.block_id;
      if( state_id != o.state_id ) return state_id < o.state_id;
      if( eqn_id   != o.eqn_id   ) return eqn_id   < o.eqn_id;
      if( dom_id   != o.dom_id   ) return dom_id   < o.dom_id;
      if( iel_lo   != o.iel_lo   ) return iel_lo   < o.iel_lo;
      if( face     != o.face     ) return face     < o.face;
      return kind < o.kind;
    }
  };

  //! @brief Interface-condition type enum (mirrors OCBase::Options).
  //! Imported here so plan structs can carry a sat_type field without
  //! a forward-declaration dependency on OCBase::Options.
  //! THE definition -- OCFESLV::Options aliases it (and re-exports its enumerators), so there is no
  //! longer a second copy to keep in sync.  (Its comment used to say it must stay in sync with
  //! OCBase::Options::InterfaceType, a path that had long since ceased to exist.)
  //! (IC_CONS = 4, the conservative flux-form condition, retired -- see OCFESLV's release notes.)
  enum InterfaceType { IC_VALUE = 0, IC_FLUX = 1, IC_UPWIND = 2, IC_AUTO = 3 };

  //! @brief Claim-kind constants (C0 value, C1 flux).
  enum { CLAIM_EXACT_C0 = 0, CLAIM_EXACT_C1 = 1, CLAIM_EXACT_ROW = CLAIM_EXACT_C0 };

  //! @brief Provenance of a receiver edge.
  //!
  //! Structural field preserved through the plan pipeline.  Used for
  //! diagnostics, IC_STRONG pivot preference, and potential future
  //! per-source penalty scaling.
  enum class InterfaceEdgeSource : int {
    PRINCIPAL_SYMBOL    = 0,   //!< Ordinary principal-symbol SAT coupling
    AUXILIARY_ALGEBRAIC = 1,   //!< Auxiliary state with algebraic (non-derivative) coupling
    GEOMETRIC_C0        = 2    //!< Geometric C0 fallback for degenerate tensor-product directions
  };

  // ------------------------------------------------------------------
  // Receiver-edge topology
  // ------------------------------------------------------------------

  //! @brief Setup-time receiver-edge classification shared by all three modes.
  //!
  //! Each candidate weak-SAT receiver edge is classified before exact trace
  //! rows and tau variables are generated.  IC_WEAK keeps weak-only edges
  //! (their contribution vanishes when the jump is zero).  IC_TRACE may
  //! introduce a free multiplier only for trace-replaceable edges.  Hard
  //! physical rows and tensor locations on a tangential physical boundary
  //! are deliberately weak-only: a free multiplier there could mask
  //! boundary/initial/constitutive residuals.
  struct InterfaceReceiverEdge {
    size_t              row_id              = std::numeric_limits<size_t>::max();
    t_InterfaceClaimKey claim_key{};
    InterfaceType       sat_type            = IC_VALUE;
    InterfaceEdgeSource source              = InterfaceEdgeSource::PRINCIPAL_SYMBOL;
    size_t              receiver_dom_id     = std::numeric_limits<size_t>::max();
    size_t              receiver_iel        = std::numeric_limits<size_t>::max();
    size_t              receiver_inode      = std::numeric_limits<size_t>::max();
    size_t              receiver_flat       = std::numeric_limits<size_t>::max(); //!< pre-boundary-collapse tensor flat used for SAT/tau injection into buf[]
    size_t              receiver_emit_flat  = std::numeric_limits<size_t>::max(); //!< post-boundary-collapse emitted residual flat used for IC_STRONG row ownership
    std::vector<std::pair<size_t,size_t>> receiver_block_el;
    double              coupling            = 0.0;
    double              orientation         = 0.0;
    bool                weak_sat_allowed    = false;
    bool                trace_tau_allowed   = false;
    bool                hard_physical_row   = false;
    bool                tangential_physical_boundary = false;
    //! @brief C2 REFUSED to derive this edge's coupling from the receiving row (the row's coefficient of the
    //! claimed state is not a nonzero constant), so `coupling` is the forced fallback, not a model quantity.
    //! Set only where the rescue fired and C2 ran; read by CRONOS_WEAK_TAU_RESCUED=2.
    bool                coupling_undefined  = false;
  };

  // ------------------------------------------------------------------
  // Executable consumer terms
  // ------------------------------------------------------------------

  //! @brief Executable weak-SAT / tau-injection consumer term.
  //!
  //! Row/claim based rather than value based.  Numeric evaluation
  //! computes trace jumps with the active arithmetic type, but
  //! topology/receiver eligibility is frozen at setup time and never
  //! recomputed during eval()/deriv().
  //!
  //! Mode semantics of individual fields:
  //!   prefactor       — used by IC_WEAK (= sigma * tau * coupling * orientation).
  //!   trace_prefactor — used by IC_TRACE and IC_STRONG (= prefactor * INTERFACE.TRACE_TAU_SCALE).
  //!   lambda_col      — column index of the tau/lambda variable (IC_TRACE, IC_STRONG).
  //!   lambda_index    — 0-based index into the trace-constraint block (IC_STRONG = lambda_col;
  //!                     IC_TRACE = lambda_col - trace_var_offset).
  //!   trace_tau_allowed — true only for receiver edges that carry an independent multiplier.
  struct FrozenWeakSatTerm {
    size_t              row_id             = std::numeric_limits<size_t>::max();
    t_InterfaceClaimKey claim_key{};
    InterfaceType       sat_type           = IC_VALUE;
    InterfaceEdgeSource source             = InterfaceEdgeSource::PRINCIPAL_SYMBOL;
    size_t              receiver_dom_id    = std::numeric_limits<size_t>::max();
    size_t              receiver_iel       = std::numeric_limits<size_t>::max();
    size_t              receiver_inode     = std::numeric_limits<size_t>::max();
    size_t              receiver_flat      = std::numeric_limits<size_t>::max(); //!< pre-boundary-collapse tensor flat used for SAT/tau injection into buf[]
    size_t              receiver_emit_flat = std::numeric_limits<size_t>::max(); //!< post-boundary-collapse emitted residual flat used for IC_STRONG row ownership
    std::vector<std::pair<size_t,size_t>> receiver_block_el;
    double              coupling           = 0.0;
    double              orientation        = 0.0;
    bool                trace_tau_allowed  = false;
    bool                uses_c1            = false;
    double              prefactor          = 0.0;  //!< IC_WEAK scalar weight (diagnostic / zero-check)
    double              trace_prefactor    = 0.0;  //!< IC_TRACE / IC_STRONG
    size_t              lambda_col         = std::numeric_limits<size_t>::max(); //!< IC_TRACE / IC_STRONG
    size_t              lambda_index       = std::numeric_limits<size_t>::max(); //!< IC_TRACE / IC_STRONG
    //! @brief IC_WEAK with multipliers only: this term's claim is imposed EXACTLY (it has a trace row and a
    //! multiplier column), so the weak evaluation injects trace_prefactor * var[lambda_col] when
    //! trace_tau_allowed, and nothing otherwise -- the IC_TRACE semantics -- instead of the penalty.
    bool                exact_claim        = false;
    //! @brief C2b: this rescued receiver term is a duplicate of a natural receiver and does not
    //! carry the multiplier (trace_tau_allowed cleared by the rule, not by the edge gate).
    bool                multiplier_silenced = false;
    //! @brief Pre-compiled sparse SAT coefficient list for IC_WEAK evaluation.
    //!
    //! Each entry (col, coeff) satisfies coeff = prefactor * w_lagrange[k].
    //! The sign convention matches the C0 jump u_lo(+1) - u_hi(-1): left-element
    //! nodes contribute positively, right-element nodes negatively.  For C1 claims
    //! the weights are additionally scaled by ±2/elem_width.
    //! Empty when prefactor == 0 (tau or coupling is zero).
    std::vector<std::pair<size_t,double>> sat_coeffs;
  };

  // ------------------------------------------------------------------
  // Block-key type and hash for two-level dispatch index
  // ------------------------------------------------------------------

  //! @brief Canonical element-block key: vector of (domain_var_id, element_index).
  //!
  //! Both plan-build and eval iterate domains in lt_FFVar order, so the
  //! entries are always sorted by domain_var_id.  This makes element-wise
  //! comparison sufficient for equality and hashing.
  typedef std::vector<std::pair<size_t,size_t>> BlockKey;

  //! @brief Hash functor for BlockKey.
  struct BlockKeyHash {
    size_t operator()( BlockKey const& k ) const noexcept
    {
      size_t h = k.size();
      for( auto const& p : k ){
        h ^= std::hash<size_t>{}( p.first  ) + 0x9e3779b9 + ( h << 6 ) + ( h >> 2 );
        h ^= std::hash<size_t>{}( p.second ) + 0x9e3779b9 + ( h << 6 ) + ( h >> 2 );
      }
      return h;
    }
  };

  //! @brief Inner map type: BlockKey -> term indices.
  typedef std::unordered_map<BlockKey,std::vector<size_t>,BlockKeyHash>
    t_BlockIndex;

  //! @brief Outer map type: row_id -> (BlockKey -> term indices).
  typedef std::unordered_map<size_t,t_BlockIndex> t_WeakSatIndex;

  //! @brief Executable exact trace constraint row (IC_TRACE / IC_STRONG).
  //!
  //! Emitted before physical residual rows.  The sparse coefficient vector
  //! holds (column_index, coefficient) pairs in the collocated variable
  //! vector, so evaluation applies only a sparse dot product without
  //! rebuilding claim maps, decoding faces, or evaluating Lagrange traces.
  struct FrozenTraceConstraintTerm {
    size_t row = std::numeric_limits<size_t>::max();
    t_InterfaceClaimKey edge{};
    std::vector<std::pair<size_t,double>> coeff;
  };

  //! @brief Stable key for a physical residual tensor entry (IC_STRONG).
  //!
  //! Used by the IC_STRONG Schur projection to identify pivot rows and
  //! projected non-pivot rows by their equation index, element-block
  //! coordinates, and flat position.
  struct FrozenPhysicalResidualRowKey {
    size_t row_id = std::numeric_limits<size_t>::max();
    std::vector<std::pair<size_t,size_t>> block_el;
    size_t flat   = std::numeric_limits<size_t>::max();
    bool operator<( FrozenPhysicalResidualRowKey const& o ) const noexcept
    {
      if( row_id   != o.row_id   ) return row_id   < o.row_id;
      if( block_el != o.block_el ) return block_el < o.block_el;
      return flat < o.flat;
    }
    bool operator==( FrozenPhysicalResidualRowKey const& o ) const noexcept
    { return row_id == o.row_id && block_el == o.block_el && flat == o.flat; }
  };

  //! @brief IC_STRONG tau-eliminated Schur projection plan.
  //!
  //! Setup selects a set of pivot rows from the tau receiver block B,
  //! inverts the pivot sub-matrix B_S, and pre-computes the correction
  //! coefficient vectors y = B_k B_S^{-1} for every non-pivot receiver
  //! row k.  Evaluation emits  R_k - y * R_S  for non-pivot rows and
  //! drops pivot rows entirely, replacing them with the exact trace rows
  //! emitted by LOOP 1.
  //!
  //! The projection coefficients reference pivot rows by dense index
  //! (0..ntrace-1) rather than by key, so the eval-time inner loop
  //! uses flat vector indexing with no map lookups.
  struct FrozenStrongTauElimination {
    bool   ready     = false;
    size_t n_trace   = 0;
    size_t rank      = 0;
    double pivot_min = 0.0;
    double pivot_max = 0.0;
    std::vector<FrozenPhysicalResidualRowKey> pivot_rows;
    std::map<FrozenPhysicalResidualRowKey,size_t> pivot_index;
    //! @brief Correction coefficients for non-pivot receiver rows.
    //!
    //! For each non-pivot row key, stores a sparse vector of
    //! (pivot_index, coefficient) pairs.  At eval time the corrected
    //! residual is:  R_k - sum_p coeff_p * pivot_value[pivot_index_p]
    //! where pivot_value[] is a dense array of pivot-row residual values
    //! collected during the LOOP 2 traversal.
    std::map<FrozenPhysicalResidualRowKey,
             std::vector<std::pair<size_t,double>>> projection;
    void clear() noexcept
    {
      ready = false; n_trace = rank = 0; pivot_min = pivot_max = 0.0;
      pivot_rows.clear(); pivot_index.clear(); projection.clear();
    }
    bool is_pivot( FrozenPhysicalResidualRowKey const& r ) const
    { return pivot_index.find( r ) != pivot_index.end(); }
  };

  // ------------------------------------------------------------------
  // Setup-time option snapshot
  // ------------------------------------------------------------------

  //! @brief Interface-plan options captured when executable terms are frozen.
  //!
  //! Public OCFESLV::options can legally be modified by callers after setup().
  //! Because FrozenWeakSatTerm::prefactor and FrozenWeakSatTerm::trace_prefactor already bake
  //! in INTERFACE.SAT_SIGMA0, sigma1 and INTERFACE.TRACE_TAU_SCALE, evaluation must reject a stale
  //! (sigma1 comes from CRONOS_SAT_SIGMA1, fixed per process, so in practice only the other two can move)
  //! plan rather than silently using coefficients from an older setup().
  //! The enum-valued options are stored as ints to keep OCPlan independent of
  //! the full OCFESLV::Options definition.
  struct FrozenOptions {
    int    imposition_type = IMPOSITION_IC_WEAK;
    int    interface_type  = IC_VALUE;
    double sat_sigma0      = 1.0;
    double sat_sigma1      = 1.0;
    double trace_tau_scale = 1.0;
  };

  //! @brief Setup-time diagnostic counters for interface-plan construction.
  //!
  //! These are not consumed by eval()/deriv(); they are frozen solely so that
  //! tests and users can audit the trace/SAT graph after setup().
  struct InterfacePlanBuildDiagnostics {
    size_t geometric_candidates = 0;          //!< Deferred geometric-C0 candidates collected before global filtering.
    size_t geometric_discarded_covered = 0;   //!< Candidates suppressed because a non-geometric edge already covered the same state/domain/interface.
    size_t geometric_emitted = 0;             //!< Candidates surviving global filtering and emitted as receiver edges.
  };

  // ------------------------------------------------------------------
  // Transactional draft (populated exclusively by OCFESLV::setup())
  // ------------------------------------------------------------------

  //! @brief Transactional staging area for a new frozen plan.
  //!
  //! OCFESLV::setup() builds all plan consumers into this local object,
  //! calls OCPlan::validate_draft(), and commits atomically via
  //! OCPlan::commit_draft().  A failed or exception-throwing build
  //! discards the draft without touching the live plan, so eval()/deriv()
  //! always consume either the previous sealed plan or none at all.
  struct InterfacePlanDraft {
    State state = State::EMPTY;
    bool  frozen = false;
    FrozenOptions options;
    std::map<t_FrozenInterfaceClaimKey,size_t>             claims;
    std::map<t_FrozenInterfaceClaimKey,std::vector<size_t>> tau_claims;
    std::vector<InterfaceReceiverEdge>                     receiver_edges;
    std::vector<FrozenWeakSatTerm>                         weak_sat_terms;
    bool                                                   weak_sat_terms_ready    = false;
    std::vector<FrozenTraceConstraintTerm>                 trace_constraints;
    bool                                                   trace_constraints_ready = false;
    std::set<FrozenPhysicalResidualRowKey>                 exact_replacement_rows;
    //! @brief IC_STRONG hybrid: continuity claims whose tau is kept EXPLICIT.
    //!
    //! These are demoted (rank-dependent receiver column) claims whose only
    //! replacement candidate is an auxiliary LINK row.  Consuming that LINK
    //! removes the aux's defining equation -> the reduced primal inherits a
    //! null direction and the aux drifts (the |dAux| corner error).  Instead of
    //! exact-replacing (consuming the LINK), these claims are kept the IC_TRACE
    //! way: the LINK is retained, an explicit tau unknown is appended to var[],
    //! and the continuity is imposed as a real row.  They are a subset of the
    //! n_trace_var tau columns and are EXCLUDED from the IC_STRONG Schur block.
    //! Empty for IC_WEAK / IC_TRACE and for IC_STRONG when no LINK sacrifice is
    //! needed, so a default plan is unaffected.
    std::set<t_FrozenInterfaceClaimKey>                    explicit_tau_claims;

    FrozenStrongTauElimination                             strong_tau_elim;
    InterfacePlanBuildDiagnostics                           diagnostics;
    size_t n_coll_eqn        = 0;
    size_t n_coll_var        = 0;
    size_t n_trace_var       = 0;
    //! @brief Count of n_trace_var taus kept explicit (IC_STRONG hybrid).
    //! Eliminated-via-Schur count = n_trace_var - n_trace_var_explicit.
    size_t n_trace_var_explicit = 0;
    size_t trace_var_offset  = 0;
    size_t n_trace_constraint = 0;
    //! @brief An IC_WEAK plan that imposes some claims exactly: n_trace_var multiplier columns after the
    //! collocation coefficients, one trace row each, held to the IC_TRACE counts by validate_draft.
    bool   weak_multipliers  = false;
    //! @brief With weak_multipliers: an exact claim that has a natural (principal-symbol) receiver able to take
    //! the multiplier feeds it ONLY there; its rescued receivers get neither multiplier nor penalty.
    bool   weak_natural_only = false;
  };

  // ------------------------------------------------------------------
  // Constructor / state query
  // ------------------------------------------------------------------

  //! @brief Default constructor — produces an empty, unsealed plan.
  OCPlan() = default;

  //! @brief True when the plan has been fully built and sealed.
  bool frozen() const noexcept { return _frozen; }

  //! @brief Lifecycle state of the plan.
  State state() const noexcept { return _state; }

  //! @brief True when the plan is sealed and evaluation may proceed.
  bool ready() const noexcept
  { return _frozen && _state == State::FROZEN && !_build_allowed; }

  //! @brief Setup-time interface-plan options associated with the frozen plan.
  FrozenOptions const& frozen_options() const noexcept { return _options; }

  // ------------------------------------------------------------------
  // Dimension query (valid only after commit)
  // ------------------------------------------------------------------

  size_t n_coll_eqn()         const noexcept { return _n_coll_eqn; }
  size_t n_coll_var()         const noexcept { return _n_coll_var; }
  size_t n_trace_var()        const noexcept { return _n_trace_var; }
  size_t n_trace_var_explicit() const noexcept { return _n_trace_var_explicit; }
  size_t trace_var_offset()   const noexcept { return _trace_var_offset; }
  size_t n_trace_constraint() const noexcept { return _n_trace_constraint; }

  //! @brief True when an IC_WEAK plan imposes some claims exactly (multiplier columns + trace rows).
  bool weak_multipliers() const noexcept { return _weak_multipliers; }

  // ------------------------------------------------------------------
  // Consumer data accessors (valid only after commit; called from _eval_eqn)
  // ------------------------------------------------------------------

  //! @brief Frozen weak-SAT / tau-injection consumer terms.
  std::vector<FrozenWeakSatTerm> const& weak_sat_terms() const noexcept
  { return _weak_sat_terms; }

  //! @brief True when the weak-SAT term list has been fully populated.
  bool weak_sat_terms_ready() const noexcept { return _weak_sat_terms_ready; }

  // ------------------------------------------------------------------
  // Dispatch indices (built by commit_draft; consumed by _eval_eqn)
  // ------------------------------------------------------------------

  //! @brief Two-level dispatch index for IC_WEAK:
  //!   outer: row_id (equation index) -> inner map
  //!   inner: BlockKey (receiver_block_el) -> vector of term indices
  //!
  //! Pre-built at commit time.  _eval_eqn LOOP 2 projects the current
  //! ndx_el to a BlockKey once per block (O(ndim)), then dispatches
  //! via outer[row_id][block_key] to the O(1) exactly-matching terms.
  //! frozen_term_matches_block is no longer on the hot path.
  t_WeakSatIndex const& weak_sat_index() const noexcept
  { return _weak_sat_index; }

  //! @brief Two-level dispatch index for IC_TRACE: same structure as
  //!        weak_sat_index but contains only trace_tau_allowed terms.
  t_WeakSatIndex const& trace_sat_index() const noexcept
  { return _trace_sat_index; }

  //! @brief Frozen exact trace constraint rows (IC_TRACE / IC_STRONG).
  std::vector<FrozenTraceConstraintTerm> const& trace_constraints() const noexcept
  { return _trace_constraints; }

  //! @brief True when the trace constraint list has been fully populated.
  bool trace_constraints_ready() const noexcept { return _trace_constraints_ready; }

  //! @brief Physical rows replaced by exact trace rows (IC_TRACE).
  std::set<FrozenPhysicalResidualRowKey> const& exact_replacement_rows() const noexcept
  { return _exact_replacement_rows; }

  //! @brief IC_STRONG hybrid: continuity claims whose tau is kept explicit.
  std::set<t_FrozenInterfaceClaimKey> const& explicit_tau_claims() const noexcept
  { return _explicit_tau_claims; }

  //! @brief IC_STRONG tau-eliminated Schur projection plan.
  FrozenStrongTauElimination const& strong_projection() const noexcept
  { return _strong_tau_elim; }

  //! @brief Frozen receiver-edge topology (used by _store_frozen_interface_plan
  //!        and by deep_copy_from remapping; not consumed directly by _eval_eqn).
  std::vector<InterfaceReceiverEdge> const& receiver_edges() const noexcept
  { return _receiver_edges; }

  //! @brief Frozen exact-trace claim map (used by _load/_store helpers).
  std::map<t_FrozenInterfaceClaimKey,size_t> const& claims_frozen() const noexcept
  { return _claims; }

  //! @brief Frozen tau-column assignment map (used by build helpers).
  std::map<t_FrozenInterfaceClaimKey,std::vector<size_t>> const& tau_claims_frozen() const noexcept
  { return _tau_claims; }

  //! @brief Frozen setup-time diagnostic counters for interface-plan construction.
  InterfacePlanBuildDiagnostics const& diagnostics() const noexcept
  { return _diagnostics; }

  // ------------------------------------------------------------------
  // Static key helpers
  // ------------------------------------------------------------------

  //! @brief Serialise a t_InterfaceClaimKey to a storable tuple key.
  static t_FrozenInterfaceClaimKey freeze_claim_key
    ( t_InterfaceClaimKey const& k ) noexcept
  {
    return t_FrozenInterfaceClaimKey{
      k.block_id, k.state_id, k.eqn_id,
      k.dom_id,   k.iel_lo,   k.face,   k.kind };
  }

  //! @brief Deserialise a storable tuple key back to t_InterfaceClaimKey.
  static t_InterfaceClaimKey unfreeze_claim_key
    ( t_FrozenInterfaceClaimKey const& k ) noexcept
  {
    return t_InterfaceClaimKey{
      std::get<0>(k), std::get<1>(k), std::get<2>(k),
      std::get<3>(k), std::get<4>(k), std::get<5>(k), std::get<6>(k) };
  }

  // ------------------------------------------------------------------
  // Plan reset (called by OCFESLV::_reset() and OCFESLV::classify_pde())
  // ------------------------------------------------------------------

  //! @brief Reset the plan to the empty state.
  //!
  //! Called by OCFESLV whenever the symbolic model or classification
  //! changes and the current plan is therefore invalidated.
  void reset() noexcept
  {
    _state                    = State::EMPTY;
    _frozen                   = false;
    _build_allowed            = false;
    _options                  = FrozenOptions{};
    _diagnostics              = InterfacePlanBuildDiagnostics{};
    _claims.clear();
    _tau_claims.clear();
    _receiver_edges.clear();
    _weak_sat_terms.clear();
    _weak_sat_terms_ready     = false;
    _weak_sat_index.clear();
    _trace_sat_index.clear();
    _trace_constraints.clear();
    _trace_constraints_ready  = false;
    _exact_replacement_rows.clear();
    _explicit_tau_claims.clear();
    _strong_tau_elim.clear();
    _n_coll_eqn               = 0;
    _n_coll_var               = 0;
    _n_trace_var              = 0;
    _n_trace_var_explicit     = 0;
    _trace_var_offset         = 0;
    _n_trace_constraint       = 0;
    _weak_multipliers         = false;
  }

  // ------------------------------------------------------------------
  // Validation (called by OCFESLV::setup() before commit)
  // ------------------------------------------------------------------

  //! @brief Validate a completed draft against the supplied imposition type.
  //!
  //! Returns true when the draft is self-consistent and safe to commit.
  //! Emits diagnostics to std::cerr and returns false otherwise.
  //! The imposition_type argument is the value of
  //! OCFESLV::options.IMPOSITION_TYPE at the time setup() was called.
  template <typename ImpositionType>
  bool validate_draft( InterfacePlanDraft const& draft,
                       ImpositionType imposition_type ) const
  {
    int const imp = static_cast<int>( imposition_type );

    if( !draft.frozen || draft.state != State::FROZEN ){
      std::cerr << "OCPlan::validate_draft ** plan build completed without "
                   "sealing the draft (state is not FROZEN)" << std::endl;
      return false;
    }
    if( !draft.weak_sat_terms_ready ){
      std::cerr << "OCPlan::validate_draft ** weak/trace SAT consumer terms "
                   "are not ready" << std::endl;
      return false;
    }
    if( !draft.trace_constraints_ready ){
      std::cerr << "OCPlan::validate_draft ** exact trace constraint rows "
                   "are not ready" << std::endl;
      return false;
    }
    // An IC_WEAK plan with multipliers is held to the IC_TRACE counts.
    bool const weak_mult = ( imp == IMPOSITION_IC_WEAK && draft.weak_multipliers );
    if( ( imp == IMPOSITION_IC_TRACE || imp == IMPOSITION_IC_STRONG || weak_mult )
     && draft.n_trace_constraint != draft.trace_constraints.size() ){
      std::cerr << "OCPlan::validate_draft ** exact trace row count mismatch: "
                   "n_trace_constraint=" << draft.n_trace_constraint
                << " rows=" << draft.trace_constraints.size() << std::endl;
      return false;
    }
    if( imp == IMPOSITION_IC_TRACE || weak_mult ){
      if( draft.n_coll_var != draft.trace_var_offset + draft.n_trace_var ){
        std::cerr << "OCPlan::validate_draft ** IC_TRACE variable count "
                     "mismatch: nvar=" << draft.n_coll_var
                  << " offset=" << draft.trace_var_offset
                  << " ntau=" << draft.n_trace_var << std::endl;
        return false;
      }
    }
    if( imp == IMPOSITION_IC_TRACE || imp == IMPOSITION_IC_STRONG || weak_mult ){
      if( draft.n_trace_constraint
       != draft.n_trace_var + draft.exact_replacement_rows.size() ){
        std::cerr << "OCPlan::validate_draft ** IC_TRACE/IC_STRONG exact/dual "
                     "count mismatch: exact_rows=" << draft.n_trace_constraint
                  << " dual_tau=" << draft.n_trace_var
                  << " exact_replacements="
                  << draft.exact_replacement_rows.size() << std::endl;
        return false;
      }
    }
    if( imp == IMPOSITION_IC_STRONG && !draft.strong_tau_elim.ready ){
      std::cerr << "OCPlan::validate_draft ** IC_STRONG tau-elimination "
                   "projection is not ready" << std::endl;
      return false;
    }
    return true;
  }

  // ------------------------------------------------------------------
  // Commit (called by OCFESLV::setup() after successful validation)
  // ------------------------------------------------------------------

  //! @brief Atomically publish a validated draft as the live sealed plan.
  //!
  //! Move-commits all draft members into the plan, resets _build_allowed
  //! to false, and seals the plan (state = FROZEN).  After this call
  //! ready() returns true and _eval_eqn() may consume the plan.
  void commit_draft( InterfacePlanDraft&& draft ) noexcept
  {
    _state                   = draft.state;
    _frozen                  = draft.frozen;
    _build_allowed           = false;
    _options                 = draft.options;
    _claims                  = std::move( draft.claims );
    _tau_claims              = std::move( draft.tau_claims );
    _receiver_edges          = std::move( draft.receiver_edges );
    _weak_sat_terms          = std::move( draft.weak_sat_terms );
    _weak_sat_terms_ready    = draft.weak_sat_terms_ready;
    _trace_constraints       = std::move( draft.trace_constraints );
    _trace_constraints_ready = draft.trace_constraints_ready;
    _exact_replacement_rows  = std::move( draft.exact_replacement_rows );
    _explicit_tau_claims     = std::move( draft.explicit_tau_claims );
    _strong_tau_elim         = std::move( draft.strong_tau_elim );
    _diagnostics             = draft.diagnostics;
    _n_coll_eqn              = draft.n_coll_eqn;
    _n_coll_var              = draft.n_coll_var;
    _n_trace_var             = draft.n_trace_var;
    _n_trace_var_explicit    = draft.n_trace_var_explicit;
    _trace_var_offset        = draft.trace_var_offset;
    _n_trace_constraint      = draft.n_trace_constraint;
    _weak_multipliers        = draft.weak_multipliers;
    // Build the dispatch indices from the newly committed term list.
    _build_weak_sat_term_indices();
  }

  // ------------------------------------------------------------------
  // Deep-copy helpers (called by OCFESLV::deep_copy_from())
  // ------------------------------------------------------------------

  //! @brief Copy the plan from src, remapping equation ids through the
  //!        supplied map produced by the DAG deep-copy operation.
  //!
  //! State, domain and block ids are preserved by FFGraph::insert, so
  //! only the eqn_id slot of t_FrozenInterfaceClaimKey needs remapping.
  //! The sentinel value numeric_limits<size_t>::max() is preserved (it
  //! marks claim-centred IC_TRACE/IC_STRONG slots that have no eqn_id).
  void copy_from( OCPlan const& src,
                  std::map<size_t,size_t> const& eqn_id_map ) noexcept
  {
    reset();
    _state          = src._state;
    _frozen         = src._frozen;
    _build_allowed  = false;
    _options        = src._options;
    _diagnostics    = src._diagnostics;
    _n_coll_eqn     = src._n_coll_eqn;
    _n_coll_var     = src._n_coll_var;
    _n_trace_var    = src._n_trace_var;
    _n_trace_var_explicit = src._n_trace_var_explicit;
    _trace_var_offset   = src._trace_var_offset;
    _n_trace_constraint = src._n_trace_constraint;
    _weak_multipliers   = src._weak_multipliers;

    auto remap_eqn_id = [&]( size_t id ) -> size_t {
      auto it = eqn_id_map.find( id );
      return it == eqn_id_map.end() ? id : it->second;
    };

    auto remap_key = [&]( t_FrozenInterfaceClaimKey const& k )
      -> t_FrozenInterfaceClaimKey
    {
      return t_FrozenInterfaceClaimKey{
        std::get<0>(k), std::get<1>(k),
        remap_eqn_id( std::get<2>(k) ),
        std::get<3>(k), std::get<4>(k), std::get<5>(k), std::get<6>(k) };
    };

    for( auto const& kv : src._claims )
      _claims[ remap_key( kv.first ) ] = kv.second;
    for( auto const& kv : src._tau_claims )
      _tau_claims[ remap_key( kv.first ) ] = kv.second;
    _explicit_tau_claims.clear();
    for( auto const& k : src._explicit_tau_claims )
      _explicit_tau_claims.insert( remap_key( k ) );

    _receiver_edges = src._receiver_edges;
    for( auto& edge : _receiver_edges ){
      auto fk = freeze_claim_key( edge.claim_key );
      edge.claim_key = unfreeze_claim_key( remap_key( fk ) );
    }

    _weak_sat_terms = src._weak_sat_terms;
    for( auto& term : _weak_sat_terms ){
      auto fk = freeze_claim_key( term.claim_key );
      term.claim_key = unfreeze_claim_key( remap_key( fk ) );
    }
    _weak_sat_terms_ready = src._weak_sat_terms_ready;

    _trace_constraints = src._trace_constraints;
    for( auto& term : _trace_constraints ){
      auto fk = freeze_claim_key( term.edge );
      term.edge = unfreeze_claim_key( remap_key( fk ) );
    }
    _trace_constraints_ready = src._trace_constraints_ready;

    _exact_replacement_rows = src._exact_replacement_rows;
    _strong_tau_elim        = src._strong_tau_elim;

    // Rebuild indices from the remapped term list — do not copy src indices
    // since row_id values are equation ids that may have been remapped above.
    _build_weak_sat_term_indices();
  }

private:

  // ------------------------------------------------------------------
  // Internal helpers
  // ------------------------------------------------------------------

  //! @brief Build two-level _weak_sat_index and _trace_sat_index from _weak_sat_terms.
  //!
  //! Called by commit_draft() and copy_from() after the term list is final.
  //! O(n_terms) time and space.  The two-level structure maps
  //!   row_id -> receiver_block_el -> [term indices]
  //! so that _eval_eqn dispatches in O(1) to the exactly-matching terms
  //! for a given (equation, element-block) pair, with no residual scan
  //! across other equations or other element interfaces.
  void _build_weak_sat_term_indices()
  {
    _weak_sat_index.clear();
    _trace_sat_index.clear();
    for( size_t i = 0; i < _weak_sat_terms.size(); ++i ){
      auto const& term = _weak_sat_terms[i];
      BlockKey const& bk = term.receiver_block_el;
      _weak_sat_index[ term.row_id ][ bk ].push_back( i );
      if( term.trace_tau_allowed )
        _trace_sat_index[ term.row_id ][ bk ].push_back( i );
    }
  }

  // ------------------------------------------------------------------
  // Live plan members (populated by commit_draft, consumed by _eval_eqn)
  // ------------------------------------------------------------------

  State  _state         = State::EMPTY;
  bool   _frozen        = false;
  bool   _build_allowed = false;   //!< True only inside OCFESLV::_build_interface_plan_draft
  FrozenOptions _options;          //!< Setup-time options baked into executable terms
  InterfacePlanBuildDiagnostics _diagnostics; //!< Setup-time diagnostic counters.

  std::map<t_FrozenInterfaceClaimKey,size_t>              _claims;
  std::map<t_FrozenInterfaceClaimKey,std::vector<size_t>> _tau_claims;
  std::vector<InterfaceReceiverEdge>                      _receiver_edges;
  std::vector<FrozenWeakSatTerm>                          _weak_sat_terms;
  bool                                                    _weak_sat_terms_ready    = false;
  std::vector<FrozenTraceConstraintTerm>                  _trace_constraints;
  bool                                                    _trace_constraints_ready = false;
  std::set<FrozenPhysicalResidualRowKey>                  _exact_replacement_rows;
  std::set<t_FrozenInterfaceClaimKey>                     _explicit_tau_claims;
  FrozenStrongTauElimination                              _strong_tau_elim;

  //! @brief Two-level dispatch: row_id → BlockKey → term indices (all terms).
  t_WeakSatIndex _weak_sat_index;
  //! @brief Two-level dispatch: row_id → BlockKey → term indices (trace_tau_allowed only).
  t_WeakSatIndex _trace_sat_index;

  size_t _n_coll_eqn        = 0;
  size_t _n_coll_var        = 0;
  size_t _n_trace_var       = 0;
  size_t _n_trace_var_explicit = 0;
  size_t _trace_var_offset  = 0;
  size_t _n_trace_constraint = 0;
  bool   _weak_multipliers   = false;

}; // class OCPlan


//============================================================================================================
// BUILDER INPUTS AND DECISIONS -- PlanInput, PlanDecisions, PlanReport
//   (was ocplan_build.hpp)
//============================================================================================================

// ocplan_build.hpp -- interface-plan decision layer: CONTRACT (phase 0 skeleton)
//
// This header defines the boundary between OCFESLV (layout, classification, receiver
// edge synthesis, numeric materialisation, assembly) and the plan DECISION layer
// (which claims are constraints, how each is realised, which physical rows are
// replaced, which taus exist).  In phase 0 the legacy decision code is moved behind
// this contract unchanged; in phase 1 the global builder implements the same
// contract; in phase 2 the legacy builder is removed.
//
// Rules of this file:
//   * no I/O, no environment variables, no OCFESLV access.  Everything the builder
//     needs arrives in PlanInput; everything it decides leaves in PlanDecisions;
//     everything it wants to say leaves in PlanReport.  A separate reporter prints.
//   * plain data.  Keys are OCPlan's own key types so that OCFESLV can pack/unpack
//     without conversion; row roles are mirrored as an int enum (EqnRole lives in
//     ocfeslv; the bridge static_assert sits in OCFESLV, as for InterfaceType).
//   * comments describe the contract as it is; history lives in the FINDING docs.

// ------------------------------------------------------------------------------
// Row roles, mirrored from OCFESLV::EqnRole (bridged by static_assert in OCFESLV).
// ------------------------------------------------------------------------------
enum class PlanRowRole : int {
  AUTO = 0, INTERIOR, INITIAL, BOUNDARY, INTERFACE, LINK, SURFACE, DIAGNOSTIC
};

// ------------------------------------------------------------------------------
// Options that reach the decision layer.  Copied from OCFESLV::Options by OCFESLV; the
// builder never sees the full option set.
// ------------------------------------------------------------------------------
struct PlanOptions
{
  enum class Imposition   : int { WEAK = 0, TRACE = 1, STRONG = 2 };

  // 4.0 (2026-10-09): the rev153c experiments -- realisation (EXACT_FIRST / TAU_FIRST), causal (KEEP_UPSTREAM /
  // SYMMETRIC), aux_implied, root_hi, face_desc -- and the builder selector (LEGACY / GLOBAL) are gone: nothing set
  // them since Phase 3 retired their variables (batch 1), so the builder always ran their defaults, which it now
  // hard-codes: exact realisation first, causal roots upstream, implied cross-direction auxiliary edges dropped,
  // trees rooted at the low DOF, faces ascending.  display_level and fault_inject_bad_keep_explicit (set, never
  // read: the legacy builder that printed with them was deleted in rev264) are gone too.
  Imposition  imposition  = Imposition::WEAK;
  bool        trace_dofs  = false;    // rev155a (shadow): non-nodal trace functionals become abstract DOFs (design item 2)
  bool        weak_multipliers = false;   // IC_WEAK realises PlanInput::rescued_claims as TAU (CRONOS_WEAK_TAU_RESCUED)
  bool        weak_natural_only = false;  // ... and a claim WITH a natural receiver feeds its multiplier only there (CRONOS_WEAK_TAU_NATURAL)
  bool        force_bad_keep_explicit = false;   // attribution testing (CRONOS_FORCE_BAD_KEEP_EXPLICIT): invert the first STRONG explicit choice so the adequacy check must catch it    // experiment: within a seam, add edges in DESCENDING transverse face order (drop low-side cycle edges)
};

// ------------------------------------------------------------------------------
// PlanInput: everything the decision layer may read.
// ------------------------------------------------------------------------------
struct PlanInput
{
  typedef OCPlan::t_InterfaceClaimKey        ClaimKey;
  typedef OCPlan::t_FrozenInterfaceClaimKey  FrozenClaimKey;
  typedef OCPlan::FrozenPhysicalResidualRowKey RowKey;
  typedef OCPlan::InterfaceReceiverEdge      ReceiverEdge;

  PlanOptions options;

  //! Claims after classification, reduction and (if any) marching collapse:
  //! internal key -> occurrence count.  Same object the legacy pipeline builds.
  std::map<ClaimKey,size_t> claims;

  //! Receiver graph: for each claim, the residual rows it may inject into, with
  //! coupling/orientation and the weak/tau/hard flags.  Same object as
  //! draft.receiver_edges; the builder reads it, never mutates it (edge flag updates
  //! are returned in PlanDecisions and applied by OCFESLV).
  std::vector<ReceiverEdge> const* receiver_edges = nullptr;

  //! Equation roles by row_id (index = row_id; LINK rows are the reduce_order()
  //! auxiliary definitions and are never replaceable).
  std::vector<PlanRowRole> row_role;

  //! Names for reporting only (index = state id / domain id).  Empty string if unknown.
  std::map<size_t,std::string> state_name;
  std::map<size_t,std::string> dom_name;

  //! Layout facts the legacy layer queried through OCFESLV helpers.
  //! trace_claim_block: (state_id, dom_id) -> block id used to normalise slot keys.
  //! Packed for every (state, dom) pair occurring in claims or receiver edges.
  std::map<std::pair<size_t,size_t>,int> trace_claim_block;

  //! Model-structural fact: does residual row_id apply a Partial in direction dom_id?
  //! (row_id, dom_id) -> bool, for every row and every domain direction.  Legacy read
  //! this by scanning the DAG subgraph inside the decision layer.
  std::map<std::pair<size_t,size_t>,bool> eqn_differentiates_in;

  //! Claims that must be backed by an explicit tau (legacy: _strong_protected_claims).
  std::set<FrozenClaimKey> protected_claims;
  //! @brief claims IC_STRONG must keep EXPLICIT (var[] columns) because the trace projector's support
  //! reaches them; filled by OCFESLV on a promotion rebuild (Options::STRONG_PROJECT).
  std::set<FrozenClaimKey> promote_explicit;

  //! Claims excluded from allocation by the two-pass DROP (legacy:
  //! _strong_dropped_claims).  Phase 0 keeps the DROP inside the legacy builder; this
  //! field exists so the GLOBAL builder can be given the same exclusions if wanted.
  std::set<FrozenClaimKey> dropped_claims;

  //! Sizes the builder must respect / extend.
  size_t n_coll_eqn       = 0;
  size_t n_coll_var       = 0;
  size_t trace_var_offset = 0;

  //! Geometry for the GLOBAL builder (phase 1; unused by LEGACY), packed by OCFESLV from
  //! the claim coefficients and the receiver edges.
  //!
  //! A DOF is identified by its var[] column.  For a C0 claim whose constraint row has
  //! exactly two unit coefficients (nodal bases with endpoint nodes), dof_lo/dof_hi are
  //! the coincident columns on the two sides of the seam; otherwise nodal=false and the
  //! GLOBAL builder leaves the claim to legacy handling (reported).
  struct DofEdge {
    ClaimKey key;                 // the claim (internal key)
    size_t   dof_lo = 0, dof_hi = 0;
    bool     nodal  = false;
    int      direction_rank = 0;  // 0 = evolution/causal direction, then spatial in registration order
    bool     causal = false;
    //! the claim's two VALUE functionals, one per side (var[] column, weight),
    //! from the same Lagrange decode as the constraint row (lo = element iel_lo at +1,
    //! hi = iel_lo+1 at -1; C1: derivative stencils scaled by 2/width).  A side with a
    //! single unit entry is a nodal DOF (that column); any other side is an abstract
    //! trace DOF, identified by the GLOBAL builder (trace_dofs) through the nodal
    //! components: node-wise equal supports are one DOF.
    std::vector<std::pair<size_t,double>> lo_support, hi_support;
  };
  std::vector<DofEdge> dof_edges;

  //! Claims the LEGACY corner rule suppressed during receiver synthesis (never in
  //! `claims`, no receiver edges), with their DOF geometry.  Shadow use only: the GLOBAL
  //! forest is checked against the legacy rule (design §3.1); in applied GLOBAL mode the
  //! rule is disabled and these arrive as ordinary claims instead.
  std::vector<DofEdge> corner_suppressed_edges;
  //! pinned-suppressed pairs (both sides pinned pointwise by the same ALL-mask
  //! IC/BC).  Never realised; under trace_dofs they JOIN the identification union-find,
  //! because the two columns are equal by the IC/BC rows.  MEASURED PDE5: U's t-corners
  //! identify only through the t=0 node, which is IC-pinned.
  std::vector<DofEdge> pinned_edges;

  //! reduce_order() auxiliaries: aux state id -> (parent state id, directions its LINK
  //! differentiates in).  Design §3.5 (implied cross-direction auxiliary edges).
  struct AuxInfo { size_t parent = 0; std::set<size_t> diff_dirs; };
  std::map<size_t,AuxInfo> aux_info;

  //! Claims suppressed upstream as pinned by an ALL-mask IC/BC: the plan never realises
  //! them, so an implication that needs them (design §3.5) does not hold there.
  std::set<ClaimKey> unrealised_claims;

  //! @brief Claims whose interface coupling had to be SUBSTITUTED (the rescue fired), so no honest penalty
  //! weight exists.  Read only when options.weak_multipliers: under IC_WEAK the builder then realises each as
  //! TAU -- an exact constraint, which needs no weight -- provided it has a tau-allowed receiver.
  std::set<ClaimKey> rescued_claims;
  //! claims whose exact trace ROW was rejected by the tensor-corner MGS rank filter,
  //! i.e. linearly DEPENDENT on the accepted rows.  For homogeneous continuity rows this means
  //! the condition is IMPLIED by the rows the plan does impose, so the claim's locus IS
  //! enforced -- the opposite of unrealised_claims (pin-suppressed = genuinely not enforced).
  //! MEASURED (OCFE_PDE3, 2026-09-07): retiring the corner dedup makes the transverse corner
  //! rows sort earlier and be accepted, which renders 32 z-direction rows at the SAME nodes
  //! dependent; rank is 952 and the forest tree 940 either way, so the two accepted sets span
  //! the same space.  Design §3.5 must not read those rejections as broken implications.
  std::set<ClaimKey> rank_implied_claims;
  //! Nodes per element per domain direction (transverse element line of a face index).
  std::map<size_t,size_t> dom_n_node;

  //! True trace/tau column coefficient per receiver edge (aligned with receiver_edges):
  //! sigma-and-width-scaled, exactly as the materialisation assembles the tau columns.
  //! The STRONG explicit scan must use THIS, not coupling*orientation (PDE20f: 20 false
  //! explicits from unscaled values).
  std::vector<double> edge_trace_coeff;

  //! Per-block drop eligibility (OCFESLV::_block_drop_eligible).  STRONG scan policy: a
  //! DEPENDENT column in an ELIGIBLE, unprotected block stays INTERNAL so the redundancy
  //! machinery can W-validate and DROP it (PDE20f); ineligible or protected -> EXPLICIT
  //! (MBC pairs, PDE7 corners: dropping loses real conditions).
  std::map<int,bool> block_drop_eligible;

  //! Owner PDE rows of a state: INTERIOR equations whose DAG depends on the state or on
  //! one of its reduce_order() auxiliaries.  state_id -> row_ids.  Used so that a
  //! primitive's claim may consume its own PDE copy at a seam node even when that row is
  //! not a receiver of the claim (reduced states: the PDE contains the auxiliary).
  std::map<size_t,std::set<size_t>> owner_rows;

  //! Collocation-set fact per (row_id, dom_id): does the equation have a row at EVERY
  //! node of that direction (true), or is one end excluded (Z_NO_LB/NO_UB/INT: false)?
  //! A seam copy of an equation that lacks an end row in d must not be consumed for a
  //! d-claim (first-order operators: continuity must ADD a condition).
  std::map<std::pair<size_t,size_t>,bool> row_full_in_dir;

  //! Every receiver row copy of every claim, tagged with the DOF it sits at (the claim's
  //! DOF on the receiver's side), whether it may be replaced by a constraint row, and the
  //! coefficient a tau on that claim would carry into it.
  struct RowAtDof {
    RowKey  row;
    size_t  dof = 0;              // column of the claim's state at the receiver node
    size_t  edge_index = 0;       // index into receiver_edges
    bool    replaceable = false;  // role INTERIOR and not hard/tangential
    bool    side_hi     = false;  // rev155a: receiver lies on the claim's hi side (iel_lo+1)
                                  //   dof == SIZE_MAX for non-nodal claims: the builder resolves
                                  //   it from the edge's side DOF (trace_dofs) or skips the row.
    bool    tau_allowed = false;
    bool    sat_allowed = false;
    double  coeff = 0.;           // coupling * orientation
  };
  std::vector<RowAtDof> rows_at_dof;
};

// ------------------------------------------------------------------------------
// PlanDecisions: everything the decision layer may write.  OCFESLV materialises
// numeric terms from these (trace constraint coefficients, weak-SAT terms, Schur
// elimination) exactly as it does today.
// ------------------------------------------------------------------------------
struct PlanDecisions
{
  typedef PlanInput::ClaimKey       ClaimKey;
  typedef PlanInput::FrozenClaimKey FrozenClaimKey;
  typedef PlanInput::RowKey         RowKey;

  //! Claims kept for allocation after the two-pass DROP (constraint rows are
  //! materialised from this map, one per occurrence, in map order).  Legacy: the
  //! input claims minus dropped ones; GLOBAL: tree edges only.
  std::map<ClaimKey,size_t>                claims_kept;

  //! Tau allocation: frozen claim key -> tau column ids (var[] indices for IC_TRACE
  //! and for the IC_STRONG kept-explicit subset; internal Schur ids otherwise).
  std::map<FrozenClaimKey,std::vector<size_t>> tau_claims;
  std::set<FrozenClaimKey>                 explicit_tau_claims;
  std::map<FrozenClaimKey,int>             keep_explicit_reason;   // legacy diagnostic; phase 2 removes

  //! Physical row copies replaced by a constraint row.
  std::set<RowKey>                         exact_replacement_rows;

  //! Receiver-edge flag updates the builder requests (index into receiver_edges).
  //! Legacy: edges made tau-inactive by the DROP.  Applied by OCFESLV.
  std::vector<std::pair<size_t,bool>>      trace_tau_allowed_updates;

  //! Counts, computed by the builder so the invariants can be checked by OCFESLV:
  //!   n_trace_constraint == n_trace_var + exact_replacement_rows.size()
  size_t n_trace_var          = 0;
  size_t n_trace_var_explicit = 0;
  size_t n_trace_constraint   = 0;
  size_t n_coll_eqn           = 0;
  size_t n_coll_var           = 0;
};

// ------------------------------------------------------------------------------
// PlanReport: what the builder wants to say.  Printed by a separate reporter.
// ------------------------------------------------------------------------------
struct PlanReport
{
  typedef PlanInput::ClaimKey ClaimKey;
  typedef PlanInput::RowKey   RowKey;

  // Legacy counters (phase 0): kept so the extraction is byte-identical in what it
  // can say; the legacy builder still prints them itself through display_level.
  size_t geom_pivot_fallback = 0;
  size_t two_pass_dropped    = 0;
  size_t edges_made_tau_inactive = 0;

  // GLOBAL builder (phase 1).
  struct Cluster {
    size_t                       id = 0;
    std::vector<size_t>          dofs;
    std::vector<ClaimKey>        edges, tree_edges, dropped_edges;
    size_t                       root = 0;
    std::map<ClaimKey,int>       realisation;      // 0 = EXACT, 1 = TAU, 2 = REDUNDANT
    std::vector<RowKey>          consumed_rows;
    size_t                       tau_block_rank = 0, tau_block_cols = 0;
    std::vector<std::pair<ClaimKey,std::pair<RowKey,double>>> tau_receivers;   // diagnostics
    size_t                       alternatives_tried = 0;
    bool                         unresolved = false;
  };
  std::vector<Cluster> clusters;
  size_t               n_unresolved = 0;

  // Shadow summary (phase 1a/1b): what the GLOBAL builder saw and would do.
  struct Shadow {
    size_t n_claims = 0, n_edges_nodal = 0, n_edges_nonnodal = 0, n_edges_c1 = 0;
    size_t n_dofs = 0, n_components = 0, n_clusters = 0, max_cluster_dofs = 0, max_component_dofs = 0;
    size_t n_rows_at_dof = 0, n_rows_replaceable = 0;
    size_t n_tree_edges = 0, n_redundant_edges = 0;
    size_t n_corner_rule = 0, n_corner_rule_forest_drops = 0, n_forest_drops_not_corner_rule = 0;
    size_t n_dropped_by_drop = 0, n_aux_implied = 0;
    size_t n_exact = 0, n_tau = 0, n_sat = 0, n_explicit = 0;
    size_t n_weak_tau = 0, n_weak_tau_norecv = 0, n_weak_standalone_sat = 0;   // IC_WEAK with multipliers
    size_t n_promoted = 0;   // rev281: IC_STRONG taus promoted to explicit for the projection
    size_t n_clusters_rank_deficient = 0, n_rank_deficit = 0;
    size_t legacy_exact = 0, legacy_tau = 0, exact_rows_common = 0;
    // rev155a trace-DOF census (trace_dofs builds only)
    size_t tdof_sides = 0, tdof_abstract = 0, tdof_rows = 0, tdof_rows_tau = 0, tdof_standalone = 0;
    size_t tdof_corner_nonnodal = 0, tdof_corner_nonnodal_forest_drops = 0, tdof_pinned_joins = 0;
    std::vector<ClaimKey> forest_dropped;   // rev155a: keys of forest-redundant edges (own and corner-suppressed)
  } shadow;
};

// ------------------------------------------------------------------------------
// Builders.  Both are pure functions of PlanInput.
// ------------------------------------------------------------------------------
// rev264: class OCPlanLegacyBuilder DELETED together with ocplan_build_legacy.hpp (954 lines).  It was
// called from exactly one place and its decisions were discarded on every plan; the corpus confirms GLOBAL's
// drop set matches it on all 5 programs that reach the redundancy-drop cycle.

class OCPlanGlobalBuilder
{
public:
  //! @brief Decide the plan.  Reads @p in (claims, seam topology, options), writes the realisation of every
  //! claim and the dropped set to @p out, and the counts and shadow statistics to @p rep.
  //! @return false only on an internal inconsistency, with the message in @p rep for the caller to raise;
  //! a model that simply has no claims yields an empty plan and true.
  static bool build( PlanInput const& in, PlanDecisions& out, PlanReport& rep );
};


//============================================================================================================
// THE BUILDER -- spanning-forest realisation (was "GLOBAL"; it is now the only one)
//   (was ocplan_build_global.hpp)
//============================================================================================================

//! @brief The interface-plan decision layer: given the claims and the seam topology, decide how each
//! continuity claim is REALISED and which redundant ones are dropped.
//!
//! WHAT IT DECIDES, per claim:
//!   EXACT  the claim replaces the receiver's row copy at the child of the seam edge;
//!   TAU    the claim becomes a trace constraint with its own multiplier column;
//!   SAT    the claim becomes a penalty term weighted by the interface coupling (IC_WEAK);
//!   DROP   the claim is redundant given the rest and is not realised at all.
//!
//! WHAT IT GUARANTEES:
//!   - the realisation is a function of the model and the claim set alone: no I/O, no environment reads,
//!     no dependence on solver state, so the same inputs always give the same plan;
//!   - within a cluster -- every row copy at one physical point, across all states, plus the edges touching
//!     them -- the tau block is rank-checked, with bounded enumeration over root and dropped-edge choices,
//!     so a plan is never committed with a tau block it knows to be deficient;
//!   - the forest per state-component is built by Kruskal in a deterministic order (evolution edges first
//!     under KEEP_UPSTREAM), so the plan does not depend on container iteration order.
//!
//! OBJECTS: a DOF is a var[] column; the seam graph is per state; a cluster is as defined above.  When
//! PlanOptions::trace_dofs is set, non-nodal and C1 sides become abstract DOFs identified through their
//! nodal components, and their receiver copies are tau-only.
//!
//! PROVENANCE (this replaced an earlier decision layer, since deleted): the two builders ran side by
//! side, this one in shadow; it became the applied builder for IC_TRACE, then for IC_STRONG, and for
//! IC_WEAK once IC_WEAK was routed through a decision layer at
//! all.  The corpus then showed the older builder's output was discarded on every one of 822 plan
//! decisions, and that its one remaining privilege -- a deferred re-setup when the redundancy-drop cycle
//! ran under this builder's plan -- produced identical drop sets and errors on all five programs that
//! reach it.  It was then deleted, with its 954-line header.

namespace ocplan_global_detail
{
  struct UnionFind {
    std::map<size_t,size_t> parent;
    size_t find( size_t x ){
      auto it = parent.find( x );
      if( it == parent.end() ){ parent[x] = x; return x; }
      size_t r = x; while( parent[r] != r ) r = parent[r];
      while( parent[x] != r ){ size_t n = parent[x]; parent[x] = r; x = n; }
      return r;
    }
    void unite( size_t a, size_t b ){ a = find(a); b = find(b); if( a != b ) parent[a] = b; }
  };
}

inline bool
OCPlanGlobalBuilder::build( PlanInput const& in, PlanDecisions& dec, PlanReport& rep )
{
  using namespace ocplan_global_detail;
  typedef PlanInput::ClaimKey ClaimKey;
  typedef PlanInput::RowKey   RowKey;
  PlanReport::Shadow& sh = rep.shadow;

  bool const is_weak   = ( in.options.imposition == PlanOptions::Imposition::WEAK );

  // ---- 1. edges (claims after the DROP), DOFs ------------------------------------
  // Occurrence counts are irrelevant here: geometry is per key; every kept key is one edge.
  // Two-pass DROP (pass 2), same rule as the legacy builder: claims in dropped_claims
  // (slot-normalised) are excluded from allocation.
  auto drop_slot_key = [&]( ClaimKey key ) -> ClaimKey {
    auto it = in.trace_claim_block.find( std::make_pair( key.state_id, key.dom_id ) );
    if( it != in.trace_claim_block.end() ) key.block_id = it->second;
    key.eqn_id = std::numeric_limits<size_t>::max();
    if( key.kind == OCPlan::CLAIM_EXACT_ROW ) key.kind = OCPlan::CLAIM_EXACT_C0;
    return key;
  };
  std::map<ClaimKey,size_t> kept;
  for( auto const& kv : in.claims ){
    if( !in.dropped_claims.empty()
     && in.dropped_claims.count( OCPlan::freeze_claim_key( drop_slot_key( kv.first ) ) ) ){ ++sh.n_dropped_by_drop; continue; }
    kept.insert( kv );
  }

  // Edge set: kept nodal C0 claims, plus (shadow only) the claims the legacy corner rule
  // suppressed, so the forest can be checked against that rule.  The latter carry no
  // receiver rows; they never become EXACT/TAU here (see 4.2), only tree/dropped.
  std::vector<PlanInput::DofEdge> all_edges = in.dof_edges;
  size_t const n_own_edges = all_edges.size();
  for( auto const& e : in.corner_suppressed_edges ) all_edges.push_back( e );
  auto is_corner_rule_edge = [&]( size_t ei ){ return ei >= n_own_edges; };

  std::vector<size_t> edge_ids;                         // indices into all_edges
  std::map<ClaimKey,size_t> idx_of_key;
  for( size_t i = 0; i < n_own_edges; ++i ) idx_of_key[ all_edges[i].key ] = i;
  std::vector<ClaimKey> standalone_tau;   // claims outside the DOF model: realised TAU as legacy does
  bool const tdof = in.options.trace_dofs;
  std::vector<size_t> trace_edge_ids;     // rev155a: non-nodal edges, DOFs assigned after the nodal components
  for( auto const& kv : kept ){
    ++sh.n_claims;
    auto it = idx_of_key.find( kv.first );
    if( it == idx_of_key.end() ){ standalone_tau.push_back( kv.first ); continue; }
    PlanInput::DofEdge const& e = all_edges[ it->second ];
    bool const has_supports = !e.lo_support.empty() && !e.hi_support.empty();
    if( kv.first.kind != OCPlan::CLAIM_EXACT_C0 ){ ++sh.n_edges_c1;
      if( tdof && has_supports ) trace_edge_ids.push_back( it->second ); else standalone_tau.push_back( kv.first ); continue; }
    if( !e.nodal ){ ++sh.n_edges_nonnodal;
      if( tdof && has_supports ) trace_edge_ids.push_back( it->second ); else standalone_tau.push_back( kv.first ); continue; }
    ++sh.n_edges_nodal; edge_ids.push_back( it->second );
  }
  for( size_t i = n_own_edges; i < all_edges.size(); ++i ){
    if( all_edges[i].nodal ){ ++sh.n_corner_rule; edge_ids.push_back( i ); }
    else if( tdof && !all_edges[i].lo_support.empty() && !all_edges[i].hi_support.empty() ){
      ++sh.n_corner_rule; ++sh.tdof_corner_nonnodal; trace_edge_ids.push_back( i ); }
  }
  sh.tdof_standalone = standalone_tau.size();
  // The standalone-TAU rule is DELETED as a designed path.  Under trace DOFs every claim
  // is a DofEdge with side supports (standalone=0 on every model measured, rev155a-e).  A claim
  // arriving here with trace DOFs on means a PACKING defect -- keep the safe TAU fallback but say
  // so loudly, ungated.
  if( tdof && !standalone_tau.empty() ){
    auto const& k0 = standalone_tau.front();
    std::cerr << "OCPlanGlobalBuilder ** NOTICE: " << standalone_tau.size()
              << " claim(s) outside the trace-DOF model (no edge or empty side supports) -- "
                 "packing defect, realised TAU as a fallback; first: state=" << k0.state_id
              << " dom=" << k0.dom_id << " iel_lo=" << k0.iel_lo << " face=" << k0.face
              << " kind=" << k0.kind << std::endl;
  }

  // ---- 2. components per state (seam graph G_S) ----------------------------------
  UnionFind comp;                                      // DOFs joined by edges of the SAME state
  std::set<size_t> dofs;
  for( size_t ei : edge_ids ){
    auto const& e = all_edges[ei];
    dofs.insert( e.dof_lo ); dofs.insert( e.dof_hi );
    comp.unite( e.dof_lo, e.dof_hi );
  }
  // rev155a: trace DOFs.  A side that is a single unit column is that column; any other
  // side is an abstract DOF keyed on its support canonicalised through the NODAL
  // components (node-wise equal supports under nodal continuity are the same trace
  // value).  Keys are computed before any trace edge joins the union-find, so the
  // identification uses nodal continuity only -- the claims that must exist for it to
  // hold.  Abstract ids live above every var[] column.
  size_t const abstract_base = std::numeric_limits<size_t>::max() / 2;
  auto is_abstract = [&]( size_t d ){ return d >= abstract_base; };
  if( tdof && !trace_edge_ids.empty() ){
    // Pinned pairs are identities (not edges): they join the identification union-find only.
    for( auto const& e : in.pinned_edges )
      if( e.nodal ){ comp.unite( e.dof_lo, e.dof_hi ); ++sh.tdof_pinned_joins; }
    typedef std::vector<std::pair<size_t,long long>> CanonKey;
    std::map<std::pair<size_t,CanonKey>,size_t> abstract_of;   // (state, canonical support) -> id
    auto side_dof = [&]( size_t state, std::vector<std::pair<size_t,double>> const& sup ) -> size_t {
      if( sup.size() == 1 && std::fabs( std::fabs( sup[0].second ) - 1. ) < 1e-12 ) return sup[0].first;
      ++sh.tdof_sides;
      CanonKey ck; ck.reserve( sup.size() );
      for( auto const& [c,w] : sup ) ck.push_back( { comp.find( c ), (long long)std::llround( w * 1e10 ) } );
      std::sort( ck.begin(), ck.end() );
      auto key = std::make_pair( state, ck );
      auto it = abstract_of.find( key );
      if( it != abstract_of.end() ) return it->second;
      size_t const id = abstract_base + abstract_of.size();
      abstract_of.emplace( key, id ); return id;
    };
    // Identification is INCREMENTAL in causal order (direction rank, then seam index):
    // the identity tau_A == tau_B at seam k needs the nodal supports joined at every
    // upstream node, and at the upstream seam node (k-1,0) that join is created only by
    // the identified trace of seam k-1 (legacy drops that x-claim as MGS-dependent).
    // MEASURED PDE6: nodal-only keys identify 4 per state (seam 1 only); causal order
    // must identify 12 (seams 1..3) -- the structural form of legacy's MGS rank filter.
    std::stable_sort( trace_edge_ids.begin(), trace_edge_ids.end(), [&]( size_t a, size_t b ){
      auto const& ea = all_edges[a]; auto const& eb = all_edges[b];
      if( ea.direction_rank != eb.direction_rank ) return ea.direction_rank < eb.direction_rank;
      if( ea.key.dom_id != eb.key.dom_id ) return ea.key.dom_id < eb.key.dom_id;
      return ea.key.iel_lo < eb.key.iel_lo; } );
    for( size_t g = 0; g < trace_edge_ids.size(); ){
      size_t h = g;
      auto const& e0 = all_edges[ trace_edge_ids[g] ];
      while( h < trace_edge_ids.size()
          && all_edges[ trace_edge_ids[h] ].key.dom_id == e0.key.dom_id
          && all_edges[ trace_edge_ids[h] ].key.iel_lo == e0.key.iel_lo ) ++h;
      for( size_t i = g; i < h; ++i ){          // keys for the whole seam first ...
        auto& e = all_edges[ trace_edge_ids[i] ];
        e.dof_lo = side_dof( e.key.state_id, e.lo_support );
        e.dof_hi = side_dof( e.key.state_id, e.hi_support );
      }
      for( size_t i = g; i < h; ++i ){          // ... then join, so later seams see it
        auto const& e = all_edges[ trace_edge_ids[i] ];
        dofs.insert( e.dof_lo ); dofs.insert( e.dof_hi );
        comp.unite( e.dof_lo, e.dof_hi );
        edge_ids.push_back( trace_edge_ids[i] );
      }
      g = h;
    }
    sh.tdof_abstract = abstract_of.size();
  }
  sh.n_dofs = dofs.size();
  std::map<size_t,std::vector<size_t>> comp_members;   // root -> dofs
  for( size_t d : dofs ) comp_members[ comp.find( d ) ].push_back( d );
  sh.n_components = comp_members.size();
  for( auto const& kv : comp_members ) sh.max_component_dofs = std::max( sh.max_component_dofs, kv.second.size() );

  // ---- 3. rows at DOFs; clusters (DOFs sharing a row copy, plus component links) ----
  std::map<size_t,std::vector<size_t>> rows_of_dof;    // dof -> indices into rows_at_dof
  std::map<RowKey,std::vector<size_t>> dofs_of_row;    // row  -> dofs
  std::map<ClaimKey,std::vector<size_t>> rows_of_claim;// claim -> indices into rows_at_dof
  // rev155a: rows of non-nodal claims arrive with dof == SIZE_MAX; under trace_dofs the
  // side DOF is resolved from the edge, otherwise the row is skipped (pre-155a behaviour).
  // A copy on an abstract side has no coincident node: nothing is redundant there, so it
  // is never replaceable (tau receiver only).
  std::vector<size_t> row_dof( in.rows_at_dof.size(), std::numeric_limits<size_t>::max() );
  std::vector<bool>   row_repl( in.rows_at_dof.size(), false );
  for( size_t i = 0; i < in.rows_at_dof.size(); ++i ){
    auto const& r = in.rows_at_dof[i];
    row_dof[i] = r.dof; row_repl[i] = r.replaceable;
    if( r.dof != std::numeric_limits<size_t>::max() ) continue;
    if( !tdof ) continue;
    auto const& edge = (*in.receiver_edges)[ r.edge_index ];
    auto it = idx_of_key.find( edge.claim_key );
    if( it == idx_of_key.end() ) continue;
    auto const& e = all_edges[ it->second ];
    if( e.lo_support.empty() || e.hi_support.empty() ) continue;
    size_t const d = r.side_hi ? e.dof_hi : e.dof_lo;
    row_dof[i] = d; row_repl[i] = r.replaceable && !is_abstract( d );
    if( is_abstract( d ) ){ ++sh.tdof_rows; if( r.tau_allowed ) ++sh.tdof_rows_tau; }
  }
  for( size_t i = 0; i < in.rows_at_dof.size(); ++i ){
    auto const& r = in.rows_at_dof[i];
    auto const& edge = (*in.receiver_edges)[ r.edge_index ];
    if( !kept.count( edge.claim_key ) ) continue;
    if( row_dof[i] == std::numeric_limits<size_t>::max() ) continue;
    rows_of_dof[ row_dof[i] ].push_back( i );
    dofs_of_row[ r.row ].push_back( row_dof[i] );
    rows_of_claim[ edge.claim_key ].push_back( i );
    ++sh.n_rows_at_dof; if( row_repl[i] ) ++sh.n_rows_replaceable;
  }
  // Replaceable row copies AT A NODE (design §2 refinement, measured on MBC1): a
  // primitive reduced in direction d has no receiver in its own PDE row for its
  // d-claim (the PDE contains the auxiliary, not the derivative), so "receiver rows of
  // the claim" under-reports what may be consumed.  Index every replaceable row copy by
  // its node key (element block, flat); a DOF's candidates are all replaceable copies at
  // the nodes of its receiver rows.
  typedef std::pair<std::vector<std::pair<size_t,size_t>>,size_t> NodeKey;
  std::map<NodeKey,std::vector<size_t>> replaceable_at_node;
  for( size_t i = 0; i < in.rows_at_dof.size(); ++i ){
    auto const& r = in.rows_at_dof[i];
    auto const& edge = (*in.receiver_edges)[ r.edge_index ];
    if( !kept.count( edge.claim_key ) || !row_repl[i] ) continue;
    if( row_dof[i] == std::numeric_limits<size_t>::max() ) continue;
    replaceable_at_node[ { r.row.block_el, r.row.flat } ].push_back( i );
  }
  auto candidates_at_dof = [&]( size_t dof ) -> std::vector<size_t> {
    std::vector<size_t> out; std::set<NodeKey> seen;
    for( size_t ri : rows_of_dof[ dof ] ){
      NodeKey nk{ in.rows_at_dof[ri].row.block_el, in.rows_at_dof[ri].row.flat };
      if( !seen.insert( nk ).second ) continue;
      auto it = replaceable_at_node.find( nk );
      if( it != replaceable_at_node.end() ) out.insert( out.end(), it->second.begin(), it->second.end() );
    }
    return out;
  };
  UnionFind clus;
  for( size_t d : dofs ) clus.unite( d, comp.find( d ) );
  for( auto const& kv : dofs_of_row )
    for( size_t k = 1; k < kv.second.size(); ++k ) clus.unite( kv.second[0], kv.second[k] );
  std::map<size_t,std::vector<size_t>> cluster_dofs;   // cluster root -> dofs
  for( size_t d : dofs ) cluster_dofs[ clus.find( d ) ].push_back( d );
  sh.n_clusters = cluster_dofs.size();
  for( auto const& kv : cluster_dofs ) sh.max_cluster_dofs = std::max( sh.max_cluster_dofs, kv.second.size() );

  // edges per cluster
  std::map<size_t,std::vector<size_t>> cluster_edges;  // cluster root -> edge ids
  for( size_t ei : edge_ids ) cluster_edges[ clus.find( all_edges[ei].dof_lo ) ].push_back( ei );

  // ---- 4. per cluster: forest, realisation, rank check -----------------------------
  // Deterministic edge order (design §3.1): direction rank, seam index, transverse key.
  auto edge_less = [&]( size_t a, size_t b ){
    auto const& ea = all_edges[a]; auto const& eb = all_edges[b];
    int ra = ea.direction_rank, rb = eb.direction_rank;
    if( ra != rb ) return ra < rb;
    if( ea.key.dom_id != eb.key.dom_id ) return ea.key.dom_id < eb.key.dom_id;
    if( ea.key.iel_lo != eb.key.iel_lo ) return ea.key.iel_lo < eb.key.iel_lo;
    if( ea.key.face   != eb.key.face   ) return ea.key.face < eb.key.face;
    return ea.key < eb.key;
  };

  std::set<RowKey> consumed;                            // global: a copy is replaced once
  std::map<std::tuple<size_t,size_t,size_t,size_t>,int> realised_at_locus;   // (state,dom,iel_lo,face) -> realisation

  std::vector<std::vector<std::pair<size_t,size_t>>> cluster_oriented;         // per cluster: (edge id, child)
  size_t cluster_id = 0;
  // IC_TRACE and an IC_WEAK plan with multipliers number their taus as var[] columns after the collocation
  // coefficients; IC_STRONG numbers them inside its Schur block.
  size_t const tau0 = ( in.options.imposition == PlanOptions::Imposition::TRACE
                     || ( is_weak && in.options.weak_multipliers ) ) ? in.trace_var_offset : size_t(0);
  size_t tau = tau0;                                    // tau column ids, legacy numbering convention
  for( auto& kv : cluster_edges ){
    PlanReport::Cluster C; C.id = cluster_id++;
    C.dofs = cluster_dofs[ kv.first ];
    std::vector<size_t> es = kv.second;
    std::sort( es.begin(), es.end(), edge_less );

    // 4.1 forest per state-component inside the cluster (Kruskal in the fixed order)
    UnionFind tree;
    // rev155d (trace_dofs): pinned pairs pre-join the forest -- a claim between two columns
    // pinned by the same ALL-mask IC/BC is implied by the pin rows (FIX1's justification), so
    // an edge closing on them is redundant.  Measured need: PDE5's two BC-pinned corner claims.
    if( tdof ) for( auto const& pe : in.pinned_edges ) if( pe.nodal ) tree.unite( pe.dof_lo, pe.dof_hi );
    std::vector<size_t> tree_edges, dropped;
    for( size_t ei : es ){
      auto const& e = all_edges[ei];
      C.edges.push_back( e.key );
      if( tree.find( e.dof_lo ) != tree.find( e.dof_hi ) ){ tree.unite( e.dof_lo, e.dof_hi ); tree_edges.push_back( ei ); }
      else dropped.push_back( ei );
    }
    for( size_t ei : tree_edges ) C.tree_edges.push_back( all_edges[ei].key );
    for( size_t ei : dropped    ){ C.dropped_edges.push_back( all_edges[ei].key ); C.realisation[ all_edges[ei].key ] = 2; }
    sh.n_tree_edges += tree_edges.size(); sh.n_redundant_edges += dropped.size();
    for( size_t ei : dropped ){
      sh.forest_dropped.push_back( all_edges[ei].key );
      if( is_corner_rule_edge( ei ) ){ ++sh.n_corner_rule_forest_drops; if( !all_edges[ei].nodal ) ++sh.tdof_corner_nonnodal_forest_drops; }
      else ++sh.n_forest_drops_not_corner_rule; }
    // Corner-rule edges that the forest KEPT have no receivers: they cannot be realised
    // here.  In shadow they are left in the tree (counted above as mismatches through
    // n_forest_drops_not_corner_rule, which then names a legacy edge the forest dropped
    // instead); in applied mode they are ordinary claims and this branch is empty.
    { std::vector<size_t> t2; for( size_t ei : tree_edges ) if( !is_corner_rule_edge( ei ) ) t2.push_back( ei ); tree_edges.swap( t2 ); }

    // 4.2 rooted realisation.  Root per component: upstream DOF if causal, else min key.
    // Orientation: BFS from the root over tree edges; child = the far endpoint.
    std::map<size_t,std::vector<std::pair<size_t,size_t>>> adj;   // dof -> (edge id, other dof)
    for( size_t ei : tree_edges ){
      auto const& e = all_edges[ei];
      adj[ e.dof_lo ].push_back( { ei, e.dof_hi } ); adj[ e.dof_hi ].push_back( { ei, e.dof_lo } );
    }
    std::set<size_t> visited;
    std::vector<std::pair<size_t,size_t>> oriented;   // (edge id, child dof)
    for( size_t d0 : C.dofs ){
      if( visited.count( d0 ) ) continue;
      // component of d0 within the tree
      std::vector<size_t> members; { std::vector<size_t> st{ d0 }; std::set<size_t> seen{ d0 };
        while( !st.empty() ){ size_t x = st.back(); st.pop_back(); members.push_back( x );
          for( auto const& [ei, y] : adj[x] ) if( seen.insert( y ).second ) st.push_back( y ); } }
      std::set<size_t> const memset_( members.begin(), members.end() );
      size_t root = std::numeric_limits<size_t>::max();
      bool has_causal = false;
      for( size_t ei : tree_edges ){ auto const& e = all_edges[ei];
        if( e.causal && memset_.count( e.dof_lo ) ){ has_causal = true; break; } }
      if( has_causal ){
        // upstream = the lo side of a causal edge OF THIS COMPONENT that is nobody's hi
        // side within the component.  (rev154f: the search must stay inside the
        // component -- scanning the whole cluster picked a sibling component's dof and
        // the BFS from an already-visited root oriented NOTHING, silently dropping every
        // further causal component of the cluster: PDE36c lost 864 of 1548 claims.)
        std::set<size_t> his;
        for( size_t ei : tree_edges ){ auto const& e = all_edges[ei];
          if( e.causal && memset_.count( e.dof_hi ) ) his.insert( e.dof_hi ); }
        for( size_t ei : tree_edges ){ auto const& e = all_edges[ei];
          if( e.causal && memset_.count( e.dof_lo ) && !his.count( e.dof_lo ) ){ root = e.dof_lo; break; } }
      }
      if( root == std::numeric_limits<size_t>::max() || !memset_.count( root ) )
        root = *std::min_element( members.begin(), members.end() );
      C.root = root;
      std::vector<size_t> q{ root }; visited.insert( root );
      for( size_t qi = 0; qi < q.size(); ++qi ){
        size_t x = q[qi];
        for( auto const& [ei, y] : adj[x] ) if( visited.insert( y ).second ){ oriented.push_back( { ei, y } ); q.push_back( y ); }
      }
    }
    // Design §3.5 is decided in a SECOND PASS (after every cluster's primitives are
    // realised), because the implication spans the auxiliary's whole transverse element
    // line, which crosses clusters.  Cross-direction auxiliary edges are deferred here.
    auto is_cross_aux = [&]( ClaimKey const& k ) -> bool {
      auto it = in.aux_info.find( k.state_id );
      return it != in.aux_info.end() && !it->second.diff_dirs.count( k.dom_id );
    };
    for( auto const& [ei, child] : oriented ){
      auto const& e = all_edges[ei];
      if( !is_weak && is_cross_aux( e.key ) ){ C.realisation[ e.key ] = 4; continue; }  // 4 = PENDING
      bool exact = false;
      // Candidate replaceable copy (design §3.2 + 2026-09-04 refinements):
      //  (A) only PRIMITIVE claims consume rows; an auxiliary's claim never does (the PDE
      //      copies among its receivers belong to its parent).
      //  (B) a copy may be consumed for a d-claim only if its equation has a row at every
      //      node in d (row_full_in_dir); one-sided sets (first-order operators) -> TAU.
      //  (C) candidates: the child's own replaceable receivers first, then the child's
      //      OWNER PDE copies at the same node (reduced states: the PDE contains the
      //      auxiliary, so it is not a receiver of the primitive's claim).
      RowKey const* cand = nullptr;
      bool upwind = false;
      for( size_t ri : rows_of_claim[ e.key ] )
        if( (*in.receiver_edges)[ in.rows_at_dof[ri].edge_index ].sat_type == OCPlan::IC_UPWIND ){ upwind = true; break; }
      bool const is_aux_state = in.aux_info.count( e.key.state_id ) > 0;
      // A claim with NO tau-allowed receiver has exactly one realisation: consume a copy
      // (legacy's geometric-C0 path; measured safe on first-order DRY_C).  A claim that
      // CAN be a tau consumes a copy only if the equation is second order in d with a
      // symmetric collocation set (row_full_in_dir); first-order operators need the
      // continuity as an ADDED condition (GAS_V, MMPDE27's PDE in t).
      bool has_tau_receiver = false;
      for( size_t ri : rows_of_claim[ e.key ] ) if( in.rows_at_dof[ri].tau_allowed ){ has_tau_receiver = true; break; }
      auto full_in_dir = [&]( RowKey const& rk ) -> bool {
        if( !has_tau_receiver ) return true;
        auto it = in.row_full_in_dir.find( std::make_pair( rk.row_id, e.key.dom_id ) );
        return it == in.row_full_in_dir.end() ? true : it->second;
      };
      if( !upwind && !is_aux_state ){
        for( size_t ri : rows_of_dof[ child ] ){
          auto const& r = in.rows_at_dof[ri];
          if( row_repl[ri] && !consumed.count( r.row ) && full_in_dir( r.row ) ){ cand = &r.row; break; }
        }
        if( !cand ){
          auto ito = in.owner_rows.find( e.key.state_id );
          if( ito != in.owner_rows.end() )
            for( size_t ri : candidates_at_dof( child ) ){
              auto const& r = in.rows_at_dof[ri];
              if( ito->second.count( r.row.row_id ) && !consumed.count( r.row ) && full_in_dir( r.row ) ){ cand = &r.row; break; }
            }
        }
      }
      if( is_weak ){
        exact = false;                                  // phase 1a: WEAK ratio rule not yet wired (SAT)
      }
      else if( upwind ){
        exact = false;                                  // first-order seam: continuity ADDS a condition (TAU); replacing a row starves the element
      }
      else{
        exact = ( cand != nullptr );
      }
      if( exact ){ consumed.insert( *cand ); C.consumed_rows.push_back( *cand ); dec.exact_replacement_rows.insert( *cand );
                   C.realisation[ e.key ] = 0; ++sh.n_exact; }
      // IC_WEAK WITH MULTIPLIERS (C2 experiment, CRONOS_WEAK_TAU_RESCUED).  A claim whose coupling had to be
      // SUBSTITUTED has no honest penalty weight; an exact constraint needs none.  So a rescued claim is
      // realised TAU -- provided it has a tau-allowed receiver, since a multiplier with no receiver would be
      // an empty column; without one it stays SAT and is counted.  Everything else under WEAK stays SAT.
      // (NOTES_20260917a: before rev273 this branch assigned TAU ids that no column or row realised.)
      else if( is_weak && in.options.weak_multipliers && in.rescued_claims.count( e.key ) ){
        if( has_tau_receiver ){ C.realisation[ e.key ] = 1; ++sh.n_tau; ++sh.n_weak_tau; }
        else                  { C.realisation[ e.key ] = 3; ++sh.n_sat; ++sh.n_weak_tau_norecv; } }
      else if( is_weak ){ C.realisation[ e.key ] = 3; ++sh.n_sat; }
      else { C.realisation[ e.key ] = 1; ++sh.n_tau; }
      realised_at_locus[ std::make_tuple( e.key.state_id, e.key.dom_id, e.key.iel_lo, e.key.face ) ] = C.realisation[ e.key ];
    }

    cluster_oriented.push_back( oriented );
    rep.clusters.push_back( C );
  }

  // ---- 5. second pass: cross-direction auxiliary edges (design §3.5, with the measured
  //         condition): implied iff the parent's continuity across the same seam is
  //         realised EXACT/TAU BY THIS PLAN at EVERY node of the auxiliary's transverse
  //         element line (the derivative stencil).  A node whose parent claim was
  //         suppressed upstream (unrealised_claims) breaks the implication.
  for( size_t ci = 0; ci < rep.clusters.size(); ++ci ){
    PlanReport::Cluster& C = rep.clusters[ci];
    for( auto const& [ei, child] : cluster_oriented[ci] ){
      auto const& e = all_edges[ei];
      auto itr = C.realisation.find( e.key );
      if( itr == C.realisation.end() || itr->second != 4 ) continue;
      auto ita = in.aux_info.find( e.key.state_id );
      bool implied = ( ita != in.aux_info.end() );
      // rev155e (knob, default = pre-155e): CRONOS_AUX_IMPLIED_EVO=0 -> Sec.3.5 does not apply to
      // an auxiliary's claims in the EVOLUTION direction.  MEASURED: blk0 AUTO TRACE, 60 such
      // drops, DETERMINED but accuracy 1.3055e-5 -> 8.18e-5; PDE3's load-bearing Sec.3.5 drops are
      // all spatial.  Row-receiver and cluster-deficit criteria were both refuted (see notes).
      // DEFAULT 0 (corpus-accepted): Sec.3.5 never applies to evolution-direction aux
      // claims; CRONOS_AUX_IMPLIED_EVO=1 restores the pre-156 behaviour for attribution.
      { static int const kEvo = 0;   // CRONOS_AUX_IMPLIED_EVO (retired 2026-10-07, WORKPLAN 3.B batch 2b)
        if( !kEvo && e.causal ) implied = false; }
      if( implied ){
        // transverse element line: 2D only (one differentiated direction); higher
        // dimensions fall back to "not implied".
        size_t nn = 0;
        if( ita->second.diff_dirs.size() == 1 ){
          auto itn = in.dom_n_node.find( *ita->second.diff_dirs.begin() );
          if( itn != in.dom_n_node.end() ) nn = itn->second;
        }
        if( !nn ) implied = false;
        else{
          size_t const e0 = ( e.key.face / nn ) * nn;
          for( size_t f = e0; f < e0 + nn && implied; ++f ){
            auto itp = realised_at_locus.find( std::make_tuple( ita->second.parent, e.key.dom_id, e.key.iel_lo, f ) );
            bool ok = ( itp != realised_at_locus.end() && ( itp->second == 0 || itp->second == 1 ) );
            // A parent claim whose exact row the corner rank filter REJECTED is not thereby
            // unenforced: the rejected row is a linear combination of rows the plan does impose,
            // so the parent's continuity still holds at that locus.  Treating the rejection as a
            // MISSING claim is a false negative, and it costs the child its aux-implied
            // realisation.  OFF by default; CRONOS_AUX_IMPLIED_RANK=1 enables the reading.
            // MEASURED (PDE3, with the corner dedup retired): 68 claims move AUX-IMPLIED -> TAU
            // and n grows 6,038 -> 6,100 when the rejection IS read as missing.
            if( !ok ){
              static int const kRankImplied = 0;   // CRONOS_AUX_IMPLIED_RANK (retired 2026-10-07, WORKPLAN 3.B batch 2b)
              if( kRankImplied && !in.rank_implied_claims.empty() ){
                ClaimKey pk = e.key;
                pk.state_id = ita->second.parent;
                pk.face     = f;
                pk.eqn_id   = std::numeric_limits<size_t>::max();
                pk.kind     = OCPlan::CLAIM_EXACT_C0;
                auto blk = in.trace_claim_block.find( std::make_pair( pk.state_id, pk.dom_id ) );
                if( blk != in.trace_claim_block.end() ) pk.block_id = blk->second;
                ok = in.rank_implied_claims.count( pk ) > 0;
              }
            }
            if( !ok ){
              // the parent may have no claim at f because it is pinned there (unrealised)
              ClaimKey pk = e.key; pk.state_id = ita->second.parent; pk.face = f;
              auto blk = in.trace_claim_block.find( std::make_pair( pk.state_id, pk.dom_id ) );
              if( blk != in.trace_claim_block.end() ) pk.block_id = blk->second;
              pk.eqn_id = std::numeric_limits<size_t>::max(); pk.kind = OCPlan::CLAIM_EXACT_C0;
              (void)pk;   // any missing parent claim at f breaks the implication
            }
            implied = ok;
          }
        }
      }
      if( implied ){ C.realisation[ e.key ] = 2; ++sh.n_aux_implied; }
      else { C.realisation[ e.key ] = 1; ++sh.n_tau; }
      realised_at_locus[ std::make_tuple( e.key.state_id, e.key.dom_id, e.key.iel_lo, e.key.face ) ] = C.realisation[ e.key ];
    }
  }

  // ---- 6. rank check per cluster (TAU edges only; receivers minus consumed) ----
  for( size_t ci = 0; ci < rep.clusters.size(); ++ci ){
    PlanReport::Cluster& C = rep.clusters[ci];
    auto const& oriented = cluster_oriented[ci];
    // 4.3 rank check of the cluster's tau block (TAU edges only; receivers minus consumed),
    //     with the structural fallback for collinear pairs (design §3.3 / §3.5): a
    //     same-direction auxiliary claim whose receiver column equals its parent's carries
    //     no independent condition at that DOF (the parent is pinned there by a
    //     non-replaceable row on both sides); it is dropped and the rank re-checked.
    if( !is_weak || in.options.weak_multipliers ){   // a weak plan has a tau block only with multipliers
      for( int pass = 0; pass < 4; ++pass ){
        std::vector<size_t> tau_edges;
        for( auto const& [ei, child] : oriented ) if( C.realisation[ all_edges[ei].key ] == 1 ) tau_edges.push_back( ei );
        if( tau_edges.empty() ) break;
        std::map<RowKey,size_t> row_index; std::vector<std::vector<std::pair<size_t,double>>> cols;
        C.tau_receivers.clear();
        for( size_t ei : tau_edges ){
          std::vector<std::pair<size_t,double>> col;
          for( size_t ri : rows_of_claim[ all_edges[ei].key ] ){
            auto const& r = in.rows_at_dof[ri];
            if( !r.tau_allowed || consumed.count( r.row ) ) continue;
            C.tau_receivers.push_back( { all_edges[ei].key, { r.row, r.coeff } } );
            auto it = row_index.find( r.row ); size_t ridx;
            if( it == row_index.end() ){ ridx = row_index.size(); row_index[ r.row ] = ridx; } else ridx = it->second;
            col.push_back( { ridx, r.coeff } );
          }
          cols.push_back( col );
        }
        arma::mat B( row_index.size() ? row_index.size() : 1, tau_edges.size(), arma::fill::zeros );
        for( size_t j = 0; j < cols.size(); ++j ) for( auto const& [ri, v] : cols[j] ) B( ri, j ) += v;
        C.tau_block_cols = tau_edges.size();
        C.tau_block_rank = row_index.empty() ? 0 : (size_t)arma::rank( B );
        if( C.tau_block_rank >= C.tau_block_cols ){ C.unresolved = false; break; }
        // (STRONG explicit selection moved to the full-tau component scan below, rev154e:
        //  per-cluster scans cannot see dependencies involving standalone taus -- PDE7.)
        // Collinear pairs are NOT dropped here (measured: each carries a genuine condition,
        // MBC1 defR 4->16 when dropped).  They are left to the faithful trace projection
        // in materialisation (the gauge fix legacy uses).  Reported as unresolved.
        bool dropped_any = false;
        if( !dropped_any ){
          ++sh.n_clusters_rank_deficient; sh.n_rank_deficit += C.tau_block_cols - C.tau_block_rank;
          C.unresolved = true; ++rep.n_unresolved;
          break;
        }
      }
    }
  }

  // Claims that get a constraint row: every edge realised EXACT or TAU (one occurrence
  // each: geometry is per key).  Redundant / implied / corner-rule edges get none.
  // Tau ids are assigned here, after every cluster is settled (legacy numbering).
  dec.claims_kept.clear();
  // Ordered tau list: standalone claims first (non-nodal trace functionals, C1), then
  // cluster TAU edges in cluster order.  (rev154a dropped standalone claims: PDE5/6/8
  // corpus regressions, 2026-09-04.)
  std::vector<ClaimKey> tau_order;
  for( auto const& k : standalone_tau ){
    // Under IC_WEAK with multipliers only the cluster edges above are imposed exactly; a standalone claim
    // stays a penalty.  Without multipliers the historical numbering is kept unchanged (a weak plan then
    // numbers a tau id no column realises -- OCFESLV reports it).
    if( is_weak && in.options.weak_multipliers ){ ++sh.n_weak_standalone_sat; continue; }
    dec.claims_kept[ k ] = 1; tau_order.push_back( k ); ++sh.n_tau; }
  for( auto const& C : rep.clusters )
    for( auto const& kv : C.realisation ){
      if( kv.second == 0 || kv.second == 1 ) dec.claims_kept[ kv.first ] = 1;
      if( kv.second == 1 ) tau_order.push_back( kv.first );
    }

  // IC_STRONG structural keep-explicit (rev154e): the Schur elimination needs the
  // ELIMINATED tau columns' receiver block at full column rank.  Scan the FULL tau set
  // -- cluster and standalone taus together, since dependencies cross that boundary
  // (PDE7: t-trace x x-claim corners) -- per connected component of the row-sharing
  // graph; columns that do not raise their component's rank stay EXPLICIT var[] columns.
  std::set<ClaimKey> explicit_keys;
  if( in.options.imposition == PlanOptions::Imposition::STRONG && !tau_order.empty() ){
    std::set<ClaimKey> tau_set( tau_order.begin(), tau_order.end() );
    std::vector<ClaimKey> const& scan_order = tau_order;
    std::map<ClaimKey,std::vector<std::pair<RowKey,double>>> recv;
    for( size_t i = 0; i < in.receiver_edges->size(); ++i ){
      auto const& edge = (*in.receiver_edges)[i];
      if( !edge.trace_tau_allowed || !tau_set.count( edge.claim_key ) ) continue;
      RowKey rk{ edge.row_id, edge.receiver_block_el, edge.receiver_emit_flat };
      if( consumed.count( rk ) ) continue;
      double const cf = ( i < in.edge_trace_coeff.size() && in.edge_trace_coeff[i] != 0. )
                      ? in.edge_trace_coeff[i] : edge.coupling * edge.orientation;
      recv[ edge.claim_key ].push_back( { rk, cf } );
    }
    // components over shared rows
    std::vector<size_t> par( scan_order.size() );
    for( size_t i = 0; i < par.size(); ++i ) par[i] = i;
    std::function<size_t(size_t)> find = [&]( size_t a ){ while( par[a] != a ){ par[a] = par[par[a]]; a = par[a]; } return a; };
    { std::map<RowKey,size_t> firstcol;
      for( size_t i = 0; i < scan_order.size(); ++i )
        for( auto const& [rk,c] : recv[ scan_order[i] ] ){
          auto it = firstcol.find( rk );
          if( it == firstcol.end() ) firstcol.emplace( rk, i );
          else { size_t a = find( i ), b = find( it->second ); if( a != b ) par[a] = b; }
        }
    }
    std::map<size_t,std::vector<size_t>> comps;
    for( size_t i = 0; i < scan_order.size(); ++i ) comps[ find( i ) ].push_back( i );
    bool forced_bad = false;
    for( auto const& [root, cols] : comps ){
      if( cols.size() < 2 ) continue;
      std::map<RowKey,size_t> ridx;
      for( size_t i : cols ) for( auto const& [rk,c] : recv[ scan_order[i] ] ) ridx.emplace( rk, ridx.size() );
      arma::mat B( ridx.size(), cols.size(), arma::fill::zeros );
      for( size_t j = 0; j < cols.size(); ++j )
        for( auto const& [rk,c] : recv[ scan_order[ cols[j] ] ] ) B( ridx[rk], j ) += c;
      if( cols.size() > 300 && (size_t)arma::rank( B ) == cols.size() ) continue;   // big and full: done
      arma::mat Bacc( B.n_rows, 0 ); size_t racc = 0, last_raising = cols.size();
      for( size_t j = 0; j < cols.size(); ++j ){
        arma::mat Bt = arma::join_rows( Bacc, B.col( j ) );
        size_t const rt = (size_t)arma::rank( Bt );
        if( rt > racc ){ Bacc = Bt; racc = rt; last_raising = j; continue; }
        ClaimKey const& dk = scan_order[ cols[j] ];
        if( in.options.force_bad_keep_explicit && !forced_bad && last_raising < cols.size() ){
          explicit_keys.insert( scan_order[ cols[ last_raising ] ] ); ++sh.n_explicit; forced_bad = true;
          continue;
        }
        // DEPENDENT column: division of labour with the redundancy machinery (rev154f).
        // Drop-ELIGIBLE block and not W-protected: leave INTERNAL -- the detection pass
        // sees it among the tau columns, W-validates the drop, and the rebuild excludes
        // it (PDE20f: 20 implied claims, def0=0).  Ineligible or protected: EXPLICIT
        // (MBC pairs, PDE7 corners: the W test proved dropping loses conditions).
        auto ite = in.block_drop_eligible.find( dk.block_id );
        bool const eligible = ( ite != in.block_drop_eligible.end() && ite->second )
                           && !in.protected_claims.count( OCPlan::freeze_claim_key( dk ) );
        if( !eligible ){ explicit_keys.insert( dk ); ++sh.n_explicit; }
      }
    }
  }

  // rev281: promoted claims are explicit whatever the per-cluster rank block decided.
  for( auto const& k : tau_order )
    if( in.promote_explicit.count( OCPlan::freeze_claim_key( k ) ) && explicit_keys.insert( k ).second ){
      ++sh.n_explicit; ++sh.n_promoted; }
  size_t n_explicit = 0;
  for( auto const& k : tau_order ){
    if( explicit_keys.count( k ) ){
      size_t const ecol = in.trace_var_offset + n_explicit++;
      dec.tau_claims[ OCPlan::freeze_claim_key( k ) ].push_back( ecol );
      dec.explicit_tau_claims.insert( OCPlan::freeze_claim_key( k ) );
      dec.keep_explicit_reason[ OCPlan::freeze_claim_key( k ) ] = 1;
    }
    else dec.tau_claims[ OCPlan::freeze_claim_key( k ) ].push_back( tau++ );
  }
  dec.n_trace_constraint = dec.claims_kept.size();

  // Receiver edges of claims that did not get a tau (EXACT, REDUNDANT, implied, dropped,
  // non-nodal) must be tau-inactive: the weak-SAT/trace materialisation requires a lambda
  // slot behind every retained tau-allowed edge (same convention as legacy).
  for( size_t i = 0; i < in.receiver_edges->size(); ++i ){
    auto const& edge = (*in.receiver_edges)[i];
    if( !edge.trace_tau_allowed ) continue;
    if( !dec.tau_claims.count( OCPlan::freeze_claim_key( edge.claim_key ) ) )
      dec.trace_tau_allowed_updates.push_back( { i, false } );
  }
  dec.n_trace_var          = ( tau - tau0 ) + n_explicit;
  dec.n_trace_var_explicit = n_explicit;
  dec.n_coll_eqn           = in.n_coll_eqn;
  dec.n_coll_var           = in.n_coll_var;
  if( in.options.imposition == PlanOptions::Imposition::TRACE ){
    dec.n_coll_var  = in.trace_var_offset + dec.n_trace_var;
    dec.n_coll_eqn += dec.n_trace_var;
  }
  else if( in.options.imposition == PlanOptions::Imposition::STRONG && n_explicit ){
    dec.n_coll_var  = in.trace_var_offset + n_explicit;
    dec.n_coll_eqn += n_explicit;
  }
  else if( is_weak && in.options.weak_multipliers && dec.n_trace_var ){
    // IC_WEAK with multipliers: the taus are var[] columns and each has its trace row, as in IC_TRACE.
    dec.n_coll_var  = in.trace_var_offset + dec.n_trace_var;
    dec.n_coll_eqn += dec.n_trace_var;
  }
  return true;
}


//============================================================================================================
// PLAN REPORTING -- the shadow and realisation lines
//   (was ocplan_report.hpp)
//============================================================================================================

//! @brief Printing for PlanReport.  THE BUILDER NEVER PRINTS -- it is a pure decision function, and this is
//! the only place its results become text.  That separation is what lets the plan be recomputed and compared
//! without side effects, and it is why the builder takes no environment reads.
//!
//! The lines it emits are the ones the sweeps parse, so their SHAPE is part of the interface:
//!   [plan/shadow] ... realisation: exact= tau= (explicit ) sat=   -- how the claims were realised
//!   [plan/shadow] ... tau-block: deficient clusters= deficit=      -- the claim-level rank prediction
//!   [plan/shadow] cluster#N ... edge state= dom= ... -> REALISATION receivers: ...
//! Changing a field name or its order breaks the corpus comparison scripts, not just a human reader.

struct OCPlanReport
{
  //! Dump every cluster touching a given state id (all edges with realisation, consumed rows).
  static void dump_state( std::ostream& os, PlanReport const& rep, PlanInput const& in, size_t sid )
  {
    auto nm = [&]( size_t s ){ auto it = in.state_name.find( s ); return it == in.state_name.end() ? std::string("?") : it->second; };
    for( auto const& C : rep.clusters ){
      bool hit = false; for( auto const& e : C.edges ) if( e.state_id == sid ) hit = true;
      if( !hit ) continue;
      os << "  [plan/dump] cluster#" << C.id << " root=" << C.root << " dofs:";
      for( size_t d : C.dofs ) os << " " << d;
      os << "  tau-block " << C.tau_block_rank << "/" << C.tau_block_cols << std::endl;
      for( auto const& e : C.edges ){
        auto it = C.realisation.find( e );
        os << "      " << nm( e.state_id ) << " dom=" << e.dom_id << " iel_lo=" << e.iel_lo << " face=" << e.face
           << " -> " << ( it == C.realisation.end() ? "?" : ( it->second == 0 ? "EXACT" : it->second == 1 ? "TAU" : it->second == 2 ? "REDUNDANT" : "SAT" ) ) << std::endl;
      }
      if( !C.consumed_rows.empty() ){
        for( auto const& r : C.consumed_rows ){ os << "      consumed row_id=" << r.row_id << " flat=" << r.flat << " blk="; for( auto const& b : r.block_el ) os << "(" << b.first << ":" << b.second << ")"; os << std::endl; }
      }
    }
  }

  //! trace-DOF census (a trace_dofs build) against the default GLOBAL build.
  static void print_tdof( std::ostream& os, PlanReport const& rep, size_t tau_common, size_t tau_only_tdof, size_t tau_only_default,
                          size_t exact_common, size_t exact_only_tdof, size_t exact_only_default )
  {
    PlanReport::Shadow const& s = rep.shadow;
    os << "  [plan/tdof] nonnodal=" << s.n_edges_nonnodal << " c1=" << s.n_edges_c1 << " standalone=" << s.tdof_standalone
       << " | abstract sides=" << s.tdof_sides << " -> dofs=" << s.tdof_abstract << " (identified " << ( s.tdof_sides > s.tdof_abstract ? s.tdof_sides - s.tdof_abstract : 0 ) << ")"
       << " rows@abstract=" << s.tdof_rows << " (tau " << s.tdof_rows_tau << ") pinned-joins=" << s.tdof_pinned_joins
       << "\n  [plan/tdof] dofs=" << s.n_dofs << " components=" << s.n_components << " (max " << s.max_component_dofs << ")"
       << " clusters=" << s.n_clusters << " (max " << s.max_cluster_dofs << ")"
       << " | forest: tree=" << s.n_tree_edges << " redundant=" << s.n_redundant_edges
       << " | corner-rule=" << s.n_corner_rule << " (nonnodal " << s.tdof_corner_nonnodal << ") forest drops " << s.n_corner_rule_forest_drops
       << " (nonnodal " << s.tdof_corner_nonnodal_forest_drops << "); drops " << s.n_forest_drops_not_corner_rule << " others"
       << "\n  [plan/tdof] realisation: exact=" << s.n_exact << " tau=" << s.n_tau << " aux-implied=" << s.n_aux_implied
       << " | tau-block deficient=" << s.n_clusters_rank_deficient << " deficit=" << s.n_rank_deficit
       << " | vs default GLOBAL: tau common=" << tau_common << " +" << tau_only_tdof << " -" << tau_only_default
       << " exact common=" << exact_common << " +" << exact_only_tdof << " -" << exact_only_default
       << ( tau_only_tdof == 0 && tau_only_default == 0 && exact_only_tdof == 0 && exact_only_default == 0 ? "  [SAME PLAN]" : "  [DIFFERS]" )
       << std::endl;
  }

  static void print_shadow( std::ostream& os, PlanReport const& rep )
  {
    PlanReport::Shadow const& s = rep.shadow;
    os << "  [plan/shadow] claims=" << s.n_claims
       << " edges: nodal=" << s.n_edges_nodal << " nonnodal=" << s.n_edges_nonnodal << " c1=" << s.n_edges_c1
       << " | dofs=" << s.n_dofs << " components=" << s.n_components << " (max " << s.max_component_dofs << ")"
       << " clusters=" << s.n_clusters << " (max " << s.max_cluster_dofs << ")"
       << " | rows@dof=" << s.n_rows_at_dof << " replaceable=" << s.n_rows_replaceable
       << "\n  [plan/shadow] drop: by-DROP=" << s.n_dropped_by_drop << " corner-rule=" << s.n_corner_rule
       << " (forest drops " << s.n_corner_rule_forest_drops << " of them; drops " << s.n_forest_drops_not_corner_rule << " others)"
       << " aux-implied=" << s.n_aux_implied
       << "\n  [plan/shadow] forest: tree=" << s.n_tree_edges << " redundant=" << s.n_redundant_edges
       << " | realisation: exact=" << s.n_exact << " tau=" << s.n_tau << " (explicit " << s.n_explicit << ") sat=" << s.n_sat
       << ( ( s.n_weak_tau || s.n_weak_tau_norecv || s.n_weak_standalone_sat )
            ? " [weak+tau: exact=" + std::to_string( s.n_weak_tau ) + " no-tau-receiver=" + std::to_string( s.n_weak_tau_norecv )
              + " standalone-kept-sat=" + std::to_string( s.n_weak_standalone_sat ) + "]" : std::string() )
       << ( s.n_promoted ? " [promoted-explicit=" + std::to_string( s.n_promoted ) + "]" : std::string() )
       << " | tau-block: deficient clusters=" << s.n_clusters_rank_deficient << " deficit=" << s.n_rank_deficit
       // rev170: name the consequence, so the number is readable where it is produced.  MEASURED
       // corpus-wide: this deficit EQUALS the IC_STRONG k' on every driver that carries one, and
       // 0 exactly where the claim drop removed the claims.  A deficit is a count of continuity
       // claims the plan IMPLIES (an algebraic state whose defining relation is pointwise in
       // states inherits continuity from its arguments); each is realised as a tau multiplier
       // that enforces nothing.  Expected, not defective -- 15 corpus drivers carry one and pass.
       // If this ever DIFFERS from the IC_STRONG k', something else is deficient and that is the
       // case worth catching.
       << ( s.n_rank_deficit ? " (implied claims; expect IC_STRONG k'=" : "" )
       << ( s.n_rank_deficit ? std::to_string( s.n_rank_deficit ) + ")" : std::string() )
       << "\n  [plan/shadow] vs legacy: legacy exact=" << s.legacy_exact << " tau=" << s.legacy_tau
       << " ; exact rows in common=" << s.exact_rows_common << std::endl;
    size_t shown = 0;
    for( auto const& C : rep.clusters ){
      if( !C.unresolved ) continue;
      if( ++shown > 8 ){ os << "  [plan/shadow] ... more deficient clusters" << std::endl; break; }
      os << "  [plan/shadow] cluster#" << C.id << " dofs=" << C.dofs.size() << " edges=" << C.edges.size()
         << " tree=" << C.tree_edges.size() << " dropped=" << C.dropped_edges.size()
         << " consumed=" << C.consumed_rows.size()
         << " tau-block " << C.tau_block_rank << "/" << C.tau_block_cols << std::endl;
      if( shown == 1 ){
        for( auto const& e : C.edges ){
          auto it = C.realisation.find( e );
          os << "      edge state=" << e.state_id << " dom=" << e.dom_id << " iel_lo=" << e.iel_lo
             << " face=" << e.face << " kind=" << e.kind << " -> "
             << ( it == C.realisation.end() ? "?" : ( it->second == 0 ? "EXACT" : it->second == 1 ? "TAU" : it->second == 2 ? "REDUNDANT" : "SAT" ) )
             << "  receivers:";
          for( auto const& rr : C.tau_receivers ) if( !( rr.first < e ) && !( e < rr.first ) )
            os << "  [row_id=" << rr.second.first.row_id << " flat=" << rr.second.first.flat << " c=" << rr.second.second << "]";
          os << std::endl;
        }
      }
    }
  }
};

} // namespace mc

#if defined(_WIN32)
# pragma pop_macro("INTERFACE")
#endif
#endif  // CRONOS__OCPLAN_HPP
