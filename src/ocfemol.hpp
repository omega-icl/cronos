// Copyright (C) Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the EPL-2.0 with GPL-2.0-or-later as a Secondary License; see the LICENSE file.
// SPDX-License-Identifier: EPL-2.0 OR GPL-2.0-or-later

#ifndef CRONOS__OCFEMOL_HPP
#define CRONOS__OCFEMOL_HPP

#include "ocfebase.hpp"

namespace mc
{

//! @brief Method of lines on orthogonal collocation: a model declared on (t, z) is collocated in z only and
////////////////////////////////////////////////////////////////////////
//! mc::OCFEMOL reads a model declared on an evolution direction t and spatial directions z as
//! M(x) d_t x = R(t, z, x, d_z x, u) (FFModel::peel_evolution), collocates the spatial directions through the
//! shared base OCFEBASE -- nodes, interface plan, imposition -- and hands the semi-discrete DAE to an
//! integrator: Mc d_t xc = Rc(t, xc), with the INITIAL rows as the Ini stream.  Phase 3 (setup): the peel, the
//! working model without t, the collocated layout and plan, the per-node DOF audit.  Phases 4-5 add the
//! residual and initial-condition operations and the DAESLV hand-over.  MoL workplan 2026-10-31, section 2.1.
////////////////////////////////////////////////////////////////////////
class OCFEMOL
: public OCFEBASE
{
public:

  //! @brief OCFEMOL's options: OCFEBASE::Options -- INTERFACE, OUTPUT, AUDIT, DETERMINACY, FATAL, the groups the
  //! collocation base reads -- plus IMPOSITION, the method-of-lines rules.  None of OCFESLV's SOLVE, SENSITIVITY
  //! or HOMOTOPY: the integrator owns those decisions here.
  struct Options
  : public OCFEBASE::Options
  {
    //! @brief What a Neumann/Robin (flux) boundary node is in the semi-discrete DAE.  After order reduction the node
    //! holds the field u and its spatial derivative D, and needs two rows.
    enum NeumannNode
    {
      NEUMANN_ALGEBRAIC = 0,   //!< OCFESLV's rule: the LINK row and the boundary row; u at the node is ALGEBRAIC (no d_t row,
                               //!< no initial condition).  Leg 1 <-> leg 3 of the gate are the SAME discrete system.
      NEUMANN_DIFFERENTIAL     //!< the PDE row and the boundary row; u at the node stays DIFFERENTIAL with its initial
                               //!< condition.  More accurate on the rod (2.3e-4 vs 7.0e-4 at 3x5); differs from OCFESLV
                               //!< at that node by spatial error.
    };
    struct t_Imposition
    {
      NeumannNode NEUMANN_NODE;   //!< default NEUMANN_ALGEBRAIC (decision 31 Oct 2026: the exactness gate first)
      bool REDUCE_HIGH_INDEX;     //!< a model of structural index > 1 along t: run FFModel's index reduction (dummy
                                  //!< derivatives, hidden constraints) on the full (t, z) model before t is taken out,
                                  //!< then collocate the reduced, index-1 system.  Default false: refuse instead (phase 8
                                  //!< of the MoL plan, offered here as an experiment).
    } IMPOSITION;

    Options
      ()
      { reset(); }

    Options& operator=
      ( Options const& opt )
      {
        OCFEBASE::Options::operator=( opt );
        IMPOSITION = opt.IMPOSITION;
        return *this;
      }

    void apply_environment_defaults
      ()
      { OCFEBASE::Options::apply_environment_defaults(); }

    void reset
      ()
      {
        OCFEBASE::Options::reset();
        IMPOSITION.NEUMANN_NODE = NEUMANN_ALGEBRAIC;
        IMPOSITION.REDUCE_HIGH_INDEX = false;
        apply_environment_defaults();
      }
  };

  Options options;

  //! @brief Constructor: OCFEBASE reads the options through a reference to this object's own.
  OCFEMOL
    ( FFGraph* dag = nullptr )
    : OCFEBASE( options )
    {
      _evoCollocated = false;   // t is the integrator's: no causal time discretisation, no evolution-input probe
      set( dag );
    }

  virtual ~OCFEMOL() = default;

  //! @brief Reset and attach a DAG (FFModel::set), with this solver's state cleared first.
  void set
    ( FFGraph* dag )
    override
    {
      _reset();
      FFModel::set( dag );
    }

  //! @brief The semi-discrete system's shape after setup: what the integrator will see.
  struct t_SemiDiscrete
  {
    size_t n_states    = 0;   //!< collocated unknowns xc: field values at nodes, plus the explicit tau multipliers (n_tau of them)
    size_t n_tau       = 0;   //!< interface multipliers minted by the plan (algebraic)
    size_t n_rows      = 0;   //!< rows of Rc after imposition
    size_t n_ini       = 0;   //!< INITIAL rows of the model (the Ini stream), before collocation
    size_t n_fields    = 0;   //!< fields (states of the model, auxiliaries included)
    size_t n_diff_fields = 0; //!< fields with a d_t somewhere: the differential ones
    size_t n_diff_cols = 0;   //!< differential UNKNOWNS: columns of xc with a mass entry (the integrator's id vector)
    size_t n_alg_rows  = 0;   //!< rows of Rc without a mass entry (plan-owned rows included)
    size_t n_ini_stream = 0;  //!< INITIAL (row, node) pairs at differential nodes: the Ini stream
    size_t n_ini_dropped = 0; //!< INITIAL (row, node) pairs at algebraic nodes: dropped, checked for consistency at t0
    double t_lo = 0., t_hi = 0.;  //!< the evolution span, from the model's t domain
    //! @brief FFModel's classification of the FULL (t, z) model, per block -- the same classifier, on the same layout,
    //! as OCFESLV's: PARABOLIC for the rod and BZ, DIFFERENTIAL_ALGEBRAIC or DIFFERENTIAL_ORDINARY for a lumped DAE or ODE.  The
    //! base's classification that follows in the report is of the spatial operator at fixed t (ELLIPTIC for the heat
    //! equation): what the plan collocates, not what the model is.
    std::map<int,EqnType> model_class;
    bool   ran = false;
  };

  //! @brief The peel this setup took along t, on the working (order-reduced) model.
  t_Peel const& peel() const { return _peel; }
  //! @brief The semi-discrete shape (valid after setup).
  t_SemiDiscrete const& semi_discrete() const { return _shape; }
  //! @brief The INITIAL rows set aside as the Ini stream (working model, order-reduced), by equation index.
  std::vector<t_Eqn> const& ini_rows() const { return _iniEqn; }

  //! @brief One entry of Mc, the mass matrix of the semi-discrete DAE  Mc(xc) d_t xc = Rc(t, xc): row @p row of Rc
  //! carries the derivative of unknown @p col with coefficient @p coeff -- a DAG expression in the model's variables,
  //! to be evaluated at the row's node (pointwise: Mc is block-diagonal over nodes).  Phase 3 records the structure
  //! and the expressions; the hand-over (phase 5) evaluates them.
  struct t_MassEntry { size_t row; size_t col; FFVar coeff; };
  //! @brief Mc's entries, in row order (valid after setup).
  std::vector<t_MassEntry> const& mass() const { return _mass; }
  //! @brief Whether unknown @p col is differential (some row carries its d_t): the integrator's id vector.
  bool is_differential( size_t col ) const { return col < _colDiff.size() && _colDiff[col]; }
  //! @brief Whether row @p row of Rc is algebraic (no mass entry).
  bool is_algebraic_row( size_t row ) const { return row < _rowAlg.size() && _rowAlg[row]; }

  //! @brief Rc(t, xc): the collocated residual, one evaluation of the interior expressions as vectors over the
  //! nodes through the OCVar environment, with t as the scalar input and the model's constants set.  @p xc has
  //! semi_discrete().n_states entries, @p rc semi_discrete().n_rows.  The mass term Mc(xc) d_t xc is not included:
  //! the integrator's residual is Mc xdot - Rc, with Mc from mass().  Phase 4 of the MoL plan.
  bool eval_residual
    ( double t, double const* xc, double* rc );

  //! @brief The Jacobian of Rc(t, xc) with respect to xc, sparse, in CSR: row r's entries are columns
  //! jac_col[jac_ptr[r] .. jac_ptr[r+1]) with values jac_val[...]; and the derivative with respect to t in @p drdt
  //! (may be null).  One compressed forward-mode pass per colour of the base's Jacobian pattern (built at setup from
  //! the dependence pass), seeded as OCFESLV seeds its own; the input columns (t, lumped inputs) are kept out of the
  //! CSR and t's reported separately.  Also fills @p rc (may be null).  Phase 4 of the MoL plan.
  bool deriv_residual
    ( double t, double const* xc, double* rc, std::vector<size_t>& jac_ptr, std::vector<size_t>& jac_col,
      std::vector<double>& jac_val, double* drdt );

  //! @brief One entry of the Ini stream: INITIAL row @p row evaluated at a node -- the fields it reads take their values
  //! from xc at @p cols (field id -> column), the spatial coordinates are @p coord (direction -> value), and the
  //! unknown it fixes is @p target (the column of the field this entry is the initial condition of).
  struct t_IniEntry
  {
    size_t row = 0;                                     //!< index into ini_rows()
    size_t target = 0;                                  //!< the differential unknown this condition fixes
    std::vector< std::pair<FFVar,size_t> > cols;        //!< (field, column of its value at the node)
    std::vector< std::pair<FFVar,double> > coord;       //!< (direction, coordinate of the node)
    bool dropped = false;                               //!< at an algebraic node: not in the stream, checked for consistency at t0
  };
  //! @brief The Ini stream (entries with dropped == false) and the dropped pairs (dropped == true), after setup.
  std::vector<t_IniEntry> const& ini_stream() const { return _iniStream; }

  //! @brief The INITIAL rows at the stream's nodes: @p ini has ini_stream().size() entries in that order (the dropped
  //! ones included, so the consistency check at t0 reads the same array).  Each is one small DAG evaluation.
  bool eval_initial
    ( double t, double const* xc, double* ini );

  //! @brief The derivative of eval_initial with respect to xc, in CSR over the entries (a few columns each).
  bool deriv_initial
    ( double t, double const* xc, double* ini, std::vector<size_t>& ptr, std::vector<size_t>& col, std::vector<double>& val );

  //! @brief The Jacobian pattern over xc alone (CSR row pointers and columns), valid after setup; the values come
  //! from deriv_residual, in the same order.
  void jacobian_pattern
    ( std::vector<size_t>& jac_ptr, std::vector<size_t>& jac_col )
    const;

protected:

  t_Peel                _peel;      //!< the peel along t (working model)
  std::vector< std::vector< std::pair<size_t,FFVar> > > _massByEqn;   //!< per row of the rewritten _mEqn: (field, coefficient)
  std::vector<t_MassEntry> _mass;    //!< Mc as data
  std::vector<t_IniEntry>  _iniStream;   //!< the Ini stream and the dropped pairs
  std::vector<bool>     _colDiff;    //!< per unknown of xc
  std::vector<bool>     _rowAlg;     //!< per row of Rc
  std::vector<t_Eqn>    _iniEqn;    //!< the INITIAL rows, kept aside from the collocated system
  double                _evoLo = 0., _evoHi = 0.;   //!< the model's t span, taken out of the working _mDom
  t_SemiDiscrete        _shape;

  //! @brief The single transfer point: the CRONOS_DISPLAY_LEVEL floor on this solver's options, then the model
  //! layer's slice of them -- FFModel::options is what the model layer reads; `options` here is the richer object.
  void _apply_options
    ()
    override
    {
      FFModel::Options::_force_display_from_env( options.DISPLAY_LEVEL );
      FFModel::options = static_cast<FFModel::Options const&>( options );
    }

  //! @brief FFModel's reset hook: the base's collocation state, then this solver's.
  void _reset
    ()
    override
    {
      _reset_model();
      _reset_colloc_state();
      _peel = t_Peel();
      _iniEqn.clear();
      _massByEqn.clear(); _mass.clear(); _colDiff.clear(); _rowAlg.clear(); _iniStream.clear();
      _evoLo = _evoHi = 0.;
      _shape = t_SemiDiscrete();
    }

  //! @brief OCBase's collocation-domain map: the working _mDom (the spatial directions, once set up).
  t_Dom const& _dom_colloc
    ()
    const
    override
    { return _mDom; }

  //! @brief Not marching, ever: t is the integrator's.
  MarchContext _march_context
    ()
    const
    override
    { return MarchContext(); }

  bool _on_setup_begin
    ()
    override
    {
      _setupStatus = SetupStatus::OK;
      _plan.reset();
      _blockSymbol.clear();
      _blockClassification.clear();
      _blockFaceData.clear();
      _idxCache.clear();
      _face_data.clear();
      _peel = t_Peel();
      _iniEqn.clear();
      _massByEqn.clear(); _mass.clear(); _colDiff.clear(); _rowAlg.clear(); _iniStream.clear();
      _shape = t_SemiDiscrete();
      return true;
    }

  //! @brief After reduce_order, before classification: the peel along t; the working model rewritten without t
  //! (INTERIOR and BOUNDARY rows replaced by their remainder R, INITIAL rows set aside, t removed from every
  //! domain set and mask, t registered as a scalar input of the collocated system); then the shared layout.
  bool _on_model_discretise
    ()
    override;

  //! @brief The shared plan over the spatial layout, then the per-node DOF audit and the [mol] report.
  bool _on_model_ready
    ()
    override;

  //! @brief The dependence set of an expression on the DAG variables it reads (leaf-seeded; needs no layout).
  FFDep _deps_of
    ( FFVar const& expr )
    const;

  //! @brief The column of field @p x at the node a row stands at, from the row's provenance record: pos_state for
  //! the element block, then the within-block tensor position in the field's direction order; a face-collocated
  //! direction resolves to the face node.  Returns false when the record does not place x (a direction x spans
  //! that the row's equation does not).
  bool _column_of
    ( FFVar const& x, t_RowNode const& rec, size_t& col )
    const;

  //! @brief Mc from the two halves of the row provenance: row r is equation k at node n, and k has mass entries
  //! (field i, coefficient) -> Mc[r, column of field i at n] = coefficient.  Fills _mass, _colDiff, _rowAlg.
  void _build_mass
    ();

  //! @brief Section 2.1 step 5, node-level: for every field, every node: a differential unknown needs an INITIAL row
  //! of its field covering the node (the Ini stream); an INITIAL row at an algebraic node is dropped, its value to be
  //! checked for consistency at t0.  Fills the stream and dropped counts; returns the number of differential unknowns
  //! WITHOUT an INITIAL row, each named -- the one genuine under-specification.
  size_t _audit_initial_coverage
    ( std::ostream& os );

  //! @brief IMPOSITION.NEUMANN_NODE == NEUMANN_DIFFERENTIAL: at a face where a boundary row reads a differential
  //! field's order-reduction auxiliary (a flux or Robin condition), the field's INTERIOR row is collocated at the face
  //! too and the auxiliary's LINK row is not -- so u at the node keeps its d_t row (and its initial condition) and the
  //! boundary row pins D.  Mask edits on the working model, before the layout.  Under NEUMANN_ALGEBRAIC (the default)
  //! nothing changes: the LINK row and the boundary row hold the node, u there is algebraic, as in OCFESLV.
  void _apply_neumann_rule
    ( FFVar const& t );
};

inline bool
OCFEMOL::_on_model_discretise
()
{
  if( !_evolution_dom_set || !_evolution_dom_var.dag() ){
    std::cerr << "OCFEMOL::setup ** the method of lines needs an evolution direction: set_evolution_domain( t )" << std::endl;
    _setupStatus = SetupStatus::EVOLUTION_DOMAIN_UNSUITABLE;
    return false;
  }
  FFVar const t = _evolution_dom_var;

  // 1. the peel on the WORKING model (order-reduced), roles resolved by now -- WITH the structural-index hook: the
  //    decomposition reads _mEqn, _mVar and the DAG only, by the equations' block ids (0 unless the model set them),
  //    so it runs here, on the full (t, z) model, before t is taken out.  An index above 1 along t refuses (phase 8).
  bool const peeled = _peel_evolution( t, _peel, false, true );

  // 1b. FFModel's classification of the FULL (t, z) model -- the layout and the classifier OCFESLV runs on the same
  //     model -- recorded as what the model IS.  (The plan's classification, later, is of the spatial operator at
  //     fixed t.)  The layout is rebuilt over the spatial directions once t is out.
  _shape.model_class.clear();
  if( _discretise_colloc_model() && _classify_pde() )
    for( auto const& [bid, cls] : _blockClassification ) _shape.model_class[bid] = cls.type;
  for( auto const& [bid, ty] : _shape.model_class )
    _disp( 1 ) << "OCFEMOL::setup [mol] model along " << t << ", block " << bid << ": " << pde_type_name( ty )
               << "  (the classification that follows is of the spatial operator at fixed " << t << ", which is what the plan collocates)\n";

  if( !peeled ){
    bool reduced = false;
    if( options.IMPOSITION.REDUCE_HIGH_INDEX && _peel.index > 1 && !_shape.model_class.empty() ){
      // EXPERIMENT (phase 8 brought forward): the full model is laid out and classified above; run the reduction
      // on it and peel again.
      _disp( 1 ) << "OCFEMOL::setup [mol] structural index " << _peel.index << " along " << t
                 << ": IMPOSITION.REDUCE_HIGH_INDEX -- reducing on the full model before collocation\n";
      _build_reduction_plan();
      if( _reduce_high_index() ){
        _disp( 1 ) << "OCFEMOL::setup [mol] reduction done: " << _mEqn.size() << " rows now; peeling again\n";
        reduced = _peel_evolution( t, _peel, false, true );
      }
    }
    if( !reduced ){
      std::cerr << "OCFEMOL::setup ** the model cannot be marched along " << t << ": " << _peel.blocker << std::endl;
      _setupStatus = SetupStatus::INCONSISTENT_MODEL;
      return false;
    }
  }

  // 2. the working model without t.  Row by row: an INTERIOR or BOUNDARY row becomes its remainder R (its d_t
  //    terms are Mc's, not Rc's); an INITIAL row leaves for the Ini stream; LINK, INTERFACE, SURFACE and
  //    DIAGNOSTIC rows stay as they are.  Every row loses its t mask.  _evoDerivRows -- the rows _defining_rows
  //    must not read as definitions -- are the rows with a mass entry, in the new numbering.
  if( options.IMPOSITION.NEUMANN_NODE == Options::NEUMANN_DIFFERENTIAL ) _apply_neumann_rule( t );

  std::vector<t_Eqn> kept;  kept.reserve( _mEqn.size() );
  std::vector<bool> is_evo_row;
  _massByEqn.clear();
  {
    std::map<size_t,t_PeelRow const*> by_eqn;
    for( auto const& r : _peel.rows ) by_eqn[r.eqn] = &r;
    std::set<size_t> ini( _peel.initial.begin(), _peel.initial.end() );
    for( size_t k = 0; k < _mEqn.size(); ++k ){
      t_Eqn e = _mEqn[k];
      e.dom.erase( t );
      if( ini.count( k ) ){ _iniEqn.push_back( std::move( e ) ); continue; }
      auto it = by_eqn.find( k );
      bool evo = false;
      std::vector< std::pair<size_t,FFVar> > mass;
      if( it != by_eqn.end() ){ e.var = it->second->rest; mass = it->second->mass; evo = !mass.empty(); }
      kept.push_back( std::move( e ) ); is_evo_row.push_back( evo ); _massByEqn.push_back( std::move( mass ) );
    }
  }
  _mEqn = std::move( kept );
  _evoDerivRows.clear();
  for( size_t k = 0; k < is_evo_row.size(); ++k ) if( is_evo_row[k] ) _evoDerivRows.insert( k );

  // t out of the domain sets and the collocation map; its span kept for the integrator
  {
    auto itd = _mDom.find( t );
    if( itd == _mDom.end() ){ _setupStatus = SetupStatus::INCONSISTENT_MODEL; return false; }
    _evoLo = itd->second.lo_dom;  _evoHi = itd->second.up_dom;
    _mDom.erase( itd );
    for( auto& [v, dm] : _mVar ) dm.erase( t );
    for( auto& [v, dm] : _mInp ) dm.erase( t );
    for( auto& [inp, disc] : _mInpDisc ) disc.erase( t );
    for( auto& f : _mFct ) f.grid.erase( t );
    for( auto& a : _auxDef ){ a.dom.erase( t ); a.diff_dom.erase( t ); }   // reduce_order recorded the auxiliaries over (t,z)
  }
  // t itself is a scalar input of the collocated system: the integrator supplies it at every evaluation
  _mInp[t] = std::set<FFVar,lt_FFVar>();

  _disp( 1 ) << "OCFEMOL::setup [mol] peel along " << t << ": " << _peel.rows.size() << " rows (INTERIOR+BOUNDARY), "
             << _peel.initial.size() << " INITIAL, " << _peel.link.size() << " LINK/INTERFACE; "
             << _peel.n_differential() << " of " << _peel.states.size() << " fields differential; structural index along "
             << t << " = " << _peel.index << "\n";

  // 3. the shared layout over the spatial directions
  return _discretise_colloc_model();
}

inline bool
OCFEMOL::_on_model_ready
()
{
  if( !_build_colloc_plan() ) return false;
  _issetup = true;

  auto const model_class = _shape.model_class;
  _shape = t_SemiDiscrete();
  _shape.model_class   = model_class;
  _shape.ran           = true;
  _shape.n_states      = _nCollVar;      // includes the explicit (kept) tau multipliers appended to var[] at plan commit
  _shape.n_tau         = _nTraceVar;     // ... which are this many of them (0 when every tau is eliminated internally)
  _shape.n_rows        = _nCollEqn;
  _shape.n_ini         = _iniEqn.size();
  _shape.n_fields      = _peel.states.size();
  _shape.n_diff_fields = _peel.n_differential();
  _shape.t_lo = _evoLo;  _shape.t_hi = _evoHi;

  _disp( 1 ) << "OCFEMOL::setup [mol] semi-discrete: " << _shape.n_states << " unknowns (" << _nCollVar - _nTraceVar << " nodal + "
             << _nTraceVar << " explicit tau), " << _shape.n_rows << " rows, " << _shape.n_ini << " INITIAL rows aside, t in ["
             << _shape.t_lo << ", " << _shape.t_hi << "]\n";

  _build_mass();
  if( !_build_jacobian_pattern() ){
    std::cerr << "OCFEMOL::setup ** the Jacobian pattern of Rc could not be built" << std::endl;
    _setupStatus = SetupStatus::INCONSISTENT_MODEL;
    return false;
  }
  _disp( 1 ) << "OCFEMOL::setup [mol] Jacobian of Rc: " << _jacNColor << " colours over " << _nCollVar + _nCollInp << " columns\n";
  _shape.n_diff_cols = 0; for( bool b : _colDiff ) _shape.n_diff_cols += b;
  _shape.n_alg_rows  = 0; for( bool b : _rowAlg )  _shape.n_alg_rows  += b;
  _disp( 1 ) << "OCFEMOL::setup [mol] Mc: " << _mass.size() << " entries; " << _shape.n_diff_cols << " differential unknowns, "
             << _shape.n_states - _shape.n_diff_cols << " algebraic; " << _shape.n_alg_rows << " algebraic rows\n";

  size_t const bad = _audit_initial_coverage( _disp( 1 ) );
  _disp( 1 ) << "OCFEMOL::setup [mol] Ini stream: " << _shape.n_ini_stream << " INITIAL (row, node) pairs at differential nodes; "
             << _shape.n_ini_dropped << " at algebraic nodes dropped (consistency at t0)\n";
  if( bad ){
    std::cerr << "OCFEMOL::setup ** " << bad << " differential unknown(s) without an initial condition (see the [mol] DOF audit)" << std::endl;
    _setupStatus = SetupStatus::INCONSISTENT_MODEL;
    return false;
  }
  return true;
}

inline bool
OCFEMOL::eval_residual
( double t, double const* xc, double* rc )
{
  if( !_shape.ran || !xc || !rc ) return false;
  // the inputs of the collocated system: t, and any lumped input of the model (zero until phase 5 wires them)
  std::vector<double> inp( _nCollInp, 0. );
  if( _nCollInp ) inp[ pos_input( _evolution_dom_var, {} ) ] = t;
  for( size_t i = 0; i < _vCst.size() && i < _vCstVal.size(); ++i ) _OCCst[i].set( _vCstVal[i] );
  _eval_eqn( rc, xc, _nCollInp? inp.data(): nullptr, _OCEqn, _OCVar, _OCInp, _OCDom, _OCCst, _wrkOC );
  return true;
}

inline void
OCFEMOL::jacobian_pattern
( std::vector<size_t>& jac_ptr, std::vector<size_t>& jac_col )
const
{
  jac_ptr.assign( 1, 0 ); jac_col.clear();
  for( size_t r = 0; r < _jacPatternEqn.size(); ++r ){
    for( size_t c : _jacPatternEqn[r] ) if( c < _nCollVar ) jac_col.push_back( c );   // the xc columns only
    jac_ptr.push_back( jac_col.size() );
  }
}

inline bool
OCFEMOL::deriv_residual
( double t, double const* xc, double* rc, std::vector<size_t>& jac_ptr, std::vector<size_t>& jac_col,
  std::vector<double>& jac_val, double* drdt )
{
  if( !_shape.ran || !_jacReady || !xc ) return false;
  size_t const ncol = _nCollVar + _nCollInp;
  // the FAD buffers, seeded by colour exactly as OCFESLV::_reset_autodiff seeds its own
  std::vector<t_FAD> Fvar( _nCollVar ), Finp( _nCollInp ), Feqn( _nCollEqn );
  for( size_t i = 0; i < _nCollVar; ++i ){ Fvar[i] = xc[i]; if( _jacNColor ) Fvar[i].diff( _jacColor[i], _jacNColor ); }
  size_t const it = _nCollInp? pos_input( _evolution_dom_var, {} ): 0;
  for( size_t i = 0; i < _nCollInp; ++i ){ Finp[i] = ( i == it? t: 0. ); if( _jacNColor ) Finp[i].diff( _jacColor[_nCollVar+i], _jacNColor ); }
  _reset_colloc( _OCFVar, _OCFInp, _OCFDom, _OCFCst );
  _OCFEqn = OCVar<t_FAD>();  _wrkOCF.clear();
  std::vector<t_FAD> Fcst( _vCstVal.begin(), _vCstVal.end() );          // the model's constants, as OCFESLV passes them
  for( size_t i = 0; i < _vCst.size() && i < Fcst.size(); ++i ) _OCFCst[i].update( &Fcst[i] );
  _eval_eqn( Feqn.data(), Fvar.empty()? nullptr: Fvar.data(), Finp.empty()? nullptr: Finp.data(),
             _OCFEqn, _OCFVar, _OCFInp, _OCFDom, _OCFCst, _wrkOCF );
  // extraction into CSR over the xc columns; t's column apart
  jac_ptr.assign( 1, 0 ); jac_col.clear(); jac_val.clear();
  for( size_t r = 0; r < _nCollEqn; ++r ){
    if( rc ) rc[r] = Feqn[r].val();
    if( drdt ) drdt[r] = 0.;
    for( size_t c : _jacPatternEqn[r] ){
      if( c >= ncol ) continue;
      double const d = _jacNColor? Feqn[r].deriv( _jacColor[c] ): 0.;
      if( c < _nCollVar ){ jac_col.push_back( c ); jac_val.push_back( d ); }
      else if( drdt && c == _nCollVar + it ) drdt[r] = d;
    }
    jac_ptr.push_back( jac_col.size() );
  }
  return true;
}

inline bool
OCFEMOL::eval_initial
( double t, double const* xc, double* ini )
{
  if( !_shape.ran || !xc || !ini ) return false;
  for( size_t k = 0; k < _iniStream.size(); ++k ){
    auto const& e = _iniStream[k];
    std::vector<FFVar> vars; std::vector<double> vals;
    for( auto const& [y, c] : e.cols ){ vars.push_back( y ); vals.push_back( xc[c] ); }
    for( auto const& [d, z] : e.coord ){ vars.push_back( d ); vals.push_back( z ); }
    if( _evolution_dom_set ){ vars.push_back( _evolution_dom_var ); vals.push_back( t ); }
    for( size_t i = 0; i < _vCst.size() && i < _vCstVal.size(); ++i ){ vars.push_back( _vCst[i] ); vals.push_back( _vCstVal[i] ); }
    std::vector<double> out( 1 );
    try{ _dag->eval( std::vector<FFVar>{ _iniEqn[e.row].var }, out, vars, vals ); } catch( ... ){ return false; }
    ini[k] = out[0];
  }
  return true;
}

inline bool
OCFEMOL::deriv_initial
( double t, double const* xc, double* ini, std::vector<size_t>& ptr, std::vector<size_t>& col, std::vector<double>& val )
{
  if( !_shape.ran || !xc ) return false;
  ptr.assign( 1, 0 ); col.clear(); val.clear();
  for( size_t k = 0; k < _iniStream.size(); ++k ){
    auto const& e = _iniStream[k];
    std::vector<FFVar> vars; std::vector<t_FAD> vals;
    for( size_t j = 0; j < e.cols.size(); ++j ){ vars.push_back( e.cols[j].first ); t_FAD f = xc[ e.cols[j].second ]; f.diff( j, e.cols.size() ); vals.push_back( f ); }
    for( auto const& [d, z] : e.coord ){ vars.push_back( d ); vals.push_back( t_FAD( z ) ); }
    if( _evolution_dom_set ){ vars.push_back( _evolution_dom_var ); vals.push_back( t_FAD( t ) ); }
    for( size_t i = 0; i < _vCst.size() && i < _vCstVal.size(); ++i ){ vars.push_back( _vCst[i] ); vals.push_back( t_FAD( _vCstVal[i] ) ); }
    std::vector<t_FAD> out( 1 );
    try{ _dag->eval( std::vector<FFVar>{ _iniEqn[e.row].var }, out, vars, vals ); } catch( ... ){ return false; }
    if( ini ) ini[k] = out[0].val();
    for( size_t j = 0; j < e.cols.size(); ++j ){ col.push_back( e.cols[j].second ); val.push_back( e.cols.size()? out[0].deriv( j ): 0. ); }
    ptr.push_back( col.size() );
  }
  return true;
}

inline FFDep
OCFEMOL::_deps_of
( FFVar const& expr )
const
{
  // the subgraph's leaves, each seeded as its own dependency: needs no layout, so it serves before the base has
  // built one (the Neumann rule runs before the layout) as well as after
  FFSubgraph sg = _dag->subgraph( 1, &expr );
  std::vector<FFVar> leaf;
  for( auto const* op : sg.l_op )
    if( op && op->type == FFOp::VAR && op->varout[0] ) leaf.push_back( *op->varout[0] );
  std::vector<FFDep> dleaf( leaf.size() ), ddep( 1 );
  for( size_t i = 0; i < leaf.size(); ++i ) dleaf[i].indep( (int)leaf[i].id().second );
  if( !leaf.empty() ){ try { _dag->eval( std::vector<FFVar>{ expr }, ddep, leaf, dleaf ); } catch( ... ){ return FFDep(); } }
  return ddep[0];
}

inline bool
OCFEMOL::_column_of
( FFVar const& x, t_RowNode const& rec, size_t& col )
const
{
  auto itv = _mVar.find( x );  auto itb = _blkVar.find( x );
  if( itv == _mVar.end() || itb == _blkVar.end() ) return false;
  auto const& var_dom = itv->second;
  auto const& blk_dim = std::get<2>( itb->second );
  std::map<FFVar,size_t,lt_FFVar> ndx_el;
  std::vector<size_t> node;  node.reserve( var_dom.size() );
  for( auto const& v : var_dom ){
    auto itr = std::find_if( rec.dims.begin(), rec.dims.end(), [&]( t_RowNode::Dim const& d ){ return d.dir_id == v.id().second; } );
    if( itr == rec.dims.end() ) return false;        // the row's equation does not span this direction of x
    ndx_el[v] = itr->iel;
    auto const& d = _mDom.at( v );
    node.push_back( itr->lim == FFDom::LB? 0: itr->lim == FFDom::UB? d.n_node - 1: itr->node );
  }
  col = pos_state( x, ndx_el );
  size_t stride = 1, k = 0;
  for( auto const& v : var_dom ){ (void)v; col += node[k] * stride; stride *= blk_dim[k]; ++k; }
  return true;
}

inline void
OCFEMOL::_build_mass
()
{
  _mass.clear();
  _colDiff.assign( _shape.n_states, false );
  _rowAlg.assign( _shape.n_rows, true );
  size_t const nrow = std::min( _rowEqnNdx.size(), _shape.n_rows );
  for( size_t r = 0; r < nrow; ++r ){
    int const k = _rowEqnNdx[r];
    if( k < 0 || (size_t)k >= _massByEqn.size() || _massByEqn[k].empty() ) continue;   // plan-owned or algebraic
    for( auto const& [i, coeff] : _massByEqn[k] ){
      size_t col = 0;
      if( !_column_of( _peel.states[i], _rowNode[r], col ) ) continue;
      _mass.push_back({ r, col, coeff });
      if( col < _colDiff.size() ) _colDiff[col] = true;
      _rowAlg[r] = false;
    }
  }
}

inline size_t
OCFEMOL::_audit_initial_coverage
( std::ostream& os )
{
  _shape.n_ini_stream = _shape.n_ini_dropped = 0;
  size_t bad = 0;
  // the INITIAL rows reading each field, with their spatial masks
  std::vector<FFDep> ini_dep;  ini_dep.reserve( _iniEqn.size() );
  for( auto const& e : _iniEqn ) ini_dep.push_back( _deps_of( e.var ) );

  for( size_t i = 0; i < _peel.states.size(); ++i ){
    FFVar const& x = _peel.states[i];
    auto itv = _mVar.find( x );  auto itb = _blkVar.find( x );
    if( itv == _mVar.end() || itb == _blkVar.end() ) continue;
    auto const& var_dom = itv->second;  auto const& blk_dim = std::get<2>( itb->second );
    if( var_dom.empty() ){                            // a lumped state (no spatial direction): one node
      size_t const col = pos_state( x, {} );
      bool const diff = col < _colDiff.size() && _colDiff[col];
      bool covered = false;
      for( size_t q = 0; q < _iniEqn.size() && !covered; ++q ) covered = ini_dep[q].dep().count( x.id().second ) > 0;
      if( covered ){
        if( diff ) ++_shape.n_ini_stream; else ++_shape.n_ini_dropped;
        for( size_t q = 0; q < _iniEqn.size(); ++q ){
          if( !ini_dep[q].dep().count( x.id().second ) ) continue;
          t_IniEntry e; e.row = q; e.target = col; e.dropped = !diff;
          for( auto const& y : _peel.states ) if( ini_dep[q].dep().count( y.id().second ) && _mVar.count( y ) && _mVar.at( y ).empty() )
            e.cols.emplace_back( y, pos_state( y, {} ) );
          _iniStream.push_back( std::move( e ) );
        }
      }
      else if( diff ){ os << "OCFEMOL::setup [mol] DOF audit: differential state " << x << " has no INITIAL row\n"; ++bad; }
      continue;
    }
    // iterate the field's nodes: element tuple x within-element node tuple
    std::vector<FFVar> dirs( var_dom.begin(), var_dom.end() );
    std::vector<size_t> nel( dirs.size() ), nnd( dirs.size() );
    for( size_t d = 0; d < dirs.size(); ++d ){ nel[d] = _mDom.at( dirs[d] ).n_elem; nnd[d] = blk_dim[d]; }
    std::vector<size_t> iel( dirs.size(), 0 ), ind( dirs.size(), 0 );
    size_t n_missing = 0;
    while( true ){
      // the column
      std::map<FFVar,size_t,lt_FFVar> ndx_el;  for( size_t d = 0; d < dirs.size(); ++d ) ndx_el[dirs[d]] = iel[d];
      size_t col = pos_state( x, ndx_el ), stride = 1;
      for( size_t d = 0; d < dirs.size(); ++d ){ col += ind[d] * stride; stride *= nnd[d]; }
      bool const diff = col < _colDiff.size() && _colDiff[col];
      // INITIAL coverage: a row reading x whose mask along every direction admits this node
      bool covered = false;
      for( size_t q = 0; q < _iniEqn.size() && !covered; ++q ){
        if( !ini_dep[q].dep().count( x.id().second ) ) continue;
        bool ok = true;
        for( size_t d = 0; d < dirs.size() && ok; ++d ){
          auto itm = _iniEqn[q].dom.find( dirs[d] );
          int const m = itm == _iniEqn[q].dom.end()? FFDom::ALL: itm->second;
          bool const atLB = iel[d] == 0 && ind[d] == 0, atUB = iel[d] + 1 == nel[d] && ind[d] + 1 == nnd[d];
          if( m == FFDom::LB ) ok = atLB;
          else if( m == FFDom::UB ) ok = atUB;
          else ok = !( ( atLB && ( m & 1 ) ) || ( atUB && ( m & 2 ) ) );
        }
        covered = ok;
      }
      if( covered ){
        if( diff ) ++_shape.n_ini_stream; else ++_shape.n_ini_dropped;
        // every INITIAL row reading x and admitting this node is an entry; its other fields' values at the same node
        for( size_t q = 0; q < _iniEqn.size(); ++q ){
          if( !ini_dep[q].dep().count( x.id().second ) ) continue;
          bool ok = true;
          for( size_t d = 0; d < dirs.size() && ok; ++d ){
            auto itm = _iniEqn[q].dom.find( dirs[d] );
            int const m = itm == _iniEqn[q].dom.end()? FFDom::ALL: itm->second;
            bool const atLB = iel[d] == 0 && ind[d] == 0, atUB = iel[d] + 1 == nel[d] && ind[d] + 1 == nnd[d];
            if( m == FFDom::LB ) ok = atLB; else if( m == FFDom::UB ) ok = atUB; else ok = !( ( atLB && ( m & 1 ) ) || ( atUB && ( m & 2 ) ) );
          }
          if( !ok ) continue;
          t_IniEntry e; e.row = q; e.target = col; e.dropped = !diff;
          for( size_t d = 0; d < dirs.size(); ++d ){
            FFDom const& dm = _mDom.at( dirs[d] );
            double const lo = dm.elem_lo( iel[d] ), hi = dm.elem_up( iel[d] ), xi = dm.lgnodes().at( ind[d] );
            e.coord.emplace_back( dirs[d], 0.5*( lo + hi ) + 0.5*( hi - lo )*xi );
          }
          for( size_t i2 = 0; i2 < _peel.states.size(); ++i2 ){
            FFVar const& y = _peel.states[i2];
            if( !ini_dep[q].dep().count( y.id().second ) ) continue;
            auto ity = _mVar.find( y ); auto itby = _blkVar.find( y );
            if( ity == _mVar.end() || itby == _blkVar.end() || ity->second != var_dom ) continue;   // same grid only
            e.cols.emplace_back( y, pos_state( y, ndx_el ) + ( col - pos_state( x, ndx_el ) ) );
          }
          _iniStream.push_back( std::move( e ) );
        }
      }
      else if( diff ) ++n_missing;
      // next node
      size_t d = 0;
      for( ; d < dirs.size(); ++d ){
        if( ++ind[d] < nnd[d] ) break;
        ind[d] = 0;
        if( ++iel[d] < nel[d] ) break;
        iel[d] = 0;
      }
      if( d == dirs.size() ) break;
    }
    if( n_missing ){
      os << "OCFEMOL::setup [mol] DOF audit: differential field " << x << " has " << n_missing << " node(s) with no INITIAL row\n";
      bad += n_missing;
    }
  }
  return bad;
}

inline void
OCFEMOL::_apply_neumann_rule
( FFVar const& t )
{
  // boundary rows (peel: BOUNDARY, no mass) that read an auxiliary of a differential field
  for( auto const& br : _peel.rows ){
    if( br.role != EqnRole::BOUNDARY || !br.mass.empty() ) continue;
    t_Eqn const& be = _mEqn[br.eqn];
    FFDep const dep = _deps_of( br.rest );
    for( auto const& a : _auxDef ){
      if( !dep.dep().count( a.aux.id().second ) ) continue;
      // the field the auxiliary differentiates, and whether it is differential
      FFDep const adep = _deps_of( a.expr );
      size_t ifield = _peel.states.size();
      for( size_t i = 0; i < _peel.states.size(); ++i )
        if( adep.dep().count( _peel.states[i].id().second ) && _peel.differential[i] ){ ifield = i; break; }
      if( ifield == _peel.states.size() ) continue;
      // the face: the boundary row's point mask along one spatial direction
      for( auto const& [dv, m] : be.dom ){
        if( dv.id() == t.id() || ( m != FFDom::LB && m != FFDom::UB ) ) continue;
        int const bit = m == FFDom::LB? 1: 2;
        // the field's INTERIOR rows: collocated at the face too
        for( auto const& ir : _peel.rows ){
          if( ir.role != EqnRole::INTERIOR ) continue;
          bool on_field = false;  for( auto const& [i, c] : ir.mass ) if( i == ifield ) on_field = true;
          if( !on_field ) continue;
          auto itm = _mEqn[ir.eqn].dom.find( dv );
          if( itm != _mEqn[ir.eqn].dom.end() && itm->second >= 0 ) itm->second &= ~bit;
        }
        // the auxiliary's LINK row: not at the face
        for( size_t k : _peel.link ){
          if( _mEqn[k].opt->role != EqnRole::LINK ) continue;
          if( !_deps_of( _mEqn[k].var ).dep().count( a.aux.id().second ) ) continue;
          auto itm = _mEqn[k].dom.find( dv );
          if( itm != _mEqn[k].dom.end() && itm->second >= 0 ) itm->second |= bit;
        }
        _disp( 1 ) << "OCFEMOL::setup [mol] NEUMANN_DIFFERENTIAL: field " << _peel.states[ifield] << " keeps its d_t row at the "
                   << ( m == FFDom::LB? "lower": "upper" ) << " face of " << dv << "; the auxiliary's LINK row leaves it\n";
      }
    }
  }
}


} // namespace mc

#endif
