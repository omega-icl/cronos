// Copyright (C) Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the EPL-2.0 with GPL-2.0-or-later as a Secondary License; see the LICENSE file.
// SPDX-License-Identifier: EPL-2.0 OR GPL-2.0-or-later

#ifndef CRONOS__OCFEMOL_HPP
#define CRONOS__OCFEMOL_HPP

#include <functional>
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
    // NEUMANN_ALGEBRAIC is OCFESLV's rule (leg 1 and leg 3 of DAESLV_molgate are then the same discrete system);
    // NEUMANN_DIFFERENTIAL is more accurate on the rod (2.3e-4 vs 7.0e-4 at 3x5) and differs from OCFESLV at that
    // node by spatial error.  Default NEUMANN_ALGEBRAIC (decision 31 Oct 2026: the exactness gate first).
    enum NeumannNode
    {
      NEUMANN_ALGEBRAIC = 0,   //!< the LINK row and the boundary row: u at the node is algebraic (no d_t row, no initial condition)
      NEUMANN_DIFFERENTIAL     //!< the PDE row and the boundary row: u at the node stays differential, with its initial condition
    };
    //! @brief How the method of lines treats a boundary node and a high-index model.
    struct t_Imposition
    {
      NeumannNode NEUMANN_NODE;   //!< What a flux boundary node is in the semi-discrete DAE (default NEUMANN_ALGEBRAIC)
      // REDUCE_HIGH_INDEX: FFModel's index reduction (dummy derivatives, hidden constraints) on the full (t, z) model
      // before t is taken out, then the reduced index-1 system collocated -- phase 8 of the MoL plan, an experiment.
      bool REDUCE_HIGH_INDEX;     //!< A model of structural index > 1 along t: reduce it first (default false: refuse it)
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
    size_t n_alg_rank  = 0;   //!< structural rank of the algebraic block (algebraic rows x algebraic unknowns): = n_alg_rows for index 1
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
  //! @brief A model input as the semi-discrete system sees it: its slots in the collocated input vector uc, in
  //! input_dofs() order (element blocks outer, nodes inner, the first direction fastest -- the order by-name values
  //! arrive in), their node coordinates at the input's own collocation, and whether the input varies in time (its
  //! declared domain held the evolution direction, recorded before setup takes t out).  t itself is time_slot().
  struct t_InputLayout
  {
    FFVar input;                                                     //!< the input, as a variable of the working model
    bool time_varying = false;
    std::vector<size_t> cols;                                        //!< its slots in uc
    std::vector< std::vector< std::pair<FFVar,double> > > nodes;     //!< each slot's coordinates
    std::vector<DofIndex> dof;                                       //!< each slot's (element, node) in the spatial directions
  };
  std::vector<t_InputLayout> const& inputs_layout() const { return _inpLayout; }
  //! @brief The collocation a time-varying input declared in the evolution direction (add_input's type and nodes per
  //! element), for the integrator's own input over t: false when it declared none (then: as the t domain).
  bool input_time_collocation
    ( FFVar const& u, FFDom::TYPE& type, size_t& n_node )
    const
    {
      if( !_evolution_dom_set ) return false;
      auto const itd = _usr._mInpDiscUsr.find( u );
      if( itd == _usr._mInpDiscUsr.cend() ) return false;
      auto const itc = itd->second.find( _evolution_dom_var );
      if( itc == itd->second.cend() ) return false;
      type = itc->second.type;  n_node = itc->second.n_node;  return true;
    }
  //! @brief Why the last setup() failed, when it did: the first refusal of the model, in words.
  std::string const& setup_error() const { return _setupError; }

  //! @brief The time part of each output (declaration order), split off at setup: the method of lines integrates in t,
  //! so an output is its SPATIAL part -- what the collocation evaluates at a given t, from (xc, uc) -- and a time part
  //! the integrator applies to it: a point in t (with its side), or an integral over t.  The spatial part is the
  //! output as declared with t taken out: a point in t removed from its point, an outermost integral or evaluation in
  //! t peeled off (re-integrated or re-evaluated over its other directions, which MC++ folds into the same node).
  struct t_OutputTime
  {
    enum Kind { POINT = 0, INTEGRAL, COMPOSED, SPATIAL, TRAJECTORY } kind = POINT;   // SPATIAL: a scalar part with no time of its own
    int    tmask = FFDom::ALL;        //!< TRAJECTORY: the output's mask in t -- one value per stage time it selects
    double tau  = 0.;                 //!< POINT: the time
    int    side = FFDom::MINUS;       //!< POINT: the side, at a time where the value may differ on either side
    //! COMPOSED (7 Nov 2026): the output is a function h of deferred values (reductions in t); each is a spatial
    //! output of its own, appended to the outputs, and the declaration composes h over their reductions
    FFVar               h;            //!< the function, in this model's dag, over the deferred inputs dvars
    std::vector<FFVar>  dvars;        //!< the deferred inputs h reads
    std::vector<size_t> parts;        //!< their spatial outputs' indices
    //! A spatial reduction (an integral over the whole of a direction, a value at a point) nested inside a
    //! NONLINEARITY in the spatial part: the collocation's assembly sums an output's per-element values at the root,
    //! exact only for what is linear in the reductions -- so each such reduction is a scalar spatial part of its own
    //! (sparts), and the spatial part is hs composed over them at every t, on the semi-discrete dag.
    FFVar               hs;           //!< the spatial part as a function of sph (empty: the rows themselves)
    std::vector<FFVar>  sph;          //!< placeholders for the nested reductions
    std::vector<size_t> sparts;       //!< the nested reductions' spatial outputs' indices
  };
  std::vector<t_OutputTime> const& outputs_time() const { return _fctTime; }
  //! @brief The number of output values (n_colloc_fct(): a point output one, a distributed one per node), in the
  //! declared outputs' order, each output's values from its row0.
  size_t n_outputs() const { return _nCollFct; }
  //! @brief The outputs' spatial parts at (t, xc, uc): @p f has n_outputs() entries.  uc as for eval_residual.
  bool eval_outputs
    ( double t, double const* xc, double const* uc, double* f );
  //! @brief Their Jacobian over [xc | uc] in CSR (columns as deriv_residual_full's): the same FAD seeding as the
  //! residual's -- the column colouring covers residual and output rows together -- through the output evaluation.
  bool deriv_outputs_full
    ( double t, double const* xc, double const* uc, double* f,
      std::vector<size_t>& ptr, std::vector<size_t>& col, std::vector<double>& val );
  //! @brief The outputs' xc columns (as jacobian_pattern's) and their input slots other than time_slot() (as
  //! input_pattern's), in CSR.
  void output_state_pattern
    ( std::vector<size_t>& ptr, std::vector<size_t>& col )
    const
    {
      auto const& pat = _nsActive? _jacPatternFctRed: _jacPatternFct;  size_t const nred = n_reduced();
      ptr.assign( 1, 0 ); col.clear();
      for( size_t r = 0; r < pat.size(); ++r ){
        for( size_t c : pat[r] ) if( c < nred ) col.push_back( c );
        ptr.push_back( col.size() );
      }
    }
  void output_input_pattern
    ( std::vector<size_t>& ptr, std::vector<size_t>& col )
    const
    {
      size_t const it = _nCollInp? pos_input( _evolution_dom_var, {} ): 0, nred = n_reduced();
      auto const& pat = _nsActive? _jacPatternFctRed: _jacPatternFct;
      ptr.assign( 1, 0 ); col.clear();
      for( size_t r = 0; r < pat.size(); ++r ){
        for( size_t c : pat[r] ) if( c >= nred && c < nred + _nCollInp && c - nred != it ) col.push_back( c - nred );
        ptr.push_back( col.size() );
      }
    }
  //! @brief Output @p k's values in the n_outputs() array: (its first, its count) -- a point output one, a distributed
  //! one per node.
  std::pair<size_t,size_t> output_rows( size_t k ) const { return { _mFct.at( k ).row0, _mFct.at( k ).nrow }; }
  //! @brief The outputs the model DECLARED: the parts of composed outputs (t_OutputTime::COMPOSED) follow them.
  size_t n_declared_outputs() const { return _nFctDeclared; }
  //! @brief This model's constants and their values, for a declaration composing an expression of this model's dag.
  std::vector<FFVar> const&  constant_vars()   const { return _vCst; }
  std::vector<double> const& constant_values() const { return _vCstVal; }
  FFGraph* dag() const { return _dag; }
  //! @brief The pattern of deriv_outputs_full.
  void output_pattern_full
    ( std::vector<size_t>& ptr, std::vector<size_t>& col )
    const
    {
      size_t const ncol = n_reduced() + _nCollInp;
      auto const& pat = _nsActive? _jacPatternFctRed: _jacPatternFct;
      ptr.assign( 1, 0 ); col.clear();
      for( size_t r = 0; r < pat.size(); ++r ){
        for( size_t c : pat[r] ) if( c < ncol ) col.push_back( c );
        ptr.push_back( col.size() );
      }
    }
  //! @brief Size of the collocated input vector uc, the evolution variable's slot included.
  size_t n_input_slots() const { return _nCollInp; }
  //! @brief The evolution variable's slot in uc.
  size_t time_slot() const { return _nCollInp? pos_input( _evolution_dom_var, {} ): 0; }

  //! @brief The columns of @p field, in MESH ORDER: element-major, node within element, the directions in the
  //! field's own order (the order a driver's own mesh loop produces).  Empty for an unknown field.
  //! A field's columns are NOT contiguous -- the layout interleaves fields element by element -- so this is the
  //! only correct way to collect them.
  std::vector<size_t> const& field_columns
    ( FFVar const& field )
    const
    { static std::vector<size_t> const none; auto it = _fieldCols.find( field ); return it == _fieldCols.end()? none: it->second; }

  //! @brief The coordinates of @p field's nodes, in the same order as field_columns: one entry per column, each the
  //! (direction, coordinate) pairs of that node.  A driver can plot or compare without rebuilding the mesh.
  std::vector< std::vector< std::pair<FFVar,double> > > const& field_nodes
    ( FFVar const& field )
    const
    { static std::vector< std::vector< std::pair<FFVar,double> > > const none; auto it = _fieldNodes.find( field );
      return it == _fieldNodes.end()? none: it->second; }

  //! @brief The field that owns unknown @p col and its reference value; false for an explicit tau multiplier.
  bool column_field
    ( size_t col, FFVar& field, double& ref )
    const
    {
      if( col >= _colOwner.size() || _colOwner[col] >= _peel.states.size() ) return false;
      field = _peel.states[ _colOwner[col] ];  ref = this->ref( field );
      return true;
    }

  //! @brief Whether row @p row of Rc is algebraic (no mass entry).
  bool is_algebraic_row( size_t row ) const { return row < _rowAlg.size() && _rowAlg[row]; }

  //! @brief Rc(t, xc): the collocated residual, one evaluation of the interior expressions as vectors over the
  //! nodes through the OCVar environment, with t as the scalar input and the model's constants set.  @p xc has
  //! semi_discrete().n_states entries, @p rc semi_discrete().n_rows.  The mass term Mc(xc) d_t xc is not included:
  //! the integrator's residual is Mc xdot - Rc, with Mc from mass().  Phase 4 of the MoL plan.
  bool eval_residual
    ( double t, double const* xc, double* rc )
    { return eval_residual( t, xc, nullptr, rc ); }
  //! @brief Rc at (t, xc, uc): @p uc is the collocated input vector (n_input_slots(); its time_slot() is overwritten
  //! with @p t), nullptr for a model without inputs.
  bool eval_residual
    ( double t, double const* xc, double const* uc, double* rc );
  //! @brief The Jacobian of Rc over [xc | uc] in CSR: column c < n_states is xc[c], column n_states + j is uc[j] -- the
  //! time slot included, so dRc/dt is column n_states + time_slot().  One FAD sweep, uc coloured with xc: what a
  //! parameter sensitivity needs (dRc/dp) comes from the same evaluation as the integrator's Jacobian.
  bool deriv_residual_full
    ( double t, double const* xc, double const* uc, double* rc,
      std::vector<size_t>& ptr, std::vector<size_t>& col, std::vector<double>& val );
  //! @brief The pattern of deriv_residual_full.
  void jacobian_pattern_full
    ( std::vector<size_t>& ptr, std::vector<size_t>& col )
    const;
  //! @brief Per row, the slots j of uc other than time_slot() it depends on, in CSR (col holds j, not n_states + j).
  void input_pattern
    ( std::vector<size_t>& ptr, std::vector<size_t>& col )
    const
    {
      size_t const it = _nCollInp? pos_input( _evolution_dom_var, {} ): 0, nred = n_reduced();
      auto const& pat = _nsActive? _jacPatternRed: _jacPatternEqn;
      ptr.assign( 1, 0 ); col.clear();
      for( size_t r = 0; r < pat.size(); ++r ){
        for( size_t c : pat[r] ) if( c >= nred && c < nred + _nCollInp && c - nred != it ) col.push_back( c - nred );
        ptr.push_back( col.size() );
      }
    }

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
    //! the inputs the row reads, each as its value AT THE NODE: a combination of its slots in uc, (slot, weight) -- one
    //! slot with weight 1 where the input's collocation is the field's, the input's Lagrange interpolant otherwise
    std::vector< std::pair< FFVar, std::vector< std::pair<size_t,double> > > > icols;
    size_t row = 0;                                     //!< index into ini_rows()
    size_t target = 0;                                  //!< the differential unknown this condition fixes
    std::vector< std::pair<FFVar,size_t> > cols;        //!< (field, column of its value at the node)
    std::vector< std::pair<FFVar,double> > coord;       //!< (direction, coordinate of the node)
    bool dropped = false;                               //!< at an algebraic node: not in the stream, checked for consistency at t0
  };
  //! @brief The Ini stream (entries with dropped == false) and the dropped pairs (dropped == true), after setup.
  std::vector<t_IniEntry> const& ini_stream() const { return _iniStream; }

  //! @brief Ini entry @p k as an ORDINARY scalar expression in the DAG of @p xc: the model's INITIAL row with each
  //! field it reads bound to that field's unknown at the entry's node, the coordinates as constants, the evolution
  //! variable as @p t and the model's constants by value -- what a hand-written semi-discrete model writes.  An
  //! integrator peels INITIAL rows one at a time by composing each with its target state; as outputs of one vectorised
  //! operation every such compose rebuilt the whole operation (measured: 415 MB of heap at 144 unknowns, all of it
  //! FFMOLINI re-insertions and their gradients).  Pointwise rows are evaluated once at t0 and gain nothing from
  //! vectorisation.
  //! The inputs the row reads are bound to @p uc, one variable per slot of the collocated input vector (its time slot
  //! unused); a row reading an input with no @p uc given is refused (std::invalid_argument), never bound to zero.
  FFVar initial_expression
    ( size_t k, std::vector<FFVar> const& xc, FFVar const& t, std::vector<FFVar> const& uc = {} )
    const
    {
      auto const& e = _iniStream.at( k );
      std::vector<FFVar> vars, vals;
      for( auto const& [y, c] : e.cols ){ vars.push_back( y ); vals.push_back( xc.at( c ) ); }
      for( auto const& [w, comb] : e.icols ){
        FFVar v( 0. );
        for( auto const& [c, a] : comb ){
          if( c >= uc.size() ) throw std::invalid_argument( "OCFEMOL::initial_expression ** an INITIAL row reads the input "
                                                           + w.name() + ": its variables must be given (uc)" );
          v = ( comb.size() == 1 && a == 1. )? uc[c]: v + a*uc[c];      // the field's collocation: the slot itself
        }
        vars.push_back( w );  vals.push_back( v );
      }
      for( auto const& [d, z] : e.coord ){ vars.push_back( d ); vals.push_back( FFVar( z ) ); }
      if( _evolution_dom_set ){ vars.push_back( _evolution_dom_var ); vals.push_back( t ); }
      for( size_t i = 0; i < _vCst.size() && i < _vCstVal.size(); ++i ){ vars.push_back( _vCst[i] ); vals.push_back( FFVar( _vCstVal[i] ) ); }
      std::vector<FFVar> out( 1 );
      _dag->eval( std::vector<FFVar>{ _iniEqn[e.row].var }, out, vars, vals );
      return out[0];
    }

  //! @brief Transitions are taken (MoL 6c): registered into the working model by FFModel::setup like the rows, and
  //! handed to the integrator in their semi-discrete form (transition_expressions), not solved here.
  bool _transitions_supported() const override { return true; }

  //! @brief The semi-discrete form of the declared transition @p k (MoL 6c): per component, at each DIFFERENTIAL node
  //! of the state it determines (the nodes that carry INITIAL rows), its left side and its right side composed at the
  //! node -- the fields they read as their columns in @p xc, the coordinates by value, the evolution variable as @p t,
  //! the constants by value.  The algebraic nodes are re-solved by their rows at tau^+ (the integrator's stop system).
  //! The determined state is the record's, or inferred as the ONE state the right side reads.  A component reading an
  //! input, or determining a lumped state with no INITIAL row here, is refused (false, @p why).
  bool transition_expressions
    ( size_t k, std::vector<FFVar> const& xc, FFVar const& t, std::vector<FFVar>& left, std::vector<FFVar>& right,
      std::string& why, std::vector<FFVar> const& uc = {} )
    const
    {
      left.clear();  right.clear();
      if( k >= var_transition().size() ){ why = "no such transition"; return false; }
      auto const& tr = var_transition()[k];                       // the WORKING copy: its expressions live in _dag
      // the states a side reads, and the INPUTS (7 Nov 2026): an input uniform in space -- a parameter, or an input over
      // t only -- is its one slot of uc, bound by the integrator at tau as any input of a transition; an input
      // distributed in space is refused (its value at the node would be an interpolant: not taken yet)
      std::set<FFVar,lt_FFVar> inputs;
      auto states_of = [&]( FFVar const& e, std::set<FFVar,lt_FFVar>& st, std::string& bad ){
        auto sg = _dag->subgraph( 1, &e );
        for( auto const* v : sg.v_indep ){
          if( _mVar.count( *v ) ){ st.insert( *v );  continue; }
          if( _evolution_dom_set && v->id() == _evolution_dom_var.id() ) continue;   // t: the integrator's
          auto iti = _mInp.find( *v );
          if( iti == _mInp.end() ) continue;
          if( !iti->second.empty() ){ bad = v->name();  continue; }
          inputs.insert( *v );
        } };
      for( size_t i = 0; i < tr.right.size() && i < tr.left.size(); ++i ){
        std::set<FFVar,lt_FFVar> sr, sl;  std::string bad;
        states_of( tr.right[i], sr, bad );  states_of( tr.left[i], sl, bad );
        if( !bad.empty() ){ why = "a transition component reads the input " + bad + ", distributed in space: not taken by the method of lines yet"; return false; }
        if( !inputs.empty() && uc.empty() ){ why = "a transition component reads an input: its semi-discrete variables must be given (uc)"; return false; }
        FFVar s;
        if( i < tr.states.size() && tr.states[i].dag() ) s = tr.states[i];
        else if( sr.size() == 1 ) s = *sr.begin();
        else { why = "a transition component determines no single state (name it: add_transition's states)"; return false; }
        auto itc = _fieldCols.find( s );  auto itn = _fieldNodes.find( s );
        if( itc == _fieldCols.end() || itn == _fieldNodes.end() ){ why = "the state " + s.name() + " has no nodes here"; return false; }
        auto column_at = [&]( FFVar const& f, std::vector< std::pair<FFVar,double> > const& coord, size_t& col ){
          auto ic = _fieldCols.find( f );  auto in = _fieldNodes.find( f );
          if( ic == _fieldCols.end() || in == _fieldNodes.end() ) return false;
          for( size_t j = 0; j < in->second.size(); ++j ){
            if( in->second[j].size() != coord.size() ) continue;
            bool same = true;
            for( size_t d = 0; d < coord.size() && same; ++d ){
              same = false;
              for( auto const& [dv, z] : in->second[j] ) if( dv.id() == coord[d].first.id() && z == coord[d].second ){ same = true; break; }
            }
            if( same ){ col = ic->second[j]; return true; }
          }
          return false; };
        size_t nnode = 0;
        for( auto const& e : _iniStream ){
          if( e.dropped ) continue;
          size_t j = itc->second.size();
          for( size_t jj = 0; jj < itc->second.size(); ++jj ) if( itc->second[jj] == e.target ){ j = jj; break; }
          if( j == itc->second.size() ) continue;                       // the entry is another state's
          auto const& coord = itn->second[j];
          std::vector<FFVar> vars, vals;
          std::set<FFVar,lt_FFVar> all( sr );  all.insert( sl.begin(), sl.end() );
          for( auto const& f : all ){
            size_t c = 0;
            if( !column_at( f, coord, c ) ){ why = "the state " + f.name() + " has no node where " + s.name() + " has one"; return false; }
            vars.push_back( f );  vals.push_back( xc.at( c ) );
          }
          for( auto const& [dv, z] : coord ){ vars.push_back( dv );  vals.push_back( FFVar( z ) ); }
          if( _evolution_dom_set ){ vars.push_back( _evolution_dom_var );  vals.push_back( t ); }
          for( auto const& w : inputs ){                             // an input uniform in space: its one slot
            size_t const c = pos_input( w, {} );
            if( c >= uc.size() ){ why = "the input " + w.name() + " of a transition has no semi-discrete variable"; return false; }
            vars.push_back( w );  vals.push_back( uc[c] );
          }
          for( size_t q = 0; q < _vCst.size() && q < _vCstVal.size(); ++q ){ vars.push_back( _vCst[q] );  vals.push_back( FFVar( _vCstVal[q] ) ); }
          std::vector<FFVar> out( 2 );
          _dag->eval( std::vector<FFVar>{ tr.left[i], tr.right[i] }, out, vars, vals );
          left.push_back( out[0] );  right.push_back( out[1] );  ++nnode;
        }
        if( !nnode ){ why = "the state " + s.name() + " has no differential node to reset"; return false; }
      }
      return !left.empty();
    }
  size_t n_transitions() const { return var_transition().size(); }
  double transition_time( size_t k ) const { return var_transition().at( k ).tau; }

  //! @brief The INITIAL rows at the stream's nodes: @p ini has ini_stream().size() entries in that order (the dropped
  //! ones included, so the consistency check at t0 reads the same array).  Each is one small DAG evaluation.
  //! @p uc: the collocated input vector, for rows that read inputs (nullptr: such a row makes the call fail).
  bool eval_initial
    ( double t, double const* xc, double* ini, double const* uc = nullptr );

  //! @brief The derivative of eval_initial with respect to xc, in CSR over the entries (a few columns each).
  bool deriv_initial
    ( double t, double const* xc, double* ini, std::vector<size_t>& ptr, std::vector<size_t>& col, std::vector<double>& val,
      double const* uc = nullptr );

  //! @brief A SEED for an integrator's own consistent initialisation, without solving anything: the base's reference
  //! initialisation (every unknown from the model's references, the order-reduction auxiliaries from their defining
  //! relations), then the differential unknowns set from the Ini stream.  The algebraic block is left to the
  //! integrator -- DAESLV's KINSOL corrector solves the INITIAL rows with the interior rows at t0, sparse, which is the
  //! consistent initialisation; OCFEMOL only makes its start good.
  bool init_seed
    ( double t, double* xc );

  //! @brief The dropped Ini entries evaluated at @p xc: an INITIAL condition at an algebraic node is not imposed, so
  //! it must already hold there.  Returns the largest |residual| over them (0 when none were dropped).
  double initial_consistency_defect
    ( double t, double const* xc );

  //! @brief The Jacobian pattern over xc alone (CSR row pointers and columns), valid after setup; the values come
  //! from deriv_residual, in the same order.
  void jacobian_pattern
    ( std::vector<size_t>& jac_ptr, std::vector<size_t>& jac_col )
    const;

protected:

  t_Peel                _peel;      //!< the peel along t (working model)
  std::vector< std::vector< std::pair<size_t,FFVar> > > _massByEqn;   //!< per row of the rewritten _mEqn: (field, coefficient)
  std::vector<t_MassEntry> _mass;    //!< Mc as data
  std::set<FFVar,lt_FFVar>     _inpTimeVarying;   //!< inputs whose declared domain held t (setup takes it out)
  std::vector<std::string>     _auditRefusals;    //!< what the INITIAL audit refused, in words
  std::vector<t_OutputTime>    _fctTime;          //!< see outputs_time()

  //! @brief The structural index of the semi-discrete system (see setup, 5b): a maximum matching of the algebraic rows
  //! to the algebraic unknowns on Rc's dependency pattern.  Full: index 1, n_alg_rank = n_alg_rows.  Deficient: index
  //! >= 2 -- the unmatched unknowns and rows are reported by name (a multiplier, or a field and its node; a constraint
  //! row reading differential unknowns only) and false is returned.
  bool _structural_index_check( std::ostream& os, std::vector<std::vector<size_t>> const& pat, size_t n, std::vector<bool> const& colDiff,
                                std::vector<bool> const& rowAlg, size_t ntau, std::vector<size_t>* badCols = nullptr, std::vector<size_t>* badRows = nullptr )
  {
    size_t const m = pat.size();
    std::vector<long> acol_of( n, -1 );  std::vector<size_t> acols;
    for( size_t j = 0; j < n; ++j ) if( j >= colDiff.size() || !colDiff[j] ){ acol_of[j] = (long)acols.size();  acols.push_back( j ); }
    std::vector<size_t> arows;  for( size_t r = 0; r < m; ++r ) if( r < rowAlg.size() && rowAlg[r] ) arows.push_back( r );
    std::vector< std::vector<size_t> > adj( arows.size() );
    for( size_t i = 0; i < arows.size(); ++i ) for( size_t c : pat[ arows[i] ] ) if( c < n && acol_of[c] >= 0 ) adj[i].push_back( (size_t)acol_of[c] );
    std::vector<long> mcol( acols.size(), -1 ), mrow( arows.size(), -1 );
    std::function<bool(size_t,std::vector<bool>&)> aug = [&]( size_t i, std::vector<bool>& seen ){
      for( size_t c : adj[i] ){ if( seen[c] ) continue;  seen[c] = true;
        if( mcol[c] < 0 || aug( (size_t)mcol[c], seen ) ){ mcol[c] = (long)i;  mrow[i] = (long)c;  return true; } }
      return false; };
    size_t rank = 0;
    for( size_t i = 0; i < arows.size(); ++i ){ std::vector<bool> seen( acols.size(), false );  if( aug( i, seen ) ) ++rank; }
    _shape.n_alg_rank = rank;
    if( rank == acols.size() && rank == arows.size() ){
      os << "OCFEMOL::setup [mol] structural index of the semi-discrete system: 1 (" << rank << " algebraic rows matched to " << rank << " algebraic unknowns)\n";
      return true;
    }
    std::ostringstream why;
    why << "the semi-discrete system has structural index > 1: " << rank << " of " << arows.size() << " algebraic rows match "
        << acols.size() << " algebraic unknowns.";
    size_t named = 0;
    for( size_t c = 0; c < acols.size(); ++c ){
      if( mcol[c] >= 0 ) continue;
      size_t const j = acols[c];
      if( badCols ) badCols->push_back( j );
      if( named >= 6 ) continue;
      if( j >= n - ntau ) why << "  Unknown " << j << " is an interface multiplier minted by the plan with no algebraic row (its rows are differential: the index-2 form).";
      else { FFVar f;  double ref = 0.;  column_field( j, f, ref );  why << "  Unknown " << j << " (" << f.name() << " at a node) has no algebraic row of its own."; }
      ++named;
    }
    for( size_t i = 0; i < arows.size(); ++i ){
      if( mrow[i] >= 0 ) continue;
      if( badRows ) badRows->push_back( arows[i] );
      if( named >= 10 ) continue;
      bool diffOnly = true;  for( size_t c : pat[ arows[i] ] ) if( c < n && acol_of[c] >= 0 ) diffOnly = false;
      why << "  Row " << arows[i] << " is an algebraic row" << ( diffOnly? " reading differential unknowns only (a constraint on them).": " that no free algebraic unknown is left for." );
      ++named;
    }
    _setupError = why.str();
    return false;
  }
  std::set<size_t>             _fctSpatial;       //!< the SPATIAL parts' indices, appended by _split_output_time
  size_t                       _nFctDeclared = 0;  //!< see n_declared_outputs()
  bool _split_output_time
    ( FFVar const& t );
  std::string                  _setupError;       //!< see setup_error()
  std::vector<t_InputLayout>   _inpLayout;        //!< see inputs_layout()
  void _build_input_layout
    ();
  std::vector<t_IniEntry>  _iniStream;   //!< the Ini stream and the dropped pairs
  std::map<FFVar,std::vector<size_t>,lt_FFVar> _fieldCols;  //!< per field: its columns in mesh order
  std::map<FFVar,std::vector< std::vector< std::pair<FFVar,double> > >,lt_FFVar> _fieldNodes;  //!< ... and their coordinates
  std::vector<size_t>      _colOwner;    //!< per column: the field's index into _peel.states (npos for a tau)

  //! @brief Fill _fieldCols and _colOwner by walking each field's (element, node) grid.
  void _build_column_map
    ();
  std::vector<bool>     _colDiff;    //!< per unknown of xc
  std::vector<bool>     _rowAlg;     //!< per row of Rc

  // ---- NULL-SPACE ELIMINATION of plan multipliers that raise the index (phase 9, step 2; 8 Nov 2026).  A multiplier
  //      lambda that the structural index check finds with no algebraic row enforces a constraint x_L - x_R = 0 between
  //      two DIFFERENTIAL copies of a node and enters only their two rows, with opposite signs: the Hessenberg index-2
  //      form.  Its null space is the merge: x_L = x_R = y, the two rows summed (lambda cancels), the constraint row and
  //      the multiplier column dropped -- exact, drift-free, index 1.  Applied as a layer: the base keeps its full
  //      layout; the integrator sees the reduced unknowns (expand on entry, contract rows and columns on exit).
  bool                  _nsActive = false;
  size_t                _nsN = 0, _nsM = 0, _nsTau = 0;   //!< reduced counts
  std::vector<size_t>   _nsCol;      //!< full column -> reduced column (npos: an eliminated multiplier)
  std::vector<size_t>   _nsRow;      //!< full row -> reduced row (npos: a dropped constraint row)
  std::vector<double>   _nsX, _nsR;  //!< work: the full state and residual
  std::vector< std::vector<size_t> > _jacPatternRed, _jacPatternFctRed;   //!< the patterns over the reduced columns (+ inputs)
  size_t n_reduced() const { return _nsActive? _nsN: _nCollVar; }
  double const* _expand( double const* y )
    { if( !_nsActive ) return y;  _nsX.assign( _nCollVar, 0. );
      for( size_t c = 0; c < _nCollVar; ++c ) if( _nsCol[c] != std::string::npos ) _nsX[c] = y[ _nsCol[c] ];
      return _nsX.data(); }
  void _contract_rows( double const* rfull, double* r ) const
    { std::fill( r, r + _nsM, 0. );  for( size_t i = 0; i < _nCollEqn; ++i ) if( _nsRow[i] != std::string::npos ) r[ _nsRow[i] ] += rfull[i]; }
  size_t _col_red( size_t c ) const { return c < _nCollVar? _nsCol[c]: _nsN + ( c - _nCollVar ); }
  void _contract_csr( std::vector<size_t>& ptr, std::vector<size_t>& col, std::vector<double>& val, size_t mred ) const
    { std::vector< std::map<size_t,double> > acc( mred );
      for( size_t r = 0; r + 1 < ptr.size(); ++r ){ size_t const rr = _nsRow[r];  if( rr == std::string::npos ) continue;
        for( size_t k = ptr[r]; k < ptr[r+1]; ++k ){ size_t const cc = _col_red( col[k] );  if( cc == std::string::npos ) continue;  acc[rr][cc] += val[k]; } }
      ptr.assign( 1, 0 ); col.clear(); val.clear();
      for( auto const& row : acc ){ for( auto const& [c, v] : row ){ col.push_back( c );  val.push_back( v ); }  ptr.push_back( col.size() ); } }
  void _contract_csr_rows_only( std::vector<size_t>& ptr, std::vector<size_t>& col, std::vector<double>& val, std::vector<size_t> const& rowmap, size_t mred ) const
    { std::vector< std::map<size_t,double> > acc( mred );
      for( size_t r = 0; r + 1 < ptr.size(); ++r ){ size_t const rr = rowmap[r];  if( rr == std::string::npos ) continue;
        for( size_t k = ptr[r]; k < ptr[r+1]; ++k ){ size_t const cc = _col_red( col[k] );  if( cc == std::string::npos ) continue;  acc[rr][cc] += val[k]; } }
      ptr.assign( 1, 0 ); col.clear(); val.clear();
      for( auto const& row : acc ){ for( auto const& [c, v] : row ){ col.push_back( c );  val.push_back( v ); }  ptr.push_back( col.size() ); } }
  bool _eliminate_multipliers( std::ostream& os, std::vector<size_t> const& badCols, std::vector<size_t> const& badRows, std::string& why );
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
      _fieldCols.clear(); _fieldNodes.clear(); _colOwner.clear();
      _evoLo = _evoHi = 0.;
      _shape = t_SemiDiscrete();
    }

  //! @brief OCBase's collocation-domain map: the working _mDom (the spatial directions, once set up).
  t_Dom const& _dom_colloc
    ()
    const
    override
    { return _mDom; }


  //! @brief The residual and the outputs are swept separately here, each on its own colouring.
  bool _colour_rows_separately() const override { return true; }

  //! @brief No element of the evolution direction is collocated: t is a scalar input, the integrator's.
  EvolutionWindow _evolution_window
    ()
    const
    override
    { EvolutionWindow w;  w.kind = EvolutionWindow::NONE;  return w; }

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
  // IC_WEAK under the method of lines takes the ENERGY-STABLE SAT form: the legacy form gives the semi-discrete operator
  // a growing mode quadratic in SAT_SIGMA0 (measured 8 Nov 2026), which OCFESLV's Radau damps and an integrator cannot
  // step over.  Switched here, before the interface plan computes its strengths (OCFEBASE::_sat_strength).
  if( options.INTERFACE.IMPOSITION == Options::IC_WEAK && options.INTERFACE.SAT_FORM == Options::SAT_LEGACY ){
    options.INTERFACE.SAT_FORM = Options::SAT_ENERGY;
    _disp( 1 ) << "OCFEMOL::setup [mol] IC_WEAK: the energy-stable SAT form (INTERFACE.SAT_FORM = SAT_ENERGY) -- the legacy"
                  " form is not stable under the method of lines\n";
  }
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
    _inpTimeVarying.clear();
    for( auto const& [v, dm] : _mInp ) if( dm.count( t ) ) _inpTimeVarying.insert( v );
    for( auto& [v, dm] : _mVar ) dm.erase( t );
    for( auto& [v, dm] : _mInp ) dm.erase( t );
    for( auto& [inp, disc] : _mInpDisc ) disc.erase( t );
    if( !_split_output_time( t ) ){ _setupStatus = SetupStatus::INCONSISTENT_MODEL; return false; }
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
  _setupError.clear();
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

  _build_column_map();
  _build_input_layout();
  _build_mass();
  if( !_build_jacobian_pattern() ){
    std::cerr << "OCFEMOL::setup ** the Jacobian pattern of Rc could not be built" << std::endl;
    _setupStatus = SetupStatus::INCONSISTENT_MODEL;
    return false;
  }
  _disp( 1 ) << "OCFEMOL::setup [mol] Jacobian of Rc: " << _jacNColorEqn << " colours over " << _nCollVar + _nCollInp
             << " columns (the outputs' " << _jacNColorFct << ", separately)\n";
  _shape.n_diff_cols = 0; for( bool b : _colDiff ) _shape.n_diff_cols += b;
  _shape.n_alg_rows  = 0; for( bool b : _rowAlg )  _shape.n_alg_rows  += b;
  _disp( 1 ) << "OCFEMOL::setup [mol] Mc: " << _mass.size() << " entries; " << _shape.n_diff_cols << " differential unknowns, "
             << _shape.n_states - _shape.n_diff_cols << " algebraic; " << _shape.n_alg_rows << " algebraic rows\n";

  // 5b. the structural index of the SEMI-DISCRETE system (the design's step 6, 8 Nov 2026): index 1 iff the algebraic
  //     rows can be matched one-to-one to the algebraic unknowns on Rc's pattern.  A deficiency names its cause -- an
  //     interface multiplier with no algebraic row, a constraint row on differential unknowns only -- which is how a
  //     multiplier minted by the plan (every value-continuity multiplier under IC_TRACE) raises the index to 2.
  _nsActive = false;  _nsN = _nCollVar;  _nsM = _nCollEqn;  _nsTau = _nTraceVar;
  { std::vector<size_t> badCols, badRows;  std::string why;
    if( !_structural_index_check( _disp( 1 ), _jacPatternEqn, _nCollVar, _colDiff, _rowAlg, _nTraceVar, &badCols, &badRows ) ){
      std::string const finding = _setupError;
      if( !_eliminate_multipliers( _disp( 1 ), badCols, badRows, why ) ){
        _setupError = finding + "  Not eliminated: " + why;
        std::cerr << "OCFEMOL::setup ** " << _setupError << std::endl;
        _setupStatus = SetupStatus::INCONSISTENT_MODEL;
        return false;
      }
      if( !_structural_index_check( _disp( 1 ), _jacPatternRed, _nsN, _colDiff, _rowAlg, _nsTau ) ){
        _setupError = "after the null-space elimination, " + _setupError;
        std::cerr << "OCFEMOL::setup ** " << _setupError << std::endl;
        _setupStatus = SetupStatus::INCONSISTENT_MODEL;
        return false;
      }
    } }
  size_t const bad = _audit_initial_coverage( _disp( 1 ) );
  _disp( 1 ) << "OCFEMOL::setup [mol] Ini stream: " << _shape.n_ini_stream << " INITIAL (row, node) pairs at differential nodes; "
             << _shape.n_ini_dropped << " at algebraic nodes dropped (consistency at t0)\n";
  if( bad ){
    if( !_auditRefusals.empty() )
      _setupError = _auditRefusals.front() + ( _auditRefusals.size() > 1? " (and " + std::to_string( _auditRefusals.size() - 1 ) + " more)": "" );
    else
      _setupError = std::to_string( bad ) + " differential unknown(s) without an initial condition (see the [mol] DOF audit)";
    std::cerr << "OCFEMOL::setup ** " << _setupError << std::endl;
    _setupStatus = SetupStatus::INCONSISTENT_MODEL;
    return false;
  }
  return true;
}

inline bool
OCFEMOL::eval_residual
( double t, double const* xc, double const* uc, double* rc )
{
  if( !_shape.ran || !xc || !rc ) return false;
  xc = _expand( xc );                                            // the full state, when multipliers were eliminated
  // the collocated input vector: the caller's values, t in its own slot
  std::vector<double> inp( _nCollInp, 0. );
  if( uc ) std::copy( uc, uc + _nCollInp, inp.begin() );
  if( _nCollInp ) inp[ pos_input( _evolution_dom_var, {} ) ] = t;
  for( size_t i = 0; i < _vCst.size() && i < _vCstVal.size(); ++i ) _OCCst[i].set( _vCstVal[i] );
  if( !_nsActive ){ _eval_eqn( rc, xc, _nCollInp? inp.data(): nullptr, _OCEqn, _OCVar, _OCInp, _OCDom, _OCCst, _wrkOC );  return true; }
  _nsR.assign( _nCollEqn, 0. );
  _eval_eqn( _nsR.data(), xc, _nCollInp? inp.data(): nullptr, _OCEqn, _OCVar, _OCInp, _OCDom, _OCCst, _wrkOC );
  _contract_rows( _nsR.data(), rc );
  return true;
}

inline bool
OCFEMOL::_eliminate_multipliers
( std::ostream& os, std::vector<size_t> const& badCols, std::vector<size_t> const& badRows, std::string& why )
{
  size_t const npos = std::string::npos, n = _nCollVar, m = _nCollEqn;
  // the Jacobian of the FULL system at the reference point: the multipliers' signs and the constraints' coefficients
  std::vector<double> xref( n, 0. ), uref( _nCollInp, 0. ), rc( m );
  for( size_t c = 0; c < n; ++c ){ FFVar f;  double ref = 0.;  if( column_field( c, f, ref ) ) xref[c] = ref; }
  std::vector<size_t> P, C;  std::vector<double> V;
  if( !deriv_residual_full( 0., xref.data(), _nCollInp? uref.data(): nullptr, rc.data(), P, C, V ) ){ why = "the Jacobian at the reference point could not be formed"; return false; }
  auto entry = [&]( size_t r, size_t c )->double{ for( size_t k = P[r]; k < P[r+1]; ++k ) if( C[k] == c ) return V[k];  return 0.; };
  // row -> its differential unknowns: a row may carry SEVERAL mass entries (the PSA's c1 balance carries F d_t q1 too),
  // so the copies a multiplier joins are named by its constraint, not by the rows' masses (8 Nov 2026)
  std::vector< std::set<size_t> > rowMass( m );  for( auto const& e : _mass ) rowMass[e.row].insert( e.col );
  // union-find over columns and rows
  std::vector<size_t> ufc( n ), ufr( m );  for( size_t c = 0; c < n; ++c ) ufc[c] = c;  for( size_t r = 0; r < m; ++r ) ufr[r] = r;
  auto find = []( std::vector<size_t>& uf, size_t i ){ while( uf[i] != i ){ uf[i] = uf[ uf[i] ];  i = uf[i]; }  return i; };
  std::vector<bool> dropCol( n, false ), dropRow( m, false );
  size_t eliminated = 0;
  for( size_t c : badCols ){
    if( c < n - _nTraceVar ){ why = "unknown " + std::to_string( c ) + " is a field node with no algebraic row, not a multiplier"; return false; }
    std::vector<size_t> rows;  for( size_t r = 0; r < m; ++r ) if( entry( r, c ) != 0. ) rows.push_back( r );
    if( rows.size() != 2 ){ why = "multiplier " + std::to_string( c ) + " enters " + std::to_string( rows.size() ) + " rows (two copies of a node expected)"; return false; }
    size_t const L = rows[0], R = rows[1];
    if( _rowAlg[L] || _rowAlg[R] || rowMass[L].empty() || rowMass[R].empty() ){ why = "multiplier " + std::to_string( c ) + " enters a row without a differential unknown"; return false; }
    double const aL = entry( L, c ), aR = entry( R, c );
    if( std::fabs( aL + aR ) > 1e-10 * std::max( std::fabs( aL ), std::fabs( aR ) ) ){ why = "multiplier " + std::to_string( c ) + " does not enter its two rows with opposite equal weights"; return false; }
    // its constraint: an unmatched algebraic row, linear, reading exactly two columns -- one differential unknown of each
    // of the multiplier's rows -- with opposite equal coefficients: those are the copies
    size_t g = npos, cL = npos, cR = npos;
    for( size_t r : badRows ){
      if( dropRow[r] || !_rowAlg[r] || !_jacLinearEqn[r] ) continue;
      std::vector<size_t> cols;  for( size_t cc : _jacPatternEqn[r] ) if( cc < n ) cols.push_back( cc );
      if( cols.size() != 2 ) continue;
      if( rowMass[L].count( cols[0] ) && rowMass[R].count( cols[1] ) ){ g = r;  cL = cols[0];  cR = cols[1];  break; }
      if( rowMass[L].count( cols[1] ) && rowMass[R].count( cols[0] ) ){ g = r;  cL = cols[1];  cR = cols[0];  break; }
    }
    if( g == npos ){ why = "no linear constraint row joining a differential unknown of each of the two rows (" + std::to_string( L ) + ", " + std::to_string( R ) + ") of multiplier " + std::to_string( c ) + " was found"; return false; }
    double const bL = entry( g, cL ), bR = entry( g, cR );
    if( std::fabs( bL + bR ) > 1e-10 * std::max( std::fabs( bL ), std::fabs( bR ) ) ){ why = "constraint row " + std::to_string( g ) + " is not a continuity x_L - x_R = 0 (coefficients not opposite and equal)"; return false; }
    ufc[ find( ufc, cR ) ] = find( ufc, cL );  ufr[ find( ufr, R ) ] = find( ufr, L );
    dropCol[c] = true;  dropRow[g] = true;  ++eliminated;
    FFVar f;  double ref = 0.;  column_field( cL, f, ref );
    os << "OCFEMOL::setup [mol] null-space elimination: multiplier " << c << " and constraint row " << g << " removed, the copies "
       << cL << " and " << cR << " of " << f.name() << " merged, their rows " << L << " and " << R << " summed\n";
  }
  if( !eliminated ){ why = "nothing to eliminate"; return false; }
  // the reduced numbering, in order of first appearance; the kept multipliers stay last
  _nsCol.assign( n, npos );  _nsRow.assign( m, npos );
  std::vector<size_t> numC( n, npos ), numR( m, npos ), repC;  _nsN = 0;  _nsM = 0;
  for( size_t c = 0; c < n; ++c ){ if( dropCol[c] ) continue;  size_t const rt = find( ufc, c );  if( numC[rt] == npos ){ numC[rt] = _nsN++;  repC.push_back( c ); }  _nsCol[c] = numC[rt]; }
  for( size_t r = 0; r < m; ++r ){ if( dropRow[r] ) continue;  size_t const rt = find( ufr, r );  if( numR[rt] == npos ) numR[rt] = _nsM++;  _nsRow[r] = numR[rt]; }
  _nsTau = _nTraceVar - eliminated;
  // the stored structures, in the reduced numbering
  { std::map< std::pair<size_t,size_t>, FFVar > acc;
    for( auto const& e : _mass ){ size_t const rr = _nsRow[e.row], cc = _nsCol[e.col];  if( rr == npos || cc == npos ) continue;
      auto it = acc.find( { rr, cc } );  if( it == acc.end() ) acc.emplace( std::make_pair( rr, cc ), e.coeff );  else it->second = it->second + e.coeff; }
    _mass.clear();  for( auto const& [k, v] : acc ){ t_MassEntry e;  e.row = k.first;  e.col = k.second;  e.coeff = v;  _mass.push_back( e ); } }
  { std::vector<bool> cd( _nsN, false ), ra( _nsM, true );
    for( size_t c = 0; c < n; ++c ) if( _nsCol[c] != npos && _colDiff[c] ) cd[ _nsCol[c] ] = true;
    for( size_t r = 0; r < m; ++r ) if( _nsRow[r] != npos && !_rowAlg[r] ) ra[ _nsRow[r] ] = false;
    _colDiff = cd;  _rowAlg = ra; }
  for( auto& [f, cols] : _fieldCols ) for( size_t& c : cols ) c = _nsCol[c];
  { std::vector<size_t> owner( _nsN, npos );  for( size_t c = 0; c < n; ++c ) if( _nsCol[c] != npos && owner[ _nsCol[c] ] == npos ) owner[ _nsCol[c] ] = _colOwner[c];  _colOwner = owner; }
  { std::vector<t_IniEntry> kept;                                  // (the stream is built AFTER this, by the audit, in the reduced numbering; this covers a rebuilt one)
    for( auto& e : _iniStream ){
      size_t const tgt = e.target;
      if( tgt < n && _nsCol[tgt] != npos && repC[ _nsCol[tgt] ] != tgt ) continue;     // the dependent copy's INITIAL row: one condition per merged unknown
      e.target = tgt < n? _nsCol[tgt]: tgt;
      for( auto& fc : e.cols ) fc.second = _nsCol[ fc.second ];
      kept.push_back( e ); }
    _iniStream = kept; }
  auto contract_pattern = [&]( std::vector<std::vector<size_t>> const& full, std::vector<std::vector<size_t>>& red, bool rows ){
    std::vector< std::set<size_t> > acc( rows? _nsM: full.size() );
    for( size_t r = 0; r < full.size(); ++r ){ size_t const rr = rows? _nsRow[r]: r;  if( rr == npos ) continue;
      for( size_t c : full[r] ){ size_t const cc = _col_red( c );  if( cc != npos ) acc[rr].insert( cc ); } }
    red.assign( acc.size(), {} );  for( size_t r = 0; r < acc.size(); ++r ) red[r].assign( acc[r].begin(), acc[r].end() ); };
  contract_pattern( _jacPatternEqn, _jacPatternRed, true );
  contract_pattern( _jacPatternFct, _jacPatternFctRed, false );
  _shape.n_states = _nsN;  _shape.n_rows = _nsM;  _shape.n_tau = _nsTau;
  _shape.n_diff_cols = 0;  for( bool b : _colDiff ) _shape.n_diff_cols += b;
  _shape.n_alg_rows = 0;   for( bool b : _rowAlg )  _shape.n_alg_rows  += b;
  _nsActive = true;
  os << "OCFEMOL::setup [mol] null-space elimination: " << eliminated << " multiplier(s) and constraint row(s) removed; " << _nsN << " unknowns ("
     << _nsTau << " explicit multipliers), " << _nsM << " rows; " << _shape.n_diff_cols << " differential, " << _nsN - _shape.n_diff_cols << " algebraic\n";
  return true;
}

inline void
OCFEMOL::jacobian_pattern
( std::vector<size_t>& jac_ptr, std::vector<size_t>& jac_col )
const
{
  auto const& pat = _nsActive? _jacPatternRed: _jacPatternEqn;  size_t const nred = n_reduced();
  jac_ptr.assign( 1, 0 ); jac_col.clear();
  for( size_t r = 0; r < pat.size(); ++r ){
    for( size_t c : pat[r] ) if( c < nred ) jac_col.push_back( c );   // the xc columns only
    jac_ptr.push_back( jac_col.size() );
  }
}

inline void
OCFEMOL::jacobian_pattern_full
( std::vector<size_t>& jac_ptr, std::vector<size_t>& jac_col )
const
{
  auto const& pat = _nsActive? _jacPatternRed: _jacPatternEqn;
  size_t const ncol = n_reduced() + _nCollInp;
  jac_ptr.assign( 1, 0 ); jac_col.clear();
  for( size_t r = 0; r < pat.size(); ++r ){
    for( size_t c : pat[r] ) if( c < ncol ) jac_col.push_back( c );
    jac_ptr.push_back( jac_col.size() );
  }
}

inline bool
OCFEMOL::deriv_residual_full
( double t, double const* xc, double const* uc, double* rc, std::vector<size_t>& jac_ptr, std::vector<size_t>& jac_col,
  std::vector<double>& jac_val )
{
  if( !_shape.ran || !_jacReady || !xc ) return false;
  bool const reduce = _nsActive;  xc = _expand( xc );
  size_t const ncol = _nCollVar + _nCollInp;
  // the FAD buffers, seeded by colour exactly as OCFESLV::_reset_autodiff seeds its own -- uc's slots coloured with xc
  std::vector<t_FAD> Fvar( _nCollVar ), Finp( _nCollInp ), Feqn( _nCollEqn );
  for( size_t i = 0; i < _nCollVar; ++i ){ Fvar[i] = xc[i]; if( _jacNColorEqn ) Fvar[i].diff( _jacColorEqn[i], _jacNColorEqn ); }
  size_t const it = _nCollInp? pos_input( _evolution_dom_var, {} ): 0;
  for( size_t i = 0; i < _nCollInp; ++i ){
    Finp[i] = ( i == it? t: ( uc? uc[i]: 0. ) );
    if( _jacNColorEqn ) Finp[i].diff( _jacColorEqn[_nCollVar+i], _jacNColorEqn );
  }
  _reset_colloc( _OCFVar, _OCFInp, _OCFDom, _OCFCst );
  _OCFEqn = OCVar<t_FAD>();  _wrkOCF.clear();
  std::vector<t_FAD> Fcst( _vCstVal.begin(), _vCstVal.end() );          // the model's constants, as OCFESLV passes them
  for( size_t i = 0; i < _vCst.size() && i < Fcst.size(); ++i ) _OCFCst[i].update( &Fcst[i] );
  _eval_eqn( Feqn.data(), Fvar.empty()? nullptr: Fvar.data(), Finp.empty()? nullptr: Finp.data(),
             _OCFEqn, _OCFVar, _OCFInp, _OCFDom, _OCFCst, _wrkOCF );
  jac_ptr.assign( 1, 0 ); jac_col.clear(); jac_val.clear();
  std::vector<double> rfull( reduce? _nCollEqn: 0 );
  for( size_t r = 0; r < _nCollEqn; ++r ){
    if( rc ){ if( reduce ) rfull[r] = Feqn[r].val(); else rc[r] = Feqn[r].val(); }
    for( size_t c : _jacPatternEqn[r] ){
      if( c >= ncol ) continue;
      jac_col.push_back( c );  jac_val.push_back( _jacNColorEqn? Feqn[r].deriv( _jacColorEqn[c] ): 0. );
    }
    jac_ptr.push_back( jac_col.size() );
  }
  if( reduce ){ _contract_csr( jac_ptr, jac_col, jac_val, _nsM );  if( rc ) _contract_rows( rfull.data(), rc ); }
  return true;
}

inline bool
OCFEMOL::eval_outputs
( double t, double const* xc, double const* uc, double* f )
{
  if( !_shape.ran || !xc || !f ) return false;
  if( !_nCollFct ) return true;
  xc = _expand( xc );
  std::vector<double> inp( _nCollInp, 0. );
  if( uc ) std::copy( uc, uc + _nCollInp, inp.begin() );
  if( _nCollInp ) inp[ pos_input( _evolution_dom_var, {} ) ] = t;
  for( size_t i = 0; i < _vCst.size() && i < _vCstVal.size(); ++i ) _OCCst[i].set( _vCstVal[i] );
  _eval_fct( f, xc, _nCollInp? inp.data(): nullptr, _OCEqn, _OCVar, _OCInp, _OCDom, _OCCst, _wrkOC );
  return true;
}

inline bool
OCFEMOL::deriv_outputs_full
( double t, double const* xc, double const* uc, double* f, std::vector<size_t>& ptr, std::vector<size_t>& col,
  std::vector<double>& val )
{
  if( !_shape.ran || !_jacReady || !xc ) return false;
  size_t const ncol = _nCollVar + _nCollInp;
  ptr.assign( 1, 0 ); col.clear(); val.clear();
  if( !_nCollFct ) return true;
  xc = _expand( xc );
  std::vector<t_FAD> Fvar( _nCollVar ), Finp( _nCollInp ), Ffct( _nCollFct );
  for( size_t i = 0; i < _nCollVar; ++i ){ Fvar[i] = xc[i]; if( _jacNColorFct ) Fvar[i].diff( _jacColorFct[i], _jacNColorFct ); }
  size_t const it = _nCollInp? pos_input( _evolution_dom_var, {} ): 0;
  for( size_t i = 0; i < _nCollInp; ++i ){
    Finp[i] = ( i == it? t: ( uc? uc[i]: 0. ) );
    if( _jacNColorFct ) Finp[i].diff( _jacColorFct[_nCollVar+i], _jacNColorFct );
  }
  _reset_colloc( _OCFVar, _OCFInp, _OCFDom, _OCFCst );
  _OCFEqn = OCVar<t_FAD>();  _wrkOCF.clear();
  std::vector<t_FAD> Fcst( _vCstVal.begin(), _vCstVal.end() );
  for( size_t i = 0; i < _vCst.size() && i < Fcst.size(); ++i ) _OCFCst[i].update( &Fcst[i] );
  _eval_fct( Ffct.data(), Fvar.empty()? nullptr: Fvar.data(), Finp.empty()? nullptr: Finp.data(),
             _OCFEqn, _OCFVar, _OCFInp, _OCFDom, _OCFCst, _wrkOCF );
  for( size_t r = 0; r < _nCollFct; ++r ){
    if( f ) f[r] = Ffct[r].val();
    for( size_t c : _jacPatternFct[r] ){
      if( c >= ncol ) continue;
      col.push_back( c );  val.push_back( _jacNColorFct? Ffct[r].deriv( _jacColorFct[c] ): 0. );
    }
    ptr.push_back( col.size() );
  }
  if( _nsActive ){ std::vector<size_t> ident( ptr.size() - 1 );  for( size_t r = 0; r < ident.size(); ++r ) ident[r] = r;  _contract_csr_rows_only( ptr, col, val, ident, ident.size() ); }
  return true;
}

inline bool
OCFEMOL::deriv_residual
( double t, double const* xc, double* rc, std::vector<size_t>& jac_ptr, std::vector<size_t>& jac_col,
  std::vector<double>& jac_val, double* drdt )
{
  // the xc columns of the full derivative, with t's column apart (one FAD code path)
  std::vector<size_t> P, C;  std::vector<double> V;
  if( !deriv_residual_full( t, xc, nullptr, rc, P, C, V ) ) return false;
  size_t const nred = n_reduced(), tc = nred + ( _nCollInp? pos_input( _evolution_dom_var, {} ): 0 );
  jac_ptr.assign( 1, 0 ); jac_col.clear(); jac_val.clear();
  for( size_t r = 0; r + 1 < P.size(); ++r ){
    if( drdt ) drdt[r] = 0.;
    for( size_t k = P[r]; k < P[r+1]; ++k ){
      if( C[k] < nred ){ jac_col.push_back( C[k] ); jac_val.push_back( V[k] ); }
      else if( drdt && _nCollInp && C[k] == tc ) drdt[r] = V[k];
    }
    jac_ptr.push_back( jac_col.size() );
  }
  return true;
}

inline bool
OCFEMOL::eval_initial
( double t, double const* xc, double* ini, double const* uc )
{
  if( !_shape.ran || !xc || !ini ) return false;
  for( size_t k = 0; k < _iniStream.size(); ++k ){
    auto const& e = _iniStream[k];
    std::vector<FFVar> vars; std::vector<double> vals;
    for( auto const& [y, c] : e.cols ){ vars.push_back( y ); vals.push_back( xc[c] ); }
    if( !e.icols.empty() && !uc ) return false;                        // reads an input: its values are needed
    for( auto const& [w, comb] : e.icols ){
      double v = 0.;  for( auto const& [c, a] : comb ) v = ( comb.size() == 1 && a == 1. )? uc[c]: v + a*uc[c];
      vars.push_back( w );  vals.push_back( v );
    }
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
( double t, double const* xc, double* ini, std::vector<size_t>& ptr, std::vector<size_t>& col, std::vector<double>& val,
  double const* uc )
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
    if( !e.icols.empty() && !uc ) return false;                        // reads an input: its values are needed
    for( auto const& [w, comb] : e.icols ){                            // held fixed
      double v = 0.;  for( auto const& [c, a] : comb ) v = ( comb.size() == 1 && a == 1. )? uc[c]: v + a*uc[c];
      vars.push_back( w );  vals.push_back( t_FAD( v ) );
    }
    std::vector<t_FAD> out( 1 );
    try{ _dag->eval( std::vector<FFVar>{ _iniEqn[e.row].var }, out, vars, vals ); } catch( ... ){ return false; }
    if( ini ) ini[k] = out[0].val();
    for( size_t j = 0; j < e.cols.size(); ++j ){ col.push_back( e.cols[j].second ); val.push_back( e.cols.size()? out[0].deriv( j ): 0. ); }
    ptr.push_back( col.size() );
  }
  return true;
}

inline bool
OCFEMOL::init_seed
( double t, double* xc )
{
  if( !_shape.ran || !xc ) return false;
  std::vector<double> inp( _nCollInp, 0. );
  if( !_nsActive ) _init_reference( xc, _nCollInp? inp.data(): nullptr, _vCstVal.empty()? nullptr: _vCstVal.data() );
  else{                                                          // the base fills the FULL layout: take the representatives
    std::vector<double> full( _nCollVar, 0. );
    _init_reference( full.data(), _nCollInp? inp.data(): nullptr, _vCstVal.empty()? nullptr: _vCstVal.data() );
    for( size_t c = 0; c < _nCollVar; ++c ) if( _nsCol[c] != std::string::npos ) xc[ _nsCol[c] ] = full[c];
  }
  std::vector<double> ini( _iniStream.size() );
  std::vector<size_t> ptr, col; std::vector<double> val;
  for( int sweep = 0; sweep < 8; ++sweep ){
    if( !deriv_initial( t, xc, ini.data(), ptr, col, val ) ) return false;
    double move = 0.;
    // every INITIAL row seeds its node -- the dropped ones too: an algebraic node is HELD by its boundary or
    // algebraic row, not by its INITIAL row, but the declared profile is where the consistent initialisation should
    // START.  From the reference value instead, a boundary row with several roots can be met on the wrong one: the
    // flux row D(u) u_z = q with D = D0 (1 + beta u) was met at 1 + beta u < 0 on the heated rod of the MoL tutorial,
    // a negative diffusivity, u(2,1) = -1.02 against 2.79 (7 Nov 2026).
    for( size_t k = 0; k < _iniStream.size(); ++k ){
      size_t const tgt = _iniStream[k].target;
      double d = 0.;  for( size_t j = ptr[k]; j < ptr[k+1]; ++j ) if( col[j] == tgt ) d = val[j];
      if( std::fabs( d ) < 1e-12 ) continue;
      double const step = -ini[k] / d;
      xc[tgt] += step;  move = std::max( move, std::fabs( step ) );
    }
    if( move < 1e-14 ) break;
  }
  return true;
}

inline double
OCFEMOL::initial_consistency_defect
( double t, double const* xc )
{
  if( !_shape.ran || !xc ) return 0.;
  std::vector<double> ini( _iniStream.size() );
  if( !eval_initial( t, xc, ini.data() ) ) return 1e30;
  double worst = 0.;
  for( size_t k = 0; k < _iniStream.size(); ++k ) if( _iniStream[k].dropped ) worst = std::max( worst, std::fabs( ini[k] ) );
  return worst;
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
OCFEMOL::_build_column_map
()
{
  _fieldCols.clear();  _fieldNodes.clear();
  _colOwner.assign( _shape.n_states, std::string::npos == 0? 0: (size_t)-1 );
  for( size_t i = 0; i < _peel.states.size(); ++i ){
    FFVar const& x = _peel.states[i];
    auto itv = _mVar.find( x );  auto itb = _blkVar.find( x );
    if( itv == _mVar.end() || itb == _blkVar.end() ) continue;
    auto const& var_dom = itv->second;  auto const& blk_dim = std::get<2>( itb->second );
    std::vector<size_t>& cols = _fieldCols[x];
    auto& nodes = _fieldNodes[x];
    if( var_dom.empty() ){                                        // lumped: one column, no coordinates
      size_t const c = pos_state( x, {} );
      cols.push_back( c );  nodes.emplace_back();  if( c < _colOwner.size() ) _colOwner[c] = i;
      continue;
    }
    std::vector<FFVar> dirs( var_dom.begin(), var_dom.end() );
    std::vector<size_t> nel( dirs.size() ), nnd( dirs.size() ), iel( dirs.size(), 0 ), ind( dirs.size(), 0 );
    for( size_t d = 0; d < dirs.size(); ++d ){ nel[d] = _mDom.at( dirs[d] ).n_elem; nnd[d] = blk_dim[d]; }
    while( true ){
      std::map<FFVar,size_t,lt_FFVar> ndx_el;  for( size_t d = 0; d < dirs.size(); ++d ) ndx_el[dirs[d]] = iel[d];
      size_t c = pos_state( x, ndx_el ), stride = 1;
      for( size_t d = 0; d < dirs.size(); ++d ){ c += ind[d]*stride; stride *= nnd[d]; }
      cols.push_back( c );  if( c < _colOwner.size() ) _colOwner[c] = i;
      nodes.emplace_back();
      for( size_t d = 0; d < dirs.size(); ++d ){
        FFDom const& dm = _mDom.at( dirs[d] );
        double const lo = dm.elem_lo( iel[d] ), hi = dm.elem_up( iel[d] ), xi = dm.lgnodes().at( ind[d] );
        nodes.back().emplace_back( dirs[d], 0.5*( lo + hi ) + 0.5*( hi - lo )*xi );
      }
      // mesh order: the LAST direction varies fastest within an element, elements likewise -- the order a driver
      // writing  for( e ) for( i )  produces, direction by direction
      size_t d = dirs.size();  bool done = false;
      while( d-- > 0 ){
        if( ++ind[d] < nnd[d] ) break;
        ind[d] = 0;
        if( ++iel[d] < nel[d] ) break;
        iel[d] = 0;
        if( d == 0 ) done = true;                 // this FIELD's grid is exhausted -- not the whole map
      }
      if( done ) break;
    }
  }
}

inline bool
OCFEMOL::_split_output_time
( FFVar const& t )
{
  _fctTime.clear();  _fctSpatial.clear();  _nFctDeclared = _mFct.size();
  FFIntegral OpI;  FFEval OpE;
  auto fail = [&]( size_t k, std::string const& why ){
    _setupError = "output " + std::to_string( k ) + ": " + why;
    std::cerr << "OCFEMOL::setup ** " << _setupError << std::endl;
    return false; };
  auto has_t_operator = [&]( FFVar const& E ) -> bool {             // a t-operator anywhere in E
    FFSubgraph const sg = _dag->subgraph( 1, &E );
    for( auto const* op : sg.l_op ){
      if( auto const* pi = mc::type_cast<FFIntegral const>( op ) ){ if( pi->Indep().expr.count( t ) ) return true; }
      else if( auto const* pe = mc::type_cast<FFEval const>( op ) ){ if( pe->Coord().count( t ) ) return true; }
      else if( auto const* pp = mc::type_cast<FFPartial const>( op ) ){ if( pp->Indep().expr.count( t ) ) return true; }
    }
    return false; };
  for( size_t k = 0; k < _mFct.size(); ++k ){
    if( _fctSpatial.count( k ) ){                              // a SPATIAL part: its rows only, no time of its own
      t_OutputTime T;  T.kind = t_OutputTime::SPATIAL;  _fctTime.push_back( T );  continue; }
    auto& f = _mFct[k];
    t_OutputTime T;  bool timed = false;
    if( f.grid.count( t ) ){
      // a TRAJECTORY (7 Nov 2026): its spatial part at every stage time its mask selects -- as the integrators take an
      // output distributed over the evolution direction; distributed over space too, one trajectory per spatial row
      if( f.point.count( t ) ) return fail( k, "it is both distributed over t and given a point in t" );
      FFVar const E = f.var;
      if( has_t_operator( E ) )
        return fail( k, "it is distributed over t and holds an operator in t -- a trajectory is a value at each stage time" );
      T.kind = t_OutputTime::TRAJECTORY;  T.tmask = f.grid.at( t );
      f.grid.erase( t );
      if( f.grid.empty() ) f.kind = FctKind::POINT;                // a scalar in space at every stage time
      _fctTime.push_back( T );
      continue;
    }
    // a point in t given with the output
    if( auto itp = f.point.find( t ); itp != f.point.end() ){
      T.kind = t_OutputTime::POINT;  T.tau = itp->second;
      if( auto its = f.side.find( t ); its != f.side.end() ){ T.side = its->second;  f.side.erase( its ); }
      f.point.erase( itp );  timed = true;
    }
    FFVar E = f.var;
    // a DEFERRED VALUE: FFModel's setup already took the reduction in t out of the output -- the output reads the
    // deferred input, and the record holds the reduction's source (the integrand, or the operand at tau) and its kind.
    // That is exactly the split here: the source is the spatial part, the kind the time part.
    if( is_deferred_input( E ) ){
      if( timed ) return fail( k, "it is a reduction in t (a deferred value) and is also given a time" );
      t_DeferredValue const* C = nullptr;
      for( auto const& c : var_deferred() ) if( c.input.id() == E.id() ){ C = &c; break; }
      if( !C ) return fail( k, "its deferred value has no record" );
      T.kind = C->accumulate? t_OutputTime::INTEGRAL: t_OutputTime::POINT;
      if( !C->accumulate ){ T.tau = C->tau;  if( C->side ) T.side = C->side; }
      E = C->source;  timed = true;
    }
    else{
      // a FUNCTION of deferred values: every deferred input it reads becomes a spatial output of its own (appended,
      // classified by the branch above on a later iteration), and the declaration composes the function over their
      // reductions -- which BASE_DE takes since 7 Nov 2026.  Its other independents must be constants.
      // the deferred inputs it reads (FFModel's captures), and the OUTERMOST reductions in t INSIDE it (not at its root:
      // that case is the outermost-operator branch below) -- under the semi-discretisation FFModel's capture leaves an
      // outer reduction whose operand holds a spatial one, so the extraction is done here
      std::vector<FFVar> dv;
      FFSubgraph const sgE = _dag->subgraph( 1, &E );
      for( auto const* v : sgE.v_indep ) if( v && is_deferred_input( *v ) ) dv.push_back( *v );
      auto is_t_reduction = [&]( FFOp const* op )->bool{
        if( auto const* pi = mc::type_cast<FFIntegral const>( op ) ) return pi->Indep().expr.count( t ) > 0;
        if( auto const* pe = mc::type_cast<FFEval const>( op ) )     return pe->Coord().count( t ) > 0;
        return false; };
      std::vector<FFVar> red;
      { FFOp const* root = E.opdef().first;
        std::vector<FFVar> cand;
        for( auto const* op : sgE.l_op ) if( op && op != root && is_t_reduction( op ) ) for( auto const* vo : op->varout ) if( vo ) cand.push_back( *vo );
        for( size_t i = 0; i < cand.size(); ++i ){                   // outermost only: not inside another candidate's operand
          bool inner = false;
          for( size_t j = 0; j < cand.size() && !inner; ++j ){
            if( i == j ) continue;
            for( auto const* vo : _dag->subgraph( 1, &cand[j] ).l_op )
              for( auto const* v : vo->varout ) if( v && v->id() == cand[i].id() && vo->varout[0]->id() != cand[j].id() ){ inner = true; break; }
          }
          if( !inner ) red.push_back( cand[i] );
        } }
      if( !dv.empty() || !red.empty() ){
        if( timed ) return fail( k, "it is a function of reductions in t and is also given a time" );
        std::vector<FFVar> ph;                                        // placeholders for the inner reductions in h
        for( size_t i = 0; i < red.size(); ++i ) ph.push_back( _dag->add_var( "fct" + std::to_string( k ) + "_red" + std::to_string( i ) ) );
        FFVar h = E;
        if( !red.empty() ) h = _dag->substitute( std::vector<FFVar>{ E }, red, ph )[0];   // node-level, as FFModel replaces its captures
        // what else h reads must be a constant: a state or an input read directly has no time here
        for( auto const* v : _dag->subgraph( 1, &h ).v_indep ){
          if( !v || is_deferred_input( *v ) ) continue;
          bool ok = false;
          for( auto const& c : _vCst ) if( c.id() == v->id() ){ ok = true; break; }
          for( auto const& q : ph )    if( q.id() == v->id() ){ ok = true; break; }
          if( !ok ) return fail( k, "it is a function of reductions in t and also reads " + v->name()
                                    + " directly -- a function of reductions may read only the reductions and constants" );
        }
        T.kind = t_OutputTime::COMPOSED;  T.h = h;
        T.dvars = dv;  T.dvars.insert( T.dvars.end(), ph.begin(), ph.end() );
        t_Fct const base = _mFct[k];     // copied BEFORE any push: f references _mFct, which the pushes may reallocate
        for( auto const& d : dv ){                                    // a captured value: its record's source, by the branch above
          t_Fct part = base;  part.var = d;  part.point.clear();  part.side.clear();
          T.parts.push_back( _mFct.size() );  _mFct.push_back( part );
        }
        for( auto const& r : red ){                                   // a reduction: the outermost-operator branch, on its own iteration
          t_Fct part = base;  part.var = r;  part.point.clear();  part.side.clear();
          T.parts.push_back( _mFct.size() );  _mFct.push_back( part );
        }
        // its own spatial part is never read: it borrows its first part's, once that part is classified (below)
        _fctTime.push_back( T );
        continue;
      }
    }
    // an outermost integral or evaluation in t
    FFOp const* pOp = E.opdef().first;
    if( pOp && pOp->varin.size() == 1 && pOp->varin[0] ){
      FFVar const operand = *pOp->varin[0];
      if( auto const* pi = mc::type_cast<FFIntegral const>( pOp ); pi && pi->Indep().expr.count( t ) ){
        if( timed ) return fail( k, "it integrates over t and is also given a time" );
        FFIntegral::t_SMon rest;
        for( auto const& [v, o] : pi->Indep().expr ) if( v.id() != t.id() ) rest += FFIntegral::t_SMon( v, o );
        E = rest.tord? OpI( operand, rest ): operand;
        T.kind = t_OutputTime::INTEGRAL;  timed = true;
      }
      else if( auto const* pe = mc::type_cast<FFEval const>( pOp ); pe && pe->Coord().count( t ) ){
        if( timed ) return fail( k, "it is evaluated at a time twice" );
        T.kind = t_OutputTime::POINT;  T.tau = pe->Coord().at( t );
        if( auto its = pe->Side().find( t ); its != pe->Side().end() ) T.side = its->second;
        std::map<FFVar,double,lt_FFVar> coord;  std::map<FFVar,int,lt_FFVar> side;  FFIntegral::t_SMon rest;
        for( auto const& [v, z] : pe->Coord() ) if( v.id() != t.id() ){ coord[v] = z;  rest += FFIntegral::t_SMon( v, 1 ); }
        for( auto const& [v, sd] : pe->Side() )  if( v.id() != t.id() ) side[v] = sd;
        E = rest.tord? OpE( std::vector<FFVar>{ operand }, rest, coord, side )[0]: operand;
        timed = true;
      }
    }
    if( !timed ) return fail( k, "it gives no time -- evaluate it at a point in t, or integrate it over t" );
    if( has_t_operator( E ) )
      return fail( k, "it holds an operator in t below its outermost one -- not taken by the method of lines" );
    // the spatial part: a reduction nested inside a nonlinearity becomes a scalar part of its own (see t_OutputTime)
    { FFOp const* root = E.opdef().first;
      auto is_space_reduction = [&]( FFOp const* op )->bool{
        if( auto const* pi = mc::type_cast<FFIntegral const>( op ) ) return !pi->Indep().expr.empty() && !pi->Indep().expr.count( t );
        if( auto const* pe = mc::type_cast<FFEval const>( op ) )     return !pe->Coord().empty() && !pe->Coord().count( t );
        return false; };
      if( root && !is_space_reduction( root ) ){
        std::vector<FFVar> cand;
        for( auto const* op : _dag->subgraph( 1, &E ).l_op ) if( op && is_space_reduction( op ) ) for( auto const* vo : op->varout ) if( vo ) cand.push_back( *vo );
        std::vector<FFVar> red;
        for( size_t i = 0; i < cand.size(); ++i ){                   // outermost only
          bool inner = false;
          for( size_t j = 0; j < cand.size() && !inner; ++j ){
            if( i == j ) continue;
            for( auto const* vo : _dag->subgraph( 1, &cand[j] ).l_op )
              for( auto const* v : vo->varout ) if( v && v->id() == cand[i].id() && vo->varout[0]->id() != cand[j].id() ){ inner = true; break; }
          }
          if( !inner ) red.push_back( cand[i] );
        }
        if( !red.empty() ){
          for( size_t i = 0; i < red.size(); ++i ) T.sph.push_back( _dag->add_var( "fct" + std::to_string( k ) + "_sred" + std::to_string( i ) ) );
          T.hs = _dag->substitute( std::vector<FFVar>{ E }, red, T.sph )[0];
          t_Fct const base = _mFct[k];   // copied BEFORE any push: f references _mFct, which the pushes may reallocate
          for( auto const& r : red ){
            t_Fct part = base;  part.var = r;  part.point.clear();  part.side.clear();
            _fctSpatial.insert( _mFct.size() );
            T.sparts.push_back( _mFct.size() );  _mFct.push_back( part );   // a reduction at its root: its rows summed over the elements
          }
          E = red[0];                                                 // a valid spatial part of its own; its rows are not read
        }
      } }
    _mFct[k].var = E;                    // by index: f may dangle after the pushes above
    _fctTime.push_back( T );
  }
  for( size_t k = 0; k < _fctTime.size() && k < _mFct.size(); ++k )     // the composed outputs: a valid spatial part of their own
    if( _fctTime[k].kind == t_OutputTime::COMPOSED && !_fctTime[k].parts.empty() ) _mFct[k].var = _mFct[ _fctTime[k].parts[0] ].var;
  return true;
}

inline void
OCFEMOL::_build_input_layout
()
{
  _inpLayout.clear();
  for( auto const& [w, dom] : _mInp ){
    if( _evolution_dom_set && w.id() == _evolution_dom_var.id() ) continue;   // t: time_slot()
    if( is_deferred_input( w ) ) continue;        // a deferred value: the outputs read its source instead (_split_output_time)
    auto const itb = _blkInp.find( w );
    if( itb == _blkInp.end() ) continue;
    t_InputLayout L;  L.input = w;  L.time_varying = _inpTimeVarying.count( w ) > 0;
    std::vector<FFVar> const dirs( dom.begin(), dom.end() );
    if( dirs.empty() ){ L.cols.push_back( pos_input( w, {} ) );  L.nodes.emplace_back();  L.dof.emplace_back();  _inpLayout.push_back( std::move( L ) );  continue; }
    size_t const nd = dirs.size();
    auto const& blk_dim = std::get<2>( itb->second );
    std::vector<size_t> nel( nd ), nnd( nd );
    std::vector< std::vector< std::vector<double> > > pos( nd );        // per direction, per element: the node positions
    for( size_t d = 0; d < nd; ++d ){
      FFDom const& dd = _input_dom( w, dirs[d] );                       // the input's own collocation
      nel[d] = dd.n_elem;  nnd[d] = blk_dim[d];
      for( size_t e = 0; e < dd.n_elem; ++e ) pos[d].push_back( dd.lgnodes( dd.elem_lo( e ), dd.elem_up( e ) ) );
    }
    // input_dofs() order: element blocks outer, nodes inner, at each level the first direction fastest
    std::vector<size_t> iel( nd, 0 );
    for( bool more_el = true; more_el; ){
      std::map<FFVar,size_t,lt_FFVar> ndx_el;  for( size_t d = 0; d < nd; ++d ) ndx_el[ dirs[d] ] = iel[d];
      size_t const base = pos_input( w, ndx_el );
      std::vector<size_t> ind( nd, 0 );
      for( bool more_nd = true; more_nd; ){
        size_t c = base, stride = 1;
        for( size_t d = 0; d < nd; ++d ){ c += ind[d]*stride;  stride *= nnd[d]; }   // within the element: first fastest
        L.cols.push_back( c );  L.nodes.emplace_back();  L.dof.emplace_back();
        for( size_t d = 0; d < nd; ++d ){
          L.nodes.back().emplace_back( dirs[d], pos[d][ iel[d] ].at( ind[d] ) );
          L.dof.back().element[ dirs[d] ] = iel[d];  L.dof.back().node[ dirs[d] ] = ind[d];
        }
        more_nd = false;  for( size_t d = 0; d < nd; ++d ){ if( ++ind[d] < nnd[d] ){ more_nd = true; break; }  ind[d] = 0; }
      }
      more_el = false;  for( size_t d = 0; d < nd; ++d ){ if( ++iel[d] < nel[d] ){ more_el = true; break; }  iel[d] = 0; }
    }
    _inpLayout.push_back( std::move( L ) );
  }
}

inline void
OCFEMOL::_build_mass
()
{
  // Mc receives the plan's row operations exactly as Rc does: an emitted IC_STRONG row is its own row minus a
  // combination of pivot rows the exact elimination dropped, so its mass is its own minus the same combination of the
  // pivots' masses.  A pivot that is the OTHER duplicate of the same field at the same physical point (a continuity
  // claim) has u_pivot = u_receiver by the continuity row, hence u_t too: its folded mass goes on the receiver's own
  // column -- u_pivot is eliminated through its constraint and stays algebraic.  Otherwise it goes on the pivot's column.
  _mass.clear();
  _colDiff.assign( _shape.n_states, false );
  _rowAlg.assign( _shape.n_rows, true );
  auto coord = [&]( t_RowNode const& rn ){                 // the physical coordinates of a row's node
    std::vector<std::pair<long,double>> c;
    for( auto const& d : rn.dims ){
      FFDom const* dm = nullptr;  for( auto const& [v, dd] : _mDom ) if( v.id().second == d.dir_id ){ dm = &dd; break; }
      if( !dm ) continue;
      double const lo = dm->elem_lo( d.iel ), hi = dm->elem_up( d.iel );
      double const xi = d.lim == FFDom::LB? -1.: d.lim == FFDom::UB? 1.: dm->lgnodes().at( d.node );
      c.emplace_back( d.dir_id, 0.5*( lo + hi ) + 0.5*( hi - lo )*xi );
    }
    return c;
  };
  auto same_point = [&]( t_RowNode const& a, t_RowNode const& b ){
    auto ca = coord( a ), cb = coord( b );
    if( ca.size() != cb.size() ) return false;
    for( size_t i = 0; i < ca.size(); ++i ) if( ca[i].first != cb[i].first || std::fabs( ca[i].second - cb[i].second ) > 1e-12 ) return false;
    return true;
  };
  size_t const nrow = std::min( _rowEqnNdx.size(), _shape.n_rows );
  // which columns carry their OWN mass row (before any fold): a pivot column without one is algebraic, pinned to the
  // receiver by its continuity row; a pivot column with one (a field with no spatial operator, hence no continuity
  // claim -- PSA's loadings) is an independent differential unknown and keeps its folded mass
  std::vector<bool> own_diff( _shape.n_states, false );
  for( size_t r = 0; r < nrow; ++r ){
    int const k = _rowEqnNdx[r];
    if( k < 0 || (size_t)k >= _massByEqn.size() ) continue;
    for( auto const& [i, coeff] : _massByEqn[k] ){ size_t col = 0; if( _column_of( _peel.states[i], _rowNode[r], col ) && col < own_diff.size() ) own_diff[col] = true; }
  }
  size_t n_folded = 0;
  for( size_t r = 0; r < nrow; ++r ){
    int const k = _rowEqnNdx[r];
    std::map<size_t,FFVar> row;                                     // column -> coefficient (merged)
    // merge: two CONSTANT coefficients are added as numbers -- the sum of two constant FFVars is a DAG node, not a constant
    auto add = [&]( size_t col, FFVar const& c ){
      auto it = row.find( col );
      if( it == row.end() ){ row.emplace( col, c ); return; }
      if( it->second.cst() && c.cst() ) it->second = FFVar( it->second.num().val() + c.num().val() );   // a plain constant: no DAG needed
      else                              it->second = it->second + c; };
    if( k >= 0 && (size_t)k < _massByEqn.size() )
      for( auto const& [i, coeff] : _massByEqn[k] ){ size_t col = 0; if( _column_of( _peel.states[i], _rowNode[r], col ) ) add( col, coeff ); }
    if( r < _rowFold.size() )
      for( auto const& f : _rowFold[r] ){
        if( f.eqn < 0 || (size_t)f.eqn >= _massByEqn.size() ) continue;
        for( auto const& [i, coeff] : _massByEqn[f.eqn] ){
          FFVar fc;  if( coeff.cst() ) fc = FFVar( -f.coeff * coeff.num().val() );  else fc = -f.coeff * coeff;
          size_t col = 0, pcol = 0;
          bool const has_pcol = _column_of( _peel.states[i], f.node, pcol );
          bool const onto_receiver = k >= 0 && has_pcol && pcol < own_diff.size() && !own_diff[pcol]
                                  && same_point( f.node, _rowNode[r] ) && _column_of( _peel.states[i], _rowNode[r], col );
          if( onto_receiver ){ add( col, fc ); ++n_folded; }
          else if( has_pcol ){ add( pcol, fc ); ++n_folded; }
        }
      }
    for( auto const& [col, c] : row ){
      if( c.cst() && c.num().val() == 0. ) continue;
      _mass.push_back({ r, col, c });
      if( col < _colDiff.size() ) _colDiff[col] = true;
      _rowAlg[r] = false;
    }
  }
  if( n_folded ) _disp( 1 ) << "OCFEMOL::setup [mol] Mc: " << n_folded << " mass entr" << ( n_folded == 1? "y": "ies" )
                            << " folded from rows the exact interface elimination dropped\n";
}

inline size_t
OCFEMOL::_audit_initial_coverage
( std::ostream& os )
{
  _shape.n_ini_stream = _shape.n_ini_dropped = 0;
  size_t bad = 0;
  _auditRefusals.clear();
  auto refuse = [&]( std::string const& m ){ os << "OCFEMOL::setup [mol] DOF audit: " << m << "\n";  _auditRefusals.push_back( m );  ++bad; };
  // An INITIAL row is evaluated at a node of the field it initialises, so an input it reads is taken at THAT node: a
  // lumped input is its one slot; a distributed one, declared over the row's directions (otherwise the model is
  // misspecified: no value at the row's nodes), is its Lagrange interpolant on the element -- in each direction the
  // weights of the input's own nodes at the field node's reference coordinate (FFDom::w_eval_lagrange), their tensor
  // product over the directions, so a linear combination of the input's slots.  This is the input as OCFESLV's
  // collocation evaluates it; at the domain's collocation it is one slot with weight 1, bound as the slot itself.
  std::set< std::pair<size_t,long> > refused_pairs;      // (INITIAL row, input) already refused
  auto over_row_dirs = [&]( FFVar const& x, FFVar const& w, std::set<FFVar,lt_FFVar> const& wdom,
                            std::set<FFVar,lt_FFVar> const& xdom ) -> bool {
    if( wdom == xdom ) return true;
    std::ostringstream m;  m << "the INITIAL row of " << x << " reads the input " << w << ", declared over other directions"
                              << " than the row's: it has no value at the row's nodes";
    refuse( m.str() );  return false; };
  // the combination of w's slots giving its value at the field node ( ndx_el, ind ), w over the field's directions dirs
  auto at_node = [&]( FFVar const& w, std::vector<FFVar> const& dirs, std::map<FFVar,size_t,lt_FFVar> const& ndx_el,
                      std::vector<size_t> const& ind, std::vector< std::pair<size_t,double> >& comb ) -> bool {
    comb.clear();
    size_t const nd = dirs.size();
    std::vector< std::vector<double> > wts( nd );
    std::vector<size_t> nni( nd );
    for( size_t d = 0; d < nd; ++d ){
      FFDom const& di = _input_dom( w, dirs[d] );  FFDom const& dx = _mDom.at( dirs[d] );
      nni[d] = di.n_node;
      if( di.type == dx.type && di.n_node == dx.n_node ){             // the field's own nodes: the unit weight
        wts[d].assign( di.n_node, 0. );  wts[d].at( ind[d] ) = 1.;  continue; }
      if( !di.w_eval_lagrange( wts[d], dx.lgnodes().at( ind[d] ), 0 ) || wts[d].size() != di.n_node ) return false;
    }
    size_t const base = pos_input( w, ndx_el );
    std::vector<size_t> j( nd, 0 );
    for( bool more = true; more; ){
      double a = 1.;  size_t c = base, stride = 1;
      for( size_t d = 0; d < nd; ++d ){ a *= wts[d][ j[d] ];  c += j[d]*stride;  stride *= nni[d]; }
      if( a != 0. ) comb.emplace_back( c, a );
      more = false;  for( size_t d = 0; d < nd; ++d ){ if( ++j[d] < nni[d] ){ more = true; break; }  j[d] = 0; }
    }
    return !comb.empty(); };
  // the INITIAL rows reading each field, with their spatial masks
  std::vector<FFDep> ini_dep;  ini_dep.reserve( _iniEqn.size() );
  for( auto const& e : _iniEqn ) ini_dep.push_back( _deps_of( e.var ) );

  for( size_t i = 0; i < _peel.states.size(); ++i ){
    FFVar const& x = _peel.states[i];
    auto itv = _mVar.find( x );  auto itb = _blkVar.find( x );
    if( itv == _mVar.end() || itb == _blkVar.end() ) continue;
    auto const& var_dom = itv->second;  auto const& blk_dim = std::get<2>( itb->second );
    if( var_dom.empty() ){                            // a lumped state (no spatial direction): one node
      size_t const colFull0 = pos_state( x, {} ), col = _nsActive? _nsCol[ colFull0 ]: colFull0;
      bool const diff = col < _colDiff.size() && _colDiff[col];
      bool covered = false;
      for( size_t q = 0; q < _iniEqn.size() && !covered; ++q ) covered = ini_dep[q].dep().count( x.id().second ) > 0;
      if( covered ){
        if( diff ) ++_shape.n_ini_stream; else ++_shape.n_ini_dropped;
        for( size_t q = 0; q < _iniEqn.size(); ++q ){
          if( !ini_dep[q].dep().count( x.id().second ) ) continue;
          t_IniEntry e; e.row = q; e.target = col; e.dropped = !diff;
          for( auto const& y : _peel.states ) if( ini_dep[q].dep().count( y.id().second ) && _mVar.count( y ) && _mVar.at( y ).empty() )
            { size_t const cy = pos_state( y, {} );  e.cols.emplace_back( y, _nsActive? _nsCol[ cy ]: cy ); }
          for( auto const& [w, wdom] : _mInp ){                         // the lumped inputs it reads
            if( ( _evolution_dom_set && w.id() == _evolution_dom_var.id() ) || !ini_dep[q].dep().count( w.id().second ) ) continue;
            if( wdom.empty() ) e.icols.emplace_back( w, std::vector< std::pair<size_t,double> >{ { pos_input( w, {} ), 1. } } );
            else over_row_dirs( x, w, wdom, var_dom );                 // a distributed input: no value without a point
          }
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
      size_t colFull = pos_state( x, ndx_el ), stride = 1;
      for( size_t d = 0; d < dirs.size(); ++d ){ colFull += ind[d] * stride; stride *= nnd[d]; }
      size_t const col = _nsActive? _nsCol[ colFull ]: colFull;            // the reduced column, when multipliers were eliminated
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
            size_t const cy = pos_state( y, ndx_el ) + ( colFull - pos_state( x, ndx_el ) );
            e.cols.emplace_back( y, _nsActive? _nsCol[ cy ]: cy );
          }
          // the inputs it reads, at the node (see over_row_dirs / at_node).  Refusals are reported once per (row, input).
          for( auto const& [w, wdom] : _mInp ){
            if( ( _evolution_dom_set && w.id() == _evolution_dom_var.id() ) || !ini_dep[q].dep().count( w.id().second ) ) continue;
            if( wdom.empty() ){ e.icols.emplace_back( w, std::vector< std::pair<size_t,double> >{ { pos_input( w, {} ), 1. } } );  continue; }
            auto const key = std::make_pair( q, w.id().second );
            if( refused_pairs.count( key ) ) continue;
            std::vector< std::pair<size_t,double> > comb;
            if( !over_row_dirs( x, w, wdom, var_dom ) ){ refused_pairs.insert( key );  continue; }
            if( !at_node( w, dirs, ndx_el, ind, comb ) ){
              std::ostringstream m;  m << "the INITIAL row of " << x << " reads the input " << w << ", whose value at the field's"
                                        << " nodes could not be interpolated";
              refuse( m.str() );  refused_pairs.insert( key );  continue; }
            e.icols.emplace_back( w, std::move( comb ) );
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
  if( _nsActive ){                                               // merged copies: one INITIAL condition per merged unknown
    std::set<size_t> seen;  std::vector<t_IniEntry> kept;  size_t dup = 0;
    for( auto& e : _iniStream ){ if( !seen.insert( e.target ).second ){ ++dup;  if( e.dropped ) --_shape.n_ini_dropped; else --_shape.n_ini_stream;  continue; }  kept.push_back( std::move( e ) ); }
    _iniStream = std::move( kept );
    if( dup ) os << "OCFEMOL::setup [mol] Ini stream: " << dup << " INITIAL condition(s) at merged copies dropped (one per merged unknown)\n";
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
