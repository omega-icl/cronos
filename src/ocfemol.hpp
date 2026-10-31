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
    size_t n_states    = 0;   //!< collocated unknowns xc: field values at nodes (+ tau multipliers, n_tau)
    size_t n_tau       = 0;   //!< interface multipliers minted by the plan (algebraic)
    size_t n_rows      = 0;   //!< rows of Rc after imposition
    size_t n_ini       = 0;   //!< INITIAL rows of the model (the Ini stream), before collocation
    size_t n_fields    = 0;   //!< fields (states of the model, auxiliaries included)
    size_t n_diff_fields = 0; //!< fields with a d_t somewhere: the differential ones
    double t_lo = 0., t_hi = 0.;  //!< the evolution span, from the model's t domain
    //! @brief The model's character ALONG t, from the peel: the highest spatial order among the differential rows.
    //! PARABOLIC (second order: a diffusion operator), HYPERBOLIC (first order: transport), ODE (no spatial
    //! derivative in a differential row), STEADY (no differential row).  This is what the base's classification
    //! does NOT say: it classifies the spatial operator at fixed t -- ELLIPTIC for the heat equation -- which is
    //! exactly what the plan collocates, and a wrong description of the model.
    enum EvolutionType { STEADY = 0, ODE, HYPERBOLIC, PARABOLIC } evolution = STEADY;
    bool   ran = false;
  };
  static char const* evolution_type_name( t_SemiDiscrete::EvolutionType e )
    { static char const* N[] = { "STEADY", "ODE", "HYPERBOLIC", "PARABOLIC" }; return N[(int)e]; }

  //! @brief The peel this setup took along t, on the working (order-reduced) model.
  t_Peel const& peel() const { return _peel; }
  //! @brief The semi-discrete shape (valid after setup).
  t_SemiDiscrete const& semi_discrete() const { return _shape; }
  //! @brief The INITIAL rows set aside as the Ini stream (working model, order-reduced), by equation index.
  std::vector<t_Eqn> const& ini_rows() const { return _iniEqn; }

protected:

  t_Peel                _peel;      //!< the peel along t (working model)
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

  //! @brief The highest spatial derivative order in @p row, on the order-reduced model: a Partial along a spatial
  //! direction whose operand is an order-reduction auxiliary (itself a first derivative) counts 2, else 1.
  size_t _spatial_order
    ( FFVar const& row )
    const
    {
      size_t order = 0;
      auto sg = _dag->subgraph( 1, &row );
      for( auto const* op : sg.l_op ){
        if( !op || op->type < FFOp::EXTERN || op->varin.size() != 1 || !op->varin[0] ) continue;
        std::string const n = op->name();
        if( !( n.size() > 9 && n.compare( 0, 8, "Partial[" ) == 0 ) ) continue;
        if( n == "Partial[" + _evolution_dom_var.name() + "]" ) continue;
        size_t o = 1;
        for( auto const& a : _auxDef ) if( a.aux.id() == op->varin[0]->id() ){ o = 2; break; }
        order = std::max( order, o );
      }
      return order;
    }

  //! @brief Section 2.1 step 5: every node of every field has exactly one row at t = 0 -- an INITIAL row for a
  //! differential field, an algebraic row of Rc for an algebraic one.  First version, on the model's row MASKS
  //! (the plan's node-level row assignment refines it in phase 3's imposition step).  Returns the number of
  //! fields whose node coverage by INITIAL rows does not match their differential status; each is named.
  size_t _audit_initial_coverage
    ( std::ostream& os )
    const;
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

  // 1. the peel on the WORKING model (order-reduced): roles are resolved by now; the index hook needs the
  //    classification's blocks and runs later
  if( !_peel_evolution( t, _peel, false, false ) ){
    std::cerr << "OCFEMOL::setup ** the model cannot be marched along " << t << ": " << _peel.blocker << std::endl;
    _setupStatus = SetupStatus::INCONSISTENT_MODEL;
    return false;
  }

  // 2. the working model without t.  Row by row: an INTERIOR or BOUNDARY row becomes its remainder R (its d_t
  //    terms are Mc's, not Rc's); an INITIAL row leaves for the Ini stream; LINK, INTERFACE, SURFACE and
  //    DIAGNOSTIC rows stay as they are.  Every row loses its t mask.  _evoDerivRows -- the rows _defining_rows
  //    must not read as definitions -- are the rows with a mass entry, in the new numbering.
  std::vector<t_Eqn> kept;  kept.reserve( _mEqn.size() );
  std::vector<bool> is_evo_row;
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
      if( it != by_eqn.end() ){ e.var = it->second->rest; evo = !it->second->mass.empty(); }
      kept.push_back( std::move( e ) ); is_evo_row.push_back( evo );
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

  {
    size_t order = 0; bool any_diff = false;
    for( auto const& r : _peel.rows ) if( !r.mass.empty() ){ any_diff = true; order = std::max( order, _spatial_order( r.rest ) ); }
    _shape.evolution = !any_diff? t_SemiDiscrete::STEADY: order >= 2? t_SemiDiscrete::PARABOLIC: order == 1? t_SemiDiscrete::HYPERBOLIC: t_SemiDiscrete::ODE;
  }
  _disp( 1 ) << "OCFEMOL::setup [mol] model along " << t << ": " << evolution_type_name( _shape.evolution )
             << " (the classification that follows is of the spatial operator at fixed " << t << ", which is what the plan collocates)\n";
  _disp( 1 ) << "OCFEMOL::setup [mol] peel along " << t << ": " << _peel.rows.size() << " rows (INTERIOR+BOUNDARY), "
             << _peel.initial.size() << " INITIAL, " << _peel.link.size() << " LINK/INTERFACE; "
             << _peel.n_differential() << " of " << _peel.states.size() << " fields differential\n";

  // 3. the shared layout over the spatial directions
  return _discretise_colloc_model();
}

inline bool
OCFEMOL::_on_model_ready
()
{
  if( !_build_colloc_plan() ) return false;
  _issetup = true;

  auto const evo = _shape.evolution;
  _shape = t_SemiDiscrete();
  _shape.evolution     = evo;
  _shape.ran           = true;
  _shape.n_states      = _nCollVar + _nTraceVar;
  _shape.n_tau         = _nTraceVar;
  _shape.n_rows        = _nCollEqn;
  _shape.n_ini         = _iniEqn.size();
  _shape.n_fields      = _peel.states.size();
  _shape.n_diff_fields = _peel.n_differential();
  _shape.t_lo = _evoLo;  _shape.t_hi = _evoHi;

  _disp( 1 ) << "OCFEMOL::setup [mol] semi-discrete: " << _shape.n_states << " unknowns (" << _nCollVar << " nodal + "
             << _nTraceVar << " tau), " << _shape.n_rows << " rows, " << _shape.n_ini << " INITIAL rows aside, t in ["
             << _shape.t_lo << ", " << _shape.t_hi << "]\n";

  size_t const bad = _audit_initial_coverage( _disp( 1 ) );
  if( bad ){
    _disp( 1 ) << "OCFEMOL::setup [mol] DOF audit: " << bad << " field(s) with initial-condition coverage not matching their"
                  " differential status (see above)\n";
    // report-only in this step: the imposition step decides what is refused
  }
  return true;
}

inline size_t
OCFEMOL::_audit_initial_coverage
( std::ostream& os )
const
{
  // Per field: is it differential (a mass entry somewhere), and do INITIAL rows cover it?  On the model's masks, a
  // field's INITIAL row is one whose expression depends on the field and whose spatial masks are ALL (or the full
  // domain).  The node-level version -- which node has which row at t = 0 -- comes with the plan's assignment.
  size_t bad = 0;
  for( size_t i = 0; i < _peel.states.size(); ++i ){
    FFVar const& x = _peel.states[i];
    bool const diff = _peel.differential[i];
    size_t n_ini = 0;
    for( auto const& e : _iniEqn ){
      FFDep dep;  // does the INITIAL row read x?
      auto sg = _dag->subgraph( 1, &e.var );
      bool reads = false;
      for( auto const* op : sg.l_op ){ if( op && op->varout.size() && op->varout[0] && op->varout[0]->id() == x.id() ){ reads = true; break; } }
      if( reads ) ++n_ini;
    }
    if( diff && n_ini == 0 ){
      os << "OCFEMOL::setup [mol] DOF audit: differential field " << x << " has no INITIAL row\n"; ++bad;
    }
    else if( !diff && n_ini ){
      os << "OCFEMOL::setup [mol] DOF audit: algebraic field " << x << " has " << n_ini << " INITIAL row(s): algebraic fields take"
            " their initial values from the collocated rows\n"; ++bad;
    }
  }
  return bad;
}

} // namespace mc

#endif
