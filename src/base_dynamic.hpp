// Copyright (C) Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the EPL-2.0 with GPL-2.0-or-later as a Secondary License; see the LICENSE file.
// SPDX-License-Identifier: EPL-2.0 OR GPL-2.0-or-later

#ifndef CRONOS__BASE_DYNAMIC_HPP
#define CRONOS__BASE_DYNAMIC_HPP

#include <map>
#include <string>
#include <vector>

#include "base_kinsol.hpp"
#include "ffmodel.hpp"

namespace mc
{
//! @brief C++ base class holding the DYNAMIC MODEL DESCRIPTION an integrator extracts from an FFModel.
////////////////////////////////////////////////////////////////////////
//! mc::BASE_DYNAMIC is a base class holding the working copy of a model as an integrator reads it:
//! the states, parameters, constants, quadratures, the (coefficient, remainder) pair of each
//! evolution row, the initial values, the quadrature integrands and the functions.  It carries NO
//! solver and NO extraction policy -- how a row is read into those members is the derived class's
//! business, and the two consumers differ: ODESLV_BASE divides each row by its coefficient to get
//! an explicit right-hand side, while a DAE consumer on IDAS wants the undivided pair as a residual.
//! Hoisted out of ODESLV_BASE when that second consumer arrived, exactly as the note there said.
////////////////////////////////////////////////////////////////////////
class BASE_DYNAMIC
: public virtual FFModel
{
 protected:

  // ---- NAMING CONVENTION, and it matters: _m* below live in the MODEL's dag (FFModel::_dag); _v* and _p* live
  // ---- in THIS class's own _dag, which _SETUP() creates.  An unqualified _dag resolves to the SOLVER's, so any
  // ---- code touching an _m* vector must say FFModel::_dag explicitly.  Mixing them does not fail loudly: the two
  // ---- graphs can carry the same indices, so the wrong one reads whatever sits there (it segfaulted inside
  // ---- set_sparse until that call was qualified).
  // ---- the model description, extracted from FFModel and still in ITS dag ---------------------------------
  // _SETUP() inserts these into this class's own _dag, exactly as it used to insert the old base class's.
  std::vector<double>                  _dT;     //!< stage times, from the evolution domain's elements -- plus
                                                //!< the time points of any function evaluated off that partition
  std::vector<FFVar>                   _mT;     //!< the evolution variable
  std::vector<FFVar>                   _mX;     //!< states
  std::vector<FFVar>                   _mP;     //!< parameters
  std::vector<FFVar>                   _mC;     //!< constants
  std::vector<FFVar>                   _mQ;     //!< quadrature variables
  std::map<std::string,std::vector<FFVar>> _mInpLevels;

  //! @brief Per time-varying input, its node positions WITHIN an element, normalised to [0,1].  Size 1 and
  //! value 0 for a piecewise-constant control, in which case the interpolant is just the nodal value.
  std::map<std::string,std::vector<double>> _mInpNodes;  //!< a time-varying input -> its per-element levels, which
                                                         //!< enter as PARAMETERS (see _extract_from_model)
  //! the coefficient multiplying each derivative (1 if explicit), PER STAGE exactly like _mDE below: one entry
  //! when no input is time-varying, else _ns entries with that stage's input levels composed in.  A coefficient
  //! can depend on a time-varying input (a control-dependent holdup, say), so it needs the same substitution the
  //! right-hand side gets -- it was a flat vector while nothing read it, which would have silently handed a DAE
  //! consumer the wrong stage's levels.
  std::vector<std::vector<FFVar>>      _mLHS;

  //! GENERAL MASS-MATRIX FORM of the evolution rows, for a consumer that wants them UNDIVIDED (a DAE on
  //! IDAS).  Row k of stage s reads   sum over _mMASS[s][k] of (coefficient * d _mX[index] / dt)  +  _mRES[s][k] = 0 .
  //! The per-row list is SPARSE -- only the derivatives actually present -- so a COUPLED mass matrix is
  //! representable here, which _mLHS (one coefficient per state) cannot do: _mLHS is the DIAGONAL projection
  //! ODESLV needs, this is the general form.  Filled only when _wants_mass_rows(), so ODESLV never pays for it.
  std::vector< std::vector< std::vector< std::pair<size_t,FFVar> > > > _mMASS;
  //! The INITIAL rows as RESIDUALS, undivided, in declaration order.  _mIC holds the same rows reduced to
  //! one value per state, which is only faithful when each row is affine in the state it pins; these are what
  //! an initial-condition SOLVE needs.  Filled only when _wants_initial_rows().
  std::vector<FFVar>                   _mICRES;

  //! The model's derivative handle for each state, d _mX[i] / dt, in extraction order.  An INITIAL row may
  //! depend on dx/dt(0) in ANY way -- h( x(0), y(0), dx/dt(0) ) = 0 is the general form -- so a consumer
  //! SUBSTITUTES these for its own dx/dt variables (FFGraph::substitute handles auxiliary nodes, which a
  //! derivative is) rather than splitting the row, which would assume it is AFFINE in each derivative.
  std::vector<FFVar>                   _mDOT;

  //! @brief Plain DAG variables standing for the state time-derivatives, and for the LEFT state at a
  //! stop.  Created ONCE with the rest of the local model, beside _mDOT, because they are part of it:
  //! the IC system is the rows with every FFPartial(.,t) node substituted for one of these, so minting
  //! them per _build_ic_system() call accumulated ~4 variables in the model's dag per restart.  They
  //! are rebuilt whenever the model is, which is the point -- a set held across a re-extraction would
  //! name variables of the PREVIOUS dag.
  std::vector<FFVar>                   _mDotVar;
  std::vector<FFVar>                   _mLftVar;

  //! The CONSISTENT-INITIALISATION system, built once at extraction and owned here because it is the same
  //! for every consumer: the INTERIOR rows at t0 together with the INITIAL rows, with every FFPartial(.,t)
  //! SUBSTITUTED for a plain variable, in its own graph.  Nothing in it is ODE- or DAE-specific -- it never
  //! asks what FORM the rows are in, only that derivative nodes become unknowns -- so neither the mass-row
  //! decomposition nor a solver's local graph is involved.  The derivative unknowns are exactly those that
  //! survive substitution: an algebraic state contributes none because its derivative appears nowhere.
  //! @brief One assembled system: its graph, rows, unknowns and the sparsity KINSOL factorises.  TWO of
  //! these are built, ONCE each, at the end of extraction -- never per call, and never per stage:
  //!   _icStart  the INTERIOR rows at t0 together with the INITIAL rows;
  //!   _icStop   the INTERIOR rows together with CONTINUITY on the differential states, for every stop.
  //! They differ only in which rows arrive, which is §3.13's point that a stage restart IS a consistent
  //! initialisation.  Neither carries a stage index: the time-varying inputs are left as the placeholders
  //! _mUpVar and BOUND per stop, because stages do not carry different equations -- they carry the same
  //! rows at a different instant.  One build, one SFAD pattern, one CSC structure, so KLU keeps its
  //! SYMBOLIC factorisation across every restart.
  struct t_ICSys
  {
    FFGraph*                           dag = nullptr;
    std::vector<FFVar>                 row;      //!< the rows, in dag
    std::vector<FFVar>                 unk;      //!< the states, then the derivatives that appear
    std::vector<size_t>                dotOf;    //!< state index of each derivative unknown
    std::vector<FFVar>                 eval;     //!< [ rows ; the SFAD Jacobian entries ]
    std::vector<FFVar>                 all;      //!< every variable, in dVar order
    std::vector<double>                dEval, dVar, wrk;
    size_t                             nnz = 0;  //!< structural nonzeros of the Jacobian
    std::vector<size_t>                jrow, jcol, at;   //!< at[e]: CSC slot of SFAD entry e
    BASE_KINSOL::Sparsity              sp;
    size_t                             nLft = 0; //!< left-state placeholders in dVar, after the derivatives
    bool                               ok = false;
    ~t_ICSys(){ delete dag; }
  };
  t_ICSys                              _icStart, _icStop;

  //! @brief The system the KINSOL callbacks are currently evaluating.  Set by _ic_solve() before the call;
  //! the callbacks are static trampolines and have nothing else to go on.
  t_ICSys*                             _icCur = nullptr;

  //! @brief The evolution rows BEFORE the per-stage input substitution, kept because the IC systems are
  //! built from them: once subst_inputs() has run, the input variable is gone from the expression and there
  //! is nothing left to parameterise.  _mMASS/_mRES stay as they are, per stage, for the integrator.
  std::vector< std::vector< std::pair<size_t,FFVar> > > _mMASSraw;
  std::vector<FFVar>                   _mRESraw;

  //! @brief The time-varying inputs, and one placeholder variable per input standing for its value at the
  //! instant the system is being solved at -- u(tau^+) at a stop, u(t0) at the start.  A jump row needing
  //! u(tau^-) as well gets a second set when transitions land; nothing reads one yet, so none is minted.
  std::vector<FFVar>                   _mTVInp, _mUpVar;

  //! @brief Element boundaries of the evolution domain, to find the element a stop falls in.
  std::vector<double>                  _mElemBnd;

  //! @brief Per time-varying input, the PARAMETER index of each of its level DOFs, so a stop can read the
  //! level values out of the parameter vector it is given.
  std::vector< std::vector<size_t> >   _mInpParNdx;

  //! Each evolution row with EVERY time derivative set to zero, per stage: _mRES[s][k].
  std::vector< std::vector<FFVar> >    _mRES;
  std::vector<std::vector<FFVar>>      _mIC;    //!< initial values (one entry: no state discontinuities)
  std::vector<std::vector<FFVar>>      _mDE;    //!< right-hand sides (one entry)
  std::vector<std::vector<FFVar>>      _mQUAD;  //!< quadrature integrands (one entry)
  std::vector<std::map<size_t,FFVar>>  _mFCT;   //!< functions, per stage: index -> expression
  std::string _extractError;                    //!< why the last extraction failed

  size_t _ns = 0;        //!< number of stages: _dT.size()-1
  size_t _nx = 0;        //!< number of states
  size_t _nx0 = 0;       //!< number of initial values
  size_t _nc = 0;        //!< number of constants
  size_t _np = 0;        //!< number of parameters
  size_t _nq = 0;        //!< number of quadratures
  size_t _nf = 0;        //!< number of state functions (ROWS of val_function(): a distributed output has one per stage time)
  //! the block [first, count] of each user output (in add_output() order) among the functions / val_function()
  std::vector< std::pair<size_t,size_t> > _blkFct;

  //! @brief POLICY HOOK.  One evolution row has been read as  L * dx_i/dt + R  by _split_in; decide what the
  //! consumer stores for state @p i.  Return false to skip the row.  The two consumers differ here and ONLY
  //! here: ODESLV_BASE refuses a vanishing L (that row is algebraic) and divides through to get an explicit
  //! right-hand side; a DAE consumer on IDAS keeps L == 0 as a constraint and stores the pair UNDIVIDED, so
  //! that the residual is  F_i = lhs[i] * dy_i/dt - de[i]  in both cases.
  virtual bool _read_evolution_row
    ( FFVar const& L, FFVar const& R, size_t i, std::vector<FFVar>& lhs, std::vector<FFVar>& de )
    = 0;

  //! @brief Whether to keep the evolution rows in their GENERAL mass-matrix form, M.dx + R, as well as
  //! divided.  Now TRUE by default: the undivided rows are what the shared consistent-initialisation
  //! system is built from, so every consumer holds them, and the divided form is being phased out.
  virtual bool _wants_mass_rows
    ()
    const
    { return true; }

  //! @brief Whether this consumer is stopped by the model's dynamic form.  The default refuses whatever
  //! dynamic_form() flags, which is what an explicit integrator wants.  A DAE consumer on IDAS accepts the
  //! two forms that exist precisely for it -- a coupled mass matrix (decoupled == false) and algebraic
  //! states -- and keeps the rest of the refusals; it decides that from the STRUCTURED fields, not the message.
  virtual bool _blocked_by_dynamic_form
    ()
    const
    { return !dynamic_form().blocker.empty(); }

  //! @brief Whether the consumer can SOLVE the INITIAL rows rather than evaluate them once.  False refuses
  //! rows implicit in each other at setup, which is right for a consumer that only evaluates _mIC; a
  //! consumer with an initial-condition solve overrides it and takes them.
  virtual bool _solves_initial_rows
    ()
    const
    { return true; }

  //! @brief Whether to keep the INITIAL rows as residuals as well as values.  Now TRUE by default: the
  //! divided form _mIC is being phased out in favour of solving the rows jointly with the dynamics
  //! residual at t0, in BOTH solvers, so the rows themselves are what every consumer should hold.  The
  //! hook stays virtual for a consumer that genuinely wants neither.
  virtual bool _wants_initial_rows
    ()
    const
    { return true; }

  //! @brief The knobs of the IC solve, from the consumer's own INITIAL group.  A hook rather than a
  //! member: the knobs belong to the solver's one `options` object, as everywhere else, and BASE_DYNAMIC
  //! has none.  The default is the plain KINSOL defaults, which is what a consumer with no INITIAL group
  //! would have got anyway.
  virtual BASE_KINSOL::Knobs _ic_knobs
    ()
    const
    { return BASE_KINSOL::Knobs(); }

  //! @brief Whether to take the explicit fast path where every INITIAL row pins one state linearly with
  //! no derivative.  There the divided form is exact and a Newton iteration buys nothing, so this is on
  //! by default; a consumer exposes it so the solve can be forced for comparison.
  virtual bool _ic_bypass_explicit
    ()
    const
    { return true; }

  //! @brief Whether every state must be given its derivative EXPLICITLY by some row.  True for an explicit
  //! integrator; a DAE consumer overrides to false, since an algebraic state has no such row by definition.
  virtual bool _require_explicit_rows
    ()
    const
    { return true; }

  //! @brief Fill the description above from this model.  False, with _extractError set, when the model is not one
  //! an integrator can take (dynamic_form().blocker says why) or a row is not solvable for the derivative it
  //! defines.  Folded in here rather than put in a class of its own: these are THIS solver's working copies, and
  //! a second consumer -- a DAESLV on IDAS, which would want the _mLHS/_mDE pair as a residual -- is the moment
  //! to hoist them, not before.
  bool _extract_from_model
    ();

 public:

  //! @brief Integrator status (was BASE_DE::STATUS -- a solver's verdict, not a model's)
  enum STATUS{
     NORMAL=0,  //!< Normal execution
     FAILURE,   //!< Integration breakdown
     FATAL      //!< Interruption due to errors in third-party libraries
  };

  //! @brief Flat index, in the parameter vector (var_parameter() order -- the vector FFODESLV's operator() takes),
  //! of the START of element block @p ndx_el of the declared input @p var: the same contract as
  //! OCFESLV::pos_input.  A time-invariant input has a single block (@p ndx_el ignored).  A distributed input
  //! has one block per element of its evolution domain, keyed in @p ndx_el by that domain variable (element 0
  //! if absent); within a block its n_node nodal DOFs are contiguous, at pos_input(...) + j.  Throws if @p var
  //! is not an input of the extracted model or the element is out of range.  Valid after setup().
  size_t pos_input
    ( FFVar const& var, std::map<FFVar,size_t,lt_FFVar> const& ndx_el = std::map<FFVar,size_t,lt_FFVar>() )
    const
    {
      auto const ndx = parameter_index( var );
      if( ndx.empty() ) throw std::runtime_error( "ODESLV::pos_input ** " + var.name() + " is not a parameter-input of the extracted model\n" );
      auto const itn = _mInpNodes.find( var.name() );
      if( itn == _mInpNodes.cend() ) return ndx[0];               // time-invariant: one block
      size_t const n_node = itn->second.size()? itn->second.size(): 1;
      size_t e = 0;
      if( !_mT.empty() ){ auto const ite = ndx_el.find( _mT[0] ); if( ite != ndx_el.cend() ) e = ite->second; }
      if( (e+1)*n_node > ndx.size() ) throw std::runtime_error( "ODESLV::pos_input ** element out of range for " + var.name() + "\n" );
      return ndx[ e*n_node ];
    }

  //! @brief The nodal DOF count of one element block of input @p var (1 for a time-invariant input).
  size_t size_input
    ( FFVar const& var )
    const
    {
      auto const itn = _mInpNodes.find( var.name() );
      return ( itn == _mInpNodes.cend() || itn->second.empty() )? 1: itn->second.size();
    }

  //! @brief Write the DOF values of declared input @p w (control_dofs order) into the flat parameter vector @p P
  //! (size np(), var_parameter() order) -- the ODESLV counterpart of OCFESLV::set_input_values.  Returns false if
  //! @p w is not a parameter-input, the value count is wrong, or @p P is too short.  Valid after setup().
  bool set_parameter_values
    ( FFVar const& w, std::vector<double> const& vals, std::vector<double>& P )
    const
    {
      auto const ndx = parameter_index( w );
      if( ndx.empty() || vals.size() != ndx.size() || P.size() < _np ) return false;
      for( size_t k=0; k<ndx.size(); ++k ) P[ ndx[k] ] = vals[k];
      return true;
    }

  //! @brief Read the DOF values of declared input @p w (control_dofs order) from the flat parameter vector @p P.
  //! Empty if @p w is not a parameter-input or @p P is too short.
  std::vector<double> get_parameter_values
    ( FFVar const& w, std::vector<double> const& P )
    const
    {
      auto const ndx = parameter_index( w );
      std::vector<double> v;
      if( ndx.empty() || P.size() < _np ) return v;
      for( size_t k : ndx ) v.push_back( P[k] );
      return v;
    }

  //! @brief Build the flat parameter vector @p P and constant vector @p c from values given BY NAME: every declared
  //! input (except a fixed one) and every constant, exactly once, each a value list in control_dofs() order or a
  //! generator.  Returns false with the reason in @p err on an unknown or duplicate entry, a wrong count, or a
  //! missing input or constant.
  bool assemble_values
    ( std::vector<FFModel::InputVal> const& vIn, std::vector<double>& P, std::vector<double>& c, std::string& err )
    const
    {
      err.clear();  P.assign( _np, 0. );
      auto const& vC = var_constant();  c.assign( vC.size(), 0. );
      std::vector<FFVar> seen;
      for( auto const& a : vIn ){
        for( auto const& u : seen ) if( u.id().second == a.input.id().second ){ err = a.input.name() + " is given more than once"; return false; }
        seen.push_back( a.input );
        int ic = -1; for( size_t k=0; k<vC.size(); ++k ) if( vC[k].id().second == a.input.id().second ){ ic = (int)k; break; }
        if( ic >= 0 ){                                    // a constant: one value
          double const v = a.gen? a.gen( DofIndex() ): ( a.vals.size() == 1? a.vals[0]: std::nan("") );
          if( !a.gen && a.vals.size() != 1 ){ err = "constant " + a.input.name() + " takes exactly one value"; return false; }
          c[ic] = v;  continue;
        }
        if( parameter_index( a.input ).empty() ){
          err = a.input.name() + " is neither a parameter-input nor a constant of the extracted model"
              + std::string( fixed_input( a.input )? " (it is fixed)": "" );
          return false;
        }
        std::vector<double> v;
        if( a.gen ){ for( auto const& d : control_dofs( a.input ) ) v.push_back( a.gen( d ) ); }
        else v = a.vals;
        if( !set_parameter_values( a.input, v, P ) ){
          err = a.input.name() + " has " + std::to_string( parameter_index( a.input ).size() ) + " DOFs, but "
              + std::to_string( v.size() ) + " values were given";
          return false;
        }
      }
      std::string missing;
      for( auto const& [w,dom] : var_declared_input() ){
        if( fixed_input( w ) ) continue;
        bool in = false; for( auto const& u : seen ) if( u.id().second == w.id().second ){ in = true; break; }
        if( !in ) missing += " " + w.name();
      }
      for( auto const& k : vC ){
        bool in = false; for( auto const& u : seen ) if( u.id().second == k.id().second ){ in = true; break; }
        if( !in ) missing += " " + k.name();
      }
      if( !missing.empty() ){ err = "every declared input and constant must be given; missing:" + missing; return false; }
      return true;
    }

  //! @brief The declared input a parameter @p p belongs to, and p's DOF index within it (element-major for
  //! a distributed input, 0 for a time-invariant one): the inverse of parameter_index().  Returns false if
  //! @p p is not a parameter of the extracted model.  Valid after setup().
  bool input_of_parameter
    ( FFVar const& p, FFVar& w, size_t& dof ) const
    {
      for( auto const& [win, dom] : var_declared_input() ){
        auto const ndx = parameter_index( win );
        for( size_t j = 0; j < ndx.size(); ++j )
          if( _mP[ ndx[j] ].id().second == p.id().second ){ w = win; dof = j; return true; }
      }
      return false;
    }

  //! @brief The parameter indices (in var_parameter() order) of the DOFs of the declared input @p w: one
  //! for a time-invariant input, its minted levels for a distributed one, empty if @p w is not an input of
  //! the extracted model.  Valid after setup().  This is the layout query fdiff_seed_all() relies on --
  //! time-invariant parameters come FIRST in the extracted order, so no positional assumption holds.
  std::vector<size_t> parameter_index
    ( FFVar const& w ) const override
    {
      std::vector<FFVar> dofs;
      auto const itl = _mInpLevels.find( w.name() );
      if( itl != _mInpLevels.cend() ) dofs = itl->second; else dofs.assign( 1, w );
      std::vector<size_t> ndx;
      for( auto const& d : dofs )
        for( size_t i=0; i<_np; ++i ) if( _mP[i].id().second == d.id().second ){ ndx.push_back( i ); break; }
      return ndx;
    }
 public:

  /** @{
   */
  //! @brief Sizes of the extracted description, and why the last extraction failed.  These report members
  //! _extract_from_model() fills, so they belong beside them rather than in one consumer: an ODE and a DAE
  //! solver answer them identically.
  //! @brief Number of state functions (was BASE_DE::nf)
  size_t nf
    ()
    const
    { return _nf; }

  //! @brief Number of states (was BASE_DE::nx)
  size_t nx
    ()
    const
    { return _nx; }

  //! @brief Rows of the consistent-initialisation system (INTERIOR at t0 + INITIAL), in _icDag.
  std::vector<FFVar> const& ic_rows
    ()
    const
    { return _icStart.row; }

  //! @brief Its unknowns: the states, then the state derivatives that appear in some row.
  std::vector<FFVar> const& ic_unknowns
    ()
    const
    { return _icStart.unk; }

  //! @brief Whether every INITIAL row pins one state EXPLICITLY (x_i - c, no other state, no derivative).
  //! Then the divided form is exact and a Newton solve buys nothing -- the common case, worth bypassing.
  bool explicit_initial_conditions
    ()
    const;

  //! @brief Reinitialise consistently at a stop at time @p tau: the evolution rows in force after it,
  //! together with CONTINUITY on the differential states, the left state @p xleft as data.  The same solve
  //! as at t0 -- a stage restart is a consistent initialisation, not a separate mechanism (plan §3.13).
  //! The system is NOT rebuilt here: it was assembled once at extraction, and all that varies between
  //! stops is @p tau and the input values bound from it, so KLU keeps its symbolic factorisation across
  //! every restart.  Takes the stop TIME rather than a stage index, because a transition's tau need not
  //! fall on a stage boundary.
  BASE_KINSOL::Report restart_ic
    ( double const tau, std::vector<double> const& xleft, std::vector<double>& x, std::vector<double>& dx,
      BASE_KINSOL::Knobs const& k, std::vector<double> const& par = std::vector<double>() )
    { if( !_icStop.ok ){ BASE_KINSOL::Report R;
        R.message = _extractError.empty()? "no restart system was built": _extractError; return R; }
      _icLeft = xleft;
      BASE_KINSOL::Report const R = _ic_solve( _icStop, tau, x, dx, k, par );
      _icLeft.clear();
      return R; }

  //! @brief Solve the consistent-initialisation system, leaving the result in @p x and @p dx.  The ONE
  //! entry point both solvers take at t0: ODESLV through _IC_D_STA, DAESLV through its setup, so the
  //! initial state and its derivative come from the same system on both sides rather than from two.
  BASE_KINSOL::Report solve_ic
    ( std::vector<double>& x, std::vector<double>& dx, BASE_KINSOL::Knobs const& k,
      std::vector<double> const& par = std::vector<double>() );

  //! @brief The model's d _mX[i] / dt handles, to substitute out of the INITIAL rows.
  std::vector<FFVar> const& state_derivatives
    ()
    const
    { return _mDOT; }

  //! @brief The INITIAL rows as residuals, undivided; empty unless the consumer asked for them.
  //! THE CONSISTENT-IC SYSTEM.  In general an INITIAL row may depend on a state's TIME DERIVATIVE at the
  //! initial point, so the solve cannot be posed over the states alone: its unknowns are (x, dx/dt) and its
  //! equations are these rows TOGETHER WITH the INTERIOR rows evaluated at t0, which are what relate dx/dt
  //! to x.  The dx/dt unknowns are only the DIFFERENTIAL ones -- an algebraic state's derivative appears in
  //! no row -- so the system is square when the number of INITIAL rows equals the number of differential
  //! states.  For the index-3 pendulum that is 4 rows against 4 differential states among 5, giving 9
  //! equations in 9 unknowns.  This is the same system IDACalcIC(IDA_YA_YDP_INIT) solves, with the id
  //! vector saying which half of x is given.
  //! CAUTION: these live in the WORKING model's dag, not necessarily FFModel::_dag -- the reduction mints
  //! its auxiliaries there -- so evaluate each through its OWN row.dag().  The two graphs can carry the
  //! same indices, which is the trap the naming-convention note above records.
  std::vector<FFVar> const& initial_rows
    ()
    const
    { return _mICRES; }

  //! @brief Number of constants
  size_t nc
    ()
    const
    { return _nc; }

  //! @brief Number of parameters (was BASE_DE::np)
  size_t np
    ()
    const
    { return _np; }

  //! @brief Number of quadratures (was BASE_DE::nq)
  size_t nq
    ()
    const
    { return _nq; }

  //! @brief Number of stages (was BASE_DE::ns)
  size_t ns
    ()
    const
    { return _ns; }

  //! @brief Why the last setup() failed to extract the model; empty after a successful one.  Usually the
  //! dynamic_form() blocker: the model is distributed, or a boundary-value problem in the evolution direction, or
  //! its derivatives are coupled through a mass matrix.
  std::string const& extract_error
    ()
    const
    { return _extractError; }
  /** @} */

 protected:

  //! @brief Residual and dense Jacobian of the IC system, for BASE_KINSOL.
  static void _ic_scatter( BASE_DYNAMIC* M, int const n, double const* z );

  static int _ic_res
    ( int const n, double const* z, double* f, void* user );
  static int _ic_jac
    ( int const n, double const* z, double* J, void* user );
  static int _ic_spjac
    ( int const n, double const* z, double* val, void* user );

  //! @brief Every variable the IC rows are evaluated against, in _dIcVar order.
  std::vector<FFVar> const& _icUnkAll
    ()
    const
    { return _icCur? _icCur->all: _icStart.all; }

  //! @brief The starting point: the declared references (virtual, so OCFESLV's sampling is honoured).
  std::vector<double> _ic_seed
    ()
    const
    { std::vector<double> sref, iref, cref, dref;
      const_cast<BASE_DYNAMIC*>(this)->_classification_reference( sref, iref, cref, dref );
      sref.resize( _nx, 0. );
      return sref; }

  //! @brief The element-local Lagrange interpolant of @p nodal at @p sc, with @p sn the node positions
  //! normalised to [0,1].  Templated so there is ONE basis in the code: subst_inputs() instantiates it on
  //! FFVar to bake the interpolant into the rows, and a stop instantiates it on double to bind u(tau).
  //! A piecewise-constant input has a single node and the interpolant is its level, as it must be.
  template <typename T>
  static T _interp_level
    ( std::vector<double> const& sn, std::vector<T> const& nodal, T const& sc )
    {
      if( sn.size() < 2 || nodal.size() < sn.size() ) return nodal.at( 0 );
      T out( 0. );
      for( size_t j = 0; j < sn.size(); ++j ){
        T Lj( 1. );
        for( size_t k = 0; k < sn.size(); ++k )
          if( k != j ) Lj = Lj * ( ( sc - sn[k] ) / ( sn[j] - sn[k] ) );
        out = out + Lj * nodal[j];
      }
      return out;
    }

  //! @brief The element @p tau falls in, on the side @p side, from the stored boundaries.  The extraction
  //! lambda of the same name, kept as a member because a stop needs it long after extraction has returned.
  size_t _elem_of
    ( double const tau, int const side )
    const
    {
      if( _mElemBnd.size() < 2 ) return 0;
      size_t const ne = _mElemBnd.size() - 1;
      size_t e = (size_t)( std::upper_bound( _mElemBnd.begin(), _mElemBnd.end(), tau ) - _mElemBnd.begin() );
      e = e? e-1: 0;  if( e >= ne ) e = ne - 1;
      if( side == FFDom::MINUS ){ if( e > 0 && std::fabs( tau - _mElemBnd[e] ) <= options.TTOL ) --e; }
      else if( e+1 < ne && std::fabs( tau - _mElemBnd[e+1] ) <= options.TTOL ) ++e;
      return e;
    }

  //! @brief Is state @p i differential, i.e. does its derivative appear in some evolution row?
  bool _differential_state
    ( size_t const i )
    const
    { for( auto const& st : _mMASS ) for( auto const& row : st )
        for( auto const& [ ix, co ] : row ) if( ix == i ) return true;
      return false; }

  //! @brief Build BOTH consistent-initialisation systems, once, at the end of extraction.
  bool _build_ic_system
    ();

  //! @brief Assemble one of them: @p at_stop picks continuity rows over the INITIAL rows.
  bool _assemble_ic
    ( t_ICSys& S, bool const at_stop );

  //! @brief Solve @p S at @p tau, with the time-varying inputs bound to their values there.
  BASE_KINSOL::Report _ic_solve
    ( t_ICSys& S, double const tau, std::vector<double>& x, std::vector<double>& dx,
      BASE_KINSOL::Knobs const& k, std::vector<double> const& par );

  //! @brief Constant values for the IC system, in _mC order.  Filled by the consumer where it already
  //! receives them -- the constants are bound per solve(), so they cannot be baked in at build time.
  std::vector<double>                  _dIcCst;

  //! @brief Left-side state values at a stop, the data of the continuity rows.  Sized _nx; only the
  //! DIFFERENTIAL entries are read, since an algebraic state is re-determined by its own rows.
  std::vector<double>                  _icLeft;

  //! @brief Read @p res as L * @p tgt + R, both as EXPRESSIONS, so a state-dependent coefficient is carried
  //! rather than refused.  L comes from differentiating a fresh leaf put in the target's place -- taking
  //! res|tgt=1 - res|tgt=0 instead mixes two substitutions into a graph that crashes on evaluation.
  bool _split_in
    ( FFVar const& res, FFVar const& tgt, FFVar& L, FFVar& R );

  //! @brief @p v, guaranteed to live in the model's dag: a constant extracted from a row carries none of its own.
  FFVar _in_model_dag
    ( FFVar const& v )
    const
    { return v.dag()? v: FFVar( FFModel::_dag, v.num().val() ); }

};

inline
bool
BASE_DYNAMIC::_extract_from_model
()
{
  _dT.clear(); _mT.clear(); _mX.clear(); _mP.clear(); _mC.clear(); _mQ.clear(); _mLHS.clear();
  _mIC.clear(); _mDE.clear(); _mQUAD.clear(); _mFCT.clear();
  _ns = _nx = _nx0 = _nc = _np = _nq = _nf = 0;
  _extractError.clear();

  if( !is_setup() ){ _extractError = "the model has not been set up"; return false; }
  // Every equation must be one ODESLV USES -- INTERIOR (it gives a derivative) or INITIAL (an initial value).  Any
  // other is REFUSED rather than silently ignored (2026-10-01: a BOUNDARY row on t = 0 used to be skipped, leaving the
  // state at 0 with solve() returning NORMAL); checked FIRST, so that the reason given is this one.  ODESLV does not
  // evaluate DIAGNOSTIC rows either.
  {
    size_t k = 0;
    for( auto const& eqn : var_equation() ){
      EqnRole const r = eqn.opt? eqn.opt->role: EqnRole::AUTO;
      if( r != EqnRole::INTERIOR && r != EqnRole::INITIAL ){
        static char const* const NAME[] = { "AUTO", "INTERIOR", "INITIAL", "BOUNDARY", "INTERFACE", "LINK", "SURFACE", "DIAGNOSTIC" };
        int const ir = (int)r;
        _extractError = "equation #" + std::to_string( k ) + " has role "
                      + ( ir >= 0 && ir < 8? NAME[ir]: "?" ) + ", which ODESLV does not use: an initial condition"
                        " has role INITIAL, the dynamics INTERIOR";
        return false;
      }
      ++k;
    }
  }
  if( _blocked_by_dynamic_form() ){ _extractError = dynamic_form().blocker; return false; }

  FFVar const& tvar = _evolution_dom_var;
  auto const itd = var_domain().find( tvar );
  if( itd == var_domain().cend() ){ _extractError = "the evolution direction has no domain"; return false; }
  FFDom const& tdom = itd->second;
  _mT.push_back( tvar );

  // The stage partition: the evolution domain's elements, PLUS the time point of every function evaluated off
  // that partition -- the solver already stops at stage boundaries and evaluates what is registered there, so a
  // function at an arbitrary tau costs nothing more than an extra boundary.
  // A point within options.TTOL of one already held is the SAME point and is dropped: the grid keeps its boundary
  // and the function is evaluated there.  A tolerant std::set comparator would be shorter and wrong -- `a < b-TTOL`
  // makes equivalence non-transitive once points cluster, which std::set does not permit -- so this scans.
  std::vector<double> stages;
  double snapped_to = 0.;
  auto add_time = [&]( double const v ) -> bool {
    for( double const s : stages )
      if( std::fabs( s - v ) <= options.TTOL ){ snapped_to = s; return false; }
    stages.push_back( v );
    return true;
  };
  // the domain's OWN element boundaries (2026-09-28: these used to be recomputed as uniform, which put the stages of
  // an FFDom( elem_bnd, ... ) domain at the wrong times)
  for( size_t k = 0; k <= tdom.n_elem; ++k )
    add_time( tdom.elem_bnd.at( k ) );
  for( auto const& C : var_deferred() ){
    if( C.accumulate || C.tau <= tdom.lo_dom || C.tau >= tdom.up_dom ) continue;
    if( !add_time( C.tau ) && options.DISPLAY_LEVEL >= 1 )
      std::cerr << "ODESLV_BASE:: function time point " << std::setprecision(16) << C.tau
                << " is within TTOL of the stage boundary " << snapped_to << std::setprecision(6)
                << "; it is evaluated there rather than splitting the stage" << std::endl;
  }
  for( auto const& tr : var_transition() )                   // a transition starts a stage (add_transition)
    if( tr.tau > tdom.lo_dom && tr.tau < tdom.up_dom ) add_time( tr.tau );
  std::sort( stages.begin(), stages.end() );
  _dT = stages;
  _ns = _dT.size() - 1;

  for( auto const& [var,dom] : var_state() ) _mX.push_back( var );
  _nx = _nx0 = _mX.size();

  // The parameters are the inputs the MODELLER declared (rev346's var_declared_input), so a deferred-value input
  // -- which setup() mints per captured output -- is never taken for a parameter: telling the two kinds apart is
  // the model's business now, not every consumer's.  The rows reference WORKING nodes, so each declared input is
  // matched to its counterpart there; for every model measured they are the SAME node, and a name match is the
  // fallback for a model that localises its inputs.
  std::vector<FFVar> tv_declared, fixed_scalar;
  for( auto const& [var,dom] : var_declared_input() ){
    FFVar const* w = nullptr;
    for( auto const& [wv,wd] : var_input() ) if( wv.id()   == var.id()   ){ w = &wv; break; }
    if( !w )
      for( auto const& [wv,wd] : var_input() ) if( wv.name() == var.name() ){ w = &wv; break; }
    if( !w ){
      _extractError = "the declared input " + var.name() + " has no counterpart in the working model";
      return false;
    }
    if( dom.empty() ){
      if( !fixed_input( *w ) ) _mP.push_back( *w ); // time-invariant: a parameter -- unless FIXED, then a
      else fixed_scalar.push_back( *w );            // constant substituted into the right-hand sides below
    }
    else              tv_declared.push_back( *w );  // time-varying: a profile, handled below
  }
  _np = _mP.size();

  _mC = var_constant();
  _nc = _mC.size();

  // A TIME-VARYING INPUT is not a parameter: an integrator needs a value, not a profile.  A PIECEWISE-CONSTANT one
  // (declared with one value per element) becomes ONE PARAMETER PER ELEMENT, and the right-hand sides become
  // stage-wise with that element's level substituted -- exactly the description the legacy wrote by hand (test2:
  // NP = 2 + NS, with RHS[k] using U[k]).  Any other profile needs its interpolant evaluated at the evolution
  // variable inside the right-hand side, which is not implemented, so it is REFUSED rather than left dangling.
  // Fixed time-invariant inputs are constants of the extracted model: substituted into EVERY expression that
  // is read from the declaration -- right-hand sides, initial conditions, functions and quadratures alike.
  std::vector<FFVar> fx_tgt, fx_rep;
  for( auto const& w : fixed_scalar ){ fx_tgt.push_back( w ); fx_rep.push_back( FFVar( (*fixed_input( w ))[0] ) ); }
  auto subst_fixed = [&]( std::vector<FFVar> const& v ){ return fx_tgt.empty()? v: FFModel::_dag->compose( v, fx_tgt, fx_rep ); };
  auto subst1 = [&]( FFVar const& f ){ if( fx_tgt.empty() ) return f; std::vector<FFVar> v{ f }; return FFModel::_dag->compose( v, fx_tgt, fx_rep )[0]; };
  _mInpLevels.clear();  _mInpNodes.clear();
  std::vector<FFVar> tv_inputs;
  for( auto const& var : tv_declared ){
    size_t n_node = 0;
    auto const iti = _usr._mInpDiscUsr.find( var );
    if( iti != _usr._mInpDiscUsr.cend() ){
      auto const itd = iti->second.find( tvar );
      if( itd != iti->second.cend() ) n_node = itd->second.n_node;
    }
    if( !n_node ){
      _extractError = "the input " + var.name() + " varies over " + tvar.name()
                    + " but declares no discretisation on it; give it one value per element to hold it"
                    + " piecewise constant, or n_node values for a polynomial profile";
      return false;
    }
    // The node FAMILY is the input's own where it declares one, the domain's otherwise -- the same
    // InpColloc{type,n_node} that FFModel::control_ndof() reads, so the DOF count and the interpolant
    // cannot disagree about the function space.
    int inp_type = (int)tdom.type;
    if( iti != _usr._mInpDiscUsr.cend() ){
      auto const itd2 = iti->second.find( tvar );
      if( itd2 != iti->second.cend() ) inp_type = (int)itd2->second.type;
    }
    std::vector<double> snode( 1, 0. );
    if( n_node > 1 ){
      BASE_OC quad;                                            // nodes only: the evaluation is symbolic below
      if( !quad.set_lgnodes( inp_type, n_node, -1., 1. ) ){
        _extractError = "cannot generate " + std::to_string( n_node ) + " collocation nodes of type "
                      + std::to_string( inp_type ) + " for the input " + var.name();
        return false;
      }
      snode = quad.lgnodes( 0., 1. );                          // element-local, in [0,1]
      if( snode.size() != n_node ){
        _extractError = "the node family for the input " + var.name() + " returned "
                      + std::to_string( snode.size() ) + " nodes, not " + std::to_string( n_node );
        return false;
      }
    }
    tv_inputs.push_back( var );
    _mInpNodes[ var.name() ] = snode;                          // documented as stored; it was not
    std::vector<FFVar> levels;
    for( size_t e = 0; e < tdom.n_elem; ++e )                  // element-major, as control_ndof() counts
      for( size_t j = 0; j < n_node; ++j ){                    // named, so the parameter vector is readable
        std::vector<double> const* fx = fixed_input( var );
        if( fx ){                                                 // FIXED: a constant level, not a parameter
          levels.push_back( FFVar( (*fx)[ e*n_node + j ] ) );
          continue;
        }
        FFVar lvl = FFModel::_dag->add_var( var.name() + "[" + std::to_string( e )
                  + ( n_node > 1? "," + std::to_string( j ): std::string() ) + "]" );
        levels.push_back( lvl );
        _mP.push_back( lvl );
      }
    _mInpLevels[ var.name() ] = levels;
    _mInpNodes [ var.name() ] = snode;
  }
  _np = _mP.size();

  FFPartial OpP;                                     // each state's row, read as L * dx/dt + R = 0
  std::vector<FFVar> rhs( _nx, FFVar( 0. ) );
  std::vector<FFVar> lhs( _nx, FFVar( 1. ) );   // kept undivided; _mLHS is filled per stage below
  std::vector<bool> got( _nx, false );
  for( auto const& eqn : var_equation() ){
    if( !eqn.opt || eqn.opt->role != EqnRole::INTERIOR ) continue;
    for( size_t i = 0; i < _nx; ++i ){
      if( got[i] ) continue;
      FFVar L, R;
      if( !_split_in( eqn.var, OpP( _mX[i], tvar ), L, R ) ) continue;
      if( !_read_evolution_row( L, R, i, lhs, rhs ) ) continue;
      got[i]   = true;
    }
  }
  // General mass-matrix form: peel each derivative off the row in turn -- _split_in returns the row with that
  // derivative SET TO ZERO, so feeding R back in exposes the next one, and what is left when every state has
  // been tried is the remainder.  Row-centric, so nothing here assumes a row belongs to one state.
  std::vector< std::vector< std::pair<size_t,FFVar> > > mass_raw;
  std::vector<FFVar> res_raw;
  if( _wants_mass_rows() ){
    for( auto const& eqn : var_equation() ){
      if( !eqn.opt || eqn.opt->role != EqnRole::INTERIOR ) continue;
      std::vector< std::pair<size_t,FFVar> > row;
      FFVar rest = eqn.var;
      for( size_t i = 0; i < _nx; ++i ){
        FFVar L, R;
        if( !_split_in( rest, OpP( _mX[i], tvar ), L, R ) ) continue;
        if( L.cst() && L.num().val() == 0. ) continue;
        row.emplace_back( i, L );
        rest = R;
      }
      mass_raw.push_back( row );
      res_raw.push_back( _in_model_dag( rest ) );
    }
  }
  if( _require_explicit_rows() )
    for( size_t i = 0; i < _nx; ++i )
      if( !got[i] ){
        _extractError = "no row gives d/d" + tvar.name() + " of " + _mX[i].name()
                      + " explicitly; this model is a DAE for the consumer and belongs to IDAS";
        return false;
      }
  if( !fixed_scalar.empty() ) rhs = subst_fixed( rhs );        // fixed time-invariant inputs: substitute once
  // Time-varying inputs, element by element (2026-09-28: shared by the dynamics, the OUTPUT evaluations and the
  // INTEGRANDS -- outputs and integrands used to keep the raw input variable, which failed at the first solve).
  auto elem_of_stage = [&]( size_t const k ){                // the element holding stage k
    double const mid = 0.5 * ( _dT[k] + _dT[k+1] );
    size_t e = (size_t)( std::upper_bound( tdom.elem_bnd.begin(), tdom.elem_bnd.end(), mid ) - tdom.elem_bnd.begin() );
    e = e? e-1: 0;  if( e >= tdom.n_elem ) e = tdom.n_elem - 1;
    return e; };
  auto elem_of_point = [&]( double const tau, int const side ){  // the element a point evaluation reads (FFDom::SIDE)
    size_t e = (size_t)( std::upper_bound( tdom.elem_bnd.begin(), tdom.elem_bnd.end(), tau ) - tdom.elem_bnd.begin() );
    e = e? e-1: 0;  if( e >= tdom.n_elem ) e = tdom.n_elem - 1;
    if( side == FFDom::MINUS ){ if( e > 0 && std::fabs( tau - tdom.elem_bnd[e] ) <= options.TTOL ) --e; }   // tau^-: ending there
    else if( e+1 < tdom.n_elem && std::fabs( tau - tdom.elem_bnd[e+1] ) <= options.TTOL ) ++e;              // tau^+: starting there
    return e; };
  auto subst_inputs = [&]( std::vector<FFVar> const& ex, size_t const elem ){
    if( tv_inputs.empty() ) return ex;
    double const te = tdom.elem_bnd.at( elem );
    double const he = tdom.elem_bnd.at( elem+1 ) - te;
    std::vector<FFVar> tgt, rep;
    for( auto const& u : tv_inputs ){
      tgt.push_back( u );
      auto const& lv = _mInpLevels[u.name()];
      auto const& sn = _mInpNodes [u.name()];
      size_t const nn = sn.size();
      if( nn < 2 ){ rep.push_back( lv[elem] ); continue; }   // piecewise constant: the level itself
      FFVar const sc = ( tvar - te ) / he;                   // element-local coordinate in [0,1]
      std::vector<FFVar> nodal;
      for( size_t j = 0; j < nn; ++j ) nodal.push_back( lv[elem*nn + j] );
      rep.push_back( _interp_level( sn, nodal, sc ) );
    }
    return FFModel::_dag->compose( ex, tgt, rep ); };
  auto subst_mass = [&]( size_t const elem ){            // the same input substitution, on the sparse rows
    std::vector<FFVar> flat;
    for( auto const& row : mass_raw ) for( auto const& [ix,co] : row ) flat.push_back( co );
    for( auto const& r : res_raw ) flat.push_back( r );
    std::vector<FFVar> const done = subst_inputs( flat, elem );
    std::vector< std::vector< std::pair<size_t,FFVar> > > M( mass_raw.size() );
    size_t n = 0;
    for( size_t k = 0; k < mass_raw.size(); ++k )
      for( auto const& [ix,co] : mass_raw[k] ) M[k].emplace_back( ix, done[n++] );
    std::vector<FFVar> R;
    for( size_t k = 0; k < res_raw.size(); ++k ) R.push_back( done[n++] );
    return std::make_pair( M, R ); };
  _mMASSraw = mass_raw;  _mRESraw = res_raw;   // the IC systems are built from these, not from _mRES
  _mTVInp   = tv_inputs;
  _mElemBnd.assign( tdom.elem_bnd.begin(), tdom.elem_bnd.end() );
  if( tv_inputs.empty() ){
    _mDE .assign( 1, rhs );
    _mLHS.assign( 1, lhs );
    _mMASS.assign( 1, mass_raw );
    _mRES .assign( 1, res_raw );
  }
  else{                                                    // one entry per stage, with that element's levels in
    _mDE .assign( _ns, std::vector<FFVar>( _nx ) );
    _mLHS.assign( _ns, std::vector<FFVar>( _nx ) );
    _mMASS.assign( _ns, mass_raw );
    _mRES .assign( _ns, res_raw );
    for( size_t k = 0; k < _ns; ++k ){
      std::vector<FFVar> const doneR = subst_inputs( rhs, elem_of_stage( k ) );
      std::vector<FFVar> const doneL = subst_inputs( lhs, elem_of_stage( k ) );
      for( size_t i = 0; i < _nx; ++i ){ _mDE[k][i] = doneR[i]; _mLHS[k][i] = doneL[i]; }
      if( _wants_mass_rows() ){ auto M = subst_mass( elem_of_stage( k ) ); _mMASS[k] = M.first; _mRES[k] = M.second; }
    }
  }



  std::vector<FFVar> ic( _nx, FFVar( 0. ) );
  _mICRES.clear();  _mDOT.clear();
  _mDotVar.clear();  _mLftVar.clear();
  if( _wants_initial_rows() ){
    _mDOT.resize( _nx );
    for( size_t i = 0; i < _nx; ++i ) _mDOT[i] = OpP( _mX[i], tvar );
  }
  for( auto const& eqn : var_equation() ){
    if( !eqn.opt || eqn.opt->role != EqnRole::INITIAL ) continue;
    // Keep the row UNDIVIDED, as _mRES keeps the evolution rows.  _split_in factorises a row as
    // L * target + R with L = d(res)/d(target), which is the row only when it is AFFINE in that target;
    // a materialised hidden constraint is not -- x^2 + y^2 - 1 split on y gives ic[y] = (1+y^2-x^2)/(2y),
    // self-referential.  The divided form LOSES the row, so an IC solve has to work from these.
    if( _wants_initial_rows() ) _mICRES.push_back( _in_model_dag( subst1( eqn.var ) ) );
    for( size_t i = 0; i < _nx; ++i ){
      FFVar L, R;
      if( !_split_in( subst1( eqn.var ), _mX[i], L, R ) ) continue;
      if( L.cst() && L.num().val() == 0. ) continue;
      ic[i] = ( L.cst() && L.num().val() == 1. )? _in_model_dag( -R ): _in_model_dag( -R / L );
    }
  }
  // An INITIAL row is IMPLICIT when its extracted value still depends on another STATE: _split_in leaves
  // the rest of the row in R, so two coupled initial conditions become ic[i] = f(x_j).  ONE evaluation
  // cannot resolve them, and the two consumers then fail differently and late -- ODESLV with a
  // missing-variable error out of the subgraph, a DAE consumer by reading whatever its working vector
  // happened to hold.  Neither names the model, so say it here, while the rows are in hand.
  {
    std::vector<FFVar> dIC;
    try{ dIC = FFModel::_dag->FAD( ic, std::vector<FFVar>( _mX.begin(), _mX.end() ) ); }
    catch( ... ){ dIC.clear(); }
    if( dIC.size() == ic.size() * _nx )
      for( size_t i = 0; i < ic.size(); ++i )
        for( size_t j = 0; j < _nx; ++j ){
          FFVar const& d = dIC[ i * _nx + j ];
          if( d.cst() && d.num().val() == 0. ) continue;
          if( _solves_initial_rows() ) break;          // the consumer resolves them; nothing to refuse
          _extractError = "the initial condition for " + _mX[i].name() + " still depends on "
                        + _mX[j].name() + ": the INITIAL rows are implicit in each other, and resolving"
                          " them needs an initial-condition SOLVE, which this solver does not do";
          return false;
        }
  }
  _mIC.assign( 1, ic );
  // TRANSITIONS (add_transition, 2026-09-29): one initial-value entry PER STAGE -- stage 0 the INITIAL equations, the
  // stage starting at a transition's tau its EXPLICIT map x_i(tau^+) = left_i(x(tau^-), p) (right_i must be the state
  // x_i itself), identity elsewhere.  ODESLV evaluates entry k with the state reached at the end of stage k-1, i.e.
  // x(tau^-), and reinitialises the integrator with the result.  Implicit maps are refused (DAESLV's business).
  std::vector< std::pair< double, std::map<size_t,FFVar> > > jumps;   // (tau, state index -> x_i(tau^+)) for PLUS evaluations
  if( !var_transition().empty() ){
    _mIC.assign( _ns, std::vector<FFVar>( _nx ) );
    _mIC[0] = ic;
    for( size_t k = 1; k < _ns; ++k ) for( size_t i = 0; i < _nx; ++i ) _mIC[k][i] = _mX[i];
    for( auto const& tr : var_transition() ){
      size_t ks = 0;
      for( size_t k = 1; k < _ns; ++k ) if( std::fabs( _dT[k] - tr.tau ) <= options.TTOL ){ ks = k; break; }
      if( !ks ){ std::ostringstream os; os << "transition at t = " << tr.tau << " is not a stage boundary"; _extractError = os.str(); return false; }
      size_t const el = elem_of_point( tr.tau, FFDom::MINUS );
      std::map<size_t,FFVar> mp;
      for( size_t c = 0; c < tr.right.size(); ++c ){
        FFVar const rw = _in_model_dag( tr.right[c] );
        size_t i = _nx;
        for( size_t j = 0; j < _nx; ++j ) if( _mX[j].id() == rw.id() ){ i = j; break; }
        if( i == _nx ){
          std::ostringstream os; os << "transition at t = " << tr.tau << ", component " << c << ": ODESLV accepts EXPLICIT maps"
            " only -- the right-hand side (evaluated at tau^+) must be a single state, x_i(tau^+) = g(x(tau^-), p)";
          _extractError = os.str(); return false;
        }
        for( auto const& [t0, m0] : jumps ) if( std::fabs( t0 - tr.tau ) <= options.TTOL && m0.count( i ) ){
          std::ostringstream os; os << "transition at t = " << tr.tau << ": state " << _mX[i].name() << " is mapped twice";
          _extractError = os.str(); return false; }
        if( mp.count( i ) ){ std::ostringstream os; os << "transition at t = " << tr.tau << ": state " << _mX[i].name() << " is mapped twice";
          _extractError = os.str(); return false; }
        FFVar const g = subst_inputs( { subst1( _in_model_dag( tr.left[c] ) ) }, el )[0];
        _mIC[ks][i] = g;  mp[i] = g;
      }
      jumps.emplace_back( tr.tau, mp );
    }
  }

  // ---- functions -------------------------------------------------------------------------------------------------
  // Every INTEGRAL record (the integral of its source over the evolution domain) is carried as a quadrature.  ODESLV
  // quadratures RESTART at every stage, so the integral over the domain is the SUM of the quadrature's values at the
  // end of every stage -- added to the function at every stage 1.._ns (before 2026-09-28 only the last stage was:
  // on a multi-stage model the integral covered the last stage alone).  Every EVALUATION record is its source at the stage boundary tau.
  // Functions are then assigned per OUTPUT, in declaration order: an output that IS one record is that record; an
  // output that COMBINES records must be AFFINE in them, g = sum_j a_j c_j + b, with coefficients a_j free of the
  // records (they may involve parameters and constants) -- it becomes ONE function whose terms a_j c_j sit in the
  // stage map at tau_j (an evaluation) or at the last stage (an integral); ODESLV sums the stage contributions of a
  // function, for the value and for the forward and adjoint sensitivities alike.  A NONLINEAR combination is refused.
  _mFCT.assign( _ns + 1, std::map<size_t,FFVar>() );
  auto const& vDef = var_deferred();
  std::vector<FFVar> quad;
  for( auto const& C : vDef ){
    if( !C.accumulate ) continue;
    quad.push_back( _in_model_dag( C.source ) );  _mQ.push_back( C.input );
  }
  // the stage whose END is tau; an evaluation at any other time cannot be represented by ODESLV (functions are
  // evaluated at stage times) and is REFUSED below -- it used to fall back to the last stage, i.e. silently to tf
  auto stage_of = [&]( double const tau )->long {
    for( size_t k = 0; k <= _ns; ++k )               // TTOL, not an exact match: a snapped point must still find it
      if( std::fabs( _dT[k] - tau ) <= options.TTOL ) return (long)k;
    return -1; };
  for( auto const& C : vDef )
    if( !C.accumulate && stage_of( C.tau ) < 0 ){
      std::ostringstream os;  os << "an output evaluates at t = " << C.tau << ", which is not a stage time: ODESLV"
        " evaluates functions at stage times only -- add a stage boundary there (one more element in the domain)";
      _extractError = os.str();  return false;
    }
  auto add_to = [&]( size_t const stg, size_t const idx, FFVar const& term ){
    auto it = _mFCT[stg].find( idx );
    if( it == _mFCT[stg].end() ) _mFCT[stg][idx] = term;  else it->second = it->second + term; };
  // the term a * (value of record C) of function idx; a == nullptr for a unit coefficient.  An INTEGRAL with a
  // non-constant coefficient (a parameter or constant, never a state or time) is carried as its OWN quadrature of
  // a * source -- a int f dt = int a f dt -- because a function in which a parameter multiplies a quadrature
  // variable is not supported by the sensitivity machinery (its parameter derivative would need the quadrature).
  auto add_record = [&]( size_t const idx, t_DeferredValue const& C, FFVar const* a ){
    if( C.accumulate && a && !a->cst() ){
      FFVar const qa( FFModel::_dag, C.input.name() + "*" );
      quad.push_back( subst1( _in_model_dag( *a ) ) * _in_model_dag( C.source ) );  _mQ.push_back( qa );
      for( size_t k = 1; k <= _ns; ++k ) add_to( k, idx, qa );          // per-stage quadrature: sum over stages
      return;
    }
    FFVar v = C.accumulate? subst1( _in_model_dag( C.input ) )
            : subst_inputs( { subst1( _in_model_dag( C.source ) ) }, elem_of_point( C.tau, C.side ) )[0];
    if( !C.accumulate && C.side == FFDom::PLUS )            // at a transition: the state JUST AFTER the jump, x(tau^+) =
      for( auto const& [t0, mp] : jumps )                   // the map applied to x(tau^-) -- evaluated at the stage end
        if( std::fabs( t0 - C.tau ) <= options.TTOL && !mp.empty() ){
          std::vector<FFVar> tgt, rep;  for( auto const& [i, g] : mp ){ tgt.push_back( _mX[i] ); rep.push_back( g ); }
          v = FFModel::_dag->compose( std::vector<FFVar>{ v }, tgt, rep )[0];
        }
    FFVar const term = !a? v: a->cst()? FFVar( a->num().val() * v ): subst1( _in_model_dag( *a ) ) * v;
    if( C.accumulate ){ for( size_t k = 1; k <= _ns; ++k ) add_to( k, idx, term ); }   // sum over stages
    else add_to( (size_t)stage_of( C.tau ), idx, term ); };
  // the records whose inputs occur in expression e (full ids)
  auto records_in = [&]( FFVar const& e ){
    std::vector<size_t> r;
    if( !e.dag() ) return r;
    auto sg = e.dag()->subgraph( 1, &e );
    for( size_t k = 0; k < vDef.size(); ++k ){
      if( !vDef[k].input.dag() ) continue;
      bool in = false;
      for( auto const* op : sg.l_op ){ if( !op ) continue;
        for( auto const* v : op->varout ) if( v && v->id() == vDef[k].input.id() ){ in = true; break; }
        if( in ) break; }
      if( in ) r.push_back( k );
    }
    return r; };
  // An output DISTRIBUTED over the evolution direction (2026-10-06, WORKPLAN 1.7) is carried at the STAGE TIMES: one
  // function per kept stage time tau_k -- the output expression's value at the END of stage k (tau_k^-; the initial
  // state at k = 0) -- in time order, the mask selecting the times: ALL all of them, LB the initial time, UB the final
  // one, ALL-LB the stage ends (not the initial time), and so on.  Each row is an ordinary per-stage function, so the
  // values and the forward and adjoint sensitivities come from the machinery that already serves the point outputs.
  auto keep_stage = [&]( int const lim, size_t const k ){
    bool const lb = ( k == 0 ), ub = ( k == _ns );
    if( lim == FFDom::ALL )                         return true;
    if( lim == FFDom::ALL - FFDom::LB )             return !lb;
    if( lim == FFDom::ALL - FFDom::UB )             return !ub;
    if( lim == FFDom::ALL - FFDom::LB - FFDom::UB ) return !lb && !ub;
    if( lim == FFDom::LB )                          return lb;
    if( lim == FFDom::UB )                          return ub;
    return false; };
  std::vector<size_t> blk_first;
  for( auto const& fct : var_output() ){
    blk_first.push_back( _nf );
    if( fct.kind == FFModel::FctKind::DISTRIBUTED ){
      if( !fct.point.empty() || fct.grid.size() != 1 || fct.grid.begin()->first.id() != _evolution_dom_var.id() ){
        _extractError = "an output of an ODESLV can be distributed over the evolution direction only (" + fct.var.name() + ")";
        return false;
      }
      if( !records_in( fct.var ).empty() ){
        _extractError = "an output distributed over the evolution direction cannot involve an evaluation or an integral ("
                      + fct.var.name() + ")";
        return false;
      }
      int const lim = fct.grid.begin()->second;  size_t nrow = 0;
      for( size_t k = 0; k <= _ns; ++k ){
        if( !keep_stage( lim, k ) ) continue;
        add_to( k, _nf++, subst_inputs( { subst1( _in_model_dag( fct.var ) ) }, elem_of_point( _dT[k], FFDom::MINUS ) )[0] );
        ++nrow;
      }
      if( !nrow ){ _extractError = "the mask of a distributed output selects no stage time (" + fct.var.name() + ")";  return false; }
      continue;
    }
    int direct = -1;
    for( size_t k = 0; k < vDef.size(); ++k )
      if( vDef[k].input.dag() && fct.var.dag() && vDef[k].input.id() == fct.var.id() ){ direct = (int)k; break; }
    if( direct >= 0 ){ add_record( _nf++, vDef[direct], nullptr ); continue; }        // the output IS one record
    std::vector<size_t> const rec = records_in( fct.var );
    if( rec.empty() ){ _mFCT[_ns][ _nf++ ] = subst1( _in_model_dag( fct.var ) ); continue; }
    // a combination of records: split g = sum_j a_j c_j + b, and require every a_j (and b) free of the records
    FFGraph* const dag = mc::type_cast<FFGraph>( fct.var.dag() );
    if( !dag ){ _extractError = "an output's DAG is not an FFGraph (" + fct.var.name() + ")"; return false; }
    std::vector<FFVar> cin;  for( size_t k : rec ) cin.push_back( vDef[k].input );
    std::vector<FFVar> A, B;
    try{
      A = dag->FAD( std::vector<FFVar>{ fct.var }, cin );
      B = dag->compose( std::vector<FFVar>{ fct.var }, cin, std::vector<FFVar>( cin.size(), FFVar( 0. ) ) );
    }
    catch( ... ){ A.clear(); }
    bool affine = ( A.size() == cin.size() && B.size() == 1 );
    for( size_t j = 0; affine && j < A.size(); ++j ) if( !A[j].cst() && !records_in( A[j] ).empty() ) affine = false;
    if( affine && !B[0].cst() && !records_in( B[0] ).empty() ) affine = false;
    if( !affine ){
      _extractError = "an output combines evaluations/integrals NONLINEARLY (" + fct.var.name() + "): ODESLV supports an"
                      " output that is one evaluation or integral, or an AFFINE combination of them";
      return false;
    }
    size_t const idx = _nf++;
    for( size_t j = 0; j < rec.size(); ++j ) add_record( idx, vDef[ rec[j] ], &A[j] );
    if( !( B[0].cst() && B[0].num().val() == 0. ) )
      add_to( _ns, idx, B[0].cst()? FFVar( B[0].num().val() ): subst1( _in_model_dag( B[0] ) ) );
  }
  blk_first.push_back( _nf );                             // the block of every output among the functions
  _blkFct.clear();
  for( size_t k = 0; k + 1 < blk_first.size(); ++k ) _blkFct.push_back( { blk_first[k], blk_first[k+1] - blk_first[k] } );
  if( tv_inputs.empty() || quad.empty() )                 // after the outputs: they may have added quadratures
    _mQUAD.assign( 1, quad );
  else{                                                    // integrands with time-varying inputs: per stage
    _mQUAD.assign( _ns, std::vector<FFVar>() );
    for( size_t k = 0; k < _ns; ++k ) _mQUAD[k] = subst_inputs( quad, elem_of_stage( k ) );
  }
  _nq = _mQ.size();
  // The IC system's placeholders are minted HERE, once, and not where _mDOT is built: they must be the
  // LAST variables added to the model's dag.  Minting them mid-extraction -- before the parameters and
  // the initial rows are finalised -- segfaults FFODESLV's embedding, which is why the original code
  // created them inside _build_ic_system().  Hoisting them out of the per-call path is right; hoisting
  // them EARLIER is not.
  _mUpVar.clear();  _mInpParNdx.clear();
  if( _wants_initial_rows() ){
    _mDotVar.resize( _nx );  _mLftVar.resize( _nx );
    for( size_t i = 0; i < _nx; ++i ){
      _mDotVar[i] = FFModel::_dag->add_var( "ic_d"   + std::to_string(i) );
      _mLftVar[i] = FFModel::_dag->add_var( "ic_lft" + std::to_string(i) );
    }
    // one placeholder per time-varying input, standing for its value at the instant being solved at, and
    // the parameter indices of that input's levels so the value can be computed there
    for( size_t j = 0; j < _mTVInp.size(); ++j ){
      _mUpVar.push_back( FFModel::_dag->add_var( "ic_u" + std::to_string(j) ) );
      _mInpParNdx.push_back( parameter_index( _mTVInp[j] ) );
    }
  }
  if( !_build_ic_system() ) return false;
  return true;
}

inline
bool
BASE_DYNAMIC::_split_in
( FFVar const& res, FFVar const& tgt, FFVar& L, FFVar& R )
{
  // Put a fresh LEAF where the target stands, then differentiate with respect to it for L and set it to zero for
  // R, so both come from ONE graph.  A LEAF target (a state, for an INITIAL row) belongs to compose(); a NON-LEAF
  // one (the derivative node) to substitute().
  FFVar d( FFModel::_dag );
  std::vector<FFVar> const dep{ res }, tar{ tgt }, rep{ d };
  std::vector<FFVar> rd;
  bool const leaf = ( tgt.id().first == FFVar::VAR );
  try{
    rd = leaf? FFModel::_dag->compose( dep, tar, rep ): FFModel::_dag->substitute( dep, tar, rep );
    if( rd.empty() ) return false;
    std::vector<FFVar> const dRdd = FFModel::_dag->FAD( rd, std::vector<FFVar>{ d } );
    if( dRdd.empty() ) return false;
    L = _in_model_dag( dRdd[0] );
    std::vector<FFVar> const r0 = FFModel::_dag->compose( rd, std::vector<FFVar>{ d },
                                                          std::vector<FFVar>{ FFVar( 0. ) } );
    if( r0.empty() ) return false;
    R = _in_model_dag( r0[0] );
  }
  catch( ... ){ return false; }
  return true;
}

inline
bool
BASE_DYNAMIC::_build_ic_system
()
{
  _icStart.ok = _icStop.ok = false;
  if( _mRESraw.empty() ) return true;
  // BOTH systems, once.  They share every row but the last block: the start takes the INITIAL rows, a stop
  // takes CONTINUITY on the differential states with the left state as data.  Built from the RAW rows, so
  // the time-varying inputs are still present to be replaced by _mUpVar and bound per stop -- that, and
  // nothing else, is what made this per-stage before.
  bool ok = _assemble_ic( _icStart, false );
  if( ok && _wants_initial_rows() ) ok = _assemble_ic( _icStop, true );
  return ok;
}

inline
bool
BASE_DYNAMIC::_assemble_ic
( t_ICSys& S, bool const at_stop )
{
  delete S.dag;  S.dag = new FFGraph;
  S.row.clear(); S.unk.clear(); S.dotOf.clear(); S.eval.clear();
  if( !at_stop && _mICRES.empty() ) return true;

  // 1. the rows, in the MODEL's graph: the undivided evolution rows, then the closing block.
  std::vector<FFVar> rows;
  for( size_t k = 0; k < _mRESraw.size(); ++k ){
    FFVar r = _mRESraw[k];
    for( auto const& [ ix, co ] : _mMASSraw[k] ) r += co * _mDOT[ix];   // M.dx + R, undivided
    rows.push_back( r );
  }
  size_t const nDyn = rows.size();
  if( !at_stop ) rows.insert( rows.end(), _mICRES.begin(), _mICRES.end() );
  else
    for( size_t i = 0; i < _nx; ++i )
      if( _differential_state( i ) ) rows.push_back( _mX[i] - _mLftVar[i] );

  // 2. the time-varying inputs become the placeholders u(tau), bound per stop.  compose(), not
  //    substitute(): an input is a plain VARIABLE, where a derivative is an auxiliary node.
  if( !_mTVInp.empty() ){
    try{ rows = FFModel::_dag->compose( rows, _mTVInp, _mUpVar ); }
    catch( ... ){ _extractError = "parameterising the time-varying inputs of the initial-condition rows failed";
                  return false; }
  }

  // 3. substitute every derivative node for a plain variable -- FFGraph::substitute takes auxiliary
  //    nodes, which a derivative is.  After this the rows are ordinary expressions.
  std::vector<FFVar> plain;
  try{ plain = FFModel::_dag->substitute( rows, _mDOT, _mDotVar ); }
  catch( ... ){ _extractError = "substituting the state derivatives in the initial-condition rows failed";
                return false; }

  // 4. which derivatives actually appear?  Those are the extra unknowns; an algebraic state has none.
  std::vector<bool> used( _nx, false );
  std::vector<FFVar> const dPh = FFModel::_dag->FAD( plain, _mDotVar );
  for( size_t r = 0; r < plain.size(); ++r )
    for( size_t i = 0; i < _nx; ++i ){
      FFVar const& d = dPh[ r*_nx + i ];
      if( !( d.cst() && d.num().val() == 0. ) ) used[i] = true;
    }

  // 5. move it all into this system's OWN graph
  std::vector<FFVar> srcV{ _mT[0] }, dstV( 1 );
  for( size_t i = 0; i < _nx; ++i ) srcV.push_back( _mX[i] );
  for( size_t j = 0; j < _np; ++j ) srcV.push_back( _mP[j] );
  for( size_t j = 0; j < _nc; ++j ) srcV.push_back( _mC[j] );
  for( auto const& v : _mUpVar ) srcV.push_back( v );
  for( size_t i = 0; i < _nx; ++i ) if( used[i] ) srcV.push_back( _mDotVar[i] );
  if( at_stop ) for( size_t i = 0; i < _nx; ++i ) srcV.push_back( _mLftVar[i] );
  dstV.assign( srcV.size(), FFVar() );
  S.dag->insert( FFModel::_dag, srcV.size(), srcV.data(), dstV.data() );
  S.all = dstV;
  S.row.assign( plain.size(), FFVar() );
  S.dag->insert( FFModel::_dag, plain.size(), plain.data(), S.row.data() );

  size_t const nU = _mUpVar.size();
  for( size_t i = 0; i < _nx; ++i ) S.unk.push_back( dstV[ 1 + i ] );
  { size_t at = 1 + _nx + _np + _nc + nU;
    for( size_t i = 0; i < _nx; ++i ) if( used[i] ){ S.unk.push_back( dstV[at++] ); S.dotOf.push_back( i ); } }
  S.nLft = at_stop? _nx: 0;

  if( S.row.size() != S.unk.size() ){
    _extractError = "the consistent-initialisation system is not square: " + std::to_string( S.row.size() )
                  + " row(s) (" + std::to_string( nDyn ) + " dynamic + "
                  + std::to_string( S.row.size() - nDyn ) + ( at_stop? " continuity": " initial" )
                  + ") for " + std::to_string( S.unk.size() ) + " unknown(s)";
    S.row.clear(); S.unk.clear();
    return true;                                   // not fatal: the consumer reports it if it needs the solve
  }
  // SFAD, not FAD: the dense form builds n^2 expressions and evaluates them all, most of them structural
  // zeros.  This mirrors what the integrator already does for its own Jacobian, and is what lets a sparse
  // linear solver be offered under KINSOL rather than only a dense one.
  S.eval = S.row;
  { auto const Js = S.dag->SFAD( S.row.size(), S.row.data(), S.unk.size(), S.unk.data() );
    S.nnz = std::get<0>( Js );
    S.jrow.assign( S.nnz, 0 ); S.jcol.assign( S.nnz, 0 );
    for( unsigned e = 0; e < S.nnz; ++e ){
      S.jrow[e] = std::get<1>( Js )[e];
      S.jcol[e] = std::get<2>( Js )[e];
    }
    S.eval.insert( S.eval.end(), std::get<3>( Js ), std::get<3>( Js ) + S.nnz );
    delete[] std::get<1>( Js ); delete[] std::get<2>( Js ); delete[] std::get<3>( Js );
    // CSC, as the integrator's KLU path builds it, plus where each SFAD entry lands in the value array
    size_t const n = S.unk.size();
    std::vector< std::pair<size_t,size_t> > key;        // (col, row) ordered
    for( size_t e = 0; e < S.nnz; ++e ) key.emplace_back( S.jcol[e], S.jrow[e] );
    std::vector<size_t> ord( S.nnz );
    for( size_t e = 0; e < S.nnz; ++e ) ord[e] = e;
    std::sort( ord.begin(), ord.end(), [&]( size_t a, size_t b ){ return key[a] < key[b]; } );
    S.sp.colptr.assign( n + 1, 0 );
    S.sp.rowidx.assign( S.nnz, 0 );
    S.at.assign( S.nnz, 0 );
    for( size_t k2 = 0; k2 < S.nnz; ++k2 ){
      size_t const e = ord[k2];
      S.sp.rowidx[k2] = (sunindextype)S.jrow[e];
      S.at[e] = k2;
      ++S.sp.colptr[ S.jcol[e] + 1 ];
    }
    for( size_t c = 0; c < n; ++c ) S.sp.colptr[c+1] += S.sp.colptr[c];
    S.sp.jac = &BASE_DYNAMIC::_ic_spjac; }
  S.dEval.assign( S.eval.size(), 0. );
  S.dVar.assign( 1 + _nx + _np + _nc + nU + S.dotOf.size() + S.nLft, 0. );
  S.ok = true;
  return true;
}

inline
bool
BASE_DYNAMIC::explicit_initial_conditions
()
const
{
  if( _mICRES.empty() ) return false;
  std::vector<FFVar> const dX = FFModel::_dag->FAD( _mICRES, _mX );
  std::vector<FFVar> const dD = FFModel::_dag->FAD( _mICRES, _mDOT );
  for( size_t r = 0; r < _mICRES.size(); ++r ){
    size_t nst = 0;
    for( size_t i = 0; i < _nx; ++i ){
      FFVar const& d = dX[ r*_nx + i ];
      if( !( d.cst() && d.num().val() == 0. ) ){ if( !d.cst() ) return false; ++nst; }
      FFVar const& e = dD[ r*_nx + i ];
      if( !( e.cst() && e.num().val() == 0. ) ) return false;      // carries a derivative
    }
    if( nst != 1 ) return false;                                   // not one state, pinned linearly
  }
  return true;
}

inline
BASE_KINSOL::Report
BASE_DYNAMIC::solve_ic
( std::vector<double>& x, std::vector<double>& dx, BASE_KINSOL::Knobs const& k,
  std::vector<double> const& par )
{
  _icLeft.clear();
  return _ic_solve( _icStart, _dT.empty()? 0.: _dT[0], x, dx, k, par );
}

inline
BASE_KINSOL::Report
BASE_DYNAMIC::_ic_solve
( t_ICSys& S, double const tau, std::vector<double>& x, std::vector<double>& dx,
  BASE_KINSOL::Knobs const& k, std::vector<double> const& par )
{
  BASE_KINSOL::Report R;
  if( S.row.empty() ){
    R.message = _extractError.empty()? "no consistent-initialisation system was built": _extractError;
    return R;
  }
  size_t const nU = _mUpVar.size();
  S.dVar[0] = tau;                      // the instant this system is being solved at, NOT always t0
  for( size_t j = 0; j < _np && j < par.size(); ++j ) S.dVar[ 1 + _nx + j ] = par[j];
  for( size_t j = 0; j < _nc && j < _dIcCst.size(); ++j ) S.dVar[ 1 + _nx + _np + j ] = _dIcCst[j];

  // Bind each time-varying input to its value AT tau.  This is the whole of what used to require a system
  // per stage: the geometry lives in computing these numbers, never in the rows.  u(tau^+) -- the element
  // STARTING at tau -- because the rows being solved are those in force after the stop.
  if( nU ){
    size_t const e = _elem_of( tau, FFDom::PLUS );
    double const te = _mElemBnd.at( e ), he = _mElemBnd.at( e+1 ) - te;
    double const sc = he > 0.? ( tau - te ) / he: 0.;
    for( size_t j = 0; j < nU; ++j ){
      auto const& sn  = _mInpNodes[ _mTVInp[j].name() ];
      auto const& ndx = _mInpParNdx[j];
      size_t const nn = sn.size();
      std::vector<double> nodal;
      if( nn < 2 ){ if( e < ndx.size() ) nodal.push_back( par.at( ndx[e] ) ); }
      else for( size_t q = 0; q < nn && e*nn+q < ndx.size(); ++q ) nodal.push_back( par.at( ndx[e*nn+q] ) );
      S.dVar[ 1 + _nx + _np + _nc + j ] = nodal.empty()? 0.: _interp_level( sn, nodal, sc );
    }
  }

  size_t const dot0 = 1 + _nx + _np + _nc + nU;
  { std::vector<double> const ref = _ic_seed();
    for( size_t i = 0; i < _nx; ++i ) S.dVar[1+i] = _icLeft.size() == _nx? _icLeft[i]: ref[i];
    for( size_t d = 0; d < S.dotOf.size(); ++d ) S.dVar[ dot0 + d ] = 0.;
    // the continuity rows read the LEFT state as data; seeding from it too is the best start there is
    for( size_t i = 0; i < S.nLft && i < _icLeft.size(); ++i )
      S.dVar[ dot0 + S.dotOf.size() + i ] = _icLeft[i]; }

  std::vector<double> z( S.unk.size() );
  for( size_t i = 0; i < _nx; ++i ) z[i] = S.dVar[1+i];
  for( size_t d = 0; d < S.dotOf.size(); ++d ) z[_nx+d] = S.dVar[ dot0 + d ];

  _icCur = &S;
  R = BASE_KINSOL::solve( (int)z.size(), &BASE_DYNAMIC::_ic_res, &BASE_DYNAMIC::_ic_jac, this, z, k, &S.sp );
  _icCur = nullptr;

  x.assign( _nx, 0. ); dx.assign( _nx, 0. );
  for( size_t i = 0; i < _nx; ++i ) x[i] = z[i];
  for( size_t d = 0; d < S.dotOf.size(); ++d ) dx[ S.dotOf[d] ] = z[_nx+d];
  return R;
}

//! Scatter the KINSOL iterate into the current system's value array: the states first, then the
//! derivative unknowns, which sit after the parameters, constants and input placeholders.
inline
void
BASE_DYNAMIC::_ic_scatter
( BASE_DYNAMIC* M, int const n, double const* z )
{
  t_ICSys& S = *M->_icCur;
  size_t const dot0 = 1 + M->_nx + M->_np + M->_nc + M->_mUpVar.size();
  for( int i = 0; i < n; ++i ) S.dVar[ i < (int)M->_nx? 1+i: dot0 + ( i - M->_nx ) ] = z[i];
}

inline
int
BASE_DYNAMIC::_ic_res
( int const n, double const* z, double* f, void* user )
{
  BASE_DYNAMIC* M = static_cast<BASE_DYNAMIC*>( user );
  t_ICSys& S = *M->_icCur;
  _ic_scatter( M, n, z );
  S.dag->eval( S.wrk, S.row.size(), S.eval.data(), S.dEval.data(),
               S.all.size(), S.all.data(), S.dVar.data() );
  for( int i = 0; i < n; ++i ) f[i] = S.dEval[i];
  return 0;
}

inline
int
BASE_DYNAMIC::_ic_jac
( int const n, double const* z, double* J, void* user )
{
  BASE_DYNAMIC* M = static_cast<BASE_DYNAMIC*>( user );
  t_ICSys& S = *M->_icCur;
  _ic_scatter( M, n, z );
  S.dag->eval( S.wrk, S.eval.size(), S.eval.data(), S.dEval.data(),
               S.all.size(), S.all.data(), S.dVar.data() );
  size_t const nr = S.row.size();
  for( int i = 0; i < n*n; ++i ) J[i] = 0.;                       // scatter the SFAD entries, column-major
  for( size_t e = 0; e < S.nnz; ++e )
    J[ S.jrow[e] + (size_t)n * S.jcol[e] ] = S.dEval[ nr + e ];
  return 0;
}

inline
int
BASE_DYNAMIC::_ic_spjac
( int const n, double const* z, double* val, void* user )
{
  BASE_DYNAMIC* M = static_cast<BASE_DYNAMIC*>( user );
  t_ICSys& S = *M->_icCur;
  _ic_scatter( M, n, z );
  S.dag->eval( S.wrk, S.eval.size(), S.eval.data(), S.dEval.data(),
               S.all.size(), S.all.data(), S.dVar.data() );
  size_t const nr = S.row.size();
  for( size_t e = 0; e < S.nnz; ++e ) val[ S.at[e] ] = S.dEval[ nr + e ];
  return 0;
}

} // end namespace mc

#endif
