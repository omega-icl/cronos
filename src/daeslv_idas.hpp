// Copyright (C) Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the EPL-2.0 with GPL-2.0-or-later as a Secondary License; see the LICENSE file.
// SPDX-License-Identifier: EPL-2.0 OR GPL-2.0-or-later

#ifndef CRONOS__DAESLV_IDAS_HPP
#define CRONOS__DAESLV_IDAS_HPP

#include <map>
#include <set>
#include <string>
#include <tuple>
#include <vector>

#include "daeslv_base.hpp"
#include "base_idas.hpp"
#include "base_kinsol.hpp"

namespace mc
{
//! @brief C++ class integrating a model's DAE form with SUNDIALS IDAS
////////////////////////////////////////////////////////////////////////
//! mc::DAESLV_IDAS is the IDAS counterpart of mc::ODESLV_CVODES.  It takes the general mass-matrix form
//! that DAESLV_BASE extracts,
//!
//!     sum over _mMASS[0][k] of ( coefficient * dy[index] )  +  _mRES[0][k]  =  0 ,
//!
//! and hands it to IDAS as the residual, with the id vector marking which unknowns are differential.
//!
//! WHAT v0 DOES NOT DO YET, each deliberate rather than overlooked:
//!   * ONE SET OF ROWS -- a model with TIME-VARYING INPUTS is refused by setup(), since its rows differ
//!     stage by stage; without them the _ns stages hold identical rows and one set covers the horizon;
//!   * DIFFERENCE-QUOTIENT Jacobian (IDAS supplies it when no IDALsJacFn is set).  The analytic one is a
//!     later step; the residual is symbolic, so it is available, just not wired;
//!   * no quadratures, no point outputs, no sensitivities.  The eight BASE_IDAS sensitivity callbacks are
//!     present and REFUSE, so phases 5/6 fill them in rather than finding a silent wrong answer;
//!   * consistent initialisation is IDACalcIC (option A of the plan).  Option C -- a one-element OCFESLV
//!     solve -- is the fallback for multi-rooted ICs and is not wired here.
//!
//! The residual is evaluated in the MODEL's dag rather than a solver-local copy, so nothing is inserted
//! or renamed; the coefficients and remainders come out of one subgraph in a single eval, and the
//! derivative values are folded in afterwards.  That keeps v0 small; a solver-local dag is what
//! ODESLV_CVODES does, and is the thing to revisit when threads or a second concurrent solve appear.
////////////////////////////////////////////////////////////////////////
class DAESLV_IDAS
: public virtual DAESLV_BASE
, public virtual BASE_IDAS
{
public:
  //! @brief Tag the DAG operation prints: "DAE[...]" and "GRAD_DAE[...]".
  static constexpr char const* FFTAG = "DAE";

  DAESLV_IDAS()
    : _nres(0), _ncoef(0)
    {}

  //! @brief Construct on the user DAG @p dag, as ODESLV_CVODES does.  Equivalent to the default
  //! constructor followed by set_dag(); it exists because consumers -- the Python binder, fdiff's
  //! factory -- construct a solver from a dag and should not have to know which solver they hold.
  DAESLV_IDAS
    ( FFGraph* dag )
    : DAESLV_IDAS()
    { set_dag( dag ); }
  virtual ~DAESLV_IDAS()
    { unregistration(); delete _dag; }

  //! @brief The solver's options: FFModel's groups (DISPLAY_LEVEL, CLASSIFY, REDUCE, FATAL, TTOL) together
  //! with the integrator's SOLVE and INITIAL.  This SHADOWS FFModel::options, exactly as OCFESLV's does, and
  //! _apply_options() is the single moment the two are reconciled -- without that override a caller's
  //! REDUCE/CLASSIFY settings would be written here and never read by the model.
  struct Options
  : public FFModel::Options
  , public BASE_IDAS::Options
  {} options;


  //! @brief Statistics of the last solve
  Stats stats;          //!< the LIVE counters, incremented by the callbacks during a solve

  //! @brief Statistics of the last plain solve / solve_fsens / solve_asens, as ODESLV keeps them:
  //! one per solve KIND, so a forward run does not overwrite what the adjoint reported.
  Stats stats_solve;
  Stats stats_fsens;
  Stats stats_asens;

  //! @brief Apply the caller's options to the model before setup() reads any of them (ocfeslv.hpp:1383).
  void _apply_options
    ()
    override
    {
      FFModel::Options::_force_display_from_env( options.DISPLAY_LEVEL );
      FFModel::options = static_cast<FFModel::Options const&>( options );
    }

  //! @brief Extract the model as a DAE and build the evaluation plan.  False with err() set.
  bool setup
    ();

  //! @brief How many of the np() parameters are per-stage INPUT LEVELS rather than model parameters.  They
  //! are passed in solve()'s parameter vector like any other.  Non-zero exactly when the model has
  //! time-varying inputs, which is also when the rows differ stage by stage and the integrator restarts.
  size_t n_input_levels
    ()
    const
    { return _vInpLvl.size(); }

  //! @brief Computes solution of parametric DAEs over the DECLARED horizon.  Mirrors
  //! ODESLV_CVODES::solve: the model's domain already fixes the horizon, so no output times are passed --
  //! the state at every stage time is left in val_state().
  STATUS solve
    ( std::vector<double> const& p, std::vector<double> const& c=std::vector<double>(),
      std::ostream& os=std::cout );

  //! @brief Computes solution of parametric DAEs
  STATUS solve
    ( double const* p, double const* c=nullptr, std::ostream& os=std::cout )
    { std::vector<double> P, C;                 // assign() rather than the range constructor: with p or c
      if( p && _np ) P.assign( p, p + _np );    // null, or the count zero, the constructor forms a pointer
      if( c && _nc ) C.assign( c, c + _nc );    // GCC cannot prove is in bounds, and warns under -Wall
      return solve( P, C, os ); }

  //! @brief Solve with the input and constant values given BY NAME, as ODESLV_CVODES does:
  //! { {u, {v1,...,vn}}, ..., {K, {k}} }.  assemble_values() validates the list and names the offender.
  STATUS solve
    ( std::vector<FFModel::InputVal> const& vIn, std::ostream& os=std::cout )
    {
      std::vector<double> P, c;  std::string err;
      if( !assemble_values( vIn, P, c, err ) ) throw std::invalid_argument( "DAESLV::solve ** " + err + "\n" );
      return solve( P, c, os );
    }

  //! @brief Solve the consistent-initialisation system.  A thin forward to the SHARED solve in
  //! BASE_DE -- the system is the same for an ODE and a DAE, so it is built and solved there; this
  //! only supplies the knobs from INITIAL.KIN and the current parameter values.
  BASE_KINSOL::Report solve_ic
    ( std::vector<double>& y, std::vector<double>& yp )
    { std::vector<double> par( _np, 0. );
      for( size_t j = 0; j < _np && 1 + 2*_nx + j < _dVAR.size(); ++j ) par[j] = _dVAR[ 1 + 2*_nx + j ];
      return BASE_DE::solve_ic( y, yp, options.INITIAL.KIN, par ); }

  //! @brief DAESLV implements transitions: a jump is a ROW of the stop system, replacing the continuity
  //! of the state it determines.  Not lifted -- the stop IS where it is applied -- so the times must
  //! still be validated, or a tau between two stage boundaries would be silently ignored (WORKPLAN 1.6).
  bool _transitions_supported
    ()
    const
    override
    { return true; }

  bool _validate_transition_tau_unlifted
    ()
    const
    override
    { return true; }

  bool _transitions_in_ic_system
    ()
    const
    override
    { return true; }

  bool _transition_derivatives_supported
    ()
    const
    override
    { return true; }

  //! @brief A tau must coincide with an ELEMENT BOUNDARY, which is where a stop already exists.  Same
  //! policy and the same words as OCFESLV, so the two solvers refuse the same models for the same reason.
  bool _accept_transition_tau
    ( double const tau, std::string& why )
    const
    override
    {
      auto const it = _issetup? _mDom.find( _evolution_dom_var )
                              : _usr._mDomUsr.find( _usr._evolution_dom_varUsr );
      auto const& M = _issetup? _mDom: _usr._mDomUsr;
      if( it == M.end() ){ why = "no evolution domain"; return false; }
      for( size_t k = 1; k < it->second.n_elem; ++k )
        if( std::fabs( it->second.elem_bnd[k] - tau ) <= 64.*DBL_EPSILON*std::max( 1., std::fabs( tau ) ) )
          return true;
      why = "tau must be an ELEMENT BOUNDARY of the evolution direction (add one there)";
      return false;
    }

  //! @brief The IC knobs, from this solver's INITIAL group.
  BASE_KINSOL::Knobs _ic_knobs
    ()
    const
    override
    { return options.INITIAL.KIN; }

  //! @brief Whether the explicit fast path may be taken; see BASE_DE's hook.
  bool _ic_bypass_explicit
    ()
    const
    override
    { return options.INITIAL.BYPASS_EXPLICIT; }

  //! @brief Nonzeros of the Jacobian IDAS factorizes -- the UNION of the dF/dy and dF/dy' patterns.
  size_t nnz_jac
    ()
    const
    { return _nnzjac; }

  //! @brief Retreive state values at stage times -- one entry for the initial point and one per stage, as
  //! ODESLV_CVODES's _xk is filled, so val_state().back() is the state at the end of the horizon.
  //! @brief The output functions, in add_output() order via function_block().  A term is accumulated at
  //! each stage time it belongs to, so a distributed output contributes once per stage.
  std::vector<double> const& val_function
    ()
    const
    { return _dF; }

  //! @brief Write the recorded STATE trajectory: one row per stage time, [ state, quadrature ].
  //! One writer per results_* container, named for it -- record_solve/record_fsens/record_asens.
  void record_solve
    ( std::ofstream& ores, unsigned const iprec=5 )
    const
    { _record( ores, results_solve, iprec ); }

  //! @brief Write the recorded SENSITIVITY trajectories, one stream per parameter.
  void record_fsens
    ( std::ofstream* osen, unsigned const iprec=5 )
    const
    { if( !osen ) return;
      for( size_t j = 0; j < results_fsens.size(); ++j ) _record( osen[j], results_fsens[j], iprec ); }

  //! @brief Write the recorded ADJOINT trajectories, one stream per output, in DESCENDING time.
  void record_asens
    ( std::ofstream* oadj, unsigned const iprec=5 )
    const
    { if( !oadj ) return;
      for( size_t r = 0; r < results_asens.size(); ++r ){
        // the container holds the order the BACKWARD sweep produced -- descending.  The FILE reads in
        // increasing time, like every other recorded trajectory, so the writer reverses rather than the
        // sweep storing against its own grain.
        std::vector<Results> asc( results_asens[r].rbegin(), results_asens[r].rend() );
        _record( oadj[r], asc, iprec );
      } }


  // ------------------------------------------------------------------------------------------
  // fdiff: the solver as a DIRECTIONAL DERIVATIVE of itself.  Mirrored from ODESLV_CVODES on
  // 2026-10-22, returning DAESLV_IDAS* as agreed.  The work is FFModel::fdiff's -- these forms
  // map the caller's parameters back to the INPUTS they belong to, seed one direction per
  // parameter, and carry constants through with their own seeds.  This is what FFDAESLV needs.
  // ------------------------------------------------------------------------------------------
  //! @brief The solver as a DIRECTIONAL DERIVATIVE of itself.  The six forms live in BASE_DE
  //! (_fdiff); these supply the factory, which is all that differed between the two solvers.
  DAESLV_IDAS* fdiff( std::vector<FFVar> const& vU, size_t const nDir,
      std::vector<std::vector<FFVar>>& vDU, std::vector<std::vector<FFVar>>& vS, std::string& err ) const
    { return _fdiff<DAESLV_IDAS>( [this]{ auto* s = new DAESLV_IDAS( _usr._dagUsr ); s->options = options; return s; }, vU, nDir, vDU, vS, err ); }
  DAESLV_IDAS* fdiff( size_t const nPar, FFVar const* pPar, std::string& err ) const
    { return _fdiff<DAESLV_IDAS>( [this]{ auto* s = new DAESLV_IDAS( _usr._dagUsr ); s->options = options; return s; }, nPar, pPar, err ); }
  DAESLV_IDAS* fdiff( size_t const nPar, FFVar const* pPar ) const
    { return _fdiff<DAESLV_IDAS>( [this]{ auto* s = new DAESLV_IDAS( _usr._dagUsr ); s->options = options; return s; }, nPar, pPar ); }
  DAESLV_IDAS* fdiff( std::vector<FFVar> const& vPar ) const
    { return _fdiff<DAESLV_IDAS>( [this]{ auto* s = new DAESLV_IDAS( _usr._dagUsr ); s->options = options; return s; }, vPar ); }
  DAESLV_IDAS* fdiff( std::vector<FFVar> const& vU, std::string& err ) const
    { return _fdiff<DAESLV_IDAS>( [this]{ auto* s = new DAESLV_IDAS( _usr._dagUsr ); s->options = options; return s; }, vU, err ); }
  DAESLV_IDAS* fdiff( std::vector<std::vector<FFVar>>& vDU, std::vector<std::vector<FFVar>>& vS,
                std::string& err ) const
    { return _fdiff<DAESLV_IDAS>( [this]{ auto* s = new DAESLV_IDAS( _usr._dagUsr ); s->options = options; return s; }, vDU, vS, err ); }

  //! @brief Adjoint values at every stage time, [stage][output][lambda], as ODESLV exposes them.
  std::vector< std::vector< std::vector<double> > > const& val_state_adjoint
    ()
    const
    { return _lk; }

  //! @brief Quadrature sensitivities at every stage time, [ time ][ direction ][ quadrature ], as
  //! ODESLV_CVODES exposes them.  Filled by solve_fsens where the state sensitivities are.
  std::vector< std::vector< std::vector<double> > > const& val_quadrature_sensitivity
    ()
    const
    { return _qpk; }

  //! @brief Gradient of the outputs w.r.t. the parameters, _nf rows of _np, after solve_fsens().
  std::vector< std::vector<double> > const& val_function_gradient
    ()
    const
    { return _fp; }

  //! @brief Forward sensitivities with the inputs given BY NAME, as ODESLV_CVODES has it.
  STATUS solve_fsens
    ( std::vector<FFModel::InputVal> const& vIn, std::ostream& os=std::cout )
    { std::vector<double> P, C;  std::string err;
      if( !assemble_values( vIn, P, C, err ) )
        throw std::invalid_argument( "DAESLV::solve_fsens ** " + err );
      return solve_fsens( P, C, os ); }

  //! @brief ... and the adjoint likewise.
  STATUS solve_asens
    ( std::vector<FFModel::InputVal> const& vIn, std::ostream& os=std::cout )
    { std::vector<double> P, C;  std::string err;
      if( !assemble_values( vIn, P, C, err ) )
        throw std::invalid_argument( "DAESLV::solve_asens ** " + err );
      return solve_asens( P, C, os ); }

  //! @brief The gradient rows of control @p var: [dof][function], its DOFs in control_dofs() order, from the last
  //! solve_fsens() or solve_asens().  Throws if @p var is not a registered control (it then has no
  //! sensitivity direction) or no sensitivity solve has been run.
  std::vector< std::vector< double > > val_function_gradient
    ( FFVar const& var )
    const
    {
      auto const b = control_block( var );
      auto const& G = val_function_gradient();
      if( !b.ndof )
        throw std::invalid_argument( "DAESLV::val_function_gradient ** " + var.name() + " is not a registered control\n" );
      if( G.size() < b.offset + b.ndof )
        throw std::invalid_argument( "DAESLV::val_function_gradient ** no sensitivity solve over the current controls\n" );
      return std::vector< std::vector< double > >( G.begin()+b.offset, G.begin()+b.offset+b.ndof );
    }

  //! @brief The gradient rows of control @p var on ONE element: [node][function], the n_node DOFs of element
  //! @p ndx_el (keyed by the evolution domain, as pos_input takes it).  A time-invariant control ignores @p ndx_el.
  std::vector< std::vector< double > > val_function_gradient
    ( FFVar const& var, std::map<FFVar,size_t,lt_FFVar> const& ndx_el )
    const
    {
      auto const all = val_function_gradient( var );
      size_t const nn = size_input( var );
      if( all.size() == 1 ) return all;                                   // time-invariant
      size_t e = 0;
      if( !_mT.empty() ){ auto const it = ndx_el.find( _mT[0] ); if( it != ndx_el.cend() ) e = it->second; }
      if( (e+1)*nn > all.size() )
        throw std::invalid_argument( "DAESLV::val_function_gradient ** element out of range for " + var.name() + "\n" );
      return std::vector< std::vector< double > >( all.begin()+e*nn, all.begin()+(e+1)*nn );
    }

  //! @brief Forward sensitivities at every stage time: [time][direction][state].  Empty unless the last
  //! call was solve_fsens().
  std::vector< std::vector< std::vector<double> > > const& val_state_sensitivity
    ()
    const
    { return _sk; }

  //! @brief Setup as an INDEPENDENT COPY of @p IVP: the same description, this solver's options, then
  //! a full setup().  The user DAG @p IVP was declared on must outlive this solver.  As ODESLV's.
  bool setup
    ( DAESLV_IDAS const& IVP )
    { std::lock_guard<std::recursive_mutex> dag_lock_( dag_mutex() );  _copy_description( IVP );
      options = IVP.options;  return setup(); }

  //! @brief Pointer forms of the sensitivity solves, which the FF operations call.  They forward to
  //! the vector forms; the sizes are np() and the number of constants.
  STATUS solve_fsens
    ( double const* p, double const* c = nullptr, std::ostream& os = std::cout )
    { return solve_fsens( std::vector<double>( p, p + np() ),
                          c? std::vector<double>( c, c + var_constant().size() ): std::vector<double>(), os ); }
  STATUS solve_asens
    ( double const* p, double const* c = nullptr, std::ostream& os = std::cout )
    { return solve_asens( std::vector<double>( p, p + np() ),
                          c? std::vector<double>( c, c + var_constant().size() ): std::vector<double>(), os ); }

  //! @brief Integrate WITH forward sensitivities w.r.t. every parameter.  The sensitivity initial data
  //! comes from ic_sens() -- the same system the state's did -- so the two cannot disagree.
  STATUS solve_fsens
    ( std::vector<double> const& p, std::vector<double> const& c = std::vector<double>(),
      std::ostream& os = std::cout )
    { _fsa = true;  STATUS const st = solve( p, c, os );  _fsa = false;  return st; }

  //! @brief Integrate and obtain the gradient of INTEGRAL output @p ifct by ADJOINT sensitivity.
  //! The result lands in val_function_gradient(), where solve_fsens() puts the forward one, so the
  //! two are directly comparable.
  //! @brief The adjoint of ONE output.  solve_asens() repeats it over every output; call it directly
  //! when only one gradient is wanted and nf is large.
  STATUS _solve_asens_one
    ( std::vector<double> const& p, size_t const ifct,
      std::vector<double> const& c = std::vector<double>(), std::ostream& os = std::cout )
    { _asa = true;  _asaF = ifct;  _asaQ = (size_t)-1;
      STATUS const st = solve( p, c, os );
      _asa = false;  _asaF = (size_t)-1;  _asaQ = (size_t)-1;  return st; }

  //! @brief The adjoint of EVERY output, in one call -- ODESLV_CVODES's signature and behaviour.
  //! One backward sweep serves one output, so this runs nf of them and accumulates their gradients
  //! into val_function_gradient(); _solve_asens_one is the single-output step it repeats.
  STATUS solve_asens
    ( std::vector<double> const& p, std::vector<double> const& c = std::vector<double>(),
      std::ostream& os = std::cout )
    { std::vector< std::vector<double> > acc;
      std::vector< std::vector< std::vector<double> > > lk;
      STATUS st = NORMAL;
      for( size_t r = 0; r < ( _nf? _nf: 1 ); ++r ){
        st = _solve_asens_one( p, r, c, os );
        if( st != NORMAL ) return st;
        if( acc.empty() ) acc.assign( _np, std::vector<double>( _nf, 0. ) );
        for( size_t j = 0; j < _np && j < _fp.size(); ++j )
          if( r < _fp[j].size() ) acc[j][r] = _fp[j][r];     // keep this output's column
        if( lk.empty() ) lk = _lk; else for( size_t k = 0; k < lk.size() && k < _lk.size(); ++k )
          for( size_t f = 0; f < lk[k].size() && f < _lk[k].size(); ++f )
            if( f == r ) lk[k][f] = _lk[k][f];
      }
      _fp.swap( acc );  _lk.swap( lk );
      return st; }

  //! @brief Quadratures at every stage time, one row per time, _nq wide.  Each row is that STAGE's
  //! increment, not a running total: a quadrature restarts at zero on every stage, and summing them is
  //! val_function()'s job.  Row 0 is zero.
  std::vector< std::vector<double> > const& val_quadrature
    ()
    const
    { return _qk; }

  std::vector< std::vector<double> > const& val_state
    ()
    const
    { return _xk; }

  //! @brief The state vector right after consistent initialisation, before any step.
  std::vector<double> const& val_state_init
    ()
    const
    { return _dX0; }

  //! @brief The derivative vector right after consistent initialisation.
  std::vector<double> const& val_deriv_init
    ()
    const
    { return _dXP0; }

protected:

  size_t _nres;                              //!< rows of the residual (== _nx)
  size_t _ncoef;                             //!< mass-matrix entries actually present
  size_t _offJy, _offJyp;                    //!< where the dF/dy and dF/dy' entries start in _dEVAL
  //! SOLVER-LOCAL dag, as ODESLV_BASE::_SETUP builds one: the model's description is INSERTED here and
  //! everything downstream -- residual, Jacobian blocks, evaluation -- happens in this graph.  Note the
  //! model's dag is never mutated: the mass coefficients and remainders are inserted FIRST and the
  //! derivative variables created HERE, so dy/dt exists only in the solver's own graph.
  FFGraph* _dag = nullptr;
  std::vector<FFVar> _pT, _pX, _pP, _pC, _pYP;    //!< the inserted variables, and dy/dt, all local
  struct StagePlan {
    std::vector<FFVar> Q;                    //!< this stage's _nq quadrature integrands, in _dag
    size_t             nFp = 0, offFp = 0;   //!< dF/dp: SFAD entries, and where they start in eval
    std::vector<int>   rFp, cFp;             //!< their (row, parameter) indices
    //!< dg/dy, dg/dy' and dg/dp for the quadrature integrands g, so a quadrature can be differentiated
    //!< by the same contraction the residual uses:  dq_j/dt = dg/dy . s_j + dg/dy' . s'_j + dg/dp_j
    size_t             nQy = 0, nQyp = 0, nQp = 0, offQy = 0, offQyp = 0, offQp = 0;
    std::vector<int>   rQy, cQy, rQyp, cQyp, rQp, cQp;
    std::vector<int>   rJy, cJy, rJyp, cJyp; //!< (row, state) of each dF/dy and dF/dy' entry.  ixJy/ixJyp
                                             //!< give where they land in the CSC value array, which the
                                             //!< Jacobian needs; the SENSITIVITY residual needs the
                                             //!< indices themselves, to contract against s and s'.
    std::vector<FFVar> F;                    //!< the residual rows of this stage
    std::vector<FFVar> eval;                 //!< [ F, dF/dy entries, dF/dy' entries ]
    std::vector<int>   ixJy, ixJyp;          //!< where each entry lands in the CSC value array
    std::vector<sunindextype> colptr, rowidx;
    size_t nJy = 0, nJyp = 0, nnz = 0, offJy = 0, offJyp = 0;
  };
  //! The consistent-IC system, assembled at setup from pieces that already exist: the stage-0 residual
  //! rows (the INTERIOR rows, which relate dx/dt to x) TOGETHER WITH the INITIAL rows, in the unknowns
  //! x and the DIFFERENTIAL dx/dt -- an algebraic state's derivative appears in no row.  Square when the
  //! INITIAL rows number the differential states, which is what initial_data() already audits.
  std::vector<StagePlan> _plan;              //!< one per stage (one entry when the rows do not change)
  std::vector<FFVar> _vYP;                   //!< unused; kept so older call sites still compile
  size_t _nJy = 0, _nJyp = 0;                //!< entry counts of the two blocks, for the active stage
  std::vector<int> _ixJy, _ixJyp;            //!< where each block entry lands in the CSC value array
  std::vector<sunindextype> _colptr, _rowidx;//!< the union pattern, CSC
  size_t _nnzjac = 0;                        //!< nonzeros of dF/dy + c_j dF/dy' for the CURRENT stage
  size_t _nnzjacmax = 0;                     //!< the worst stage's count, which the matrix is sized on
  SUNMatrix _sun_mat = nullptr;              //!< kept so a restart can reinitialize KLU on it
  SUNLinearSolver _sun_ls = nullptr;
  std::vector<FFVar>  _vInpLvl;              //!< the per-stage input LEVEL variables, in declaration order
  std::vector< std::pair<size_t,size_t> > _coefAt;  //!< (row, state) of each coefficient, in eval order
  //! @brief The function-evaluation context.  Its OWN variable list rather than _vVAR's, because a
  //! function term may reference a QUADRATURE (an INTEGRAL capture is carried as one) and the residual
  //! path has no quadrature variables in scope -- and because _vVAR is the hot path of every residual
  //! and Jacobian call, which this is not.
  std::vector<FFVar>  _pQv;                  //!< quadrature variables (_mQ), in _dag
  std::vector<FFVar>  _vFVAR;                //!< [ t, states, quadratures, parameters, constants ]
  std::vector<double> _dFVAR;                //!< their values
  std::vector< std::vector<FFVar> > _vFCT;   //!< function terms, one entry per stage time
  std::vector<double> _dF;                   //!< the accumulated function values, val_function() order
  std::vector< std::vector<FFVar> > _vFJac;  //!< df/d(x,q,p) per stage time, and its index pairs
  std::vector< std::vector<int> >   _rFJ, _cFJ;
  std::vector<double> _dFGrad;               //!< the accumulated gradient, _nf x _np row-major (internal)
  std::vector< std::vector<double> > _fp;    //!< the same as [ parameter ][ function ], the shape
                                             //!< ODESLV exposes and the FF operations index

  size_t              _nFp = 0, _offFp = 0;  //!< this stage's dF/dp block within _vEVAL
  std::vector<int>    _rFp, _cFp, _rJy, _cJy, _rJyp, _cJyp;
  size_t              _nQy = 0, _nQyp = 0, _nQp = 0, _offQy = 0, _offQyp = 0, _offQp = 0;
  std::vector<int>    _rQy, _cQy, _rQyp, _cQyp, _rQp, _cQp;

  std::vector<FFVar>  _vQUAD;                //!< this stage's quadrature integrands, in _dag
  std::vector<double> _dQUAD;                //!< their values
  bool                _fsa = false;          //!< this solve carries forward sensitivities
  bool                _asa = false;          //!< ... or runs the adjoint instead
  size_t              _asaF = (size_t)-1;    //!< which OUTPUT this backward problem is the adjoint of
  size_t              _asaQ = (size_t)-1;    //!< the quadrature it reads, if it is an integral output
  size_t              _asaWqCol = 0;         //!< the stage time the quadrature coefficient sits at
  bool                _asaHasX = false;      //!< the output has a TERMINAL part, so lambda(T) is seeded
  double              _asaWq = 1.;           //!< df/dq: the weight the running term enters with, which
                                             //!< is 1 for a bare integral and need not be in general
  std::vector<double> _asaHx;                //!< (df/dx)(T) -- the terminal seed
  std::vector< std::vector< std::vector<double> > > _sk;   //!< sensitivities at every stage time

  std::vector< std::vector< std::vector<double> > > _qpk;  //!< quadrature sensitivities per stage time
  std::vector< std::vector< std::vector<double> > > _lk;   //!< [stage][output][lambda], as ODESLV's
  std::vector< std::vector<double> > _xpk;   //!< state derivatives at every stage time (LEFT limit)
  std::vector< std::vector<double> > _xplus, _xpplus;   //!< the stop solutions (x+, dx+) at each stage start
  std::vector< std::vector<double> > _qk;    //!< quadratures at every stage time, val_quadrature() order

 public:
  //! @brief Recorded state trajectory of the last solve, for record().  Derivatives are NOT recorded:
  //! the format matches ODESLV's exactly, [ state, quadrature ] per row.
  std::vector< Results > results_solve;
  //! @brief Recorded sensitivity trajectories, one vector per parameter, as ODESLV indexes by direction.
  std::vector< std::vector< Results > > results_fsens;
  //! @brief Recorded adjoint trajectories, one vector per OUTPUT, in DESCENDING time -- the order the
  //! backward sweep produces them, which is what ODESLV stores too.
  std::vector< std::vector< Results > > results_asens;
 protected:
  std::vector< std::vector<double> > _fk;    //!< function values, per stage, val_function() order

  std::vector<FFVar>  _vEVAL;                //!< [ every coefficient, then every remainder ]
  std::vector<FFVar>  _vVAR;                 //!< [ t, states..., parameters... ]
  std::vector<double> _dVAR;                 //!< values of _vVAR
  std::vector<double> _dEVAL;                //!< values of _vEVAL
  std::vector<double> _dWRK;                 //!< eval work array
  std::vector<FFVar>  _vIC;                  //!< initial values, as expressions
  std::vector< std::vector<double> > _xk;    //!< state at the initial point and at each stage time
  std::vector<double> _dX, _dX0, _dXP0;      //!< state now, and (state, derivative) after IDACalcIC

  //! @brief The DAE residual: F = M(t,y,p) * yp + R(t,y,p).
  int IDRES__
    ( sunrealtype tt, N_Vector yy, N_Vector yp, N_Vector rr, void* user_data )
    override;

  //! @brief The quadrature integrands of the current stage.  IDAS hands this BOTH yy and yp, where the
  //! CVODES counterpart gets only y, so an integrand may depend on the state derivative -- strictly more
  //! than ODESLV can express.  Nothing uses that yet; the point is that the mass matrix costs nothing
  //! here, because the integrand is evaluated, never solved.
  int IDQUAD__
    ( sunrealtype tt, N_Vector yy, N_Vector yp, N_Vector rrQ, void* user_data )
    override
    { if( !_nq ) return -1;
      double const* y = N_VGetArrayPointer( yy ); double const* ypv = N_VGetArrayPointer( yp );
      _dVAR[0] = tt;
      for( size_t i = 0; i < _nx; ++i ){ _dVAR[1+i] = y[i]; _dVAR[1+_nx+i] = ypv[i]; }
      _dag->eval( _dWRK, _vQUAD.size(), _vQUAD.data(), _dQUAD.data(),
                  _vVAR.size(), _vVAR.data(), _dVAR.data() );
      double* q = N_VGetArrayPointer( rrQ );
      for( size_t i = 0; i < _nq; ++i ) q[i] = _dQUAD[i];
      return 0; }

  //! @brief v0 uses the IDAS difference-quotient Jacobian; this is never installed.
  //! @brief The analytic Jacobian  dF/dy + c_j dF/dy' .  With F = M(t,y,p) y' + R(t,y,p),
  //!     dF_k/dy_j  =  sum_i ( dM_ki/dy_j ) y'_i  +  dR_k/dy_j        and   dF_k/dy'_i = M_ki ,
  //! so every piece comes from the same symbolic rows the residual uses -- nothing is re-derived.
  int IDJAC__
    ( sunrealtype t, sunrealtype c_j, N_Vector yy, N_Vector yp, N_Vector rr, SUNMatrix Jac,
      void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3 )
    override;

  //! @brief Adjoint sensitivity is phase 6; refusing beats a silent wrong gradient.
  //! @brief The ADJOINT residual.  Read off IDAS's own idasAkzoNob_ASAi_dns example rather than
  //! derived:
  //!     resB = (dF/dy')^T . lambda'  -  (dF/dy)^T . lambda  +  (dg/dy)^T
  //! Its row 6 -- the ALGEBRAIC one -- carries no lambda' term at all, which is the general case
  //! working itself out: where dF/dy' has a zero row the derivative term drops and that row of the
  //! adjoint is algebraic in lambda, inheriting the index structure exactly as IDASetId does forward.
  //! Jyp and Jy are the SAME stage-plan blocks the state residual and the sensitivity residual use.
  int IDRESB__
    ( sunrealtype tt, N_Vector yy, N_Vector yp, N_Vector yyB, N_Vector ypB, N_Vector rrB, void* )
    override
    {
      _eval_rows( tt, N_VGetArrayPointer( yy ), N_VGetArrayPointer( yp ) );
      double const* lam = N_VGetArrayPointer( yyB );
      double const* lap = N_VGetArrayPointer( ypB );
      double* r = N_VGetArrayPointer( rrB );
      for( size_t i = 0; i < _nx; ++i ) r[i] = 0.;
      for( size_t e = 0; e < _nJyp; ++e ) r[ _cJyp[e] ] += _dEVAL[ _offJyp + e ] * lap[ _rJyp[e] ];
      for( size_t e = 0; e < _nJy;  ++e ) r[ _cJy[e]  ] -= _dEVAL[ _offJy  + e ] * lam[ _rJy[e]  ];
      // (dg/dy)^T for the objective's integrand: the quadrature this backward problem is the adjoint of
      if( _asaQ < _nq )
        for( size_t e = 0; e < _nQy; ++e )
          if( (size_t)_rQy[e] == _asaQ ) r[ _cQy[e] ] += _asaWq * _dEVAL[ _offQy + e ];
      return 0;
    }
  //! @brief The gradient quadrature of the backward problem:  dqB/dp_j = (dg/dp_j) - (dF/dp_j)^T lambda.
  //! dF/dp is the same SFAD block the forward sensitivity residual uses as its inhomogeneous term; this
  //! is that block transposed, which is the whole of what ASA costs once FSA exists.
  int IDQUADB__
    ( sunrealtype tt, N_Vector yy, N_Vector yp, N_Vector yyB, N_Vector, N_Vector rrQB, void* )
    override
    {
      _eval_rows( tt, N_VGetArrayPointer( yy ), N_VGetArrayPointer( yp ) );
      double const* lam = N_VGetArrayPointer( yyB );
      double* q = N_VGetArrayPointer( rrQB );
      for( size_t j = 0; j < _np; ++j ) q[j] = 0.;
      for( size_t e = 0; e < _nFp; ++e ) q[ _cFp[e] ] += _dEVAL[ _offFp + e ] * lam[ _rFp[e] ];
      if( _asaQ < _nq )
        for( size_t e = 0; e < _nQp; ++e )
          if( (size_t)_rQp[e] == _asaQ ) q[ _cQp[e] ] += _asaWq * _dEVAL[ _offQp + e ];
      return 0;
    }
  //! @brief Adjoint Jacobian; phase 6, refuses for now.
  int IDJACB__
    ( sunrealtype, sunrealtype, N_Vector, N_Vector, N_Vector, N_Vector, N_Vector, SUNMatrix, void*,
      N_Vector, N_Vector, N_Vector )
    override
    { return -1; }

  //! @brief Forward sensitivity is phase 5.
  //! @brief The FORWARD SENSITIVITY residual, one per direction:
  //!     F_s,j = dF/dy . s_j  +  dF/dy' . s'_j  +  dF/dp_j
  //! All three blocks are SFAD entries of the SAME stage plan the state residual uses, evaluated in one
  //! pass -- so the sensitivity equations cannot drift from the equations they differentiate, which is
  //! the failure mode this module has already paid for once on the initial conditions.
  int IDRESF__
    ( int Ns, sunrealtype tt, N_Vector yy, N_Vector yp, N_Vector, N_Vector* yyS,
      N_Vector* ypS, N_Vector* resvalS, void*, N_Vector, N_Vector, N_Vector )
    override
    {
      _eval_rows( tt, N_VGetArrayPointer( yy ), N_VGetArrayPointer( yp ) );
      size_t const nr = _nx;
      for( int j = 0; j < Ns; ++j ){
        double const* sj  = N_VGetArrayPointer( yyS[j] );
        double const* spj = N_VGetArrayPointer( ypS[j] );
        double* rj = N_VGetArrayPointer( resvalS[j] );
        for( size_t i = 0; i < nr; ++i ) rj[i] = 0.;
        for( size_t e = 0; e < _nJy;  ++e ) rj[ _rJy[e]  ] += _dEVAL[ _offJy  + e ] * sj [ _cJy[e]  ];
        for( size_t e = 0; e < _nJyp; ++e ) rj[ _rJyp[e] ] += _dEVAL[ _offJyp + e ] * spj[ _cJyp[e] ];
        for( size_t e = 0; e < _nFp;  ++e ) if( _cFp[e] == j ) rj[ _rFp[e] ] += _dEVAL[ _offFp + e ];
      }
      return 0;
    }
  //! @brief The QUADRATURE sensitivity integrands:  dq_j/dt = dg/dy . s_j + dg/dy' . s'_j + dg/dp_j.
  //! The same contraction as the sensitivity residual, over the same stage plan -- so a quadrature's
  //! sensitivity is differentiated from the integrand the quadrature itself uses, never from a copy.
  int IDQUADF__
    ( int Ns, sunrealtype tt, N_Vector yy, N_Vector yp, N_Vector* yyS, N_Vector* ypS, N_Vector,
      N_Vector* rhsQS, void*, N_Vector, N_Vector, N_Vector )
    override
    {
      if( !_nq ) return -1;
      double const* y = N_VGetArrayPointer( yy ); double const* ypv = N_VGetArrayPointer( yp );
      _dVAR[0] = tt;
      for( size_t i = 0; i < _nx; ++i ){ _dVAR[1+i] = y[i]; _dVAR[1+_nx+i] = ypv[i]; }
      _dag->eval( _dWRK, _vEVAL.size(), _vEVAL.data(), _dEVAL.data(),
                  _vVAR.size(), _vVAR.data(), _dVAR.data() );
      for( int j = 0; j < Ns; ++j ){
        double const* sj  = N_VGetArrayPointer( yyS[j] );
        double const* spj = N_VGetArrayPointer( ypS[j] );
        double* rj = N_VGetArrayPointer( rhsQS[j] );
        for( size_t i = 0; i < _nq; ++i ) rj[i] = 0.;
        for( size_t e = 0; e < _nQy;  ++e ) rj[ _rQy[e]  ] += _dEVAL[ _offQy  + e ] * sj [ _cQy[e]  ];
        for( size_t e = 0; e < _nQyp; ++e ) rj[ _rQyp[e] ] += _dEVAL[ _offQyp + e ] * spj[ _cQyp[e] ];
        for( size_t e = 0; e < _nQp;  ++e ) if( _cQp[e] == j ) rj[ _rQp[e] ] += _dEVAL[ _offQp + e ];
      }
      return 0;
    }

  //! @brief Evaluate the rows and their y-derivatives at (t, y, p) into _dEVAL.
  void _eval_rows
    ( double const t, double const* y, double const* yp );

  //! @brief Accumulate the function terms of stage time @p k.
  void _eval_fct
    ( size_t const k, double const t, double const* y, double const* q,
      std::vector< std::vector<double> > const& sfct = std::vector< std::vector<double> >(),
      std::vector< std::vector<double> > const& sqfct = std::vector< std::vector<double> >() );

  //! @brief Build the evaluation plan for stage @p s: the residual as ONE expression in (t, y, y', p), and
  //! its two sparse Jacobian blocks by SFAD, mirroring ODESLV_BASE's _pJAC = _dag->SFAD( nx, RHS, nx, X ).
  bool _build_rows
    ( size_t const s );

  //! @brief Build the solver-local dag and every stage's residual + Jacobian plan.
  bool _build_plans
    ();

  //! @brief Mass-matrix entries of stage @p s.
  size_t _ncoefOf
    ( size_t const s )
    const
    { size_t n = 0; for( auto const& row : _mMASS[s] ) n += row.size(); return n; }
};

//! @brief C trampolines for the ADJOINT residual and its gradient quadrature.
inline
int
MC_IDRESB__
( sunrealtype tt, N_Vector yy, N_Vector yp, N_Vector yyB, N_Vector ypB, N_Vector rrB, void* ud )
{
  return ( BASE_IDAS::PTR_BASE_IDAS->*PTR_IDRESB )( tt, yy, yp, yyB, ypB, rrB, ud );
}

inline
int
MC_IDQUADB__
( sunrealtype tt, N_Vector yy, N_Vector yp, N_Vector yyB, N_Vector ypB, N_Vector rrQB, void* ud )
{
  return ( BASE_IDAS::PTR_BASE_IDAS->*PTR_IDQUADB )( tt, yy, yp, yyB, ypB, rrQB, ud );
}

//! @brief C trampoline for the quadrature-sensitivity integrands, reaching the solver's IDQUADF__.
inline
int
MC_IDQUADF__
( int Ns, sunrealtype tt, N_Vector yy, N_Vector yp, N_Vector* yyS, N_Vector* ypS, N_Vector rrQ,
  N_Vector* rhsQS, void* user_data, N_Vector t1, N_Vector t2, N_Vector t3 )
{
  return ( BASE_IDAS::PTR_BASE_IDAS->*PTR_IDQUADF )( Ns, tt, yy, yp, yyS, ypS, rrQ, rhsQS,
                                                     user_data, t1, t2, t3 );
}

//! @brief C trampoline for the forward-sensitivity residual, reaching the solver's IDRESF__.
inline
int
MC_IDRESF__
( int Ns, sunrealtype tt, N_Vector yy, N_Vector yp, N_Vector resval, N_Vector* yyS,
  N_Vector* ypS, N_Vector* resvalS, void* user_data, N_Vector t1, N_Vector t2, N_Vector t3 )
{
  return ( BASE_IDAS::PTR_BASE_IDAS->*PTR_IDRESF )( Ns, tt, yy, yp, resval, yyS, ypS, resvalS,
                                                    user_data, t1, t2, t3 );
}

//! @brief C trampoline for the quadrature integrands, reaching the registered solver's IDQUAD__.
inline
int
MC_IDQUAD__
( sunrealtype tt, N_Vector yy, N_Vector yp, N_Vector rrQ, void* user_data )
{
  return ( BASE_IDAS::PTR_BASE_IDAS->*PTR_IDQUAD )( tt, yy, yp, rrQ, user_data );
}

//! @brief C trampoline for the residual: IDAS calls this, which reaches the registered solver's IDRES__.
//! Named as ODESLV_CVODES names its own (MC_CVRHS__, MC_CVJAC__).
inline
int
MC_IDRES__
( sunrealtype tt, N_Vector yy, N_Vector yp, N_Vector rr, void* user_data )
{
  return ( BASE_IDAS::PTR_BASE_IDAS->*PTR_IDRES )( tt, yy, yp, rr, user_data );
}

//! @brief C trampoline for the analytic Jacobian.
inline
int
MC_IDJAC__
( sunrealtype tt, sunrealtype c_j, N_Vector yy, N_Vector yp, N_Vector rr, SUNMatrix Jac, void* user_data,
  N_Vector tmp1, N_Vector tmp2, N_Vector tmp3 )
{
  return ( BASE_IDAS::PTR_BASE_IDAS->*PTR_IDJAC )( tt, c_j, yy, yp, rr, Jac, user_data, tmp1, tmp2, tmp3 );
}

inline
bool
DAESLV_IDAS::_build_plans
()
{
  delete _dag; _dag = new FFGraph;
  _pT.assign( 1, FFVar() );   _dag->insert( FFModel::_dag, 1,   _mT.data(), _pT.data() );
  _pX.assign( _nx, FFVar() ); _dag->insert( FFModel::_dag, _nx, _mX.data(), _pX.data() );
  _pP.assign( _np, FFVar() );
  if( _np ) _dag->insert( FFModel::_dag, _np, _mP.data(), _pP.data() );
  _pC.assign( _nc, FFVar() );
  if( _nc ) _dag->insert( FFModel::_dag, _nc, _mC.data(), _pC.data() );
  // the function context: quadrature variables, then the terms per stage time, exactly as ODESLV
  // assembles _vFCT -- the shared extraction already split every output into INTEGRAL and POINT parts
  _pQv.assign( _nq, FFVar() );
  if( _nq ) _dag->insert( FFModel::_dag, _nq, _mQ.data(), _pQv.data() );
  _vFVAR.clear();
  _vFVAR.push_back( _pT[0] );
  for( auto const& v : _pX  ) _vFVAR.push_back( v );
  for( auto const& v : _pQv ) _vFVAR.push_back( v );
  for( auto const& v : _pP  ) _vFVAR.push_back( v );
  for( auto const& v : _pC  ) _vFVAR.push_back( v );
  _dFVAR.assign( _vFVAR.size(), 0. );
  _dFGrad.assign( _nf * _np, 0. );
  _vFCT.clear();
  if( _nf ){
    if( _mFCT.size() == 1 ){                 // not distributed: the terms belong to the final time
      _vFCT.assign( _ns+1, std::vector<FFVar>() );
      _vFCT[_ns].assign( _nf, FFVar( 0. ) );
      _dag->insert( FFModel::_dag, _mFCT[0], _vFCT[_ns] );
    }
    else{
      _vFCT.assign( _mFCT.size(), std::vector<FFVar>( _nf, FFVar( 0. ) ) );
      size_t k = 0;
      for( auto const& mFCTk : _mFCT ) _dag->insert( FFModel::_dag, mFCTk, _vFCT[k++] );
    }
    // df/d(x, q, p) per stage time: the gradient of an output is the chain rule over what it reads --
    // the states, the quadratures it was built from, and the parameters it mentions directly.
    _vFJac.assign( _vFCT.size(), std::vector<FFVar>() );
    _rFJ.assign( _vFCT.size(), std::vector<int>() );
    _cFJ.assign( _vFCT.size(), std::vector<int>() );
    std::vector<FFVar> wrt;
    for( auto const& v : _pX  ) wrt.push_back( v );
    for( auto const& v : _pQv ) wrt.push_back( v );
    for( auto const& v : _pP  ) wrt.push_back( v );
    for( size_t k = 0; k < _vFCT.size(); ++k ){
      if( _vFCT[k].empty() || wrt.empty() ) continue;
      auto J = _dag->SFAD( _vFCT[k].size(), _vFCT[k].data(), wrt.size(), wrt.data() );
      size_t const n = std::get<0>(J);
      _rFJ[k].resize(n); _cFJ[k].resize(n);
      for( unsigned e = 0; e < n; ++e ){ _rFJ[k][e] = (int)std::get<1>(J)[e]; _cFJ[k][e] = (int)std::get<2>(J)[e]; }
      _vFJac[k].assign( std::get<3>(J), std::get<3>(J) + n );
      delete[] std::get<1>(J); delete[] std::get<2>(J); delete[] std::get<3>(J);
    }
  }

  _pYP.clear();
  for( size_t i = 0; i < _nx; ++i ) _pYP.push_back( _dag->add_var( "d" + _mX[i].name() + "/dt" ) );

  _vVAR.clear();
  _vVAR.push_back( _pT[0] );
  for( auto const& v : _pX  ) _vVAR.push_back( v );
  for( auto const& v : _pYP ) _vVAR.push_back( v );
  for( auto const& v : _pP  ) _vVAR.push_back( v );
  for( auto const& v : _pC  ) _vVAR.push_back( v );   // constants are variables of the eval, with fixed values
  _dVAR.assign( _vVAR.size(), 0. );

  _plan.assign( _mMASS.size(), StagePlan() );
  _nnzjacmax = 0;
  for( size_t sg = 0; sg < _mMASS.size(); ++sg ){
    StagePlan& P = _plan[sg];
    // insert this stage's coefficients and remainders, then build F HERE
    std::vector<FFVar> M( _ncoefOf( sg ) ), R( _mRES[sg].size() );
    std::vector<FFVar> Msrc;
    std::vector< std::pair<size_t,size_t> > at;
    for( size_t k = 0; k < _mMASS[sg].size(); ++k )
      for( auto const& [ ix, co ] : _mMASS[sg][k] ){ at.emplace_back( k, ix ); Msrc.push_back( co ); }
    if( !Msrc.empty() ) _dag->insert( FFModel::_dag, Msrc.size(), Msrc.data(), M.data() );
    _dag->insert( FFModel::_dag, _mRES[sg].size(), const_cast<FFVar*>( _mRES[sg].data() ), R.data() );
    // the quadrature integrands of this stage, from the SAME per-stage forms ODESLV uses: INTEGRAL
    // captures were turned into quadratures by the shared extraction, so there is nothing to re-derive
    if( _nq ){
      auto const& mq = _mQUAD[ sg < _mQUAD.size()? sg: _mQUAD.size()-1 ];
      P.Q.assign( _nq, FFVar() );
      _dag->insert( FFModel::_dag, _nq, const_cast<FFVar*>( mq.data() ), P.Q.data() );
    }
    P.F = R;
    for( size_t c = 0; c < at.size(); ++c ) P.F[ at[c].first ] += M[c] * _pYP[ at[c].second ];

    std::tuple<unsigned,unsigned*,unsigned*,FFVar*> Jy, Jyp, Fp;
    bool haveFp = false;
    try{
      Jy  = _dag->SFAD( P.F.size(), P.F.data(), _nx, _pX.data()  );
      Jyp = _dag->SFAD( P.F.size(), P.F.data(), _nx, _pYP.data() );
      // dF/dp, the inhomogeneous term of the sensitivity residual
      //     F_s,j = dF/dy . s_j + dF/dy' . s'_j + dF/dp_j
      // SFAD like the other two: a parameter touches few rows, and the dense form would build
      // rows x np expressions to find that out.
      if( _np ){ Fp = _dag->SFAD( P.F.size(), P.F.data(), _np, _pP.data() ); haveFp = true; }
    }
    catch( ... ){ _extractError = "SFAD failed on the DAE residual"; return false; }

    std::set< std::pair<size_t,size_t> > pat;                       // (col, row)
    for( unsigned e = 0; e < std::get<0>(Jy);  ++e ) pat.emplace( std::get<2>(Jy)[e],  std::get<1>(Jy)[e]  );
    for( unsigned e = 0; e < std::get<0>(Jyp); ++e ) pat.emplace( std::get<2>(Jyp)[e], std::get<1>(Jyp)[e] );
    P.nnz = pat.size();
    P.colptr.assign( _nx + 1, 0 );
    P.rowidx.clear();
    std::map< std::pair<size_t,size_t>, int > where;
    for( auto const& [ c, r ] : pat ){ where[{c,r}] = (int)P.rowidx.size(); P.rowidx.push_back( (sunindextype)r ); }
    for( auto const& [ c, r ] : pat ) ++P.colptr[ c + 1 ];
    for( size_t c = 0; c < _nx; ++c ) P.colptr[c+1] += P.colptr[c];
    P.nJy = std::get<0>(Jy);   P.ixJy.resize( P.nJy );  P.rJy.resize( P.nJy );  P.cJy.resize( P.nJy );
    for( unsigned e = 0; e < P.nJy;  ++e ){
      P.ixJy[e] = where[{ std::get<2>(Jy)[e], std::get<1>(Jy)[e] }];
      P.rJy[e]  = (int)std::get<1>(Jy)[e];  P.cJy[e] = (int)std::get<2>(Jy)[e];
    }
    P.nJyp = std::get<0>(Jyp); P.ixJyp.resize( P.nJyp ); P.rJyp.resize( P.nJyp ); P.cJyp.resize( P.nJyp );
    for( unsigned e = 0; e < P.nJyp; ++e ){
      P.ixJyp[e] = where[{ std::get<2>(Jyp)[e], std::get<1>(Jyp)[e] }];
      P.rJyp[e]  = (int)std::get<1>(Jyp)[e]; P.cJyp[e] = (int)std::get<2>(Jyp)[e];
    }

    P.eval = P.F;
    P.offJy  = P.eval.size();  P.eval.insert( P.eval.end(), std::get<3>(Jy),  std::get<3>(Jy)  + P.nJy  );
    P.offJyp = P.eval.size();  P.eval.insert( P.eval.end(), std::get<3>(Jyp), std::get<3>(Jyp) + P.nJyp );
    if( haveFp ){
      P.nFp = std::get<0>(Fp);
      P.rFp.resize( P.nFp ); P.cFp.resize( P.nFp );
      for( unsigned e = 0; e < P.nFp; ++e ){ P.rFp[e] = (int)std::get<1>(Fp)[e]; P.cFp[e] = (int)std::get<2>(Fp)[e]; }
      P.offFp = P.eval.size();
      P.eval.insert( P.eval.end(), std::get<3>(Fp), std::get<3>(Fp) + P.nFp );
      delete[] std::get<1>(Fp); delete[] std::get<2>(Fp); delete[] std::get<3>(Fp);
    }
    // the quadrature integrands differentiated the same three ways, appended to the same eval vector
    // so one pass values the residual, its Jacobian, and everything the sensitivities need
    if( _nq && !P.Q.empty() ){
      auto app = [&]( std::tuple<unsigned,unsigned*,unsigned*,FFVar*> J, size_t& n, size_t& off,
                      std::vector<int>& r, std::vector<int>& c ){
        n = std::get<0>(J);
        r.resize(n); c.resize(n);
        for( unsigned e = 0; e < n; ++e ){ r[e] = (int)std::get<1>(J)[e]; c[e] = (int)std::get<2>(J)[e]; }
        off = P.eval.size();
        P.eval.insert( P.eval.end(), std::get<3>(J), std::get<3>(J) + n );
        delete[] std::get<1>(J); delete[] std::get<2>(J); delete[] std::get<3>(J); };
      app( _dag->SFAD( P.Q.size(), P.Q.data(), _nx, _pX.data()  ), P.nQy,  P.offQy,  P.rQy,  P.cQy  );
      app( _dag->SFAD( P.Q.size(), P.Q.data(), _nx, _pYP.data() ), P.nQyp, P.offQyp, P.rQyp, P.cQyp );
      if( _np ) app( _dag->SFAD( P.Q.size(), P.Q.data(), _np, _pP.data() ), P.nQp, P.offQp, P.rQp, P.cQp );
    }
    for( auto* t : { &Jy, &Jyp } ){ delete[] std::get<1>(*t); delete[] std::get<2>(*t); delete[] std::get<3>(*t); }
    if( P.nnz > _nnzjacmax ) _nnzjacmax = P.nnz;
  }
  // The IC system.  Its rows live in the SAME local dag as everything else, so no cross-graph insert.
  return _build_rows( 0 );
}

inline
bool
DAESLV_IDAS::_build_rows
( size_t const s )
{
  size_t const ss = ( s < _plan.size() )? s: _plan.size() - 1;
  StagePlan const& P = _plan[ss];
  _vEVAL = P.eval;
  _vQUAD = P.Q;  _dQUAD.assign( _nq, 0. );
  _nFp = P.nFp;  _offFp = P.offFp;  _rFp = P.rFp;  _cFp = P.cFp;
  _rJy = P.rJy;  _cJy = P.cJy;  _rJyp = P.rJyp;  _cJyp = P.cJyp;
  _nQy = P.nQy; _offQy = P.offQy; _rQy = P.rQy; _cQy = P.cQy;
  _nQyp = P.nQyp; _offQyp = P.offQyp; _rQyp = P.rQyp; _cQyp = P.cQyp;
  _nQp = P.nQp; _offQp = P.offQp; _rQp = P.rQp; _cQp = P.cQp;
  _offJy = P.offJy;  _offJyp = P.offJyp;
  _nJy = P.nJy;      _nJyp = P.nJyp;
  _ixJy = P.ixJy;    _ixJyp = P.ixJyp;
  _colptr = P.colptr; _rowidx = P.rowidx;
  _nnzjac = P.nnz;
  _dEVAL.assign( _vEVAL.size(), 0. );
  return true;
}

inline
int
DAESLV_IDAS::IDJAC__
( sunrealtype t, sunrealtype c_j, N_Vector yy, N_Vector yp, N_Vector, SUNMatrix Jac, void*, N_Vector,
  N_Vector, N_Vector )
{
  _eval_rows( t, N_VGetArrayPointer( yy ), N_VGetArrayPointer( yp ) );
  std::vector<double> v( _nnzjac, 0. );
  for( size_t e = 0; e < _nJy;  ++e ) v[ _ixJy[e]  ] += _dEVAL[ _offJy  + e ];
  for( size_t e = 0; e < _nJyp; ++e ) v[ _ixJyp[e] ] += c_j * _dEVAL[ _offJyp + e ];
  if( SUNMatGetID( Jac ) == SUNMATRIX_SPARSE ){
    sunindextype* cp = SUNSparseMatrix_IndexPointers( Jac );
    sunindextype* ri = SUNSparseMatrix_IndexValues( Jac );
    sunrealtype*  d  = SUNSparseMatrix_Data( Jac );
    for( size_t c = 0; c <= _nx; ++c ) cp[c] = _colptr[c];
    for( size_t e = 0; e < _nnzjac; ++e ){ ri[e] = _rowidx[e]; d[e] = v[e]; }
  }
  else{
    SUNMatZero( Jac );
    for( size_t c = 0; c < _nx; ++c )
      for( sunindextype e = _colptr[c]; e < _colptr[c+1]; ++e )
        SM_ELEMENT_D( Jac, _rowidx[e], c ) = v[e];
  }
  ++stats.numJAC;
  return 0;
}

inline
bool
DAESLV_IDAS::setup
()
{
  if( !FFModel::setup() ){ _extractError = "FFModel::setup() failed"; return false; }
  if( !_extract_dae() ) return false;
  _nres = n_residual();
  // the sensitivity directions, from the registered controls -- as ODESLV resolves them after its own
  // extraction.  The FF operations check sensitivity_index().size() against n_control_dof().
  if( !_resolve_sensitivity_index() ) return false;

  if( !_build_plans() ) return false;
  _vInpLvl.clear();
  for( auto const& [ name, levels ] : _mInpLevels )
    for( auto const& lv : levels ) _vInpLvl.push_back( lv );
  _vIC.assign( _mIC[0].size(), FFVar() );
  _dag->insert( FFModel::_dag, _mIC[0].size(), const_cast<FFVar*>( _mIC[0].data() ), _vIC.data() );
  _dX.assign( _nx, 0. );
  _dX0.assign( _nx, 0. );
  _dXP0.assign( _nx, 0. );
  registration();
  return true;
}

//! Accumulate the function terms of stage time @p k.  The eval ADDS (the trailing unit scalar), because a
//! function can take contributions at several stage times -- a distributed output has one per stage.
inline
void
DAESLV_IDAS::_eval_fct
( size_t const k, double const t, double const* y, double const* q,
  std::vector< std::vector<double> > const& sfct, std::vector< std::vector<double> > const& sqfct )
{
  if( !_nf || k >= _vFCT.size() || _vFCT[k].empty() ) return;
  _dFVAR[0] = t;
  for( size_t i = 0; i < _nx; ++i ) _dFVAR[1+i] = y[i];
  for( size_t i = 0; i < _nq; ++i ) _dFVAR[1+_nx+i] = q? q[i]: 0.;
  for( size_t j = 0; j < _np; ++j ) _dFVAR[1+_nx+_nq+j] = _dVAR[ 1 + 2*_nx + j ];
  for( size_t j = 0; j < _nc; ++j ) _dFVAR[1+_nx+_nq+_np+j] = _dVAR[ 1 + 2*_nx + _np + j ];
  static double const one = 1.;
  _dag->eval( _nf, _vFCT[k].data(), _dF.data(), _vFVAR.size(), _vFVAR.data(), _dFVAR.data(), &one );

  // ... and the GRADIENT, by the same accumulation: df/dp_j = df/dx . s_j + df/dq . sq_j + df/dp_j
  if( _fsa && _np && k < _vFJac.size() && !_vFJac[k].empty() ){
    std::vector<double> dj( _vFJac[k].size(), 0. );
    _dag->eval( _vFJac[k].size(), _vFJac[k].data(), dj.data(),
                _vFVAR.size(), _vFVAR.data(), _dFVAR.data() );
    for( size_t e = 0; e < dj.size(); ++e ){
      size_t const r = (size_t)_rFJ[k][e], c = (size_t)_cFJ[k][e];
      if( r >= _nf ) continue;
      for( size_t j = 0; j < _np; ++j ){
        double w = 0.;
        if( c < _nx )                 w = ( j < sfct.size() && c < sfct[j].size() )? sfct[j][c]: 0.;
        else if( c < _nx + _nq )      w = ( j < sqfct.size() && (c-_nx) < sqfct[j].size() )? sqfct[j][c-_nx]: 0.;
        else if( c - _nx - _nq == j ) w = 1.;
        _dFGrad[ r*_np + j ] += dj[e] * w;
      }
    }
  }
}

inline
void
DAESLV_IDAS::_eval_rows
( double const t, double const* y, double const* yp )
{
  _dVAR[0] = t;
  for( size_t i = 0; i < _nx; ++i ) _dVAR[1+i] = y[i];
  for( size_t i = 0; i < _nx; ++i ) _dVAR[1+_nx+i] = yp[i];
  _dag->eval( _dWRK, _vEVAL.size(), _vEVAL.data(), _dEVAL.data(),
              _vVAR.size(), _vVAR.data(), _dVAR.data() );
}

inline
int
DAESLV_IDAS::IDRES__
( sunrealtype tt, N_Vector yy, N_Vector yp, N_Vector rr, void* )
{
  sunrealtype const* y  = N_VGetArrayPointer( yy );
  sunrealtype const* dy = N_VGetArrayPointer( yp );
  sunrealtype*       r  = N_VGetArrayPointer( rr );
  _eval_rows( tt, y, dy );
  for( size_t k = 0; k < _nres; ++k ) r[k] = _dEVAL[k];
  ++stats.numRES;
  return 0;
}

inline
DAESLV_IDAS::STATUS
DAESLV_IDAS::solve
( std::vector<double> const& p, std::vector<double> const& c, std::ostream& os )
{
  if( _vEVAL.empty() ){ _extractError = "solve() before a successful setup()"; return FATAL; }
  if( p.size() != _np ){
    _extractError = std::to_string( p.size() ) + " parameter value(s) for " + std::to_string( _np )
                  + " parameter(s)";
    return FATAL;
  }
  if( c.size() != _nc ){
    _extractError = std::to_string( c.size() ) + " constant value(s) for " + std::to_string( _nc );
    return FATAL;
  }
  for( size_t j = 0; j < _np; ++j ) _dVAR[ 1 + 2*_nx + j ] = p[j];
  for( size_t j = 0; j < _nc; ++j ) _dVAR[ 1 + 2*_nx + _np + j ] = c[j];
  _dIcCst.assign( c.begin(), c.begin() + _nc );   // the IC system is evaluated in its own graph
  _xk.clear();  _xk.reserve( _ns + 1 );

  // Initial values: the INITIAL rows, which may depend on the parameters AND on the states -- a
  // materialised hidden constraint becomes ic[i] = -R/L with R carrying other states.  Seed from the
  // declared references first, so the result does not depend on what the previous call left in _dVAR.
  // NOTE this is ONE PASS: when the rows are implicit in each other it is a starting point for
  // IDACalcIC, not a solution of them.  See the plan's IC-solve item.
  std::vector<double> y0( _nx, 0. );
  // Stage 0's rows are built by setup(), but the stage loop leaves them at the LAST stage, so a second
  // solve() on the same object started stage 0 with the previous run's final input levels -- measured at
  // 3.4e-02 on an 8-element model, reproducibly.  Rebuilt HERE, before the IDAS memory and its linear
  // solver are sized from the rows, not after.  Invisible to a single-stage model, which is why the
  // idempotence check in DAESLV_test1 passed throughout: with one stage the rows never change.
  if( _mMASS.size() > 1 && !_build_rows( 0 ) ){ _extractError = "_build_rows(0)"; return FAILURE; }

  _dVAR[0] = _dT[0];
  { // _classification_reference() is the canonical reference point and is VIRTUAL -- OCFESLV overrides it
    // to sample the collocated field -- so take the seed from it rather than re-deriving one here.
    std::vector<double> sref, iref, cref, dref;
    _classification_reference( sref, iref, cref, dref );
    for( size_t i = 0; i < _nx; ++i ){ _dVAR[1+i] = i < sref.size()? sref[i]: 0.; _dVAR[1+_nx+i] = 0.; } }
  // The DIVIDED form is only usable where each INITIAL row is affine in the state it pins AND free of
  // derivatives: _split_in leaves anything else in the remainder, so ic[x] can carry z and dx/dt and
  // cannot be evaluated at all.  With a solve available the divided form is not needed -- the reference
  // seed goes straight to solve_ic, which takes the rows as they are.
  if( options.INITIAL.SOLVER != BASE_IDAS::Options::t_Initial::KINSOL )
    _dag->eval( _dWRK, _vIC.size(), _vIC.data(), y0.data(),
                _vVAR.size(), _vVAR.data(), _dVAR.data() );
  else
    for( size_t i = 0; i < _nx; ++i ) y0[i] = _dVAR[1+i];      // the reference seed, already in _dVAR

  // With INITIAL.SOLVER = KINSOL the one-pass evaluation above is only the SEED: solve_ic() then resolves
  // the INITIAL rows together with the interior rows at t0, which is what rows implicit in each other
  // need -- one evaluation cannot resolve them, and _mIC's divided form cannot even represent them.
  if( options.INITIAL.SOLVER == BASE_IDAS::Options::t_Initial::KINSOL ){
    for( size_t i = 0; i < _nx; ++i ) _dVAR[1+i] = y0[i];
    std::vector<double> yv, ypv;
    BASE_KINSOL::Report const r = solve_ic( yv, ypv );
    if( r.status == BASE_KINSOL::FAILED ){ _extractError = "initial-condition solve: " + r.message; return FAILURE; }
    if( r.status == BASE_KINSOL::RANK_DEFICIENT )
      os << "DAESLV_IDAS ** " << r.message << std::endl;   // a warning: the point is usable, just inexact
    y0 = yv;
    _dXP0 = ypv;
  }

  _init_stats( stats );
  SUNContext ctx = nullptr;
  if( SUNContext_Create( SUN_COMM_NULL, &ctx ) < 0 ){ _extractError = "SUNContext_Create failed"; return FATAL; }
  N_Vector yy = N_VNew_Serial( (sunindextype)_nx, ctx );
  N_Vector yp = N_VNew_Serial( (sunindextype)_nx, ctx );
  N_Vector id = N_VNew_Serial( (sunindextype)_nx, ctx );
  for( size_t i = 0; i < _nx; ++i ){
    N_VGetArrayPointer( yy )[i] = y0[i];
    N_VGetArrayPointer( yp )[i] = 0.;
    N_VGetArrayPointer( id )[i] = _mID[i];
  }

  bool ok = true;
  std::string why;
  N_Vector Nq = nullptr;                       // the quadrature vector, null when the model has none
  N_Vector* yS = nullptr; N_Vector* ypS = nullptr;   // the sensitivity vectors, null without FSA
  N_Vector* yQS = nullptr;                           // and the quadrature sensitivities
  void* mem = IDACreate( ctx );
  // mirrors ODESLV_CVODES's switch, minus DIAG, which IDAS does not have
#if !defined( CRONOS__WITH_KLU )
  if( options.SOLVE.LINSOL == BASE_IDAS::Options::SPARSE ){
    std::cerr << "DAESLV_IDAS ** LINSOL = SPARSE requires KLU, which this build does not include"
              << "; choose DENSE or DENSEDQ, or use a build with KLU" << std::endl;
    N_VDestroy( yy ); N_VDestroy( yp ); N_VDestroy( id ); if( Nq ) N_VDestroy( Nq );
  if( yS ) N_VDestroyVectorArray( yS, (int)_np );
  if( yQS ) N_VDestroyVectorArray( yQS, (int)_np );
  if( ypS ) N_VDestroyVectorArray( ypS, (int)_np );
    if( yS ) N_VDestroyVectorArray( yS, (int)_np );
    if( yQS ) N_VDestroyVectorArray( yQS, (int)_np );
    if( ypS ) N_VDestroyVectorArray( ypS, (int)_np ); if( Nq ) N_VDestroy( Nq );
    IDAFree( &mem ); SUNContext_Free( &ctx );
    _extractError = "LINSOL = SPARSE without KLU";
    return FATAL;
  }
#endif
  SUNMatrix A = nullptr;
  SUNLinearSolver LS = nullptr;
#if defined( CRONOS__WITH_KLU )
  if( options.SOLVE.LINSOL == BASE_IDAS::Options::SPARSE ){
    A  = SUNSparseMatrix( (sunindextype)_nx, (sunindextype)_nx, (sunindextype)_nnzjacmax, CSC_MAT, ctx );
    LS = SUNLinSol_KLU( yy, A, ctx );
  }
  else
#endif
  {
    A  = SUNDenseMatrix( (sunindextype)_nx, (sunindextype)_nx, ctx );
    LS = SUNLinSol_Dense( yy, A, ctx );
  }
  _sun_mat = A; _sun_ls = LS;
  reregistration();
  do{
    if( IDAInit( mem, MC_IDRES__, _dT[0], yy, yp ) != IDA_SUCCESS ){ ok = false; why = "IDAInit"; break; }
    IDASStolerances( mem, options.SOLVE.RTOL, options.SOLVE.ATOL );

    // Forward sensitivities.  yS0/ypS0 come from ic_sens(), which solved
    //     ( dF/dz ) . ( dz/dp ) = - dF/dp
    // against the Jacobian the state IC solve had just factorised -- so the sensitivity initial data
    // is the derivative OF the initial data, not of some other object that resembles it.
    if( _fsa && _np ){
      std::vector< std::vector<double> > dxdp, ddxdp;
      if( !ic_sens( dxdp, ddxdp ) ){ ok = false; why = "ic_sens for the sensitivity initial data"; break; }
      yS  = N_VCloneVectorArray( (int)_np, yy );
      ypS = N_VCloneVectorArray( (int)_np, yp );
      if( !yS || !ypS ){ ok = false; why = "N_VCloneVectorArray"; break; }
      for( size_t j = 0; j < _np; ++j )
        for( size_t i = 0; i < _nx; ++i ){
          N_VGetArrayPointer( yS[j]  )[i] = j < dxdp.size()?  dxdp[j][i]:  0.;
          N_VGetArrayPointer( ypS[j] )[i] = j < ddxdp.size()? ddxdp[j][i]: 0.;
        }
      int const ism = ( options.SOLVE.FSACORR == BASE_IDAS::Options::SIMULTANEOUS )? IDA_SIMULTANEOUS: IDA_STAGGERED;
      if( IDASensInit( mem, (int)_np, ism, MC_IDRESF__, yS, ypS ) != IDA_SUCCESS ){
        ok = false; why = "IDASensInit"; break; }
      if( IDASensEEtolerances( mem ) != IDA_SUCCESS ){ ok = false; why = "IDASensEEtolerances"; break; }
      IDASetSensErrCon( mem, options.SOLVE.FSAERR );
    }

    // Quadratures, mirroring ODESLV's CVodeQuad* block.  Ordering is deliberate and was measured, not
    // assumed: IDAQuadInit AFTER IDAInit and the tolerances, and the ReInit at a restart AFTER
    // IDAReInit, because both allocate against the memory the preceding call sized.
    if( _nq ){
      Nq = N_VNew_Serial( (sunindextype)_nq, ctx );
      for( size_t i = 0; i < _nq; ++i ) N_VGetArrayPointer( Nq )[i] = 0.;
      if( IDAQuadInit( mem, MC_IDQUAD__, Nq ) != IDA_SUCCESS ){ ok = false; why = "IDAQuadInit"; break; }
      IDASetQuadErrCon( mem, options.SOLVE.QERR );
      if( options.SOLVE.QERR
       && IDAQuadSStolerances( mem, options.SOLVE.RTOL, options.SOLVE.ATOL ) != IDA_SUCCESS ){
        ok = false; why = "IDAQuadSStolerances"; break; }
      // Quadrature SENSITIVITIES, after IDASensInit (which allocated the sensitivity vectors this
      // reads) and after IDAQuadInit (which sized the quadrature this differentiates).  The same
      // ordering rule as everywhere else here: init against the memory the preceding call sized.
      if( _fsa && _np ){
        yQS = N_VCloneVectorArray( (int)_np, Nq );
        if( !yQS ){ ok = false; why = "N_VCloneVectorArray for the quadrature sensitivities"; break; }
        for( size_t j = 0; j < _np; ++j ) N_VConst( 0., yQS[j] );
        if( IDAQuadSensInit( mem, MC_IDQUADF__, yQS ) != IDA_SUCCESS ){
          ok = false; why = "IDAQuadSensInit"; break; }
        IDASetQuadSensErrCon( mem, options.SOLVE.QERRS );
        if( options.SOLVE.QERRS && IDAQuadSensEEtolerances( mem ) != IDA_SUCCESS ){
          ok = false; why = "IDAQuadSensEEtolerances"; break; }
      }
    }
    if( IDASetLinearSolver( mem, LS, A ) != IDA_SUCCESS ){ ok = false; why = "IDASetLinearSolver"; break; }
    IDASetId( mem, id );
    if( options.SOLVE.LINSOL != BASE_IDAS::Options::DENSEDQ )
      IDASetJacFn( mem, MC_IDJAC__ );         // analytic; DENSEDQ leaves IDAS its difference quotient
    IDASetSuppressAlg( mem, options.SOLVE.SUPPRESSALG );
    IDASetMaxOrd( mem, options.SOLVE.MAXORD );
    if( options.SOLVE.NMAX ) IDASetMaxNumSteps( mem, (long)options.SOLVE.NMAX );
    if( options.SOLVE.HMAX > 0. ) IDASetMaxStep( mem, options.SOLVE.HMAX );
    if( options.SOLVE.H0 > 0. ) IDASetInitStep( mem, options.SOLVE.H0 );
    IDASetMaxNumStepsIC( mem, options.INITIAL.NMAXIC );
    IDASetMaxNumJacsIC( mem, options.INITIAL.MAXJACIC );
    IDASetMaxNumItersIC( mem, options.INITIAL.MAXNLIIC );
    IDASetLineSearchOffIC( mem, options.INITIAL.LSOFFIC );
    if( options.INITIAL.STEPTOLIC > 0. ) IDASetStepToleranceIC( mem, options.INITIAL.STEPTOLIC );

    if( options.INITIAL.CALCIC != BASE_IDAS::Options::NONE ){
      double t1 = options.INITIAL.TOUT1;
      if( t1 == 0. ) t1 = _dT.size() > 1? _dT[1]: _dT[0] + 1.;
      int const icopt = ( options.INITIAL.CALCIC == BASE_IDAS::Options::YA_YDP )? IDA_YA_YDP_INIT: IDA_Y_INIT;
      if( IDACalcIC( mem, icopt, t1 ) != IDA_SUCCESS ){ ok = false; why = "IDACalcIC"; break; }
      IDAGetConsistentIC( mem, yy, yp );
      ++stats.numIC;
    }
    for( size_t i = 0; i < _nx; ++i ){
      _dX0[i]  = N_VGetArrayPointer( yy )[i];
      _dXP0[i] = N_VGetArrayPointer( yp )[i];
    }
    if( options.OUTPUT.RESRECORD ){
    }
    _xk.push_back( _dX0 );                       // the initial point, as ODESLV_CVODES records it
    results_solve.clear();  results_fsens.clear();  results_asens.clear();
    bool const rec = options.OUTPUT.RESRECORD;     // opt-in, as in ODESLV: recording costs a copy
    if( rec && ( _fsa || _asa ) ){ results_fsens.resize( _np );  results_asens.resize( _nf ); }
    _xpk.clear();  _xpk.push_back( _dXP0 );
    _xplus.assign( _ns, std::vector<double>() );  _xpplus.assign( _ns, std::vector<double>() );
    _qk.clear();  _qk.push_back( std::vector<double>( _nq, 0. ) );   // quadratures start at zero
    _sk.clear();  _qpk.clear();
    if( _fsa && yS ){
      std::vector< std::vector<double> > s0( _np, std::vector<double>( _nx, 0. ) );
      for( size_t j = 0; j < _np; ++j )
        for( size_t i = 0; i < _nx; ++i ) s0[j][i] = N_VGetArrayPointer( yS[j] )[i];
      _sk.push_back( s0 );
      if( rec ) for( size_t j = 0; j < _np; ++j )
        results_fsens[j].push_back( Results( _dT[0], _nx, s0[j].data(),
                                             _nq, yQS? N_VGetArrayPointer( yQS[j] ): nullptr ) );
    }
    if( rec ) results_solve.push_back( Results( _dT[0], _nx, _dX0.data(), _nq, _qk.back().data() ) );
    _dF.assign( _nf, 0. );                                           // the function accumulator
    _dFGrad.assign( _nf * _np, 0. );
    { std::vector< std::vector<double> > s0, sq0;
      if( _fsa && yS ){
        s0.assign( _np, std::vector<double>( _nx, 0. ) );
        for( size_t j = 0; j < _np; ++j )
          for( size_t i = 0; i < _nx; ++i ) s0[j][i] = N_VGetArrayPointer( yS[j] )[i];
        sq0.assign( _np, std::vector<double>( _nq, 0. ) );            // quadratures start at zero
      }
      _eval_fct( 0, _dT[0], _dX0.data(), _qk.back().data(), s0, sq0 ); }   // the t0 contribution

    // The rows differ stage by stage exactly when an input is time-varying -- _mMASS then holds _ns
    // entries rather than one.  At such a boundary the BDF history describes a problem that no longer
    // exists (ODESLV's _stage_restart rule), AND, unlike an ODE, the state itself is no longer
    // consistent: the DIFFERENTIAL components carry over, but the ALGEBRAIC ones are whatever the new
    // constraints say and generally JUMP.  So the restart is IDAReInit followed by IDACalcIC in
    // IDA_YA_YDP_INIT, which holds the differential components and re-solves the algebraic ones and the
    // derivatives -- the split the id vector already carries.  Merely re-initialising would resume from
    // an inconsistent point, and with SUPPRESSALG on that first step can look plausible rather than fail.
    // A stop exists wherever the rows change OR a transition was declared: without the second clause a
    // jump on a model with no time-varying input was never applied, the loop never restarting there.
    // ASA is armed BEFORE any step: IDAAdjInit'd memory must be driven by IDASolveF, never by a plain
    // IDASolve.  Multi-stage backward solving mirrors ODESLV's per-stage scheme -- re-integrate each
    // stage forward on its own checkpoints, then adjoin across it -- and is NOT done yet, so a model
    // with more than one stage is REFUSED rather than silently given the last stage's adjoint.
    bool const restarts = ( _mMASS.size() > 1 ) || !var_transition().empty();

    // ASA is armed BEFORE any step: IDAAdjInit'd memory must be driven by IDASolveF, never by a plain
    // IDASolve.  The restriction is RESTARTS, not stages: a stage boundary that neither changes the
    // rows nor carries a transition is not reinitialised, so the forward sweep runs continuously and
    // IDASolveF checkpoints the whole horizon -- the backward solve crosses it with nothing to do.
    // Where the forward DOES reinitialise, the checkpoints do not span the boundary and the backward
    // sweep must go stage by stage, re-integrating each on its own checkpoints as ODESLV does, with
    // the stop transpose applied between them.  That is not written, so it is REFUSED.
    if( _asa ){
      // Which kind of output is this?  Its own Jacobian says: a column in the QUADRATURE block makes it
      // an integral (a running g_y term, lambda(T) = 0); a column in the STATE block at the final stage
      // time makes it a point value (no running term, lambda(T) seeded from df/dx there).  Detected
      // rather than declared, so a caller cannot mislabel it.
      _asaQ = (size_t)-1;  _asaHx.assign( _nx, 0. );
      if( _asaF >= _nf ){ ok = false; why = "solve_asens: no such output"; break; }
      { bool anyQ = false, anyX = false, anyT = false;
        size_t const kT = _vFCT.empty()? 0: _vFCT.size()-1;
        for( size_t k = 0; k < _rFJ.size(); ++k )
          for( size_t e = 0; e < _rFJ[k].size(); ++e ){
            if( (size_t)_rFJ[k][e] != _asaF ) continue;
            size_t const col = (size_t)_cFJ[k][e];
            if( col >= _nx && col < _nx + _nq ){ anyQ = true; _asaQ = col - _nx; _asaWqCol = k; }
            else if( col < _nx ){ anyX = true;  if( k == kT ) anyT = true; }
          }
        _asaHasX = anyT;      // the TERMINAL seed; interior terms arrive as impulses in the sweep
        (void)anyX;
        if( !anyQ && !anyX ){ ok = false; why = "solve_asens: that output depends on neither the states"
                                               " at any stage time nor a quadrature"; break; }
      }
      if( IDAAdjInit( mem, (int)options.SOLVE.ASACHKPT, IDA_HERMITE ) != IDA_SUCCESS ){
        ok = false; why = "IDAAdjInit"; break; }
    }
    for( size_t sg = 0; sg < _ns && ok; ++sg ){
      double const ta = _dT[sg], tb = _dT[sg+1];
      if( sg && restarts ){
        if( !_build_rows( sg ) ){ ok = false; why = "_build_rows"; break; }
        if( IDAReInit( mem, ta, yy, yp ) != IDA_SUCCESS ){ ok = false; why = "IDAReInit"; break; }
#if defined( CRONOS__WITH_KLU )
        // the stage's rows changed, so the Jacobian's nonzero count or structure may have too: force a new
        // SYMBOLIC and numeric factorization, exactly as ODESLV_CVODES does at its own restarts
        if( options.SOLVE.LINSOL == BASE_IDAS::Options::SPARSE
         && SUNLinSol_KLUReInit( _sun_ls, _sun_mat, (sunindextype)_nnzjacmax, 2 ) ){
          ok = false; why = "SUNLinSol_KLUReInit"; break;
        }
#endif
        if( options.INITIAL.SOLVER == BASE_IDAS::Options::t_Initial::KINSOL ){
          // A stage boundary is a consistent initialisation: the stage's rows together with CONTINUITY on
          // the differential states, the left state as data.  Same solve as t0 (plan §3.13), and unlike
          // IDACalcIC it reports WHY it failed -- a degenerate configuration is named, not just flagged.
          std::vector<double> const xL( N_VGetArrayPointer( yy ), N_VGetArrayPointer( yy ) + _nx );
          std::vector<double> xp, dxp, par( _np, 0. );
          for( size_t j = 0; j < _np; ++j ) par[j] = _dVAR[ 1 + 2*_nx + j ];
          std::vector<double> const dxL( N_VGetArrayPointer( yp ), N_VGetArrayPointer( yp ) + _nx );
          if( _fsa && yS ){        // (s^-, s'^-): the jump law reads x^-, so its sensitivity reads s^-
            std::vector< std::vector<double> > sL( _np, std::vector<double>( _nx, 0. ) ),
                                               sdL( _np, std::vector<double>( _nx, 0. ) );
            for( size_t j = 0; j < _np; ++j )
              for( size_t i = 0; i < _nx; ++i ){
                sL[j][i]  = N_VGetArrayPointer( yS[j]  )[i];
                sdL[j][i] = N_VGetArrayPointer( ypS[j] )[i];
              }
            set_ic_left_sens( sL, sdL );
          }
          BASE_KINSOL::Report const R = restart_ic( ta, xL, xp, dxp, options.INITIAL.KIN, par, dxL );
          _xplus[sg] = xp;  _xpplus[sg] = dxp;
          if( R.status == BASE_KINSOL::FAILED ){ ok = false; why = "restart_ic: " + R.message; break; }
          if( R.status == BASE_KINSOL::RANK_DEFICIENT )
            os << "DAESLV_IDAS ** at the stage boundary: " << R.message << std::endl;
          for( size_t i = 0; i < _nx; ++i ){
            N_VGetArrayPointer( yy )[i] = xp[i];
            N_VGetArrayPointer( yp )[i] = dxp[i];
          }
          if( IDAReInit( mem, ta, yy, yp ) != IDA_SUCCESS ){ ok = false; why = "IDAReInit after restart_ic"; break; }
          ++stats.numIC;
          // The sensitivities are reinitialised from the SAME stop system the state just used, by the
          // same linear solve: ic_sens() linearises about whatever restart_ic() converged to, so the
          // closing rows that applied there -- continuity, or a declared jump -- apply to s and s'
          // too, with no separate derivation of the jump law.  Carrying BOTH s^- and s'^- is what
          // plan 3.13 requires, and the stop system has both because dx/dt is one of its unknowns.
          if( _fsa && yS ){
            std::vector< std::vector<double> > dxdp, ddxdp;
            if( !ic_sens( dxdp, ddxdp ) ){ ok = false; why = "ic_sens at the stop"; break; }
            for( size_t j = 0; j < _np; ++j )
              for( size_t i = 0; i < _nx; ++i ){
                N_VGetArrayPointer( yS[j]  )[i] = j < dxdp.size()?  dxdp[j][i]:  0.;
                N_VGetArrayPointer( ypS[j] )[i] = j < ddxdp.size()? ddxdp[j][i]: 0.;
              }
            int const ism2 = ( options.SOLVE.FSACORR == BASE_IDAS::Options::SIMULTANEOUS )? IDA_SIMULTANEOUS: IDA_STAGGERED;
            if( IDASensReInit( mem, ism2, yS, ypS ) != IDA_SUCCESS ){ ok = false; why = "IDASensReInit"; break; }
          }
        }
        else if( options.INITIAL.CALCIC != BASE_IDAS::Options::NONE ){
          int const icopt = ( options.INITIAL.CALCIC == BASE_IDAS::Options::YA_YDP )? IDA_YA_YDP_INIT: IDA_Y_INIT;
          if( IDACalcIC( mem, icopt, tb ) != IDA_SUCCESS ){ ok = false; why = "IDACalcIC at a stage boundary"; break; }
          IDAGetConsistentIC( mem, yy, yp );
          ++stats.numIC;
        }
      }
      // The quadrature is reset to ZERO at every stage, as ODESLV's _IC_D_QUAD does, so a quadrature is a
      // PER-STAGE increment and the accumulation across stages is the FUNCTION's job: extraction puts one
      // term per stage time in _mFCT, and _eval_fct adds them.  Carrying the accumulator across instead
      // double-counts every earlier stage -- measured as exactly Q(tau) too much on a two-stage model.
      // Gated on the stage index ALONE, not on whether the rows change: a stage boundary exists wherever
      // a POINT capture put one, and the quadrature must restart there too.  Gating it on `restarts` left
      // the accumulator running on a model whose stages come from an interior tau rather than an input.
      if( ok && _nq && sg ){
        for( size_t i = 0; i < _nq; ++i ) N_VGetArrayPointer( Nq )[i] = 0.;
        if( IDAQuadReInit( mem, Nq ) != IDA_SUCCESS ){ ok = false; why = "IDAQuadReInit"; break; }
        // the quadrature sensitivity is the derivative of a PER-STAGE increment, so it restarts at
        // zero with it; the accumulation across stages is val_function_gradient()'s, as the value's
        // accumulation is val_function()'s
        if( _fsa && yQS ){
          for( size_t j = 0; j < _np; ++j ) N_VConst( 0., yQS[j] );
          if( IDAQuadSensReInit( mem, yQS ) != IDA_SUCCESS ){ ok = false; why = "IDAQuadSensReInit"; break; }
        }
      }

      // Stop AT the boundary rather than stepping past it and interpolating back, as ODESLV does with
      // CVodeSetStopTime.  Without it IDAS leaves its internal time beyond tb, so zeroing the quadrature
      // at the next stage zeroes it at that later instant and the interval [tb, tn] is lost from the
      // stage increment -- measured as 4.6e-03 on a two-stage model, i.e. as if stage 2 began at 0.4069
      // rather than 0.4.  The state was unaffected, being interpolated back to tb, which is why only the
      // quadrature showed it.
      // The STAGE-START row, at the same t as the previous stage's last: _record writes a blank line
      // when it sees a repeated time, which separates the blocks -- and across a declared jump the two
      // rows ARE the discontinuity, x(tau^-) then x(tau^+).  ODESLV pushes the same pair.
      if( ok && rec && sg ){
        results_solve.push_back( Results( ta, _nx, N_VGetArrayPointer( yy ),
                                          _nq, _nq? N_VGetArrayPointer( Nq ): nullptr ) );
        if( _fsa && yS ){
          sunrealtype ts0 = ta;
          if( IDAGetSens( mem, &ts0, yS ) == IDA_SUCCESS )
            for( size_t j = 0; j < _np; ++j )
              results_fsens[j].push_back( Results( ta, _nx, N_VGetArrayPointer( yS[j] ),
                                                   _nq, yQS? N_VGetArrayPointer( yQS[j] ): nullptr ) );
        }
      }

      // RESRECORD intermediate points: the solver is told to STOP at each of them, as ODESLV does with
      // CVodeSetStopTime, rather than being run to the boundary and interpolated back afterwards -- so a
      // recorded row is a point the integrator actually landed on.  The last stop is forced exactly onto
      // tb, as ODESLV forces its own, so no accumulation of TSTEP can drift off the boundary.
      sunrealtype tret = ta;
      if( ok && tret < tb ){
        unsigned const NSTEP = ( rec && options.OUTPUT.RESRECORD )? options.OUTPUT.RESRECORD: 1;
        double const TSTEP = ( tb - ta ) / (double)NSTEP;
        double TSTOP = ta + TSTEP;
        for( unsigned st = 0; st < NSTEP && ok; ++st, TSTOP += TSTEP ){
          if( st+1 == NSTEP ) TSTOP = tb;
          if( IDASetStopTime( mem, TSTOP ) != IDA_SUCCESS ){ ok = false; why = "IDASetStopTime"; break; }
          // IDA_TSTOP_RETURN, not IDA_SUCCESS, is what a solve that lands ON the stop time returns.
          int ncheck = 0;
          int const fl = _asa? IDASolveF( mem, TSTOP, &tret, yy, yp, IDA_NORMAL, &ncheck )
                             : IDASolve ( mem, TSTOP, &tret, yy, yp, IDA_NORMAL );
          if( fl != IDA_SUCCESS && fl != IDA_TSTOP_RETURN ){ ok = false; why = "IDASolve to the stage boundary"; break; }
          if( !rec || st+1 == NSTEP ) continue;        // the boundary row is pushed below, beside _xk
          if( _nq && IDAGetQuadDky( mem, TSTOP, 0, Nq ) != IDA_SUCCESS ){ ok = false; why = "IDAGetQuadDky"; break; }
          results_solve.push_back( Results( TSTOP, _nx, N_VGetArrayPointer( yy ),
                                            _nq, _nq? N_VGetArrayPointer( Nq ): nullptr ) );
          if( _fsa && yS ){
            sunrealtype ts2 = TSTOP;
            if( IDAGetSens( mem, &ts2, yS ) != IDA_SUCCESS ){ ok = false; why = "IDAGetSens"; break; }
            sunrealtype tq3 = TSTOP;
            if( yQS && IDAGetQuadSens( mem, &tq3, yQS ) != IDA_SUCCESS ){ ok = false; why = "IDAGetQuadSens"; break; }
            for( size_t j = 0; j < _np; ++j )
              results_fsens[j].push_back( Results( TSTOP, _nx, N_VGetArrayPointer( yS[j] ),
                                                   _nq, yQS? N_VGetArrayPointer( yQS[j] ): nullptr ) );
          }
        }
      }
      if( ok ){
        _xk.push_back( std::vector<double>( N_VGetArrayPointer( yy ), N_VGetArrayPointer( yy ) + _nx ) );
        _xpk.push_back( std::vector<double>( N_VGetArrayPointer( yp ), N_VGetArrayPointer( yp ) + _nx ) );
        // the recorded row, beside the vector it mirrors.  Hooking this to a neighbouring line
        // instead is how the first attempt recorded only t0.
        if( rec ) results_solve.push_back( Results( tret, _nx, N_VGetArrayPointer( yy ),
                                                   _nq, _nq? _qk.back().data(): nullptr ) );
      }
      if( ok && _nq ){
        // Dky rather than IDAGetQuad, as ODESLV uses CVodeGetQuadDky: with the stop time set the two
        // agree, but Dky asks for the value AT tret instead of at whatever the last internal step was,
        // which is what the caller means.
        if( IDAGetQuadDky( mem, tret, 0, Nq ) != IDA_SUCCESS ){ ok = false; why = "IDAGetQuadDky"; break; }
        if( _fsa && yQS ){
          sunrealtype tq2 = tret;
          if( IDAGetQuadSens( mem, &tq2, yQS ) != IDA_SUCCESS ){ ok = false; why = "IDAGetQuadSens"; break; }
        }
        _qk.push_back( std::vector<double>( N_VGetArrayPointer( Nq ), N_VGetArrayPointer( Nq ) + _nq ) );
      }
      else if( ok ) _qk.push_back( std::vector<double>( _nq, 0. ) );
      if( ok && _fsa && yS ){
        if( IDAGetSens( mem, &tret, yS ) != IDA_SUCCESS ){ ok = false; why = "IDAGetSens"; break; }
        std::vector< std::vector<double> > sk( _np, std::vector<double>( _nx, 0. ) );
        for( size_t j = 0; j < _np; ++j )
          for( size_t i = 0; i < _nx; ++i ) sk[j][i] = N_VGetArrayPointer( yS[j] )[i];
        _sk.push_back( sk );
        { std::vector< std::vector<double> > qk( _np, std::vector<double>( _nq, 0. ) );
          if( yQS ) for( size_t j = 0; j < _np; ++j )
            for( size_t m = 0; m < _nq; ++m ) qk[j][m] = N_VGetArrayPointer( yQS[j] )[m];
          _qpk.push_back( qk ); }
        // the BOUNDARY row; the RESRECORD loop above records only the intermediate ones
        if( rec ) for( size_t j = 0; j < _np; ++j )
          results_fsens[j].push_back( Results( tret, _nx, sk[j].data(),
                                               _nq, yQS? N_VGetArrayPointer( yQS[j] ): nullptr ) );
        // (the stale note below is what the loop does NOT cover)
      }
      if( ok ){
        std::vector< std::vector<double> > sf, sqf;
        if( _fsa && yS ){
          sf.assign( _np, std::vector<double>( _nx, 0. ) );
          for( size_t j = 0; j < _np; ++j )
            for( size_t i = 0; i < _nx; ++i ) sf[j][i] = N_VGetArrayPointer( yS[j] )[i];
          sqf.assign( _np, std::vector<double>( _nq, 0. ) );
          if( yQS ) for( size_t j = 0; j < _np; ++j )
            for( size_t i = 0; i < _nq; ++i ) sqf[j][i] = N_VGetArrayPointer( yQS[j] )[i];
        }
        _eval_fct( sg+1, tret, N_VGetArrayPointer( yy ), _qk.back().data(), sf, sqf );
      }
    }
    if( !ok ) break;
    for( size_t i = 0; i < _nx; ++i ) _dX[i] = N_VGetArrayPointer( yy )[i];
    long n = 0;
    IDAGetNumSteps( mem, &n );    stats.numSteps = (unsigned)n;
    IDAGetNumJacEvals( mem, &n ); stats.numJAC   = (unsigned)n;
  } while( false );

  // ---- the BACKWARD sweep.  lambda(T) = 0 for an objective that is a pure quadrature: there is no
  // terminal term to seed it with.  One backward problem per adjoined output; a single stage here, so
  // no re-integration is needed -- the forward sweep's checkpoints still cover it.
  if( ok && _asa ){
    int indexB = -1;
    N_Vector yB = N_VNew_Serial( (sunindextype)_nx, ctx );
    N_Vector ypB = N_VNew_Serial( (sunindextype)_nx, ctx );
    N_Vector qB = N_VNew_Serial( (sunindextype)std::max<size_t>( _np, 1 ), ctx );
    SUNMatrix AB = nullptr; SUNLinearSolver LSB = nullptr;
    N_VConst( 0., yB ); N_VConst( 0., ypB ); N_VConst( 0., qB );
    double const tEnd = _dT.back(), tBeg = _dT.front();
    // df/dq, the weight the running term enters with.  1 for a bare integral output, and read from the
    // same Jacobian that classified the output rather than assumed.  Evaluated at the final stage time:
    // a coefficient that is itself state-dependent AND sits at an earlier stage time would want
    // re-evaluating there, which no output in the corpus does.
    _asaWq = 0.;
    if( _asaQ < _nq && _asaWqCol < _vFJac.size() && !_vFJac[_asaWqCol].empty() ){
      std::vector<double> dj( _vFJac[_asaWqCol].size(), 0. );
      _dag->eval( _vFJac[_asaWqCol].size(), _vFJac[_asaWqCol].data(), dj.data(),
                  _vFVAR.size(), _vFVAR.data(), _dFVAR.data() );
      for( size_t e = 0; e < dj.size(); ++e )
        if( (size_t)_rFJ[_asaWqCol][e] == _asaF && (size_t)_cFJ[_asaWqCol][e] == _nx + _asaQ )
          _asaWq += dj[e];
    }
    do{
      if( IDACreateB( mem, &indexB ) != IDA_SUCCESS ){ ok = false; why = "IDACreateB"; break; }
      if( IDAInitB( mem, indexB, MC_IDRESB__, tEnd, yB, ypB ) != IDA_SUCCESS ){
        ok = false; why = "IDAInitB"; break; }
      if( IDASStolerancesB( mem, indexB, options.SOLVE.RTOLB, options.SOLVE.ATOLB ) != IDA_SUCCESS ){
        ok = false; why = "IDASStolerancesB"; break; }
      AB = SUNDenseMatrix( (sunindextype)_nx, (sunindextype)_nx, ctx );
      LSB = SUNLinSol_Dense( yB, AB, ctx );
      if( !AB || !LSB || IDASetLinearSolverB( mem, indexB, LSB, AB ) != IDA_SUCCESS ){
        ok = false; why = "IDASetLinearSolverB"; break; }
      if( _np && IDAQuadInitB( mem, indexB, MC_IDQUADB__, qB ) != IDA_SUCCESS ){
        ok = false; why = "IDAQuadInitB"; break; }
      // lambda(T) = 0 is the right VALUE but not a consistent STATE: resB = Jyp^T lambda' - Jy^T lambda
      // + g_y, so lambda'(T) must satisfy Jyp^T lambda'(T) = -g_y there.  Seeding it at zero made
      // IDASolveB fail outright rather than drift, which is the good failure.  IDACalcICB solves it,
      // and needs the backward id vector: adjoint row i is DIFFERENTIAL exactly where column i of Jyp
      // has an entry, the transpose of the forward partition.
      // A POINT output has no running integrand; its whole influence enters through the TERMINAL
      // condition  Jyp^T lambda(T) = ( df/dx )( T ).  _dFVAR still holds the final stage time from the
      // forward sweep's last _eval_fct, so the function Jacobian is evaluated where it belongs.
      if( _asaHasX && !_vFJac.empty() ){
        size_t const kT = _vFJac.size()-1;
        if( !_vFJac[kT].empty() ){
          std::vector<double> dj( _vFJac[kT].size(), 0. );
          _dag->eval( _vFJac[kT].size(), _vFJac[kT].data(), dj.data(),
                      _vFVAR.size(), _vFVAR.data(), _dFVAR.data() );
          for( size_t e = 0; e < dj.size(); ++e )
            if( (size_t)_rFJ[kT][e] == _asaF && (size_t)_cFJ[kT][e] < _nx )
              _asaHx[ _cFJ[kT][e] ] += dj[e];
        }
        // Jyp^T lambda = h_x, densely, at the final time.  Singular where Jyp has zero rows, which is
        // exactly where lambda is algebraic -- IDACalcICB resolves those below.
        _eval_rows( tEnd, N_VGetArrayPointer( yy ), N_VGetArrayPointer( yp ) );
        std::vector<double> Mt( _nx*_nx, 0. );
        for( size_t e = 0; e < _nJyp; ++e ) Mt[ (size_t)_cJyp[e] + _nx*(size_t)_rJyp[e] ] += _dEVAL[ _offJyp + e ];
        SUNMatrix AT = SUNDenseMatrix( (sunindextype)_nx, (sunindextype)_nx, ctx );
        N_Vector bT = N_VNew_Serial( (sunindextype)_nx, ctx );
        SUNLinearSolver LST = SUNLinSol_Dense( bT, AT, ctx );
        bool okT = ( AT && bT && LST );
        if( okT ){
          for( size_t c = 0; c < _nx; ++c )
            for( size_t r = 0; r < _nx; ++r ) SM_ELEMENT_D( AT, r, c ) = Mt[ r + _nx*c ];
          for( size_t r = 0; r < _nx; ++r ) N_VGetArrayPointer( bT )[r] = _asaHx[r];
          okT = ( SUNLinSolSetup( LST, AT ) == 0 ) && ( SUNLinSolSolve( LST, AT, bT, bT, 0. ) == 0 );
        }
        if( okT ) for( size_t i = 0; i < _nx; ++i ) N_VGetArrayPointer( yB )[i] = N_VGetArrayPointer( bT )[i];
        if( LST ) SUNLinSolFree( LST );
        if( AT ) SUNMatDestroy( AT );
        if( bT ) N_VDestroy( bT );
        if( !okT ){ ok = false; why = "the terminal adjoint condition Jyp^T lambda(T) = df/dx is singular"; break; }
        if( IDAReInitB( mem, indexB, tEnd, yB, ypB ) != IDA_SUCCESS ){ ok = false; why = "IDAReInitB with the terminal seed"; break; }
      }

      { N_Vector idB = N_VNew_Serial( (sunindextype)_nx, ctx );
        N_VConst( 0., idB );
        for( size_t e = 0; e < _nJyp; ++e ) N_VGetArrayPointer( idB )[ _cJyp[e] ] = 1.;
        int const rc = IDASetIdB( mem, indexB, idB );
        N_VDestroy( idB );
        if( rc != IDA_SUCCESS ){ ok = false; why = "IDASetIdB"; break; }
        if( IDACalcICB( mem, indexB, tBeg, yy, yp ) != IDA_SUCCESS ){
          ok = false; why = "IDACalcICB (the backward initial condition)"; break; }
        if( IDAGetConsistentICB( mem, indexB, yB, ypB ) != IDA_SUCCESS ){
          ok = false; why = "IDAGetConsistentICB"; break; }
      }
      if( _np && options.SOLVE.QERRB ){
        IDASetQuadErrConB( mem, indexB, SUNTRUE );
        if( IDAQuadSStolerancesB( mem, indexB, options.SOLVE.RTOLB, options.SOLVE.ATOLB ) != IDA_SUCCESS ){
          ok = false; why = "IDAQuadSStolerancesB"; break; }
      }
      // Backward STAGE BY STAGE, so an output term at an INTERIOR stage time can be picked up.  A
      // measurement at tau gives lambda an impulse there:  Jyp^T lambda(tau^-) = Jyp^T lambda(tau^+)
      // + (df/dx)(tau).  ODESLV expresses the same thing as one Hamiltonian  f_k(x) + lambda^T.map
      // differentiated in one go; here the map term is the stop system rather than a divided map, so
      // the two parts are taken separately -- this is the measurement part, which is what parameter
      // estimation against time-series data needs and which needs no restart machinery.
      sunrealtype tb2 = tBeg;
      if( _lk.size() != _ns+1 ) _lk.assign( _ns+1, std::vector<std::vector<double>>( _nf ) );
      if( _asaF < _nf )                       // lambda(T), at the TERMINAL -- as ODESLV fills _lk[_ns]
        _lk[_ns][_asaF].assign( N_VGetArrayPointer( yB ), N_VGetArrayPointer( yB ) + _nx );
      for( size_t kk = _ns; kk > 0 && ok; --kk ){
        if( kk < _ns ){
          // the impulse at this stage time, from the output's own Jacobian evaluated THERE
          std::vector<double> hx( _nx, 0. );
          if( kk < _vFJac.size() && !_vFJac[kk].empty() ){
            _dFVAR[0] = _dT[kk];
            for( size_t i = 0; i < _nx && kk < _xk.size(); ++i ) _dFVAR[1+i] = _xk[kk][i];
            for( size_t m = 0; m < _nq && kk < _qk.size(); ++m ) _dFVAR[1+_nx+m] = _qk[kk][m];
            std::vector<double> dj( _vFJac[kk].size(), 0. );
            _dag->eval( _vFJac[kk].size(), _vFJac[kk].data(), dj.data(),
                        _vFVAR.size(), _vFVAR.data(), _dFVAR.data() );
            for( size_t e = 0; e < dj.size(); ++e )
              if( (size_t)_rFJ[kk][e] == _asaF && (size_t)_cFJ[kk][e] < _nx ) hx[ _cFJ[kk][e] ] += dj[e];
          }
          bool any = false;
          for( double v : hx ) if( v != 0. ){ any = true; break; }
          if( any ){
            _eval_rows( _dT[kk], _xk[kk].data(), _dXP0.data() );
            std::vector<double> Mt( _nx*_nx, 0. );
            for( size_t e = 0; e < _nJyp; ++e )
              Mt[ (size_t)_cJyp[e] + _nx*(size_t)_rJyp[e] ] += _dEVAL[ _offJyp + e ];
            SUNMatrix AK = SUNDenseMatrix( (sunindextype)_nx, (sunindextype)_nx, ctx );
            N_Vector bK = N_VNew_Serial( (sunindextype)_nx, ctx );
            SUNLinearSolver LSK = SUNLinSol_Dense( bK, AK, ctx );
            bool okK = ( AK && bK && LSK );
            if( okK ){
              for( size_t c = 0; c < _nx; ++c )
                for( size_t r = 0; r < _nx; ++r ) SM_ELEMENT_D( AK, r, c ) = Mt[ r + _nx*c ];
              for( size_t r = 0; r < _nx; ++r ) N_VGetArrayPointer( bK )[r] = hx[r];
              okK = ( SUNLinSolSetup( LSK, AK ) == 0 ) && ( SUNLinSolSolve( LSK, AK, bK, bK, 0. ) == 0 );
            }
            if( okK )
              for( size_t i = 0; i < _nx; ++i )
                N_VGetArrayPointer( yB )[i] += N_VGetArrayPointer( bK )[i];
            if( LSK ) SUNLinSolFree( LSK );
            if( AK ) SUNMatDestroy( AK );
            if( bK ) N_VDestroy( bK );
            if( !okK ){ ok = false; why = "the interior adjoint impulse Jyp^T dlambda = df/dx is singular"; break; }
          }
          // A RESTARTING boundary: the forward solved a stop system here, so the co-state must cross it
          // by that system's TRANSPOSE -- v^- from v^+ -- and the stop contributes its own -dF/dp^T mu
          // to the gradient.  Then lambda^- solves Jyp^-^T lambda^- = v^-.  ODESLV does the same through
          // _CC_D_SEN, the transpose of its divided map; here the map is a system, so it is a solve.
          bool const restarts_bwd = ( _mMASS.size() > 1 ) || !var_transition().empty();
          bool const stopped_here = restarts_bwd && kk < _xplus.size() && !_xplus[kk].empty();
          if( stopped_here ){
            // v^+ = Jyp^+^T lambda^+ at the start of stage kk, plus whatever impulse was just added
            _eval_rows( _dT[kk], _xplus[kk].data(), _xpplus[kk].data() );
            std::vector<double> vplus( _nx, 0. );
            double const* l = N_VGetArrayPointer( yB );
            for( size_t e = 0; e < _nJyp; ++e ) vplus[ _cJyp[e] ] += _dEVAL[ _offJyp + e ] * l[ _rJyp[e] ];
            std::vector<double> vminus, dvminus, dGdp( _np, 0. );
            std::vector<double> parv( _np, 0. );
            for( size_t j = 0; j < _np; ++j ) parv[j] = _dVAR[ 1 + 2*_nx + j ];
            if( !ic_adjoint( _dT[kk], _xk[kk], _xpk[kk], _xplus[kk], _xpplus[kk], vplus, parv,
                             vminus, dvminus, dGdp ) ){ ok = false; why = "ic_adjoint at a stop"; break; }
            for( size_t j = 0; j < _np; ++j ) N_VGetArrayPointer( qB )[j] += dGdp[j];
            // lambda^- from v^-: the rows of stage kk-1 at its END, where Jyp^- is evaluated
            if( !_build_rows( kk-1 ) ){ ok = false; why = "_build_rows for the adjoint stop"; break; }
            _eval_rows( _dT[kk], _xk[kk].data(), _xpk[kk].data() );
            std::vector<double> Mt( _nx*_nx, 0. );
            for( size_t e = 0; e < _nJyp; ++e ) Mt[ (size_t)_cJyp[e] + _nx*(size_t)_rJyp[e] ] += _dEVAL[ _offJyp + e ];
            SUNMatrix AS = SUNDenseMatrix( (sunindextype)_nx, (sunindextype)_nx, ctx );
            N_Vector bS = N_VNew_Serial( (sunindextype)_nx, ctx );
            SUNLinearSolver LSS = SUNLinSol_Dense( bS, AS, ctx );
            bool okS = ( AS && bS && LSS );
            if( okS ){
              for( size_t c = 0; c < _nx; ++c )
                for( size_t r = 0; r < _nx; ++r ) SM_ELEMENT_D( AS, r, c ) = Mt[ r + _nx*c ];
              for( size_t r = 0; r < _nx; ++r ) N_VGetArrayPointer( bS )[r] = vminus[r];
              okS = ( SUNLinSolSetup( LSS, AS ) == 0 ) && ( SUNLinSolSolve( LSS, AS, bS, bS, 0. ) == 0 );
            }
            if( okS ) for( size_t i = 0; i < _nx; ++i ) N_VGetArrayPointer( yB )[i] = N_VGetArrayPointer( bS )[i];
            if( LSS ) SUNLinSolFree( LSS );
            if( AS ) SUNMatDestroy( AS );
            if( bS ) N_VDestroy( bS );
            if( !okS ){ ok = false; why = "Jyp^T lambda^- = v^- at a stop is singular"; break; }
            // The forward sweep's checkpoints do not span a reinitialisation.  Re-integrate THIS stage
            // forward on its own, as ODESLV does, so the backward solve across it has checkpoints.
            if( IDAAdjReInit( mem ) != IDA_SUCCESS ){ ok = false; why = "IDAAdjReInit"; break; }
            for( size_t i = 0; i < _nx; ++i ){
              N_VGetArrayPointer( yy )[i] = kk-1 < _xplus.size() && !_xplus[kk-1].empty()? _xplus[kk-1][i]: _xk[kk-1][i];
              N_VGetArrayPointer( yp )[i] = kk-1 < _xpplus.size() && !_xpplus[kk-1].empty()? _xpplus[kk-1][i]: _xpk[kk-1][i];
            }
            if( IDAReInit( mem, _dT[kk-1], yy, yp ) != IDA_SUCCESS ){ ok = false; why = "IDAReInit for the re-integration"; break; }
            // the rows changed with _build_rows: a sparse factorisation must be redone on the new
            // pattern, exactly as the forward restart does after ITS _build_rows
#if defined( CRONOS__WITH_KLU )
            if( options.SOLVE.LINSOL == BASE_IDAS::Options::SPARSE
             && SUNLinSol_KLUReInit( _sun_ls, _sun_mat, (sunindextype)_nnzjacmax, 2 ) ){
              ok = false; why = "SUNLinSol_KLUReInit (re-integration)"; break; }
#endif
            { int nck = 0; sunrealtype tr = _dT[kk-1];
              if( IDASetStopTime( mem, _dT[kk] ) != IDA_SUCCESS ){ ok = false; why = "IDASetStopTime (re-integration)"; break; }
              int const fl = IDASolveF( mem, _dT[kk], &tr, yy, yp, IDA_NORMAL, &nck );
              if( fl != IDA_SUCCESS && fl != IDA_TSTOP_RETURN ){ ok = false; why = "IDASolveF (re-integration)"; break; } }
            any = true;      // lambda changed: the backward problem must be reinitialised below
          }
          // Reinitialise ONLY where lambda actually jumped.  Doing it at every stage time restarts the
          // backward BDF history for nothing and cost three orders: 5e-08 where the single-sweep
          // version had 8e-11.  qB is passed back so the accumulated gradient carries across.
          if( any ){
            if( IDAReInitB( mem, indexB, _dT[kk], yB, ypB ) != IDA_SUCCESS ){
              ok = false; why = "IDAReInitB at an interior stage time"; break; }
            if( _np && IDAQuadReInitB( mem, indexB, qB ) != IDA_SUCCESS ){
              ok = false; why = "IDAQuadReInitB"; break; }
          }
        }
        // The backward sweep sub-samples by RESRECORD exactly as the forward one does, stopping at
        // each recorded time rather than interpolating, and pushes the STAGE-END row first so the
        // repeated t separates the blocks -- the same pair that shows a discontinuity forward.
        bool const brec = options.OUTPUT.RESRECORD && _asaF < results_asens.size();
        if( brec ){
          if( _np && IDAGetQuadB( mem, indexB, &tb2, qB ) != IDA_SUCCESS ){ ok = false; why = "IDAGetQuadB"; break; }
          results_asens[_asaF].push_back( Results( _dT[kk], _nx, N_VGetArrayPointer( yB ),
                                                   _np, _np? N_VGetArrayPointer( qB ): nullptr ) );
        }
        { unsigned const NB = options.OUTPUT.RESRECORD? options.OUTPUT.RESRECORD: 1;
          double const TB = ( _dT[kk] - _dT[kk-1] ) / (double)NB;
          double tgt = _dT[kk] - TB;
          for( unsigned st = 0; st < NB && ok; ++st, tgt -= TB ){
            if( st+1 == NB ) tgt = _dT[kk-1];           // land exactly on the stage start
            if( IDASolveB( mem, tgt, IDA_NORMAL ) < 0 ){ ok = false; why = "IDASolveB"; break; }
            if( IDAGetB( mem, indexB, &tb2, yB, ypB ) != IDA_SUCCESS ){ ok = false; why = "IDAGetB"; break; }
            if( _np && IDAGetQuadB( mem, indexB, &tb2, qB ) != IDA_SUCCESS ){ ok = false; why = "IDAGetQuadB"; break; }
            // the t0 row waits: qB there is only the RUNNING part of the gradient, and the
            // initial-condition term is added below -- recording it early is what left ODESLV's
            // d f/da = 0.5596 against a column of zeros here, for a parameter that acts only
            // through x(0)
            if( brec )
              results_asens[_asaF].push_back( Results( tgt, _nx, N_VGetArrayPointer( yB ),
                                                       _np, _np? N_VGetArrayPointer( qB ): nullptr ) );
            if( st+1 == NB && _asaF < _nf )          // the stage boundary, as ODESLV records it
              _lk[kk-1][_asaF].assign( N_VGetArrayPointer( yB ), N_VGetArrayPointer( yB ) + _nx );
          }
        }
        if( !ok ) break;
      }
      if( !ok ) break;
      // the gradient of the adjoined output lands in the same place the forward one does, so the two
      // are directly comparable -- which is the gate
      // The quadrature covers only what the parameters do INSIDE the horizon.  A parameter that enters
      // through the initial data -- x(0) = a, say -- has no dF/dp anywhere and would read as having no
      // influence at all.  Its contribution is the adjoint of what ic_sens() supplies forward:
      //     dG/dp += ( Jyp^T lambda )( t0 ) . dx0/dp
      // Jyp^T rather than lambda alone because the mass matrix maps the adjoint onto the state, and
      // dx0/dp from the SAME system the state's initial condition came from.
      // a point output may also mention the parameters DIRECTLY
      std::vector<double> dp( _np, 0. );
      if( !_vFJac.empty() ){
        size_t const kT = _vFJac.size()-1;
        if( !_vFJac[kT].empty() ){
          std::vector<double> dj( _vFJac[kT].size(), 0. );
          _dag->eval( _vFJac[kT].size(), _vFJac[kT].data(), dj.data(),
                      _vFVAR.size(), _vFVAR.data(), _dFVAR.data() );
          for( size_t e = 0; e < dj.size(); ++e ){
            size_t const col = (size_t)_cFJ[kT][e];
            if( (size_t)_rFJ[kT][e] == _asaF && col >= _nx + _nq && col - _nx - _nq < _np )
              dp[ col - _nx - _nq ] += dj[e];
          }
        }
      }
      std::vector<double> ic( _np, 0. );   // one per PARAMETER -- sized by _nx it aliased across them
      { std::vector< std::vector<double> > dxdp, ddxdp;
        // the START system, explicitly: ic_sens() linearises whatever was solved LAST, which after a
        // restart is the last STOP -- and dx+/dp at a jump is not dx0/dp.  Measured: on x(1+) = a x(1-)
        // the term came out as lambda(t0).x(1-) = 0.1075 on top of the right 0.1653.
        _icLast = &_icStart;
        if( ic_sens( dxdp, ddxdp ) && !dxdp.empty() ){
          _eval_rows( tBeg, _dX0.data(), _dXP0.data() );
          std::vector<double> Ml( _nx, 0. );
          double const* l0 = N_VGetArrayPointer( yB );
          for( size_t e = 0; e < _nJyp; ++e ) Ml[ _cJyp[e] ] += _dEVAL[ _offJyp + e ] * l0[ _rJyp[e] ];
          for( size_t j = 0; j < _np && j < dxdp.size(); ++j ){
            double acc = 0.;
            for( size_t i = 0; i < _nx && i < dxdp[j].size(); ++i ) acc += Ml[i] * dxdp[j][i];
            ic[j] = acc;
          }
        } }
      _dFGrad.assign( _nf * _np, 0. );
      // qB becomes the WHOLE gradient, not just its running part, so the recorded lambda_q at t0 IS
      // what val_function_gradient() reports rather than a piece of it
      for( size_t j = 0; j < _np; ++j ) N_VGetArrayPointer( qB )[j] += ic[j] + dp[j];
      // A SECOND row at t0, carrying the whole gradient.  It repeats t0, so _record gives it a block
      // of its own at the head of the file -- ODESLV writes the same pair, the summary row then the
      // trajectory's own t0 row with the running part alone.
      if( options.OUTPUT.RESRECORD && _asaF < results_asens.size() )
        results_asens[_asaF].push_back( Results( tBeg, _nx, N_VGetArrayPointer( yB ),
                                                 _np, _np? N_VGetArrayPointer( qB ): nullptr ) );
      for( size_t j = 0; j < _np; ++j )
        for( size_t r = 0; r < _nf; ++r )
          if( blk_fct( r ).first == r )
            _dFGrad[ r*_np + j ] = N_VGetArrayPointer( qB )[j];
    } while( false );
    if( LSB ) SUNLinSolFree( LSB );
    if( AB ) SUNMatDestroy( AB );
    N_VDestroy( yB ); N_VDestroy( ypB ); N_VDestroy( qB );
  }

  // The forward linear solver is freed HERE, after the backward sweep: the re-integration of a stage
  // drives the FORWARD problem again (IDAReInit + IDASolveF) and needs it.  Freeing it before the
  // backward block -- as it was -- left every backward call working against a dead linear solver whose
  // memory the backward vectors then reused: an int 1 (their length) over its vtable, then 0.0, then
  // x(1^-) itself.  Found by tracing _sun_mat->ops, not by reading.
  SUNLinSolFree( LS );
  SUNMatDestroy( A );
  // val_function_gradient() is [ parameter ][ function ], ODESLV's shape -- the FF operations read
  // gradient()[ op input i ][ output k ] and both solvers must agree on which index is which
  _fp.assign( _np, std::vector<double>( _nf, 0. ) );
  if( _dFGrad.size() == _nf*_np )
    for( size_t r = 0; r < _nf; ++r ) for( size_t j = 0; j < _np; ++j ) _fp[j][r] = _dFGrad[ r*_np + j ];

  // file the counters under the kind of solve this was
  if     ( _asa ) stats_asens = stats;
  else if( _fsa ) stats_fsens = stats;
  else            stats_solve = stats;

  N_VDestroy( yy ); N_VDestroy( yp ); N_VDestroy( id );
  IDAFree( &mem );
  SUNContext_Free( &ctx );
  _final_stats( stats );
  if( !ok ){ _extractError = why + " failed"; os << "DAESLV_IDAS ** " << _extractError << std::endl; }
  return ok? NORMAL: FAILURE;
}

} // end namespace mc

#endif
