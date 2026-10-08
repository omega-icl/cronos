// Copyright (C) Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

#ifndef CRONOS__ODESLV_BASE_HPP
#define CRONOS__ODESLV_BASE_HPP

#undef  CRONOS__ODESLV_BASE_DEBUG

#include <limits>
#include <utility>
#include <cstdint>
#include <new>
#include <sstream>
#include <stdexcept>
#include <cassert>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <vector>
#include <list>

#include "ffmodel.hpp"

namespace mc
{

//! @brief Allocate an array of n FFVar (2026-10-03).  The bound makes the size computation of new[] provably free of
//! overflow: GCC (with LTO, -Walloc-size-larger-than) otherwise warns on its compiler-generated overflow branch,
//! which calls operator new[] with SIZE_MAX.  An impossible size throws std::bad_array_new_length, as new[] would.
inline FFVar*
odeslv_new_vars
( size_t const n )
{
  if( n > static_cast<size_t>( PTRDIFF_MAX ) / sizeof( FFVar ) - 1 ) throw std::bad_array_new_length();
  return new FFVar[n];
}
//! @brief C++ base class for computing solutions of parametric ODEs using continuous-time real-valued integration.
////////////////////////////////////////////////////////////////////////
//! mc::ODESLV_BASE is a C++ base class for computing solutions of
//! parametric ordinary differential equation (ODEs) using
//! continuous-time real-valued integration.
////////////////////////////////////////////////////////////////////////
class ODESLV_BASE
: public virtual FFModel
{
 public:
  /** @defgroup ODESLV Continuous-time real-valued integration of parametric ODEs
   *  @{
   */
  //! @brief Default constructor
  ODESLV_BASE();

  //! @brief Virtual destructor
  virtual ~ODESLV_BASE();

  //! @brief Integration results at a given time instant
  struct Results
  {
    //! @brief Constructors
    Results
      ( double const& tk, unsigned const nx1, double const* x1,
        unsigned const nx2=0, double const* x2=nullptr ):
      t( tk )
      { x.resize( nx1+nx2 );
        unsigned int ix=0;
        for( ; ix<nx1; ix++ )      x[ix] = x1[ix];
        for( ; ix<nx1+nx2; ix++ )  x[ix] = x2[ix-nx1]; }
    Results
      ( Results const& res ):
      t( res.t ), x( res.x )
      {}
    //! @brief Time point
    double t;
    //! @brief Solution point
    std::vector<double> x;
  };
  /** @} */

  //! @brief Why the last setup() failed to extract the model; empty after a successful one.  Usually the
  //! dynamic_form() blocker: the model is distributed, or a boundary-value problem in the evolution direction, or
  //! its derivatives are coupled through a mass matrix.
  std::string const& extract_error
    ()
    const
    { return _extractError; }

  //! @brief Set the RHS Jacobian sparsity size (was BASE_DE::set_sparse)
  bool set_sparse
    ();

  //! @brief Display intermediate results (was BASE_DE::_print_interm)
  template<typename U>
  static void _print_interm
    ( size_t const nx, U const* x, std::string const& var, std::ostream& os=std::cout );

  //! @brief Display intermediate results
  template<typename U>
  static void _print_interm
    ( double const& t, size_t const nx, U const* x, std::string const& var, std::ostream& os=std::cout );

  //! @brief Display intermediate results
  template<typename U, typename V>
  static void _print_interm
    ( double const& t, size_t const nx, U const* x, V const& r, std::string const& var,
      std::ostream& os=std::cout );

  //! @brief ODESLV implements EXPLICIT transitions (add_transition): x_i(tau^+) = g(x(tau^-), p)
  bool _transitions_supported() const override { return true; }

  //! @brief Take the model declaration of @p IVP (_usr and the control flags, shallow on the same user DAG, which
  //! is not owned) so that a following setup() yields an independent solver.  Everything extraction derives --
  //! dimensions, _ndxSEN, minted levels, the stage partition -- is recomputed by that setup(), not copied.
  //! (2026-09-29: the legacy description path -- set_dag/set_time/set_state/set_parameter/set_constant/
  //! set_differential/set_initial/set_quadrature/set_function -- is retired: a model is declared through FFModel.)
  void _copy_description
    ( ODESLV_BASE const& IVP )
    {
      FFModel::set( IVP._usr._dagUsr );
      _usr = IVP._usr;
      _decisionFlag = IVP._decisionFlag;
    }

  //! @brief Integrator status (was BASE_DE::STATUS -- a solver's verdict, not a model's)
  enum STATUS{
     NORMAL=0,  //!< Normal execution
     FAILURE,   //!< Integration breakdown
     FATAL      //!< Interruption due to errors in third-party libraries
  };

  //! @brief Number of parameters (was BASE_DE::np)
  //! @brief The parameter vector the solver actually integrates against: the time-invariant inputs and,
  //! for each time-varying input, its minted per-element nodal values -- in _mP order, the order np()
  //! counts and the sensitivity directions index.  Exposed because fdiff() needs FFVars to differentiate
  //! against and the minted levels are reachable nowhere else.
  std::vector<FFVar> const& var_parameter
    () const
    { return _mP; }

  //! @brief The parameter indices the sensitivity runs over -- FFModel's canonical control order, so
  //! direction j is DOF j of the control vector.  Every parameter when no control is registered.
  std::vector<size_t> const& sensitivity_index
    () const override
    { return _ndxSEN; }

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

  //! @brief ODESLV substitutes fixed-input values at extraction (setup), not per call.
  bool _fixed_inputs_applied_per_call
    () const override
    { return false; }

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

  size_t np
    ()
    const
    { return _np; }

  //! @brief Revision of the ODESLV headers (odeslv_base.hpp, odeslvs_base.hpp, odeslv_cvodes.hpp, odeslvs_cvodes.hpp,
  //! which change together): bump it with any change to one of them.  Since 2026-10-07, so that a sweep log can show
  //! which ODESLV headers it ran (ODESLV_revision prints it).
  static constexpr char const* HEADER_ID
    = "odeslv  rev1  2026-10-07";

  //! @brief The revision of the ODESLV headers (HEADER_ID), e.g. for a bug report
  static char const* revision() { return HEADER_ID; }

  //! @brief Number of state functions (was BASE_DE::nf)
  size_t nf
    ()
    const
    { return _nf; }

  //! @brief The block of output @p ndx (in add_output() order) in val_function(): (first, count).  A scalar output
  //! takes ONE value; an output DISTRIBUTED over the evolution direction takes one value per kept STAGE TIME
  //! (ODESLV has no collocation nodes), in time order.  (max, 0) if setup() has not run or @p ndx is out of range --
  //! the convention of OCFESLV::blk_fct, so a script reads the outputs of either solver alike (WORKPLAN 2.2, 1.7).
  std::pair<size_t,size_t> blk_fct
    ( size_t const ndx )
    const
    {
      return this->is_setup() && ndx < _blkFct.size()? _blkFct[ndx]
                                                     : std::make_pair( std::numeric_limits<size_t>::max(), size_t(0) );
    }

  //! @brief Number of states (was BASE_DE::nx)
  size_t nx
    ()
    const
    { return _nx; }

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

protected:
  // ---- what the old base class provided, now this class's own --------------------------------------------
  size_t _ns = 0;        //!< number of stages: _dT.size()-1
  size_t _nx = 0;        //!< number of states
  size_t _nx0 = 0;       //!< number of initial values
  size_t _nc = 0;        //!< number of constants
  size_t _np = 0;        //!< number of parameters
  //! @brief The parameter indices the sensitivity runs over, in FFModel's CANONICAL CONTROL ORDER, so
  //! direction j is DOF j of the control vector and control_columns() slices the gradient.  Empty registry
  //! -> 0.._np-1, i.e. every parameter, which is the behaviour before controls were honoured.
  std::vector<size_t> _ndxSEN;
  size_t _nsen = 0;      //!< number of sensitivity DIRECTIONS (_ndxSEN.size()); _np counts PARAMETERS
  bool _resolve_sensitivity_index();

  //! @brief The control registry changed.  If the model is already extracted, re-resolve the sensitivity
  //! directions: that reads only extraction outputs (_mP, _mInpLevels), and CVODES sizes its sensitivity memory
  //! per solve, so no re-extraction is needed.  Lets a DAG operation reset the registry after setup().
  void _on_controls_changed
    () override
    { if( _dag ) _resolve_sensitivity_index(); }

  size_t _nq = 0;        //!< number of quadratures
  size_t _nf = 0;        //!< number of state functions (ROWS of val_function(): a distributed output has one per stage time)
  //! the block [first, count] of each user output (in add_output() order) among the functions / val_function()
  std::vector< std::pair<size_t,size_t> > _blkFct;
  double _t = 0.;        //!< current time
  size_t _istg = 0;      //!< current stage
  size_t _nnzjac = 0;    //!< nonzeros in the Jacobian

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
  std::vector<FFVar>                   _mLHS;   //!< the coefficient multiplying each derivative (1 if explicit)
  std::vector<std::vector<FFVar>>      _mIC;    //!< initial values (one entry: no state discontinuities)
  std::vector<std::vector<FFVar>>      _mDE;    //!< right-hand sides (one entry)
  std::vector<std::vector<FFVar>>      _mQUAD;  //!< quadrature integrands (one entry)
  std::vector<std::map<size_t,FFVar>>  _mFCT;   //!< functions, per stage: index -> expression
  std::string _extractError;                    //!< why the last extraction failed

  //! @brief Fill the description above from this model.  False, with _extractError set, when the model is not one
  //! an integrator can take (dynamic_form().blocker says why) or a row is not solvable for the derivative it
  //! defines.  Folded in here rather than put in a class of its own: these are THIS solver's working copies, and
  //! a second consumer -- a DAESLV on IDAS, which would want the _mLHS/_mDE pair as a residual -- is the moment
  //! to hoist them, not before.
  bool _extract_from_model
    ();

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

  //! @brief local copy of DAG
  FFGraph* _dag;

  //! @brief local copy of initial value functions
  std::vector<std::vector<FFVar>> _vIC;

  //! @brief local copy of right-hand-side functions
  std::vector<std::vector<FFVar>> _vRHS;

  //! @brief local copy of quadrature functions
  std::vector<std::vector<FFVar>> _vQUAD;

  //! @brief local copy of output functions
  std::vector<std::vector<FFVar>> _vFCT;

  //! @brief local copy of constants
  FFVar* _pC;

  //! @brief local copy of parameters
  FFVar* _pP;

  //! @brief local copy of time/independent variable
  FFVar* _pT;

  //! @brief local copy of state variables
  FFVar* _pX;

  //! @brief local copy of quadrature variables
  FFVar* _pQ;

  //! @brief Sugraph of ODE RHS function
  FFSubgraph _opRHS;

  //! @brief Subgraph of ODE RHS Jacobian
  FFSubgraph _opJAC;

  //! @brief Subgraph of quadrature RHS function
  FFSubgraph _opQUAD;

  //! @brief const pointer to RHS function in current stage of ODE system
  FFVar const* _pRHS;

  //! @brief const pointer to quadrature integrand in current stage of ODE system
  FFVar const* _pQUAD;

  //! @brief const pointer to IC function in current stage of ODE system
  FFVar const* _pIC;

  //! @brief sparse representation of RHS Jacobian in current stage of ODE system
  std::tuple< unsigned, unsigned*, unsigned*, FFVar* > _pJAC;

  //! @brief sparse representation of RHS Jacobian in current stage of ODE system: i-th entry is the index in data where the first non-zero matrix entry of the i-th column is stored (length NEQ + 1), last entry is number of non-zeros
  std::vector< size_t > _pJACCOLNDX;

  //! @brief number of variables for DAG evaluation
  unsigned _nVAR;

  //! @brief number of variables for DAG evaluation (without quadratures)
  unsigned _nVAR0;

  //! @brief array of variables for DAG evaluation
  FFVar* _pVAR;

  //! @brief array of variable values for DAG evaluation
  double* _DVAR;

  //! @brief pointer to time value **DO NOT FREE**
  double* _Dt;

  //! @brief pointer to state values **DO NOT FREE**
  double* _Dx;

  //! @brief pointer to parameter values **DO NOT FREE**
  double* _Dp;

  //! @brief pointer to quadrature values **DO NOT FREE**
  double* _Dq;

  //! @brief vector to hold function values
  std::vector<double> _Df;

  //! @brief vector to hold Jacobian evaluation results
  std::vector<double> _DJAC;

  //! @brief storage vector for DAG evaluation
  std::vector<double> _DWRK;

  //! @brief Function setting up local DAG
  bool _SETUP
    ();

  //! @brief Function setting up DAG of IVP
  bool _SETUP
    ( ODESLV_BASE const& IVP );

  //! @brief Function converting integrator array to internal format
  template <typename REALTYPE>
  static void _vec2D
    ( REALTYPE const* vec, unsigned const n, double* d );

  //! @brief Function converting integrator array to internal format
  template <typename REALTYPE>
  static void _D2vec
    ( double const* d, unsigned const n, REALTYPE* vec );

  //! @brief Function to initialize state integration
  bool _INI_D_STA
    ( double const* p, double const* c );

  //! @brief Function to finalize state integration
  bool _END_D_STA
    ();

  //! @brief Function to retreive state bounds
  template <typename REALTYPE>
  void _GET_D_STA
    ( REALTYPE const* x, REALTYPE const* q );

  //! @brief Function to set state/quadrature at initial time
  bool _IC_D_SET
    ();

  //! @brief Function to initialize quaratures
  template <typename REALTYPE>
  bool _IC_D_QUAD
    ( REALTYPE* vec );

  //! @brief Function to initialize states
  template <typename REALTYPE>
  bool _IC_D_STA
    ( double const& t, REALTYPE* vec );

  //! @brief Function to reset state at intermediate time
  bool _CC_D_SET
    ( unsigned const iIC );

  //! @brief Function to reinitialize state at intermediate time
  template <typename REALTYPE>
  bool _CC_D_STA
    ( double const& t, REALTYPE* vec );

  //! @brief Function to set RHS and QUAD pointers
  bool _RHS_D_SET
    ( unsigned const iRHS, unsigned const iQUAD );

  //! @brief Function to set RHS and QUAD pointers
  bool _RHS_D_SET
    ();

  //! @brief Function to calculate the ODE RHS
  template <typename REALTYPE>
  bool _RHS_D_STA
    ( double const& t, REALTYPE const* x, REALTYPE* xdot );

  //! @brief Function to calculate the ODE quadrature
  template <typename REALTYPE>
  bool _RHS_D_QUAD
    ( double const& t, REALTYPE const* x, REALTYPE* qdot );

  //! @brief Function to calculate the ODE Jacobian
  template <typename REALTYPE, typename INDEXTYPE>
  bool _JAC_D_STA
    ( double const& t, REALTYPE const* x, REALTYPE* jac, INDEXTYPE* ptrs, INDEXTYPE* vals );

  //! @brief Function to calculate the ODE Jacobian
  template <typename REALTYPE>
  bool _JAC_D_STA
    ( double const& t, REALTYPE const* x, REALTYPE** jac );

  //! @brief Function to calculate the functions at intermediate/end point
  bool _FCT_D_STA
    ( unsigned const iFCT, double const& t );

  //! @brief Record results in file <a>bndrec</a>, with accuracy of <a>iprec</a> digits
  template <typename VRES>
  static void _record
    ( std::ofstream& ofile, VRES const& bnd, unsigned const iprec=5 );

  //! @brief Block default compiler methods
  ODESLV_BASE( ODESLV_BASE const& ) = delete;
  ODESLV_BASE& operator=( ODESLV_BASE const& ) = delete;
};

inline 
ODESLV_BASE::ODESLV_BASE
()
: _dag(nullptr),
  _pC(nullptr), _pP(nullptr), _pT(nullptr), _pX(nullptr), _pQ(nullptr),
  _pJAC(0,nullptr,nullptr,nullptr),
  _nVAR(0), _nVAR0(0), _pVAR(nullptr),
  _DVAR(nullptr), _Dt(nullptr), _Dx(nullptr), _Dp(nullptr), _Dq(nullptr)
{}

inline
ODESLV_BASE::~ODESLV_BASE
()
{
  /* DO NOT FREE _pRHS, _pQUAD, _pIC */
  delete[] std::get<1>(_pJAC);  std::get<1>(_pJAC) = nullptr;
  delete[] std::get<2>(_pJAC);  std::get<2>(_pJAC) = nullptr;
  delete[] std::get<3>(_pJAC);  std::get<3>(_pJAC) = nullptr;
  delete[] _pVAR;
  delete[] _DVAR;
  delete[] _pX;
  delete[] _pQ;
  delete[] _pP;
  delete[] _pC;
  delete   _pT;
  delete   _dag;
}

inline
bool
ODESLV_BASE::_SETUP
()
{
  if( !_extract_from_model() ) return false;
  // The control registry -> parameter indices (without it _ndxSEN would be empty and _nsen zero, silently
  // disabling the sensitivity solve).
  if( !_resolve_sensitivity_index() ) return false;
  delete _dag; _dag = new FFGraph;
#ifdef CRONOS__ODESLV_BASE_DEBUG
  std::cout << "ODESLV_BASE:: Original DAG: " << FFModel::_dag << std::endl;
  std::cout << "ODESLV_BASE:: Copied DAG:   " << _dag << std::endl;
#endif

  delete[] _pC; _pC = nullptr;
  if( _nc ){
    _pC  = odeslv_new_vars( _nc );
    _dag->insert( FFModel::_dag, _nc, _mC.data(), _pC );
  }

  delete[] _pP; _pP = nullptr;
  if( _np ){
    _pP  = odeslv_new_vars( _np );
    _dag->insert( FFModel::_dag, _np, _mP.data(), _pP );
  }

  delete _pT; _pT = nullptr;
  if( _mT.size() ){
    _pT  = new FFVar;
    _dag->insert( FFModel::_dag, 1, _mT.data(), _pT );
  }

  delete[] _pX; _pX = nullptr;
  if( _nx ){
    _pX  = odeslv_new_vars( _nx );
    _dag->insert( FFModel::_dag, _nx, _mX.data(), _pX );
  }

  delete[] _pQ; _pQ = nullptr;
  if( _nq ){
    _pQ  = odeslv_new_vars( _nq );
    _dag->insert( FFModel::_dag, _nq, _mQ.data(), _pQ );
  }

  _vIC.assign( _mIC.size(), std::vector<FFVar>(_nx0) );
  size_t k = 0;
  for( auto const& vICk : _mIC )
    _dag->insert( FFModel::_dag, _nx0, vICk.data(), _vIC[k++].data() );

  _vRHS.assign( _mDE.size(), std::vector<FFVar>(_nx) );
  k = 0;
  for( auto const& vRHSk : _mDE )
    _dag->insert( FFModel::_dag, _nx, vRHSk.data(), _vRHS[k++].data() );

  _vQUAD.assign( _mQUAD.size(), std::vector<FFVar>(_nq) );
  k = 0;
  for( auto const& vQUADk : _mQUAD )
    _dag->insert( FFModel::_dag, _nq, vQUADk.data(), _vQUAD[k++].data() );

  if( _mFCT.size() == 1 ){
    _vFCT.assign( _ns+1, std::vector<FFVar>() ); // empty vectors except last stage
    _vFCT[_ns].resize( _nf );
    _dag->insert( FFModel::_dag, _mFCT[0], _vFCT[_ns] );
  }
  else{// if( _mFCT.size() == _ns+1 ){
    _vFCT.assign( _mFCT.size(), std::vector<FFVar>(_nf, 0. ) );
    k = 0;
    for( auto const& mFCTk : _mFCT )
      _dag->insert( FFModel::_dag, mFCTk, _vFCT[k++] );
  }

  return true;
}

inline
bool
ODESLV_BASE::_SETUP
( ODESLV_BASE const& IVP )
{
  delete _dag; _dag = new FFGraph;
#ifdef CRONOS__ODESLV_BASE_DEBUG
  std::cout << "ODESLV_BASE:: IVP DAG:   " << IVP._dag << std::endl;
#endif
  
  delete _pT; _pT = nullptr;
  if( IVP._pT ){
    _pT  = new FFVar;
    _dag->insert( IVP._dag, 1, IVP._pT, _pT );
  }

  delete[] _pX; _pX = nullptr;
  if( _nx ){
    _pX  = odeslv_new_vars( _nx );
    _dag->insert( IVP._dag, _nx, IVP._pX, _pX );
  }

  delete[] _pQ; _pQ = nullptr;
  if( _nq ){
    _pQ  = odeslv_new_vars( _nq );
    _dag->insert( IVP._dag, _nq, IVP._pQ, _pQ );
  }

  delete[] _pC; _pC = nullptr;
  if( _nc ){
    _pC  = odeslv_new_vars( _nc );
    _dag->insert( IVP._dag, _nc, IVP._pC, _pC );
  }

  delete[] _pP; _pP = nullptr;
  if( _np ){
    _pP  = odeslv_new_vars( _np );
    _dag->insert( IVP._dag, _np, IVP._pP, _pP );
  }

  _vIC.assign( IVP._vIC.size(), std::vector<FFVar>(_nx0) );
  size_t k = 0;
  for( auto const& vICk : IVP._vIC )
    _dag->insert( IVP._dag, _nx0, vICk.data(), _vIC[k++].data() );
    
  _vQUAD.assign( IVP._vQUAD.size(), std::vector<FFVar>(_nq) );
  k = 0;
  for( auto const& vQUADk : IVP._vQUAD )
    _dag->insert( IVP._dag, _nq, vQUADk.data(), _vQUAD[k++].data() );

  _vRHS.assign( IVP._mDE.size(), std::vector<FFVar>(_nx) );
  k = 0;
  for( auto const& vRHSk : IVP._vRHS )
    _dag->insert( IVP._dag, _nx, vRHSk.data(), _vRHS[k++].data() );

  _vFCT.assign( IVP._vFCT.size(), std::vector<FFVar>(_nf, 0. ) );
  k = 0;
  for( auto const& vFCTk : IVP._vFCT )
    _dag->insert( IVP._dag, vFCTk, _vFCT[k++] );
  
  return true;
}

template <typename VRES> 
inline
void
ODESLV_BASE::_record
( std::ofstream& ofile, VRES const& res, unsigned const iprec )
{
  if( !ofile ) return;

  // Specify format
  ofile << std::right << std::scientific << std::setprecision(iprec);

  // Record computed states at stage times
  auto it = res.begin(), it0 = it;
  for( ; it != res.end(); ++it ){
    if( it != res.begin() && it->t == it0->t )
      ofile << std::endl;
    ofile << std::setw(iprec+9) << it->t;
    for( auto const& xi : it->x )
      ofile << std::setw(iprec+9) << xi;
    ofile << std::endl;
    it0 = it;
  }
}

template <typename REALTYPE>
inline
void
ODESLV_BASE::_vec2D
( REALTYPE const* vec, unsigned const n, double* d )
{
  for( unsigned i=0; i<n; i++  ) d[i] = vec[i];
  return;
}

template <typename REALTYPE>
inline
void
ODESLV_BASE::_D2vec
( double const* d, unsigned const n, REALTYPE* vec )
{
  for( unsigned i=0; i<n; i++  ) vec[i] = d[i];
  return;
}

inline
bool
ODESLV_BASE::_INI_D_STA
( double const* p, double const* c )
{
  // Set constants
  for( unsigned ic=0; c && ic<_nc; ++ic ) _pC[ic].set( c[ic] );

  // Size and set DAG evaluation arrays
  _nVAR0 = _nx+_np+1;
  _nVAR  = _nVAR0+_nq;
  delete[] _pVAR; _pVAR = odeslv_new_vars( _nVAR );
  delete[] _DVAR; _DVAR = new double[_nVAR];
  for( unsigned ix=0; ix<_nx; ix++ ) _pVAR[ix] = _pX[ix];
  for( unsigned ip=0; ip<_np; ip++ ) _pVAR[_nx+ip] = _pP[ip];
  _pVAR[_nx+_np] = (_pT? *_pT: 0. );
  for( unsigned iq=0; iq<_nq; iq++ ) _pVAR[_nx+_np+1+iq] = _pQ?_pQ[iq]:0.;

  _Dx = _DVAR;
  _Dp = _Dx + _nx;
  _Dt = _Dp + _np;
  _Dq = _Dt + 1;
  for( unsigned ip=0; ip<_np; ip++ ) _Dp[ip] = p[ip];

  _Df.assign( _nf, 0. );

  return true;
}

inline
bool
ODESLV_BASE::_END_D_STA
()
{
  // Set constants
  for( unsigned ic=0; ic<_nc; ++ic ) _pC[ic].unset();

  return true;
}

template <typename REALTYPE>
inline
void
ODESLV_BASE::_GET_D_STA
( REALTYPE const* x, REALTYPE const* q )
{
  _vec2D( x, _nx, _Dx );
  if( q ) _vec2D( q, _nq, _Dq );
}

inline
bool
ODESLV_BASE::_IC_D_SET
()
{
  //std::cout << "ENTERING _IC_D_SET\n";
  //std::cout << _vIC.size() << " " << _nx0 << " " << _nx << std::endl;
  if( !_vIC.size() || _nx0 != _nx ) return false;
  _pIC = _vIC.at(0).data();
  return true;
}

template <typename REALTYPE>
inline
bool
ODESLV_BASE::_IC_D_STA
( double const& t, REALTYPE* x )
{
  *_Dt = t; // current time
  _dag->eval( _nx, _pIC, (double*)x, _np+1, _pVAR+_nx, _DVAR+_nx );
  return true;
}

template <typename REALTYPE>
inline
bool
ODESLV_BASE::_IC_D_QUAD
( REALTYPE* q )
{
  for( unsigned iq=0; iq<_nq; iq++ ) q[iq] = 0.;
  return true;
}

inline
bool
ODESLV_BASE::_CC_D_SET
( unsigned const iIC )
{
  if( _vIC.size() <= iIC || _nx0 != _nx ) return false;
  _pIC = _vIC.at( iIC ).data();
  return true;
}

template <typename REALTYPE>
inline
bool
ODESLV_BASE::_CC_D_STA
( double const& t, REALTYPE* x )
{
  *_Dt = t; // current time
  _vec2D( x, _nx, _Dx ); // current state
  _dag->eval( _nx, _pIC, (double*)x, _nVAR0, _pVAR, _DVAR );
  return true;
}

inline
bool
ODESLV_BASE::_RHS_D_SET
( unsigned const iRHS, unsigned const iQUAD )
{
  if( _vRHS.size() <= iRHS ) return false;
  _pRHS = _vRHS.at( iRHS ).data();

  if( _nq && _vQUAD.size() <= iQUAD ) return false;
  _pQUAD = _nq? _vQUAD.at( iQUAD ).data(): nullptr;

  // Generate Jacobian using sparse forward AD
  delete[] std::get<1>(_pJAC); delete[] std::get<2>(_pJAC); delete[] std::get<3>(_pJAC);
  _pJAC = _dag->SFAD( _nx, _pRHS, _nx, _pX ); // Jacobian in sparse format, ordered columnwise
  _pJACCOLNDX.resize( _nx+1 );
  for( unsigned ie=0, ic=0; ie<std::get<0>(_pJAC); ++ie ){
#ifdef CRONOS__ODESLV_BASE_DEBUG
    std::cout << "  JAC[" << std::get<1>(_pJAC)[ie] << ", " << std::get<2>(_pJAC)[ie] << "]" << std::endl;
#endif
    for( ; std::get<2>(_pJAC)[ie] >= ic; ++ic ){
      _pJACCOLNDX[ic] = ie;
#ifdef CRONOS__ODESLV_BASE_DEBUG
      std::cout << "  JACCOLNDX[" << ic << "] = " << ie << std::endl;
#endif
    }
  }
  _pJACCOLNDX[_nx] = std::get<0>(_pJAC);
#ifdef CRONOS__ODESLV_BASE_DEBUG
  std::cout << "  JACCOLNDX[" << _nx << "] = " << std::get<0>(_pJAC) << std::endl;
  std::cout << "PAUSED - <1> TO CONTINUE"; int dum; std::cin >> dum;
#endif

  return _RHS_D_SET();
}

inline
bool
ODESLV_BASE::_RHS_D_SET
()
{
  _opRHS  = _dag->subgraph( _nx, _pRHS );

  if( _pQUAD ) _opQUAD = _dag->subgraph( _nq, _pQUAD );

  _opJAC = _dag->subgraph( std::get<0>(_pJAC), std::get<3>(_pJAC) );
  _DJAC.resize( std::get<0>(_pJAC) );

  return true;
}

template <typename REALTYPE>
inline
bool
ODESLV_BASE::_RHS_D_STA
( double const& t, REALTYPE const* x, REALTYPE* xdot )
{
  if( !_pRHS ) return false;
  *_Dt = t; // current time
  _vec2D( x, _nx, _Dx ); // current state
  _dag->eval( _opRHS, _DWRK, _nx, _pRHS, (double*)xdot, _nVAR0, _pVAR, _DVAR );
  return true;
}

template <typename REALTYPE>
inline
bool
ODESLV_BASE::_RHS_D_QUAD
( double const& t, REALTYPE const* x, REALTYPE* qdot )
{
  if( !_pQUAD ) return false;
  // 2026-09-28: load the callback's OWN (t, x) -- relying on the last RHS call made the Jacobian / quadrature
  // stale whenever another evaluation came in between (the forward replay inside CVodeB: wrong J, replay overrun)
  *_Dt = t;  _vec2D( x, _nx, _Dx );
  _dag->eval( _opQUAD, _DWRK, _nq, _pQUAD, (double*)qdot, _nVAR0, _pVAR, _DVAR );
  return true;
}

template <typename REALTYPE>
inline
bool
ODESLV_BASE::_JAC_D_STA
( double const& t, REALTYPE const* x, REALTYPE** jac )
{
  // 2026-09-28: load the callback's OWN (t, x) -- relying on the last RHS call made the Jacobian / quadrature
  // stale whenever another evaluation came in between (the forward replay inside CVodeB: wrong J, replay overrun)
  *_Dt = t;  _vec2D( x, _nx, _Dx );
  _dag->eval( _opJAC, _DWRK, std::get<0>(_pJAC), std::get<3>(_pJAC),
              _DJAC.data(), _nVAR0, _pVAR, _DVAR );
  for( unsigned ie=0; ie<std::get<0>(_pJAC); ++ie ){
    jac[std::get<2>(_pJAC)[ie]][std::get<1>(_pJAC)[ie]] = _DJAC[ie];
#ifdef CRONOS__ODESLV_BASE_DEBUG
    std::cout << "  jac[" << std::get<1>(_pJAC)[ie] << ", "
              << std::get<2>(_pJAC)[ie] << "] = " << _DJAC[ie] << std::endl;
#endif
  }
  return true;
}

#if defined( CRONOS__WITH_KLU )
template <typename REALTYPE, typename INDEXTYPE>
inline
bool
ODESLV_BASE::_JAC_D_STA
( double const& t, REALTYPE const* x, REALTYPE* jac, INDEXTYPE* ptrs, INDEXTYPE* vals )
{
  // 2026-09-28: load the callback's OWN (t, x) -- relying on the last RHS call made the Jacobian / quadrature
  // stale whenever another evaluation came in between (the forward replay inside CVodeB: wrong J, replay overrun)
  *_Dt = t;  _vec2D( x, _nx, _Dx );
  _dag->eval( _opJAC, _DWRK, std::get<0>(_pJAC), std::get<3>(_pJAC),
              (double*)jac, _nVAR0, _pVAR, _DVAR );
  for( unsigned ie=0; ie<std::get<0>(_pJAC); ++ie ){
    vals[ie] = (INDEXTYPE)std::get<1>(_pJAC)[ie];
#ifdef CRONOS__ODESLV_BASE_DEBUG
    std::cout << "  jac[" << ie << "] = " << jac[ie] << std::endl;
    std::cout << "  vals[" << ie << "] = " << vals[ie] << std::endl;
#endif
  }
  for( unsigned ic=0; ic<=_nx; ++ic ){
    ptrs[ic] = (INDEXTYPE) _pJACCOLNDX[ic];
#ifdef CRONOS__ODESLV_BASE_DEBUG
    std::cout << "  ptrs[" << ic << "] = " << ptrs[ic] << std::endl;
#endif
  }
  return true;
}
#endif

inline
bool
ODESLV_BASE::_FCT_D_STA
( unsigned const iFCT, double const& t )
{
  if( !_nf || _vFCT.at( iFCT ).empty() ) return true; // nothing to do if no function

  *_Dt = t; // current time
  FFVar const* pFCT = _vFCT.at( iFCT ).data();
  //std::cout << "evaluating functions @" << t << std::endl;
  static double const one = 1.;
  _dag->eval( _nf, pFCT, _Df.data(), _nVAR, _pVAR, _DVAR, &one );//iFCT? &one: nullptr );
  
  return true;
}


//! @brief Resolve FFModel's control registry into parameter indices of _mP.  A control that is distributed
//! over the evolution variable contributes one parameter per element (_mInpLevels), in element order, which
//! is the same block FFModel::control_ndof() counts -- the two are cross-checked here rather than assumed.
inline bool
ODESLV_BASE::_resolve_sensitivity_index
()
{
  _ndxSEN.clear();
  if( controls().empty() ){                       // nothing registered: differentiate w.r.t. everything
    _ndxSEN.resize( _np );
    for( size_t i=0; i<_np; ++i ) _ndxSEN[i] = i;
    _nsen = _np;
    return true;
  }
  for( auto const& [u,spec] : controls() ){       // canonical control order: offsets ascending
    std::vector<FFVar> dofs;
    auto const itl = _mInpLevels.find( u.name() );
    if( itl != _mInpLevels.cend() ) dofs = itl->second;
    else                            dofs.assign( 1, u );
    if( dofs.size() != spec.ndof ){
      _extractError = "the control " + u.name() + " has " + std::to_string( spec.ndof )
                    + " DOFs in its declared function space but " + std::to_string( dofs.size() )
                    + " parameters in the extracted model";
      return false;
    }
    for( auto const& d : dofs ){
      bool found = false;
      for( size_t i=0; i<_np; ++i )
        if( _mP[i].id().second == d.id().second ){ _ndxSEN.push_back( i ); found = true; break; }
      if( !found ){
        _extractError = "the registered control " + u.name()
                      + " is not a parameter of the extracted model";
        return false;
      }
    }
  }
  _nsen = _ndxSEN.size();
  return true;
}

inline
bool
ODESLV_BASE::_split_in
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
ODESLV_BASE::_extract_from_model
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
  if( !dynamic_form().blocker.empty() ){ _extractError = dynamic_form().blocker; return false; }

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
  _mLHS.assign( _nx, FFVar( 1. ) );
  std::vector<bool> got( _nx, false );
  for( auto const& eqn : var_equation() ){
    if( !eqn.opt || eqn.opt->role != EqnRole::INTERIOR ) continue;
    for( size_t i = 0; i < _nx; ++i ){
      if( got[i] ) continue;
      FFVar L, R;
      if( !_split_in( eqn.var, OpP( _mX[i], tvar ), L, R ) ) continue;
      if( L.cst() && L.num().val() == 0. ) continue;
      _mLHS[i] = L;
      rhs[i]   = ( L.cst() && L.num().val() == 1. )? _in_model_dag( -R ): _in_model_dag( -R / L );
      got[i]   = true;
    }
  }
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
      FFVar expr( 0. );
      for( size_t j = 0; j < nn; ++j ){
        FFVar Lj( 1. );
        for( size_t k = 0; k < nn; ++k )
          if( k != j ) Lj *= ( sc - sn[k] ) / ( sn[j] - sn[k] );
        expr += Lj * lv[elem*nn + j];
      }
      rep.push_back( expr );
    }
    return FFModel::_dag->compose( ex, tgt, rep ); };
  if( tv_inputs.empty() )
    _mDE.assign( 1, rhs );
  else{                                                    // one entry per stage, with that element's levels in
    _mDE.assign( _ns, std::vector<FFVar>( _nx ) );
    for( size_t k = 0; k < _ns; ++k ){
      std::vector<FFVar> const done = subst_inputs( rhs, elem_of_stage( k ) );
      for( size_t i = 0; i < _nx; ++i ) _mDE[k][i] = done[i];
    }
  }

  std::vector<FFVar> ic( _nx, FFVar( 0. ) );
  for( auto const& eqn : var_equation() ){
    if( !eqn.opt || eqn.opt->role != EqnRole::INITIAL ) continue;
    for( size_t i = 0; i < _nx; ++i ){
      FFVar L, R;
      if( !_split_in( subst1( eqn.var ), _mX[i], L, R ) ) continue;
      if( L.cst() && L.num().val() == 0. ) continue;
      ic[i] = ( L.cst() && L.num().val() == 1. )? _in_model_dag( -R ): _in_model_dag( -R / L );
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
  return true;
}

template <typename U>
inline
void
ODESLV_BASE::_print_interm
( double const& t, size_t const nx, U const* x, std::string const& var,
  std::ostream& os )
{
  os << " @t = " << std::scientific << std::setprecision(6)
                 << std::left << t << " :" << std::endl;
  _print_interm( nx, x, var, os );
  return;
}

template <typename U, typename V>
inline
void
ODESLV_BASE::_print_interm
( double const& t, size_t const nx, U const* x, V const& r,
  std::string const& var, std::ostream& os )
{
  os << " @t = " << std::scientific << std::setprecision(6)
                 << std::left << t << " :" << std::endl;
  _print_interm( nx, x, var, os );
  os << " " << "R" << var << " =" << r << std::endl;
  return;
}

template <typename U>
inline
void
ODESLV_BASE::_print_interm
( size_t const nx, U const* x, std::string const& var, std::ostream& os )
{
  if( !x ) return;
  for( size_t ix=0; ix<nx; ix++ )
    os << " " << var << "[" << ix << "] = " << x[ix] << std::endl;
  return;
}


inline
bool
ODESLV_BASE::set_sparse
()
{
  _nnzjac = 0;
  // ODE systems only.  The old "_nd != _nx" half of this guard is gone: the extraction already refuses a model
  // whose rows are not explicitly solvable for their derivatives, so every extracted state IS differential.
  if( _nx0 != _nx ) return false;

  // Intermediates
  std::vector<FFVar> vVAR = _mX;
  vVAR.insert( vVAR.end(), _mP.cbegin(), _mP.cend() );
  vVAR.push_back( _mT.size()? _mT.front(): 0. );
  size_t const nVAR = vVAR.size();
  std::vector<FFDep> depSTA( nVAR, 0 ), depRHS( _nx );
  for( size_t ix=0; ix<_nx; ++ix ) depSTA[ix].indep( ix );

  // RHS dependencies
  for( size_t is=0; is<_ns; ++is ){
    size_t const pos_rhs  = ( _mDE.size()<=1?  0: is );
    FFVar const* pRHS  = ( _mDE.empty()? nullptr: _mDE[pos_rhs].data() );
    if( !pRHS ) return false;
    // FFModel::_dag, NOT _dag: the staging vectors live in the MODEL's graph, while _dag is this solver's own
    // working copy.  The two can carry the same indices, so evaluating one graph's expressions against the other
    // reads whatever happens to sit there -- which is the segfault valgrind caught inside this function.
    if( !_nc ) FFModel::_dag->eval( _nx, pRHS, depRHS.data(), nVAR, vVAR.data(), depSTA.data() );
    else       FFModel::_dag->eval( _nx, pRHS, depRHS.data(), nVAR, vVAR.data(), depSTA.data(),
                                    _nc, _mC.data(), std::vector<FFDep>(_nc, 0 ).data() );
    size_t nnz = 0;
    for( size_t ix=0; ix<_nx; ix++ ){
#ifdef CRONOS__DEBUG__BASE_DE
      std::cout << "RHS[" << is << "][" << ix << "]: " << depRHS[ix] << std::endl;
#endif
      nnz += depRHS[ix].dep().size();
    }

    if( _nnzjac < nnz ) _nnzjac = nnz;
#ifdef CRONOS__DEBUG__BASE_DE
    std::cout << "NNZ: " << nnz << "   MAX: " << _nnzjac << std::endl;
#endif
  }

  return true;
}

// 2026-09-28: over-alignment guard (GCC wrong-code with over-aligned virtual bases)
static_assert( alignof( ODESLV_BASE ) <= alignof( void* ), "ODESLV_BASE is a VIRTUAL base of the CRONOS solvers and must not be over-aligned (alignof > 8): GCC (6 to at least 16) emits aligned vector stores in base-object constructors assuming the full alignment, but a virtual-base subobject is only placed at its non-virtual alignment -> SIGSEGV at -O2/-O3 (see gccvb/pr_vbase_align.cpp). Keep over-aligned members (Armadillo/Eigen fixed-size types, alignas) behind a pointer, as FFModel::_pClassification does." );

} // end namescape mc

#endif

