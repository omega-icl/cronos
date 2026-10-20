// Copyright (C) Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the EPL-2.0 with GPL-2.0-or-later as a Secondary License; see the LICENSE file.
// SPDX-License-Identifier: EPL-2.0 OR GPL-2.0-or-later

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
#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <list>

#include "ffmodel.hpp"
#include "base_dynamic.hpp"

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
: public virtual BASE_DYNAMIC
, public virtual FFModel
{
 protected:
  //! @brief ODE policy for one evolution row read as L * dx/dt + R: REFUSE a vanishing coefficient
  //! (that row is algebraic and the model belongs to IDAS), and DIVIDE so the integrator gets an
  //! explicit right-hand side.  Byte-for-byte what _extract_from_model did inline before the hook.
  bool _read_evolution_row
    ( FFVar const& L, FFVar const& R, size_t i, std::vector<FFVar>& lhs, std::vector<FFVar>& de )
    override
    { if( L.cst() && L.num().val() == 0. ) return false;
      lhs[i] = L;
      de[i]  = ( L.cst() && L.num().val() == 1. )? _in_model_dag( -R ): _in_model_dag( -R / L );
      return true; }

  //! @brief 4.10 (2026-10-09): a POINT-FORM output, add_output( f, {t}, {tau}, {side} ), becomes the evaluation
  //! FFEval( f, t, tau, side ) before the evolution-direction reduction, so it is served like the expression form
  //! add_output( OpE( f, t, tau, side ) ).  Until now ODESLV ignored an output's point: it read f at the FINAL time
  //! (a state: silently the wrong value) or failed on a time-varying input ("missing variable").
  bool _on_before_reduction
    ()
    override
    {
      if( !_evolution_dom_var.dag() ) return true;
      FFEval OpE;
      for( auto& fct : _mFct ){
        auto const it = fct.point.find( _evolution_dom_var );
        if( it == fct.point.end() ) continue;
        if( fct.point.size() != 1 ){
          _extractError = "a point output fixes a coordinate other than the evolution direction (" + fct.var.name()
                        + "): ODESLV has no other direction";
          return false;
        }
        auto const is = fct.side.find( _evolution_dom_var );
        int const side = ( is != fct.side.end() && is->second == FFDom::PLUS )? (int)FFDom::PLUS: (int)FFDom::MINUS;
        fct.var = OpE( fct.var, _evolution_dom_var, it->second, side );
        fct.point.clear();  fct.side.clear();
      }
      return true;
    }

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

  //! @brief ODESLV substitutes fixed-input values at extraction (setup), not per call.
  bool _fixed_inputs_applied_per_call
    () const override
    { return false; }

  //! @brief Revision of the ODESLV headers (odeslv_base.hpp, odeslvs_base.hpp, odeslv_cvodes.hpp, odeslvs_cvodes.hpp,
  //! which change together): bump it with any change to one of them.  Since 2026-10-07, so that a sweep log can show
  //! which ODESLV headers it ran (ODESLV_revision prints it).
  static constexpr char const* HEADER_ID
    = "odeslv  rev1  2026-10-07";

  //! @brief The revision of the ODESLV headers (HEADER_ID), e.g. for a bug report
  static char const* revision() { return HEADER_ID; }

protected:
  // ---- what the old base class provided, now this class's own --------------------------------------------
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

  double _t = 0.;        //!< current time
  size_t _istg = 0;      //!< current stage
  size_t _nnzjac = 0;    //!< nonzeros in the Jacobian





  //! @brief local copy of DAG
  FFGraph* _dag;

  //! @brief local copy of initial value functions
  //! @brief Sensitivities of the consistent initial state, d x_i(0) / d p_j, from the SHARED system --
  //! filled only where the initial condition was SOLVED rather than divided.  Empty on the explicit
  //! fast path, where the divided form IS the initial condition and differentiating it is exact.
  std::vector<std::vector<double>> _dICSens, _dICSensDot;

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
  // and give the shared IC system its own copy: it is evaluated in _icDag, which these .set() values
  // do not reach
  _dIcCst.assign( _nc, 0. );
  for( unsigned ic=0; c && ic<_nc; ++ic ) _dIcCst[ic] = c[ic];

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
  // The DIVIDED form is exact only where every INITIAL row pins one state linearly with no derivative,
  // and the bypass is a KNOB (INITIAL.BYPASS_EXPLICIT) rather than a fixed rule, so the shared solve can
  // be forced on an explicit model and the two routes compared on one model.
  // Otherwise ic[i] carries other states, or a derivative, and cannot be evaluated at all -- so the rows
  // go to the SHARED consistent-initialisation solve in BASE_DYNAMIC instead, which poses them jointly
  // with the dynamics residual at t0 over (x, dx/dt).  The explicit case stays on the fast path: a Newton
  // solve buys nothing there, and it is the common one.
  if( !( _ic_bypass_explicit() && explicit_initial_conditions() ) && !ic_rows().empty() ){
    std::vector<double> xv, dxv, par( _np, 0. );
    // _Dp, not an offset of my own: _INI_D_STA lays _DVAR out as [ x, p, t, q ], so the parameters
    // start at _nx and NOT at _nx+1.  Reading them one slot late gave the IC solve zeros for every
    // parameter -- a wrong ANSWER, not a failure.  Unreachable from the corpus until now, since the
    // only driver with implicit initial rows declares no parameters.
    for( size_t j = 0; j < _np; ++j ) par[j] = _Dp[j];
    BASE_KINSOL::Report const R = BASE_DYNAMIC::solve_ic( xv, dxv, _ic_knobs(), par );
    if( R.status == BASE_KINSOL::FAILED ) return false;
    // A rank-deficient root is a DIAGNOSIS, not a failure: the iterate is accurate only to sqrt(FTOL)
    // there, so say which it is rather than returning it as if it were converged, exactly as DAESLV does
    // at its own stops.  Dropping it here was an asymmetry between the two consumers of one solve.
    if( R.status == BASE_KINSOL::RANK_DEFICIENT )
      std::cerr << "ODESLV ** at the initial time: " << R.message << std::endl;
    for( size_t i = 0; i < _nx; ++i ) ((double*)x)[i] = xv[i];
    // and the SENSITIVITY initial conditions from the same system, by the linear solve against the
    // Jacobian just factorised.  The divided form _vIC does not represent these rows, so differentiating
    // it -- which is what _IC_SET_FSA does -- gives a wrong answer here rather than a refusal.
    if( !BASE_DYNAMIC::ic_sens( _dICSens, _dICSensDot ) ){ _dICSens.clear(); _dICSensDot.clear(); }
    return true;
  }
  _dICSens.clear();  _dICSensDot.clear();       // explicit: the divided form is exact, leave FSA alone
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

