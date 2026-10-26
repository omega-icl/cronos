// Copyright (C) Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the EPL-2.0 with GPL-2.0-or-later as a Secondary License; see the LICENSE file.
// SPDX-License-Identifier: EPL-2.0 OR GPL-2.0-or-later

#ifndef CRONOS__ODESLV_BASE_HPP
#define CRONOS__ODESLV_BASE_HPP

#undef  CRONOS__ODESLV_BASE_DEBUG

#include <chrono>
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
#include "base_de.hpp"

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
: public virtual BASE_DE
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


  //! @brief Number of parameters (was BASE_DE::np)


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


  //! @brief Block default compiler methods
  ODESLV_BASE( ODESLV_BASE const& ) = delete;
  ODESLV_BASE& operator=( ODESLV_BASE const& ) = delete;

  // ------------------------------------------------------------------------------------------
  // SENSITIVITY ANALYSIS.  Was class ODESLV_BASE until 2026-10-21.  It existed only to extend
  // this one; 33 of 34 drivers named the derived class, and the split cost a second set of
  // evaluation arrays that SHADOWED these -- renamed _pVARS, _DVARS, ... before this merge.
  // ------------------------------------------------------------------------------------------
public:
  /** @ingroup ODESLV
   *  @{
   */
  //! @brief Default constructor

  //! @brief Virtual destructor
  /** @} */

protected:

  //! @brief The sensitivity direction _IC_SET_FSA last prepared, for _IC_D_SEN to index ic_sens() by.
  unsigned _isenIC = 0;


  //using ODESLV_BASE::_pIC;

  //! @brief size of sensitivity variables
  unsigned _ny;

  //! @brief local copy of sensitivity variables
  FFVar* _pY;
  
  //! @brief size of sensitivity quadrature variables
  unsigned _nyq;

  //! @brief local copy of sensitivity quadrature variables
  FFVar* _pYQ;

  //! @brief array of subgraphs of sensitivity/adjoint RHS functions
  std::vector<FFSubgraph> _opSARHS;

  //! @brief subgraph of sensitivity/adjoint RHS Jacobian 
  FFSubgraph _opSAJAC;
  //! @brief The ADJOINT Jacobian pattern, OWNED here (2026-09-28): entry e is (row, column, node) of the adjoint
  //! matrix -- the forward entry transposed, its node negated when neg -- in CSR order, with its row pointers.  The
  //! FORWARD _pJAC / _pJACCOLNDX are never modified (they used to be swapped in place: see _RHS_SET_ASA).
  std::vector<unsigned> _JACB_row, _JACB_col;
  std::vector<FFVar>    _JACB_node;
  std::vector<size_t>   _JACB_rowptr;
  //! @brief Scratch for the dense ADJOINT Jacobian values -- its own buffer, NOT ODESLV_BASE::_DJAC, which the
  //! forward Jacobian (_JAC_D_STA, inherited) writes: the two used to share it (2026-09-28 audit)
  std::vector<double>   _DJACB;

  //! @brief array of subgraphs of sensitivity/adjoint quadrature functions
  std::vector<FFSubgraph> _opSAQUAD;

  //! @brief const pointer to RHS function in current stage of ODE system
  FFVar const* _pRHSS;

  //! @brief const pointer to quadrature integrand in current stage of ODE system
  FFVar const* _pQUADS;

  //! @brief vector of const pointers to sensitivity/adjoint RHS function in current stage of ODE system
  std::vector<FFVar const*> _vSARHS;

  //! @brief vector of const pointers to sensitivity/adjoint quadrature integrand in current stage of ODE system
  std::vector<FFVar const*> _vSAQUAD;

  //! @brief const pointer to adjoint TC function in current stage of ODE system
  FFVar const* _pSACFCT;

  //! @brief pointer to adjoint discontinuity function in current stage of ODE system
  FFVar* _pSAFCT;

  //! @brief number of variables for DAG evaluation
  unsigned _nVARS;

  //! @brief number of variables for DAG evaluation (without quadratures)
  unsigned _nVAR0S;

  //! @brief pointer to variables for DAG evaluation
  FFVar* _pVARS;

  //! @brief pointer of variable values for DAG evaluation
  double* _DVARS;

  //! @brief pointer to time **DO NOT FREE**
  double* _DtS;

  //! @brief pointer to sensitivity/adjoint values **DO NOT FREE**
  double* _Dy;

  //! @brief pointer to state values **DO NOT FREE**
  double* _DxS;

  //! @brief pointer to parameter values **DO NOT FREE**
  double* _DpS;

  //! @brief pointer to quadrature values **DO NOT FREE**
  double* _DqS;

  //! @brief pointer to sensitivity/adjoint quadrature values **DO NOT FREE**
  double* _Dyq;

  //! @brief vector to hold function derivatives
  std::vector<double> _Dfp;

  //! @brief storage vector for sensitivity/adjoint DAG evaluation
  std::vector<double> _DWRKS;

  //! @brief Function setting up local DAG
  bool _SETUP_SEN
    ();

  //! @brief Function setting up DAG of IVP
  bool _SETUP_SEN
    ( ODESLV_BASE const& IVP );

  //! @brief Function to initialize sensitivity for parameter <a>isen</a>
  bool _IC_SET_FSA
    ( unsigned const isen );

  //! @brief Function to add initial state contribution to adjoint quadrature values
  bool _IC_SET_ASA
    ( unsigned const ifct ); // ();

  //! @brief Function to reinitialize sensitivity at stage times for parameter <a>isen</a>
  bool _CC_SET_FSA
    ( unsigned const pos_ic, unsigned const isen );

  //! @brief Function to reinitialize adjoint at stage times for function <a>ifct</a>
  bool _CC_SET_ASA
    ( unsigned const pos_ic, unsigned const pos_fct, unsigned const ifct );

  //! @brief Function to initial adjoint at terminal time for function <a>ifct</a>
  bool _TC_SET_ASA
    ( unsigned const pos_fct, unsigned const ifct );

  //! @brief Function to set sensitivity RHS pointer
  bool _RHS_SET_FSA
    ( unsigned const iRHS, unsigned const iQUAD );

  //! @brief Function to set adjoint RHS pointer
  bool _RHS_SET_ASA
    ( unsigned const iRHS, unsigned const iQUAD,
      unsigned const pos_fct, bool const neg=true );

  //! @brief Function to initialize sensitivity/adjoint values
  bool _INI_D_SEN
    ( double const* p, unsigned const nf, unsigned const nyq );

  //! @brief Function to retreive sensitivity/adjoint values
  template <typename REALTYPE>
  void _GET_D_SEN
    ( REALTYPE const* y, unsigned const nyq, REALTYPE const* yq );

  //! @brief Function to retreive state and sensitivity/adjoint values
  template <typename REALTYPE>
  void _GET_D_SEN
    ( REALTYPE const* x, REALTYPE const* y, REALTYPE const* q,
      unsigned const nyq, REALTYPE const* yq );

  //! @brief Function to initialize sensitivity/adjoint values at terminal time
  template <typename REALTYPE>
  bool _TC_D_SEN
    ( double const& t, REALTYPE const* x, REALTYPE* y );

  //! @brief Function to initialize adjoint quadrature values at terminal time
  template <typename REALTYPE>
  bool _TC_D_QUAD_ASA
   ( REALTYPE* yq ); 

//  //! @brief Function to add initial state contribution to sensitivity values
  template <typename REALTYPE>
  bool _IC_D_SEN
    ( double const& t, REALTYPE* y );

  //! @brief Function to add initial state contributionvto adjoint values
  template <typename REALTYPE>
  bool _IC_D_SEN
    ( double const& t, REALTYPE const* x, REALTYPE const* y );

  //! @brief Function to add initial state contribution to adjoint quadrature values
  template <typename REALTYPE>
  bool _IC_D_QUAD_ASA
    ( REALTYPE* yq );

  //! @brief Function to set sensitivity/adjoint transitions at stage times
  bool _CC_D_SET
    ();

  //! @brief Function to transition sensitivity/adjoint values at stage times
  template <typename REALTYPE>
  bool _CC_D_SEN
    ( double const& t, REALTYPE const* x, REALTYPE* y );

  //! @brief Function to transition adjoint quadrature values at stage times
  template <typename REALTYPE>
  bool _CC_D_QUAD_ASA
    ( REALTYPE* yq );

  //! @brief Function to set sensitivity/adjoint RHS pointer
  bool _RHS_D_SET_SEN
    ( unsigned const nf, unsigned const nyq );

  //! @brief Function to calculate the RHS of sensitivity/adjoint ODEs
  template <typename REALTYPE>
  bool _RHS_D_SEN
  ( double const& t, REALTYPE const* x, REALTYPE const* y, REALTYPE* ydot,
    unsigned const ifct );

  //! @brief Function to calculate the Jacobian RHS of sensitivity/adjoint ODEs
  template <typename REALTYPE>
  bool _JAC_D_SEN
    ( double const& t, REALTYPE const* x, REALTYPE const* y, REALTYPE** jac );

  //! @brief Function to calculate the Jacobian RHS of sensitivity/adjoint ODEs
  template <typename REALTYPE, typename INDEXTYPE>
  bool _JAC_D_SEN
    ( double const& t, REALTYPE const* x, REALTYPE const* y, REALTYPE* jac,
      INDEXTYPE* ptrs, INDEXTYPE* vals );

  //! @brief Function to calculate the RHS of sensitivity/adjoint quadrature ODEs
  template <typename REALTYPE>
  bool _RHS_D_QUAD
    ( unsigned const nyq, REALTYPE* qdot, unsigned const ifct );

  //! @brief Function to calculate the function sensitivities at intermediate/end point
  bool _FCT_D_SEN
    ( unsigned const pos_fct, unsigned const isen, double const& t );

  //! @brief Block default compiler methods

};

inline 
ODESLV_BASE::ODESLV_BASE
()
: _dag(nullptr),
  _pC(nullptr), _pP(nullptr), _pT(nullptr), _pX(nullptr), _pQ(nullptr),
  _pJAC(0,nullptr,nullptr,nullptr),
  _nVAR(0), _nVAR0(0), _pVAR(nullptr),
  _DVAR(nullptr), _Dt(nullptr), _Dx(nullptr), _Dp(nullptr), _Dq(nullptr),
  _ny(0), _pY(nullptr), _nyq(0), _pYQ(nullptr),
  _pRHSS(nullptr), _pQUADS(nullptr), _pSACFCT(nullptr), _pSAFCT(nullptr),
  _nVARS(0), _nVAR0S(0), _pVARS(nullptr)
{
  _DVARS = _DtS = _DpS = _Dy = _DxS = _DqS = _Dyq = nullptr;
}

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

  // merged from ~ODESLV_BASE, 2026-10-21
  delete[] _pVARS;
  delete[] _DVARS;
  /* DO NOT FREE _pRHSS, _pQUADS */
  for( auto& rhs  : _vSARHS  ) delete[] rhs;
  for( auto& quad : _vSAQUAD ) delete[] quad;
  delete[] _pSACFCT;
  delete[] _pSAFCT;
  delete[] _pY;
  delete[] _pYQ;
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
  // go to the SHARED consistent-initialisation solve in BASE_DE instead, which poses them jointly
  // with the dynamics residual at t0 over (x, dx/dt).  The explicit case stays on the fast path: a Newton
  // solve buys nothing there, and it is the common one.
  if( !( _ic_bypass_explicit() && explicit_initial_conditions() ) && !ic_rows().empty() ){
    std::vector<double> xv, dxv, par( _np, 0. );
    // _Dp, not an offset of my own: _INI_D_STA lays _DVAR out as [ x, p, t, q ], so the parameters
    // start at _nx and NOT at _nx+1.  Reading them one slot late gave the IC solve zeros for every
    // parameter -- a wrong ANSWER, not a failure.  Unreachable from the corpus until now, since the
    // only driver with implicit initial rows declares no parameters.
    for( size_t j = 0; j < _np; ++j ) par[j] = _Dp[j];
    BASE_KINSOL::Report const R = BASE_DE::solve_ic( xv, dxv, _ic_knobs(), par );
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
    if( !BASE_DE::ic_sens( _dICSens, _dICSensDot ) ){ _dICSens.clear(); _dICSensDot.clear(); }
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
  // Row by row, each row against only the states it depends on (BASE_DE::_sfad_by_row): the one-call
  // SFAD is O(nx^2) however sparse the result, and this runs at EVERY solve -- 1.8 s at 1600 states,
  // against a 0.44 s integration (2026-10-25).  Column-ordered, as the pointer loop below requires.
  _pJAC = _sfad_by_row( _dag, std::vector<FFVar>( _pRHS, _pRHS + _nx ), std::vector<FFVar>( _pX, _pX + _nx ) );
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


inline
bool
ODESLV_BASE::_SETUP_SEN
()
{
  if( !ODESLV_BASE::_SETUP() ) return false;
  
  delete[] _pY; _pY = nullptr;
  _ny = _nx;
  if( _ny ){
    _pY  = odeslv_new_vars( _ny );
    for( unsigned iy=0; iy<_ny; ++iy ) _pY[iy].set( _dag );
  }
//  std::cerr << "_pY[" << _ny << "]: " << _pY << std::endl;

  delete[] _pYQ; _pYQ = nullptr;
  _nyq = ( _nsen<_nq? _nq: _nsen );
  if( _nyq ){
    _pYQ  = odeslv_new_vars( _nyq );
    for( unsigned iyq=0; iyq<_nyq; ++iyq ) _pYQ[iyq].set( _dag );
  }
//  std::cerr << "_pYQ[" << _nyq << "]: " << _pYQ << std::endl;
  
  return true;
}

inline
bool
ODESLV_BASE::_SETUP_SEN
( ODESLV_BASE const& IVP )
{
  if( !ODESLV_BASE::_SETUP( IVP) ) return false;

  delete[] _pY; _pY = nullptr;
  _ny = _nx;
  if( _ny ){
    _pY  = odeslv_new_vars( _ny );
    for( unsigned iy=0; iy<_ny; ++iy ) _pY[iy].set( _dag );
  }
//  std::cerr << "_pY[" << _ny << "]: " << _pY << std::endl;

  delete[] _pYQ; _pYQ = nullptr;
  _nyq = ( _nsen<_nq? _nq: _nsen );
  if( _nyq ){
    _pYQ  = odeslv_new_vars( _nyq );
    for( unsigned iyq=0; iyq<_nyq; ++iyq ) _pYQ[iyq].set( _dag );
  }
//  std::cerr << "_pYQ[" << _nyq << "]: " << _pYQ << std::endl;
  
  return true;
}

inline
bool
ODESLV_BASE::_INI_D_SEN
( double const* p, unsigned const nvec, unsigned const nyq )
{
  // Size and set DAG evaluation arrays
  if( nyq > _nyq ) return false;
  _nVAR0S = _ny + _nx + _np + 1;
  _nVARS  = _nVAR0S + _nq + nyq;
  delete[] _pVARS; _pVARS = odeslv_new_vars( _nVARS );
  delete[] _DVARS; _DVARS = new double[_nVARS];
  delete[] _pSAFCT;  _pSAFCT  = odeslv_new_vars( _ny+_np+1+_nq );
  for( auto it=_vSARHS.begin(); it!=_vSARHS.end(); ++it ){ delete[] *it; *it=0; }
  for( auto it=_vSAQUAD.begin(); it!=_vSAQUAD.end(); ++it ){ delete[] *it; *it=0; }
  _vSARHS.resize(nvec); _vSAQUAD.resize(nvec);

  for( unsigned ix=0; ix<_nx; ix++ ) _pVARS[ix] = _pY[ix];
  for( unsigned ix=0; ix<_nx; ix++ ) _pVARS[_ny+ix] = _pX[ix];
  for( unsigned ip=0; ip<_np; ip++ ) _pVARS[_ny+_nx+ip] = _pP[ip];
  _pVARS[2*_nx+_np] = (_pT? *_pT: 0. );
  for( unsigned iq=0; iq<_nq; iq++ ) _pVARS[_ny+_nx+_np+1+iq] = _pQ?_pQ[iq]:0.;
  for( unsigned iyq=0; iyq<nyq; iyq++ ) _pVARS[_ny+_nx+_np+1+_nq+iyq] = _pYQ?_pYQ[iyq]:0.;
  _Dy = _DVARS;
  _DxS = _Dy + _ny;
  _DpS = _DxS + _nx;
  _DtS = _DpS + _np;
  _DqS = _DtS + 1;
  _Dyq = _DqS + _nq;
  for( unsigned ip=0; ip<_np; ip++ ) _DpS[ip] = p[ip];

  _Dfp.assign( _nf*_nsen, 0. );

  return true;
}

template <typename REALTYPE>
inline
void
ODESLV_BASE::_GET_D_SEN
( REALTYPE const* y, unsigned const nyq, REALTYPE const* yq )
{
  _vec2D( y, _ny, _Dy );
  if( yq ) _vec2D( yq, nyq, _Dyq );
}

template <typename REALTYPE>
inline
void
ODESLV_BASE::_GET_D_SEN
( REALTYPE const* x, REALTYPE const* y, REALTYPE const* q,
  unsigned const nyq, REALTYPE const* yq )
{
  _vec2D( x, _nx, _DxS );
  _vec2D( y, _ny, _Dy );
  if( q )  _vec2D( q, _nq, _DqS );
  if( yq ) _vec2D( yq, nyq, _Dyq );
}

inline
bool
ODESLV_BASE::_IC_SET_FSA
( unsigned const isen )
{
  FFVar const* pIC = _vIC.at( 0 ).data();
#ifdef CRONOS__ODESLVS_BASE_DEBUG
  auto sgIC = _dag->subgraph( _nx, pIC );
  //_dag->output( sgIC, " pIC" );
  std::vector<FFExpr> exIC = FFExpr::subgraph( _dag, sgIC ); 
  for( size_t i=0; i<_nx; ++i )
    std::cout << "IC[" << i << "] = " << exIC[i] << std::endl;
  std::cout << "P[" << _ndxSEN[isen] << "] = " << _pVARS[_ny+_nx+_ndxSEN[isen]] << std::endl;
  {std::cout << "Enter <1> to continue"; int dum; std::cin >> dum;}
#endif

  _isenIC = isen;                               // which direction _IC_D_SEN is about to be asked for
  delete[] _pSACFCT; _pSACFCT = _dag->FAD( _nx, pIC, 1, _pVARS+_ny+_nx+_ndxSEN[isen] );
  pIC = _pSACFCT;
  return true;
}
/*
inline
bool
ODESLV_BASE::_IC_SET_ASA
()
{
  FFVar const* pIC = _vIC.at( 0 );
  FFVar pHAM( 0. );
  for( unsigned ix=0; ix<_nx; ix++ ) pHAM += _pY[ix] * pIC[ix];
#ifndef CRONOS__ODESLVS_USE_BAD
  std::vector<FFVar> Psen( _nsen );   // the selected parameters, in control order
  for( unsigned is=0; is<_nsen; ++is ) Psen[is] = _pP[_ndxSEN[is]];
  delete[] _pSACFCT; _pSACFCT = _dag->FAD( 1, &pHAM, _nsen, Psen.data() );
#else
  std::vector<FFVar> Psen( _nsen );
  for( unsigned is=0; is<_nsen; ++is ) Psen[is] = _pP[_ndxSEN[is]];
  delete[] _pSACFCT; _pSACFCT = _dag->BAD( 1, &pHAM, _nsen, Psen.data() );
#endif
  return true;
}
*/
inline
bool
ODESLV_BASE::_IC_SET_ASA
( unsigned const ifct )
{
  FFVar const* pIC = _vIC.at( 0 ).data();
  FFVar pHAM = (!_vFCT.at( 0 ).empty()? _vFCT.at( 0 )[ifct]: 0.);
  for( unsigned ix=0; ix<_nx; ix++ ) pHAM += _pY[ix] * pIC[ix];
  //_dag->output( _dag->subgraph( 1, &pHAM ) );
#ifndef CRONOS__ODESLVS_USE_BAD
  delete[] _pSACFCT; _pSACFCT = _dag->FAD( 1, &pHAM, _nx+_np, _pVARS+_ny );
#else
  delete[] _pSACFCT; _pSACFCT = _dag->BAD( 1, &pHAM, _nx+_np, _pVARS+_ny );
#endif
  bool df0dx_nonzero = false;
  for( unsigned ix=0; ix<_nx && !df0dx_nonzero; ix++ )
    if( !( _pSACFCT[ix].cst() && _pSACFCT[ix].num().val() == 0. ) ) df0dx_nonzero = true;
  for( unsigned is=0; is<_nsen; is++ ){
    _pSAFCT[is] = _pSACFCT[_nx+_ndxSEN[is]];          // d( f0 + y.IC )/dp at a fixed initial state
    // + the chain through the initial state, (df0/dx) . d(x0)/dp: d(x0)/dp is the derivative of the INITIAL-CONDITION
    // expression w.r.t. the parameter.  It used to multiply by the initial VALUE x0 instead, so an output at the
    // initial time had a wrong adjoint gradient (1.0 for x0 = 1, whatever the parameter) -- WORKPLAN 1.8, 2026-10-06.
    // Only a function at stage 0 makes df0/dx non-zero: skip the derivative of the initial condition otherwise.
    if( df0dx_nonzero ){
      FFVar* dIC = _dag->FAD( _nx, pIC, 1, _pP + _ndxSEN[is] );
      for( unsigned ix=0; ix<_nx; ix++ ) _pSAFCT[is] += _pSACFCT[ix] * dIC[ix];
      delete[] dIC;
    }
  }

  return true;
}

inline
bool
ODESLV_BASE::_CC_SET_FSA
( unsigned const pos_ic, unsigned const isen )
{
  FFVar const* pIC = _vIC.at( pos_ic ).data();
  for( unsigned iy=0; iy<_ny; iy++ ) _pSAFCT[iy] = _pVARS[iy];
  for( unsigned ip=0; ip<_np; ip++ ) _pSAFCT[_ny+ip] = (ip==_ndxSEN[isen]? 1.: 0.);
  delete[] _pSACFCT; _pSACFCT = _dag->DFAD( _nx, pIC, _nx+_np, _pVARS+_ny, _pSAFCT );
  for( unsigned iy=0; iy<_ny; iy++ ) _pSAFCT[iy] = _pSACFCT[iy];

  return true;
}

inline
bool
ODESLV_BASE::_CC_SET_ASA
( unsigned const pos_ic, unsigned const pos_fct, unsigned const ifct )
{
  FFVar const* pIC = _vIC.at( pos_ic ).data();
  FFVar pHAM = (!_vFCT.at(pos_fct).empty()? _vFCT.at(pos_fct)[ifct]: 0.);
  //FFVar pHAM = (pos_fct? _vFCT.at(pos_fct-1)[ifct]: 0.);
#ifdef CRONOS__ODESLVS_BASE_DEBUG
  std::cout << "pos_ic: " << pos_ic << std::endl;
  std::cout << "pos_fct: " << pos_fct << std::endl;
#endif
  for( unsigned ix=0; ix<_nx; ix++ ) pHAM += _pY[ix] * (pos_ic? pIC[ix]: _pX[ix] );
  //_dag->output( _dag->subgraph( 1, &pHAM ) );
  // w.r.t. the states, then the SELECTED parameters in control order -- as _IC_SET_ASA does, and as _CC_D_QUAD_ASA
  // reads it (_nsen entries after the _ny state ones).  2026-09-29: this used ALL _np parameters, so with a registered
  // subset of controls the stage-transition increments landed on the wrong gradient entries.
  std::vector<FFVar> W( _pVARS+_ny, _pVARS+_ny+_nx );
  for( unsigned is=0; is<_nsen; ++is ) W.push_back( _pP[_ndxSEN[is]] );
#ifndef CRONOS__ODESLVS_USE_BAD
  delete[] _pSACFCT; _pSACFCT = _dag->FAD( 1, &pHAM, W.size(), W.data() );
#else
  delete[] _pSACFCT; _pSACFCT = _dag->BAD( 1, &pHAM, W.size(), W.data() );
#endif
  for( unsigned iy=0; iy<_nx; iy++ )
    _pSAFCT[iy] = _pSACFCT[iy];

  return true;
}

inline
bool
ODESLV_BASE::_TC_SET_ASA
( unsigned const pos_fct, unsigned const ifct )
{
  FFVar const* pFCT = _vFCT.at(pos_fct).data();
  //FFVar const* _pIC = _vFCT.at(pos_fct);//+ifct;

#ifndef CRONOS__ODESLVS_USE_BAD
  delete[] _pSACFCT; _pSACFCT = _dag->FAD( 1, pFCT+ifct, _nx+_np, _pVARS+_nx );
#else
  delete[] _pSACFCT; _pSACFCT = _dag->BAD( 1, pFCT+ifct, _nx+_np, _pVARS+_nx );
#endif

  return true;
}

inline
bool
ODESLV_BASE::_RHS_SET_FSA
( unsigned const iRHS, unsigned const iQUAD )
{
  if( _vRHS.size() <= iRHS ) return false; 
  if( _nq && _vQUAD.size() <= iQUAD ) return false;

  _pRHSS =  _vRHS.at( iRHS ).data();
#ifdef CRONOS__ODESLVS_BASE_DEBUG
  std::ostringstream ofilename;
  ofilename << "vRHS.dot";
  std::ofstream ofile( ofilename.str(), std::ios_base::out );
  _dag->dot_script( _nx, _pRHSS, ofile );
  ofile.close();
#endif
  _pQUADS = _nq? _vQUAD.at( iQUAD ).data(): nullptr;

  // Set sensitivity ODEs using directional derivatives
  for( unsigned iy=0; iy<_nx; iy++ ) _pSAFCT[iy] = _pVARS[iy];
  for( unsigned is=0; is<_nsen; is++ ){
    unsigned const ip = _ndxSEN[is];
    for( unsigned jp=0; jp<_np; jp++ ) _pSAFCT[_nx+jp] = (ip==jp? 1.: 0.);
    delete[] _vSARHS[is];  _vSARHS[is]  = _dag->DFAD( _nx, _pRHSS, _nx+_np, _pVARS+_nx, _pSAFCT );
#ifdef CRONOS__ODESLVS_BASE_DEBUG
    std::ostringstream ofilename;
    ofilename << "vSARHS" << is << ".dot";
    std::ofstream ofile( ofilename.str(), std::ios_base::out );
    _dag->dot_script( _nx, _vSARHS[ip], ofile );
    ofile.close();
#endif
    if( !_nq ) continue;
    delete[] _vSAQUAD[is]; _vSAQUAD[is] = _dag->DFAD( _nq, _pQUADS, _nx+_np, _pVARS+_nx, _pSAFCT );

  }

  return true;
}

inline
bool
ODESLV_BASE::_RHS_SET_ASA
( unsigned const iRHS, unsigned const iQUAD,
  unsigned const iFCT, bool const neg )
{
  if( _vRHS.size() <= iRHS ) return false; 
  if( _nq && _vQUAD.size() <= iQUAD ) return false;

  _pRHSS = _vRHS.at( iRHS ).data();
  FFVar pHAM( 0. );
  for( unsigned ix=0; ix<_nx; ix++ ){
    if( neg ) pHAM -= _pY[ix] * _pRHSS[ix];
    else      pHAM += _pY[ix] * _pRHSS[ix];
  }
  _pQUADS  = _nq? _vQUAD.at( iQUAD ).data(): nullptr;
  std::vector<FFVar> vHAM( _nf, pHAM );
  
  if( _nq && !_vFCT.at(iFCT).empty() ){
    FFVar const* pFCT = _vFCT.at(iFCT).data();
    //FFVar const* pFCT = _vFCT.at(pos_fct).data();
    for( unsigned ifct=0; ifct<_nf; ifct++ ){
#ifndef CRONOS__ODESLVS_USE_BAD
      delete[] _pSACFCT; _pSACFCT = _dag->FAD( 1, pFCT+ifct, _nq, _pQ );
#else
      delete[] _pSACFCT; _pSACFCT = _dag->BAD( 1, pFCT+ifct, _nq, _pQ );
#endif
      for( unsigned iq=0; iq<_nq; iq++ ){
        if( !_pSACFCT[iq].cst() ) return false; // quadrature appears nonlinearly in function
        if( neg ) vHAM[ifct] -= _pQUADS[iq] * _pSACFCT[iq];
        else      vHAM[ifct] += _pQUADS[iq] * _pSACFCT[iq];
      }
    }
  }
#ifdef CRONOS__ODESLVS_BASE_DEBUG
  _dag->output( _dag->subgraph( _nf, vHAM.data() ) );
#endif

  for( unsigned ifct=0; ifct<_nf; ifct++ ){
    std::vector<FFVar> Psen( _nsen );   // the selected parameters, in control order
    for( unsigned is=0; is<_nsen; ++is ) Psen[is] = _pP[_ndxSEN[is]];
#ifndef CRONOS__ODESLVS_USE_BAD
    delete[] _vSARHS[ifct];  _vSARHS[ifct]  = _dag->FAD( 1, vHAM.data()+ifct, _nx, _pX   );
    delete[] _vSAQUAD[ifct]; _vSAQUAD[ifct] = _dag->FAD( 1, vHAM.data()+ifct, _nsen, Psen.data() );
#else
    delete[] _vSARHS[ifct];  _vSARHS[ifct]  = _dag->BAD( 1, vHAM.data()+ifct, _nx, _pX   );
    delete[] _vSAQUAD[ifct]; _vSAQUAD[ifct] = _dag->BAD( 1, vHAM.data()+ifct, _nsen, Psen.data() );
#endif
  }

  // The ADJOINT Jacobian pattern, built in its OWN storage (2026-09-28).  This used to std::swap the row and column
  // arrays of the FORWARD pattern _pJAC in place, negate its nodes, and overwrite _pJACCOLNDX with adjoint row
  // pointers -- so every FORWARD Jacobian assembled during CVodeB (the checkpoint replays in CVAdataStore) came out
  // TRANSPOSED: a wrong Newton matrix, a replay that departed from the forward pass, and CVODES' dt_mem overrun.
  // (The dense adjoint writer then transposed back, so the dense ADJOINT matrix was J, not J^T.)
  // SFAD orders the forward entries COLUMNWISE, i.e. by the adjoint's ROW: CSR order of the transpose, which is
  // what _INI_ASA declares (CSR_MAT).
  unsigned const nnzJ = std::get<0>(_pJAC);
  _JACB_row.resize( nnzJ );  _JACB_col.resize( nnzJ );  _JACB_node.resize( nnzJ );
  _JACB_rowptr.assign( _nx+1, 0 );
  unsigned ir = 0;
  for( unsigned ie=0; ie<nnzJ; ++ie ){
    _JACB_row[ie]  = std::get<2>(_pJAC)[ie];     // forward column -> adjoint row
    _JACB_col[ie]  = std::get<1>(_pJAC)[ie];     // forward row    -> adjoint column
    _JACB_node[ie] = neg? -std::get<3>(_pJAC)[ie]: std::get<3>(_pJAC)[ie];
    for( ; _JACB_row[ie] >= ir; ++ir ) _JACB_rowptr[ir] = ie;
#ifdef CRONOS__ODESLVS_BASE_DEBUG
    std::cout << "  JACB[" << _JACB_row[ie] << ", " << _JACB_col[ie] << "]" << std::endl;
#endif
  }
  for( ; ir<=_nx; ++ir ) _JACB_rowptr[ir] = nnzJ;  // empty trailing rows, and the end pointer

  return true;
}

template <typename REALTYPE>
inline
bool
ODESLV_BASE::_TC_D_SEN
( double const& t, REALTYPE const* x, REALTYPE* y )
{
  if( !_pSACFCT ){
    for( unsigned iy=0; iy<_ny; ++iy )
      y[iy] = 0.;
    return true;
  }
  
  *_DtS = t; // current time
  _vec2D( x, _nx, _DxS );
  _dag->eval( _ny, _pSACFCT, (double*)y, _nx+_np+1, _pVARS+_ny, _DVARS+_ny );
  return true;
}

template <typename REALTYPE>
inline
bool
ODESLV_BASE::_TC_D_QUAD_ASA
( REALTYPE* yq )
{
  if( !_pSACFCT ){
    for( unsigned iy=0; iy<_nsen; ++iy )
      yq[iy] = 0.;
    return true;
  }

  _dag->eval( _nsen, _pSACFCT+_ny, (double*)yq, _nx+_np+1, _pVARS+_ny, _DVARS+_ny );
  return true;
}

template <typename REALTYPE>
inline
bool
ODESLV_BASE::_IC_D_SEN
( double const& t, REALTYPE* y )
{
#ifdef CRONOS__ODESLVS_BASE_DEBUG
  auto sgSACFCT = _dag->subgraph( _nx, _pSACFCT );
  //_dag->output( sgSACFCT ), " _pSACFCT" );
  std::vector<FFExpr> exSACFCT = FFExpr::subgraph( _dag, sgSACFCT ); 
  for( size_t i=0; i<_nx; ++i )
    std::cout << "SACFCT[" << i << "] = " << exSACFCT[i] << std::endl;
  for( size_t i=0; i<_np+1; ++i )
    std::cout << "P[" << i << "] = " << _pVARS[_ny+_nx+i] << std::endl;
  {std::cout << "Enter <1> to continue"; int dum; std::cin >> dum;}
#endif

  *_DtS = t; // current time
  // Where the initial condition was SOLVED, its sensitivities came from the same system (ic_sens) and
  // _pSACFCT -- the FAD of the divided form -- is not the derivative of it.  Prefer the solved pair.
  if( !_dICSens.empty() ){
    size_t const j = _ndxSEN[_isenIC];
    if( j < _dICSens.size() && _dICSens[j].size() >= _nx ){
      for( unsigned i = 0; i < _nx; ++i ) ((double*)y)[i] = _dICSens[j][i];
      return true;
    }
  }
  _dag->eval( _nx, _pSACFCT, (double*)y, _np+1, _pVARS+_ny+_nx, _DVARS+_ny+_nx );
  return true;
}

template <typename REALTYPE>
inline
bool
ODESLV_BASE::_IC_D_SEN
( double const& t, REALTYPE const* x, REALTYPE const* y )
{
  *_DtS = t; // current time
  _vec2D( x, _nx, _DxS );
  _vec2D( y, _ny, _Dy );
  return true;
}

template <typename REALTYPE>
inline
bool
ODESLV_BASE::_IC_D_QUAD_ASA
( REALTYPE* yq )
{
  _vec2D( yq, _nsen, _Dyq );
  static double const one = 1.;
  _dag->eval( _nsen, _pSAFCT, (double*)yq, _nVAR0S, _pVARS, _DVARS, &one );
  _vec2D( yq, _nsen, _Dyq );
  return true;
}

template <typename REALTYPE>
inline
bool
ODESLV_BASE::_CC_D_SEN
( double const& t, REALTYPE const* x, REALTYPE* y )
{
  *_DtS = t; // current time
  _vec2D( x, _nx, _DxS );
  _vec2D( y, _ny, _Dy );
  _dag->eval( _ny, _pSAFCT, (double*)y, _nVAR0S, _pVARS, _DVARS );
  return true;
}

template <typename REALTYPE>
inline
bool
ODESLV_BASE::_CC_D_QUAD_ASA
( REALTYPE* yq )
{
  _vec2D( yq, _nsen, _Dyq );
  static double const one = 1.;
  _dag->eval( _nsen, _pSACFCT+_ny, (double*)yq, _nVAR0S, _pVARS, _DVARS, &one );
  return true;
}

inline
bool
ODESLV_BASE::_RHS_D_SET_SEN
( unsigned const nf, unsigned const nyq )
{
  _opSARHS.resize( nf );
  for( unsigned ifct=0; ifct<nf; ifct++ )
    _opSARHS[ifct]  = _dag->subgraph( _ny, _vSARHS[ifct] );

  if( nyq ){
    _opSAQUAD.resize( nf );
    for( unsigned ifct=0; ifct<nf; ifct++ )
      _opSAQUAD[ifct] = _dag->subgraph( nyq, _vSAQUAD[ifct] );
  }

  _opSAJAC = _dag->subgraph( _JACB_node.size(), _JACB_node.data() );
  _DJACB.resize( _JACB_node.size() );

  return true;
}

template <typename REALTYPE>
inline
bool
ODESLV_BASE::_RHS_D_SEN
( double const& t, REALTYPE const* x, REALTYPE const* y, REALTYPE* ydot,
  unsigned const ifct )
{
  if( !_vSARHS.size() || !_vSARHS[ifct] ) return false; // **error** ADJRHS not defined
  *_DtS = t; // set current time
  _vec2D( x, _nx, _DxS ); // set current state bounds
  _vec2D( y, _ny, _Dy ); // set current sensitivity/adjoint bounds
#ifdef CRONOS__ODESLVS_BASE_DEBUG
  _dag->output( _opSARHS[ifct] );
#endif
  _dag->eval( _opSARHS[ifct], _DWRKS, _ny, _vSARHS[ifct], (double*)ydot,
               _nVAR0S, _pVARS, _DVARS );
  return true;
}

template <typename REALTYPE>
inline
bool
ODESLV_BASE::_JAC_D_SEN
( double const& t, REALTYPE const* x, REALTYPE const* y, REALTYPE** jac )
{
  // 2026-09-28: load the callback's OWN (t, x) -- relying on the last RHS call made the Jacobian / quadrature
  // stale whenever another evaluation came in between (the forward replay inside CVodeB: wrong J, replay overrun)
  *_DtS = t;  _vec2D( x, _nx, _DxS );  _vec2D( y, _ny, _Dy );
  _dag->eval( _opSAJAC, _DWRKS, _JACB_node.size(), _JACB_node.data(),
               _DJACB.data(), _nVAR0S, _pVARS, _DVARS );
  for( unsigned ie=0; ie<_JACB_node.size(); ++ie ){
    jac[_JACB_col[ie]][_JACB_row[ie]] = _DJACB[ie];     // SM_COLS_D: jac[column][row]
#ifdef CRONOS__ODESLVS_BASE_DEBUG
    std::cout << "  jacB[" << _JACB_row[ie] << ", "
              << _JACB_col[ie] << "] = " << _DJACB[ie] << std::endl;
#endif
  }
  return true;
}

template <typename REALTYPE, typename INDEXTYPE>
inline
bool
ODESLV_BASE::_JAC_D_SEN
( double const& t, REALTYPE const* x, REALTYPE const* y, REALTYPE* jac,
  INDEXTYPE* ptrs, INDEXTYPE* vals )
{
  // 2026-09-28: load the callback's OWN (t, x) -- relying on the last RHS call made the Jacobian / quadrature
  // stale whenever another evaluation came in between (the forward replay inside CVodeB: wrong J, replay overrun)
  *_DtS = t;  _vec2D( x, _nx, _DxS );  _vec2D( y, _ny, _Dy );
  _dag->eval( _opSAJAC, _DWRKS, _JACB_node.size(), _JACB_node.data(),
               (double*)jac, _nVAR0S, _pVARS, _DVARS );
  for( unsigned ie=0; ie<_JACB_node.size(); ++ie ){
    vals[ie] = (INDEXTYPE)_JACB_col[ie];             // CSR index values are COLUMNS
#ifdef CRONOS__ODESLVS_BASE_DEBUG
    std::cout << "  jac[" << ie << "] = " << jac[ie] << std::endl;
    std::cout << "  vals[" << ie << "] = " << vals[ie] << std::endl;
#endif
  }
  for( unsigned ic=0; ic<=_nx; ++ic ){
    ptrs[ic] = (INDEXTYPE) _JACB_rowptr[ic];         // CSR ROW pointers of the adjoint
#ifdef CRONOS__ODESLVS_BASE_DEBUG
    std::cout << "  ptrs[" << ic << "] = " << ptrs[ic] << std::endl;
#endif
  }
  return true;
}

template <typename REALTYPE>
inline
bool
ODESLV_BASE::_RHS_D_QUAD
( unsigned const nyq, REALTYPE* yqdot, unsigned const ifct )
{
  if( !_vSAQUAD.size() || !_vSAQUAD[ifct] ) return false;
  // _DVARS (t, x, y) is loaded by the CALLER, which has them (CVQUADB__, the sensitivity-quadrature callback)
  _dag->eval( _opSAQUAD[ifct], _DWRKS, nyq, _vSAQUAD[ifct], (double*)yqdot,
               _nVAR0S, _pVARS, _DVARS );
  return true;
}

inline
bool
ODESLV_BASE::_FCT_D_SEN
( unsigned const iFCT, unsigned const isen, double const& t )
{
  if( !_nf || _vFCT.at( iFCT ).empty() ) return true; // nothing to do if no function

  *_DtS = t; // set current time
  FFVar const* pFCT = _vFCT.at( iFCT ).data();
  //_pIC = _vFCT.at( pos_fct ).data();
  for( unsigned iy=0; iy<_ny; iy++ )   _pSAFCT[iy] = _pY[iy];
  for( unsigned ip=0; ip<_np+1; ip++ ) _pSAFCT[_ny+ip] = (ip==_ndxSEN[isen]? 1.: 0.); // includes time
  for( unsigned iq=0; iq<_nq; iq++ )   _pSAFCT[_nx+_np+1+iq] = _pYQ[iq];
  delete[] _pSACFCT;
  _pSACFCT = _dag->DFAD( _nf, pFCT, _nx+_np+1+_nq, _pVARS+_ny, _pSAFCT );
  static double const one = 1.;
  _dag->eval( _nf, _pSACFCT, _Dfp.data()+isen*_nf, _nVARS, _pVARS, _DVARS, &one );//iFCT? &one: nullptr );

  return true;
}

// 2026-09-28: over-alignment guard (GCC wrong-code with over-aligned virtual bases)
static_assert( alignof( ODESLV_BASE ) <= alignof( void* ), "ODESLV_BASE is a VIRTUAL base of the CRONOS solvers and must not be over-aligned (alignof > 8): GCC (6 to at least 16) emits aligned vector stores in base-object constructors assuming the full alignment, but a virtual-base subobject is only placed at its non-virtual alignment -> SIGSEGV at -O2/-O3 (see gccvb/pr_vbase_align.cpp). Keep over-aligned members (Armadillo/Eigen fixed-size types, alignas) behind a pointer, as FFModel::_pClassification does." );

//! @brief ODESLV_BASE was merged into ODESLV_BASE on 2026-10-21; the alias keeps consumers compiling.
using ODESLVS_BASE = ODESLV_BASE;

} // end namescape mc

#endif

