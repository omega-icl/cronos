// Copyright (C) Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the EPL-2.0 with GPL-2.0-or-later as a Secondary License; see the LICENSE file.
// SPDX-License-Identifier: EPL-2.0 OR GPL-2.0-or-later

#ifndef CRONOS__ODESLVS_BASE_HPP
#define CRONOS__ODESLVS_BASE_HPP

#undef  CRONOS__ODESLVS_BASE_DEBUG

#include "ffexpr.hpp"
#include "odeslv_base.hpp"

namespace mc
{
//! @brief C++ class computing solutions of parametric ODEs with forward/adjoint sensitivity analysis using continuous-time real-valued integration.
////////////////////////////////////////////////////////////////////////
//! mc::ODESLVS_BASE is a C++ base class for computing solutions of
//! parametric ordinary differential equation (ODEs) with forward/
//! adjoint sensitivity analysis using continuous-time real-valued
//! integration.
////////////////////////////////////////////////////////////////////////
class ODESLVS_BASE
: public virtual ODESLV_BASE
, public virtual FFModel
{
public:
  /** @ingroup ODESLV
   *  @{
   */
  //! @brief Default constructor
  ODESLVS_BASE();

  //! @brief Virtual destructor
  virtual ~ODESLVS_BASE();
  /** @} */

protected:
  using ODESLV_BASE::_nx;
  using ODESLV_BASE::_nx0;
  using ODESLV_BASE::_np;
  using ODESLV_BASE::_ndxSEN;
  using ODESLV_BASE::_nsen;
  using ODESLV_BASE::_nq;
  using ODESLV_BASE::_nf;
  using ODESLV_BASE::_t;
  using ODESLV_BASE::_istg;

  using ODESLV_BASE::_dag;
  using ODESLV_BASE::_vIC;
  using ODESLV_BASE::_vRHS;
  using ODESLV_BASE::_vQUAD;
  using ODESLV_BASE::_vFCT;
  using ODESLV_BASE::_pT;
  using ODESLV_BASE::_pX;
  using ODESLV_BASE::_pQ;
  using ODESLV_BASE::_pP;

  using ODESLV_BASE::_vec2D;
  //using ODESLV_BASE::_pIC;
  using ODESLV_BASE::_pJAC;
  using ODESLV_BASE::_pJACCOLNDX;

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
  FFVar const* _pRHS;

  //! @brief const pointer to quadrature integrand in current stage of ODE system
  FFVar const* _pQUAD;

  //! @brief vector of const pointers to sensitivity/adjoint RHS function in current stage of ODE system
  std::vector<FFVar const*> _vSARHS;

  //! @brief vector of const pointers to sensitivity/adjoint quadrature integrand in current stage of ODE system
  std::vector<FFVar const*> _vSAQUAD;

  //! @brief const pointer to adjoint TC function in current stage of ODE system
  FFVar const* _pSACFCT;

  //! @brief pointer to adjoint discontinuity function in current stage of ODE system
  FFVar* _pSAFCT;

  //! @brief number of variables for DAG evaluation
  unsigned _nVAR;

  //! @brief number of variables for DAG evaluation (without quadratures)
  unsigned _nVAR0;

  //! @brief pointer to variables for DAG evaluation
  FFVar* _pVAR;

  //! @brief pointer of variable values for DAG evaluation
  double* _DVAR;

  //! @brief pointer to time **DO NOT FREE**
  double* _Dt;

  //! @brief pointer to sensitivity/adjoint values **DO NOT FREE**
  double* _Dy;

  //! @brief pointer to state values **DO NOT FREE**
  double* _Dx;

  //! @brief pointer to parameter values **DO NOT FREE**
  double* _Dp;

  //! @brief pointer to quadrature values **DO NOT FREE**
  double* _Dq;

  //! @brief pointer to sensitivity/adjoint quadrature values **DO NOT FREE**
  double* _Dyq;

  //! @brief vector to hold function derivatives
  std::vector<double> _Dfp;

  //! @brief storage vector for sensitivity/adjoint DAG evaluation
  std::vector<double> _DWRK;

  //! @brief Function setting up local DAG
  bool _SETUP
    ();

  //! @brief Function setting up DAG of IVP
  bool _SETUP
    ( ODESLVS_BASE const& IVP );

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
  bool _RHS_D_SET
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
  ODESLVS_BASE( ODESLVS_BASE const& ) = delete;
  ODESLVS_BASE& operator=( ODESLVS_BASE const& ) = delete;
};

inline
ODESLVS_BASE::ODESLVS_BASE
()
: ODESLV_BASE(),
  _ny(0), _pY(nullptr), _nyq(0), _pYQ(nullptr),
  _pRHS(nullptr),  _pQUAD(nullptr), _pSACFCT(nullptr), _pSAFCT(nullptr),
  _nVAR(0), _nVAR0(0), _pVAR(nullptr)
{
  _DVAR = _Dt = _Dp = _Dy = _Dx = _Dq = _Dyq = nullptr;
}

inline
ODESLVS_BASE::~ODESLVS_BASE
()
{
  delete[] _pVAR;
  delete[] _DVAR;
  /* DO NOT FREE _pRHS, _pQUAD */
  for( auto& rhs  : _vSARHS  ) delete[] rhs;
  for( auto& quad : _vSAQUAD ) delete[] quad;
  delete[] _pSACFCT;
  delete[] _pSAFCT;
  delete[] _pY;
  delete[] _pYQ;
}

inline
bool
ODESLVS_BASE::_SETUP
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
ODESLVS_BASE::_SETUP
( ODESLVS_BASE const& IVP )
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
ODESLVS_BASE::_INI_D_SEN
( double const* p, unsigned const nvec, unsigned const nyq )
{
  // Size and set DAG evaluation arrays
  if( nyq > _nyq ) return false;
  _nVAR0 = _ny + _nx + _np + 1;
  _nVAR  = _nVAR0 + _nq + nyq;
  delete[] _pVAR; _pVAR = odeslv_new_vars( _nVAR );
  delete[] _DVAR; _DVAR = new double[_nVAR];
  delete[] _pSAFCT;  _pSAFCT  = odeslv_new_vars( _ny+_np+1+_nq );
  for( auto it=_vSARHS.begin(); it!=_vSARHS.end(); ++it ){ delete[] *it; *it=0; }
  for( auto it=_vSAQUAD.begin(); it!=_vSAQUAD.end(); ++it ){ delete[] *it; *it=0; }
  _vSARHS.resize(nvec); _vSAQUAD.resize(nvec);

  for( unsigned ix=0; ix<_nx; ix++ ) _pVAR[ix] = _pY[ix];
  for( unsigned ix=0; ix<_nx; ix++ ) _pVAR[_ny+ix] = _pX[ix];
  for( unsigned ip=0; ip<_np; ip++ ) _pVAR[_ny+_nx+ip] = _pP[ip];
  _pVAR[2*_nx+_np] = (_pT? *_pT: 0. );
  for( unsigned iq=0; iq<_nq; iq++ ) _pVAR[_ny+_nx+_np+1+iq] = _pQ?_pQ[iq]:0.;
  for( unsigned iyq=0; iyq<nyq; iyq++ ) _pVAR[_ny+_nx+_np+1+_nq+iyq] = _pYQ?_pYQ[iyq]:0.;
  _Dy = _DVAR;
  _Dx = _Dy + _ny;
  _Dp = _Dx + _nx;
  _Dt = _Dp + _np;
  _Dq = _Dt + 1;
  _Dyq = _Dq + _nq;
  for( unsigned ip=0; ip<_np; ip++ ) _Dp[ip] = p[ip];

  _Dfp.assign( _nf*_nsen, 0. );

  return true;
}

template <typename REALTYPE>
inline
void
ODESLVS_BASE::_GET_D_SEN
( REALTYPE const* y, unsigned const nyq, REALTYPE const* yq )
{
  _vec2D( y, _ny, _Dy );
  if( yq ) _vec2D( yq, nyq, _Dyq );
}

template <typename REALTYPE>
inline
void
ODESLVS_BASE::_GET_D_SEN
( REALTYPE const* x, REALTYPE const* y, REALTYPE const* q,
  unsigned const nyq, REALTYPE const* yq )
{
  _vec2D( x, _nx, _Dx );
  _vec2D( y, _ny, _Dy );
  if( q )  _vec2D( q, _nq, _Dq );
  if( yq ) _vec2D( yq, nyq, _Dyq );
}

inline
bool
ODESLVS_BASE::_IC_SET_FSA
( unsigned const isen )
{
  FFVar const* pIC = _vIC.at( 0 ).data();
#ifdef CRONOS__ODESLVS_BASE_DEBUG
  auto sgIC = _dag->subgraph( _nx, pIC );
  //_dag->output( sgIC, " pIC" );
  std::vector<FFExpr> exIC = FFExpr::subgraph( _dag, sgIC ); 
  for( size_t i=0; i<_nx; ++i )
    std::cout << "IC[" << i << "] = " << exIC[i] << std::endl;
  std::cout << "P[" << _ndxSEN[isen] << "] = " << _pVAR[_ny+_nx+_ndxSEN[isen]] << std::endl;
  {std::cout << "Enter <1> to continue"; int dum; std::cin >> dum;}
#endif

  delete[] _pSACFCT; _pSACFCT = _dag->FAD( _nx, pIC, 1, _pVAR+_ny+_nx+_ndxSEN[isen] );
  pIC = _pSACFCT;
  return true;
}
/*
inline
bool
ODESLVS_BASE::_IC_SET_ASA
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
ODESLVS_BASE::_IC_SET_ASA
( unsigned const ifct )
{
  FFVar const* pIC = _vIC.at( 0 ).data();
  FFVar pHAM = (!_vFCT.at( 0 ).empty()? _vFCT.at( 0 )[ifct]: 0.);
  for( unsigned ix=0; ix<_nx; ix++ ) pHAM += _pY[ix] * pIC[ix];
  //_dag->output( _dag->subgraph( 1, &pHAM ) );
#ifndef CRONOS__ODESLVS_USE_BAD
  delete[] _pSACFCT; _pSACFCT = _dag->FAD( 1, &pHAM, _nx+_np, _pVAR+_ny );
#else
  delete[] _pSACFCT; _pSACFCT = _dag->BAD( 1, &pHAM, _nx+_np, _pVAR+_ny );
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
ODESLVS_BASE::_CC_SET_FSA
( unsigned const pos_ic, unsigned const isen )
{
  FFVar const* pIC = _vIC.at( pos_ic ).data();
  for( unsigned iy=0; iy<_ny; iy++ ) _pSAFCT[iy] = _pVAR[iy];
  for( unsigned ip=0; ip<_np; ip++ ) _pSAFCT[_ny+ip] = (ip==_ndxSEN[isen]? 1.: 0.);
  delete[] _pSACFCT; _pSACFCT = _dag->DFAD( _nx, pIC, _nx+_np, _pVAR+_ny, _pSAFCT );
  for( unsigned iy=0; iy<_ny; iy++ ) _pSAFCT[iy] = _pSACFCT[iy];

  return true;
}

inline
bool
ODESLVS_BASE::_CC_SET_ASA
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
  std::vector<FFVar> W( _pVAR+_ny, _pVAR+_ny+_nx );
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
ODESLVS_BASE::_TC_SET_ASA
( unsigned const pos_fct, unsigned const ifct )
{
  FFVar const* pFCT = _vFCT.at(pos_fct).data();
  //FFVar const* _pIC = _vFCT.at(pos_fct);//+ifct;

#ifndef CRONOS__ODESLVS_USE_BAD
  delete[] _pSACFCT; _pSACFCT = _dag->FAD( 1, pFCT+ifct, _nx+_np, _pVAR+_nx );
#else
  delete[] _pSACFCT; _pSACFCT = _dag->BAD( 1, pFCT+ifct, _nx+_np, _pVAR+_nx );
#endif

  return true;
}

inline
bool
ODESLVS_BASE::_RHS_SET_FSA
( unsigned const iRHS, unsigned const iQUAD )
{
  if( _vRHS.size() <= iRHS ) return false; 
  if( _nq && _vQUAD.size() <= iQUAD ) return false;

  _pRHS =  _vRHS.at( iRHS ).data();
#ifdef CRONOS__ODESLVS_BASE_DEBUG
  std::ostringstream ofilename;
  ofilename << "vRHS.dot";
  std::ofstream ofile( ofilename.str(), std::ios_base::out );
  _dag->dot_script( _nx, _pRHS, ofile );
  ofile.close();
#endif
  _pQUAD = _nq? _vQUAD.at( iQUAD ).data(): nullptr;

  // Set sensitivity ODEs using directional derivatives
  for( unsigned iy=0; iy<_nx; iy++ ) _pSAFCT[iy] = _pVAR[iy];
  for( unsigned is=0; is<_nsen; is++ ){
    unsigned const ip = _ndxSEN[is];
    for( unsigned jp=0; jp<_np; jp++ ) _pSAFCT[_nx+jp] = (ip==jp? 1.: 0.);
    delete[] _vSARHS[is];  _vSARHS[is]  = _dag->DFAD( _nx, _pRHS, _nx+_np, _pVAR+_nx, _pSAFCT );
#ifdef CRONOS__ODESLVS_BASE_DEBUG
    std::ostringstream ofilename;
    ofilename << "vSARHS" << is << ".dot";
    std::ofstream ofile( ofilename.str(), std::ios_base::out );
    _dag->dot_script( _nx, _vSARHS[ip], ofile );
    ofile.close();
#endif
    if( !_nq ) continue;
    delete[] _vSAQUAD[is]; _vSAQUAD[is] = _dag->DFAD( _nq, _pQUAD, _nx+_np, _pVAR+_nx, _pSAFCT );

  }

  return true;
}

inline
bool
ODESLVS_BASE::_RHS_SET_ASA
( unsigned const iRHS, unsigned const iQUAD,
  unsigned const iFCT, bool const neg )
{
  if( _vRHS.size() <= iRHS ) return false; 
  if( _nq && _vQUAD.size() <= iQUAD ) return false;

  _pRHS = _vRHS.at( iRHS ).data();
  FFVar pHAM( 0. );
  for( unsigned ix=0; ix<_nx; ix++ ){
    if( neg ) pHAM -= _pY[ix] * _pRHS[ix];
    else      pHAM += _pY[ix] * _pRHS[ix];
  }
  _pQUAD  = _nq? _vQUAD.at( iQUAD ).data(): nullptr;
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
        if( neg ) vHAM[ifct] -= _pQUAD[iq] * _pSACFCT[iq];
        else      vHAM[ifct] += _pQUAD[iq] * _pSACFCT[iq];
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
ODESLVS_BASE::_TC_D_SEN
( double const& t, REALTYPE const* x, REALTYPE* y )
{
  if( !_pSACFCT ){
    for( unsigned iy=0; iy<_ny; ++iy )
      y[iy] = 0.;
    return true;
  }
  
  *_Dt = t; // current time
  _vec2D( x, _nx, _Dx );
  _dag->eval( _ny, _pSACFCT, (double*)y, _nx+_np+1, _pVAR+_ny, _DVAR+_ny );
  return true;
}

template <typename REALTYPE>
inline
bool
ODESLVS_BASE::_TC_D_QUAD_ASA
( REALTYPE* yq )
{
  if( !_pSACFCT ){
    for( unsigned iy=0; iy<_nsen; ++iy )
      yq[iy] = 0.;
    return true;
  }

  _dag->eval( _nsen, _pSACFCT+_ny, (double*)yq, _nx+_np+1, _pVAR+_ny, _DVAR+_ny );
  return true;
}

template <typename REALTYPE>
inline
bool
ODESLVS_BASE::_IC_D_SEN
( double const& t, REALTYPE* y )
{
#ifdef CRONOS__ODESLVS_BASE_DEBUG
  auto sgSACFCT = _dag->subgraph( _nx, _pSACFCT );
  //_dag->output( sgSACFCT ), " _pSACFCT" );
  std::vector<FFExpr> exSACFCT = FFExpr::subgraph( _dag, sgSACFCT ); 
  for( size_t i=0; i<_nx; ++i )
    std::cout << "SACFCT[" << i << "] = " << exSACFCT[i] << std::endl;
  for( size_t i=0; i<_np+1; ++i )
    std::cout << "P[" << i << "] = " << _pVAR[_ny+_nx+i] << std::endl;
  {std::cout << "Enter <1> to continue"; int dum; std::cin >> dum;}
#endif

  *_Dt = t; // current time
  _dag->eval( _nx, _pSACFCT, (double*)y, _np+1, _pVAR+_ny+_nx, _DVAR+_ny+_nx );
  return true;
}

template <typename REALTYPE>
inline
bool
ODESLVS_BASE::_IC_D_SEN
( double const& t, REALTYPE const* x, REALTYPE const* y )
{
  *_Dt = t; // current time
  _vec2D( x, _nx, _Dx );
  _vec2D( y, _ny, _Dy );
  return true;
}

template <typename REALTYPE>
inline
bool
ODESLVS_BASE::_IC_D_QUAD_ASA
( REALTYPE* yq )
{
  _vec2D( yq, _nsen, _Dyq );
  static double const one = 1.;
  _dag->eval( _nsen, _pSAFCT, (double*)yq, _nVAR0, _pVAR, _DVAR, &one );
  _vec2D( yq, _nsen, _Dyq );
  return true;
}

template <typename REALTYPE>
inline
bool
ODESLVS_BASE::_CC_D_SEN
( double const& t, REALTYPE const* x, REALTYPE* y )
{
  *_Dt = t; // current time
  _vec2D( x, _nx, _Dx );
  _vec2D( y, _ny, _Dy );
  _dag->eval( _ny, _pSAFCT, (double*)y, _nVAR0, _pVAR, _DVAR );
  return true;
}

template <typename REALTYPE>
inline
bool
ODESLVS_BASE::_CC_D_QUAD_ASA
( REALTYPE* yq )
{
  _vec2D( yq, _nsen, _Dyq );
  static double const one = 1.;
  _dag->eval( _nsen, _pSACFCT+_ny, (double*)yq, _nVAR0, _pVAR, _DVAR, &one );
  return true;
}

inline
bool
ODESLVS_BASE::_RHS_D_SET
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
ODESLVS_BASE::_RHS_D_SEN
( double const& t, REALTYPE const* x, REALTYPE const* y, REALTYPE* ydot,
  unsigned const ifct )
{
  if( !_vSARHS.size() || !_vSARHS[ifct] ) return false; // **error** ADJRHS not defined
  *_Dt = t; // set current time
  _vec2D( x, _nx, _Dx ); // set current state bounds
  _vec2D( y, _ny, _Dy ); // set current sensitivity/adjoint bounds
#ifdef CRONOS__ODESLVS_BASE_DEBUG
  _dag->output( _opSARHS[ifct] );
#endif
  _dag->eval( _opSARHS[ifct], _DWRK, _ny, _vSARHS[ifct], (double*)ydot,
               _nVAR0, _pVAR, _DVAR );
  return true;
}

template <typename REALTYPE>
inline
bool
ODESLVS_BASE::_JAC_D_SEN
( double const& t, REALTYPE const* x, REALTYPE const* y, REALTYPE** jac )
{
  // 2026-09-28: load the callback's OWN (t, x) -- relying on the last RHS call made the Jacobian / quadrature
  // stale whenever another evaluation came in between (the forward replay inside CVodeB: wrong J, replay overrun)
  *_Dt = t;  _vec2D( x, _nx, _Dx );  _vec2D( y, _ny, _Dy );
  _dag->eval( _opSAJAC, _DWRK, _JACB_node.size(), _JACB_node.data(),
               _DJACB.data(), _nVAR0, _pVAR, _DVAR );
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
ODESLVS_BASE::_JAC_D_SEN
( double const& t, REALTYPE const* x, REALTYPE const* y, REALTYPE* jac,
  INDEXTYPE* ptrs, INDEXTYPE* vals )
{
  // 2026-09-28: load the callback's OWN (t, x) -- relying on the last RHS call made the Jacobian / quadrature
  // stale whenever another evaluation came in between (the forward replay inside CVodeB: wrong J, replay overrun)
  *_Dt = t;  _vec2D( x, _nx, _Dx );  _vec2D( y, _ny, _Dy );
  _dag->eval( _opSAJAC, _DWRK, _JACB_node.size(), _JACB_node.data(),
               (double*)jac, _nVAR0, _pVAR, _DVAR );
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
ODESLVS_BASE::_RHS_D_QUAD
( unsigned const nyq, REALTYPE* yqdot, unsigned const ifct )
{
  if( !_vSAQUAD.size() || !_vSAQUAD[ifct] ) return false;
  // _DVAR (t, x, y) is loaded by the CALLER, which has them (CVQUADB__, the sensitivity-quadrature callback)
  _dag->eval( _opSAQUAD[ifct], _DWRK, nyq, _vSAQUAD[ifct], (double*)yqdot,
               _nVAR0, _pVAR, _DVAR );
  return true;
}

inline
bool
ODESLVS_BASE::_FCT_D_SEN
( unsigned const iFCT, unsigned const isen, double const& t )
{
  if( !_nf || _vFCT.at( iFCT ).empty() ) return true; // nothing to do if no function

  *_Dt = t; // set current time
  FFVar const* pFCT = _vFCT.at( iFCT ).data();
  //_pIC = _vFCT.at( pos_fct ).data();
  for( unsigned iy=0; iy<_ny; iy++ )   _pSAFCT[iy] = _pY[iy];
  for( unsigned ip=0; ip<_np+1; ip++ ) _pSAFCT[_ny+ip] = (ip==_ndxSEN[isen]? 1.: 0.); // includes time
  for( unsigned iq=0; iq<_nq; iq++ )   _pSAFCT[_nx+_np+1+iq] = _pYQ[iq];
  delete[] _pSACFCT;
  _pSACFCT = _dag->DFAD( _nf, pFCT, _nx+_np+1+_nq, _pVAR+_ny, _pSAFCT );
  static double const one = 1.;
  _dag->eval( _nf, _pSACFCT, _Dfp.data()+isen*_nf, _nVAR, _pVAR, _DVAR, &one );//iFCT? &one: nullptr );

  return true;
}

// 2026-09-28: over-alignment guard (GCC wrong-code with over-aligned virtual bases)
static_assert( alignof( ODESLVS_BASE ) <= alignof( void* ), "ODESLVS_BASE is a VIRTUAL base of the CRONOS solvers and must not be over-aligned (alignof > 8): GCC (6 to at least 16) emits aligned vector stores in base-object constructors assuming the full alignment, but a virtual-base subobject is only placed at its non-virtual alignment -> SIGSEGV at -O2/-O3 (see gccvb/pr_vbase_align.cpp). Keep over-aligned members (Armadillo/Eigen fixed-size types, alignas) behind a pointer, as FFModel::_pClassification does." );

} // end namescape mc

#endif

