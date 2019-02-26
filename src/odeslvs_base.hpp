// Copyright (C) 2019 Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

#ifndef MC__ODESLVS_BASE_HPP
#define MC__ODESLVS_BASE_HPP

#undef  MC__ODESLVS_BASE_DEBUG

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
class ODESLVS_BASE:
  public virtual BASE_DE,
  public virtual ODESLV_BASE
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
  //! @brief array of subgraphs of sensitivity/adjoint RHS functions
  std::vector<FFSubgraph> _opSARHS;

  //! @brief subgraph of sensitivity/adjoint RHS Jacobian 
  FFSubgraph _opSAJAC;

  //! @brief array of subgraphs of sensitivity/adjoint quadrature functions
  std::vector<FFSubgraph> _opSAQUAD;

  //! @brief const pointer to RHS function in current stage of ODE system
  const FFVar* _pRHS;

  //! @brief const pointer to quadrature integrand in current stage of ODE system
  const FFVar* _pQUAD;

  //! @brief vector of const pointers to sensitivity/adjoint RHS function in current stage of ODE system
  std::vector<const FFVar*> _vSARHS;

  //! @brief vector of const pointers to sensitivity/adjoint quadrature integrand in current stage of ODE system
  std::vector<const FFVar*> _vSAQUAD;

  //! @brief const pointer to adjoint TC function in current stage of ODE system
  const FFVar* _pSACFCT;

  //! @brief pointer to adjoint discontinuity function in current stage of ODE system
  FFVar* _pSAFCT;

  //! @brief number of variables for DAG evaluation
  unsigned _nVAR;

  //! @brief number of variables for DAG evaluation (without quadratures)
  unsigned _nVAR0;

  //! @brief pointer to variables for DAG evaluation
  FFVar* _pVAR;

  //! @brief pointer of variable values for DAG evaluation
  double *_DVAR;

  //! @brief pointer to time **DO NOT FREE**
  double *_Dt;

  //! @brief pointer to sensitivity/adjoint values **DO NOT FREE**
  double *_Dy;

  //! @brief pointer to state values **DO NOT FREE**
  double *_Dx;

  //! @brief pointer to parameter values **DO NOT FREE**
  double *_Dp;

  //! @brief pointer to quadrature values **DO NOT FREE**
  double *_Dq;

  //! @brief pointer to sensitivity/adjoint quadrature values **DO NOT FREE**
  double *_Dyq;

  //! @brief vector to hold function derivatives
  std::vector<double> _Dfp;

  //! @brief storage vector for sensitivity/adjoint DAG evaluation
  std::vector<double> _DWRK;

  //! @brief Function to initialize sensitivity for parameter <a>isen</a>
  bool _IC_SET_FSA
    ( const unsigned isen );

  //! @brief Function to add initial state contribution to adjoint quadrature values
  bool _IC_SET_ASA
    ();

  //! @brief Function to reinitialize sensitivity at stage times for parameter <a>isen</a>
  bool _CC_SET_FSA
    ( const unsigned pos_ic, unsigned isen );

  //! @brief Function to reinitialize adjoint at stage times for function <a>ifct</a>
  bool _CC_SET_ASA
    ( const unsigned pos_fct, unsigned ifct );

  //! @brief Function to initial adjoint at terminal time for function <a>ifct</a>
  bool _TC_SET_ASA
    ( const unsigned pos_fct, const unsigned ifct );

  //! @brief Function to set sensitivity RHS pointer
  bool _RHS_SET_FSA
    ( const unsigned iRHS, const unsigned iQUAD );

  //! @brief Function to set adjoint RHS pointer
  bool _RHS_SET_ASA
    ( const unsigned iRHS, const unsigned iQUAD,
      const unsigned pos_fct, const bool neg=true );

  //! @brief Function to initialize sensitivity/adjoint values
  bool _INI_D_SEN
    ( const double *p, const unsigned nf, const unsigned nyq );

  //! @brief Function to retreive sensitivity/adjoint values
  template <typename REALTYPE> void _GET_D_SEN
    ( const REALTYPE*y, const unsigned nyq, const REALTYPE*yq );

  //! @brief Function to retreive state and sensitivity/adjoint values
  template <typename REALTYPE> void _GET_D_SEN
    ( const REALTYPE*x, const REALTYPE*y, const REALTYPE*q,
      const unsigned nyq, const REALTYPE*yq );

  //! @brief Function to initialize sensitivity/adjoint values at terminal time
  template <typename REALTYPE> bool _TC_D_SEN
    ( const double t, const REALTYPE *x, REALTYPE *y );

  //! @brief Function to initialize adjoint quadrature values at terminal time
  template <typename REALTYPE> bool _TC_D_QUAD_ASA
   ( REALTYPE*yq ); 

  //! @brief Function to add initial state contribution to sensitivity/adjoint quadrature values
  template <typename REALTYPE> bool _IC_D_SEN
    ( const double t, const REALTYPE*x, const REALTYPE*y );

  //! @brief Function to add initial state contribution to adjoint quadrature values
  template <typename REALTYPE> bool _IC_D_QUAD_ASA
    ( REALTYPE*yq );

  //! @brief Function to set sensitivity/adjoint transitions at stage times
  bool _CC_D_SET
    ();

  //! @brief Function to transition sensitivity/adjoint values at stage times
  template <typename REALTYPE> bool _CC_D_SEN
    ( const double t, const REALTYPE *x, REALTYPE *y );

  //! @brief Function to transition adjoint quadrature values at stage times
  template <typename REALTYPE> bool _CC_D_QUAD_ASA
    ( REALTYPE*yq );

  //! @brief Function to set sensitivity/adjoint RHS pointer
  bool _RHS_D_SET
    ( const unsigned nf, const unsigned nyq );

  //! @brief Function to calculate the RHS of sensitivity/adjoint ODEs
  template <typename REALTYPE> bool _RHS_D_SEN
  (  double t, const REALTYPE*x, const REALTYPE*y, REALTYPE*ydot,
     const unsigned ifct );

  //! @brief Function to calculate the Jacobian RHS of sensitivity/adjoint ODEs
  template <typename REALTYPE> bool _JAC_D_SEN
    ( double t, const REALTYPE* x, const REALTYPE* y, REALTYPE**jac );

  //! @brief Function to calculate the RHS of sensitivity/adjoint quadrature ODEs
  template <typename REALTYPE> bool _RHS_D_QUAD
    ( const unsigned nyq, REALTYPE*qdot, const unsigned ifct );

  //! @brief Function to calculate the function sensitivities at intermediate/end point
  bool _FCT_D_SEN
    ( const unsigned pos_fct, const unsigned isen, const double t );

  //! @brief Private methods to block default compiler methods
  ODESLVS_BASE(const ODESLVS_BASE&);
  ODESLVS_BASE& operator=(const ODESLVS_BASE&);
};

ODESLVS_BASE::ODESLVS_BASE
()
: BASE_DE(), ODESLV_BASE(), _pRHS(0),  _pQUAD(0), _pSACFCT(0),
  _pSAFCT(0), _nVAR(0), _nVAR0(0), _pVAR(0)
{
  _DVAR = _Dt = _Dp = _Dy = _Dx = _Dq = _Dyq = 0;
}

ODESLVS_BASE::~ODESLVS_BASE
()
{
  delete[] _pVAR;
  delete[] _DVAR;
  /* DO NOT FREE _pRHS, _pQUAD */
  for( auto it=_vSARHS.begin(); it!=_vSARHS.end(); ++it ) delete[] *it;
  for( auto it=_vSAQUAD.begin(); it!=_vSAQUAD.end(); ++it ) delete[] *it;
  delete[] _pSACFCT;
  delete[] _pSAFCT;
}

bool
ODESLVS_BASE::_INI_D_SEN
( const double*p, const unsigned nf, const unsigned nyq )
{
  // Size and set DAG evaluation arrays
  BASE_DE::set_sensitivity( _nx, nyq );
  _nVAR0 = _ny + _nx + _np + 1;
  _nVAR  = _nVAR0 + _nq + nyq;
  delete[] _pVAR; _pVAR = new FFVar[_nVAR];
  delete[] _DVAR; _DVAR = new double[_nVAR];
  delete[] _pSAFCT;  _pSAFCT  = new FFVar[_ny+_np+1+_nq];
  for( auto it=_vSARHS.begin(); it!=_vSARHS.end(); ++it ){ delete[] *it; *it=0; }
  for( auto it=_vSAQUAD.begin(); it!=_vSAQUAD.end(); ++it ){ delete[] *it; *it=0; }
  _vSARHS.resize(nf); _vSAQUAD.resize(nf);

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

  _Dfp.resize( _nf*_np );

  return true;
}

template <typename REALTYPE>
inline
void
ODESLVS_BASE::_GET_D_SEN
( const REALTYPE*y, const unsigned nyq, const REALTYPE*yq )
{
  _vec2D( y, _ny, _Dy );
  if( yq ) _vec2D( yq, nyq, _Dyq );
}

template <typename REALTYPE>
inline
void
ODESLVS_BASE::_GET_D_SEN
( const REALTYPE*x, const REALTYPE*y, const REALTYPE*q,
  const unsigned nyq, const REALTYPE*yq )
{
  _vec2D( x, _nx, _Dx );
  _vec2D( y, _ny, _Dy );
  if( q )  _vec2D( q, _nq, _Dq );
  if( yq ) _vec2D( yq, nyq, _Dyq );
}

bool
ODESLVS_BASE::_IC_SET_FSA
( const unsigned isen )
{
  _pIC = _vIC.at(0);
  delete[] _pSACFCT; _pSACFCT = _pDAG->FAD( _nx, _pIC, 1, _pVAR+_ny+_nx+isen);
  _pIC = _pSACFCT;
  return true;
}

bool
ODESLVS_BASE::_IC_SET_ASA
()
{
  const FFVar* pIC = _vIC.at(0);
  FFVar pHAM( 0. );
  for( unsigned ix=0; ix<_nx; ix++ ) pHAM += _pY[ix] * pIC[ix];
#ifndef MC__ODEBNDS_GSL_USE_BAD
  delete[] _pSACFCT; _pSACFCT = _pDAG->FAD( 1, &pHAM, _np, _pP );
#else
  delete[] _pSACFCT; _pSACFCT = _pDAG->BAD( 1, &pHAM, _np, _pP );
#endif
  return true;
}

bool
ODESLVS_BASE::_CC_SET_FSA
( const unsigned pos_ic, const unsigned isen )
{
  _pIC = _vIC.at( pos_ic );
  for( unsigned iy=0; iy<_ny; iy++ )   _pSAFCT[iy] = _pVAR[iy];
  for( unsigned ip=0; ip<_np; ip++ ) _pSAFCT[_ny+ip] = (ip==isen? 1.: 0.);
  delete[] _pSACFCT; _pSACFCT = _pDAG->DFAD( _nx, _pIC, _nx+_np, _pVAR+_ny, _pSAFCT );
  for( unsigned iy=0; iy<_ny; iy++ )   _pSAFCT[iy] = _pSACFCT[iy];

  return true;
}

bool
ODESLVS_BASE::_CC_SET_ASA
( const unsigned pos_fct, const unsigned ifct )
{
  _pIC = _vFCT.at(pos_fct-1)+ifct;
#ifndef MC__ODEBNDS_GSL_USE_BAD
  delete[] _pSACFCT; _pSACFCT = _pDAG->FAD( 1, _pIC, _nx+_np, _pVAR+_ny );
#else
  delete[] _pSACFCT; _pSACFCT = _pDAG->BAD( 1, _pIC, _nx+_np, _pVAR+_ny );
#endif
  for( unsigned iy=0; iy<_nx; iy++ )
    _pSAFCT[iy] = _pY[iy] + _pSACFCT[iy];

  return true;
}

bool
ODESLVS_BASE::_TC_SET_ASA
( const unsigned pos_fct, const unsigned ifct )
{
  const FFVar* _pIC = _vFCT.at(pos_fct)+ifct;

#ifndef MC__ODEBNDS_GSL_USE_BAD
  delete[] _pSACFCT; _pSACFCT = _pDAG->FAD( 1, _pIC, _nx+_np, _pVAR+_nx );
#else
  delete[] _pSACFCT; _pSACFCT = _pDAG->BAD( 1, _pIC, _nx+_np, _pVAR+_nx );
#endif

  return true;
}

bool
ODESLVS_BASE::_RHS_SET_FSA
( const unsigned iRHS, const unsigned iQUAD )
{
  if( _vRHS.size() <= iRHS ) return false; 
  if( _nq && _vQUAD.size() <= iQUAD ) return false;

  _pRHS =  _vRHS.at( iRHS );
#ifdef MC__ODESLVS_BASE_DEBUG
  std::ostringstream ofilename;
  ofilename << "vRHS.dot";
  std::ofstream ofile( ofilename.str(), std::ios_base::out );
  _pDAG->dot_script( _nx, _pRHS, ofile );
  ofile.close();
#endif
  _pQUAD  = _nq? _vQUAD.at( iQUAD ): 0;

  // Set sensitivity ODEs using directional derivatives
  for( unsigned iy=0; iy<_nx; iy++ ) _pSAFCT[iy] = _pVAR[iy];
  for( unsigned ip=0; ip<_np; ip++ ){
    for( unsigned jp=0; jp<_np; jp++ ) _pSAFCT[_nx+jp] = (ip==jp? 1.: 0.);
    delete[] _vSARHS[ip];  _vSARHS[ip]  = _pDAG->DFAD( _nx, _pRHS, _nx+_np, _pVAR+_nx, _pSAFCT );
#ifdef MC__ODESLVS_BASE_DEBUG
    std::ostringstream ofilename;
    ofilename << "vSARHS" << ip << ".dot";
    std::ofstream ofile( ofilename.str(), std::ios_base::out );
    _pDAG->dot_script( _nx, _vSARHS[ip], ofile );
    ofile.close();
#endif
    if( !_nq ) continue;
    delete[] _vSAQUAD[ip]; _vSAQUAD[ip] = _pDAG->DFAD( _nq, _pQUAD, _nx+_np, _pVAR+_nx, _pSAFCT );

  }

  return true;
}

bool
ODESLVS_BASE::_RHS_SET_ASA
( const unsigned iRHS, const unsigned iQUAD,
  const unsigned pos_fct, const bool neg )
{
  if( _vRHS.size() <= iRHS ) return false; 
  if( _nq && _vQUAD.size() <= iQUAD ) return false;

  _pRHS = _vRHS.at( iRHS );
  FFVar pHAM( 0. );
  for( unsigned ix=0; ix<_nx; ix++ ){
    if( neg ) pHAM -= _pVAR[ix] * _pRHS[ix];
    else      pHAM += _pVAR[ix] * _pRHS[ix];
  }
  _pQUAD  = _nq? _vQUAD.at( iQUAD ): 0;
  std::vector<FFVar> vHAM( _nf, pHAM );
  const FFVar* pFCT = _vFCT.at(pos_fct);
  for( unsigned ifct=0; ifct<_nf; ifct++ ){
#ifndef MC__ODEBNDS_GSL_USE_BAD
    delete[] _pSACFCT; _pSACFCT = _nq? _pDAG->FAD( 1, pFCT+ifct, _nq, _pQ ): 0;
#else
    delete[] _pSACFCT; _pSACFCT = _nq? _pDAG->BAD( 1, pFCT+ifct, _nq, _pQ ): 0;
#endif
    for( unsigned iq=0; iq<_nq; iq++ ){
      if( !_pSACFCT[iq].cst() ) return false; // quadrature appears nonlinearly in function
      if( neg ) vHAM[ifct] -= _pQUAD[iq] * _pSACFCT[iq];
      else      vHAM[ifct] += _pQUAD[iq] * _pSACFCT[iq];
    }
  }

  for( unsigned ifct=0; ifct<_nf; ifct++ ){
#ifndef MC__ODEBNDS_GSL_USE_BAD
    delete[] _vSARHS[ifct];  _vSARHS[ifct]  = _pDAG->FAD( 1, vHAM.data()+ifct, _nx, _pX   );
    delete[] _vSAQUAD[ifct]; _vSAQUAD[ifct] = _pDAG->FAD( 1, vHAM.data()+ifct, _np, _pP );
#else
    delete[] _vSARHS[ifct];  _vSARHS[ifct]  = _pDAG->BAD( 1, vHAM.data()+ifct, _nx, _pX   );
    delete[] _vSAQUAD[ifct]; _vSAQUAD[ifct] = _pDAG->BAD( 1, vHAM.data()+ifct, _np, _pP );
#endif
  }

  delete[] std::get<1>(_pJAC); delete[] std::get<2>(_pJAC); delete[] std::get<3>(_pJAC);
  _pJAC = _pDAG->SFAD( _nx, _vSARHS[0], _ny, _pY ); // Jacobian in sparse format

  return true;
}

template <typename REALTYPE>
inline
bool
ODESLVS_BASE::_TC_D_SEN
( const double t, const REALTYPE*x, REALTYPE*y )
{
  *_Dt = t; // current time
  _vec2D( x, _nx, _Dx );
  _pDAG->eval( _ny, _pSACFCT, (double*)y, _nx+_np+1, _pVAR+_ny, _DVAR+_ny );
  return true;
}

template <typename REALTYPE>
inline
bool
ODESLVS_BASE::_TC_D_QUAD_ASA
( REALTYPE*yq )
{
  _pDAG->eval( _np, _pSACFCT+_ny, (double*)yq, _nx+_np+1, _pVAR+_ny, _DVAR+_ny );
  return true;
}

template <typename REALTYPE>
inline
bool
ODESLVS_BASE::_IC_D_SEN
( const double t, const REALTYPE*x, const REALTYPE*y )
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
( REALTYPE*yq )
{
  _vec2D( yq, _np, _Dyq );
  _pDAG->eval( _np, _pSACFCT, (double*)yq, _nVAR0, _pVAR, _DVAR, true );
  _vec2D( yq, _np, _Dyq );
  return true;
}

template <typename REALTYPE>
inline
bool
ODESLVS_BASE::_CC_D_SEN
( const double t, const REALTYPE*x, REALTYPE*y )
{
  *_Dt = t; // current time
  _vec2D( x, _nx, _Dx );
  _vec2D( y, _ny, _Dy );
  _pDAG->eval( _ny, _pSAFCT, (double*)y, _nVAR0, _pVAR, _DVAR );
  return true;
}

template <typename REALTYPE>
inline
bool
ODESLVS_BASE::_CC_D_QUAD_ASA
( REALTYPE*yq )
{
  _vec2D( yq, _np, _Dyq );
  _pDAG->eval( _np, _pSACFCT+_ny, (double*)yq, _nVAR0, _pVAR, _DVAR, true );
  return true;
}

bool
ODESLVS_BASE::_RHS_D_SET
( const unsigned nf, const unsigned nyq )
{
  unsigned nWRK = 0;

  _opSARHS.resize( nf );
  for( unsigned ifct=0; ifct<nf; ifct++ ){
    _opSARHS[ifct]  = _pDAG->subgraph( _ny, _vSARHS[ifct] );
    if( nWRK < _opSARHS[ifct].l_op.size()  ) nWRK = _opSARHS[ifct].l_op.size();
  }
  if( nyq ){
    _opSAQUAD.resize( nf );
    for( unsigned ifct=0; ifct<nf; ifct++ ){
      _opSAQUAD[ifct] = _pDAG->subgraph( nyq, _vSAQUAD[ifct] );
      if( nWRK < _opSAQUAD[ifct].l_op.size() ) nWRK = _opSAQUAD[ifct].l_op.size();
    }
  }

  //_opSAJAC.clear();
  _opSAJAC = _pDAG->subgraph( std::get<0>(_pJAC), std::get<3>(_pJAC) );
  _DJAC.resize( std::get<0>(_pJAC) );
  if( nWRK < _opSAJAC.l_op.size() )  nWRK = _opSAJAC.l_op.size();

  _DWRK.reserve( nWRK );
  return true;
}

template <typename REALTYPE>
inline
bool
ODESLVS_BASE::_RHS_D_SEN
( double t, const REALTYPE* x, const REALTYPE* y, REALTYPE* ydot,
  const unsigned ifct )
{
  if( !_vSARHS.size() || !_vSARHS[ifct] ) return false; // **error** ADJRHS not defined
  *_Dt = t; // set current time
  _vec2D( x, _nx, _Dx ); // set current state bounds
  _vec2D( y, _ny, _Dy ); // set current sensitivity/adjoint bounds
#ifdef MC__ODESLVS_BASE_DEBUG
  _pDAG->output( _opSARHS[ifct] );
#endif
  _pDAG->eval( _opSARHS[ifct], _DWRK, _ny, _vSARHS[ifct], (double*)ydot,
               _nVAR0, _pVAR, _DVAR );
  return true;
}

template <typename REALTYPE>
inline
bool
ODESLVS_BASE::_JAC_D_SEN
( double t, const REALTYPE* x, const REALTYPE* y, REALTYPE** jac )
{
  //*_Dt = t; // current time
  //_vec2D( x, _nx, _Dx ); // current state
  //_vec2D( y, _ny, _Dy ); // set current adjoint bounds
  _pDAG->eval( _opSAJAC, _DWRK, std::get<0>(_pJAC), std::get<3>(_pJAC),
               _DJAC.data(), _nVAR0, _pVAR, _DVAR );
  for( unsigned ie=0; ie<std::get<0>(_pJAC); ++ie ){
    jac[std::get<1>(_pJAC)[ie]][std::get<2>(_pJAC)[ie]] = _DJAC[ie];
#ifdef MC__ODESLVS_BASE_DEBUG
    std::cout << "  jac[" << std::get<1>(_pJAC)[ie] << ", "
              << std::get<2>(_pJAC)[ie] << "] = " << _DJAC[ie] << std::endl;
#endif
  }
  return true;
}

template <typename REALTYPE>
inline
bool
ODESLVS_BASE::_RHS_D_QUAD
( const unsigned nyq, REALTYPE*yqdot, const unsigned ifct )
{
  if( !_vSAQUAD.size() || !_vSAQUAD[ifct] ) return false;
  _pDAG->eval( _opSAQUAD[ifct], _DWRK, nyq, _vSAQUAD[ifct], (double*)yqdot,
               _nVAR0, _pVAR, _DVAR );
  return true;
}

inline
bool
ODESLVS_BASE::_FCT_D_SEN
( const unsigned iFCT, const unsigned isen, const double t )
{
  if( !_nf ) return true;
  *_Dt = t; // set current time
  _pIC = _vFCT.at( iFCT );
  for( unsigned iy=0; iy<_ny; iy++ )   _pSAFCT[iy] = _pY[iy];
  for( unsigned ip=0; ip<_np+1; ip++ ) _pSAFCT[_ny+ip] = (ip==isen? 1.: 0.); // includes time
  for( unsigned iq=0; iq<_nq; iq++ )   _pSAFCT[_nx+_np+1+iq] = _pYQ[iq];
  delete[] _pSACFCT; _pSACFCT = _pDAG->DFAD( _nf, _pIC, _nx+_np+1+_nq, _pVAR+_ny, _pSAFCT );
  _pDAG->eval( _nf, _pSACFCT, _Dfp.data()+isen*_nf, _nVAR, _pVAR, _DVAR, iFCT?true:false );

  return true;
}

} // end namescape mc

#endif

