// Copyright (C) 2016 Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

#ifndef MC__ODESLV_BASE_HPP
#define MC__ODESLV_BASE_HPP

#undef  MC__ODESLV_BASE_DEBUG

#include <stdexcept>
#include <cassert>
#include <cmath>
#include <fstream>
#include <vector>
#include <list>
#include <sys/time.h>

#include "base_de.hpp"

namespace mc
{
//! @brief C++ base class for computing solutions of parametric ODEs using continuous-time real-valued integration.
////////////////////////////////////////////////////////////////////////
//! mc::ODESLV_BASE is a C++ base class for computing solutions of
//! parametric ordinary differential equation (ODEs) using
//! continuous-time real-valued integration.
////////////////////////////////////////////////////////////////////////
class ODESLV_BASE:
  public virtual BASE_DE
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
      ( const double tk, const unsigned nx1, const double*x1,
        const unsigned nx2=0, const double*x2=0 ):
      t( tk ), nx( nx1+nx2 )
      { X = new double[nx];
        unsigned int ix=0;
        for( ; ix<nx1; ix++ ) X[ix] = x1[ix];
        for( ; ix<nx; ix++ )  X[ix] = x2[ix-nx1]; }
    Results
      ( const Results&res ):
      t( res.t ), nx( res.nx )
      { X = new double[nx];
        for( unsigned int ix=0; ix<nx; ix++ ) X[ix] = res.X[ix]; }
    //! @brief Destructor
    ~Results()
      { delete[] X; }
    //! @brief Time point
    double t;
    //! @brief Solution dimension
    unsigned int nx;
    //! @brief Solution point
    double* X;
  };
  /** @} */

 protected:
  //! @brief list of operations in RHS evaluation
  std::list<const FFOp*> _opRHS;

  //! @brief list of operations in RHS Jacobian
  std::list<const FFOp*> _opJAC;

  //! @brief list of operations in quadrature evaluation
  std::list<const FFOp*> _opQUAD;

  //! @brief const pointer to RHS function in current stage of ODE system
  const FFVar* _pRHS;

  //! @brief sparse representation of RHS Jacobian in current stage of ODE system
  std::tuple< unsigned, const unsigned*, const unsigned*, const FFVar* > _pJAC;

  //! @brief const pointer to quadrature integrand in current stage of ODE system
  const FFVar* _pQUAD;

  //! @brief const pointer to IC function in current stage of ODE system
  const FFVar* _pIC;

  //! @brief number of variables for DAG evaluation
  unsigned _nVAR;

  //! @brief number of variables for DAG evaluation (without quadratures)
  unsigned _nVAR0;

  //! @brief pointer to variables for DAG evaluation
  FFVar* _pVAR;

  //! @brief pointer of variable values for DAG evaluation
  double *_DVAR;

  //! @brief pointer to state time **DO NOT FREE**
  double *_Dt;

  //! @brief pointer to state interval bounds **DO NOT FREE**
  double *_Dx;

  //! @brief pointer to parameter interval bounds **DO NOT FREE**
  double *_Dp;

  //! @brief quadrature values **DO NOT FREE**
  double *_Dq;

  //! @brief function values
  double *_Df;

  //! @brief preallocated array for holding Jacobian evaluation results
  double* _DJAC;

  //! @brief preallocated array for DAG evaluation
  double* _DWRK;

  //! @brief Function converting integrator array to internal format
  template <typename REALTYPE> static void _vec2D
    ( const REALTYPE*vec, const unsigned n, double*d );

  //! @brief Function to initialize state integration
  bool _INI_D_STA
    ( const double*p );

  //! @brief Function to retreive state bounds
  template <typename REALTYPE> void _GET_D_STA
    ( const REALTYPE*x, const REALTYPE*q );

  //! @brief Function to set state/quadrature at initial time
  bool _IC_D_SET
    ();

  //! @brief Function to initialize quaratures
  template <typename REALTYPE> bool _IC_D_QUAD
    ( REALTYPE*vec );

  //! @brief Function to initialize states
  template <typename REALTYPE> bool _IC_D_STA
    ( const double t, REALTYPE*vec );

  //! @brief Function to reset state at intermediate time
  bool _CC_D_SET
    ( const unsigned iIC );

  //! @brief Function to reinitialize state at intermediate time
  template <typename REALTYPE> bool _CC_D_STA
    ( const double t, REALTYPE*vec );

  //! @brief Function to set RHS and QUAD pointers
  bool _RHS_D_SET
    ( const unsigned iRHS, const unsigned iQUAD );

  //! @brief Function to set RHS and QUAD pointers
  bool _RHS_D_SET
    ();

  //! @brief Function to calculate the ODE RHS
  template <typename REALTYPE> bool _RHS_D_STA
    ( double t, const REALTYPE*x, REALTYPE*xdot );

  //! @brief Function to calculate the ODE quadrature
  template <typename REALTYPE> bool _RHS_D_QUAD
    ( double t, const REALTYPE*x, REALTYPE*qdot );

  //! @brief Function to calculate the ODE Jacobian
  template <typename REALTYPE> bool _JAC_D_STA
    ( double t, const REALTYPE*x, REALTYPE**jac );

  //! @brief Function to calculate the functions at intermediate/end point
  bool _FCT_D_STA
    ( const unsigned iFCT, const double t );

  //! @brief Record results in file <a>bndrec</a>, with accuracy of <a>iprec</a> digits
  template <typename VRES> static void _record
    ( std::ofstream&ofile, const VRES&bnd, const unsigned iprec=5 );

  //! @brief Private methods to block default compiler methods
  ODESLV_BASE(const ODESLV_BASE&);
  ODESLV_BASE& operator=(const ODESLV_BASE&);
};

inline
ODESLV_BASE::ODESLV_BASE
()
: BASE_DE(), _pRHS(0), _pJAC(0,0,0,0), _pQUAD(0), _pIC(0), _nVAR(0), _nVAR0(0),
  _pVAR(0), _DVAR(0), _Dt(0), _Dx(0), _Dp(0), _Dq(0), _Df(0), _DJAC(0), _DWRK(0)
{}

inline
ODESLV_BASE::~ODESLV_BASE
()
{
  /* DO NOT FREE _pRHS, _pQUAD, _pIC */
  delete[] _DWRK;
  delete[] std::get<1>(_pJAC);  std::get<1>(_pJAC) = 0;
  delete[] std::get<2>(_pJAC);  std::get<2>(_pJAC) = 0;
  delete[] std::get<3>(_pJAC);  std::get<3>(_pJAC) = 0;
  delete[] _Df;
  delete[] _DJAC;
  delete[] _pVAR;
  delete[] _DVAR;   // **DO NOT FREE _Dt, _Dx, _Dp, _Dq**
}

template <typename REALTYPE>
inline void
ODESLV_BASE::_vec2D
( const REALTYPE*vec, const unsigned n, double*d )
{
  for( unsigned i=0; i<n; i++  ) d[i] = vec[i];
  return;
}

inline bool
ODESLV_BASE::_INI_D_STA
( const double*p )
{
  // Size and set DAG evaluation arrays
  _nVAR0 = _nx+_np+1;
  _nVAR  = _nVAR0+_nq;
  delete[] _pVAR; _pVAR = new FFVar[_nVAR];
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

  delete[] _Df; _Df = _nf? new double[_nf]: 0;

  return true;
}

template <typename REALTYPE> inline void
ODESLV_BASE::_GET_D_STA
( const REALTYPE*x, const REALTYPE*q )
{
  _vec2D( x, _nx, _Dx );
  if( q ) _vec2D( q, _nq, _Dq );
}

inline bool
ODESLV_BASE::_IC_D_SET
()
{
  if( !_vIC.size() || _nx0 != _nx ) return false;
  _pIC = _vIC.at(0);
  return true;
}

template <typename REALTYPE> inline bool
ODESLV_BASE::_IC_D_STA
( const double t, REALTYPE*x )
{
  *_Dt = t; // current time
  _pDAG->eval( _nx, _pIC, (double*)x, _np+1, _pVAR+_nx, _DVAR+_nx );
  return true;
}

template <typename REALTYPE> inline bool
ODESLV_BASE::_IC_D_QUAD
( REALTYPE*q )
{
  for( unsigned iq=0; iq<_nq; iq++ ) q[iq] = 0.;
  return true;
}

inline bool
ODESLV_BASE::_CC_D_SET
( const unsigned iIC )
{
  if( _vIC.size() <= iIC || _nx0 != _nx ) return false;
  _pIC = _vIC.at( iIC );
  return true;
}

template <typename REALTYPE> inline bool
ODESLV_BASE::_CC_D_STA
( const double t, REALTYPE*x )
{
  *_Dt = t; // current time
  _vec2D( x, _nx, _Dx ); // current state
  _pDAG->eval( _nx, _pIC, (double*)x, _nVAR0, _pVAR, _DVAR );
  return true;
}

inline bool
ODESLV_BASE::_RHS_D_SET
( const unsigned iRHS, const unsigned iQUAD )
{
  if( _vRHS.size() <= iRHS ) return false;
  _pRHS = _vRHS.at( iRHS );

  if( _nq && _vQUAD.size() <= iQUAD ) return false;
  _pQUAD = _nq? _vQUAD.at( iQUAD ): 0;

  delete[] std::get<1>(_pJAC); delete[] std::get<2>(_pJAC); delete[] std::get<3>(_pJAC);
  _pJAC = _pDAG->SFAD( _nx, _pRHS, _nx, _pX ); // Jacobian in sparse format

  return _RHS_D_SET();
}

inline bool
ODESLV_BASE::_RHS_D_SET
()
{
  unsigned opmax = 0;

  _opRHS.clear();
  _opRHS  = _pDAG->subgraph( _nx, _pRHS );
  if( _opRHS.size() > opmax )  opmax = _opRHS.size();

  _opQUAD.clear();
  if( _pQUAD ) _opQUAD = _pDAG->subgraph( _nq, _pQUAD );
  if( _opQUAD.size() > opmax ) opmax = _opQUAD.size();

  _opJAC.clear();
  _opJAC = _pDAG->subgraph( std::get<0>(_pJAC), std::get<3>(_pJAC) );
  delete[] _DJAC; _DJAC = new double[ std::get<0>(_pJAC) ];
  if( _opJAC.size() > opmax )  opmax = _opJAC.size();

  delete[] _DWRK; _DWRK = new double[ opmax ];
  //std::cout << "size for: " << opmax << std::endl;
  return true;
}

template <typename REALTYPE> inline bool
ODESLV_BASE::_RHS_D_STA
( double t, const REALTYPE*x, REALTYPE*xdot )
{
  if( !_pRHS ) return false;
  *_Dt = t; // current time
  _vec2D( x, _nx, _Dx ); // current state
  _pDAG->eval( _opRHS, _DWRK, _nx, _pRHS, (double*)xdot, _nVAR0, _pVAR, _DVAR );
  return true;
}

template <typename REALTYPE> inline bool
ODESLV_BASE::_RHS_D_QUAD
( double t, const REALTYPE*x, REALTYPE*qdot )
{
  if( !_pQUAD ) return false;
  _pDAG->eval( _opQUAD, _DWRK, _nq, _pQUAD, (double*)qdot, _nVAR0, _pVAR, _DVAR );
  return true;
}

template <typename REALTYPE> inline bool
ODESLV_BASE::_JAC_D_STA
( double t, const REALTYPE*x, REALTYPE**jac )
{
  //*_Dt = t; // current time
  //_vec2D( x, _nx, _Dx ); // current state
  //std::cout << "need for: " << _opJAC.size() << std::endl;
  _pDAG->eval( _opJAC, _DWRK, std::get<0>(_pJAC), std::get<3>(_pJAC), _DJAC, _nVAR0, _pVAR, _DVAR );
  for( unsigned ie=0; ie<std::get<0>(_pJAC); ++ie ){
    jac[std::get<2>(_pJAC)[ie]][std::get<1>(_pJAC)[ie]] = _DJAC[ie];
#ifdef MC__ODESLV_BASE_DEBUG
    std::cout << "  jac[" << std::get<1>(_pJAC)[ie] << ", "
              << std::get<2>(_pJAC)[ie] << "] = " << _DJAC[ie] << std::endl;
#endif
  }
  return true;
}

inline bool
ODESLV_BASE::_FCT_D_STA
( const unsigned iFCT, const double t )
{
  if( !_nf ) return true;
  *_Dt = t; // current time
  const FFVar* pFCT = _vFCT.at( iFCT );
  _pDAG->eval( _nf, pFCT, _Df, _nVAR, _pVAR, _DVAR, iFCT?true:false );
  return true;
}

template <typename VRES> inline void
ODESLV_BASE::_record
( std::ofstream&ofile, const VRES&res, const unsigned iprec )
{
  if( !ofile ) return;

  // Specify format
  ofile << std::right << std::scientific << std::setprecision(iprec);

  // Record computed states at stage times
  auto it = res.begin();
  for( ; it != res.end(); ++it ){
    ofile << std::setw(iprec+9) << it->t;
    for( unsigned ix=0; ix<it->nx; ix++ )
      ofile << std::setw(iprec+9) << it->X[ix];
    ofile << std::endl;
  }
}

} // end namescape mc

#endif

