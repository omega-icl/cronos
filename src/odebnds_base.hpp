// Copyright (C) 2015 Benoit Chachuat & Nikola Peric, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

#ifndef MC__ODEBNDS_BASE_HPP
#define MC__ODEBNDS_BASE_HPP

#undef  MC__ODEBNDS_BASE_DINEQI_DEBUG
#undef  MC__ODEBNDS_BASE_DINEQPM_DEBUG
#undef  MC__ODEBNDS_BASE_MVXP_USE

#include <stdexcept>
#include <cassert>
#include <cmath>
#include <fstream>
#include <vector>
#include <sys/time.h>

#include "ellipsoid.hpp"
#include "odebnd_base.hpp"

namespace mc
{
//! @brief C++ base class for computing enclosures of the reachable set of parametric ODEs and their sensitivities using continuous-time set-valued integration.
////////////////////////////////////////////////////////////////////////
//! mc::ODEBNDS_BASE is a C++ base class for computing enclosures of the
//! reachable set of parametric ordinary differential equations (ODEs)
//! and their sensitivities using continuous-time set-valued integration.
////////////////////////////////////////////////////////////////////////
template <typename T, typename PMT, typename PVT>
class ODEBNDS_BASE:
  public virtual BASE_DE,
  public virtual ODEBND_BASE<T,PMT,PVT>
{
  typedef Ellipsoid E;
  typedef BASE_DE::STATUS STATUS;

protected:
  using ODEBND_BASE<T,PMT,PVT>::NORMAL; 
  using ODEBND_BASE<T,PMT,PVT>::FAILURE;
  using ODEBND_BASE<T,PMT,PVT>::FATAL;

  using ODEBND_BASE<T,PMT,PVT>::_npar;
  using ODEBND_BASE<T,PMT,PVT>::_nVAR;
  using ODEBND_BASE<T,PMT,PVT>::_pVAR;

  using ODEBND_BASE<T,PMT,PVT>::_IVAR;
  using ODEBND_BASE<T,PMT,PVT>::_It;
  using ODEBND_BASE<T,PMT,PVT>::_pRHS;
  using ODEBND_BASE<T,PMT,PVT>::_pQUAD;
  using ODEBND_BASE<T,PMT,PVT>::_PMIC;
  using ODEBND_BASE<T,PMT,PVT>::_Q;
  using ODEBND_BASE<T,PMT,PVT>::_Er;
  using ODEBND_BASE<T,PMT,PVT>::_Ir;
  using ODEBND_BASE<T,PMT,PVT>::_pref;
  using ODEBND_BASE<T,PMT,PVT>::_Ip;
  using ODEBND_BASE<T,PMT,PVT>::_B;
  using ODEBND_BASE<T,PMT,PVT>::_xref;
  using ODEBND_BASE<T,PMT,PVT>::_Ixdot;
  using ODEBND_BASE<T,PMT,PVT>::_diam;
  using ODEBND_BASE<T,PMT,PVT>::_vec2I;
  using ODEBND_BASE<T,PMT,PVT>::_vec2E;
  using ODEBND_BASE<T,PMT,PVT>::_ep2x;
  using ODEBND_BASE<T,PMT,PVT>::_E2vec;
  using ODEBND_BASE<T,PMT,PVT>::_I2vec;
  using ODEBND_BASE<T,PMT,PVT>::_IC_I_ELL;
  using ODEBND_BASE<T,PMT,PVT>::_CC_I_ELL;
  using ODEBND_BASE<T,PMT,PVT>::_QUAD_I;
  using ODEBND_BASE<T,PMT,PVT>::_RHS_I_DI;
  using ODEBND_BASE<T,PMT,PVT>::_RHS_I_ELL;

  using ODEBND_BASE<T,PMT,PVT>::_PMenv;
  using ODEBND_BASE<T,PMT,PVT>::_PMVAR;
  using ODEBND_BASE<T,PMT,PVT>::_PMt;
  using ODEBND_BASE<T,PMT,PVT>::_PMx;
  using ODEBND_BASE<T,PMT,PVT>::_PMxdot;
  using ODEBND_BASE<T,PMT,PVT>::_PMp;
  using ODEBND_BASE<T,PMT,PVT>::_vec2PMI;
  using ODEBND_BASE<T,PMT,PVT>::_vec2PME;
  using ODEBND_BASE<T,PMT,PVT>::_PMI2vec;
  using ODEBND_BASE<T,PMT,PVT>::_PME2vec;
  using ODEBND_BASE<T,PMT,PVT>::_IC_PM_ELL;
  using ODEBND_BASE<T,PMT,PVT>::_CC_PM_ELL;
  using ODEBND_BASE<T,PMT,PVT>::_QUAD_PM;
  using ODEBND_BASE<T,PMT,PVT>::_RHS_PM_DI;
  using ODEBND_BASE<T,PMT,PVT>::_RHS_PM_ELL0;
  using ODEBND_BASE<T,PMT,PVT>::_RHS_PM_ELL1;
  using ODEBND_BASE<T,PMT,PVT>::_RHS_PM_ELL2;
  using ODEBND_BASE<T,PMT,PVT>::_RHS_PM_ELL;
  using ODEBND_BASE<T,PMT,PVT>::_e2x;

public:
  //! @brief Default constructor
  ODEBNDS_BASE();

  //! @brief Virtual destructor
  virtual ~ODEBNDS_BASE();

protected:
  //! @brief list of operations in adjoint RHS evaluation
  std::list<const FFOp*> _opADJRHS;

  //! @brief array of list of operations in individual adjoint RHS evaluations
  std::list<const FFOp*> *_opADJRHSi;

  //! @brief list of operations in adjoint RHS Jacobian
  std::list<const FFOp*> _opADJJAC;

  //! @brief list of operations in adjoint quadrature evaluation
  std::list<const FFOp*> _opADJQUAD;

  //! @brief list of operations in adjoint TC evaluation
  std::list<const FFOp*> _opADJTC;

  //! @brief list of operations in adjoint TC Jacobian
  std::list<const FFOp*> _opADJDTC;

  //! @brief list of operations in adjoint quadrature TC evaluation
  std::list<const FFOp*> _opTCQUAD;

  //! @brief preallocated array for evaluation of adjoint RHS function in T arithmetic
  T* _IADJRHS;

  //! @brief preallocated array for evaluation of adjoint RHS Jacobian in T arithmetic
  T* _IADJJAC;

  //! @brief preallocated array for evaluation of adjoint TC Jacobian in T arithmetic
  T* _IADJDTC;

  //! @brief preallocated array for evaluation of adjoint RHS function in PM arithmetic
  PVT* _PMADJRHS;

  //! @brief preallocated array for evaluation of adjoint RHS Jacobian in PM arithmetic
  PVT* _PMADJJAC;

  //! @brief preallocated array for evaluation of adjoint TC function in PM arithmetic
  PVT* _PMADJTC;

  //! @brief preallocated array for evaluation of adjoint TC Jacobian in PM arithmetic
  PVT* _PMADJDTC;

  //! @brief const pointer to adjoint RHS function in current stage of ODE system
  const FFVar* _pADJRHS;

  //! @brief const pointer to adjoint RHS Jacobian in current stage of ODE system
  const FFVar* _pADJJAC;

  //! @brief const pointer to adjoint quadrature integrand in current stage of ODE system
  const FFVar* _pADJQUAD;

  //! @brief const pointer to adjoint TC function in current stage of ODE system
  const FFVar* _pADJTC;

  //! @brief pointer to adjoint discontinuity function in current stage of ODE system
  FFVar* _pADJCC;

  //! @brief const pointer to adjoint TC Jacobian in current stage of ODE system
  const FFVar* _pADJDTC;

  //! @brief pointer to adjoint interval bounds **DO NOT FREE**
  T *_Iy;

  //! @brief pointer to adjoint parameter (states and parameters) interval bounds **DO NOT FREE**
  T *_Iz;

  //! @brief adjoint interval bounds time derivatives
  T *_Iydot;

  //! @brief adjoint lower bound time derivatives
  double *_yLdot;

  //! @brief adjoint upper bound time derivatives
  double *_yUdot;

  //! @brief pointer to adjoint quadrature interval bounds **DO NOT FREE**
  T *_Iyq;

  //! @brief adjoint quadrature derivative interval bounds
  T *_Iyqdot;

  //! @brief adjoint polynomial model **DO NOT FREE**
  PVT *_PMy;

  //! @brief adjoint parameter (states-parameters) polynomial model **DO NOT FREE**
  PVT *_PMz;

  //! @brief adjoint derivative polynomial model
  PVT *_PMydot;

  //! @brief adjoint PM remainder lower bound time derivatives
  double *_RyLdot;

  //! @brief adjoint PM remainder upper bound time derivatives
  double *_RyUdot;

  //! @brief adjoint quadrature polynomial model **DO NOT FREE**
  PVT *_PMyq;

  //! @brief adjoint quadrature derivative polynomial model
  PVT *_PMyqdot;

  //! @brief adjoint quadrature PM remainder radius time derivatives
  double *_Ryqdot;

  //! @brief adjoint parameter reference
  double *_zref;

  //! @brief adjoint reference
  double *_yref;

  //! @brief linear transformation A matrix (adjoint system)
  double *_Ay;

  //! @brief linear transformation B matrix (adjoint system)
  double *_By;

  //! @brief linear transformed adjoint interval bounds
  T *_Idy;

  //! @brief adjoint RHS Jacobian interval bounds
  T *_Idgdy;

  //! @brief linear transformed adjoint ellipsoidal bounds
  E _Edy;

  //! @brief shape matrix (lower triangular) in ellipsoidal bounds (adjoint system)
  double *_Qy;

  //! @brief adjoint reference time derivatives
  double *_yrefdot;

  //! @brief linear transformation B matrix time derivatives (adjoint system)
  double *_Bydot;

  //! @brief rotated adjoint interval bounds time derivatives
  T *_Idydot;

  //! @brief shape matrix time directives in ellipsoidal bounds (adjoint system)
  double *_Qydot;

  //! @brief polynomial model environment for mean-value theorem in (Y,Z)
  PMT *_MVYZenv;

  //! @brief rotated state polynomial model
  PVT *_MVYZr;

  //! @brief rotated adjoint polynomial model
  PVT *_MVYZd;

  //! @brief adjoint RHS polynomial model
  PVT *_MVYZf;

  //! @brief adjoint RHS Jacobian polynomial model
  PVT *_MVYZdgdy;

  //! @brief polynomial models for variables in mean-value theorem (adjoint system)
  PVT *_MVYZVAR;

  //! @brief pointer to time polynomial models (adjoint system) **DO NOT FREE**
  PVT *_MVYZt;

  //! @brief pointer to adjoint polynomial models **DO NOT FREE**
  PVT *_MVYZy;

  //! @brief pointer to parameter polynomial models (adjoint system) **DO NOT FREE**
  PVT *_MVYZz;

  //! @brief pointer to parameter polynomial models (adjoint system) **DO NOT FREE**
  PVT *_MVYZp;

  //! @brief polynomial model environment for mean-value theorem in Z
  PMT *_MVZenv;

  //! @brief pointer to rotated state polynomial models (adjoint initialization)
  PVT *_MVZd;

  //! @brief pointer to adjoint polynomial models (adjoint initialization)
  PVT *_MVZy;

  //! @brief polynomial models for variables in mean-value theorem (adjoint initialization)
  PVT *_MVZVAR;

  //! @brief pointer to parameter polynomial models (adjoint initialization) **DO NOT FREE**
  PVT *_MVZz;

  //! @brief pointer to parameter polynomial models (adjoint initialization) **DO NOT FREE**
  PVT *_MVZp;

  //! @brief pointer to time polynomial models (adjoint initialization) **DO NOT FREE**
  PVT *_MVZt;

  //! @brief number of parameters in adjoint ODE system
  unsigned _nz;

  //! @brief Static function converting rotated states into original coordinates
  template<typename U> static void _ep2x
    ( const unsigned ny, const unsigned nx, const unsigned np, const U*d,
      const double*xref, const U*x, const double*Bx, const double*pref,
      const U*p, const double*Bp, const double*yref, U*y );

  //! @brief Static function converting ellipsoids into interval bounds
  static void _ep2x
    ( const unsigned ny, const unsigned nx, const unsigned np, const double*Q,
      E&Ed, T*Id, const double*xref, const T*Ix, const double*Bx,
      const double*pref, const T*Ip, const double*Bp, const double*yref, T*Iy );

  //! @brief Static function converting GSL array to ellipsoidal bounds
  template <typename REALTYPE> static void _vec2E
    ( const REALTYPE*vec, const unsigned ny, const unsigned nx, const unsigned np,
      double*Q, E&Ed, T*Id, const double*xref, const T*Ix, double*Bx,
      const double*pref, const T*Ip, double*Bp, double*yref, T*Iy );

  //! @brief Function to initialize adjoint interval bounding
  template <typename OPT> bool _INI_I_ADJ
    ( const OPT &options, const unsigned np, const T *Ip );

  //! @brief Function to initialize adjoint interval bounds
  template <typename REALTYPE, typename OPT> bool _TC_I_ADJ
    ( const OPT &options, const double t, REALTYPE *vec,
      unsigned pos_fct, unsigned ifct );

  //! @brief Function to initialize quarature interval bounds
  template <typename REALTYPE, typename OPT> bool _TC_I_QUAD
   ( const OPT&options, const double t, REALTYPE*vec );

  //! @brief Function to add initial state contribution to function derivatives
  bool _IC_I_ADJ
    ();

  //! @brief Function to reinitialize adjoint bounds after discontinuity
  template <typename REALTYPE, typename OPT> bool _CC_I_ADJ
    ( const OPT &options, const double t, REALTYPE *vec,
      unsigned pos_fct, unsigned ifct );

  //! @brief Function to reinitialize adjoint bounds after discontinuity
  template <typename REALTYPE, typename OPT> bool _CC_I_QUAD
    ( const OPT&options, const double t, REALTYPE*vec );

  //! @brief Set adjoint RHS pointer and corresponding Jacobian
  template <typename OPT> bool _SET_I_ADJ
    ( const OPT &options, unsigned iADJRHS, unsigned iQUAD,
      unsigned pos_fct, unsigned ifct );

  //! @brief Function to calculate the adjoint ODEs RHS values in interval arithmetic
  template <typename REALTYPE, typename OPT> bool _RHS_I_ADJ
  (  const OPT &options, double t, const REALTYPE* y, REALTYPE* ydot,
     REALTYPE* vec_sta  );

  //! @brief Function to calculate the RHS of auxiliary ODEs in interval arithmetic
  template <typename REALTYPE, typename OPT> bool _RHS_I_QUAD
    ( const OPT&options, double t, const REALTYPE*y, REALTYPE*qdot,
      REALTYPE* vec_sta, const bool bndinit=true );

  //! @brief Function to initialize adjoint polynomial models
  template <typename OPT> bool _INI_PM_ADJ
    ( const OPT &options, const PVT *PMp );

  //! @brief Function to initialize adjoint polynomial model
  template <typename REALTYPE, typename OPT> bool _TC_PM_ADJ
    ( const OPT &options, const double t, REALTYPE *vec,  unsigned pos_fct, unsigned ifct );

  //! @brief Function to initialize quadrature polynomial model
  template <typename REALTYPE, typename OPT> bool _TC_PM_QUAD
    ( const OPT &options, REALTYPE *vec );

  //! @brief Function to add initial state contribution to function derivatives
  bool _IC_PM_ADJ
    ();

  //! @brief Function to reinitialize state polynomial bounds
  template <typename REALTYPE, typename OPT> bool _CC_PM_ADJ
    ( const OPT &options, const double t, REALTYPE *vec, unsigned pos_fct, unsigned ifct );

  //! @brief Function to set adjoint RHS pointer and corresponding polynomial
  template <typename OPT> bool _SET_PM_ADJ
    ( const OPT &options, unsigned iADJRHS, unsigned iQUAD, unsigned pos_fct, unsigned ifct );

  //! @brief Function to calculate the RHS of adjoint ODEs in polynomial mode arithmetic
  template <typename REALTYPE, typename OPT> bool _RHS_PM_ADJ
  (  const OPT &options, double t, const REALTYPE* y, REALTYPE* ydot, REALTYPE* vec_sta  );

  //! @brief Function to calculate the RHS of auxiliary ODEs in polynomial mode arithmetic
  template <typename REALTYPE, typename OPT> bool _RHS_PM_QUAD
    ( const OPT&options, double t, const REALTYPE*y, REALTYPE*qdot,
      const bool bndinit=true );

  //! @brief Private methods to block default compiler methods
  ODEBNDS_BASE(const ODEBNDS_BASE&);
  ODEBNDS_BASE& operator=(const ODEBNDS_BASE&);
};

template <typename T, typename PMT, typename PVT> inline
ODEBNDS_BASE<T,PMT,PVT>::ODEBNDS_BASE
()
: ODEBND_BASE<T,PMT,PVT>(), _pADJRHS(0), _pADJJAC(0), _pADJQUAD(0), _pADJTC(0),
  _pADJCC(0), _pADJDTC(0)
{
  // Initialize adjoint arrays
  _opADJRHSi = 0;
  _IADJRHS = _IADJJAC = _Iy = _Iz  = _IADJDTC = 0;
  _PMADJRHS = _PMADJJAC = _PMy = _PMz = _PMydot = _PMyq = _PMyqdot =
  _PMADJTC = _PMADJDTC = 0;

  // Initialize parameterization arrays
  _zref = _yref = _yrefdot = 0;
  _Ay = _By = _Qy = _Bydot = _Qydot = 0;
  _Idy = _Idydot = _Iydot = _Iyqdot = _Idgdy = 0;
  _yLdot = _yUdot = _RyLdot = _RyUdot = _Ryqdot = 0;

  // Initialize Taylor model environments
  _MVYZenv = _MVZenv = 0;
  _MVYZd = _MVYZr = _MVYZf = _MVYZdgdy = _MVYZVAR = _MVYZy = _MVYZz = _MVYZp = _MVYZt = 0;
  _MVZd  = _MVZy  = _MVZVAR  = _MVZz  = _MVZp = _MVZt = 0;

  // Initialize ellipsoidal calculus
  E::options.PSDCHK = false;
}

template <typename T, typename PMT, typename PVT> inline
ODEBNDS_BASE<T,PMT,PVT>::~ODEBNDS_BASE
()
{
  // ** DO NOT DELETE _pADJRHS **
  delete[] _opADJRHSi;
  delete[] _IADJRHS;
  delete[] _PMADJRHS;
  delete[] _IADJJAC;
  delete[] _PMADJJAC;
  delete[] _pADJJAC;
  delete[] _pADJRHS;
  delete[] _pADJCC;
  delete[] _pADJTC;
  delete[] _pADJDTC;
  delete[] _IADJDTC;
  delete[] _PMADJTC;
  delete[] _PMADJDTC;

  // Free adjoint arrays
  delete[] _Iydot;
  delete[] _yLdot;
  delete[] _yUdot;
  delete[] _Iyqdot;
  delete[] _PMydot;
  delete[] _RyLdot;
  delete[] _RyUdot;
  delete[] _PMyqdot;
  delete[] _Ryqdot;
      
  // Free linear transformation arrays
  delete[] _zref;
  delete[] _yref;
  delete[] _Ay;
  delete[] _By;
  delete[] _Qy;
  delete[] _yrefdot;
  delete[] _Bydot;
  delete[] _Idy;
  delete[] _Idydot;
  delete[] _Qydot;
  delete[] _Idgdy;
  // ** DO NOT DELETE _MVYZy, _MVYZz, _MVYZp, _MVYZt **
  delete   _MVYZenv;
  delete[] _MVYZd;
  delete[] _MVYZr;
  delete[] _MVYZf;
  delete[] _MVYZdgdy;
  delete[] _MVYZVAR;
  // ** DO NOT DELETE _MVZz, _MVZp, _MVZt **
  delete   _MVZenv;
  delete[] _MVZd;
  delete[] _MVZy;
  delete[] _MVZVAR;
}

template <typename T, typename PMT, typename PVT>
template <typename U> inline void
ODEBNDS_BASE<T,PMT,PVT>::_ep2x
( const unsigned ny, const unsigned nx, const unsigned np, const U*d,
  const double*xref, const U*x, const double*Bx, const double*pref,
  const U*p, const double*Bp, const double*yref, U*y )
{
  for( unsigned iy=0; iy<ny; iy++ ){   
    y[iy] = yref[iy] + d[iy];
    for( unsigned jp=0; jp<np; jp++ ) {y[iy] += ( p[jp] - (pref?pref[jp]:0.) ) * Bp[jp*ny+iy];}
    for( unsigned jx=0; jx<nx; jx++ ) {y[iy] += ( x[jx] - (xref?xref[jx]:0.) ) * Bx[jx*ny+iy];}
  }
  return;
}

template <typename T, typename PMT, typename PVT>
inline void
ODEBNDS_BASE<T,PMT,PVT>::_ep2x
( const unsigned ny, const unsigned nx, const unsigned np, const double*Q,
  E&Ed, T*Id, const double*xref, const T*Ix, const double*Bx,
  const double*pref, const T*Ip, const double*Bp, const double*yref, T*Iy )
{
  Ed.set( ny, Q );
  for( unsigned iy=0; iy<ny; iy++ ) {Id[iy] = T( Ed.l(iy), Ed.u(iy) );}
  _ep2x( ny, nx, np, Id, xref, Ix, Bx, pref, Ip, Bp, yref, Iy );

  return;
}

template <typename T, typename PMT, typename PVT>
template <typename REALTYPE>
inline void
ODEBNDS_BASE<T,PMT,PVT>::_vec2E
( const REALTYPE*vec, const unsigned ny, const unsigned nx, const unsigned np,
  double*Q, E&Ed, T*Id, const double*xref, const T*Ix, double*Bx,
  const double*pref, const T*Ip, double*Bp, double*yref, T*Iy )
{
  unsigned ivec = 0;
  for( unsigned iy=0; iy<ny; iy++ ){yref[iy] = vec[ivec++];}
  for( unsigned iQ=0; iQ<ny*(ny+1)/2; iQ++ ){Q[iQ] = vec[ivec++];}
  for( unsigned iB=0; iB<ny*nx; iB++ ) {Bx[iB] = vec[ivec++];}
  for( unsigned iB=0; iB<ny*np; iB++ ) {Bp[iB] = vec[ivec++];}

  return _ep2x( ny, nx, np, Q, Ed, Id, xref, Ix, Bx, pref, Ip, Bp, yref, Iy );
}


template <typename T, typename PMT, typename PVT> 
template <typename OPT> inline bool
ODEBNDS_BASE<T,PMT,PVT>::_INI_I_ADJ
( const OPT&options, const unsigned np, const T* Ip )
{ 
  // Update effective number of parameters
  // (possibly larger than _np if lifting is used)
  _npar = np;

  // Size and set DAG evaluation arrays
  _nz = _nx + _npar;
  _nVAR = _nx + _nz + 1 + _nq + _npar;
  delete[] _pVAR; _pVAR = new FFVar[_nVAR];
  delete[] _IVAR; _IVAR = new T[_nVAR];
  delete[] _pADJCC;  _pADJCC  = new FFVar[_nx];

  BASE_DE::set_adjoint();
  for( unsigned ix=0; ix<_nx; ix++ ) _pVAR[ix] = _pY[ix];
  for( unsigned ix=0; ix<_nx; ix++ ) _pVAR[_nx + ix] = _pX[ix];
  for( unsigned ip=0; ip<_npar; ip++ ) _pVAR[_nx + _nx + ip] = _pP[ip];
  _pVAR[_nx + _nz] = (_pT? *_pT: 0. );
  _Iy = _IVAR;
  _Iz = _Iy + _nx;
  _Ip = _Iz + _nx;
  _It = _Ip + _npar;
  _Iyq = _It + 1;
  for( unsigned ip=0; ip<_npar; ip++ ) _Ip[ip] = Ip[ip];

  // Reset _MVYZVenv and related variables
  if( _MVYZenv && ( _MVYZenv->nord() != options.ORDMIT 
                 || _MVYZenv->nvar() != _nx+_nz ) ){
    delete[] _MVYZf;   _MVYZf = 0;
    delete[] _MVYZd;   _MVYZd = 0;
    delete[] _MVYZr;   _MVYZr = 0;
    delete[] _MVYZVAR; _MVYZVAR = 0;
    delete   _MVYZenv; _MVYZenv = 0;
  }

  // Reset _MVZVenv and related variables
  if( _MVZenv && ( _MVZenv->nord() != 1 
                || _MVZenv->nvar() != _nz ) ){
    delete[] _MVZy;   _MVZy = 0;
    delete[] _MVZd;   _MVZd = 0;
    delete[] _MVZVAR; _MVZVAR = 0;
    delete   _MVZenv; _MVZenv = 0;
  }

  // Set parameterization variables
  delete[] _Iyqdot; _Iyqdot = new T[_npar];
  switch( options.WRAPMIT){
  case OPT::NONE:
    delete[] _Iydot;   _Iydot = new T[_nx];
    break;
  case OPT::DINEQ:
    delete[] _yLdot;   _yLdot = new double[_nx];
    delete[] _yUdot;   _yUdot = new double[_nx];
    delete[] _Iydot;   _Iydot = new T[_nx];
    break;
  case OPT::ELLIPS:
  default:
    if( !_zref )    _zref     = new double[_nz];
    for( unsigned ip=0; ip<_npar; ip++ )
      _pref[ip] = _zref[_nx+ip] = Op<T>::mid( _Ip[ip] );
    if( !_yref )    _yref     = new double[_nx];
    if( !_Ay )	    _Ay       = new double[_nx*_nx];
    if( !_By )	    _By       = new double[_nx*_nz];
    if( !_Qy )	    _Qy       = new double[_nx*(_nx+1)/2];
    if( !_Idy )      _Idy     = new T[_nx];
    if( !_yrefdot ) _yrefdot  = new double[_nx];
    if( !_Bydot )    _Bydot   = new double[_nx*_nz];
    if( !_Qydot )    _Qydot   = new double[_nx*(_nx+1)/2];
    if( !_Idydot )   _Idydot  = new T[_nx];
    if( !_MVYZenv ) _MVYZenv  = new PMT( _nx+_nz, options.ORDMIT );
    _MVYZenv->options = options.PMOPT;
    if( !_MVYZd )   _MVYZd    = new PVT[_nx];
    if( !_MVYZr )   _MVYZr    = new PVT[_nx];
    if( !_MVYZf )   _MVYZf    = new PVT[_nx];
    if( !_MVYZVAR ) _MVYZVAR  = new PVT[_nVAR];
    _MVYZy = _MVYZVAR;
    _MVYZz = _MVYZy + _nx;
    _MVYZp = _MVYZz + _nx;
    _MVYZt = _MVYZp + _npar;
    if( !_MVZenv )  _MVZenv   = new PMT( _nz, 1 );
    _MVZenv->options = options.PMOPT;
    if( !_MVZd )    _MVZd     = new PVT[_nx];
    if( !_MVZy )    _MVZy     = new PVT[_nx];
    if( !_MVZVAR )  _MVZVAR   = new PVT[_nz+1];
    _MVZz = _MVZVAR;
    _MVZp = _MVZz + _nx;
    _MVZt = _MVZp + _npar;
    break;
  }

  return true;
}

template <typename T, typename PMT, typename PVT>
template <typename REALTYPE, typename OPT> inline bool
ODEBNDS_BASE<T,PMT,PVT>::_TC_I_ADJ
( const OPT&options, const double t, REALTYPE*vec,
  unsigned pos_fct, unsigned ifct )
{
  const FFVar* pFCT = _vFCT.at(pos_fct)+ifct;
#ifndef MC__ODEBNDS_BASE_USE_BAD
  delete[] _pADJTC; _pADJTC = _pDAG->FAD( 1, pFCT, _nz, _pVAR+_nx );
#else
  delete[] _pADJTC; _pADJTC = _pDAG->BAD( 1, pFCT, _nz, _pVAR+_nx );
#endif
  _opADJTC  = _pDAG->subgraph( _nx, _pADJTC );
  *_It = -t; // current time

  switch( options.WRAPMIT){

  case OPT::NONE:
  case OPT::DINEQ:
    _pDAG->eval( _opADJTC, _nx, _pADJTC, _Iy, _nz+1, _pVAR+_nx, _IVAR+_nx );
    _I2vec( _nx, _Iy, vec );
    break;

  case OPT::ELLIPS:
  default:
    *_MVZt = -t;
    for( unsigned ip=0; ip<_npar; ip++ )
      _MVZp[ip].set( _MVZenv, _nx+ip, _Ip[ip] );
    for( unsigned jx=0; jx<_nx; jx++ )
      _MVZd[jx].set( _MVZenv, jx, _Ir[jx] );
    _ep2x( _nx, _npar, _MVZd, _pref, _MVZp, _B, _zref, _MVZz );
    _IC_I_ELL( _pDAG, _opADJTC, _nx, _pADJTC, _nz+1, _pVAR+_nx, _MVZVAR, _MVZy,
               _yref, _Qy, _nz, _By, _Iy );
    _E2vec( _nx, _nz, _yref, _Qy, _By, vec );
    break;
  }

  return true;
}

template <typename T, typename PMT, typename PVT>
template <typename REALTYPE, typename OPT> inline bool
ODEBNDS_BASE<T,PMT,PVT>::_TC_I_QUAD
( const OPT&options, const double t, REALTYPE*vec )
{
  *_It = -t; // current time
  _opTCQUAD = _pDAG->subgraph( _npar, _pADJTC+_nx );
  _pDAG->eval( _opTCQUAD, _npar, _pADJTC+_nx, _Iyq, _nz+1, _pVAR+_nx, _IVAR+_nx );
  _I2vec( _npar, _Iyq, vec );
  return true;
}

template <typename T, typename PMT, typename PVT> inline bool
ODEBNDS_BASE<T,PMT,PVT>::_IC_I_ADJ
()
{
  const FFVar* pIC = _vIC.at(0);
  FFVar pHAM( 0. );
  for( unsigned ix=0; ix<_nx; ix++ ) pHAM += _pVAR[ix] * pIC[ix];
#ifndef MC__ODEBNDS_GSL_USE_BAD
  delete[] _pADJTC; _pADJTC = _pDAG->FAD( 1, &pHAM, _npar, _pVAR+_nx+_nx );
#else
  delete[] _pADJTC; _pADJTC = _pDAG->BAD( 1, &pHAM, _npar, _pVAR+_nx+_nx );
#endif

  // Add initial state contribution to derivative bounds
  _pDAG->eval( _npar, _pADJTC, _Iyq, _nVAR-_npar, _pVAR, _IVAR, true );

  return true;
}

template <typename T, typename PMT, typename PVT> 
template <typename REALTYPE, typename OPT> inline bool
ODEBNDS_BASE<T,PMT,PVT>::_CC_I_ADJ
( const OPT&options, const double t, REALTYPE*vec, unsigned pos_fct, unsigned ifct )
{
  const FFVar* pFCT = _vFCT.at(pos_fct-1)+ifct;
#ifndef MC__ODEBNDS_GSL_USE_BAD
  delete[] _pADJTC; _pADJTC = _pDAG->FAD( 1, pFCT, _nz, _pVAR+_nx );
#else
  delete[] _pADJTC; _pADJTC = _pDAG->BAD( 1, pFCT, _nz, _pVAR+_nx );
#endif
  for( unsigned iy=0; iy<_nx; iy++ )
    _pADJCC[iy] = _pY[iy] + _pADJTC[iy];

  switch( options.WRAPMIT){
  case OPT::NONE:
  case OPT::DINEQ:
    _pDAG->eval( _nx, _pADJCC, _Ixdot, _nVAR-_npar, _pVAR, _IVAR );
    _I2vec( _nx, _Ixdot, vec );
    break;

  case OPT::ELLIPS:
  default:
    for( unsigned jx=0; jx<_nx; jx++ )
      _MVYZr[jx].set( _MVYZenv, _nx+jx, _Ir[jx] );
    for( unsigned jp=0; jp<_npar; jp++ )
      _MVYZp[jp].set( _MVYZenv, _nx+_nx+jp, _Ip[jp] );
    for( unsigned jy=0; jy<_nx; jy++ )
      _MVYZd[jy].set( _MVYZenv, jy, _Idy[jy] );
    *_MVYZt = -t; // current time
    _ep2x( _nx, _npar, _MVYZr, _pref, _MVYZp, _B, _zref, _MVYZz );
    _ep2x( _nx, _nx, _npar, _MVYZd, 0, _MVYZr, _By, _pref, _MVYZp, _By+_nx*_nx,
           _yref, _MVYZy );
    _opADJTC  = _pDAG->subgraph( _nx, _pADJCC );
    delete[] _PMIC; _PMIC = new PVT[_opADJTC.size()];
    _CC_I_ELL( _pDAG, _opADJTC, _PMIC, _nx, _pADJCC, _nVAR-_npar, _pVAR, _MVYZVAR,
               _MVYZf, _Edy, _Ay, _nz, _yrefdot, _Bydot, _Idydot, _Qydot,
               options.QTOL, machprec() );
    _E2vec( _nx, _nz, _yrefdot, _Qydot, _Bydot, vec );
    break;
  }

  return true;
}

template <typename T, typename PMT, typename PVT> 
template <typename REALTYPE, typename OPT> inline bool
ODEBNDS_BASE<T,PMT,PVT>::_CC_I_QUAD
( const OPT&options, const double t, REALTYPE*vec )
{
  *_It = -t; // current time
  _pDAG->eval( _npar, _pADJTC+_nx, _Iyq, _nVAR-_npar, _pVAR, _IVAR, true );
  _I2vec( _npar, _Iyq, vec );

  return true;
}

template <typename T, typename PMT, typename PVT>
template <typename OPT> inline bool
ODEBNDS_BASE<T,PMT,PVT>::_SET_I_ADJ
( const OPT &options, unsigned iADJRHS, unsigned iQUAD, unsigned pos_fct, unsigned ifct )
{
  if( _vRHS.size() <= iADJRHS ) return false; 
  if( _nq && _vQUAD.size() <= iQUAD ) return false;

  FFVar pHAM( 0. );
  _pRHS = _vRHS.at( iADJRHS );
  for( unsigned ix=0; ix<_nx; ix++ ) pHAM += _pVAR[ix] * _pRHS[ix];
  _pQUAD  = _nq? _vQUAD.at( iQUAD ): 0;
  const FFVar* pFCT = _vFCT.at(pos_fct)+ifct;
#ifndef MC__ODEBNDS_GSL_USE_BAD
  delete[] _pADJTC; _pADJTC = _nq? _pDAG->FAD( 1, pFCT, _nq, _pQ ): 0;
#else
  delete[] _pADJTC; _pADJTC = _nq? _pDAG->BAD( 1, pFCT, _nq, _pQ ): 0;
#endif
  for( unsigned iq=0; iq<_nq; iq++ ){
    if( !_pADJTC[iq].cst() ) return false; // quadrature appears nonlinearly in function
    pHAM += _pQUAD[iq] * _pADJTC[iq];
  }

  delete[] _pADJRHS;
#ifndef MC__ODEBNDS_GSL_USE_BAD
  _pADJRHS = _pDAG->FAD( 1, &pHAM, _nx+_npar, _pVAR+_nx );
#else
  _pADJRHS = _pDAG->BAD( 1, &pHAM, _nx+_npar, _pVAR+_nx );
#endif
  _opADJRHS = _pDAG->subgraph( _nx, _pADJRHS );
  _pADJQUAD = _pADJRHS + _nx;
  _opADJQUAD = _pDAG->subgraph(_npar, _pADJQUAD);
  const unsigned opmax = _opADJRHS.size()>_opADJQUAD.size()?_opADJRHS.size():_opADJQUAD.size();

  delete[] _IADJRHS;   _IADJRHS = 0;
  delete[] _opADJRHSi; _opADJRHSi = 0;
  delete[] _PMADJRHS;  _PMADJRHS = 0;
  switch( options.WRAPMIT){
  case OPT::NONE:
    _IADJRHS = new T[ opmax ];
    break;
  case OPT::DINEQ:
    _opADJRHSi = new std::list<const FFOp*>[_nx];
    for( unsigned ix=0; ix<_nx; ix++ )
      _opADJRHSi[ix] = _pDAG->subgraph( 1, _pADJRHS+ix );
    _IADJRHS = new T[ opmax ];
    break;
  case OPT::ELLIPS:
  default:
    _IADJRHS = new T[ _opADJQUAD.size() ];
    _PMADJRHS = new PVT[_opADJRHS.size()];
    break;
  }

  return true;
}

template <typename T, typename PMT, typename PVT>
template <typename REALTYPE, typename OPT> inline bool
ODEBNDS_BASE<T,PMT,PVT>::_RHS_I_ADJ
( const OPT &options, double t, const REALTYPE* y, REALTYPE* ydot,
  REALTYPE* vec_sta )
{
  if( !_pADJRHS ) return false; // **error** ADJRHS not defined

  switch( options.WRAPMIT){
  case OPT::NONE:
    _vec2I( vec_sta, _nx, _Iz );
    _vec2I( y, _nx, _Iy );// set current adjoint bounds
    *_It = -t; // set current time
    _pDAG->eval( _opADJRHS, _IADJRHS, _nx, _pADJRHS, _Iydot, _nVAR-_npar, _pVAR,
                 _IVAR );
    _I2vec( _nx, _Iydot, ydot );
    break;  
   
  case OPT::DINEQ:
    _vec2I( vec_sta, _nx, _Iz );
    _vec2I( y, _nx, _Iy );  // set current adjoint bounds
    *_It = -t; // set current time
    _RHS_I_DI( _pDAG, _opADJRHSi, _IADJRHS, _nx, _pADJRHS, _nVAR-_npar, _pVAR,
               _IVAR, _Iydot, _yLdot, _yUdot);
    _I2vec( _nx, _yLdot, _yUdot, ydot );
    break;  
   
  case OPT::ELLIPS:
  default:
    _vec2E( vec_sta, _nx, _npar, _Q, _Er, _Ir, _pref, _Ip, _B, _zref, _Iz );
    for( unsigned jx=0; jx<_nx; jx++ )
      _MVYZr[jx].set( _MVYZenv, _nx+jx, _Ir[jx] );
    for( unsigned jp=0; jp<_npar; jp++ )
      _MVYZp[jp].set( _MVYZenv, _nx+_nx+jp, _Ip[jp] );
    _ep2x( _nx, _npar, _MVYZr, _pref, _MVYZp, _B, _zref, _MVYZz );

    // Adjoint bounds
    _vec2E( y, _nx, _nx, _npar, _Qy, _Edy, _Idy, 0, _Ir, _By, _pref, _Ip,
            _By+_nx*_nx, _yref, _Iy );
    for( unsigned jy=0; jy<_nx; jy++ ) { _MVYZd[jy].set( _MVYZenv, jy, _Idy[jy] ); }
    _ep2x( _nx, _nx, _npar, _MVYZd, 0, _MVYZr, _By, _pref, _MVYZp, _By+_nx*_nx,
           _yref, _MVYZy );

    // Adjoint derivative bounds
    *_MVYZt = -t; // Current time
    _RHS_I_ELL( _pDAG, _opADJRHS, _PMADJRHS, _nx, _pADJRHS, _nVAR-_npar, 
                _pVAR, _MVYZVAR, _MVYZf, _Qy, _Ay, _nz, _yrefdot, _Bydot,
                _Idydot, _Qydot, options.QTOL, machprec() );
    _E2vec( _nx, _nz, _yrefdot, _Qydot, _Bydot, ydot );
    break;
  }
  return true;
}

template <typename T, typename PMT, typename PVT>
template <typename REALTYPE, typename OPT> inline bool
ODEBNDS_BASE<T,PMT,PVT>::_RHS_I_QUAD
( const OPT&options, double t, const REALTYPE*y, REALTYPE*qdot,
  REALTYPE* vec_sta, const bool bndinit )
{
  if( !_pADJQUAD ) return false;

  switch( options.WRAPMIT){
  case OPT::NONE:
  case OPT::DINEQ:
    if( !bndinit ) break;
    _vec2I( vec_sta, _nx, _Iz );
    *_It = -t; // set current time
    _vec2I( y, _nx, _Iy );// set current adjoint bounds
    break;

  case OPT::ELLIPS:
  default:
    if( !bndinit ) break;
    _vec2E( vec_sta, _nx, _npar, _Q, _Er, _Ir, _pref, _Ip, _B, _zref, _Iz );
    _vec2E( y, _nx, _nx, _npar, _Qy, _Edy, _Idy, 0, _Ir, _By, _pref, _Ip, 
            _By+_nx*_nx, _yref, _Iy ); // set current adj bounds
    *_It = -t; // set current time
    break;
  }

  _QUAD_I( _pDAG, _opADJQUAD, _IADJRHS, _npar, _pADJQUAD, _nVAR-_npar, _pVAR,
           _IVAR, _Iyqdot );
  _I2vec( _npar, _Iyqdot, qdot );
  return true;
}

template <typename T, typename PMT, typename PVT>
template <typename OPT> inline bool
ODEBNDS_BASE<T,PMT,PVT>::_INI_PM_ADJ
( const OPT &options, const PVT* PMp )
{
  // Check polynomial model compatibility and size
  // Not sure if this should be p or p-x (i.e. z)
  unsigned kp=_npar;
  for( unsigned ip=0; ip<_npar && kp==_npar; ip++ )
    if( PMp[ip].env() ) kp = ip;
  if( kp==_npar || PMp[kp].env()->nvar()!=_npar ) return false;
  _PMenv = PMp[kp].env();

  // Size and set DAG evaluation arrays
  _nz = _nx + _npar;
  _nVAR = _nx + _nz + 1 + _nq + _npar;
  delete[] _pVAR; _pVAR = new FFVar[_nVAR];
  delete[] _IVAR; _IVAR = new T[_nVAR];
  delete[] _PMVAR; _PMVAR = new PVT[_nVAR];
  delete[] _pADJCC; _pADJCC = new FFVar[_nx];
  BASE_DE::set_adjoint();
  for( unsigned ix=0; ix<_nx; ix++ ) _pVAR[ix] = _pY[ix];
  for( unsigned ix=0; ix<_nx; ix++ ) _pVAR[_nx + ix] = _pX[ix];
  for( unsigned ip=0; ip<_npar; ip++ ) _pVAR[_nx + _nx +ip] = _pP[ip];
  _pVAR[_nx + _nz] = (_pT? *_pT: 0. );
  _Iy = _Iz = _Ip = _Iyq = 0; // not used
  for( unsigned ip=0; ip<_npar; ip++ ) _PMVAR[_nx + _nx +ip] = PMp[ip];
  _PMy = _PMVAR;
  _PMz = _PMy + _nx;
  _PMp = _PMz + _nx;
  _PMyq = _PMp + _npar+1; // +1 for time
  
  // Reset _MVXPenv and related variables
  unsigned MVYZsize = ( options.ORDMIT<_PMenv->nord()? options.ORDMIT: _PMenv->nord() ); 
  if( _MVYZenv && ( _MVYZenv->nord() != MVYZsize || _MVYZenv->nvar() != _nx+_nz ) ){
    delete[] _MVYZf;    _MVYZf = 0;
    delete[] _MVYZdgdy; _MVYZdgdy = 0;
    delete[] _MVYZd;    _MVYZd = 0;
    delete[] _MVYZVAR;  _MVYZVAR = _MVYZy = _MVYZz = 0; 
    delete   _MVYZenv;  _MVYZenv = 0;
  }

  // Set parameterization variables
  delete[] _PMydot; _PMydot = new PVT[_nx];
  delete[] _PMyqdot; _PMyqdot = _npar? new PVT[_npar]: 0;
  delete[] _Ryqdot;  _Ryqdot  = _npar? new double[_npar]: 0;
  switch( options.WRAPMIT){
  case OPT::NONE:
    break;
  case OPT::DINEQ:
    delete[] _RyLdot; _RyLdot = new double[_nx];
    delete[] _RyUdot; _RyUdot = new double[_nx];
    break;
  case OPT::ELLIPS:
  default:
    if( !_yref )     _yref      = new double[_nx];
    if( !_Ay )       _Ay        = new double[_nx*_nx];
    if( !_Qy )       _Qy        = new double[_nx*(_nx+1)/2];
    if( !_Idy )      _Idy       = new T[_nx];
    if( !_Qydot )    _Qydot     = new double[_nx*(_nx+1)/2];
    if( !_Idydot )   _Idydot    = new T[_nx];
    if( !_Idgdy )    _Idgdy     = new T[_nx*_nx];
    if( !_MVYZenv ) _MVYZenv    = new PMT( _nx+_nz, MVYZsize );
    _MVYZenv->options = _PMenv->options;
    if( !_MVYZd )    _MVYZd     = new PVT[_nx];
    if( !_MVYZr )    _MVYZr     = new PVT[_nx];
    if( !_MVYZf )    _MVYZf     = new PVT[_nx];
    if( !_MVYZdgdy ) _MVYZdgdy  = new PVT[_nx*_nx];
    if( !_MVYZVAR )  _MVYZVAR   = new PVT[_nVAR];
    _MVYZy = _MVYZVAR;
    _MVYZz = _MVYZy + _nx;
    _MVYZp = _MVYZz + _nx;
    for( unsigned ip=0; ip<_npar; ip++ ){ _MVYZp[ip].set( _MVYZenv, ip, _PMp[ip].B() ); }
    _MVYZt = _MVYZp + _npar;
    if( !_MVZd )     _MVZd      = new PVT[_nx];
    if( !_MVZy )     _MVZy      = new PVT[_nx];
    if( !_MVZVAR )   _MVZVAR    = new PVT[_nz+1];
    _MVZz = _MVZVAR;
    _MVZp = _MVZz + _nx;
    break;
  }

  return true;
}

template <typename T, typename PMT, typename PVT>
template <typename REALTYPE, typename OPT> inline bool
ODEBNDS_BASE<T,PMT,PVT>::_TC_PM_ADJ
( const OPT &options, const double t, REALTYPE*vec, unsigned pos_fct, unsigned ifct )
{
  const FFVar* pFCT = _vFCT.at(pos_fct)+ifct;
#ifndef MC__ODEBNDS_GSL_USE_BAD
  delete[] _pADJTC; _pADJTC = _pDAG->FAD( 1, pFCT, _nz, _pVAR+_nx );
#else
  delete[] _pADJTC; _pADJTC = _pDAG->BAD( 1, pFCT, _nz, _pVAR+_nx );
#endif
  _opADJTC  = _pDAG->subgraph( _nx, _pADJTC );
  delete[] _PMADJTC; _PMADJTC = new PVT[_opADJTC.size()];

  switch( options.WRAPMIT){
  case OPT::NONE:
  case OPT::DINEQ:
    _pDAG->eval( _opADJTC, _PMADJTC, _nx, _pADJTC, _PMy, _nz+1, _pVAR+_nx, _PMVAR+_nx );
    _PMI2vec( _PMenv, _nx, _PMy, vec );
    //??_PMI2vec( _PMenv, _nx, _PMy, _npar, _PMyq, _vec_adj+_pos_adj );
    break;

  case OPT::ELLIPS:
  default:
    _IC_PM_ELL( _pDAG, _opADJTC, _nx, _pADJTC, _nz+1, _pVAR+_nx, _PMz, _PMy, _Qy, _Edy, _Idy );
    _PME2vec( _PMenv, _nx, _PMy, _Qy, vec );
    //??_PME2vec( _PMenv, _nx, _PMy, _Qy, _npar, _PMyq, _vec_adj+_pos_adj );
    break;
  }

  return true;
}

template <typename T, typename PMT, typename PVT>
template <typename REALTYPE, typename OPT> inline bool
ODEBNDS_BASE<T,PMT,PVT>::_TC_PM_QUAD
( const OPT &options, REALTYPE*vec )
{
  if( !_vQUAD.size() || !_nq ) return true;
  _pDAG->eval( _npar, _pADJTC+_nx, _PMyq, _nz+1, _pVAR+_nx, _PMVAR+_nx );
  _PMI2vec( _PMenv, _npar, _PMyq, vec, true );
  return true;
}

template <typename T, typename PMT, typename PVT>
template <typename REALTYPE, typename OPT> inline bool
ODEBNDS_BASE<T,PMT,PVT>::_CC_PM_ADJ
( const OPT &options, const double t, REALTYPE*vec, unsigned pos_fct, unsigned ifct )
{
  const FFVar* pFCT = _vFCT.at(pos_fct-1)+ifct;
#ifndef MC__ODEBNDS_GSL_USE_BAD
  delete[] _pADJTC; _pADJTC = _pDAG->FAD( 1, pFCT, _nz, _pVAR+_nx );
#else
  delete[] _pADJTC; _pADJTC = _pDAG->BAD( 1, pFCT, _nz, _pVAR+_nx );
#endif
  for( unsigned iy=0; iy<_nx; iy++ ){ _pADJCC[iy] = _pY[iy] + _pADJTC[iy]; }
  _opADJTC = _pDAG->subgraph( _nx, _pADJCC );
  _opADJDTC.clear();
  delete[] _PMADJTC;  _PMADJTC  = 0;
  delete[] _PMADJDTC; _PMADJDTC = 0;
  delete[] _IADJDTC;  _IADJDTC  = 0;

  switch( options.WRAPMIT){
  case OPT::NONE:
  case OPT::DINEQ:
    _PMADJTC = new PVT[_opADJTC.size()];
    _pDAG->eval( _opADJTC, _PMADJTC, _nx, _pADJCC, _PMxdot, _nVAR-_npar, _pVAR, _PMVAR );
    _pDAG->eval( _npar, _pADJTC+_nx, _PMyq, _nVAR-_npar, _pVAR, _PMVAR, true );
    _PMI2vec( _PMenv, _nx, _PMxdot, vec );
    //??_PMI2vec( _PMenv, _nx, _PMxdot, _npar, _PMyq, _vec_adj+_pos_adj );
    break;

  case OPT::ELLIPS:
  default:
    for( unsigned jx=0; jx<_nx; jx++ ){ _MVYZr[jx].set( _MVYZenv, _nx+jx, _Ir[jx] ); }
    _vec2PME( vec, _PMenv, _nx, _PMy, _Qy, _Edy, _Idy );
    //??_vec2PME( _vec_adj+_pos_adj, _PMenv, _nx, _PMy, _Qy, _Edy, _Idy, _npar, _PMyq ); // current adjoint/quadrature polynomial model
    _PMADJTC = new PVT[_opADJTC.size()];

    // In this variant a bound on the Jacobian matrix is computed and the
    // linear part is taken as the mid-point of this matrix
    if( !options.ORDMIT ){
      for( unsigned ix=0; ix<_nx; ix++ ){
        _IVAR[ix] = _PMy[ix].bound(); // set current adjoint bounds
        ( _PMVAR[ix].center() ).set( T(0.) ); // cancel remainder term
      }
      _pADJDTC = _pDAG->FAD( _nx, _pADJCC, _nx, _pVAR ); // This is the identity matrix!!
      _opADJDTC = _pDAG->subgraph( _nx*_nx, _pADJDTC );
      _IADJDTC = new T[ _opADJDTC.size() ];
      _RHS_PM_ELL0( _pDAG, _opADJTC, _PMIC, _opADJDTC, _IADJDTC, _nx, _pADJCC,
                    _pADJDTC, _nVAR-_npar, _pVAR, _PMVAR, _IVAR, _PMxdot, _Idgdy,
                    _Idy, _Ay, _Idydot );
    }
    // In this variant a polynomial model of the Jacobian matrix is computed and the
    // linear part is taken as the mid-point of this matrix
    else if( _PMenv->nord() > _MVYZenv->nord() ){
      for( unsigned jx=0; jx<_nx; jx++ ){
        _MVYZd[jx].set( _MVYZenv, _npar+jx, _Idy[jx] );
        _MVYZy[jx].set( _MVYZenv ).set( _PMy[jx].center().set( T(0.) ), true );
      }
      _e2x( _nx, _MVYZd, _MVYZy, false );
      *_MVYZt = t; // set current time
      _pADJDTC = _pDAG->FAD( _nx, _pADJCC, _nx, _pVAR );
      _opADJDTC = _pDAG->subgraph( _nx*_nx, _pADJDTC );
      _PMADJDTC = new PVT[ _opADJDTC.size() ];

      _RHS_PM_ELL1( _pDAG, _opADJTC, _PMIC, _opADJDTC, _PMADJDTC, _nx, _pADJCC,
                    _pADJDTC, _nVAR-_npar, _pVAR, _PMVAR, _MVYZVAR, _PMxdot,
                    _MVYZdgdy, _Idy, _Ay, _Idydot );
    }
    // In this variant a polynomial model in the joint adjoint-parameter and
    // of the same order as the parameter polynomial model is computed
    else{

      for( unsigned jx=0; jx<_nx; jx++ ){
        _MVYZd[jx].set( _MVYZenv, _npar+jx, _Idy[jx] );
        _MVYZy[jx].set( _MVYZenv ).set( _PMy[jx].center().set( T(0.) ), true );
      }
      _e2x( _nx, _MVYZd, _MVYZy, false );
      *_MVYZt = t; // set current time   
      _RHS_PM_ELL2( _pDAG, _opADJTC, _PMIC, _nx, _pADJRHS, _nVAR-_npar, _pVAR,
                    _MVYZVAR, _PMenv, _PMxdot, _MVYZf, _npar, _Idy, _Ay, _Idydot );
    }

    _CC_PM_ELL( _nx, _Edy, _Ay, _Idydot, _Qydot, options.QTOL, machprec() );
    _pDAG->eval( _npar, _pADJTC+_nx, _PMyq, _nVAR-_npar, _pVAR, _PMVAR, true );
    //_RHS_PM_ELL( _nx, _Qy, _Ay, _Idydot, _Qydot, options.QTOL, machprec() );
    _PME2vec( _PMenv, _nx, _PMxdot, _Qydot, vec );
    //??_PME2vec( _PMenv, _nx, _PMxdot, _Qydot, _npar, _PMyq, _vec_adj+_pos_adj ); 
    break;
  }

  return true;
}

template <typename T, typename PMT, typename PVT> inline bool
ODEBNDS_BASE<T,PMT,PVT>::_IC_PM_ADJ
()
{
  const FFVar* pIC = _vIC.at(0);
  FFVar pHAM( 0. );
  for( unsigned ix=0; ix<_nx; ix++ ) pHAM += _pVAR[ix] * pIC[ix];
#ifndef MC__ODEBNDS_GSL_USE_BAD
  delete[] _pADJTC; _pADJTC = _pDAG->FAD( 1, &pHAM, _npar, _pVAR+_nx+_nx );
#else
  delete[] _pADJTC; _pADJTC = _pDAG->BAD( 1, &pHAM, _npar, _pVAR+_nx+_nx );
#endif

  // Add initial state contribution to derivative bounds
  _pDAG->eval( _npar, _pADJTC, _PMyq, _nVAR-_npar, _pVAR, _PMVAR, true );

  return true;
}

template <typename T, typename PMT, typename PVT>
template <typename OPT> inline bool
ODEBNDS_BASE<T,PMT,PVT>::_SET_PM_ADJ
( const OPT &options, unsigned iADJRHS, unsigned iQUAD, unsigned pos_fct, unsigned ifct )
{
  if( _vRHS.size() <= iADJRHS ) return false; 
  if( _nq && _vQUAD.size() <= iQUAD ) return false;

  FFVar pHAM( 0. );
  _pRHS = _vRHS.at( iADJRHS );
  for( unsigned ix=0; ix<_nx; ix++ ) pHAM += _pVAR[ix] * _pRHS[ix];
  _pQUAD  = _nq? _vQUAD.at( iQUAD ): 0;
  const FFVar* pFCT = _vFCT.at(pos_fct)+ifct;
#ifndef MC__ODEBNDS_GSL_USE_BAD
  delete[] _pADJTC; _pADJTC = _nq? _pDAG->FAD( 1, pFCT, _nq, _pQ ): 0;
#else
  delete[] _pADJTC; _pADJTC = _nq? _pDAG->BAD( 1, pFCT, _nq, _pQ ): 0;
#endif
#ifdef MC__ODEBNDS_GSL_DINEQI_DEBUG
  for( unsigned iq=0; iq<_nq; iq++ )
    std::cout << "  dfdq(" << _ifct << "," << iq << "):" << _pADJTC[iq] << std::endl;
#endif

  for( unsigned iq=0; iq<_nq; iq++ ){
    if( !_pADJTC[iq].cst() ) return false; // quadrature appears nonlinearly in function
    pHAM += _pQUAD[iq] * _pADJTC[iq];
  }

  delete[] _pADJRHS;
#ifndef MC__ODEBNDS_GSL_USE_BAD
  _pADJRHS = _pDAG->FAD( 1, &pHAM, _nx+_npar, _pVAR+_nx );
#else
  _pADJRHS = _pDAG->BAD( 1, &pHAM, _nx+_npar, _pVAR+_nx );
#endif
  _opADJRHS = _pDAG->subgraph( _nx, _pADJRHS );
  _pADJQUAD = _pADJRHS + _nx;
  _opADJQUAD = _pDAG->subgraph(_npar, _pADJQUAD);
  const unsigned opmax = _opADJRHS.size()>_opADJQUAD.size()?_opADJRHS.size():_opADJQUAD.size();

  delete[] _IADJRHS;   _IADJRHS = 0;
  delete[] _PMADJRHS;  _PMADJRHS = 0;
  delete[] _opADJRHSi; _opADJRHSi = 0;
  delete[] _IADJJAC;   _IADJJAC = 0;
  delete[] _PMADJJAC;  _PMADJJAC = 0;
  delete[] _pADJJAC;   _pADJJAC = 0; _opADJJAC.clear();

  _IADJRHS = new T[opmax];
  _PMADJRHS = new PVT[opmax];

  switch( options.WRAPMIT){
  case OPT::NONE:
    break;
  case OPT::DINEQ:
    _opADJRHSi = new std::list<const FFOp*>[_nx];
    for( unsigned ix=0; ix<_nx; ix++ )
      _opADJRHSi[ix] = _pDAG->subgraph( 1, _pADJRHS+ix );
    break;
  case OPT::ELLIPS:
  default:
    _pADJJAC = _pDAG->FAD( _nx, _pADJRHS, _nx, _pVAR );
    _opADJJAC = _pDAG->subgraph( _nx*_nx, _pADJJAC );
    if( !options.ORDMIT )
      _IADJJAC = new T[ _opADJJAC.size() ];
    else if( _PMenv->nord() > _MVYZenv->nord() )
      _PMADJJAC = new PVT[ _opADJJAC.size() ];
    break;
  }

  return true;
}

template <typename T, typename PMT, typename PVT>
template <typename REALTYPE, typename OPT> inline bool
ODEBNDS_BASE<T,PMT,PVT>::_RHS_PM_ADJ
( const OPT &options, double t, const REALTYPE* y, REALTYPE* ydot, REALTYPE* vec_sta )
{
  if( !_pADJRHS ) return false;

  // interpolated state NOW SET IN EITHER GSL OR SUNDIALS

  switch( options.WRAPMIT){
  case OPT::NONE:

    _vec2PMI( vec_sta, _PMenv, _nx, _PMz );
    _vec2PMI( y, _PMenv, _nx, _PMy );// set current adjoint bounds
    _QUAD_PM( _pDAG, _opADJQUAD, _PMADJRHS, _npar, _pADJQUAD, _nVAR-_npar, _pVAR,
              _PMVAR, _PMyqdot );
    _pDAG->eval( _opADJRHS, _PMADJRHS, _nx, _pADJRHS, _PMydot, _nVAR-_npar, _pVAR,
                 _PMVAR );
    _PMI2vec( _PMenv, _nx, _PMydot, ydot );
    //_PMI2vec( _PMenv, _npar, _PMyqdot, ydot+_offset_quad);
    //??_PMI2vec( _PMenv, _nx, _PMydot, _npar, _PMyqdot, ydot );
    break;  
  case OPT::DINEQ:
    _vec2PMI( vec_sta, _PMenv, _nx, _PMz );
    _vec2PMI( y, _PMenv, _nx, _PMy );// set current adjoint bounds
    _RHS_PM_DI( _pDAG, _opADJRHSi, _PMADJRHS, _nx, _pADJRHS, _nVAR-_npar, _pVAR,
                _PMVAR, _PMydot, _RyLdot, _RyUdot);
    _PMI2vec( _PMenv, _nx, _PMydot, _RyLdot, _RyUdot, ydot );
    //??_PMI2vec( _PMenv, _nx, _PMydot, _RyLdot, _RyUdot, _npar, _PMyqdot, ydot );
    break;  
  case OPT::ELLIPS:
  default:
    _vec2PME( vec_sta, _PMenv, _nx, _PMz, _Q, _Er, _Ir );
    for( unsigned jx=0; jx<_nx; jx++ ){ _MVYZr[jx].set( _MVYZenv, _nx+jx, _Ir[jx] ); }
    
    _vec2PME( y, _PMenv, _nx, _PMy, _Qy, _Edy, _Idy ); // set current adjoint polynomial model

    // Adjoint derivative bounds
    _QUAD_PM( _pDAG, _opADJQUAD, _PMADJRHS, _npar, _pADJQUAD, _nVAR-_npar, _pVAR, _PMVAR,
              _PMyqdot );

    // In this variant a bound on the Jacobian matrix is computed and the
    // linear part is taken as the mid-point of this matrix
    if( !options.ORDMIT ){
      for( unsigned ix=0; ix<_nx; ix++ ){
        _IVAR[ix] = _PMy[ix].bound(); // set current adjoint bounds
        ( _PMVAR[ix].center() ).set( T(0.) ); // cancel remainder term
      }
      _RHS_PM_ELL0( _pDAG, _opADJRHS, _PMADJRHS, _opADJJAC, _IADJJAC, _nx, _pADJRHS,
                    _pADJJAC, _nVAR-_npar, _pVAR, _PMVAR, _IVAR, _PMydot, _Idgdy,
                    _Idy, _Ay, _Idydot );
    }

    // In this variant a polynomial model of the Jacobian matrix is computed and the
    // linear part is taken as the mid-point of this matrix
    else if( _PMenv->nord() > _MVYZenv->nord() ){
      for( unsigned jx=0; jx<_nx; jx++ ){
        _MVYZd[jx].set( _MVYZenv, _npar+jx, _Idy[jx] );
        _MVYZy[jx].set( _MVYZenv ).set( _PMy[jx].center().set( T(0.) ), true );
      }
      _e2x( _nx, _MVYZd, _MVYZy, false );
      *_MVYZt = t; // set current time
      //std::cout << "_RHS_PM_ELL1" << std::endl;
      _RHS_PM_ELL1( _pDAG, _opADJRHS, _PMADJRHS, _opADJJAC, _PMADJJAC, _nx, _pADJRHS,
                    _pADJJAC, _nVAR-_npar, _pVAR, _PMVAR, _MVYZVAR, _PMydot,
                    _MVYZdgdy, _Idy, _Ay, _Idydot );
    }

    // In this variant a polynomial model in the joint adjoint-parameter and
    // of the same order as the parameter polynomial model is computed
    else{
      for( unsigned jx=0; jx<_nx; jx++ ){
        _MVYZd[jx].set( _MVYZenv, _npar+jx, _Idy[jx] );
        _MVYZy[jx].set( _MVYZenv ).set( _PMy[jx].center().set( T(0.) ), true );
      }
      _e2x( _nx, _MVYZd, _MVYZy, false );
      *_MVYZt = t; // set current time 
      _RHS_PM_ELL2( _pDAG, _opADJRHS, _PMADJRHS, _nx, _pADJRHS, _nVAR-_npar, _pVAR,
                    _MVYZVAR, _PMenv, _PMydot, _MVYZf, _npar, _Idy, _Ay, _Idydot );
    }

    _RHS_PM_ELL( _nx, _Qy, _Ay, _Idydot, _Qydot, options.QTOL, machprec() );
    _PME2vec( _PMenv, _nx, _PMydot, _Qydot, ydot);
    //??_PME2vec( _PMenv, _nx, _PMydot, _Qydot, _npar, _PMyqdot, ydot );
    break;
  }

  return true;
}

template <typename T, typename PMT, typename PVT>
template <typename REALTYPE, typename OPT> inline bool
ODEBNDS_BASE<T,PMT,PVT>::_RHS_PM_QUAD
( const OPT&options, double t, const REALTYPE*y, REALTYPE*qdot, const bool reinit )
{
  if( !_pQUAD ) return false;

  switch( options.WRAPMIT){
  case OPT::NONE:
  case OPT::DINEQ:
    if( !reinit ) break;
    _vec2PMI( y, _PMenv, _nx, _PMy );// set current adjoint bounds
    break;
  case OPT::ELLIPS:
  default:
    if( !reinit ) break;
    *_MVYZt = t; // set current time 
    _vec2PME( y, _PMenv, _nx, _PMy, _Qy, _Edy, _Idy ); // set current adjoint polynomial model
    break;
  }

  _QUAD_PM( _pDAG, _opADJQUAD, _PMADJRHS, _npar, _pADJQUAD, _nVAR-_npar, _pVAR,
              _PMVAR, _PMyqdot );

  // Whether or not to ignore the remainder
  if( !options.PMNOREM ) {_PMI2vec( _PMenv, _npar, _PMyqdot, qdot, true);}
  else {_PMI2vec( _PMenv, _npar, _PMyqdot, 0, qdot );}

  return true;  
}

} // end namescape mc

#endif

















