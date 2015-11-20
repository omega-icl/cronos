// Copyright (C) 2015 Benoit Chachuat & Nikola Peric, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

#ifndef MC__ODEBNDS_BASE_HPP
#define MC__ODEBNDS_BASE_HPP

#undef  MC__ODEBNDS_BASE_DINEQI_DEBUG
#undef  MC__ODEBNDS_BASE_DINEQPM_DEBUG
#undef  MC__ODEBNDS_BASE_MVYZ_USE

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
  using ODEBND_BASE<T,PMT,PVT>::_RHS_PM_ELL0;
  using ODEBND_BASE<T,PMT,PVT>::_RHS_PM_ELL1;
  using ODEBND_BASE<T,PMT,PVT>::_RHS_PM_ELL2;
  using ODEBND_BASE<T,PMT,PVT>::_RHS_PM_ELL;
  using ODEBND_BASE<T,PMT,PVT>::_e2x;
  using ODEBND_BASE<T,PMT,PVT>::_print_interm;

public:
  //! @brief Default constructor
  ODEBNDS_BASE();

  //! @brief Virtual destructor
  virtual ~ODEBNDS_BASE();

protected:
  //! @brief array of list of operations in adjoint RHS evaluation
  std::list<const FFOp*>* _opADJRHS;

  //! @brief array of list of operations in adjoint RHS Jacobian
  std::list<const FFOp*> *_opADJJAC;

  //! @brief array of list of operations in adjoint quadrature evaluation
  std::list<const FFOp*>* _opADJQUAD;

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
  PVT *_MVYZg;

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
    ( const double t );

  //! @brief Function to reinitialize adjoint bounds after discontinuity
  template <typename REALTYPE, typename OPT> bool _CC_I_ADJ
    ( const OPT &options, const double t, REALTYPE *vec,
      unsigned pos_fct, unsigned ifct );

  //! @brief Function to reinitialize adjoint quadrature bounds after discontinuity
  template <typename REALTYPE, typename OPT> bool _CC_I_QUAD
    ( const OPT&options, const double t, REALTYPE*vec );

  //! @brief Set adjoint RHS pointer and corresponding Jacobian
  template <typename OPT> bool _SET_I_ADJ
    ( const OPT &options, unsigned iADJRHS, unsigned iQUAD,
      unsigned pos_fct, const bool neg=true );

  //! @brief Function to calculate the adjoint ODEs RHS values in interval arithmetic
  template <typename REALTYPE, typename OPT> bool _RHS_I_ADJ
  (  const OPT &options, double t, const REALTYPE* y, REALTYPE* ydot,
     REALTYPE* x, const unsigned ifct, const bool neg=false  );

  //! @brief Function to calculate the RHS of auxiliary ODEs in interval arithmetic
  template <typename REALTYPE, typename OPT> bool _RHS_I_QUAD
    ( const OPT&options, double t, const REALTYPE*y, REALTYPE*qdot,
      REALTYPE* x, const unsigned ifct, const bool bndinit=true );

  //! @brief Function to initialize adjoint polynomial models
  template <typename OPT> bool _INI_PM_ADJ
    ( const OPT &options, const unsigned np, const PVT *PMp );

  //! @brief Function to initialize adjoint polynomial model
  template <typename REALTYPE, typename OPT> bool _TC_PM_ADJ
    ( const OPT &options, const double t, REALTYPE *vec,  unsigned pos_fct, unsigned ifct );

  //! @brief Function to initialize quadrature polynomial model
  template <typename REALTYPE, typename OPT> bool _TC_PM_QUAD
    ( const OPT &options, const double t, REALTYPE *vec );

  //! @brief Function to add initial state contribution to function derivatives
  bool _IC_PM_ADJ
    ( const double t );

  //! @brief Function to reinitialize adjoint polynomial bounds after discontinuity
  template <typename REALTYPE, typename OPT> bool _CC_PM_ADJ
    ( const OPT &options, const double t, REALTYPE *vec, unsigned pos_fct, unsigned ifct );

  //! @brief Function to reinitialize adjoint quadrature polynomial bounds after discontinuity
  template <typename REALTYPE, typename OPT> bool _CC_PM_QUAD
    ( const OPT&options, const double t, REALTYPE*vec );

  //! @brief Function to set adjoint RHS pointer and corresponding polynomial
  template <typename OPT> bool _SET_PM_ADJ
    ( const OPT &options, unsigned iADJRHS, unsigned iQUAD,
      unsigned pos_fct, const bool neg=true );

  //! @brief Function to calculate the RHS of adjoint ODEs in polynomial model arithmetic
  template <typename REALTYPE, typename OPT> bool _RHS_PM_ADJ
    ( const OPT &options, double t, const REALTYPE* y, REALTYPE* ydot,
      REALTYPE* x, const unsigned ifct, const bool neg=false );

  //! @brief Function to calculate the RHS of auxiliary ODEs in polynomial model arithmetic
  template <typename REALTYPE, typename OPT> bool _RHS_PM_QUAD
    ( const OPT&options, double t, const REALTYPE*y, REALTYPE*qdot,
      REALTYPE* x, const unsigned ifct, const bool bndinit=true );

  //! @brief Private methods to block default compiler methods
  ODEBNDS_BASE(const ODEBNDS_BASE&);
  ODEBNDS_BASE& operator=(const ODEBNDS_BASE&);
};

template <typename T, typename PMT, typename PVT> inline
ODEBNDS_BASE<T,PMT,PVT>::ODEBNDS_BASE
()
: BASE_DE(), ODEBND_BASE<T,PMT,PVT>(), _pADJRHS(0), _pADJJAC(0),
  _pADJQUAD(0), _pADJTC(0), _pADJCC(0), _pADJDTC(0)
{
  // Initialize adjoint arrays
  _opADJRHS = _opADJQUAD = _opADJJAC = 0;
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
  _MVYZd = _MVYZr = _MVYZg = _MVYZdgdy = _MVYZVAR = _MVYZy = _MVYZz = _MVYZp = _MVYZt = 0;
  _MVZd  = _MVZy  = _MVZVAR  = _MVZz  = _MVZp = _MVZt = 0;
}

template <typename T, typename PMT, typename PVT> inline
ODEBNDS_BASE<T,PMT,PVT>::~ODEBNDS_BASE
()
{
  delete[] _opADJRHS;
  delete[] _opADJQUAD;
  delete[] _opADJJAC;
  delete[] _IADJRHS;
  delete[] _PMADJRHS;
  delete[] _IADJJAC;
  delete[] _PMADJJAC;
  delete[] _pADJJAC;
  delete[] _pADJRHS;
  delete[] _pADJQUAD;
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
  delete[] _MVYZg;
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
    delete[] _MVYZg;   _MVYZg = 0;
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
    delete[] _zref;    _zref    = new double[_nz];
    for( unsigned ip=0; ip<_npar; ip++ )
      _pref[ip] = _zref[_nx+ip] = Op<T>::mid( _Ip[ip] );
    delete[] _yref;    _yref    = new double[_nx];
    delete[] _Ay;      _Ay      = new double[_nx*_nx];
    delete[] _By;      _By      = new double[_nx*_nz];
    delete[] _Qy;      _Qy      = new double[_nx*(_nx+1)/2];
    delete[] _Idy;     _Idy     = new T[_nx];
    delete[] _yrefdot; _yrefdot = new double[_nx];
    delete[] _Bydot;   _Bydot   = new double[_nx*_nz];
    delete[] _Qydot;   _Qydot   = new double[_nx*(_nx+1)/2];
    delete[] _Idydot;  _Idydot  = new T[_nx];

    delete   _MVYZenv; _MVYZenv = new PMT( _nx+_nz, options.ORDMIT );
    _MVYZenv->options = options.PMOPT;
    delete[] _MVYZd;   _MVYZd   = new PVT[_nx];
    delete[] _MVYZr;   _MVYZr   = new PVT[_nx];
    delete[] _MVYZg;   _MVYZg   = new PVT[_nx];
    delete[] _MVYZVAR; _MVYZVAR = new PVT[_nVAR];
    _MVYZy = _MVYZVAR;
    _MVYZz = _MVYZy + _nx;
    _MVYZp = _MVYZz + _nx;
    _MVYZt = _MVYZp + _npar;
    for( unsigned jp=0; jp<_npar; jp++ )
      _MVYZp[jp].set( _MVYZenv, _nx+_nx+jp, _Ip[jp] );

    delete   _MVZenv;  _MVZenv  = new PMT( _nz, 1 );
    _MVZenv->options = options.PMOPT;
    delete[] _MVZd;    _MVZd    = new PVT[_nx];
    delete[] _MVZy;    _MVZy    = new PVT[_nx];
    delete[] _MVZVAR;  _MVZVAR  = new PVT[_nz+1];
    _MVZz = _MVZVAR;
    _MVZp = _MVZz + _nx;
    _MVZt = _MVZp + _npar;
    for( unsigned ip=0; ip<_npar; ip++ )
      _MVZp[ip].set( _MVZenv, _nx+ip, _Ip[ip] );
    break;
  }

  return true;
}

template <typename T, typename PMT, typename PVT>
template <typename REALTYPE, typename OPT> inline bool
ODEBNDS_BASE<T,PMT,PVT>::_TC_I_ADJ
( const OPT&options, const double t, REALTYPE*y,
  unsigned pos_fct, unsigned ifct )
{
  const FFVar* pFCT = _vFCT.at(pos_fct)+ifct;
#ifndef MC__ODEBNDS_BASE_USE_BAD
  delete[] _pADJTC; _pADJTC = _pDAG->FAD( 1, pFCT, _nz, _pVAR+_nx );
#else
  delete[] _pADJTC; _pADJTC = _pDAG->BAD( 1, pFCT, _nz, _pVAR+_nx );
#endif
  _opADJTC  = _pDAG->subgraph( _nx, _pADJTC );
  *_It = t; // current time

  switch( options.WRAPMIT){

  case OPT::NONE:
  case OPT::DINEQ:
    _pDAG->eval( _opADJTC, _nx, _pADJTC, _Iy, _nz+1, _pVAR+_nx, _IVAR+_nx );
    _I2vec( _nx, _Iy, y );
    break;

  case OPT::ELLIPS:
  default:
    *_MVZt = t;
    for( unsigned jx=0; jx<_nx; jx++ )
      _MVZd[jx].set( _MVZenv, jx, _Ir[jx] );
    _ep2x( _nx, _npar, _MVZd, _pref, _MVZp, _B, _zref, _MVZz );
    _pDAG->eval( _opADJTC, _nx, _pADJTC, _MVZy, _nz+1, _pVAR+_nx, _MVZVAR+_nx );
    _IC_I_ELL( _nx, _MVZy, _yref, _Qy, _nz, _By, _Iy );
    _E2vec( _nx, _nz, _yref, _Qy, _By, y );
    break;
  }

  return true;
}

template <typename T, typename PMT, typename PVT>
template <typename REALTYPE, typename OPT> inline bool
ODEBNDS_BASE<T,PMT,PVT>::_TC_I_QUAD
( const OPT&options, const double t, REALTYPE*yq )
{
  *_It = t; // current time
  _opTCQUAD = _pDAG->subgraph( _npar, _pADJTC+_nx );
  _pDAG->eval( _opTCQUAD, _npar, _pADJTC+_nx, _Iyq, _nz+1, _pVAR+_nx, _IVAR+_nx );
  _I2vec( _npar, _Iyq, yq );
  return true;
}

template <typename T, typename PMT, typename PVT> inline bool
ODEBNDS_BASE<T,PMT,PVT>::_IC_I_ADJ
( const double t )
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
  *_It = t; // current time
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
    *_It = t; // current time
    _pDAG->eval( _nx, _pADJCC, _Ixdot, _nVAR-_npar, _pVAR, _IVAR );
    _I2vec( _nx, _Ixdot, vec );
    break;

  case OPT::ELLIPS:
  default:
    for( unsigned jx=0; jx<_nx; jx++ )
      _MVYZr[jx].set( _MVYZenv, _nx+jx, _Ir[jx] );
    for( unsigned jy=0; jy<_nx; jy++ )
      _MVYZd[jy].set( _MVYZenv, jy, _Idy[jy] );
    *_MVYZt = t; // current time
    _ep2x( _nx, _npar, _MVYZr, _pref, _MVYZp, _B, _zref, _MVYZz );
    _ep2x( _nx, _nx, _npar, _MVYZd, 0, _MVYZr, _By, _pref, _MVYZp, _By+_nx*_nx,
           _yref, _MVYZy );
    _opADJTC  = _pDAG->subgraph( _nx, _pADJCC );
    delete[] _PMIC; _PMIC = new PVT[_opADJTC.size()];
    _pDAG->eval( _opADJTC, _PMIC, _nx, _pADJCC, _MVYZg, _nVAR-_npar, _pVAR, _MVYZVAR );
    _CC_I_ELL( _nx, _MVYZg, _Edy, _Ay, _nz, _yrefdot, _Bydot, _Idydot, _Qydot,
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
  *_It = t; // current time
  //_pDAG->eval( _npar, _pADJTC+_nx, _Iyq, _nVAR-_npar, _pVAR, _IVAR, true );
  _pDAG->eval( _npar, _pADJTC+_nx, _Iyq, _nVAR-_npar, _pVAR, _IVAR, true );
  _I2vec( _npar, _Iyq, vec );

  return true;
}

template <typename T, typename PMT, typename PVT>
template <typename OPT> inline bool
ODEBNDS_BASE<T,PMT,PVT>::_SET_I_ADJ
( const OPT &options, unsigned iADJRHS, unsigned iQUAD, unsigned pos_fct,
  const bool neg )
{
  if( _vRHS.size() <= iADJRHS ) return false; 
  if( _nq && _vQUAD.size() <= iQUAD ) return false;

  _pRHS = _vRHS.at( iADJRHS );
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
    delete[] _pADJTC; _pADJTC = _nq? _pDAG->FAD( 1, pFCT+ifct, _nq, _pQ ): 0;
#else
    delete[] _pADJTC; _pADJTC = _nq? _pDAG->BAD( 1, pFCT+ifct, _nq, _pQ ): 0;
#endif
    for( unsigned iq=0; iq<_nq; iq++ ){
      if( !_pADJTC[iq].cst() ) return false; // quadrature appears nonlinearly in function
      if( neg ) vHAM[ifct] -= _pQUAD[iq] * _pADJTC[iq];
      else      vHAM[ifct] += _pQUAD[iq] * _pADJTC[iq];
    }
  }

  delete[] _pADJRHS;
  delete[] _pADJQUAD;
#ifndef MC__ODEBNDS_GSL_USE_BAD
  _pADJRHS  = _pDAG->FAD( _nf, vHAM.data(), _nx,   _pVAR+_nx   );
  _pADJQUAD = _pDAG->FAD( _nf, vHAM.data(), _npar, _pVAR+2*_nx );
#else
  _pADJRHS  = _pDAG->BAD( _nf, vHAM.data(), _nx,   _pVAR+_nx   );
  _pADJQUAD = _pDAG->BAD( _nf, vHAM.data(), _npar, _pVAR+2*_nx );
#endif

  delete[] _opADJRHS;  _opADJRHS = 0;
  delete[] _opADJQUAD; _opADJQUAD = 0;
  delete[] _IADJRHS;   _IADJRHS = 0;
  delete[] _PMADJRHS;  _PMADJRHS = 0;
  unsigned Iopmax = 0, PMopmax = 0;

  switch( options.WRAPMIT){
  case OPT::NONE:
    _opADJRHS  = new std::list<const FFOp*>[_nf];
    _opADJQUAD = new std::list<const FFOp*>[_nf];
    for( unsigned ifct=0; ifct<_nf; ifct++ ){
      _opADJRHS[ifct]  = _pDAG->subgraph( _nx,   _pADJRHS+ifct*_nx );
      _opADJQUAD[ifct] = _pDAG->subgraph( _npar, _pADJQUAD+ifct*_npar );
      if( Iopmax < _opADJRHS[ifct].size()  ) Iopmax = _opADJRHS[ifct].size();
      if( Iopmax < _opADJQUAD[ifct].size() ) Iopmax = _opADJQUAD[ifct].size();
    }
    break;

  case OPT::DINEQ:
    _opADJRHS = new std::list<const FFOp*>[_nf*_nx];
    for( unsigned ifx=0; ifx<_nf*_nx; ifx++ ){
      _opADJRHS[ifx] = _pDAG->subgraph( 1, _pADJRHS+ifx );
      if( Iopmax < _opADJRHS[ifx].size()  ) Iopmax = _opADJRHS[ifx].size();
    }
    _opADJQUAD = new std::list<const FFOp*>[_nf];
    for( unsigned ifct=0; ifct<_nf; ifct++ ){
      _opADJQUAD[ifct] = _pDAG->subgraph( _npar, _pADJQUAD+ifct*_npar );
      if( Iopmax < _opADJQUAD[ifct].size() ) Iopmax = _opADJQUAD[ifct].size();
    }
    break;

  case OPT::ELLIPS:
  default:
    _opADJRHS  = new std::list<const FFOp*>[_nf];
    _opADJQUAD = new std::list<const FFOp*>[_nf];
    for( unsigned ifct=0; ifct<_nf; ifct++ ){
      _opADJRHS[ifct]  = _pDAG->subgraph( _nx,   _pADJRHS+ifct*_nx );
      _opADJQUAD[ifct] = _pDAG->subgraph( _npar, _pADJQUAD+ifct*_npar );
      if( PMopmax < _opADJRHS[ifct].size() ) PMopmax = _opADJRHS[ifct].size();
      if( Iopmax < _opADJQUAD[ifct].size() ) Iopmax  = _opADJQUAD[ifct].size();
    }
    break;
  }
  // Intermediate arrays in DAG evaluation
  _IADJRHS   = Iopmax?  new T[ Iopmax ]: 0;
  _PMADJRHS  = PMopmax? new PVT[ PMopmax ]: 0;

  return true;
}

template <typename T, typename PMT, typename PVT>
template <typename REALTYPE, typename OPT> inline bool
ODEBNDS_BASE<T,PMT,PVT>::_RHS_I_ADJ
( const OPT &options, double t, const REALTYPE* y, REALTYPE* ydot,
  REALTYPE*x, const unsigned ifct, const bool neg )
{
  if( !_pADJRHS ) return false; // **error** ADJRHS not defined

  switch( options.WRAPMIT){
  case OPT::NONE:
    _vec2I( x, _nx, _Iz ); // set current state bounds
    _vec2I( y, _nx, _Iy ); // set current adjoint bounds
    *_It = t; // set current time
    _pDAG->eval( _opADJRHS[ifct], _IADJRHS, _nx, _pADJRHS+ifct*_nx, _Iydot,
                 _nVAR-_npar, _pVAR, _IVAR );
    _I2vec( _nx, _Iydot, ydot );
    break;  
   
  case OPT::DINEQ:
    _vec2I( x, _nx, _Iz );  // set current state bounds
    _vec2I( y, _nx, _Iy );  // set current adjoint bounds
    *_It = t; // set current time
    for( unsigned ix=0; ix<_nx; ix++ ){
      T Iyi = _IVAR[ix];
      for( unsigned up=0; up<2; up++ ){ // separate lower/upper bounding subproblems
        _IVAR[ix] = up? Op<T>::u( Iyi ): Op<T>::l( Iyi );
        _pDAG->eval( _opADJRHS[ifct*_nx+ix], _IADJRHS, 1, _pADJRHS+ifct*_nx+ix,
                     _Iydot+ix, _nVAR-_npar, _pVAR, _IVAR );
        if( up && !neg ) _yUdot[ix] = Op<T>::u( _Iydot[ix] );
        else if( up )    _yUdot[ix] = Op<T>::l( _Iydot[ix] );
        else if( !neg )  _yLdot[ix] = Op<T>::l( _Iydot[ix] );
        else             _yLdot[ix] = Op<T>::u( _Iydot[ix] );
      }
      _IVAR[ix] = Iyi;
    }
    _I2vec( _nx, _yLdot, _yUdot, ydot );
    break;  

  case OPT::ELLIPS:
  default:
    _vec2E( x, _nx, _npar, _Q, _Er, _Ir, _pref, _Ip, _B, _zref, _Iz ); // set current state bounds
    _vec2E( y, _nx, _nx, _npar, _Qy, _Edy, _Idy, 0, _Ir, _By, _pref, _Ip,
            _By+_nx*_nx, _yref, _Iy );  // set current adjoint bounds

    // Setup polynomial model expansion of adjoint RHS
    *_MVYZt = t; // Current time
    for( unsigned jx=0; jx<_nx; jx++ )
      _MVYZr[jx].set( _MVYZenv, _nx+jx, _Ir[jx] );
    for( unsigned jy=0; jy<_nx; jy++ )
      _MVYZd[jy].set( _MVYZenv, jy, _Idy[jy] );
    _ep2x( _nx, _npar, _MVYZr, _pref, _MVYZp, _B, _zref, _MVYZz );
    _ep2x( _nx, _nx, _npar, _MVYZd, 0, _MVYZr, _By, _pref, _MVYZp, _By+_nx*_nx,
           _yref, _MVYZy );

    // Construct the adjoint ellipsoidal remainder derivatives
    _pDAG->eval( _opADJRHS[ifct], _PMADJRHS, _nx, _pADJRHS+ifct*_nx, _MVYZg,
                 _nVAR-_npar, _pVAR, _MVYZVAR );
    if( options.QSCALE )
      _RHS_I_ELL( _nx, _MVYZg, _Qy, _Ay, _nz, _yrefdot, _Bydot, _Idydot,
                  _Qydot, options.QTOL, machprec(), _Iy );
    else
      _RHS_I_ELL( _nx, _MVYZg, _Qy, _Ay, _nz, _yrefdot, _Bydot, _Idydot,
                  _Qydot, options.QTOL, machprec() );
    _E2vec( _nx, _nz, _yrefdot, _Qydot, _Bydot, ydot );
    break;
  }
  return true;
}

template <typename T, typename PMT, typename PVT>
template <typename REALTYPE, typename OPT> inline bool
ODEBNDS_BASE<T,PMT,PVT>::_RHS_I_QUAD
( const OPT&options, double t, const REALTYPE*y, REALTYPE*qdot,
  REALTYPE*x, const unsigned ifct, const bool bndinit )
{
  if( !_pADJQUAD ) return false;

  switch( options.WRAPMIT){
  case OPT::NONE:
  case OPT::DINEQ:
    if( !bndinit ) break;
    _vec2I( x, _nx, _Iz );
    _vec2I( y, _nx, _Iy );// set current adjoint bounds
    break;

  case OPT::ELLIPS:
  default:
    if( !bndinit ) break;
    _vec2E( x, _nx, _npar, _Q, _Er, _Ir, _pref, _Ip, _B, _zref, _Iz );
    _vec2E( y, _nx, _nx, _npar, _Qy, _Edy, _Idy, 0, _Ir, _By, _pref, _Ip, 
            _By+_nx*_nx, _yref, _Iy ); // set current adj bounds
    break;
  }

  *_It = t; // set current time
  _QUAD_I( _pDAG, _opADJQUAD[ifct], _IADJRHS, _npar, _pADJQUAD+ifct*_npar,
           _nVAR-_npar, _pVAR, _IVAR, _Iyqdot );
  _I2vec( _npar, _Iyqdot, qdot );
  return true;
}

template <typename T, typename PMT, typename PVT>
template <typename OPT> inline bool
ODEBNDS_BASE<T,PMT,PVT>::_INI_PM_ADJ
( const OPT &options, const unsigned np, const PVT* PMp )
{
  // Update effective number of parameters
  // (possibly larger than _np if lifting is used)
  _npar = np;

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
  _Iy = _IVAR;
  _Iz = _Iy + _nx;
  _Ip = _Iz + _nx;
  _It = _Ip + _npar;
  _Iyq = 0; // not used
  for( unsigned ip=0; ip<_npar; ip++ ) _Ip[ip] = PMp[ip].bound();
  _PMy = _PMVAR;
  _PMz = _PMy + _nx;
  _PMp = _PMz + _nx;
  _PMt = _PMp + _npar;
  _PMyq = _PMt + 1;
  for( unsigned ip=0; ip<_npar; ip++ ) _PMp[ip] = PMp[ip];
  
  // Reset _MVYZenv and related variables
  unsigned MVYZsize = ( options.ORDMIT<_PMenv->nord()? options.ORDMIT: _PMenv->nord() ); 
#ifdef MC__ODEBNDS_BASE_MVYZ_USE
  unsigned MVYZdim  = _nx+_npar;
#else
  unsigned MVYZdim  = ( options.ORDMIT<_PMenv->nord()? _npar: _nx+_npar ); 
#endif
  if( _MVYZenv && ( _MVYZenv->nord() != MVYZsize
                 || _MVYZenv->nvar() != MVYZdim  ) ){
    delete[] _MVYZg;    _MVYZg = 0;
    delete[] _MVYZdgdy; _MVYZdgdy = 0;
    delete[] _MVYZd;    _MVYZd = 0;
    delete[] _MVYZVAR;  _MVYZVAR = _MVYZy = _MVYZz = _MVYZp = _MVYZt = 0; 
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

    if( !_MVYZenv ) _MVYZenv    = new PMT( MVYZdim, MVYZsize );
    _MVYZenv->options = _PMenv->options;
    if( !_MVYZd )    _MVYZd     = new PVT[_nx];
    if( !_MVYZg )    _MVYZg     = new PVT[_nx];
    if( !_MVYZdgdy ) _MVYZdgdy  = new PVT[_nx*_nx];
    if( !_MVYZVAR )  _MVYZVAR   = new PVT[_nVAR];
    _MVYZy = _MVYZVAR;
    _MVYZz = _MVYZy + _nx;
    _MVYZp = _MVYZz + _nx;
    _MVYZt = _MVYZp + _npar;
    for( unsigned ip=0; ip<_npar; ip++ )
      _MVYZp[ip].set( _MVYZenv, ip, _PMp[ip].B() );
    break;
  }

  return true;
}

template <typename T, typename PMT, typename PVT>
template <typename REALTYPE, typename OPT> inline bool
ODEBNDS_BASE<T,PMT,PVT>::_TC_PM_ADJ
( const OPT &options, const double t, REALTYPE*y, unsigned pos_fct, unsigned ifct )
{
  const FFVar* pFCT = _vFCT.at(pos_fct)+ifct;
#ifndef MC__ODEBNDS_GSL_USE_BAD
  delete[] _pADJTC; _pADJTC = _pDAG->FAD( 1, pFCT, _nz, _pVAR+_nx );
#else
  delete[] _pADJTC; _pADJTC = _pDAG->BAD( 1, pFCT, _nz, _pVAR+_nx );
#endif

  _opADJTC  = _pDAG->subgraph( _nx, _pADJTC );
  *_PMt = t; // current time
  _pDAG->eval( _opADJTC, _nx, _pADJTC, _PMy, _nz+1, _pVAR+_nx, _PMVAR+_nx );

  switch( options.WRAPMIT){

  case OPT::NONE:
    // Whether or not to ignore the remainder
    if( !options.PMNOREM )
      _PMI2vec( _PMenv, _nx, _PMy, y, true );
    else
      _PMI2vec( _PMenv, _nx, _PMy, 0, y );
    break;

  case OPT::DINEQ:
    // Whether or not to ignore the remainder
    if( !options.PMNOREM )
      _PMI2vec( _PMenv, _nx, _PMy, y );
    else
      _PMI2vec( _PMenv, _nx, _PMy, 0, y );
    break;

  case OPT::ELLIPS:
  default:
    _IC_PM_ELL( _nx, _PMy, _Qy, _Edy, _Idy );
    // Whether or not to ignore the remainder
    if( !options.PMNOREM )
      _PME2vec( _PMenv, _nx, _PMy, _Qy, y );
    else
      _PME2vec( _PMenv, _nx, _PMy, 0, y );
    break;
  }

  return true;
}

template <typename T, typename PMT, typename PVT>
template <typename REALTYPE, typename OPT> inline bool
ODEBNDS_BASE<T,PMT,PVT>::_TC_PM_QUAD
( const OPT &options, const double t, REALTYPE*yq )
{
  *_PMt = t; // current time
  _opTCQUAD = _pDAG->subgraph( _npar, _pADJTC+_nx );
  _pDAG->eval( _opTCQUAD, _npar, _pADJTC+_nx, _PMyq, _nz+1, _pVAR+_nx, _PMVAR+_nx );

  // Whether or not to ignore the adjoint remainder
  if( !options.PMNOREM )
    _PMI2vec( _PMenv, _npar, _PMyq, yq, true ); // centered
  else
    _PMI2vec( _PMenv, _npar, _PMyq, 0, yq );

  return true;
}

template <typename T, typename PMT, typename PVT>
template <typename REALTYPE, typename OPT> inline bool
ODEBNDS_BASE<T,PMT,PVT>::_CC_PM_ADJ
( const OPT &options, const double t, REALTYPE*y, unsigned pos_fct, unsigned ifct )
{
  const FFVar* pFCT = _vFCT.at(pos_fct-1)+ifct;
#ifndef MC__ODEBNDS_GSL_USE_BAD
  delete[] _pADJTC; _pADJTC = _pDAG->FAD( 1, pFCT, _nz, _pVAR+_nx );
#else
  delete[] _pADJTC; _pADJTC = _pDAG->BAD( 1, pFCT, _nz, _pVAR+_nx );
#endif
  for( unsigned iy=0; iy<_nx; iy++ )
    _pADJCC[iy] = _pY[iy] + _pADJTC[iy];

  _opADJTC = _pDAG->subgraph( _nx, _pADJCC );
#ifdef MC__ODEBNDS_BASE_DINEQPM_DEBUG
  std::cout << "number of operations in adjoint transition: " << _opADJTC.size() << std::endl;
  { int dum; std::cin >> dum; }
#endif

  delete[] _PMADJTC;  _PMADJTC  = 0;
  _opADJDTC.clear();
  delete[] _PMADJDTC; _PMADJDTC = 0;
  delete[] _IADJDTC;  _IADJDTC  = 0;

  switch( options.WRAPMIT){
  case OPT::NONE:
    *_PMt = t; // current time
    _pDAG->eval( _nx, _pADJCC, _PMydot, _nVAR-_npar, _pVAR, _PMVAR );
    // Whether or not to ignore the adjoint remainder
    if( !options.PMNOREM )
      _PMI2vec( _PMenv, _npar, _PMydot, y, true );
    else
      _PMI2vec( _PMenv, _npar, _PMydot, 0, y );
    break;

  case OPT::DINEQ:
    *_PMt = t; // current time
    _pDAG->eval( _nx, _pADJCC, _PMydot, _nVAR-_npar, _pVAR, _PMVAR );
    // Whether or not to ignore the adjoint remainder
    if( !options.PMNOREM )
      _PMI2vec( _PMenv, _npar, _PMydot, y );
    else
      _PMI2vec( _PMenv, _npar, _PMydot, 0, y );
    break;

  case OPT::ELLIPS:
  default:
    *_PMt = t; // current time
#ifdef MC__ODEBNDS_BASE_DINEQPM_DEBUG
    std::cout << "Function:" << ifct << std::endl;
    _print_interm( t, _nx+_npar, _PMz, _Er, "PMz Intermediate", std::cerr );
    _print_interm( t, _nx, _PMy, _Edy, "PMy Intermediate", std::cerr );
#endif

    // In this variant a bound on the Jacobian matrix is computed and the
    // linear part is taken as the mid-point of this matrix
    if( !options.ORDMIT ){
      *_It = t; // set current time
      for( unsigned ix=0; ix<_nx; ix++ ){
        _Iz[ix] = _PMz[ix].bound(); // set current state bounds
        _Iy[ix] = _PMy[ix].bound(); // set current adjoint bounds
        ( _PMy[ix].center() ).set( T(0.) ); // cancel remainder term
      }
      _PMADJTC = new PVT[_opADJTC.size()];
      _pADJDTC = _pDAG->FAD( _nx, _pADJCC, _nx, _pVAR ); // This is the identity matrix!!
      _opADJDTC = _pDAG->subgraph( _nx*_nx, _pADJDTC );
      _IADJDTC = new T[ _opADJDTC.size() ];
      _pDAG->eval( _opADJTC, _PMADJTC, _nx, _pADJCC, _PMydot, _nVAR-_npar,
                   _pVAR, _PMVAR );
      _pDAG->eval( _opADJDTC, _IADJDTC, _nx*_nx, _pADJDTC, _Idgdy, _nVAR-_npar,
                   _pVAR, _IVAR );
      _RHS_PM_ELL0( _nx, _PMydot, _Idgdy, _Idy, _Ay, _Idydot );
    }

    // In this variant a polynomial model of the Jacobian matrix is computed and the
    // linear part is taken as the mid-point of this matrix
    else if( _PMenv->nord() > _MVYZenv->nord() ){
      // Setup polynomial model expansion of adjoint RHS
      *_MVYZt = t; // time
      for( unsigned jx=0; jx<_nx; jx++ ) // state bounds
        _MVYZz[jx].set( _MVYZenv ).set( _PMz[jx].center(), true );
#ifdef MC__ODEBNDS_BASE_MVYZ_USE
      for( unsigned jx=0; jx<_nx; jx++ ){ // adjoint bounds
        _MVYZd[jx].set( _MVYZenv, _npar+jx, _Idy[jx] );
        _MVYZy[jx].set( _MVYZenv ).set( _PMy[jx].center().set( T(0.) ), true );
      }
      _e2x( _nx, _MVYZd, _MVYZy, false );
#else
      for( unsigned jx=0; jx<_nx; jx++ ){
        _MVYZy[jx].set( _MVYZenv ).set( _PMy[jx].center(), true );
        _PMy[jx].set( T(0.) );
      }
#endif
      _PMADJTC = new PVT[_opADJTC.size()];
      _pADJDTC = _pDAG->FAD( _nx, _pADJCC, _nx, _pVAR ); // This is the identity matrix!!
      _opADJDTC = _pDAG->subgraph( _nx*_nx, _pADJDTC );
      _IADJDTC = new T[ _opADJDTC.size() ];
      _pDAG->eval( _opADJTC, _PMADJTC, _nx, _pADJCC, _PMydot, _nVAR-_npar,
                   _pVAR, _PMVAR );
      _pDAG->eval( _opADJDTC, _PMADJDTC, _nx*_nx, _pADJDTC, _MVYZdgdy, _nVAR-_npar,
                   _pVAR, _MVYZVAR );
      _RHS_PM_ELL1( _nx, _PMydot, _MVYZdgdy, _Idy, _Ay, _Idydot );
    }

    // In this variant a polynomial model in the joint adjoint-parameter and
    // of the same order as the parameter polynomial model is computed
    else{
      //{std::cout << "here!\n";  int dum; std::cin >> dum; }
      // Setup polynomial model expansion of adjoint RHS
      *_MVYZt = t; // set current time 
      for( unsigned jx=0; jx<_nx; jx++ ) // state bounds
        _MVYZz[jx].set( _MVYZenv ).set( _PMz[jx].center(), true );
      for( unsigned jx=0; jx<_nx; jx++ ){ // adjoint bounds
        _MVYZd[jx].set( _MVYZenv, _npar+jx, _Idy[jx] );
        _MVYZy[jx].set( _MVYZenv ).set( _PMy[jx].center().set( T(0.) ), true );
      }
      _e2x( _nx, _MVYZd, _MVYZy, false );
      _PMADJTC = new PVT[_opADJTC.size()];
      _pDAG->eval( _opADJTC, _PMADJTC, _nx, _pADJCC, _MVYZg, _nVAR-_npar,
                   _pVAR, _MVYZVAR );
      _RHS_PM_ELL2( _nx, _PMenv, _PMydot, _MVYZg, _npar, _Idy, _Ay, _Idydot );
    }

    // Whether or not to ignore the remainder
    if( !options.PMNOREM ){
      _CC_PM_ELL( _nx, _Edy, _Ay, _Idydot, _Qydot, options.QTOL, machprec() );
      _PME2vec( _PMenv, _nx, _PMydot, _Qydot, y );
    }
    else
      _PME2vec( _PMenv, _nx, _PMydot, 0, y );
    break;
  }

  return true;
}

template <typename T, typename PMT, typename PVT> 
template <typename REALTYPE, typename OPT> inline bool
ODEBNDS_BASE<T,PMT,PVT>::_CC_PM_QUAD
( const OPT&options, const double t, REALTYPE*yq )
{
  *_PMt = t; // current time
  _pDAG->eval( _npar, _pADJTC+_nx, _PMyq, _nVAR-_npar, _pVAR, _PMVAR, true );

  // Whether or not to ignore the adjoint remainder
  if( !options.PMNOREM )
    _PMI2vec( _PMenv, _npar, _PMyq, yq, true );
  else
    _PMI2vec( _PMenv, _npar, _PMyq, 0, yq );

  return true;
}

template <typename T, typename PMT, typename PVT> inline bool
ODEBNDS_BASE<T,PMT,PVT>::_IC_PM_ADJ
( const double t )
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
  *_PMt = t; // current time
  _pDAG->eval( _npar, _pADJTC, _PMyq, _nVAR-_npar, _pVAR, _PMVAR, true );

  return true;
}

template <typename T, typename PMT, typename PVT>
template <typename OPT> inline bool
ODEBNDS_BASE<T,PMT,PVT>::_SET_PM_ADJ
( const OPT &options, unsigned iADJRHS, unsigned iQUAD, unsigned pos_fct,
  const bool neg )
{
  if( _vRHS.size() <= iADJRHS ) return false; 
  if( _nq && _vQUAD.size() <= iQUAD ) return false;

  _pRHS = _vRHS.at( iADJRHS );
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
    delete[] _pADJTC; _pADJTC = _nq? _pDAG->FAD( 1, pFCT+ifct, _nq, _pQ ): 0;
#else
    delete[] _pADJTC; _pADJTC = _nq? _pDAG->BAD( 1, pFCT+ifct, _nq, _pQ ): 0;
#endif
    for( unsigned iq=0; iq<_nq; iq++ ){
      if( !_pADJTC[iq].cst() ) return false; // quadrature appears nonlinearly in function
      if( neg ) vHAM[ifct] -= _pQUAD[iq] * _pADJTC[iq];
      else      vHAM[ifct] += _pQUAD[iq] * _pADJTC[iq];
    }
  }

  delete[] _pADJRHS;
  delete[] _pADJQUAD;
#ifndef MC__ODEBNDS_GSL_USE_BAD
  _pADJRHS  = _pDAG->FAD( _nf, vHAM.data(), _nx,   _pVAR+_nx   );
  _pADJQUAD = _pDAG->FAD( _nf, vHAM.data(), _npar, _pVAR+2*_nx );
#else
  _pADJRHS  = _pDAG->BAD( _nf, vHAM.data(), _nx,   _pVAR+_nx   );
  _pADJQUAD = _pDAG->BAD( _nf, vHAM.data(), _npar, _pVAR+2*_nx );
#endif
  delete[] _pADJJAC;   _pADJJAC = 0;
  delete[] _opADJJAC;  _opADJJAC = 0;
  delete[] _IADJJAC;   _IADJJAC = 0;
  delete[] _PMADJJAC;  _PMADJJAC = 0;

  delete[] _opADJRHS;  _opADJRHS = 0;
  delete[] _opADJQUAD; _opADJQUAD = 0;
  delete[] _IADJRHS;   _IADJRHS = 0;
  delete[] _PMADJRHS;  _PMADJRHS = 0;
  unsigned Iopmax = 0, PMopmax = 0, JACopmax = 0;

  switch( options.WRAPMIT){
  case OPT::NONE:
    _opADJRHS  = new std::list<const FFOp*>[_nf];
    _opADJQUAD = new std::list<const FFOp*>[_nf];
    for( unsigned ifct=0; ifct<_nf; ifct++ ){
      _opADJRHS[ifct]  = _pDAG->subgraph( _nx,   _pADJRHS+ifct*_nx );
      _opADJQUAD[ifct] = _pDAG->subgraph( _npar, _pADJQUAD+ifct*_npar );
      if( PMopmax < _opADJRHS[ifct].size()  ) PMopmax = _opADJRHS[ifct].size();
      if( PMopmax < _opADJQUAD[ifct].size() ) PMopmax = _opADJQUAD[ifct].size();
    }
    break;

  case OPT::DINEQ:
    _opADJRHS = new std::list<const FFOp*>[_nf*_nx];
    for( unsigned ifx=0; ifx<_nf*_nx; ifx++ ){
      _opADJRHS[ifx] = _pDAG->subgraph( 1, _pADJRHS+ifx );
      if( PMopmax < _opADJRHS[ifx].size()  ) PMopmax = _opADJRHS[ifx].size();
    }
    _opADJQUAD = new std::list<const FFOp*>[_nf];
    for( unsigned ifct=0; ifct<_nf; ifct++ ){
      _opADJQUAD[ifct] = _pDAG->subgraph( _npar, _pADJQUAD+ifct*_npar );
      if( PMopmax < _opADJQUAD[ifct].size() ) PMopmax = _opADJQUAD[ifct].size();
    }
    break;

  case OPT::ELLIPS:
  default:
    _opADJRHS  = new std::list<const FFOp*>[_nf];
    _opADJQUAD = new std::list<const FFOp*>[_nf];
    for( unsigned ifct=0; ifct<_nf; ifct++ ){
      _opADJRHS[ifct]  = _pDAG->subgraph( _nx,   _pADJRHS+ifct*_nx );
      _opADJQUAD[ifct] = _pDAG->subgraph( _npar, _pADJQUAD+ifct*_npar );
      if( PMopmax < _opADJRHS[ifct].size() ) PMopmax = _opADJRHS[ifct].size();
      if( Iopmax < _opADJQUAD[ifct].size() ) Iopmax  = _opADJQUAD[ifct].size();
    }
    if( options.ORDMIT && _PMenv->nord() <= _MVYZenv->nord() ) break;
    _opADJJAC  = new std::list<const FFOp*>[_nf];
    _pADJJAC = _pDAG->FAD( _nx*_nf, _pADJRHS, _nx, _pVAR );
    for( unsigned ifct=0; ifct<_nf; ifct++ ){
      //_pADJJAC = _pDAG->FAD( _nx, _pADJRHS[ifct], _nx, _pVAR );
      _opADJJAC[ifct] = _pDAG->subgraph( _nx*_nx, _pADJJAC+ifct*_nx*_nx );
      if( JACopmax < _opADJJAC[ifct].size() ) JACopmax = _opADJJAC[ifct].size();
    }
    if( !options.ORDMIT )
      _IADJJAC  = JACopmax?  new T[ JACopmax ]: 0;
    else if( _PMenv->nord() > _MVYZenv->nord() )
      _PMADJJAC = JACopmax?  new PVT[ JACopmax ]: 0;
  }

  // Intermediate arrays in DAG evaluation
  _IADJRHS   = Iopmax?  new T[ Iopmax ]: 0;
  _PMADJRHS  = PMopmax? new PVT[ PMopmax ]: 0;

  return true;
}

template <typename T, typename PMT, typename PVT>
template <typename REALTYPE, typename OPT> inline bool
ODEBNDS_BASE<T,PMT,PVT>::_RHS_PM_ADJ
( const OPT&options, double t, const REALTYPE*y, REALTYPE*ydot,
  REALTYPE*x, const unsigned ifct, const bool neg )
{
  if( !_pADJRHS ) return false;

  switch( options.WRAPMIT){
  case OPT::NONE:
    *_PMt = t; // current time
    _vec2PMI( x, _PMenv, _nx, _PMz, true );// set current state bounds
    _vec2PMI( y, _PMenv, _nx, _PMy, true );// set current adjoint bounds
    _pDAG->eval( _opADJRHS[ifct], _PMADJRHS, _nx, _pADJRHS+ifct*_nx, _PMydot,
                 _nVAR-_npar, _pVAR, _PMVAR );
    if( !options.PMNOREM )
      _PMI2vec( _PMenv, _nx, _PMydot, ydot, true );
    else
      _PMI2vec( _PMenv, _nx, _PMydot, 0, ydot );
    break;
  
  case OPT::DINEQ:
    *_PMt = t; // current time
    _vec2PMI( x, _PMenv, _nx, _PMz );// set current state bounds
    _vec2PMI( y, _PMenv, _nx, _PMy );// set current adjoint bounds
    for( unsigned ix=0; ix<_nx; ix++ ){
      T Ryi = _PMVAR[ix].remainder();
      for( unsigned up=0; up<2; up++ ){ // separate lower/upper bounding subproblems
        if( up ) _PMVAR[ix].set( Op<T>::u( Ryi ) );
        else     _PMVAR[ix].set( Op<T>::l( Ryi ) );
        _pDAG->eval( _opADJRHS[ifct*_nx+ix], _PMADJRHS, 1, _pADJRHS+ifct*_nx+ix,
                     _PMydot+ix, _nVAR-_npar, _pVAR, _PMVAR );
        if( up && !neg ) _RyUdot[ix] = Op<T>::u( _PMydot[ix].remainder() );
        else if( up )    _RyUdot[ix] = Op<T>::l( _PMydot[ix].remainder() );
        else if( !neg )  _RyLdot[ix] = Op<T>::l( _PMydot[ix].remainder() );
        else             _RyLdot[ix] = Op<T>::u( _PMydot[ix].remainder() );
      }
      _PMVAR[ix].set( Ryi );
    }
    //_RHS_PM_DI( _pDAG, _opADJRHS+ifct*_nx, _PMADJRHS, _nx, _pADJRHS+ifct*_nx,
    //            _nVAR-_npar, _pVAR, _PMVAR, _PMydot, _RyLdot, _RyUdot, neg );
    if( !options.PMNOREM )
      _PMI2vec( _PMenv, _nx, _PMydot, _RyLdot, _RyUdot, ydot );
    else
      _PMI2vec( _PMenv, _nx, _PMydot, 0, ydot );
#ifdef MC__ODEBND_BASE_DINEQPM_DEBUG
    _print_interm( t, _nx, _PMydot, "PMydot Intermediate", std::cerr );
    _print_interm( t, _nx, _RyLdot, "RyLdot Intermediate", std::cerr );
    _print_interm( t, _nx, _RyUdot, "RyUdot Intermediate", std::cerr );
    { std::cout << "--paused--"; int dum; std::cin >> dum; }
#endif
    break;

  case OPT::ELLIPS:
  default:
    *_PMt = t; // current time
    _vec2PME( x, _PMenv, _nx, _PMz, _Q, _Er, _Ir );// set current state bounds
#ifdef MC__ODEBND_BASE_DINEQPM_DEBUG
    std::cout << "Function:" << ifct << std::endl;
    _print_interm( t, _nx+_npar, _PMz, _Er, "PMz Intermediate", std::cerr );
#endif
    _vec2PME( y, _PMenv, _nx, _PMy, _Qy, _Edy, _Idy ); // set current adjoint bounds
#ifdef MC__ODEBND_BASE_DINEQPM_DEBUG
    _print_interm( t, _nx, _PMy, _Edy, "PMy Intermediate", std::cerr );
#endif

    // In this variant a bound on the Jacobian matrix is computed and the
    // linear part is taken as the mid-point of this matrix
    if( !options.ORDMIT ){
      *_It = t; // set current time
      for( unsigned ix=0; ix<_nx; ix++ ){
        _Iz[ix] = _PMz[ix].bound(); // set current state bounds
        _Iy[ix] = _PMy[ix].bound(); // set current adjoint bounds
        ( _PMy[ix].center() ).set( T(0.) ); // cancel remainder term
      }
      _pDAG->eval( _opADJRHS[ifct], _PMADJRHS, _nx, _pADJRHS+ifct*_nx, _PMydot,
                   _nVAR-_npar, _pVAR, _PMVAR );
      _pDAG->eval( _opADJJAC[ifct], _IADJJAC, _nx*_nx, _pADJJAC+ifct*_nx*_nx,
                   _Idgdy, _nVAR-_npar, _pVAR, _IVAR );
      _RHS_PM_ELL0( _nx, _PMydot, _Idgdy, _Idy, _Ay, _Idydot );
    }

    // In this variant a polynomial model of the Jacobian matrix is computed and the
    // linear part is taken as the mid-point of this matrix
    else if( _PMenv->nord() > _MVYZenv->nord() ){
      // Setup polynomial model expansion of adjoint RHS
      *_MVYZt = t; // time
      for( unsigned jx=0; jx<_nx; jx++ ) // state bounds
        _MVYZz[jx].set( _MVYZenv ).set( _PMz[jx].center(), true );
#ifdef MC__ODEBNDS_BASE_MVYZ_USE
      for( unsigned jx=0; jx<_nx; jx++ ){ // adjoint bounds
        _MVYZd[jx].set( _MVYZenv, _npar+jx, _Idy[jx] );
        _MVYZy[jx].set( _MVYZenv ).set( _PMy[jx].center().set( T(0.) ), true );
      }
      _e2x( _nx, _MVYZd, _MVYZy, false );
#else
      for( unsigned jx=0; jx<_nx; jx++ ){
        _MVYZy[jx].set( _MVYZenv ).set( _PMy[jx].center(), true );
        _PMy[jx].set( T(0.) );
      }
#endif
      _pDAG->eval( _opADJRHS[ifct], _PMADJRHS, _nx, _pADJRHS+ifct*_nx, _PMydot,
                   _nVAR-_npar, _pVAR, _PMVAR );
      _pDAG->eval( _opADJJAC[ifct], _PMADJJAC, _nx*_nx, _pADJJAC+ifct*_nx*_nx,
                   _MVYZdgdy, _nVAR-_npar, _pVAR, _MVYZVAR );
      _RHS_PM_ELL1( _nx, _PMydot, _MVYZdgdy, _Idy, _Ay, _Idydot );
    }

    // In this variant a polynomial model in the joint adjoint-parameter and
    // of the same order as the parameter polynomial model is computed
    else{
      //{std::cout << "here!\n";  int dum; std::cin >> dum; }
      // Setup polynomial model expansion of adjoint RHS
      *_MVYZt = t; // set current time 
      for( unsigned jx=0; jx<_nx; jx++ ) // state bounds
        _MVYZz[jx].set( _MVYZenv ).set( _PMz[jx].center(), true );
      for( unsigned jx=0; jx<_nx; jx++ ){ // adjoint bounds
        _MVYZd[jx].set( _MVYZenv, _npar+jx, _Idy[jx] );
        _MVYZy[jx].set( _MVYZenv ).set( _PMy[jx].center().set( T(0.) ), true );
      }
      _e2x( _nx, _MVYZd, _MVYZy, false );
      _pDAG->eval( _opADJRHS[ifct], _PMADJRHS, _nx, _pADJRHS+ifct*_nx, _MVYZg,
                   _nVAR-_npar, _pVAR, _MVYZVAR );
      _RHS_PM_ELL2( _nx, _PMenv, _PMydot, _MVYZg, _npar, _Idy, _Ay, _Idydot );
    }

    // Construct the ellipsoidal remainder derivatives
    if( options.QSCALE )
      _RHS_PM_ELL( _nx, _Qy, _Ay, _Idydot, _Qydot, options.QTOL, machprec(), _PMy );
    else
      _RHS_PM_ELL( _nx, _Qy, _Ay, _Idydot, _Qydot, options.QTOL, machprec() );

    // Whether or not to ignore the adjoint remainder
    if( !options.PMNOREM )
      _PME2vec( _PMenv, _nx, _PMydot, _Qydot, ydot );
    else
      _PME2vec( _PMenv, _nx, _PMydot, 0, ydot );
#ifdef MC__ODEBND_BASE_DINEQPM_DEBUG
    _print_interm( t, _nx, _PMydot, E(_nx,_Qydot), "PMydot Intermediate", std::cerr );
    { std::cout << "--paused--"; int dum; std::cin >> dum; }
#endif
    break;
  }

  return true;
}

template <typename T, typename PMT, typename PVT>
template <typename REALTYPE, typename OPT> inline bool
ODEBNDS_BASE<T,PMT,PVT>::_RHS_PM_QUAD
( const OPT&options, double t, const REALTYPE*y, REALTYPE*yqdot, 
  REALTYPE*x, const unsigned ifct, const bool bndinit )
{
  if( !_pADJQUAD ) return false;

  switch( options.WRAPMIT){
  case OPT::NONE:
    if( !bndinit ) break;
    _vec2PMI( x, _PMenv, _nx, _PMz, true );// set current state bounds
    _vec2PMI( y, _PMenv, _nx, _PMy, true );// set current adjoint bounds
    break;

  case OPT::DINEQ:
    if( !bndinit ) break;
    _vec2PMI( x, _PMenv, _nx, _PMz );// set current state bounds
    _vec2PMI( y, _PMenv, _nx, _PMy );// set current adjoint bounds
    break;

  case OPT::ELLIPS:
  default:
    if( !bndinit ) break;
    _vec2PME( x, _PMenv, _nx, _PMz, _Q, _Er, _Ir );// set current state bounds
    _vec2PME( y, _PMenv, _nx, _PMy, _Qy, _Edy, _Idy );// set current adjoint bounds
    break;
  }

  *_PMt = t; // current time
  _QUAD_PM( _pDAG, _opADJQUAD[ifct], _PMADJRHS, _npar, _pADJQUAD+ifct*_npar,
            _nVAR-_npar, _pVAR, _PMVAR, _PMyqdot );

  // Whether or not to ignore the adjoint remainder
  if( !options.PMNOREM )
    _PMI2vec( _PMenv, _npar, _PMyqdot, yqdot, true );
  else
    _PMI2vec( _PMenv, _npar, _PMyqdot, 0, yqdot );

  return true;  
}

} // end namescape mc

#endif

















