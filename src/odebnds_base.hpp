// Copyright (C) 2015 Benoit Chachuat & Nikola Peric, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

#ifndef MC__ODEBNDS_BASE_HPP
#define MC__ODEBNDS_BASE_HPP

#undef  MC__ODEBNDS_BASE_DINEQI_DEBUG
#undef  MC__ODEBNDS_BASE_DINEQPM_DEBUG
#undef  MC__ODEBNDS_BASE_MVYXP_USE

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
  using ODEBND_BASE<T,PMT,PVT>::_A;
  using ODEBND_BASE<T,PMT,PVT>::_B;
  using ODEBND_BASE<T,PMT,PVT>::_xref;
  using ODEBND_BASE<T,PMT,PVT>::_Ixdot;
  using ODEBND_BASE<T,PMT,PVT>::_Idfdx;
  using ODEBND_BASE<T,PMT,PVT>::_MVXPenv;
  using ODEBND_BASE<T,PMT,PVT>::_MVXPVAR;
  using ODEBND_BASE<T,PMT,PVT>::_MVXPt;
  using ODEBND_BASE<T,PMT,PVT>::_MVXPx;
  using ODEBND_BASE<T,PMT,PVT>::_MVXPp;
  using ODEBND_BASE<T,PMT,PVT>::_MVXPd;
  using ODEBND_BASE<T,PMT,PVT>::_MVXPf;
  using ODEBND_BASE<T,PMT,PVT>::_MVXPdfdx;
  using ODEBND_BASE<T,PMT,PVT>::_diam;
  using ODEBND_BASE<T,PMT,PVT>::_vec2I;
  using ODEBND_BASE<T,PMT,PVT>::_vec2E;
  using ODEBND_BASE<T,PMT,PVT>::_ep2x;
  using ODEBND_BASE<T,PMT,PVT>::_E2vec;
  using ODEBND_BASE<T,PMT,PVT>::_I2vec;
  using ODEBND_BASE<T,PMT,PVT>::_QUAD_I_ELL;
  using ODEBND_BASE<T,PMT,PVT>::_ndxLT;

  using ODEBND_BASE<T,PMT,PVT>::_PMenv;
  using ODEBND_BASE<T,PMT,PVT>::_PMVAR;
  using ODEBND_BASE<T,PMT,PVT>::_PMt;
  using ODEBND_BASE<T,PMT,PVT>::_PMp;
  using ODEBND_BASE<T,PMT,PVT>::_PMx;
  using ODEBND_BASE<T,PMT,PVT>::_vec2PMI;
  using ODEBND_BASE<T,PMT,PVT>::_vec2PME;
  using ODEBND_BASE<T,PMT,PVT>::_PMI2vec;
  using ODEBND_BASE<T,PMT,PVT>::_PME2vec;
  using ODEBND_BASE<T,PMT,PVT>::_QUAD_PM;
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

  //! @brief pointer to adjoint parameter (states) interval bounds **DO NOT FREE**
  T *_Ix;

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

  //! @brief adjoint quadrature lower bound time derivatives
  double *_yqLdot;

  //! @brief adjoint quadrature upper bound time derivatives
  double *_yqUdot;

  //! @brief adjoint polynomial model **DO NOT FREE**
  PVT *_PMy;

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

  //! @brief linear transformation A matrix (adjoint system) (to multiply states)
  double *_Ax;

  //! @brief linear transformation B matrix (adjoint system)
  double *_By;

  //! @brief linear transformation B matrix (adjoint quadrature)
  double *_Byq;

  //! @brief linear transformed adjoint interval bounds
  T *_Idy;

  //! @brief linear transformed adjoint quadrature bounds
  T *_Idyq;

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

  //! @brief linear transformation B matrix time derivatives (adjoint quadrature)
  double *_Byqdot;

  //! @brief rotated adjoint interval bounds time derivatives
  T *_Idydot;

  //! @brief rotated adjoint quadrature bounds time derivatives
  T *_Idyqdot;

  //! @brief shape matrix time directives in ellipsoidal bounds (adjoint system)
  double *_Qydot;

  //! @brief polynomial model environment for mean-value theorem in (Y,X,P)
  PMT *_MVYXPenv;

  //! @brief rotated scaled state polynomial model
  PVT *_MVYXPw;

  //! @brief rotated unscaled state polynomial model
  PVT *_MVYXPr;

  //! @brief rotated adjoint polynomial model
  PVT *_MVYXPd;

  //! @brief adjoint RHS polynomial model
  PVT *_MVYXPg;

  //! @brief adjoint RHS Jacobian polynomial model
  PVT *_MVYXPdgdy;

  //! @brief adjoint quadrature RHS polynomial model
  PVT *_MVYXPyqdot;

  //! @brief polynomial models for variables in mean-value theorem (adjoint system)
  PVT *_MVYXPVAR;

  //! @brief pointer to time polynomial models (adjoint system) **DO NOT FREE**
  PVT *_MVYXPt;

  //! @brief pointer to adjoint polynomial models **DO NOT FREE**
  PVT *_MVYXPy;

  //! @brief pointer to state polynomial models (adjoint system) **DO NOT FREE**
  PVT *_MVYXPx;

  //! @brief pointer to parameter polynomial models (adjoint system) **DO NOT FREE**
  PVT *_MVYXPp;

  //! @brief polynomial model environment for mean-value theorem in Z
  //PMT *_MVXPenv;

  //! @brief pointer to rotated state polynomial models (adjoint initialization)
  //PVT *_MVXPd;

  //! @brief pointer to adjoint polynomial models (adjoint initialization)
  PVT *_MVXPy;

  //! @brief pointer to adjoint quadrature polynomial models (adjoint initialization)
  PVT *_MVXPyq;

  //! @brief polynomial models for variables in mean-value theorem (adjoint initialization)
  //PVT *_MVXPVAR;

  //! @brief pointer to state polynomial models (adjoint initialization) **DO NOT FREE**
  //PVT *_MVXPx; //can go -> MVXPx

  //! @brief pointer to parameter polynomial models (adjoint initialization) **DO NOT FREE**
  //PVT *_MVXPp;

  //! @brief pointer to time polynomial models (adjoint initialization) **DO NOT FREE**
  //PVT *_MVXPt;

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

  //! @brief Function to reinitialize adjoint interval bounds after discontinuity
  template <typename REALTYPE, typename OPT> bool _CC_I_ADJ
    ( const OPT &options, const double t, REALTYPE *vec,
      unsigned pos_fct, unsigned ifct );

  //! @brief Function to reinitialize adjoint interval bounds w/ ellipsoidal bounds 
  static void _CC_I_ELL
    ( const unsigned ny, PVT*MVYXPg, T*Iy, const E&Edy, double*Ay, const E&Edx,
      double*Ax, const unsigned np, double*yref, double*By, T*Idy, double*Qy,
      const double QTOL, const double EPS );

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

  //! @brief Static function to calculate the RHS of adjoint ODEs in interval arithmetic w/ ellipsoidal contractor
  static void _RHS_I_ELL
    ( const unsigned nx, PVT*MVYXPg, const double*Qy, double*Ay, double *Ax,
      const unsigned np, double*yrefdot, double*Bydot, T*Idydot, double*Qydot,
      const double QTOL, const double EPS, const bool neg, const T*W=0 );

  //! @brief Static function to calculate the RHS of adjoint ODEs w/ ellipsoidal contractor
  template <typename U> static void _RHS_ELL
    ( const unsigned nx, const double*Qy, const double*Ay, const double *Ax,
      const T*Idydot, double*Qydot, const double QTOL, const double EPS,
      const bool neg, const U*W );

  //! @brief Function to calculate the RHS of auxiliary ODEs in interval arithmetic
  template <typename REALTYPE, typename OPT> bool _RHS_I_QUAD
    ( const OPT&options, double t, const REALTYPE*y, REALTYPE*qdot,
      REALTYPE* x, const unsigned ifct, const bool neg=false );

  //! @brief Function to initialize adjoint polynomial models
  template <typename OPT> bool _INI_PM_ADJ
    ( const OPT &options, const unsigned np, const PVT *PMp );

  //! @brief Function to initialize adjoint polynomial model
  template <typename REALTYPE, typename OPT> bool _TC_PM_ADJ
    ( const OPT &options, const double t, const REALTYPE *x, REALTYPE *y,
      unsigned pos_fct, unsigned ifct );

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

  //! @brief Static function to calculate the RHS of auxiliary adjoint ODEs in polynomial model arithmetic w/ ellipsoidal contractor - approximation using mean-value theorem and interval analysis
  static void _RHS_PM_ELL0
    ( const unsigned nx, PVT*PMg, T*Idgdy, const T*Ix, double*Ax,
      const T*Iy, double*Ay, T*Idydot );

  //! @brief Static function to calculate the RHS of auxiliary adjoint ODEs in polynomial model arithmetic w/ ellipsoidal contractor - approximation using mean-value theorem and PM arithmetic
  static void _RHS_PM_ELL1
    ( const unsigned nx, PVT*PMg, PVT*MVYXPdgdy, const T*Ix, double*Ax,
      const T*Iy, double*Ay, T*Idydot );

  //! @brief Static function to calculate the RHS of auxiliary adjoint ODEs in polynomial model arithmetic w/ ellipsoidal contractor - joint polynomial model in states and parameters
  static void _RHS_PM_ELL2
    ( const unsigned nx, PMT*PMenv, PVT*PMg, PVT*MVYXPg, const unsigned np,
      const T*Ix, double*Ax, const T*Iy, double*Ay, T*Idydot );

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
  _IADJRHS = _IADJJAC = _Iy = _Ix  = _IADJDTC = 0;
  _PMADJRHS = _PMADJJAC = _PMy = _PMydot = _PMyq = _PMyqdot =
  _PMADJTC = _PMADJDTC = 0;

  // Initialize parameterization arrays
  _zref = _yref = _yrefdot = 0;
  _Ay = _Ax = _By = _Byq = _Qy = _Bydot = _Byqdot = _Qydot = 0;
  _Idy = _Idydot = _Iydot = _Iyqdot = _Idgdy = _Idyq = _Idyqdot = 0;
  _yLdot = _yUdot = _yqLdot = _yqUdot = _RyLdot = _RyUdot = _Ryqdot = 0;

  // Initialize Taylor model environments
  _MVYXPenv = 0;//_MVXPenv = 0;
  _MVYXPd = _MVYXPr = _MVYXPw = _MVYXPg = _MVYXPdgdy = _MVYXPyqdot = _MVYXPVAR =
  _MVYXPy = _MVYXPx = _MVYXPp = _MVYXPt = 0;
  _MVXPy  = _MVXPyq = 0; //_MVXPd  = _MVXPVAR  = _MVXPx  = _MVXPp = _MVXPt = 0;
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
  delete[] _yqLdot;
  delete[] _yqUdot;
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
  delete[] _Ax;
  delete[] _By;
  delete[] _Byq;
  delete[] _Qy;
  delete[] _yrefdot;
  delete[] _Bydot;
  delete[] _Byqdot;
  delete[] _Idy;
  delete[] _Idyq;
  delete[] _Idydot;
  delete[] _Idyqdot;
  delete[] _Qydot;
  delete[] _Idgdy;
  // ** DO NOT DELETE _MVYXPy, _MVYXPx, _MVYXPp, _MVYXPt **
  delete   _MVYXPenv;
  delete[] _MVYXPd;
  delete[] _MVYXPr;
  delete[] _MVYXPw;
  delete[] _MVYXPg;
  delete[] _MVYXPdgdy;
  delete[] _MVYXPVAR;
  delete[] _MVYXPyqdot;
  // ** DO NOT DELETE _MVXPx, _MVXPp, _MVXPt **
  //delete   _MVXPenv;
  //delete[] _MVXPd;
  delete[] _MVXPy;
  delete[] _MVXPyq;
  //delete[] _MVXPVAR;
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
  _nVAR = _nx + _nx + _npar + 1 + _nq + _npar;
  delete[] _pVAR; _pVAR = new FFVar[_nVAR];
  delete[] _IVAR; _IVAR = new T[_nVAR];
  delete[] _pADJCC;  _pADJCC  = new FFVar[_nx];

  BASE_DE::set_adjoint();
  for( unsigned ix=0; ix<_nx; ix++ ) _pVAR[ix] = _pY[ix];
  for( unsigned ix=0; ix<_nx; ix++ ) _pVAR[_nx+ix] = _pX[ix];
  for( unsigned ip=0; ip<_npar; ip++ ) _pVAR[2*_nx+ip] = _pP[ip];
  _pVAR[2*_nx+_npar] = (_pT? *_pT: 0. );
  _Iy = _IVAR;
  _Ix = _Iy + _nx;
  _Ip = _Ix + _nx;
  _It = _Ip + _npar;
  _Iyq = _It + 1;
  for( unsigned ip=0; ip<_npar; ip++ ) _Ip[ip] = Ip[ip];

  // Reset _MVYZVenv and related variables
  if( _MVYXPenv && ( _MVYXPenv->nord() != options.ORDMIT 
                  || _MVYXPenv->nvar() != 2*_nx+_npar ) ){
    delete[] _MVYXPg;   _MVYXPg = 0;
    delete[] _MVYXPd;   _MVYXPd = 0;
    delete[] _MVYXPr;   _MVYXPr = 0;
    delete[] _MVYXPw;   _MVYXPw = 0;
    delete[] _MVYXPVAR; _MVYXPVAR = _MVYXPy = _MVYXPx = _MVYXPp = _MVYXPt = 0;
    delete[] _MVYXPyqdot; _MVYXPyqdot = 0;
    delete   _MVYXPenv; _MVYXPenv = 0;
  }
/*
  // Reset _MVZVenv and related variables
  if( _MVXPenv && ( _MVXPenv->nord() != 1 
                 || _MVXPenv->nvar() != _nx+_npar ) ){
    delete[] _MVXPy;   _MVXPy = 0;
    delete[] _MVXPd;   _MVXPd = 0;
    delete[] _MVXPVAR; _MVXPVAR = 0;_MVYXPyqdot
    delete   _MVXPenv; _MVXPenv = 0;
  }
*/
  // Set parameterization variables
  delete[] _Iyqdot;   _Iyqdot = new T[_npar];
  delete[] _yqLdot;   _yqLdot = new double[_npar];
  delete[] _yqUdot;   _yqUdot = new double[_npar];

  switch( options.WRAPMIT){
  case OPT::NONE:
  case OPT::DINEQ:
    delete[] _yLdot;   _yLdot = new double[_nx];
    delete[] _yUdot;   _yUdot = new double[_nx];
    delete[] _Iydot;   _Iydot = new T[_nx];
    break;

  case OPT::ELLIPS:
  default:
    for( unsigned ip=0; ip<_npar; ip++ )
      _pref[ip] = Op<T>::mid( _Ip[ip] );
    delete[] _yref;    _yref    = new double[_nx];
    delete[] _Ay;      _Ay      = new double[_nx*_nx];
    delete[] _Ax;      _Ax      = new double[_nx*_nx];
    delete[] _By;      _By      = new double[_nx*_npar];
    delete[] _Byq;     _Byq     = new double[_npar*_npar];
    delete[] _Qy;      _Qy      = new double[_nx*(_nx+1)/2];
    delete[] _Idy;     _Idy     = new T[_nx];
    delete[] _Idyq;    _Idyq    = new T[_nx];
    delete[] _yrefdot; _yrefdot = new double[_nx];
    delete[] _Bydot;   _Bydot   = new double[_nx*_npar];
    delete[] _Byqdot;  _Byqdot  = new double[_npar*_npar];
    delete[] _Qydot;   _Qydot   = new double[_nx*(_nx+1)/2];
    delete[] _Idydot;  _Idydot  = new T[_nx];
    delete[] _Idyqdot; _Idyqdot = new T[_npar];

    if( !_MVYXPenv )  _MVYXPenv = new PMT( 2*_nx+_npar, options.ORDMIT );
    _MVYXPenv->options = options.PMOPT;
    if( !_MVYXPd )    _MVYXPd   = new PVT[_nx];
    if( !_MVYXPr )    _MVYXPr   = new PVT[_nx];
    if( !_MVYXPw )    _MVYXPw   = new PVT[_nx];
    if( !_MVYXPg )    _MVYXPg   = new PVT[_nx];
    if( !_MVYXPyqdot ) _MVYXPyqdot = new PVT[_npar];
    if( !_MVYXPVAR )  _MVYXPVAR = new PVT[_nVAR];
    _MVYXPy = _MVYXPVAR;
    _MVYXPx = _MVYXPy + _nx;
    _MVYXPp = _MVYXPx + _nx;
    _MVYXPt = _MVYXPp + _npar;
    for( unsigned jp=0; jp<_npar; jp++ )
      _MVYXPp[jp].set( _MVYXPenv, _nx+_nx+jp, _Ip[jp] );

    //delete   _MVXPenv;  _MVXPenv  = new PMT( _nx+_npar, 1 );
    //_MVXPenv->options = options.PMOPT;
    //delete[] _MVXPd;    _MVXPd    = new PVT[_nx];
    delete[] _MVXPy;    _MVXPy    = new PVT[_nx];
    delete[] _MVXPyq;   _MVXPyq   = new PVT[_npar];
    //delete[] _MVXPVAR;  _MVXPVAR  = new PVT[_nx+_npar+1];
    //_MVXPx = _MVXPVAR;
    //_MVXPp = _MVXPx + _nx;
    //_MVXPt = _MVXPp + _npar;
    //for( unsigned ip=0; ip<_npar; ip++ )
    //  _MVXPp[ip].set( _MVXPenv, _nx+ip, _Ip[ip] );
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
  delete[] _pADJTC; _pADJTC = _pDAG->FAD( 1, pFCT, _nx+_npar, _pVAR+_nx );
#else
  delete[] _pADJTC; _pADJTC = _pDAG->BAD( 1, pFCT, _nx+_npar, _pVAR+_nx );
#endif
  _opADJTC  = _pDAG->subgraph( _nx, _pADJTC );

  switch( options.WRAPMIT){

  case OPT::NONE:
  case OPT::DINEQ:
    *_It = t; // current time
    _pDAG->eval( _opADJTC, _nx, _pADJTC, _Iy, _nx+_npar+1, _pVAR+_nx, _IVAR+_nx );
    _I2vec( _nx, _Iy, y );
    break;

  case OPT::ELLIPS:
  default:
#ifdef MC__ODEBNDS_BASE_DINEQI_DEBUG
    std::cout << "@t=" << t << std::endl;
    //E::options.PSDCHK = true;
    for( unsigned ix=0; ix<_nx; ix++ )
      std::cout << "xref[" << ix << "] = " << _xref[ix] << std::endl;
    for( unsigned ix=0; ix<_nx; ix++ ){
      std::cout << "B[" << ix << ",#] = ";
      for( unsigned jp=0; jp<_npar; jp++ )
        std::cout << _B[jp*_nx+ix] << "  ";
      std::cout << std::endl;
    }
    std::cout << "Er = " << _Er << std::endl;
#endif
    *_MVXPt = t;
    for( unsigned jx=0; jx<_nx; jx++ )
      _MVXPd[jx].set( _MVXPenv, jx, _Ir[jx] );
    _ep2x( _nx, _npar, _MVXPd, _pref, _MVXPp, _B, _xref, _MVXPx );
    _pDAG->eval( _opADJTC, _nx, _pADJTC, _MVXPy, _nx+_npar+1, _pVAR+_nx, _MVXPVAR );
#ifdef MC__ODEBNDS_BASE_DINEQI_DEBUG
    for( unsigned ix=0; ix<_nx; ix++ )
      std::cout << "MVXPy[" << ix << "] = " << _MVXPy[ix] << std::endl;
#endif
    ODEBND_BASE<T,PMT,PVT>::_CC_I_ELL( _nx, _MVXPy, _Iy, _Er, _Ay, _npar, _yref, _By,
      _Idy, _Qy, options.QTOL, machprec() );
#ifdef MC__ODEBNDS_BASE_DINEQI_DEBUG
    for( unsigned ix=0; ix<_nx; ix++ )
      std::cout << "yref[" << ix << "] = " << _yref[ix] << std::endl;
    for( unsigned ix=0; ix<_nx; ix++ ){
      std::cout << "By[" << ix << ",#] = ";
      for( unsigned jp=0; jp<_npar; jp++ )
        std::cout << _By[jp*_nx+ix] << "  ";
      std::cout << std::endl;
    }
    for( unsigned ix=0; ix<_nx; ix++ )
      std::cout << "Iy[" << ix << ",#] = " << _Iy[ix] << std::endl;
    std::cout << "Edy = " << E(_nx,_Qy) << std::endl;
#endif
    _E2vec( _nx, _npar, _yref, _Qy, _By, y );
#ifdef MC__ODEBNDS_BASE_DINEQI_DEBUG
    for( unsigned i=0; i<_nx*(1+_npar)+_nx*(_nx+1)/2; i++ )
      std::cout << "yvec[" << i << "] = " << y[i] << std::endl;
    { int dum; std::cin >> dum; }
#endif
    break;
  }

  return true;
}

template <typename T, typename PMT, typename PVT>
template <typename REALTYPE, typename OPT> inline bool
ODEBNDS_BASE<T,PMT,PVT>::_TC_I_QUAD
( const OPT&options, const double t, REALTYPE*yq )
{
  _opTCQUAD = _pDAG->subgraph( _npar, _pADJTC+_nx );
  switch( options.WRAPMIT){

  case OPT::NONE:
  case OPT::DINEQ:
    _pDAG->eval( _opTCQUAD, _npar, _pADJTC+_nx, _Iyq, _nx+_npar+1, _pVAR+_nx, _IVAR+_nx );
    _I2vec( _npar, _Iyq, yq );
    break;

  case OPT::ELLIPS:
  default:
    _pDAG->eval( _opTCQUAD, _npar, _pADJTC+_nx, _MVXPyq, _nx+_npar+1, _pVAR+_nx, _MVXPVAR );
    _QUAD_I_ELL( _npar, _nx, _MVXPyq, _Byq, _Idyq, _Iyq );
    _I2vec( _npar, _npar, _Byq, _Idyq, yq );
    break;
  }

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
  //std::cout << "In _CC_I_ADJ function" << std::endl;
  const FFVar* pFCT = _vFCT.at(pos_fct-1)+ifct;
#ifndef MC__ODEBNDS_GSL_USE_BAD
  delete[] _pADJTC; _pADJTC = _pDAG->FAD( 1, pFCT, _nx+_npar, _pVAR+_nx );
#else
  delete[] _pADJTC; _pADJTC = _pDAG->BAD( 1, pFCT, _nx+_npar, _pVAR+_nx );
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
#ifdef MC__ODEBNDS_BASE_DINEQI_DEBUG
    std::cout << "@t=" << t << std::endl;
    //E::options.PSDCHK = true;
    for( unsigned ix=0; ix<_nx; ix++ )
      std::cout << "xref[" << ix << "] = " << _xref[ix] << std::endl;
    for( unsigned ix=0; ix<_nx; ix++ ){
      std::cout << "B[" << ix << ",#] = ";
      for( unsigned jp=0; jp<_npar; jp++ )
        std::cout << _B[jp*_nx+ix] << "  ";
      std::cout << std::endl;
    }
    std::cout << "Er = " << _Er << std::endl;
#endif
#ifdef MC__ODEBNDS_BASE_DINEQI_DEBUG
    for( unsigned i=0; i<_nx*(1+_npar)+_nx*(_nx+1)/2+2*_npar; i++ )
      std::cout << "yvec[" << i << "] = " << vec[i] << std::endl;
    //return true;
    for( unsigned ix=0; ix<_nx; ix++ )
      std::cout << "yref[" << ix << "] = " << _yref[ix] << std::endl;
    for( unsigned ix=0; ix<_nx; ix++ ){
      std::cout << "By[" << ix << ",#] = ";
      for( unsigned jp=0; jp<_npar; jp++ )
        std::cout << _By[jp*_nx+ix] << "  ";
      std::cout << std::endl;
    }
    std::cout << "Edy = " << _Edy << std::endl;
    //{ int dum; std::cin >> dum; }
#endif
    // Setup polynomial model expansion of adjoint RHS
    *_MVYXPt = t; // Current time
    for( unsigned jx=0; jx<_nx; jx++ )
      _MVYXPr[jx].set( _MVYXPenv, _nx+jx, _Ir[jx] );
    _ep2x( _nx, _npar, _MVYXPr, _pref, _MVYXPp, _B, _xref, _MVYXPx );
    for( unsigned jy=0; jy<_nx; jy++ )
      _MVYXPd[jy].set( _MVYXPenv, jy, _Idy[jy] );
    _ep2x( _nx, _npar, _MVYXPd, _pref, _MVYXPp, _By, _yref, _MVYXPy );
#ifdef MC__ODEBNDS_BASE_DINEQI_DEBUG
    std::cout << "function #" << ifct << std::endl;
    for( unsigned ix=0; ix<_nx; ix++ )
      std::cout << "MVYXPd[ " << ix << "] = " << _MVYXPd[ix] << std::endl;
    for( unsigned ix=0; ix<_nx; ix++ )
      std::cout << "MVYXPy[ " << ix << "] = " << _MVYXPy[ix] << std::endl;
    for( unsigned ix=0; ix<_nx; ix++ )
      std::cout << "MVYXPr[ " << ix << "] = " << _MVYXPr[ix] << std::endl;
    for( unsigned ix=0; ix<_nx; ix++ )
      std::cout << "MVYXPx[ " << ix << "] = " << _MVYXPx[ix] << std::endl;
    for( unsigned ip=0; ip<_npar; ip++ )
      std::cout << "MVYXPp[ " << ip << "] = " << _MVYXPp[ip] << std::endl;
#endif
    _opADJTC  = _pDAG->subgraph( _nx, _pADJCC );
    delete[] _PMIC; _PMIC = new PVT[_opADJTC.size()];
    _pDAG->eval( _opADJTC, _PMIC, _nx, _pADJCC, _MVYXPg, _nVAR-_npar, _pVAR, _MVYXPVAR );
    _CC_I_ELL( _nx, _MVYXPg, _Iy, _Edy, _Ay, _Er, _A, _npar, _yrefdot, _Bydot, _Idydot,
      _Qydot, options.QTOL, machprec() );
    _E2vec( _nx, _npar, _yrefdot, _Qydot, _Bydot, vec );
#ifdef MC__ODEBNDS_BASE_DINEQI_DEBUG
    for( unsigned i=0; i<_nx*(1+_npar)+_nx*(_nx+1)/2+2*_npar; i++ )
      std::cout << "yvec[" << i << "] = " << vec[i] << std::endl;
    //{ int dum; std::cin >> dum; }
    for( unsigned ix=0; ix<_nx; ix++ )
      std::cout << "yref[" << ix << "] = " << _yrefdot[ix] << std::endl;
    for( unsigned ix=0; ix<_nx; ix++ ){
      std::cout << "By[" << ix << ",#] = ";
      for( unsigned jp=0; jp<_npar; jp++ )
        std::cout << _Bydot[jp*_nx+ix] << "  ";
      std::cout << std::endl;
    }
    std::cout << "Edy = " << E(_nx,_Qydot) << std::endl;
    { int dum; std::cin >> dum; }
#endif
    break;
  }

  return true;
}

template <typename T, typename PMT, typename PVT>
inline void
ODEBNDS_BASE<T,PMT,PVT>::_CC_I_ELL
( const unsigned nx, PVT*MVYXPg, T*Iy, const E&Edy, double*Ay, const E&Edx, double*Ax,
  const unsigned np, double*yref, double*By, T*Idy, double*Qy,
  const double QTOL, const double EPS )
{
  // Extract time derivatives of constant, linear and remainder parts
  for( unsigned ix=0; ix<nx; ix++ ){
    Iy[ix] = MVYXPg[ix].bound();
#ifdef MC__ODEBNDS_BASE_DINEQI_DEBUG
    std::cout << "MVYXPg[" << ix << "] = " << MVYXPg[ix] << std::endl;
#endif
    for( unsigned jx=0; jx<nx; jx++ )
      Ay[ix+jx*nx] = MVYXPg[ix].linear(jx,true);
#ifdef MC__ODEBNDS_BASE_DINEQI_DEBUG
    std::cout << "Ay[" << ix << ",#] = ";
    for( unsigned jx=0; jx<nx; jx++ )
      std::cout << Ay[ix+jx*nx] << "  ";
    std::cout << std::endl;
#endif
    for( unsigned jx=0; jx<nx; jx++ )
      Ax[ix+jx*nx] = MVYXPg[ix].linear(nx+jx,true);
#ifdef MC__ODEBNDS_BASE_DINEQI_DEBUG
    std::cout << "Ax[" << ix << ",#] = ";
    for( unsigned jx=0; jx<nx; jx++ )
      std::cout << Ax[ix+jx*nx] << "  ";
    std::cout << std::endl;
#endif
    for( unsigned jp=0; jp<np; jp++ )
      By[ix+jp*nx] = MVYXPg[ix].linear(2*nx+jp,true);
#ifdef MC__ODEBNDS_BASE_DINEQI_DEBUG
    std::cout << "By[" << ix << ",#] = ";
    for( unsigned jp=0; jp<np; jp++ )
      std::cout << By[ix+jp*nx] << "  ";
    std::cout << std::endl;
    std::cout << "MVYXPg[" << ix << "] = " << MVYXPg[ix] << std::endl;
#endif
    T Rgi = MVYXPg[ix].B();
    yref[ix] = Op<T>::mid(Rgi);
    Idy[ix] = Rgi - Op<T>::mid(Rgi);
#ifdef MC__ODEBNDS_BASE_DINEQI_DEBUG
    std::cout << "yref[" << ix << "] = " << yref[ix] << std::endl;
    std::cout << "Idy[" << ix << "] = " << Idy[ix]
              << " : " << Op<T>::mid(Idy[ix]) << std::endl;
#endif
  }
  CPPL::dgematrix matAy(nx,nx), matAx(nx,nx);
  for( unsigned ix=0; ix<nx; ix++ )
    for( unsigned jx=0; jx<nx; jx++ ){
      matAy(ix,jx) = Ay[ix+jx*nx];
      matAx(ix,jx) = Ax[ix+jx*nx];
    }
#ifdef MC__ODEBNDS_BASE_DINEQI_DEBUG
  std::cout << "Ay·Edy =" << mtimes(Edy,matAy) << std::endl;
  std::cout << "Ax·Edx =" << mtimes(Edx,matAx) << std::endl;
#endif
  E Eg = minksum_ea( mtimes(Edy,matAy), mtimes(Edx,matAx), EPS );
#ifdef MC__ODEBNDS_BASE_DINEQI_DEBUG
  std::cout << "Ay·Edy + Ax·Edx =" << Eg << std::endl;
#endif
  Eg = minksum_ea( Eg, Idy, QTOL, EPS );
#ifdef MC__ODEBNDS_BASE_DINEQI_DEBUG
  std::cout << "Eg =" << Eg << std::endl;
#endif
  for( unsigned jx=0; jx<nx; jx++ )
    for( unsigned ix=jx; ix<nx; ix++ )
      Qy[_ndxLT(ix,jx,nx)] = Eg.Q(ix,jx);
#ifdef MC__ODEBNDS_BASE_DINEQI_DEBUG
  for( unsigned ix=0; ix<nx; ix++ ){
    std::cout << "Qy[" << ix << ",#] = ";
    for( unsigned jx=0; jx<=ix; jx++ )
      std::cout << Qy[_ndxLT(ix,jx,nx)] << "  ";
    std::cout << std::endl;
  }
  { int dum; std::cin >> dum; }
  // throw(0);
#endif
}

template <typename T, typename PMT, typename PVT> 
template <typename REALTYPE, typename OPT> inline bool
ODEBNDS_BASE<T,PMT,PVT>::_CC_I_QUAD
( const OPT&options, const double t, REALTYPE*yq )
{
  _opTCQUAD = _pDAG->subgraph( _npar, _pADJTC+_nx );
  switch( options.WRAPMIT){

  case OPT::NONE:
  case OPT::DINEQ:
    _pDAG->eval( _opTCQUAD, _npar, _pADJTC+_nx, _Iyq, _nVAR-_npar, _pVAR, _IVAR, true );
    _I2vec( _npar, _Iyq, yq );
    break;

  case OPT::ELLIPS:
  default:
    // Construct the adjoint quadrature derivatives
    _pDAG->eval( _opTCQUAD, _npar, _pADJTC+_nx, _MVXPyq, _nVAR-_npar, _pVAR, _MVXPVAR, true );
    _QUAD_I_ELL( _npar, _nx, _MVXPyq, _Byq, _Idyq, _Iyq );
    _I2vec( _npar, _npar, _Byq, _Idyq, yq );
    break;
  }

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
      if( PMopmax < _opADJRHS[ifct].size() )  PMopmax = _opADJRHS[ifct].size();
      if( PMopmax < _opADJQUAD[ifct].size() ) PMopmax = _opADJQUAD[ifct].size();
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
  REALTYPE* x, const unsigned ifct, const bool neg )
{
  if( !_pADJRHS ) return false; // **error** ADJRHS not defined

  switch( options.WRAPMIT){
  case OPT::NONE:
    _vec2I( x, _nx, _Ix ); // set current state bounds
    _vec2I( y, _nx, _Iy ); // set current adjoint bounds
    *_It = t; // set current time
    _pDAG->eval( _opADJRHS[ifct], _IADJRHS, _nx, _pADJRHS+ifct*_nx, _Iydot,
                 _nVAR-_npar, _pVAR, _IVAR );
    if( !neg )
      _I2vec( _nx, _Iydot, ydot );
    else{
      for( unsigned ix=0; ix<_nx; ix++ ){
        _yUdot[ix] = Op<T>::l( _Iydot[ix] );
        _yLdot[ix] = Op<T>::u( _Iydot[ix] );
      }
      _I2vec( _nx, _yLdot, _yUdot, ydot );
    }
    break;  
   
  case OPT::DINEQ:
    _vec2I( x, _nx, _Ix );  // set current state bounds
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
    _vec2E( x, _nx, _npar, _Q, _Er, _Ir, _pref, _Ip, _B, _xref, _Ix ); // set current state bounds
    _vec2E( y, _nx, _npar, _Qy, _Edy, _Idy, _pref, _Ip, _By, _yref, _Iy); // set current adjoint bounds
#ifdef MC__ODEBNDS_BASE_DINEQI_DEBUG
    std::cout << "@t=" << t << std::endl;
    //E::options.PSDCHK = true;
    for( unsigned ix=0; ix<_nx; ix++ )
      std::cout << "xref[" << ix << "] = " << _xref[ix] << std::endl;
    for( unsigned ix=0; ix<_nx; ix++ ){
      std::cout << "B[" << ix << ",#] = ";
      for( unsigned jp=0; jp<_npar; jp++ )
        std::cout << _B[jp*_nx+ix] << "  ";
      std::cout << std::endl;
    }
    std::cout << "Er = " << _Er << std::endl;
    for( unsigned i=0; i<_nx*(1+_npar)+_nx*(_nx+1)/2; i++ )
      std::cout << "yvec[" << i << "] = " << y[i] << std::endl;
    for( unsigned ix=0; ix<_nx; ix++ )
      std::cout << "yref[" << ix << "] = " << _yref[ix] << std::endl;
    for( unsigned ix=0; ix<_nx; ix++ ){
      std::cout << "By[" << ix << ",#] = ";
      for( unsigned jp=0; jp<_npar; jp++ )
        std::cout << _By[jp*_nx+ix] << "  ";
      std::cout << std::endl;
    }
    std::cout << "Edy = " << _Edy << std::endl;
    //{ int dum; std::cin >> dum; }
#endif
    // Setup polynomial model expansion of adjoint RHS
    *_MVYXPt = t; // Current time
    for( unsigned jx=0; jx<_nx; jx++ ){
      _MVYXPw[jx].set( _MVYXPenv, _nx+jx, T(-1.,1.) );
      _MVYXPr[jx] = _MVYXPw[jx] * 0.5*Op<T>::diam(_Ir[jx]);
    }
    _ep2x( _nx, _npar, _MVYXPr, _pref, _MVYXPp, _B, _xref, _MVYXPx );
    for( unsigned jy=0; jy<_nx; jy++ )
      _MVYXPd[jy].set( _MVYXPenv, jy, _Idy[jy] );
    _ep2x( _nx, _npar, _MVYXPd, _pref, _MVYXPp, _By, _yref, _MVYXPy );
#ifdef MC__ODEBNDS_BASE_DINEQI_DEBUG
    std::cout << "function #" << ifct << std::endl;
    for( unsigned ix=0; ix<_nx; ix++ )
      std::cout << "MVYXPy[ " << ix << "] = " << _MVYXPy[ix] << std::endl;
    for( unsigned ix=0; ix<_nx; ix++ )
      std::cout << "MVYXPx[ " << ix << "] = " << _MVYXPx[ix] << std::endl;
    for( unsigned ip=0; ip<_npar; ip++ )
      std::cout << "MVYXPp[ " << ip << "] = " << _MVYXPp[ip] << std::endl;
#endif
    // Construct the adjoint ellipsoidal remainder derivatives
    _pDAG->eval( _opADJRHS[ifct], _PMADJRHS, _nx, _pADJRHS+ifct*_nx, _MVYXPg,
                 _nVAR-_npar, _pVAR, _MVYXPVAR );
    if( options.QSCALE )
      _RHS_I_ELL( _nx, _MVYXPg, _Qy, _Ay, _Ax, _npar, _yrefdot, _Bydot, _Idydot, _Qydot, 
                    options.QTOL, machprec(), neg, _Iy );
    else
      _RHS_I_ELL( _nx, _MVYXPg, _Qy, _Ay, _Ax, _npar, _yrefdot, _Bydot, _Idydot, _Qydot, 
                    options.QTOL, machprec(), neg );
    _E2vec( _nx, _npar, _yrefdot, _Qydot, _Bydot, ydot );
    break;
  }
  return true;
}

template <typename T, typename PMT, typename PVT>
inline void
ODEBNDS_BASE<T,PMT,PVT>::_RHS_I_ELL
( const unsigned nx, PVT*MVYXPg, const double*Qy, double*Ay, double *Ax,
  const unsigned np, double*yrefdot, double*Bydot, T*Idydot, double*Qydot,
  const double QTOL, const double EPS, const bool neg, const T*W )
{
  // Extract time derivatives of constant, linear and remainder parts
  // Set reference and linear block RHS
  for( unsigned ix=0; ix<nx; ix++ ){
#ifdef MC__ODEBNDS_BASE_DINEQI_DEBUG
    std::cout << "MVYXPg[" << ix << "] = " << MVYXPg[ix] << std::endl;
#endif
    for( unsigned jx=0; jx<nx; jx++ )
      Ay[ix+jx*nx] = MVYXPg[ix].linear(jx,true);
#ifdef MC__ODEBNDS_BASE_DINEQI_DEBUG
    std::cout << "Ay[" << ix << ",#] = ";
    for( unsigned jx=0; jx<nx; jx++ )
      std::cout << Ay[ix+jx*nx] << "  ";
    std::cout << std::endl;
#endif
    for( unsigned jx=0; jx<nx; jx++ )
      Ax[ix+jx*nx] = MVYXPg[ix].linear(nx+jx,true);
#ifdef MC__ODEBNDS_BASE_DINEQI_DEBUG
    std::cout << "Ax[" << ix << ",#] = ";
    for( unsigned jx=0; jx<nx; jx++ )
      std::cout << Ax[ix+jx*nx] << "  ";
    std::cout << std::endl;
#endif
    for( unsigned jp=0; jp<np; jp++ )
      Bydot[ix+jp*nx] = MVYXPg[ix].linear(2*nx+jp,true);
#ifdef MC__ODEBNDS_BASE_DINEQI_DEBUG
    std::cout << "Bydot[" << ix << ",#] = ";
    for( unsigned jp=0; jp<np; jp++ )
      std::cout << Bydot[ix+jp*nx] << "  ";
    std::cout << std::endl;
    std::cout << "MVYXPg[" << ix << "] = " << MVYXPg[ix] << std::endl;
#endif
    T Rgi = MVYXPg[ix].B();
    yrefdot[ix] = Op<T>::mid(Rgi);
    Idydot[ix] = Rgi - Op<T>::mid(Rgi);
#ifdef MC__ODEBNDS_BASE_DINEQI_DEBUG
    std::cout << "xrefdot[" << ix << "] = " << xrefdot[ix] << std::endl;
    std::cout << "Iddot[" << ix << "] = " << Iddot[ix]
              << " : " << Op<T>::mid(Iddot[ix]) << std::endl;
#endif
  }

  return _RHS_ELL( nx, Qy, Ay, Ax, Idydot, Qydot, QTOL, EPS, neg, W );
}

template <typename T, typename PMT, typename PVT>
template <typename U>
inline void
ODEBNDS_BASE<T,PMT,PVT>::_RHS_ELL
( const unsigned nx, const double*Qy, const double*Ay, const double *Ax,
  const T*Idydot, double*Qydot, const double QTOL, const double EPS,
  const bool neg, const U*W )
{
  // Construct trajectories kappa and eta
  const double WTOL = 1e-8;//EPS*1e2;
  double trQ = 0.;
  for( unsigned ix=0; ix<nx; ix++ ){
    double sqr_wi = W? sqr(Op<U>::abs(W[ix])+WTOL): 1.;
    trQ += ( Qy[_ndxLT(ix,ix,nx)]>EPS? Qy[_ndxLT(ix,ix,nx)]/sqr_wi: EPS );
  }
  double sumkappa = 0., sumeta = 0., AxTAx[nx];
  const double srqt_trQ = (trQ>0? std::sqrt( trQ ): 0.) + QTOL;
  for( unsigned ix=0; ix<nx; ix++ ){
    double wi = W? Op<U>::abs(W[ix])+WTOL: 1.;
    //sumkappa += .1;
    sumkappa += Op<T>::diam( Idydot[ix] ) / ( 2. * wi * srqt_trQ );
    AxTAx[ix] = 0.;
    for( unsigned jx=0; jx<nx; jx++)
      AxTAx[ix] += Ax[jx+ix*nx] * Ax[jx+ix*nx];
    if( AxTAx[ix] < EPS ) AxTAx[ix] = EPS;
    sumeta += std::sqrt(AxTAx[ix]) / ( wi * srqt_trQ );
#ifdef MC__ODEBNDS_BASE_DINEQPM_DEBUG
    std::cout << "eta[" << ix << "] = "
              << std::sqrt(AxTAx[ix]) / ( wi * srqt_trQ )
              << std::endl;
    std::cout << "kappa[" << ix << "] = "
              << Op<T>::diam( Idydot[ix] ) / ( 2. * wi * srqt_trQ )
              << std::endl;
#endif
  }

  // Set ellipsoidal remainder RHS
  double pm = neg? -1.: 1.;
  for( unsigned jx=0; jx<nx; jx++ ){
    for( unsigned ix=jx; ix<nx; ix++ ){
      Qydot[_ndxLT(ix,jx,nx)] = pm * ( sumkappa + sumeta ) * Qy[_ndxLT(ix,jx,nx)];
      for( unsigned kx=0; kx<nx; kx++ ){
        double wk = W? Op<U>::abs(W[kx])+WTOL: 1.,
               wksqtrQ = wk * ( std::sqrt( trQ ) + QTOL );
        Qydot[_ndxLT(ix,jx,nx)] += Qy[_ndxLT(ix,kx,nx)] * Ay[jx+kx*nx]
          + Ay[ix+kx*nx] * Qy[_ndxLT(kx,jx,nx)]
          + pm * Ax[ix+kx*nx] * Ax[jx+kx*nx] / std::sqrt(AxTAx[kx]) * wksqtrQ;
      }
    }
    double wj = W? Op<U>::abs(W[jx])+WTOL: 1.;
    Qydot[_ndxLT(jx,jx,nx)] += pm * Op<T>::diam( Idydot[jx] ) / 2. * wj * srqt_trQ;
    //Qydot[_ndxLT(jx,jx,nx)] += pm * sqr( Op<T>::diam( Idydot[jx] ) / 2. ) / 0.1;
  }
#ifdef MC__ODEBNDS_BASE_DINEQI_DEBUG
  for( unsigned ix=0; ix<nx; ix++ ){
    std::cout << "Qydot[" << ix << ",#] = ";
    for( unsigned jx=0; jx<=ix; jx++ )
      std::cout << Qydot[_ndxLT(ix,jx,nx)] << "  ";
    std::cout << std::endl;
  }
  E Eydot( nx, Qydot, yrefdot );
  std::cout << "Eydot =" << Eydot << std::endl;
  { int dum; std::cin >> dum; }
#endif
}

template <typename T, typename PMT, typename PVT>
template <typename REALTYPE, typename OPT> inline bool
ODEBNDS_BASE<T,PMT,PVT>::_RHS_I_QUAD
( const OPT&options, double t, const REALTYPE*y, REALTYPE*qdot,
  REALTYPE*x, const unsigned ifct, const bool neg )
{
  if( !_pADJQUAD ) return false;

  switch( options.WRAPMIT){
  case OPT::NONE:
  case OPT::DINEQ:
    _pDAG->eval( _opADJQUAD[ifct], _IADJRHS, _npar, _pADJQUAD+ifct*_npar,
                 _Iyqdot, _nVAR-_npar, _pVAR, _IVAR );
    if( !neg )
      _I2vec( _npar, _Iyqdot, qdot );
    else{
      for( unsigned ip=0; ip<_npar; ip++ ){
        _yqUdot[ip] = Op<T>::l( _Iyqdot[ip] );
        _yqLdot[ip] = Op<T>::u( _Iyqdot[ip] );
      }
      _I2vec( _npar, _yqLdot, _yqUdot, qdot );
    }
    break;

  case OPT::ELLIPS:
  default:
    _pDAG->eval( _opADJQUAD[ifct], _PMADJRHS, _npar, _pADJQUAD+ifct*_npar,
                 _MVYXPyqdot, _nVAR-_npar, _pVAR, _MVYXPVAR );
    _QUAD_I_ELL( _npar, 2*_nx, _MVYXPyqdot, _Byqdot, _Idyqdot, _Iyqdot );
    if( !neg )
      _I2vec( _npar, _npar, _Byqdot, _Idyqdot, qdot );
    else{
      for( unsigned ip=0; ip<_npar; ip++ ){
        _yqUdot[ip] = Op<T>::l( _Idyqdot[ip] );
        _yqLdot[ip] = Op<T>::u( _Idyqdot[ip] );
      }
      _I2vec( _npar, _npar, _Byqdot, _yqLdot, _yqUdot, qdot );
    }
    break;
  }

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
  _nVAR = _nx + _nx + _npar + 1 + _nq + _npar;
  delete[] _pVAR; _pVAR = new FFVar[_nVAR];
  delete[] _IVAR; _IVAR = new T[_nVAR];
  delete[] _PMVAR; _PMVAR = new PVT[_nVAR];
  delete[] _pADJCC; _pADJCC = new FFVar[_nx];

  BASE_DE::set_adjoint();
  for( unsigned ix=0; ix<_nx; ix++ ) _pVAR[ix] = _pY[ix];
  for( unsigned ix=0; ix<_nx; ix++ ) _pVAR[_nx+ix] = _pX[ix];
  for( unsigned ip=0; ip<_npar; ip++ ) _pVAR[2*_nx+ip] = _pP[ip];
  _pVAR[2*_nx+_npar] = (_pT? *_pT: 0. );
  _Iy = _IVAR;
  _Ix = _Iy + _nx;
  _Ip = _Ix + _nx;
  _It = _Ip + _npar;
  _Iyq = 0; // not used
  for( unsigned ip=0; ip<_npar; ip++ ) _Ip[ip] = PMp[ip].bound();
  _PMy = _PMVAR;
  _PMx = _PMy + _nx;
  _PMp = _PMx + _nx;
  _PMt = _PMp + _npar;
  _PMyq = _PMt + 1;
  for( unsigned ip=0; ip<_npar; ip++ ) _PMp[ip] = PMp[ip];
  
  // Reset _MVYXPenv and related variables
  unsigned MVYZsize = ( options.ORDMIT<_PMenv->nord()? options.ORDMIT: _PMenv->nord() ); 
#ifdef MC__ODEBNDS_BASE_MVYXP_USE
  unsigned MVYZdim  = 2*_nx+_npar;
#else
  unsigned MVYZdim  = ( options.ORDMIT<_PMenv->nord()? _nx+_npar: 2*_nx+_npar ); 
#endif
  if( _MVYXPenv && ( _MVYXPenv->nord() != MVYZsize
                  || _MVYXPenv->nvar() != MVYZdim  ) ){
    delete[] _MVYXPg;    _MVYXPg = 0;
    delete[] _MVYXPdgdy; _MVYXPdgdy = 0;
    delete[] _MVYXPw;    _MVYXPw = 0;
    delete[] _MVYXPr;    _MVYXPr = 0;
    delete[] _MVYXPd;    _MVYXPd = 0;
    delete[] _MVYXPVAR;  _MVYXPVAR = _MVYXPy = _MVYXPx = _MVYXPp = _MVYXPt = 0; 
    delete   _MVYXPenv;  _MVYXPenv = 0;
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
    delete[] _yref;     _yref     = new double[_nx];
    delete[] _Ay;       _Ay       = new double[_nx*_nx];
    delete[] _Ax;       _Ax       = new double[_nx*_nx];
    delete[] _Qy;       _Qy       = new double[_nx*(_nx+1)/2];
    delete[] _Idy;      _Idy      = new T[_nx];
    delete[] _Qydot;    _Qydot    = new double[_nx*(_nx+1)/2];
    delete[] _Idydot;   _Idydot   = new T[_nx];
    delete[] _Idgdy;    _Idgdy    = new T[2*_nx*_nx];

    if( !_MVYXPenv ) _MVYXPenv    = new PMT( MVYZdim, MVYZsize );
    _MVYXPenv->options = _PMenv->options;
    if( !_MVYXPw )    _MVYXPw     = new PVT[_nx];
    if( !_MVYXPr )    _MVYXPr     = new PVT[_nx];
    if( !_MVYXPd )    _MVYXPd     = new PVT[_nx];
    if( !_MVYXPg )    _MVYXPg     = new PVT[_nx];
    if( !_MVYXPdgdy ) _MVYXPdgdy  = new PVT[2*_nx*_nx];
    if( !_MVYXPVAR )  _MVYXPVAR   = new PVT[_nVAR];
    _MVYXPy = _MVYXPVAR;
    _MVYXPx = _MVYXPy + _nx;
    _MVYXPp = _MVYXPx + _nx;
    _MVYXPt = _MVYXPp + _npar;
    for( unsigned ip=0; ip<_npar; ip++ )
      _MVYXPp[ip].set( _MVYXPenv, ip, _PMp[ip].B() );
    break;
  }

  return true;
}

template <typename T, typename PMT, typename PVT>
template <typename REALTYPE, typename OPT> inline bool
ODEBNDS_BASE<T,PMT,PVT>::_TC_PM_ADJ
( const OPT &options, const double t, const REALTYPE*x, REALTYPE*y,
  unsigned pos_fct, unsigned ifct )
{
  const FFVar* pFCT = _vFCT.at(pos_fct)+ifct;
#ifndef MC__ODEBNDS_GSL_USE_BAD
  delete[] _pADJTC; _pADJTC = _pDAG->FAD( 1, pFCT, _nx+_npar, _pVAR+_nx );
#else
  delete[] _pADJTC; _pADJTC = _pDAG->BAD( 1, pFCT, _nx+_npar, _pVAR+_nx );
#endif
  _opADJTC = _pDAG->subgraph( _nx, _pADJTC );
  delete[] _pADJDTC; _pADJDTC = 0;
  delete[] _PMADJTC; _PMADJTC = new PVT[_opADJTC.size()];
  _opADJDTC.clear();
  delete[] _PMADJDTC; _PMADJDTC = 0;
  delete[] _IADJDTC;  _IADJDTC  = 0;

  switch( options.WRAPMIT){
  case OPT::NONE:
    _vec2PMI( x, _PMenv, _nx, _PMx, true );
    *_PMt = t; // current time
    _pDAG->eval( _opADJTC, _PMADJTC, _nx, _pADJTC, _PMy, _nx+_npar+1, _pVAR+_nx, _PMVAR+_nx );
    // Whether or not to ignore the remainder
    if( !options.PMNOREM )
      _PMI2vec( _PMenv, _nx, _PMy, y, true );
    else
      _PMI2vec( _PMenv, _nx, _PMy, 0, y );
    break;

  case OPT::DINEQ:
    _vec2PMI( x, _PMenv, _nx, _PMx );
    *_PMt = t; // current time
    _pDAG->eval( _opADJTC, _PMADJTC, _nx, _pADJTC, _PMy, _nx+_npar+1, _pVAR+_nx, _PMVAR+_nx );
    // Whether or not to ignore the remainder
    if( !options.PMNOREM )
      _PMI2vec( _PMenv, _nx, _PMy, y );
    else
      _PMI2vec( _PMenv, _nx, _PMy, 0, y );
    break;
   
  case OPT::ELLIPS:
  default:{
    _vec2PME( x, _PMenv, _nx, _PMx, _Q, _Er, _Ir );
    *_PMt = t; // current time   

    // In this variant a bound on the Jacobian matrix is computed and the
    // linear part is taken as the mid-point of this matrix
    if( !options.ORDMIT ){
      *_It = t; // current time
      for( unsigned ix=0; ix<_nx; ix++ )
        _PMx[ix].center().set( T(0.) ); // cancel remainder term
      _pADJDTC = _pDAG->FAD( _nx, _pADJTC, _nx, _pVAR+_nx );
      _opADJDTC = _pDAG->subgraph( _nx*_nx, _pADJDTC );
      _IADJDTC = new T[ _opADJDTC.size() ];
      _pDAG->eval( _opADJTC, _PMADJTC, _nx, _pADJTC, _PMy, _nx+_npar+1, _pVAR+_nx, _PMVAR+_nx );
      _pDAG->eval( _opADJDTC, _IADJDTC, _nx*_nx, _pADJDTC, _Idfdx, _nx+_npar+1, _pVAR+_nx, _IVAR+_nx );
      ODEBND_BASE<T,PMT,PVT>::_RHS_PM_ELL0( _nx, _PMy, _Idfdx, _Ir, _Ax, _Idy );
    }

    // In this variant a polynomial model of the Jacobian matrix is computed and the
    // linear part is taken as the mid-point of this matrix
    else if( _PMenv->nord() >= _MVXPenv->nord() ){
      *_MVXPt = t; // current time
#ifdef MC__ODEBND_BASE_MVXP_USE
      for( unsigned jx=0; jx<_nx; jx++ ){
        _MVXPd[jx].set( _MVXPenv, _npar+jx, _Ir[jx] );
        _MVXPx[jx].set( _MVXPenv ).set( _PMx[jx].center().set( T(0.) ), true );
      }
      _e2x( _nx, _MVXPd, _MVXPx, false );
#else
      for( unsigned jx=0; jx<_nx; jx++ ){
        _MVXPx[jx].set( _MVXPenv ).set( _PMx[jx].center(), true );
        _PMx[jx].set( T(0.) );
      }
#endif
      _pADJDTC = _pDAG->FAD( _nx, _pADJTC, _nx, _pVAR+_nx );
      _opADJDTC = _pDAG->subgraph( _nx*_nx, _pADJDTC );
      _PMADJDTC = new PVT[ _opADJDTC.size() ];
      _pDAG->eval( _opADJTC, _PMADJTC, _nx, _pADJTC, _PMy, _nx+_npar+1, _pVAR+_nx, _PMVAR+_nx );
      _pDAG->eval( _opADJDTC, _PMADJDTC, _nx*_nx, _pADJDTC, _MVXPdfdx, _nx+_npar+1, _pVAR+_nx, _MVXPVAR );
      ODEBND_BASE<T,PMT,PVT>::_RHS_PM_ELL1( _nx, _PMy, _MVXPdfdx, _Ir, _Ax, _Idy );
    }

    // In this variant a polynomial model in the joint state-parameter and
    // of the same order as the parameter polynomial model is computed
    else{
      *_MVXPt = t; // current time   
      for( unsigned jx=0; jx<_nx; jx++ ){
        _MVXPd[jx].set( _MVXPenv, _npar+jx, _Ir[jx] );
        _MVXPx[jx].set( _MVXPenv ).set( _PMx[jx].center().set( T(0.) ), true );
      }
      _e2x( _nx, _MVXPd, _MVXPx, false );
      _pDAG->eval( _opADJTC, _nx, _pADJTC, _MVXPf, _nx+_npar+1, _pVAR+_nx, _MVXPVAR );
      ODEBND_BASE<T,PMT,PVT>::_RHS_PM_ELL2( _nx, _PMenv, _PMy, _MVXPf, _npar, _Ir, _Ax, _Idy );
    }

    // Whether or not to ignore the remainder
    if( !options.PMNOREM ){
      ODEBND_BASE<T,PMT,PVT>::_CC_PM_ELL( _nx, _Er, _Ax, _Idy, _Qy, options.QTOL, machprec() );
      _PME2vec( _PMenv, _nx, _PMy, _Qy, y );
    }
    else
      _PME2vec( _PMenv, _nx, _PMy, 0, y );
    break;
   }
  }

  return true;
}

template <typename T, typename PMT, typename PVT>
template <typename REALTYPE, typename OPT> inline bool
ODEBNDS_BASE<T,PMT,PVT>::_TC_PM_QUAD
( const OPT &options, const double t, REALTYPE*yq )
{
  // THIS FUNCTION MUST BE CALLED AFTER _TC_PM_ADJ
  *_PMt = t; // current time
  _opTCQUAD = _pDAG->subgraph( _npar, _pADJTC+_nx );
  _pDAG->eval( _opTCQUAD, _npar, _pADJTC+_nx, _PMyq, _nx+_npar+1, _pVAR+_nx, _PMVAR+_nx );

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
  delete[] _pADJTC; _pADJTC = _pDAG->FAD( 1, pFCT, _nx+_npar, _pVAR+_nx );
#else
  delete[] _pADJTC; _pADJTC = _pDAG->BAD( 1, pFCT, _nx+_npar, _pVAR+_nx );
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
    _print_interm( t, _nx+_npar, _PMx, _Er, "PMz Intermediate", std::cerr );
    _print_interm( t, _nx, _PMy, _Edy, "PMy Intermediate", std::cerr );
#endif
    return true;

    // In this variant a bound on the Jacobian matrix is computed and the
    // linear part is taken as the mid-point of this matrix
    if( !options.ORDMIT ){
      *_It = t; // set current time
      for( unsigned ix=0; ix<_nx; ix++ ){
        _Ix[ix] = _PMx[ix].bound(); // set current state bounds
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
      ODEBND_BASE<T,PMT,PVT>::_RHS_PM_ELL0( _nx, _PMydot, _Idgdy, _Idy, _Ay, _Idydot );
    }

    // In this variant a polynomial model of the Jacobian matrix is computed and the
    // linear part is taken as the mid-point of this matrix
    else if( _PMenv->nord() > _MVYXPenv->nord() ){
      // Setup polynomial model expansion of adjoint RHS
      *_MVYXPt = t; // time
      for( unsigned jx=0; jx<_nx; jx++ ) // state bounds
        _MVYXPx[jx].set( _MVYXPenv ).set( _PMx[jx].center(), true );
#ifdef MC__ODEBNDS_BASE_MVYXP_USE
      for( unsigned jx=0; jx<_nx; jx++ ){ // adjoint bounds
        _MVYXPd[jx].set( _MVYXPenv, _npar+jx, _Idy[jx] );
        _MVYXPy[jx].set( _MVYXPenv ).set( _PMy[jx].center().set( T(0.) ), true );
      }
      _e2x( _nx, _MVYXPd, _MVYXPy, false );
#else
      for( unsigned jx=0; jx<_nx; jx++ ){
        _MVYXPy[jx].set( _MVYXPenv ).set( _PMy[jx].center(), true );
        _PMy[jx].set( T(0.) );
      }
#endif
      _PMADJTC = new PVT[_opADJTC.size()];
      _pADJDTC = _pDAG->FAD( _nx, _pADJCC, _nx, _pVAR ); // This is the identity matrix!!
      _opADJDTC = _pDAG->subgraph( _nx*_nx, _pADJDTC );
      _IADJDTC = new T[ _opADJDTC.size() ];
      _pDAG->eval( _opADJTC, _PMADJTC, _nx, _pADJCC, _PMydot, _nVAR-_npar,
                   _pVAR, _PMVAR );
      _pDAG->eval( _opADJDTC, _PMADJDTC, _nx*_nx, _pADJDTC, _MVYXPdgdy, _nVAR-_npar,
                   _pVAR, _MVYXPVAR );
      ODEBND_BASE<T,PMT,PVT>::_RHS_PM_ELL1( _nx, _PMydot, _MVYXPdgdy, _Idy, _Ay, _Idydot );
    }

    // In this variant a polynomial model in the joint adjoint-parameter and
    // of the same order as the parameter polynomial model is computed
    else{
      //{std::cout << "here!\n";  int dum; std::cin >> dum; }
      // Setup polynomial model expansion of adjoint RHS
      *_MVYXPt = t; // set current time 
      for( unsigned jx=0; jx<_nx; jx++ ) // state bounds
        _MVYXPx[jx].set( _MVYXPenv ).set( _PMx[jx].center(), true );
      for( unsigned jx=0; jx<_nx; jx++ ){ // adjoint bounds
        _MVYXPd[jx].set( _MVYXPenv, _npar+jx, _Idy[jx] );
        _MVYXPy[jx].set( _MVYXPenv ).set( _PMy[jx].center().set( T(0.) ), true );
      }
      _e2x( _nx, _MVYXPd, _MVYXPy, false );
      _PMADJTC = new PVT[_opADJTC.size()];
      _pDAG->eval( _opADJTC, _PMADJTC, _nx, _pADJCC, _MVYXPg, _nVAR-_npar,
                   _pVAR, _MVYXPVAR );
      ODEBND_BASE<T,PMT,PVT>::_RHS_PM_ELL2( _nx, _PMenv, _PMydot, _MVYXPg, _npar, _Idy, _Ay, _Idydot );
    }

    // Whether or not to ignore the remainder
    if( !options.PMNOREM ){
      ODEBND_BASE<T,PMT,PVT>::_CC_PM_ELL( _nx, _Edy, _Ay, _Idydot, _Qydot, options.QTOL, machprec() );
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
    if( options.ORDMIT && _PMenv->nord() <= _MVYXPenv->nord() ) break;
    _opADJJAC  = new std::list<const FFOp*>[_nf];
    _pADJJAC = _pDAG->FAD( _nx*_nf, _pADJRHS, 2*_nx, _pVAR );
    for( unsigned ifct=0; ifct<_nf; ifct++ ){
      //_pADJJAC = _pDAG->FAD( 2*_nx, _pADJRHS[ifct], 2*_nx, _pVAR );
      _opADJJAC[ifct] = _pDAG->subgraph( 2*_nx*_nx, _pADJJAC+ifct*2*_nx*_nx );
      if( JACopmax < _opADJJAC[ifct].size() ) JACopmax = _opADJJAC[ifct].size();
    }
    if( !options.ORDMIT )
      _IADJJAC  = JACopmax?  new T[ JACopmax ]: 0;
    else if( _PMenv->nord() > _MVYXPenv->nord() )
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
    _vec2PMI( x, _PMenv, _nx, _PMx, true );// set current state bounds
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
    _vec2PMI( x, _PMenv, _nx, _PMx );// set current state bounds
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
    if( !options.PMNOREM )
      _PMI2vec( _PMenv, _nx, _PMydot, _RyLdot, _RyUdot, ydot );
    else
      _PMI2vec( _PMenv, _nx, _PMydot, 0, ydot );
#ifdef MC__ODEBNDS_BASE_DINEQPM_DEBUG
    _print_interm( t, _nx, _PMydot, "PMydot Intermediate", std::cerr );
    _print_interm( t, _nx, _RyLdot, "RyLdot Intermediate", std::cerr );
    _print_interm( t, _nx, _RyUdot, "RyUdot Intermediate", std::cerr );
    { std::cout << "--paused--"; int dum; std::cin >> dum; }
#endif
    break;

  case OPT::ELLIPS:
  default:
    *_PMt = t; // current time
    _vec2PME( x, _PMenv, _nx, _PMx, _Q, _Er, _Ir );// set current state bounds
#ifdef MC__ODEBNDS_BASE_DINEQPM_DEBUG
    std::cout << "Function:" << ifct << std::endl;
    _print_interm( t, _nx, _PMx, _Er, "PMx Intermediate", std::cerr );
#endif
    _vec2PME( y, _PMenv, _nx, _PMy, _Qy, _Edy, _Idy ); // set current adjoint bounds
#ifdef MC__ODEBNDS_BASE_DINEQPM_DEBUG
    _print_interm( t, _nx, _PMy, _Edy, "PMy Intermediate", std::cerr );
#endif

    // In this variant a bound on the Jacobian matrix is computed and the
    // linear part is taken as the mid-point of this matrix
    if( !options.ORDMIT ){
      *_It = t; // set current time
      for( unsigned ix=0; ix<_nx; ix++ ){
        _Ix[ix] = _PMx[ix].bound(); // set current state bounds
        _Iy[ix] = _PMy[ix].bound(); // set current adjoint bounds
        ( _PMx[ix].center() ).set( T(0.) ); // cancel state remainder term
        ( _PMy[ix].center() ).set( T(0.) ); // cancel adjoint remainder term
      }
      _pDAG->eval( _opADJRHS[ifct], _PMADJRHS, _nx, _pADJRHS+ifct*_nx, _PMydot,
                   _nVAR-_npar, _pVAR, _PMVAR );
      _pDAG->eval( _opADJJAC[ifct], _IADJJAC, 2*_nx*_nx, _pADJJAC+ifct*2*_nx*_nx,
                   _Idgdy, _nVAR-_npar, _pVAR, _IVAR );
      _RHS_PM_ELL0( _nx, _PMydot, _Idgdy, _Ir, _Ax, _Idy, _Ay, _Idydot );
    }

    // In this variant a polynomial model of the Jacobian matrix is computed and the
    // linear part is taken as the mid-point of this matrix
    else if( _PMenv->nord() > _MVYXPenv->nord() ){
      // Setup polynomial model expansion of adjoint RHS
      *_MVYXPt = t; // time
//#ifdef MC__ODEBNDS_BASE_MVYXP_USE
      for( unsigned jx=0; jx<_nx; jx++ ){ // state bounds
        _MVYXPw[jx].set( _MVYXPenv, _npar+jx, T(-1.,1.) );
        _MVYXPr[jx] = _MVYXPw[jx] * 0.5*Op<T>::diam(_Ir[jx]);
        _MVYXPx[jx].set( _MVYXPenv ).set( _PMx[jx].center().set( T(0.) ), true );
      }
      _e2x( _nx, _MVYXPr, _MVYXPx, false );
//#else
//      for( unsigned jx=0; jx<_nx; jx++ ){
//        _MVYXPx[jx].set( _MVYXPenv ).set( _PMx[jx].center(), true );
//        _PMx[jx].set( T(0.) );
//      }
//#endif
#ifdef MC__ODEBNDS_BASE_DINEQPM_DEBUG
      _print_interm( t, _nx, _MVYXPx, "MVYXPx Intermediate", std::cerr );
#endif
#ifdef MC__ODEBNDS_BASE_MVYXP_USE
      for( unsigned jx=0; jx<_nx; jx++ ){ // adjoint bounds
        _MVYXPd[jx].set( _MVYXPenv, _npar+_nx+jx, _Idy[jx] );
        _MVYXPy[jx].set( _MVYXPenv ).set( _PMy[jx].center().set( T(0.) ), true );
      }
      _e2x( _nx, _MVYXPd, _MVYXPy, false );
#else
      for( unsigned jx=0; jx<_nx; jx++ ){
        _MVYXPy[jx].set( _MVYXPenv ).set( _PMy[jx].center(), true );
        _PMy[jx].set( T(0.) );
      }
#endif
#ifdef MC__ODEBNDS_BASE_DINEQPM_DEBUG
      _print_interm( t, _nx, _MVYXPy, "MVYXPy Intermediate", std::cerr );
#endif
      _pDAG->eval( _opADJRHS[ifct], _PMADJRHS, _nx, _pADJRHS+ifct*_nx, _PMydot,
                   _nVAR-_npar, _pVAR, _PMVAR );
      _pDAG->eval( _opADJJAC[ifct], _PMADJJAC, 2*_nx*_nx, _pADJJAC+ifct*2*_nx*_nx,
                   _MVYXPdgdy, _nVAR-_npar, _pVAR, _MVYXPVAR );
      _RHS_PM_ELL1( _nx, _PMydot, _MVYXPdgdy, _Ir, _Ax, _Idy, _Ay, _Idydot );
    }

    // In this variant a polynomial model in the joint adjoint-parameter and
    // of the same order as the parameter polynomial model is computed
    else{
      *_MVYXPt = t; // set current time 
      for( unsigned jx=0; jx<_nx; jx++ ){ // state bounds
        _MVYXPw[jx].set( _MVYXPenv, _npar+jx, T(-1.,1.) );
        _MVYXPr[jx] = _MVYXPw[jx] * 0.5*Op<T>::diam(_Ir[jx]);
        _MVYXPx[jx].set( _MVYXPenv ).set( _PMx[jx].center().set( T(0.) ), true );
      }
      _e2x( _nx, _MVYXPr, _MVYXPx, false );
      for( unsigned jx=0; jx<_nx; jx++ ){ // adjoint bounds
        _MVYXPd[jx].set( _MVYXPenv, _npar+_nx+jx, _Idy[jx] );
        _MVYXPx[jx].set( _MVYXPenv ).set( _PMy[jx].center().set( T(0.) ), true );
      }
      _e2x( _nx, _MVYXPd, _MVYXPy, false );
      _pDAG->eval( _opADJRHS[ifct], _PMADJRHS, _nx, _pADJRHS+ifct*_nx, _MVYXPg,
                   _nVAR-_npar, _pVAR, _MVYXPVAR );
      _RHS_PM_ELL2( _nx, _PMenv, _PMydot, _MVYXPg, _npar, _Ir, _Ax, _Idy, _Ay, _Idydot );
    }

    // Construct the ellipsoidal remainder derivatives
    if( options.QSCALE )
      _RHS_ELL( _nx, _Qy, _Ay, _Ax, _Idydot, _Qydot, options.QTOL, machprec(), neg, _PMy );
    else
      _RHS_ELL( _nx, _Qy, _Ay, _Ax, _Idydot, _Qydot, options.QTOL, machprec(), neg, (PVT*)0 );

    // Whether or not to ignore the adjoint remainder
    if( !options.PMNOREM )
      _PME2vec( _PMenv, _nx, _PMydot, _Qydot, ydot );
    else
      _PME2vec( _PMenv, _nx, _PMydot, 0, ydot );
#ifdef MC__ODEBNDS_BASE_DINEQPM_DEBUG
    _print_interm( t, _nx, _PMydot, E(_nx,_Qydot), "PMydot Intermediate", std::cerr );
    { std::cout << "--paused--"; int dum; std::cin >> dum; }
#endif
    break;
  }

  return true;
}

template <typename T, typename PMT, typename PVT>
inline void
ODEBNDS_BASE<T,PMT,PVT>::_RHS_PM_ELL0
( const unsigned nx, PVT*PMg, T*Idgdy, const T*Idx, double*Ax,
  const T*Idy, double*Ay, T*Idydot )
{
  for( unsigned ix=0; ix<nx; ix++ ){
    for( unsigned jx=0; jx<nx; jx++ ){
      Ay[ix+jx*nx] = Op<T>::mid( Idgdy[ix*2*nx+jx] );
      Ax[ix+jx*nx] = Op<T>::mid( Idgdy[ix*2*nx+nx+jx] );
    }
    Idydot[ix] = PMg[ix].remainder();
    for( unsigned jx=0; jx<nx; jx++ )
      Idydot[ix] += ( Idgdy[ix*2*nx+nx+jx] - Ax[ix+jx*nx] ) * Idx[jx]
                  + ( Idgdy[ix*2*nx+jx] - Ay[ix+jx*nx] ) * Idy[jx];
#ifdef MC__ODEBNDS_BASE_DINEQPM_DEBUG
    std::cout << "Idydot[" << ix << "] = " << Idydot[ix]
              << " : " << Op<T>::mid(Idydot[ix]) << std::endl;
#endif
  }
}

template <typename T, typename PMT, typename PVT>
inline void
ODEBNDS_BASE<T,PMT,PVT>::_RHS_PM_ELL1
( const unsigned nx, PVT*PMg, PVT*MVYXPdgdy, const T*Idx, double*Ax,
  const T*Idy, double*Ay, T*Idydot )
{
  for( unsigned ix=0; ix<nx; ix++ ){
#ifdef MC__ODEBNDS_BASE_DINEQPM_DEBUG
    for( unsigned jx=0; jx<2*nx; jx++ )
      std::cout << "MVYXPdgdy[" << ix << "," << jx << "] = " << MVYXPdgdy[ix*nx+jx] << std::endl;
#endif
    // Extract constant coefficients and set to 0
    for( unsigned jx=0; jx<nx; jx++ ){
      Ay[ix+jx*nx] = ( MVYXPdgdy[ix*2*nx+jx].center() ).constant( true );
      Ax[ix+jx*nx] = ( MVYXPdgdy[ix*2*nx+nx+jx].center() ).constant( true );
    }
#ifdef MC__ODEBNDS_BASE_DINEQPM_DEBUG
    std::cout << "Ax[" << ix << ",#] = ";
    for( unsigned jx=0; jx<nx; jx++ )
      std::cout << Ax[ix+jx*nx] << "  ";
    std::cout << "Ay[" << ix << ",#] = ";
    for( unsigned jx=0; jx<nx; jx++ )
      std::cout << Ay[ix+jx*nx] << "  ";
    std::cout << std::endl;
    for( unsigned jx=0; jx<2*nx; jx++ )
      std::cout << "MVYXPdgdy[" << ix << "," << jx << "] = " << MVYXPdgdy[ix*nx+jx] << std::endl;
#endif
    // Bound remaining terms
    Idydot[ix] = PMg[ix].remainder();
    for( unsigned jx=0; jx<nx; jx++ )
      Idydot[ix] += MVYXPdgdy[ix*2*nx+nx+jx].bound() * Idx[jx]
                  + MVYXPdgdy[ix*2*nx+jx].bound() * Idy[jx];
#ifdef MC__ODEBND_BASES_DINEQPM_DEBUG
    std::cout << "Idydot[" << ix << "] = " << Izdot[ix]
              << " : " << Op<T>::mid(Idydot[ix]) << std::endl;
#endif
  }
}

template <typename T, typename PMT, typename PVT>
inline void
ODEBNDS_BASE<T,PMT,PVT>::_RHS_PM_ELL2
( const unsigned nx, PMT*PMenv, PVT*PMg, PVT*MVYXPg, const unsigned np,
  const T*Idx, double*Ax, const T*Idy, double*Ay, T*Idydot )
{
  for( unsigned ix=0; ix<nx; ix++ ){
#ifdef MC__ODEBNDS_BASE_DINEQPM_DEBUG
    std::cout << "MVYXPg[" << ix << "] = " << MVYXPg[ix] << std::endl;
#endif
    // Extract polynomial model in P and set to 0
    MVYXPg[ix].get( PMg[ix].set(PMenv), true );
#ifdef MC__ODEBNDS_BASE_DINEQPM_DEBUG
    std::cout << "PMydot[" << ix << "] =" << PMg[ix] << std::endl;
    std::cout << "MVYXPg[" << ix << "] = " << MVYXPg[ix] << std::endl;
    { int dum; std::cin >> dum; }
#endif
    // Extract linear part in X and set to 0
    for( unsigned jx=0; jx<nx; jx++ ){
      Ax[ix+jx*nx] = MVYXPg[ix].linear( np+jx, true );
      Ay[ix+jx*nx] = MVYXPg[ix].linear( np+nx+jx, true );
    }
#ifdef MC__ODEBNDS_BASE_DINEQPM_DEBUG
    std::cout << "Ax[" << ix << ",#] = ";
    for( unsigned jx=0; jx<nx; jx++ )
      std::cout << Ax[ix+jx*nx] << "  ";
    std::cout << "Ay[" << ix << ",#] = ";
    for( unsigned jx=0; jx<nx; jx++ )
      std::cout << Ay[ix+jx*nx] << "  ";
    std::cout << std::endl;
    std::cout << "MVYXPg[" << ix << "] = " << MVYXPg[ix] << std::endl;
#endif
    // Bound remaining terms
    Idydot[ix] = MVYXPg[ix].bound();
#ifdef MC__ODEBNDS_BASE_DINEQPM_DEBUG
    std::cout << "Idydot[" << ix << "] = " << Idydot[ix]
              << " : " << Op<T>::mid(Idydot[ix]) << std::endl;
#endif
  }
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
    _vec2PMI( x, _PMenv, _nx, _PMx, true );// set current state bounds
    _vec2PMI( y, _PMenv, _nx, _PMy, true );// set current adjoint bounds
    break;

  case OPT::DINEQ:
    if( !bndinit ) break;
    _vec2PMI( x, _PMenv, _nx, _PMx );// set current state bounds
    _vec2PMI( y, _PMenv, _nx, _PMy );// set current adjoint bounds
    break;

  case OPT::ELLIPS:
  default:
    if( !bndinit ) break;
    _vec2PME( x, _PMenv, _nx, _PMx, _Q, _Er, _Ir );// set current state bounds
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

















