// Copyright (C) 2015 Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

#ifndef MC__ODEBND_BASE_HPP
#define MC__ODEBND_BASE_HPP

#undef  MC__ODEBND_BASE_DINEQI_DEBUG
#undef  MC__ODEBND_BASE_DINEQPM_DEBUG
#undef  MC__ODEBND_BASE_MVXP_USE

#include <stdexcept>
#include <cassert>
#include <cmath>
#include <fstream>
#include <vector>
#include <sys/time.h>

#include "ellipsoid.hpp"

// *** TO DO
// - Detect block structure and use it in ellipsoidal approach
// - Implement scaling in ellipsoidal approach

namespace mc
{
//! @brief C++ base class for computing enclosures of the reachable set of parametric ODEs using continuous-time set-valued integration.
////////////////////////////////////////////////////////////////////////
//! mc::ODEBND_BASE is a C++ base class for computing enclosures of the
//! reachable set of parametric ordinary differential equations
//! (ODEs) using continuous-time set-valued integration.
////////////////////////////////////////////////////////////////////////
template <typename T, typename PMT, typename PVT>
class ODEBND_BASE:
  public virtual BASE_DE
{
 typedef Ellipsoid E;

 public:
  /** @defgroup ODEBND_BASE Continuous-time set-valued integration of parametric ODEs
   *  @{
   */
  //! @brief Default constructor
  ODEBND_BASE();

  //! @brief Virtual destructor
  virtual ~ODEBND_BASE();

 protected:
  //! @brief list of operations in RHS evaluation
  std::list<const FFOp*> _opRHS;

  //! @brief array of list of operations in individual RHS evaluations
  std::list<const FFOp*> *_opRHSi;

  //! @brief list of operations in RHS Jacobian
  std::list<const FFOp*> _opJAC;

  //! @brief list of operations in quadrature evaluation
  std::list<const FFOp*> _opQUAD;

  //! @brief list of operations in IC evalution
  std::list<const FFOp*> _opIC;

  //! @brief list of operations in IC Jacobian
  std::list<const FFOp*> _opDIC;

  //! @brief preallocated array for evaluation of RHS function in T arithmetic
  T* _IRHS;

  //! @brief preallocated array for evaluation of RHS Jacobian in T arithmetic
  T* _IJAC;

  //! @brief preallocated array for evaluation of IC function in T arithmetic
  T* _IIC;

  //! @brief preallocated array for evaluation of IC Jacobian in T arithmetic
  T* _IDIC;

  //! @brief preallocated array for evaluation of RHS function in PM arithmetic
  PVT* _PMRHS;

  //! @brief preallocated array for evaluation of RHS Jacobian in PM arithmetic
  PVT* _PMJAC;

  //! @brief preallocated array for evaluation of IC function in PM arithmetic
  PVT* _PMIC;

  //! @brief preallocated array for evaluation of IC Jacobian in PM arithmetic
  PVT* _PMDIC;

  //! @brief const pointer to RHS function in current stage of ODE system
  const FFVar* _pRHS;

  //! @brief const pointer to RHS Jacobian in current stage of ODE system
  const FFVar* _pJAC;

  //! @brief const pointer to quadrature integrand in current stage of ODE system
  const FFVar* _pQUAD;

  //! @brief const pointer to IC function in current stage of ODE system
  const FFVar* _pIC;

  //! @brief const pointer to IC Jacobian in current stage of ODE system
  const FFVar* _pDIC;

  //! @brief number of effective parameters (possibly larger than _np due to lifting)
  unsigned _npar;

  //! @brief number of variables for DAG evaluation
  unsigned _nVAR;

  //! @brief pointer to variables for DAG evaluation
  FFVar* _pVAR;

  //! @brief pointer of variable bounds for DAG evaluation
  T *_IVAR;

  //! @brief pointer to adjoint time **DO NOT FREE**
  T *_It;

  //! @brief pointer to state interval bounds **DO NOT FREE**
  T *_Ix;

  //! @brief pointer to parameter interval bounds **DO NOT FREE**
  T *_Ip;

  //! @brief state interval bounds time derivatives
  T *_Ixdot;

  //! @brief state lower bound time derivatives
  double *_xLdot;

  //! @brief state upper bound time derivatives
  double *_xUdot;

  //! @brief quadrature interval bounds **DO NOT FREE**
  T *_Iq;

  //! @brief quadrature derivative interval bounds
  T *_Iqdot;

  //! @brief polynomial model environment
  PMT *_PMenv;

  //! @brief variable polynomial models for DAG evaluation
  PVT *_PMVAR;

  //! @brief time polynomial model **DO NOT FREE**
  PVT *_PMt;

  //! @brief state polynomial model **DO NOT FREE**
  PVT *_PMx;

  //! @brief parameter polynomial model **DO NOT FREE**
  PVT *_PMp;

  //! @brief state derivative polynomial model
  PVT *_PMxdot;

  //! @brief state PM remainder lower bound time derivatives
  double *_RxLdot;

  //! @brief state PM remainder upper bound time derivatives
  double *_RxUdot;

  //! @brief quadrature polynomial model **DO NOT FREE**
  PVT *_PMq;

  //! @brief quadrature derivative polynomial model
  PVT *_PMqdot;

  //! @brief quadrature PM remainder radius time derivatives
  double *_radRqdot;

  //! @brief parameter reference
  double *_pref;

  //! @brief state reference
  double *_xref;

  //! @brief linear transformation A matrix
  double *_A;

  //! @brief linear transformation B matrix
  double *_B;

  //! @brief linear transformed state interval bounds
  T *_Ir;

  //! @brief RHS Jacobian interval bounds
  T *_Idfdx;

  //! @brief linear transformed state ellipsoidal bounds
  E _Er;

  //! @brief shape matrix (lower triangular) in ellipsoidal bounds
  double *_Q;

  //! @brief state reference time derivatives
  double *_xrefdot;

  //! @brief linear transformation B matrix time derivatives
  double *_Bdot;

  //! @brief rotated state interval bounds time derivatives
  T *_Irdot;

  //! @brief shape matrix time derivatives in ellipsoidal bounds
  double *_Qdot;

  //! @brief polynomial model environment for mean-value theorem in (X,P)
  PMT *_MVXPenv;

  //! @brief rotated state polynomial model
  PVT *_MVXPd;

  //! @brief RHS polynomial model
  PVT *_MVXPf;

  //! @brief RHS Jacobian polynomial model
  PVT *_MVXPdfdx;

  //! @brief polynomial models for variables in mean-value theorem
  PVT *_MVXPVAR;

  //! @brief pointer to time polynomial model **DO NOT FREE**
  PVT *_MVXPt;

  //! @brief pointer to state polynomial models **DO NOT FREE**
  PVT *_MVXPx;

  //! @brief pointer to parameter polynomial models **DO NOT FREE**
  PVT *_MVXPp;

  //! @brief polynomial model environment for linear transformation of intial value
  PMT *_MVPenv;

  //! @brief array of parameter polynomial model for initial value
  PVT *_MVPp;

  //! @brief array of initial state polynomial model for initial value
  PVT *_MVPx;

  //! @brief Function to initialize state interval bounding
  template <typename OPT> bool _INI_I_STA
    ( const OPT&options, const unsigned np, const T*Ip, const unsigned ns );

  //! @brief Static function converting interval bounds to GSL array
  template <typename REALTYPE> static void _I2vec
    ( const unsigned nx, const double*xL, const double*xU, REALTYPE*vec );

  //! @brief Static function converting interval bounds to GSL array
  template <typename REALTYPE> static void _I2vec
    ( const unsigned nx, const T*Ix, REALTYPE*vec, const bool centered=false );

  //! @brief Static function converting ellipsoidal bounds to GSL array
  template <typename REALTYPE> static void _E2vec
    ( const unsigned nx, const unsigned np, const double*xref, const double*Qx,
      const double*Bx, REALTYPE*vec );

  //! @brief Function converting GSL array to interval bounds
  template <typename REALTYPE> static void _vec2I
    ( const REALTYPE*vec, const unsigned nx, T*Ix, const bool centered=false );

  //! @brief Static function converting GSL array to ellipsoidal bounds
  template <typename REALTYPE> static void _vec2E
    ( const REALTYPE*vec, const unsigned nx, const unsigned np, double*Q,
      E&Ed, T*Id, const double*pref, const T*Ip, double*B, double*xref,
      T*Ix );

  //! @brief Static function converting rotated states into original coordinates
  template<typename U> static void _ep2x
    ( const unsigned nx, const unsigned np, const U*d, const double*pref,
      const U*p, const double*B, const double*xref, U*x );

  //! @brief Static function converting ellipsoids into interval bounds
  static void _ep2x
    ( const unsigned nx, const unsigned np, const double*Q, E&Ed, T*Id,
      const double*pref, const T*Ip, const double*B, const double*xref, T*Ix );

  //! @brief Function to initialize quarature interval bounds
  template <typename REALTYPE, typename OPT> bool _IC_I_QUAD
    ( const OPT&options, REALTYPE*vec );

  //! @brief Function to initialize state interval bounds
  template <typename REALTYPE, typename OPT> bool _IC_I_STA
    ( const OPT&options, REALTYPE*vec );

  //! @brief Function to initialize state interval bounds w/ ellipsoidal bounds 
  static void _IC_I_ELL
    ( FFGraph*DAG, std::list<const FFOp*>&opIC, const unsigned nx, const FFVar*pIC,
      const unsigned nVAR, const FFVar*pVAR, PVT*MVPVAR, PVT*MVPx, double*xref,
      double*Q, const unsigned np, double*B, T*Ix );

  //! @brief Function to reinitialize state interval bounds
  template <typename REALTYPE, typename OPT> bool _CC_I_STA
    ( const OPT&options, const unsigned iIC, const double t, REALTYPE*vec );

  //! @brief Function to reinitialize state interval bounds w/ ellipsoidal bounds 
  static void _CC_I_ELL
    ( FFGraph*DAG, std::list<const FFOp*>&opf, PVT*PMf,
      const unsigned nx, const FFVar*pf, const unsigned nVAR,
      const FFVar*pVAR, PVT*MVXPVAR, PVT*MVXPf, const E&Ed, double*Af,
      const unsigned np, double*reff, double*Bf, T*Idf, double*Qf,
      const double QTOL, const double EPS );

  //! @brief Function to set RHS and QUAD pointers
  template <typename OPT> bool _SET_I_STA
    ( const OPT&options, const unsigned iRHS, const unsigned iQUAD );

  //! @brief Function to calculate the RHS of auxiliary ODEs in interval arithmetic
  template <typename REALTYPE, typename OPT> bool _RHS_I_STA
    ( const OPT&options, double t, const REALTYPE*x, REALTYPE*xdot );

  //! @brief Static function to calculate the RHS of auxiliary ODEs in interval arithmetic w/o contractor
  static void _RHS_I_NONE
    ( FFGraph*DAG, std::list<const FFOp*>&opRHS, T*IRHS,
      const unsigned nx, const FFVar*pRHS, const unsigned nVAR,
      const FFVar*pVAR, T*IVAR, T*Ixdot );

  //! @brief Static function to calculate the RHS of auxiliary ODEs in interval arithmetic w/ differential inequality contractor
  static void _RHS_I_DI
    ( FFGraph*DAG, std::list<const FFOp*>*opRHSi, T*IRHS,
      const unsigned nx, const FFVar*pRHS, const unsigned nVAR,
      const FFVar*pVAR, T*IVAR, T*Ixdot, double*xLdot, double*xUdot );

  //! @brief Static function to calculate the RHS of auxiliary ODEs in interval arithmetic w/ ellipsoidal contractor
  static void _RHS_I_ELL
    ( FFGraph*DAG, std::list<const FFOp*>&opRHS, PVT*PMRHS,
      const unsigned nx, const FFVar*pRHS, const unsigned nVAR,
      const FFVar*pVAR, PVT*MVXPVAR, PVT*MVXPf, const double*Q, double*A,
      const unsigned np, double*xrefdot, double*Bdot, T*Iddot, double*Qdot,
      const double QTOL, const double EPS );

  //! @brief Function to calculate the RHS of auxiliary ODEs in interval arithmetic
  template <typename REALTYPE, typename OPT> bool _RHS_I_QUAD
    ( const OPT&options, double t, const REALTYPE*x, REALTYPE*qdot,
      const bool bndinit=true );

  //! @brief Static function to calculate the quadratures in interval arithmetic
  static void _QUAD_I
    ( FFGraph*DAG, std::list<const FFOp*>&opQUAD, T*IQUAD,
      const unsigned nq, const FFVar*pQUAD, const unsigned nVAR,
      const FFVar*pVAR, T*IVAR, T*Iqdot );

  //! @brief Function to calculate the Jacobian of auxiliary ODEs in interval arithmetic
  template <typename REALTYPE, typename OPT> bool _JAC_I_STA
    ( const OPT&options, double t, const REALTYPE*x, REALTYPE*jac, REALTYPE*xdot );

  //! @brief Function to calculate the functions at intermediate/end point
  bool _FCT_I_STA
    ( const unsigned iFCT, const double t, T*If );

  //! @brief Function converting polynomial model with interval remainder to GSL array
  template <typename REALTYPE> static void _PMI2vec
    ( const PMT*PMenv, const unsigned nx, const PVT*PMx, const double*RxL,
      const double*RxU, REALTYPE*vec );

  //! @brief Function converting polynomial model with interval remainder to GSL array
  template <typename REALTYPE> static void _PMI2vec
    ( const PMT*PMenv, const unsigned nx, const PVT*PMx, const T*IRx,
      REALTYPE*vec );

  //! @brief Function converting polynomial model with interval remainder to GSL array
  template <typename REALTYPE> static void _PMI2vec
    ( const PMT*PMenv, const unsigned nx, const PVT*PMx, REALTYPE*vec,
      const bool centered=false );

  //! @brief Function converting polynomial model with ellipsoidal remainder to GSL array
  template <typename REALTYPE> static void _PME2vec
    ( const PMT*PMenv, const unsigned nx, const PVT*PMx, const double*Qx,
      REALTYPE*vec );

  //! @brief Function converting GSL array to polynomial model with interval remainder
  template <typename REALTYPE> static void _vec2PMI
    ( const REALTYPE*vec, PMT*PMenv, const unsigned nx, PVT*PMx,
      const bool centered=false );

  //! @brief Function converting GSL array to polynomial model with ellipsoidal remainder
  template <typename REALTYPE> static void _vec2PME
    ( const REALTYPE*vec, PMT*PMenv, const unsigned nx, PVT*PMx, double*Qr,
      E&Er, T*Ir );

  //! @brief Function converting rotated interval remainder bound back into original coordinates
  template <typename U> static void _e2x
    ( const unsigned nx, const U*d, U*x, const bool reinit=true );

  //! @brief Function converting ellipsoidal remainder bound into interval remainder bound
  static void _e2x
    ( const unsigned nx, const double*Qr, E&Er, T*Ir );

  //! @brief Function to initialize GSL for state polynomial models
  template <typename OPT> bool _INI_PM_STA
    ( const OPT&options, const unsigned np, const PVT*PMp, const unsigned ns );

  //! @brief Function to initialize quarature polynomial models
  template <typename REALTYPE, typename OPT> bool _IC_PM_QUAD
    ( const OPT&options, REALTYPE*vec );

  //! @brief Function to initialize state polynomial models
  template <typename REALTYPE, typename OPT> bool _IC_PM_STA
    ( const OPT&options, REALTYPE*vec );

  //! @brief Function to initialize state polynomial model w/ ellipsoidal remainder
  static void _IC_PM_ELL
    ( FFGraph*DAG, std::list<const FFOp*>&opIC, const unsigned nx, const FFVar*pIC,
      const unsigned np, const FFVar*pP, PVT*PMp, PVT*PMx, double*Qr, E&Er, T*Ir );

  //! @brief Function to reinitialize state polynomial bounds
  template <typename REALTYPE, typename OPT> bool _CC_PM_STA
    ( const OPT&options, const unsigned iIC, const double t, REALTYPE*vec );

  //! @brief Function to reinitialize state polynomial model w/ ellipsoidal remainder
  static void _CC_PM_ELL
    ( const unsigned nx, const E&Exr, const double*Afr, const T*Ifr,
      double*Qfr, const double QTOL, const double EPS );

  //! @brief Function to set RHS and QUAD pointers
  template <typename OPT> bool _SET_PM_STA
    ( const OPT&options, const unsigned iRHS, const unsigned iQUAD );

  //! @brief Function to calculate the RHS of auxiliary ODEs in polynomial mode arithmetic
  template <typename REALTYPE, typename OPT> bool _RHS_PM_STA
    ( const OPT&options, double t, const REALTYPE*x, REALTYPE*xdot );

  //! @brief Static function to calculate the RHS of auxiliary ODEs in polynomial model arithmetic w/o contractor
  static void _RHS_PM_NONE
    ( FFGraph*DAG, std::list<const FFOp*>&opRHS, PVT*PMRHS,
      const unsigned nx, const FFVar*pRHS, const unsigned nVAR,
      const FFVar*pVAR, PVT*PMVAR, PVT*PMxdot );

  //! @brief Static function to calculate the RHS of auxiliary ODEs in polynomial model arithmetic w/ differential inequality contractor
  static void _RHS_PM_DI
    ( FFGraph*DAG, std::list<const FFOp*>*opRHSi, PVT*PMRHS,
      const unsigned nx, const FFVar*pRHS, const unsigned nVAR,
      const FFVar*pVAR, PVT*PMVAR, PVT*PMxdot, double*rLdot, double*rUdot );

  //! @brief Static function to calculate the RHS of auxiliary ODEs in polynomial model arithmetic w/ ellipsoidal contractor - approximation using mean-value theorem and interval analysis
  static void _RHS_PM_ELL0
    ( FFGraph*DAG, std::list<const FFOp*>&opRHS, PVT*PMRHS,
      std::list<const FFOp*>&opJAC, T*IJAC, const unsigned nx,
      const FFVar*pRHS, const FFVar*pJAC, const unsigned nVAR,
      const FFVar*pVAR, PVT*PMVAR, T*IVAR, PVT*PMf, T*Idfdx,
      const T*Ir, double*Ar, T*Irdot );

  //! @brief Static function to calculate the RHS of auxiliary ODEs in polynomial model arithmetic w/ ellipsoidal contractor - approximation using mean-value theorem and PM arithmetic
  static void _RHS_PM_ELL1
    ( FFGraph*DAG, std::list<const FFOp*>&opRHS, PVT*PMRHS,
      std::list<const FFOp*>&opJAC, PVT*PMJAC, const unsigned nx,
      const FFVar*pRHS, const FFVar*pJAC, const unsigned nVAR,
      const FFVar*pVAR, PVT*PMVAR, PVT*MVXPVAR, PVT*PMf, PVT*MVXPdfdx,
      const T*Ir, double*Ar, T*Irdot );

  //! @brief Static function to calculate the RHS of auxiliary ODEs in polynomial model arithmetic w/ ellipsoidal contractor - joint polynomial model in states and parameters
  static void _RHS_PM_ELL2
    ( FFGraph*DAG, std::list<const FFOp*>&opRHS, PVT*PMRHS,
      const unsigned nx, const FFVar*pRHS, const unsigned nVAR, const FFVar*pVAR,
      PVT*MVXPVAR, PMT*PMenv, PVT*PMf, PVT*MVXPf, const unsigned np, const T*Ir,
      double*Ar, T*Irdot );

  //! @brief Static function to calculate the RHS of auxiliary ODEs in polynomial model arithmetic - ellipsoidal contractor
  static void _RHS_PM_ELL
    ( const unsigned nx, const double*Qr, const double*Ar, const T*Irdot,
      double*Qrdot, const double QTOL, const double EPS );

  //! @brief Function to calculate the RHS of auxiliary ODEs in polynomial mode arithmetic
  template <typename REALTYPE, typename OPT> bool _RHS_PM_QUAD
    ( const OPT&options, double t, const REALTYPE*x, REALTYPE*qdot,
      const bool bndinit=true );

  //! @brief Static function to calculate the quadratures in polynomial model arithmetic
  static void _QUAD_PM
    ( FFGraph*DAG, std::list<const FFOp*>&opQUAD, PVT*PMQUAD,
      const unsigned nq, const FFVar*pQUAD, const unsigned nVAR,
      const FFVar*pVAR, PVT*PMVAR, PVT*PMqdot );

  //! @brief Function to calculate the Jacobian of auxiliary ODEs in polynomial mode arithmetic
  template <typename REALTYPE, typename OPT> bool _JAC_PM_STA
    ( const OPT&options, double t, const REALTYPE*x, REALTYPE*jac, REALTYPE*xdot );

  //! @brief Function to calculate the functions at intermediate/end point
  bool _FCT_PM_STA
    ( const unsigned iFCT, const double t, PVT*PMf );

  //! @brief Computes Hausdorff distance between interval enclosure and actual reachable set of parametric ODEs, using parameter sampling
  template <typename ODEBND, typename ODESLV> inline bool _hausdorff
    ( const unsigned ns, const double*tk, const T*Ip, double**Hxk,
      ODEBND&dineq, ODESLV&traj, const unsigned nsamp,
      std::ostream&os=std::cout );

  //! @brief Computes Hausdorff distance between polynomial model remainder enclosure and actual remainder function range, using parameter sampling
  template <typename ODEBND, typename ODESLV> inline bool _hausdorff
    ( const unsigned ns, const double*tk, const PVT*PMp, double**Hxk,
      ODEBND&dineq, ODESLV&traj, const unsigned nsamp,
      std::ostream&os=std::cout );

  //! @brief Function to bound the remainder function relative to a polynomial model at sampling points -- return value is status
  template<typename ODESLV> bool _remainders
    ( const unsigned ns, const double*tk, const unsigned np, const T*Ip,
      const PVT*const*PMxk, T**Rxk, ODESLV&traj, const unsigned nsamp,
      std::ostream&os=std::cout );

  //! @brief Recrusive function computing bounds on errors between solutions of IVP in ODEs and polynomial approximant using sampling
  template<typename ODESLV> bool _remainders
    ( const unsigned ns, const double*tk, const unsigned np, const T*Ip,
      const PVT*const*PMxk, T**Rxk, ODESLV&traj, const unsigned nsamp,
      unsigned* vsamp, const unsigned ip, double*p, double**xk,
      std::ostream&os );

  //! @brief Function to display intermediate results
  template<typename U> static void _print_interm
    ( const unsigned nx, const U*x, const std::string&var, std::ostream&os=std::cout );

  //! @brief Function to display intermediate results
  template<typename U> static void _print_interm
    ( const double t, const unsigned nx, const U*x, const std::string&var,
      std::ostream&os=std::cout );

  //! @brief Function to display intermediate results
  template<typename U, typename V> static void _print_interm
    ( const double t, const unsigned nx, const U*x, const V&r, const std::string&var,
      std::ostream&os=std::cout );

  //! @brief Position in symmetric matrix stored in lower triangular form
  static unsigned _ndxLT
    ( const unsigned i, const unsigned j, const unsigned n )
    { return( i<j? _ndxLT(j,i,n): i+j*n-j*(j+1)/2 ); }

  //! @brief Function computing set diameter (max-norm component-wise)
  static double _diam
    ( const unsigned nx, T*X );

  //! @brief Function computing set diameter (max-norm component-wise)
  static double _diam
    ( const unsigned nx, PVT*X );

  //! @brief Function computing Hausdorff distance between intervals
  template <typename U> static double _dH
    ( const U&X, const U&Y );

  //! @brief Record results in file <a>bndrec</a>, with accuracy of <a>iprec</a> digits
  template <typename VRES> static void _record
    ( std::ofstream&ofile, const VRES&bnd, const unsigned iprec=5 );

  //! @brief Private methods to block default compiler methods
  ODEBND_BASE(const ODEBND_BASE&);
  ODEBND_BASE& operator=(const ODEBND_BASE&);
};

template <typename T, typename PMT, typename PVT>
inline
ODEBND_BASE<T,PMT,PVT>::ODEBND_BASE
()
: BASE_DE(), _pRHS(0), _pJAC(0), _pQUAD(0), _pIC(0), _pDIC(0), _npar(0),
  _nVAR(0), _pVAR(0)
{
  // Initalize state/parameter arrays
  _opRHSi = 0;
  _IRHS = _IJAC = _IIC = _IDIC = _IVAR = _Ip = _Ix = _Ixdot = _Iq = _Iqdot = 0;
  _PMRHS = _PMJAC = _PMIC = _PMDIC = _PMVAR = _PMp = _PMx = _PMxdot = _PMq = _PMqdot = 0;
  _xLdot = _xUdot = _RxLdot = _RxUdot = _radRqdot = 0;
  _PMenv = 0;

  // Initialize parameterization arrays
  _pref = _xref = _xrefdot = 0;
  _A = _B = _Bdot = _Q = _Qdot = 0;
  _Ir = _Idfdx = _Irdot = 0;

  // Initialize polynomial model environments
  _MVXPenv = _MVPenv = 0;
  _MVXPd = _MVXPf = _MVXPdfdx = _MVXPVAR = _MVXPx = _MVXPp = _MVPp = _MVPx = 0;

  // Initialize ellipsoidal calculus
  E::options.PSDCHK = false;
}

template <typename T, typename PMT, typename PVT>
inline
ODEBND_BASE<T,PMT,PVT>::~ODEBND_BASE
()
{
  /* DO NOT FREE _pRHS */
  delete[] _opRHSi;
  delete[] _IRHS;
  delete[] _PMRHS;
  delete[] _IJAC;
  delete[] _PMJAC;
  delete[] _pJAC;
  delete[] _IIC;
  delete[] _PMIC;
  delete[] _IDIC;
  delete[] _PMDIC;
  delete[] _pDIC;
  delete[] _pVAR;

  // Free state/quadrature arrays -- Do *NOT* delete _Ip _PMp _PMenv
  delete[] _IVAR;   // **DO NOT DELETE _Ix, _Ip, _Iq**
  delete[] _PMVAR;  // **DO NOT DELETE _PMx, _PMp, _PMq**
  delete[] _PMxdot;
  delete[] _RxLdot;
  delete[] _RxUdot;
  delete[] _PMqdot;
  delete[] _radRqdot;
  delete[] _Ixdot;
  delete[] _xLdot;
  delete[] _xUdot;
  delete[] _Iqdot;

  // Free linear transformation arrays
  delete[] _pref;
  delete[] _xref;
  delete[] _A;
  delete[] _B;
  delete[] _Ir;
  delete[] _Idfdx;
  delete[] _Q;
  delete[] _xrefdot;
  delete[] _Bdot;
  delete[] _Irdot;
  delete[] _Qdot;
  delete[] _MVXPd;
  delete[] _MVXPf;
  delete[] _MVXPdfdx;
  delete[] _MVXPVAR;  // **DO NOT DELETE _MVXPx, _MVXPp**
  delete   _MVXPenv;
  delete[] _MVPp;
  delete[] _MVPx;
  delete   _MVPenv;
}

template <typename T, typename PMT, typename PVT>
template <typename REALTYPE>
inline void
ODEBND_BASE<T,PMT,PVT>::_vec2I
( const REALTYPE*vec, const unsigned nx, T*Ix, const bool centered )
{
  unsigned ivec = 0;
  for( unsigned ix=0; ix<nx; ix++ ){
    if( !centered ){
      Ix[ix] = T( vec[ivec], vec[ivec+1] );
      ivec += 2;
    }
    else{
      Ix[ix] = T( -vec[ivec], vec[ivec] );
      ivec++;
    }
  }
  return;
}

template <typename T, typename PMT, typename PVT>
template <typename U> inline void
ODEBND_BASE<T,PMT,PVT>::_ep2x
( const unsigned nx, const unsigned np, const U*d, const double*pref,
  const U*p, const double*B, const double*xref, U*x )
{
  for( unsigned ix=0; ix<nx; ix++ ){   
    x[ix] = xref[ix] + d[ix];
    for( unsigned jp=0; jp<np; jp++ )
      x[ix] += ( p[jp] - pref[jp] ) * B[jp*nx+ix];
  }
  return;
}

template <typename T, typename PMT, typename PVT>
inline void
ODEBND_BASE<T,PMT,PVT>::_ep2x
( const unsigned nx, const unsigned np, const double*Q, E&Ed, T*Id,
  const double*pref, const T*Ip, const double*B, const double*xref, T*Ix )
{
  Ed.set( nx, Q );
#ifdef MC__ODEBND_BASE_DINEQI_DEBUG
  std::cout << "Ed =" << Ed << std::endl;
#endif
  for( unsigned ix=0; ix<nx; ix++ )
    Id[ix] = T( Ed.l(ix), Ed.u(ix) );
  _ep2x( nx, np, Id, pref, Ip, B, xref, Ix );

#ifdef MC__ODEBND_BASE_DINEQI_DEBUG
  for( unsigned ix=0; ix<nx; ix++ )
    std::cout << "Id[" << ix << "] = " << Id[ix] << std::endl;
  for( unsigned ix=0; ix<nx; ix++ )
    std::cout << "Ix[" << ix << "] = " << Ix[ix] << std::endl;
  { int dum; std::cin >> dum; }
#endif
  return;
}

template <typename T, typename PMT, typename PVT>
template <typename REALTYPE>
inline void
ODEBND_BASE<T,PMT,PVT>::_vec2E
( const REALTYPE*vec, const unsigned nx, const unsigned np, double*Q,
  E&Ed, T*Id, const double*pref, const T*Ip, double*B, double*xref,
  T*Ix )
{
  unsigned ivec = 0;
  for( unsigned ix=0; ix<nx; ix++ )
    xref[ix] = vec[ivec++];
  for( unsigned iQ=0; iQ<nx*(nx+1)/2; iQ++ )
    Q[iQ] = vec[ivec++];
  for( unsigned iB=0; iB<nx*np; iB++ )
    B[iB] = vec[ivec++];
  return _ep2x( nx, np, Q, Ed, Id, pref, Ip, B, xref, Ix );
}

template <typename T, typename PMT, typename PVT>
template <typename REALTYPE>
inline void
ODEBND_BASE<T,PMT,PVT>::_I2vec
( const unsigned nx, const double*xL, const double*xU,
  REALTYPE*vec )
{
  unsigned ivec = 0;
  for( unsigned ix=0; ix<nx; ix++ ){
    vec[ivec++] = xL[ix];
    vec[ivec++] = xU[ix];
  }
  return;
}

template <typename T, typename PMT, typename PVT>
template <typename REALTYPE>
inline void
ODEBND_BASE<T,PMT,PVT>::_I2vec
( const unsigned nx, const T*Ix, REALTYPE*vec, const bool centered )
{
  unsigned ivec = 0;
  if( !centered ){
    for( unsigned ix=0; ix<nx; ix++ ){
      vec[ivec++] = Op<T>::l( Ix[ix] );
      vec[ivec++] = Op<T>::u( Ix[ix] );
    }
  }
  else
    for( unsigned ix=0; ix<nx; ix++ )
      vec[ivec++] = 0.5*Op<T>::diam( Ix[ix] );
  return;
}

template <typename T, typename PMT, typename PVT>
template <typename REALTYPE>
inline void
ODEBND_BASE<T,PMT,PVT>::_E2vec
( const unsigned nx, const unsigned np, const double*xref, const double*Q,
  const double*B, REALTYPE*vec )
{
  unsigned ivec = 0;
  for( unsigned ix=0; ix<nx; ix++ )
    vec[ivec++] = xref[ix];
  for( unsigned iQ=0; iQ<nx*(nx+1)/2; iQ++ )
    vec[ivec++] = Q[iQ];
  for( unsigned iB=0; iB<nx*np; iB++ )
    vec[ivec++] = B[iB];
  return;
}

template <typename T, typename PMT, typename PVT>
template <typename REALTYPE, typename OPT> inline bool
ODEBND_BASE<T,PMT,PVT>::_IC_I_QUAD
( const OPT&options, REALTYPE*vec )
{
  if( !_vQUAD.size() || !_nq ) return true;
  for( unsigned iq=0; iq<_nq; iq++ ) _Iq[iq] = 0.;
  _I2vec( _nq, _Iq, vec );
  return true;
}

template <typename T, typename PMT, typename PVT>
template <typename REALTYPE, typename OPT> inline bool
ODEBND_BASE<T,PMT,PVT>::_IC_I_STA
( const OPT&options, REALTYPE*vec )
{
  if( !_vIC.size() || _nx0 != _nx ) return false;
  _pIC = _vIC.at(0);
  _opIC = _pDAG->subgraph( _nx, _pIC );

  switch( options.WRAPMIT){

  case OPT::NONE:
  case OPT::DINEQ:
    _pDAG->eval( _opIC, _nx, _pIC, _Ix, _npar, _pVAR+_nx, _Ip );
    _I2vec( _nx, _Ix, vec );
    break;

  case OPT::ELLIPS:
  default:
    for( unsigned ip=0; ip<_npar; ip++ ){
      _pref[ip] = Op<T>::mid( _Ip[ip] );
      _MVPp[ip].set( _MVPenv, ip, _Ip[ip] );
    }
    _IC_I_ELL( _pDAG, _opIC, _nx, _pIC, _npar, _pVAR+_nx, _MVPp, _MVPx,
               _xref, _Q, _npar, _B, _Ix );
    _E2vec( _nx, _npar, _xref, _Q, _B, vec );
    break;
  }

  return true;
}

template <typename T, typename PMT, typename PVT>
inline void
ODEBND_BASE<T,PMT,PVT>::_IC_I_ELL
( FFGraph*DAG, std::list<const FFOp*>&opIC, const unsigned nx, const FFVar*pIC,
  const unsigned nVAR, const FFVar*pVAR, PVT*MVPVAR, PVT*MVPx, double*xref,
  double*Q, const unsigned np, double*B, T*Ix )
{
  DAG->eval( opIC, nx, pIC, MVPx, nVAR, pVAR, MVPVAR );

  double norm1R = 0.;
  for( unsigned ix=0; ix<nx; ix++ )
    norm1R += Op<T>::diam( MVPx[ix].remainder() ) / 2.;
  for( unsigned ix=0, iQ=0; ix<nx; ix++ ){
    xref[ix] = MVPx[ix].constant();
    for( unsigned jp=0; jp<np; jp++ ) B[jp*nx+ix] = MVPx[ix].linear( jp );
    Ix[ix] = MVPx[ix].bound();
    Q[iQ++] = norm1R * Op<T>::diam( MVPx[ix].remainder() ) / 2.;
    for( unsigned jx=ix+1; jx<nx; jx++ ) Q[iQ++] = 0.;
  }
#ifdef MC__ODEBND_BASE_DINEQI_DEBUG
  std::cout << "@t0" << std::endl;
  for( unsigned ix=0; ix<nx; ix++ )
    std::cout << "xref[" << ix << "] = " << xref[ix] << std::endl;
  for( unsigned ix=0; ix<nx; ix++ ){
    std::cout << "B[" << ix << ",#] = ";
    for( unsigned jp=0; jp<np; jp++ )
      std::cout << B[jp*nx+ix] << "  ";
    std::cout << std::endl;
  }
  E Ex0( nx, Q, xref );
  std::cout << "Ex0 =" << Ex0 << std::endl;
  { int dum; std::cin >> dum; }
#endif
}

template <typename T, typename PMT, typename PVT>
template <typename REALTYPE, typename OPT> inline bool
ODEBND_BASE<T,PMT,PVT>::_CC_I_STA
( const OPT&options, const unsigned iIC, const double t, REALTYPE*vec )
{
  if( _vIC.size() <= iIC || _nx0 != _nx ) return false;
  _pIC = _vIC.at( iIC );
  _opIC = _pDAG->subgraph( _nx, _pIC );
  delete[] _IIC;   _IIC = 0;
  delete[] _PMIC;  _PMIC = 0;

  switch( options.WRAPMIT){
  case OPT::NONE:
  case OPT::DINEQ:{
    _vec2I( vec, _nx, _Ix ); // current state bounds
    *_It = t; // current time
    _IIC = new T[_opIC.size()];
    _pDAG->eval( _opIC, _IIC, _nx, _pIC, _Ixdot, _nVAR-_nq, _pVAR, _IVAR );
    _I2vec( _nx, _Ixdot, vec );
    break;
   }
  case OPT::ELLIPS:
  default:{
    _vec2E( vec, _nx, _npar, _Q, _Er, _Ir, _pref, _Ip, _B, _xref, _Ix ); // current state enclosure
    for( unsigned jx=0; jx<_nx; jx++ )
      _MVXPd[jx].set( _MVXPenv, jx, _Ir[jx] );
    for( unsigned jp=0; jp<_npar; jp++ )
      _MVXPp[jp].set( _MVXPenv, _nx+jp, _Ip[jp] );
    _ep2x( _nx, _npar, _MVXPd, _pref, _MVXPp, _B, _xref, _MVXPx );
    *_MVXPt = t; // current time
    _PMIC = new PVT[_opIC.size()];
    _CC_I_ELL( _pDAG, _opIC, _PMIC, _nx, _pIC, _nVAR-_nq, _pVAR, _MVXPVAR,
       _MVXPf, _Er, _A, _npar, _xrefdot, _Bdot, _Irdot, _Qdot, options.QTOL,
       machprec() );
    _E2vec( _nx, _npar, _xrefdot, _Qdot, _Bdot, vec );
    break;
   }
  }

  return true;
}

template <typename T, typename PMT, typename PVT>
inline void
ODEBND_BASE<T,PMT,PVT>::_CC_I_ELL
( FFGraph*DAG, std::list<const FFOp*>&opf, PVT*PMf,
  const unsigned nx, const FFVar*pf, const unsigned nVAR,
  const FFVar*pVAR, PVT*MVXPVAR, PVT*MVXPf, const E&Ed, double*Af,
  const unsigned np, double*reff, double*Bf, T*Idf, double*Qf,
  const double QTOL, const double EPS )
{
#ifdef MC__ODEBND_BASE_DINEQI_DEBUG
  for( unsigned ix=0; ix<nVAR; ix++ )
    std::cout << "MVXPVAR[" << ix << "] = " << MVXPVAR[ix] << std::endl;
#endif

  // Compute polynomial expansion of function
  DAG->eval( opf, PMf, nx, pf, MVXPf, nVAR, pVAR, MVXPVAR );
#ifdef MC__ODEBND_BASE_DINEQI_DEBUG
  DAG->output( opf );

#endif

  // Extract time derivatives of constant, linear and remainder parts
  for( unsigned ix=0; ix<nx; ix++ ){
#ifdef MC__ODEBND_BASE_DINEQI_DEBUG
    std::cout << "MVXPf[" << ix << "] = " << MVXPf[ix] << std::endl;
#endif
    for( unsigned jx=0; jx<nx; jx++ )
      Af[ix+jx*nx] = MVXPf[ix].linear(jx,true);
#ifdef MC__ODEBND_BASE_DINEQI_DEBUG
    std::cout << "Af[" << ix << ",#] = ";
    for( unsigned jx=0; jx<nx; jx++ )
      std::cout << Af[ix+jx*nx] << "  ";
    std::cout << std::endl;
#endif
    for( unsigned jp=0; jp<np; jp++ )
      Bf[ix+jp*nx] = MVXPf[ix].linear(nx+jp,true);
#ifdef MC__ODEBND_BASE_DINEQI_DEBUG
    std::cout << "Bf[" << ix << ",#] = ";
    for( unsigned jp=0; jp<np; jp++ )
      std::cout << Bf[ix+jp*nx] << "  ";
    std::cout << std::endl;
#endif
#ifdef MC__ODEBND_BASE_DINEQI_DEBUG
    std::cout << "MVXPf[" << ix << "] = " << MVXPf[ix] << std::endl;
#endif
    T Rfi = MVXPf[ix].B();
    reff[ix] = Op<T>::mid(Rfi);
    Idf[ix] = Rfi - Op<T>::mid(Rfi);
#ifdef MC__ODEBND_BASE_DINEQI_DEBUG
    std::cout << "reff[" << ix << "] = " << reff[ix] << std::endl;
    std::cout << "Idf[" << ix << "] = " << Idf[ix]
              << " : " << Op<T>::mid(Idf[ix]) << std::endl;
#endif
  }
  CPPL::dgematrix matAf(nx,nx);
  for( unsigned ix=0; ix<nx; ix++ )
    for( unsigned jx=0; jx<nx; jx++ )
      matAf(ix,jx) = Af[ix+jx*nx];
  E Ef = minksum_ea( mtimes(Ed,matAf), Idf, QTOL, EPS );
  for( unsigned jx=0; jx<nx; jx++ )
    for( unsigned ix=jx; ix<nx; ix++ )
      Qf[_ndxLT(ix,jx,nx)] = Ef.Q(ix,jx);
#ifdef MC__ODEBND_BASE_DINEQI_DEBUG
  for( unsigned ix=0; ix<nx; ix++ ){
    std::cout << "Qf[" << ix << ",#] = ";
    for( unsigned jx=0; jx<=ix; jx++ )
      std::cout << Qf[_ndxLT(ix,jx,nx)] << "  ";
    std::cout << std::endl;
  }
  std::cout << "Ef =" << Ef << std::endl;
  { int dum; std::cin >> dum; }
    throw(0);
#endif
}

template <typename T, typename PMT, typename PVT>
template <typename REALTYPE, typename OPT> inline bool
ODEBND_BASE<T,PMT,PVT>::_RHS_I_STA
( const OPT&options, double t, const REALTYPE*x, REALTYPE*xdot )
{
  if( !_pRHS ) return false;

  switch( options.WRAPMIT){
  case OPT::NONE:
    _vec2I( x, _nx, _Ix );   // set current state bounds
    *_It = t; // set current time
    _RHS_I_NONE( _pDAG, _opRHS, _IRHS, _nx, _pRHS, _nVAR-_nq, _pVAR, _IVAR, _Ixdot );
    _I2vec( _nx, _Ixdot, xdot );
    return true;
   
  case OPT::DINEQ:
    if( !_opRHSi ) return false;
    _vec2I( x, _nx, _Ix );   // set current state bounds
    *_It = t; // set current time
    _RHS_I_DI( _pDAG, _opRHSi, _IRHS, _nx, _pRHS, _nVAR-_nq, _pVAR, _IVAR, _Ixdot,
               _xLdot, _xUdot );
    _I2vec( _nx, _xLdot, _xUdot, xdot );
    return true;
   
  case OPT::ELLIPS:
  default:{
    _vec2E( x, _nx, _npar, _Q, _Er, _Ir, _pref, _Ip, _B, _xref, _Ix ); // set current state enclosure
#ifdef MC__ODEBND_BASE_DINEQI_DEBUG
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
    //{ int dum; std::cin >> dum; }
#endif

    // Compute polynomial expansion of ODE's RHS
    for( unsigned jx=0; jx<_nx; jx++ )
      _MVXPd[jx].set( _MVXPenv, jx, _Ir[jx] );
    for( unsigned jp=0; jp<_npar; jp++ )
      _MVXPp[jp].set( _MVXPenv, _nx+jp, _Ip[jp] );
    _ep2x( _nx, _npar, _MVXPd, _pref, _MVXPp, _B, _xref, _MVXPx );
#ifdef MC__ODEBND_BASE_DINEQI_DEBUG
    for( unsigned ix=0; ix<_nx; ix++ )
      std::cout << "MVXPx[ " << ix << "] = " << _MVXPx[ix] << std::endl;
    { int dum; std::cin >> dum; }
#endif
    *_MVXPt = t; // set current time
    _RHS_I_ELL( _pDAG, _opRHS, _PMRHS, _nx, _pRHS, _nVAR-_nq, _pVAR, _MVXPVAR,
       _MVXPf, _Q, _A, _npar, _xrefdot, _Bdot, _Irdot, _Qdot, options.QTOL,
       machprec() );
    _E2vec( _nx, _npar, _xrefdot, _Qdot, _Bdot, xdot );
#ifdef MC__ODEBND_BASE_DINEQI_DEBUG
    //E::options.PSDCHK = true;
    for( unsigned ix=0; ix<_nx; ix++ )
      std::cout << "xrefdot[" << ix << "] = " << _xrefdot[ix] << std::endl;
    for( unsigned ix=0; ix<_nx; ix++ ){
      std::cout << "Bdot[" << ix << ",#] = ";
      for( unsigned jp=0; jp<_npar; jp++ )
        std::cout << _Bdot[jp*_nx+ix] << "  ";
      std::cout << std::endl;
    }
    std::cout << "Er = " << E(_nx,_Qdot) << std::endl;
    { int dum; std::cin >> dum; }
#endif
    
    return true;
   }
  }
}

template <typename T, typename PMT, typename PVT>
inline void
ODEBND_BASE<T,PMT,PVT>::_RHS_I_DI
( FFGraph*DAG, std::list<const FFOp*>*opRHSi, T*IRHS,
  const unsigned nx, const FFVar*pRHS, const unsigned nVAR,
  const FFVar*pVAR, T*IVAR, T*Ixdot, double*xLdot, double*xUdot )
{
  for( unsigned ix=0; ix<nx; ix++ ){
    T Ixi = IVAR[ix];
    IVAR[ix] = Op<T>::l( Ixi );
    DAG->eval( opRHSi[ix], IRHS, 1, pRHS+ix, Ixdot+ix, nVAR, pVAR, IVAR );
    xLdot[ix] = Op<T>::l( Ixdot[ix] );
    IVAR[ix] = Op<T>::u( Ixi );
    DAG->eval( opRHSi[ix], IRHS, 1, pRHS+ix, Ixdot+ix, nVAR, pVAR, IVAR );
    xUdot[ix] = Op<T>::u( Ixdot[ix] );
    IVAR[ix] = Ixi;
  }
}

template <typename T, typename PMT, typename PVT>
inline void
ODEBND_BASE<T,PMT,PVT>::_RHS_I_ELL
( FFGraph*DAG, std::list<const FFOp*>&opRHS, PVT*PMRHS,
  const unsigned nx, const FFVar*pRHS, const unsigned nVAR,
  const FFVar*pVAR, PVT*MVXPVAR, PVT*MVXPf, const double*Q, double*A,
  const unsigned np, double*xrefdot, double*Bdot, T*Iddot, double*Qdot,
  const double QTOL, const double EPS )
{
  // Compute polynomial expansion of ODE's RHS
  DAG->eval( opRHS, PMRHS, nx, pRHS, MVXPf, nVAR, pVAR, MVXPVAR );

  // Extract time derivatives of constant, linear and remainder parts
  double trQ = 0.;
  for( unsigned ix=0; ix<nx; ix++ )
    trQ += ( Q[_ndxLT(ix,ix,nx)]>0? Q[_ndxLT(ix,ix,nx)]: EPS );
  double sumkappa = 0.;
  for( unsigned ix=0; ix<nx; ix++ ){
#ifdef MC__ODEBND_BASE_DINEQI_DEBUG
    std::cout << "MVXPf[" << ix << "] = " << MVXPf[ix] << std::endl;
#endif
    for( unsigned jx=0; jx<nx; jx++ )
      A[ix+jx*nx] = MVXPf[ix].linear(jx,true);
#ifdef MC__ODEBND_BASE_DINEQI_DEBUG
    std::cout << "A[" << ix << ",#] = ";
    for( unsigned jx=0; jx<nx; jx++ )
      std::cout << A[ix+jx*nx] << "  ";
    std::cout << std::endl;
#endif
    for( unsigned jp=0; jp<np; jp++ )
      Bdot[ix+jp*nx] = MVXPf[ix].linear(nx+jp,true);
#ifdef MC__ODEBND_BASE_DINEQI_DEBUG
    std::cout << "Bdot[" << ix << ",#] = ";
    for( unsigned jp=0; jp<np; jp++ )
      std::cout << Bdot[ix+jp*nx] << "  ";
    std::cout << std::endl;
#endif
#ifdef MC__ODEBND_BASE_DINEQI_DEBUG
    std::cout << "MVXPf[" << ix << "] = " << MVXPf[ix] << std::endl;
#endif
    T Rfi = MVXPf[ix].B();
    xrefdot[ix] = Op<T>::mid(Rfi);
    Iddot[ix] = Rfi - Op<T>::mid(Rfi);
#ifdef MC__ODEBND_BASE_DINEQI_DEBUG
    std::cout << "xrefdot[" << ix << "] = " << xrefdot[ix] << std::endl;
    std::cout << "Iddot[" << ix << "] = " << Iddot[ix]
              << " : " << Op<T>::mid(Iddot[ix]) << std::endl;
    std::cout << "kappa[" << ix << "] = "
              << Op<T>::diam( Iddot[ix] ) / 2.
               / ( std::sqrt( trQ ) + QTOL ) << std::endl;
#endif
    sumkappa += ( Op<T>::diam( Iddot[ix] ) / 2. )
              / ( std::sqrt( trQ ) + QTOL );
  }

  for( unsigned jx=0; jx<nx; jx++ ){
    for( unsigned ix=jx; ix<nx; ix++ ){
      Qdot[_ndxLT(ix,jx,nx)] = sumkappa * Q[_ndxLT(ix,jx,nx)];
      for( unsigned kx=0; kx<nx; kx++ )
        Qdot[_ndxLT(ix,jx,nx)] += Q[_ndxLT(ix,kx,nx)] * A[jx+kx*nx]
                                + A[ix+kx*nx] * Q[_ndxLT(kx,jx,nx)];
    }
    Qdot[_ndxLT(jx,jx,nx)] += ( Op<T>::diam( Iddot[jx] ) / 2. )
                            * ( std::sqrt( trQ ) + QTOL );
  }
#ifdef MC__ODEBND_BASE_DINEQI_DEBUG
  for( unsigned ix=0; ix<nx; ix++ ){
    std::cout << "Qdot[" << ix << ",#] = ";
    for( unsigned jx=0; jx<=ix; jx++ )
      std::cout << Qdot[_ndxLT(ix,jx,nx)] << "  ";
    std::cout << std::endl;
  }
  E Exdot( nx, Qdot, xrefdot );
  std::cout << "Exdot =" << Exdot << std::endl;
  { int dum; std::cin >> dum; }
#endif
}

template <typename T, typename PMT, typename PVT>
template <typename REALTYPE, typename OPT> inline bool
ODEBND_BASE<T,PMT,PVT>::_RHS_I_QUAD
( const OPT&options, double t, const REALTYPE*x, REALTYPE*qdot,
  const bool bndinit )
{
  if( !_pQUAD ) return false;

  switch( options.WRAPMIT){
  case OPT::NONE:
  case OPT::DINEQ:
    if( !bndinit ) break;
    _vec2I( x, _nx, _Ix );   // set current state bounds
    break;

  case OPT::ELLIPS:
  default:
    if( !bndinit ) break;
    _vec2E( x, _nx, _npar, _Q, _Er, _Ir, _pref, _Ip, _B, _xref, _Ix ); // set current state enclosure
    break;
  }

  *_It = t; // set current time
  _QUAD_I( _pDAG, _opQUAD, _IRHS, _nq, _pQUAD, _nVAR-_nq, _pVAR, _IVAR, _Iqdot );
  _I2vec( _nq, _Iqdot, qdot );
    
   return true;
}

template <typename T, typename PMT, typename PVT>
inline void
ODEBND_BASE<T,PMT,PVT>::_QUAD_I
( FFGraph*DAG, std::list<const FFOp*>&opQUAD, T*IQUAD,
  const unsigned nq, const FFVar*pQUAD, const unsigned nVAR,
  const FFVar*pVAR, T*IVAR, T*Iqdot )
{
  if( !nq ) return;
  DAG->eval( opQUAD, IQUAD, nq, pQUAD, Iqdot, nVAR, pVAR, IVAR );
}

template <typename T, typename PMT, typename PVT>
inline void
ODEBND_BASE<T,PMT,PVT>::_RHS_I_NONE
( FFGraph*DAG, std::list<const FFOp*>&opRHS, T*IRHS,
  const unsigned nx, const FFVar*pRHS, const unsigned nVAR,
  const FFVar*pVAR, T*IVAR, T*Ixdot )
{
  DAG->eval( opRHS, IRHS, nx, pRHS, Ixdot, nVAR, pVAR, IVAR );
}

template <typename T, typename PMT, typename PVT>
template <typename REALTYPE, typename OPT> inline bool
ODEBND_BASE<T,PMT,PVT>::_JAC_I_STA
( const OPT&options, double t, const REALTYPE*x, REALTYPE*jac, REALTYPE*xdot )
{
  // Jacobian not (yet) implemented
  return false;
}

template <typename T, typename PMT, typename PVT>
inline bool
ODEBND_BASE<T,PMT,PVT>::_FCT_I_STA
( const unsigned iFCT, const double t, T*If )
{
  if( !_nf || !If ) return true;

  *_It = t; // set current time
  const FFVar* pFCT = _vFCT.at( iFCT );
  _pDAG->eval( _nf, pFCT, If, _nVAR, _pVAR, _IVAR, iFCT?true:false );
  return true;
}

template <typename T, typename PMT, typename PVT>
template <typename OPT> inline bool
ODEBND_BASE<T,PMT,PVT>::_INI_I_STA
( const OPT&options, const unsigned np, const T*Ip, const unsigned ns )
{
  // Update effective number of parameters
  // (possibly larger than _np if lifting is used)
  _npar = np;

  // Size and set DAG evaluation arrays
  _nVAR = _nx+_npar+1+_nq;
  delete[] _pVAR; _pVAR = new FFVar[_nVAR];
  delete[] _IVAR; _IVAR = new T[_nVAR];
  for( unsigned ix=0; ix<_nx; ix++ ) _pVAR[ix] = _pX[ix];
  for( unsigned ip=0; ip<_npar; ip++ ) _pVAR[_nx+ip] = _pP[ip];
  _pVAR[_nx+_npar] = (_pT? *_pT: 0. );
  for( unsigned iq=0; iq<_nq; iq++ ) _pVAR[_nx+_npar+1+iq] = _pQ?_pQ[iq]:0.;
  for( unsigned ip=0; ip<_npar; ip++ ) _IVAR[_nx+ip] = Ip[ip];
  _Ix = _IVAR;
  _Ip = _Ix + _nx;
  _It = _Ip + _npar;
  _Iq = _It + 1;

  // Reset _MVXPVenv and related variables
  if( _MVXPenv && ( _MVXPenv->nord() != options.ORDMIT 
                 || _MVXPenv->nvar() != _nx+_npar ) ){
    delete[] _MVXPf;   _MVXPf = 0;
    delete[] _MVXPd;   _MVXPd = 0;
    delete[] _MVXPVAR; _MVXPVAR = 0;
    delete   _MVXPenv; _MVXPenv = 0;
  }

  // Reset _MVPVenv and related variables
  if( _MVPenv && ( _MVPenv->nord() != 1 
                || _MVPenv->nvar() != _npar ) ){
    delete[] _MVXPx;   _MVPx = 0;
    delete[] _MVXPp;   _MVPp = 0;
    delete   _MVPenv; _MVPenv = 0;
  }

  // Set parameterization variables
  delete[] _Iqdot; _Iqdot = new T[_nq];

  switch( options.WRAPMIT){
  case OPT::DINEQ:
    delete[] _xLdot;   _xLdot = new double[_nx];
    delete[] _xUdot;   _xUdot = new double[_nx];
  case OPT::NONE:
    delete[] _Ixdot;   _Ixdot = new T[_nx];
    break;

  case OPT::ELLIPS:
  default:
    delete[] _pref;    _pref    = new double[_npar];
    delete[] _xref;    _xref    = new double[_nx];
    delete[] _A;       _A       = new double[_nx*_nx];
    delete[] _B;       _B       = new double[_nx*_npar];
    delete[] _Q;       _Q       = new double[_nx*(_nx+1)/2];
    delete[] _Ir;      _Ir      = new T[_nx];
    delete[] _xrefdot; _xrefdot = new double[_nx];
    delete[] _Bdot;    _Bdot    = new double[_nx*_npar];
    delete[] _Qdot;    _Qdot    = new double[_nx*(_nx+1)/2];
    delete[] _Irdot;   _Irdot   = new T[_nx];
    delete   _MVXPenv; _MVXPenv = new PMT( _nx+_npar, options.ORDMIT );
    _MVXPenv->options = options.PMOPT;
    delete[] _MVXPd;   _MVXPd   = new PVT[_nx];
    delete[] _MVXPf;   _MVXPf   = new PVT[_nx];
    delete[] _MVXPVAR; _MVXPVAR = new PVT[_nVAR-_nq];
    _MVXPx = _MVXPVAR;
    _MVXPp = _MVXPx + _nx;
    _MVXPt = _MVXPp +1;
    delete   _MVPenv;  _MVPenv  = new PMT( _npar, 1 );
    _MVPenv->options = options.PMOPT;
    delete[] _MVPp;    _MVPp    = new PVT[_npar];
    delete[] _MVPx;    _MVPx    = new PVT[_nx];
    break;
  }

  return true;
}

template <typename T, typename PMT, typename PVT>
template <typename OPT> inline bool
ODEBND_BASE<T,PMT,PVT>::_SET_I_STA
( const OPT&options, const unsigned iRHS, const unsigned iQUAD )
{
  if( _vRHS.size() <= iRHS ) return false;
  if( _nq && _vQUAD.size() <= iQUAD ) return false;

  _pRHS   = _vRHS.at( iRHS );
  _opRHS  = _pDAG->subgraph( _nx, _pRHS );
  _pQUAD  = _nq? _vQUAD.at( iQUAD ): 0;
  _opQUAD = _pDAG->subgraph( _nq, _pQUAD );
  const unsigned opmax = _opRHS.size()>_opQUAD.size()?_opRHS.size():_opQUAD.size();

  delete[] _IRHS;   _IRHS = 0;
  delete[] _PMRHS;  _PMRHS = 0;
  delete[] _opRHSi; _opRHSi = 0;
  switch( options.WRAPMIT){
  case OPT::NONE:
    _IRHS = new T[ opmax ];
    break;
  case OPT::DINEQ:
    _opRHSi = new std::list<const FFOp*>[_nx];
    for( unsigned ix=0; ix<_nx; ix++ )
      _opRHSi[ix] = _pDAG->subgraph( 1, _pRHS+ix );
    _IRHS = new T[ opmax ];
    break;   
  case OPT::ELLIPS:
  default:
    _IRHS = new T[ _opQUAD.size() ];
    _PMRHS = new PVT[ _opRHS.size() ];
    break;
  }

  delete[] _pJAC;
  _pJAC = 0;
  _opJAC.clear();

  return true;
}

template <typename T, typename PMT, typename PVT>
template <typename REALTYPE>
inline void
ODEBND_BASE<T,PMT,PVT>::_vec2PMI
( const REALTYPE*vec, PMT*PMenv, const unsigned nx, PVT*PMx,
  const bool centered )
{
  unsigned ivec = 0;
  for( unsigned ix=0; ix<nx; ix++  ){
    PMx[ix].set( PMenv );
    PMx[ix].set( vec+ivec );
    ivec += PMenv->nmon();
    if( !centered ){
      PMx[ix].set( T( vec[ivec], vec[ivec+1] ) );
      ivec += 2;
    }
    else{
      PMx[ix].set( T( -vec[ivec], vec[ivec] ) );
      ivec++;
    }
  }
  return;
}

template <typename T, typename PMT, typename PVT>
template <typename REALTYPE>
inline void
ODEBND_BASE<T,PMT,PVT>::_vec2PME
( const REALTYPE*vec, PMT*PMenv, const unsigned nx, PVT*PMx, double*Qr,
  E&Er, T*Ir )
{
  unsigned ivec=0;
  for( unsigned ix=0; ix<nx; ix++  ){
    PMx[ix].set( PMenv );
    PMx[ix].set( vec+ivec );
    ivec += PMenv->nmon();
  }
  for( unsigned iQ=0; iQ<(nx*(nx+1))/2; iQ++ )
    Qr[iQ] = vec[ivec++];
  _e2x( nx, Qr, Er, Ir );
  for( unsigned ix=0; ix<nx; ix++  )
    PMx[ix].set( Ir[ix] );
  return;
}

template <typename T, typename PMT, typename PVT>
template <typename U> inline void
ODEBND_BASE<T,PMT,PVT>::_e2x
( const unsigned nx, const U*d, U*x, const bool reinit )
{
  for( unsigned ix=0; ix<nx; ix++ ){   
    if( reinit ) x[ix] = 0.;
    x[ix] += d[ix];
  }
  return;
}

template <typename T, typename PMT, typename PVT>
inline void
ODEBND_BASE<T,PMT,PVT>::_e2x
( const unsigned nx, const double*Qr, E&Er, T*Ir )
{
  Er.set( nx, Qr );
#ifdef MC__ODEBND_BASE_DINEQPM_DEBUG
  std::cout << "Er =" << Er << std::endl;
#endif

  for( unsigned ix=0; ix<nx; ix++ )
    Ir[ix] = T( Er.l(ix), Er.u(ix) );

#ifdef MC__ODEBND_BASE_DINEQPM_DEBUG
  for( unsigned ix=0; ix<nx; ix++ )
    std::cout << "Ir[" << ix << "] = " << Ir[ix] << std::endl;
  { int dum; std::cin >> dum; }
#endif
  return;
}

template <typename T, typename PMT, typename PVT>
template <typename REALTYPE>
inline void
ODEBND_BASE<T,PMT,PVT>::_PMI2vec
( const PMT*PMenv, const unsigned nx, const PVT*PMx, REALTYPE*vec,
  const bool centered )
{
  unsigned ivec=0;
  for( unsigned ix=0; ix<nx; ix++ ){
    std::pair<unsigned, const double*> PMcoef = PMx[ix].coefmon();
    unsigned imon = 0;
    for( ; imon<PMcoef.first; imon++ )
      vec[ivec++] = PMcoef.second[imon];
    for( ; imon<PMenv->nmon(); imon++ )
      vec[ivec++] = 0.;
    if( !centered ){
      vec[ivec++] = Op<T>::l( PMx[ix].remainder() );
      vec[ivec++] = Op<T>::u( PMx[ix].remainder() );
    }
    else{
      vec[ivec++] = 0.5*Op<T>::diam( PMx[ix].remainder() );
    }
  }
  return;
}

template <typename T, typename PMT, typename PVT>
template <typename REALTYPE>
inline void
ODEBND_BASE<T,PMT,PVT>::_PMI2vec
( const PMT*PMenv, const unsigned nx, const PVT*PMx, const double*RxL,
  const double*RxU, REALTYPE*vec )
{
  unsigned ivec=0;
  for( unsigned ix=0; ix<nx; ix++ ){
    std::pair<unsigned, const double*> PMcoef = PMx[ix].coefmon();
    unsigned imon = 0;
    for( ; imon<PMcoef.first; imon++ )
      vec[ivec++] = PMcoef.second[imon];
    for( ; imon<PMenv->nmon(); imon++ )
      vec[ivec++] = 0.;
    vec[ivec++] = ( RxL? RxL[ix]: 0. );
    vec[ivec++] = ( RxU? RxU[ix]: 0. );
  }
  return;
}

template <typename T, typename PMT, typename PVT>
template <typename REALTYPE>
inline void
ODEBND_BASE<T,PMT,PVT>::_PMI2vec
( const PMT*PMenv, const unsigned nx, const PVT*PMx, const T*Rx,
  REALTYPE*vec )
{
  unsigned ivec=0;
  for( unsigned ix=0; ix<nx; ix++ ){
    std::pair<unsigned, const double*> PMcoef = PMx[ix].coefmon();
    unsigned imon = 0;
    for( ; imon<PMcoef.first; imon++ )
      vec[ivec++] = PMcoef.second[imon];
    for( ; imon<PMenv->nmon(); imon++ )
      vec[ivec++] = 0.;
    vec[ivec++] = ( Rx? Op<T>::l( Rx[ix] ): 0. );
    vec[ivec++] = ( Rx? Op<T>::u( Rx[ix] ): 0. );
  }
  return;
}

template <typename T, typename PMT, typename PVT>
template <typename REALTYPE>
inline void
ODEBND_BASE<T,PMT,PVT>::_PME2vec
( const PMT*PMenv, const unsigned nx, const PVT*PMx, const double*QRx,
  REALTYPE*vec )
{
  unsigned ivec=0;
  for( unsigned ix=0; ix<nx; ix++ ){
    std::pair<unsigned, const double*> PMcoef = PMx[ix].coefmon();
    unsigned imon = 0;
    for( ; imon<PMcoef.first; imon++ )
      vec[ivec++] = PMcoef.second[imon];
    for( ; imon<PMenv->nmon(); imon++ )
      vec[ivec++] = 0.;
  }
  for( unsigned iQ=0; iQ<nx*(nx+1)/2; iQ++ )
    vec[ivec++] = ( QRx? QRx[iQ]: 0. );
  return;
}

template <typename T, typename PMT, typename PVT>
template <typename REALTYPE, typename OPT> inline bool
ODEBND_BASE<T,PMT,PVT>::_IC_PM_QUAD
( const OPT&options, REALTYPE*vec )
{
  if( !_vQUAD.size() || !_nq ) return true;
  for( unsigned iq=0; iq<_nq; iq++ ) _PMq[iq] = 0.;
  _PMI2vec( _PMenv, _nq, _PMq, vec, true );
  return true;
}

template <typename T, typename PMT, typename PVT>
template <typename REALTYPE, typename OPT> inline bool
ODEBND_BASE<T,PMT,PVT>::_IC_PM_STA
( const OPT&options, REALTYPE*vec )
{
  if( !_vIC.size() || _nx0 != _nx ) return false;
  _pIC = _vIC.at(0);
  _opIC = _pDAG->subgraph( _nx, _pIC );

  switch( options.WRAPMIT){

  case OPT::NONE:
  case OPT::DINEQ:
    _pDAG->eval( _opIC, _nx, _pIC, _PMx, _npar, _pVAR+_nx, _PMp );
    // Whether or not to ignore the remainder
    if( !options.PMNOREM )
      _PMI2vec( _PMenv, _nx, _PMx, vec, false );
    else
      _PMI2vec( _PMenv, _nx, _PMx, 0, vec );
    break;

  case OPT::ELLIPS:
  default:
    _IC_PM_ELL( _pDAG, _opIC, _nx, _pIC, _npar, _pVAR+_nx, _PMp, _PMx, _Q, _Er, _Ir );
    // Whether or not to ignore the remainder
    if( !options.PMNOREM )
      _PME2vec( _PMenv, _nx, _PMx, _Q, vec );
    else
      _PME2vec( _PMenv, _nx, _PMx, 0, vec );
    break;
  }
  return true;
}

template <typename T, typename PMT, typename PVT>
inline void
ODEBND_BASE<T,PMT,PVT>::_IC_PM_ELL
( FFGraph*DAG, std::list<const FFOp*>&opIC, const unsigned nx, const FFVar*pIC,
  const unsigned np, const FFVar*pP, PVT*PMp, PVT*PMx, double*Qr, E&Er, T*Ir )
{
  DAG->eval( opIC, nx, pIC, PMx, np, pP, PMp );

  double norm1R = 0.;
  for( unsigned ix=0; ix<nx; ix++ )
    norm1R += Op<T>::diam( PMx[ix].remainder() ) / 2.;
  for( unsigned ix=0, iQ=0; ix<nx; ix++ ){
    Qr[iQ++] = norm1R * Op<T>::diam( PMx[ix].remainder() ) / 2.;
    for( unsigned jx=ix+1; jx<nx; jx++ ) Qr[iQ++] = 0.;
    Ir[ix] = PMx[ix].remainder();
  }
  Er.set( nx, Qr );
#ifdef MC__ODEBND_BASE_DINEQPM_DEBUG
  std::cout << "ERx0 =" << Er << std::endl;
#endif
}

template <typename T, typename PMT, typename PVT>
template <typename REALTYPE, typename OPT> inline bool
ODEBND_BASE<T,PMT,PVT>::_CC_PM_STA
( const OPT&options, const unsigned iIC, const double t, REALTYPE*vec )
{
  if( _vIC.size() <= iIC || _nx0 != _nx ) return false;
  _pIC = _vIC.at( iIC );
  _opIC = _pDAG->subgraph( _nx, _pIC );
  delete[] _IIC;   _IIC = 0;
  delete[] _PMIC;  _PMIC = 0;
  delete[] _IDIC;  _IDIC = 0;
  delete[] _PMDIC; _PMDIC = 0;
  delete[] _pDIC;  _pDIC = 0; _opDIC.clear();

  switch( options.WRAPMIT){
  case OPT::NONE:
  case OPT::DINEQ:
    _vec2PMI( vec, _PMenv, _nx, _PMx ); // current state/quadrature polynomial model
    *_PMt = t; // current time
    _PMIC = new PVT[_opIC.size()];
    _pDAG->eval( _opIC, _PMIC, _nx, _pIC, _PMxdot, _nVAR-_nq, _pVAR, _PMVAR );
    // Whether or not to ignore the remainder
    if( !options.PMNOREM )
      _PMI2vec( _PMenv, _nx, _PMxdot, vec, false );
    else
      _PMI2vec( _PMenv, _nx, _PMxdot, 0, vec );
    break;
   
  case OPT::ELLIPS:
  default:{
    _vec2PME( vec, _PMenv, _nx, _PMx, _Q, _Er, _Ir ); // current state/quadrature polynomial model
    *_PMt = t; // current time   
    _PMIC = new PVT[_opIC.size()];

    // In this variant a bound on the Jacobian matrix is computed and the
    // linear part is taken as the mid-point of this matrix
    if( !options.ORDMIT ){
      for( unsigned ix=0; ix<_nx; ix++ ){
        _IVAR[ix] = _PMx[ix].bound(); // set current state bounds
        _PMVAR[ix].center().set( T(0.) ); // cancel remainder term
      }
      *_It = t; // current time
      _pDIC = _pDAG->FAD( _nx, _pIC, _nx, _pVAR );
      _opDIC = _pDAG->subgraph( _nx*_nx, _pDIC );
      _IDIC = new T[ _opDIC.size() ];
      _RHS_PM_ELL0( _pDAG, _opIC, _PMIC, _opDIC, _IDIC, _nx, _pIC,
                    _pDIC, _nVAR-_nq, _pVAR, _PMVAR, _IVAR, _PMxdot, _Idfdx,
                    _Ir, _A, _Irdot );
    }

    // In this variant a polynomial model of the Jacobian matrix is computed and the
    // linear part is taken as the mid-point of this matrix
    else if( _PMenv->nord() >= _MVXPenv->nord() ){
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
      *_MVXPt = t; // current time
      _pDIC = _pDAG->FAD( _nx, _pIC, _nx, _pVAR );
      _opDIC = _pDAG->subgraph( _nx*_nx, _pDIC );
      _PMDIC = new PVT[ _opDIC.size() ];
      _RHS_PM_ELL1( _pDAG, _opIC, _PMIC, _opDIC, _PMDIC, _nx, _pIC,
                    _pDIC, _nVAR-_nq, _pVAR, _PMVAR, _MVXPVAR, _PMxdot,
                    _MVXPdfdx, _Ir, _A, _Irdot );
    }

    // In this variant a polynomial model in the joint state-parameter and
    // of the same order as the parameter polynomial model is computed
    else{
      for( unsigned jx=0; jx<_nx; jx++ ){
        _MVXPd[jx].set( _MVXPenv, _npar+jx, _Ir[jx] );
        _MVXPx[jx].set( _MVXPenv ).set( _PMx[jx].center().set( T(0.) ), true );
      }
      _e2x( _nx, _MVXPd, _MVXPx, false );
      *_MVXPt = t; // current time   
      _RHS_PM_ELL2( _pDAG, _opIC, _PMIC, _nx, _pRHS, _nVAR-_nq, _pVAR,
                    _MVXPVAR, _PMenv, _PMxdot, _MVXPf, _npar, _Ir, _A, _Irdot );
    }

    // Whether or not to ignore the remainder
    if( !options.PMNOREM ){
      _CC_PM_ELL( _nx, _Er, _A, _Irdot, _Qdot, options.QTOL, machprec() );
      _PME2vec( _PMenv, _nx, _PMxdot, _Qdot, vec );
    }
    else
      _PME2vec( _PMenv, _nx, _PMxdot, 0, vec );
    break;
   }
  }

  return true;
}

template <typename T, typename PMT, typename PVT>
inline void
ODEBND_BASE<T,PMT,PVT>::_CC_PM_ELL
( const unsigned nx, const E&Exr, const double*Afr, const T*Ifr,
  double*Qfr, const double QTOL, const double EPS )
{
  // Set shape matrix of discontinuity
  E Efr;
  if( Afr ){
    CPPL::dgematrix matAfr(nx,nx);
    for( unsigned ix=0; ix<nx; ix++ )
      for( unsigned jx=0; jx<nx; jx++ )
        matAfr(ix,jx) = Afr[ix+jx*nx];
    Efr = minksum_ea( mtimes(Exr,matAfr), Ifr, QTOL, EPS );
  }
  else
    Efr = minksum_ea( Exr, Ifr, QTOL, EPS );    
  for( unsigned jx=0; jx<nx; jx++ )
    for( unsigned ix=jx; ix<nx; ix++ )
      Qfr[_ndxLT(ix,jx,nx)] = Efr.Q(ix,jx);

#ifdef MC__ODEBND_BASE_DINEQPM_DEBUG
  for( unsigned ix=0; ix<nx; ix++ ){
    std::cout << "Qfr[" << ix << ",#] = ";
    for( unsigned jx=0; jx<=ix; jx++ )
      std::cout << Qfr[_ndxLT(ix,jx,nx)] << "  ";
    std::cout << std::endl;
  }
  std::cout << "Efr =" << Efr << std::endl;
  { int dum; std::cin >> dum; }
#endif
}

template <typename T, typename PMT, typename PVT>
template <typename REALTYPE, typename OPT> inline bool
ODEBND_BASE<T,PMT,PVT>::_RHS_PM_STA
( const OPT&options, double t, const REALTYPE*x, REALTYPE*xdot )
{
  if( !_pRHS ) return false;

  switch( options.WRAPMIT){
  case OPT::NONE:
    _vec2PMI( x, _PMenv, _nx, _PMx );   // set current state polynomial model
#ifdef MC__ODEBND_BASE_DINEQPM_DEBUG
    _print_interm( t, _nx, _PMx, "PMx Intermediate", std::cerr );
#endif
    *_PMt = t; // set current time
    _QUAD_PM( _pDAG, _opQUAD, _PMRHS, _nq, _pQUAD, _nVAR-_nq, _pVAR, _PMVAR,
              _PMqdot );
    _RHS_PM_NONE( _pDAG, _opRHS, _PMRHS, _nx, _pRHS, _nVAR-_nq, _pVAR, _PMVAR,
                  _PMxdot );
    // Whether or not to ignore the remainder
    if( !options.PMNOREM )
      _PMI2vec( _PMenv, _nx, _PMxdot, xdot, false );
    else
      _PMI2vec( _PMenv, _nx, _PMxdot, 0, xdot );
#ifdef MC__ODEBND_BASE_DINEQPM_DEBUG
    _print_interm( t, _nx, _PMxdot, "PMxdot Intermediate", std::cerr );
    {int dum; std::cin >> dum;}
#endif
    return true;  

  case OPT::DINEQ:
    _vec2PMI( x, _PMenv, _nx, _PMx );   // set current state polynomial model
    for( unsigned ix=0; ix<_nx; ix++ ) _PMx[ix].center();
#ifdef MC__ODEBND_BASE_DINEQPM_DEBUG
    _print_interm( t, _nx, _PMx, "PMx Intermediate", std::cerr );
#endif
    *_PMt = t; // set current time
    _QUAD_PM( _pDAG, _opQUAD, _PMRHS, _nq, _pQUAD, _nVAR-_nq, _pVAR, _PMVAR,
              _PMqdot );
    _RHS_PM_DI( _pDAG, _opRHSi, _PMRHS, _nx, _pRHS, _nVAR-_nq, _pVAR, _PMVAR,
                _PMxdot, _RxLdot, _RxUdot );
    // Whether or not to ignore the remainder
    if( !options.PMNOREM )
      _PMI2vec( _PMenv, _nx, _PMxdot, _RxLdot, _RxUdot, xdot );
    else
      _PMI2vec( _PMenv, _nx, _PMxdot, 0, xdot );
#ifdef MC__ODEBND_BASE_DINEQPM_DEBUG
    _print_interm( t, _nx, _PMxdot, "PMxdot Intermediate", std::cerr );
    {int dum; std::cin >> dum;}
#endif
    return true;  
   
  case OPT::ELLIPS:
  default:{
    _vec2PME( x, _PMenv, _nx, _PMx, _Q, _Er, _Ir ); // set current state polynomial model
    *_PMt = t; // set current time   
#ifdef MC__ODEBND_BASE_DINEQPM_DEBUG
    _print_interm( t, _nx, _PMx, _Er, "PMx Intermediate", std::cerr );
#endif
    _QUAD_PM( _pDAG, _opQUAD, _PMRHS, _nq, _pQUAD, _nVAR-_nq, _pVAR, _PMVAR,
              _PMqdot );

    // In this variant a bound on the Jacobian matrix is computed and the
    // linear part is taken as the mid-point of this matrix
    if( !options.ORDMIT ){
      for( unsigned ix=0; ix<_nx; ix++ ){
        _IVAR[ix] = _PMx[ix].bound(); // set current state bounds
        _PMVAR[ix].center().set( T(0.) ); // cancel remainder term
      }
      *_It = t; // set current time
      _RHS_PM_ELL0( _pDAG, _opRHS, _PMRHS, _opJAC, _IJAC, _nx, _pRHS,
                    _pJAC, _nVAR-_nq, _pVAR, _PMVAR, _IVAR, _PMxdot, _Idfdx,
                    _Ir, _A, _Irdot );
    }

    // In this variant a polynomial model of the Jacobian matrix is computed and the
    // linear part is taken as the mid-point of this matrix
    else if( _PMenv->nord() >= _MVXPenv->nord() ){
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
#ifdef MC__ODEBND_BASE_DINEQPM_DEBUG
      _print_interm( t, _nx, _MVXPx, "MVXPx Intermediate", std::cerr );
#endif
      *_MVXPt = t; // set current time
      _RHS_PM_ELL1( _pDAG, _opRHS, _PMRHS, _opJAC, _PMJAC, _nx, _pRHS,
                    _pJAC, _nVAR-_nq, _pVAR, _PMVAR, _MVXPVAR, _PMxdot,
                    _MVXPdfdx, _Ir, _A, _Irdot );
    }

    // In this variant a polynomial model in the joint state-parameter and
    // of the same order as the parameter polynomial model is computed
    else{
      for( unsigned jx=0; jx<_nx; jx++ ){
        _MVXPd[jx].set( _MVXPenv, _npar+jx, _Ir[jx] );
        _MVXPx[jx].set( _MVXPenv ).set( _PMx[jx].center().set( T(0.) ), true );
      }
      _e2x( _nx, _MVXPd, _MVXPx, false );
#ifdef MC__ODEBND_BASE_DINEQPM_DEBUG
      _print_interm( t, _nx, _MVXPx, "MVXPx Intermediate", std::cerr );
#endif
      *_MVXPt = t; // set current time   
      _RHS_PM_ELL2( _pDAG, _opRHS, _PMRHS, _nx, _pRHS, _nVAR-_nq, _pVAR,
                    _MVXPVAR, _PMenv, _PMxdot, _MVXPf, _npar, _Ir, _A, _Irdot );
    }

    _RHS_PM_ELL( _nx, _Q, _A, _Irdot, _Qdot, options.QTOL, machprec() );
    // Whether or not to ignore the remainder
    if( !options.PMNOREM )
      _PME2vec( _PMenv, _nx, _PMxdot, _Qdot, xdot );
    else
      _PME2vec( _PMenv, _nx, _PMxdot, 0, xdot );
#ifdef MC__ODEBND_BASE_DINEQPM_DEBUG
    _print_interm( t, _nx, _PMxdot, "PMxdot Intermediate", std::cerr );
    {int dum; std::cin >> dum;}
#endif
    return true;
   }
  }
}

template <typename T, typename PMT, typename PVT>
inline void
ODEBND_BASE<T,PMT,PVT>::_RHS_PM_NONE
( FFGraph*DAG, std::list<const FFOp*>&opRHS, PVT*PMRHS,
  const unsigned nx, const FFVar*pRHS, const unsigned nVAR,
  const FFVar*pVAR, PVT*PMVAR, PVT*PMxdot )
{
  DAG->eval( opRHS, PMRHS, nx, pRHS, PMxdot, nVAR, pVAR, PMVAR );
}

template <typename T, typename PMT, typename PVT>
inline void
ODEBND_BASE<T,PMT,PVT>::_RHS_PM_DI
( FFGraph*DAG, std::list<const FFOp*>*opRHSi, PVT*PMRHS,
  const unsigned nx, const FFVar*pRHS, const unsigned nVAR,
  const FFVar*pVAR, PVT*PMVAR, PVT*PMxdot, double*rLdot, double*rUdot )
{
  for( unsigned ix=0; ix<nx; ix++ ){
    T Rxi = PMVAR[ix].remainder();
    PMVAR[ix].set( Op<T>::l( Rxi ) );
    DAG->eval( opRHSi[ix], PMRHS, 1, pRHS+ix, PMxdot+ix, nVAR, pVAR, PMVAR );
    rLdot[ix] = Op<T>::l( PMxdot[ix].remainder() );
    PMVAR[ix].set( Op<T>::u( Rxi ) );
    DAG->eval( opRHSi[ix], PMRHS, 1, pRHS+ix, PMxdot+ix, nVAR, pVAR, PMVAR );
    rUdot[ix] = Op<T>::u( PMxdot[ix].remainder() );
    PMVAR[ix].set( Rxi );
  }
}

template <typename T, typename PMT, typename PVT>
inline void
ODEBND_BASE<T,PMT,PVT>::_RHS_PM_ELL
( const unsigned nx, const double*Qr, const double*Ar, const T*Irdot,
  double*Qrdot, const double QTOL, const double EPS )
{
  // Set dynamics of shape matrix
  double trQ = 0., sumkappa = 0.;
  for( unsigned ix=0; ix<nx; ix++ )
    trQ += ( Qr[_ndxLT(ix,ix,nx)]>0? Qr[_ndxLT(ix,ix,nx)]: EPS );
  const double srqt_trQ = (trQ>0? std::sqrt( trQ ): 0.) + QTOL;
  for( unsigned ix=0; ix<nx; ix++ ){
#ifdef MC__ODEBND_BASE_DINEQPM_DEBUG
      std::cout << "kappa[" << ix << "] = "
                << Op<T>::diam( Irdot[ix] ) / ( 2. * srqt_trQ ) << std::endl;
#endif
      sumkappa += Op<T>::diam( Irdot[ix] ) / ( 2. * srqt_trQ );
  }

  for( unsigned jx=0; jx<nx; jx++ ){
    for( unsigned ix=jx; ix<nx; ix++ ){
      Qrdot[_ndxLT(ix,jx,nx)] = sumkappa * Qr[_ndxLT(ix,jx,nx)];
      for( unsigned kx=0; kx<nx; kx++ )
        Qrdot[_ndxLT(ix,jx,nx)] += Qr[_ndxLT(ix,kx,nx)] * Ar[jx+kx*nx]
                                 + Ar[ix+kx*nx] * Qr[_ndxLT(kx,jx,nx)];
    }
    Qrdot[_ndxLT(jx,jx,nx)] += Op<T>::diam( Irdot[jx] ) / 2. * srqt_trQ;
  }

#ifdef MC__ODEBND_BASE_DINEQPM_DEBUG
  for( unsigned ix=0; ix<nx; ix++ )
    std::cout << "Irdot[" << ix << "] =" << Irdot[ix] << std::endl;
  E Erdot( nx, Qrdot );
  std::cout << "Erdot =" << Erdot << std::endl;
  { int dum; std::cin >> dum; }
#endif
}

template <typename T, typename PMT, typename PVT>
inline void
ODEBND_BASE<T,PMT,PVT>::_RHS_PM_ELL0
( FFGraph*DAG, std::list<const FFOp*>&opRHS, PVT*PMRHS,
  std::list<const FFOp*>&opJAC, T*IJAC, const unsigned nx,
  const FFVar*pRHS, const FFVar*pJAC, const unsigned nVAR,
  const FFVar*pVAR, PVT*PMVAR, T*IVAR, PVT*PMf, T*Idfdx,
  const T*Ir, double*Ar, T*Irdot )
{
  DAG->eval( opRHS, PMRHS, nx, pRHS, PMf, nVAR, pVAR, PMVAR );
  DAG->eval( opJAC, IJAC, nx*nx, pJAC, Idfdx, nVAR, pVAR, IVAR );

  for( unsigned ix=0; ix<nx; ix++ ){
    for( unsigned jx=0; jx<nx; jx++ )
      Ar[ix+jx*nx] = Op<T>::mid( Idfdx[ix*nx+jx] );
    Irdot[ix] = PMf[ix].remainder();
    for( unsigned jx=0; jx<nx; jx++ )
      Irdot[ix] += ( Idfdx[ix*nx+jx] - Ar[ix+jx*nx] ) * Ir[jx];
#ifdef MC__ODEBND_BASE_DINEQPM_DEBUG
    std::cout << "Irdot[" << ix << "] = " << Irdot[ix]
              << " : " << Op<T>::mid(Irdot[ix]) << std::endl;
#endif
  }
}

template <typename T, typename PMT, typename PVT>
inline void
ODEBND_BASE<T,PMT,PVT>::_RHS_PM_ELL1
( FFGraph*DAG, std::list<const FFOp*>&opRHS, PVT*PMRHS,
  std::list<const FFOp*>&opJAC, PVT*PMJAC, const unsigned nx,
  const FFVar*pRHS, const FFVar*pJAC, const unsigned nVAR,
  const FFVar*pVAR, PVT*PMVAR, PVT*MVXPVAR, PVT*PMf, PVT*MVXPdfdx,
  const T*Ir, double*Ar, T*Irdot )
{
  DAG->eval( opRHS, PMRHS, nx, pRHS, PMf, nVAR, pVAR, PMVAR );
  DAG->eval( opJAC, PMJAC, nx*nx, pJAC, MVXPdfdx, nVAR, pVAR, MVXPVAR );

  for( unsigned ix=0; ix<nx; ix++ ){
#ifdef MC__ODEBND_BASE_DINEQPM_DEBUG
    for( unsigned jx=0; jx<nx; jx++ )
      std::cout << "PMdfdx[" << ix << "," << jx << "] = " << MVXPdfdx[ix*nx+jx] << std::endl;
#endif
    // Extract constant coefficients and set to 0
    for( unsigned jx=0; jx<nx; jx++ )
      Ar[ix+jx*nx] = ( MVXPdfdx[ix*nx+jx].center() ).constant( true );
#ifdef MC__ODEBND_BASE_DINEQPM_DEBUG
    std::cout << "Ar[" << ix << ",#] = ";
    for( unsigned jx=0; jx<nx; jx++ )
      std::cout << Ar[ix+jx*nx] << "  ";
    std::cout << std::endl;
    for( unsigned jx=0; jx<nx; jx++ )
      std::cout << "PMdfdx[" << ix << "," << jx << "] = " << MVXPdfdx[ix*nx+jx] << std::endl;
#endif
    // Bound remaining terms
    Irdot[ix] = PMf[ix].remainder();
    for( unsigned jx=0; jx<nx; jx++ )
      Irdot[ix] += MVXPdfdx[ix*nx+jx].bound() * Ir[jx];
#ifdef MC__ODEBND_BASE_DINEQPM_DEBUG
    std::cout << "Iddot[" << ix << "] = " << Irdot[ix]
              << " : " << Op<T>::mid(Irdot[ix]) << std::endl;
#endif
  }
}

template <typename T, typename PMT, typename PVT>
inline void
ODEBND_BASE<T,PMT,PVT>::_RHS_PM_ELL2
( FFGraph*DAG, std::list<const FFOp*>&opRHS, PVT*PMRHS,
  const unsigned nx, const FFVar*pRHS, const unsigned nVAR, const FFVar*pVAR,
  PVT*MVXPVAR, PMT*PMenv, PVT*PMf, PVT*MVXPf, const unsigned np, const T*Ir,
  double*Ar, T*Irdot )
{
  DAG->eval( opRHS, PMRHS, nx, pRHS, MVXPf, nVAR, pVAR, MVXPVAR );

  for( unsigned ix=0; ix<nx; ix++ ){
#ifdef MC__ODEBND_BASE_DINEQPM_DEBUG
    std::cout << "MVXPf[" << ix << "] = " << MVXPf[ix] << std::endl;
#endif
    // Extract polynomial model in P and set to 0
    MVXPf[ix].get( PMf[ix].set(PMenv), true );
#ifdef MC__ODEBND_BASE_DINEQPM_DEBUG
    std::cout << "PMxdot[" << ix << "] =" << PMf[ix] << std::endl;
    std::cout << "MVXPf[" << ix << "] = " << MVXPf[ix] << std::endl;
    { int dum; std::cin >> dum; }
#endif
    // Extract linear part in X and set to 0
    for( unsigned jx=0; jx<nx; jx++ )
      Ar[ix+jx*nx] = MVXPf[ix].linear( np+jx, true );
#ifdef MC__ODEBND_BASE_DINEQPM_DEBUG
    std::cout << "Ar[" << ix << ",#] = ";
    for( unsigned jx=0; jx<nx; jx++ )
      std::cout << Ar[ix+jx*nx] << "  ";
    std::cout << std::endl;
    std::cout << "MVXPf[" << ix << "] = " << MVXPf[ix] << std::endl;
#endif
    // Bound remaining terms
    Irdot[ix] = MVXPf[ix].bound();
#ifdef MC__ODEBND_BASE_DINEQPM_DEBUG
    std::cout << "Irdot[" << ix << "] = " << Irdot[ix]
              << " : " << Op<T>::mid(Irdot[ix]) << std::endl;
#endif
  }
}

template <typename T, typename PMT, typename PVT>
template <typename REALTYPE, typename OPT> inline bool
ODEBND_BASE<T,PMT,PVT>::_RHS_PM_QUAD
( const OPT&options, double t, const REALTYPE*x, REALTYPE*qdot,
  const bool reinit )
{
  if( !_pQUAD ) return false;

  switch( options.WRAPMIT){
  case OPT::NONE:
  case OPT::DINEQ:
    if( !reinit ) break;
    _vec2PMI( x, _PMenv, _nx, _PMx );   // set current state polynomial model
#ifdef MC__ODEBND_BASE_DINEQPM_DEBUG
    _print_interm( t, _nx, _PMx, "PMx Intermediate", std::cerr );
#endif
    break;
   
  case OPT::ELLIPS:
  default:
    if( !reinit ) break;
    _vec2PME( x, _PMenv, _nx, _PMx, _Q, _Er, _Ir ); // set current state polynomial model
#ifdef MC__ODEBND_BASE_DINEQPM_DEBUG
    _print_interm( t, _nx, _PMx, _Er, "PMx Intermediate", std::cerr );
#endif
    break;
  }

  *_PMt = t; // set current time
  _QUAD_PM( _pDAG, _opQUAD, _PMRHS, _nq, _pQUAD, _nVAR-_nq, _pVAR, _PMVAR,
            _PMqdot );
#ifdef MC__ODEBND_BASE_DINEQPM_DEBUG
    _print_interm( t, _nq, _PMqdot, "PMqdot Intermediate", std::cerr );
    {int dum; std::cin >> dum;}
#endif

  // Whether or not to ignore the remainder
  if( !options.PMNOREM )
    _PMI2vec( _PMenv, _nq, _PMqdot, qdot, true );
  else
    _PMI2vec( _PMenv, _nq, _PMqdot, 0, qdot );

  return true;  
}

template <typename T, typename PMT, typename PVT>
inline void
ODEBND_BASE<T,PMT,PVT>::_QUAD_PM
( FFGraph*DAG, std::list<const FFOp*>&opQUAD, PVT*PMQUAD,
  const unsigned nq, const FFVar*pQUAD, const unsigned nVAR,
  const FFVar*pVAR, PVT*PMVAR, PVT*PMqdot )
{
  if( !nq ) return;
  DAG->eval( opQUAD, PMQUAD, nq, pQUAD, PMqdot, nVAR, pVAR, PMVAR );
  for( unsigned iq=0; iq<nq; iq++ ) PMqdot[iq].center();
}

template <typename T, typename PMT, typename PVT>
template <typename REALTYPE, typename OPT> inline bool
ODEBND_BASE<T,PMT,PVT>::_JAC_PM_STA
( const OPT&options, double t, const REALTYPE*x, REALTYPE*jac, REALTYPE*xdot )
{
  return false;
}

template <typename T, typename PMT, typename PVT> inline bool
ODEBND_BASE<T,PMT,PVT>::_FCT_PM_STA
( const unsigned iFCT, const double t, PVT*PMf )
{
  if( !_nf || !PMf ) return true;

  *_PMt = t; // set current time
  const FFVar* pFCT = _vFCT.at( iFCT );
  _pDAG->eval( _nf, pFCT, PMf, _nVAR, _pVAR, _PMVAR, iFCT?true:false );
  return true;
}

template <typename T, typename PMT, typename PVT>
template <typename OPT> inline bool
ODEBND_BASE<T,PMT,PVT>::_INI_PM_STA
( const OPT&options, const unsigned np, const PVT* PMp, const unsigned ns )
{
  // Update effective number of parameters
  // (possibly larger than _np due to lifting)
  _npar = np;

  // Check polynomial model compatibility and size
  unsigned kp=_npar;
  for( unsigned ip=0; ip<_npar && kp==_npar; ip++ )
    if( PMp[ip].env() ) kp = ip;
  if( kp==_npar || PMp[kp].env()->nvar()!=_npar ) return false;
  _PMenv = PMp[kp].env();  

  // Size and set DAG evaluation arrays
  _nVAR = _nx+_npar+1+_nq;
  delete[] _pVAR; _pVAR = new FFVar[_nVAR];
  delete[] _IVAR; _IVAR = new T[_nVAR];
  delete[] _PMVAR; _PMVAR = new PVT[_nVAR];
  for( unsigned ix=0; ix<_nx; ix++ ) _pVAR[ix] = _pX[ix];
  for( unsigned ip=0; ip<_npar; ip++ ) _pVAR[_nx+ip] = _pP[ip];
  _pVAR[_nx+_npar] = (_pT? *_pT: 0. );
  for( unsigned iq=0; iq<_nq; iq++ ) _pVAR[_nx+_npar+1+iq] = _pQ?_pQ[iq]:0.;
  for( unsigned ip=0; ip<_npar; ip++ ) _IVAR[_nx+ip] = PMp[ip].bound();
  _Ix = _IVAR;
  _Ip = _Ix + _nx;
  _It = _Ip + _npar;
  _Iq = _It + 1;
  for( unsigned ip=0; ip<_npar; ip++ ) _PMVAR[_nx+ip] = PMp[ip];
  _PMx = _PMVAR;
  _PMp = _PMx + _nx;
  _PMt = _PMp + _npar;
  _PMq = _PMt + 1;

  // Reset _MVXPVenv and related variables
  unsigned MVXPsize = ( options.ORDMIT<_PMenv->nord()? options.ORDMIT: _PMenv->nord() ); 
#ifdef MC__ODEBND_BASE_MVXP_USE
  unsigned MVXPdim  = _nx+_npar; 
#else
  unsigned MVXPdim  = ( options.ORDMIT<_PMenv->nord()? _npar: _nx+_npar ); 
#endif
  if( _MVXPenv && ( _MVXPenv->nord() != MVXPsize || _MVXPenv->nvar() != MVXPdim ) ){
    delete[] _MVXPf;    _MVXPf = 0;
    delete[] _MVXPdfdx; _MVXPdfdx = 0;
    delete[] _MVXPd;    _MVXPd = 0;
    delete[] _MVXPVAR;  _MVXPVAR = _MVXPx = _MVXPp = 0; 
    delete   _MVXPenv;  _MVXPenv = 0;
  }

  // Size state/quadrature derivative arrays
  delete[] _PMxdot; _PMxdot = new PVT[_nx];
  delete[] _PMqdot; _PMqdot = _nq? new PVT[_nq]: 0;
  delete[] _radRqdot;  _radRqdot  = _nq? new double[_nq]: 0;

  switch( options.WRAPMIT){
  case OPT::NONE:
    break;

  case OPT::DINEQ:
    delete[] _RxLdot; _RxLdot = new double[_nx];
    delete[] _RxUdot; _RxUdot = new double[_nx];
    break;

  case OPT::ELLIPS:
    delete[] _xref;     _xref     = new double[_nx];
    delete[] _A;        _A        = new double[_nx*_nx];
    delete[] _Q;        _Q        = new double[_nx*(_nx+1)/2];
    delete[] _Ir;       _Ir       = new T[_nx];
    delete[] _Qdot;     _Qdot     = new double[_nx*(_nx+1)/2];
    delete[] _Irdot;    _Irdot    = new T[_nx];
    delete[] _Idfdx;    _Idfdx    = new T[_nx*_nx];
    delete   _MVXPenv;  _MVXPenv  = new PMT( MVXPdim, MVXPsize );
    _MVXPenv->options = _PMenv->options;
    delete[] _MVXPd;    _MVXPd    = new PVT[_nx];
    delete[] _MVXPf;    _MVXPf    = new PVT[_nx];
    delete[] _MVXPdfdx; _MVXPdfdx = new PVT[_nx*_nx];
    delete[] _MVXPVAR;  _MVXPVAR  = new PVT[_nVAR];
    _MVXPx = _MVXPVAR;
    _MVXPp = _MVXPx + _nx;
    _MVXPt = _MVXPp + _npar;
    for( unsigned ip=0; ip<_npar; ip++ )
      _MVXPp[ip].set( _MVXPenv, ip, _PMp[ip].B() );
    break;
  }

  return true;
}

template <typename T, typename PMT, typename PVT>
template <typename OPT> inline bool
ODEBND_BASE<T,PMT,PVT>::_SET_PM_STA
( const OPT&options, const unsigned iRHS, const unsigned iQUAD )
{
  if( _vRHS.size() <= iRHS ) return false;
  if( _nq && _vQUAD.size() <= iQUAD ) return false;

  _pRHS   = _vRHS.at( iRHS );
  _opRHS  = _pDAG->subgraph( _nx, _pRHS );
  _pQUAD  = _nq? _vQUAD.at( iQUAD ): 0;
  _opQUAD = _pDAG->subgraph( _nq, _pQUAD );
  const unsigned opmax = _opRHS.size()>_opQUAD.size()?_opRHS.size():_opQUAD.size();

  delete[] _IRHS;   _IRHS = 0;
  delete[] _PMRHS;  _PMRHS = 0;
  delete[] _opRHSi; _opRHSi = 0;
  delete[] _IJAC;   _IJAC = 0;
  delete[] _PMJAC;  _PMJAC = 0;
  delete[] _pJAC;   _pJAC = 0; _opJAC.clear();

  _PMRHS = new PVT[ opmax ];

  switch( options.WRAPMIT){
  case OPT::NONE: default:
    break;

  case OPT::DINEQ:
    _opRHSi = new std::list<const FFOp*>[_nx];
    for( unsigned ix=0; ix<_nx; ix++ )
      _opRHSi[ix] = _pDAG->subgraph( 1, _pRHS+ix );
    break;

  case OPT::ELLIPS:
    _pJAC = _pDAG->FAD( _nx, _pRHS, _nx, _pVAR );
    _opJAC = _pDAG->subgraph( _nx*_nx, _pJAC );
    if( !options.ORDMIT )
      _IJAC = new T[ _opJAC.size() ];
    else if( _PMenv->nord() >= _MVXPenv->nord() )
      _PMJAC = new PVT[ _opJAC.size() ];
    break;
  }

  return true;
}

template <typename T, typename PMT, typename PVT>
template <typename ODEBND, typename ODESLV> inline bool
ODEBND_BASE<T,PMT,PVT>::_hausdorff
( const unsigned ns, const double*tk, const T*Ip, double**Hxk,
  ODEBND&dineq, ODESLV&traj, const unsigned nsamp, std::ostream&os )
{
  int DISPLAY_ODEBND = dineq.options.DISPLAY;
  int DISPLAY_ODESLV = traj.options.DISPLAY;
  dineq.options.DISPLAY = traj.options.DISPLAY = 0;

  // Compute exact bounds 
  T** Ixk0 = new T*[ns+1];
  for( unsigned is=0; is<ns+1; is++ ) Ixk0[is] = new T[_nx];
  if( traj.bounds( ns, tk, Ip, Ixk0, 0, 0, nsamp, os ) != ODEBND::NORMAL ){
    for( unsigned is=0; is<ns+1; is++ ) delete[] Ixk0[is];
    delete[] Ixk0;
    return false;
  }

  // Compute approximate bounds
  T** Ixk = new T*[ns+1];
  for( unsigned is=0; is<ns+1; is++ ) Ixk[is] = new T[_nx];
  try{ dineq.bounds( ns, tk, Ip, Ixk, 0, 0, os ); }
  catch(...){;}
  unsigned nsf = dineq.final_stage();

  dineq.options.DISPLAY = DISPLAY_ODEBND;
  traj.options.DISPLAY = DISPLAY_ODESLV;
  for( unsigned is=0; is<ns+1; is++ ){
    for( unsigned ix=0; ix<_nx; ix++ )
      Hxk[is][ix] = is<=nsf? _dH( Ixk[is][ix], Ixk0[is][ix] ): 0./0.;
    if( dineq.options.DISPLAY >= 1 )
      _print_interm( tk[is], _nx, Hxk[is], "dH", os );
  }

  for( unsigned is=0; is<ns+1; is++ ) delete[] Ixk0[is];
  delete[] Ixk0;
  for( unsigned is=0; is<ns+1; is++ ) delete[] Ixk[is];
  delete[] Ixk;

  return true;
}

template <typename T, typename PMT, typename PVT>
template <typename ODEBND, typename ODESLV> inline bool
ODEBND_BASE<T,PMT,PVT>::_hausdorff
( const unsigned ns, const double*tk, const PVT*PMp, double**Hxk,
  ODEBND&dineq, ODESLV&traj, const unsigned nsamp, std::ostream&os )
{
  int DISPLAY_ODEBND = dineq.options.DISPLAY;
  int DISPLAY_ODESLV = traj.options.DISPLAY;
  dineq.options.DISPLAY = traj.options.DISPLAY = 0;

  // Compute approximate bounds
  PVT** PMxk = new PVT*[ns+1];
  for( unsigned is=0; is<ns+1; is++ ) PMxk[is] = new PVT[_nx];
  try{ dineq.bounds( ns, tk, PMp, PMxk, 0, 0, os ); }
  catch(...){;}
  unsigned nsf = dineq.final_stage();

  // Compute remainder bounds 
  T* Ip = new T[_npar];
  for( unsigned ip=0; ip<_npar; ip++ ) Ip[ip] = PMp[ip].B();
  T** Rxk = new T*[ns+1];
  for( unsigned is=0; is<ns+1; is++ ) Rxk[is] = new T[_nx];
  _remainders( nsf, tk, _npar, Ip, PMxk, Rxk, traj, nsamp, os );

  dineq.options.DISPLAY = DISPLAY_ODEBND;
  traj.options.DISPLAY = DISPLAY_ODESLV;
  for( unsigned is=0; is<ns+1; is++ ){
    for( unsigned ix=0; ix<_nx; ix++ )
      Hxk[is][ix] = is<=nsf? _dH( PMxk[is][ix].R(), Rxk[is][ix] ): 0./0.;
    if( dineq.options.DISPLAY >= 1 ){
      _print_interm( tk[is], _nx, Hxk[is], "dH", os );
    }
  }

  for( unsigned is=0; is<ns+1; is++ ) delete[] Rxk[is];
  delete[] Rxk;
  for( unsigned is=0; is<ns+1; is++ ) delete[] PMxk[is];
  delete[] PMxk;
  delete[] Ip;

  return NORMAL;
}

template <typename T, typename PMT, typename PVT>
template<typename ODESLV> inline bool
ODEBND_BASE<T,PMT,PVT>::_remainders
( const unsigned ns, const double*tk, const unsigned np, const T*Ip,
  const PVT*const*PMxk, T**Rxk, ODESLV&traj, const unsigned nsamp,
  std::ostream&os )
{
   // Initialization of sampled bounds at parameter lower bound
  double *p = new double[np];
  for( unsigned ip=0; ip<np; ip++ )
    p[ip] = Op<T>::l(Ip[ip]);
  double **xk = new double*[ns+1];
  for( unsigned is=0; is<=ns; is++ )
    xk[is] = new double[_nx];
  typename ODESLV::STATUS flag = traj.states( ns, tk, p, xk, 0, 0, os );
  if( flag != ODESLV::NORMAL || nsamp <= 1 ){
    delete[] p;
    for( unsigned is=0; is<=ns; is++ ) delete[] xk[is]; delete[] xk;
    return( flag==ODESLV::NORMAL? true: false );
  }   
  for( unsigned is=0; is<=ns; is++ )
    for( unsigned ix=0; ix<_nx; ix++ )
      Rxk[is][ix] = xk[is][ix] - PMxk[is][ix].polynomial( p );
  
  // Start sampling process
  unsigned* vsamp = new unsigned[np];
  bool flag2 = _remainders( ns, tk, np, Ip, PMxk, Rxk, traj, nsamp, vsamp,
                            0, p, xk, os );
  
  // Clean-up
  delete[] p;
  for( unsigned is=0; is<=ns; is++ ) delete[] xk[is]; delete[] xk;
  delete[] vsamp;
  
  return flag2;
}

template <typename T, typename PMT, typename PVT>
template<typename ODESLV> inline bool
ODEBND_BASE<T,PMT,PVT>::_remainders
( const unsigned ns, const double*tk, const unsigned np, const T*Ip,
  const PVT*const*PMxk, T**Rxk, ODESLV&traj, const unsigned nsamp,
  unsigned* vsamp, const unsigned ip, double*p, double**xk,
  std::ostream&os )
{
  typename ODESLV::STATUS flag = ODESLV::NORMAL;

  // Update bounds for all sampling points
  for( unsigned isamp=0; isamp<nsamp; isamp++ ){
    vsamp[ip] = isamp;

    // Continue recursive call
    if( ip+1 < np ){
      if( !_remainders( ns, tk, np, Ip, PMxk, Rxk, traj, nsamp, vsamp,
          ip+1, p, xk, os ) ) return false;
      continue;
    }

    // Update bounds for current point
#ifdef MC__ODEBND_SUNDIALS_DINEQPM_DEBUG
    std::cout << "Sample: ";
#endif
    for( unsigned ip=0; ip<np; ip++ ){
      p[ip] = Op<T>::l( Ip[ip] ) + vsamp[ip]/(nsamp-1.) * Op<T>::diam( Ip[ip] );
#ifdef MC__ODEBND_SUNDIALS_DINEQPM_DEBUG
      std::cout << p[ip] << "  ";
#endif
    }
#ifdef MC__ODEBND_SUNDIALS_DINEQPM_DEBUG
    std::cout << std::endl;
#endif
    flag = traj.states( ns, tk, p, xk, 0, 0, os );
    if( flag != ODESLV::NORMAL ) return false;
    for( unsigned is=0; is<=ns; is++ )
      for( unsigned ix=0; ix<_nx; ix++ )
        Rxk[is][ix] = Op<T>::hull( xk[is][ix]-PMxk[is][ix].polynomial(p),
                                   Rxk[is][ix] );
  }

  return true;
}  

template <typename T, typename PMT, typename PVT>
inline double
ODEBND_BASE<T,PMT,PVT>::_diam
( const unsigned nx, T*X )
{
  double diam = 0.;
  for( unsigned ix=0; X && ix<nx; ix++ )
    diam = std::max( diam, Op<T>::diam( X[ix] ) );
  return diam;
}

template <typename T, typename PMT, typename PVT>
inline double
ODEBND_BASE<T,PMT,PVT>::_diam
( const unsigned nx, PVT*X )
{
  double diam = 0.;
  for( unsigned ix=0; X && ix<nx; ix++ )
    diam = std::max( diam, Op<T>::diam( X[ix].R() ) );
  return diam;
}

template <typename T, typename PMT, typename PVT>
template <typename U> inline double
ODEBND_BASE<T,PMT,PVT>::_dH
( const U&X, const U&Y )
{
  return std::max( std::fabs(Op<U>::l(X)-Op<U>::l(Y)),
                   std::fabs(Op<U>::u(X)-Op<U>::u(Y)) );
}

template <typename T, typename PMT, typename PVT>
template <typename U> inline void
ODEBND_BASE<T,PMT,PVT>::_print_interm
( const double t, const unsigned nx, const U*x, const std::string&var,
  std::ostream&os )
{
  os << " @t = " << std::scientific << std::setprecision(4)
                 << std::left << t << " :" << std::endl;
  _print_interm( nx, x, var, os );
  return;
}

template <typename T, typename PMT, typename PVT>
template <typename U, typename V> inline void
ODEBND_BASE<T,PMT,PVT>::_print_interm
( const double t, const unsigned nx, const U*x, const V&r,
  const std::string&var, std::ostream&os )
{
  os << " @t = " << std::scientific << std::setprecision(4)
                 << std::left << t << " :" << std::endl;
  _print_interm( nx, x, var, os );
  os << " " << "R" << var.c_str() << " =" << r << std::endl;
  return;
}

template <typename T, typename PMT, typename PVT>
template <typename U> inline void
ODEBND_BASE<T,PMT,PVT>::_print_interm
( const unsigned nx, const U*x, const std::string&var, std::ostream&os )
{
  if( !x ) return;
  for( unsigned ix=0; ix<nx; ix++ )
    os << " " << var.c_str() << "[" << ix << "] = " << x[ix] << std::endl;
  return;
}

template <typename T, typename PMT, typename PVT>
template <typename VRES> inline void
ODEBND_BASE<T,PMT,PVT>::_record
( std::ofstream&ofile, const VRES&bnd, const unsigned iprec )
{
  if( !ofile ) return;

  // Specify format
  ofile << std::right << std::scientific << std::setprecision(iprec);

  // Record computed interval bounds at stage times
  auto it = bnd.begin();
  for( ; it != bnd.end(); ++it ){
    ofile << std::setw(iprec+9) << it->t;
    for( unsigned ix=0; ix<it->nx; ix++ )
      ofile << std::setw(iprec+9) << mc::Op<T>::l( it->X[ix] )
            << std::setw(iprec+9) << mc::Op<T>::u( it->X[ix] );
    ofile << std::endl;
  }
}

} // end namescape mc

#endif

