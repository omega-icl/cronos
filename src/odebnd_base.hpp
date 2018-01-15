// Copyright (C) 2015 Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

#ifndef MC__ODEBND_BASE_HPP
#define MC__ODEBND_BASE_HPP

#undef  MC__ODEBND_BASE_DINEQI_DEBUG
#undef  MC__ODEBND_BASE_DINEQPM_DEBUG

#include <stdexcept>
#include <cassert>
#include <cmath>
#include <fstream>
#include <vector>
#include <list>
#include <sys/time.h>

#include "ellipsoid.hpp"
#include "base_de.hpp"

// *** TO DO
// - Detect block structure and use it in ellipsoidal approach
// - Implement scaling in ellipsoidal approach -> DONE

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

  //! @brief Integration results at a given time instant
  struct Results
  {
    //! @brief Constructors
    Results
      ( const double tk, const unsigned nx1, const double*x1,
        const unsigned nx2=0, const double*x2=0 ):
      t( tk ), nx( nx1+nx2 )
      { X = new T[nx];
        for( unsigned ix=0; ix<nx; ix++ ) X[ix] = ix<nx1? x1[ix]: x2[ix-nx1]; }
    Results
      ( const double tk, const unsigned nx1, const T*Ix1,
        const unsigned nx2=0, const T*Ix2=0 ):
      t( tk ), nx( nx1+nx2 )
      { X = new T[nx];
        for( unsigned ix=0; ix<nx; ix++ ) X[ix] = ix<nx1? Ix1[ix]: Ix2[ix-nx1]; }
    Results
      ( const double tk, const unsigned nx1, const PVT*PMx1,
        const unsigned nx2=0, const PVT*PMx2=0 ):
      t( tk ), nx( nx1+nx2 )
      { X = new T[nx];
        for( unsigned ix=0; ix<nx; ix++ ) X[ix] = ix<nx1? PMx1[ix].B(): PMx2[ix-nx1].B(); }
    Results
      ( const Results&res ):
      t( res.t ), nx( res.nx )
      { X = new T[nx];
        for( unsigned ix=0; ix<nx; ix++ ) X[ix] = res.X[ix]; }
    //! @brief Destructor
    ~Results()
      { delete[] X; }
    //! @brief Time point
    double t;
    //! @brief Solution dimension
    unsigned nx;
    //! @brief Solution bounds
    T* X;
  };

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

  //! @brief storage vector for DAG evaluation in T arithmetic
  std::vector<T> _IWK;

  //! @brief storage vector for DAG evaluation in PM arithmetic
  std::vector<PVT> _PMWK;

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

  //! @brief const pointer to state function in current stage of ODE system
  const FFVar* _pFCT;

  //! @brief number of effective parameters (possibly larger than _np due to lifting)
  unsigned _npar;

  //! @brief number of variables for DAG evaluation
  unsigned _nVAR;

  //! @brief number of variables for DAG evaluation (without quadratures)
  unsigned _nVAR0;

  //! @brief pointer to variables for DAG evaluation
  FFVar* _pVAR;

  //! @brief pointer of variable bounds for DAG evaluation
  T *_IVAR;

  //! @brief pointer to state time **DO NOT FREE**
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

  //! @brief function bounds
  T *_If;

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

  //! @brief function polynomial model
  PVT *_PMf;

  //! @brief quadrature PM remainder radius time derivatives
  double *_radRqdot;

  //! @brief parameter reference
  double *_pref;

  //! @brief state reference
  double *_xref;

  //! @brief linear transformation A matrix for state equations
  double *_A;

  //! @brief linear transformation B matrix for state equations
  double *_B;

  //! @brief linear transformation B matrix for quadrature equations
  double *_Bq;

  //! @brief linear transformation B matrix for state functions
  double *_Bfct;

  //! @brief linear transformed state interval bounds
  T *_Ir;

  //! @brief linear transformed state quadrature bounds
  T *_Irq;

  //! @brief linear transformed state function bounds
  T *_Irfct;

  //! @brief RHS Jacobian interval bounds
  T *_Idfdx;

  //! @brief linear transformed state ellipsoidal bounds
  E _Er;

  //! @brief shape matrix (lower triangular) in ellipsoidal bounds
  double *_Q;

  //! @brief state reference time derivatives
  double *_xrefdot;

  //! @brief linear transformation B matrix time derivatives for state equations
  double *_Bdot;

  //! @brief linear transformation B matrix time derivatives for quadrature equations
  double *_Bqdot;

  //! @brief rotated state interval bounds time derivatives for state equations
  T *_Irdot;

  //! @brief rotated state interval bounds time derivatives for quadrature equations
  T *_Irqdot;

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

  //! @brief RHS polynomial model for quadrature equations
  PVT *_MVXPqdot;

  //! @brief RHS polynomial model for state functions
  PVT *_MVXPfct;

  //! @brief polynomial models for variables in mean-value theorem
  PVT *_MVXPVAR;

  //! @brief pointer to time polynomial model **DO NOT FREE**
  PVT *_MVXPt;

  //! @brief pointer to state polynomial models **DO NOT FREE**
  PVT *_MVXPx;

  //! @brief pointer to parameter polynomial models **DO NOT FREE**
  PVT *_MVXPp;

  //! @brief pointer to quadrature polynomial models **DO NOT FREE**
  PVT *_MVXPq;

  //! @brief polynomial model environment for linear transformation of intial value
  PMT *_MVPenv;

  //! @brief array of initial state polynomial model for initial value
  PVT *_MVPx;

  //! @brief array of parameter polynomial model for initial value
  PVT *_MVPp;

  //! @brief pointer to time polynomial model for initial value **DO NOT FREE**
  PVT *_MVPt;

  //! @brief absolute tolerance for shape matrix of ellipsoidal remainder
  double _QTOLx;

  //! @brief Function to set IC pointer
  bool _IC_SET
    ( const unsigned iIC );

  //! @brief Function to set RHS and QUAD pointers
  bool _RHS_SET
    ( const unsigned iRHS, const unsigned iQUAD );

  //! @brief Static function converting interval bounds to array
  template <typename REALTYPE> static void _I2vec
    ( const unsigned nx, const double*xL, const double*xU, REALTYPE*vec );

  //! @brief Static function converting interval bounds to array
  template <typename REALTYPE> static void _I2vec
    ( const unsigned nx, const T*Ix, REALTYPE*vec, const bool centered=false );

  //! @brief Static function converting interval bounds to array
  template <typename REALTYPE> static void _I2vec
    ( const unsigned nx, const unsigned np, const double*B, const double*dL,
      const double*dU, REALTYPE*vec );

  //! @brief Static function converting interval bounds to array
  template <typename REALTYPE> static void _I2vec
    ( const unsigned nx, const unsigned np, const double*B, const T*Id,
      REALTYPE*vec, const bool centered=false );

  //! @brief Static function converting ellipsoidal bounds to array
  template <typename REALTYPE> static void _E2vec
    ( const unsigned nx, const unsigned np, const double*xref, const double*Qx,
      const double*Bx, REALTYPE*vec );

  //! @brief Function converting GSL array to interval bounds
  template <typename REALTYPE> static void _vec2I
    ( const REALTYPE*vec, const unsigned nx, T*Ix, const bool centered=false );

  //! @brief Function converting GSL array to interval bounds
  template <typename REALTYPE> static void _vec2I
    ( const REALTYPE*vec, const unsigned nx, const unsigned np, const double*pref,
      const T*Ip, double*B, T*Id, T*Ix, const bool centered=false );

  //! @brief Static function converting GSL array to ellipsoidal bounds
  template <typename REALTYPE> static void _vec2E
    ( const REALTYPE*vec, const unsigned nx, const unsigned np, double*Q,
      E&Ed, T*Id, const double*pref, const T*Ip, double*B, double*xref,
      T*Ix, const bool regPSD=false );

  //! @brief Static function converting rotated states into original coordinates
  template<typename U, typename V> static void _ep2x
    ( const unsigned nx, const unsigned np, const V*d, const double*pref,
      const U*p, const double*B, const double*xref, U*x );

  //! @brief Static function converting ellipsoids into interval bounds
  static void _ep2x
    ( const unsigned nx, const unsigned np, double*Q, E&Ed, T*Id,
      const double*pref, const T*Ip, const double*B, const double*xref, T*Ix,
      const bool regPSD );
      
  //! @brief Function to initialize state interval bounding
  template <typename OPT> bool _INI_I_STA
    ( const OPT&options, const unsigned np, const T*Ip, const double QTOL );

  //! @brief Function to retreive state bounds
  template <typename REALTYPE, typename OPT> void _GET_I_STA
    ( const OPT &options, const REALTYPE*x, const REALTYPE*q );

  //! @brief Function to set IC pointer
  template <typename OPT> bool _IC_I_SET
    ( const OPT &options );

  //! @brief Function to initialize quarature interval bounds
  template <typename REALTYPE, typename OPT> bool _IC_I_QUAD
    ( const OPT&options, REALTYPE*vec );

  //! @brief Function to initialize state interval bounds
  template <typename REALTYPE, typename OPT> bool _IC_I_STA
    ( const OPT&options, const double t, REALTYPE*vec );

  //! @brief Function to initialize state interval bounds w/ ellipsoidal bounds 
  static void _IC_I_ELL
    ( const unsigned nx, PVT*MVPx, double*xref, double*Q,
      const unsigned np, double*B, T*Ix, T*Ir, E&Er );

  //! @brief Function to set CC pointer
  template <typename OPT> bool _CC_I_SET
    ( const OPT &options, const unsigned iIC );

  //! @brief Function to reinitialize state bounds at intermediate time
  template <typename REALTYPE, typename OPT> bool _CC_I_STA
    ( const OPT&options, const double t, REALTYPE*vec );

  //! @brief Function to reinitialize state bounds w/ ellipsoidal bounder 
  static void _CC_I_ELL
    ( const unsigned nx, PVT*MVXPf, T*If, const E&Ed, double*Af,
      const unsigned np, double*reff, double*Bf, T*Idf, double*Qf,
      const double QTOL, const double EPS );

  //! @brief Function to set RHS and QUAD pointers
  template <typename OPT> bool _RHS_I_SET
    ( const OPT&options, const unsigned iRHS, const unsigned iQUAD );

  //! @brief Function to calculate the RHS of auxiliary ODEs in interval arithmetic
  template <typename REALTYPE, typename OPT> bool _RHS_I_STA
    ( const OPT&options, double t, const REALTYPE*x, REALTYPE*xdot );

  //! @brief Static function to calculate the RHS of auxiliary ODEs in interval arithmetic w/ ellipsoidal contractor
  static void _RHS_I_ELL
    ( const unsigned nx, PVT*MVXPf, const double*Q, double*A,
      const unsigned np, double*xrefdot, double*Bdot, T*Iddot, double*Qdot,
      const double QTOL, const double EPS, const double QSCALE, const T*W=0 );

  //! @brief Static function to calculate the RHS of auxiliary ODEs w/ ellipsoidal contractor
  template <typename U>
  static bool _RHS_ELL
    ( const unsigned nx, const double*Qr, const double*Ar, const T*Irdot,
      double*Qrdot, const double QTOL, const double EPS, const double QSCALE,
      const U*W );

  //! @brief Function to calculate the RHS of auxiliary ODEs in interval arithmetic
  template <typename REALTYPE, typename OPT> bool _RHS_I_QUAD
    ( const OPT&options, double t, const REALTYPE*x, REALTYPE*qdot );

  //! @brief Static function to calculate the RHS of quadratures in interval arithmetic w/ ellipsoidal contractor
  static void _QUAD_I_ELL
    ( const unsigned nq, const unsigned np, const unsigned offset, PVT*MVYXPg,
      double*Bqdot, T*Idqdot, T*Iqdot );

  //! @brief Function to calculate the Jacobian of auxiliary ODEs in interval arithmetic
  template <typename REALTYPE, typename OPT> bool _JAC_I_STA
    ( const OPT&options, double t, const REALTYPE*x, REALTYPE*jac, REALTYPE*xdot );

  //! @brief Function to calculate the functions at intermediate/end point
  bool _FCT_I_STA
    ( const unsigned pos_fct, const double t );

  //! @brief Function to calculate the functions at intermediate/end point
  template <typename REALTYPE, typename OPT> bool _FCT_I_STA
    ( const OPT&options, const unsigned pos_fct, const double t, REALTYPE*fct );

  //! @brief Function converting polynomial model with interval remainder to GSL array
  template <typename REALTYPE> static void _PMI2vec
    ( const PMT*PMenv, const unsigned nx, const PVT*PMx, const std::pair<double*,
      double*>&Rx, REALTYPE*vec );

  //! @brief Function converting polynomial model with interval remainder to GSL array
  template <typename REALTYPE> static void _PMI2vec
    ( const PMT*PMenv, const unsigned nx, const PVT*PMx, const T*IRx,
      REALTYPE*vec );

  //! @brief Function converting polynomial model with interval remainder to GSL array
  template <typename REALTYPE> static void _PMI2vec
    ( const PMT*PMenv, const unsigned nx, const PVT*PMx, REALTYPE*vec,
      const bool centered=false, const bool neg=false );

  //! @brief Function converting polynomial model with ellipsoidal remainder to GSL array
  template <typename REALTYPE> static void _PME2vec
    ( const PMT*PMenv, const unsigned nx, const PVT*PMx, const double*Qx,
      REALTYPE*vec );

  //! @brief Function converting GSL array to polynomial model with interval remainder
  template <typename REALTYPE> static void _vec2PMI
    ( const REALTYPE*vec, PMT*PMenv, const unsigned nx, PVT*PMx,
      const bool centered=false );

  //! @brief Function converting GSL array to polynomial model with ellipsoidal remainder
  template <typename REALTYPE> static bool _vec2PME
    ( const REALTYPE*vec, PMT*PMenv, const unsigned nx, PVT*PMx, double*Qr,
      E&Er, T*Ir, const bool regPSD=false );

  //! @brief Function converting rotated interval remainder bound back into original coordinates
  template <typename U> static void _e2x
    ( const unsigned nx, const U*d, U*x, const bool reinit=true );

  //! @brief Function converting rotated interval remainder bound back into original coordinates
  template <typename U> static void _e2x
    ( const unsigned nx, const T*w, const U*d, U*x, const bool reinit=true );

  //! @brief Function converting rotated interval remainder bound back into original coordinates
  template <typename U> static void _e2x
    ( const unsigned nx, const CPPL::dsymatrix&w, const U*d, U*x,
      const bool reinit=true );

  //! @brief Function converting ellipsoidal remainder bound into interval remainder bound
  static bool _e2x
    ( const unsigned nx, double*Qr, E&Er, T*Ir, const bool regPSD );

  //! @brief Function to set-valued integration for state polynomial models
  bool _INI_PM_STA
    ( const unsigned np, const PVT*PMp );

  //! @brief Function to set-valued integration for state polynomial models
  template <typename OPT> bool _INI_PM_STA
    ( const OPT&options, const unsigned np, const PVT*PMp, const double QTOL );
      
  //! @brief Function to retreive state polynomial models
  template <typename REALTYPE, typename OPT> void _GET_PM_STA
    ( const OPT&options, const REALTYPE*x, const REALTYPE*q );

  //! @brief Function to set state polynomial models at initial time
  template <typename OPT> bool _IC_PM_SET
    ( const OPT &options );

  //! @brief Function to initialize quarature polynomial models
  void _IC_PM_QUAD
    ();

  //! @brief Function to initialize quarature polynomial models
  template <typename REALTYPE, typename OPT> bool _IC_PM_QUAD
    ( const OPT&options, REALTYPE*vec );

  //! @brief Function to initialize state polynomial models
  void _IC_PM_STA
    ( const double t );

  //! @brief Function to initialize state polynomial models
  template <typename REALTYPE, typename OPT> bool _IC_PM_STA
    ( const OPT&options, const double t, REALTYPE*vec );

  //! @brief Function to initialize state polynomial model w/ ellipsoidal remainder
  static void _IC_PM_ELL
    ( const unsigned nx, PVT*PMx, double*Qr, E&Er, T*Ir, double Qtol );

  //! @brief Function to reset state polynomial models at intermediate time
  template <typename OPT> bool _CC_PM_SET
    ( const OPT &options, const unsigned iIC );

  //! @brief Function to reinitialize state polynomial models at intermediate time
  void _CC_PM_STA
    ( const double t );

  //! @brief Function to reinitialize state polynomial models at intermediate time
  template <typename REALTYPE, typename OPT> bool _CC_PM_STA
    ( const OPT&options, const double t, REALTYPE*vec );

  //! @brief Function to reinitialize state polynomial model w/ ellipsoidal remainder
  static void _CC_PM_ELL
    ( const unsigned nx, const E&Exr, const double*Afr, const T*Ifr,
      PVT*PMf, double*Qfr, const double QTOL, const double EPS );

  //! @brief Function to set RHS and QUAD pointers
  template <typename OPT> bool _RHS_PM_SET
    ( const OPT&options, const unsigned iRHS, const unsigned iQUAD );

  //! @brief Function to calculate the RHS of auxiliary ODEs in polynomial mode arithmetic
  template <typename REALTYPE, typename OPT> int _RHS_PM_STA
    ( const OPT&options, double t, const REALTYPE*x, REALTYPE*xdot );

  //! @brief Static function to calculate the RHS of auxiliary ODEs in polynomial model arithmetic w/ ellipsoidal contractor - approximation using mean-value theorem and interval analysis
  static void _RHS_PM_ELL0
    ( const unsigned nx, PVT*PMf, T*Idfdx, const T*Ir, double*Ar, T*Irdot );

  //! @brief Static function to calculate the RHS of auxiliary ODEs in polynomial model arithmetic w/ ellipsoidal contractor - approximation using mean-value theorem and PM arithmetic
  static void _RHS_PM_ELL1
    ( const unsigned nx, PVT*PMf, PVT*MVXPdfdx, const T*Ir, double*Ar, T*Irdot );

  //! @brief Static function to calculate the RHS of auxiliary ODEs in polynomial model arithmetic w/ ellipsoidal contractor - joint polynomial model in states and parameters
  static void _RHS_PM_ELL2
    ( const unsigned nx, PMT*PMenv, PVT*PMf, PVT*MVXPf, const unsigned np,
      const T*Ir, double*Ar, T*Irdot );

  //! @brief Function to calculate the RHS of auxiliary ODEs in polynomial mode arithmetic
  template <typename REALTYPE, typename OPT> bool _RHS_PM_QUAD
    ( const OPT&options, double t, const REALTYPE*x, REALTYPE*qdot,
      const bool bndinit=true );

  //! @brief Static function to calculate the quadratures in polynomial model arithmetic
  static void _QUAD_PM
    ( FFGraph*DAG, std::list<const FFOp*>&opQUAD, std::vector<PVT>&PMQUAD,
      const unsigned nq, const FFVar*pQUAD, const unsigned nVAR,
      const FFVar*pVAR, PVT*PMVAR, PVT*PMqdot, const bool append=false );

  //! @brief Function to calculate the Jacobian of auxiliary ODEs in polynomial mode arithmetic
  template <typename REALTYPE, typename OPT> bool _JAC_PM_STA
    ( const OPT&options, double t, const REALTYPE*x, REALTYPE*jac, REALTYPE*xdot );

  //! @brief Function to calculate the functions at intermediate/end point
  bool _FCT_PM_STA
    ( const unsigned pos_fct, const double t );

  //! @brief Computes inner bounds approximation using parameter sampling
  template <typename ODESLV> inline bool _bounds
    ( const T*Ip, T**Ixk, T*If, ODESLV&traj, const unsigned nsamp,
      std::vector<Results>&results, std::ostream&os );

  //! @brief Recursive function computing bounds on solutions of IVP in ODEs using sampling
  template <typename ODESLV> inline bool _sampling
    ( const T*Ip, T**Ixk, T*If, ODESLV&traj, const unsigned nsamp,
      unsigned* vsamp, const unsigned ipar, double*p, double**xk,
      double*f, std::vector<Results>&results, std::ostream&os );

  //! @brief Computes Hausdorff distance between interval enclosure and actual reachable set of parametric ODEs, using parameter sampling
  template <typename ODEBND> inline bool _hausdorff
    ( const T*Ip, double**Hxk, double*Hf, ODEBND&dineq, const unsigned nsamp,
      std::ostream&os=std::cout );

  //! @brief Computes Hausdorff distance between polynomial model remainder enclosure and actual remainder function range, using parameter sampling
  template <typename ODEBND, typename ODESLV> inline bool _hausdorff
    ( const PVT*PMp, double**Hxk, double*Hf, ODEBND&dineq, ODESLV&traj,
      const unsigned nsamp, std::ostream&os=std::cout );

  //! @brief Function to bound the remainder function relative to a polynomial model at sampling points
  template<typename ODESLV> bool _remainders
    ( const unsigned ns, const unsigned np, const T*Ip, const PVT*const*PMxk,
      const PVT*PMf, T**Rxk, T*Rf, ODESLV&traj, const unsigned nsamp,
      std::ostream&os=std::cout );

  //! @brief Recrusive function computing bounds on errors between solutions of IVP in ODEs and polynomial approximant using sampling
  template<typename ODESLV> bool _remainders
    ( const unsigned ns, const unsigned np, const T*Ip, const PVT*const*PMxk,
      const PVT*PMf, T**Rxk, T*Rf, ODESLV&traj, const unsigned nsamp,
      unsigned* vsamp, const unsigned ip, double*p, double**xk, double*f,
      std::ostream&os );

  //! @brief Position in symmetric matrix stored in lower triangular form
  static unsigned _ndxLT
    ( const unsigned i, const unsigned j, const unsigned n )
    { return( i<j? _ndxLT(j,i,n): i+j*n-j*(j+1)/2 ); }

  //! @brief Scaling fuction
  template <typename U> static double _scaling
    ( const unsigned ix, const U*W, const double WMAX, const double EPS,
      const double QSCALE );

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
: BASE_DE(), _pRHS(0), _pJAC(0), _pQUAD(0), _pIC(0), _pDIC(0), _pFCT(0),
  _npar(0), _nVAR(0), _nVAR0(0), _pVAR(0)
{
  // Initalize state/parameter arrays
  _opRHSi = 0;
  _IVAR = _It = _Ip = _Ix = _Ixdot =
  _Iq = _Iqdot = _If = 0;
  _PMVAR = _PMt = _PMp = _PMx = _PMxdot =
  _PMq = _PMqdot = _PMf = 0;
  _xLdot = _xUdot = _RxLdot = _RxUdot = _radRqdot = 0;
  _PMenv = 0;

  // Initialize parameterization arrays
  _pref = _xref = _xrefdot = 0;
  _A = _B = _Bq = _Bdot = _Bqdot = _Bfct = _Q = _Qdot = 0;
  _Ir = _Irq = _Idfdx = _Irdot = _Irqdot = _Irfct = 0;

  // Initialize polynomial model environments
  _MVXPenv = _MVPenv = 0;
  _MVXPd = _MVXPf = _MVXPdfdx = _MVXPqdot = _MVXPfct = _MVXPVAR =
  _MVXPt = _MVXPx = _MVXPp = _MVXPq = _MVPx = _MVPp = _MVPt = 0;

  // Initialize ellipsoidal calculus
  E::options.PSDCHK = false;
}

template <typename T, typename PMT, typename PVT>
inline
ODEBND_BASE<T,PMT,PVT>::~ODEBND_BASE
()
{
  /* DO NOT FREE _pRHS, _pQUAD, _pIC, _pFCT */
  delete[] _opRHSi;
  delete[] _pJAC;
  delete[] _pDIC;
  delete[] _pVAR;

  // Free state/quadrature arrays -- Do *NOT* delete _Ip _PMp _PMenv
  delete[] _IVAR;   // **DO NOT DELETE _Ix, _Ip, _Iq**
  delete[] _PMVAR;  // **DO NOT DELETE _PMx, _PMp, _PMq**
  delete[] _PMxdot;
  delete[] _RxLdot;
  delete[] _RxUdot;
  delete[] _PMqdot;
  delete[] _PMf;
  delete[] _radRqdot;
  delete[] _Ixdot;
  delete[] _If;
  delete[] _xLdot;
  delete[] _xUdot;
  delete[] _Iqdot;

  // Free linear transformation arrays
  delete[] _pref;
  delete[] _xref;
  delete[] _A;
  delete[] _B;
  delete[] _Bq;
  delete[] _Ir;
  delete[] _Irq;
  delete[] _Idfdx;
  delete[] _Q;
  delete[] _xrefdot;
  delete[] _Bdot;
  delete[] _Bqdot;
  delete[] _Bfct;
  delete[] _Irdot;
  delete[] _Irqdot;
  delete[] _Irfct;
  delete[] _Qdot;
  delete[] _MVXPd;
  delete[] _MVXPf;
  delete[] _MVXPdfdx;
  delete[] _MVXPqdot;
  delete[] _MVXPfct;
  delete[] _MVXPVAR;  // **DO NOT DELETE _MVXPx, _MVXPp**
  delete   _MVXPenv;
  delete[] _MVPp;
  delete[] _MVPx;
  delete   _MVPenv;
}

template <typename T, typename PMT, typename PVT>
inline bool
ODEBND_BASE<T,PMT,PVT>::_IC_SET
( const unsigned iIC )
{
  if( _vIC.size() <= iIC || _nx0 != _nx ) return false;
  _pIC = _vIC.at( iIC );
  _opIC.clear();
  if( _pIC ) _opIC = _pDAG->subgraph( _nx, _pIC );

  return true;
}

template <typename T, typename PMT, typename PVT>
inline bool
ODEBND_BASE<T,PMT,PVT>::_RHS_SET
( const unsigned iRHS, const unsigned iQUAD )
{
  if( _vRHS.size() <= iRHS ) return false;
  _pRHS = _vRHS.at( iRHS );
  _opRHS.clear();
  if( _pRHS ) _opRHS  = _pDAG->subgraph( _nx, _pRHS );

  if( _nq && _vQUAD.size() <= iQUAD ) return false;
  _pQUAD = _nq? _vQUAD.at( iQUAD ): 0;
  _opQUAD.clear();
  if( _pQUAD ) _opQUAD = _pDAG->subgraph( _nq, _pQUAD );

  return true;
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
template <typename REALTYPE>
inline void
ODEBND_BASE<T,PMT,PVT>::_vec2I
( const REALTYPE*vec, const unsigned nx, const unsigned np, const double*pref,
  const T*Ip, double*B, T*Id, T*Ix, const bool centered )
{
  unsigned ivec = 0;
  for( unsigned iB=0; iB<nx*np; iB++ )
    B[iB] = vec[ivec++];
  for( unsigned ix=0; ix<nx; ix++ ){
    if( !centered ){
      Id[ix] = T( vec[ivec], vec[ivec+1] );
      ivec += 2;
    }
    else{
      Id[ix] = T( -vec[ivec], vec[ivec] );
      ivec++;
    }
  }
  _ep2x( nx, np, Id, pref, Ip, B, 0, Ix );

  return;
}

template <typename T, typename PMT, typename PVT>
template <typename U, typename V> inline void
ODEBND_BASE<T,PMT,PVT>::_ep2x
( const unsigned nx, const unsigned np, const V*d, const double*pref,
  const U*p, const double*B, const double*xref, U*x )
{
  for( unsigned ix=0; ix<nx; ix++ ){   
    x[ix] = 0;
    if( xref ) x[ix] += xref[ix];
    if( d )    x[ix] += d[ix];
    for( unsigned jp=0; jp<np; jp++ )
      x[ix] += ( p[jp] - pref[jp] ) * B[jp*nx+ix];
  }
  return;
}

template <typename T, typename PMT, typename PVT>
inline void
ODEBND_BASE<T,PMT,PVT>::_ep2x
( const unsigned nx, const unsigned np, double*Q, E&Ed, T*Id,
  const double*pref, const T*Ip, const double*B, const double*xref, T*Ix,
  const bool regPSD )
{
  Ed.set( nx, Q );
#ifdef MC__ODEBND_BASE_DINEQI_DEBUG
  std::cout << "Ed =" << Ed << std::endl;
#endif
  if( regPSD && !Ed.psdQ() ){
    double lmin = Ed.eigQ().first(0);
#ifdef MC__ODEBND_BASE_DINEQPM_DEBUG
    std::cout << "lmin =" << lmin << std::endl;
#endif
    for( unsigned ix=0; ix<nx; ix++ )
      Q[_ndxLT(ix,ix,nx)] = Ed.Q(ix,ix) -= lmin;
#ifdef MC__ODEBND_BASE_DINEQPM_DEBUG
    std::cout << "Ed =" << Ed << std::endl;
    int dum; std::cin >> dum;
#endif
  }

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
  T*Ix, const bool regPSD )
{
  unsigned ivec = 0;
  for( unsigned ix=0; ix<nx; ix++ )
    xref[ix] = vec[ivec++];
  for( unsigned iQ=0; iQ<nx*(nx+1)/2; iQ++ )
    Q[iQ] = vec[ivec++];
  for( unsigned iB=0; iB<nx*np; iB++ )
    B[iB] = vec[ivec++];
  return _ep2x( nx, np, Q, Ed, Id, pref, Ip, B, xref, Ix, regPSD );
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
ODEBND_BASE<T,PMT,PVT>::_I2vec
( const unsigned nx, const unsigned np, const double*B, const double*dL,
  const double*dU, REALTYPE*vec )
{
  unsigned ivec = 0;
  for( unsigned iB=0; iB<nx*np; iB++ )
    vec[ivec++] = B[iB];
  for( unsigned ix=0; ix<nx; ix++ ){
    vec[ivec++] = dL[ix];
    vec[ivec++] = dU[ix];
  }
  return;
}

template <typename T, typename PMT, typename PVT>
template <typename REALTYPE>
inline void
ODEBND_BASE<T,PMT,PVT>::_I2vec
( const unsigned nx, const unsigned np, const double*B, const T*Id,
  REALTYPE*vec, const bool centered )
{
  unsigned ivec = 0;
  for( unsigned iB=0; iB<nx*np; iB++ )
    vec[ivec++] = B[iB];
  if( !centered ){
    for( unsigned ix=0; ix<nx; ix++ ){
      vec[ivec++] = Op<T>::l( Id[ix] );
      vec[ivec++] = Op<T>::u( Id[ix] );
    }
  }
  else
    for( unsigned ix=0; ix<nx; ix++ )
      vec[ivec++] = 0.5*Op<T>::diam( Id[ix] );
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
template <typename OPT> inline bool
ODEBND_BASE<T,PMT,PVT>::_INI_I_STA
( const OPT&options, const unsigned np, const T*Ip, const double QTOL )
{
  // Update effective number of parameters
  // (possibly larger than _np if lifting is used)
  _npar = np;

  // Size and set DAG evaluation arrays
  _nVAR0 = _nx+_npar+1;
  _nVAR  = _nVAR0+_nq;
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

  // Reset _MVXPenv and related variables
  unsigned ordmit = options.ORDMIT<0? -options.ORDMIT: options.ORDMIT;
  if( _MVXPenv && ( _MVXPenv->nord() != ordmit 
                 || _MVXPenv->nvar() != _nx+_npar ) ){
    delete[] _MVXPf;   _MVXPf = 0;
    delete[] _MVXPd;   _MVXPd = 0;
    delete[] _MVXPVAR; _MVXPVAR = 0;
    delete   _MVXPenv; _MVXPenv = 0;
  }

  // Reset _MVPenv and related variables
  if( _MVPenv && ( _MVPenv->nord() != 1 
                || _MVPenv->nvar() != _npar ) ){
    delete[] _MVPx;   _MVPx = 0;
    delete[] _MVPp;   _MVPp = 0;
    delete   _MVPenv; _MVPenv = 0;
  }

  // Set parameterization variables
  delete[] _Iqdot; _Iqdot = _nq? new T[_nq]: 0;
  delete[] _If; _If = _nf? new T[_nf]: 0;

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
    delete[] _Bq;      _Bq      = _nq? new double[_nq*_npar]: 0;
    delete[] _Q;       _Q       = new double[_nx*(_nx+1)/2];
    delete[] _Ir;      _Ir      = new T[_nx];
    delete[] _Irq;     _Irq     = _nq? new T[_nq]: 0;
    delete[] _xrefdot; _xrefdot = new double[_nx];
    delete[] _Bdot;    _Bdot    = new double[_nx*_npar];
    delete[] _Bqdot;   _Bqdot   = _nq? new double[_nq*_npar]: 0;
    delete[] _Bfct;    _Bfct    = _nf? new double[_nf*_npar]: 0;
    delete[] _Qdot;    _Qdot    = new double[_nx*(_nx+1)/2];
    delete[] _Irdot;   _Irdot   = new T[_nx];
    delete[] _Irqdot;  _Irqdot  = _nq? new T[_nq]: 0;
    delete[] _Irfct;   _Irfct   = _nf? new T[_nf]: 0;

    delete   _MVXPenv; _MVXPenv = new PMT( _nx+_npar, ordmit );
    _MVXPenv->options = options.PMOPT;
    delete[] _MVXPd;   _MVXPd   = new PVT[_nx];
    delete[] _MVXPf;   _MVXPf   = new PVT[_nx];
    delete[] _MVXPqdot;_MVXPqdot= _nq? new PVT[_nq]: 0;
    delete[] _MVXPfct; _MVXPfct = _nf? new PVT[_nf]: 0;
    delete[] _MVXPVAR; _MVXPVAR = new PVT[_nVAR];
    _MVXPx = _MVXPVAR;
    _MVXPp = _MVXPx + _nx;
    _MVXPt = _MVXPp + _npar;
    _MVXPq = _MVXPt + 1;
    for( unsigned jp=0; jp<_npar; jp++ )
      _MVXPp[jp].set( _MVXPenv, _nx+jp, _Ip[jp] );

    delete   _MVPenv;  _MVPenv  = new PMT( _npar, 1 );
    _MVPenv->options = options.PMOPT;
    delete[] _MVPx;    _MVPx    = new PVT[_nx];
    delete[] _MVPp;    _MVPp    = new PVT[_npar+1];
    _MVPt = _MVPp + _npar;
    for( unsigned ip=0; ip<_npar; ip++ ){
      _pref[ip] = Op<T>::mid( _Ip[ip] );
      _MVPp[ip].set( _MVPenv, ip, _Ip[ip] );
    }
    _QTOLx = QTOL;
    break;
  }

  return true;
}

template <typename T, typename PMT, typename PVT>
template <typename REALTYPE, typename OPT> inline void
ODEBND_BASE<T,PMT,PVT>::_GET_I_STA
( const OPT &options, const REALTYPE*x, const REALTYPE*q )
{
  switch( options.WRAPMIT){
  case OPT::NONE:
  case OPT::DINEQ:
    _vec2I( x, _nx, _Ix );
    if( q ) _vec2I( q, _nq, _Iq );
    break;

  case OPT::ELLIPS:
    default:
    _vec2E( x, _nx, _npar, _Q, _Er, _Ir, _pref, _Ip, _B, _xref, _Ix );
    if( q ) _vec2I( q, _nq, _npar, _pref, _Ip, _Bq, _Irq, _Iq );
    break;
  }
}

template <typename T, typename PMT, typename PVT>
template <typename OPT> inline bool
ODEBND_BASE<T,PMT,PVT>::_IC_I_SET
( const OPT &options )
{
  return _IC_SET( 0 );
}

template <typename T, typename PMT, typename PVT>
template <typename REALTYPE, typename OPT> inline bool
ODEBND_BASE<T,PMT,PVT>::_IC_I_STA
( const OPT&options, const double t, REALTYPE*vec )
{
  switch( options.WRAPMIT){

  case OPT::NONE:
  case OPT::DINEQ:
    *_It = t; // current time
    _pDAG->eval( _opIC, _IWK, _nx, _pIC, _Ix, _npar+1, _pVAR+_nx, _Ip );
    _I2vec( _nx, _Ix, vec );
    break;

  case OPT::ELLIPS:
  default:
    *_MVPt = t; // current time
    _pDAG->eval( _opIC, _PMWK, _nx, _pIC, _MVPx, _npar+1, _pVAR+_nx, _MVPp );
    _IC_I_ELL( _nx, _MVPx, _xref, _Q, _npar, _B, _Ix, _Ir, _Er );
    _E2vec( _nx, _npar, _xref, _Q, _B, vec );
    break;
  }

  return true;
}

template <typename T, typename PMT, typename PVT>
inline void
ODEBND_BASE<T,PMT,PVT>::_IC_I_ELL
( const unsigned nx, PVT*MVPx, double*xref, double*Q,
  const unsigned np, double*B, T*Ix, T*Ir, E&Er )
{
  for( unsigned ix=0; ix<nx; ix++ ){
    Ix[ix] = MVPx[ix].bound();
    for( unsigned jp=0; jp<np; jp++ )
      B[jp*nx+ix] = MVPx[ix].linear( jp, true );
    Ir[ix] = MVPx[ix].bound();
  }
  Er.set( nx, Ir );
  for( unsigned ix=0, iQ=0; ix<nx; ix++ ){
    xref[ix] = Er.c(ix);
    for( unsigned jx=ix; jx<nx; jx++ )
      Q[iQ++] = Er.Q(ix,jx);
  }
#ifdef MC__ODEBND_BASE_DINEQI_DEBUG
  std::cout << "@t0" << std::endl;
  for( unsigned ix=0; ix<nx; ix++ ){
    std::cout << "B[" << ix << ",#] = ";
    for( unsigned jp=0; jp<np; jp++ )
      std::cout << B[jp*nx+ix] << "  ";
    std::cout << std::endl;
  }
  std::cout << "Er =" << Er << std::endl;
  { int dum; std::cin >> dum; }
#endif
}

template <typename T, typename PMT, typename PVT>
template <typename REALTYPE, typename OPT> inline bool
ODEBND_BASE<T,PMT,PVT>::_IC_I_QUAD
( const OPT&options, REALTYPE*vec )
{
  switch( options.WRAPMIT){

  case OPT::NONE:
  case OPT::DINEQ:
    for( unsigned iq=0; iq<_nq; iq++ ) _Iq[iq] = 0.;
    _I2vec( _nq, _Iq, vec );
    break;

  case OPT::ELLIPS:
  default:
    for( unsigned iq=0; iq<_nq; iq++ ){
      _Irq[iq] = _Iq[iq] = 0.;
      for( unsigned jp=0; jp<_npar; jp++ )
        _Bq[jp*_nq+iq] = 0.;
    }
    _I2vec( _nq, _npar, _Bq, _Irq, vec );
    break;
  }

  return true;
}

template <typename T, typename PMT, typename PVT>
template <typename OPT> inline bool
ODEBND_BASE<T,PMT,PVT>::_CC_I_SET
( const OPT &options, const unsigned iIC )
{
  if( !_IC_SET( iIC ) ) return false;
  return true;
}

template <typename T, typename PMT, typename PVT>
template <typename REALTYPE, typename OPT> inline bool
ODEBND_BASE<T,PMT,PVT>::_CC_I_STA
( const OPT&options, const double t, REALTYPE*vec )
{
  switch( options.WRAPMIT){
  case OPT::NONE:
  case OPT::DINEQ:{
    _vec2I( vec, _nx, _Ix ); // current state bounds
    *_It = t; // current time
    _pDAG->eval( _opIC, _IWK, _nx, _pIC, _Ixdot, _nVAR0, _pVAR, _IVAR );
    _I2vec( _nx, _Ixdot, vec );
    break;
   }

  case OPT::ELLIPS:
  default:{
    _vec2E( vec, _nx, _npar, _Q, _Er, _Ir, _pref, _Ip, _B, _xref, _Ix ); // current state enclosure
    for( unsigned jx=0; jx<_nx; jx++ )
      _MVXPd[jx].set( _MVXPenv, jx, _Ir[jx] );
    _ep2x( _nx, _npar, _MVXPd, _pref, _MVXPp, _B, _xref, _MVXPx );
    *_MVXPt = t; // current time
    _pDAG->eval( _opIC, _PMWK, _nx, _pIC, _MVXPf, _nVAR0, _pVAR, _MVXPVAR );
    _CC_I_ELL( _nx, _MVXPf, _Ix, _Er, _A, _npar, _xrefdot, _Bdot, _Irdot, _Qdot,
               _QTOLx, machprec() );
    _E2vec( _nx, _npar, _xrefdot, _Qdot, _Bdot, vec );
    break;
   }
  }

  return true;
}

template <typename T, typename PMT, typename PVT>
inline void
ODEBND_BASE<T,PMT,PVT>::_CC_I_ELL
( const unsigned nx, PVT*MVXPf, T*If, const E&Ed, double*Af,
  const unsigned np, double*reff, double*Bf, T*Idf, double*Qf,
  const double QTOL, const double EPS )
{
  for( unsigned ix=0; ix<nx; ix++ ){
    If[ix] = MVXPf[ix].bound();
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
template <typename OPT> inline bool
ODEBND_BASE<T,PMT,PVT>::_RHS_I_SET
( const OPT&options, const unsigned iRHS, const unsigned iQUAD )
{
  if( !_RHS_SET( iRHS, iQUAD ) ) return false;

  delete[] _opRHSi; _opRHSi = 0;
  switch( options.WRAPMIT){
  case OPT::DINEQ:
    _opRHSi = new std::list<const FFOp*>[_nx];
    for( unsigned ix=0; ix<_nx; ix++ )
      if( _pRHS ) _opRHSi[ix] = _pDAG->subgraph( 1, _pRHS+ix );
    break;   
  default:
    break;
  }

  delete[] _pJAC;
  _pJAC = 0;
  _opJAC.clear();

  return true;
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
    _pDAG->eval( _opRHS, _IWK, _nx, _pRHS, _Ixdot, _nVAR0, _pVAR, _IVAR );
    _I2vec( _nx, _Ixdot, xdot );
    return true;
   
  case OPT::DINEQ:
    if( !_opRHSi ) return false;
    _vec2I( x, _nx, _Ix );   // set current state bounds
    *_It = t; // set current time
    for( unsigned ix=0; ix<_nx; ix++ ){
      T Ixi = _IVAR[ix];
      for( unsigned up=0; up<2; up++ ){ // separate lower/upper bounding subproblems
        _IVAR[ix] = up? Op<T>::u( Ixi ): Op<T>::l( Ixi );
        _pDAG->eval( _opRHSi[ix], _IWK, 1, _pRHS+ix, _Ixdot+ix, _nVAR0,
                     _pVAR, _IVAR );
        if( up ) _xUdot[ix] = Op<T>::u( _Ixdot[ix] );
        else     _xLdot[ix] = Op<T>::l( _Ixdot[ix] );
      }
      _IVAR[ix] = Ixi;
    }
    _I2vec( _nx, _xLdot, _xUdot, xdot );
    return true;
   
  case OPT::ELLIPS:
  default:{
    _vec2E( x, _nx, _npar, _Q, _Er, _Ir, _pref, _Ip, _B, _xref, _Ix, true ); // current state enclosure
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

    // Setup polynomial model expansion of RHS
    *_MVXPt = t; // set current time
    for( unsigned jx=0; jx<_nx; jx++ )
      _MVXPd[jx].set( _MVXPenv, jx, _Ir[jx] );
    _ep2x( _nx, _npar, _MVXPd, _pref, _MVXPp, _B, _xref, _MVXPx );
#ifdef MC__ODEBND_BASE_DINEQI_DEBUG
    for( unsigned ix=0; ix<_nx; ix++ )
      std::cout << "MVXPx[ " << ix << "] = " << _MVXPx[ix] << std::endl;
    //{ int dum; std::cin >> dum; }
#endif

    // Construct the ellipsoidal remainder derivatives
    _pDAG->eval( _opRHS, _PMWK, _nx, _pRHS, _MVXPf, _nVAR0, _pVAR, _MVXPVAR );
    _RHS_I_ELL( _nx, _MVXPf, _Q, _A, _npar, _xrefdot, _Bdot, _Irdot,
       _Qdot, _QTOLx, machprec(), options.QSCALE, _Ir );
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
ODEBND_BASE<T,PMT,PVT>::_RHS_I_ELL
( const unsigned nx, PVT*MVXPf, const double*Q, double*A,
  const unsigned np, double*xrefdot, double*Bdot, T*Iddot, double*Qdot,
  const double QTOL, const double EPS, const double QSCALE, const T*W )
{
  // Extract time derivatives of constant, linear and remainder parts
  // Set reference and linear block RHS
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
#endif
  }

  _RHS_ELL( nx, Q, A, Iddot, Qdot, QTOL, EPS, QSCALE, W );
}

template <typename T, typename PMT, typename PVT>
template <typename U>
inline double
ODEBND_BASE<T,PMT,PVT>::_scaling
( const unsigned ix, const U*W, const double WMAX, const double EPS,
  const double QSCALE )
{
  if( !W || QSCALE <= 0. || WMAX <= EPS ) return 1.;
  double wi = Op<U>::abs(W[ix]);
  return wi/WMAX + std::sqrt(QSCALE);
  //return wi/WMAX < std::sqrt(QSCALE)? std::sqrt(QSCALE): wi/WMAX;
}

template <typename T, typename PMT, typename PVT>
template <typename U>
inline bool
ODEBND_BASE<T,PMT,PVT>::_RHS_ELL
( const unsigned nx, const double*Qx, const double*Ax, const T*Idxdot,
  double*Qxdot, const double QTOL, const double EPS, const double QSCALE,
  const U*W )
{
#ifdef MC__ODEBND_BASE_DINEQPM_DEBUG
  std::cout << "Ex =" << E(nx,Qx) << std::endl;
#endif

 // Set dynamics of shape matrix
  double trQW = 0., WMAX = 0.;
  for( unsigned ix=0; ix<nx; ix++ )
    if( Op<U>::abs(W[ix]) > WMAX ) WMAX = Op<U>::abs(W[ix]);
  for( unsigned ix=0; ix<nx; ix++ ){
    double sqr_wi = sqr(_scaling(ix,W,WMAX,EPS,QSCALE));
    trQW += ( Qx[_ndxLT(ix,ix,nx)]>EPS? Qx[_ndxLT(ix,ix,nx)]/sqr_wi: EPS );
  }
  double sumkappa = 0.;
  const double srqt_trQW = (trQW>0? std::sqrt( trQW ): 0.) + QTOL;
  for( unsigned ix=0; ix<nx; ix++ ){
    double wi = _scaling(ix,W,WMAX,EPS,QSCALE);
#ifdef MC__ODEBND_BASE_DINEQPM_DEBUG
    std::cout << "kappa[" << ix << "] = "
              << Op<T>::diam( Idxdot[ix] ) / ( 2. * wi * srqt_trQW ) << std::endl;
#endif
    sumkappa += Op<T>::diam( Idxdot[ix] ) / ( 2. * wi * srqt_trQW );
  }

  for( unsigned jx=0; jx<nx; jx++ ){
    for( unsigned ix=jx; ix<nx; ix++ ){
      Qxdot[_ndxLT(ix,jx,nx)] = sumkappa * Qx[_ndxLT(ix,jx,nx)];
      for( unsigned kx=0; kx<nx; kx++ )
        Qxdot[_ndxLT(ix,jx,nx)] += Qx[_ndxLT(ix,kx,nx)] * Ax[jx+kx*nx]
                                 + Ax[ix+kx*nx] * Qx[_ndxLT(kx,jx,nx)];
    }
    double wj = _scaling(jx,W,WMAX,EPS,QSCALE);
    Qxdot[_ndxLT(jx,jx,nx)] += Op<T>::diam( Idxdot[jx] ) / 2. * wj * srqt_trQW;
  }

#ifdef MC__ODEBND_BASE_DINEQPM_DEBUG
  for( unsigned ix=0; ix<nx; ix++ )
    std::cout << "Idxdot[" << ix << "] =" << Idxdot[ix] << std::endl;
  E Exdot( nx, Qxdot );
  std::cout << "Exdot =" << Exdot << std::endl;
  { int dum; std::cin >> dum; }
#endif
  return true;
}

template <typename T, typename PMT, typename PVT>
template <typename REALTYPE, typename OPT> inline bool
ODEBND_BASE<T,PMT,PVT>::_RHS_I_QUAD
( const OPT&options, double t, const REALTYPE*x, REALTYPE*qdot )
{
  if( !_pQUAD ) return false;

  switch( options.WRAPMIT){
  case OPT::NONE:
  case OPT::DINEQ:
    _pDAG->eval( _opQUAD, _IWK, _nq, _pQUAD, _Iqdot, _nVAR0, _pVAR, _IVAR );
    _I2vec( _nq, _Iqdot, qdot );
    break;

  case OPT::ELLIPS:
  default:
    _pDAG->eval( _opQUAD, _PMWK, _nq, _pQUAD, _MVXPqdot, _nVAR0, _pVAR, _MVXPVAR );
    _QUAD_I_ELL( _nq, _npar, _nx, _MVXPqdot, _Bqdot, _Irqdot, _Iqdot );
    _I2vec( _nq, _npar, _Bqdot, _Irqdot, qdot );
    break;
  }

   return true;
}

template <typename T, typename PMT, typename PVT>
inline void
ODEBND_BASE<T,PMT,PVT>::_QUAD_I_ELL
( const unsigned nq, const unsigned np, const unsigned offset, PVT*MVYXPg,
  double*Bqdot, T*Idqdot, T*Iqdot )
{
  // Extract time derivatives of constant, linear and remainder parts
  // Set reference and linear block RHS
  for( unsigned iq=0; iq<nq; iq++ ){
#ifdef MC__ODEBND_BASE_DINEQI_DEBUG
    std::cout << "MVYXPg[" << iq << "] = " << MVYXPg[iq] << std::endl;
#endif
    Iqdot[iq] = MVYXPg[iq].B();
#ifdef MC__ODEBND_BASE_DINEQI_DEBUG
    std::cout << "Iqdot[" << iq << "] = " << Iqdot[iq]
              << " : " << Op<T>::mid(Iqdot[iq]) << std::endl;
#endif
    for( unsigned jp=0; jp<np; jp++ )
      Bqdot[iq+jp*nq] = MVYXPg[iq].linear(offset+jp,true);
#ifdef MC__ODEBND_BASE_DINEQI_DEBUG
    std::cout << "Bqdot[" << iq << ",#] = ";
    for( unsigned jp=0; jp<np; jp++ )
      std::cout << Bqdot[iq+jp*nq] << "  ";
    std::cout << std::endl;
    std::cout << "MVYXPg[" << iq << "] = " << MVYXPg[iq] << std::endl;
#endif
    Idqdot[iq] = MVYXPg[iq].B();
#ifdef MC__ODEBND_BASE_DINEQI_DEBUG
    std::cout << "Idqdot[" << iq << "] = " << Idqdot[iq]
              << " : " << Op<T>::mid(Idqdot[iq]) << std::endl;
#endif
  }
}

template <typename T, typename PMT, typename PVT>
template <typename REALTYPE, typename OPT> inline bool
ODEBND_BASE<T,PMT,PVT>::_JAC_I_STA
( const OPT&options, double t, const REALTYPE*x, REALTYPE*jac, REALTYPE*xdot )
{
  // Jacobian not (yet) implemented
  return false;
}

template <typename T, typename PMT, typename PVT> inline bool
ODEBND_BASE<T,PMT,PVT>::_FCT_I_STA
( const unsigned pos_fct, const double t )
{
  if( !_nf ) return true;

  *_It = t; // set current time
  const FFVar* pFCT = _vFCT.at( pos_fct );
#ifdef MC__ODEBNDS_BASE_DINEQI_DEBUG
    for( unsigned j=0; j<_nq; j++ )
      std::cout << "Iq[" << j << "] = " << _Iq[j] << std::endl;
#endif
  _pDAG->eval( _nf, pFCT, _If, _nVAR, _pVAR, _IVAR, pos_fct?true:false );
  return true;
}

template <typename T, typename PMT, typename PVT>
template <typename REALTYPE, typename OPT> inline bool
ODEBND_BASE<T,PMT,PVT>::_FCT_I_STA
( const OPT&options, const unsigned pos_fct, const double t,
  REALTYPE*fct )
{
  if( !_nf || !fct ) return true;
  _pFCT = _vFCT.at( pos_fct );

  switch( options.WRAPMIT){
  case OPT::NONE:
  case OPT::DINEQ:
    *_It = t; // set current time
    _pDAG->eval( _nf, _pFCT, _If, _nVAR, _pVAR, _IVAR, pos_fct?true:false );
    _I2vec( _nf, _If, fct );
    break;

  case OPT::ELLIPS:
  default:
    *_MVXPt = t; // set current time
    for( unsigned jx=0; jx<_nx; jx++ )
      _MVXPd[jx].set( _MVXPenv, jx, _Ir[jx] );
    _ep2x( _nx, _npar, _MVXPd, _pref, _MVXPp, _B, _xref, _MVXPx );
    _ep2x( _nq, _npar, _Irq, _pref, _MVXPp, _Bq, 0, _MVXPq );
#ifdef MC__ODEBND_BASE_DINEQI_DEBUG
    for( unsigned j=0; j<_nq; j++ )
      std::cout << "MVXPq[" << j << "] = " << _MVXPq[j] << std::endl;
#endif
    if( pos_fct ){
      _vec2I( fct, _nf, _npar, _pref, _Ip, _Bfct, _Irfct, _If );
      _ep2x( _nf, _npar, _Irfct, _pref, _MVXPp, _Bfct, 0, _MVXPfct );
      _pDAG->eval( _nf, _pFCT, _MVXPfct, _nVAR, _pVAR, _MVXPVAR, true );
    }
    else
      _pDAG->eval( _nf, _pFCT, _MVXPfct, _nVAR, _pVAR, _MVXPVAR );
#ifdef MC__ODEBND_BASE_DINEQI_DEBUG
    for( unsigned j=0; j<_nf; j++ )
      std::cout << "MVXPfct[" << j << "] = " << _MVXPfct[j] << std::endl;
#endif
    _QUAD_I_ELL( _nf, _npar, _nx, _MVXPfct, _Bfct, _Irfct, _If );
    _I2vec( _nf, _npar, _Bfct, _Irfct, fct );
    break;
  }

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
inline bool
ODEBND_BASE<T,PMT,PVT>::_vec2PME
( const REALTYPE*vec, PMT*PMenv, const unsigned nx, PVT*PMx, double*Qr,
  E&Er, T*Ir, const bool regPSD )
{
  unsigned ivec=0;
  for( unsigned ix=0; ix<nx; ix++  ){
    PMx[ix].set( PMenv );
    PMx[ix].set( vec+ivec );
    ivec += PMenv->nmon();
  }
  for( unsigned iQ=0; iQ<(nx*(nx+1))/2; iQ++ )
    Qr[iQ] = vec[ivec++];
  if( !_e2x( nx, Qr, Er, Ir, regPSD ) ) return false;
  for( unsigned ix=0; ix<nx; ix++  )
    PMx[ix].set( Ir[ix] );
  return true;
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
template <typename U> inline void
ODEBND_BASE<T,PMT,PVT>::_e2x
( const unsigned nx, const T*w, const U*d, U*x, const bool reinit )
{
  for( unsigned ix=0; ix<nx; ix++ ){   
    if( reinit ) x[ix] = 0.;
    x[ix] += d[ix] * Op<T>::diam(w[ix]) / 2.;
  }
  return;
}

template <typename T, typename PMT, typename PVT>
template <typename U> inline void
ODEBND_BASE<T,PMT,PVT>::_e2x
( const unsigned nx, const CPPL::dsymatrix&w, const U*d, U*x,
  const bool reinit )
{
  for( unsigned ix=0; ix<nx; ix++ ){   
    if( reinit ) x[ix] = 0.;
    for( unsigned jx=0; jx<nx; jx++ )
      x[ix] += w(ix,jx) * d[jx];
  }
  return;
}

template <typename T, typename PMT, typename PVT>
inline bool
ODEBND_BASE<T,PMT,PVT>::_e2x
( const unsigned nx, double*Qr, E&Er, T*Ir, const bool regPSD )
{
  Er.set( nx, Qr );
#ifdef MC__ODEBND_BASE_DINEQPM_DEBUG
  std::cout << "Er =" << Er << std::endl;
#endif
/*
  if( regPSD && !Er.psdQ() ){
//#ifdef MC__ODEBND_BASE_DINEQPM_DEBUG
    std::cout << "Er not a p.s.d. matrix; lmin =" << Er.eigQ().first(0) << std::endl;
//#endif
    return false;
  }
*/
  if( regPSD && !Er.psdQ() ){
    auto eigdec = Er.eigQ();
#ifdef MC__ODEBND_BASE_DINEQPM_DEBUG
    std::cout << "Er not a p.s.d. matrix; lmin =" << eigdec.first(0) << std::endl;
#endif
    for( unsigned i=0; eigdec.first(i)<0 && i<nx; i++ ){
      auto Qcor = eigdec.first(i) * ( eigdec.second.col(i) * t( eigdec.second.col(i) ) );
      for( unsigned ix=0; ix<nx; ix++ )
        for( unsigned jx=0; jx<=ix; jx++ )
          Qr[_ndxLT(ix,jx,nx)] = Er.Q(ix,jx) -= Qcor(ix,jx);
#ifdef MC__ODEBND_BASE_DINEQPM_DEBUG
      Er.reset_aux();
      std::cout << "Er =" << Er << std::endl;
      std::cout << "lmin =" << Er.eigQ().first(0) << std::endl;
      int dum; std::cin >> dum;
#endif
    }
/*
    double lmin = Er.eigQ().first(0);
//#ifdef MC__ODEBND_BASE_DINEQPM_DEBUG
    std::cout << "Er not a p.s.d. matrix; lmin =" << lmin << std::endl;
//#endif
    for( unsigned ix=0; ix<nx; ix++ )
      Qr[_ndxLT(ix,ix,nx)] = Er.Q(ix,ix) -= lmin;
#ifdef MC__ODEBND_BASE_DINEQPM_DEBUG
    std::cout << "Er =" << Er << std::endl;
    int dum; std::cin >> dum;
#endif
*/
  }

  for( unsigned ix=0; ix<nx; ix++ )
    Ir[ix] = T( Er.l(ix), Er.u(ix) );

#ifdef MC__ODEBND_BASE_DINEQPM_DEBUG
  for( unsigned ix=0; ix<nx; ix++ )
    std::cout << "Ir[" << ix << "] = " << Ir[ix] << std::endl;
  { int dum; std::cin >> dum; }
#endif
  return true;
}

template <typename T, typename PMT, typename PVT>
template <typename REALTYPE>
inline void
ODEBND_BASE<T,PMT,PVT>::_PMI2vec
( const PMT*PMenv, const unsigned nx, const PVT*PMx, REALTYPE*vec,
  const bool centered, const bool neg )
{
  unsigned ivec=0;
  for( unsigned ix=0; ix<nx; ix++ ){
    //if( centered ) PMx[ix].center();
    std::pair<unsigned, const double*> PMcoef = PMx[ix].coefmon();
    unsigned imon = 0;
    for( ; imon<PMcoef.first; imon++ )
      vec[ivec++] = PMcoef.second[imon];
    for( ; imon<PMenv->nmon(); imon++ )
      vec[ivec++] = 0.;
    if( !centered ){
      vec[ivec++] = neg? Op<T>::u( PMx[ix].remainder() ): Op<T>::l( PMx[ix].remainder() );
      vec[ivec++] = neg? Op<T>::l( PMx[ix].remainder() ): Op<T>::u( PMx[ix].remainder() );
    }
    else{
      vec[ivec++] = neg? -0.5*Op<T>::diam( PMx[ix].remainder() ): 0.5*Op<T>::diam( PMx[ix].remainder() );
    }
  }
  return;
}

template <typename T, typename PMT, typename PVT>
template <typename REALTYPE>
inline void
ODEBND_BASE<T,PMT,PVT>::_PMI2vec
( const PMT*PMenv, const unsigned nx, const PVT*PMx, const std::pair<double*,
  double*>&Rx, REALTYPE*vec )
{
  unsigned ivec=0;
  for( unsigned ix=0; ix<nx; ix++ ){
    std::pair<unsigned, const double*> PMcoef = PMx[ix].coefmon();
    unsigned imon = 0;
    for( ; imon<PMcoef.first; imon++ )
      vec[ivec++] = PMcoef.second[imon];
    for( ; imon<PMenv->nmon(); imon++ )
      vec[ivec++] = 0.;
    vec[ivec++] = ( Rx.first?  Rx.first[ix]:  0. );
    vec[ivec++] = ( Rx.second? Rx.second[ix]: 0. );
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
inline bool
ODEBND_BASE<T,PMT,PVT>::_INI_PM_STA
( const unsigned np, const PVT* PMp )
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
  _nVAR0 = _nx+_npar+1;
  _nVAR  = _nVAR0+_nq;
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


  return true;
}

template <typename T, typename PMT, typename PVT>
template <typename OPT> inline bool
ODEBND_BASE<T,PMT,PVT>::_INI_PM_STA
( const OPT&options, const unsigned np, const PVT* PMp, const double QTOL )
{
  if( !_INI_PM_STA( np, PMp ) ) return false;

  // Reset _MVXPVenv and related variables
  unsigned ordmit = options.ORDMIT<0? -options.ORDMIT: options.ORDMIT;
  unsigned MVXPsize = ( options.ORDMIT<(int)_PMenv->nord()? ordmit: _PMenv->nord() ); 
  unsigned MVXPdim  = ( options.ORDMIT<0? _npar: _nx+_npar  );
  if( _MVXPenv && ( _MVXPenv->nord() != MVXPsize || _MVXPenv->nvar() != MVXPdim ) ){
    delete[] _MVXPf;    _MVXPf = 0;
    delete[] _MVXPdfdx; _MVXPdfdx = 0;
    delete[] _MVXPd;    _MVXPd = 0;
    delete[] _MVXPVAR;  _MVXPVAR = _MVXPx = _MVXPp = 0; 
    delete   _MVXPenv;  _MVXPenv = 0;
  }

  // Size state/quadrature derivative arrays
  delete[] _PMxdot;   _PMxdot = new PVT[_nx];
  delete[] _PMqdot;   _PMqdot = _nq? new PVT[_nq]: 0;
  delete[] _PMf;      _PMf = _nf? new PVT[_nf]: 0;
  delete[] _radRqdot; _radRqdot  = _nq? new double[_nq]: 0;

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
    delete[] _MVXPVAR;  _MVXPVAR  = new PVT[_nVAR0];
    _MVXPx = _MVXPVAR;
    _MVXPp = _MVXPx + _nx;
    _MVXPt = _MVXPp + _npar;
    for( unsigned ip=0; ip<_npar; ip++ )
      _MVXPp[ip].set( _MVXPenv, ip, _PMp[ip].B() );

    _QTOLx = QTOL;
    break;
  }

  return true;
}

template <typename T, typename PMT, typename PVT>
template <typename REALTYPE, typename OPT> inline void
ODEBND_BASE<T,PMT,PVT>::_GET_PM_STA
( const OPT &options, const REALTYPE*x, const REALTYPE*q )
{
  switch( options.WRAPMIT){
   case OPT::NONE:
    _vec2PMI( x, _PMenv, _nx, _PMx, true );
    break;

   case OPT::DINEQ:
    _vec2PMI( x, _PMenv, _nx, _PMx, false );
    break;

   case OPT::ELLIPS:
   default:
    _vec2PME( x, _PMenv, _nx, _PMx, _Q, _Er, _Ir );
    break;
  }

  if( q ) _vec2PMI( q, _PMenv, _nq, _PMq, true );
}

template <typename T, typename PMT, typename PVT>
template <typename OPT> inline bool
ODEBND_BASE<T,PMT,PVT>::_IC_PM_SET
( const OPT &options )
{
  return _IC_SET( 0 );
}

template <typename T, typename PMT, typename PVT>
inline void
ODEBND_BASE<T,PMT,PVT>::_IC_PM_STA
( const double t )
{
  *_PMt = t; // current time
  _pDAG->eval( _opIC, _PMWK, _nx, _pIC, _PMx, _npar+1, _pVAR+_nx, _PMp );
}

template <typename T, typename PMT, typename PVT>
template <typename REALTYPE, typename OPT> inline bool
ODEBND_BASE<T,PMT,PVT>::_IC_PM_STA
( const OPT&options, const double t, REALTYPE*vec )
{
  _IC_PM_STA( t );

  switch( options.WRAPMIT){
  case OPT::NONE:
    // Whether or not to ignore the remainder
    if( !options.PMNOREM )
      _PMI2vec( _PMenv, _nx, _PMx, vec, true );
    else
      _PMI2vec( _PMenv, _nx, _PMx, 0, vec );
    break;

  case OPT::DINEQ:
    // Whether or not to ignore the remainder
    if( !options.PMNOREM )
      _PMI2vec( _PMenv, _nx, _PMx, vec, false );
    else
      _PMI2vec( _PMenv, _nx, _PMx, 0, vec );
    break;

  case OPT::ELLIPS:
  default:
    _IC_PM_ELL( _nx, _PMx, _Q, _Er, _Ir, _QTOLx );
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
( const unsigned nx, PVT*PMx, double*Qr, E&Er, T*Ir, double QTOL )
{
  double norm1R = 0.;
  for( unsigned ix=0; ix<nx; ix++ )
    norm1R += Op<T>::diam( PMx[ix].remainder() ) / 2.;
  for( unsigned ix=0, iQ=0; ix<nx; ix++ ){
    Qr[iQ++] = norm1R * Op<T>::diam( PMx[ix].remainder() ) / 2. + 1e2*QTOL;
    for( unsigned jx=ix+1; jx<nx; jx++ ) Qr[iQ++] = 0.;
    Ir[ix] = PMx[ix].remainder();
  }
  Er.set( nx, Qr );
#ifdef MC__ODEBND_BASE_DINEQPM_DEBUG
  std::cout << "ERx0 =" << Er << std::endl;
#endif
}

template <typename T, typename PMT, typename PVT>
inline void
ODEBND_BASE<T,PMT,PVT>::_IC_PM_QUAD
()
{
  for( unsigned iq=0; iq<_nq; iq++ ) _PMq[iq] = 0.;
}

template <typename T, typename PMT, typename PVT>
template <typename REALTYPE, typename OPT> inline bool
ODEBND_BASE<T,PMT,PVT>::_IC_PM_QUAD
( const OPT&options, REALTYPE*vec )
{
  if( !_vQUAD.size() || !_nq ) return true;
  _IC_PM_QUAD();
  _PMI2vec( _PMenv, _nq, _PMq, vec, true );
  return true;
}

template <typename T, typename PMT, typename PVT>
template <typename OPT> inline bool
ODEBND_BASE<T,PMT,PVT>::_CC_PM_SET
( const OPT &options, const unsigned iIC )
{
  if( !_IC_SET( iIC ) ) return false;
  delete[] _pDIC;  _pDIC = 0; _opDIC.clear();

  switch( options.WRAPMIT){
  case OPT::NONE:
  case OPT::DINEQ:
    break;

  case OPT::ELLIPS:
  default:
    if( !options.ORDMIT ){
      _pDIC = _pDAG->FAD( _nx, _pIC, _nx, _pVAR );
      _opDIC = _pDAG->subgraph( _nx*_nx, _pDIC );
    }

    else if( options.ORDMIT < 0 || _MVXPenv->nord() <= _PMenv->nord() ){
      _pDIC = _pDAG->FAD( _nx, _pIC, _nx, _pVAR );
      _opDIC = _pDAG->subgraph( _nx*_nx, _pDIC );
    }

    break;
  }

  return true;
}

template <typename T, typename PMT, typename PVT>
inline void
ODEBND_BASE<T,PMT,PVT>::_CC_PM_STA
( const double t )
{
  *_PMt = t; // current time
  _pDAG->eval( _opIC, _PMWK, _nx, _pIC, _PMxdot, _nVAR0, _pVAR, _PMVAR );
  for( unsigned ix=0; ix<_nx; ix++ ) _PMx[ix] = _PMxdot[ix];
}

template <typename T, typename PMT, typename PVT>
template <typename REALTYPE, typename OPT> inline bool
ODEBND_BASE<T,PMT,PVT>::_CC_PM_STA
( const OPT&options, const double t, REALTYPE*x )
{
  switch( options.WRAPMIT){
  case OPT::NONE:
    *_PMt = t; // current time
    _vec2PMI( x, _PMenv, _nx, _PMx, true ); // current state/quadrature polynomial model
    _pDAG->eval( _opIC, _PMWK, _nx, _pIC, _PMxdot, _nVAR0, _pVAR, _PMVAR );
    // Whether or not to ignore the remainder
    if( !options.PMNOREM )
      _PMI2vec( _PMenv, _nx, _PMxdot, x, true );
    else
      _PMI2vec( _PMenv, _nx, _PMxdot, 0, x );
    break;
   
  case OPT::DINEQ:
    *_PMt = t; // current time
    _vec2PMI( x, _PMenv, _nx, _PMx, false ); // current state/quadrature polynomial model
    _pDAG->eval( _opIC, _PMWK, _nx, _pIC, _PMxdot, _nVAR0, _pVAR, _PMVAR );
    // Whether or not to ignore the remainder
    if( !options.PMNOREM )
      _PMI2vec( _PMenv, _nx, _PMxdot, x, false );
    else
      _PMI2vec( _PMenv, _nx, _PMxdot, 0, x );
    break;
   
  case OPT::ELLIPS:
  default:{
    *_PMt = t; // current time   
    _vec2PME( x, _PMenv, _nx, _PMx, _Q, _Er, _Ir, true ); // current state/quadrature polynomial model

    // In this variant a bound on the Jacobian matrix is computed and the
    // linear part is taken as the mid-point of this matrix
    if( !options.ORDMIT ){
      *_It = t; // current time
      for( unsigned ix=0; ix<_nx; ix++ ){
        _Ix[ix] = _PMx[ix].bound(); // set current state bounds
        _PMx[ix].center().set( T(0.) ); // cancel remainder term
      }
      _pDAG->eval( _opIC, _PMWK, _nx, _pIC, _PMxdot, _nVAR0, _pVAR, _PMVAR );
      _pDAG->eval( _opDIC, _IWK, _nx*_nx, _pDIC, _Idfdx, _nVAR0, _pVAR, _IVAR );
      _RHS_PM_ELL0( _nx, _PMxdot, _Idfdx, _Ir, _A, _Irdot );
    }

    // In this variant a polynomial model of the Jacobian matrix is computed and the
    // linear part is taken as the mid-point of this matrix - reduced space
    else if( options.ORDMIT < 0 ){
      *_MVXPt = t; // current time
      for( unsigned jx=0; jx<_nx; jx++ ){
        _MVXPx[jx].set( _MVXPenv ).set( _PMx[jx].center(), true );
        _PMx[jx].set( T(0.) );
      }
      _pDAG->eval( _opIC, _PMWK, _nx, _pIC, _PMxdot, _nVAR0, _pVAR, _PMVAR );
      _pDAG->eval( _opDIC, _PMWK, _nx*_nx, _pDIC, _MVXPdfdx, _nVAR0, _pVAR, _MVXPVAR );
      _RHS_PM_ELL1( _nx, _PMxdot, _MVXPdfdx, _Ir, _A, _Irdot );
    }

    // In this variant a polynomial model of the Jacobian matrix is computed and the
    // linear part is taken as the mid-point of this matrix - full space
    else if( _MVXPenv->nord() <= _PMenv->nord() ){
      *_MVXPt = t; // current time
      for( unsigned jx=0; jx<_nx; jx++ ){
        _MVXPd[jx].set( _MVXPenv, _npar+jx, _Ir[jx] );
        _MVXPx[jx].set( _MVXPenv ).set( _PMx[jx].center().set( T(0.) ), true );
      }
      _e2x( _nx, _MVXPd, _MVXPx, false );
      _pDAG->eval( _opIC, _PMWK, _nx, _pIC, _PMxdot, _nVAR0, _pVAR, _PMVAR );
      _pDAG->eval( _opDIC, _PMWK, _nx*_nx, _pDIC, _MVXPdfdx, _nVAR0, _pVAR, _MVXPVAR );
      _RHS_PM_ELL1( _nx, _PMxdot, _MVXPdfdx, _Ir, _A, _Irdot );
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
      _pDAG->eval( _opIC, _PMWK, _nx, _pIC, _MVXPf, _nVAR0, _pVAR, _MVXPVAR );
      _RHS_PM_ELL2( _nx, _PMenv, _PMxdot, _MVXPf, _npar, _Ir, _A, _Irdot );
    }

    // Whether or not to ignore the remainder
    if( !options.PMNOREM ){
      _CC_PM_ELL( _nx, _Er, _A, _Irdot, _PMxdot, _Qdot, _QTOLx, machprec() );
      _PME2vec( _PMenv, _nx, _PMxdot, _Qdot, x );
    }
    else
      _PME2vec( _PMenv, _nx, _PMxdot, 0, x );
    break;
   }
  }

  return true;
}

template <typename T, typename PMT, typename PVT>
inline void
ODEBND_BASE<T,PMT,PVT>::_CC_PM_ELL
( const unsigned nx, const E&Exr, const double*Afr, const T*Ifr,
  PVT*PMf, double*Qfr, const double QTOL, const double EPS )
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
  for( unsigned jx=0; jx<nx; jx++ ){
    Qfr[_ndxLT(jx,jx,nx)] = Efr.Q(jx,jx) + 1e2*QTOL;
    for( unsigned ix=jx+1; ix<nx; ix++ )
      Qfr[_ndxLT(ix,jx,nx)] = Efr.Q(ix,jx);
    PMf[jx].set( T(Efr.l(jx),Efr.u(jx)) );
  }

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
template <typename REALTYPE, typename OPT> inline int
ODEBND_BASE<T,PMT,PVT>::_RHS_PM_STA
( const OPT&options, double t, const REALTYPE*x, REALTYPE*xdot )
{
  if( !_pRHS ) return -1;

  switch( options.WRAPMIT ){
  case OPT::NONE:
    _vec2PMI( x, _PMenv, _nx, _PMx, true );   // set current state polynomial model
#ifdef MC__ODEBND_BASE_DINEQPM_DEBUG
    _print_interm( t, _nx, _PMx, "PMx Intermediate", std::cerr );
#endif
    *_PMt = t; // set current time
    _pDAG->eval( _opRHS, _PMWK, _nx, _pRHS, _PMxdot, _nVAR0, _pVAR, _PMVAR );
    // Whether or not to ignore the remainder
    if( !options.PMNOREM )
      _PMI2vec( _PMenv, _nx, _PMxdot, xdot, true );
    else
      _PMI2vec( _PMenv, _nx, _PMxdot, 0, xdot );
#ifdef MC__ODEBND_BASE_DINEQPM_DEBUG
    _print_interm( t, _nx, _PMxdot, "PMxdot Intermediate", std::cerr );
    {int dum; std::cin >> dum;}
#endif
    return 0;  

  case OPT::DINEQ:
    _vec2PMI( x, _PMenv, _nx, _PMx, false );   // set current state polynomial model
    for( unsigned ix=0; ix<_nx; ix++ ) _PMx[ix].center();
#ifdef MC__ODEBND_BASE_DINEQPM_DEBUG
    _print_interm( t, _nx, _PMx, "PMx Intermediate", std::cerr );
#endif
    *_PMt = t; // set current time
    for( unsigned ix=0; ix<_nx; ix++ ){
      T Rxi = _PMVAR[ix].remainder();
      for( unsigned up=0; up<2; up++ ){ // separate lower/upper bounding subproblems
        if( up ) _PMVAR[ix].set( Op<T>::u( Rxi ) );
        else     _PMVAR[ix].set( Op<T>::l( Rxi ) );
        _pDAG->eval( _opRHSi[ix], _PMWK, 1, _pRHS+ix, _PMxdot+ix, _nVAR0,
                     _pVAR, _PMVAR );
        if( up ) _RxUdot[ix] = Op<T>::u( _PMxdot[ix].remainder() );
        else     _RxLdot[ix] = Op<T>::l( _PMxdot[ix].remainder() );
      }
      _PMVAR[ix].set( Rxi );
    }
    // Whether or not to ignore the remainder
    if( !options.PMNOREM )
      _PMI2vec( _PMenv, _nx, _PMxdot, std::make_pair(_RxLdot, _RxUdot), xdot );
    else
      _PMI2vec( _PMenv, _nx, _PMxdot, 0, xdot );
#ifdef MC__ODEBND_BASE_DINEQPM_DEBUG
    _print_interm( t, _nx, _PMxdot, "PMxdot Intermediate", std::cerr );
    {int dum; std::cin >> dum;}
#endif
    return 0;  
   
  case OPT::ELLIPS:
  default:{
    if( !_vec2PME( x, _PMenv, _nx, _PMx, _Q, _Er, _Ir, true ) ) return 1; // set current state polynomial model
    *_PMt = t; // set current time   
#ifdef MC__ODEBND_BASE_DINEQPM_DEBUG
    _print_interm( t, _nx, _PMx, _Er, "PMx Intermediate", std::cerr );
    std::cout << "l =" << _Er.eigQ().first(0) << " ; " << _Er.eigQ().first(_nx-1) << std::endl;
#endif

    // In this variant a bound on the Jacobian matrix is computed and the
    // linear part is taken as the mid-point of this matrix
    if( !options.ORDMIT ){
      *_It = t; // set current time
      for( unsigned ix=0; ix<_nx; ix++ ){
        _IVAR[ix] = _PMx[ix].bound(); // set current state bounds
        _PMVAR[ix].center().set( T(0.) ); // cancel remainder term
      }
      _pDAG->eval( _opRHS, _PMWK, _nx, _pRHS, _PMxdot, _nVAR0, _pVAR, _PMVAR );
      _pDAG->eval( _opJAC, _IWK, _nx*_nx, _pJAC, _Idfdx, _nVAR0, _pVAR, _IVAR );
      _RHS_PM_ELL0( _nx, _PMxdot, _Idfdx, _Ir, _A, _Irdot );
    }

    // In this variant a polynomial model of the Jacobian matrix is computed and the
    // linear part is taken as the mid-point of this matrix
    else if( options.ORDMIT < 0 ){
      *_MVXPt = t; // set current time
      for( unsigned jx=0; jx<_nx; jx++ ){
        _MVXPx[jx].set( _MVXPenv ).set( _PMx[jx].center(), true );
        _PMx[jx].set( T(0.) );
      }
#ifdef MC__ODEBND_BASE_DINEQPM_DEBUG
      _print_interm( t, _nx, _MVXPx, "MVXPx Intermediate", std::cerr );
#endif
      _pDAG->eval( _opRHS, _PMWK, _nx, _pRHS, _PMxdot, _nVAR0, _pVAR, _PMVAR );
      _pDAG->eval( _opJAC, _PMWK, _nx*_nx, _pJAC, _MVXPdfdx, _nVAR0, _pVAR, _MVXPVAR );
      _RHS_PM_ELL1( _nx, _PMxdot, _MVXPdfdx, _Ir, _A, _Irdot );
    }

    // In this variant a polynomial model of the Jacobian matrix is computed and the
    // linear part is taken as the mid-point of this matrix
    else if( _MVXPenv->nord() <= _PMenv->nord() ){
      *_MVXPt = t; // set current time
      for( unsigned jx=0; jx<_nx; jx++ ){
        _MVXPd[jx].set( _MVXPenv, _npar+jx, _Ir[jx] );
        _MVXPx[jx].set( _MVXPenv ).set( _PMx[jx].center().set( T(0.) ), true );
      }
      _e2x( _nx, _MVXPd, _MVXPx, false );
#ifdef MC__ODEBND_BASE_DINEQPM_DEBUG
      _print_interm( t, _nx, _MVXPx, "MVXPx Intermediate", std::cerr );
#endif
      _pDAG->eval( _opRHS, _PMWK, _nx, _pRHS, _PMxdot, _nVAR0, _pVAR, _PMVAR );
      _pDAG->eval( _opJAC, _PMWK, _nx*_nx, _pJAC, _MVXPdfdx, _nVAR0, _pVAR, _MVXPVAR );
      _RHS_PM_ELL1( _nx, _PMxdot, _MVXPdfdx, _Ir, _A, _Irdot );
    }

    // In this variant a polynomial model in the joint state-parameter and
    // of the same order as the parameter polynomial model is computed
    else{
      *_MVXPt = t; // set current time   
      for( unsigned jx=0; jx<_nx; jx++ ){
        _MVXPd[jx].set( _MVXPenv, _npar+jx, _Ir[jx] );
        _MVXPx[jx].set( _MVXPenv ).set( _PMx[jx].center().set( T(0.) ), true );
      }
      _e2x( _nx, _MVXPd, _MVXPx, false );
#ifdef MC__ODEBND_BASE_DINEQPM_DEBUG
      _print_interm( t, _nx, _MVXPx, "MVXPx Intermediate", std::cerr );
#endif
      _pDAG->eval( _opRHS, _PMWK, _nx, _pRHS, _MVXPf, _nVAR0, _pVAR, _MVXPVAR );
      _RHS_PM_ELL2( _nx, _PMenv, _PMxdot, _MVXPf, _npar, _Ir, _A, _Irdot );
    }

    // Construct the ellipsoidal remainder derivatives
    if( !_RHS_ELL( _nx, _Q, _A, _Irdot, _Qdot, _QTOLx, machprec(), options.QSCALE, _Ir ) )
      return 1.;

    // Whether or not to ignore the remainder
    if( !options.PMNOREM )
      _PME2vec( _PMenv, _nx, _PMxdot, _Qdot, xdot );
    else
      _PME2vec( _PMenv, _nx, _PMxdot, 0, xdot );
#ifdef MC__ODEBND_BASE_DINEQPM_DEBUG
    _Er = E(_nx,_Qdot);
    _print_interm( t, _nx, _PMxdot, _Er, "PMxdot Intermediate", std::cerr );
    std::cout << "l =" << _Er.eigQ().first(0) << " ; " << _Er.eigQ().first(_nx-1) << std::endl;
    { std::cout << "--paused--"; int dum; std::cin >> dum; }
#endif
    return 0;
   }
  }
}

template <typename T, typename PMT, typename PVT>
inline void
ODEBND_BASE<T,PMT,PVT>::_RHS_PM_ELL0
( const unsigned nx, PVT*PMf, T*Idfdx, const T*Ir, double*Ar, T*Irdot )
{
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
( const unsigned nx, PVT*PMf, PVT*MVXPdfdx, const T*Ir, double*Ar, T*Irdot )
{
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
( const unsigned nx, PMT*PMenv, PVT*PMf, PVT*MVXPf, const unsigned np,
  const T*Ir, double*Ar, T*Irdot )
{
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
    if( !reinit ) break;
    _vec2PMI( x, _PMenv, _nx, _PMx, true );   // set current state polynomial model
#ifdef MC__ODEBND_BASE_DINEQPM_DEBUG
    _print_interm( t, _nx, _PMx, "PMx Intermediate", std::cerr );
#endif
    break;
   
  case OPT::DINEQ:
    if( !reinit ) break;
    _vec2PMI( x, _PMenv, _nx, _PMx, false );   // set current state polynomial model
#ifdef MC__ODEBND_BASE_DINEQPM_DEBUG
    _print_interm( t, _nx, _PMx, "PMx Intermediate", std::cerr );
#endif
    break;
   
  case OPT::ELLIPS:
  default:
    if( !reinit ) break;
    _vec2PME( x, _PMenv, _nx, _PMx, _Q, _Er, _Ir, true ); // set current state polynomial model
#ifdef MC__ODEBND_BASE_DINEQPM_DEBUG
    _print_interm( t, _nx, _PMx, _Er, "PMx Intermediate", std::cerr );
#endif
    break;
  }

  *_PMt = t; // set current time
  _QUAD_PM( _pDAG, _opQUAD, _PMWK, _nq, _pQUAD, _nVAR0, _pVAR, _PMVAR,
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
( FFGraph*DAG, std::list<const FFOp*>&opQUAD, std::vector<PVT>&PMQUAD,
  const unsigned nq, const FFVar*pQUAD, const unsigned nVAR,
  const FFVar*pVAR, PVT*PMVAR, PVT*PMqdot, const bool append )
{
  if( !nq ) return;
  DAG->eval( opQUAD, PMQUAD, nq, pQUAD, PMqdot, nVAR, pVAR, PMVAR, append );
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
( const unsigned pos_fct, const double t )
{
  if( !_nf ) return true;

  *_PMt = t; // set current time
  const FFVar* pFCT = _vFCT.at( pos_fct );
#ifdef MC__ODEBNDS_BASE_DINEQPM_DEBUG
    for( unsigned j=0; j<_nq; j++ )
      std::cout << "PMq[" << j << "] = " << _PMq[j] << std::endl;
#endif
  _pDAG->eval( _nf, pFCT, _PMf, _nVAR, _pVAR, _PMVAR, pos_fct?true:false );
  return true;
}

template <typename T, typename PMT, typename PVT>
template <typename OPT> inline bool
ODEBND_BASE<T,PMT,PVT>::_RHS_PM_SET
( const OPT&options, const unsigned iRHS, const unsigned iQUAD )
{
  if( !_RHS_SET( iRHS, iQUAD ) ) return false;

  delete[] _opRHSi; _opRHSi = 0;
  delete[] _pJAC;   _pJAC = 0; _opJAC.clear();

  switch( options.WRAPMIT){
  case OPT::DINEQ:
    _opRHSi = new std::list<const FFOp*>[_nx];
    for( unsigned ix=0; ix<_nx; ix++ )
      if( _pRHS ) _opRHSi[ix] = _pDAG->subgraph( 1, _pRHS+ix );
    break;

  case OPT::ELLIPS:
    if( _pRHS ) _pJAC = _pDAG->FAD( _nx, _pRHS, _nx, _pVAR );
    if( _pJAC ) _opJAC = _pDAG->subgraph( _nx*_nx, _pJAC );
    break;

  default:
    break;
  }

  return true;
}

template <typename T, typename PMT, typename PVT>
template <typename ODESLV> inline bool
ODEBND_BASE<T,PMT,PVT>::_bounds
( const T*Ip, T**Ixk, T*If, ODESLV&traj, const unsigned nsamp, 
  std::vector<Results>&results, std::ostream&os )
{
  int DISPLAY_SAVE = traj.options.DISPLAY;
  traj.options.DISPLAY = 0;

  // Initialization of sampled bounds at parameter lower bound
  double *p = new double[_np];
  for( unsigned ip=0; ip<_np; ip++ )
    p[ip] = Op<T>::l(Ip[ip]);
  double **xk = Ixk? new double*[_nsmax+1]: 0;
  for( unsigned is=0; Ixk && is<=_nsmax; is++ ){
    if( !Ixk[is] ) Ixk[is] = new T[_nx+_nq];
    xk[is] = new double[_nx+_nq];
  }
  double *f = If? new double[_nf]: 0;
  STATUS stat = traj.states( p, xk, f, os );
  if( stat != NORMAL || nsamp <= 1 ){
    delete[] p;
    delete[] f;
    for( unsigned is=0; is<=_nsmax; is++ ) delete[] xk[is];
    delete[] xk;
    return false;
  }   
  for( unsigned is=0; Ixk && is<=_nsmax; is++ )
    for( unsigned ix=0; ix<_nx+_nq; ix++ )
      Ixk[is][ix] = xk[is][ix];
  for( unsigned ifn=0; If && ifn<_nf; ifn++ )
    If[ifn] = f[ifn];

  // Initialize result vector with current trajectory
  results.clear();
  auto it=traj.results_sta.begin();
  for( ; traj.options.RESRECORD && it!=traj.results_sta.end(); ++it )
    results.push_back( Results( it->t, it->nx, it->X ) );

  // Start sampling process
  unsigned* vsamp = new unsigned[_np];
  bool flag = _sampling( Ip, Ixk, If, traj, nsamp, vsamp, 0, p, xk, f, results, os );
  traj.options.DISPLAY = DISPLAY_SAVE;

  // Clean-up
  delete[] p; delete[] f;
  for( unsigned is=0; xk && is<=_nsmax; is++ ) delete[] xk[is]; delete[] xk;
  delete[] vsamp;
  
  return flag;
}

template <typename T, typename PMT, typename PVT>
template <typename ODESLV> inline bool
ODEBND_BASE<T,PMT,PVT>::_sampling
( const T*Ip, T**Ixk, T*If, ODESLV&traj, const unsigned nsamp,
  unsigned* vsamp, const unsigned ipar, double*p, double**xk,
  double*f, std::vector<Results>&results, std::ostream&os )
{
  // Update bounds for all sampling points
  for( unsigned isamp=0; isamp<nsamp; isamp++ ){
    vsamp[ipar] = isamp;

    // Continue recursive call
    if( ipar+1 < _np ){
      if( !_sampling( Ip, Ixk, If, traj, nsamp, vsamp, ipar+1, p, xk, f, results, os ) )
        return false;
      continue;
    }

    // Update bounds for current point
#ifdef MC__ODEBND_BASE_SAMPLE_DEBUG
    std::cout << "Sample: ";
#endif
    for( unsigned ip=0; ip<_np; ip++ ){
      p[ip] = Op<T>::l( Ip[ip] ) + vsamp[ip]/(nsamp-1.) * Op<T>::diam( Ip[ip] );
#ifdef MC__ODEBND_BASE_SAMPLE_DEBUG
      std::cout << p[ip] << "  ";
#endif
    }
#ifdef MC__ODEBND_BASE_SAMPLE_DEBUG
    std::cout << std::endl;
#endif
    typename ODESLV::STATUS flag = traj.states( p, xk, f, os );
    if( flag != ODESLV::NORMAL ) return flag;
    for( unsigned is=0; Ixk && is<=_nsmax; is++ )
      for( unsigned ix=0; ix<_nx+_nq; ix++ )
        Ixk[is][ix] = Op<T>::hull( xk[is][ix], Ixk[is][ix] );
    for( unsigned ifn=0; If && ifn<_nf; ifn++ )
      If[ifn] = Op<T>::hull( f[ifn], If[ifn] );

    // Update result vector with current trajectory
    auto it=traj.results_sta.begin();
    for( unsigned k=0; traj.options.RESRECORD && it!=traj.results_sta.end(); ++it, k++ )
      for( unsigned i=0; i<it->nx; i++ )
        results[k].X[i] = Op<T>::hull( it->X[i], results[k].X[i] );
  }

  return true;
}  

template <typename T, typename PMT, typename PVT>
template <typename ODEBND> inline bool
ODEBND_BASE<T,PMT,PVT>::_hausdorff
( const T*Ip, double**Hxk, double*Hf, ODEBND&dineq, const unsigned nsamp,
  std::ostream&os )
{
  int DISPLAY_ODEBND = dineq.options.DISPLAY;
  dineq.options.DISPLAY = 0;

  // Compute inner bounds 
  T** Ixk0 = new T*[_nsmax+1];
  for( unsigned is=0; is<_nsmax+1; is++ ) Ixk0[is] = new T[_nx+_nq];
  T* If0 = _nf? new T[_nf]: 0;
  if( dineq.bounds( nsamp, Ip, Ixk0, If0, os ) != ODEBND::NORMAL ){
    for( unsigned is=0; is<_nsmax+1; is++ ) delete[] Ixk0[is];
    delete[] Ixk0;
    delete[] If0;
    return false;
  }

  // Compute outer bounds
  T** Ixk = new T*[_nsmax+1];
  for( unsigned is=0; is<_nsmax+1; is++ ) Ixk[is] = new T[_nx+_nq];
  T* If = _nf? new T[_nf]: 0;
  try{ dineq.bounds( Ip, Ixk, If, os ); }
  catch(...){;}
  unsigned nsf = dineq.final_stage();

  dineq.options.DISPLAY = DISPLAY_ODEBND;
  double **Hxk_ = Hxk? Hxk: new double*[_nsmax+1];
  for( unsigned is=0; is<=_nsmax; is++ ){
    if( !Hxk_[is] ) Hxk_[is] = new double[_nx+_nq];
    for( unsigned i=0; i<_nx+_nq; i++ )
      Hxk_[is][i] = is<=nsf? _dH( Ixk[is][i], Ixk0[is][i] ): 0./0.;
    if( dineq.options.DISPLAY >= 1 )
      _print_interm( _dT[is], _nx, Hxk_[is], "dHx", os );
      _print_interm( _nq, Hxk_[is]+_nx, "dHq", os );
  }
  double *Hf_ = Hf? Hf: new double[_nf];
  for( unsigned i=0; i<_nf; i++ )
    Hf_[i] = nsf==_nsmax? _dH( If[i], If0[i] ): 0./0.;
  if( dineq.options.DISPLAY >= 1 && _nf )
    _print_interm( _nf, Hf_, "dHf", os );

  for( unsigned is=0; is<_nsmax+1; is++ ) delete[] Ixk0[is];
  delete[] Ixk0;
  delete[] If0;
  for( unsigned is=0; is<_nsmax+1; is++ ) delete[] Ixk[is];
  delete[] Ixk;
  delete[] If;
  for( unsigned is=0; !Hxk && is<_nsmax+1; is++ ) delete[] Hxk_[is];
  if( !Hxk ) delete[] Hxk_;
  if( !Hf )  delete[] Hf_;

  return true;
}

template <typename T, typename PMT, typename PVT>
template <typename ODEBND, typename ODESLV> inline bool
ODEBND_BASE<T,PMT,PVT>::_hausdorff
( const PVT*PMp, double**Hxk, double*Hf, ODEBND&dineq, ODESLV&traj,
  const unsigned nsamp, std::ostream&os )
{
  int DISPLAY_ODEBND = dineq.options.DISPLAY;
  int DISPLAY_ODESLV = traj.options.DISPLAY;
  dineq.options.DISPLAY = traj.options.DISPLAY = 0;

  // Compute inner bounds
  PVT** PMxk = new PVT*[_nsmax+1];
  for( unsigned is=0; is<_nsmax+1; is++ ){
    PMxk[is] = new PVT[_nx+_nq];
    for( unsigned ix=0; ix<_nx+_nq; ix++ )
      PMxk[is][ix] = 0./0.;
  }
  PVT* PMf = _nf? new PVT[_nf]: 0;
  try{ dineq.bounds( PMp, PMxk, PMf, os ); }
  catch(...){;}
  unsigned nsf = dineq.final_stage();
  traj.set_time( nsf, _dT.data(), _pT );

  // Compute remainder (outer) bounds 
  T* Ip = new T[_npar];
  for( unsigned ip=0; ip<_npar; ip++ ) Ip[ip] = PMp[ip].B();
  T** Rxk = new T*[_nsmax+1];
  for( unsigned is=0; is<_nsmax+1; is++ ) Rxk[is] = new T[_nx+_nq];
  T* Rf = _nf? new T[_nf]: 0;
  bool flag = _remainders( nsf, _npar, Ip, PMxk, PMf, Rxk, Rf, traj, nsamp, os );

  dineq.options.DISPLAY = DISPLAY_ODEBND;
  traj.options.DISPLAY = DISPLAY_ODESLV;
  double **Hxk_ = Hxk? Hxk: new double*[_nsmax+1];
  for( unsigned is=0; is<=_nsmax; is++ ){
    if( !Hxk_[is] ) Hxk_[is] = new double[_nx+_nq];
    for( unsigned ix=0; ix<_nx+_nq; ix++ )
      Hxk_[is][ix] = is<=nsf? _dH( PMxk[is][ix].R(), Rxk[is][ix] ): 0./0.;
    if( dineq.options.DISPLAY >= 1 ){
      _print_interm( _dT[is], _nx, Hxk_[is], "dHx", os );
      _print_interm( _nq, Hxk_[is]+_nx, "dHq", os );
    }
  }
  double *Hf_ = Hf? Hf: new double[_nf];
  for( unsigned i=0; i<_nf; i++ )
    Hf_[i] = nsf==_nsmax? _dH( PMf[i].R(), Rf[i] ): 0./0.;
  if( dineq.options.DISPLAY >= 1 && _nf )
    _print_interm( _nf, Hf_, "dHf", os );

  for( unsigned is=0; is<_nsmax+1; is++ ) delete[] Rxk[is];
  delete[] Rxk;
  delete[] Rf;
  for( unsigned is=0; is<_nsmax+1; is++ ) delete[] PMxk[is];
  delete[] PMxk;
  delete[] PMf;
  delete[] Ip;
  for( unsigned is=0; !Hxk && is<_nsmax+1; is++ ) delete[] Hxk_[is];
  if( !Hxk ) delete[] Hxk_;
  if( !Hf )  delete[] Hf_;

  return flag;
}

template <typename T, typename PMT, typename PVT>
template<typename ODESLV> inline bool
ODEBND_BASE<T,PMT,PVT>::_remainders
( const unsigned ns, const unsigned np, const T*Ip, const PVT*const*PMxk,
  const PVT*PMf, T**Rxk, T*Rf, ODESLV&traj, const unsigned nsamp,
  std::ostream&os )
{
   // Initialization of sampled bounds at parameter lower bound
  double *p = new double[np];
  for( unsigned ip=0; ip<np; ip++ )
    p[ip] = Op<T>::l(Ip[ip]);
  double **xk = new double*[ns+1];
  for( unsigned is=0; is<=ns; is++ )
    xk[is] = new double[_nx+_nq];
  double *f = _nf? new double[_nf]:0;
  typename ODESLV::STATUS flag = traj.states( p, xk, f, os );
  if( flag != ODESLV::NORMAL || nsamp <= 1 ){
    delete[] p; delete[] f;
    for( unsigned is=0; is<=ns; is++ ) delete[] xk[is]; delete[] xk;
    return( flag==ODESLV::NORMAL? true: false );
  }   
  for( unsigned is=0; is<=ns; is++ )
    for( unsigned i=0; i<_nx+_nq; i++ )
      Rxk[is][i] = xk[is][i] - PMxk[is][i].polynomial( p );
  for( unsigned i=0; i<_nf; i++ )
    Rf[i] = f[i] - PMf[i].polynomial(p);

  // Start sampling process
  unsigned* vsamp = new unsigned[np];
  bool flag2 = _remainders( ns, np, Ip, PMxk, PMf, Rxk, Rf, traj, nsamp,
                            vsamp, 0, p, xk, f, os );
  
  // Clean-up
  delete[] p; delete[] f;
  for( unsigned is=0; is<=ns; is++ ) delete[] xk[is]; delete[] xk;
  delete[] vsamp;
  
  return flag2;
}

template <typename T, typename PMT, typename PVT>
template<typename ODESLV> inline bool
ODEBND_BASE<T,PMT,PVT>::_remainders
( const unsigned ns, const unsigned np, const T*Ip, const PVT*const*PMxk,
  const PVT*PMf, T**Rxk, T*Rf, ODESLV&traj, const unsigned nsamp, unsigned* vsamp,
  const unsigned ip, double*p, double**xk, double*f, std::ostream&os )
{
  typename ODESLV::STATUS flag = ODESLV::NORMAL;

  // Update bounds for all sampling points
  for( unsigned isamp=0; isamp<nsamp; isamp++ ){
    vsamp[ip] = isamp;

    // Continue recursive call
    if( ip+1 < np ){
      if( !_remainders( ns, np, Ip, PMxk, PMf, Rxk, Rf, traj, nsamp, vsamp,
                        ip+1, p, xk, f, os ) ) return false;
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
    flag = traj.states( p, xk, f, os );
    if( flag != ODESLV::NORMAL ) return false;
    for( unsigned is=0; is<=ns; is++ )
      for( unsigned i=0; i<_nx+_nq; i++ )
        Rxk[is][i] = Op<T>::hull( xk[is][i]-PMxk[is][i].polynomial(p), Rxk[is][i] );
    for( unsigned i=0; i<_nf; i++ )
      Rf[i] = Op<T>::hull( f[i]-PMf[i].polynomial(p), Rf[i] );
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
template <typename VRES> inline void
ODEBND_BASE<T,PMT,PVT>::_record
( std::ofstream&ofile, const VRES&bnd, const unsigned iprec )
{
  if( !ofile ) return;

  // Specify format
  ofile << std::right << std::scientific << std::setprecision(iprec);

  // Record computed interval bounds at stage times
  auto it = bnd.begin(), it0 = it;
  for( ; it != bnd.end(); ++it ){
    if( it != bnd.begin() && it->t == it0->t )
      ofile << std::endl;
    ofile << std::setw(iprec+9) << it->t;
    for( unsigned ix=0; ix<it->nx; ix++ )
      ofile << std::setw(iprec+9) << mc::Op<T>::l( it->X[ix] )
            << std::setw(iprec+9) << mc::Op<T>::u( it->X[ix] );
    ofile << std::endl;
    it0 = it;
  }
}

} // end namescape mc

#endif

