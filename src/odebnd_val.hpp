// Copyright (C) 2012-2014 Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

#ifndef MC__ODEBND_VAL_HPP
#define MC__ODEBND_VAL_HPP

#undef  MC__ODEBND_VAL_DEBUG
#undef  MC__ODEBND_VAL_DEBUG_INVARIANT
#undef  MC__ODEBND_VAL_STEPSIZE_HSTAB

#include <stdexcept>
#include <cassert>
#include <cmath>
#include <fstream>
#include <vector>
#include <sys/time.h>

#include "ellipsoid.hpp"
#include "tmodel.hpp"
#include "cmodel.hpp"

#include "odeslv_gsl.hpp"

namespace mc
{
//! @brief C++ class computing (validated) enclosures of the reachable set of parametric ODEs using Taylor series expansion and polynomial models with convex remainder terms
////////////////////////////////////////////////////////////////////////
//! mc::ODEBND_VAL is a C++ class that computes enclosures of the
//! reachable set of parametric ordinary differential equations
//! (ODEs) using MC++. It implements a validated method based on Taylor
//! series expansion in time of the ODE solutions, whereby polynomial
//! models with interval or ellipsoidal remainders are used to enable
//! high-order convergence. The use of ellipsoidal remainders enables
//! stability of the enclosures for asymptotically stable ODE systems
//! when the parameter host is sufficiently small.
////////////////////////////////////////////////////////////////////////
template <typename T, typename PMT=mc::TModel<T>, typename PVT=mc::TVar<T> >
class ODEBND_VAL: public BASE_DE                
{
private:

  typedef Ellipsoid E;

  //! @brief number of independent variables in DAG
  unsigned int _nVAR;

  //! @brief pointer to independent variables in DAG
  FFVar* _pVAR;

  //! @brief pointer to T bounds of independent variables (state, parameter, time)
  T *_IVAR;

  //! @brief Polynomial model environment
  PMT *_PMenv;

  //! @brief pointer to PVT bounds of independent variables (state, parameter, time)
  PVT *_PMVAR;

  //! @brief internal polynomial model environment for ODE right-hand side
  PMT *_MVenv;

  //! @brief pointer to PVT bounds of independent variables (state, parameter, time)
  PVT *_MVVAR;


  //! @brief state interval bounds
  T *_Ix;

  //! @brief state polynomial model
  PVT *_PMx;

  //! @brief state ellipsoidal remainder bounds in polynomial model
  E _Ex;

  //! @brief state ellipsoidal remainder bounds in polynomial model at expansion point
  E _Ex_TE;

  //! @brief state interval remainder bounds in polynomial model
  T *_Rx;

  //! @brief state interval remainder bounds in polynomial model at expansion point
  T *_Rx_TE;

  //! @brief RHS derivatives w.r.t. state
  CPPL::dgematrix _Ax;

  //! @brief Invariant remainder part
  T *_Rinv;

  //! @brief Invariant remainder enclosure
  T *_Ninv;

  //! @brief Invariant linear part
  CPPL::dgematrix _Ainv;


  //! @brief preallocated array for DAG evaluation in T arithmetic
  T* _IDAG;

  //! @brief preallocated array for DAG evaluation in PVT arithmetic
  PVT* _PMDAG;


  //! @brief list of operations in RHS Taylor coefficients
  std::list<const FFOp*> _opTRHS;

  //! @brief const pointer to RHS Taylor coefficients in DAG
  const FFVar* _pTRHS;

  //! @brief pointer to PVT bounds of RHS Taylor coefficients
  PVT *_PMTRHS;


  //! @brief list of operations in RHS Taylor coefficient derivatives
  std::list<const FFOp*> _opFTRHS;

  //! @brief const pointer to RHS Taylor coefficient derivatives in DAG
  const FFVar* _pFTRHS;

  //! @brief pointer to T bounds of RHS Taylor coefficient derivatives
  T *_IFTRHS;

  //! @brief pointer to PVT bounds of RHS Taylor coefficient derivatives
  PVT *_PMFTRHS;


  //! @brief list of operations in RHS highest-order Taylor coefficient
  std::list<const FFOp*> _opTRHSval;

  //! @brief pointer to T bounds of RHS highest-order Taylor coefficient
  T *_ITRHSval;

  //! @brief pointer to PVT bounds of RHS highest-order Taylor coefficient
  PVT *_PMTRHSval;


  //! @brief list of operations in initial conditions
  std::list<const FFOp*> _opIC;

  //! @brief const pointer to IC function in DAG
  const FFVar* _pIC;


  //! @brief list of operations in invariants
  std::list<const FFOp*> _opINV;

  //! @brief const pointer to invariants in DAG
  const FFVar* _pINV;

  //! @brief pointer to PVT bounds of invariants
  PVT *_PMINV;


  //! @brief list of operations in invariant derivatives
  std::list<const FFOp*> _opFINV;

  //! @brief const pointer to invariant derivatives in DAG
  const FFVar* _pFINV;

  //! @brief pointer to T bounds of invariant derivatives
  T *_IFINV;

  //! @brief pointer to PVT bounds of invariant derivatives
  PVT *_PMFINV;

  //! @brief current time
  double _t;

  //! @brief current stage
  unsigned _istg;

  //! @brief stepsize
  double _h;

  //! @brief maximal valid stepsize
  double _hmax;

  //! @brief Maximal propagation time
  double _tmax;

  //! @brief Time for predictor validation
  double _tval;

  //! @brief ODE solver based on GSL
  ODESLV_GSL<T>* _ODESLV_GSL;

public:
  /** @defgroup ODEBND_VAL Discretized (validated) set-valued integration of parametric ODEs
   *  @{
   */
  //! @brief Default constructor
  ODEBND_VAL();

  //! @brief Default destructor
  virtual ~ODEBND_VAL();

  //! @brief Integration results at a given time instant
  struct Results
  {
    //! @brief Constructors
    Results
      ( const double tk, const unsigned int nxk, const T*Ixk ):
      t( tk ), nx( nxk )
      { X = new T[nx];
        for( unsigned int ix=0; ix<nx; ix++ ) X[ix] = Ixk[ix]; }
    Results
      ( const double tk, const unsigned int nxk, const PVT*PMxk ):
      t( tk ), nx( nxk )
      { X = new T[nx];
        for( unsigned int ix=0; ix<nx; ix++ ) X[ix] = PMxk[ix].B(); }
    Results
      ( const Results&res ):
      t( res.t ), nx( res.nx )
      { X = new T[nx];
        for( unsigned int ix=0; ix<nx; ix++ ) X[ix] = res.X[ix]; }
    //! @brief Destructor
    ~Results()
      { delete[] X; }
    //! @brief Time instant
    double t;
    //! @brief Solution dimension
    unsigned int nx;
    //! @brief Solution bounds
    T* X;
  };

  //! @brief Integrator options
  struct Options
  {
    //! @brief Constructor
    Options():
      TSORDER(7), WRAPMIT(ELLIPS), ORDMIT(1), DMAX(1e20), HMIN(1e-8), HMAX(1e8),
      HREDUC(0.8), HSTAB(true), TOL(1e-7), SCALING(true), ATOL(1e-10),
      QTOL(machprec()), TSTOP(true), PMVALID(false), USEINV(true),
      ODESLVOPT(typename ODESLV_GSL<T>::Options()), DISPLAY(1), RESRECORD(false)
      {}
    //! @brief Enumeration of wrapping mitigation strategies
    enum WRAPPING_STRATEGY{
      NONE=0,		//!< No wrapping mitigation
      ELLIPS		//!< Ellipsoidal contractor with linear preconditioning [Default]
    };
    //! @brief Taylor series expansion order (Default: 7)
    unsigned int TSORDER;
    //! @brief Wrapping mitigation strategy
    WRAPPING_STRATEGY WRAPMIT;
    //! @brief Order of wrapping mitigation strategy (Default: 1)
    unsigned int ORDMIT;
    //! @brief Maximum enclosure diameter, \f$D_{\rm max}\f$ (Default: 1e20)
    double DMAX;
    //! @brief Minimum step-size, \f$h_{\rm min}\f$ (Default: 1e-8)
    double HMIN;
    //! @brief Maximum step-size, \f$h_{\rm max}\f$ (Default: 1e8)
    double HMAX;
    //! @brief Reduction factor for step-size validation between 0-1 (default: 0.8)
    double HREDUC;
    //! @brief Whether or not to consider stepsize dependence in Taylor truncation term (default: true)
    double HSTAB;
    //! @brief Tolerance on truncation error for step size selection (Default: 1e-7)
    double TOL;
    //! @brief Whether or not to adjust scaling for state components (default: true)
    bool SCALING;
    //! @brief Absolute tolerance for scaling (Default: 1e-10)
    double ATOL;
    //! @brief Tolerance when dividing by trace of shape matrix in ellipsoidal bounds (Default: mc::machprec())
    double QTOL;
    //! @brief Whether to stop and reinitialize the integrator at time steps (default: true)
    bool TSTOP;
    //! @brief Whether or not to use polynomial models for stepsize validation (Default: false)
    bool PMVALID;
    //! @brief Whether or not to use the specified invariants for bounds contraction (default: true)
    bool USEINV;
     //! @brief Options of non-validated ODE solver for reachable set approximation
    typename ODESLV_GSL<T>::Options ODESLVOPT;
  //! @brief Display level (default: 1)
    int DISPLAY;
    //! @brief Whether or not to record results (default: false)
    bool RESRECORD;
  } options;

  //! @brief Structure storing integration statistics
  struct Stats
  {
    //! @brief Constructor
    Stats():
      cputime(0.), numSteps(0), numTRHS_I(0), numFTRHS_I(0), numTRHS_PM(0), numFTRHS_PM(0)
      {}
    //! @brief Constructor
    void reset()
      { cputime = 0.;
        numSteps = numTRHS_I = numFTRHS_I = numTRHS_PM = numFTRHS_PM = 0; }

    //! @brief CPU time
    double cputime;
    //! @brief Step count
    unsigned long numSteps;
    //! @brief Number of RHS Taylor expansion evaluations in T arithmetic
    unsigned long numTRHS_I;
    //! @brief Number of RHS Taylor expansion derivative evaluations in T arithmetic
    unsigned long numFTRHS_I;
    //! @brief Number of RHS Taylor expansion evaluations in PM arithmetic
    unsigned long numTRHS_PM;
    //! @brief Number of RHS Taylor expansion derivative evaluations in PM arithmetic
    unsigned long numFTRHS_PM;
  };

  //! @brief Statistics for state bounds integration
  Stats stats_traj;

  //! @brief Structure for setting up storing the solver exceptions
  class Exceptions
  {
  public:
    //! @brief Enumeration type for ODEBND_VAL exception handling
    enum TYPE{
      REINIT=1,	//!< Error due to state reinitialization at intermediate stage
      HREDUC,	//!< Error due to an invalid stepsize reduction parameter
      UNDEF=-33	//!< Error due to calling a function/feature not yet implemented
    };
    //! @brief Constructor for error <a>ierr</a>
    Exceptions( TYPE ierr ) : _ierr( ierr ){}
    //! @brief Inline function returning the error flag
    int ierr(){ return _ierr; }
    //! @brief Inline function returning the error description
    std::string what(){
      switch( _ierr ){
      case REINIT:
        return "ODEBND_VAL::Exceptions  State reinitialization at intermediate stage not allowed";
      case HREDUC:
        return "ODEBND_VAL::Exceptions  Invalid stepsize reduction parameter";
      case UNDEF: default:
        return "ODEBND_VAL::Exceptions  Calling a feature not yet implemented";
      }
    }
  private:
    TYPE _ierr;
  };

  //! @brief Computes polynomial model enclosure of reachable set of parametric ODEs
  STATUS bounds
    ( const unsigned int ns, const double*tk, const PVT*PMp, PVT**PMxk,
      E*ERxk=0, std::ostream&os=std::cout );

  //! @brief Computes Hausdorff distance between polynomial model remainder enclosure and actual remainder function range, using parameter sampling
  STATUS hausdorff
    ( const unsigned int ns, const double*tk, const PVT*PMp, double**Hxk,
      const unsigned int nsamp, std::ostream&os=std::cout );

  //! @brief Record results in file <a>bndrec</a>, with accuracy of <a>iprec</a> digits
  void record
    ( std::ofstream&bndrec, const unsigned int iprec=5 ) const;

  //! @brief Return value of final time reached
  double final_time() const
    { return _t; };
  /** @} */

private:
  //! @brief Function to initialize integrator statistics
  static void _init_stats
    ( Stats&stats );

  //! @brief Function to finalize integrator statistics
  static void _final_stats
    ( Stats&stats );

  //! @brief Function to display integrator statistics
  static void _print_stats
    ( const Stats&stats, std::ostream&os=std::cout );

  //! @brief Vector storing interval bound results (upon request only)
  std::vector< Results > _results;

  //! @brief Function doing scaling
  template <typename U> double _scale
    ( U&R, const T&X ) const;

  //! @brief Function computing scaled max-norm
  template <typename U> double _smaxnorm
    ( U*R, const T*X ) const;

  //! @brief Function computing scaled 1-norm
  template <typename U> double _snorm
    ( U*R, const T*X ) const;

  //! @brief Function computing set diameter
  template <typename U> double _diam
    ( U*X ) const;

  //! @brief Compute stepsize using Taylor remainder
  template <typename U> double _stepsize
    ( const unsigned q, U*UTRHS, const T*X, const bool hdep ) const;

  //!@brief Select and validate stepsize for current predictor
  double _validation
    ( const unsigned q, PVT*PMTRHS, const T*X );

  //! @brief Function to compute the Taylor-based predictor for state component <a>ix</a> within a (possibly interval valued) step-size <a>h</a> in T arithmetic
  template <typename U, typename V> static U _predictor
    ( const V&h, const unsigned q, const unsigned n, const U*U_TE,
      const unsigned i );

  //! @brief Function to compute the Taylor-based predictor for derivative component <a>jx</a> of state component <a>ix</a> within a (possibly interval valued) step-size <a>h</a> in T arithmetic
  template <typename U, typename V> static U _predictor
    ( const V&h, const unsigned q, const unsigned n, const U*U_FTE,
      const unsigned i, const unsigned j );

  //! @brief Function to prepare for state bound propagation
  bool _prepare
    ( const PVT*PMp );

  //! @brief Function to initialize state polynomial models
  bool _init
    ( const unsigned k, const double tk, PVT*PMxk, E*Exk );

  //! @brief Function to set/size RHS evaluation
  bool _resize
    ( const unsigned q, const unsigned int iRHS );

  //! @brief Function to propagate enclosures in PMT arithmetic
  STATUS _propagate
    ( const double tnext, std::ostream&os );

  //! @brief Function to propagate enclosures in PMT arithmetic with interval remainder
  STATUS _propagate_INT
    ( const double tnext, std::ostream&os );

  //! @brief Function to propagate enclosures in PMT arithmetic with ellipsoidal remainder
  STATUS _propagate_ELL
    ( const double tnext, std::ostream&os );

  //! @brief Function to propagate enclosures in PMT arithmetic with ellipsoidal remainder -- approximation using mean-value theorem and interval analysis
  STATUS _propagate_ELL0
    ( const unsigned q, const double tnext, std::ostream&os );

  //! @brief Function to propagate enclosures in PMT arithmetic with ellipsoidal remainder -- approximation using mean-value theorem and polynomial model arithmetic
  STATUS _propagate_ELL1
    ( const unsigned q, const double tnext, std::ostream&os );

  //! @brief Function to propagate enclosures in PMT arithmetic with ellipsoidal remainder -- joint polynomial model in states and parameters
  STATUS _propagate_ELL2
    ( const unsigned q, const double tnext, std::ostream&os );

  //! @brief Function to contract ellipsoidal remainder enclosure in PMT arithmetic using ODE invariants
  STATUS _invariant_ELL
    ( std::ostream&os );

  //! @brief Function to contract ellipsoidal remainder enclosure in PMT arithmetic using ODE invariants -- approximation using mean-value theorem and interval analysis
  STATUS _invariant_ELL0
    ( const PVT*Px, const T*Rx, std::ostream&os );

  //! @brief Function to contract ellipsoidal remainder enclosure in PMT arithmetic using ODE invariants -- approximation using mean-value theorem and polynomial model arithmetic
  STATUS _invariant_ELL1
    ( const PVT*Px, const T*Rx, std::ostream&os );

  //! @brief Function to contract ellipsoidal remainder enclosure in PMT arithmetic using ODE invariants -- joint polynomial model in states and parameters
  STATUS _invariant_ELL2
    ( const PVT*Px, const T*Rx, std::ostream&os );

  //! @brief Function to update the enclosures in PMT arithmetic with ellipsoidal remainder
  STATUS _update_ELL
    ( std::ostream&os );

  //! @brief Function to bound the remainder function relative to a polynomial model at sampling points -- return value is status
  STATUS _remainders
    ( const unsigned int ns, const double*tk, const T*Ip, const PVT*const*PMxk,
      T**Rxk, const unsigned int nsamp, std::ostream&os=std::cout );

  //! @brief Recrusive function computing bounds on errors between solutions of IVP in ODEs and polynomial approximant using sampling
  STATUS _remainders
    ( const unsigned int ns, const double*tk, const T*Ip, const PVT*const*PMxk,
      T**Rxk, const unsigned int nsamp, unsigned int* vsamp,
      const unsigned int ip, double*p, double**xk, std::ostream&os );

  //! @brief Function to display intermediate results
  template<typename U> void _print_interm
    ( const double t, const U*x, const std::string&var, std::ostream&os=std::cout ) const;

  //! @brief Function to display intermediate results
  template<typename U, typename V> void _print_interm
    ( const double t, const U*x, const V&r, const std::string&var, std::ostream&os=std::cout ) const;

  //! @brief Recursive function computing the factorial of an integer
  static double _factorial
    ( const unsigned int n )
    {
      return( n==1? 1: n*_factorial(n-1) );
    }

  //! @brief Private methods to block default compiler methods
  ODEBND_VAL(const ODEBND_VAL&);
  ODEBND_VAL& operator=(const ODEBND_VAL&);
};

template <typename T, typename PMT, typename PVT>
inline
ODEBND_VAL<T,PMT,PVT>::ODEBND_VAL
() : BASE_DE()
{
  // Initalize DAG computation
  _nVAR = 0;
  _pVAR = 0;
  _PMVAR = _MVVAR = _PMDAG = 0;
  _IVAR = _IDAG = 0;

  // Initalize dynamic system computation
  _pTRHS = _pFTRHS = _pIC = _pINV = _pFINV = 0;
  _PMenv = _MVenv = 0;
  _Ix = _Rx = _Rx_TE = _Rinv = _Ninv = 0;
  _PMx = 0;
  _IFTRHS = _ITRHSval = _IFINV = 0;
  _PMTRHS = _PMFTRHS = _PMTRHSval = _PMINV = _PMFINV = 0;

  // Initialize ellipsoidal calculus
#ifdef MC__ODEBND_VAL_DEBUG
  E::options.PSDCHK = true;
#else
  E::options.PSDCHK = false;
#endif

  // Initalize GSL
  _ODESLV_GSL = new ODESLV_GSL<T>();
}

template <typename T, typename PMT, typename PVT>
inline
ODEBND_VAL<T,PMT,PVT>::~ODEBND_VAL
()
{
  delete[] _pVAR;
  delete[] _IVAR;
  delete[] _PMVAR;
  delete[] _MVVAR;
  delete[] _IDAG;
  delete[] _PMDAG;

  // Do *NOT* free _PMenv
  delete _MVenv;
  delete[] _PMx;
  delete[] _Ix;
  delete[] _Rx;
  delete[] _Rx_TE;
  delete[] _Rinv;
  delete[] _Ninv;

  delete[] _pTRHS;
  delete[] _PMTRHS;
  delete[] _pFTRHS;
  delete[] _PMFTRHS;
  delete[] _IFTRHS;
  delete[] _PMTRHSval;
  delete[] _ITRHSval;
  // Do *NOT* free _pIC;
  // Do *NOT* free _pINV;
  delete[] _PMINV;
  delete[] _pFINV;
  delete[] _PMFINV;
  delete[] _IFINV;

  delete _ODESLV_GSL;
}

template <typename T, typename PMT, typename PVT>
template <typename U> inline double
ODEBND_VAL<T,PMT,PVT>::_scale
( U&R, const T&X ) const
{
  return( options.SCALING?
          Op<U>::abs(R)/(2.*options.TOL*Op<T>::abs(X)+options.ATOL):
	  Op<U>::abs(R)/options.TOL );
}

template <typename T, typename PMT, typename PVT>
template <typename U> inline double
ODEBND_VAL<T,PMT,PVT>::_smaxnorm
( U*R, const T*X ) const
{
  double maxnorm = 0.;
  for( unsigned int ix=0; ix<_nx; ix++ )
    maxnorm = std::max( maxnorm, _scale(R[ix],X[ix]) );
  return maxnorm;
}

template <typename T, typename PMT, typename PVT>
template <typename U> inline double
ODEBND_VAL<T,PMT,PVT>::_snorm
( U*R, const T*X ) const
{
  double snorm = 0.;
  for( unsigned int ix=0; ix<_nx; ix++ )
    snorm += _scale(R[ix],X[ix]);
  return snorm/(double)_nx;
}

template <typename T, typename PMT, typename PVT>
template <typename U> inline double
ODEBND_VAL<T,PMT,PVT>::_diam
( U*X ) const
{
  double diam = 0.;
  for( unsigned int ix=0; ix<_nx; ix++ )
    diam = std::max( diam, Op<U>::diam(X[ix]) );
  return diam;
}

template <typename T, typename PMT, typename PVT>
template <typename U, typename V> inline U
ODEBND_VAL<T,PMT,PVT>::_predictor
( const V&h, const unsigned q, const unsigned n, const U*U_TE,
  const unsigned i )
{
#ifdef MC__ODEBND_VAL_DEBUG
  for( unsigned k=0; k<q; q++ )
    std::cout << "ODEBND_VAL::_predictor *** f[" << q << "](" << i << ") = "
              << U_TE[q*n+i] << std::endl;
#endif
  // Compute predictor of xi(_t+h)
  unsigned pTEi = q*n+i;
  U Xi( U_TE[pTEi] );
  pTEi -= n;
  for( unsigned k=0; k<q; k++, pTEi-=n )
    { Xi *= h; Xi += U_TE[pTEi]; }
  return Xi;

}

template <typename T, typename PMT, typename PVT>
template <typename U, typename V> inline U
ODEBND_VAL<T,PMT,PVT>::_predictor
( const V&h, const unsigned q, const unsigned n, const U*U_FTE,
  const unsigned i, const unsigned j )
{
  // Compute predictor of dxi(_t+h)/dx
  unsigned pFTEij = q*n*n+i*n+j;
  U DjXi( U_FTE[pFTEij] );
  pFTEij -= n*n;
  for( unsigned k=0; k<q; k++, pFTEij-=n*n )
    { DjXi *= h; DjXi += U_FTE[pFTEij]; };
  return DjXi;
}

template <typename T, typename PMT, typename PVT>
template <typename U> inline double
ODEBND_VAL<T,PMT,PVT>::_stepsize
( const unsigned q, U*UTRHS, const T*X, const bool hdep ) const
{
#ifdef MC__ODEBND_VAL_DEBUG_STEPSIZE
  std::cout << "UTRHS[q+1]:";
  for( unsigned int ix=0; ix<_nx; ix++ )
    std::cout << "  " << UTRHS[ix];
  std::cout << std::endl << "factorial: " << _factorial(q+1) << std::endl;
#endif
  return hdep? std::pow( _smaxnorm(UTRHS,X)*_factorial(q+1) + machprec(), -1./q )
             : std::pow( _smaxnorm(UTRHS,X)*_factorial(q+1) + machprec(), -1./(q+1) );
}

template <typename T, typename PMT, typename PVT>
inline double
ODEBND_VAL<T,PMT,PVT>::_validation
( const unsigned q, PVT*PMTRHS, const T*X )
{
  if( options.HREDUC >= 1. || options.HREDUC <=0. ) throw Exceptions( Exceptions::HREDUC );

  // Stepsize validatinon
#ifdef MC__ODEBND_VAL_DEBUG_STEPSIZE
  for( unsigned k=1; k<=q; k++ ){
    double h0 = _stepsize( k, PMTRHS+k*_nx, X, options.HSTAB );
    std::cout << k << "  " << h0 << std::endl;
  }
  { int dum; std::cin >> dum; }
#endif
  double h = _stepsize( q, PMTRHS+q*_nx, X, options.HSTAB ) * options.HREDUC, hval;
  for( ; true; h *= options.HREDUC ){

    // Check minimum stepsize criterion
    if( h < options.HMIN ) return h;

    // Check predictor validity on full step [0,h]
    if( options.PMVALID ){ // polynomial model-based validation
      for( unsigned ix=0; ix<_nx; ix++ ){
        _PMVAR[ix] = _predictor( T(0,1)*h, q, _nx, PMTRHS, ix ).center();
        _PMVAR[ix] += options.HSTAB? h*options.TOL*T(-1,1): options.TOL*T(-1,1);
      }
      _PMVAR[_nx+_np] = _t; // Keep time when stepsize was computed
      _pDAG->eval( _opTRHSval, _PMDAG, _nx, _pTRHS+(q+1)*_nx, _PMTRHSval, _nVAR, _pVAR, _PMVAR );
      ++stats_traj.numTRHS_PM;
      hval = _stepsize(q, _PMTRHSval, X, options.HSTAB );
    }
    else{ // Interval-based validation
      for( unsigned ix=0; ix<_nx; ix++ ){
        _IVAR[ix] = _predictor( T(0,1)*h, q, _nx, PMTRHS, ix ).bound();
        _IVAR[ix] += options.HSTAB? h*options.TOL*T(-1,1): options.TOL*T(-1,1);
      }
      _IVAR[_nx+_np] = _t; // Keep time when stepsize was computed
      _pDAG->eval( _opTRHSval, _IDAG, _nx, _pTRHS+(q+1)*_nx, _ITRHSval, _nVAR, _pVAR, _IVAR );
      ++stats_traj.numTRHS_I;
      hval = _stepsize( q, _ITRHSval, X, options.HSTAB );
    }
#ifdef MC__ODEBND_VAL_DEBUG_STEPSIZE
  std::cout << "current step: h=" << h << " <=? " << hval << std::endl;
#endif
    if( h < hval ) return h;
  }
}

template <typename T, typename PMT, typename PVT>
inline typename ODEBND_VAL<T,PMT,PVT>::STATUS
ODEBND_VAL<T,PMT,PVT>::_propagate_INT
( const double tnext, std::ostream&os )
{
  const unsigned q = options.TSORDER;

  // Taylor series expansion in PVT arithmetic up to order TSORDER+1
  for( unsigned ix=0; ix<_nx; ix++ ) _PMVAR[ix] = _PMx[ix];
  _PMVAR[_nx+_np] = _t;
  _pDAG->eval( _opTRHS, _PMDAG, (q+2)*_nx, _pTRHS, _PMTRHS, _nVAR, _pVAR, _PMVAR );
  ++stats_traj.numTRHS_PM;

  // Stepsize selection and validation
  _h = _validation( q, _PMTRHS, _Ix );
  if( _h < options.HMIN ) return FAILURE;
  if( _h > options.HMAX ) _h = options.HMAX;
  // Modify stepsize if stage time is exceeded
  if( _t+_h > tnext ) _h = tnext-_t;
  ++stats_traj.numSteps;

  // Enclosure propagation
  for( unsigned int ix=0; ix<_nx; ix++ ){
    _PMx[ix] = _predictor( _h, q, _nx, _PMTRHS, ix );
    // Truncation error
#ifdef MC__ODEBND_VAL_STEPSIZE_HSTAB
    _PMx[ix] += _h*options.TOL*T(-1.,1.);
#else
    _PMx[ix] += options.TOL*T(-1.,1.);
#endif
#ifdef MC__ODEBND_VAL_DEBUG
    std::cout << "PMx(@" << _t+_h << ")[" << ix << "] = " << _PMx[ix] << std::endl;
#endif
    _Ix[ix] = _PMx[ix].bound();
  }
  _Ex.set();

  // Check enclosure diameter before returning
  return( _diam(_Ix) < options.DMAX? NORMAL: FAILURE );
}

template <typename T, typename PMT, typename PVT>
inline typename ODEBND_VAL<T,PMT,PVT>::STATUS
ODEBND_VAL<T,PMT,PVT>::_propagate_ELL0
( const unsigned q, const double tnext, std::ostream&os )
{
  // Stepsize validation
  if( options.TSTOP || isequal( _t, _tmax ) ){

    _Ex_TE = _Ex; // Copy of _Ex at expansion point for updates later on
    for( unsigned int ix=0; ix<_nx; ix++ ){
      _Rx_TE[ix] = _PMx[ix].remainder(); // Copy of polynomial model remainder at expansion time for updates later on
      _PMx[ix].center().set( T(0.) ); // Keep centered polynomial part
    }

    // Taylor series expansion in PVT arithmetic up to order TSORDER+1
    for( unsigned ix=0; ix<_nx; ix++ ) _PMVAR[ix] = _PMx[ix];
    _PMVAR[_nx+_np] = _t; // Keep time when stepsize was computed
    _pDAG->eval( _opTRHS, _PMDAG, (q+2)*_nx, _pTRHS, _PMTRHS, _nVAR, _pVAR, _PMVAR );
    ++stats_traj.numTRHS_PM;

    // Taylor series expansion and differentiation in T arithmetic up to order TSORDER
    for( unsigned ix=0; ix<_nx; ix++ ) _IVAR[ix] = _Ix[ix];
    _IVAR[_nx+_np] = _t;
    _pDAG->eval( _opFTRHS, _IDAG, (q+1)*_nx*_nx, _pFTRHS, _IFTRHS, _nVAR, _pVAR, _IVAR );
    ++stats_traj.numFTRHS_I;

    // Stepsize selection and validation
    _hmax = _h = _validation( q, _PMTRHS, _Ix );
    if( _hmax < options.HMIN ) return FAILURE;
    if( _hmax > options.HMAX ) _hmax = _h = options.HMAX;
    // Modify stepsize if stage time is exceeded
    _tmax = _t+_hmax;
    // Keep time when stepsize was computed
    _tval = _t; 
    ++stats_traj.numSteps;
#ifdef MC__ODEBND_VAL_DEBUG
    std::cout << "\nCurrent step: hmax=" << _hmax << "  tmax=" << _tmax << std::endl;
    { int dum; std::cin >> dum; }
#endif
  }

  // By-pass stepsize validation if previous integration step can be
  // carried out further and ODE RHS is continuous
  else{
    _h = _hmax;
  }

  // Limit stepsize so as to stop at stage time tnext
  if( _tmax > tnext ) _h = tnext-_tval; // Do not use _t directly as it may be different from the validation time
#ifdef MC__ODEBND_VAL_DEBUG
  std::cout << "\nCurrent step: t=" << _t << "  h=" << _h << "  tnext=" << tnext
            << "  hmax=" << _hmax << "  tmax=" << _tmax << std::endl;
#endif

  // Enclosure propagation
  for( unsigned int ix=0; ix<_nx; ix++ ){
    _PMx[ix] = _predictor( _h, q, _nx, _PMTRHS, ix ).center();
#ifdef MC__ODEBND_VAL_DEBUG
    std::cout << "_PMx[" << ix << "] =" << _PMx[ix];
#endif
    // Remainder of polynomial model propagation of state polynomial approximant
    _Rx[ix] = _PMx[ix].remainder();
    // Linear part for ellipsoidal remainder update
    for( unsigned int jx=0; jx<_nx; jx++ ){
      T IAx = _predictor( _h, q, _nx, _IFTRHS, ix, jx );
      _Ax(ix,jx) = Op<T>::mid( IAx );
      _Rx[ix] += ( IAx - _Ax(ix,jx) ) * _Rx_TE[jx];
    }
    // Truncation error
#ifdef MC__ODEBND_VAL_STEPSIZE_HSTAB
    _Rx[ix] += _h*options.TOL*T(-1.,1.);
#else
    _Rx[ix] += options.TOL*T(-1.,1.);
#endif
#ifdef MC__ODEBND_VAL_DEBUG
    std::cout << "_Rx[" << ix << "] =" << _Rx[ix]
              << "  diam(_Rx[" << ix << "]) =" << Op<T>::diam(_Rx[ix])
              << "  midp(_Rx[" << ix << "]) =" << Op<T>::mid(_Rx[ix])
              << std::endl;
#endif
  }
#ifdef MC__ODEBND_VAL_DEBUG
  std::cout << "_Ax(@" << _tval << ") =\n" << _Ax << std::endl;
  std::cout << "_Ex(@" << _tval << ") =" << _Ex_TE << std::endl;
#endif

  // Correct stepsize in case stepsize selection/validation was by-passed
  _h -= (_t-_tval); 
#ifdef MC__ODEBND_VAL_DEBUG
  std::cout << "\nAdjusted step: h=" << _h << "  tval=" << _tval << std::endl;
#endif

  return NORMAL;
}

template <typename T, typename PMT, typename PVT>
inline typename ODEBND_VAL<T,PMT,PVT>::STATUS
ODEBND_VAL<T,PMT,PVT>::_propagate_ELL1
( const unsigned q, const double tnext, std::ostream&os )
{
  // Stepsize validation
  if( options.TSTOP || isequal( _t, _tmax ) ){

    _Ex_TE = _Ex; // Copy of _Ex at expansion point for updates later on
    for( unsigned int ix=0; ix<_nx; ix++ ){
      _Rx_TE[ix] = _PMx[ix].remainder(); // Copy of polynomial model remainder at expansion time for updates later on
      _PMx[ix].center().set( T(0.) ); // get centered polynomial part
    }

    // Taylor series expansion in PVT arithmetic up to order TSORDER+1
    for( unsigned ix=0; ix<_nx; ix++ ) _PMVAR[ix] = _PMx[ix];
    _PMVAR[_nx+_np] = _t; // Keep time when stepsize was computed
    _pDAG->eval( _opTRHS, _PMDAG, (q+2)*_nx, _pTRHS, _PMTRHS, _nVAR, _pVAR, _PMVAR );
    ++stats_traj.numTRHS_PM;

    // Taylor series expansion and differentiation in PVT arithmetic up to order TSORDER
    for( unsigned ix=0; ix<_nx; ix++ ){
      PVT MVri( _MVenv, _np+ix, _Rx_TE[ix] );
      _MVVAR[ix].set( _MVenv ).set( _PMx[ix], true );
      _MVVAR[ix] += MVri;
    }
    _MVVAR[_nx+_np] = _t;
    _pDAG->eval( _opFTRHS, _PMDAG, (q+1)*_nx*_nx, _pFTRHS, _PMFTRHS, _nVAR, _pVAR, _MVVAR );
    ++stats_traj.numFTRHS_PM;

    // Stepsize selection and validation
    _hmax = _h = _validation( q, _PMTRHS, _Ix );
    if( _hmax < options.HMIN ) return FAILURE;
    if( _hmax > options.HMAX ) _hmax = _h = options.HMAX;
    // Modify stepsize if stage time is exceeded
    _tmax = _t+_hmax;
    // Keep time when stepsize was computed
    _tval = _t; 
    ++stats_traj.numSteps;
#ifdef MC__ODEBND_VAL_DEBUG
    std::cout << "\nCurrent step: hmax=" << _hmax << "  tmax=" << _tmax << std::endl;
    //{ int dum; std::cin >> dum; }
#endif
  }

  // By-pass stepsize validation if previous integration step can be
  // carried out further and ODE RHS is continuous
  else{
    _h = _hmax;
  }

  // Limit stepsize so as to stop at stage time tnext
  // Do not use _t since might be different from Taylor expansion time
  if( _tmax > tnext ) _h = tnext-_tval;
#ifdef MC__ODEBND_VAL_DEBUG
  std::cout << "\nCurrent step: t=" << _t << "  h=" << _h << "  tnext=" << tnext
            << "  hmax=" << _hmax << "  tmax=" << _tmax << std::endl;
#endif

  // Enclosure propagation
  for( unsigned int ix=0; ix<_nx; ix++ ){
    _PMx[ix] = _predictor( _h, q, _nx, _PMTRHS, ix ).center();
#ifdef MC__ODEBND_VAL_DEBUG
    std::cout << "_PMx[" << ix << "] =" << _PMx[ix];
#endif
    // Remainder of polynomial model propagation of state polynomial approximant
    _Rx[ix] = _PMx[ix].remainder();
    // Linear part for ellipsoidal remainder update
    for( unsigned int jx=0; jx<_nx; jx++ ){
      PVT MVfx = _predictor( _h, q, _nx, _PMFTRHS, ix, jx ).center();
      _Ax(ix,jx) = MVfx.constant();
      _Rx[ix] += ( MVfx.bound() - _Ax(ix,jx) ) * _Rx_TE[jx];
    }
    // Truncation error
#ifdef MC__ODEBND_VAL_STEPSIZE_HSTAB
    _Rx[ix] += _h*options.TOL*T(-1.,1.);
#else
    _Rx[ix] += options.TOL*T(-1.,1.);
#endif
#ifdef MC__ODEBND_VAL_DEBUG
    std::cout << "_Rx[" << ix << "] =" << _Rx[ix]
              << "  diam(_Rx[" << ix << "]) =" << Op<T>::diam(_Rx[ix])
              << "  midp(_Rx[" << ix << "]) =" << Op<T>::mid(_Rx[ix])
              << std::endl;
#endif
  }
#ifdef MC__ODEBND_VAL_DEBUG
  std::cout << "_Ax(@" << _tval << ") =\n" << _Ax << std::endl;
  std::cout << "_Ex(@" << _tval << ") =" << _Ex_TE << std::endl;
#endif

  // Correct stepsize in case stepsize selection/validation was by-passed
  _h -= (_t-_tval); 
#ifdef MC__ODEBND_VAL_DEBUG
  std::cout << "\nAdjusted step: h=" << _h << "  tval=" << _tval << std::endl;
#endif

  return NORMAL;
}
template <typename T, typename PMT, typename PVT>
inline typename ODEBND_VAL<T,PMT,PVT>::STATUS
ODEBND_VAL<T,PMT,PVT>::_propagate_ELL2
( const unsigned q, const double tnext, std::ostream&os )
{
  // Stepsize validation
  if( options.TSTOP || isequal( _t, _tmax ) ){

    _Ex_TE = _Ex; // Copy of _Ex at expansion point for updates later on
    for( unsigned int ix=0; ix<_nx; ix++ ){
      _Rx[ix] = _PMx[ix].center().remainder(); // get centered remainder
      _PMx[ix].set( T(0.) ); // cancel remainder term to get polynomial part
#ifdef MC__ODEBND_VAL_DEBUG
      std::cout << "_PMx[" << ix << "] =" << _PMx[ix];
      std::cout << "_Rx[" << ix << "] =" << _Rx[ix] << std::endl;
#endif
    }

    // Taylor series expansion in PVT arithmetic (joint parameter-state) up to order TSORDER+1
   for( unsigned ix=0; ix<_nx; ix++ ){
      PVT MVri( _MVenv, _np+ix, _Rx[ix] );
      _MVVAR[ix].set( _MVenv ).set( _PMx[ix], true );
      _MVVAR[ix] += MVri;
    }
    _MVVAR[_nx+_np] = _t;
    _pDAG->eval( _opTRHS, _PMDAG, (q+2)*_nx, _pTRHS, _PMTRHS, _nVAR, _pVAR, _MVVAR );
    ++stats_traj.numTRHS_PM;

    // Stepsize selection and validation
    _hmax = _h = _validation( q, _PMTRHS, _Ix );
    if( _hmax < options.HMIN ) return FAILURE;
    if( _hmax > options.HMAX ) _hmax = _h = options.HMAX;
    // Modify stepsize if stage time is exceeded
    _tmax = _t+_hmax;
    // Keep time when stepsize was computed
    _tval = _t; 
    ++stats_traj.numSteps;
  }

  // By-pass stepsize validation if previous integration step can be
  // carried out further and ODE RHS is continuous
  else{
#ifdef MC__ODEBND_VAL_DEBUG
    std::cout << "Expansion by-passed: t=" << _t << "  tmax=" << _tmax << "  hmax=" << _hmax << std::endl;
#endif
    _h = _hmax;
  }

  // Limit stepsize so as to stop at stage time tnext
  // Do not use _t since might be different from Taylor expansion time
  if( _tmax > tnext ) _h = tnext-_tval;
#ifdef MC__ODEBND_VAL_DEBUG
  std::cout << "\nCurrent step: t=" << _t << "  h=" << _h << "  tnext=" << tnext
            << "  hmax=" << _hmax << "  tmax=" << _tmax << std::endl;
#endif

  // Enclosure propagation
  for( unsigned int ix=0; ix<_nx; ix++ ){
    PVT MVXh = _predictor( _h, q, _nx, _PMTRHS, ix ).center();
#ifdef MC__ODEBND_VAL_DEBUG
    std::cout << "MVXh[" << ix << "] =" << MVXh;
#endif
    // Extract state polynomial approximant in p and set to 0
    MVXh.get( _PMx[ix].set(_PMenv), true );
#ifdef MC__ODEBND_VAL_DEBUG
    std::cout << "_PMx[" << ix << "] =" << _PMx[ix];
    std::cout << "MVXh[" << ix << "] =" << MVXh;
#endif
    // Extract first-order terms and set to 0
    for( unsigned int jx=0; jx<_nx; jx++ )
      _Ax(ix,jx) = MVXh.linear(_np+jx,true);
    // Remainder is higher-order terms *plus* truncation error
#ifdef MC__ODEBND_VAL_STEPSIZE_HSTAB
    _Rx[ix] = MVXh.bound() + _h*options.TOL*T(-1.,1.);
#else
    _Rx[ix] = MVXh.bound() + options.TOL*T(-1.,1.);
#endif
#ifdef MC__ODEBND_VAL_DEBUG
    std::cout << "_Rx[" << ix << "] =" << _Rx[ix]
              << "  diam(_Rx[" << ix << "]) =" << Op<T>::diam(_Rx[ix])
              << "  midp(_Rx[" << ix << "]) =" << Op<T>::mid(_Rx[ix])
              << std::endl;
#endif
  }
#ifdef MC__ODEBND_VAL_DEBUG
  std::cout << "_Ax(@" << _tval << ") =\n" << _Ax << std::endl;
  std::cout << "_Ex(@" << _tval << ") =" << _Ex_TE << std::endl;
#endif

  // Correct stepsize in case stepsize selection/validation was by-passed
  _h -= (_t-_tval); 
#ifdef MC__ODEBND_VAL_DEBUG
  std::cout << "\nAdjusted step: h=" << _h << "  tval=" << _tval << std::endl;
#endif

  return NORMAL;
}

template <typename T, typename PMT, typename PVT>
inline typename ODEBND_VAL<T,PMT,PVT>::STATUS
ODEBND_VAL<T,PMT,PVT>::_invariant_ELL0
( const PVT*Px, const T*Rx, std::ostream&os )
{
  // Compute invariants in PVT arithmetic
  for( unsigned i=0; i<_nx; i++ ) _PMVAR[i] = Px[i];
  _PMVAR[_nx+_np] = _t;
  _pDAG->eval( _opINV, _PMDAG, _ni, _pINV, _PMINV, _nVAR, _pVAR, _PMVAR );

  // Compute invariant derivatives in T arithmetic
  for( unsigned i=0; i<_nx; i++ ) _IVAR[i] = _PMx[i].bound() + _Rx[i];
  _IVAR[_nx+_np] = _t;
  _pDAG->eval( _opFINV, _IDAG, _ni*_nx, _pFINV, _IFINV, _nVAR, _pVAR, _IVAR );

  // Compute invariant directions and remainders
  for( unsigned i=0, ij=0; i<_ni; i++ ){
    _Rinv[i] = _PMINV[i].remainder();
    for( unsigned j=0; j<_nx; j++, ij++ ){
      _Ainv(i,j) = Op<T>::mid(_IFINV[ij]); // row-wise storage in IFINV
      _Rinv[i] += ( _IFINV[ij] - _Ainv(i,j) ) * Rx[j];
    }
#ifdef MC__ODEBND_VAL_DEBUG
    std::cout << "_Rinv[" << i << "] =" << _Rinv[i]
              << "  diam(_Rinv[" << i << "]) =" << Op<T>::diam(_Rinv[i])
              << "  midp(_Rinv[" << i << "]) =" << Op<T>::mid(_Rinv[i])
              << std::endl;
#endif
  }
#ifdef MC__ODEBND_VAL_DEBUG
  std::cout << "_Ainv(@" << _t << ") =\n" << _Ainv << std::endl;
  { int dum; std::cin >> dum; }
#endif

  return NORMAL;
}

template <typename T, typename PMT, typename PVT>
inline typename ODEBND_VAL<T,PMT,PVT>::STATUS
ODEBND_VAL<T,PMT,PVT>::_invariant_ELL1
( const PVT*Px, const T*Rx, std::ostream&os )
{
  // Compute invariants in PVT arithmetic
  for( unsigned i=0; i<_nx; i++ ) _PMVAR[i] = Px[i];
  _PMVAR[_nx+_np] = _t;
  _pDAG->eval( _opINV, _PMDAG, _ni, _pINV, _PMINV, _nVAR, _pVAR, _PMVAR );

  // Compute invariant derivatives in PVT arithmetic
  for( unsigned i=0; i<_nx; i++ ){
    PVT MVri( _MVenv, _np+i, Rx[i] );
    _MVVAR[i].set( _MVenv ).set( Px[i], true );
    _MVVAR[i] += MVri;
#ifdef MC__ODEBND_VAL_DEBUG
    std::cout << "Px[" << i << "] = " << Px[i] << std::endl;
    std::cout << "_MVx[" << i << "] = " << _MVVAR[i] << std::endl;
#endif
  }
  _MVVAR[_nx+_np] = _t;
  _pDAG->eval( _opFINV, _PMDAG, _ni*_nx, _pFINV, _PMFINV, _nVAR, _pVAR, _MVVAR );

  for( unsigned i=0, ij=0; i<_ni; i++ ){
    _Rinv[i] = _PMINV[i].remainder();
    for( unsigned j=0; j<_nx; j++, ij++ ){
#ifdef MC__ODEBND_VAL_DEBUG
      std::cout << "_PMFINV[" << i << "," << j << "] =" << _PMFINV[ij]
                << std::endl;
#endif
      _Ainv(i,j) = _PMFINV[ij].center().constant();
      _Rinv[i] += ( _PMFINV[ij].bound() - _Ainv(i,j) ) * Rx[j];
    }
#ifdef MC__ODEBND_VAL_DEBUG
    std::cout << "_Rinv[" << i << "] =" << _Rinv[i]
              << "  diam(_Rinv[" << i << "]) =" << Op<T>::diam(_Rinv[i])
              << "  midp(_Rinv[" << i << "]) =" << Op<T>::mid(_Rinv[i])
              << std::endl;
#endif
  }
#ifdef MC__ODEBND_VAL_DEBUG
  std::cout << "_Ainv(@" << _t << ") =\n" << _Ainv << std::endl;
  { int dum; std::cin >> dum; }
#endif

  return NORMAL;
}

template <typename T, typename PMT, typename PVT>
inline typename ODEBND_VAL<T,PMT,PVT>::STATUS
ODEBND_VAL<T,PMT,PVT>::_invariant_ELL2
( const PVT*Px, const T*Rx, std::ostream&os )
{
  // Compute invariants in PVT arithmetic (joint X,P expansion)
  for( unsigned i=0; i<_nx; i++ ){
    PVT MVri( _MVenv, _np+i, Rx[i] );
    _MVVAR[i].set( _MVenv ).set( Px[i], true );
    _MVVAR[i] += MVri;
#ifdef MC__ODEBND_VAL_DEBUG
    std::cout << "Px[" << i << "] = " << Px[i] << std::endl;
    std::cout << "_MVx[" << i << "] = " << _MVVAR[i] << std::endl;
#endif
  }
  _MVVAR[_nx+_np] = _t;
  _pDAG->eval( _opINV, _PMDAG, _ni, _pINV, _PMINV, _nVAR, _pVAR, _MVVAR );

  for( unsigned int i=0; i<_ni; i++ ){
    // Extract state polynomial approximant in p and set to 0
    PVT Pinv; _PMINV[i].get( Pinv.set(_PMenv), true ); //THIS SHOULD NORMALLY BE ZERO!!!
    // Extract first-order terms and set to 0
    for( unsigned int j=0; j<_nx; j++ )
      _Ainv(i,j) = _PMINV[i].linear( _np+j, true );
    _Rinv[i] = _PMINV[i].bound();
#ifdef MC__ODEBND_VAL_DEBUG
    std::cout << "_Rinv[" << i << "] =" << _Rinv[i]
              << "  diam(_Rinv[" << i << "]) =" << Op<T>::diam(_Rinv[i])
              << "  midp(_Rinv[" << i << "]) =" << Op<T>::mid(_Rinv[i])
              << std::endl;
#endif
  }
#ifdef MC__ODEBND_VAL_DEBUG
  std::cout << "_Ainv(@" << _t << ") =\n" << _Ainv << std::endl;
  { int dum; std::cin >> dum; }
#endif

  return NORMAL;
}

template <typename T, typename PMT, typename PVT>
inline typename ODEBND_VAL<T,PMT,PVT>::STATUS
ODEBND_VAL<T,PMT,PVT>::_invariant_ELL
( std::ostream&os )
{
  if( !options.USEINV || !_ni ) return NORMAL;
  
  // Update polynomial model and derived bounds
  for( unsigned int ix=0; ix<_nx; ix++ ){
    _PMx[ix].set( T(0.) ); // polynomial part
    _Rx[ix] = T( _Ex.l(ix), _Ex.u(ix) ); // centered remainder
  }

  // In this variant a bound on the Jacobian matrix of the invariant is computed
  // and the linear part is taken as the mid-point of this matrix
  if( !options.ORDMIT ){
    if( _invariant_ELL0( _PMx, _Rx, os ) != NORMAL ) return FAILURE;
  }

  // In this variant a polynomial model of the Jacobian matrix of the invariant is
  // computed and the linear part is taken as the mid-point of this matrix
  else if( options.ORDMIT < _PMenv->nord() ){
    if( _invariant_ELL1( _PMx, _Rx, os ) != NORMAL ) return FAILURE;
  }

  // In this variant a polynomial model of the invariant in the joint state-parameter
  // and of the same order as the parameter polynomial model is computed
  else{
    if( _invariant_ELL2( _PMx, _Rx, os ) != NORMAL ) return FAILURE;
  }

  // Contract ellipsoidal remainder
  for( unsigned int ii=0; ii<_ni; ii++ ){
    // Try exact intersection with hyperplane
    E Ex_red;
    try{
      Ex_red = hpintersection( _Ex, std::make_pair( t(_Ainv.row(ii)), 0. ) );
#ifdef MC__ODEBND_VAL_DEBUG_INVARIANT
      std::cout << "_Ex (intersection) =\n" << Ex_red << std::endl;
#endif
      // Nonlinear invariant case <-- NOT YET TESTED EXTENSIVELY-->
      if( Op<T>::diam(_Rinv[ii]) > machprec() ){
        for( unsigned int jx=0; jx<_nx; jx++ )
          _Ninv[jx] = _Ainv(ii,jx)/(_Ainv.row(ii)%_Ainv.row(ii))*_Rinv[ii];
        Ex_red = minksum_ea( Ex_red, _Ninv, options.QTOL );
#ifdef MC__ODEBND_VAL_DEBUG_INVARIANT
        std::cout << "_Ex (minkowski sum) =\n" << Ex_red << std::endl;
#endif
      }
#ifdef MC__ODEBND_VAL_DEBUG_INVARIANT
      int dum; std::cin >> dum;
#endif
    }
    catch( E::Exceptions &eObj){
      continue;
    }
    if( Ex_red.trQ() < _Ex.trQ() ){
#ifdef MC__ODEBND_VAL_DEBUG_INVARIANT
      std::cout << "_Ex (initial) =\n" << _Ex << std::endl;
      std::cout << "_Ex (updated) =\n" << Ex_red << std::endl;
      int dum; std::cin >> dum;
#endif
      _Ex = Ex_red;
    }
  }

  return NORMAL;
}

template <typename T, typename PMT, typename PVT>
inline typename ODEBND_VAL<T,PMT,PVT>::STATUS
ODEBND_VAL<T,PMT,PVT>::_update_ELL
( std::ostream&os )
{
  try{
    //_Ex = mtimes(_Ex_TE,_Ax);
    _Ex = minksum_ea( mtimes(_Ex_TE,_Ax), _Rx, options.QTOL );
    for( unsigned int ix=0; ix<_nx; ix++ ) _PMx[ix] += _Ex_TE.c(ix);
    _Ex = _Ex.O();
#ifdef MC__ODEBND_VAL_DEBUG
    std::cout << "_Ex(" << _t+_h << ") =" << _Ex << std::endl;
    {int dum; std::cin >> dum;}
#endif
  }
  catch( E::Exceptions &eObj){
    return FAILURE;
  }
#ifdef MC__ODEBND_VAL_DEBUG
  for( unsigned int ix=0; ix<_nx; ix++ )
    std::cout << "_Rx[" << ix << " =\n" << T( _Ex.l(ix), _Ex.u(ix) ) << std::endl;
#endif

  // Contract ellipsoidal remainder using ODE invariants
  if( _invariant_ELL( os ) != NORMAL ) return FAILURE;

  // Update polynomial model and derived bounds
  for( unsigned int ix=0; ix<_nx; ix++ ){
    _PMx[ix].set( T( _Ex.l(ix), _Ex.u(ix) ) );
    _Ix[ix] =_PMx[ix].bound();
#ifdef MC__ODEBND_VAL_DEBUG
    std::cout << "_Rx[" << ix << " =\n" << T( _Ex.l(ix), _Ex.u(ix) ) << std::endl;
    std::cout << "_Ix[" << ix << " (bef.) = " << _Ix[ix] << std::endl;
#endif
  }

  return NORMAL;
}

template <typename T, typename PMT, typename PVT>
inline typename ODEBND_VAL<T,PMT,PVT>::STATUS
ODEBND_VAL<T,PMT,PVT>::_propagate_ELL
( const double tnext, std::ostream&os )
{
  // In this variant a bound on the Jacobian matrix is computed and the
  // linear part is taken as the mid-point of this matrix
  if( !options.ORDMIT ){
    if( _propagate_ELL0( options.TSORDER, tnext, os ) != NORMAL ) return FAILURE;
  }
  // In this variant a polynomial model of the Jacobian matrix is computed and the
  // linear part is taken as the mid-point of this matrix
  else if( options.ORDMIT < _PMenv->nord() ){
    if( _propagate_ELL1( options.TSORDER, tnext, os ) != NORMAL ) return FAILURE;
  }
  // In this variant a polynomial model in the joint state-parameter and
  // of the same order as the parameter polynomial model is computed
  else{
    if( _propagate_ELL2( options.TSORDER, tnext, os ) != NORMAL ) return FAILURE;
  }

  // Compute ellispoidal remainder
  if( _update_ELL( os ) != NORMAL ) return FAILURE;

  // Check enclosure diameter before returning
  return( _diam(_Ix) < options.DMAX? NORMAL: FAILURE );
}

template <typename T, typename PMT, typename PVT>
inline typename ODEBND_VAL<T,PMT,PVT>::STATUS
ODEBND_VAL<T,PMT,PVT>::_propagate
( const double tnext, std::ostream&os )
{
  // Propagate validated enclosure until tnext
  STATUS flag = NORMAL;
  for( _h = 0; flag == NORMAL && tnext > _t; _t += _h ){
    switch( options.WRAPMIT){
      case Options::NONE:
        flag = _propagate_INT( tnext, os ); break;

      case Options::ELLIPS:
      default:
        flag = _propagate_ELL( tnext, os ); break;
    }
  }
  return flag;
}

template <typename T, typename PMT, typename PVT> inline bool
ODEBND_VAL<T,PMT,PVT>::_resize
( const unsigned q, const unsigned iRHS )
{
  if( _vRHS.size() <= iRHS ) return false;

  delete[] _pTRHS; _pTRHS = _pDAG->TAD( q+1, _nx, _vRHS[iRHS], _nx, _pVAR, _pT );
  _opTRHS = _pDAG->subgraph( (q+2)*_nx, _pTRHS );
  delete[] _PMTRHS; _PMTRHS = new PVT[(q+2)*_nx];
  unsigned nPMDAG = _opTRHS.size();

  _opTRHSval = _pDAG->subgraph( _nx, _pTRHS+(q+1)*_nx );
  delete[] _PMTRHSval; _PMTRHSval = options.PMVALID? new PVT[_nx]: 0;
  delete[] _ITRHSval; _ITRHSval = options.PMVALID? 0: new T[_nx];
  unsigned nIDAG = options.PMVALID? 0: _opTRHSval.size();

  delete[] _pFTRHS; _pFTRHS = 0;
  _opFTRHS.clear();
  delete[] _PMFTRHS; _PMFTRHS = 0;
  delete[] _IFTRHS;  _IFTRHS  = 0;

  _pINV = 0;
  _opINV.clear();
  delete[] _PMINV; _PMINV = 0;

  delete[] _pFINV; _pFINV = 0;
  _opFINV.clear();
  delete[] _PMFINV; _PMFINV = 0;
  delete[] _IFINV;  _IFINV  = 0;

  switch( options.WRAPMIT){
  case Options::NONE:
    break;

  case Options::ELLIPS: default:
    if( options.USEINV && _vINV.size() && _ni ){   // INV
      _pINV = _vINV[iRHS];
      _opINV = _pDAG->subgraph( _ni, _pINV );
      _PMINV = new PVT[_ni];
      if( nPMDAG < _opINV.size() ) nPMDAG = _opINV.size();
    }
    if( !options.ORDMIT ){                         // ELL0
      _pFTRHS = _pDAG->FAD( (q+1)*_nx, _pTRHS, _nx, _pVAR );
      _opFTRHS = _pDAG->subgraph( (q+1)*_nx*_nx, _pFTRHS );
      _IFTRHS  = new T[(q+1)*_nx*_nx];
      if( nIDAG < _opFTRHS.size() ) nIDAG = _opFTRHS.size();
      if( options.USEINV && _vINV.size() && _ni ){
        _pFINV = _pDAG->FAD( _ni, _pINV, _nx, _pVAR );
        _opFINV = _pDAG->subgraph( _ni*_nx, _pFINV );
        _IFINV  = new T[_ni*_nx];
        if( nIDAG < _opFINV.size() ) nIDAG = _opFINV.size();
      }
    }
    else if( _PMenv->nord() > _MVenv->nord() ){    // ELL1
      _pFTRHS = _pDAG->FAD( (q+1)*_nx, _pTRHS, _nx, _pVAR );
      _opFTRHS = _pDAG->subgraph( (q+1)*_nx*_nx, _pFTRHS );
      _PMFTRHS = new PVT[(q+1)*_nx*_nx];
      if( nPMDAG < _opFTRHS.size() ) nPMDAG = _opFTRHS.size();
      if( options.USEINV && _vINV.size() && _ni ){
        _pFINV = _pDAG->FAD( _ni, _pINV, _nx, _pVAR );
        _opFINV = _pDAG->subgraph( _ni*_nx, _pFINV );
        _PMFINV = new PVT[_ni*_nx];
        if( nPMDAG < _opFINV.size() ) nPMDAG = _opFINV.size();
      }
    }
                                                   // ELL2
    break;
  }

  delete[] _PMDAG; _PMDAG = new PVT[nPMDAG];
  delete[] _IDAG;  _IDAG  = new T[nIDAG];

  return true;
}

template <typename T, typename PMT, typename PVT>
inline bool
ODEBND_VAL<T,PMT,PVT>::_prepare
( const PVT* PMp )
{
  // Check polynomial model compatibility and size
  unsigned int kp=_np;
  for( unsigned int ip=0; ip<_np && kp==_np; ip++ )
    if( PMp[ip].env() ) kp = ip;
  if( kp==_np || PMp[kp].env()->nvar()!=_np ) return false;
  _PMenv = PMp[kp].env();  

  // Internal polynomial model environment reset
  unsigned int MVord = ( options.ORDMIT<_PMenv->nord()? options.ORDMIT: _PMenv->nord() ); 
  if( !_MVenv || _MVenv->nord() != MVord || _MVenv->nvar() != _nx+_np ){
    delete _MVenv; _MVenv = new PMT( _nx+_np, MVord );
    _MVenv->options = _PMenv->options;
  }

  // Size and set DAG evaluation arrays
  if( _nVAR != _nx+_np+1 ){
    _nVAR = _nx+_np+1;
    delete[] _pVAR; _pVAR = new FFVar[_nVAR];
    delete[] _IVAR; _IVAR = new T[_nVAR];
    delete[] _PMVAR; _PMVAR = new PVT[_nVAR];
    delete[] _MVVAR; _MVVAR = new PVT[_nVAR];
  }
  for( unsigned ix=0; ix<_nx; ix++ )
    _pVAR[ix] = _pX[ix];
  for( unsigned ip=0; ip<_np; ip++ ){
    _pVAR[_nx+ip] = _pP[ip];
    _IVAR[_nx+ip]  = PMp[ip].bound();
    _PMVAR[_nx+ip] = PMp[ip];
    _MVVAR[_nx+ip].set( _MVenv, ip, _IVAR[_nx+ip] );
  }
  _pVAR[_nx+_np] = (_pT? *_pT: 0. );

  // State bounds
  if( !_Ix )  _Ix  = new T[_nx];
  if( !_PMx ) _PMx = new PVT[_nx];

  // Size environments and arrays
  switch( options.WRAPMIT){
  case Options::NONE:
    break;

  case Options::ELLIPS:
  default:
    if( _Ax.n != _nx ){
      _Ax.resize(_nx,_nx);
      delete[] _Rx;    _Rx    = new T[_nx];
      delete[] _Rx_TE; _Rx_TE = new T[_nx];
    }
    if( options.USEINV && _ni && ( _Ainv.n != _nx || _Ainv.m != _ni ) ){
      _Ainv.resize(_ni,_nx);
      delete[] _Rinv;  _Rinv = new T[_ni];
      delete[] _Ninv;  _Ninv = new T[_ni];
    }
    break;
  }

  return true;
}

template <typename T, typename PMT, typename PVT>
inline bool
ODEBND_VAL<T,PMT,PVT>::_init
( const unsigned k, const double tk, PVT*PMxk, E*Exk )
{
  if( !_vIC.size() ) return false;
  if( k && _vIC.size()>1 ) throw Exceptions( Exceptions::REINIT );
  if( k ) return true;

  _t = _tmax = tk;
  _pIC = _vIC.at(0);
  _opIC = _pDAG->subgraph( _nx, _pIC );
 
  switch( options.WRAPMIT){

  case Options::NONE:
    _pDAG->eval( _opIC, _nx, _pIC, _PMx, _np, _pVAR+_nx, _PMVAR+_nx );
    for( unsigned int ix=0; ix<_nx; ix++ )
      _Ix[ix] = _PMx[ix].center().bound();
    _Ex.set();
    break;

  case Options::ELLIPS:
  default:
    _pDAG->eval( _opIC, _nx, _pIC, _PMx, _np, _pVAR+_nx, _PMVAR+_nx );
    for( unsigned int ix=0; ix<_nx; ix++ ){
      _PMx[ix].center();
      _Rx[ix] = _PMx[ix].remainder();
      _Ix[ix] = _PMx[ix].bound();
#ifdef MC__ODEBND_VAL_DEBUG
      std::cout << "PMx0[" << ix << "] =" << _PMx[ix] << std::endl;
#endif
    }
    _Ex.set( _nx, _Rx );
#ifdef MC__ODEBND_VAL_DEBUG
      std::cout << "Ex0 =" << _Ex << std::endl;
#endif
    // THERE MIGHT BE A BETTER WAY OF INITIALIZING IF AN ELLIPSOIDAL
    // REMAINDER BOUND IS KNOWN FOR THE PARAMETERS
    break;
  }

  // Copy internal variables back into arguments
  for( unsigned ix=0; ix<_nx; ix++ ) PMxk[ix] = _PMx[ix];
  if( Exk ) *Exk = _Ex;
  return true;
}

//! @fn template <typename T, typename PMT, typename PVT> inline typename ODEBND_VAL<T,PMT,PVT>::STATUS ODEBND_VAL<T,PMT,PVT>::bounds(
//! const unsigned int ntk, const double*tk, const PVT*PMp, PVT**PMxk,
//! E*Exk=0, std::ostream&os=std::cout )
//!
//! This function computes an interval enclosure of the reachable set of 
//! the parametric ODEs defined in IVP using equally spaced samples:
//!   - <a>ntk</a> [input] number of time stages
//!   - <a>tk</a>  [input] stage times, including the initial time
//!   - <a>PMp</a> [input] polynomial model parameter enclosure
//!   - <a>PMxk</a> [output] polynomial model state enclosures at stage times
//!   - <a>Exk</a> [output] ellipsoidal remainders in polynomial model state enclosures at stage times (default: NULL)
//!   - <a>os</a>  [input/output] output stream (default: std::cout)
//! .
//! The return value is the status.
template <typename T, typename PMT, typename PVT>
inline typename ODEBND_VAL<T,PMT,PVT>::STATUS
ODEBND_VAL<T,PMT,PVT>::bounds
( const unsigned int ntk, const double*tk, const PVT*PMp, PVT**PMxk,
  E*Exk, std::ostream&os )
{

  // Initialize validated ODE bounding
  if( !_prepare( PMp ) ) return FATAL;
  _results.clear();  
  _init_stats( stats_traj );

  try{
    // Initial state bounds
    if( !_init( 0, tk[0], PMxk[0], Exk ) ) return FATAL;
    if( options.DISPLAY >= 1 )
      switch( options.WRAPMIT ){
        case Options::NONE:   _print_interm( tk[0], _PMx, "x", os );      break;
        case Options::ELLIPS: _print_interm( tk[0], _PMx, _Ex, "x", os ); break;
      }
    // Record initial results
    if( options.RESRECORD )
      _results.push_back( Results( tk[0], _nx, _PMx ) );

    // Integrate ODEs through each stage
    for( _istg=0; _istg<ntk; _istg++ ){

      // update list of operations in RHS and INV
      const unsigned int iRHS = ( _vRHS.size()<2? 0: _istg );
      if( (!_istg || iRHS) && !_resize( options.TSORDER, iRHS ) ) return FATAL;

      if( _istg && !_init( _istg, tk[_istg], PMxk[_istg], Exk+_istg ) ) return FATAL;
      if( _propagate( tk[_istg+1], os ) != NORMAL ){
        _final_stats( stats_traj );
        if( options.DISPLAY >= 1 ) _print_stats( stats_traj, os );
	return FAILURE;
      }
      for( unsigned int ix=0; ix<_nx; ix++ ) PMxk[_istg+1][ix] = _PMx[ix];
      if( Exk ) Exk[_istg+1] = _Ex;
      if( options.DISPLAY >= 1 )
        switch( options.WRAPMIT ){
          case Options::NONE:   _print_interm( tk[_istg+1], _PMx, "x", os );      break;
          case Options::ELLIPS: _print_interm( tk[_istg+1], _PMx, _Ex, "x", os ); break;
        }
      // Record intermediate results
      if( options.RESRECORD ) _results.push_back( Results( tk[_istg+1], _nx, _PMx ) );
    }
  }
  catch(...){
    _final_stats( stats_traj );
    if( options.DISPLAY >= 1 ) _print_stats( stats_traj, os );
    return FAILURE;
  }

  _final_stats( stats_traj );
  if( options.DISPLAY >= 1 ) _print_stats( stats_traj, os );
  return NORMAL;
}

template <typename T, typename PMT, typename PVT> inline typename ODEBND_VAL<T,PMT,PVT>::STATUS
ODEBND_VAL<T,PMT,PVT>::_remainders
( const unsigned int ns, const double*tk, const T*Ip, const PVT*const*PMxk,
  T**Rxk, const unsigned int nsamp, std::ostream&os )
{
  STATUS flag = NORMAL;
  
  // Initialization of sampled bounds at parameter lower bound
  double *p = new double[_np];
  for( unsigned int ip=0; ip<_np; ip++ )
    p[ip] = Op<T>::l(Ip[ip]);
  double **xk = new double*[ns+1];
  for( unsigned int is=0; is<=ns; is++ )
    xk[is] = new double[_nx];
  flag = _ODESLV_GSL->states( ns, tk, p, xk, 0, os );
  if( flag != NORMAL || nsamp <= 1 ){
    delete[] p;
    for( unsigned int is=0; is<=ns; is++ ) delete[] xk[is]; delete[] xk;
    return flag;
  }   
  for( unsigned int is=0; is<=ns; is++ )
    for( unsigned int ix=0; ix<_nx; ix++ )
      Rxk[is][ix] = xk[is][ix] - PMxk[is][ix].polynomial( p );
  
  // Start sampling process
  unsigned int* vsamp = new unsigned int[_np];
  flag = _remainders( ns, tk, Ip, PMxk, Rxk, nsamp, vsamp, 0, p, xk, os );
  
  // Clean-up
  delete[] p;
  for( unsigned int is=0; is<=ns; is++ ) delete[] xk[is]; delete[] xk;
  delete[] vsamp;
  
  return flag;
}

template <typename T, typename PMT, typename PVT> inline typename ODEBND_VAL<T,PMT,PVT>::STATUS
ODEBND_VAL<T,PMT,PVT>::_remainders
( const unsigned int ns, const double*tk, const T*Ip, const PVT*const*PMxk,
  T**Rxk, const unsigned int nsamp, unsigned int* vsamp, const unsigned int ip,
  double*p, double**xk, std::ostream&os )
{
  STATUS flag = NORMAL;

  // Update bounds for all sampling points
  for( unsigned int isamp=0; isamp<nsamp; isamp++ ){
    vsamp[ip] = isamp;

    // Continue recursive call
    if( ip+1 < _np ){
      flag = _remainders( ns, tk, Ip, PMxk, Rxk, nsamp, vsamp, ip+1, p, xk, os );
      if( flag != NORMAL ) return flag;
      continue;
    }

    // Update bounds for current point
#ifdef MC__ODEBND_GSL_SAMPLE_DEBUG
    std::cout << "Sample: ";
#endif
    for( unsigned int ip=0; ip<_np; ip++ ){
      p[ip] = Op<T>::l( Ip[ip] ) + vsamp[ip]/(nsamp-1.) * Op<T>::diam( Ip[ip] );
#ifdef MC__ODEBND_GSL_SAMPLE_DEBUG
      std::cout << p[ip] << "  ";
#endif
    }
#ifdef MC__ODEBND_GSL_SAMPLE_DEBUG
    std::cout << std::endl;
#endif
    flag = _ODESLV_GSL->states( ns, tk, p, xk, 0, os );
    if( flag != NORMAL ) return flag;
    for( unsigned int is=0; is<=ns; is++ )
      for( unsigned int ix=0; ix<_nx; ix++ )
        Rxk[is][ix] = Op<T>::hull( xk[is][ix]-PMxk[is][ix].polynomial(p),
                                   Rxk[is][ix] );
  }

  return flag;
}  

//! @fn template <typename T, typename PMT, typename PVT> inline typename ODEBND_VAL<T,PMT,PVT>::STATUS ODEBND_VAL<T,PMT,PVT>::hausdorff(
//! const unsigned int ns, const double*tk, const PVT*PMp, double**Hxk,
//! const unsigned int nsamp, std::ostream&os )
//!
//! This function computes the Hausdorff distance between the polynomial model
//! remainder and the actual (sampled) range of the remainder function
//! in projection onto each variable and for each stage time
//! remainder and the actual range of the remainder function:
//!   - <a>ns</a> [input] number of time stages
//!   - <a>tk</a> [input] stage times, including the initial time
//!   - <a>PMp</a> [input] polynomial model of parameter set
//!   - <a>Hxk</a> [output] Hausdorff distance between the polynomial model
//!     remainder and the actual (sampled) range of the remainder function,
//!     at stage times
//!   - <a>nsamp</a> [input] number of samples for each parameter
//!   - <a>os</a> [input] output stream (default: std::cout)

template <typename T, typename PMT, typename PVT> inline typename ODEBND_VAL<T,PMT,PVT>::STATUS
ODEBND_VAL<T,PMT,PVT>::hausdorff
( const unsigned int ns, const double*tk, const PVT*PMp, double**Hxk,
  const unsigned int nsamp, std::ostream&os )
{
  int DISPLAY_ODEBND = options.DISPLAY;
  int DISPLAY_ODESLV = options.ODESLVOPT.DISPLAY;
  options.DISPLAY = 0;
  options.ODESLVOPT.DISPLAY = 0;
  _ODESLV_GSL->set( *this ); // in base class BASE_ODE
  _ODESLV_GSL->options = options.ODESLVOPT;
  _ODESLV_GSL->options.RESRECORD = false;

  // Compute approximate bounds
  PVT** PMxk = new PVT*[ns+1];
  for( unsigned int is=0; is<ns+1; is++ ) PMxk[is] = new PVT[_nx];
  try{ bounds( ns, tk, PMp, PMxk, 0, os ); }
  catch(...){;}
  unsigned int nsf = _istg;
  for( unsigned int is=nsf; is<ns; is++ )
    for( unsigned int ix=0; ix<_nx; ix++ )
      PMxk[is+1][ix] = T( -1e20, 1e20 );

  // Compute remainder bounds 
  T* Ip = new T[_np];
  for( unsigned int ip=0; ip<_np; ip++ ) Ip[ip] = PMp[ip].B();
  T** Rxk = new T*[ns+1];
  for( unsigned int is=0; is<ns+1; is++ ) Rxk[is] = new T[_nx];
  _remainders( nsf, tk, Ip, PMxk, Rxk, nsamp, os );
  for( unsigned int is=nsf; is<ns; is++ )
    for( unsigned int ix=0; ix<_nx; ix++ )
      Rxk[is+1][ix] = T( 0. );

  // Compute Hausdorff metric at each time step
  struct loc{ static double hausdorff
    ( const T&Ix, const T&Ix0 )
    { return std::max( std::fabs(Op<T>::l(Ix)-Op<T>::l(Ix0)),
                       std::fabs(Op<T>::u(Ix)-Op<T>::u(Ix0)) ); }
  };

  options.DISPLAY = DISPLAY_ODEBND;
  options.ODESLVOPT.DISPLAY = DISPLAY_ODESLV;
  for( unsigned int is=0; is<ns+1; is++ ){
    for( unsigned int ix=0; ix<_nx; ix++ )
      Hxk[is][ix] = loc::hausdorff( PMxk[is][ix].R(), Rxk[is][ix] );
    if( options.DISPLAY >= 1 ){
      _print_interm( tk[is], Hxk[is], "dH", os );
    }
  }

  for( unsigned int is=0; is<ns+1; is++ ) delete[] Rxk[is];
  delete[] Rxk;
  for( unsigned int is=0; is<ns+1; is++ ) delete[] PMxk[is];
  delete[] PMxk;
  delete[] Ip;

  return NORMAL;
}

template <typename T, typename PMT, typename PVT>
template <typename U> inline void
ODEBND_VAL<T,PMT,PVT>::_print_interm
( const double t, const U*x, const std::string&var, std::ostream&os ) const
{
  os << " @t = " << std::scientific << std::setprecision(4)
                 << std::left << t << " :" << std::endl;
  for( unsigned int ix=0; ix<_nx; ix++ )
    os << " " << var.c_str() << "[" << ix << "] = " << x[ix] << std::endl;
  return;
}

template <typename T, typename PMT, typename PVT>
template <typename U, typename V> inline void
ODEBND_VAL<T,PMT,PVT>::_print_interm
( const double t, const U*x, const V&r, const std::string&var, std::ostream&os ) const
{
  os << " @t = " << std::scientific << std::setprecision(4)
                 << std::left << t << " :" << std::endl;
  for( unsigned int ix=0; ix<_nx; ix++ )
    os << " " << var.c_str() << "[" << ix << "] = " << x[ix] << std::endl;
  os << " " << "R" << var.c_str() << " =" << r << std::endl;
  return;
}

template <typename T, typename PMT, typename PVT>
inline void
ODEBND_VAL<T,PMT,PVT>::record
( std::ofstream&bndrec, const unsigned int iprec ) const
{
  if( !bndrec ) return;

  // Specify format
  bndrec << std::right << std::scientific << std::setprecision(iprec);

  // Record computed interval bounds at stage times
  typename std::vector< Results >::const_iterator it = _results.begin();
  for( ; it != _results.end(); ++it ){
    bndrec << std::setw(iprec+9) << (*it).t;
    for( unsigned int ix=0; ix<_nx; ix++ )
      bndrec << std::setw(iprec+9) << mc::Op<T>::l( (*it).X[ix] )
             << std::setw(iprec+9) << mc::Op<T>::u( (*it).X[ix] );
    bndrec << std::endl;
  }
}

template <typename T, typename PMT, typename PVT>
inline void
ODEBND_VAL<T,PMT,PVT>::_init_stats
( Stats&stats )
{
  // Initialize statistics
  stats.reset();
  timeval time;
  gettimeofday(&time, 0) ;
  stats.cputime = - time.tv_sec - time.tv_usec*1e-6;
}

template <typename T, typename PMT, typename PVT>
inline void
ODEBND_VAL<T,PMT,PVT>::_final_stats
( Stats&stats )
{
  // Get final CPU time
  timeval time;
  gettimeofday(&time, 0);
  stats.cputime += time.tv_sec + time.tv_usec*1e-6;
}

template <typename T, typename PMT, typename PVT>
inline void
ODEBND_VAL<T,PMT,PVT>::_print_stats
( const Stats&stats, std::ostream&os )
{
  // Statistics
  os << " No STEPS  " << stats.numSteps
     << std::endl
     << " No EVALATIONS"
     << "   TRHS:  " << stats.numTRHS_I << " (IA)  " << stats.numTRHS_PM << " (PM)"
     << "   FTRHS: " << stats.numFTRHS_I << " (IA)  " << stats.numFTRHS_PM << " (PM)"
     << std::endl
     << " CPU TIME (SEC)     " << std::fixed << std::left
     << std::setprecision(5) << stats.cputime  << std::endl;
  return;
}

} // end namescape mc
#endif
