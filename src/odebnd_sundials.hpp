// Copyright (C) 2015 Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

#ifndef MC__ODEBND_SUNDIALS_HPP
#define MC__ODEBND_SUNDIALS_HPP

#undef  MC__ODEBND_SUNDIALS_DINEQI_DEBUG
#undef  MC__ODEBND_SUNDIALS_DINEQPM_DEBUG

#include "odebnd_base.hpp"
#include "base_sundials.hpp"

#include "tmodel.hpp"
#include "cmodel.hpp"

namespace mc
{
//! @brief C++ class computing enclosures of the reachable set of parametric ODEs using continuous-time set-valued integration and SUNDIALS-CVODE code.
////////////////////////////////////////////////////////////////////////
//! mc::ODEBND_SUNDIALS is a C++ class that computes enclosures of the
//! reachable set of parametric ordinary differential equations
//! (ODEs) using continuous-time set-valued integration. It implements
//! the methods of differential inequalities, whereby polynomial models
//! with interval or ellipsoidal remainders are used to enable high-
//! order convergence. The use of ellipsoidal remainders enables
//! stability of the enclosures for asymptotically stable ODE systems
//! when the parameter host is sufficiently small. The numerical
//! integrator is CVODE in SUNDIALS.
////////////////////////////////////////////////////////////////////////
template <typename T, typename PMT=mc::TModel<T>, typename PVT=mc::TVar<T> >
class ODEBND_SUNDIALS: public ODEBND_BASE<T,PMT,PVT>, public virtual BASE_SUNDIALS
{
 typedef Ellipsoid E;
 typedef BASE_DE::STATUS STATUS;
 using ODEBND_BASE<T,PMT,PVT>::NORMAL;
 using ODEBND_BASE<T,PMT,PVT>::FAILURE;
 using ODEBND_BASE<T,PMT,PVT>::FATAL;

 using ODEBND_BASE<T,PMT,PVT>::_nx;
 using ODEBND_BASE<T,PMT,PVT>::_np;
 using ODEBND_BASE<T,PMT,PVT>::_nq;
 using ODEBND_BASE<T,PMT,PVT>::_nf;
 using ODEBND_BASE<T,PMT,PVT>::_vRHS;
 using ODEBND_BASE<T,PMT,PVT>::_vQUAD;
 using ODEBND_BASE<T,PMT,PVT>::_vIC;
 using ODEBND_BASE<T,PMT,PVT>::_vFCT;

 using ODEBND_BASE<T,PMT,PVT>::_Q;
 using ODEBND_BASE<T,PMT,PVT>::_Er;
 using ODEBND_BASE<T,PMT,PVT>::_Ir;
 using ODEBND_BASE<T,PMT,PVT>::_pref;
 using ODEBND_BASE<T,PMT,PVT>::_Ip;
 using ODEBND_BASE<T,PMT,PVT>::_B;
 using ODEBND_BASE<T,PMT,PVT>::_xref;
 using ODEBND_BASE<T,PMT,PVT>::_Ix;
 using ODEBND_BASE<T,PMT,PVT>::_Iq;
 using ODEBND_BASE<T,PMT,PVT>::_vec2I;
 using ODEBND_BASE<T,PMT,PVT>::_vec2E;
 using ODEBND_BASE<T,PMT,PVT>::_IC_I_STA;
 using ODEBND_BASE<T,PMT,PVT>::_IC_I_QUAD;
 using ODEBND_BASE<T,PMT,PVT>::_CC_I_STA;
 using ODEBND_BASE<T,PMT,PVT>::_SET_I_STA;
 using ODEBND_BASE<T,PMT,PVT>::_RHS_I_STA;
 using ODEBND_BASE<T,PMT,PVT>::_RHS_I_QUAD;
 using ODEBND_BASE<T,PMT,PVT>::_JAC_I_STA;
 using ODEBND_BASE<T,PMT,PVT>::_FCT_I_STA;

 using ODEBND_BASE<T,PMT,PVT>::_PMenv;
 using ODEBND_BASE<T,PMT,PVT>::_PMx;
 using ODEBND_BASE<T,PMT,PVT>::_PMp;
 using ODEBND_BASE<T,PMT,PVT>::_PMq;
 using ODEBND_BASE<T,PMT,PVT>::_vec2PMI;
 using ODEBND_BASE<T,PMT,PVT>::_vec2PME;
 using ODEBND_BASE<T,PMT,PVT>::_IC_PM_STA;
 using ODEBND_BASE<T,PMT,PVT>::_IC_PM_QUAD;
 using ODEBND_BASE<T,PMT,PVT>::_CC_PM_STA;
 using ODEBND_BASE<T,PMT,PVT>::_SET_PM_STA;
 using ODEBND_BASE<T,PMT,PVT>::_RHS_PM_STA;
 using ODEBND_BASE<T,PMT,PVT>::_RHS_PM_QUAD;
 using ODEBND_BASE<T,PMT,PVT>::_JAC_PM_STA;
 using ODEBND_BASE<T,PMT,PVT>::_FCT_PM_STA;

 using ODEBND_BASE<T,PMT,PVT>::_diam;
 using ODEBND_BASE<T,PMT,PVT>::_dH;
 using ODEBND_BASE<T,PMT,PVT>::_hausdorff;
 using ODEBND_BASE<T,PMT,PVT>::_remainders;
 using ODEBND_BASE<T,PMT,PVT>::_print_interm;

typedef int (*CVRhsFn)( realtype t, N_Vector y, N_Vector ydot, void *user_data );

 private:
  //! @brief Pointer to the CVODE memory block
  void *_cv_mem;

  //! @brief Return flag for SUNDIALS methods
  int _cv_flag;

 protected:
  //! @brief N_Vector object holding current states
  N_Vector _Nx;

  //! @brief N_Vector object holding current state quadratures
  N_Vector _Nq;

  //! @brief stepsize
  double _h;

  //! @brief position in _vIC
  unsigned _pos_ic;

  //! @brief position in _vRHS
  unsigned _pos_rhs;

  //! @brief position in _vQUAD
  unsigned _pos_quad;
  
  //! @brief position in _vFCT
  unsigned _pos_fct;

  //! @brief Static const member for mesh type
  static const int PMOFFSET = 10;

  //! @brief static pointer to class
  static ODEBND_SUNDIALS<T,PMT,PVT> *_pODEBND;

public:
  /** @defgroup ODEBND_SUNDIALS Continuous-time set-valued integration of parametric ODEs
   *  @{
   */
  //! @brief Default constructor
  ODEBND_SUNDIALS
    ();

  //! @brief Virtual destructor
  virtual ~ODEBND_SUNDIALS
    ();

  //! @brief Integrator options
  struct Options: public BASE_SUNDIALS::Options
  {
    //! @brief Constructor
    Options():
      BASE_SUNDIALS::Options(), WRAPMIT(ELLIPS), ORDMIT(2), PMOPT(typename PMT::Options()),
      QTOL(machprec()), PMNOREM(false), USEINV(true), DMAX(1e20), DISPLAY(1), 
      RESRECORD(false)
      {}
    //! @brief Assignment operator
    Options& operator=
      ( Options&options ){
        BASE_SUNDIALS::Options::operator=(options);
        WRAPMIT   = options.WRAPMIT;
        ORDMIT    = options.ORDMIT;
        PMOPT     = options.PMOPT;
        QTOL      = options.QTOL;
        PMNOREM   = options.PMNOREM;
        USEINV    = options.USEINV;
        DMAX      = options.DMAX;
        DISPLAY   = options.DISPLAY;
        RESRECORD = options.RESRECORD;
        return *this;
      }
    //! @brief Enumeration of wrapping mitigation strategies
    enum WRAPPING_STRATEGY{
      NONE=0,		//!< No wrapping mitigation
      DINEQ,		//!< Differential inequality contractor
      ELLIPS		//!< Ellipsoidal contractor with linear preconditioning [Default]
    };
    //! @brief Wrapping mitigation strategy
    WRAPPING_STRATEGY WRAPMIT;
    //! @brief Order of wrapping mitigation strategy (Default: 2)
    unsigned ORDMIT;
    //! @brief Options of internal polynomial model arithmetic for wrapping mitigation
    typename PMT::Options PMOPT;
    //! @brief Tolerance when dividing by trace of shape matrix in ellipsoidal bounds (Default: machprec())
    double QTOL;
    //! @brief Whether or not to cancel the remainder term in polynomial models (non-validated integration; Default: false)
    bool PMNOREM;
    //! @brief Whether or not to use the specified invariants for bounds contraction (default: true)
    bool USEINV;
    //! @brief Maximum enclosure diameter, \f$D_{\rm max}\f$ (Default: 1e20)
    double DMAX;
    //! @brief Display level (default: 1)
    int DISPLAY;
    //! @brief Whether or not to record results (default: false)
    bool RESRECORD;
  } options;

  //! @brief Structure for setting up storing the solver exceptions
  class Exceptions
  {
  public:
    //! @brief Enumeration type for ODEBND_SUNDIALS exception handling
    enum TYPE{
      INTERN=-33	//!< Internal error
    };
    //! @brief Constructor for error <a>ierr</a>
    Exceptions( TYPE ierr ) : _ierr( ierr ){}
    //! @brief Inline function returning the error flag
    int ierr(){ return _ierr; }
    //! @brief Inline function returning the error description
    std::string what(){
      switch( _ierr ){
      case INTERN: default:
        return "ODEBND_SUNDIALS::Exceptions  Internal error";
       }
    }
  private:
    TYPE _ierr;
  };

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

  //! @brief Vector storing interval bound results (upon request only)
  std::vector< Results > results_sta;

  //! @brief Statistics for state bounds integration
  Stats stats_sta;

  //! @brief Computes interval enclosure of reachable set of parametric ODEs
  STATUS bounds
    ( const unsigned ns, const double*tk, const T*Ip, T**Ixk=0,
      T*Iq=0, T*If=0, std::ostream&os=std::cout );

  //! @brief Computes Hausdorff distance between interval enclosure and actual reachable set of parametric ODEs, using parameter sampling
  template <typename ODESLV> STATUS hausdorff
    ( const unsigned ns, const double*tk, const T*Ip, double**Hxk,
      ODESLV&traj, const unsigned nsamp, std::ostream&os=std::cout );

  //! @brief Computes polynomial model enclosure of reachable set of parametric ODEs
  STATUS bounds
    ( const unsigned ns, const double*tk, const PVT*PMp, PVT**PMxk=0,
      PVT*PMq=0, PVT*PMf=0, std::ostream&os=std::cout );

  //! @brief Computes Hausdorff distance between polynomial model remainder enclosure and actual remainder function range, using parameter sampling
  template <typename ODESLV> STATUS hausdorff
    ( const unsigned ns, const double*tk, const PVT*PMp, double**Hxk,
      ODESLV&traj, const unsigned nsamp, std::ostream&os=std::cout );

  //! @brief Record results in file <a>bndrec</a>, with accuracy of <a>iprec</a> digits
  void record
    ( std::ofstream&bndrec, const unsigned iprec=5 ) const
    { return ODEBND_BASE<T,PMT,PVT>::_record( bndrec, results_sta, iprec ); }
  /** @} */

protected:

  //! @brief Function to initialize CVode memory block
  bool _INI_CVODE
    ( CVRhsFn MC_CVRHS, CVRhsFn MC_CVQUAD );

  //! @brief Function to reinitialize CVode memory block
  bool _CC_CVODE
    ();

  //! @brief Function to finalize state bounding
  void _END_STA
    ();

  //! @brief Function to initialize state interval bounding
  bool _INI_I_STA
    ( const unsigned np, const T*Ip, const unsigned ns );

  //! @brief Static wrapper to function computing the DINEQs RHS in interval arithmetic
  static int MC_CVRHSI__
    ( realtype t, N_Vector Nx, N_Vector Nxdot, void *user_data );

  //! @brief Static wrapper to function computing the quadratures RHS in interval arithmetic
  static int MC_CVQUADI__
    ( realtype t, N_Vector Nx, N_Vector Nqdot, void *user_data );

  //! @brief Function to initialize GSL for state polynomial models
  bool _INI_PM_STA
    ( const unsigned np, const PVT*PMp, const unsigned ns );

  //! @brief Static wrapper to function computing the DINEQs RHS in polynomial model arithmetic
  static int MC_CVRHSPM__
    ( realtype t, N_Vector Nx, N_Vector Nxdot, void *user_data );

  //! @brief Static wrapper to function computing the quadratures RHS in polynomial model arithmetic
  static int MC_CVQUADPM__
    ( realtype t, N_Vector Nx, N_Vector Nqdot, void *user_data );

  //! @brief Private methods to block default compiler methods
  ODEBND_SUNDIALS(const ODEBND_SUNDIALS&);
  ODEBND_SUNDIALS& operator=(const ODEBND_SUNDIALS&);
};

template <typename T, typename PMT, typename PVT>
 ODEBND_SUNDIALS<T,PMT,PVT>* ODEBND_SUNDIALS<T,PMT,PVT>::_pODEBND = 0;

template <typename T, typename PMT, typename PVT> inline
ODEBND_SUNDIALS<T,PMT,PVT>::ODEBND_SUNDIALS
()
: BASE_SUNDIALS(), ODEBND_BASE<T,PMT,PVT>(), _cv_mem(0), _cv_flag(0),
  _Nx(0), _Nq(0)
{}

template <typename T, typename PMT, typename PVT> inline
ODEBND_SUNDIALS<T,PMT,PVT>::~ODEBND_SUNDIALS
()
{
  if( _Nx ) N_VDestroy_Serial( _Nx );
  if( _Nq ) N_VDestroy_Serial( _Nq );
  if( _cv_mem ) CVodeFree( &_cv_mem );
}

template <typename T, typename PMT, typename PVT> inline bool
ODEBND_SUNDIALS<T,PMT,PVT>::_INI_CVODE
( CVRhsFn MC_CVRHS, CVRhsFn MC_CVQUAD )
{
  // Reset CVode memory block
  if( _cv_mem ) CVodeFree( &_cv_mem );
  switch( options.INTMETH ){
   // Non-stiff case: create CVode memory block for the ADAMS method
   // and use of a functional iteration
   case Options::MSADAMS: default:
    _cv_mem = CVodeCreate( CV_ADAMS, CV_FUNCTIONAL );
    if( _check_cv_flag((void *)_cv_mem, "CVodeCreate", 0) ) return false;
    break;

   // Stiff case: create CVode memory block for the BDF method
   // and use of a Newton iteration
   case Options::MSBDF:
    _cv_mem = CVodeCreate( CV_BDF, CV_NEWTON );
    if( _check_cv_flag((void *)_cv_mem, "CVodeCreate", 0) ) return false;
    break;
  }

  // Initialize CVode memory and specify right hand side function,
  // current time _t, and current state _Nx
  _cv_flag = CVodeInit( _cv_mem, MC_CVRHS, _t, _Nx );
  if( _check_cv_flag(&_cv_flag, "CVodeInit", 1) ) return false;

   // Specify the Jacobian approximation and linear solver
  if( options.INTMETH == Options::MSBDF ){
    switch( options.JACAPPROX ){
     case Options::CV_DIAG: default:
       _cv_flag = CVDiag( _cv_mem );
       if( _check_cv_flag(&_cv_flag, "CVDiag", 1)) return false;
       break;

     case Options::CV_LAPACKDENSE:
       //_cv_flag = CVLapackDense( _cv_mem, NV_LENGTH_S( _Nx ) );
       //if( _check_cv_flag(&_cv_flag, "CVLapackDense", 1)) return false;
       //break;

     case Options::CV_DENSE:
       _cv_flag = CVDense( _cv_mem, NV_LENGTH_S( _Nx ) );
       if( _check_cv_flag(&_cv_flag, "CVDense", 1)) return false;
       break;
    }
    _cv_flag = CVDlsSetDenseJacFn( _cv_mem, NULL );
    if ( _check_cv_flag(&_cv_flag, "CVDlsSetDenseJacFn", 1) ) return false;
  }

  // Specify the relative and absolute tolerances for states
  _cv_flag = CVodeSStolerances( _cv_mem, options.RTOL, options.ATOL );
  if( _check_cv_flag(&_cv_flag, "CVodeSStolerances", 1) ) return false;

  // Set maximum number of error test failures
  _cv_flag = CVodeSetMaxErrTestFails( _cv_mem, options.MAXFAIL );
  if ( _check_cv_flag(&_cv_flag, "CVodeSetMaxErrTestFails", 1) ) return false;

  // Specify minimum stepsize
  _cv_flag = CVodeSetMinStep( _cv_mem, options.HMIN>0.? options.HMIN:0. );
  if( _check_cv_flag(&_cv_flag, "CVodeSetMinStep", 1) ) return false;

  // Specify maximum stepsize
  _cv_flag = CVodeSetMaxStep( _cv_mem, options.HMAX>0.? options.HMAX: 0. );
  if( _check_cv_flag(&_cv_flag, "CVodeSetMaxStep", 1) ) return false;

  // Specify maximum number of steps between two stage times
  _cv_flag = CVodeSetMaxNumSteps( _cv_mem, options.NMAX );
  if( _check_cv_flag(&_cv_flag, "CVodeSetMaxNumSteps", 1) ) return false;

  // Initialize the integrator memory for the quadrature variables
  if( !_Nq ) return true;
  _cv_flag = CVodeQuadInit( _cv_mem, MC_CVQUAD, _Nq );
  if( _check_cv_flag(&_cv_flag, "CVodeQuadInit", 1) ) return false;

  // Specify whether or not to perform error control on quadrature
  _cv_flag = CVodeSetQuadErrCon( _cv_mem, options.QERR );
  if( _check_cv_flag(&_cv_flag, "CVodeSetQuadErrCon", 1) ) return false;

  // Specify the relative and absolute tolerances for quadratures
  _cv_flag = CVodeQuadSStolerances( _cv_mem, options.RTOL, options.ATOL );
  if( _check_cv_flag(&_cv_flag, "CVodeQuadSStolerances", 1) ) return false;

  return true;
}

template <typename T, typename PMT, typename PVT> inline bool
ODEBND_SUNDIALS<T,PMT,PVT>::_CC_CVODE
()
{
  // Reinitialize CVode memory block for current time _t and
  // current state _Nx
  _cv_flag = CVodeReInit( _cv_mem, _t, _Nx );
  if( _check_cv_flag(&_cv_flag, "CVodeReInit", 1) ) return false;

  // Reinitialize CVode memory block for current quarature _Nq
  if( !_Nq ) return true;
  _cv_flag = CVodeQuadReInit( _cv_mem, _Nq );
  if( _check_cv_flag(&_cv_flag, "CVodeQuadReInit", 1) ) return false;

  return true;
}

template <typename T, typename PMT, typename PVT> inline void
ODEBND_SUNDIALS<T,PMT,PVT>::_END_STA()
{
  // Get final CPU time
  _final_stats( stats_sta );
}

template <typename T, typename PMT, typename PVT> inline bool
ODEBND_SUNDIALS<T,PMT,PVT>::_INI_I_STA
( const unsigned np, const T*Ip, const unsigned ns )
{
  // Initialize bound propagation
  if( !ODEBND_BASE<T,PMT,PVT>::_INI_I_STA( options, np, Ip, ns ) )
    return false;

  // Set SUNDIALS state/quadrature arrays
  unsigned cv_Nx_size = 0, cv_Nq_size = 2*_nq;
  switch( options.WRAPMIT){
  case Options::NONE:
  case Options::DINEQ:
    cv_Nx_size += 2*_nx;
    break;
  case Options::ELLIPS:
  default:
    cv_Nx_size += _nx*(1+np)+_nx*(_nx+1)/2;
    break;
  }
  if( !_Nx || cv_Nx_size != NV_LENGTH_S( _Nx ) ){
    if( _Nx ) N_VDestroy_Serial( _Nx );
    _Nx  = N_VNew_Serial( cv_Nx_size );
  }
  if( !_Nq || cv_Nq_size != NV_LENGTH_S( _Nq ) ){
    if( _Nq ) N_VDestroy_Serial( _Nq );
    _Nq  = cv_Nq_size? N_VNew_Serial( cv_Nq_size ): 0;
  }

  // Reset result record and statistics
  results_sta.clear();
  _init_stats( stats_sta );

  return true;
}

template <typename T, typename PMT, typename PVT> inline int
ODEBND_SUNDIALS<T,PMT,PVT>::MC_CVRHSI__
( realtype t, N_Vector y, N_Vector ydot, void *user_data )
{
  ODEBND_SUNDIALS<T,PMT,PVT> *pODEBND = ODEBND_SUNDIALS<T,PMT,PVT>::_pODEBND;
  bool flag = pODEBND->_RHS_I_STA( pODEBND->options, t, NV_DATA_S( y ),
    NV_DATA_S( ydot ) );
  ODEBND_SUNDIALS<T,PMT,PVT>::_pODEBND = pODEBND;
  pODEBND->stats_sta.numRHS++;
  return( flag? 0: -1 );
}

template <typename T, typename PMT, typename PVT> inline int
ODEBND_SUNDIALS<T,PMT,PVT>::MC_CVQUADI__
( realtype t, N_Vector y, N_Vector qdot, void *user_data )
{
  ODEBND_SUNDIALS<T,PMT,PVT> *pODEBND = ODEBND_SUNDIALS<T,PMT,PVT>::_pODEBND;
  bool flag = pODEBND->_RHS_I_QUAD( pODEBND->options, t, NV_DATA_S( y ),
    NV_DATA_S( qdot ) );
  ODEBND_SUNDIALS<T,PMT,PVT>::_pODEBND = pODEBND;
  //pODEBND->stats_sta.numRHS++;
  return( flag? 0: -1 );
}

//! @fn template <typename T, typename PMT, typename PVT> inline typename ODEBND_SUNDIALS<T,PMT,PVT>::STATUS ODEBND_SUNDIALS<T,PMT,PVT>::bounds(
//! const unsigned ns, const double*tk, const T*Ip, T**Ixk=0, T*Iq=0, T*If=0,
//! std::ostream&os=std::cout )
//!
//! This function computes an interval enclosure of the reachable set of 
//! the parametric ODEs defined in IVP using equally spaced samples:
//!   - <a>ns</a> [input] number of time stages
//!   - <a>tk</a> [input] stage times, including the initial time
//!   - <a>Ip</a> [input] interval parameter set
//!   - <a>Ixk</a> [output] interval state enclosures at stage times (default: NULL)
//!   - <a>Iq</a> [output] interval quadrature enclosures at final time (default: NULL)
//!   - <a>If</a> [output] interval function enclosures (default: NULL)
//!   - <a>os</a> [input] output stream (default: std::cout)
//! .
//! The return value is the status.
template <typename T, typename PMT, typename PVT>
inline typename ODEBND_SUNDIALS<T,PMT,PVT>::STATUS
ODEBND_SUNDIALS<T,PMT,PVT>::bounds
( const unsigned ns, const double*tk, const T*Ip, T**Ixk, T*Iq, T*If,
  std::ostream&os )
{
  // Check size
  if( !tk || !Ixk || !Ip || (_nf && !If) ) return FATAL;

  try{
    // Initialize trajectory integration
    if( !_INI_I_STA( _np, Ip, ns ) ) return FATAL;

    // Bounds on initial states/quadratures
    _t = tk[0];
    if( !_IC_I_STA( options, NV_DATA_S( _Nx ) )
     || (_Nq && !_IC_I_QUAD( options, NV_DATA_S( _Nq )) ) )
      { _END_STA(); return FATAL; }
    if( options.DISPLAY >= 1 ){
      _print_interm( _t, _nx, _Ix, "x", os );
      _print_interm( _nq, _Iq, "q", os );
    }
    if( Ixk && !Ixk[0] ) Ixk[0] = new T[_nx];
    for( unsigned ix=0; Ixk[0] && ix<_nx; ix++ ) Ixk[0][ix] = _Ix[ix];

    // Record initial results
    if( options.RESRECORD )
      results_sta.push_back( Results( tk[0], _nx, Ixk[0] ) );

    // Integrate ODEs through each stage using SUNDIALS
    _pODEBND = this;
    if( !_INI_CVODE( MC_CVRHSI__, MC_CVQUADI__ ) )
      { _END_STA(); return FATAL; }

    for( _istg=0; _istg<ns; _istg++ ){
      // Bounds on state discontinuities (if any) at stage times
      // and integrator reinitialization (if applicable)
      _pos_ic = ( _vIC.size()>=ns? _istg:0 );
      if( _pos_ic
       && ( !_CC_I_STA( options, _pos_ic, _t, NV_DATA_S( _Nx ) )
         || !_CC_CVODE() ) )
        { _END_STA(); return FAILURE; }

      // update list of operations in RHS, JAC and QUAD
      _pos_rhs  = ( _vRHS.size()<=1? 0: _istg );
      _pos_quad = ( _vQUAD.size()<=1? 0: _istg );
      if( (!_istg || _pos_rhs || _pos_quad)
        && !_SET_I_STA( options, _pos_rhs, _pos_quad ) )
        { _END_STA(); return FATAL; }

      // integrate till end of time stage
      _cv_flag = CVodeSetStopTime( _cv_mem, tk[_istg+1] );
      if( _check_cv_flag(&_cv_flag, "CVodeSetStopTime", 1) )
        { _END_STA(); return FATAL; }
      while( _t < tk[_istg+1] ){
        _cv_flag = CVode( _cv_mem, tk[_istg+1], _Nx, &_t, CV_ONE_STEP );
        if( _check_cv_flag(&_cv_flag, "CVode", 1)
         || (options.NMAX && stats_sta.numSteps > options.NMAX)
         || _diam( _nx, _Ix ) > options.DMAX )
          throw Exceptions( Exceptions::INTERN );
        stats_sta.numSteps++;
      }

      // Bounds on intermediate states
      switch( options.WRAPMIT){
      case Options::NONE:
      case Options::DINEQ:
        _vec2I( NV_DATA_S( _Nx ), _nx, _Ix );
        break;
      case Options::ELLIPS:
      default:
        _vec2E( NV_DATA_S( _Nx ), _nx, _np, _Q, _Er, _Ir, _pref, _Ip, _B, _xref, _Ix );
        break;
      }

      // Bounds on intermediate quadratures
      if( _nq ){
        _cv_flag = CVodeGetQuad( _cv_mem, &_t, _Nq );
        if( _check_cv_flag(&_cv_flag, "CVodeGetQuad", 1) )
          { _END_STA(); return FATAL; }
        _vec2I( NV_DATA_S( _Nq ), _nq, _Iq );
      }

      // Keep track/display/record stage results
      if( options.DISPLAY >= 1 ){
        _print_interm( _t, _nx, _Ix, "x", os );
        _print_interm( _nq, _Iq, "q", os );
      }
      if( Ixk && !Ixk[_istg+1] ) Ixk[_istg+1] = new T[_nx];
      for( unsigned ix=0; Ixk[_istg+1] && ix<_nx; ix++ ) Ixk[_istg+1][ix] = _Ix[ix];
      if( options.RESRECORD )
        results_sta.push_back( Results( tk[_istg+1], _nx, Ixk[_istg+1] ) );

      // Add intermediate function terms
      _pos_fct = ( _vFCT.size()>=ns? _istg:0 );
      if( (_vFCT.size()>=ns || _istg==ns-1) && !_FCT_I_STA( _pos_fct, _t, If ) )
        { _END_STA(); return FATAL; }
    }

    // Bounds on final quadratures and functions
    for( unsigned iq=0; Iq && iq<_nq; iq++ ) Iq[iq] = _Iq[iq];
    if( options.DISPLAY >= 1 ) _print_interm( _nf, If, "f", os );
  }
  catch(...){
    _END_STA();
    if( options.DISPLAY >= 1 ) _print_stats( stats_sta, os );
    return FAILURE;
  }

  _END_STA();
  if( options.DISPLAY >= 1 ) _print_stats( stats_sta, os );
  return NORMAL;
}

//! @fn template <typename T, typename PMT, typename PVT> template<typename ODESLV> inline typename ODEBND_SUNDIALS<T,PMT,PVT>::STATUS ODEBND_SUNDIALS<T,PMT,PVT>::hausdorff(
//! const unsigned ns, const double*tk, const T*Ip, double**Hxk,
//! ODESLV*traj, const unsigned nsamp, std::ostream&os )
//!
//! This function computes the Hausdorff distance between the interval enclosure
//! and the exact reachable set projected onto each variable:
//!   - <a>ns</a> [input] number of time stages
//!   - <a>tk</a> [input] stage times, including the initial time
//!   - <a>Ip</a> [input] interval parameter set
//!   - <a>Hxk</a> [output] Hausdorff distance between the interval enclosure
//!     and the exact reachable set projected onto each variable, at stage times
//!   - <a>traj</a> [input] trajectory ODE integrator
//!   - <a>nsamp</a> [input] number of samples for each parameter
//!   - <a>os</a> [input] output stream
//! .
//! The return value is the status.
template <typename T, typename PMT, typename PVT> template<typename ODESLV> inline typename ODEBND_SUNDIALS<T,PMT,PVT>::STATUS
ODEBND_SUNDIALS<T,PMT,PVT>::hausdorff
( const unsigned ns, const double*tk, const T*Ip, double**Hxk,
  ODESLV&traj, const unsigned nsamp, std::ostream&os )
{
  return _hausdorff( ns, tk, Ip, Hxk, *this, traj, nsamp, os )? NORMAL: FAILURE;
}

template <typename T, typename PMT, typename PVT> inline int
ODEBND_SUNDIALS<T,PMT,PVT>::MC_CVRHSPM__
( realtype t, N_Vector y, N_Vector ydot, void *user_data )
{
  ODEBND_SUNDIALS<T,PMT,PVT> *pODEBND = ODEBND_SUNDIALS<T,PMT,PVT>::_pODEBND;
  bool flag = pODEBND->_RHS_PM_STA( pODEBND->options, t, NV_DATA_S( y ),
    NV_DATA_S( ydot ) );
  ODEBND_SUNDIALS<T,PMT,PVT>::_pODEBND = pODEBND;
  pODEBND->stats_sta.numRHS++;
  return( flag? 0: -1 );
}

template <typename T, typename PMT, typename PVT> inline int
ODEBND_SUNDIALS<T,PMT,PVT>::MC_CVQUADPM__
( realtype t, N_Vector y, N_Vector qdot, void *user_data )
{
  ODEBND_SUNDIALS<T,PMT,PVT> *pODEBND = ODEBND_SUNDIALS<T,PMT,PVT>::_pODEBND;
  bool flag = pODEBND->_RHS_PM_QUAD( pODEBND->options, t, NV_DATA_S( y ),
    NV_DATA_S( qdot ) );
  ODEBND_SUNDIALS<T,PMT,PVT>::_pODEBND = pODEBND;
  //pODEBND->stats_sta.numRHS++;
  return( flag? 0: -1 );
}

template <typename T, typename PMT, typename PVT> inline bool
ODEBND_SUNDIALS<T,PMT,PVT>::_INI_PM_STA
( const unsigned np, const PVT* PMp, const unsigned ns )
{
  // Initialize bound propagation
  if( !ODEBND_BASE<T,PMT,PVT>::_INI_PM_STA( options, np, PMp, ns ) )
    return false;

  // Set SUNDIALS state/quadrature arrays
  unsigned cv_Nx_size = _PMenv->nmon()*_nx, cv_Nq_size = (_PMenv->nmon()+1)*_nq;
  switch( options.WRAPMIT){
  case Options::NONE:
    cv_Nx_size += _nx;
    break;
  case Options::DINEQ:
    cv_Nx_size += 2*_nx;
    break;
  case Options::ELLIPS:
  default:
    cv_Nx_size += _nx*(_nx+1)/2;
    break;
  }
  if( !_Nx || cv_Nx_size != NV_LENGTH_S( _Nx ) ){
    if( _Nx ) N_VDestroy_Serial( _Nx );
    _Nx  = N_VNew_Serial( cv_Nx_size );
  }
  if( !_Nq || cv_Nq_size != NV_LENGTH_S( _Nq ) ){
    if( _Nq ) N_VDestroy_Serial( _Nq );
    _Nq  = cv_Nq_size? N_VNew_Serial( cv_Nq_size ): 0;
  }

  // Reset result record and statistics
  results_sta.clear();
  _init_stats( stats_sta );

  return true;
}

//! @fn template <typename T, typename PMT, typename PVT> inline typename ODEBND_SUNDIALS<T,PMT,PVT>::STATUS ODEBND_SUNDIALS<T,PMT,PVT>::bounds(
//! const unsigned ns, const double*tk, const PVT*PMp, PVT**PMxk=0, 
//! PVT*PMq=0, PVT*PMf=0, std::ostream&os=std::cout )
//!
//! This function computes an enclosure of the reachable set of the parametric ODEs
//! using propagation of polynomial models with convex remainders (intervals, ellipsoids):
//!   - <a>ns</a> [input] number of time stages
//!   - <a>tk</a> [input] stage times, including the initial time
//!   - <a>PMp</a> [input] polynomial model of parameter set
//!   - <a>PMxk</a> [output] polynomial model of state enclosures at stage times (default: NULL)
//!   - <a>PMq</a> [output] polynomial model of quadrature variables (default: NULL)
//!   - <a>PMf</a> [output] polynomial model of state/quadrature functionals (default: NULL)
//!   - <a>os</a> [input] output stream (default: std::cout)
//! .
//! The return value is the status.
template <typename T, typename PMT, typename PVT>
inline typename ODEBND_SUNDIALS<T,PMT,PVT>::STATUS
ODEBND_SUNDIALS<T,PMT,PVT>::bounds
( const unsigned ns, const double*tk, const PVT*PMp, PVT**PMxk,
  PVT*PMq, PVT*PMf, std::ostream&os )
{
  // Check arguments
  if( !tk || !PMxk || !PMp || (_nf && !PMf) ) return FATAL;

  try{
    // Initialize trajectory integration
    if( !_INI_PM_STA( _np, PMp, ns ) ) return FATAL;

    // Bounds on initial states/quadratures
    _t = tk[0];
    if( !_IC_PM_STA( options, NV_DATA_S( _Nx ) )
     || (_Nq && !_IC_PM_QUAD( options, NV_DATA_S( _Nq )) ) )
      { _END_STA(); return FATAL; }
    if( options.DISPLAY >= 1 ){
      _print_interm( _t, _nx, _PMx, "x", os );
      _print_interm( _nq, _PMq, "q", os );
    }
    if( PMxk && !PMxk[0] ) PMxk[0] = new PVT[_nx];
    for( unsigned ix=0; PMxk[0] && ix<_nx; ix++ ) PMxk[0][ix] = _PMx[ix];

    // Record initial results
    if( options.RESRECORD )
      results_sta.push_back( Results( tk[0], _nx, PMxk[0] ) );

    // Integrate ODEs through each stage using SUNDIALS
    _pODEBND = this;
    if( !_INI_CVODE( MC_CVRHSPM__, MC_CVQUADPM__ ) )
      { _END_STA(); return FATAL; }

    for( _istg=0; _istg<ns; _istg++ ){
      // Bounds on state discontinuities (if any) at stage times
      // and integrator reinitialization (if applicable)
      _pos_ic = ( _vIC.size()>=ns? _istg:0 );
      if( _pos_ic
       && ( !_CC_PM_STA( options, _pos_ic, _t, NV_DATA_S( _Nx ) )
         || !_CC_CVODE() ) )
        { _END_STA(); return FAILURE; }

      // update list of operations in RHS, JAC and QUAD
      _pos_rhs  = ( _vRHS.size()<=1? 0: _istg );
      _pos_quad = ( _vQUAD.size()<=1? 0: _istg );
      if( (!_istg || _pos_rhs || _pos_quad)
        && !_SET_PM_STA( options, _pos_rhs, _pos_quad ) )
        { _END_STA(); return FATAL; }

      // integrate till end of time stage
      _cv_flag = CVodeSetStopTime( _cv_mem, tk[_istg+1] );
      if( _check_cv_flag(&_cv_flag, "CVodeSetStopTime", 1) )
        { _END_STA(); return FATAL; }
      while( _t < tk[_istg+1] ){
        _cv_flag = CVode( _cv_mem, tk[_istg+1], _Nx, &_t, CV_ONE_STEP );
        if( _check_cv_flag(&_cv_flag, "CVode", 1)
         || (options.NMAX && stats_sta.numSteps > options.NMAX)
         || _diam( _nx, _PMx ) > options.DMAX
         || _diam( _nx, _Ir  ) > options.DMAX )
          throw Exceptions( Exceptions::INTERN );
        stats_sta.numSteps++;
      }

      // Bounds on intermediate states
      switch( options.WRAPMIT){
      case Options::NONE:
      case Options::DINEQ:
        _vec2PMI( NV_DATA_S( _Nx ), _PMenv, _nx, _PMx, false );
        break;
      case Options::ELLIPS:
      default:
        _vec2PME( NV_DATA_S( _Nx ), _PMenv, _nx, _PMx, _Q, _Er, _Ir );
        //std::cout << _Er.eigQ().first;
        break;
      }

      // Bounds on intermediate quadratures
      if( _nq ){
        _cv_flag = CVodeGetQuad( _cv_mem, &_t, _Nq );
        if( _check_cv_flag(&_cv_flag, "CVodeGetQuad", 1) )
          { _END_STA(); return FATAL; }
        _vec2PMI( NV_DATA_S( _Nq ), _PMenv, _nq, _PMq, true );
      }

      // Keep track/display/record stage results
      if( options.DISPLAY >= 1 ){
        _print_interm( _t, _nx, _PMx, "x", os );
        _print_interm( _nq, _PMq, "q", os );
      }
      if( PMxk && !PMxk[_istg+1] ) PMxk[_istg+1] = new PVT[_nx];
      for( unsigned ix=0; PMxk[_istg+1] && ix<_nx; ix++ ) PMxk[_istg+1][ix] = _PMx[ix];
      if( options.RESRECORD )
        results_sta.push_back( Results( tk[_istg+1], _nx, PMxk[_istg+1] ) );

      // Add intermediate function terms
      _pos_fct = ( _vFCT.size()>=ns? _istg:0 );
      if( (_vFCT.size()>=ns || _istg==ns-1) && !_FCT_PM_STA( _pos_fct, _t, PMf ) )
        { _END_STA(); return FATAL; }
    }

    // Bounds on final quadratures and functions
    for( unsigned iq=0; PMq && iq<_nq; iq++ ) PMq[iq] = _PMq[iq];
    if( options.DISPLAY >= 1 ) _print_interm( _nf, PMf, "f", os );
  }
  catch(...){
    _END_STA();
    if( options.DISPLAY >= 1 ) _print_stats( stats_sta, os );
    return FAILURE;
  }

  _END_STA();
  if( options.DISPLAY >= 1 ) _print_stats( stats_sta, os );
  return NORMAL;
}

//! @fn template <typename T, typename PMT, typename PVT> template<typename ODESLV> inline typename ODEBND_SUNDIALS<T,PMT,PVT>::STATUS ODEBND_SUNDIALS<T,PMT,PVT>::hausdorff(
//! const unsigned ns, const double*tk, const PVT*PMp, double**Hxk,
//! ODESLV&traj, const unsigned nsamp, std::ostream&os )
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
//!   - <a>traj</a> [input] trajectory ODE integrator
//!   - <a>nsamp</a> [input] number of samples for each parameter
//!   - <a>os</a> [input] output stream (default: std::cout)

template <typename T, typename PMT, typename PVT> template<typename ODESLV> inline typename ODEBND_SUNDIALS<T,PMT,PVT>::STATUS
ODEBND_SUNDIALS<T,PMT,PVT>::hausdorff
( const unsigned ns, const double*tk, const PVT*PMp, double**Hxk,
  ODESLV&traj, const unsigned nsamp, std::ostream&os )
{
  return _hausdorff( ns, tk, PMp, Hxk, *this, traj, nsamp, os )? NORMAL: FAILURE;
}

} // end namescape mc

#endif

