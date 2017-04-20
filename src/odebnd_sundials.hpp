// Copyright (C) 2015 Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

#ifndef MC__ODEBND_SUNDIALS_HPP
#define MC__ODEBND_SUNDIALS_HPP

#undef  MC__ODEBND_SUNDIALS_DINEQI_DEBUG
#undef  MC__ODEBND_SUNDIALS_DINEQPM_DEBUG

#include "odebnd_base.hpp"
#include "base_sundials.hpp"
#include "odeslv_sundials.hpp"

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
class ODEBND_SUNDIALS:
  public virtual BASE_DE,
  public virtual BASE_SUNDIALS,
  public virtual ODEBND_BASE<T,PMT,PVT>
{
  typedef Ellipsoid E;
  typedef BASE_DE::STATUS STATUS;
  typedef int (*CVRhsFn)( realtype t, N_Vector y, N_Vector ydot, void *user_data );

 protected:
  using ODEBND_BASE<T,PMT,PVT>::_diam;
  using ODEBND_BASE<T,PMT,PVT>::_hausdorff;
  using ODEBND_BASE<T,PMT,PVT>::_bounds;
  using ODEBND_BASE<T,PMT,PVT>::_print_interm;
  using ODEBND_BASE<T,PMT,PVT>::_PMenv;

  using ODEBND_BASE<T,PMT,PVT>::_Ix;
  using ODEBND_BASE<T,PMT,PVT>::_Ir;
  using ODEBND_BASE<T,PMT,PVT>::_Er;
  using ODEBND_BASE<T,PMT,PVT>::_PMx;
  using ODEBND_BASE<T,PMT,PVT>::_GET_I_STA;
  using ODEBND_BASE<T,PMT,PVT>::_IC_I_SET;
  using ODEBND_BASE<T,PMT,PVT>::_IC_I_STA;
  using ODEBND_BASE<T,PMT,PVT>::_IC_I_QUAD;
  using ODEBND_BASE<T,PMT,PVT>::_CC_I_SET;
  using ODEBND_BASE<T,PMT,PVT>::_CC_I_STA;
  using ODEBND_BASE<T,PMT,PVT>::_RHS_I_SET;
  using ODEBND_BASE<T,PMT,PVT>::_RHS_I_STA;
  using ODEBND_BASE<T,PMT,PVT>::_RHS_I_QUAD;
  using ODEBND_BASE<T,PMT,PVT>::_JAC_I_STA;
  using ODEBND_BASE<T,PMT,PVT>::_FCT_I_STA;
  using ODEBND_BASE<T,PMT,PVT>::_GET_PM_STA;
  using ODEBND_BASE<T,PMT,PVT>::_IC_PM_SET;
  using ODEBND_BASE<T,PMT,PVT>::_IC_PM_STA;
  using ODEBND_BASE<T,PMT,PVT>::_IC_PM_QUAD;
  using ODEBND_BASE<T,PMT,PVT>::_CC_PM_SET;
  using ODEBND_BASE<T,PMT,PVT>::_CC_PM_STA;
  using ODEBND_BASE<T,PMT,PVT>::_RHS_PM_SET;
  using ODEBND_BASE<T,PMT,PVT>::_RHS_PM_STA;
  using ODEBND_BASE<T,PMT,PVT>::_RHS_PM_QUAD;
  using ODEBND_BASE<T,PMT,PVT>::_JAC_PM_STA;
  using ODEBND_BASE<T,PMT,PVT>::_FCT_PM_STA;

 protected:
  //! @brief Pointer to the CVODE memory block
  void *_cv_mem;

  //! @brief Return flag for SUNDIALS methods
  int _cv_flag;

  //! @brief N_Vector object holding current state parameterization
  N_Vector _Nx;

  //! @brief N_Vector object holding current state quadrature parameterization
  N_Vector _Nq;

  //! @brief N_Vector object holding current state function parameterization
  N_Vector _Nf;
  
  //! @brief N_Vector object holding absolute tolerances
  N_Vector _NTOLx;

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

  //! @brief checkpoints for adjoint integration
  int _nchk;

  //! @brief state parameterizations at time stages for adjoint integration
  std::vector< std::vector<realtype> > _vec_sta;

  //! @brief static pointer to class
  static ODEBND_SUNDIALS<T,PMT,PVT> *pODEBND;

  //! @brief static pointer to local integrator class
  ODESLV_SUNDIALS pODESLV;

public:
  /** @defgroup ODEBND Continuous-time set-valued integration of parametric ODEs
   *  @{
   */
  typedef typename ODEBND_BASE<T,PMT,PVT>::Results Results;

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
      BASE_SUNDIALS::Options(), WRAPMIT(ELLIPS), ORDMIT(-2), PMOPT(typename PMT::Options()),
      QSCALE(1e-5), PMNOREM(false), USEINV(true), DMAX(1e20), DISPLAY(1), 
      RESRECORD(false), ODESLV(typename ODESLV_SUNDIALS::Options())
      {}
    //! @brief Assignment operator
    template <typename OPT> Options& operator=
      ( OPT&options ){
        BASE_SUNDIALS::Options::operator=(options);
        WRAPMIT   = options.WRAPMIT;
        ORDMIT    = options.ORDMIT;
        PMOPT     = options.PMOPT;
        QSCALE    = options.QSCALE;
        PMNOREM   = options.PMNOREM;
        USEINV    = options.USEINV;
        DMAX      = options.DMAX;
        DISPLAY   = options.DISPLAY;
        RESRECORD = options.RESRECORD;
        ODESLV    = options.ODESLV;
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
    //! @brief Option for wrapping mitigation strategy with ellipsoidal bounder (Default: -2)
    int ORDMIT;
    //! @brief Options of internal polynomial model arithmetic for wrapping mitigation
    typename PMT::Options PMOPT;
    //! @brief Tolerance for applying ellipsoidal remainder scaling (Default: machprec())
    double QSCALE;
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
    //! @brief Options of real-valued integrator for bounds sampling
    typename ODESLV_SUNDIALS::Options ODESLV;
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

  //! @brief Vector storing interval bound results (upon request only)
  std::vector< Results > results_sta;

  //! @brief Statistics for state bounds integration
  Stats stats_sta;

  //! @brief Computes interval enclosure of reachable set of parametric ODEs
  STATUS bounds
    ( const T*Ip, T**Ixk, T*If, std::ostream&os=std::cout );

  //! @brief Computes Hausdorff distance between interval enclosure and actual reachable set of parametric ODEs, using parameter sampling
  STATUS hausdorff
    ( const T*Ip, double**Hxk, double*Hf, const unsigned nsamp, std::ostream&os=std::cout );

  //! @brief Computes polynomial model enclosure of reachable set of parametric ODEs
  STATUS bounds
    ( const PVT*PMp, PVT**PMxk, PVT*PMf, std::ostream&os=std::cout );

  //! @brief Computes polynomial model enclosure of reachable set of parametric ODEs
  STATUS bounds
    ( const PVT*PMp, PVT**PMxk, E*ERxk, PVT*PMf, std::ostream&os=std::cout );

  //! @brief Computes Hausdorff distance between polynomial model remainder enclosure and actual remainder function range, using parameter sampling
  STATUS hausdorff
    ( const PVT*PMp, double**Hxk, double*Hf, const unsigned nsamp, std::ostream&os=std::cout );

 //! @brief Compute approximate interval enclosure of reachable set of parametric ODEs using parameter sampling
  STATUS bounds
    ( const T*Ip, T**Ixk, T*If, const unsigned nsamp, std::ostream&os=std::cout );

  //! @brief Record results in file <a>bndrec</a>, with accuracy of <a>iprec</a> digits
  void record
    ( std::ofstream&bndrec, const unsigned iprec=5 ) const
    { return ODEBND_BASE<T,PMT,PVT>::_record( bndrec, results_sta, iprec ); }
  /** @} */

protected:

  //! @brief Function to initialize CVode memory block
  virtual bool _INI_CVODE
    ( CVRhsFn MC_CVRHS, CVRhsFn MC_CVQUAD );

  //! @brief Function to reinitialize CVode memory block
  bool _CC_CVODE_STA
    ();

  //! @brief Function to reinitialize CVode memory block
  bool _CC_CVODE_QUAD
    ();

  //! @brief Function to finalize state bounding
  void _END_STA
    ();

  //! @brief Function to initialize state interval bounding
  bool _INI_I_STA
    ( const unsigned np, const T*Ip );

  //! @brief Function to initialize state interval bounding
  bool _INI_I_STA
    ( const unsigned np );

  //! @brief Static wrapper to function computing the DINEQs RHS in interval arithmetic
  static int MC_CVRHSI__
    ( realtype t, N_Vector Nx, N_Vector Nxdot, void *user_data );

  //! @brief Static wrapper to function computing the quadratures RHS in interval arithmetic
  static int MC_CVQUADI__
    ( realtype t, N_Vector Nx, N_Vector Nqdot, void *user_data );

  //! @brief Propagate state/quadrature interval bounds forward in time through every time stages
  virtual STATUS _bounds
    ( const T*Ip, T**Ixk, T*If, const bool store, std::ostream&os );

  //! @brief Function to initialize GSL for state polynomial models
  bool _INI_PM_STA
    ( const unsigned np, const PVT*PMp );

  //! @brief Function to initialize GSL for state polynomial models
  bool _INI_PM_STA
    ( const unsigned np );

  //! @brief Static wrapper to function computing the DINEQs RHS in polynomial model arithmetic
  static int MC_CVRHSPM__
    ( realtype t, N_Vector Nx, N_Vector Nxdot, void *user_data );

  //! @brief Static wrapper to function computing the quadratures RHS in polynomial model arithmetic
  static int MC_CVQUADPM__
    ( realtype t, N_Vector Nx, N_Vector Nqdot, void *user_data );

  //! @brief Propagate state/quadrature polynomial models forward in time through every time stages
  virtual STATUS _bounds
    ( const PVT*PMp, PVT**PMxk, E*ERxk, PVT*PMf, const bool store, std::ostream&os );
 
  //! @brief Private methods to block default compiler methods
  ODEBND_SUNDIALS(const ODEBND_SUNDIALS&);
  ODEBND_SUNDIALS& operator=(const ODEBND_SUNDIALS&);
};

template <typename T, typename PMT, typename PVT>
 ODEBND_SUNDIALS<T,PMT,PVT>* ODEBND_SUNDIALS<T,PMT,PVT>::pODEBND = 0;

template <typename T, typename PMT, typename PVT> inline
ODEBND_SUNDIALS<T,PMT,PVT>::ODEBND_SUNDIALS
()
: _cv_mem(0), _cv_flag(0), _Nx(0), _Nq(0), _Nf(0), _NTOLx(0)
{}

template <typename T, typename PMT, typename PVT> inline
ODEBND_SUNDIALS<T,PMT,PVT>::~ODEBND_SUNDIALS
()
{
  if( _Nx )    N_VDestroy_Serial( _Nx );
  if( _Nq )    N_VDestroy_Serial( _Nq );
  if( _Nf )    N_VDestroy_Serial( _Nf );
  if( _NTOLx ) N_VDestroy_Serial( _NTOLx );
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
  if( options.INTMETH == Options::MSBDF ){ // Applies to Newton iteration only
    switch( options.JACAPPROX ){
     case Options::CV_DIAG: default:
       _cv_flag = CVDiag( _cv_mem );
       if( _check_cv_flag(&_cv_flag, "CVDiag", 1)) return false;
       break;

     case Options::CV_LAPACKDENSE:
       //_cv_flag = CVLapackDense( _cv_mem, NV_LENGTH_S( _Nx ) );
       //if( _check_cv_flag(&_cv_flag, "CVLapackDense", 1)) return false;
       //_cv_flag = CVDlsSetDenseJacFn( _cv_mem, NULL );
       //if ( _check_cv_flag(&_cv_flag, "CVDlsSetDenseJacFn", 1) ) return false;
       //break;

     case Options::CV_DENSE:
       _cv_flag = CVDense( _cv_mem, NV_LENGTH_S( _Nx ) );
       if( _check_cv_flag(&_cv_flag, "CVDense", 1)) return false;
       _cv_flag = CVDlsSetDenseJacFn( _cv_mem, NULL );
       if ( _check_cv_flag(&_cv_flag, "CVDlsSetDenseJacFn", 1) ) return false;
       break;
    }
  }
  
  // Specify the relative and absolute tolerances for states, with
  // different tolerances for entries of ellipsoid shape matrix ODEs
  if( (options.WRAPMIT == Options::ELLIPS) && (options.ETOL < options.ATOL) ){
    for(unsigned i=0; i<NV_LENGTH_S(_NTOLx); i++ )
      NV_Ith_S(_NTOLx,i) = i<NV_LENGTH_S(_NTOLx)-_nx*(_nx+1)/2? options.ATOL: options.ETOL;
    _cv_flag = CVodeSVtolerances( _cv_mem, options.RTOL, _NTOLx );
    if( _check_cv_flag(&_cv_flag, "CVodeSVtolerances", 1) ) return false;   
  }
  else{
    _cv_flag = CVodeSStolerances( _cv_mem, options.RTOL, options.ATOL );
    if( _check_cv_flag(&_cv_flag, "CVodeSStolerances", 1) ) return false;
  }

  // Set maximum number of nonlinear solver iterations
  _cv_flag = CVodeSetMaxNonlinIters( _cv_mem, options.MAXCORR );
  if ( _check_cv_flag(&_cv_flag, "CVodeSetMaxNonlinIters", 1) ) return false;

  // Set maximum number of error test failures
  _cv_flag = CVodeSetMaxErrTestFails( _cv_mem, options.MAXFAIL );
  if ( _check_cv_flag(&_cv_flag, "CVodeSetMaxErrTestFails", 1) ) return false;

  // Set maximum number of convergence test failures
  _cv_flag = CVodeSetMaxConvFails( _cv_mem, options.MAXFAIL );
  if ( _check_cv_flag(&_cv_flag, "CVodeSetMaxConvFails", 1) ) return false;

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

  // Allocate memory for adjoint sensitivity analysis (if applicable)
  //_cv_flag = CVodeAdjInit( _cv_mem, options.ASACHKPT, options.ASAINTERP );
  //if( _check_cv_flag(&_cv_flag, "CVodeAdjInit", 1) ) return false;

  return true;
}

template <typename T, typename PMT, typename PVT> inline bool
ODEBND_SUNDIALS<T,PMT,PVT>::_CC_CVODE_STA
()
{
  // Reinitialize CVode memory block for current time _t and current state _Nx
  _cv_flag = CVodeReInit( _cv_mem, _t, _Nx );
  if( _check_cv_flag(&_cv_flag, "CVodeReInit", 1) ) return false;

  return true;
}

template <typename T, typename PMT, typename PVT> inline bool
ODEBND_SUNDIALS<T,PMT,PVT>::_CC_CVODE_QUAD
()
{
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
( const unsigned np, const T*Ip )
{
  // Initialize bound propagation
  if( !ODEBND_BASE<T,PMT,PVT>::_INI_I_STA( options, np, Ip, options.ETOL )
   || !ODEBND_SUNDIALS<T,PMT,PVT>::_INI_I_STA( np ) )
    return false;
  return true;
}

template <typename T, typename PMT, typename PVT> inline bool
ODEBND_SUNDIALS<T,PMT,PVT>::_INI_I_STA
( const unsigned np )
{
  // Set SUNDIALS state/quadrature arrays
  unsigned cv_Nx_size, cv_Nq_size, cv_Nf_size;
  switch( options.WRAPMIT){
  case Options::NONE:
  case Options::DINEQ:
    cv_Nx_size = 2*_nx;
    cv_Nq_size = 2*_nq;
    cv_Nf_size = 2*_nf;
    break;
  case Options::ELLIPS:
  default:
    cv_Nx_size = _nx*(1+np)+_nx*(_nx+1)/2;
    cv_Nq_size = _nq*(2+np);
    cv_Nf_size = _nf*(2+np);
    break;
  }
  if( !_Nx || cv_Nx_size != NV_LENGTH_S( _Nx ) ){
    if( _Nx ) N_VDestroy_Serial( _Nx );
    _Nx = N_VNew_Serial( cv_Nx_size );
  }
  if( !_Nq || cv_Nq_size != NV_LENGTH_S( _Nq ) ){
    if( _Nq ) N_VDestroy_Serial( _Nq );
    _Nq = cv_Nq_size? N_VNew_Serial( cv_Nq_size ): 0;
  }
  if( !_Nf || cv_Nf_size != NV_LENGTH_S( _Nf ) ){
    if( _Nf ) N_VDestroy_Serial( _Nf );
    _Nf = cv_Nf_size? N_VNew_Serial( cv_Nf_size ): 0;
  }
  if( !_NTOLx || cv_Nx_size != NV_LENGTH_S( _NTOLx ) ){
    if( _NTOLx ) N_VDestroy_Serial( _NTOLx );
      _NTOLx = N_VNew_Serial( cv_Nx_size );
  }
  // Initialize state parameterization at time stages
  _vec_sta.clear();

  // Reset result record and statistics
  results_sta.clear();
  _init_stats( stats_sta );

  return true;
}

template <typename T, typename PMT, typename PVT> inline int
ODEBND_SUNDIALS<T,PMT,PVT>::MC_CVRHSI__
( realtype t, N_Vector y, N_Vector ydot, void *user_data )
{
  ODEBND_SUNDIALS<T,PMT,PVT> *pODEBND = ODEBND_SUNDIALS<T,PMT,PVT>::pODEBND;
#ifdef MC__ODEBNDS_SUNDIALS_DEBUG
  std::cout << "@t=" << t << "\nx:\n";
  for( unsigned i=0; i<NV_LENGTH_S( y ); i++ ) std::cout << NV_Ith_S( y, i ) << std::endl;
#endif
  bool flag = pODEBND->_RHS_I_STA( pODEBND->options, t, NV_DATA_S( y ),
    NV_DATA_S( ydot ) );
#ifdef MC__ODEBNDS_SUNDIALS_DEBUG
  std::cout << "xdot:\n";
  for( unsigned i=0; i<NV_LENGTH_S( ydot ); i++ ) std::cout << NV_Ith_S( ydot, i ) << std::endl;
  { int dum; std::cin >> dum; }
#endif
  ODEBND_SUNDIALS<T,PMT,PVT>::pODEBND = pODEBND;
  pODEBND->stats_sta.numRHS++;
  return( flag? 0: -1 );
}

template <typename T, typename PMT, typename PVT> inline int
ODEBND_SUNDIALS<T,PMT,PVT>::MC_CVQUADI__
( realtype t, N_Vector y, N_Vector qdot, void *user_data )
{
  ODEBND_SUNDIALS<T,PMT,PVT> *pODEBND = ODEBND_SUNDIALS<T,PMT,PVT>::pODEBND;
  bool flag = pODEBND->_RHS_I_QUAD( pODEBND->options, t, NV_DATA_S( y ),
    NV_DATA_S( qdot ) );
  ODEBND_SUNDIALS<T,PMT,PVT>::pODEBND = pODEBND;
  return( flag? 0: -1 );
}

template <typename T, typename PMT, typename PVT>
inline typename ODEBND_SUNDIALS<T,PMT,PVT>::STATUS
ODEBND_SUNDIALS<T,PMT,PVT>::_bounds
( const T*Ip, T**Ixk, T*If, const bool store, std::ostream&os )
{
  // Check size
  if( !Ip || !Ixk || (_nf && !If) ) return FATAL;

  try{
    // Initialize trajectory integration
    if( !_INI_I_STA( _np, Ip ) ) return FATAL;
    _t = _dT[0];

    // Bounds on initial states/quadratures
    if( !_IC_I_SET( options )
     || !_IC_I_STA( options, _t, NV_DATA_S( _Nx ) )
     || (_Nq && !_IC_I_QUAD( options, NV_DATA_S( _Nq ) ) ) )
      { _END_STA(); return FATAL; }
    if( Ixk && !Ixk[0] ) Ixk[0] = new T[_nx];
    for( unsigned ix=0; Ixk[0] && ix<_nx; ix++ ) Ixk[0][ix] = _Ix[ix];

    // Store full state at initial time
    if( store ){
      realtype*vsta = NV_DATA_S(_Nx);
      unsigned lsta = NV_LENGTH_S(_Nx);
      _vec_sta.push_back( std::vector<realtype>( vsta, vsta+lsta ) );
    }

    // Display & record initial results
    if( options.DISPLAY >= 1 ){
      _print_interm( _t, _nx, Ixk[0], "x", os );
    }
    if( options.RESRECORD )
      results_sta.push_back( Results( _dT[0], _nx, Ixk[0] ) );

    // Integrate ODEs through each stage using SUNDIALS
    pODEBND = this;
    if( !_INI_CVODE( MC_CVRHSI__, MC_CVQUADI__ ) )
      { _END_STA(); return FATAL; }

    for( _istg=0; _istg<_nsmax; _istg++ ){
      // Bounds on state discontinuities (if any) at stage times
      // and integrator reinitialization (if applicable)
      _pos_ic = ( _vIC.size()>=_nsmax? _istg:0 );
      if( _pos_ic
       && ( !_CC_I_SET( options, _pos_ic )
         || !_CC_I_STA( options, _t, NV_DATA_S( _Nx ) )
         || !_CC_CVODE_STA() ) )
        { _END_STA(); return FAILURE; }
      if( _istg 
       && ( (_Nq && !_IC_I_QUAD( options, NV_DATA_S( _Nq ) ) ) // quadrature reinitialization
         || !_CC_CVODE_QUAD() ) )
        { _END_STA(); return FAILURE; }

      // update list of operations in RHS, JAC and QUAD
      _pos_rhs  = ( _vRHS.size()<=1? 0: _istg );
      _pos_quad = ( _vQUAD.size()<=1? 0: _istg );
      if( (!_istg || _pos_rhs || _pos_quad)
        && !_RHS_I_SET( options, _pos_rhs, _pos_quad ) )
        { _END_STA(); return FATAL; }

      // integrate till end of time stage
      _cv_flag = CVodeSetStopTime( _cv_mem, _dT[_istg+1] );
      if( _check_cv_flag(&_cv_flag, "CVodeSetStopTime", 1) )
        { _END_STA(); return FATAL; }
      while( _t < _dT[_istg+1] ){
        if( !store )
          _cv_flag = CVode( _cv_mem, _dT[_istg+1], _Nx, &_t, CV_ONE_STEP );
        else
          _cv_flag = CVodeF( _cv_mem, _dT[_istg+1], _Nx, &_t, CV_ONE_STEP, &_nchk );
        if( _check_cv_flag(&_cv_flag, store?"CVodeF":"CVode", 1)
         || (options.NMAX && stats_sta.numSteps > options.NMAX)
         || _diam( _nx, _Ix ) > options.DMAX )
          throw Exceptions( Exceptions::INTERN );
        stats_sta.numSteps++;
      }

      // Store full state at stage time
      if( store ){
        realtype*vsta = NV_DATA_S(_Nx);
        unsigned lsta = NV_LENGTH_S(_Nx);
        _vec_sta.push_back( std::vector<realtype>( vsta, vsta+lsta ) );
      }

      // Bounds on intermediate states and quadratures
      if( _nq ){
        _cv_flag = CVodeGetQuad( _cv_mem, &_t, _Nq );
        if( _check_cv_flag(&_cv_flag, "CVodeGetQuad", 1) )
          { _END_STA(); return FATAL; }
      }
      _GET_I_STA( options, NV_DATA_S(_Nx), _nq && _Nq? NV_DATA_S(_Nq): 0 );

      // Keep track/display/record stage results
      if( Ixk && !Ixk[_istg+1] ) Ixk[_istg+1] = new T[_nx];
      for( unsigned ix=0; Ixk[_istg+1] && ix<_nx; ix++ ) Ixk[_istg+1][ix] = _Ix[ix];
      if( options.DISPLAY >= 1 )
        _print_interm( _t, _nx, Ixk[_istg+1], "x", os );
      if( options.RESRECORD )
        results_sta.push_back( Results( _dT[_istg+1], _nx, Ixk[_istg+1] ) );

      // Add intermediate function terms
      _GET_I_STA( options, NV_DATA_S(_Nx), _nq && _Nq? NV_DATA_S(_Nq): 0 );
      _pos_fct = ( _vFCT.size()>=_nsmax? _istg:0 );
      if( _nf
       && (_vFCT.size()>=_nsmax || _istg==_nsmax-1)
       && !_FCT_I_STA( options, _pos_fct, _t, NV_DATA_S( _Nf ), If ) )
        { _END_STA(); return FATAL; }
    }

    // Bounds on final quadratures and functions
    if( options.DISPLAY >= 1 ) _print_interm( _nf, If, "f", os );
  }
  catch(...){
    _END_STA();
    if( options.DISPLAY >= 1 ){
      os << " ABORT TIME  " << std::scientific << std::left
                            << std::setprecision(5) << _t << std::endl;
      _print_stats( stats_sta, os );
    }
    return FAILURE;
  }

  _END_STA();
  if( options.DISPLAY >= 1 ) _print_stats( stats_sta, os );
  return NORMAL;
}

//! @fn template <typename T, typename PMT, typename PVT> inline typename ODEBND_SUNDIALS<T,PMT,PVT>::STATUS ODEBND_SUNDIALS<T,PMT,PVT>::bounds(
//! const T*Ip, T**Ixk, T*If, std::ostream&os=std::cout )
//!
//! This function computes an interval enclosure of the reachable set of 
//! the parametric ODEs defined in IVP using equally spaced samples:
//!   - <a>Ip</a> [input] interval parameter set
//!   - <a>Ixk</a> [output] interval state enclosures at stage times
//!   - <a>If</a> [output] interval function enclosures
//!   - <a>os</a> [input] output stream [default: std::cout]
//! .
//! The return value is the status.
template <typename T, typename PMT, typename PVT>
inline typename ODEBND_SUNDIALS<T,PMT,PVT>::STATUS
ODEBND_SUNDIALS<T,PMT,PVT>::bounds
( const T*Ip, T**Ixk, T*If, std::ostream&os )
{
  return _bounds( Ip, Ixk, If, false, os );
}

//! @fn template <typename T, typename PMT, typename PVT> inline typename ODEBND_SUNDIALS<T,PMT,PVT>::STATUS ODEBND_SUNDIALS<T,PMT,PVT>::hausdorff(
//! const T*Ip, double**Hxk, const unsigned nsamp, std::ostream&os )
//!
//! This function computes the Hausdorff distance between the interval enclosure
//! and the exact reachable set projected onto each variable:
//!   - <a>Ip</a>    [input]  interval parameter set
//!   - <a>Hxk</a>   [output] Hausdorff distance between the interval enclosure
//!                           and the exact reachable set projected onto each
//!                           variable, at stage times
//!   - <a>nsamp</a> [input]  number of samples for each parameter
//!   - <a>os</a>    [input]  output stream
//! .
//! The return value is the status.
template <typename T, typename PMT, typename PVT>
inline typename ODEBND_SUNDIALS<T,PMT,PVT>::STATUS
ODEBND_SUNDIALS<T,PMT,PVT>::hausdorff
( const T*Ip, double**Hxk, double*Hf, const unsigned nsamp, std::ostream&os )
{
  return _hausdorff( Ip, Hxk, Hf, *this, nsamp, os )? NORMAL: FAILURE;
}

template <typename T, typename PMT, typename PVT> inline bool
ODEBND_SUNDIALS<T,PMT,PVT>::_INI_PM_STA
( const unsigned np, const PVT* PMp )
{
  // Initialize bound propagation
  if( !ODEBND_BASE<T,PMT,PVT>::_INI_PM_STA( options, np, PMp, options.ETOL )
   || !ODEBND_SUNDIALS<T,PMT,PVT>::_INI_PM_STA( np ) )
    return false;
  return true;
}

template <typename T, typename PMT, typename PVT> inline bool
ODEBND_SUNDIALS<T,PMT,PVT>::_INI_PM_STA
( const unsigned np )
{
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
  if( !_NTOLx || cv_Nx_size != NV_LENGTH_S( _NTOLx ) ){
    if( _NTOLx ) N_VDestroy_Serial( _NTOLx );
      _NTOLx = N_VNew_Serial( cv_Nx_size );
  }
  // Initialize state parameterization at time stages
  _vec_sta.clear();

  // Reset result record and statistics
  results_sta.clear();
  _init_stats( stats_sta );

  return true;
}

template <typename T, typename PMT, typename PVT> inline int
ODEBND_SUNDIALS<T,PMT,PVT>::MC_CVRHSPM__
( realtype t, N_Vector y, N_Vector ydot, void *user_data )
{
  ODEBND_SUNDIALS<T,PMT,PVT> *pODEBND = ODEBND_SUNDIALS<T,PMT,PVT>::pODEBND;
#ifdef MC__ODEBNDS_SUNDIALS_DEBUG
  std::cout << "@t=" << t << "\nx:\n";
  for( unsigned i=0; i<NV_LENGTH_S( y ); i++ ) std::cout << NV_Ith_S( y, i ) << std::endl;
#endif
  int flag = pODEBND->_RHS_PM_STA( pODEBND->options, t, NV_DATA_S( y ), NV_DATA_S( ydot ) );
#ifdef MC__ODEBNDS_SUNDIALS_DEBUG
  std::cout << "xdot:\n";
  for( unsigned i=0; i<NV_LENGTH_S( ydot ); i++ ) std::cout << NV_Ith_S( ydot, i ) << std::endl;
  { int dum; std::cin >> dum; }
#endif
  ODEBND_SUNDIALS<T,PMT,PVT>::pODEBND = pODEBND;
  pODEBND->stats_sta.numRHS++;
  if( _diam( pODEBND->_nx, pODEBND->_PMx ) > pODEBND->options.DMAX ) return -1;
  return flag;
}

template <typename T, typename PMT, typename PVT> inline int
ODEBND_SUNDIALS<T,PMT,PVT>::MC_CVQUADPM__
( realtype t, N_Vector y, N_Vector qdot, void *user_data )
{
  ODEBND_SUNDIALS<T,PMT,PVT> *pODEBND = ODEBND_SUNDIALS<T,PMT,PVT>::pODEBND;
  bool flag = pODEBND->_RHS_PM_QUAD( pODEBND->options, t, NV_DATA_S( y ),
    NV_DATA_S( qdot ) );
#ifdef MC__ODEBNDS_SUNDIALS_DEBUG
  std::cout << "qdot:\n";
  for( unsigned i=0; i<NV_LENGTH_S( qdot ); i++ ) std::cout << NV_Ith_S( qdot, i ) << std::endl;
  { int dum; std::cin >> dum; }
#endif
  ODEBND_SUNDIALS<T,PMT,PVT>::pODEBND = pODEBND;
  //pODEBND->stats_sta.numRHS++;
  return( flag? 0: -1 );
}

template <typename T, typename PMT, typename PVT>
inline typename ODEBND_SUNDIALS<T,PMT,PVT>::STATUS
ODEBND_SUNDIALS<T,PMT,PVT>::_bounds
( const PVT*PMp, PVT**PMxk, E*ERxk, PVT*PMf, const bool store, std::ostream&os )
{
  // Check arguments
  if( !PMxk || !PMp || (_nf && !PMf) ) return FATAL;

  try{
    // Initialize trajectory integration
    if( !_INI_PM_STA( _np, PMp ) ) return FATAL;
    _t = _dT[0];

    // Bounds on initial states/quadratures
    if( !_IC_PM_SET( options )
     || !_IC_PM_STA( options, _t, NV_DATA_S( _Nx ) )
     || ( _Nq && !_IC_PM_QUAD( options, NV_DATA_S( _Nq ) ) ) )
      { _END_STA(); return FATAL; }
    if( PMxk && !PMxk[0] ) PMxk[0] = new PVT[_nx];
    for( unsigned ix=0; PMxk[0] && ix<_nx; ix++ ) PMxk[0][ix] = _PMx[ix];
    if( options.WRAPMIT == Options::ELLIPS && ERxk ) ERxk[0] = _Er;

    // Store full state at initial time
    if( store ){
      realtype*vsta = NV_DATA_S(_Nx);
      unsigned lsta = NV_LENGTH_S(_Nx);
      _vec_sta.push_back( std::vector<realtype>( vsta, vsta+lsta ) );
    }

    // Display & record initial results
    if( options.DISPLAY >= 1 )
      _print_interm( _t, _nx, PMxk[0], "x", os );
    if( options.RESRECORD )
      results_sta.push_back( Results( _dT[0], _nx, PMxk[0] ) );

    // Integrate ODEs through each stage using SUNDIALS
    pODEBND = this;
    if( !_INI_CVODE( MC_CVRHSPM__, MC_CVQUADPM__ ) )
      { _END_STA(); return FATAL; }

    for( _istg=0; _istg<_nsmax; _istg++ ){
      // Bounds on state discontinuities (if any) at stage times
      // and integrator reinitialization (if applicable)
      _pos_ic = ( _vIC.size()>=_nsmax? _istg:0 );
      if( _pos_ic
       && ( !_CC_PM_SET( options, _pos_ic )
         || !_CC_PM_STA( options, _t, NV_DATA_S( _Nx ) )
         || !_CC_CVODE_STA() ) )
        { _END_STA(); return FAILURE; }
      if( _istg 
       && ( (_Nq && !_IC_PM_QUAD( options, NV_DATA_S( _Nq ) ) ) // quadrature reinitialization
         || !_CC_CVODE_QUAD() ) )
        { _END_STA(); return FAILURE; }

      // update list of operations in RHS, JAC and QUAD
      _pos_rhs  = ( _vRHS.size()<=1? 0: _istg );
      _pos_quad = ( _vQUAD.size()<=1? 0: _istg );
      if( (!_istg || _pos_rhs || _pos_quad)
        && !_RHS_PM_SET( options, _pos_rhs, _pos_quad ) )
        { _END_STA(); return FATAL; }

      // integrate till end of time stage
      _cv_flag = CVodeSetStopTime( _cv_mem, _dT[_istg+1] );
      if( _check_cv_flag(&_cv_flag, "CVodeSetStopTime", 1) )
        { _END_STA(); return FATAL; }
      while( _t < _dT[_istg+1] ){
        if( !store )
          _cv_flag = CVode( _cv_mem, _dT[_istg+1], _Nx, &_t, CV_ONE_STEP );
        else
          _cv_flag = CVodeF( _cv_mem, _dT[_istg+1], _Nx, &_t, CV_ONE_STEP, &_nchk );
        if( _check_cv_flag(&_cv_flag, store?"CVodeF":"CVode", 1)
         || (options.NMAX && stats_sta.numSteps > options.NMAX)
         || _diam( _nx, _PMx ) > options.DMAX
         || _diam( _nx, _Ir  ) > options.DMAX )
          throw Exceptions( Exceptions::INTERN );
        stats_sta.numSteps++;
      }

      // Store full state at stage time
      if( store ){
        realtype*vsta = NV_DATA_S(_Nx);
        unsigned lsta = NV_LENGTH_S(_Nx);
        _vec_sta.push_back( std::vector<realtype>( vsta, vsta+lsta ) );
      }

      // Bounds on intermediate states and quadratures
      if( _nq ){
        _cv_flag = CVodeGetQuad( _cv_mem, &_t, _Nq );
        if( _check_cv_flag(&_cv_flag, "CVodeGetQuad", 1) )
          { _END_STA(); return FATAL; }
      }
     _GET_PM_STA( options, NV_DATA_S(_Nx), _nq && _Nq? NV_DATA_S(_Nq): 0 );

      // Keep track/display/record stage results
      if( PMxk && !PMxk[_istg+1] ) PMxk[_istg+1] = new PVT[_nx];
      for( unsigned ix=0; PMxk[_istg+1] && ix<_nx; ix++ ) PMxk[_istg+1][ix] = _PMx[ix];
      if( options.WRAPMIT == Options::ELLIPS && ERxk ) ERxk[_istg+1] = _Er;
      if( options.DISPLAY >= 1 )
        _print_interm( _t, _nx, PMxk[_istg+1], "x", os );
      if( options.RESRECORD )
        results_sta.push_back( Results( _dT[_istg+1], _nx, PMxk[_istg+1] ) );

      // Add intermediate function terms
      _pos_fct = ( _vFCT.size()>=_nsmax? _istg:0 );
      if( _nf
       && (_vFCT.size()>=_nsmax || _istg==_nsmax-1)
       && !_FCT_PM_STA( _pos_fct, _t, PMf ) )
        { _END_STA(); return FATAL; }
    }

    // Bounds on final quadratures and functions
    if( options.DISPLAY >= 1 ) _print_interm( _nf, PMf, "f", os );
  }
  catch(...){
    _END_STA();
    if( options.DISPLAY >= 1 ){
      os << " ABORT TIME  " << std::scientific << std::left
                            << std::setprecision(5) << _t << std::endl;
      _print_stats( stats_sta, os );
    }
    return FAILURE;
  }

  _END_STA();
  if( options.DISPLAY >= 1 ) _print_stats( stats_sta, os );
  return NORMAL;
}

//! @fn template <typename T, typename PMT, typename PVT> inline typename ODEBND_SUNDIALS<T,PMT,PVT>::STATUS ODEBND_SUNDIALS<T,PMT,PVT>::bounds(
//! const PVT*PMp, PVT**PMxk, PVT*PMf, std::ostream&os=std::cout )
//!
//! This function computes an enclosure of the reachable set of the parametric ODEs
//! using propagation of polynomial models with convex remainders (intervals, ellipsoids):
//!   - <a>PMp</a>  [input]  polynomial model of parameter set
//!   - <a>PMxk</a> [output] polynomial model of state varibles at stage times
//!   - <a>PMf</a>  [output] polynomial model of state-dependent functions
//!   - <a>os</a>   [input]  output stream [default: std::cout]
//! .
//! The return value is the status.
template <typename T, typename PMT, typename PVT>
inline typename ODEBND_SUNDIALS<T,PMT,PVT>::STATUS
ODEBND_SUNDIALS<T,PMT,PVT>::bounds
( const PVT*PMp, PVT**PMxk, PVT*PMf, std::ostream&os )
{
  return _bounds( PMp, PMxk, 0, PMf, false, os );
}

//! @fn template <typename T, typename PMT, typename PVT> inline typename ODEBND_SUNDIALS<T,PMT,PVT>::STATUS ODEBND_SUNDIALS<T,PMT,PVT>::bounds(
//! const PVT*PMp, PVT**PMxk, E*ERxk, PVT*PMf, std::ostream&os=std::cout )
//!
//! This function computes an enclosure of the reachable set of the parametric ODEs
//! using propagation of polynomial models with convex remainders (intervals, ellipsoids):
//!   - <a>PMp</a>  [input]  polynomial model of parameter set
//!   - <a>PMxk</a> [output] polynomial model of state variables at stage times
//!   - <a>ERxk</a> [output] ellipsoidal remainder of state variables at stage times (only if ellipsoidal bounder is selected)
//!   - <a>PMf</a>  [output] polynomial model of state-dependent functions
//!   - <a>os</a>   [input]  output stream [default: std::cout]
//! .
//! The return value is the status.
template <typename T, typename PMT, typename PVT>
inline typename ODEBND_SUNDIALS<T,PMT,PVT>::STATUS
ODEBND_SUNDIALS<T,PMT,PVT>::bounds
( const PVT*PMp, PVT**PMxk, E*ERxk, PVT*PMf, std::ostream&os )
{
  return _bounds( PMp, PMxk, ERxk, PMf, false, os );
}

//! @fn template <typename T, typename PMT, typename PVT> inline typename ODEBND_SUNDIALS<T,PMT,PVT>::STATUS ODEBND_SUNDIALS<T,PMT,PVT>::hausdorff(
//! const PVT*PMp, double**Hxk, double*Hf, const unsigned nsamp, std::ostream&os=std::cout )
//!
//! This function computes the Hausdorff distance between the polynomial model
//! remainder and the actual (sampled) range of the remainder function
//! in projection onto each variable and for each stage time
//! remainder and the actual range of the remainder function:
//!   - <a>PMp</a>   [input]  polynomial model of parameter set
//!   - <a>Hxk</a>   [output] Hausdorff distance between the polynomial model
//!                           remainder and the actual (sampled) range of the
//!                           remainder term for states at stage times
//!   - <a>Hf</a>    [output] Hausdorff distance between the polynomial model
//!                           remainder and the actual (sampled) range of the
//!                           remainder term for state-dependent functions
//!   - <a>nsamp</a> [input]  number of samples for each parameter
//!   - <a>os</a>    [input]  output stream [default: std::cout]

template <typename T, typename PMT, typename PVT>
inline typename ODEBND_SUNDIALS<T,PMT,PVT>::STATUS
ODEBND_SUNDIALS<T,PMT,PVT>::hausdorff
( const PVT*PMp, double**Hxk, double*Hf, const unsigned nsamp, std::ostream&os )
{
  pODESLV.set( *this );
  pODESLV.options = options.ODESLV;
  return _hausdorff( PMp, Hxk, Hf, *this, pODESLV, nsamp, os )?
         NORMAL: FAILURE;
}

//! @fn template <typename T, typename PMT, typename PVT> inline typename ODEBND_SUNDIALS<T,PMT,PVT>::STATUS ODEBND_SUNDIALS<T,PMT,PVT>::bounds(
//! const T*Ip, T**Ixk, T*If, const unsigned nsamp, std::ostream&os=std::cout )
//!
//! This function computes projections of an inner-approximation enclosure of
//! the reachable set of the parametric ODEs using sampling and continuous-time
//! integration:
//!   - <a>Ip</a>    [input] interval parameter set
//!   - <a>Ixk</a>   [output] approximate interval state enclosures at stage times
//!   - <a>If</a>    [output] approximate state-dependent function enclosures
//!   - <a>nsamp</a> [input] number of samples for each parameter
//!   - <a>os</a>    [input] output stream
//! .
//! The return value is the status.
template <typename T, typename PMT, typename PVT>
inline typename ODEBND_SUNDIALS<T,PMT,PVT>::STATUS
ODEBND_SUNDIALS<T,PMT,PVT>::bounds
( const T*Ip, T**Ixk, T*If, const unsigned nsamp, std::ostream&os )
{
  // Sample inner approximation
  STATUS flag = NORMAL;
  pODESLV.set( *this );
  pODESLV.options = options.ODESLV;
  if( !_bounds( Ip, Ixk, If, pODESLV, nsamp, os ) )
    flag = FAILURE;

  // Display results
  if( options.DISPLAY >= 1 ){
    for( unsigned is=0; Ixk && is<=_nsmax; is++ )
      _print_interm( _dT[is], _nx, Ixk[is], "x", os );
    if( If ) _print_interm( _nf, If, "f", os );
  }

  // Record intermediate results
  results_sta.clear();
  if( options.RESRECORD )
    for( unsigned is=0; Ixk && is<=_nsmax; is++ )
      results_sta.push_back( Results( _dT[is], _nx, Ixk[is] ) );
  
  return flag;
}

} // end namescape mc

#endif

