// Copyright (C) 2015 Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

#ifndef MC__ODEBND_GSL_HPP
#define MC__ODEBND_GSL_HPP

#undef  MC__ODEBND_GSL_DINEQI_DEBUG
#undef  MC__ODEBND_GSL_DINEQPM_DEBUG

#include "base_de.hpp"
#include "odebnd_base.hpp"
#include "base_gsl.hpp"
#include "mesh_gsl.hpp"

#include "tmodel.hpp"
#include "cmodel.hpp"

namespace mc
{
//! @brief C++ class computing enclosures of the reachable set of parametric ODEs using continuous-time set-valued integration.
////////////////////////////////////////////////////////////////////////
//! mc::ODEBND_GSL is a C++ class that computes enclosures of the
//! reachable set of parametric ordinary differential equations
//! (ODEs) using continuous-time set-valued integration. It implements
//! the methods of differential inequalities, whereby polynomial models
//! with interval or ellipsoidal remainders are used to enable high-
//! order convergence. The use of ellipsoidal remainders enables
//! stability of the enclosures for asymptotically stable ODE systems
//! when the parameter host is sufficiently small. The numerical
//! integrator is gsl_odeiv2 in GSL.
////////////////////////////////////////////////////////////////////////
template <typename T, typename PMT=mc::TModel<T>, typename PVT=mc::TVar<T> >
class ODEBND_GSL:
  public virtual BASE_DE,
  public virtual BASE_GSL,
  public virtual ODEBND_BASE<T,PMT,PVT>
{
 protected:
  typedef Ellipsoid E;

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
  using ODEBND_BASE<T,PMT,PVT>::_hausdorff;
  using ODEBND_BASE<T,PMT,PVT>::_remainders;
  using ODEBND_BASE<T,PMT,PVT>::_print_interm;

  template <typename U, typename PMU, typename PVU>
  friend int MC_GSLRHSI__
    ( double t, const double* x, double* xdot, void* user_data );
  template <typename U, typename PMU, typename PVU>
  friend int MC_GSLJACI__
    ( double t, const double* x, double* jac, double* xdot, void* user_data );
  template <typename U, typename PMU, typename PVU>
  friend int MC_GSLRHSPM__
    ( double t, const double* x, double* xdot, void* user_data );
  template <typename U, typename PMU, typename PVU>
  friend int MC_GSLJACPM__
    ( double t, const double* x, double* jac, double* xdot, void* user_data );

 private:
  //! @brief GSL data type for bounding ODE system
  gsl_odeiv2_system _sys_sta;

  //! @brief GSL driver for bounding ODE system
  gsl_odeiv2_driver *_driver_sta;

 protected:
  //! @brief full GSL state
  double *_vec_sta;

  //! @brief pointer to quadratures in full GSL state
  double *_vec_quad;

  //! @brief offset to quadratures in full GSL state
  unsigned _offset_quad;

  //! @brief stepsize
  double _h;

  //! @brief Mesh storing state bound parameterizations 
  MESH_GSL _mesh_sta;

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
  static ODEBND_GSL<T,PMT,PVT> *_pODEBND;

public:
  typedef typename ODEBND_BASE<T,PMT,PVT>::Results Results;

  /** @defgroup ODEBND_GSL Continuous-time set-valued integration of parametric ODEs
   *  @{
   */
  //! @brief Default constructor
  ODEBND_GSL();

  //! @brief Virtual destructor
  virtual ~ODEBND_GSL();

  //! @brief Integrator options
  struct Options: public BASE_GSL::Options
  {
    //! @brief Constructor
    Options():
      BASE_GSL::Options(), WRAPMIT(ELLIPS), ORDMIT(2), PMOPT(typename PMT::Options()),
      QTOL(machprec()), QSCALE(true), PMNOREM(false), USEINV(true), MESHPREALLOC(0),
      DMAX(1e20), DISPLAY(1), RESRECORD(false)
      {}
    //! @brief Assignment operator
    Options& operator=
      ( Options&options ){
        BASE_GSL::Options::operator=(options);
        WRAPMIT   = options.WRAPMIT;
        ORDMIT    = options.ORDMIT;
        PMOPT     = options.PMOPT;
        QTOL      = options.QTOL;
        QSCALE    = options.QSCALE;
        PMNOREM   = options.PMNOREM;
        USEINV    = options.USEINV;
        MESHPREALLOC = options.MESHPREALLOC;
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
    //! @brief Whether or not to scale the ellispoidal remainder (Default: true)
    double QSCALE;
    //! @brief Whether or not to cancel the remainder term in polynomial models (non-validated integration; Default: false)
    bool PMNOREM;
    //! @brief Whether or not to use the specified invariants for bounds contraction (Default: true)
    bool USEINV;
    //! @brief Preallocated mesh size (Default: 0)
    bool MESHPREALLOC;
    //! @brief Maximum enclosure diameter, \f$D_{\rm max}\f$ (Default: 1e20)
    double DMAX;
    //! @brief Display level (Default: 1)
    int DISPLAY;
    //! @brief Whether or not to record results (Default: false)
    bool RESRECORD;
  } options;

  //! @brief Structure for setting up storing the solver exceptions
  class Exceptions
  {
  public:
    //! @brief Enumeration type for ODEBND_GSL exception handling
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
        return "ODEBND_GSL::Exceptions  Internal error";
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
    ( const unsigned ns, const double*tk, const T*Ip, T**Ixk=0,
      T*Iq=0, T*If=0, std::ostream&os=std::cout );

  //! @brief Computes polynomial model enclosure of reachable set of parametric ODEs
  STATUS bounds
    ( const unsigned ns, const double*tk, const PVT*PMp, PVT**PMxk=0,
      PVT*PMq=0, PVT*PMf=0, std::ostream&os=std::cout );

  //! @brief Computes Hausdorff distance between interval enclosure and actual reachable set of parametric ODEs, using parameter sampling
  template<typename ODESLV> STATUS hausdorff
    ( const unsigned ns, const double*tk, const T*Ip, double**Hxk,
      ODESLV&traj, const unsigned nsamp, std::ostream&os=std::cout );

  //! @brief Computes Hausdorff distance between polynomial model remainder enclosure and actual remainder function range, using parameter sampling
  template<typename ODESLV> STATUS hausdorff
    ( const unsigned ns, const double*tk, const PVT*PMp, double**Hxk,
      ODESLV&traj, const unsigned nsamp, std::ostream&os=std::cout );

  //! @brief Record results in file <a>bndrec</a>, with accuracy of <a>iprec</a> digits
  void record
    ( std::ofstream&bndrec, const unsigned iprec=5 ) const
    { return ODEBND_BASE<T,PMT,PVT>::_record( bndrec, results_sta, iprec ); }
  /** @} */

protected:

  //! @brief Function to initialize GSL driver
  void _INI_GSL
    ( gsl_odeiv2_system &sys, gsl_odeiv2_driver *&driver );

  //! @brief Function to finalize state bounding
  void _END_STA();

  //! @brief Function to initialize state interval bounding
  bool _INI_I_STA
    ( const unsigned np, const T*Ip, const unsigned ns );

  //! @brief Static wrapper to function to calculate the DINEQs RHS values
  static int MC_GSLRHSI__
    ( double t, const double* x, double* xdot, void* user_data );

  //! @brief Static wrapper to function to calculate the DINEQs RHS derivatives
  static int MC_GSLJACI__
    ( double t, const double* x, double* jac, double* xdot, void* user_data );

  //! @brief Propagate state/quadrature interval bounds forward in time through every time stages
  STATUS _bounds
    ( const unsigned ns, const double*tk, const T*Ip, T**Ixk,
      T*Iq, T*If, const bool store, std::ostream&os );

  //! @brief Function to initialize GSL for state polynomial models
  bool _INI_PM_STA
    ( const unsigned np, const PVT*PMp, const unsigned ns );

  //! @brief Static wrapper to function to calculate the DINEQ-PMs RHS values
  static int MC_GSLRHSPM__
    ( double t, const double*x, double*xdot, void*user_data );

  //! @brief Static wrapper to function to calculate the DINEQ-PMs RHS derivatives
  static int MC_GSLJACPM__
    ( double t, const double*x, double*jac, double*xdot, void*user_data );

  //! @brief Propagate state/quadrature polynomial models forward in time through every time stages
  STATUS _bounds
    ( const unsigned ns, const double*tk, const PVT*PMp, PVT**PMxk,
      PVT*PMq, PVT*PMf, const bool store, std::ostream&os );

  //! @brief Private methods to block default compiler methods
  ODEBND_GSL(const ODEBND_GSL&);
  ODEBND_GSL& operator=(const ODEBND_GSL&);
};

template <typename T, typename PMT, typename PVT>
 ODEBND_GSL<T,PMT,PVT>* ODEBND_GSL<T,PMT,PVT>::_pODEBND = 0;

template <typename T, typename PMT, typename PVT> inline
ODEBND_GSL<T,PMT,PVT>::ODEBND_GSL
()
: BASE_DE(), BASE_GSL(), ODEBND_BASE<T,PMT,PVT>(), _driver_sta(0),
  _vec_sta(0), _vec_quad(0)
{}

template <typename T, typename PMT, typename PVT> inline
ODEBND_GSL<T,PMT,PVT>::~ODEBND_GSL
()
{
  delete[] _vec_sta;
  if( _driver_sta ) gsl_odeiv2_driver_free( _driver_sta );
}

template <typename T, typename PMT, typename PVT> inline void
ODEBND_GSL<T,PMT,PVT>::_INI_GSL
( gsl_odeiv2_system &sys, gsl_odeiv2_driver *&driver )
{
  // Set GSL driver
  if( driver ) gsl_odeiv2_driver_free( driver );
  switch( options.INTMETH ){
  case Options::RK8PD:
    driver = gsl_odeiv2_driver_alloc_y_new( &sys, gsl_odeiv2_step_rk8pd,
      options.H0, options.ATOL, options.RTOL );
    break;
  case Options::MSADAMS:
    driver = gsl_odeiv2_driver_alloc_y_new( &sys, gsl_odeiv2_step_msadams,
      options.H0, options.ATOL, options.RTOL );
    break;
  case Options::MSBDF:
    driver = gsl_odeiv2_driver_alloc_y_new( &sys, gsl_odeiv2_step_msbdf,
      options.H0, options.ATOL, options.RTOL );
    break;
  case Options::RKF45: default:
    driver = gsl_odeiv2_driver_alloc_y_new( &sys, gsl_odeiv2_step_rkf45,
      options.H0, options.ATOL, options.RTOL );
    break;
  }
  gsl_odeiv2_driver_set_hmin( driver, options.HMIN );  
  gsl_odeiv2_driver_set_nmax( driver, options.NMAX );  
}

template <typename T, typename PMT, typename PVT> inline void
ODEBND_GSL<T,PMT,PVT>::_END_STA()
{
  // Get final CPU time
  _final_stats( stats_sta );
}

template <typename T, typename PMT, typename PVT> inline int
ODEBND_GSL<T,PMT,PVT>::MC_GSLRHSI__
( double t, const double* x, double* xdot, void* user_data )
{
  ODEBND_GSL<T,PMT,PVT> *pODEBND = ODEBND_GSL<T,PMT,PVT>::_pODEBND;
  bool flag = pODEBND->_RHS_I_STA( pODEBND->options, t, x, xdot );
  if( flag && pODEBND->_nq ){
    double* qdot = xdot + pODEBND->_offset_quad;
    flag = pODEBND->_RHS_I_QUAD( pODEBND->options, t, x, qdot );
  }
  pODEBND->stats_sta.numRHS++;
  ODEBND_GSL<T,PMT,PVT>::_pODEBND = pODEBND;
  return( flag? GSL_SUCCESS: GSL_EBADFUNC );
}

template <typename T, typename PMT, typename PVT> inline int
ODEBND_GSL<T,PMT,PVT>::MC_GSLJACI__
( double t, const double* x, double* jac, double* xdot, void* user_data )
{
  ODEBND_GSL<T,PMT,PVT> *pODEBND = ODEBND_GSL<T,PMT,PVT>::_pODEBND;
  bool flag = pODEBND->_JAC_I_STA( pODEBND->options, t, x, jac, xdot );
  pODEBND->stats_sta.numJAC++;
  ODEBND_GSL<T,PMT,PVT>::_pODEBND = pODEBND;
  return( flag? GSL_SUCCESS: GSL_EBADFUNC );
}

template <typename T, typename PMT, typename PVT> inline bool
ODEBND_GSL<T,PMT,PVT>::_INI_I_STA
( const unsigned np, const T*Ip, const unsigned ns )
{
  // Initialize bound propagation
  if( !ODEBND_BASE<T,PMT,PVT>::_INI_I_STA( options, np, Ip, ns ) )
    return false;

  // Set GSL driver
  _sys_sta.function = MC_GSLRHSI__;
  _sys_sta.jacobian = MC_GSLJACI__;
  _sys_sta.params = 0;
  _sys_sta.dimension = 2*_nq;
  switch( options.WRAPMIT){
  case Options::NONE:
  case Options::DINEQ:
    _sys_sta.dimension += 2*_nx;
    break;
  case Options::ELLIPS:
  default:
    _sys_sta.dimension += _nx*(1+np)+_nx*(_nx+1)/2;
    break;
  }
  _INI_GSL( _sys_sta, _driver_sta );
  delete [] _vec_sta;
  _vec_sta  = new double[ _sys_sta.dimension ];
  _offset_quad = _sys_sta.dimension - 2*_nq;
  _vec_quad = _vec_sta + _offset_quad;

  // Initialize mesh
  if( !_mesh_sta.set( ns, _offset_quad,
    options.MESHPREALLOC, options.WRAPMIT ) ) return false;

  // Reset result record and statistics
  results_sta.clear();
  _init_stats( stats_sta );

  return true;
}

template <typename T, typename PMT, typename PVT>
inline typename ODEBND_GSL<T,PMT,PVT>::STATUS
ODEBND_GSL<T,PMT,PVT>::_bounds
( const unsigned ns, const double*tk, const T*Ip, T**Ixk, T*Iq, T*If,
  const bool store, std::ostream&os )
{
  // Check size
  if( !tk || !Ixk || !Ip ) return FATAL;

  try{
    // Initialize trajectory integration with GSL
    _INI_I_STA( _np, Ip, ns );

    // Bounds on initial states/quadratures
    _t = tk[0];
    if( !_IC_I_STA( options, _vec_sta )
     || !_IC_I_QUAD( options, _vec_quad ) )
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

    // Integrate ODEs through each stage using GSL
    _h = options.H0;
    _pODEBND = this;

    for( _istg=0; _istg<ns; _istg++ ){
      // Bounds on state discontinuities (if any) at stage times
      // and integrator reinitialization (if applicable)
      _pos_ic = ( _vIC.size()>=ns? _istg:0 );
      if( _pos_ic ){
        if( !_CC_I_STA( options, _pos_ic, _t, _vec_sta ) )
          { _END_STA(); return FATAL; }
        gsl_odeiv2_driver_reset( _driver_sta );
      }

      // update list of operations in RHS, JAC and QUAD
      _pos_rhs  = ( _vRHS.size()<=1? 0: _istg );
      _pos_quad = ( _vQUAD.size()<=1? 0: _istg );
      if( (!_istg || _pos_rhs || _pos_quad)
        && !_SET_I_STA( options, _pos_rhs, _pos_quad ) )
        { _END_STA(); return FATAL; }

      // Store mesh point at time stage
      if( store && !_mesh_sta.add( _istg, _t, _vec_sta, true ) )
        { _END_STA(); return FATAL; }

      // integrate till end of time stage
      while( _t < tk[_istg+1] ){
        if( gsl_odeiv2_evolve_apply( _driver_sta->e, _driver_sta->c, _driver_sta->s,
            &_sys_sta, &_t, tk[_istg+1], &_h, _vec_sta ) != GSL_SUCCESS
         || _h < options.HMIN
         || (options.NMAX && stats_sta.numSteps > options.NMAX)
         || _diam( _nx, _Ix ) > options.DMAX )
          throw Exceptions( Exceptions::INTERN );
        stats_sta.numSteps++;
        if( options.HMAX > 0 && _h > options.HMAX ) _h = options.HMAX;

        // Store mesh point at intermediate time
        if( store && !_mesh_sta.add( _istg, _t, _vec_sta, false ) )
          { _END_STA(); return FATAL; }
      }

      // Bounds on intermediate states/quadratures
      switch( options.WRAPMIT){
      case Options::NONE:
      case Options::DINEQ:
        _vec2I( _vec_sta, _nx, _Ix );
        break;
      case Options::ELLIPS:
      default:
        _vec2E( _vec_sta, _nx, _np, _Q, _Er, _Ir, _pref, _Ip, _B, _xref, _Ix );
        break;
      }
      _vec2I( _vec_quad, _nq, _Iq );
      if( options.DISPLAY >= 1 ){
        _print_interm( _t, _nx, _Ix, "x", os );
        _print_interm( _nq, _Iq, "q", os );
      }
      if( Ixk && !Ixk[_istg+1] ) Ixk[_istg+1] = new T[_nx];
      for( unsigned ix=0; Ixk[_istg+1] && ix<_nx; ix++ ) Ixk[_istg+1][ix] = _Ix[ix];

      // Record intermediate results
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

//! @fn template <typename T, typename PMT, typename PVT> inline typename ODEBND_GSL<T,PMT,PVT>::STATUS ODEBND_GSL<T,PMT,PVT>::bounds(
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
//!   - <a>If</a> [output] interval functional enclosures (default: NULL)
//!   - <a>os</a> [input] output stream (default: std::cout)
//! .
//! The return value is the status.
template <typename T, typename PMT, typename PVT>
inline typename ODEBND_GSL<T,PMT,PVT>::STATUS
ODEBND_GSL<T,PMT,PVT>::bounds
( const unsigned ns, const double*tk, const T*Ip, T**Ixk, T*Iq, T*If,
  std::ostream&os )
{
  return _bounds( ns, tk, Ip, Ixk, Iq, If, false, os );
}

//! @fn template <typename T, typename PMT, typename PVT> template<typename ODESLV> inline typename ODEBND_GSL<T,PMT,PVT>::STATUS ODEBND_GSL<T,PMT,PVT>::hausdorff(
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
template <typename T, typename PMT, typename PVT> template<typename ODESLV> inline typename ODEBND_GSL<T,PMT,PVT>::STATUS
ODEBND_GSL<T,PMT,PVT>::hausdorff
( const unsigned ns, const double*tk, const T*Ip, double**Hxk,
  ODESLV&traj, const unsigned nsamp, std::ostream&os )
{
  return _hausdorff( ns, tk, Ip, Hxk, *this, traj, nsamp, os )? NORMAL: FAILURE;
}

template <typename T, typename PMT, typename PVT> inline int
ODEBND_GSL<T,PMT,PVT>::MC_GSLRHSPM__
( double t, const double* x, double* xdot, void* user_data )
{
  ODEBND_GSL<T,PMT,PVT> *pODEBND = ODEBND_GSL<T,PMT,PVT>::_pODEBND;
  bool flag = pODEBND->_RHS_PM_STA( pODEBND->options, t, x, xdot );
  if( flag && pODEBND->_nq ){
    double* qdot = xdot + pODEBND->_offset_quad;
    flag = pODEBND->_RHS_PM_QUAD( pODEBND->options, t, x, qdot );
  }
  pODEBND->stats_sta.numRHS++;
  ODEBND_GSL<T,PMT,PVT>::_pODEBND = pODEBND;
  return( flag? GSL_SUCCESS: GSL_EBADFUNC );
}

template <typename T, typename PMT, typename PVT> inline int
ODEBND_GSL<T,PMT,PVT>::MC_GSLJACPM__
( double t, const double* x, double* jac, double* xdot, void* user_data )
{
  ODEBND_GSL<T,PMT,PVT> *pODEBND = ODEBND_GSL<T,PMT,PVT>::_pODEBND;
  bool flag = pODEBND->_JAC_PM_STA( pODEBND->options, t, x, jac, xdot );
  pODEBND->stats_sta.numJAC++;
  ODEBND_GSL<T,PMT,PVT>::_pODEBND = pODEBND;
  return( flag? GSL_SUCCESS: GSL_EBADFUNC );
}

template <typename T, typename PMT, typename PVT> inline bool
ODEBND_GSL<T,PMT,PVT>::_INI_PM_STA
( const unsigned np, const PVT* PMp, const unsigned ns )
{

  // Initialize bound propagation
  if( !ODEBND_BASE<T,PMT,PVT>::_INI_PM_STA( options, np, PMp, ns ) )
    return false;

  // Define ODE system in GSL format
  _sys_sta.function = MC_GSLRHSPM__;
  _sys_sta.jacobian = MC_GSLJACPM__;
  _sys_sta.params = 0;
  _sys_sta.dimension = _PMenv->nmon()*(_nx+_nq)+_nq;
  switch( options.WRAPMIT){
  case Options::NONE:
    _sys_sta.dimension += _nx;
    break;
  case Options::DINEQ:
    _sys_sta.dimension += 2*_nx;
    break;
  case Options::ELLIPS:
  default:
    _sys_sta.dimension += _nx*(_nx+1)/2;
    break;
  }
  _INI_GSL( _sys_sta, _driver_sta );
  delete [] _vec_sta;
  _vec_sta = new double[ _sys_sta.dimension ];
  _offset_quad = _sys_sta.dimension - _PMenv->nmon()*_nq-_nq;
  _vec_quad = _vec_sta + _offset_quad;

  // Initialize mesh
  if( !_mesh_sta.set( ns, _offset_quad,
    options.MESHPREALLOC, PMOFFSET+(int)options.WRAPMIT ) ) return false;

  // Reset result record and statistics
  results_sta.clear();  
  _init_stats( stats_sta );

  return true;
}

template <typename T, typename PMT, typename PVT>
inline typename ODEBND_GSL<T,PMT,PVT>::STATUS
ODEBND_GSL<T,PMT,PVT>::_bounds
( const unsigned ns, const double*tk, const PVT*PMp, PVT**PMxk,
  PVT*PMq, PVT*PMf, const bool store, std::ostream&os )
{
  // Check arguments
  if( !tk || !PMxk || !PMp || (_nf && !PMf) ) return FATAL;

  try{
    // Initialize trajectory integration with GSL
    if( !_INI_PM_STA( _np, PMp, ns ) ) return FATAL;

    // Bounds on initial states/quadratures
    _t = tk[0];
    if( !_IC_PM_STA( options, _vec_sta )
     || !_IC_PM_QUAD( options, _vec_quad ) )
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

    // Integrate ODEs through each stage using GSL
    _h = options.H0;
    _pODEBND = this;

    for( _istg=0; _istg<ns; _istg++ ){
      // Bounds on state discontinuities (if any) at stage times
      // and integrator reinitialization (if applicable)
      _pos_ic = ( _vIC.size()>=ns? _istg:0 );
      if( _pos_ic ){
        if( !_CC_PM_STA( options, _pos_ic, _t, _vec_sta ) )
          { _END_STA(); return FATAL; }
        gsl_odeiv2_driver_reset( _driver_sta );
      }

      // update list of operations in RHS, JAC and QUAD
      _pos_rhs = ( _vRHS.size()<=1? 0: _istg );
      _pos_quad = ( _vQUAD.size()<=1? 0: _istg );
      if( (!_istg || _pos_rhs || _pos_quad)
        && !_SET_PM_STA( options, _pos_rhs, _pos_quad ) )
        { _END_STA(); return FATAL; }

      // Store mesh point at time stage
      if( store && !_mesh_sta.add( _istg, _t, _vec_sta, true ) )
        { _END_STA(); return FATAL; }

      // integrate till end of time stage
      while( _t < tk[_istg+1] ){
        if( gsl_odeiv2_evolve_apply( _driver_sta->e, _driver_sta->c, _driver_sta->s,
            &_sys_sta, &_t, tk[_istg+1], &_h, _vec_sta ) != GSL_SUCCESS
         || _h < options.HMIN
         || (options.NMAX && stats_sta.numSteps > options.NMAX)
         || _diam( _nx, _PMx ) > options.DMAX
         || _diam( _nx, _Ir  ) > options.DMAX )
          throw Exceptions( Exceptions::INTERN );
        stats_sta.numSteps++;
        if( options.HMAX > 0 && _h > options.HMAX ) _h = options.HMAX;

        // Store mesh point at intermediate time
        if( store && !_mesh_sta.add( _istg, _t, _vec_sta, false ) )
          { _END_STA(); return FATAL; }
      }

      // Bounds on intermediate states/quadratures
      switch( options.WRAPMIT){
      case Options::NONE:
        _vec2PMI( _vec_sta,  _PMenv, _nx, _PMx, true );
        break;
      case Options::DINEQ:
        _vec2PMI( _vec_sta,  _PMenv, _nx, _PMx, false );
        break;
      case Options::ELLIPS:
      default:
        _vec2PME( _vec_sta, _PMenv, _nx, _PMx, _Q, _Er, _Ir );
        //std::cout << _Er.eigQ().first;
        break;
      }
      _vec2PMI( _vec_quad, _PMenv, _nq, _PMq, true );

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

//! @fn template <typename T, typename PMT, typename PVT> inline typename ODEBND_GSL<T,PMT,PVT>::STATUS ODEBND_GSL<T,PMT,PVT>::bounds(
//! const unsigned ns, const double*tk, const PVT*PMp, PVT**PMxk=0, 
//! PVT*PMq=0, PVT*PMf=0, std::ostream&os=std::cout )
//!
//! This function computes an enclosure of the reachable set of the parametric ODEs
//! using propagation of polynomial models with convex remainders (intervals, ellipsoids):
//!   - <a>ns</a> [input] number of time stages
//!   - <a>tk</a> [input] stage times, including the initial time
//!   - <a>PMp</a> [input] parameter polynomial models
//!   - <a>PMxk</a> [output] state polynomial models at stage times (default: NULL)
//!   - <a>PMq</a> [output] quadrature polynomial models (default: NULL)
//!   - <a>PMf</a> [output] functional polynomial models (default: NULL)
//!   - <a>os</a> [input] output stream (default: std::cout)
//! .
//! The return value is the status.
template <typename T, typename PMT, typename PVT>
inline typename ODEBND_GSL<T,PMT,PVT>::STATUS
ODEBND_GSL<T,PMT,PVT>::bounds
( const unsigned ns, const double*tk, const PVT*PMp, PVT**PMxk,
  PVT*PMq, PVT*PMf, std::ostream&os )
{
  return _bounds( ns, tk, PMp, PMxk, PMq, PMf, false, os );
}

//! @fn template <typename T, typename PMT, typename PVT> template<typename ODESLV> inline typename ODEBND_GSL<T,PMT,PVT>::STATUS ODEBND_GSL<T,PMT,PVT>::hausdorff(
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

template <typename T, typename PMT, typename PVT> template<typename ODESLV> inline typename ODEBND_GSL<T,PMT,PVT>::STATUS
ODEBND_GSL<T,PMT,PVT>::hausdorff
( const unsigned ns, const double*tk, const PVT*PMp, double**Hxk,
  ODESLV&traj, const unsigned nsamp, std::ostream&os )
{
  return _hausdorff( ns, tk, PMp, Hxk, *this, traj, nsamp, os )? NORMAL: FAILURE;
}

} // end namescape mc

#endif

