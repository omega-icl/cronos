// Copyright (C) 2015 Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

#ifndef MC__IODEBND_SUNDIALS_HPP
#define MC__IODEBND_SUNDIALS_HPP

#undef  MC__IODEBND_SUNDIALS_DINEQI_DEBUG
#undef  MC__IODEBND_SUNDIALS_DINEQPM_DEBUG

#include "iodebnd_base.hpp"
#include "odebnd_sundials.hpp"

namespace mc
{
//! @brief C++ class computing enclosures of the reachable set of parametric IODEs or DAEs using continuous-time set-valued integration and the numerical integrator SUNDIALS-CVODE.
////////////////////////////////////////////////////////////////////////
//! mc::IODEBND_SUNDIALS is a C++ class that computes enclosures of the
//! reachable set of parametric implicit ordinary differential equations
//! (IODEs) or differential-algebraic equations (DAEs) using
//! continuous-time set-valued integration. The method for DAEs
//! considers the underlying ODEs. Moreover, it implements the method
//! of differential inequalities, whereby polynomial models
//! with interval or ellipsoidal remainders are used to enable high-
//! order convergence. The use of ellipsoidal remainders enables
//! stability of the enclosures for asymptotically stable ODE systems
//! when the parameter host is sufficiently small. The numerical
//! integrator is CVODE in SUNDIALS.
////////////////////////////////////////////////////////////////////////
template <typename T, typename PMT=mc::TModel<T>, typename PVT=mc::TVar<T> >
class IODEBND_SUNDIALS:
  public virtual BASE_DE,
  public virtual BASE_SUNDIALS,
  public virtual IODEBND_BASE<T,PMT,PVT>,
  protected virtual ODEBND_SUNDIALS<T,PMT,PVT>
{
 typedef BASE_DE::STATUS STATUS;
 typedef int (*CVRhsFn)( realtype t, N_Vector y, N_Vector ydot, void *user_data );

 protected:
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
  using IODEBND_BASE<T,PMT,PVT>::_IC_I_STA;
  using ODEBND_BASE<T,PMT,PVT>::_IC_I_QUAD;
  using IODEBND_BASE<T,PMT,PVT>::_CC_I_STA;
  using ODEBND_BASE<T,PMT,PVT>::_SET_I_STA;
  using IODEBND_BASE<T,PMT,PVT>::_RHS_I_STA;
  using ODEBND_BASE<T,PMT,PVT>::_RHS_I_QUAD;
  using ODEBND_BASE<T,PMT,PVT>::_JAC_I_STA;
  using ODEBND_BASE<T,PMT,PVT>::_FCT_I_STA;
  using IODEBND_BASE<T,PMT,PVT>::_URHS_STA;

  using ODEBND_BASE<T,PMT,PVT>::_PMenv;
  using ODEBND_BASE<T,PMT,PVT>::_PMx;
  using ODEBND_BASE<T,PMT,PVT>::_PMp;
  using ODEBND_BASE<T,PMT,PVT>::_PMq;
  using ODEBND_BASE<T,PMT,PVT>::_vec2PMI;
  using ODEBND_BASE<T,PMT,PVT>::_vec2PME;
  using IODEBND_BASE<T,PMT,PVT>::_IC_PM_STA;
  using ODEBND_BASE<T,PMT,PVT>::_IC_PM_QUAD;
  using IODEBND_BASE<T,PMT,PVT>::_CC_PM_STA;
  using ODEBND_BASE<T,PMT,PVT>::_SET_PM_STA;
  using IODEBND_BASE<T,PMT,PVT>::_RHS_PM_STA;
  using ODEBND_BASE<T,PMT,PVT>::_RHS_PM_QUAD;
  using ODEBND_BASE<T,PMT,PVT>::_JAC_PM_STA;
  using ODEBND_BASE<T,PMT,PVT>::_FCT_PM_STA;

  using ODEBND_BASE<T,PMT,PVT>::_diam;
  using ODEBND_BASE<T,PMT,PVT>::_hausdorff;
  using ODEBND_BASE<T,PMT,PVT>::_remainders;
  using ODEBND_BASE<T,PMT,PVT>::_print_interm;

  using ODEBND_SUNDIALS<T,PMT,PVT>::_cv_mem;
  using ODEBND_SUNDIALS<T,PMT,PVT>::_cv_flag;
  using ODEBND_SUNDIALS<T,PMT,PVT>::_Nx;
  using ODEBND_SUNDIALS<T,PMT,PVT>::_Nq;
  using ODEBND_SUNDIALS<T,PMT,PVT>::_nchk;
  using ODEBND_SUNDIALS<T,PMT,PVT>::_vec_sta;
  using ODEBND_SUNDIALS<T,PMT,PVT>::_pos_rhs;
  using ODEBND_SUNDIALS<T,PMT,PVT>::_pos_quad;
  using ODEBND_SUNDIALS<T,PMT,PVT>::_pos_fct;

  using ODEBND_SUNDIALS<T,PMT,PVT>::_INI_CVODE;
  using ODEBND_SUNDIALS<T,PMT,PVT>::_CC_CVODE;
  using ODEBND_SUNDIALS<T,PMT,PVT>::_END_STA;

  //! @brief static pointer to class
  static IODEBND_SUNDIALS<T,PMT,PVT> *_pIODEBND;

 public:
  typedef typename ODEBND_BASE<T,PMT,PVT>::Results Results;
  typedef typename ODEBND_SUNDIALS<T,PMT,PVT>::Exceptions Exceptions;
  using ODEBND_SUNDIALS<T,PMT,PVT>::results_sta;
  using ODEBND_SUNDIALS<T,PMT,PVT>::stats_sta;

  /** @defgroup IODEBND_SUNDIALS Continuous-time set-valued integration of parametric IODEs or DAEs
   *  @{
   */
  //! @brief Default constructor
  IODEBND_SUNDIALS();

  //! @brief Virtual destructor
  virtual ~IODEBND_SUNDIALS();

  //! @brief Integrator options
  struct Options: public ODEBND_SUNDIALS<T,PMT,PVT>::Options
  {
    //! @brief Constructor
    Options():
      ODEBND_SUNDIALS<T,PMT,PVT>::Options(),
      ICBNDOPT(typename AEBND<T,PMT,PVT>::Options(0)),
      RHSBNDOPT(typename AEBND<T,PMT,PVT>::Options(0)), RHSNUMER(true)
      {}
    //! @brief Assignment operator
    Options& operator=
      ( Options&options ){
        ODEBND_SUNDIALS<T,PMT,PVT>::Options::operator=(options);
        ICBNDOPT   = options.AEBNDOPT;
        RHSBNDOPT  = options.RHSBNDOPT;
        RHSNUMER   = options.RHSNUMER;
        return *this;
      }
    //! @brief Options of AE bounder for initial conditions
    typename AEBND<T,PMT,PVT>::Options ICBNDOPT;
    //! @brief Options of AE bounder for right-hand side evaluation
    typename AEBND<T,PMT,PVT>::Options RHSBNDOPT;
    //! @brief Whether to compute the underlying ODE numerically (true) or symbolically (false)
    unsigned int RHSNUMER;
  } options;

  //! @brief Propagate state/quadrature interval bounds forward in time through every time stages
  STATUS bounds
    ( const unsigned ns, const double*tk, const T*Ip, T**Ixk,
      T*Iq=0, T*If=0, std::ostream&os=std::cout )
      { ODEBND_SUNDIALS<T,PMT,PVT>::options = options;
        return _bounds( ns, tk, Ip, Ixk, Iq, If, false, os); }

  //! @brief Propagate state/quadrature polynomial models forward in time through every time stages
  STATUS bounds
    ( const unsigned ns, const double*tk, const PVT*PMp, PVT**PMxk,
      PVT*PMq=0, PVT*PMf=0, std::ostream&os=std::cout )
      { ODEBND_SUNDIALS<T,PMT,PVT>::options = options;
        return _bounds( ns, tk, PMp, PMxk, PMq, PMf, false, os); }

  //! @brief Record state bounds in file <a>obndsta</a>, with accuracy of <a>iprec</a> digits
  void record
    ( std::ofstream&obndsta, const unsigned iprec=5 ) const
    { this->ODEBND_SUNDIALS<T,PMT,PVT>::record( obndsta, iprec ); }
  /** @} */

protected:
  //! @brief Function to initialize state interval bounding
  bool _INI_I_STA
    ( const unsigned np, const T*Ip, const unsigned ns );

  //! @brief Static wrapper to function computing the DINEQs RHS in interval arithmetic
  static int MC_CVRHSI__
    ( realtype t, N_Vector Nx, N_Vector Nxdot, void *user_data );

  //! @brief Static wrapper to function computing the quadratures RHS in interval arithmetic
  static int MC_CVQUADI__
    ( realtype t, N_Vector Nx, N_Vector Nqdot, void *user_data );

  //! @brief Propagate state/quadrature interval bounds forward in time through every time stages
  virtual STATUS _bounds
    ( const unsigned ns, const double*tk, const T*Ip, T**Ixk,
      T*Iq, T*If, const bool store, std::ostream&os );

  //! @brief Function to initialize SUNDIALS for state polynomial models
  bool _INI_PM_STA
    ( const unsigned np, const PVT*PMp, const unsigned ns );

  //! @brief Static wrapper to function computing the DINEQs RHS in polynomial model arithmetic
  static int MC_CVRHSPM__
    ( realtype t, N_Vector Nx, N_Vector Nxdot, void *user_data );

  //! @brief Static wrapper to function computing the quadratures RHS in polynomial model arithmetic
  static int MC_CVQUADPM__
    ( realtype t, N_Vector Nx, N_Vector Nqdot, void *user_data );

  //! @brief Propagate state/quadrature polynomial models forward in time through every time stages
  virtual STATUS _bounds
    ( const unsigned ns, const double*tk, const PVT*PMp, PVT**PMxk,
      PVT*PMq, PVT*PMf, const bool store, std::ostream&os );

  //! @brief Private methods to block default compiler methods
  IODEBND_SUNDIALS(const IODEBND_SUNDIALS&);
  IODEBND_SUNDIALS& operator=(const IODEBND_SUNDIALS&);
};

template <typename T, typename PMT, typename PVT>
 IODEBND_SUNDIALS<T,PMT,PVT>* IODEBND_SUNDIALS<T,PMT,PVT>::_pIODEBND = 0;

template <typename T, typename PMT, typename PVT> inline
IODEBND_SUNDIALS<T,PMT,PVT>::IODEBND_SUNDIALS
()
: BASE_DE(), BASE_SUNDIALS(), IODEBND_BASE<T,PMT,PVT>(), ODEBND_SUNDIALS<T,PMT,PVT>()
{}

template <typename T, typename PMT, typename PVT> inline
IODEBND_SUNDIALS<T,PMT,PVT>::~IODEBND_SUNDIALS
()
{}

template <typename T, typename PMT, typename PVT> inline bool
IODEBND_SUNDIALS<T,PMT,PVT>::_INI_I_STA
( const unsigned np, const T*Ip, const unsigned ns )
{
  // Initialize bound propagation
  if( !IODEBND_BASE<T,PMT,PVT>::_INI_I_STA( options, np, Ip, ns )
   || !ODEBND_SUNDIALS<T,PMT,PVT>::_INI_I_STA( np, ns ) )
    return false;
  return true;
}

template <typename T, typename PMT, typename PVT> inline int
IODEBND_SUNDIALS<T,PMT,PVT>::MC_CVRHSI__
( realtype t, N_Vector y, N_Vector ydot, void *user_data )
{
  IODEBND_SUNDIALS<T,PMT,PVT> *pIODEBND = IODEBND_SUNDIALS<T,PMT,PVT>::_pIODEBND;
  bool flag = pIODEBND->_RHS_I_STA( pIODEBND->options, t, NV_DATA_S( y ),
    NV_DATA_S( ydot ) );
  IODEBND_SUNDIALS<T,PMT,PVT>::_pIODEBND = pIODEBND;
  pIODEBND->stats_sta.numRHS++;
  return( flag? 0: -1 );
}

template <typename T, typename PMT, typename PVT> inline int
IODEBND_SUNDIALS<T,PMT,PVT>::MC_CVQUADI__
( realtype t, N_Vector y, N_Vector qdot, void *user_data )
{
  IODEBND_SUNDIALS<T,PMT,PVT> *pIODEBND = IODEBND_SUNDIALS<T,PMT,PVT>::_pIODEBND;
  bool flag = pIODEBND->_RHS_I_QUAD( pIODEBND->options, t, NV_DATA_S( y ),
    NV_DATA_S( qdot ) );
  IODEBND_SUNDIALS<T,PMT,PVT>::_pIODEBND = pIODEBND;
  return( flag? 0: -1 );
}

template <typename T, typename PMT, typename PVT>
inline typename IODEBND_SUNDIALS<T,PMT,PVT>::STATUS
IODEBND_SUNDIALS<T,PMT,PVT>::_bounds
( const unsigned ns, const double*tk, const T*Ip, T**Ixk, T*Iq, T*If,
  const bool store, std::ostream&os )
{
  // Check size
  if( !tk || !Ixk || !Ip || (_nf && !If) ) return FATAL;

  try{
    // Initialize trajectory integration
    if( !_INI_I_STA( _np, Ip, ns ) ) return FATAL;

    // Bounds on initial states/quadratures
    _t = tk[0];
    if( !_IC_I_STA( options, _t, NV_DATA_S( _Nx ) )
     || (_Nq && !_IC_I_QUAD( options, NV_DATA_S( _Nq )) ) )
      { _END_STA(); return FATAL; }
    if( options.DISPLAY >= 1 ){
      _print_interm( _t, _nx, _Ix, "x", os );
      _print_interm( _nq, _Iq, "q", os );
    }
    if( Ixk && !Ixk[0] ) Ixk[0] = new T[_nx];
    for( unsigned ix=0; Ixk[0] && ix<_nx; ix++ ) Ixk[0][ix] = _Ix[ix];

    // Store full state at initial time
    if( store ){
      realtype*vsta = NV_DATA_S(_Nx);
      unsigned lsta = NV_LENGTH_S(_Nx);
      _vec_sta.push_back( std::vector<realtype>( vsta, vsta+lsta ) );
    }

    // Record initial results
    if( options.RESRECORD )
      results_sta.push_back( Results( tk[0], _nx, Ixk[0] ) );

    // Integrate ODEs through each stage using SUNDIALS
    _pIODEBND = this;
    if( !_INI_CVODE( MC_CVRHSI__, MC_CVQUADI__ ) )
      { _END_STA(); return FATAL; }

    for( _istg=0; _istg<ns; _istg++ ){
      // Bounds on state discontinuities (if any) at stage times
      // and integrator reinitialization (if applicable)
      _pos_rhs = ( _vAE.size()<=1? 0: _istg );
      if( _pos_rhs
       && ( !_CC_I_STA( options, _pos_rhs, _t, NV_DATA_S( _Nx ) )
         || !_CC_CVODE() ) )
        { _END_STA(); return FAILURE; }

      // update list of operations in RHS, JAC and QUAD
      _pos_rhs  = ( _vRHS.size()<=1? 0: _istg );
      _pos_quad = ( _vQUAD.size()<=1? 0: _istg );
      if( (!_istg || _pos_rhs || _pos_quad)
       && (!_URHS_STA( options, _pos_rhs, _pos_quad )
        || !_SET_I_STA( options )) )
        { _END_STA(); return FATAL; }

      // integrate till end of time stage
      _cv_flag = CVodeSetStopTime( _cv_mem, tk[_istg+1] );
      if( _check_cv_flag(&_cv_flag, "CVodeSetStopTime", 1) )
        { _END_STA(); return FATAL; }
      while( _t < tk[_istg+1] ){
        if( !store )
          _cv_flag = CVode( _cv_mem, tk[_istg+1], _Nx, &_t, CV_ONE_STEP );
        else
          _cv_flag = CVodeF( _cv_mem, tk[_istg+1], _Nx, &_t, CV_ONE_STEP, &_nchk );
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

template <typename T, typename PMT, typename PVT> inline bool
IODEBND_SUNDIALS<T,PMT,PVT>::_INI_PM_STA
( const unsigned np, const PVT*PMp, const unsigned ns )
{
  // Initialize bound propagation
  if( !IODEBND_BASE<T,PMT,PVT>::_INI_PM_STA( options, np, PMp, ns )
   || !ODEBND_SUNDIALS<T,PMT,PVT>::_INI_PM_STA( np, ns ) )
    return false;
  return true;
}

template <typename T, typename PMT, typename PVT> inline int
IODEBND_SUNDIALS<T,PMT,PVT>::MC_CVRHSPM__
( realtype t, N_Vector y, N_Vector ydot, void *user_data )
{
  IODEBND_SUNDIALS<T,PMT,PVT> *pIODEBND = IODEBND_SUNDIALS<T,PMT,PVT>::_pIODEBND;
  bool flag = pIODEBND->_RHS_PM_STA( pIODEBND->options, t, NV_DATA_S( y ),
    NV_DATA_S( ydot ) );
  IODEBND_SUNDIALS<T,PMT,PVT>::_pIODEBND = pIODEBND;
  pIODEBND->stats_sta.numRHS++;
  return( flag? 0: -1 );
}

template <typename T, typename PMT, typename PVT> inline int
IODEBND_SUNDIALS<T,PMT,PVT>::MC_CVQUADPM__
( realtype t, N_Vector y, N_Vector qdot, void *user_data )
{
  IODEBND_SUNDIALS<T,PMT,PVT> *pIODEBND = IODEBND_SUNDIALS<T,PMT,PVT>::_pIODEBND;
  bool flag = pIODEBND->_RHS_PM_QUAD( pIODEBND->options, t, NV_DATA_S( y ),
    NV_DATA_S( qdot ) );
  IODEBND_SUNDIALS<T,PMT,PVT>::_pIODEBND = pIODEBND;
  return( flag? 0: -1 );
}

template <typename T, typename PMT, typename PVT>
inline typename IODEBND_SUNDIALS<T,PMT,PVT>::STATUS
IODEBND_SUNDIALS<T,PMT,PVT>::_bounds
( const unsigned ns, const double*tk, const PVT*PMp, PVT**PMxk,
  PVT*PMq, PVT*PMf, const bool store, std::ostream&os )
{
  // Check arguments
  if( !tk || !PMxk || !PMp || (_nf && !PMf) ) return FATAL;

  try{
    // Initialize trajectory integration
    if( !_INI_PM_STA( _np, PMp, ns ) ) return FATAL;

    // Bounds on initial states/quadratures
    _t = tk[0];
    if( !_IC_PM_STA( options, _t, NV_DATA_S( _Nx ) )
     || (_Nq && !_IC_PM_QUAD( options, NV_DATA_S( _Nq )) ) )
      { _END_STA(); return FATAL; }
    if( options.DISPLAY >= 1 ){
      _print_interm( _t, _nx, _PMx, "x", os );
      _print_interm( _nq, _PMq, "q", os );
    }
    if( PMxk && !PMxk[0] ) PMxk[0] = new PVT[_nx];
    for( unsigned ix=0; PMxk[0] && ix<_nx; ix++ ) PMxk[0][ix] = _PMx[ix];

    // Store full state at initial time
    if( store ){
      realtype*vsta = NV_DATA_S(_Nx);
      unsigned lsta = NV_LENGTH_S(_Nx);
      _vec_sta.push_back( std::vector<realtype>( vsta, vsta+lsta ) );
    }

    // Record initial results
    if( options.RESRECORD )
      results_sta.push_back( Results( tk[0], _nx, PMxk[0] ) );

    // Integrate ODEs through each stage using SUNDIALS
    _pIODEBND = this;
    if( !_INI_CVODE( MC_CVRHSPM__, MC_CVQUADPM__ ) )
      { _END_STA(); return FATAL; }

    for( _istg=0; _istg<ns; _istg++ ){
      // Bounds on state discontinuities (if any) at stage times
      // and integrator reinitialization (if applicable)
      _pos_rhs = ( _vAE.size()<=1? 0: _istg );
      if( _pos_rhs
       && ( !_CC_PM_STA( options, _pos_rhs, _t, NV_DATA_S( _Nx ) )
         || !_CC_CVODE() ) )
        { _END_STA(); return FAILURE; }

      // update list of operations in RHS, JAC and QUAD
      _pos_rhs  = ( _vRHS.size()<=1? 0: _istg );
      _pos_quad = ( _vQUAD.size()<=1? 0: _istg );
      if( (!_istg || _pos_rhs || _pos_quad)
       && (!_URHS_STA( options, _pos_rhs, _pos_quad )
        || !_SET_PM_STA( options )) )
        { _END_STA(); return FATAL; }

      // integrate till end of time stage
      _cv_flag = CVodeSetStopTime( _cv_mem, tk[_istg+1] );
      if( _check_cv_flag(&_cv_flag, "CVodeSetStopTime", 1) )
        { _END_STA(); return FATAL; }
      while( _t < tk[_istg+1] ){
        if( !store )
          _cv_flag = CVode( _cv_mem, tk[_istg+1], _Nx, &_t, CV_ONE_STEP );
        else
          _cv_flag = CVodeF( _cv_mem, tk[_istg+1], _Nx, &_t, CV_ONE_STEP, &_nchk );
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

      // Bounds on intermediate states
      switch( options.WRAPMIT){
      case Options::NONE:
        _vec2PMI( NV_DATA_S( _Nx ), _PMenv, _nx, _PMx, true );
        break;
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

} // end namescape mc

#endif

