// Copyright (C) 2015 Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

#ifndef MC__IODEBND_GSL_HPP
#define MC__IODEBND_GSL_HPP

#undef  MC__IODEBND_GSL_DINEQI_DEBUG
#undef  MC__IODEBND_GSL_DINEQPM_DEBUG

#include "iodebnd_base.hpp"
#include "odebnd_gsl.hpp"

namespace mc
{
//! @brief C++ class computing enclosures of the reachable set of parametric IODEs or DAEs using continuous-time set-valued integration and the numerical integrator in GSL.
////////////////////////////////////////////////////////////////////////
//! mc::IODEBND_GSL is a C++ class that computes enclosures of the
//! reachable set of parametric implicit ordinary differential equations
//! (IODEs) or differential-algebraic equations (DAEs) using
//! continuous-time set-valued integration. The method for DAEs
//! considers the underlying ODEs. Moreover, it implements the method
//! of differential inequalities, whereby polynomial models
//! with interval or ellipsoidal remainders are used to enable high-
//! order convergence. The use of ellipsoidal remainders enables
//! stability of the enclosures for asymptotically stable ODE systems
//! when the parameter host is sufficiently small. The numerical
//! integrator is gsl_odeiv2 in GSL.
////////////////////////////////////////////////////////////////////////
template <typename T, typename PMT=mc::TModel<T>, typename PVT=mc::TVar<T> >
class IODEBND_GSL:
  public virtual BASE_DE,
  public virtual BASE_GSL,
  public virtual IODEBND_BASE<T,PMT,PVT>,
  protected virtual ODEBND_GSL<T,PMT,PVT>
{
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

  using ODEBND_GSL<T,PMT,PVT>::_sys_sta;
  using ODEBND_GSL<T,PMT,PVT>::_driver_sta;
  using ODEBND_GSL<T,PMT,PVT>::_vec_sta;
  using ODEBND_GSL<T,PMT,PVT>::_vec_quad;
  using ODEBND_GSL<T,PMT,PVT>::_offset_quad;
  using ODEBND_GSL<T,PMT,PVT>::_h;
  using ODEBND_GSL<T,PMT,PVT>::_mesh_sta;
  using ODEBND_GSL<T,PMT,PVT>::_pos_ic;
  using ODEBND_GSL<T,PMT,PVT>::_pos_rhs;
  using ODEBND_GSL<T,PMT,PVT>::_pos_quad;
  using ODEBND_GSL<T,PMT,PVT>::_pos_fct;
  using ODEBND_GSL<T,PMT,PVT>::PMOFFSET;

  using ODEBND_GSL<T,PMT,PVT>::_INI_GSL;
  using ODEBND_GSL<T,PMT,PVT>::_END_STA;

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

  //! @brief static pointer to class
  static IODEBND_GSL<T,PMT,PVT> *_pIODEBND;

 public:
  typedef typename ODEBND_BASE<T,PMT,PVT>::Results Results;
  typedef typename ODEBND_GSL<T,PMT,PVT>::Exceptions Exceptions;
  using ODEBND_GSL<T,PMT,PVT>::results_sta;
  using ODEBND_GSL<T,PMT,PVT>::stats_sta;

  /** @defgroup IODEBND_GSL Continuous-time set-valued integration of parametric IODEs or DAEs
   *  @{
   */
  //! @brief Default constructor
  IODEBND_GSL();

  //! @brief Virtual destructor
  virtual ~IODEBND_GSL();

  //! @brief Integrator options
  struct Options: public ODEBND_GSL<T,PMT,PVT>::Options
  {
    //! @brief Constructor
    Options():
      ODEBND_GSL<T,PMT,PVT>::Options(),
      ICBNDOPT(typename AEBND<T,PMT,PVT>::Options(0)),
      RHSBNDOPT(typename AEBND<T,PMT,PVT>::Options(0)), RHSNUMER(true)
      {}
    //! @brief Assignment operator
    Options& operator=
      ( Options&options ){
        ODEBND_GSL<T,PMT,PVT>::Options::operator=(options);
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
      { ODEBND_GSL<T,PMT,PVT>::options = options;
        return _bounds( ns, tk, Ip, Ixk, Iq, If, false, os); }

  //! @brief Propagate state/quadrature polynomial models forward in time through every time stages
  STATUS bounds
    ( const unsigned ns, const double*tk, const PVT*PMp, PVT**PMxk,
      PVT*PMq=0, PVT*PMf=0, std::ostream&os=std::cout )
      { ODEBND_GSL<T,PMT,PVT>::options = options;
        return _bounds( ns, tk, PMp, PMxk, PMq, PMf, false, os); }

  //! @brief Record state bounds in file <a>obndsta</a>, with accuracy of <a>iprec</a> digits
  void record
    ( std::ofstream&obndsta, const unsigned iprec=5 ) const
    { this->ODEBND_GSL<T,PMT,PVT>::record( obndsta, iprec ); }
  /** @} */

 protected:
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
  virtual STATUS _bounds
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
  virtual STATUS _bounds
    ( const unsigned ns, const double*tk, const PVT*PMp, PVT**PMxk,
      PVT*PMq, PVT*PMf, const bool store, std::ostream&os );

  //! @brief Private methods to block default compiler methods
  IODEBND_GSL(const IODEBND_GSL&);
  IODEBND_GSL& operator=(const IODEBND_GSL&);
};

template <typename T, typename PMT, typename PVT>
 IODEBND_GSL<T,PMT,PVT>* IODEBND_GSL<T,PMT,PVT>::_pIODEBND = 0;

template <typename T, typename PMT, typename PVT> inline
IODEBND_GSL<T,PMT,PVT>::IODEBND_GSL
()
: BASE_DE(), BASE_GSL(), IODEBND_BASE<T,PMT,PVT>(), ODEBND_GSL<T,PMT,PVT>()
{}

template <typename T, typename PMT, typename PVT> inline
IODEBND_GSL<T,PMT,PVT>::~IODEBND_GSL
()
{}

template <typename T, typename PMT, typename PVT> inline bool
IODEBND_GSL<T,PMT,PVT>::_INI_I_STA
( const unsigned np, const T*Ip, const unsigned ns )
{
  // Initialize bound propagation
  ODEBND_GSL<T,PMT,PVT>::_sys_sta.function = MC_GSLRHSI__;
  ODEBND_GSL<T,PMT,PVT>::_sys_sta.jacobian = MC_GSLJACI__;
  if( !IODEBND_BASE<T,PMT,PVT>::_INI_I_STA( options, np, Ip, ns )
   || !ODEBND_GSL<T,PMT,PVT>::_INI_I_STA( np, ns ) )
    return false;
  return true;
}

template <typename T, typename PMT, typename PVT> inline int
IODEBND_GSL<T,PMT,PVT>::MC_GSLRHSI__
( double t, const double* x, double* xdot, void* user_data )
{
  IODEBND_GSL<T,PMT,PVT> *pIODEBND = IODEBND_GSL<T,PMT,PVT>::_pIODEBND;
  bool flag = pIODEBND->_RHS_I_STA( pIODEBND->options, t, x, xdot );
  if( flag && pIODEBND->_nq ){
    double* qdot = xdot + pIODEBND->_offset_quad;
    flag = pIODEBND->_RHS_I_QUAD( pIODEBND->options, t, x, qdot );
  }
  pIODEBND->stats_sta.numRHS++;
  IODEBND_GSL<T,PMT,PVT>::_pIODEBND = pIODEBND;
  return( flag? GSL_SUCCESS: GSL_EBADFUNC );
}

template <typename T, typename PMT, typename PVT> inline int
IODEBND_GSL<T,PMT,PVT>::MC_GSLJACI__
( double t, const double* x, double* jac, double* xdot, void* user_data )
{
  IODEBND_GSL<T,PMT,PVT> *pIODEBND = IODEBND_GSL<T,PMT,PVT>::_pIODEBND;
  bool flag = pIODEBND->_JAC_I_STA( pIODEBND->options, t, x, jac, xdot );
  pIODEBND->stats_sta.numJAC++;
  IODEBND_GSL<T,PMT,PVT>::_pIODEBND = pIODEBND;
  return( flag? GSL_SUCCESS: GSL_EBADFUNC );
}

template <typename T, typename PMT, typename PVT>
inline typename IODEBND_GSL<T,PMT,PVT>::STATUS
IODEBND_GSL<T,PMT,PVT>::_bounds
( const unsigned ns, const double*tk, const T*Ip, T**Ixk, T*Iq, T*If,
  const bool store, std::ostream&os )
{
  // Check size
  if( !tk || !Ixk || !Ip || (_nf && !If) ) return FATAL;

  try{
    // Initialize trajectory integration with GSL
    _INI_I_STA( _np, Ip, ns );

    // Bounds on initial states/quadratures
    _t = tk[0];
    if( !_IC_I_STA( options, _t, _vec_sta )
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
    _pIODEBND = this;

    for( _istg=0; _istg<ns; _istg++ ){
      // Bounds on state discontinuities (if any) at stage times
      // and integrator reinitialization (if applicable)
      _pos_rhs = ( _vAE.size()<=1? 0: _istg );
      if( _pos_rhs ){
        if( !_CC_I_STA( options, _pos_rhs, _t, _vec_sta ) )
          { _END_STA(); return FATAL; }
        gsl_odeiv2_driver_reset( _driver_sta );
      }

      // update list of operations in RHS, JAC and QUAD
      _pos_rhs  = ( _vRHS.size()<=1? 0: _istg );
      _pos_quad = ( _vQUAD.size()<=1? 0: _istg );
      if( (!_istg || _pos_rhs || _pos_quad)
       && (!_URHS_STA( options, _pos_rhs, _pos_quad )
        || !_SET_I_STA( options )) )
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

template <typename T, typename PMT, typename PVT> inline bool
IODEBND_GSL<T,PMT,PVT>::_INI_PM_STA
( const unsigned np, const PVT*PMp, const unsigned ns )
{
  // Initialize bound propagation
  ODEBND_GSL<T,PMT,PVT>::_sys_sta.function = MC_GSLRHSPM__;
  ODEBND_GSL<T,PMT,PVT>::_sys_sta.jacobian = MC_GSLJACPM__;
  if( !IODEBND_BASE<T,PMT,PVT>::_INI_PM_STA( options, np, PMp, ns )
   || !ODEBND_GSL<T,PMT,PVT>::_INI_PM_STA( np, ns ) )
    return false;
  return true;
}

template <typename T, typename PMT, typename PVT> inline int
IODEBND_GSL<T,PMT,PVT>::MC_GSLRHSPM__
( double t, const double* x, double* xdot, void* user_data )
{
  IODEBND_GSL<T,PMT,PVT> *pIODEBND = IODEBND_GSL<T,PMT,PVT>::_pIODEBND;
  bool flag = pIODEBND->_RHS_PM_STA( pIODEBND->options, t, x, xdot );
  if( flag && pIODEBND->_nq ){
    double* qdot = xdot + pIODEBND->_offset_quad;
    flag = pIODEBND->_RHS_PM_QUAD( pIODEBND->options, t, x, qdot );
  }
  pIODEBND->stats_sta.numRHS++;
  IODEBND_GSL<T,PMT,PVT>::_pIODEBND = pIODEBND;
  return( flag? GSL_SUCCESS: GSL_EBADFUNC );
}

template <typename T, typename PMT, typename PVT> inline int
IODEBND_GSL<T,PMT,PVT>::MC_GSLJACPM__
( double t, const double* x, double* jac, double* xdot, void* user_data )
{
  IODEBND_GSL<T,PMT,PVT> *pIODEBND = IODEBND_GSL<T,PMT,PVT>::_pIODEBND;
  bool flag = pIODEBND->_JAC_PM_STA( pIODEBND->options, t, x, jac, xdot );
  pIODEBND->stats_sta.numJAC++;
  IODEBND_GSL<T,PMT,PVT>::_pIODEBND = pIODEBND;
  return( flag? GSL_SUCCESS: GSL_EBADFUNC );
}

template <typename T, typename PMT, typename PVT>
inline typename IODEBND_GSL<T,PMT,PVT>::STATUS
IODEBND_GSL<T,PMT,PVT>::_bounds
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
    if( !_IC_PM_STA( options, _t, _vec_sta )
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
    _pIODEBND = this;

    for( _istg=0; _istg<ns; _istg++ ){
      // Bounds on state discontinuities (if any) at stage times
      // and integrator reinitialization (if applicable)
      _pos_rhs = ( _vAE.size()<=1? 0: _istg );
      if( _pos_rhs ){
        if( !_CC_PM_STA( options, _pos_rhs, _t, _vec_sta ) )
          { _END_STA(); return FATAL; }
        gsl_odeiv2_driver_reset( _driver_sta );
      }

      // update list of operations in RHS, JAC and QUAD
      _pos_rhs = ( _vRHS.size()<=1? 0: _istg );
      _pos_quad = ( _vQUAD.size()<=1? 0: _istg );
      if( (!_istg || _pos_rhs || _pos_quad)
       && (!_URHS_STA( options, _pos_rhs, _pos_quad )
        || !_SET_PM_STA( options )) )
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

} // end namescape mc

#endif

