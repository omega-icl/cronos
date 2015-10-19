// Copyright (C) 2012-2014 Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

#ifndef MC__ODESLV_GSL_HPP
#define MC__ODESLV_GSL_HPP

#undef  MC__ODESLV_GSL_SAMPLE_DEBUG
#undef  MC__ODESLV_GSL_DEBUG_INTERP
#undef  MC__ODESLV_GSL_DEBUG_ADJOINT
#define MC__ODESLV_GSL_CHECK

#include <stdexcept>
#include <cassert>
#include <cmath>
#include <fstream>
#include <vector>
#include <sys/time.h>

#include "base_de.hpp"
#include "base_gsl.hpp"
#include "mesh_gsl.hpp"

namespace mc
{
//! @brief C++ class computing solutions of parametric ODEs using non-validated integration.
////////////////////////////////////////////////////////////////////////
//! mc::ODESLV_GSL is a C++ class that computes solutions of parametric
//! ordinary differential equations (ODEs) using GSL.
////////////////////////////////////////////////////////////////////////
template <typename T>
class ODESLV_GSL: public virtual BASE_DE, public virtual BASE_GSL
{
private:
  //! @brief GSL data type for ODE integration
  gsl_odeiv2_system _sys_traj;

  //! @brief GSL driver for ODE integration
  gsl_odeiv2_driver *_driver_traj;

protected:
  //! @brief full GSL state
  double *_vec_sta;

  //! @brief stepsize
  double _h;

  //! @brief number of variables for DAG evaluation
  unsigned _nVAR;

  //! @brief pointer to variables for DAG evaluation
  FFVar* _pVAR;

  //! @brief pointer to variable values for DAG evaluation
  double *_dVAR;

  //! @brief list of operations in RHS evaluation
  std::list<const FFOp*> _opRHS;

  //! @brief list of operations in RHS Jacobian
  std::list<const FFOp*> _opJAC;

  //! @brief list of operations in IC evalution
  std::list<const FFOp*> _opIC;

  //! @brief pointer to RHS function in current stage of ODE system
  FFVar* _pRHS;

  //! @brief preallocated array for evaluation of RHS function
  double* _dRHS;

  //! @brief const pointer to RHS Jacobian in current stage of ODE system
  const FFVar* _pJAC;

  //! @brief preallocated array for evaluation of RHS Jacobian
  double* _dJAC;

  //! @brief const pointer to IC function in current stage of ODE system
  const FFVar* _pIC;

  //! @brief Mesh storing state bound parameterizations 
  MESH_GSL _mesh_sta;

public:
  /** @defgroup ODESLV_GSL Real-valued (non-validated) integration of parametric ODEs
   *  @{
   */
  //! @brief Default class constructor
  ODESLV_GSL()
    : BASE_GSL(), _driver_traj(0), _vec_sta(0), _nVAR(0), _pVAR(0), _dVAR(0),
      _pRHS(0), _dRHS(0), _pJAC(0), _dJAC(0), _pIC(0)
    {}

  //! @brief Default destructor
  virtual ~ODESLV_GSL()
    {    
      if( _driver_traj ) gsl_odeiv2_driver_free( _driver_traj );
      delete[] _vec_sta;
      delete[] _dVAR;
      delete[] _dRHS;
      delete[] _dJAC;
      delete[] _pVAR;
      /* DO NOT FREE _pIC */
      delete[] _pRHS;
      delete[] _pJAC;
    }

  //! @brief Integrator options
  struct Options: public BASE_GSL::Options
  {
    //! @brief Constructor
    Options():
      BASE_GSL::Options(), MESHPREALLOC(0), DISPLAY(0), RESRECORD(false)
      {}
    //! @brief Assignment operator
    template <typename U> Options& operator=
      ( U&options ){
        BASE_GSL::Options::operator=(options);
        MESHPREALLOC = options.MESHPREALLOC;
        DISPLAY      = options.DISPLAY;
        RESRECORD    = options.RESRECORD;
        return *this;
      }
    //! @brief Preallocated mesh size (default: 0)
    bool MESHPREALLOC;
    //! @brief Display level
    int DISPLAY;
    //! @brief Whether or not to record results (default: false)
    bool RESRECORD;
  } options;

  //! @brief Structure for setting up storing the solver exceptions
  class Exceptions
  {
  public:
    //! @brief Enumeration type for NLPSBB exception handling
    enum TYPE{
      UNDEF=-33	//!< Error due to calling a function/feature not yet implemented
    };
    //! @brief Constructor for error <a>ierr</a>
    Exceptions( TYPE ierr ) : _ierr( ierr ){}
    //! @brief Inline function returning the error flag
    int ierr(){ return _ierr; }
    //! @brief Inline function returning the error description
    std::string what(){
      switch( _ierr ){
      case UNDEF: default:
        return "ODESLV_GSL::Exceptions  Error due to calling a not yet implemented feature";
      }
    }
  private:
    TYPE _ierr;
  };

  //! @brief Store integration bounds at a given time instant
  struct Results
  {
    //! @brief Constructors
    Results
      ( const double tk, const unsigned nxk, const T*Ixk ):
      t( tk ), nx( nxk )
      { X = new T[nx];
        for( unsigned ix=0; ix<nx; ix++ ) X[ix] = Ixk[ix]; }
    Results
      ( const Results&res ):
      t( res.t ), nx( res.nx )
      { X = new T[nx];
        for( unsigned ix=0; ix<nx; ix++ ) X[ix] = res.X[ix]; }
    //! @brief Destructor
    virtual ~Results()
      { delete[] X; }
    //! @brief Time instant
    double t;
    //! @brief Solution dimension
    unsigned nx;
    //! @brief Solution bounds
    T* X;
  };

  //! @brief Statistics for state integration
  Stats stats_sta;

  //! @brief Integrate trajectory of parametric ODEs
  STATUS states
    ( const unsigned ns, const double*tk, const double*p, double**xk,
      double*q, double*f, std::ostream&os=std::cout );

  //! @brief Compute approximate interval enclosure of reachable set of parametric ODEs using parameter sampling
  STATUS bounds
    ( const unsigned ns, const double*tk, const T*Ip, T**Ixk,
      T*Iq, T*If, const unsigned nsamp, std::ostream&os=std::cout );

  //! @brief Record state bounds in file <a>bndsta</a>, with accuracy of <a>iprec</a> digits
  void record
    ( std::ofstream&bndsta, const unsigned iprec=5 ) const;

  //! @brief static pointer to class
  static ODESLV_GSL<T> *pODESLV_GSL;
  /** @} */

protected:
  //! @brief Vector storing interval state bounds (see Options::RESRECORD)
  std::vector< Results > _results_sta;

  //! @brief Function to finalize statistics for GSL
  void _END_STA
    ( Stats &stats );

  //! @brief Function to initialize GSL numerical integration driver
  void _INI_GSL
    ( gsl_odeiv2_system &sys, gsl_odeiv2_driver *&driver );

  //! @brief Initialize GSL for state integration
  bool _INI_STA
    ( const double* p, const unsigned ns );

  //! @brief Evaluate initial conditions
  bool _IC_STA
   ();

  //! @brief Evaluate continuity conditions
  bool _CC_STA
    ();

  //! @brief Evaluate the functions at intermediate/end point
  bool _FCT_STA
    ( const unsigned iFCT, double*f );

  //! @brief Set RHS pointer and corresponding Jacobian
  bool _RHS_STA_SET
    ( const unsigned iRHS );

  //! @brief Integrate state trajectories on current stages
  STATUS _states_traj
    ( const double tf, const bool store );

  //! @brief Integrate state trajectories on every time stages
  STATUS _states
    ( const unsigned ns, const double*tk, const double*p, double**xk,
      double*q, double*f, const bool store, std::ostream&os );

  //! @brief Static wrapper to function to calculate the ODEs RHS values
  static int MC_GSLRHS__
    ( double t, const double* x, double* xdot, void* user_data );

  //! @brief Function to calculate the ODEs RHS values
  int _RHS_STA
    ( double t, const double* x, double* xdot, void* user_data );

  //! @brief Static wrapper to function to calculate the ODEs RHS derivatives
  static int MC_GSLJAC__
    ( double t, const double* x, double* jac, double* xdot, void* user_data );

  //! @brief Function to calculate the ODEs RHS derivatives
  int _JAC_STA
    ( double t, const double* x, double* jac, double* xdot, void* user_data );

  //! @brief Recursive function computing bounds on solutions of IVP in ODEs using sampling
  STATUS _states
    ( const unsigned ns, const double*tk, const T*Ip, T**Ixk, T*Iq, T*If,
      const unsigned nsamp, unsigned* vsamp, const unsigned ip,
      double*p, double**xk, double*q, double*f, std::ostream&os );

  //! @brief Function to display intermediate results
  template<typename U> static void _print_interm
    ( const double t, const unsigned nx, const U*x, const std::string&var,
      std::ostream&os=std::cout );

  //! @brief Function to display intermediate results
  template<typename U> static void _print_interm
    ( const unsigned nx, const U*x, const std::string&var,
      std::ostream&os=std::cout );

private:
  //! @brief Private methods to block default compiler methods
  ODESLV_GSL(const ODESLV_GSL&);
  ODESLV_GSL& operator=(const ODESLV_GSL&);
};

template <typename T>
ODESLV_GSL<T>* ODESLV_GSL<T>::pODESLV_GSL = 0;

template <typename T> inline bool
ODESLV_GSL<T>::_INI_STA
( const double*p, const unsigned ns )
{
  // Define ODE system in GSL format
  _sys_traj.function = MC_GSLRHS__;
  _sys_traj.jacobian = MC_GSLJAC__;
  _sys_traj.dimension = _nx+_nq;
  _sys_traj.params = 0;
  
  // Set GSL drivers for ODE integration
  _INI_GSL( _sys_traj, _driver_traj );

  // Initialize mesh
  if( !_mesh_sta.set( ns, _nx, options.MESHPREALLOC ) )
    return false;

  // Size and set DAG evaluation arrays
  if( _nVAR != _nx+_np+1+_nq ){
    _nVAR = _nx+_np+1+_nq;
    delete[] _pVAR; _pVAR = new FFVar[_nVAR];
    delete[] _dVAR; _dVAR = new double[_nVAR];
  }
  for( unsigned ix=0; ix<_nx; ix++ ) _pVAR[ix] = _pX[ix];
  for( unsigned ip=0; ip<_np; ip++ ) _pVAR[_nx+ip] = _pP[ip];
  _pVAR[_nx+_np] = (_pT? *_pT: 0. );
  for( unsigned iq=0; iq<_nq; iq++ ) _pVAR[_nx+_np+1+iq] = _pQ?_pQ[iq]:0.;
  for( unsigned ip=0; ip<_np; ip++ ) _dVAR[_nx+ip] = p[ip];

  // Initialize statistics
  _init_stats( stats_sta );

  return true;
}

template <typename T> inline void
ODESLV_GSL<T>::_INI_GSL
( gsl_odeiv2_system &sys, gsl_odeiv2_driver *&driver )
{
  // Set GSL numerical integration driver
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

  // Set state storage
  delete [] _vec_sta;
  _vec_sta = new double[ sys.dimension ];

  return;
}

template <typename T> inline void
ODESLV_GSL<T>::_END_STA
( Stats& stats )
{
  // Get final CPU time
  _final_stats( stats );
}

template <typename T> inline int
ODESLV_GSL<T>::_RHS_STA
( double t, const double* x, double* xdot, void* user_data )
{
  stats_sta.numRHS++; // increment RHS counter
  if( !_pRHS ) return GSL_EBADFUNC;
  for( unsigned i=0; i<_nx; i++ ) _dVAR[i] = x[i]; // set current state values
  _dVAR[_nx+_np] = t; // set current time
  _pDAG->eval( _opRHS, _dRHS, _nx+_nq, _pRHS, xdot, _nVAR-_nq, _pVAR, _dVAR );
  return GSL_SUCCESS;  
}

template <typename T> inline int
ODESLV_GSL<T>::MC_GSLRHS__
( double t, const double* x, double* xdot, void* user_data )
{
  ODESLV_GSL<T> *pODESLV_GSL = ODESLV_GSL<T>::pODESLV_GSL;
  int flag = pODESLV_GSL->_RHS_STA( t, x, xdot, user_data );
  ODESLV_GSL<T>::pODESLV_GSL = pODESLV_GSL;
  return flag;
}

template <typename T> inline int
ODESLV_GSL<T>::_JAC_STA
( double t, const double* x, double* jac, double* xdot, void* user_data )
{
  stats_sta.numJAC++; // increment JAC counter
  if( !_pRHS || !_pJAC ) return GSL_EBADFUNC;
  for( unsigned i=0; i<_nx; i++ ) _dVAR[i] = x[i]; // set current state values
  _dVAR[_nx+_np] = t; // set current time
  _pDAG->eval( _opRHS, _dRHS, _nx+_nq, _pRHS, xdot, _nVAR-_nq, _pVAR, _dVAR );
  _pDAG->eval( _opJAC, _dJAC, (_nx+_nq)*(_nx+_nq), _pJAC, jac, _nVAR-_nq, _pVAR, _dVAR );
  return GSL_SUCCESS;
}

template <typename T> inline int
ODESLV_GSL<T>::MC_GSLJAC__
( double t, const double* x, double* jac, double* xdot, void* user_data )
{
  ODESLV_GSL<T> *pODESLV_GSL = ODESLV_GSL<T>::pODESLV_GSL;
  int flag = pODESLV_GSL->_JAC_STA( t, x, jac, xdot, user_data );
  ODESLV_GSL<T>::pODESLV_GSL = pODESLV_GSL;
  return flag;
}

template <typename T> inline bool
ODESLV_GSL<T>::_IC_STA
()
{
  if( !_vIC.size() || _nx0 != _nx ) return false;
  try{
    _pIC = _vIC.at(0);
    _opIC = _pDAG->subgraph( _nx, _pIC );
    _dVAR[_nx+_np] = _t; // set initial time
    _pDAG->eval( _opIC, _nx, _pIC, _vec_sta, _np+1, _pVAR+_nx, _dVAR+_nx );
    for( unsigned iq=0; iq<_nq; iq++ ) _vec_sta[_nx+iq] = 0.;
  }
  catch(...){
    return false;
  }
  return true;
}

template <typename T> inline bool
ODESLV_GSL<T>::_CC_STA
()
{
  if( !_istg || _vIC.size() == 1 ) return true;
  if( _nx0 != _nx ) return false;
  try{
    _pIC = _vIC.at(_istg);
    _opIC = _pDAG->subgraph( _nx, _pIC );
    for( unsigned i=0; i<_nx; i++ ) _dVAR[i] = _vec_sta[i]; // set current state values
    _dVAR[_nx+_np] = _t; // set current time
    _pDAG->eval( _opIC, _nx, _pIC, _vec_sta, _nVAR-_nq, _pVAR, _dVAR );
    // Reset ODE solver - needed to account for discontinuity
    gsl_odeiv2_driver_reset( _driver_traj );
  }
  catch(...){
    return false;
  }
  return true;
}

template <typename T> inline bool
ODESLV_GSL<T>::_FCT_STA
( const unsigned iFCT, double*f )
{
  if( !_nf || !f ) return true;
  try{
    for( unsigned i=0; i<_nx; i++ ) _dVAR[i] = _vec_sta[i]; // set current states
    _dVAR[_nx+_np] = _t; // set current time
    for( unsigned i=0; i<_nq; i++ ) _dVAR[_nx+_np+1+i] = _vec_sta[_nx+i]; // set current quadratures
    const FFVar* pFCT = _vFCT.at( iFCT );
    _pDAG->eval( _nf, pFCT, f, _nVAR, _pVAR, _dVAR, iFCT?true:false );
  }
  catch(...){
    return false;
  }
  return true;
}

template <typename T> inline bool
ODESLV_GSL<T>::_RHS_STA_SET
( const unsigned iRHS )
{
  if( _vRHS.size() <= iRHS ) return false;
  if( _nq && _vQUAD.size() <= iRHS ) return false;

  delete[] _pRHS; _pRHS = new FFVar[_nx+_nq];
  for( unsigned i=0; i<_nx; i++ ) _pRHS[i] = _vRHS.at( iRHS )[i];
  for( unsigned i=0; i<_nq; i++ ) _pRHS[_nx+i] = _vQUAD.at( iRHS )[i];
  _opRHS = _pDAG->subgraph( _nx+_nq, _pRHS );
  delete[] _dRHS; _dRHS = new double[ _opRHS.size() ];

  switch( options.INTMETH ){
  case Options::MSADAMS:
  case Options::MSBDF:{
    FFVar* pXQ = new FFVar[_nx+_nq];
    for( unsigned i=0; i<_nx; i++ ) pXQ[i] = _pX[i];
    for( unsigned i=0; i<_nq; i++ ) pXQ[_nx+i] = _pQ?_pQ[i]:0.;
    delete[] _pJAC; _pJAC  = _pDAG->FAD( _nx+_nq, _pRHS, _nx+_nq, pXQ );
    delete[] pXQ;
    _opJAC = _pDAG->subgraph( (_nx+_nq)*(_nx+_nq), _pJAC );
    delete[] _dJAC; _dJAC  = new double[ _opJAC.size() ];
    break;
   }
  case Options::RK8PD:
  case Options::RKF45:
  default:
    delete[] _pJAC; _pJAC = 0;
    _opJAC.clear();
    delete[] _dJAC; _dJAC = 0;
    break;
  }

  return true;
}

template <typename T> inline typename ODESLV_GSL<T>::STATUS
ODESLV_GSL<T>::_states_traj
( const double tf, const bool store )
{
  // Store initial point
  if( store ) _mesh_sta.add( _istg, _t, _vec_sta, true );

  // integrate till end of time stage
  while( _t < tf ){
    if( gsl_odeiv2_evolve_apply( _driver_traj->e, _driver_traj->c, _driver_traj->s,
        &_sys_traj, &_t, tf, &_h, _vec_sta ) != GSL_SUCCESS
     || _h < options.HMIN
     || (options.NMAX && stats_sta.numSteps > options.NMAX) ) return FAILURE;
    stats_sta.numSteps++;
    if( options.HMAX > 0 && _h > options.HMAX ) _h = options.HMAX;

    // Store intermediate point
    if( store ) _mesh_sta.add( _istg, _t, _vec_sta, false );
  }
  return NORMAL;
}

template <typename T> inline typename ODESLV_GSL<T>::STATUS
ODESLV_GSL<T>::_states
( const unsigned ns, const double*tk, const double*p, double**xk,
  double*q, double*f, const bool store, std::ostream&os )
{
  try{
    // Initialize trajectory integration with GSL
    if( !_INI_STA( p, ns ) )
      { _END_STA( stats_sta ); return FATAL; }

    // States initial conditions
    _t = tk[0];
    if( !_IC_STA() )
      { _END_STA( stats_sta ); return FATAL; } 
    if( options.DISPLAY >= 1 ){
      _print_interm( _t, _nx, _vec_sta, "x", os );
      _print_interm( _nq, _vec_sta+_nx, "q", os );
    }
    if( xk && !xk[0] ) xk[0] = new double[_nx];
    for( unsigned ix=0; xk[0] && ix<_nx; ix++ ) xk[0][ix] = _vec_sta[ix];

    // Integrate ODEs through each stage using GSL
    _h = options.H0;
    pODESLV_GSL = this;

    for( _istg=0; _istg<ns; _istg++ ){

      // Reinitialize states at intermediate times
      if( !_CC_STA() )
        { _END_STA( stats_sta ); return FATAL; } 

      // Update list of operations in RHS and JAC
      const unsigned iRHS = ( _vRHS.size()<2? 0: _istg );
      if( (!_istg || iRHS) && !_RHS_STA_SET( iRHS ) )
        { _END_STA( stats_sta ); return FATAL; }

      // Integrate until end of time stage
      if( _states_traj( tk[_istg+1], store ) == FAILURE )
        { _END_STA( stats_sta ); return FAILURE; }
      if( options.DISPLAY >= 1 ){
        _print_interm( _t, _nx, _vec_sta, "x", os );
        _print_interm( _nq, _vec_sta+_nx, "q", os );
      }
      if( xk && !xk[_istg+1] ) xk[_istg+1] = new double[_nx];
      for( unsigned ix=0; xk[_istg+1] && ix<_nx; ix++ ) xk[_istg+1][ix] = _vec_sta[ix];

      // Add intermediate function terms
      const unsigned iFCT = ( _vFCT.size()>=ns? _istg:0 );
      if( (_vFCT.size()>=ns || _istg==ns-1) && !_FCT_STA( iFCT, f ) )
        { _END_STA( stats_sta ); return FATAL; }
    }

    // Final quadrature and function values
    for( unsigned iq=0; q && iq<_nq; iq++ ) q[iq] = _vec_sta[_nx+iq];
    if( options.DISPLAY >= 1 ) _print_interm( _nf, f, "f", os );
  }
  catch(...){
    _END_STA( stats_sta );
    if( options.DISPLAY >= 1 ) _print_stats( stats_sta, os );
    return FAILURE;
  }

  _END_STA( stats_sta);
  if( options.DISPLAY >= 1 ) _print_stats( stats_sta, os );
  return NORMAL;
}

template <typename T> inline typename ODESLV_GSL<T>::STATUS
ODESLV_GSL<T>::_states
( const unsigned ns, const double*tk, const T*Ip, T**Ixk, T*Iq, T*If,
  const unsigned nsamp, unsigned* vsamp, const unsigned ipar,
  double*p, double**xk, double*q, double*f, std::ostream&os )
{
  STATUS flag = NORMAL;

  // Update bounds for all sampling points
  for( unsigned isamp=0; isamp<nsamp; isamp++ ){
    vsamp[ipar] = isamp;

    // Continue recursive call
    if( ipar+1 < _np ){
      flag = _states( ns, tk, Ip, Ixk, Iq, If, nsamp, vsamp, ipar+1, p, xk, q, f, os );
      if( flag != NORMAL ) return flag;
      continue;
    }

    // Update bounds for current point
#ifdef MC__ODESLV_GSL_SAMPLE_DEBUG
    std::cout << "Sample: ";
#endif
    for( unsigned ip=0; ip<_np; ip++ ){
      p[ip] = Op<T>::l( Ip[ip] ) + vsamp[ip]/(nsamp-1.) * Op<T>::diam( Ip[ip] );
#ifdef MC__ODESLV_GSL_SAMPLE_DEBUG
      std::cout << p[ip] << "  ";
#endif
    }
#ifdef MC__ODESLV_GSL_SAMPLE_DEBUG
    std::cout << std::endl;
#endif
    flag = states( ns, tk, p, xk, q, f, os );
    if( flag != NORMAL ) return flag;
    for( unsigned is=0; Ixk && is<=ns; is++ )
      for( unsigned ix=0; ix<_nx; ix++ )
        Ixk[is][ix] = Op<T>::hull( xk[is][ix], Ixk[is][ix] );
    for( unsigned iq=0; Iq && iq<_nq; iq++ )
      Iq[iq] = Op<T>::hull( q[iq], Iq[iq] );
    for( unsigned ifn=0; If && ifn<_nf; ifn++ )
      If[ifn] = Op<T>::hull( f[ifn], If[ifn] );
  }

  return flag;
}  

//! @fn template <typename T> inline typename ODESLV_GSL<T>::STATUS ODESLV_GSL<T>::states(
//! const unsigned ns, const double*tk, const double*p, double**xk,
//! double*q, double*f, std::ostream&os )
//!
//! This function computes the solution of the parametric ODEs defined in IVP:
//!   - <a>ns</a> [input]  number of time stages
//!   - <a>tk</a> [input]  stage times, including the initial time
//!   - <a>p</a>  [input]  parameter values
//!   - <a>xk</a> [output] state values at stage times
//!   - <a>q</a>  [output] quadrature values at final time
//!   - <a>f</a>  [output] function values
//!   - <a>os</a> [input]  output stream
//!   .
//! The return value is the status.
template <typename T> inline typename ODESLV_GSL<T>::STATUS
ODESLV_GSL<T>::states
( const unsigned ns, const double*tk, const double*p, double**xk,
  double*q, double*f, std::ostream&os )
{
  return _states( ns, tk, p, xk, q, f, false, os );
}

//! @fn template <typename T> inline typename ODESLV_GSL<T>::STATUS ODESLV_GSL<T>::bounds(
//! const unsigned ns, const double*tk, const T*Ip, T**Ixk, T*q, T*f,
//! const unsigned nsamp, std::ostream&os )
//!
//! This function computes an approximate interval enclosure of the
//! reachable set of the parametric ODEs defined in IVP using equally
//! spaced samples:
//!   - <a>ns</a> [input] number of time stages
//!   - <a>tk</a> [input] stage times, including the initial time
//!   - <a>Ip</a> [input] interval parameter set
//!   - <a>Ixk</a> [output] approximate interval state enclosures at stage times
//!   - <a>Iq</a>  [output] approximate quadrature enclosures at final time
//!   - <a>If</a>  [output] approximate function enclosures
//!   - <a>nsamp</a> [input] number of samples for each parameter
//!   - <a>os</a> [input] output stream
//! .
//! The return value is the status.
template <typename T> inline typename ODESLV_GSL<T>::STATUS
ODESLV_GSL<T>::bounds
( const unsigned ns, const double*tk, const T*Ip, T**Ixk,
  T*Iq, T*If, const unsigned nsamp, std::ostream&os )
{
  int DISPLAY_SAVE = options.DISPLAY;
  options.DISPLAY = 0;
  STATUS flag = NORMAL;
  
  // Initialization of sampled bounds at parameter lower bound
  double *p = new double[_np];
  for( unsigned ip=0; ip<_np; ip++ )
    p[ip] = Op<T>::l(Ip[ip]);
  double **xk = Ixk? new double*[ns+1]: 0;
  for( unsigned is=0; Ixk && is<=ns; is++ ){
    if( !Ixk[is] ) Ixk[is] = new T[ns];
    xk[is] = new double[_nx];
  }
  double *q = Iq? new double[_nq]: 0;
  double *f = If? new double[_nf]: 0;
  flag = states( ns, tk, p, xk, q, f, os );
  if( flag != NORMAL || nsamp <= 1 ){
    delete[] p; delete[] f;
    for( unsigned is=0; is<=ns; is++ ) delete[] xk[is]; delete[] xk;
    return flag;
  }   
  for( unsigned is=0; Ixk && is<=ns; is++ )
    for( unsigned ix=0; ix<_nx; ix++ )
      Ixk[is][ix] = xk[is][ix];
  for( unsigned iq=0; Iq && iq<_nq; iq++ )
    Iq[iq] = q[iq];
  for( unsigned ifn=0; If && ifn<_nf; ifn++ )
    If[ifn] = f[ifn];

  // Start sampling process
  unsigned* vsamp = new unsigned[_np];
  flag = _states( ns, tk, Ip, Ixk, Iq, If, nsamp, vsamp, 0, p, xk, q, f, os );

  // Display results
  options.DISPLAY = DISPLAY_SAVE;
  if( options.DISPLAY >= 1 ){
    for( unsigned is=0; Ixk && is<=ns; is++ )
      _print_interm( tk[is], _nx, Ixk[is], "x", os );
    if( Iq ) _print_interm( _nq, Iq, "q", os );
    if( If ) _print_interm( _nf, If, "f", os );
  }

  // Record intermediate results
  _results_sta.clear();
  if( options.RESRECORD )
    for( unsigned is=0; Ixk && is<=ns; is++ )
      _results_sta.push_back( Results( tk[is], _nx, Ixk[is] ) );
  
  // Clean-up
  delete[] p; delete[] q; delete[] f;
  for( unsigned is=0; xk && is<=ns; is++ ) delete[] xk[is]; delete[] xk;
  delete[] vsamp;
  
  return flag;
}

template <typename T> template<typename U> inline void
ODESLV_GSL<T>::_print_interm
( const double t, const unsigned nx, const U*x, const std::string&var,
  std::ostream&os )
{
  os << " @t = " << std::scientific << std::setprecision(6)
                 << std::left << t << " :" << std::endl;
  for( unsigned ix=0; ix<nx; ix++ )
    os << " " << var.c_str() << "[" << ix << "] = " << x[ix] << std::endl;
  return;
}

template <typename T> template<typename U> inline void
ODESLV_GSL<T>::_print_interm
( const unsigned nx, const U*x, const std::string&var, std::ostream&os )
{
  if( !nx || !x ) return;
  os << std::scientific << std::setprecision(6) << std::left;
  for( unsigned ix=0; ix<nx; ix++ )
    os << " " << var.c_str() << "[" << ix << "] = " << x[ix] << std::endl;
  return;
}

template <typename T> inline void
ODESLV_GSL<T>::record
( std::ofstream&bndsta, const unsigned iprec ) const
{
  if( !bndsta ) return;

  // Specify format
  bndsta << std::right << std::scientific << std::setprecision(iprec);

  // Record computed state interval bounds at stage times
  typename std::vector< Results >::const_iterator it = _results_sta.begin();
  for( ; it != _results_sta.end(); ++it ){
    bndsta << std::setw(iprec+9) << (*it).t;
    for( unsigned ix=0; ix<(*it).nx; ix++ )
      bndsta << std::setw(iprec+9) << mc::Op<T>::l( (*it).X[ix] )
             << std::setw(iprec+9) << mc::Op<T>::u( (*it).X[ix] );
    bndsta << std::endl;
  }
}

} // end namescape mc

#endif

