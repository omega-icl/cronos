// Copyright (C) 2015 Benoit Chachuat & Nikola Peric, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

#ifndef MC__ODEBNDS_GSL_HPP
#define MC__ODEBNDS_GSL_HPP

#include "odebnds_base_NP.hpp"
#include "odebnd_gsl.hpp"

#define MC__ODEBNDS_GSL_USE_BAD

//! @TODO
//! - handle state quadratures in adjoint sensitivity ==> OK
//! - handle state discontinuities in adjoint sensitivity 
//! - propagation of polynomial models in adjoint sensitivity


//! @brief C++ class computing enclosures of the reachable set of parametric ODEs and adjoint sensitivity analysis using GSL and MC++.
////////////////////////////////////////////////////////////////////////
//! mc::ODEBNDS_GSL is a C++ class that computes enclosures of the
//! reachable set of parametric ordinary differential equations
//! (ODEs) using GSL and MC++. It implements the methods of differential 
//! inequalities, whereby polynomial models with interval or ellipsoidal
//! remainders are used to enable high-order convergence. In addition,
//! the class computes enclosures on the adjoint system of equations,
//! which can be used to calculate the derivatives of state-dependent
//! functionals. The use of ellipsoidal remainders enables stability
//! of the enclosures for asymptotically stable ODE systems when the
//! parameter host is sufficiently small.
////////////////////////////////////////////////////////////////////////
namespace mc
{
template < typename T, typename PMT=mc::TModel<T>, typename PVT=mc::TVar<T> >
class ODEBNDS_GSL:
  public virtual BASE_DE,
  public virtual BASE_GSL,
  public virtual ODEBNDS_BASE<T,PMT,PVT>,
  public virtual ODEBND_GSL<T,PMT,PVT>
{
  
  typedef Ellipsoid E;
  typedef BASE_DE::STATUS STATUS;
  using ODEBND_BASE<T,PMT,PVT>::NORMAL; 
  using ODEBND_BASE<T,PMT,PVT>::FAILURE;
  using ODEBND_BASE<T,PMT,PVT>::FATAL;

  using ODEBND_BASE<T,PMT,PVT>::_Q;
  using ODEBND_BASE<T,PMT,PVT>::_Er;
  using ODEBND_BASE<T,PMT,PVT>::_Ir;
  using ODEBND_BASE<T,PMT,PVT>::_pref;
  using ODEBND_BASE<T,PMT,PVT>::_Ip;
  using ODEBND_BASE<T,PMT,PVT>::_xref;
  using ODEBND_BASE<T,PMT,PVT>::_B;
  using ODEBND_BASE<T,PMT,PVT>::_diam;
  using ODEBND_BASE<T,PMT,PVT>::_vec2I;
  using ODEBND_BASE<T,PMT,PVT>::_vec2E;
  using ODEBND_BASE<T,PMT,PVT>::_PMenv;
  using ODEBND_BASE<T,PMT,PVT>::_PMq;
  using ODEBND_BASE<T,PMT,PVT>::_vec2PMI;
  using ODEBND_BASE<T,PMT,PVT>::_vec2PME;
  using ODEBND_BASE<T,PMT,PVT>::_print_interm;

  using ODEBNDS_BASE<T,PMT,PVT>::_Ix;
  using ODEBNDS_BASE<T,PMT,PVT>::_Iy;
  using ODEBNDS_BASE<T,PMT,PVT>::_Idy;
  using ODEBNDS_BASE<T,PMT,PVT>::_Iyq;
  using ODEBNDS_BASE<T,PMT,PVT>::_Qy;
  using ODEBNDS_BASE<T,PMT,PVT>::_By;
  using ODEBNDS_BASE<T,PMT,PVT>::_Edy;
  using ODEBNDS_BASE<T,PMT,PVT>::_yref;
  using ODEBNDS_BASE<T,PMT,PVT>::_PMz;
  using ODEBNDS_BASE<T,PMT,PVT>::_PMy;
  using ODEBNDS_BASE<T,PMT,PVT>::_PMyq;
  using ODEBNDS_BASE<T,PMT,PVT>::_vec2E;
  using ODEBNDS_BASE<T,PMT,PVT>::_IC_I_ADJ;
  using ODEBNDS_BASE<T,PMT,PVT>::_TC_I_ADJ;
  using ODEBNDS_BASE<T,PMT,PVT>::_TC_I_QUAD;
  using ODEBNDS_BASE<T,PMT,PVT>::_CC_I_ADJ;
  using ODEBNDS_BASE<T,PMT,PVT>::_CC_I_QUAD;
  using ODEBNDS_BASE<T,PMT,PVT>::_SET_I_ADJ;
  using ODEBNDS_BASE<T,PMT,PVT>::_RHS_I_ADJ;
  using ODEBNDS_BASE<T,PMT,PVT>::_RHS_I_QUAD;
  using ODEBNDS_BASE<T,PMT,PVT>::_IC_PM_ADJ;
  using ODEBNDS_BASE<T,PMT,PVT>::_TC_PM_ADJ;
  using ODEBNDS_BASE<T,PMT,PVT>::_TC_PM_QUAD;
  using ODEBNDS_BASE<T,PMT,PVT>::_CC_PM_ADJ;
  using ODEBNDS_BASE<T,PMT,PVT>::_CC_PM_QUAD;
  using ODEBNDS_BASE<T,PMT,PVT>::_SET_PM_ADJ;
  using ODEBNDS_BASE<T,PMT,PVT>::_RHS_PM_ADJ;
  using ODEBNDS_BASE<T,PMT,PVT>::_RHS_PM_QUAD;

  using ODEBND_GSL<T,PMT,PVT>::_vec_sta;
  using ODEBND_GSL<T,PMT,PVT>::_mesh_sta;
  using ODEBND_GSL<T,PMT,PVT>::_pos_rhs;
  using ODEBND_GSL<T,PMT,PVT>::_pos_quad;
  using ODEBND_GSL<T,PMT,PVT>::_pos_fct;
  using ODEBND_GSL<T,PMT,PVT>::_offset_quad;
  using ODEBND_GSL<T,PMT,PVT>::_init_stats;
  using ODEBND_GSL<T,PMT,PVT>::_final_stats;
  using ODEBND_GSL<T,PMT,PVT>::_print_stats;

 template <typename U, typename PMU, typename PVU>
 friend int MC_GSLADJRHSI__
   ( double t, const double* y, double* ydot, void* user_data );

 template <typename U, typename PMU, typename PVU>
 friend int MC_GSLADJRHSPM__
   ( double t, const double* y, double* ydot, void* user_data );

private:
  //! @brief GSL drivers for adjoint ODE integration
  std::vector<gsl_odeiv2_driver*> _driver_adj;

  //! @brief GSL data type for adjoint ODE integration
  gsl_odeiv2_system _sys_adj;

protected:

  //! @brief current function
  unsigned _ifct;

  //! @brief vector storing stepsize during adjoint integration
  std::vector<double> _h_adj; // In GSL and SUNDIALS headers seperately

  //! @brief array for adjoint GSL integration
  double *_vec_adj; // In GSL and SUNDIALS headers seperately
  
  //! @brief position in _vec_adj
  unsigned _pos_adj; // In GSL and SUNDIALS headers seperately

  //! @brief static pointer to class
  static ODEBNDS_GSL<T,PMT,PVT> *pODEBNDS;
  
public: // EVERYTHING to be called/used from outside (mostly functions and a few variables)
  //typedef BASE_GSL::STATUS STATUS;
  typedef BASE_GSL::Stats Stats;
  typedef typename ODEBND_BASE<T,PMT,PVT>::Results Results;
  typedef typename ODEBND_GSL<T,PMT,PVT>::Exceptions Exceptions;

  //! @brief Default constructor
  ODEBNDS_GSL();

  //! @brief Virtual destructor
  virtual ~ODEBNDS_GSL();

  //! @brief Integrator options
  struct Options: public ODEBND_GSL<T,PMT,PVT>::Options
  {
    //! @brief Constructor
    Options():
      ODEBND_GSL<T,PMT,PVT>::Options(), INTERPMETH(MESH_GSL::CSPLINE)
      {}
    //! @brief Assignment operator
    template <typename U> Options& operator=
      ( U&options ){
        ODEBND_GSL<T,PMT,PVT>::Options::operator=(options);
        INTERPMETH   = options.INTERPMETH;
        return *this;
      }
    //! @brief Numerical interpolation method
    MESH_GSL::INTERPOLATION_METHOD INTERPMETH;
  } options;

  //! @brief Statistics for adjoint integration
  Stats stats_adj;

  //! @brief Vector storing interval adjoint bounds (see Options::RESRECORD)
  std::vector< Results > _results_adj;

  //! @brief Propagate state/quadrature interval bounds forward in time through every time stages
  STATUS bounds
    ( const unsigned ns, const double*tk, const T*Ip, T**Ixk,
      T*Iq, T*If, std::ostream&os=std::cout )
      { ODEBND_GSL<T,PMT,PVT>::options = options;
        return ODEBND_GSL<T,PMT,PVT>::_bounds( ns, tk, Ip, Ixk, Iq, If, false, os); }

  //! @brief Propagate state/quadrature polynomial models forward in time through every time stages
  STATUS bounds
    ( const unsigned ns, const double*tk, const PVT*PMp, PVT**PMxk,
      PVT*PMq, PVT*PMf, std::ostream&os=std::cout )
      { ODEBND_GSL<T,PMT,PVT>::options = options;
        return ODEBND_GSL<T,PMT,PVT>::_bounds( ns, tk, PMp, PMxk, PMq, PMf, false, os); }

  //! @brief Propagate state and adjoint interval bounds forward and backward in time through every time stages
  STATUS bounds_ASA
    ( const unsigned ns, const double*tk, const T*Ip, T**Ixk,
      T*Iq, T*If, T**Ilk, T*Idf, std::ostream&os=std::cout );

  //! @brief Propagate state and adjoint polynomial models forward and backward in time through every time stages
  STATUS bounds_ASA
    ( const unsigned ns, const double*tk, const PVT*PMp, PVT**PMxk,
      PVT*PMq, PVT*PMf, PVT**PMlk, PVT*PMdf, std::ostream&os=std::cout );

  //! @brief Record state and sensitivity bounds in files <a>obndsta</a> and <a>obndsa</a>, with accuracy of <a>iprec</a> digits
  void record
    ( std::ofstream&obndsta, std::ofstream&obndsa, const unsigned iprec=5 ) const
    { this->ODEBND_GSL<T,PMT,PVT>::record( obndsta, iprec );
      this->ODEBND_BASE<T,PMT,PVT>::_record( obndsa, _results_adj, iprec ); }

  //! @brief Record state bounds in file <a>obndsta</a>, with accuracy of <a>iprec</a> digits
  void record
    ( std::ofstream&obndsta, const unsigned iprec=5 ) const
    { this->ODEBND_GSL<T,PMT,PVT>::record( obndsta, iprec ); }

private:

  //! @brief Static wrapper to function to calculate the adjoint DINEQ RHS values
  static int MC_GSLADJRHSI__
    ( double t, const double* y, double* ydot, void* user_data );

  //! @brief Function to initialize adjoint interval bounding
  bool _INI_I_ADJ //split
    ( const unsigned np, const T *Ip );

  //! @brief Static wrapper to function to calculate the adjoint DINEQ-PM RHS values
  static int MC_GSLADJRHSPM__
    ( double t, const double* y, double* ydot, void* user_data );

  //! @brief Function to initialize GSL for adjoint polynomial models
  bool _INI_PM_ADJ //split
    ( const unsigned np, const PVT*PMp );

  //! @brief Function to initialize GSL numerical integration drivers
  void _INI_GSL
    ( gsl_odeiv2_system &sys, std::vector<gsl_odeiv2_driver*>&driver );

  //! @brief Function to finalize adjoint bounding
  void _END_ADJ();

};

template <typename T, typename PMT, typename PVT>
 ODEBNDS_GSL<T,PMT,PVT>* ODEBNDS_GSL<T,PMT,PVT>::pODEBNDS = 0;

template <typename T, typename PMT, typename PVT> inline
ODEBNDS_GSL<T,PMT,PVT>::ODEBNDS_GSL
()
: BASE_DE(), BASE_GSL(), ODEBNDS_BASE<T,PMT,PVT>(), ODEBND_GSL<T,PMT,PVT>(),
  _ifct(0), _vec_adj(0)
{}

template <typename T, typename PMT, typename PVT> inline
ODEBNDS_GSL<T,PMT,PVT>::~ODEBNDS_GSL
()
{
  // Free GSL array
  delete[] _vec_adj;
  
  // Free GSL arrays
  for( unsigned i=0; i<_driver_adj.size(); i++ )
    {if( _driver_adj[i] )  gsl_odeiv2_driver_free( _driver_adj[i] );}
}

template <typename T, typename PMT, typename PVT> inline void
ODEBNDS_GSL<T,PMT,PVT>::_INI_GSL
( gsl_odeiv2_system &sys, std::vector<gsl_odeiv2_driver*>&driver )
{
  // Reset GSL numerical integration drivers
  for( unsigned i=0; i<driver.size(); i++ )
    gsl_odeiv2_driver_free( driver[i] );
  driver.clear();

  for( unsigned i=0; i<_nf; i++ ){
    switch( options.INTMETH ){
    case Options::RK8PD:
      driver.push_back( gsl_odeiv2_driver_alloc_y_new( &sys,
        gsl_odeiv2_step_rk8pd, options.H0, options.ATOL, options.RTOL ) );
      break;
    case Options::MSADAMS:
      driver.push_back( gsl_odeiv2_driver_alloc_y_new( &sys,
        gsl_odeiv2_step_msadams, options.H0, options.ATOL, options.RTOL ) );
      break;
    case Options::MSBDF:
      driver.push_back( gsl_odeiv2_driver_alloc_y_new( &sys,
        gsl_odeiv2_step_msbdf, options.H0, options.ATOL, options.RTOL ) );
      break;
    case Options::RKF45: default:
      driver.push_back( gsl_odeiv2_driver_alloc_y_new( &sys,
        gsl_odeiv2_step_rkf45, options.H0, options.ATOL, options.RTOL ) );
      break;
    }
    gsl_odeiv2_driver_set_hmin( driver[i], options.HMIN );  
    gsl_odeiv2_driver_set_nmax( driver[i], options.NMAX );  
  }
  return;
}

template <typename T, typename PMT, typename PVT> inline void
ODEBNDS_GSL<T,PMT,PVT>::_END_ADJ()
{
  // Get final CPU time
  _final_stats( stats_adj );
}

template <typename T, typename PMT, typename PVT> inline bool
ODEBNDS_GSL<T,PMT,PVT>::_INI_I_ADJ
( const unsigned np, const T* Ip )
{ 
  // Initialize bound propagation
  if( !ODEBNDS_BASE<T,PMT,PVT>::_INI_I_ADJ( options, np, Ip ) )
    return false;

  // Define adjoint ODE system in GSL format
  _sys_adj.function = MC_GSLADJRHSI__;
  _sys_adj.params = 0;
  _sys_adj.dimension = 2*np;
  switch( options.WRAPMIT){
  case Options::NONE:
  case Options::DINEQ:
    _sys_adj.dimension += 2*_nx;
    break;
  case Options::ELLIPS:
  default:
    _sys_adj.dimension += _nx*(1+np)+_nx*(_nx+1)/2;
    break;
  }

  // Set GSL drivers for adjoint ODE integration
  _INI_GSL( _sys_adj, _driver_adj );
  delete [] _vec_adj;
  _vec_adj  = new double[ _sys_adj.dimension*_nf ];
  _offset_quad = _sys_adj.dimension - 2*np;

  // Reset result record and statistics
  _results_adj.clear();
  _init_stats( stats_adj );

  return true;
}

template <typename T, typename PMT, typename PVT> inline int
ODEBNDS_GSL<T,PMT,PVT>::MC_GSLADJRHSI__
( double t, const double* y, double* ydot, void* user_data )
{
  ODEBNDS_GSL<T,PMT,PVT> *pODEBNDS = ODEBNDS_GSL<T,PMT,PVT>::pODEBNDS;
  if( !pODEBNDS->_mesh_sta.eval( pODEBNDS->_istg, -t, pODEBNDS->_vec_sta ) )
    { return GSL_EBADFUNC; } // set interpolated state
  bool flag = pODEBNDS->_RHS_I_ADJ( pODEBNDS->options, -t, y, ydot,
                                    pODEBNDS->_vec_sta, pODEBNDS->_ifct );
  if( flag ){
    double* qdot = ydot + pODEBNDS->_offset_quad;
    flag = pODEBNDS->_RHS_I_QUAD( pODEBNDS->options, -t, y, qdot,
                                  pODEBNDS->_vec_sta, pODEBNDS->_ifct );
  }
  pODEBNDS->stats_adj.numRHS++; // increment RHS counter
  ODEBNDS_GSL<T,PMT,PVT>::pODEBNDS = pODEBNDS;
  return ( flag? GSL_SUCCESS: GSL_EBADFUNC );
}

//! @fn template <typename T, typename PMT, typename PVT> inline typename ODEBND_GSL<T,PMT,PVT>::STATUS ODEBND_GSL<T,PMT,PVT>::bounds_ASA
//!( const unsigned ns, const double*tk, const T*Ip, T**Ixk,
//!  T*Iq, T*If, T**Ilk, T*Idf, std::ostream&os=std::cout )
//!
//! This function computes an interval enclosure on the functional defined
//! in F together with its derivatives as well as an interval enclosure of the reachable set of 
//! the parametric ODEs defined in IVP using equally spaced samples:
//!   - <a>ns</a> [input] number of time stages
//!   - <a>tk</a> [input] stage times, including the initial time
//!   - <a>Ip</a> [input] parameters interval parameter set
//!   - <a>Ixk</a> [output] states interval enclosures at stage times
//!   - <a>If</a> [output] functions interval enclosures at stage times
//!   - <a>Idf</a> [output] function derivatives interval enclosures at stage times
//!   - <a>os</a> [input] output stream
//! .
//! The return value is the status.
template <typename T, typename PMT, typename PVT>
inline typename ODEBNDS_GSL<T,PMT,PVT>::STATUS
ODEBNDS_GSL<T,PMT,PVT>::bounds_ASA
( const unsigned ns, const double*tk, const T*Ip, T**Ixk, T*Iq, T*If,
  T**Ilk, T*Idf, std::ostream&os )
{
  // Compute state bounds and store intermediate results in _mesh_sta
  STATUS flag = NORMAL;
  ODEBND_GSL<T,PMT,PVT>::options = options;
  flag = ODEBND_GSL<T,PMT,PVT>::_bounds( ns, tk, Ip, Ixk, Iq, If, true, os);
  if( flag != NORMAL ) return flag;

  // Nothing to do if no functions are defined
  if( !_nf ) return NORMAL;

  // Check size
  if( !Ilk || !Idf ) return FATAL;

  try{
    // Initialize adjoint bound integration using GSL
    if( !_INI_I_ADJ( _np, Ip )) return FATAL;

    // Interpolate state mesh in final stage
    _t = -tk[ns];
    if( !_mesh_sta.interp( ns, options.INTERPMETH ) )
      { _END_ADJ(); return FAILURE; }

    // Bounds on terminal states/quadratures
    if( !_mesh_sta.eval( ns, tk[ns], _vec_sta ) )
      { _END_ADJ(); return FAILURE; }
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

    // Bounds on terminal adjoints/quadratures
    if( Ilk && !Ilk[ns] ) Ilk[ns] = new T[_nx*_nf];
    for( _ifct=_pos_adj=0; _ifct < _nf; _ifct++, _pos_adj+=_sys_adj.dimension ){
      _pos_fct = ( _vFCT.size()>=ns? ns-1:0 );
      if( !_TC_I_ADJ( options, -_t, _vec_adj+_pos_adj, _pos_fct, _ifct )
       || !_TC_I_QUAD( options, -_t, _vec_adj+_pos_adj+_offset_quad ) )
        { _END_ADJ(); return FATAL; }
      for( unsigned iy=0; Ilk[ns] && iy<_nx; iy++ )
        Ilk[ns][_ifct*_nx+iy] = _Iy[iy];
    }

    // Display & record adjoint terminal results
    if( options.DISPLAY >= 1 ){
      _print_interm( tk[ns], _nx*_nf, Ilk[ns], "l", os );
      _print_interm( _np, _Iyq, "q", os );
    }
    if( options.RESRECORD )
      _results_adj.push_back( Results( tk[ns], _nf*_nx, Ilk[ns] ) );

    // Integrate adjoint ODEs through each stage using GSL
    _h_adj.assign( _nf, options.H0 );
    pODEBNDS = this;
    for( _istg=ns; _istg>0; _istg-- ){

      // Interpolate state mesh in current stage
      if( _istg<ns && !_mesh_sta.interp( _istg, options.INTERPMETH ) )
        { _END_ADJ(); return FAILURE; }

      // Update list of operations in RHSADJ and QUADADJ
      _pos_rhs  = ( _vRHS.size() <=1? 0: _istg-1 );
      _pos_quad = ( _vQUAD.size()<=1? 0: _istg-1 );
      if( !_SET_I_ADJ( options, _pos_rhs, _pos_quad, _pos_fct, false ) )
        { _END_ADJ(); return FATAL; }

      // Integrate backward through current stage for each function
      for( _ifct=_pos_adj=0; _ifct < _nf; _ifct++, _pos_adj+=_sys_adj.dimension ){
        _t = -tk[_istg];
#ifdef MC__ODEBNDS_GSL_DINEQI_DEBUG
        std::cout << "@t=" << -_t << std::endl;
        std::cout << "_pos_adj = " << _pos_adj << std::endl;
        for( unsigned iy=0; iy<_sys_adj.dimension; iy++ )
          std::cout << "_vec_adj[" << _pos_adj+iy << "] = "
                    << _vec_adj[_pos_adj+iy] << std::endl;
        { int dum; std::cin >> dum; }
#endif

        // Propagate bounds backward to previous stage time
        while( _t < -tk[_istg-1] ){
          if( gsl_odeiv2_evolve_apply( _driver_adj[_ifct]->e, _driver_adj[_ifct]->c,
              _driver_adj[_ifct]->s, &_sys_adj, &_t, -tk[_istg-1], &_h_adj[_ifct],
              _vec_adj+_pos_adj ) != GSL_SUCCESS
           || _h_adj[_ifct] < options.HMIN
           || (options.NMAX && stats_adj.numSteps > options.NMAX)
           || _diam(_nx, _Iy) > options.DMAX )
            throw Exceptions( Exceptions::INTERN );
          stats_adj.numSteps++;    
          if( options.HMAX > 0 && _h_adj[_ifct] > options.HMAX ) _h_adj[_ifct] = options.HMAX;
        }
#ifdef MC__ODEBNDS_GSL_DINEQI_DEBUG
        std::cout << "@t=" << -_t << std::endl;
        for( unsigned iy=0; iy<_sys_adj.dimension; iy++ )
          std::cout << "_vec_adj[" << _pos_adj+iy << "] = " << _vec_adj[_pos_adj+iy] << std::endl;
        { int dum; std::cin >> dum; }
#endif

        // Bounds on states/adjoints/quadratures at stage time
        if( !_mesh_sta.eval( _istg, tk[_istg-1], _vec_sta ) )
          { _END_ADJ(); return FATAL; }
        switch( options.WRAPMIT){
        case Options::NONE:
        case Options::DINEQ:
          _vec2I( _vec_sta, _nx, _Ix );
          _vec2I( _vec_adj+_pos_adj, _nx, _Iy);
          break;
        case Options::ELLIPS:
        default:
          _vec2E( _vec_sta, _nx, _np, _Q, _Er, _Ir, _pref, _Ip, _B, _xref, _Ix );
          _vec2E( _vec_adj+_pos_adj, _nx, _np, _Qy, _Edy, _Idy, _pref, _Ip, _By, _yref, _Iy); 
          break;
        }
        _vec2I( _vec_adj+_pos_adj+_offset_quad, _np, _Iyq);
        if( Ilk && !Ilk[_istg-1] ) Ilk[_istg-1] = new T[_nx*_nf];
        for( unsigned iy=0; Ilk[_istg-1] && iy<_nx; iy++ )
          Ilk[_istg-1][_ifct*_nx+iy] = _Iy[iy];

        // Add function contribution to adjoint bounds (discontinuities)
        if( _istg>1  ){
          _pos_fct = ( _vFCT.size()>=ns? _istg-1:0 );
          if( _pos_fct
           && !_CC_I_ADJ( options, -_t, _vec_adj+_pos_adj, _pos_fct, _ifct )
           && !_CC_I_QUAD( options, -_t, _vec_adj+_pos_adj+_offset_quad ) )
            { _END_ADJ(); return FATAL; }
            // Reset ODE solver - needed in case of discontinuity (used to be in _CC_I_ADJ)
            gsl_odeiv2_driver_reset( _driver_adj[_ifct] );
        }

        // Add initial state contribution to derivative bounds
        else if( Idf ){
          if( !_IC_I_ADJ( -_t ) ){ _END_ADJ(); return FATAL; }
          for( unsigned iq=0; iq<_np; iq++ )
            Idf[_ifct*_np+iq] = _Iyq[iq];
        }
      }

      // Display & record adjoint intermediate results
      if( options.DISPLAY >= 1 ){
        _print_interm( tk[_istg-1], _nf*_nx, Ilk[_istg-1], "l", os );
        _print_interm( _np, _Iyq, "q", os );
      }
      if( options.RESRECORD )
         _results_adj.push_back( Results( tk[_istg-1], _nf*_nx, Ilk[_istg-1] ) );
    }

    if( options.DISPLAY >= 1 )
      _print_interm( _nf*_np, Idf, "df", os );
  }
  catch(...){
    _END_ADJ();
    if( options.DISPLAY >= 1 ) _print_stats( stats_adj, os );
    return FAILURE;
  }
  _END_ADJ();
  if( options.DISPLAY >= 1 ) _print_stats( stats_adj, os );

  return NORMAL;
}

template <typename T, typename PMT, typename PVT> inline bool
ODEBNDS_GSL<T,PMT,PVT>::_INI_PM_ADJ
( const unsigned np, const PVT* PMp )
{
  // Initialize bound propagation
  if( !ODEBNDS_BASE<T,PMT,PVT>::_INI_PM_ADJ( options, np, PMp ) )
    return false;

  // Define Adjoint system in GSL format
  _sys_adj.function = MC_GSLADJRHSPM__;
  _sys_adj.params = 0;
  switch( options.WRAPMIT){
  case Options::NONE:
    _sys_adj.dimension = (_PMenv->nmon()+1)*(_nx+np);
    break;
  case Options::DINEQ:
    _sys_adj.dimension = _PMenv->nmon()*(_nx+np) + 2*_nx + np;
    break;
  case Options::ELLIPS:
  default:
    _sys_adj.dimension = _PMenv->nmon()*(_nx+np) + _nx*(_nx+1)/2 + np;
    break;
  }

  // Set GSL drivers for adjoint ODE integration
  _INI_GSL(_sys_adj, _driver_adj);
  delete[] _vec_adj;
  _vec_adj = new double[_sys_adj.dimension*_nf];
  _offset_quad = _sys_adj.dimension - _PMenv->nmon()*np-np;

  // Reset result record and statistics
  _results_adj.clear();
  _init_stats( stats_adj );

  return true;
}

template <typename T, typename PMT, typename PVT> inline int
ODEBNDS_GSL<T,PMT,PVT>::MC_GSLADJRHSPM__
( double t, const double* y, double* ydot, void* user_data )
{
  ODEBNDS_GSL<T,PMT,PVT> *pODEBNDS = ODEBNDS_GSL<T,PMT,PVT>::pODEBNDS;
  if( !pODEBNDS->_mesh_sta.eval( pODEBNDS->_istg, -t, pODEBNDS->_vec_sta ) )
    { return GSL_EBADFUNC; } // set interpolated state
  bool flag = pODEBNDS->_RHS_PM_ADJ( pODEBNDS->options, -t, y, ydot,
                                     pODEBNDS->_vec_sta, pODEBNDS->_ifct );
  if( flag ){
    double* qdot = ydot + pODEBNDS->_offset_quad;
    flag = pODEBNDS->_RHS_PM_QUAD( pODEBNDS->options, -t, y, qdot,
                                   pODEBNDS->_vec_sta, pODEBNDS->_ifct );
  }
  pODEBNDS->stats_adj.numRHS++; // increment RHS counter
  ODEBNDS_GSL<T,PMT,PVT>::pODEBNDS = pODEBNDS;
  return ( flag? GSL_SUCCESS: GSL_EBADFUNC );
}

//! @fn template <typename T, typename PMT, typename PVT> inline typename ODEBNDS_GSL<T,PMT,PVT>::STATUS ODEBNDS_GSL<T,PMT,PVT>::bounds_ASA(
//! const unsigned ns, const double*tk, const PVT*PMp, PVT**PMxk=0, 
//! PVT*PMq=0, PVT*PMf=0, PVT**PMlk=0, PVT*PMdf=0, std::ostream&os=std::cout )
//!
//! This function computes an enclosure on the functional defined in F together with its
//! derivatives, as well as an enclosure of the reachable set of the parametric ODEs defined in
//! IVP using equally spaced samples, using the propagation of polynomial models with convex
//! remainders (intervals, ellipsoids):
//!  - <a>ns</a> [input] number of time stages
//!  - <a>tk</a> [input] stage times, including the initial time
//!  - <a>PMp</a> [input] polynomial model of parameter set
//!  - <a>PMxk</a> [output] polynomial model of state enclosures at stage times (default: NULL)
//!  - <a>PMq</a> [output] polynomial model of quadrature variables (default: NULL)
//!  - <a>PMf</a> [output] polynomial model of state/quadrature functionals (default: NULL)
//!  - <a>PMlk</a> [output] polynomial model of adjoint enclosures at stage times (default: NULL)
//!  - <a>PMdf</a> [output] polynomial model of functional derivatives (default: NULL)
//!  - <a>os</a> [input] output stream (default: std::cout)
//! .
//! The return value is the status.
template <typename T, typename PMT, typename PVT>
inline typename ODEBNDS_GSL<T,PMT,PVT>::STATUS
ODEBNDS_GSL<T,PMT,PVT>::bounds_ASA
( const unsigned ns, const double*tk, const PVT*PMp, PVT**PMxk,
  PVT*PMq, PVT*PMf, PVT**PMlk, PVT*PMdf, std::ostream&os )
{
  // Compute state bounds and store intermediate results in _mesh_sta
  STATUS flag = NORMAL;
  ODEBND_GSL<T,PMT,PVT>::options = options;
  flag = ODEBND_GSL<T,PMT,PVT>::_bounds( ns, tk, PMp, PMxk, PMq, PMf, true, os);
  if( flag != NORMAL ) return flag;

  // Nothing to do if no functions are defined
  if( !_nf ) return NORMAL;

  // Check size
  if( !PMlk || !PMdf ) return FATAL;

  try{
    // Initialize adjoint bound integration using GSL
    if( !_INI_PM_ADJ( _np, PMp ) ) return FATAL;
    _t = -tk[ns];

    // Interpolate state mesh in final stage
    if( !_mesh_sta.interp( ns, options.INTERPMETH ) )
      { _END_ADJ(); return FAILURE; }

    // Bounds on terminal states/quadratures
    if( !_mesh_sta.eval( ns, tk[ns], _vec_sta ) )
      { _END_ADJ(); return FAILURE; }
    switch( options.WRAPMIT){
    case Options::NONE:
      _vec2PMI( _vec_sta, _PMenv, _nx, _PMz, true );
      break;
    case Options::DINEQ:
      _vec2PMI( _vec_sta, _PMenv, _nx, _PMz );
      break;
    case Options::ELLIPS:
    default:
      _vec2PME( _vec_sta, _PMenv, _nx, _PMz, _Q, _Er, _Ir );
      break;
    }

    // Bounds on terminal adjoints/quadratures
    if( PMlk && !PMlk[ns] ) PMlk[ns] = new PVT[_nx*_nf];
    for( _ifct=_pos_adj=0; _ifct < _nf; _ifct++, _pos_adj+=_sys_adj.dimension ){
      _pos_fct = ( _vFCT.size()>=ns? ns-1:0 );
      if( !_TC_PM_ADJ( options, -_t, _vec_adj+_pos_adj, _pos_fct, _ifct )
       || !_TC_PM_QUAD( options, -_t, _vec_adj+_pos_adj+_offset_quad ) )
        { _END_ADJ(); return FATAL; }
      for( unsigned iy=0; PMlk[ns] && iy<_nx; iy++ )
        PMlk[ns][_ifct*_nx+iy] = _PMy[iy];
    }

    // Display & record adjoint terminal results
    if( options.DISPLAY >= 1 ){
      _print_interm( tk[ns], _nx*_nf, PMlk[ns], "l", os );
      _print_interm( _np, _PMyq, "q", os );
    }
    if( options.RESRECORD )
      _results_adj.push_back( Results( tk[ns], _nf*_nx, PMlk[ns] ) );

    // Integrate adjoint ODEs through each stage using GSL
    _h_adj.assign( _nf, options.H0 );
    pODEBNDS = this;
    for( _istg=ns; _istg>0; _istg-- ){

      // Interpolate state mesh in current stage
      if(_istg<ns && !_mesh_sta.interp(_istg, options.INTERPMETH) )
        {_END_ADJ(); return FAILURE;}

      // Update list of operations in RHSADJ and QUADADJ
      _pos_rhs  = ( _vRHS.size() <=1? 0: _istg-1 );
      _pos_quad = ( _vQUAD.size()<=1? 0: _istg-1 );
      if( !_SET_PM_ADJ( options, _pos_rhs, _pos_quad, _pos_fct, false ) )
        { _END_ADJ(); return FATAL; }

      // Integrate backward through current stage for each function
      for( _ifct=_pos_adj=0; _ifct < _nf; _ifct++, _pos_adj+=_sys_adj.dimension ){
        _t = -tk[_istg];
#ifdef MC__ODEBNDS_GSL_DINEQI_DEBUG
        std::cout << "@t=" << -_t << std::endl;
        std::cout << "Function Number " << _ifct << std::endl;
        for( unsigned iy=0; iy<_sys_adj.dimension; iy++ )
          std::cout << "_vec_adj[" << _pos_adj+iy << "] = "<<_vec_adj[_pos_adj+iy]<<std::endl;
        { int dum; std::cin >> dum; }
#endif

        // Propagate bounds backward to previous stage time
        while( _t < -tk[_istg-1] ){
          if( gsl_odeiv2_evolve_apply( _driver_adj[_ifct]->e, _driver_adj[_ifct]->c,
              _driver_adj[_ifct]->s, &_sys_adj, &_t, -tk[_istg-1], &_h_adj[_ifct],
              _vec_adj+_pos_adj ) != GSL_SUCCESS
           || _h_adj[_ifct] < options.HMIN
           || (options.NMAX && stats_adj.numSteps > options.NMAX)
           || _diam(_nx, _PMy) > options.DMAX )
            throw Exceptions( Exceptions::INTERN );
          stats_adj.numSteps++;
          if( options.HMAX > 0 && _h_adj[_ifct] > options.HMAX ) _h_adj[_ifct] = options.HMAX;
        }
#ifdef MC__ODEBNDS_GSL_DINEQI_DEBUG
        std::cout << "@t=" << -_t << std::endl;
        for( unsigned iy=0; iy<_sys_adj.dimension; iy++ )
          std::cout << "_vec_adj[" << _pos_adj+iy << "] = " << _vec_adj[_pos_adj+iy] << std::endl;
        { int dum; std::cin >> dum; }
#endif

        // Bounds on states/adjoints/quadratures at stage time
        if( !_mesh_sta.eval( _istg, tk[_istg-1], _vec_sta ) )
          { _END_ADJ(); return FATAL; }
        switch( options.WRAPMIT){
        case Options::NONE:
          _vec2PMI( _vec_sta, _PMenv, _nx, _PMz, true );
          _vec2PMI( _vec_adj+_pos_adj, _PMenv, _nx, _PMy, true );
          break;
        case Options::DINEQ:
          _vec2PMI( _vec_sta, _PMenv, _nx, _PMz );
          _vec2PMI( _vec_adj+_pos_adj, _PMenv, _nx, _PMy );
          break;
        case Options::ELLIPS:
        default:
          _vec2PME( _vec_sta, _PMenv, _nx, _PMz, _Q, _Er, _Ir );
          _vec2PME( _vec_adj+_pos_adj, _PMenv, _nx, _PMy, _Qy, _Edy, _Idy );
          break;
        }
        _vec2PMI( _vec_adj+_pos_adj+_offset_quad, _PMenv, _np, _PMyq, true );
        if( PMlk && !PMlk[_istg-1] ) PMlk[_istg-1] = new PVT[_nx*_nf];
        for( unsigned iy=0; PMlk[_istg-1] && iy<_nx; iy++)
          PMlk[_istg-1][_ifct*_nx+iy] =_PMy[iy];

        // Add function contribution to adjoint bounds (discontinuities)
        if( _istg>1  ){
          _pos_fct = ( _vFCT.size()>=ns? _istg-1:0 );
          if( _pos_fct
           && !_CC_PM_ADJ( options, -_t, _vec_adj+_pos_adj, _pos_fct, _ifct )
           && !_CC_PM_QUAD( options, -_t, _vec_adj+_pos_adj+_offset_quad ) )
            { _END_ADJ(); return FATAL; }
          // Reset ODE solver - needed in case of discontinuity (used to be in _CC_PM_ADJ)
          gsl_odeiv2_driver_reset( _driver_adj[_ifct] );
        }

        // Add initial state contribution to derivative bounds
        else if( PMdf ){
          if( !_IC_PM_ADJ( -_t ) ){ _END_ADJ(); return FATAL; }
          for( unsigned iq=0; iq<_np; iq++ )
            PMdf[_ifct*_np+iq] = _PMyq[iq];
        }
      }

      // Display & record adjoint intermediate results
      if( options.DISPLAY >= 1 ){
        _print_interm( tk[_istg-1], _nf*_nx, PMlk[_istg-1], "l", os );
        _print_interm( _np, _PMyq, "q", os );
      }
      if( options.RESRECORD )
        _results_adj.push_back( Results( tk[_istg-1], _nf*_nx, PMlk[_istg-1] ) );
    }

    if( options.DISPLAY >= 1 )
      _print_interm( _nf*_np, PMdf, "df", os );
  }
  catch(...){
    _END_ADJ();
    if( options.DISPLAY >= 1 ) _print_stats( stats_adj, os );
    return FAILURE;
  }

  _END_ADJ();
  if( options.DISPLAY >= 1 ) {_print_stats( stats_adj, os );}

  return NORMAL;
}

} // end namescape mc

#endif
