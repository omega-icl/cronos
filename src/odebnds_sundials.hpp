// Copyright (C) 2015 Benoit Chachuat & Nikola Peric, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

#ifndef MC__ODEBNDS_SUNDIALS_HPP
#define MC__ODEBNDS_SUNDIALS_HPP

#undef  MC__ODEBNDS_SUNDIALS_DINEQI_DEBUG
#undef  MC__ODEBNDS_SUNDIALS_DINEQPM_DEBUG

#include "odebnds_base.hpp"
#include "odebnd_sundials.hpp"

#define MC__ODEBNDS_GSL_USE_BAD

//! @brief C++ class computing enclosures of the reachable set of parametric ODEs and adjoint sensitivity analysis using SUNDIALS and MC++.
////////////////////////////////////////////////////////////////////////
//! mc::ODEBNDS_SUNDIALS is a C++ class that computes enclosures of the
//! reachable set of parametric ordinary differential equations (ODEs)
//! using SUNDIALS and MC++. It implements the methods of differential
//! inequalities, whereby polynomial models with interval or
//! ellipsoidalremainders are used to enable high-order convergence. 
//! In addition,the class computes enclosures on the adjoint system of 
//! equations, which can be used to calculate the derivatives of 
//! state-dependent functionals. The use of ellipsoidal remainders 
//! enables stability of the enclosures for asymptotically stable ODE 
//! systems when the parameter host is sufficiently small.
////////////////////////////////////////////////////////////////////////

namespace mc
{

template < typename T, typename PMT=mc::TModel<T>, typename PVT=mc::TVar<T> >
class ODEBNDS_SUNDIALS: public virtual ODEBNDS_BASE<T,PMT,PVT>, 
                        public virtual ODEBND_SUNDIALS<T,PMT,PVT>, public virtual BASE_SUNDIALS
{
private:
  typedef Ellipsoid E;
  typedef BASE_DE::STATUS STATUS;
  typedef int (*CVADJRhsFn)( realtype t, N_Vector x, N_Vector y, N_Vector ydot, void *user_data );

  using ODEBND_BASE<T,PMT,PVT>::NORMAL; 
  using ODEBND_BASE<T,PMT,PVT>::FAILURE;
  using ODEBND_BASE<T,PMT,PVT>::FATAL;
  using ODEBND_BASE<T,PMT,PVT>::_nx;
  using ODEBND_BASE<T,PMT,PVT>::_np;
  using ODEBND_BASE<T,PMT,PVT>::_nf;
  using ODEBND_BASE<T,PMT,PVT>::_vRHS;
  using ODEBND_BASE<T,PMT,PVT>::_vQUAD;
  using ODEBND_BASE<T,PMT,PVT>::_vIC;
  using ODEBND_BASE<T,PMT,PVT>::_vFCT;
  using ODEBND_BASE<T,PMT,PVT>::_t;
  using ODEBND_BASE<T,PMT,PVT>::_istg;
  using ODEBND_BASE<T,PMT,PVT>::_Q;
  using ODEBND_BASE<T,PMT,PVT>::_Er;
  using ODEBND_BASE<T,PMT,PVT>::_Ir;
  using ODEBND_BASE<T,PMT,PVT>::_pref;
  using ODEBND_BASE<T,PMT,PVT>::_Ip;
  using ODEBND_BASE<T,PMT,PVT>::_B;
  using ODEBND_BASE<T,PMT,PVT>::_diam;
  using ODEBND_BASE<T,PMT,PVT>::_vec2I;
  using ODEBND_BASE<T,PMT,PVT>::_vec2E;
  using ODEBND_BASE<T,PMT,PVT>::_PMenv;
  using ODEBND_BASE<T,PMT,PVT>::_PMq;
  using ODEBND_BASE<T,PMT,PVT>::_vec2PMI;
  using ODEBND_BASE<T,PMT,PVT>::_vec2PME;
  using ODEBND_BASE<T,PMT,PVT>::_print_interm;

  using ODEBNDS_BASE<T,PMT,PVT>::_nz;
  using ODEBNDS_BASE<T,PMT,PVT>::_Iz;
  using ODEBNDS_BASE<T,PMT,PVT>::_zref;
  using ODEBNDS_BASE<T,PMT,PVT>::_Iy;
  using ODEBNDS_BASE<T,PMT,PVT>::_TC_I_ADJ;
  using ODEBNDS_BASE<T,PMT,PVT>::_TC_I_QUAD;
  using ODEBNDS_BASE<T,PMT,PVT>::_CC_I_ADJ;
  using ODEBNDS_BASE<T,PMT,PVT>::_SET_I_ADJ;
  using ODEBNDS_BASE<T,PMT,PVT>::_RHS_I_ADJ;
  using ODEBNDS_BASE<T,PMT,PVT>::_RHS_I_QUAD;

  using ODEBND_SUNDIALS<T,PMT,PVT>::_cv_mem;
  using ODEBND_SUNDIALS<T,PMT,PVT>::_cv_flag;
  using ODEBND_SUNDIALS<T,PMT,PVT>::_Nx;
  using ODEBND_SUNDIALS<T,PMT,PVT>::_Nq;
  using ODEBND_SUNDIALS<T,PMT,PVT>::_pos_rhs;
  using ODEBND_SUNDIALS<T,PMT,PVT>::_pos_quad;
  using ODEBND_SUNDIALS<T,PMT,PVT>::_pos_fct;
  using ODEBND_SUNDIALS<T,PMT,PVT>::_init_stats;
  using ODEBND_SUNDIALS<T,PMT,PVT>::_final_stats;
  using ODEBND_SUNDIALS<T,PMT,PVT>::_print_stats;

protected:
  //! @brief current function
  unsigned _ifct;

  //! @brief N_Vector object holding current adjoints
  N_Vector _Ny;
  
  //! @brief CVodeS index
  int indexB;

  //! @brief Size of adjoint bounding system
  unsigned int cv_Ny_size;

  //! @brief Size of quadrature bounding system
  unsigned int cv_Nq_size;

  //! @brief position in _Ny
  unsigned _pos_adj;

  //! @brief position in _Nq
  unsigned _pos_adjquad; 

  //! @brief static pointer to class
  static ODEBNDS_SUNDIALS<T,PMT,PVT> *pODEBNDS;

public:
  typedef BASE_SUNDIALS::Stats Stats;
  typedef typename ODEBND_SUNDIALS<T,PMT,PVT>::Results Results;
  typedef typename ODEBND_SUNDIALS<T,PMT,PVT>::Exceptions Exceptions;

  //! @brief Default constructor
  ODEBNDS_SUNDIALS
    ();

  //! @brief Virtual destructor
  virtual ~ODEBNDS_SUNDIALS
    ();

  //! @brief Integrator options
  struct Options: public ODEBND_SUNDIALS<T,PMT,PVT>::Options
  {
    //! @brief Constructor
    Options():
      ODEBND_SUNDIALS<T,PMT,PVT>::Options()
      {}
    //! @brief Assignment operator
    template <typename U> Options& operator=
      ( U&options ){
        ODEBND_SUNDIALS<T,PMT,PVT>::Options::operator=(options);
        return *this;
      }
  } options;

  //! @brief Statistics for adjoint integration
  Stats stats_adj;

  //! @brief Vector storing interval adjoint bounds (see Options::RESRECORD)
  std::vector< Results > _results_adj;

  //! @brief Propagate state and adjoint interval bounds forward and backward in time through every time stages
  STATUS bounds_ASA
    ( const unsigned ns, const double*tk, const T*Ip, T**Ixk, T*Iq, T*If, T**Ilk, T*Idf,
      std::ostream&os=std::cout );

private:
  //! @brief Function to initialize CVodeS memory block
  bool _INI_CVODES
    ( CVADJRhsFn MC_CVADJRHS, CVADJRhsFn MC_CVADJQUAD );

  //! @brief Function to reinitialize CVodeS memory block
  bool _CC_CVODES
    ();

  //! @brief Function to finalize adjoint bounding
  void _END_ADJ
    ();

  //! @brief Function to initialize adjoint interval bounding
  bool _INI_I_ADJ
    ( const T *Ip );

  //! @brief Static wrapper to function computing the adjoint DINEQ RHS values
  static int MC_CVADJRHSI__
    ( realtype t, N_Vector x, N_Vector y, N_Vector ydot, void *user_data );

  //! @brief Static wrapper to function computing the quadratures RHS in interval arithmetic
  static int MC_CVADJQUADI__
    ( realtype t, N_Vector x, N_Vector y, N_Vector qdot, void *user_data );

  //! @brief Private methods to block default compiler methods
  ODEBNDS_SUNDIALS(const ODEBNDS_SUNDIALS&);
  ODEBNDS_SUNDIALS& operator=(const ODEBNDS_SUNDIALS&);
};

template <typename T, typename PMT, typename PVT>
 ODEBNDS_SUNDIALS<T,PMT,PVT>* ODEBNDS_SUNDIALS<T,PMT,PVT>::pODEBNDS = 0;

template <typename T, typename PMT, typename PVT> inline
ODEBNDS_SUNDIALS<T,PMT,PVT>::ODEBNDS_SUNDIALS
()
: BASE_SUNDIALS(), ODEBNDS_BASE<T,PMT,PVT>(), ODEBND_SUNDIALS<T,PMT,PVT>(), _ifct(0), _Ny(0)
{}

template <typename T, typename PMT, typename PVT> inline
ODEBNDS_SUNDIALS<T,PMT,PVT>::~ODEBNDS_SUNDIALS
()
{
  if( _Ny ) N_VDestroy_Serial( _Ny );
}

template <typename T, typename PMT, typename PVT> inline bool
ODEBNDS_SUNDIALS<T,PMT,PVT>::_INI_CVODES
( CVADJRhsFn MC_CVADJRHS, CVADJRhsFn MC_CVADJQUAD )
{
  switch( options.INTMETH ){
   // Non-stiff case: create CVodeS memory block for the ADAMS method
   // and use of a functional iteration
   case Options::MSADAMS: default:
    _cv_flag = CVodeCreateB( _cv_mem, CV_ADAMS, CV_FUNCTIONAL, &indexB );
    if( _check_cv_flag(&_cv_flag, "CVodeCreateB", 1) ) return false;
    break;

   // Stiff case: create CVodeS memory block for the BDF method
   // and use of a Newton iteration
   case Options::MSBDF:
    _cv_flag = CVodeCreateB( _cv_mem, CV_BDF, CV_NEWTON, &indexB );
    if( _check_cv_flag(&_cv_flag, "CVodeCreateB", 1) ) return false;
    break;
  }

  // Initialize CVodeS memory and specify right hand side function,
  // current time _t, and current adjoint _Ny
  _cv_flag = CVodeInitB( _cv_mem, indexB, MC_CVADJRHS, _t, _Ny );
  if( _check_cv_flag(&_cv_flag, "CVodeInitB", 1) ) return false;

  // Specify the Jacobian approximation and linear solver
  if( options.INTMETH == Options::MSBDF ){
    switch( options.JACAPPROX ){
     case Options::CV_DIAG: default:
       _cv_flag = CVDiagB( _cv_mem, indexB );
       if( _check_cv_flag(&_cv_flag, "CVDiagB", 1)) return false;
       break;

     case Options::CV_LAPACKDENSE:
       //_cv_flag = CVLapackDenseB( _cv_mem, indexB, NV_LENGTH_S( _Ny ) );
       //if( _check_cv_flag(&_cv_flag, "CVLapackDenseB", 1)) return false;
       //break;

     case Options::CV_DENSE:
       _cv_flag = CVDenseB( _cv_mem, indexB, NV_LENGTH_S( _Ny ) );
       if( _check_cv_flag(&_cv_flag, "CVDenseB", 1)) return false;
       break;
    }
    _cv_flag = CVDlsSetDenseJacFnB( _cv_mem, indexB, NULL );
    if ( _check_cv_flag(&_cv_flag, "CVDlsSetDenseJacFnB", 1) ) return false;
  }

  // Specify the relative and absolute tolerances for states
  _cv_flag = CVodeSStolerancesB( _cv_mem, indexB, options.RTOL, options.ATOL );
  if( _check_cv_flag(&_cv_flag, "CVodeSStolerancesB", 1) ) return false;

  // Set maximum number of error test failures
  //_cv_flag = CVodeSetMaxErrTestFailsB( _cv_mem, indexB, options.MAXFAIL );
  //if ( _check_cv_flag(&_cv_flag, "CVodeSetMaxErrTestFailsB", 1) ) return false;

  // Specify minimum stepsize
  _cv_flag = CVodeSetMinStepB( _cv_mem, indexB, options.HMIN>0.? options.HMIN:0. );
  if( _check_cv_flag(&_cv_flag, "CVodeSetMinStepB", 1) ) return false;

  // Specify maximum stepsize
  _cv_flag = CVodeSetMaxStepB( _cv_mem, indexB, options.HMAX>0.? options.HMAX: 0. );
  if( _check_cv_flag(&_cv_flag, "CVodeSetMaxStepB", 1) ) return false;

  // Specify maximum number of steps between two stage times
  _cv_flag = CVodeSetMaxNumStepsB( _cv_mem, indexB, options.NMAX );
  if( _check_cv_flag(&_cv_flag, "CVodeSetMaxNumStepsB", 1) ) return false;

  // Initialize the integrator memory for the quadrature variables
  if( !_Nq ) return true;
  _cv_flag = CVodeQuadInitB( _cv_mem, indexB, MC_CVADJQUAD, _Nq );
  if( _check_cv_flag(&_cv_flag, "CVodeQuadInitB", 1) ) return false;

  // Specify whether or not to perform error control on quadrature
  _cv_flag = CVodeSetQuadErrConB( _cv_mem, indexB, options.QERR );
  if( _check_cv_flag(&_cv_flag, "CVodeSetQuadErrConB", 1) ) return false;

  // Specify the relative and absolute tolerances for quadratures
  _cv_flag = CVodeQuadSStolerancesB( _cv_mem, indexB, options.RTOL, options.ATOL );
  if( _check_cv_flag(&_cv_flag, "CVodeQuadSStolerancesB", 1) ) return false;

  return true;
}

template <typename T, typename PMT, typename PVT> inline bool
ODEBNDS_SUNDIALS<T,PMT,PVT>::_CC_CVODES
()
{
  // Reinitialize CVodeS memory block for current time _t and
  // current adjoint _Ny
  _cv_flag = CVodeReInitB( _cv_mem, indexB, _t, _Ny );
  if( _check_cv_flag(&_cv_flag, "CVodeReInitB", 1) ) return false;

  // Reinitialize CVodeS memory block for current quarature _Nq
  if( !_Nq ) return true;
  _cv_flag = CVodeQuadReInitB( _cv_mem, indexB, _Nq );
  if( _check_cv_flag(&_cv_flag, "CVodeQuadReInitB", 1) ) return false;

  return true;
}

template <typename T, typename PMT, typename PVT> inline void
ODEBNDS_SUNDIALS<T,PMT,PVT>::_END_ADJ()
{
  // Get final CPU time
  _final_stats( stats_adj );
}

template <typename T, typename PMT, typename PVT> inline bool
ODEBNDS_SUNDIALS<T,PMT,PVT>::_INI_I_ADJ
( const T* Ip )
{
  // Initialize bound propagation
  if( !ODEBNDS_BASE<T,PMT,PVT>::_INI_I_ADJ( options, Ip ) ) return false;

  // Set SUNDIALS adjoint/quadrature arrays
  cv_Nq_size = 2*_np;
  _nz = _nx + _np;
  switch( options.WRAPMIT){
  case Options::NONE:
  case Options::DINEQ:
    cv_Ny_size = 2*_nx;
    break;
  case Options::ELLIPS:
  default:
    cv_Ny_size = _nx*(1+_nz)+_nx*(_nx+1)/2;
    break;
  }
  if( !_Ny || cv_Ny_size != NV_LENGTH_S( _Ny ) ){
    if( _Ny ) N_VDestroy_Serial( _Ny );
    _Ny  = N_VNew_Serial( cv_Ny_size*_nf );
  }
  if( !_Nq || cv_Nq_size != NV_LENGTH_S( _Nq ) ){
    if( _Nq ) N_VDestroy_Serial( _Nq );
    _Nq  = N_VNew_Serial( cv_Nq_size*_nf );
  }

  // Reset result record and statistics
  _results_adj.clear();
  _init_stats( stats_adj );

  return true;
}

template <typename T, typename PMT, typename PVT> inline int
ODEBNDS_SUNDIALS<T,PMT,PVT>::MC_CVADJRHSI__
( realtype t, N_Vector x, N_Vector y, N_Vector ydot, void *user_data )
{
  ODEBNDS_SUNDIALS<T,PMT,PVT> *pODEBNDS = ODEBNDS_SUNDIALS<T,PMT,PVT>::pODEBNDS;
  bool flag = pODEBNDS->_RHS_I_ADJ( pODEBNDS->options, t, NV_DATA_S( y ), NV_DATA_S( ydot ),
                                    NV_DATA_S( x ) );
  ODEBNDS_SUNDIALS<T,PMT,PVT>::pODEBNDS = pODEBNDS;
  pODEBNDS->stats_adj.numRHS++;
  return( flag? 0: -1 );
}

template <typename T, typename PMT, typename PVT> inline int
ODEBNDS_SUNDIALS<T,PMT,PVT>::MC_CVADJQUADI__
( realtype t, N_Vector x, N_Vector y, N_Vector qdot, void *user_data )
{
  ODEBNDS_SUNDIALS<T,PMT,PVT> *pODEBNDS = ODEBNDS_SUNDIALS<T,PMT,PVT>::pODEBNDS;
  bool flag = pODEBNDS->_RHS_I_QUAD( pODEBNDS->options, t, NV_DATA_S( y ), NV_DATA_S( qdot ) );
  ODEBNDS_SUNDIALS<T,PMT,PVT>::pODEBNDS = pODEBNDS;
  return( flag? 0: -1 );
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
inline typename ODEBNDS_SUNDIALS<T,PMT,PVT>::STATUS
ODEBNDS_SUNDIALS<T,PMT,PVT>::bounds_ASA
( const unsigned ns, const double*tk, const T*Ip, T**Ixk, T*Iq, T*If, T**Ilk, T*Idf,
  std::ostream&os )
{
  // Compute state bounds and store intermediate results
  //std::cout << "Before state bounding:\n" << ODEBND_BASE<T,PMT,PVT>::_IVAR[0] << std::endl;
  STATUS flag = NORMAL;
  ODEBND_SUNDIALS<T,PMT,PVT>::options = options;
  flag = ODEBND_SUNDIALS<T,PMT,PVT>::bounds( ns, tk, Ip, Ixk, Iq, If, true, os);
  if( flag != NORMAL ) return flag;

  // Nothing to do if no functions are defined
  if( !_nf ) return NORMAL;

  // Check size
  if( !Ilk || !Idf ) return FATAL;

  try{
    // Initialize adjoint bound integration
    if( !_INI_I_ADJ( Ip )) return FATAL;

    // No explicit interpolation...

    // Bounds on terminal states/quadratures
    switch( options.WRAPMIT){
    case Options::NONE:
    case Options::DINEQ:
      _vec2I( NV_DATA_S( _Nx ), _nx, _Iz );
      break;
    case Options::ELLIPS:
    default:
      _vec2E( NV_DATA_S( _Nx ), _nx, _np, _Q, _Er, _Ir, _pref, _Ip, _B, _zref, _Iz );
      break;
    }

    // Set time
    _t = tk[ns];

    // Bounds on terminal adjoints/quadratures
    if( Ilk && !Ilk[ns] ) Ilk[ns] = new T[_nx*_nf];
    for( _ifct=_pos_adj=_pos_adjquad=0; _ifct < _nf; _ifct++, 
         _pos_adj+=cv_Ny_size, _pos_adjquad+=cv_Nq_size  ){
      _pos_fct = ( _vFCT.size()>=ns? ns-1:0 );
      if( !_TC_I_ADJ( options, _t, NV_DATA_S( _Ny )+_pos_adj, _pos_fct, _ifct)
       || ( _Nq && !_TC_I_QUAD( options, NV_DATA_S( _Nq )+_pos_adjquad ) ) )
        { _END_ADJ(); return FATAL; }
      for( unsigned iy=0; Ilk[ns] && iy<_nx; iy++ ) {Ilk[ns][_ifct*_nx+iy] = _Iy[iy];}
    }
    // Display & record adjoint terminal results
    if( options.DISPLAY >= 1 )
      _print_interm( tk[ns], _nx*_nf, Ilk[ns], "l", os );
    if( options.RESRECORD )
      _results_adj.push_back( Results( tk[ns], _nf*_nx, Ilk[ns] ) );

    // Integrate adjoint ODEs through each stage using SUNDIALS
    pODEBNDS = this;
    if( !_INI_CVODES( MC_CVADJRHSI__, MC_CVADJQUADI__ ) ) { _END_ADJ(); return FATAL;}

    for( _istg=ns; _istg>0; _istg-- ){
      // Integrate backward through current stage for each function
      if( Ilk && !Ilk[_istg-1] ) Ilk[_istg-1] = new T[_nx*_nf];

      for( _ifct=_pos_adj=_pos_adjquad=0; _ifct < _nf; _ifct++, 
           _pos_adj+=cv_Ny_size, _pos_adjquad+=cv_Nq_size ){

        _t = tk[_istg];

        // Update list of operations in RHSADJ and QUADADJ
        _pos_rhs  = ( _vRHS.size() <=1? 0: _istg-1 );
        _pos_quad = ( _vQUAD.size()<=1? 0: _istg-1 );
        if( !_SET_I_ADJ( options, _pos_rhs, _pos_quad, _pos_fct, _ifct ) )
          { _END_ADJ(); return FATAL; }

        // Propagate bounds backward to previous stage time
        //_cv_flag = CVodeB( _cv_mem, tk[_istg-1], CV_NORMAL );
        //if( _check_cv_flag(&_cv_flag, "CVodeB", 1) ) { _END_ADJ(); return FATAL; }

      }

    }





  }
  catch(...){
    _END_ADJ();
    if( options.DISPLAY >= 1 ) _print_stats( stats_adj, os );
    return FAILURE;
  }

  _END_ADJ();
  if( options.DISPLAY >= 1 ) {_print_stats( stats_adj, os );}
  std::cout << "At the end of bounds" << std::endl;
  return NORMAL;
}

} // end namescape mc

#endif




















