// Copyright (C) 2015 Benoit Chachuat & Nikola Peric, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

#ifndef MC__ODEBNDS_SUNDIALS_HPP
#define MC__ODEBNDS_SUNDIALS_HPP

#undef  MC__ODEBNDS_SUNDIALS_DINEQI_DEBUG
#undef  MC__ODEBNDS_SUNDIALS_DINEQPM_DEBUG

#include "odebnds_base.hpp"
#include "odebnd_sundials.hpp"

#define MC__ODEBNDS_SUNDIALS_USE_BAD

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
class ODEBNDS_SUNDIALS:
  public virtual BASE_DE,
  public virtual BASE_SUNDIALS,
  public virtual ODEBNDS_BASE<T,PMT,PVT>, 
  public virtual ODEBND_SUNDIALS<T,PMT,PVT>
{
 protected:
  typedef Ellipsoid E;
  typedef BASE_DE::STATUS STATUS;
  typedef int (*CVRhsFn)( realtype t, N_Vector y, N_Vector ydot, void *user_data );
  typedef int (*CVADJRhsFn)( realtype t, N_Vector x, N_Vector y, N_Vector ydot, void *user_data );

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
  using ODEBND_BASE<T,PMT,PVT>::_PMp;
  using ODEBND_BASE<T,PMT,PVT>::_vec2PMI;
  using ODEBND_BASE<T,PMT,PVT>::_vec2PME;
  using ODEBND_BASE<T,PMT,PVT>::_print_interm;

  using ODEBNDS_BASE<T,PMT,PVT>::_Ix;
  using ODEBNDS_BASE<T,PMT,PVT>::_Iy;
  using ODEBNDS_BASE<T,PMT,PVT>::_Iyq;
  using ODEBNDS_BASE<T,PMT,PVT>::_Idy;
  using ODEBNDS_BASE<T,PMT,PVT>::_Idyq;
  using ODEBNDS_BASE<T,PMT,PVT>::_Qy;
  using ODEBNDS_BASE<T,PMT,PVT>::_By;
  using ODEBNDS_BASE<T,PMT,PVT>::_Byq;
  using ODEBNDS_BASE<T,PMT,PVT>::_Edy;
  using ODEBNDS_BASE<T,PMT,PVT>::_yref;
  using ODEBNDS_BASE<T,PMT,PVT>::_PMx;
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

  using ODEBND_SUNDIALS<T,PMT,PVT>::_cv_mem;
  using ODEBND_SUNDIALS<T,PMT,PVT>::_cv_flag;
  using ODEBND_SUNDIALS<T,PMT,PVT>::_Nx;
  using ODEBND_SUNDIALS<T,PMT,PVT>::_vec_sta;
  using ODEBND_SUNDIALS<T,PMT,PVT>::_pos_rhs;
  using ODEBND_SUNDIALS<T,PMT,PVT>::_pos_quad;
  using ODEBND_SUNDIALS<T,PMT,PVT>::_pos_fct;
  using ODEBND_SUNDIALS<T,PMT,PVT>::_init_stats;
  using ODEBND_SUNDIALS<T,PMT,PVT>::_final_stats;
  using ODEBND_SUNDIALS<T,PMT,PVT>::_print_stats;

 protected:
  //! @brief current function
  unsigned _ifct;

  //! @brief number of function in N_vector arrays
  unsigned _nfct;

  //! @brief N_Vector object holding current adjoints
  N_Vector *_Ny;

  //! @brief N_Vector object holding current quadratures
  N_Vector *_Nyq;
  
  //! @brief pointer to array holding identifiers of the backward problems
  int* _indexB;
  
  //! @brief pointer to array holding identifiers of the backward problems
  unsigned* _iusrB;

  //! @brief position in _Ny
  unsigned _pos_adj;

  //! @brief position in _Nq
  unsigned _pos_adjquad; 

  //! @brief static pointer to class
  static ODEBNDS_SUNDIALS<T,PMT,PVT> *pODEBNDS;

 public:
  typedef BASE_SUNDIALS::Stats Stats;
  typedef typename ODEBND_BASE<T,PMT,PVT>::Results Results;
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

  //! @brief Propagate state/quadrature interval bounds forward in time through every time stages
  STATUS bounds
    ( const unsigned ns, const double*tk, const T*Ip, T**Ixk,
      T*Iq, T*If, std::ostream&os=std::cout )
      { ODEBND_SUNDIALS<T,PMT,PVT>::options = options;
        return ODEBND_SUNDIALS<T,PMT,PVT>::_bounds( ns, tk, Ip, Ixk, Iq, If, false, os); }

  //! @brief Propagate state/quadrature polynomial models forward in time through every time stages
  STATUS bounds
    ( const unsigned ns, const double*tk, const PVT*PMp, PVT**PMxk,
      PVT*PMq, PVT*PMf, std::ostream&os=std::cout )
      { ODEBND_SUNDIALS<T,PMT,PVT>::options = options;
        return ODEBND_SUNDIALS<T,PMT,PVT>::_bounds( ns, tk, PMp, PMxk, PMq, PMf, false, os); }

  //! @brief Propagate state and adjoint interval bounds forward and backward in time through every time stages
  STATUS bounds_ASA
    ( const unsigned ns, const double*tk, const T*Ip, T**Ixk, T*Iq, T*If, T**Ilk, T*Idf,
      std::ostream&os=std::cout );

  //! @brief Propagate state and adjoint polynomial models forward and backward in time through every time stages
  STATUS bounds_ASA
    ( const unsigned ns, const double*tk, const PVT*PMp, PVT**PMxk, PVT*PMq, PVT*PMf,
      PVT**PMlk, PVT*PMdf, std::ostream&os=std::cout );

  //! @brief Record state and sensitivity bounds in files <a>obndsta</a> and <a>obndsa</a>, with accuracy of <a>iprec</a> digits
  void record
    ( std::ofstream&obndsta, std::ofstream&obndsa, const unsigned iprec=5 ) const
    { this->ODEBND_SUNDIALS<T,PMT,PVT>::record( obndsta, iprec );
      this->ODEBND_BASE<T,PMT,PVT>::_record( obndsa, _results_adj, iprec ); }

  //! @brief Record state and sensitivity bounds in files <a>obndsta</a> and <a>obndsa</a>, with accuracy of <a>iprec</a> digits
  void record
    ( std::ofstream&obndsta, const unsigned iprec=5 ) const
    { this->ODEBND_SUNDIALS<T,PMT,PVT>::record( obndsta, iprec ); }

 private:
  //! @brief Function to initialize CVode memory block (virtual)
  virtual bool _INI_CVODE
    ( CVRhsFn MC_CVRHS, CVRhsFn MC_CVQUAD );

  //! @brief Function to initialize CVodeS memory block
  bool _INI_CVODES
    ( CVADJRhsFn MC_CVADJRHS, CVADJRhsFn MC_CVADJQUAD, const unsigned ifct,
      int&indexB, unsigned&iusrB );

  //! @brief Function to reinitialize CVodeS memory block
  bool _CC_CVODES
    ( const unsigned ifct, const int indexB );

  //! @brief Function to finalize adjoint bounding
  void _END_ADJ
    ();

  //! @brief Function to initialize adjoint interval bounding
  bool _INI_I_ADJ
    ( const unsigned np, const T *Ip );

  //! @brief Static wrapper to function computing the adjoint DINEQ RHS values
  static int MC_CVADJRHSI__
    ( realtype t, N_Vector x, N_Vector y, N_Vector ydot, void *user_data );

  //! @brief Static wrapper to function computing the quadratures RHS in interval arithmetic
  static int MC_CVADJQUADI__
    ( realtype t, N_Vector x, N_Vector y, N_Vector qdot, void *user_data );

  //! @brief Function to initialize adjoint polynomial models
  bool _INI_PM_ADJ
    ( const unsigned np, const PVT *PMp );

  //! @brief Static wrapper to function to calculate the adjoint DINEQ-PM RHS values
  static int MC_CVADJRHSPM__
    ( realtype t, N_Vector x, N_Vector y, N_Vector ydot, void *user_data );

  //! @brief Static wrapper to function computing the PM-quadratures RHS values
  static int MC_CVADJQUADPM__
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
: BASE_SUNDIALS(), ODEBNDS_BASE<T,PMT,PVT>(), ODEBND_SUNDIALS<T,PMT,PVT>(),
  _ifct(0), _nfct(0), _Ny(0), _Nyq(0), _indexB(0), _iusrB(0)
{}

template <typename T, typename PMT, typename PVT> inline
ODEBNDS_SUNDIALS<T,PMT,PVT>::~ODEBNDS_SUNDIALS
()
{
  if( _Ny )   N_VDestroyVectorArray_Serial( _Ny,  _nfct );
  if( _Nyq )  N_VDestroyVectorArray_Serial( _Nyq, _nfct );
  delete[] _indexB;
  delete[] _iusrB;
}

template <typename T, typename PMT, typename PVT> inline bool
ODEBNDS_SUNDIALS<T,PMT,PVT>::_INI_CVODE
( CVRhsFn MC_CVRHS, CVRhsFn MC_CVQUAD )
{
  // Call _INI_CVODE in ODEBND_SUNDIALS
  this->ODEBND_SUNDIALS<T,PMT,PVT>::_INI_CVODE( MC_CVRHS, MC_CVQUAD );
  
  // Allocate memory for adjoint integration
  _cv_flag = CVodeAdjInit( _cv_mem, options.ASACHKPT, options.ASAINTERP );
  if( _check_cv_flag(&_cv_flag, "CVodeAdjInit", 1) ) return false;

  // Reinitialize adjoint holding vectors
  if( _nfct != _nf ){
    if( _Ny )   N_VDestroyVectorArray_Serial( _Ny,  _nfct );
    if( _Nyq )  N_VDestroyVectorArray_Serial( _Nyq, _nfct );
    delete[] _indexB;
    delete[] _iusrB;
    _nfct = _nf;
    _Ny  = N_VCloneVectorArray_Serial( _nfct, _Nx );
    _Nyq = N_VCloneVectorArray_Serial( _nfct, _Nx );
    _indexB = new int[_nfct];
    _iusrB = new unsigned[_nfct];
  }
  
  return true;
}

template <typename T, typename PMT, typename PVT> inline bool
ODEBNDS_SUNDIALS<T,PMT,PVT>::_INI_CVODES
( CVADJRhsFn MC_CVADJRHS, CVADJRhsFn MC_CVADJQUAD, const unsigned ifct,
  int&indexB, unsigned&iusrB )
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

  // Initialize CVodeS memory and specify the adjoint RHS function,
  // terminal time _t, and terminal adjoint _Ny
  _cv_flag = CVodeInitB( _cv_mem, indexB, MC_CVADJRHS, _t, _Ny[ifct] );
  if( _check_cv_flag(&_cv_flag, "CVodeInitB", 1) ) return false;

  // Specify the user_data to pass the function index corresponding to
  // the RHS functuinrelative and absolute tolerancesor states
  iusrB = _ifct;
  _cv_flag = CVodeSetUserDataB( _cv_mem, indexB, &iusrB );
  if( _check_cv_flag(&_cv_flag, "CVodeSetUserDataB", 1) ) return false;

  // Specify the Jacobian approximation and linear solver
  if( options.INTMETH == Options::MSBDF ){
    switch( options.JACAPPROX ){
     case Options::CV_DIAG: default:
       _cv_flag = CVDiagB( _cv_mem, indexB );
       if( _check_cv_flag(&_cv_flag, "CVDiagB", 1)) return false;
       break;

     case Options::CV_LAPACKDENSE:
       //_cv_flag = CVLapackDenseB( _cv_mem, indexB, NV_LENGTH_S( _Ny[ifct] ) );
       //if( _check_cv_flag(&_cv_flag, "CVLapackDenseB", 1)) return false;
       //break;

     case Options::CV_DENSE:
       _cv_flag = CVDenseB( _cv_mem, indexB, NV_LENGTH_S( _Ny[ifct] ) );
       if( _check_cv_flag(&_cv_flag, "CVDenseB", 1)) return false;
       break;
    }
    _cv_flag = CVDlsSetDenseJacFnB( _cv_mem, indexB, NULL );
    if ( _check_cv_flag(&_cv_flag, "CVDlsSetDenseJacFnB", 1) ) return false;
  }

  // Specify the relative and absolute tolerances for states
  _cv_flag = CVodeSStolerancesB( _cv_mem, indexB, options.RTOLB, options.ATOLB );
  if( _check_cv_flag(&_cv_flag, "CVodeSStolerancesB", 1) ) return false;

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
  if( !_Nyq ) return true;
  _cv_flag = CVodeQuadInitB( _cv_mem, indexB, MC_CVADJQUAD, _Nyq[ifct] );
  if( _check_cv_flag(&_cv_flag, "CVodeQuadInitB", 1) ) return false;

  // Specify whether or not to perform error control on quadrature
  _cv_flag = CVodeSetQuadErrConB( _cv_mem, indexB, options.QERRB );
  if( _check_cv_flag(&_cv_flag, "CVodeSetQuadErrConB", 1) ) return false;

  // Specify the relative and absolute tolerances for quadratures
  _cv_flag = CVodeQuadSStolerancesB( _cv_mem, indexB, options.RTOLB, options.ATOLB );
  if( _check_cv_flag(&_cv_flag, "CVodeQuadSStolerancesB", 1) ) return false;

  return true;
}

template <typename T, typename PMT, typename PVT> inline bool
ODEBNDS_SUNDIALS<T,PMT,PVT>::_CC_CVODES
( const unsigned ifct, const int indexB )
{
  // Reinitialize CVodeS memory block for current time _t and
  // current adjoint _Ny
  _cv_flag = CVodeReInitB( _cv_mem, indexB, _t, _Ny[ifct] );
  if( _check_cv_flag(&_cv_flag, "CVodeReInitB", 1) ) return false;

  // Reinitialize CVodeS memory block for current quarature _Nq
  if( !_Nyq ) return true;
  _cv_flag = CVodeQuadReInitB( _cv_mem, indexB, _Nyq[ifct] );
  if( _check_cv_flag(&_cv_flag, "CVodeQuadReInitB", 1) ) return false;

  return true;
}

template <typename T, typename PMT, typename PVT> inline bool
ODEBNDS_SUNDIALS<T,PMT,PVT>::_INI_I_ADJ
( const unsigned np, const T*Ip )
{
  // Initialize bound propagation
  if( !ODEBNDS_BASE<T,PMT,PVT>::_INI_I_ADJ( options, np, Ip ) )
    return false;

  // Set SUNDIALS adjoint/quadrature arrays
  unsigned cv_Ny_size, cv_Nyq_size;
  switch( options.WRAPMIT){
  case Options::NONE:
  case Options::DINEQ:
    cv_Ny_size  = 2*_nx;
    cv_Nyq_size = 2*np;
    break;
  case Options::ELLIPS:
  default:
    cv_Ny_size  = _nx*(1+np)+_nx*(_nx+1)/2;
    cv_Nyq_size = _np*(2+np);
    break;
  }
  for( unsigned i=0; i<_nf; i++ ){
    if( !_Ny[i] || cv_Ny_size != NV_LENGTH_S( _Ny[i] ) ){
      if( _Ny[i] ) N_VDestroy_Serial( _Ny[i] );
      _Ny[i] = N_VNew_Serial( cv_Ny_size );
    }
    if( !_Nyq[i] || cv_Nyq_size != NV_LENGTH_S( _Nyq[i] ) ){
      if( _Nyq[i] ) N_VDestroy_Serial( _Nyq[i] );
      _Nyq[i] = N_VNew_Serial( cv_Nyq_size );
    }
  }

  // Reset result record and statistics
  _results_adj.clear();
  _init_stats( stats_adj );

  return true;
}

template <typename T, typename PMT, typename PVT> inline void
ODEBNDS_SUNDIALS<T,PMT,PVT>::_END_ADJ()
{
  // Get final CPU time
  _final_stats( stats_adj );
}

template <typename T, typename PMT, typename PVT> inline int
ODEBNDS_SUNDIALS<T,PMT,PVT>::MC_CVADJRHSI__
( realtype t, N_Vector x, N_Vector y, N_Vector ydot, void *user_data )
{
  ODEBNDS_SUNDIALS<T,PMT,PVT> *pODEBNDS = ODEBNDS_SUNDIALS<T,PMT,PVT>::pODEBNDS;
  pODEBNDS->_ifct = *static_cast<unsigned*>( user_data );
  bool flag = pODEBNDS->_RHS_I_ADJ( pODEBNDS->options, t, NV_DATA_S( y ),
    NV_DATA_S( ydot ), NV_DATA_S( x ), pODEBNDS->_ifct, true );
  ODEBNDS_SUNDIALS<T,PMT,PVT>::pODEBNDS = pODEBNDS;
  pODEBNDS->stats_adj.numRHS++;
  return( (flag && _diam( pODEBNDS->_nx, pODEBNDS->_Iy ) < pODEBNDS->options.DMAX)? 0: -1 );
}

template <typename T, typename PMT, typename PVT> inline int
ODEBNDS_SUNDIALS<T,PMT,PVT>::MC_CVADJQUADI__
( realtype t, N_Vector x, N_Vector y, N_Vector qdot, void *user_data )
{
  ODEBNDS_SUNDIALS<T,PMT,PVT> *pODEBNDS = ODEBNDS_SUNDIALS<T,PMT,PVT>::pODEBNDS;
  pODEBNDS->_ifct = *static_cast<unsigned*>( user_data );
  bool flag = pODEBNDS->_RHS_I_QUAD( pODEBNDS->options, t, NV_DATA_S( y ),
    NV_DATA_S( qdot ), NV_DATA_S( x ), pODEBNDS->_ifct, true );
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
( const unsigned ns, const double*tk, const T*Ip, T**Ixk, T*Iq, T*If,
  T**Ilk, T*Idf, std::ostream&os )
{
  // Compute state bounds and store intermediate results
  STATUS flag = NORMAL;
  ODEBND_SUNDIALS<T,PMT,PVT>::options = options;
  flag = ODEBND_SUNDIALS<T,PMT,PVT>::_bounds( ns, tk, Ip, Ixk, Iq, If, true, os);
  if( flag != NORMAL ) return flag;

  // Nothing to do if no functions are defined
  if( !_nf ) return NORMAL;

  // Check size
  if( !Ilk || !Idf ) return FATAL;

  try{
    // Initialize adjoint bound integration
    if( !_INI_I_ADJ( _np, Ip )) return FATAL;
    _t = tk[ns];

    // Bounds on terminal states/quadratures
    //for( unsigned i=0; i<_vec_sta[ns].size(); i++ )
    //  NV_DATA_S(_Nx)[i] = _vec_sta[ns][i];
    switch( options.WRAPMIT){
    case Options::NONE:
    case Options::DINEQ:
      _vec2I( _vec_sta[ns].data(), _nx, _Ix );
      break;
    case Options::ELLIPS:
    default:
      _vec2E( _vec_sta[ns].data(), _nx, _np, _Q, _Er, _Ir, _pref, _Ip, _B, _xref, _Ix );
      break;
    }

    // Bounds on terminal adjoints/quadratures
    if( Ilk && !Ilk[ns] ) Ilk[ns] = new T[_nx*_nf];
    for( _ifct=0; _ifct < _nf; _ifct++ ){
      _pos_fct = ( _vFCT.size()>=ns? ns-1:0 );
      if( !_TC_I_ADJ( options, _t, NV_DATA_S(_Ny[_ifct]), _pos_fct, _ifct)
       || ( _Nyq && _Nyq[_ifct] && !_TC_I_QUAD( options, _t, NV_DATA_S(_Nyq[_ifct]) ) ) )
        { _END_ADJ(); return FATAL; }
      for( unsigned iy=0; Ilk[ns] && iy<_nx; iy++ )
        Ilk[ns][_ifct*_nx+iy] = _Iy[iy];
      for( unsigned iq=0; iq<_np; iq++ )
        Idf[_ifct*_np+iq] = _Iyq[iq];
    }

    // Display & record adjoint terminal results
    if( options.DISPLAY >= 1 ){
      _print_interm( tk[ns], _nf*_nx, Ilk[ns], "l", os );//_nf
      //_print_interm( _nf*_np, Idf, "df", os );
    }
    if( options.RESRECORD )
      _results_adj.push_back( Results( tk[ns], _nx*_nf, Ilk[ns] ) );//_nf

    // Initialization of adjoint integration
    for( _ifct=0; _ifct < _nf; _ifct++ )
      if( !_INI_CVODES( MC_CVADJRHSI__, MC_CVADJQUADI__, _ifct,
        _indexB[_ifct], _iusrB[_ifct] ) )
        { _END_ADJ(); return FATAL;}

    // Integrate adjoint ODEs through each stage using SUNDIALS
    pODEBNDS = this;
    for( _istg=ns; _istg>0; _istg-- ){

      // Update list of operations in RHSADJ and QUADADJ
      _pos_rhs  = ( _vRHS.size() <=1? 0: _istg-1 );
      _pos_quad = ( _vQUAD.size()<=1? 0: _istg-1 );
      if( !_SET_I_ADJ( options, _pos_rhs, _pos_quad, _pos_fct ) )
        { _END_ADJ(); return FATAL; }

      // Propagate bounds backward to previous stage time
      //_cv_flag = CVodeSetStopTimeB( _cv_mem, tk[_istg-1] );
      //if( _check_flag(&_cv_flag, "CVodeSetStopTimeB", 1) )
      //  { _END_ADJ(); return FATAL; }
      _cv_flag = CVodeB( _cv_mem, tk[_istg-1], CV_NORMAL );
      if( _check_cv_flag(&_cv_flag, "CVodeB", 1) )
        { _END_ADJ(); return FATAL; }
      _t = tk[_istg-1];

      // Bounds on states/adjoints/quadratures at stage time
      //stats_adj.numSteps = 0; 
      for( _ifct=0; _ifct < _nf; _ifct++ ){
        //for( unsigned i=0; i<_vec_sta[_istg-1].size(); i++ )
        //  NV_DATA_S(_Nx)[i] = _vec_sta[_istg-1][i];
        _cv_flag = CVodeGetB( _cv_mem, _indexB[_ifct], &_t, _Ny[_ifct]);
        if( _check_cv_flag( &_cv_flag, "CVodeGetB", 1) )
          { _END_ADJ(); return FATAL; }
        _cv_flag = CVodeGetQuadB( _cv_mem, _indexB[_ifct], &_t, _Nyq[_ifct]);
        if( _check_cv_flag( &_cv_flag, "CVodeGetQuadB", 1) )
          { _END_ADJ(); return FATAL; }
        void *cv_memB = CVodeGetAdjCVodeBmem(_cv_mem, _indexB[_ifct] );
        long int nstpB;
        _cv_flag = CVodeGetNumSteps( cv_memB, &nstpB );
        stats_adj.numSteps += nstpB;
        switch( options.WRAPMIT){
        case Options::NONE:
        case Options::DINEQ:
          _vec2I( _vec_sta[_istg-1].data(), _nx, _Ix );
          _vec2I( NV_DATA_S(_Ny[_ifct]), _nx, _Iy);
          _vec2I( NV_DATA_S( _Nyq[_ifct] ), _np, _Iyq );
          break;
        case Options::ELLIPS:
        default:
          _vec2E( _vec_sta[_istg-1].data(), _nx, _np, _Q, _Er, _Ir, _pref, _Ip, _B, _xref, _Ix );
          _vec2E( NV_DATA_S(_Ny[_ifct]), _nx, _np, _Qy, _Edy, _Idy, _pref, _Ip, _By, _yref, _Iy );
#ifdef MC__ODEBNDS_GSL_DINEQI_DEBUG
          std::cout << "Edy: " << _Edy;
#endif
          _vec2I( NV_DATA_S( _Nyq[_ifct] ), _np, _np, _pref, _Ip, _Byq, _Idyq, _Iyq );
          break;
        }
        if( Ilk && !Ilk[_istg-1] ) Ilk[_istg-1] = new T[_nx*_nf];
        for( unsigned iy=0; Ilk[_istg-1] && iy<_nx; iy++ )
          Ilk[_istg-1][_ifct*_nx+iy] = _Iy[iy];

        // Add function contribution to adjoint bounds (discontinuities)
        if( _istg > 1  ){
          _pos_fct = ( _vFCT.size()>=ns? _istg-1:0 );
          if( _pos_fct
           && !_CC_I_ADJ( options, _t, NV_DATA_S(_Ny[_ifct]), _pos_fct, _ifct )
           && !_CC_I_QUAD( options, _t, NV_DATA_S(_Nyq[_ifct]) ) )
            { _END_ADJ(); return FATAL; }

          // Reset ODE solver - needed in case of discontinuity
          if( !_CC_CVODES( _ifct, _indexB[_ifct] ) )
            { _END_ADJ(); return FATAL; }
        }
        // Add initial state contribution to derivative bounds
        else if( !_IC_I_ADJ( _t ) ){ _END_ADJ(); return FATAL; }
        for( unsigned iq=0; iq<_np; iq++ )
          Idf[_ifct*_np+iq] = _Iyq[iq];
      }

      // Display & record adjoint intermediate results
      if( options.DISPLAY >= 1 ){
        _print_interm( tk[_istg-1], _nf*_nx, Ilk[_istg-1], "l", os );
        //_print_interm( _nf*_np, Idf, "df", os );
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
ODEBNDS_SUNDIALS<T,PMT,PVT>::_INI_PM_ADJ
( const unsigned np, const PVT *PMp )
{

  // Initialize bound propagation
  if( !ODEBNDS_BASE<T,PMT,PVT>::_INI_PM_ADJ( options, np, PMp ) )
    return false;

  // Set SUNDIALS adjoint/quadrature arrays
  unsigned cv_Ny_size;
  switch( options.WRAPMIT){
  case Options::NONE:
    cv_Ny_size = (_PMenv->nmon()+1)*_nx;
    break;
  case Options::DINEQ:
    cv_Ny_size = _PMenv->nmon()*_nx + 2*_nx ;
    break;
  case Options::ELLIPS:
  default:
    cv_Ny_size = _PMenv->nmon()*_nx + _nx*(_nx+1)/2;
    break;
  }
  unsigned cv_Nyq_size = (_PMenv->nmon()+1)*_np;

  for( unsigned i=0; i<_nf; i++ ){
    if( !_Ny[i] || cv_Ny_size != NV_LENGTH_S( _Ny[i] ) ){
      if( _Ny[i] ) N_VDestroy_Serial( _Ny[i] );
      _Ny[i] = N_VNew_Serial( cv_Ny_size );
    }
    if( !_Nyq[i] || cv_Nyq_size != NV_LENGTH_S( _Nyq[i] ) ){
      if( _Nyq[i] ) N_VDestroy_Serial( _Nyq[i] );
      _Nyq[i] = N_VNew_Serial( cv_Nyq_size );
    }
  }

  // Reset result record and statistics
  _results_adj.clear();
  _init_stats( stats_adj );

  return true;
}

template <typename T, typename PMT, typename PVT> inline int
ODEBNDS_SUNDIALS<T,PMT,PVT>::MC_CVADJRHSPM__
( realtype t, N_Vector x, N_Vector y, N_Vector ydot, void *user_data )
{
  ODEBNDS_SUNDIALS<T,PMT,PVT> *pODEBNDS = ODEBNDS_SUNDIALS<T,PMT,PVT>::pODEBNDS;
  pODEBNDS->_ifct = *static_cast<unsigned*>( user_data );
  bool flag = pODEBNDS->_RHS_PM_ADJ( pODEBNDS->options, t, NV_DATA_S( y ),
    NV_DATA_S( ydot ), NV_DATA_S( x ), pODEBNDS->_ifct, true );
  ODEBNDS_SUNDIALS<T,PMT,PVT>::pODEBNDS = pODEBNDS;
  pODEBNDS->stats_adj.numRHS++;
  return( (flag && _diam( pODEBNDS->_nx, pODEBNDS->_PMy ) < pODEBNDS->options.DMAX)? 0: -1 );
}

template <typename T, typename PMT, typename PVT> inline int
ODEBNDS_SUNDIALS<T,PMT,PVT>::MC_CVADJQUADPM__
( realtype t, N_Vector x, N_Vector y, N_Vector qdot, void *user_data )
{
  ODEBNDS_SUNDIALS<T,PMT,PVT> *pODEBNDS = ODEBNDS_SUNDIALS<T,PMT,PVT>::pODEBNDS;
  pODEBNDS->_ifct = *static_cast<unsigned*>( user_data );
  bool flag = pODEBNDS->_RHS_PM_QUAD( pODEBNDS->options, t, NV_DATA_S( y ),
    NV_DATA_S( qdot ), NV_DATA_S( x ), pODEBNDS->_ifct );
  ODEBNDS_SUNDIALS<T,PMT,PVT>::pODEBNDS = pODEBNDS;
  return( flag? 0: -1 );
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
inline typename ODEBNDS_SUNDIALS<T,PMT,PVT>::STATUS
ODEBNDS_SUNDIALS<T,PMT,PVT>::bounds_ASA
( const unsigned ns, const double*tk, const PVT*PMp, PVT**PMxk, PVT*PMq, PVT*PMf, PVT**PMlk,
  PVT*PMdf, std::ostream&os )
{
  // Compute state bounds and store intermediate results
  STATUS flag = NORMAL;
  ODEBND_SUNDIALS<T,PMT,PVT>::options = options;
  flag = ODEBND_SUNDIALS<T,PMT,PVT>::_bounds( ns, tk, PMp, PMxk, PMq, PMf, true, os);
  if( flag != NORMAL ) return flag;

  // Nothing to do if no functions are defined
  if( !_nf ) return NORMAL;

  // Check size
  if( !PMlk || !PMdf ) return FATAL;

  try{
    // Initialize adjoint bound integration
    if( !_INI_PM_ADJ( _np, PMp )) return FATAL;
    _t = tk[ns];

    // Bounds on terminal adjoints/quadratures
    if( PMlk && !PMlk[ns] ) PMlk[ns] = new PVT[_nx*_nf];
    for( _ifct=0; _ifct < _nf; _ifct++ ){
      _pos_fct = ( _vFCT.size()>=ns? ns-1:0 );
      if( !_TC_PM_ADJ( options, _t, _vec_sta[ns].data(), NV_DATA_S(_Ny[_ifct]), _pos_fct, _ifct )
       || ( _Nyq && _Nyq[_ifct] && !_TC_PM_QUAD( options, _t, NV_DATA_S(_Nyq[_ifct]) ) ) )
        { _END_ADJ(); return FATAL; }
      for( unsigned iy=0; PMlk[ns] && iy<_nx; iy++ )
        PMlk[ns][_ifct*_nx+iy] = _PMy[iy];
      for( unsigned iq=0; iq<_np; iq++ )
        PMdf[_ifct*_np+iq] = _PMyq[iq];
    }
    // Display & record adjoint terminal results
    if( options.DISPLAY >= 1 ){
      _print_interm( tk[ns], _nx*_nf, PMlk[ns], "l", os );
      //_print_interm( _np, PMdf, "df", os );
    }
    if( options.RESRECORD )
      _results_adj.push_back( Results( tk[ns], _nf*_nx, PMlk[ns] ) );

    // Initialization of adjoint integration
    for( _ifct=0; _ifct < _nf; _ifct++ )
      if( !_INI_CVODES( MC_CVADJRHSPM__, MC_CVADJQUADPM__, _ifct,
        _indexB[_ifct], _iusrB[_ifct] ) )
        { _END_ADJ(); return FATAL;}

    // Integrate adjoint ODEs through each stage using SUNDIALS
    pODEBNDS = this;
    for( _istg=ns; _istg>0; _istg-- ){

      // Update list of operations in RHSADJ and QUADADJ
      _pos_rhs  = ( _vRHS.size() <=1? 0: _istg-1 );
      _pos_quad = ( _vQUAD.size()<=1? 0: _istg-1 );
      if( !_SET_PM_ADJ( options, _pos_rhs, _pos_quad, _pos_fct ) )
        { _END_ADJ(); return FATAL; }

      // Propagate bounds backward to previous stage time
      _cv_flag = CVodeB( _cv_mem, tk[_istg-1], CV_NORMAL );
      if( _check_cv_flag(&_cv_flag, "CVodeB", 1) )
        { _END_ADJ(); return FATAL; }
      _t = tk[_istg-1];

      // Bounds on states/adjoints/quadratures at stage time
      //stats_adj.numSteps = 0; 
      for( _ifct=0; _ifct < _nf; _ifct++ ){
        _cv_flag = CVodeGetB( _cv_mem, _indexB[_ifct], &_t, _Ny[_ifct]);
        if( _check_cv_flag( &_cv_flag, "CVodeGetB", 1) )
          { _END_ADJ(); return FATAL; }
        _cv_flag = CVodeGetQuadB( _cv_mem, _indexB[_ifct], &_t, _Nyq[_ifct]);
        if( _check_cv_flag( &_cv_flag, "CVodeGetQuadB", 1) )
          { _END_ADJ(); return FATAL; }
        void *cv_memB = CVodeGetAdjCVodeBmem(_cv_mem, _indexB[_ifct] );
        long int nstpB;
        _cv_flag = CVodeGetNumSteps( cv_memB, &nstpB );
        stats_adj.numSteps += nstpB;
        switch( options.WRAPMIT){
        case Options::NONE:
          _vec2PMI( _vec_sta[_istg-1].data(), _PMenv, _nx, _PMx, true );
          _vec2PMI( NV_DATA_S(_Ny[_ifct]), _PMenv, _nx, _PMy, true );
          break;
        case Options::DINEQ:
          _vec2PMI( _vec_sta[_istg-1].data(), _PMenv, _nx, _PMx );
          _vec2PMI( NV_DATA_S(_Ny[_ifct]), _PMenv, _nx, _PMy );
          break;
        case Options::ELLIPS:
        default:
          _vec2PME( _vec_sta[_istg-1].data(), _PMenv, _nx, _PMx, _Q, _Er, _Ir );
          _vec2PME( NV_DATA_S(_Ny[_ifct]), _PMenv, _nx, _PMy, _Qy, _Edy, _Idy );
          break;
        }
        _vec2PMI( NV_DATA_S( _Nyq[_ifct] ), _PMenv, _np, _PMyq, true );
        if( PMlk && !PMlk[_istg-1] ) PMlk[_istg-1] = new PVT[_nx*_nf];
        for( unsigned iy=0; PMlk[_istg-1] && iy<_nx; iy++)
          PMlk[_istg-1][_ifct*_nx+iy] =_PMy[iy];

        // Add function contribution to adjoint bounds (discontinuities)
        if( _istg > 1  ){
          _pos_fct = ( _vFCT.size()>=ns? _istg-1:0 );
          if( _pos_fct 
           && !_CC_PM_ADJ( options, _t, NV_DATA_S(_Ny[_ifct]), _pos_fct, _ifct )
           && !_CC_PM_QUAD( options, _t, NV_DATA_S(_Nyq[_ifct]) ) )
            { _END_ADJ(); return FATAL; }

          // Reset ODE solver - needed in case of discontinuity
          if( !_CC_CVODES( _ifct, _indexB[_ifct] ) )
            { _END_ADJ(); return FATAL; }
        }
        // Add initial state contribution to derivative bounds
        else if( !_IC_PM_ADJ( _t ) ){ _END_ADJ(); return FATAL; }
        for( unsigned iq=0; iq<_np; iq++ )
          PMdf[_ifct*_np+iq] = _PMyq[iq];
      }

      // Display & record adjoint intermediate results
      if( options.DISPLAY >= 1 ){
        _print_interm( tk[_istg-1], _nf*_nx, PMlk[_istg-1], "l", os );
        //_print_interm( _np, PMdf, "df", os );
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




















