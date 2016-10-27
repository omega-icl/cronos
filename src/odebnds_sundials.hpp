// Copyright (C) 2015- Benoit Chachuat & Nikola Peric, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

#ifndef MC__ODEBNDS_SUNDIALS_HPP
#define MC__ODEBNDS_SUNDIALS_HPP

#undef  MC__ODEBNDS_SUNDIALS_DINEQI_DEBUG
#undef  MC__ODEBNDS_SUNDIALS_DINEQPM_DEBUG

#include "odebnds_base.hpp"
#include "odebnd_sundials.hpp"
#include "odeslvs_sundials.hpp"

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
  typedef int (*CVSENRhs1Fn)( int Ns, realtype t, N_Vector x, N_Vector xdot, int is,
    N_Vector y, N_Vector ydot, void *user_data, N_Vector tmp1, N_Vector tmp2 );
  typedef int (*CVSENQuadFn)( int Ns, realtype t, N_Vector x, N_Vector *y, N_Vector qdot,
    N_Vector *qSdot, void *user_data, N_Vector tmp1, N_Vector tmp2 );

  using ODEBND_BASE<T,PMT,PVT>::NORMAL; 
  using ODEBND_BASE<T,PMT,PVT>::FAILURE;
  using ODEBND_BASE<T,PMT,PVT>::FATAL;

  using ODEBND_BASE<T,PMT,PVT>::_diam;
  using ODEBND_BASE<T,PMT,PVT>::_print_interm;
  using ODEBND_BASE<T,PMT,PVT>::_Ix;
  using ODEBND_BASE<T,PMT,PVT>::_PMenv;
  using ODEBND_BASE<T,PMT,PVT>::_PMx;

  using ODEBNDS_BASE<T,PMT,PVT>::_Ir;
  using ODEBNDS_BASE<T,PMT,PVT>::_Er;
  using ODEBNDS_BASE<T,PMT,PVT>::_Iy;
  using ODEBNDS_BASE<T,PMT,PVT>::_Iyq;
  using ODEBNDS_BASE<T,PMT,PVT>::_Edy;
  using ODEBNDS_BASE<T,PMT,PVT>::_PMy;
  using ODEBNDS_BASE<T,PMT,PVT>::_PMyq;
  using ODEBNDS_BASE<T,PMT,PVT>::_IC_SET_FSA;
  using ODEBNDS_BASE<T,PMT,PVT>::_CC_SET_FSA;
  using ODEBNDS_BASE<T,PMT,PVT>::_CC_SET_ASA;
  using ODEBNDS_BASE<T,PMT,PVT>::_RHS_SET_FSA;
  using ODEBNDS_BASE<T,PMT,PVT>::_RHS_SET_ASA;
  using ODEBNDS_BASE<T,PMT,PVT>::_GET_I_SEN;
  using ODEBNDS_BASE<T,PMT,PVT>::_INI_I_SEN;
  using ODEBNDS_BASE<T,PMT,PVT>::_IC_I_SET_ASA;
  using ODEBNDS_BASE<T,PMT,PVT>::_IC_I_SEN;
  using ODEBNDS_BASE<T,PMT,PVT>::_IC_I_QUAD_ASA;
  using ODEBNDS_BASE<T,PMT,PVT>::_TC_I_SET_ASA;
  using ODEBNDS_BASE<T,PMT,PVT>::_TC_I_SEN;
  using ODEBNDS_BASE<T,PMT,PVT>::_TC_I_QUAD_ASA;
  using ODEBNDS_BASE<T,PMT,PVT>::_CC_I_SET;
  using ODEBNDS_BASE<T,PMT,PVT>::_CC_I_SEN;
  using ODEBNDS_BASE<T,PMT,PVT>::_CC_I_QUAD_ASA;
  using ODEBNDS_BASE<T,PMT,PVT>::_RHS_I_SET;
  using ODEBNDS_BASE<T,PMT,PVT>::_RHS_I_SEN;
  using ODEBNDS_BASE<T,PMT,PVT>::_RHS_I_QUAD;
  using ODEBNDS_BASE<T,PMT,PVT>::_FCT_I_SEN;
  using ODEBNDS_BASE<T,PMT,PVT>::_GET_PM_SEN;
  using ODEBNDS_BASE<T,PMT,PVT>::_INI_PM_SEN;
  using ODEBNDS_BASE<T,PMT,PVT>::_IC_PM_SET_ASA;
  using ODEBNDS_BASE<T,PMT,PVT>::_IC_PM_SEN;
  using ODEBNDS_BASE<T,PMT,PVT>::_IC_PM_QUAD_ASA;
  using ODEBNDS_BASE<T,PMT,PVT>::_TC_PM_SET_ASA;
  using ODEBNDS_BASE<T,PMT,PVT>::_TC_PM_SEN;
  using ODEBNDS_BASE<T,PMT,PVT>::_TC_PM_QUAD_ASA;
  using ODEBNDS_BASE<T,PMT,PVT>::_CC_PM_SET;
  using ODEBNDS_BASE<T,PMT,PVT>::_CC_PM_SEN;
  using ODEBNDS_BASE<T,PMT,PVT>::_CC_PM_QUAD_ASA;
  using ODEBNDS_BASE<T,PMT,PVT>::_RHS_PM_SET;
  using ODEBNDS_BASE<T,PMT,PVT>::_RHS_PM_SEN;
  using ODEBNDS_BASE<T,PMT,PVT>::_RHS_PM_QUAD;
  using ODEBNDS_BASE<T,PMT,PVT>::_FCT_PM_SEN;
  using ODEBNDS_BASE<T,PMT,PVT>::_bounds_ASA;
  using ODEBNDS_BASE<T,PMT,PVT>::_bounds_FSA;

  using ODEBND_SUNDIALS<T,PMT,PVT>::_cv_mem;
  using ODEBND_SUNDIALS<T,PMT,PVT>::_cv_flag;
  using ODEBND_SUNDIALS<T,PMT,PVT>::_Nx;
  using ODEBND_SUNDIALS<T,PMT,PVT>::_Nq;
  using ODEBND_SUNDIALS<T,PMT,PVT>::_Nf;
  using ODEBND_SUNDIALS<T,PMT,PVT>::_vec_sta;
  using ODEBND_SUNDIALS<T,PMT,PVT>::_pos_ic;
  using ODEBND_SUNDIALS<T,PMT,PVT>::_pos_rhs;
  using ODEBND_SUNDIALS<T,PMT,PVT>::_pos_quad;
  using ODEBND_SUNDIALS<T,PMT,PVT>::_pos_fct;
  using ODEBND_SUNDIALS<T,PMT,PVT>::_init_stats;
  using ODEBND_SUNDIALS<T,PMT,PVT>::_final_stats;
  using ODEBND_SUNDIALS<T,PMT,PVT>::_print_stats;
  using ODEBND_SUNDIALS<T,PMT,PVT>::_END_STA;
  using ODEBND_SUNDIALS<T,PMT,PVT>::stats_sta;
  using ODEBND_SUNDIALS<T,PMT,PVT>::results_sta;

public:
  using ODEBND_SUNDIALS<T,PMT,PVT>::options;

 protected:
  //! @brief current function index
  unsigned _ifct;

  //! @brief current parameter sensitivity index
  unsigned _isen;

  //! @brief size of N_vector arrays
  unsigned _nvec;

  //! @brief N_Vector object holding current sensitivity/adjoint parameterizations
  N_Vector *_Ny;

  //! @brief N_Vector object holding current quadrature sensitivity parameterizations
  N_Vector *_Nyq;
  
  //! @brief N_Vector object holding current state function derivative parameterizations
  N_Vector *_Nfp;

  //! @brief N_Vector object holding absolute tolerances
  N_Vector _NTOLy;

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

  //! @brief static pointer to local integrator class
  ODESLVS_SUNDIALS pODESLVS;

 public:
  /** @defgroup ODEBNDS Continuous-time set-valued integration of parametric ODEs with sensitivity analysis
   *  @{
   */
  typedef BASE_SUNDIALS::Stats Stats;
  typedef typename ODEBND_BASE<T,PMT,PVT>::Results Results;
  typedef typename ODEBND_SUNDIALS<T,PMT,PVT>::Exceptions Exceptions;
  typedef typename ODEBND_SUNDIALS<T,PMT,PVT>::Options Options;

  //! @brief Default constructor
  ODEBNDS_SUNDIALS
    ();

  //! @brief Virtual destructor
  virtual ~ODEBNDS_SUNDIALS
    ();

  //! @brief Statistics for adjoint integration
  Stats stats_sen;

  //! @brief Vector storing interval adjoint bounds (see Options::RESRECORD)
  std::vector< Results > results_sen;

  //! @brief Propagate state/quadrature interval bounds forward in time through every time stages
  STATUS bounds
    ( const unsigned ns, const double*tk, const T*Ip, T**Ixk,
      T*If, std::ostream&os=std::cout )
      { return ODEBND_SUNDIALS<T,PMT,PVT>::bounds( ns, tk, Ip, Ixk, If, os); }

  //! @brief Propagate state/quadrature polynomial models forward in time through every time stages
  STATUS bounds
    ( const unsigned ns, const double*tk, const PVT*PMp, PVT**PMxk,
      PVT*PMf, std::ostream&os=std::cout )
      { return ODEBND_SUNDIALS<T,PMT,PVT>::bounds( ns, tk, PMp, PMxk, PMf, os); }

  //! @brief Propagate state/quadrature polynomial models forward in time through every time stages
  STATUS bounds
    ( const unsigned ns, const double*tk, const PVT*PMp, PVT**PMxk,
      E*ERxk, PVT*PMf, std::ostream&os=std::cout )
      { return ODEBND_SUNDIALS<T,PMT,PVT>::bounds( ns, tk, PMp, PMxk, ERxk, PMf, os); }

 //! @brief Propagate state and sensitivity interval bounds forward in time through every time stages
  STATUS bounds_FSA
    ( const unsigned ns, const double*tk, const T*Ip, T**Ixk, T*If, T**Ixpk,
      T*Ifp, std::ostream&os=std::cout );

 //! @brief Propagate state and sensitivity polynomial models forward in time through every time stages
  STATUS bounds_FSA
    ( const unsigned ns, const double*tk, const PVT*PMp, PVT**PMxk, PVT*PMf,
      PVT**PMxpk, PVT*PMfp, std::ostream&os=std::cout );

 //! @brief Propagate state and sensitivity polynomial models forward in time through every time stages
  STATUS bounds_FSA
    ( const unsigned ns, const double*tk, const PVT*PMp, PVT**PMxk, E*ERxk,
      PVT*PMf, PVT**PMxpk, E**ERxpk, PVT*PMfp, std::ostream&os=std::cout );

 //! @brief Compute approximate state and sensitivity interval bounds forward in time using sampling through every time stages
  STATUS bounds_FSA
    ( const unsigned ns, const double*tk, const T*Ip, T**Ixk, T*If, T**Ixpk,
      T*Ifp, const unsigned nsamp, std::ostream&os=std::cout );

  //! @brief Propagate state and adjoint interval bounds forward and backward in time through every time stages
  STATUS bounds_ASA
    ( const unsigned ns, const double*tk, const T*Ip, T**Ixk, T*If, T**Ilk, T*Ifp,
      std::ostream&os=std::cout );

  //! @brief Propagate state and adjoint polynomial models forward and backward in time through every time stages
  STATUS bounds_ASA
    ( const unsigned ns, const double*tk, const PVT*PMp, PVT**PMxk, PVT*PMf,
      PVT**PMlk, PVT*PMfp, std::ostream&os=std::cout );

  //! @brief Propagate state and adjoint polynomial models forward and backward in time through every time stages
  STATUS bounds_ASA
    ( const unsigned ns, const double*tk, const PVT*PMp, PVT**PMxk, E*ERxk,
      PVT*PMf, PVT**PMlk, E**ERlk, PVT*PMfp, std::ostream&os=std::cout );

 //! @brief Compute approximate state and adjoint interval bounds forward and backward in time using sampling through every time stages
  STATUS bounds_ASA
    ( const unsigned ns, const double*tk, const T*Ip, T**Ixk, T*If, T**Ilk,
      T*Ifp, const unsigned nsamp, std::ostream&os=std::cout );


  //! @brief Record state and sensitivity bounds in files <a>obndsta</a> and <a>obndsa</a>, with accuracy of <a>iprec</a> digits
  void record
    ( std::ofstream&obndsta, std::ofstream&obndsa, const unsigned iprec=5 ) const
    { this->ODEBND_SUNDIALS<T,PMT,PVT>::record( obndsta, iprec );
      this->ODEBND_BASE<T,PMT,PVT>::_record( obndsa, results_sen, iprec ); }

  //! @brief Record state and sensitivity bounds in files <a>obndsta</a> and <a>obndsa</a>, with accuracy of <a>iprec</a> digits
  void record
    ( std::ofstream&obndsta, const unsigned iprec=5 ) const
    { this->ODEBND_SUNDIALS<T,PMT,PVT>::record( obndsta, iprec ); }
  /** @} */

 private:
  //! @brief Function to initialize CVode memory block (virtual)
  virtual bool _INI_CVODE
    ( CVRhsFn MC_CVRHS, CVRhsFn MC_CVQUAD );

  //! @brief Function to initialize CVodeS memory block for forward sensitivity
  bool _INI_CVODES
    ( CVSENRhs1Fn MC_CVSENRHS, CVSENQuadFn MC_CVSENQUAD );

  //! @brief Function to initialize CVodeS memory block for adjoint sensitivity
  bool _INI_CVODES
    ( CVADJRhsFn MC_CVADJRHS, CVADJRhsFn MC_CVADJQUAD, const unsigned ifct,
      int&indexB, unsigned&iusrB );

  //! @brief Function to reinitialize CVodeS memory block for forward sensitivity
  bool _CC_CVODES_FSA
    ();

  //! @brief Function to reinitialize CVodeS memory block for forward sensitivity
  bool _CC_CVODES_QUAD
    ();

  //! @brief Function to reinitialize CVodeS memory block for adjoint sensitivity
  bool _CC_CVODES_ASA
    ( const unsigned ifct, const int indexB );

  //! @brief Function to finalize sensitivity/adjoint bounding
  void _END_SEN
    ();

  //! @brief Function to initialize adjoint interval bounding
  bool _INI_I_ASA
    ( const unsigned np, const T *Ip );

  //! @brief Static wrapper to function computing the adjoint DINEQ RHS
  static int MC_CVASARHSI__
    ( realtype t, N_Vector x, N_Vector y, N_Vector ydot, void *user_data );

  //! @brief Static wrapper to function computing the adjoint quadrature RHS
  static int MC_CVASAQUADI__
    ( realtype t, N_Vector x, N_Vector y, N_Vector qdot, void *user_data );

  //! @brief Function to initialize adjoint polynomial models
  bool _INI_PM_ASA
    ( const unsigned np, const PVT *PMp );

  //! @brief Static wrapper to function to calculate the adjoint DINEQ-PM RHS
  static int MC_CVASARHSPM__
    ( realtype t, N_Vector x, N_Vector y, N_Vector ydot, void *user_data );

  //! @brief Static wrapper to function computing the adjoint PM-quadrature RHS
  static int MC_CVASAQUADPM__
    ( realtype t, N_Vector x, N_Vector y, N_Vector qdot, void *user_data );

  //! @brief Function to initialize adjoint polynomial models
  bool _INI_I_FSA
    ( const unsigned np, const T *Ip );

  //! @brief Static wrapper to function to calculate the sensitivity DINEQ RHS
  static int MC_CVFSARHSI__
    ( int Ns, realtype t, N_Vector x, N_Vector xdot, int is, N_Vector y,
      N_Vector ydot, void *user_data, N_Vector tmp1, N_Vector tmp2 );

  //! @brief Static wrapper to function computing the sensitivity quadrature RHS
  static int MC_CVFSAQUADI__
    ( int Ns, realtype t, N_Vector x, N_Vector *y, N_Vector qdot, N_Vector *qSdot, 
      void *user_data, N_Vector tmp1, N_Vector tmp2 );

  //! @brief Function to initialize adjoint polynomial models
  bool _INI_PM_FSA
    ( const unsigned np, const PVT *PMp );

  //! @brief Static wrapper to function to calculate the sensitivity DINEQ-PM RHS
  static int MC_CVFSARHSPM__
    ( int Ns, realtype t, N_Vector x, N_Vector xdot, int is, N_Vector y,
      N_Vector ydot, void *user_data, N_Vector tmp1, N_Vector tmp2 );

  //! @brief Static wrapper to function computing the sensitivity PM-quadrature RHS
  static int MC_CVFSAQUADPM__
    ( int Ns, realtype t, N_Vector x, N_Vector *y, N_Vector qdot, N_Vector *qSdot, 
      void *user_data, N_Vector tmp1, N_Vector tmp2 );

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
  _ifct(0), _isen(0), _nvec(0), _Ny(0), _Nyq(0), _Nfp(0), _NTOLy(0),
  _indexB(0), _iusrB(0), pODESLVS()
{}

template <typename T, typename PMT, typename PVT> inline
ODEBNDS_SUNDIALS<T,PMT,PVT>::~ODEBNDS_SUNDIALS
()
{
  if( _Ny )    N_VDestroyVectorArray_Serial( _Ny,  _nvec );
  if( _Nyq )   N_VDestroyVectorArray_Serial( _Nyq, _nvec );
  if( _Nfp )   N_VDestroyVectorArray_Serial( _Nfp, _nvec );
  if( _NTOLy ) N_VDestroy_Serial( _NTOLy );
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
  //{std::cout<<"here!\n"; int dum; std::cin>>dum;}
  delete[] _indexB; _indexB = new int[_nf];
  delete[] _iusrB;  _iusrB  = new unsigned[_nf];

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
       //_cv_flag = CVDlsSetDenseJacFnB( _cv_mem, indexB, NULL );
       //if ( _check_cv_flag(&_cv_flag, "CVDlsSetDenseJacFnB", 1) ) return false;
       //break;

     case Options::CV_DENSE:
       _cv_flag = CVDenseB( _cv_mem, indexB, NV_LENGTH_S( _Ny[ifct] ) );
       if( _check_cv_flag(&_cv_flag, "CVDenseB", 1)) return false;
       _cv_flag = CVDlsSetDenseJacFnB( _cv_mem, indexB, NULL );
       if ( _check_cv_flag(&_cv_flag, "CVDlsSetDenseJacFnB", 1) ) return false;
       break;
    }
  }

  // Specify the relative and absolute tolerances for adjoints, with
  // different tolerances for entries of ellipsoid shape matrix
  if( (options.WRAPMIT == Options::ELLIPS) && (options.ETOLB < options.ATOLB) ){
    for(unsigned i=0; i<NV_LENGTH_S(_NTOLy); i++ )
      NV_Ith_S(_NTOLy,i) = i<NV_LENGTH_S(_NTOLy)-_nx*(_nx+1)/2? options.ATOLB: options.ETOLB;
    _cv_flag = CVodeSVtolerancesB( _cv_mem, indexB, options.RTOL, _NTOLy );
    if( _check_cv_flag(&_cv_flag, "CVodeSVtolerancesB", 1) ) return false;  
  }
  else{
    _cv_flag = CVodeSStolerancesB( _cv_mem, indexB, options.RTOLB, options.ATOLB );
    if( _check_cv_flag(&_cv_flag, "CVodeSStolerancesB", 1) ) return false;
  }

  // Set maximum number of error test failures
  //_cv_flag = CVodeSetMaxErrTestFailsB( _cv_mem, indexB, options.MAXFAIL );
  //if ( _check_cv_flag(&_cv_flag, "CVodeSetMaxErrTestFailsB", 1) ) return false;

  // Set maximum number of error test failures
  //_cv_flag = CVodeSetMaxConvFailsB( _cv_mem, indexB, options.MAXFAIL );
  //if ( _check_cv_flag(&_cv_flag, "CVodeSetMaxConvFailsB", 1) ) return false;

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
ODEBNDS_SUNDIALS<T,PMT,PVT>::_INI_CVODES
( CVSENRhs1Fn MC_CVSENRHS, CVSENQuadFn MC_CVSENQUAD )
{
  // Allocate memory for sensitivity integration
  _cv_flag = CVodeSensInit1( _cv_mem, _np, options.FSACORR, MC_CVSENRHS, _Ny );
  if( _check_cv_flag(&_cv_flag, "CVodeSensInit1", 1) ) return false;

  // Specify absolute and relative tolerances for sensitivities
  if( options.AUTOTOLS ){
    _cv_flag = CVodeSensEEtolerances( _cv_mem );
    if( _check_cv_flag(&_cv_flag, "CVodeSensEEtolerances", 1)) return false;
  }
  else if( (options.WRAPMIT == Options::ELLIPS) && (options.ETOLS < options.ATOLS) ){
    N_Vector* ATOLS = N_VCloneVectorArray_Serial( _nvec, _NTOLy );
    for(unsigned jp=0; jp<_np; jp++)
      for(unsigned iy=0; iy<NV_LENGTH_S(ATOLS[jp]); iy++ )
        NV_Ith_S(ATOLS[jp],iy) = iy<NV_LENGTH_S(ATOLS[jp])-_nx*(_nx+1)/2? options.ATOLS: options.ETOLS;
    _cv_flag = CVodeSensSVtolerances( _cv_mem, options.RTOLS, ATOLS );
    N_VDestroyVectorArray_Serial( ATOLS, _nvec );
    if( _check_cv_flag(&_cv_flag, "CVodeSensSVtolerances", 1)) return false;     
  }
  else{
    realtype ATOLS[_np];
    for( unsigned int ip=0; ip<_np; ip++ ) ATOLS[ip] = options.ATOLS;
    _cv_flag = CVodeSensSStolerances( _cv_mem, options.RTOLS, ATOLS );
    if( _check_cv_flag(&_cv_flag, "CVodeSensSStolerances", 1)) return false;
  }

  // Specify the error control strategy for sensitivity variables
  _cv_flag = CVodeSetSensErrCon( _cv_mem, options.FSAERR );
  if( _check_cv_flag(&_cv_flag, "CVodeSetSensErrCon", 1) ) return false;

  // Specify problem parameter information for sensitivity calculations
  _cv_flag = CVodeSetSensParams( _cv_mem, 0, 0, 0 );
  if( _check_cv_flag(&_cv_flag, "CVodeSetSensParams", 1) ) return false;

  // Initialize integrator memory for quadratures
  if( !_nq ) return true;
  _cv_flag = CVodeQuadSensInit( _cv_mem, MC_CVSENQUAD, _Nyq );
  if( _check_cv_flag(&_cv_flag, "CVodeQuadSensInit", 1) ) return false;
  
  // Specify whether or not to perform error control on quadrature
  _cv_flag = CVodeSetQuadSensErrCon( _cv_mem, options.QERRS );
  if( _check_cv_flag(&_cv_flag, "CVodeSetQuadSensErrCon", 1) ) return false;
  
  // Specify absolute and relative tolerances for quadratures
  if( options.AUTOTOLS ){
    _cv_flag = CVodeQuadSensEEtolerances( _cv_mem );
    if( _check_cv_flag(&_cv_flag, "CVodeQuadSensEEtolerances", 1) ) return false;
  }
  else{
    realtype ATOLS[_np];
    for( unsigned int ip=0; ip<_np; ip++ ) ATOLS[ip] = options.ATOLS;
    _cv_flag = CVodeQuadSensSStolerances(_cv_mem, options.RTOLS, ATOLS);
    if( _check_cv_flag(&_cv_flag, "CVodeQuadSensSStolerances", 1) ) return false;
  }

  return true;
}

template <typename T, typename PMT, typename PVT> inline bool
ODEBNDS_SUNDIALS<T,PMT,PVT>::_CC_CVODES_ASA
( const unsigned ifct, const int indexB )
{
  // Reinitialize CVodeS memory block for current time _t and adjoint _Ny
  _cv_flag = CVodeReInitB( _cv_mem, indexB, _t, _Ny[ifct] );
  if( _check_cv_flag(&_cv_flag, "CVodeReInitB", 1) ) return false;

  // Reinitialize CVodeS memory block for current adjoint quarature _Nyq
  if( !_np ) return true;
  _cv_flag = CVodeQuadReInitB( _cv_mem, indexB, _Nyq[ifct] );
  if( _check_cv_flag(&_cv_flag, "CVodeQuadReInitB", 1) ) return false;

  return true;
}

template <typename T, typename PMT, typename PVT> inline bool
ODEBNDS_SUNDIALS<T,PMT,PVT>::_CC_CVODES_FSA
()
{
  // Reinitialize CVodeS memory block for current sensitivity _Ny
  _cv_flag = CVodeSensReInit( _cv_mem, options.FSACORR, _Ny );
  if( _check_cv_flag(&_cv_flag, "CVodeSensReInit", 1) ) return false;

  return true;
}

template <typename T, typename PMT, typename PVT> inline bool
ODEBNDS_SUNDIALS<T,PMT,PVT>::_CC_CVODES_QUAD
()
{
  // Reinitialize CVode memory block for current sensitivity quarature _Nyq
  if( !_nq ) return true;
  _cv_flag = CVodeQuadSensReInit( _cv_mem, _Nyq );
  if( _check_cv_flag(&_cv_flag, "CVodeQuadSensReInit", 1) ) return false;

  return true;
}

template <typename T, typename PMT, typename PVT> inline bool
ODEBNDS_SUNDIALS<T,PMT,PVT>::_INI_I_ASA
( const unsigned np, const T*Ip )
{
  // Initialize bound propagation
  if( !_INI_I_SEN( options, np, Ip, _nf, np, options.ETOLB ) )
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
  if( _nvec != _nf ){
    if( _Ny )   N_VDestroyVectorArray_Serial( _Ny,  _nvec );
    if( _Nyq )  N_VDestroyVectorArray_Serial( _Nyq, _nvec );
     _nvec = _nf;
    _Ny  = N_VCloneVectorArray_Serial( _nvec, _Nx );
    _Nyq = N_VCloneVectorArray_Serial( _nvec, _Nx );
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
  if( !_NTOLy || cv_Ny_size != NV_LENGTH_S( _NTOLy ) )
    if( _NTOLy ) N_VDestroy_Serial( _NTOLy );
      _NTOLy = N_VNew_Serial( cv_Ny_size );

  // Reset result record and statistics
  results_sen.clear();
  _init_stats( stats_sen );

  return true;
}

template <typename T, typename PMT, typename PVT> inline void
ODEBNDS_SUNDIALS<T,PMT,PVT>::_END_SEN()
{
  // Get final CPU time
  _final_stats( stats_sen );
}

template <typename T, typename PMT, typename PVT> inline int
ODEBNDS_SUNDIALS<T,PMT,PVT>::MC_CVASARHSI__
( realtype t, N_Vector x, N_Vector y, N_Vector ydot, void *user_data )
{
  ODEBNDS_SUNDIALS<T,PMT,PVT> *pODEBNDS = ODEBNDS_SUNDIALS<T,PMT,PVT>::pODEBNDS;
  pODEBNDS->_ifct = *static_cast<unsigned*>( user_data );
  bool flag = pODEBNDS->_RHS_I_SEN( pODEBNDS->options, t, NV_DATA_S( x ), NV_DATA_S( y ),
    NV_DATA_S( ydot ), pODEBNDS->_ifct, true );
  ODEBNDS_SUNDIALS<T,PMT,PVT>::pODEBNDS = pODEBNDS;
  pODEBNDS->stats_sen.numRHS++;
  return( (flag && _diam( pODEBNDS->_nx, pODEBNDS->_Iy ) < pODEBNDS->options.DMAX)? 0: -1 );
}

template <typename T, typename PMT, typename PVT> inline int
ODEBNDS_SUNDIALS<T,PMT,PVT>::MC_CVASAQUADI__
( realtype t, N_Vector x, N_Vector y, N_Vector qdot, void *user_data )
{
  ODEBNDS_SUNDIALS<T,PMT,PVT> *pODEBNDS = ODEBNDS_SUNDIALS<T,PMT,PVT>::pODEBNDS;
  pODEBNDS->_ifct = *static_cast<unsigned*>( user_data );
  bool flag = pODEBNDS->_RHS_I_QUAD( pODEBNDS->options, pODEBNDS->_np, NV_DATA_S( qdot ),
    pODEBNDS->_ifct, true );
  ODEBNDS_SUNDIALS<T,PMT,PVT>::pODEBNDS = pODEBNDS;
  return( flag? 0: -1 );
}

//! @fn template <typename T, typename PMT, typename PVT> inline typename ODEBNDS_SUNDIALS<T,PMT,PVT>::STATUS ODEBNDS_SUNDIALS<T,PMT,PVT>::bounds_ASA(
//! const unsigned ns, const double*tk, const T*Ip, T**Ixk, 
//! T*If, T**Ilk, T*Ifp, std::ostream&os=std::cout )
//!
//! This function computes an enclosure of the reachable set of the parametric ODEs
//! using continuous-time propagation of convex sets (intervals, ellipsoids),
//! together with adjoint sensitivity bounds:
//!  - <a>ns</a>    [input]  number of time stages
//!  - <a>tk</a>    [input]  stage times, including the initial time
//!  - <a>Ip</a>    [input]  interval enclosure of parameter set
//!  - <a>Ixk</a>   [output] interval enclosure of state variables at stage times
//!  - <a>If</a>    [output] interval enclosure of state functions
//!  - <a>Ilk</a>   [output] interval enclosure of adjoint variables at stage times
//!  - <a>Ifp</a>   [output] interval enclosure of state function derivatives
//!  - <a>os</a>    [input]  output stream [default: std::cout]
//! .
//! The return value is the status.
template <typename T, typename PMT, typename PVT>
inline typename ODEBNDS_SUNDIALS<T,PMT,PVT>::STATUS
ODEBNDS_SUNDIALS<T,PMT,PVT>::bounds_ASA
( const unsigned ns, const double*tk, const T*Ip, T**Ixk, T*If,
  T**Ilk, T*Ifp, std::ostream&os )
{
  // Compute state bounds and store intermediate results
  STATUS flag = NORMAL;
  flag = ODEBND_SUNDIALS<T,PMT,PVT>::_bounds( ns, tk, Ip, Ixk, If, true, os);
  if( flag != NORMAL ) return flag;

  // Nothing to do if no functions are defined
  if( !_nf || !_np ) return NORMAL;

  // Check size
  if( !Ilk || (_nf && !Ifp) ) return FATAL;

  try{
    // Initialize adjoint bound integration
    if( !_INI_I_ASA( _np, Ip )) return FATAL;
    _t = tk[ns];

    // Bounds on terminal adjoints/quadratures
    if( Ilk && !Ilk[ns] ) Ilk[ns] = new T[_nx*_nf];
    for( _ifct=0; _ifct < _nf; _ifct++ ){
      _pos_fct = ( _vFCT.size()>=ns? ns-1:0 );
      if( !_TC_I_SET_ASA( options, _pos_fct, _ifct )
       || !_TC_I_SEN( options, _t, _vec_sta[ns].data(), NV_DATA_S(_Ny[_ifct]) )
       || ( _Nyq && _Nyq[_ifct] && !_TC_I_QUAD_ASA( options, NV_DATA_S(_Nyq[_ifct]) ) ) )
        { _END_SEN(); return FATAL; }
      for( unsigned iy=0; Ilk[ns] && iy<_nx; iy++ )
        Ilk[ns][_ifct*_nx+iy] = _Iy[iy];
      for( unsigned iq=0; iq<_np; iq++ )
        Ifp[iq*_nf+_ifct] = _Iyq[iq];
    }

    // Display & record adjoint terminal results
    if( options.DISPLAY >= 1 )
      _print_interm( tk[ns], _nf*_nx, Ilk[ns], "l", os );//_nf
    if( options.RESRECORD )
      results_sen.push_back( Results( tk[ns], _nx*_nf, Ilk[ns] ) );//_nf

    // Initialization of adjoint integration
    for( _ifct=0; _ifct < _nf; _ifct++ )
      if( !_INI_CVODES( MC_CVASARHSI__, MC_CVASAQUADI__, _ifct,
        _indexB[_ifct], _iusrB[_ifct] ) )
        { _END_SEN(); return FATAL;}

    // Integrate adjoint ODEs through each stage using SUNDIALS
    pODEBNDS = this;
    for( _istg=ns; _istg>0; _istg-- ){

      // Update list of operations in RHSADJ and QUADADJ
      _pos_rhs  = ( _vRHS.size() <=1? 0: _istg-1 );
      _pos_quad = ( _vQUAD.size()<=1? 0: _istg-1 );
      if( !_RHS_SET_ASA( _pos_rhs, _pos_quad, _pos_fct )
       || !_RHS_I_SET( options, _nf, _np ) )
        { _END_SEN(); return FATAL; }

      // Propagate bounds backward to previous stage time
      _cv_flag = CVodeB( _cv_mem, tk[_istg-1], CV_NORMAL );
      if( _check_cv_flag(&_cv_flag, "CVodeB", 1) )
        { _END_SEN(); return FATAL; }
      _t = tk[_istg-1];

      // Bounds on states/adjoints/quadratures at stage time
      //stats_sen.numSteps = 0; 
      for( _ifct=0; _ifct < _nf; _ifct++ ){
        _cv_flag = CVodeGetB( _cv_mem, _indexB[_ifct], &_t, _Ny[_ifct]);
        if( _check_cv_flag( &_cv_flag, "CVodeGetB", 1) )
          { _END_SEN(); return FATAL; }
        _cv_flag = CVodeGetQuadB( _cv_mem, _indexB[_ifct], &_t, _Nyq[_ifct]);
#ifdef MC__ODEBNDS_SUNDIALS_DEBUG
          for( unsigned iy=0; iy<NV_LENGTH_S(_Ny[_ifct]); iy++ )
            std::cout << "_Ny" << _ifct << "[iy] = " << NV_Ith_S(_Ny[_ifct],iy) << std::endl;
#endif
        if( _check_cv_flag( &_cv_flag, "CVodeGetQuadB", 1) )
          { _END_SEN(); return FATAL; }
#ifdef MC__ODEBNDS_SUNDIALS_DEBUG
          for( unsigned iy=0; iy<NV_LENGTH_S(_Nyq[_ifct]); iy++ )
            std::cout << "_Nyq" << _ifct << "[iy] = " << NV_Ith_S(_Nyq[_ifct],iy) << std::endl;
#endif
        void *cv_memB = CVodeGetAdjCVodeBmem(_cv_mem, _indexB[_ifct] );
        long int nstpB;
        _cv_flag = CVodeGetNumSteps( cv_memB, &nstpB );
        stats_sen.numSteps += nstpB;

        // Add function contribution to adjoint bounds (discontinuities)
        if( _istg > 1  ){
          _pos_fct = ( _vFCT.size()>=ns? _istg-1:0 );
          if( _pos_fct
           && ( !_CC_SET_ASA( _pos_fct, _ifct )
             || !_CC_I_SET( options )
             || !_CC_I_SEN( options, _t, _vec_sta[_istg-1].data(), NV_DATA_S(_Ny[_ifct]) )
             || ( _Nyq && _Nyq[_ifct] && !_CC_I_QUAD_ASA( options, NV_DATA_S(_Nyq[_ifct]) ) ) ) )
            { _END_SEN(); return FATAL; }
          else if( !_pos_fct )
             _GET_I_SEN( options, NV_DATA_S(_Ny[_ifct]), _np, _Nyq && _Nyq[_ifct]? NV_DATA_S(_Nyq[_ifct]): 0 );
#ifdef MC__ODEBNDS_SUNDIALS_DEBUG
          for( unsigned iy=0; iy<NV_LENGTH_S(_Ny[_ifct]); iy++ )
            std::cout << "_Ny" << _ifct << "[iy] = " << NV_Ith_S(_Ny[_ifct],iy) << std::endl;
          for( unsigned iy=0; iy<NV_LENGTH_S(_Nyq[_ifct]); iy++ )
            std::cout << "_Nyq" << _ifct << "[iy] = " << NV_Ith_S(_Nyq[_ifct],iy) << std::endl;
#endif

          // Reset ODE solver - needed in case of discontinuity
          if( !_CC_CVODES_ASA( _ifct, _indexB[_ifct] ) )
            { _END_SEN(); return FATAL; }
        }
        // Add initial state contribution to derivative bounds
        else if( !_IC_I_SET_ASA( options )
              || !_IC_I_SEN( options, _t, _vec_sta[_istg-1].data(), NV_DATA_S(_Ny[_ifct]) )
              || ( _Nyq && _Nyq[_ifct] && !_IC_I_QUAD_ASA( options, NV_DATA_S(_Nyq[_ifct]) ) ) )
          { _END_SEN(); return FATAL; }

        // Keep track of results at stage times
        if( Ilk && !Ilk[_istg-1] ) Ilk[_istg-1] = new T[_nx*_nf];
        for( unsigned iy=0; Ilk[_istg-1] && iy<_nx; iy++ )
          Ilk[_istg-1][_ifct*_nx+iy] = _Iy[iy];
        for( unsigned iq=0; iq<_np; iq++ )
          Ifp[iq*_nf+_ifct] = _Iyq[iq];
      }

      // Display & record adjoint intermediate results
      if( options.DISPLAY >= 1 ){
        _print_interm( tk[_istg-1], _nf*_nx, Ilk[_istg-1], "l", os );
        //_print_interm( _nf*_np, Ifp, "fp", os );
      }
      if( options.RESRECORD )
        results_sen.push_back( Results( tk[_istg-1], _nf*_nx, Ilk[_istg-1] ) );
    }
    if( options.DISPLAY >= 1 )
      _print_interm( _nf*_np, Ifp, "fp", os );
  }
  catch(...){
    _END_SEN();
    if( options.DISPLAY >= 1 ) _print_stats( stats_sen, os );
    return FAILURE;
  }
  _END_SEN();
  if( options.DISPLAY >= 1 ) _print_stats( stats_sen, os );

  return NORMAL;
}

template <typename T, typename PMT, typename PVT> inline bool
ODEBNDS_SUNDIALS<T,PMT,PVT>::_INI_I_FSA
( const unsigned np, const T *Ip )
{
  // Initialize bound propagation
  if( !_INI_I_SEN( options, np, Ip, np, _nq, options.ETOLS ) )
    return false;

  // Set SUNDIALS sensitivity/quadrature arrays
  unsigned cv_Ny_size, cv_Nyq_size, cv_Nfp_size;
  switch( options.WRAPMIT){
  case Options::NONE:
  case Options::DINEQ:
    cv_Ny_size  = 2*_nx;
    cv_Nyq_size = 2*_nq;
    cv_Nfp_size = 2*_nf;
    break;
  case Options::ELLIPS:
  default:
    cv_Ny_size  = _nx*(1+np)+_nx*(_nx+1)/2;
    cv_Nyq_size = _nq*(2+np);
    cv_Nfp_size = _nf*(2+np);
    break;
  }
  if( _nvec != np ){
    if( _Ny )    N_VDestroyVectorArray_Serial( _Ny,  _nvec ); _Ny = 0;
    if( _Nyq )   N_VDestroyVectorArray_Serial( _Nyq, _nvec ); _Nyq = 0;
    if( _Nfp )   N_VDestroyVectorArray_Serial( _Nfp, _nvec ); _Nfp = 0;
    _nvec = np;
    _Ny    = N_VCloneVectorArray_Serial( np, _Nx );
    if( _nq ) _Nyq = N_VCloneVectorArray_Serial( np, _Nx );
    if( _nf ) _Nfp = N_VCloneVectorArray_Serial( np, _Nx );
  }
  for( unsigned i=0; i<np; i++ ){
    if( !_Ny[i] || cv_Ny_size != NV_LENGTH_S( _Ny[i] ) ){
      if( _Ny[i] ) N_VDestroy_Serial( _Ny[i] );
      _Ny[i] = N_VNew_Serial( cv_Ny_size );
    }
    if( _Nyq && (!_Nyq[i] || cv_Nyq_size != NV_LENGTH_S( _Nyq[i]) ) ){
      if( _Nyq[i] ) N_VDestroy_Serial( _Nyq[i] );
      _Nyq[i] = cv_Nyq_size? N_VNew_Serial( cv_Nyq_size ): 0;
    }
    if( _Nfp && (!_Nfp[i] || cv_Nfp_size != NV_LENGTH_S( _Nfp[i]) ) ){
      if( _Nfp[i] ) N_VDestroy_Serial( _Nfp[i] );
      _Nfp[i] = cv_Nfp_size? N_VNew_Serial( cv_Nfp_size ): 0;
    }
  }
  if( !_NTOLy || cv_Ny_size != NV_LENGTH_S( _NTOLy ) )
    if( _NTOLy ) N_VDestroy_Serial( _NTOLy );
      _NTOLy = N_VNew_Serial( cv_Ny_size );

  // Reset result record and statistics
  results_sen.clear();
  _init_stats( stats_sen );

  return true;
}

template <typename T, typename PMT, typename PVT> inline int
ODEBNDS_SUNDIALS<T,PMT,PVT>::MC_CVFSARHSI__
( int Ns, realtype t, N_Vector x, N_Vector xdot, int is, N_Vector y,
  N_Vector ydot, void *user_data, N_Vector tmp1, N_Vector tmp2 )
{
  ODEBNDS_SUNDIALS<T,PMT,PVT> *pODEBNDS = ODEBNDS_SUNDIALS<T,PMT,PVT>::pODEBNDS;
#ifdef MC__ODEBNDS_SUNDIALS_DEBUG
  std::cout << "@t=" << t << "\nx:\n";
  for( unsigned i=0; i<NV_LENGTH_S( x ); i++ ) std::cout << NV_Ith_S( x, i ) << std::endl;
  std::cout << "y:\n";
  for( unsigned i=0; i<NV_LENGTH_S( y ); i++ ) std::cout << NV_Ith_S( y, i ) << std::endl;
#endif
  bool flag = pODEBNDS->_RHS_I_SEN( pODEBNDS->options, t, NV_DATA_S( x ), NV_DATA_S( y ),
    NV_DATA_S( ydot ), is );
#ifdef MC__ODEBNDS_SUNDIALS_DEBUG
  std::cout << "ydot:\n";
  for( unsigned i=0; i<NV_LENGTH_S( ydot ); i++ ) std::cout << NV_Ith_S( ydot, i ) << std::endl;
  { int dum; std::cin >> dum; }
#endif
  ODEBNDS_SUNDIALS<T,PMT,PVT>::pODEBNDS = pODEBNDS;
  pODEBNDS->stats_sen.numRHS++;
  return( (flag && _diam( pODEBNDS->_nx, pODEBNDS->_Iy ) < pODEBNDS->options.DMAX)? 0: -1 );
}

template <typename T, typename PMT, typename PVT> inline int
ODEBNDS_SUNDIALS<T,PMT,PVT>::MC_CVFSAQUADI__
( int Ns, realtype t, N_Vector x, N_Vector *y, N_Vector qdot, N_Vector *qSdot, 
  void *user_data, N_Vector tmp1, N_Vector tmp2 )
{
  ODEBNDS_SUNDIALS<T,PMT,PVT> *pODEBNDS = ODEBNDS_SUNDIALS<T,PMT,PVT>::pODEBNDS;
#ifdef MC__ODEBNDS_SUNDIALS_DEBUG
  std::cout << "@t=" << t << "\nx:\n";
  for( unsigned i=0; i<NV_LENGTH_S( x ); i++ ) std::cout << NV_Ith_S( x, i ) << std::endl;
  std::cout << "qdot:\n";
  for( unsigned i=0; i<NV_LENGTH_S( qdot ); i++ ) std::cout << NV_Ith_S( qdot, i ) << std::endl;
#endif
  bool flag = true;
  for( int is=0; is<Ns && flag; is++ ){
    pODEBNDS->_GET_I_SEN( pODEBNDS->options, NV_DATA_S(x), NV_DATA_S(y[is]), (realtype*)0, 0, (realtype*)0 );
#ifdef MC__ODEBNDS_SUNDIALS_DEBUG
    std::cout << "y:\n";
    for( unsigned i=0; i<NV_LENGTH_S( y[is] ); i++ ) std::cout << NV_Ith_S( y[is], i ) << std::endl;
#endif
    flag = pODEBNDS->_RHS_I_QUAD( pODEBNDS->options, pODEBNDS->_nq, NV_DATA_S( qSdot[is] ), is );
#ifdef MC__ODEBNDS_SUNDIALS_DEBUG
    std::cout << "qSdot:\n";
    for( unsigned i=0; i<NV_LENGTH_S( qSdot[is] ); i++ ) std::cout << NV_Ith_S( qSdot[is], i ) << std::endl;
    { int dum; std::cin >> dum; }
#endif
  }
  ODEBNDS_SUNDIALS<T,PMT,PVT>::pODEBNDS = pODEBNDS;
  return( flag? 0: -1 );
}

//! @fn template <typename T, typename PMT, typename PVT> inline typename ODEBNDS_SUNDIALS<T,PMT,PVT>::STATUS ODEBNDS_SUNDIALS<T,PMT,PVT>::bounds_FSA(
//! const unsigned ns, const double*tk, const T*Ip, T**Ixk, 
//! T*If, T**Ixpk, T*Ifp, std::ostream&os=std::cout )
//!
//! This function computes an enclosure of the reachable set of the parametric ODEs
//! using continuous-time propagation of convex sets (intervals, ellipsoids),
//! together with forward sensitivity bounds:
//!  - <a>ns</a>    [input]  number of time stages
//!  - <a>tk</a>    [input]  stage times, including the initial time
//!  - <a>Ip</a>   [input]   interval enclosure of parameter set
//!  - <a>Ixk</a>  [output]  interval enclosure of state variables at stage times
//!  - <a>If</a>   [output]  interval enclosure of state functions
//!  - <a>Ixpk</a> [output]  interval enclosure of state sensitivity variables at stage times
//!  - <a>Ifp</a>  [output]  interval enclosure of state function derivatives
//!  - <a>os</a>    [input]  output stream [default: std::cout]
//! .
//! The return value is the status.
template <typename T, typename PMT, typename PVT>
inline typename ODEBNDS_SUNDIALS<T,PMT,PVT>::STATUS
ODEBNDS_SUNDIALS<T,PMT,PVT>::bounds_FSA
( const unsigned ns, const double*tk, const T*Ip, T**Ixk, T*If,
  T**Ixpk, T*Ifp, std::ostream&os )
{
  if( !_np ) return bounds( ns, tk, Ip, Ixk, If, os );

  // Check arguments
  if( !tk || (_np && !Ip) || !Ixk || !Ixpk || ( _nf && (!If || !Ifp) ) ) return FATAL;

  try{
    // Initialize trajectory integration
    if( !ODEBND_SUNDIALS<T,PMT,PVT>::_INI_I_STA( _np, Ip ) 
     || !_INI_I_FSA( _np, Ip ) ) return FATAL;
    _t = tk[0];

    // Bounds on initial states/quadratures
    if( !ODEBND_BASE<T,PMT,PVT>::_IC_I_SET( options )
     || !ODEBND_BASE<T,PMT,PVT>::_IC_I_STA( options, _t, NV_DATA_S( _Nx ) )
     || ( _Nq && !ODEBND_BASE<T,PMT,PVT>::_IC_I_QUAD( options, NV_DATA_S( _Nq ) ) ) )
      { _END_STA(); _END_SEN(); return FATAL; }
    if( Ixk && !Ixk[0] ) Ixk[0] = new T[_nx];
    for( unsigned ix=0; Ixk[0] && ix<_nx; ix++ ) Ixk[0][ix] = _Ix[ix];

    // Bounds on initial state/quadrature sensitivities
    if( Ixpk && !Ixpk[0] ) Ixpk[0] = new T[_nx*_np];
    for( _isen=0; _isen<_np; _isen++ ){
      if( !_IC_SET_FSA( _isen )
       || !ODEBND_BASE<T,PMT,PVT>::_IC_I_STA( options, _t, NV_DATA_S(_Ny[_isen]) )
       || ( _Nyq && _Nyq[_isen]
         && !ODEBND_BASE<T,PMT,PVT>::_IC_I_QUAD( options, NV_DATA_S(_Nyq[_isen]) ) ) ) 
        { _END_STA(); _END_SEN(); return FATAL; }
      for( unsigned iy=0; Ixpk[0] && iy<_nx; iy++ ) Ixpk[0][_isen*_nx+iy] = _Ix[iy];
#ifdef MC__ODEBNDS_SUNDIALS_DEBUG
      std::cout << "qS[" << _isen << "]: " << _Nyq[_isen] << "(" << NV_LENGTH_S( _Nyq[_isen] ) << ")\n";
      for( unsigned i=0; i<NV_LENGTH_S( _Nyq[_isen] ); i++ ) std::cout << NV_Ith_S( _Nyq[_isen], i ) << std::endl;
      { int dum; std::cin >> dum; }
#endif
    }

    // Display & record initial results
    if( options.DISPLAY >= 1 ){
      _print_interm( _t, _nx, Ixk[0], "x", os );
      _print_interm( _nx*_np, Ixpk[0], "xp", os );
    }
    if( options.RESRECORD ){
      results_sta.push_back( Results( _t, _nx, Ixk[0] ) );
      results_sen.push_back( Results( _t, _nx*_np, Ixpk[0] ) );
    }

    // Integrate ODEs through each stage using SUNDIALS
    ODEBND_SUNDIALS<T,PMT,PVT>::pODEBND = pODEBNDS = this;
    if( !ODEBND_SUNDIALS<T,PMT,PVT>::_INI_CVODE
          ( ODEBND_SUNDIALS<T,PMT,PVT>::MC_CVRHSI__,
            ODEBND_SUNDIALS<T,PMT,PVT>::MC_CVQUADI__ )
     || !_INI_CVODES( MC_CVFSARHSI__, MC_CVFSAQUADI__ ) )
      { _END_STA(); _END_SEN(); return FATAL; }

    for( _istg=0; _istg<ns; _istg++ ){
      // Bounds on state discontinuities (if any) at stage times
      // and integrator reinitialization (if applicable)
      _pos_ic = ( _vIC.size()>=ns? _istg:0 );
      if( _pos_ic
       && ( !ODEBND_BASE<T,PMT,PVT>::_CC_I_SET( options, _pos_ic )
         || !ODEBND_BASE<T,PMT,PVT>::_CC_I_STA( options, _t, NV_DATA_S( _Nx ) )
         || !ODEBND_SUNDIALS<T,PMT,PVT>::_CC_CVODE_STA() ) )
        { _END_STA(); _END_SEN(); return FATAL; }
      if( _istg 
       && ( (_Nq && !ODEBND_BASE<T,PMT,PVT>::_IC_I_QUAD( options, NV_DATA_S( _Nq ) ) ) // quadrature reinitialization
         || !ODEBND_SUNDIALS<T,PMT,PVT>::_CC_CVODE_QUAD() ) )
        { _END_STA(); _END_SEN(); return FATAL; }
      for( _isen=0; _isen<_np; _isen++ ){
        if( _pos_ic
         && ( !_CC_SET_FSA( _pos_ic, _isen )
           || !_CC_I_SET( options )
           || !_CC_I_SEN( options, _t, NV_DATA_S( _Nx ), NV_DATA_S(_Ny[_isen]) )
           || ( !_isen && !_CC_CVODES_FSA() ) ) )
          { _END_STA(); _END_SEN(); return FATAL; }
        if( _istg
         && ( ( _Nyq && _Nyq[_isen] && !ODEBND_BASE<T,PMT,PVT>::_IC_I_QUAD( options, NV_DATA_S(_Nyq[_isen]) ) ) //quadrature sensitivity reinitialization
           || ( !_isen && !_CC_CVODES_QUAD() ) ) )
            { _END_STA(); _END_SEN(); return FATAL; }
#ifdef MC__ODEBNDS_SUNDIALS_DEBUG
        for( unsigned iy=0; iy<NV_LENGTH_S(_Ny[_isen]); iy++ )
          std::cout << "_Ny" << _isen << "[iy] = " << NV_Ith_S(_Ny[_isen],iy) << std::endl;
        for( unsigned iy=0; iy<NV_LENGTH_S(_Nyq[_isen]); iy++ )
          std::cout << "_Nyq" << _isen << "[iy] = " << NV_Ith_S(_Nyq[_isen],iy) << std::endl;
#endif
      }

      // update list of operations in RHS, JAC, QUAD, RHSFSA and QUADFSA
      _pos_rhs  = ( _vRHS.size()<=1?  0: _istg );
      _pos_quad = ( _vQUAD.size()<=1? 0: _istg );
      if( (!_istg || _pos_rhs || _pos_quad)
        && ( !ODEBND_BASE<T,PMT,PVT>::_RHS_I_SET( options, _pos_rhs, _pos_quad )
          || !_RHS_SET_FSA( _pos_rhs, _pos_quad )
          || !_RHS_I_SET( options, _np, _nq ) ) )
        { _END_STA(); _END_SEN(); return FATAL; }

      // integrate till end of time stage
      _cv_flag = CVodeSetStopTime( _cv_mem, tk[_istg+1] );
      if( _check_cv_flag(&_cv_flag, "CVodeSetStopTime", 1) )
        { _END_STA(); return FATAL; }
      while( _t < tk[_istg+1] ){
        _cv_flag = CVode( _cv_mem, tk[_istg+1], _Nx, &_t, CV_ONE_STEP );
        if( _check_cv_flag(&_cv_flag, "CVode", 1)
         || (options.NMAX && stats_sta.numSteps > options.NMAX)
         || _diam( _nx, _Ix  ) > options.DMAX )
          throw Exceptions( Exceptions::INTERN );
        stats_sta.numSteps++;
        stats_sen.numSteps++;
      }

      // Bounds on intermediate states, quadratures and their sensitivities
      if( _nq ){
        for( unsigned ip=0; ip<_np; ip++ ){
          _cv_flag = CVodeGetQuadSens1(_cv_mem, &_t, ip, _Nyq[ip]);
          if( _check_cv_flag( &_cv_flag, "CVodeGetQuadSens", 1) )
            { _END_STA(); _END_SEN(); return FATAL; }
#ifdef MC__ODEBNDS_SUNDIALS_DEBUG
          std::cout << "qS[" << ip << "]: " << _Nyq[ip] << "(" << NV_LENGTH_S( _Nyq[ip] ) << ")\n";
          for( unsigned i=0; i<NV_LENGTH_S( _Nyq[ip] ); i++ ) std::cout << NV_Ith_S( _Nyq[ip], i ) << std::endl;
          { int dum; std::cin >> dum; }
#endif
        }
        _cv_flag = CVodeGetQuad( _cv_mem, &_t, _Nq );
        if( _check_cv_flag(&_cv_flag, "CVodeGetQuad", 1) )
          { _END_STA(); _END_SEN(); return FATAL; }
      }
      _cv_flag = CVodeGetSens(_cv_mem, &_t, _Ny);
      if( _check_cv_flag( &_cv_flag, "CVodeGetSens", 1) )
         { _END_STA(); _END_SEN(); return FATAL; }

      // Add intermediate function terms and derivatives
      _pos_fct = ( _vFCT.size()>=ns? _istg:0 );
      ODEBND_SUNDIALS<T,PMT,PVT>::_GET_I_STA( options, NV_DATA_S(_Nx), _nq && _Nq? NV_DATA_S(_Nq): 0 );
      if( (_vFCT.size()>=ns || _istg==ns-1)
       && !ODEBND_BASE<T,PMT,PVT>::_FCT_I_STA( options, _pos_fct, _t, NV_DATA_S( _Nf ), If ) )
        { _END_STA(); _END_SEN(); return FATAL; }
      if( Ixk && !Ixk[_istg+1] ) Ixk[_istg+1] = new T[_nx];
      for( unsigned ix=0; ix<_nx; ix++ ) Ixk[_istg+1][ix] = _Ix[ix];
      for( _isen=0; _isen<_np; _isen++ ){
        _GET_I_SEN( options, NV_DATA_S(_Nx), NV_DATA_S(_Ny[_isen]), _nq && _Nq? NV_DATA_S(_Nq): 0,
                    _nq, _nq && _Nyq[_isen]? NV_DATA_S(_Nyq[_isen]): 0 );
        if( (_vFCT.size()>=ns || _istg==ns-1)
         && !ODEBND_BASE<T,PMT,PVT>::_FCT_I_SEN( options, _pos_fct, _isen, _t,
                    NV_DATA_S( _Nfp[_isen] ), Ifp+_isen*_nf ) )
          { _END_STA(); _END_SEN(); return FATAL; }
        if( Ixpk && !Ixpk[_istg+1] ) Ixpk[_istg+1] = new T[_nx*_np];
        for( unsigned iy=0; iy<_nx; iy++ ) Ixpk[_istg+1][_isen*_nx+iy] = _Iy[iy];
      }

      // Display & record stage results
      if( options.DISPLAY >= 1 ){
        _print_interm( _t, _nx, Ixk[_istg+1], "x", os );
        _print_interm( _nx*_np, Ixpk[_istg+1], "xp", os );
      }
      if( options.RESRECORD ){
        results_sta.push_back( Results( tk[_istg+1], _nx, Ixk[_istg+1] ) );
        results_sen.push_back( Results( tk[_istg+1], _nx*_np, Ixpk[_istg+1] ) );
      }
    }

    // Bounds on final quadratures and functions
    if( options.DISPLAY >= 1 ){
      _print_interm( _nf, If, "f", os );
      _print_interm( _nf*_np, Ifp, "fp", os );
    }
  }
  catch(...){
    _END_STA();
    if( options.DISPLAY >= 1 ){
      os << " ABORT TIME  " << std::scientific << std::left
                            << std::setprecision(5) << _t << std::endl;
      _print_stats( stats_sta, os );
      _print_stats( stats_sen, os );
    }
    return FAILURE;
  }

  _END_STA();
  if( options.DISPLAY >= 1 ) _print_stats( stats_sta, os );
  return NORMAL;
}

template <typename T, typename PMT, typename PVT> inline bool
ODEBNDS_SUNDIALS<T,PMT,PVT>::_INI_PM_ASA
( const unsigned np, const PVT *PMp )
{
  // Initialize bound propagation
  if( !_INI_PM_SEN( options, np, PMp, _nf, np, options.ETOLB ) )
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
  if( _nvec != _nf ){
    if( _Ny )   N_VDestroyVectorArray_Serial( _Ny,  _nvec );
    if( _Nyq )  N_VDestroyVectorArray_Serial( _Nyq, _nvec );
    _nvec = _nf;
    _Ny  = N_VCloneVectorArray_Serial( _nvec, _Nx );
    _Nyq = N_VCloneVectorArray_Serial( _nvec, _Nx );
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
  if( !_NTOLy || cv_Ny_size != NV_LENGTH_S( _NTOLy ) )
    if( _NTOLy ) N_VDestroy_Serial( _NTOLy );
      _NTOLy = N_VNew_Serial( cv_Ny_size );

  // Reset result record and statistics
  results_sen.clear();
  _init_stats( stats_sen );

  return true;
}

template <typename T, typename PMT, typename PVT> inline int
ODEBNDS_SUNDIALS<T,PMT,PVT>::MC_CVASARHSPM__
( realtype t, N_Vector x, N_Vector y, N_Vector ydot, void *user_data )
{
  ODEBNDS_SUNDIALS<T,PMT,PVT> *pODEBNDS = ODEBNDS_SUNDIALS<T,PMT,PVT>::pODEBNDS;
  pODEBNDS->_ifct = *static_cast<unsigned*>( user_data );
  int flag = pODEBNDS->_RHS_PM_SEN( pODEBNDS->options, t, NV_DATA_S( x ), NV_DATA_S( y ),
    NV_DATA_S( ydot ), pODEBNDS->_ifct, true );
  ODEBNDS_SUNDIALS<T,PMT,PVT>::pODEBNDS = pODEBNDS;
  pODEBNDS->stats_sen.numRHS++;
  if( _diam( pODEBNDS->_nx, pODEBNDS->_PMy ) > pODEBNDS->options.DMAX ) return -1;
  return flag;
}

template <typename T, typename PMT, typename PVT> inline int
ODEBNDS_SUNDIALS<T,PMT,PVT>::MC_CVASAQUADPM__
( realtype t, N_Vector x, N_Vector y, N_Vector qdot, void *user_data )
{
  ODEBNDS_SUNDIALS<T,PMT,PVT> *pODEBNDS = ODEBNDS_SUNDIALS<T,PMT,PVT>::pODEBNDS;
  pODEBNDS->_ifct = *static_cast<unsigned*>( user_data );
  bool flag = pODEBNDS->_RHS_PM_QUAD( pODEBNDS->options, pODEBNDS->_np,
    NV_DATA_S( qdot ), pODEBNDS->_ifct, true );
  ODEBNDS_SUNDIALS<T,PMT,PVT>::pODEBNDS = pODEBNDS;
  return( flag? 0: -1 );
}

//! @fn template <typename T, typename PMT, typename PVT> inline typename ODEBNDS_SUNDIALS<T,PMT,PVT>::STATUS ODEBNDS_SUNDIALS<T,PMT,PVT>::bounds_ASA(
//! const unsigned ns, const double*tk, const PVT*PMp, PVT**PMxk, 
//! PVT*PMf, PVT**PMlk, PVT*PMfp, std::ostream&os=std::cout )
//!
//! This function computes an enclosure of the reachable set of the parametric ODEs
//! using continuous-time propagation of polynomial models with convex remainders
//! (intervals, ellipsoids), together with forward sensitivity bounds:
//!  - <a>ns</a>    [input]  number of time stages
//!  - <a>tk</a>    [input]  stage times, including the initial time
//!  - <a>PMp</a>   [input]  polynomial model of parameter set
//!  - <a>PMxk</a>  [output] polynomial model of state variables at stage times
//!  - <a>PMf</a>   [output] polynomial model of state functions
//!  - <a>PMlk</a>  [output] polynomial model of adjoint variables at stage times
//!  - <a>PMfp</a>  [output] polynomial model of state function derivatives
//!  - <a>os</a>    [input]  output stream [default: std::cout]
//! .
//! The return value is the status.
template <typename T, typename PMT, typename PVT>
inline typename ODEBNDS_SUNDIALS<T,PMT,PVT>::STATUS
ODEBNDS_SUNDIALS<T,PMT,PVT>::bounds_ASA
( const unsigned ns, const double*tk, const PVT*PMp, PVT**PMxk, PVT*PMf,
  PVT**PMlk, PVT*PMfp, std::ostream&os )
{
  return bounds_ASA( ns, tk, PMp, PMxk, 0, PMf, PMlk, 0, PMfp, os );
}

//! @fn template <typename T, typename PMT, typename PVT> inline typename ODEBNDS_SUNDIALS<T,PMT,PVT>::STATUS ODEBNDS_SUNDIALS<T,PMT,PVT>::bounds_ASA(
//! const unsigned ns, const double*tk, const PVT*PMp, PVT**PMxk, E*ERxk, 
//! PVT*PMf, PVT**PMlk, E**ERlk, PVT*PMfp, std::ostream&os=std::cout )
//!
//! This function computes an enclosure of the reachable set of the parametric ODEs
//! using continuous-time propagation of polynomial models with convex remainders
//! (intervals, ellipsoids), together with adjoint sensitivity bounds:
//!  - <a>ns</a>    [input]  number of time stages
//!  - <a>tk</a>    [input]  stage times, including the initial time
//!  - <a>PMp</a>   [input]  polynomial model of parameter set
//!  - <a>PMxk</a>  [output] polynomial model of state variables at stage times
//!  - <a>ERxk</a>  [output] ellipsoidal remainder of state variables at stage times (only if ellipsoidal bounder is selected)
//!  - <a>PMf</a>   [output] polynomial model of state functions
//!  - <a>PMlk</a>  [output] polynomial model of adjoint variables at stage times
//!  - <a>ERxpk</a> [output] ellipsoidal remainder of adjoint variables at stage times (only if ellipsoidal bounder is selected)
//!  - <a>PMfp</a>  [output] polynomial model of state function derivatives
//!  - <a>os</a>    [input]  output stream [default: std::cout]
//! .
//! The return value is the status.
template <typename T, typename PMT, typename PVT>
inline typename ODEBNDS_SUNDIALS<T,PMT,PVT>::STATUS
ODEBNDS_SUNDIALS<T,PMT,PVT>::bounds_ASA
( const unsigned ns, const double*tk, const PVT*PMp, PVT**PMxk, E*ERxk,
  PVT*PMf, PVT**PMlk, E**ERlk, PVT*PMfp, std::ostream&os )
{
  // Compute state bounds and store intermediate results
  STATUS flag = NORMAL;
  flag = ODEBND_SUNDIALS<T,PMT,PVT>::_bounds( ns, tk, PMp, PMxk, ERxk, PMf, true, os);
  if( flag != NORMAL ) return flag;

  // Nothing to do if no functions are defined
  if( !_nf || !_np ) return NORMAL;

  // Check size
  if( !PMlk || (_nf && !PMfp) ) return FATAL;

  try{
    // Initialize adjoint bound integration
    if( !_INI_PM_ASA( _np, PMp )) return FATAL;
    _t = tk[ns];

    // Bounds on terminal adjoints/quadratures
    if( PMlk && !PMlk[ns] ) PMlk[ns] = new PVT[_nx*_nf];
    _pos_fct = ( _vFCT.size()>=ns? ns-1:0 );
    for( _ifct=0; _ifct < _nf; _ifct++ ){
      if( !_TC_PM_SET_ASA( options, _pos_fct, _ifct )
       || !_TC_PM_SEN( options, _t, _vec_sta[ns].data(), NV_DATA_S(_Ny[_ifct]) )
       || ( _Nyq && _Nyq[_ifct] && !_TC_PM_QUAD_ASA( options, NV_DATA_S(_Nyq[_ifct]) ) ) )
        { _END_SEN(); return FATAL; }
      for( unsigned iy=0; PMlk[ns] && iy<_nx; iy++ )
        PMlk[ns][_ifct*_nx+iy] = _PMy[iy];
      if( options.WRAPMIT == Options::ELLIPS && ERlk && ERlk[0] ) ERlk[0][_ifct] = _Edy;
      for( unsigned iq=0; iq<_np; iq++ )
        PMfp[iq*_nf+_ifct] = _PMyq[iq];
    }
    // Display & record adjoint terminal results
    if( options.DISPLAY >= 1 )
      _print_interm( tk[ns], _nx*_nf, PMlk[ns], "l", os );
    if( options.RESRECORD )
      results_sen.push_back( Results( tk[ns], _nf*_nx, PMlk[ns] ) );

    // Initialization of adjoint integration
    for( _ifct=0; _ifct < _nf; _ifct++ )
      if( !_INI_CVODES( MC_CVASARHSPM__, MC_CVASAQUADPM__, _ifct,
        _indexB[_ifct], _iusrB[_ifct] ) )
        { _END_SEN(); return FATAL;}

    // Integrate adjoint ODEs through each stage using SUNDIALS
    pODEBNDS = this;
    for( _istg=ns; _istg>0; _istg-- ){

      // Update list of operations in RHSADJ and QUADADJ
      _pos_rhs  = ( _vRHS.size() <=1? 0: _istg-1 );
      _pos_quad = ( _vQUAD.size()<=1? 0: _istg-1 );
      if( !_RHS_SET_ASA( _pos_rhs, _pos_quad, _pos_fct ) 
       || !_RHS_PM_SET( options, _nf, _np ) )
        { _END_SEN(); return FATAL; }

      // Propagate bounds backward to previous stage time
      _cv_flag = CVodeB( _cv_mem, tk[_istg-1], CV_NORMAL );
      if( _check_cv_flag(&_cv_flag, "CVodeB", 1) )
        { _END_SEN(); return FATAL; }
      _t = tk[_istg-1];

      // Bounds on states/adjoints/quadratures at stage time
      //stats_sen.numSteps = 0; 
      for( _ifct=0; _ifct < _nf; _ifct++ ){
        _cv_flag = CVodeGetB( _cv_mem, _indexB[_ifct], &_t, _Ny[_ifct]);
        if( _check_cv_flag( &_cv_flag, "CVodeGetB", 1) )
          { _END_SEN(); return FATAL; }
#ifdef MC__ODEBNDS_SUNDIALS_DEBUG
        for( unsigned iy=0; iy<NV_LENGTH_S(_Ny[_ifct]); iy++ )
          std::cout << "_Ny" << _ifct << "[iy] = " << NV_Ith_S(_Ny[_ifct],iy) << std::endl;
#endif
        _cv_flag = CVodeGetQuadB( _cv_mem, _indexB[_ifct], &_t, _Nyq[_ifct]);
        if( _check_cv_flag( &_cv_flag, "CVodeGetQuadB", 1) )
          { _END_SEN(); return FATAL; }
#ifdef MC__ODEBNDS_SUNDIALS_DEBUG
        for( unsigned iy=0; iy<NV_LENGTH_S(_Nyq[_ifct]); iy++ )
          std::cout << "_Nyq" << _ifct << "[iy] = " << NV_Ith_S(_Nyq[_ifct],iy) << std::endl;
#endif
        void *cv_memB = CVodeGetAdjCVodeBmem(_cv_mem, _indexB[_ifct] );
        long int nstpB;
        _cv_flag = CVodeGetNumSteps( cv_memB, &nstpB );
        stats_sen.numSteps += nstpB;
        //{ int dum; std::cin >> dum; }

        // Add function contribution to adjoint bounds (discontinuities)
        if( _istg > 1  ){
          _pos_fct = ( _vFCT.size()>=ns? _istg-1:0 );
          if( _pos_fct 
           && ( !_CC_SET_ASA( _pos_fct, _ifct )
             || !_CC_PM_SET( options )
             || !_CC_PM_SEN( options, _t, _vec_sta[_istg-1].data(), NV_DATA_S(_Ny[_ifct]) )
             || ( _Nyq && _Nyq[_ifct] && !_CC_PM_QUAD_ASA( options, NV_DATA_S(_Nyq[_ifct]) ) ) ) )
            { _END_SEN(); return FATAL; }
          else if( !_pos_fct )
             _GET_PM_SEN( options, NV_DATA_S(_Ny[_ifct]), _np, _Nyq && _Nyq[_ifct]? NV_DATA_S(_Nyq[_ifct]): 0 );
#ifdef MC__ODEBNDS_SUNDIALS_DEBUG
          for( unsigned iy=0; iy<NV_LENGTH_S(_Ny[_ifct]); iy++ )
            std::cout << "_Ny" << _ifct << "[iy] = " << NV_Ith_S(_Ny[_ifct],iy) << std::endl;
          for( unsigned iy=0; iy<NV_LENGTH_S(_Nyq[_ifct]); iy++ )
            std::cout << "_Nyq" << _ifct << "[iy] = " << NV_Ith_S(_Nyq[_ifct],iy) << std::endl;
#endif
          // Reset ODE solver - needed in case of discontinuity
          if( !_CC_CVODES_ASA( _ifct, _indexB[_ifct] ) )
            { _END_SEN(); return FATAL; }
        }
        // Add initial state contribution to derivative bounds
        else if( !_IC_PM_SET_ASA( options )
              || !_IC_PM_SEN( options, _t, _vec_sta[_istg-1].data(), NV_DATA_S(_Ny[_ifct]) )
              || ( _Nyq && _Nyq[_ifct] && !_IC_PM_QUAD_ASA( options, NV_DATA_S(_Nyq[_ifct]) ) ) )
          { _END_SEN(); return FATAL; }

        // Keep track of results at stage times
        if( PMlk && !PMlk[_istg-1] ) PMlk[_istg-1] = new PVT[_nx*_nf];
        for( unsigned iy=0; PMlk[_istg-1] && iy<_nx; iy++)
          PMlk[_istg-1][_ifct*_nx+iy] =_PMy[iy];
        if( options.WRAPMIT == Options::ELLIPS && ERlk && ERlk[_istg-1] ) ERlk[_istg-1][_ifct] = _Edy;
        for( unsigned iq=0; iq<_np; iq++ )
          PMfp[iq*_nf+_ifct] = _PMyq[iq];
      }

      // Display & record adjoint intermediate results
      if( options.DISPLAY >= 1 )
        _print_interm( tk[_istg-1], _nf*_nx, PMlk[_istg-1], "l", os );
      if( options.RESRECORD )
        results_sen.push_back( Results( tk[_istg-1], _nf*_nx, PMlk[_istg-1] ) );
      //{ int dum; std::cin >> dum; }
    }

    if( options.DISPLAY >= 1 )
      _print_interm( _nf*_np, PMfp, "fp", os );
  }
  catch(...){
    _END_SEN();
    if( options.DISPLAY >= 1 ) _print_stats( stats_sen, os );
    return FAILURE;
  }

  _END_SEN();
  if( options.DISPLAY >= 1 ) {_print_stats( stats_sen, os );}

  return NORMAL;
}

template <typename T, typename PMT, typename PVT> inline bool
ODEBNDS_SUNDIALS<T,PMT,PVT>::_INI_PM_FSA
( const unsigned np, const PVT *PMp )
{
  // Initialize bound propagation
  if( !_INI_PM_SEN( options, np, PMp, np, _nq, options.ETOLS ) )
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
  unsigned cv_Nyq_size = (_PMenv->nmon()+1)*_nq,
           cv_Nfp_size = (_PMenv->nmon()+1)*_nf;
  if( _nvec != np ){
    if( _Ny )    N_VDestroyVectorArray_Serial( _Ny,  _nvec ); _Ny = 0;
    if( _Nyq )   N_VDestroyVectorArray_Serial( _Nyq, _nvec ); _Nyq = 0;
    if( _Nfp )   N_VDestroyVectorArray_Serial( _Nfp, _nvec ); _Nfp = 0;
    _nvec = np;
    _Ny    = N_VCloneVectorArray_Serial( np, _Nx );
    if( _nq ) _Nyq = N_VCloneVectorArray_Serial( np, _Nx );
    if( _nf ) _Nfp = N_VCloneVectorArray_Serial( np, _Nx );
  }
  for( unsigned i=0; i<np; i++ ){
    if( !_Ny[i] || cv_Ny_size != NV_LENGTH_S( _Ny[i] ) ){
      if( _Ny[i] ) N_VDestroy_Serial( _Ny[i] );
      _Ny[i] = N_VNew_Serial( cv_Ny_size );
    }
    if( _Nyq && (!_Nyq[i] || cv_Nyq_size != NV_LENGTH_S( _Nyq[i]) ) ){
      if( _Nyq[i] ) N_VDestroy_Serial( _Nyq[i] );
      _Nyq[i] = cv_Nyq_size? N_VNew_Serial( cv_Nyq_size ): 0;
    }
    if( _Nfp && (!_Nfp[i] || cv_Nfp_size != NV_LENGTH_S( _Nfp[i]) ) ){
      if( _Nfp[i] ) N_VDestroy_Serial( _Nfp[i] );
      _Nfp[i] = cv_Nfp_size? N_VNew_Serial( cv_Nfp_size ): 0;
    }
  }
  if( !_NTOLy || cv_Ny_size != NV_LENGTH_S( _NTOLy ) )
    if( _NTOLy ) N_VDestroy_Serial( _NTOLy );
      _NTOLy = N_VNew_Serial( cv_Ny_size );

  // Reset result record and statistics
  results_sen.clear();
  _init_stats( stats_sen );

  return true;
}

template <typename T, typename PMT, typename PVT> inline int
ODEBNDS_SUNDIALS<T,PMT,PVT>::MC_CVFSARHSPM__
( int Ns, realtype t, N_Vector x, N_Vector xdot, int is, N_Vector y,
  N_Vector ydot, void *user_data, N_Vector tmp1, N_Vector tmp2 )
{
  ODEBNDS_SUNDIALS<T,PMT,PVT> *pODEBNDS = ODEBNDS_SUNDIALS<T,PMT,PVT>::pODEBNDS;
#ifdef MC__ODEBNDS_SUNDIALS_DEBUG
  std::cout << "@t=" << t << "\nx:\n";
  for( unsigned i=0; i<NV_LENGTH_S( x ); i++ ) std::cout << NV_Ith_S( x, i ) << std::endl;
  std::cout << "y:\n";
  for( unsigned i=0; i<NV_LENGTH_S( y ); i++ ) std::cout << NV_Ith_S( y, i ) << std::endl;
#endif
  int flag = pODEBNDS->_RHS_PM_SEN( pODEBNDS->options, t, NV_DATA_S( x ), NV_DATA_S( y ),
    NV_DATA_S( ydot ), is );
#ifdef MC__ODEBNDS_SUNDIALS_DEBUG
  std::cout << "ydot:\n";
  for( unsigned i=0; i<NV_LENGTH_S( ydot ); i++ ) std::cout << NV_Ith_S( ydot, i ) << std::endl;
  { int dum; std::cin >> dum; }
#endif
  ODEBNDS_SUNDIALS<T,PMT,PVT>::pODEBNDS = pODEBNDS;
  pODEBNDS->stats_sen.numRHS++;
  if( _diam( pODEBNDS->_nx, pODEBNDS->_PMy ) > pODEBNDS->options.DMAX ) return -1;
  return flag;
}

template <typename T, typename PMT, typename PVT> inline int
ODEBNDS_SUNDIALS<T,PMT,PVT>::MC_CVFSAQUADPM__
( int Ns, realtype t, N_Vector x, N_Vector *y, N_Vector qdot, N_Vector *qSdot, 
  void *user_data, N_Vector tmp1, N_Vector tmp2 )
{
  ODEBNDS_SUNDIALS<T,PMT,PVT> *pODEBNDS = ODEBNDS_SUNDIALS<T,PMT,PVT>::pODEBNDS;
#ifdef MC__ODEBNDS_SUNDIALS_DEBUG
  std::cout << "@t=" << t << "\nx:\n";
  for( unsigned i=0; i<NV_LENGTH_S( x ); i++ ) std::cout << NV_Ith_S( x, i ) << std::endl;
  std::cout << "qdot:\n";
  for( unsigned i=0; i<NV_LENGTH_S( qdot ); i++ ) std::cout << NV_Ith_S( qdot, i ) << std::endl;
#endif
  bool flag = true;
  for( int is=0; is<Ns && flag; is++ ){
    pODEBNDS->_GET_PM_SEN( pODEBNDS->options, NV_DATA_S(x), NV_DATA_S(y[is]), (realtype*)0, 0, (realtype*)0 );
#ifdef MC__ODEBNDS_SUNDIALS_DEBUG
    std::cout << "y:\n";
    for( unsigned i=0; i<NV_LENGTH_S( y[is] ); i++ ) std::cout << NV_Ith_S( y[is], i ) << std::endl;
#endif
    flag = pODEBNDS->_RHS_PM_QUAD( pODEBNDS->options, pODEBNDS->_nq, NV_DATA_S( qSdot[is] ), is );
#ifdef MC__ODEBNDS_SUNDIALS_DEBUG
    std::cout << "qSdot:\n";
    for( unsigned i=0; i<NV_LENGTH_S( qSdot[is] ); i++ ) std::cout << NV_Ith_S( qSdot[is], i ) << std::endl;
    { int dum; std::cin >> dum; }
#endif
  }
  ODEBNDS_SUNDIALS<T,PMT,PVT>::pODEBNDS = pODEBNDS;
  return( flag? 0: -1 );
}

//! @fn template <typename T, typename PMT, typename PVT> inline typename ODEBNDS_SUNDIALS<T,PMT,PVT>::STATUS ODEBNDS_SUNDIALS<T,PMT,PVT>::bounds_FSA(
//! const unsigned ns, const double*tk, const PVT*PMp, PVT**PMxk, 
//! PVT*PMf, PVT**PMxpk, PVT*PMfp, std::ostream&os=std::cout )
//!
//! This function computes an enclosure of the reachable set of the parametric ODEs
//! using continuous-time propagation of polynomial models with convex remainders
//! (intervals, ellipsoids), together with forward sensitivity bounds:
//!  - <a>ns</a>    [input]  number of time stages
//!  - <a>tk</a>    [input]  stage times, including the initial time
//!  - <a>PMp</a>   [input]  polynomial model of parameter set
//!  - <a>PMxk</a>  [output] polynomial model of state variables at stage times
//!  - <a>PMf</a>   [output] polynomial model of state functions
//!  - <a>PMxpk</a> [output] polynomial model of state sensitivity variables at stage times
//!  - <a>PMfp</a>  [output] polynomial model of state function derivatives
//!  - <a>os</a>    [input]  output stream [default: std::cout]
//! .
//! The return value is the status.
template <typename T, typename PMT, typename PVT>
inline typename ODEBNDS_SUNDIALS<T,PMT,PVT>::STATUS
ODEBNDS_SUNDIALS<T,PMT,PVT>::bounds_FSA
( const unsigned ns, const double*tk, const PVT*PMp, PVT**PMxk, PVT*PMf,
  PVT**PMxpk, PVT*PMfp, std::ostream&os )
{
  return bounds_FSA( ns, tk, PMp, PMxk, 0, PMf, PMxpk, 0, PMfp, os );
}

//! @fn template <typename T, typename PMT, typename PVT> inline typename ODEBNDS_SUNDIALS<T,PMT,PVT>::STATUS ODEBNDS_SUNDIALS<T,PMT,PVT>::bounds_FSA(
//! const unsigned ns, const double*tk, const PVT*PMp, PVT**PMxk, 
//! PVT*PMf, PVT**PMxpk, PVT*PMfp, std::ostream&os=std::cout )
//!
//! This function computes an enclosure of the reachable set of the parametric ODEs
//! using continuous-time propagation of polynomial models with convex remainders
//! (intervals, ellipsoids), together with forward sensitivity bounds:
//!  - <a>ns</a>    [input]  number of time stages
//!  - <a>tk</a>    [input]  stage times, including the initial time
//!  - <a>PMp</a>   [input]  polynomial model of parameter set
//!  - <a>PMxk</a>  [output] polynomial model of state variables at stage times
//!  - <a>ERxk</a>  [output] ellipsoidal remainder of state variables at stage times (only if ellipsoidal bounder is selected)
//!  - <a>PMf</a>   [output] polynomial model of state functions
//!  - <a>PMxpk</a> [output] polynomial model of state sensitivity variables at stage times
//!  - <a>ERxpk</a> [output] ellipsoidal remainder of state sensitivity variables at stage times (only if ellipsoidal bounder is selected)
//!  - <a>PMfp</a>  [output] polynomial model of state function derivatives
//!  - <a>os</a>    [input]  output stream [default: std::cout]
//! .
//! The return value is the status.
template <typename T, typename PMT, typename PVT>
inline typename ODEBNDS_SUNDIALS<T,PMT,PVT>::STATUS
ODEBNDS_SUNDIALS<T,PMT,PVT>::bounds_FSA
( const unsigned ns, const double*tk, const PVT*PMp, PVT**PMxk, E*ERxk,
  PVT*PMf, PVT**PMxpk, E**ERxpk, PVT*PMfp, std::ostream&os )
{
  // Check arguments
  if( !tk || !PMp || !PMxk || !PMxpk || ( _nf && (!PMf || !PMfp) ) ) return FATAL;

  try{
    // Initialize trajectory integration
    if( !ODEBND_SUNDIALS<T,PMT,PVT>::_INI_PM_STA( _np, PMp ) 
     || !_INI_PM_FSA( _np, PMp ) ) return FATAL;
    _t = tk[0];

    // Bounds on initial states/quadratures
    if( !ODEBND_BASE<T,PMT,PVT>::_IC_PM_SET( options )
     || !ODEBND_BASE<T,PMT,PVT>::_IC_PM_STA( options, _t, NV_DATA_S( _Nx ) )
     || ( _Nq && !ODEBND_BASE<T,PMT,PVT>::_IC_PM_QUAD( options, NV_DATA_S( _Nq ) ) ) )
      { _END_STA(); _END_SEN(); return FATAL; }
    if( PMxk && !PMxk[0] ) PMxk[0] = new PVT[_nx];
    for( unsigned ix=0; PMxk[0] && ix<_nx; ix++ ) PMxk[0][ix] = _PMx[ix];
    if( options.WRAPMIT == Options::ELLIPS && ERxk ) ERxk[0] = _Er;

    // Bounds on initial state/quadrature sensitivities
    if( PMxpk && !PMxpk[0] ) PMxpk[0] = new PVT[_nx*_np];
    for( _isen=0; _isen<_np; _isen++ ){
      if( !_IC_SET_FSA( _isen )
       || !ODEBND_BASE<T,PMT,PVT>::_IC_PM_STA( options, _t, NV_DATA_S(_Ny[_isen]) )
       || ( _Nyq && _Nyq[_isen]
         && !ODEBND_BASE<T,PMT,PVT>::_IC_PM_QUAD( options, NV_DATA_S(_Nyq[_isen]) ) ) ) 
        { _END_STA(); _END_SEN(); return FATAL; }
      for( unsigned iy=0; PMxpk[0] && iy<_nx; iy++ ) PMxpk[0][_isen*_nx+iy] = _PMx[iy];
#ifdef MC__ODEBNDS_SUNDIALS_DEBUG
      std::cout << "qS[" << _isen << "]: " << _Nyq[_isen] << "(" << NV_LENGTH_S( _Nyq[_isen] ) << ")\n";
      for( unsigned i=0; i<NV_LENGTH_S( _Nyq[_isen] ); i++ ) std::cout << NV_Ith_S( _Nyq[_isen], i ) << std::endl;
      { int dum; std::cin >> dum; }
#endif
      if( options.WRAPMIT == Options::ELLIPS && ERxpk && ERxpk[0] ) ERxpk[0][_isen] = _Er;
    }

    // Display & record initial results
    if( options.DISPLAY >= 1 ){
      _print_interm( _t, _nx, PMxk[0], "x", os );
      _print_interm( _nx*_np, PMxpk[0], "xp", os );
    }
    if( options.RESRECORD ){
      results_sta.push_back( Results( _t, _nx, PMxk[0] ) );
      results_sen.push_back( Results( _t, _nx*_np, PMxpk[0] ) );
    }

    // Integrate ODEs through each stage using SUNDIALS
    ODEBND_SUNDIALS<T,PMT,PVT>::pODEBND = pODEBNDS = this;
    if( !ODEBND_SUNDIALS<T,PMT,PVT>::_INI_CVODE
          ( ODEBND_SUNDIALS<T,PMT,PVT>::MC_CVRHSPM__,
            ODEBND_SUNDIALS<T,PMT,PVT>::MC_CVQUADPM__ )
     || !_INI_CVODES( MC_CVFSARHSPM__, MC_CVFSAQUADPM__ ) )
      { _END_STA(); _END_SEN(); return FATAL; }

    for( _istg=0; _istg<ns; _istg++ ){
      // Bounds on state discontinuities (if any) at stage times
      // and integrator reinitialization (if applicable)
      _pos_ic = ( _vIC.size()>=ns? _istg:0 );
      if( _pos_ic
       && ( !ODEBND_BASE<T,PMT,PVT>::_CC_PM_SET( options, _pos_ic )
         || !ODEBND_BASE<T,PMT,PVT>::_CC_PM_STA( options, _t, NV_DATA_S( _Nx ) )
         || !ODEBND_SUNDIALS<T,PMT,PVT>::_CC_CVODE_STA() ) )
        { _END_STA(); _END_SEN(); return FATAL; }
      if( _istg 
       && ( (_Nq && !ODEBND_BASE<T,PMT,PVT>::_IC_PM_QUAD( options, NV_DATA_S( _Nq ) ) ) // quadrature reinitialization
         || !ODEBND_SUNDIALS<T,PMT,PVT>::_CC_CVODE_QUAD() ) )
        { _END_STA(); _END_SEN(); return FATAL; }
      for( _isen=0; _isen<_np; _isen++ ){
        if( _pos_ic
         && ( !_CC_SET_FSA( _pos_ic, _isen )
           || !_CC_PM_SET( options )
           || !_CC_PM_SEN( options, _t, NV_DATA_S( _Nx ), NV_DATA_S(_Ny[_isen]) )
           || ( _isen==_np-1 && !_CC_CVODES_FSA() ) ) )
          { _END_STA(); _END_SEN(); return FATAL; }
        if( _istg
         && ( ( _Nyq && _Nyq[_isen] && !ODEBND_BASE<T,PMT,PVT>::_IC_PM_QUAD( options, NV_DATA_S(_Nyq[_isen]) ) ) //quadrature sensitivity reinitialization
           || ( _isen==_np-1 && !_CC_CVODES_QUAD() ) ) )
            { _END_STA(); _END_SEN(); return FATAL; }
#ifdef MC__ODEBNDS_SUNDIALS_DEBUG
        std::cout << "Entering Stage #" << _istg << std::endl;
        for( unsigned iy=0; iy<NV_LENGTH_S(_Ny[_isen]); iy++ )
          std::cout << "_Ny" << _isen << "[iy] = " << NV_Ith_S(_Ny[_isen],iy) << std::endl;
        for( unsigned iy=0; _nq && iy<NV_LENGTH_S(_Nyq[_isen]); iy++ )
          std::cout << "_Nyq" << _isen << "[iy] = " << NV_Ith_S(_Nyq[_isen],iy) << std::endl;
#endif
      }

      // update list of operations in RHS, QUAD, RHSFSA and QUADFSA
      _pos_rhs  = ( _vRHS.size()<=1?  0: _istg );
      _pos_quad = ( _vQUAD.size()<=1? 0: _istg );
      if( (!_istg || _pos_rhs || _pos_quad)
        && ( !ODEBND_BASE<T,PMT,PVT>::_RHS_PM_SET( options, _pos_rhs, _pos_quad )
          || !_RHS_SET_FSA( _pos_rhs, _pos_quad )
          || !_RHS_PM_SET( options, _np, _nq ) ) )
        { _END_STA(); _END_SEN(); return FATAL; }

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
        stats_sen.numSteps++;
      }

      // Bounds on intermediate states, quadratures and their sensitivities
      if( _nq ){
        for( unsigned ip=0; ip<_np; ip++ ){
          _cv_flag = CVodeGetQuadSens1(_cv_mem, &_t, ip, _Nyq[ip]);
          if( _check_cv_flag( &_cv_flag, "CVodeGetQuadSens", 1) )
            { _END_STA(); _END_SEN(); return FATAL; }
#ifdef MC__ODEBNDS_SUNDIALS_DEBUG
          std::cout << "qS[" << ip << "]: " << _Nyq[ip] << "(" << NV_LENGTH_S( _Nyq[ip] ) << ")\n";
          for( unsigned i=0; i<NV_LENGTH_S( _Nyq[ip] ); i++ ) std::cout << NV_Ith_S( _Nyq[ip], i ) << std::endl;
          { int dum; std::cin >> dum; }
#endif
        }
        _cv_flag = CVodeGetQuad( _cv_mem, &_t, _Nq );
        if( _check_cv_flag(&_cv_flag, "CVodeGetQuad", 1) )
          { _END_STA(); _END_SEN(); return FATAL; }
      }
      _cv_flag = CVodeGetSens(_cv_mem, &_t, _Ny);
      if( _check_cv_flag( &_cv_flag, "CVodeGetSens", 1) )
         { _END_STA(); _END_SEN(); return FATAL; }

      // Add intermediate function terms and derivatives
      _pos_fct = ( _vFCT.size()>=ns? _istg:0 );
      ODEBND_SUNDIALS<T,PMT,PVT>::_GET_PM_STA( options, NV_DATA_S(_Nx), _nq && _Nq? NV_DATA_S(_Nq): 0 );
      if( (_vFCT.size()>=ns || _istg==ns-1)
       && !ODEBND_BASE<T,PMT,PVT>::_FCT_PM_STA( _pos_fct, _t, PMf ) )
        { _END_STA(); _END_SEN(); return FATAL; }
      if( PMxk && !PMxk[_istg+1] ) PMxk[_istg+1] = new PVT[_nx];
      for( unsigned ix=0; ix<_nx; ix++ ) PMxk[_istg+1][ix] = _PMx[ix];
      if( options.WRAPMIT == Options::ELLIPS && ERxk ) ERxk[_istg+1] = _Er;
      for( _isen=0; _isen<_np; _isen++ ){
        _GET_PM_SEN( options, NV_DATA_S(_Nx), NV_DATA_S(_Ny[_isen]), _nq && _Nq? NV_DATA_S(_Nq): 0,
                     _nq, _nq && _Nyq[_isen]? NV_DATA_S(_Nyq[_isen]): 0 );
        if( (_vFCT.size()>=ns || _istg==ns-1)
         && !_FCT_PM_SEN( _pos_fct, _isen, _t, PMfp+_isen*_nf ) )
          { _END_STA(); _END_SEN(); return FATAL; }
        if( PMxpk && !PMxpk[_istg+1] ) PMxpk[_istg+1] = new PVT[_nx*_np];
        for( unsigned iy=0; iy<_nx; iy++ ) PMxpk[_istg+1][_isen*_nx+iy] = _PMy[iy];
        if( options.WRAPMIT == Options::ELLIPS && ERxpk && ERxpk[_istg+1] ) ERxpk[_istg+1][_isen] = _Edy;
      }

      // Display & record stage results
      if( options.DISPLAY >= 1 ){
        _print_interm( _t, _nx, PMxk[_istg+1], "x", os );
        _print_interm( _nx*_np, PMxpk[_istg+1], "xp", os );
      }
      if( options.RESRECORD ){
        results_sta.push_back( Results( tk[_istg+1], _nx, PMxk[_istg+1] ) );
        results_sen.push_back( Results( tk[_istg+1], _nx*_np, PMxpk[_istg+1] ) );
      }
    }

    // Bounds on final quadratures and functions
    if( options.DISPLAY >= 1 ){
      _print_interm( _nf, PMf, "f", os );
      _print_interm( _nf*_np, PMfp, "fp", os );
    }
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

//! @fn template <typename T, typename PMT, typename PVT> inline typename ODEBNDS_SUNDIALS<T,PMT,PVT>::STATUS ODEBNDS_SUNDIALS<T,PMT,PVT>::bounds_FSA(
//! const unsigned ns, const double*tk, const T*Ip, T**Ixk, 
//! T*If, T**Ixpk, T*Ifp, const unsigned nsamp,
//! std::ostream&os=std::cout )
//!
//! This function computes projections of an inner-approximation enclosure of
//! the reachable set of the parametric ODEs using sampling and continuous-time
//! integration, together with forward sensitivity information:
//!  - <a>ns</a>    [input]  number of time stages
//!  - <a>tk</a>    [input]  stage times, including the initial time
//!  - <a>Ip</a>    [input]  interval enclosure of parameter set
//!  - <a>Ixk</a>   [output] interval enclosure of state variables at stage times
//!  - <a>If</a>    [output] interval enclosure of state functions
//!  - <a>Ixpk</a>  [output] interval enclosure of state sensitivity variables at stage times
//!  - <a>Ifp</a>   [output] interval enclosure of state function derivatives
//!  - <a>nsamp</a> [input]  number of samples for each parameter
//!  - <a>os</a>    [input]  output stream [default: std::cout]
//! .
//! The return value is the status.
template <typename T, typename PMT, typename PVT>
inline typename ODEBNDS_SUNDIALS<T,PMT,PVT>::STATUS
ODEBNDS_SUNDIALS<T,PMT,PVT>::bounds_FSA
( const unsigned ns, const double*tk, const T*Ip, T**Ixk,
  T*If, T**Ixpk, T*Ifp, const unsigned nsamp,
  std::ostream&os )
{
  // Sample inner approximation
  STATUS flag = NORMAL;
  pODESLVS.set( *this );
  pODESLVS.options = options.ODESLV;
  if( !_bounds_FSA( ns, tk, Ip, Ixk, If, Ixpk, Ifp, pODESLVS, nsamp, os ) )
    flag = FAILURE;

  // Display results
  if( options.DISPLAY >= 1 ){
    for( unsigned is=0; Ixk && is<=ns; is++ ){
      _print_interm( tk[is], _nx, Ixk[is], "x", os );
      _print_interm( _nx*_np, Ixpk[is], "xp", os );
    }
    if( If )  _print_interm( _nf, If, "f", os );
    if( Ifp ) _print_interm( _nf*_np, Ifp, "fp", os );
  }

  // Record intermediate results
  results_sta.clear();
  results_sen.clear();
  if( options.RESRECORD ){
    for( unsigned is=0; Ixk && is<=ns; is++ )
      results_sta.push_back( Results( tk[is], _nx, Ixk[is] ) );
    for( unsigned is=0; Ixpk && is<=ns; is++ )
      results_sen.push_back( Results( tk[is], _nx*_np, Ixpk[is] ) );
  }
  
  return flag;
}

//! @fn template <typename T, typename PMT, typename PVT> inline typename ODEBNDS_SUNDIALS<T,PMT,PVT>::STATUS ODEBNDS_SUNDIALS<T,PMT,PVT>::bounds_ASA(
//! const unsigned ns, const double*tk, const T*Ip, T**Ixk, 
//! T*If, T**Ilk, T*Ifp, const unsigned nsamp, std::ostream&os=std::cout )
//!
//! This function computes projections of an inner-approximation enclosure of
//! the reachable set of the parametric ODEs using sampling and continuous-time
//! integration, together with forward sensitivity information:
//!  - <a>ns</a>    [input]  number of time stages
//!  - <a>tk</a>    [input]  stage times, including the initial time
//!  - <a>Ip</a>    [input]  interval enclosure of parameter set
//!  - <a>Ixk</a>   [output] interval enclosure of state variables at stage times
//!  - <a>If</a>    [output] interval enclosure of state functions
//!  - <a>Ilk</a>   [output] interval enclosure of adjoint variables at stage times
//!  - <a>Ifp</a>   [output] interval enclosure of state function derivatives
//!  - <a>nsamp</a> [input]  number of samples for each parameter
//!  - <a>os</a>    [input]  output stream [default: std::cout]
//! .
//! The return value is the status.
template <typename T, typename PMT, typename PVT>
inline typename ODEBNDS_SUNDIALS<T,PMT,PVT>::STATUS
ODEBNDS_SUNDIALS<T,PMT,PVT>::bounds_ASA
( const unsigned ns, const double*tk, const T*Ip, T**Ixk,
  T*If, T**Ilk, T*Ifp, const unsigned nsamp,
  std::ostream&os )
{
  // Sample inner approximation
  STATUS flag = NORMAL;
  pODESLVS.set( *this );
  pODESLVS.options = options.ODESLV;
  if( !_bounds_ASA( ns, tk, Ip, Ixk, If, Ilk, Ifp, pODESLVS, nsamp, os ) )
    flag = FAILURE;

  // Display results
  if( options.DISPLAY >= 1 ){
    for( unsigned is=0; Ixk && is<=ns; is++ )
      _print_interm( tk[is], _nx, Ixk[is], "x", os );
    if( If )  _print_interm( _nf, If, "f", os );
    for( unsigned is=ns; Ilk && is<=ns; is-- )
      _print_interm( tk[is], _nx*_nf, Ilk[is], "l", os );
    if( Ifp ) _print_interm( _nf*_np, Ifp, "fp", os );
  }

  // Record intermediate results
  results_sta.clear();
  results_sen.clear();
  if( options.RESRECORD ){
    for( unsigned is=0; Ixk && is<=ns; is++ )
      results_sta.push_back( Results( tk[is], _nx, Ixk[is] ) );
    for( unsigned is=0; Ilk && is<=ns; is++ )
      results_sen.push_back( Results( tk[is], _nx*_nf, Ilk[is] ) );
  }

  return flag;
}

} // end namescape mc

#endif




















