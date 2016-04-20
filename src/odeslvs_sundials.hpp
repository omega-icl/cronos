// Copyright (C) 2015- Benoit Chachuat & Nikola Peric, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

#ifndef MC__ODESLVS_SUNDIALS_HPP
#define MC__ODESLVS_SUNDIALS_HPP

#undef  MC__ODESLVS_SUNDIALS_DEBUG

#include "odeslvs_base.hpp"
#include "odeslv_sundials.hpp"

#define MC__ODESLVS_SUNDIALS_USE_BAD

namespace mc
{
//! @brief C++ class computing solutions of parametric ODEs with forward/adjoint sensitivity analysis capability using SUNDIALS and MC++.
////////////////////////////////////////////////////////////////////////
//! mc::ODESLV_SUNDIALS is a C++ class for solution of IVPs in ODEs
//! with forward/adjoint sensitivity analysis capability using the code
//! CVODES in SUNDIALS and MC++.
////////////////////////////////////////////////////////////////////////
class ODESLVS_SUNDIALS:
  public virtual BASE_DE,
  //public virtual BASE_SUNDIALS,
  public virtual ODESLVS_BASE, 
  public virtual ODESLV_SUNDIALS
{
 protected:
  typedef int (*CVRhsFn)( realtype t, N_Vector y, N_Vector ydot, void *user_data );
  typedef int (*CVADJRhsFn)( realtype t, N_Vector x, N_Vector y, N_Vector ydot, void *user_data );
  typedef int (*CVSENRhs1Fn)( int Ns, realtype t, N_Vector x, N_Vector xdot, int is,
    N_Vector y, N_Vector ydot, void *user_data, N_Vector tmp1, N_Vector tmp2 );
  typedef int (*CVSENQuadFn)( int Ns, realtype t, N_Vector x, N_Vector *y, N_Vector qdot,
    N_Vector *qSdot, void *user_data, N_Vector tmp1, N_Vector tmp2 );
  typedef int (*CVDlsDenseJacFn)( long int N, realtype t, N_Vector y, N_Vector fy,
     DlsMat Jac, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3 );
  typedef int (*CVDlsDenseJacFnB)( long int NeqB, realtype t, N_Vector y, N_Vector yB,
     N_Vector fyB, DlsMat JacB, void *user_dataB, N_Vector tmp1B, N_Vector tmp2B,
     N_Vector tmp3B );

 protected:
  //! @brief current function index
  unsigned _ifct;

  //! @brief current parameter sensitivity index
  unsigned _isen;

  //! @brief size of N_vector arrays
  unsigned _nvec;

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
  static ODESLVS_SUNDIALS *pODESLVS;

 public:
  /** @ingroup ODESLV
   *  @{
   */
  //! @brief Default constructor
  ODESLVS_SUNDIALS
    ();

  //! @brief Virtual destructor
  virtual ~ODESLVS_SUNDIALS
    ();

  //! @brief Statistics for sensitivity/adjoint integration
  Stats stats_sen;

  //! @brief Vector storing interval adjoint bounds (see Options::RESRECORD)
  std::vector< Results > results_sen;

 //! @brief Propagate states and state-sensitivities forward in time through every time stages
  STATUS states_FSA
    ( const unsigned ns, const double*tk, const double*p, double**xk,
      double*f, double**xpk, double*fp, std::ostream&os=std::cout );

  //! @brief Propagate states and adjoints forward and backward in time through every time stages
  STATUS states_ASA
    ( const unsigned ns, const double*tk, const double*p, double**xk,
      double*f, double**lk, double*fp, std::ostream&os=std::cout );

  //! @brief Record state and sensitivity bounds in files <a>obndsta</a> and <a>obndsa</a>, with accuracy of <a>iprec</a> digits
  void record
    ( std::ofstream&obndsta, std::ofstream&obndsa, const unsigned iprec=5 ) const
    { this->ODESLV_SUNDIALS::record( obndsta, iprec );
      this->ODESLV_BASE::_record( obndsa, results_sen, iprec ); }

  //! @brief Record state and sensitivity bounds in files <a>obndsta</a> and <a>obndsa</a>, with accuracy of <a>iprec</a> digits
  void record
    ( std::ofstream&obndsta, const unsigned iprec=5 ) const
    { ODESLV_SUNDIALS::record( obndsta, iprec ); }
  /** @} */

 private:
  //! @brief Function to initialize CVodes memory block (virtual)
  virtual bool _INI_CVODE
    ( CVRhsFn MC_CVRHS, CVRhsFn MC_CVQUAD, CVDlsDenseJacFn MC_CVJAC );

  //! @brief Function to initialize CVodeS memory block for forward sensitivity
  bool _INI_CVODES
    ( CVSENRhs1Fn MC_CVSENRHS, CVSENQuadFn MC_CVSENQUAD );

  //! @brief Function to initialize CVodeS memory block for adjoint sensitivity
  bool _INI_CVODES
    ( CVADJRhsFn MC_CVADJRHS, CVADJRhsFn MC_CVADJQUAD, CVDlsDenseJacFnB MC_CVADJJAC,
      const unsigned ifct, int&indexB, unsigned&iusrB );

  //! @brief Function to reinitialize CVodeS memory block for forward sensitivity
  bool _CC_CVODES_FSA
    ();

  //! @brief Function to reinitialize CVodeS memory block for forward quadrature sensitivity
  bool _CC_CVODES_QUAD
    ();

  //! @brief Function to reinitialize CVodeS memory block for adjoint sensitivity
  bool _CC_CVODES_ASA
    ( const unsigned ifct, const int indexB );

  //! @brief Function to finalize sensitivity/adjoint bounding
  void _END_SEN
    ();

  //! @brief Function to initialize adjoint sensitivity analysis
  bool _INI_ASA
    ( const double *p );

  //! @brief Static wrapper to function computing the adjoint RHS
  static int MC_CVASARHSD__
    ( realtype t, N_Vector x, N_Vector y, N_Vector ydot, void *user_data );

  //! @brief Static wrapper to function computing the adjoint quadrature RHS
  static int MC_CVASAQUADD__
    ( realtype t, N_Vector x, N_Vector y, N_Vector qdot, void *user_data );

  //! @brief Static wrapper to function to calculate the adjoint RHS Jacobian
  static int MC_CVASAJACD__
    ( long int NeqB, realtype t, N_Vector y, N_Vector yB, N_Vector fyB, DlsMat JacB,
      void *user_dataB, N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B );

  //! @brief Function to initialize forward sensitivity analysis
  bool _INI_FSA
    ( const double *p );

  //! @brief Static wrapper to function to calculate the sensitivity RHS
  static int MC_CVFSARHSD__
    ( int Ns, realtype t, N_Vector x, N_Vector xdot, int is, N_Vector y,
      N_Vector ydot, void *user_data, N_Vector tmp1, N_Vector tmp2 );

  //! @brief Static wrapper to function computing the sensitivity quadrature RHS
  static int MC_CVFSAQUADD__
    ( int Ns, realtype t, N_Vector x, N_Vector *y, N_Vector qdot, N_Vector *qSdot, 
      void *user_data, N_Vector tmp1, N_Vector tmp2 );

  //! @brief Static wrapper to function to calculate the sensitivity RHS Jacobian
  static int MC_CVFSAJACD__
    ( long int NeqB, realtype t, N_Vector y, N_Vector yB, N_Vector fyB, DlsMat JacB,
      void *user_dataB, N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B );

  //! @brief Private methods to block default compiler methods
  ODESLVS_SUNDIALS(const ODESLVS_SUNDIALS&);
  ODESLVS_SUNDIALS& operator=(const ODESLVS_SUNDIALS&);
};

ODESLVS_SUNDIALS* ODESLVS_SUNDIALS::pODESLVS = 0;

inline
ODESLVS_SUNDIALS::ODESLVS_SUNDIALS
()
//: BASE_SUNDIALS(), 
: ODESLVS_BASE(), ODESLV_SUNDIALS(),
  _ifct(0), _isen(0), _nvec(0), _Ny(0), _Nyq(0), _indexB(0), _iusrB(0)
{}

inline
ODESLVS_SUNDIALS::~ODESLVS_SUNDIALS
()
{
  if( _Ny )   N_VDestroyVectorArray_Serial( _Ny,  _nvec );
  if( _Nyq )  N_VDestroyVectorArray_Serial( _Nyq, _nvec );
  delete[] _indexB;
  delete[] _iusrB;
}

inline bool
ODESLVS_SUNDIALS::_INI_CVODE
( CVRhsFn MC_CVRHS, CVRhsFn MC_CVQUAD, CVDlsDenseJacFn MC_CVJAC )
{
  // Call _INI_CVODE in ODEBND_SUNDIALS
  this->ODESLV_SUNDIALS::_INI_CVODE( MC_CVRHS, MC_CVQUAD, MC_CVJAC );

  // Allocate memory for adjoint integration
  _cv_flag = CVodeAdjInit( _cv_mem, options.ASACHKPT, options.ASAINTERP );
  if( _check_cv_flag(&_cv_flag, "CVodeAdjInit", 1) ) return false;

  // Reinitialize adjoint holding vectors
  delete[] _indexB; _indexB = new int[_nf];
  delete[] _iusrB;  _iusrB  = new unsigned[_nf];

  return true;
}

inline bool
ODESLVS_SUNDIALS::_INI_CVODES
( CVADJRhsFn MC_CVADJRHS, CVADJRhsFn MC_CVADJQUAD, CVDlsDenseJacFnB MC_CVADJJAC,
  const unsigned ifct, int&indexB, unsigned&iusrB )
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
  // the RHS function
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
       _cv_flag = CVDlsSetDenseJacFnB( _cv_mem, indexB, MC_CVADJJAC );
       if ( _check_cv_flag(&_cv_flag, "CVDlsSetDenseJacFnB", 1) ) return false;
       break;
    }
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

inline bool
ODESLVS_SUNDIALS::_INI_CVODES
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

inline bool
ODESLVS_SUNDIALS::_CC_CVODES_ASA
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

inline bool
ODESLVS_SUNDIALS::_CC_CVODES_FSA
()
{
  // Reinitialize CVodeS memory block for current sensitivity _Ny
  _cv_flag = CVodeSensReInit( _cv_mem, options.FSACORR, _Ny );
  if( _check_cv_flag(&_cv_flag, "CVodeSensReInit", 1) ) return false;

  return true;
}

inline bool
ODESLVS_SUNDIALS::_CC_CVODES_QUAD
()
{
  // Reinitialize CVode memory block for current sensitivity quarature _Nyq
  if( !_nq ) return true;
  _cv_flag = CVodeQuadSensReInit( _cv_mem, _Nyq );
  if( _check_cv_flag(&_cv_flag, "CVodeQuadSensReInit", 1) ) return false;

  return true;
}

inline void
ODESLVS_SUNDIALS::_END_SEN()
{
  // Get final CPU time
  _final_stats( stats_sen );
}

inline bool
ODESLVS_SUNDIALS::_INI_ASA
( const double*p )
{
  // Initialize bound propagation
  if( !_INI_D_SEN( p, _nf, _np ) )
    return false;

  // Set SUNDIALS adjoint/quadrature arrays
  if( _nvec != _nf ){
    if( _Ny )   N_VDestroyVectorArray_Serial( _Ny,  _nvec );
    if( _Nyq )  N_VDestroyVectorArray_Serial( _Nyq, _nvec );
    _nvec = _nf;
    _Ny  = N_VCloneVectorArray_Serial( _nvec, _Nx );
    _Nyq = N_VCloneVectorArray_Serial( _nvec, _Nx );
  }
  for( unsigned i=0; i<_nf; i++ ){
    if( !_Ny[i] || NV_LENGTH_S( _Ny[i] ) != _ny ){
      if( _Ny[i] ) N_VDestroy_Serial( _Ny[i] );
      _Ny[i] = N_VNew_Serial( _ny );
    }
    if( !_Nyq[i] || NV_LENGTH_S( _Nyq[i] ) != _np ){
      if( _Nyq[i] ) N_VDestroy_Serial( _Nyq[i] );
      _Nyq[i] = N_VNew_Serial( _np );
    }
  }

  // Reset result record and statistics
  results_sen.clear();
  _init_stats( stats_sen );

  return true;
}

inline int
ODESLVS_SUNDIALS::MC_CVASARHSD__
( realtype t, N_Vector x, N_Vector y, N_Vector ydot, void *user_data )
{
  ODESLVS_SUNDIALS *pODESLVS = ODESLVS_SUNDIALS::pODESLVS;
  pODESLVS->_ifct = *static_cast<unsigned*>( user_data );
#ifdef MC__ODESLVS_SUNDIALS_DEBUG
  std::cout << "@t=" << t << "\nx:\n";
  for( unsigned i=0; i<NV_LENGTH_S( x ); i++ ) std::cout << NV_Ith_S( x, i ) << std::endl;
  std::cout << "y" << pODESLVS->_ifct << ":\n";
  for( unsigned i=0; i<NV_LENGTH_S( y ); i++ ) std::cout << NV_Ith_S( y, i ) << std::endl;
#endif
  bool flag = pODESLVS->_RHS_D_SEN( t, NV_DATA_S( x ), NV_DATA_S( y ),
    NV_DATA_S( ydot ), pODESLVS->_ifct );
#ifdef MC__ODESLVS_SUNDIALS_DEBUG
  std::cout << "ydot" << pODESLVS->_ifct << ":\n";
  for( unsigned i=0; i<NV_LENGTH_S( ydot ); i++ ) std::cout << NV_Ith_S( ydot, i ) << std::endl;
  { int dum; std::cin >> dum; }
#endif
  ODESLVS_SUNDIALS::pODESLVS = pODESLVS;
  pODESLVS->stats_sen.numRHS++;
  return( flag? 0: -1 );
}

inline int
ODESLVS_SUNDIALS::MC_CVASAJACD__
( long int NeqB, realtype t, N_Vector x, N_Vector y, N_Vector fy, DlsMat JacB,
  void *user_dataB, N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B )
{
  ODESLVS_SUNDIALS *pODESLVS = ODESLVS_SUNDIALS::pODESLVS;
  bool flag = pODESLVS->_JAC_D_SEN( t, NV_DATA_S( x ), NV_DATA_S( y ), JacB->cols );
  ODESLVS_SUNDIALS::pODESLVS = pODESLVS;
  pODESLVS->stats_sen.numJAC++; // increment JAC counter
  return( flag? 0: -1 );
}

inline int
ODESLVS_SUNDIALS::MC_CVASAQUADD__
( realtype t, N_Vector x, N_Vector y, N_Vector qdot, void *user_data )
{
  ODESLVS_SUNDIALS *pODESLVS = ODESLVS_SUNDIALS::pODESLVS;
  pODESLVS->_ifct = *static_cast<unsigned*>( user_data );
  bool flag = pODESLVS->_RHS_D_QUAD( pODESLVS->_np, NV_DATA_S( qdot ),
    pODESLVS->_ifct );
  ODESLVS_SUNDIALS::pODESLVS = pODESLVS;
  return( flag? 0: -1 );
}

//! @fn inline typename ODESLVS_SUNDIALS::STATUS ODESLVS_SUNDIALS::states_ASA
//!( const unsigned ns, const double*tk, const double*p, double**xk,
//!  double*f, double**lk, double*fp, std::ostream&os=std::cout )
//!
//! This function computes a solution to the parametric ODEs with adjoint
//! sensitivity analysis:
//!  - <a>ns</a> [input]  number of time stages
//!  - <a>tk</a> [input]  stage times, including the initial time
//!  - <a>p</a>  [input]  parameter values
//!  - <a>xk</a> [output] state values at stage times
//!  - <a>f</a>  [output] function values
//!  - <a>lk</a> [output] adjoint values at stage times
//!  - <a>fp</a> [output] function derivatives
//!  - <a>os</a> [input]  output stream [default: std::cout]
//! .
//! The return value is the status.
inline typename ODESLVS_SUNDIALS::STATUS
ODESLVS_SUNDIALS::states_ASA
( const unsigned ns, const double*tk, const double*p, double**xk,
  double*f, double**lk, double*fp, std::ostream&os )
{
  // Compute state bounds and store intermediate results
  STATUS flag = NORMAL;
  flag = ODESLV_SUNDIALS::_states( ns, tk, p, xk, f, true, os);
  if( flag != NORMAL ) return flag;

  // Nothing to do if no functions are defined
  if( !_nf ) return NORMAL;

  // Check size
  if( !lk || !fp ) return FATAL;

  try{
    // Initialize adjoint integration
    if( !_INI_ASA( p )) return FATAL;
    _t = tk[ns];

    // Bounds on terminal adjoints/quadratures
    if( lk && !lk[ns] ) lk[ns] = new double[_nx*_nf];
    for( _ifct=0; _ifct < _nf; _ifct++ ){
      _pos_fct = ( _vFCT.size()>=ns? ns-1:0 );
      if( !_TC_SET_ASA( _pos_fct, _ifct )
       || !_TC_D_SEN( _t, _vec_sta[ns].data(), NV_DATA_S(_Ny[_ifct]) )
       || ( _Nyq && _Nyq[_ifct] && !_TC_D_QUAD_ASA( NV_DATA_S(_Nyq[_ifct]) ) ) )
        { _END_SEN(); return FATAL; }
      _GET_D_SEN( NV_DATA_S(_Ny[_ifct]), _np, _Nyq? NV_DATA_S(_Nyq[_ifct]): 0 );
      for( unsigned iy=0; lk[ns] && iy<_ny; iy++ )
        lk[ns][_ifct*_nx+iy] = _Dy[iy];
      for( unsigned iq=0; iq<_np; iq++ )
        fp[_ifct*_np+iq] = _Dyq[iq];
    }

    // Display & record adjoint terminal results
    if( options.DISPLAY >= 1 )
      _print_interm( tk[ns], _nf*_nx, lk[ns], "l", os );
    if( options.RESRECORD )
      results_sen.push_back( Results( tk[ns], _nx*_nf, lk[ns] ) );//_nf

    // Initialization of adjoint integration
    for( _ifct=0; _ifct < _nf; _ifct++ )
      if( !_INI_CVODES( MC_CVASARHSD__, MC_CVASAQUADD__, MC_CVASAJACD__,
        _ifct, _indexB[_ifct], _iusrB[_ifct] ) )
        { _END_SEN(); return FATAL;}

    // Integrate adjoint ODEs through each stage using SUNDIALS
    pODESLVS = this;
    for( _istg=ns; _istg>0; _istg-- ){

      // Update list of operations in RHSADJ and QUADADJ
      _pos_rhs  = ( _vRHS.size() <=1? 0: _istg-1 );
      _pos_quad = ( _vQUAD.size()<=1? 0: _istg-1 );
      if( !_RHS_SET_ASA( _pos_rhs, _pos_quad, _pos_fct )
       || !_RHS_D_SET( _nf, _np ) )
        { _END_SEN(); return FATAL; }

      // Propagate bounds backward to previous stage time
      _cv_flag = CVodeB( _cv_mem, tk[_istg-1], CV_NORMAL );
      if( _check_cv_flag(&_cv_flag, "CVodeB", 1) )
        { _END_SEN(); return FATAL; }
      _t = tk[_istg-1];

      // Bounds on states/adjoints/quadratures at stage time
      //stats_sen.numSteps = 0; 
#ifdef MC__ODESLVS_SUNDIALS_DEBUG
      for( unsigned ix=0; ix<_nx; ix++ )
        std::cout << "_vec_sta = " << _vec_sta[_istg-1][ix] << std::endl;
#endif
      for( _ifct=0; _ifct < _nf; _ifct++ ){
        _cv_flag = CVodeGetB( _cv_mem, _indexB[_ifct], &_t, _Ny[_ifct]);
        if( _check_cv_flag( &_cv_flag, "CVodeGetB", 1) )
          { _END_SEN(); return FATAL; }
        _cv_flag = CVodeGetQuadB( _cv_mem, _indexB[_ifct], &_t, _Nyq[_ifct]);
#ifdef MC__ODESLVS_SUNDIALS_DEBUG
          for( unsigned iy=0; iy<_ny; iy++ )
            std::cout << "_Ny" << _ifct << "[" << iy << "] = " << NV_Ith_S(_Ny[_ifct],iy) << std::endl;
#endif
        if( _check_cv_flag( &_cv_flag, "CVodeGetQuadB", 1) )
          { _END_SEN(); return FATAL; }
#ifdef MC__ODESLVS_SUNDIALS_DEBUG
          for( unsigned ip=0; ip<_np; ip++ )
            std::cout << "_Nyq" << _ifct << "[" << ip << "] = " << NV_Ith_S(_Nyq[_ifct],ip) << std::endl;
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
             || !_CC_D_SEN( _t, _vec_sta[_istg-1].data(), NV_DATA_S(_Ny[_ifct]) )
             || ( _Nyq && _Nyq[_ifct] && !_CC_D_QUAD_ASA( NV_DATA_S(_Nyq[_ifct]) ) ) ) )
            { _END_SEN(); return FATAL; }
          //else if( !_pos_fct )
             _GET_D_SEN( NV_DATA_S(_Ny[_ifct]), _np, _Nyq && _Nyq[_ifct]? NV_DATA_S(_Nyq[_ifct]): 0 );
#ifdef MC__ODESLVS_SUNDIALS_DEBUG
          for( unsigned iy=0; iy<_ny; iy++ )
            std::cout << "_Ny" << _ifct << "[" << iy << "] = " << NV_Ith_S(_Ny[_ifct],iy) << std::endl;
          for( unsigned ip=0; ip<_np; ip++ )
            std::cout << "_Nyq" << _ifct << "[" << ip << "] = " << NV_Ith_S(_Nyq[_ifct],ip) << std::endl;
#endif

          // Reset ODE solver - needed in case of discontinuity
          if( !_CC_CVODES_ASA( _ifct, _indexB[_ifct] ) )
            { _END_SEN(); return FATAL; }
        }
        // Add initial state contribution to derivative bounds
        else if( !_IC_SET_ASA()
              || !_IC_D_SEN( _t, _vec_sta[_istg-1].data(), NV_DATA_S(_Ny[_ifct]) )
              || ( _Nyq && _Nyq[_ifct] && !_IC_D_QUAD_ASA( NV_DATA_S(_Nyq[_ifct]) ) ) )
          { _END_SEN(); return FATAL; }

        // Keep track of results at stage times
        if( lk && !lk[_istg-1] ) lk[_istg-1] = new double[_nx*_nf];
        for( unsigned iy=0; lk[_istg-1] && iy<_ny; iy++ )
          lk[_istg-1][_ifct*_nx+iy] = _Dy[iy];
        for( unsigned iq=0; iq<_np; iq++ )
          fp[_ifct*_np+iq] = _Dyq[iq];
      }

      // Display & record adjoint intermediate results
      if( options.DISPLAY >= 1 )
        _print_interm( tk[_istg-1], _nf*_nx, lk[_istg-1], "l", os );
      if( options.RESRECORD )
        results_sen.push_back( Results( tk[_istg-1], _nf*_nx, lk[_istg-1] ) );
    }
    if( options.DISPLAY >= 1 )
      _print_interm( _nf*_np, fp, "fp", os );
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

inline bool
ODESLVS_SUNDIALS::_INI_FSA
( const double *p )
{
  // Initialize bound propagation
  if( !_INI_D_SEN( p, _np, _nq ) )
    return false;

  // Set SUNDIALS sensitivity/quadrature arrays
  if( _nvec != _np ){
    if( _Ny )   N_VDestroyVectorArray_Serial( _Ny,  _nvec ); _Ny = 0;
    if( _Nyq )  N_VDestroyVectorArray_Serial( _Nyq, _nvec ); _Nyq = 0;
    _nvec = _np;
    _Ny  = N_VCloneVectorArray_Serial( _np, _Nx );
    if( _nq ) _Nyq = N_VCloneVectorArray_Serial( _np, _Nq );
  }
  for( unsigned i=0; i<_np; i++ ){
    if( !_Ny[i] || NV_LENGTH_S( _Ny[i] ) != _ny ){
      if( _Ny[i] ) N_VDestroy_Serial( _Ny[i] );
      _Ny[i] = N_VNew_Serial( _ny );
    }
    if( _Nyq && (!_Nyq[i] || NV_LENGTH_S( _Nyq[i]) ) != _nq ){
      if( _Nyq[i] ) N_VDestroy_Serial( _Nyq[i] );
      _Nyq[i] = _nq? N_VNew_Serial( _nq ): 0;
    }
  }

  // Reset result record and statistics
  results_sen.clear();
  _init_stats( stats_sen );

  return true;
}

inline int
ODESLVS_SUNDIALS::MC_CVFSARHSD__
( int Ns, realtype t, N_Vector x, N_Vector xdot, int is, N_Vector y,
  N_Vector ydot, void *user_data, N_Vector tmp1, N_Vector tmp2 )
{
  ODESLVS_SUNDIALS *pODESLVS = ODESLVS_SUNDIALS::pODESLVS;
#ifdef MC__ODESLVS_SUNDIALS_DEBUG
  std::cout << "@t=" << t << "\nx:\n";
  for( unsigned i=0; i<NV_LENGTH_S( x ); i++ ) std::cout << NV_Ith_S( x, i ) << std::endl;
  std::cout << "y:\n";
  for( unsigned i=0; i<NV_LENGTH_S( y ); i++ ) std::cout << NV_Ith_S( y, i ) << std::endl;
#endif
  bool flag = pODESLVS->_RHS_D_SEN( t, NV_DATA_S( x ), NV_DATA_S( y ),
    NV_DATA_S( ydot ), is );
#ifdef MC__ODESLVS_SUNDIALS_DEBUG
  std::cout << "ydot:\n";
  for( unsigned i=0; i<NV_LENGTH_S( ydot ); i++ ) std::cout << NV_Ith_S( ydot, i ) << std::endl;
  { int dum; std::cin >> dum; }
#endif
  ODESLVS_SUNDIALS::pODESLVS = pODESLVS;
  pODESLVS->stats_sen.numRHS++;
  return( flag? 0: -1 );
}

inline int
ODESLVS_SUNDIALS::MC_CVFSAQUADD__
( int Ns, realtype t, N_Vector x, N_Vector *y, N_Vector qdot, N_Vector *qSdot, 
  void *user_data, N_Vector tmp1, N_Vector tmp2 )
{
  ODESLVS_SUNDIALS *pODESLVS = ODESLVS_SUNDIALS::pODESLVS;
#ifdef MC__ODESLVS_SUNDIALS_DEBUG
  std::cout << "@t=" << t << "\nx:\n";
  for( unsigned i=0; i<NV_LENGTH_S( x ); i++ ) std::cout << NV_Ith_S( x, i ) << std::endl;
  std::cout << "qdot:\n";
  for( unsigned i=0; i<NV_LENGTH_S( qdot ); i++ ) std::cout << NV_Ith_S( qdot, i ) << std::endl;
#endif
  bool flag = true;
  for( int is=0; is<Ns && flag; is++ ){
    pODESLVS->_GET_D_SEN( NV_DATA_S(x), NV_DATA_S(y[is]), (realtype*)0, 0, (realtype*)0 );
#ifdef MC__ODESLVS_SUNDIALS_DEBUG
    std::cout << "y:\n";
    for( unsigned i=0; i<NV_LENGTH_S( y[is] ); i++ ) std::cout << NV_Ith_S( y[is], i ) << std::endl;
#endif
    flag = pODESLVS->_RHS_D_QUAD( pODESLVS->_nq, NV_DATA_S( qSdot[is] ), is );
#ifdef MC__ODESLVS_SUNDIALS_DEBUG
    std::cout << "qSdot:\n";
    for( unsigned i=0; i<NV_LENGTH_S( qSdot[is] ); i++ ) std::cout << NV_Ith_S( qSdot[is], i ) << std::endl;
    { int dum; std::cin >> dum; }
#endif
  }
  ODESLVS_SUNDIALS::pODESLVS = pODESLVS;
  return( flag? 0: -1 );
}

//! @fn inline typename ODESLVS_SUNDIALS::STATUS ODESLVS_SUNDIALS::states_FSA(
//! const unsigned ns, const double*tk, const double*p, double**xk, 
//! double*f, double**xpk, double*fp, std::ostream&os=std::cout )
//!
//! This function computes a solution to the parametric ODEs with forward
//! sensitivity analysis:
//!  - <a>ns</a>  [input]  number of time stages
//!  - <a>tk</a>  [input]  stage times, including the initial time
//!  - <a>p</a>   [input]  parameter values
//!  - <a>xk</a>  [output] state values at stage times
//!  - <a>f</a>   [output] function values
//!  - <a>xpk</a> [output] state-sensitivity values at stage times
//!  - <a>fp</a>  [output] function derivatives
//!  - <a>os</a>  [input]  output stream [default: std::cout]
//! .
//! The return value is the status.
inline typename ODESLVS_SUNDIALS::STATUS
ODESLVS_SUNDIALS::states_FSA
( const unsigned ns, const double*tk, const double*p, double**xk, 
  double*f, double**xpk, double*fp, std::ostream&os )
{
  // Check arguments
  if( !tk || !p || !xk || !xpk || ( _nf && (!f || !fp) ) ) return FATAL;

  try{
    // Initialize trajectory integration
    if( !ODESLV_SUNDIALS::_INI_STA( p ) 
     || !_INI_FSA( p ) ) return FATAL;
    _t = tk[0];

    // Bounds on initial states/quadratures
    if( !ODESLV_BASE::_IC_D_SET()
     || !ODESLV_BASE::_IC_D_STA( _t, NV_DATA_S( _Nx ) )
     || ( _Nq && !ODESLV_BASE::_IC_D_QUAD( NV_DATA_S( _Nq ) ) ) )
      { _END_STA(); _END_SEN(); return FATAL; }
    _GET_D_STA( NV_DATA_S(_Nx), _nq && _Nq? NV_DATA_S(_Nq): 0 );
    if( xk && !xk[0] ) xk[0] = new double[_nx];
    for( unsigned ix=0; xk[0] && ix<_nx; ix++ ) xk[0][ix] = ODESLV_BASE::_Dx[ix];

    // Bounds on initial state/quadrature sensitivities
    if( xpk && !xpk[0] ) xpk[0] = new double[_nx*_np];
    for( _isen=0; _isen<_np; _isen++ ){
      if( !_IC_SET_FSA( _isen )
       || !ODESLV_BASE::_IC_D_STA( _t, NV_DATA_S(_Ny[_isen]) )
       || ( _Nyq && _Nyq[_isen] && !ODESLV_BASE::_IC_D_QUAD( NV_DATA_S(_Nyq[_isen]) ) ) ) 
        { _END_STA(); _END_SEN(); return FATAL; }
      _GET_D_SEN( NV_DATA_S(_Ny[_isen]), _np, _nq && _Nyq? NV_DATA_S(_Nyq[_isen]): 0 );
      for( unsigned iy=0; xpk[0] && iy<_ny; iy++ ) xpk[0][_isen*_ny+iy] = _Dy[iy];
#ifdef MC__ODESLVS_SUNDIALS_DEBUG
      std::cout << "qS[" << _isen << "]: " << _Nyq[_isen] << "(" << NV_LENGTH_S( _Nyq[_isen] ) << ")\n";
      for( unsigned i=0; i<NV_LENGTH_S( _Nyq[_isen] ); i++ ) std::cout << NV_Ith_S( _Nyq[_isen], i ) << std::endl;
      { int dum; std::cin >> dum; }
#endif
    }

    // Display & record initial results
    if( options.DISPLAY >= 1 ){
      _print_interm( _t, _nx, xk[0], "x", os );
      _print_interm( _nx*_np, xpk[0], "xp", os );
    }
    if( options.RESRECORD ){
      results_sta.push_back( Results( _t, _nx, xk[0] ) );
      results_sen.push_back( Results( _t, _nx*_np, xpk[0] ) );
    }

    // Integrate ODEs through each stage using SUNDIALS
    ODESLV_SUNDIALS::pODESLV = pODESLVS = this;
    if( !ODESLV_SUNDIALS::_INI_CVODE( MC_CVRHSD__, MC_CVQUADD__, MC_CVJACD__ )
     || !_INI_CVODES( MC_CVFSARHSD__, MC_CVFSAQUADD__ ) )
      { _END_STA(); _END_SEN(); return FATAL; }

    for( _istg=0; _istg<ns; _istg++ ){
      // Bounds on state discontinuities (if any) at stage times
      // and integrator reinitialization (if applicable)
      _pos_ic = ( _vIC.size()>=ns? _istg:0 );
      if( _pos_ic
       && ( !ODESLV_BASE::_CC_D_SET( _pos_ic )
         || !ODESLV_BASE::_CC_D_STA( _t, NV_DATA_S( _Nx ) )
         || !ODESLV_SUNDIALS::_CC_CVODE_STA() ) )
        { _END_STA(); _END_SEN(); return FAILURE; }
      if( _istg 
       && ( ( _Nq && !ODESLV_BASE::_IC_D_QUAD( NV_DATA_S( _Nq ) ) ) // quadrature reinitialization
         || !ODESLV_SUNDIALS::_CC_CVODE_QUAD() ) )
        { _END_STA(); _END_SEN(); return FAILURE; }
      for( _isen=0; _isen<_np; _isen++ ){
        if( _pos_ic
         && ( !_CC_SET_FSA( _pos_ic, _isen )
           || !_CC_D_SEN( _t, NV_DATA_S( _Nx ), NV_DATA_S(_Ny[_isen]) )
           || ( !_isen && !_CC_CVODES_FSA() ) ) )
            { _END_STA(); _END_SEN(); return FATAL; }
        if( _istg
         && ( ( _Nyq && _Nyq[_isen] && !ODESLV_BASE::_IC_D_QUAD( NV_DATA_S(_Nyq[_isen]) ) ) //quadrature sensitivity reinitialization
           || ( !_isen && !_CC_CVODES_QUAD() ) ) )
            { _END_STA(); _END_SEN(); return FATAL; }
#ifdef MC__ODESLVS_SUNDIALS_DEBUG
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
        && ( !ODESLV_BASE::_RHS_D_SET( _pos_rhs, _pos_quad )
          || !_RHS_SET_FSA( _pos_rhs, _pos_quad )
          || !_RHS_D_SET( _np, _nq ) ) )
        { _END_STA(); _END_SEN(); return FATAL; }

      // integrate till end of time stage
      _cv_flag = CVodeSetStopTime( _cv_mem, tk[_istg+1] );
      if( _check_cv_flag(&_cv_flag, "CVodeSetStopTime", 1) )
        { _END_STA(); return FATAL; }
      while( _t < tk[_istg+1] ){
        _cv_flag = CVode( _cv_mem, tk[_istg+1], _Nx, &_t, CV_ONE_STEP );
        if( _check_cv_flag(&_cv_flag, "CVode", 1)
         || (options.NMAX && stats_sta.numSteps > options.NMAX) )
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
#ifdef MC__ODESLVS_SUNDIALS_DEBUG
          std::cout << "qS[" << ip << "]: " << _Nyq[ip] << "(" << NV_LENGTH_S( _Nyq[ip] ) << ")\n";
          for( unsigned i=0; i<NV_LENGTH_S( _Nyq[ip] ); i++ ) std::cout << NV_Ith_S( _Nyq[ip], i ) << std::endl;
          { int dum; std::cin >> dum; }
#endif
        }
        _cv_flag = CVodeGetQuad( _cv_mem, &_t, _Nq );
        if( _check_cv_flag(&_cv_flag, "CVodeGetQuad", 1) )
          { _END_STA(); _END_SEN(); return FATAL; }
      }
      _cv_flag = CVodeGetSens(_cv_mem, &_t, _Ny );
      if( _check_cv_flag( &_cv_flag, "CVodeGetSens", 1) )
         { _END_STA(); _END_SEN(); return FATAL; }

      // Add intermediate function terms and derivatives
      _pos_fct = ( _vFCT.size()>=ns? _istg:0 );
      ODESLV_SUNDIALS::_GET_D_STA( NV_DATA_S(_Nx), _nq && _Nq? NV_DATA_S(_Nq): 0 );
      if( (_vFCT.size()>=ns || _istg==ns-1)
       && !ODESLV_BASE::_FCT_D_STA( _pos_fct, _t, f ) )
        { _END_STA(); _END_SEN(); return FATAL; }
      if( xk && !xk[_istg+1] ) xk[_istg+1] = new double[_nx];
      for( unsigned ix=0; ix<_nx; ix++ ) xk[_istg+1][ix] = ODESLV_BASE::_Dx[ix];
      for( _isen=0; _isen<_np; _isen++ ){
        _GET_D_SEN( NV_DATA_S(_Nx), NV_DATA_S(_Ny[_isen]), _nq && _Nq? NV_DATA_S(_Nq): 0,
                    _nq, _nq && _Nyq[_isen]? NV_DATA_S(_Nyq[_isen]): 0 );
        if( (_vFCT.size()>=ns || _istg==ns-1)
         && !_FCT_D_SEN( _pos_fct, _isen, _t, fp+_isen*_nf ) )
          { _END_STA(); _END_SEN(); return FATAL; }
        if( xpk && !xpk[_istg+1] ) xpk[_istg+1] = new double[_nx*_np];
        for( unsigned iy=0; iy<_nx; iy++ ) xpk[_istg+1][_isen*_nx+iy] = _Dy[iy];
      }

      // Display & record stage results
      if( options.DISPLAY >= 1 ){
        _print_interm( _t, _nx, xk[_istg+1], "x", os );
        _print_interm( _nx*_np, xpk[_istg+1], "xp", os );
      }
      if( options.RESRECORD ){
        results_sta.push_back( Results( tk[_istg+1], _nx, xk[_istg+1] ) );
        results_sen.push_back( Results( tk[_istg+1], _nx*_np, xpk[_istg+1] ) );
      }
    }

    // Bounds on final quadratures and functions
    if( options.DISPLAY >= 1 ){
      _print_interm( _nf, f, "f", os );
      _print_interm( _nf*_np, fp, "fp", os );
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

} // end namescape mc

#endif

