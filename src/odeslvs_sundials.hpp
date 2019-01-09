// Copyright (C) 2015- Benoit Chachuat & Nikola Peric, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

#ifndef MC__ODESLVS_SUNDIALS_HPP
#define MC__ODESLVS_SUNDIALS_HPP

#undef  MC__ODESLVS_SUNDIALS_DEBUG

#include <sstream>
#include "odeslvs_base.hpp"
#include "odeslv_sundials.hpp"

#define MC__ODESLVS_SUNDIALS_USE_BAD
#undef  MC__ODESLVS_SUNDIALS_DEBUG

namespace mc
{
//! @brief C++ class computing solutions of parametric ODEs with forward/adjoint sensitivity analysis capability using SUNDIALS and MC++.
////////////////////////////////////////////////////////////////////////
//! mc::ODESLV_SUNDIALS is a C++ class for solution of IVPs in ODEs
//! with forward/adjoint sensitivity analysis capability using the code
//! CVODES in SUNDIALS and MC++.
////////////////////////////////////////////////////////////////////////
class ODESLVS_SUNDIALS:
  public virtual ODESLV_SUNDIALS,
  public virtual BASE_DE,
  public virtual ODESLVS_BASE
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

  //! @brief Vector storing adjoint/sensitivity trajectyories (see Options::RESRECORD)
  std::vector< std::vector< Results > > results_sen;

 //! @brief Propagate states and state-sensitivities forward in time through every time stages
  STATUS states_FSA
    ( const double*p, double**xk=0, double*f=0, double**xpk=0, double*fp=0, std::ostream&os=std::cout );

  //! @brief Propagate states and adjoints forward and backward in time through every time stages
  STATUS states_ASA
    ( const double*p, double**xk=0, double*f=0, double**lk=0, double*fp=0, std::ostream&os=std::cout );

  //! @brief Record state and sensitivity trajectories in files <a>obndsta</a> and <a>obndsa</a>, with accuracy of <a>iprec</a> digits
  void record
    ( std::ofstream&obndsta, std::ofstream*obndsen, const unsigned iprec=5 ) const
    { this->ODESLV_SUNDIALS::record( obndsta, iprec );
      for( unsigned isen=0; isen<results_sen.size(); ++isen )
        this->ODESLV_BASE::_record( obndsen[isen], results_sen[isen], iprec ); }

  //! @brief Record state trajectories in files <a>obndsta</a>, with accuracy of <a>iprec</a> digits
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
: //ODESLVS_BASE(), ODESLV_SUNDIALS(),
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
  results_sen.resize( _nf );
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
//!  - <a>p</a>  [input]  parameter values
//!  - <a>xk</a> [output] state + quadrature values at stage times
//!  - <a>f</a>  [output] function values
//!  - <a>lk</a> [output] adjoint + adjoint quadrature values at stage times
//!  - <a>fp</a> [output] function derivatives
//!  - <a>os</a> [input]  output stream [default: std::cout]
//! .
//! The return value is the status.
inline typename ODESLVS_SUNDIALS::STATUS
ODESLVS_SUNDIALS::states_ASA
( const double*p, double**xk, double*f, double**lk, double*fp, std::ostream&os )
{
  // Compute state bounds and store intermediate results
  STATUS flag = NORMAL;
  flag = ODESLV_SUNDIALS::_states( p, xk, f, true, os);
  if( flag != NORMAL ) return flag;

  // Nothing to do if no functions are defined
  if( !_nf ) return NORMAL;

  try{
    // Initialize adjoint integration
    if( !_INI_ASA( p )) return FATAL;
    _t = _dT[_nsmax];
    const unsigned NSTEP = options.RESRECORD? options.RESRECORD: 1;

    // Terminal adjoint & quadrature values
    if( lk && !lk[_nsmax] ) lk[_nsmax] = new double[(_nx+_np)*_nf];
    _pos_fct = ( _vFCT.size()>=_nsmax? _nsmax-1:0 );
    for( _ifct=0; _ifct < _nf; _ifct++ ){
      if( !_TC_SET_ASA( _pos_fct, _ifct )
       || !_TC_D_SEN( _t, _vec_sta[_nsmax].data(), NV_DATA_S(_Ny[_ifct]) )
       || ( _Nyq && _Nyq[_ifct] && !_TC_D_QUAD_ASA( NV_DATA_S(_Nyq[_ifct]) ) ) )
        { _END_SEN(); return FATAL; }
      _GET_D_SEN( NV_DATA_S(_Ny[_ifct]), _np, _Nyq? NV_DATA_S(_Nyq[_ifct]): 0 );
      for( unsigned iq=0; iq<_np; iq++ )
        _Dfp[iq*_nf+_ifct] = _Dyq[iq];
      // Display / record / return adjoint terminal values
      if( options.DISPLAY >= 1 ){
        std::ostringstream ol; ol << " l[" << _ifct << "]";
        if( !_ifct ) _print_interm( _dT[_nsmax], _nx, _Dy, ol.str(), os );
        else        _print_interm( _nx, _Dy, ol.str(), os );
        std::ostringstream oq; oq << " qp[" << _ifct << "]";
        _print_interm( _np, _Dyq, oq.str(), os );
      }
      if( options.RESRECORD )
        results_sen[_ifct].push_back( Results( _t, _nx, _Dy, _np, _Dyq ) );
      for( unsigned iy=0; lk && iy<_ny+_np; iy++ )
        lk[_nsmax][_ifct*(_nx+_np)+iy] = iy<_ny? _Dy[iy]: _Dyq[iy-_ny];
    }

    // Initialization of adjoint integration
    for( _ifct=0; _ifct < _nf; _ifct++ )
      if( !_INI_CVODES( MC_CVASARHSD__, MC_CVASAQUADD__, MC_CVASAJACD__,
        _ifct, _indexB[_ifct], _iusrB[_ifct] ) )
        { _END_SEN(); return FATAL;}

    // Integrate adjoint ODEs through each stage using SUNDIALS
    pODESLVS = this;
    for( _istg=_nsmax; _istg>0; _istg-- ){

      // Update list of operations in RHSADJ and QUADADJ
      _pos_rhs  = ( _vRHS.size() <=1? 0: _istg-1 );
      _pos_quad = ( _vQUAD.size()<=1? 0: _istg-1 );
      if( !ODESLV_BASE::_RHS_D_SET( _pos_rhs, _pos_quad )
       || !_RHS_SET_ASA( _pos_rhs, _pos_quad, _pos_fct )
       || !_RHS_D_SET( _nf, _np ) )
        { _END_SEN(); return FATAL; }

      // Propagate bounds backward to previous stage time
      const double TSTEP = ( _t - _dT[_istg-1] ) / NSTEP;
      double TSTOP = _t-TSTEP;
      for( unsigned k=0; k<NSTEP; k++, TSTOP-=TSTEP ){
        if( k+1 == NSTEP ) TSTOP = _dT[_istg-1];
        _cv_flag = CVodeB( _cv_mem, TSTOP, CV_NORMAL );
        if( _check_cv_flag(&_cv_flag, "CVodeB", 1) )
          { _END_SEN(); return FATAL; }

        // intermediate record
        if( options.RESRECORD ){
          for( _ifct=0; _ifct < _nf; _ifct++ ){
            _cv_flag = CVodeGetB( _cv_mem, _indexB[_ifct], &_t, _Ny[_ifct]);
            if( _check_cv_flag( &_cv_flag, "CVodeGetB", 1) )
              { _END_SEN(); return FATAL; }
#ifdef MC__ODESLVS_SUNDIALS_DEBUG
            std::cout << "Adjoint #" << _ifct << ": " << _t << std::endl;
            for( unsigned iy=0; iy<_ny; iy++ )
              std::cout << "_Ny" << _ifct << "[" << iy << "] = " << NV_Ith_S(_Ny[_ifct],iy) << std::endl;
#endif
            _cv_flag = CVodeGetQuadB( _cv_mem, _indexB[_ifct], &_t, _Nyq[_ifct]);
            if( _check_cv_flag( &_cv_flag, "CVodeGetQuadB", 1) )
              { _END_SEN(); return FATAL; }
#ifdef MC__ODESLVS_SUNDIALS_DEBUG
            for( unsigned ip=0; ip<_np; ip++ )
              std::cout << "_Nyq" << _ifct << "[" << ip << "] = " << NV_Ith_S(_Nyq[_ifct],ip) << std::endl;
#endif
            results_sen[_ifct].push_back( Results( _t, _nx, NV_DATA_S(_Ny[_ifct]), _np, _Nyq && _Nyq[_ifct]? NV_DATA_S(_Nyq[_ifct]): 0 ) );
          }
        }
      }
      for( _ifct=0; _ifct < _nf; _ifct++ ){
        void *cv_memB = CVodeGetAdjCVodeBmem(_cv_mem, _indexB[_ifct] );
        long int nstpB;
        _cv_flag = CVodeGetNumSteps( cv_memB, &nstpB );
        stats_sen.numSteps += nstpB;
#ifdef MC__ODESLVS_SUNDIALS_DEBUG
        std::cout << "Number of steps for adjoint #" << _ifct << ": " 
                  << nstpB << std::endl;
#endif
      }

      // states/adjoints/quadratures at stage time
#ifdef MC__ODESLVS_SUNDIALS_DEBUG
      for( unsigned ix=0; ix<_nx; ix++ )
        std::cout << "_vec_sta = " << _vec_sta[_istg-1][ix] << std::endl;
#endif
      if( lk && !lk[_istg-1] ) lk[_istg-1] = new double[(_nx+_np)*_nf];
      for( _ifct=0; _ifct < _nf; _ifct++ ){
        _cv_flag = CVodeGetB( _cv_mem, _indexB[_ifct], &_t, _Ny[_ifct]);
        if( _check_cv_flag( &_cv_flag, "CVodeGetB", 1) )
          { _END_SEN(); return FATAL; }
#ifdef MC__ODESLVS_SUNDIALS_DEBUG
        for( unsigned iy=0; iy<_ny; iy++ )
          std::cout << "_Ny" << _ifct << "[" << iy << "] = " << NV_Ith_S(_Ny[_ifct],iy) << std::endl;
#endif
        _cv_flag = CVodeGetQuadB( _cv_mem, _indexB[_ifct], &_t, _Nyq[_ifct]);
        if( _check_cv_flag( &_cv_flag, "CVodeGetQuadB", 1) )
          { _END_SEN(); return FATAL; }
#ifdef MC__ODESLVS_SUNDIALS_DEBUG
        for( unsigned ip=0; ip<_np; ip++ )
          std::cout << "_Nyq" << _ifct << "[" << ip << "] = " << NV_Ith_S(_Nyq[_ifct],ip) << std::endl;
#endif
        // Add function contribution to adjoint bounds (discontinuities)
        if( _istg > 1  ){
          _pos_fct = ( _vFCT.size()>=_nsmax? _istg-1:0 );
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

        // Display / record / return adjoint terminal values
        if( options.DISPLAY >= 1 ){
          std::ostringstream ol; ol << " l[" << _ifct << "]";
          if( !_ifct ) _print_interm( _dT[_istg-1], _nx, _Dy, ol.str(), os );
          else         _print_interm( _nx, _Dy, ol.str(), os );
          std::ostringstream oq; oq << " qp[" << _ifct << "]";
          _print_interm( _np, _Dyq, oq.str(), os );
        }
        if( options.RESRECORD )
          results_sen[_ifct].push_back( Results( _t, _nx, _Dy, _np, _Dyq ) );
        for( unsigned iy=0; lk && iy<_ny+_np; iy++ )
          lk[_istg-1][_ifct*(_nx+_np)+iy] = iy<_ny? _Dy[iy]: _Dyq[iy-_ny];

        // Keep track of function derivatives
        for( unsigned iq=0; iq<_np; iq++ ) _Dfp[iq*_nf+_ifct] = _Dyq[iq];
      }
    }

    // Display / return function derivatives
    for( unsigned i=0; fp && i<_nf*_np; i++ ) fp[i] = _Dfp[i];
    if( options.DISPLAY >= 1 ){
      for( unsigned iq=0; iq<_np; iq++ ){
        std::ostringstream ofp; ofp << " fp[" << iq << "]";
        _print_interm( _nf, _Dfp.data()+iq*_nf, ofp.str(), os );
      }
    }
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
    if( _Ny )   N_VDestroyVectorArray_Serial( _Ny,  _nvec );
    _Ny = 0;
    if( _Nyq )  N_VDestroyVectorArray_Serial( _Nyq, _nvec );
    _Nyq = 0;
    _nvec = _np;
    _Ny  = N_VCloneVectorArray_Serial( _np, _Nx );
    if( _nq ) _Nyq = N_VCloneVectorArray_Serial( _np, _Nq );
  }
  for( unsigned i=0; i<_np; i++ ){
    if( !_Ny[i] || NV_LENGTH_S( _Ny[i] ) != _nx ){
      if( _Ny[i] ) N_VDestroy_Serial( _Ny[i] );
      _Ny[i] = N_VNew_Serial( _nx );
    }
    if( _Nyq && (!_Nyq[i] || NV_LENGTH_S( _Nyq[i]) ) != _nq ){
      if( _Nyq[i] ) N_VDestroy_Serial( _Nyq[i] );
      _Nyq[i] = _nq? N_VNew_Serial( _nq ): 0;
    }
  }

  // Reset result record and statistics
  results_sen.clear();
  results_sen.resize( _np );
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
  std::cout << "y[" << is << "]:\n";
  for( unsigned i=0; i<NV_LENGTH_S( y ); i++ ) std::cout << NV_Ith_S( y, i ) << std::endl;
#endif
  bool flag = pODESLVS->_RHS_D_SEN( t, NV_DATA_S( x ), NV_DATA_S( y ),
    NV_DATA_S( ydot ), is );
#ifdef MC__ODESLVS_SUNDIALS_DEBUG
  std::cout << "ydot[" << is << "]:\n";
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
//! const double*p, double**xk, double*f, double**xpk, double*fp, std::ostream&os=std::cout )
//!
//! This function computes a solution to the parametric ODEs with forward
//! sensitivity analysis:
//!  - <a>p</a>   [input]  parameter values
//!  - <a>xk</a>  [output] state & quadrature values at stage times
//!  - <a>f</a>   [output] function values
//!  - <a>xpk</a> [output] state- & quadrature-sensitivity values at stage times
//!  - <a>fp</a>  [output] function derivatives
//!  - <a>os</a>  [input]  output stream [default: std::cout]
//! .
//! The return value is the status.
inline typename ODESLVS_SUNDIALS::STATUS
ODESLVS_SUNDIALS::states_FSA
( const double*p, double**xk, double*f, double**xpk, double*fp, std::ostream&os )
{
  // Check arguments
  if( !p ) return FATAL;

  try{
    // Initialize trajectory integration
    if( !ODESLV_SUNDIALS::_INI_STA( p ) 
     || !_INI_FSA( p ) ) return FATAL;
    _t = _dT[0];
    const unsigned NSTEP = options.RESRECORD? options.RESRECORD: 1;

    // Initial state/quadrature values
    if( xk && !xk[0] ) xk[0] = new double[_nx+_nq];
    if( !ODESLV_BASE::_IC_D_SET()
     || !ODESLV_BASE::_IC_D_STA( _t, NV_DATA_S( _Nx ) )
     || ( _Nq && !ODESLV_BASE::_IC_D_QUAD( NV_DATA_S( _Nq ) ) ) )
      { _END_STA(); _END_SEN(); return FATAL; }
    _GET_D_STA( NV_DATA_S(_Nx), _nq && _Nq? NV_DATA_S(_Nq): 0 );

    // Display / record / return initial results
    if( options.DISPLAY >= 1 ){
      _print_interm( _t, _nx, ODESLV_BASE::_Dx, " x", os );
      _print_interm( _nq, ODESLV_BASE::_Dq, " q", os );
    }
    if( options.RESRECORD )
      results_sta.push_back( Results( _t, _nx, NV_DATA_S(_Nx), _nq, _nq? NV_DATA_S(_Nq): 0 ) );
    for( unsigned ix=0; xk && ix<_nx+_nq; ix++ )
      xk[0][ix] = ix<_nx? ODESLV_BASE::_Dx[ix]: ODESLV_BASE::_Dq[ix-_nx];

    // Initial state/quadrature sensitivities
    if( xpk && !xpk[0] ) xpk[0] = new double[(_nx+_nq)*_np];
    for( _isen=0; _isen<_np; _isen++ ){
      if( !_IC_SET_FSA( _isen )
       || !ODESLV_BASE::_IC_D_STA( _t, NV_DATA_S(_Ny[_isen]) )
       || ( _Nyq && _Nyq[_isen] && !ODESLV_BASE::_IC_D_QUAD( NV_DATA_S(_Nyq[_isen]) ) ) ) 
        { _END_STA(); _END_SEN(); return FATAL; }
      _GET_D_SEN( NV_DATA_S(_Ny[_isen]), _nq, _nq && _Nyq? NV_DATA_S(_Nyq[_isen]): 0 );

      // Display / record / return initial results
      if( options.DISPLAY >= 1 ){
        std::ostringstream oxp; oxp << " xp[" << _isen << "]";
        _print_interm( _nx, _Dy, oxp.str(), os );
        std::ostringstream oqp; oqp << " qp[" << _isen << "]";
        _print_interm( _nq, _Dyq, oqp.str(), os );
      }
      if( options.RESRECORD )
        results_sen[_isen].push_back( Results( _t, _nx, NV_DATA_S(_Ny[_isen]), _nq, _nq? NV_DATA_S(_Nyq[_isen]):0 ) );
      for( unsigned ix=0; xpk && ix<_nx+_nq; ix++ )
        xpk[0][(_nx+_nq)*_isen+ix] = ix<_nx? _Dy[ix]: _Dyq[ix-_nx];
    }

    // Integrate ODEs through each stage using SUNDIALS
    ODESLV_SUNDIALS::pODESLV = pODESLVS = this;
    if( !ODESLV_SUNDIALS::_INI_CVODE( MC_CVRHSD__, MC_CVQUADD__, MC_CVJACD__ )
     || !_INI_CVODES( MC_CVFSARHSD__, MC_CVFSAQUADD__ ) )
      { _END_STA(); _END_SEN(); return FATAL; }

    for( _istg=0; _istg<_nsmax; _istg++ ){
      // Bounds on state discontinuities (if any) at stage times
      // and integrator reinitialization (if applicable)
      _pos_ic = ( _vIC.size()>=_nsmax? _istg:0 );
      if( _pos_ic
       && ( !ODESLV_BASE::_CC_D_SET( _pos_ic )
         || !ODESLV_BASE::_CC_D_STA( _t, NV_DATA_S( _Nx ) )
         || !ODESLV_SUNDIALS::_CC_CVODE_STA() ) )
        { _END_STA(); _END_SEN(); return FAILURE; }
      if( _istg 
       && ( ( _Nq && !ODESLV_BASE::_IC_D_QUAD( NV_DATA_S( _Nq ) ) ) // quadrature reinitialization
         || !ODESLV_SUNDIALS::_CC_CVODE_QUAD() ) )
        { _END_STA(); _END_SEN(); return FAILURE; }
      if( options.RESRECORD )
        results_sta.push_back( Results( _t, _nx, NV_DATA_S(_Nx), _nq, _nq? NV_DATA_S(_Nq): 0 ) );
      for( _isen=0; _isen<_np; _isen++ ){
        if( _pos_ic
         && ( !_CC_SET_FSA( _pos_ic, _isen )
           || !_CC_D_SEN( _t, NV_DATA_S( _Nx ), NV_DATA_S(_Ny[_isen]) )
           || ( _isen==_np-1 && !_CC_CVODES_FSA() ) ) )
            { _END_STA(); _END_SEN(); return FATAL; }
        if( _istg
         && ( ( _Nyq && _Nyq[_isen] && !ODESLV_BASE::_IC_D_QUAD( NV_DATA_S(_Nyq[_isen]) ) ) //quadrature sensitivity reinitialization
           || ( _isen==_np-1 && !_CC_CVODES_QUAD() ) ) )
            { _END_STA(); _END_SEN(); return FATAL; }
#ifdef MC__ODESLVS_SUNDIALS_DEBUG
        for( unsigned iy=0; iy<NV_LENGTH_S(_Ny[_isen]); iy++ )
          std::cout << "_Ny" << _isen << "[" << iy << "] = " << NV_Ith_S(_Ny[_isen],iy) << std::endl;
        for( unsigned iy=0; _nq && iy<NV_LENGTH_S(_Nyq[_isen]); iy++ )
          std::cout << "_Nyq" << _isen << "[" << iy << "] = " << NV_Ith_S(_Nyq[_isen],iy) << std::endl;
#endif
        if( options.RESRECORD )
          results_sen[_isen].push_back( Results( _t, _nx, NV_DATA_S(_Ny[_isen]), _nq, _nq? NV_DATA_S(_Nyq[_isen]): 0 ) );
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
      _cv_flag = CVodeSetStopTime( _cv_mem, _dT[_istg+1] );
      if( _check_cv_flag(&_cv_flag, "CVodeSetStopTime", 1) )
        { _END_STA(); return FATAL; }
/*
      while( _t < _dT[_istg+1] ){
        _cv_flag = CVode( _cv_mem, _dT[_istg+1], _Nx, &_t, CV_ONE_STEP );
        if( _check_cv_flag(&_cv_flag, "CVode", 1)
         || (options.NMAX && stats_sen.numSteps > options.NMAX) )
          throw Exceptions( Exceptions::INTERN );
        stats_sta.numSteps++;
        stats_sen.numSteps++;
*/
      const double TSTEP = ( _dT[_istg+1] - _t ) / NSTEP;
      double TSTOP = _t+TSTEP;
      for( unsigned k=0; k<NSTEP; k++, TSTOP+=TSTEP ){
        if( k+1 == NSTEP ) TSTOP = _dT[_istg+1];
        _cv_flag = CVode( _cv_mem, TSTOP, _Nx, &_t, CV_NORMAL );
        if( _check_cv_flag(&_cv_flag, "CVode", 1) )
         //|| (options.NMAX && stats_sen.numSteps > options.NMAX) )
          throw Exceptions( Exceptions::INTERN );

        // intermediate record
        if( options.RESRECORD ){
          if( _nq ){
            _cv_flag = CVodeGetQuad( _cv_mem, &_t, _Nq );
            if( _check_cv_flag(&_cv_flag, "CVodeGetQuad", 1) )
              { _END_STA(); return FATAL; }
          }
          results_sta.push_back( Results( _t, _nx, NV_DATA_S(_Nx), _nq, _nq? NV_DATA_S(_Nq): 0 ) );
          for( _isen=0; _isen<_np; _isen++ ){
            _cv_flag = CVodeGetSens1(_cv_mem, &_t, _isen, _Ny[_isen] );
            if( _check_cv_flag( &_cv_flag, "CVodeGetSens", 1) )
             { _END_STA(); _END_SEN(); return FATAL; }
            if( _nq ){
              _cv_flag = CVodeGetQuadSens1(_cv_mem, &_t, _isen, _Nyq[_isen]);
              if( _check_cv_flag( &_cv_flag, "CVodeGetQuadSens", 1) )
                { _END_STA(); _END_SEN(); return FATAL; }
            }
            results_sen[_isen].push_back( Results( _t, _nx, NV_DATA_S(_Ny[_isen]), _nq, _nq? NV_DATA_S(_Nyq[_isen]): 0 ) );
          }
        }
      }

      // Intermediate states and quadratures
      if( xk && !xk[_istg+1] ) xk[_istg+1] = new double[_nx+_nq];
      if( _nq ){
        _cv_flag = CVodeGetQuad( _cv_mem, &_t, _Nq );
        if( _check_cv_flag(&_cv_flag, "CVodeGetQuad", 1) )
          { _END_STA(); return FATAL; }
      }
      ODESLV_SUNDIALS::_GET_D_STA( NV_DATA_S(_Nx), _nq && _Nq? NV_DATA_S(_Nq): 0 );
      // Display / return stage results
      for( unsigned ix=0; xk && ix<_nx+_nq; ix++ )
        xk[_istg+1][ix] = ix<_nx? ODESLV_SUNDIALS::_Dx[ix]: ODESLV_SUNDIALS::_Dq[ix-_nx];
      if( options.DISPLAY >= 1 ){
        _print_interm( _t, _nx, ODESLV_BASE::_Dx, " x", os );
        _print_interm( _nq, ODESLV_BASE::_Dq, " q", os );
      }
      // Add intermediate function terms
      _pos_fct = ( _vFCT.size()>=_nsmax? _istg:0 );
      if( (_vFCT.size()>=_nsmax || _istg==_nsmax-1)
       && !ODESLV_BASE::_FCT_D_STA( _pos_fct, _t ) )
        { _END_STA(); _END_SEN(); return FATAL; }

      // Intermediate state and quadrature sensitivities
      if( xpk && !xpk[_istg+1] ) xpk[_istg+1] = new double[(_nx+_nq)*_np];
      for( _isen=0; _isen<_np; _isen++ ){
        _cv_flag = CVodeGetSens1(_cv_mem, &_t, _isen, _Ny[_isen] );
        if( _check_cv_flag( &_cv_flag, "CVodeGetSens", 1) )
         { _END_STA(); _END_SEN(); return FATAL; }
        if( _nq ){
          _cv_flag = CVodeGetQuadSens1(_cv_mem, &_t, _isen, _Nyq[_isen]);
          if( _check_cv_flag( &_cv_flag, "CVodeGetQuadSens", 1) )
            { _END_STA(); _END_SEN(); return FATAL; }
        }
        _GET_D_SEN( NV_DATA_S(_Nx), NV_DATA_S(_Ny[_isen]), _nq && _Nq? NV_DATA_S(_Nq): 0,
                    _nq, _nq && _Nyq[_isen]? NV_DATA_S(_Nyq[_isen]): 0 );
        // Display / return stage results
        if( options.DISPLAY >= 1 ){
          std::ostringstream oxp; oxp << " xp[" << _isen << "]";
          _print_interm( _nx, _Dy, oxp.str(), os );
          std::ostringstream oqp; oqp << " qp[" << _isen << "]";
          _print_interm( _nq, _Dyq, oqp.str(), os );
        }
        for( unsigned ix=0; xpk && ix<_nx+_nq; ix++ )
          xpk[_istg+1][(_nx+_nq)*_isen+ix] = ix<_nx? _Dy[ix]: _Dyq[ix-_nx];
        // Add intermediate function derivative terms
        if( (_vFCT.size()>=_nsmax || _istg==_nsmax-1)
         && !_FCT_D_SEN( _pos_fct, _isen, _t ) )
          { _END_STA(); _END_SEN(); return FATAL; }
      }
    }

    // Display / return function values and derivatives
    for( unsigned i=0; f && i<_nf; i++ ) f[i] = _Df[i];
    for( unsigned i=0; fp && i<_nf*_np; i++ ) fp[i] = _Dfp[i];
    if( options.DISPLAY >= 1 ){
      _print_interm( _nf, _Df.data(), " f", os );
      for( unsigned iq=0; iq<_np; iq++ ){
        std::ostringstream ofp; ofp << " fp[" << iq << "]";
        _print_interm( _nf, _Dfp.data()+iq*_nf, ofp.str(), os );
      }
    }
  }
  catch(...){
    _END_STA(); _END_SEN();
    long int nstp;
    _cv_flag = CVodeGetNumSteps( _cv_mem, &nstp );
    stats_sta.numSteps += nstp;
    stats_sen.numSteps += nstp;
    if( options.DISPLAY >= 1 ) _print_stats( stats_sen, os );
    return FAILURE;
  }

  long int nstp;
  _cv_flag = CVodeGetNumSteps( _cv_mem, &nstp );
  stats_sta.numSteps += nstp;
  stats_sen.numSteps += nstp;
#ifdef MC__ODESLVS_SUNDIALS_DEBUG
  std::cout << "number of steps: " << nstp << std::endl;
#endif

  _END_STA(); _END_SEN();
  if( options.DISPLAY >= 1 ) _print_stats( stats_sen, os );
  return NORMAL;
}

} // end namescape mc

#endif

