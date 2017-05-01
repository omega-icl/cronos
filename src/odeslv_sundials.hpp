// Copyright (C) 2015 Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

#ifndef MC__ODESLV_SUNDIALS_HPP
#define MC__ODESLV_SUNDIALS_HPP

#undef  MC__ODESLV_SUNDIALS_DEBUG

#include "base_sundials.hpp"
#include "odeslv_base.hpp"

namespace mc
{
//! @brief C++ class computing solutions of parametric ODEs using SUNDIALS and MC++.
////////////////////////////////////////////////////////////////////////
//! mc::ODESLV_SUNDIALS is a C++ class for solution of IVPs in ODEs
//! using the code CVODES in SUNDIALS and MC++.
////////////////////////////////////////////////////////////////////////
class ODESLV_SUNDIALS:
  public virtual BASE_SUNDIALS,
  public virtual ODESLV_BASE,
  public virtual BASE_DE
{
  typedef int (*CVRhsFn)( realtype t, N_Vector y, N_Vector ydot, void *user_data );
  typedef int (*CVDlsDenseJacFn)( long int N, realtype t, N_Vector y, N_Vector fy,
     DlsMat Jac, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3 );

 protected:
  //! @brief Pointer to the CVODE memory block
  void *_cv_mem;

  //! @brief Return flag for SUNDIALS methods
  int _cv_flag;

  //! @brief N_Vector object holding current states
  N_Vector _Nx;

  //! @brief N_Vector object holding current state quadratures
  N_Vector _Nq;

  //! @brief stepsize
  double _h;

  //! @brief position in _vIC
  unsigned _pos_ic;

  //! @brief position in _vRHS
  unsigned _pos_rhs;

  //! @brief position in _vQUAD
  unsigned _pos_quad;
  
  //! @brief position in _vFCT
  unsigned _pos_fct;

  //! @brief checkpoints for adjoint integration
  int _nchk;

  //! @brief state parameterizations at time stages for adjoint integration
  std::vector< std::vector<realtype> > _vec_sta;

  //! @brief static pointer to class
  static ODESLV_SUNDIALS *pODESLV;

public:
  /** @ingroup ODESLV
   *  @{
   */
  //! @brief Default constructor
  ODESLV_SUNDIALS
    ();

  //! @brief Virtual destructor
  virtual ~ODESLV_SUNDIALS
    ();

  //! @brief Integrator options
  struct Options:
    public BASE_SUNDIALS::Options
  {
    //! @brief Constructor
    Options():
      BASE_SUNDIALS::Options(), DISPLAY(1), RESRECORD(false)
      { JACAPPROX = CV_DENSE; }
    //! @brief Assignment operator
    template <typename OPT> Options& operator=
      ( OPT&options ){
        BASE_SUNDIALS::Options::operator=(options);
        DISPLAY   = options.DISPLAY;
        RESRECORD = options.RESRECORD;
        return *this;
      }
    //! @brief Display level (default: 1)
    int DISPLAY;
    //! @brief Whether or not to record results (default: false)
    unsigned RESRECORD;
  } options;

  //! @brief Structure for setting up storing the solver exceptions
  class Exceptions
  {
  public:
    //! @brief Enumeration type for ODESLV_SUNDIALS exception handling
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
        return "ODESLV_SUNDIALS::Exceptions  Internal error";
       }
    }
  private:
    TYPE _ierr;
  };

  //! @brief Vector storing results (upon request only)
  std::vector< Results > results_sta;

  //! @brief Statistics for state integration
  Stats stats_sta;

  //! @brief Computes solution of parametric ODEs
  STATUS states
    ( const double*p, double**xk=0, double*f=0, std::ostream&os=std::cout );

  //! @brief Record results in file <a>ores</a>, with accuracy of <a>iprec</a> digits
  void record
    ( std::ofstream&ores, const unsigned iprec=5 ) const
    { return _record( ores, results_sta, iprec ); }
  /** @} */

protected:

  //! @brief Function to initialize CVode memory block
  virtual bool _INI_CVODE
    ( CVRhsFn MC_CVRHS, CVRhsFn MC_CVQUAD, CVDlsDenseJacFn MC_CVJAC );

  //! @brief Function to reinitialize CVode memory block
  bool _CC_CVODE_STA
    ();

  //! @brief Function to reinitialize CVode memory block
  bool _CC_CVODE_QUAD
    ();

  //! @brief Function to finalize state integration
  void _END_STA
    ();

  //! @brief Function to initialize state integration
  bool _INI_STA
    ( const double*p );

  //! @brief Function to initialize state integration
  bool _INI_STA
    ();

  //! @brief Static wrapper to function to calculate the ODEs RHS derivatives
  static int MC_CVRHSD__
    ( realtype t, N_Vector Nx, N_Vector Nxdot, void *user_data );

  //! @brief Static wrapper to function to calculate the quadrature RHS derivatives
  static int MC_CVQUADD__
    ( realtype t, N_Vector Nx, N_Vector Nqdot, void *user_data );

  //! @brief Static wrapper to function to calculate the ODEs RHS Jacobian
  static int MC_CVJACD__
    ( long int N, realtype t, N_Vector y, N_Vector fy, DlsMat Jac, void *user_data,
      N_Vector tmp1, N_Vector tmp2, N_Vector tmp3 );

  //! @brief Solve parametric ODEs forward in time through every time stages
  virtual STATUS _states
    ( const double*p, double**xk, double*f, const bool store, std::ostream&os );

  //! @brief Private methods to block default compiler methods
  ODESLV_SUNDIALS(const ODESLV_SUNDIALS&);
  ODESLV_SUNDIALS& operator=(const ODESLV_SUNDIALS&);
};

ODESLV_SUNDIALS* ODESLV_SUNDIALS::pODESLV = 0;

inline
ODESLV_SUNDIALS::ODESLV_SUNDIALS
()
: //BASE_DE(), BASE_SUNDIALS(), ODESLV_BASE(),
  _cv_mem(0), _cv_flag(0),
  _Nx(0), _Nq(0)
{}

inline
ODESLV_SUNDIALS::~ODESLV_SUNDIALS
()
{
  if( _Nx ) N_VDestroy_Serial( _Nx );
  if( _Nq ) N_VDestroy_Serial( _Nq );
  if( _cv_mem ) CVodeFree( &_cv_mem );
}

inline bool
ODESLV_SUNDIALS::_INI_CVODE
( CVRhsFn MC_CVRHS, CVRhsFn MC_CVQUAD, CVDlsDenseJacFn MC_CVJAC )
{
  // Reset CVode memory block
  if( _cv_mem ) CVodeFree( &_cv_mem );
  switch( options.INTMETH ){
   // Non-stiff case: create CVode memory block for the ADAMS method
   // and use of a functional iteration
   case Options::MSADAMS: default:
    _cv_mem = CVodeCreate( CV_ADAMS, CV_FUNCTIONAL );
    if( _check_cv_flag((void *)_cv_mem, "CVodeCreate", 0) ) return false;
    break;

   // Stiff case: create CVode memory block for the BDF method
   // and use of a Newton iteration
   case Options::MSBDF:
    _cv_mem = CVodeCreate( CV_BDF, CV_NEWTON );
    if( _check_cv_flag((void *)_cv_mem, "CVodeCreate", 0) ) return false;
    break;
  }

  // Initialize CVode memory and specify right hand side function,
  // current time _t, and current state _Nx
  _cv_flag = CVodeInit( _cv_mem, MC_CVRHS, _t, _Nx );
  if( _check_cv_flag(&_cv_flag, "CVodeInit", 1) ) return false;

   // Specify the Jacobian approximation and linear solver
  if( options.INTMETH == Options::MSBDF ){ // Applies to Newton iteration only
    switch( options.JACAPPROX ){
     case Options::CV_DIAG: default:
       _cv_flag = CVDiag( _cv_mem );
       if( _check_cv_flag(&_cv_flag, "CVDiag", 1)) return false;
       break;

     case Options::CV_LAPACKDENSE:
       //_cv_flag = CVLapackDense( _cv_mem, NV_LENGTH_S( _Nx ) );
       //if( _check_cv_flag(&_cv_flag, "CVLapackDense", 1)) return false;
       //_cv_flag = CVDlsSetDenseJacFn( _cv_mem, NULL );
       //if ( _check_cv_flag(&_cv_flag, "CVDlsSetDenseJacFn", 1) ) return false;
       //break;

     case Options::CV_DENSE:
       _cv_flag = CVDense( _cv_mem, NV_LENGTH_S( _Nx ) );
       if( _check_cv_flag(&_cv_flag, "CVDense", 1)) return false;
       _cv_flag = CVDlsSetDenseJacFn( _cv_mem, MC_CVJAC );
       if ( _check_cv_flag(&_cv_flag, "CVDlsSetDenseJacFn", 1) ) return false;
       break;
    }
  }

  // Specify the relative and absolute tolerances for states
  _cv_flag = CVodeSStolerances( _cv_mem, options.RTOL, options.ATOL );
  if( _check_cv_flag(&_cv_flag, "CVodeSStolerances", 1) ) return false;

  // Set maximum number of error test failures
  _cv_flag = CVodeSetMaxErrTestFails( _cv_mem, options.MAXFAIL );
  if ( _check_cv_flag(&_cv_flag, "CVodeSetMaxErrTestFails", 1) ) return false;

  // Set maximum number of error test failures
  _cv_flag = CVodeSetMaxConvFails( _cv_mem, options.MAXFAIL );
  if ( _check_cv_flag(&_cv_flag, "CVodeSetMaxConvFails", 1) ) return false;

  // Specify minimum stepsize
  _cv_flag = CVodeSetMinStep( _cv_mem, options.HMIN>0.? options.HMIN:0. );
  if( _check_cv_flag(&_cv_flag, "CVodeSetMinStep", 1) ) return false;

  // Specify maximum stepsize
  _cv_flag = CVodeSetMaxStep( _cv_mem, options.HMAX>0.? options.HMAX: 0. );
  if( _check_cv_flag(&_cv_flag, "CVodeSetMaxStep", 1) ) return false;

  // Specify maximum number of steps between two stage times
  _cv_flag = CVodeSetMaxNumSteps( _cv_mem, options.NMAX );
  if( _check_cv_flag(&_cv_flag, "CVodeSetMaxNumSteps", 1) ) return false;

  // Initialize the integrator memory for the quadrature variables
  if( !_Nq ) return true;
  _cv_flag = CVodeQuadInit( _cv_mem, MC_CVQUAD, _Nq );
  if( _check_cv_flag(&_cv_flag, "CVodeQuadInit", 1) ) return false;

  // Specify whether or not to perform error control on quadrature
  _cv_flag = CVodeSetQuadErrCon( _cv_mem, options.QERR );
  if( _check_cv_flag(&_cv_flag, "CVodeSetQuadErrCon", 1) ) return false;

  // Specify the relative and absolute tolerances for quadratures
  _cv_flag = CVodeQuadSStolerances( _cv_mem, options.RTOL, options.ATOL );
  if( _check_cv_flag(&_cv_flag, "CVodeQuadSStolerances", 1) ) return false;

  // Allocate memory for adjoint sensitivity analysis (if applicable)
  //_cv_flag = CVodeAdjInit( _cv_mem, options.ASACHKPT, options.ASAINTERP );
  //if( _check_cv_flag(&_cv_flag, "CVodeAdjInit", 1) ) return false;

  return true;
}

inline bool
ODESLV_SUNDIALS::_CC_CVODE_STA
()
{
  // Reinitialize CVode memory block for current time _t and current state _Nx
  _cv_flag = CVodeReInit( _cv_mem, _t, _Nx );
  if( _check_cv_flag(&_cv_flag, "CVodeReInit", 1) ) return false;

  return true;
}

inline bool
ODESLV_SUNDIALS::_CC_CVODE_QUAD
()
{
  // Reinitialize CVode memory block for current quarature _Nq
  if( !_Nq ) return true;
  _cv_flag = CVodeQuadReInit( _cv_mem, _Nq );
  if( _check_cv_flag(&_cv_flag, "CVodeQuadReInit", 1) ) return false;

  return true;
}

inline void
ODESLV_SUNDIALS::_END_STA()
{
  // Get final CPU time
  _final_stats( stats_sta );
}

inline bool
ODESLV_SUNDIALS::_INI_STA
( const double*p )
{
  // Initialize bound propagation
  if( !_INI_D_STA( p ) || !_INI_STA() )
    return false;
  return true;
}

inline bool
ODESLV_SUNDIALS::_INI_STA
()
{
  // Set SUNDIALS state/quadrature arrays
  if( !_Nx || NV_LENGTH_S( _Nx ) != _nx ){
    if( _Nx ) N_VDestroy_Serial( _Nx );
    _Nx  = N_VNew_Serial( _nx );
  }
  if( !_Nq || NV_LENGTH_S( _Nq ) != _nq ){
    if( _Nq ) N_VDestroy_Serial( _Nq );
    _Nq  = _nq? N_VNew_Serial( _nq ): 0;
  }

  // Initialize state parameterization at time stages
  _vec_sta.clear();

  // Reset result record and statistics
  results_sta.clear();
  _init_stats( stats_sta );

  return true;
}

inline int
ODESLV_SUNDIALS::MC_CVRHSD__
( realtype t, N_Vector y, N_Vector ydot, void *user_data )
{
  ODESLV_SUNDIALS *pODESLV = ODESLV_SUNDIALS::pODESLV;
#ifdef MC__ODESLV_SUNDIALS_DEBUG
  std::cout << "@t=" << t << "\nx:\n";
  for( unsigned i=0; i<NV_LENGTH_S( y ); i++ ) std::cout << NV_Ith_S( y, i ) << std::endl;
#endif
  bool flag = pODESLV->_RHS_D_STA( t, NV_DATA_S( y ), NV_DATA_S( ydot ) );
#ifdef MC__ODESLV_SUNDIALS_DEBUG
  std::cout << "xdot:\n";
  for( unsigned i=0; i<NV_LENGTH_S( ydot ); i++ ) std::cout << NV_Ith_S( ydot, i ) << std::endl;
  { int dum; std::cin >> dum; }
#endif
  ODESLV_SUNDIALS::pODESLV = pODESLV;
  pODESLV->stats_sta.numRHS++;
  return( flag? 0: -1 );
}

inline int
ODESLV_SUNDIALS::MC_CVQUADD__
( realtype t, N_Vector y, N_Vector qdot, void *user_data )
{
  ODESLV_SUNDIALS *pODESLV = ODESLV_SUNDIALS::pODESLV;
  bool flag = pODESLV->_RHS_D_QUAD( t, NV_DATA_S( y ), NV_DATA_S( qdot ) );
  ODESLV_SUNDIALS::pODESLV = pODESLV;
  return( flag? 0: -1 );
}

inline int
ODESLV_SUNDIALS::MC_CVJACD__
( long int N, realtype t, N_Vector y, N_Vector ydot, DlsMat Jac,
  void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3 )
{
  ODESLV_SUNDIALS *pODESLV = ODESLV_SUNDIALS::pODESLV;
  bool flag = pODESLV->_JAC_D_STA( t, NV_DATA_S( y ), Jac->cols );
  ODESLV_SUNDIALS::pODESLV = pODESLV;
  pODESLV->stats_sta.numJAC++; // increment JAC counter
  return( flag? 0: -1 );
}

inline typename ODESLV_SUNDIALS::STATUS
ODESLV_SUNDIALS::_states
( const double*p, double**xk, double*f, const bool store, std::ostream&os )
{
  // Check size
  if( !p ) return FATAL;

  try{
    // Initialize trajectory integration
    if( !_INI_STA( p ) ) return FATAL;
    _t = _dT[0];
    const unsigned NSTEP = options.RESRECORD? options.RESRECORD: 1;

    // Initial state/quadrature values
    if( !_IC_D_SET()
     || !_IC_D_STA( _t, NV_DATA_S( _Nx ) )
     || (_Nq && !_IC_D_QUAD( NV_DATA_S( _Nq ) ) ) )
      { _END_STA(); return FATAL; }
    _GET_D_STA( NV_DATA_S(_Nx), _nq && _Nq? NV_DATA_S(_Nq): 0 );

    // Store full state at initial time
    if( store ){
      realtype*vsta = NV_DATA_S(_Nx);
      unsigned lsta = NV_LENGTH_S(_Nx);
      _vec_sta.push_back( std::vector<realtype>( vsta, vsta+lsta ) );
    }

    // Display / record / return initial results
    if( options.DISPLAY >= 1 ){
      _print_interm( _t, _nx, _Dx, " x", os );
      _print_interm( _nq, _Dq, " q", os );
    }
    if( options.RESRECORD )
      results_sta.push_back( Results( _t, _nx, NV_DATA_S(_Nx), _nq, _nq? NV_DATA_S(_Nq): 0 ) );
    if( xk ){
      if( !xk[0] ) xk[0] = new double[_nx+_nq];
      for( unsigned ix=0; ix<_nx+_nq; ix++ )
        xk[0][ix] = ix<_nx? _Dx[ix]: _Dq[ix-_nx];
    }

    // Integrate ODEs through each stage using SUNDIALS
    pODESLV = this;
    if( !_INI_CVODE( MC_CVRHSD__, MC_CVQUADD__, MC_CVJACD__ ) )
      { _END_STA(); return FATAL; }

    for( _istg=0; _istg<_nsmax; _istg++ ){
      // State discontinuities (if any) at stage times
      // and integrator reinitialization (if applicable)
      _pos_ic = ( _vIC.size()>=_nsmax? _istg:0 );
      if( _pos_ic
       && ( !_CC_D_SET( _pos_ic )
         || !_CC_D_STA( _t, NV_DATA_S( _Nx ) )
         || !_CC_CVODE_STA() ) )
        { _END_STA(); return FAILURE; }
      if( _istg 
       && ( ( _Nq && !_IC_D_QUAD( NV_DATA_S( _Nq ) ) ) // quadrature reinitialization
         || !_CC_CVODE_QUAD() ) )
        { _END_STA(); return FAILURE; }
      if( options.RESRECORD )
        results_sta.push_back( Results( _t, _nx, NV_DATA_S(_Nx), _nq, _nq? NV_DATA_S(_Nq): 0 ) );

      // update list of operations in RHS, JAC and QUAD
      _pos_rhs  = ( _vRHS.size()<=1? 0: _istg );
      _pos_quad = ( _vQUAD.size()<=1? 0: _istg );
      if( (!_istg || _pos_rhs || _pos_quad)
        && !_RHS_D_SET( _pos_rhs, _pos_quad ) )
        { _END_STA(); return FATAL; }

      // integrate till end of time stage
      _cv_flag = CVodeSetStopTime( _cv_mem, _dT[_istg+1] );
      if( _check_cv_flag(&_cv_flag, "CVodeSetStopTime", 1) )
        { _END_STA(); return FATAL; }
/*
      while( _t < _dT[_istg+1] ){
        if( !store )
          _cv_flag = CVode( _cv_mem, _dT[_istg+1], _Nx, &_t, CV_ONE_STEP );
        else
          _cv_flag = CVodeF( _cv_mem, _dT[_istg+1], _Nx, &_t, CV_ONE_STEP, &_nchk );
        if( _check_cv_flag(&_cv_flag, store?"CVodeF":"CVode", 1)
         || (options.NMAX && stats_sta.numSteps > options.NMAX) )
          throw Exceptions( Exceptions::INTERN );
        stats_sta.numSteps++;
      }
*/
      const double TSTEP = ( _dT[_istg+1] - _t  ) / NSTEP;
      double TSTOP = _t+TSTEP;
      for( unsigned k=0; k<NSTEP; k++, TSTOP+=TSTEP ){
        if( k+1 == NSTEP ) TSTOP = _dT[_istg+1];
        if( !store )
          _cv_flag = CVode( _cv_mem, TSTOP, _Nx, &_t, CV_NORMAL );
        else
          _cv_flag = CVodeF( _cv_mem, TSTOP, _Nx, &_t, CV_NORMAL, &_nchk );
        if( _check_cv_flag(&_cv_flag, store?"CVodeF":"CVode", 1) )
         //|| (options.NMAX && stats_sta.numSteps > options.NMAX) )
          throw Exceptions( Exceptions::INTERN );

        // intermediate record
        if( options.RESRECORD ){
          if( _nq ){
            _cv_flag = CVodeGetQuad( _cv_mem, &_t, _Nq );
            if( _check_cv_flag(&_cv_flag, "CVodeGetQuad", 1) )
              { _END_STA(); return FATAL; }
          }
          results_sta.push_back( Results( _t, _nx, NV_DATA_S(_Nx), _nq, _nq? NV_DATA_S(_Nq): 0 ) );
        }
      }

      // Store full state at stage time
      if( store ){
        realtype*vsta = NV_DATA_S(_Nx);
        unsigned lsta = NV_LENGTH_S(_Nx);
        _vec_sta.push_back( std::vector<realtype>( vsta, vsta+lsta ) );
      }

      // Intermediate states and quadratures
      if( _nq ){
        _cv_flag = CVodeGetQuad( _cv_mem, &_t, _Nq );
        if( _check_cv_flag(&_cv_flag, "CVodeGetQuad", 1) )
          { _END_STA(); return FATAL; }
      }
     _GET_D_STA( NV_DATA_S(_Nx), _nq && _Nq? NV_DATA_S(_Nq): 0 );

      // Display / return stage results
      if( xk ){
        if( !xk[_istg+1] ) xk[_istg+1] = new double[_nx+_nq];
        unsigned ix = 0;
        for( ; ix<_nx; ix++ ) xk[_istg+1][ix] = _Dx[ix];
        for( ; ix<_nx+_nq; ix++ ) xk[_istg+1][ix] = _Dq[ix-_nx];
      }
      if( options.DISPLAY >= 1 ){
        _print_interm( _t, _nx, _Dx, " x", os );
        _print_interm( _nq, _Dq, " q", os );
      }

      // Add intermediate function terms
      _pos_fct = ( _vFCT.size()>=_nsmax? _istg:0 );
      if( (_vFCT.size()>=_nsmax || _istg==_nsmax-1) && !_FCT_D_STA( _pos_fct, _t ) )
        { _END_STA(); return FATAL; }
    }
#ifdef MC__ODESLV_SUNDIALS_DEBUG
    if( store ) std::cout << "number of checkpoints: " << _nchk << std::endl;
#endif

    // Display / return function values
    for( unsigned i=0; f && i<_nf; i++ ) f[i] = _Df[i];
    if( options.DISPLAY >= 1 ) _print_interm( _nf, _Df, " f", os );
  }
  catch(...){
    _END_STA();
    long int nstp;
    _cv_flag = CVodeGetNumSteps( _cv_mem, &nstp );
    stats_sta.numSteps += nstp;
    if( options.DISPLAY >= 1 ) _print_stats( stats_sta, os );
    return FAILURE;
  }

  long int nstp;
  _cv_flag = CVodeGetNumSteps( _cv_mem, &nstp );
  stats_sta.numSteps += nstp;
#ifdef MC__ODESLV_SUNDIALS_DEBUG
  std::cout << "number of steps: " << nstp << std::endl;
#endif

  _END_STA();
  if( options.DISPLAY >= 1 ) _print_stats( stats_sta, os );
  return NORMAL;
}

//! @fn inline typename ODESLV_SUNDIALS::STATUS ODESLV_SUNDIALS::states(
//! const double*p, double**xk, double*q, double*f, std::ostream&os=std::cout )
//!
//! This function computes a solution to the parametric ODEs:
//!   - <a>p</a>  [input]  parameter values
//!   - <a>xk</a> [output] state & quadrature values at stage times
//!   - <a>f</a>  [output] function values
//!   - <a>os</a> [input]  output stream [default: std::cout]
//! .
//! The return value is the status.
inline typename ODESLV_SUNDIALS::STATUS
ODESLV_SUNDIALS::states
( const double*p, double**xk, double*f, std::ostream&os )
{
  return _states( p, xk, f, false, os );
}

} // end namescape mc

#endif

