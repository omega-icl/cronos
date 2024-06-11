// Copyright (C) Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

#ifndef MC__ODESLV_CVODES_HPP
#define MC__ODESLV_CVODES_HPP

#undef  MC__ODESLV_CVODES_DEBUG

#include "base_cvodes.hpp"
#include "odeslv_base.hpp"

namespace mc
{
//! @brief C++ class computing solutions of parametric ODEs using SUNDIALS/CVODES and MC++.
////////////////////////////////////////////////////////////////////////
//! mc::ODESLV_CVODES is a C++ class for solution of IVPs in ODEs
//! using the code CVODES in SUNDIALS and MC++.
////////////////////////////////////////////////////////////////////////
template <typename... ExtOps>
class ODESLV_CVODES
: public virtual BASE_CVODES
, public virtual ODESLV_BASE<ExtOps...>
, public virtual BASE_DE<ExtOps...>
{
 protected:
 
  using ODESLV_BASE<ExtOps...>::_print_interm;
  using ODESLV_BASE<ExtOps...>::_record;

  using ODESLV_BASE<ExtOps...>::_nsmax;
  using ODESLV_BASE<ExtOps...>::_istg;
  using ODESLV_BASE<ExtOps...>::_t;
  using ODESLV_BASE<ExtOps...>::_dT;
  using ODESLV_BASE<ExtOps...>::_nx;
  using ODESLV_BASE<ExtOps...>::_Dx;
  using ODESLV_BASE<ExtOps...>::_nq;
  using ODESLV_BASE<ExtOps...>::_Dq;
  using ODESLV_BASE<ExtOps...>::_nf;
  using ODESLV_BASE<ExtOps...>::_Df;
  using ODESLV_BASE<ExtOps...>::_vIC;
  using ODESLV_BASE<ExtOps...>::_vRHS;
  using ODESLV_BASE<ExtOps...>::_vQUAD;
  using ODESLV_BASE<ExtOps...>::_vFCT;

  using ODESLV_BASE<ExtOps...>::_INI_D_STA;
  using ODESLV_BASE<ExtOps...>::_GET_D_STA;
  using ODESLV_BASE<ExtOps...>::_IC_D_SET;
  using ODESLV_BASE<ExtOps...>::_IC_D_STA;
  using ODESLV_BASE<ExtOps...>::_IC_D_QUAD;
  using ODESLV_BASE<ExtOps...>::_CC_D_SET;
  using ODESLV_BASE<ExtOps...>::_CC_D_STA;
  using ODESLV_BASE<ExtOps...>::_RHS_D_SET;
  using ODESLV_BASE<ExtOps...>::_RHS_D_STA;
  using ODESLV_BASE<ExtOps...>::_RHS_D_QUAD;
  using ODESLV_BASE<ExtOps...>::_JAC_D_STA;
  using ODESLV_BASE<ExtOps...>::_FCT_D_STA;

  //! @brief Pointer to the CVode memory block
  void* _cv_mem;

  //! @brief Return flag for CVode methods
  int _cv_flag;

  //! @brief SUNMatrix for use in linear solves
  SUNMatrix _sun_mat;

  //! @brief SUNLinearSolver object for use by CVode
  SUNLinearSolver _sun_ls;

  //! @brief SUNNonlinearSolver object for use by CVode
  SUNNonlinearSolver _sun_nls;

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
  std::vector< std::vector< sunrealtype > > _vec_sta;

public:

  using BASE_DE<ExtOps...>::set;

  /** @ingroup ODESLV
   *  @{
   */
  typedef typename BASE_DE<ExtOps...>::STATUS STATUS;
  typedef typename ODESLV_BASE<ExtOps...>::Results Results;
  using BASE_DE<ExtOps...>::np;
  using BASE_DE<ExtOps...>::nf;

  //! @brief Default constructor
  ODESLV_CVODES
    ();

  //! @brief Virtual destructor
  virtual ~ODESLV_CVODES
    ();

  //! @brief Integrator options
  struct Options:
   public BASE_CVODES::Options
   {
    //! @brief Constructor
    Options():
      BASE_CVODES::Options(),
      DISPLAY(1), RESRECORD(0)
      {}
      //{ JACAPPROX = CV_DENSE; }
    //! @brief Assignment operator
    template <typename OPT>
    Options& operator=
      ( OPT const& options )
      {
        BASE_CVODES::Options::operator=( options );
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
    //! @brief Enumeration type for ODESLV_CVODES exception handling
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
        return "ODESLV_CVODES<ExtOps...>::Exceptions  Internal error";
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
    ( double const* p, double** xk=nullptr, double* f=nullptr,
      std::ostream& os=std::cout );

  //! @brief Setup local copy of parametric ODEs
  bool setup
    ()
    { return ODESLV_BASE<ExtOps...>::_SETUP(); }

  //! @brief Setup local copy of parametric ODEs based on IVP
  bool setup
    ( ODESLV_CVODES<ExtOps...> const& IVP )
    { return ODESLV_BASE<ExtOps...>::_SETUP( IVP ); }

  //! @brief Record results in file <a>ores</a>, with accuracy of <a>iprec</a> digits
  void record
    ( std::ofstream& ores, unsigned const iprec=5 )
    const
    { return _record( ores, results_sta, iprec ); }
  /** @} */

protected:

  //! @brief Function to initialize CVode memory block
  virtual bool _INI_CVODE
    ();

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
    ( double const* p );

  //! @brief Function to initialize state integration
  bool _INI_STA
    ();

  //! @brief Static wrapper to function to calculate the ODEs RHS derivatives
  static int MC_CVRHS__
    ( sunrealtype t, N_Vector Nx, N_Vector Nxdot, void* user_data );

  //! @brief Function to calculate the ODEs RHS derivatives
  virtual int CVRHS__
    ( sunrealtype t, N_Vector Nx, N_Vector Nxdot, void* user_data );

  //! @brief Static wrapper to function to calculate the quadrature RHS derivatives
  static int MC_CVQUAD__
    ( sunrealtype t, N_Vector Nx, N_Vector Nqdot, void* user_data );

  //! @brief Function to calculate the quadrature RHS derivatives
  virtual int CVQUAD__
    ( sunrealtype t, N_Vector Nx, N_Vector Nqdot, void* user_data );

  //! @brief Static wrapper to function to calculate the ODEs RHS Jacobian
  static int MC_CVJAC__
    ( sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix Jac, void *user_data,
      N_Vector tmp1, N_Vector tmp2, N_Vector tmp3 );

  //! @brief Function to calculate the ODEs RHS Jacobian
  virtual int CVJAC__
    ( sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix Jac, void *user_data,
      N_Vector tmp1, N_Vector tmp2, N_Vector tmp3 );

  //! @brief Solve parametric ODEs forward in time through every time stages
  virtual STATUS _states
    ( double const* p, double** xk, double* f, bool const store, std::ostream& os );

private:

  //! @brief Virtual function to calculate the adjoint ODEs RHS derivatives
  virtual int CVRHSB__
    ( sunrealtype t, N_Vector x, N_Vector y, N_Vector ydot, void* user_data )
    { throw Exceptions( Exceptions::INTERN ); }

  //! @brief Virtual function to calculate the adjoint quadrature RHS derivatives
  virtual int CVQUADB__
    ( sunrealtype t, N_Vector x, N_Vector y, N_Vector qdot, void* user_data )
    { throw Exceptions( Exceptions::INTERN ); }

  //! @brief Virtual function to calculate the adjoint ODEs RHS Jacobian
  virtual int CVJACB__
    ( sunrealtype t, N_Vector x, N_Vector y, N_Vector fB, SUNMatrix JacB,
      void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3 )
    { throw Exceptions( Exceptions::INTERN ); }

  //! @brief Virtual function to calculate the sensitivity ODEs RHS derivatives
  virtual int CVRHSF__
    ( int Ns, sunrealtype t, N_Vector x, N_Vector xdot, int is, N_Vector y,
      N_Vector ydot, void* user_data, N_Vector tmp1, N_Vector tmp2 )
    { throw Exceptions( Exceptions::INTERN ); }

  //! @brief Virtual function to calculate the sensitivity quadrature RHS derivatives
  virtual int CVQUADF__
    ( int Ns, sunrealtype t, N_Vector x, N_Vector* y, N_Vector qdot, N_Vector* qSdot, 
      void *user_data, N_Vector tmp1, N_Vector tmp2 )
    { throw Exceptions( Exceptions::INTERN ); }

  //! @brief Virtual function to calculate the sensitivity ODEs RHS Jacobian
  virtual int CVJACF__
    ( sunrealtype t, N_Vector x, N_Vector y, N_Vector fS, SUNMatrix JacS,
      void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3 )
    { throw Exceptions( Exceptions::INTERN ); }

  //! @brief Private methods to block default compiler methods
  ODESLV_CVODES( ODESLV_CVODES<ExtOps...> const& ) = delete;
  ODESLV_CVODES<ExtOps...>& operator=( ODESLV_CVODES<ExtOps...> const& ) = delete;
};

template <typename... ExtOps>
ODESLV_CVODES<ExtOps...>::ODESLV_CVODES
()
: BASE_CVODES(), BASE_DE<ExtOps...>(), ODESLV_BASE<ExtOps...>(),
  _cv_mem(nullptr), _cv_flag(0), _sun_mat(nullptr), _sun_ls(nullptr), _sun_nls(nullptr),
  _Nx(nullptr), _Nq(nullptr)
{}

template <typename... ExtOps>
ODESLV_CVODES<ExtOps...>::~ODESLV_CVODES
()
{
  if( _Nx ) N_VDestroy( _Nx );
  if( _Nq ) N_VDestroy( _Nq );
  if( _cv_mem )  CVodeFree( &_cv_mem );         /* Free CVODES memory */
  if( _sun_nls ) SUNNonlinSolFree( _sun_nls );  /* Free the nonlinear solver memory */
  if( _sun_ls )  SUNLinSolFree( _sun_ls );      /* Free the linear solver memory */
  if( _sun_mat ) SUNMatDestroy( _sun_mat );     /* Free the matrix memory */

}

template <typename... ExtOps>
bool
ODESLV_CVODES<ExtOps...>::_INI_CVODE
()
{
  // Reset CVode memory block
  if( _cv_mem ) CVodeFree( &_cv_mem );

  // Create CVode memory block for the ADAMS or BDF method
  _cv_mem = CVodeCreate( options.INTMETH, sunctx );
  if( _check_cv_flag( (void*)_cv_mem, "CVodeCreate", 0 ) ) return false;

//  // Specify error output
//  if( options.DISPLAY < 0 )
//    _cv_flag = ( SUNLogger_SetErrorFilename( _cv_mem, NULL )
//              || SUNLogger_SetWarningFilename( _cv_mem, NULL ) );
////  else{
////    _cv_flag = SUNLogger_SetErrorFilename( _cv_mem, stderr );
////    _cv_flag = SUNLogger_SetWarningFilename( _cv_mem, stderr );
////  }
//  if( _check_cv_flag( &_cv_flag, "SUNLogger_SetErrorFilename", 1 ) ) return false;

  // Initialize CVode memory and specify right hand side function,
  // current time _t, and current state _Nx
  _cv_flag = CVodeInit( _cv_mem, MC_CVRHS__, _t, _Nx );
  if( _check_cv_flag( &_cv_flag, "CVodeInit", 1 ) ) return false;

  // Specify the nonlinear solver
  if( _sun_nls ){ SUNNonlinSolFree( _sun_nls );  _sun_nls = nullptr; } /* Free the nonlinear solver memory */
  if( _sun_ls ) { SUNLinSolFree( _sun_ls );      _sun_ls  = nullptr; } /* Free the linear solver memory */
  if( _sun_mat ){ SUNMatDestroy( _sun_mat );     _sun_mat = nullptr; } /* Free the matrix memory */
  switch( options.NLINSOL ){
   // Fixed point nonlinear solver
   case Options::FIXEDPOINT:
    _sun_nls = SUNNonlinSol_FixedPoint( _Nx, 0, sunctx );
    if( _check_cv_flag( (void *)_sun_nls, "SUNNonlinSol_FixedPoint", 0 ) ) return false;
    _cv_flag = CVodeSetNonlinearSolver( _cv_mem, _sun_nls );
    if( _check_cv_flag( &_cv_flag, "CVodeSetNonlinearSolver", 1 ) ) return false;
    break;
   
   // Newton nonlinear solver
   case Options::NEWTON:
    _sun_nls = SUNNonlinSol_Newton( _Nx, sunctx );
    if( _check_cv_flag( (void *)_sun_nls, "SUNNonlinSol_Newton", 0 ) ) return false;
    _cv_flag = CVodeSetNonlinearSolver( _cv_mem, _sun_nls );
    if( _check_cv_flag( &_cv_flag, "CVodeSetNonlinearSolver", 1 ) ) return false;

    // Specify the linear solver and Jacobian approximation
    switch( options.LINSOL ){
     case Options::DIAG: default:
       _cv_flag = CVDiag( _cv_mem );
       if( _check_cv_flag( &_cv_flag, "CVDiag", 1) ) return false;
       break;

     // Dense Jacobian
     case Options::DENSE:
     case Options::DENSEDQ:
       // Create dense SUNMatrix for use in linear solves
       _sun_mat = SUNDenseMatrix( NV_LENGTH_S( _Nx ), NV_LENGTH_S( _Nx ), sunctx );
       if( _check_cv_flag( (void*)_sun_mat, "SUNDenseMatrix", 0 ) ) return false;
       // Create dense SUNLinearSolver object for use by CVode
       _sun_ls = SUNLinSol_Dense( _Nx, _sun_mat, sunctx );
       if( _check_cv_flag( (void *)_sun_ls, "SUNLinSol_Dense", 0 ) ) return false;
       // Attach the matrix and linear solver
       _cv_flag = CVodeSetLinearSolver( _cv_mem, _sun_ls, _sun_mat );
       if( _check_cv_flag( &_cv_flag, "CVodeSetLinearSolver", 1 ) ) return false;
       // Set the user-supplied Jacobian routine Jac
       _cv_flag = CVodeSetJacFn( _cv_mem, options.LINSOL==Options::DENSE? MC_CVJAC__: nullptr );
       if ( _check_cv_flag( &_cv_flag, "CVodeSetJacFn", 1 ) ) return false;
       break;
/*
     case Options::SPARSE:
       // Create sparse SUNMatrix for use in linear solves
       _sun_mat = SUNSparseMatrix( NV_LENGTH_S( _Nx ), NV_LENGTH_S( _Nx ), NV_LENGTH_S( _Nx )*NV_LENGTH_S( _Nx ), CSR_MAT, sunctx);
       if( _check_cv_flag( (void*)_sun_mat, "SUNSSarseMatrix", 0 ) ) return false;
       // Create sparse SUNLinearSolver object for use by CVode
       _sun_ls = SUNLinSol_KLU( _Nx, _sun_mat, sunctx );
       if( _check_cv_flag( (void *)_sun_ls, "SUNLinSol_Dense", 0 ) ) return false;
       // Function SUNLinSol_KLUReInit(SUNLinearSolver S, SUNMatrix A, sunindextype nnz, int reinit_type)
       // needed to reinitialize memory and flag for a new factorization (symbolic and numeric) to be conducted at
       // the next solver setup call. This routine is useful in the cases where the number of nonzeroes has changed or if
       // the structure of the linear system has changed which would require a new symbolic (and numeric factorization).
       // Attach the matrix and linear solver
       _cv_flag = CVodeSetLinearSolver( _cv_mem, _sun_ls, _sun_mat );
       if( _check_cv_flag( &_cv_flag, "CVodeSetLinearSolver", 1 ) ) return false;
       // Set the user-supplied Jacobian routine Jac
       _cv_flag = CVodeSetJacFn( _cv_mem, MC_CVJAC__ );
       if ( _check_cv_flag( &_cv_flag, "CVodeSetJacFn", 1 ) ) return false;
       break;
*/
    }
  }

  // Specify the relative and absolute tolerances for states
  _cv_flag = CVodeSStolerances( _cv_mem, options.RTOL, options.ATOL );
  if( _check_cv_flag( &_cv_flag, "CVodeSStolerances", 1 ) ) return false;

  // Set maximum number of error test failures
  _cv_flag = CVodeSetMaxErrTestFails( _cv_mem, options.MAXFAIL );
  if ( _check_cv_flag( &_cv_flag, "CVodeSetMaxErrTestFails", 1 ) ) return false;

  // Set maximum number of error test failures
  _cv_flag = CVodeSetMaxConvFails( _cv_mem, options.MAXFAIL );
  if ( _check_cv_flag( &_cv_flag, "CVodeSetMaxConvFails", 1 ) ) return false;

  // Specify minimum stepsize
  _cv_flag = CVodeSetMinStep( _cv_mem, options.HMIN>0.? options.HMIN:0. );
  if( _check_cv_flag( &_cv_flag, "CVodeSetMinStep", 1 ) ) return false;

  // Specify maximum stepsize
  _cv_flag = CVodeSetMaxStep( _cv_mem, options.HMAX>0.? options.HMAX: 0. );
  if( _check_cv_flag( &_cv_flag, "CVodeSetMaxStep", 1 ) ) return false;

  // Specify maximum number of steps between two stage times
  _cv_flag = CVodeSetMaxNumSteps( _cv_mem, options.NMAX );
  if( _check_cv_flag( &_cv_flag, "CVodeSetMaxNumSteps", 1 ) ) return false;

  // Initialize the integrator memory for the quadrature variables
  if( !_Nq ) return true;
  _cv_flag = CVodeQuadInit( _cv_mem, MC_CVQUAD__, _Nq );
  if( _check_cv_flag( &_cv_flag, "CVodeQuadInit", 1 ) ) return false;

  // Specify whether or not to perform error control on quadrature
  _cv_flag = CVodeSetQuadErrCon( _cv_mem, options.QERR );
  if( _check_cv_flag( &_cv_flag, "CVodeSetQuadErrCon", 1 ) ) return false;

  // Specify the relative and absolute tolerances for quadratures
  _cv_flag = CVodeQuadSStolerances( _cv_mem, options.RTOL, options.ATOL );
  if( _check_cv_flag( &_cv_flag, "CVodeQuadSStolerances", 1 ) ) return false;

  return true;
}

template <typename... ExtOps>
bool
ODESLV_CVODES<ExtOps...>::_CC_CVODE_STA
()
{
  // Reinitialize CVode memory block for current time _t and current state _Nx
  _cv_flag = CVodeReInit( _cv_mem, _t, _Nx );
  if( _check_cv_flag( &_cv_flag, "CVodeReInit", 1 ) ) return false;

  return true;
}

template <typename... ExtOps>
bool
ODESLV_CVODES<ExtOps...>::_CC_CVODE_QUAD
()
{
  // Reinitialize CVode memory block for current quarature _Nq
  if( !_Nq ) return true;
  _cv_flag = CVodeQuadReInit( _cv_mem, _Nq );
  if( _check_cv_flag( &_cv_flag, "CVodeQuadReInit", 1 ) ) return false;

  return true;
}

template <typename... ExtOps>
void
ODESLV_CVODES<ExtOps...>::_END_STA()
{
  // Get final CPU time
  _final_stats( stats_sta );
}

template <typename... ExtOps>
bool
ODESLV_CVODES<ExtOps...>::_INI_STA
( double const* p )
{
  // Initialize bound propagation
  if( !_INI_D_STA( p ) || !_INI_STA() )
    return false;
  return true;
}

template <typename... ExtOps>
bool
ODESLV_CVODES<ExtOps...>::_INI_STA
()
{
  // Set SUNDIALS state/quadrature arrays
  if( !_Nx || NV_LENGTH_S( _Nx ) != _nx ){
    if( _Nx ) N_VDestroy( _Nx );
    _Nx  = N_VNew_Serial( _nx, sunctx );
  }
  if( !_Nq || NV_LENGTH_S( _Nq ) != _nq ){
    if( _Nq ) N_VDestroy_Serial( _Nq );
    _Nq  = _nq? N_VNew_Serial( _nq, sunctx ): nullptr;
  }

  // Initialize state parameterization at time stages
  _vec_sta.clear();

  // Reset result record and statistics
  results_sta.clear();
  _init_stats( stats_sta );

  return true;
}

template <typename... ExtOps>
inline
int
ODESLV_CVODES<ExtOps...>::MC_CVRHS__
( sunrealtype t, N_Vector y, N_Vector ydot, void *user_data )
{
#ifdef MC__BASE_CVODES_CHECK
  //std::cout << "BASE_CVODES::PTR_BASE_CVODES: " << BASE_CVODES::PTR_BASE_CVODES
  //          << "  PTR_CVRHS: " << PTR_CVRHS << std::endl;
  assert( BASE_CVODES::PTR_BASE_CVODES != nullptr && PTR_CVRHS != nullptr);
#endif
  return (BASE_CVODES::PTR_BASE_CVODES->*PTR_CVRHS)( t, y, ydot, user_data );
}

template <typename... ExtOps>
inline
int
ODESLV_CVODES<ExtOps...>::CVRHS__
( sunrealtype t, N_Vector y, N_Vector ydot, void *user_data )
{
#ifdef MC__ODESLV_CVODES_DEBUG
  std::cout << "@t=" << t << "\nx:\n";
  for( unsigned i=0; i<NV_LENGTH_S( y ); i++ ) std::cout << NV_Ith_S( y, i ) << std::endl;
#endif
  bool flag = _RHS_D_STA( t, NV_DATA_S( y ), NV_DATA_S( ydot ) );
#ifdef MC__ODESLV_CVODES_DEBUG
  std::cout << "xdot:\n";
  for( unsigned i=0; i<NV_LENGTH_S( ydot ); i++ ) std::cout << NV_Ith_S( ydot, i ) << std::endl;
  { int dum; std::cin >> dum; }
#endif
  stats_sta.numRHS++;
  return( flag? 0: -1 );
}

template <typename... ExtOps>
inline
int
ODESLV_CVODES<ExtOps...>::MC_CVQUAD__
( sunrealtype t, N_Vector y, N_Vector qdot, void *user_data )
{
#ifdef MC__BASE_CVODES_CHECK
  //std::cout << "BASE_CVODES::PTR_BASE_CVODES: " << BASE_CVODES::PTR_BASE_CVODES
  //          << "  PTR_CVQUAD: " << PTR_CVQUAD << std::endl;
  assert( BASE_CVODES::PTR_BASE_CVODES != nullptr && PTR_CVQUAD != nullptr);
#endif
  return (BASE_CVODES::PTR_BASE_CVODES->*PTR_CVQUAD)( t, y, qdot, user_data );
}

template <typename... ExtOps>
inline
int
ODESLV_CVODES<ExtOps...>::CVQUAD__
( sunrealtype t, N_Vector y, N_Vector qdot, void *user_data )
{
  bool flag = _RHS_D_QUAD( t, NV_DATA_S( y ), NV_DATA_S( qdot ) );
  return( flag? 0: -1 );
}

template <typename... ExtOps>
inline
int
ODESLV_CVODES<ExtOps...>::MC_CVJAC__
( sunrealtype t, N_Vector y, N_Vector ydot, SUNMatrix Jac, void *user_data,
  N_Vector tmp1, N_Vector tmp2, N_Vector tmp3 )
{
#ifdef MC__BASE_CVODES_CHECK
  //std::cout << "BASE_CVODES::PTR_BASE_CVODES: " << BASE_CVODES::PTR_BASE_CVODES
  //          << "  PTR_CVJAC: " << PTR_CVJAC << std::endl;
  assert( BASE_CVODES::PTR_BASE_CVODES != nullptr && PTR_CVJAC != nullptr);
#endif
  return (BASE_CVODES::PTR_BASE_CVODES->*PTR_CVJAC)( t, y, ydot, Jac, user_data, tmp1, tmp2, tmp3 );
}

template <typename... ExtOps>
inline
int
ODESLV_CVODES<ExtOps...>::CVJAC__
( sunrealtype t, N_Vector y, N_Vector ydot, SUNMatrix Jac, void *user_data,
  N_Vector tmp1, N_Vector tmp2, N_Vector tmp3 )
{
  //std::cerr << "Entering: ODESLV_CVODES<ExtOps...>::CVJAC__\n";
  bool flag = false;
  //std::cerr << "Option: " << options.LINSOL << "\n";
  switch( options.LINSOL ){
   case Options::DIAG:
   case Options::DENSEDQ:
   case Options::SPARSE:
   default:
    flag = false;
    break;
   case Options::DENSE:
    flag = _JAC_D_STA( t, NV_DATA_S( y ), SM_COLS_D(Jac) );
    break;
#if 0
   case Options::SPARSE:
    flag = _JAC_D_STA( t, NV_DATA_S( y ), SM_DATA_S(Jac), SM_INDEXPTRS_S(Jac), SM_INDEXVALS_S(Jac) );
    break;
#endif
  }
  stats_sta.numJAC++; // increment JAC counter
  return( flag? 0: -1 );
}

template <typename... ExtOps>
typename ODESLV_CVODES<ExtOps...>::STATUS
ODESLV_CVODES<ExtOps...>::_states
( double const* p, double** xk, double* f, bool const store, std::ostream& os )
{
  // Check size
  //if( !p ) return STATUS::FATAL;

  try{
    // Initialize trajectory integration
    if( !_INI_STA( p ) ) return STATUS::FATAL;
    _t = _dT[0];
    const unsigned NSTEP = options.RESRECORD? options.RESRECORD: 1;

    // Initial state/quadrature values
    if( !_IC_D_SET()
     || !_IC_D_STA( _t, NV_DATA_S( _Nx ) )
     || (_Nq && !_IC_D_QUAD( NV_DATA_S( _Nq ) ) ) )
      { _END_STA(); return STATUS::FATAL; }
    _GET_D_STA( NV_DATA_S(_Nx), _nq && _Nq? NV_DATA_S(_Nq): 0 );

    // Store full state at initial time
    if( store ){
      sunrealtype*vsta = NV_DATA_S(_Nx);
      unsigned lsta = NV_LENGTH_S(_Nx);
      _vec_sta.push_back( std::vector<sunrealtype>( vsta, vsta+lsta ) );
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
    if( !_INI_CVODE() )
      { _END_STA(); return STATUS::FATAL; }

    for( _istg=0; _istg<_nsmax; _istg++ ){
      // State discontinuities (if any) at stage times
      // and integrator reinitialization (if applicable)
      _pos_ic = ( _vIC.size()>=_nsmax? _istg:0 );
      if( _pos_ic
       && ( !_CC_D_SET( _pos_ic )
         || !_CC_D_STA( _t, NV_DATA_S( _Nx ) )
         || !_CC_CVODE_STA() ) )
        { _END_STA(); return STATUS::FAILURE; }
      if( _istg 
       && !_CC_CVODE_QUAD() )
       //&& ( ( _Nq && !_IC_D_QUAD( NV_DATA_S( _Nq ) ) ) // quadrature reinitialization
       //  || !_CC_CVODE_QUAD() ) )
        { _END_STA(); return STATUS::FAILURE; }
      if( options.RESRECORD )
        results_sta.push_back( Results( _t, _nx, NV_DATA_S(_Nx), _nq, _nq? NV_DATA_S(_Nq): 0 ) );

      // update list of operations in RHS, JAC and QUAD
      _pos_rhs  = ( _vRHS.size()<=1? 0: _istg );
      _pos_quad = ( _vQUAD.size()<=1? 0: _istg );
      if( (!_istg || _pos_rhs || _pos_quad)
        && !_RHS_D_SET( _pos_rhs, _pos_quad ) )
        { _END_STA(); return STATUS::FATAL; }

      // integrate till end of time stage
      _cv_flag = CVodeSetStopTime( _cv_mem, _dT[_istg+1] );
      if( _check_cv_flag( &_cv_flag, "CVodeSetStopTime", 1 ) )
        { _END_STA(); return STATUS::FATAL; }
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
        if( _check_cv_flag( &_cv_flag, store?"CVodeF":"CVode", 1 ) )
         //|| (options.NMAX && stats_sta.numSteps > options.NMAX) )
          throw Exceptions( Exceptions::INTERN );

        // intermediate record
        if( options.RESRECORD ){
          if( _nq ){
            _cv_flag = CVodeGetQuad( _cv_mem, &_t, _Nq );
            if( _check_cv_flag(&_cv_flag, "CVodeGetQuad", 1) )
              { _END_STA(); return STATUS::FATAL; }
          }
          results_sta.push_back( Results( _t, _nx, NV_DATA_S(_Nx), _nq, _nq? NV_DATA_S(_Nq): 0 ) );
        }
      }

      // Store full state at stage time
      if( store ){
        sunrealtype*vsta = NV_DATA_S(_Nx);
        unsigned lsta = NV_LENGTH_S(_Nx);
        _vec_sta.push_back( std::vector<sunrealtype>( vsta, vsta+lsta ) );
      }

      // Intermediate states and quadratures
      if( _nq ){
        _cv_flag = CVodeGetQuad( _cv_mem, &_t, _Nq );
        if( _check_cv_flag(&_cv_flag, "CVodeGetQuad", 1) )
          { _END_STA(); return STATUS::FATAL; }
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
        { _END_STA(); return STATUS::FATAL; }
    }
#ifdef MC__ODESLV_CVODES_DEBUG
    if( store ) std::cout << "number of checkpoints: " << _nchk << std::endl;
#endif

    // Display / return function values
    for( unsigned i=0; f && i<_nf; i++ ) f[i] = _Df[i];
    if( options.DISPLAY >= 1 ){
      _print_interm( _nf, _Df.data(), " f", os );
      // Print final statistics to the screen
      //os << "\nFinal Statistics:\n";
      //_cv_flag = CVodePrintAllStats( _cv_mem, stdout, SUN_OUTPUTFORMAT_TABLE);
    }
  }
  catch(...){
    _END_STA();
    long int nstp;
    _cv_flag = CVodeGetNumSteps( _cv_mem, &nstp );
    stats_sta.numSteps += nstp;
    if( options.DISPLAY >= 1 ) _print_stats( stats_sta, os );
    return STATUS::FAILURE;
  }

  long int nstp;
  _cv_flag = CVodeGetNumSteps( _cv_mem, &nstp );
  stats_sta.numSteps += nstp;
#ifdef MC__ODESLV_CVODES_DEBUG
  std::cout << "number of steps: " << nstp << std::endl;
#endif

  _END_STA();
  if( options.DISPLAY >= 1 ) _print_stats( stats_sta, os );
  return STATUS::NORMAL;
}

//! @fn template <typename... ExtOps> inline typename ODESLV_CVODES<ExtOps...>::STATUS ODESLV_CVODES<ExtOps...>::states(
//! double const* p, double** xk, double* q, double* f, std::ostream& os=std::cout )
//!
//! This function computes a solution to the parametric ODEs:
//!   - <a>p</a>  [input]  parameter values
//!   - <a>xk</a> [output] state & quadrature values at stage times
//!   - <a>f</a>  [output] function values
//!   - <a>os</a> [input]  output stream [default: std::cout]
//! .
//! The return value is the status.
template <typename... ExtOps>
typename ODESLV_CVODES<ExtOps...>::STATUS
ODESLV_CVODES<ExtOps...>::states
( double const* p, double** xk, double* f, std::ostream& os )
{
  registration();
  STATUS flag = _states( p, xk, f, false, os );
  unregistration();
  return flag;
}

} // end namescape mc

#endif

