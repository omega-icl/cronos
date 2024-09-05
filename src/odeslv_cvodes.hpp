// Copyright (C) Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

#ifndef CRONOS__ODESLV_CVODES_HPP
#define CRONOS__ODESLV_CVODES_HPP

#undef  CRONOS__ODESLV_CVODES_DEBUG

#include "base_cvodes.hpp"
#include "odeslv_base.hpp"

namespace mc
{
//! @brief C++ class computing solutions of parametric ODEs using SUNDIALS/CVODES and MC++.
////////////////////////////////////////////////////////////////////////
//! mc::ODESLV_CVODES is a C++ class for solution of IVPs in ODEs
//! using the code CVODES in SUNDIALS and MC++.
////////////////////////////////////////////////////////////////////////
class ODESLV_CVODES
: public virtual BASE_CVODES
, public virtual ODESLV_BASE
, public virtual BASE_DE
{
 protected:
 
  using ODESLV_BASE::_print_interm;
  using ODESLV_BASE::_record;

  using ODESLV_BASE::_nsmax;
  using ODESLV_BASE::_istg;
  using ODESLV_BASE::_t;
  using ODESLV_BASE::_dT;
  using ODESLV_BASE::_nx;
  using ODESLV_BASE::_Dx;
  using ODESLV_BASE::_nq;
  using ODESLV_BASE::_Dq;
  using ODESLV_BASE::_nf;
  using ODESLV_BASE::_Df;
  using ODESLV_BASE::_vIC;
  using ODESLV_BASE::_vRHS;
  using ODESLV_BASE::_vQUAD;
  using ODESLV_BASE::_vFCT;
  using ODESLV_BASE::_nnzjac;

  using ODESLV_BASE::_INI_D_STA;
  using ODESLV_BASE::_END_D_STA;
  using ODESLV_BASE::_GET_D_STA;
  using ODESLV_BASE::_IC_D_SET;
  using ODESLV_BASE::_IC_D_STA;
  using ODESLV_BASE::_IC_D_QUAD;
  using ODESLV_BASE::_CC_D_SET;
  using ODESLV_BASE::_CC_D_STA;
  using ODESLV_BASE::_RHS_D_SET;
  using ODESLV_BASE::_RHS_D_STA;
  using ODESLV_BASE::_RHS_D_QUAD;
  using ODESLV_BASE::_JAC_D_STA;
  using ODESLV_BASE::_FCT_D_STA;

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

  //! @brief state values at stage times
  std::vector< std::vector< double > > _xk;

  //! @brief quadrature values at stage times
  std::vector< std::vector< double > > _qk;

  //! @brief function values
  std::vector< double > _f;

public:

  using BASE_DE::set;

  /** @ingroup ODESLV
   *  @{
   */
  typedef typename BASE_DE::STATUS STATUS;
  typedef typename ODESLV_BASE::Results Results;
  using BASE_DE::np;
  using BASE_DE::nf;

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
    //! @brief Whether or not to record results (default: 0)
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
        return "ODESLV_CVODES::Exceptions  Internal error";
       }
    }
   private:
    TYPE _ierr;
   };

  //! @brief Vector storing results (upon request only)
  std::vector< Results > results_state;

  //! @brief Statistics for state integration
  Stats stats_state;

  //! @brief Computes solution of parametric ODEs
  STATUS solve_state
    ( std::vector<double> const& p, std::vector<double> const& c=std::vector<double>() , std::ostream& os=std::cout );

  //! @brief Computes solution of parametric ODEs
  STATUS solve_state
    ( double const*p, double const*c=nullptr, std::ostream& os=std::cout );

  //! @brief Setup local copy of parametric ODEs
  bool setup
    ()
    { return ODESLV_BASE::_SETUP(); }

  //! @brief Setup local copy of parametric ODEs based on IVP
  bool setup
    ( ODESLV_CVODES const& IVP )
    { return ODESLV_BASE::_SETUP( IVP ); }

  //! @brief Record results in file <a>ores</a>, with accuracy of <a>iprec</a> digits
  void record
    ( std::ofstream& ores, unsigned const iprec=5 )
    const
    { return _record( ores, results_state, iprec ); }
    
  //! @brief Retreive state values at stage times
  std::vector< std::vector< double > > const& val_state
    ()
    const
    { return _xk; }

  //! @brief Retreive quadrature values at stage times
  std::vector< std::vector< double > > const& val_quadrature
    ()
    const
    { return _qk; }

  //! @brief Retreive function values at stage times
  std::vector< double > const& val_function
    ()
    const
    { return _f; }
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
    ( double const* p, double const* c );

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
  STATUS _states
    ( double const* p, double const* c, bool const store, std::ostream& os );

  //! @brief Solve parametric ODEs forward in time in current time stage
  STATUS _states_stage
    ( unsigned istg, double& t, N_Vector& Nx, N_Vector& Nq, bool const reinit, 
      bool const store, bool const record, std::ostream& os );

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
  ODESLV_CVODES( ODESLV_CVODES const& ) = delete;
  ODESLV_CVODES& operator=( ODESLV_CVODES const& ) = delete;
};

ODESLV_CVODES::ODESLV_CVODES
()
: BASE_CVODES(), BASE_DE(), ODESLV_BASE(),
  _cv_mem(nullptr), _cv_flag(0), _sun_mat(nullptr), _sun_ls(nullptr), _sun_nls(nullptr),
  _Nx(nullptr), _Nq(nullptr)
{}

ODESLV_CVODES::~ODESLV_CVODES
()
{
  if( _Nx ) N_VDestroy( _Nx );
  if( _Nq ) N_VDestroy( _Nq );
  if( _cv_mem )  CVodeFree( &_cv_mem );         /* Free CVODES memory */
  if( _sun_nls ) SUNNonlinSolFree( _sun_nls );  /* Free the nonlinear solver memory */
  if( _sun_ls )  SUNLinSolFree( _sun_ls );      /* Free the linear solver memory */
  if( _sun_mat ) SUNMatDestroy( _sun_mat );     /* Free the matrix memory */

}

bool
ODESLV_CVODES::_INI_CVODE
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
       _sun_mat = SUNDenseMatrix( _nx, _nx, sunctx );
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
       
#if defined( CRONOS__WITH_KLU )
     // Sparse Jacobian
     case Options::SPARSE:
       // Create sparse SUNMatrix for use in linear solves
       if( !this->set_sparse() ) return false;
       _sun_mat = SUNSparseMatrix( _nx, _nx, _nnzjac, CSC_MAT, sunctx );
       if( _check_cv_flag( (void*)_sun_mat, "SUNSparseMatrix", 0 ) ) return false;
       // Create sparse SUNLinearSolver object for use by CVode
       _sun_ls = SUNLinSol_KLU( _Nx, _sun_mat, sunctx );
       if( _check_cv_flag( (void *)_sun_ls, "SUNLinSol_Dense", 0 ) ) return false;
       // Attach the matrix and linear solver
       _cv_flag = CVodeSetLinearSolver( _cv_mem, _sun_ls, _sun_mat );
       if( _check_cv_flag( &_cv_flag, "CVodeSetLinearSolver", 1 ) ) return false;
       // Set the user-supplied Jacobian routine Jac
       _cv_flag = CVodeSetJacFn( _cv_mem, MC_CVJAC__ );
       if ( _check_cv_flag( &_cv_flag, "CVodeSetJacFn", 1 ) ) return false;
       break;
#endif
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

bool
ODESLV_CVODES::_CC_CVODE_STA
()
{
  // Reinitialize CVode memory block for current time _t and current state _Nx
  _cv_flag = CVodeReInit( _cv_mem, _t, _Nx );
  if( _check_cv_flag( &_cv_flag, "CVodeReInit", 1 ) ) return false;

#if defined( CRONOS__WITH_KLU )
  switch( options.LINSOL ){
    case Options::SPARSE:
      // Function SUNLinSol_KLUReInit(SUNLinearSolver S, SUNMatrix A, sunindextype nnz, int reinit_type)
      // needed to reinitialize memory and flag for a new factorization (symbolic and numeric) to be conducted at
      // the next solver setup call. This routine is useful in the cases where the number of nonzeroes has changed or if
      // the structure of the linear system has changed which would require a new symbolic (and numeric factorization).
      // std::cout << "Calling SUNLinSol_KLUReInit" << std::endl;
      _cv_flag = SUNLinSol_KLUReInit( _sun_ls, _sun_mat, _nnzjac, 2 );
      if( _cv_flag ) return false;
      break;
    default:
      break;
  }
#endif

  return true;
}

bool
ODESLV_CVODES::_CC_CVODE_QUAD
()
{

  // Reinitialize CVode memory block for current quarature _Nq
  if( !_Nq ) return true;
  _cv_flag = CVodeQuadReInit( _cv_mem, _Nq );
  if( _check_cv_flag( &_cv_flag, "CVodeQuadReInit", 1 ) ) return false;

  return true;
}

void
ODESLV_CVODES::_END_STA()
{
  // Unset constants
  _END_D_STA();

  // Get final CPU time
  _final_stats( stats_state );
}

bool
ODESLV_CVODES::_INI_STA
( double const* p, double const* c )
{
  // Initialize bound propagation
  if( !_INI_D_STA( p, c ) || !_INI_STA() )
    return false;
  return true;
}

bool
ODESLV_CVODES::_INI_STA
()
{
  // Set SUNDIALS state/quadrature arrays
  if( !_Nx || NV_LENGTH_S( _Nx ) != (sunindextype)_nx ){
    if( _Nx ) N_VDestroy( _Nx );
    _Nx  = N_VNew_Serial( _nx, sunctx );
  }
  if( !_Nq || NV_LENGTH_S( _Nq ) != (sunindextype)_nq ){
    if( _Nq ) N_VDestroy( _Nq );
    _Nq  = _nq? N_VNew_Serial( _nq, sunctx ): nullptr;
  }

  // reset at time stages
  _xk.clear(); _xk.reserve(_nsmax);
  _qk.clear(); _qk.reserve(_nsmax);
  _f.clear();  _f.reserve(_nf);

  // Reset result record and statistics
  results_state.clear();
  _init_stats( stats_state );

  return true;
}

inline
int
ODESLV_CVODES::MC_CVRHS__
( sunrealtype t, N_Vector y, N_Vector ydot, void *user_data )
{
#ifdef CRONOS__BASE_CVODES_CHECK
  //std::cout << "BASE_CVODES::PTR_BASE_CVODES: " << BASE_CVODES::PTR_BASE_CVODES
  //          << "  PTR_CVRHS: " << PTR_CVRHS << std::endl;
  assert( BASE_CVODES::PTR_BASE_CVODES != nullptr && PTR_CVRHS != nullptr);
#endif
  return (BASE_CVODES::PTR_BASE_CVODES->*PTR_CVRHS)( t, y, ydot, user_data );
}

inline
int
ODESLV_CVODES::CVRHS__
( sunrealtype t, N_Vector y, N_Vector ydot, void *user_data )
{
#ifdef CRONOS__ODESLV_CVODES_DEBUG
  std::cout << std::scientific << std::setprecision(6) << t;
  for( unsigned i=0; i<NV_LENGTH_S( y ); i++ ) std::cout << "  " << NV_Ith_S( y, i );
#endif
  bool flag = _RHS_D_STA( t, NV_DATA_S( y ), NV_DATA_S( ydot ) );
#ifdef CRONOS__ODESLV_CVODES_DEBUG
  for( unsigned i=0; i<NV_LENGTH_S( ydot ); i++ ) std::cout << "  " << NV_Ith_S( ydot, i );
  std::cout << std::endl;
  { int dum; std::cin >> dum; }
#endif
  stats_state.numRHS++;
  return( flag? 0: -1 );
}

inline
int
ODESLV_CVODES::MC_CVQUAD__
( sunrealtype t, N_Vector y, N_Vector qdot, void *user_data )
{
#ifdef CRONOS__BASE_CVODES_CHECK
  //std::cout << "BASE_CVODES::PTR_BASE_CVODES: " << BASE_CVODES::PTR_BASE_CVODES
  //          << "  PTR_CVQUAD: " << PTR_CVQUAD << std::endl;
  assert( BASE_CVODES::PTR_BASE_CVODES != nullptr && PTR_CVQUAD != nullptr);
#endif
  return (BASE_CVODES::PTR_BASE_CVODES->*PTR_CVQUAD)( t, y, qdot, user_data );
}

inline
int
ODESLV_CVODES::CVQUAD__
( sunrealtype t, N_Vector y, N_Vector qdot, void *user_data )
{
  bool flag = _RHS_D_QUAD( t, NV_DATA_S( y ), NV_DATA_S( qdot ) );
  return( flag? 0: -1 );
}

inline
int
ODESLV_CVODES::MC_CVJAC__
( sunrealtype t, N_Vector y, N_Vector ydot, SUNMatrix Jac, void *user_data,
  N_Vector tmp1, N_Vector tmp2, N_Vector tmp3 )
{
#ifdef CRONOS__BASE_CVODES_CHECK
  //std::cout << "BASE_CVODES::PTR_BASE_CVODES: " << BASE_CVODES::PTR_BASE_CVODES
  //          << "  PTR_CVJAC: " << PTR_CVJAC << std::endl;
  assert( BASE_CVODES::PTR_BASE_CVODES != nullptr && PTR_CVJAC != nullptr);
#endif
  return (BASE_CVODES::PTR_BASE_CVODES->*PTR_CVJAC)( t, y, ydot, Jac, user_data, tmp1, tmp2, tmp3 );
}

inline
int
ODESLV_CVODES::CVJAC__
( sunrealtype t, N_Vector y, N_Vector ydot, SUNMatrix Jac, void *user_data,
  N_Vector tmp1, N_Vector tmp2, N_Vector tmp3 )
{
  //std::cerr << "Entering: ODESLV_CVODES::CVJAC__\n";
  bool flag = false;
  //std::cerr << "Option: " << options.LINSOL << "\n";
  switch( options.LINSOL ){
   case Options::DIAG:
   case Options::DENSEDQ:
   default:
    flag = false;
    break;

   case Options::DENSE:
    flag = _JAC_D_STA( t, NV_DATA_S( y ), SM_COLS_D(Jac) );
    break;

#if defined( CRONOS__WITH_KLU )
   case Options::SPARSE:
    flag = _JAC_D_STA( t, NV_DATA_S( y ), SUNSparseMatrix_Data(Jac),
                       SUNSparseMatrix_IndexPointers(Jac), SUNSparseMatrix_IndexValues(Jac) );
    break;
#endif
  }
  stats_state.numJAC++; // increment JAC counter
  return( flag? 0: -1 );
}

typename ODESLV_CVODES::STATUS
ODESLV_CVODES::_states_stage
( unsigned istg, double& t, N_Vector& Nx, N_Vector& Nq, bool const reinit, bool const store,
  bool const record, std::ostream& os )
{
  // State discontinuities (if any) at stage times
  // and integrator reinitialization (if applicable)
  _pos_ic = ( _vIC.size()>=_nsmax? istg:0 );
  if( _pos_ic && ( !_CC_D_SET( _pos_ic )
                || !_CC_D_STA( t, NV_DATA_S( Nx ) )
                || !_CC_CVODE_STA() ) )
    { _END_STA(); return STATUS::FAILURE; }
  else if( !istg && reinit && !_CC_CVODE_STA() )
    { _END_STA(); return STATUS::FAILURE; }
  //if( istg && !_CC_CVODE_QUAD() )
  if( ( istg || reinit )
   && ( ( Nq && !_IC_D_QUAD( NV_DATA_S( Nq ) ) ) // quadrature reinitialization
     || !_CC_CVODE_QUAD() ) )
    { _END_STA(); return STATUS::FAILURE; }
  if( record )
    results_state.push_back( Results( t, _nx, NV_DATA_S(Nx), _nq, _nq? NV_DATA_S(Nq): nullptr ) );

  // update list of operations in RHS, JAC and QUAD
  _pos_rhs  = ( _vRHS.size()<=1? 0: istg );
  _pos_quad = ( _vQUAD.size()<=1? 0: istg );
  if( (!istg || _pos_rhs || _pos_quad)
    && !_RHS_D_SET( _pos_rhs, _pos_quad ) )
    { _END_STA(); return STATUS::FATAL; }

  // integrate till end of time stage
  _cv_flag = CVodeSetStopTime( _cv_mem, _dT[istg+1] );
  if( _check_cv_flag( &_cv_flag, "CVodeSetStopTime", 1 ) )
    { _END_STA(); return STATUS::FATAL; }

  unsigned const NSTEP = options.RESRECORD? options.RESRECORD: 1;
  double const TSTEP = ( _dT[istg+1] - t ) / NSTEP;
  double TSTOP = t + TSTEP;
  for( unsigned k=0; k<NSTEP; k++, TSTOP+=TSTEP ){

    if( k+1 == NSTEP ) TSTOP = _dT[istg+1];
    if( !store )
      _cv_flag = CVode( _cv_mem, TSTOP, Nx, &t, CV_NORMAL );
    else
      _cv_flag = CVodeF( _cv_mem, TSTOP, Nx, &t, CV_NORMAL, &_nchk );
    if( _check_cv_flag( &_cv_flag, store?"CVodeF":"CVode", 1 ) )
     //|| (options.NMAX && stats_state.numSteps > options.NMAX) )
      throw Exceptions( Exceptions::INTERN );

    // intermediate record
    if( record ){
      if( _nq ){
        _cv_flag = CVodeGetQuad( _cv_mem, &t, Nq );
        if( _check_cv_flag(&_cv_flag, "CVodeGetQuad", 1) )
          { _END_STA(); return STATUS::FATAL; }
      }
      results_state.push_back( Results( t, _nx, NV_DATA_S(Nx), _nq, _nq? NV_DATA_S(Nq): nullptr ) );
    }
  }

  return STATUS::NORMAL;
}

typename ODESLV_CVODES::STATUS
ODESLV_CVODES::_states
( double const* p, double const* c, bool const store, std::ostream& os )
{
  try{
    // Initialize trajectory integration
    if( !_INI_STA( p, c ) ) return STATUS::FATAL;
    _t = _dT[0];

    // Initial state/quadrature values
    if( !_IC_D_SET()
     || !_IC_D_STA( _t, NV_DATA_S( _Nx ) )
     || (_Nq && !_IC_D_QUAD( NV_DATA_S( _Nq ) ) ) )
      { _END_STA(); return STATUS::FATAL; }
    _GET_D_STA( NV_DATA_S(_Nx), _nq && _Nq? NV_DATA_S(_Nq): nullptr );

//    // Store full state at initial time
//    if( store ){
//      sunrealtype*vsta = NV_DATA_S(_Nx);
//      unsigned lsta = NV_LENGTH_S(_Nx);
//      _vec_sta.push_back( std::vector<double>( vsta, vsta+lsta ) );
//    }

    // Display / record / return initial results
    _xk.push_back( std::vector<double>( _Dx, _Dx+_nx ) );
    if( _nq ) _qk.push_back( std::vector<double>( _Dq, _Dq+_nq ) );
    if( options.DISPLAY >= 1 ){
      _print_interm( _t, _nx, _Dx, " x", os );
      _print_interm( _nq, _Dq, " q", os );
    }
//    if( options.RESRECORD )
//      results_state.push_back( Results( _t, _nx, NV_DATA_S(_Nx), _nq, _nq? NV_DATA_S(_Nq): nullptr ) );

    // Integrate ODEs through each stage using SUNDIALS
    if( !_INI_CVODE() )
      { _END_STA(); return STATUS::FATAL; }

    for( _istg=0; _istg<_nsmax; _istg++ ){

      // Integrate states over stage
      _states_stage( _istg, _t, _Nx, _Nq, false, store, options.RESRECORD, os );

//      // Store full state at stage time
//      if( store ){
//        sunrealtype*vsta = NV_DATA_S(_Nx);
//        unsigned lsta = NV_LENGTH_S(_Nx);
//        _vec_sta.push_back( std::vector<double>( vsta, vsta+lsta ) );
//      }

      // Intermediate states and quadratures
      if( _nq ){
        _cv_flag = CVodeGetQuad( _cv_mem, &_t, _Nq );
        if( _check_cv_flag(&_cv_flag, "CVodeGetQuad", 1) )
          { _END_STA(); return STATUS::FATAL; }
      }
     _GET_D_STA( NV_DATA_S(_Nx), _nq && _Nq? NV_DATA_S(_Nq): 0 );

      // Display / return stage results
      _xk.push_back( std::vector<double>( _Dx, _Dx+_nx ) );
      if( _nq ) _qk.push_back( std::vector<double>( _Dq, _Dq+_nq ) );
      if( options.DISPLAY >= 1 ){
        _print_interm( _t, _nx, _Dx, " x", os );
        _print_interm( _nq, _Dq, " q", os );
      }

      // Add intermediate function terms
      _pos_fct = ( _vFCT.size()>=_nsmax? _istg:0 );
      if( (_vFCT.size()>=_nsmax || _istg==_nsmax-1) && !_FCT_D_STA( _pos_fct, _t ) )
        { _END_STA(); return STATUS::FATAL; }
    }
#ifdef CRONOS__ODESLV_CVODES_DEBUG
    if( store ) std::cout << "number of checkpoints: " << _nchk << std::endl;
#endif

    // Display / return function values
    _f = _Df;
    if( options.DISPLAY >= 1 ){
      _print_interm( _nf, _Df.data(), " f", os );
      // Print final statistics to the screen
      // os << "\nFinal Statistics:\n";
      // _cv_flag = CVodePrintAllStats( _cv_mem, stdout, SUN_OUTPUTFORMAT_TABLE);
    }
  }
  catch(...){
    _END_STA();
    long int nstp;
    _cv_flag = CVodeGetNumSteps( _cv_mem, &nstp );
    stats_state.numSteps += nstp;
    if( options.DISPLAY >= 1 ) _print_stats( stats_state, os );
    //std::cout << "failed status: " << STATUS::FAILURE << std::endl;
    return STATUS::FAILURE;
  }

  long int nstp;
  _cv_flag = CVodeGetNumSteps( _cv_mem, &nstp );
  stats_state.numSteps += nstp;
#ifdef CRONOS__ODESLV_CVODES_DEBUG
  std::cout << "number of steps: " << nstp << std::endl;
#endif

  _END_STA();
  if( options.DISPLAY >= 1 ) _print_stats( stats_state, os );
  //std::cout << "normal status: " << STATUS::NORMAL << std::endl;
  return STATUS::NORMAL;
}

//! @fn inline typename ODESLV_CVODES::STATUS ODESLV_CVODES::solve_state(
//! std::vector<double> const& p, std::vector<double> const& c=std::vector<double>(), std::ostream& os=std::cout )
//!
//! This function computes a solution to the parametric ODEs:
//!   - <a>p</a>  [input]  parameter values
//!   - <a>c</a>  [input]  constant values
//!   - <a>os</a> [input/ouptut]  output stream [default: std::cout]
//! .
//! The return value is the status.
typename ODESLV_CVODES::STATUS
ODESLV_CVODES::solve_state
( std::vector<double> const& p, std::vector<double> const& c, std::ostream& os )
{
  registration();
  STATUS flag = _states( p.data(), c.data(), false, os );
  unregistration();
  return flag;
}

//! @fn inline typename ODESLV_CVODES::STATUS ODESLV_CVODES::solve_state(
//! double const* p, double const* c=nullptr, std::ostream& os=std::cout )
//!
//! This function computes a solution to the parametric ODEs:
//!   - <a>p</a>  [input]  parameter values
//!   - <a>c</a>  [input]  constant values
//!   - <a>os</a> [input/output]  output stream [default: std::cout]
//! .
//! The return value is the status.
typename ODESLV_CVODES::STATUS
ODESLV_CVODES::solve_state
( double const* p, double const* c, std::ostream& os )
{
  registration();
  STATUS flag = _states( p, c, false, os );
  unregistration();
  return flag;
}

} // end namescape mc

#endif

