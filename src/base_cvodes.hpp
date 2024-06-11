// Copyright (C) Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

#ifndef MC__BASE_CVODES_HPP
#define MC__BASE_CVODES_HPP

#include <iostream>
#include <iomanip>
#include <cassert>
#include <thread>
#include <sys/time.h>

#include <cvodes/cvodes.h>              /* prototypes for CVODE fcts., consts.  */
#include <nvector/nvector_serial.h>     /* access to serial N_Vector            */
//#include <sundials/sundials_math.h>     /* definition of SUNRabs, SUNRexp, etc. */
#include <sunmatrix/sunmatrix_dense.h>  /* access to dense SUNMatrix            */
//#include <sunmatrix/sunmatrix_sparse.h> /* access to sparse SUNMatrix           */
#include <sunlinsol/sunlinsol_dense.h>  /* access to dense SUNLinearSolver      */
//#include <sunlinsol/sunlinsol_klu.h>    /* access to sparse SUNLinearSolver     */
#include <cvodes/cvodes_diag.h>         /* access to CVDIAG linear solver       */
#include "sunnonlinsol/sunnonlinsol_newton.h"     /* access to the newton SUNNonlinearSolver      */
#include "sunnonlinsol/sunnonlinsol_fixedpoint.h" /* access to the fixed point SUNNonlinearSolver */

#include "base_sundials.hpp"

namespace mc
{
//! @brief C++ base class for computing solutions of parametric parametric differential-algebraic equations using SUNDIALS
////////////////////////////////////////////////////////////////////////
//! mc::BASE_CVODES is a C++ base class for computing solutions of
//! parametric differential-algebraic equations (DAEs) using SUNDIALS
////////////////////////////////////////////////////////////////////////
class BASE_CVODES
: public virtual BASE_SUNDIALS
{
public:
  /** @ingroup ODESLV_CVODES
   *  @ingroup ODEBND_CVODES
   *  @{
   */
  //! @brief Default class constructor
  BASE_CVODES()
    : _is_registered( false )
    {}

  //! @brief Class destructor
  virtual ~BASE_CVODES()
    {}
    
  //! @brief Function registering thread
  void registration
    ();

  //! @brief Function unregistering thread
  void unregistration
    ();

  //! @brief GSL options
  struct Options
  {
    //! @brief Constructor
    Options():
      INTMETH(MSBDF), NLINSOL(NEWTON), LINSOL(DIAG), 
      H0(0e0), HMIN(1e-12), HMAX(0e0), NMAX(2000),
      RTOL(1e-8), ATOL(1e-8), ETOL(1e-16),
      QERR(true), MAXFAIL(10), MAXCORR(5),
      AUTOTOLS(false), RTOLS(1e-8), ATOLS(1e-8), ETOLS(1e-16),
      FSACORR(STAGGERED), FSAERR(true), QERRS(true),
      RTOLB(1e-8), ATOLB(1e-8), ETOLB(1e-16), QERRB(true),
      ASAINTERP(HERMITE), ASACHKPT(2000)//, RTOLFD(1e-3), ATOLFD(1e-3), CENFD(true)
      {}
    //! @brief Assignment operator
    template <typename U> Options& operator=
      ( U&options ){
        INTMETH   = options.INTMETH;
        NLINSOL   = options.NLINSOL;
        LINSOL    = options.LINSOL;
        H0        = options.H0;
        HMIN      = options.HMIN;
        HMAX      = options.HMAX;
        NMAX      = options.NMAX;
        RTOL      = options.RTOL;
        ATOL      = options.ATOL;
        ETOL      = options.ETOL;
        QERR      = options.QERR;
        MAXFAIL   = options.MAXFAIL;
        MAXCORR   = options.MAXCORR;
        AUTOTOLS  = options.AUTOTOLS;
        RTOLS     = options.RTOLS;
        ATOLS     = options.ATOLS;
        ETOLS     = options.ETOLS;
        FSACORR   = options.FSACORR;
        FSAERR    = options.FSAERR;
        QERRS     = options.QERRS;
        RTOLB     = options.RTOLB;
        ATOLB     = options.ATOLB;
        ETOLB     = options.ETOLB;
        QERRB     = options.QERRB;
        ASAINTERP = options.ASAINTERP;
        ASACHKPT  = options.ASACHKPT;
        //RTOLFD    = options.RTOLFD;
        //ATOLFD    = options.ATOLFD;
        //CENFD     = options.CENFD;
        return *this;
      }
    //! @brief Enumeration type for FSA method
    enum FSA_STRATEGY{
      SIMULTANEOUS=CV_SIMULTANEOUS, //!< Simultaneous state/sensitivity correction
      STAGGERED=CV_STAGGERED,       //!< Simultaneous sensitivity corrections after state corrections
      STAGGERED1=CV_STAGGERED1      //!< Sequential sensitivity corrections after state corrections
    };
    //! @brief Enumeration type for ASA method
    enum ASA_STRATEGY{
      HERMITE=CV_HERMITE,	//!< Cubic Hermite interpolation
      POLYNOMIAL=CV_POLYNOMIAL	//!< Variable degree polynomial interpolation
    };
    //! @brief Enumeration of numerical integration algorithms
    enum INTEGRATION_METHOD{
      MSADAMS=CV_ADAMS,	//!< Variable-coefficient linear multistep Adams method (non-stiff systems)
      MSBDF=CV_BDF	//!< Variable-coefficient linear multistep backward differentiation formula (BDF) method (stiff systems)
    };
    //! @brief Enumeration of nonlinear solver strategies
    enum NONLINEAR_SOLVER{
      FIXEDPOINT=0,	//!< Fixed point nonlinear solver
      NEWTON		//!< Newton nonlinear solver
    };
    //! @brief Enumeration of linear solver strategies (within Newton nonlinear solver)
    enum LINEAR_SOLVER{
      DIAG=0,	//!< Approximate diagonal Jacobian formed by way of a difference quotient
      DENSE,	//!< Use analytic dense Jacobian and internal direct dense linear algebra functions
      DENSEDQ,	//!< Use approximate dense Jacobian by way of a difference quotient and internal direct dense linear algebra functions
      SPARSE	//!< Use analytic sparse Jacobian and use of internal direct dense linear algebra functions
    };
    //! @brief Numerical integration method [Default: MSADAMS]
    INTEGRATION_METHOD INTMETH;
    //! @brief Nonlinear solver method [Default: FIXEDPOINT]
    NONLINEAR_SOLVER NLINSOL;
    //! @brief Linear solver method and Jacobian approximation [Default: DIAG]
    LINEAR_SOLVER LINSOL;
    //! @brief Initial step-size [Default: 0e0 (auto)]
    double H0;
    //! @brief Minimum step-size [Default: 0e0]
    double HMIN;
    //! @brief Maximum step-size [Default: 0e0 (+inf)]
    double HMAX;
    //! @brief Maximum number of steps in a time stage [Default: 0 (500)]
    unsigned int NMAX;
    //! @brief Relative (scalar) integration tolerance
    sunrealtype RTOL;
    //! @brief Absolute (scalar) integration tolerance
    sunrealtype ATOL;
    //! @brief Absolute (scalar) integration tolerance for ellipsoid shape matrix
    sunrealtype ETOL;    
    //! @brief Whether or not (state) quadrature error control is performed?
    bool QERR;
    //! @brief Maximum number of error test failures per step
    int MAXFAIL;
    //! @brief Maximum  number of nonlinear solver iterations per step
    int MAXCORR;
    //! @brief Whether integration tolerances for FS are to be set automatically?
    bool AUTOTOLS;
    //! @brief Relative (scalar) integration tolerance for FSA
    sunrealtype RTOLS;
    //! @brief Absolute (scalar) integration tolerance for FSA
    sunrealtype ATOLS;
    //! @brief Absolute (scalar) integration tolerance for ellipsoid shape matrix in FSA
    sunrealtype ETOLS;      
    //! @brief FSA correction strategy
    FSA_STRATEGY FSACORR;
    //! @brief Whether or not FSA error control is performed?
    bool FSAERR;
    //! @brief Whether or not FSA quadrature error control is performed?
    bool QERRS;
    //! @brief Relative (scalar) integration tolerance for ASA
    sunrealtype RTOLB;
    //! @brief Absolute (scalar) integration tolerance for ASA
    sunrealtype ATOLB;
    //! @brief Absolute (scalar) integration tolerance for ellipsoid shape matrix in ASA
    sunrealtype ETOLB;    
    //! @brief Whether or not ASA quadrature error control is performed?
    bool QERRB;
    //! @brief ASA interpolation strategy
    ASA_STRATEGY ASAINTERP;
    //! @brief Number of steps between each check point for ASA
    int ASACHKPT;
    //! @brief Relative tolerance for finite differences
    double RTOLFD;
    //! @brief Absolute tolerance for finite differences
    double ATOLFD;
    //! @brief Whether or not centered finite differences are used?
    bool CENFD;
  } options;

  //! @brief Structure storing integration statistics
  struct Stats
  {
    //! @brief Constructor
    Stats():
      cputime(0.), numSteps(0), numRHS(0), numJAC(0)
      {}
    //! @brief Constructor
    void reset()
      { cputime = 0.; numSteps = numRHS = numJAC = 0; }

    //! @brief CPU time
    double cputime;
    //! @brief Number of integration steps
    unsigned int numSteps;
    //! @brief Number of right-hand side (RHS) evaluations
    unsigned int numRHS;
    //! @brief Number of Jacobian evaluations
    unsigned int numJAC;
  };

  //! @brief Thread-safe static pointer to self
  thread_local static BASE_CVODES* PTR_BASE_CVODES;
  /** @} */

protected:
  //! @brief Flag indicating whether the thread is registered
  bool _is_registered;

  //! @brief Pure virtual function to calculate the ODEs RHS derivatives
  virtual int CVRHS__
    ( sunrealtype t, N_Vector Nx, N_Vector Nxdot, void* user_data )
    = 0;

  //! @brief Pure virtual function to calculate the quadrature RHS derivatives
  virtual int CVQUAD__
    ( sunrealtype t, N_Vector Nx, N_Vector Nqdot, void* user_data )
    = 0;

  //! @brief Pure virtual function to calculate the ODEs RHS Jacobian
  virtual int CVJAC__
    ( sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix Jac, void *user_data,
      N_Vector tmp1, N_Vector tmp2, N_Vector tmp3 )
    = 0;

  //! @brief Pure virtual function to calculate the adjoint ODEs RHS derivatives
  virtual int CVRHSB__
    ( sunrealtype t, N_Vector x, N_Vector y, N_Vector ydot, void* user_data )
    = 0;

  //! @brief Pure virtual function to calculate the adjoint quadrature RHS derivatives
  virtual int CVQUADB__
    ( sunrealtype t, N_Vector x, N_Vector y, N_Vector qdot, void* user_data )
    = 0;

  //! @brief Pure virtual function to calculate the adjoint ODEs RHS Jacobian
  virtual int CVJACB__
    ( sunrealtype t, N_Vector x, N_Vector y, N_Vector fB, SUNMatrix JacB,
      void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3 )
    = 0;

  //! @brief Pure virtual function to calculate the sensitivity ODEs RHS derivatives
  virtual int CVRHSF__
    ( int Ns, sunrealtype t, N_Vector x, N_Vector xdot, int is, N_Vector y,
      N_Vector ydot, void* user_data, N_Vector tmp1, N_Vector tmp2 )
    = 0;

  //! @brief Pure virtual function to calculate the sensitivity quadrature RHS derivatives
  virtual int CVQUADF__
    ( int Ns, sunrealtype t, N_Vector x, N_Vector* y, N_Vector qdot, N_Vector* qSdot, 
      void *user_data, N_Vector tmp1, N_Vector tmp2 )
    = 0;

  //! @brief Function to initialize CVODES statistics
  static void _init_stats
    ( Stats& stats );

  //! @brief Function to finalize CVODES statistics
  static void _final_stats
    ( Stats& stats );

  //! @brief Function to display integration statistics
  static void _print_stats
    ( Stats const& stats, std::ostream& os=std::cout );

  //! @brief Function to check function return values with CVODE
  static bool _check_cv_flag
    ( void *flagvalue, std::string const& funcname, int const opt );

  //! @brief Function to display CVODE statistics
  static void _print_stats_cvode
    ( void* cvode_mem, bool const sens=false, std::ostream& os=std::cout );

  //! @brief Private methods to block default compiler methods
  BASE_CVODES( BASE_CVODES const& ) = delete;
  BASE_CVODES& operator=( BASE_CVODES const& ) = delete;
};

thread_local BASE_CVODES* BASE_CVODES::PTR_BASE_CVODES = nullptr;
thread_local int (BASE_CVODES::*PTR_CVRHS)  ( sunrealtype, N_Vector, N_Vector, void* ) = nullptr;
thread_local int (BASE_CVODES::*PTR_CVQUAD) ( sunrealtype, N_Vector, N_Vector, void* ) = nullptr;
thread_local int (BASE_CVODES::*PTR_CVJAC)  ( sunrealtype, N_Vector, N_Vector, SUNMatrix, void*, N_Vector, N_Vector, N_Vector ) = nullptr;
thread_local int (BASE_CVODES::*PTR_CVRHSB) ( sunrealtype, N_Vector, N_Vector, N_Vector, void* ) = nullptr;
thread_local int (BASE_CVODES::*PTR_CVQUADB)( sunrealtype, N_Vector, N_Vector, N_Vector, void* ) = nullptr;
thread_local int (BASE_CVODES::*PTR_CVJACB) ( sunrealtype, N_Vector, N_Vector, N_Vector, SUNMatrix, void*, N_Vector, N_Vector, N_Vector ) = nullptr;
thread_local int (BASE_CVODES::*PTR_CVRHSF) ( int, sunrealtype, N_Vector, N_Vector, int, N_Vector, N_Vector, void*, N_Vector, N_Vector ) = nullptr;
thread_local int (BASE_CVODES::*PTR_CVQUADF)( int, sunrealtype, N_Vector, N_Vector*, N_Vector, N_Vector*, void*, N_Vector, N_Vector ) = nullptr;

void REG_BASE_CVODES
( BASE_CVODES* CV,
  int (BASE_CVODES::*CVRhsFn)  ( sunrealtype, N_Vector, N_Vector, void* ),
  int (BASE_CVODES::*CVQuadFn) ( sunrealtype, N_Vector, N_Vector, void* ),
  int (BASE_CVODES::*CVJacFn)  ( sunrealtype, N_Vector, N_Vector, SUNMatrix, void*, N_Vector, N_Vector, N_Vector ),
  int (BASE_CVODES::*CVRhsFnB) ( sunrealtype, N_Vector, N_Vector, N_Vector, void* ),
  int (BASE_CVODES::*CVQuadFnB)( sunrealtype, N_Vector, N_Vector, N_Vector, void* ),
  int (BASE_CVODES::*CVJacFnB) ( sunrealtype, N_Vector, N_Vector, N_Vector, SUNMatrix, void*, N_Vector, N_Vector, N_Vector ),
  int (BASE_CVODES::*CVRhsFnF) ( int, sunrealtype, N_Vector, N_Vector, int, N_Vector, N_Vector, void*, N_Vector, N_Vector ),
  int (BASE_CVODES::*CVQuadFnF)( int, sunrealtype, N_Vector, N_Vector*, N_Vector, N_Vector*, void*, N_Vector, N_Vector ) )
{
#ifdef MC__BASE_CVODES_CHECK
  assert( BASE_CVODES::PTR_BASE_CVODES == nullptr && PTR_CVRHS == nullptr );
#endif
  BASE_CVODES::PTR_BASE_CVODES = CV;
  PTR_CVRHS   = CVRhsFn;
  PTR_CVQUAD  = CVQuadFn;
  PTR_CVJAC   = CVJacFn;
  PTR_CVRHSB  = CVRhsFnB;
  PTR_CVQUADB = CVQuadFnB;
  PTR_CVJACB  = CVJacFnB;
  PTR_CVRHSF  = CVRhsFnF;
  PTR_CVQUADF = CVQuadFnF;
}

inline
void
BASE_CVODES::registration
()
{
  if( _is_registered ) return;
  REG_BASE_CVODES( this,
                   &BASE_CVODES::CVRHS__,
                   &BASE_CVODES::CVQUAD__,
                   &BASE_CVODES::CVJAC__,
                   &BASE_CVODES::CVRHSB__,
                   &BASE_CVODES::CVQUADB__,
                   &BASE_CVODES::CVJACB__,
                   &BASE_CVODES::CVRHSF__,
                   &BASE_CVODES::CVQUADF__
                 );
  _is_registered = true;
}

void UNREG_BASE_CVODES
()
{
#ifdef MC__BASE_CVODES_CHECK
  assert( BASE_CVODES::PTR_BASE_CVODES != nullptr
       && PTR_CVRHS   != nullptr
       && PTR_CVQUAD  != nullptr
       && PTR_CVJAC   != nullptr
       && PTR_CVRHSB  != nullptr
       && PTR_CVQUADB != nullptr
       && PTR_CVJACB  != nullptr
       && PTR_CVRHSF  != nullptr
       && PTR_CVQUADF != nullptr );
#endif
  BASE_CVODES::PTR_BASE_CVODES = nullptr;
  PTR_CVRHS   = nullptr;
  PTR_CVQUAD  = nullptr;
  PTR_CVJAC   = nullptr;
  PTR_CVRHSB  = nullptr;
  PTR_CVQUADB = nullptr;
  PTR_CVJACB  = nullptr;
  PTR_CVRHSF  = nullptr;
  PTR_CVQUADF = nullptr;
}

inline
void
BASE_CVODES::unregistration
()
{
  if( !_is_registered ) return;
  UNREG_BASE_CVODES();
  _is_registered = false;
}

inline
void
BASE_CVODES::_init_stats
( Stats& stats )
{
  // Initialize statistics
  stats.reset();
  timeval time;
  gettimeofday(&time, 0) ;
  stats.cputime = - time.tv_sec - time.tv_usec*1e-6;
}

inline
void
BASE_CVODES::_final_stats
( Stats& stats )
{
  // Get final CPU time
  timeval time;
  gettimeofday(&time, 0);
  stats.cputime += time.tv_sec + time.tv_usec*1e-6;
}

inline
void
BASE_CVODES::_print_stats
( Stats const& stats, std::ostream& os )
{
  // Statistics
  os << " No STEPS    " << stats.numSteps
     << std::endl
     << " No EVALATIONS" << "   RHS: " << stats.numRHS
                         << "   JAC: " << stats.numJAC
     << std::endl
     << " CPU TIME (SEC)     " << std::fixed << std::left
                               << std::setprecision(5) << stats.cputime
     << std::endl << std::endl;
  return;
}

inline
bool
BASE_CVODES::_check_cv_flag
( void* flagvalue, std::string const& funcname, int const opt )
{
  // Check if SUNDIALS function returned NULL pointer - no memory allocated
  if( opt == 0 && flagvalue == nullptr ){
    std::cerr << "\nSUNDIALS_ERROR: " << funcname
              << "() failed - returned NULL pointer\n\n";
    return true;
  }

  // Check if flag < 0
  else if( opt == 1 ){
    int *errflag = (int *) flagvalue;
    if( *errflag < 0 ){
      std::cerr << "\nSUNDIALS_ERROR: " << funcname
                << "() failed with flag = " << *errflag << "\n\n";
      return true;
    }
  }

  // Check if function returned NULL pointer - no memory allocated
  else if( opt == 2 && flagvalue == nullptr ){
    std::cerr << "\nMEMORY_ERROR: " << funcname
              << "() failed - returned NULL pointer\n\n";
    return true;
  }

  return false;
}

inline
void
BASE_CVODES::_print_stats_cvode
( void* cvode_mem, bool const sens, std::ostream& os )
{
  long int nst, nfe, nsetups, nje, nfeLS, nni, ncfn, netf;

  int flag = CVodeGetNumSteps(cvode_mem, &nst);
  if( _check_cv_flag(&flag, "CVodeGetNumSteps", 1) ) return;
  flag = CVodeGetNumRhsEvals(cvode_mem, &nfe);
  if( _check_cv_flag(&flag, "CVodeGetNumRhsEvals", 1) ) return;
  flag = CVodeGetNumLinSolvSetups(cvode_mem, &nsetups);
  if( _check_cv_flag(&flag, "CVodeGetNumLinSolvSetups", 1) ) return;
  flag = CVodeGetNumErrTestFails(cvode_mem, &netf);
  if( _check_cv_flag(&flag, "CVodeGetNumErrTestFails", 1) ) return;
  flag = CVodeGetNumNonlinSolvIters(cvode_mem, &nni);
  if( _check_cv_flag(&flag, "CVodeGetNumNonlinSolvIters", 1) ) return;
  flag = CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);
  if( _check_cv_flag(&flag, "CVodeGetNumNonlinSolvConvFails", 1) ) return;

  flag = CVodeGetNumJacEvals(cvode_mem, &nje);
  if( _check_cv_flag(&flag, "CVDlsGetNumJacEvals", 1) ) return;
  flag = CVodeGetNumRhsEvals(cvode_mem, &nfeLS);
  if( _check_cv_flag(&flag, "CVDlsGetNumRhsEvals", 1) ) return;

  os << "\nFinal Statistics:\n"
     << "   nst = " << nst << "   nfe  = " << nfe << "   nsetups = " << nsetups
     << "   nfeLS = " << nfeLS << "   nje = " << nje << std::endl
     << "   nni = " << nni << "   ncfn = " << ncfn << "   netf = " << netf
     << std::endl;

  if( sens ){
    long int nfSe, nfeS, nsetupsS, nniS, ncfnS, netfS;

    flag = CVodeGetSensNumRhsEvals(cvode_mem, &nfSe);
    if( _check_cv_flag(&flag, "CVodeGetSensNumRhsEvals", 1) ) return;
    flag = CVodeGetNumRhsEvalsSens(cvode_mem, &nfeS);
    if( _check_cv_flag(&flag, "CVodeGetNumRhsEvalsSens", 1) ) return;
    flag = CVodeGetSensNumLinSolvSetups(cvode_mem, &nsetupsS);
    if( _check_cv_flag(&flag, "CVodeGetSensNumLinSolvSetups", 1) ) return;
    flag = CVodeGetSensNumErrTestFails(cvode_mem, &netfS);
    if( _check_cv_flag(&flag, "CVodeGetSensNumErrTestFails", 1) ) return;
    flag = CVodeGetSensNumNonlinSolvIters(cvode_mem, &nniS);
    if( _check_cv_flag(&flag, "CVodeGetSensNumNonlinSolvIters", 1) ) return;
    flag = CVodeGetSensNumNonlinSolvConvFails(cvode_mem, &ncfnS);
    if( _check_cv_flag(&flag, "CVodeGetSensNumNonlinSolvConvFails", 1) ) return;

    os << "   nfSe = " << nfSe << "   nfeS  = " << nfeS
       << "   nsetupsS = " << nsetupsS << "   netfS = " << netfS
       << "   nniS = " << nniS << "   ncfnS = " << ncfnS << std::endl
       << std::endl;
  }
}

} // end namescape mc

#endif

