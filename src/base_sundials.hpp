// Copyright (C) 2019 Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

#ifndef MC__BASE_SUNDIALS_HPP
#define MC__BASE_SUNDIALS_HPP

#include <iostream>
#include <iomanip>
#include <sys/time.h>

#include "cvodes.h"
#include "nvector_serial.h"
#include "cvodes_dense.h"
//#include "cvodes_lapack.h"
#include "cvodes_diag.h"
#include "sundials_dense.h"
#include "sundials_types.h"

namespace mc
{
//! @brief C++ base class for computing solutions of parametric DAEs using GSL
////////////////////////////////////////////////////////////////////////
//! mc::BASE_SUNDIALS is a C++ base class for computing solutions of
//! parametric differential-algebraic equations (DAEs) using SUNDIALS
////////////////////////////////////////////////////////////////////////
class BASE_SUNDIALS
{
public:
  /** @ingroup ODESLV_SUNDIALS
   *  @ingroup ODEBND_SUNDIALS
   *  @{
   */
  //! @brief Default class constructor
  BASE_SUNDIALS()
    {}

  //! @brief Class destructor
  virtual ~BASE_SUNDIALS()
    {}

  //! @brief GSL options
  struct Options
  {
    //! @brief Constructor
    Options():
      INTMETH(MSADAMS), H0(0e0), HMIN(1e-12), HMAX(0e0), NMAX(2000), RTOL(1e-8), ATOL(1e-8),
      ETOL(1e-16), QERR(true), MAXFAIL(10), MAXCORR(5), JACAPPROX(CV_DIAG),
      AUTOTOLS(false), RTOLS(1e-8), ATOLS(1e-8), ETOLS(1e-16), FSACORR(STAGGERED),
      FSAERR(true), QERRS(true), RTOLB(1e-8), ATOLB(1e-8), ETOLB(1e-16), QERRB(true),
      ASAINTERP(HERMITE), ASACHKPT(2000), RTOLFD(1e-3), ATOLFD(1e-3), CENFD(true)
      {}
    //! @brief Assignment operator
    template <typename U> Options& operator=
      ( U&options ){
        INTMETH   = options.INTMETH;
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
        JACAPPROX = options.JACAPPROX;
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
        RTOLFD    = options.RTOLFD;
        ATOLFD    = options.ATOLFD;
        CENFD     = options.CENFD;
        return *this;
      }
    //! @brief Enumeration type for FSA method
    enum FSA_STRATEGY{
      SIMULTANEOUS=CV_SIMULTANEOUS,//!< Simultaneous state/sensitivity correction
      STAGGERED=CV_STAGGERED,	//!< Simultaneous sensitivity corrections after state corrections
      STAGGERED1=CV_STAGGERED1	//!< Sequential sensitivity corrections after state corrections
    };
    //! @brief Enumeration type for ASA method
    enum ASA_STRATEGY{
      HERMITE=CV_HERMITE,	//!< Cubic Hermite interpolation
      POLYNOMIAL=CV_POLYNOMIAL	//!< Variable degree polynomial interpolation
    };
    //! @brief Enumeration of numerical integration algorithms
    enum INTEGRATION_METHOD{
      MSADAMS=2,	//!< Variable-coefficient linear multistep Adams method (non-stiff systems)
      MSBDF		//!< Variable-coefficient linear multistep backward differentiation formula (BDF) method (stiff systems)
    };
    //! @brief Enumeration type for FSA method
    enum JAC_STRATEGY{
      CV_DIAG=0,	//!< Approximate diagonal Jacobian formed by way of a difference quotient
      CV_DENSE,		//!< Approximate dense Jacobian and use of internal direct dense linear algebra functions
      CV_LAPACKDENSE	//!< Approximate dense Jacobian and use of LAPACK dense linear algebra functions
    };
    //! @brief Numerical integration method [Default: MSADAMS]
    INTEGRATION_METHOD INTMETH;
    //! @brief Initial step-size [Default: 0e0 (auto)]
    double H0;
    //! @brief Minimum step-size [Default: 0e0]
    double HMIN;
    //! @brief Maximum step-size [Default: 0e0 (+inf)]
    double HMAX;
    //! @brief Maximum number of steps in a time stage [Default: 0 (500)]
    unsigned int NMAX;
    //! @brief Relative (scalar) integration tolerance
    realtype RTOL;
    //! @brief Absolute (scalar) integration tolerance
    realtype ATOL;
    //! @brief Absolute (scalar) integration tolerance for ellipsoid shape matrix
    realtype ETOL;    
    //! @brief Whether or not (state) quadrature error control is performed?
    bool QERR;
    //! @brief Maximum number of error test failures per step
    int MAXFAIL;
    //! @brief Maximum  number of nonlinear solver iterations per step
    int MAXCORR;
    //! @brief Jacobian approximation method [Default: DIAG]
    JAC_STRATEGY JACAPPROX;
    //! @brief Whether integration tolerances for FS are to be set automatically?
    bool AUTOTOLS;
    //! @brief Relative (scalar) integration tolerance for FSA
    realtype RTOLS;
    //! @brief Absolute (scalar) integration tolerance for FSA
    realtype ATOLS;
    //! @brief Absolute (scalar) integration tolerance for ellipsoid shape matrix in FSA
    realtype ETOLS;      
    //! @brief FSA correction strategy
    FSA_STRATEGY FSACORR;
    //! @brief Whether or not FSA error control is performed?
    bool FSAERR;
    //! @brief Whether or not FSA quadrature error control is performed?
    bool QERRS;
    //! @brief Relative (scalar) integration tolerance for ASA
    realtype RTOLB;
    //! @brief Absolute (scalar) integration tolerance for ASA
    realtype ATOLB;
    //! @brief Absolute (scalar) integration tolerance for ellipsoid shape matrix in ASA
    realtype ETOLB;    
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
  /** @} */

protected:
  //! @brief Function to initialize GSL statistics
  static void _init_stats
    ( Stats&stats );

  //! @brief Function to finalize GSL statistics
  static void _final_stats
    ( Stats&stats );

  //! @brief Function to display integration statistics
  static void _print_stats
    ( const Stats&stats, std::ostream&os=std::cout );

  //! @brief Function to check function return values with CVODE
  static bool _check_cv_flag
    ( void *flagvalue, std::string funcname, int opt );

  //! @brief Function to display CVODE statistics
  static void _print_stats_cvode
    ( void *cvode_mem, bool sens=false, std::ostream&os=std::cout );

  //! @brief Private methods to block default compiler methods
  BASE_SUNDIALS(const BASE_SUNDIALS&);
  BASE_SUNDIALS& operator=(const BASE_SUNDIALS&);
};

inline
void
BASE_SUNDIALS::_init_stats
( Stats&stats )
{
  // Initialize statistics
  stats.reset();
  timeval time;
  gettimeofday(&time, 0) ;
  stats.cputime = - time.tv_sec - time.tv_usec*1e-6;
}

inline
void
BASE_SUNDIALS::_final_stats
( Stats&stats )
{
  // Get final CPU time
  timeval time;
  gettimeofday(&time, 0);
  stats.cputime += time.tv_sec + time.tv_usec*1e-6;
}

inline
void
BASE_SUNDIALS::_print_stats
( const Stats&stats, std::ostream&os )
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
BASE_SUNDIALS::_check_cv_flag
( void *flagvalue, std::string funcname, int opt )
{
  // Check if SUNDIALS function returned NULL pointer - no memory allocated
  if (opt == 0 && flagvalue == NULL) {
    std::cerr << "\nSUNDIALS_ERROR: " << funcname
              << "() failed - returned NULL pointer\n\n";
    return(true); }

  // Check if flag < 0
  else if (opt == 1) {
    int *errflag = (int *) flagvalue;
    if (*errflag < 0) {
      std::cerr << "\nSUNDIALS_ERROR: " << funcname
                << "() failed with flag = " << *errflag << "\n\n";
      return(true); }}

  // Check if function returned NULL pointer - no memory allocated
  else if (opt == 2 && flagvalue == NULL) {
    std::cerr << "\nMEMORY_ERROR: " << funcname
              << "() failed - returned NULL pointer\n\n";
    return(true); }

  return(false);
}

inline
void
BASE_SUNDIALS::_print_stats_cvode
( void *cvode_mem, bool sens, std::ostream&os )
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

  flag = CVDlsGetNumJacEvals(cvode_mem, &nje);
  if( _check_cv_flag(&flag, "CVDlsGetNumJacEvals", 1) ) return;
  flag = CVDlsGetNumRhsEvals(cvode_mem, &nfeLS);
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

