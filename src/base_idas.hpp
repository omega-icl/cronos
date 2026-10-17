// Copyright (C) Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the EPL-2.0 with GPL-2.0-or-later as a Secondary License; see the LICENSE file.
// SPDX-License-Identifier: EPL-2.0 OR GPL-2.0-or-later

#ifndef CRONOS__BASE_IDAS_HPP
#define CRONOS__BASE_IDAS_HPP

#include <iostream>
#include <iomanip>
#include <cassert>
#include <thread>
#include <chrono>     // elapsed time of an integration (portable; was gettimeofday)

#include <idas/idas.h>                  /* prototypes for IDAS fcts., consts.   */
#include <idas/idas_ls.h>               /* access to the IDAS linear solver i/f */
#include <nvector/nvector_serial.h>     /* access to serial N_Vector            */
#include <sunmatrix/sunmatrix_dense.h>  /* access to dense SUNMatrix            */
#if defined (CRONOS__WITH_KLU)
  #include <sunmatrix/sunmatrix_sparse.h> /* access to sparse SUNMatrix           */
  #include <sunlinsol/sunlinsol_klu.h>    /* access to sparse SUNLinearSolver     */
#endif
#include <sunlinsol/sunlinsol_dense.h>  /* access to dense SUNLinearSolver      */
#include "sunnonlinsol/sunnonlinsol_newton.h"     /* access to the newton SUNNonlinearSolver      */
#include "sunnonlinsol/sunnonlinsol_fixedpoint.h" /* access to the fixed point SUNNonlinearSolver */

#include "base_sundials.hpp"

namespace mc
{
//! @brief C++ base class for computing solutions of parametric differential-algebraic equations using SUNDIALS IDAS
////////////////////////////////////////////////////////////////////////
//! mc::BASE_IDAS is the IDAS counterpart of mc::BASE_CVODES: options, statistics, thread
//! registration and flag checking for the DAE integrators.  It mirrors BASE_CVODES field for
//! field WHERE THE TWO SOLVERS AGREE, and departs from it where IDAS has no counterpart --
//! every departure is commented, because a silently-dropped option is worse than an absent one:
//!
//!   * NO `DIAG` linear solver.  CVDiag has no IDAS analogue at all, and it is the CVODES
//!     DEFAULT, so the IDAS default is DENSE instead.
//!   * NO `MSADAMS`.  IDAS is BDF-only, so INTEGRATION_METHOD admits MSBDF alone and an attempt
//!     to ask for Adams is a COMPILE error rather than a value quietly ignored at run time.
//!   * NO `STAGGERED1` FSA strategy.  IDAS offers IDA_SIMULTANEOUS and IDA_STAGGERED only.
//!   * EXTRA options with no CVODES counterpart: SUPPRESSALG, and the consistent-initialisation
//!     block (CALCIC, TOUT1, NMAXIC, MAXJACIC, MAXNLIIC, LSOFFIC, STEPTOLIC) driving IDACalcIC.
////////////////////////////////////////////////////////////////////////
class BASE_IDAS
: public virtual BASE_SUNDIALS
{
public:
  /** @{
   */
  //! @brief Default class constructor
  BASE_IDAS()
    : _is_registered( false )
    {}

  //! @brief Class destructor
  virtual ~BASE_IDAS()
    {}

  //! @brief Function registering thread
  void registration
    ();

  //! @brief Function re-registering thread
  void reregistration
    ();

  //! @brief Function unregistering thread
  void unregistration
    ();

  //! @brief IDAS options
  struct Options
  {
    enum FSA_STRATEGY{
      SIMULTANEOUS=IDA_SIMULTANEOUS, //!< Simultaneous state/sensitivity correction
      STAGGERED=IDA_STAGGERED        //!< Sensitivity corrections after state corrections
    };
    enum ASA_STRATEGY{
      HERMITE=IDA_HERMITE,	//!< Cubic Hermite interpolation
      POLYNOMIAL=IDA_POLYNOMIAL	//!< Variable degree polynomial interpolation
    };
    enum INTEGRATION_METHOD{
      MSBDF=0	//!< Variable-coefficient linear multistep backward differentiation formula (BDF) method
    };
    enum NONLINEAR_SOLVER{
      FIXEDPOINT=0,	//!< Fixed point nonlinear solver
      NEWTON		//!< Newton nonlinear solver
    };
    enum LINEAR_SOLVER{
      DENSE=0,	//!< Use analytic dense Jacobian and internal direct dense linear algebra functions
      DENSEDQ,	//!< Use approximate dense Jacobian by way of a difference quotient and internal direct dense linear algebra
      SPARSE	//!< Use analytic sparse Jacobian and the KLU sparse direct solver.  Requires a build with KLU
		//!< (CRONOS__WITH_KLU); otherwise the integrator REFUSES it at setup
    };
    enum CALCIC_MODE{
      NONE=0,	//!< Take the supplied (y0, yp0) as consistent; do not call IDACalcIC
      YA_YDP,	//!< IDA_YA_YDP_INIT: solve for the algebraic components of y and for y', given the
		//!< differential components of y.  Requires the id vector (IDASetId)
      Y		//!< IDA_Y_INIT: solve for y given y', holding y' fixed
    };

    //! @brief Everything the integration takes: method, linear and nonlinear solver, step control, tolerances,
    //! and the sensitivity and adjoint settings.  Grouped as OCFESLV groups its SOLVE, which has no
    //! SENS/ADJOINT substructures either.
    struct t_Solve
    {
      //! @brief Enumeration type for FSA method.  NOTE: IDAS has NO STAGGERED1 (CVODES does)
      //! @brief Enumeration type for ASA method
      //! @brief Enumeration of numerical integration algorithms.  NOTE: IDAS is BDF-ONLY -- there is
      //! no MSADAMS here on purpose, so asking for Adams does not compile rather than being ignored
      //! @brief Enumeration of nonlinear solver strategies
      //! @brief Enumeration of linear solver strategies (within Newton nonlinear solver).  NOTE: there
      //! is NO DIAG here -- CVDiag has no IDAS counterpart, which is why the default is DENSE
      //! @brief Enumeration of consistent-initialisation modes, passed to IDACalcIC
      //! @brief Numerical integration method [Default: MSBDF -- the only one IDAS has]
      INTEGRATION_METHOD INTMETH = MSBDF;
      //! @brief Nonlinear solver method [Default: NEWTON]
      NONLINEAR_SOLVER NLINSOL = NEWTON;
      //! @brief Linear solver method and Jacobian approximation [Default: DENSE -- CVODES says DIAG,
      //! which IDAS does not have]
      LINEAR_SOLVER LINSOL = DENSE;
      //! @brief Initial step-size [Default: 0e0 (auto)]
      double H0 = 0e0;
      //! @brief Minimum step-size [Default: 0e0]
      double HMIN = 1e-12;
      //! @brief Maximum step-size [Default: 0e0 (+inf)]
      double HMAX = 0e0;
      //! @brief Maximum number of steps in a time stage [Default: 0 (2000)]
      unsigned int NMAX = 2000;
      //! @brief Maximum BDF order, 1 to 5 [Default: 5] (IDASetMaxOrd)
      int MAXORD = 5;
      //! @brief Relative (scalar) integration tolerance
      sunrealtype RTOL = 1e-7;
      //! @brief Absolute (scalar) integration tolerance
      sunrealtype ATOL = 1e-9;
      //! @brief Absolute (scalar) integration tolerance for ellipsoid shape matrix
      sunrealtype ETOL = 1e-16;
      //! @brief Whether or not (state) quadrature error control is performed?
      bool QERR = true;
      //! @brief Maximum number of error test failures per step
      int MAXFAIL = 10;
      //! @brief Maximum number of nonlinear solver iterations per step
      int MAXCORR = 5;
      //! @brief Whether to EXCLUDE the algebraic variables from the local error test [Default: false]
      //! (IDASetSuppressAlg).  No CVODES counterpart -- there are no algebraic variables there
      bool SUPPRESSALG = false;
      //! @brief Whether integration tolerances for FS are to be set automatically?
      bool AUTOTOLS = false;
      //! @brief Relative (scalar) integration tolerance for FSA
      sunrealtype RTOLS = 1e-7;
      //! @brief Absolute (scalar) integration tolerance for FSA
      sunrealtype ATOLS = 1e-9;
      //! @brief Absolute (scalar) integration tolerance for ellipsoid shape matrix in FSA
      sunrealtype ETOLS = 1e-14;
      //! @brief FSA correction strategy
      FSA_STRATEGY FSACORR = STAGGERED;
      //! @brief Whether or not FSA error control is performed?
      bool FSAERR = true;
      //! @brief Whether or not FSA quadrature error control is performed?
      bool QERRS = true;
      //! @brief Relative (scalar) integration tolerance for ASA
      sunrealtype RTOLB = 1e-7;
      //! @brief Absolute (scalar) integration tolerance for ASA
      sunrealtype ATOLB = 1e-9;
      //! @brief Absolute (scalar) integration tolerance for ellipsoid shape matrix in ASA
      sunrealtype ETOLB = 1e-16;
      //! @brief Whether or not ASA quadrature error control is performed?
      bool QERRB = true;
      //! @brief ASA interpolation strategy
      ASA_STRATEGY ASAINTERP = HERMITE;
      //! @brief Number of steps between each check point for ASA
      int ASACHKPT = 2000;
    } SOLVE;

    //! @brief Consistent initialisation: the IDACalcIC controls, and where a KINSOL initial-condition solve's
    //! own options belong.
    struct t_Initial
    {
      //! @brief Consistent-initialisation mode [Default: YA_YDP] (IDACalcIC)
      CALCIC_MODE CALCIC = YA_YDP;
      //! @brief First output time handed to IDACalcIC, which only uses its DIRECTION and rough
      //! magnitude [Default: 0e0 = let the integrator pick from the stage length]
      double TOUT1 = 0e0;
      //! @brief Maximum number of steps in the IC calculation [Default: 5] (IDASetMaxNumStepsIC)
      int NMAXIC = 5;
      //! @brief Maximum number of Jacobian re-evaluations in the IC calculation [Default: 4] (IDASetMaxNumJacsIC)
      int MAXJACIC = 4;
      //! @brief Maximum number of Newton iterations per IC step [Default: 10] (IDASetMaxNumItersIC)
      int MAXNLIIC = 10;
      //! @brief Whether to switch the IC line search OFF [Default: false] (IDASetLineSearchOffIC)
      bool LSOFFIC = false;
      //! @brief Step tolerance of the IC calculation [Default: 0e0 (auto)] (IDASetStepToleranceIC)
      double STEPTOLIC = 0e0;
    } INITIAL;

    //! @brief Assignment operator
    template <typename U> Options& operator=
      ( U& options ){ SOLVE = options.SOLVE; INITIAL = options.INITIAL; return *this; }
  };

  //! @brief Structure storing integration statistics
  struct Stats
  {
    //! @brief Constructor
    Stats():
      cputime(0.), numSteps(0), numRES(0), numJAC(0), numIC(0)
      {}
    //! @brief Constructor
    void reset()
      { cputime = 0.; numSteps = numRES = numJAC = numIC = 0; }

    //! @brief CPU time
    double cputime;
    //! @brief Number of integration steps
    unsigned int numSteps;
    //! @brief Number of residual evaluations (the DAE counterpart of CVODES' numRHS)
    unsigned int numRES;
    //! @brief Number of Jacobian evaluations
    unsigned int numJAC;
    //! @brief Number of consistent-initialisation solves (no CVODES counterpart)
    unsigned int numIC;
  };

  //! @brief Thread-safe static pointer to self
  thread_local static BASE_IDAS* PTR_BASE_IDAS;
  /** @} */

protected:
  //! @brief Flag indicating whether the thread is registered
  bool _is_registered;

  //! @brief Pure virtual function to calculate the DAE residual
  virtual int IDRES__
    ( sunrealtype tt, N_Vector yy, N_Vector yp, N_Vector rr, void* user_data )
    = 0;

  //! @brief Pure virtual function to calculate the quadrature RHS derivatives
  virtual int IDQUAD__
    ( sunrealtype tres, N_Vector yy, N_Vector yp, N_Vector rrQ, void* user_data )
    = 0;

  //! @brief Pure virtual function to calculate the DAE residual Jacobian dF/dy + c_j dF/dy'
  virtual int IDJAC__
    ( sunrealtype t, sunrealtype c_j, N_Vector y, N_Vector yp, N_Vector r, SUNMatrix Jac,
      void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3 )
    = 0;

  //! @brief Pure virtual function to calculate the adjoint DAE residual
  virtual int IDRESB__
    ( sunrealtype tt, N_Vector yy, N_Vector yp, N_Vector yyB, N_Vector ypB, N_Vector rrB,
      void* user_dataB )
    = 0;

  //! @brief Pure virtual function to calculate the adjoint quadrature RHS derivatives
  virtual int IDQUADB__
    ( sunrealtype tt, N_Vector yy, N_Vector yp, N_Vector yyB, N_Vector ypB, N_Vector rhsvalBQ,
      void* user_dataB )
    = 0;

  //! @brief Pure virtual function to calculate the adjoint DAE residual Jacobian
  virtual int IDJACB__
    ( sunrealtype tt, sunrealtype c_jB, N_Vector yy, N_Vector yp, N_Vector yyB, N_Vector ypB,
      N_Vector rrB, SUNMatrix JacB, void* user_dataB, N_Vector tmp1B, N_Vector tmp2B,
      N_Vector tmp3B )
    = 0;

  //! @brief Pure virtual function to calculate the sensitivity DAE residuals
  virtual int IDRESF__
    ( int Ns, sunrealtype t, N_Vector yy, N_Vector yp, N_Vector resval, N_Vector* yyS,
      N_Vector* ypS, N_Vector* resvalS, void* user_data, N_Vector tmp1, N_Vector tmp2,
      N_Vector tmp3 )
    = 0;

  //! @brief Pure virtual function to calculate the sensitivity quadrature RHS derivatives
  virtual int IDQUADF__
    ( int Ns, sunrealtype t, N_Vector yy, N_Vector yp, N_Vector* yyS, N_Vector* ypS,
      N_Vector rrQ, N_Vector* rhsvalQS, void* user_data, N_Vector yytmp, N_Vector yptmp,
      N_Vector tmpQS )
    = 0;

  //! @brief Function to initialize IDAS statistics
  static void _init_stats
    ( Stats& stats );

  //! @brief Function to finalize IDAS statistics
  static void _final_stats
    ( Stats& stats );

  //! @brief Function to display integration statistics
  static void _print_stats
    ( Stats const& stats, std::ostream& os=std::cout );

  //! @brief Function to check function return values with IDAS
  static bool _check_ida_flag
    ( void *flagvalue, std::string const& funcname, int const opt );

  //! @brief Function to display IDAS statistics
  static void _print_stats_idas
    ( void* ida_mem, bool const sens=false, std::ostream& os=std::cout );

  //! @brief Private methods to block default compiler methods
  BASE_IDAS( BASE_IDAS const& ) = delete;
  BASE_IDAS& operator=( BASE_IDAS const& ) = delete;
};

inline thread_local BASE_IDAS* BASE_IDAS::PTR_BASE_IDAS = nullptr;
inline thread_local int (BASE_IDAS::*PTR_IDRES)  ( sunrealtype, N_Vector, N_Vector, N_Vector, void* ) = nullptr;
inline thread_local int (BASE_IDAS::*PTR_IDQUAD) ( sunrealtype, N_Vector, N_Vector, N_Vector, void* ) = nullptr;
inline thread_local int (BASE_IDAS::*PTR_IDJAC)  ( sunrealtype, sunrealtype, N_Vector, N_Vector, N_Vector, SUNMatrix, void*, N_Vector, N_Vector, N_Vector ) = nullptr;
inline thread_local int (BASE_IDAS::*PTR_IDRESB) ( sunrealtype, N_Vector, N_Vector, N_Vector, N_Vector, N_Vector, void* ) = nullptr;
inline thread_local int (BASE_IDAS::*PTR_IDQUADB)( sunrealtype, N_Vector, N_Vector, N_Vector, N_Vector, N_Vector, void* ) = nullptr;
inline thread_local int (BASE_IDAS::*PTR_IDJACB) ( sunrealtype, sunrealtype, N_Vector, N_Vector, N_Vector, N_Vector, N_Vector, SUNMatrix, void*, N_Vector, N_Vector, N_Vector ) = nullptr;
inline thread_local int (BASE_IDAS::*PTR_IDRESF) ( int, sunrealtype, N_Vector, N_Vector, N_Vector, N_Vector*, N_Vector*, N_Vector*, void*, N_Vector, N_Vector, N_Vector ) = nullptr;
inline thread_local int (BASE_IDAS::*PTR_IDQUADF)( int, sunrealtype, N_Vector, N_Vector, N_Vector*, N_Vector*, N_Vector, N_Vector*, void*, N_Vector, N_Vector, N_Vector ) = nullptr;

inline
void
REG_BASE_IDAS
( BASE_IDAS* ID,
  int (BASE_IDAS::*IDResFn)  ( sunrealtype, N_Vector, N_Vector, N_Vector, void* ),
  int (BASE_IDAS::*IDQuadFn) ( sunrealtype, N_Vector, N_Vector, N_Vector, void* ),
  int (BASE_IDAS::*IDJacFn)  ( sunrealtype, sunrealtype, N_Vector, N_Vector, N_Vector, SUNMatrix, void*, N_Vector, N_Vector, N_Vector ),
  int (BASE_IDAS::*IDResFnB) ( sunrealtype, N_Vector, N_Vector, N_Vector, N_Vector, N_Vector, void* ),
  int (BASE_IDAS::*IDQuadFnB)( sunrealtype, N_Vector, N_Vector, N_Vector, N_Vector, N_Vector, void* ),
  int (BASE_IDAS::*IDJacFnB) ( sunrealtype, sunrealtype, N_Vector, N_Vector, N_Vector, N_Vector, N_Vector, SUNMatrix, void*, N_Vector, N_Vector, N_Vector ),
  int (BASE_IDAS::*IDResFnF) ( int, sunrealtype, N_Vector, N_Vector, N_Vector, N_Vector*, N_Vector*, N_Vector*, void*, N_Vector, N_Vector, N_Vector ),
  int (BASE_IDAS::*IDQuadFnF)( int, sunrealtype, N_Vector, N_Vector, N_Vector*, N_Vector*, N_Vector, N_Vector*, void*, N_Vector, N_Vector, N_Vector ) )
{
#ifdef CRONOS__BASE_IDAS_CHECK
  assert( BASE_IDAS::PTR_BASE_IDAS == nullptr && PTR_IDRES == nullptr );
#endif
  BASE_IDAS::PTR_BASE_IDAS = ID;
  PTR_IDRES   = IDResFn;
  PTR_IDQUAD  = IDQuadFn;
  PTR_IDJAC   = IDJacFn;
  PTR_IDRESB  = IDResFnB;
  PTR_IDQUADB = IDQuadFnB;
  PTR_IDJACB  = IDJacFnB;
  PTR_IDRESF  = IDResFnF;
  PTR_IDQUADF = IDQuadFnF;
}

inline
void
BASE_IDAS::registration
()
{
  if( _is_registered ) return;
  REG_BASE_IDAS( this,
                 &BASE_IDAS::IDRES__,
                 &BASE_IDAS::IDQUAD__,
                 &BASE_IDAS::IDJAC__,
                 &BASE_IDAS::IDRESB__,
                 &BASE_IDAS::IDQUADB__,
                 &BASE_IDAS::IDJACB__,
                 &BASE_IDAS::IDRESF__,
                 &BASE_IDAS::IDQUADF__
               );
  _is_registered = true;
}

inline
void
BASE_IDAS::reregistration
()
{
  assert( _is_registered );
  REG_BASE_IDAS( this,
                 &BASE_IDAS::IDRES__,
                 &BASE_IDAS::IDQUAD__,
                 &BASE_IDAS::IDJAC__,
                 &BASE_IDAS::IDRESB__,
                 &BASE_IDAS::IDQUADB__,
                 &BASE_IDAS::IDJACB__,
                 &BASE_IDAS::IDRESF__,
                 &BASE_IDAS::IDQUADF__
               );
}

inline
void
UNREG_BASE_IDAS
()
{
#ifdef CRONOS__BASE_IDAS_CHECK
  assert( BASE_IDAS::PTR_BASE_IDAS != nullptr
       && PTR_IDRES   != nullptr
       && PTR_IDQUAD  != nullptr
       && PTR_IDJAC   != nullptr
       && PTR_IDRESB  != nullptr
       && PTR_IDQUADB != nullptr
       && PTR_IDJACB  != nullptr
       && PTR_IDRESF  != nullptr
       && PTR_IDQUADF != nullptr );
#endif
  BASE_IDAS::PTR_BASE_IDAS = nullptr;
  PTR_IDRES   = nullptr;
  PTR_IDQUAD  = nullptr;
  PTR_IDJAC   = nullptr;
  PTR_IDRESB  = nullptr;
  PTR_IDQUADB = nullptr;
  PTR_IDJACB  = nullptr;
  PTR_IDRESF  = nullptr;
  PTR_IDQUADF = nullptr;
}

inline
void
BASE_IDAS::unregistration
()
{
  if( !_is_registered ) return;
  UNREG_BASE_IDAS();
  _is_registered = false;
}

inline
void
BASE_IDAS::_init_stats
( Stats& stats )
{
  // Initialize statistics
  stats.reset();
  stats.cputime = - std::chrono::duration<double>( std::chrono::steady_clock::now().time_since_epoch() ).count();
}

inline
void
BASE_IDAS::_final_stats
( Stats& stats )
{
  // Get final CPU time
  stats.cputime += std::chrono::duration<double>( std::chrono::steady_clock::now().time_since_epoch() ).count();
}

inline
void
BASE_IDAS::_print_stats
( Stats const& stats, std::ostream& os )
{
  // Statistics -- SAME layout as BASE_CVODES::_print_stats, with RES for RHS and the IC count
  // appended, so the two registers read alike in a driver that runs both
  os << " No STEPS    " << stats.numSteps
     << std::endl
     << " No EVALATIONS" << "   RES: " << stats.numRES
                         << "   JAC: " << stats.numJAC
                         << "   IC: "  << stats.numIC
     << std::endl
     << " CPU TIME (SEC)     " << std::fixed << std::left
                               << std::setprecision(5) << stats.cputime
     << std::endl << std::endl;
  return;
}

inline
bool
BASE_IDAS::_check_ida_flag
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
BASE_IDAS::_print_stats_idas
( void* ida_mem, bool const sens, std::ostream& os )
{
  long int nst, nre, nsetups, nje, nreLS, nni, ncfn, netf;

  int flag = IDAGetNumSteps(ida_mem, &nst);
  if( _check_ida_flag(&flag, "IDAGetNumSteps", 1) ) return;
  flag = IDAGetNumResEvals(ida_mem, &nre);
  if( _check_ida_flag(&flag, "IDAGetNumResEvals", 1) ) return;
  flag = IDAGetNumLinSolvSetups(ida_mem, &nsetups);
  if( _check_ida_flag(&flag, "IDAGetNumLinSolvSetups", 1) ) return;
  flag = IDAGetNumErrTestFails(ida_mem, &netf);
  if( _check_ida_flag(&flag, "IDAGetNumErrTestFails", 1) ) return;
  flag = IDAGetNumNonlinSolvIters(ida_mem, &nni);
  if( _check_ida_flag(&flag, "IDAGetNumNonlinSolvIters", 1) ) return;
  flag = IDAGetNumNonlinSolvConvFails(ida_mem, &ncfn);
  if( _check_ida_flag(&flag, "IDAGetNumNonlinSolvConvFails", 1) ) return;

  flag = IDAGetNumJacEvals(ida_mem, &nje);
  if( _check_ida_flag(&flag, "IDAGetNumJacEvals", 1) ) return;
  flag = IDAGetNumLinResEvals(ida_mem, &nreLS);
  if( _check_ida_flag(&flag, "IDAGetNumLinResEvals", 1) ) return;

  os << "\nFinal Statistics:\n"
     << "   nst = " << nst << "   nre  = " << nre << "   nsetups = " << nsetups
     << "   nreLS = " << nreLS << "   nje = " << nje << std::endl
     << "   nni = " << nni << "   ncfn = " << ncfn << "   netf = " << netf
     << std::endl;

  if( sens ){
    long int nfSe, nreS, nsetupsS, nniS, ncfnS, netfS;

    flag = IDAGetSensNumResEvals(ida_mem, &nfSe);
    if( _check_ida_flag(&flag, "IDAGetSensNumResEvals", 1) ) return;
    flag = IDAGetNumResEvalsSens(ida_mem, &nreS);
    if( _check_ida_flag(&flag, "IDAGetNumResEvalsSens", 1) ) return;
    flag = IDAGetSensNumLinSolvSetups(ida_mem, &nsetupsS);
    if( _check_ida_flag(&flag, "IDAGetSensNumLinSolvSetups", 1) ) return;
    flag = IDAGetSensNumErrTestFails(ida_mem, &netfS);
    if( _check_ida_flag(&flag, "IDAGetSensNumErrTestFails", 1) ) return;
    flag = IDAGetSensNumNonlinSolvIters(ida_mem, &nniS);
    if( _check_ida_flag(&flag, "IDAGetSensNumNonlinSolvIters", 1) ) return;
    flag = IDAGetSensNumNonlinSolvConvFails(ida_mem, &ncfnS);
    if( _check_ida_flag(&flag, "IDAGetSensNumNonlinSolvConvFails", 1) ) return;

    os << "   nfSe = " << nfSe << "   nreS  = " << nreS
       << "   nsetupsS = " << nsetupsS << "   netfS = " << netfS
       << "   nniS = " << nniS << "   ncfnS = " << ncfnS << std::endl
       << std::endl;
  }
}

// 2026-09-28: over-alignment guard (GCC wrong-code with over-aligned virtual bases)
static_assert( alignof( BASE_IDAS ) <= alignof( void* ), "BASE_IDAS is a VIRTUAL base of the CRONOS solvers and must not be over-aligned (alignof > 8): GCC (6 to at least 16) emits aligned vector stores in base-object constructors assuming the full alignment, but a virtual-base subobject is only placed at its non-virtual alignment -> SIGSEGV at -O2/-O3 (see gccvb/pr_vbase_align.cpp). Keep over-aligned members (Armadillo/Eigen fixed-size types, alignas) behind a pointer, as FFModel::_pClassification does." );

} // end namescape mc

#endif
