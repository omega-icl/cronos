// Copyright (C) Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the EPL-2.0 with GPL-2.0-or-later as a Secondary License; see the LICENSE file.
// SPDX-License-Identifier: EPL-2.0 OR GPL-2.0-or-later

#ifndef CRONOS__BASE_KINSOL_HPP
#define CRONOS__BASE_KINSOL_HPP

#include <cmath>
#include <string>
#include <vector>

#include <kinsol/kinsol.h>
#include <kinsol/kinsol_ls.h>
#include <nvector/nvector_serial.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunlinsol/sunlinsol_dense.h>
#if defined (CRONOS__WITH_KLU)
  #include <sunmatrix/sunmatrix_sparse.h>
  #include <sunlinsol/sunlinsol_klu.h>
#endif

#include "base_sundials.hpp"

namespace mc
{
//! @brief C++ base class wrapping SUNDIALS KINSOL for the consistent initial-condition solve
////////////////////////////////////////////////////////////////////////
//! mc::BASE_KINSOL is the corrector behind an initial-condition SOLVE: the INITIAL rows together with
//! the INTERIOR rows at t0, in the unknowns (x, dx/dt).  It is shared -- an ODE consumer with implicit
//! initial rows wants it as much as a DAE consumer does -- and it deliberately owns NO options object:
//! its knobs are a plain struct a solver embeds in its own INITIAL group, so there is one `options` per
//! solver, as everywhere else.
//!
//! WHAT IT KNOWS THAT PLAIN KINSOL DOES NOT (plan §3.12).  Newton converges LINEARLY at a rank-deficient
//! root -- the error halves per iteration, so ~24 iterations reach ~6e-8 and then stall -- and no
//! corrector fixes that; it needs rank-revealing linear algebra.  The pendulum released horizontally is
//! such a root BY CONSTRUCTION (the position constraint's gradient collapses onto the row pinning x), and
//! the DOF COUNT passes while the Jacobian is singular.  So a stall is REPORTED as a rank diagnosis, not
//! retried, and the caller decides whether to hand off to a rank-revealing path.
//!
//! Also handled explicitly: starting AT a root.  KINSOL's linear-solver setup fails (-11) on a zero
//! residual, and a reinitialisation after a discontinuity frequently starts at or near consistency.
////////////////////////////////////////////////////////////////////////
class BASE_KINSOL
: public virtual BASE_SUNDIALS
{
public:
  //! @brief The knobs of the IC solve, to be EMBEDDED in a solver's INITIAL option group
  struct Knobs
  {
    //! @brief stopping tolerance on the scaled residual norm
    double FTOL = 1e-10;
    //! @brief stopping tolerance on the scaled step -- what stops a slowly converging iteration
    double STEPTOL = 1e-14;
    //! @brief iteration cap
    int MAXIT = 200;
    //! @brief use the line-search globalisation (KIN_LINESEARCH) rather than plain Newton (KIN_NONE).
    //! Made no measurable difference on the §3.12 systems; kept because it costs nothing where it does not help
    bool LINESEARCH = true;
    //! @brief Jacobian refresh interval: 1 = every iteration (KINSetMaxSetupCalls)
    int MAXSETUP = 1;
    //! @brief use the sparse KLU linear solver (needs CRONOS__WITH_KLU) rather than dense LU
    bool SPARSE = false;
    //! @brief iteration count above which a CONVERGED solve is read as linearly convergent, i.e. the
    //! rank-deficiency signature.  Newton at a simple root is quadratic -- a handful of iterations from any
    //! sane guess (4 on the §3.12 full-rank system from far away) -- while at a double root it halves the
    //! error each step.  The residual is NO guide: there F ~ y^2, so |F| meets FTOL while the SOLUTION is
    //! only accurate to sqrt(FTOL).  Measured: 17 iterations to |F| = 1e-10 with |y| ~ 1e-5.
    int SLOW_ITERS = 12;
  };

  //! @brief How the solve ended
  enum Status {
    OK = 0,          //!< converged to FTOL
    AT_ROOT,         //!< the starting point already satisfied FTOL: nothing was solved
    RANK_DEFICIENT,  //!< converged LINEARLY (or stalled): the Jacobian is rank-deficient at this root and the
                     //!< solution is accurate only to about sqrt(FTOL), whatever the residual says
    FAILED           //!< KINSOL returned an error that is not the plateau signature
  };

  //! @brief What happened, for the caller's diagnosis
  struct Report
  {
    Status status = FAILED;
    int kin_flag = 0;          //!< KINSOL's own return code
    long iterations = 0;
    double fnorm = 0.;         //!< final scaled residual norm
    std::string message;
  };

  //! @brief Residual callback: F(y) into f, both of size n.  Returns 0 on success.
  typedef int (*ResFn)( int n, double const* y, double* f, void* user );
  //! @brief Dense Jacobian callback, column-major J[i + n*j] = dF_i/dy_j.  Returns 0 on success.
  typedef int (*JacFn)( int n, double const* y, double* J, void* user );

  //! @brief Solve F(y) = 0 from the starting point in @p y, which holds the solution on return
  static Report solve
    ( int const n, ResFn res, JacFn jac, void* user, std::vector<double>& y, Knobs const& k );

protected:
  BASE_KINSOL() {}
  virtual ~BASE_KINSOL() {}

private:
  struct Shim { ResFn res; JacFn jac; void* user; int n; std::vector<double> Jbuf; };
  static int _res( N_Vector uu, N_Vector ff, void* p )
    { Shim* s = static_cast<Shim*>( p );
      return s->res( s->n, N_VGetArrayPointer( uu ), N_VGetArrayPointer( ff ), s->user ); }
  static int _jac( N_Vector uu, N_Vector, SUNMatrix J, void* p, N_Vector, N_Vector )
    { Shim* s = static_cast<Shim*>( p );
      int const rc = s->jac( s->n, N_VGetArrayPointer( uu ), s->Jbuf.data(), s->user );
      if( rc ) return rc;
      for( int j = 0; j < s->n; ++j ) for( int i = 0; i < s->n; ++i ) SM_ELEMENT_D( J, i, j ) = s->Jbuf[ i + s->n*j ];
      return 0; }
};

inline
BASE_KINSOL::Report
BASE_KINSOL::solve
( int const n, ResFn res, JacFn jac, void* user, std::vector<double>& y, Knobs const& k )
{
  Report R;
  // AT the root already?  KINSOL's setup fails (-11) on a zero residual, and a reinitialisation after a
  // discontinuity often starts consistent; this is the common case, not the odd one.
  { std::vector<double> f( n );
    if( res( n, y.data(), f.data(), user ) == 0 ){
      double nrm = 0.; for( double const v : f ) nrm = std::max( nrm, std::fabs( v ) );
      if( nrm <= k.FTOL ){ R.status = AT_ROOT; R.fnorm = nrm; R.message = "already consistent"; return R; }
    } }

  SUNContext ctx = nullptr;
  if( SUNContext_Create( SUN_COMM_NULL, &ctx ) < 0 ){ R.message = "SUNContext_Create failed"; return R; }
  N_Vector uu = N_VNew_Serial( n, ctx ), sc = N_VNew_Serial( n, ctx );
  N_VConst( 1., sc );
  for( int i = 0; i < n; ++i ) N_VGetArrayPointer( uu )[i] = y[i];
  Shim shim{ res, jac, user, n, std::vector<double>( (size_t)n * n ) };

  void* mem = KINCreate( ctx );
  KINInit( mem, _res, uu );
  KINSetUserData( mem, &shim );
  KINSetFuncNormTol( mem, k.FTOL );
  KINSetScaledStepTol( mem, k.STEPTOL );
  KINSetNumMaxIters( mem, k.MAXIT );
  KINSetMaxSetupCalls( mem, k.MAXSETUP );
  SUNMatrix A = nullptr; SUNLinearSolver LS = nullptr;
#if defined (CRONOS__WITH_KLU)
  if( k.SPARSE ){ A = SUNSparseMatrix( n, n, n*n, CSC_MAT, ctx ); LS = SUNLinSol_KLU( uu, A, ctx ); }
  else
#endif
  { A = SUNDenseMatrix( n, n, ctx ); LS = SUNLinSol_Dense( uu, A, ctx ); }
  KINSetLinearSolver( mem, LS, A );
  if( jac && !k.SPARSE ) KINSetJacFn( mem, _jac );        // dense analytic; sparse analytic is a later step

  R.kin_flag = KINSol( mem, uu, k.LINESEARCH? KIN_LINESEARCH: KIN_NONE, sc, sc );
  KINGetNumNonlinSolvIters( mem, &R.iterations );
  KINGetFuncNorm( mem, &R.fnorm );
  for( int i = 0; i < n; ++i ) y[i] = N_VGetArrayPointer( uu )[i];

  bool const converged = ( R.kin_flag == KIN_SUCCESS || R.kin_flag == KIN_INITIAL_GUESS_OK );
  if( converged && R.iterations <= k.SLOW_ITERS ){ R.status = OK; R.message = "converged"; }
  else if( R.kin_flag == KIN_LINSOLV_NO_RECOVERY || R.kin_flag == KIN_LSETUP_FAIL ){
    // The LU setup itself failed: the Jacobian is SINGULAR at this point, not merely ill-conditioned.
    // That is the same rank deficiency as the plateau below, seen earlier and more plainly.
    R.status = RANK_DEFICIENT;
    R.message = "the linear solve failed (" + std::to_string( R.kin_flag ) + "): the initial-condition"
                " Jacobian is SINGULAR at this configuration -- the free initial data is degenerate here";
  }
  else if( converged || R.kin_flag == KIN_STEP_LT_STPTOL || R.kin_flag == KIN_MAXITER_REACHED ){
    // the §3.12 signature: the error has been HALVING rather than squaring.  At a double root F ~ y^2, so a
    // residual that meets FTOL says nothing about the solution, which is accurate only to ~sqrt(FTOL).
    // The Jacobian is rank-deficient at this root; no corrector does better, and the caller should know.
    R.status = RANK_DEFICIENT;
    R.message = std::string( converged? "converged LINEARLY": "stalled" ) + " (" + std::to_string( R.iterations )
              + " iterations, |F| = " + std::to_string( R.fnorm ) + "): the initial-condition Jacobian is"
                " rank-deficient at this configuration -- the free initial data is degenerate here, and the"
                " solution is accurate only to about sqrt(FTOL)";
  }
  else { R.status = FAILED; R.message = "KINSol returned " + std::to_string( R.kin_flag ); }

  SUNLinSolFree( LS ); SUNMatDestroy( A ); N_VDestroy( uu ); N_VDestroy( sc );
  KINFree( &mem ); SUNContext_Free( &ctx );
  return R;
}

} // end namespace mc

#endif
