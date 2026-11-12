// Copyright (C) Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the EPL-2.0 with GPL-2.0-or-later as a Secondary License; see the LICENSE file.
// SPDX-License-Identifier: EPL-2.0 OR GPL-2.0-or-later

#ifndef CRONOS__SUNLINSOL_CRONOS_HPP
#define CRONOS__SUNLINSOL_CRONOS_HPP

// Direct sparse linear solvers for SUNDIALS' matrix-based solvers (KINSOL here) beyond the KLU SUNDIALS wraps itself:
// SuperLU (serial, the one Armadillo uses), UMFPACK and SPQR, each a custom SUNLinearSolver on a CSC SUNMatrix_Sparse --
// `setup` factorises, `solve` back-substitutes (SPQR: `setup` takes the values, `solve` factorises and solves at once).  Each backend is compiled where its flag is set:
//   CRONOS__WITH_SUPERLU   SuperLU's C API (slu_ddefs.h; -lsuperlu) -- set with CRONOS_WITH_SUPERLU (CMake), FLAG_SUPERLU (make)
//   CRONOS__WITH_UMFPACK   UMFPACK (umfpack.h; -lumfpack)
//   CRONOS__WITH_SPQR      SuiteSparseQR (SuiteSparseQR.hpp; -lspqr -lcholmod)
// SunLinSolCronos::available() says which are.  (10 Nov 2026, the KINSOL extension of OCFESLV, K1b.)

#include <vector>
#include <string>
#include <cstring>
#include <cmath>
#include <algorithm>
#include <sundials/sundials_linearsolver.h>
#include <sundials/sundials_context.h>
#include <sunmatrix/sunmatrix_sparse.h>
#include <nvector/nvector_serial.h>

#if defined(CRONOS__WITH_SUPERLU)
// ARMA_USE_SUPERLU is set by Armadillo's own configuration, read only once <armadillo> is included: include it first,
// so the choice below sees it (otherwise SuperLU's headers, included globally here, would empty Armadillo's copy)
#  if __has_include(<armadillo>)
#    include <armadillo>
#  endif
#  if defined(ARMA_USE_SUPERLU)
// Armadillo declares SuperLU's types and the functions used here inside arma::superlu, under the SAME include guards as
// SuperLU's headers -- including those as well would leave one of the two empty.  So Armadillo's declarations, its
// arma_wrapper() naming (the C symbol, or Armadillo's wrapper of it), and SuperMatrix stores filled here directly.
namespace cronos_slu = arma::superlu;
#    define CRONOS_SLU_FN( f ) arma::arma_wrapper( f )     // the functions in namespace arma, the types in arma::superlu
#  else
extern "C" {
#    include <slu_ddefs.h>          // SuperLU's own include directory, as CMake and makeoptions.mk give it
}
namespace cronos_slu {
  using ::SuperMatrix;  using ::NCformat;  using ::DNformat;  using ::superlu_options_t;  using ::SuperLUStat_t;  using ::GlobalLU_t;
  using ::SLU_NC;  using ::SLU_DN;  using ::SLU_D;  using ::SLU_GE;  using ::NOTRANS;  using ::COLAMD;
}
#    define CRONOS_SLU_FN( f ) ::f
#  endif
#endif
#if defined(CRONOS__WITH_UMFPACK)
#include <umfpack.h>
#endif
#if defined(CRONOS__WITH_SPQR)
#include <SuiteSparseQR.hpp>
#endif

namespace mc
{

struct SunLinSolCronos
{
  //! SUPERLU_SPQR: sparse LU, falling back to the rank-revealing QR for a setup whose LU reports a singular matrix --
  //! the tiers of OCFESLV's LM under SOLVE_SPQR (LU for the Newton step, QR when it cannot)
  enum Backend { SUPERLU = 0, UMFPACK, SPQR, SUPERLU_SPQR };

  //! @brief Whether @p b was compiled in
  static bool available( Backend b )
    {
      switch( b ){
#if defined(CRONOS__WITH_SUPERLU)
        case SUPERLU: return true;
#endif
#if defined(CRONOS__WITH_UMFPACK)
        case UMFPACK: return true;
#endif
#if defined(CRONOS__WITH_SPQR)
        case SPQR:    return true;
#endif
#if defined(CRONOS__WITH_SUPERLU) || defined(CRONOS__WITH_SPQR)
        case SUPERLU_SPQR: return true;       // with one of the two only, that one
#endif
        default:      return false;
      }
    }
  static char const* name( Backend b )
    { return b == SUPERLU? "SuperLU": b == UMFPACK? "UMFPACK": b == SPQR? "SPQR": "SuperLU/SPQR"; }

  //! @brief The solver on @p A's CSC pattern, @p b's backend; nullptr if @p b is not available or @p A is not CSC
  static SUNLinearSolver create( N_Vector y, SUNMatrix A, Backend b, SUNContext ctx );

  struct Content
  {
    Backend      be = SUPERLU;
    int          last = 0;
    sunindextype n = 0;
    bool         factored = false;
    bool         lu_ready = false, use_qr = false;   // SuperLU's factors held; this setup fell back to (or is) the QR
#if defined(CRONOS__WITH_SUPERLU)
    std::vector<int>          sl_colptr, sl_rowidx, perm_c, perm_r, etree;
    std::vector<double>       sl_val;
    cronos_slu::NCformat      sl_nc;          // A's store, ours (not SuperLU's to free)
    cronos_slu::SuperMatrix   sl_A, sl_AC, sl_L, sl_U;
    cronos_slu::GlobalLU_t    sl_Glu;
    bool                      sl_have_perm = false;
#endif
#if defined(CRONOS__WITH_UMFPACK)
    std::vector<int64_t> um_colptr, um_rowidx;
    void*                um_sym = nullptr;
    void*                um_num = nullptr;
#endif
#if defined(CRONOS__WITH_SPQR)
    cholmod_common                         sq_cc;
    bool                                   sq_started = false;
    cholmod_sparse*                        sq_A = nullptr;
    std::vector<double>                    sq_r, sq_c;      // the row and column scales of the matrix SPQR is given
#endif
  };

  static Content* C( SUNLinearSolver S ) { return static_cast<Content*>( S->content ); }

  static void release( Content* c )
    {
      if( !c ) return;
#if defined(CRONOS__WITH_SUPERLU)
      if( c->lu_ready ){                                // L, U and AC are SuperLU's; A's store is ours
        CRONOS_SLU_FN( Destroy_SuperNode_Matrix )( &c->sl_L );  CRONOS_SLU_FN( Destroy_CompCol_Matrix )( &c->sl_U );
        CRONOS_SLU_FN( Destroy_CompCol_Permuted )( &c->sl_AC );
        c->lu_ready = false;
      }
#endif
#if defined(CRONOS__WITH_UMFPACK)
      if( c->um_num ){ umfpack_dl_free_numeric( &c->um_num );  c->um_num = nullptr; }
#endif
#if defined(CRONOS__WITH_SPQR)
#endif
      c->factored = false;
    }

  static SUNLinearSolver_Type gettype( SUNLinearSolver ) { return SUNLINEARSOLVER_DIRECT; }
  static SUNLinearSolver_ID   getid  ( SUNLinearSolver ) { return SUNLINEARSOLVER_CUSTOM; }
  static SUNErrCode           initialize( SUNLinearSolver S ) { C( S )->last = 0;  return SUN_SUCCESS; }
  static sunindextype         lastflag( SUNLinearSolver S ) { return C( S )->last; }

#if defined(CRONOS__WITH_SPQR)
  //! @brief A's values (and, once, its pattern) into the cholmod matrix SPQR solves from
  static bool spqr_store( Content* c, sunindextype n, sunindextype nnz, sunindextype const* cp, sunindextype const* ri, double const* v )
    {
      if( !c->sq_started ){ cholmod_l_start( &c->sq_cc );  c->sq_started = true; }
      if( !c->sq_A || (sunindextype)c->sq_A->ncol != n || (sunindextype)c->sq_A->nzmax != nnz ){
        if( c->sq_A ) cholmod_l_free_sparse( &c->sq_A, &c->sq_cc );
        c->sq_A = cholmod_l_allocate_sparse( n, n, nnz, 1, 1, 0, CHOLMOD_REAL, &c->sq_cc );
        if( !c->sq_A ) return false;
        std::memcpy( c->sq_A->p, cp, sizeof( int64_t )*( n + 1 ) );
        std::memcpy( c->sq_A->i, ri, sizeof( int64_t )*nnz );
      }
      // EQUILIBRATED (11 Nov 2026): rows by 1/max|a_ij|, then columns likewise.  SPQR's rank tolerance is relative to
      // the largest column norm; on the unscaled Jacobian of a multi-physics model it discarded columns small only by
      // their units, the step became a minimum-norm one, not Newton's, and KINSOL stalled (OCFE_MBC5c at 1.1e-5).
      c->sq_r.assign( n, 0. );  c->sq_c.assign( n, 0. );
      for( sunindextype j = 0; j < n; ++j ) for( sunindextype e = cp[j]; e < cp[j+1]; ++e ) c->sq_r[ ri[e] ] = std::max( c->sq_r[ ri[e] ], std::fabs( v[e] ) );
      for( sunindextype i = 0; i < n; ++i ) c->sq_r[i] = c->sq_r[i] > 0.? 1./c->sq_r[i]: 1.;
      for( sunindextype j = 0; j < n; ++j ){
        double cm = 0.;  for( sunindextype e = cp[j]; e < cp[j+1]; ++e ) cm = std::max( cm, std::fabs( v[e] )*c->sq_r[ ri[e] ] );
        c->sq_c[j] = cm > 0.? 1./cm: 1.;
      }
      double* x = static_cast<double*>( c->sq_A->x );
      for( sunindextype j = 0; j < n; ++j ) for( sunindextype e = cp[j]; e < cp[j+1]; ++e ) x[e] = v[e]*c->sq_r[ ri[e] ]*c->sq_c[j];
      return true;
    }
  //! @brief x = the minimum-2-norm solution of A x = b, SPQR factorising and solving in one call
  static bool spqr_solve( Content* c, sunindextype n, double* xd, double const* bd )
    {
      cholmod_dense* B = cholmod_l_allocate_dense( n, 1, n, CHOLMOD_REAL, &c->sq_cc );
      double* bs = static_cast<double*>( B->x );
      for( sunindextype i = 0; i < n; ++i ) bs[i] = bd[i]*c->sq_r[i];                       // (R A C)(C^-1 x) = R b
      cholmod_dense* X = SuiteSparseQR_min2norm<double>( SPQR_ORDERING_DEFAULT, SPQR_DEFAULT_TOL, c->sq_A, B, &c->sq_cc );
      bool const ok = X && c->sq_cc.status == CHOLMOD_OK;
      if( ok ){ double const* xs = static_cast<double const*>( X->x );  for( sunindextype j = 0; j < n; ++j ) xd[j] = xs[j]*c->sq_c[j]; }
      cholmod_l_free_dense( &B, &c->sq_cc );  if( X ) cholmod_l_free_dense( &X, &c->sq_cc );
      return ok;
    }
#endif

  //! @brief Factorise A (its values, on the pattern fixed at the first call)
  static int setup( SUNLinearSolver S, SUNMatrix A )
    {
      Content* c = C( S );
      sunindextype const n = SUNSparseMatrix_Columns( A );
      sunindextype const* cp = SUNSparseMatrix_IndexPointers( A );
      sunindextype const* ri = SUNSparseMatrix_IndexValues( A );
      double const* v = SUNSparseMatrix_Data( A );
      sunindextype const nnz = cp[n];
      release( c );
      c->n = n;
      switch( c->be ){
#if defined(CRONOS__WITH_SUPERLU)
        case SUPERLU: case SUPERLU_SPQR: {
          c->use_qr = false;
          if( c->sl_colptr.size() != (size_t)n + 1 ){               // the pattern, once: 32-bit indices for SuperLU
            c->sl_colptr.assign( cp, cp + n + 1 );  c->sl_rowidx.assign( ri, ri + nnz );
            c->perm_c.assign( n, 0 );  c->perm_r.assign( n, 0 );  c->etree.assign( n, 0 );  c->sl_have_perm = false;
          }
          c->sl_val.assign( v, v + nnz );
          c->sl_nc.nnz = (int)nnz;  c->sl_nc.nzval = c->sl_val.data();  c->sl_nc.rowind = c->sl_rowidx.data();  c->sl_nc.colptr = c->sl_colptr.data();
          c->sl_A.Stype = cronos_slu::SLU_NC;  c->sl_A.Dtype = cronos_slu::SLU_D;  c->sl_A.Mtype = cronos_slu::SLU_GE;
          c->sl_A.nrow = (int)n;  c->sl_A.ncol = (int)n;  c->sl_A.Store = &c->sl_nc;
          cronos_slu::superlu_options_t opts;  CRONOS_SLU_FN( set_default_options )( &opts );  opts.ColPerm = cronos_slu::COLAMD;
          if( !c->sl_have_perm ){ CRONOS_SLU_FN( get_perm_c )( 3, &c->sl_A, c->perm_c.data() );  c->sl_have_perm = true; }   // COLAMD, the pattern's
          CRONOS_SLU_FN( sp_preorder )( &opts, &c->sl_A, c->perm_c.data(), c->etree.data(), &c->sl_AC );
          cronos_slu::SuperLUStat_t stat;  CRONOS_SLU_FN( StatInit )( &stat );
          int info = 0;
          CRONOS_SLU_FN( dgstrf )( &opts, &c->sl_AC, CRONOS_SLU_FN( sp_ienv )( 2 ), CRONOS_SLU_FN( sp_ienv )( 1 ), c->etree.data(), nullptr, 0,
                                   c->perm_c.data(), c->perm_r.data(), &c->sl_L, &c->sl_U, &c->sl_Glu, &stat, &info );
          CRONOS_SLU_FN( StatFree )( &stat );
          c->factored = true;  c->lu_ready = true;
          if( info == 0 ) break;
#if defined(CRONOS__WITH_SPQR)
          if( c->be == SUPERLU_SPQR && info > 0 ){                   // singular: this setup is the QR's
            release( c );  c->use_qr = true;
            if( !spqr_store( c, n, nnz, cp, ri, v ) ){ c->last = -1;  return SUN_ERR_EXT_FAIL; }
            c->factored = true;
            break;
          }
#endif
          c->last = (int)info;  return info > 0? SUNLS_PACKAGE_FAIL_REC: SUN_ERR_EXT_FAIL;
        }
#elif defined(CRONOS__WITH_SPQR)
        case SUPERLU_SPQR: {                                         // no SuperLU in this build: the QR alone
          c->use_qr = true;
          if( !spqr_store( c, n, nnz, cp, ri, v ) ){ c->last = -1;  return SUN_ERR_EXT_FAIL; }
          c->factored = true;
          break;
        }
#endif
#if defined(CRONOS__WITH_UMFPACK)
        case UMFPACK: {
          if( c->um_colptr.size() != (size_t)n + 1 ){
            c->um_colptr.assign( cp, cp + n + 1 );  c->um_rowidx.assign( ri, ri + nnz );
            if( c->um_sym ){ umfpack_dl_free_symbolic( &c->um_sym );  c->um_sym = nullptr; }
          }
          if( !c->um_sym && umfpack_dl_symbolic( n, n, c->um_colptr.data(), c->um_rowidx.data(), v, &c->um_sym, nullptr, nullptr ) != UMFPACK_OK ){
            c->last = -1;  return SUN_ERR_EXT_FAIL; }
          int const st = umfpack_dl_numeric( c->um_colptr.data(), c->um_rowidx.data(), v, c->um_sym, &c->um_num, nullptr, nullptr );
          c->factored = true;
          if( st != UMFPACK_OK ){ c->last = st;  return st == UMFPACK_WARNING_singular_matrix? SUNLS_PACKAGE_FAIL_REC: SUN_ERR_EXT_FAIL; }
          break;
        }
#endif
#if defined(CRONOS__WITH_SPQR)
        case SPQR: {              // the values are taken: SPQR factorises and solves in ONE call, in solve() -- Q is
          c->use_qr = true;         // never stored (as OCFESLV's LM calls it)
          if( !spqr_store( c, n, nnz, cp, ri, v ) ){ c->last = -1;  return SUN_ERR_EXT_FAIL; }
          c->factored = true;
          break;
        }
#endif
        default: c->last = -1;  return SUN_ERR_EXT_FAIL;
      }
      c->last = 0;
      return SUN_SUCCESS;
    }

  //! @brief x = A^-1 b, from the factorisation of the last setup
  static int solve( SUNLinearSolver S, SUNMatrix A, N_Vector x, N_Vector b, sunrealtype )
    {
      Content* c = C( S );
      if( !c->factored ){ c->last = -1;  return SUN_ERR_EXT_FAIL; }
      double* xd = N_VGetArrayPointer( x );
      double const* bd = N_VGetArrayPointer( b );
      sunindextype const n = c->n;
      switch( c->be ){
#if defined(CRONOS__WITH_SUPERLU)
        case SUPERLU: case SUPERLU_SPQR: {
#if defined(CRONOS__WITH_SPQR)
          if( c->use_qr ){ if( !spqr_solve( c, n, xd, bd ) ){ c->last = c->sq_cc.status;  return SUNLS_PACKAGE_FAIL_REC; }  break; }
#endif
          if( xd != bd ) std::memcpy( xd, bd, sizeof( double )*n );
          cronos_slu::DNformat dn;  dn.lda = (int)n;  dn.nzval = xd;                     // B's store, ours
          cronos_slu::SuperMatrix B;  B.Stype = cronos_slu::SLU_DN;  B.Dtype = cronos_slu::SLU_D;  B.Mtype = cronos_slu::SLU_GE;
          B.nrow = (int)n;  B.ncol = 1;  B.Store = &dn;
          cronos_slu::SuperLUStat_t stat;  CRONOS_SLU_FN( StatInit )( &stat );
          int info = 0;
          CRONOS_SLU_FN( dgstrs )( cronos_slu::NOTRANS, &c->sl_L, &c->sl_U, c->perm_c.data(), c->perm_r.data(), &B, &stat, &info );
          CRONOS_SLU_FN( StatFree )( &stat );
          if( info != 0 ){ c->last = info;  return SUN_ERR_EXT_FAIL; }
          break;
        }
#endif
#if defined(CRONOS__WITH_UMFPACK)
        case UMFPACK: {
          std::vector<double> rhs( bd, bd + n );
          int const st = umfpack_dl_solve( UMFPACK_A, c->um_colptr.data(), c->um_rowidx.data(), SUNSparseMatrix_Data( A ),
                                           xd, rhs.data(), c->um_num, nullptr, nullptr );
          if( st != UMFPACK_OK ){ c->last = st;  return SUNLS_PACKAGE_FAIL_REC; }
          break;
        }
#endif
#if defined(CRONOS__WITH_SPQR)
        case SPQR:
#if !defined(CRONOS__WITH_SUPERLU)
        case SUPERLU_SPQR:
#endif
        {
          if( !spqr_solve( c, n, xd, bd ) ){ c->last = c->sq_cc.status;  return SUNLS_PACKAGE_FAIL_REC; }
          break;
        }
#endif
        default: c->last = -1;  return SUN_ERR_EXT_FAIL;
      }
      c->last = 0;
      return SUN_SUCCESS;
    }

  static SUNErrCode free_( SUNLinearSolver S )
    {
      if( !S ) return SUN_SUCCESS;
      Content* c = C( S );
      if( c ){
        release( c );
#if defined(CRONOS__WITH_UMFPACK)
        if( c->um_sym ) umfpack_dl_free_symbolic( &c->um_sym );
#endif
#if defined(CRONOS__WITH_SPQR)
        if( c->sq_A ) cholmod_l_free_sparse( &c->sq_A, &c->sq_cc );
        if( c->sq_started ) cholmod_l_finish( &c->sq_cc );
#endif
        delete c;  S->content = nullptr;
      }
      SUNLinSolFreeEmpty( S );
      return SUN_SUCCESS;
    }
};

inline SUNLinearSolver
SunLinSolCronos::create
( N_Vector, SUNMatrix A, Backend b, SUNContext ctx )
{
  if( !available( b ) || !A || SUNMatGetID( A ) != SUNMATRIX_SPARSE || SUNSparseMatrix_SparseType( A ) != CSC_MAT ) return nullptr;
  SUNLinearSolver S = SUNLinSolNewEmpty( ctx );
  if( !S ) return nullptr;
  S->ops->gettype    = gettype;
  S->ops->getid      = getid;
  S->ops->initialize = initialize;
  S->ops->setup      = setup;
  S->ops->solve      = solve;
  S->ops->lastflag   = lastflag;
  S->ops->free       = free_;
  Content* c = new Content;  c->be = b;
  S->content = c;
  return S;
}

} // namespace mc

#endif
