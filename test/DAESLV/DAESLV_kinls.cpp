////////////////////////////////////////////////////////////////////////
// DAESLV_kinls -- the sparse direct solvers KINSOL can factorise with (10 Nov 2026, the KINSOL extension of OCFESLV,
// K1b): BASE_KINSOL::Knobs::SPARSE_SOLVER = KLU (SUNDIALS' own, the default), or SuperLU, UMFPACK, SPQR as CRONOS'
// own SUNLinearSolvers (src/sunlinsol_cronos.hpp), each where its build flag is set -- and SUPERLU_SPQR, sparse LU with
// the QR for a setup whose LU is singular (the tiers of OCFESLV's LM under SOLVE_SPQR).
//   (A) 2-D Bratu, -lap u = 5 e^u (5-point differences, 1600 unknowns), through BASE_KINSOL: every compiled backend
//       converges to KLU's root (to 1e-12), in KLU's iteration count;
//   (B) DAESLV under the method of lines (the Peclet-10 rod, its boundary nodes algebraic): the consistent
//       initialisation with each backend (options.INITIAL.KIN.SPARSE_SOLVER) gives KLU's outputs to 1e-10.
////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <iomanip>
#include <cmath>
#include "daeslv_idas.hpp"
using namespace mc;

static int nfail = 0;
static void check( std::string const& w, bool ok ){ std::cout << "  " << std::left << std::setw(96) << w << ( ok? "PASS": "FAIL" ) << "\n"; if( !ok ) ++nfail; }
static char const* NM[5] = { "KLU", "SuperLU", "UMFPACK", "SPQR", "SuperLU/SPQR" };
static bool compiled( int b ){ return b == 0 || SunLinSolCronos::available( (SunLinSolCronos::Backend)( b - 1 ) ); }

static int M = 40;
static int bratu( int n, double const* y, double* f, void* ){
  double const h = 1./( M + 1 ), ih2 = 1./( h*h );
  for( int i = 0; i < M; ++i ) for( int j = 0; j < M; ++j ){
    int const k = i*M + j;
    double const l = i > 0? y[k-M]: 0., r = i < M-1? y[k+M]: 0., d = j > 0? y[k-1]: 0., u = j < M-1? y[k+1]: 0.;
    f[k] = ( 4.*y[k] - l - r - d - u )*ih2 - 5.*std::exp( y[k] ); }
  return 0; }
static std::vector<sunindextype> CP, RI;
static int bratu_jac( int n, double const* y, double* val, void* ){
  double const h = 1./( M + 1 ), ih2 = 1./( h*h );
  for( int c = 0; c < n; ++c ) for( sunindextype e = CP[c]; e < CP[c+1]; ++e ) val[e] = (int)RI[e] == c? 4.*ih2 - 5.*std::exp( y[c] ): -ih2;
  return 0; }

static bool rod( int backend, std::vector<double>& F ){
  FFGraph G; DAESLV_IDAS S( &G ); FFPartial OpP; FFEval OpE; FFIntegral OpI;
  FFVar t = G.add_var("t"), z = G.add_var("z"), u = G.add_var("u");
  S.add_domain( t, FFDom( 0., 1., 2, FFDom::LGR, 3 ) ); S.add_domain( z, FFDom( 0., 1., 6, FFDom::LGL, 5 ) );
  S.set_evolution_domain( t ); S.add_state( u, {t,z} ); S.update_ref( u, 0.5 );
  int const TI = FFDom::ALL - FFDom::LB, ZI = FFDom::ALL - FFDom::LB - FFDom::UB;
  DAESLV_IDAS::EqnOptions io( DAESLV_IDAS::EqnRole::INTERIOR, 0 ), ib( DAESLV_IDAS::EqnRole::BOUNDARY, 0 ), ii( DAESLV_IDAS::EqnRole::INITIAL, 0 );
  S.add_equation( OpP(u,t) + OpP(u,z) - 0.1*OpP(OpP(u,z),z), {t,z}, {TI, ZI}, io );
  S.add_equation( u - 1., {t,z}, {TI, FFDom::LB}, ib );  S.add_equation( OpP(u,z), {t,z}, {TI, FFDom::UB}, ib );
  S.add_equation( u - 0.5 - 0.5*(1.-z)*(1.-z), {t,z}, {FFDom::LB, FFDom::ALL}, ii );
  S.add_output( OpE( u, {{t,1},{z,1}}, {{t,1.},{z,1.}} ) );  S.add_output( OpI( u, z ), {t}, {1.} );
  S.options.INTERFACE.IMPOSITION = OCFEMOL::Options::IC_STRONG;
  S.options.INITIAL.KIN.SPARSE_SOLVER = (BASE_KINSOL::Knobs::SparseSolver)backend;
  S.options.SOLVE.RTOL = S.options.SOLVE.ATOL = 1e-11;  S.options.SOLVE.VERBOSE = 0;  S.options.DISPLAY_LEVEL = 0;
  if( !S.setup() || S.solve( std::vector<FFModel::InputVal>{} ) != DAESLV_IDAS::NORMAL ) return false;
  F = S.val_function();  return true;
}

int main()
{
  std::cout << std::scientific << std::setprecision(3);
  std::cout << "(A) 2-D Bratu through BASE_KINSOL, " << M*M << " unknowns: each compiled backend against KLU\n";
  { int const n = M*M;  CP.assign( n + 1, 0 );  RI.clear();
    for( int c = 0; c < n; ++c ){ int const i = c / M, j = c % M;
      if( i > 0 ) RI.push_back( c - M );  if( j > 0 ) RI.push_back( c - 1 );  RI.push_back( c );
      if( j < M-1 ) RI.push_back( c + 1 );  if( i < M-1 ) RI.push_back( c + M );  CP[c+1] = RI.size(); }
    BASE_KINSOL::Sparsity sp;  sp.colptr = CP;  sp.rowidx = RI;  sp.jac = bratu_jac;
    std::vector<double> ref;  long itref = 0;
    for( int b = 0; b < 5; ++b ){
      if( !compiled( b ) ){ std::cout << "    " << NM[b] << ": not in this build\n";  continue; }
      BASE_KINSOL::Knobs k;  k.SPARSE = true;  k.SPARSE_SOLVER = (BASE_KINSOL::Knobs::SparseSolver)b;  k.FTOL = 1e-9;  k.MAXIT = 50;
      std::vector<double> y( n, 0. );
      auto r = BASE_KINSOL::solve( n, bratu, nullptr, nullptr, y, k, &sp );
      bool const ok = r.status == BASE_KINSOL::OK || r.status == BASE_KINSOL::AT_ROOT;
      double d = 0.;  if( b == 0 ){ ref = y;  itref = r.iterations; } else for( int i = 0; i < n; ++i ) d = std::max( d, std::fabs( y[i] - ref[i] ) );
      std::cout << "    " << NM[b] << ": " << r.iterations << " iterations, max |u - u_KLU| " << d << "\n";
      check( std::string( "(A) " ) + NM[b] + ": converged to KLU's root to 1e-12, in KLU's iteration count", ok && !ref.empty() && d < 1e-12 && r.iterations == itref );
    }
  }
  std::cout << "(B) DAESLV under MoL, the consistent initialisation with each backend (INITIAL.KIN.SPARSE_SOLVER)\n";
  { std::vector<double> Fr;  bool const okr = rod( 0, Fr );
    check( "(B) KLU: solves", okr );
    for( int b = 1; b < 5; ++b ){
      if( !compiled( b ) ){ std::cout << "    " << NM[b] << ": not in this build\n";  continue; }
      std::vector<double> F;  bool const ok = rod( b, F );
      double d = 0.;  if( ok && okr ) for( size_t i = 0; i < F.size(); ++i ) d = std::max( d, std::fabs( F[i] - Fr[i] )/std::max( 1., std::fabs( Fr[i] ) ) );
      std::cout << "    " << NM[b] << ": outputs against KLU's " << d << "\n";
      check( std::string( "(B) " ) + NM[b] + ": KLU's outputs to 1e-10", ok && okr && d < 1e-10 );
    }
  }
  std::cout << "==========  DAESLV_kinls: " << ( nfail? std::to_string( nfail ) + " FAILED": "ALL PASS" ) << "  ==========\n";
  return nfail? 1: 0;
}
