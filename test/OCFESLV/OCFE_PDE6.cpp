// Copyright (C) Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the EPL-2.0 with GPL-2.0-or-later as a Secondary License; see the LICENSE file.
// SPDX-License-Identifier: EPL-2.0 OR GPL-2.0-or-later
//
// OCFE_PDE6_solve2.cpp
// --------------------
// Two independent PDE blocks in a single OCFESLV, refactored from the legacy
// test4_multiblock.cpp into the OCFE_PDEx_solve2.cpp driver format (config
// macros, Par struct, exact-profile functions, the generic LM/Newton harness,
// a per-mode run_mode() and an IC_WEAK/IC_TRACE/IC_STRONG acceptance sweep).
//
//   Block 0: parabolic heat equation in conservative first-order form
//     T_t - q_x = 0,      q - a*T_x = 0
//     T(0,x) = Ts + T0*sin(pi/2 * x/xf)
//     T(t,0) = Ts,        q(t,xf) = 0
//
//   Block 1: hyperbolic linear advection
//     u_t + c*u_x = 0
//     u(0,x) = u0*sin(2*pi * x/xf)
//     u(t,0) = -u0*sin(2*pi*c*t/xf)            (c > 0)
//
// Unlike the manufactured PDEx cases, the reference profiles here are the genuine
// analytical solutions (sin/exp), which are NOT polynomials.  The collocation
// solution therefore differs from the reference by spectral truncation error, so
// the reference-profile residual is a small-but-nonzero truncation level (not
// ~1e-14) and the exact-match tolerances are spectral, not machine.
//
// Purpose: exercise block-level classification and per-block IC_AUTO dispatch
// (parabolic block 0 -> IC_VALUE; hyperbolic block 1 -> IC_UPWIND on spatial
// interfaces; time interfaces stay IC_VALUE), and confirm Jacobian correctness
// via the finite-difference checks in test_deriv_utils.hpp -- including at a
// state with a deliberate element-local jump in the hyperbolic field u.
//
// Build knobs (all overridable with -D):
//   -DTEST_HA_NEL_T   axial(time) finite elements        (default 2)
//   -DTEST_HA_NEL_X   spatial finite elements            (default 2)
//   -DTEST_HA_NT      nodes per time element (LGR)       (default 8)
//   -DTEST_HA_NX      nodes per space element (LGL)      (default 8)
//   -DTEST_HA_INIT_MODE   0 exact, 1 perturbed exact, 2 constant fields (default 1)
//   -DTEST_HA_MAXIT / -DTEST_HA_SOLVE_TOL / -DTEST_HA_SAT_SIGMA0
//   -DTEST_HA_CLASSIFY    run/report classify_pde()      (default 1)
//   -DTEST_HA_DERIV_CHECK run the FD Jacobian checks      (default 1)

#include <algorithm>
#include <cctype>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <string>
#include <vector>
#include <fstream>

#include <armadillo>

#ifndef OCFE_OCFESLV_HEADER
#define OCFE_OCFESLV_HEADER "ocfeslv.hpp"
#endif

//#define MC__OCFESLV_INTERFACE_DECISION_V2
//#define MC__OCFESLV_STRONG_EXPLICIT_TAU
//#define MC__OCFESLV_STRONG_DUMMY_DIAG

#include OCFE_OCFESLV_HEADER

#include "test_deriv_utils.hpp"

#ifndef TEST_HA_NEL_T
#define TEST_HA_NEL_T 5
#endif
#ifndef TEST_HA_NEL_X
#define TEST_HA_NEL_X 5
#endif
#ifndef TEST_HA_NT
#define TEST_HA_NT 5
#endif
#ifndef TEST_HA_NX
#define TEST_HA_NX 5
#endif
#ifndef TEST_HA_INIT_MODE
#define TEST_HA_INIT_MODE 1
#endif
#ifndef TEST_HA_MAXIT
#define TEST_HA_MAXIT 40
#endif
#ifndef TEST_HA_SOLVE_TOL
#define TEST_HA_SOLVE_TOL 1e-9
#endif
#ifndef TEST_HA_SAT_SIGMA0
#define TEST_HA_SAT_SIGMA0 10.0
#endif
#ifndef TEST_HA_CLASSIFY
#define TEST_HA_CLASSIFY 1
#endif
#ifndef TEST_HA_DERIV_CHECK
#define TEST_HA_DERIV_CHECK 1
#endif
#ifndef TEST_HA_EXACT_TOL
#define TEST_HA_EXACT_TOL 5e-3
#endif
// Element-interface duplicate-spread gate (refinement-based).  A continuity that
// is explicitly imposed is machine-tight, but one that is *implied* (a soundly
// dropped/redundant claim, e.g. PDE6 q-flux) or *weak* (IC_WEAK SAT) holds only
// to the discretization truncation, which MR.ref_res (residual of the exact
// profile on the discrete system) measures.  So gate the spread against
// max(absolute floor, factor*ref_res): the original strict floor still applies
// to well-resolved explicit continuity, while truncation-level implied/weak
// continuity passes, and the bound tightens automatically under refinement (a
// genuine O(1) discontinuity still fails at any resolution).
#ifndef TEST_HA_DUP_ABS_FLOOR
#define TEST_HA_DUP_ABS_FLOOR 1e-6
#endif
#ifndef TEST_HA_DUP_TRUNC_FACTOR
#define TEST_HA_DUP_TRUNC_FACTOR 3.0
#endif
#ifndef TEST_HA_OUTPUT_PREFIX
#define TEST_HA_OUTPUT_PREFIX "OCFE_PDE6_solve"
#endif

using namespace mc;

static constexpr double kPi = 3.14159265358979323846264338327950288;

// ---------------------------------------------------------------------------
// Model parameters and analytical solution
// ---------------------------------------------------------------------------
struct Par
{
  double a   = 0.1;   // thermal diffusivity
  double T0  = 10.0;  // initial temperature amplitude
  double Ts  = 5.0;   // boundary/surface temperature
  double c   = 0.5;   // advection speed (>0)
  double u0  = 1.0;   // advection amplitude
  double xf  = 1.0;   // spatial extent
  double tf  = 1.0;   // time horizon
};

static double T_exact( double t, double x, Par const& p )
{
  double const k = kPi/2.0/p.xf;
  return p.Ts + p.T0 * std::exp( -p.a*k*k*t ) * std::sin( k*x );
}
static double q_exact( double t, double x, Par const& p )
{
  double const k = kPi/2.0/p.xf;
  return p.a * p.T0 * k * std::exp( -p.a*k*k*t ) * std::cos( k*x );
}
static double u_exact( double t, double x, Par const& p )
{
  return p.u0 * std::sin( 2.0*kPi/p.xf * ( x - p.c*t ) );
}

// Reference value for a collocation state, keyed by state name and node coords.
// The node coordinate order matches the domain order used in add_state below:
// every state lives on {t,x}, so xy = { t, x }.
static double exact_value_for_state( std::string const& nm, std::vector<double> const& xy, Par const& p )
{
  double const t = xy[0], x = xy[1];
  switch( nm.empty() ? '?' : nm[0] ){
    case 'T': return T_exact( t, x, p );
    case 'q': return q_exact( t, x, p );
    case 'u': return u_exact( t, x, p );
    default:  return 0.0;
  }
}
static double constant_initial_value_for_state( std::string const& nm, Par const& p )
{
  switch( nm.empty() ? '?' : nm[0] ){
    case 'T': return p.Ts;   // flat field at the surface temperature
    case 'q': return 0.0;
    case 'u': return 0.0;
    default:  return 0.0;
  }
}

// ---------------------------------------------------------------------------
// Generic residual / Jacobian / LM-Newton harness  (shared with OCFE_PDEx)
// ---------------------------------------------------------------------------
static double max_abs( std::vector<double> const& r )
{ double m=0.; for( double v: r ) m=std::max(m,std::abs(v)); return m; }

static void print_residuals( std::string const& label, std::vector<double> const& r )
{
  double sum=0.; for( double v: r ) sum += std::abs(v);
  std::cout << std::left << std::setw(46) << label
            << " n=" << std::setw(5) << r.size()
            << " max|r|=" << std::scientific << std::setprecision(4) << max_abs(r)
            << " mean|r|=" << (r.empty()?0.:sum/r.size()) << "\n";
}
static void print_largest_residuals( std::string const& label, std::vector<double> const& r, size_t nprint=10 )
{
  std::vector<size_t> idx(r.size());
  std::iota(idx.begin(),idx.end(),size_t(0));
  std::partial_sort(idx.begin(),idx.begin()+std::min(nprint,idx.size()),idx.end(),
    [&](size_t a,size_t b){return std::abs(r[a])>std::abs(r[b]);});
  std::cout << label << " largest residual rows:\n";
  for( size_t k=0; k<std::min(nprint,idx.size()); ++k ){
    size_t const i=idx[k];
    std::cout << "  row " << std::setw(6) << i << "  r=" << std::scientific << std::setprecision(8) << r[i] << "\n";
  }
}
static bool check_close( std::string const& label, double val, double tol )
{
  bool ok = std::isfinite(val) && val <= tol;
  std::cout << std::left << std::setw(50) << label
            << " value=" << std::scientific << std::setprecision(6) << val
            << " tol=" << tol << "  " << (ok?"PASS":"FAIL") << "\n";
  return ok;
}

#ifdef MC__OCFESLV_ACOND_PROBE
// Dense Jacobian assembly (sparse AD pattern + values via deriv(), densified).
// Retained only for the conditioning probe below; the production solve uses the
// sparse path inside OCFESLV::solve().
static bool dense_jacobian( OCFESLV& oc, std::vector<double> const& x, arma::mat& J )
{
  size_t const nEqn=oc.n_colloc_eqn(), nVar=oc.n_colloc_sta();
  std::vector<size_t> nnz(nEqn,0);
  std::vector<std::vector<size_t>> col(nEqn);
  std::vector<size_t*> colptr(nEqn,nullptr);
  if( !oc.deriv(nnz.data(),nullptr) ) return false;
  for( size_t i=0; i<nEqn; ++i ){ col[i].assign(nnz[i],0); if(nnz[i]) colptr[i]=col[i].data(); }
  if( !oc.deriv(nnz.data(),colptr.data()) ) return false;
  size_t nnzsum=0; for( size_t v: nnz ) nnzsum += v;
  std::vector<double> grad(nnzsum,0.0);
  if( !oc.deriv(grad.data(),nullptr,x.data(),nullptr,nullptr) ) return false;
  J.zeros(nEqn,nVar);
  size_t off=0;
  for( size_t i=0; i<nEqn; ++i ){
    for( size_t k=0; k<nnz[i]; ++k ) if( col[i][k] < nVar ) J(i,col[i][k]) = grad[off+k];
    off += nnz[i];
  }
  return true;
}
// A-invertibility probe (Schur Option B gate): is the (1,1) physical block A of
// the IC_TRACE saddle [A Bt; B C] invertible, so S = B A^{-1} Bt exists?  The
// suite is linear => J is state-independent (evaluate at x=0).  The physical/
// trace split is at trace_var_offset = n_colloc_sta() - n_colloc_trace().  For
// IC_STRONG the taus are eliminated (n_colloc_trace()==0), so J is already the
// reduced operator and there is no explicit A/B block to split out of var[].
static void probe_saddle_conditioning( OCFESLV& oc, std::string const& tag, int attempt )
{
  size_t const nEqn = oc.n_colloc_eqn();
  size_t const nVar = oc.n_colloc_sta();
  size_t const nTr  = oc.n_colloc_trace();
  std::cerr << "[Acond " << tag << " attempt=" << attempt << "] nEqn=" << nEqn
            << " nVar=" << nVar << " nTrace=" << nTr
            << " trace_var_offset=" << ( nVar>=nTr ? nVar-nTr : size_t(0) ) << "\n";
  if( nEqn != nVar ){ std::cerr << "  non-square system; skipping\n"; return; }

  std::vector<double> x0( nVar, 0.0 );
  arma::mat J;
  if( !dense_jacobian( oc, x0, J ) ){ std::cerr << "  dense_jacobian failed\n"; return; }

  auto report = [&]( char const* nm, arma::mat const& M ){
    if( M.n_rows==0 || M.n_cols==0 ){ std::cerr << "    " << nm << ": empty\n"; return; }
    arma::vec sv;
    if( !arma::svd( sv, M ) ){ std::cerr << "    " << nm << ": svd failed\n"; return; }
    double const smax = sv.n_elem ? sv.max() : 0.0;
    double const smin = sv.n_elem ? sv.min() : 0.0;
    double const tol  = (double)std::max(M.n_rows,M.n_cols) * smax
                      * std::numeric_limits<double>::epsilon();
    arma::uword rnk = 0; for( arma::uword i=0;i<sv.n_elem;++i) if( sv(i)>tol ) ++rnk;
    arma::uword const mn = (arma::uword)std::min(M.n_rows,M.n_cols);
    std::cerr << "    " << nm << ": " << M.n_rows << "x" << M.n_cols
              << "  sigma_max=" << smax << "  sigma_min=" << smin
              << "  cond=" << ( smin>0 ? smax/smin
                                       : std::numeric_limits<double>::infinity() )
              << "  numrank=" << rnk << "/" << mn
              << ( rnk==mn ? "  [full rank]" : "  [RANK DEFICIENT]" ) << "\n";
  };

  report( "J  (full system, must be nonsingular)", J );
  if( nTr>0 && nTr<nVar ){
    size_t const tvo = nVar - nTr;
    arma::mat const A  = J.submat( 0,   0,   tvo-1,  tvo-1 );
    arma::mat const B  = J.submat( tvo, 0,   nVar-1, tvo-1 );
    arma::mat const Bt = J.submat( 0,   tvo, tvo-1,  nVar-1 );
    arma::mat const C  = J.submat( tvo, tvo, nVar-1, nVar-1 );
    report( "A  (1,1 physical; Schur needs this INVERTIBLE)", A );
    report( "B  (2,1 continuity rows)", B );
    double const ntb = Bt.n_elem ? arma::abs(Bt).max() : 0.0;
    double const ntc = C.n_elem  ? arma::abs(C).max()  : 0.0;
    std::cerr << "    ||Bt(1,2)||_max=" << ntb << "  ||C(2,2)||_max=" << ntc
              << "  (pure saddle => C ~ 0)\n";
  } else {
    std::cerr << "    (IC_STRONG: taus eliminated; J is the reduced operator, "
                 "no explicit A/B split in var[])\n";
  }
}
#endif

static bool eval_residual( OCFESLV& oc, std::vector<double> const& x, std::vector<double>& r )
{ std::fill(r.begin(),r.end(),0.0); return oc.eval(r.data(),nullptr,x.data(),nullptr,nullptr); }


// ---------------------------------------------------------------------------
// Per-state exact-error reporting and element-interface duplicate spread
// ---------------------------------------------------------------------------
struct StateExactErrors { double eT=0.0, eq=0.0, eu=0.0, emax=0.0; };

static StateExactErrors print_variable_exact_errors
( OCFESLV const& oc, std::vector<double> const& var, Par const& p, std::string const& label )
{
  std::cout << "\nVariable comparison against exact profiles (" << label << "):\n";
  StateExactErrors E;
  size_t off=0;
  for( auto const& st : oc.states_colloc() ){
    auto nodes = oc.node_colloc(st);
    std::string const nm = st.name();
    double emax=0.0;
    for( size_t i=0; i<nodes.size(); ++i )
      emax = std::max( emax, std::abs( var[off+i] - exact_value_for_state(nm,nodes[i],p) ) );
    std::cout << "  " << std::setw(10) << nm
              << " max|num-exact|=" << std::scientific << std::setprecision(6) << emax << "\n";
    switch( nm.empty()?'?':nm[0] ){
      case 'T': E.eT = emax; break;
      case 'q': E.eq = emax; break;
      case 'u': E.eu = emax; break;
      default: break;
    }
    E.emax = std::max( E.emax, emax );
    off += nodes.size();
  }
  return E;
}

static double max_duplicate_spread
( OCFESLV const& oc, std::vector<double> const& var, std::string const& label )
{
  std::cout << "Element-interface duplicate-node spreads (" << label << "):\n";
  double max_pair=0.0;
  size_t off=0;
  for( auto const& st : oc.states_colloc() ){
    auto nodes = oc.node_colloc(st);
    std::map<std::vector<long long>, std::pair<double,double>> g; // key -> (lo,hi)
    for( size_t i=0; i<nodes.size(); ++i ){
      std::vector<long long> key; key.reserve(nodes[i].size());
      for( double cc: nodes[i] ) key.push_back( (long long)std::llround(cc*1.0e12) );
      double const v=var[off+i];
      auto it=g.find(key);
      if( it==g.end() ) g.emplace(std::move(key),std::make_pair(v,v));
      else { it->second.first=std::min(it->second.first,v); it->second.second=std::max(it->second.second,v); }
    }
    double sp=0.0; for( auto const& kv: g ) sp=std::max(sp,kv.second.second-kv.second.first);
    max_pair=std::max(max_pair,sp);
    std::cout << "  " << std::setw(10) << st.name()
              << " max_pair_spread=" << std::scientific << std::setprecision(6) << sp << "\n";
    off += nodes.size();
  }
  return max_pair;
}

// ---------------------------------------------------------------------------
// One imposition mode
// ---------------------------------------------------------------------------
struct ModeResult
{
  std::string name;
  bool        ok        = false;
  bool        solved    = false;
  bool        deriv_ok  = false;
  size_t      nVar      = 0;
  size_t      nEqn      = 0;
  size_t      nTrace    = 0;
  double      ref_res   = 0.0;
  double      final_res = std::numeric_limits<double>::infinity();
  double      max_spread= 0.0;
  double      eT=0.0, eq=0.0, eu=0.0, emax=0.0;
  std::string pde_type  = "n/a";
};

// Physical collocation-node coordinates of a domain (element nodes concatenated;
// element-interface nodes appear once per adjoining element, as in PDE2/3/4).
static std::vector<double> physical_nodes( FFDom const& dom )
{
  std::vector<double> xs;
  for( size_t ie=0; ie<dom.n_elem; ++ie ){
    auto xe = dom.lgnodes( dom.elem_lo(ie), dom.elem_up(ie) );
    xs.insert( xs.end(), xe.begin(), xe.end() );
  }
  return xs;
}
// Evaluate a collocated state at a physical (t,x) point from the solution vector.
static double interp( OCFESLV const& oc, FFVar const& v,
                      std::map<FFVar,double,lt_FFVar> const& pt,
                      std::vector<double> const& var )
{ return oc.eval_colloc<double>( v, pt, var.data(), nullptr, nullptr ); }

static ModeResult run_mode( OCFESLV::Options::ImpositionType imp, std::string const& strimp )
{
  ModeResult MR; MR.name = strimp;
  bool ok = true;
  Par p;

  std::cout << "\n========== two-block heat(parabolic)+advection(hyperbolic) test ==========\n";
  std::cout << "imposition: " << strimp
            << ", finite elements: t=" << TEST_HA_NEL_T << " x=" << TEST_HA_NEL_X
            << ", nodes/element: t=" << TEST_HA_NT << " x=" << TEST_HA_NX
            << ", init_mode=" << TEST_HA_INIT_MODE << "\n";

  // --- symbolic model -------------------------------------------------------
  FFGraph DAG;
  FFVar t  = DAG.add_var("t");
  FFVar x  = DAG.add_var("x");
  FFVar T  = DAG.add_var("T(t,x)");
  FFVar q  = DAG.add_var("q(t,x)");
  FFVar u  = DAG.add_var("u(t,x)");

  FFPartial OpP;

  // Model constants fold into the DAG as compile-time doubles (OCFE_PDEx idiom),
  // so no separate constant array needs threading through eval/deriv.
  double const av=p.a, cv=p.c, T0v=p.T0, Tsv=p.Ts, u0v=p.u0, xf=p.xf;

  // Block 0: conservative first-order heat system (parabolic).
  //FFVar HEAT_PDE  = OpP(T,t) - OpP(av*OpP(T,x),x);
  FFVar HEAT_PDE  = OpP(T,t) - OpP(q,x);
  FFVar HEAT_LINK = q - av*OpP(T,x);
  FFVar HEAT_INI  = T - Tsv - T0v*sin( kPi/2.0/xf * x );
  FFVar HEAT_BCL  = T - Tsv;
  FFVar HEAT_BCU  = q;

  // Block 1: linear advection (hyperbolic).
  FFVar ADV_PDE = OpP(u,t) + cv*OpP(u,x);
  FFVar ADV_INI = u - u0v*sin( 2.0*kPi/xf * x );
  FFVar ADV_BC  = u + u0v*sin( 2.0*kPi*cv/xf * t );

  OCFESLV oc(&DAG);
  oc.add_domain( t, FFDom(0., p.tf, TEST_HA_NEL_T, FFDom::LGR, TEST_HA_NT) );
  oc.add_domain( x, FFDom(0., p.xf, TEST_HA_NEL_X, FFDom::LGL, TEST_HA_NX) );
  oc.add_state( T, {t,x} );
  oc.add_state( q, {t,x} );
  oc.add_state( u, {t,x} );

  oc.update_ref( T, [&](OCFESLV::t_Coord const& cd){ return T_exact(cd.at(t),cd.at(x),p); } );
  oc.update_ref( q, [&](OCFESLV::t_Coord const& cd){ return q_exact(cd.at(t),cd.at(x),p); } );
  oc.update_ref( u, [&](OCFESLV::t_Coord const& cd){ return u_exact(cd.at(t),cd.at(x),p); } );

  oc.set_evolution_domain( t );

  // Block 0 (T,q) -> PARABOLIC.
  oc.add_equation( HEAT_PDE,  {t,x}, {FFDom::ALL-FFDom::LB, FFDom::ALL-FFDom::LB-FFDom::UB},
                   OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( HEAT_LINK, {t,x}, {FFDom::ALL,           FFDom::ALL},
                   OCFESLV::EqnOptions( OCFESLV::EqnRole::LINK, 0 ) );
  oc.add_equation( HEAT_INI,  {t,x}, {FFDom::LB,            FFDom::ALL},
                   OCFESLV::EqnOptions( OCFESLV::EqnRole::INITIAL, 0 ) );
  oc.add_equation( HEAT_BCL,  {t,x}, {FFDom::ALL-FFDom::LB, FFDom::LB},
                   OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );
  oc.add_equation( HEAT_BCU,  {t,x}, {FFDom::ALL-FFDom::LB, FFDom::UB},
                   OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );

  // Block 1 (u) -> HYPERBOLIC.
  oc.add_equation( ADV_PDE, {t,x}, {FFDom::ALL-FFDom::LB, FFDom::ALL-FFDom::LB},
                   OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 1 ) );
  oc.add_equation( ADV_INI, {t,x}, {FFDom::LB,            FFDom::ALL},
                   OCFESLV::EqnOptions( OCFESLV::EqnRole::INITIAL, 1 ) );
  oc.add_equation( ADV_BC,  {t,x}, {FFDom::ALL-FFDom::LB, FFDom::LB},
                   OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 1 ) );

  oc.options.INTERFACE.TYPE  = OCFESLV::Options::IC_AUTO;
  oc.options.INTERFACE.IMPOSITION = imp;
  // rev317: this driver closes its outflow faces itself, by extending its PDE masks to them.  The automatic
  // closure now appends only the per-face deficit, so those rows are no longer duplicated and the automatic
  // closure stays on (CRONOS_AUTO_HYP_CLOSURE, environment-only since rev318).  (Until rev317 it had to be switched off here, or each outflow face was closed twice.)
  oc.options.INTERFACE.SAT_SIGMA0      = TEST_HA_SAT_SIGMA0;
  // Classification (parabolic block 0 -> IC_VALUE, hyperbolic block 1 ->
  // IC_UPWIND) runs inside setup() when CLASSIFY is enabled; do NOT call
  // classify_pde() afterwards -- it deliberately clears the setup flag and would
  // require a second setup() before any node/eval query.
  oc.options.CLASSIFY.MODE = TEST_HA_CLASSIFY ? OCFESLV::Options::CLASS_AUTO
                                         : OCFESLV::Options::CLASS_NONE;
  // This driver hand is for monolithic only
  oc.options.SOLVE.MARCHING  = false;

  // item 11: the principal-symbol block-eligibility gate decides the continuity
  // drop / keep-explicit set entirely at setup() (IC_TRACE un-flag; IC_STRONG
  // setup-internal two-pass), so no driver-level re-derive loop is needed.  The
  // post-solve verify_interface_drop() below is a passive assertion only.
  size_t nVar=0, nEqn=0, nTrace=0;
  std::vector<double> xv;
  bool solved=false;
  {

  if( !oc.setup() ){
    std::cerr << "ERROR: OCFESLV::setup() failed for " << strimp << "\n";
    MR.ok = false; return MR;
  }
  std::cout << oc;

#ifdef MC__OCFESLV_ACOND_PROBE
  probe_saddle_conditioning( oc, strimp, 0 );
#endif

#if TEST_HA_CLASSIFY
  auto const& cls = oc.pde_type();
  MR.pde_type = OCFESLV::pde_type_name( cls.type );
  std::cout << "PDE type: " << MR.pde_type
            << "  At_singular=" << (cls.At_singular?"yes":"no")
            << "  parabolic=" << (cls.parabolic_structure_detected?"yes":"no")
            << "  evolution_hyperbolic=" << (cls.evolution_hyperbolic?"yes":"no") << "\n";

  // Per-block IC_AUTO dispatch dump.  The exact-coupling resolver folds IC_TRACE
  // into the "weak-like" branch, and IC_WEAK uses the SAT (weak) path, so only
  // IC_STRONG exposes IC_UPWIND for the hyperbolic block; the evolution (t)
  // interface is always causal IC_VALUE.  This reports the resolution actually
  // used by the current imposition mode.
  {
    auto iftype_name = []( OCFESLV::Options::InterfaceType it )->const char* {
      switch( it ){
        case OCFESLV::Options::IC_VALUE:  return "IC_VALUE";
        case OCFESLV::Options::IC_FLUX:   return "IC_FLUX";
        case OCFESLV::Options::IC_UPWIND: return "IC_UPWIND";
        case OCFESLV::Options::IC_AUTO:   return "IC_AUTO";
      }
      return "?";
    };
    bool const weak_path = ( imp == OCFESLV::Options::IC_WEAK );
    auto const& bclass = oc.block_classification();
    std::cout << "Per-block IC_AUTO dispatch (" << strimp << ", path="
              << (weak_path?"SAT/weak":"exact") << "):\n";
    if( bclass.empty() )
      std::cout << "  (no per-block classification map; resolving via aggregate type "
                << MR.pde_type << ")\n";
    int const block_ids[2] = { 0, 1 };       // block 0 = heat, block 1 = advection
    for( int bid : block_ids ){
      auto it = bclass.find( bid );
      OCFESLV::t_Classify const& bc = ( it != bclass.end() ) ? it->second : cls;
      std::cout << "  block " << bid << " (" << (bid==0?"heat/parabolic":"advection/hyperbolic")
                << "): type=" << OCFESLV::pde_type_name(bc.type)
                << " evolution_hyperbolic=" << (bc.evolution_hyperbolic?"yes":"no")
                << " parabolic=" << (bc.parabolic_structure_detected?"yes":"no")
                << " At_singular=" << (bc.At_singular?"yes":"no")
                << (it==bclass.end()?"  [aggregate]":"") << "\n";
      std::cout << "      t-interface (evolution) -> "
                << iftype_name( oc.resolved_interface_type( bid, t, OCFESLV::EqnRole::INTERIOR, weak_path ) )
                << "   |   x-interface (spatial) -> "
                << iftype_name( oc.resolved_interface_type( bid, x, OCFESLV::EqnRole::INTERIOR, weak_path ) )
                << "\n";
    }
  }
#endif

  nVar=oc.n_colloc_sta();
  nEqn=oc.n_colloc_eqn();
  nTrace=oc.n_colloc_trace();
  MR.nVar=nVar; MR.nEqn=nEqn; MR.nTrace=nTrace;
  std::cout << "nVar=" << nVar << " nEqn=" << nEqn << " nTrace=" << nTrace
            << " square=" << (nVar==nEqn?"yes":"NO") << "\n";
  ok &= (nVar==nEqn);

  // --- build the analytical reference vector (ordinary states; taus -> 0) ----
  std::vector<double> xExact; xExact.reserve(nVar);
  for( auto const& st : oc.states_colloc() ){
    std::string const nm = st.name();
    for( auto const& xy : oc.node_colloc(st) )
      xExact.push_back( exact_value_for_state(nm,xy,p) );
  }
  size_t const nState = (nVar>=nTrace) ? nVar-nTrace : nVar;
  if( xExact.size() != nState ){
    std::cerr << "ERROR: reference vector size mismatch: ordinary states=" << xExact.size()
              << " expected=" << nState << " (nVar=" << nVar << " nTrace=" << nTrace << ")\n";
    MR.ok = false; return MR;
  }
  if( nTrace ){
    std::cout << "Trace/tau variables appended: " << nTrace << " (initialised to zero)\n";
    xExact.resize(nVar,0.0);
  }

  std::vector<double> res(nEqn,0.0);
  if( !eval_residual(oc,xExact,res) ){
    std::cerr << "ERROR: reference residual evaluation failed\n"; MR.ok=false; return MR;
  }
  print_residuals("Reference-profile residual (truncation)", res);
  MR.ref_res = max_abs(res);
  print_variable_exact_errors( oc, xExact, p, "analytical reference" );
  max_duplicate_spread( oc, xExact, "analytical reference" );

  // --- initial guess --------------------------------------------------------
#if TEST_HA_INIT_MODE == 0
  xv = xExact;
  std::cout << "Initialisation: analytical reference profile\n";
#elif TEST_HA_INIT_MODE == 1
  xv = xExact;
  for( size_t i=0; i<nState; ++i )
    xv[i] += 1.0e-3*std::sin(0.37*double(i+1))*std::max(1.0,std::abs(xv[i]));
  std::cout << "Initialisation: analytical reference plus small perturbation\n";
#elif TEST_HA_INIT_MODE == 2
  xv.reserve(nVar);
  for( auto const& st : oc.states_colloc() ){
    std::string const nm=st.name();
    for( size_t i=0, N=oc.node_colloc(st).size(); i<N; ++i )
      xv.push_back( constant_initial_value_for_state(nm,p) );
  }
  xv.resize(nVar,0.0);
  std::cout << "Initialisation: constant primitive fields, zero trace variables\n";
#else
#error "TEST_HA_INIT_MODE must be 0, 1, or 2"
#endif

  if( !eval_residual(oc,xv,res) ){ MR.ok=false; return MR; }
  print_residuals("Initial residual", res);

  // item 12: nonlinear solve via OCFESLV::solve() (equilibrated LM/Newton, sparse AD).
  oc.options.SOLVE.MAX_ITER = TEST_HA_MAXIT;
  oc.options.SOLVE.RES_TOL  = TEST_HA_SOLVE_TOL;
  OCFESLV::SolveReport const srep = oc.solve( xv.data() );
  solved = srep.converged;
  if( !solved )
    std::cerr << "OCFESLV::solve did not converge: final|r|=" << srep.final_residual
              << " after " << srep.iterations << " it ("
              << srep.newton_steps << " Newton / " << srep.lm_steps << " LM)\n";
  MR.solved = solved; ok &= solved;
  if( !eval_residual(oc,xv,res) ){ MR.ok=false; return MR; }
  print_residuals("Final residual", res);
  MR.final_res = max_abs(res);

  StateExactErrors const E = print_variable_exact_errors( oc, xv, p, "final "+strimp+" solution" );
  MR.eT=E.eT; MR.eq=E.eq; MR.eu=E.eu; MR.emax=E.emax;
  MR.max_spread = max_duplicate_spread( oc, xv, "final "+strimp+" solution" );

  // Passive post-solve continuity assertion: item 11 fixes the drop/keep-explicit
  // set at setup(), so verify must not trigger a re-derive.  An OVERDROP here means
  // a mis-classified block and is a hard FAIL.
  if( oc.verify_interface_drop( xv.data() ) == OCFESLV::InterfaceDropStatus::OVERDROP ){
    std::cerr << "** PDE6: interface-drop verification FAILED (over-drop)\n";
    ok = false;
  }
  }   // end single-pass build/solve block

  // --- derivative (Jacobian) checks -----------------------------------------
#if TEST_HA_DERIV_CHECK
  bool deriv_ok = true;
  mc_test::DerivCheckOptions deriv_opt; deriv_opt.max_columns = 48;
  deriv_ok &= mc_test::check_oc_derivatives( oc, strimp+" @ solution", xv, nullptr, nullptr, deriv_opt );

  // Hyperbolic-interface probe: inject an element-local jump into u at the
  // shared spatial node x=xf/2 (perturb only every second duplicate occurrence,
  // plus everything to its right) and re-check the Jacobian at that state.
  std::vector<double> xPert = xv;
  size_t off=0;
  double const xmid = 0.5*p.xf;
  for( auto const& st : oc.states_colloc() ){
    auto nodes = oc.node_colloc(st);
    if( !st.name().empty() && st.name()[0]=='u' ){
      size_t dup=0;
      for( size_t i=0; i<nodes.size(); ++i ){
        double const xc = nodes[i][1];
        bool perturb = ( xc > xmid + 1e-12 );
        if( std::abs(xc-xmid) < 1e-12 ){ ++dup; if( dup%2==0 ) perturb=true; }
        if( perturb ) xPert[off+i] += 1e-2;
      }
    }
    off += nodes.size();
  }
  deriv_ok &= mc_test::check_oc_derivatives( oc, strimp+" @ local-u-jump state", xPert, nullptr, nullptr, deriv_opt );
  MR.deriv_ok = deriv_ok;
  ok &= deriv_ok;
#else
  MR.deriv_ok = true;
#endif

  // --- per-mode gates -------------------------------------------------------
  std::cout << "\nGates (" << strimp << "):\n";
  ok &= check_close("square system (nVar-nEqn)", std::abs(double(nVar)-double(nEqn)), 0.0);
  ok &= check_close("solve converged (final |r|)", MR.final_res, TEST_HA_SOLVE_TOL );
  double const dup_tol = std::max( double(TEST_HA_DUP_ABS_FLOOR),
                                   double(TEST_HA_DUP_TRUNC_FACTOR) * MR.ref_res );
  ok &= check_close("element-interface duplicate spread (vs truncation floor)",
                    MR.max_spread, dup_tol );
  ok &= check_close("max |T-T_exact| (spectral)", MR.eT, TEST_HA_EXACT_TOL );
  ok &= check_close("max |q-q_exact| (spectral)", MR.eq, TEST_HA_EXACT_TOL );
  ok &= check_close("max |u-u_exact| (spectral)", MR.eu, TEST_HA_EXACT_TOL );
#if TEST_HA_DERIV_CHECK
  ok &= check_close("Jacobian finite-difference checks", MR.deriv_ok?0.0:1.0, 0.0 );
#endif

  // --- solution / plot file output (matches PDE2/3/4 convention) ------------
  // One file per mode: gnuplot pm3d-friendly (blank line between t-blocks),
  // columns: t x T q u T_exact q_exact u_exact.
  {
    std::string suffix = strimp;
    auto const pos = suffix.find("IC_");
    if( pos != std::string::npos ) suffix = suffix.substr(pos+3);
    for( auto& ch : suffix ) ch = char(std::tolower((unsigned char)ch));
    std::string const prefix = std::string(TEST_HA_OUTPUT_PREFIX) + "_" + suffix;

    auto tNodes = physical_nodes( oc.var_domain().at(t) );
    auto xNodes = physical_nodes( oc.var_domain().at(x) );
    std::ofstream out( prefix + std::string(".out") );
    out << "# t x T q u T_exact q_exact u_exact\n";
    for( double tt : tNodes ){
      for( double xx : xNodes ){
        std::map<FFVar,double,lt_FFVar> pt{{t,tt},{x,xx}};
        out << std::setprecision(16)
            << tt << " " << xx << " "
            << interp(oc,T,pt,xv) << " "
            << interp(oc,q,pt,xv) << " "
            << interp(oc,u,pt,xv) << " "
            << T_exact(tt,xx,p) << " "
            << q_exact(tt,xx,p) << " "
            << u_exact(tt,xx,p) << "\n";
      }
      out << "\n";
    }
    std::cout << "Wrote plot file: " << prefix << ".out\n";
  }

  MR.ok = ok;
  std::cout << "Two-block heat+advection test (" << strimp << "): " << (ok?"PASS":"FAIL") << "\n";
  return MR;
}

int main()
{
  std::vector<ModeResult> results;
  results.push_back( run_mode( OCFESLV::Options::IC_WEAK,   "IC_WEAK"   ) );
  results.push_back( run_mode( OCFESLV::Options::IC_TRACE,  "IC_TRACE"  ) );
  results.push_back( run_mode( OCFESLV::Options::IC_STRONG, "IC_STRONG" ) );

  std::cout << "\n============== PDE6 two-block (parabolic+hyperbolic) sweep ==============\n";
  std::cout << std::left  << std::setw(11) << "mode"
            << std::right << std::setw(8)  << "nTrace"
            << std::setw(13) << "ref|r|"
            << std::setw(13) << "final|r|"
            << std::setw(12) << "dup-spread"
            << std::setw(11) << "|dT|"
            << std::setw(11) << "|dq|"
            << std::setw(11) << "|du|"
            << std::setw(8)  << "deriv"
            << std::setw(9)  << "result" << "\n";
  bool all_ok = true;
  for( auto const& r : results ){
    all_ok &= r.ok;
    std::cout << std::left  << std::setw(11) << r.name
              << std::right << std::scientific << std::setprecision(3)
              << std::setw(8)  << r.nTrace
              << std::setw(13) << r.ref_res
              << std::setw(13) << r.final_res
              << std::setw(12) << r.max_spread
              << std::setw(11) << r.eT
              << std::setw(11) << r.eq
              << std::setw(11) << r.eu
              << std::setw(8)  << (r.deriv_ok?"ok":"BAD")
              << std::setw(9)  << (r.ok?"PASS":"FAIL") << "\n";
  }
  std::cout << "========================================================================\n";
  std::cout << "Two-block parabolic/hyperbolic IC_AUTO sweep over IC_WEAK, IC_TRACE, IC_STRONG: "
            << (all_ok?"PASS":"FAIL") << "\n";
  return all_ok ? 0 : 1;
}
