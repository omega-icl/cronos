// Copyright (C) Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the EPL-2.0 with GPL-2.0-or-later as a Secondary License; see the LICENSE file.
// SPDX-License-Identifier: EPL-2.0 OR GPL-2.0-or-later
//
// OCFE_PDE7_solve2.cpp
// --------------------
// A genuinely characteristic 2x2 hyperbolic system, refactored from the legacy
// test7.cpp into the OCFE_PDEx_solve2.cpp driver format (config macros, Par
// struct, exact-profile functions, the generic LM/Newton harness, a per-mode
// run_mode() and an IC_WEAK/IC_TRACE/IC_STRONG acceptance sweep + plot output).
//
//   Symmetric wave / acoustic system (single hyperbolic block 0):
//     p_t + v_x = 0,
//     v_t + p_x = 0
//   principal matrices  A_t = I,  A_x = [ 0 1 ; 1 0 ]  (eigenvalues +1, -1).
//
// Because A_x has BOTH a +1 and a -1 eigenvalue, this is genuinely bidirectional:
// a true upwind SAT must split A_x into incoming/outgoing characteristic parts
//     A_x^+ = 1/2 [ 1 1 ; 1 1 ],   A_x^- = 1/2 [ -1 1 ; 1 -1 ],
// and contributes to BOTH PDE rows on BOTH sides of every internal x-interface.
// This is the distinguishing feature relative to OCFE_PDE6's scalar advection
// block (A_x = c, one-directional), and the reason the model is a much harder
// stress test of the characteristic interface machinery.
//
// Analytic (manufactured) solution -- NOT polynomial, so the reference-profile
// residual is a spectral-truncation level and the exact-match tolerances are
// spectral, not machine:
//     p(t,x) = sin(2*pi*(x - t)) + cos(2*pi*(x + t)),
//     v(t,x) = sin(2*pi*(x - t)) - cos(2*pi*(x + t)).
//
// Characteristic inflow data (one incoming characteristic per physical boundary):
//     x=0  : p + v = -2 sin(2*pi*t)        (right-going characteristic enters),
//     x=xf : p - v =  2 cos(2*pi*(xf + t)) (left-going characteristic enters).
//
// WELL-POSEDNESS CAVEAT (read before interpreting the sweep)
// ----------------------------------------------------------
// The legacy test7.cpp was deliberately an IC_WEAK + IC_UPWIND *residual* test:
// it never solved and never asserted squareness, dropping BOTH PDE rows at BOTH
// physical x-faces and closing the faces with the two scalar characteristic BCs
// alone.  That collocation is non-square under a real solve -- the outgoing
// characteristic at each x-face is left unpinned, so the system is short by
// (2 faces x #t-collocation nodes) rows (deficit 30 on the 2x3 / 8x8 grid).
//
// This driver therefore uses the asymmetric characteristic closure throughout:
// keep ONE scalar PDE row per physical x-face so the BC pins the incoming
// characteristic and a PDE row pins the outgoing one.  With A_x=[[0,1],[1,0]] the
// right-going characteristic is w+=p+v and the left-going is w-=p-v, so
//   x=0  : incoming w+ (BCL),  outgoing w- closed by PDEp  (PDEp keeps x=LB)
//   x=xf : incoming w- (BCR),  outgoing w+ closed by PDEv  (PDEv keeps x=UB).
// Since PDEp=(char+ + char-)/2 and PDEv=(char+ - char-)/2, each PDE row combined
// with its face BC closes the outgoing characteristic.  The layout is square
// (nEqn == nVar) by construction and is shared by all three imposition modes:
// IC_WEAK + IC_UPWIND is the legacy-validated characteristic-SAT reference, while
// IC_TRACE / IC_STRONG exercise the same closure under exact (tau/Schur) interface
// imposition.
//
// Build knobs (all overridable with -D):
//   -DTEST_W_NEL_T  axial(time) finite elements   (default 2)
//   -DTEST_W_NEL_X  spatial finite elements        (default 3)
//   -DTEST_W_NT     nodes per time element  (LGR)  (default 8)
//   -DTEST_W_NX     nodes per space element (LGL)  (default 8)
//   -DTEST_W_INIT_MODE  0 exact, 1 perturbed exact, 2 zero fields (default 1)
//   -DTEST_W_MAXIT / -DTEST_W_SOLVE_TOL / -DTEST_W_SAT_SIGMA0
//   -DTEST_W_CLASSIFY    run/report classify_pde()  (default 1)
//   -DTEST_W_DERIV_CHECK run the FD Jacobian checks  (default 1)

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

#ifndef TEST_W_NEL_T
#define TEST_W_NEL_T 2
#endif
#ifndef TEST_W_NEL_X
#define TEST_W_NEL_X 3
#endif
#ifndef TEST_W_NT
#define TEST_W_NT 8
#endif
#ifndef TEST_W_NX
#define TEST_W_NX 8
#endif
#ifndef TEST_W_INIT_MODE
#define TEST_W_INIT_MODE 1
#endif
#ifndef TEST_W_MAXIT
#define TEST_W_MAXIT 40
#endif
#ifndef TEST_W_SOLVE_TOL
#define TEST_W_SOLVE_TOL 1e-9
#endif
#ifndef TEST_W_SAT_SIGMA0
// Characteristic interface SAT penalty multiplier.  The baseline (1.0) under-
// penalizes the weak interface continuity for this bidirectional upwind system,
// leaving IC_WEAK loose (interface dup-spread ~8.4e-3, |p-exact|~2.2e-2).  10.0
// tightens it (dup-spread ~4.6e-5, |p-exact|~1.4e-3 -> PASS) with no effect on
// the exact modes; consistent with PDE6's TEST_HA_SAT_SIGMA0=10.0.
#define TEST_W_SAT_SIGMA0 10.0
#endif
#ifndef TEST_W_CLASSIFY
#define TEST_W_CLASSIFY 1
#endif
#ifndef TEST_W_DERIV_CHECK
#define TEST_W_DERIV_CHECK 1
#endif
#ifndef TEST_W_EXACT_TOL
#define TEST_W_EXACT_TOL 5e-3
#endif
// Element-interface duplicate-spread gate (refinement-based) -- see the PDE6
// driver for the rationale.  Explicit continuity is machine-tight; implied
// (soundly dropped) or weak (IC_WEAK SAT) continuity holds only to truncation,
// which MR.ref_res measures, so gate against max(absolute floor, factor*ref_res).
#ifndef TEST_W_DUP_ABS_FLOOR
#define TEST_W_DUP_ABS_FLOOR 1e-6
#endif
#ifndef TEST_W_DUP_TRUNC_FACTOR
#define TEST_W_DUP_TRUNC_FACTOR 3.0
#endif
#ifndef TEST_W_OUTPUT_PREFIX
#define TEST_W_OUTPUT_PREFIX "OCFE_PDE7_solve"
#endif

#define TEST_LAP_SPQR

using namespace mc;

static constexpr double kPi = 3.14159265358979323846264338327950288;

// ---------------------------------------------------------------------------
// Model parameters and analytical solution
// ---------------------------------------------------------------------------
struct Par
{
  double xf = 1.0;   // spatial extent
  double tf = 1.0;   // time horizon
};

static double p_exact( double t, double x, Par const& )
{
  return std::sin( 2.0*kPi*( x - t ) ) + std::cos( 2.0*kPi*( x + t ) );
}
static double v_exact( double t, double x, Par const& )
{
  return std::sin( 2.0*kPi*( x - t ) ) - std::cos( 2.0*kPi*( x + t ) );
}

// Reference value for a collocation state, keyed by state name and node coords.
// Node coordinate order matches the {t,x} domain order used in add_state below.
static double exact_value_for_state( std::string const& nm, std::vector<double> const& xy, Par const& p )
{
  double const t = xy[0], x = xy[1];
  switch( nm.empty() ? '?' : nm[0] ){
    case 'p': return p_exact( t, x, p );
    case 'v': return v_exact( t, x, p );
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
struct StateExactErrors { double ep=0.0, ev=0.0, emax=0.0; };

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
      case 'p': E.ep = emax; break;
      case 'v': E.ev = emax; break;
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
  double      ep=0.0, ev=0.0, emax=0.0;
  std::string pde_type  = "n/a";
};

// Physical collocation-node coordinates of a domain (element nodes concatenated;
// element-interface nodes appear once per adjoining element, as in PDE2/3/4/6).
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

  std::cout << "\n========== 2x2 characteristic wave system (hyperbolic) ==========\n";
  std::cout << "imposition: " << strimp << " + IC_UPWIND"
            << ", finite elements: t=" << TEST_W_NEL_T << " x=" << TEST_W_NEL_X
            << ", nodes/element: t=" << TEST_W_NT << " x=" << TEST_W_NX
            << ", init_mode=" << TEST_W_INIT_MODE << "\n";

  // --- symbolic model -------------------------------------------------------
  FFGraph DAG;
  FFVar t  = DAG.add_var("t");
  FFVar x  = DAG.add_var("x");
  FFVar p_ = DAG.add_var("p(t,x)");
  FFVar v_ = DAG.add_var("v(t,x)");

  FFPartial OpP;
  double const xf = p.xf;
  double const w  = 2.0*kPi;           // angular wavenumber (one period over [0,1])

  // Single hyperbolic block 0: symmetric 2x2 wave system.
  FFVar WAVE_PDEp = OpP(p_,t) + OpP(v_,x);
  FFVar WAVE_PDEv = OpP(v_,t) + OpP(p_,x);
  FFVar WAVE_INIp = p_ - ( sin( w*x ) + cos( w*x ) );
  FFVar WAVE_INIv = v_ - ( sin( w*x ) - cos( w*x ) );
  // Incoming right-going characteristic at x=0:  p + v = -2 sin(2*pi*t).
  FFVar WAVE_BCL  = p_ + v_ + 2.0*sin( w*t );
  // Incoming left-going characteristic at x=xf:  p - v =  2 cos(2*pi*(xf+t)).
  FFVar WAVE_BCR  = p_ - v_ - 2.0*cos( w*( xf + t ) );

  OCFESLV oc(&DAG);
  oc.add_domain( t, FFDom(0., p.tf, TEST_W_NEL_T, FFDom::LGR, TEST_W_NT) );
  oc.add_domain( x, FFDom(0., p.xf, TEST_W_NEL_X, FFDom::LGL, TEST_W_NX) );
  oc.add_state( p_, {t,x} );
  oc.add_state( v_, {t,x} );

  oc.update_ref( p_, [&](OCFESLV::t_Coord const& cd){ return p_exact(cd.at(t),cd.at(x),p); } );
  oc.update_ref( v_, [&](OCFESLV::t_Coord const& cd){ return v_exact(cd.at(t),cd.at(x),p); } );

  oc.set_evolution_domain( t );

  // Both PDE rows are in the same hyperbolic block 0.  PDE rows exclude the
  // initial face (t=LB).  Asymmetric characteristic closure on the spatial faces:
  // keep ONE scalar PDE row per physical x-face so the OUTGOING characteristic
  // there is pinned by a PDE row while the BC pins the incoming one.  With
  // A_x=[[0,1],[1,0]] the right-going characteristic is w+=p+v, the left-going is
  // w-=p-v:
  //   x=0  : incoming w+ (BCL), outgoing w- -> keep PDEp at x=LB (drop x=UB)
  //   x=xf : incoming w- (BCR), outgoing w+ -> keep PDEv at x=UB (drop x=LB)
  // PDEp=(char+ + char-)/2 with char+ pinned by BCL closes char-; PDEv=(char+ -
  // char-)/2 with char- pinned by BCR closes char+.  Square by construction
  // (nEqn == nVar) for all three imposition modes.
  oc.add_equation( WAVE_PDEp, {t,x}, {FFDom::ALL-FFDom::LB, FFDom::ALL-FFDom::UB},
                   OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( WAVE_PDEv, {t,x}, {FFDom::ALL-FFDom::LB, FFDom::ALL-FFDom::LB},
                   OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( WAVE_INIp, {t,x}, {FFDom::LB,            FFDom::ALL},
                   OCFESLV::EqnOptions( OCFESLV::EqnRole::INITIAL, 0 ) );
  oc.add_equation( WAVE_INIv, {t,x}, {FFDom::LB,            FFDom::ALL},
                   OCFESLV::EqnOptions( OCFESLV::EqnRole::INITIAL, 0 ) );
  oc.add_equation( WAVE_BCL,  {t,x}, {FFDom::ALL-FFDom::LB, FFDom::LB},
                   OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );
  oc.add_equation( WAVE_BCR,  {t,x}, {FFDom::ALL-FFDom::LB, FFDom::UB},
                   OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );

  // First-order system already -- no order reduction.  IC_UPWIND requests the
  // genuine characteristic A^+/A^- split SAT (the point of this model).
  oc.options.REDUCE.ORDER    = OCFESLV::Options::RED_NONE;
  oc.options.INTERFACE.TYPE  = OCFESLV::Options::IC_UPWIND;
  oc.options.INTERFACE.IMPOSITION = imp;
  // rev317: this driver closes its outflow faces itself, by extending its PDE masks to them.  The automatic
  // closure now appends only the per-face deficit, so those rows are no longer duplicated and the automatic
  // closure stays on (CRONOS_AUTO_HYP_CLOSURE, environment-only since rev318).  (Until rev317 it had to be switched off here, or each outflow face was closed twice.)
  oc.options.INTERFACE.SAT_SIGMA0      = TEST_W_SAT_SIGMA0;
  oc.options.CLASSIFY.MODE = TEST_W_CLASSIFY ? OCFESLV::Options::CLASS_AUTO
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

#if TEST_W_CLASSIFY
  auto const& cls = oc.pde_type();
  MR.pde_type = OCFESLV::pde_type_name( cls.type );
  std::cout << "PDE type: " << MR.pde_type
            << "  At_singular=" << (cls.At_singular?"yes":"no")
            << "  parabolic=" << (cls.parabolic_structure_detected?"yes":"no")
            << "  evolution_hyperbolic=" << (cls.evolution_hyperbolic?"yes":"no") << "\n";

  // Single-block IC dispatch report.  For IC_WEAK the SAT (weak) path is used
  // with the characteristic A^+/A^- split; IC_TRACE folds the spatial interface
  // into IC_VALUE; IC_STRONG exposes the IC_UPWIND characteristic exact rows.
  // The evolution (t) interface is always causal IC_VALUE.
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
    auto it0 = bclass.find( 0 );
    OCFESLV::t_Classify const& bc = ( it0 != bclass.end() ) ? it0->second : cls;
    std::cout << "Block 0 (wave/hyperbolic) IC dispatch (" << strimp << ", path="
              << (weak_path?"SAT/weak":"exact") << "):\n";
    std::cout << "  type=" << OCFESLV::pde_type_name(bc.type)
              << " evolution_hyperbolic=" << (bc.evolution_hyperbolic?"yes":"no")
              << (it0==bclass.end()?"  [aggregate]":"") << "\n";
    std::cout << "      t-interface (evolution) -> "
              << iftype_name( oc.resolved_interface_type( 0, t, OCFESLV::EqnRole::INTERIOR, weak_path ) )
              << "   |   x-interface (spatial) -> "
              << iftype_name( oc.resolved_interface_type( 0, x, OCFESLV::EqnRole::INTERIOR, weak_path ) )
              << "\n";
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
#if TEST_W_INIT_MODE == 0
  xv = xExact;
  std::cout << "Initialisation: analytical reference profile\n";
#elif TEST_W_INIT_MODE == 1
  xv = xExact;
  for( size_t i=0; i<nState; ++i )
    xv[i] += 1.0e-3*std::sin(0.37*double(i+1))*std::max(1.0,std::abs(xv[i]));
  std::cout << "Initialisation: analytical reference plus small perturbation\n";
#elif TEST_W_INIT_MODE == 2
  xv.assign(nVar,0.0);
  std::cout << "Initialisation: zero fields, zero trace variables\n";
#else
#error "TEST_W_INIT_MODE must be 0, 1, or 2"
#endif

  if( !eval_residual(oc,xv,res) ){ MR.ok=false; return MR; }
  print_residuals("Initial residual", res);

  // item 12: nonlinear solve via OCFESLV::solve() (equilibrated LM/Newton, sparse AD).
  oc.options.SOLVE.MAX_ITER = TEST_W_MAXIT;
  oc.options.SOLVE.RES_TOL  = TEST_W_SOLVE_TOL;
#if defined(CRONOS__WITH_SPQR) && defined(TEST_LAP_SPQR)
  // item 13: route the augmented (IC_TRACE/IC_STRONG) solve through sparse
  // rank-revealing QR instead of the JtJ normal equations.  Requires the header
  // built with -DCRONOS__WITH_SPQR and SuiteSparse linked.
  oc.options.SOLVE.FACTORIZATION = OCFESLV::Options::SOLVE_SPQR;
#endif

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
  MR.ep=E.ep; MR.ev=E.ev; MR.emax=E.emax;
  MR.max_spread = max_duplicate_spread( oc, xv, "final "+strimp+" solution" );

  // Passive post-solve continuity assertion: item 11 fixes the drop/keep-explicit
  // set at setup(), so verify must not trigger a re-derive.  An OVERDROP here means
  // a mis-classified block and is a hard FAIL.
  if( oc.verify_interface_drop( xv.data() ) == OCFESLV::InterfaceDropStatus::OVERDROP ){
    std::cerr << "** PDE7: interface-drop verification FAILED (over-drop)\n";
    ok = false;
  }
  }   // end single-pass build/solve block

  // --- derivative (Jacobian) checks -----------------------------------------
#if TEST_W_DERIV_CHECK
  bool deriv_ok = true;
  mc_test::DerivCheckOptions deriv_opt; deriv_opt.max_columns = 48;
  deriv_ok &= mc_test::check_oc_derivatives( oc, strimp+" @ solution", xv, nullptr, nullptr, deriv_opt );

  // Characteristic-SAT probe (legacy test7 signature): a pure p-jump across the
  // x-element interfaces.  The perturbation is constant within each x-element
  // (zero local x-derivative), so its only effect is the jump at every x
  // interface -- which the true A^+/A^- characteristic SAT must feel on BOTH
  // PDE rows and BOTH sides of each interface.  Re-check the Jacobian there.
  std::vector<double> xPert = xv;
  size_t off=0;
  double const elw = p.xf / double(TEST_W_NEL_X);
  for( auto const& st : oc.states_colloc() ){
    auto nodes = oc.node_colloc(st);
    if( !st.name().empty() && st.name()[0]=='p' ){
      for( size_t i=0; i<nodes.size(); ++i ){
        double const xc = nodes[i][1];
        int iel = int( std::floor( ( xc + 1e-9 ) / elw ) );
        if( iel < 0 ) iel = 0;
        if( iel >= TEST_W_NEL_X ) iel = TEST_W_NEL_X-1;
        xPert[off+i] += 1.0e-2 * double(iel);
      }
    }
    off += nodes.size();
  }
  deriv_ok &= mc_test::check_oc_derivatives( oc, strimp+" @ p-jump state", xPert, nullptr, nullptr, deriv_opt );
  MR.deriv_ok = deriv_ok;
  ok &= deriv_ok;
#else
  MR.deriv_ok = true;
#endif

  // --- per-mode gates -------------------------------------------------------
  std::cout << "\nGates (" << strimp << "):\n";
  ok &= check_close("square system (nVar-nEqn)", std::abs(double(nVar)-double(nEqn)), 0.0);
  ok &= check_close("solve converged (final |r|)", MR.final_res, TEST_W_SOLVE_TOL );
  double const dup_tol = std::max( double(TEST_W_DUP_ABS_FLOOR),
                                   double(TEST_W_DUP_TRUNC_FACTOR) * MR.ref_res );
  ok &= check_close("element-interface duplicate spread (vs truncation floor)",
                    MR.max_spread, dup_tol );
  ok &= check_close("max |p-p_exact| (spectral)", MR.ep, TEST_W_EXACT_TOL );
  ok &= check_close("max |v-v_exact| (spectral)", MR.ev, TEST_W_EXACT_TOL );
#if TEST_W_DERIV_CHECK
  ok &= check_close("Jacobian finite-difference checks", MR.deriv_ok?0.0:1.0, 0.0 );
#endif

  // --- solution / plot file output (matches PDE2/3/4/6 convention) ----------
  // One file per mode: gnuplot pm3d-friendly (blank line between t-blocks),
  // columns: t x p v p_exact v_exact.
  {
    std::string suffix = strimp;
    auto const pos = suffix.find("IC_");
    if( pos != std::string::npos ) suffix = suffix.substr(pos+3);
    for( auto& ch : suffix ) ch = char(std::tolower((unsigned char)ch));
    std::string const prefix = std::string(TEST_W_OUTPUT_PREFIX) + "_" + suffix;

    auto tNodes = physical_nodes( oc.var_domain().at(t) );
    auto xNodes = physical_nodes( oc.var_domain().at(x) );
    std::ofstream out( prefix + std::string(".out") );
    out << "# t x p v p_exact v_exact\n";
    for( double tt : tNodes ){
      for( double xx : xNodes ){
        std::map<FFVar,double,lt_FFVar> pt{{t,tt},{x,xx}};
        out << std::setprecision(16)
            << tt << " " << xx << " "
            << interp(oc,p_,pt,xv) << " "
            << interp(oc,v_,pt,xv) << " "
            << p_exact(tt,xx,p) << " "
            << v_exact(tt,xx,p) << "\n";
      }
      out << "\n";
    }
    std::cout << "Wrote plot file: " << prefix << ".out\n";
  }

  MR.ok = ok;
  std::cout << "2x2 wave system test (" << strimp << "): " << (ok?"PASS":"FAIL") << "\n";
  return MR;
}

int main()
{
  std::vector<ModeResult> results;
  results.push_back( run_mode( OCFESLV::Options::IC_WEAK,   "IC_WEAK"   ) );
  results.push_back( run_mode( OCFESLV::Options::IC_TRACE,  "IC_TRACE"  ) );
  results.push_back( run_mode( OCFESLV::Options::IC_STRONG, "IC_STRONG" ) );

  std::cout << "\n============== PDE7 2x2 characteristic wave sweep ==============\n";
  std::cout << std::left  << std::setw(11) << "mode"
            << std::right << std::setw(8)  << "nTrace"
            << std::setw(13) << "ref|r|"
            << std::setw(13) << "final|r|"
            << std::setw(12) << "dup-spread"
            << std::setw(11) << "|dp|"
            << std::setw(11) << "|dv|"
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
              << std::setw(11) << r.ep
              << std::setw(11) << r.ev
              << std::setw(8)  << (r.deriv_ok?"ok":"BAD")
              << std::setw(9)  << (r.ok?"PASS":"FAIL") << "\n";
  }
  std::cout << "================================================================\n";
  std::cout << "2x2 characteristic-wave IC_UPWIND sweep over IC_WEAK, IC_TRACE, IC_STRONG: "
            << (all_ok?"PASS":"FAIL") << "\n";
  // IC_WEAK + IC_UPWIND is the legacy-validated configuration; a FAIL in
  // IC_TRACE/IC_STRONG most likely reflects the characteristic-closure
  // well-posedness caveat documented in the header, not a solver defect.
  return all_ok ? 0 : 1;
}
