// ===========================================================================
// OCFE_PDE10_solve2.cpp
//
// Three-direction (t,x,y) scalar-hyperbolic test -- the 3D corner test under
// genuine UPWIND.  A single first-order advection block, mixed orientation:
//
//   d_t u + c d_x u + d d_y u = 0,   c = +0.6 (inflow x=0),  d = -0.4 (inflow y=L)
//
// PURPOSE.  The hyperbolic sibling of the PDE5/PDE9 reduced-aux corner work.
// Both spatial directions are first-order (no second derivative, no reduced
// auxiliary), so _resolve_interface_type takes the UPWIND branch for x and y.
// This exercises the *primitive* tensor-corner dedup under genuine UPWIND at 3D
// edges/corners -- the exact-row/donor corner path, a SIBLING of the reduced-aux
// receiver-graph suppression (the reduced-aux corner probe is inert here: with
// no aux links, reduced_flux_state_continuity is false everywhere).
//
// A SINGLE block suffices: the corner machinery is per-block and there is no
// cross-block coupling, so two decoupled blocks would only duplicate the test.
// Instead the one block uses MIXED orientation -- forward in x, backward in y --
// so the x-receiver and y-receiver at every 3D edge/corner come from opposite
// sides.  That is the adversarial case for the ordered-ownership corner dedup,
// which a same-orientation (uniform-inflow) block never produces.  The element-
// interface duplicate-spread gate below is the detector: if the primitive corner
// dedup fails to span a 3D edge/corner the spread blows past the truncation
// floor, exactly as the reduced auxes do in PDE5.
//
// The exact solution is a travelling plane wave u0*sin(k(x-c t)+k(y-d t)),
// which solves the 2D advection exactly for either sign of the speeds; inflow
// Dirichlet data on the two upwind faces (x=0, y=L), initial data on the whole
// spatial grid.
//
// Build flags (gate is default-ON in the header; the tag probe is opt-in):
//   -DMC__OCFESLV_STRONG_EXPLICIT_TAU -DMC__OCFESLV_INTERFACE_DECISION_V2
//
// Knobs (compile-time -D overrides):
//   -DTEST_SA_NEL_T   evolution(time) finite elements   (default 3)
//   -DTEST_SA_NEL_X   x finite elements                 (default 3)
//   -DTEST_SA_NEL_Y   y finite elements                 (default 3)
//   -DTEST_SA_NT      nodes per time element  (LGR)     (default 6)
//   -DTEST_SA_NX      nodes per x element     (LGL)     (default 6)
//   -DTEST_SA_NY      nodes per y element     (LGL)     (default 6)
// ===========================================================================

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
//#define MC__OCFESLV_SYMBOL_TAG_PROBE      // print per-block symbol fingerprint + per-claim flags

#include OCFE_OCFESLV_HEADER

#include "test_deriv_utils.hpp"

#ifndef TEST_SA_NEL_T
#define TEST_SA_NEL_T 3
#endif
#ifndef TEST_SA_NEL_X
#define TEST_SA_NEL_X 3
#endif
#ifndef TEST_SA_NEL_Y
#define TEST_SA_NEL_Y 3
#endif
#ifndef TEST_SA_NT
#define TEST_SA_NT 6
#endif
#ifndef TEST_SA_NX
#define TEST_SA_NX 6
#endif
#ifndef TEST_SA_NY
#define TEST_SA_NY 6
#endif
#ifndef TEST_SA_INIT_MODE
#define TEST_SA_INIT_MODE 1
#endif
#ifndef TEST_SA_MAXIT
#define TEST_SA_MAXIT 40
#endif
#ifndef TEST_SA_SOLVE_TOL
#define TEST_SA_SOLVE_TOL 1e-9
#endif
#ifndef TEST_SA_SAT_SIGMA0
#define TEST_SA_SAT_SIGMA0 10.0
#endif
#ifndef TEST_SA_CLASSIFY
#define TEST_SA_CLASSIFY 1
#endif
#ifndef TEST_SA_DERIV_CHECK
#define TEST_SA_DERIV_CHECK 1
#endif
#ifndef TEST_SA_EXACT_TOL
#define TEST_SA_EXACT_TOL 5e-3
#endif
// Element-interface duplicate-spread gate (refinement-based, same rationale as
// PDE6): explicit continuity is machine-tight; implied (dropped/redundant) or
// weak (SAT) continuity holds only to truncation, which MR.ref_res measures.
#ifndef TEST_SA_DUP_ABS_FLOOR
#define TEST_SA_DUP_ABS_FLOOR 1e-6
#endif
#ifndef TEST_SA_DUP_TRUNC_FACTOR
#define TEST_SA_DUP_TRUNC_FACTOR 3.0
#endif
#ifndef TEST_SA_OUTPUT_PREFIX
#define TEST_SA_OUTPUT_PREFIX "OCFE_PDE10_solve"
#endif

#define TEST_LAP_SPQR

using namespace mc;

static constexpr double kPi = 3.14159265358979323846264338327950288;

// ---------------------------------------------------------------------------
// Model parameters and analytical solutions
// ---------------------------------------------------------------------------
struct Par
{
  double c =  0.6;    // x-advection speed (>0, forward; inflow x=0)
  double d = -0.4;    // y-advection speed (<0, backward; inflow y=L)
                      // c+d != 0 so the plane wave genuinely evolves in time
                      // (c=-d would collapse the phase to k(x+y), steady).
  double u0 = 1.0;    // amplitude
  double xf = 1.0;    // x extent
  double yf = 1.0;    // y extent
  double tf = 1.0;    // time horizon
};

static double u_exact( double t, double x, double y, Par const& p )
{
  double const k = 2.0*kPi/p.xf;
  return p.u0 * std::sin( k*( x - p.c*t ) + k*( y - p.d*t ) );
}

// Reference value for a collocation state, keyed by state name and node coords.
// Node coordinate order matches the domain order used in add_state: (t,x,y).
static double exact_value_for_state( std::string const& nm, std::vector<double> const& xyz, Par const& p )
{
  double const t = xyz[0], x = xyz[1], y = xyz.size()>2 ? xyz[2] : 0.0;
  switch( nm.empty()?'?':nm[0] ){
    case 'u': return u_exact( t, x, y, p );
    default:  return 0.0;
  }
}
static double constant_initial_value_for_state( std::string const&, Par const& )
{ return 0.0; }

// ---------------------------------------------------------------------------
// Generic residual / Jacobian / LM-Newton harness (shared idiom with OCFE_PDEx)
// ---------------------------------------------------------------------------
static double max_abs( std::vector<double> const& r )
{ double m=0.; for( double v: r ) m=std::max(m,std::abs(v)); return m; }

static void print_residuals( std::string const& label, std::vector<double> const& r )
{
  double sum=0.; for( double v: r ) sum += std::abs(v);
  std::cout << std::left << std::setw(46) << label
            << std::right << " n=" << std::setw(5) << r.size()
            << "  max|r|=" << std::scientific << std::setprecision(4) << max_abs(r)
            << " mean|r|=" << (r.empty()?0.0:sum/double(r.size())) << "\n";
}

static bool check_close( std::string const& label, double val, double tol )
{
  bool ok = std::isfinite(val) && val <= tol;
  std::cout << std::left << std::setw(50) << label
            << std::right << " value=" << std::scientific << std::setprecision(6) << val
            << " tol=" << tol << "  " << (ok?"PASS":"FAIL") << "\n";
  return ok;
}

static bool eval_residual( OCFESLV& oc, std::vector<double> const& x, std::vector<double>& r )
{ std::fill(r.begin(),r.end(),0.0); return oc.eval(r.data(),nullptr,x.data(),nullptr,nullptr); }


// ---------------------------------------------------------------------------
// Comparison / spread helpers
// ---------------------------------------------------------------------------
struct StateExactErrors { double eu=0.0, ev=0.0, emax=0.0; };

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
      case 'u': E.eu = emax; break;
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
    std::map<std::vector<long long>, std::pair<double,double>> g;
    for( size_t i=0; i<nodes.size(); ++i ){
      std::vector<long long> key; key.reserve(nodes[i].size());
      for( double cc: nodes[i] ) key.push_back( (long long)std::llround(cc*1.0e12) );
      double const vv=var[off+i];
      auto it=g.find(key);
      if( it==g.end() ) g.emplace(std::move(key),std::make_pair(vv,vv));
      else { it->second.first=std::min(it->second.first,vv); it->second.second=std::max(it->second.second,vv); }
    }
    double sp=0.0; for( auto const& kv: g ) sp=std::max(sp,kv.second.second-kv.second.first);
    max_pair=std::max(max_pair,sp);
    std::cout << "  " << std::setw(10) << st.name()
              << " max_pair_spread=" << std::scientific << std::setprecision(6) << sp << "\n";
    off += nodes.size();
  }
  return max_pair;
}

static std::vector<double> physical_nodes( FFDom const& dom )
{
  std::vector<double> xs;
  for( size_t ie=0; ie<dom.n_elem; ++ie ){
    auto xe = dom.lgnodes( dom.elem_lo(ie), dom.elem_up(ie) );
    xs.insert( xs.end(), xe.begin(), xe.end() );
  }
  return xs;
}
static double interp( OCFESLV const& oc, FFVar const& v,
                      std::map<FFVar,double,lt_FFVar> const& pt,
                      std::vector<double> const& var )
{ return oc.eval_colloc<double>( v, pt, var.data(), nullptr, nullptr ); }

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
  double      eu=0.0, ev=0.0, emax=0.0;
  std::string pde_type  = "n/a";
};

// ---------------------------------------------------------------------------
// One imposition mode
// ---------------------------------------------------------------------------
static ModeResult run_mode( OCFESLV::Options::ImpositionType imp, std::string const& strimp )
{
  ModeResult MR; MR.name = strimp;
  bool ok = true;
  Par p;

  std::cout << "\n========== three-direction single-block scalar advection (mixed: fwd-x, bwd-y, hyperbolic) ==========\n";
  std::cout << "imposition: " << strimp
            << ", finite elements: t=" << TEST_SA_NEL_T << " x=" << TEST_SA_NEL_X
            << " y=" << TEST_SA_NEL_Y
            << ", nodes/element: t=" << TEST_SA_NT << " x=" << TEST_SA_NX
            << " y=" << TEST_SA_NY
            << ", init_mode=" << TEST_SA_INIT_MODE << "\n";

  // --- symbolic model -------------------------------------------------------
  FFGraph DAG;
  FFVar t  = DAG.add_var("t");
  FFVar x  = DAG.add_var("x");
  FFVar y  = DAG.add_var("y");
  FFVar u  = DAG.add_var("u(t,x,y)");

  FFPartial OpP;

  double const cv=p.c, dv=p.d, u0v=p.u0;
  double const k2 = 2.0*kPi/p.xf;

  // Travelling plane wave: exact solution of the 2D advection.  Keeping x,y
  // symbolic, it evaluates to the correct inflow/initial data on whichever face
  // it is collocated, so the one expression serves both INI and the inflow BCs.
  FFVar UE = u0v*sin( k2*( x - cv*t ) + k2*( y - dv*t ) );   // u_exact(t,x,y)

  // Single block: forward in x (c>0, inflow x=LB), backward in y (d<0, inflow
  // y=UB).  The mixed orientation makes the x- and y-receivers at every 3D
  // edge/corner come from OPPOSITE sides -- the adversarial case for the
  // tensor-corner ownership/spanning under UPWIND.
  FFVar ADV_PDE = OpP(u,t) + cv*OpP(u,x) + dv*OpP(u,y);
  FFVar ADV_REF = u - UE;

  OCFESLV oc(&DAG);
  oc.add_domain( t, FFDom(0., p.tf, TEST_SA_NEL_T, FFDom::LGR, TEST_SA_NT) );
  oc.add_domain( x, FFDom(0., p.xf, TEST_SA_NEL_X, FFDom::LGL, TEST_SA_NX) );
  oc.add_domain( y, FFDom(0., p.yf, TEST_SA_NEL_Y, FFDom::LGL, TEST_SA_NY) );
  oc.add_state( u, {t,x,y} );

  oc.update_ref( u, [&](OCFESLV::t_Coord const& cd){ return u_exact(cd.at(t),cd.at(x),cd.at(y),p); } );

  oc.set_evolution_domain( t );

  int const T_INT = FFDom::ALL - FFDom::LB;   // t>0 (post-initial)

  // Node partition for t>0:  x=LB (all y)   -> x-inflow BC;
  //                          x>LB, y=UB     -> y-inflow BC;
  //                          x>LB, y<UB     -> interior PDE (outflow faces
  //                          x=UB and y=LB carry the PDE, no BC).
  //                          t=LB           -> initial data (all x,y).
  oc.add_equation( ADV_PDE, {t,x,y}, {T_INT, FFDom::ALL-FFDom::LB, FFDom::ALL-FFDom::UB},
                   OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( ADV_REF, {t,x,y}, {FFDom::LB, FFDom::ALL, FFDom::ALL},
                   OCFESLV::EqnOptions( OCFESLV::EqnRole::INITIAL, 0 ) );
  oc.add_equation( ADV_REF, {t,x,y}, {T_INT, FFDom::LB, FFDom::ALL},
                   OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );   // x-inflow (x=LB)
  oc.add_equation( ADV_REF, {t,x,y}, {T_INT, FFDom::ALL-FFDom::LB, FFDom::UB},
                   OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );   // y-inflow (x>LB, y=UB)

  oc.options.INTERFACE.TYPE  = OCFESLV::Options::IC_AUTO;
  oc.options.INTERFACE.IMPOSITION = imp;
  // rev317: this driver closes its outflow faces itself, by extending its PDE masks to them.  The automatic
  // closure now appends only the per-face deficit, so those rows are no longer duplicated and the automatic
  // closure stays on (CRONOS_AUTO_HYP_CLOSURE, environment-only since rev318).  (Until rev317 it had to be switched off here, or each outflow face was closed twice.)
  oc.options.INTERFACE.SAT_SIGMA0      = TEST_SA_SAT_SIGMA0;
  oc.options.CLASSIFY.MODE = TEST_SA_CLASSIFY ? OCFESLV::Options::CLASS_AUTO
                                         : OCFESLV::Options::CLASS_NONE;

  // item 11: the principal-symbol block-eligibility gate decides the continuity
  // drop / keep-explicit set entirely at setup() (IC_TRACE un-flag; IC_STRONG
  // setup-internal two-pass).  For these scalar (eig=1) blocks the gate stays a
  // no-op; no driver-level re-derive loop is needed and the post-solve
  // verify_interface_drop() below is a passive assertion only.
  size_t nVar=0, nEqn=0, nTrace=0;
  std::vector<double> xv;
  bool solved=false;
  {

  if( !oc.setup() ){
    std::cerr << "ERROR: OCFESLV::setup() failed for " << strimp << "\n";
    MR.ok = false; return MR;
  }
  std::cout << oc;

#if TEST_SA_CLASSIFY
  auto const& cls = oc.pde_type();
  MR.pde_type = OCFESLV::pde_type_name( cls.type );
  std::cout << "PDE type: " << MR.pde_type
            << "  At_singular=" << (cls.At_singular?"yes":"no")
            << "  parabolic=" << (cls.parabolic_structure_detected?"yes":"no")
            << "  evolution_hyperbolic=" << (cls.evolution_hyperbolic?"yes":"no") << "\n";
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
    int const block_ids[1] = { 0 };
    for( int bid : block_ids ){
      auto it = bclass.find( bid );
      OCFESLV::t_Classify const& bc = ( it != bclass.end() ) ? it->second : cls;
      std::cout << "  block " << bid << " (mixed: fwd-x, bwd-y)"
                << ": type=" << OCFESLV::pde_type_name(bc.type)
                << " evolution_hyperbolic=" << (bc.evolution_hyperbolic?"yes":"no")
                << " parabolic=" << (bc.parabolic_structure_detected?"yes":"no")
                << " At_singular=" << (bc.At_singular?"yes":"no")
                << (it==bclass.end()?"  [aggregate]":"") << "\n";
      std::cout << "      t-interface (evolution) -> "
                << iftype_name( oc.resolved_interface_type( bid, t, OCFESLV::EqnRole::INTERIOR, weak_path ) )
                << "   |   x-interface -> "
                << iftype_name( oc.resolved_interface_type( bid, x, OCFESLV::EqnRole::INTERIOR, weak_path ) )
                << "   |   y-interface -> "
                << iftype_name( oc.resolved_interface_type( bid, y, OCFESLV::EqnRole::INTERIOR, weak_path ) )
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

  // --- analytical reference vector (ordinary states; taus -> 0) -------------
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
#if TEST_SA_INIT_MODE == 0
  xv = xExact;
  std::cout << "Initialisation: analytical reference profile\n";
#elif TEST_SA_INIT_MODE == 1
  xv = xExact;
  for( size_t i=0; i<nState; ++i )
    xv[i] += 1.0e-3*std::sin(0.37*double(i+1))*std::max(1.0,std::abs(xv[i]));
  std::cout << "Initialisation: analytical reference plus small perturbation\n";
#elif TEST_SA_INIT_MODE == 2
  xv.reserve(nVar);
  for( auto const& st : oc.states_colloc() ){
    std::string const nm=st.name();
    for( size_t i=0, N=oc.node_colloc(st).size(); i<N; ++i )
      xv.push_back( constant_initial_value_for_state(nm,p) );
  }
  xv.resize(nVar,0.0);
  std::cout << "Initialisation: constant fields, zero trace variables\n";
#else
#error "TEST_SA_INIT_MODE must be 0, 1, or 2"
#endif

  if( !eval_residual(oc,xv,res) ){ MR.ok=false; return MR; }
  print_residuals("Initial residual", res);

  // This driver hand is for monolithic only
  oc.options.SOLVE.MARCHING = false;
  // item 12: nonlinear solve via OCFESLV::solve() (equilibrated LM/Newton, sparse AD).
  oc.options.SOLVE.MAX_ITER = TEST_SA_MAXIT;
  oc.options.SOLVE.RES_TOL  = TEST_SA_SOLVE_TOL;
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
  MR.eu=E.eu; MR.emax=E.emax;
  MR.max_spread = max_duplicate_spread( oc, xv, "final "+strimp+" solution" );

  // Passive post-solve continuity assertion: item 11 fixes the drop/keep-explicit
  // set at setup(), so verify must not trigger a re-derive.  An OVERDROP here means
  // a mis-classified block and is a hard FAIL.
  if( oc.verify_interface_drop( xv.data() ) == OCFESLV::InterfaceDropStatus::OVERDROP ){
    std::cerr << "** PDE10: interface-drop verification FAILED (over-drop)\n";
    ok = false;
  }
  }   // end single-pass build/solve block

  // --- derivative (Jacobian) checks -----------------------------------------
#if TEST_SA_DERIV_CHECK
  bool deriv_ok = true;
  mc_test::DerivCheckOptions deriv_opt; deriv_opt.max_columns = 48;
  deriv_ok &= mc_test::check_oc_derivatives( oc, strimp+" @ solution", xv, nullptr, nullptr, deriv_opt );
  MR.deriv_ok = deriv_ok;
  ok &= deriv_ok;
#else
  MR.deriv_ok = true;
#endif

  // --- per-mode gates -------------------------------------------------------
  std::cout << "\nGates (" << strimp << "):\n";
  ok &= check_close("square system (nVar-nEqn)", std::abs(double(nVar)-double(nEqn)), 0.0);
  ok &= check_close("solve converged (final |r|)", MR.final_res, TEST_SA_SOLVE_TOL );
  double const dup_tol = std::max( double(TEST_SA_DUP_ABS_FLOOR),
                                   double(TEST_SA_DUP_TRUNC_FACTOR) * MR.ref_res );
  ok &= check_close("element-interface duplicate spread (vs truncation floor)",
                    MR.max_spread, dup_tol );
  ok &= check_close("max |u-u_exact| (spectral)", MR.eu, TEST_SA_EXACT_TOL );
#if TEST_SA_DERIV_CHECK
  ok &= check_close("Jacobian finite-difference checks", MR.deriv_ok?0.0:1.0, 0.0 );
#endif

  // --- plot file output -----------------------------------------------------
  {
    std::string suffix = strimp;
    auto const pos = suffix.find("IC_");
    if( pos != std::string::npos ) suffix = suffix.substr(pos+3);
    for( auto& ch : suffix ) ch = char(std::tolower((unsigned char)ch));
    std::string const prefix = std::string(TEST_SA_OUTPUT_PREFIX) + "_" + suffix;

    auto tNodes = physical_nodes( oc.var_domain().at(t) );
    auto xNodes = physical_nodes( oc.var_domain().at(x) );
    auto yNodes = physical_nodes( oc.var_domain().at(y) );
    std::ofstream out( prefix + std::string(".out") );
    out << "# t x y u u_exact\n";
    for( double tt : tNodes ){
      for( double xx : xNodes ){
        for( double yy : yNodes ){
          std::map<FFVar,double,lt_FFVar> pt{{t,tt},{x,xx},{y,yy}};
          out << std::setprecision(16)
              << tt << " " << xx << " " << yy << " "
              << interp(oc,u,pt,xv) << " "
              << u_exact(tt,xx,yy,p) << "\n";
        }
      }
      out << "\n";
    }
    std::cout << "Wrote plot file: " << prefix << ".out\n";
  }

  MR.ok = ok;
  std::cout << "Three-direction single-block scalar-advection test (" << strimp << "): " << (ok?"PASS":"FAIL") << "\n";
  return MR;
}

int main()
{
  std::vector<ModeResult> results;
  results.push_back( run_mode( OCFESLV::Options::IC_WEAK,   "IC_WEAK"   ) );
  results.push_back( run_mode( OCFESLV::Options::IC_TRACE,  "IC_TRACE"  ) );
  results.push_back( run_mode( OCFESLV::Options::IC_STRONG, "IC_STRONG" ) );

  std::cout << "\n============== PDE10 three-direction single-block scalar-advection (mixed fwd-x/bwd-y) sweep ==============\n";
  std::cout << std::left  << std::setw(11) << "mode"
            << std::right << std::setw(8)  << "nTrace"
            << std::setw(13) << "ref|r|"
            << std::setw(13) << "final|r|"
            << std::setw(12) << "dup-spread"
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
              << std::setw(11) << r.eu
              << std::setw(8)  << (r.deriv_ok?"ok":"BAD")
              << std::setw(9)  << (r.ok?"PASS":"FAIL") << "\n";
  }
  std::cout << "==========================================================================\n";
  std::cout << "Single-block scalar-advection (mixed fwd-x/bwd-y, 3D) sweep over IC_WEAK, IC_TRACE, IC_STRONG: "
            << (all_ok?"PASS":"FAIL") << "\n";
  return all_ok ? 0 : 1;
}
