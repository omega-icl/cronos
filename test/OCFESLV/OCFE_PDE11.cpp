// ===========================================================================
// OCFE_PDE11_solve2.cpp
//
// First single-direction ORDER>2 test -- a 3rd-order (dispersive / linearized-
// KdV) operator, the oracle for the deferred item-11 Part-2 order>2 read.
//
//   d_t u + d_xxx u = f(t,x)      on (t,x), t evolution, x spatial
//
// PURPOSE (phase 1a -- reduction + solve through the current guard).  The read
// presently GUARDS order>2 (rd.order_gt2 -> block-level fallback), so this PDE
// should solve via the fallback; the point of phase 1a is to (a) exercise
// reduce_order at chain depth 2 for the FIRST time (D1=u_x, D2=u_xx; LINK chain
// D1=d_x u, D2=d_x D1), which the suite has never tested beyond depth 1, and
// (b) surface the boundary-condition mechanics of a reduced 3rd-order operator.
//
// BC NOTE (the structural unknown this build resolves).  After reduction the
// block carries 3 states (u,D1,D2); each node takes 2 auto-collocated LINK rows
// + 1 user row, so a 1-D x-domain offers only TWO boundary-node slots (x=LB,
// x=UB) yet a 3rd-order operator needs THREE x-conditions for well-posedness
// (two value BCs leave a 1-D dispersive null space).  The third condition (a
// derivative BC at x=LB) therefore shares the x=LB node with the value BC, which
// over-fills that node by one unless the framework's redundancy handling drops
// the now-redundant LINK1 there (D1=d_x u and D1=UE_x both pin D1, consistently).
// Whether that resolves cleanly is exactly what this first build tells us; the
// BC block below is written to be easy to retune from the build's square/rank
// report.
//
// Build flags:  -DMC__OCFESLV_INTERFACE_DECISION_V2  (+ -DCRONOS__WITH_SPQR for QR)
//   optional probe: -DMC__OCFESLV_SYMBOL_BALANCE_PROBE to print the per-direction
//   read table (should show order>2 flagged, deferring to block-level).
//
// Knobs (-D overrides):
//   -DTEST_KDV_NEL_T / _NEL_X      finite elements    (default 3 / 3)
//   -DTEST_KDV_NT   / _NX          nodes per element  (default 6 / 7)
//   -DTEST_KDV_K  -DTEST_KDV_TF    wavenumber / final time
// ===========================================================================

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <string>
#include <vector>

#include <armadillo>

#ifndef OCFE_OCFESLV_HEADER
#define OCFE_OCFESLV_HEADER "ocfeslv.hpp"
#endif
#include OCFE_OCFESLV_HEADER

using namespace mc;

#ifndef TEST_KDV_NEL_T
#define TEST_KDV_NEL_T 3
#endif
#ifndef TEST_KDV_NEL_X
#define TEST_KDV_NEL_X 3
#endif
#ifndef TEST_KDV_NT
#define TEST_KDV_NT 6
#endif
#ifndef TEST_KDV_NX
#define TEST_KDV_NX 12   // per-element degree: order-3 chain double-differentiates
                         // D_xx via the aux chain, so low degree (e.g. 7) is
                         // under-resolved; 12 converges with comfortable margin.
#endif
#ifndef TEST_KDV_TF
#define TEST_KDV_TF 0.5
#endif
#ifndef TEST_KDV_XF
#define TEST_KDV_XF 1.0
#endif
#ifndef TEST_KDV_MAXIT
#define TEST_KDV_MAXIT 40
#endif
#ifndef TEST_KDV_SOLVE_TOL
#define TEST_KDV_SOLVE_TOL 1e-9
#endif
#ifndef TEST_KDV_EXACT_TOL
#define TEST_KDV_EXACT_TOL 1e-6
#endif

#define TEST_KDV_SPQR

// Third (derivative) boundary condition u_x(x=LB)=UE_x.  Always imposed: a
// 3rd-order operator needs three spatial conditions, and the value+derivative
// pair at x=LB is the over-determining face that the order>2 boundary closure
// resolves by displacing the redundant LINK1 (validated against the single-
// element oracle and shown to converge spectrally on the multi-element grid).

static double const kPi = 3.14159265358979323846;

struct Par {
  double tf  = TEST_KDV_TF;
  double xf  = TEST_KDV_XF;
  double k   = 2.0*kPi;     // wavenumber
  double phi = 0.7;         // phase (keeps u and its derivatives nonzero at the ends)
  double a   = 0.3;         // linear-in-time growth
  double u0  = 1.0;         // amplitude
};

// Manufactured solution and the x-derivatives we need for the derivative BC and
// the aux comparison: u = u0 sin(k x + phi) (1 + a t).
static double U_exact ( double t, double x, Par const& p ){ return p.u0*std::sin(p.k*x+p.phi)*(1.0+p.a*t); }
static double Ux_exact( double t, double x, Par const& p ){ return p.u0*p.k*std::cos(p.k*x+p.phi)*(1.0+p.a*t); }
static double Uxx_exact(double t, double x, Par const& p ){ return -p.u0*p.k*p.k*std::sin(p.k*x+p.phi)*(1.0+p.a*t); }

static char const* iftype_name( OCFESLV::Options::InterfaceType t ){
  switch( t ){
    case OCFESLV::Options::IC_VALUE:  return "VALUE";
    case OCFESLV::Options::IC_UPWIND: return "UPWIND";
    case OCFESLV::Options::IC_FLUX:   return "FLUX";
    case OCFESLV::Options::IC_AUTO:   return "AUTO";
    default:                        return "?";
  }
}
static char const* imp_name( OCFESLV::Options::ImpositionType t ){
  switch( t ){
    case OCFESLV::Options::IC_WEAK:   return "IC_WEAK";
    case OCFESLV::Options::IC_TRACE:  return "IC_TRACE";
    case OCFESLV::Options::IC_STRONG: return "IC_STRONG";
    default:                        return "IC_?";
  }
}

static double max_abs( std::vector<double> const& v ){
  double m=0.0; for( double x: v ) m = std::max(m,std::fabs(x)); return m;
}
static bool check_close( std::string const& label, double value, double tol ){
  bool ok = ( value <= tol );
  std::cout << std::left << std::setw(40) << label
            << " value=" << std::scientific << std::setprecision(6) << value
            << " tol=" << tol << "  " << (ok?"PASS":"FAIL") << "\n";
  return ok;
}

// Per-state interface continuity probe: spread (max-min) among co-located
// (duplicate) collocation nodes.  At interior element interfaces a state's
// duplicate nodes share coordinates; a non-zero spread is a continuity jump.
// Lets us see WHICH state of the reduced chain (u / Dx_u / Dx_Dx_u) jumps.
struct DuplicateSpread { double max_pair=0.0; size_t max_mult=0; };

static DuplicateSpread duplicate_node_spread
( OCFESLV const& oc, FFVar const& st, std::vector<double> const& var, size_t off )
{
  struct Accum { double lo, hi; size_t count; };
  std::map< std::vector<long long>, Accum > groups;
  auto nodes = oc.node_colloc(st);
  for( size_t i=0; i<nodes.size(); ++i ){
    std::vector<long long> key; key.reserve(nodes[i].size());
    for( double c: nodes[i] ) key.push_back( static_cast<long long>( std::llround(c*1.0e12) ) );
    double const v = var[off+i];
    auto it = groups.find(key);
    if( it == groups.end() ) groups.emplace( std::move(key), Accum{v,v,1} );
    else { it->second.lo=std::min(it->second.lo,v); it->second.hi=std::max(it->second.hi,v); ++it->second.count; }
  }
  DuplicateSpread out;
  for( auto const& kv: groups ){
    out.max_mult = std::max(out.max_mult, kv.second.count);
    if( kv.second.count >= 2 ) out.max_pair = std::max(out.max_pair, kv.second.hi-kv.second.lo);
  }
  return out;
}

static void print_duplicate_spreads
( OCFESLV const& oc, std::vector<double> const& var )
{
  std::cout << "Element-interface duplicate-node spreads (per state):\n";
  size_t off=0;
  for( auto const& st: oc.states_colloc() ){
    DuplicateSpread const d = duplicate_node_spread(oc,st,var,off);
    std::cout << "  " << std::left << std::setw(16) << st.name()
              << " max_pair_spread=" << std::scientific << std::setprecision(6) << d.max_pair
              << "  max_multiplicity=" << d.max_mult << "\n";
    off += oc.node_colloc(st).size();
  }
}

struct ModeResult { std::string name; bool ok=false; size_t nVar=0,nEqn=0,nTrace=0;
                    bool solved=false; double final_res=0., eU=0.; };

static ModeResult run_mode( OCFESLV::Options::ImpositionType imp,
                            std::string const& suffix, Par const& p )
{
  (void)suffix;
  ModeResult R; R.name = imp_name(imp); bool ok=true;

  std::cout << "\n========== 3rd-order dispersive (KdV-type) test ==========\n";
  std::cout << "imposition: " << R.name
            << ", finite elements: t=" << TEST_KDV_NEL_T << " x=" << TEST_KDV_NEL_X
            << ", nodes/element: t=" << TEST_KDV_NT << " x=" << TEST_KDV_NX << "\n";
  std::cout << "PDE: d_t u + d_xxx u = f,  u_exact = u0 sin(k x + phi)(1 + a t),"
            << " k=" << p.k << "\n";

  FFGraph DAG;
  FFVar t = DAG.add_var("t");
  FFVar x = DAG.add_var("x");
  FFVar u = DAG.add_var("u(t,x)");
  FFPartial OpP;

  // Manufactured solution and the forcing, in CLOSED FORM.  OpP is applied only
  // to the state u (so reduce_order turns OpP(u,{x,3}) into the D1/D2 chain); the
  // manufactured expression UE is never differentiated with OpP -- doing so would
  // create a placeholder derivative-equation in only the domain variables (t,x)
  // with no state, which the domain-consistency check rejects.  Exact terms:
  //   UE      = u0 sin(k x + phi)(1 + a t)
  //   dUE/dx  = u0 k cos(k x + phi)(1 + a t)
  //   f = dUE/dt + d3UE/dx3 = u0 a sin(k x + phi) - u0 k^3 cos(k x + phi)(1 + a t)
  FFVar sinx = sin( p.k*x + p.phi );
  FFVar cosx = cos( p.k*x + p.phi );
  FFVar UE  = p.u0*sinx*( 1.0 + p.a*t );
  FFVar UEx = p.u0*p.k*cosx*( 1.0 + p.a*t );
  FFVar FE  = p.u0*p.a*sinx - p.u0*(p.k*p.k*p.k)*cosx*( 1.0 + p.a*t );

  FFVar PDE   = OpP(u,t) + OpP(u,{x,3}) - FE;   // d_t u + d_xxx u - f
  FFVar BCVAL = u - UE;                         // Dirichlet value
  FFVar BCDX  = OpP(u,x) - UEx;                 // first-derivative (Neumann) condition

  OCFESLV oc(&DAG);
  oc.add_domain( t, FFDom(0., p.tf,  TEST_KDV_NEL_T, FFDom::LGR, TEST_KDV_NT) );
  oc.add_domain( x, FFDom(0., p.xf,  TEST_KDV_NEL_X, FFDom::LGL, TEST_KDV_NX) );
  oc.add_state( u, {t,x} );
  oc.update_ref( u, [&]( OCFESLV::t_Coord const& c ){ return U_exact(c.at(t),c.at(x),p); } );

  OCFESLV::EqnOptions int_opt( OCFESLV::EqnRole::INTERIOR, 0 );
  OCFESLV::EqnOptions ini_opt( OCFESLV::EqnRole::INITIAL,  0 );
  OCFESLV::EqnOptions bnd_opt( OCFESLV::EqnRole::BOUNDARY, 0 );

  int const T_INT = FFDom::ALL - FFDom::LB;
  int const X_INT = FFDom::ALL - FFDom::LB - FFDom::UB;

  // Interior PDE on x-interior, t>0.
  oc.add_equation( PDE, {t,x}, {T_INT, X_INT}, int_opt );
  // Initial data on the whole spatial grid.
  oc.add_equation( BCVAL, {t,x}, {FFDom::LB, FFDom::ALL}, ini_opt );
  // x-boundary conditions for t>0.  Three conditions for the 3rd-order operator:
  // value at both ends + first-derivative at x=LB (2 at LB, 1 at UB -- a standard
  // well-posed split for d_xxx).  The value+derivative pair at x=LB is the node
  // that over-fills if the framework does not drop the redundant LINK1 there.
  oc.add_equation( BCVAL, {t,x}, {T_INT, FFDom::LB}, bnd_opt );   // u(x=LB) = UE
  oc.add_equation( BCDX,  {t,x}, {T_INT, FFDom::LB}, bnd_opt );   // u_x(x=LB) = UE_x  (3rd condition; closure displaces LINK1)
  oc.add_equation( BCVAL, {t,x}, {T_INT, FFDom::UB}, bnd_opt );   // u(x=UB) = UE

  oc.set_evolution_domain( t );
  oc.options.REDUCE.ORDER    = OCFESLV::Options::RED_FULL;
  oc.options.CLASSIFY.MODE        = OCFESLV::Options::CLASS_AUTO;
  oc.options.INTERFACE.TYPE  = OCFESLV::Options::IC_AUTO;
  oc.options.INTERFACE.IMPOSITION = imp;
  oc.options.INTERFACE.SAT_SIGMA0      = 10.0;
  oc.options.DISPLAY_LEVEL   = 1;

  if( !oc.setup() ){
    std::cerr << "ERROR: OCFESLV::setup() failed for " << R.name << "\n";
    R.ok=false; return R;
  }

  size_t const nVar=oc.n_colloc_sta();
  size_t const nEqn=oc.n_colloc_eqn();
  size_t const nTrace=oc.n_colloc_trace();
  R.nVar=nVar; R.nEqn=nEqn; R.nTrace=nTrace;
  std::cout << "nVar=" << nVar << " nEqn=" << nEqn << " nTrace=" << nTrace
            << " square=" << (nVar==nEqn?"yes":"no") << "\n";
  ok &= check_close("square system (nVar==nEqn)", nVar==nEqn?0.0:1.0, 0.0);

  // Diagnostics: the reduction chain (does depth-2 produce two auxes?), the
  // block classification, and the resolver verdict for x (order>2 => the read
  // guards and the block-level fallback governs).
  std::cout << "States after setup (primitive + reduce_order aux chain):";
  for( auto const& st: oc.states_colloc() ) std::cout << ' ' << st.name();
  std::cout << "\nAuxiliary states introduced: "
            << (oc.states_colloc().size()>1?oc.states_colloc().size()-1:0)
            << " (expect 2: D1=u_x, D2=u_xx)\n";
  {
    auto const& cls = oc.pde_type();
    std::cout << "PDE type: " << OCFESLV::pde_type_name(cls.type)
              << "  evolution_hyperbolic=" << (cls.evolution_hyperbolic?"yes":"no")
              << "  parabolic=" << (cls.parabolic_structure_detected?"yes":"no") << "\n";
    bool const weak_path = ( imp != OCFESLV::Options::IC_STRONG );
    std::cout << "Resolved interface type:  t(evolution) -> "
              << iftype_name( oc.resolved_interface_type( 0, t, OCFESLV::EqnRole::INTERIOR, weak_path ) )
              << "   |   x -> "
              << iftype_name( oc.resolved_interface_type( 0, x, OCFESLV::EqnRole::INTERIOR, weak_path ) )
              << "\n";
  }

  // Initial guess: exact primitive, auxes left at zero (the LINK chain drives
  // them to u_x, u_xx during the solve).  Phase 1a does not yet seed the aux
  // exact values -- the printed aux names above tell us the chain naming for the
  // next iteration's reference-residual / aux-error gates.
  std::vector<double> var(nVar,0.0);
  {
    size_t off=0, sidx=0;
    for( auto const& st: oc.states_colloc() ){
      auto nodes = oc.node_colloc(st);          // each node: coords in domain order
      bool const is_primitive = ( sidx == 0 );  // primitive is first; reduce_order auxes follow
      for( size_t i=0; i<nodes.size(); ++i )
        var[off+i] = is_primitive ? U_exact( nodes[i][0], nodes[i][1], p ) : 0.0;
      off += nodes.size(); ++sidx;
    }
  }

  // This driver hand is for monolithic only
  oc.options.SOLVE.MARCHING = false;
  oc.options.SOLVE.MAX_ITER = TEST_KDV_MAXIT;
  oc.options.SOLVE.RES_TOL  = TEST_KDV_SOLVE_TOL;
#if defined(CRONOS__WITH_SPQR) && defined(TEST_KDV_SPQR)
  oc.options.SOLVE.FACTORIZATION = OCFESLV::Options::SOLVE_SPQR;
#endif
  OCFESLV::SolveReport const srep = oc.solve( var.data() );
  R.solved = srep.converged;
  if( !R.solved )
    std::cerr << "OCFESLV::solve did not converge: final|r|=" << srep.final_residual
              << " after " << srep.iterations << " it\n";
  ok &= R.solved;

  std::vector<double> res(nEqn,0.0);
  oc.eval(res.data(),nullptr,var.data(),nullptr,nullptr);
  R.final_res = max_abs(res);
  ok &= check_close("final max residual", R.final_res, TEST_KDV_SOLVE_TOL);

  // Primitive accuracy against the manufactured solution.
  {
    double eU=0.0; size_t off=0, sidx=0;
    for( auto const& st: oc.states_colloc() ){
      auto nodes = oc.node_colloc(st);
      if( sidx == 0 )                            // primitive
        for( size_t i=0; i<nodes.size(); ++i )
          eU = std::max( eU, std::fabs( var[off+i] - U_exact(nodes[i][0],nodes[i][1],p) ) );
      off += nodes.size(); ++sidx;
    }
    R.eU = eU;
    ok &= check_close("max |u - u_exact|", R.eU, TEST_KDV_EXACT_TOL);
  }

  // Per-state error breakdown (which field of the reduced chain is off) and
  // per-state interface continuity (which field jumps across element faces).
  {
    size_t off=0, sidx=0;
    for( auto const& st: oc.states_colloc() ){
      auto nodes = oc.node_colloc(st);
      double e=0.0;
      for( size_t i=0; i<nodes.size(); ++i ){
        double const ex = (sidx==0)? U_exact (nodes[i][0],nodes[i][1],p)
                        : (sidx==1)? Ux_exact(nodes[i][0],nodes[i][1],p)
                                   : Uxx_exact(nodes[i][0],nodes[i][1],p);
        e = std::max( e, std::fabs( var[off+i] - ex ) );
      }
      std::cout << "  max|" << std::left << std::setw(12) << (st.name()+" - exact")
                << "|=" << std::scientific << std::setprecision(6) << e << "\n";
      off += nodes.size(); ++sidx;
    }
    print_duplicate_spreads( oc, var );
  }

  std::cout << "3rd-order dispersive test (" << R.name << "): " << (ok?"PASS":"FAIL") << "\n";
  R.ok = ok;
  return R;
}

int main()
{
  Par p;
  std::vector<ModeResult> results;
  results.push_back( run_mode(OCFESLV::Options::IC_WEAK,   "weak",   p) );
  results.push_back( run_mode(OCFESLV::Options::IC_TRACE,  "trace",  p) );
  results.push_back( run_mode(OCFESLV::Options::IC_STRONG, "strong", p) );

  std::cout << "\n==================== PDE11 (3rd-order) sweep summary ====================\n";
  std::cout << std::left << std::setw(12) << "mode"
            << std::setw(10) << "nTrace" << std::setw(16) << "final|r|"
            << std::setw(16) << "|u-exact|" << "result\n";
  bool all_ok=true;
  for( auto const& r: results ){
    std::cout << std::left << std::setw(12) << r.name
              << std::setw(10) << r.nTrace
              << std::scientific << std::setprecision(3)
              << std::setw(16) << r.final_res
              << std::setw(16) << r.eU
              << (r.ok?"PASS":"FAIL") << "\n";
    all_ok &= r.ok;
  }
  std::cout << "========================================================================\n";
  std::cout << "3rd-order dispersive sweep over IC_WEAK, IC_TRACE, IC_STRONG: "
            << (all_ok?"PASS":"FAIL") << "\n";
  return all_ok ? 0 : 1;
}
