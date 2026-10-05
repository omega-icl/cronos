// ===========================================================================
// OCFE_PDE13_solve2.cpp
//
// Scalar linear ADVECTION-DIFFUSION oracle -- the first stage of the PSA
// case study.  It isolates the two pieces of interface-condition machinery
// that the multi-physics PSA model stresses but that the suite has never
// actually exercised:
//
//   (1) the advection-diffusion CORRECTOR (the k_h vs k* = adv/diff
//       singular-perturbation guard in _symbol_balance_probe).  In every
//       prior test it has run INERT -- PDE5 is pure diffusion (k* = 0).
//       Here the advection term makes k* = v/D finite, so a Peclet sweep
//       drives the read across the mesh threshold and the corrector should
//       transition from inert (diffusion-dominated) to FIRED (advection-
//       dominated), shifting the z interface read FLUX -> VALUE/UPWIND.
//
//   (2) the DIRECTIONAL boundary read -- a Danckwerts inflow flux (Robin,
//       value+derivative in one equation) at z=LB plus a Neumann outflow
//       condition at z=UB, the canonical well-posed split for a parabolic
//       transport operator and the exact BC structure of the PSA composition
//       and energy balances (sec. 8-9 of the PSA model).
//
//       d_t c + v d_z c - D d_zz c = f(t,z)      on (t,z)
//
//       t evolution,  z spatial.  v > 0 (inflow at z=0), D > 0.
//
// This is deliberately ORDER 2 in z (one reduce_order aux Dz_c = c_z, chain
// depth 1, the order>2 boundary closure stays inert) and LINEAR in c (Newton
// is exact in one step), so convergence is never the variable under test --
// what we are watching is the SETUP: the corrector read, the resolved z
// interface type, squareness, and interface continuity.  The manufactured
// forcing makes c_exact the exact solution at every Peclet number, so accuracy
// is expected to hold across the whole sweep; the observable change is the
// corrector row in the probe table, not the error.
//
// Build flags:  -DMC__OCFESLV_INTERFACE_DECISION_V2  (+ -DCRONOS__WITH_SPQR for QR,
//   + -DCRONOS__WITH_UMFPACK for the IC_STRONG band path).  STRONGLY RECOMMENDED:
//   -DMC__OCFESLV_SYMBOL_BALANCE_PROBE  so the per-direction read table prints --
//   that table (the (c,z) row's k*, k_h, corrector columns) is the whole point.
//
// Knobs (-D overrides):
//   -DTEST_AD_NEL_T / _NEL_Z       finite elements    (default 3 / 3)
//   -DTEST_AD_NT   / _NZ           nodes per element  (default 6 / 12)
//   -DTEST_AD_V    -DTEST_AD_D     advection / diffusion coefficients
//   -DTEST_AD_K  -DTEST_AD_TF  -DTEST_AD_ZF   wavenumber / final time / length
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

#ifndef TEST_AD_NEL_T
#define TEST_AD_NEL_T 3
#endif
#ifndef TEST_AD_NEL_Z
#define TEST_AD_NEL_Z 3
#endif
#ifndef TEST_AD_NT
#define TEST_AD_NT 6
#endif
#ifndef TEST_AD_NZ
#define TEST_AD_NZ 12
#endif
#ifndef TEST_AD_V
#define TEST_AD_V 1.0      // advection velocity (v > 0 => inflow at z=LB)
#endif
#ifndef TEST_AD_D
#define TEST_AD_D 0.02     // diffusion (default Pe = v*zf/D = 50)
#endif
#ifndef TEST_AD_TF
#define TEST_AD_TF 0.5
#endif
#ifndef TEST_AD_ZF
#define TEST_AD_ZF 1.0
#endif
#ifndef TEST_AD_MAXIT
#define TEST_AD_MAXIT 40
#endif
#ifndef TEST_AD_SOLVE_TOL
#define TEST_AD_SOLVE_TOL 1e-9
#endif
#ifndef TEST_AD_EXACT_TOL
#define TEST_AD_EXACT_TOL 1e-6
#endif

#define TEST_AD_SPQR

static double const kPi = 3.14159265358979323846;

struct Par {
  double tf  = TEST_AD_TF;
  double zf  = TEST_AD_ZF;
  double v   = TEST_AD_V;    // advection velocity
  double D   = TEST_AD_D;    // diffusion coefficient
  double k   = 2.0*kPi;      // wavenumber
  double phi = 0.7;          // phase (keeps c and c_z nonzero at the ends)
  double a   = 0.3;          // linear-in-time growth
  double c0  = 1.0;          // amplitude
};

// Manufactured solution and the z-derivative needed for the aux comparison,
// the Danckwerts/Neumann BC data and the initial guess:
//   c = c0 sin(k z + phi) (1 + a t).
static double C_exact ( double t, double z, Par const& p ){ return p.c0*std::sin(p.k*z+p.phi)*(1.0+p.a*t); }
static double Cz_exact( double t, double z, Par const& p ){ return p.c0*p.k*std::cos(p.k*z+p.phi)*(1.0+p.a*t); }

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
// (duplicate) collocation nodes.  A non-zero spread at an interior interface
// is a continuity jump; lets us see whether c or its aux Dz_c jumps.
struct DuplicateSpread { double max_pair=0.0; size_t max_mult=0; };

static DuplicateSpread duplicate_node_spread
( OCFESLV const& oc, FFVar const& st, std::vector<double> const& var, size_t off )
{
  struct Accum { double lo, hi; size_t count; };
  std::map< std::vector<long long>, Accum > groups;
  auto nodes = oc.node_colloc(st);
  for( size_t i=0; i<nodes.size(); ++i ){
    std::vector<long long> key; key.reserve(nodes[i].size());
    for( double cc: nodes[i] ) key.push_back( static_cast<long long>( std::llround(cc*1.0e12) ) );
    double const vv = var[off+i];
    auto it = groups.find(key);
    if( it == groups.end() ) groups.emplace( std::move(key), Accum{vv,vv,1} );
    else { it->second.lo=std::min(it->second.lo,vv); it->second.hi=std::max(it->second.hi,vv); ++it->second.count; }
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
                    bool solved=false; double final_res=0., eC=0.;
                    double Pe=0., Dval=0.; };

static ModeResult run_mode( OCFESLV::Options::ImpositionType imp,
                            std::string const& suffix, Par const& p )
{
  (void)suffix;
  ModeResult R; R.name = imp_name(imp); bool ok=true;
  R.Dval = p.D; R.Pe = p.v*p.zf/p.D;

  std::cout << "\n========== scalar advection-diffusion test ==========\n";
  std::cout << "imposition: " << R.name
            << ", finite elements: t=" << TEST_AD_NEL_T << " z=" << TEST_AD_NEL_Z
            << ", nodes/element: t=" << TEST_AD_NT << " z=" << TEST_AD_NZ << "\n";
  std::cout << "PDE: d_t c + v d_z c - D d_zz c = f,  v=" << p.v << " D=" << p.D
            << " (Pe=v*zf/D=" << R.Pe << "),  c_exact = c0 sin(k z + phi)(1 + a t),"
            << " k=" << p.k << "\n";

  FFGraph DAG;
  FFVar t = DAG.add_var("t");
  FFVar z = DAG.add_var("z");
  FFVar c = DAG.add_var("c(t,z)");
  FFPartial OpP;

  // Manufactured solution, forcing and BC data, all in CLOSED FORM.  OpP is
  // applied only to the state c (so reduce_order turns OpP(c,{z,2}) into the
  // single aux Dz_c, and reuse_aux folds the advection OpP(c,z) onto the same
  // aux); the manufactured expressions are never differentiated with OpP --
  // doing so would create a derivative-equation in only the domain variables
  // (t,z) with no state, which the domain-consistency check (Z33) rejects.
  //   c_exact = c0 sin(k z + phi)(1 + a t)
  //   c_z     = c0 k cos(k z + phi)(1 + a t)
  //   f = d_t c + v d_z c - D d_zz c
  //     = c0 a sin(.) + v c0 k (1+at) cos(.) + D c0 k^2 (1+at) sin(.)
  FFVar sinz = sin( p.k*z + p.phi );
  FFVar cosz = cos( p.k*z + p.phi );
  FFVar CE  = p.c0*sinz*( 1.0 + p.a*t );                       // c_exact
  FFVar FE  = p.c0*p.a*sinz
            + p.v*p.c0*p.k*cosz*( 1.0 + p.a*t )                // v c_z
            + p.D*p.c0*(p.k*p.k)*sinz*( 1.0 + p.a*t );         // -D c_zz = +D k^2 sin

  // Danckwerts inflow flux RHS at z=LB:  g_in(t) = v c_exact(t,0) - D c_z(t,0).
  double const sphi = std::sin(p.phi), cphi = std::cos(p.phi);
  FFVar GIN = ( p.v*p.c0*sphi - p.D*p.c0*p.k*cphi )*( 1.0 + p.a*t );
  // Neumann outflow RHS at z=UB:  h_out(t) = c_z(t,L).
  double const ckL = std::cos(p.k*p.zf + p.phi);
  FFVar HOUT = p.c0*p.k*ckL*( 1.0 + p.a*t );

  FFVar PDE    = OpP(c,t) + p.v*OpP(c,z) - p.D*OpP(c,{z,2}) - FE;  // d_t c + v d_z c - D d_zz c - f
  FFVar BCINIT = c - CE;                                           // initial value c(0,z)=c_exact
  FFVar BCROB  = p.v*c - p.D*OpP(c,z) - GIN;                       // Danckwerts inflow flux at z=LB
  FFVar BCNEU  = OpP(c,z) - HOUT;                                  // Neumann outflow at z=UB

  OCFESLV oc(&DAG);
  oc.add_domain( t, FFDom(0., p.tf, TEST_AD_NEL_T, FFDom::LGR, TEST_AD_NT) );
  oc.add_domain( z, FFDom(0., p.zf, TEST_AD_NEL_Z, FFDom::LGL, TEST_AD_NZ) );
  oc.add_state( c, {t,z} );
  oc.update_ref( c, [&]( OCFESLV::t_Coord const& crd ){ return C_exact(crd.at(t),crd.at(z),p); } );

  OCFESLV::EqnOptions int_opt( OCFESLV::EqnRole::INTERIOR, 0 );
  OCFESLV::EqnOptions ini_opt( OCFESLV::EqnRole::INITIAL,  0 );
  OCFESLV::EqnOptions bnd_opt( OCFESLV::EqnRole::BOUNDARY, 0 );

  int const T_INT = FFDom::ALL - FFDom::LB;
  int const Z_INT = FFDom::ALL - FFDom::LB - FFDom::UB;

  // Interior advection-diffusion PDE on z-interior, t>0.
  oc.add_equation( PDE,    {t,z}, {T_INT, Z_INT},         int_opt );
  // Initial data on the whole spatial grid.
  oc.add_equation( BCINIT, {t,z}, {FFDom::LB, FFDom::ALL}, ini_opt );
  // Two spatial conditions for the 2nd-order operator: Danckwerts inflow flux
  // (Robin) at z=LB + Neumann outflow at z=UB.  The Robin condition carries the
  // value term (v c), so the system has no constant null space.
  oc.add_equation( BCROB,  {t,z}, {T_INT, FFDom::LB},     bnd_opt );  // v c - D c_z = g_in  (inflow)
  oc.add_equation( BCNEU,  {t,z}, {T_INT, FFDom::UB},     bnd_opt );  // c_z = h_out         (outflow)

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

  // Diagnostics: the reduction chain (depth-1 -> one aux Dz_c), the block
  // classification, and the RESOLVED z interface read -- the corrector's
  // effect on z shows up here (FLUX when diffusion-dominated, shifting toward
  // VALUE/UPWIND when advection dominates at high Pe).
  std::cout << "States after setup (primitive + reduce_order aux chain):";
  for( auto const& st: oc.states_colloc() ) std::cout << ' ' << st.name();
  std::cout << "\nAuxiliary states introduced: "
            << (oc.states_colloc().size()>1?oc.states_colloc().size()-1:0)
            << " (expect 1: D1=c_z)\n";
  {
    auto const& cls = oc.pde_type();
    std::cout << "PDE type: " << OCFESLV::pde_type_name(cls.type)
              << "  evolution_hyperbolic=" << (cls.evolution_hyperbolic?"yes":"no")
              << "  parabolic=" << (cls.parabolic_structure_detected?"yes":"no") << "\n";
    bool const weak_path = ( imp != OCFESLV::Options::IC_STRONG );
    std::cout << "Resolved interface type:  t(evolution) -> "
              << iftype_name( oc.resolved_interface_type( 0, t, OCFESLV::EqnRole::INTERIOR, weak_path ) )
              << "   |   z -> "
              << iftype_name( oc.resolved_interface_type( 0, z, OCFESLV::EqnRole::INTERIOR, weak_path ) )
              << "\n";
  }

  // Initial guess: exact primitive, aux left at zero (the LINK drives Dz_c to
  // c_z during the solve).  The problem is LINEAR, so Newton is exact in one
  // step regardless of the seed -- this just keeps the seed honest.
  std::vector<double> var(nVar,0.0);
  {
    size_t off=0, sidx=0;
    for( auto const& st: oc.states_colloc() ){
      auto nodes = oc.node_colloc(st);
      bool const is_primitive = ( sidx == 0 );
      for( size_t i=0; i<nodes.size(); ++i )
        var[off+i] = is_primitive ? C_exact( nodes[i][0], nodes[i][1], p ) : 0.0;
      off += nodes.size(); ++sidx;
    }
  }

  oc.options.SOLVE.MAX_ITER = TEST_AD_MAXIT;
  oc.options.SOLVE.RES_TOL  = TEST_AD_SOLVE_TOL;
#if defined(CRONOS__WITH_SPQR) && defined(TEST_AD_SPQR)
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
  ok &= check_close("final max residual", R.final_res, TEST_AD_SOLVE_TOL);

  // Primitive accuracy against the manufactured solution.
  {
    double eC=0.0; size_t off=0, sidx=0;
    for( auto const& st: oc.states_colloc() ){
      auto nodes = oc.node_colloc(st);
      if( sidx == 0 )
        for( size_t i=0; i<nodes.size(); ++i )
          eC = std::max( eC, std::fabs( var[off+i] - C_exact(nodes[i][0],nodes[i][1],p) ) );
      off += nodes.size(); ++sidx;
    }
    R.eC = eC;
    ok &= check_close("max |c - c_exact|", R.eC, TEST_AD_EXACT_TOL);
  }

  // Per-state error breakdown (primitive vs aux c_z) and per-state interface
  // continuity (which field jumps across element faces).
  {
    size_t off=0, sidx=0;
    for( auto const& st: oc.states_colloc() ){
      auto nodes = oc.node_colloc(st);
      double e=0.0;
      for( size_t i=0; i<nodes.size(); ++i ){
        double const ex = (sidx==0)? C_exact (nodes[i][0],nodes[i][1],p)
                                   : Cz_exact(nodes[i][0],nodes[i][1],p);
        e = std::max( e, std::fabs( var[off+i] - ex ) );
      }
      std::cout << "  max|" << std::left << std::setw(12) << (st.name()+" - exact")
                << "|=" << std::scientific << std::setprecision(6) << e << "\n";
      off += nodes.size(); ++sidx;
    }
    print_duplicate_spreads( oc, var );
  }

  std::cout << "scalar advection-diffusion test (" << R.name << "): " << (ok?"PASS":"FAIL") << "\n";
  R.ok = ok;
  return R;
}

int main()
{
  bool all_ok = true;

  // ---------------------------------------------------------------------- //
  // Section A -- correctness across imposition modes at the default Peclet. //
  // ---------------------------------------------------------------------- //
  Par p;
  std::vector<ModeResult> modeA;
  modeA.push_back( run_mode(OCFESLV::Options::IC_WEAK,   "weak",   p) );
  modeA.push_back( run_mode(OCFESLV::Options::IC_TRACE,  "trace",  p) );
  modeA.push_back( run_mode(OCFESLV::Options::IC_STRONG, "strong", p) );

  std::cout << "\n============== PDE13 mode sweep (Pe=" << p.v*p.zf/p.D
            << ") summary ==============\n";
  std::cout << std::left << std::setw(12) << "mode"
            << std::setw(10) << "nTrace" << std::setw(16) << "final|r|"
            << std::setw(16) << "|c-exact|" << "result\n";
  for( auto const& r: modeA ){
    std::cout << std::left << std::setw(12) << r.name
              << std::setw(10) << r.nTrace
              << std::scientific << std::setprecision(3)
              << std::setw(16) << r.final_res
              << std::setw(16) << r.eC
              << (r.ok?"PASS":"FAIL") << "\n";
    all_ok &= r.ok;
  }
  std::cout << "============================================================\n";

  // ---------------------------------------------------------------------- //
  // Section B -- Peclet transition sweep at IC_WEAK.  Watch the (c,z) row of //
  // the _symbol_balance_probe table across these runs: k* = v/D grows, and   //
  // the corrector should move from inert (low Pe) to FIRED (high Pe).  The   //
  // solution stays exact throughout (manufactured forcing), so the error     //
  // column is a control, not the observable.                                 //
  // ---------------------------------------------------------------------- //
  double const Dlist[] = { 1.0, 0.1, 0.01, 0.001 };   // Pe = 1, 10, 100, 1000
  std::vector<ModeResult> sweepB;
  for( double Dval : Dlist ){
    Par pD = p; pD.D = Dval;
    sweepB.push_back( run_mode(OCFESLV::Options::IC_WEAK, "weak", pD) );
  }

  std::cout << "\n============== PDE13 Peclet sweep (IC_WEAK) summary ==============\n";
  std::cout << std::left << std::setw(12) << "D"
            << std::setw(12) << "Pe=v*zf/D" << std::setw(14) << "k*=v/D"
            << std::setw(16) << "final|r|" << std::setw(16) << "|c-exact|"
            << "result\n";
  for( auto const& r: sweepB ){
    std::cout << std::left << std::scientific << std::setprecision(3)
              << std::setw(12) << r.Dval
              << std::setw(12) << r.Pe
              << std::setw(14) << (p.v/r.Dval)
              << std::setw(16) << r.final_res
              << std::setw(16) << r.eC
              << (r.ok?"PASS":"FAIL") << "\n";
    all_ok &= r.ok;
  }
  std::cout << "=================================================================\n";
  std::cout << "(corrector transition is read from the (c,z) probe row per run, "
               "not from this table)\n";

  std::cout << "\nscalar advection-diffusion oracle (modes + Peclet sweep): "
            << (all_ok?"PASS":"FAIL") << "\n";
  return all_ok ? 0 : 1;
}
