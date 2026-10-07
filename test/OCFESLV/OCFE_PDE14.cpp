// ===========================================================================
// OCFE_PDE14_solve2.cpp
//
// Off-diagonal first-order HYPERBOLIC oracle -- stage 2 of the PSA case study.
// Where PDE13 (scalar advection-diffusion) is diagonal -- one state, one
// governing equation, so a per-state read suffices and VALUE happens to be the
// right interface answer at every Peclet number -- this driver is the minimal
// system where a per-state read CANNOT recover the interface structure.
//
//   d_t c + a d_z u = f_c(t,z)        on (t,z)
//   d_t u + a d_z c = f_u(t,z)
//
//   t evolution,  z spatial.  The z principal symbol is a*[[0,1],[1,0]] --
//   ZERO on the diagonal, eigenvalues +-a, characteristic variables w+ = c+u
//   (speed +a) and w- = c-u (speed -a) travelling in OPPOSITE directions.
//
// This is the PSA continuity/momentum analog: the equation is first-order in
// BOTH states, no single state owns it, and the well-posed conditions live on
// the characteristic combinations c+-u, not on c or u individually --
//   w+ = c+u  incoming at z=LB  (prescribe there)
//   w- = c-u  incoming at z=UB  (prescribe there)
// A per-state read finds each state's own z-derivative only in the OTHER
// equation (the diagonal speed is zero), so it has no upwind direction to read
// and cannot split the conditions between the ends.  Recovering this requires
// the joint eigenstructure of the 2x2 block -- a per-equation, per-direction,
// sign-aware characteristic read.
//
// OBSERVED (gap-finder run, TEST_HYP_MANUAL_CLOSURE=0):
//   Setup SUCCEEDS and the framework correctly classifies the block
//   EVOL_HYPERBOLIC and reads (c,z)/(u,z) as order-1 UPWIND -- the read and
//   classifier recognise the structure.  But the system is NOT square:
//   nEqn-nVar = -2*(t-interior nodes) = -34 on the 3x3 / 6x10 grid, because a
//   first-order block has no reduce_order LINK to close the boundary and the
//   PDE is collocated on the z-interior, so the OUTGOING characteristic at each
//   end (w- at z=LB, w+ at z=UB) is left unequationed.  Interface continuity is
//   fine (IC_TRACE/IC_STRONG dup-spreads ~1e-15) but the 34-dim null space is
//   unconstrained, so the solve reaches |r|~0 at the WRONG solution
//   (|c-exact|~0.87).  Diagnosis: the gap is the hyperbolic BOUNDARY CLOSURE,
//   not the read, the classifier, or interface continuity.
//
// DEFAULT (framework auto-closure, TEST_HYP_MANUAL_CLOSURE=0): the header's
//   auto-closure (enabled by default) projects the collocated PDEs onto the
//   outgoing eigenvector at each hyperbolic boundary, closing the 34, going
//   square, and recovering the manufactured solution.  TEST_HYP_MANUAL_CLOSURE=1
//   instead supplies that outgoing-characteristic closure by hand (see below) as
//   a cross-check oracle, and disables the framework closure (so the two never
//   double-close) -- the reference for what the auto-closure must reproduce.
// The driver is written to be easy to retune from the square/rank report.
//
// Build flags:  auto-closure and the incoming-BC guard are now ON by default in
//   the header (opt out with -DMC__OCFESLV_NO_AUTO_HYP_CLOSURE / _NO_HYP_BC_GUARD);
//   V2 interface decision is also the default.
//   + -DCRONOS__WITH_SPQR for QR, + -DCRONOS__WITH_UMFPACK for IC_STRONG.  Recommended:
//   -DMC__OCFESLV_SYMBOL_BALANCE_PROBE  to print the per-direction read table
//   (watch the (c,z)/(u,z) rows -- a per-state read shows order-1 with no
//   coherent characteristic direction).
//
// Knobs (-D overrides):
//   -DTEST_HYP_NEL_T / _NEL_Z      finite elements    (default 3 / 3)
//   -DTEST_HYP_NT   / _NZ          nodes per element  (default 6 / 10)
//   -DTEST_HYP_A                   off-diagonal coupling / wave speed (default 1)
//   -DTEST_HYP_K  -DTEST_HYP_TF  -DTEST_HYP_ZF
// ===========================================================================

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <map>
#include <string>
#include <vector>

#include <armadillo>

// Closure source.  The framework auto-closure is the default (enabled in the
// header).  TEST_HYP_MANUAL_CLOSURE selects the driver's hand-written cross-check
// oracle instead; when on, we suppress the framework closure here BEFORE the
// header is included, so the two never both emit rows (double closure).
// Default 0 = use the validated framework auto-closure.
#ifndef TEST_HYP_MANUAL_CLOSURE
#define TEST_HYP_MANUAL_CLOSURE 0
#endif
#if TEST_HYP_MANUAL_CLOSURE && !defined(MC__OCFESLV_NO_AUTO_HYP_CLOSURE)
#define MC__OCFESLV_NO_AUTO_HYP_CLOSURE
#endif

#define TEST_HYP_SPQR

#ifndef OCFE_OCFESLV_HEADER
#define OCFE_OCFESLV_HEADER "ocfeslv.hpp"
#endif
#include OCFE_OCFESLV_HEADER

using namespace mc;

#ifndef TEST_HYP_NEL_T
#define TEST_HYP_NEL_T 3
#endif
#ifndef TEST_HYP_NEL_Z
#define TEST_HYP_NEL_Z 3
#endif
#ifndef TEST_HYP_NT
#define TEST_HYP_NT 6
#endif
#ifndef TEST_HYP_NZ
#define TEST_HYP_NZ 10
#endif
#ifndef TEST_HYP_A
#define TEST_HYP_A 1.0     // off-diagonal coupling = wave speed (characteristics +-a)
#endif
#ifndef TEST_HYP_TF
#define TEST_HYP_TF 0.5
#endif
#ifndef TEST_HYP_ZF
#define TEST_HYP_ZF 1.0
#endif
#ifndef TEST_HYP_MAXIT
#define TEST_HYP_MAXIT 40
#endif
#ifndef TEST_HYP_SOLVE_TOL
#define TEST_HYP_SOLVE_TOL 1e-9
#endif
#ifndef TEST_HYP_EXACT_TOL
#define TEST_HYP_EXACT_TOL 1e-6
#endif
// (TEST_HYP_MANUAL_CLOSURE is resolved before the header include above, where it
// also disables the framework auto-closure to avoid double closure.  Default 0 =
// framework auto-closure; =1 selects the hand-written manual oracle.)

static double const kPi = 3.14159265358979323846;

struct Par {
  double tf  = TEST_HYP_TF;
  double zf  = TEST_HYP_ZF;
  double a   = TEST_HYP_A;    // coupling / wave speed
  double k   = 2.0*kPi;       // wavenumber
  double Ac  = 1.0;           // amplitude of c
  double Au  = 0.6;           // amplitude of u (distinct so c+-u are nontrivial)
  double pc  = 0.7;           // phase of c
  double pu  = 0.3;           // phase of u
  double alp = 0.3;           // linear-in-time growth
};

// Manufactured solution:
//   c = Ac sin(k z + pc)(1 + alp t),  u = Au sin(k z + pu)(1 + alp t).
static double C_exact( double t, double z, Par const& p ){ return p.Ac*std::sin(p.k*z+p.pc)*(1.0+p.alp*t); }
static double U_exact( double t, double z, Par const& p ){ return p.Au*std::sin(p.k*z+p.pu)*(1.0+p.alp*t); }

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
// (duplicate) collocation nodes -- a non-zero spread at an interior interface
// is a continuity jump.  For this hyperbolic system the physically meaningful
// continuity is on c+-u, but a per-state jump in c or u still signals trouble.
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
                    bool square=false, solved=false, rejected_bc=false;
                    double final_res=0., eC=0., eU=0.; };

static ModeResult run_mode( OCFESLV::Options::ImpositionType imp,
                            std::string const& suffix, Par const& p,
                            bool force_fwd = false )
{
  (void)suffix;
  ModeResult R; R.name = imp_name(imp); bool ok=true;

  // Propagation/inflow geometry from the sign of the wave speed a:
  //   a>0: w+ (speed +a) enters at z=LB, w- (speed -a) enters at z=UB.
  //   a<0: it flips -- w- enters at z=LB, w+ enters at z=UB.
  // force_fwd pins the a>0 convention regardless of sign: a deliberately
  // MIS-PLACED control that prescribes the OUTGOING characteristic at each end,
  // leaving the incoming one unconstrained -> square by count but rank-deficient
  // -> wrong (null-space) answer.  It demonstrates the answer tracks the sign.
  bool const fwd = force_fwd ? true : ( p.a >= 0.0 );

  std::cout << "\n========== off-diagonal first-order hyperbolic test ==========\n";
  std::cout << "imposition: " << R.name
            << ", finite elements: t=" << TEST_HYP_NEL_T << " z=" << TEST_HYP_NEL_Z
            << ", nodes/element: t=" << TEST_HYP_NT << " z=" << TEST_HYP_NZ << "\n";
  std::cout << "PDE: d_t c + a d_z u = f_c ; d_t u + a d_z c = f_u,  a=" << p.a
            << " (characteristics c+-u at speeds +-a),  k=" << p.k << "\n";
  std::cout << "placement: " << (fwd?"FORWARD (w+ in@LB, w- in@UB)":"REVERSE (w- in@LB, w+ in@UB)")
            << (force_fwd && p.a<0.0 ? "  [MIS-PLACED control: expect wrong answer]" : "")
            << "\n";

  FFGraph DAG;
  FFVar t = DAG.add_var("t");
  FFVar z = DAG.add_var("z");
  FFVar c = DAG.add_var("c(t,z)");
  FFVar u = DAG.add_var("u(t,z)");
  FFPartial OpP;

  // Manufactured solution, forcing and characteristic-BC data in CLOSED FORM.
  // OpP is applied only to states (c,u); the manufactured expressions are never
  // differentiated with OpP (that would create a state-free derivative-equation
  // the domain-consistency check rejects).
  //   f_c = d_t c + a d_z u = Ac alp sin(.c) + a Au k cos(.u)(1+alp t)
  //   f_u = d_t u + a d_z c = Au alp sin(.u) + a Ac k cos(.c)(1+alp t)
  FFVar sinc = sin( p.k*z + p.pc );
  FFVar cosc = cos( p.k*z + p.pc );
  FFVar sinu = sin( p.k*z + p.pu );
  FFVar cosu = cos( p.k*z + p.pu );
  FFVar CE = p.Ac*sinc*( 1.0 + p.alp*t );
  FFVar UE = p.Au*sinu*( 1.0 + p.alp*t );
  FFVar FC = p.Ac*p.alp*sinc + p.a*p.Au*p.k*cosu*( 1.0 + p.alp*t );
  FFVar FU = p.Au*p.alp*sinu + p.a*p.Ac*p.k*cosc*( 1.0 + p.alp*t );

  // Incoming-characteristic boundary data:
  //   w+ = c+u,  w- = c-u, each evaluated at both ends (z=LB, z=UB).
  //   Which one is prescribed (incoming) vs propagated (outgoing) at each end
  //   is set by the sign of a (see fwd above).
  // All four end/characteristic values are needed once placement can flip:
  double const wpLB0 = p.Ac*std::sin(p.pc)          + p.Au*std::sin(p.pu);          // w+ at z=LB
  double const wmLB0 = p.Ac*std::sin(p.pc)          - p.Au*std::sin(p.pu);          // w- at z=LB
  double const wpUB0 = p.Ac*std::sin(p.k*p.zf+p.pc) + p.Au*std::sin(p.k*p.zf+p.pu); // w+ at z=UB
  double const wmUB0 = p.Ac*std::sin(p.k*p.zf+p.pc) - p.Au*std::sin(p.k*p.zf+p.pu); // w- at z=UB
  FFVar wpLB = wpLB0*( 1.0 + p.alp*t );
  FFVar wmLB = wmLB0*( 1.0 + p.alp*t );
  FFVar wpUB = wpUB0*( 1.0 + p.alp*t );
  FFVar wmUB = wmUB0*( 1.0 + p.alp*t );

  FFVar PDEC = OpP(c,t) + p.a*OpP(u,z) - FC;   // d_t c + a d_z u - f_c
  FFVar PDEU = OpP(u,t) + p.a*OpP(c,z) - FU;   // d_t u + a d_z c - f_u
  FFVar ICC  = c - CE;                         // initial c(0,z)
  FFVar ICU  = u - UE;                         // initial u(0,z)

  // Characteristic combinations:
  //   w+ = c+u transports at speed +a ;  w- = c-u transports at speed -a.
  //   transport eqns:  WPz = (PDEC+PDEU) = d_t w+ + a d_z w+ - (f_c+f_u)
  //                    WMz = (PDEC-PDEU) = d_t w- - a d_z w- - (f_c-f_u)
  FFVar WPz = PDEC + PDEU;                      // w+ transport
  FFVar WMz = PDEC - PDEU;                      // w- transport

  // Sign-aware placement: at each end prescribe the INCOMING characteristic
  // (its value from the manufactured solution) and collocate the OUTGOING
  // characteristic's transport.  Reverse-sign swaps which is which.
  FFVar BC_LB  = fwd ? ( c + u - wpLB ) : ( c - u - wmLB );  // incoming @LB: w+ (fwd) / w- (rev)
  FFVar CLO_LB = fwd ? WMz              : WPz;               // outgoing @LB: w- (fwd) / w+ (rev)
  FFVar BC_UB  = fwd ? ( c - u - wmUB ) : ( c + u - wpUB );  // incoming @UB: w- (fwd) / w+ (rev)
  FFVar CLO_UB = fwd ? WPz              : WMz;               // outgoing @UB: w+ (fwd) / w- (rev)

  OCFESLV oc(&DAG);
  oc.add_domain( t, FFDom(0., p.tf, TEST_HYP_NEL_T, FFDom::LGR, TEST_HYP_NT) );
  oc.add_domain( z, FFDom(0., p.zf, TEST_HYP_NEL_Z, FFDom::LGL, TEST_HYP_NZ) );
  oc.add_state( c, {t,z} );
  oc.add_state( u, {t,z} );
  oc.update_ref( c, [&]( OCFESLV::t_Coord const& crd ){ return C_exact(crd.at(t),crd.at(z),p); } );
  oc.update_ref( u, [&]( OCFESLV::t_Coord const& crd ){ return U_exact(crd.at(t),crd.at(z),p); } );

  OCFESLV::EqnOptions int_opt( OCFESLV::EqnRole::INTERIOR, 0 );
  OCFESLV::EqnOptions ini_opt( OCFESLV::EqnRole::INITIAL,  0 );
  OCFESLV::EqnOptions bnd_opt( OCFESLV::EqnRole::BOUNDARY, 0 );
  // Outgoing-characteristic closure rows: interior-operator collocations (so
  // they receive SAT / collocate like the PDE) but EXCLUDED from the principal
  // symbol (classify=false).  They are not governing-evolution equations, so
  // leaving them in classification makes the block symbol non-square (nEqn>nState)
  // and classify() bails to UNDETERMINED before the eigenvalue analysis.
  OCFESLV::EqnOptions clo_opt( OCFESLV::EqnRole::INTERIOR, 0, OCFESLV::Options::IC_AUTO,
                             /*classify=*/false, /*sat=*/true );

  int const T_INT = FFDom::ALL - FFDom::LB;
  int const Z_INT = FFDom::ALL - FFDom::LB - FFDom::UB;

  // Two first-order PDEs on the z-interior, t>0.  Initial data for both states.
  // One incoming-characteristic condition at each end, and (default) the
  // matching outgoing-characteristic transport.  Without the outgoing closure
  // the count is short by 2*(t-interior nodes); with it mis-placed (force_fwd
  // for a<0) the count is right but the rank is not -> wrong answer.
  oc.add_equation( PDEC,  {t,z}, {T_INT, Z_INT},         int_opt );
  oc.add_equation( PDEU,  {t,z}, {T_INT, Z_INT},         int_opt );
  oc.add_equation( ICC,   {t,z}, {FFDom::LB, FFDom::ALL}, ini_opt );
  oc.add_equation( ICU,   {t,z}, {FFDom::LB, FFDom::ALL}, ini_opt );
  oc.add_equation( BC_LB, {t,z}, {T_INT, FFDom::LB},     bnd_opt );  // incoming characteristic @LB
  oc.add_equation( BC_UB, {t,z}, {T_INT, FFDom::UB},     bnd_opt );  // incoming characteristic @UB
#if TEST_HYP_MANUAL_CLOSURE
  oc.add_equation( CLO_LB, {t,z}, {T_INT, FFDom::LB},    clo_opt );  // outgoing characteristic @LB
  oc.add_equation( CLO_UB, {t,z}, {T_INT, FFDom::UB},    clo_opt );  // outgoing characteristic @UB
#else
  (void)CLO_LB; (void)CLO_UB;
#endif

  oc.set_evolution_domain( t );
  oc.options.REDUCE.ORDER    = OCFESLV::Options::RED_FULL;   // inert: no high-order spatial
  oc.options.CLASSIFY.MODE        = OCFESLV::Options::CLASS_AUTO;
  oc.options.INTERFACE.TYPE  = OCFESLV::Options::IC_AUTO;
  oc.options.INTERFACE.IMPOSITION = imp;
  oc.options.INTERFACE.SAT_SIGMA0      = 10.0;
  oc.options.DISPLAY_LEVEL   = 1;
#if TEST_HYP_MANUAL_CLOSURE
  // rev317/rev318: the MANUAL ORACLE compares a hand-written closure against the automatic one, so the automatic
  // closure must be OFF here.  (It used to be disabled through
  // MC__OCFESLV_NO_AUTO_HYP_CLOSURE, which no header has honoured for some time: this variant had been
  // double-closed and failing until rev317.)
  // An OPTION again since 2026-10-07 (AUTO.HYP_CLOSURE, WORKPLAN 3.C); the environment variable is retired.
  oc.options.AUTO.HYP_CLOSURE = false;
#endif

  if( !oc.setup() ){
    std::cerr << "ERROR: OCFESLV::setup() failed for " << R.name
              << ": " << OCFESLV::setup_status_str( oc.setup_status() ) << "\n";
    R.rejected_bc = ( oc.setup_status() == OCFESLV::SetupStatus::HYP_INCOMING_BC
                   || oc.setup_status() == OCFESLV::SetupStatus::HYP_BC_MISDIRECTED );
    R.ok=false; return R;
  }

  size_t const nVar=oc.n_colloc_sta();
  size_t const nEqn=oc.n_colloc_eqn();
  size_t const nTrace=oc.n_colloc_trace();
  R.nVar=nVar; R.nEqn=nEqn; R.nTrace=nTrace; R.square=(nVar==nEqn);
  std::cout << "nVar=" << nVar << " nEqn=" << nEqn << " nTrace=" << nTrace
            << " square=" << (R.square?"yes":"no") << "\n";
  if( !R.square )
    std::cout << "  [count] nEqn-nVar=" << (long long)nEqn-(long long)nVar
              << " (without outgoing closure: short by 2*(t-interior) = one missing"
                 " condition per (boundary end x t>0 node))\n";
  ok &= check_close("square system (nVar==nEqn)", R.square?0.0:1.0, 0.0);

  // Diagnostics: no reduce_order auxes expected (first-order); the block
  // classification (hyperbolic?) and the RESOLVED z interface read -- a coupled
  // characteristic read should give UPWIND on z, a per-state read has no
  // coherent direction (zero diagonal speed).
  std::cout << "States after setup:";
  for( auto const& st: oc.states_colloc() ) std::cout << ' ' << st.name();
  std::cout << "\nAuxiliary states introduced: "
            << (oc.states_colloc().size()>2?oc.states_colloc().size()-2:0)
            << " (expect 0: first-order system, no high-order spatial)\n";
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

  // Initial guess: zero (perturbed from exact) -- the system is LINEAR, so a
  // square non-singular setup recovers the manufactured solution in one Newton
  // step; failure to recover it flags an ill-posed/under-determined setup.
  std::vector<double> var(nVar,0.0);

  // This driver hand is for monolithic only
  oc.options.SOLVE.MARCHING = false;
  oc.options.SOLVE.MAX_ITER = TEST_HYP_MAXIT;
  oc.options.SOLVE.RES_TOL  = TEST_HYP_SOLVE_TOL;
#if defined(CRONOS__WITH_SPQR) && defined(TEST_HYP_SPQR)
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
  ok &= check_close("final max residual", R.final_res, TEST_HYP_SOLVE_TOL);

  // Per-state accuracy against the manufactured solution.
  {
    double eC=0.0, eU=0.0; size_t off=0, sidx=0;
    for( auto const& st: oc.states_colloc() ){
      auto nodes = oc.node_colloc(st);
      for( size_t i=0; i<nodes.size(); ++i ){
        double const ex = (sidx==0)? C_exact(nodes[i][0],nodes[i][1],p)
                                   : U_exact(nodes[i][0],nodes[i][1],p);
        double const e = std::fabs( var[off+i] - ex );
        if( sidx==0 ) eC = std::max(eC,e); else eU = std::max(eU,e);
      }
      off += nodes.size(); ++sidx;
    }
    R.eC=eC; R.eU=eU;
    ok &= check_close("max |c - c_exact|", R.eC, TEST_HYP_EXACT_TOL);
    ok &= check_close("max |u - u_exact|", R.eU, TEST_HYP_EXACT_TOL);
  }

  print_duplicate_spreads( oc, var );

  std::cout << "off-diagonal hyperbolic test (" << R.name << "): " << (ok?"PASS":"FAIL") << "\n";
  R.ok = ok;
  return R;
}

int main()
{
  double const amag = std::fabs( (double)TEST_HYP_A );
  Par pf; pf.a = +amag;   // forward propagation (w+ rightward)
  Par pr; pr.a = -amag;   // reverse propagation (w+ leftward)

  auto run3 = []( Par const& p, bool force_fwd, std::vector<ModeResult>& out ){
    out.push_back( run_mode(OCFESLV::Options::IC_WEAK,   "weak",   p, force_fwd) );
    out.push_back( run_mode(OCFESLV::Options::IC_TRACE,  "trace",  p, force_fwd) );
    out.push_back( run_mode(OCFESLV::Options::IC_STRONG, "strong", p, force_fwd) );
  };

  // A: forward sign, sign-correct placement.   Expect PASS (baseline).
  // B: reverse sign, sign-correct placement.   Expect PASS (closure follows sign).
  // C: reverse sign, MIS-PLACED (forward) BCs.  Expect WRONG ANSWER (square but
  //    rank-deficient) -- the demonstration that the answer tracks the sign.
  std::vector<ModeResult> A, B, C;
  run3( pf, /*force_fwd=*/false, A );
  run3( pr, /*force_fwd=*/false, B );
  C.push_back( run_mode(OCFESLV::Options::IC_WEAK, "weak", pr, /*force_fwd=*/true) );

  auto print_block = []( char const* title, std::vector<ModeResult> const& v ){
    std::cout << "\n============== " << title << " ==============\n";
    std::cout << std::left << std::setw(12) << "mode"
              << std::setw(8)  << "square" << std::setw(10) << "nTrace"
              << std::setw(15) << "final|r|" << std::setw(13) << "|c-exact|"
              << std::setw(13) << "|u-exact|" << "result\n";
    for( auto const& r: v )
      std::cout << std::left << std::setw(12) << r.name
                << std::setw(8)  << (r.square?"yes":"no")
                << std::setw(10) << r.nTrace
                << std::scientific << std::setprecision(3)
                << std::setw(15) << r.final_res
                << std::setw(13) << r.eC
                << std::setw(13) << r.eU
                << (r.rejected_bc ? "REJECTED (expected)" : (r.ok?"PASS":"FAIL")) << "\n";
  };

  print_block( "A: forward sign, correct placement (a=+|a|)", A );
  print_block( "B: reverse sign, correct placement (a=-|a|)", B );
  print_block( "C: reverse sign, MIS-PLACED control (expect WRONG)", C );

  bool all_ok = true;
  for( auto const& r: A ) all_ok &= r.ok;   // baseline must pass
  for( auto const& r: B ) all_ok &= r.ok;   // sign-correct reverse must pass

  // C is a control: it must NOT recover the solution.  Two correct outcomes:
  //   * guard ON (default): the mis-placed inflow is REJECTED at setup
  //     (HYP_BC_MISDIRECTED for a direction reject, or HYP_INCOMING_BC for a
  //     count mismatch) -- the framework refuses the ill-posed placement.
  //   * guard OFF (-DMC__OCFESLV_NO_HYP_BC_GUARD): setup stays square but the
  //     answer is visibly wrong, proving it depends on the characteristic sign
  //     and not merely on the equation count.
  bool const demo_ok = !C.empty() &&
       ( C[0].rejected_bc || ( C[0].square && C[0].eC > 1e-3 ) );
  std::cout << "\nDirection-matters demonstration: reverse-sign system, forward-placed BCs -> "
            << ( C.empty() ? "?"
               : C[0].rejected_bc ? "REJECTED at setup (mis-directed inflow)"
               : C[0].square ? "square but wrong" : "under-determined" )
            << "  |c-exact|=" << std::scientific << std::setprecision(3)
            << (C.empty()?0.0:C[0].eC)
            << "  => " << (demo_ok?"direction is load-bearing (as expected)":"inconclusive")
            << "\n";
  all_ok &= demo_ok;

  std::cout << "\n====================================================================\n";
  std::cout << "off-diagonal hyperbolic oracle (sign sweep + placement control): "
            << (all_ok?"PASS":"FAIL") << "\n";
  return all_ok ? 0 : 1;
}
