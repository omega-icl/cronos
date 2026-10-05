// ===========================================================================
// OCFE_PDE15_solve2.cpp
//
// NON-SYMMETRIC off-diagonal first-order HYPERBOLIC oracle -- stage 3 of the
// PSA case study, the symmetry-breaking sibling of PDE14.
//
// PDE14 used the SYMMETRIC symbol a*[[0,1],[1,0]]: there the left and right
// characteristic eigenvectors coincide, so the distinction between "which
// equation combination transports the outgoing characteristic" (a ROW-space /
// LEFT-eigenvector question) and "which state direction the incoming BC may
// couple to" (a COLUMN-space / RIGHT-eigenvector question) is invisible -- both
// answers are c+-u.  PSA's (pressure,velocity) continuity/momentum block is
// NOT symmetric, so this driver breaks the symmetry to exercise exactly that
// distinction before the framework meets it on the real model.
//
//   d_t c +     d_z u = f_c(t,z)        on (t,z)
//   d_t u + 4   d_z c = f_u(t,z)
//
//   z principal symbol  A_z = a*[[0,1],[4,0]]  (off-diagonals 1:4, a=speed scale)
//   eigenvalues  +-2a,
//   RIGHT eigenvectors  r_+ = (1, 2),  r_- = (1,-2)     (column space of A_z)
//   LEFT  eigenvectors  l_+ = (2, 1),  l_- = (2,-1)     (row space of A_z)
//   characteristic vars w_+ = 2c+u (speed +2a),  w_- = 2c-u (speed -2a),
//   travelling in OPPOSITE directions:
//       w_+ = 2c+u  incoming at z=LB  (prescribe there)
//       w_- = 2c-u  incoming at z=UB  (prescribe there)
//
//   The well-posed BC coupling is the LEFT eigenvector (the characteristic
//   combination l . (c,u)); the outgoing-characteristic CLOSURE is also a LEFT-
//   eigenvector combination of the governing equations (= rowspace of A_z@end).
//   Because left != right here, this driver is the first that can tell a
//   rowspace-based closure/guard apart from a colspace-based one.
//
// WHAT PDE15 VALIDATES (two framework paths PDE14's symmetry could not test):
//
//  1. AUTO-CLOSURE on a non-symmetric symbol.  The generator builds each
//     outgoing-characteristic row from rowspace(A_z) at the end (V-columns of
//     the SVD), i.e. the LEFT eigenvector -- here 2*PDEC-PDEU at z=LB (w_-) and
//     2*PDEC+PDEU at z=UB (w_+).  With the default build (auto-closure ON),
//     sections A/B PASS, confirming the
//     closure combines equations by the LEFT eigenvector, not the right.  (On
//     PDE14 a colspace closure would have passed too; here only rowspace does.)
//
//  2. INCOMING-BC GUARD direction test.  The guard checks the BC's algebraic
//     coupling d(BC)/d(state) lies in the INCOMING characteristic subspace.
//     That subspace is spanned by the LEFT eigenvector (l_+ = (2,1) at LB).  A
//     guard that projects onto colspace(A_+) = U-columns = the RIGHT eigenvector
//     (1,2) will WRONGLY REJECT the correct coupling (2,1), because (2,1) is not
//     in span{(1,2)}.  With the guard ON (default), if A/B were REJECTED
//     (setup_status HYP_INCOMING_BC) this driver had caught the colspace/rowspace
//     bug -- the fix (now in) projects onto rowspace (V-columns), matching the
//     closure.  With the fix, A/B pass with the guard on and C (mis-placed) is
//     correctly rejected.
//
// GAP-FINDER (TEST_HYP_MANUAL_CLOSURE=0, no auto): identical structural gap to
//   PDE14 -- a first-order block has no reduce_order LINK, the PDEs collocate on
//   the z-interior, so the OUTGOING characteristic at each end is left
//   unequationed and nEqn-nVar = -2*(t-interior nodes).  Continuity is fine; the
//   null space is unconstrained and the solve lands on a wrong solution.
//
// MANUAL ORACLE (TEST_HYP_MANUAL_CLOSURE=1, which auto-disables the framework
//   closure before the header include): the driver supplies the LEFT-eigenvector
//   outgoing closure by hand, the system goes square, and the manufactured
//   solution is recovered across all three IC modes for both wave-speed signs --
//   the cross-check that the framework auto-closure (the default) reproduces.
//
// Recommended probe: -DMC__OCFESLV_SYMBOL_BALANCE_PROBE prints the per-direction
//   read table and the A_+/A_- projectors -- for this NON-symmetric symbol the
//   projectors are non-symmetric (rank-1 onto the eigenspaces), visibly unlike
//   PDE14's symmetric 0.5[[1,1],[1,1]].
//
// Build flags:  (V2 interface decision is the default in the header now)
//   + -DCRONOS__WITH_SPQR for QR, + -DCRONOS__WITH_UMFPACK for IC_STRONG.
//
// Knobs (-D overrides):
//   -DTEST_HYP_NEL_T / _NEL_Z      finite elements    (default 3 / 3)
//   -DTEST_HYP_NT   / _NZ          nodes per element  (default 6 / 10)
//   -DTEST_HYP_A                   signed wave-speed scale (eigenvalues +-2A) (default 1)
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
// oracle instead; when it is on we suppress the framework closure here, BEFORE the
// header is included, so the two never both emit rows (double closure -> over-
// determined).  Default 0 = use the validated framework auto-closure.
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
#define TEST_HYP_A 1.0     // signed wave-speed scale; eigenvalues +-2A (NON-symmetric symbol)
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
// Manual characteristic boundary closure.  DEFAULT ON: supplies the OUTGOING-
// characteristic transport equation at each z-end (w- at z=LB, w+ at z=UB) as a
// LEFT-eigenvector combination of the governing PDEs, so the system goes square
// and recovers the manufactured solution.  Build with -DTEST_HYP_MANUAL_CLOSURE=0
// to reproduce the gap-finder result: the outgoing characteristic is left
// unequationed at each boundary node, giving nEqn-nVar = -2*(t-interior nodes)
// (= -34 on the default 3x3 / 6x10 grid) and an under-determined solve.
// (TEST_HYP_MANUAL_CLOSURE is resolved before the header include above, where it
// also disables the framework auto-closure to avoid double closure.)

static double const kPi = 3.14159265358979323846;

struct Par {
  double tf  = TEST_HYP_TF;
  double zf  = TEST_HYP_ZF;
  double a   = TEST_HYP_A;    // signed wave-speed scale (eigenvalues +-2a)
  double g   = 4.0;           // ratio coeff: A_z = a[[0,1],[g,0]]; g!=1 => NON-SYMMETRIC
  double sg  = 2.0;           // sqrt(g): left eigvecs (sg,+-1), right eigvecs (1,+-sg)
  double k   = 2.0*kPi;       // wavenumber
  double Ac  = 1.0;           // amplitude of c
  double Au  = 0.6;           // amplitude of u (distinct so sg*c+-u are nontrivial)
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
// continuity is on the characteristics sg*c+-u, but a per-state jump in c or u
// still signals trouble.
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

  // Propagation/inflow geometry from the sign of the wave-speed scale a:
  //   a>0: w+ (speed +2a) enters at z=LB, w- (speed -2a) enters at z=UB.
  //   a<0: it flips -- w- enters at z=LB, w+ enters at z=UB.
  // force_fwd pins the a>0 convention regardless of sign: a deliberately
  // MIS-PLACED control that prescribes the OUTGOING characteristic at each end,
  // leaving the incoming one unconstrained -> square by count but rank-deficient
  // -> wrong (null-space) answer.  It demonstrates the answer tracks the sign.
  bool const fwd = force_fwd ? true : ( p.a >= 0.0 );

  std::cout << "\n====== NON-SYMMETRIC off-diagonal first-order hyperbolic test ======\n";
  std::cout << "imposition: " << R.name
            << ", finite elements: t=" << TEST_HYP_NEL_T << " z=" << TEST_HYP_NEL_Z
            << ", nodes/element: t=" << TEST_HYP_NT << " z=" << TEST_HYP_NZ << "\n";
  std::cout << "PDE: d_t c + a d_z u = f_c ; d_t u + g a d_z c = f_u,  a=" << p.a
            << " g=" << p.g << " (eigenvalues +-" << p.sg << "a, characteristics "
            << p.sg << "c+-u),  k=" << p.k << "\n";
  std::cout << "left eigvecs (" << p.sg << ",+-1) != right eigvecs (1,+-" << p.sg
            << ")  [NON-SYMMETRIC: closure & BC couple via LEFT/rowspace]\n";
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
  //   f_c = d_t c +   a d_z u = Ac alp sin(.c) +   a Au k cos(.u)(1+alp t)
  //   f_u = d_t u + g a d_z c = Au alp sin(.u) + g a Ac k cos(.c)(1+alp t)
  FFVar sinc = sin( p.k*z + p.pc );
  FFVar cosc = cos( p.k*z + p.pc );
  FFVar sinu = sin( p.k*z + p.pu );
  FFVar cosu = cos( p.k*z + p.pu );
  FFVar CE = p.Ac*sinc*( 1.0 + p.alp*t );
  FFVar UE = p.Au*sinu*( 1.0 + p.alp*t );
  FFVar FC = p.Ac*p.alp*sinc +     p.a*p.Au*p.k*cosu*( 1.0 + p.alp*t );
  FFVar FU = p.Au*p.alp*sinu + p.g*p.a*p.Ac*p.k*cosc*( 1.0 + p.alp*t );

  // Incoming-characteristic boundary data:
  //   w+ = sg c + u,  w- = sg c - u, each evaluated at both ends (z=LB, z=UB).
  //   Which one is prescribed (incoming) vs propagated (outgoing) at each end
  //   is set by the sign of a (see fwd above).
  // All four end/characteristic values are needed once placement can flip:
  double const wpLB0 = p.sg*p.Ac*std::sin(p.pc)          + p.Au*std::sin(p.pu);          // w+ at z=LB
  double const wmLB0 = p.sg*p.Ac*std::sin(p.pc)          - p.Au*std::sin(p.pu);          // w- at z=LB
  double const wpUB0 = p.sg*p.Ac*std::sin(p.k*p.zf+p.pc) + p.Au*std::sin(p.k*p.zf+p.pu); // w+ at z=UB
  double const wmUB0 = p.sg*p.Ac*std::sin(p.k*p.zf+p.pc) - p.Au*std::sin(p.k*p.zf+p.pu); // w- at z=UB
  FFVar wpLB = wpLB0*( 1.0 + p.alp*t );
  FFVar wmLB = wmLB0*( 1.0 + p.alp*t );
  FFVar wpUB = wpUB0*( 1.0 + p.alp*t );
  FFVar wmUB = wmUB0*( 1.0 + p.alp*t );

  FFVar PDEC = OpP(c,t) +     p.a*OpP(u,z) - FC;   // d_t c +   a d_z u - f_c
  FFVar PDEU = OpP(u,t) + p.g*p.a*OpP(c,z) - FU;   // d_t u + g a d_z c - f_u
  FFVar ICC  = c - CE;                             // initial c(0,z)
  FFVar ICU  = u - UE;                             // initial u(0,z)

  // Characteristic combinations (LEFT eigenvectors l_+ = (sg,1), l_- = (sg,-1)):
  //   w+ = sg c + u transports at speed +sg a ;  w- = sg c - u at speed -sg a.
  //   transport eqns:  WPz = (sg PDEC + PDEU) = d_t w+ + sg a d_z w+ - (sg f_c+f_u)
  //                    WMz = (sg PDEC - PDEU) = d_t w- - sg a d_z w- - (sg f_c-f_u)
  FFVar WPz = p.sg*PDEC + PDEU;                    // w+ transport (rowspace l_+)
  FFVar WMz = p.sg*PDEC - PDEU;                    // w- transport (rowspace l_-)

  // Sign-aware placement: at each end prescribe the INCOMING characteristic
  // (its value from the manufactured solution, coupled via the LEFT eigenvector)
  // and collocate the OUTGOING characteristic's transport.  Reverse-sign swaps
  // which is which.
  FFVar BC_LB  = fwd ? ( p.sg*c + u - wpLB ) : ( p.sg*c - u - wmLB );  // incoming @LB: w+ (fwd) / w- (rev)
  FFVar CLO_LB = fwd ? WMz                   : WPz;                    // outgoing @LB: w- (fwd) / w+ (rev)
  FFVar BC_UB  = fwd ? ( p.sg*c - u - wmUB ) : ( p.sg*c + u - wpUB );  // incoming @UB: w- (fwd) / w+ (rev)
  FFVar CLO_UB = fwd ? WPz                   : WMz;                    // outgoing @UB: w+ (fwd) / w- (rev)

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
  // closure must be OFF here.  Environment-only since rev318, hence setenv.  (It used to be disabled through
  // MC__OCFESLV_NO_AUTO_HYP_CLOSURE, which no header has honoured for some time: this variant had been
  // double-closed and failing until rev317.)
  setenv( "CRONOS_AUTO_HYP_CLOSURE", "0", 1 );
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

  std::cout << "non-symmetric hyperbolic test (" << R.name << "): " << (ok?"PASS":"FAIL") << "\n";
  R.ok = ok;
  return R;
}

int main()
{
  // Self-documenting build config -- removes the "which flags did I compile
  // with?" ambiguity from the log (closure source + guard state).
  std::cout << "PDE15 build config: closure="
#if TEST_HYP_MANUAL_CLOSURE && defined(MC__OCFESLV_AUTO_HYP_CLOSURE)
            << "MANUAL+AUTO(double!)"
#elif TEST_HYP_MANUAL_CLOSURE
            << "MANUAL"
#elif defined(MC__OCFESLV_AUTO_HYP_CLOSURE)
            << "AUTO"
#else
            << "NONE(gap-finder)"
#endif
            << "  guard="
            << "ON"
            << "\n";

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
  //    (With -DMC__OCFESLV_HYP_BC_GUARD this run is instead correctly REJECTED.)
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
  print_block( "C: reverse sign, MIS-PLACED control (expect WRONG / rejected)", C );

  bool all_ok = true;
  for( auto const& r: A ) all_ok &= r.ok;   // baseline must pass
  for( auto const& r: B ) all_ok &= r.ok;   // sign-correct reverse must pass

  // C is a control: it must NOT recover the solution.  Two correct outcomes,
  // depending on the guard:
  //   * guard ON (default): the mis-placed inflow is REJECTED at setup
  //     (HYP_INCOMING_BC) -- the framework refuses the ill-posed placement.
  //   * guard OFF (-DMC__OCFESLV_NO_HYP_BC_GUARD): setup is square but the answer
  //     is wrong (~1e-4, 1000s x the A/B machine accuracy) because the auto-
  //     closure is correct and only the mis-placed BC corrupts the solution.
  // Either way, direction is load-bearing.
  bool const demo_ok = !C.empty() &&
       ( C[0].rejected_bc || ( C[0].square && C[0].eC > 1.0e-5 ) );
  std::cout << "\nDirection-matters demonstration: reverse-sign system, forward-placed BCs -> "
            << ( C.empty() ? "?"
               : C[0].rejected_bc ? "REJECTED at setup (HYP_INCOMING_BC)"
               : C[0].square ? "square but wrong" : "under-determined" )
            << "  |c-exact|=" << std::scientific << std::setprecision(3)
            << (C.empty()?0.0:C[0].eC)
            << "  => " << (demo_ok?"direction is load-bearing (as expected)":"inconclusive")
            << "\n";
  all_ok &= demo_ok;

  std::cout << "\n====================================================================\n";
  std::cout << "non-symmetric hyperbolic oracle (sign sweep + placement control): "
            << (all_ok?"PASS":"FAIL") << "\n";
  return all_ok ? 0 : 1;
}
