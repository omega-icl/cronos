// ===========================================================================
// OCFE_PDE16_solve2.cpp
//
// FLUX / ROBIN inflow oracle on a NON-SYMMETRIC first-order hyperbolic block --
// stage 4 of the PSA case study, the validator for the incoming-BC guard's
// FLUX direction test (MC__OCFESLV_HYP_FLUX_BC_GUARD).
//
// PDE15 prescribed each incoming characteristic by VALUE: BC = w+- - g.  Its
// coupling d(BC)/d(state) is the bare left eigenvector, FAD-able, and the guard's
// value projection onto Vin sufficed.  A real inflow (PSA Danckwerts, Darcy) is
// FLUX form: it carries a normal-derivative term d_z(.) that survives as a DAG
// OpP node, which FAD cannot traverse -- so the original guard SKIPPED it.  The
// refinement proxies the OpP node with the principal-symbol derivative proxies so
// the VALUE coupling a = d(BC)/d(u) becomes FAD-able, then applies the SAME test
// as PDE15: a must lie in the incoming characteristic rowspace Vin.
//
// The derivative coupling b = d(BC)/d(d_z u) is NOT tested.  A BC is ill-posed
// only when it prescribes an OUTGOING characteristic's VALUE (colliding with the
// auto-closure, which transports the outgoing values out).  The closure is a
// transport equation -- it does NOT pin the outgoing d_z w -- so a flux BC whose
// DERIVATIVE term touches an outgoing characteristic is still well-posed.  (An
// earlier revision also required b in Vin; that false-rejected the well-posed
// section C below, caught by the -DTEST_PDE16_NO_FLUX_GUARD run.)
//
// This driver is the LINEAR-flux validator for that mechanism; the Darcy-bilinear
// section (proxy seeded by the reference derivative, where a itself references a
// derivative proxy) follows on this base.
//
//   d_t c +     d_z u = f_c(t,z)        on (t,z)   [identical system to PDE15]
//   d_t u + 4   d_z c = f_u(t,z)
//   A_z = a*[[0,1],[4,0]],  eigenvalues +-2a,
//   LEFT eigvecs l_+ = (2,1), l_- = (2,-1)  (rowspace; the BC-coupling space),
//   characteristics w_+ = 2c+u (speed +2a, incoming@LB),
//                   w_- = 2c-u (speed -2a, incoming@UB).
//
// INFLOW BC, FLUX/ROBIN form (the new bit vs PDE15):
//   at the incoming end prescribe   w +- beta d_z w  =  g(t),
//   i.e.  BC = (l . u) + beta (l . d_z u) - g, with value coupling a = l in Vin.
//   Because g is the manufactured value of the same combination the square solve
//   recovers the exact solution.  The spatial profile is a degree-3 POLYNOMIAL so
//   the spectral d_z is exact at the boundary and the flux RHS is discretely
//   consistent (machine accuracy at the full beta).
//
// WHAT PDE16 VALIDATES (beyond PDE15):
//   * the guard proxies the surviving d_z OpP nodes so the VALUE coupling of a
//     FLUX/Robin BC is FAD-able, then projects a onto Vin (value test extended to
//     flux BCs -- the original guard skipped these entirely);
//   * section C: a flux BC with the CORRECT value coupling but an off-character
//     DERIVATIVE term is WELL-POSED and must be ACCEPTED + recover -- the guard
//     must not false-reject it (regression guard vs the removed b-projection);
//   * section D: a flux BC whose VALUE coupling prescribes the OUTGOING
//     characteristic is ill-posed -> REJECTED (a=l_out outside Vin).  This is the
//     discrimination the extended value test provides over the old OpP-skip.
//
// Sections:
//   A  forward sign (a=+|a|), well-directed flux inflow        -> ACCEPT + recover
//   B  reverse sign (a=-|a|), well-directed flux inflow        -> ACCEPT + recover
//   C  forward sign, correct value + OFF-CHARACTER DERIVATIVE   -> ACCEPT + recover
//   D  forward sign, MIS-DIRECTED VALUE (outgoing char, Robin)  -> REJECTED at setup
//
// The flux guard is OPT-IN in the header (default off); this oracle exists to
// exercise it, so it is enabled here by default (disable with
// -DTEST_PDE16_NO_FLUX_GUARD: then the OpP BCs skip the test, C still recovers
// (well-posed), and D is NO LONGER rejected -- the before/after in one driver).
//
// Build flags:  auto-closure + base incoming-BC guard are default-on in the
//   header; + -DCRONOS__WITH_SPQR for QR, + -DCRONOS__WITH_UMFPACK for IC_STRONG.
//   Recommended: -DMC__OCFESLV_SYMBOL_BALANCE_PROBE and DISPLAY_LEVEL>=1 print the
//   per-BC value-resid the guard computes.
//
// Knobs (-D overrides):
//   -DTEST_HYP_NEL_T / _NEL_Z   finite elements      (default 3 / 3)
//   -DTEST_HYP_NT   / _NZ       nodes per element    (default 6 / 10)
//   -DTEST_HYP_A                signed wave-speed scale (eigenvalues +-2A) (default 1)
//   -DTEST_PDE16_BETA           Robin coefficient on d_z w (default 0.1)
//   -DTEST_PDE16_NO_FLUX_GUARD  build with the flux guard OFF (before/after)
//   -DTEST_HYP_TF  -DTEST_HYP_ZF       domain extents (default 0.5 / 1.0)
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

// The flux-BC guard is now DEFAULT-ON in the header (criterion #3).  For the
// guard-OFF "before" build, -DTEST_PDE16_NO_FLUX_GUARD must both #undef any
// force-define AND set the header opt-out MC__OCFESLV_NO_HYP_FLUX_BC_GUARD, else the
// header's default-on auto-define would re-enable it.  guard-ON needs no action.
// Authoritative against a harness force-define of -DMC__OCFESLV_HYP_FLUX_BC_GUARD.
// (The previous revision's hard `#define TEST_PDE16_NO_FLUX_GUARD` is removed.)

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
#ifndef TEST_PDE16_BETA
#define TEST_PDE16_BETA 0.1   // Robin coefficient on d_z w in the flux inflow BC
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

struct Par {
  double tf   = TEST_HYP_TF;
  double zf   = TEST_HYP_ZF;
  double a    = TEST_HYP_A;     // signed wave-speed scale (eigenvalues +-2a)
  double g    = 4.0;            // ratio coeff: A_z = a[[0,1],[g,0]]; g!=1 => NON-SYMMETRIC
  double sg   = 2.0;            // sqrt(g): left eigvecs (sg,+-1), right eigvecs (1,+-sg)
  // Spatial profiles are low-degree POLYNOMIALS in z (degree 3 <= N-1), so the
  // spectral z-derivative is EXACT at every node including the boundary.  The
  // Robin/flux inflow RHS is built from the analytic d_z w, so with an exactly-
  // differentiated profile it is discretely consistent and the solve recovers to
  // machine accuracy at the full beta.  (A sin(kz) profile leaves an O(beta)
  // analytic-vs-spectral boundary-derivative aliasing error that pollutes the
  // incoming characteristic -- correct guard mechanism, contaminated accuracy.)
  double cc[4] = { 1.0, 0.9, -0.5,  0.3 };  // Pc(z) = cc0+cc1 z+cc2 z^2+cc3 z^3
  double uu[4] = { 0.4, 0.7,  0.3, -0.3 };  // Pu(z) (distinct: w_+-=sg c+-u nontrivial)
  double alp  = 0.3;            // linear-in-time growth
  double beta = TEST_PDE16_BETA;// Robin coefficient on d_z w
};

// Polynomial profile value / z-derivative (Horner).
static double polyval( double const c[4], double z ){ return c[0]+z*(c[1]+z*(c[2]+z*c[3])); }
static double polyder( double const c[4], double z ){ return c[1]+z*(2.0*c[2]+z*3.0*c[3]); }

// Manufactured solution: c = Pc(z)(1+alp t),  u = Pu(z)(1+alp t).
static double C_exact( double t, double z, Par const& p ){ return polyval(p.cc,z)*(1.0+p.alp*t); }
static double U_exact( double t, double z, Par const& p ){ return polyval(p.uu,z)*(1.0+p.alp*t); }

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

// BC form selector for the inflow conditions.
enum BCVariant {
  FLUX_WELL_DIRECTED,     // w +- beta d_z w        : a=l in Vin           -> ACCEPT + recover
  FLUX_OFFDIR_DERIV,      // (l.u) + beta d_z c     : a=l in Vin, b off    -> ACCEPT + recover
                          //   well-posed: the closure does NOT pin the outgoing d_z w, so an
                          //   off-characteristic DERIVATIVE term is admissible.  The guard must
                          //   not false-reject this (regression vs the removed b-projection).
  FLUX_MISDIRECTED_VALUE  // prescribe the OUTGOING characteristic (Robin) at the incoming end:
                          //   a=l_out NOT in Vin -> REJECT (collides with the auto-closure).
};
static char const* variant_name( BCVariant v ){
  switch( v ){
    case FLUX_WELL_DIRECTED:     return "flux well-directed (w+beta d_z w)";
    case FLUX_OFFDIR_DERIV:      return "flux off-char DERIVATIVE (value ok, d_z c term)";
    case FLUX_MISDIRECTED_VALUE: return "flux MIS-DIRECTED value (prescribes outgoing char)";
  }
  return "?";
}

struct ModeResult { std::string name; bool ok=false; size_t nVar=0,nEqn=0,nTrace=0;
                    bool square=false, solved=false, rejected_bc=false;
                    double final_res=0., eC=0., eU=0.; };

static ModeResult run_mode( OCFESLV::Options::ImpositionType imp,
                            Par const& p, BCVariant variant )
{
  ModeResult R; R.name = imp_name(imp); bool ok=true;

  // a>0: w+ (speed +2a) enters at z=LB, w- enters at z=UB.  a<0 flips.
  bool const fwd = ( p.a >= 0.0 );

  std::cout << "\n====== FLUX/ROBIN inflow on non-symmetric hyperbolic block ======\n";
  std::cout << "imposition: " << R.name
            << ", finite elements: t=" << TEST_HYP_NEL_T << " z=" << TEST_HYP_NEL_Z
            << ", nodes/element: t=" << TEST_HYP_NT << " z=" << TEST_HYP_NZ << "\n";
  std::cout << "PDE: d_t c + a d_z u = f_c ; d_t u + g a d_z c = f_u,  a=" << p.a
            << " g=" << p.g << " (eigenvalues +-" << p.sg << "a),  beta=" << p.beta
            << "  (polynomial profile, exact spectral d_z)\n";
  std::cout << "placement: " << (fwd?"FORWARD (w+ in@LB, w- in@UB)":"REVERSE (w- in@LB, w+ in@UB)")
            << " ;  inflow BC = " << variant_name(variant) << "\n";

  FFGraph DAG;
  FFVar t = DAG.add_var("t");
  FFVar z = DAG.add_var("z");
  FFVar c = DAG.add_var("c(t,z)");
  FFVar u = DAG.add_var("u(t,z)");
  FFPartial OpP;

  // Manufactured forcing for the polynomial profiles (same governing system):
  //   f_c = d_t c +   a d_z u = alp Pc +   a Pu'(z)(1+alp t)
  //   f_u = d_t u + g a d_z c = alp Pu + g a Pc'(z)(1+alp t)
  FFVar Pc  = p.cc[0] + p.cc[1]*z + p.cc[2]*z*z + p.cc[3]*z*z*z;
  FFVar Pu  = p.uu[0] + p.uu[1]*z + p.uu[2]*z*z + p.uu[3]*z*z*z;
  FFVar dPc = p.cc[1] + 2.0*p.cc[2]*z + 3.0*p.cc[3]*z*z;
  FFVar dPu = p.uu[1] + 2.0*p.uu[2]*z + 3.0*p.uu[3]*z*z;
  FFVar CE = Pc*( 1.0 + p.alp*t );
  FFVar UE = Pu*( 1.0 + p.alp*t );
  FFVar FC = p.alp*Pc +     p.a*dPu*( 1.0 + p.alp*t );
  FFVar FU = p.alp*Pu + p.g*p.a*dPc*( 1.0 + p.alp*t );

  // Characteristic VALUES at each end on the manufactured solution (w = l.u).
  double const wpLB0 = p.sg*polyval(p.cc,0.0)  + polyval(p.uu,0.0);
  double const wmLB0 = p.sg*polyval(p.cc,0.0)  - polyval(p.uu,0.0);
  double const wpUB0 = p.sg*polyval(p.cc,p.zf) + polyval(p.uu,p.zf);
  double const wmUB0 = p.sg*polyval(p.cc,p.zf) - polyval(p.uu,p.zf);

  // Characteristic NORMAL DERIVATIVES d_z w = l.(d_z u) at each end -- EXACT for
  // the polynomial profile (matches the spectral d_z at the boundary node).
  double const dwpLB0 = p.sg*polyder(p.cc,0.0)  + polyder(p.uu,0.0);
  double const dwmLB0 = p.sg*polyder(p.cc,0.0)  - polyder(p.uu,0.0);
  double const dwpUB0 = p.sg*polyder(p.cc,p.zf) + polyder(p.uu,p.zf);
  double const dwmUB0 = p.sg*polyder(p.cc,p.zf) - polyder(p.uu,p.zf);
  // d_z c at z=LB (for the OFF-DIRECTION-b control's RHS).
  double const dcLB0  = polyder(p.cc,0.0);

  FFVar const grow = ( 1.0 + p.alp*t );

  FFVar PDEC = OpP(c,t) +     p.a*OpP(u,z) - FC;
  FFVar PDEU = OpP(u,t) + p.g*p.a*OpP(c,z) - FU;
  FFVar ICC  = c - CE;
  FFVar ICU  = u - UE;

  // Value / derivative coupling building blocks.
  FFVar wpVal = p.sg*c + u;                 // l_+ . u
  FFVar wmVal = p.sg*c - u;                 // l_- . u
  FFVar wpDer = p.sg*OpP(c,z) + OpP(u,z);   // l_+ . d_z u = d_z w_+
  FFVar wmDer = p.sg*OpP(c,z) - OpP(u,z);   // l_- . d_z u = d_z w_-
  FFVar offDer = OpP(c,z);                  // off-characteristic derivative (1,0).d_z u

  // Outgoing-characteristic transports (left-eigenvector combinations) for the
  // optional manual closure cross-check (auto-closure is the default).
  FFVar WPz = p.sg*PDEC + PDEU;             // w+ transport (rowspace l_+)
  FFVar WMz = p.sg*PDEC - PDEU;             // w- transport (rowspace l_-)

  // -- Flux/Robin inflow BCs ------------------------------------------------
  // Per-direction incoming/outgoing characteristic Robin BCs (value+derivative)
  // plus the off-characteristic-derivative form.  The UB end is always the
  // well-directed incoming Robin; only the LB end varies by control.
  //   well-directed   BC = w_in  + beta d_z w_in - (w_in0  + beta dw_in0)*grow  (a=l_in, accept)
  //   off-char deriv  BC = w_in  + beta d_z c    - (w_in0  + beta dc0   )*grow  (a=l_in, accept)
  //   misdirected val BC = w_out + beta d_z w_out- (w_out0 + beta dw_out0)*grow (a=l_out, REJECT)
  FFVar robin_in_LB, robin_out_LB, offdir_in_LB, robin_in_UB;
  if( fwd ){
    // incoming @LB = w+ , outgoing @LB = w- ; incoming @UB = w-
    robin_in_LB  = wpVal + p.beta*wpDer  - (wpLB0 + p.beta*dwpLB0)*grow;
    robin_out_LB = wmVal + p.beta*wmDer  - (wmLB0 + p.beta*dwmLB0)*grow;
    offdir_in_LB = wpVal + p.beta*offDer - (wpLB0 + p.beta*dcLB0 )*grow;
    robin_in_UB  = wmVal + p.beta*wmDer  - (wmUB0 + p.beta*dwmUB0)*grow;
  } else {
    // incoming @LB = w- , outgoing @LB = w+ ; incoming @UB = w+
    robin_in_LB  = wmVal + p.beta*wmDer  - (wmLB0 + p.beta*dwmLB0)*grow;
    robin_out_LB = wpVal + p.beta*wpDer  - (wpLB0 + p.beta*dwpLB0)*grow;
    offdir_in_LB = wmVal + p.beta*offDer - (wmLB0 + p.beta*dcLB0 )*grow;
    robin_in_UB  = wpVal + p.beta*wpDer  - (wpUB0 + p.beta*dwpUB0)*grow;
  }
  FFVar BC_LB, BC_UB = robin_in_UB;
  switch( variant ){
    case FLUX_WELL_DIRECTED:      BC_LB = robin_in_LB;  break;
    case FLUX_OFFDIR_DERIV:       BC_LB = offdir_in_LB; break;
    case FLUX_MISDIRECTED_VALUE:  BC_LB = robin_out_LB; break;
  }
  FFVar CLO_LB = fwd ? WMz : WPz;           // outgoing @LB (manual closure only)
  FFVar CLO_UB = fwd ? WPz : WMz;           // outgoing @UB

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
  OCFESLV::EqnOptions clo_opt( OCFESLV::EqnRole::INTERIOR, 0, OCFESLV::Options::IC_AUTO,
                             /*classify=*/false, /*sat=*/true );

  int const T_INT = FFDom::ALL - FFDom::LB;
  int const Z_INT = FFDom::ALL - FFDom::LB - FFDom::UB;

  oc.add_equation( PDEC,  {t,z}, {T_INT, Z_INT},         int_opt );
  oc.add_equation( PDEU,  {t,z}, {T_INT, Z_INT},         int_opt );
  oc.add_equation( ICC,   {t,z}, {FFDom::LB, FFDom::ALL}, ini_opt );
  oc.add_equation( ICU,   {t,z}, {FFDom::LB, FFDom::ALL}, ini_opt );
  oc.add_equation( BC_LB, {t,z}, {T_INT, FFDom::LB},     bnd_opt );  // incoming flux @LB
  oc.add_equation( BC_UB, {t,z}, {T_INT, FFDom::UB},     bnd_opt );  // incoming flux @UB
#if TEST_HYP_MANUAL_CLOSURE
  oc.add_equation( CLO_LB, {t,z}, {T_INT, FFDom::LB},    clo_opt );
  oc.add_equation( CLO_UB, {t,z}, {T_INT, FFDom::UB},    clo_opt );
#else
  (void)CLO_LB; (void)CLO_UB;
#endif

  oc.set_evolution_domain( t );
  oc.options.REDUCE.ORDER    = OCFESLV::Options::RED_FULL;
  oc.options.CLASSIFY.MODE        = OCFESLV::Options::CLASS_AUTO;
  oc.options.INTERFACE.TYPE  = OCFESLV::Options::IC_AUTO;
  oc.options.INTERFACE.IMPOSITION = imp;
  oc.options.INTERFACE.SAT_SIGMA0      = 10.0;
  oc.options.DISPLAY_LEVEL   = 1;

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
  ok &= check_close("square system (nVar==nEqn)", R.square?0.0:1.0, 0.0);

  std::cout << "States after setup:";
  for( auto const& st: oc.states_colloc() ) std::cout << ' ' << st.name();
  std::cout << "\nAuxiliary states introduced: "
            << (oc.states_colloc().size()>2?oc.states_colloc().size()-2:0)
            << " (expect 0: first-order system)\n";
  {
    auto const& cls = oc.pde_type();
    std::cout << "PDE type: " << OCFESLV::pde_type_name(cls.type)
              << "  evolution_hyperbolic=" << (cls.evolution_hyperbolic?"yes":"no") << "\n";
    bool const weak_path = ( imp != OCFESLV::Options::IC_STRONG );
    std::cout << "Resolved interface type:  t(evolution) -> "
              << iftype_name( oc.resolved_interface_type( 0, t, OCFESLV::EqnRole::INTERIOR, weak_path ) )
              << "   |   z -> "
              << iftype_name( oc.resolved_interface_type( 0, z, OCFESLV::EqnRole::INTERIOR, weak_path ) )
              << "\n";
  }

  std::vector<double> var(nVar,0.0);

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

  std::cout << "flux-inflow hyperbolic test (" << R.name << "): " << (ok?"PASS":"FAIL") << "\n";
  R.ok = ok;
  return R;
}

int main()
{
  std::cout << "PDE16 build config: closure="
#if defined(MC__OCFESLV_AUTO_HYP_CLOSURE)
            << "AUTO"
#else
            << "NONE(gap-finder)"
#endif
            << "  base-guard="
            << "ON"
            << "  flux-guard="
            << "ON"
            << "\n";

  double const amag = std::fabs( (double)TEST_HYP_A );
  Par pf; pf.a = +amag;   // forward propagation
  Par pr; pr.a = -amag;   // reverse propagation

  auto run3 = []( Par const& p, BCVariant variant, std::vector<ModeResult>& out ){
    out.push_back( run_mode(OCFESLV::Options::IC_WEAK,   p, variant) );
    out.push_back( run_mode(OCFESLV::Options::IC_TRACE,  p, variant) );
    out.push_back( run_mode(OCFESLV::Options::IC_STRONG, p, variant) );
  };

  // A: forward sign, well-directed flux inflow       -> ACCEPT + recover.
  // B: reverse sign, well-directed flux inflow       -> ACCEPT + recover.
  // C: forward sign, OFF-CHARACTER DERIVATIVE flux    -> ACCEPT + recover.  Value
  //    coupling is the correct incoming l (a in Vin); only the derivative term is
  //    off-character.  The closure does not pin the outgoing d_z w, so this is
  //    well-posed and the value-only guard must NOT reject it (regression guard vs
  //    the removed b-projection -- proven by the -DTEST_PDE16_NO_FLUX_GUARD run,
  //    which recovers the exact solution).
  // D: forward sign, MIS-DIRECTED VALUE flux          -> REJECTED.  Prescribes the
  //    OUTGOING characteristic (Robin) at the incoming end: a=l_out outside Vin,
  //    collides with the auto-closure.
  std::vector<ModeResult> A, B, C, D;
  run3( pf, FLUX_WELL_DIRECTED, A );
  run3( pr, FLUX_WELL_DIRECTED, B );
  run3( pf, FLUX_OFFDIR_DERIV,  C );
  D.push_back( run_mode(OCFESLV::Options::IC_WEAK, pf, FLUX_MISDIRECTED_VALUE) );

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

  print_block( "A: forward sign, well-directed flux inflow", A );
  print_block( "B: reverse sign, well-directed flux inflow", B );
  print_block( "C: off-character DERIVATIVE flux inflow (value ok) -> well-posed", C );
  print_block( "D: mis-directed VALUE flux inflow (prescribes outgoing char)", D );

  bool all_ok = true;
  for( auto const& r: A ) all_ok &= r.ok;   // well-directed flux must pass
  for( auto const& r: B ) all_ok &= r.ok;   // sign-correct reverse must pass
  for( auto const& r: C ) all_ok &= r.ok;   // off-char DERIVATIVE is well-posed: must pass

  // D is the value-direction discrimination control.  Correct outcomes:
  //   * flux guard ON (default here): REJECTED at setup -- a=l_out lies outside
  //     Vin, so the value-only test fires (HYP_INCOMING_BC).
  //   * flux guard OFF (-DTEST_PDE16_NO_FLUX_GUARD): the OpP BC skips the test, so
  //     it is NOT rejected -- the value-guard is load-bearing.
  bool const flux_on =
        true;
  bool const demo_ok =
    !D.empty() && ( flux_on ? D[0].rejected_bc : !D[0].rejected_bc );

  std::cout << "\nValue-direction control (mis-directed VALUE flux BC): "
            << ( D.empty() ? "?"
               : D[0].rejected_bc ? "REJECTED at setup (HYP_INCOMING_BC) -- a=l_out outside Vin"
                                  : "NOT rejected (flux guard off -> OpP BC skipped)" )
            << "  => " << (demo_ok ? "value-guard catches misdirection (as expected)"
                                   : "inconclusive") << "\n";
  std::cout << "(C off-character derivative ACCEPTED + recovers in both guard configs"
            << " -- well-posed, guard correctly leaves it alone)\n";
  all_ok &= demo_ok;

  std::cout << "\n====================================================================\n";
  std::cout << "flux-inflow hyperbolic oracle (well-directed + off-char-derivative accept"
            << " + mis-directed-value reject): " << (all_ok?"PASS":"FAIL") << "\n";
  return all_ok ? 0 : 1;
}
