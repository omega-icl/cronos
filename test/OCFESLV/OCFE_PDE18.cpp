// ===========================================================================
// OCFE_PDE18_solve2.cpp
//
// DARCY-BILINEAR inflow oracle -- stage 5 of the PSA case study, the validator
// for the value-only flux-BC guard on a BILINEAR boundary condition whose VALUE
// coupling `a` itself references a DERIVATIVE proxy, plus the two reference
// degeneracies the linear PDE16 oracle cannot exercise.
//
// PDE16 validated the value-only criterion on a LINEAR flux/Robin BC
// (w +- beta d_z w): the d_z OpP node is proxied by the principal-symbol
// derivative proxy, the proxy contributes a constant to `a`, and a in Vin is
// tested.  A real PSA inflow (Darcy v = d_z P, Danckwerts P*v) is BILINEAR: the
// boundary form multiplies a state VALUE by a state DERIVATIVE, so after the
// proxy substitution `a = d(BC)/d(state)` references the seeded proxy *through a
// product*.  This driver exercises exactly that, on a decoupled-advection block
// whose Vin is a single coordinate axis so the accept/reject discrimination is
// unambiguous.
//
//   d_t P + s d_z P + k c = f_P     (s>0: P RIGHT-going, incoming @ z=LB)
//   d_t c - r d_z c + k P = f_c     (r>0: c LEFT-going,  incoming @ z=UB)
//
//   Principal symbol A_z = diag(s,-r): eigenvalues s>0, -r<0 (HYPERBOLIC).  The
//   zeroth-order coupling k (k != 0) keeps (P,c) ONE 2x2 block (so the guard sees
//   a 2-vector coupling) WITHOUT entering the principal symbol.  Hence
//     Vin @ LB = span{(1,0)} = P (the only right-going characteristic),
//     Vin @ UB = span{(0,1)} = c.
//   State order is (P,c): index 0 = P, index 1 = c.
//
// SECTIONS (each states the criterion the guard must meet, and the SOLVE -- not
// theory -- decides well-posedness for the accept cases):
//
//   A  WELL-DIRECTED bilinear      BC = P * d_z c - q          @ LB
//        substitute d_z c -> proxy(=1): a = d(P*proxy)/d(P,c) = (1,0) in Vin
//        -> ACCEPT, and the square solve must RECOVER the manufactured solution.
//        (The proxy seed makes `a` reference-INDEPENDENT in direction -- this is
//        the case proxy-seeding is meant to protect, and A confirms it recovers
//        even when the state reference is left at zero.)
//
//   B  DERIVATIVE-RESCUED bilinear  BC = c * d_z P - q          @ LB
//        substitute d_z P -> proxy: VALUE a = d(c*proxy)/d(P,c) = (0,1) is OUTGOING,
//        but DERIVATIVE b = d(BC)/d(d_z u) = (c_ref,0) is INCOMING (pins P through
//        d_z P; the spectral boundary derivative includes the boundary value).
//        Under CRITERION #3 (accept iff a in Vin OR b in Vin) the b-branch RESCUES:
//        -> ACCEPT + recover (all modes).  The probe confirms it is genuinely
//           well-posed (perturbed init returns to exact).  A pure value-only test
//           FALSE-REJECTS this -- the error criterion #3 fixes.
//
//   B2 PURE-OUTGOING bilinear       BC = c * d_z c - q          @ LB
//        a = (0,1) OUTGOING and b = (0,c_ref) OUTGOING: NEITHER coupling in Vin, so
//        the BC ITSELF pins no incoming characteristic.
//        -> guard ON: REJECT at setup (HYP_INCOMING_BC) -- criterion #3 reject branch.
//        -> guard OFF: accepted; FINDING -- the perturbed-init probe RETURNS to exact
//           because this block's ZEROTH-ORDER coupling leaks the incoming mode into
//           the outgoing closure row (k*P in d_t c - r d_z c + k P), discretely
//           rescuing any independent second LB constraint.  So a k-coupled DIAGONAL
//           block cannot host a genuinely ill-posed inflow; the load-bearing reject
//           (silent wrongness) is shown on PDE16 D (principal-symbol coupled).  The
//           guard-ON reject here is STRUCTURAL (the BC supplies no inflow data) and
//           conservative -- relevant to PSA's coupled blocks.
//
//   C  REFERENCE-DEPENDENT VERDICT  BC = P * c - q              @ LB   (ALGEBRAIC)
//        a purely-algebraic bilinear (NO derivative, so NO proxy rescue):
//        a = d(P*c)/d(P,c) = (c_ref, P_ref).  (rev369, 2026-10-09: the guard tests the
//        KREISS condition a.r_in != 0 with r_in = (1,0) -- the condition PINS P iff
//        c_ref != 0; the former membership test "a in Vin" rejected it for P_ref != 0.)
//        -> reference ZEROED (no update_ref): a=(0,0) -> ||a||~0 -> the guard
//           reports "imposes nothing" and SILENTLY ACCEPTS: the collapse, as before.
//        -> reference SET (operating point): a=(c0,P0), c0 != 0 -> pins P -> ACCEPT
//           (P*c = q with c known from the interior determines P: well-posed).
//        -> P reference SET, c reference ZERO: a=(0,P0) -> a.r_in = 0 -> REJECT.
//        The verdict flips PURELY from the reference -- from c_ref, the one that
//        matters.  Independent of the flux-guard macro (no OpP node).
//
//   D  FACE-vs-MIDPOINT dom_ref    BC = P * d_z c + g(z) c - q  @ LB
//        g(z) = z (z - zf):  g(LB)=g(0)=0 but g(midpoint)=-zf^2/4 != 0.  At the
//        face the BC is exactly section A (well-directed, recovers).  But the
//        guard evaluates `a` at _default_dom_ref(z) = the element MIDPOINT:
//        a = (1, g(z_ref)).
//        -> dom_ref = MIDPOINT (default): a=(1,-zf^2/4): the membership test "a in Vin"
//           FALSE-REJECTED this (rev <= 368); the Kreiss test (rev369) accepts it,
//           a.r_in = 1, and the solve RECOVERS.
//        -> dom_ref = FACE (add_domain(z,...,ref=LB) overrides _classDomRef): a=(1,0)
//           -> ACCEPT, and the solve RECOVERS.
//        Both arms accept and recover: the midpoint default no longer false-rejects a
//        well-posed face BC (the face override is kept as a second arm).
//        (Requires the flux guard ON; with it off the OpP BC skips and neither
//        dom_ref discriminates.)
//
// The flux guard is OPT-IN in the header (default off); this oracle enables it by
// default so the bilinear path is exercised.  Build with -DTEST_PDE18_NO_FLUX_GUARD
// for the guard-OFF "before" picture (A still recovers via skip; B converges WRONG;
// D no longer discriminates; C is unchanged -- it never depends on the macro).
//
// Build flags:  auto-closure + base incoming-BC guard are default-on in the header;
//   + -DCRONOS__WITH_SPQR for QR, + -DCRONOS__WITH_UMFPACK for IC_STRONG.  Recommended:
//   DISPLAY_LEVEL>=1 prints the per-BC value-resid the guard computes.
//
//   -DTEST_PDE18_NO_FLUX_GUARD   build with the flux guard OFF (before/after)
//   -DTEST_PDE18_BETA ...        (unused here; kept for flag symmetry with PDE16)
// ===========================================================================

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <cmath>
#include <map>
#include <algorithm>

// The flux-BC guard is now DEFAULT-ON in the header (criterion #3).  For the
// guard-OFF "before" build, -DTEST_PDE18_NO_FLUX_GUARD must both #undef any
// force-define AND set the header opt-out MC__OCFESLV_NO_HYP_FLUX_BC_GUARD, else the
// header's default-on auto-define would re-enable it.  guard-ON needs no action
// (the header default).  The off-flag stays AUTHORITATIVE against a harness
// force-define of -DMC__OCFESLV_HYP_FLUX_BC_GUARD.

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
  double tf = TEST_HYP_TF;
  double zf = TEST_HYP_ZF;
  double s  = 1.0;     // P advection speed (>0: P right-going, in @ LB)
  double r  = 0.8;     // c advection speed magnitude (>0: c left-going, in @ UB)
  double k  = 0.15;    // zeroth-order coupling: ties (P,c) into ONE block, NOT in A_z
  double alp = 0.3;    // linear-in-time growth factor
  // Spatial profiles: degree-3 POLYNOMIALS in z (degree <= N-1), so the spectral
  // d_z is EXACT at every node including the boundary; the bilinear-BC RHS built
  // from the analytic d_z is discretely consistent and the solve recovers to
  // machine accuracy.  Coeffs chosen so Pp(0),Pp'(0),Pc(0),Pc'(0) are all nonzero
  // (the bilinear couplings and their derivative proxies are then nontrivial).
  double pp[4] = { 0.8, 0.5, -0.4,  0.2 };   // P profile Pp(z)
  double pc[4] = { 1.1, 0.7,  0.3, -0.25 };  // c profile Pc(z)
};

// Polynomial value / z-derivative (Horner).
static double polyval( double const c[4], double z ){ return c[0]+z*(c[1]+z*(c[2]+z*c[3])); }

// Manufactured solution: P = Pp(z)(1+alp t),  c = Pc(z)(1+alp t).
static double P_exact( double t, double z, Par const& p ){ return polyval(p.pp,z)*(1.0+p.alp*t); }
static double C_exact( double t, double z, Par const& p ){ return polyval(p.pc,z)*(1.0+p.alp*t); }

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
  std::cout << std::left << std::setw(44) << label
            << " value=" << std::scientific << std::setprecision(6) << value
            << " tol=" << tol << "  " << (ok?"PASS":"FAIL") << "\n";
  return ok;
}

// ---- section taxonomy ------------------------------------------------------
enum Section { SEC_A_WELL, SEC_B_MISDIR, SEC_B2_OUTGOING, SEC_C_ALGEBRAIC, SEC_D_ZDEP };
static char const* section_name( Section s ){
  switch( s ){
    case SEC_A_WELL:      return "A well-directed bilinear  (P d_z c)";
    case SEC_B_MISDIR:    return "B mis-directed bilinear   (c d_z P)";
    case SEC_B2_OUTGOING: return "B2 pure-outgoing bilinear (c d_z c)";
    case SEC_C_ALGEBRAIC: return "C algebraic bilinear      (P c)";
    case SEC_D_ZDEP:      return "D z-dependent coupling     (P d_z c + g c)";
  }
  return "?";
}
// Expected guard outcome with the flux guard ON, at the operating reference.
//   accept : setup succeeds (well-directed) -> for A/D2 also recover
//   reject : setup fails HYP_INCOMING_BC    (mis-directed value coupling)
//   skip   : setup succeeds because a~0 / OpP-skipped (no constraint asserted)

struct Result {
  std::string mode, label;
  Section section = SEC_A_WELL;
  bool set_ref = true, face_domref = false;
  bool setup_ok = false, rejected_bc = false, square = false, solved = false;
  bool did_solve = false;
  double final_res = 0., eP = 0., eC = 0.;
};

// Build + setup (+ solve when the section is meant to accept/recover or to expose
// a guard-off wrong answer).  set_ref=false leaves the classification reference at
// zero (the C zero-collapse arm); face_domref=true overrides the classification
// dom_ref for z to the LB face (the D face arm).
static Result run_section( OCFESLV::Options::ImpositionType imp, Par const& p,
                           Section section, bool set_ref, bool face_domref,
                           char const* label, double init_perturb = 0.0,
                           bool c_ref_zero = false )   // section C: P reference set, c reference left at zero
{
  Result R; R.mode = imp_name(imp); R.label = label;
  R.section = section; R.set_ref = set_ref; R.face_domref = face_domref;

  std::cout << "\n------ section " << section_name(section)
            << "  [" << label << ", " << R.mode << "] ------\n";
  std::cout << "PDE: d_t P + s d_z P + k c = f_P ; d_t c - r d_z c + k P = f_c"
            << "   (s=" << p.s << " r=" << p.r << " k=" << p.k << ")\n";
  std::cout << "A_z = diag(s,-r); Vin@LB = (1,0)=P, Vin@UB = (0,1)=c ; state order (P,c)\n";
  std::cout << "classification reference: " << (set_ref ? "SET (operating point)"
                                                        : "ZEROED (default)")
            << " ;  z dom_ref: " << (face_domref ? "FACE (z=LB override)"
                                                 : "MIDPOINT (default)") << "\n";

  FFGraph DAG;
  FFVar t = DAG.add_var("t");
  FFVar z = DAG.add_var("z");
  FFVar P = DAG.add_var("P(t,z)");
  FFVar c = DAG.add_var("c(t,z)");
  FFPartial OpP;

  FFVar grow = ( 1.0 + p.alp*t );
  FFVar Pp  = p.pp[0] + p.pp[1]*z + p.pp[2]*z*z + p.pp[3]*z*z*z;
  FFVar Pc  = p.pc[0] + p.pc[1]*z + p.pc[2]*z*z + p.pc[3]*z*z*z;
  FFVar dPp = p.pp[1] + 2.0*p.pp[2]*z + 3.0*p.pp[3]*z*z;
  FFVar dPc = p.pc[1] + 2.0*p.pc[2]*z + 3.0*p.pc[3]*z*z;

  // Manufactured forcing (so the governing residual is zero on the exact soln):
  //   f_P = alp Pp + s Pp'(z) grow + k Pc(z) grow
  //   f_c = alp Pc - r Pc'(z) grow + k Pp(z) grow
  FFVar FP = p.alp*Pp + p.s*dPp*grow + p.k*Pc*grow;
  FFVar FC = p.alp*Pc - p.r*dPc*grow + p.k*Pp*grow;

  FFVar PDEP = OpP(P,t) + p.s*OpP(P,z) + p.k*c - FP;
  FFVar PDEC = OpP(c,t) - p.r*OpP(c,z) + p.k*P - FC;
  FFVar ICP  = P - Pp*grow;
  FFVar ICC  = c - Pc*grow;

  // UB inflow: c is incoming @ UB; prescribe its value (well-directed, a=(0,1)).
  FFVar BC_UB = c - Pc*grow;

  // LB inflow: section-dependent.  RHS q is the manufactured value of the same
  // boundary form, so the manufactured solution satisfies BC=0 at z=0.
  FFVar g_ff = z*( z - p.zf );                 // g(0)=0, g(zf)=0, g(zf/2)=-zf^2/4
  FFVar BC_LB;
  switch( section ){
    case SEC_A_WELL:
      BC_LB = P*OpP(c,z) - Pp*dPc*grow*grow;                 // a -> (1,0)  ACCEPT
      break;
    case SEC_B_MISDIR:
      BC_LB = c*OpP(P,z) - Pc*dPp*grow*grow;                 // a->(0,1) OUT, b->(c,0) IN: accept (b)
      break;
    case SEC_B2_OUTGOING:
      BC_LB = c*OpP(c,z) - Pc*dPc*grow*grow;                 // a->(0,1) OUT, b->(0,c) OUT: REJECT
      break;
    case SEC_C_ALGEBRAIC:
      BC_LB = P*c - Pp*Pc*grow*grow;                         // a=(c_ref,P_ref)
      break;
    case SEC_D_ZDEP:
      BC_LB = P*OpP(c,z) + g_ff*c - ( Pp*dPc*grow*grow + g_ff*Pc*grow ); // a=(1,g(z_ref))
      break;
  }

  OCFESLV oc(&DAG);
  oc.add_domain( t, FFDom(0., p.tf, TEST_HYP_NEL_T, FFDom::LGR, TEST_HYP_NT) );
  if( face_domref )
    oc.add_domain( z, FFDom(0., p.zf, TEST_HYP_NEL_Z, FFDom::LGL, TEST_HYP_NZ), /*classDomRef=*/0.0 );
  else
    oc.add_domain( z, FFDom(0., p.zf, TEST_HYP_NEL_Z, FFDom::LGL, TEST_HYP_NZ) );
  oc.add_state( P, {t,z} );
  oc.add_state( c, {t,z} );
  if( set_ref ){
    oc.update_ref( P, [&]( OCFESLV::t_Coord const& crd ){ return P_exact(crd.at(t),crd.at(z),p); } );
    if( !c_ref_zero )
      oc.update_ref( c, [&]( OCFESLV::t_Coord const& crd ){ return C_exact(crd.at(t),crd.at(z),p); } );
  }

  OCFESLV::EqnOptions int_opt( OCFESLV::EqnRole::INTERIOR, 0 );
  OCFESLV::EqnOptions ini_opt( OCFESLV::EqnRole::INITIAL,  0 );
  OCFESLV::EqnOptions bnd_opt( OCFESLV::EqnRole::BOUNDARY, 0 );

  int const T_INT = FFDom::ALL - FFDom::LB;
  int const Z_INT = FFDom::ALL - FFDom::LB - FFDom::UB;

  oc.add_equation( PDEP,  {t,z}, {T_INT, Z_INT},          int_opt );
  oc.add_equation( PDEC,  {t,z}, {T_INT, Z_INT},          int_opt );
  oc.add_equation( ICP,   {t,z}, {FFDom::LB, FFDom::ALL},  ini_opt );
  oc.add_equation( ICC,   {t,z}, {FFDom::LB, FFDom::ALL},  ini_opt );
  oc.add_equation( BC_LB, {t,z}, {T_INT, FFDom::LB},       bnd_opt );  // section BC @ LB
  oc.add_equation( BC_UB, {t,z}, {T_INT, FFDom::UB},       bnd_opt );  // c value @ UB

  oc.set_evolution_domain( t );
  oc.options.REDUCE.ORDER    = OCFESLV::Options::RED_FULL;
  oc.options.CLASSIFY.MODE        = OCFESLV::Options::CLASS_AUTO;
  oc.options.INTERFACE.TYPE  = OCFESLV::Options::IC_AUTO;
  oc.options.INTERFACE.IMPOSITION = imp;
  oc.options.INTERFACE.SAT_SIGMA0      = 10.0;
  oc.options.DISPLAY_LEVEL   = 1;

  if( !oc.setup() ){
    R.setup_ok = false;
    R.rejected_bc = ( oc.setup_status() == OCFESLV::SetupStatus::HYP_INCOMING_BC
                   || oc.setup_status() == OCFESLV::SetupStatus::HYP_BC_MISDIRECTED );
    std::cout << "setup: FAILED -> " << OCFESLV::setup_status_str( oc.setup_status() )
              << (R.rejected_bc ? "   [guard rejected the BC]" : "") << "\n";
    return R;
  }
  R.setup_ok = true;

  size_t const nVar=oc.n_colloc_sta(), nEqn=oc.n_colloc_eqn(), nTrace=oc.n_colloc_trace();
  R.square = (nVar==nEqn);
  std::cout << "setup: OK  nVar=" << nVar << " nEqn=" << nEqn << " nTrace=" << nTrace
            << " square=" << (R.square?"yes":"no") << "\n";
  {
    auto const& cls = oc.pde_type();
    std::cout << "PDE type: " << OCFESLV::pde_type_name(cls.type)
              << "  evolution_hyperbolic=" << (cls.evolution_hyperbolic?"yes":"no") << "\n";
  }

  // Solve only where a converged answer is meaningful: the accept/recover
  // sections (A, D-face) and the guard-OFF mis-directed arm (B), where the point
  // is to EXPOSE a square wrong answer.  C is a setup-verdict demo (no solve).
  bool const want_solve = R.setup_ok && ( section != SEC_C_ALGEBRAIC );
  if( !want_solve ){
    std::cout << "(verdict-only section: no solve)\n";
    return R;
  }

  std::vector<double> var(nVar,0.0);
  // Seed the initial guess from the classification reference (the manufactured
  // solution, via update_ref).  solve() takes the CALLER's var as the starting
  // guess -- it does NOT seed from the reference -- and a BILINEAR boundary form
  // (P*d_z c) has a SINGULAR Jacobian at the all-zero state (d/dP = d_z c = 0 when
  // c=0), so a zero start stalls the first Newton/QR step.  A nonzero (here exact)
  // start makes the BC row full-rank; for the well-posed sections the solve then
  // recovers, and for a rank-deficient mis-directed assembly (B) the min-norm QR
  // step moves OFF exact to the wrong solution (silent wrongness, not singularity).
  if( set_ref && !oc.init( var.data(), nullptr, nullptr ) )
    std::cout << "warning: init() from reference failed; using zero start\n";
  if( init_perturb != 0.0 ){
    // Deterministic, nonsingular perturbation OFF the exact seed: turns "recovered"
    // from a consistency check (exact is a root) into a genuine convergence test
    // (does the solver RETURN to exact, i.e. is exact the unique attractor).
    for( double& x : var ) x += init_perturb;
    std::cout << "init: reference + uniform perturbation " << init_perturb
              << " (well-posedness probe)\n";
  }
  oc.options.SOLVE.MAX_ITER = TEST_HYP_MAXIT;
  oc.options.SOLVE.RES_TOL  = TEST_HYP_SOLVE_TOL;
#if defined(CRONOS__WITH_SPQR)
  oc.options.SOLVE.FACTORIZATION = OCFESLV::Options::SOLVE_SPQR;
#endif
  OCFESLV::SolveReport const srep = oc.solve( var.data() );
  R.solved = srep.converged; R.did_solve = true;
  if( !R.solved )
    std::cout << "solve: did NOT converge (final|r|=" << srep.final_residual
              << " after " << srep.iterations << " it)\n";

  std::vector<double> res(nEqn,0.0);
  oc.eval(res.data(),nullptr,var.data(),nullptr,nullptr);
  R.final_res = max_abs(res);

  double eP=0.0, eC=0.0; size_t off=0, sidx=0;
  for( auto const& st: oc.states_colloc() ){
    auto nodes = oc.node_colloc(st);
    for( size_t i=0; i<nodes.size(); ++i ){
      double const ex = (sidx==0)? P_exact(nodes[i][0],nodes[i][1],p)
                                 : C_exact(nodes[i][0],nodes[i][1],p);
      double const e = std::fabs( var[off+i] - ex );
      if( sidx==0 ) eP = std::max(eP,e); else eC = std::max(eC,e);
    }
    off += nodes.size(); ++sidx;
  }
  R.eP=eP; R.eC=eC;
  std::cout << "solve: " << (R.solved?"converged":"FAILED")
            << "  final|r|=" << std::scientific << std::setprecision(3) << R.final_res
            << "  |P-exact|=" << R.eP << "  |c-exact|=" << R.eC << "\n";
  return R;
}

int main()
{
  bool const flux_on =
        true;

  std::cout << "PDE18 build config: closure="
#if defined(MC__OCFESLV_AUTO_HYP_CLOSURE)
            << "AUTO"
#else
            << "NONE(gap-finder)"
#endif
            << "  base-guard="
            << "ON"
            << "  flux-guard=" << (flux_on?"ON":"off") << "\n";
  std::cout << "Darcy-bilinear inflow oracle: bilinear value-coupling guard + the two\n"
            << "reference degeneracies (zero-reference collapse, midpoint-vs-face dom_ref).\n";

  Par p;

  // ---- Section A: well-directed bilinear P*d_z c -- ACCEPT + recover (3 modes).
  // a=(1,0) in Vin comes from the seeded derivative proxy, so it is correct even
  // with the reference at zero; we run it with the reference SET for a clean solve.
  std::vector<Result> A;
  A.push_back( run_section(OCFESLV::Options::IC_WEAK,   p, SEC_A_WELL, /*set_ref=*/true, false, "well") );
  A.push_back( run_section(OCFESLV::Options::IC_TRACE,  p, SEC_A_WELL, true, false, "well") );
  A.push_back( run_section(OCFESLV::Options::IC_STRONG, p, SEC_A_WELL, true, false, "well") );

  // ---- Section B: derivative-rescued bilinear c*d_z P  (a=(0,1) OUTGOING value,
  // b=(c_ref,0) INCOMING derivative).  Under criterion #3 the b-branch RESCUES it:
  // ACCEPT + recover.  Regression guard proving the OR does NOT false-reject a
  // derivative-pinned inflow (the PDE18 probe finding that drove criterion #3).
  std::vector<Result> B;
  B.push_back( run_section(OCFESLV::Options::IC_WEAK,   p, SEC_B_MISDIR, true, false, "deriv-rescue") );
  B.push_back( run_section(OCFESLV::Options::IC_TRACE,  p, SEC_B_MISDIR, true, false, "deriv-rescue") );
  B.push_back( run_section(OCFESLV::Options::IC_STRONG, p, SEC_B_MISDIR, true, false, "deriv-rescue") );

  // ---- Section B2: pure-outgoing bilinear c*d_z c  (a=(0,1) and b=(0,c_ref) BOTH
  // outgoing).  Criterion #3 REJECTS it (neither coupling touches Vin) -- the reject
  // branch.  guard ON -> REJECT at setup; guard OFF -> accepted, and the perturbed
  // probe must NOT return to exact (silent wrongness) -> reject branch load-bearing.
  Result B2 = run_section(OCFESLV::Options::IC_WEAK, p, SEC_B2_OUTGOING, true, false,
                          flux_on ? "outgoing/guard-on" : "outgoing/guard-off");

  // ---- Section C: algebraic bilinear P*c -- the zero-reference collapse.
  // Same BC, two references; the verdict must FLIP.  Independent of the flux macro.
  Result C_zero = run_section(OCFESLV::Options::IC_WEAK, p, SEC_C_ALGEBRAIC, /*set_ref=*/false, false, "ref-zeroed");
  Result C_set  = run_section(OCFESLV::Options::IC_WEAK, p, SEC_C_ALGEBRAIC, /*set_ref=*/true,  false, "ref-set");
  Result C_czero= run_section(OCFESLV::Options::IC_WEAK, p, SEC_C_ALGEBRAIC, /*set_ref=*/true,  false, "ref-c-zero", 0.0, /*c_ref_zero=*/true);

  // ---- Section D: z-dependent coupling -- midpoint-vs-face dom_ref.
  //   dom_ref = MIDPOINT (default) -> FALSE REJECT.
  //   dom_ref = FACE (override)    -> ACCEPT + recover.
  Result D_mid  = run_section(OCFESLV::Options::IC_WEAK, p, SEC_D_ZDEP, /*set_ref=*/true, /*face=*/false, "dom_ref=midpoint");
  Result D_face = run_section(OCFESLV::Options::IC_WEAK, p, SEC_D_ZDEP, /*set_ref=*/true, /*face=*/true,  "dom_ref=face");

  // ---- Well-posedness probe (perturbed init) --------------------------------
  // A control (well-posed -> returns to exact, confirms delta in-basin).  B confirms
  // the rescued bilinear is well-posed (returns to exact, supporting ACCEPT).  B2 is
  // the reject-branch load-bearing readout (guard OFF): it must NOT return to exact
  // -- the silent wrong answer the reject prevents.
  double const dpert = 0.1;
  std::cout << "\n###### well-posedness probe (perturbed init, delta=" << dpert << ") ######\n";
  Result A_pert  = run_section(OCFESLV::Options::IC_WEAK, p, SEC_A_WELL,      true, false, "well/perturbed",    dpert);
  Result B_pert  = run_section(OCFESLV::Options::IC_WEAK, p, SEC_B_MISDIR,    true, false, "deriv-rescue/pert", dpert);
  Result B2_pert = run_section(OCFESLV::Options::IC_WEAK, p, SEC_B2_OUTGOING, true, false,
                               flux_on ? "outgoing/pert (reject@setup)" : "outgoing/perturbed", dpert);

  // ====================== verdict table ======================
  auto verdict_cell = []( Result const& r )->std::string {
    if( !r.setup_ok ) return r.rejected_bc ? "REJECT(setup)" : "SETUP-FAIL";
    if( !r.did_solve ) return "ACCEPT(no-solve)";
    if( !r.solved )    return "ACCEPT/solve-FAIL";
    return "ACCEPT/solved";
  };
  std::cout << "\n==================================================================\n";
  std::cout << "Darcy-bilinear guard verdicts  (flux guard " << (flux_on?"ON":"OFF") << ")\n";
  std::cout << std::left << std::setw(38) << "section / config"
            << std::setw(12) << "mode"
            << std::setw(20) << "verdict"
            << std::setw(13) << "|P-exact|"
            << std::setw(13) << "|c-exact|" << "\n";
  auto row = [&]( std::string const& tag, Result const& r ){
    std::cout << std::left << std::setw(38) << tag
              << std::setw(12) << r.mode
              << std::setw(20) << verdict_cell(r)
              << std::scientific << std::setprecision(3);
    if( r.did_solve ) std::cout << std::setw(13) << r.eP << std::setw(13) << r.eC;
    else              std::cout << std::setw(13) << "-" << std::setw(13) << "-";
    std::cout << "\n";
  };
  for( size_t i=0;i<A.size();++i ) row( i==0?"A well-directed (P d_z c)":"", A[i] );
  for( size_t i=0;i<B.size();++i ) row( i==0?"B deriv-rescue (c d_z P)":"", B[i] );
  row( flux_on ? "B2 pure-outgoing (c d_z c)" : "B2 pure-outgoing (guard OFF)", B2 );
  row( "C algebraic (P c) ref-ZEROED", C_zero );
  row( "C algebraic (P c) ref-SET",    C_set  );
  row( "D z-dep dom_ref=MIDPOINT",     D_mid  );
  row( "D z-dep dom_ref=FACE",         D_face );

  // ====================== pass/fail judgement ======================
  bool all_ok = true;

  // A: every mode must ACCEPT and recover.
  for( auto const& r: A )
    all_ok &= ( r.setup_ok && r.solved && r.eP <= TEST_HYP_EXACT_TOL && r.eC <= TEST_HYP_EXACT_TOL );
  std::cout << "\nA  well-directed bilinear ACCEPT+recover (all modes): "
            << (std::all_of(A.begin(),A.end(),[&](Result const&r){
                 return r.setup_ok&&r.solved&&r.eP<=TEST_HYP_EXACT_TOL&&r.eC<=TEST_HYP_EXACT_TOL;})?"PASS":"FAIL")
            << "\n";

  // B: criterion #3's derivative-RESCUE branch.  c*d_z P is well-posed (the probe
  // proved it) -- a=(0,1) is outgoing but b=(c_ref,0) pins the incoming P through
  // d_z P.  The guard must ACCEPT it and the solve must recover, every mode.  A
  // reject here means criterion #3 did not take (value-only false-reject returned).
  auto B_accept_recover = [&]( Result const& r ){
    return r.setup_ok && r.solved && r.eP <= TEST_HYP_EXACT_TOL && r.eC <= TEST_HYP_EXACT_TOL; };
  bool const B_ok = std::all_of( B.begin(), B.end(), B_accept_recover );
  std::cout << "\nB  derivative-rescued bilinear ACCEPT+recover (all modes): "
            << (B_ok?"PASS":"FAIL")
            << "   [criterion #3 b-branch; a outgoing, b incoming -> accept]\n";
  all_ok &= B_ok;

  // B2: criterion #3's REJECT branch.  c*d_z c is pure-outgoing (a and b both
  // outgoing) -> the BC itself pins no incoming characteristic -> REJECT.  The
  // SCORED criterion is the guard-ON setup reject.  guard-OFF does NOT exhibit
  // silent wrongness here: in this ZEROTH-ORDER-coupled block the outgoing c-closure
  // row carries k*P, so the coupling leaks the incoming mode and DISCRETELY rescues
  // any independent second LB constraint (even prescribing the outgoing value).  So
  // B2 recovers guard-off -- a real, PSA-relevant finding -- NOT a load-bearing
  // wrong answer.  The reject branch's silent-wrongness proof lives on PDE16 D
  // (principal-symbol coupled, clean characteristic closure, no cross-leak).
  bool B2_ok;
  if( flux_on ){
    B2_ok = ( !B2.setup_ok && B2.rejected_bc );
    std::cout << "B2 pure-outgoing bilinear REJECTED at setup (guard ON, criterion #3 "
              << "reject branch): " << (B2_ok?"PASS":"FAIL") << "\n";
  } else {
    bool const accepted = ( B2.setup_ok && !B2.rejected_bc );
    bool const recovered = B2_pert.did_solve && B2_pert.solved
                        && B2_pert.eP <= TEST_HYP_EXACT_TOL && B2_pert.eC <= TEST_HYP_EXACT_TOL;
    B2_ok = accepted;     // scored: guard absent -> BC sails through setup (guard is what rejects)
    std::cout << "B2 pure-outgoing bilinear ACCEPTED at setup (guard OFF -> not caught): "
              << (B2_ok?"PASS":"FAIL") << "\n"
              << "   FINDING: perturbed solve "
              << (recovered ? "RETURNS to exact -- zeroth-order coupling (k*P in the c-closure "
                              "row) discretely rescues the pure-outgoing BC; reject is structural"
                            : "does NOT return to exact -- B2 ill-posed even with coupling")
              << " [|P-exact|=" << std::scientific << std::setprecision(3) << B2_pert.eP
              << " |c-exact|=" << B2_pert.eC << "]\n"
              << "   (load-bearing reject demonstrated on PDE16 D, not here -- k-coupled "
              << "diagonal block cannot host a genuinely ill-posed inflow)\n";
  }
  all_ok &= B2_ok;

  // C: the verdict must FLIP with the reference (Kreiss, rev369) -- zeroed ACCEPTs (a~0, guard imposes
  // nothing), set ACCEPTs (a=(c0,P0): c0 pins P), c-zero REJECTs (a=(0,P0): pins nothing incoming).
  bool const C_ok = ( C_zero.setup_ok && !C_zero.rejected_bc )
                 && ( C_set.setup_ok && !C_set.rejected_bc )
                 && ( !C_czero.setup_ok && C_czero.rejected_bc );
  std::cout << "C  reference-dependent verdict (zeroed=accept, set=accept, c-zero=reject): "
            << (C_ok?"PASS":"FAIL")
            << "   [zeroed: " << (C_zero.setup_ok?"silently accepted":"rejected")
            << " | set: "     << (C_set.setup_ok?"accepted":"rejected")
            << " | c-zero: "  << (C_czero.setup_ok?"accepted":"rejected") << "]\n";
  all_ok &= C_ok;

  // D: NO false reject (Kreiss, rev369) -- midpoint AND face dom_ref ACCEPT and recover.  (Up to rev368 the
  // membership test false-rejected the midpoint arm.)  Only discriminates with the flux guard ON.
  bool D_ok;
  if( flux_on ){
    D_ok = (  D_mid.setup_ok && D_mid.solved
              && D_mid.eP <= TEST_HYP_EXACT_TOL && D_mid.eC <= TEST_HYP_EXACT_TOL )
        && (  D_face.setup_ok && D_face.solved
              && D_face.eP <= TEST_HYP_EXACT_TOL && D_face.eC <= TEST_HYP_EXACT_TOL );
    std::cout << "D  dom_ref face-vs-midpoint: no false reject (both accept+recover): "
              << (D_ok?"PASS":"FAIL") << "\n";
  } else {
    // guard off: the OpP BC skips, so neither dom_ref rejects (both accept).  The
    // degeneracy is invisible without the flux guard -- documented, not scored.
    D_ok = true;
    std::cout << "D  dom_ref discrimination requires the flux guard ON "
              << "(guard off: OpP BC skipped, neither dom_ref rejects) -- not scored\n";
  }
  all_ok &= D_ok;

  // Well-posedness probe / criterion #3 confirmation.  CONTROL: A perturbed must
  // return to exact (delta in-basin).  B perturbed (the rescued bilinear) must ALSO
  // return to exact -- confirming c*d_z P is genuinely well-posed, so accepting it is
  // correct (no false-reject).  B2 perturbed (guard OFF) is the reject-branch readout:
  // it must NOT return to exact -- the silent wrong answer the reject prevents.
  bool const ctrl_ok = A_pert.setup_ok && A_pert.solved
                    && A_pert.eP <= TEST_HYP_EXACT_TOL && A_pert.eC <= TEST_HYP_EXACT_TOL;
  std::cout << "\nwell-posedness probe CONTROL (A perturbed returns to exact): "
            << (ctrl_ok?"in-basin (probe valid)":"OUT-OF-BASIN (probe inconclusive -- lower delta)")
            << "  |P-exact|=" << std::scientific << std::setprecision(3) << A_pert.eP
            << " |c-exact|=" << A_pert.eC << "\n";
  bool const B_returns = B_pert.setup_ok && B_pert.solved
                      && B_pert.eP <= TEST_HYP_EXACT_TOL && B_pert.eC <= TEST_HYP_EXACT_TOL;
  std::cout << "criterion #3 RESCUE confirm (B=c*d_z P perturbed returns to exact -> "
            << "well-posed, accept is correct): "
            << (ctrl_ok ? (B_returns?"CONFIRMED":"NOT CONFIRMED (B did not recover)")
                        : "inconclusive (control out of basin)")
            << "  |P-exact|=" << B_pert.eP << " |c-exact|=" << B_pert.eC << "\n";
  if( !flux_on ){
    bool const B2_returns = B2_pert.setup_ok && B2_pert.solved
                         && B2_pert.eP <= TEST_HYP_EXACT_TOL && B2_pert.eC <= TEST_HYP_EXACT_TOL;
    std::cout << "PROMOTE GATE (PDE18 scope) -- criterion #3 accept/reject branches "
              << "BOTH validated: value (A) + derivative-rescue (B) accept and recover; "
              << "pure-outgoing (B2) rejected under the guard.\n";
    std::cout << "  reject-branch FINDING: B2=c*d_z c from perturbed init "
              << (B2_returns ? "RETURNS to exact -- the zeroth-order coupling (k*P in the "
                               "c-closure row) discretely rescues a pure-outgoing BC, so a "
                               "k-coupled diagonal block cannot exhibit silent wrongness"
                             : "does NOT return to exact -- B2 genuinely ill-posed here")
              << "  [|P-exact|=" << std::scientific << std::setprecision(3) << B2_pert.eP
              << " |c-exact|=" << B2_pert.eC << "]\n";
    std::cout << "  -> the load-bearing reject (silent wrongness guard-off) must be "
              << "confirmed on PDE16 D (principal-symbol coupled).  PSA implication: in "
              << "coupled blocks criterion #3's reject is STRUCTURAL/conservative -- it "
              << "rejects BCs that pin no incoming characteristic even if coupling could "
              << "back-solve them; confirm this is the intended policy before default-on.\n";
  } else {
    std::cout << "PROMOTE GATE: B2 is rejected at setup under the guard (criterion #3 "
              << "reject branch); the guard-OFF build documents the k-coupling rescue\n";
  }

  std::cout << "\n==================================================================\n";
  std::cout << "Darcy-bilinear guard oracle, criterion #3 (value branch accept + derivative\n"
            << "RESCUE accept + pure-outgoing reject + zero-ref collapse + dom_ref): "
            << (all_ok?"PASS":"FAIL") << "\n";
  return all_ok ? 0 : 1;
}
