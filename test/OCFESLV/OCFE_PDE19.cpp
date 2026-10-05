// ===========================================================================
// OCFE_PDE19_solve2.cpp
//
// PSA case study, step 2: the ALGEBRAIC-momentum coupled first-order (c,u)
// block -- the pressure/velocity subsystem in its cleanest isolating form.
//
// PDE15 (stage 3) already validated the NON-SYMMETRIC off-diagonal read, the
// rowspace AUTO-closure, and the incoming-BC direction guard -- but with BOTH
// states evolution (d_t c, d_t u).  PSA's real (P,u) block is NOT all-evolution:
// the total molar balance is differential in the advected density, while Ergun
// momentum is ALGEBRAIC (no d_t) -- a differential-algebraic (DAE) pair whose
// z-principal symbol is off-diagonal.  No existing driver declares an algebraic
// state inside a coupled first-order block; this driver is the first.
//
// Linear-surrogate model (gradual step; static-per-step flow direction):
//
//   continuity (c differential):  d_t c + d_z(u c) = f_c        on (t,z)
//   momentum   (u ALGEBRAIC):     d_z c + R u      = g          on (t,z)   [no d_t u]
//
//   The bilinear advection d_z(u c) linearises about (u0,c0) to
//       u0 d_z dc + c0 d_z du                                   (off-diagonal),
//   and the algebraic momentum contributes d_z dc (with u algebraic).  The
//   z-principal symbol is therefore
//       A_z = [[u0, c0],
//              [ 1,  0 ]]        det = -c0 < 0
//   so the eigenvalues  lam_+- = ( u0 +- sqrt(u0^2 + 4 c0) ) / 2  are ALWAYS one
//   positive and one negative: a FIXED 1+1 incoming split (one condition at each
//   end), with the eigenVECTORS rotating as u0 changes sign.  The "upwind end
//   flips with sgn(u)" intuition is the scalar-advection view; the coupled-block
//   COUNT does not flip -- which is itself a thing to confirm (reference-
//   robustness probe: 1+1 invariant under the fwd/rev reference swap).
//
//   A per-(state,dir) scan cannot recover this: the momentum row differentiates
//   c, not u, so no diagonal "equation i owns state i" assignment exists -- the
//   joint eigenstructure of A_z is required, and the Hungarian continuity-slot
//   matching must cope with a 2-eq/2-state block that is not diagonally
//   assignable.
//
// WHAT THIS DRIVER PROBES (first cut: IC_WEAK, forward + reverse):
//   1. Does the framework accept an algebraic state (u: no d_t, no IC) inside a
//      coupled first-order block, and classify the (c,u) block correctly?
//   2. Does it read the off-diagonal 1+1 z-split and place one incoming
//      condition at each end (criterion #3 accepts the physical c-value BCs)?
//   3. Does the AUTO outgoing-characteristic closure fire on the algebraic block,
//      leaving a SQUARE system?
//   4. fwd (u>0) and rev (u<0) as two STATIC setups: both well-posed, recover the
//      manufactured solution, reference-robustness reports the 1+1 split invariant.
//
//   This is DIAGNOSTIC-FIRST: classification / squareness / BC directionality /
//   robustness are printed BEFORE the solve, so a setup that needs iteration
//   still tells us how the framework reads the new structure.  IC_TRACE/IC_STRONG,
//   nonlinear Ergun + ideal-gas EOS, and dynamic mid-step reversal are deferred
//   follow-ups once this WEAK static cut is confirmed.
//
// EXPLORATORY ASSUMPTIONS (flagged; the first log adjudicates):
//   - u is marked algebraic purely by OMITTING d_t u and its IC (no explicit
//     "algebraic state" API used).
//   - the algebraic momentum is collocated on {ALL t, z-interior} (the constraint
//     holds at every t, including t=LB), while continuity is on {t>LB, z-interior}.
//   - physical c-value BCs at BOTH ends (feed density at inflow, back-pressure at
//     outflow via the linear EOS P=c); each is expected well-directed because the
//     incoming right-eigenvector has a non-zero c-component at either end.
//   - the framework AUTO-closure (default) supplies the outgoing rows.
//
// Build:  + -DCRONOS__WITH_SPQR (QR), + -DCRONOS__WITH_UMFPACK (IC_STRONG, later).
//   -DMC__OCFESLV_REFERENCE_ROBUSTNESS_PROBE  recommended for this first run.
//   Knobs:  -DTEST_HYP_NEL_T/_NEL_Z (3/3), -DTEST_HYP_NT/_NZ (6/10),
//           -DTEST_HYP_R (Darcy resistance, default 1.5), -DTEST_HYP_TF/_ZF.
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
#ifndef TEST_HYP_R
#define TEST_HYP_R 1.5      // linear Darcy resistance:  d_z c + R u = g
#endif

// ---- LINK-detection probe (scaffold, opt-in) --------------------------------
// LOAD-BEARING configuration (default): the momentum equation `d_z c + R u = g`
// declares u as a GENUINE algebraic state (no d_t, no IC, EqnRole::INTERIOR), and
// the framework's algebraic-boundary closure (MC__OCFESLV_AUTO_ALG_CLOSURE) supplies
// the z=LB/UB constraint rows that the {ALL t, z-interior} collocation omits --
// closing the block to square exactly as the M4 twin in the PDE20 corpus does.
//
// PROBE (opt-in, -DTEST_PDE19_LINK_MOMENTUM=1): HAND-TAG the momentum equation as a
// first-order LINK (the role reduce_order would assign to an auxiliary = derivative
// + LINK equation) to measure how far the LINK/DESCRIPTOR->PARABOLIC machinery gets
// when the link is user-supplied (u is NOT in the _auxDef registry).  This is a
// diagnostic scaffold only -- it is NOT the load-bearing path: an EqnRole::LINK is
// treated as an order-reduction auxiliary and is therefore (correctly) skipped by
// the algebraic-boundary closure, leaving the block 36-short.  The eventual goal is
// to DETECT such linking algebraic equations internally; until then, the genuine-
// algebraic declaration above IS the supported route.
#ifndef TEST_PDE19_LINK_MOMENTUM
#define TEST_PDE19_LINK_MOMENTUM 0
#endif
#ifndef TEST_HYP_TF
#define TEST_HYP_TF 0.5
#endif
#ifndef TEST_HYP_ZF
#define TEST_HYP_ZF 0.5
#endif

// ---------------------------------------------------------------------------
// Manufactured solution (polynomials in z, linear-in-t growth) so the spectral
// d_z is EXACT and recovery is at machine accuracy -- the observable is the
// SETUP (classification / placement / squareness), not aliasing error.
//   C(z)   = c0 + c1 z + c2 z^2 + c3 z^3        (kept > 0 on [0,zf])
//   Up(z)  = u0 + u1 z + u2 z^2                 (kept > 0 on [0,zf]; signed by dir)
//   c_exact(t,z) = C(z)  (1 + alp t)
//   u_exact(t,z) = dir Up(z) (1 + alp t)        dir = +1 forward, -1 reverse
// ---------------------------------------------------------------------------
struct Par {
  double tf  = TEST_HYP_TF;
  double zf  = TEST_HYP_ZF;
  double R   = TEST_HYP_R;
  double alp = 0.3;                      // linear-in-time growth
  double c0=1.0, c1=0.4, c2=-0.3, c3=0.1;   // C(z): C(0)=1, C(zf=.5)=1.1375 > 0
  double u0=0.8, u1=0.3, u2=-0.2;           // Up(z): Up(0)=0.8, Up(.5)=0.9 > 0
};

static double Cpoly ( double z, Par const& p ){ return p.c0 + p.c1*z + p.c2*z*z + p.c3*z*z*z; }
static double Uppoly( double z, Par const& p ){ return p.u0 + p.u1*z + p.u2*z*z; }
static double C_exact( double t, double z, Par const& p ){ return Cpoly(z,p)*(1.0+p.alp*t); }
static double U_exact( double t, double z, Par const& p, double dir ){ return dir*Uppoly(z,p)*(1.0+p.alp*t); }

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

struct ModeResult { std::string name; bool ok=true; bool rejected_bc=false;
                    bool square=false; size_t nVar=0,nEqn=0,nTrace=0;
                    double final_res=0., eC=0., eU=0.; };

// dir = +1 forward (u>0, inflow at z=LB), -1 reverse (u<0, inflow at z=UB).
static ModeResult run_step( OCFESLV::Options::ImpositionType imp,
                            std::string const& dir_tag, Par const& p, double dir )
{
  ModeResult R; R.name = std::string(OCFESLV::interface_imposition_name(imp)) + " / " + dir_tag;
  bool ok=true;

  FFGraph DAG;
  FFVar t = DAG.add_var("t");
  FFVar z = DAG.add_var("z");
  FFVar c = DAG.add_var("c(t,z)");
  FFVar u = DAG.add_var("u(t,z)");
  FFPartial OpP;

  // z-polynomial building blocks as DAG expressions.
  FFVar Cexpr  = p.c0 + p.c1*z + p.c2*z*z + p.c3*z*z*z;
  FFVar Czexpr = p.c1 + 2.0*p.c2*z + 3.0*p.c3*z*z;
  FFVar Upexpr = p.u0 + p.u1*z + p.u2*z*z;
  FFVar Upzexpr= p.u1 + 2.0*p.u2*z;
  FFVar tfac   = 1.0 + p.alp*t;

  FFVar CE = Cexpr*tfac;                 // c_exact
  FFVar UE = dir*Upexpr*tfac;            // u_exact (signed)

  // f_c = d_t c + d_z(u c) = alp C + dir (1+alp t)^2 ( Up' C + Up C' )
  FFVar FC = p.alp*Cexpr + dir*tfac*tfac*( Upzexpr*Cexpr + Upexpr*Czexpr );
  // g   = d_z c + R u     = (1+alp t)( C' + R dir Up )
  FFVar G  = tfac*( Czexpr + p.R*dir*Upexpr );

  FFVar PDEC = OpP(c,t) + u*OpP(c,z) + c*OpP(u,z) - FC;   // continuity (bilinear, off-diagonal)
  FFVar PDEU = OpP(c,z) + p.R*u - G;                      // momentum (ALGEBRAIC, linear Darcy)
  FFVar ICC  = c - CE;                                    // initial c (only c -- u algebraic)
  FFVar BCC  = c - CE;                                    // c-value BC (well-directed at either end)

  OCFESLV oc(&DAG);
  oc.add_domain( t, FFDom(0., p.tf, TEST_HYP_NEL_T, FFDom::LGR, TEST_HYP_NT) );
  oc.add_domain( z, FFDom(0., p.zf, TEST_HYP_NEL_Z, FFDom::LGL, TEST_HYP_NZ) );
  oc.add_state( c, {t,z} );
  oc.add_state( u, {t,z} );
  oc.update_ref( c, [&]( OCFESLV::t_Coord const& crd ){ return C_exact(crd.at(t),crd.at(z),p); } );
  oc.update_ref( u, [&]( OCFESLV::t_Coord const& crd ){ return U_exact(crd.at(t),crd.at(z),p,dir); } );

  OCFESLV::EqnOptions int_opt( OCFESLV::EqnRole::INTERIOR, 0 );
  OCFESLV::EqnOptions ini_opt( OCFESLV::EqnRole::INITIAL,  0 );
  OCFESLV::EqnOptions bnd_opt( OCFESLV::EqnRole::BOUNDARY, 0 );

  int const T_INT = FFDom::ALL - FFDom::LB;
  int const Z_INT = FFDom::ALL - FFDom::LB - FFDom::UB;

  // continuity (differential c): t-interior, z-interior.
  oc.add_equation( PDEC, {t,z}, {T_INT,     Z_INT},     int_opt );
  // momentum (algebraic u): ALL t (the constraint holds at every t, incl. LB), z-interior.
  // PROBE: optionally hand-tag it as a first-order LINK (the role reduce_order would
  // assign) to measure the parabolic-closure path for a user-supplied linking eqn.
#if TEST_PDE19_LINK_MOMENTUM
  OCFESLV::EqnOptions mom_opt( OCFESLV::EqnRole::LINK, 0 );
#else
  OCFESLV::EqnOptions mom_opt = int_opt;
#endif
  oc.add_equation( PDEU, {t,z}, {FFDom::ALL, Z_INT},    mom_opt );
  // initial data for c ONLY (u is algebraic -> no IC).
  oc.add_equation( ICC,  {t,z}, {FFDom::LB,  FFDom::ALL}, ini_opt );
  // physical c-value BC at each end (1+1 incoming); AUTO closure supplies outgoing.
  oc.add_equation( BCC,  {t,z}, {T_INT,     FFDom::LB}, bnd_opt );
  oc.add_equation( BCC,  {t,z}, {T_INT,     FFDom::UB}, bnd_opt );

  oc.set_evolution_domain( t );
  oc.options.REDUCE.ORDER    = OCFESLV::Options::RED_FULL;
  oc.options.CLASSIFY.MODE        = OCFESLV::Options::CLASS_AUTO;
  oc.options.INTERFACE.TYPE  = OCFESLV::Options::IC_AUTO;
  oc.options.INTERFACE.IMPOSITION = imp;
  // Genuine-algebraic momentum (TEST_PDE19_LINK_MOMENTUM=0): the algebraic-state
  // boundary closure is load-bearing here.  Was -DMC__OCFESLV_AUTO_ALG_CLOSURE,
  // now a runtime Option.
  oc.options.INTERFACE.SAT_SIGMA0      = 10.0;
  oc.options.DISPLAY_LEVEL   = 1;
  oc.options.SOLVE.MAX_ITER  = 50;
  oc.options.SOLVE.RES_TOL   = 1e-9;
#ifdef CRONOS__WITH_SPQR
  oc.options.SOLVE.FACTORIZATION = OCFESLV::Options::SOLVE_SPQR;
#endif

#ifdef MC__OCFESLV_REFERENCE_ROBUSTNESS_PROBE
  // Two references straddling the fwd/rev sign: the per-face incoming/outgoing
  // split must be invariant (1+1) -- this is the count-invariance check.
  oc.add_classification_reference_sample( { {c, p.c0}, {u,  0.8} } );
  oc.add_classification_reference_sample( { {c, p.c0}, {u, -0.8} } );
#endif

  std::cout << "\n------ PDE19 (c,u) algebraic off-diagonal  [" << R.name << "] ------\n";
  std::cout << "PDE: d_t c + d_z(u c) = f_c ;  d_z c + R u = g   (u ALGEBRAIC), R="
            << p.R << "  dir=" << (dir>0?"+1 (fwd, u>0)":"-1 (rev, u<0)") << "\n";

  bool const setup_ok = oc.setup();

  // ---- DIAGNOSTIC-FIRST: report the classification/placement read regardless ----
  {
    auto const pt = oc.pde_type();
    std::cout << "pde_type: " << OCFESLV::pde_type_name( pt.type )
              << "  evolution_hyperbolic=" << (pt.evolution_hyperbolic?"yes":"no") << "\n";
  }
  std::cout << "setup: " << (setup_ok?"OK":"FAILED")
            << "  status=" << OCFESLV::setup_status_str( oc.setup_status() ) << "\n";
  R.rejected_bc = ( oc.setup_status() == OCFESLV::SetupStatus::HYP_INCOMING_BC
                 || oc.setup_status() == OCFESLV::SetupStatus::HYP_BC_MISDIRECTED );

  if( !setup_ok ){
    std::cout << "  -> setup did not complete; see status above (this is the structure-learning signal).\n";
    R.ok=false; return R;
  }

  size_t const nVar=oc.n_colloc_sta(), nEqn=oc.n_colloc_eqn(), nTrace=oc.n_colloc_trace();
  R.nVar=nVar; R.nEqn=nEqn; R.nTrace=nTrace; R.square=(nVar==nEqn);
  std::cout << "nVar=" << nVar << " nEqn=" << nEqn << " nTrace=" << nTrace
            << " square=" << (R.square?"yes":"no") << "\n";
  ok &= check_close( "square system (nVar==nEqn)", R.square?0.0:1.0, 0.0 );

#ifdef MC__OCFESLV_REFERENCE_ROBUSTNESS_PROBE
  bool const rob = oc.reference_robustness_ok();
  std::cout << "reference-robustness (1+1 split invariant under u-sign swap): "
            << (rob?"OK":"VIOLATED") << "\n";
  ok &= rob;
#endif

  // ---- solve + manufactured-accuracy recovery ----
  std::vector<double> var( nVar, 0.0 );
  if( !oc.init( var.data(), nullptr, nullptr ) ){
    std::cout << "  init() failed\n"; R.ok=false; return R;
  }
  auto rep = oc.solve( var.data() );
  std::cout << "solve: converged=" << (rep.converged?"yes":"no")
            << " iters=" << rep.iterations
            << " final_residual=" << std::scientific << std::setprecision(3) << rep.final_residual << "\n";

  // residual via eval (mirror PDE15), used for the PASS check.
  std::vector<double> res( nEqn, 0.0 );
  oc.eval( res.data(), nullptr, var.data(), nullptr, nullptr );
  R.final_res = max_abs( res );
  ok &= check_close( "final max residual", R.final_res, 1e-7 );

  // per-state error vs manufactured (state order follows add_state: 0=c, 1=u).
  {
    double eC=0.0, eU=0.0; size_t off=0, sidx=0;
    for( auto const& st : oc.states_colloc() ){
      auto nodes = oc.node_colloc( st );
      for( size_t i=0;i<nodes.size();++i ){
        double const ex = (sidx==0) ? C_exact(nodes[i][0],nodes[i][1],p)
                                    : U_exact(nodes[i][0],nodes[i][1],p,dir);
        double const e = std::fabs( var[off+i] - ex );
        if( sidx==0 ) eC=std::max(eC,e); else eU=std::max(eU,e);
      }
      off += nodes.size(); ++sidx;
    }
    R.eC=eC; R.eU=eU;
    ok &= check_close( "max |c - c_exact|", eC, 1e-6 );
    ok &= check_close( "max |u - u_exact|", eU, 1e-6 );
  }

  R.ok=ok;
  std::cout << "[" << R.name << "]: " << (ok?"PASS":"FAIL") << "\n";
  return R;
}

int main()
{
  std::cout << "PDE19 build config: closure=AUTO  base-guard="
            << "ON"
            << "  flux-guard="
            << "ON"
            << "  robustness-probe="
#ifdef MC__OCFESLV_REFERENCE_ROBUSTNESS_PROBE
            << "ON"
#else
            << "off"
#endif
            << "  momentum-role="
#if TEST_PDE19_LINK_MOMENTUM
            << "LINK (probe)"
#else
            << "INTERIOR (genuine-algebraic, load-bearing)"
#endif
            << "\n";

  Par p;
  bool all_ok = true;

  // First cut: IC_WEAK only (direct solve, no interface augmentation), both
  // flow directions as static per-step setups.
  std::cout << "\n========== PDE19: algebraic-momentum (c,u) off-diagonal block ==========\n";
  ModeResult const fwd = run_step( OCFESLV::Options::IC_WEAK, "forward", p, +1.0 );
  ModeResult const rev = run_step( OCFESLV::Options::IC_WEAK, "reverse", p, -1.0 );
  all_ok &= fwd.ok && rev.ok;

  std::cout << "\n==================================================================\n";
  std::cout << "PDE19 algebraic-momentum off-diagonal oracle (IC_WEAK, fwd+rev): "
            << (all_ok?"PASS":"FAIL") << "\n";
  return all_ok ? 0 : 1;
}
