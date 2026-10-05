// OCFE_SYMBOLCASES1.cpp  ---  principal-symbol classification cases
// =========================================================================
// PURPOSE (2026-09-08).  Two EqnType values are UNEXERCISED corpus-wide:
// MIXED and WEAK_HYPERBOLIC occur 0 times in 1410 classification lines
// (NOTES_20260908k).  STAGE D1/D2 propose to retire/rename them, and that is
// correct either way -- but if the PREDICATES behind them never fire, the
// detection they represent has never been validated and may simply be wrong.
// This driver settles that with models built to hit each branch on purpose.
//
//   u_t + A u_x = s(t,z),   u = (u1,u2),   A constant 2x2
//
// The characteristic speeds are the eigenvalues of A, so A selects the branch:
//   DIAG     A = [[1,0],[0,2]]        real, distinct      -> EVOL_HYPERBOLIC
//   JORDAN   A = [[1,1],[0,1]]        real, DEFECTIVE     -> WEAK_HYPERBOLIC?
//   NEARJOR  A = [[1,1],[1e-8,1]]     real, near-defective-> tolerance probe
//   COMPLEX  A = [[0,-1],[1,0]]       eigenvalues +-i     -> MIXED?
//
// MANUFACTURED SOLUTION, so every case has an exact answer regardless of
// well-posedness:  u1 = 1 + x + t + x^2/2,  u2 = 2 + 2x - t + 3t^2/10,
// with s := u_t + A u_x evaluated from those expressions.  Both are quadratic
// and lie in the collocation space at n_nd >= 3, so a well-posed case should
// reproduce them to round-off.
//
// WHAT IS GATED, AND WHAT IS NOT.
//   GATED: the classification FLAGS -- evolution_hyperbolic and
//     weak_hyperbolic -- because those are stable across STAGE D1/D2, which
//     move WEAK_HYPERBOLIC from the type axis to the flag axis.  A driver
//     gated on the type NAME would have to be edited by the very change it is
//     meant to protect.
//   REPORTED, NOT GATED: the type name, and the solve.  The COMPLEX case is
//     an ill-posed evolution problem (elliptic in the time-like direction);
//     its discrete system may or may not converge and that is not a defect of
//     the interface plan.  Convergence is printed so the reader can see it,
//     never asserted.
// =========================================================================

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <cmath>

#define MC__OCFESLV_SYMBOL_TAG_PROBE

#ifndef OCFE_OCFESLV_HEADER
#define OCFE_OCFESLV_HEADER "ocfeslv.hpp"
#endif
#include OCFE_OCFESLV_HEADER
using namespace mc;

struct Case {
  char const* name;
  double a11, a12, a21, a22;
  char const* expect;      // documented expectation, reported not gated
  bool  expect_evol_hyp;   // GATED
  bool  expect_weak_hyp;   // GATED
};

static Case const CASES[] = {
  { "DIAG   ",  1.0,  0.0, 0.0, 1.0*2.0, "EVOL_HYPERBOLIC (real, distinct)",      true,  false },
  { "JORDAN ",  1.0,  1.0, 0.0, 1.0,     "WEAK_HYPERBOLIC (real, defective)",     true,  true  },
  { "NEARJOR",  1.0,  1.0, 1e-8, 1.0,    "near-defective: tolerance probe",       true,  false },
  { "COMPLEX",  0.0, -1.0, 1.0, 0.0,     "MIXED (complex speeds)",                false, false }
};

static inline double u1M( double t, double x ){ return 1.0 + x + t + 0.5*x*x; }
static inline double u2M( double t, double x ){ return 2.0 + 2.0*x - t + 0.3*t*t; }

static bool run_case( Case const& C, size_t nelt, size_t nelz, size_t n_nd, bool& gate_ok )
{
  std::cout << "\n---- case " << C.name << "  A=[[" << C.a11 << "," << C.a12 << "],["
            << C.a21 << "," << C.a22 << "]]   expect: " << C.expect << " ----\n";

  FFGraph DAG;
  FFVar t = DAG.add_var("t"), x = DAG.add_var("x");
  FFVar u1 = DAG.add_var("u1(t,x)"), u2 = DAG.add_var("u2(t,x)");
  FFPartial OpP;

  // manufactured sources: s = u_t + A u_x, from the exact profiles above
  FFVar const u1_t = 1.0,        u1_x = 1.0 + x;
  FFVar const u2_t = -1.0 + 0.6*t, u2_x = 2.0;
  FFVar const s1 = u1_t + C.a11*u1_x + C.a12*u2_x;
  FFVar const s2 = u2_t + C.a21*u1_x + C.a22*u2_x;

  FFVar E1 = OpP(u1,t) + C.a11*OpP(u1,x) + C.a12*OpP(u2,x) - s1;
  FFVar E2 = OpP(u2,t) + C.a21*OpP(u1,x) + C.a22*OpP(u2,x) - s2;
  FFVar IC1 = u1 - ( 1.0 + x + 0.5*x*x );
  FFVar IC2 = u2 - ( 2.0 + 2.0*x );
  FFVar BC1 = u1 - ( 1.0 + t );
  FFVar BC2 = u2 - ( 2.0 - t + 0.3*t*t );

  OCFESLV oc( &DAG );
  oc.add_domain( t, FFDom( 0., 1., nelt, FFDom::LGL, n_nd ) );
  oc.add_domain( x, FFDom( 0., 1., nelz, FFDom::LGL, n_nd ) );
  oc.add_state( u1, {t,x} );
  oc.add_state( u2, {t,x} );
  oc.update_ref( u1, [&](OCFESLV::t_Coord const& c){ return u1M(c.at(t),c.at(x)); } );
  oc.update_ref( u2, [&](OCFESLV::t_Coord const& c){ return u2M(c.at(t),c.at(x)); } );
  oc.set_evolution_domain( t );

  int const T_NO_LB = FFDom::ALL - FFDom::LB;
  int const X_NO_LB = FFDom::ALL - FFDom::LB;
  oc.add_equation( E1,  {t,x}, {T_NO_LB,   X_NO_LB},   OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( E2,  {t,x}, {T_NO_LB,   X_NO_LB},   OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( IC1, {t,x}, {FFDom::LB, FFDom::ALL}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INITIAL,  0 ) );
  oc.add_equation( IC2, {t,x}, {FFDom::LB, FFDom::ALL}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INITIAL,  0 ) );
  oc.add_equation( BC1, {t,x}, {T_NO_LB,   FFDom::LB},  OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );
  oc.add_equation( BC2, {t,x}, {T_NO_LB,   FFDom::LB},  OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );

  oc.options.REDUCE.ORDER    = OCFESLV::Options::RED_FULL;
  oc.options.CLASSIFY.MODE        = OCFESLV::Options::CLASS_AUTO;
  oc.options.INTERFACE.TYPE  = OCFESLV::Options::IC_AUTO;
  oc.options.INTERFACE.IMPOSITION = OCFESLV::Options::IC_WEAK;
  oc.options.DISPLAY_LEVEL   = 1;

  if( !oc.setup() ){ std::cout << "  setup FAILED (no classification available)\n"; gate_ok = false; return false; }

  OCFESLV::t_Classify const& cls = oc.pde_type();
  std::cout << "  type=" << OCFESLV::pde_type_name( cls.type )
            << "  evolution_hyperbolic=" << cls.evolution_hyperbolic
            << "  weak_hyperbolic=" << cls.weak_hyperbolic
            << "  max_imag_eig=" << std::scientific << std::setprecision(2) << cls.max_imag_eig
            << "  imag_tol=" << cls.imag_tol << "\n";

  bool const flags_ok = ( cls.evolution_hyperbolic == C.expect_evol_hyp )
                     && ( cls.weak_hyperbolic      == C.expect_weak_hyp );
  std::cout << "  FLAG GATE: evolution_hyperbolic/weak_hyperbolic expected "
            << C.expect_evol_hyp << "/" << C.expect_weak_hyp
            << "  ->  " << ( flags_ok ? "PASS" : "** FAIL **" ) << "\n";
  if( !flags_ok ) gate_ok = false;

  // solve is REPORTED, never gated (see header)
  std::vector<double> v, ip;
  if( oc.init( v, ip, nullptr ) ){
    std::vector<double> xv = v;
    for( size_t i = 0; i < xv.size(); ++i ) xv[i] += 1e-3*std::sin( 0.7*double(i) );
    OCFESLV::SolveReport rep = oc.solve( xv.data(), nullptr, nullptr );
    double err = 0.;
    double const sg[3] = { 0.25, 0.5, 0.75 };
    for( double ts : sg ) for( double xs : sg ){
      OCFESLV::t_Coord p; p[t] = ts; p[x] = xs;
      err = std::max( err, std::fabs( oc.eval_colloc<double>( u1, p, xv.data(), nullptr, nullptr ) - u1M(ts,xs) ) );
      err = std::max( err, std::fabs( oc.eval_colloc<double>( u2, p, xv.data(), nullptr, nullptr ) - u2M(ts,xs) ) );
    }
    std::cout << "  solve (reported, NOT gated): conv=" << ( rep.converged ? "y" : "n" )
              << " it=" << rep.iterations << " |r|=" << rep.final_residual
              << "  max|u-u*|=" << err << "\n";
  }
  else std::cout << "  init failed (reported, not gated)\n";
  return true;
}

int main()
{
  std::cout << "==================================================================\n"
            << "  Principal-symbol classification cases: u_t + A u_x = s\n"
            << "  Built to exercise the branches the corpus never reaches:\n"
            << "  defective eigenvectors (weak hyperbolic) and complex speeds.\n"
            << "  GATED on the FLAGS (evolution_hyperbolic, weak_hyperbolic),\n"
            << "  which survive STAGE D1/D2; the type NAME is reported only.\n"
            << "==================================================================\n";
  bool gate_ok = true;
  for( Case const& C : CASES ) run_case( C, 2, 2, 4, gate_ok );
  std::cout << "\n  READ THIS AS:\n"
            << "    all four FLAG GATEs PASS   -> the detection works; MIXED and\n"
            << "        WEAK_HYPERBOLIC are unexercised because the CORPUS lacks\n"
            << "        such models, and STAGE D1/D2 are pure tidying.\n"
            << "    JORDAN weak_hyperbolic=0   -> the defective-eigenvector test\n"
            << "        never fires; the WEAK_HYPERBOLIC branch is dead code and\n"
            << "        the detection, not just the label, needs work.\n"
            << "    COMPLEX evolution_hyperbolic=1 -> complex speeds are being\n"
            << "        reported as hyperbolic, which would be a real defect.\n"
            << "    NEARJOR is a tolerance probe: either verdict is informative,\n"
            << "        but it must not differ from JORDAN by accident -- compare\n"
            << "        max_imag_eig against imag_tol on both rows.\n"
            << "  OCFE_SYMBOLCASES1: " << ( gate_ok ? "PASS" : "FAIL" ) << "\n";
  return gate_ok ? 0 : 1;
}
