// OCFE_M8_coupled_solve2.cpp  ---  Stage-2 GATE: coupled index-2 (TWO witnesses)
// ===========================================================================
// The first corpus case that exercises the BIPARTITE MATCHING of Pantelides
// (not just the differentiation engine).  M2/M3 are single-deficiency single
// chains -- one witness, one constraint -- so the assignment "which constraint
// pins which witness" is trivial.  Here there are TWO hidden algebraic
// witnesses and TWO constraints, and only ONE valid assignment:
//
//   states:   x1,x2 dynamic ;  y1,y2 algebraic (HIDDEN: absent from constraints)
//   DIFFx1:   x1' - y1            = 0
//   DIFFx2:   x2' - y2            = 0
//   g1   :    x1 - a(t)           = 0    (1 diff -> y1 - a'      ; reaches y1 only)
//   g2   :    x1 + x2 - c(t)      = 0    (1 diff -> y1 + y2 - c' ; reaches y1 AND y2)
//   ICx1 :    x1 - a(t)           = 0    @ LB   (consistency)
//   ICx2 :    x1 + x2 - c(t)      = 0    @ LB   (consistency, pins x2(0) given x1(0))
//
//   a(t) = 1 + 0.5 t - 0.3 t^2 ,  a'(t) = 0.5 - 0.6 t
//   c(t) = 2 + 0.4 t + 0.2 t^2 ,  c'(t) = 0.4 + 0.4 t
//
// Manufactured exact solution:
//   x1 = a ,  x2 = c - a ,  y1 = a' ,  y2 = c' - a'
//      = a ,  1 - 0.1t + 0.5t^2 ,  0.5-0.6t ,  -0.1 + 1.0 t
//
// WHY MATCHING IS REQUIRED: y2 is reachable ONLY through g2; y1 through either
// g1 or g2.  The unique valid assignment is y1<-g1, y2<-g2.  A greedy/arbitrary
// pick that assigns g2 to y1 strands y2 (no constraint left), and the current
// header's `constraints.front()` assigns g1 to BOTH -> y2 never pinned.
//
// GATE SEMANTICS (ADVERSARIAL seed x1=x2=1, y1=y2=0):
//   * front()-plan (current header): assign1 differentiates g1 -> y1-a' (pins y1,
//     replaces g1 interior); assign2 ALSO targets g1, whose interior row is now
//     gone -> "interior constraint not found"; g2 untouched -> y2@LB unpinned.
//     Expect BASELINE: max|y2-exact| concentrated at t=0.
//   * augmenting-path plan (after the decision-layer upgrade): y1<-g1, y2<-g2;
//     reduced rows y1-a'=0 and y1+y2-c'=0 pin both witnesses -> REDUCED:PASS.
// ===========================================================================

#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <string>

#ifndef OCFE_OCFESLV_HEADER
#define OCFE_OCFESLV_HEADER "ocfeslv.hpp"
#endif
#include OCFE_OCFESLV_HEADER

#ifndef TEST_DAE_REDUCE
#define TEST_DAE_REDUCE RED_MAIN
#endif

using namespace mc;

static double const TF = 0.5;
static size_t const NEL_T = 3, NT = 6;
static inline double a_e ( double t ){ return 1.0 + 0.5*t - 0.3*t*t; }
static inline double ap_e( double t ){ return 0.5 - 0.6*t; }
static inline double c_e ( double t ){ return 2.0 + 0.4*t + 0.2*t*t; }
static inline double cp_e( double t ){ return 0.4 + 0.4*t; }
static inline double x2_e( double t ){ return c_e(t) - a_e(t); }     // 1 - 0.1t + 0.5t^2
static inline double y1_e( double t ){ return ap_e(t); }             // a'
static inline double y2_e( double t ){ return cp_e(t) - ap_e(t); }   // c' - a' = -0.1 + t

int main()
{
  std::cout << "================================================================\n";
  std::cout << "  Stage-2 GATE: coupled index-2 (TWO witnesses) -- matching test\n";
  std::cout << "  x1'=y1 ; x2'=y2 ; 0=x1-a ; 0=x1+x2-c\n";
  std::cout << "  (exact: x1=a, x2=c-a, y1=a', y2=c'-a')\n";
  std::cout << "================================================================\n";

  FFGraph DAG;
  FFVar t  = DAG.add_var( "t" );
  FFVar x1 = DAG.add_var( "x1(t)" );
  FFVar x2 = DAG.add_var( "x2(t)" );
  FFVar y1 = DAG.add_var( "y1(t)" );        // ALGEBRAIC, HIDDEN
  FFVar y2 = DAG.add_var( "y2(t)" );        // ALGEBRAIC, HIDDEN
  FFPartial OpP;

  FFVar Ae = 1.0 + 0.5*t - 0.3*t*t;         // a(t)
  FFVar Ce = 2.0 + 0.4*t + 0.2*t*t;         // c(t)
  FFVar DIFFx1 = OpP( x1, t ) - y1;         // x1' - y1 = 0
  FFVar DIFFx2 = OpP( x2, t ) - y2;         // x2' - y2 = 0
  FFVar g1     = x1 - Ae;                    // 0 = x1 - a        (reaches y1)
  FFVar g2     = x1 + x2 - Ce;               // 0 = x1 + x2 - c   (reaches y1,y2)
  FFVar ICx1   = x1 - Ae;                    // consistency @ t=0
  FFVar ICx2   = x1 + x2 - Ce;               // consistency @ t=0

  OCFESLV oc( &DAG );
  oc.add_domain( t, FFDom( 0., TF, NEL_T, FFDom::LGR, NT ) );
  oc.add_state( x1, {t} );
  oc.add_state( x2, {t} );
  oc.add_state( y1, {t} );
  oc.add_state( y2, {t} );
  // ADVERSARIAL seed (non-exact): forces unpinned witnesses to reveal themselves.
  oc.update_ref( x1, []( OCFESLV::t_Coord const& ){ return 1.0; } );
  oc.update_ref( x2, []( OCFESLV::t_Coord const& ){ return 1.0; } );
  oc.update_ref( y1, []( OCFESLV::t_Coord const& ){ return 0.0; } );
  oc.update_ref( y2, []( OCFESLV::t_Coord const& ){ return 0.0; } );

  OCFESLV::EqnOptions io( OCFESLV::EqnRole::INTERIOR, 0 ), ii( OCFESLV::EqnRole::INITIAL, 0 );
  int const T_INT = FFDom::ALL - FFDom::LB;
  oc.add_equation( DIFFx1, {t}, {T_INT},      io );
  oc.add_equation( DIFFx2, {t}, {T_INT},      io );
  oc.add_equation( g1,     {t}, {FFDom::ALL}, io );
  oc.add_equation( g2,     {t}, {FFDom::ALL}, io );
  // 2026-09-23: index 2 through TWO constraints, so two hidden levels against two differential states (x1,x2):
  // this model has NO free initial data -- both ICs were consistency conditions.  REDUCE.HIDDEN_IC materialises
  // g1 and g2 at t=0, which is also the stronger matching test: each hidden level must pin its own witness.
  (void)ICx1; (void)ICx2;
  oc.set_evolution_domain( t );

  oc.options.REDUCE.ORDER     = OCFESLV::Options::TEST_DAE_REDUCE;
  oc.options.REDUCE.HIDDEN_IC = true;       // the hidden levels are rows; the model declares the free data (none)
  oc.options.CLASSIFY.MODE         = OCFESLV::Options::CLASS_AUTO;
  oc.options.INTERFACE.TYPE   = OCFESLV::Options::IC_AUTO;
  oc.options.INTERFACE.IMPOSITION  = OCFESLV::Options::IC_WEAK;
  oc.options.DISPLAY_LEVEL    = 1;          // show the structural probe (index 2, witnesses {y1,y2})
  oc.options.SOLVE.MAX_ITER   = 50;
  oc.options.SOLVE.RES_TOL    = 1e-9;
#if defined(CRONOS__WITH_SPQR)
  oc.options.SOLVE.FACTORIZATION = OCFESLV::Options::SOLVE_SPQR;
#endif

  bool setup_ok = false, threw = false;
  try { setup_ok = oc.setup(); }
  catch( ... ) { threw = true; setup_ok = false; }
  if( threw ){
    std::cout << "\n  setup: THREW\n  ==> BASELINE: coupled index-2 detected, setup threw\n";
    return 0;
  }

  std::cout << "\n  setup: " << (setup_ok?"OK":"FAIL") << "\n";
  for( auto const& [bid, cls] : oc.block_classification() )
    std::cout << "  block " << bid << ": differential_index=" << cls.differential_index
              << "  pde_type=" << OCFESLV::pde_type_name( cls.type ) << "\n";

  auto const& plan = oc.reduction_plan();
  std::cout << "  reduction_plan: assigns=" << plan.assigns.size()
            << " max_index=" << plan.max_index
            << " resolved=" << (plan.resolved?"yes":"no") << "\n";
  for( auto const& asg : plan.assigns )
    std::cout << "    assign: block=" << asg.block_id << " n_diff=" << asg.n_diff << "\n";

  if( !setup_ok ){ std::cout << "  ==> BASELINE: setup did not complete\n"; return 0; }

  size_t const nVar = oc.n_colloc_sta(), nEqn = oc.n_colloc_eqn();
  bool const square = ( nVar == nEqn );
  std::cout << "  nVar=" << nVar << " nEqn=" << nEqn << " square=" << (square?"yes":"no") << "\n";
  if( !square ){ std::cout << "  ==> BASELINE: detected but NOT square\n"; return 0; }

  std::vector<double> varInit, inpInit;
  if( !oc.init( varInit, inpInit, nullptr ) ){ std::cerr << "  init failed\n"; return 1; }
  std::vector<double> xv = varInit;
  OCFESLV::SolveReport rep = oc.solve( xv.data(), nullptr, nullptr );
  std::cout << "  solve: converged=" << (rep.converged?"yes":"no")
            << " iters=" << rep.iterations << " final|r|=" << std::scientific
            << std::setprecision(3) << rep.final_residual << "\n";
  if( rep.iterations == 0 )
    std::cout << "  !! WARNING: 0 iterations -- seed was already a fixed point; test is VACUOUS\n";
  if( !rep.converged ){
    std::cout << "  ==> BASELINE: square but solve did not converge\n";
    return 0;
  }

  auto eval = [&]( FFVar const& V, double tt ){
    OCFESLV::t_Coord pt; pt[t]=tt;
    return oc.eval_colloc<double>( V, pt, xv.data(), nullptr, nullptr ); };

  double e_x1=0,e_x2=0,e_y1=0,e_y2=0,e_y2_at=0;
  int const NS=101;
  for( int k=0; k<=NS; ++k ){ double tt=double(k)/NS*TF;
    e_x1=std::max(e_x1,std::fabs(eval(x1,tt)-a_e(tt)));
    e_x2=std::max(e_x2,std::fabs(eval(x2,tt)-x2_e(tt)));
    e_y1=std::max(e_y1,std::fabs(eval(y1,tt)-y1_e(tt)));
    double d2=std::fabs(eval(y2,tt)-y2_e(tt)); if(d2>e_y2){ e_y2=d2; e_y2_at=tt; } }
  double const e_y1_LB=std::fabs(eval(y1,0.0)-y1_e(0.0));
  double const e_y2_LB=std::fabs(eval(y2,0.0)-y2_e(0.0));
  std::cout << "  verify vs exact:  max|x1-a|=" << std::scientific << std::setprecision(3) << e_x1
            << "  max|x2-(c-a)|=" << e_x2 << "\n"
            << "                    max|y1-a'|=" << e_y1 << " (LB " << e_y1_LB << ")"
            << "  max|y2-(c'-a')|=" << e_y2 << " (at t=" << std::fixed << std::setprecision(3) << e_y2_at
            << ", LB " << std::scientific << std::setprecision(3) << e_y2_LB << ")\n";

  bool const pass = ( e_x1<1e-8 && e_x2<1e-8 && e_y1<1e-8 && e_y2<1e-8 );
  if( pass ){
    std::cout << "  ==> REDUCED:PASS (coupled index-2 -> index-1; both witnesses pinned; matching correct)\n";
    return 0;
  }
  std::cout << "  ==> BASELINE: coupled index-2 NOT fully reduced "
            << "(witness(es) unpinned -- assignment requires augmenting-path matching, "
               "not constraints.front())\n";
  return 0;   // measurement, not a failure
}
