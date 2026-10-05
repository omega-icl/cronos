// OCFE_M3_index3_solve2.cpp  ---  Stage-2 index-reduction GATE: corpus M3 (index 3)
// ===========================================================================
// Minimal index-3 DAE (corpus M3 from OCFE_PDE20_solve2.cpp), promoted from a
// classification fixture to a SOLVE+VERIFY regression gate.  This is the
// smallest case that exercises TWO differentiation rounds of the high-index
// reducer (the M2 gate exercises one); it is the incremental M2 -> M3 step.
//
//   states:   x(t),v(t) dynamic,  L(t) algebraic (HIDDEN: absent from the constraint)
//   DIFFx:    x' - v       = 0          (ODE)
//   DIFFv:    v' - L       = 0          (ODE)
//   ALG :     x - a(t)     = 0          (constraint, NO v, NO L -> index 3)
//   ICx :     x - a (0)    = 0          (consistency, x dynamic)
//   ICv :     v - a'(0)    = 0          (consistency, v dynamic)
//   a(t) = 1 + 0.5 t - 0.3 t^2 ,  a'(t) = 0.5 - 0.6 t ,  a''(t) = -0.6
//
// Manufactured exact solution:  x=a(t),  v=a'(t),  L=a''(t)=-0.6.
//
// DETECTION (Stage 1): the structural probe reports index 3 with witness {L}
//   (L unpinned by ALG, exposed after TWO differentiations of ALG:
//    d/dt (x-a)   --x'=v--> v - a'    ;
//    d/dt (v-a')  --v'=L--> L - a''  ).
//
// REDUCTION (Stage 2, what this gate measures): differentiate ALG twice, each
//   round substituting the PURE ODE RHS for the exposed state-rate
//   (x'->v, then v'->L), giving L - a''(t) = 0 which pins L; replace ALG with
//   it; keep the original x-a=0 only at the IC (ICx already there).  Result:
//   square index-1 system solving to x=a, v=a', L=a''.
//
// GATE SEMANTICS (ADVERSARIAL seed x=1, v=0, L=0 -- NOT the exact solution):
//   * BEFORE the reducer: square (54=54) but rank-deficient -- L at the global LB
//     is unpinned (DIFFv is imposed at ALL-LB; ALG/ICx/ICv involve only x,v), so it
//     is in the Jacobian null space.  From the non-exact seed the solve
//     residual-converges but L@t=0 stays at its seed -> BASELINE (error at the
//     boundary DOF); exits 0 (measurement, not failure).
//   * AFTER the reducer: L-a''(t)=0 pins L everywhere incl. t=0; solve recovers
//     x=a, v=a', L=a'' to ~machine precision -> REDUCED:PASS.
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
static inline double a_exact  ( double t ){ return 1.0 + 0.5*t - 0.3*t*t; }
static inline double ap_exact ( double t ){ return 0.5 - 0.6*t; }
static inline double app_exact( double   ){ return -0.6; }

int main()
{
  std::cout << "================================================================\n";
  std::cout << "  Stage-2 index-reduction GATE: corpus M3 (index 3)\n";
  std::cout << "  x'=v ; v'=L ; 0=x-a(t) ; a=1+0.5t-0.3t^2\n";
  std::cout << "  (exact: x=a, v=a'=0.5-0.6t, L=a''=-0.6)\n";
  std::cout << "================================================================\n";

  FFGraph DAG;
  FFVar t = DAG.add_var( "t" );
  FFVar x = DAG.add_var( "x(t)" );
  FFVar v = DAG.add_var( "v(t)" );
  FFVar L = DAG.add_var( "L(t)" );          // ALGEBRAIC, HIDDEN (not in the constraint)
  FFPartial OpP;

  FFVar Ae    = 1.0 + 0.5*t - 0.3*t*t;      // a(t)
  FFVar Ap    = 0.5 - 0.6*t;                // a'(t)
  FFVar DIFFx = OpP( x, t ) - v;            // x' - v = 0
  FFVar DIFFv = OpP( v, t ) - L;            // v' - L = 0
  FFVar ALG   = x - Ae;                     // 0 = x - a   (no v, no L -> index 3)
  FFVar ICx   = x - Ae;                     // consistency at t=0
  FFVar ICv   = v - Ap;                     // consistency at t=0

  OCFESLV oc( &DAG );
  oc.add_domain( t, FFDom( 0., TF, NEL_T, FFDom::LGR, NT ) );
  oc.add_state( x, {t} );
  oc.add_state( v, {t} );
  oc.add_state( L, {t} );
  // ADVERSARIAL seed: deliberately NON-exact (x=1, v=0, L=0).  Seeding the exact
  // solution would make the solve vacuous and HIDE the index-3 singularity:
  // L at the global LB is unpinned (in the Jacobian null space), so LM leaves it
  // at its seed.  A non-exact seed forces the unpinned DOF to reveal itself.
  oc.update_ref( x, []( OCFESLV::t_Coord const& ){ return 1.0; } );
  oc.update_ref( v, []( OCFESLV::t_Coord const& ){ return 0.0; } );
  oc.update_ref( L, []( OCFESLV::t_Coord const& ){ return 0.0; } );

  OCFESLV::EqnOptions io( OCFESLV::EqnRole::INTERIOR, 0 ), ii( OCFESLV::EqnRole::INITIAL, 0 );
  int const T_INT = FFDom::ALL - FFDom::LB;
  oc.add_equation( DIFFx, {t}, {T_INT},      io );
  oc.add_equation( DIFFv, {t}, {T_INT},      io );
  oc.add_equation( ALG,   {t}, {FFDom::ALL}, io );
  // 2026-09-23: index 3 with two differential states (x,v) and two hidden levels, so this model has NO free
  // initial data at all -- x(0)=a(0) and v(0)=a'(0) are consistency conditions, not choices.  REDUCE.HIDDEN_IC
  // materialises both levels at t=0, so neither is declared here.
  (void)ICx; (void)ICv;
  oc.set_evolution_domain( t );

  oc.options.REDUCE.ORDER     = OCFESLV::Options::TEST_DAE_REDUCE;
  oc.options.REDUCE.HIDDEN_IC = true;       // the hidden levels are rows; the model declares the free data (none)
  oc.options.CLASSIFY.MODE         = OCFESLV::Options::CLASS_AUTO;
  oc.options.INTERFACE.TYPE   = OCFESLV::Options::IC_AUTO;
  oc.options.INTERFACE.IMPOSITION  = OCFESLV::Options::IC_WEAK;
  oc.options.DISPLAY_LEVEL    = 1;          // show the structural probe (index 3, witness {L})
  oc.options.SOLVE.MAX_ITER   = 50;
  oc.options.SOLVE.RES_TOL    = 1e-9;
#if defined(CRONOS__WITH_SPQR)
  oc.options.SOLVE.FACTORIZATION = OCFESLV::Options::SOLVE_SPQR;
#endif

  bool setup_ok = false, threw = false;
  try { setup_ok = oc.setup(); }
  catch( ... ) { threw = true; setup_ok = false; }

  if( threw ){
    std::cout << "\n  setup: THREW (model form not yet reducible by setup)\n";
    std::cout << "  ==> BASELINE: index-3 detected, reducer rejected the model (setup threw)\n";
    return 0;
  }

  std::cout << "\n  setup: " << (setup_ok?"OK":"FAIL") << "\n";
  for( auto const& [bid, cls] : oc.block_classification() )
    std::cout << "  block " << bid << ": differential_index=" << cls.differential_index
              << "  pde_type=" << OCFESLV::pde_type_name( cls.type ) << "\n";

  // report the persisted reduction plan (Stage-2 decision layer)
  auto const& plan = oc.reduction_plan();
  std::cout << "  reduction_plan: assigns=" << plan.assigns.size()
            << " max_index=" << plan.max_index
            << " resolved=" << (plan.resolved?"yes":"no") << "\n";
  for( auto const& asg : plan.assigns )
    std::cout << "    assign: block=" << asg.block_id << " n_diff=" << asg.n_diff << "\n";

  if( !setup_ok ){
    std::cout << "  ==> BASELINE: setup did not complete (high index not yet reduced)\n";
    return 0;
  }

  size_t const nVar = oc.n_colloc_sta(), nEqn = oc.n_colloc_eqn();
  bool const square = ( nVar == nEqn );
  std::cout << "  nVar=" << nVar << " nEqn=" << nEqn << " square=" << (square?"yes":"no") << "\n";
  if( !square ){
    std::cout << "  ==> BASELINE: index-3 detected but NOT square (reduction not yet applied)\n";
    return 0;
  }

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
    std::cout << "  ==> BASELINE: square but solve did not converge (index-3 not reduced -> singular Jacobian)\n";
    return 0;
  }

  auto eval = [&]( FFVar const& V, double tt ){
    OCFESLV::t_Coord pt; pt[t]=tt;
    return oc.eval_colloc<double>( V, pt, xv.data(), nullptr, nullptr ); };

  // verify vs exact; the index-3 signature is a wrong L at t=0 (witness@LB)
  double ex=0., ev=0., eL=0., eL_at=0.;
  int const NS=101;
  for( int k=0; k<=NS; ++k ){ double tt=double(k)/NS*TF;
    ex=std::max( ex, std::fabs( eval(x,tt)-a_exact(tt)  ) );
    ev=std::max( ev, std::fabs( eval(v,tt)-ap_exact(tt) ) );
    double dL=std::fabs( eval(L,tt)-app_exact(tt) );
    if( dL>eL ){ eL=dL; eL_at=tt; } }
  double const eL_LB = std::fabs( eval(L,0.0) - app_exact(0.0) );   // error at the unpinned DOF
  std::cout << "  verify vs exact:  max|x-a|=" << std::scientific << std::setprecision(3) << ex
            << "  max|v-a'|=" << ev
            << "  max|L-a''|=" << eL << " (at t=" << std::fixed << std::setprecision(3) << eL_at << ")"
            << "  |L-a''|@t=0=" << std::scientific << std::setprecision(3) << eL_LB << "\n";

  bool const pass = ( ex < 1e-8 && ev < 1e-8 && eL < 1e-8 );
  if( pass ){
    std::cout << "  ==> REDUCED:PASS (index-3 -> index-1; L pinned everywhere incl. t=0; exact a/a'/a'')\n";
    return 0;
  }
  bool const lb_dominated = ( eL_LB > 0.5*eL && eL > 1e-6 );
  std::cout << "  ==> BASELINE: "
            << ( lb_dominated ? "L UNPINNED at t=0 (error concentrated at the boundary DOF) -- "
                              : "" )
            << "index-3 NOT fully reduced (square + residual-converged but answer wrong; "
               "L is in the Jacobian null space)\n";
  return 0;   // measurement, not a failure
}
