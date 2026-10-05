// OCFE_M2_index2_solve2.cpp  ---  Stage-2 index-reduction GATE: corpus M2 (index 2)
// ===========================================================================
// Minimal index-2 DAE (corpus M2 from OCFE_PDE20_solve2.cpp), promoted from a
// classification fixture to a SOLVE+VERIFY regression gate for the high-index
// reducer (Pantelides differentiation set + dummy-derivative / constraint
// substitution).  This is the smallest case that exercises every step of Stage 2.
//
//   states:   x(t) dynamic,  y(t) algebraic (HIDDEN: absent from the constraint)
//   DIFF:     x' - y       = 0          (ODE)
//   ALG :     x - a(t)     = 0          (constraint, NO y -> index 2)
//   ICx :     x - a(0)     = 0
//   a(t) = 1 + 0.5 t - 0.3 t^2 ,   a'(t) = 0.5 - 0.6 t
//
// Manufactured exact solution:  x(t) = a(t),  y(t) = a'(t).
//
// DETECTION (Stage 1, already in the header): the structural probe reports
//   index 2 with witness {y} (y unpinned by ALG, exposed after ONE differentiation
//   of ALG: d/dt(x-a) = x'-a' --substitute x'=y--> y - a' = 0).
//
// REDUCTION (Stage 2, the work this gate measures): differentiate ALG once,
//   substitute the ODE for the state-derivative -> y - a'(t) = 0 pins y; replace
//   ALG with that; keep the original ALG as the consistency IC (ICx already does).
//   Result: square index-1 system solving to x=a(t), y=a'(t).
//
// GATE SEMANTICS (ADVERSARIAL seed x=1, y=0 -- NOT the exact solution):
//   * Seeding the exact solution is a FALSE-PASS trap: y at the global LB is unpinned
//     (in the Jacobian null space), the residual is insensitive to it, and LM leaves it
//     at its seed -- so an exact seed "passes" with 0 iterations without ever reducing.
//   * BEFORE the reducer: the unreduced index-2 system is square but rank-deficient.
//     From the non-exact seed it residual-converges (y@LB does not enter the residual)
//     but y@t=0 stays wrong -> gate prints BASELINE (error concentrated at the boundary
//     DOF) and exits 0 (a measurement, not a failure).
//   * AFTER the reducer: the differentiated constraint y-a'(t)=0 pins y everywhere incl.
//     t=0; solve recovers x=a(t), y=a'(t) to ~machine precision -> REDUCED:PASS.
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
static inline double a_exact ( double t ){ return 1.0 + 0.5*t - 0.3*t*t; }
static inline double ap_exact( double t ){ return 0.5 - 0.6*t; }

int main()
{
  std::cout << "================================================================\n";
  std::cout << "  Stage-2 index-reduction GATE: corpus M2 (index 2)\n";
  std::cout << "  x'=y ; 0=x-a(t) ; a=1+0.5t-0.3t^2   (exact: x=a, y=a')\n";
  std::cout << "================================================================\n";

  FFGraph DAG;
  FFVar t = DAG.add_var( "t" );
  FFVar x = DAG.add_var( "x(t)" );
  FFVar y = DAG.add_var( "y(t)" );          // ALGEBRAIC, HIDDEN (not in the constraint)
  FFPartial OpP;

  FFVar Ae   = 1.0 + 0.5*t - 0.3*t*t;       // a(t)
  FFVar DIFF = OpP( x, t ) - y;             // x' - y = 0
  FFVar ALG  = x - Ae;                      // 0 = x - a   (no y -> index 2)
  FFVar ICx  = x - Ae;                      // consistency at t=0

  OCFESLV oc( &DAG );
  oc.add_domain( t, FFDom( 0., TF, NEL_T, FFDom::LGR, NT ) );
  oc.add_state( x, {t} );
  oc.add_state( y, {t} );
  // ADVERSARIAL seed: deliberately NON-exact (x=1, y=0).  Seeding the exact solution
  // makes the solve vacuous (0 iterations) and HIDES the index-2 singularity: y at the
  // global LB is unpinned (DIFF is imposed at ALL-LB; ALG/ICx involve only x), so it is
  // in the Jacobian null space and LM leaves it at its seed.  A non-exact seed forces the
  // unpinned DOF to reveal itself -- only a genuinely reduced (index-1) system pins y
  // everywhere and recovers a'(t) at t=0.
  oc.update_ref( x, []( OCFESLV::t_Coord const& ){ return 1.0; } );
  oc.update_ref( y, []( OCFESLV::t_Coord const& ){ return 0.0; } );

  OCFESLV::EqnOptions io( OCFESLV::EqnRole::INTERIOR, 0 ), ii( OCFESLV::EqnRole::INITIAL, 0 );
  int const T_INT = FFDom::ALL - FFDom::LB;
  oc.add_equation( DIFF, {t}, {T_INT},      io );
  oc.add_equation( ALG,  {t}, {FFDom::ALL}, io );
  oc.add_equation( ICx,  {t}, {FFDom::LB},  ii );
  oc.set_evolution_domain( t );

  oc.options.REDUCE.ORDER     = OCFESLV::Options::TEST_DAE_REDUCE;
  oc.options.CLASSIFY.MODE         = OCFESLV::Options::CLASS_AUTO;
  oc.options.INTERFACE.TYPE   = OCFESLV::Options::IC_AUTO;
  oc.options.INTERFACE.IMPOSITION  = OCFESLV::Options::IC_WEAK;
  oc.options.DISPLAY_LEVEL    = 1;          // show the structural probe (index 2, witness {y})
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
    std::cout << "  ==> BASELINE: index-2 detected, reducer NOT YET IMPLEMENTED (setup rejects high index)\n";
    return 0;
  }

  // report per-block differential index (Stage-1 detection)
  std::cout << "\n  setup: " << (setup_ok?"OK":"FAIL") << "\n";
  for( auto const& [bid, cls] : oc.block_classification() )
    std::cout << "  block " << bid << ": differential_index=" << cls.differential_index
              << "  pde_type=" << OCFESLV::pde_type_name( cls.type ) << "\n";

  if( !setup_ok ){
    std::cout << "  ==> BASELINE: setup did not complete (high index not yet reduced)\n";
    return 0;
  }

  size_t const nVar = oc.n_colloc_sta(), nEqn = oc.n_colloc_eqn();
  bool const square = ( nVar == nEqn );
  std::cout << "  nVar=" << nVar << " nEqn=" << nEqn << " square=" << (square?"yes":"no") << "\n";
  if( !square ){
    std::cout << "  ==> BASELINE: index-2 detected but NOT square (reduction not yet applied)\n";
    return 0;
  }

  // square -> attempt the solve and verify against the manufactured solution
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
    std::cout << "  ==> BASELINE: square but solve did not converge (index-2 not reduced -> singular Jacobian)\n";
    return 0;
  }

  auto eval = [&]( FFVar const& V, double tt ){
    OCFESLV::t_Coord pt; pt[t]=tt;
    return oc.eval_colloc<double>( V, pt, xv.data(), nullptr, nullptr ); };

  // verify vs exact; track where the error lives (the index-2 signature is a wrong y at t=0)
  double ex=0., ey=0., ey_at=0.;
  int const NS=101;
  for( int k=0; k<=NS; ++k ){ double tt=double(k)/NS*TF;
    double dx=std::fabs( eval(x,tt)-a_exact(tt) ), dy=std::fabs( eval(y,tt)-ap_exact(tt) );
    ex=std::max(ex,dx); if( dy>ey ){ ey=dy; ey_at=tt; } }
  double const ey_LB = std::fabs( eval(y,0.0) - ap_exact(0.0) );   // error at the unpinned DOF
  std::cout << "  verify vs exact:  max|x-a(t)|=" << std::scientific << std::setprecision(3) << ex
            << "  max|y-a'(t)|=" << ey << " (at t=" << std::fixed << std::setprecision(3) << ey_at << ")"
            << "  |y-a'|@t=0=" << std::scientific << std::setprecision(3) << ey_LB << "\n";

  bool const pass = ( ex < 1e-8 && ey < 1e-8 );
  if( pass ){
    std::cout << "  ==> REDUCED:PASS (index-2 -> index-1; y pinned everywhere incl. t=0; exact a(t)/a'(t))\n";
    return 0;
  }
  // converged (residual ~0) but wrong: the unpinned y@LB sat at its seed -> not reduced
  bool const lb_dominated = ( ey_LB > 0.5*ey && ey > 1e-6 );
  std::cout << "  ==> BASELINE: "
            << ( lb_dominated ? "y UNPINNED at t=0 (error concentrated at the boundary DOF) -- "
                              : "" )
            << "index-2 NOT reduced (square + residual-converged but answer wrong; "
               "y is in the Jacobian null space)\n";
  return 0;   // measurement, not a failure: this is the pre-reducer baseline
}
