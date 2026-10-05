// OCFE_PDE20e.cpp  ---  Stage-2 GATE: initial data of a REDUCED high-index model
// ===========================================================================
// REWRITTEN 2026-09-23 for the free-data convention (REDUCE.HIDDEN_IC).
//
// The previous version gated the OLD convention: the modeller declared one IC per differential
// variable and the reducer auto-synthesized a witness-only IC through the re-pivot path.  That
// path is superseded: _reduce_high_index now materialises the HIDDEN constraint levels at the
// evolution boundary, so the model declares only the data it may actually CHOOSE.
//
// Corpus M2 (index 2):  x(t) differential, y(t) algebraic (hidden);  a(t)=1+0.5t-0.3t^2
//   DIFF : x' - y = 0                 (interior)
//   ALG  : x - a  = 0   --reduce-->   y - a'(t) = 0   (pins y at ALL nodes)
//   hidden level 0 (x - a) is materialised at t=0.
// ONE differential state, ONE hidden level  =>  FREE = 0: this model admits NO initial data.
//
// CASES
//   (A) nothing declared          -> the level pins x@LB; balance 0; audit clean; solve exact
//   (B) one IC declared (x@LB)    -> SURPLUS 1: the count now SEES over-specified initial data,
//                                    which the old convention could not (it was "the" correct IC)
//   (C) y pinned twice, x@LB gap, WITH THE SWITCH OFF -> genuine SQUARE rank deficiency; the audit
//                                    still catches it.  The switch is off deliberately: with the
//                                    levels materialised this model is merely non-square, so the
//                                    square-but-singular path -- which is what the audit exists for --
//                                    would no longer be exercised at all.
//   (D) as (C) with the opt-out   -> warns and proceeds, answer wrong
//   (E) marched, nothing declared -> the hidden level is re-imposed at EACH window LB (it is an
//                                    invariant, not data), so the solve stays exact across windows
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

struct CaseResult
{
  bool   setup_ok    = false;
  bool   status_inconsistent = false;   // setup_status()==REDUCED_DOF_INCONSISTENT
  bool   reduced     = false;
  bool   audit_ran   = false;
  bool   audit_square= false;
  size_t audit_rows  = 0, audit_cols = 0, audit_rank = 0, audit_def = 0;
  bool   audit_ok    = false;
  bool   solved      = false;           // a solve was attempted (setup succeeded)
  bool   solve_conv  = false;
  double err_x = 0., err_y = 0.;
  std::string balance; bool balanced = false;
  size_t declared = 0, free = 0, hidden = 0;
};

enum Spec { FREEDATA, OVERSPEC, DEFICIENT };

static CaseResult run_case( Spec spec, bool fatal, bool marching = false, bool hidden = true )
{
  CaseResult R;

  FFGraph DAG;
  FFVar t = DAG.add_var( "t" );
  FFVar x = DAG.add_var( "x(t)" );
  FFVar y = DAG.add_var( "y(t)" );          // ALGEBRAIC, HIDDEN
  FFPartial OpP;

  FFVar Ae    = 1.0 + 0.5*t - 0.3*t*t;      // a(t)
  FFVar Ape   = 0.5 - 0.6*t;                // a'(t)
  FFVar DIFF  = OpP( x, t ) - y;            // x' - y = 0
  FFVar ALG   = x - Ae;                     // 0 = x - a   (no y -> index 2)
  FFVar ICx   = x - Ae;                     // consistency IC for x  (anchors x@LB)
  FFVar ICy   = y - Ape;                    // witness IC: pins y@LB
  FFVar ICy2  = 2.0*y - 2.0*Ape;            // redundant 2nd y@LB pin (same root, distinct row)

  OCFESLV oc( &DAG );
  oc.add_domain( t, FFDom( 0., TF, NEL_T, FFDom::LGR, NT ) );
  oc.add_state( x, {t} );
  oc.add_state( y, {t} );
  oc.update_ref( x, []( OCFESLV::t_Coord const& ){ return 1.0; } );  // adversarial seed
  oc.update_ref( y, []( OCFESLV::t_Coord const& ){ return 0.0; } );

  OCFESLV::EqnOptions io( OCFESLV::EqnRole::INTERIOR, 0 ), ii( OCFESLV::EqnRole::INITIAL, 0 );
  int const T_INT = FFDom::ALL - FFDom::LB;
  oc.add_equation( DIFF, {t}, {T_INT}, io );                  // x' - y = 0   (interior)
  // CORRECT/AUTOFILL: the algebraic constraint spans the evolution-LB, so the
  // reducer can pin x@t=0 from it (auto-synthesizing the IC when only a witness
  // IC is supplied).  DEFICIENT: the constraint is restricted to t>0, which
  // suppresses that synthesis and leaves a genuine x@LB gap.
  oc.add_equation( ALG, {t}, { spec==DEFICIENT ? T_INT : (int)FFDom::ALL }, io );
  switch( spec ){
    case FREEDATA:                                           // nothing: FREE = 0, the hidden
      (void)ICx; (void)ICy; (void)ICy2;                      // level pins x@LB
      break;
    case OVERSPEC:                                           // one IC too many: the level
      oc.add_equation( ICx, {t}, {FFDom::LB}, ii );          // already pins x@LB
      break;
    case DEFICIENT:                                          // y pinned twice, x@LB gap
      oc.add_equation( ICy,  {t}, {FFDom::LB}, ii );
      oc.add_equation( ICy2, {t}, {FFDom::LB}, ii );
      break;
  }
  oc.set_evolution_domain( t );

  oc.options.REDUCE.ORDER     = OCFESLV::Options::TEST_DAE_REDUCE;
  oc.options.CLASSIFY.MODE         = OCFESLV::Options::CLASS_AUTO;
  oc.options.INTERFACE.TYPE   = OCFESLV::Options::IC_AUTO;
  oc.options.INTERFACE.IMPOSITION  = OCFESLV::Options::IC_WEAK;
  oc.options.DISPLAY_LEVEL    = 1;
  oc.options.SOLVE.MAX_ITER   = 50;
  oc.options.SOLVE.RES_TOL    = 1e-9;
  oc.options.FATAL.REDUCED_DOF = fatal;   // <-- the knob under test
  oc.options.REDUCE.HIDDEN_IC  = hidden;  // on: the model declares the free data, levels are rows
  oc.options.SOLVE.MARCHING    = marching;
#if defined(CRONOS__WITH_SPQR)
  oc.options.SOLVE.FACTORIZATION = OCFESLV::Options::SOLVE_SPQR;
#endif

  try { R.setup_ok = oc.setup(); }
  catch( ... ) { R.setup_ok = false; }

  // The audit + reduction plan are populated DURING setup, before any fatal
  // return, so read them unconditionally (valid even when setup() returned false).
  R.status_inconsistent =
      ( oc.setup_status() == OCFESLV::SetupStatus::REDUCED_DOF_INCONSISTENT );
  R.reduced = !oc.reduction_plan().empty();
  OCFESLV::t_DofAudit const& A = oc.reduced_dof_audit();
  R.audit_ran    = A.ran;
  R.audit_square = A.square;
  R.audit_rows   = A.rows;
  R.audit_cols   = A.cols;
  R.audit_rank   = A.rank;
  R.audit_def    = A.deficiency;
  R.audit_ok     = A.ok();
  R.balance      = oc.dof_balance().str;
  R.balanced     = oc.dof_balance().balanced;
  R.declared     = oc.initial_data().declared;
  R.free         = oc.initial_data().free;
  R.hidden       = oc.initial_data().hidden;

  // solve only if setup succeeded (a rejected model has no valid plan to solve)
  if( R.setup_ok ){
    std::vector<double> varInit, inpInit;
    if( oc.init( varInit, inpInit, nullptr ) ){
      R.solved = true;
      std::vector<double> xv = varInit;
      OCFESLV::SolveReport rep = oc.solve( xv.data(), nullptr, nullptr );
      R.solve_conv = rep.converged;
      auto eval = [&]( FFVar const& V, double tt ){
        OCFESLV::t_Coord pt; pt[t]=tt;
        return oc.eval_colloc<double>( V, pt, xv.data(), nullptr, nullptr ); };
      int const NS=101;
      for( int k=0;k<=NS;++k ){ double tt=double(k)/NS*TF;
        R.err_x=std::max(R.err_x,std::fabs(eval(x,tt)-a_exact(tt)));
        R.err_y=std::max(R.err_y,std::fabs(eval(y,tt)-ap_exact(tt))); }
    }
  }
  return R;
}

static void report( char const* tag, CaseResult const& R )
{
  std::cout << "\n  [" << tag << "]\n";
  std::cout << "    setup=" << (R.setup_ok?"OK":"FALSE")
            << "  status_inconsistent=" << (R.status_inconsistent?"yes":"no")
            << "  reduced=" << (R.reduced?"yes":"no") << "\n";
  std::cout << "    initial data: declared=" << R.declared << " hidden=" << R.hidden
            << " free=" << R.free << "   balance: rows-unknowns=" << R.balance
            << ( R.balanced? "  (balanced)": "  (NOT balanced)" ) << "\n";
  std::cout << "    audit: ran=" << (R.audit_ran?"yes":"no")
            << " rows=" << R.audit_rows << " cols=" << R.audit_cols
            << " rank=" << R.audit_rank << " deficiency=" << R.audit_def
            << " square=" << (R.audit_square?"yes":"no")
            << "  ok()=" << (R.audit_ok?"yes":"NO") << "\n";
  if( R.solved )
    std::cout << "    solve: converged=" << (R.solve_conv?"yes":"no")
              << "  max|x-a|=" << std::scientific << std::setprecision(3) << R.err_x
              << "  max|y-a'|=" << R.err_y << "\n";
  else
    std::cout << "    solve: (not attempted -- setup rejected the model)\n";
}

int main()
{
  std::cout << "================================================================\n";
  std::cout << "  Stage-2 GATE: post-index-reduction IC/BC DOF audit (corpus M2)\n";
  std::cout << "  reducer auto-completes a witness-only IC; audit still rejects a\n";
  std::cout << "  genuine square rank-deficiency (opt-out tested)\n";
  std::cout << "================================================================\n";

  std::cout << "\n---- (A) nothing declared: FREE = 0, the hidden level pins x@LB ----\n";
  CaseResult A = run_case( FREEDATA, /*fatal=*/true );
  report( "A free data (none)", A );

  std::cout << "\n---- (B) one IC declared: the count should see a SURPLUS ----------\n";
  CaseResult B = run_case( OVERSPEC, /*fatal=*/true );
  report( "B over-specified", B );

  std::cout << "\n---- (C) RANK-DEFICIENT, switch OFF: square but singular, fatal -\n";
  CaseResult C = run_case( DEFICIENT, /*fatal=*/true, /*marching=*/false, /*hidden=*/false );
  report( "C deficient, fatal", C );

  std::cout << "\n---- (D) RANK-DEFICIENT, switch OFF, fatal=false: warn + proceed -\n";
  CaseResult D = run_case( DEFICIENT, /*fatal=*/false, /*marching=*/false, /*hidden=*/false );
  report( "D deficient, warn", D );

  std::cout << "\n---- (E) MARCHED, nothing declared: the level is an INVARIANT ----\n";
  CaseResult E = run_case( FREEDATA, /*fatal=*/true, /*marching=*/true );
  report( "E free data, marched", E );

  // verdicts
  bool const A_ok = A.setup_ok && !A.status_inconsistent && A.reduced && A.balanced &&
                    A.free==0 && A.declared==0 && A.hidden==1 &&
                    A.audit_ran && A.audit_ok && A.audit_def==0 &&
                    A.solved && A.err_x < 1e-8 && A.err_y < 1e-8;
  // B: the count sees it (+1) AND the existing DOF audit rejects it as non-square -- the check the
  // old convention could not make, since one IC per differential variable WAS the correct declaration.
  bool const B_surplus = ( B.balance == "1" ) && !B.balanced && B.declared==1 && B.free==0
                         && !B.setup_ok && !B.audit_square && !B.solved;
  bool const C_rejected = !C.setup_ok && C.status_inconsistent && C.reduced &&
                          C.audit_ran && C.audit_square && C.audit_def==1 && !C.solved;
  bool const D_warned = D.setup_ok && !D.status_inconsistent && D.reduced &&
                        D.audit_ran && D.audit_def>=1 && D.solved &&
                        !( D.err_x < 1e-8 && D.err_y < 1e-8 );
  bool const E_ok = E.setup_ok && E.balanced && E.solved && E.solve_conv &&
                    E.err_x < 1e-8 && E.err_y < 1e-8;

  std::cout << "\n==================== verdict ====================\n";
  std::cout << "  (A) nothing declared -> level pins x@LB, balance 0, exact         : " << (A_ok?"yes":"NO") << "\n";
  std::cout << "  (B) one IC declared  -> reported SURPLUS of 1                     : " << (B_surplus?"yes":"NO") << "\n";
  std::cout << "  (C) square+singular  -> deficiency 1, setup REJECTED (switch off) : " << (C_rejected?"yes":"NO") << "\n";
  std::cout << "  (D) square+singular  -> opt-out warns + proceeds, answer wrong    : " << (D_warned?"yes":"NO") << "\n";
  std::cout << "  (E) marched          -> invariant re-imposed per window, exact    : " << (E_ok?"yes":"NO") << "\n";
  bool const pass = A_ok && B_surplus && C_rejected && D_warned && E_ok;
  std::cout << "\n  ==> " << ( pass
        ? "INITIAL-DATA:PASS (free-data convention: the levels pin the model, a declared IC is a surplus,"
          " a genuine deficiency is still caught, and marching re-imposes the invariant)"
        : "INITIAL-DATA:FAIL (did not behave as expected -- read the per-case report above)" ) << "\n";
  return pass ? 0 : 1;
}
