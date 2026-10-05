// OCFE_PDE6_blk0.cpp
// ===========================================================================
//  PDE6 BLOCK 0 -- FULL FACTORIAL EQUIVALENCE HARNESS
//    formulation { MANUAL, AUTO } x solve { mono, march }
//  x imposition  { IC_WEAK, IC_TRACE, IC_STRONG } x backend { SuperLU, SPQR }
//                                                              = 24 cells
// ===========================================================================
//
//  PURPOSE: this driver exists to be FIXED AGAINST.  MANUAL and AUTO are two ways
//  of writing the SAME parabolic problem; the target is that they agree in every
//  cell.  Today they do not, and the harness localises exactly where.
//
//      MANUAL:  T_t - q_x = 0,  q - a*T_x = 0     q declared via add_state
//      AUTO  :  T_t - (a*T_x)_x = 0               auxiliary built by reduce_order
//
//  Same physics, mesh, IC, BCs and perturbed start.
//
//  MEASURED (rev123, NEL_T=NEL_X=5, IC_STRONG, mono)
//      MANUAL SuperLU  |r|=4.4253e-13  err=4.8179e-04
//      MANUAL SPQR     |r|=4.3343e-13  err=2.9510e+00    ** 6125x apart **
//      AUTO   SuperLU  |r|=3.8414e-13  err=1.3055e-05
//      AUTO   SPQR     |r|=3.8414e-13  err=1.3055e-05    identical
//
//  Both converge to ~4e-13.  MANUAL's ANSWER DEPENDS ON THE LINEAR SOLVER.  The
//  cause is a determinacy defect, and it is NOT the k family:
//      PLAN METRIC: rank(B) = rank(C) = n_trace_var,  k = 0
//      Js 1250x1250 cond=3.60e+16 rank_deficiency=4  (AUTO: 1.60e+02, 0)
//      right-null ~100% in the STATES, 0.1% in tau
//      *** CONSISTENT but PRIMAL UNDETERMINED ***
//  rank_deficiency = NEL_X-1 (2,3,4,5 at NEL_X 3,4,5,6), CONSTANT in NEL_T.
//
//  THE BACKEND PAIR IS THE DETERMINACY TEST.  A non-unique primal shows as two
//  solvers reaching the SAME RESIDUAL and DIFFERENT ANSWERS -- cheap, no SVD, no
//  audit flag.  Disagreement is PROOF; agreement is evidence.
//
//  THE CENSUS REFUTED THE OBVIOUS EXPLANATION.  reduce_order's auxiliary is claimed
//  MORE, not less, and the difference is in the EVOLUTION direction:
//                      MANUAL      AUTO
//      t / T               92        92
//      t / aux     -- none --       100      <-- the whole difference
//      x / T               84        80
//      x / aux            100        84
//      total              276       356
//  q is algebraic in t, so no t-claim is minted on it; Daux3_T gets 100.  MORE
//  constraints, not fewer, is what makes AUTO determinate.
//
//  Also: the redundancy detector FIRES on MANUAL (rebuild at phase
//  'detect_redundant_continuity_claims'), drops 4 claims 276->272, and leaves a
//  rank deficiency of exactly 4.  Whether they are the SAME 4 is NOT established
//  and is worth knowing before a fix is written.
//
//  USING IT WHILE FIXING OCFESLV
//    1. every EQUIV row should read EQUIV; every DETERM row should read DETERM
//    2. exit 0 iff both, in all cells -- that is the target
//    3. IC_WEAK is the WITHIN-RUN CONTROL: it imposes no exact continuity, so a
//       failure there is the model or the harness, not the interface plan.  Do not
//       chase IC_TRACE/IC_STRONG until IC_WEAK is clean.
//
//  For the numbers behind a failing cell (dense SVD per cell -- narrow with --only):
//    CRONOS_AUDIT_SPECTRUM=1 ./OCFE_PDE6_blk0 --only strong 2>&1 | grep -E "rank_def|PRIMAL"
//    CRONOS_AUDIT_EVOFACE=1  ./OCFE_PDE6_blk0 --only strong 2>&1 | grep -E "claim census|evoface\]     "
//
//  USAGE
//    ./OCFE_PDE6_blk0 [--nelx N] [--nelt N] [--tol R] [--only weak|trace|strong] [-v]
//      --legend   print each cell's equation list -- what reduce_order BUILT, with roles
//      --drop     post-solve over-drop gate -- NAMES the claims the detector removed
//      --linklast register MANUAL's LINK row LAST, as reduce_order does.  The domain
//                 specs are identical between the forms (rev125), so registration
//                 ORDER is the one structural difference left standing.
//
//    the two open questions, and the flag for each:
//      "which 4 claims were dropped?"          -> --drop --only strong
//      "why 100 claims apart if identical?"    -> --legend --only weak
//
// ===========================================================================

#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <limits>

#include OCFE_OCFESLV_HEADER

using namespace mc;

static constexpr double kPi = 3.14159265358979323846264338327950288;

struct Par
{
  double a  = 0.1;    // thermal diffusivity
  double T0 = 10.0;   // initial amplitude
  double Ts = 5.0;    // surface temperature
  double xf = 1.0;
  double tf = 1.0;
};

// Analytical solution of the heat equation with T(t,0)=Ts, T_x(t,xf)=0.
static double T_exact( double t, double x, Par const& p )
{
  double const k = kPi/2.0/p.xf;
  return p.Ts + p.T0 * std::exp( -p.a*k*k*t ) * std::sin( k*x );
}
static double q_exact( double t, double x, Par const& p )   // q = a * T_x
{
  double const k = kPi/2.0/p.xf;
  return p.a * p.T0 * k * std::exp( -p.a*k*k*t ) * std::cos( k*x );
}

static size_t NELT = 5, NELX = 5, NT = 5, NX = 5;
static int    g_disp = 0;
static int    g_only = -1;          // -1 = all three impositions
static bool   g_legend = false;     // --legend: print each cell's equation list
static bool   g_drop   = false;     // --drop  : post-solve over-drop gate, names the claims
static bool   g_linklast = false;   // --linklast: register MANUAL's LINK row LAST, as reduce_order does
static bool   g_bcuderiv = false;   // --bcuderiv: write MANUAL's upper BC as a DERIVATIVE, as AUTO must

// Two cells agree when their max-errors match to this RELATIVE tolerance.
//
// NOT round-off.  MANUAL and AUTO are genuinely DIFFERENT discretisations -- they
// carry different state sets ({T,q} vs {T,Daux3_T}) -- so equality at 1e-15 is the
// wrong expectation.  Measured when both are determinate: MANUAL/march 1.3057e-05,
// AUTO/march 1.3057e-05, AUTO/mono 1.3055e-05 -- four digits.  The defect targeted
// here is 37x (and 6125x across backends).  1e-3 separates those by three orders
// in both directions.
static double AGREE_REL = 1.0e-3;

static char const* IMPNAME[3] = { "IC_WEAK", "IC_TRACE", "IC_STRONG" };
static OCFESLV::Options::ImpositionType imp_of( int i )
{
  return i == 0 ? OCFESLV::Options::IC_WEAK
       : i == 1 ? OCFESLV::Options::IC_TRACE
                : OCFESLV::Options::IC_STRONG;
}

struct Run
{
  bool   ran = false, setup_ok = false, square = false, conv = false, threw = false;
  size_t nVar = 0, nEqn = 0, nTrace = 0;
  double eT = std::numeric_limits<double>::infinity();
  double resid = std::numeric_limits<double>::infinity();
};

static Run run_cell( bool manual, bool march, int impi, bool spqr )
{
  Run R;
  Par p;

  // The backend is a run-time option in rev113+, so both columns come from ONE
  // binary and one build.  A rebuild between them would reintroduce exactly the
  // stale-object ambiguity this corpus has been bitten by.
#if !defined(CRONOS__WITH_SPQR)
  if( spqr ) return R;   // ran == false -> reported n/a, NEVER folded into a pass
#endif
  R.ran = true;

  FFGraph DAG;
  FFVar t = DAG.add_var( "t" );
  FFVar x = DAG.add_var( "x" );
  FFVar T = DAG.add_var( "T(t,x)" );
  FFVar q = DAG.add_var( "q(t,x)" );

  FFPartial OpP;
  double const av = p.a, T0v = p.T0, Tsv = p.Ts, xf = p.xf;

  OCFESLV oc( &DAG );
  oc.options.DISPLAY_LEVEL = g_disp;
  oc.add_domain( t, FFDom( 0., p.tf, NELT, FFDom::LGR, NT ) );
  oc.add_domain( x, FFDom( 0., p.xf, NELX, FFDom::LGL, NX ) );

  oc.add_state( T, { t, x } );
  oc.update_ref( T, [&]( OCFESLV::t_Coord const& cd ){ return T_exact( cd.at(t), cd.at(x), p ); } );
  if( manual ){
    oc.add_state( q, { t, x } );
    oc.update_ref( q, [&]( OCFESLV::t_Coord const& cd ){ return q_exact( cd.at(t), cd.at(x), p ); } );
  }
  oc.set_evolution_domain( t );

  FFVar HEAT_PDE  = manual ? ( OpP( T, t ) - OpP( q, x ) )
                           : ( OpP( T, t ) - OpP( av*OpP( T, x ), x ) );
  FFVar HEAT_LINK = q - av*OpP( T, x );
  FFVar HEAT_INI  = T - Tsv - T0v*sin( kPi/2.0/xf * x );
  FFVar HEAT_BCL  = T - Tsv;
  // Upper BC: zero flux -- AND A CONFOUND, now switchable.
  //
  // AUTO has no q, so its BC must go on the derivative.  MANUAL has q and by default
  // uses it.  Mathematically the same condition; STRUCTURALLY NOT:
  //     MANUAL default:  q       -> a VALUE condition on a state
  //     AUTO          :  a*T_x   -> a DERIVATIVE condition on T
  //
  // The reduction keeps those in SEPARATE per-(block,direction,face) maps:
  //     value_prims  -> root prims carrying a VALUE condition
  //     deriv_auxes  -> auxes carrying a DERIVATIVE condition
  //     derivbc_eqn  -> aux id to the wEqn index of its derivative BC
  // so it is not cosmetic; it lands in bookkeeping the reduction consults.
  //
  // The comment here previously read "this is the ONLY place the two formulations
  // differ outside HEAT_PDE", which was wrong twice: it is a SECOND difference, and it
  // is not incidental.
  //
  // --bcuderiv writes MANUAL's BC as a derivative too, leaving the PDE form as the
  // SOLE difference.
  FFVar HEAT_BCU  = ( manual && !g_bcuderiv ) ? q : ( av*OpP( T, x ) );

  // --linklast: the ONE structural difference still standing between the two forms.
  //
  // rev125 printed the domain specs and they are IDENTICAL -- every row, both forms:
  //     e INTERIOR |t] x |x|   e LINK [t] x [x]   e INITIAL |t x [x]
  //     e BOUNDARY |t] x |x    e BOUNDARY |t] x x|
  // Same five roles, same five bound masks, same states, same partials.  Only the
  // POSITION of LINK differs: MANUAL registers it second, reduce_order appends it
  // last.  If moving it makes MANUAL determinate, the 100-claim gap is registration
  // order and nothing else.
  //
  // State order is NOT a confound: MANUAL declares T then q, and reduce_order creates
  // Daux3_T after T, so the auxiliary is second in both.
  auto add_link = [&](){
    oc.add_equation( HEAT_LINK, { t, x }, { FFDom::ALL, FFDom::ALL },
                     OCFESLV::EqnOptions( OCFESLV::EqnRole::LINK, 0 ) );
  };

  oc.add_equation( HEAT_PDE, { t, x }, { FFDom::ALL-FFDom::LB, FFDom::ALL-FFDom::LB-FFDom::UB },
                   OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  if( manual && !g_linklast ) add_link();
  oc.add_equation( HEAT_INI, { t, x }, { FFDom::LB, FFDom::ALL },
                   OCFESLV::EqnOptions( OCFESLV::EqnRole::INITIAL, 0 ) );
  oc.add_equation( HEAT_BCL, { t, x }, { FFDom::ALL-FFDom::LB, FFDom::LB },
                   OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );
  oc.add_equation( HEAT_BCU, { t, x }, { FFDom::ALL-FFDom::LB, FFDom::UB },
                   OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );
  if( manual && g_linklast ) add_link();      // AUTO's order: PDE, INI, BCL, BCU, LINK

  oc.options.INTERFACE.TYPE  = OCFESLV::Options::IC_AUTO;
  oc.options.INTERFACE.IMPOSITION = imp_of( impi );
  oc.options.INTERFACE.SAT_SIGMA0      = 10.0;
  oc.options.CLASSIFY.MODE        = OCFESLV::Options::CLASS_AUTO;
  // MONO vs MARCH is the third axis.  The free direction lives at an INTERIOR
  // t-INTERFACE (the right-null dominant columns are all at t=0.2, the first one,
  // across every x-element), and marching solves one evolution window at a time --
  // so there is no interior t-interface inside the system being factorised.
  //
  // PREDICTION: marching removes the MANUAL deficiency, and MANUAL/march AGREES
  // across backends while MANUAL/mono does not.  If marching does NOT fix it, the
  // free direction is not the t-interface after all and the reading above is wrong.
  oc.options.SOLVE.MARCHING  = march;
  // THE BACKEND IS SET ON THE OPTION, NOT THROUGH THE ENVIRONMENT.
  //
  // CRONOS_SOLVE_FACTORIZATION arrived in rev113.  Run this harness against a header
  // that predates it -- ocfeslv_causalflux9, say -- and setenv() is silently ignored,
  // BOTH columns run SuperLU, and the report proudly announces 12/12 DETERM.  That is
  // exactly what happened: every causalflux9 cell reproduced rev125's SuperLU value,
  // and the RESIDUALS were bit-identical across the two "backends" (6.2769e-11 in
  // both), which two different factorisations never are -- under rev125 the same cell
  // gives 6.2769e-11 and 2.1775e-11.
  //
  // Assigning the option directly cannot fail quietly: a header without the field does
  // not COMPILE, which is the loud failure this needs.  The env var is also cleared so
  // a stale one in the shell cannot override the assignment.
  // SOLVE_SPQR only EXISTS under CRONOS__WITH_SPQR, so the reference needs the same
  // guard -- and that is a feature: a header lacking the field, or lacking the enum
  // value, fails to COMPILE.  Silent fallback is impossible in either direction.
  unsetenv( "CRONOS_SOLVE_FACTORIZATION" );
#if defined(CRONOS__WITH_SPQR)
  oc.options.SOLVE.FACTORIZATION = spqr ? OCFESLV::Options::SOLVE_SPQR
                                        : OCFESLV::Options::SOLVE_SUPERLU;
#else
  oc.options.SOLVE.FACTORIZATION = OCFESLV::Options::SOLVE_SUPERLU;   // spqr already returned
#endif
  oc.options.SOLVE.RES_TOL   = 1.0e-9;
  // MANUAL supplies its own auxiliary, so order reduction has nothing to do.
  // AUTO relies on it -- that is the whole comparison.
  // RED_FULL, not RED_MAIN: the point is to let the framework build the auxiliary
  // for the second-order term, which is what removed the deficiency on full PDE6.
  // If RED_MAIN gives a different rank, that difference is itself the finding and
  // this line is the knob to turn.
  oc.options.REDUCE.ORDER    = manual ? OCFESLV::Options::RED_NONE
                                      : OCFESLV::Options::RED_FULL;

  try { R.setup_ok = oc.setup(); }
  catch( ... ) { R.threw = true; return R; }
  if( !R.setup_ok ) return R;

  // --legend answers "why 100 claims apart?".  print_equation_legend() is public and
  // needs no audit flag; until now it only appeared when tier-0 went singular, which
  // is precisely the case AUTO does NOT hit -- so the AUTO equation list, the one
  // that would say what reduce_order actually built and with what ROLE, was never
  // visible.  Measured census, and the roles are SWAPPED rather than merely different:
  //     MANUAL  q    : FULL 100 in x,   NONE in t
  //     AUTO    aux  : FULL 100 in t,   84 in x
  // q = a*T_x has no time derivative, so MANUAL minting no t-claim on it is correct.
  // AUTO's auxiliary is meant to BE that quantity and gets the full 100 t-claims, as
  // if it were a differentiated state.  T moves too: x/T goes 84 -> 80.
  if( g_legend ){
    std::cerr << "\n===== equation legend: " << ( manual ? "MANUAL" : "AUTO" )
              << "  " << IMPNAME[impi] << ( march ? "  march" : "  mono" )
              << ( spqr ? "  SPQR" : "  SuperLU" ) << " =====" << std::endl;
    try { oc.print_equation_legend( std::cerr ); } catch( ... ) {}
  }

  try {
    R.nVar   = oc.n_colloc_sta();
    R.nEqn   = oc.n_colloc_eqn();
    R.nTrace = oc.n_colloc_trace();
  } catch( ... ) {}
  R.square = ( R.nVar && R.nVar == R.nEqn );

  std::vector<double> varInit, inpInit;
  if( !oc.init( varInit, inpInit, nullptr ) ) return R;
  std::vector<double> xv = varInit;
  // Perturb off the exact profile, as PDE6 does (TEST_HA_INIT_MODE 1).  Starting
  // ON the solution would hide a determinacy defect entirely: with nothing to
  // solve for, every solver returns the point it was handed.
  for( size_t i = 0; i < xv.size(); ++i )
    xv[i] += 1.0e-3*std::sin( 0.37*double(i+1) )*std::max( 1.0, std::fabs( xv[i] ) );

  OCFESLV::SolveReport const rep = oc.solve( xv.data(), nullptr, nullptr );
  R.conv  = rep.converged;
  R.resid = rep.final_residual;

  // --drop answers "which 4 claims were dropped, and were they safe to drop?".
  // The redundancy detector FIRES on MANUAL (plan rebuild at phase
  // 'detect_redundant_continuity_claims'), removes 4 claims (census 276 -> tau 272),
  // and leaves a rank deficiency of exactly 4.  Whether those are the SAME 4 is the
  // open question, and this gate is what can answer it: it re-evaluates each DROPPED
  // claim's continuity on the converged solution and returns the keys of any whose
  // continuity does NOT in fact hold.
  //
  // rel_tol must sit ABOVE the truncation floor and BELOW a genuine jump.  1e-2 is
  // the documented default; the max|T-T*| spread here runs 1.3e-05 to 6.2e+01, so a
  // real violation is far above it and truncation far below.
  if( g_drop && R.conv ){
    // The overdropped-set overload takes t_FrozenInterfaceClaimKey, which is
    // PROTECTED -- a driver cannot name the claims.  So pass nullptr and rely on the
    // one-line summary verify_interface_drop always prints (worst relative violation,
    // claims checked).  Naming them belongs on the HEADER side, where the type is
    // visible: CRONOS_AUDIT_DROPSET=1 in rev124 does that.
    try {
      auto const st = oc.verify_interface_drop( xv.data(), 1.0e-2, nullptr );
      std::cerr << "  [drop] " << ( manual ? "MANUAL" : "AUTO" ) << "  " << IMPNAME[impi]
                << ( march ? "  march" : "  mono" ) << ( spqr ? "  SPQR" : "  SuperLU" )
                << "  status="
                << ( st == OCFESLV::InterfaceDropStatus::OK       ? "OK"
                   : st == OCFESLV::InterfaceDropStatus::OVERDROP ? "OVERDROP"
                                                                : "UNDERDROP" )
                << "   (counts and worst violation are in the line above, printed by"
                   " verify_interface_drop itself; CRONOS_AUDIT_DROPSET=1 names them)"
                << std::endl;
    }
    catch( ... ) { std::cerr << "  [drop] verify_interface_drop threw" << std::endl; }
  }

  // Error on T only: it is the state BOTH formulations carry, so it is the one
  // number that is comparable across them.  q exists only in the manual form.
  {
    double emax = 0.;
    size_t const NS = 21;
    for( size_t i = 0; i <= NS; ++i )
      for( size_t j = 0; j <= NS; ++j ){
        double const tt = p.tf*double(i)/double(NS), xx = p.xf*double(j)/double(NS);
        OCFESLV::t_Coord pt; pt[t] = tt; pt[x] = xx;
        double v = 0.;
        try { v = oc.eval_colloc<double>( T, pt, xv.data(), nullptr, nullptr ); }
        catch( ... ) { continue; }
        emax = std::max( emax, std::fabs( v - T_exact( tt, xx, p ) ) );
      }
    R.eT = emax;
  }
  return R;
}


// ---------------------------------------------------------------------------
static bool ok( Run const& R )
{ return R.ran && R.setup_ok && R.square && R.conv && std::isfinite( R.eT ); }

// 0 = agree, 1 = disagree, -1 = NOT COMPARABLE.
// Not-comparable is never folded into agree: a skipped SPQR build or a failed
// setup must not read as a pass.
static int cmp( Run const& A, Run const& B )
{
  if( !ok( A ) || !ok( B ) ) return -1;
  double const s = std::max( std::fabs( A.eT ), std::fabs( B.eT ) );
  if( s <= 0. ) return 0;
  return ( std::fabs( A.eT - B.eT ) / s < AGREE_REL ) ? 0 : 1;
}
static char const* verd( int c, char const* yes, char const* no )
{ return c < 0 ? "n/a" : ( c == 0 ? yes : no ); }

static std::string stat_of( Run const& R )
{
  if( !R.ran )      return "n/a (no SPQR build)";
  if( R.threw )     return "THREW";
  if( !R.setup_ok ) return "setup failed";
  if( !R.square )   return "NOT SQUARE";
  if( !R.conv )     return "not converged";
  return "solved";
}
static std::string num( double v, bool valid )
{
  if( !valid ) return "      --     ";
  std::ostringstream o; o << std::scientific << std::setprecision(4) << v;
  return o.str();
}
static std::string ratio_of( Run const& A, Run const& B, int c )
{
  if( c < 0 || A.eT <= 0. || B.eT <= 0. ) return "-";
  std::ostringstream o; o << std::fixed << std::setprecision(2)
                          << ( std::max(A.eT,B.eT)/std::min(A.eT,B.eT) ) << "x";
  return o.str();
}

int main( int argc, char* argv[] )
{
  for( int i = 1; i < argc; ++i ){
    if( !std::strcmp( argv[i], "-v" ) ){ g_disp = 1; continue; }
    if( !std::strcmp( argv[i], "--legend" ) ){ g_legend = true; continue; }
    if( !std::strcmp( argv[i], "--drop"   ) ){ g_drop   = true; continue; }
    if( !std::strcmp( argv[i], "--linklast" ) ){ g_linklast = true; continue; }
    if( !std::strcmp( argv[i], "--bcuderiv" ) ){ g_bcuderiv = true; continue; }
    if( !std::strcmp( argv[i], "--nelx" ) && i+1 < argc ){ NELX = (size_t)std::atoi( argv[++i] ); continue; }
    if( !std::strcmp( argv[i], "--nelt" ) && i+1 < argc ){ NELT = (size_t)std::atoi( argv[++i] ); continue; }
    if( !std::strcmp( argv[i], "--tol"  ) && i+1 < argc ){ AGREE_REL = std::atof( argv[++i] ); continue; }
    if( !std::strcmp( argv[i], "--only" ) && i+1 < argc ){
      std::string const w = argv[++i];
      g_only = ( w=="weak" ? 0 : w=="trace" ? 1 : w=="strong" ? 2 : -1 );
      if( g_only < 0 ){
        std::cerr << "OCFE_PDE6_blk0 ** --only '" << w
                  << "' is not weak|trace|strong." << std::endl;
        return 2;
      }
      continue;
    }
    // REJECT anything unrecognised.  This loop used to fall through silently, and an
    // ignored --bcuderiv produced a "comparison" of default against default that looked
    // like a clean negative result.  A switch that does nothing must say so.
    std::cerr << "OCFE_PDE6_blk0 ** unrecognised argument '" << argv[i] << "'.\n"
                 "  known: --nelx N --nelt N --tol R --only weak|trace|strong\n"
                 "         --legend --drop --linklast --bcuderiv -v" << std::endl;
    return 2;
  }

  static Run R[2][2][3][2];          // [manual][march][imposition][spqr]
  for( int m = 0; m < 2; ++m )
   for( int r = 0; r < 2; ++r )
    for( int i = 0; i < 3; ++i ){
      if( g_only >= 0 && i != g_only ) continue;
      for( int b = 0; b < 2; ++b )
        R[m][r][i][b] = run_cell( m == 0, r == 1, i, b == 1 );
    }

  std::cout
    << "\n================================================================================\n"
    << "  PDE6 BLOCK 0 -- equivalence harness   MANUAL vs AUTO\n"
    << "    MANUAL:  T_t - q_x = 0,  q - a*T_x = 0    (q via add_state)\n"
    << "    AUTO  :  T_t - (a*T_x)_x = 0              (auxiliary via reduce_order)\n"
    << "  NEL_T=" << NELT << " NEL_X=" << NELX << "  agree_rel="
    << std::scientific << std::setprecision(1) << AGREE_REL << "  solve_tol=1e-09\n"
    << "  TARGET: every DETERM row DETERM, every EQUIV row EQUIV.  Exit 0 iff both.\n"
    << "  header: " << OCFESLV::HEADER_ID << "\n"
    << "  MANUAL LINK registered " << ( g_linklast ? "LAST (--linklast, matching reduce_order)"
                                                   : "SECOND (default)" ) << "\n"
    << "================================================================================\n\n"
    << "  ---- ALL CELLS ---------------------------------------------------------------\n"
    << "  imposition  solve  backend   form     nVar nTrace      final|r|      max|T-T*|  status\n";

  for( int i = 0; i < 3; ++i ){
    if( g_only >= 0 && i != g_only ) continue;
    for( int r = 0; r < 2; ++r )
      for( int b = 0; b < 2; ++b )
        for( int m = 0; m < 2; ++m ){
          Run const& C = R[m][r][i][b];
          std::cout << "  " << std::left << std::setw(12) << IMPNAME[i]
                    << std::setw(7)  << ( r ? "march" : "mono" )
                    << std::setw(10) << ( b ? "SPQR" : "SuperLU" )
                    << std::setw(9)  << ( m ? "AUTO" : "MANUAL" )
                    << std::right << std::setw(5) << C.nVar << std::setw(7) << C.nTrace
                    << "  " << num( C.resid, ok(C) ) << "  " << num( C.eT, ok(C) )
                    << "  " << stat_of( C ) << "\n";
        }
    std::cout << "\n";
  }

  std::cout
    << "  ---- DETERMINACY  (SuperLU vs SPQR, SAME formulation) ------------------------\n"
    << "  disagreement is PROOF the primal is not unique: same equations, same\n"
    << "  residual, different answer.  agreement is evidence, not proof.\n"
    << "  imposition  solve  form        SuperLU          SPQR       ratio  verdict\n";
  int det_ok=0, det_bad=0, det_na=0;
  for( int i = 0; i < 3; ++i ){
    if( g_only >= 0 && i != g_only ) continue;
    for( int r = 0; r < 2; ++r )
      for( int m = 0; m < 2; ++m ){
        Run const& A = R[m][r][i][0]; Run const& B = R[m][r][i][1];
        int const c = cmp( A, B );
        c < 0 ? ++det_na : ( c ? ++det_bad : ++det_ok );
        std::cout << "  " << std::left << std::setw(12) << IMPNAME[i]
                  << std::setw(7) << ( r ? "march" : "mono" )
                  << std::setw(9) << ( m ? "AUTO" : "MANUAL" )
                  << num( A.eT, ok(A) ) << "  " << num( B.eT, ok(B) )
                  << "  " << std::left << std::setw(7) << ratio_of( A, B, c )
                  << verd( c, "DETERM", "** UNDETERMINED **" ) << "\n";
      }
  }

  std::cout
    << "\n  ---- EQUIVALENCE  (MANUAL vs AUTO, SAME cell) --------------------------------\n"
    << "  imposition  solve  backend      MANUAL           AUTO       ratio  verdict\n";
  int eq_ok=0, eq_bad=0, eq_na=0;
  for( int i = 0; i < 3; ++i ){
    if( g_only >= 0 && i != g_only ) continue;
    for( int r = 0; r < 2; ++r )
      for( int b = 0; b < 2; ++b ){
        Run const& A = R[0][r][i][b]; Run const& B = R[1][r][i][b];
        int const c = cmp( A, B );
        c < 0 ? ++eq_na : ( c ? ++eq_bad : ++eq_ok );
        std::cout << "  " << std::left << std::setw(12) << IMPNAME[i]
                  << std::setw(7) << ( r ? "march" : "mono" )
                  << std::setw(10) << ( b ? "SPQR" : "SuperLU" )
                  << num( A.eT, ok(A) ) << "  " << num( B.eT, ok(B) )
                  << "  " << std::left << std::setw(7) << ratio_of( A, B, c )
                  << verd( c, "EQUIV", "** DIFFER **" ) << "\n";
      }
  }

  bool const clean = ( det_bad==0 && eq_bad==0 && det_na==0 && eq_na==0 );
  std::cout
    << "\n================================================================================\n"
    << "  DETERMINACY : " << det_ok << " DETERM, " << det_bad << " UNDETERMINED, " << det_na << " n/a\n"
    << "  EQUIVALENCE : " << eq_ok  << " EQUIV, "  << eq_bad  << " DIFFER, "       << eq_na  << " n/a\n"
    << "  PDE6_blk0: " << ( clean ? "PASS" : "FAIL" ) << "\n"
    << "================================================================================\n"
    << "  READ IN THIS ORDER:\n"
    << "   1. IC_WEAK rows.  No exact continuity is imposed there, so a failure means\n"
    << "      the MODEL or this harness is wrong, not the interface plan.  Nothing\n"
    << "      downstream is interpretable until IC_WEAK is clean.\n"
    << "   2. DETERMINACY.  '** UNDETERMINED **' is proof of a non-unique primal.\n"
    << "   3. EQUIVALENCE.  '** DIFFER **' with BOTH sides determinate means the two\n"
    << "      formulations really are different discretisations.  With one side\n"
    << "      undetermined it only means you are comparing against a coin flip --\n"
    << "      fix determinacy FIRST, then re-read equivalence.\n"
    << "  'n/a' counts as FAIL.  A skipped SPQR build (no -DCRONOS__WITH_SPQR) or a\n"
    << "  failed setup must never look like agreement.\n"
    << "\n"
    << "  A PASS does NOT mean the plan is right -- only that these cells agree at one\n"
    << "  mesh.  The deficiency measured NEL_X-1 and was CONSTANT in NEL_T, so widen\n"
    << "  with --nelx/--nelt before believing a fix.\n"
    << "================================================================================\n";

  // Exit code, unlike the earlier probe version of this driver.  The purpose has
  // changed from "measure what is happening" to "give a fix something to converge
  // against", and a harness meant to be fixed against needs a gate.
  return clean ? 0 : 1;
}
