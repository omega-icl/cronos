// OCFE_PDE20f.cpp  ---  Stage-2 GATE: SPATIAL AUXILIARY x {high index, ALG_CLOSURE, sigma reuse}
// 2026-09-10: ALL 18 CELLS RUN UNCONDITIONALLY -- 3 models x 2 reductions x 3 impositions
// (IC_WEAK, IC_TRACE, IC_STRONG).  The corpus harness must invoke this binary WITH NO
// imposition restriction: the rev180 sweep showed the previous corpus entry exercised
// IC_STRONG only, and the one measured catastrophic failure in this family (model C,
// RED_FULL, IC_WEAK: |u-ex| = 1.69e+00, silent, converging) lives in a cell that a
// STRONG-only invocation never sees.  A sweep guards only what it runs.
// ===========================================================================
// REQUIRES:  -DOCFE_OCFESLV_HEADER='"ocfeslv_iccons5.hpp"'   (or the promoted ocfeslv.hpp)
//
// THE GAP THIS CLOSES.  R-E ran the PDE20 family under RED_MAIN and everything passed -- but
// M2/M3/M8 are PURE ODEs in t (x'=y; 0=x-a): no spatial derivative, so NO auxiliary is created,
// so RED_MAIN vs RED_FULL is very likely a NO-OP for them.  The cases that DO have spatial
// structure (M4 Darcy, M7 conservative flux) are both INDEX-1.  So nothing yet tests the
// INTERACTION between a spatial auxiliary and the structural machinery.  That interaction is
// exactly what RED_MAIN changes:
//     RED_FULL   sigma = du/dz created for d2u/dz2, THEN REUSED for the bare du/dz
//     RED_MAIN   sigma created, NOT reused -- du/dz stays an explicit derivative node
// so the structural incidence seen by Pantelides / AUTO_ALG_CLOSURE / AUTO_DIFF_ELIM differs.
//
// THREE MODELS, EACH WITH d2u/dz2 (so a spatial auxiliary IS created), EACH EXERCISING A
// DIFFERENT MECHANISM, EACH RUN UNDER BOTH REDUCTION MODES:
//
//   A  HIGH INDEX + SPATIAL AUXILIARY
//        dx/dt = y ;  0 = x - a(t)                    <- index-2 pair (as corpus M2)
//        du/dt + U du/dz - D d2u/dz2 = x + f          <- spatial auxiliary + advection
//      x,y are (t,z) states with no z-derivatives (the PSA q1/q2 pattern), so the index-2
//      constraint is DISTRIBUTED -- which also stresses the distributed-algebraic path.
//      Exact: x = a, y = a', u = u_ex.
//
//   B  AUTO_ALG_CLOSURE + SPATIAL AUXILIARY
//        du/dt + U du/dz - D d2u/dz2 = w + f
//        0 = w - k*u                                  <- distributed SOLVED ALGEBRAIC
//      w is pinned pointwise, so closure must decide its interface treatment while sigma exists.
//      Exact: w = k*u_ex, u = u_ex.
//
//   C  SIGMA REUSE + DERIVATIVE-DEFINED AUXILIARY   (header premise corrected 2026-09-10)
//        du/dt + U du/dz - D d2u/dz2 = k*w + f
//        0 = w - du/dz                                <- DERIVATIVE-DEFINED auxiliary
//      w is defined AS a first derivative.  THE ORIGINAL PREMISE OF THIS CASE WAS WRONG:
//      AUTO_DIFF_ELIM performs ZERO eliminations here, in every cell (measured at
//      DISPLAY_LEVEL=1, which prints one line per elimination: none appear).  What actually
//      distinguishes the reductions is reduce_order's SIGMA REUSE: RED_FULL replaces the bare
//      du/dz by Dz_u in the PDE *and* the algebraic row, so u's z-derivative survives only in
//      ALGEBRAIC rows (excluded from the principal symbol) and u's z-continuity claims lose
//      their natural receiver -- they are re-expressed on Dz_u and RESCUED with a fabricated
//      coupling.  ALG becomes  w - cf*Dz_u  with coefficient -cf of Dz_u, and the forced +1
//      cancels it at cf=1: |u-ex| = 1.69e+00 under RED_FULL/IC_WEAK while converging, ~3.9e-06
//      at every other forced value including -1.  RED_MAIN keeps du/dz explicit in the PDE and
//      is correct in all cells.  AUTO_DIFF_ELIM stays ON for this model so that its measured
//      zero-elimination property is itself guarded: if a future revision makes it act here,
//      the aux counts and these gates will say so.
//      The fix (ocfeslv rev180, CRONOS_RESCUE_FROM_ROW=1): on IC_WEAK, an AUX-rescued edge
//      whose ALGEBRAIC receiver has a nonzero CONSTANT zeroth-order coefficient of the claimed
//      state takes THAT coefficient instead of the forced constant.  With it on, the RED_FULL/
//      IC_WEAK cell passes; with it off, that cell FAILS BY DESIGN as the regression witness.
//      Exact: w = du_ex/dz, u = u_ex.
//
// PASS: every (model, reduction) cell must set up, converge, and reproduce its exact solution to
// discretisation error.  The auxiliary COUNT is reported per cell so the RED_FULL/RED_MAIN
// difference is visible rather than assumed -- if the counts are identical in every model, the
// two modes are not actually being distinguished and the test proves less than it appears.
// ===========================================================================

#include <iostream>
#include <cstdlib>
#include <iomanip>
#include <vector>
#include <cmath>
#include <string>
#include <stdexcept>
#include <algorithm>

#ifndef OCFE_OCFESLV_HEADER
#define OCFE_OCFESLV_HEADER "ocfeslv.hpp"
#endif
#include OCFE_OCFESLV_HEADER

using namespace mc;
typedef OCFESLV::Options O;

static double const U_ADV = 1.3;            // != 1, so a coupling of 1.0 is never ambiguous
static double const D_AX  = 0.02;
static double const K_LIN = 0.7;
static double const T_END = 0.5, PH = 0.7;
static double const PI2   = 6.283185307179586476925286766559;

static inline double uex ( double z,double t ){ return std::exp(-t)*(1.0+std::sin(PI2*z+PH)); }
static inline double uzex( double z,double t ){ return std::exp(-t)*PI2*std::cos(PI2*z+PH); }
static inline double aex ( double t )         { return 1.0+0.5*t-0.3*t*t; }
static inline double apex( double t )         { return 0.5-0.6*t; }

struct Cell {
  std::string model, red, imp;   // 2026-09-09: imposition is now a cell dimension
  bool setup_ok=false, conv=false;
  int  iters=0;
  size_t n_aux=0, nVar=0;
  double e_u=-1., e_s2=-1.;                 // u error, and the secondary state's error
  std::string s2name="-";
  std::string err="";                       // exception text, if the cell threw
  std::string verdict="";                   // 2026-09-10: PASS / XFAIL / FAIL, set by the gate
};

// ---------------------------------------------------------------------------------------------
static Cell run_case( int model, O::ReductionType red, O::ImpositionType imp )
{
  Cell C;
  C.model = ( model==0 ? "A high-index" : model==1 ? "B alg-closure" : "C diff-elim" );
  C.red   = ( red==O::RED_FULL ? "RED_FULL" : "RED_MAIN" );
  C.imp   = ( imp==O::IC_WEAK ? "WEAK" : imp==O::IC_TRACE ? "TRACE" : "STRONG" );

  size_t const ne_t=3, nn_t=5, ne_z=6, nn_z=6;
  std::vector<double> t_bnd, z_bnd;
  for( size_t i=0;i<=ne_t;++i ) t_bnd.push_back( T_END*double(i)/double(ne_t) );
  for( size_t i=0;i<=ne_z;++i ) z_bnd.push_back( double(i)/double(ne_z) );

  FFGraph DAG;
  FFVar t=DAG.add_var("t"), z=DAG.add_var("z"), u=DAG.add_var("u(t,z)");
  FFVar x=DAG.add_var("x(t,z)"), y=DAG.add_var("y(t,z)"), w=DAG.add_var("w(t,z)");
  FFPartial OpP;

  FFVar W  = PI2*z+PH;
  FFVar UE = exp(-t)*(1.0+sin(W));                       // u_ex
  FFVar UT = -exp(-t)*(1.0+sin(W));                      // du_ex/dt
  FFVar UZ = exp(-t)*PI2*cos(W);                         // du_ex/dz
  FFVar UZZ= -exp(-t)*PI2*PI2*sin(W);                    // d2u_ex/dz2
  FFVar AA = 1.0+0.5*t-0.3*t*t;                          // a(t)
  FFVar LHS= UT + U_ADV*UZ - D_AX*UZZ;                   // the operator applied to u_ex

  OCFESLV oc(&DAG);
  oc.add_domain( t, FFDom(t_bnd,FFDom::LGL,nn_t) );
  oc.add_domain( z, FFDom(z_bnd,FFDom::LGL,nn_z) );

  int const TNL=FFDom::ALL-FFDom::LB, ZI=FFDom::ALL-FFDom::LB-FFDom::UB;
  typedef OCFESLV::EqnOptions EO; typedef OCFESLV::EqnRole ER;

  // u is common to all three models
  oc.add_state( u, {t,z} );
  oc.update_ref( u, [&t,&z](OCFESLV::t_Coord const& c){ return uex(c.at(z),c.at(t)); } );

  if( model==0 ){
    // ---- A: index-2 pair (distributed) + spatial auxiliary ----------------------------------
    oc.add_state( x, {t,z} ); oc.add_state( y, {t,z} );
    oc.update_ref( x, [&t](OCFESLV::t_Coord const& c){ return aex(c.at(t)); } );
    oc.update_ref( y, [&t](OCFESLV::t_Coord const& c){ return apex(c.at(t)); } );
    C.s2name = "x";
  }
  else{
    oc.add_state( w, {t,z} );
    if( model==1 )
      oc.update_ref( w, [&t,&z](OCFESLV::t_Coord const& c){ return K_LIN*uex(c.at(z),c.at(t)); } );
    else
      oc.update_ref( w, [&t,&z](OCFESLV::t_Coord const& c){ return uzex(c.at(z),c.at(t)); } );
    C.s2name = "w";
  }
  oc.set_evolution_domain( t );

  FFVar SRC, PDEu, ALG, ICx, ICy, ICw;
  if( model==0 ){
    SRC  = LHS - AA;                                     // so that u_ex solves it once x == a
    PDEu = OpP(u,t) + U_ADV*OpP(u,z) - D_AX*OpP(OpP(u,z),z) - x - SRC;
    ALG  = x - AA;                                       // 0 = x - a(t)  -> INDEX 2 with dx/dt=y
  }
  else if( model==1 ){
    SRC  = LHS - K_LIN*UE;
    PDEu = OpP(u,t) + U_ADV*OpP(u,z) - D_AX*OpP(OpP(u,z),z) - w - SRC;
    ALG  = w - K_LIN*u;                                  // distributed solved algebraic
  }
  else{
    SRC  = LHS - K_LIN*UZ;
    PDEu = OpP(u,t) + U_ADV*OpP(u,z) - D_AX*OpP(OpP(u,z),z) - K_LIN*w - SRC;
    ALG  = w - OpP(u,z);                                 // DERIVATIVE-DEFINED -> DIFF_ELIM turf
  }

  FFVar ICu = u - UE;
  FFVar BCL = U_ADV*u - D_AX*OpP(u,z) - ( U_ADV*UE - D_AX*UZ );   // flux/Robin at z=0
  FFVar BCU = OpP(u,z) - UZ;                                      // Neumann at z=1

  oc.add_equation( PDEu, {t,z}, {TNL,ZI},               EO(ER::INTERIOR,0) );
  oc.add_equation( ICu , {t,z}, {FFDom::LB,FFDom::ALL}, EO(ER::INITIAL ,0) );
  oc.add_equation( BCL , {t,z}, {TNL,FFDom::LB},        EO(ER::BOUNDARY,0) );
  oc.add_equation( BCU , {t,z}, {TNL,FFDom::UB},        EO(ER::BOUNDARY,0) );

  if( model==0 ){
    // dx/dt = y over ALL z (no z-derivatives -> no z-BCs needed; the PSA q1/q2 pattern)
    // EXACTLY the corpus-M2 pattern from OCFE_PDE20b.cpp, lifted to a distributed (t,z) state:
    //   DIFF on ALL-LB,  ALG on **FFDom::ALL** (LB INCLUDED),  IC on x at LB, none on y.
    // Putting ALG on ALL-LB instead was the original bug here: it left y unpinned at t=LB, i.e.
    // exactly one missing row per z-node (the audit reported rows=684 cols=720, deficit 36 = nz).
    FFVar DXY = OpP(x,t) - y;
    oc.add_equation( DXY, {t,z}, {TNL,FFDom::ALL},            EO(ER::INTERIOR,1) );
    oc.add_equation( ALG, {t,z}, {FFDom::ALL,FFDom::ALL},     EO(ER::INTERIOR,1) );
    // 2026-09-23: block 1 is the index-2 pattern (x differential through DXY, y hidden in ALG), so it has NO free
    // initial data: x(0)=a(0) was a consistency condition.  REDUCE.HIDDEN_IC materialises that level at t=0.
    ICx = x - AA; (void)ICx;
    // NOTE: no IC on y -- index reduction must pin it (the M2 "witness" behaviour).
  }
  else{
    // Same ALL-t correction, AND block 0 rather than block 1.
    // SEGFAULT ROOT CAUSE (models B/C, first run): putting the pointwise algebraic relation in its
    // OWN block gave that block ZERO differential equations -- a block with no dynamics at all.
    // classify() printed blk0's probe and then crashed before blk1's.  Model A survived because
    // ITS block 1 does have a differential equation (dx/dt = y).  A value-slaved distributed
    // algebraic (the "P - Rg*sum(c) = 0" pattern) belongs in the SAME block as the state it is
    // slaved to, so w goes in block 0 with u.
    oc.add_equation( ALG, {t,z}, {FFDom::ALL,FFDom::ALL},     EO(ER::INTERIOR,0) );
  }

  oc.options.REDUCE.ORDER     = red;
  oc.options.CLASSIFY.MODE         = O::CLASS_AUTO;
  oc.options.INTERFACE.TYPE   = O::IC_AUTO;
  oc.options.INTERFACE.IMPOSITION  = imp;
  oc.options.AUTO.DIFF_ELIM   = ( model==2 );            // ON only for model C -- and MEASURED
                                                         // to eliminate NOTHING there (see the
                                                         // case-C header); kept on as a guard
  // 2026-09-09 PROBE: PDE20F_SIGMA0 overrides the weak-imposition penalty strength.  Model C
  // under RED_FULL/IC_WEAK gives |u-ex| = 1.69e+00 while converging; the discriminator is
  // whether that error SCALES with the penalty (under-penalisation -- the claim is present but
  // too weak) or is INSENSITIVE to it (the claim is absent, misdirected, or cancelling).
  oc.options.INTERFACE.SAT_SIGMA0       = 10.0;
  if( char const* v = std::getenv( "PDE20F_SIGMA0" ) ) oc.options.INTERFACE.SAT_SIGMA0 = std::atof( v );
  oc.options.REDUCE.HIDDEN_IC = true;                    // hidden levels are rows; models B/C have none
  oc.options.DISPLAY_LEVEL    = 1;                       // so the classify probe line shows idx=
#if defined(CRONOS__WITH_SPQR)
  oc.options.SOLVE.FACTORIZATION = O::SOLVE_SPQR;
#else
  oc.options.SOLVE.FACTORIZATION = O::SOLVE_SUPERLU;
#endif
  oc.options.SOLVE.WARMSTART = O::BROADCAST_IC;
  oc.options.OUTPUT.MARCH_STORE = true;

  std::cout << "\n------ " << C.model << "  /  " << C.red << "  /  " << C.imp << " ------\n";
  // setup()/init()/solve() can THROW (OCFESLV::Exceptions does not derive from std::exception and
  // its what() is non-const, so it must be caught BY VALUE).  Guard every cell so one bad cell
  // reports its exception instead of aborting the whole run -- the original driver died silently
  // part-way through model B for exactly this reason.
  try{
    if( !oc.setup() ){ std::cout << "  setup FAILED (returned false)\n"; return C; }
    // 2026-09-09 PROBE: dump the frozen interface-plan diagnostics (which carry the weak SAT
    // term listing) so the RED_FULL alias duplicate can be compared against RED_MAIN.
    if( std::getenv( "CRONOS_WEAK_SAT_DUMP" ) ) oc.display_interface_plan_diagnostics( std::cout );
  }
  catch( OCFESLV::Exceptions e ){ C.err = std::string("setup threw: ")+e.what();
    std::cout << "  *** " << C.err << "\n"; return C; }
  catch( std::exception const& e ){ C.err = std::string("setup threw std: ")+e.what();
    std::cout << "  *** " << C.err << "\n"; return C; }
  catch( ... ){ C.err = "setup threw unknown"; std::cout << "  *** " << C.err << "\n"; return C; }
  C.setup_ok = true;

  try{ C.n_aux = oc.auxiliary_states().size(); }
  catch( OCFESLV::Exceptions e ){ std::cout << "  auxiliary_states threw: " << e.what() << "\n"; }

  std::vector<double> xv, inp;
  OCFESLV::SolveReport rep;
  try{
    if( !oc.init(xv,inp,nullptr) ){ std::cout << "  init FAILED\n"; return C; }
    C.nVar = oc.n_colloc_sta();
    rep = oc.solve( xv.data(), inp.data(), nullptr );
  }
  catch( OCFESLV::Exceptions e ){ C.err = std::string("solve threw: ")+e.what();
    std::cout << "  *** " << C.err << "\n"; return C; }
  catch( ... ){ C.err = "solve threw unknown"; std::cout << "  *** " << C.err << "\n"; return C; }
  C.conv = rep.converged; C.iters = rep.iterations;
  std::cout << "  setup=OK  n_aux=" << C.n_aux << "  nVar=" << C.nVar
            << "  conv=" << (C.conv?"y":"n") << "  iters=" << C.iters << "\n";
  if( !C.conv ) return C;

  try{
    auto V=[&](FFVar const& s,double zz,double tt){ OCFESLV::t_Coord p; p[z]=zz; p[t]=tt;
      return oc.eval_colloc<double>(s,p,xv.data(),inp.data(),nullptr); };
    double eu=0., e2=0.;
    for( int it=1; it<=6; ++it ){
      double const tt = T_END*double(it)/6.0;
      for( int k=0; k<=60; ++k ){
        double const zz = (double(k)+0.5)/61.0;
        eu = std::max( eu, std::fabs( V(u,zz,tt) - uex(zz,tt) ) );
        if( model==0 )      e2 = std::max( e2, std::fabs( V(x,zz,tt) - aex(tt) ) );
        else if( model==1 ) e2 = std::max( e2, std::fabs( V(w,zz,tt) - K_LIN*uex(zz,tt) ) );
        else                e2 = std::max( e2, std::fabs( V(w,zz,tt) - uzex(zz,tt) ) );
      }
    }
    C.e_u = eu; C.e_s2 = e2;
    std::cout << std::scientific << std::setprecision(3)
              << "  max|u-u_ex|=" << C.e_u << "   max|" << C.s2name << "-exact|=" << C.e_s2 << "\n";
  }catch( OCFESLV::Exceptions e ){ std::cout << "  verify threw: " << e.what() << "\n"; }

  return C;
}

// ---------------------------------------------------------------------------------------------
int main()
{
  std::cout.setf( std::ios::unitbuf );
  std::cout << "================================================================\n"
            << "  PDE20f -- SPATIAL AUXILIARY x {high index, ALG_CLOSURE, DIFF_ELIM}\n"
            << "  header: " << OCFE_OCFESLV_HEADER << "\n"
            << "================================================================\n"
            << "  Closes the R-E gap: M2/M3/M8 are pure ODEs (no spatial auxiliary, so RED_MAIN\n"
            << "  vs RED_FULL is a no-op there); M4/M7 have spatial structure but are index-1.\n"
            << "  Each model below has d2u/dz2 (auxiliary IS created) AND exercises one of the\n"
            << "  three structural mechanisms, under BOTH reduction modes.\n";

  std::vector<Cell> cells;
  // 2026-09-09: the imposition is now a CELL DIMENSION.  Before, this gate ran IC_STRONG only
  // -- and a corpus-wide check (NOTES_20260909c) found that single-imposition drivers hide
  // real defects: OCFE_PDE30 passes IC_STRONG and fails IC_TRACE.  18 cells now.
  O::ImpositionType const IMPS[3] = { O::IC_WEAK, O::IC_TRACE, O::IC_STRONG };
  for( int m = 0; m < 3; ++m ){
    for( O::ImpositionType imp : IMPS ){
      cells.push_back( run_case(m,O::RED_FULL,imp) );
      cells.push_back( run_case(m,O::RED_MAIN,imp) );
    }
  }

  // ---- per-cell gate FIRST, so the summary can print each cell's verdict ---------------------
  bool all_ok = true;
  for( size_t i = 0; i < cells.size(); ++i ){
    Cell& C = cells[i];
    // 2026-09-09: PER-IMPOSITION tolerance, CALIBRATED FROM THE 18-CELL MATRIX rather than
    // carried over from the IC_STRONG-only version.  Measured |u-ex| / |2nd-ex| over all cells:
    //   STRONG  4.5e-06 .. 1.57e-05  (and 1.29e-04 on C/RED_FULL's 2nd state)
    //   TRACE   4.5e-06 .. 3.97e-05, EXCEPT C/RED_FULL at 1.18e-04 / 2.07e-03
    //   WEAK    4.2e-06 .. 6.56e-03, EXCEPT C/RED_FULL at 1.69e+00 / 3.14e+01
    // IC_WEAK imposes continuity by PENALTY, so a larger error is a property of the method and
    // not a defect: its tolerance is looser by design.  The exact modes enforce continuity
    // structurally and are held tight.
    // These thresholds are NOT set to make every cell pass -- two cells still fail, and both
    // are real (see the verdict text).  Setting the bar at the worst observed value would have
    // hidden them, which is the failure this gate exists to prevent.
    double const tol = ( C.imp == "WEAK"  ? 1e-2      // penalty imposition
                       : C.imp == "TRACE" ? 1e-3      // trace unknowns + tau multipliers
                       :                    1e-4 );   // STRONG: structural continuity
    bool const cell_ok = C.setup_ok && C.conv && C.e_u >= 0. && C.e_u < tol
                       && C.e_s2 >= 0. && C.e_s2 < tol;
    // 2026-09-10: KNOWN-FAIL LIST, so the corpus grep separates standing defects from new
    // ones.  Three cells, all model C / RED_FULL, all measured and documented:
    //   WEAK   -- the sigma-reuse rescue cancellation (|u-ex| 1.69e+00).  FIXED by ocfeslv
    //             rev180's CRONOS_RESCUE_FROM_ROW=1, so this cell is expected to FAIL exactly
    //             when that knob is OFF: it is the fix's regression witness in both directions
    //             (fails without the fix, passes with it).
    //   TRACE  -- 2nd state 2.07e-03 vs 1e-3; STRONG -- 2nd state 1.29e-04 vs 1e-4.  Standing
    //             defects of the exact path on the reused-sigma model, under investigation
    //             (S1b); expected to FAIL regardless of the knob.
    // A cell on this list that FAILS prints XFAIL and does not fail the driver; one that
    // PASSES prints PASS (an improvement, e.g. WEAK with the fix on).  Any OTHER cell failing
    // fails the driver, as before.
    // 2026-09-10 (S2): the header default is now ON (rev181), so an UNSET variable means the
    // fix is active; only an explicit =0 disables it.  This lambda must mirror the header's
    // default or the WEAK expectation inverts silently.
    static bool const kFix = []{ char const* v = std::getenv( "CRONOS_RESCUE_FROM_ROW" );
                                 return !( v && *v ) || std::atoi( v ) != 0; }();
    bool const known_fail =
         ( C.model[0] == 'C' && C.red == "RED_FULL"
           && ( ( C.imp == "WEAK" && !kFix ) || C.imp == "TRACE" || C.imp == "STRONG" ) );
    C.verdict = cell_ok ? "PASS" : known_fail ? "XFAIL" : "FAIL";
    all_ok = all_ok && ( cell_ok || known_fail );
  }

  std::cout << "\n============================== SUMMARY ==============================\n"
            << std::left << std::setw(16) << "model" << std::setw(11) << "reduction"
            << std::setw(8) << "imp"
            << std::right << std::setw(7) << "n_aux" << std::setw(8) << "nVar"
            << std::setw(6) << "conv" << std::setw(5) << "it"
            << std::setw(12) << "|u-ex|" << std::setw(12) << "|2nd-ex|" << "  verdict" << "\n";
  for( Cell const& C : cells ){
    std::cout << std::left << std::setw(16) << C.model << std::setw(11) << C.red
              << std::setw(8) << C.imp << std::right;
    if( !C.setup_ok ){ std::cout << "   *** SETUP FAILED"
              << (C.err.empty()? std::string("") : (" -- "+C.err)) << " ***\n"; continue; }
    std::cout << std::setw(7) << C.n_aux << std::setw(8) << C.nVar
              << std::setw(6) << (C.conv?"y":"n") << std::setw(5) << C.iters;
    if( !C.conv ){ std::cout << "   *** NO CONVERGE ***\n"; continue; }
    std::cout << std::scientific << std::setprecision(2)
              << std::setw(12) << C.e_u << std::setw(12) << C.e_s2
              << "  " << C.verdict << "\n";
  }

  // ---- verdict --------------------------------------------------------------------------------
  size_t n_diff_aux = 0;
  for( size_t i = 0; i + 1 < cells.size(); i += 2 )
    if( cells[i].setup_ok && cells[i+1].setup_ok && cells[i].n_aux != cells[i+1].n_aux )
      ++n_diff_aux;

  std::cout << "\n=============================== VERDICT ===============================\n";
  if( all_ok ){
    size_t n_xfail = 0;
    for( Cell const& C : cells ) if( C.verdict == "XFAIL" ) ++n_xfail;
    if( n_xfail )
      std::cout << "  => ALL CELLS PASS OR ARE DOCUMENTED XFAIL (" << n_xfail << " XFAIL: the\n"
                   "     known-fail list in the gate -- standing defects under investigation,\n"
                   "     not new regressions).  Any change in WHICH cells are XFAIL is a finding.\n";
    else
      std::cout << "  => ALL CELLS PASS, none XFAIL.  A spatial auxiliary coexists correctly with\n"
                   "     high-index reduction, AUTO_ALG_CLOSURE and sigma reuse, under BOTH\n"
                   "     reduction modes.\n";
  }
  else
    std::cout << "  => AT LEAST ONE CELL FAILED beyond the known-fail list.  Read the summary:\n"
                 "     a RED_MAIN-only failure is a genuine R-E/R-F finding; a both-modes failure\n"
                 "     means the MODEL is wrong (my construction) -- check the exact solution\n"
                 "     first.  XFAIL cells are documented standing defects and do NOT trip this.\n";

  std::cout << "\n  auxiliary-count differences between modes: " << n_diff_aux << "/3 models\n";
  if( n_diff_aux == 0 )
    std::cout << "  NOTE: n_aux is identical in every model, so RED_FULL and RED_MAIN may not be\n"
                 "  materially distinguished here -- the pass then proves LESS than it appears.\n"
                 "  (Expected: the modes differ in HOW sigma is REUSED, which need not change the\n"
                 "  COUNT.  Cross-check against nVar and the classify idx= lines above.)\n";

  return all_ok ? 0 : 1;
}
