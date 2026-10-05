// OCFE_PDE20g.cpp  ---  Stage-2 GATE: the SIGMA-REUSE FAMILY of OCFE_PDE20f case C
//
// PDE20f case C found one catastrophic, silent failure: an auxiliary DEFINED through a first
// derivative, under RED_FULL, has that derivative replaced by the reused sigma in its own
// defining row, which (i) turns the row into an apparent identity in the auxiliary and sigma,
// (ii) moves u's z-continuity onto the auxiliary's RESCUED edges, and (iii) lets the forced
// rescue coupling +1 cancel the row's own coefficient of sigma.  PDE20f tests exactly ONE point
// of that family: coefficient -1, linear.  This driver walks the neighbourhood, one axis at a
// time, so the fix (ocfeslv rev180 S1c; rev182 R1+R2 under CRONOS_REUSE_AWARE=1) is tested
// where its assumptions are STRESSED rather than where they were derived:
//
//   G1  0 = w - du/dz           coefficient -1  (the PDE20f point: regression witness)
//   G2  0 = w - 2 du/dz         coefficient -2  (the forced +1 does NOT cancel -2: predicted to
//                                pass WITHOUT the fix, and the fix must not break it; also
//                                settles the cf=2 discrepancy recorded 2026-09-10h)
//   G3  0 = w + 0.5 du/dz       coefficient +0.5 (positive partner: +1 reinforces; c0 = +0.5)
//   G4  0 = w - du/dz - u       mixed bare terms: parent u bare beside the alias -- the R2
//                                alias test must fire on Dz_u and not be confused by u
//   G5  0 = w - EPS*(du/dz)^2   NON-LINEAR in the derivative, EPS = 0.05 by default.  MEASURED
//                                2x2 (2026-09-10): {u_z has zeros | strictly monotone reference}
//                                x {EPS=1 | EPS=0.05}.  EPS=1 never converges properly in either
//                                reference (K*w ~ 30 swamps the linear operator); EPS=0.05
//                                converges in 8 Newton steps to u ~ 4e-06 in BOTH references.
//                                The zeros of u_z -- where d(ALG)/d(Dz_u) = -2*EPS*Dz_u vanishes
//                                -- are HARMLESS to convergence here.  PDE20G_BETA (reference
//                                shift) and PDE20G_EPS remain as runtime knobs for that 2x2.
//                                With the pair OFF, RED_FULL/TRACE is NOT equivalent to RED_MAIN
//                                (2.5x-39x on the second state); with CRONOS_REUSE_AWARE=1 all six
//                                cells are equivalent to the digit -- so R2 (alias-aware drop)
//                                alone suffices on the non-constant path; R1 is inert there by
//                                design (c0 non-constant -> S1c falls back).  c0 is
//                                NON-CONSTANT, which is G5's point:
//                                S1c must fall back (counted) and the forced rescue stands
//
// Cells: 5 models x 2 reductions x 3 impositions = 30.  Same discretisation, exact solution,
// tolerances and verdict logic as PDE20f.  KNOWN-FAIL LIST: EMPTY on first delivery -- every
// XFAIL entry must be EARNED by a measurement and an explanation, as PDE20f's were.  Read the
// first run as a map, then decide.
//
// Exact: u = exp(-t)(1+sin(2 pi z + phase)); w_ex follows from the defining relation.
#include <iostream>
#include <cstdlib>
#include <iomanip>
#include <vector>
#include <cmath>
#include <string>
#include <sstream>
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
// G5 experiment knobs (runtime, so the 2x2 {zero of u_z} x {magnitude of K*w} runs without rebuilds):
//   PDE20G_BETA   reference shift u_ex + BETA*z  (0: u_z has zeros; >2pi: strictly monotone)
//   PDE20G_EPS    nonlinearity scale, ALG = w - EPS*(du/dz)^2  (1: K*w ~ 30; 0.05: K*w ~ 1)
//   PDE20G_MODELS comma-separated model indices to run, e.g. "4" for G5 only
static double knob( char const* name, double def ){ char const* v = std::getenv( name ); return ( v && *v ) ? std::atof( v ) : def; }
static double const BETA_G5 = knob( "PDE20G_BETA", 0.0 );   // MEASURED harmless: default OFF
static double const EPS_G5  = knob( "PDE20G_EPS",  0.05 );  // MEASURED decisive: default 0.05
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
  double wscale=1.;                         // max|w_ex| over the sampled nodes (relative bar)
  double e_u=-1., e_s2=-1.;                 // u error, and the secondary state's error
  std::string s2name="-";
  std::string err="";                       // exception text, if the cell threw
  std::string verdict="";                   // 2026-09-10: PASS / XFAIL / FAIL, set by the gate
};

// ---------------------------------------------------------------------------------------------
static Cell run_case( int model, O::ReductionType red, O::ImpositionType imp, size_t nn_z_override = 0 )
{
  Cell C;
  static char const* NAMES[5] = { "G1 coef -1", "G2 coef -2", "G3 coef +.5", "G4 mixed u", "G5 nonlin" };
  C.model = NAMES[model];
  if( nn_z_override ) C.model += " p" + std::to_string( nn_z_override );   // the p-refined convergence cell
  C.red   = ( red==O::RED_FULL ? "RED_FULL" : "RED_MAIN" );
  C.imp   = ( imp==O::IC_WEAK ? "WEAK" : imp==O::IC_TRACE ? "TRACE" : "STRONG" );

  size_t const ne_t=3, nn_t=5, ne_z=6, nn_z=( nn_z_override ? nn_z_override : 6 );
  std::vector<double> t_bnd, z_bnd;
  for( size_t i=0;i<=ne_t;++i ) t_bnd.push_back( T_END*double(i)/double(ne_t) );
  for( size_t i=0;i<=ne_z;++i ) z_bnd.push_back( double(i)/double(ne_z) );

  FFGraph DAG;
  FFVar t=DAG.add_var("t"), z=DAG.add_var("z"), u=DAG.add_var("u(t,z)");
  FFVar x=DAG.add_var("x(t,z)"), y=DAG.add_var("y(t,z)"), w=DAG.add_var("w(t,z)");
  FFPartial OpP;

  FFVar W  = PI2*z+PH;
  double const BETA = ( model==4 ? BETA_G5 : 0.0 );      // G5 only: strictly monotone reference
  FFVar UE = exp(-t)*(1.0+sin(W)) + BETA*z;              // u_ex
  FFVar UT = -exp(-t)*(1.0+sin(W));                      // du_ex/dt
  FFVar UZ = exp(-t)*PI2*cos(W) + BETA;                  // du_ex/dz
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
  oc.update_ref( u, [&t,&z,BETA](OCFESLV::t_Coord const& c){ return uex(c.at(z),c.at(t)) + BETA*c.at(z); } );

  oc.add_state( w, {t,z} );
  oc.update_ref( w, [&t,&z,model,BETA](OCFESLV::t_Coord const& c){
      double const uz = uzex(c.at(z),c.at(t)) + BETA, ue = uex(c.at(z),c.at(t)) + BETA*c.at(z);
      switch( model ){ case 0: return uz; case 1: return 2.0*uz; case 2: return -0.5*uz;
                       case 3: return uz + ue; default: return EPS_G5*uz*uz; } } );
  C.s2name = "w";
  oc.set_evolution_domain( t );

  FFVar SRC, PDEu, ALG;
  FFVar WEX = ( model==0 ? UZ : model==1 ? 2.0*UZ : model==2 ? -0.5*UZ
              : model==3 ? UZ + UE : EPS_G5*UZ*UZ );   // w_ex as a DAG expression
  SRC  = LHS - K_LIN*WEX;                              // u_ex is exact once w == w_ex
  PDEu = OpP(u,t) + U_ADV*OpP(u,z) - D_AX*OpP(OpP(u,z),z) - K_LIN*w - SRC;
  switch( model ){
    case 0:  ALG = w - OpP(u,z);               break;
    case 1:  ALG = w - 2.0*OpP(u,z);           break;
    case 2:  ALG = w + 0.5*OpP(u,z);           break;
    case 3:  ALG = w - OpP(u,z) - u;           break;
    default: ALG = w - EPS_G5*OpP(u,z)*OpP(u,z); break;
  }
  FFVar ICu = u - UE;
  FFVar BCL = U_ADV*u - D_AX*OpP(u,z) - ( U_ADV*UE - D_AX*UZ );   // flux/Robin at z=0
  FFVar BCU = OpP(u,z) - UZ;                                      // Neumann at z=1

  oc.add_equation( PDEu, {t,z}, {TNL,ZI},               EO(ER::INTERIOR,0) );
  oc.add_equation( ICu , {t,z}, {FFDom::LB,FFDom::ALL}, EO(ER::INITIAL ,0) );
  oc.add_equation( BCL , {t,z}, {TNL,FFDom::LB},        EO(ER::BOUNDARY,0) );
  oc.add_equation( BCU , {t,z}, {TNL,FFDom::UB},        EO(ER::BOUNDARY,0) );

  oc.add_equation( ALG, {t,z}, {FFDom::ALL,FFDom::ALL},     EO(ER::INTERIOR,0) );

  oc.options.REDUCE.ORDER     = red;
  oc.options.CLASSIFY.MODE         = O::CLASS_AUTO;
  oc.options.INTERFACE.TYPE   = O::IC_AUTO;
  oc.options.INTERFACE.IMPOSITION  = imp;
  oc.options.AUTO.DIFF_ELIM   = true;                    // ON for all: measured inert on PDE20f C, and
                                                         // to eliminate NOTHING there (see the
                                                         // case-C header); kept on as a guard
  // 2026-09-09 PROBE: PDE20F_SIGMA0 overrides the weak-imposition penalty strength.  Model C
  // under RED_FULL/IC_WEAK gives |u-ex| = 1.69e+00 while converging; the discriminator is
  // whether that error SCALES with the penalty (under-penalisation -- the claim is present but
  // too weak) or is INSENSITIVE to it (the claim is absent, misdirected, or cancelling).
  oc.options.INTERFACE.SAT_SIGMA0       = 10.0;
  if( char const* v = std::getenv( "PDE20F_SIGMA0" ) ) oc.options.INTERFACE.SAT_SIGMA0 = std::atof( v );
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
        eu = std::max( eu, std::fabs( V(u,zz,tt) - ( uex(zz,tt) + BETA*zz ) ) );
        { double const uz = uzex(zz,tt) + BETA, ue = uex(zz,tt) + BETA*zz;
          double const wex = ( model==0 ? uz : model==1 ? 2.0*uz : model==2 ? -0.5*uz
                             : model==3 ? uz + ue : EPS_G5*uz*uz );  // MUST match WEX above
          e2 = std::max( e2, std::fabs( V(w,zz,tt) - wex ) );
          C.wscale = std::max( C.wscale, std::fabs( wex ) ); }
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
  std::vector<int> models = {0,1,2,3,4};
  if( char const* v = std::getenv( "PDE20G_MODELS" ) ){ models.clear(); std::string t; std::istringstream is(v);
    while( std::getline( is, t, ',' ) ) if( !t.empty() ) models.push_back( std::atoi( t.c_str() ) ); }
  for( int m : models ){
    for( O::ImpositionType imp : IMPS ){
      cells.push_back( run_case(m,O::RED_FULL,imp) );
      cells.push_back( run_case(m,O::RED_MAIN,imp) );
    }
  }
  // ---- 2026-09-17: THE G4 CONVERGENCE CELLS (one p-refinement, RED_MAIN, each imposition) --------------------
  // G4's RED_MAIN evaluates a bare one-sided du/dz beside the auxiliary in w's row.  Under the receiver rule
  // (rev277 WEAK, rev278/279 exact modes) its second-state error at the base mesh (6 elements x 6 nodes) is the
  // HONEST 1.93e-03 -- above the exact-mode bar, which was calibrated while the rescued injection into the LINK
  // row was tying du/dz to the accurate auxiliary (a coarse-mesh artefact whose size is non-monotone in the
  // rescue's fabricated weight).  Refinement decided which it is (NOTES_20260917l): w at the z-interfaces,
  // rule on, RED_MAIN --  h: 2.05e-03 (6 el) -> 5.7e-05 (12) -> 4.3e-06 (24);  p: 2.9e-06 (8 nodes) ->
  // 2.77e-06 (10) -- a spectral collapse onto the t-discretisation floor that every variant shares.  So the
  // number is a discretisation error at that mesh, not a defect, and the assertion that decides it is
  // convergence: RED_MAIN with 8 nodes per z-element must put w below 1e-5 in every imposition.  The base-mesh
  // G4 cells keep their honest values with a G4-specific second-state bar (below).
  bool const refine_g4 = std::find( models.begin(), models.end(), 3 ) != models.end();
  if( refine_g4 )
    for( O::ImpositionType imp : IMPS ) cells.push_back( run_case( 3, O::RED_MAIN, imp, 8 ) );

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
    // 2026-09-10 calibration: the ABSOLUTE bar on |2nd-ex| is what PDE20f uses, and it is right
    // there (max|w_ex| ~ 2 pi).  Across this family w_ex ranges from 0.5*u_z (G3) to 2*u_z (G2)
    // to u_z+u (G4), and an absolute 1e-4 under STRONG flagged G2 and G4 in BOTH reductions at
    // the SAME value -- i.e. not a reuse effect but the bar.  So the second state is gated
    // RELATIVE to max|w_ex| (scale >= 1 so small-w models are not loosened).  Absolute
    // errors stay in the table so G1 remains directly comparable to PDE20f case C.
    double s2bar = tol * std::max( 1.0, C.wscale );
    bool const g4_refined = ( C.model == "G4 mixed u p8" );
    if( g4_refined ) s2bar = 1e-5;                     // the convergence assertion (see the cell comment above)
    else if( C.model == "G4 mixed u" && C.red == "RED_MAIN" )
      s2bar = std::max( s2bar, 5e-3 );                 // the honest base-mesh value, 1.93e-03; bound by p8 above
    bool const cell_ok = C.setup_ok && C.conv && C.e_u >= 0. && C.e_u < tol
                       && C.e_s2 >= 0. && C.e_s2 < s2bar;
    bool const known_fail = false;
    C.verdict = cell_ok ? "PASS" : known_fail ? "XFAIL" : "FAIL";
    all_ok = all_ok && ( cell_ok || known_fail );
  }

  // ---- THE FAMILY INVARIANT: RED_FULL must be EQUIVALENT to RED_MAIN, cell by cell ----------
  // Sigma reuse is a REPRESENTATION choice; it must not change the answer.  Measured
  // (NOTES_20260910m): without the S1b pair it does, on the exact path (G1 TRACE 48x, G2
  // TRACE ~1000x, G1/G3 STRONG 2-3x); with CRONOS_REUSE_AWARE=1 every pair is equivalent to
  // the digit.  A non-equivalent pair is the reuse defect: XFAIL when the pair is OFF
  // (documented, fixed by rev182), FAIL when it is ON.
  // rev183 mirrors the header default: UNSET means the pair is ON; only an explicit =0 disables
  // it.  Without this the equivalence expectation would invert silently after the flip.
  static bool const kPair = []{ char const* v = std::getenv( "CRONOS_REUSE_AWARE" );
                                return !( v && *v ) || std::atoi( v ) != 0; }();
  double const EQ_FACTOR = 2.0;
  std::vector<std::string> eqv;   // one entry per (model, imp), in cell order
  for( size_t i = 0; i + 1 < cells.size(); i += 2 ){
    Cell& RF = cells[i]; Cell& RM = cells[i+1];
    if( RF.red != "RED_FULL" || RM.red != "RED_MAIN" || RF.model != RM.model ) break;   // the p8 cells follow the pairs
    bool ok = RF.conv && RM.conv && RF.e_u >= 0. && RM.e_u >= 0.;
    double ru = 0., r2 = 0.;
    if( ok ){ ru = RF.e_u / std::max( RM.e_u, 1e-300 ); r2 = RF.e_s2 / std::max( RM.e_s2, 1e-300 ); }
    // 2026-09-17 (rev277): the gate is TWO-SIDED.  Until now it only caught RED_FULL being worse than RED_MAIN,
    // so a 0.025x second-state ratio printed EQUIV.  Equivalence means both ratios in [1/EQ_FACTOR, EQ_FACTOR].
    auto within = []( double r, double f ){ return r <= f && r >= 1.0/f; };
    bool const equiv = ok && within( ru, EQ_FACTOR ) && within( r2, EQ_FACTOR );
    // The ONE expected non-equivalence: G4 under IC_WEAK.  RED_MAIN keeps the bare one-sided du/dz beside the
    // auxiliary in w's row, and at interface nodes that derivative is accurate only to ~2e-03 (NOTES_20260917f);
    // RED_FULL reuses the auxiliary and is ~40x better.  The default's smaller gap came from the rescued
    // injection into the LINK row, whose effect is non-monotone in a fabricated number (coupling 0.3/1/7 ->
    // 2.7e-03/4.7e-04/1.3e-03).  Expected, documented, and not a defect of either reduction.
    // The exception is a MODEL fact, not a mode fact: RED_MAIN evaluates a bare one-sided du/dz at interface
    // nodes in w's row whatever the imposition.  Measured under the receiver rule (2026-09-17): WEAK 2.47e-02x,
    // TRACE and STRONG 2.44e-02x -- the three impositions agree once none of them injects into the LINK
    // duplicate.  Without the rule the exact modes print ~1.00x, and that agreement is the artefact.
    bool const expected_gap = ( RF.model == "G4 mixed u" && ok && within( ru, EQ_FACTOR ) && r2 < 1.0 );
    std::ostringstream os;
    os << std::left << std::setw(13) << RF.model << std::setw(8) << RF.imp << std::scientific
       << std::setprecision(2) << "  u: " << ru << "x   2nd: " << r2 << "x   "
       << ( equiv ? "EQUIV" : expected_gap ? "NOT EQUIV (expected: RED_MAIN bare du/dz at interface nodes)"
          : kPair ? "** NOT EQUIV **  FAIL" : "** NOT EQUIV **  XFAIL (reuse defect; rev182 pair OFF)" );
    eqv.push_back( os.str() );
    if( !equiv && !expected_gap ){ RF.verdict = kPair ? "FAIL" : "XFAIL"; if( kPair ) all_ok = false; }
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
  for( size_t i = 0; i + 1 < cells.size(); i += 2 ){
    if( cells[i].red != "RED_FULL" || cells[i+1].red != "RED_MAIN" ) break;   // pairs only (the p8 cells follow)
    if( cells[i].setup_ok && cells[i+1].setup_ok && cells[i].n_aux != cells[i+1].n_aux )
      ++n_diff_aux; }

  std::cout << "\n---- EQUIVALENCE  (RED_FULL vs RED_MAIN, same model & imposition; factor "
            << EQ_FACTOR << ") ----\n";
  for( auto const& l : eqv ) std::cout << "  " << l << "\n";
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
