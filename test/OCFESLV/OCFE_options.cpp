// OCFE_options.cpp -- the Options / EqnOptions MECHANISMS, tested directly (2026-09-21).
//
// Why this driver exists.  Every failure in the rev308-rev310 refactor was invisible to the gate used at the
// time (OCFE_PDE6_blk0, a physics model) and surfaced only in a corpus sweep or by accident:
//   - OCFESLV::deep_copy_from carried FFModel's options but not the solver's, so a COPIED environment kept the
//     source's frozen interface plan beside default-constructed solver options -- 8 programs aborted;
//   - OCFE_exceptions_api read stored options by value after they became a shared_ptr, and did not compile.
// Both are failures of MECHANISM, not of numerics.  This driver checks the mechanisms themselves, on a tiny
// linear model, so each refactor step can be gated in well under a second.
//
// What it asserts (each item names the failure or design decision it guards):
//   O1  requested vs applied: a model option set on oc.options reaches FFModel::options only at setup();
//   O2  post-setup mutation of a model option does nothing until the next setup() (the COPY semantics
//       chosen over a pointer in rev308, which nothing else asserts);
//   O3  copy constructor, operator= and deep_copy_from carry BOTH the solver's and the model's options
//       (the rev309 bug);
//   O4  the eight retired knobs are environment-only and default to their measured values;
//   O5  CRONOS_DISPLAY_LEVEL sets both layers;
//   O6  the frozen-plan guard trips when IMPOSITION_TYPE changes after setup();
//   O8  rev314's fourteen environment-only knobs keep their rev313 defaults (verified field-for-accessor
//       against rev313 separately; asserted here so the driver guards it from now on);
//   O9  rev316: SAT_SIGMA1 and ALG_CLOSURE environment-only -- sigma1 defaults to 1.0 through the rev51 path
//       (read once per process), ALG_CLOSURE defaults on and follows CRONOS_AUTO_ALG_CLOSURE;
//   O10 rev318: HYP_CLOSURE environment-only -- default on, CRONOS_AUTO_HYP_CLOSURE=0 turns it off;
//   O7  rev313's enum move: SKIP_NT_ON_TRACE keeps its in-class default after leaving FFModel::Options (it had
//       NO reset() line, so a lost initializer would have left it indeterminate), and the collapsed
//       InterfaceType is ONE type -- OCFESLV's and OCPlan's enumerators are the same values;
//   E1  scalar AND vector add_equation store every equation (the vector overloads call add_equation on this,
//       which is why it had to become virtual);
//   E2  role-derived defaults of the stored options;
//   E3  stored options read back through var_equation() by pointer (the rev310 storage change);
//   E4  a copied environment OWNS its equation options -- mutating the copy's leaves the source's untouched.
//       rev310's copy loop shared the shared_ptr, coupling the two; harmless only while nothing mutated in
//       place, and step B's in-place normaliser would have armed it (fixed by EqnOptions::clone(), rev311);
//   E5  an equation MINTED BY THE MODEL (a reduce_order LINK row) carries the solver's options type -- the
//       case the rev312 factory exists for, since those rows never pass through add_equation;
//   E6  a COPY's stored options keep the solver's type and field values -- the clone() override contract.
#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <cmath>
#include <string>
#include <cstdlib>
#include <type_traits>

#ifndef OCFE_OCFESLV_HEADER
#define OCFE_OCFESLV_HEADER "ocfeslv.hpp"
#endif
#include OCFE_OCFESLV_HEADER
using namespace mc;

static int g_pass=0, g_fail=0;
static void check( std::string const& nm, bool ok, std::string const& detail = "" )
{
  (ok?g_pass:g_fail)++;
  std::cout << "  " << std::left << std::setw(66) << nm << std::right << (ok?"PASS":"FAIL")
            << ( detail.empty() ? "" : "   " + detail ) << "\n";
}

struct Vars { FFVar z, u, w; };

// A linear two-state problem: u'' = 1 on (0,1) with u(0)=0, u(1)=1, and an algebraic w = 2u.
static void build( FFGraph& DAG, OCFESLV& oc, Vars& V, bool vector_form = false )
{
  V.z = DAG.add_var( "z" );
  V.u = DAG.add_var( "u(z)" );
  V.w = DAG.add_var( "w(z)" );
  FFPartial OpP;
  oc.add_domain( V.z, FFDom( 0., 1., 2, FFDom::LGL, 5 ) );
  oc.add_state ( V.u, { V.z } );
  oc.add_state ( V.w, { V.z } );
  int const Z_INT = FFDom::ALL - FFDom::LB - FFDom::UB;
  typedef OCFESLV::EqnOptions EO; typedef OCFESLV::EqnRole ER;
  if( !vector_form ){
    oc.add_equation( OpP( OpP( V.u, V.z ), V.z ) - 1.0, { V.z }, { Z_INT },     EO( ER::INTERIOR, 0 ) );
    oc.add_equation( V.w - 2.0*V.u,                     { V.z }, { FFDom::ALL }, EO( ER::INTERIOR, 0 ) );
  }
  else{
    // the VECTOR overload: FFModel loops `add_equation( e, vDom, vLim, opt )` on itself, which is exactly the
    // call a non-virtual override would have silently missed
    std::vector<FFVar> const vE = { OpP( OpP( V.u, V.z ), V.z ) - 1.0 };
    oc.add_equation( vE, { V.z }, { Z_INT }, EO( ER::INTERIOR, 0 ) );
    std::vector<FFVar> const vA = { V.w - 2.0*V.u };
    oc.add_equation( vA, { V.z }, { FFDom::ALL }, EO( ER::INTERIOR, 0 ) );
  }
  oc.add_equation( V.u,       { V.z }, { FFDom::LB }, EO( ER::BOUNDARY, 0 ) );
  oc.add_equation( V.u - 1.0, { V.z }, { FFDom::UB }, EO( ER::BOUNDARY, 0 ) );
}

static void configure( OCFESLV& oc, OCFESLV::Options::ImpositionType imp )
{
  oc.options.DISPLAY_LEVEL   = 0;
  oc.options.REDUCE.ORDER    = OCFESLV::Options::RED_FULL;
  oc.options.CLASSIFY.MODE        = OCFESLV::Options::CLASS_AUTO;
  oc.options.INTERFACE.IMPOSITION = imp;
}

int main()
{
  std::cout << "================================================================\n"
            << "  Options / EqnOptions mechanisms          header: " << OCFE_OCFESLV_HEADER << "\n"
            << "================================================================\n";

  // ---- O1 requested vs applied ------------------------------------------------------------------------------
  std::cout << "\n---- O1: requested vs applied ----\n";
  {
    FFGraph DAG; OCFESLV oc( &DAG ); Vars V; build( DAG, oc, V );
    configure( oc, OCFESLV::Options::IC_TRACE );
    oc.options.REDUCE.ORDER = OCFESLV::Options::RED_MAIN;          // a MODEL option, set on the solver's object
    FFModel const& fm = oc;
    bool const before = ( fm.options.REDUCE.ORDER != OCFESLV::Options::RED_MAIN );
    bool const ok = oc.setup();
    bool const after  = ( fm.options.REDUCE.ORDER == OCFESLV::Options::RED_MAIN );
    check( "O1 a model option set on oc.options is NOT applied before setup()", before );
    check( "O1 setup() succeeds", ok );
    check( "O1 ...and IS applied to FFModel::options by setup()", after );
  }

  // ---- O2 post-setup mutation ------------------------------------------------------------------------------
  std::cout << "\n---- O2: a model option changed after setup() waits for the next setup() ----\n";
  {
    FFGraph DAG; OCFESLV oc( &DAG ); Vars V; build( DAG, oc, V );
    configure( oc, OCFESLV::Options::IC_TRACE );
    oc.setup();
    FFModel const& fm = oc;
    oc.options.CLASSIFY.MODE = OCFESLV::Options::CLASS_NONE;             // requested, after setup
    bool const unchanged = ( fm.options.CLASSIFY.MODE == OCFESLV::Options::CLASS_AUTO );
    oc.setup();
    bool const now = ( fm.options.CLASSIFY.MODE == OCFESLV::Options::CLASS_NONE );
    check( "O2 the applied value is unchanged until the next setup()", unchanged );
    check( "O2 ...and changes at the next setup()", now );
  }

  // ---- O3 copies carry both halves --------------------------------------------------------------------------
  std::cout << "\n---- O3: copy ctor, operator= and deep_copy_from carry BOTH option objects ----\n";
  {
    FFGraph DAG; OCFESLV src( &DAG ); Vars V; build( DAG, src, V );
    configure( src, OCFESLV::Options::IC_STRONG );
    src.options.SOLVE.MAX_ITER = 37;                                 // a SOLVER option
    src.options.REDUCE.ORDER   = OCFESLV::Options::RED_FULL;         // a MODEL option
    src.options.REDUCE.HIDDEN_IC = true;                             // rev337, a MODEL option
    src.options.TTOL             = 5.0e-7;                           // rev345, a MODEL option
    bool const ok = src.setup();
    check( "O3 source setup() succeeds (IC_STRONG)", ok );
    auto same = []( OCFESLV const& a, OCFESLV const& b ){
      FFModel const& fa = a; FFModel const& fb = b;
      // NOTE: these are REPRESENTATIVES, not the whole object -- a field added to Options and forgotten in
      // operator= is caught here only once it is listed.  REDUCE.HIDDEN_IC (rev337) and TTOL (rev345) were both
      // added without a check and passed this driver unchanged; anything new belongs on this list.
      return a.options.INTERFACE.IMPOSITION == b.options.INTERFACE.IMPOSITION
          && a.options.SOLVE.MAX_ITER    == b.options.SOLVE.MAX_ITER
          && fa.options.REDUCE.ORDER     == fb.options.REDUCE.ORDER
          && fa.options.REDUCE.HIDDEN_IC == fb.options.REDUCE.HIDDEN_IC
          && fa.options.TTOL             == fb.options.TTOL; };
    { OCFESLV c( src );                        check( "O3 copy constructor carries both halves", same( c, src ) ); }
    { FFGraph D2; OCFESLV a( &D2 ); a = src;   check( "O3 operator= carries both halves",        same( a, src ) ); }
    { OCFESLV d; bool const r = d.deep_copy_from( src );
      check( "O3 deep_copy_from succeeds and carries both halves", r && same( d, src ) ); }
    // the rev309 failure mode, directly: a copy must EVALUATE without tripping the frozen-plan guard
    { OCFESLV c( src ); std::vector<double> xv, inp;
      bool const ev = c.init( xv, inp, nullptr ) && c.solve( xv.data(), inp.data(), nullptr ).converged;
      check( "O3 a copy of an IC_STRONG setup solves (the rev309 guard failure)", ev ); }
  }

  // ---- O11 the newest model options: defaults, and application at setup ---------------------------------------
  std::cout << "\n---- O11: REDUCE.HIDDEN_IC (rev337) and TTOL (rev345) ----\n";
  {
    FFModel::Options fresh;
    check( "O11 REDUCE.HIDDEN_IC defaults to false", fresh.REDUCE.HIDDEN_IC == false );
    check( "O11 TTOL defaults to 1e-9",              fresh.TTOL == 1.0e-9,
           [&]{ std::ostringstream os; os << std::scientific << std::setprecision(3) << "got " << fresh.TTOL;
                return os.str(); }() );

    FFGraph DAG; OCFESLV oc( &DAG ); Vars V; build( DAG, oc, V );
    configure( oc, OCFESLV::Options::IC_WEAK );
    oc.options.TTOL             = 2.5e-6;
    oc.options.REDUCE.HIDDEN_IC = true;
    FFModel const& mdl = oc;
    bool const before = ( mdl.options.TTOL != 2.5e-6 );          // not applied before setup(), as O1 established
    bool const ok     = oc.setup();
    bool const after  = ( mdl.options.TTOL == 2.5e-6 && mdl.options.REDUCE.HIDDEN_IC == true );
    check( "O11 they are not applied before setup()", before );
    check( "O11 setup() succeeds",                    ok );
    check( "O11 ...and both reach FFModel::options",  after );
  }

  // ---- O4 retired knobs -------------------------------------------------------------------------------------
  std::cout << "\n---- O4: the eight retired knobs are environment-only, measured defaults ----\n";
  {
    char const* keys[] = { "CRONOS_RESCUE_C2", "CRONOS_WEAK_TAU_RESCUED", "CRONOS_WEAK_TAU_NATURAL",
                           "CRONOS_EXACT_NATURAL", "CRONOS_EXACT_NATURAL_STRONG", "CRONOS_STRONG_PROJECT",
                           "CRONOS_WEAK_NATURAL_PENALTY", "CRONOS_DROP_WDECIDE" };
    for( auto k : keys ) unsetenv( k );
    check( "O4 RESCUE_C2 defaults to true",            OCFESLV::_knob_RESCUE_C2() == true );
    check( "O4 WEAK_TAU_RESCUED defaults to 0",        OCFESLV::_knob_WEAK_TAU_RESCUED() == 0 );
    check( "O4 WEAK_TAU_NATURAL defaults to false",    OCFESLV::_knob_WEAK_TAU_NATURAL() == false );
    check( "O4 EXACT_NATURAL defaults to 2",           OCFESLV::_knob_EXACT_NATURAL() == 2 );
    check( "O4 EXACT_NATURAL_STRONG defaults to true", OCFESLV::_knob_EXACT_NATURAL_STRONG() == true );
    check( "O4 STRONG_PROJECT defaults to true",       OCFESLV::_knob_STRONG_PROJECT() == true );
    check( "O4 WEAK_NATURAL_PENALTY defaults to 2",    OCFESLV::_knob_WEAK_NATURAL_PENALTY() == 2 );
    check( "O4 DROP_WDECIDE defaults to true",         OCFESLV::_knob_DROP_WDECIDE() == true );
    setenv( "CRONOS_EXACT_NATURAL", "0", 1 );
    check( "O4 ...and follow the environment when set (EXACT_NATURAL=0)", OCFESLV::_knob_EXACT_NATURAL() == 0 );
    unsetenv( "CRONOS_EXACT_NATURAL" );
  }

  // ---- O5 DISPLAY_LEVEL on both layers (rev314: resolved at READ time) ----------------------------------------
  // rev308 added SOLVE_DISPLAY_LEVEL and never read it.  The obvious wiring -- copy DISPLAY_LEVEL at reset() --
  // would have been WRONG: drivers set DISPLAY_LEVEL AFTER construction, so the copy would be stale and every
  // such driver's solver output would silently change.  So it defaults to -1 = follow, resolved when read.
  std::cout << "\n---- O5: the solver's display level follows the model's unless set ----\n";
  {
    FFGraph DAG; OCFESLV oc( &DAG ); Vars V; build( DAG, oc, V );
    check( "O5 SOLVE_DISPLAY_LEVEL defaults to -1 (follow)", oc.options.SOLVE.DISPLAY_LEVEL == -1 );
    oc.options.DISPLAY_LEVEL = 3;                                   // set AFTER construction, as drivers do
    check( "O5 ...so the solver follows a DISPLAY_LEVEL set after construction", oc._solve_display_level() == 3,
           "solver level=" + std::to_string( oc._solve_display_level() ) );
    oc.options.SOLVE.DISPLAY_LEVEL = 0;
    check( "O5 ...and an explicit SOLVE_DISPLAY_LEVEL makes it independent",
           oc._solve_display_level() == 0 && oc.options.DISPLAY_LEVEL == 3 );
  }

  // ---- O6 the frozen-plan guard -----------------------------------------------------------------------------
  std::cout << "\n---- O6: the frozen-plan guard trips when IMPOSITION_TYPE changes after setup() ----\n";
  {
    FFGraph DAG; OCFESLV oc( &DAG ); Vars V; build( DAG, oc, V );
    configure( oc, OCFESLV::Options::IC_STRONG );
    oc.setup();
    oc.options.INTERFACE.IMPOSITION = OCFESLV::Options::IC_WEAK;          // changed without re-running setup()
    std::vector<double> xv, inp; bool threw = false, solved = false;
    try{ solved = oc.init( xv, inp, nullptr ) && oc.solve( xv.data(), inp.data(), nullptr ).converged; }
    catch( ... ){ threw = true; }
    check( "O6 an evaluation with a changed imposition does NOT silently solve", threw || !solved );
  }

  // ---- O7 the enum move ------------------------------------------------------------------------------------
  std::cout << "\n---- O7: the enum move (rev313) ----\n";
  {
    unsetenv( "CRONOS_SKIP_NT_ON_TRACE" );
    // (rev314: SKIP_NT_ON_TRACE is environment-only now; its default is asserted through the accessor in O8.)
    bool const one_type = std::is_same<OCFESLV::Options::InterfaceType, OCPlan::InterfaceType>::value;
    check( "O7 InterfaceType is ONE type (OCFESLV aliases OCPlan's)", one_type );
    check( "O7 ...and its enumerators agree (IC_VALUE..IC_AUTO = 0..3; IC_CONS retired at rev320)",
           OCFESLV::Options::IC_VALUE == 0 && OCFESLV::Options::IC_FLUX == 1 && OCFESLV::Options::IC_UPWIND == 2
        && OCFESLV::Options::IC_AUTO == 3 );
    // rev319: ImpositionType collapsed into OCPlan's enum, as InterfaceType was at rev313
    check( "O7 ImpositionType is ONE type (OCFESLV aliases OCPlan's, rev319)",
           std::is_same<OCFESLV::Options::ImpositionType, OCPlan::ImpositionSentinel>::value );
    check( "O7 ...and its enumerators agree (IC_WEAK..IC_TRACE = 0..2)",
           OCFESLV::Options::IC_WEAK == 0 && OCFESLV::Options::IC_STRONG == 1 && OCFESLV::Options::IC_TRACE == 2 );
  }

  // ---- O8 rev314's fourteen environment-only knobs ------------------------------------------------------------
  std::cout << "\n---- O8: the fourteen rev314 knobs, environment-only, at their measured defaults ----\n";
  {
    for( auto k : { "CRONOS_LINK_PAIR_SUPPRESS", "CRONOS_NO_PIN_CLOSURE", "CRONOS_SENS_TAU_REG", "CRONOS_SENS_MINNORM_PINV",
                    "CRONOS_CHAIN_SYMBOL", "CRONOS_CLAIM_C0", "CRONOS_AUDIT_SETUP_ONLY", "CRONOS_AUDIT_BETA",
                    "CRONOS_AUDIT_SENS", "CRONOS_DUP_SPREAD", "CRONOS_DUP_SPREAD_WARN",
                    "CRONOS_DETERMINACY_GAP_MIN", "CRONOS_DETERMINACY_MAX_DENSE" } ) unsetenv( k );
    check( "O8 LINK_PAIR off, PIN_CLOSURE on",   !OCFESLV::_knob_INTERFACE_SUPPRESS_LINK_PAIR() && OCFESLV::_knob_INTERFACE_SUPPRESS_PIN_CLOSURE() );
    check( "O8 SKIP_NT_ON_TRACE = NT_SKIP_ALWAYS", OCFESLV::_knob_SKIP_NT_ON_TRACE() == OCFESLV::Options::NT_SKIP_ALWAYS );
    check( "O8 SENS_TAU_REG 0, SENS_MINNORM_PINV off", OCFESLV::_knob_SENS_TAU_REG() == 0.0 && !OCFESLV::_knob_SENS_MINNORM_PINV() );
    check( "O8 CHAIN_SYMBOL on (model side)",     FFModel::_knob_CHAIN_SYMBOL() );
    check( "O8 CLAIM_FROM_C0 off",                !OCFESLV::_knob_CLAIM_FROM_C0() );
    check( "O8 the four instrument audits off",   !OCFESLV::_knob_AUDIT_SETUP_ONLY() && !OCFESLV::_knob_AUDIT_BETA()
                                                && !OCFESLV::_knob_AUDIT_SENS_NULL() && !OCFESLV::_knob_AUDIT_DUP_SPREAD() );
    check( "O8 DUP_SPREAD_WARN 1e-10",            OCFESLV::_knob_DUP_SPREAD_WARN() == 1e-10 );
    check( "O8 DETERMINACY GAP_MIN 1e2, MAX_DENSE 3000", OCFESLV::_knob_DETERMINACY_GAP_MIN() == 1e2
                                                && OCFESLV::_knob_DETERMINACY_MAX_DENSE() == 3000 );
    setenv( "CRONOS_NO_PIN_CLOSURE", "1", 1 );
    check( "O8 ...and the inverted CRONOS_NO_PIN_CLOSURE still inverts", !OCFESLV::_knob_INTERFACE_SUPPRESS_PIN_CLOSURE() );
    unsetenv( "CRONOS_NO_PIN_CLOSURE" );
  }

  // ---- O9 rev316 ---------------------------------------------------------------------------------------------
  std::cout << "\n---- O9: SAT_SIGMA1 and ALG_CLOSURE are environment-only (rev316) ----\n";
  {
    FFGraph DAG; OCFESLV oc( &DAG ); Vars V; build( DAG, oc, V );
    // CRONOS_SAT_SIGMA1 is read ONCE per process (rev51), so only its default can be asserted in-process
    check( "O9 sigma1 defaults to 1.0 (CRONOS_SAT_SIGMA1 unset at start)",
           std::getenv( "CRONOS_SAT_SIGMA1" ) != nullptr || oc._eff_sat_sigma1() == 1.0 );
    unsetenv( "CRONOS_AUTO_ALG_CLOSURE" );
    check( "O9 ALG_CLOSURE defaults on",                   FFModel::_knob_AUTO_ALG_CLOSURE() );
    setenv( "CRONOS_AUTO_ALG_CLOSURE", "0", 1 );
    check( "O9 ...and CRONOS_AUTO_ALG_CLOSURE=0 turns it off", !FFModel::_knob_AUTO_ALG_CLOSURE() );
    unsetenv( "CRONOS_AUTO_ALG_CLOSURE" );
  }

  // ---- O10 rev318 --------------------------------------------------------------------------------------------
  std::cout << "\n---- O10: HYP_CLOSURE is environment-only (rev318) ----\n";
  {
    unsetenv( "CRONOS_AUTO_HYP_CLOSURE" );
    check( "O10 HYP_CLOSURE defaults on",                        FFModel::_knob_AUTO_HYP_CLOSURE() );
    setenv( "CRONOS_AUTO_HYP_CLOSURE", "0", 1 );
    check( "O10 ...and CRONOS_AUTO_HYP_CLOSURE=0 turns it off",  !FFModel::_knob_AUTO_HYP_CLOSURE() );
    unsetenv( "CRONOS_AUTO_HYP_CLOSURE" );
  }

  // ---- E1 scalar and vector add_equation --------------------------------------------------------------------
  std::cout << "\n---- E1: scalar AND vector add_equation store every equation ----\n";
  {
    FFGraph D1; OCFESLV s( &D1 ); Vars V1; build( D1, s, V1, false );
    FFGraph D2; OCFESLV v( &D2 ); Vars V2; build( D2, v, V2, true );
    check( "E1 scalar form stores 4 equations", s.var_equation().size() == 4,
           "n=" + std::to_string( s.var_equation().size() ) );
    check( "E1 vector form stores the same 4", v.var_equation().size() == 4,
           "n=" + std::to_string( v.var_equation().size() ) );
  }

  // ---- E2/E3 role-derived defaults, read back by pointer -----------------------------------------------------
  std::cout << "\n---- E2/E3: role-derived defaults, read back through var_equation() ----\n";
  {
    FFGraph DAG; OCFESLV oc( &DAG ); Vars V; build( DAG, oc, V );
    auto const& eq = oc.var_equation();
    bool all_ptr = true;
    for( auto const& e : eq ) if( !e.opt ) all_ptr = false;
    check( "E3 every stored equation holds non-null options", all_ptr );
    if( eq.size() == 4 && all_ptr ){
      check( "E2 INTERIOR participates in classification", eq[0].opt->participate_in_classification );
      check( "E2 BOUNDARY does not participate",            !eq[2].opt->participate_in_classification );
      check( "E2 BOUNDARY is donated for state continuity", OCFESLV::_sopt( *eq[2].opt ).donate_for_state_continuity );
      check( "E2 INTERIOR receives SAT",                    OCFESLV::_sopt( *eq[0].opt ).receive_sat );
      check( "E3 block_id reads back as supplied",          eq[0].opt->block_id == 0 );
    }
  }

  // ---- E4 copies OWN their equation options ------------------------------------------------------------
  std::cout << "\n---- E4: a copy owns its equation options (no shared state with the source) ----\n";
  {
    FFGraph DAG; OCFESLV src( &DAG ); Vars V; build( DAG, src, V );
    configure( src, OCFESLV::Options::IC_TRACE );
    bool const ok = src.setup();
    check( "E4 source setup() succeeds", ok );
    OCFESLV cp( src );
    FFModel const& fs = src; FFModel const& fc = cp;
    // the WORKING equations (the ones _deep_copy_model copies), reached through the model's table
    auto const& es = fs.var_equation(); auto const& ec = fc.var_equation();
    bool distinct = ( es.size() == ec.size() && !es.empty() );
    for( size_t k = 0; distinct && k < es.size(); ++k )
      if( !es[k].opt || !ec[k].opt || es[k].opt.get() == ec[k].opt.get() ) distinct = false;
    check( "E4 source and copy hold DISTINCT options objects", distinct );
    if( distinct ){
      int const before = es[0].opt->block_id;
      const_cast<FFModel::EqnOptions&>( *ec[0].opt ).block_id = before + 7;   // mutate the COPY in place
      check( "E4 mutating the copy's options leaves the source's unchanged", es[0].opt->block_id == before,
             "source block_id=" + std::to_string( es[0].opt->block_id ) );
    }
  }

  // ---- E5 model-minted rows carry the solver's type -------------------------------------------------------
  std::cout << "\n---- E5: equations minted by the model itself carry the solver's options type ----\n";
  {
    FFGraph DAG; OCFESLV oc( &DAG ); Vars V; build( DAG, oc, V );
    configure( oc, OCFESLV::Options::IC_TRACE );        // RED_FULL: reduce_order mints LINK rows for u''
    bool const ok = oc.setup();
    FFModel const& fm = oc;
    size_t nlink = 0, nall = 0, nsolver = 0;
    for( auto const& e : fm.var_equation() ){
      ++nall;
      if( e.opt && dynamic_cast<OCFESLV::EqnOptions const*>( e.opt.get() ) ) ++nsolver;
      if( e.opt && e.opt->role == OCFESLV::EqnRole::LINK ) ++nlink;
    }
    check( "E5 setup() succeeds", ok );
    check( "E5 reduce_order minted at least one LINK row", nlink > 0, "nlink=" + std::to_string( nlink ) );
    check( "E5 EVERY stored equation carries OCFESLV::EqnOptions", nall > 0 && nsolver == nall,
           std::to_string( nsolver ) + "/" + std::to_string( nall ) );
  }

  // ---- E6 a copy keeps the solver's type and values ------------------------------------------------------
  std::cout << "\n---- E6: a copy keeps the solver's options type and field values (clone override) ----\n";
  {
    FFGraph DAG; OCFESLV src( &DAG ); Vars V; build( DAG, src, V );
    configure( src, OCFESLV::Options::IC_TRACE );
    src.setup();
    OCFESLV cp( src );
    FFModel const& fs = src; FFModel const& fc = cp;
    auto const& es = fs.var_equation(); auto const& ec = fc.var_equation();
    bool typed = ( es.size() == ec.size() && !es.empty() ), same = typed;
    for( size_t k = 0; typed && k < ec.size(); ++k ){
      auto const* a = dynamic_cast<OCFESLV::EqnOptions const*>( es[k].opt.get() );
      auto const* b = dynamic_cast<OCFESLV::EqnOptions const*>( ec[k].opt.get() );
      if( !a || !b ){ typed = false; break; }
      if( a->receive_sat != b->receive_sat || a->donate_for_state_continuity != b->donate_for_state_continuity
       || a->interface_type != b->interface_type ) same = false;
    }
    check( "E6 every copied equation still carries OCFESLV::EqnOptions", typed );
    check( "E6 ...with the source's solver-field values", typed && same );
  }

  std::cout << "\n  " << g_pass << " passed, " << g_fail << " failed\n"
            << "OCFE_options: " << ( g_fail ? "FAIL" : "PASS" ) << "\n";
  return g_fail ? 1 : 0;
}
