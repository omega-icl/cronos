// OCFE_DAE5sd.cpp
// ===========================================================================
//  DAE5 IN THE ROW-SPACE VARIABLES OF A_e -- does squaring the core remove k?
// ===========================================================================
//
//  WHAT THIS DRIVER IS FOR
//  -----------------------
//  OCFE_DAE5 is the same physical problem written on {x1,x2}:
//
//      SUM  :  dx1/dt + dx2/dt = -(x1 + x2)     ONE row, TWO derivatives
//      SPLIT:  x1 - x2        = D0 exp(-2t)     algebraic
//
//  so A_e = [1 1], a 1x2 RECTANGULAR differential core.  Measured consequences
//  (rev119, IC_STRONG, --nomarch):
//
//      m = 2*(nel_t - 1) continuity claims -- one per STATE per interior interface
//      k = nel_t - 1     = 3 / 7 / 15 at nel_t = 4 / 8 / 16
//
//  k is exactly one per interface because the principal-symbol coupling of BOTH
//  claims at an interface is A_d = [1 1]: the two columns are COLLINEAR, so B
//  cannot carry an independent multiplier for each.  Downstream of that: dead
//  multiplier columns, a numerically singular Js, tier-0 Newton refused, LM
//  stalling, and a least-squares floor of 1.94e-07 / 6.94e-09 / 2.32e-10 against
//  a fixed SOLVE_RES_TOL of 1e-09 -- so conv=no at nel_t 4 and 8, and DAE5 FAILs
//  on gates (a) and (c) while its ANSWERS are correct at every mesh.
//
//  THE REFORMULATION
//  -----------------
//  A_e = [1 1].  Its ROW SPACE is spanned by (1,1), so the combination the row
//  actually differentiates is
//
//      s = x1 + x2         (differential)
//      d = x1 - x2         (the complement, pinned by SPLIT)
//
//  and the model becomes
//
//      SUMs  :  ds/dt = -s              1x1, SQUARE
//      SPLITd:  d     = D0 exp(-2t)     algebraic
//
//  NOTE what is NOT needed: no constraint is differentiated.  SUM already
//  contains d(x1+x2)/dt = ds/dt, so this is a CHANGE OF VARIABLES, not index
//  reduction -- no drift, no extra initial condition, no differentiation of an
//  algebraic row.
//
//  WHY THE ROW SPACE AND NOT "eliminate x2"
//  ----------------------------------------
//  Substituting x2 = x1 - D0 exp(-2t) into SUM also gives a square 1x1 system,
//  but it requires differentiating SPLIT, and -- more importantly -- the CHOICE
//  of which state to eliminate is arbitrary.  A standalone replica of this
//  kernel measured a factor-6 accuracy spread across such choices
//  (1.13e-07 for the symbol-row combination, 1.33e-07 for keeping both,
//  6.66e-07 and 4.93e-07 for eliminating x1 or x2 respectively).  The row space
//  of A_e is basis-independent: there is nothing to choose.
//
//  WHAT THIS DRIVER MEASURES
//  -------------------------
//  With ONE differential state there is one claim per interior interface and
//  nothing for it to be collinear with, so the PREDICTION is
//
//      m = nel_t - 1,  k = 0,  no dead columns, no floor, Newton at every mesh
//
//  and accuracy NO WORSE than DAE5's 2.56e-07 / 2.07e-08 / 1.51e-09.
//
//  The error is measured on x1 and x2 -- reconstructed as (s+d)/2 and (s-d)/2 --
//  NOT on s and d, so the number is directly comparable with OCFE_DAE5's.
//  Measuring s and d instead would be measuring a different quantity and would
//  make the comparison meaningless.
//
//  HOW TO READ THE RESULT
//  ----------------------
//    k = 0 and conv=yes at every mesh, err <= DAE5's
//        -> squaring the core removes the kernel AT SOURCE.  The rectangular
//           shape was the cause, not the interface plan.
//    k = 0 but err noticeably WORSE than DAE5's
//        -> the change of variables costs accuracy; the over-determined
//           least-squares solution on {x1,x2} was using information this
//           formulation discards.  That would need understanding before any
//           automation.
//    k > 0
//        -> the collinearity is not what the rectangular core produces, and the
//           whole account above is wrong.
//
//  SCOPE, stated plainly: DAE5 is the ONLY driver in the 139-driver corpus with a
//  rectangular core (MMPDE27 is 2x2, PDE9 and bisect2 are 1x1, all others report
//  NOT).  This validates the idea on the founding problem and on nothing else.
//
//  USAGE
//    ./OCFE_DAE5sd [--nomarch] [--weak|--trace] [--tol R] [-v] [nel_t ...]
//
// ===========================================================================

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <limits>
#include <sstream>

#include OCFE_OCFESLV_HEADER

using mc::FFGraph;
using mc::FFVar;
using mc::FFPartial;
using mc::OCFESLV;
using mc::FFDom;

static double const TF = 1.0;
static double const S0 = 2.0;      // s(0) = x1(0) + x2(0)
static double const D0 = 0.5;      // d(0) = x1(0) - x2(0)

static size_t NND_T  = 5;
static double RESTOL = 1.0e-9;
static bool   g_trace = false, g_weak = false, NOMARCH = false;
static int    g_disp  = 0;

static double const XTOL   = 1.0e-8;
static double const ORDTOL = 3.0;

static double s_exact ( double t ){ return S0 * std::exp( -t ); }
static double d_exact ( double t ){ return D0 * std::exp( -2.0 * t ); }
static double x1_exact( double t ){ return 0.5 * ( s_exact(t) + d_exact(t) ); }
static double x2_exact( double t ){ return 0.5 * ( s_exact(t) - d_exact(t) ); }

struct Result
{
  bool   setup_ok = false, square = false, conv = false, threw = false;
  double err = std::numeric_limits<double>::infinity();
  size_t nVar = 0;
};

// ---------------------------------------------------------------------------
static Result run( size_t nel_t, int display )
{
  Result R;

  FFGraph DAG;
  FFVar t = DAG.add_var( "t" );
  FFVar s = DAG.add_var( "s(t)" );      // x1 + x2 -- the row space of A_e = [1 1]
  FFVar d = DAG.add_var( "d(t)" );      // x1 - x2 -- the complement

  FFPartial OpP;   // an OBJECT, not a free function (see OCFE_DAE5.cpp)

  // SUMs differentiates exactly ONE state, so A_e is 1x1 and SQUARE.  This is the
  // whole difference from OCFE_DAE5; everything else below is unchanged.
  FFVar SUMs   = OpP( s, t ) + s;

  // SPLITd pins the complement.  Algebraic -> routed to vAlgEqn, as SPLIT was.
  FFVar SPLITd = d - D0 * exp( -2.0 * t );

  // ONE initial condition, on the ONE differential state.  OCFE_DAE5 needs a note
  // here about which of x1/x2 to close because SPLIT ties them; with s and d
  // separated there is no such choice -- d carries no derivative and needs no IC.
  FFVar IC_S = s - S0;

  OCFESLV oc( &DAG );
  oc.options.DISPLAY_LEVEL = display;
  oc.add_domain( t, FFDom( 0., TF, nel_t, FFDom::LGR, NND_T ) );
  oc.add_state ( s, { t } );
  oc.add_state ( d, { t } );

  oc.reset_evolution_domain();

  OCFESLV::EqnOptions interior( OCFESLV::EqnRole::INTERIOR, 0 );
  OCFESLV::EqnOptions initial ( OCFESLV::EqnRole::INITIAL,  0 );

  int const T_NO_LB = FFDom::ALL - FFDom::LB;

  oc.add_equation( SUMs,   { t }, { T_NO_LB    }, interior );
  oc.add_equation( SPLITd, { t }, { FFDom::ALL }, interior );
  oc.add_equation( IC_S,   { t }, { FFDom::LB  }, initial  );

  oc.options.REDUCE.ORDER    = OCFESLV::Options::RED_NONE;
  oc.options.INTERFACE.TYPE  = OCFESLV::Options::IC_AUTO;
  oc.options.INTERFACE.IMPOSITION = g_weak  ? OCFESLV::Options::IC_WEAK
                             : g_trace ? OCFESLV::Options::IC_TRACE
                                       : OCFESLV::Options::IC_STRONG;
  oc.options.SOLVE.RES_TOL   = RESTOL;
  oc.options.CLASSIFY.MODE        = OCFESLV::Options::CLASS_AUTO;
  if( NOMARCH ) oc.options.SOLVE.MARCHING = false;

  try { R.setup_ok = oc.setup(); }
  catch( ... ) { R.threw = true; return R; }
  if( !R.setup_ok ) return R;

  size_t nEqn = 0;
  try { R.nVar = oc.n_colloc_sta(); nEqn = oc.n_colloc_eqn(); } catch( ... ) {}
  R.square = ( R.nVar && R.nVar == nEqn );

  std::vector<double> varInit, inpInit;
  if( !oc.init( varInit, inpInit, nullptr ) ) return R;
  std::vector<double> xv = varInit;

  OCFESLV::SolveReport rep = oc.solve( xv.data(), nullptr, nullptr );
  R.conv = rep.converged;

  {
    // Reconstruct x1, x2 from (s,d) and measure THEIR error, so the number is
    // directly comparable with OCFE_DAE5's.  Same 41 sample points.
    double emax = 0.;
    size_t const NS = 41;
    for( size_t k = 0; k <= NS; ++k ){
      double const tt = TF * double( k ) / double( NS );
      OCFESLV::t_Coord pt; pt[t] = tt;
      double sv = 0., dv = 0.;
      try { sv = oc.eval_colloc<double>( s, pt, xv.data(), nullptr, nullptr );
            dv = oc.eval_colloc<double>( d, pt, xv.data(), nullptr, nullptr ); }
      catch( ... ) { continue; }
      double const a = 0.5 * ( sv + dv );      // x1
      double const b = 0.5 * ( sv - dv );      // x2
      emax = std::max( emax, std::max( std::fabs( a - x1_exact(tt) ),
                                       std::fabs( b - x2_exact(tt) ) ) );
    }
    R.err = emax;
  }

  return R;
}

// ---------------------------------------------------------------------------
int main( int argc, char* argv[] )
{
  std::vector<size_t> nelList;

  for( int i = 1; i < argc; ++i ){
    if( !std::strcmp( argv[i], "--trace"   ) ){ g_trace = true; continue; }
    if( !std::strcmp( argv[i], "--weak"    ) ){ g_weak  = true; continue; }
    if( !std::strcmp( argv[i], "--nomarch" ) ){ NOMARCH = true; continue; }
    if( !std::strcmp( argv[i], "-v"        ) ){ g_disp  = 1;    continue; }
    if( !std::strcmp( argv[i], "--verbose" ) ){ g_disp  = 1;    continue; }
    if( !std::strcmp( argv[i], "--tol" ) && i+1 < argc ){ RESTOL = std::atof( argv[++i] ); continue; }
    if( argv[i][0] != '-' ){
      double const v = std::atof( argv[i] );
      if( v > 0. ) nelList.push_back( (size_t)v );
      continue;
    }
  }
  if( nelList.empty() ) nelList = { 4, 8, 16 };

  std::cout
    << "================================================================\n"
    << "  DAE5sd -- the SAME problem on the ROW SPACE of A_e\n"
    << "    SUMs  :  ds/dt = -s                    ONE row, ONE derivative\n"
    << "    SPLITd:  d     = " << D0 << " exp(-2t)             algebraic\n"
    << "  s = x1 + x2,  d = x1 - x2;  error reported on x1,x2 = (s+d)/2, (s-d)/2\n"
    << "  exact: s = " << S0 << " e^-t,  d = " << D0 << " e^-2t\n"
    << "  mode=" << ( g_weak ? "IC_WEAK" : g_trace ? "IC_TRACE" : "IC_STRONG" )
    << ( NOMARCH ? "  --nomarch (SOLVE_MARCHING=false: no collapse)" : "" )
    << "  tol=" << RESTOL << "\n"
    << "  PREDICTION vs OCFE_DAE5: k = 0 (was nel_t-1), conv=yes at EVERY mesh,\n"
    << "  err <= 2.56e-07 / 2.07e-08 / 1.51e-09.  Read [promote-scan] k= and the\n"
    << "  [classify] line: expect sym=1x1, not 1x2.\n"
    << "================================================================\n"
    << "  nel_t  nVar   square  conv   max|x-x*|   order   status\n";

  bool all_ok  = true;
  bool fine_ok = false;
  double prev_h = 0., prev_e = 0., worst_ord = 1e30;
  size_t nord  = 0;

  for( size_t n : nelList ){
    Result const R = run( n, g_disp );
    bool const solved = R.setup_ok && R.square && R.conv && std::isfinite( R.err );
    all_ok  = all_ok && solved;
    fine_ok = solved && R.err < XTOL;

    double ord = 0.; bool have_ord = false;
    double const h = TF / double( n );
    if( solved && prev_h > 0. && prev_e > 0. && R.err > 0. ){
      ord = std::log( prev_e / R.err ) / std::log( prev_h / h );
      have_ord = true;
      worst_ord = std::min( worst_ord, ord );
      ++nord;
    }
    if( solved ){ prev_h = h; prev_e = R.err; }

    std::cout << "  " << std::left << std::setw(5) << n
              << std::setw(7) << R.nVar
              << std::setw(8) << ( R.square ? "yes" : "no" )
              << std::setw(7) << ( R.conv   ? "yes" : "no" )
              << std::right << std::scientific << std::setprecision(2)
              << std::setw(11) << R.err << "   "
              << std::left << std::setw(7)
              << ( have_ord ? ( std::ostringstream() << std::fixed
                                << std::setprecision(2) << ord ).str() : std::string("-") )
              << ( R.threw ? "THREW" : solved ? "solved" : "[--]" ) << "\n";
  }

  bool const ord_ok = ( nord == 0 ) || ( worst_ord >= ORDTOL );

  std::cout
    << "  ----------------------------------------------------------------\n"
    << "    (a) every mesh solved            : " << ( all_ok  ? "yes" : "NO" ) << "\n"
    << "    (b) finest mesh |x-x*| < " << XTOL << " : " << ( fine_ok ? "yes" : "NO" ) << "\n"
    << "    (c) worst observed order >= " << ORDTOL << "  : " << ( ord_ok ? "yes" : "NO" );
  if( nord ) std::cout << "  (worst " << std::fixed << std::setprecision(2) << worst_ord
                       << " over " << nord << " refinement(s))";
  std::cout << "\n================================================================\n"
            << "  DAE5sd: " << ( ( all_ok && fine_ok && ord_ok ) ? "PASS" : "FAIL" ) << "\n"
            << "  Compare with ./OCFE_DAE5 --nomarch on the SAME meshes.  The point of\n"
            << "  this driver is the DIFFERENCE, not the verdict: if DAE5sd passes where\n"
            << "  DAE5 fails, with k=0 and no worse error, the rectangular core was the\n"
            << "  cause and the fix belongs at the model level, not in the claim logic.\n"
            << "================================================================\n";

  return ( all_ok && fine_ok && ord_ok ) ? 0 : 1;
}
