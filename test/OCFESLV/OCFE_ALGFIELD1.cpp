// OCFE_ALGFIELD1.cpp
// ===========================================================================
//  A DISTRIBUTED, PURELY ALGEBRAIC FIELD -- NO DERIVATIVE OPERATOR ANYWHERE
// ===========================================================================
//
//  WHY THIS DRIVER EXISTS
//  ----------------------
//  The corpus covers the LUMPED algebraic case thoroughly -- OCFE_AE0/AE1/AE2 are
//  zero-domain, add_state(x,{}), no collocation at all.  The DISTRIBUTED algebraic
//  case is not covered: a state y(t,xi) living on the collocation mesh, defined
//  POINTWISE, with no partial derivative of any state appearing in any equation.
//
//  THIS IS THE SHARPEST AVAILABLE TEST OF THE INTERFACE MACHINERY.
//
//  A distributed algebraic field has NO REASON to need a continuity claim.  The
//  value at each duplicate node is fixed pointwise by the SAME equation on both
//  sides of a seam, so continuity is IMPLIED BY THE PHYSICAL ROWS -- exactly the
//  case rev66's own tau-block note describes:
//
//      "when the continuity rows are IMPLIED by the physical rows the conditions
//       hold identically, the kernel is a gauge freedom and the primal answer is
//       correct (measured: PDE2 ... max|r|=3e-10)"
//
//  So the outcome is binary, and needs no oracle beyond y itself, which is known
//  in closed form:
//
//    ALL THREE MODES AGREE, claims suppressed or harmless
//        -> the redundancy detection works and the machinery has no quarrel with
//           algebraic distributed fields.  That EXONERATES the interface plan for
//           this class and points the MESHONLY/MMPDE investigation back at the
//           interaction with DIFFERENTIATION.
//
//    IC_TRACE / IC_STRONG DIFFER FROM IC_WEAK
//        -> a defect with NO POSSIBLE PHYSICAL EXCUSE.  There is nothing here that
//           can legitimately be discontinuous, no front, no layer, no under-resolved
//           feature, no moving mesh, no ALE term, no lifting alias, no LINK.  Every
//           confound that has dogged this arc is absent by construction.
//
//  WHAT IT ISOLATES THAT MESHONLY CANNOT
//  -------------------------------------
//  MESHONLY still carries x_t and THREE spatial derivative aliases (xg = x_xi,
//  gg = g_xi, and g = M*xg feeding both).  Its exact modes converge to a TANGLED
//  mesh and its residual anomaly sits on the g/gg pair at xi-seams.  Here there is
//  no derivative at all, so if the exact modes still misbehave the problem is in
//  the interface plan ALONE and not in any interaction with differentiation.
//
//  MODEL
//  -----
//    y(t,xi) distributed on (t,xi), and ONE equation, registered on ALL nodes:
//
//        ALG :  y^3 + y - h(t,xi) = 0 ,     h = ( s^3 + s ),  s = sin(2 pi xi)(1+t)
//
//    so the exact solution is y*(t,xi) = sin(2 pi xi) * (1+t), in closed form, and
//    the cubic keeps Newton meaningful rather than making the system linear.
//    y^3 + y is strictly monotone in y, so the root is UNIQUE for every (t,xi) --
//    there is no second branch for an exact mode to fall into, which removes the
//    multiple-root explanation that MESHONLY's tangled field left open.
//
//    NO initial condition and NO boundary condition are registered: an algebraic
//    field has no evolution to initialise and no characteristic to feed.  If setup()
//    demands one, that is itself a finding and is reported rather than worked around.
//
//  USAGE
//    ./OCFE_ALGFIELD1 [--weak|--trace] [--evolve] [--tol R] [--nelt N] [--verbose]
//                     [nel_xi ...]
//    Every bare numeric is a nel_xi.  There is deliberately no Pe-style rule here:
//    MMPDE26's "first bare numeric >= 10 is Pe" caused three separate confusions in
//    one session.
// ===========================================================================

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <limits>

#ifndef OCFE_OCFESLV_HEADER
  #define OCFE_OCFESLV_HEADER "ocfeslv.hpp"
#endif
#include OCFE_OCFESLV_HEADER

using namespace mc;

static double const TF  = 0.5;
static double const PI2 = 6.283185307179586;

static size_t NND_XI = 8, NEL_T = 8, NND_T = 5;
static double RESTOL = 1.0e-9;
static double SIGMA0 = 1.0;
static bool   g_trace = false, g_weak = false, EVOLVE = false;
static int    g_disp  = 0;

// GATE (corrected).  The first version used YTOL = 1e-9 on the reasoning that "a
// pointwise algebraic root should be exact to round-off".  THAT WAS WRONG: y is
// represented in a degree-(NND_XI-1) polynomial basis per element, so the SOLVED field
// carries interpolation error exactly like any other collocated state -- the algebraic
// closure is satisfied AT the nodes (max|r| ~ 1.9e-11 in every run) while y between them
// is a polynomial approximation of a sinusoid.
//
// MEASURED: |y-y*| = 4.37e-08, 2.26e-10, 5.80e-12 at nel_xi = 4, 8, 16 -- observed order
// 7.6 then 5.3, i.e. spectral, as it should be for a smooth field.  Only the coarsest
// mesh failed the old gate, and it failed for the right reason.
//
// The absolute bound below is a SANITY check sized to pass the coarsest mesh in the
// default sweep.  The REAL gate is the convergence rate, checked across the sweep in
// main(): a correct spectral discretisation must cut the error by a large factor per
// refinement, and an absolute threshold cannot distinguish "accurate" from "accurate
// because the mesh happened to be fine enough".
static double const YTOL      = 1.0e-6;   // absolute sanity bound
static double const YRATE_MIN = 4.0;      // minimum error reduction per mesh doubling

// exact solution and the forcing built from it
static double y_exact( double tt, double qq )
  { return std::sin( PI2 * qq ) * ( 1.0 + tt ); }

struct Result
{
  bool   setup_ok = false, square = false, conv = false, threw = false;
  double errY = std::numeric_limits<double>::infinity();
  double spread = -1.0;
  size_t nVar = 0, nTau = 0;
};

// ---------------------------------------------------------------------------
static Result run( size_t nel_xi, int display )
{
  Result R;

  FFGraph DAG;
  FFVar t  = DAG.add_var( "t" );
  FFVar xi = DAG.add_var( "xi" );
  FFVar y  = DAG.add_var( "y(t,xi)" );

  // h(t,xi) = s^3 + s with s = sin(2 pi xi)(1+t), written in the DAG so the
  // forcing is exactly the value that makes y* the root.
  FFVar const sv  = sin( PI2 * xi ) * ( 1.0 + t );
  FFVar const hv  = sv * sv * sv + sv;
  FFVar ALG = y * y * y + y - hv;          // POINTWISE.  No OpP anywhere.

  OCFESLV oc( &DAG );
  oc.options.DISPLAY_LEVEL = display;
  oc.add_domain( t,  FFDom( 0., TF, NEL_T,  FFDom::LGR, NND_T  ) );
  oc.add_domain( xi, FFDom( 0., 1.0, nel_xi, FFDom::LGL, NND_XI ) );
  oc.add_state ( y, { t, xi } );

  if( EVOLVE ) oc.set_evolution_domain( t );
  else         oc.reset_evolution_domain();

  OCFESLV::EqnOptions int_opt( OCFESLV::EqnRole::INTERIOR, 0 );

  // Registered on ALL nodes in BOTH domains: an algebraic field is determined
  // everywhere by the same closure, including at every element face.  No IC, no BC.
  oc.add_equation( ALG, { t, xi }, { FFDom::ALL, FFDom::ALL }, int_opt );

  oc.options.REDUCE.ORDER    = OCFESLV::Options::RED_MAIN;
  oc.options.INTERFACE.TYPE  = OCFESLV::Options::IC_AUTO;
  oc.options.INTERFACE.IMPOSITION = g_weak  ? OCFESLV::Options::IC_WEAK
                             : g_trace ? OCFESLV::Options::IC_TRACE
                                       : OCFESLV::Options::IC_STRONG;
  oc.options.INTERFACE.SAT_SIGMA0      = SIGMA0;
  oc.options.SOLVE.RES_TOL   = RESTOL;
  oc.options.CLASSIFY.MODE        = OCFESLV::Options::CLASS_AUTO;

  try { R.setup_ok = oc.setup(); }
  catch( ... ) { R.threw = true; return R; }
  if( !R.setup_ok ) return R;

  size_t nEqn = 0;
  try { R.nVar = oc.n_colloc_sta(); R.nTau = oc.n_colloc_trace();
        nEqn = oc.n_colloc_eqn(); } catch( ... ) {}
  R.square = ( R.nVar && R.nVar == nEqn );

  std::vector<double> varInit, inpInit;
  if( !oc.init( varInit, inpInit, nullptr ) ) return R;
  std::vector<double> xv = varInit;

  OCFESLV::SolveReport rep = oc.solve( xv.data(), nullptr, nullptr );
  R.conv = rep.converged;

  // |y - y*| over a sampling lattice, and the duplicate-node spread at xi-seams.
  // The spread is the DIRECT question: at a seam the two copies of y solve the SAME
  // scalar cubic with the SAME forcing, so any difference between them is produced
  // entirely by the imposition machinery.
  {
    double emax = 0., smax = 0.;
    size_t const NT_S = 11, PPE = 6;
    for( size_t k = 0; k <= NT_S; ++k ){
      double const tt = TF * double( k ) / double( NT_S );
      for( size_t e = 0; e < nel_xi; ++e )
        for( size_t sI = 0; sI <= PPE; ++sI ){
          double const qq = ( double( e ) + double( sI ) / double( PPE ) ) / double( nel_xi );
          OCFESLV::t_Coord pt; pt[t] = tt; pt[xi] = qq;
          double yn = 0.;
          try { yn = oc.eval_colloc<double>( y, pt, xv.data(), nullptr, nullptr ); }
          catch( ... ) { continue; }
          emax = std::max( emax, std::fabs( yn - y_exact( tt, qq ) ) );
        }
      // seam spread: evaluate either side of each interior xi-element boundary
      for( size_t e = 1; e < nel_xi; ++e ){
        double const qs = double( e ) / double( nel_xi );
        OCFESLV::t_Coord lo, hi;
        lo[t] = tt; lo[xi] = qs - 1e-12;
        hi[t] = tt; hi[xi] = qs + 1e-12;
        double a = 0., b = 0.;
        try { a = oc.eval_colloc<double>( y, lo, xv.data(), nullptr, nullptr );
              b = oc.eval_colloc<double>( y, hi, xv.data(), nullptr, nullptr ); }
        catch( ... ) { continue; }
        smax = std::max( smax, std::fabs( a - b ) );
      }
    }
    R.errY = emax; R.spread = smax;
  }

  return R;
}

// ---------------------------------------------------------------------------
int main( int argc, char* argv[] )
{
  std::vector<size_t> nelList;

  for( int i = 1; i < argc; ++i ){
    if( !std::strcmp( argv[i], "--small"  ) ){ NND_XI = 5; NEL_T = 2; NND_T = 3; continue; }
    if( !std::strcmp( argv[i], "--trace"  ) ){ g_trace = true; continue; }
    if( !std::strcmp( argv[i], "--weak"   ) ){ g_weak  = true; continue; }
    if( !std::strcmp( argv[i], "--evolve" ) ){ EVOLVE  = true; continue; }
    if( !std::strcmp( argv[i], "--verbose") ){ g_disp  = 1;    continue; }
    if( !std::strcmp( argv[i], "-v"       ) ){ g_disp  = 1;    continue; }
    if( !std::strcmp( argv[i], "--nelt"   ) && i+1 < argc ){ NEL_T  = (size_t)std::atol( argv[++i] ); continue; }
    if( !std::strcmp( argv[i], "--tol"    ) && i+1 < argc ){ RESTOL = std::atof( argv[++i] ); continue; }
    if( !std::strcmp( argv[i], "--sigma0" ) && i+1 < argc ){ SIGMA0 = std::atof( argv[++i] ); continue; }
    if( argv[i][0] != '-' ){
      double const v = std::atof( argv[i] );
      if( v > 0. ) nelList.push_back( (size_t)v );
      continue;
    }
  }
  if( nelList.empty() ) nelList = { 4, 8, 16 };

  std::cout
    << "================================================================\n"
    << "  ALGFIELD1 -- distributed PURELY ALGEBRAIC field, no derivatives\n"
    << "  y^3 + y = h(t,xi),  h built so y* = sin(2 pi xi)(1+t) exactly.\n"
    << "  ONE state, ONE equation on ALL nodes.  No IC, no BC, no OpP.\n"
    << "  y^3+y is strictly monotone, so the pointwise root is UNIQUE --\n"
    << "  there is no second branch for an exact mode to fall into.\n"
    << "  mode=" << ( g_weak ? "IC_WEAK" : g_trace ? "IC_TRACE" : "IC_STRONG" )
    << "   cfg=" << ( EVOLVE ? "evolve" : "monolithic" )
    << "   sigma0=" << SIGMA0 << " tol=" << RESTOL << "\n"
    << "  t-mesh " << NEL_T << "x" << NND_T << "  xi order " << NND_XI << "\n"
    << "  gates: square + conv + |y-y*| < " << YTOL
    << "  AND error ratio >= " << YRATE_MIN << " per refinement\n"
    << "  A seam spread here is produced ENTIRELY by the imposition machinery:\n"
    << "  both copies of y solve the SAME scalar cubic with the SAME forcing.\n"
    << "================================================================\n"
    << "  nel   nVar    nTau   square  conv   |y-y*|      seam spread  status\n";

  bool all_ok = true;
  std::vector<double> errs;
  std::vector<size_t> nels;
  for( size_t n : nelList ){
    Result const R = run( n, g_disp );
    bool const ok = R.setup_ok && R.square && R.conv
                 && std::isfinite( R.errY ) && R.errY < YTOL;
    all_ok = all_ok && ok;
    if( std::isfinite( R.errY ) ){ errs.push_back( R.errY ); nels.push_back( n ); }
    std::cout << "  " << std::setw(4) << n
              << std::setw(8) << R.nVar
              << std::setw(7) << R.nTau
              << std::setw(9) << ( R.square ? "yes" : "NO" )
              << std::setw(7) << ( R.conv ? "yes" : "no" )
              << std::setw(12) << std::scientific << std::setprecision(2) << R.errY
              << std::setw(13) << R.spread
              << std::setw(9) << ( R.threw ? "THREW" : !R.setup_ok ? "SETUP" : ok ? "[OK]" : "[--]" )
              << std::defaultfloat << "\n";
  }

  // CONVERGENCE GATE.  More discriminating than the absolute bound: it fails a solve that
  // is accurate only because the mesh is fine, and it fails a discretisation that has
  // stopped converging even while every individual value looks small.
  if( errs.size() >= 2 ){
    std::cout << "  convergence:\n";
    for( size_t i = 1; i < errs.size(); ++i ){
      double const ratio = errs[i] > 0. ? errs[i-1] / errs[i]
                                        : std::numeric_limits<double>::infinity();
      double const order = ( errs[i] > 0. && nels[i] > nels[i-1] )
                         ? std::log( ratio ) / std::log( double(nels[i]) / double(nels[i-1]) )
                         : std::numeric_limits<double>::infinity();
      bool const rate_ok = !( ratio < YRATE_MIN );
      all_ok = all_ok && rate_ok;
      std::cout << "    nel " << nels[i-1] << " -> " << nels[i]
                << "   ratio=" << std::scientific << std::setprecision(2) << ratio
                << "   order=" << std::fixed << std::setprecision(1) << order
                << ( rate_ok ? "   ok" : "   *** BELOW MINIMUM RATE ***" )
                << std::defaultfloat << "\n";
    }
  }

  std::cout << "================================================================\n"
            << "  ALGFIELD1: " << ( all_ok ? "PASS" : "FAIL" ) << "\n"
            << "  ALL THREE MODES AGREE  -> the interface machinery has no quarrel\n"
            << "    with distributed algebraic fields; the MESHONLY/MMPDE anomaly\n"
            << "    then lies in the interaction with DIFFERENTIATION, not in the\n"
            << "    plan alone.\n"
            << "  EXACT MODES DIFFER     -> a defect with no physical excuse: no\n"
            << "    front, no layer, no moving mesh, no ALE term, no lifting alias,\n"
            << "    no LINK, and a unique pointwise root.  Every confound of this\n"
            << "    arc is absent by construction.\n"
            << "  SETUP failure          -> the framework declines a distributed\n"
            << "    algebraic field; that is itself the answer to the question.\n"
            << "================================================================\n";
  return all_ok ? 0 : 1;
}
