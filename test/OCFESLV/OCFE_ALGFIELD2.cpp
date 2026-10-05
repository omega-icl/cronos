// OCFE_ALGFIELD2.cpp
// ===========================================================================
//  TWO DISTRIBUTED ALGEBRAIC STATES COUPLED BY A STATE x STATE PRODUCT
// ===========================================================================
//
//  WHY THIS DRIVER EXISTS
//  ----------------------
//  ALGFIELD1 established that the framework handles a SINGLE distributed algebraic
//  field correctly: no claims emitted (nTau=0), all three imposition modes
//  bit-identical, |y-y*| converging at order ~7.  That exonerated the interface
//  plan for a state whose closure is pointwise and SCALAR.
//
//  The moving-mesh models are not scalar-pointwise.  They carry PRODUCTS:
//
//      MESHONLY   DEF  :  g  = M(t,xi) * xg        KNOWN coefficient x STATE
//      MMPDE26    QDEF :  q  * x_xi   = u_xi       STATE x STATE  ("multiply form")
//
//  The second is the structure this arc kept circling.  Its conditioning depends
//  entirely on x_xi staying away from zero -- and MESHONLY's exact modes converge
//  to a TANGLED mesh, min(xg) = -4.77e+02, where x_xi CHANGES SIGN.  A product
//  closure p*a = w has a unique solution for p only while a is bounded away from 0;
//  as a -> 0 the map is singular and near it the Jacobian is ill-conditioned.
//
//  This driver reproduces exactly that product, with NO derivative anywhere, so the
//  product and its conditioning can be studied without a mesh, a front, an ALE term,
//  a lifting alias, or a differential operator of any kind.
//
//  MODEL
//  -----
//    two distributed states a(t,xi), p(t,xi); two equations, both on ALL nodes:
//
//      DEF  :  a^3 + a - ( a*^3 + a* )     = 0     <- a alone, strictly monotone
//      LINK :  p * a - ( p* * a* )         = 0     <- STATE x STATE product
//
//    with        s(t,xi) = sin(2 pi xi) * (1+t) / (1+TF)      so |s| <= 1
//                a*      = AMIN + 1 + s                       so a* in [AMIN, AMIN+2]
//                p*      = cos(2 pi xi) * (1+t)
//
//    DEF fixes a uniquely (z^3+z is strictly monotone, so there is no second branch),
//    and LINK then fixes p uniquely PROVIDED a != 0.  --amin sets how close a* comes
//    to zero, and therefore how ill-conditioned the product closure is:
//
//        --amin 1.0   well separated, cond ~ 1        (the untangled analogue)
//        --amin 1e-2  a* dips to 1e-2                 (near-singular)
//        --amin 1e-4  a* dips to 1e-4                 (severely so)
//        --amin 0     a* TOUCHES zero -- LINK is singular there by construction, and
//                     this is the ALGEBRAIC analogue of a tangled mesh (x_xi -> 0)
//
//  WHAT EACH OUTCOME MEANS
//  -----------------------
//    all three modes agree at every amin
//        -> the interface machinery is indifferent to product closures too, and the
//           MESHONLY/MMPDE anomaly needs the DIFFERENTIAL structure, not the product.
//
//    modes agree at amin=1 but diverge as amin -> 0
//        -> the product's conditioning is what the exact modes are sensitive to.  That
//           would connect directly to the tangling: MESHONLY's exact modes converge to
//           min(xg) < 0, and a sign change requires passing through zero.
//
//    modes differ even at amin=1
//        -> a defect in handling a state x state closure per se, with no conditioning
//           excuse available.
//
//  Gates are resolution-aware (see ALGFIELD1's correction): the absolute bound is a
//  sanity check and the REAL gate is the convergence rate across the sweep.  a and p
//  are smooth, so both must fall spectrally.
//
//  USAGE
//    ./OCFE_ALGFIELD2 [--weak|--trace] [--evolve] [--amin A] [--a0 V] [--showinit]
//                     [--tol R] [--nelt N] [-v] [nel_xi ...]
//
//  --showinit  report what init() actually produced for the a and p blocks (ranges,
//              located via pos_state/node_colloc, not assumed)
//  --a0 V      reseed a at V, leaving p at init()'s value.  d(LINK)/dp = a, so an a
//              block containing 0 makes every p column FREE rather than
//              one-pivot-from-free, and the Jacobian is singular at the first iterate.
// ===========================================================================

#include <iostream>
#include <iomanip>
#include <vector>
#include <map>            // ALGFIELD2: pos_state takes std::map<FFVar,size_t,lt_FFVar>
#include <utility>        // std::pair for the block lookup
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
static double AMIN   = 1.0;      // --amin : minimum of a*, i.e. distance from singular
static double A0SEED = std::numeric_limits<double>::quiet_NaN();  // --a0 <v> : seed a at v
static bool   SHOWINIT = false;  // --showinit : report what init() actually produced
static double RESTOL = 1.0e-9;
static double SIGMA0 = 1.0;
static bool   g_trace = false, g_weak = false, EVOLVE = false;
static int    g_disp  = 0;

static double const ATOL      = 1.0e-6;   // absolute sanity bound
static double const ARATE_MIN = 4.0;      // minimum error reduction per refinement

static double s_of  ( double tt, double qq )
  { return std::sin( PI2 * qq ) * ( 1.0 + tt ) / ( 1.0 + TF ); }
static double a_exact( double tt, double qq ) { return AMIN + 1.0 + s_of( tt, qq ); }
static double p_exact( double tt, double qq )
  { return std::cos( PI2 * qq ) * ( 1.0 + tt ); }

struct Result
{
  bool   setup_ok = false, square = false, conv = false, threw = false;
  double errA = std::numeric_limits<double>::infinity();
  double errP = std::numeric_limits<double>::infinity();
  double amin_solved = std::numeric_limits<double>::infinity();
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
  FFVar a  = DAG.add_var( "a(t,xi)" );
  FFVar p  = DAG.add_var( "p(t,xi)" );

  // exact profiles as DAG expressions, so the forcing is exactly consistent
  FFVar const sv  = sin( PI2 * xi ) * ( 1.0 + t ) / ( 1.0 + TF );
  FFVar const av  = AMIN + 1.0 + sv;                       // a*
  FFVar const pv  = cos( PI2 * xi ) * ( 1.0 + t );         // p*

  FFVar DEF  = a * a * a + a - ( av * av * av + av );      // a alone, monotone
  FFVar LINK = p * a - pv * av;                            // STATE x STATE product

  OCFESLV oc( &DAG );
  oc.options.DISPLAY_LEVEL = display;
  oc.add_domain( t,  FFDom( 0., TF,  NEL_T,  FFDom::LGR, NND_T  ) );
  oc.add_domain( xi, FFDom( 0., 1.0, nel_xi, FFDom::LGL, NND_XI ) );
  oc.add_state ( a, { t, xi } );
  oc.add_state ( p, { t, xi } );

  if( EVOLVE ) oc.set_evolution_domain( t );
  else         oc.reset_evolution_domain();

  OCFESLV::EqnOptions int_opt( OCFESLV::EqnRole::INTERIOR, 0 );
  oc.add_equation( DEF,  { t, xi }, { FFDom::ALL, FFDom::ALL }, int_opt );
  oc.add_equation( LINK, { t, xi }, { FFDom::ALL, FFDom::ALL }, int_opt );

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

  // --------------------------------------------------------------------------
  // INITIAL GUESS AND THE MULTIPLY-FORM SINGULARITY
  // --------------------------------------------------------------------------
  // OBSERVED, every run: "tier-0 Newton singular on the augmented interface system;
  // using LM steps", with 1280 of 2560 columns reported singly covered, all p.
  //
  // p appears in exactly ONE equation, so its Jacobian column has ONE entry,
  // d(LINK)/dp = a, and the system is block-triangular:
  //
  //       [ dDEF/da     0      ]
  //       [ dLINK/da  diag(a)  ]
  //
  // Singly-covered is therefore STRUCTURALLY NORMAL for a state defined by a single
  // closure -- "one pivot from free" is accurate, and that pivot is a.  But if a is
  // seeded at ZERO the diag(a) block is EXACTLY zero, every p column is free rather
  // than one-pivot-from-free, and the Jacobian is singular at the first iterate.
  //
  // INFERRED, not yet measured: the reported initial residual 2.9933e+01 matches
  // |DEF| at a=0, which is a*^3 + a* = 3^3 + 3 = 30 at amin=1, to three digits.
  // --showinit prints what init() ACTUALLY produced so the inference is replaced by a
  // measurement -- this arc has had a long run of plausible inferences turn out wrong.
  //
  // --a0 <v> then reseeds a (leaving p at init()'s value).  If tier-0 succeeds at
  // a0=1 and fails at a0=0, the mechanism is confirmed, and it is a general statement
  // about MULTIPLY-FORM closures rather than anything specific to meshes: MMPDE26's
  // q*x_xi = u_xi is singular wherever its multiplicand vanishes -- at a zero initial
  // guess, and at a TANGLED mesh where x_xi crosses zero.
  size_t const n_sta_blk = oc.n_colloc_sta() - oc.n_colloc_trace();
  if( SHOWINIT || std::isfinite( A0SEED ) ){
    // Locate the a and p blocks by measurement, not by assuming a layout.
    auto blockof = [&]( FFVar const& V )->std::pair<size_t,size_t>{
      std::map<FFVar,size_t,lt_FFVar> ndx0;
      auto const& vs = oc.var_state();
      auto const it = vs.find( V );
      if( it != vs.end() ) for( auto const& d : it->second ) ndx0[d] = 0;
      size_t off = 0;
      try { off = oc.pos_state( V, ndx0 ); } catch( ... ) { return { 0, 0 }; }
      return { off, oc.node_colloc( V ).size() };
    };
    auto const ba = blockof( a ), bp = blockof( p );

    if( SHOWINIT ){
      auto rng = [&]( std::pair<size_t,size_t> b )->std::pair<double,double>{
        if( !b.second ) return { 0., 0. };
        double lo = xv[b.first], hi = xv[b.first];
        for( size_t i = 0; i < b.second && b.first + i < xv.size(); ++i ){
          lo = std::min( lo, xv[b.first+i] ); hi = std::max( hi, xv[b.first+i] ); }
        return { lo, hi }; };
      auto const ra = rng( ba ), rp = rng( bp );
      std::cerr << "  [showinit] n_sta=" << n_sta_blk
                << "  a: off=" << ba.first << " n=" << ba.second
                << " range=[" << std::scientific << std::setprecision(3)
                << ra.first << "," << ra.second << "]"
                << "  p: off=" << bp.first << " n=" << bp.second
                << " range=[" << rp.first << "," << rp.second << "]"
                << std::defaultfloat << "\n"
                << "  [showinit] d(LINK)/dp = a, so an a-range containing 0 makes the"
                   " p block of the Jacobian singular at the first iterate.\n";
    }

    if( std::isfinite( A0SEED ) ){
      if( !ba.second )
        std::cerr << "  [a0seed] REFUSED: could not locate the a block; not reseeding.\n";
      else{
        for( size_t i = 0; i < ba.second && ba.first + i < xv.size(); ++i )
          xv[ba.first+i] = A0SEED;
        std::cerr << "  [a0seed] a seeded at " << A0SEED << " over " << ba.second
                  << " coefficient(s); p left at init() value.\n";
      }
    }
  }

  OCFESLV::SolveReport rep = oc.solve( xv.data(), nullptr, nullptr );
  R.conv = rep.converged;

  {
    double ea = 0., ep = 0., amn = std::numeric_limits<double>::infinity(), smax = 0.;
    size_t const NT_S = 11, PPE = 6;
    for( size_t k = 0; k <= NT_S; ++k ){
      double const tt = TF * double( k ) / double( NT_S );
      for( size_t e = 0; e < nel_xi; ++e )
        for( size_t sI = 0; sI <= PPE; ++sI ){
          double const qq = ( double( e ) + double( sI ) / double( PPE ) ) / double( nel_xi );
          OCFESLV::t_Coord pt; pt[t] = tt; pt[xi] = qq;
          double an = 0., pn = 0.;
          try { an = oc.eval_colloc<double>( a, pt, xv.data(), nullptr, nullptr );
                pn = oc.eval_colloc<double>( p, pt, xv.data(), nullptr, nullptr ); }
          catch( ... ) { continue; }
          ea  = std::max( ea, std::fabs( an - a_exact( tt, qq ) ) );
          ep  = std::max( ep, std::fabs( pn - p_exact( tt, qq ) ) );
          amn = std::min( amn, an );          // how close the SOLVED a came to zero
        }
      // seam spread on p: at a seam both copies solve the SAME product closure with the
      // SAME forcing, so any difference is produced entirely by the imposition machinery.
      for( size_t e = 1; e < nel_xi; ++e ){
        double const qs = double( e ) / double( nel_xi );
        OCFESLV::t_Coord lo, hi;
        lo[t] = tt; lo[xi] = qs - 1e-12;
        hi[t] = tt; hi[xi] = qs + 1e-12;
        double u = 0., v = 0.;
        try { u = oc.eval_colloc<double>( p, lo, xv.data(), nullptr, nullptr );
              v = oc.eval_colloc<double>( p, hi, xv.data(), nullptr, nullptr ); }
        catch( ... ) { continue; }
        smax = std::max( smax, std::fabs( u - v ) );
      }
    }
    R.errA = ea; R.errP = ep; R.amin_solved = amn; R.spread = smax;
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
    if( !std::strcmp( argv[i], "-v"       ) ){ g_disp  = 1;    continue; }
    if( !std::strcmp( argv[i], "--verbose") ){ g_disp  = 1;    continue; }
    if( !std::strcmp( argv[i], "--amin"   ) && i+1 < argc ){ AMIN   = std::atof( argv[++i] ); continue; }
    if( !std::strcmp( argv[i], "--a0" ) && i+1 < argc ){ A0SEED = std::atof( argv[++i] ); continue; }
    if( !std::strcmp( argv[i], "--showinit" ) ){ SHOWINIT = true; continue; }
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
    << "  ALGFIELD2 -- TWO distributed algebraic states, STATE x STATE product\n"
    << "    DEF  :  a^3 + a = a*^3 + a*      (a alone, strictly monotone -> unique)\n"
    << "    LINK :  p * a   = p* * a*        (the multiply-form structure)\n"
    << "  a* = " << AMIN << " + 1 + sin(2 pi xi)(1+t)/(1+tf)   so min a* = " << AMIN << "\n"
    << "  p* = cos(2 pi xi)(1+t).   NO derivative operator anywhere.\n"
    << "  mode=" << ( g_weak ? "IC_WEAK" : g_trace ? "IC_TRACE" : "IC_STRONG" )
    << "   cfg=" << ( EVOLVE ? "evolve" : "monolithic" )
    << "   sigma0=" << SIGMA0 << " tol=" << RESTOL << "\n"
    << "  t-mesh " << NEL_T << "x" << NND_T << "  xi order " << NND_XI << "\n"
    << "  gates: square + conv + |a-a*|,|p-p*| < " << ATOL
    << "  AND error ratio >= " << ARATE_MIN << " per refinement\n"
    << "  --showinit reports what init() produced; --a0 <v> reseeds a.  d(LINK)/dp = a,\n"
    << "  so a seeded at 0 makes every p column FREE, not one-pivot-from-free.\n"
    << "  --amin sets how close the product closure comes to SINGULAR: p is determined\n"
    << "  by p*a=w only while a is bounded away from 0.  This is the algebraic analogue\n"
    << "  of a TANGLED mesh -- MESHONLY's exact modes reach min(xg) = -4.77e+02, and a\n"
    << "  sign change must pass through zero.\n"
    << "================================================================\n"
    << "  nel   nVar   nTau  square conv   |a-a*|      |p-p*|      min a     seam sprd  status\n";

  bool all_ok = true;
  std::vector<double> errs; std::vector<size_t> nels;
  for( size_t n : nelList ){
    Result const R = run( n, g_disp );
    double const worst = std::max( R.errA, R.errP );
    bool const ok = R.setup_ok && R.square && R.conv
                 && std::isfinite( worst ) && worst < ATOL;
    all_ok = all_ok && ok;
    if( std::isfinite( worst ) ){ errs.push_back( worst ); nels.push_back( n ); }
    std::cout << "  " << std::setw(4) << n
              << std::setw(7) << R.nVar
              << std::setw(6) << R.nTau
              << std::setw(7) << ( R.square ? "yes" : "NO" )
              << std::setw(6) << ( R.conv ? "yes" : "no" )
              << std::setw(12) << std::scientific << std::setprecision(2) << R.errA
              << std::setw(12) << R.errP
              << std::setw(11) << R.amin_solved
              << std::setw(12) << R.spread
              << std::setw(8) << ( R.threw ? "THREW" : !R.setup_ok ? "SETUP" : ok ? "[OK]" : "[--]" )
              << std::defaultfloat << "\n";
  }

  // Convergence gate -- more discriminating than the absolute bound, which cannot tell
  // "accurate" from "accurate because the mesh happened to be fine enough".
  if( errs.size() >= 2 ){
    std::cout << "  convergence (worst of a,p):\n";
    for( size_t i = 1; i < errs.size(); ++i ){
      double const ratio = errs[i] > 0. ? errs[i-1] / errs[i]
                                        : std::numeric_limits<double>::infinity();
      double const order = ( errs[i] > 0. && nels[i] > nels[i-1] )
                         ? std::log( ratio ) / std::log( double(nels[i]) / double(nels[i-1]) )
                         : std::numeric_limits<double>::infinity();
      bool const rate_ok = !( ratio < ARATE_MIN );
      all_ok = all_ok && rate_ok;
      std::cout << "    nel " << nels[i-1] << " -> " << nels[i]
                << "   ratio=" << std::scientific << std::setprecision(2) << ratio
                << "   order=" << std::fixed << std::setprecision(1) << order
                << ( rate_ok ? "   ok" : "   *** BELOW MINIMUM RATE ***" )
                << std::defaultfloat << "\n";
    }
  }

  std::cout << "================================================================\n"
            << "  ALGFIELD2: " << ( all_ok ? "PASS" : "FAIL" ) << "\n"
            << "  Compare against ALGFIELD1 (one state, scalar closure), where all three\n"
            << "  modes were bit-identical with nTau=0.  Sweep --amin to separate:\n"
            << "    modes agree at every amin      -> the machinery is indifferent to\n"
            << "      product closures; the MESHONLY/MMPDE anomaly needs DIFFERENTIATION.\n"
            << "    modes agree at amin=1 but not\n"
            << "      as amin -> 0                 -> conditioning of the product is what\n"
            << "      the exact modes are sensitive to, which connects directly to the\n"
            << "      tangling (a sign change in x_xi must pass through zero).\n"
            << "    modes differ at amin=1         -> a defect in state x state closures\n"
            << "      per se, with no conditioning excuse available.\n"
            << "================================================================\n";
  return all_ok ? 0 : 1;
}
