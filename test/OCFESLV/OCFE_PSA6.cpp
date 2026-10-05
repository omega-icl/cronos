// OCFE_PSA6.cpp  ---  FLUX-CONSERVATION ARC, STEP 1, on the PSA STEP-INLET oracle
// ===========================================================================
// Base: OCFE_PSA5 (rung-9 PSA, marching, graded t-mesh, manufactured inlet-composition
// step at t* = 2.0, controls b01 + T_ic, outputs Inv1 / Eff1).  PSA5 already establishes
// that the graded march CONVERGES through the step and that forward == adjoint
// sensitivity across it.  PSA6 asks the Step-1 question the arc actually turns on:
//
//     AT EQUAL COST, HOW FAR DOES MESH GRADING ALONE CLOSE THE DISCRETE FLUX BALANCE
//     AS THE FRONT SHARPENS -- AND DOES THE REMAINING DEFECT CONTAMINATE GRADIENTS?
//
// Why PSA5 and not PDE37
// ----------------------
// PDE37 sharpens the feed at t = 0, where the Danckwerts BC meets the IC.  That
// conflates TWO defects: a corner incompatibility (the BC demands flux while the IC
// demands zero) and the strong-form conservation error.  PSA5's step sits at t* = 2,
// INTERIOR to the horizon, so the front is a clean interior feature and the grading
// question is asked in isolation.  PSA5 is also already marching, already on the
// non-uniform FFDom ctor, and already has registered controls -- so the conservation
// defect can be DIFFERENTIATED, which is the thing PDE37 could not do.
//
// NOTE, and it matters for reading the table: PSA5's feed is c_A at t=0 against a ZERO
// initial bed, so there is ALSO an unsmoothed START-UP step at the t=0 corner.  It is
// left in place by default and MEASURED rather than hidden -- the per-window drift
// separates it from the t* front.  SMOOTH_START turns on a startup ramp that is exactly
// IC-compatible (r(0)=0, r'(0)=O(sech^2 4)), which isolates the interior front.
//
// Axes
// ----
//   (1) eps    step width: {0.2, 0.1, 0.05, 0.025, 0.0125}   (PSA5 ships eps = 0.05)
//   (2) G      how many of the ne_tot t-elements are packed into [t*-W, t*+W],
//              W = 4*eps.  G = 0 is the uniform-mesh baseline.  ne_tot is HELD FIXED,
//              so every run in the sweep has the same element count -- equal cost.
//   (3) startup smoothing OFF (default) / ON.
//
// Because W is PROPORTIONAL to eps, the graded columns are a CONSTANT-RESOLUTION family:
// the element width in the zone is 8*eps/G and the LGL first-gap fraction for n_nd = 6
// is (1 - 0.7650553239)/2 = 0.1174725, so
//        dt1/eps = 8*0.1174725/G      ->   G=2: 0.470     G=4: 0.235   (eps-INDEPENDENT)
// whereas the uniform column has dt1 fixed and dt1/eps blowing up as eps falls.
// THE TEST IS THEREFORE FLATNESS IN eps, not a fitted collapse:
//   * flat along G=2 and G=4  -> the defect is PURE RESOLUTION; grading solves it and
//     Step 2 need only target CONSERVATION (flux-form telescoping), not oscillation;
//   * still degrading at pinned dt1/eps -> the defect is the strong-form conservation
//     error itself and option 2 is mandatory, not optional;
//   * mb flat but undershoot not (or vice versa) -> the two pathologies are decoupled
//     and Step 2 must be scoped to whichever survives.
//
// Conservation diagnostic (the point of the driver)
// -------------------------------------------------
// Continuous balance, per species, with Danckwerts in / zero-gradient out:
//     Inv_i(t) - Inv_i(0) = Feed_i(t) - Eff_i(t)
//     Inv_i(t)  = \int_0^1 (c_i + F q_i)(z,t) dz     Eff_i(t) = \int_0^t U c_i(1,s) ds
//     Feed_i(t) = \int_0^t U c_i,feed(s) ds          (analytic integrand, see FeedInt)
// Reported three ways:
//   mb_KPI    -- from the framework's own outputs at T_end (Inv_i + Eff_i - Feed_i).
//                This is the number an OPTIMIZER would see.
//   mb_GL     -- Gauss-Legendre re-quadrature of the collocated solution, element by
//                element, EXACT for the piecewise polynomial (GL-10 is exact to degree
//                19; the solution is degree <= n_nd-1 = 5 per element).  Isolates the
//                conservation defect OF THE DISCRETE SOLUTION -- exactly what flux-form
//                DG telescoping would drive to machine precision.
//   d_k       -- PER-MARCH-WINDOW increment of that defect.  Under marching the windows
//                are chained exactly (terminal -> IC transfer is an identity), so the
//                global defect is the SUM of per-window quadrature defects.  d_k says
//                WHICH window loses the mass: the t=0 corner, the t* front, or all of
//                them uniformly.  A uniform d_k would mean a background quadrature
//                error; a spike at the front window means the front is the problem.
// If mb_KPI and mb_GL disagree, that is a framework-quadrature finding, not a scheme
// finding, and it must be resolved before any DG work is scoped.
//
// Gradient contamination (the question PDE37 could not ask)
// ---------------------------------------------------------
// MB_i = Inv_i + Eff_i - Feed_i, and Feed_i does not depend on the controls, so
//     dMB_i/dp = dInv_i/dp + dEff_i/dp
// is just the SUM OF TWO ROWS of the already-validated reduced Jacobian -- free, and
// exercising the "instrumented and differentiable mass balance" the arc assumes.  If
// |dMB/dp| is at the level of |dF/dp| itself, the conservation defect is control-
// sensitive and WILL corrupt a reduced-space optimizer's gradients; if it is orders
// below, the defect is a benign offset and optimization can proceed ahead of the DG work.
// Run at the sharpest eps for each G only (one extra fsens solve apiece).
//
// Root-identity guard (carried over from the homotopy arc)
// -------------------------------------------------------
// "Converged" is not evidence of correctness on these models.  Each config builds its
// own OCFESLV from the analytic reference -- no warm start is ever carried across the
// grading axis, where the same vector would mean different physical points.  Every
// graded run then reports max_t |c1(1,t) - c1(1,t)|_{G=0} at the SAME eps: O(1e-3) is a
// refinement gain, O(1) is a branch change.  That column also charges the honest PRICE
// of grading (the off-zone elements coarsen when G rises at fixed ne_tot).
//
// Outputs: *_summary.out, *_window.out (per-window drift), *_outlet.out, *.gp
// ===========================================================================

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <array>
#include <cmath>
#include <string>
#include <algorithm>

#ifndef OCFE_OCFESLV_HEADER
#define OCFE_OCFESLV_HEADER "ocfeslv.hpp"
#endif
#include OCFE_OCFESLV_HEADER

using namespace mc;

#ifndef PSA6_OUT_PREFIX
#define PSA6_OUT_PREFIX "OCFE_PSA6"
#endif

// ---------------------------------------------------------------------------
// physical model -- IDENTICAL to OCFE_PSA5 (rung-9 / PDE36c).  Do not retune here:
// this is a NUMERICS diagnostic and must stay directly comparable to the PSA5 log.
// ---------------------------------------------------------------------------
static double const U_VEL = 1.0, D_ax = 0.1, F_ph = 0.5, k_ldf = 2.0;
static double const qs1 = 1.0, qs2 = 1.0, b02_L = 1.0;
static double const beta1 = 2.0, beta2 = 1.0, T0_ref = 1.0;
static double const Cp_e = 1.0, G_cv = 1.0, lam = 0.1, dH1 = 1.0, dH2 = 0.5, hw = 0.5, Tw = 1.0;
static double const T_end = 5.0;
static double const t_step = 2.0;                       // interior step location
static double const cA1 = 0.2, cB1 = 0.5;               // c1 feed: steps UP   at t_step
static double const cA2 = 0.5, cB2 = 0.2;               // c2 feed: steps DOWN at t_step
static double const b01_nom = 3.0;

static double const START_KC = 4.0;                     // startup ramp centre = START_KC*eps0
static double const ZONE_KW  = 4.0;                     // refinement half-width W = ZONE_KW*eps
static size_t const NE_TOT   = 9;                       // t-elements, FIXED across the grading axis
static size_t const NTS      = 501;                     // common outlet-history grid

static inline double tsamp( size_t kk ){ return double(kk)/double(NTS-1)*T_end; }

static inline double q1star_g( double c1g, double c2g, double b01 )
{ return qs1*b01*c1g/( 1.0 + b01*c1g + b02_L*c2g ); }
static inline double q2star_g( double c1g, double c2g, double b01 )
{ return qs2*b02_L*c2g/( 1.0 + b01*c1g + b02_L*c2g ); }

// ---------------------------------------------------------------------------
// feed programs.  The composition step is a centred tanh of width eps.  The optional
// startup ramp is SHIFTED AND NORMALISED so that r(0) = 0 EXACTLY and r(inf) = 1
// EXACTLY -- i.e. it is compatible with the zero IC in value AND (to O(sech^2 4)) in
// derivative, which is what makes it isolate the interior front instead of trading one
// corner artefact for another.  A tanh merely centred at 0 would give r(0) = 1/2
// against a zero bed and would be WORSE than the unsmoothed step.
// ---------------------------------------------------------------------------
static inline double lncosh( double x )
{ x = std::fabs(x); return x + std::log1p( std::exp(-2.0*x) ) - std::log(2.0); }

struct Feed
{
  double eps = 0.05;
  bool   smooth_start = false;
  double eps0 = 0.05;

  double startup( double t ) const
  {
    if( !smooth_start ) return 1.0;
    double const t0 = START_KC*eps0, Th = std::tanh( START_KC );
    return ( std::tanh( (t-t0)/eps0 ) + Th )/( 1.0 + Th );
  }
  double step( double t ) const { return 0.5*( 1.0 + std::tanh( (t-t_step)/eps ) ); }
  double c1( double t ) const { return startup(t)*( cA1 + (cB1-cA1)*step(t) ); }
  double c2( double t ) const { return startup(t)*( cA2 + (cB2-cA2)*step(t) ); }

  //! closed-form \int_0^t c_i,feed  -- valid ONLY without startup smoothing (the
  //! product of two tanh factors has no elementary primitive).  Used as a self-test
  //! of the numeric integrator below, never as the working path.
  double closed_int( int i, double t ) const
  {
    double const A = ( i==1? cA1 : cA2 ), B = ( i==1? cB1-cA1 : cB2-cA2 );
    double const S = 0.5*( t + eps*( lncosh( (t-t_step)/eps ) - lncosh( t_step/eps ) ) );
    return A*t + B*S;
  }
};

//! Composite Gauss-Legendre cumulative integral of an ANALYTIC scalar feed on a panel
//! grid fine enough to resolve every feature (panel <= min(eps,eps0)/4).  The integrand
//! is smooth, so GL-10 per panel is accurate to ~1e-15 -- twelve orders below the ~1e-3
//! signal being measured.  The referee must not inject error into what it is refereeing.
static double const GLX[10] = { -0.9739065285171717, -0.8650633666889845, -0.6794095682990244,
                                -0.4333953941292472, -0.1488743389816312,  0.1488743389816312,
                                 0.4333953941292472,  0.6794095682990244,  0.8650633666889845,
                                 0.9739065285171717 };
static double const GLW[10] = {  0.0666713443086881,  0.1494513491505806,  0.2190863625159820,
                                 0.2692667193099963,  0.2955242247147529,  0.2955242247147529,
                                 0.2692667193099963,  0.2190863625159820,  0.1494513491505806,
                                 0.0666713443086881 };

struct FeedInt
{
  Feed                fd;
  double              h = 0.;
  std::vector<double> cum1, cum2;      // cumulative integrals at panel boundaries

  void build( Feed const& f )
  {
    fd = f;
    double const feat = std::min( f.eps, f.smooth_start? f.eps0 : f.eps );
    size_t np = (size_t)std::ceil( T_end/( feat/4.0 ) );
    np = std::min<size_t>( std::max<size_t>( np, 200 ), 200000 );
    h = T_end/double(np);
    cum1.assign( np+1, 0. ); cum2.assign( np+1, 0. );
    for( size_t k=0; k<np; ++k ){
      double const a=double(k)*h, hm=0.5*h, mid=a+hm;
      double s1=0., s2=0.;
      for( int g=0; g<10; ++g ){
        double const s=mid+hm*GLX[g], w=hm*GLW[g];
        s1 += w*fd.c1(s); s2 += w*fd.c2(s); }
      cum1[k+1]=cum1[k]+s1; cum2[k+1]=cum2[k]+s2;
    }
  }
  //! \int_0^t c_i,feed : whole panels from the table + one exact partial panel
  double operator()( int i, double t ) const
  {
    if( t <= 0. ) return 0.;
    double const tt = std::min( t, T_end );
    size_t k = (size_t)std::floor( tt/h );
    if( k >= cum1.size()-1 ) k = cum1.size()-2;
    double v = ( i==1? cum1[k] : cum2[k] );
    double const a=double(k)*h, b=tt, hm=0.5*(b-a), mid=0.5*(a+b);
    for( int g=0; g<10; ++g ){
      double const s=mid+hm*GLX[g], w=hm*GLW[g];
      v += w*( i==1? fd.c1(s) : fd.c2(s) ); }
    return v;
  }
};

// ---------------------------------------------------------------------------
// three-zone t-mesh: [0,t*-W] | [t*-W,t*+W] with G elements | [t*+W,T_end].
// ne_tot is conserved, so raising G COARSENS the off-zone elements -- that is the price
// of grading at fixed cost, and it is charged to the outdev column, not hidden.
// G even puts an element boundary exactly AT t*, where LGL clusters nodes -- so the
// steepest point of the front sits on the finest part of the stencil.
// ---------------------------------------------------------------------------
static std::vector<double> mesh_step( double L, double ts, size_t ne_tot, size_t G, double W )
{
  std::vector<double> b;
  double const lo = ts - W, hi = ts + W;
  if( G == 0 || G + 2 > ne_tot || lo <= 0. || hi >= L ){       // uniform fallback
    for( size_t i=0; i<=ne_tot; ++i ) b.push_back( L*double(i)/double(ne_tot) );
    return b;
  }
  size_t const nrest = ne_tot - G;
  double const La = lo, Lc = L - hi;
  size_t nA = (size_t)std::llround( double(nrest)*La/(La+Lc) );
  if( nA < 1 ) nA = 1;
  if( nA > nrest-1 ) nA = nrest-1;
  size_t const nC = nrest - nA;
  for( size_t i=0; i<nA; ++i ) b.push_back( La*double(i)/double(nA) );
  for( size_t i=0; i<G;  ++i ) b.push_back( lo + 2.0*W*double(i)/double(G) );
  for( size_t i=0; i<=nC; ++i ) b.push_back( hi + Lc*double(i)/double(nC) );
  b.front() = 0.; b.back() = L;
  return b;
}

//! LGL first-node gap fraction for n nodes (n=6 -> 0.1174725); used for reporting only
static double lgl_gap_frac( size_t n )
{
  if( n == 6 ){ double const x = std::sqrt( ( 7.0 + 2.0*std::sqrt(7.0) )/21.0 ); return 0.5*(1.0-x); }
  double const PI = std::acos(-1.0);                     // CGL-like fallback estimate
  return 0.5*( 1.0 - std::cos( PI/double(n-1) ) );
}

// ---------------------------------------------------------------------------
struct Cfg {
  double eps = 0.05;
  size_t G   = 0;
  bool   smooth_start = false;
  double eps0 = 0.05;
  bool   want_sens = false;
  std::string tag;
};

struct SR {
  Cfg    cfg;
  bool   converged = false;
  int    iters = 0;
  size_t nVar = 0, nwin = 0;
  double dt1 = 0., ratio = 0.;
  double mb_kpi1 = 1., mb_kpi2 = 1.;
  double mb_gl1  = 1., mb_gl2  = 1.;
  double dwin_max1 = 0., dwin_max2 = 0.;      // largest single-window defect
  size_t kwin_max1 = 0;                       // which window that was
  double under1 = 0., under2 = 0., over1 = 0., over2 = 0.;
  double dMB1_b01 = 0., dMB1_inf = 0.;        // gradient contamination (if want_sens)
  double dF1_b01  = 0., dF1_inf  = 0.;        // scale to compare it against
  double outdev = -1.;
  std::vector<double> c1out;                  // c1(1,t) on the common grid
  std::vector<double> wt, wd1, wd2;           // per-window boundary t and defect increments
};

// ---------------------------------------------------------------------------
static SR run_case( Cfg const& cfg )
{
  SR R; R.cfg = cfg;
  size_t const n_el = 5, n_nd = 6;
  FFDom::TYPE const coltype = FFDom::LGL;

  Feed fd; fd.eps = cfg.eps; fd.smooth_start = cfg.smooth_start; fd.eps0 = cfg.eps0;
  FeedInt FI; FI.build( fd );

  double const W = ZONE_KW*cfg.eps;
  std::vector<double> const t_bnd = mesh_step( T_end, t_step, NE_TOT, cfg.G, W );
  size_t kz = 0; for( size_t i=0; i+1<t_bnd.size(); ++i )      // the element starting at t*-W
    if( cfg.G && std::fabs( t_bnd[i] - (t_step-W) ) < 1e-12 ){ kz = i; break; }
  double const hzone = cfg.G ? ( 2.0*W/double(cfg.G) ) : ( T_end/double(NE_TOT) );
  R.dt1   = hzone*lgl_gap_frac( n_nd );
  R.ratio = R.dt1/cfg.eps;

  std::cout << "\n---- eps=" << std::fixed << std::setprecision(5) << cfg.eps
            << "  G=" << cfg.G << "  start=" << (cfg.smooth_start? "smooth":"step")
            << "  dt1=" << std::scientific << std::setprecision(3) << R.dt1
            << "  dt1/eps=" << std::fixed << std::setprecision(3) << R.ratio
            << "  " << cfg.tag << " ----\n     t-mesh:";
  for( double b : t_bnd ) std::cout << " " << std::setprecision(4) << b;
  std::cout << "   (zone starts at element " << kz << ")\n";

  // self-test of the feed integrator against the closed form (no-startup case only)
  if( !cfg.smooth_start ){
    double const e1 = std::fabs( FI(1,T_end) - fd.closed_int(1,T_end) );
    double const e2 = std::fabs( FI(2,T_end) - fd.closed_int(2,T_end) );
    std::cout << "     feed-integrator self-test vs closed form: "
              << std::scientific << std::setprecision(2) << std::max(e1,e2)
              << ( std::max(e1,e2) < 1e-11 ? "  OK\n" : "  ** SUSPECT **\n" );
  }

  // ---------------- model (== OCFE_PSA5) ----------------
  FFGraph DAG;
  FFVar t  = DAG.add_var( "t" );
  FFVar z  = DAG.add_var( "z" );
  FFVar c1 = DAG.add_var( "c1(t,z)" );
  FFVar c2 = DAG.add_var( "c2(t,z)" );
  FFVar q1 = DAG.add_var( "q1(t,z)" );
  FFVar q2 = DAG.add_var( "q2(t,z)" );
  FFVar T  = DAG.add_var( "T(t,z)"  );
  FFVar c1_ic = DAG.add_var( "c1_ic(z)" );
  FFVar c2_ic = DAG.add_var( "c2_ic(z)" );
  FFVar q1_ic = DAG.add_var( "q1_ic(z)" );
  FFVar q2_ic = DAG.add_var( "q2_ic(z)" );
  FFVar T_ic  = DAG.add_var( "T_ic(z)"  );
  FFVar b01   = DAG.add_var( "b01" );

  FFPartial  OpP;
  FFIntegral OpI;

  // DAG-side feed, matched term-for-term to Feed::c1 / Feed::c2
  // Built WITHOUT a ternary over FFVar: a "0.0*t + 1.0" constant arm would leave a
  // degenerate t-dependence in the DAG, and a spurious evolution-variable dependence in
  // a BC is exactly what bit PDE38 ("DOMAIN VARIABLE t USED IN OUTPUT").  Branch instead.
  double const s0 = START_KC*cfg.eps0, Sh = std::tanh( START_KC );
  FFVar sstep  = 0.5*( 1.0 + tanh( ( t - t_step )/cfg.eps ) );
  FFVar c1feed = cA1 + ( cB1 - cA1 )*sstep;
  FFVar c2feed = cA2 + ( cB2 - cA2 )*sstep;
  if( cfg.smooth_start ){
    FFVar rstart = ( tanh( ( t - s0 )/cfg.eps0 ) + Sh )/( 1.0 + Sh );
    c1feed = rstart*c1feed;
    c2feed = rstart*c2feed;
  }

  FFVar b1 = b01  *exp( beta1*( 1.0/T - 1.0/T0_ref ) );
  FFVar b2 = b02_L*exp( beta2*( 1.0/T - 1.0/T0_ref ) );
  FFVar den = 1.0 + b1*c1 + b2*c2;
  FFVar q1star = qs1*b1*c1/den;
  FFVar q2star = qs2*b2*c2/den;

  FFVar CONT1 = OpP( c1, t ) + U_VEL*OpP( c1, z ) - D_ax*OpP( OpP( c1, z ), z ) + F_ph*OpP( q1, t );
  FFVar CONT2 = OpP( c2, t ) + U_VEL*OpP( c2, z ) - D_ax*OpP( OpP( c2, z ), z ) + F_ph*OpP( q2, t );
  FFVar LDF1  = OpP( q1, t ) - k_ldf*( q1star - q1 );
  FFVar LDF2  = OpP( q2, t ) - k_ldf*( q2star - q2 );
  FFVar ENE_T = Cp_e*OpP( T, t ) + G_cv*OpP( T, z ) - lam*OpP( OpP( T, z ), z )
              - F_ph*( dH1*OpP( q1, t ) + dH2*OpP( q2, t ) ) + hw*( T - Tw );

  FFVar IC_c1 = c1 - c1_ic, IC_c2 = c2 - c2_ic;
  FFVar IC_q1 = q1 - q1_ic, IC_q2 = q2 - q2_ic;
  FFVar IC_T  = T  - T_ic;

  FFVar BC_L1 = U_VEL*c1 - D_ax*OpP( c1, z ) - U_VEL*c1feed;
  FFVar BC_L2 = U_VEL*c2 - D_ax*OpP( c2, z ) - U_VEL*c2feed;
  FFVar BC_LT = G_cv*T  - lam*OpP( T,  z ) - G_cv*T0_ref;
  FFVar BC_U1 = OpP( c1, z ), BC_U2 = OpP( c2, z ), BC_UT = OpP( T, z );

  FFVar Inv1 = OpI( c1 + F_ph*q1, z );    // fct 0 (terminal)
  FFVar Inv2 = OpI( c2 + F_ph*q2, z );    // fct 1 (terminal)
  FFVar Eff1 = OpI( U_VEL*c1, t );        // fct 2 (accumulated over the march)
  FFVar Eff2 = OpI( U_VEL*c2, t );        // fct 3 (accumulated over the march)

  OCFESLV oc( &DAG );
  oc.add_domain( t, FFDom( t_bnd, coltype, n_nd ) );        // explicit-boundary (non-uniform) ctor
  oc.add_domain( z, FFDom( 0., 1.0, n_el, coltype, n_nd ) );
  oc.add_state ( c1, {t,z} );  oc.add_state ( c2, {t,z} );
  oc.add_state ( q1, {t,z} );  oc.add_state ( q2, {t,z} );
  oc.add_state ( T,  {t,z} );
  oc.add_input ( b01, b01_nom, /*is_decision=*/true );

  // update_ref LIFETIME: these callables are RETAINED and invoked during init().  c1g/c2g
  // are captured BY VALUE into the q1/q2 lambdas -- capturing a local helper lambda by
  // reference dangles and crashes far from the cause.
  auto c1g = [&t,&z,fd]( OCFESLV::t_Coord const& cr ){ return fd.c1( cr.at(t) )*( 1.0 - 0.5*cr.at(z) ); };
  auto c2g = [&t,&z,fd]( OCFESLV::t_Coord const& cr ){ return fd.c2( cr.at(t) )*( 1.0 - 0.5*cr.at(z) ); };
  oc.update_ref( c1, c1g );
  oc.update_ref( c2, c2g );
  oc.update_ref( q1, [c1g,c2g]( OCFESLV::t_Coord const& cr ){ return q1star_g( c1g(cr), c2g(cr), b01_nom ); } );
  oc.update_ref( q2, [c1g,c2g]( OCFESLV::t_Coord const& cr ){ return q2star_g( c1g(cr), c2g(cr), b01_nom ); } );
  oc.update_ref( T,  [c1g,c2g]( OCFESLV::t_Coord const& cr ){
    return T0_ref + F_ph*( dH1*q1star_g(c1g(cr),c2g(cr),b01_nom)
                         + dH2*q2star_g(c1g(cr),c2g(cr),b01_nom) )/Cp_e; } );

  oc.add_input ( c1_ic, {z} );  oc.add_input ( c2_ic, {z} );
  oc.add_input ( q1_ic, {z} );  oc.add_input ( q2_ic, {z} );
  oc.add_input ( T_ic, {z}, []( OCFESLV::t_Coord const& ){ return T0_ref; }, /*is_decision=*/true );
  oc.update_ref( c1_ic, []( OCFESLV::t_Coord const& ){ return 0.0; } );
  oc.update_ref( c2_ic, []( OCFESLV::t_Coord const& ){ return 0.0; } );
  oc.update_ref( q1_ic, []( OCFESLV::t_Coord const& ){ return 0.0; } );
  oc.update_ref( q2_ic, []( OCFESLV::t_Coord const& ){ return 0.0; } );

  oc.set_evolution_domain( t );

  int const T_NO_LB = FFDom::ALL - FFDom::LB;
  int const Z_INT   = FFDom::ALL - FFDom::LB - FFDom::UB;
  oc.add_equation( CONT1, {t,z}, {T_NO_LB,   Z_INT},      OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( CONT2, {t,z}, {T_NO_LB,   Z_INT},      OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( LDF1,  {t,z}, {T_NO_LB,   FFDom::ALL}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( LDF2,  {t,z}, {T_NO_LB,   FFDom::ALL}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( ENE_T, {t,z}, {T_NO_LB,   Z_INT},      OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( IC_c1, {t,z}, {FFDom::LB, FFDom::ALL}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INITIAL,  0 ) );
  oc.add_equation( IC_c2, {t,z}, {FFDom::LB, FFDom::ALL}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INITIAL,  0 ) );
  oc.add_equation( IC_q1, {t,z}, {FFDom::LB, FFDom::ALL}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INITIAL,  0 ) );
  oc.add_equation( IC_q2, {t,z}, {FFDom::LB, FFDom::ALL}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INITIAL,  0 ) );
  oc.add_equation( IC_T,  {t,z}, {FFDom::LB, FFDom::ALL}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INITIAL,  0 ) );
  oc.add_equation( BC_L1, {t,z}, {T_NO_LB,   FFDom::LB},  OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );
  oc.add_equation( BC_L2, {t,z}, {T_NO_LB,   FFDom::LB},  OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );
  oc.add_equation( BC_LT, {t,z}, {T_NO_LB,   FFDom::LB},  OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );
  oc.add_equation( BC_U1, {t,z}, {T_NO_LB,   FFDom::UB},  OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );
  oc.add_equation( BC_U2, {t,z}, {T_NO_LB,   FFDom::UB},  OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );
  oc.add_equation( BC_UT, {t,z}, {T_NO_LB,   FFDom::UB},  OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );
  oc.add_output( Inv1, {t}, {T_end} );   // fct 0
  oc.add_output( Inv2, {t}, {T_end} );   // fct 1
  oc.add_output( Eff1, {z}, {1.0}   );   // fct 2
  oc.add_output( Eff2, {z}, {1.0}   );   // fct 3

  oc.options.REDUCE.ORDER     = OCFESLV::Options::RED_FULL;
  oc.options.CLASSIFY.MODE         = OCFESLV::Options::CLASS_AUTO;
  oc.options.INTERFACE.TYPE   = OCFESLV::Options::IC_AUTO;
  oc.options.INTERFACE.IMPOSITION  = OCFESLV::Options::IC_STRONG;
  oc.options.INTERFACE.SAT_SIGMA0       = 10.0;
  oc.options.DISPLAY_LEVEL    = 0;
#if defined(CRONOS__WITH_SPQR)
  oc.options.SOLVE.FACTORIZATION = OCFESLV::Options::SOLVE_SPQR;
#else
  oc.options.SOLVE.FACTORIZATION = OCFESLV::Options::SOLVE_SUPERLU;
#endif
  oc.options.SOLVE.WARMSTART     = OCFESLV::Options::BROADCAST_IC;   // marching, as PSA5
  oc.options.OUTPUT.MARCH_STORE = true;                        // needed by eval_solution()

  if( !oc.setup() ){ std::cerr << "  setup() FAILED\n"; return R; }
  std::vector<double> xv, inp;
  if( !oc.init( xv, inp, nullptr ) ){ std::cerr << "  init() FAILED\n"; return R; }

  R.nVar = oc.n_colloc_sta();
  R.nwin = oc.n_march_steps();
  size_t const ncf = oc.n_colloc_fct(), ncd = oc.n_control_dof();

  OCFESLV::SolveReport const rep = oc.solve( xv.data(), inp.data(), nullptr );
  R.converged = rep.converged; R.iters = rep.iterations;
  std::cout << "  nVar=" << R.nVar << " marching=" << (oc.is_marching()?"yes":"no")
            << " windows=" << R.nwin << " ncf=" << ncf << " ncd=" << ncd
            << " converged=" << (rep.converged?"yes":"no")
            << " iters=" << rep.iterations << " final|r|=" << std::scientific
            << std::setprecision(3) << rep.final_residual << "\n";
  if( !rep.converged ) return R;

  // ---- (i) the KPI balance -- what an optimizer sees ----
  std::vector<double> const Fref = oc.val_functions();
  double const Feed1T = U_VEL*FI(1,T_end), Feed2T = U_VEL*FI(2,T_end);
  if( Fref.size() >= 4 ){
    R.mb_kpi1 = std::fabs( Fref[0] + Fref[2] - Feed1T )/std::max(Feed1T,1e-30);
    R.mb_kpi2 = std::fabs( Fref[1] + Fref[3] - Feed2T )/std::max(Feed2T,1e-30);
  }

  // ---- field reads MUST happen here: eval_solution() serves the MOST RECENT solve,
  //      and the fsens / FD solves below would overwrite the stored trajectory ----
  auto S = [&]( FFVar const& v, double zz, double tt ){
    OCFESLV::t_Coord pt; pt[z]=zz; pt[t]=tt; return oc.eval_solution( v, pt ); };

  // outlet history (root-identity guard)
  R.c1out.resize( NTS );
  for( size_t kk=0; kk<NTS; ++kk ) R.c1out[kk] = S( c1, 1.0, tsamp(kk) );

  // ---- (ii) exact GL re-quadrature -> per-window conservation defect ----
  std::vector<double> zb( n_el+1 );                              // z-elements (uniform)
  for( size_t i=0; i<=n_el; ++i ) zb[i] = double(i)/double(n_el);
  auto inv_at = [&]( double tt, double& I1, double& I2 ){
    I1=0.; I2=0.;
    for( size_t ie=0; ie+1<zb.size(); ++ie ){
      double const a=zb[ie], b=zb[ie+1], hm=0.5*(b-a), mid=0.5*(a+b);
      for( int g=0; g<10; ++g ){
        double const zz=mid+hm*GLX[g], w=hm*GLW[g];
        I1 += w*( S(c1,zz,tt) + F_ph*S(q1,zz,tt) );
        I2 += w*( S(c2,zz,tt) + F_ph*S(q2,zz,tt) ); } } };

  R.wt = t_bnd; R.wd1.assign( t_bnd.size()-1, 0. ); R.wd2.assign( t_bnd.size()-1, 0. );
  double Iprev1, Iprev2; inv_at( t_bnd[0], Iprev1, Iprev2 );
  double eff1=0., eff2=0., cum1=0., cum2=0.;
  for( size_t k=0; k+1<t_bnd.size(); ++k ){
    double const a=t_bnd[k], b=t_bnd[k+1];
    // ONE GL-10 panel per window is EXACT, not approximate: within a window the marched
    // solution is a polynomial of degree <= n_nd-1 = 5 in t, and GL-10 integrates degree
    // 19.  Sub-panelling would only cost evaluations without adding a digit.
    double const hm=0.5*(b-a), mid=0.5*(a+b);
    double de1=0., de2=0.;
    for( int g=0; g<10; ++g ){
      double const tt=mid+hm*GLX[g], w=hm*GLW[g];
      de1 += w*U_VEL*S(c1,1.0,tt); de2 += w*U_VEL*S(c2,1.0,tt); }
    eff1 += de1; eff2 += de2;
    double I1,I2; inv_at( b, I1, I2 );
    double const dF1 = U_VEL*( FI(1,b) - FI(1,a) ), dF2 = U_VEL*( FI(2,b) - FI(2,a) );
    R.wd1[k] = ( (I1-Iprev1) - ( dF1 - de1 ) )/std::max(Feed1T,1e-30);
    R.wd2[k] = ( (I2-Iprev2) - ( dF2 - de2 ) )/std::max(Feed2T,1e-30);
    if( std::fabs(R.wd1[k]) > R.dwin_max1 ){ R.dwin_max1 = std::fabs(R.wd1[k]); R.kwin_max1 = k; }
    R.dwin_max2 = std::max( R.dwin_max2, std::fabs(R.wd2[k]) );
    cum1 += R.wd1[k]; cum2 += R.wd2[k];
    Iprev1=I1; Iprev2=I2;
  }
  R.mb_gl1 = std::fabs(cum1); R.mb_gl2 = std::fabs(cum2);

  // ---- over/undershoot: the feed brackets the solution, so any excursion outside
  //      [min(cA,cB), max(cA,cB)] is a numerical artefact of the front ----
  {
    double m1=1e30,M1=-1e30,m2=1e30,M2=-1e30;
    int const NZ=40, NT=400;
    auto scan = [&]( double t_lo, double t_hi ){
      for( int iz=0; iz<=NZ; ++iz ){ double const zz=double(iz)/double(NZ);
        for( int it=0; it<=NT; ++it ){
          double const tt = t_lo + (t_hi-t_lo)*double(it)/double(NT);
          if( tt < 0. || tt > T_end ) continue;
          double const v1=S(c1,zz,tt), v2=S(c2,zz,tt);
          m1=std::min(m1,v1); M1=std::max(M1,v1);
          m2=std::min(m2,v2); M2=std::max(M2,v2); } } };
    scan( 0., T_end );                                   // global pass
    scan( std::max(0.,t_step-8.0*W), std::min(T_end,t_step+8.0*W) );   // front pass
    scan( 0., std::min( T_end, 8.0*std::max(cfg.eps0,cfg.eps) ) );     // t=0 corner pass
    R.under1 = std::min( 0.0, m1 );                      // any excursion below 0
    R.under2 = std::min( 0.0, m2 );
    R.over1  = std::max( 0.0, M1 - std::max(cA1,cB1) );  // any excursion above the feed
    R.over2  = std::max( 0.0, M2 - std::max(cA2,cB2) );
  }

  std::cout << std::scientific << std::setprecision(3)
            << "  massbal(KPI): c1=" << R.mb_kpi1 << " c2=" << R.mb_kpi2
            << "  |  massbal(GL-exact): c1=" << R.mb_gl1 << " c2=" << R.mb_gl2 << "\n"
            << "  worst window: k=" << R.kwin_max1 << " on [" << std::fixed << std::setprecision(3)
            << t_bnd[R.kwin_max1] << "," << t_bnd[R.kwin_max1+1] << "] defect="
            << std::scientific << std::setprecision(3) << R.dwin_max1
            << "   overshoot c1=" << R.over1 << " undershoot c1=" << R.under1 << "\n";

  // ---- (iii) gradient contamination -- LAST, since it overwrites the trajectory ----
  if( cfg.want_sens && ncd > 0 && ncf >= 4 ){
    std::vector<double> xvf( xv ), inpf( inp );
    if( oc.solve_fsens( xvf.data(), inpf.data(), nullptr ) ){
      std::vector<double> const& J = oc.sens_jacobian();
      auto const& C = oc.controls();
      size_t const ib01 = C.at( b01 ).offset;
      double gi=0., fi=0.;
      for( size_t j=0; j<ncd; ++j ){
        double const g = J[0*ncd+j] + J[2*ncd+j];            // dInv1/dp + dEff1/dp
        gi = std::max( gi, std::fabs(g) );
        fi = std::max( fi, std::fabs( J[0*ncd+j] ) );
      }
      R.dMB1_b01 = ( J[0*ncd+ib01] + J[2*ncd+ib01] )/std::max(Feed1T,1e-30);
      R.dMB1_inf = gi/std::max(Feed1T,1e-30);
      R.dF1_b01  = J[0*ncd+ib01]/std::max(Feed1T,1e-30);
      R.dF1_inf  = fi/std::max(Feed1T,1e-30);
      std::cout << "  gradient contamination: |dMB1/dp|_inf=" << R.dMB1_inf
                << "  vs |dInv1/dp|_inf=" << R.dF1_inf
                << "   ratio=" << ( R.dF1_inf>0? R.dMB1_inf/R.dF1_inf : 0. ) << "\n";
    }
    else std::cerr << "  solve_fsens FAILED (gradient contamination not measured)\n";
  }

  return R;
}

// ---------------------------------------------------------------------------
int main()
{
  std::cout << "================================================================\n"
            << "  PSA6 -- flux-conservation STEP 1 on the PSA5 step-inlet oracle\n"
            << "  graded mesh at EQUAL COST vs a sharpening INTERIOR front\n"
            << "================================================================\n";

  std::vector<double> const epss = { 0.2, 0.1, 0.05, 0.025, 0.0125 };
  std::vector<size_t> const Gs   = { 0, 2, 4 };

  std::vector<SR> rows;
  for( size_t G : Gs )
    for( double eps : epss ){
      Cfg c; c.eps=eps; c.G=G; c.tag="sweep";
      c.want_sens = ( eps == epss.back() );        // gradient probe at the sharpest eps
      rows.push_back( run_case( c ) );
    }

  // startup-isolation probe: same two configs with the t=0 corner step smoothed away.
  // If the drift collapses here, the corner -- not the interior front -- dominates.
  for( size_t G : { size_t(0), size_t(4) } ){
    Cfg c; c.eps=epss.back(); c.G=G; c.smooth_start=true; c.eps0=epss.back();
    c.tag="probe:smooth-start";
    rows.push_back( run_case( c ) );
  }

  // root-identity guard: each run vs the SAME-eps, SAME-startup, G=0 baseline
  for( SR& R : rows ){
    if( R.cfg.G==0 ){ R.outdev = 0.; continue; }
    for( SR const& B : rows )
      if( B.cfg.G==0 && B.cfg.eps==R.cfg.eps && B.cfg.smooth_start==R.cfg.smooth_start
       && B.converged && R.converged && B.c1out.size()==R.c1out.size() ){
        double m=0.; for( size_t k=0;k<R.c1out.size();++k )
          m=std::max(m,std::fabs(R.c1out[k]-B.c1out[k]));
        R.outdev=m; break; }
  }

  // ---------------- the table ----------------
  std::cout << "\n=============== PSA6 STEP-1 TABLE (equal cost: " << NE_TOT
            << " t-elements throughout) ===============\n";
  std::cout << std::left
            << std::setw(7) << "start" << std::setw(4) << "G" << std::setw(9) << "eps"
            << std::right
            << std::setw(11) << "dt1" << std::setw(9) << "dt1/eps" << std::setw(5) << "it"
            << std::setw(11) << "mb_KPI" << std::setw(11) << "mb_GL"
            << std::setw(11) << "worstwin" << std::setw(5) << "k*"
            << std::setw(11) << "under1" << std::setw(11) << "outdev" << "\n";
  for( SR const& R : rows ){
    std::cout << std::left
              << std::setw(7) << (R.cfg.smooth_start? "smooth":"step")
              << std::setw(4) << R.cfg.G
              << std::fixed << std::setprecision(4) << std::setw(9) << R.cfg.eps
              << std::right << std::scientific << std::setprecision(2) << std::setw(11) << R.dt1
              << std::fixed << std::setprecision(3) << std::setw(9) << R.ratio
              << std::setw(5) << R.iters;
    if( !R.converged ){ std::cout << "   *** DIVERGED ***\n"; continue; }
    std::cout << std::scientific << std::setprecision(2)
              << std::setw(11) << R.mb_kpi1 << std::setw(11) << R.mb_gl1
              << std::setw(11) << R.dwin_max1 << std::setw(5) << R.kwin_max1
              << std::setw(11) << R.under1 << std::setw(11) << R.outdev << "\n";
  }

  std::cout << "\n--- gradient contamination at the sharpest eps ---\n"
            << std::left << std::setw(7) << "start" << std::setw(4) << "G"
            << std::right << std::setw(13) << "|dMB1/dp|inf" << std::setw(13) << "|dInv1/dp|inf"
            << std::setw(11) << "ratio" << "\n";
  for( SR const& R : rows ){
    if( !R.cfg.want_sens || !R.converged ) continue;
    std::cout << std::left << std::setw(7) << (R.cfg.smooth_start? "smooth":"step")
              << std::setw(4) << R.cfg.G << std::right << std::scientific << std::setprecision(3)
              << std::setw(13) << R.dMB1_inf << std::setw(13) << R.dF1_inf
              << std::setw(11) << ( R.dF1_inf>0? R.dMB1_inf/R.dF1_inf : 0. ) << "\n";
  }

  std::cout << "\n  READING THE TABLE\n"
               "  * dt1/eps is pinned along each G>0 column by construction, so the test is\n"
               "    FLATNESS IN eps, not a fitted collapse.\n"
               "  * mb_KPI vs mb_GL: if they differ, the defect is framework QUADRATURE, not\n"
               "    the scheme -- fix that before scoping any DG work.\n"
               "  * k* names the window that loses the mass.  k*=0 means the t=0 corner step\n"
               "    dominates (compare the smooth-start probe); k* at the zone means the\n"
               "    interior front does; k* scattered means background quadrature error.\n"
               "  * outdev is the root-identity guard: O(1e-3) refinement, O(1) branch change.\n"
               "  * ratio ~ 1 means the conservation defect moves with the controls and WILL\n"
               "    corrupt reduced-space gradients; ratio << 1 means it is a benign offset.\n";

  // ---------------- data files ----------------
  std::string const pre = std::string(PSA6_OUT_PREFIX);
  { std::ofstream f( pre+"_summary.out" );
    f << "# start G eps dt1 dt1_over_eps iters mb_KPI1 mb_GL1 mb_KPI2 mb_GL2 worstwin kstar under1 outdev conv\n";
    for( SR const& R : rows )
      f << (R.cfg.smooth_start?"smooth":"step") << " " << R.cfg.G << " "
        << std::setprecision(8) << R.cfg.eps << " " << R.dt1 << " " << R.ratio << " "
        << R.iters << " " << R.mb_kpi1 << " " << R.mb_gl1 << " "
        << R.mb_kpi2 << " " << R.mb_gl2 << " " << R.dwin_max1 << " " << R.kwin_max1 << " "
        << R.under1 << " " << R.outdev << " " << (R.converged?1:0) << "\n"; }

  { std::ofstream f( pre+"_window.out" );
    f << "# per-window conservation defect: one block per config\n";
    for( SR const& R : rows ){
      if( !R.converged ) continue;
      f << "\n\n# start=" << (R.cfg.smooth_start?"smooth":"step")
        << " G=" << R.cfg.G << " eps=" << R.cfg.eps << "\n";
      f << "# t_lo t_hi d1 d2\n";
      for( size_t k=0; k<R.wd1.size(); ++k )
        f << std::setprecision(8) << R.wt[k] << " " << R.wt[k+1] << " "
          << R.wd1[k] << " " << R.wd2[k] << "\n"; } }

  { std::ofstream f( pre+"_outlet.out" );
    f << "# t  c1(1,t) : sharpest eps, one column per G (step startup)\n";
    std::vector<SR const*> sel;
    for( SR const& R : rows )
      if( !R.cfg.smooth_start && R.cfg.eps==epss.back() && R.converged ) sel.push_back(&R);
    for( size_t kk=0; kk<NTS; ++kk ){
      f << std::setprecision(8) << tsamp(kk);
      for( SR const* R : sel ) f << " " << R->c1out[kk];
      f << "\n"; } }

  { std::ofstream gp( pre+".gp" );
    gp << "set datafile commentschars '#'\n"
          "set terminal pngcairo size 1000,700 enhanced font 'Helvetica,12'\n"
          "set grid\n\n"
          "set output '" << pre << "_flatness.png'\n"
          "set title 'conservation defect vs eps -- FLAT along a graded column = pure resolution'\n"
          "set xlabel 'eps'; set ylabel '|mass-balance defect|'; set logscale xy; set key top right\n"
          "plot '" << pre << "_summary.out' u ($2==0?$3:1/0):($8>1e-16?$8:1e-16) w lp lw 2 pt 7 t 'G=0 (uniform)', \\\n"
          "     '' u ($2==2?$3:1/0):($8>1e-16?$8:1e-16) w lp lw 2 pt 5 t 'G=2', \\\n"
          "     '' u ($2==4?$3:1/0):($8>1e-16?$8:1e-16) w lp lw 2 pt 9 t 'G=4'\n"
          "unset logscale\n\n"
          "set output '" << pre << "_window.png'\n"
          "set title 'per-window conservation defect: WHERE the mass is lost'\n"
          "set xlabel 't (window lower edge)'; set ylabel 'defect'; set key outside\n"
          "plot for [i=0:4] '" << pre << "_window.out' index i u 1:3 w steps lw 2 t 'block '.i\n\n"
          "set output '" << pre << "_outlet.png'\n"
          "set title 'outlet c_1(1,t): must be grading-INVARIANT (root-identity guard)'\n"
          "set xlabel 't'; set ylabel 'c_1(1,t)'; set key bottom right\n"
          "plot for [i=2:4] '" << pre << "_outlet.out' u 1:i w l lw 2 t 'col '.i\n"; }

  bool all_conv=true; for( SR const& R : rows ) all_conv &= R.converged;
  std::cout << "\n  wrote " << pre << "_{summary,window,outlet}.out  and  " << pre << ".gp\n"
            << "  Overall: " << ( all_conv ? "ALL CONVERGED" : "SOME DID NOT CONVERGE" ) << "\n";
  return all_conv ? 0 : 1;
}
