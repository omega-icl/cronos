// OCFE_PSA7.cpp  ---  IS THE MESH-RATIO FLOOR REMOVABLE BY MESH DESIGN ALONE?
// ===========================================================================
// Follow-on to OCFE_PSA6, same oracle (PSA5 step-inlet PSA, marching, controls b01+T_ic).
// PSA6 established four things, and this driver exists because of the fourth:
//
//   1. mb_KPI == mb_GL in all 17 rows -> the framework's OpI quadrature is EXACT; the
//      defect lives entirely in the discrete solution.  (Re-checked here as a regression.)
//   2. At G=0 the defect is 100% localised: eight windows at ~1e-6, one window straddling
//      the step at 1.2e-2.  Background strong-form conservation is ALREADY 1e-6.
//   3. Grading at equal cost kills the front defect: 1.24e-2 -> 8.5e-6 at G=2 (x1500).
//   4. *** MORE GRADING IS WORSE. ***  G=4 floors at ~1.1e-4 -- 14x worse than G=2 -- and
//      k* moves off the front onto window 6 = [2.05, 3.033], the FIRST COARSE WINDOW AFTER
//      THE ZONE.  That window alone is 1.10e-4 of the 1.18e-4 total.  Element-size jump out
//      of the zone: 39x at G=4 versus 15x at G=2, where the same window contributes 8.8e-6.
//
// THE QUESTION
// ------------
// Is that floor a MESH-RATIO artefact -- removable by not slamming a 0.025-wide element
// into a 0.98-wide one -- or is it an irreducible property of the strong-form operator?
// The answer decides how much of the flux-form DG case survives PSA6:
//   * removable by mesh design  -> the remaining DG justification is narrow and specific:
//     fronts whose LOCATION is not known a priori (a real PSA shock, a cyclic process where
//     t* moves), since grading presupposes knowing where to grade.  Different verification
//     ladder than the one scoped in 0r, and it should be re-scoped as such.
//   * NOT removable -> telescoping is exact regardless of element-size ratio, and option 2
//     is back on its original footing.
//
// THE MESH LAW (the only thing that changes from PSA6)
// ----------------------------------------------------
// ABRUPT  == PSA6 verbatim: G uniform elements on [t*-W, t*+W], W = 4*eps, uniform fill on
//            each side.  Reproduced EXACTLY (checked: G=4, eps=0.0125 gives
//            0 0.975 1.95 1.975 2 2.025 2.05 3.033 4.017 5, PSA6's mesh), so the abrupt rows
//            are a CROSS-DRIVER REGRESSION -- they must reproduce PSA6's mb_GL1 to the digit.
// GEOM    == same zone, but each side fills geometrically outward from the zone width h0.
//            Given a side length L and n elements, the growth ratio rho solves
//                 h0 * rho * (rho^n - 1)/(rho - 1) = L
//            exactly (bisection, monotone).  The left/right element split is then chosen to
//            MINIMISE the achieved ratio.  So rho is a MEASURED OUTCOME of the budget, not a
//            requested knob -- a requested-ratio law is not monotone in the achieved ratio
//            (a greedy fill leaves a runt element at the outer end and the ratio blows up
//            there instead), which would have made the experiment meaningless.
//
// TWO ARMS -- "removable at all" is a different question from "removable for free"
// ---------------------------------------------------------------------------------
//   ARM A  budget = 9 elements (PSA6's cost), ABRUPT vs GEOM, G in {2,3,4}, eps in
//          {0.05, 0.0125}.  Geometric fill cuts rho at equal cost:
//              G=4, eps=0.0125:  rho 39.33 -> 8.35        G=2:  14.75 -> 3.00
//          If d_post falls with it, the floor is a ratio artefact and partly free.
//   ARM B  budget in {12,16,20}, GEOM, G in {2,4}, eps = 0.0125.  rho falls monotonically
//          (G=4: 8.35 -> 2.99 -> 1.97 -> 1.61).  This is the "removable at all" arm, and
//          the element count is the honest price to weigh against building DG.
//
// WHAT IS MEASURED
// ----------------
// The per-window defect is DECOMPOSED, because the whole point is that the mass moves from
// one place to another:
//     d_zone  sum |d_k| over the refinement-zone windows      (the front itself)
//     d_post  |d_k| on the first window AFTER the zone         (PSA6's culprit)
//     d_pre   |d_k| on the last window BEFORE the zone
//     d_bg    sum |d_k| over everything else                   (the ~1e-6 background)
// Headline: d_post versus rho.  If it collapses on rho, the floor is mesh design.
// Retained regressions: mb_KPI == mb_GL (quadrature layer), and under1 pinned at
// -1.6e-3..-1.9e-3 -- PSA6 showed it is invariant to every t-direction knob, hence a
// z-direction artefact; the z-mesh is untouched here, so it MUST NOT move.  If it does,
// something else changed and the comparison is void.
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

#ifndef PSA7_OUT_PREFIX
#define PSA7_OUT_PREFIX "OCFE_PSA7"
#endif

// ---------- model: IDENTICAL to PSA5/PSA6 (rung-9 / PDE36c).  Do not retune. ----------
static double const U_VEL = 1.0, D_ax = 0.1, F_ph = 0.5, k_ldf = 2.0;
static double const qs1 = 1.0, qs2 = 1.0, b02_L = 1.0;
static double const beta1 = 2.0, beta2 = 1.0, T0_ref = 1.0;
static double const Cp_e = 1.0, G_cv = 1.0, lam = 0.1, dH1 = 1.0, dH2 = 0.5, hw = 0.5, Tw = 1.0;
static double const T_end = 5.0, t_step = 2.0;
static double const cA1 = 0.2, cB1 = 0.5, cA2 = 0.5, cB2 = 0.2;
static double const b01_nom = 3.0;
static double const ZONE_KW = 4.0;               // zone half-width W = ZONE_KW*eps
static size_t const NTS = 501;

static inline double tsamp( size_t kk ){ return double(kk)/double(NTS-1)*T_end; }
static inline double q1star_g( double c1g, double c2g, double b01 )
{ return qs1*b01*c1g/( 1.0 + b01*c1g + b02_L*c2g ); }
static inline double q2star_g( double c1g, double c2g, double b01 )
{ return qs2*b02_L*c2g/( 1.0 + b01*c1g + b02_L*c2g ); }

// ---------- feed: plain step, no startup smoothing.  PSA6 settled that PSA5's unsmoothed
// startup is NOT a numerical problem: the discontinuity sits exactly on the t=0 corner
// node, which T_NO_LB excludes from the BC, so the scheme never samples it.  (PSA6's
// smooth-start probe merely manufactured a NEW unresolved front in window 0, d1 going
// -1.2e-5 -> +5.3e-3, leaving every other window bit-identical.)  Dropped here. ----------
static inline double lncosh( double x )
{ x = std::fabs(x); return x + std::log1p( std::exp(-2.0*x) ) - std::log(2.0); }

struct Feed
{
  double eps = 0.05;
  double stepf( double t ) const { return 0.5*( 1.0 + std::tanh( (t-t_step)/eps ) ); }
  double c1( double t ) const { return cA1 + (cB1-cA1)*stepf(t); }
  double c2( double t ) const { return cA2 + (cB2-cA2)*stepf(t); }
  //! exact \int_0^t c_i,feed
  double integ( int i, double t ) const
  {
    double const A = ( i==1? cA1 : cA2 ), B = ( i==1? cB1-cA1 : cB2-cA2 );
    return A*t + B*0.5*( t + eps*( lncosh( (t-t_step)/eps ) - lncosh( t_step/eps ) ) );
  }
};

// ---------- Gauss-Legendre 10 (exact to degree 19; the solution is degree <= n_nd-1 = 5
// per element, so ONE panel per element is exact, not approximate) ----------
static double const GLX[10] = { -0.9739065285171717, -0.8650633666889845, -0.6794095682990244,
                                -0.4333953941292472, -0.1488743389816312,  0.1488743389816312,
                                 0.4333953941292472,  0.6794095682990244,  0.8650633666889845,
                                 0.9739065285171717 };
static double const GLW[10] = {  0.0666713443086881,  0.1494513491505806,  0.2190863625159820,
                                 0.2692667193099963,  0.2955242247147529,  0.2955242247147529,
                                 0.2692667193099963,  0.2190863625159820,  0.1494513491505806,
                                 0.0666713443086881 };

// ---------------------------------------------------------------------------
// mesh law
// ---------------------------------------------------------------------------
static double series( double h0, double rho, size_t n )
{
  if( n == 0 ) return 0.;
  if( std::fabs( rho - 1.0 ) < 1e-12 ) return h0*double(n);
  return h0*rho*( std::pow( rho, double(n) ) - 1.0 )/( rho - 1.0 );
}
//! growth ratio rho >= 1 such that n geometric elements starting at h0*rho exactly span L
static double solve_ratio( double h0, double L, size_t n )
{
  if( n == 0 || !(L > 0.) ) return 1.0;
  if( series( h0, 1.0, n ) >= L ) return 1.0;             // budget ample: uniform is finer
  double a = 1.0, b = 2.0;
  while( series( h0, b, n ) < L && b < 1e8 ) b *= 2.0;
  for( int i=0; i<200; ++i ){ double const m = 0.5*(a+b); if( series(h0,m,n) < L ) a=m; else b=m; }
  return 0.5*(a+b);
}
static std::vector<double> geom_widths( double h0, double L, size_t n, double rho )
{
  std::vector<double> w;
  if( n == 0 ) return w;
  if( rho <= 1.0 + 1e-12 ){ w.assign( n, L/double(n) ); return w; }
  for( size_t k=0; k<n; ++k ) w.push_back( h0*std::pow( rho, double(k+1) ) );
  double s=0.; for( double x : w ) s += x;
  for( double& x : w ) x *= L/s;
  return w;
}

struct MeshOut { std::vector<double> bnd; size_t kz_lo=0, kz_hi=0, nL=0, nR=0; double rho=1., h0=0.; };

static MeshOut build_mesh( double L, double ts, double W, size_t G, size_t budget, bool geom )
{
  MeshOut M;
  double const lo=ts-W, hi=ts+W, La=lo, Lc=L-hi;
  M.h0 = 2.0*W/double(G);
  size_t const avail = ( budget > G+1 ? budget-G : 2 );
  size_t bL=1, bR=avail-1; double best=1e300, rL=1., rR=1.;
  for( size_t nl=1; nl<avail; ++nl ){                     // split chosen to MINIMISE the ratio
    size_t const nr = avail-nl;
    double const a = geom? solve_ratio(M.h0,La,nl) : (La/double(nl))/M.h0;
    double const b = geom? solve_ratio(M.h0,Lc,nr) : (Lc/double(nr))/M.h0;
    double const m = std::max( std::max(a,1.0/std::max(a,1e-30)),
                               std::max(b,1.0/std::max(b,1e-30)) );
    if( m < best ){ best=m; bL=nl; bR=nr; rL=a; rR=b; }
  }
  std::vector<double> const wl = geom? geom_widths(M.h0,La,bL,rL) : std::vector<double>(bL,La/double(bL));
  std::vector<double> const wr = geom? geom_widths(M.h0,Lc,bR,rR) : std::vector<double>(bR,Lc/double(bR));
  M.nL=bL; M.nR=bR;
  M.bnd.push_back( 0. );
  for( size_t i=wl.size(); i-->0; ) M.bnd.push_back( M.bnd.back()+wl[i] );
  M.kz_lo = M.bnd.size()-1;
  for( size_t i=0; i<G; ++i ) M.bnd.push_back( M.bnd.back()+M.h0 );
  M.kz_hi = M.bnd.size()-1;
  for( size_t i=0; i<wr.size(); ++i ) M.bnd.push_back( M.bnd.back()+wr[i] );
  M.bnd.back() = L;
  M.rho = 1.;
  for( size_t k=0; k+2<M.bnd.size(); ++k ){
    double const a=M.bnd[k+1]-M.bnd[k], b=M.bnd[k+2]-M.bnd[k+1];
    M.rho = std::max( M.rho, std::max(a/b,b/a) ); }
  return M;
}

// ---------------------------------------------------------------------------
struct Cfg { double eps=0.05; size_t G=4, budget=9; bool geom=false; bool want_sens=false; };

struct SR {
  Cfg    cfg;
  MeshOut M;
  bool   converged=false;
  int    iters=0;
  size_t nVar=0, nwin=0;
  double mb_kpi1=1., mb_gl1=1., mb_kpi2=1., mb_gl2=1.;
  double d_zone=0., d_post=0., d_pre=0., d_bg=0.;
  double under1=0.;
  double dMB1_inf=0., dF1_inf=0.;
  double outdev=-1.;
  std::vector<double> c1out, wd1, wd2;
};

// ---------------------------------------------------------------------------
static SR run_case( Cfg const& cfg )
{
  SR R; R.cfg = cfg;
  size_t const n_el = 5, n_nd = 6;
  FFDom::TYPE const coltype = FFDom::LGL;

  Feed fd; fd.eps = cfg.eps;
  double const W = ZONE_KW*cfg.eps;
  R.M = build_mesh( T_end, t_step, W, cfg.G, cfg.budget, cfg.geom );
  std::vector<double> const& t_bnd = R.M.bnd;

  std::cout << "\n---- eps=" << std::fixed << std::setprecision(5) << cfg.eps
            << "  G=" << cfg.G << "  " << (cfg.geom? "GEOM  ":"ABRUPT")
            << "  budget=" << cfg.budget << "  ne=" << t_bnd.size()-1
            << "  rho=" << std::setprecision(3) << R.M.rho
            << "  nL/nR=" << R.M.nL << "/" << R.M.nR
            << "  zone=[" << R.M.kz_lo << "," << R.M.kz_hi << ") ----\n     t-mesh:";
  for( double b : t_bnd ) std::cout << " " << std::setprecision(4) << b;
  std::cout << "\n";

  // referee self-test: closed-form feed integral vs GL re-quadrature of the same integrand
  { double s1=0.; size_t const np=4000; double const hh=T_end/double(np);
    for( size_t k=0;k<np;++k ){ double const a=double(k)*hh, hm=0.5*hh, mid=a+hm;
      for( int g=0; g<10; ++g ) s1 += hm*GLW[g]*fd.c1( mid+hm*GLX[g] ); }
    double const e = std::fabs( s1 - fd.integ(1,T_end) );
    std::cout << "     feed-integral self-test: " << std::scientific << std::setprecision(2) << e
              << ( e < 1e-11 ? "  OK\n" : "  ** SUSPECT **\n" ); }

  // ---------------- model (== PSA5/PSA6) ----------------
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

  FFVar sstep  = 0.5*( 1.0 + tanh( ( t - t_step )/cfg.eps ) );
  FFVar c1feed = cA1 + ( cB1 - cA1 )*sstep;
  FFVar c2feed = cA2 + ( cB2 - cA2 )*sstep;

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

  FFVar Inv1 = OpI( c1 + F_ph*q1, z );    // fct 0
  FFVar Inv2 = OpI( c2 + F_ph*q2, z );    // fct 1
  FFVar Eff1 = OpI( U_VEL*c1, t );        // fct 2
  FFVar Eff2 = OpI( U_VEL*c2, t );        // fct 3

  OCFESLV oc( &DAG );
  oc.add_domain( t, FFDom( t_bnd, coltype, n_nd ) );
  oc.add_domain( z, FFDom( 0., 1.0, n_el, coltype, n_nd ) );   // z-mesh UNTOUCHED (see under1)
  oc.add_state ( c1, {t,z} );  oc.add_state ( c2, {t,z} );
  oc.add_state ( q1, {t,z} );  oc.add_state ( q2, {t,z} );
  oc.add_state ( T,  {t,z} );
  oc.add_input ( b01, b01_nom, /*is_decision=*/true );

  // update_ref LIFETIME: retained and invoked during init().  FFVars by reference (they
  // outlive init() here), everything else BY VALUE -- a local helper lambda captured by
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
  oc.add_output( Inv1, {t}, {T_end} );
  oc.add_output( Inv2, {t}, {T_end} );
  oc.add_output( Eff1, {z}, {1.0}   );
  oc.add_output( Eff2, {z}, {1.0}   );

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
  oc.options.SOLVE.WARMSTART        = OCFESLV::Options::BROADCAST_IC;
  oc.options.OUTPUT.MARCH_STORE = true;

  if( !oc.setup() ){ std::cerr << "  setup() FAILED\n"; return R; }
  std::vector<double> xv, inp;
  if( !oc.init( xv, inp, nullptr ) ){ std::cerr << "  init() FAILED\n"; return R; }

  R.nVar = oc.n_colloc_sta(); R.nwin = oc.n_march_steps();
  size_t const ncf = oc.n_colloc_fct(), ncd = oc.n_control_dof();

  OCFESLV::SolveReport const rep = oc.solve( xv.data(), inp.data(), nullptr );
  R.converged = rep.converged; R.iters = rep.iterations;
  std::cout << "  nVar=" << R.nVar << " windows=" << R.nwin << " ncf=" << ncf << " ncd=" << ncd
            << " converged=" << (rep.converged?"yes":"no") << " iters=" << rep.iterations
            << " final|r|=" << std::scientific << std::setprecision(3) << rep.final_residual << "\n";
  if( !rep.converged ) return R;

  std::vector<double> const Fref = oc.val_functions();
  double const Feed1T = U_VEL*fd.integ(1,T_end), Feed2T = U_VEL*fd.integ(2,T_end);
  if( Fref.size() >= 4 ){
    R.mb_kpi1 = std::fabs( Fref[0] + Fref[2] - Feed1T )/std::max(Feed1T,1e-30);
    R.mb_kpi2 = std::fabs( Fref[1] + Fref[3] - Feed2T )/std::max(Feed2T,1e-30);
  }

  // Field reads MUST precede the fsens solve: eval_solution() serves the MOST RECENT solve.
  auto S = [&]( FFVar const& v, double zz, double tt ){
    OCFESLV::t_Coord pt; pt[z]=zz; pt[t]=tt; return oc.eval_solution( v, pt ); };

  R.c1out.resize( NTS );
  for( size_t kk=0; kk<NTS; ++kk ) R.c1out[kk] = S( c1, 1.0, tsamp(kk) );

  std::vector<double> zb( n_el+1 );
  for( size_t i=0; i<=n_el; ++i ) zb[i] = double(i)/double(n_el);
  auto inv_at = [&]( double tt, double& I1, double& I2 ){
    I1=0.; I2=0.;
    for( size_t ie=0; ie+1<zb.size(); ++ie ){
      double const a=zb[ie], b=zb[ie+1], hm=0.5*(b-a), mid=0.5*(a+b);
      for( int g=0; g<10; ++g ){
        double const zz=mid+hm*GLX[g], w=hm*GLW[g];
        I1 += w*( S(c1,zz,tt) + F_ph*S(q1,zz,tt) );
        I2 += w*( S(c2,zz,tt) + F_ph*S(q2,zz,tt) ); } } };

  size_t const nw = t_bnd.size()-1;
  R.wd1.assign( nw, 0. ); R.wd2.assign( nw, 0. );
  double Ip1, Ip2; inv_at( t_bnd[0], Ip1, Ip2 );
  double cum1=0., cum2=0.;
  for( size_t k=0; k<nw; ++k ){
    double const a=t_bnd[k], b=t_bnd[k+1], hm=0.5*(b-a), mid=0.5*(a+b);
    double de1=0., de2=0.;
    for( int g=0; g<10; ++g ){                        // ONE panel per window: EXACT
      double const tt=mid+hm*GLX[g], w=hm*GLW[g];
      de1 += w*U_VEL*S(c1,1.0,tt); de2 += w*U_VEL*S(c2,1.0,tt); }
    double I1,I2; inv_at( b, I1, I2 );
    double const dF1 = U_VEL*( fd.integ(1,b) - fd.integ(1,a) );
    double const dF2 = U_VEL*( fd.integ(2,b) - fd.integ(2,a) );
    R.wd1[k] = ( (I1-Ip1) - ( dF1 - de1 ) )/std::max(Feed1T,1e-30);
    R.wd2[k] = ( (I2-Ip2) - ( dF2 - de2 ) )/std::max(Feed2T,1e-30);
    cum1 += R.wd1[k]; cum2 += R.wd2[k];
    Ip1=I1; Ip2=I2;
  }
  R.mb_gl1 = std::fabs(cum1); R.mb_gl2 = std::fabs(cum2);

  // decomposition: zone / immediately-after / immediately-before / background
  for( size_t k=0; k<nw; ++k ){
    double const v = std::fabs( R.wd1[k] );
    if( k >= R.M.kz_lo && k < R.M.kz_hi )        R.d_zone += v;
    else if( k == R.M.kz_hi )                    R.d_post  = v;
    else if( R.M.kz_lo>0 && k == R.M.kz_lo-1 )   R.d_pre   = v;
    else                                         R.d_bg   += v;
  }

  // under1 regression: PSA6 showed this is a z-artefact, invariant to every t-knob
  { double m1=1e30; int const NZ=30, NT=250;
    auto scan = [&]( double lo, double hi ){
      for( int iz=0; iz<=NZ; ++iz ){ double const zz=double(iz)/double(NZ);
        for( int it=0; it<=NT; ++it ){
          double const tt = lo + (hi-lo)*double(it)/double(NT);
          if( tt<0. || tt>T_end ) continue;
          m1 = std::min( m1, S(c1,zz,tt) ); } } };
    scan( 0., T_end );
    scan( std::max(0.,t_step-8.0*W), std::min(T_end,t_step+8.0*W) );
    R.under1 = std::min( 0.0, m1 ); }

  std::cout << std::scientific << std::setprecision(3)
            << "  mb_KPI1=" << R.mb_kpi1 << "  mb_GL1=" << R.mb_gl1
            << "   d_zone=" << R.d_zone << " d_post=" << R.d_post
            << " d_pre=" << R.d_pre << " d_bg=" << R.d_bg
            << "   under1=" << R.under1 << "\n";

  if( cfg.want_sens && ncd > 0 && ncf >= 4 ){
    std::vector<double> xvf( xv ), inpf( inp );
    if( oc.solve_fsens( xvf.data(), inpf.data(), nullptr ) ){
      std::vector<double> const& J = oc.sens_jacobian();
      double gi=0., fi=0.;
      for( size_t j=0; j<ncd; ++j ){
        gi = std::max( gi, std::fabs( J[0*ncd+j] + J[2*ncd+j] ) );
        fi = std::max( fi, std::fabs( J[0*ncd+j] ) ); }
      R.dMB1_inf = gi/std::max(Feed1T,1e-30);
      R.dF1_inf  = fi/std::max(Feed1T,1e-30);
      std::cout << "  gradient contamination: |dMB1/dp|inf=" << R.dMB1_inf
                << " vs |dInv1/dp|inf=" << R.dF1_inf
                << "  ratio=" << ( R.dF1_inf>0? R.dMB1_inf/R.dF1_inf : 0. ) << "\n";
    }
    else std::cerr << "  solve_fsens FAILED\n";
  }
  return R;
}

// ---------------------------------------------------------------------------
int main()
{
  std::cout << "================================================================\n"
            << "  PSA7 -- is the PSA6 mesh-ratio floor removable by mesh design?\n"
            << "  ARM A: equal cost, abrupt vs geometric.  ARM B: cost sweep.\n"
            << "================================================================\n";

  std::vector<Cfg> cfgs;
  // ARM A -- equal cost (9 elements, PSA6's budget); the ABRUPT rows are the PSA6 regression
  for( double eps : { 0.05, 0.0125 } )
    for( size_t G : { size_t(2), size_t(3), size_t(4) } )
      for( int gm=0; gm<2; ++gm ){
        Cfg c; c.eps=eps; c.G=G; c.budget=9; c.geom=(gm==1);
        c.want_sens = ( G==4 && eps==0.0125 );
        cfgs.push_back(c); }
  // ARM B -- cost sweep at the sharpest eps
  for( size_t G : { size_t(2), size_t(4) } )
    for( size_t bud : { size_t(12), size_t(16), size_t(20) } ){
      Cfg c; c.eps=0.0125; c.G=G; c.budget=bud; c.geom=true;
      c.want_sens = ( G==4 );
      cfgs.push_back(c); }

  std::vector<SR> rows;
  for( Cfg const& c : cfgs ) rows.push_back( run_case(c) );

  // root-identity guard: geometric vs the ABRUPT run at the same (eps, G) and budget 9
  for( SR& R : rows ){
    if( !R.cfg.geom ){ R.outdev = 0.; continue; }
    for( SR const& B : rows )
      if( !B.cfg.geom && B.cfg.G==R.cfg.G && B.cfg.eps==R.cfg.eps
       && B.converged && R.converged && B.c1out.size()==R.c1out.size() ){
        double m=0.; for( size_t k=0;k<R.c1out.size();++k )
          m=std::max(m,std::fabs(R.c1out[k]-B.c1out[k]));
        R.outdev=m; break; }
  }

  std::cout << "\n======================== PSA7 TABLE ========================\n"
            << std::left
            << std::setw(9) << "eps" << std::setw(4) << "G" << std::setw(8) << "mesh"
            << std::setw(5) << "ne"
            << std::right
            << std::setw(9) << "rho" << std::setw(5) << "it"
            << std::setw(11) << "mb_KPI1" << std::setw(11) << "mb_GL1"
            << std::setw(11) << "d_zone" << std::setw(11) << "d_post"
            << std::setw(11) << "d_pre" << std::setw(11) << "d_bg"
            << std::setw(11) << "under1" << std::setw(11) << "outdev" << "\n";
  for( SR const& R : rows ){
    std::cout << std::left << std::fixed << std::setprecision(4) << std::setw(9) << R.cfg.eps
              << std::setw(4) << R.cfg.G << std::setw(8) << (R.cfg.geom?"geom":"abrupt")
              << std::setw(5) << (R.M.bnd.size()-1)
              << std::right << std::setprecision(2) << std::setw(9) << R.M.rho
              << std::setw(5) << R.iters;
    if( !R.converged ){ std::cout << "   *** DIVERGED ***\n"; continue; }
    std::cout << std::scientific << std::setprecision(2)
              << std::setw(11) << R.mb_kpi1 << std::setw(11) << R.mb_gl1
              << std::setw(11) << R.d_zone  << std::setw(11) << R.d_post
              << std::setw(11) << R.d_pre   << std::setw(11) << R.d_bg
              << std::setw(11) << R.under1  << std::setw(11) << R.outdev << "\n";
  }

  std::cout << "\n--- gradient contamination (G=4) ---\n"
            << std::left << std::setw(9) << "eps" << std::setw(8) << "mesh" << std::setw(5) << "ne"
            << std::right << std::setw(9) << "rho" << std::setw(13) << "|dMB1/dp|"
            << std::setw(13) << "|dInv1/dp|" << std::setw(11) << "ratio" << "\n";
  for( SR const& R : rows ){
    if( !R.cfg.want_sens || !R.converged ) continue;
    std::cout << std::left << std::fixed << std::setprecision(4) << std::setw(9) << R.cfg.eps
              << std::setw(8) << (R.cfg.geom?"geom":"abrupt") << std::setw(5) << (R.M.bnd.size()-1)
              << std::right << std::setprecision(2) << std::setw(9) << R.M.rho
              << std::scientific << std::setprecision(3)
              << std::setw(13) << R.dMB1_inf << std::setw(13) << R.dF1_inf
              << std::setw(11) << ( R.dF1_inf>0? R.dMB1_inf/R.dF1_inf : 0. ) << "\n"; }

  std::cout << "\n  READING THE TABLE\n"
               "  * REGRESSION FIRST: the abrupt rows at (G=2,4) x (eps=0.05,0.0125) must\n"
               "    reproduce PSA6's mb_GL1 -- 6.9e-5 / 8.5e-6 / 7.2e-5 / 1.18e-4 -- to the\n"
               "    digit, and mb_KPI1 must equal mb_GL1 everywhere.  If either fails, the\n"
               "    model drifted and nothing below means anything.\n"
               "  * under1 must stay pinned at -1.6e-3..-1.9e-3 (z-artefact, z-mesh untouched).\n"
               "  * HEADLINE: d_post vs rho.  If d_post collapses on rho while d_zone and d_bg\n"
               "    stay put, the PSA6 floor is a MESH-RATIO artefact and mesh design fixes it.\n"
               "  * If d_post refuses to fall below ~1e-5 even at rho -> 1.6 (ARM B), the floor\n"
               "    is intrinsic to the strong form and flux-form telescoping is the remedy.\n"
               "  * ARM B's ne column is the PRICE of the fix -- weigh it against building DG.\n";

  std::string const pre = std::string(PSA7_OUT_PREFIX);
  { std::ofstream f( pre+"_summary.out" );
    f << "# eps G geom ne rho iters mb_KPI1 mb_GL1 mb_KPI2 mb_GL2 d_zone d_post d_pre d_bg under1 outdev nVar conv\n";
    for( SR const& R : rows )
      f << std::setprecision(8) << R.cfg.eps << " " << R.cfg.G << " " << (R.cfg.geom?1:0) << " "
        << R.M.bnd.size()-1 << " " << R.M.rho << " " << R.iters << " "
        << R.mb_kpi1 << " " << R.mb_gl1 << " " << R.mb_kpi2 << " " << R.mb_gl2 << " "
        << R.d_zone << " " << R.d_post << " " << R.d_pre << " " << R.d_bg << " "
        << R.under1 << " " << R.outdev << " " << R.nVar << " " << (R.converged?1:0) << "\n"; }

  { std::ofstream f( pre+"_window.out" );
    f << "# per-window defect; one block per config\n";
    for( SR const& R : rows ){
      if( !R.converged ) continue;
      f << "\n\n# eps=" << R.cfg.eps << " G=" << R.cfg.G << " " << (R.cfg.geom?"geom":"abrupt")
        << " ne=" << R.M.bnd.size()-1 << " rho=" << R.M.rho
        << " zone=[" << R.M.kz_lo << "," << R.M.kz_hi << ")\n# t_lo t_hi d1 d2\n";
      for( size_t k=0; k<R.wd1.size(); ++k )
        f << std::setprecision(8) << R.M.bnd[k] << " " << R.M.bnd[k+1] << " "
          << R.wd1[k] << " " << R.wd2[k] << "\n"; } }

  { std::ofstream gp( pre+".gp" );
    gp << "set datafile commentschars '#'\n"
          "set terminal pngcairo size 1000,700 enhanced font 'Helvetica,12'\n"
          "set grid\n\n"
          "set output '" << pre << "_dpost.png'\n"
          "set title 'post-zone defect vs achieved mesh ratio -- collapse => mesh artefact'\n"
          "set xlabel 'rho (max neighbouring element-size ratio)'; set ylabel 'd_{post}'\n"
          "set logscale xy; set key top left\n"
          "plot '" << pre << "_summary.out' u ($3==0?$5:1/0):($12>1e-16?$12:1e-16) "
          "w p pt 7 ps 1.4 t 'abrupt', \\\n"
          "     '' u ($3==1?$5:1/0):($12>1e-16?$12:1e-16) w p pt 5 ps 1.4 t 'geometric'\n"
          "unset logscale\n\n"
          "set output '" << pre << "_cost.png'\n"
          "set title 'ARM B: total defect vs element budget (the price of the fix)'\n"
          "set xlabel 'ne'; set ylabel 'mb_{GL1}'; set logscale y; set key top right\n"
          "plot '" << pre << "_summary.out' u ($2==2&&$3==1?$4:1/0):($8>1e-16?$8:1e-16) "
          "w lp lw 2 pt 7 t 'G=2', \\\n"
          "     '' u ($2==4&&$3==1?$4:1/0):($8>1e-16?$8:1e-16) w lp lw 2 pt 5 t 'G=4'\n"
          "unset logscale\n"; }

  bool all=true; for( SR const& R : rows ) all &= R.converged;
  std::cout << "\n  wrote " << pre << "_{summary,window}.out and " << pre << ".gp\n"
            << "  Overall: " << ( all ? "ALL CONVERGED" : "SOME DID NOT CONVERGE" ) << "\n";
  return all ? 0 : 1;
}
