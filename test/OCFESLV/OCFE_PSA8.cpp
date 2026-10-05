// OCFE_PSA8.cpp  ---  SETTLING THE DG QUESTION: the z-direction half
// ===========================================================================
// PSA6 and PSA7 dissolved the CONSERVATION motivation for flux-form DG, but they only
// ever varied the t-mesh; the z-mesh was held at uniform 5x6 as a control.  Since 0r
// proposes DG **in z**, two of its three motivations were never tested at all.  This
// driver tests both, on the same oracle (PSA5 step-inlet PSA, marching, b01+T_ic).
//
// The t-mesh is FROZEN at PSA7's best configuration -- G=2, geometric fill, 20 elements,
// which put the global balance at 1.6e-6, the background floor.  So whatever this driver
// measures is a z-direction property, not t-contamination.
//
// ---------------------------------------------------------------------------
// QUESTION A -- is the undershoot a RESOLUTION artefact or an OPERATOR artefact?
// ---------------------------------------------------------------------------
// under1 sat inert at -1.6e-3..-1.9e-3 through every t-knob in PSA6 and PSA7 (and moved
// only when PSA7's tight-budget geometric fill left a 2.32-wide TERMINAL element, i.e. by
// far-field t-coarsening, not by anything near the front).  0r's Gibbs diagnosis makes a
// SHARP, FALSIFIABLE prediction about which kind of refinement should fix it:
//
//     p-refinement (raise n_nd) only NARROWS the wiggle at fixed amplitude;
//     h-refinement (raise n_el) is what CONVERGES.
//
// So h and p are swept SEPARATELY here, at comparable cost, and the two hypotheses give
// opposite signatures:
//   * under1 -> 0 under h-refinement, flat under p-refinement  => Gibbs/resolution.  Mesh
//     fixes it; upwind interface dissipation buys nothing; the oscillation motivation for
//     DG is dead and the whole 0r case is dead with it.
//   * under1 plateaus under BOTH                                => it is the strong-form
//     advection operator itself.  Upwinding is the textbook remedy and DG has a live case
//     that has nothing to do with conservation.
// The LOCATION (z*,t*) of the minimum is reported -- it should have been all along.  If it
// sits at the z=0 Danckwerts layer (thickness D/U = 0.1, unresolved by a 0.2-wide first
// element) that is a boundary-layer resolution story and the graded-z row settles it; if it
// sits at the front or downstream of it, it is an advection story.
//
// ---------------------------------------------------------------------------
// QUESTION B -- would telescoping actually buy anything HERE?
// ---------------------------------------------------------------------------
// DG's guarantee is that element K sees +F_hat and K+1 sees -F_hat at a shared face, so
// summing local balances telescopes exactly.  But under RED_FULL + IC_STRONG this framework
// already carries BOTH c and its order-reduction auxiliary Dz_c as states with continuity
// claimed at interfaces -- so the total axial flux
//         N = U*c - D*dc/dz
// may ALREADY be single-valued at every z-face, in which case the telescoping DG advertises
// is structurally present and the refactor adds nothing at the face level.  That has never
// been checked.  It is checked here directly, and it is the single most decisive number in
// the driver:
//
//   FACE JUMPS.  At each interior z-face, the one-sided limits of c, Dz_c and N are formed
//   INDEPENDENTLY from the left and right element polynomials by 4-point Richardson
//   extrapolation to the face (samples at delta, 2delta, 3delta, 4delta inside each element;
//   v(0) = 4v0 - 6v1 + 4v2 - v3, error O(delta^4)).  Then:
//     jump -> 0 as delta -> 0  => flux already single-valued; DG's telescoping is ALREADY
//                                 satisfied and option 2 adds nothing but exactness of the
//                                 PER-ELEMENT balance;
//     jump -> nonzero limit    => the flux is genuinely double-valued, the sum does NOT
//                                 telescope, and option 2 has its original footing back.
//
//   PER-ELEMENT BALANCE.  For each z-element over the whole horizon (zero IC):
//     defect_j = \int_{zj}^{zj+1}(c+F q)(z,T) dz - [ \int_0^T N(zj,t) dt - \int_0^T N(zj+1,t) dt ]
//   Reported as max_j|defect_j| and sum_j|defect_j|.  Note the SIGNED sum telescopes by
//   construction here (each face flux is computed ONCE and shared), so it carries no
//   information and is not the headline -- max_j|defect_j| is.  This is the local analogue
//   of the PSA7 finding that a 1.6e-6 global balance sat on ~5.6e-5 of cancelling per-window
//   defects; DG would make each element exact, and this says what that is worth.
//
// Two in-situ validations of eval_colloc_deriv (0PDE5 flags it as "inherently ~one order
// coarser", so it must not be trusted blind): at z=0 the Danckwerts BC forces N = U*c_feed,
// known ANALYTICALLY, and at z=1 the outflow BC forces dc/dz = 0 so N = U*c(1,t).  Both are
// computed the hard way and differenced against the exact value.
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

#ifndef PSA8_OUT_PREFIX
#define PSA8_OUT_PREFIX "OCFE_PSA8"
#endif

// ---------- model: IDENTICAL to PSA5/6/7 (rung-9 / PDE36c) ----------
static double const U_VEL = 1.0, D_ax = 0.1, F_ph = 0.5, k_ldf = 2.0;
static double const qs1 = 1.0, qs2 = 1.0, b02_L = 1.0;
static double const beta1 = 2.0, beta2 = 1.0, T0_ref = 1.0;
static double const Cp_e = 1.0, G_cv = 1.0, lam = 0.1, dH1 = 1.0, dH2 = 0.5, hw = 0.5, Tw = 1.0;
static double const T_end = 5.0, t_step = 2.0;
static double const cA1 = 0.2, cB1 = 0.5, cA2 = 0.5, cB2 = 0.2;
static double const b01_nom = 3.0;
static double const ZONE_KW = 4.0;
static double const Z_LAYER = D_ax/U_VEL;          // Danckwerts inlet-layer thickness = 0.1

static inline double q1star_g( double a, double b, double p ){ return qs1*p*a/( 1.0 + p*a + b02_L*b ); }
static inline double q2star_g( double a, double b, double p ){ return qs2*b02_L*b/( 1.0 + p*a + b02_L*b ); }
static inline double lncosh( double x )
{ x = std::fabs(x); return x + std::log1p( std::exp(-2.0*x) ) - std::log(2.0); }

struct Feed {
  double eps = 0.0125;
  double stepf( double t ) const { return 0.5*( 1.0 + std::tanh( (t-t_step)/eps ) ); }
  double c1( double t ) const { return cA1 + (cB1-cA1)*stepf(t); }
  double c2( double t ) const { return cA2 + (cB2-cA2)*stepf(t); }
  double integ( int i, double t ) const {
    double const A=(i==1?cA1:cA2), B=(i==1?cB1-cA1:cB2-cA2);
    return A*t + B*0.5*( t + eps*( lncosh((t-t_step)/eps) - lncosh(t_step/eps) ) ); }
};

static double const GLX[10] = { -0.9739065285171717, -0.8650633666889845, -0.6794095682990244,
                                -0.4333953941292472, -0.1488743389816312,  0.1488743389816312,
                                 0.4333953941292472,  0.6794095682990244,  0.8650633666889845,
                                 0.9739065285171717 };
static double const GLW[10] = {  0.0666713443086881,  0.1494513491505806,  0.2190863625159820,
                                 0.2692667193099963,  0.2955242247147529,  0.2955242247147529,
                                 0.2692667193099963,  0.2190863625159820,  0.1494513491505806,
                                 0.0666713443086881 };

// ---------- PSA7's frozen t-mesh: G=2, geometric fill, budget 20 ----------
static double series( double h0, double rho, size_t n )
{ if( n==0 ) return 0.;
  if( std::fabs(rho-1.0) < 1e-12 ) return h0*double(n);
  return h0*rho*( std::pow(rho,double(n)) - 1.0 )/( rho - 1.0 ); }
static double solve_ratio( double h0, double L, size_t n )
{ if( n==0 || !(L>0.) ) return 1.0;
  if( series(h0,1.0,n) >= L ) return 1.0;
  double a=1.0,b=2.0; while( series(h0,b,n) < L && b < 1e8 ) b*=2.0;
  for( int i=0;i<200;++i ){ double const m=0.5*(a+b); if( series(h0,m,n)<L ) a=m; else b=m; }
  return 0.5*(a+b); }
static std::vector<double> geom_widths( double h0, double L, size_t n, double rho )
{ std::vector<double> w; if( n==0 ) return w;
  if( rho <= 1.0+1e-12 ){ w.assign(n,L/double(n)); return w; }
  for( size_t k=0;k<n;++k ) w.push_back( h0*std::pow(rho,double(k+1)) );
  double s=0.; for(double x:w) s+=x; for(double& x:w) x*=L/s; return w; }

static std::vector<double> t_mesh( double eps, size_t G, size_t budget )
{
  double const W=ZONE_KW*eps, lo=t_step-W, hi=t_step+W, h0=2.0*W/double(G);
  double const La=lo, Lc=T_end-hi;
  size_t const avail = ( budget>G+1 ? budget-G : 2 );
  size_t bL=1,bR=avail-1; double best=1e300, rL=1., rR=1.;
  for( size_t nl=1; nl<avail; ++nl ){
    size_t const nr=avail-nl;
    double const a=solve_ratio(h0,La,nl), b=solve_ratio(h0,Lc,nr);
    double const m=std::max(std::max(a,1.0/std::max(a,1e-30)),std::max(b,1.0/std::max(b,1e-30)));
    if( m<best ){ best=m; bL=nl; bR=nr; rL=a; rR=b; } }
  std::vector<double> const wl=geom_widths(h0,La,bL,rL), wr=geom_widths(h0,Lc,bR,rR);
  std::vector<double> b; b.push_back(0.);
  for( size_t i=wl.size(); i-->0; ) b.push_back( b.back()+wl[i] );
  for( size_t i=0;i<G;++i ) b.push_back( b.back()+h0 );
  for( size_t i=0;i<wr.size();++i ) b.push_back( b.back()+wr[i] );
  b.back()=T_end; return b;
}

//! z-mesh: uniform, or two-zone graded fine over the Danckwerts layer [0, 2*Z_LAYER]
static std::vector<double> z_mesh( size_t ne, size_t nfine )
{
  std::vector<double> b;
  if( nfine == 0 || nfine+1 > ne ){
    for( size_t i=0;i<=ne;++i ) b.push_back( double(i)/double(ne) );
    return b; }
  double const Wz = 2.0*Z_LAYER;
  for( size_t i=0;i<nfine;++i )   b.push_back( Wz*double(i)/double(nfine) );
  size_t const nc = ne-nfine;
  for( size_t i=0;i<=nc;++i )     b.push_back( Wz + (1.0-Wz)*double(i)/double(nc) );
  b.front()=0.; b.back()=1.; return b;
}

//! 4-point Richardson extrapolation to the face from ONE side (error O(d^4))
static inline double extrap4( double v0, double v1, double v2, double v3 )
{ return 4.0*v0 - 6.0*v1 + 4.0*v2 - v3; }

// ---------------------------------------------------------------------------
struct Cfg { size_t ne_z=5, nn_z=6, nfine_z=0; double eps=0.0125; std::string tag; };

struct SR {
  Cfg cfg;
  bool converged=false; int iters=0; size_t nVar=0, nwin=0;
  double mb_kpi1=1.;
  double under1=0., under2=0., zmin1=-1., tmin1=-1.;
  double over1=0.;
  double jump_c=0., jump_dz=0., jump_N=0.;          // face jumps at the SMALLEST delta
  double jump_c_r=0., jump_N_r=0.;                  // ratio jump(4d)/jump(d): ~4^4 => ->0
  double loc_max=0., loc_sum=0.;                    // per-z-element balance
  double bc0_err=0., bc1_err=0.;                    // eval_colloc_deriv validations
  std::vector<double> zd;                           // per-element defects
};

// ---------------------------------------------------------------------------
static SR run_case( Cfg const& cfg )
{
  SR R; R.cfg = cfg;
  size_t const n_nd_t = 6;
  FFDom::TYPE const coltype = FFDom::LGL;
  Feed fd; fd.eps = cfg.eps;

  std::vector<double> const t_bnd = t_mesh( cfg.eps, 2, 20 );      // PSA7 best: frozen
  std::vector<double> const z_bnd = z_mesh( cfg.ne_z, cfg.nfine_z );

  std::cout << "\n---- z " << cfg.ne_z << "x" << cfg.nn_z
            << ( cfg.nfine_z? "  GRADED" : "  uniform" )
            << "  eps=" << std::fixed << std::setprecision(5) << cfg.eps
            << "  (t-mesh frozen: G=2 geom ne=20)  " << cfg.tag << " ----\n     z-mesh:";
  for( double b : z_bnd ) std::cout << " " << std::setprecision(4) << b;
  std::cout << "\n";

  FFGraph DAG;
  FFVar t=DAG.add_var("t"), z=DAG.add_var("z");
  FFVar c1=DAG.add_var("c1(t,z)"), c2=DAG.add_var("c2(t,z)");
  FFVar q1=DAG.add_var("q1(t,z)"), q2=DAG.add_var("q2(t,z)"), T=DAG.add_var("T(t,z)");
  FFVar c1_ic=DAG.add_var("c1_ic(z)"), c2_ic=DAG.add_var("c2_ic(z)");
  FFVar q1_ic=DAG.add_var("q1_ic(z)"), q2_ic=DAG.add_var("q2_ic(z)");
  FFVar T_ic=DAG.add_var("T_ic(z)"), b01=DAG.add_var("b01");

  FFPartial  OpP;
  FFIntegral OpI;

  FFVar sstep  = 0.5*( 1.0 + tanh( ( t - t_step )/cfg.eps ) );
  FFVar c1feed = cA1 + ( cB1 - cA1 )*sstep;
  FFVar c2feed = cA2 + ( cB2 - cA2 )*sstep;

  FFVar b1 = b01  *exp( beta1*( 1.0/T - 1.0/T0_ref ) );
  FFVar b2 = b02_L*exp( beta2*( 1.0/T - 1.0/T0_ref ) );
  FFVar den = 1.0 + b1*c1 + b2*c2;
  FFVar q1star = qs1*b1*c1/den, q2star = qs2*b2*c2/den;

  FFVar CONT1 = OpP( c1, t ) + U_VEL*OpP( c1, z ) - D_ax*OpP( OpP( c1, z ), z ) + F_ph*OpP( q1, t );
  FFVar CONT2 = OpP( c2, t ) + U_VEL*OpP( c2, z ) - D_ax*OpP( OpP( c2, z ), z ) + F_ph*OpP( q2, t );
  FFVar LDF1  = OpP( q1, t ) - k_ldf*( q1star - q1 );
  FFVar LDF2  = OpP( q2, t ) - k_ldf*( q2star - q2 );
  FFVar ENE_T = Cp_e*OpP( T, t ) + G_cv*OpP( T, z ) - lam*OpP( OpP( T, z ), z )
              - F_ph*( dH1*OpP( q1, t ) + dH2*OpP( q2, t ) ) + hw*( T - Tw );
  FFVar IC_c1=c1-c1_ic, IC_c2=c2-c2_ic, IC_q1=q1-q1_ic, IC_q2=q2-q2_ic, IC_T=T-T_ic;
  FFVar BC_L1 = U_VEL*c1 - D_ax*OpP( c1, z ) - U_VEL*c1feed;
  FFVar BC_L2 = U_VEL*c2 - D_ax*OpP( c2, z ) - U_VEL*c2feed;
  FFVar BC_LT = G_cv*T  - lam*OpP( T, z ) - G_cv*T0_ref;
  FFVar BC_U1 = OpP( c1, z ), BC_U2 = OpP( c2, z ), BC_UT = OpP( T, z );
  FFVar Inv1 = OpI( c1 + F_ph*q1, z ), Inv2 = OpI( c2 + F_ph*q2, z );
  FFVar Eff1 = OpI( U_VEL*c1, t ),     Eff2 = OpI( U_VEL*c2, t );

  OCFESLV oc( &DAG );
  oc.add_domain( t, FFDom( t_bnd, coltype, n_nd_t ) );
  oc.add_domain( z, FFDom( z_bnd, coltype, cfg.nn_z ) );
  oc.add_state ( c1, {t,z} ); oc.add_state( c2, {t,z} );
  oc.add_state ( q1, {t,z} ); oc.add_state( q2, {t,z} ); oc.add_state( T, {t,z} );
  oc.add_input ( b01, b01_nom, /*is_decision=*/true );

  // update_ref LIFETIME: retained and invoked during init() -- FFVars by reference (they
  // outlive init() here), everything else BY VALUE.
  auto c1g = [&t,&z,fd]( OCFESLV::t_Coord const& cr ){ return fd.c1(cr.at(t))*(1.0-0.5*cr.at(z)); };
  auto c2g = [&t,&z,fd]( OCFESLV::t_Coord const& cr ){ return fd.c2(cr.at(t))*(1.0-0.5*cr.at(z)); };
  oc.update_ref( c1, c1g );
  oc.update_ref( c2, c2g );
  oc.update_ref( q1, [c1g,c2g]( OCFESLV::t_Coord const& cr ){ return q1star_g(c1g(cr),c2g(cr),b01_nom); } );
  oc.update_ref( q2, [c1g,c2g]( OCFESLV::t_Coord const& cr ){ return q2star_g(c1g(cr),c2g(cr),b01_nom); } );
  oc.update_ref( T,  [c1g,c2g]( OCFESLV::t_Coord const& cr ){
    return T0_ref + F_ph*( dH1*q1star_g(c1g(cr),c2g(cr),b01_nom)
                         + dH2*q2star_g(c1g(cr),c2g(cr),b01_nom) )/Cp_e; } );
  oc.add_input ( c1_ic, {z} ); oc.add_input( c2_ic, {z} );
  oc.add_input ( q1_ic, {z} ); oc.add_input( q2_ic, {z} );
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
  oc.add_output( Inv1, {t}, {T_end} ); oc.add_output( Inv2, {t}, {T_end} );
  oc.add_output( Eff1, {z}, {1.0}   ); oc.add_output( Eff2, {z}, {1.0}   );

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

  OCFESLV::SolveReport const rep = oc.solve( xv.data(), inp.data(), nullptr );
  R.converged = rep.converged; R.iters = rep.iterations;
  std::cout << "  nVar=" << R.nVar << " windows=" << R.nwin
            << " converged=" << (rep.converged?"yes":"no") << " iters=" << rep.iterations
            << " final|r|=" << std::scientific << std::setprecision(3) << rep.final_residual << "\n";
  if( !rep.converged ) return R;

  std::vector<double> const Fref = oc.val_functions();
  double const Feed1T = U_VEL*fd.integ(1,T_end);
  if( Fref.size() >= 4 )
    R.mb_kpi1 = std::fabs( Fref[0] + Fref[2] - Feed1T )/std::max(Feed1T,1e-30);

  auto V  = [&]( FFVar const& v, double zz, double tt ){
    OCFESLV::t_Coord pt; pt[z]=zz; pt[t]=tt;
    return oc.eval_colloc<double>( v, pt, xv.data(), inp.data(), nullptr ); };
  auto DZ = [&]( FFVar const& v, double zz, double tt ){
    OCFESLV::t_Coord pt; pt[z]=zz; pt[t]=tt;
    return oc.eval_colloc_deriv<double>( v, pt, z, 1, xv.data(), inp.data(), nullptr ); };
  auto Nflux = [&]( FFVar const& v, double zz, double tt ){ return U_VEL*V(v,zz,tt) - D_ax*DZ(v,zz,tt); };

  // ---- QUESTION A: undershoot WITH ITS LOCATION ----
  // two-pass: coarse global scan, then a local refine box around the argmin.  Cheaper than
  // a uniformly fine scan AND better localised -- which matters, since the whole point is
  // to find out WHICH feature the undershoot belongs to.
  { double m1=1e30, M1=-1e30, m2=1e30;
    auto scan = [&]( double z0, double z1, double s0, double s1, int NZ, int NT ){
      for( int iz=0; iz<=NZ; ++iz ){
        double const zz = z0 + (z1-z0)*double(iz)/double(NZ);
        if( zz<0. || zz>1. ) continue;
        for( int it=0; it<=NT; ++it ){
          double const tt = s0 + (s1-s0)*double(it)/double(NT);
          if( tt<0. || tt>T_end ) continue;
          double const v1=V(c1,zz,tt), v2=V(c2,zz,tt);
          if( v1 < m1 ){ m1=v1; R.zmin1=zz; R.tmin1=tt; }
          M1=std::max(M1,v1); m2=std::min(m2,v2); } } };
    scan( 0., 1., 0., T_end, 120, 240 );
    double const dz=2.0/120.0, dt=2.0*T_end/240.0;
    scan( R.zmin1-dz, R.zmin1+dz, R.tmin1-dt, R.tmin1+dt, 60, 60 );
    R.under1=std::min(0.0,m1); R.under2=std::min(0.0,m2);
    R.over1 =std::max(0.0,M1-std::max(cA1,cB1)); }

  // ---- eval_colloc_deriv validation against the two BCs ----
  { double e0=0., e1=0.;
    for( int k=0;k<=40;++k ){
      double const tt = T_end*double(k)/40.0;
      if( tt <= 0. ) continue;
      e0 = std::max( e0, std::fabs( Nflux(c1,0.0,tt) - U_VEL*fd.c1(tt) ) );
      e1 = std::max( e1, std::fabs( DZ(c1,1.0,tt) ) ); }
    R.bc0_err=e0; R.bc1_err=e1; }

  // ---- QUESTION B1: FACE JUMPS, one-sided limits formed independently ----
  { double const scale = std::max(cB1,cA1);
    for( int pass=0; pass<2; ++pass ){          // pass 0: delta ; pass 1: 4*delta
      double jc=0., jd=0., jn=0.;
      for( size_t f=1; f+1<z_bnd.size(); ++f ){
        double const hL = z_bnd[f]-z_bnd[f-1], hR = z_bnd[f+1]-z_bnd[f];
        double const d  = ( pass? 0.04 : 0.01 )*std::min(hL,hR);
        for( int k=0; k<=15; ++k ){
          double const tt = T_end*double(k)/15.0;
          if( tt<=0. ) continue;
          double cl[4],cr[4],dl[4],dr[4];
          for( int s=0;s<4;++s ){
            double const off = double(s+1)*d;
            cl[s]=V (c1,z_bnd[f]-off,tt); cr[s]=V (c1,z_bnd[f]+off,tt);
            dl[s]=DZ(c1,z_bnd[f]-off,tt); dr[s]=DZ(c1,z_bnd[f]+off,tt); }
          double const CL=extrap4(cl[0],cl[1],cl[2],cl[3]), CR=extrap4(cr[0],cr[1],cr[2],cr[3]);
          double const DL=extrap4(dl[0],dl[1],dl[2],dl[3]), DR=extrap4(dr[0],dr[1],dr[2],dr[3]);
          jc=std::max(jc,std::fabs(CR-CL));
          jd=std::max(jd,std::fabs(DR-DL));
          jn=std::max(jn,std::fabs( (U_VEL*CR-D_ax*DR) - (U_VEL*CL-D_ax*DL) )); } }
      if( pass==0 ){ R.jump_c=jc/scale; R.jump_dz=jd; R.jump_N=jn/scale; }
      else { R.jump_c_r = ( R.jump_c>0? (jc/scale)/R.jump_c : 0. );
             R.jump_N_r = ( R.jump_N>0? (jn/scale)/R.jump_N : 0. ); } } }

  // ---- QUESTION B2: PER-Z-ELEMENT BALANCE over the whole horizon ----
  { size_t const nz = z_bnd.size()-1;
    std::vector<double> face( nz+1, 0. );        // \int_0^T N(z_f,t) dt : computed ONCE per face
    for( size_t f=0; f<=nz; ++f ){
      double acc=0.;
      for( size_t k=0; k+1<t_bnd.size(); ++k ){
        double const a=t_bnd[k], b=t_bnd[k+1], hm=0.5*(b-a), mid=0.5*(a+b);
        for( int g=0; g<10; ++g ){
          double const tt=mid+hm*GLX[g], w=hm*GLW[g];
          acc += w*( f==0 ? U_VEL*fd.c1(tt)                       // Danckwerts BC: exact
                   : f==nz? U_VEL*V(c1,1.0,tt)                    // outflow BC: dc/dz = 0
                   :        Nflux(c1,z_bnd[f],tt) ); } }
      face[f]=acc; }
    R.zd.assign( nz, 0. );
    for( size_t j=0; j<nz; ++j ){
      double const a=z_bnd[j], b=z_bnd[j+1], hm=0.5*(b-a), mid=0.5*(a+b);
      double inv=0.;
      for( int g=0; g<10; ++g ){
        double const zz=mid+hm*GLX[g], w=hm*GLW[g];
        inv += w*( V(c1,zz,T_end) + F_ph*V(q1,zz,T_end) ); }
      R.zd[j] = ( inv - ( face[j] - face[j+1] ) )/std::max(Feed1T,1e-30);
      R.loc_max = std::max( R.loc_max, std::fabs(R.zd[j]) );
      R.loc_sum += std::fabs(R.zd[j]); } }

  std::cout << std::scientific << std::setprecision(3)
            << "  under1=" << R.under1 << " at (z*,t*)=(" << std::fixed << std::setprecision(4)
            << R.zmin1 << "," << R.tmin1 << ")   under2=" << std::scientific << R.under2
            << "  over1=" << R.over1 << "\n"
            << "  face jumps:  |[c]|=" << R.jump_c << " (ratio@4d " << std::fixed
            << std::setprecision(2) << R.jump_c_r << ")   |[dc/dz]|=" << std::scientific << R.jump_dz
            << "   |[N]|=" << R.jump_N << " (ratio@4d " << std::fixed << std::setprecision(2)
            << R.jump_N_r << ")\n"
            << "  per-element: max=" << std::scientific << std::setprecision(3) << R.loc_max
            << " sum|.|=" << R.loc_sum << "   mb_KPI1=" << R.mb_kpi1
            << "   BCcheck: z=0 " << R.bc0_err << "  z=1 " << R.bc1_err << "\n";
  return R;
}

// ---------------------------------------------------------------------------
int main()
{
  std::cout << "================================================================\n"
            << "  PSA8 -- settling the DG question in the direction 0r proposed it (z)\n"
            << "  A: is the undershoot h-convergent (Gibbs) or a plateau (operator)?\n"
            << "  B: are the z-face fluxes ALREADY single-valued (telescoping free)?\n"
            << "================================================================\n";

  std::vector<Cfg> cfgs;
  for( size_t ne : { size_t(5), size_t(8), size_t(12), size_t(16), size_t(24) } ){
    Cfg c; c.ne_z=ne; c.nn_z=6; c.tag="h-refine"; cfgs.push_back(c); }
  for( size_t nn : { size_t(8), size_t(10) } ){
    Cfg c; c.ne_z=5; c.nn_z=nn; c.tag="p-refine"; cfgs.push_back(c); }
  { Cfg c; c.ne_z=8; c.nn_z=6; c.nfine_z=3; c.tag="z-graded (inlet layer)"; cfgs.push_back(c); }
  { Cfg c; c.ne_z=5; c.nn_z=6; c.eps=0.05; c.tag="eps control (ties to PSA7)"; cfgs.push_back(c); }

  std::vector<SR> rows;
  for( Cfg const& c : cfgs ) rows.push_back( run_case(c) );

  std::cout << "\n============================== PSA8 TABLE ==============================\n"
            << std::left << std::setw(6) << "ne_z" << std::setw(6) << "nn_z"
            << std::setw(7) << "grad" << std::setw(9) << "eps" << std::setw(8) << "nVar"
            << std::right
            << std::setw(11) << "under1" << std::setw(9) << "z*" << std::setw(9) << "t*"
            << std::setw(11) << "|[c]|" << std::setw(11) << "|[N]|" << std::setw(7) << "r@4d"
            << std::setw(11) << "loc_max" << std::setw(11) << "mb_KPI1" << "\n";
  for( SR const& R : rows ){
    std::cout << std::left << std::setw(6) << R.cfg.ne_z << std::setw(6) << R.cfg.nn_z
              << std::setw(7) << (R.cfg.nfine_z? "yes":"no")
              << std::fixed << std::setprecision(4) << std::setw(9) << R.cfg.eps
              << std::setw(8) << R.nVar << std::right;
    if( !R.converged ){ std::cout << "   *** DIVERGED ***\n"; continue; }
    std::cout << std::scientific << std::setprecision(2) << std::setw(11) << R.under1
              << std::fixed << std::setprecision(4) << std::setw(9) << R.zmin1
              << std::setw(9) << R.tmin1
              << std::scientific << std::setprecision(2)
              << std::setw(11) << R.jump_c << std::setw(11) << R.jump_N
              << std::fixed << std::setprecision(2) << std::setw(7) << R.jump_N_r
              << std::scientific << std::setprecision(2)
              << std::setw(11) << R.loc_max << std::setw(11) << R.mb_kpi1 << "\n";
  }

  std::cout << "\n  READING THE TABLE\n"
               "  A -- under1 across the h-rows (ne_z 5..24 at nn_z=6) vs the p-rows (nn_z 8,10\n"
               "       at ne_z=5).  0r's Gibbs diagnosis predicts h CONVERGES and p only NARROWS\n"
               "       at fixed amplitude.  If h drives under1 to zero, the oscillation motivation\n"
               "       for DG is dead.  If BOTH plateau, it is the strong-form advection operator\n"
               "       and upwinding has a live case independent of conservation.\n"
               "       (z*,t*) says which feature it is: z* ~ 0 with t* small is the Danckwerts\n"
               "       inlet layer (thickness D/U = 0.1) -- then the z-graded row should kill it.\n"
               "  B -- r@4d is the decisive number.  The jump is re-measured at 4x the offset:\n"
               "       r ~ 256 means the 'jump' is pure extrapolation truncation, i.e. the flux is\n"
               "       CONTINUOUS and DG's telescoping is ALREADY satisfied by IC_STRONG on c and\n"
               "       on the Dz_c order-reduction auxiliary -- option 2 then adds nothing at the\n"
               "       face level and only makes the PER-ELEMENT balance exact.\n"
               "       r ~ 1 means the jump is a genuine limit: the flux is double-valued, the sum\n"
               "       does NOT telescope, and option 2 recovers its original justification.\n"
               "       NOISE FLOOR: extrap4 amplifies roundoff ~15x, so jumps below ~1e-13 are noise\n"
               "       and r is then erratic rather than ~256 -- that ALSO reads as 'continuous'.\n"
               "  B2 - loc_max is what DG would zero.  Compare it against mb_KPI1: if loc_max is\n"
               "       orders ABOVE the global balance, the elements are trading mass and only the\n"
               "       total is accurate -- exactly the PSA7 cancellation picture, but in z.\n"
               "  Validations: BCcheck z=0 (Danckwerts flux vs analytic) and z=1 (dc/dz vs 0) bound\n"
               "  the eval_colloc_deriv error; if they are not small, the flux numbers are not\n"
               "  trustworthy and nothing in column B means anything.\n";

  std::string const pre = std::string(PSA8_OUT_PREFIX);
  { std::ofstream f( pre+"_summary.out" );
    f << "# ne_z nn_z graded eps nVar iters under1 z_min t_min under2 over1 jump_c jump_dz jump_N "
         "jump_c_r jump_N_r loc_max loc_sum mb_KPI1 bc0 bc1 conv\n";
    for( SR const& R : rows )
      f << R.cfg.ne_z << " " << R.cfg.nn_z << " " << (R.cfg.nfine_z?1:0) << " "
        << std::setprecision(8) << R.cfg.eps << " " << R.nVar << " " << R.iters << " "
        << R.under1 << " " << R.zmin1 << " " << R.tmin1 << " " << R.under2 << " " << R.over1 << " "
        << R.jump_c << " " << R.jump_dz << " " << R.jump_N << " "
        << R.jump_c_r << " " << R.jump_N_r << " " << R.loc_max << " " << R.loc_sum << " "
        << R.mb_kpi1 << " " << R.bc0_err << " " << R.bc1_err << " " << (R.converged?1:0) << "\n"; }

  { std::ofstream f( pre+"_zdefect.out" );
    f << "# per-z-element balance defect; one block per config\n";
    for( SR const& R : rows ){
      if( !R.converged ) continue;
      std::vector<double> const zb = z_mesh( R.cfg.ne_z, R.cfg.nfine_z );
      f << "\n\n# ne_z=" << R.cfg.ne_z << " nn_z=" << R.cfg.nn_z
        << " graded=" << (R.cfg.nfine_z?1:0) << " eps=" << R.cfg.eps << "\n# z_lo z_hi defect\n";
      for( size_t j=0; j<R.zd.size(); ++j )
        f << std::setprecision(8) << zb[j] << " " << zb[j+1] << " " << R.zd[j] << "\n"; } }

  { std::ofstream gp( pre+".gp" );
    gp << "set datafile commentschars '#'\n"
          "set terminal pngcairo size 1000,700 enhanced font 'Helvetica,12'\n"
          "set grid\n\n"
          "set output '" << pre << "_under.png'\n"
          "set title 'undershoot vs z-resolution: h-convergence (Gibbs) or plateau (operator)?'\n"
          "set xlabel 'ne_z (h-refinement, nn_z=6)'; set ylabel '|under1|'; set logscale y\n"
          "plot '" << pre << "_summary.out' u ($2==6&&$3==0?$1:1/0):($7<0?-$7:1e-16) "
          "w lp lw 2 pt 7 t 'h-refinement'\n"
          "unset logscale\n\n"
          "set output '" << pre << "_zdefect.png'\n"
          "set title 'per-z-element balance defect (what DG would zero)'\n"
          "set xlabel 'z'; set ylabel 'defect'; set key outside\n"
          "plot for [i=0:4] '" << pre << "_zdefect.out' index i u 1:3 w steps lw 2 t 'block '.i\n"; }

  bool all=true; for( SR const& R : rows ) all &= R.converged;
  std::cout << "\n  wrote " << pre << "_{summary,zdefect}.out and " << pre << ".gp\n"
            << "  Overall: " << ( all ? "ALL CONVERGED" : "SOME DID NOT CONVERGE" ) << "\n";
  return all ? 0 : 1;
}
