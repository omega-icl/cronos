// OCFE_PDE38_solve2.cpp  ---  PSA rung 9c: inlet-sharpness sweep on the 7-STATE INDEX-1 structure
// ===========================================================================
// The rung-9b Heaviside sweep, but on the FULL index-1 momentum structure of
// rungs 6-8: velocity u slaved by DARCY to the gradient of a value-slaved total
// pressure P = Rg T (C_base + c1 + c2).  An inert carrier C_base keeps P > 0 on the
// regenerated bed (else P=0 -> u ill-defined).  States: c1,c2,q1,q2,T,u,P (7);
// probe should read  5 dyn {c1,c2,q1,q2,T} / 5 alg {u,P,Dz_c1,Dz_c2,Dz_T},
// INDEX 1 [parabolic], with the value-slaved drop flagging P (and u kept).
//
//   mass_i: dci/dt + d(u ci)/dz - D d2ci/dz2 + F dqi/dt              = 0
//   LDF_i:  dqi/dt - k ( qi*(c1,c2,T) - qi )                         = 0
//   energy: Cp dT/dt + G dT/dz - lam d2T/dz2 - F(dH1 dq1/dt+dH2 dq2/dt)+hw(T-Tw)=0
//   MOMENTUM (default DARCY):     (u - U0) + (kappa/mu) dP/dz = 0      [u = U0 - kappa dP/dz]
//     (-DMC__ERGUN toggle):       dP/dz + a(u-U0) + b(u-U0)|u-U0| = 0  [stays INDEX 1]
//   EOS (value-slaved):           P - Rg T ( C_base + c1 + c2 ) = 0
//
// PHYSICAL REGIME.  P is value-slaved by the EOS; u is the algebraic Darcy velocity about a base
// throughput U0:  u = U0 - kappa dP/dz.  The base flow U0 removes the t=0 startup singularity of the
// pure front-driven closure (uniform bed -> dP/dz=0 -> u=0 -> the sharp feed has no carrier velocity,
// which stalls Newton at the inlet corner); with U0>0 the bed has a sustained throughput that the
// EOS-pressure gradient modulates (~30-50% across the front at kappa=0.1).  Set U0_base=0 to recover
// the pure front-driven regime.  The index-1 structure is unchanged either way (u derivative-defined
// through dP/dz -> kept; P value-slaved -> drop fires).  A genuinely sustained COMPRESSIBLE flow with
// velocity from total continuity + EOS-density is INDEX 2 (Pantelides + dummy-derivatives -- separate
// milestone).  Ergun keeps the index at 1 in THIS closure; it only adds nonlinearity (toggle to check).
//
// Sweep tau in {0.2,0.1,0.05,0.025} (feed ci=c0_i(1-e^{-t/tau}) -> Heaviside), CGL 6x8
// (finest near-0 gap ~0.041: 0.2/0.1 resolved, 0.05 marginal, 0.025 under-resolved),
// tau-continuation, plus a finer-grid recovery run at the sharpest tau.  Output: gnuplot
// data (feed/inlet-c2/outlet-c2/inlet-velocity/summary) + OCFE_PDE38.gp.
// ===========================================================================

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>   // std::sort, std::max/min -- ranked coincidence groups
#include <limits>      // std::numeric_limits -- NaN order marker
#include <map>         // coincidence grouping by quantised coordinate
#include <utility>     // std::pair
#include <string>

#include "fflin.hpp"
// Configurable header guard.  A hardcoded #include "ocfeslv.hpp" silently ignores
// -DOCFE_OCFESLV_HEADER, so the driver compiles against whatever ocfeslv.hpp happens to be
// while the build log claims otherwise.  This is a documented trap in this corpus.
#ifndef OCFE_OCFESLV_HEADER
#define OCFE_OCFESLV_HEADER "ocfeslv.hpp"
#endif
#include OCFE_OCFESLV_HEADER

using namespace mc;

#define MC__ERGUN

#ifndef PDE38_OUT_PREFIX
#define PDE38_OUT_PREFIX "OCFE_PDE38"
#endif

// momentum / pressure
static double const Kperm = 0.1, Rg = 1.0, C_base = 1.0;   // inert carrier -> P >= Rg T0 C_base > 0
static double const U0_base = 1.0;   // base (pressure-driven) throughput: u = U0 - Kperm dP/dz.
// U0_base removes the t=0 startup singularity of the pure front-driven closure (uniform bed -> dP/dz=0
// -> u=0 -> feed cannot be carried in). Set U0_base=0 to recover the pure front-driven regime. The
// index-1 structure is unchanged: u is still derivative-defined through dP/dz (kept), P value-slaved.
#if defined(MC__ERGUN)
static double const Avisc = 1.0/Kperm, Bins = 1.0, ureg = 1.0e-3;  // Ergun (Darcy when Bins=0)
#endif
// transport / kinetics
static double const D_ax = 0.1, F_ph = 0.5, k_ldf = 2.0;
static double const qs1 = 1.0, qs2 = 1.0, b01_L = 3.0, b02_L = 1.0;
static double const beta1 = 2.0, beta2 = 1.0, T0_ref = 1.0;
static double const Cp_e = 1.0, G_cv = 1.0, lam = 0.1, dH1 = 1.0, dH2 = 0.5, hw = 0.5, Tw = 1.0;
static double const c0_1 = 0.5, c0_2 = 0.5, T_end = 5.0;

static size_t const NTS = 401;
static inline double tsamp( size_t kk ){ return double(kk)/double(NTS-1)*T_end; }
static inline double q1star_T0( double c1g, double c2g )
{ return qs1*b01_L*c1g/( 1.0 + b01_L*c1g + b02_L*c2g ); }
static inline double q2star_T0( double c1g, double c2g )
{ return qs2*b02_L*c2g/( 1.0 + b01_L*c1g + b02_L*c2g ); }

static size_t g_nonconv = 0;   // 2026-09-09: non-converged solves in cells NOT fed into `rows`

struct SR {
  double tau=0.; size_t n_el=0, n_nd=0; bool converged=false; int iters=0; double final_r=0.;
  double mb1=1., mb2=1., overshoot=0., undershoot=0., umax=0.;
  std::vector<double> c2out, c2inlet, uinlet;
  // --- imposition-mode comparison -----------------------------------------------------
  OCFESLV::Options::ImpositionType imp = OCFESLV::Options::IC_STRONG;
  char const* imp_name = "IC_STRONG";
  double dup_worst = -1.;      //!< worst per-state duplicate-node spread over ALL states
  double dup_u     = -1.;      //!< spread of u(t,z) alone -- the state under investigation
  double dup_claimed_max = -1.;//!< worst spread over the CLAIMED primitive states c1,c2,T,P
};

// The three modes, swept for every tau.
// ---------------------------------------------------------------------------------------
//  Where does the spread live?  --  t-seam vs z-interface
// ---------------------------------------------------------------------------------------
// u(t,z) is the ONLY claimed state in this model that does not reach round-off:
//   IC_WEAK 3.449e-06 / IC_TRACE 1.570e-08 / IC_STRONG 1.570e-08, while c1, c2, T, P,
//   Dz_c1, Dz_c2 and Dz_T are all EXACTLY 0.
// The exact modes improve it ~220x over weak, so the claim is being ENFORCED and merely
// falls short -- unlike the scalar_nest case, which turned out to be a reader bug and where
// the modes agreed to seven digits precisely because nothing was being measured.
//
// TWO CANDIDATE OWNERS, and they are not the same team:
//
//   (a) z-INTERFACE.  u is the only ALGEBRAIC claimed state: u = U0 - Kperm dP/dz (ERGUN
//       here, so nonlinear in u).  P is C0 but NOT C1, so dP/dz genuinely jumps across an
//       element interface and the momentum row at each of the two coincident nodes solves
//       for u from its own one-sided gradient.  The continuity claim and the momentum row
//       then pull against each other and 1.6e-08 is the compromise.  If so the spread is
//       TRUNCATION-limited and must fall under z-refinement -- not an interface-plan defect.
//       Consistent with P itself being exactly 0: P = Rg T (C_base+c1+c2) is a VALUE function
//       of C0 states and inherits continuity exactly.
//
//   (b) t-SEAM.  This driver MARCHES (six windows).  A window boundary is an IC TRANSFER,
//       not a trace claim: differential states transfer exactly, an algebraic state is
//       recomputed in the new window from its own equation.  If the spread sits only at
//       t-seams then the interface plan is not involved at all and this belongs to the
//       moving-mesh-under-marching line (roadmap C3), which flags exactly this gap.
//
// Printing the COORDINATES of the worst group separates them in one run.  Do not infer the
// answer from the magnitude; infer it from where the nodes are.
static void dump_worst_group( OCFESLV const& oc, FFVar const& V, std::vector<double> const& xv,
                              size_t ne_t, double t_end, size_t ne_z, double z_end,
                              char const* mode_name, size_t max_groups = 6 )
{
  // pos_state() does ndx_el.at(d) for EVERY domain of the state, so an EMPTY index map throws
  // for a distributed state -- pin every domain to element 0 instead.  (Two revisions were
  // lost to the empty-map idiom, which is valid only for scalars.)
  std::map<FFVar,size_t,lt_FFVar> ndx0;
  {
    auto const& vs = oc.var_state();
    auto const it = vs.find( V );
    if( it != vs.end() ) for( auto const& d : it->second ) ndx0[d] = 0;
  }
  size_t const off = oc.pos_state( V, ndx0 );

  std::vector<std::vector<double>> const nodes = oc.node_colloc( V );
  std::map<std::vector<long long>,std::vector<size_t>> groups;
  for( size_t i=0; i<nodes.size(); ++i ){
    std::vector<long long> key;
    for( double c : nodes[i] ) key.push_back( (long long)std::llround( c * 1.0e12 ) );
    groups[key].push_back( i );
  }

  // Rank coincidence groups by spread so the worst offenders print first.
  std::vector<std::pair<double,std::vector<size_t>>> ranked;
  for( auto const& kv : groups ){
    if( kv.second.size() < 2 ) continue;
    double lo = 1e300, hi = -1e300;
    for( size_t i : kv.second ){
      if( off+i >= xv.size() ) continue;
      lo = std::min( lo, xv[off+i] ); hi = std::max( hi, xv[off+i] );
    }
    ranked.emplace_back( hi-lo, kv.second );
  }
  std::sort( ranked.begin(), ranked.end(),
             []( auto const& A, auto const& B ){ return A.first > B.first; } );

  std::cout << "  [where] " << mode_name << "  " << V.name()
            << "  base_off=" << off << "  n_node=" << nodes.size()
            << "  coincidence_groups=" << ranked.size() << "\n";

  // Attribute each group by COORDINATE against the known element boundaries, not by index
  // stride: the flat ordering is an implementation detail and guessing at it is how the last
  // three instrument bugs happened.  Elements are uniform, so the interior boundaries are
  // exactly k*t_end/ne_t and k*z_end/ne_z.
  auto on_interior_boundary = []( double c, size_t ne, double span ){
    if( ne < 2 ) return false;
    for( size_t k = 1; k < ne; ++k )
      if( std::fabs( c - double(k)*span/double(ne) ) < 1e-9*std::max(1.0,span) ) return true;
    return false;
  };
  size_t n_t = 0, n_z = 0, n_both = 0, n_neither = 0;
  double worst_t = 0., worst_z = 0., worst_both = 0.;
  for( auto const& r : ranked ){
    if( r.second.empty() ) continue;
    std::vector<double> const& c0 = nodes[ r.second.front() ];
    bool const bt = ( c0.size() > 0 ) && on_interior_boundary( c0[0], ne_t, t_end );
    bool const bz = ( c0.size() > 1 ) && on_interior_boundary( c0[1], ne_z, z_end );
    if     ( bt && bz ){ ++n_both;    worst_both = std::max( worst_both, r.first ); }
    else if( bt )      { ++n_t;       worst_t    = std::max( worst_t,    r.first ); }
    else if( bz )      { ++n_z;       worst_z    = std::max( worst_z,    r.first ); }
    else                 ++n_neither;
  }
  std::cout << "  [where]   t-seam only groups=" << n_t     << " worst=" << std::scientific
            << std::setprecision(3) << worst_t    << "\n"
            << "  [where]   z-face only groups=" << n_z     << " worst=" << worst_z    << "\n"
            << "  [where]   corner (both)  groups=" << n_both << " worst=" << worst_both << "\n"
            << "  [where]   unattributed   groups=" << n_neither
            << "   (nonzero here means the coordinate ordering is not (t,z) -- check before"
               " reading the rest)\n";

  size_t shown = 0;
  for( auto const& r : ranked ){
    if( shown++ >= max_groups ) break;
    std::cout << "  [where]   spread=" << std::scientific << std::setprecision(6) << r.first
              << "  members=" << r.second.size() << "\n";
    for( size_t i : r.second ){
      double const val = ( off+i < xv.size() ) ? xv[off+i]
                                               : std::numeric_limits<double>::quiet_NaN();
      std::cout << "  [where]     i=" << std::setw(5) << i
                << " var[" << std::setw(6) << (off+i) << "]  coord=(";
      for( size_t c=0; c<nodes[i].size(); ++c )
        std::cout << ( c? "," : "" ) << std::fixed << std::setprecision(6) << nodes[i][c];
      std::cout << ")  value=" << std::scientific << std::setprecision(10) << val << "\n";
    }
  }
}

struct ImpMode { OCFESLV::Options::ImpositionType t; char const* name; };
static ImpMode const IMP_MODES[3] = {
  { OCFESLV::Options::IC_WEAK,   "IC_WEAK"   },
  { OCFESLV::Options::IC_TRACE,  "IC_TRACE"  },
  { OCFESLV::Options::IC_STRONG, "IC_STRONG" }
};

static SR run_sharp( double tau, size_t ne_t, size_t nn_t, size_t ne_z, size_t nn_z,
                     std::vector<double>& warm, std::string const& label,
                     ImpMode const& mode = IMP_MODES[2], bool want_where = false )
{
  SR R; R.tau=tau; R.n_el=ne_t; R.n_nd=nn_t;   // n_el/n_nd hold the t-grid
  R.imp = mode.t; R.imp_name = mode.name;
  std::cout << "\n---- tau=" << std::fixed << std::setprecision(5) << tau
            << "  t-grid " << ne_t << "x" << nn_t << " z-grid " << ne_z << "x" << nn_z
            << "  imposition=" << mode.name
            << "  " << label << " ----\n";

  FFGraph DAG;
  FFVar t  = DAG.add_var( "t" );
  FFVar z  = DAG.add_var( "z" );
  FFVar c1 = DAG.add_var( "c1(t,z)" );
  FFVar c2 = DAG.add_var( "c2(t,z)" );
  FFVar q1 = DAG.add_var( "q1(t,z)" );
  FFVar q2 = DAG.add_var( "q2(t,z)" );
  FFVar T  = DAG.add_var( "T(t,z)" );
  FFVar u  = DAG.add_var( "u(t,z)" );
  FFVar P  = DAG.add_var( "P(t,z)" );
  std::vector<FFVar> c12{ c1, c2 };

  FFPartial  OpP;
  FFIntegral OpI;
  FFLin      OpSum;

  FFVar c1feed = c0_1*( 1.0 - exp( -t/tau ) );
  FFVar c2feed = c0_2*( 1.0 - exp( -t/tau ) );
  FFVar b1 = b01_L*exp( beta1*( 1.0/T - 1.0/T0_ref ) );
  FFVar b2 = b02_L*exp( beta2*( 1.0/T - 1.0/T0_ref ) );
  FFVar den = 1.0 + b1*c1 + b2*c2;
  FFVar q1star = qs1*b1*c1/den;
  FFVar q2star = qs2*b2*c2/den;

  FFVar CONT1 = OpP( c1, t ) + OpP( u*c1, z ) - D_ax*OpP( OpP( c1, z ), z ) + F_ph*OpP( q1, t );
  FFVar CONT2 = OpP( c2, t ) + OpP( u*c2, z ) - D_ax*OpP( OpP( c2, z ), z ) + F_ph*OpP( q2, t );
  FFVar LDF1  = OpP( q1, t ) - k_ldf*( q1star - q1 );
  FFVar LDF2  = OpP( q2, t ) - k_ldf*( q2star - q2 );
  FFVar ENE_T = Cp_e*OpP( T, t ) + G_cv*OpP( T, z ) - lam*OpP( OpP( T, z ), z )
              - F_ph*( dH1*OpP( q1, t ) + dH2*OpP( q2, t ) ) + hw*( T - Tw );
#if defined(MC__ERGUN)
  FFVar du = u - U0_base;
  FFVar MOM = OpP( P, z ) + Avisc*du + Bins*du*sqrt( du*du + ureg*ureg ); // Ergun on (u-U0) (INDEX 1)
#else
  FFVar MOM = ( u - U0_base ) + Kperm*OpP( P, z );                        // Darcy: u = U0 - Kperm dP/dz
#endif
  FFVar EOS = P - Rg*T*( OpSum( c12, 1., C_base ) );                            // value-slaved total pressure
  //FFVar EOS = P - Rg*T*( C_base + c1 + c2 );                            // value-slaved total pressure
  FFVar IC_c1 = c1, IC_c2 = c2, IC_q1 = q1, IC_q2 = q2, IC_T = T - T0_ref;
  FFVar BC_L1 = u*c1 - D_ax*OpP( c1, z ) - u*c1feed;     // Danckwerts inlet (variable u)
  FFVar BC_L2 = u*c2 - D_ax*OpP( c2, z ) - u*c2feed;
  FFVar BC_LT = G_cv*T - lam*OpP( T, z ) - G_cv*T0_ref;
  FFVar BC_U1 = OpP( c1, z ), BC_U2 = OpP( c2, z ), BC_UT = OpP( T, z );

  FFVar Inv1 = OpI( c1 + F_ph*q1, z ), Inv2 = OpI( c2 + F_ph*q2, z );
  // feed/effluent fluxes (variable u) are computed NUMERICALLY in post-processing, NOT as OpI-over-t
  // outputs: OpI(u*c,t) / OpI(u*c_feed,t) leave a residual t-dependence the output checker rejects
  // ("DOMAIN VARIABLE t USED IN OUTPUT ... NO EVALUATION VALUE"). Inv (OpI-over-z) is unaffected.

  OCFESLV oc( &DAG );
  oc.add_domain( t, FFDom( 0., T_end, ne_t, FFDom::CGL, nn_t ) );
  oc.add_domain( z, FFDom( 0., 1.0,   ne_z, FFDom::CGL, nn_z ) );
  oc.add_state ( c1, {t,z} ); oc.add_state( c2, {t,z} );
  oc.add_state ( q1, {t,z} ); oc.add_state( q2, {t,z} ); oc.add_state( T, {t,z} );
  oc.add_state ( u,  {t,z} ); oc.add_state( P,  {t,z} );
  oc.update_ref( c1, [&]( OCFESLV::t_Coord const& cr ){
    double zz=cr.at(z), tt=cr.at(t); return c0_1*(1.0-std::exp(-tt/tau))*(1.0-0.5*zz); } );
  oc.update_ref( c2, [&]( OCFESLV::t_Coord const& cr ){
    double zz=cr.at(z), tt=cr.at(t); return c0_2*(1.0-std::exp(-tt/tau))*(1.0-0.5*zz); } );
  oc.update_ref( q1, [&]( OCFESLV::t_Coord const& cr ){
    double zz=cr.at(z), tt=cr.at(t); double f=(1.0-std::exp(-tt/tau))*(1.0-0.5*zz);
    return q1star_T0( c0_1*f, c0_2*f ); } );
  oc.update_ref( q2, [&]( OCFESLV::t_Coord const& cr ){
    double zz=cr.at(z), tt=cr.at(t); double f=(1.0-std::exp(-tt/tau))*(1.0-0.5*zz);
    return q2star_T0( c0_1*f, c0_2*f ); } );
  oc.update_ref( T,  [&]( OCFESLV::t_Coord const& cr ){
    double zz=cr.at(z), tt=cr.at(t); double f=(1.0-std::exp(-tt/tau))*(1.0-0.5*zz);
    return T0_ref + F_ph*( dH1*q1star_T0(c0_1*f,c0_2*f) + dH2*q2star_T0(c0_1*f,c0_2*f) )/Cp_e; } );
  oc.update_ref( u,  [&]( OCFESLV::t_Coord const& cr ){
    double tt=cr.at(t); return U0_base + 0.5*Kperm*Rg*(c0_1+c0_2)*(1.0-std::exp(-tt/tau)); } ); // U0 + Darcy estimate
  oc.update_ref( P,  [&]( OCFESLV::t_Coord const& cr ){
    double zz=cr.at(z), tt=cr.at(t); double f=(1.0-std::exp(-tt/tau))*(1.0-0.5*zz);
    double Tg=T0_ref+F_ph*(dH1*q1star_T0(c0_1*f,c0_2*f)+dH2*q2star_T0(c0_1*f,c0_2*f))/Cp_e;
    return Rg*Tg*( C_base + c0_1*f + c0_2*f ); } );
  oc.set_evolution_domain( t );

  int const T_NO_LB = FFDom::ALL - FFDom::LB;
  int const Z_INT   = FFDom::ALL - FFDom::LB - FFDom::UB;
  oc.add_equation( CONT1, {t,z}, {T_NO_LB,   Z_INT},      OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( CONT2, {t,z}, {T_NO_LB,   Z_INT},      OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( LDF1,  {t,z}, {T_NO_LB,   FFDom::ALL}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( LDF2,  {t,z}, {T_NO_LB,   FFDom::ALL}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( ENE_T, {t,z}, {T_NO_LB,   Z_INT},      OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( MOM,   {t,z}, {FFDom::ALL,FFDom::ALL}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( EOS,   {t,z}, {FFDom::ALL,FFDom::ALL}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
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
  oc.add_output( Inv1,   {t}, {T_end} );  oc.add_output( Inv2,   {t}, {T_end} );

  oc.options.REDUCE.ORDER     = OCFESLV::Options::RED_FULL;
  oc.options.CLASSIFY.MODE         = OCFESLV::Options::CLASS_AUTO;
  oc.options.INTERFACE.TYPE   = OCFESLV::Options::IC_AUTO;
  oc.options.INTERFACE.IMPOSITION  = mode.t;
  oc.options.INTERFACE.SAT_SIGMA0       = 10.0;
  oc.options.DISPLAY_LEVEL    = ( label.rfind("first",0)==0 ? 1 : 0 );   // show probe once
#if defined(CRONOS__WITH_SPQR)
  oc.options.SOLVE.FACTORIZATION = OCFESLV::Options::SOLVE_SPQR;
#else
  oc.options.SOLVE.FACTORIZATION = OCFESLV::Options::SOLVE_SUPERLU;
#endif

  if( !oc.setup() ){ std::cerr << "ERROR: setup failed\n"; return R; }
  size_t const nVar=oc.n_colloc_sta(), nEqn=oc.n_colloc_eqn(), nFct=oc.n_colloc_fct();
  if( nVar!=nEqn ){ std::cerr << "ERROR: not square (" << nVar << "/" << nEqn << ")\n"; return R; }

  std::vector<double> varInit, inpInit;
  if( !oc.init( varInit, inpInit, nullptr ) ){ std::cerr << "ERROR: init failed\n"; return R; }
  std::vector<double> xv = ( warm.size()==varInit.size() ? warm : varInit );
  OCFESLV::SolveReport rep = oc.solve( xv.data(), nullptr, nullptr );
  R.converged=rep.converged;
  // 2026-09-09: MEASURED -- with the claim drop unlocked this driver reports 7 non-converged
  // marching windows and still exits 0, because only the tau-sweep `rows` reach the verdict:
  // the MODE-CMP and Z-REFINE cells are printed and discarded.  Count every cell here.
  if( !R.converged ) ++g_nonconv; R.iters=rep.iterations; R.final_r=rep.final_residual;
  std::cout << "  nVar=" << nVar << " converged=" << (rep.converged?"yes":"no")
            << " iters=" << rep.iterations << " final|r|=" << std::scientific << std::setprecision(3)
            << rep.final_residual << "\n";

  // ---- duplicate-node spreads, per state, straight from the header -------------------
  // Called explicitly rather than relying on CRONOS_DUP_SPREAD so the comparison is part
  // of the driver's own output and cannot be lost by a sweep-script parsing change.
  //
  // WHAT THIS IS TESTING.  The rev43 corpus sweep found u(t,z) at 1.6e-08..3.0e-08 under
  // IC_STRONG with conv=yes, while EVERY other claimed state sat at EXACTLY 0.  u carries a
  // retained continuity claim ([claims] states WITH a retained claim: c1 c2 T u P Dz_c1
  // Dz_c2 Dz_T), so this is not a state that was never claimed -- it is a claim that is not
  // being delivered.  u is also the only state defined THROUGH a derivative of a
  // value-slaved quantity: u = U0 - Kperm dP/dz with P = Rg T (C_base+c1+c2).
  //
  // The discriminating question is whether the same claim is honoured under IC_WEAK and
  // IC_TRACE.  If u is clean under TRACE and dirty under STRONG, the defect is in the
  // Schur/keep-explicit path and it is Stage 2 work.  If it is dirty in all three, the
  // claim generation for derivative-defined states predates this arc entirely.
  {
    std::string const lab = std::string("PDE38 tau=") + std::to_string(tau) + " " + mode.name;
    R.dup_worst = oc.report_duplicate_node_spreads( xv, lab.c_str(), &std::cout );

    auto sp = [&]( FFVar const& V ){
      size_t off = 0;
      for( auto const& st : oc.states_colloc() ){
        if( st.id() == V.id() )
          return oc.duplicate_node_spread( st, xv, off ).max_pair_spread;
        off += oc.node_colloc( st ).size();
      }
      return -1.0;
    };
    R.dup_u = sp( u );
    R.dup_claimed_max = std::max( std::max( sp(c1), sp(c2) ), std::max( sp(T), sp(P) ) );
    std::cout << "  [mode-cmp] " << mode.name
              << "  u=" << std::scientific << std::setprecision(3) << R.dup_u
              << "  max(c1,c2,T,P)=" << R.dup_claimed_max
              << "  worst_any=" << R.dup_worst << "\n";

    // Locate the u spread.  Also dump c1 as a CONTROL: it is claimed and at exactly 0, so
    // its groups must attribute the same way while showing no spread anywhere.  If c1's
    // attribution looks different from u's, the attribution itself is suspect.
    if( want_where ){
      dump_worst_group( oc, u,  xv, ne_t, T_end, ne_z, 1.0, mode.name );
      dump_worst_group( oc, c1, xv, ne_t, T_end, ne_z, 1.0, mode.name, 2 );
    }
  }

  if( !rep.converged ) return R;
  if( warm.size()==varInit.size() || warm.empty() ) warm = xv;

  auto eval = [&]( FFVar const& V, double zz, double tt ){
    OCFESLV::t_Coord pt; pt[z]=zz; pt[t]=tt;
    return oc.eval_colloc<double>( V, pt, xv.data(), nullptr, nullptr ); };

  R.c2out.resize(NTS); R.c2inlet.resize(NTS); R.uinlet.resize(NTS);
  for( size_t kk=0; kk<NTS; ++kk ){ double tt=tsamp(kk);
    R.c2out[kk]=eval(c2,1.0,tt); R.c2inlet[kk]=eval(c2,0.0,tt); R.uinlet[kk]=eval(u,0.0,tt); }

  double cmax=-1e30, cmin=1e30; R.umax=-1e30;
  int const NZ=60, NTf=600;
  for( int iz=0; iz<=NZ; ++iz ){ double zz=double(iz)/NZ;
    for( int it=0; it<=NTf; ++it ){ double tt=double(it)/NTf*1.0;
      double v1=eval(c1,zz,tt), v2=eval(c2,zz,tt);
      cmax=std::max(cmax,std::max(v1,v2)); cmin=std::min(cmin,std::min(v1,v2)); } }
  for( int iz=0; iz<=40; ++iz ) for( int it=0; it<=40; ++it )
    R.umax=std::max(R.umax, eval(u,double(iz)/40.0,double(it)/40.0*T_end));
  R.overshoot=cmax-std::max(c0_1,c0_2); R.undershoot=cmin;

  std::vector<double> res(nEqn,0.), fct(nFct,0.);
  if( !oc.eval( res.data(), fct.data(), xv.data(), nullptr, nullptr ) ){ std::cerr << "ERROR: output eval\n"; return R; }
  double inv1=fct[oc.row_fct(0)-nEqn], inv2=fct[oc.row_fct(1)-nEqn];     // Inv_i = int_z (c_i+F q_i)|_Tend
  // variable-u feed/effluent fluxes by trapezoidal quadrature over the t-sample grid:
  //   FeedF_i = int_0^Tend u(0,t) c_feed,i(t) dt   (inlet Danckwerts flux = u c_feed)
  //   EfflF_i = int_0^Tend u(1,t) c_i(1,t) dt      (outlet; dispersive flux ~0 since d_z c=0 there)
  double ff1=0.,ff2=0.,ef1=0.,ef2=0., dt=T_end/double(NTS-1);
  for( size_t kk=0; kk<NTS; ++kk ){
    double tt=tsamp(kk), w=(kk==0||kk==NTS-1)?0.5:1.0;
    double u0=eval(u,0.0,tt), u1=eval(u,1.0,tt);
    double cf1=c0_1*(1.0-std::exp(-tt/tau)), cf2=c0_2*(1.0-std::exp(-tt/tau));
    ff1+=w*u0*cf1; ff2+=w*u0*cf2;
    ef1+=w*u1*eval(c1,1.0,tt); ef2+=w*u1*eval(c2,1.0,tt);
  }
  ff1*=dt; ff2*=dt; ef1*=dt; ef2*=dt;
  R.mb1=std::fabs(inv1-(ff1-ef1))/std::max(std::fabs(ff1),1e-30);
  R.mb2=std::fabs(inv2-(ff2-ef2))/std::max(std::fabs(ff2),1e-30);

  std::cout << std::scientific << std::setprecision(3)
            << "  massbal: c1=" << R.mb1 << " c2=" << R.mb2
            << "  umax=" << R.umax << "  near-inlet overshoot=" << R.overshoot
            << " undershoot=" << R.undershoot << "\n";
  return R;
}

int main()
{
  std::cout << "================================================================\n";
#if defined(MC__ERGUN)
  std::cout << "  PSA rung 9c: inlet-sharpness sweep, 7-state INDEX-1 (ERGUN momentum)\n";
#else
  std::cout << "  PSA rung 9c: inlet-sharpness sweep, 7-state INDEX-1 (DARCY momentum)\n";
#endif
  std::cout << "  value-slaved P=Rg T(C_base+c1+c2), C_base=" << C_base
            << "; u = U0 - Kperm dP/dz, U0=" << U0_base << " Kperm=" << Kperm << "\n";
  std::cout << "================================================================\n";

  size_t const NET=6, NNT=8;                     // t-grid: resolves the sharp inlet (near-0 gap ~0.041)
  size_t const NEZ=5, NNZ=6;                      // z-grid: coarse; spatial front is smooth
  std::vector<double> taus = { 0.2, 0.1, 0.05, 0.025 };

  // ======================================================================================
  //  IMPOSITION-MODE COMPARISON  (the reason this driver was instrumented)
  // ======================================================================================
  // Same model, same mesh, same tau -- only IMPOSITION_TYPE changes.  Each mode gets its own
  // cold start (no warm vector shared across modes): a warm start from another mode's
  // solution would let the comparison inherit that mode's interface behaviour, which is the
  // one thing being measured.
  //
  // tau=0.2 is the best-resolved case in the sweep (near-0 gap ~0.041), so all three modes
  // should converge and any spread difference is attributable to the imposition, not to a
  // marginal solve.  Check conv= before reading any number here.
  {
    std::cout << "\n================================================================\n"
              << "  IMPOSITION-MODE COMPARISON at tau=" << taus.front()
              << "  (u is claimed continuous in all three modes)\n"
              << "================================================================\n";
    SR cmp[3];
    for( int m=0; m<3; ++m ){
      std::vector<double> warm_m;                       // cold start, deliberately
      cmp[m] = run_sharp( taus.front(), NET, NNT, NEZ, NNZ, warm_m,
                          "MODE-CMP", IMP_MODES[m], /*want_where=*/true );
    }
    std::cout << "\n  ---- duplicate-node spread by imposition mode ----\n"
              << "  " << std::left << std::setw(12) << "mode"
              << std::setw(7) << "conv"
              << std::right << std::setw(13) << "u(t,z)"
              << std::setw(15) << "max(c1,c2,T,P)"
              << std::setw(13) << "worst_any" << "\n";
    for( int m=0; m<3; ++m )
      std::cout << "  " << std::left << std::setw(12) << cmp[m].imp_name
                << std::setw(7) << ( cmp[m].converged ? "yes" : "NO" )
                << std::right << std::scientific << std::setprecision(3)
                << std::setw(13) << cmp[m].dup_u
                << std::setw(15) << cmp[m].dup_claimed_max
                << std::setw(13) << cmp[m].dup_worst << "\n";
    std::cout << "\n  READ THIS AS:\n"
              << "    u dirty under IC_STRONG only      -> defect is in the Schur/keep-explicit\n"
              << "                                         path; Stage 2 work.\n"
              << "    u dirty in all three modes        -> claim generation for derivative-\n"
              << "                                         defined states; predates this arc.\n"
              << "    u dirty under WEAK+STRONG, clean\n"
              << "      under TRACE                     -> the projection is doing real work and\n"
              << "                                         Stage 2 would inherit the benefit.\n"
              << "    max(c1,c2,T,P) should be ~1e-16 or below in the exact modes in every case;\n"
              << "    if it is not, the finding is broader than u and this table is not the\n"
              << "    right frame for it.\n";
  }

  // ======================================================================================
  //  z-REFINEMENT SWEEP  --  is u's spread truncation-limited or structural?
  // ======================================================================================
  // The hypothesis under test: u = U0 - Kperm dP/dz is the only ALGEBRAIC claimed state, P is
  // C0 but not C1, so dP/dz jumps across a z-element interface and the momentum row at each
  // coincident node solves for u from its own one-sided gradient.  If that is what 1.6e-08
  // is, it is TRUNCATION and must fall as the z-grid is refined.  If it holds roughly
  // constant, it is structural and the interface plan owns it.
  //
  // Refining ne_z at fixed order is the h-sweep; the observable is the CONVERGENCE ORDER, not
  // any single value.  IC_STRONG only -- the mode comparison above already established that
  // TRACE and STRONG agree to four digits, and IC_WEAK is a different (penalty) mechanism
  // whose spread is truncation-limited by construction and would only confuse the fit.
  //
  // Read the c1/c2/T/P column first: those claims are at EXACTLY 0 and must STAY there at
  // every refinement.  If they move, the run is not measuring what this sweep assumes.
  {
    std::cout << "\n================================================================\n"
              << "  z-REFINEMENT SWEEP at tau=" << taus.front() << ", IC_STRONG\n"
              << "================================================================\n";
    std::vector<size_t> const nez_list = { 4, 5, 8, 10, 16 };
    std::vector<double> uu, hh;
    std::cout << "  " << std::left << std::setw(8) << "ne_z"
              << std::setw(8) << "conv"
              << std::right << std::setw(14) << "u(t,z)"
              << std::setw(16) << "max(c1,c2,T,P)"
              << std::setw(12) << "order" << "\n";
    for( size_t nez : nez_list ){
      std::vector<double> warm_r;                      // cold start at every resolution
      SR const R = run_sharp( taus.front(), NET, NNT, nez, NNZ, warm_r,
                              "Z-REFINE", IMP_MODES[2] );
      double order = std::numeric_limits<double>::quiet_NaN();
      if( R.converged && R.dup_u > 0. && !uu.empty() && uu.back() > 0. )
        order = std::log( uu.back()/R.dup_u ) / std::log( double(nez)/hh.back() );
      if( R.converged ){ uu.push_back( R.dup_u ); hh.push_back( double(nez) ); }
      std::cout << "  " << std::left << std::setw(8) << nez
                << std::setw(8) << ( R.converged ? "yes" : "NO" )
                << std::right << std::scientific << std::setprecision(3)
                << std::setw(14) << R.dup_u
                << std::setw(16) << R.dup_claimed_max;
      if( std::isnan(order) ) std::cout << std::setw(12) << "-";
      else                    std::cout << std::setw(12) << std::fixed
                                        << std::setprecision(2) << order;
      std::cout << "\n";
    }
    std::cout << "\n  READ THIS AS:\n"
              << "    order ~ p or better, u falling steadily  -> TRUNCATION.  u's continuity is\n"
              << "        limited by the C1 jump in P, the claim is doing its job, and this is\n"
              << "        NOT an interface-plan defect.  Close the item.\n"
              << "    order ~ 0, u roughly flat                -> STRUCTURAL.  Something is\n"
              << "        capping the enforcement independently of resolution; the interface\n"
              << "        plan owns it and it needs its own investigation.\n"
              << "    max(c1,c2,T,P) leaves 0 at any row       -> STOP.  The other claims were\n"
              << "        exact at ne_z=5; if they are not exact here the comparison is void.\n"
              << "    non-convergence at the finer grids       -> read nothing from those rows\n"
              << "        (conv before accuracy), and note the basin, not the spread.\n";
  }

  std::vector<SR> rows; rows.reserve(taus.size());
  std::vector<double> warm;
  for( size_t i=0;i<taus.size();++i )
    rows.push_back( run_sharp( taus[i], NET, NNT, NEZ, NNZ, warm, i==0? "first(sweep)" : "sweep" ) );

  std::vector<double> warm2;
  SR rec = run_sharp( taus.back(), 10, 8, NEZ, NNZ, warm2, "RECOVERY (finer t-grid)" );

  std::string pre = std::string(PDE38_OUT_PREFIX);
  size_t const nt = rows.size();
  { std::ofstream f( pre+"_feed.out" );
    f << "# t  c2feed for tau ="; for(double tau:taus) f << " " << tau; f << "\n";
    for( size_t kk=0; kk<NTS; ++kk ){ double tt=tsamp(kk); f << std::setprecision(8) << tt;
      for( double tau : taus ) f << " " << c0_2*(1.0-std::exp(-tt/tau)); f << "\n"; } }
  auto colfile = [&]( std::string const& suffix, std::vector<double> SR::* mem ){
    std::ofstream f( pre+suffix );
    f << "# t  per tau ="; for(double tau:taus) f << " " << tau; f << "\n";
    for( size_t kk=0; kk<NTS; ++kk ){ f << std::setprecision(8) << tsamp(kk);
      for( size_t j=0;j<nt;++j ) f << " " << (rows[j].converged? (rows[j].*mem)[kk] : NAN); f << "\n"; } };
  colfile( "_inlet.out",  &SR::c2inlet );
  colfile( "_outlet.out", &SR::c2out );
  colfile( "_vel.out",    &SR::uinlet );
  { std::ofstream f( pre+"_summary.out" );
    f << "# tau  massbal1 massbal2  overshoot |undershoot| umax  iters conv\n";
    for( SR const& R : rows )
      f << std::setprecision(8) << R.tau << " " << R.mb1 << " " << R.mb2 << " "
        << std::max(R.overshoot,0.0) << " " << std::max(-R.undershoot,0.0) << " " << R.umax << " "
        << R.iters << " " << (R.converged?1:0) << "\n";
    f << "# recovery (tau=" << rec.tau << " on " << rec.n_el << "x" << rec.n_nd << "):\n";
    f << std::setprecision(8) << rec.tau << " " << rec.mb1 << " " << rec.mb2 << " "
      << std::max(rec.overshoot,0.0) << " " << std::max(-rec.undershoot,0.0) << " " << rec.umax << " "
      << rec.iters << " " << (rec.converged?1:0) << "\n"; }

  { std::ofstream gp( pre+".gp" );
    gp << "set datafile commentschars '#'\n"
          "set terminal pngcairo size 1000,700 enhanced font 'Helvetica,12'\nset grid\n\n";
    auto overlay = [&]( std::ofstream& g, std::string const& file ){
      for( size_t j=0;j<nt;++j )
        g << (j? ", \\\n     ":"") << "'" << pre << file << "' u 1:" << (2+(int)j)
          << " w l lw 2 t 'tau=" << taus[j] << "'"; g << "\n"; };
    gp << "set output '" << pre << "_feed.png'\n"
          "set title 'Inlet feed programs (-> Heaviside)'\n"
          "set xlabel 't'; set ylabel 'c_{2,feed}'; set xrange [0:0.8]; set key bottom right\n"
          "plot "; overlay( gp, "_feed.out" ); gp << "set xrange [*:*]\n\n";
    gp << "set output '" << pre << "_inlet.png'\n"
          "set title 'Near-inlet c_2(0,t): sharpening + Gibbs ringing as tau->0'\n"
          "set xlabel 't'; set ylabel 'c_2(0,t)'; set xrange [0:1.5]; set key bottom right\n"
          "plot "; overlay( gp, "_inlet.out" ); gp << "set xrange [*:*]\n\n";
    gp << "set output '" << pre << "_outlet.png'\n"
          "set title 'Outlet c_2(1,t): index-1 front-driven breakthrough'\n"
          "set xlabel 't'; set ylabel 'c_2(1,t)'; set key bottom right\n"
          "plot "; overlay( gp, "_outlet.out" ); gp << "\n";
    gp << "set output '" << pre << "_velocity.png'\n"
          "set title 'Darcy inlet velocity u(0,t): front/pressurisation-driven, decays as bed fills'\n"
          "set xlabel 't'; set ylabel 'u(0,t)'; set key top right\n"
          "plot "; overlay( gp, "_vel.out" ); gp << "\n";
    gp << "set output '" << pre << "_summary.png'\n"
          "set title 'Resolution diagnostics vs tau (7-state index-1)'\n"
          "set xlabel 'tau'; set ylabel 'magnitude'; set logscale xy; set key top right\n"
          "plot '" << pre << "_summary.out' every ::0::" << (nt-1) << " u 1:($3>1e-16?$3:1e-16) w lp lw 2 pt 7 t 'massbal comp2', \\\n"
          "     '' every ::0::" << (nt-1) << " u 1:($4>1e-16?$4:1e-16) w lp lw 2 pt 5 t 'near-inlet overshoot', \\\n"
          "     '' every ::0::" << (nt-1) << " u 1:($5>1e-16?$5:1e-16) w lp lw 2 pt 9 t '|undershoot|'\n"
          "unset logscale\n"; }

  std::cout << "\n==================== rung 9c summary (7-state index-1 Heaviside sweep) ====================\n";
  std::cout << std::left << std::setw(10) << "tau" << std::setw(8) << "conv" << std::setw(7) << "iters"
            << std::right << std::setw(12) << "massbal2" << std::setw(12) << "overshoot"
            << std::setw(12) << "undershoot" << std::setw(10) << "umax" << "\n";
  for( SR const& R : rows )
    std::cout << std::left << std::fixed << std::setprecision(4) << std::setw(10) << R.tau
              << std::setw(8) << (R.converged?"yes":"no") << std::setw(7) << R.iters
              << std::right << std::scientific << std::setprecision(2) << std::setw(12) << R.mb2
              << std::setw(12) << std::max(R.overshoot,0.0) << std::setw(12) << R.undershoot
              << std::setw(10) << R.umax << "\n";
  std::cout << "  recovery tau=" << std::fixed << std::setprecision(4) << rec.tau << " on "
            << rec.n_el << "x" << rec.n_nd << ": overshoot " << std::scientific << std::setprecision(2)
            << std::max(rows.back().overshoot,0.0) << " -> " << std::max(rec.overshoot,0.0)
            << " , |undershoot| " << std::max(-rows.back().undershoot,0.0) << " -> " << std::max(-rec.undershoot,0.0) << "\n";

  bool all_conv=true; for( SR const& R : rows ) all_conv &= R.converged;
  std::cout << "\n  wrote " << pre << "_{feed,inlet,outlet,vel,summary}.out  and  " << pre << ".gp\n";
  std::cout << "  visualise:  gnuplot " << pre << ".gp\n";
  std::cout << "  Overall: " << ( all_conv && rec.converged ? "ALL CONVERGED" : "SOME DID NOT CONVERGE" ) << "\n";
  // 2026-09-09: g_nonconv covers the MODE-CMP and Z-REFINE cells that never enter `rows`.
  if( g_nonconv )
    std::cout << "  NOTE: " << g_nonconv << " solve(s) did not converge, including cells outside"
                 " the tau sweep -- these were previously printed but not gated.\n";
  return ( all_conv && rec.converged && g_nonconv == 0 ) ? 0 : 1;
}
