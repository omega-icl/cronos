// OCFE_PDE37_solve2.cpp   ---  PSA rung 9b: INLET-SHARPNESS (Heaviside-limit) sweep
// ===========================================================================
// Variant of rung 9 (PDE36): hold the physical binary bed fixed and sweep the
// inlet ramp time-constant tau DOWNWARD so the feed step  ci,feed = c0_i(1-e^{-t/tau})
// sharpens toward a HEAVISIDE step (tau -> 0).  Two things are being probed:
//
//   (physics)   the packed bed is a LOW-PASS FILTER: as the inlet sharpens, the
//               OUTLET breakthrough curve converges to a tau-independent limiting
//               shape set by D, k and the retardation -- the column smooths the step.
//   (numerics)  the t-collocation must resolve the inlet transient.  CGL nodes
//               cluster at element ends, so the first t-element refines near t=0;
//               once tau drops below the local node spacing the NEAR-INLET region
//               develops Gibbs ringing (overshoot > c0, undershoot < 0), which the
//               bed then filters out before the outlet.
//
// Sweep: tau in {0.2, 0.1, 0.05, 0.025, 0.0125} on an ASYMMETRIC CGL grid -- t fine (8x8, finest
// near-0 gap ~0.031: 0.2/0.1 resolved, 0.05 marginal, 0.025/0.0125 under-resolved) to resolve the
// sharp inlet, z coarse (5x6) since the spatial front is smooth (PDE36 passed at z=5x6).  This keeps
// the Heaviside threshold (a t-only property) while cutting the variable count ~2-3x vs 8x8/8x8.
// tau-continuation warm-starts each solve from the previous (larger) tau.  The RECOVERY run refines
// the t-grid ONLY (12x8 / z 5x6) so it stays tractable.  Heaviside feed-throughput
// limit:  Feed_i -> U c0_i T_end  as tau -> 0.
//
// Output for gnuplot: *_feed.out (feed programs, approach to Heaviside),
//   *_inlet.out (near-inlet response c2(0,t) -> ringing), *_outlet.out
//   (c2(1,t) -> convergence to limiting breakthrough), *_summary.out
//   (massbal / overshoot / undershoot vs tau), and OCFE_PDE37.gp.
// ===========================================================================

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>

#ifndef OCFE_OCFESLV_HEADER
#define OCFE_OCFESLV_HEADER "ocfeslv.hpp"
#endif
#include OCFE_OCFESLV_HEADER

using namespace mc;

#ifndef PDE37_OUT_PREFIX
#define PDE37_OUT_PREFIX "OCFE_PDE37"
#endif

// physical model (identical to rung 9 / PDE36)
static double const U_VEL = 1.0, D_ax = 0.1, F_ph = 0.5, k_ldf = 2.0;
static double const qs1 = 1.0, qs2 = 1.0, b01_L = 3.0, b02_L = 1.0;
static double const beta1 = 2.0, beta2 = 1.0, T0_ref = 1.0;
static double const Cp_e = 1.0, G_cv = 1.0, lam = 0.1, dH1 = 1.0, dH2 = 0.5, hw = 0.5, Tw = 1.0;
static double const c0_1 = 0.5, c0_2 = 0.5, T_end = 5.0;

static size_t const NTS = 401;                 // common history sample grid
static inline double tsamp( size_t kk ){ return double(kk)/double(NTS-1)*T_end; }
static inline double feed_throughput( double c0, double tau )
{ return U_VEL*c0*( T_end - tau*( 1.0 - std::exp( -T_end/tau ) ) ); }
static inline double q1star_T0( double c1g, double c2g )
{ return qs1*b01_L*c1g/( 1.0 + b01_L*c1g + b02_L*c2g ); }
static inline double q2star_T0( double c1g, double c2g )
{ return qs2*b02_L*c2g/( 1.0 + b01_L*c1g + b02_L*c2g ); }

struct SR {
  double tau=0.; size_t n_el=0, n_nd=0; bool converged=false; int iters=0; double final_r=0.;
  double mb1=1., mb2=1., overshoot=0., undershoot=0., Tmax=0., c2bt=0.;
  std::vector<double> c2out, c2inlet;          // c2(1,t_k), c2(0,t_k) on the common grid
};

static SR run_sharp( double tau, size_t ne_t, size_t nn_t, size_t ne_z, size_t nn_z,
                     std::vector<double>& warm, std::string const& label )
{
  SR R; R.tau=tau; R.n_el=ne_t; R.n_nd=nn_t;   // n_el/n_nd hold the t-grid (the sharpness-relevant one)
  std::cout << "\n---- tau=" << std::fixed << std::setprecision(5) << tau
            << "  t-grid " << ne_t << "x" << nn_t << " z-grid " << ne_z << "x" << nn_z
            << "  " << label << " ----\n";

  FFGraph DAG;
  FFVar t  = DAG.add_var( "t" );
  FFVar z  = DAG.add_var( "z" );
  FFVar c1 = DAG.add_var( "c1(t,z)" );
  FFVar c2 = DAG.add_var( "c2(t,z)" );
  FFVar q1 = DAG.add_var( "q1(t,z)" );
  FFVar q2 = DAG.add_var( "q2(t,z)" );
  FFVar T  = DAG.add_var( "T(t,z)" );

  FFPartial  OpP;
  FFIntegral OpI;

  FFVar c1feed = c0_1*( 1.0 - exp( -t/tau ) );
  FFVar c2feed = c0_2*( 1.0 - exp( -t/tau ) );
  FFVar b1 = b01_L*exp( beta1*( 1.0/T - 1.0/T0_ref ) );
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
  FFVar IC_c1 = c1, IC_c2 = c2, IC_q1 = q1, IC_q2 = q2, IC_T = T - T0_ref;
  FFVar BC_L1 = U_VEL*c1 - D_ax*OpP( c1, z ) - U_VEL*c1feed;
  FFVar BC_L2 = U_VEL*c2 - D_ax*OpP( c2, z ) - U_VEL*c2feed;
  FFVar BC_LT = G_cv*T - lam*OpP( T, z ) - G_cv*T0_ref;
  FFVar BC_U1 = OpP( c1, z ), BC_U2 = OpP( c2, z ), BC_UT = OpP( T, z );

  FFVar Inv1 = OpI( c1 + F_ph*q1, z ), Inv2 = OpI( c2 + F_ph*q2, z );
  FFVar Eff1 = OpI( U_VEL*c1, t ),     Eff2 = OpI( U_VEL*c2, t );

  OCFESLV oc( &DAG );
  oc.add_domain( t, FFDom( 0., T_end, ne_t, FFDom::CGL, nn_t ) );
  oc.add_domain( z, FFDom( 0., 1.0,   ne_z, FFDom::CGL, nn_z ) );
  oc.add_state ( c1, {t,z} ); oc.add_state( c2, {t,z} );
  oc.add_state ( q1, {t,z} ); oc.add_state( q2, {t,z} ); oc.add_state( T, {t,z} );
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
  oc.add_output( Eff1, {z}, {1.0} );
  oc.add_output( Eff2, {z}, {1.0} );

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

  if( !oc.setup() ){ std::cerr << "ERROR: setup failed\n"; return R; }
  R.n_el=ne_t; R.n_nd=nn_t;
  size_t const nVar=oc.n_colloc_sta(), nEqn=oc.n_colloc_eqn(), nFct=oc.n_colloc_fct();
  if( nVar!=nEqn ){ std::cerr << "ERROR: not square\n"; return R; }

  std::vector<double> varInit, inpInit;
  if( !oc.init( varInit, inpInit, nullptr ) ){ std::cerr << "ERROR: init failed\n"; return R; }
  std::vector<double> xv = ( warm.size()==varInit.size() ? warm : varInit );   // tau-continuation
  OCFESLV::SolveReport rep = oc.solve( xv.data(), nullptr, nullptr );
  R.converged=rep.converged; R.iters=rep.iterations; R.final_r=rep.final_residual;
  std::cout << "  nVar=" << nVar << " converged=" << (rep.converged?"yes":"no")
            << " iters=" << rep.iterations << " final|r|=" << std::scientific << std::setprecision(3)
            << rep.final_residual << "\n";
  if( !rep.converged ) return R;
  if( warm.size()==varInit.size() || warm.empty() ) warm = xv;                 // save for next (same grid)

  auto eval = [&]( FFVar const& V, double zz, double tt ){
    OCFESLV::t_Coord pt; pt[z]=zz; pt[t]=tt;
    return oc.eval_colloc<double>( V, pt, xv.data(), nullptr, nullptr ); };

  // histories on the common grid (outlet z=1, near-inlet z=0)
  R.c2out.resize(NTS); R.c2inlet.resize(NTS);
  for( size_t kk=0; kk<NTS; ++kk ){
    double tt=tsamp(kk);
    R.c2out[kk]   = eval(c2,1.0,tt);
    R.c2inlet[kk] = eval(c2,0.0,tt);
  }
  R.c2bt=0.; for( size_t kk=0;kk<NTS;++kk ){ double tt=tsamp(kk);
    if( R.c2bt==0. && R.c2out[kk]>0.5*c0_2 && kk>0 ){ R.c2bt=tt; break; } }

  // overshoot/undershoot: FINE scan of the near-inlet region t in [0,1], z in [0,1]
  double cmax=-1e30, cmin=1e30;
  int const NZ=60, NTf=600;
  for( int iz=0; iz<=NZ; ++iz ){ double zz=double(iz)/NZ;
    for( int it=0; it<=NTf; ++it ){ double tt=double(it)/NTf*1.0;
      double v1=eval(c1,zz,tt), v2=eval(c2,zz,tt);
      cmax=std::max(cmax,std::max(v1,v2)); cmin=std::min(cmin,std::min(v1,v2)); } }
  R.overshoot = cmax - std::max(c0_1,c0_2);     // >0 means above feed (Gibbs overshoot)
  R.undershoot= cmin;                            // <0 means below zero (Gibbs undershoot)
  // Tmax over full domain
  R.Tmax=-1e30; for( int iz=0;iz<=40;++iz ) for( int it=0;it<=40;++it )
    R.Tmax=std::max(R.Tmax, eval(T,double(iz)/40.0,double(it)/40.0*T_end));

  // species mass balance
  std::vector<double> res(nEqn,0.), fct(nFct,0.);
  if( !oc.eval( res.data(), fct.data(), xv.data(), nullptr, nullptr ) ){ std::cerr << "ERROR: output eval\n"; return R; }
  double inv1=fct[oc.row_fct(0)-nEqn], inv2=fct[oc.row_fct(1)-nEqn];
  double eff1=fct[oc.row_fct(2)-nEqn], eff2=fct[oc.row_fct(3)-nEqn];
  double f1=feed_throughput(c0_1,tau), f2=feed_throughput(c0_2,tau);
  R.mb1=std::fabs(inv1-(f1-eff1))/std::max(f1,1e-30);
  R.mb2=std::fabs(inv2-(f2-eff2))/std::max(f2,1e-30);

  std::cout << std::scientific << std::setprecision(3)
            << "  massbal: c1=" << R.mb1 << " c2=" << R.mb2
            << "  near-inlet overshoot=" << R.overshoot << " undershoot=" << R.undershoot
            << "  c2 bt50=" << std::fixed << std::setprecision(3) << R.c2bt << "\n";
  return R;
}

int main()
{
  std::cout << "================================================================\n";
  std::cout << "  PSA rung 9b: inlet-sharpness sweep tau -> 0 (Heaviside limit)\n";
  std::cout << "  binary bed (PDE36 model); feed ci=c0_i(1-e^{-t/tau}); sweep tau down\n";
  std::cout << "================================================================\n";

  size_t const NET=8, NNT=8;                     // t-grid: resolves the sharp inlet (near-0 gap ~0.031)
  size_t const NEZ=5, NNZ=6;                      // z-grid: coarse; spatial front is smooth (PDE36 passed at 5x6)
  std::vector<double> taus = { 0.2, 0.1, 0.05, 0.025, 0.0125 };

  std::vector<SR> rows; rows.reserve(taus.size());
  std::vector<double> warm;                       // tau-continuation across the (fixed-grid) sweep
  for( double tau : taus ) rows.push_back( run_sharp( tau, NET, NNT, NEZ, NNZ, warm, "sweep" ) );

  // resolution-recovery: sharpest tau on a finer grid (fresh start)
  std::vector<double> warm2;
  // resolution-recovery: refine the t-grid ONLY at the sharpest tau (z stays coarse) -> tractable
  SR rec = run_sharp( taus.back(), 12, 8, NEZ, NNZ, warm2, "RECOVERY (finer t-grid)" );

  // ---- write gnuplot data ----
  std::string pre = std::string(PDE37_OUT_PREFIX);
  size_t const nt = rows.size();
  { std::ofstream f( pre+"_feed.out" );
    f << "# t  c2feed for tau ="; for(double tau:taus) f << " " << tau; f << "\n";
    for( size_t kk=0; kk<NTS; ++kk ){ double tt=tsamp(kk); f << std::setprecision(8) << tt;
      for( double tau : taus ) f << " " << c0_2*(1.0-std::exp(-tt/tau)); f << "\n"; } }
  { std::ofstream f( pre+"_inlet.out" );
    f << "# t  c2(0,t) for tau ="; for(double tau:taus) f << " " << tau; f << "\n";
    for( size_t kk=0; kk<NTS; ++kk ){ f << std::setprecision(8) << tsamp(kk);
      for( size_t j=0;j<nt;++j ) f << " " << (rows[j].converged? rows[j].c2inlet[kk] : NAN); f << "\n"; } }
  { std::ofstream f( pre+"_outlet.out" );
    f << "# t  c2(1,t) for tau ="; for(double tau:taus) f << " " << tau; f << "\n";
    for( size_t kk=0; kk<NTS; ++kk ){ f << std::setprecision(8) << tsamp(kk);
      for( size_t j=0;j<nt;++j ) f << " " << (rows[j].converged? rows[j].c2out[kk] : NAN); f << "\n"; } }
  { std::ofstream f( pre+"_summary.out" );
    f << "# tau  massbal1  massbal2  overshoot  |undershoot|  iters  conv  (last row = recovery)\n";
    for( SR const& R : rows )
      f << std::setprecision(8) << R.tau << " " << R.mb1 << " " << R.mb2 << " "
        << std::max(R.overshoot,0.0) << " " << std::max(-R.undershoot,0.0) << " "
        << R.iters << " " << (R.converged?1:0) << "\n";
    f << "# recovery (tau=" << rec.tau << " on " << rec.n_el << "x" << rec.n_nd << "):\n";
    f << std::setprecision(8) << rec.tau << " " << rec.mb1 << " " << rec.mb2 << " "
      << std::max(rec.overshoot,0.0) << " " << std::max(-rec.undershoot,0.0) << " "
      << rec.iters << " " << (rec.converged?1:0) << "\n"; }

  // ---- emit gnuplot script ----
  { std::ofstream gp( pre+".gp" );
    gp << "set datafile commentschars '#'\n"
          "set terminal pngcairo size 1000,700 enhanced font 'Helvetica,12'\n"
          "set grid\n\n";
    auto overlay = [&]( std::ofstream& g, std::string const& file, int yidx_base ){
      for( size_t j=0;j<nt;++j ){
        g << (j? ", \\\n     ":"") << "'" << pre << file << "' u 1:" << (2+(int)j)
          << " w l lw 2 t 'tau=" << taus[j] << "'"; (void)yidx_base; }
      g << "\n"; };
    // 1. feed programs -> approach to Heaviside
    gp << "set output '" << pre << "_feed.png'\n"
          "set title 'Inlet feed programs c_{2,feed}(t) = c0 (1 - e^{-t/tau})  (-> Heaviside)'\n"
          "set xlabel 't'; set ylabel 'c_{2,feed}'; set xrange [0:0.8]; set key bottom right\n"
          "plot "; overlay( gp, "_feed.out", 0 ); gp << "set xrange [*:*]\n\n";
    // 2. near-inlet response -> ringing
    gp << "set output '" << pre << "_inlet.png'\n"
          "set title 'Near-inlet response c_2(z=0,t): sharpening + Gibbs ringing as tau->0'\n"
          "set xlabel 't'; set ylabel 'c_2(0,t)'; set xrange [0:1.5]; set key bottom right\n"
          "plot "; overlay( gp, "_inlet.out", 0 ); gp << "set xrange [*:*]\n\n";
    // 3. outlet -> convergence to limiting breakthrough
    gp << "set output '" << pre << "_outlet.png'\n"
          "set title 'Outlet breakthrough c_2(z=1,t): converges to a tau-independent limit'\n"
          "set xlabel 't'; set ylabel 'c_2(1,t)'; set key bottom right\n"
          "plot "; overlay( gp, "_outlet.out", 0 ); gp << "\n";
    // 4. summary vs tau (log-log)
    gp << "set output '" << pre << "_summary.png'\n"
          "set title 'Resolution diagnostics vs tau (smaller tau = sharper step)'\n"
          "set xlabel 'tau'; set ylabel 'magnitude'; set logscale xy; set key top right\n"
          "plot '" << pre << "_summary.out' every ::0::" << (nt-1) << " u 1:($2>1e-16?$2:1e-16) w lp lw 2 pt 7 t 'massbal comp1', \\\n"
          "     '' every ::0::" << (nt-1) << " u 1:($4>1e-16?$4:1e-16) w lp lw 2 pt 5 t 'near-inlet overshoot', \\\n"
          "     '' every ::0::" << (nt-1) << " u 1:($5>1e-16?$5:1e-16) w lp lw 2 pt 9 t '|undershoot|'\n"
          "unset logscale\n"; }

  // ---- console summary table ----
  std::cout << "\n==================== rung 9b summary (Heaviside sweep) ====================\n";
  std::cout << std::left << std::setw(10) << "tau" << std::setw(8) << "conv"
            << std::setw(7) << "iters" << std::right << std::setw(13) << "massbal2"
            << std::setw(13) << "overshoot" << std::setw(13) << "undershoot"
            << std::setw(9) << "bt50" << "\n";
  for( SR const& R : rows )
    std::cout << std::left << std::fixed << std::setprecision(4) << std::setw(10) << R.tau
              << std::setw(8) << (R.converged?"yes":"no") << std::setw(7) << R.iters
              << std::right << std::scientific << std::setprecision(3) << std::setw(13) << R.mb2
              << std::setw(13) << std::max(R.overshoot,0.0) << std::setw(13) << R.undershoot
              << std::fixed << std::setprecision(3) << std::setw(9) << R.c2bt << "\n";
  std::cout << "  recovery tau=" << std::fixed << std::setprecision(4) << rec.tau
            << " on " << rec.n_el << "x" << rec.n_nd << ": overshoot "
            << std::scientific << std::setprecision(3) << std::max(rows.back().overshoot,0.0)
            << " -> " << std::max(rec.overshoot,0.0)
            << " , |undershoot| " << std::max(-rows.back().undershoot,0.0)
            << " -> " << std::max(-rec.undershoot,0.0)
            << "  (refinement " << ( (std::max(rec.overshoot,0.0)<std::max(rows.back().overshoot,0.0)*0.7
                                      || std::max(-rec.undershoot,0.0)<std::max(-rows.back().undershoot,0.0)*0.7)
                                     ? "reduces ringing -> resolution effect" : "inconclusive" ) << ")\n";

  bool all_conv=true; for( SR const& R : rows ) all_conv &= R.converged;
  std::cout << "\n  wrote " << pre << "_{feed,inlet,outlet,summary}.out  and  " << pre << ".gp\n";
  std::cout << "  visualise:  gnuplot " << pre << ".gp\n";
  std::cout << "  Overall: " << ( all_conv && rec.converged ? "ALL CONVERGED (sweep complete)" : "SOME DID NOT CONVERGE" ) << "\n";
  return ( all_conv && rec.converged ) ? 0 : 1;
}
