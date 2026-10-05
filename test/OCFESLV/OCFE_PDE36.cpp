// OCFE_PDE36_solve2.cpp   ---  PSA case study, RUNG 9 (PHYSICAL binary breakthrough sweep)
// ===========================================================================
// First PHYSICAL multi-component run (no MMS).  A regenerated bed is fed a binary
// mixture; the strongly-adsorbed component (1) is retarded and breaks through
// LATE, while the weakly-adsorbed component (2) breaks through EARLY and ROLLS UP
// (its outlet concentration overshoots its feed value as component 1 displaces it
// from the adsorbed phase).  The adsorption exotherm drives a thermal wave; wall
// cooling (hw) relaxes it.  Competitive van't Hoff Langmuir couples everything.
//
//   mass_i: dci/dt + d(U ci)/dz - D d2ci/dz2 + F dqi/dt              = 0   (i=1,2)
//   LDF_i:  dqi/dt - k ( qi*(c1,c2,T) - qi )                         = 0
//   energy: Cp dT/dt + G dT/dz - lam d2T/dz2 - F(dH1 dq1/dt+dH2 dq2/dt)
//                                                       + hw (T-Tw)  = 0
//   b_i(T) = b0_i exp( beta_i ( 1/T - 1/T0 ) );  qi* = qs_i b_i ci /(1+b1 c1+b2 c2)
//   IC: c1=c2=q1=q2=0, T=T0 (regenerated, clean bed).
//   inlet: Danckwerts mass  ci,feed = c0_i (1 - e^{-t/tau}); thermal feed at T0.
//   outlet: dci/dz = dT/dz = 0.   Constant interstitial velocity U.
//
// Validation (no exact solution -- physical + numerical diagnostics):
//   * per-species mass balance  int(ci+F qi)|_T dz = Feed_i - Effluent_i  (closes);
//   * outlet curves c1(1,t), c2(1,t), T(1,t) self-converge under p-refinement;
//   * physical sanity: c2 roll-up (outlet overshoot), c1 breaks through after c2,
//     thermal wave rises above T0 and stays bounded.
// Output for gnuplot:  *_field.out (z,t fields), *_outlet.out (breakthrough
//   histories + purity), *_snapshots.out (spatial profiles), and OCFE_PDE36.gp
//   (ready-to-run script producing PNGs).
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

#ifndef PDE36_OUT_PREFIX
#define PDE36_OUT_PREFIX "OCFE_PDE36"
#endif

// transport / kinetics
static double const U_VEL = 1.0, D_ax = 0.1, F_ph = 0.5, k_ldf = 2.0;
// competitive van't Hoff Langmuir: component 1 strong (b01=3), component 2 weak (b02=1)
static double const qs1 = 1.0, qs2 = 1.0, b01_L = 3.0, b02_L = 1.0;
static double const beta1 = 2.0, beta2 = 1.0, T0_ref = 1.0;
// energy
static double const Cp_e = 1.0, G_cv = 1.0, lam = 0.1, dH1 = 1.0, dH2 = 0.5, hw = 0.5, Tw = 1.0;
// feed program / horizon
static double const c0_1 = 0.5, c0_2 = 0.5, tau_in = 0.15, T_end = 5.0;

static inline double c1feed_d( double t ){ return c0_1*( 1.0 - std::exp( -t/tau_in ) ); }
static inline double c2feed_d( double t ){ return c0_2*( 1.0 - std::exp( -t/tau_in ) ); }
static inline double feed_throughput( double c0 )
{ return U_VEL*c0*( T_end - tau_in*( 1.0 - std::exp( -T_end/tau_in ) ) ); }
// competitive Langmuir at T0 (b(T0)=b0) -- used only for the Newton initial guess
static inline double q1star_T0( double c1g, double c2g )
{ return qs1*b01_L*c1g/( 1.0 + b01_L*c1g + b02_L*c2g ); }
static inline double q2star_T0( double c1g, double c2g )
{ return qs2*b02_L*c2g/( 1.0 + b01_L*c1g + b02_L*c2g ); }

struct BT {
  std::string tag;
  bool   converged=false, ok=false;
  size_t nVar=0; int iters=0;
  double c1max_out=0., c2max_out=0., Tmax=0., Tmin=0.;
  double t1_bt=0., t2_bt=0.;                 // 50%-feed breakthrough times (z=1)
  double mb1=1., mb2=1.;                      // relative species mass-balance error
  std::vector<double> c1o, c2o, To;           // outlet histories at fixed sample times
};

static std::vector<double> sample_times()
{ std::vector<double> ts; for( int k=0;k<=40;++k ) ts.push_back( double(k)*T_end/40.0 ); return ts; }

static BT run_BT( FFDom::TYPE coltype, std::string const& cname,
                  size_t n_el, size_t n_nd, bool write_files )
{
  BT R;
  R.tag = cname + " n_el=" + std::to_string(n_el) + " n_nd=" + std::to_string(n_nd);
  std::cout << "\n---- binary breakthrough  " << R.tag << " ----\n";

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

  FFVar c1feed = c0_1*( 1.0 - exp( -t/tau_in ) );
  FFVar c2feed = c0_2*( 1.0 - exp( -t/tau_in ) );
  FFVar b1 = b01_L*exp( beta1*( 1.0/T - 1.0/T0_ref ) );
  FFVar b2 = b02_L*exp( beta2*( 1.0/T - 1.0/T0_ref ) );
  FFVar den = 1.0 + b1*c1 + b2*c2;
  FFVar q1star = qs1*b1*c1/den;
  FFVar q2star = qs2*b2*c2/den;

  // mass balances (conservative flux form; constant U so d(U ci)/dz = U dci/dz, mass-conservative)
  FFVar CONT1 = OpP( c1, t ) + U_VEL*OpP( c1, z ) - D_ax*OpP( OpP( c1, z ), z ) + F_ph*OpP( q1, t );
  FFVar CONT2 = OpP( c2, t ) + U_VEL*OpP( c2, z ) - D_ax*OpP( OpP( c2, z ), z ) + F_ph*OpP( q2, t );
  FFVar LDF1  = OpP( q1, t ) - k_ldf*( q1star - q1 );
  FFVar LDF2  = OpP( q2, t ) - k_ldf*( q2star - q2 );
  FFVar ENE_T = Cp_e*OpP( T, t ) + G_cv*OpP( T, z ) - lam*OpP( OpP( T, z ), z )
              - F_ph*( dH1*OpP( q1, t ) + dH2*OpP( q2, t ) ) + hw*( T - Tw );
  FFVar IC_c1 = c1;
  FFVar IC_c2 = c2;
  FFVar IC_q1 = q1;
  FFVar IC_q2 = q2;
  FFVar IC_T  = T - T0_ref;
  FFVar BC_L1 = U_VEL*c1 - D_ax*OpP( c1, z ) - U_VEL*c1feed;   // Danckwerts inlet
  FFVar BC_L2 = U_VEL*c2 - D_ax*OpP( c2, z ) - U_VEL*c2feed;
  FFVar BC_LT = G_cv*T - lam*OpP( T, z ) - G_cv*T0_ref;        // thermal feed at T0
  FFVar BC_U1 = OpP( c1, z );                                  // zero-gradient outlets
  FFVar BC_U2 = OpP( c2, z );
  FFVar BC_UT = OpP( T,  z );

  // species inventories at final time and cumulative effluents (for mass balance)
  FFVar Inv1 = OpI( c1 + F_ph*q1, z );
  FFVar Inv2 = OpI( c2 + F_ph*q2, z );
  FFVar Eff1 = OpI( U_VEL*c1, t );
  FFVar Eff2 = OpI( U_VEL*c2, t );

  OCFESLV oc( &DAG );
  oc.add_domain( t, FFDom( 0., T_end, n_el, coltype, n_nd ) );
  oc.add_domain( z, FFDom( 0., 1.0,   n_el, coltype, n_nd ) );
  oc.add_state ( c1, {t,z} );
  oc.add_state ( c2, {t,z} );
  oc.add_state ( q1, {t,z} );
  oc.add_state ( q2, {t,z} );
  oc.add_state ( T,  {t,z} );
  // physically-motivated smooth Newton initial guess (decaying feed into the bed)
  oc.update_ref( c1, [&]( OCFESLV::t_Coord const& cr ){
    double zz=cr.at(z), tt=cr.at(t); return c1feed_d(tt)*(1.0-0.5*zz); } );
  oc.update_ref( c2, [&]( OCFESLV::t_Coord const& cr ){
    double zz=cr.at(z), tt=cr.at(t); return c2feed_d(tt)*(1.0-0.5*zz); } );
  oc.update_ref( q1, [&]( OCFESLV::t_Coord const& cr ){
    double zz=cr.at(z), tt=cr.at(t);
    return q1star_T0( c1feed_d(tt)*(1.0-0.5*zz), c2feed_d(tt)*(1.0-0.5*zz) ); } );
  oc.update_ref( q2, [&]( OCFESLV::t_Coord const& cr ){
    double zz=cr.at(z), tt=cr.at(t);
    return q2star_T0( c1feed_d(tt)*(1.0-0.5*zz), c2feed_d(tt)*(1.0-0.5*zz) ); } );
  oc.update_ref( T,  [&]( OCFESLV::t_Coord const& cr ){
    double zz=cr.at(z), tt=cr.at(t);
    double c1g=c1feed_d(tt)*(1.0-0.5*zz), c2g=c2feed_d(tt)*(1.0-0.5*zz);
    return T0_ref + F_ph*( dH1*q1star_T0(c1g,c2g) + dH2*q2star_T0(c1g,c2g) )/Cp_e; } );
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
  oc.add_output( Eff1, {z}, {1.0} );     // fct 2
  oc.add_output( Eff2, {z}, {1.0} );     // fct 3

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
  R.nVar = oc.n_colloc_sta();
  size_t const nEqn = oc.n_colloc_eqn(), nFct = oc.n_colloc_fct();
  if( R.nVar != nEqn ){ std::cerr << "ERROR: not square (" << R.nVar << "/" << nEqn << ")\n"; return R; }

  std::vector<double> varInit, inpInit;
  if( !oc.init( varInit, inpInit, nullptr ) ){ std::cerr << "ERROR: init failed\n"; return R; }
  std::vector<double> xv = varInit;
  OCFESLV::SolveReport rep = oc.solve( xv.data(), nullptr, nullptr );
  R.converged = rep.converged; R.iters = rep.iterations;
  std::cout << "  nVar=" << R.nVar << "  converged=" << (rep.converged?"yes":"no")
            << " iters=" << rep.iterations
            << " final|r|=" << std::scientific << std::setprecision(3) << rep.final_residual << "\n";
  if( !rep.converged ) return R;

  // outlet histories at fixed sample times (z=1)
  auto eval = [&]( FFVar const& V, double zz, double tt ){
    OCFESLV::t_Coord pt; pt[z]=zz; pt[t]=tt;
    return oc.eval_colloc<double>( V, pt, xv.data(), nullptr, nullptr ); };
  auto ts = sample_times();
  R.c1o.reserve(ts.size()); R.c2o.reserve(ts.size()); R.To.reserve(ts.size());
  for( double tt : ts ){
    R.c1o.push_back( eval(c1,1.0,tt) );
    R.c2o.push_back( eval(c2,1.0,tt) );
    R.To .push_back( eval(T ,1.0,tt) );
  }

  // field ranges + breakthrough times (first crossing of 50% feed at z=1) + roll-up
  R.Tmin=1e30; R.Tmax=-1e30;
  for( int iz=0;iz<=40;++iz ) for( int it=0;it<=40;++it ){
    double Tv = eval( T, double(iz)/40.0, double(it)/40.0*T_end );
    R.Tmin=std::min(R.Tmin,Tv); R.Tmax=std::max(R.Tmax,Tv);
  }
  int const NT=400;
  double prev1=0., prev2=0.;
  for( int kk=0; kk<=NT; ++kk ){
    double tt=double(kk)/NT*T_end;
    double v1=eval(c1,1.0,tt), v2=eval(c2,1.0,tt);
    R.c1max_out=std::max(R.c1max_out,v1); R.c2max_out=std::max(R.c2max_out,v2);
    if( R.t1_bt==0. && v1>0.5*c0_1 && kk>0 ) R.t1_bt=tt;
    if( R.t2_bt==0. && v2>0.5*c0_2 && kk>0 ) R.t2_bt=tt;
    prev1=v1; prev2=v2;
  }
  (void)prev1; (void)prev2;

  // species mass balance:  Inv_i(T_end) ?= Feed_i - Eff_i
  auto const& fct = oc.val_functions();   // window-summed evolution integrals; correct for both modes
  if( fct.size() < nFct ){ std::cerr << "ERROR: output functionals unavailable\n"; return R; }
  double inv1=fct[oc.row_fct(0)-nEqn], inv2=fct[oc.row_fct(1)-nEqn];
  double eff1=fct[oc.row_fct(2)-nEqn], eff2=fct[oc.row_fct(3)-nEqn];
  double feed1=feed_throughput(c0_1), feed2=feed_throughput(c0_2);
  R.mb1=std::fabs( inv1-(feed1-eff1) )/std::max(feed1,1e-30);
  R.mb2=std::fabs( inv2-(feed2-eff2) )/std::max(feed2,1e-30);

  std::cout << std::scientific << std::setprecision(4)
            << "  c1max_out=" << R.c1max_out << " c2max_out=" << R.c2max_out
            << " (c0=" << c0_2 << ")  rollup=" << (R.c2max_out>c0_2*1.02?"yes":"no") << "\n"
            << "  bt50: t1(strong)=" << R.t1_bt << " t2(weak)=" << R.t2_bt
            << "  ordered=" << (R.t1_bt>R.t2_bt?"yes":"no") << "\n"
            << "  massbal_rel: comp1=" << R.mb1 << " comp2=" << R.mb2
            << "  Tmax=" << R.Tmax << " (dT=" << (R.Tmax-T0_ref) << ")\n";

  bool csane = (R.c1max_out>-0.02)&&(R.c2max_out>-0.02)&&(R.c1max_out<c0_1+0.6)&&(R.c2max_out<c0_2+0.6);
  bool Tsane = (R.Tmax>T0_ref+0.02)&&(R.Tmax<T0_ref+3.0)&&(R.Tmin>T0_ref-0.05);
  R.ok = R.converged && (R.mb1<1e-2)&&(R.mb2<1e-2) && csane && Tsane;

  if( write_files ){
    std::string pre = std::string(PDE36_OUT_PREFIX)+"_"+cname;
    // 1) full field (gnuplot pm3d: blank line between z-blocks)
    { std::ofstream of( pre+"_field.out" );
      of << "# z  t  c1  c2  q1  q2  T   (pm3d: splot u 1:2:N)\n";
      int const NG=61;
      for( int iz=0; iz<NG; ++iz ){ double zz=double(iz)/(NG-1);
        for( int it=0; it<NG; ++it ){ double tt=double(it)/(NG-1)*T_end;
          of << std::setprecision(8) << zz << " " << tt << " "
             << eval(c1,zz,tt) << " " << eval(c2,zz,tt) << " "
             << eval(q1,zz,tt) << " " << eval(q2,zz,tt) << " " << eval(T,zz,tt) << "\n"; }
        of << "\n"; } }
    // 2) outlet breakthrough histories + purity (z=1)
    { std::ofstream oo( pre+"_outlet.out" );
      oo << "# t  c1(1,t)  c2(1,t)  T(1,t)  q1(1,t)  q2(1,t)  y1=c1/(c1+c2)  c1feed  c2feed\n";
      for( int kk=0;kk<=400;++kk ){ double tt=double(kk)/400.0*T_end;
        double v1=eval(c1,1.0,tt), v2=eval(c2,1.0,tt), den2=v1+v2;
        oo << std::setprecision(8) << tt << " " << v1 << " " << v2 << " "
           << eval(T,1.0,tt) << " " << eval(q1,1.0,tt) << " " << eval(q2,1.0,tt) << " "
           << ( den2>1e-12 ? v1/den2 : 0.0 ) << " " << c1feed_d(tt) << " " << c2feed_d(tt) << "\n"; } }
    // 3) spatial snapshots at several times  (DOUBLE blank line per block: gnuplot `index` needs it)
    { std::ofstream os( pre+"_snapshots.out" );
      os << "# z  c1  c2  q1  q2  T   (index blocks 0..4 at t = 1,2,3,4,5)\n";
      for( double tt : { 1.0,2.0,3.0,4.0,5.0 } ){ os << "# t=" << tt << "\n";
        for( int kk=0;kk<=120;++kk ){ double zz=double(kk)/120.0;
          os << std::setprecision(8) << zz << " "
             << eval(c1,zz,tt) << " " << eval(c2,zz,tt) << " "
             << eval(q1,zz,tt) << " " << eval(q2,zz,tt) << " " << eval(T,zz,tt) << "\n"; }
        os << "\n\n"; } }
    // 4) ready-to-run gnuplot script (produces PNGs from the three data files)
    { std::ofstream gp( std::string(PDE36_OUT_PREFIX)+".gp" );
      gp << "# gnuplot " << PDE36_OUT_PREFIX << ".gp   ->  PNGs of the binary breakthrough\n"
            "# data written by OCFE_PDE36_solve2 (tag " << cname << ")\n"
            "set datafile commentschars '#'\n"
            "set terminal pngcairo size 1000,700 enhanced font 'Helvetica,12'\n"
            "pre = '" << pre << "'\n\n"
            "# --- 1. outlet breakthrough curves (mass + thermal) ---\n"
            "set output '" << PDE36_OUT_PREFIX << "_breakthrough.png'\n"
            "set title 'Binary breakthrough at column outlet (z=1)'\n"
            "set xlabel 't'; set ylabel 'concentration c_i(1,t)'\n"
            "set y2label 'temperature T(1,t)'; set ytics nomirror; set y2tics\n"
            "set grid; set key top left\n"
            "plot pre.'_outlet.out' u 1:2 w l lw 2 lc rgb 'blue'   t 'c_1 (strong)', \\\n"
            "     ''             u 1:3 w l lw 2 lc rgb 'forest-green' t 'c_2 (weak, roll-up)', \\\n"
            "     ''             u 1:8 w l dt 2 lc rgb 'blue'   t 'c_1,feed', \\\n"
            "     ''             u 1:9 w l dt 2 lc rgb 'forest-green' t 'c_2,feed', \\\n"
            "     ''             u 1:4 axes x1y2 w l lw 2 lc rgb 'red' t 'T(1,t)'\n\n"
            "# --- 2. outlet purity of component 1 ---\n"
            "set output '" << PDE36_OUT_PREFIX << "_purity.png'\n"
            "set title 'Outlet mole fraction of component 1'\n"
            "unset y2label; unset y2tics; set ytics mirror; set ylabel 'y_1 = c_1/(c_1+c_2)'\n"
            "set yrange [0:1]\n"
            "plot pre.'_outlet.out' u 1:7 w l lw 2 lc rgb 'purple' t 'y_1(1,t)'\n"
            "unset yrange\n\n"
            "# --- 3. concentration / temperature fields (pm3d) ---\n"
            "set pm3d map; set xlabel 'z'; set ylabel 't'; unset key\n"
            "set output '" << PDE36_OUT_PREFIX << "_field_c1.png'\n"
            "set title 'c_1(z,t)  (strong, retarded)'; splot pre.'_field.out' u 1:2:3 notitle\n"
            "set output '" << PDE36_OUT_PREFIX << "_field_c2.png'\n"
            "set title 'c_2(z,t)  (weak, early + roll-up)'; splot pre.'_field.out' u 1:2:4 notitle\n"
            "set output '" << PDE36_OUT_PREFIX << "_field_T.png'\n"
            "set title 'T(z,t)  (thermal wave)'; splot pre.'_field.out' u 1:2:7 notitle\n"
            "unset pm3d\n\n"
            "# --- 4. spatial snapshots (profiles vs z at t=1..5) ---\n"
            "set output '" << PDE36_OUT_PREFIX << "_snapshots.png'\n"
            "set title 'Spatial profiles vs z (t = 1,2,3,4,5)'\n"
            "set xlabel 'z'; set ylabel 'concentration'; set key top right; set grid\n"
            "plot for [i=0:4] pre.'_snapshots.out' index i u 1:2 w l lw 2 t sprintf('c_1 t=%d',i+1), \\\n"
            "     for [i=0:4] pre.'_snapshots.out' index i u 1:3 w l dt 2 t sprintf('c_2 t=%d',i+1)\n"; }
    std::cout << "  wrote " << pre << "_{field,outlet,snapshots}.out  and  "
              << PDE36_OUT_PREFIX << ".gp\n";
    std::cout << "  visualise:  gnuplot " << PDE36_OUT_PREFIX << ".gp   (produces "
              << PDE36_OUT_PREFIX << "_*.png)\n";
  }

  std::cout << "  " << R.tag << ": " << (R.ok?"PASS":"FAIL") << "\n";
  return R;
}

static double curve_diff( std::vector<double> const& a, std::vector<double> const& b )
{
  if( a.size()!=b.size() || a.empty() ) return 1e30;
  double m=0.; for( size_t i=0;i<a.size();++i ) m=std::max(m,std::fabs(a[i]-b[i]));
  return m;
}

int main()
{
  std::cout << "================================================================\n";
  std::cout << "  PSA rung 9: PHYSICAL binary breakthrough (competitive, non-isothermal)\n";
  std::cout << "  U=" << U_VEL << " D=" << D_ax << " F=" << F_ph << " k=" << k_ldf
            << "  b01=" << b01_L << " b02=" << b02_L << " (selectivity " << b01_L/b02_L << ")\n";
  std::cout << "  dH1=" << dH1 << " dH2=" << dH2 << " hw=" << hw
            << "  c0=(" << c0_1 << "," << c0_2 << ") tau=" << tau_in << " T_end=" << T_end << "\n";
  std::cout << "================================================================\n";

  BT coarse = run_BT( FFDom::CGL, "CGL", 5, 6, false );  // self-convergence reference
  BT fine   = run_BT( FFDom::CGL, "CGL", 5, 8, true  );  // primary + gnuplot files

  double d1 = curve_diff( coarse.c1o, fine.c1o );
  double d2 = curve_diff( coarse.c2o, fine.c2o );
  double dT = curve_diff( coarse.To,  fine.To  );

  std::cout << "\n==================== rung 9 summary ====================\n";
  std::cout << std::scientific << std::setprecision(3);
  std::cout << "  coarse(5x6): conv=" << (coarse.converged?"y":"n")
            << " mb1=" << coarse.mb1 << " mb2=" << coarse.mb2 << " dT=" << (coarse.Tmax-T0_ref) << "\n";
  std::cout << "  fine  (5x8): conv=" << (fine.converged?"y":"n")
            << " mb1=" << fine.mb1 << " mb2=" << fine.mb2 << " dT=" << (fine.Tmax-T0_ref) << "\n";
  std::cout << "  self-conv  max|c1_6-c1_8|=" << d1 << " max|c2|=" << d2 << " max|T|=" << dT << "\n";
  std::cout << "  roll-up(c2 outlet overshoot)=" << (fine.c2max_out>c0_2*1.02?"yes":"no")
            << "  breakthrough order t1>t2=" << (fine.t1_bt>fine.t2_bt?"yes":"no") << "\n";

  bool ok = fine.ok && coarse.converged && (d1<2e-2)&&(d2<2e-2)&&(dT<2e-2);
  std::cout << "\n  Overall: " << ( ok ? "ALL PASS (resolved, conservative, physical)" : "SOME FAILED" ) << "\n";
  return ok ? 0 : 1;
}
