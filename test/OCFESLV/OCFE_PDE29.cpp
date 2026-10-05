// OCFE_PDE29_solve2.cpp   ---  PSA case study, RUNG 5b (non-isothermal breakthrough)
// ===========================================================================
// PHYSICAL non-isothermal adsorption sweep: rung 4's breakthrough + rung 5's
// energy balance and van't Hoff isotherm.  Regenerated bed at T0, cold feed at
// T0; the adsorption exotherm drives a THERMAL WAVE that travels ahead of the
// (retarded) concentration front, and b(T) lowers capacity in the hot zone.
//
//   mass:  dc/dt + u dc/dz - D d2c/dz2 + F dq/dt              = 0
//   LDF:   dq/dt - k ( qs b(T) c/(1+b(T) c) - q )             = 0
//   energy:Cp dT/dt + G dT/dz - lam d2T/dz2 - dH F dq/dt + hw (T-Tw) = 0
//   b(T) = b0 exp( beta ( 1/T - 1/T0 ) )
//   IC: c=q=0, T=T0 (regenerated).  inlet: Danckwerts mass (c_feed=c0(1-e^{-t/tau})),
//   thermal Danckwerts (T_feed=T0).  outlet: dc/dz=dT/dz=0.
//
// Validation (energy discretisation already MMS-proven in rung 5):
//   * concentration mass balance  int(c+F q)|_T dz = feed - effluent (converging);
//   * outlet curves c(1,t), T(1,t) self-converge under refinement;
//   * thermal-wave sanity: T rises above T0, stays bounded, relaxes behind front.
// Output: full c/q/T field over (z,t), outlet histories, spatial snapshots.
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

#ifndef PDE29_OUT_PREFIX
#define PDE29_OUT_PREFIX "OCFE_PDE29"
#endif

static double const U_VEL = 1.0, D_ax = 0.1, F_ph = 0.5, k_ldf = 2.0;
static double const qs_L = 1.0, b0_L = 1.0, beta_vh = 2.0, T0_ref = 1.0;
static double const Cp_e = 1.0, G_cv = 1.0, lam = 0.1, dH = 2.0, hw = 0.5, Tw = 1.0;
static double const c0_feed = 1.0, tau_in = 0.10, T_end = 3.0;

static inline double cfeed_exact( double t ){ return c0_feed*( 1.0 - std::exp( -t/tau_in ) ); }
static inline double qstar_T0( double cg ){ return qs_L*b0_L*cg/( 1.0 + b0_L*cg ); } // b(T0)=b0
static inline double feed_throughput()
{ return U_VEL*c0_feed*( T_end - tau_in*( 1.0 - std::exp( -T_end/tau_in ) ) ); }

struct NI {
  std::string tag;
  bool   converged=false, ok=false;
  size_t nVar=0;
  int    iters=0;
  double cmin=0., cmax=0., Tmin=0., Tmax=0., inv_T=0., effl=0., feed=0., massbal_rel=1.0;
  std::vector<double> cout_t, Tout_t;   // c(1,t_k), T(1,t_k) at fixed sample times
};

static std::vector<double> sample_times()
{ std::vector<double> ts; for( int k=0;k<=20;++k ) ts.push_back( double(k)*T_end/20.0 ); return ts; }

static NI run_NI( FFDom::TYPE coltype, std::string const& cname,
                  size_t n_el, size_t n_nd, bool write_files )
{
  NI R;
  R.tag = cname + " n_el=" + std::to_string(n_el) + " n_nd=" + std::to_string(n_nd);
  std::cout << "\n---- non-isothermal breakthrough  " << R.tag << " ----\n";

  FFGraph DAG;
  FFVar t = DAG.add_var( "t" );
  FFVar z = DAG.add_var( "z" );
  FFVar c = DAG.add_var( "c(t,z)" );
  FFVar q = DAG.add_var( "q(t,z)" );
  FFVar T = DAG.add_var( "T(t,z)" );

  FFPartial  OpP;
  FFIntegral OpI;

  FFVar cfeed = c0_feed*( 1.0 - exp( -t/tau_in ) );
  FFVar bT    = b0_L*exp( beta_vh*( 1.0/T - 1.0/T0_ref ) );
  FFVar qstar = qs_L*bT*c/( 1.0 + bT*c );

  FFVar PDE_c = OpP( c, t ) + U_VEL*OpP( c, z ) - D_ax*OpP( OpP( c, z ), z ) + F_ph*OpP( q, t );
  FFVar LDF_q = OpP( q, t ) - k_ldf*( qstar - q );
  FFVar ENE_T = Cp_e*OpP( T, t ) + G_cv*OpP( T, z ) - lam*OpP( OpP( T, z ), z )
              - dH*F_ph*OpP( q, t ) + hw*( T - Tw );
  FFVar IC_c  = c;
  FFVar IC_q  = q;
  FFVar IC_T  = T - T0_ref;
  FFVar BC_Lc = U_VEL*c - D_ax*OpP( c, z ) - U_VEL*cfeed;
  FFVar BC_Uc = OpP( c, z );
  FFVar BC_LT = G_cv*T - lam*OpP( T, z ) - G_cv*T0_ref;   // cold feed at T0
  FFVar BC_UT = OpP( T, z );

  FFVar Inv = OpI( c + F_ph*q, z );
  FFVar Eff = OpI( U_VEL*c, t );

  OCFESLV oc( &DAG );
  oc.add_domain( t, FFDom( 0., T_end, n_el, coltype, n_nd ) );
  oc.add_domain( z, FFDom( 0., 1.0,   n_el, coltype, n_nd ) );
  oc.add_state ( c, {t,z} );
  oc.add_state ( q, {t,z} );
  oc.add_state ( T, {t,z} );
  oc.update_ref( c, [&]( OCFESLV::t_Coord const& cr ){
    double zz=cr.at(z), tt=cr.at(t); return cfeed_exact(tt)*(1.0-0.5*zz); } );
  oc.update_ref( q, [&]( OCFESLV::t_Coord const& cr ){
    double zz=cr.at(z), tt=cr.at(t); return qstar_T0( cfeed_exact(tt)*(1.0-0.5*zz) ); } );
  oc.update_ref( T, [&]( OCFESLV::t_Coord const& cr ){
    double zz=cr.at(z), tt=cr.at(t); double qg=qstar_T0( cfeed_exact(tt)*(1.0-0.5*zz) );
    return T0_ref + dH*F_ph*qg/Cp_e; } );
  oc.set_evolution_domain( t );

  int const T_NO_LB = FFDom::ALL - FFDom::LB;
  int const Z_INT   = FFDom::ALL - FFDom::LB - FFDom::UB;

  oc.add_equation( PDE_c, {t,z}, {T_NO_LB,   Z_INT},      OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( LDF_q, {t,z}, {T_NO_LB,   FFDom::ALL}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( ENE_T, {t,z}, {T_NO_LB,   Z_INT},      OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( IC_c,  {t,z}, {FFDom::LB, FFDom::ALL}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INITIAL,  0 ) );
  oc.add_equation( IC_q,  {t,z}, {FFDom::LB, FFDom::ALL}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INITIAL,  0 ) );
  oc.add_equation( IC_T,  {t,z}, {FFDom::LB, FFDom::ALL}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INITIAL,  0 ) );
  oc.add_equation( BC_Lc, {t,z}, {T_NO_LB,   FFDom::LB},  OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );
  oc.add_equation( BC_Uc, {t,z}, {T_NO_LB,   FFDom::UB},  OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );
  oc.add_equation( BC_LT, {t,z}, {T_NO_LB,   FFDom::LB},  OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );
  oc.add_equation( BC_UT, {t,z}, {T_NO_LB,   FFDom::UB},  OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );

  oc.add_output( Inv, {t}, {T_end} );
  oc.add_output( Eff, {z}, {1.0} );

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
  // SOLVE_MARCHING defaults true (auto-march when an evolution domain + a differential state exist);
  // set = false to force a monolithic solve.  Both modes now read identically via the trajectory-
  // routed eval_colloc() (fields) and val_functions() (output functionals).
  // oc.options.SOLVE.MARCHING = false;

  if( !oc.setup() ){ std::cerr << "ERROR: setup failed\n"; return R; }
  R.nVar = oc.n_colloc_sta();
  size_t const nEqn = oc.n_colloc_eqn(), nFct = oc.n_colloc_fct();
  if( R.nVar != nEqn ){ std::cerr << "ERROR: not square\n"; return R; }

  std::vector<double> varInit, inpInit;
  if( !oc.init( varInit, inpInit, nullptr ) ){ std::cerr << "ERROR: init failed\n"; return R; }
  std::vector<double> xv = varInit;
  OCFESLV::SolveReport rep = oc.solve( xv.data(), nullptr, nullptr );
  R.converged = rep.converged; R.iters = rep.iterations;
  std::cout << "  nVar=" << R.nVar << "  [B] converged=" << (rep.converged?"yes":"no")
            << " iters=" << rep.iterations
            << " final|r|=" << std::scientific << std::setprecision(3) << rep.final_residual << "\n";
  if( !rep.converged ) return R;

  auto ts = sample_times();
  R.cout_t.reserve(ts.size()); R.Tout_t.reserve(ts.size());
  for( double tt : ts ){
    OCFESLV::t_Coord pt; pt[z]=1.0; pt[t]=tt;
    R.cout_t.push_back( oc.eval_colloc<double>( c, pt, xv.data(), nullptr, nullptr ) );
    R.Tout_t.push_back( oc.eval_colloc<double>( T, pt, xv.data(), nullptr, nullptr ) );
  }
  R.cmin=1e30; R.cmax=-1e30; R.Tmin=1e30; R.Tmax=-1e30;
  for( int iz=0;iz<=40;++iz ) for( int it=0;it<=40;++it ){
    OCFESLV::t_Coord pt; pt[z]=double(iz)/40.0; pt[t]=double(it)/40.0*T_end;
    double cv = oc.eval_colloc<double>( c, pt, xv.data(), nullptr, nullptr );
    double Tv = oc.eval_colloc<double>( T, pt, xv.data(), nullptr, nullptr );
    R.cmin=std::min(R.cmin,cv); R.cmax=std::max(R.cmax,cv);
    R.Tmin=std::min(R.Tmin,Tv); R.Tmax=std::max(R.Tmax,Tv);
  }

  // Output functionals from the framework's stored result (val_functions()) -- correct for BOTH
  // monolithic and marching.  Under marching the evolution-direction integral Eff = int_0^T u c dt
  // is the sum over windows accumulated during the march (see _valFct); a re-eval on xv would see
  // only the LAST window and undercount the effluent -> the mass-balance failure.  This is the
  // functional analogue of the eval_colloc read consolidation.
  auto const& fct = oc.val_functions();
  if( fct.size() < nFct ){ std::cerr << "ERROR: output functionals unavailable\n"; return R; }
  R.inv_T = fct[ oc.row_fct(0)-nEqn ];
  R.effl  = fct[ oc.row_fct(1)-nEqn ];
  R.feed  = feed_throughput();
  R.massbal_rel = std::fabs( R.inv_T - ( R.feed - R.effl ) ) / std::max(R.feed,1e-30);

  std::cout << std::scientific << std::setprecision(4)
            << "  c-range=[" << R.cmin << "," << R.cmax << "]  T-range=[" << R.Tmin << "," << R.Tmax << "]\n"
            << "  massbal_rel=" << R.massbal_rel << "  peak dT=" << (R.Tmax-T0_ref) << "\n";

  bool csane = (R.cmin > -0.05) && (R.cmax < c0_feed+0.05) && (R.cout_t.back() > 0.3*c0_feed);
  bool Tsane = (R.Tmax > T0_ref+0.02) && (R.Tmax < T0_ref+2.0) && (R.Tmin > T0_ref-0.05);
  R.ok = R.converged && (R.massbal_rel < 3e-3) && csane && Tsane;

  if( write_files ){
    // Full c/q/T field over (z,t): gnuplot pm3d (blank line between z-blocks).
    std::string ff = std::string(PDE29_OUT_PREFIX)+"_"+cname+"_field.out";
    std::ofstream of( ff );
    of << "# z  t  c  q  T   (splot u 1:2:5 w pm3d for the thermal wave)\n";
    int const NG=41;
    for( int iz=0; iz<NG; ++iz ){
      double zz=double(iz)/(NG-1);
      for( int it=0; it<NG; ++it ){
        double tt=double(it)/(NG-1)*T_end; OCFESLV::t_Coord pt; pt[z]=zz; pt[t]=tt;
        of << std::setprecision(8) << zz << " " << tt << " "
           << oc.eval_colloc<double>( c, pt, xv.data(), nullptr, nullptr ) << " "
           << oc.eval_colloc<double>( q, pt, xv.data(), nullptr, nullptr ) << " "
           << oc.eval_colloc<double>( T, pt, xv.data(), nullptr, nullptr ) << "\n";
      }
      of << "\n";
    }
    // Outlet histories: mass and thermal breakthrough at z=1.
    std::string fo = std::string(PDE29_OUT_PREFIX)+"_"+cname+"_outlet.out";
    std::ofstream oo( fo );
    oo << "# t  c(1,t)  T(1,t)  q(1,t)  c_feed(t)\n";
    for( int kk=0;kk<=200;++kk ){
      double tt=double(kk)/200.0*T_end; OCFESLV::t_Coord pt; pt[z]=1.0; pt[t]=tt;
      oo << std::setprecision(8) << tt << " "
         << oc.eval_colloc<double>( c, pt, xv.data(), nullptr, nullptr ) << " "
         << oc.eval_colloc<double>( T, pt, xv.data(), nullptr, nullptr ) << " "
         << oc.eval_colloc<double>( q, pt, xv.data(), nullptr, nullptr ) << " " << cfeed_exact(tt) << "\n";
    }
    // Spatial snapshots: c/q/T vs z at several times (moving wave).
    std::string fs = std::string(PDE29_OUT_PREFIX)+"_"+cname+"_snapshots.out";
    std::ofstream os( fs );
    os << "# z  c  q  T   (blocks at t = 0.5,1.0,1.5,2.0,2.5,3.0)\n";
    for( double tt : { 0.5,1.0,1.5,2.0,2.5,3.0 } ){
      os << "# t=" << tt << "\n";
      for( int kk=0;kk<=120;++kk ){
        double zz=double(kk)/120.0; OCFESLV::t_Coord pt; pt[z]=zz; pt[t]=tt;
        os << std::setprecision(8) << zz << " "
           << oc.eval_colloc<double>( c, pt, xv.data(), nullptr, nullptr ) << " "
           << oc.eval_colloc<double>( q, pt, xv.data(), nullptr, nullptr ) << " "
           << oc.eval_colloc<double>( T, pt, xv.data(), nullptr, nullptr ) << "\n";
      }
      os << "\n";
    }
    std::cout << "  wrote " << ff << " , " << fo << " , " << fs << "\n";
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
  std::cout << "  PSA rung 5b: non-isothermal breakthrough + thermal wave\n";
  std::cout << "  u=" << U_VEL << " D=" << D_ax << " lam=" << lam << " F=" << F_ph
            << " k=" << k_ldf << " dH=" << dH << " hw=" << hw << " beta=" << beta_vh
            << "  c0=" << c0_feed << " tau=" << tau_in << " T=" << T_end << "\n";
  std::cout << "================================================================\n";

  NI coarse = run_NI( FFDom::CGL, "CGL", 5, 6, false );  // convergence reference
  NI fine   = run_NI( FFDom::CGL, "CGL", 5, 8, true  );  // primary + wave files

  double dc = curve_diff( coarse.cout_t, fine.cout_t );
  double dT = curve_diff( coarse.Tout_t, fine.Tout_t );

  std::cout << "\n==================== rung 5b summary ====================\n";
  std::cout << std::scientific << std::setprecision(3);
  std::cout << "  coarse(5x6): conv=" << (coarse.converged?"y":"n")
            << " massbal=" << coarse.massbal_rel << " peakdT=" << (coarse.Tmax-T0_ref) << "\n";
  std::cout << "  fine  (5x8): conv=" << (fine.converged?"y":"n")
            << " massbal=" << fine.massbal_rel << " peakdT=" << (fine.Tmax-T0_ref) << "\n";
  std::cout << "  self-conv  max|c6-c8|=" << dc << "  max|T6-T8|=" << dT << "\n";

  bool ok = fine.ok && coarse.converged && (dc < 5e-3) && (dT < 5e-3);
  std::cout << "\n  Overall: " << ( ok ? "ALL PASS (thermal wave resolved & conservative)" : "SOME FAILED" ) << "\n";
  return ok ? 0 : 1;
}
