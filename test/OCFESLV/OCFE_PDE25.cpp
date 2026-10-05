// OCFE_PDE25_solve2.cpp   ---  PSA case study, RUNG 2 (validation + diagnostic)
// ===========================================================================
// Isothermal, single-component packed bed with LINEAR-DARCY velocity.
//   gas:  dc/dt + u dc/dz + c du/dz - D d2c/dz2 + F dq/dt = s_c(z,t)
//   LDF:  dq/dt - k ( q*(c) - q )                          = s_q(z,t)
//   mom:  u + K dc/dz = 0   (algebraic, index-1, u slaved to dc/dz)
//
// Manufactured: c = 1 + (Az/2)(1-z)^2 + Ct t,  q = Q0+Qz z+Qt t,  u = K Az (1-z).
//
// NONLINEAR STRUCTURE / MULTIPLE ROOTS.  Slaving u = -K c_z turns the convection
// into  -K (c_z)^2 - K c c_zz : a quadratic-gradient (Hamilton-Jacobi-like) term,
// so the discrete system has a SECOND root ~0.44 from the manufactured one (both
// satisfy the residual to machine zero).  The manufactured root's basin shrinks
// on finer/LGL grids; an aggressive initial guess can land in the spurious basin
// under IC_WEAK.  Because BOTH roots are continuous, the SAT penalties (SAT_SIGMA0,
// SAT_SIGMA1) are zero at either root and cannot select between them -- penalty
// tuning does NOT fix it.  IC_STRONG / IC_TRACE are STRUCTURAL (Schur-eliminated /
// exact continuity rows): they reshape the Newton system, tighten the feasible set,
// and recover the physical root robustly.  Hence the validation sweep uses WEAK
// only where it is robust (CGL) and STRONG/TRACE for LGL; the diagnostics localise
// the basin vs perturbation and confirm SAT_SIGMA0 is (near-)inert.
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

#ifndef PDE25_OUT_PREFIX
#define PDE25_OUT_PREFIX "OCFE_PDE25"
#endif

static double const D_ax = 0.1, F_ph = 0.5, k_ldf = 2.0, qs_L = 1.0, b_L = 1.0;
static double const K_dar = 2.0;
static double const Az = 0.5, Ct = 0.2;
static double const Q0 = 0.4, Qz = 0.1, Qt = 0.3;

static inline double cM_exact( double z, double t ){ return 1.0 + 0.5*Az*(1.0-z)*(1.0-z) + Ct*t; }
static inline double qM_exact( double z, double t ){ return Q0 + Qz*z + Qt*t; }
static inline double uM_exact( double z, double   ){ return K_dar*Az*(1.0-z); }

struct Cfg {
  FFDom::TYPE coltype = FFDom::CGL;
  std::string cname   = "CGL";
  bool        conservative = false;
  OCFESLV::Options::ImpositionType imp = OCFESLV::Options::IC_WEAK;
  std::string impname = "WEAK";
  int         fact    = 0;            // 0=SUPERLU, 2=SPQR (needs -DCRONOS__WITH_SPQR)
  std::string factname = "SUPERLU";
  OCFESLV::Options::ReductionType red = OCFESLV::Options::RED_FULL;   // 2026-09-19: the matrix's second axis
  std::string redname = "RED_FULL";
  size_t      n_el = 3, n_nd = 6;
  double      sigma0 = 1.0, sigma1 = 1.0;
  double      pert   = 0.05;
  bool        write_file = false;
};

struct Result {
  std::string tag;
  size_t nVar=0, nEqn=0, nTrace=0;
  bool   square=false, converged=false, ok=false;
  int    iters=0;
  double sigma0=0., sigma1=0., pert=0., finalr=0., ce=0., qe=0., ue=0.;
  std::vector<double> xv;
};

static Result run_psa2( Cfg const& cfg )
{
  Result R;
  R.tag = cfg.cname + "/" + (cfg.conservative?"cons":"exp ") + "/" + cfg.impname + "/" + cfg.factname
        + ( cfg.red == OCFESLV::Options::RED_MAIN ? "/RED_MAIN" : "" );
  R.sigma0 = cfg.sigma0; R.sigma1 = cfg.sigma1; R.pert = cfg.pert;
  std::cout << "\n---- " << R.tag
            << "  [s0=" << cfg.sigma0 << " s1=" << cfg.sigma1 << " pert=" << cfg.pert << "] ----\n";

  FFGraph DAG;
  FFVar t = DAG.add_var( "t" );
  FFVar z = DAG.add_var( "z" );
  FFVar c = DAG.add_var( "c(t,z)" );
  FFVar q = DAG.add_var( "q(t,z)" );
  FFVar u = DAG.add_var( "u(t,z)" );

  FFPartial OpP;

  FFVar cMan  = 1.0 + 0.5*Az*(1.0-z)*(1.0-z) + Ct*t;
  FFVar qMan  = Q0 + Qz*z + Qt*t;
  FFVar uMan  = K_dar*Az*(1.0-z);
  FFVar dzc   = -Az*(1.0-z);
  FFVar qstarM = qs_L*b_L*cMan / ( 1.0 + b_L*cMan );
  FFVar s_c   = Ct + uMan*dzc + cMan*( -K_dar*Az ) - D_ax*Az + F_ph*Qt;
  FFVar s_q   = Qt - k_ldf*( qstarM - qMan );
  FFVar g_in  = K_dar*Az*( 1.0 + 0.5*Az + Ct*t ) + D_ax*Az;

  FFVar conv = cfg.conservative ? OpP( u*c, z )
                                : ( u*OpP( c, z ) + c*OpP( u, z ) );
  FFVar PDE_c = OpP( c, t ) + conv - D_ax*OpP( OpP( c, z ), z ) + F_ph*OpP( q, t ) - s_c;
  FFVar LDF_q = OpP( q, t ) - k_ldf*( qs_L*b_L*c/( 1.0 + b_L*c ) - q ) - s_q;
  FFVar MOM   = u + K_dar*OpP( c, z );
  FFVar IC_c  = c - ( 1.0 + 0.5*Az*(1.0-z)*(1.0-z) );
  FFVar IC_q  = q - ( Q0 + Qz*z );
  FFVar BC_L  = u*c - D_ax*OpP( c, z ) - g_in;
  FFVar BC_U  = OpP( c, z );

  OCFESLV oc( &DAG );
  oc.add_domain( t, FFDom( 0., 1., cfg.n_el, cfg.coltype, cfg.n_nd ) );
  oc.add_domain( z, FFDom( 0., 1., cfg.n_el, cfg.coltype, cfg.n_nd ) );
  oc.add_state ( c, {t,z} );
  oc.add_state ( q, {t,z} );
  oc.add_state ( u, {t,z} );
  oc.update_ref( c, [&]( OCFESLV::t_Coord const& cr ){ return cM_exact( cr.at(z), cr.at(t) ); } );
  oc.update_ref( q, [&]( OCFESLV::t_Coord const& cr ){ return qM_exact( cr.at(z), cr.at(t) ); } );
  oc.update_ref( u, [&]( OCFESLV::t_Coord const& cr ){ return uM_exact( cr.at(z), cr.at(t) ); } );
  oc.set_evolution_domain( t );

  int const T_NO_LB = FFDom::ALL - FFDom::LB;
  int const Z_INT   = FFDom::ALL - FFDom::LB - FFDom::UB;

  oc.add_equation( PDE_c, {t,z}, {T_NO_LB,   Z_INT},      OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( LDF_q, {t,z}, {T_NO_LB,   FFDom::ALL}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( MOM,   {t,z}, {FFDom::ALL,FFDom::ALL}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( IC_c,  {t,z}, {FFDom::LB, FFDom::ALL}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INITIAL,  0 ) );
  oc.add_equation( IC_q,  {t,z}, {FFDom::LB, FFDom::ALL}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INITIAL,  0 ) );
  oc.add_equation( BC_L,  {t,z}, {T_NO_LB,   FFDom::LB},  OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );
  oc.add_equation( BC_U,  {t,z}, {T_NO_LB,   FFDom::UB},  OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );

  oc.options.REDUCE.ORDER     = cfg.red;
  oc.options.CLASSIFY.MODE         = OCFESLV::Options::CLASS_AUTO;
  oc.options.INTERFACE.TYPE   = OCFESLV::Options::IC_AUTO;
  oc.options.INTERFACE.IMPOSITION  = cfg.imp;
  oc.options.INTERFACE.SAT_SIGMA0       = cfg.sigma0;
  oc.options.DISPLAY_LEVEL    = 0;
#if defined(CRONOS__WITH_SPQR)
  oc.options.SOLVE.FACTORIZATION = ( cfg.fact==2 ) ? OCFESLV::Options::SOLVE_SPQR
                                                   : OCFESLV::Options::SOLVE_SUPERLU;
#else
  oc.options.SOLVE.FACTORIZATION = OCFESLV::Options::SOLVE_SUPERLU;
#endif

  if( !oc.setup() ){ std::cerr << "ERROR: setup failed (" << R.tag << ")\n"; return R; }

  R.nVar = oc.n_colloc_sta();
  R.nEqn = oc.n_colloc_eqn();
  R.nTrace = oc.n_colloc_trace();
  R.square = ( R.nVar == R.nEqn );
  std::cout << "  type=" << OCFESLV::pde_type_name( oc.pde_type().type )
            << " nVar=" << R.nVar << " nEqn=" << R.nEqn << " nTrace=" << R.nTrace
            << " square=" << (R.square?"yes":"NO") << "\n";

  std::vector<double> varInit, inpInit;
  if( !oc.init( varInit, inpInit, nullptr ) ){
    std::cerr << "ERROR: init failed (" << R.tag << ")\n"; return R;
  }
  std::vector<double> res( R.nEqn, 0. );
  if( oc.eval( res.data(), nullptr, varInit.data(), nullptr, nullptr ) ){
    double rm=0.; for( double v : res ) rm = std::max( rm, std::fabs(v) );
    std::cout << "  [A] residual at manufactured exact: max|r|="
              << std::scientific << std::setprecision(4) << rm
              << "  " << ( rm<1e-9 ? "PASS" : "(not exact)" ) << "\n";
  }

  std::vector<double> xv = varInit;
  for( size_t i=0;i<xv.size();++i ) xv[i] += cfg.pert*std::sin( 0.7*double(i) + 0.2 );
  OCFESLV::SolveReport rep = oc.solve( xv.data(), nullptr, nullptr );
  R.converged = rep.converged; R.iters = rep.iterations; R.finalr = rep.final_residual;
  R.xv = xv;
  std::cout << "  [B] solve: converged=" << (rep.converged?"yes":"no")
            << " iters=" << rep.iterations
            << " final|r|=" << std::scientific << std::setprecision(4) << rep.final_residual << "\n";

  if( rep.converged ){
    double const sg[4] = { 0.15, 0.35, 0.65, 0.85 };
    for( double zs : sg ) for( double ts : sg ){
      OCFESLV::t_Coord pt; pt[z]=zs; pt[t]=ts;
      R.ce = std::max( R.ce, std::fabs( oc.eval_colloc<double>( c, pt, xv.data(), nullptr, nullptr ) - cM_exact(zs,ts) ) );
      R.qe = std::max( R.qe, std::fabs( oc.eval_colloc<double>( q, pt, xv.data(), nullptr, nullptr ) - qM_exact(zs,ts) ) );
      R.ue = std::max( R.ue, std::fabs( oc.eval_colloc<double>( u, pt, xv.data(), nullptr, nullptr ) - uM_exact(zs,ts) ) );
    }
    std::cout << std::scientific << std::setprecision(4)
              << "  [C] err c=" << R.ce << " q=" << R.qe << " u=" << R.ue << "\n";
  }
  R.ok = R.square && R.converged && R.ce<1e-7 && R.qe<1e-7 && R.ue<1e-7;

  if( cfg.write_file && rep.converged ){
    std::string const fn = std::string(PDE25_OUT_PREFIX) + "_" + cfg.cname + "_" + cfg.impname + ".out";
    std::ofstream out( fn );
    out << "# z t c c_exact q q_exact u u_exact   (" << R.tag << ")\n";
    int const NG = 11;
    for( int iz=0; iz<NG; ++iz ){
      double const zz = double(iz)/double(NG-1);
      for( int it=0; it<NG; ++it ){
        double const tt = double(it)/double(NG-1);
        OCFESLV::t_Coord pt; pt[z]=zz; pt[t]=tt;
        out << std::setprecision(16)
            << zz << " " << tt << " "
            << oc.eval_colloc<double>( c, pt, xv.data(), nullptr, nullptr ) << " " << cM_exact(zz,tt) << " "
            << oc.eval_colloc<double>( q, pt, xv.data(), nullptr, nullptr ) << " " << qM_exact(zz,tt) << " "
            << oc.eval_colloc<double>( u, pt, xv.data(), nullptr, nullptr ) << " " << uM_exact(zz,tt) << "\n";
      }
      out << "\n";
    }
    std::cout << "  wrote " << fn << "\n";
  }

  std::cout << "  " << R.tag << ": " << (R.ok?"PASS":"FAIL") << "\n";
  return R;
}

int main()
{
  using O = OCFESLV::Options;
  std::cout << "================================================================\n";
  std::cout << "  PSA rung 2: validation (robust impositions) + basin/SAT diagnostics\n";
  std::cout << "  Default imposition = IC_WEAK.  n_el=3 n_nd=6, pert=0.05.\n";
  std::cout << "================================================================\n";

  // ============================ VALIDATION SWEEP =============================
  // WEAK is used only where robust (CGL); LGL uses STRUCTURAL impositions
  // (TRACE/STRONG), which recover the physical root regardless of initial guess.
  std::vector<Result> V;

  Cfg base; base.write_file = true;                         // CGL/WEAK
  Result rbase = run_psa2( base ); V.push_back( rbase );

  Cfg ccons = base; ccons.conservative = true; ccons.write_file = false;  // conservative form
  Result rcons = run_psa2( ccons ); V.push_back( rcons );
  if( rbase.converged && rcons.converged && rbase.xv.size()==rcons.xv.size() ){
    double dmax=0.; for( size_t i=0;i<rbase.xv.size();++i ) dmax = std::max( dmax, std::fabs(rbase.xv[i]-rcons.xv[i]) );
    std::cout << "\n[cons vs exp] max|dx| = " << std::scientific << std::setprecision(4) << dmax
              << "  -> " << ( dmax<1e-8 ? "EQUAL to solver tol" : "DIFFERENT" )
              << "  (bitwise diff is FAD term order; FP add non-associative)\n";
  }

  Cfg ctr = base; ctr.imp = O::IC_TRACE;  ctr.impname = "TRACE";  ctr.write_file = true;
  V.push_back( run_psa2( ctr ) );
  Cfg cst = base; cst.imp = O::IC_STRONG; cst.impname = "STRONG"; cst.write_file = true;
  V.push_back( run_psa2( cst ) );

#if defined(CRONOS__WITH_SPQR)
  Cfg cqr = base; cqr.fact = 2; cqr.factname = "SPQR";
  V.push_back( run_psa2( cqr ) );
#else
  std::cout << "\n[SPQR] skipped: build with -DCRONOS__WITH_SPQR (+ -lspqr -lcholmod -lsuitesparseconfig).\n";
#endif

  // LGL with structural imposition -> robust where LGL/WEAK is not.
  Cfg lt = base; lt.coltype = FFDom::LGL; lt.cname = "LGL"; lt.imp = O::IC_TRACE;  lt.impname = "TRACE";
  V.push_back( run_psa2( lt ) );
  Cfg ls = base; ls.coltype = FFDom::LGL; ls.cname = "LGL"; ls.imp = O::IC_STRONG; ls.impname = "STRONG";
  V.push_back( run_psa2( ls ) );

  std::cout << "\n==================== validation summary ====================\n";
  std::cout << std::left << std::setw(26) << "config"
            << std::right << std::setw(7) << "nVar" << std::setw(7) << "nTrc"
            << std::setw(7) << "it" << std::setw(12) << "final|r|"
            << std::setw(12) << "c err" << std::setw(12) << "u err"
            << std::setw(8) << "result" << "\n";
  bool all = true;
  for( auto const& r : V ){
    all &= r.ok;
    std::cout << std::left << std::setw(26) << r.tag
              << std::right << std::setw(7) << r.nVar << std::setw(7) << r.nTrace
              << std::setw(7) << r.iters << std::scientific << std::setprecision(3)
              << std::setw(12) << r.finalr << std::setw(12) << r.ce << std::setw(12) << r.ue
              << std::setw(8) << (r.ok?"PASS":"FAIL") << "\n";
  }

  // ===================== DIAGNOSTIC 1: basin vs perturbation =================
  // LGL/WEAK only: recovers the physical root for small perturbations, escapes to
  // the spurious root as the initial guess worsens.  This is the multiple-root /
  // basin signature, independent of penalties.
  std::cout << "\n=========== diagnostic 1: LGL/WEAK basin vs perturbation ===========\n";
  std::cout << std::left << std::setw(10) << "pert"
            << std::right << std::setw(7) << "it" << std::setw(13) << "final|r|"
            << std::setw(13) << "c err" << std::setw(10) << "result" << "\n";
  std::vector<Result> D1;
  for( double pe : { 0.002, 0.01, 0.05, 0.10 } ){
    Cfg d; d.coltype = FFDom::LGL; d.cname = "LGL"; d.pert = pe;
    D1.push_back( run_psa2( d ) );
  }
  for( auto const& r : D1 )
    std::cout << std::left << std::setw(10) << std::fixed << std::setprecision(3) << r.pert
              << std::right << std::setw(7) << r.iters << std::scientific << std::setprecision(3)
              << std::setw(13) << r.finalr << std::setw(13) << r.ce
              << std::setw(10) << (r.ok?"PASS":"FAIL") << "\n";

  // ===================== DIAGNOSTIC 2: SAT_SIGMA0 strength ===================
  // At the failing point (LGL/WEAK, pert=0.05), does a stronger C0 penalty move
  // the basin?  Expected near-inert: both roots are continuous, so the penalty is
  // zero at either and cannot select between them.
  std::cout << "\n=========== diagnostic 2: LGL/WEAK @ pert=0.05, SAT_SIGMA0 sweep ===========\n";
  std::cout << std::left << std::setw(10) << "sigma0"
            << std::right << std::setw(7) << "it" << std::setw(13) << "final|r|"
            << std::setw(13) << "c err" << std::setw(10) << "result" << "\n";
  std::vector<Result> D2;
  for( double s0 : { 1.0, 10.0, 100.0, 1000.0 } ){
    Cfg d; d.coltype = FFDom::LGL; d.cname = "LGL"; d.pert = 0.05; d.sigma0 = s0;
    D2.push_back( run_psa2( d ) );
  }
  for( auto const& r : D2 )
    std::cout << std::left << std::setw(10) << std::fixed << std::setprecision(1) << r.sigma0
              << std::right << std::setw(7) << r.iters << std::scientific << std::setprecision(3)
              << std::setw(13) << r.finalr << std::setw(13) << r.ce
              << std::setw(10) << (r.ok?"PASS":"FAIL") << "\n";

  // ---------------------------------------------------------------------------------------------------------
  // 2026-09-19 -- THE SYSTEMATIC MATRIX: 3 impositions x 2 reductions, one cell each, on the CGL/expanded form.
  //   Corpus census: only five drivers covered all three impositions AND both reductions, and PDE25 -- the
  //   hidden-c0 witness (rev296), one of only three models carrying that population -- covered RED_FULL only.
  //   A cell EXPECTED to differ gets a documented XFAIL naming the mechanism, never a loosened bar (the
  //   OCFE_scalar_nest4 discipline); a cell that unexpectedly passes is FLAGGED so the fix is noticed.
  // ---------------------------------------------------------------------------------------------------------
  std::cout << "\n---- SYSTEMATIC MATRIX (CGL, expanded form): imposition x reduction ----\n";
  {
    struct MCell { char const* imp; OCFESLV::Options::ImpositionType it;
                   char const* red; OCFESLV::Options::ReductionType rt; };
    MCell const cells[6] = {
      { "WEAK",   OCFESLV::Options::IC_WEAK,   "RED_FULL", OCFESLV::Options::RED_FULL },
      { "TRACE",  OCFESLV::Options::IC_TRACE,  "RED_FULL", OCFESLV::Options::RED_FULL },
      { "STRONG", OCFESLV::Options::IC_STRONG, "RED_FULL", OCFESLV::Options::RED_FULL },
      { "WEAK",   OCFESLV::Options::IC_WEAK,   "RED_MAIN", OCFESLV::Options::RED_MAIN },
      { "TRACE",  OCFESLV::Options::IC_TRACE,  "RED_MAIN", OCFESLV::Options::RED_MAIN },
      { "STRONG", OCFESLV::Options::IC_STRONG, "RED_MAIN", OCFESLV::Options::RED_MAIN } };
    bool matrix_ok = true;
    std::cout << std::left << std::setw(10) << "IMPOSITION" << std::setw(11) << "REDUCTION"
              << std::setw(7) << "CONV" << std::setw(8) << "nTRACE" << std::setw(13) << "final|r|"
              << std::setw(13) << "|c-c*|" << "VERDICT\n";
    for( auto const& mc : cells ){
      Cfg c; c.write_file = false; c.imp = mc.it; c.impname = mc.imp; c.red = mc.rt; c.redname = mc.red;
      Result r = run_psa2( c );
      bool const xfail = false;   // none expected: PDE25's rows keep their derivative under RED_MAIN
      std::string verdict = r.ok ? "PASS" : ( xfail ? "XFAIL (documented)" : "FAIL" );
      if( r.ok && xfail ) verdict = "PASS (unexpected: the documented defect is gone?)";
      if( !r.ok && !xfail ) matrix_ok = false;
      std::cout << std::left << std::setw(10) << mc.imp << std::setw(11) << mc.red
                << std::setw(7) << ( r.converged ? "yes" : "NO" ) << std::setw(8) << r.nTrace
                << std::scientific << std::setprecision(3) << std::setw(13) << r.finalr
                << std::setw(13) << r.ce << verdict << "\n";
    }
    std::cout << "  MATRIX: " << ( matrix_ok ? "PASS" : "FAIL" ) << "\n";
    all &= matrix_ok;
  }

  std::cout << "\n  Overall (validation): " << ( all ? "ALL PASS" : "SOME FAILED" ) << "\n";
  return all ? 0 : 1;
}
