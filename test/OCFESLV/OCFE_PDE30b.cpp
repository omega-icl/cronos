// OCFE_PDE30b_diag.cpp  ---  PSA rung 6 FAIL diagnosis
// ===========================================================================
// Same explicit-P / EOS / Darcy / continuity model as PDE30, run as a 2-axis
// sweep to localize the ~1e-3 recovery floor:
//
//   AXIS 1 (resolution): n_nd = 6,8,10 and n_el = 3,5 at IC_STRONG.
//     If the floor is DISCRETIZATION it falls with refinement.  The manufactured
//     c is degree-2 in z (exactly representable), [A] is machine-zero, so the
//     prediction is a FLAT (or rising, via operator conditioning) floor.
//
//   AXIS 2 (imposition): IC_WEAK / IC_TRACE / IC_STRONG at fixed 3x6.
//     If the singularity is the IC_STRONG redundancy from the implied P-continuity
//     (P=Rg c), a non-eliminating imposition (WEAK/TRACE) recovers where STRONG
//     does not.
//
// Reports [A] (manufactured residual, = discretization error) and [C] (recovery
// error) per cell.  No physics here -- pure structural diagnosis.
// ===========================================================================

#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <string>

#ifndef OCFE_OCFESLV_HEADER
#define OCFE_OCFESLV_HEADER "ocfeslv.hpp"
#endif
#include OCFE_OCFESLV_HEADER
using namespace mc;

static double const Kperm = 0.5, Rg = 1.0;
static double const D_ax = 0.1, F_ph = 0.5, k_ldf = 2.0, qs_L = 1.0, b_L = 1.0;
static double const Ac = 0.3, Cc = 0.2, Q0 = 0.4, Qz = 0.1, Qt = 0.3;
static double const kappa_g = Kperm*Rg;

static size_t g_nonconv = 0;   // 2026-09-09: non-converged solves across all rows

static inline double cM( double z, double t ){ return 1.0 + Ac*(1.0-z)*(1.0-z) + Cc*t; }
static inline double qM( double z, double   ){ return Q0 + Qz*z; }            // ref @ t handled below
static inline double qMt( double z, double t ){ return Q0 + Qz*z + Qt*t; }
static inline double PM( double z, double t ){ return Rg*cM(z,t); }
static inline double uM( double z, double   ){ return 2.0*kappa_g*Ac*(1.0-z); }

struct Cell { double Aerr=0., cerr=0., uerr=0., Perr=0.; bool conv=false; int iters=0; double finalr=0.; };

static Cell run( FFDom::TYPE coltype, int imposition, size_t n_el, size_t n_nd )
{
  Cell R;
  FFGraph DAG;
  FFVar t=DAG.add_var("t"), z=DAG.add_var("z");
  FFVar c=DAG.add_var("c(t,z)"), q=DAG.add_var("q(t,z)"), u=DAG.add_var("u(t,z)"), P=DAG.add_var("P(t,z)");
  FFPartial OpP;

  FFVar cMan = 1.0 + Ac*(1.0-z)*(1.0-z) + Cc*t;
  FFVar qMan = Q0 + Qz*z + Qt*t;
  FFVar qstarMan = qs_L*b_L*cMan/( 1.0 + b_L*cMan );
  FFVar s_c = Cc - 2.0*kappa_g*Ac*( 1.0 + 3.0*Ac*(1.0-z)*(1.0-z) + Cc*t ) - D_ax*( 2.0*Ac ) + F_ph*Qt;
  FFVar s_q = Qt - k_ldf*( qstarMan - qMan );
  FFVar g_c_in = 2.0*kappa_g*Ac*( 1.0 + Ac + Cc*t ) + 2.0*D_ax*Ac;

  FFVar qstar = qs_L*b_L*c/( 1.0 + b_L*c );
  FFVar CONT  = OpP(c,t) + u*OpP(c,z) + c*OpP(u,z) - D_ax*OpP(OpP(c,z),z) + F_ph*OpP(q,t) - s_c;
  FFVar LDF_q = OpP(q,t) - k_ldf*( qstar - q ) - s_q;
  FFVar DARCY = u + Kperm*OpP(P,z);
  FFVar EOS   = P - Rg*c;
  FFVar IC_c  = c - ( 1.0 + Ac*(1.0-z)*(1.0-z) );
  FFVar IC_q  = q - ( Q0 + Qz*z );
  FFVar BC_Lc = u*c - D_ax*OpP(c,z) - g_c_in;
  FFVar BC_Uc = OpP(c,z);

  OCFESLV oc( &DAG );
  oc.add_domain( t, FFDom( 0., 1., n_el, coltype, n_nd ) );
  oc.add_domain( z, FFDom( 0., 1., n_el, coltype, n_nd ) );
  oc.add_state( c, {t,z} ); oc.add_state( q, {t,z} ); oc.add_state( u, {t,z} ); oc.add_state( P, {t,z} );
  oc.update_ref( c, [&](OCFESLV::t_Coord const& cr){ return cM(cr.at(z),cr.at(t)); } );
  oc.update_ref( q, [&](OCFESLV::t_Coord const& cr){ return qMt(cr.at(z),cr.at(t)); } );
  oc.update_ref( u, [&](OCFESLV::t_Coord const& cr){ return uM(cr.at(z),cr.at(t)); } );
  oc.update_ref( P, [&](OCFESLV::t_Coord const& cr){ return PM(cr.at(z),cr.at(t)); } );
  oc.set_evolution_domain( t );

  int const T_NO_LB = FFDom::ALL - FFDom::LB;
  int const Z_INT   = FFDom::ALL - FFDom::LB - FFDom::UB;
  oc.add_equation( CONT,  {t,z}, {T_NO_LB,   Z_INT},      OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( LDF_q, {t,z}, {T_NO_LB,   FFDom::ALL}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( DARCY, {t,z}, {FFDom::ALL,FFDom::ALL}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( EOS,   {t,z}, {FFDom::ALL,FFDom::ALL}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( IC_c,  {t,z}, {FFDom::LB, FFDom::ALL}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INITIAL,  0 ) );
  oc.add_equation( IC_q,  {t,z}, {FFDom::LB, FFDom::ALL}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INITIAL,  0 ) );
  oc.add_equation( BC_Lc, {t,z}, {T_NO_LB,   FFDom::LB},  OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );
  oc.add_equation( BC_Uc, {t,z}, {T_NO_LB,   FFDom::UB},  OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );

  oc.options.REDUCE.ORDER     = OCFESLV::Options::RED_FULL;
  oc.options.CLASSIFY.MODE         = OCFESLV::Options::CLASS_AUTO;
  oc.options.INTERFACE.TYPE   = OCFESLV::Options::IC_AUTO;
  oc.options.INTERFACE.IMPOSITION  = static_cast<OCFESLV::Options::ImpositionType>( imposition );
  oc.options.INTERFACE.SAT_SIGMA0       = 10.0;
  oc.options.DISPLAY_LEVEL    = 0;

  if( !oc.setup() ){ std::cerr << "  setup failed\n"; return R; }
  size_t const nVar = oc.n_colloc_sta(), nEqn = oc.n_colloc_eqn();
  if( nVar != nEqn ){ std::cerr << "  not square\n"; return R; }

  std::vector<double> varInit, inpInit;
  if( !oc.init( varInit, inpInit, nullptr ) ){ std::cerr << "  init failed\n"; return R; }
  std::vector<double> res( nEqn, 0. );
  if( oc.eval( res.data(), nullptr, varInit.data(), nullptr, nullptr ) ){
    for( double v : res ) R.Aerr = std::max( R.Aerr, std::fabs(v) );
  }
  std::vector<double> xv = varInit;
  for( size_t i=0;i<xv.size();++i ) xv[i] += 0.05*std::sin( 0.7*double(i)+0.2 );
  OCFESLV::SolveReport rep = oc.solve( xv.data(), nullptr, nullptr );
  R.conv=rep.converged;
  // 2026-09-09: record non-convergence globally.  MEASURED: with the claim drop unlocked this
  // driver reports 3 non-converged marching windows and still exits 0, because main() returned
  // 0 unconditionally and no caller inspected R.conv.  A driver that cannot fail cannot gate.
  if( !rep.converged ) g_nonconv++; R.iters=rep.iterations; R.finalr=rep.final_residual;

  double const sg[4]={0.15,0.35,0.65,0.85};
  for( double zs : sg ) for( double ts : sg ){
    OCFESLV::t_Coord pt; pt[z]=zs; pt[t]=ts;
    R.cerr=std::max(R.cerr,std::fabs(oc.eval_colloc<double>(c,pt,xv.data(),nullptr,nullptr)-cM(zs,ts)));
    R.uerr=std::max(R.uerr,std::fabs(oc.eval_colloc<double>(u,pt,xv.data(),nullptr,nullptr)-uM(zs,ts)));
    R.Perr=std::max(R.Perr,std::fabs(oc.eval_colloc<double>(P,pt,xv.data(),nullptr,nullptr)-PM(zs,ts)));
  }
  return R;
}

static char const* imp_name( int imp )
{
  if( imp==OCFESLV::Options::IC_WEAK )   return "WEAK ";
  if( imp==OCFESLV::Options::IC_STRONG ) return "STRONG";
  if( imp==OCFESLV::Options::IC_TRACE )  return "TRACE";
  return "?";
}

static void row( FFDom::TYPE ct, char const* ctn, int imp, size_t nel, size_t nnd )
{
  Cell R = run( ct, imp, nel, nnd );
  std::cout << "  " << ctn << "  " << std::setw(6) << imp_name(imp)
            << "  n_el=" << nel << " n_nd=" << nnd
            << "  [A]=" << std::scientific << std::setprecision(2) << R.Aerr
            << "  conv=" << (R.conv?"y":"n") << " it=" << R.iters
            << " |r|=" << R.finalr
            << "  [C] c=" << R.cerr << " u=" << R.uerr << " P=" << R.Perr
            << "  " << ( (R.cerr<1e-7&&R.uerr<1e-7&&R.Perr<1e-7) ? "PASS" : "fail" ) << "\n";
}

int main()
{
  std::cout << "================================================================\n";
  std::cout << "  PSA rung 6 FAIL diagnosis: resolution x imposition\n";
  std::cout << "  explicit P + EOS (P=Rg c) is exact-rep => [A] machine-zero;\n";
  std::cout << "  if [C] floor is discretization it falls with n_nd/n_el.\n";
  std::cout << "================================================================\n";

  std::cout << "\n--- AXIS 1: resolution (LGL, IC_STRONG) -- does the floor fall? ---\n";
  row( FFDom::LGL, "LGL", OCFESLV::Options::IC_STRONG, 3, 6  );
  row( FFDom::LGL, "LGL", OCFESLV::Options::IC_STRONG, 3, 8  );
  row( FFDom::LGL, "LGL", OCFESLV::Options::IC_STRONG, 3, 10 );
  row( FFDom::LGL, "LGL", OCFESLV::Options::IC_STRONG, 5, 6  );
  row( FFDom::LGL, "LGL", OCFESLV::Options::IC_STRONG, 8, 6  );

  std::cout << "\n--- AXIS 2: imposition (LGL, 3x6) -- does STRONG alone fail? ---\n";
  row( FFDom::LGL, "LGL", OCFESLV::Options::IC_WEAK,   3, 6 );
  row( FFDom::LGL, "LGL", OCFESLV::Options::IC_TRACE,  3, 6 );
  row( FFDom::LGL, "LGL", OCFESLV::Options::IC_STRONG, 3, 6 );

  std::cout << "\n  Read: flat/rising [C] across AXIS 1 => not discretization.\n";
  std::cout << "        WEAK/TRACE PASS where STRONG fails => IC_STRONG P-continuity redundancy.\n";

  // 2026-09-09: the verdict.  Every row above is a real solve; a non-converged one is a
  // failure whatever the accuracy columns say.
  if( g_nonconv ){
    std::cout << "\n  Overall: " << g_nonconv << " SOLVE(S) DID NOT CONVERGE -- FAIL\n";
    return 1;
  }
  std::cout << "\n  Overall: all solves converged -- PASS\n";
  return 0;
}
