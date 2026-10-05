// OCFE_PDE30_solve2.cpp   ---  PSA case study, RUNG 6 (pressure-driven momentum, staged)
// ===========================================================================
// First pressure-momentum rung: an explicit PRESSURE state P with an ideal-gas
// EOS, Darcy's law in PRESSURE, and the total continuity.  Isothermal, single
// component, MMS-validated -- the smallest form that introduces the
// P / EOS / Darcy / continuity machinery in isolation.
//
//   continuity: dc/dt + d(u c)/dz - D d2c/dz2 + F dq/dt = s_c
//   LDF:        dq/dt - k ( qs b c/(1+b c) - q )         = s_q
//   Darcy:      u + (kappa/mu) dP/dz = 0          (algebraic; velocity from dP/dz)
//   EOS:        P - Rg c = 0                       (algebraic; isothermal ideal gas)
//
// For an isothermal single component P = Rg c, so pressure is tied to c by the
// EOS -- but this is the correct physics, not a degeneracy: the transient
// pressure is real (P evolves with c through continuity) and u is driven by dP/dz
// through the explicit P state, exactly the formulation that becomes genuinely
// independent once a second component or temperature is added.  The manufactured
// u and P are DETERMINED by c via Darcy and the EOS (those two residuals hold by
// construction; only continuity and LDF carry sources).
//
// Open question the log answers: is the u <-> P <-> c chain INDEX 1 (algebraic
// u,P with spatial derivatives from collocation, as for rung 2's velocity), or
// does it push to higher index (the Pantelides trigger)?
//
// Manufactured:  c = 1 + Ac(1-z)^2 + Cc t,  q = Q0 + Qz z + Qt t,
//                P = Rg c,  u = 2 (kappa/mu) Rg Ac (1-z)  [ = -(kappa/mu) dP/dz ].
// Imposition: IC_STRONG (the PSA default).  P=Rg c is a value-slaved algebraic
// state (derivative-free defining constraint), so its C0-continuity claim is
// redundant and was making the STRONG/TRACE tau receiver block singular.  The
// value-slaved eligibility-gate fix (drop the implied claim in eligible blocks)
// removes that singularity, so IC_STRONG now recovers the exact root to ~1e-13.
// (The earlier IC_WEAK + SIGMA0=10 workaround, diagnosed in PDE30b_diag, is
// retired.)
// ===========================================================================

#include <iostream>
#include <cstdlib>
#include <iomanip>
#include <vector>
#include <cmath>
#include <string>

#ifndef OCFE_OCFESLV_HEADER
#define OCFE_OCFESLV_HEADER "ocfeslv.hpp"
#endif
#include OCFE_OCFESLV_HEADER

using namespace mc;

static double const Kperm = 0.5, Rg = 1.0;             // kappa/mu, and RT (lumped)
static double const D_ax = 0.1, F_ph = 0.5, k_ldf = 2.0, qs_L = 1.0, b_L = 1.0;
static double const Ac = 0.3, Cc = 0.2;
static double const Q0 = 0.4, Qz = 0.1, Qt = 0.3;
static double const kappa_g = Kperm*Rg;                 // u_man = 2 kappa_g Ac (1-z)

static inline double cM_exact( double z, double t ){ return 1.0 + Ac*(1.0-z)*(1.0-z) + Cc*t; }
static inline double qM_exact( double z, double t ){ return Q0 + Qz*z + Qt*t; }
static inline double PM_exact( double z, double t ){ return Rg*cM_exact(z,t); }
static inline double uM_exact( double z, double   ){ return 2.0*kappa_g*Ac*(1.0-z); }

static bool run_psa6( FFDom::TYPE coltype, std::string const& cname )
{
  std::cout << "\n================================================================\n";
  std::cout << "  PSA rung 6  pressure-driven Darcy (explicit P + EOS)  " << cname << "\n";
  std::cout << "================================================================\n";

  size_t const n_el = 3, n_nd = 6;

  FFGraph DAG;
  FFVar t = DAG.add_var( "t" );
  FFVar z = DAG.add_var( "z" );
  FFVar c = DAG.add_var( "c(t,z)" );
  FFVar q = DAG.add_var( "q(t,z)" );
  FFVar u = DAG.add_var( "u(t,z)" );
  FFVar P = DAG.add_var( "P(t,z)" );

  FFPartial OpP;

  FFVar cMan = 1.0 + Ac*(1.0-z)*(1.0-z) + Cc*t;
  FFVar qMan = Q0 + Qz*z + Qt*t;
  FFVar qstarMan = qs_L*b_L*cMan/( 1.0 + b_L*cMan );

  // sources (continuity convection term: d(u_man c_man)/dz = -2 kappa_g Ac (1 + 3 Ac (1-z)^2 + Cc t))
  FFVar s_c = Cc - 2.0*kappa_g*Ac*( 1.0 + 3.0*Ac*(1.0-z)*(1.0-z) + Cc*t ) - D_ax*( 2.0*Ac ) + F_ph*Qt;
  FFVar s_q = Qt - k_ldf*( qstarMan - qMan );
  FFVar g_c_in = 2.0*kappa_g*Ac*( 1.0 + Ac + Cc*t ) + 2.0*D_ax*Ac;   // u_man(0) c_man(0) - D dz c_man(0)

  FFVar qstar = qs_L*b_L*c/( 1.0 + b_L*c );
  FFVar CONT  = OpP( c, t ) + u*OpP( c, z ) + c*OpP( u, z ) - D_ax*OpP( OpP( c, z ), z ) + F_ph*OpP( q, t ) - s_c;
  FFVar LDF_q = OpP( q, t ) - k_ldf*( qstar - q ) - s_q;
  FFVar DARCY = u + Kperm*OpP( P, z );                 // algebraic
  FFVar EOS   = P - Rg*c;                              // algebraic
  FFVar IC_c  = c - ( 1.0 + Ac*(1.0-z)*(1.0-z) );
  FFVar IC_q  = q - ( Q0 + Qz*z );
  FFVar BC_Lc = u*c - D_ax*OpP( c, z ) - g_c_in;       // Danckwerts inlet
  FFVar BC_Uc = OpP( c, z );                            // zero-gradient outlet

  OCFESLV oc( &DAG );
  oc.add_domain( t, FFDom( 0., 1., n_el, coltype, n_nd ) );
  oc.add_domain( z, FFDom( 0., 1., n_el, coltype, n_nd ) );
  oc.add_state ( c, {t,z} );
  oc.add_state ( q, {t,z} );
  oc.add_state ( u, {t,z} );
  oc.add_state ( P, {t,z} );
  oc.update_ref( c, [&]( OCFESLV::t_Coord const& cr ){ return cM_exact( cr.at(z), cr.at(t) ); } );
  oc.update_ref( q, [&]( OCFESLV::t_Coord const& cr ){ return qM_exact( cr.at(z), cr.at(t) ); } );
  oc.update_ref( u, [&]( OCFESLV::t_Coord const& cr ){ return uM_exact( cr.at(z), cr.at(t) ); } );
  oc.update_ref( P, [&]( OCFESLV::t_Coord const& cr ){ return PM_exact( cr.at(z), cr.at(t) ); } );
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
  // IC_STRONG (the PSA default).  For an isothermal single component, P=Rg c is a
  // value-slaved algebraic state: P's interface-continuity claims are an exact
  // scalar multiple of c's.  Pre-fix, the tau-based impositions (STRONG/TRACE)
  // built a frozen trace receiver block the redundancy pre-pass flagged but did
  // not de-singularize -> tier-0 singular -> LM floor ~1e-3 (see PDE30b_diag).
  // The value-slaved eligibility-gate fix drops that implied claim in this
  // (eligible, parabolic) block, removing the singularity, so IC_STRONG now
  // recovers the exact polynomial root to ~1e-13.  sigma0=10 is retained for
  // basin robustness (the slaved u makes convection quadratic-gradient, as in
  // rung 2); it is no longer load-bearing for the singularity.
  oc.options.INTERFACE.IMPOSITION  = OCFESLV::Options::IC_STRONG;
  // 2026-09-09: PSA_IMP overrides the imposition so the gated solve can be run in the
  // modes this driver does not otherwise test.  MEASURED: OCFE_PDE30 gates IC_STRONG
  // only and fails in IC_TRACE with the claim drop unlocked -- a real defect the corpus
  // could not see.  Env-only; unset leaves the driver's own choice untouched.
  if( char const* v = std::getenv( "PSA_IMP" ) ){
    std::string const m_( v );
    if     ( m_ == "WEAK"   ) oc.options.INTERFACE.IMPOSITION = OCFESLV::Options::IC_WEAK;
    else if( m_ == "TRACE"  ) oc.options.INTERFACE.IMPOSITION = OCFESLV::Options::IC_TRACE;
    else if( m_ == "STRONG" ) oc.options.INTERFACE.IMPOSITION = OCFESLV::Options::IC_STRONG;
  }
  oc.options.INTERFACE.SAT_SIGMA0       = 10.0;
  oc.options.DISPLAY_LEVEL    = 1;   // show classification/index for this frontier rung

  if( !oc.setup() ){ std::cerr << "ERROR: setup failed\n"; return false; }

  size_t const nVar = oc.n_colloc_sta(), nEqn = oc.n_colloc_eqn();
  std::cout << "  type=" << OCFESLV::pde_type_name( oc.pde_type().type )
            << " nVar=" << nVar << " nEqn=" << nEqn << " nTrace=" << oc.n_colloc_trace()
            << " square=" << (nVar==nEqn?"yes":"NO") << "\n";
  if( nVar != nEqn ){ std::cerr << "ERROR: not square\n"; return false; }

  std::vector<double> varInit, inpInit;
  if( !oc.init( varInit, inpInit, nullptr ) ){ std::cerr << "ERROR: init failed\n"; return false; }

  std::vector<double> res( nEqn, 0. );
  if( oc.eval( res.data(), nullptr, varInit.data(), nullptr, nullptr ) ){
    double rm=0.; for( double v : res ) rm = std::max( rm, std::fabs(v) );
    std::cout << "  [A] residual at manufactured exact: max|r|="
              << std::scientific << std::setprecision(4) << rm
              << "  " << ( rm<1e-9 ? "PASS" : "(not exact)" ) << "\n";
  }

  std::vector<double> xv = varInit;
  for( size_t i=0;i<xv.size();++i ) xv[i] += 0.05*std::sin( 0.7*double(i) + 0.2 );
  OCFESLV::SolveReport rep = oc.solve( xv.data(), nullptr, nullptr );
  std::cout << "  [B] solve: converged=" << (rep.converged?"yes":"no")
            << " iters=" << rep.iterations
            << " final|r|=" << std::scientific << std::setprecision(4) << rep.final_residual << "\n";
  if( !rep.converged ){ std::cerr << "ERROR: solve did not converge\n"; return false; }

  double ce=0., qe=0., ue=0., Pe=0.;
  double const sg[4] = { 0.15, 0.35, 0.65, 0.85 };
  for( double zs : sg ) for( double ts : sg ){
    OCFESLV::t_Coord pt; pt[z]=zs; pt[t]=ts;
    ce = std::max( ce, std::fabs( oc.eval_colloc<double>( c, pt, xv.data(), nullptr, nullptr ) - cM_exact(zs,ts) ) );
    qe = std::max( qe, std::fabs( oc.eval_colloc<double>( q, pt, xv.data(), nullptr, nullptr ) - qM_exact(zs,ts) ) );
    ue = std::max( ue, std::fabs( oc.eval_colloc<double>( u, pt, xv.data(), nullptr, nullptr ) - uM_exact(zs,ts) ) );
    Pe = std::max( Pe, std::fabs( oc.eval_colloc<double>( P, pt, xv.data(), nullptr, nullptr ) - PM_exact(zs,ts) ) );
  }
  std::cout << std::scientific << std::setprecision(4)
            << "  [C] err c=" << ce << " q=" << qe << " u=" << ue << " P=" << Pe << "\n";

  bool ok = ( ce<1e-7 && qe<1e-7 && ue<1e-7 && Pe<1e-7 );
  std::cout << "  PSA6 " << cname << ": " << (ok?"PASS":"FAIL") << "\n";
  return ok;
}

int main()
{
  std::cout << "================================================================\n";
  std::cout << "  PSA rung 6: pressure-driven Darcy momentum (staged, index-1)\n";
  std::cout << "  explicit P state + ideal-gas EOS + Darcy + total continuity, MMS\n";
  std::cout << "  Kperm=" << Kperm << " Rg=" << Rg << "  IC_STRONG, SIGMA0=10\n";
  std::cout << "================================================================\n";
  bool ok = true;
  ok &= run_psa6( FFDom::CGL, "CGL" );
  ok &= run_psa6( FFDom::LGL, "LGL" );
  std::cout << "\n  Overall: " << ( ok ? "ALL PASS" : "SOME FAILED" ) << "\n";
  return ok ? 0 : 1;
}
