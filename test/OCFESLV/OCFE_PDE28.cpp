// OCFE_PDE28_solve2.cpp   ---  PSA case study, RUNG 5 (non-isothermal, MMS)
// ===========================================================================
// Adds the ENERGY BALANCE and a temperature-dependent (van't Hoff) isotherm to
// the constant-velocity bed, closing the thermal loop
//     T -> b(T) -> q*(c,T) -> heat of adsorption -> T.
//
//   mass:  dc/dt + u dc/dz - D d2c/dz2 + F dq/dt              = s_c
//   LDF:   dq/dt - k ( qs b(T) c/(1+b(T) c) - q )             = s_q
//   energy:Cp dT/dt + G dT/dz - lam d2T/dz2 - dH F dq/dt + hw (T-Tw) = s_E
//   b(T) = b0 exp( beta ( 1/T - 1/T0 ) )      (exothermic: b decreases with T)
//
// The energy equation is a second parabolic transport equation (T is a regular
// dynamic state, like c) -- no new DAE-index difficulty.  Validated by the method
// of manufactured solutions to machine precision; analytic sources force the
// chosen polynomial fields, and both c and T use a (1-z)^2 profile so their
// zero-gradient outlets hold exactly.  Standard PSA config: IC_STRONG, SIGMA0=10.
//
// Manufactured:  c = 1 + Ac(1-z)^2 + Cc t,  T = 1 + At(1-z)^2 + Ctt t,
//                q = Q0 + Qz z + Qt t.
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

// transport / kinetics
static double const U_VEL = 1.0, D_ax = 0.1, F_ph = 0.5, k_ldf = 2.0;
// isotherm (van't Hoff Langmuir)
static double const qs_L = 1.0, b0_L = 1.0, beta_vh = 1.0, T0_ref = 1.0;
// energy
static double const Cp_e = 1.0, G_cv = 1.0, lam = 0.1, dH = 0.5, hw = 0.2, Tw = 1.0;
// manufactured coefficients
static double const Ac = 0.3, Cc = 0.2;
static double const At = 0.2, Ctt = 0.15;
static double const Q0 = 0.4, Qz = 0.1, Qt = 0.3;

static inline double cM_exact( double z, double t ){ return 1.0 + Ac*(1.0-z)*(1.0-z) + Cc*t; }
static inline double qM_exact( double z, double t ){ return Q0 + Qz*z + Qt*t; }
static inline double TM_exact( double z, double t ){ return 1.0 + At*(1.0-z)*(1.0-z) + Ctt*t; }

static bool run_psa5( FFDom::TYPE coltype, std::string const& cname )
{
  std::cout << "\n================================================================\n";
  std::cout << "  PSA rung 5  non-isothermal (energy + van't Hoff)  " << cname << "\n";
  std::cout << "================================================================\n";

  size_t const n_el = 3, n_nd = 6;

  FFGraph DAG;
  FFVar t = DAG.add_var( "t" );
  FFVar z = DAG.add_var( "z" );
  FFVar c = DAG.add_var( "c(t,z)" );
  FFVar q = DAG.add_var( "q(t,z)" );
  FFVar T = DAG.add_var( "T(t,z)" );

  FFPartial OpP;

  // Manufactured fields and the analytic sources / inlet fluxes.
  FFVar cMan = 1.0 + Ac*(1.0-z)*(1.0-z) + Cc*t;
  FFVar TMan = 1.0 + At*(1.0-z)*(1.0-z) + Ctt*t;
  FFVar qMan = Q0 + Qz*z + Qt*t;
  FFVar bTMan    = b0_L*exp( beta_vh*( 1.0/TMan - 1.0/T0_ref ) );
  FFVar qstarMan = qs_L*bTMan*cMan/( 1.0 + bTMan*cMan );

  // derivatives of the manufactured fields:
  //   dc/dt=Cc, dc/dz=-2Ac(1-z), d2c/dz2=2Ac ;  dT/dt=Ctt, dT/dz=-2At(1-z), d2T/dz2=2At
  FFVar s_c = Cc + U_VEL*( -2.0*Ac*(1.0-z) ) - D_ax*( 2.0*Ac ) + F_ph*Qt;
  FFVar s_q = Qt - k_ldf*( qstarMan - qMan );
  FFVar s_E = Cp_e*Ctt + G_cv*( -2.0*At*(1.0-z) ) - lam*( 2.0*At ) - dH*F_ph*Qt + hw*( TMan - Tw );
  FFVar g_c_in = U_VEL*( 1.0 + Ac + Cc*t ) + 2.0*D_ax*Ac;   // u c(0,t) - D dz c(0,t)
  FFVar g_T_in = G_cv*( 1.0 + At + Ctt*t ) + 2.0*lam*At;     // G T(0,t) - lam dz T(0,t)

  // Residuals (state T enters the isotherm -> the non-isothermal coupling).
  FFVar bT    = b0_L*exp( beta_vh*( 1.0/T - 1.0/T0_ref ) );
  FFVar qstar = qs_L*bT*c/( 1.0 + bT*c );
  FFVar PDE_c = OpP( c, t ) + U_VEL*OpP( c, z ) - D_ax*OpP( OpP( c, z ), z ) + F_ph*OpP( q, t ) - s_c;
  FFVar LDF_q = OpP( q, t ) - k_ldf*( qstar - q ) - s_q;
  FFVar ENE_T = Cp_e*OpP( T, t ) + G_cv*OpP( T, z ) - lam*OpP( OpP( T, z ), z )
              - dH*F_ph*OpP( q, t ) + hw*( T - Tw ) - s_E;
  FFVar IC_c  = c - ( 1.0 + Ac*(1.0-z)*(1.0-z) );
  FFVar IC_q  = q - ( Q0 + Qz*z );
  FFVar IC_T  = T - ( 1.0 + At*(1.0-z)*(1.0-z) );
  FFVar BC_Lc = U_VEL*c - D_ax*OpP( c, z ) - g_c_in;
  FFVar BC_Uc = OpP( c, z );
  FFVar BC_LT = G_cv*T - lam*OpP( T, z ) - g_T_in;
  FFVar BC_UT = OpP( T, z );

  OCFESLV oc( &DAG );
  oc.add_domain( t, FFDom( 0., 1., n_el, coltype, n_nd ) );
  oc.add_domain( z, FFDom( 0., 1., n_el, coltype, n_nd ) );
  oc.add_state ( c, {t,z} );
  oc.add_state ( q, {t,z} );
  oc.add_state ( T, {t,z} );
  oc.update_ref( c, [&]( OCFESLV::t_Coord const& cr ){ return cM_exact( cr.at(z), cr.at(t) ); } );
  oc.update_ref( q, [&]( OCFESLV::t_Coord const& cr ){ return qM_exact( cr.at(z), cr.at(t) ); } );
  oc.update_ref( T, [&]( OCFESLV::t_Coord const& cr ){ return TM_exact( cr.at(z), cr.at(t) ); } );
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

  oc.options.REDUCE.ORDER     = OCFESLV::Options::RED_FULL;
  oc.options.CLASSIFY.MODE         = OCFESLV::Options::CLASS_AUTO;
  oc.options.INTERFACE.TYPE   = OCFESLV::Options::IC_AUTO;
  oc.options.INTERFACE.IMPOSITION  = OCFESLV::Options::IC_STRONG;
  oc.options.INTERFACE.SAT_SIGMA0       = 10.0;
  oc.options.DISPLAY_LEVEL    = 0;

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

  double ce=0., qe=0., Te=0.;
  double const sg[4] = { 0.15, 0.35, 0.65, 0.85 };
  for( double zs : sg ) for( double ts : sg ){
    OCFESLV::t_Coord pt; pt[z]=zs; pt[t]=ts;
    ce = std::max( ce, std::fabs( oc.eval_colloc<double>( c, pt, xv.data(), nullptr, nullptr ) - cM_exact(zs,ts) ) );
    qe = std::max( qe, std::fabs( oc.eval_colloc<double>( q, pt, xv.data(), nullptr, nullptr ) - qM_exact(zs,ts) ) );
    Te = std::max( Te, std::fabs( oc.eval_colloc<double>( T, pt, xv.data(), nullptr, nullptr ) - TM_exact(zs,ts) ) );
  }
  std::cout << std::scientific << std::setprecision(4)
            << "  [C] err c=" << ce << " q=" << qe << " T=" << Te << "\n";

  bool ok = ( ce<1e-7 && qe<1e-7 && Te<1e-7 );
  std::cout << "  PSA5 " << cname << ": " << (ok?"PASS":"FAIL") << "\n";
  return ok;
}

int main()
{
  std::cout << "================================================================\n";
  std::cout << "  PSA rung 5: non-isothermal (energy balance + van't Hoff), MMS\n";
  std::cout << "  IC_STRONG, SIGMA0=10.  Coupling beta=" << beta_vh
            << " dH=" << dH << " hw=" << hw << "\n";
  std::cout << "================================================================\n";
  bool ok = true;
  ok &= run_psa5( FFDom::CGL, "CGL" );
  ok &= run_psa5( FFDom::LGL, "LGL" );
  std::cout << "\n  Overall: " << ( ok ? "ALL PASS" : "SOME FAILED" ) << "\n";
  return ok ? 0 : 1;
}
