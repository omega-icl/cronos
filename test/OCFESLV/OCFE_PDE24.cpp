// OCFE_PDE24_solve2.cpp   ---  PSA case study, RUNG 1
// ===========================================================================
// Isothermal, single-component, CONSTANT-velocity packed bed, single step.
// The PSA "skeleton" with momentum switched off: gas-phase mass balance
// (axial dispersion + convection + adsorption sink) coupled to a Linear
// Driving Force (LDF) loading equation with a Langmuir equilibrium.
//
//   gas:  dc/dt + u dc/dz - D d2c/dz2 + F dq/dt = s_c(z,t)
//   LDF:  dq/dt - k ( q*(c) - q )               = s_q(z,t),   q* = qs b c/(1+b c)
//
//   with u, D, F, k, qs, b CONSTANT (rung 1: no momentum coupling, u fixed).
//   F = (1-eps)/eps * rho_p  (lumped phase ratio).
//
// Boundary / initial:
//   inlet  z=0 : Danckwerts flux   u c - D dc/dz = g_in(t)      (BOUNDARY)
//   outlet z=L : zero gradient     dc/dz = 0                    (BOUNDARY)
//   initial t=0: c = c0(z), q = q0(z)                           (INITIAL)
//   q carries NO spatial BC (LDF is local in z): its dt-ODE holds at every z.
//
// VALIDATION: method of manufactured solutions.  Pick low-degree polynomial
//   c(z,t), q(z,t) (exactly representable on the collocation grid), derive the
//   source terms s_c, s_q and the Danckwerts RHS g_in analytically, and write
//   them as DAG expressions in (z,t).  Because the manufactured fields lie in
//   the collocation space, the discrete solution equals them exactly and the
//   solver recovers them to round-off -- even though the Langmuir source is
//   rational in c (it is evaluated exactly at every node).
//
//   c(z,t) = 1 + Az (z - z^2/2) + Ct t   ->  dc/dz = Az(1-z),  dc/dz(L=1)=0
//   q(z,t) = Q0 + Qz z + Qt t
//
// MODELLING RULE honoured: every flux term is constant * bare state-partial;
//   no conservative OpP(product) form (which the residual-Jacobian FAD rejects).
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

// ---- physical / numerical parameters --------------------------------------
static double const u_vel = 1.0;   // interstitial velocity
static double const D_ax  = 0.1;   // axial dispersion
static double const F_ph  = 0.5;   // phase ratio (1-eps)/eps * rho_p
static double const k_ldf = 2.0;   // LDF mass-transfer coefficient
static double const qs_L  = 1.0;   // Langmuir saturation loading
static double const b_L   = 1.0;   // Langmuir affinity

// ---- manufactured solution ------------------------------------------------
static double const Az = 0.5, Ct = 0.2;            // c(z,t)
static double const Q0 = 0.4, Qz = 0.1, Qt = 0.3;  // q(z,t)

static inline double cM_exact( double z, double t ){ return 1.0 + Az*( z - 0.5*z*z ) + Ct*t; }
static inline double qM_exact( double z, double t ){ return Q0 + Qz*z + Qt*t; }

static bool run_psa1( FFDom::TYPE coltype, std::string const& name )
{
  std::cout << "\n================================================================\n";
  std::cout << "  PSA rung 1  isothermal single-component breakthrough  " << name << "\n";
  std::cout << "================================================================\n";

  size_t const n_el_t = 2, n_nd_t = 4, n_el_z = 2, n_nd_z = 4;

  FFGraph DAG;
  FFVar t = DAG.add_var( "t" );
  FFVar z = DAG.add_var( "z" );
  FFVar c = DAG.add_var( "c(t,z)" );
  FFVar q = DAG.add_var( "q(t,z)" );

  FFPartial OpP;

  // Manufactured fields and the analytic sources, all as DAG expressions in (z,t).
  FFVar cMan  = 1.0 + Az*( z - 0.5*z*z ) + Ct*t;
  FFVar qMan  = Q0 + Qz*z + Qt*t;
  FFVar qstarM = qs_L*b_L*cMan / ( 1.0 + b_L*cMan );

  // s_c = dc/dt + u dc/dz - D d2c/dz2 + F dq/dt   (dc/dz=Az(1-z), d2c/dz2=-Az)
  FFVar s_c   = Ct + u_vel*Az*( 1.0 - z ) + D_ax*Az + F_ph*Qt;
  // s_q = dq/dt - k ( q*(cMan) - qMan )
  FFVar s_q   = Qt - k_ldf*( qstarM - qMan );
  // Danckwerts inlet RHS: u c(0,t) - D dc/dz(0) = u(1+Ct t) - D Az
  FFVar g_in  = u_vel*( 1.0 + Ct*t ) - D_ax*Az;

  // Residuals (every flux term is constant * bare state-partial).
  FFVar PDE_c = OpP( c, t ) + u_vel*OpP( c, z ) - D_ax*OpP( OpP( c, z ), z )
              + F_ph*OpP( q, t ) - s_c;
  FFVar LDF_q = OpP( q, t ) - k_ldf*( qs_L*b_L*c/( 1.0 + b_L*c ) - q ) - s_q;
  FFVar IC_c  = c - ( 1.0 + Az*( z - 0.5*z*z ) );    // c(z,0) = cM(z,0)
  FFVar IC_q  = q - ( Q0 + Qz*z );                   // q(z,0) = qM(z,0)
  FFVar BC_L  = u_vel*c - D_ax*OpP( c, z ) - g_in;   // Danckwerts inlet @ z=0
  FFVar BC_U  = OpP( c, z );                         // zero-gradient outlet @ z=1

  OCFESLV oc( &DAG );
  oc.add_domain( t, FFDom( 0., 1., n_el_t, coltype, n_nd_t ) );
  oc.add_domain( z, FFDom( 0., 1., n_el_z, coltype, n_nd_z ) );
  oc.add_state ( c, {t,z} );
  oc.add_state ( q, {t,z} );
  oc.update_ref( c, [&]( OCFESLV::t_Coord const& cr ){ return cM_exact( cr.at(z), cr.at(t) ); } );
  oc.update_ref( q, [&]( OCFESLV::t_Coord const& cr ){ return qM_exact( cr.at(z), cr.at(t) ); } );
  oc.set_evolution_domain( t );

  int const T_NO_LB = FFDom::ALL - FFDom::LB;
  int const Z_INT   = FFDom::ALL - FFDom::LB - FFDom::UB;

  oc.add_equation( PDE_c, {t,z}, {T_NO_LB,   Z_INT},      OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( LDF_q, {t,z}, {T_NO_LB,   FFDom::ALL}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( IC_c,  {t,z}, {FFDom::LB, FFDom::ALL}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INITIAL,  0 ) );
  oc.add_equation( IC_q,  {t,z}, {FFDom::LB, FFDom::ALL}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INITIAL,  0 ) );
  oc.add_equation( BC_L,  {t,z}, {T_NO_LB,   FFDom::LB},  OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );
  oc.add_equation( BC_U,  {t,z}, {T_NO_LB,   FFDom::UB},  OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );

  oc.options.REDUCE.ORDER     = OCFESLV::Options::RED_FULL;
  oc.options.CLASSIFY.MODE         = OCFESLV::Options::CLASS_AUTO;
  oc.options.INTERFACE.TYPE   = OCFESLV::Options::IC_AUTO;
  oc.options.INTERFACE.IMPOSITION  = OCFESLV::Options::IC_WEAK;
  oc.options.INTERFACE.SAT_SIGMA0       = 1.0;
  oc.options.DISPLAY_LEVEL    = 1;

  if( !oc.setup() ){ std::cerr << "ERROR: setup failed (" << name << ")\n"; return false; }
  std::cout << oc;

  size_t const nVar = oc.n_colloc_sta(), nEqn = oc.n_colloc_eqn();
  std::cout << "\nPDE type: " << OCFESLV::pde_type_name( oc.pde_type().type )
            << "  nVar=" << nVar << " nEqn=" << nEqn
            << " square=" << ( nVar==nEqn ? "yes":"NO" ) << "\n";
  bool ok = ( nVar == nEqn );

  std::vector<double> varInit, inpInit;
  if( !oc.init( varInit, inpInit, nullptr ) || varInit.size()!=nVar ){
    std::cerr << "ERROR: init failed (" << name << ")\n"; return false;
  }

  // Residual at the manufactured solution should be ~0 (it is the exact discrete root).
  std::vector<double> res( nEqn, 0. );
  if( !oc.eval( res.data(), nullptr, varInit.data(), nullptr, nullptr ) ){
    std::cerr << "ERROR: eval (" << name << ")\n"; return false;
  }
  double resmax = 0.; for( double v : res ) resmax = std::max( resmax, std::fabs(v) );
  std::cout << "[A] residual at manufactured exact (init): max|r|="
            << std::scientific << std::setprecision(4) << resmax
            << "  " << ( resmax<1e-9 ? "PASS" : "(init not exact)" ) << "\n";

  // Solve from a perturbed start.
  std::vector<double> xv = varInit;
  for( size_t i=0;i<nVar;++i ) xv[i] += 0.05*std::sin( 0.7*double(i) + 0.2 );
  OCFESLV::SolveReport rep = oc.solve( xv.data(), nullptr, nullptr );
  std::cout << "[B] solve(): converged=" << (rep.converged?"yes":"no")
            << " iters=" << rep.iterations
            << " init|r|=" << std::scientific << std::setprecision(4) << rep.initial_residual
            << " final|r|=" << rep.final_residual << "\n";

  // Recovery of both fields at sample points.
  double cerr_max = 0., qerr_max = 0., zc=0., tc=0., zq=0., tq=0.;
  if( rep.converged ){
    double const sg[4] = { 0.15, 0.35, 0.65, 0.85 };
    for( double zs : sg ) for( double ts : sg ){
      OCFESLV::t_Coord pt; pt[z]=zs; pt[t]=ts;
      double const ce = oc.eval_colloc<double>( c, pt, xv.data(), nullptr, nullptr );
      double const qe = oc.eval_colloc<double>( q, pt, xv.data(), nullptr, nullptr );
      double const dce = std::fabs( ce - cM_exact(zs,ts) );
      double const dqe = std::fabs( qe - qM_exact(zs,ts) );
      if( dce>cerr_max ){ cerr_max=dce; zc=zs; tc=ts; }
      if( dqe>qerr_max ){ qerr_max=dqe; zq=zs; tq=ts; }
    }
  }
  bool const c_ok = ( cerr_max < 1e-7 ), q_ok = ( qerr_max < 1e-7 );
  std::cout << std::defaultfloat
            << "[C] eval_colloc(c) vs manufactured: max err="
            << std::scientific << std::setprecision(4) << cerr_max
            << " at (z=" << std::fixed << std::setprecision(2) << zc << ",t=" << tc << ")  "
            << ( c_ok ? "PASS" : "FAIL" ) << "\n";
  std::cout << std::defaultfloat
            << "[C] eval_colloc(q) vs manufactured: max err="
            << std::scientific << std::setprecision(4) << qerr_max
            << " at (z=" << std::fixed << std::setprecision(2) << zq << ",t=" << tq << ")  "
            << ( q_ok ? "PASS" : "FAIL" ) << "\n";

  bool const pass = ok && rep.converged && c_ok && q_ok;
  std::cout << std::defaultfloat << "PSA1 " << name << ": " << ( pass ? "PASS" : "FAIL" );
  if( pass ) std::cout << "  -> isothermal LDF/Langmuir breakthrough skeleton recovered.";
  std::cout << "\n";
  return pass;
}

int main()
{
  std::cout << "================================================================\n";
  std::cout << "  PSA rung 1: isothermal, single-component, constant velocity\n";
  std::cout << "  gas balance + LDF + Langmuir, Danckwerts inlet (MMS)\n";
  std::cout << "================================================================\n";

  bool all = true;
  all &= run_psa1( FFDom::CGL, "CGL" );
  all &= run_psa1( FFDom::LGL, "LGL" );

  std::cout << "\n  Overall: " << ( all ? "ALL PASS" : "SOME FAILED" ) << "\n";
  return all ? 0 : 1;
}
