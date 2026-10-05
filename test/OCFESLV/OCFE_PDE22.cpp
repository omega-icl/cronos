// Copyright (C) Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// OCFE_PDE22_solve2.cpp
// ---------
// Per-input continuity SPATIAL oracles:
//   A) coefficient discontinuity in z (differential state, flux continuity)
//   B) algebraic state fed by a z-discontinuous input (must stay discontinuous)
//
// The discontinuous quantity is supplied as an input carrying piecewise-constant data
// (it jumps at the z-join).  Each section runs twice: once with the input on the state
// z-grid (matching) and once with an n_node=1 LG input (one Gauss node per z-element,
// the genuine piecewise-constant coefficient).  The coarse case requires the OCVar
// operator*= finer-grid union patch: when the coarse coefficient multiplies a state
// gradient that RED_FULL folds into a flux aux (D*dz(u)), the product must lift the
// coarse operand UP to the state grid rather than collapsing onto the coarse grid.
// Input vectors order block-major over z: flat index i -> z-element iel = i / inodes.
//
// ---- A: coefficient discontinuity ------------------------------------------------
//   Steady 1-D:  -d/dz( D(z) du/dz ) = 0,   z in [0,1],  D = 1 (z<1/2), 2 (z>1/2)
//   BCs:         u(0)=0,  u(1)=0.75
//   Manufactured (flux-continuous, value-continuous, slope-kinked):
//     u = z                 on [0,1/2]   (u'=1,   flux D u' = 1)
//     u = 1/2 + (z-1/2)/2   on [1/2,1]   (u'=1/2, flux D u' = 2*1/2 = 1)
//   Flux D u' = 1 is continuous; u' kinks 1 -> 1/2.  The default must impose
//   FLUX continuity (C0 of D*du/dz), not GRADIENT continuity (C0 of du/dz), to
//   recover this.  Whether RED_FULL folds the input D into the flux aux is the
//   discriminator (read the LINK row in the dump).
//
// ---- B: algebraic state under a z-discontinuous input ----------------------------
//   (t,z):  dc/dt - d/dz(D0 dc/dz) = s_c        c differential, smooth
//           a - c - w(z)           = 0          a ALGEBRAIC (no d_t, no d_z, no IC)
//   c = 1 + t + z^2 (smooth);  w = 0.3 (z<1/2), 0.7 (z>1/2);  a = c + w (jumps in z).
//   Per the propagation rule, a (algebraic, fed by a z-discontinuous input) is
//   discontinuous in z: the default must emit NO z-C0 claim for a, letting it jump,
//   while c stays continuous.

#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <string>
#include <map>
#include <algorithm>

#ifndef OCFE_OCFESLV_HEADER
#define OCFE_OCFESLV_HEADER "ocfeslv.hpp"
#endif
#include OCFE_OCFESLV_HEADER

using namespace mc;

static const char* drop_status_name( OCFESLV::InterfaceDropStatus s )
{
  switch( s ){
    case OCFESLV::InterfaceDropStatus::OK:        return "OK";
    case OCFESLV::InterfaceDropStatus::OVERDROP:  return "OVERDROP (dropped a NEEDED claim)";
    case OCFESLV::InterfaceDropStatus::UNDERDROP: return "UNDERDROP";
  }
  return "?";
}

// ================================================================================
//  Section A: coefficient discontinuity (differential state, flux continuity)
// ================================================================================
static inline double uA_exact( double z )
{ return ( z < 0.5 ) ? z : 0.5 + 0.5*( z - 0.5 ); }

static bool run_coef_disc( FFDom::TYPE coltype, std::string const& name, bool coarse )
{
  std::string const tag = name + ( coarse ? "  [coarse n_node=1 input]" : "  [matching-grid input]" );
  std::cout << "\n================================================================\n";
  std::cout << "  A  coefficient discontinuity  " << tag << "\n";
  std::cout << "================================================================\n";

  size_t const n_el_z = 2, n_nd_z = 4;
  size_t const inodes = coarse ? 1 : n_nd_z;   // input nodes per z-element

  FFGraph DAG;
  FFVar z    = DAG.add_var( "z" );
  FFVar u    = DAG.add_var( "u(z)" );
  FFVar D    = DAG.add_var( "D(z)" );
  FFVar Uout = DAG.add_var( "Uout" );

  FFPartial OpP;
  FFVar PDE  = OpP( D * OpP( u, z ), z );          // -d/dz(D du/dz)=0 (sign irrelevant, =0)
  FFVar BC0  = u;                                  // u(0)=0
  FFVar BC1  = u - Uout;                           // u(1)=0.75

  OCFESLV oc( &DAG );
  oc.add_domain( z, FFDom( 0., 1., n_el_z, coltype, n_nd_z ) );
  oc.add_state ( u, {z} );
  oc.update_ref( u, [&]( OCFESLV::t_Coord const& c ){ return uA_exact( c.at(z) ); } );
  if( coarse )
    oc.add_input ( D, {z}, FFDom::LG, 1, std::optional<double>( 1.5 ) ); // 1 LG node/element
  else
    oc.add_input ( D, {z}, std::optional<double>( 1.5 ) );               // matching state z-grid
  // D carries piecewise-CONSTANT data (jumps at the join); set below.
  oc.set_constant( {Uout}, {0.75} );
  oc.reset_evolution_domain();

  int const Z_INT = FFDom::ALL - FFDom::LB - FFDom::UB;
  oc.add_equation( PDE, {z}, {Z_INT},     OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( BC0, {z}, {FFDom::LB}, OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );
  oc.add_equation( BC1, {z}, {FFDom::UB}, OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );

  oc.options.REDUCE.ORDER    = OCFESLV::Options::RED_FULL;
  oc.options.CLASSIFY.MODE        = OCFESLV::Options::CLASS_AUTO;
  oc.options.INTERFACE.TYPE  = OCFESLV::Options::IC_AUTO;
  oc.options.INTERFACE.IMPOSITION = OCFESLV::Options::IC_WEAK;
  oc.options.INTERFACE.SAT_SIGMA0      = 1.0;
  oc.options.DISPLAY_LEVEL   = 1;

  if( !oc.setup() ){ std::cerr << "ERROR: setup failed (A)\n"; return false; }
  std::cout << oc;

  size_t const nVar = oc.n_colloc_sta(), nEqn = oc.n_colloc_eqn(), nInp = oc.n_colloc_inp();
  std::cout << "\nPDE type: " << OCFESLV::pde_type_name( oc.pde_type().type )
            << "  nVar=" << nVar << " nEqn=" << nEqn << " nInp=" << nInp
            << " square=" << ( nVar==nEqn ? "yes":"NO" ) << "\n";
  bool ok = ( nVar == nEqn );

  double const cval = 0.75;
  std::vector<double> varInit, inpInit;
  if( !oc.init( varInit, inpInit, &cval ) || varInit.size()!=nVar ){
    std::cerr << "ERROR: init failed (A)\n"; return false;
  }

  // Discontinuous coefficient input data: D=1 on el0, D=2 on el1.
  auto vn_D = oc.node_colloc( D );
  if( vn_D.size() != nInp || nInp != n_el_z*inodes ){
    std::cerr << "ERROR: D input layout unexpected (nodes=" << vn_D.size()
              << " nInp=" << nInp << ")\n"; return false;
  }
  std::vector<double> Dinp( nInp, 0. );
  std::cout << "\nInput D(z) nodes (" << inodes << " per z-element):\n";
  for( size_t i=0;i<nInp;++i ){
    size_t const iel = i / inodes;
    Dinp[i] = ( iel==0 ) ? 1.0 : 2.0;
    std::cout << "  i=" << std::setw(2) << i << "  el=" << iel
              << "  z=" << std::fixed << std::setprecision(6) << vn_D[i][0]
              << "  D=" << Dinp[i] << "\n";
  }
  std::cout << std::defaultfloat;

  std::vector<double> res( nEqn, 0. );
  if( !oc.eval( res.data(), nullptr, varInit.data(), Dinp.data(), &cval ) ){
    std::cerr << "ERROR: eval (A)\n"; return false;
  }
  double resmax = 0.; for( double v : res ) resmax = std::max( resmax, std::fabs(v) );
  std::cout << "[A] residual at manufactured exact (init): max|r|="
            << std::scientific << std::setprecision(4) << resmax
            << "  " << ( resmax<1e-9 ? "PASS" : "(init not exact)" ) << "\n";

  std::vector<double> xv = varInit;
  for( size_t i=0;i<nVar;++i ) xv[i] += 0.05*std::sin(0.7*double(i)+0.2);
  OCFESLV::SolveReport rep = oc.solve( xv.data(), Dinp.data(), &cval );
  std::cout << "[B] solve(): converged=" << (rep.converged?"yes":"no")
            << " iters=" << rep.iterations
            << " init|r|=" << std::scientific << std::setprecision(4) << rep.initial_residual
            << " final|r|=" << rep.final_residual << "\n";
  bool drop_ok = false;
  if( rep.converged ){
    auto st = oc.verify_and_protect( xv.data() );
    drop_ok = ( st == OCFESLV::InterfaceDropStatus::OK );
    std::cout << "    verify_interface_drop: " << drop_status_name( st ) << "\n";
  }

  // Layout-agnostic recovery + explicit slope-kink probe at the z-join.
  double const zs[] = { 0.0, 0.15, 0.3, 0.45, 0.49, 0.51, 0.6, 0.8, 1.0 };
  double uerr = 0., zb = 0.;
  for( double zi : zs ){
    OCFESLV::t_Coord pt; pt[z]=zi;
    double const ue = oc.eval_colloc<double>( u, pt, xv.data(), Dinp.data(), &cval );
    double const e  = std::fabs( ue - uA_exact(zi) );
    if( e>uerr ){ uerr=e; zb=zi; }
  }
  bool const u_ok = ( uerr < 1e-7 );
  std::cout << std::defaultfloat
            << "[C] eval_colloc(u) vs manufactured: max err="
            << std::scientific << std::setprecision(4) << uerr
            << " at z=" << std::fixed << std::setprecision(4) << zb
            << "  " << ( u_ok ? "PASS" : "FAIL" ) << "\n";

  // Slope on each side of the join (finite difference of the interpolant).
  auto uz = [&]( double zi ){ OCFESLV::t_Coord pt; pt[z]=zi;
    return oc.eval_colloc<double>( u, pt, xv.data(), Dinp.data(), &cval ); };
  double const h = 1e-4;
  double const slopeL = ( uz(0.5-h) - uz(0.5-3.*h) ) / (2.*h);
  double const slopeR = ( uz(0.5+3.*h) - uz(0.5+h) ) / (2.*h);
  std::cout << "    join slopes: left~" << std::fixed << std::setprecision(5) << slopeL
            << " (exact 1.0)   right~" << slopeR << " (exact 0.5)   flux L="
            << 1.0*slopeL << "  flux R=" << 2.0*slopeR << " (both exact 1.0)\n";
  std::cout << std::defaultfloat;

  bool const pass = ok && rep.converged && drop_ok && u_ok;
  std::cout << "TEST A " << tag << ": " << ( pass ? "PASS" : "FAIL" );
  if( pass ) std::cout << "  -> flux continuity recovered across a discontinuous coefficient.";
  else       std::cout << "  -> coefficient discontinuity NOT handled by default (Hook-2 demotion needed).";
  std::cout << "\n";
  return pass;
}

// ================================================================================
//  Section B: algebraic state fed by a z-discontinuous input
// ================================================================================
static inline double cB_exact( double z, double t ){ return 1. + t + z*z; }
static inline double wB_value( double z ){ return ( z < 0.5 ) ? 0.3 : 0.7; }
static inline double aB_exact( double z, double t ){ return cB_exact(z,t) + wB_value(z); }

static bool run_algebraic( FFDom::TYPE coltype, std::string const& name, bool coarse )
{
  std::string const tag = name + ( coarse ? "  [coarse n_node=1 input]" : "  [matching-grid input]" );
  std::cout << "\n================================================================\n";
  std::cout << "  B  algebraic state + z-discontinuous input  " << tag << "\n";
  std::cout << "================================================================\n";

  double const D0 = 0.1;
  size_t const n_el_t = 3, n_nd_t = 4, n_el_z = 2, n_nd_z = 4;
  size_t const inodes = coarse ? 1 : n_nd_z;   // input nodes per z-element

  FFGraph DAG;
  FFVar t   = DAG.add_var( "t" );
  FFVar z   = DAG.add_var( "z" );
  FFVar c   = DAG.add_var( "c(t,z)" );
  FFVar a   = DAG.add_var( "a(t,z)" );
  FFVar w   = DAG.add_var( "w(z)" );
  FFVar D0c = DAG.add_var( "D0" );

  FFPartial OpP;
  FFVar s_c   = 1.0 - 2.0*D0c;                          // dc/dt - d/dz(D0 dc/dz) = 1 - 2D0
  FFVar PDE_C = OpP( c, t ) - OpP( D0c*OpP( c, z ), z ) - s_c;
  FFVar ALG_A = a - c - w;                              // a = c + w (algebraic, pointwise)
  FFVar IC_C  = c - ( 1.0 + z*z );
  FFVar BC_L  = OpP( c, z );                            // dc/dz(0,t)=0
  FFVar BC_U  = OpP( c, z ) - 2.0;                      // dc/dz(L,t)=2

  OCFESLV oc( &DAG );
  oc.add_domain( t, FFDom( 0., 1., n_el_t, coltype, n_nd_t ) );
  oc.add_domain( z, FFDom( 0., 1., n_el_z, coltype, n_nd_z ) );
  oc.add_state ( c, {t,z} );
  oc.add_state ( a, {t,z} );                            // ALGEBRAIC: no d_t, no IC
  oc.update_ref( c, [&]( OCFESLV::t_Coord const& cr ){ return cB_exact( cr.at(z), cr.at(t) ); } );
  oc.update_ref( a, [&]( OCFESLV::t_Coord const& cr ){ return aB_exact( cr.at(z), cr.at(t) ); } );
  if( coarse )
    oc.add_input ( w, {z}, FFDom::LG, 1, std::optional<double>( 0.5 ) ); // 1 LG node/element
  else
    oc.add_input ( w, {z}, std::optional<double>( 0.5 ) );               // matching state z-grid
  oc.set_constant( {D0c}, {D0} );
  oc.set_evolution_domain( t );

  int const T_NO_LB = FFDom::ALL - FFDom::LB;
  int const Z_INT   = FFDom::ALL - FFDom::LB - FFDom::UB;

  oc.add_equation( PDE_C, {t,z}, {T_NO_LB,   Z_INT},     OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( ALG_A, {t,z}, {FFDom::ALL, FFDom::ALL},OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( IC_C,  {t,z}, {FFDom::LB,  FFDom::ALL},OCFESLV::EqnOptions( OCFESLV::EqnRole::INITIAL,  0 ) );
  oc.add_equation( BC_L,  {t,z}, {T_NO_LB,   FFDom::LB}, OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );
  oc.add_equation( BC_U,  {t,z}, {T_NO_LB,   FFDom::UB}, OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );

  oc.options.REDUCE.ORDER    = OCFESLV::Options::RED_FULL;
  oc.options.CLASSIFY.MODE        = OCFESLV::Options::CLASS_AUTO;
  oc.options.INTERFACE.TYPE  = OCFESLV::Options::IC_AUTO;
  oc.options.INTERFACE.IMPOSITION = OCFESLV::Options::IC_WEAK;
  oc.options.INTERFACE.SAT_SIGMA0      = 1.0;
  oc.options.DISPLAY_LEVEL   = 1;

  if( !oc.setup() ){ std::cerr << "ERROR: setup failed (B)\n"; return false; }
  std::cout << oc;

  size_t const nVar = oc.n_colloc_sta(), nEqn = oc.n_colloc_eqn(), nInp = oc.n_colloc_inp();
  std::cout << "\nPDE type: " << OCFESLV::pde_type_name( oc.pde_type().type )
            << "  nVar=" << nVar << " nEqn=" << nEqn << " nInp=" << nInp
            << " square=" << ( nVar==nEqn ? "yes":"NO" ) << "\n";
  std::cout << "States after reduction:";
  for( auto const& st : oc.states_colloc() ) std::cout << " " << st;
  std::cout << "\n";
  bool ok = ( nVar == nEqn );

  double const cval = D0;
  std::vector<double> varInit, inpInit;
  if( !oc.init( varInit, inpInit, &cval ) || varInit.size()!=nVar ){
    std::cerr << "ERROR: init failed (B)\n"; return false;
  }

  auto vn_w = oc.node_colloc( w );
  if( vn_w.size() != nInp || nInp != n_el_z*inodes ){
    std::cerr << "ERROR: w input layout unexpected (nodes=" << vn_w.size()
              << " nInp=" << nInp << ")\n"; return false;
  }
  std::vector<double> winp( nInp, 0. );
  std::cout << "\nInput w(z) nodes (" << inodes << " per z-element):\n";
  for( size_t i=0;i<nInp;++i ){
    size_t const iel = i / inodes;
    winp[i] = ( iel==0 ) ? 0.3 : 0.7;
    std::cout << "  i=" << std::setw(2) << i << "  el=" << iel
              << "  z=" << std::fixed << std::setprecision(6) << vn_w[i][0]
              << "  w=" << winp[i] << "\n";
  }
  std::cout << std::defaultfloat;

  std::vector<double> res( nEqn, 0. );
  if( !oc.eval( res.data(), nullptr, varInit.data(), winp.data(), &cval ) ){
    std::cerr << "ERROR: eval (B)\n"; return false;
  }
  double resmax = 0.; for( double v : res ) resmax = std::max( resmax, std::fabs(v) );
  std::cout << "[A] residual at manufactured exact (init): max|r|="
            << std::scientific << std::setprecision(4) << resmax
            << "  " << ( resmax<1e-9 ? "PASS" : "(init not exact at the a-join)" ) << "\n";

#if defined(CRONOS__WITH_SPQR)
  oc.options.SOLVE.FACTORIZATION = OCFESLV::Options::SOLVE_SPQR;
#endif

  std::vector<double> xv = varInit;
  for( size_t i=0;i<nVar;++i ) xv[i] += 0.05*std::sin(0.6*double(i)+0.1);
  OCFESLV::SolveReport rep = oc.solve( xv.data(), winp.data(), &cval );
  std::cout << "[B] solve(): converged=" << (rep.converged?"yes":"no")
            << " iters=" << rep.iterations
            << " init|r|=" << std::scientific << std::setprecision(4) << rep.initial_residual
            << " final|r|=" << rep.final_residual << "\n";
  bool drop_ok = false;
  if( rep.converged ){
    auto st = oc.verify_and_protect( xv.data() );
    drop_ok = ( st == OCFESLV::InterfaceDropStatus::OK );
    std::cout << "    verify_interface_drop: " << drop_status_name( st ) << "\n";
  }

  // c recovery (continuous) and a recovery (must jump in z) on a sample grid.
  double const zs[] = { 0.1, 0.3, 0.4, 0.6, 0.7, 0.9 };
  double const ts[] = { 0.1, 0.4, 0.7, 0.9 };
  double cerr=0., aerr=0.;
  for( double zi : zs ) for( double ti : ts ){
    OCFESLV::t_Coord pt; pt[z]=zi; pt[t]=ti;
    double const ce = oc.eval_colloc<double>( c, pt, xv.data(), winp.data(), &cval );
    double const ae = oc.eval_colloc<double>( a, pt, xv.data(), winp.data(), &cval );
    cerr = std::max( cerr, std::fabs( ce - cB_exact(zi,ti) ) );
    aerr = std::max( aerr, std::fabs( ae - aB_exact(zi,ti) ) );
  }
  bool const c_ok = ( cerr < 1e-7 ), a_ok = ( aerr < 1e-7 );
  std::cout << "[C] eval_colloc(c) continuous: max err=" << std::scientific << std::setprecision(4)
            << cerr << "  " << ( c_ok ? "PASS" : "FAIL" ) << "\n";
  std::cout << "[C] eval_colloc(a) jumps in z:  max err=" << aerr
            << "  " << ( a_ok ? "PASS" : "FAIL" ) << "\n";

  // Explicit jump probe at the z-join (t=0.5): a(0.4) uses w_L, a(0.6) uses w_R.
  OCFESLV::t_Coord pL; pL[z]=0.4; pL[t]=0.5;
  OCFESLV::t_Coord pR; pR[z]=0.6; pR[t]=0.5;
  double const aL = oc.eval_colloc<double>( a, pL, xv.data(), winp.data(), &cval );
  double const aR = oc.eval_colloc<double>( a, pR, xv.data(), winp.data(), &cval );
  std::cout << std::defaultfloat
            << "    a(z=0.4,t=0.5)=" << std::setprecision(6) << aL
            << " (exact " << aB_exact(0.4,0.5) << ")   a(z=0.6,t=0.5)=" << aR
            << " (exact " << aB_exact(0.6,0.5) << ")   jump~"
            << std::fabs(aR-aL) << " (exact 0.4)\n";

  bool const pass = ok && rep.converged && drop_ok && c_ok && a_ok;
  std::cout << "TEST B " << tag << ": " << ( pass ? "PASS" : "FAIL" );
  if( pass ) std::cout << "  -> algebraic a left discontinuous (no z-C0 claim); c continuous.";
  else       std::cout << "  -> algebraic discontinuity NOT handled by default.";
  std::cout << "\n";
  return pass;
}

int main()
{
  std::cout << "================================================================\n";
  std::cout << "  PDE22: spatial input-continuity oracles\n";
  std::cout << "    A) discontinuous coefficient D(z) (flux continuity)\n";
  std::cout << "    B) algebraic state a fed by z-discontinuous input w(z)\n";
  std::cout << "================================================================\n";

  bool all = true;
  // matching-grid inputs: patch must stay INERT (these already passed pre-patch)
  all &= run_coef_disc( FFDom::CGL, "CGL", false );
  all &= run_coef_disc( FFDom::LGL, "LGL", false );
  all &= run_algebraic( FFDom::CGL, "CGL", false );
  all &= run_algebraic( FFDom::LGL, "LGL", false );
  // coarse n_node=1 inputs: the case that threw pre-patch; the operator*= finer-grid
  // union must now lift the coarse coefficient up to the state grid.
  all &= run_coef_disc( FFDom::CGL, "CGL", true );
  all &= run_coef_disc( FFDom::LGL, "LGL", true );
  all &= run_algebraic( FFDom::CGL, "CGL", true );
  all &= run_algebraic( FFDom::LGL, "LGL", true );

  std::cout << "\n================================================================\n";
  std::cout << "  Overall: " << ( all ? "ALL PASS" : "SOME FAILED" ) << "\n";
  std::cout << "================================================================\n";
  return all ? 0 : 1;
}
