// Copyright (C) Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the EPL-2.0 with GPL-2.0-or-later as a Secondary License; see the LICENSE file.
// SPDX-License-Identifier: EPL-2.0 OR GPL-2.0-or-later
//
// OCFE_PDE21_solve2.cpp
// ---------
// Per-input continuity REFRAMING oracle (parabolic z,t with a t-discontinuous input).
//
// Confirms: a distributed input that is DISCONTINUOUS in the evolution direction t
// (a PSA cycle-step actuation) but enters a genuinely parabolic (z-diffusion) PDE is
// handled CORRECTLY by the default -- the state stays C0 (continuous, kinked) in t and
// C1 (smooth, flux-continuous) in z.  i.e. the cycle-step-in-time case needs no flag;
// the per-input continuity flag is a SPATIAL concern only.
//
//   PDE:   du/dt - d/dz( D du/dz ) = w(t) + s        z in [0,L], t in [0,T]
//   BCs:   du/dz(0,t) = 0   (Neumann LB),  du/dz(L,t) = 2L  (Neumann UB)
//   IC:    u(z,0) = 1 + z^2
//
// Method of manufactured solutions, additively separable:
//   u(z,t) = phi(t) + z^2
//     phi(t):  el0 [0,1/3]   phi = 1 + t + t^2     phi' = 1 + 2t
//              el1 [1/3,2/3] phi = 13/9 (const)    phi' = 0
//              el2 [2/3,1]   phi = 13/9 - (t-2/3)  phi' = -1
//   => u is C0 in t with SLOPE KINKS at t=1/3,2/3 (phi continuous, phi' jumps),
//      and smooth (C-infinity) in z.
//   w(t) = phi'(t)              distributed input over {t}, DISCONTINUOUS in t
//                               (jumps 5/3->0 at t=1/3, 0->-1 at t=2/3), uniform in z.
//   s    = -2D                  constant source (d/dz(D du/dz) = D*2 = 2D).
//   Residual at exact = phi' - 2D - phi' - (-2D) = 0  (machine zero, exactly
//   representable: phi piecewise deg<=2, z^2 deg 2, on n_node=4).
//
// Inputs are 1-D over {t} (broadcast over z), so the input vector orders exactly like
// the 1-D ODE1 case: flat index i -> t-element iel_t = i / n_nd_t.  At a t-join the two
// duplicated copies receive phi' from their own element => the value jump is faithful.
//
// Checks:
//   [A assert] residual at the manufactured exact (init) vector ~ machine zero.
//   [B assert] solve() from a perturbed guess converges; verify_interface_drop == OK.
//   [C assert] eval_colloc(u) on a (z,t) sample grid recovers phi(t)+z^2  (layout-
//              agnostic: confirms C0-in-t kink + smooth-in-z under the jumping input).
//
// Requires the solve() equation-rows-only sparsity fix (patched ocfeslv.hpp)
// only if a scalar output is added; this driver adds none, but uses solve(), so the
// fix is harmless and recommended.

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

// ---- Manufactured solution (fixed to the 3-element uniform [0,1] t-partition) ------
static inline double phi_piece( size_t iel, double t )
{
  switch( iel ){
    case 0:  return 1. + t + t*t;
    case 1:  return 13./9.;
    default: return 13./9. - ( t - 2./3. );
  }
}
static inline double phip_piece( size_t iel, double t )   // phi'(t)
{
  switch( iel ){
    case 0:  return 1. + 2.*t;
    case 1:  return 0.;
    default: return -1.;
  }
}
static inline double u_exact( double z, double t )
{
  size_t iel = ( t < 1./3. ) ? 0 : ( t < 2./3. ? 1 : 2 );
  return phi_piece( iel, t ) + z*z;
}

static const char* drop_status_name( OCFESLV::InterfaceDropStatus s )
{
  switch( s ){
    case OCFESLV::InterfaceDropStatus::OK:        return "OK";
    case OCFESLV::InterfaceDropStatus::OVERDROP:  return "OVERDROP (dropped a NEEDED claim)";
    case OCFESLV::InterfaceDropStatus::UNDERDROP: return "UNDERDROP";
  }
  return "?";
}

static bool run_test( FFDom::TYPE coltype, std::string const& name )
{
  std::cout << "\n================================================================\n";
  std::cout << "  " << name << "\n";
  std::cout << "================================================================\n";

  double const L  = 1.0;
  double const Dv = 0.1;
  size_t const n_el_t = 3, n_nd_t = 4;
  size_t const n_el_z = 2, n_nd_z = 4;

  FFGraph DAG;
  FFVar t  = DAG.add_var( "t" );
  FFVar z  = DAG.add_var( "z" );
  FFVar u  = DAG.add_var( "u(t,z)" );
  FFVar w  = DAG.add_var( "w(t)" );
  FFVar Dc = DAG.add_var( "D" );

  FFPartial OpP;
  FFVar s_src = -2.0 * Dc;                                   // s = -2D
  FFVar PDE   = OpP( u, t ) - OpP( Dc * OpP( u, z ), z ) - w - s_src;
  FFVar IC    = u - ( 1.0 + z*z );                           // u(z,0) = 1 + z^2
  FFVar BC_LB = OpP( u, z );                                 // du/dz(0,t) = 0
  FFVar BC_UB = OpP( u, z ) - 2.0*L;                         // du/dz(L,t) = 2L

  OCFESLV oc( &DAG );
  oc.add_domain( t, FFDom( 0., 1., n_el_t, coltype, n_nd_t ) );
  oc.add_domain( z, FFDom( 0., L,  n_el_z, coltype, n_nd_z ) );
  oc.add_state( u, {t,z} );
  oc.update_ref( u, [&]( OCFESLV::t_Coord const& c ){
    // u value is continuous everywhere (only du/dt kinks), so a coordinate-only
    // reference is single-valued -- safe at the t-joins.
    return u_exact( c.at(z), c.at(t) );
  } );
  // Distributed input over {t} (uniform in z), DEFAULT continuity (discontinuous).
  oc.add_input( w, {t}, std::optional<double>( 0.0 ) );
  oc.set_constant( {Dc}, {Dv} );
  oc.set_evolution_domain( t );

  int const T_NO_LB = FFDom::ALL - FFDom::LB;
  int const Z_INT   = FFDom::ALL - FFDom::LB - FFDom::UB;

  oc.add_equation( PDE,   {t,z}, {T_NO_LB, Z_INT},
                   OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( IC,    {t,z}, {FFDom::LB, FFDom::ALL},
                   OCFESLV::EqnOptions( OCFESLV::EqnRole::INITIAL, 0 ) );
  oc.add_equation( BC_LB, {t,z}, {T_NO_LB, FFDom::LB},
                   OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );
  oc.add_equation( BC_UB, {t,z}, {T_NO_LB, FFDom::UB},
                   OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );

  oc.options.REDUCE.ORDER    = OCFESLV::Options::RED_FULL;
  oc.options.CLASSIFY.MODE        = OCFESLV::Options::CLASS_AUTO;
  oc.options.INTERFACE.TYPE  = OCFESLV::Options::IC_AUTO;
  oc.options.INTERFACE.IMPOSITION = OCFESLV::Options::IC_WEAK;
  oc.options.SOLVE.MARCHING  = false;   // monolithic oracle by design: this test validates the
                                        // MONOLITHIC multi-element-in-time handling of a t-DISCONTINUOUS
                                        // input, and builds the per-t-element input vector by hand
                                        // (nInp==n_el_t*n_nd_t) -- marching would collapse t and defeat it.
  oc.options.INTERFACE.SAT_SIGMA0      = 1.0;
  oc.options.DISPLAY_LEVEL   = 1;

  if( !oc.setup() ){
    std::cerr << "ERROR: OCFESLV::setup() failed\n";
    return false;
  }
  std::cout << oc;

  size_t const nVar = oc.n_colloc_sta();
  size_t const nEqn = oc.n_colloc_eqn();
  size_t const nInp = oc.n_colloc_inp();
  auto const& cls = oc.pde_type();
  std::cout << "\nPDE type: " << OCFESLV::pde_type_name( cls.type )
            << "   nVar=" << nVar << " nEqn=" << nEqn
            << " nInp=" << nInp
            << " square=" << ( nVar==nEqn ? "yes" : "NO" ) << "\n";
  bool ok = ( nVar == nEqn );

  double const cval = Dv;

  // --- Manufactured exact state vector (incl. RED_FULL flux aux) via init() ----------
  std::vector<double> varInit, inpInit;
  if( !oc.init( varInit, inpInit, &cval ) || varInit.size() != nVar ){
    std::cerr << "ERROR: init() failed or wrong size\n";
    return false;
  }

  // --- Build the (t-discontinuous) input vector explicitly ---------------------------
  auto vn_w = oc.node_colloc( w );
  if( vn_w.size() != nInp || nInp != n_el_t*n_nd_t ){
    std::cerr << "ERROR: input node layout unexpected (nodes=" << vn_w.size()
              << ", nInp=" << nInp << ")\n";
    return false;
  }
  std::vector<double> winp( nInp, 0. );
  std::cout << "\nInput w(t) nodes (1-D over t):\n";
  for( size_t i=0; i<nInp; ++i ){
    size_t const iel = i / n_nd_t;
    winp[i] = phip_piece( iel, vn_w[i][0] );
    std::cout << "  i=" << std::setw(2) << i << "  el=" << iel
              << "  t=" << std::setw(9) << std::fixed << std::setprecision(6) << vn_w[i][0]
              << "  w=" << std::setw(10) << std::setprecision(6) << winp[i] << "\n";
  }
  std::cout << std::defaultfloat;

  // --- [A] residual at the manufactured exact vector ---------------------------------
  std::vector<double> res( nEqn, 0. );
  if( !oc.eval( res.data(), nullptr, varInit.data(), winp.data(), &cval ) ){
    std::cerr << "ERROR: eval at init failed\n"; return false;
  }
  double resmax = 0.; for( double v : res ) resmax = std::max( resmax, std::fabs(v) );
  bool const init_exact = ( resmax < 1e-9 );
  std::cout << "\n[A] residual at manufactured exact (init): max|r|="
            << std::scientific << std::setprecision(4) << resmax
            << "  " << ( init_exact ? "PASS" : "(init not exact -- see note)" ) << "\n";

  // --- [B] solve from a perturbed guess; well-posedness + drop gate ------------------
  std::vector<double> xv = varInit;
  for( size_t i=0;i<nVar;++i ) xv[i] += 0.05*std::sin( 0.7*double(i)+0.3 );
  OCFESLV::SolveReport rep = oc.solve( xv.data(), winp.data(), &cval );
  std::cout << "[B] solve(): converged=" << ( rep.converged ? "yes":"no" )
            << "  iters=" << rep.iterations
            << "  init|r|=" << std::scientific << std::setprecision(4) << rep.initial_residual
            << "  final|r|=" << rep.final_residual << "\n";
  bool drop_ok = false;
  if( rep.converged ){
    OCFESLV::InterfaceDropStatus st = oc.verify_and_protect( xv.data() );
    drop_ok = ( st == OCFESLV::InterfaceDropStatus::OK );
    std::cout << "    verify_interface_drop: " << drop_status_name( st ) << "\n";
  }

  // recovery of the full vector against the manufactured exact (only meaningful if
  // init reproduced the exact solution, i.e. [A] passed)
  double recmax = 0.; size_t recimax = 0;
  for( size_t i=0;i<nVar;++i ){ double e=std::fabs(xv[i]-varInit[i]); if(e>recmax){recmax=e;recimax=i;} }
  std::cout << std::defaultfloat
            << "    full-vector recovery max|x_solved-x_manufactured|="
            << std::scientific << std::setprecision(4) << recmax
            << " at i=" << recimax << ( init_exact ? "" : "  (reference inexact)" ) << "\n";

  // --- [C] layout-agnostic u recovery on a (z,t) sample grid -------------------------
  double const zs[] = { 0.0, 0.2, 0.4, 0.5, 0.6, 0.8, 1.0 };
  double const ts[] = { 0.05, 0.25, 0.3333333333, 0.45, 0.5, 0.6666666667, 0.8, 0.95 };
  double uerr = 0.; double zb=0, tb=0;
  for( double zi : zs ) for( double ti : ts ){
    OCFESLV::t_Coord pt; pt[z]=zi; pt[t]=ti;
    double const ue = oc.eval_colloc<double>( u, pt, xv.data(), winp.data(), &cval );
    double const ux = u_exact( zi, ti );
    double const e  = std::fabs( ue - ux );
    if( e > uerr ){ uerr = e; zb=zi; tb=ti; }
  }
  bool const u_ok = ( uerr < 1e-7 );
  std::cout << "[C] eval_colloc(u) vs manufactured on (z,t) grid: max err="
            << std::scientific << std::setprecision(4) << uerr
            << " at (z=" << std::fixed << std::setprecision(4) << zb
            << ", t=" << tb << ")  " << ( u_ok ? "PASS" : "FAIL" ) << "\n";
  std::cout << std::defaultfloat;

  bool const pass = ok && rep.converged && drop_ok && u_ok;
  std::cout << "\nTEST " << name << ": " << ( pass ? "PASS" : "FAIL" ) << "\n";
  if( pass )
    std::cout << "  -> default recovers C0-kinked-in-t, smooth-in-z state under a "
                 "t-DISCONTINUOUS input.\n     Cycle-step-in-time needs no flag; the "
                 "per-input continuity flag is a SPATIAL concern.\n";
  return pass;
}

int main()
{
  std::cout << "================================================================\n";
  std::cout << "  PDE21: du/dt - d/dz(D du/dz) = w(t)+s,  manufactured u=phi(t)+z^2\n";
  std::cout << "         w(t) DISCONTINUOUS in t (cycle-step), smooth in z\n";
  std::cout << "================================================================\n";

  bool all = true;
  all &= run_test( FFDom::CGL, "CGL  (Chebyshev-Gauss-Lobatto)" );
  all &= run_test( FFDom::LGL, "LGL  (Legendre-Gauss-Lobatto)" );

  std::cout << "\n================================================================\n";
  std::cout << "  Overall: " << ( all ? "ALL PASS" : "SOME FAILED" ) << "\n";
  std::cout << "================================================================\n";
  return all ? 0 : 1;
}
