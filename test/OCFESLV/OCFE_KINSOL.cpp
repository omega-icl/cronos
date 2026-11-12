// OCFE_KINSOL.cpp  ---  KINSOL as OCFESLV's nonlinear solver (SOLVE.NLSOLVER = KINSOL; 10 Nov 2026, K2)
// ===========================================================================
// SOLVE.NLSOLVER = KINSOL solves the collocation system with SUNDIALS KINSOL through BASE_KINSOL, on solve()'s own
// residual and assembly -- so a monolithic solve, every marching window and every homotopy step alike -- factorising
// with SOLVE.FACTORIZATION (SUPERLU, KLU, UMFPACK, SPQR sparse; LAPACK dense), scaled by the equilibration, and
// accepting on RES_TOL as LM does (KIN.FTOL = 0 derives it).  Measured first (K0): the same roots as LM, KINSOL 15/15
// once FTOL follows RES_TOL and the rows are equilibrated.
//   (A) KINSOL = LM, every factorisation: the nonlinear-diffusion rod (monolithic and MARCHING), steady Bratu, the PSA
//       column (DAESLV_molsens's declaration) -- the solution and an output to round-off;
//       -- and the rod under IC_TRACE (the plan eliminates the multipliers before the solve: square in the states)
//       and IC_WEAK;
//   (B) the refusal, by name: SOLVE_KLU under LM;
//   (C) forward sensitivities after a KINSOL solve = LM's;
//   (E) a SINGULAR START (11 Nov 2026, from OCFE_ALGFIELD2): a^3 + a = a*^3 + a*, p a = p* a*, both fields started
//       at zero, so d(p a)/dp = a = 0 and the Jacobian's p block is zero at the first iterate -- KINSOL/SUPERLU takes
//       the QR's minimum-norm step there (SuperLU alone failed) and converges to a* and p*, as LM;
//   (F) a STRUCTURALLY rank-deficient system the caller opted past (OCFE_PDE20e's case D: y pinned twice, x free at
//       t0, FATAL.REDUCED_DOF = false): KINSOL takes the QR directly -- SuperLU, handed it, called BLAS with an invalid
//       dimension ('DGEMV parameter number 2') -- and returns without a BLAS error;
//   (D) homotopy under KINSOL (K3): Bratu continued in lambda from 0.5 to 3.4 (the turning point ~3.51), from zero --
//       every continuation step re-enters solve(), so takes the selected solver: KINSOL = LM on the final solution.
#include <iostream>
#include <iomanip>
#include <cmath>
#include <functional>
#include "ocfeslv.hpp"
#include "daeslv_idas.hpp"
using namespace mc;

static int nfail = 0;
static void check( std::string const& w, bool ok ){ std::cout << "  " << std::left << std::setw(96) << w << ( ok? "PASS": "FAIL" ) << "\n"; if( !ok ) ++nfail; }

static double const U_VEL = 1.0, D_ax = 0.1, F_ph = 0.5, k_ldf = 2.0;
static double const qs1 = 1.0, qs2 = 1.0, b02 = 1.0, beta1 = 2.0, beta2 = 1.0, T0_ref = 1.0;
static double const Cp_e = 1.0, G_cv = 1.0, lam = 0.1, dH1 = 1.0, dH2 = 0.5, hw = 0.5, Tw = 1.0;
static double const c0_1 = 0.5, c0_2 = 0.5, tau_in = 0.15, T_END = 5.0;

struct Fields { FFVar t, z, c1, c2, q1, q2, T, b01; };

// the collocation options, the same for both solvers -- OCFESLV holds them itself, DAESLV for its semi-discretisation
static void colloc_options( OCFESLV& oc )
{ oc.options.INTERFACE.IMPOSITION = OCFESLV::Options::IC_STRONG;  oc.options.INTERFACE.SAT_SIGMA0 = 10.0; }
static void colloc_options( DAESLV_IDAS& S )
{ S.options.INTERFACE.IMPOSITION = OCFEMOL::Options::IC_STRONG;  S.options.INTERFACE.SAT_SIGMA0 = 10.0; }

template <typename OC>
static void declare_psa( OC& oc, FFGraph& G, Fields& F, size_t NEZ, size_t NNZ, size_t NELT, size_t NODT )
{
  F.t = G.add_var("t"); F.z = G.add_var("z");
  F.c1 = G.add_var("c1"); F.c2 = G.add_var("c2"); F.q1 = G.add_var("q1"); F.q2 = G.add_var("q2"); F.T = G.add_var("T");
  F.b01 = G.add_var("b01");
  FFVar const &t = F.t, &z = F.z, &c1 = F.c1, &c2 = F.c2, &q1 = F.q1, &q2 = F.q2, &T = F.T, &b01 = F.b01;
  FFPartial OpP;  FFIntegral OpI;
  FFVar c1feed = c0_1*( 1.0 - exp( -t/tau_in ) ), c2feed = c0_2*( 1.0 - exp( -t/tau_in ) );
  FFVar b1 = b01*exp( beta1*( 1.0/T - 1.0/T0_ref ) ), b2 = b02*exp( beta2*( 1.0/T - 1.0/T0_ref ) );
  FFVar den = 1.0 + b1*c1 + b2*c2, q1star = qs1*b1*c1/den, q2star = qs2*b2*c2/den;
  FFVar CONT1 = OpP( c1, t ) + U_VEL*OpP( c1, z ) - D_ax*OpP( OpP( c1, z ), z ) + F_ph*OpP( q1, t );
  FFVar CONT2 = OpP( c2, t ) + U_VEL*OpP( c2, z ) - D_ax*OpP( OpP( c2, z ), z ) + F_ph*OpP( q2, t );
  FFVar LDF1  = OpP( q1, t ) - k_ldf*( q1star - q1 );
  FFVar LDF2  = OpP( q2, t ) - k_ldf*( q2star - q2 );
  FFVar ENE_T = Cp_e*OpP( T, t ) + G_cv*OpP( T, z ) - lam*OpP( OpP( T, z ), z )
              - F_ph*( dH1*OpP( q1, t ) + dH2*OpP( q2, t ) ) + hw*( T - Tw );
  oc.add_domain( t, FFDom( 0., T_END, NELT, FFDom::LGR, NODT ) );
  oc.add_domain( z, FFDom( 0., 1.0,   NEZ,  FFDom::LGL, NNZ  ) );
  oc.set_evolution_domain( t );
  oc.add_state( c1, {t,z} ); oc.update_ref( c1, 0.25 );  oc.add_state( c2, {t,z} ); oc.update_ref( c2, 0.25 );
  oc.add_state( q1, {t,z} ); oc.update_ref( q1, 0.25 );  oc.add_state( q2, {t,z} ); oc.update_ref( q2, 0.25 );
  oc.add_state( T,  {t,z} ); oc.update_ref( T, T0_ref );
  oc.add_input( b01, {} );
  int const TI = FFDom::ALL - FFDom::LB, ZI = FFDom::ALL - FFDom::LB - FFDom::UB;
  typename OC::EqnOptions io( OC::EqnRole::INTERIOR, 0 ), ib( OC::EqnRole::BOUNDARY, 0 ), ii( OC::EqnRole::INITIAL, 0 );
  oc.add_equation( CONT1, {t,z}, {TI, ZI}, io );  oc.add_equation( CONT2, {t,z}, {TI, ZI}, io );
  oc.add_equation( LDF1,  {t,z}, {TI, FFDom::ALL}, io );  oc.add_equation( LDF2, {t,z}, {TI, FFDom::ALL}, io );
  oc.add_equation( ENE_T, {t,z}, {TI, ZI}, io );
  oc.add_equation( c1, {t,z}, {FFDom::LB, FFDom::ALL}, ii );  oc.add_equation( c2, {t,z}, {FFDom::LB, FFDom::ALL}, ii );
  oc.add_equation( q1, {t,z}, {FFDom::LB, FFDom::ALL}, ii );  oc.add_equation( q2, {t,z}, {FFDom::LB, FFDom::ALL}, ii );
  oc.add_equation( T - T0_ref, {t,z}, {FFDom::LB, FFDom::ALL}, ii );
  oc.add_equation( U_VEL*c1 - D_ax*OpP( c1, z ) - U_VEL*c1feed, {t,z}, {TI, FFDom::LB}, ib );
  oc.add_equation( U_VEL*c2 - D_ax*OpP( c2, z ) - U_VEL*c2feed, {t,z}, {TI, FFDom::LB}, ib );
  oc.add_equation( G_cv*T - lam*OpP( T, z ) - G_cv*T0_ref,      {t,z}, {TI, FFDom::LB}, ib );
  oc.add_equation( OpP( c1, z ), {t,z}, {TI, FFDom::UB}, ib );
  oc.add_equation( OpP( c2, z ), {t,z}, {TI, FFDom::UB}, ib );
  oc.add_equation( OpP( T,  z ), {t,z}, {TI, FFDom::UB}, ib );
  oc.add_output( OpI( c1 + F_ph*q1, z ), {t}, {T_END} );          // Inv1
  oc.add_output( OpI( U_VEL*c1, t ),     {z}, {1.0}   );          // Eff1
  colloc_options( oc );
  oc.options.DISPLAY_LEVEL = 0;
}


static char const* FN[5] = { "SUPERLU", "KLU", "UMFPACK", "SPQR", "LAPACK" };
static OCFESLV::Options::SolveFactorization FZ[5] = { OCFESLV::Options::SOLVE_SUPERLU, OCFESLV::Options::SOLVE_KLU, OCFESLV::Options::SOLVE_UMFPACK,
                                                      OCFESLV::Options::SOLVE_SPQR, OCFESLV::Options::SOLVE_LAPACK };
// one model: LM, then KINSOL with each factorisation, from the same initialisation
static void compare( std::string const& name, std::function<void(OCFESLV&, FFGraph&, std::vector<double>&)> declare, bool march, int nfz = 5,
                     OCFESLV::Options::ImpositionType imp = OCFESLV::Options::IC_STRONG ){
  std::vector<double> ref;  double fref = NAN;
  for( int v = -1; v < nfz; ++v ){
    FFGraph G;  OCFESLV oc( &G );  std::vector<double> inpv;  declare( oc, G, inpv );
    oc.options.INTERFACE.IMPOSITION = imp;  oc.options.SOLVE.MARCHING = march;  oc.options.SOLVE.RES_TOL = 1e-9;  oc.options.DISPLAY_LEVEL = 0;
    if( v >= 0 ){ oc.options.SOLVE.NLSOLVER = OCFESLV::Options::KINSOL;  oc.options.SOLVE.FACTORIZATION = FZ[v]; }
    bool ok = oc.setup();
    std::vector<double> var, inp;
    if( ok ){ oc.init( var, inp );  for( size_t i = 0; i < inpv.size() && i < inp.size(); ++i ) inp[i] = inpv[i]; }
    auto rep = ok? oc.solve( var.data(), inp.empty()? nullptr: inp.data(), nullptr ): OCFESLV::SolveReport();
    ok = ok && rep.converged;
    double const f = ok && !oc.val_function().empty()? oc.val_function()[0]: NAN;
    if( v < 0 ){ ref = var;  fref = f;  check( "(A) " + name + ": LM converges", ok );  continue; }
    double d = NAN;  if( ok ){ d = 0.; for( size_t i = 0; i < var.size(); ++i ) d = std::max( d, std::fabs( var[i] - ref[i] )/std::max( 1., std::fabs( ref[i] ) ) ); }
    std::cout << "    KINSOL/" << std::left << std::setw(8) << FN[v] << std::right << rep.iterations << " it, max|r| " << std::scientific << std::setprecision(1)
              << rep.final_residual << ", solution vs LM " << d << ", output vs LM " << std::fabs( f - fref )/std::max( 1e-300, std::fabs( fref ) ) << "\n";
    check( "(A) " + name + ": KINSOL/" + FN[v] + " = LM, the solution to 1e-10 and the output to 1e-10", ok && d < 1e-10 && std::fabs( f - fref ) <= 1e-10*std::max( 1., std::fabs( fref ) ) );
  }
}

int main()
{
  FFPartial OpP;  FFIntegral OpI;
  auto rod = [&]( OCFESLV& oc, FFGraph& G, std::vector<double>& ){
    FFVar t = G.add_var("t"), z = G.add_var("z"), u = G.add_var("u");
    oc.add_domain( t, FFDom( 0., 0.2, 4, FFDom::LGR, 4 ) );  oc.add_domain( z, FFDom( 0., 1., 4, FFDom::LGL, 5 ) );
    oc.set_evolution_domain( t );  oc.add_state( u, {t,z} );  oc.update_ref( u, 0.5 );
    int const TI = FFDom::ALL - FFDom::LB, ZI = FFDom::ALL - FFDom::LB - FFDom::UB;
    OCFESLV::EqnOptions io( OCFESLV::EqnRole::INTERIOR, 0 ), ib( OCFESLV::EqnRole::BOUNDARY, 0 ), ii( OCFESLV::EqnRole::INITIAL, 0 );
    oc.add_equation( OpP(u,t) - OpP( ( 1. + u*u )*OpP(u,z), z ), {t,z}, {TI, ZI}, io );
    oc.add_equation( u, {t,z}, {TI, FFDom::LB}, ib );  oc.add_equation( u, {t,z}, {TI, FFDom::UB}, ib );
    oc.add_equation( u - sin( M_PI*z ), {t,z}, {FFDom::LB, FFDom::ALL}, ii );
    oc.add_output( OpI( OpI( u, z ), t ) ); };
  auto bratu = [&]( OCFESLV& oc, FFGraph& G, std::vector<double>& ){
    FFVar z = G.add_var("z"), u = G.add_var("u");
    oc.add_domain( z, FFDom( 0., 1., 4, FFDom::LGL, 6 ) );  oc.add_state( u, {z} );  oc.update_ref( u, 0.1 );
    OCFESLV::EqnOptions io( OCFESLV::EqnRole::INTERIOR, 0 ), ib( OCFESLV::EqnRole::BOUNDARY, 0 );
    oc.add_equation( -OpP( OpP(u,z), z ) - 3.*exp( u ), {z}, { FFDom::ALL - FFDom::LB - FFDom::UB }, io );
    oc.add_equation( u, {z}, {FFDom::LB}, ib );  oc.add_equation( u, {z}, {FFDom::UB}, ib );
    oc.add_output( OpI( u, z ) ); };
  Fields PF;
  auto psa = [&]( OCFESLV& oc, FFGraph& G, std::vector<double>& inp ){ declare_psa( oc, G, PF, 3, 5, 10, 3 );  inp = { 3.0 }; };

  std::cout << "(A) KINSOL = LM, every factorisation\n";
  compare( "the nonlinear-diffusion rod, monolithic", rod, false );
  compare( "the nonlinear-diffusion rod, MARCHING", rod, true );
  compare( "steady Bratu", bratu, false );
  compare( "the PSA column, monolithic", psa, false, 2 );     // SuperLU and KLU (3600 unknowns)
  compare( "the rod under IC_TRACE (the plan eliminates its multipliers)", rod, false, 2, OCFESLV::Options::IC_TRACE );
  compare( "the rod under IC_WEAK (the legacy SAT)", rod, false, 2, OCFESLV::Options::IC_WEAK );

  std::cout << "(B) the refusals, by name\n";
  { FFGraph G;  OCFESLV oc( &G );  std::vector<double> iv;  rod( oc, G, iv );
    oc.options.SOLVE.MARCHING = false;  oc.options.SOLVE.FACTORIZATION = OCFESLV::Options::SOLVE_KLU;  oc.options.DISPLAY_LEVEL = 0;
    bool const su = oc.setup();  std::vector<double> var, inp;  if( su ) oc.init( var, inp );
    auto rep = su? oc.solve( var.data(), inp.empty()? nullptr: inp.data(), nullptr ): OCFESLV::SolveReport();
    check( "(B) SOLVE_KLU under LEVENBERG_MARQUARDT: refused (not converged, named on stderr)", su && !rep.converged && rep.iterations == 0 ); }


  std::cout << "(C) forward sensitivities after a KINSOL solve = LM's (the rod with a source p e^-t)\n";
  { std::vector< std::vector<double> > Gr[2];  bool ok = true;
    for( int v = 0; v < 2; ++v ){
      FFGraph G;  OCFESLV oc( &G );
      FFVar t = G.add_var("t"), z = G.add_var("z"), u = G.add_var("u"), p = G.add_var("p");
      oc.add_domain( t, FFDom( 0., 0.5, 3, FFDom::LGR, 4 ) );  oc.add_domain( z, FFDom( 0., 1., 3, FFDom::LGL, 5 ) );
      oc.set_evolution_domain( t );  oc.add_state( u, {t,z} );  oc.update_ref( u, 1. );  oc.add_input( p, {} );
      int const TI = FFDom::ALL - FFDom::LB, ZI = FFDom::ALL - FFDom::LB - FFDom::UB;
      OCFESLV::EqnOptions io( OCFESLV::EqnRole::INTERIOR, 0 ), ib( OCFESLV::EqnRole::BOUNDARY, 0 ), ii( OCFESLV::EqnRole::INITIAL, 0 );
      oc.add_equation( OpP(u,t) - OpP( ( 1. + 0.5*u*u )*OpP(u,z), z ) - p*exp( -t ), {t,z}, {TI, ZI}, io );
      oc.add_equation( OpP(u,z), {t,z}, {TI, FFDom::LB}, ib );  oc.add_equation( u - 1., {t,z}, {TI, FFDom::UB}, ib );
      oc.add_equation( u - 1., {t,z}, {FFDom::LB, FFDom::ALL}, ii );
      oc.add_output( OpI( OpI( u*u, z ), t ) );
      oc.options.SOLVE.MARCHING = false;  oc.options.SOLVE.RES_TOL = 1e-10;  oc.options.DISPLAY_LEVEL = 0;
      if( v ) oc.options.SOLVE.NLSOLVER = OCFESLV::Options::KINSOL;
      if( !oc.setup() || oc.solve_fsens( std::vector<FFModel::InputVal>{ { p, { 0.7 } } } ) != OCFESLV::NORMAL ){ ok = false;  break; }
      Gr[v] = oc.val_function_gradient();
    }
    double d = NAN;  if( ok && !Gr[0].empty() && !Gr[1].empty() ) d = std::fabs( Gr[1][0][0] - Gr[0][0][0] )/std::max( 1e-300, std::fabs( Gr[0][0][0] ) );
    std::cout << "    d(output)/dp: LM " << std::scientific << std::setprecision(6) << ( ok? Gr[0][0][0]: NAN ) << "  KINSOL " << ( ok? Gr[1][0][0]: NAN ) << "  rel " << std::setprecision(1) << d << "\n";
    check( "(C) solve_fsens after KINSOL = after LM to 1e-8", ok && d < 1e-8 ); }

  std::cout << "(D) homotopy under KINSOL: Bratu, lambda continued from 0.5 to 3.4, from zero\n";
  { std::vector<double> X[3];  bool ok[3] = { false, false, false };  int its[3] = { 0, 0, 0 };
    for( int v = 0; v < 3; ++v ){
      FFGraph G;  OCFESLV oc( &G );
      FFVar z = G.add_var("z"), u = G.add_var("u"), lam = G.add_var("lam");
      oc.add_domain( z, FFDom( 0., 1., 4, FFDom::LGL, 6 ) );  oc.add_state( u, {z} );  oc.update_ref( u, 0. );  oc.add_input( lam, {} );
      OCFESLV::EqnOptions io( OCFESLV::EqnRole::INTERIOR, 0 ), ib( OCFESLV::EqnRole::BOUNDARY, 0 );
      oc.add_equation( -OpP( OpP(u,z), z ) - lam*exp( u ), {z}, { FFDom::ALL - FFDom::LB - FFDom::UB }, io );
      oc.add_equation( u, {z}, {FFDom::LB}, ib );  oc.add_equation( u, {z}, {FFDom::UB}, ib );
      oc.options.SOLVE.RES_TOL = 1e-10;  oc.options.DISPLAY_LEVEL = 0;
      if( v ){ oc.options.SOLVE.NLSOLVER = OCFESLV::Options::KINSOL;  oc.options.SOLVE.FACTORIZATION = v == 1? OCFESLV::Options::SOLVE_SUPERLU: OCFESLV::Options::SOLVE_KLU; }
      if( !oc.setup() ) continue;
      oc.add_homotopy( lam, 0.5, 3.4 );  oc.options.HOMOTOPY.ENABLE = true;  oc.options.HOMOTOPY.SOLVE_S0 = true;  oc.options.HOMOTOPY.STEP_CAP = 0.2;
      std::vector<double> var, inp;  oc.init( var, inp );  std::fill( var.begin(), var.end(), 0. );  if( !inp.empty() ) inp[0] = 3.4;
      auto rep = oc.solve( var.data(), inp.empty()? nullptr: inp.data(), nullptr );
      ok[v] = rep.converged;  its[v] = rep.iterations;  X[v] = var;
    }
    double d1 = NAN, d2 = NAN, umax = 0.;
    if( ok[0] ){ for( double x : X[0] ) umax = std::max( umax, x );
      if( ok[1] ){ d1 = 0.; for( size_t i = 0; i < X[0].size(); ++i ) d1 = std::max( d1, std::fabs( X[1][i] - X[0][i] ) ); }
      if( ok[2] ){ d2 = 0.; for( size_t i = 0; i < X[0].size(); ++i ) d2 = std::max( d2, std::fabs( X[2][i] - X[0][i] ) ); } }
    std::cout << "    max state value " << std::scientific << std::setprecision(4) << umax << " (u and its auxiliary d_z u);  iterations LM " << its[0] << ", KINSOL/SuperLU " << its[1]
              << ", KINSOL/KLU " << its[2] << ";  vs LM " << std::setprecision(1) << d1 << ", " << d2 << "\n";
    check( "(D) LM converges along the homotopy", ok[0] );
    check( "(D) KINSOL/SUPERLU and KINSOL/KLU converge along it to LM's solution, to 1e-10", ok[1] && ok[2] && d1 < 1e-10 && d2 < 1e-10 ); }

  std::cout << "(E) a singular start: a^3 + a = a*^3 + a*, p a = p* a*, both started at zero\n";
  { double err[2] = { NAN, NAN };  bool ok[2] = { false, false };
    for( int v = 0; v < 2; ++v ){
      FFGraph G;  OCFESLV oc( &G );
      FFVar xi = G.add_var("xi"), a = G.add_var("a"), p = G.add_var("p");
      oc.add_domain( xi, FFDom( 0., 1., 4, FFDom::LGL, 5 ) );  oc.add_state( a, {xi} );  oc.add_state( p, {xi} );
      oc.update_ref( a, 0. );  oc.update_ref( p, 0. );
      FFVar const av = 1. + xi, pv = 2. - xi;
      OCFESLV::EqnOptions io( OCFESLV::EqnRole::INTERIOR, 0 );
      oc.add_equation( a*a*a + a - ( av*av*av + av ), {xi}, {FFDom::ALL}, io );
      oc.add_equation( p*a - pv*av, {xi}, {FFDom::ALL}, io );
      oc.add_output( OpI( ( a - av )*( a - av ) + ( p - pv )*( p - pv ), xi ) );
      oc.options.SOLVE.RES_TOL = 1e-10;  oc.options.DISPLAY_LEVEL = 0;
      if( v ) oc.options.SOLVE.NLSOLVER = OCFESLV::Options::KINSOL;   // SOLVE_SUPERLU, the default factorisation
      if( !oc.setup() ) continue;
      std::vector<double> var, inp;  oc.init( var, inp );
      auto rep = oc.solve( var.data(), inp.empty()? nullptr: inp.data(), nullptr );
      ok[v] = rep.converged;  if( ok[v] && !oc.val_function().empty() ) err[v] = std::sqrt( std::fabs( oc.val_function()[0] ) );
    }
    std::cout << "    L2 error against a*, p*: LM " << std::scientific << std::setprecision(1) << err[0] << ", KINSOL/SUPERLU " << err[1] << "\n";
    check( "(E) LM converges from the singular start", ok[0] && err[0] < 1e-8 );
    check( "(E) KINSOL/SUPERLU converges from it too (the QR's minimum-norm step at the singular LU), to 1e-8", ok[1] && err[1] < 1e-8 ); }

  std::cout << "(F) a structurally rank-deficient system, the audit opted past: KINSOL takes the QR\n";
  { FFGraph G;  OCFESLV oc( &G );                // OCFE_PDE20e's case D: x' = y, 0 = x - a(t) on t > 0 only, y pinned twice
    FFVar t = G.add_var("t"), x = G.add_var("x"), y = G.add_var("y");
    FFVar const a = 1. + 0.5*t - 0.3*t*t, ap = 0.5 - 0.6*t;
    oc.add_domain( t, FFDom( 0., 0.5, 3, FFDom::LGR, 4 ) );
    oc.add_state( x, {t} );  oc.add_state( y, {t} );  oc.update_ref( x, 1. );  oc.update_ref( y, 0. );
    OCFESLV::EqnOptions io( OCFESLV::EqnRole::INTERIOR, 0 ), ii( OCFESLV::EqnRole::INITIAL, 0 );
    int const T_INT = FFDom::ALL - FFDom::LB;
    oc.add_equation( OpP(x,t) - y, {t}, {T_INT}, io );
    oc.add_equation( x - a, {t}, {T_INT}, io );
    oc.add_equation( y - ap, {t}, {FFDom::LB}, ii );  oc.add_equation( 2.*y - 2.*ap, {t}, {FFDom::LB}, ii );
    oc.set_evolution_domain( t );
    oc.options.REDUCE.ORDER = OCFESLV::Options::RED_MAIN;  oc.options.INTERFACE.IMPOSITION = OCFESLV::Options::IC_WEAK;
    oc.options.REDUCE.HIDDEN_IC = false;  oc.options.FATAL.REDUCED_DOF = false;  oc.options.SOLVE.MARCHING = false;
    oc.options.SOLVE.NLSOLVER = OCFESLV::Options::KINSOL;  oc.options.DISPLAY_LEVEL = 0;
    bool const su = oc.setup();
    auto const& au = oc.reduced_dof_audit();
    std::vector<double> var, inp;  if( su ) oc.init( var, inp );
    OCFESLV::SolveReport rep;  if( su ) rep = oc.solve( var.data(), inp.empty()? nullptr: inp.data(), nullptr );
    std::cout << "    setup " << su << ", audit ran " << au.ran << " rank " << au.rank << "/" << au.rows << " deficiency " << au.deficiency
              << ";  KINSOL factorised with " << ( rep.linear_solver.empty()? "-": rep.linear_solver ) << "\n";
    check( "(F) the case is what it is meant to be: set up, the audit's structural deficiency > 0", su && au.ran && au.deficiency > 0 );
    check( "(F) KINSOL took the QR directly (SuperLU passed BLAS an invalid dimension on such a matrix)", su && rep.linear_solver == "SPQR" ); }

  std::cout << "==========  OCFE_KINSOL: " << ( nfail? std::to_string( nfail ) + " FAILED": "ALL PASS" ) << "  ==========\n";
  return nfail? 1: 0;
}
