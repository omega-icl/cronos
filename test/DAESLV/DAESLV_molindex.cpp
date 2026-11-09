////////////////////////////////////////////////////////////////////////
// DAESLV_molindex -- phase 9 of the MoL plan (8 Nov 2026): the index of the SEMI-DISCRETE system when the interface plan
// mints multipliers, and the null-space elimination that brings it back to 1.
//
//   (A) the rod u_t = D u_zz, u(t,0) = 0, u_z(t,1) = 0, u0 = sin(pi z/2), D = 0.1 p: u(1,1) = exp(-0.1 p (pi/2)^2).
//       Under IC_TRACE every value-continuity multiplier enters only the two interface copies' differential rows and its
//       constraint reads differential unknowns only -- Hessenberg index 2 (KINSOL's IC Jacobian singular before this
//       work).  The null-space elimination merges the copies and sums their rows: the merged operator's spectrum is
//       IC_STRONG's to every digit, and so is the solution.  Value, dJ/dp forward and adjoint against the closed form
//       (spatial error ~1e-5 on 3 elements of LGL(5): the gate is the IC_STRONG = IC_TRACE agreement at tolerance, and
//       both at the discretisation's distance from the closed form).
//   (C) the PSA column of DAESLV_molsens (5 fields, convective-dispersive with a Danckwerts inlet, LDF, an energy
//       balance, a nonlinear isotherm) under IC_TRACE through MoL: 6 value multipliers eliminated, the c1 and c2 balances'
//       rows carrying a second mass entry (F d_t q) -- the case that showed the copies must be named by the constraint,
//       not by the rows' masses.  = IC_STRONG to 4e-11 / 4e-10, and to OCFESLV at order 6 as IC_STRONG is.
//   (B) model G of DAESLV_molsetup under IC_STRONG, on 3 elements: its 6 kept-explicit multipliers (three per interface)
//       are pinned by algebraic rows (index 1 already) and must be left alone -- n_tau stays 6 -- and the solution is the initial profile
//       exactly (u_t = u_zz + lam = 0): u(1, 0.5) = 0.25, lam = 2 at every node.
////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <iomanip>
#include <cmath>
#include "ocfeslv.hpp"
#include "daeslv_idas.hpp"
using namespace mc;

static int nfail = 0;
static void check( std::string const& w, bool ok ){ std::cout << "  " << std::left << std::setw(96) << w << ( ok? "PASS": "FAIL" ) << "\n"; if( !ok ) ++nfail; }
static double rel( double a, double b ){ return std::fabs( a - b ) / std::max( 1e-300, std::fabs( b ) ); }

// ---- the PSA column of DAESLV_molsens (copied verbatim, the imposition made a switch), for case (C)
static size_t G_NST = 0, G_NTAU = 0;
static double const U_VEL = 1.0, D_ax = 0.1, F_ph = 0.5, k_ldf = 2.0;
static double const qs1 = 1.0, qs2 = 1.0, b02 = 1.0, beta1 = 2.0, beta2 = 1.0, T0_ref = 1.0;
static double const Cp_e = 1.0, G_cv = 1.0, lam = 0.1, dH1 = 1.0, dH2 = 0.5, hw = 0.5, Tw = 1.0;
static double const c0_1 = 0.5, c0_2 = 0.5, tau_in = 0.15, T_END = 5.0;

struct Fields { FFVar t, z, c1, c2, q1, q2, T, b01; };

// the collocation options, the same for both solvers -- OCFESLV holds them itself, DAESLV for its semi-discretisation
static void colloc_options( OCFESLV& oc )
{ oc.options.INTERFACE.IMPOSITION = OCFESLV::Options::IC_STRONG;  oc.options.INTERFACE.SAT_SIGMA0 = 10.0; }
static int G_IMP = 0;   // 0 IC_STRONG, 1 IC_TRACE: the MoL imposition of case (C)
static void colloc_options( DAESLV_IDAS& S )
{ S.options.INTERFACE.IMPOSITION = G_IMP? OCFEMOL::Options::IC_TRACE: OCFEMOL::Options::IC_STRONG;  S.options.INTERFACE.SAT_SIGMA0 = 10.0; }

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

static void tight( DAESLV_IDAS& S )
{
  S.options.SOLVE.RTOL  = 1e-10;  S.options.SOLVE.ATOL  = 1e-12;
  S.options.SOLVE.RTOLB = 1e-10;  S.options.SOLVE.ATOLB = 1e-12;
  S.options.SOLVE.VERBOSE = 0;    S.options.DISPLAY_LEVEL = 0;
}

// DAESLV on the PSA column at b01: mode 0 value, 1 forward, 2 adjoint; f = (Inv1, Eff1), g = their d/db01
static bool psa_daeslv( double b01v, int mode, std::vector<double>& f, std::vector<double>& g, std::string& why )
{
  FFGraph G;  DAESLV_IDAS S( &G );  Fields F;
  declare_psa( S, G, F, 3, 5, 5, 3 );
  tight( S );
  if( !S.setup() ){ why = "setup: " + S.extract_error(); return false; }
  if( mode == 0 ){ auto const& SD = S.semidiscrete_model()->semi_discrete();  G_NST = SD.n_states;  G_NTAU = SD.n_tau; }
  std::vector<FFModel::InputVal> const in{ { F.b01, { b01v } } };
  auto const st = mode == 0? S.solve( in ): mode == 1? S.solve_fsens( in ): S.solve_asens( in );
  if( st != DAESLV_IDAS::NORMAL ){ why = "solve: " + S.extract_error(); return false; }
  f = S.val_function();
  if( mode ){ auto const& Gr = S.val_function_gradient();            // [ parameter ][ function ]
              if( Gr.size() != 1 || Gr[0].size() != 2 ){ why = "gradient shape"; return false; }
              g = { Gr[0][0], Gr[0][1] }; }
  return f.size() == 2;
}

// OCFESLV on the same declaration, 10 time elements of order NODT -- DAESLV_molpsa's grid, where the fields agree with
// MoL to 1e-6 at order 6 (one element of order 6 does not resolve the breakthrough: 3.8% on Eff1, measured)
static bool psa_ocfeslv( double b01v, size_t NODT, std::vector<double>& f, std::vector<double>& g, std::string& why )
{
  FFGraph G;  OCFESLV oc( &G );  Fields F;
  declare_psa( oc, G, F, 3, 5, 10, NODT );
  oc.options.SOLVE.MARCHING = false;
#if defined(CRONOS__WITH_SPQR)
  oc.options.SOLVE.FACTORIZATION = OCFESLV::Options::SOLVE_SPQR;
#endif
  if( !oc.setup() ){ why = "setup"; return false; }
  std::vector<double> xv, inp;
  if( !oc.init( xv, inp, nullptr ) ){ why = "init"; return false; }
  oc.set_input_values( F.b01, { b01v }, inp.data() );
  oc.register_control( F.b01 );
  if( !oc.solve_fsens( xv.data(), inp.data(), nullptr ) ){ why = "solve_fsens"; return false; }
  auto const& J = oc.sens_jacobian();  f = oc.sens_functions();
  if( f.size() != 2 || J.size() != 2 ){ why = "shape"; return false; }
  g = { J[0], J[1] };                                              // n_fct x n_control_dof, one control DOF
  return true;
}


int main()
{
  std::cout << std::scientific << std::setprecision(4);
  std::cout << "(A) the rod under IC_STRONG and IC_TRACE: value and dJ/dp, forward and adjoint\n";
  double const P = 1.3, TEND = 1., ex = std::exp( -0.1*P*0.25*M_PI*M_PI*TEND ), dex = -0.1*0.25*M_PI*M_PI*TEND*ex;
  double F[2] = { 0., 0. }, Gf[2] = { 0., 0. }, Ga[2] = { 0., 0. };  size_t ntau[2] = { 0, 0 }, nst[2] = { 0, 0 };
  for( int imp = 0; imp < 2; ++imp ){
    FFGraph G; DAESLV_IDAS S( &G ); FFPartial OpP; FFEval OpE;
    FFVar t = G.add_var("t"), z = G.add_var("z"), u = G.add_var("u"), p = G.add_var("p");
    S.add_domain( t, FFDom( 0., TEND, 2, FFDom::LGR, 3 ) ); S.add_domain( z, FFDom( 0., 1., 3, FFDom::LGL, 5 ) );
    S.set_evolution_domain( t ); S.add_state( u, {t,z} ); S.update_ref( u, 0.5 ); S.add_input( p, {} ); S.update_ref( p, 1. );
    int const TI = FFDom::ALL - FFDom::LB, ZI = FFDom::ALL - FFDom::LB - FFDom::UB;
    DAESLV_IDAS::EqnOptions io( DAESLV_IDAS::EqnRole::INTERIOR, 0 ), ib( DAESLV_IDAS::EqnRole::BOUNDARY, 0 ), ii( DAESLV_IDAS::EqnRole::INITIAL, 0 );
    S.add_equation( OpP(u,t) - 0.1*p*OpP(OpP(u,z),z), {t,z}, {TI, ZI}, io );
    S.add_equation( u, {t,z}, {TI, FFDom::LB}, ib ); S.add_equation( OpP(u,z), {t,z}, {TI, FFDom::UB}, ib );
    S.add_equation( u - sin( 0.5*M_PI*z ), {t,z}, {FFDom::LB, FFDom::ALL}, ii );
    S.add_output( OpE( u, {{t,1},{z,1}}, {{t,TEND},{z,1.}} ) );
    S.options.INTERFACE.IMPOSITION = imp? OCFEMOL::Options::IC_TRACE: OCFEMOL::Options::IC_STRONG; S.options.INTERFACE.SAT_SIGMA0 = 10.;
    S.options.SOLVE.RTOL = S.options.SOLVE.ATOL = 1e-10; S.options.SOLVE.RTOLB = S.options.SOLVE.ATOLB = 1e-10; S.options.SOLVE.VERBOSE = 0; S.options.DISPLAY_LEVEL = 0;
    char const* name = imp? "IC_TRACE": "IC_STRONG";
    if( !S.setup() ){ check( std::string( name ) + ": setup -- " + S.extract_error(), false ); continue; }
    ntau[imp] = S.semidiscrete_model()->semi_discrete().n_tau;  nst[imp] = S.semidiscrete_model()->semi_discrete().n_states;
    std::vector<FFModel::InputVal> in{ { p, { P } } };
    bool ok = S.solve( in ) == DAESLV_IDAS::NORMAL;  if( ok ) F[imp] = S.val_function()[0];
    ok = ok && S.solve_fsens( in ) == DAESLV_IDAS::NORMAL;  if( ok ) Gf[imp] = S.val_function_gradient()[0][0];
    ok = ok && S.solve_asens( in ) == DAESLV_IDAS::NORMAL;  if( ok ) Ga[imp] = S.val_function_gradient()[0][0];
    std::cout << "    " << name << ": " << nst[imp] << " unknowns, " << ntau[imp] << " explicit multipliers;  u(1,1) " << F[imp] << " (closed form " << ex
              << ")  dJ/dp fwd " << Gf[imp] << " adj " << Ga[imp] << " (closed form " << dex << ")\n";
    check( std::string( name ) + ": set up and solved, value, forward and adjoint", ok );
  }
  check( "IC_STRONG: no multipliers; IC_TRACE: the 2 value multipliers eliminated, the 2 flux multipliers kept, the same unknown count",
         ntau[0] == 0 && ntau[1] == 2 && nst[0] == nst[1] );
  check( "IC_TRACE = IC_STRONG: value to 1e-8, forward and adjoint gradients to 1e-6", rel( F[1], F[0] ) < 1e-8 && rel( Gf[1], Gf[0] ) < 1e-6 && rel( Ga[1], Ga[0] ) < 1e-6 );
  check( "both = closed form to the discretisation's 1e-4: value and gradients", rel( F[0], ex ) < 1e-4 && rel( F[1], ex ) < 1e-4 && rel( Gf[1], dex ) < 1e-4 && rel( Ga[1], dex ) < 1e-4 );

  std::cout << "(B) model G under IC_STRONG: kept-explicit multipliers left alone\n";
  {
    FFGraph G; DAESLV_IDAS S( &G ); FFPartial OpP; FFEval OpE;
    FFVar t = G.add_var("t"), z = G.add_var("z"), u = G.add_var("u"), lam = G.add_var("lam");
    S.add_domain( t, FFDom( 0., 1., 2, FFDom::LGR, 3 ) ); S.add_domain( z, FFDom( 0., 1., 3, FFDom::LGL, 4 ) );
    S.set_evolution_domain( t );
    S.add_state( u, {t,z} ); S.update_ref( u, 0.25 ); S.add_state( lam, {t,z} ); S.update_ref( lam, 2. );
    int const TI = FFDom::ALL - FFDom::LB, ZI = FFDom::ALL - FFDom::LB - FFDom::UB;
    DAESLV_IDAS::EqnOptions io( DAESLV_IDAS::EqnRole::INTERIOR, 0 ), ib( DAESLV_IDAS::EqnRole::BOUNDARY, 0 ), ii( DAESLV_IDAS::EqnRole::INITIAL, 0 );
    S.add_equation( OpP(u,t) - OpP(OpP(u,z),z) - lam, {t,z}, {TI, ZI}, io );
    S.add_equation( OpP(OpP(u,z),z) + lam,            {t,z}, {TI, FFDom::ALL}, io );
    S.add_equation( u, {t,z}, {TI, FFDom::LB}, ib ); S.add_equation( u, {t,z}, {TI, FFDom::UB}, ib );
    S.add_equation( u - z*(1.-z), {t,z}, {FFDom::LB, FFDom::ALL}, ii );
    S.add_output( OpE( u, {{t,1},{z,1}}, {{t,1.},{z,0.5}} ) );  S.add_output( OpE( lam, {{t,1},{z,1}}, {{t,1.},{z,0.5}} ) );
    S.options.INTERFACE.IMPOSITION = OCFEMOL::Options::IC_STRONG;
    S.options.SOLVE.RTOL = S.options.SOLVE.ATOL = 1e-10; S.options.SOLVE.VERBOSE = 0; S.options.DISPLAY_LEVEL = 0;
    bool ok = S.setup();
    size_t const nt = ok? S.semidiscrete_model()->semi_discrete().n_tau: 0;
    ok = ok && S.solve( std::vector<FFModel::InputVal>{} ) == DAESLV_IDAS::NORMAL;
    double const uv = ok? S.val_function()[0]: 0., lv = ok? S.val_function()[1]: 0.;
    std::cout << "    " << nt << " explicit multipliers kept;  u(1, 0.5) " << uv << " (0.25)  lam(1, 0.5) " << lv << " (2)" << ( ok? "": "  -- " + S.extract_error() ) << "\n";
    check( "model G: set up and solved with its kept-explicit multipliers untouched (6: three per interface)", ok && nt == 6 );
    check( "model G: u(1, 0.5) = 0.25 and lam = 2, to 1e-8 (the initial profile, exactly)", ok && rel( uv, 0.25 ) < 1e-8 && rel( lv, 2. ) < 1e-8 );
  }
  std::cout << "(C) the PSA column of DAESLV_molsens under IC_STRONG and IC_TRACE through MoL: Inv1, Eff1 and d/db01, forward and adjoint\n";
  {
    std::vector<double> f[2][3], g[2][3];  size_t nst[2] = { 0, 0 }, ntau[2] = { 0, 0 };  bool ok[2] = { true, true };
    for( G_IMP = 0; G_IMP < 2; ++G_IMP ){
      std::string why;
      for( int mode = 0; mode < 3 && ok[G_IMP]; ++mode ) ok[G_IMP] = psa_daeslv( 3.0, mode, f[G_IMP][mode], g[G_IMP][mode], why );
      nst[G_IMP] = G_NST;  ntau[G_IMP] = G_NTAU;
      char const* name = G_IMP? "IC_TRACE": "IC_STRONG";
      if( ok[G_IMP] ) std::cout << "    " << name << ": " << nst[G_IMP] << " unknowns, " << ntau[G_IMP] << " explicit multipliers;  Inv1 " << f[G_IMP][0][0] << " Eff1 " << f[G_IMP][0][1]
                                << "  d/db01 fwd " << g[G_IMP][1][0] << " " << g[G_IMP][1][1] << "  adj " << g[G_IMP][2][0] << " " << g[G_IMP][2][1] << "\n";
      check( std::string( name ) + ": set up and solved, value, forward and adjoint" + ( ok[G_IMP]? std::string(): " -- " + why ), ok[G_IMP] );
    }
    std::vector<double> fo, go;  std::string why;  bool const oko = psa_ocfeslv( 3.0, 6, fo, go, why );
    check( "OCFESLV at time order 6 (IC_STRONG), the reference" + ( oko? std::string(): " -- " + why ), oko );
    if( ok[0] && ok[1] ){
      auto gap = []( std::vector<double> const& x, std::vector<double> const& y ){ double m = 0.; for( size_t i = 0; i < x.size() && i < y.size(); ++i ) m = std::max( m, rel( x[i], y[i] ) ); return m; };
      double const ev = gap( f[1][0], f[0][0] ), ef = gap( g[1][1], g[0][1] ), ea = gap( g[1][2], g[0][2] ), eaf = gap( g[1][2], g[1][1] );
      std::cout << "    IC_TRACE vs IC_STRONG: values " << ev << "  forward " << ef << "  adjoint " << ea << ";  IC_TRACE adjoint vs forward " << eaf << "\n";
      check( "IC_TRACE: the 6 value multipliers (c1, c2, T at 2 interfaces) eliminated, the 6 flux multipliers kept, the same unknown count",
             ntau[0] == 0 && ntau[1] == 6 && nst[0] == nst[1] );
      check( "IC_TRACE = IC_STRONG: values to 1e-9, forward and adjoint gradients to 1e-8; adjoint = forward to 1e-8", ev < 1e-9 && ef < 1e-8 && ea < 1e-8 && eaf < 1e-8 );
      if( oko ){
        double const ov = gap( f[1][0], fo ), og = gap( g[1][1], go );
        std::cout << "    IC_TRACE vs OCFESLV(6): values " << ov << "  gradients " << og << "  (DAESLV_molsens records 4.6e-5 / 5.6e-6 for IC_STRONG)\n";
        check( "IC_TRACE = OCFESLV at order 6 to 1e-4 (values) and 1e-5 (gradients), as IC_STRONG", ov < 1e-4 && og < 1e-5 );
      }
    }
  }
  std::cout << "==========  DAESLV_molindex: " << ( nfail? std::to_string( nfail ) + " FAILED": "ALL PASS" ) << "  ==========\n";
  return nfail? 1: 0;
}
