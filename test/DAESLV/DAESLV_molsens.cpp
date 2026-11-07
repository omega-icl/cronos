// DAESLV_molsens -- MoL phase 6a's gate, with 6b: a model distributed in space declared ONCE with inputs and outputs,
// solved by DAESLV (the method of lines) and by OCFESLV (collocating t and z); values and sensitivities.
//
// (A) the PSA column of OCFE_PDE36c / DAESLV_molpsa, the isotherm affinity b01 an INPUT, its two outputs: the bed
//     inventory Inv1 = int (c1 + F q1) dz at T_END and the cumulative effluent Eff1 = int U c1 dt at z = 1.  DAESLV's
//     forward and adjoint gradients against each other and against a central difference of DAESLV's own values; values
//     and forward gradients against OCFESLV on 10 time elements of order 6 (DAESLV_molpsa: the fields agree to 1e-6).
// (B) a CONSERVED inventory: u_t = alpha u_zz + p e^-t, zero flux at both ends, u0 = p2 + cos(pi z), J = int u dz at T.
//     The diffusion conserves J, so the costate is (nearly) constant while the parameter quadrature ~ e^-t varies --
//     the case that stalled IDAS's order selection (DAESLV_quadacc).  A uniform field is in the discrete diffusion
//     operator's kernel, so du/dp = 1 - e^-t at every node and dJ/dp = 1 - e^-1, dJ/dp2 = 1 EXACTLY on any mesh, to the
//     integration tolerance; p2 enters through an INITIAL row (initial sensitivities through input-binding rows), and
//     the two inputs check that the gradient's rows are the declared inputs, in order.
// (C) an input in a DIRICHLET condition -- an ALGEBRAIC node's value at t0 depends on it: u_t = alpha u_zz, u(0,t) = p,
//     u_z(1,t) = 0, u0 = p + A sin(pi z / 2), J = int u dz at T.  A uniform field satisfies both conditions and is in the
//     operator's kernel, so dJ/dp = 1 exactly; J is linear in A, so a central difference with a LARGE step has no
//     truncation error.  Whether D1 (the algebraic derivatives at t0) was needed for sensitivities: it is not (measured
//     4 Nov 2026, here and in A, whose Danckwerts inlet is algebraic and time-varying at t0).
#include <iostream>
#include <iomanip>
#include <cmath>
#include "ocfeslv.hpp"
#include "daeslv_idas.hpp"
using namespace mc;

static int nfail = 0;
static void check( std::string const& w, bool ok )
{ std::cout << "  " << std::left << std::setw(86) << w << ( ok? "PASS": "FAIL" ) << "\n"; if( !ok ) ++nfail; }
static double rel( double a, double b ){ return std::fabs( a - b )/std::max( 1e-300, std::fabs( b ) ); }

// OCFE_PDE36c's constants, the isotherm affinity at its nominal value
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

// (B) the conserved-inventory rod: J = int u dz at T, inputs p (a source) and p2 (the initial level)
static bool rod_daeslv( int mode, double& J, double& dJdp, double& dJdp2, std::string& why )
{
  double const ALPHA = 0.1, TEND = 1., P = 0.7, P2 = 1.3;
  FFGraph G;  DAESLV_IDAS S( &G );
  FFVar t = G.add_var("t"), z = G.add_var("z"), u = G.add_var("u"), p = G.add_var("p"), p2 = G.add_var("p2");
  S.add_domain( t, FFDom( 0., TEND, 2, FFDom::LGR, 3 ) );
  S.add_domain( z, FFDom( 0., 1., 3, FFDom::LGL, 5 ) );
  S.set_evolution_domain( t );
  S.add_state( u, {t,z} ); S.update_ref( u, 1. );
  S.add_input( p, {} );  S.add_input( p2, {} );
  FFPartial OpP;  FFIntegral OpI;
  int const TI = FFDom::ALL - FFDom::LB, ZI = FFDom::ALL - FFDom::LB - FFDom::UB;
  DAESLV_IDAS::EqnOptions io( DAESLV_IDAS::EqnRole::INTERIOR, 0 ), ib( DAESLV_IDAS::EqnRole::BOUNDARY, 0 ), ii( DAESLV_IDAS::EqnRole::INITIAL, 0 );
  S.add_equation( OpP(u,t) - OpP( ALPHA*OpP(u,z), z ) - p*exp( -t ), {t,z}, {TI, ZI},                io );
  S.add_equation( OpP(u,z),                                        {t,z}, {TI, FFDom::LB},         ib );
  S.add_equation( OpP(u,z),                                        {t,z}, {TI, FFDom::UB},         ib );
  S.add_equation( u - ( p2 + cos( M_PI*z ) ),                      {t,z}, {FFDom::LB, FFDom::ALL}, ii );
  S.add_output( OpI( u, z ), {t}, {TEND} );
  tight( S );
  if( !S.setup() ){ why = "setup: " + S.extract_error(); return false; }
  std::vector<FFModel::InputVal> const in{ { p, { P } }, { p2, { P2 } } };
  auto const st = mode == 0? S.solve( in ): mode == 1? S.solve_fsens( in ): S.solve_asens( in );
  if( st != DAESLV_IDAS::NORMAL ){ why = "solve: " + S.extract_error(); return false; }
  J = S.val_function().at( 0 );
  if( mode ){ auto const& Gr = S.val_function_gradient();
              if( Gr.size() != 2 || Gr[0].size() != 1 ){ why = "gradient shape"; return false; }
              dJdp = Gr[0][0];  dJdp2 = Gr[1][0]; }
  return true;
}

// (C) the Dirichlet rod
static bool dirichlet_daeslv( int mode, double P, double A, double& J, double& gp, double& gA, std::string& why )
{
  double const ALPHA = 0.1, TEND = 1.;
  FFGraph G;  DAESLV_IDAS S( &G );
  FFVar t = G.add_var("t"), z = G.add_var("z"), u = G.add_var("u"), p = G.add_var("p"), a = G.add_var("a");
  S.add_domain( t, FFDom( 0., TEND, 2, FFDom::LGR, 3 ) );
  S.add_domain( z, FFDom( 0., 1., 3, FFDom::LGL, 5 ) );
  S.set_evolution_domain( t );
  S.add_state( u, {t,z} ); S.update_ref( u, 1. );
  S.add_input( p, {} );  S.add_input( a, {} );
  FFPartial OpP;  FFIntegral OpI;
  int const TI = FFDom::ALL - FFDom::LB, ZI = FFDom::ALL - FFDom::LB - FFDom::UB;
  DAESLV_IDAS::EqnOptions io( DAESLV_IDAS::EqnRole::INTERIOR, 0 ), ib( DAESLV_IDAS::EqnRole::BOUNDARY, 0 ), ii( DAESLV_IDAS::EqnRole::INITIAL, 0 );
  S.add_equation( OpP(u,t) - OpP( ALPHA*OpP(u,z), z ), {t,z}, {TI, ZI},                io );
  S.add_equation( u - p,                               {t,z}, {TI, FFDom::LB},         ib );
  S.add_equation( OpP(u,z),                            {t,z}, {TI, FFDom::UB},         ib );
  S.add_equation( u - ( p + a*sin( 0.5*M_PI*z ) ),     {t,z}, {FFDom::LB, FFDom::ALL}, ii );
  S.add_output( OpI( u, z ), {t}, {TEND} );
  tight( S );
  if( !S.setup() ){ why = "setup: " + S.extract_error(); return false; }
  std::vector<FFModel::InputVal> const in{ { p, { P } }, { a, { A } } };
  auto const st = mode == 0? S.solve( in ): mode == 1? S.solve_fsens( in ): S.solve_asens( in );
  if( st != DAESLV_IDAS::NORMAL ){ why = "solve: " + S.extract_error(); return false; }
  J = S.val_function().at( 0 );
  if( mode ){ auto const& Gr = S.val_function_gradient();
              if( Gr.size() != 2 ){ why = "gradient shape"; return false; }
              gp = Gr[0].at( 0 );  gA = Gr[1].at( 0 ); }
  return true;
}

int main()
{
  std::cout << std::scientific << std::setprecision(3);
  // ---- (A)
  std::cout << "(A) PSA column, b01 an input; outputs Inv1(T), Eff1(z = 1)\n";
  double const B = 3.0, H = 1e-4;
  std::vector<double> f0, gF, gA, fp, fm, dummy, fo, go;  std::string why;
  bool const v  = psa_daeslv( B, 0, f0, dummy, why );
  check( "(A) DAESLV: the PDE model with an input and two outputs, solved by name" + ( v? std::string(): " -- " + why ), v );
  bool const fw = psa_daeslv( B, 1, fp, gF, why );
  check( "(A) DAESLV forward sensitivities by name" + ( fw? std::string(): " -- " + why ), fw );
  bool const ad = psa_daeslv( B, 2, fp, gA, why );
  check( "(A) DAESLV adjoint sensitivities by name" + ( ad? std::string(): " -- " + why ), ad );
  bool const fd = psa_daeslv( B + H, 0, fp, dummy, why ) && psa_daeslv( B - H, 0, fm, dummy, why );
  // OCFESLV's time error is the gap: a ladder in its time order, as DAESLV_molpsa's gate (a convergence, not a tolerance)
  std::vector<std::vector<double>> fL( 3 ), gL( 3 );  bool ocL[3];
  for( int k = 0; k < 3; ++k ) ocL[k] = psa_ocfeslv( B, 4 + k, fL[k], gL[k], why );
  bool const oc = ocL[0] && ocL[1] && ocL[2];
  if( oc ){ fo = fL[2];  go = gL[2]; }
  check( "(A) OCFESLV at time orders 4, 5, 6: values and forward sensitivities" + ( oc? std::string(): " -- " + why ), oc );
  if( v && fw && fd ){                                  // each comparison needs only the solves it reads
    for( int i = 0; i < 2; ++i ){
      double const gd = ( fp[i] - fm[i] )/( 2.*H );
      std::cout << "    " << ( i? "Eff1": "Inv1" ) << ": DAESLV " << f0[i] << ( oc? "  OCFESLV " + std::to_string( fo[i] ): std::string() )
                << "   d/db01: fwd " << gF[i] << "  adj " << ( ad? std::to_string( gA[i] ): std::string( "--" ) ) << "  FD " << gd
                << ( oc? "  OCFESLV " + std::to_string( go[i] ): std::string() ) << "\n";
    }
    double eFA = 0., eFD = 0., eV = 0., eG = 0.;
    for( int i = 0; i < 2; ++i ){
      double const gd = ( fp[i] - fm[i] )/( 2.*H );
      eFD = std::max( eFD, rel( gF[i], gd ) );
      if( ad ) eFA = std::max( eFA, rel( gA[i], gF[i] ) );
      if( oc ){ eV = std::max( eV, rel( f0[i], fo[i] ) );  eG = std::max( eG, rel( gF[i], go[i] ) ); }
    }
    std::cout << "    max rel: adj-fwd " << ( ad? eFA: NAN ) << "  fwd-FD " << eFD << "  values DAESLV-OCFESLV " << eV << "  gradients DAESLV-OCFESLV " << eG << "\n";
    // the central difference of two ADAPTIVE solves carries their uncorrelated global errors, ~1e-9 at RTOL 1e-10,
    // over H = 1e-4: noise ~1e-5 by construction, so 1e-5 sat at the floor (2.6e-6, 5.5e-6, 1.0e-5 across changes that
    // left the forward sensitivity, the adjoint and the OCFESLV comparison unchanged -- 7 Nov 2026).  A wrong
    // sensitivity is wrong at order one; the sharp checks are adjoint = forward and the convergence to OCFESLV.
    check( "(A) DAESLV forward = central difference of its own values, to 1e-4 (the FD's noise ~1e-5)", eFD <= 1e-4 );
    bool conv = oc;
    if( oc ){
      std::cout << std::setprecision(10);
      for( int i = 0; i < 2; ++i ){
        std::cout << "    " << ( i? "Eff1": "Inv1" ) << "  DAESLV " << f0[i] << "   |DAESLV - OCFESLV| / |.|, orders 4 5 6:  values";
        for( int k = 0; k < 3; ++k ) std::cout << " " << std::setprecision(2) << rel( f0[i], fL[k][i] );
        std::cout << "   gradients";
        for( int k = 0; k < 3; ++k ) std::cout << " " << rel( gF[i], gL[k][i] );
        std::cout << "\n" << std::setprecision(10);
        for( int k = 0; k < 2; ++k ) conv = conv && rel( f0[i], fL[k][i] ) >= 3.*rel( f0[i], fL[k+1][i] );
      }
      std::cout << std::setprecision(3);
    }
    check( "(A) values: the gap to OCFESLV falls >= 3x per time order (it is OCFESLV's time error)", conv );
    check( "(A) at order 6: values and gradients DAESLV = OCFESLV, to 1e-4",   oc && eV <= 1e-4 && eG <= 1e-4 );
    check( "(A) DAESLV adjoint = forward, to 1e-6",                            ad && eFA <= 1e-6 );
  }
  // ---- (B)
  std::cout << "(B) conserved inventory: u_t = alpha u_zz + p e^-t, zero flux, u0 = p2 + cos(pi z), J = int u dz at T\n";
  double const P = 0.7, P2 = 1.3, e1 = 1. - std::exp( -1. ), Jex = P2 + P*e1;
  double J0 = 0., Jf = 0., Ja = 0., dpF = 0., dp2F = 0., dpA = 0., dp2A = 0., d0 = 0.;
  bool const rv = rod_daeslv( 0, J0, d0, d0, why );
  check( "(B) value solve" + ( rv? std::string(): " -- " + why ), rv );
  bool const rf = rod_daeslv( 1, Jf, dpF, dp2F, why );
  check( "(B) forward sensitivities" + ( rf? std::string(): " -- " + why ), rf );
  bool const ra = rod_daeslv( 2, Ja, dpA, dp2A, why );
  check( "(B) adjoint sensitivities: the constant-costate case completes" + ( ra? std::string(): " -- " + why ), ra );
  if( rv && rf && ra ){
    std::cout << "    J " << J0 << " (closed form " << Jex << ")   dJ/dp: fwd " << dpF << " adj " << dpA << " (exact " << e1
              << ")   dJ/dp2: fwd " << dp2F << " adj " << dp2A << " (exact 1)\n";
    check( "(B) J = p2 + p (1 - e^-1), to the spatial error (1e-6)",          rel( J0, Jex ) <= 1e-6 );
    check( "(B) dJ/dp = 1 - e^-1 and dJ/dp2 = 1, forward, to 1e-8",          rel( dpF, e1 ) <= 1e-8 && rel( dp2F, 1. ) <= 1e-8 );
    check( "(B) dJ/dp = 1 - e^-1 and dJ/dp2 = 1, adjoint, to 1e-8",          rel( dpA, e1 ) <= 1e-8 && rel( dp2A, 1. ) <= 1e-8 );
  }
  // ---- (C)
  std::cout << "(C) an input in a Dirichlet condition: an algebraic node at t0 depends on it\n";
  { double const PC = 0.4, AC = 1.2, HC = 0.1;
    double Jc = 0., Jp = 0., Jm = 0., gpF = 0., gAF = 0., gpA = 0., gAA = 0., d = 0.;
    bool const okv = dirichlet_daeslv( 0, PC, AC, Jc, d, d, why ) && dirichlet_daeslv( 0, PC, AC + HC, Jp, d, d, why )
                  && dirichlet_daeslv( 0, PC, AC - HC, Jm, d, d, why );
    bool const okf = okv && dirichlet_daeslv( 1, PC, AC, Jc, gpF, gAF, why );
    bool const oka = okv && dirichlet_daeslv( 2, PC, AC, Jc, gpA, gAA, why );
    check( "(C) value, forward and adjoint solves" + ( okv && okf && oka? std::string(): " -- " + why ), okv && okf && oka );
    if( okv && okf && oka ){
      double const fdA = ( Jp - Jm )/( 2.*HC );
      std::cout << "    dJ/dp: fwd " << gpF << " adj " << gpA << " (exact 1)   dJ/dA: fwd " << gAF << " adj " << gAA << " (FD " << fdA << ")\n";
      check( "(C) dJ/dp = 1, forward and adjoint, to 1e-9",                 rel( gpF, 1. ) <= 1e-9 && rel( gpA, 1. ) <= 1e-9 );
      check( "(C) dJ/dA = central difference (J linear in A), forward and adjoint, to 1e-9", rel( gAF, fdA ) <= 1e-9 && rel( gAA, fdA ) <= 1e-9 );
    }
  }
  std::cout << "==========  DAESLV_molsens: " << ( nfail? std::to_string( nfail ) + " FAILED": std::string( "ALL PASS" ) ) << "  ==========\n";
  return nfail? 1: 0;
}
