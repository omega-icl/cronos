// OCFE_PDE41.cpp  ---  a SPATIAL reduction nested inside a nonlinearity of an output, on OCFESLV (9 Nov 2026)
// ===========================================================================
// The collocation evaluates an output element by element and sums at the root -- exact for anything linear in a
// reduction, wrong for a nonlinearity of one: (int u dz)^2 came out as sum_e (int_e u dz)^2, a THIRD of the value on
// three equal elements, for every output a user would write as a ratio of spatial integrals (a bed loading fraction).
// OCFESLV::_lift_nested_reductions (FFModel::_on_declaration) now lifts each such reduction into a hidden algebraic
// state J over the directions it does not consume, pinned by J - R = 0 at every node (linear: exact), and the output
// reads J -- mirroring OCFEMOL's spatial parts under the method of lines and BASE_DE's lift of nonlinear outputs.
//   (A) the zero-flux rod with u0 = 1 and source p e^-t: u = 1 + p (1 - e^-t) uniform, J = int u dz = u exactly on
//       any mesh -- int J^2 dt, J(T)^2 and the ratio J / int 2u dz against closed forms; int J dt (linear) untouched
//       and the state count says what was lifted;
//   (C) a constant or a parameter inside a spatial integral at an instant (int (1 + u) dz, int (u + p) dz at T): the
//       row section 4a assembles element by element counted its field-free part once -- g was evaluated with the
//       integrated FIELDS zeroed, which left int 1 dz standing on one element (1 - 1/n_el short); g is now the residual
//       with its integral NODES replaced by 0 (9 Nov 2026) -- and a ratio with a constant inside (lifted);
//   (D) a reduction consuming EVERY direction of its fields -- the steady BVP -u'' = 1, u(0) = u(1) = 0, int u dz = 1/12
//       exactly at degree >= 2: lifted into a LUMPED hidden state pinned by a domain-less row (10 Nov 2026; refused
//       before), (int u dz)^2 and int u dz / int 2u dz to round-off;
//   (B) the DAESLV_molsens PSA column: the bed loading fraction int q1 dz / int (q1 + q2) dz at the end, OCFESLV
//       (lifted) against DAESLV under the method of lines (composed), one declaration.
#include <iostream>
#include <iomanip>
#include <cmath>
#include "ocfeslv.hpp"
#include "daeslv_idas.hpp"
using namespace mc;

static int nfail = 0;
static void check( std::string const& w, bool ok ){ std::cout << "  " << std::left << std::setw(96) << w << ( ok? "PASS": "FAIL" ) << "\n"; if( !ok ) ++nfail; }
static double rel( double a, double b ){ return std::fabs( a - b ) / std::max( 1e-300, std::fabs( b ) ); }

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


int main()
{
  std::cout << std::scientific << std::setprecision(4);
  std::cout << "(A) the rod: nested reductions against closed forms\n";
  {
    double const P = 0.7, TEND = 1., e1 = std::exp( -TEND ), e2 = std::exp( -2.*TEND );
    double const Be = (1+P)*(1+P) - 2*P*(1+P)*(1-e1) + P*P*(1-e2)/2., JT = 1.+P*(1.-e1), Ae = 1.+P*e1;
    FFGraph G; OCFESLV S( &G ); FFPartial OpP; FFIntegral OpI; FFEval OpE;
    FFVar t = G.add_var("t"), z = G.add_var("z"), u = G.add_var("u"), p = G.add_var("p");
    S.add_domain( t, FFDom( 0., TEND, 4, FFDom::LGR, 6 ) ); S.add_domain( z, FFDom( 0., 1., 3, FFDom::LGL, 5 ) );
    S.set_evolution_domain( t ); S.add_state( u, {t,z} ); S.update_ref( u, 1. ); S.add_input( p, {} );
    int const TI = FFDom::ALL - FFDom::LB, ZI = FFDom::ALL - FFDom::LB - FFDom::UB;
    OCFESLV::EqnOptions io( OCFESLV::EqnRole::INTERIOR, 0 ), ib( OCFESLV::EqnRole::BOUNDARY, 0 ), ii( OCFESLV::EqnRole::INITIAL, 0 );
    S.add_equation( OpP(u,t) - OpP( 0.1*OpP(u,z), z ) - p*exp( -t ), {t,z}, {TI, ZI}, io );
    S.add_equation( OpP(u,z), {t,z}, {TI, FFDom::LB}, ib ); S.add_equation( OpP(u,z), {t,z}, {TI, FFDom::UB}, ib );
    S.add_equation( u - 1., {t,z}, {FFDom::LB, FFDom::ALL}, ii );
    FFVar J = OpI( u, z );
    S.add_output( OpI( J*J, t ) );  S.add_output( OpE( J*J, t, TEND ) );  S.add_output( OpE( J / OpI( 2.*u, z ), t, TEND ) );  S.add_output( OpI( J, t ) );
    S.options.INTERFACE.IMPOSITION = OCFESLV::Options::IC_STRONG; S.options.INTERFACE.SAT_SIGMA0 = 10.;
    S.options.SOLVE.MARCHING = false; S.options.SOLVE.RES_TOL = 1e-12; S.options.DISPLAY_LEVEL = 0;
    bool const su = S.setup();
    size_t nhid = 0;  for( auto const& kv : S.var_state() ) if( kv.first.name().rfind( "nest_", 0 ) == 0 ) ++nhid;   // the hidden states, by name
    bool const ok = su && S.solve( std::vector<FFModel::InputVal>{ { p, { P } } } ) == OCFESLV::NORMAL;
    double const F0 = ok? S.val_function()[0]: NAN, F1 = ok? S.val_function()[1]: NAN, F2 = ok? S.val_function()[2]: NAN, F3 = ok? S.val_function()[3]: NAN;
    std::cout << "    hidden states " << nhid << ";  int J^2 dt " << F0 << " (" << Be << ")  J(T)^2 " << F1 << " (" << JT*JT << ")  ratio " << F2 << " (0.5)  int J dt " << F3 << " (" << Ae << ")\n";
    check( "(A) sets up and solves; the two distinct reductions lifted once each (J = int u dz, int 2u dz)", ok && nhid == 2 );
    check( "(A) int J^2 dt and J(T)^2 = closed forms to 1e-8 (a third of them before)", ok && rel( F0, Be ) < 1e-8 && rel( F1, JT*JT ) < 1e-8 );
    check( "(A) the ratio of two spatial integrals = 1/2 to 1e-10", ok && rel( F2, 0.5 ) < 1e-10 );
    check( "(A) the linear int J dt = closed form to 1e-8 (the element-wise assembly, exact as before)", ok && rel( F3, Ae ) < 1e-8 );
  }
  std::cout << "(C) constants and parameters inside a spatial integral at an instant, 3 elements\n";
  {
    double const P = 0.7, TEND = 1., JT = 1.+P*(1.-std::exp(-TEND));
    FFGraph G; OCFESLV S( &G ); FFPartial OpP; FFIntegral OpI; FFEval OpE;
    FFVar t = G.add_var("t"), z = G.add_var("z"), u = G.add_var("u"), p = G.add_var("p");
    S.add_domain( t, FFDom( 0., TEND, 4, FFDom::LGR, 6 ) ); S.add_domain( z, FFDom( 0., 1., 3, FFDom::LGL, 5 ) );
    S.set_evolution_domain( t ); S.add_state( u, {t,z} ); S.update_ref( u, 1. ); S.add_input( p, {} );
    int const TI = FFDom::ALL - FFDom::LB, ZI = FFDom::ALL - FFDom::LB - FFDom::UB;
    OCFESLV::EqnOptions io( OCFESLV::EqnRole::INTERIOR, 0 ), ib( OCFESLV::EqnRole::BOUNDARY, 0 ), ii( OCFESLV::EqnRole::INITIAL, 0 );
    S.add_equation( OpP(u,t) - OpP( 0.1*OpP(u,z), z ) - p*exp( -t ), {t,z}, {TI, ZI}, io );
    S.add_equation( OpP(u,z), {t,z}, {TI, FFDom::LB}, ib ); S.add_equation( OpP(u,z), {t,z}, {TI, FFDom::UB}, ib );
    S.add_equation( u - 1., {t,z}, {FFDom::LB, FFDom::ALL}, ii );
    S.add_output( OpE( OpI( 1. + u, z ), t, TEND ) );  S.add_output( OpE( OpI( u + p, z ), t, TEND ) );
    S.add_output( OpE( OpI( u, z ) / OpI( 1. + u, z ), t, TEND ) );
    S.options.INTERFACE.IMPOSITION = OCFESLV::Options::IC_STRONG; S.options.SOLVE.MARCHING = false; S.options.SOLVE.RES_TOL = 1e-12; S.options.DISPLAY_LEVEL = 0;
    bool const ok = S.setup() && S.solve( std::vector<FFModel::InputVal>{ { p, { P } } } ) == OCFESLV::NORMAL;
    double const F0 = ok? S.val_function()[0]: NAN, F1 = ok? S.val_function()[1]: NAN, F2 = ok? S.val_function()[2]: NAN;
    std::cout << "    int (1+u) dz " << F0 << " (" << 1.+JT << ")  int (u+p) dz " << F1 << " (" << JT+P << ")  ratio " << F2 << " (" << JT/(1.+JT) << ")\n";
    check( "(C) int (1 + u) dz and int (u + p) dz at T = closed forms to 1e-8 (2/3 of the constant missing before)", ok && rel( F0, 1.+JT ) < 1e-8 && rel( F1, JT+P ) < 1e-8 );
    check( "(C) int u dz / int (1 + u) dz at T, lifted = closed form to 1e-8", ok && rel( F2, JT/(1.+JT) ) < 1e-8 );
  }
  std::cout << "(D) a reduction consuming every direction: the steady BVP -u'' = 1, lumped hidden states\n";
  {
    FFGraph G; OCFESLV S( &G ); FFPartial OpP; FFIntegral OpI;
    FFVar z = G.add_var("z"), u = G.add_var("u");
    S.add_domain( z, FFDom( 0., 1., 3, FFDom::LGL, 5 ) );  S.add_state( u, {z} );  S.update_ref( u, 0.1 );
    OCFESLV::EqnOptions io( OCFESLV::EqnRole::INTERIOR, 0 ), ib( OCFESLV::EqnRole::BOUNDARY, 0 );
    S.add_equation( -OpP( OpP(u,z), z ) - 1., {z}, { FFDom::ALL - FFDom::LB - FFDom::UB }, io );
    S.add_equation( u, {z}, { FFDom::LB }, ib );  S.add_equation( u, {z}, { FFDom::UB }, ib );
    S.add_output( OpI( u, z ) * OpI( u, z ) );  S.add_output( OpI( u, z ) / OpI( 2.*u, z ) );  S.add_output( OpI( u, z ) );
    S.options.INTERFACE.IMPOSITION = OCFESLV::Options::IC_STRONG;  S.options.SOLVE.RES_TOL = 1e-12;  S.options.DISPLAY_LEVEL = 0;
    bool const su = S.setup();
    bool const ok = su && S.solve( std::vector<FFModel::InputVal>{} ) == OCFESLV::NORMAL;
    size_t nlump = 0;  for( auto const& kv : S.var_state() ) if( kv.first.name().rfind( "nest_", 0 ) == 0 && kv.second.empty() ) ++nlump;
    double const F0 = ok? S.val_function()[0]: NAN, F1 = ok? S.val_function()[1]: NAN, F2 = ok? S.val_function()[2]: NAN;
    std::cout << "    lumped hidden states " << nlump << ";  (int u)^2 " << F0 << " (" << 1./144. << ")  ratio " << F1 << " (0.5)  int u dz " << F2 << " (" << 1./12. << ")\n";
    check( "(D) sets up and solves; the two reductions lifted into LUMPED hidden states (refused before)", ok && nlump == 2 );
    check( "(D) (int u dz)^2 = 1/144 and the ratio = 1/2 to 1e-12; the linear int u dz = 1/12 untouched", ok && rel( F0, 1./144. ) < 1e-12 && rel( F1, 0.5 ) < 1e-12 && rel( F2, 1./12. ) < 1e-12 );
  }
  std::cout << "(B) the PSA column: the bed loading fraction int q1 dz / int (q1 + q2) dz at the end, OCFESLV vs DAESLV (MoL)\n";
  {
    double fo = NAN, fd = NAN;  std::string why;
    { FFGraph G; OCFESLV oc( &G ); Fields F;  declare_psa( oc, G, F, 3, 5, 10, 6 );
      FFIntegral OpI;  oc.add_output( OpI( F.q1, F.z ) / OpI( F.q1 + F.q2, F.z ), {F.t}, {T_END} );
      colloc_options( oc );  oc.options.SOLVE.RES_TOL = 1e-10;  oc.options.SOLVE.MARCHING = false;  oc.options.DISPLAY_LEVEL = 0;
      if( oc.setup() && oc.solve( std::vector<FFModel::InputVal>{ { F.b01, { 3.0 } } } ) == OCFESLV::NORMAL ) fo = oc.val_function().back(); else why += " OCFESLV failed;"; }
    { FFGraph G; DAESLV_IDAS S( &G ); Fields F;  declare_psa( S, G, F, 3, 5, 5, 3 );
      FFIntegral OpI;  S.add_output( OpI( F.q1, F.z ) / OpI( F.q1 + F.q2, F.z ), {F.t}, {T_END} );
      colloc_options( S );  S.options.SOLVE.RTOL = 1e-10;  S.options.SOLVE.ATOL = 1e-12;  S.options.SOLVE.VERBOSE = 0;  S.options.DISPLAY_LEVEL = 0;
      if( S.setup() && S.solve( std::vector<FFModel::InputVal>{ { F.b01, { 3.0 } } } ) == DAESLV_IDAS::NORMAL ) fd = S.val_function().back(); else why += " DAESLV failed: " + S.extract_error(); }
    std::cout << "    loading fraction: OCFESLV " << fo << "  DAESLV " << fd << "  rel gap " << rel( fo, fd ) << why << "\n";
    check( "(B) both solve the column with the ratio as an output", std::isfinite( fo ) && std::isfinite( fd ) );
    check( "(B) OCFESLV (lifted) = DAESLV under MoL (composed) to 1e-4 -- the reference's time discretisation", std::isfinite( fo ) && std::isfinite( fd ) && rel( fo, fd ) < 1e-4 );
  }
  std::cout << "==========  OCFE_PDE41: " << ( nfail? std::to_string( nfail ) + " FAILED": "ALL PASS" ) << "  ==========\n";
  return nfail? 1: 0;
}
