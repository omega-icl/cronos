// DAESLV_molpsa -- MoL on an adsorption column: the PSA marching model of OCFE_PDE36c (which OCFE_PSA2/5 build on),
// discretised in z by OCFEMOL, exported to DAESLV_IDAS, against OCFESLV collocating t and z.
//
// Five fields: c1, c2 (advection-dispersion, Pe = 10, Danckwerts), q1, q2 (LDF towards a temperature-dependent
// competitive Langmuir isotherm -- an ODE at EVERY node, the column ends included), T (energy, with heats of
// adsorption and wall loss).  The mass balances carry mass on TWO fields (c_t + F q_t) and the energy balance on
// three: the first models in the gate where Mc is not one entry per row.
//
// The check is the gate's same-system check: OCFESLV's time collocation refined in ORDER must drive |OCFESLV - MoL|
// down spectrally, all five fields, at the final time.  Then the whole-solve cost of the two, at matched accuracy --
// the measurement the residual-kernel decision (MoL workplan, section 1) has been waiting for.  Timings print unless
// CRONOS_SNAPSHOT is set.
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include "daeslv_idas.hpp"
#include "ffocfemol.hpp"
#include "ocfeslv.hpp"

using namespace mc;
static int n_pass = 0, n_fail = 0;
static void check( std::string const& w, bool ok )
{ ok? ++n_pass: ++n_fail; std::cout << "    " << std::left << std::setw(66) << w << ( ok? "PASS": "FAIL" ) << "\n"; }
static void check_lt( std::string const& w, double got, double bound )
{
  bool const ok = ( got == got ) && got < bound;  ok? ++n_pass: ++n_fail;
  std::cout << "    " << std::left << std::setw(54) << w << std::right << std::scientific << std::setprecision(2)
            << std::setw(10) << got << " < " << std::setw(8) << bound << "  " << ( ok? "PASS": "FAIL" ) << std::defaultfloat << "\n";
}

// OCFE_PDE36c's constants, the isotherm affinity at its nominal value
static double const U_VEL = 1.0, D_ax = 0.1, F_ph = 0.5, k_ldf = 2.0;
static double const qs1 = 1.0, qs2 = 1.0, b01 = 3.0, b02 = 1.0, beta1 = 2.0, beta2 = 1.0, T0_ref = 1.0;
static double const Cp_e = 1.0, G_cv = 1.0, lam = 0.1, dH1 = 1.0, dH2 = 0.5, hw = 0.5, Tw = 1.0;
static double const c0_1 = 0.5, c0_2 = 0.5, tau_in = 0.15, T_END = 5.0;

struct Fields { FFVar t, z, c1, c2, q1, q2, T; };

template <typename OC>
static void declare_psa( OC& oc, FFGraph& G, Fields& F, size_t NEZ, size_t NNZ, size_t NELT, size_t NODT )
{
  F.t = G.add_var("t"); F.z = G.add_var("z");
  F.c1 = G.add_var("c1"); F.c2 = G.add_var("c2"); F.q1 = G.add_var("q1"); F.q2 = G.add_var("q2"); F.T = G.add_var("T");
  FFVar const &t = F.t, &z = F.z, &c1 = F.c1, &c2 = F.c2, &q1 = F.q1, &q2 = F.q2, &T = F.T;
  FFPartial OpP;
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
  oc.options.INTERFACE.IMPOSITION = OC::Options::IC_STRONG;
  oc.options.INTERFACE.SAT_SIGMA0 = 10.0;
  oc.options.DISPLAY_LEVEL = 0;
}

// the five fields at T_END, field-major (c1, c2, q1, q2, T), at the z coordinates @p zs
struct Result { std::vector<double> u; std::vector<double> zs; size_t n_states = 0, n_diff = 0; int index = -1;
                double setup_ms = 0., solve_ms = 0.; bool ok = false; };

static Result leg3( size_t NEZ, size_t NNZ, double rtol, double atol )
{
  Result R;
  auto const t0 = std::chrono::steady_clock::now();
  FFGraph G; OCFEMOL mol( &G ); Fields F;
  declare_psa( mol, G, F, NEZ, NNZ, 1, 2 );
  if( !mol.setup() ){ std::cout << "      leg3: OCFEMOL::setup failed\n"; return R; }
  R.n_states = mol.semi_discrete().n_states;  R.n_diff = mol.semi_discrete().n_diff_cols;
  FFGraph H; DAESLV_IDAS S; std::vector<FFVar> xc; FFVar tt; std::string err;
  if( !export_semidiscrete( mol, H, S, FFDom( 0., T_END, 5, FFDom::LGR, 3 ), xc, tt, &err ) ){ std::cout << "      leg3: " << err << "\n"; return R; }
  std::vector<double> x0( R.n_states );
  if( !mol.init_seed( 0., x0.data() ) ) return R;
  for( size_t i = 0; i < x0.size(); ++i ) S.update_ref( xc[i], x0[i] );
  FFEval OpE;
  for( FFVar const* f : { &F.c1, &F.c2, &F.q1, &F.q2, &F.T } )
    for( size_t c : mol.field_columns( *f ) ) S.add_output( OpE( xc[c], tt, T_END ) );
  for( auto const& nd : mol.field_nodes( F.c1 ) ) R.zs.push_back( nd.empty()? 0.: nd[0].second );
  S.options.SOLVE.RTOL = rtol;  S.options.SOLVE.ATOL = atol;  S.options.SOLVE.NMAX = 200000;
  S.options.SOLVE.VERBOSE = 0;  S.options.DISPLAY_LEVEL = 0;
  if( !S.setup() ){ std::cout << "      leg3: DAESLV::setup failed: " << S.extract_error() << "\n"; return R; }
  R.index = S.structural_index( 0, tt ).index;
  auto const t1 = std::chrono::steady_clock::now();
  std::vector<double> P( S.np(), 0. );
  if( S.solve( P ) != DAESLV_IDAS::NORMAL ){ std::cout << "      leg3: solve failed: " << S.extract_error() << "\n"; return R; }
  auto const t2 = std::chrono::steady_clock::now();
  R.u = S.val_function();
  R.ok = ( R.u.size() == 5*R.zs.size() );
  R.setup_ms = std::chrono::duration<double,std::milli>( t1 - t0 ).count();
  R.solve_ms = std::chrono::duration<double,std::milli>( t2 - t1 ).count();
  return R;
}

static Result leg1( size_t NEZ, size_t NNZ, size_t NELT, size_t NODT, std::vector<double> const& zs )
{
  Result R;
  auto const t0 = std::chrono::steady_clock::now();
  FFGraph G; OCFESLV oc( &G ); Fields F;
  declare_psa( oc, G, F, NEZ, NNZ, NELT, NODT );
  oc.options.SOLVE.MARCHING = false;  oc.options.SOLVE.RES_TOL = 1e-12;
  if( !oc.setup() ){ std::cout << "      leg1: setup failed\n"; return R; }
  R.n_states = oc.n_colloc_sta();
  std::vector<double> xv, inp;
  if( !oc.init( xv, inp, nullptr ) ){ std::cout << "      leg1: init failed\n"; return R; }
  auto const t1 = std::chrono::steady_clock::now();
  auto const rep = oc.solve( xv.data(), inp.empty()? nullptr: inp.data(), nullptr );
  auto const t2 = std::chrono::steady_clock::now();
  if( !rep.converged ){ std::cout << "      leg1: solve did not converge\n"; return R; }
  for( FFVar const* f : { &F.c1, &F.c2, &F.q1, &F.q2, &F.T } )
    for( double z : zs ){ std::map<FFVar,double,lt_FFVar> pt{ { F.t, T_END }, { F.z, z } }; R.u.push_back( oc.eval_solution( *f, pt ) ); }
  R.zs = zs;  R.ok = true;
  R.setup_ms = std::chrono::duration<double,std::milli>( t1 - t0 ).count();
  R.solve_ms = std::chrono::duration<double,std::milli>( t2 - t1 ).count();
  return R;
}

static double maxdiff( std::vector<double> const& a, std::vector<double> const& b )
{
  if( a.size() != b.size() || a.empty() ) return 1e30;
  double d = 0.;  for( size_t i = 0; i < a.size(); ++i ) d = std::max( d, std::fabs( a[i] - b[i] ) );
  return d;
}

int main()
{
  std::cout.setf( std::ios::unitbuf );
  bool const show_time = std::getenv( "CRONOS_SNAPSHOT" ) == nullptr;
  std::cout << "\n--- DAESLV_molpsa: the PSA column through MoL, against OCFESLV ----------------------------\n";
  size_t const NEZ = 6, NNZ = 5;
  std::cout << "\n  z: " << NEZ << " LGL elements x " << NNZ << " nodes;  t in [0, " << T_END << "]\n\n";

  auto const L3 = leg3( NEZ, NNZ, 1e-11, 1e-13 );
  check( "PSA leg 3: OCFEMOL set up, exported and integrated", L3.ok );
  if( !L3.ok ){ std::cout << "\n==========  DAESLV_molpsa: " << n_pass << " passed, " << n_fail << " failed -- FAIL  ==========\n\n"; return 1; }
  std::cout << "      semi-discrete: " << L3.n_states << " unknowns, " << L3.n_diff << " differential;  the exported DAE's index along t: "
            << L3.index << "\n";
  check( "PSA: the exported DAE has structural index <= 1", L3.index >= 0 && L3.index <= 1 );
  // a physical sanity check: by T = 5 the feed (0.5 each) has broken through, the bed is loaded, T is near ambient
  double cmax = 0., qmin = 1e30; for( size_t k = 0; k < L3.zs.size(); ++k ){ cmax = std::max( cmax, L3.u[k] ); qmin = std::min( qmin, L3.u[2*L3.zs.size()+k] ); }
  std::cout << "      at T = " << T_END << ": max c1 = " << std::setprecision(4) << cmax << ", min q1 = " << qmin << std::defaultfloat << "\n";
  check( "PSA: c1 in [0, feed], q1 > 0 at T (breakthrough, loaded bed)", cmax > 0.1 && cmax <= 0.5 + 1e-6 && qmin > 0. );

  std::cout << "\n  OCFESLV refined in time ORDER (10 LGR time elements), all five fields at T:\n";
  std::vector<double> d;  std::vector<Result> L1s;
  for( size_t nodt : { 3u, 4u, 5u, 6u } ){
    auto const L1 = leg1( NEZ, NNZ, 10, nodt, L3.zs );  L1s.push_back( L1 );
    d.push_back( L1.ok? maxdiff( L1.u, L3.u ): 1e30 );
    std::cout << "      NODT = " << nodt << "   (" << std::setw(5) << L1.n_states << " unknowns)   |OCFESLV - MoL| = "
              << std::scientific << std::setprecision(3) << d.back() << std::defaultfloat << "\n";
  }
  check_lt( "PSA: the two converge spectrally in OCFESLV's time order: >= 5x per order", 5.*d[3], d[2] );
  check_lt( "PSA: they agree to 1e-6 at order 6",                                         d[3], 1e-6 );

  if( show_time ){
    std::cout << "\n  whole-solve cost (setup | solve, ms):\n";
    std::cout << "      MoL (OCFEMOL + DAESLV, rtol 1e-11):  " << std::fixed << std::setprecision(1) << std::setw(8) << L3.setup_ms << " | " << std::setw(8) << L3.solve_ms << "\n";
    for( size_t i = 0; i < L1s.size(); ++i )
      std::cout << "      OCFESLV, NODT " << ( i + 3 ) << ":                     " << std::setw(8) << L1s[i].setup_ms << " | " << std::setw(8) << L1s[i].solve_ms << "\n";
    std::cout << std::defaultfloat;
  }

  std::cout << "\n==========  DAESLV_molpsa: " << n_pass << " passed, " << n_fail << " failed -- "
            << ( n_fail? "FAIL": "ALL PASS" ) << "  ==========\n\n";
  return n_fail? 1: 0;
}
