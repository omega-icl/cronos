// ============================================================================
//  OCFE_SYMDIFF1.cpp  --  FFOCFESLV SYMBOLIC differentiation (SYMDIFF) and
//                         FIXED inputs (fix_input) on a time-dependent model
//                         with a distributed control, a distributed fixed
//                         input and a constant; monolithic and marching.
//
//  MODEL (1 state, t in [0,T], N_EL elements, LGR collocation)
//        dx/dt = -k x + g u(t) + d(t) ,   x(0) = x0
//    k     scalar input       -> map 1 (numerical gradient), flat
//    u(t)  distributed input, piecewise-constant per element (n_node=1)
//                             -> map 1, by GENERATOR per element
//    x0    scalar input       -> map 2, DAG variable
//    g     CONSTANT           -> map 2, DAG variable
//    d(t)  distributed disturbance, FIXED by fix_input: the MODEL supplies
//          its values; it is NOT mapped (FFOCFESLV exempts fixed inputs)
//    outputs  G0 = x(T),  G1 = x(T/2)   (T/2 is an element boundary)
//
//  ORACLE (independent of the collocation)
//    With u, d constant on each element the exact solution is the recurrence
//        x_{e+1} = c_e/k + ( x_e - c_e/k ) exp(-k h),   c_e = g u_e + d_e ,
//    and the 7 derivatives [k, u_0..u_3, x0, g] are central differences of it.
//
//  TESTS (per mode: monolithic, marching)
//    F1  init() writes d's FIXED values into its slots
//    F2  solve with garbage in d's slots == oracle (the model's values win)
//    V1  FFOCFESLV two-map value == oracle ; d left unmapped is accepted
//    N1  numeric route: dF/d(k,u) == oracle ; x0 and g columns ZERO (contract)
//    S1  SYMDIFF={x0,g}: map-2 input and CONSTANT symbolically == oracle ;
//        k and u columns ZERO (contract)
//    S2  SYMDIFF=all 7: full Jacobian == oracle ; and == numeric on k,u
//    F3  SHALLOW op: fix_input(d) again AFTER the op was built -> the SAME op
//        now returns the oracle at the new d (fixed values apply per call)
//  CROSS-MODE
//    XM  SYMDIFF Jacobian monolithic ~ marching
//
//  Build:
//    g++ -std=c++17 <suite flags> -DOCFE_OCFESLV_HEADER='"ocfeslv.hpp"' \
//        OCFE_SYMDIFF1.cpp -o OCFE_SYMDIFF1 <libs>
// ============================================================================

#include <array>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#ifndef OCFE_OCFESLV_HEADER
#define OCFE_OCFESLV_HEADER "ocfeslv.hpp"
#endif
#include OCFE_OCFESLV_HEADER
#include "ffocfe.hpp"

using namespace mc;

static double const T_end = 2.0;
static size_t const N_EL  = 4;
static size_t const N_ND  = 8;                               // discretisation error well below the oracle tolerances
static double const k_nom = 0.8, x0_nom = 1.0, g_nom = 1.2;
static std::vector<double> const u_nom { 0.3, 0.6, 0.2, 0.5 };
static std::vector<double> const d_fix { 0.1, -0.05, 0.2, 0.0 };
static std::vector<double> const d_new { 0.4, 0.4, 0.4, 0.4 };

static int g_pass = 0, g_fail = 0;
static void check_close( std::string const& name, double got, double want, double tol )
{
  bool const ok = std::fabs( got - want ) <= tol;
  ( ok ? g_pass : g_fail )++;
  std::cout << "  " << std::left << std::setw(52) << name << std::right << std::scientific << std::setprecision(3)
            << " |got-want|=" << std::fabs( got - want ) << " tol=" << tol << "  " << ( ok ? "PASS" : "FAIL" ) << "\n";
}
static void check_true( std::string const& name, bool ok )
{
  ( ok ? g_pass : g_fail )++;
  std::cout << "  " << std::left << std::setw(52) << name << "  " << ( ok ? "PASS" : "FAIL" ) << "\n";
}

// ---- oracle: exact element recurrence; p = [ k, u_0..u_{N-1}, x0, g ] -------------------------------------------
static std::array<double,2> oracle( std::vector<double> const& p, std::vector<double> const& d )
{
  double const k = p[0], x0 = p[N_EL+1], g = p[N_EL+2], h = T_end / N_EL;
  double x = x0;  std::array<double,2> out{ 0., 0. };
  for( size_t e = 0; e < N_EL; ++e ){
    double const c = g * p[1+e] + d[e];
    x = c/k + ( x - c/k ) * std::exp( -k*h );
    if( e+1 == N_EL/2 ) out[1] = x;
  }
  out[0] = x;
  return out;
}
static std::vector<std::array<double,2>> oracle_jac( std::vector<double> const& p, std::vector<double> const& d )
{
  std::vector<std::array<double,2>> J( p.size() );                     // J[col][fct]
  for( size_t i = 0; i < p.size(); ++i ){
    double const s = 1e-6 * std::max( 1., std::fabs( p[i] ) );
    std::vector<double> pp( p ), pm( p ); pp[i] += s; pm[i] -= s;
    auto const fp = oracle( pp, d ), fm = oracle( pm, d );
    for( size_t j = 0; j < 2; ++j ) J[i][j] = ( fp[j] - fm[j] ) / ( 2.*s );
  }
  return J;
}

struct Cell { bool ok = false; std::vector<double> Jsym; };           // Jsym[fct*7 + col]

static Cell run_mode( bool marching )
{
  Cell R;
  std::cout << "\n================ " << ( marching ? "MARCHING" : "MONOLITHIC" ) << " ================\n";
  FFGraph DAG;
  FFVar t  = DAG.add_var( "t" ),  x = DAG.add_var( "x(t)" );
  FFVar k  = DAG.add_var( "k" ),  u = DAG.add_var( "u(t)" ), d = DAG.add_var( "d(t)" );
  FFVar x0 = DAG.add_var( "x0" ), g = DAG.add_var( "g" );

  OCFESLV oc( &DAG );
  oc.add_domain( t, FFDom( 0., T_end, N_EL, FFDom::LGR, N_ND ) );
  oc.set_evolution_domain( t );
  oc.add_state( x, {t} );
  oc.add_input( k,  {} );
  oc.add_input( u,  {t}, FFDom::LGR, 1 );                              // piecewise constant per element
  oc.add_input( d,  {t}, FFDom::LGR, 1 );
  oc.add_input( x0, {} );
  oc.FFModel::set_constant( { g }, { g_nom } );
  oc.update_ref( x, []( OCFESLV::t_Coord const& ){ return x0_nom; } );

  FFPartial OpP;
  oc.add_equation( OpP( x, t ) + k*x - g*u - d, {t}, { FFDom::ALL - FFDom::LB }, OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( x - x0,                      {t}, { FFDom::LB },              OCFESLV::EqnOptions( OCFESLV::EqnRole::INITIAL,  0 ) );
  oc.add_output( x, {t}, { T_end } );                                  // G0
  oc.add_output( x, {t}, { 0.5*T_end } );                              // G1

  check_true( "fix_input( d, 4 values ) accepted before setup", oc.fix_input( d, d_fix ) );
  oc.options.SOLVE.MARCHING = marching;
  oc.options.SOLVE.MAX_ITER = 60;
  oc.options.SOLVE.RES_TOL  = 1.0e-12;
  oc.options.DISPLAY_LEVEL  = 0;
  if( !oc.setup() ){ std::cerr << "  setup() FAILED: " << OCFESLV::setup_status_str( oc.setup_status() ) << "\n"; return R; }
  check_true( "MODE is_marching()==requested", oc.is_marching() == marching );

  std::vector<double> p0{ k_nom };  for( double v : u_nom ) p0.push_back( v );  p0.push_back( x0_nom );  p0.push_back( g_nom );
  auto const Fo = oracle( p0, d_fix );
  auto const Jo = oracle_jac( p0, d_fix );
  double const tolV = 1e-7, tolJ = 1e-6;                               // spectral + oracle FD

  // ---- F1, F2: the model supplies d ------------------------------------------------------------------------------
  std::vector<double> xv, inp;
  oc.init( xv, inp, nullptr );
  auto const dval = oc.get_input_values( d, inp.data() );
  double ed = 0.; for( size_t e = 0; e < N_EL && e < dval.size(); ++e ) ed = std::max( ed, std::fabs( dval[e] - d_fix[e] ) );
  check_close( "F1 init() writes d's fixed values", ed, 0., 0. );
  oc.set_input_values( k,  { k_nom },  inp.data() );
  oc.set_input_values( u,  u_nom,      inp.data() );
  oc.set_input_values( x0, { x0_nom }, inp.data() );
  oc.set_input_values( d,  std::vector<double>( N_EL, 99. ), inp.data() );   // garbage: must be ignored
  std::vector<double> xs( xv );
  auto const rep = oc.solve( xs.data(), inp.data(), &g_nom );
  auto const Fs  = oc.val_functions();
  check_true ( "F2 solve converged", rep.converged );
  check_close( "F2 G0 with garbage in d's slots == oracle", Fs[0], Fo[0], tolV );
  check_close( "F2 G1 with garbage in d's slots == oracle", Fs[1], Fo[1], tolV );

  // ---- the reduced operation: map 1 {k, u}, map 2 {x0, g}; d unmapped ----------------------------------------------
  FFGraph rdag;
  FFVar pk = rdag.add_var( "pk" ), px0 = rdag.add_var( "px0" ), pg = rdag.add_var( "pg" );
  std::vector<FFVar> pu( N_EL );  for( size_t e = 0; e < N_EL; ++e ) pu[e] = rdag.add_var( "pu" + std::to_string( e ) );
  FFOCFESLV op;  std::vector<FFVar> F;
  try{
    F = op( { { k, std::vector<FFVar>{ pk } },
              { u, [&]( FFModel::DofIndex const& di ){ return pu.at( di.element.at( t ) ); } } },
            { { x0, std::vector<FFVar>{ px0 } }, { g, std::vector<FFVar>{ pg } } },
            &oc, FFOCFESLV::COPY, marching ? "SYM1_march" : "SYM1_mono" );
  } catch( std::exception& e ){ check_true( std::string( "V1 two-map embedding threw: " ) + e.what(), false ); return R; }
  check_true( "V1 two-map op built with d unmapped (fixed)", F.size() == 2 && op.n_input() == 7 );

  std::vector<FFVar>  vX { pk };  for( auto const& v : pu ) vX.push_back( v );  vX.push_back( px0 );  vX.push_back( pg );
  std::vector<double> vXv( p0 );                                        // same order as the oracle's p
  size_t const nX = vX.size();
  std::vector<double> Fv( 2 );  rdag.eval( F, Fv, vX, vXv );
  check_close( "V1 op value G0 == oracle", Fv[0], Fo[0], tolV );
  check_close( "V1 op value G1 == oracle", Fv[1], Fo[1], tolV );

  auto jac = [&]( std::vector<FFVar> const& sym ){
    FFOCFESLV::options.SYMDIFF = sym;
    auto const dF = rdag.FAD( F, vX );  std::vector<double> v( dF.size() );  rdag.eval( dF, v, vX, vXv );
    FFOCFESLV::options.SYMDIFF.clear();  return v; };                 // v[fct*nX + col]
  try{
    // ---- N1: numeric route (FFGradOCFESLV) --------------------------------------------------------------------
    auto const dN = jac( {} );
    double en = 0., zn = 0.;
    for( size_t j = 0; j < 2; ++j ){
      for( size_t c = 0; c <= N_EL; ++c ) en = std::max( en, std::fabs( dN[j*nX+c] - Jo[c][j] ) );
      for( size_t c = N_EL+1; c < nX; ++c ) zn = std::max( zn, std::fabs( dN[j*nX+c] ) );
    }
    check_close( "N1 numeric dG/d(k,u) == oracle", en, 0., tolJ );
    check_close( "N1 numeric dG/d(x0,g) ZERO by contract", zn, 0., 0. );
    // ---- S1: SYMDIFF on a map-2 input and a constant ------------------------------------------------------------
    auto const dS = jac( { px0, pg } );
    double es = 0., zs = 0.;
    for( size_t j = 0; j < 2; ++j ){
      for( size_t c = N_EL+1; c < nX; ++c ) es = std::max( es, std::fabs( dS[j*nX+c] - Jo[c][j] ) );
      for( size_t c = 0; c <= N_EL; ++c ) zs = std::max( zs, std::fabs( dS[j*nX+c] ) );
    }
    check_close( "S1 SYMDIFF={x0,g} dG/d(x0,g) == oracle", es, 0., tolJ );
    check_close( "S1 SYMDIFF={x0,g} dG/d(k,u) ZERO by contract", zs, 0., 0. );
    // ---- S2: SYMDIFF on everything --------------------------------------------------------------------------------
    auto const dA = jac( vX );
    double ea = 0., ex = 0.;
    for( size_t j = 0; j < 2; ++j )
      for( size_t c = 0; c < nX; ++c ){
        ea = std::max( ea, std::fabs( dA[j*nX+c] - Jo[c][j] ) );
        if( c <= N_EL ) ex = std::max( ex, std::fabs( dA[j*nX+c] - dN[j*nX+c] ) );
      }
    check_close( "S2 SYMDIFF=all: full Jacobian == oracle", ea, 0., tolJ );
    check_close( "S2 SYMDIFF == numeric on k,u (same discretisation)", ex, 0., 1e-9 );
    R.Jsym = dA;
  } catch( std::exception& e ){ check_true( std::string( "Jacobians threw: " ) + e.what(), false ); return R; }

  // ---- F3: a SHALLOW op sees fix_input changed AFTER it was built -------------------------------------------------
  {
    FFGraph sdag;
    FFVar sk = sdag.add_var( "sk" ), sx0 = sdag.add_var( "sx0" ), sg = sdag.add_var( "sg" );
    std::vector<FFVar> su( N_EL );  for( size_t e = 0; e < N_EL; ++e ) su[e] = sdag.add_var( "su" + std::to_string( e ) );
    FFOCFESLV sop;
    auto const FS = sop( { { k, std::vector<FFVar>{ sk } }, { u, su } }, { { x0, std::vector<FFVar>{ sx0 } }, { g, std::vector<FFVar>{ sg } } },
                         &oc, FFOCFESLV::SHALLOW );
    std::vector<FFVar> sX { sk };  for( auto const& v : su ) sX.push_back( v );  sX.push_back( sx0 );  sX.push_back( sg );
    std::vector<double> FvS( 2 );
    check_true( "F3 fix_input( d ) again after setup, on the shared solver", oc.fix_input( d, d_new ) );
    sdag.eval( FS, FvS, sX, vXv );
    auto const Fn = oracle( p0, d_new );
    check_close( "F3 SAME op now returns G0 at the new d == oracle", FvS[0], Fn[0], tolV );
    check_close( "F3 SAME op now returns G1 at the new d == oracle", FvS[1], Fn[1], tolV );
  }
  R.ok = true;
  return R;
}

int main()
{
  std::cout << "================================================================\n"
            << "  OCFE_SYMDIFF1 : FFOCFESLV SYMDIFF + fix_input on a time-dependent model\n"
            << "  dx/dt = -k x + g u(t) + d(t), x(0)=x0 ; u control, d FIXED, g constant\n"
            << "================================================================\n";
  Cell const mono  = run_mode( false );
  Cell const march = run_mode( true );

  std::cout << "\n---------------- cross-mode ----------------\n";
  if( mono.ok && march.ok && mono.Jsym.size() == march.Jsym.size() ){
    double e = 0.; for( size_t i = 0; i < mono.Jsym.size(); ++i ) e = std::max( e, std::fabs( mono.Jsym[i] - march.Jsym[i] ) );
    check_close( "XM SYMDIFF Jacobian monolithic ~ marching", e, 0., 1e-6 );
  }
  else check_true( "XM both modes completed", false );

  std::cout << "\n============================================================\n"
            << "  RESULT: " << g_pass << " passed, " << g_fail << " failed  -- " << ( g_fail == 0 ? "ALL PASS" : "SOME FAILED" ) << "\n"
            << "============================================================\n";
  return g_fail == 0 ? 0 : 1;
}
