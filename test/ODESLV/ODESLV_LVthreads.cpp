// ODESLV_LVthreads.cpp -- repurposed from test0b.cpp (Lotka-Volterra, threaded solver copies).
// =============================================================================================
//   dx0/dt = p x0 (1 - x1),   dx1/dt = p x1 (x0 - 1),   t in [0,10]
// Everything is DECLARED through FFModel (no set_time / set_function).  Functions: x0 x1 at tf; the integral of x1
// over the WHOLE domain (OpI) -- test0b's second function (its q(T1) + q(tf) was a slip in the original); and, in
// the two-stage model, x0^2 SUMMED over the stage times T1 and tf (evaluations at different times, one function).
// The state is continuous across stages (test0b's transition x1 -> x1 - 0.5 at stage 1 is omitted -- see
// KNOWN_ISSUE_ASA_per_stage_IC.md).
//
// PART A -- ODESLVS, 2 stages (T1 = 5): the integral spans both stages.
//   NTH copies made by setup( LV ), each solving FORWARD and ADJOINT sensitivities at its own p in its own thread.
//     A1  threaded values and gradients BITWISE EQUAL to a serial rerun (thread safety of the copies)
//     A2  forward sensitivity == adjoint
//     A3  ORACLE: x0 x1 (tf) and int x1 dt of the 2-stage model == the single-stage model (independent of the stage split)
// PART B -- FFODESLV on the single-stage model, same two outputs.
//     B1  NTH FFODESLV ops, one DAG each (COPY), evaluated in threads == serial
//     B2  numeric (FFGradODESLV) == SYMDIFF (fdiff), and == central FD through the op, per thread's p
//
// Build: g++ -std=c++17 ... -pthread ODESLV_LVthreads.cpp

#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <thread>
#include <vector>
#include "ffode.hpp"

static int npass = 0, nfail = 0;
static void check( bool c, std::string const& w ){ std::cout << "  " << ( c? "PASS  ": "FAIL  " ) << w << "\n"; c? ++npass: ++nfail; }
static double rel( double a, double b ){ double s = std::fabs( b ); return std::fabs( a - b ) / ( s > 1e-12? s: 1. ); }

unsigned const NTH = 8;
double const t0 = 0., tf = 10., T1 = 5.;
double const pL = 2.95, dp = 0.1 / ( NTH - 1. );

//! @brief The declared Lotka-Volterra model on [t0,tf], 2 stages (@p staged) or 1; outputs x0 x1 (tf), int x1 dt,
//! and (staged) x0^2 (T1) + x0^2 (tf).
struct LV { mc::FFGraph DAG; mc::ODESLVS_CVODES IVP{ &DAG }; mc::FFVar t, x0, x1, p; };
static bool build( LV& M, bool staged )
{
  auto& G = M.DAG;  auto& I = M.IVP;
  M.t = G.add_var( "t" );  M.x0 = G.add_var( "x0(t)" );  M.x1 = G.add_var( "x1(t)" );  M.p = G.add_var( "p" );
  mc::FFPartial OpP;  mc::FFEval OpE;  mc::FFIntegral OpI;
  I.add_domain( M.t, mc::FFDom( t0, tf, staged? 2: 1, mc::FFDom::LGR, 4 ) );
  I.add_state( M.x0, {M.t} );  I.add_state( M.x1, {M.t} );
  I.add_input( M.p );
  I.set_evolution_domain( M.t );  I.update_ref( M.x0, 1.2 );  I.update_ref( M.x1, 1.1 );
  int const T_INT = mc::FFDom::ALL - mc::FFDom::LB;
  auto const io = mc::FFModel::EqnOptions( mc::FFModel::EqnRole::INTERIOR, 0 );
  auto const ii = mc::FFModel::EqnOptions( mc::FFModel::EqnRole::INITIAL,  0 );
  I.add_equation( OpP( M.x0, M.t ) - M.p * M.x0 * ( 1. - M.x1 ), {M.t}, {T_INT},         io );
  I.add_equation( OpP( M.x1, M.t ) - M.p * M.x1 * ( M.x0 - 1. ), {M.t}, {T_INT},         io );
  I.add_equation( M.x0 - 1.2,                                    {M.t}, {mc::FFDom::LB}, ii );
  I.add_equation( M.x1 - ( 1.1 + 0.01 * M.p ),                   {M.t}, {mc::FFDom::LB}, ii );
  I.add_output( OpE( M.x0 * M.x1, M.t, tf ) );
  I.add_output( OpI( std::vector<mc::FFVar>{ M.x1 }, { M.t } )[0] );   // int x1 dt over the whole domain
  if( staged ) I.add_output( OpE( sqr( M.x0 ), M.t, T1 ) + OpE( sqr( M.x0 ), M.t, tf ) );   // x0^2 at the stage times
  I.options.INTMETH = mc::BASE_CVODES::Options::MSBDF;  I.options.NLINSOL = mc::BASE_CVODES::Options::NEWTON;
  I.options.LINSOL  = mc::BASE_CVODES::Options::DENSE;  I.options.NMAX = 2000;  I.options.DISPLAY = 0;
  I.options.ATOL = I.options.ATOLS = I.options.ATOLB = 1e-10;
  I.options.RTOL = I.options.RTOLS = I.options.RTOLB = 1e-10;
  return I.setup();
}

struct Res { std::vector<double> F; std::vector<std::vector<double>> Gf, Ga; bool ok = false; };
static void solve_all( mc::ODESLVS_CVODES& ivp, std::vector<double> const& p, Res& r )
{
  r.ok = ivp.solve_fsens( p ) == mc::ODESLVS_CVODES::STATUS::NORMAL;
  r.F = ivp.val_function();  r.Gf = ivp.val_function_gradient();
  r.ok = r.ok && ivp.solve_asens( p ) == mc::ODESLVS_CVODES::STATUS::NORMAL;
  r.Ga = ivp.val_function_gradient();
}

int main()
{
  std::cout << std::scientific << std::setprecision(6)
            << "================================================================\n"
            << "  ODESLV_LVthreads: Lotka-Volterra, threaded ODESLVS copies and FFODESLV ops\n"
            << "================================================================\n";
  std::vector<double> pv( NTH );  for( unsigned i = 0; i < NTH; ++i ) pv[i] = pL + i*dp;
  std::vector<std::vector<double>> pvv( NTH );  for( unsigned i = 0; i < NTH; ++i ) pvv[i] = { pv[i] };

  std::cout << "\n--- PART A: ODESLVS, declared, 2 stages: an integral spanning both, x0^2 summed at the stage times\n";
  LV M;  check( build( M, true ) && M.IVP.nf() == 3, "A0 declared two-stage model sets up, 3 functions" );
  { std::vector<mc::ODESLVS_CVODES> civ( NTH );  std::vector<Res> rt( NTH ), rs( NTH );  std::vector<std::thread> th( NTH );
    bool setup_ok = true;
    for( unsigned i = 0; i < NTH; ++i ) setup_ok = civ[i].setup( M.IVP ) && setup_ok;
    check( setup_ok, "A0 " + std::to_string( NTH ) + " copies set up by setup( LV )" );
    for( unsigned i = 0; i < NTH; ++i ) th[i] = std::thread( solve_all, std::ref( civ[i] ), std::cref( pvv[i] ), std::ref( rt[i] ) );
    for( auto& t : th ) t.join();
    bool bit = true, allok = true;
    for( unsigned i = 0; i < NTH; ++i ){
      mc::ODESLVS_CVODES fresh;  fresh.setup( M.IVP );  solve_all( fresh, pvv[i], rs[i] );
      allok = allok && rt[i].ok && rs[i].ok;
      bit = bit && rt[i].F == rs[i].F && rt[i].Gf == rs[i].Gf && rt[i].Ga == rs[i].Ga;
    }
    check( allok, "A0 every threaded and serial solve returned NORMAL" );
    check( bit, "A1 threaded F, forward and adjoint gradients BITWISE EQUAL to a serial rerun" );
    double wfa = 0.;
    for( unsigned i = 0; i < NTH; ++i ) for( size_t k = 0; k < rt[i].F.size(); ++k ) wfa = std::max( wfa, rel( rt[i].Ga[0][k], rt[i].Gf[0][k] ) );
    { std::ostringstream o; o << "A2 forward == adjoint, 2 stages (worst rel " << wfa << ")"; check( allok && wfa < 1e-6, o.str() ); }
  }

  std::cout << "\n--- PART B: FFODESLV on the declared single-stage model\n";
  LV S;  check( build( S, false ), "B0 declared single-stage model sets up" );
  { // ORACLE for Part A: x0 x1 (tf) and int x1 dt do not depend on the stage split -- the single-stage model (where the
    // integral needs no summing over stages) must give the same values and gradients as the 2-stage model
    Res r2, r1;  mc::ODESLVS_CVODES c2, c1;  c2.setup( M.IVP );  c1.setup( S.IVP );
    solve_all( c2, pvv[0], r2 );  solve_all( c1, pvv[0], r1 );
    double wv = 0., wg = 0.;
    for( size_t k = 0; k < 2; ++k ){ wv = std::max( wv, rel( r2.F[k], r1.F[k] ) ); wg = std::max( wg, rel( r2.Gf[0][k], r1.Gf[0][k] ) ); }
    std::ostringstream o; o << "A3 ORACLE: 2-stage x0x1(tf), int x1 == single-stage (value " << wv << ", gradient " << wg << ")";
    check( r1.ok && r2.ok && wv < 1e-7 && wg < 1e-6, o.str() ); }
  struct Op { mc::FFGraph NLP; mc::FFVar pp; std::vector<mc::FFVar> F; mc::FFODESLV op; };
  std::vector<Op> ops( NTH );
  for( unsigned i = 0; i < NTH; ++i ){
    ops[i].pp = ops[i].NLP.add_var( "pp" );
    ops[i].F  = ops[i].op( { { S.p, std::vector<mc::FFVar>{ ops[i].pp } } }, &S.IVP, mc::FFODESLV::COPY );
  }
  auto value_and_numeric = []( Op& o, double p, std::vector<double>& F, std::vector<double>& dF ){
    F.assign( 2, 0. );  o.NLP.eval( o.F, F, std::vector<mc::FFVar>{ o.pp }, std::vector<double>{ p } );
    auto const d = o.NLP.FAD( o.F, std::vector<mc::FFVar>{ o.pp } );
    dF.assign( d.size(), 0. );  o.NLP.eval( d, dF, std::vector<mc::FFVar>{ o.pp }, std::vector<double>{ p } ); };
  mc::FFODESLV::options.SYMDIFF.clear();                                // numeric route in the threads
  std::vector<std::vector<double>> Ft( NTH ), dFt( NTH ), Fs( NTH ), dFs( NTH );
  { std::vector<std::thread> th( NTH );
    for( unsigned i = 0; i < NTH; ++i ) th[i] = std::thread( value_and_numeric, std::ref( ops[i] ), pv[i], std::ref( Ft[i] ), std::ref( dFt[i] ) );
    for( auto& t : th ) t.join(); }
  bool bit = true;
  for( unsigned i = 0; i < NTH; ++i ){ value_and_numeric( ops[i], pv[i], Fs[i], dFs[i] ); bit = bit && Ft[i] == Fs[i] && dFt[i] == dFs[i]; }
  check( bit, "B1 " + std::to_string( NTH ) + " ops evaluated in threads: value and numeric derivative BITWISE EQUAL to serial" );
  double wsn = 0., wfd = 0.;
  for( unsigned i = 0; i < NTH; ++i ){
    mc::FFODESLV::options.SYMDIFF = { ops[i].pp };
    auto const d = ops[i].NLP.FAD( ops[i].F, std::vector<mc::FFVar>{ ops[i].pp } );
    std::vector<double> dS( d.size() );  ops[i].NLP.eval( d, dS, std::vector<mc::FFVar>{ ops[i].pp }, std::vector<double>{ pv[i] } );
    mc::FFODESLV::options.SYMDIFF.clear();
    double const h = 1e-6 * pv[i];  std::vector<double> Fp( 2 ), Fm( 2 );
    ops[i].NLP.eval( ops[i].F, Fp, std::vector<mc::FFVar>{ ops[i].pp }, std::vector<double>{ pv[i] + h } );
    ops[i].NLP.eval( ops[i].F, Fm, std::vector<mc::FFVar>{ ops[i].pp }, std::vector<double>{ pv[i] - h } );
    for( size_t k = 0; k < 2; ++k ){ wsn = std::max( wsn, rel( dS[k], dFs[i][k] ) ); wfd = std::max( wfd, rel( dS[k], ( Fp[k] - Fm[k] ) / ( 2.*h ) ) ); }
  }
  { std::ostringstream o; o << "B2 SYMDIFF (fdiff) == numeric (FFGradODESLV), every p (worst rel " << wsn << ")"; check( wsn < 1e-6, o.str() ); }
  { std::ostringstream o; o << "B2 SYMDIFF == central finite differences through the op (worst rel " << wfd << ")"; check( wfd < 1e-5, o.str() ); }

  std::cout << "\n  ODESLV_LVthreads: " << npass << " passed, " << nfail << " failed -- " << ( nfail? "FAILURES": "ALL PASS" ) << "\n";
  return nfail? 1: 0;
}
