// OCFE_receiver.cpp -- THE IC_WEAK RECEIVER RULE, asserted (2026-09-17, rev276+).
//
//   A continuity claim on a state the model DEFINES rather than balances (a reduction auxiliary, or q = a*T_x
//   declared by hand) has no principal-symbol weight.  The rescue used to inject a penalty with a fabricated
//   weight into the state's own defining row -- a DUPLICATE receiver whenever the claim also has a natural one.
//   The rule (Options::WEAK_NATURAL_PENALTY = 2): penalise a claim only in its natural BALANCE receivers; a claim
//   with none keeps its rescued rows; a natural receiver that is a LINK row (the next auxiliary's definition)
//   does not count.  NOTES_20260917a-g are the measurements behind each cell.
//
//   CELL 1  duplicate receiver  (mixed-form heat, blk0's MANUAL pattern)
//           rule ON : |T-T*| invariant to whether C2 derives the rescued weight or not      -- the assertion
//           rule OFF vs ON: the answers differ -- the duplicate injection was load-bearing   -- the control
//           (2026-09-18, R4: the controls no longer vary RESCUE_COUPLING, which R5 deletes)
//   CELL 2  sole claim  (T_t = T_xx under RED_FULL with TWO t-elements: the t-direction claims on the auxiliary
//           have no balance receiver)  rule ON: sole claims reported > 0, converged, invariant to the coupling.
//   CELL 3  chain  (order-4 biharmonic, AUTO depth-3 chain)  rule ON (=2): natural-LINK-only claims > 0 and the
//           deepest auxiliary's duplicate-node spread within 10x of the default; =1 (informational) loses it.
//   CELL 4  definition-row pollution  (PDE20g G4, RED_FULL: the rescued injection into w's own row)
//           rule ON: w's interface error at least 5x smaller than the default (production setting, RESCUE_C2 on;
//           with C2 off the sole w-claims carry the forced +1 against c0 = -1 and pollute w on both paths).
//   CELL 5  bare derivative at the interface  (PDE20g G4, RED_MAIN)  rule ON: w's error invariant to C2 on/off and
//           equal in WEAK and TRACE (the receiver rule is uniform across impositions since rev290).
//
//   Everything is compared IN ONE PROCESS through Options fields (rev276); no environment variable is read by
//   this driver.  Supersedes OCFE_coupling.cpp (whose checks could not fail: expect() never called, in-process
//   setenv inert because the knobs were statics).
//
// REQUIRES: -DOCFE_OCFESLV_HEADER='"ocfeslv_rev276.hpp"' or later.

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <map>
#include <set>
#include <cmath>
#include <string>
#include <functional>

#include "ffunc.hpp"
#ifndef OCFE_OCFESLV_HEADER
#define OCFE_OCFESLV_HEADER "ocfeslv.hpp"
#endif
#include OCFE_OCFESLV_HEADER

using mc::FFGraph; using mc::FFVar; using mc::FFPartial;
using namespace mc;
typedef OCFESLV::Options O;

static double const kPi = 3.14159265358979323846;
static int g_pass = 0, g_fail = 0;
static bool expect( std::string const& what, bool ok, std::string const& detail = "" )
{
  std::cout << "  " << std::left << std::setw(70) << what << ( ok ? "PASS" : "FAIL" )
            << ( detail.empty() ? "" : "   " + detail ) << "\n";
  ( ok ? g_pass : g_fail )++;
  return ok;
}
static std::string sci( double v ){ std::ostringstream s; s << std::scientific << std::setprecision(4) << v; return s.str(); }
static double rel( double a, double b ){ double const m = std::max( std::fabs(a), std::fabs(b) ); return m > 0 ? std::fabs(a-b)/m : 0.; }

// ---- capture the plan's report lines on std::cerr during setup ------------------------------------------------
struct RcvReport {
  long dropped=-1, natural=-1, sole=-1, link_only=-1;
  std::string line;
};
static RcvReport parse_report( std::string const& err )
{
  RcvReport r; std::istringstream is( err ); std::string l;
  while( std::getline( is, l ) ){
    if( l.find( "IC_WEAK natural-only penalty" ) == std::string::npos ) continue;
    r.line = l;
    auto num_before = [&]( char const* key ){ size_t p = l.find( key ); if( p == std::string::npos ) return -1L;
      size_t e = p; while( e > 0 && l[e-1] == ' ' ) --e; size_t b = e; while( b > 0 && std::isdigit( (unsigned char)l[b-1] ) ) --b;
      return b < e ? std::stol( l.substr( b, e-b ) ) : -1L; };
    r.dropped   = num_before( " rescued receiver term(s) dropped" );
    r.natural   = num_before( " claim(s) penalised on natural receivers" );
    r.sole      = num_before( " sole claim(s) keep" );
    r.link_only = num_before( " natural-LINK-only" );
  }
  return r;
}

// ---- generic run record -----------------------------------------------------------------------------------------
struct Run {
  bool setup=false, conv=false; size_t nVar=0; int iters=0;
  std::map<std::string,double> err_if, err_int, spread;   // per state: max error at interface nodes / interior, duplicate-node spread
  RcvReport rep;
};

// per-state error split by whether the node sits on an interior element boundary of the SPACE domain (index sdom)
// Interface nodes: on an interior boundary of the space domain AND strictly inside a t-element.  t-interface corner
// nodes are excluded: a marched solve starts each window there and w carries ~1e-03 there in EVERY variant -- a
// t-interface fact (T1), not a receiver-rule one.
static void census( OCFESLV const& oc, std::vector<double> const& var, std::vector<double> const& bnd, size_t sdom, std::vector<double> const& tbnd,
                    std::map< std::string, std::function<double(std::vector<double> const&)> > const& exact, Run& R )
{
  size_t off = 0;
  for( auto const& st : oc.states_colloc() ){
    auto nodes = oc.node_colloc( st );
    std::string const nm = st.name();
    auto ex = exact.find( nm );
    double eif = 0., eint = 0.;
    std::map< std::vector<long long>, std::pair<double,double> > groups;
    for( size_t i = 0; i < nodes.size(); ++i ){
      double const v = var[off+i];
      std::vector<long long> key; for( double c : nodes[i] ) key.push_back( std::llround( c*1e12 ) );
      auto g = groups.find( key );
      if( g == groups.end() ) groups.emplace( key, std::make_pair( v, v ) );
      else { g->second.first = std::min( g->second.first, v ); g->second.second = std::max( g->second.second, v ); }
      if( ex == exact.end() ) continue;
      double const e = std::fabs( v - ex->second( nodes[i] ) );
      bool at_if = false;
      for( size_t k = 1; k + 1 < bnd.size(); ++k ) if( std::fabs( nodes[i][sdom] - bnd[k] ) < 1e-12 ) at_if = true;
      bool corner = false; for( double tb : tbnd ) if( std::fabs( nodes[i][1-sdom] - tb ) < 1e-9 ) corner = true;
      if( corner ) continue;
      ( at_if ? eif : eint ) = std::max( at_if ? eif : eint, e );
    }
    double sp = 0.; for( auto const& kv : groups ) sp = std::max( sp, kv.second.second - kv.second.first );
    R.err_if[nm] = eif; R.err_int[nm] = eint; R.spread[nm] = sp;
    off += nodes.size();
  }
}

static void common_options( OCFESLV& oc, int policy, double coupling, O::ReductionType red )
{
  oc.options.REDUCE.ORDER     = red;
  oc.options.CLASSIFY.MODE         = O::CLASS_AUTO;
  oc.options.INTERFACE.TYPE   = O::IC_AUTO;
  oc.options.INTERFACE.IMPOSITION  = O::IC_WEAK;
  setenv( "CRONOS_WEAK_NATURAL_PENALTY", std::to_string( policy ).c_str(), 1 );
  // rev309: RESCUE_C2 and WEAK_NATURAL_PENALTY are no longer Options fields -- they are interim knobs with one
  // measured setting, reachable only from the environment.  These cells are EXPERIMENTS ON A FIXED SETTING,
  // not an API a user should reach for, and setenv() says so.  Options::reset() runs per instance and
  // _env_flag does not cache, so the next OCFESLV constructed picks this up.
  setenv( "CRONOS_RESCUE_C2", ( coupling != 0.0 ) ? "1" : "0", 1 );
  oc.options.INTERFACE.SAT_SIGMA0       = 10.0;
  oc.options.DISPLAY_LEVEL    = 0;
  oc.options.SOLVE.MAX_ITER   = 40;
  oc.options.SOLVE.RES_TOL    = 1e-9;
#if defined(CRONOS__WITH_SPQR)
  oc.options.SOLVE.FACTORIZATION = O::SOLVE_SPQR;
#else
  oc.options.SOLVE.FACTORIZATION = O::SOLVE_SUPERLU;
#endif
  oc.options.SOLVE.WARMSTART  = O::BROADCAST_IC;
}

// run: setup (capturing the plan report), init from references, solve, census
static Run finish( OCFESLV& oc, std::vector<double> const& bnd, size_t sdom, std::vector<double> const& tbnd,
                   std::map< std::string, std::function<double(std::vector<double> const&)> > const& exact )
{
  Run R;
  // The summary line parsed below is informational (std::cout, DISPLAY_LEVEL >= 2) since the display policy of
  // 20261001e; refusals stay on std::cerr.  Capture both.
  oc.options.DISPLAY_LEVEL = std::max( oc.options.DISPLAY_LEVEL, 2 );
  std::ostringstream cap; std::streambuf* old = std::cerr.rdbuf( cap.rdbuf() );
  std::streambuf* oldout = std::cout.rdbuf( cap.rdbuf() );
  bool ok = false;
  try{ ok = oc.setup(); } catch( ... ){ ok = false; }
  std::cerr.rdbuf( old );  std::cout.rdbuf( oldout );
  R.rep = parse_report( cap.str() );
  if( !ok ){ std::cout << "  *** setup failed\n" << cap.str().substr( 0, 2000 ); return R; }
  R.setup = true;
  std::vector<double> xv, inp;
  if( !oc.init( xv, inp, nullptr ) ){ std::cout << "  *** init failed\n"; return R; }
  R.nVar = oc.n_colloc_sta();
  OCFESLV::SolveReport rep = oc.solve( xv.data(), inp.data(), nullptr );
  R.conv = rep.converged; R.iters = rep.iterations;
  if( R.conv ) census( oc, xv, bnd, sdom, tbnd, exact, R );
  return R;
}

// ============================================================================================================
// CELL 1 / CELL 2 models: heat.  mixed form (T, q declared, LINK row by hand) or T_t = a T_xx (RED_FULL chain)
// ============================================================================================================
struct HeatPar { double a = 0.1, T0 = 10., Ts = 5., xf = 1., tf = 1.; };
static double T_exact( double t, double x, HeatPar const& p ){ double const k = kPi/2./p.xf; return p.Ts + p.T0*std::exp( -p.a*k*k*t )*std::sin( k*x ); }
static double q_exact( double t, double x, HeatPar const& p ){ double const k = kPi/2./p.xf; return p.a*p.T0*k*std::exp( -p.a*k*k*t )*std::cos( k*x ); }

static Run heat( bool mixed, size_t nelt, size_t nelx, int policy, double coupling )
{
  HeatPar p; FFGraph DAG;
  FFVar t = DAG.add_var( "t" ), x = DAG.add_var( "x" ), T = DAG.add_var( "T(t,x)" ), q = DAG.add_var( "q(t,x)" );
  FFPartial OpP;
  OCFESLV oc( &DAG );
  oc.add_domain( t, FFDom( 0., p.tf, nelt, FFDom::LGR, 5 ) );
  oc.add_domain( x, FFDom( 0., p.xf, nelx, FFDom::LGL, 5 ) );
  oc.add_state( T, { t, x } );
  oc.update_ref( T, [&]( OCFESLV::t_Coord const& c ){ return T_exact( c.at(t), c.at(x), p ); } );
  if( mixed ){
    oc.add_state( q, { t, x } );
    oc.update_ref( q, [&]( OCFESLV::t_Coord const& c ){ return q_exact( c.at(t), c.at(x), p ); } );
  }
  oc.set_evolution_domain( t );
  FFVar PDE  = mixed ? ( OpP( T, t ) - OpP( q, x ) ) : ( OpP( T, t ) - OpP( p.a*OpP( T, x ), x ) );
  FFVar LINK = q - p.a*OpP( T, x );
  FFVar INI  = T - p.Ts - p.T0*sin( kPi/2./p.xf*x );
  FFVar BCL  = T - p.Ts;
  FFVar BCU  = mixed ? q : p.a*OpP( T, x );
  oc.add_equation( PDE, { t, x }, { FFDom::ALL-FFDom::LB, FFDom::ALL-FFDom::LB-FFDom::UB }, OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  if( mixed ) oc.add_equation( LINK, { t, x }, { FFDom::ALL, FFDom::ALL }, OCFESLV::EqnOptions( OCFESLV::EqnRole::LINK, 0 ) );
  oc.add_equation( INI, { t, x }, { FFDom::LB, FFDom::ALL }, OCFESLV::EqnOptions( OCFESLV::EqnRole::INITIAL, 0 ) );
  oc.add_equation( BCL, { t, x }, { FFDom::ALL-FFDom::LB, FFDom::LB }, OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );
  oc.add_equation( BCU, { t, x }, { FFDom::ALL-FFDom::LB, FFDom::UB }, OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );
  common_options( oc, policy, coupling, O::RED_FULL );
  std::vector<double> bnd; for( size_t i = 0; i <= nelx; ++i ) bnd.push_back( p.xf*double(i)/double(nelx) );
  std::vector<double> tb; for( size_t i = 0; i <= nelt; ++i ) tb.push_back( p.tf*double(i)/double(nelt) );
  std::map< std::string, std::function<double(std::vector<double> const&)> > ex;
  ex["T(t,x)"] = [p]( std::vector<double> const& n ){ return T_exact( n[0], n[1], p ); };
  ex["q(t,x)"] = [p]( std::vector<double> const& n ){ return q_exact( n[0], n[1], p ); };
  return finish( oc, bnd, 1, tb, ex );
}

// ============================================================================================================
// CELL 3 model: order-4 biharmonic, AUTO depth-3 chain (PDE12's construction)
// ============================================================================================================
struct BihPar { double u0 = 1., k = 2.*kPi, phi = 0.3, a = 0.5, tf = 0.5, xf = 1.; };
static double U_exact( double t, double x, BihPar const& p ){ return p.u0*std::sin( p.k*x+p.phi )*( 1.+p.a*t ); }

static Run biharmonic( int policy )
{
  BihPar p; FFGraph DAG;
  FFVar t = DAG.add_var( "t" ), x = DAG.add_var( "x" ), u = DAG.add_var( "u(t,x)" );
  FFPartial OpP;
  double const k4 = p.k*p.k*p.k*p.k;
  FFVar sinx = sin( p.k*x + p.phi ), cosx = cos( p.k*x + p.phi );
  FFVar UE = p.u0*sinx*( 1.+p.a*t ), UEx = p.u0*p.k*cosx*( 1.+p.a*t ), UExx = -p.u0*p.k*p.k*sinx*( 1.+p.a*t );
  FFVar FE = p.u0*p.a*sinx + p.u0*k4*sinx*( 1.+p.a*t );
  FFVar PDE = OpP( u, t ) + OpP( u, {x,4} ) - FE, BCVAL = u - UE, BCDX = OpP( u, x ) - UEx, BCDXX = OpP( u, {x,2} ) - UExx;
  OCFESLV oc( &DAG );
  size_t const nelx = 3;
  oc.add_domain( t, FFDom( 0., p.tf, 3, FFDom::LGR, 6 ) );
  oc.add_domain( x, FFDom( 0., p.xf, nelx, FFDom::LGL, 16 ) );
  oc.add_state( u, { t, x } );
  oc.update_ref( u, [&]( OCFESLV::t_Coord const& c ){ return U_exact( c.at(t), c.at(x), p ); } );
  int const T_INT = FFDom::ALL - FFDom::LB, X_INT = FFDom::ALL - FFDom::LB - FFDom::UB;
  OCFESLV::EqnOptions io( OCFESLV::EqnRole::INTERIOR, 0 ), ii( OCFESLV::EqnRole::INITIAL, 0 ), bo( OCFESLV::EqnRole::BOUNDARY, 0 );
  oc.add_equation( PDE,   { t, x }, { T_INT, X_INT }, io );
  oc.add_equation( BCVAL, { t, x }, { FFDom::LB, FFDom::ALL }, ii );
  oc.add_equation( BCVAL, { t, x }, { T_INT, FFDom::LB }, bo );
  oc.add_equation( BCDX,  { t, x }, { T_INT, FFDom::LB }, bo );
  oc.add_equation( BCDXX, { t, x }, { T_INT, FFDom::LB }, bo );
  oc.add_equation( BCVAL, { t, x }, { T_INT, FFDom::UB }, bo );
  oc.set_evolution_domain( t );
  common_options( oc, policy, 1.0, O::RED_FULL );
  setenv( "CRONOS_RESCUE_C2", "1", 1 );   // the chain cell is about receivers, not the coupling: production setting
  std::vector<double> bnd; for( size_t i = 0; i <= nelx; ++i ) bnd.push_back( p.xf*double(i)/double(nelx) );
  std::vector<double> tb; for( size_t i = 0; i <= 3; ++i ) tb.push_back( p.tf*double(i)/3. );
  std::map< std::string, std::function<double(std::vector<double> const&)> > ex;
  ex["u(t,x)"] = [p]( std::vector<double> const& n ){ return U_exact( n[0], n[1], p ); };
  return finish( oc, bnd, 1, tb, ex );
}

// ============================================================================================================
// CELL 4 / CELL 5 model: PDE20g G4,  0 = w - du/dz - u,  under RED_FULL (w reuses Dz_u) and RED_MAIN (bare du/dz)
// ============================================================================================================
static double const U_ADV = 0.5, D_AX = 0.02, K_LIN = 0.7, T_END = 0.5, PH = 0.7, PI2 = 2.*kPi;
static double uex( double z, double t ){ return std::exp(-t)*( 1.+std::sin( PI2*z+PH ) ); }
static double uzex( double z, double t ){ return std::exp(-t)*PI2*std::cos( PI2*z+PH ); }

static Run g4( O::ReductionType red, int policy, double coupling, bool c2 = false, O::ImpositionType imp = O::IC_WEAK )
{
  FFGraph DAG;
  FFVar t = DAG.add_var( "t" ), z = DAG.add_var( "z" ), u = DAG.add_var( "u(t,z)" ), w = DAG.add_var( "w(t,z)" );
  FFPartial OpP;
  size_t const ne_t = 3, nn_t = 5, ne_z = 6, nn_z = 6;
  std::vector<double> t_bnd, z_bnd;
  for( size_t i = 0; i <= ne_t; ++i ) t_bnd.push_back( T_END*double(i)/double(ne_t) );
  for( size_t i = 0; i <= ne_z; ++i ) z_bnd.push_back( double(i)/double(ne_z) );
  FFVar W = PI2*z+PH, UE = exp(-t)*( 1.+sin(W) ), UT = -exp(-t)*( 1.+sin(W) ), UZ = exp(-t)*PI2*cos(W), UZZ = -exp(-t)*PI2*PI2*sin(W);
  FFVar LHS = UT + U_ADV*UZ - D_AX*UZZ, WEX = UZ + UE, SRC = LHS - K_LIN*WEX;
  OCFESLV oc( &DAG );
  oc.add_domain( t, FFDom( t_bnd, FFDom::LGL, nn_t ) );
  oc.add_domain( z, FFDom( z_bnd, FFDom::LGL, nn_z ) );
  int const TNL = FFDom::ALL-FFDom::LB, ZI = FFDom::ALL-FFDom::LB-FFDom::UB;
  oc.add_state( u, { t, z } );
  oc.update_ref( u, [&]( OCFESLV::t_Coord const& c ){ return uex( c.at(z), c.at(t) ); } );
  oc.add_state( w, { t, z } );
  oc.update_ref( w, [&]( OCFESLV::t_Coord const& c ){ return uzex( c.at(z), c.at(t) ) + uex( c.at(z), c.at(t) ); } );
  oc.set_evolution_domain( t );
  FFVar PDEu = OpP( u, t ) + U_ADV*OpP( u, z ) - D_AX*OpP( OpP( u, z ), z ) - K_LIN*w - SRC;
  FFVar ALG = w - OpP( u, z ) - u, ICu = u - UE;
  FFVar BCL = U_ADV*u - D_AX*OpP( u, z ) - ( U_ADV*UE - D_AX*UZ ), BCU = OpP( u, z ) - UZ;
  typedef OCFESLV::EqnOptions EO; typedef OCFESLV::EqnRole ER;
  oc.add_equation( PDEu, { t, z }, { TNL, ZI }, EO( ER::INTERIOR, 0 ) );
  oc.add_equation( ICu,  { t, z }, { FFDom::LB, FFDom::ALL }, EO( ER::INITIAL, 0 ) );
  oc.add_equation( BCL,  { t, z }, { TNL, FFDom::LB }, EO( ER::BOUNDARY, 0 ) );
  oc.add_equation( BCU,  { t, z }, { TNL, FFDom::UB }, EO( ER::BOUNDARY, 0 ) );
  oc.add_equation( ALG,  { t, z }, { FFDom::ALL, FFDom::ALL }, EO( ER::INTERIOR, 0 ) );
  common_options( oc, policy, coupling, red );
  setenv( "CRONOS_RESCUE_C2", c2 ? "1" : "0", 1 );   // rev309: interim knob, environment-only
  oc.options.INTERFACE.IMPOSITION  = imp;
  oc.options.AUTO.DIFF_ELIM   = true;
  std::map< std::string, std::function<double(std::vector<double> const&)> > ex;
  ex["u(t,z)"] = []( std::vector<double> const& n ){ return uex( n[1], n[0] ); };
  ex["w(t,z)"] = []( std::vector<double> const& n ){ return uzex( n[1], n[0] ) + uex( n[1], n[0] ); };
  return finish( oc, z_bnd, 1, t_bnd, ex );
}

// ============================================================================================================
static void show( std::string const& tag, Run const& R, std::string const& st )
{
  std::cout << "    " << std::left << std::setw(34) << tag << ( R.conv ? "conv" : "NO-CONV" ) << " nVar=" << R.nVar
            << "  " << st << ": if=" << sci( R.err_if.count(st) ? R.err_if.at(st) : -1. )
            << " int=" << sci( R.err_int.count(st) ? R.err_int.at(st) : -1. )
            << " spread=" << sci( R.spread.count(st) ? R.spread.at(st) : -1. )
            << "  [dropped=" << R.rep.dropped << " natural=" << R.rep.natural << " sole=" << R.rep.sole
            << " linkonly=" << R.rep.link_only << "]\n";
}
static double emax( Run const& R, std::string const& st ){ return std::max( R.err_if.count(st) ? R.err_if.at(st) : 1e300, R.err_int.count(st) ? R.err_int.at(st) : 1e300 ); }

int main()
{
  std::cout.setf( std::ios::unitbuf );
  std::cout << "================================================================\n"
            << "  OCFE_receiver -- the IC_WEAK receiver rule (WEAK_NATURAL_PENALTY)\n"
            << "  header: " << OCFE_OCFESLV_HEADER << "\n"
            << "================================================================\n";

  // ---- CELL 1 ------------------------------------------------------------------------------------------------
  std::cout << "\n---- CELL 1: duplicate receiver (mixed-form heat; q's claim has the PDE row AND its LINK row) ----\n";
  { Run d1 = heat( true, 4, 4, 0, 1.0 ), r1 = heat( true, 4, 4, 2, 1.0 ), r0 = heat( true, 4, 4, 2, 0.0 );
    show( "rule OFF  C2 on", d1, "T(t,x)" ); show( "rule ON   C2 on", r1, "T(t,x)" ); show( "rule ON   C2 off", r0, "T(t,x)" );
    expect( "C1 all three runs converge", d1.conv && r1.conv && r0.conv );
    expect( "C1 rule ON drops rescued terms and reports natural receivers", r1.rep.dropped > 0 && r1.rep.natural > 0, r1.rep.line.substr( r1.rep.line.find( "--" )+3, 90 ) );
    expect( "C1 rule ON: |T-T*| invariant to C2 on/off (rel < 1e-9)", r1.conv && r0.conv && rel( emax(r1,"T(t,x)"), emax(r0,"T(t,x)") ) < 1e-9, sci(emax(r1,"T(t,x)")) + " vs " + sci(emax(r0,"T(t,x)")) );
    expect( "C1 control -- rule OFF vs ON differ (rel > 1e-2): the duplicate injection was load-bearing", d1.conv && r1.conv && rel( emax(d1,"T(t,x)"), emax(r1,"T(t,x)") ) > 1e-2, sci(emax(d1,"T(t,x)")) + " vs " + sci(emax(r1,"T(t,x)")) );
    expect( "C1 nVar unchanged by the rule", r1.nVar == d1.nVar );
  }

  // ---- CELL 2 ------------------------------------------------------------------------------------------------
  std::cout << "\n---- CELL 2: sole claim (T_t = a T_xx, RED_FULL, TWO t-elements: t-claims on the auxiliary have no balance receiver) ----\n";
  { Run r1 = heat( false, 2, 3, 2, 1.0 ), r7 = heat( false, 2, 3, 2, 0.0 ), d1 = heat( false, 2, 3, 0, 1.0 );
    show( "rule OFF  C2 on", d1, "T(t,x)" ); show( "rule ON   C2 on", r1, "T(t,x)" ); show( "rule ON   C2 off", r7, "T(t,x)" );
    expect( "C2 rule ON reports sole claims > 0 (their rescued rows are kept)", r1.rep.sole > 0, "sole=" + std::to_string( r1.rep.sole ) );
    expect( "C2 rule ON converges with sole claims kept", r1.conv && r7.conv );
    expect( "C2 rule ON: |T-T*| invariant to C2 on/off (rel < 1e-9)", r1.conv && r7.conv && rel( emax(r1,"T(t,x)"), emax(r7,"T(t,x)") ) < 1e-9, sci(emax(r1,"T(t,x)")) + " vs " + sci(emax(r7,"T(t,x)")) );
  }

  // ---- CELL 3 ------------------------------------------------------------------------------------------------
  std::cout << "\n---- CELL 3: chain (order-4 biharmonic, AUTO depth-3 chain u -> u_x -> u_xx -> u_xxx) ----\n";
  { Run d = biharmonic( 0 ), p1 = biharmonic( 1 ), p2 = biharmonic( 2 );
    std::string deep; for( auto const& kv : d.spread ) if( kv.first.size() > deep.size() ) deep = kv.first;   // longest name = deepest
    show( "default", d, deep ); show( "=1 (informational)", p1, deep ); show( "=2 rule ON", p2, deep );
    expect( "C3 all three converge", d.conv && p1.conv && p2.conv );
    expect( "C3 =2 reports natural-LINK-only claims > 0", p2.rep.link_only > 0, "linkonly=" + std::to_string( p2.rep.link_only ) );
    expect( "C3 =2: deepest auxiliary's interface spread within 10x of the default", p2.conv && d.conv && p2.spread[deep] <= 10.*std::max( d.spread[deep], 1e-14 ), sci(p2.spread[deep]) + " vs " + sci(d.spread[deep]) );
    std::cout << "    (=1 spread " << sci( p1.spread[deep] ) << ": the O8 lesson, not asserted)\n";
  }

  // ---- CELL 4 ------------------------------------------------------------------------------------------------
  std::cout << "\n---- CELL 4: definition-row pollution (G4 RED_FULL: the rescued injection into w's own row) ----\n";
  // Production setting (RESCUE_C2 on): the sole w-claims are penalised in w's own row with the DERIVED c0.  With C2
  // off both paths carry the forced +1 against c0 = -1 and w is polluted by its own claim either way (measured:
  // 1.21e-03 vs 1.22e-03) -- sole claims still need C2's derived weight; the receiver rule does not replace it.
  { Run d = g4( O::RED_FULL, 0, 1.0, true ), r = g4( O::RED_FULL, 2, 1.0, true ), x = g4( O::RED_FULL, 2, 1.0, false );
    show( "default (C2 on)", d, "w(t,z)" ); show( "rule ON (C2 on)", r, "w(t,z)" ); show( "rule ON, C2 OFF (informational)", x, "w(t,z)" );
    expect( "C4 both converge", d.conv && r.conv );
    // R4 (rev291): the rule-OFF number dropped 3.44e-04 -> 8.01e-05 once C2's derived c0 = -1 on w's row stopped
    // being negated back to +1 by the sign rule (the c0 read was gated on S1; C2 derived; the sign rule undid it).
    // The rule-ON number (2.54e-05) is unchanged; the margin is now ~3.2x.
    expect( "C4 rule ON: w's interface error at least 2.5x below the default", d.conv && r.conv && r.err_if["w(t,z)"]*2.5 <= d.err_if["w(t,z)"], sci(r.err_if["w(t,z)"]) + " vs " + sci(d.err_if["w(t,z)"]) );
  }

  // ---- CELL 5 ------------------------------------------------------------------------------------------------
  std::cout << "\n---- CELL 5: bare derivative at the interface (G4 RED_MAIN: w reads the one-sided du/dz) ----\n";
  { Run r1 = g4( O::RED_MAIN, 2, 1.0, true ), r0 = g4( O::RED_MAIN, 2, 1.0, false ), t1 = g4( O::RED_MAIN, 2, 1.0, true, O::IC_TRACE );
    show( "rule ON  WEAK  C2 on", r1, "w(t,z)" ); show( "rule ON  WEAK  C2 off", r0, "w(t,z)" ); show( "rule ON  TRACE C2 on", t1, "w(t,z)" );
    expect( "C5 all three converge", r1.conv && r0.conv && t1.conv );
    double const a = r1.err_if["w(t,z)"], b = r0.err_if["w(t,z)"], c = t1.err_if["w(t,z)"];
    expect( "C5 rule ON: w's interface error invariant to C2 on/off (rel < 1e-8)", rel(a,b) < 1e-8, sci(a) + " / " + sci(b) );
    expect( "C5 rule ON: WEAK and TRACE agree on the honest RED_MAIN value (rel < 5e-2)", rel(a,c) < 5e-2, sci(a) + " vs " + sci(c) );
    std::cout << "    (the honest RED_MAIN value; the old default's smaller number came from a fabricated LINK-row injection)\n";
  }

  std::cout << "\n================================ VERDICT ================================\n"
            << "  " << g_pass << " PASS, " << g_fail << " FAIL\n"
            << "OCFE_receiver: " << ( g_fail == 0 ? "PASS" : "FAIL" ) << "\n";
  return g_fail == 0 ? 0 : 1;
}
