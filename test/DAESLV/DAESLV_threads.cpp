// DAESLV_threads -- FFDAESLV under veval on SEVERAL threads.
// =====================================================================================================
// A DAG operation makes one solver per thread: FFBaseDE's copy constructor does `new SOLVER` then
// `setup( *Op._pODESLV )`, an independent copy of the description.  That setup EXTRACTS into the user
// DAG, which the threads share, so it has to be serialised -- FFModel::dag_mutex() exists for exactly
// that, and its own note says setups are locked while solves are not, each solver evaluating on its
// own working DAG afterwards.
//
// ODESLV_CVODES::setup() takes that lock.  DAESLV_IDAS::setup() did not, so a multi-threaded veval of
// an FFDAESLV raced in extraction (reported from the tutorial's scenario sweep, 2026-10-24: fine at
// MAXTHREAD=1, throwing above it, while the same sweep over an FFODESLV was fine).
//
// This driver runs the same sweep on one thread and on many, and requires the results to AGREE --
// not merely to avoid throwing, since a race that corrupts an extraction can equally produce wrong
// numbers quietly.
#include <cmath>
#include <iomanip>
#include <iostream>
#include <thread>
#include "ffde.hpp"

using namespace mc;

static int n_pass = 0, n_fail = 0;
static void check_bool( std::string const& w, bool ok )
{
  ok? ++n_pass: ++n_fail;
  std::cout << "    " << std::left << std::setw( 46 ) << w << ( ok? "  PASS\n": "  FAIL\n" ) << std::flush;
}

// the model, embedded and swept over scenarios.  `dae` picks the index-2 vented reactor, which only
// DAESLV can take; the ODESLV control needs an ODE, so it gets the same reactor with the vent rate
// given explicitly instead of by a constraint.
template <class SOLVER, class FFOP>
static bool sweep( std::string const& tag, unsigned const nthread, std::vector<double>& out,
                   bool const dae, int const policy )
{
  FFPartial OpP;  FFEval OpE;
  int const TI = FFDom::ALL - FFDom::LB;
  auto const INTER = FFModel::EqnRole::INTERIOR;
  auto const INIT  = FFModel::EqnRole::INITIAL;

  FFGraph G;  SOLVER M;  M.set_dag( &G );
  FFVar t = G.add_var("t"), nA = G.add_var("nA(t)"), nB = G.add_var("nB(t)"), q = G.add_var("q(t)"),
        k = G.add_var("k"), n0 = G.add_var("nA0"), nT = G.add_var("nT");
  M.add_domain( t, FFDom( 0., 1., 4, FFDom::LGR, 3 ) );
  M.add_state( nA, {t} );  M.add_state( nB, {t} );
  if( dae ) M.add_state( q, {t} );        // the ODE control has no algebraic state at all
  M.add_input( k );  M.add_input( n0 );
  M.set_constant( {nT} );
  M.set_evolution_domain( t );
  M.update_ref( nA, 0.5 );  M.update_ref( nB, 0.5 );
  if( dae ) M.update_ref( q, 0.1 );
  if( dae ){                              // index 2: the vent rate is fixed by the constraint
    M.add_equation( OpP( nA, t ) + k*nA + q*nA,    {t}, {TI}, INTER );
    M.add_equation( OpP( nB, t ) - 2.*k*nA + q*nB, {t}, {TI}, INTER );
    //M.add_equation( nA + nB - 1.0,                 {t}, {FFDom::ALL}, INTER );
    M.add_equation( nA + nB - nT,                  {t}, {FFDom::ALL}, INTER );
  }
  else{                                   // the ODE control: q = k nA substituted in
    M.add_equation( OpP( nA, t ) + k*nA + (k*nA)*nA,    {t}, {TI}, INTER );
    M.add_equation( OpP( nB, t ) - 2.*k*nA + (k*nA)*nB, {t}, {TI}, INTER );
  }
  M.add_equation( nA - n0,                       {t}, {FFDom::LB}, INIT );
  if( !dae ) M.add_equation( nB - 0.2,           {t}, {FFDom::LB}, INIT );
  M.add_output( OpE( nA, t, 1. ) );
  M.options.SOLVE.RTOL = M.options.SOLVE.ATOL = 1e-10;
  M.options.DISPLAY_LEVEL = 0;
  M.options.SOLVE.VERBOSE = 0;
  if( !M.setup() ){ std::cout << "    " << tag << " setup failed: " << M.extract_error() << "\n"; return false; }

  FFOP Op;
  //FFVar Pk = G.add_var("Pk"), Pn = G.add_var("Pn");
  //std::vector<FFModel::InputArg> vArg{ { k, std::vector<FFVar>{ Pk } },
  //                                     { n0, std::vector<FFVar>{ Pn } } };
  std::vector<FFModel::InputArg> vArg{ { k,  std::vector<FFVar>{ k } },
                                       { n0, std::vector<FFVar>{ n0 } } };
  std::vector<FFModel::InputArg> vCst{ { nT, std::vector<FFVar>{ nT } } };
  // COPY, not SHALLOW.  SHALLOW leaves _ownODESLV false, so FFBaseDE's copy constructor SHARES the
  // solver and never calls setup() -- the very path dag_mutex() protects.  A threaded test on a
  // SHALLOW op therefore passes vacuously, which is what this driver did when first written.
  //std::vector<FFVar> vF = Op( vArg, &M, policy );
  std::vector<FFVar> vF = Op( vArg, vCst, &M, policy );
  if( vF.empty() ){ std::cout << "    " << tag << ": the embedding returned nothing\n"; return false; }

  std::vector< std::vector<double> > sample;
  for( int i = 0; i < 24; ++i )
    sample.push_back( { 0.5, 0.4 } );//0.5 + 0.1*i, 0.4 + 0.02*i } );

  G.options.MAXTHREAD = nthread;
  std::vector< std::vector<double> > res;
  //G.veval( vF, res, std::vector<FFVar>{ Pk, Pn }, sample );
  G.veval( vF, res, std::vector<FFVar>{ k, n0 }, sample, std::vector<FFVar>{ nT }, std::vector<double>{ 1.0 } );
  out.clear();
  for( auto const& r : res ) for( double v : r ) out.push_back( v );
  return true;
}

int main()
{
  // A race cannot be observed on one core, and a green result there would be worthless: say so.
  unsigned const ncore = std::thread::hardware_concurrency();
  std::cout << "\n    hardware_concurrency = " << ncore << "\n";
  if( ncore < 2 )
    std::cout << "    *** ONE CORE: MAXTHREAD=0 resolves to a single thread, so the agreement below\n"
                 "    *** is one thread against one thread and proves NOTHING about concurrency.\n";

  // ---------------------------------------------------------------- the ODE operation, as a control
  std::cout << "\n--- FFODESLV with COPY: the control, already fine on many threads\n";
  { std::vector<double> one, many;
    std::cout << "starting single threaded run\n";
    bool const a = sweep<ODESLV_CVODES,FFODESLV>( "ODESLV 1", 1, one, false, FFODESLV::COPY );
    std::cout << "starting multiple threaded run\n";
    bool const b = sweep<ODESLV_CVODES,FFODESLV>( "ODESLV N", 0, many, false, FFODESLV::COPY );   // 0 = every thread
    check_bool( "both sweeps ran", a && b );
    if( a && b ){
      double worst = 0.;
      for( size_t i = 0; i < one.size() && i < many.size(); ++i )
        worst = std::max( worst, std::fabs( one[i] - many[i] ) );
      std::cout << "    " << std::left << std::setw( 46 ) << "max |1 thread - N threads|"
                << std::scientific << std::setprecision(2) << worst << "\n" << std::flush;
      check_bool( "same sizes", one.size() == many.size() );
      check_bool( "and the SAME numbers", worst < 1e-9 );
    } }

  // ---------------------------------------------------------------- the DAE operation
  std::cout << "\n--- FFDAESLV with COPY: one setup() PER THREAD, which is what raced\n";
  { std::vector<double> one, many;
    std::cout << "starting single threaded run\n";
    bool const a = sweep<DAESLV_IDAS,FFDAESLV>( "DAESLV 1", 1, one, true, FFDAESLV::COPY );
    std::cout << "starting multiple threaded run\n";
    bool const b = sweep<DAESLV_IDAS,FFDAESLV>( "DAESLV N", 0, many, true, FFDAESLV::COPY );
    check_bool( "both sweeps ran", a && b );
    if( a && b ){
      double worst = 0.;
      for( size_t i = 0; i < one.size() && i < many.size(); ++i )
        worst = std::max( worst, std::fabs( one[i] - many[i] ) );
      std::cout << "    " << std::left << std::setw( 46 ) << "max |1 thread - N threads|"
                << std::scientific << std::setprecision(2) << worst << "\n" << std::flush;
      check_bool( "same sizes", one.size() == many.size() );
      check_bool( "and the SAME numbers", worst < 1e-9 );
    } }

  std::cout << "\n==========  DAESLV_threads: " << n_pass << " passed, " << n_fail << " failed -- "
            << ( n_fail? "FAILURES": "ALL PASS" ) << "  ==========\n\n";
  return n_fail? 1: 0;
}
