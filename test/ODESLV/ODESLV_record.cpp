// ODESLV_record -- the recorded trajectories and the three writers, mirroring DAESLV_record.
// =====================================================================================================
// record_solve / record_fsens / record_asens, one per results_* container.  record_asens is the one
// worth having: results_asens has been filled by every adjoint solve and, until the three-way split,
// no writer existed for it -- so this driver is the first thing to read it.
//
//   x' = -k x , x(0) = a , on [0,1];  f = int_0^1 x dt
//   x(t) = a e^{-k t}, so every recorded state row is checkable against the closed form.
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include "odeslv_cvodes.hpp"

using namespace mc;

static int n_pass = 0, n_fail = 0;
static void check( std::string const& w, double got, double want, double tol )
{
  double const d = std::fabs( got - want );
  ( d <= tol )? ++n_pass: ++n_fail;
  std::cout << "    " << std::left << std::setw( 40 ) << w << std::right
            << " got=" << std::scientific << std::setprecision(10) << got
            << " want=" << want << " |d|=" << std::setprecision(2) << d
            << ( d <= tol? "  PASS\n": "  FAIL\n" );
}
static void check_bool( std::string const& w, bool ok )
{
  ok? ++n_pass: ++n_fail;
  std::cout << "    " << std::left << std::setw( 40 ) << w << ( ok? "  PASS\n": "  FAIL\n" );
}
static size_t count_lines( std::string const& f )
{
  std::ifstream in( f ); size_t n = 0; std::string l;
  while( std::getline( in, l ) ) if( !l.empty() ) ++n;
  return n;
}

int main()
{
  FFPartial OpP;  FFIntegral OpI;
  int const TI = FFDom::ALL - FFDom::LB;
  auto io = FFModel::EqnOptions( FFModel::EqnRole::INTERIOR, 0 );
  auto ii = FFModel::EqnOptions( FFModel::EqnRole::INITIAL,  0 );

  FFGraph G;  ODESLV_CVODES I( &G );
  FFVar t = G.add_var("t"), x = G.add_var("x(t)"), k = G.add_var("k"), a = G.add_var("a");
  I.add_domain( t, FFDom( 0., 1., 4, FFDom::LGR, 4 ) );
  I.add_state( x, {t} );  I.add_input( k );  I.add_input( a );
  I.set_evolution_domain( t );  I.update_ref( x, 1. );
  I.add_equation( OpP( x, t ) + k*x, {t}, {TI}, io );
  I.add_equation( x - a,             {t}, {FFDom::LB}, ii );
  I.add_output( OpI( std::vector<FFVar>{ x }, { t } )[0] );
  I.options.SOLVE.RTOL  = I.options.SOLVE.ATOL  = 1e-11;
  I.options.SOLVE.RTOLB = I.options.SOLVE.ATOLB = 1e-11;
  I.options.OUTPUT.RESRECORD = 1;     // recording is OPT-IN in both solvers
  I.options.DISPLAY_LEVEL = 0;
  // ODESLV fills results_* only when RESRECORD is set -- its recording is OPT-IN, where DAESLV's is
  // unconditional.  Without this every container comes back empty while the solve prints normally.
  I.options.OUTPUT.RESRECORD = 1;
  if( !I.setup() ){ std::cout << "setup failed: " << I.extract_error() << "\n"; return 1; }

  double const K = 1.3, A = 2.1;
  size_t const jk = I.parameter_index( k )[0], ja = I.parameter_index( a )[0];
  std::vector<double> P( I.np(), 0. );  P[jk] = K;  P[ja] = A;

  // ------------------------------------------------------------------ A: the STATE trajectory
  { std::cout << "\n--- A. results_solve and record_solve()\n";
    I.solve( P );
    auto const& R = I.results_solve;
    std::cout << "    rows=" << R.size() << "  width=" << ( R.empty()? 0: R[0].x.size() )
              << "  (nx=" << I.nx() << ", nq=" << I.nq() << ")\n";
    check_bool( "rows were recorded", !R.empty() );
    if( !R.empty() ){
      check( "row width = nx + nq", (double)R[0].x.size(), (double)(I.nx()+I.nq()), 0. );
      check( "t ASCENDS: first row", R.front().t, 0., 1e-12 );
      check( "              last row", R.back().t, 1., 1e-9 );
      check( "x(0) = a",      R.front().x[0], A, 1e-9 );
      check( "x(1) = a e^-k", R.back().x[0],  A*std::exp(-K), 1e-7 );
      bool allok = true;
      for( auto const& r : R ) if( std::fabs( r.x[0] - A*std::exp(-K*r.t) ) > 1e-6 ) allok = false;
      check_bool( "EVERY row matches a e^{-k t}", allok );
      std::ofstream os( "ODESLV_record_sta.dat" );  I.record_solve( os );  os.close();
      check( "record_solve() wrote one line per row",
             (double)count_lines( "ODESLV_record_sta.dat" ), (double)R.size(), 0. );
    } }

  // ------------------------------------------------------------------ B: the SENSITIVITY trajectories
  { std::cout << "\n--- B. results_fsens and record_fsens()\n";
    I.solve_fsens( P );
    auto const& S = I.results_fsens;
    std::cout << "    directions=" << S.size() << "  rows=" << ( S.empty()? 0: S[0].size() )
              << "  width=" << ( S.empty()||S[0].empty()? 0: S[0][0].x.size() ) << "\n";
    check_bool( "one vector per sensitivity direction", !S.empty() );
    if( !S.empty() && !S[ja].empty() ){
      // ODESLV records [ dy, dyq ]: the sensitivity state AND its quadrature
      check( "row width = nx + nq", (double)S[ja][0].x.size(), (double)(I.nx()+I.nq()), 0. );
      check( "d x(0)/da = 1",       S[ja].front().x[0], 1.,               1e-9 );
      check( "d x(0)/dk = 0",       S[jk].front().x[0], 0.,               1e-9 );
      check( "d x(1)/da = e^-k",    S[ja].back().x[0],  std::exp(-K),     1e-7 );
      check( "d x(1)/dk = -a e^-k", S[jk].back().x[0], -A*std::exp(-K),   1e-7 );
      std::ofstream osen[2];
      osen[0].open( "ODESLV_record_sen0.dat" );  osen[1].open( "ODESLV_record_sen1.dat" );
      I.record_fsens( osen );  osen[0].close();  osen[1].close();
      check( "record_fsens() wrote one line per row",
             (double)count_lines( "ODESLV_record_sen0.dat" ), (double)S[0].size(), 0. );
    } }

  // ------------------------------------------------------------------ C: the ADJOINT trajectory
  { std::cout << "\n--- C. results_asens and record_asens() -- no writer existed before the split\n";
    I.solve_asens( P );
    auto const& Ad = I.results_asens;
    std::cout << "    outputs=" << Ad.size() << "  rows=" << ( Ad.empty()? 0: Ad[0].size() )
              << "  width=" << ( Ad.empty()||Ad[0].empty()? 0: Ad[0][0].x.size() ) << "\n";
    check_bool( "one vector per output function", !Ad.empty() );
    if( !Ad.empty() && !Ad[0].empty() ){
      check( "row width = nx + nsen", (double)Ad[0][0].x.size(), (double)(I.nx()+I.np()), 0. );
      check_bool( "t DESCENDS, as the backward sweep produces it",
                  Ad[0].front().t > Ad[0].back().t );
      check( "the sweep reached t0", Ad[0].back().t, 0., 1e-9 );
      // lambda(t0) carries dG/dx(0); with x(0)=a that is d f/da = (1-e^-k)/k
      check( "lambda(t0) = d f / da", Ad[0].back().x[0], ( 1.-std::exp(-K) )/K, 1e-6 );
      std::ofstream oadj[1];  oadj[0].open( "ODESLV_record_adj.dat" );
      I.record_asens( oadj );  oadj[0].close();
            { std::ifstream in( "ODESLV_record_adj.dat" );
        double t0 = 0., t1 = 0., v = 0.; bool first = true, asc = true;
        std::string line;
        while( std::getline( in, line ) ){
          if( line.empty() ) continue;
          std::istringstream is( line ); is >> v;
          if( first ){ t0 = v; first = false; } else if( v < t1 ) asc = false;
          t1 = v;
        }
        check_bool( "the FILE reads in INCREASING time", asc && t1 > t0 ); }
      check( "record_asens() wrote one line per row",
             (double)count_lines( "ODESLV_record_adj.dat" ), (double)Ad[0].size(), 0. );
    } }

  // ------------------------------------------------------------------ D: the new accessors
  { std::cout << "\n--- D. val_state / val_state_init, as DAESLV exposes them\n";
    auto const& X = I.val_state();
    check_bool( "val_state() returns the stage-time states", !X.empty() );
    if( !X.empty() ){
      check( "x at the last stage time = a e^-k", X.back().at(0), A*std::exp(-K), 1e-7 );
      check( "val_state_init() = x(0) = a",       I.val_state_init().at(0), A, 1e-9 );
    } }

  std::cout << "\n==========  ODESLV_record: " << n_pass << " passed, " << n_fail << " failed -- "
            << ( n_fail? "FAILURES": "ALL PASS" ) << "  ==========\n\n";
  return n_fail? 1: 0;
}
