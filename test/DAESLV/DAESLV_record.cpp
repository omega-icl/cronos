// DAESLV_record -- the recorded trajectories and record(), mirroring ODESLV.
// =====================================================================================================
// results_solve / results_fsens / results_asens are filled where _xk, _sk and the backward lambda
// already flow, and _record (hoisted into BASE_DE so both solvers share ONE writer) puts them on
// a stream.  Derivatives are deliberately NOT recorded: the row format is ODESLV's, [ state, quad ].
//
//   x' = -k x , x(0) = a , on [0,1] in 4 stages;  f = int_0^1 x dt
//   x(t) = a e^{-k t}, so every recorded row is checkable against the closed form.
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include "daeslv_idas.hpp"

using namespace mc;

static int n_pass = 0, n_fail = 0;
static void check( std::string const& w, double got, double want, double tol )
{
  double const d = std::fabs( got - want );
  ( d <= tol )? ++n_pass: ++n_fail;
  std::cout << "    " << std::left << std::setw( 38 ) << w << std::right
            << " got=" << std::scientific << std::setprecision(10) << got
            << " want=" << want << " |d|=" << std::setprecision(2) << d
            << ( d <= tol? "  PASS\n": "  FAIL\n" );
}
static void check_bool( std::string const& w, bool ok )
{
  ok? ++n_pass: ++n_fail;
  std::cout << "    " << std::left << std::setw( 38 ) << w << ( ok? "  PASS\n": "  FAIL\n" );
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

  FFGraph G; DAESLV_IDAS M; M.set_dag( &G );
  FFVar t = G.add_var("t"), x = G.add_var("x(t)"), k = G.add_var("k"), a = G.add_var("a");
  M.add_domain( t, FFDom( 0., 1., 4, FFDom::LGR, 4 ) );
  M.add_state( x, {t} ); M.add_input( k ); M.add_input( a );
  M.set_evolution_domain( t ); M.update_ref( x, 1. );
  M.add_equation( OpP( x, t ) + k*x, {t}, {TI}, io );
  M.add_equation( x - a,             {t}, {FFDom::LB}, ii );
  M.add_output( OpI( std::vector<FFVar>{ x }, { t } )[0] );
  M.options.SOLVE.RTOL  = M.options.SOLVE.ATOL  = 1e-11;
  M.options.SOLVE.RTOLB = M.options.SOLVE.ATOLB = 1e-11;
  unsigned const NREC = 10;
  M.options.OUTPUT.RESRECORD = NREC;   // opt-in, and NREC recorded points per stage
  M.options.DISPLAY_LEVEL = 0;
  if( !M.setup() ){ std::cout << "setup failed: " << M.extract_error() << "\n"; return 1; }

  double const K = 1.3, A = 2.1;
  size_t const jk = M.parameter_index( k )[0], ja = M.parameter_index( a )[0];
  std::vector<double> P( M.np(), 0. ); P[jk] = K; P[ja] = A;

  // ------------------------------------------------------------------ A: the STATE trajectory
  { std::cout << "\n--- A. results_solve, and record() to a stream\n";
    M.solve( P );
    auto const& R = M.results_solve;
    std::cout << "    rows=" << R.size() << "  width=" << ( R.empty()? 0: R[0].x.size() )
              << "  (stages=" << M.ns() << ", nx=" << M.nx() << ", nq=" << M.nq() << ")\n";
    // one row per stage time PLUS a stage-START row at every interior boundary (the repeated t that
    // _record turns into a block separator, and that shows a discontinuity across a jump)
    check( "rows = stage times + block starts", (double)R.size(), (double)( 1 + M.ns()*NREC + M.ns()-1 ), 0. );
    check( "row width = nx + nq",      (double)R[0].x.size(), (double)(M.nx()+M.nq()), 0. );
    check( "sensitivity rows carry the quadrature too", 0., 0., 0. );
    check( "t ASCENDS: first row",     R.front().t, 0., 0. );
    check( "              last row",   R.back().t,  1., 1e-12 );
    check( "x(0) = a",                 R.front().x[0], A, 1e-9 );
    check( "x(1) = a e^-k",            R.back().x[0],  A*std::exp(-K), 1e-8 );
    // every row against the closed form
    bool allok = true;
    for( auto const& r : R ) if( std::fabs( r.x[0] - A*std::exp(-K*r.t) ) > 1e-7 ) allok = false;
    check_bool( "EVERY row matches a e^{-k t}", allok );

    std::ofstream os( "DAESLV_record_sta.dat" );
    M.record_solve( os );  os.close();
    check( "record_solve() wrote one line per row", (double)count_lines( "DAESLV_record_sta.dat" ), (double)R.size(), 0. ); }

  // ------------------------------------------------------------------ B: the SENSITIVITY trajectories
  { std::cout << "\n--- B. results_fsens, one vector per parameter\n";
    M.solve_fsens( P );
    auto const& S = M.results_fsens;
    check( "one vector per parameter", (double)S.size(), (double)M.np(), 0. );
    check( "rows = stage times + block starts", (double)S[ja].size(), (double)( 1 + M.ns()*NREC + M.ns()-1 ), 0. );
    // d x / da = e^{-k t}, d x / dk = -t a e^{-k t}
    check( "d x(1)/da = e^-k",    S[ja].back().x[0],  std::exp(-K),      1e-7 );
    check( "d x(1)/dk = -a e^-k", S[jk].back().x[0], -A*std::exp(-K),    1e-7 );
    check( "d x(0)/da = 1",       S[ja].front().x[0], 1.,                1e-9 );
    check( "d x(0)/dk = 0",       S[jk].front().x[0], 0.,                1e-9 );
    std::ofstream osen[2]; osen[0].open( "DAESLV_record_sen0.dat" ); osen[1].open( "DAESLV_record_sen1.dat" );
    std::ofstream os( "DAESLV_record_sta2.dat" );
    M.record_solve( os );  M.record_fsens( osen );  os.close(); osen[0].close(); osen[1].close();
    check_bool( "record_solve + record_fsens wrote both",
                count_lines( "DAESLV_record_sta2.dat" ) == S[0].size()
             && count_lines( "DAESLV_record_sen0.dat" ) == S[0].size() ); }

  // ------------------------------------------------------------------ C: the ADJOINT trajectory
  { std::cout << "\n--- C. results_asens, in DESCENDING time as ODESLV stores it\n";
    M._solve_asens_one( P, 0 );
    auto const& A_ = M.results_asens;
    check( "one vector per output", (double)A_.size(), (double)M.nf(), 0. );
    check_bool( "the adjoint sweep recorded rows", !A_[0].empty() );
    if( !A_[0].empty() ){
      std::cout << "    rows=" << A_[0].size() << "  t: " << std::scientific << std::setprecision(3)
                << A_[0].front().t << " -> " << A_[0].back().t << "\n";
      check_bool( "t DESCENDS, as the backward sweep produces it",
                  A_[0].front().t > A_[0].back().t );
      check( "the sweep reached t0", A_[0].back().t, 0., 1e-12 );
      check( "adjoint rows = NREC per stage + block ends",
             (double)A_[0].size(), (double)( M.ns()*NREC + M.ns() + 1 ), 0. );   // +1: the t0 SUMMARY row
      std::ofstream oadj[1]; oadj[0].open( "DAESLV_record_adj.dat" );
      M.record_asens( oadj );  oadj[0].close();
            { std::ifstream in( "DAESLV_record_adj.dat" );
        double t0 = 0., t1 = 0., v = 0.; bool first = true, asc = true;
        std::string line;
        while( std::getline( in, line ) ){
          if( line.empty() ) continue;
          std::istringstream is( line ); is >> v;
          if( first ){ t0 = v; first = false; } else if( v < t1 ) asc = false;
          t1 = v;
        }
        check_bool( "the FILE reads in INCREASING time", asc && t1 > t0 ); }
      // the recorded lambda_q at t0 must BE the gradient the solver reports -- tying the trajectory
      // to the number, so the two cannot drift apart again
      { auto const& Gd = M.val_function_gradient();
        check( "recorded lambda_q(t0) = d f/dk", A_[0].back().x.at(1+jk), Gd[jk][0], 1e-9 );
        check( "recorded lambda_q(t0) = d f/da", A_[0].back().x.at(1+ja), Gd[ja][0], 1e-9 ); }
      // val_state_adjoint(): [stage][output][lambda], as ODESLV exposes it
      { auto const& L = M.val_state_adjoint();
        check( "val_state_adjoint: one entry per stage+1", (double)L.size(), (double)(M.ns()+1), 0. );
        check( "  one per output",      (double)L.front().size(), (double)M.nf(), 0. );
        check( "  lambda(t0) = d f/da", L.front().at(0).at(0), M.val_function_gradient()[ja][0], 1e-9 );
        check( "  lambda(T)  = 0",      L.back().at(0).at(0), 0., 1e-9 ); }
      check( "record_asens() wrote one line per row",
             (double)count_lines( "DAESLV_record_adj.dat" ), (double)A_[0].size(), 0. );
    } }

  std::cout << "\n==========  DAESLV_record: " << n_pass << " passed, " << n_fail << " failed -- "
            << ( n_fail? "FAILURES": "ALL PASS" ) << "  ==========\n\n";
  return n_fail? 1: 0;
}
