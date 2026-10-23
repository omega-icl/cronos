// DAESLV_asensmulti -- the adjoint against the forward on a model with SEVERAL outputs.
// =====================================================================================================
// check_daeslv.py disagrees here and no C++ driver did, because every existing adjoint case adjoins
// ONE output.  solve_asens() now loops the outputs (ODESLV's behaviour, adopted when ifct was dropped
// from its signature), and this is the first thing to exercise that loop.
//
//   x' = -k x , x(0) = 1 , on [0,1] in 4 stages
//   outputs: x at each stage time, then x(1) again -- 6 in all, so the loop runs six sweeps
//   x(tau) = e^{-k tau} , d x(tau)/dk = -tau e^{-k tau}
#include <cmath>
#include <iomanip>
#include <iostream>
#include "daeslv_idas.hpp"

using namespace mc;

static int n_pass = 0, n_fail = 0;
static void check( std::string const& w, double got, double want, double tol )
{
  double const d = std::fabs( got - want );
  ( d <= tol )? ++n_pass: ++n_fail;
  std::cout << "    " << std::left << std::setw( 34 ) << w << std::right
            << " got=" << std::scientific << std::setprecision(10) << got
            << " want=" << want << " |d|=" << std::setprecision(2) << d
            << ( d <= tol? "  PASS\n": "  FAIL\n" );
}

int main()
{
  FFPartial OpP;  FFEval OpE;
  int const TI = FFDom::ALL - FFDom::LB;
  auto io = FFModel::EqnOptions( FFModel::EqnRole::INTERIOR, 0 );
  auto ii = FFModel::EqnOptions( FFModel::EqnRole::INITIAL,  0 );

  FFGraph G;  DAESLV_IDAS M;  M.set_dag( &G );
  FFVar t = G.add_var("t"), x = G.add_var("x(t)"), k = G.add_var("k");
  M.add_domain( t, FFDom( 0., 1., 4, FFDom::LGR, 3 ) );
  M.add_state( x, {t} );  M.add_input( k );
  M.set_evolution_domain( t );  M.update_ref( x, 1. );
  M.add_equation( OpP( x, t ) + k*x, {t}, {TI}, io );
  M.add_equation( x - 1.,            {t}, {FFDom::LB}, ii );
  double const TAU[5] = { 0., 0.25, 0.5, 0.75, 1. };
  for( double tv : TAU ) M.add_output( OpE( x, t, tv ) );
  M.add_output( OpE( x, t, 1. ) );                 // a sixth, repeating the last
  M.options.SOLVE.RTOL  = M.options.SOLVE.ATOL  = 1e-11;
  M.options.SOLVE.RTOLB = M.options.SOLVE.ATOLB = 1e-11;
  M.options.DISPLAY_LEVEL = 0;
  if( !M.setup() ){ std::cout << "setup failed: " << M.extract_error() << "\n"; return 1; }
  std::cout << "\n--- nf=" << M.nf() << " outputs, np=" << M.np() << "\n";

  double const K = 0.8;
  std::vector<double> P( M.np(), 0. );  size_t const jk = M.parameter_index( k )[0];  P[jk] = K;

  M.solve_fsens( P );
  std::vector< std::vector<double> > const Gf = M.val_function_gradient();
  M.solve_asens( P );
  auto const& Ga = M.val_function_gradient();

  std::cout << "    fwd:";
  for( size_t r = 0; r < M.nf(); ++r ) std::cout << " " << std::scientific << std::setprecision(4)
            << Gf.at(jk).at( M.blk_fct(r).first );
  std::cout << "\n    adj:";
  for( size_t r = 0; r < M.nf(); ++r ) std::cout << " " << Ga.at(jk).at( M.blk_fct(r).first );
  std::cout << "\n";

  for( size_t r = 0; r < 5; ++r ){
    size_t const b = M.blk_fct(r).first;
    double const want = -TAU[r]*std::exp(-K*TAU[r]);
    check( "fwd d x(" + std::to_string(TAU[r]).substr(0,4) + ")/dk", Gf.at(jk).at(b), want, 1e-6 );
    check( "adj d x(" + std::to_string(TAU[r]).substr(0,4) + ")/dk", Ga.at(jk).at(b), want, 1e-6 );
  }
  { size_t const b = M.blk_fct(5).first;
    check( "fwd the repeated output", Gf.at(jk).at(b), -std::exp(-K), 1e-6 );
    check( "adj the repeated output", Ga.at(jk).at(b), -std::exp(-K), 1e-6 ); }

  std::cout << "\n==========  DAESLV_asensmulti: " << n_pass << " passed, " << n_fail << " failed -- "
            << ( n_fail? "FAILURES": "ALL PASS" ) << "  ==========\n\n";
  return n_fail? 1: 0;
}
