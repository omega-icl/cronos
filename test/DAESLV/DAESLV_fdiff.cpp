// DAESLV_fdiff -- the solver as a DIRECTIONAL DERIVATIVE of itself.
// =====================================================================================================
// fdiff() returns a NEW DAESLV_IDAS whose model is the derivative of this one along the given
// parameters.  Solving THAT solver gives the same numbers solve_fsens() gives, by a different route:
// fsens differentiates the integration, fdiff differentiates the MODEL and integrates the result.
// Agreement between the two is the gate, and it is what FFDAESLV will rest on.
//
//   x' = -k x , x(0) = a ;  f = int_0^1 x dt
//   d f/da = (1-e^-k)/k ,  d f/dk = a[ k e^-k - (1-e^-k) ]/k^2
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

int main()
{
  FFPartial OpP;  FFIntegral OpI;
  int const TI = FFDom::ALL - FFDom::LB;
  auto io = FFModel::EqnOptions( FFModel::EqnRole::INTERIOR, 0 );
  auto ii = FFModel::EqnOptions( FFModel::EqnRole::INITIAL,  0 );

  FFGraph G;  DAESLV_IDAS M;  M.set_dag( &G );
  FFVar t = G.add_var("t"), x = G.add_var("x(t)"), k = G.add_var("k"), a = G.add_var("a");
  M.add_domain( t, FFDom( 0., 1., 1, FFDom::LGR, 5 ) );
  M.add_state( x, {t} );  M.add_input( k );  M.add_input( a );
  M.set_evolution_domain( t );  M.update_ref( x, 1. );
  M.add_equation( OpP( x, t ) + k*x, {t}, {TI}, io );
  M.add_equation( x - a,             {t}, {FFDom::LB}, ii );
  M.add_output( OpI( std::vector<FFVar>{ x }, { t } )[0] );
  M.options.SOLVE.RTOL = M.options.SOLVE.ATOL = 1e-11;
  M.options.DISPLAY_LEVEL = 0;
  if( !M.setup() ){ std::cout << "setup failed: " << M.extract_error() << "\n"; return 1; }

  double const K = 1.3, A = 2.1, e = std::exp(-K);
  size_t const jk = M.parameter_index( k )[0], ja = M.parameter_index( a )[0];
  std::vector<double> P( M.np(), 0. );  P[jk] = K;  P[ja] = A;

  // the reference: forward sensitivities of the ORIGINAL solver, already gated elsewhere
  M.solve_fsens( P );
  std::vector<double> const Gref = M.val_function_gradient();
  size_t const b = M.blk_fct(0).first;
  std::cout << "\n--- reference: solve_fsens on the original model\n";
  check( "d f/da = (1-e^-k)/k", Gref.at( b*M.np()+ja ), ( 1.-e )/K, 1e-7 );
  check( "d f/dk",              Gref.at( b*M.np()+jk ), A*( K*e - (1.-e) )/(K*K), 1e-7 );

  // ------------------------------------------------------------------ A: fdiff over both parameters
  { std::cout << "\n--- A. fdiff( {k,a} ): a NEW solver whose model is the derivative\n";
    std::string err;
    std::vector<FFVar> vU{ k, a };
    DAESLV_IDAS* S = M.fdiff( vU, err );
    check_bool( "fdiff returned a solver", S != nullptr );
    if( !S ){ std::cout << "    err: " << err << "\n"; }
    else{
      check_bool( "it is a DAESLV_IDAS", dynamic_cast<DAESLV_IDAS*>( S ) != nullptr );
      if( !S->setup() ){ std::cout << "    setup of the derivative model failed: " << S->extract_error() << "\n"; ++n_fail; }
      else{
        std::cout << "    derivative model: np=" << S->np() << " nf=" << S->nf()
                  << " nx=" << S->nx() << "\n";
        check_bool( "it carries more outputs than the original", S->nf() > M.nf() );
        // The derivative model's parameters are the ORIGINALS plus one DIRECTION input per
        // parameter differentiated.  Every one needs a value -- a missing d1[k] is what the dag
        // reports, not a defaulted zero.  They are filled by NAME, as ODESLV_test8 does.
        // var_declared_input() is the generic accessor (ODESLV spells it var_parameter); each input
        // may carry several DOFs, so values go in by parameter_index().
        std::vector<double> Q( S->np(), 0. );
        size_t ndir = 0;
        for( auto const& [w, dom] : S->var_declared_input() ){
          std::string const n = w.name();
          double v = 0.;
          if      ( n == "k" ) v = K;
          else if ( n == "a" ) v = A;
          else { ++ndir; v = 0.; }           // a direction input, seeded by fdiff through fix_input
          for( size_t j : S->parameter_index( w ) ) if( j < Q.size() ) Q[j] = v;
        }
        std::cout << "    inputs: " << S->var_declared_input().size()
                  << " (" << ndir << " direction input(s)), np=" << S->np() << "\n";
        if( S->solve( Q ) != DAESLV_IDAS::NORMAL ){ std::cout << "    solve failed\n"; ++n_fail; }
        else{
          auto const& F = S->val_function();
          std::cout << "    derivative outputs:";
          for( size_t i = 0; i < F.size() && i < 6; ++i ) std::cout << " " << std::scientific
                    << std::setprecision(6) << F[i];
          std::cout << "\n";
          // the directional derivatives appear among the outputs: look for the two we know
          bool got_da = false, got_dk = false;
          for( double v : F ){
            if( std::fabs( v - ( 1.-e )/K ) < 1e-6 ) got_da = true;
            if( std::fabs( v - A*( K*e - (1.-e) )/(K*K) ) < 1e-6 ) got_dk = true;
          }
          check_bool( "d f/da appears among the derivative outputs", got_da );
          check_bool( "d f/dk appears among the derivative outputs", got_dk );
        }
      }
      delete S;
    } }

  std::cout << "\n==========  DAESLV_fdiff: " << n_pass << " passed, " << n_fail << " failed -- "
            << ( n_fail? "FAILURES": "ALL PASS" ) << "  ==========\n\n";
  return n_fail? 1: 0;
}
