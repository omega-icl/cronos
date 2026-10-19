// DAESLV_outputs -- INTEGRAL captures carried as quadratures by IDAS.
// =====================================================================================================
//  A. plain ODE, integrand of the state:   x' = -x, x(0)=1,  Q = int_0^1 x dt = 1 - e^{-1}
//  B. MASS MATRIX, so the quadrature is checked where the rows are NOT in divided form:
//        (2+t) x' + x = 0 , x(0)=1  =>  x(t) = 2/(2+t) ,  Q = int_0^2 x dt = 2 ln 2
//  C. ACROSS RESTARTS: u piecewise constant over 4 elements, Q = int_0^2 u dt, exact by construction.
//     A quadrature accumulates from t0 and must CARRY across a stage boundary, not restart at it.
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
  std::cout << "    " << std::left << std::setw( 44 ) << w << std::right
            << " got=" << std::scientific << std::setprecision(10) << got
            << " want=" << want << " |d|=" << std::setprecision(2) << d
            << ( d <= tol? "  PASS\n": "  FAIL\n" );
}

int main()
{
  FFPartial OpP;  FFEval OpE;  FFIntegral OpI;
  int const TI = FFDom::ALL - FFDom::LB;

  { std::cout << "\n--- A. x' = -x, x(0) = 1;  Q = int_0^1 x dt = 1 - e^-1\n";
    FFGraph G; DAESLV_IDAS M; M.set_dag( &G );
    FFVar t = G.add_var("t"), x = G.add_var("x(t)");
    M.add_domain( t, FFDom( 0., 1., 1, FFDom::LGR, 4 ) );
    M.add_state( x, {t} ); M.set_evolution_domain( t ); M.update_ref( x, 1. );
    auto io = FFModel::EqnOptions( FFModel::EqnRole::INTERIOR, 0 );
    auto ii = FFModel::EqnOptions( FFModel::EqnRole::INITIAL,  0 );
    M.add_equation( OpP( x, t ) + x, {t}, {TI}, io );
    M.add_equation( x - 1.,          {t}, {FFDom::LB}, ii );
    M.add_output( OpI( std::vector<FFVar>{ x }, { t } )[0] );
    M.options.SOLVE.RTOL = M.options.SOLVE.ATOL = 1e-10;
    M.options.DISPLAY_LEVEL = 0;
    if( !M.setup() ){ std::cout << "    setup failed: " << M.extract_error() << "\n"; return 1; }
    std::cout << "    nq=" << M.nq() << " nf=" << M.nf() << "\n";
    check( "the INTEGRAL became a quadrature", (double)M.nq(), 1., 0. );
    M.solve( {} );
    check( "Q(1) = 1 - e^-1", M.val_quad().back()[0], 1. - std::exp(-1.), 1e-8 );
    check( "Q(0) = 0",        M.val_quad().front()[0], 0.,                0. ); }

  { std::cout << "\n--- B. (2+t) x' + x = 0, x(0) = 1;  Q = int_0^2 x dt = 2 ln 2\n";
    FFGraph G; DAESLV_IDAS M; M.set_dag( &G );
    FFVar t = G.add_var("t"), x = G.add_var("x(t)");
    M.add_domain( t, FFDom( 0., 2., 1, FFDom::LGR, 5 ) );
    M.add_state( x, {t} ); M.set_evolution_domain( t ); M.update_ref( x, 1. );
    auto io = FFModel::EqnOptions( FFModel::EqnRole::INTERIOR, 0 );
    auto ii = FFModel::EqnOptions( FFModel::EqnRole::INITIAL,  0 );
    M.add_equation( ( 2. + t ) * OpP( x, t ) + x, {t}, {TI}, io );   // NON-TRIVIAL mass coefficient
    M.add_equation( x - 1.,                       {t}, {FFDom::LB}, ii );
    M.add_output( OpI( std::vector<FFVar>{ x }, { t } )[0] );
    M.options.SOLVE.RTOL = M.options.SOLVE.ATOL = 1e-10;
    M.options.DISPLAY_LEVEL = 0;
    if( !M.setup() ){ std::cout << "    setup failed: " << M.extract_error() << "\n"; return 1; }
    M.solve( {} );
    check( "x(2) = 1/2",        M.val_state().back()[0], 0.5,                1e-8 );
    check( "Q(2) = 2 ln 2",     M.val_quad().back()[0],  2.*std::log(2.),    1e-7 ); }

  { std::cout << "\n--- C. u piecewise constant over 4 elements;  Q = int_0^2 u dt, across 3 restarts\n";
    size_t const NE = 4;
    FFGraph G; DAESLV_IDAS M; M.set_dag( &G );
    FFVar t = G.add_var("t"), x = G.add_var("x(t)"), u = G.add_var("u");
    M.add_domain( t, FFDom( 0., 2., NE, FFDom::LGR, 4 ) );
    M.add_state( x, {t} ); M.add_input( u, {t}, FFDom::LGR, 1 );
    M.set_evolution_domain( t ); M.update_ref( x, 1. );
    auto io = FFModel::EqnOptions( FFModel::EqnRole::INTERIOR, 0 );
    auto ii = FFModel::EqnOptions( FFModel::EqnRole::INITIAL,  0 );
    M.add_equation( OpP( x, t ) + x - u, {t}, {TI}, io );
    M.add_equation( x - 1.,              {t}, {FFDom::LB}, ii );
    M.add_output( OpI( std::vector<FFVar>{ u }, { t } )[0] );
    M.options.SOLVE.RTOL = M.options.SOLVE.ATOL = 1e-10;
    M.options.DISPLAY_LEVEL = 0;
    if( !M.setup() ){ std::cout << "    setup failed: " << M.extract_error() << "\n"; return 1; }
    std::vector<double> P;  double exact = 0.;
    for( size_t e = 0; e < NE; ++e ){ P.push_back( 1. + (double)e ); exact += P.back() * 2./(double)NE; }
    M.solve( P );
    std::cout << "    numIC=" << M.stats.numIC << "  rows of val_quad()=" << M.val_quad().size() << "\n";
    // The OUTPUT is the sum over stages, which is val_function()'s job; val_quad() holds the per-stage
    // increments it is summed from.  Checking both pins the division of labour, not just the total.
    check( "the output integral over 4 stages", M.val_function()[ M.blk_fct(0).first ], exact, 1e-7 );
    for( size_t e = 0; e < NE; ++e )
      check( "  stage " + std::to_string(e+1) + " increment = u_e h",
             M.val_quad().at(e+1)[0], P[e] * 2./(double)NE, 1e-7 );
    M.solve( P );
    check( "solve() idempotent on the output", M.val_function()[ M.blk_fct(0).first ], exact, 1e-7 ); }

  { std::cout << "\n--- D. POINT captures: a value at the end, at an INTERIOR tau, and mixed with an integral\n";
    FFGraph G; DAESLV_IDAS M; M.set_dag( &G );
    FFVar t = G.add_var("t"), x = G.add_var("x(t)");
    M.add_domain( t, FFDom( 0., 1., 1, FFDom::LGR, 4 ) );
    M.add_state( x, {t} ); M.set_evolution_domain( t ); M.update_ref( x, 1. );
    auto io = FFModel::EqnOptions( FFModel::EqnRole::INTERIOR, 0 );
    auto ii = FFModel::EqnOptions( FFModel::EqnRole::INITIAL,  0 );
    M.add_equation( OpP( x, t ) + x, {t}, {TI}, io );
    M.add_equation( x - 1.,          {t}, {FFDom::LB}, ii );
    M.add_output( OpE( x, t, 1.  ) );                                  // x(1) = e^-1
    M.add_output( OpE( x, t, 0.4 ) );                                  // x(0.4) = e^-0.4, an INTERIOR tau
    M.add_output( OpE( x, t, 1. ) + OpI( std::vector<FFVar>{ x }, { t } )[0] );  // mixed POINT + INTEGRAL
    M.add_output( OpI( std::vector<FFVar>{ x }, { t } )[0] );          // the SAME integral, on its own
    M.options.SOLVE.RTOL = M.options.SOLVE.ATOL = 1e-10;
    M.options.DISPLAY_LEVEL = 0;
    if( !M.setup() ){ std::cout << "    setup failed: " << M.extract_error() << "\n"; return 1; }
    std::cout << "    nq=" << M.nq() << " nf=" << M.nf() << " stages=" << M.ns() << "\n";
    check( "an interior tau created its own stage", (double)M.ns(), 2., 0. );
    M.solve( {} );
    auto const& F = M.val_function();
    if( F.size() < 3 ){ std::cout << "    val_function() is " << F.size() << " long\n"; ++n_fail; }
    else{
      check( "f0 = x(1) = e^-1",            F[ M.blk_fct(0).first ], std::exp(-1.),  1e-8 );
      check( "f1 = x(0.4) = e^-0.4",        F[ M.blk_fct(1).first ], std::exp(-0.4), 1e-8 );
      check( "f2 = x(1) + int_0^1 x dt",    F[ M.blk_fct(2).first ],
             std::exp(-1.) + 1. - std::exp(-1.), 1e-7 );
      check( "f3 = int_0^1 x dt alone",     F[ M.blk_fct(3).first ], 1. - std::exp(-1.), 1e-7 );
      std::cout << "    stage times / per-stage quadrature increments:\n";
      for( size_t k = 0; k < M.val_quad().size(); ++k )
        std::cout << "      k=" << k << "  q=" << std::scientific << std::setprecision(10)
                  << M.val_quad()[k][0] << "\n";
      std::cout << "      exact stage 1 int_0^0.4 = " << 1.-std::exp(-0.4)
                << "   stage 2 int_0.4^1 = " << std::exp(-0.4)-std::exp(-1.) << "\n";
    }
    M.solve( {} );
    check( "solve() idempotent on f0", M.val_function()[ M.blk_fct(0).first ], std::exp(-1.), 1e-8 ); }

  std::cout << "\n==========  DAESLV_outputs: " << n_pass << " passed, " << n_fail << " failed -- "
            << ( n_fail? "FAILURES": "ALL PASS" ) << "  ==========\n\n";
  return n_fail? 1: 0;
}
