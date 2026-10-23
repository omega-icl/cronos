
// =====================================================================================================
// The sensitivity residual is  F_s,j = dF/dy . s_j + dF/dy' . s'_j + dF/dp_j,  all three blocks SFAD
// entries of the SAME stage plan the state residual uses.  The initial data comes from ic_sens(), the
// same system the state's came from -- so neither the equations nor their initial conditions can drift
// from what they differentiate.
//
//   A. ODE, two parameters:  x' = -k x ,  x(0) = a
//      x(t) = a e^{-k t} ,  dx/da = e^{-k t} ,  dx/dk = -a t e^{-k t}
//
//   B. INDEX-1 DAE with a mass matrix and an algebraic state:
//      (2+p) y' + y = 0 ,  0 = z - p y ,  y(0) = 1
//      y(t) = exp( -t/(2+p) ) ,  dy/dp = t y /(2+p)^2 ,  z = p y ,  dz/dp = y + p dy/dp
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
  std::cout << "    " << std::left << std::setw( 32 ) << w << std::right
            << " got=" << std::scientific << std::setprecision(10) << got
            << " want=" << want << " |d|=" << std::setprecision(2) << d
            << ( d <= tol? "  PASS\n": "  FAIL\n" );
}

int main()
{
  FFPartial OpP;
  int const TI = FFDom::ALL - FFDom::LB;
  auto io = FFModel::EqnOptions( FFModel::EqnRole::INTERIOR, 0 );
  auto ii = FFModel::EqnOptions( FFModel::EqnRole::INITIAL,  0 );

  // ------------------------------------------------------------------ E: ADJOINT vs FORWARD
  //   x' = -k x , x(0) = a ,  f = int_0^1 x dt = a(1-e^-k)/k
  //   df/da = (1-e^-k)/k ,  df/dk = a[ k e^-k - (1-e^-k) ]/k^2
  // The forward gradient is already gated against this closed form to ~1e-11 (case D), so it is a
  // reference the adjoint can be held to as well as the analytic answer.
  { std::cout << "\n--- E. ADJOINT gradient of a quadrature, vs forward and vs the closed form\n";
    FFGraph G; DAESLV_IDAS M; M.set_dag( &G );
    FFIntegral OpI;
    FFVar t = G.add_var("t"), x = G.add_var("x(t)"), k = G.add_var("k"), a = G.add_var("a");
    M.add_domain( t, FFDom( 0., 1., 1, FFDom::LGR, 5 ) );
    M.add_state( x, {t} ); M.add_input( k ); M.add_input( a );
    M.set_evolution_domain( t ); M.update_ref( x, 1. );
    M.add_equation( OpP( x, t ) + k*x, {t}, {TI}, io );
    M.add_equation( x - a,             {t}, {FFDom::LB}, ii );
    M.add_output( OpI( std::vector<FFVar>{ x }, { t } )[0] );
    M.options.SOLVE.RTOL  = M.options.SOLVE.ATOL  = 1e-11;
    M.options.SOLVE.RTOLB = M.options.SOLVE.ATOLB = 1e-11;   // the BACKWARD pass has its own
    M.options.DISPLAY_LEVEL = 0;
    M.options.OUTPUT.RESRECORD = 10;
    if( !M.setup() ){ std::cout << "    setup failed: " << M.extract_error() << "\n"; return 1; }
    double const K = 1.3, A = 2.1, e = std::exp(-K);
    size_t const jk = M.parameter_index( k )[0], ja = M.parameter_index( a )[0];
    std::vector<double> P( M.np(), 0. ); P[jk] = K; P[ja] = A;

    if( M.solve_fsens( P ) != DAESLV_IDAS::NORMAL ){ std::cout << "    solve_fsens failed\n"; return 1; }
    auto const Gf = M.val_function_gradient();

    if( M._solve_asens_one( P, 0 ) != DAESLV_IDAS::NORMAL ){
      std::cout << "    solve_asens failed: " << M.extract_error() << "\n"; ++n_fail; }
    else{
      auto const& Ga = M.val_function_gradient();
      size_t const b = M.blk_fct(0).first;
      std::cout << "    adjoint  d/da=" << std::scientific << std::setprecision(10) << Ga[ja][b]
                << "   d/dk=" << Ga[jk][b] << "\n";
      check( "adjoint d f/da = (1-e^-k)/k", Ga[ja][b], (1.-e)/K, 1e-7 );
      check( "adjoint d f/dk",              Ga[jk][b], A*( K*e - (1.-e) )/(K*K), 1e-7 );
      check( "adjoint == forward on d/da",  Ga[ja][b], Gf[ja][b], 1e-7 );
      check( "adjoint == forward on d/dk",  Ga[jk][b], Gf[jk][b], 1e-7 );
    } }

  // ------------------------------------------------------------------ E: ADJOINT vs FORWARD
  //   x' = -k x , x(0) = a ,  f = int_0^1 x dt = a(1-e^-k)/k
  //   df/da = (1-e^-k)/k ,  df/dk = a[ k e^-k - (1-e^-k) ]/k^2
  // The forward gradient is already gated against this closed form to ~1e-11 (case D), so it is a
  // reference the adjoint can be held to as well as the analytic answer.
  { std::cout << "\n--- E. ADJOINT gradient of a quadrature, vs forward and vs the closed form\n";
    FFGraph G; DAESLV_IDAS M; M.set_dag( &G );
    FFIntegral OpI;
    FFVar t = G.add_var("t"), x = G.add_var("x(t)"), k = G.add_var("k"), a = G.add_var("a");
    M.add_domain( t, FFDom( 0., 1., 1, FFDom::LGR, 5 ) );
    M.add_state( x, {t} ); M.add_input( k ); M.add_input( a );
    M.set_evolution_domain( t ); M.update_ref( x, 1. );
    M.add_equation( OpP( x, t ) + k*x, {t}, {TI}, io );
    M.add_equation( x - a,             {t}, {FFDom::LB}, ii );
    M.add_output( OpI( std::vector<FFVar>{ x }, { t } )[0] );
    M.options.SOLVE.RTOL  = M.options.SOLVE.ATOL  = 1e-11;
    M.options.SOLVE.RTOLB = M.options.SOLVE.ATOLB = 1e-11;   // the BACKWARD pass has its own
    M.options.DISPLAY_LEVEL = 0;
    if( !M.setup() ){ std::cout << "    setup failed: " << M.extract_error() << "\n"; return 1; }
    double const K = 1.3, A = 2.1, e = std::exp(-K);
    size_t const jk = M.parameter_index( k )[0], ja = M.parameter_index( a )[0];
    std::vector<double> P( M.np(), 0. ); P[jk] = K; P[ja] = A;

    if( M.solve_fsens( P ) != DAESLV_IDAS::NORMAL ){ std::cout << "    solve_fsens failed\n"; return 1; }
    auto const Gf = M.val_function_gradient();

    if( M._solve_asens_one( P, 0 ) != DAESLV_IDAS::NORMAL ){
      std::cout << "    solve_asens failed: " << M.extract_error() << "\n"; ++n_fail; }
    else{
      auto const& Ga = M.val_function_gradient();
      size_t const b = M.blk_fct(0).first;
      std::cout << "    adjoint  d/da=" << std::scientific << std::setprecision(10) << Ga[ja][b]
                << "   d/dk=" << Ga[jk][b] << "\n";
      check( "adjoint d f/da = (1-e^-k)/k", Ga[ja][b], (1.-e)/K, 1e-7 );
      check( "adjoint d f/dk",              Ga[jk][b], A*( K*e - (1.-e) )/(K*K), 1e-7 );
      check( "adjoint == forward on d/da",  Ga[ja][b], Gf[ja][b], 1e-7 );
      check( "adjoint == forward on d/dk",  Ga[jk][b], Gf[jk][b], 1e-7 );
    } }

  // ------------------------------------------------------------------ G: adjoint of a POINT output
  //   f = x(1) = a e^-k  -- no running integrand at all, so the whole influence enters through the
  //   TERMINAL condition  Jyp^T lambda(T) = ( df/dx )( T ).  df/da = e^-k, df/dk = -a e^-k.
  { std::cout << "\n--- G. ADJOINT of a POINT output (terminal seed, no running term)\n";
    FFGraph G; DAESLV_IDAS M; M.set_dag( &G );
    FFEval OpE;  FFIntegral OpI;
    FFVar t = G.add_var("t"), x = G.add_var("x(t)"), k = G.add_var("k"), a = G.add_var("a");
    M.add_domain( t, FFDom( 0., 1., 1, FFDom::LGR, 5 ) );
    M.add_state( x, {t} ); M.add_input( k ); M.add_input( a );
    M.set_evolution_domain( t ); M.update_ref( x, 1. );
    M.add_equation( OpP( x, t ) + k*x, {t}, {TI}, io );
    M.add_equation( x - a,             {t}, {FFDom::LB}, ii );
    M.add_output( OpE( x, t, 1. ) );
    M.add_output( OpE( x, t, 1. ) + OpI( std::vector<FFVar>{ x }, { t } )[0] );   // MIXED
    M.options.SOLVE.RTOL  = M.options.SOLVE.ATOL  = 1e-11;
    M.options.SOLVE.RTOLB = M.options.SOLVE.ATOLB = 1e-11;
    M.options.DISPLAY_LEVEL = 0;
    if( !M.setup() ){ std::cout << "    setup failed: " << M.extract_error() << "\n"; return 1; }
    double const K = 1.3, A = 2.1, e = std::exp(-K);
    size_t const jk = M.parameter_index( k )[0], ja = M.parameter_index( a )[0];
    std::vector<double> P( M.np(), 0. ); P[jk] = K; P[ja] = A;
    M.solve_fsens( P );
    auto const Gf = M.val_function_gradient();
    if( M._solve_asens_one( P, 0 ) != DAESLV_IDAS::NORMAL ){
      std::cout << "    solve_asens failed: " << M.extract_error() << "\n"; ++n_fail; }
    else{
      auto const& Ga = M.val_function_gradient();
      size_t const b = M.blk_fct(0).first;
      check( "point-output adjoint d f/da = e^-k",    Ga[ja][b], e,    1e-7 );
      check( "point-output adjoint d f/dk = -a e^-k", Ga[jk][b], -A*e, 1e-7 );
      check( "  == forward on d/da", Ga[ja][b], Gf[ja][b], 1e-7 );
      check( "  == forward on d/dk", Ga[jk][b], Gf[jk][b], 1e-7 );
    }
    // the MIXED output: a terminal seed AND a running term, superposed
    if( M._solve_asens_one( P, 1 ) != DAESLV_IDAS::NORMAL ){
      std::cout << "    mixed solve_asens failed: " << M.extract_error() << "\n"; ++n_fail; }
    else{
      auto const& Gm = M.val_function_gradient();
      size_t const b1 = M.blk_fct(1).first;
      check( "MIXED adjoint d f/da", Gm[ja][b1], e + (1.-e)/K, 1e-7 );
      check( "MIXED adjoint d f/dk", Gm[jk][b1], -A*e + A*( K*e - (1.-e) )/(K*K), 1e-7 );
      check( "  == forward on d/da", Gm[ja][b1], Gf[ja][b1], 1e-7 );
      check( "  == forward on d/dk", Gm[jk][b1], Gf[jk][b1], 1e-7 );
    } }

  // ------------------------------------------------------------------ D: across a DECLARED JUMP
  // x' = -k x on [0,2] with x(1+) = a x(1-), so x(2) = a e^-2k.  The co-state crosses the stop by
  // the TRANSPOSE of the system the forward solved there (ic_adjoint).
  { std::cout << "\n--- D. ADJOINT across a DECLARED JUMP\n";
    FFGraph G; DAESLV_IDAS M; M.set_dag( &G );
    FFEval OpE;
    FFVar t = G.add_var("t"), x = G.add_var("x(t)"), k = G.add_var("k"), a = G.add_var("a");
    M.add_domain( t, FFDom( 0., 2., 2, FFDom::LGR, 4 ) );
    M.add_state( x, {t} ); M.add_input( k ); M.add_input( a );
    M.set_evolution_domain( t ); M.update_ref( x, 1. );
    M.add_equation( OpP( x, t ) + k*x, {t}, {TI}, io );
    M.add_equation( x - 1.,            {t}, {FFDom::LB}, ii );
    M.add_transition( a*x, x, t, 1. );
    M.add_output( OpE( x, t, 2. ) );
    M.options.SOLVE.RTOL  = M.options.SOLVE.ATOL  = 1e-11;
    M.options.SOLVE.RTOLB = M.options.SOLVE.ATOLB = 1e-11;
    M.options.DISPLAY_LEVEL = 0;
    if( !M.setup() ){ std::cout << "    setup failed: " << M.extract_error() << "\n"; return 1; }
    double const K = 0.9, A = 1.6, e2 = std::exp(-2.*K);
    size_t const jk = M.parameter_index( k )[0], ja = M.parameter_index( a )[0];
    std::vector<double> P( M.np(), 0. ); P[jk] = K; P[ja] = A;
    M.solve_fsens( P );
    auto const Gf = M.val_function_gradient();
    if( M._solve_asens_one( P, 0 ) != DAESLV_IDAS::NORMAL ){
      std::cout << "    solve_asens failed: " << M.extract_error() << "\n"; ++n_fail; }
    else{
      auto const& Ga = M.val_function_gradient();  size_t const b = M.blk_fct(0).first;
      check( "d x(2) / da = e^-2k",     Ga[ja][b], e2,       1e-7 );
      check( "d x(2) / dk = -2a e^-2k", Ga[jk][b], -2.*A*e2, 1e-7 );
      check( "  == forward on d/da",    Ga[ja][b], Gf[ja][b], 1e-7 );
      check( "  == forward on d/dk",    Ga[jk][b], Gf[jk][b], 1e-7 );
    } }

  // ------------------------------------------------------------------ E: INTERIOR measurement times
  // The parameter-estimation shape: measurements at 0.3, 0.6 and 1.0, each an IMPULSE on lambda at
  // its own stage time, plus their sum whose terms sit at three stage times at once.
  { std::cout << "\n--- E. ADJOINT at INTERIOR measurement times (time-series data)\n";
    FFGraph G; DAESLV_IDAS M; M.set_dag( &G );
    FFEval OpE;
    FFVar t = G.add_var("t"), x = G.add_var("x(t)"), k = G.add_var("k"), a = G.add_var("a");
    M.add_domain( t, FFDom( 0., 1., 1, FFDom::LGR, 5 ) );
    M.add_state( x, {t} ); M.add_input( k ); M.add_input( a );
    M.set_evolution_domain( t ); M.update_ref( x, 1. );
    M.add_equation( OpP( x, t ) + k*x, {t}, {TI}, io );
    M.add_equation( x - a,             {t}, {FFDom::LB}, ii );
    double const TAU[3] = { 0.3, 0.6, 1.0 };
    for( double tv : TAU ) M.add_output( OpE( x, t, tv ) );
    M.add_output( OpE( x, t, TAU[0] ) + OpE( x, t, TAU[1] ) + OpE( x, t, TAU[2] ) );
    M.options.SOLVE.RTOL  = M.options.SOLVE.ATOL  = 1e-11;
    M.options.SOLVE.RTOLB = M.options.SOLVE.ATOLB = 1e-11;
    M.options.DISPLAY_LEVEL = 0;
    if( !M.setup() ){ std::cout << "    setup failed: " << M.extract_error() << "\n"; return 1; }
    double const K = 1.3, A = 2.1;
    size_t const jk = M.parameter_index( k )[0], ja = M.parameter_index( a )[0];
    std::vector<double> P( M.np(), 0. ); P[jk] = K; P[ja] = A;
    M.solve_fsens( P );
    auto const Gf = M.val_function_gradient();
    double sa = 0., sk = 0.;
    for( int m = 0; m < 3; ++m ){
      double const ex = std::exp( -K*TAU[m] );
      sa += ex;  sk += -TAU[m]*A*ex;
      if( M._solve_asens_one( P, (size_t)m ) != DAESLV_IDAS::NORMAL ){
        std::cout << "    adjoint at tau=" << TAU[m] << " failed\n"; ++n_fail; continue; }
      auto const& Ga = M.val_function_gradient();  size_t const b = M.blk_fct( m ).first;
      check( "d x(" + std::to_string(TAU[m]).substr(0,3) + ") / da", Ga[ja][b], ex, 1e-7 );
      check( "d x(" + std::to_string(TAU[m]).substr(0,3) + ") / dk", Ga[jk][b], -TAU[m]*A*ex, 1e-7 );
    }
    if( M._solve_asens_one( P, 3 ) != DAESLV_IDAS::NORMAL ){
      std::cout << "    adjoint of the sum failed\n"; ++n_fail; }
    else{
      auto const& Ga = M.val_function_gradient();  size_t const b3 = M.blk_fct(3).first;
      check( "d (sum of the three) / da", Ga[ja][b3], sa, 1e-7 );
      check( "d (sum of the three) / dk", Ga[jk][b3], sk, 1e-7 );
      check( "  == forward on d/da", Ga[ja][b3], Gf[ja][b3], 1e-7 );
      check( "  == forward on d/dk", Ga[jk][b3], Gf[jk][b3], 1e-7 );
    } }

  std::cout << "\n==========  DAESLV_asens: " << n_pass << " passed, " << n_fail << " failed -- "
            << ( n_fail? "FAILURES": "ALL PASS" ) << "  ==========\n\n";
  return n_fail? 1: 0;
}
