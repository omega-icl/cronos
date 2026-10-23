// DAESLV_fsens -- FORWARD sensitivities through IDAS.  The adjoint lives in DAESLV_asens.
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

  // ------------------------------------------------------------------ A
  { std::cout << "\n--- A. x' = -k x, x(0) = a;  sensitivities w.r.t. BOTH parameters\n";
    FFGraph G; DAESLV_IDAS M; M.set_dag( &G );
    FFVar t = G.add_var("t"), x = G.add_var("x(t)"), k = G.add_var("k"), a = G.add_var("a");
    M.add_domain( t, FFDom( 0., 1., 1, FFDom::LGR, 4 ) );
    M.add_state( x, {t} ); M.add_input( k ); M.add_input( a );
    M.set_evolution_domain( t ); M.update_ref( x, 1. );
    M.add_equation( OpP( x, t ) + k*x, {t}, {TI}, io );
    M.add_equation( x - a,             {t}, {FFDom::LB}, ii );
    M.options.SOLVE.RTOL = M.options.SOLVE.ATOL = 1e-11;
    M.options.DISPLAY_LEVEL = 0;
    if( !M.setup() ){ std::cout << "    setup failed: " << M.extract_error() << "\n"; return 1; }
    double const K = 1.3, A = 2.1;
    std::vector<double> P( M.np(), 0. );
    size_t const jk = M.parameter_index( k )[0], ja = M.parameter_index( a )[0];
    P[jk] = K;  P[ja] = A;
    if( M.solve_fsens( P ) != DAESLV_IDAS::NORMAL ){ std::cout << "    solve_fsens failed\n"; return 1; }
    double const xT = A*std::exp(-K);
    check( "x(1) = a e^-k", M.val_state().back()[0], xT, 1e-8 );
    auto const& S = M.val_state_sensitivity();
    std::cout << "    val_state_sensitivity(): " << S.size() << " time(s) x " << S[0].size() << " direction(s)\n";
    check( "d x(0) / da = 1", S.front().at(ja).at(0), 1., 1e-9 );
    check( "d x(0) / dk = 0", S.front().at(jk).at(0), 0., 1e-9 );
    check( "d x(1) / da = e^-k",      S.back().at(ja).at(0),  std::exp(-K), 1e-7 );
    check( "d x(1) / dk = -a e^-k",   S.back().at(jk).at(0), -A*std::exp(-K), 1e-7 ); }

  // ------------------------------------------------------------------ B
  { std::cout << "\n--- B. (2+p) y' + y = 0, 0 = z - p y, y(0) = 1  (mass matrix, algebraic state)\n";
    FFGraph G; DAESLV_IDAS M; M.set_dag( &G );
    FFVar t = G.add_var("t"), y = G.add_var("y(t)"), z = G.add_var("z(t)"), p = G.add_var("p");
    M.add_domain( t, FFDom( 0., 1., 1, FFDom::LGR, 5 ) );
    M.add_state( y, {t} ); M.add_state( z, {t} ); M.add_input( p );
    M.set_evolution_domain( t ); M.update_ref( y, 1. ); M.update_ref( z, 1. );
    M.add_equation( ( 2.+p )*OpP( y, t ) + y, {t}, {TI}, io );
    M.add_equation( z - p*y,                  {t}, {FFDom::ALL}, io );
    M.add_equation( y - 1.,                   {t}, {FFDom::LB}, ii );
    M.options.SOLVE.RTOL = M.options.SOLVE.ATOL = 1e-11;
    M.options.DISPLAY_LEVEL = 0;
    if( !M.setup() ){ std::cout << "    setup failed: " << M.extract_error() << "\n"; return 1; }
    check( "the model is DOF-balanced", M.dof_balance().balanced? 1.: 0., 1., 0. );
    double const PV = 0.7;
    size_t const jp = M.parameter_index( p )[0];
    std::vector<double> P( M.np(), 0. ); P[jp] = PV;
    if( M.solve_fsens( P ) != DAESLV_IDAS::NORMAL ){ std::cout << "    solve_fsens failed\n"; return 1; }
    double const tau = 2.+PV, yT = std::exp( -1./tau ), dyT = yT/( tau*tau );
    check( "y(1) = exp(-1/(2+p))", M.val_state().back()[0], yT, 1e-8 );
    check( "z(1) = p y(1)",        M.val_state().back()[1], PV*yT, 1e-8 );
    auto const& S = M.val_state_sensitivity();
    // at t = 0: y is pinned, so dy/dp = 0 and dz/dp = y(0) = 1
    check( "d y(0) / dp = 0", S.front().at(jp).at(0), 0., 1e-9 );
    check( "d z(0) / dp = 1", S.front().at(jp).at(1), 1., 1e-9 );
    check( "d y(1) / dp = t y /(2+p)^2", S.back().at(jp).at(0), dyT, 1e-7 );
    check( "d z(1) / dp = y + p dy/dp",  S.back().at(jp).at(1), yT + PV*dyT, 1e-7 ); }

  // ------------------------------------------------------------------ C: ACROSS A STOP
  // x' = -k x on [0,2] with a DECLARED JUMP x(1+) = a x(1-).  Exact:
  //     x(1-) = e^-k ,  x(1+) = a e^-k ,  x(2) = a e^-2k
  //     d x(2)/da = e^-2k        (the jump law is linear in a)
  //     d x(2)/dk = -2 a e^-2k
  // The sensitivity must be REINITIALISED at the stop from the same system the state used, or it
  // carries the pre-jump value across and d/da comes out zero.
  { std::cout << "\n--- C. a DECLARED JUMP at t = 1:  x(1+) = a x(1-)\n";
    FFGraph G; DAESLV_IDAS M; M.set_dag( &G );
    FFVar t = G.add_var("t"), x = G.add_var("x(t)"), k = G.add_var("k"), a = G.add_var("a");
    M.add_domain( t, FFDom( 0., 2., 2, FFDom::LGR, 4 ) );          // the boundary at t = 1 is the stop
    M.add_state( x, {t} ); M.add_input( k ); M.add_input( a );
    M.set_evolution_domain( t ); M.update_ref( x, 1. );
    M.add_equation( OpP( x, t ) + k*x, {t}, {TI}, io );
    M.add_equation( x - 1.,            {t}, {FFDom::LB}, ii );
    M.add_transition( a*x, x, t, 1. );
    M.options.SOLVE.RTOL = M.options.SOLVE.ATOL = 1e-11;
    M.options.DISPLAY_LEVEL = 0;
    if( !M.setup() ){ std::cout << "    setup failed: " << M.extract_error() << "\n"; return 1; }
    double const K = 0.9, A = 1.6;
    size_t const jk = M.parameter_index( k )[0], ja = M.parameter_index( a )[0];
    std::vector<double> P( M.np(), 0. ); P[jk] = K; P[ja] = A;
    if( M.solve_fsens( P ) != DAESLV_IDAS::NORMAL ){ std::cout << "    solve_fsens failed\n"; return 1; }
    auto const& S = M.val_state_sensitivity();
    std::cout << "    numIC=" << M.stats.numIC << "   sensitivity rows=" << S.size() << "\n";
    double const e1 = std::exp(-K), e2 = std::exp(-2.*K);
    check( "x(1-) = e^-k",          M.val_state().at(1)[0], e1,    1e-8 );
    check( "x(2)  = a e^-2k",       M.val_state().back()[0], A*e2, 1e-7 );
    // at the LEFT of the stop the jump has not applied yet: d x(1-)/da = 0
    check( "d x(1-) / da = 0",      S.at(1).at(ja).at(0), 0.,        1e-8 );
    check( "d x(1-) / dk = -e^-k",  S.at(1).at(jk).at(0), -e1,       1e-7 );
    check( "d x(2) / da = e^-2k",   S.back().at(ja).at(0), e2,       1e-7 );
    check( "d x(2) / dk = -2a e^-2k", S.back().at(jk).at(0), -2.*A*e2, 1e-7 );

    // The adjoint across this same jump is DAESLV_asens case D.
  }


  // ------------------------------------------------------------------ D: QUADRATURE and OUTPUT gradients
  //   x' = -k x , x(0) = a , on [0,1]
  //   f0 = x(1)                   = a e^-k        df0/da = e^-k ,  df0/dk = -a e^-k
  //   f1 = int_0^1 x dt           = a(1-e^-k)/k   df1/da = (1-e^-k)/k
  //                                               df1/dk = a[ k e^-k - (1-e^-k) ]/k^2
  //   f2 = x(1) + int_0^1 x dt    -- MIXED, so the chain rule must pick up both the state
  //                                  sensitivity and the QUADRATURE sensitivity
  { std::cout << "\n--- D. gradients of a POINT output, a QUADRATURE, and their sum\n";
    FFGraph G; DAESLV_IDAS M; M.set_dag( &G );
    FFEval OpE;  FFIntegral OpI;
    FFVar t = G.add_var("t"), x = G.add_var("x(t)"), k = G.add_var("k"), a = G.add_var("a");
    M.add_domain( t, FFDom( 0., 1., 1, FFDom::LGR, 5 ) );
    M.add_state( x, {t} ); M.add_input( k ); M.add_input( a );
    M.set_evolution_domain( t ); M.update_ref( x, 1. );
    M.add_equation( OpP( x, t ) + k*x, {t}, {TI}, io );
    M.add_equation( x - a,             {t}, {FFDom::LB}, ii );
    M.add_output( OpE( x, t, 1. ) );
    M.add_output( OpI( std::vector<FFVar>{ x }, { t } )[0] );
    M.add_output( OpE( x, t, 1. ) + OpI( std::vector<FFVar>{ x }, { t } )[0] );
    M.options.SOLVE.RTOL = M.options.SOLVE.ATOL = 1e-11;
    M.options.DISPLAY_LEVEL = 0;
    if( !M.setup() ){ std::cout << "    setup failed: " << M.extract_error() << "\n"; return 1; }
    double const K = 1.3, A = 2.1;
    size_t const jk = M.parameter_index( k )[0], ja = M.parameter_index( a )[0];
    std::vector<double> P( M.np(), 0. ); P[jk] = K; P[ja] = A;
    if( M.solve_fsens( P ) != DAESLV_IDAS::NORMAL ){ std::cout << "    solve_fsens failed\n"; return 1; }
    double const e = std::exp(-K), I = A*(1.-e)/K;
    auto const& F = M.val_function();
    auto const& Gd = M.val_function_gradient();
    std::cout << "    nf=" << M.nf() << " nq=" << M.nq() << "  gradient entries=" << Gd.size() << "\n";
    size_t const b0 = M.blk_fct(0).first, b1 = M.blk_fct(1).first, b2 = M.blk_fct(2).first;
    check( "f0 = x(1)",            F.at(b0), A*e,     1e-8 );
    check( "f1 = int x dt",        F.at(b1), I,       1e-7 );
    check( "f2 = f0 + f1",         F.at(b2), A*e + I, 1e-7 );
    check( "d f0 / da = e^-k",     Gd[ja][b0], e,        1e-7 );
    check( "d f0 / dk = -a e^-k",  Gd[jk][b0], -A*e,     1e-7 );
    check( "d f1 / da = (1-e^-k)/k", Gd[ja][b1], (1.-e)/K, 1e-7 );
    check( "d f1 / dk", Gd[jk][b1], A*( K*e - (1.-e) )/(K*K), 1e-7 );
    check( "d f2 / da = sum",      Gd[ja][b2], e + (1.-e)/K, 1e-7 );
    check( "d f2 / dk = sum",      Gd[jk][b2], -A*e + A*( K*e - (1.-e) )/(K*K), 1e-7 ); }

  // ------------------------------------------------------------------ H: INTERMEDIATE point terms
  // Parameter estimation reads the states at measurement times, not just at the end.  Three point
  // outputs at tau = 0.3, 0.6 and 1.0, plus their sum -- the shape of a least-squares residual.
  //   x(tau) = a e^{-k tau} ,  dx/da = e^{-k tau} ,  dx/dk = -tau a e^{-k tau}
  { std::cout << "\n--- H. FORWARD gradients of point outputs at INTERMEDIATE times\n";
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
    M.options.SOLVE.RTOLB = M.options.SOLVE.ATOLB = 1e-11;   // the backward pass has its OWN; left at
    M.options.DISPLAY_LEVEL = 0;                               // default it sat at 1e-07 and looked like a
                                                               // defect in the per-stage sweep -- twice now
    if( !M.setup() ){ std::cout << "    setup failed: " << M.extract_error() << "\n"; return 1; }
    std::cout << "    stages=" << M.ns() << " nf=" << M.nf() << "\n";
    check( "each measurement time made a stage", (double)M.ns(), 3., 0. );
    double const K = 1.3, A = 2.1;
    size_t const jk = M.parameter_index( k )[0], ja = M.parameter_index( a )[0];
    std::vector<double> P( M.np(), 0. ); P[jk] = K; P[ja] = A;
    if( M.solve_fsens( P ) != DAESLV_IDAS::NORMAL ){ std::cout << "    solve_fsens failed\n"; ++n_fail; }
    else{
      auto const& Gd = M.val_function_gradient();
      double sa = 0., sk = 0.;
      for( int m = 0; m < 3; ++m ){
        double const ex = std::exp( -K*TAU[m] );
        size_t const b = M.blk_fct( m ).first;
        check( "d x(" + std::to_string(TAU[m]).substr(0,3) + ") / da", Gd[ja][b], ex, 1e-7 );
        check( "d x(" + std::to_string(TAU[m]).substr(0,3) + ") / dk", Gd[jk][b], -TAU[m]*A*ex, 1e-7 );
        sa += ex;  sk += -TAU[m]*A*ex;
      }
      size_t const b3 = M.blk_fct( 3 ).first;
      check( "d (sum of the three) / da", Gd[ja][b3], sa, 1e-7 );
      check( "d (sum of the three) / dk", Gd[jk][b3], sk, 1e-7 );

    } }

  std::cout << "\n==========  DAESLV_fsens: " << n_pass << " passed, " << n_fail << " failed -- "
            << ( n_fail? "FAILURES": "ALL PASS" ) << "  ==========\n\n";
  return n_fail? 1: 0;
}
