// DAESLV_inputs -- constants, and the two kinds of input, against closed forms.
// =====================================================================================================
//  A  CONSTANT          y' = -k y , y(0)=1 , k a constant      -> y(T) = exp(-k T)
//     plus the refusal when the wrong number of constant values is handed to solve().
//  B  LUMPED input      the same model with k an input carrying ONE level: nothing varies along t, so
//     the rows are the same in every stage and the integrator must NOT restart (numIC == 1).
//  C  DISTRIBUTED input y' = -y + u with u = t, piecewise LINEAR (2 nodes per element, so u = t is
//     EXACT)  -> y = t - 1 + exp(-t).  The other drivers only ever use one node per element, so the
//     multi-node level packing is exercised here and nowhere else.
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
static void check_bool( std::string const& w, bool ok )
{
  ok? ++n_pass: ++n_fail;
  std::cout << "    " << std::left << std::setw( 44 ) << w << std::right
            << ( ok? "                                      PASS\n": "                                      FAIL\n" );
}

int main()
{
  double const K = 0.7, T_END = 2.0;

  // ---------------------------------------------------------------- A  constant
  { std::cout << "\n--- A. CONSTANT k in the residual: y' = -k y, y(0) = 1, k = " << K << "\n";
    FFGraph G; DAESLV_IDAS M; M.set_dag( &G );
    FFVar t=G.add_var("t"), y=G.add_var("y(t)"), k=G.add_var("k");
    FFPartial OpP; int const TI = FFDom::ALL - FFDom::LB;
    M.add_domain( t, FFDom( 0., T_END, 2, FFDom::LGR, 4 ) );
    M.add_state( y, {t} );
    M.set_constant( { k } );
    M.set_evolution_domain( t ); M.update_ref( y, 1. );
    auto io=FFModel::EqnOptions(FFModel::EqnRole::INTERIOR,0), ii=FFModel::EqnOptions(FFModel::EqnRole::INITIAL,0);
    M.add_equation( OpP(y,t) + k*y, {t}, {TI},        io );
    M.add_equation( y - 1.,         {t}, {FFDom::LB}, ii );
    M.options.DISPLAY_LEVEL = 0;
    M.options.SOLVE.RTOL = M.options.SOLVE.ATOL = 1e-10;
    if( !M.setup() ){ std::cout << "    setup failed: " << M.extract_error() << "\n"; ++n_fail; }
    else{
      std::cout << "    nc=" << M.nc() << " np=" << M.np() << "\n";
      if( M.solve( {}, { K } ) != DAESLV_IDAS::NORMAL ){
        std::cout << "    solve failed: " << M.extract_error() << "\n"; ++n_fail; }
      else
        check( "y(T) = exp(-k T)", M.val_state().back()[0], std::exp( -K*T_END ), 1e-8 );
      // the POINTER overload, which nothing else runs
      double const cp[1] = { K };
      if( M.solve( nullptr, cp ) != DAESLV_IDAS::NORMAL ){
        std::cout << "    pointer-form solve failed: " << M.extract_error() << "\n"; ++n_fail; }
      else
        check( "pointer-form solve agrees", M.val_state().back()[0], std::exp( -K*T_END ), 1e-8 );
      // and the refusal: a wrong number of constant values must not be integrated through
      check_bool( "wrong-length c is REFUSED",
                  M.solve( {}, { K, K } ) == DAESLV_IDAS::FATAL );
    } }

  // ---------------------------------------------------------------- B  lumped input
  { std::cout << "\n--- B. LUMPED input k (one level, nothing varies along t)\n";
    FFGraph G; DAESLV_IDAS M; M.set_dag( &G );
    FFVar t=G.add_var("t"), y=G.add_var("y(t)"), k=G.add_var("k");
    FFPartial OpP; int const TI = FFDom::ALL - FFDom::LB;
    M.add_domain( t, FFDom( 0., T_END, 2, FFDom::LGR, 4 ) );
    M.add_state( y, {t} );
    M.add_input( k );                                   // LUMPED: no domain
    M.set_evolution_domain( t ); M.update_ref( y, 1. );
    auto io=FFModel::EqnOptions(FFModel::EqnRole::INTERIOR,0), ii=FFModel::EqnOptions(FFModel::EqnRole::INITIAL,0);
    M.add_equation( OpP(y,t) + k*y, {t}, {TI},        io );
    M.add_equation( y - 1.,         {t}, {FFDom::LB}, ii );
    M.options.DISPLAY_LEVEL = 0;
    M.options.SOLVE.RTOL = M.options.SOLVE.ATOL = 1e-10;
    if( !M.setup() ){ std::cout << "    setup failed: " << M.extract_error() << "\n"; ++n_fail; }
    else{
      std::cout << "    np=" << M.np() << " input levels=" << M.n_input_levels() << "\n";
      if( M.solve( std::vector<double>{ K } ) != DAESLV_IDAS::NORMAL ){
        std::cout << "    solve failed: " << M.extract_error() << "\n"; ++n_fail; }
      else{
        check( "y(T) = exp(-k T), k as a lumped input", M.val_state().back()[0], std::exp( -K*T_END ), 1e-8 );
        check_bool( "no restart: the rows never change (numIC == 1)", M.stats.numIC == 1 );
      } } }

  // ---------------------------------------------------------------- C  distributed, piecewise linear
  { std::cout << "\n--- C. DISTRIBUTED input, piecewise LINEAR: y' = -y + u, u = t, y(0) = 0\n";
    FFGraph G; DAESLV_IDAS M; M.set_dag( &G );
    FFVar t=G.add_var("t"), y=G.add_var("y(t)"), u=G.add_var("u(t)");
    FFPartial OpP; int const TI = FFDom::ALL - FFDom::LB;
    M.add_domain( t, FFDom( 0., T_END, 2, FFDom::LGR, 4 ) );
    M.add_state( y, {t} );
    M.add_input( u, {t}, FFDom::LGL, 2 );               // TWO nodes per element -> piecewise linear
    M.set_evolution_domain( t ); M.update_ref( y, 0. );
    auto io=FFModel::EqnOptions(FFModel::EqnRole::INTERIOR,0), ii=FFModel::EqnOptions(FFModel::EqnRole::INITIAL,0);
    M.add_equation( OpP(y,t) + y - u, {t}, {TI},        io );
    M.add_equation( y - 0.,           {t}, {FFDom::LB}, ii );
    M.options.DISPLAY_LEVEL = 0;
    M.options.SOLVE.RTOL = M.options.SOLVE.ATOL = 1e-10;
    if( !M.setup() ){ std::cout << "    setup failed: " << M.extract_error() << "\n"; ++n_fail; }
    else{
      size_t const nl = M.n_input_levels();
      std::cout << "    input levels=" << nl << " (2 elements x 2 nodes)\n";
      // u = t sampled at the element nodes: 0,1 on [0,1] and 1,2 on [1,2], shared or not
      std::vector<double> lv;
      if( nl == 4 ) lv = { 0., 1., 1., 2. };
      else if( nl == 3 ) lv = { 0., 1., 2. };
      else { std::cout << "    unexpected level count\n"; ++n_fail; }
      if( !lv.empty() ){
        if( M.solve( lv, std::vector<double>{} ) != DAESLV_IDAS::NORMAL ){
          std::cout << "    solve failed: " << M.extract_error() << "\n"; ++n_fail; }
        else
          check( "y(T) = T - 1 + exp(-T)", M.val_state().back()[0],
                 T_END - 1. + std::exp( -T_END ), 1e-7 );
      } } }

  // ---------------------------------------------------------------- D  implicit INITIAL rows
  { std::cout << "\n--- D. IMPLICIT initial rows, resolved by the KINSOL corrector\n";
    // x0' = -x0 , x1' = -2 x1 , with the initial conditions given IMPLICITLY:
    //     x0(0) + x1(0) - 3 = 0   and   x0(0) - x1(0) - 1 = 0   ->   x0(0)=2 , x1(0)=1
    // Each row split on its own state leaves the OTHER in the remainder, so ic[x0] = 3 - x1 and
    // ic[x1] = x0 - 1: one evaluation cannot resolve them.  ODESLV refuses such a model; DAESLV
    // solves the rows.  Closed form x0 = 2 e^-t , x1 = e^-2t .
    FFGraph G; DAESLV_IDAS M; M.set_dag( &G );
    FFVar t=G.add_var("t"), x0=G.add_var("x0(t)"), x1=G.add_var("x1(t)");
    FFPartial OpP; int const TI = FFDom::ALL - FFDom::LB;
    M.add_domain( t, FFDom( 0., 1., 2, FFDom::LGR, 4 ) );
    M.add_state( x0, {t} ); M.add_state( x1, {t} );
    M.set_evolution_domain( t ); M.update_ref( x0, 0.5 ); M.update_ref( x1, 0.5 );  // NOT consistent
    auto io=FFModel::EqnOptions(FFModel::EqnRole::INTERIOR,0), ii=FFModel::EqnOptions(FFModel::EqnRole::INITIAL,0);
    M.add_equation( OpP(x0,t) + x0,    {t}, {TI},        io );
    M.add_equation( OpP(x1,t) + 2.*x1, {t}, {TI},        io );
    M.add_equation( x0 + x1 - 3.,      {t}, {FFDom::LB}, ii );
    M.add_equation( x0 - x1 - 1.,      {t}, {FFDom::LB}, ii );
    M.options.DISPLAY_LEVEL = 0;
    M.options.SOLVE.RTOL = M.options.SOLVE.ATOL = 1e-10;
    M.options.INITIAL.SOLVER = BASE_IDAS::Options::t_Initial::KINSOL;
    check_bool( "implicit INITIAL rows are ACCEPTED", M.setup() );
    if( M.setup() ){
      std::vector<double> y, yp;
      auto const R = M.solve_ic( y, yp );
      std::cout << "    solve_ic: " << R.message << " (" << R.iterations << " iterations)\n";
      check_bool( "the IC system solves", R.status == BASE_KINSOL::OK );
      check( "x0(0) = 2 recovered from the implicit rows", y[0], 2., 1e-12 );
      check( "x1(0) = 1 recovered from the implicit rows", y[1], 1., 1e-12 );
      if( M.solve( std::vector<double>{} ) != DAESLV_IDAS::NORMAL ){
        std::cout << "    solve failed: " << M.extract_error() << "\n"; ++n_fail; }
      else{
        check( "x0(1) = 2 e^-1", M.val_state().back()[0], 2.*std::exp(-1.), 1e-8 );
        check( "x1(1) = e^-2",   M.val_state().back()[1], std::exp(-2.),    1e-8 );
      } } }

  // ---------------------------------------------------------------- E  general h(x, y, dx/dt)
  { std::cout << "\n--- E. a GENERAL initial condition h( x(0), y(0), dx/dt(0) ) = 0\n";
    // x' = -x z  (differential x) , 0 = z - x^2  (algebraic z) , and ONE initial row coupling all three:
    //     x(0) + 2 z(0) + dx/dt(0) - 2 = 0
    // With x(0)=a the DAE forces z=a^2 and dx/dt=-a^3, so the row reads a + 2a^2 - a^3 = 2, root a = 1,
    // giving x(0)=1, z(0)=1, dx/dt(0)=-1.  Nothing here is solvable row by row: the algebraic state, the
    // derivative and the differential state are determined TOGETHER, which is the whole point of posing
    // the consistent-initialisation problem jointly over (x, y, dx/dt).  The coefficients matter: with
    // x + z + dx/dt - 1 the row's gradient at the solution is [1,1,1], IDENTICAL to the first interior
    // row's, and the joint Jacobian is singular -- which the solver correctly diagnosed when this case
    // was first written that way.
    FFGraph G; DAESLV_IDAS M; M.set_dag( &G );
    FFVar t=G.add_var("t"), x=G.add_var("x(t)"), z=G.add_var("z(t)");
    FFPartial OpP; int const TI = FFDom::ALL - FFDom::LB;
    M.add_domain( t, FFDom( 0., 1., 2, FFDom::LGR, 4 ) );
    M.add_state( x, {t} ); M.add_state( z, {t} );
    M.set_evolution_domain( t ); M.update_ref( x, 0.4 ); M.update_ref( z, 0.2 );   // NOT consistent
    auto io=FFModel::EqnOptions(FFModel::EqnRole::INTERIOR,0), ii=FFModel::EqnOptions(FFModel::EqnRole::INITIAL,0);
    M.add_equation( OpP(x,t) + x*z,          {t}, {TI},         io );
    M.add_equation( z - x*x,                 {t}, {FFDom::ALL}, io );
    M.add_equation( x + 2.*z + OpP(x,t) - 2., {t}, {FFDom::LB}, ii );   // h(x, y, dx/dt)
    M.options.DISPLAY_LEVEL = 0;
    M.options.SOLVE.RTOL = M.options.SOLVE.ATOL = 1e-10;
    M.options.INITIAL.SOLVER = BASE_IDAS::Options::t_Initial::KINSOL;
    check_bool( "a row in x, y AND dx/dt is accepted", M.setup() );
    if( M.setup() ){
      std::vector<double> y, yp;
      auto const R = M.solve_ic( y, yp );
      std::cout << "    solve_ic: " << R.message << " (" << R.iterations << " iterations)\n";
      check_bool( "the joint system solves", R.status == BASE_KINSOL::OK );
      check( "x(0) = 1",     y[0],  1., 1e-10 );
      check( "z(0) = 1",     y[1],  1., 1e-10 );
      check( "dx/dt(0) = -1", yp[0], -1., 1e-10 );
      // the same solve through the SPARSE (KLU) linear solver: same pattern the integrator's Jacobian
      // uses, and it must land on the same point
      M.options.INITIAL.KIN.SPARSE = true;
      std::vector<double> ys, yps;
      auto const Rs = M.solve_ic( ys, yps );
      check_bool( "SPARSE (KLU) IC solve converges", Rs.status == BASE_KINSOL::OK );
      if( Rs.status == BASE_KINSOL::OK ){
        check( "SPARSE agrees with DENSE on x(0)",     ys[0],  y[0],  1e-14 );
        check( "SPARSE agrees with DENSE on dx/dt(0)", yps[0], yp[0], 1e-14 );
      }
      M.options.INITIAL.KIN.SPARSE = false;
      if( M.solve( std::vector<double>{} ) == DAESLV_IDAS::NORMAL )
        check( "x(1) = (1+2t)^(-1/2)", M.val_state().back()[0], 1./std::sqrt(3.), 1e-8 );
      else { std::cout << "    solve failed: " << M.extract_error() << "\n"; ++n_fail; }
    } }

  std::cout << "\n==========  DAESLV_inputs: " << n_pass << " passed, " << n_fail << " failed -- "
            << ( n_fail? "FAILURES": "ALL PASS" ) << "  ==========\n\n";
  return n_fail? 1: 0;
}
