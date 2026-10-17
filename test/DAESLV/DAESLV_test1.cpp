// DAESLV_test1 -- DAE integration against closed forms, and the ODE cross-check against ODESLV.
// =====================================================================================================
//  A  semi-explicit index-1 DAE   y' = -y z , 0 = z - y^2 , y(0)=1   vs the closed form
//  B  harmonic oscillator          x0' = x1 , x1' = -x0              DAESLV vs ODESLV vs closed form
//  C  COUPLED mass matrix          ODESLV cannot take it; DAESLV integrates it implicitly and must agree
//     with ODESLV on the EXPLICIT form of the same system.
#include <cmath>
#include <iomanip>
#include <iostream>
#include "daeslv_idas.hpp"
#include "odeslvs_cvodes.hpp"

using namespace mc;

//! Exposes the two callbacks so the ANALYTIC Jacobian can be checked against finite differences of the
//! residual itself -- agreement between solver backends is weak evidence, since a wrong Jacobian usually
//! still converges, just more slowly.
struct JacProbe : public DAESLV_IDAS {
  int res( double t, N_Vector y, N_Vector yp, N_Vector r ){ return IDRES__( t, y, yp, r, nullptr ); }
  int jac( double t, double cj, N_Vector y, N_Vector yp, SUNMatrix J )
    { return IDJAC__( t, cj, y, yp, nullptr, J, nullptr, nullptr, nullptr, nullptr ); }
};

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
static double const T_END = 1.0;

int main()
{
  // ---------------------------------------------------------------- A
  { std::cout << "\n--- A. semi-explicit index-1 DAE: y' = -y z, 0 = z - y^2, y(0) = 1\n";
    FFGraph G; DAESLV_IDAS M; M.set_dag( &G );
    FFVar t = G.add_var("t"), y = G.add_var("y(t)"), z = G.add_var("z(t)");
    FFPartial OpP;  int const TI = FFDom::ALL - FFDom::LB;
    M.add_domain( t, FFDom( 0., T_END, 2, FFDom::LGR, 4 ) );
    M.add_state( y, {t} ); M.add_state( z, {t} );
    M.set_evolution_domain( t ); M.update_ref( y, 1. ); M.update_ref( z, 1. );
    auto io = FFModel::EqnOptions( FFModel::EqnRole::INTERIOR, 0 );
    auto ii = FFModel::EqnOptions( FFModel::EqnRole::INITIAL,  0 );
    M.add_equation( OpP(y,t) + y*z, {t}, {TI},          io );
    M.add_equation( z - y*y,        {t}, {FFDom::ALL},  io );
    M.add_equation( y - 1.,         {t}, {FFDom::LB},   ii );
    M.options.DISPLAY_LEVEL = 0;
    M.options.SOLVE.RTOL = M.options.SOLVE.ATOL = 1e-10;
    if( !M.setup() ){ std::cout << "    setup failed: " << M.extract_error() << "\n"; ++n_fail; }
    else{
      std::cout << "    rows=" << M.n_residual() << " differential=" << M.n_differential()
                << " algebraic=" << M.n_algebraic() << "\n";
      if( M.solve( {} ) != DAESLV_IDAS::NORMAL ){ std::cout << "    solve failed: " << M.extract_error() << "\n"; ++n_fail; }
      else{
        check( "IDACalcIC: z(0) = y(0)^2 = 1", M.val_state_init()[1], 1., 1e-9 );
        check( "IDACalcIC: y'(0) = -1",        M.val_deriv_init()[0], -1., 1e-9 );
        check( "y(1) = (1+2t)^(-1/2)", M.val_state().back()[0], 1./std::sqrt(3.), 1e-8 );
        check( "z(1) = (1+2t)^(-1)",   M.val_state().back()[1], 1./3.,            1e-8 );
        // solve() must not depend on what the previous call left behind: the INITIAL rows can carry
        // states, and evaluating them against stale values once gave a different answer each run.
        std::vector<double> const first = M.val_state().back();
        M.solve( {} );
        check( "solve() is idempotent on y", M.val_state().back()[0], first[0], 0. );
        check( "solve() is idempotent on z", M.val_state().back()[1], first[1], 0. );
        // the same model through the KINSOL corrector: it resolves the INITIAL rows TOGETHER WITH the
        // interior rows at t0, so it must land on the same consistent point IDACalcIC does.
        M.options.INITIAL.SOLVER = BASE_IDAS::Options::t_Initial::KINSOL;
        if( M.solve( {} ) != DAESLV_IDAS::NORMAL ){
          std::cout << "    KINSOL path failed: " << M.extract_error() << "\n"; ++n_fail; }
        else{
          check( "KINSOL corrector agrees with IDACalcIC on y", M.val_state().back()[0], first[0], 1e-9 );
          check( "KINSOL corrector agrees with IDACalcIC on z", M.val_state().back()[1], first[1], 1e-9 );
        }
        M.options.INITIAL.SOLVER = BASE_IDAS::Options::t_Initial::IDACALCIC;
        std::cout << "    steps=" << M.stats.numSteps << " res=" << M.stats.numRES
                  << " jac=" << M.stats.numJAC << "\n";
      } } }

  // ---------------------------------------------------------------- B
  { std::cout << "\n--- B. harmonic oscillator x0'=x1, x1'=-x0  (DAESLV vs ODESLV vs closed form)\n";
    double dae0 = 0., dae1 = 0., ode0 = 0., ode1 = 0.;
    { FFGraph G; DAESLV_IDAS M; M.set_dag( &G );
      FFVar t=G.add_var("t"), x0=G.add_var("x0(t)"), x1=G.add_var("x1(t)");
      FFPartial OpP; int const TI = FFDom::ALL - FFDom::LB;
      M.add_domain( t, FFDom( 0., T_END, 2, FFDom::LGR, 4 ) );
      M.add_state( x0, {t} ); M.add_state( x1, {t} );
      M.set_evolution_domain( t ); M.update_ref( x0, 0. ); M.update_ref( x1, 1. );
      auto io=FFModel::EqnOptions(FFModel::EqnRole::INTERIOR,0), ii=FFModel::EqnOptions(FFModel::EqnRole::INITIAL,0);
      M.add_equation( OpP(x0,t) - x1, {t}, {TI}, io );
      M.add_equation( OpP(x1,t) + x0, {t}, {TI}, io );
      M.add_equation( x0 - 0., {t}, {FFDom::LB}, ii );
      M.add_equation( x1 - 1., {t}, {FFDom::LB}, ii );
      M.options.DISPLAY_LEVEL = 0; M.options.SOLVE.RTOL=M.options.SOLVE.ATOL=1e-10;
      if( M.setup() && M.solve( {} ) == DAESLV_IDAS::NORMAL ){ dae0=M.val_state().back()[0]; dae1=M.val_state().back()[1]; }
      else { std::cout << "    DAESLV failed: " << M.extract_error() << "\n"; ++n_fail; } }
    { FFGraph G; ODESLVS_CVODES A; A.set_dag( &G );
      FFVar t=G.add_var("t"), x0=G.add_var("x0(t)"), x1=G.add_var("x1(t)");
      FFPartial OpP; FFEval OpE; int const TI = FFDom::ALL - FFDom::LB;
      A.add_domain( t, FFDom( 0., T_END, 2, FFDom::LGR, 4 ) );
      A.add_state( x0, {t} ); A.add_state( x1, {t} );
      A.set_evolution_domain( t ); A.update_ref( x0, 0. ); A.update_ref( x1, 1. );
      auto io=FFModel::EqnOptions(FFModel::EqnRole::INTERIOR,0), ii=FFModel::EqnOptions(FFModel::EqnRole::INITIAL,0);
      A.add_equation( OpP(x0,t) - x1, {t}, {TI}, io );
      A.add_equation( OpP(x1,t) + x0, {t}, {TI}, io );
      A.add_equation( x0 - 0., {t}, {FFDom::LB}, ii );
      A.add_equation( x1 - 1., {t}, {FFDom::LB}, ii );
      A.add_output( OpE( x0, t, T_END ) ); A.add_output( OpE( x1, t, T_END ) );
      A.options.SOLVE.RTOL=A.options.SOLVE.ATOL=1e-10; A.options.SOLVE.VERBOSE = false;
      if( A.setup() ){ A.solve( std::vector<double>{} ); ode0=A.val_function()[0]; ode1=A.val_function()[1]; }
      else { std::cout << "    ODESLV setup failed\n"; ++n_fail; } }
    check( "DAESLV x0(1) = sin(1)", dae0, std::sin(T_END), 1e-8 );
    check( "DAESLV x1(1) = cos(1)", dae1, std::cos(T_END), 1e-8 );
    check( "ODESLV x0(1) = sin(1)", ode0, std::sin(T_END), 1e-8 );
    check( "DAESLV == ODESLV on x0", dae0, ode0, 1e-8 );
    check( "DAESLV == ODESLV on x1", dae1, ode1, 1e-8 ); }

  // ---------------------------------------------------------------- C
  { std::cout << "\n--- C. coupled mass matrix (DAESLV implicit) vs ODESLV on the inverted form\n";
    double d0=0., d1=0., o0=0., o1=0.;
    { FFGraph G; DAESLV_IDAS M; M.set_dag( &G );
      FFVar t=G.add_var("t"), x0=G.add_var("x0(t)"), x1=G.add_var("x1(t)");
      FFPartial OpP; int const TI = FFDom::ALL - FFDom::LB;
      M.add_domain( t, FFDom( 0., T_END, 2, FFDom::LGR, 4 ) );
      M.add_state( x0, {t} ); M.add_state( x1, {t} );
      M.set_evolution_domain( t ); M.update_ref( x0, 1. ); M.update_ref( x1, 2. );
      auto io=FFModel::EqnOptions(FFModel::EqnRole::INTERIOR,0), ii=FFModel::EqnOptions(FFModel::EqnRole::INITIAL,0);
      M.add_equation( 2.*OpP(x0,t) + 1.*OpP(x1,t) - ( x0 + 2.*x1 ), {t}, {TI}, io );
      M.add_equation( 1.*OpP(x0,t) + 3.*OpP(x1,t) - ( 3.*x0 - x1 ), {t}, {TI}, io );
      M.add_equation( x0 - 1., {t}, {FFDom::LB}, ii );
      M.add_equation( x1 - 2., {t}, {FFDom::LB}, ii );
      M.options.DISPLAY_LEVEL = 0; M.options.SOLVE.RTOL=M.options.SOLVE.ATOL=1e-10;
      if( M.setup() && M.solve( {} ) == DAESLV_IDAS::NORMAL ){ d0=M.val_state().back()[0]; d1=M.val_state().back()[1]; }
      else { std::cout << "    DAESLV failed: " << M.extract_error() << "\n"; ++n_fail; } }
    { FFGraph G; ODESLVS_CVODES A; A.set_dag( &G );
      FFVar t=G.add_var("t"), x0=G.add_var("x0(t)"), x1=G.add_var("x1(t)");
      FFPartial OpP; FFEval OpE; int const TI = FFDom::ALL - FFDom::LB;
      A.add_domain( t, FFDom( 0., T_END, 2, FFDom::LGR, 4 ) );
      A.add_state( x0, {t} ); A.add_state( x1, {t} );
      A.set_evolution_domain( t ); A.update_ref( x0, 1. ); A.update_ref( x1, 2. );
      auto io=FFModel::EqnOptions(FFModel::EqnRole::INTERIOR,0), ii=FFModel::EqnOptions(FFModel::EqnRole::INITIAL,0);
      A.add_equation( OpP(x0,t) - 1.4*x1,            {t}, {TI}, io );   // M^-1 applied by hand
      A.add_equation( OpP(x1,t) - ( x0 - 0.8*x1 ),   {t}, {TI}, io );
      A.add_equation( x0 - 1., {t}, {FFDom::LB}, ii );
      A.add_equation( x1 - 2., {t}, {FFDom::LB}, ii );
      A.add_output( OpE( x0, t, T_END ) ); A.add_output( OpE( x1, t, T_END ) );
      A.options.SOLVE.RTOL=A.options.SOLVE.ATOL=1e-10; A.options.SOLVE.VERBOSE = false;
      if( A.setup() ){ A.solve( std::vector<double>{} ); o0=A.val_function()[0]; o1=A.val_function()[1]; }
      else { std::cout << "    ODESLV setup failed\n"; ++n_fail; } }
    check( "implicit == explicit on x0", d0, o0, 1e-7 );
    check( "implicit == explicit on x1", d1, o1, 1e-7 ); }


  std::cout << "\n==========  DAESLV_test1: " << n_pass << " passed, " << n_fail << " failed -- "
            << ( n_fail? "FAILURES": "ALL PASS" ) << "  ==========\n\n";
  return n_fail? 1: 0;
}
