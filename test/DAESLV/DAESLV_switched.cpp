// DAESLV_switched -- piecewise-constant inputs: stage restarts with consistent reinitialisation.
// =====================================================================================================
//  D  y' = -y + u , 0 = z - u y with u = 1 then 2.  y is continuous across the switch, the ALGEBRAIC z
//     JUMPS, and numIC == 2 shows the restart happened rather than the integrator stepping through.
//  E  the same model through DENSE, DENSEDQ and SPARSE (KLU), which must agree.
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
  // ---------------------------------------------------------------- D
  { std::cout << "\n--- D. PIECEWISE-CONSTANT input: y' = -y + u, 0 = z - u y, u = 1 then 2\n";
    FFGraph G; DAESLV_IDAS M; M.set_dag( &G );
    FFVar t=G.add_var("t"), y=G.add_var("y(t)"), z=G.add_var("z(t)"), u=G.add_var("u(t)");
    FFPartial OpP; int const TI = FFDom::ALL - FFDom::LB;
    M.add_domain( t, FFDom( 0., 2., 2, FFDom::LGR, 4 ) );     // two elements -> the switch sits at t = 1
    M.add_state( y, {t} ); M.add_state( z, {t} );
    M.add_input( u, {t}, FFDom::LGL, 1 );                      // ONE node per element = piecewise constant
    M.set_evolution_domain( t ); M.update_ref( y, 0. ); M.update_ref( z, 0. );
    auto io=FFModel::EqnOptions(FFModel::EqnRole::INTERIOR,0), ii=FFModel::EqnOptions(FFModel::EqnRole::INITIAL,0);
    M.add_equation( OpP(y,t) + y - u, {t}, {TI},         io );
    M.add_equation( z - u*y,          {t}, {FFDom::ALL}, io );
    M.add_equation( y - 0.,           {t}, {FFDom::LB},  ii );
    M.FFModel::options.DISPLAY_LEVEL = 0;
    M.options.RTOL = M.options.ATOL = 1e-10;
    if( !M.setup() ){ std::cout << "    setup failed: " << M.extract_error() << "\n"; ++n_fail; }
    else{
      std::cout << "    input levels=" << M.n_input_levels() << " differential=" << M.n_differential()
                << " algebraic=" << M.n_algebraic() << "\n";
      if( !M.solve( { 1., 2. }, { 2. } ) ){ std::cout << "    solve failed: " << M.extract_error() << "\n"; ++n_fail; }
      else{
        double const y1 = 1. - std::exp(-1.);
        double const y2 = 2. + ( y1 - 2. ) * std::exp(-1.);
        check( "y(2) = 2 + (y(1)-2)e^-1",  M.val_state()[0], y2,      1e-8 );
        check( "z(2) = u y(2) with u = 2", M.val_state()[1], 2.*y2,   1e-8 );
        check( "z tracks the NEW level, not the old", M.val_state()[1] > 1.5*y2, true*1., 0.5 );
        std::cout << "    consistent initialisations = " << M.stats.numIC
                  << "  (1 at t0 + 1 per restarting stage boundary)\n";
        check( "restarted at the switch (numIC == 2)", (double)M.stats.numIC, 2., 0. );
        std::cout << "    steps=" << M.stats.numSteps << " res=" << M.stats.numRES
                  << " jac=" << M.stats.numJAC << " (analytic)\n";
      } } }

  // ---------------------------------------------------------------- E
  { std::cout << "\n--- E. same switched model through DENSE, DENSEDQ and SPARSE (KLU)\n";
    double ref0 = 0., ref1 = 0.; size_t nnz = 0;
    for( int mode = 0; mode < 3; ++mode ){
      FFGraph G; DAESLV_IDAS M; M.set_dag( &G );
      FFVar t=G.add_var("t"), y=G.add_var("y(t)"), z=G.add_var("z(t)"), u=G.add_var("u(t)");
      FFPartial OpP; int const TI = FFDom::ALL - FFDom::LB;
      M.add_domain( t, FFDom( 0., 2., 2, FFDom::LGR, 4 ) );
      M.add_state( y, {t} ); M.add_state( z, {t} );
      M.add_input( u, {t}, FFDom::LGL, 1 );
      M.set_evolution_domain( t ); M.update_ref( y, 0. ); M.update_ref( z, 0. );
      auto io=FFModel::EqnOptions(FFModel::EqnRole::INTERIOR,0), ii=FFModel::EqnOptions(FFModel::EqnRole::INITIAL,0);
      M.add_equation( OpP(y,t) + y - u, {t}, {TI},         io );
      M.add_equation( z - u*y,          {t}, {FFDom::ALL}, io );
      M.add_equation( y - 0.,           {t}, {FFDom::LB},  ii );
      M.FFModel::options.DISPLAY_LEVEL = 0;
      M.options.RTOL = M.options.ATOL = 1e-10;
      M.options.LINSOL = ( mode == 0 )? BASE_IDAS::Options::DENSE
                       : ( mode == 1 )? BASE_IDAS::Options::DENSEDQ
                                      : BASE_IDAS::Options::SPARSE;
      char const* nm = ( mode == 0 )? "DENSE   ": ( mode == 1 )? "DENSEDQ ": "SPARSE  ";
      if( !M.setup() || !M.solve( { 1., 2. }, { 2. } ) ){
        std::cout << "    " << nm << " failed: " << M.extract_error() << "\n"; ++n_fail; continue; }
      if( mode == 0 ){ ref0 = M.val_state()[0]; ref1 = M.val_state()[1]; nnz = M.nnz_jac();
        std::cout << "    Jacobian nonzeros (union of dF/dy and dF/dy') = " << nnz << " of "
                  << M.id().size()*M.id().size() << "\n"; }
      std::cout << "    " << nm << " y(2)=" << std::setprecision(12) << M.val_state()[0]
                << "  jac evals=" << M.stats.numJAC << "\n";
      check( std::string(nm) + "agrees with DENSE on y(2)", M.val_state()[0], ref0, 1e-7 );
      check( std::string(nm) + "agrees with DENSE on z(2)", M.val_state()[1], ref1, 1e-7 );
    } }


  std::cout << "\n==========  DAESLV_switched: " << n_pass << " passed, " << n_fail << " failed -- "
            << ( n_fail? "FAILURES": "ALL PASS" ) << "  ==========\n\n";
  return n_fail? 1: 0;
}
