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
    M.options.DISPLAY_LEVEL = 0;
    M.options.SOLVE.RTOL = M.options.SOLVE.ATOL = 1e-10;
    if( !M.setup() ){ std::cout << "    setup failed: " << M.extract_error() << "\n"; ++n_fail; }
    else{
      std::cout << "    input levels=" << M.n_input_levels() << " differential=" << M.n_differential()
                << " algebraic=" << M.n_algebraic() << "\n";
      // by NAME, as ODESLV_CVODES allows: assemble_values() checks the list and names any offender
      if( M.solve( { { u, { 1., 2. } } } ) != DAESLV_IDAS::NORMAL ){
        std::cout << "    solve failed: " << M.extract_error() << "\n"; ++n_fail; }
      else{
        double const y1 = 1. - std::exp(-1.);
        double const y2 = 2. + ( y1 - 2. ) * std::exp(-1.);
        check( "y(2) = 2 + (y(1)-2)e^-1",  M.val_state().back()[0], y2,      1e-8 );
        check( "z(2) = u y(2) with u = 2", M.val_state().back()[1], 2.*y2,   1e-8 );
        check( "z tracks the NEW level, not the old", M.val_state().back()[1] > 1.5*y2, true*1., 0.5 );
        std::cout << "    consistent initialisations = " << M.stats.numIC
                  << "  (1 at t0 + 1 per restarting stage boundary)\n";
        check( "restarted at the switch (numIC == 2)", (double)M.stats.numIC, 2., 0. );

        // What the restart COMPUTES, not merely that it happened.  Across the switch the differential y
        // is continuous while its DERIVATIVE jumps by exactly the input step: y' = -y + u, so
        // y'(1-) = -y(1) + 1 and y'(1+) = -y(1) + 2, differing by 1.  Nothing asserted this before, so a
        // reinitialisation that produced a plausible-but-wrong post-boundary state would have passed.
        { std::vector<double> xL{ y1, 1.*y1 }, xp, dxp;          // the left state at the switch
          BASE_KINSOL::Knobs kk;
          auto const R = M.restart_ic( 1, xL, xp, dxp, kk, { 1., 2. } );
          check( "restart_ic resolves the boundary", R.status == BASE_KINSOL::OK? 1.: 0., 1., 0. );
          if( R.status == BASE_KINSOL::OK ){
            check( "y continuous across the switch",        xp[0],  y1,          1e-12 );
            check( "z jumps to the NEW level, u+ y+",       xp[1],  2.*y1,       1e-12 );
            check( "dy/dt(1+) = -y(1) + 2",                 dxp[0], -y1 + 2.,    1e-12 );
            check( "the derivative JUMPS by the input step", dxp[0] - ( -y1 + 1. ), 1., 1e-12 );
          } }

        // and the INTEGRATOR driven through that same path must give the same trajectory as IDACalcIC
        { FFGraph G2; DAESLV_IDAS M2; M2.set_dag( &G2 );
          FFVar t2=G2.add_var("t"), y2v=G2.add_var("y(t)"), z2=G2.add_var("z(t)"), u2=G2.add_var("u(t)");
          FFPartial P2; int const TI2 = FFDom::ALL - FFDom::LB;
          M2.add_domain( t2, FFDom( 0., 2., 2, FFDom::LGR, 4 ) );
          M2.add_state( y2v, {t2} ); M2.add_state( z2, {t2} );
          M2.add_input( u2, {t2}, FFDom::LGL, 1 );
          M2.set_evolution_domain( t2 ); M2.update_ref( y2v, 1. ); M2.update_ref( z2, 1. );
          auto o2=FFModel::EqnOptions(FFModel::EqnRole::INTERIOR,0), i2=FFModel::EqnOptions(FFModel::EqnRole::INITIAL,0);
          M2.add_equation( P2(y2v,t2) + y2v - u2, {t2}, {TI2},        o2 );
          M2.add_equation( z2 - u2*y2v,           {t2}, {FFDom::ALL}, o2 );
          M2.add_equation( y2v - 0.,              {t2}, {FFDom::LB},  i2 );
          M2.options.DISPLAY_LEVEL = 0;
          M2.options.SOLVE.RTOL = M2.options.SOLVE.ATOL = 1e-10;
          M2.options.INITIAL.SOLVER = BASE_IDAS::Options::t_Initial::KINSOL;
          if( M2.setup() && M2.solve( { 1., 2. } ) == DAESLV_IDAS::NORMAL ){
            check( "KINSOL restart gives the same y(2) as IDACalcIC", M2.val_state().back()[0], y2, 1e-8 );
            check( "KINSOL restart gives the same z(2)",              M2.val_state().back()[1], 2.*y2, 1e-8 );
            check( "and still restarted once (numIC == 2)", (double)M2.stats.numIC, 2., 0. );
          }
          else { std::cout << "    KINSOL-restart solve failed: " << M2.extract_error() << "\n"; ++n_fail; } }
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
      M.options.DISPLAY_LEVEL = 0;
      M.options.SOLVE.RTOL = M.options.SOLVE.ATOL = 1e-10;
      M.options.SOLVE.LINSOL = ( mode == 0 )? BASE_IDAS::Options::DENSE
                       : ( mode == 1 )? BASE_IDAS::Options::DENSEDQ
                                      : BASE_IDAS::Options::SPARSE;
      char const* nm = ( mode == 0 )? "DENSE   ": ( mode == 1 )? "DENSEDQ ": "SPARSE  ";
      if( !M.setup() || M.solve( { 1., 2. } ) != DAESLV_IDAS::NORMAL ){
        std::cout << "    " << nm << " failed: " << M.extract_error() << "\n"; ++n_fail; continue; }
      if( mode == 0 ){ ref0 = M.val_state().back()[0]; ref1 = M.val_state().back()[1]; nnz = M.nnz_jac();
        std::cout << "    Jacobian nonzeros (union of dF/dy and dF/dy') = " << nnz << " of "
                  << M.id().size()*M.id().size() << "\n"; }
      std::cout << "    " << nm << " y(2)=" << std::setprecision(12) << M.val_state().back()[0]
                << "  jac evals=" << M.stats.numJAC << "\n";
      check( std::string(nm) + "agrees with DENSE on y(2)", M.val_state().back()[0], ref0, 1e-7 );
      check( std::string(nm) + "agrees with DENSE on z(2)", M.val_state().back()[1], ref1, 1e-7 );
    } }


  std::cout << "\n==========  DAESLV_switched: " << n_pass << " passed, " << n_fail << " failed -- "
            << ( n_fail? "FAILURES": "ALL PASS" ) << "  ==========\n\n";
  return n_fail? 1: 0;
}
