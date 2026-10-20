// DAESLV_icsens -- the SENSITIVITY initial conditions from the shared system.
// =====================================================================================================
// ic_sens() solves ( dF/dz ) . ( dz/dp ) = - dF/dp against the Jacobian the state solve just used.
// The model below has INITIAL rows implicit in each other and carrying the parameter, which is exactly
// the case the divided form cannot represent:
//
//   x0' = -x0 ,  x1' = -2 x1
//   x0(0) + x1(0) - 3 p = 0      =>  x0(0) = 2p ,  x1(0) = p
//   x0(0) - x1(0) -   p = 0          d x0(0)/dp = 2 ,  d x1(0)/dp = 1
//                                    d x0'(0)/dp = -2 ,  d x1'(0)/dp = -2
//
// Everything is linear in p, so the expected values are exact and the comparison has no tolerance
// argument.  The EXPLICIT declaration of the same data is the control: both routes must agree.
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

static void run( char const* label, bool implicit )
{
  std::cout << "\n--- " << label << "\n";
  FFGraph G; DAESLV_IDAS M; M.set_dag( &G );
  FFPartial OpP;
  int const TI = FFDom::ALL - FFDom::LB;
  auto io = FFModel::EqnOptions( FFModel::EqnRole::INTERIOR, 0 );
  auto ii = FFModel::EqnOptions( FFModel::EqnRole::INITIAL,  0 );
  FFVar t = G.add_var("t"), x0 = G.add_var("x0(t)"), x1 = G.add_var("x1(t)"), p = G.add_var("p");
  M.add_domain( t, FFDom( 0., 1., 1, FFDom::LGR, 4 ) );
  M.add_state( x0, {t} ); M.add_state( x1, {t} ); M.add_input( p );
  M.set_evolution_domain( t );
  M.update_ref( x0, 2. ); M.update_ref( x1, 1. );
  M.add_equation( OpP( x0, t ) + x0,    {t}, {TI}, io );
  M.add_equation( OpP( x1, t ) + 2.*x1, {t}, {TI}, io );
  if( implicit ){
    M.add_equation( x0 + x1 - 3.*p, {t}, {FFDom::LB}, ii );
    M.add_equation( x0 - x1 -    p, {t}, {FFDom::LB}, ii );
  }
  else{
    M.add_equation( x0 - 2.*p, {t}, {FFDom::LB}, ii );
    M.add_equation( x1 -    p, {t}, {FFDom::LB}, ii );
  }
  M.options.SOLVE.RTOL = M.options.SOLVE.ATOL = 1e-11;
  M.options.INITIAL.BYPASS_EXPLICIT = false;      // force the SOLVE on both, so they are comparable
  M.options.DISPLAY_LEVEL = 0;
  if( !M.setup() ){ std::cout << "    setup failed: " << M.extract_error() << "\n"; ++n_fail; return; }

  double const P = 1.7;
  std::vector<double> xv, dxv;
  BASE_KINSOL::Knobs k;  k.FTOL = 1e-12;
  auto const R = M.BASE_DYNAMIC::solve_ic( xv, dxv, k, std::vector<double>{ P } );
  if( R.status == BASE_KINSOL::FAILED ){ std::cout << "    solve_ic failed: " << R.message << "\n"; ++n_fail; return; }
  check( "x0(0) = 2p", xv.at(0), 2.*P, 1e-9 );
  check( "x1(0) =  p", xv.at(1),    P, 1e-9 );

  std::vector< std::vector<double> > dxdp, ddxdp;
  if( !M.ic_sens( dxdp, ddxdp ) ){ std::cout << "    ic_sens failed\n"; ++n_fail; return; }
  std::cout << "    ic_sens: " << dxdp.size() << " parameter(s) x " << dxdp[0].size() << " state(s)\n";
  check( "d x0(0) / dp = 2",  dxdp[0].at(0),  2., 1e-9 );
  check( "d x1(0) / dp = 1",  dxdp[0].at(1),  1., 1e-9 );
  // x0' = -x0 so d x0'(0)/dp = -2 ; x1' = -2 x1 so d x1'(0)/dp = -2
  check( "d x0'(0) / dp = -2", ddxdp[0].at(0), -2., 1e-9 );
  check( "d x1'(0) / dp = -2", ddxdp[0].at(1), -2., 1e-9 );
}

int main()
{
  std::cout << "\n=== ic_sens: ( dF/dz ) . ( dz/dp ) = - dF/dp\n";
  run( "CONTROL: the same data declared EXPLICITLY", false );
  run( "the data declared IMPLICITLY (coupled rows)", true );
  // ------------------------------------------------------------------ HIGH INDEX
  // Do the initial sensitivities account for the HIDDEN constraint sensitivities?  They must, and the
  // mechanism is that REDUCE.HIDDEN_IC materialises every hidden level as an INITIAL row -- so the
  // levels are ROWS of the system ic_sens() differentiates, not an afterthought bolted onto it.
  //
  //   x' = z + u ,  0 = x - p        INDEX 2.  Differentiating pins z: z = -u.
  //                                  The HIDDEN level is x - p = 0, and it CARRIES the parameter.
  //   so x(0) = p, dx/dp = 1 EXACTLY, and dz/dp = 0.
  //
  // The discriminating part is dx/dp.  x is pinned by the hidden level ALONE -- no other row mentions
  // it without a derivative -- so if the hidden level were absent from the system, x would not be
  // determined and dx/dp could not come out 1.
  { std::cout << "\n--- INDEX-2: the hidden level x - p = 0 carries the parameter\n";
    FFGraph G; DAESLV_IDAS M; M.set_dag( &G );
    FFPartial OpP;
    int const TI2 = FFDom::ALL - FFDom::LB;
    auto io2 = FFModel::EqnOptions( FFModel::EqnRole::INTERIOR, 0 );
    FFVar t = G.add_var("t"), x = G.add_var("x(t)"), z = G.add_var("z(t)"),
          u = G.add_var("u"), pp = G.add_var("p");
    M.add_domain( t, FFDom( 0., 1., 1, FFDom::LGR, 4 ) );
    M.add_state( x, {t} ); M.add_state( z, {t} );
    M.add_input( u ); M.add_input( pp );
    M.set_evolution_domain( t );
    M.update_ref( x, 1. ); M.update_ref( z, 0. );
    M.add_equation( OpP( x, t ) - z - u, {t}, {TI2}, io2 );
    M.add_equation( x - pp,              {t}, {FFDom::ALL}, io2 );     // INDEX 2
    M.options.SOLVE.RTOL = M.options.SOLVE.ATOL = 1e-11;
    M.options.DISPLAY_LEVEL = 0;
    if( !M.setup() ){ std::cout << "    setup failed: " << M.extract_error() << "\n"; ++n_fail; }
    else{
      std::cout << "    max_index=" << M.reduction_plan().max_index
                << "  hidden levels kept=" << M.hidden_levels().size()
                << "  ic_rows=" << M.ic_rows().size() << " ic_unknowns=" << M.ic_unknowns().size() << "\n";
      check( "a hidden level was materialised", (double)M.hidden_levels().size(), 1., 0. );
      double const U = 0.3, PV = 2.4;
      std::vector<double> xv, dxv;
      BASE_KINSOL::Knobs k;  k.FTOL = 1e-12;
      auto const R = M.BASE_DYNAMIC::solve_ic( xv, dxv, k, std::vector<double>{ U, PV } );
      if( R.status == BASE_KINSOL::FAILED ){ std::cout << "    solve_ic failed: " << R.message << "\n"; ++n_fail; }
      else{
        check( "x(0) = p (from the hidden level)", xv.at(0), PV, 1e-9 );
        check( "z(0) = -u", xv.at(1), -U, 1e-9 );
        std::vector< std::vector<double> > dxdp, ddxdp;
        if( !M.ic_sens( dxdp, ddxdp ) ){ std::cout << "    ic_sens failed\n"; ++n_fail; }
        else{
          size_t const ju = M.parameter_index( u )[0], jp = M.parameter_index( pp )[0];
          check( "d x(0) / dp = 1  (HIDDEN level)", dxdp.at(jp).at(0), 1., 1e-9 );
          check( "d z(0) / dp = 0",                 dxdp.at(jp).at(1), 0., 1e-9 );
          check( "d x(0) / du = 0",                 dxdp.at(ju).at(0), 0., 1e-9 );
          check( "d z(0) / du = -1",                dxdp.at(ju).at(1), -1., 1e-9 );
        }
      }
    }
  }

  std::cout << "\n==========  DAESLV_icsens: " << n_pass << " passed, " << n_fail << " failed -- "
            << ( n_fail? "FAILURES": "ALL PASS" ) << "  ==========\n\n";
  return n_fail? 1: 0;
}
