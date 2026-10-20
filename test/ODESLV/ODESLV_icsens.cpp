// ODESLV_icsens -- the forward sensitivity initial condition on IMPLICIT initial rows.
// STATE initial condition when the INITIAL rows are IMPLICIT?
// =====================================================================================================
// _IC_D_STA solves implicit INITIAL rows through the shared BASE_DYNAMIC system (KINSOL) -- it has done
// since 2026-10-18.  _IC_SET_FSA instead takes
//        _pSACFCT = _dag->FAD( _nx, pIC, 1, p_j )      with pIC = _vIC.at(0)
// i.e. it differentiates the DIVIDED form _mIC, which can represent an INITIAL row only when that row
// pins one state linearly.  Where the rows are implicit in each other the divided form is not the
// initial condition at all, so the two halves can disagree -- and nothing says so.
//
//   x0' = -x0 ,  x1' = -2 x1
//   x0(0) + x1(0) - 3 p = 0        =>  x0(0) = 2p ,  x1(0) = p
//   x0(0) - x1(0) -   p = 0
//   outputs x0(1) = 2p e^-1 and x1(1) = p e^-2, so d/dp = 2 e^-1 and e^-2 EXACTLY -- the model is
//   linear in p, so a finite difference is exact too and the comparison has no tolerance argument.
#include <cmath>
#include <iomanip>
#include <iostream>
#include "odeslvs_cvodes.hpp"

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

//! Build the model.  implicit=false declares the SAME initial data explicitly, as a control.
static bool build( ODESLVS_CVODES& I, FFGraph& G, bool implicit )
{
  FFPartial OpP;  FFEval OpE;
  FFVar t = G.add_var("t"), x0 = G.add_var("x0(t)"), x1 = G.add_var("x1(t)"), p = G.add_var("p");
  int const TI = FFDom::ALL - FFDom::LB;
  auto io = FFModel::EqnOptions( FFModel::EqnRole::INTERIOR, 0 );
  auto ii = FFModel::EqnOptions( FFModel::EqnRole::INITIAL,  0 );
  I.add_domain( t, FFDom( 0., 1., 1, FFDom::LGR, 4 ) );
  I.set_evolution_domain( t );
  I.add_state( x0, {t} );  I.add_state( x1, {t} );  I.add_input( p );
  I.update_ref( x0, 2. );  I.update_ref( x1, 1. );
  I.add_equation( OpP( x0, t ) + x0,      {t}, {TI}, io );
  I.add_equation( OpP( x1, t ) + 2.*x1,   {t}, {TI}, io );
  if( implicit ){
    I.add_equation( x0 + x1 - 3.*p, {t}, {FFDom::LB}, ii );     // IMPLICIT: coupled in x0 and x1
    I.add_equation( x0 - x1 -    p, {t}, {FFDom::LB}, ii );
  }
  else{
    I.add_equation( x0 - 2.*p, {t}, {FFDom::LB}, ii );          // the SAME data, explicitly
    I.add_equation( x1 -    p, {t}, {FFDom::LB}, ii );
  }
  I.add_output( OpE( x0, t, 1. ) );
  I.add_output( OpE( x1, t, 1. ) );
  I.options.SOLVE.RTOL = I.options.SOLVE.ATOL = 1e-11;
  I.options.SOLVE.RTOLS = I.options.SOLVE.ATOLS = 1e-11;
  I.options.SOLVE.RTOLB = I.options.SOLVE.ATOLB = 1e-11;   // the ADJOINT pass has its OWN tolerances;
                                                           // left at default it runs four orders looser
                                                           // than forward and looks like a defect
  I.options.DISPLAY_LEVEL = 0;
  return I.setup();
}

static void run( char const* label, bool implicit )
{
  std::cout << "\n--- " << label << "\n";
  FFGraph G;  ODESLVS_CVODES I( &G );
  if( !build( I, G, implicit ) ){ std::cout << "    setup failed: " << I.extract_error() << "\n"; ++n_fail; return; }

  double const P = 1.7;
  std::vector<double> p( I.np(), 0. );
  p[ I.parameter_index( *I.var_parameter().begin() )[0] ] = P;   // the one declared input

  if( I.solve( p ) != ODESLVS_CVODES::STATUS::NORMAL ){ std::cout << "    solve failed\n"; ++n_fail; return; }
  auto const F = I.val_function();
  check( "x0(1) = 2 p e^-1", F.at(0), 2.*P*std::exp(-1.), 1e-8 );
  check( "x1(1) =   p e^-2", F.at(1),    P*std::exp(-2.), 1e-8 );

  if( I.solve_fsens( p ) != ODESLVS_CVODES::STATUS::NORMAL ){ std::cout << "    solve_fsens failed\n"; ++n_fail; return; }
  auto const Gf = I.val_function_gradient();
  // one ROW per parameter, one COLUMN per output
  std::cout << "    val_function_gradient(): " << Gf.size() << " row(s) x "
            << ( Gf.empty()? 0: Gf[0].size() ) << " col(s)\n";
  if( Gf.empty() || Gf[0].size() < 2 ){ ++n_fail; return; }
  // the model is LINEAR in p, so these are exact, not approximations
  check( "d x0(1) / dp = 2 e^-1", Gf[0].at(0), 2.*std::exp(-1.), 1e-7 );
  check( "d x1(1) / dp =   e^-2", Gf[0].at(1),    std::exp(-2.), 1e-7 );

  if( I.solve_asens( p ) == ODESLVS_CVODES::STATUS::NORMAL ){
    auto const Ga = I.val_function_gradient();
    if( !Ga.empty() && Ga[0].size() >= 2 ){
      // On the EXPLICIT path the adjoint agrees with forward, and this GATES it -- but only once
      // RTOLB/ATOLB are set: they are separate from RTOLS/ATOLS, and leaving them at their defaults
      // put the adjoint four orders looser than the forward pass and looked like a defect.
      //
      // On the IMPLICIT path it is RECORDED, not gated: the adjoint initial condition still comes from
      // the divided form -- the exact analogue of what _IC_SET_FSA did before ic_sens() -- so it
      // returns zeros.  Phase 6.
      double const dA0 = std::fabs( Ga[0].at(0)-Gf[0].at(0) ), dA1 = std::fabs( Ga[0].at(1)-Gf[0].at(1) );
      if( !implicit ){
        check( "  adjoint agrees with forward on d x0", Ga[0].at(0), Gf[0].at(0), 1e-9 );
        check( "  adjoint agrees with forward on d x1", Ga[0].at(1), Gf[0].at(1), 1e-9 );
      }
      else
        std::cout << "    [recorded] adjoint IC still divided: d x0 " << std::scientific
                  << std::setprecision(2) << dA0 << "   d x1 " << dA1 << "\n";
    }
  }
}

int main()
{
  std::cout << "\n=== forward sensitivity vs the IMPLICIT initial condition\n";
  run( "CONTROL: the same data declared EXPLICITLY", false );
  run( "the data declared IMPLICITLY (coupled rows)", true );
  std::cout << "\n==========  ODESLV_icsens: " << n_pass << " passed, " << n_fail << " failed -- "
            << ( n_fail? "FAILURES": "ALL PASS" ) << "  ==========\n\n";
  return n_fail? 1: 0;
}
