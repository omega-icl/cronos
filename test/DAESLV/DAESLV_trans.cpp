// DAESLV_trans -- state transitions (add_transition) applied at a stop.
// =====================================================================================================
// A transition replaces the CONTINUITY of the states it determines with its own rows, in the same stop
// system the restart already solves.  That is the whole mechanism: a jump is a row, not a mechanism of
// its own, which is why an IMPLICIT map costs nothing here where _mIC (ODESLV's divided path) refuses it.
//
//   A. EXPLICIT jump on a pure ODE:  x' = -x on [0,2], x(0)=1, with x(1^+) = 2 x(1^-).
//      Exact: x(t) = e^-t on [0,1), then 2e^-1 e^-(t-1) on (1,2], so x(2) = 2 e^-2.
//   B. IMPLICIT jump on an index-1 DAE: the post-jump state is given by a row that does NOT solve for
//      it explicitly -- x(1^+)^2 = x(1^-), i.e. the jump is the positive root.  _mIC cannot express it.
//   C. REFUSAL: a tau that is not an element boundary must be refused at setup, not ignored.
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
{ ok? ++n_pass: ++n_fail;
  std::cout << "    " << std::left << std::setw( 44 ) << w << std::right << ( ok? "  PASS\n": "  FAIL\n" ); }

int main()
{
  FFPartial OpP;
  int const TI = FFDom::ALL - FFDom::LB;
  auto io = FFModel::EqnOptions( FFModel::EqnRole::INTERIOR, 0 );
  auto ii = FFModel::EqnOptions( FFModel::EqnRole::INITIAL,  0 );

  // ----------------------------------------------------------------- A: explicit jump
  { std::cout << "\n--- A. x' = -x on [0,2], x(0)=1, EXPLICIT jump x(1+) = 2 x(1-)\n";
    FFGraph G; DAESLV_IDAS M; M.set_dag( &G );
    FFVar t = G.add_var("t"), x = G.add_var("x(t)");
    M.add_domain( t, FFDom( 0., 2., 2, FFDom::LGR, 4 ) );      // the boundary at t = 1 is the stop
    M.add_state( x, {t} ); M.set_evolution_domain( t ); M.update_ref( x, 1. );
    M.add_equation( OpP( x, t ) + x, {t}, {TI}, io );
    M.add_equation( x - 1.,          {t}, {FFDom::LB}, ii );
    M.add_transition( 2.*x, x, t, 1. );                        // left(tau^-) = right(tau^+)
    M.options.SOLVE.RTOL = M.options.SOLVE.ATOL = 1e-10;
    M.options.DISPLAY_LEVEL = 0;
    bool const ok = M.setup();
    std::cout << "    setup=" << ok << ( ok? std::string(): "  :: " + M.extract_error() ) << "\n";
    check_bool( "setup ACCEPTS a transition at a boundary", ok );
    if( ok ){
      M.solve( std::vector<double>{} );
      auto const& X = M.val_state();
      std::cout << "    numIC=" << M.stats.numIC << "   x at the stage times:";
      for( auto const& r : X ) std::cout << " " << std::scientific << std::setprecision(6) << r[0];
      std::cout << "\n";
      check( "x(1-) = e^-1 (the LEFT limit is recorded)", X.at(1)[0], std::exp(-1.), 1e-8 );
      check( "x(2)  = 2 e^-2 (the jump was applied)",     X.back()[0], 2.*std::exp(-2.), 1e-7 );
    }
  }

  // ----------------------------------------------------------------- B: implicit jump
  { std::cout << "\n--- B. IMPLICIT jump: x(1+)^2 = x(1-), which _mIC cannot express\n";
    FFGraph G; DAESLV_IDAS M; M.set_dag( &G );
    FFVar t = G.add_var("t"), x = G.add_var("x(t)");
    M.add_domain( t, FFDom( 0., 2., 2, FFDom::LGR, 4 ) );
    M.add_state( x, {t} ); M.set_evolution_domain( t ); M.update_ref( x, 1. );
    M.add_equation( OpP( x, t ) + x, {t}, {TI}, io );
    M.add_equation( x - 1.,          {t}, {FFDom::LB}, ii );
    M.add_transition( x, x*x, t, 1. );                         // left(tau^-) = right(tau^+)^2
    M.options.SOLVE.RTOL = M.options.SOLVE.ATOL = 1e-10;
    M.options.DISPLAY_LEVEL = 0;
    bool const ok = M.setup();
    std::cout << "    setup=" << ok << ( ok? std::string(): "  :: " + M.extract_error() ) << "\n";
    check_bool( "setup ACCEPTS an IMPLICIT jump", ok );
    if( ok ){
      M.solve( std::vector<double>{} );
      auto const& X = M.val_state();
      double const xm = std::exp(-1.), xp = std::sqrt( xm );   // positive root
      std::cout << "    numIC=" << M.stats.numIC << "   x(1-)=" << std::scientific << X.at(1)[0]
                << "   expected x(1+)=" << xp << "\n";
      check( "x(2) = sqrt(e^-1) e^-1", X.back()[0], xp*std::exp(-1.), 1e-7 );
    }
  }

  // ----------------------------------------------------------------- C: refusal off a boundary
  { std::cout << "\n--- C. a tau that is NOT an element boundary must be REFUSED\n";
    FFGraph G; DAESLV_IDAS M; M.set_dag( &G );
    FFVar t = G.add_var("t"), x = G.add_var("x(t)");
    M.add_domain( t, FFDom( 0., 2., 2, FFDom::LGR, 4 ) );      // boundaries at 0, 1, 2 only
    M.add_state( x, {t} ); M.set_evolution_domain( t ); M.update_ref( x, 1. );
    M.add_equation( OpP( x, t ) + x, {t}, {TI}, io );
    M.add_equation( x - 1.,          {t}, {FFDom::LB}, ii );
    M.add_transition( 2.*x, x, t, 0.4 );                       // 0.4 is inside element 0
    M.options.DISPLAY_LEVEL = 0;
    bool const ok = M.setup();
    std::cout << "    setup=" << ok << "  (refusal is the PASS here)\n";
    check_bool( "a tau off a boundary is refused, not ignored", !ok );
  }

  // ----------------------------------------------------------------- D: the LEFT LIMITS u^- and dx^-
  { std::cout << "\n--- D. the jump reads u(1-), NOT u(1+)\n";
    FFGraph G; DAESLV_IDAS M; M.set_dag( &G );
    FFVar t = G.add_var("t"), x = G.add_var("x(t)"), u = G.add_var("u");
    M.add_domain( t, FFDom( 0., 2., 2, FFDom::LGR, 4 ) );
    M.add_state( x, {t} ); M.add_input( u, {t}, FFDom::LGR, 1 );   // u = 1 on [0,1], 5 on [1,2]
    M.set_evolution_domain( t ); M.update_ref( x, 1. );
    M.add_equation( OpP( x, t ) + x - u, {t}, {TI}, io );
    M.add_equation( x - 1.,              {t}, {FFDom::LB}, ii );
    // x(1+) = x(1-) + u(1-).  x' = -x + u with x(0) = u = 1 gives x == 1 on [0,1], so x(1+) = 2 and
    // x(2) = 5 - 3 e^-1.  Had u^+ = 5 been bound instead, x(1+) = 6 and x(2) = 5 + e^-1 -- far apart,
    // so this cannot pass by accident.
    //
    // x'(tau^-) is NOT exercised: _check_transitions refuses a derivative operator in a transition
    // expression ("must be POINTWISE"), so no jump can read it today.  The placeholders are not minted
    // for that reason -- see base_dynamic.hpp.
    M.add_transition( x + u, x, t, 1. );
    M.options.SOLVE.RTOL = M.options.SOLVE.ATOL = 1e-10;
    M.options.DISPLAY_LEVEL = 0;
    bool const ok = M.setup();
    std::cout << "    setup=" << ok << ( ok? std::string(): "  :: " + M.extract_error() ) << "\n";
    check_bool( "a jump may read u(tau-)", ok );
    if( ok ){
      M.solve( std::vector<double>{ 1., 5. } );
      auto const& X = M.val_state();
      std::cout << "    x(1-)=" << std::scientific << std::setprecision(6) << X.at(1)[0]
                << "   so x(1+) should be 2 exactly\n";
      // on [1,2] with u = 5: x(t) = 5 + ( x(1+) - 5 ) e^-(t-1), so x(2) = 5 - 3 e^-1
      check( "x(2) = 5 - 3 e^-1  (u^- bound, not u^+)", X.back()[0], 5. - 3.*std::exp(-1.), 1e-7 );
    }
  }

  // ----------------------------------------------------------------- E: naming the post-jump state
  { std::cout << "\n--- E. the OPTIONAL post-jump state argument\n";
    auto build = [&]( DAESLV_IDAS& M, FFGraph& G, int mode ){
      FFVar t = G.add_var("t"), x = G.add_var("x(t)"), y = G.add_var("y(t)");
      M.set_dag( &G );
      M.add_domain( t, FFDom( 0., 2., 2, FFDom::LGR, 4 ) );
      M.add_state( x, {t} ); M.add_state( y, {t} );
      M.set_evolution_domain( t ); M.update_ref( x, 1. ); M.update_ref( y, 1. );
      M.add_equation( OpP( x, t ) + x, {t}, {TI}, io );
      M.add_equation( OpP( y, t ) + y, {t}, {TI}, io );
      M.add_equation( x - 1., {t}, {FFDom::LB}, ii );
      M.add_equation( y - 1., {t}, {FFDom::LB}, ii );
      // BOTH components depend on both x and y at tau^+, so the matching is legal but AMBIGUOUS
      std::vector<FFVar> L{ 2.*x, 3.*y }, R{ x + 0.*y, y + 0.*x };
      if( mode == 0 ) M.add_transition( L, R, t, 1. );
      else if( mode == 1 ) M.add_transition( L, R, t, 1., std::vector<FFVar>{ x, y } );
      else                 M.add_transition( L, R, t, 1., std::vector<FFVar>{ x, x } );
      M.options.SOLVE.RTOL = M.options.SOLVE.ATOL = 1e-10;
      M.options.DISPLAY_LEVEL = 0;
    };
    { FFGraph G; DAESLV_IDAS M; build( M, G, 1 );
      bool const ok = M.setup();
      check_bool( "naming {x,y} is accepted", ok );
      if( ok ){ M.solve( std::vector<double>{} );
        check( "x(2) = 2 e^-2 with the named matching", M.val_state().back()[0], 2.*std::exp(-2.), 1e-7 );
        check( "y(2) = 3 e^-2 with the named matching", M.val_state().back()[1], 3.*std::exp(-2.), 1e-7 ); } }
    { FFGraph G; DAESLV_IDAS M; build( M, G, 2 );
      std::cout << "    (naming x twice -- a refusal is the PASS)\n";
      check_bool( "naming the SAME state twice is refused", !M.setup() ); }
  }

  // ----------------------------------------------------------------- F: a derivative in the jump law
  { std::cout << "\n--- F. the jump reads x'(1-), the LEFT-limit derivative\n";
    FFGraph G; DAESLV_IDAS M; M.set_dag( &G );
    FFVar t = G.add_var("t"), x = G.add_var("x(t)"), u = G.add_var("u");
    M.add_domain( t, FFDom( 0., 2., 2, FFDom::LGR, 4 ) );
    M.add_state( x, {t} ); M.add_input( u, {t}, FFDom::LGR, 1 );     // u = 1 on [0,1], 5 on [1,2]
    M.set_evolution_domain( t ); M.update_ref( x, 2. );
    M.add_equation( OpP( x, t ) + x - u, {t}, {TI}, io );
    M.add_equation( x - 2.,              {t}, {FFDom::LB}, ii );
    // x(0)=2, u=1 on [0,1]  =>  x(t) = 1 + e^-t, so x(1-) = 1 + e^-1 and x'(1-) = -e^-1.
    // The jump is x(1+) = -x'(1-) = e^-1.  Were x'(1+) used instead, the law would read
    // x(1+) = -(-x(1+)+5) = x(1+)-5, which has NO solution -- so this cannot pass by the wrong side.
    M.add_transition( -OpP( x, t ), x, t, 1. );
    M.options.SOLVE.RTOL = M.options.SOLVE.ATOL = 1e-10;
    M.options.DISPLAY_LEVEL = 0;
    bool const ok = M.setup();
    std::cout << "    setup=" << ok << ( ok? std::string(): "  :: " + M.extract_error() ) << "\n";
    check_bool( "a jump may read x'(tau-)", ok );
    if( ok ){
      M.solve( std::vector<double>{ 1., 5. } );
      auto const& X = M.val_state();
      double const xm = 1. + std::exp(-1.), xp = std::exp(-1.);      // x(1-), then x(1+) = -x'(1-)
      std::cout << "    x(1-)=" << std::scientific << std::setprecision(6) << X.at(1)[0]
                << " (exact " << xm << ")   x(1+) should be " << xp << "\n";
      check( "x(1-) = 1 + e^-1", X.at(1)[0], xm, 1e-8 );
      // on [1,2] with u = 5: x(2) = 5 + ( x(1+) - 5 ) e^-1
      check( "x(2) = 5 + (e^-1 - 5) e^-1", X.back()[0], 5. + ( xp - 5. )*std::exp(-1.), 1e-7 );
    }
  }

  // ----------------------------------------------------------------- G: a derivative in the WRONG direction
  { std::cout << "\n--- G. differentiating along a NON-evolution direction must be refused\n";
    FFGraph G; DAESLV_IDAS M; M.set_dag( &G );
    FFVar t = G.add_var("t"), z = G.add_var("z"), x = G.add_var("x(t,z)");
    M.add_domain( t, FFDom( 0., 2., 2, FFDom::LGR, 4 ) );
    M.add_domain( z, FFDom( 0., 1., 1, FFDom::LGR, 3 ) );
    M.add_state( x, {t, z} ); M.set_evolution_domain( t ); M.update_ref( x, 1. );
    M.add_equation( OpP( x, t ) + x, {t, z}, {TI, FFDom::ALL}, io );
    M.add_equation( x - 1.,          {t, z}, {FFDom::LB, FFDom::ALL}, ii );
    M.add_transition( OpP( x, z ), x, t, 1. );                        // along z, not t
    M.options.DISPLAY_LEVEL = 0;
    std::cout << "    (a refusal is the PASS here)\n";
    check_bool( "a derivative along another direction is refused", !M.setup() );
  }

  std::cout << "\n==========  DAESLV_trans: " << n_pass << " passed, " << n_fail << " failed -- "
            << ( n_fail? "FAILURES": "ALL PASS" ) << "  ==========\n\n";
  return n_fail? 1: 0;
}
