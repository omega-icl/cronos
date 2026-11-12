// DAESLV_ffdae -- FFDAESLV: the DAE solver as a DAG operation.
// =====================================================================================================
// FFDAESLV is FFDE<DAESLV_IDAS>, the same template FFODESLV is an instantiation of.  So the test is
// the one FFODESLV answers: embed the solver in a DAG, evaluate the DAG, differentiate it -- and get
// what the solver gives directly.
//
//   x' = -k x , x(0) = a ;  f = int_0^1 x dt = a(1-e^-k)/k
#include <cmath>
#include <iomanip>
#include <iostream>
#include "ffde.hpp"

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
static void check_bool( std::string const& w, bool ok )
{
  ok? ++n_pass: ++n_fail;
  std::cout << "    " << std::left << std::setw( 40 ) << w << ( ok? "  PASS\n": "  FAIL\n" );
}

int main()
{
  FFPartial OpP;  FFIntegral OpI;
  int const TI = FFDom::ALL - FFDom::LB;
  auto io = FFModel::EqnOptions( FFModel::EqnRole::INTERIOR, 0 );
  auto ii = FFModel::EqnOptions( FFModel::EqnRole::INITIAL,  0 );

  FFGraph G;  DAESLV_IDAS M;  M.set_dag( &G );
  FFVar t = G.add_var("t"), x = G.add_var("x(t)"), k = G.add_var("k"), a = G.add_var("a");
  M.add_domain( t, FFDom( 0., 1., 1, FFDom::LGR, 5 ) );
  M.add_state( x, {t} );  M.add_input( k );  M.add_input( a );
  M.set_evolution_domain( t );  M.update_ref( x, 1. );
  M.add_equation( OpP( x, t ) + k*x, {t}, {TI}, io );
  M.add_equation( x - a,             {t}, {FFDom::LB}, ii );
  M.add_output( OpI( std::vector<FFVar>{ x }, { t } )[0] );
  M.options.SOLVE.RTOL = M.options.SOLVE.ATOL = 1e-11;
  M.options.DISPLAY_LEVEL = 0;
  if( !M.setup() ){ std::cout << "setup failed: " << M.extract_error() << "\n"; return 1; }

  double const K = 1.3, A = 2.1, e = std::exp(-K);
  double const F  = A*( 1.-e )/K, dFda = ( 1.-e )/K, dFdk = A*( K*e - (1.-e) )/(K*K);

  // ------------------------------------------------------------------ A: embed and evaluate
  std::cout << "\n--- A. FFDAESLV embedded in the DAG, evaluated\n";
  FFDAESLV OpDAE;
  std::cout << "    op name tag: " << OpDAE.name() << "\n";
  check_bool( "the op reports a DAE tag", OpDAE.name().rfind( "DAE[", 0 ) == 0 );
  std::vector<FFVar> vPar{ k, a };
  // an InputArg pairs a model INPUT with the DAG variable(s) that supply it -- here the same FFVar
  std::vector<FFModel::InputArg> vArg{ { k, std::vector<FFVar>{ k } }, { a, std::vector<FFVar>{ a } } };
  std::vector<FFVar> vF = OpDAE( vArg, &M, FFDAESLV::SHALLOW );
  check_bool( "the embedding returned outputs", !vF.empty() );
  if( vF.empty() ) return 1;
  std::vector<double> dP{ K, A }, dF( M.nf(), 0. );
  G.eval( vF.size(), vF.data(), dF.data(), vPar.size(), vPar.data(), dP.data() );
  check( "DAG eval = int_0^1 x dt", dF[0], F, 1e-8 );

  // ------------------------------------------------------------------ B: differentiate the DAG
  std::cout << "\n--- B. the DAG differentiated through the operation\n";
  std::vector<FFVar> const dFdP = G.FAD( vF, vPar );
  std::vector<double> dG( dFdP.size(), 0. );
  G.eval( dFdP.size(), dFdP.data(), dG.data(), vPar.size(), vPar.data(), dP.data() );
  check( "d f/dk through the DAG", dG[0], dFdk, 1e-7 );
  check( "d f/da through the DAG", dG[1], dFda, 1e-7 );

  // the same numbers the solver gives directly, by solve_fsens
  M.solve_fsens( dP );
  auto const& Gs = M.val_function_gradient();
  size_t const b = M.blk_fct(0).first, jk = M.parameter_index( k )[0], ja = M.parameter_index( a )[0];
  check( "  == solve_fsens on d/dk", dG[0], Gs[jk][b], 1e-9 );
  check( "  == solve_fsens on d/da", dG[1], Gs[ja][b], 1e-9 );

  // ------------------------------------------------------------------ C: a model with a TRANSITION, SYMDIFF
  // x(0.5+) = x(0.5-) + d: the numeric route takes the forward sensitivities across the jump; the symbolic route
  // (SYMDIFF) is fdiff's, whose product carries the transition and its derivative per direction (refused until 12 Nov 2026)
  std::cout << "\n--- C. a model with a transition: numeric and symbolic (SYMDIFF) derivatives\n";
  { FFGraph GT;  DAESLV_IDAS MT;  MT.set_dag( &GT );
    FFVar tt = GT.add_var("t"), xt = GT.add_var("x(t)"), kt = GT.add_var("k"), at = GT.add_var("a"), dt = GT.add_var("d");
    FFEval OpE;
    MT.add_domain( tt, FFDom( 0., 1., 2, FFDom::LGR, 5 ) );
    MT.add_state( xt, {tt} );  MT.add_input( kt );  MT.add_input( at );  MT.add_input( dt );
    MT.set_evolution_domain( tt );  MT.update_ref( xt, 1. );
    MT.add_equation( OpP( xt, tt ) + kt*xt, {tt}, {TI}, io );
    MT.add_equation( xt - at,               {tt}, {FFDom::LB}, ii );
    MT.add_transition( xt + dt, xt, tt, 0.5 );
    MT.add_output( OpE( xt, tt, 1. ) );
    MT.options.SOLVE.RTOL = MT.options.SOLVE.ATOL = 1e-11;  MT.options.DISPLAY_LEVEL = 0;
    bool const su = MT.setup();
    check_bool( "setup of the model with a transition", su );
    if( su ){
      double const DV = 0.7;  std::vector<FFVar> vP{ kt, at, dt };  std::vector<double> dPT{ K, A, DV };
      std::vector<FFModel::InputArg> vA{ { kt, std::vector<FFVar>{ kt } }, { at, std::vector<FFVar>{ at } }, { dt, std::vector<FFVar>{ dt } } };
      std::vector<FFVar> vFT = OpDAE( vA, &MT, FFDAESLV::SHALLOW );
      // x(1) = ( a e^{-k/2} + d ) e^{-k/2}
      double const eh = std::exp( -0.5*K ), dFdd = eh;
      FFDAESLV::options.SYMDIFF.clear();
      std::vector<FFVar> const dFT = GT.FAD( vFT, vP );
      std::vector<double> dGT( dFT.size(), 0. );
      GT.eval( dFT.size(), dFT.data(), dGT.data(), vP.size(), vP.data(), dPT.data() );
      check( "numeric d x(1)/dd through the DAG", dGT[2], dFdd, 1e-7 );
      std::vector<double> dST( 3, 1e30 );  bool sym = true;
      FFDAESLV::options.SYMDIFF = vP;
      try{ std::vector<FFVar> const dFS = GT.FAD( vFT, vP );
           GT.eval( dFS.size(), dFS.data(), dST.data(), vP.size(), vP.data(), dPT.data() ); }
      catch( std::exception const& e ){ sym = false;  std::cout << "    SYMDIFF threw: " << e.what() << "\n"; }
      FFDAESLV::options.SYMDIFF.clear();
      check_bool( "SYMDIFF differentiates across the transition (was refused)", sym );
      check( "symbolic d x(1)/dk", dST[0], -A*eh*eh - 0.5*DV*eh, 1e-7 );
      check( "symbolic d x(1)/da", dST[1], eh*eh,                1e-7 );
      check( "symbolic d x(1)/dd", dST[2], dFdd,                 1e-7 );
      check( "  == the numeric route on d/dk", dST[0], dGT[0], 1e-8 );
    } }

  std::cout << "\n==========  DAESLV_ffdae: " << n_pass << " passed, " << n_fail << " failed -- "
            << ( n_fail? "FAILURES": "ALL PASS" ) << "  ==========\n\n";
  return n_fail? 1: 0;
}
