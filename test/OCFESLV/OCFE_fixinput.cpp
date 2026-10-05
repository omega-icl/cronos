// test13_ff.cpp -- FFOCFESLV TWO-MAP contract on a closed-form model (AE1 + a non-control input c).
//   states x1,x2; inputs a,b (map 1) and c (map 2):  F1: x1 - a = 0,  F2: x1 x2 - (b + c) = 0
//   outputs G0 = x1 + x2 = a + (b+c)/a,  G1 = x1 x2 = b + c
//   dG0/da = 1-(b+c)/a^2, dG0/db = 1/a, dG1/da = 0, dG1/db = 1;  map-2 column (c) ZERO by contract
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <functional>
#include "ocfeslv.hpp"
#include "ffocfe.hpp"
using namespace mc;
static int npass = 0, nfail = 0;
static void check( bool c, std::string const& w ){ std::cout << "  " << (c? "PASS  ": "FAIL  ") << w << "\n"; c? ++npass: ++nfail; }
static double rel( double a, double b ){ double s=std::fabs(b); return std::fabs(a-b)/( s>1e-12? s: 1. ); }
struct M { FFGraph DAG; FFVar x1, x2, a, b, c; OCFESLV* oc = nullptr; };
static void build( M& m, std::vector<double> const& pre_fix = {} )
{
  m.x1 = m.DAG.add_var("x1"); m.x2 = m.DAG.add_var("x2"); m.a = m.DAG.add_var("a"); m.b = m.DAG.add_var("b"); m.c = m.DAG.add_var("c");
  m.oc = new OCFESLV( &m.DAG );
  m.oc->add_state( m.x1, {} ); m.oc->add_state( m.x2, {} );
  m.oc->add_input( m.a, {} ); m.oc->add_input( m.b, {} ); m.oc->add_input( m.c, {} );
  m.oc->update_ref( m.x1, 1.0 ); m.oc->update_ref( m.x2, 1.0 );
  OCFESLV::EqnOptions alg( OCFESLV::EqnRole::INTERIOR, 0 );
  m.oc->add_equation( m.x1 - m.a, {}, {}, alg );
  m.oc->add_equation( m.x1 * m.x2 - ( m.b + m.c ), {}, {}, alg );
  m.oc->add_output( m.x1 + m.x2, std::vector<FFVar>{}, std::vector<double>{} );
  m.oc->add_output( m.x1 * m.x2, std::vector<FFVar>{}, std::vector<double>{} );
  m.oc->options.SOLVE.MAX_ITER = 60; m.oc->options.SOLVE.RES_TOL = 1e-12; m.oc->options.DISPLAY_LEVEL = 0;
  if( !pre_fix.empty() ) m.oc->fix_input( m.c, pre_fix );
  m.oc->setup();
}
int main()
{
  std::cout << std::scientific << std::setprecision(6)
            << "================================================================\n"
            << "  test15_ff: OCFESLV honours fix_input -- the MODEL supplies c's values\n"
            << "================================================================\n";
  double const A = 2., B = 6.;
  auto G0 = [&]( double c ){ return A + (B+c)/A; };  auto G1 = [&]( double c ){ return B + c; };
  for( int when : { 0, 1 } ){
    std::string const tag = when? "fixed AFTER setup": "fixed BEFORE setup";
    std::cout << "\n--- " << tag << " (c = 0.5)\n";
    M m; if( when ){ build( m ); m.oc->fix_input( m.c, { 0.5 } ); } else build( m, { 0.5 } );
    check( m.oc->fixed_input( m.c ) != nullptr, tag + ": c is fixed" );
    std::vector<double> xv, inp; m.oc->init( xv, inp, nullptr );
    check( m.oc->get_input_values( m.c, inp.data() ) == std::vector<double>{ 0.5 }, tag + ": init() writes c's fixed value into its slot" );
    m.oc->set_input_values( m.a, {A}, inp.data() ); m.oc->set_input_values( m.b, {B}, inp.data() );
    m.oc->set_input_values( m.c, {999.}, inp.data() );                   // the caller's value must be IGNORED
    std::vector<double> x1( xv ); auto const rep = m.oc->solve( x1.data(), inp.data(), nullptr );
    auto const F = m.oc->val_functions();
    { std::ostringstream o; o << tag << ": solve with 999 in c's slot == closed form at c=0.5 (worst rel "
        << std::max( rel(F[0],G0(.5)), rel(F[1],G1(.5)) ) << ")";
      check( rep.converged && rel(F[0],G0(.5)) < 1e-9 && rel(F[1],G1(.5)) < 1e-9, o.str() ); }
    check( m.oc->get_input_values( m.c, inp.data() ) == std::vector<double>{ 999. }, tag + ": the caller's inp itself is left untouched" );
    m.oc->register_control( m.a ); m.oc->register_control( m.b );
    for( int adj : { 0, 1 } ){
      std::vector<double> x2( xv );
      bool const ok = adj? m.oc->solve_asens( x2.data(), inp.data(), nullptr ): m.oc->solve_fsens( x2.data(), inp.data(), nullptr );
      auto const& J = m.oc->sens_jacobian();
      double const w = ok && J.size() >= 4? std::max( std::max( rel( J[0], 1.-(B+.5)/(A*A) ), rel( J[1], 1./A ) ),
                                                     std::max( std::fabs( J[2] ), rel( J[3], 1. ) ) ): 1.;
      std::ostringstream o; o << tag << ( adj? ": solve_asens": ": solve_fsens" ) << " dG/d(a,b) == closed form at c=0.5 (worst " << w << ")";
      check( w < 1e-8, o.str() );
    }
    if( when ){                                                          // re-fix after setup: applied at the next call
      m.oc->fix_input( m.c, { 0.25 } );
      std::vector<double> x3( xv ); m.oc->solve( x3.data(), inp.data(), nullptr );
      check( rel( m.oc->val_functions()[1], G1(.25) ) < 1e-9, tag + ": fixed again at 0.25 -- the next solve uses it" );
    }
    FFGraph rdag; FFVar pa = rdag.add_var("pa"), pb = rdag.add_var("pb");
    FFOCFESLV op; std::vector<FFVar> Fo;
    try{ Fo = op( { { m.a, {pa} }, { m.b, {pb} } }, m.oc, FFOCFESLV::COPY ); }
    catch( std::exception& e ){ check( false, tag + ": FFOCFESLV one-map threw: " + e.what() ); delete m.oc; continue; }
    double const cf = when? .25: .5;
    std::vector<double> Fv( 2 ); rdag.eval( Fo, Fv, std::vector<FFVar>{ pa, pb }, std::vector<double>{ A, B } );
    { std::ostringstream o; o << tag << ": FFOCFESLV one-map {a,b}, c left to the model == closed form (worst rel "
        << std::max( rel(Fv[0],G0(cf)), rel(Fv[1],G1(cf)) ) << ")";
      check( rel(Fv[0],G0(cf)) < 1e-9 && rel(Fv[1],G1(cf)) < 1e-9, o.str() ); }
    delete m.oc;
  }
  std::cout << "\n--- a model with NOTHING fixed is untouched\n";
  { M m; build( m ); std::vector<double> xv, inp; m.oc->init( xv, inp, nullptr );
    m.oc->set_input_values( m.a, {A}, inp.data() ); m.oc->set_input_values( m.b, {B}, inp.data() ); m.oc->set_input_values( m.c, {0.75}, inp.data() );
    m.oc->solve( xv.data(), inp.data(), nullptr );
    check( rel( m.oc->val_functions()[1], G1(.75) ) < 1e-9, "unfixed c: the caller's value is used" );
    delete m.oc; }
  std::cout << "\n  test15_ff: " << npass << " passed, " << nfail << " failed -- " << (nfail? "FAILURES":"ALL PASS") << "\n";
  return nfail? 1: 0;
}
