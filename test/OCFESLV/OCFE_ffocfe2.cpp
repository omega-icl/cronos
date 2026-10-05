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
static void build( M& m )
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
  m.oc->setup();
}
int main()
{
  std::cout << std::scientific << std::setprecision(6)
            << "================================================================\n"
            << "  test13_ff: FFOCFESLV TWO-MAP form (map 1 {a,b}, map 2 {c})\n"
            << "================================================================\n";
  double const A = 2., B = 6., C = 0.5;
  double const G0 = A + (B+C)/A, G1 = B + C, J00 = 1.-(B+C)/(A*A), J01 = 1./A;
  for( int policy : { (int)FFOCFESLV::SHALLOW, (int)FFOCFESLV::COPY } ){
    std::string const pn = policy? "COPY": "SHALLOW";
    for( int lit : { 1, 0 } ){
      std::string const tag = pn + ( lit? " / c LITERAL": " / c DAG VARIABLE" );
      std::cout << "\n--- " << tag << "\n";
      M m; build( m );
      m.oc->register_control( m.a );                                   // the caller's registry: {a} only
      FFGraph rdag; FFVar pa = rdag.add_var("pa"), pb = rdag.add_var("pb"), pcv = rdag.add_var("pc");
      FFVar const cin = lit? FFVar( C ): pcv;
      FFOCFESLV op;  std::vector<FFVar> F;
      try{ F = op( { { m.a, {pa} }, { m.b, {pb} } }, { { m.c, {cin} } }, m.oc, policy, "t13" ); }
      catch( std::exception& e ){ check( false, tag + " embedding threw: " + e.what() ); delete m.oc; continue; }
      check( F.size() == 2 && op.n_input() == 3, tag + ": 2 outputs, 3 op inputs [a b | c]" );
      check( m.oc->controls().size() == ( policy? 1u: 2u ), tag + ( policy? ": the caller's registry {a} was RESTORED": ": the shared registry is now map 1 {a,b}" ) );
      std::vector<FFVar> vX{ pa, pb }; std::vector<double> vXv{ A, B };
      if( !lit ){ vX.push_back( pcv ); vXv.push_back( C ); }
      std::vector<double> Fv( 2 );
      try{ rdag.eval( F, Fv, vX, vXv ); } catch( std::exception& e ){ check( false, tag + " DAG eval threw: " + e.what() ); delete m.oc; continue; }
      { std::ostringstream o; o << tag << ": DAG value == closed form (worst rel " << std::max( rel(Fv[0],G0), rel(Fv[1],G1) ) << ")";
        check( rel(Fv[0],G0) < 1e-9 && rel(Fv[1],G1) < 1e-9, o.str() ); }
      std::vector<double> in3{ A, B, C }, Fd( 2 );                       // direct eval: [ map 1 | map 2 ]
      op.eval( 2u, Fd.data(), 3u, in3.data(), nullptr );
      { std::ostringstream o; o << tag << ": direct eval on [a b c] == closed form (worst rel " << std::max( rel(Fd[0],G0), rel(Fd[1],G1) ) << ")";
        check( rel(Fd[0],G0) < 1e-9 && rel(Fd[1],G1) < 1e-9, o.str() ); }
      if( !lit ){                                                         // c's value flows through the DAG
        vXv[2] = 2.*C; rdag.eval( F, Fv, vX, vXv );
        check( rel( Fv[1], B + 2.*C ) < 1e-9, tag + ": changing c in the DAG changes the value, same op" );
        vXv[2] = C;
      }
      auto const dF = rdag.FAD( F, vX ); std::vector<double> dFv( dF.size() ); rdag.eval( dF, dFv, vX, vXv );
      size_t const nX = vX.size();
      double w = std::max( std::max( rel( dFv[0], J00 ), rel( dFv[1], J01 ) ), std::max( std::fabs( dFv[nX+0] ), rel( dFv[nX+1], 1. ) ) );
      { std::ostringstream o; o << tag << ": dG/d(a,b) numeric (FFGradOCFESLV) == closed form (worst " << w << ")"; check( w < 1e-8, o.str() ); }
      if( !lit ){ double z = std::max( std::fabs( dFv[2] ), std::fabs( dFv[nX+2] ) );
        std::ostringstream o; o << tag << ": dG/dc (map 2) ZERO by contract (max |.| " << z << ")"; check( z == 0., o.str() ); }
      if( !policy ){
        m.oc->clear_controls(); bool threw = false; std::string msg;
        try{ rdag.eval( F, Fv, vX, vXv ); } catch( std::exception& e ){ threw = true; msg = e.what(); }
        check( threw && msg.find("registry changed") != std::string::npos, tag + ": evaluation refused after the caller changed the registry" );
      }
      delete m.oc;
    }
  }
  std::cout << "\n--- refusals and the one-map form\n";
  { M m; build( m ); FFGraph rdag; FFVar pa = rdag.add_var("pa"), pb = rdag.add_var("pb");
    auto refused = [&]( std::function<void()> f, std::string const& frag ){
      try{ f(); } catch( std::exception& e ){ return std::string( e.what() ).find( frag ) != std::string::npos; } return false; };
    FFOCFESLV op;
    check( refused( [&]{ op( { { m.a, {pa} }, { m.b, {pb} } }, m.oc ); }, "missing: c" ), "one-map form with c unmapped: refused, naming c" );
    check( refused( [&]{ op( { { m.a, {pa} } }, { { m.a, {pb} }, { m.b, {pb} }, { m.c, {pb} } }, m.oc ); }, "more than once" ), "an input in both maps: refused" );
    check( refused( [&]{ op( { { m.a, std::vector<FFVar>{ pa, pb } } }, { { m.b, {pb} }, { m.c, {pb} } }, m.oc ); }, "a has 1 DOFs" ), "a wrong DOF count: refused" );
    check( m.oc->controls().empty(), "refusals left the registry untouched" );
    FFOCFESLV op1; auto F1 = op1( 1, { { m.a, {pa} }, { m.b, {pb} } }, { { m.c, { FFVar( C ) } } }, m.oc, FFOCFESLV::COPY );
    std::vector<FFVar> one{ F1 }; std::vector<double> v1( 1 ); rdag.eval( one, v1, std::vector<FFVar>{ pa, pb }, std::vector<double>{ A, B } );
    check( rel( v1[0], G1 ) < 1e-9, "two-map, iDep=1: the single output == G1" );
    delete m.oc; }
  std::cout << "\n  test13_ff: " << npass << " passed, " << nfail << " failed -- " << (nfail? "FAILURES":"ALL PASS") << "\n";
  return nfail? 1: 0;
}
