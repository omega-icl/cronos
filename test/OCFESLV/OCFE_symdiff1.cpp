// test16_ff.cpp -- OCFESLV::fdiff: the forward-differentiated model, legacy layout f[i*nf+j] = dF_j/d(dir_i).
//  (1) algebraic, closed form: x1 - a = 0, x1 x2 - (b + K c) = 0; G0 = x1 + x2, G1 = x1 x2; K a CONSTANT.
//      directions (a,0), (c,0), (K,0): an input, another input, a constant.
//  (2) dynamic: dx/dt = -k x, x(0) = x0; G0 = x(T), G1 = x(T/2).  The product (differentiate, then collocate) is
//      checked against solve_fsens on the original (collocate, then differentiate) and against the closed form.
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include "ocfeslv.hpp"
#include "ffocfe.hpp"
using namespace mc;
static int npass = 0, nfail = 0;
static void check( bool c, std::string const& w ){ std::cout << "  " << (c? "PASS  ": "FAIL  ") << w << "\n"; c? ++npass: ++nfail; }
static double rel( double a, double b ){ double s=std::fabs(b); return std::fabs(a-b)/( s>1e-12? s: 1. ); }
int main()
{
  std::cout << std::scientific << std::setprecision(6)
            << "================================================================\n"
            << "  test16_ff: OCFESLV::fdiff -- product outputs f[i*nf+j] = dF_j/d(dir_i)\n"
            << "================================================================\n";
  { std::cout << "\n--- (1) algebraic, closed form; directions a, c (inputs) and K (constant)\n";
    FFGraph DAG; FFVar x1 = DAG.add_var("x1"), x2 = DAG.add_var("x2"), a = DAG.add_var("a"), b = DAG.add_var("b"),
                         c = DAG.add_var("c"), K = DAG.add_var("K");
    OCFESLV oc( &DAG );
    oc.add_state( x1, {} ); oc.add_state( x2, {} ); oc.add_input( a, {} ); oc.add_input( b, {} ); oc.add_input( c, {} );
    oc.FFModel::set_constant( { K }, { 1.5 } );
    oc.update_ref( x1, 1.0 ); oc.update_ref( x2, 1.0 );
    OCFESLV::EqnOptions alg( OCFESLV::EqnRole::INTERIOR, 0 );
    oc.add_equation( x1 - a, {}, {}, alg ); oc.add_equation( x1 * x2 - ( b + K * c ), {}, {}, alg );
    oc.add_output( x1 + x2, std::vector<FFVar>{}, std::vector<double>{} ); oc.add_output( x1 * x2, std::vector<FFVar>{}, std::vector<double>{} );
    oc.options.SOLVE.MAX_ITER = 60; oc.options.SOLVE.RES_TOL = 1e-12; oc.options.DISPLAY_LEVEL = 0;
    oc.setup();
    double const A = 2., B = 6., C = 0.5, Kv = 1.5, S = B + Kv*C;
    double const J[3][2] = { { 1. - S/(A*A), 0. }, { Kv/A, Kv }, { C/A, C } };     // [dir][fct]: d/da, d/dc, d/dK
    std::string err;
    OCFESLV* P = oc.fdiff( { { a, 0 }, { c, 0 }, { K, 0 } }, err );
    check( P != nullptr, "fdiff along (a, c, K) returns a product" + ( err.empty()? std::string(): " -- " + err ) );
    if( P ){
      bool const ok = P->setup();
      check( ok && P->n_colloc_fct() == 6, "product sets up with 3 x 2 = 6 outputs" );
      std::vector<double> xv, inp; P->init( xv, inp, nullptr );
      P->set_input_values( a, {A}, inp.data() ); P->set_input_values( b, {B}, inp.data() ); P->set_input_values( c, {C}, inp.data() );
      auto const rep = P->solve( xv.data(), inp.data(), &Kv );
      auto const f = P->val_functions();
      double w = 0.; for( size_t i=0;i<3;++i ) for( size_t j=0;j<2;++j ) w = std::max( w, rel( f[i*2+j], J[i][j] ) );
      std::ostringstream o; o << "f[i*nf+j] == closed-form dG_j/d(a, c, K)  (worst rel " << w << ")";
      check( rep.converged && f.size() == 6 && w < 1e-9, o.str() );
      delete P;
    }
    check( oc.fdiff( { { x1, 0 } }, err ) == nullptr && err.find("neither") != std::string::npos, "a state as direction: refused" );
    check( oc.fdiff( { { a, 1 } }, err ) == nullptr && err.find("out of range") != std::string::npos, "DOF out of range: refused" );
  }
  { std::cout << "\n--- (2) dynamic dx/dt = -k x: product vs solve_fsens (discrete) and vs closed form\n";
    double const T = 2., k0 = 0.7, x00 = 1.0;
    FFGraph DAG; FFVar t = DAG.add_var("t"), x = DAG.add_var("x(t)"), k = DAG.add_var("k"), x0 = DAG.add_var("x0");
    OCFESLV oc( &DAG );
    oc.add_domain( t, FFDom( 0., T, 6, FFDom::LGR, 5 ) ); oc.set_evolution_domain( t );
    oc.add_state( x, {t} ); oc.add_input( k, {} ); oc.add_input( x0, {} );
    oc.update_ref( x, []( OCFESLV::t_Coord const& ){ return 1.0; } );                // a FUNCTION reference
    FFPartial OpP;
    oc.add_equation( OpP( x, t ) + k * x, {t}, { FFDom::ALL - FFDom::LB }, OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
    oc.add_equation( x - x0, {t}, { FFDom::LB }, OCFESLV::EqnOptions( OCFESLV::EqnRole::INITIAL, 0 ) );
    oc.add_output( x, {t}, { T } ); oc.add_output( x, {t}, { 0.5*T } );
    oc.options.SOLVE.MAX_ITER = 60; oc.options.SOLVE.RES_TOL = 1e-12; oc.options.DISPLAY_LEVEL = 0; oc.options.SOLVE.MARCHING = false;
    oc.setup();
    std::vector<double> xv, inp; oc.init( xv, inp, nullptr );
    oc.set_input_values( k, {k0}, inp.data() ); oc.set_input_values( x0, {x00}, inp.data() );
    oc.register_control( k ); oc.register_control( x0 );
    std::vector<double> xs( xv ); oc.solve_fsens( xs.data(), inp.data(), nullptr );
    auto const Jd = oc.sens_jacobian();                                              // [fct*ncd + ctrl]
    size_t const ik = oc.controls().at( k ).offset, ix = oc.controls().at( x0 ).offset;
    std::string err;
    OCFESLV* P = oc.fdiff( { { k, 0 }, { x0, 0 } }, err );
    check( P != nullptr, "fdiff along (k, x0) returns a product" + ( err.empty()? std::string(): " -- " + err ) );
    if( P ){
      check( P->setup() && P->n_colloc_fct() == 4, "product sets up with 2 x 2 = 4 outputs" );
      std::vector<double> pv, pinp; P->init( pv, pinp, nullptr );
      P->set_input_values( k, {k0}, pinp.data() ); P->set_input_values( x0, {x00}, pinp.data() );
      auto const rep = P->solve( pv.data(), pinp.data(), nullptr );
      auto const f = P->val_functions();
      double wd = 0., wa = 0.;
      double const an[2][2] = { { -T*x00*std::exp(-k0*T), -0.5*T*x00*std::exp(-k0*0.5*T) }, { std::exp(-k0*T), std::exp(-k0*0.5*T) } };
      for( size_t j=0;j<2;++j ){ wd = std::max( wd, rel( f[0*2+j], Jd[j*2+ik] ) ); wd = std::max( wd, rel( f[1*2+j], Jd[j*2+ix] ) );
                                  wa = std::max( wa, rel( f[0*2+j], an[0][j] ) );  wa = std::max( wa, rel( f[1*2+j], an[1][j] ) ); }
      { std::ostringstream o; o << "product == solve_fsens on the original (worst rel " << wd << ")"; check( rep.converged && wd < 1e-9, o.str() ); }
      { std::ostringstream o; o << "product == closed form, spectral accuracy (worst rel " << wa << ")"; check( wa < 1e-5, o.str() ); }
      delete P;
    }
  }
  { std::cout << "\n--- (3) FFOCFESLV SYMDIFF: map 1 {a,b}, map 2 {c, K} as DAG variables\n";
    FFGraph DAG; FFVar x1 = DAG.add_var("x1"), x2 = DAG.add_var("x2"), a = DAG.add_var("a"), b = DAG.add_var("b"),
                         c = DAG.add_var("c"), K = DAG.add_var("K");
    OCFESLV oc( &DAG );
    oc.add_state( x1, {} ); oc.add_state( x2, {} ); oc.add_input( a, {} ); oc.add_input( b, {} ); oc.add_input( c, {} );
    oc.FFModel::set_constant( { K }, { 1.5 } );
    oc.update_ref( x1, 1.0 ); oc.update_ref( x2, 1.0 );
    OCFESLV::EqnOptions alg( OCFESLV::EqnRole::INTERIOR, 0 );
    oc.add_equation( x1 - a, {}, {}, alg ); oc.add_equation( x1 * x2 - ( b + K * c ), {}, {}, alg );
    oc.add_output( x1 + x2, std::vector<FFVar>{}, std::vector<double>{} ); oc.add_output( x1 * x2, std::vector<FFVar>{}, std::vector<double>{} );
    oc.options.SOLVE.MAX_ITER = 60; oc.options.SOLVE.RES_TOL = 1e-12; oc.options.DISPLAY_LEVEL = 0;
    oc.setup();
    double const A = 2., B = 6., C = 0.5, Kv = 1.5, S = B + Kv*C;
    // closed form [output][a, b, c, K]
    double const J[2][4] = { { 1. - S/(A*A), 1./A, Kv/A, C/A }, { 0., 1., Kv, C } };
    FFGraph rdag; FFVar pa = rdag.add_var("pa"), pb = rdag.add_var("pb"), pc = rdag.add_var("pc"), pK = rdag.add_var("pK");
    FFOCFESLV op;
    std::vector<FFVar> F = op( { { a, {pa} }, { b, {pb} } }, { { c, {pc} }, { K, {pK} } }, &oc, FFOCFESLV::COPY );
    std::vector<FFVar> vX{ pa, pb, pc, pK }; std::vector<double> vXv{ A, B, C, Kv };
    auto run = [&]( std::vector<FFVar> const& sym ){
      FFOCFESLV::options.SYMDIFF = sym;
      auto const dF = rdag.FAD( F, vX ); std::vector<double> d( dF.size() ); rdag.eval( dF, d, vX, vXv );
      FFOCFESLV::options.SYMDIFF.clear(); return d; };                              // d[j*4 + col]
    try{
      auto const dN = run( {} );
      double wn = 0., zn = 0.; for( size_t j=0;j<2;++j ){ for( size_t u=0;u<2;++u ) wn = std::max( wn, std::fabs( dN[j*4+u] - J[j][u] ) );
                                                          for( size_t u=2;u<4;++u ) zn = std::max( zn, std::fabs( dN[j*4+u] ) ); }
      { std::ostringstream o; o << "numeric: dG/d(a,b) closed form (" << wn << "), map-2 columns zero (" << zn << ")"; check( wn < 1e-9 && zn == 0., o.str() ); }
      auto const dC = run( { pc, pK } );
      double wc = 0., zc = 0.; for( size_t j=0;j<2;++j ){ for( size_t u=2;u<4;++u ) wc = std::max( wc, std::fabs( dC[j*4+u] - J[j][u] ) );
                                                          for( size_t u=0;u<2;++u ) zc = std::max( zc, std::fabs( dC[j*4+u] ) ); }
      { std::ostringstream o; o << "SYMDIFF={c, K}: map-2 input and CONSTANT differentiated symbolically (" << wc << "); a, b columns zero (" << zc << ")";
        check( wc < 1e-9 && zc == 0., o.str() ); }
      auto const dA = run( { pa, pb, pc, pK } );
      double wa = 0.; for( size_t j=0;j<2;++j ) for( size_t u=0;u<4;++u ) wa = std::max( wa, std::fabs( dA[j*4+u] - J[j][u] ) );
      { std::ostringstream o; o << "SYMDIFF=all four: the full Jacobian, closed form (" << wa << ")"; check( wa < 1e-9, o.str() ); }
    } catch( std::exception& e ){ check( false, std::string( "FFOCFESLV SYMDIFF threw: " ) + e.what() ); }
  }
  std::cout << "\n  test16_ff: " << npass << " passed, " << nfail << " failed -- " << (nfail? "FAILURES":"ALL PASS") << "\n";
  return nfail? 1: 0;
}
