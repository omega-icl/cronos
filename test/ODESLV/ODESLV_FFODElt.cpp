// FFODE_lt.cpp -- identity of FFODESLV operations in the DAG (FFBaseODESLV::lt, 2026-10-01).
// Calls that differ only by the OUTPUT INDEX (idep) must share ONE operation; embeddings of the same solver over the
// same variables that differ by the MAP SPLIT (one-map vs two-map) or the POLICY (COPY vs SHALLOW) must be distinct
// operations, each with its own derivatives.  Before 2026-10-01 FFODESLV defined no lt: the DAG compared only the
// operands and the solver pointer, so the second of a one-map / two-map pair returned the first's node and
// derivatives.  Model: x' = -p q x, x(0) = 1, outputs x(0.5) and x(1); closed-form derivatives.
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>
#include "ffode.hpp"

static int npass = 0, nfail = 0;
static void check( bool c, std::string const& w ){ std::printf( "  %s  %s\n", c? "PASS": "FAIL", w.c_str() ); c? ++npass: ++nfail; }
static bool same( std::vector<mc::FFVar> const& a, std::vector<mc::FFVar> const& b )
{ if( a.size() != b.size() ) return false; for( size_t i=0; i<a.size(); ++i ) if( a[i].id() != b[i].id() ) return false; return true; }

double const Pv = 0.7, Qv = 1.;
struct Model {
  mc::FFVar t, x, p, q;
  mc::ODESLVS_CVODES* S;
  Model( mc::FFGraph& G, int k ){
    mc::FFPartial OpP; mc::FFEval OpE;
    t = G.add_var( "t"+std::to_string(k) ); x = G.add_var( "x"+std::to_string(k) );
    p = G.add_var( "p"+std::to_string(k) ); q = G.add_var( "q"+std::to_string(k) );
    S = new mc::ODESLVS_CVODES( &G );
    S->add_domain( t, mc::FFDom( 0., 1., 2, mc::FFDom::LGR, 3 ) );  S->set_evolution_domain( t );
    S->add_state( x, {t} );  S->update_ref( x, 1. );  S->add_input( p );  S->add_input( q );
    int const T_INT = mc::FFDom::ALL - mc::FFDom::LB;
    S->add_equation( OpP( x, t ) + p * q * x, {t}, {T_INT}, mc::FFModel::EqnOptions( mc::FFModel::EqnRole::INTERIOR, 0 ) );
    S->add_equation( x - 1.,                  {t}, {mc::FFDom::LB}, mc::FFModel::EqnOptions( mc::FFModel::EqnRole::INITIAL, 0 ) );
    S->add_output( OpE( x, t, .5 ) );  S->add_output( OpE( x, t, 1. ) );
    S->options.RTOL = S->options.ATOL = S->options.RTOLS = S->options.ATOLS = S->options.RTOLB = S->options.ATOLB = 1e-10;
    S->options.DISPLAY = 0;
    S->setup();
  }
  ~Model(){ delete S; }
};
typedef std::vector<mc::FFModel::InputArg> t_Map;
// d y_k / d (P,Q) at (Pv,Qv), row-major [k][i]
static std::vector<double> jac( mc::FFGraph& G, std::vector<mc::FFVar> const& y, std::vector<mc::FFVar> const& X ){
  auto const dF = G.FAD( y, X );  std::vector<double> v( dF.size() );  G.eval( dF, v, X, std::vector<double>{ Pv, Qv } );  return v; }
static bool close( std::vector<double> const& a, std::vector<double> const& b ){
  if( a.size() != b.size() ) return false; for( size_t i=0; i<a.size(); ++i ) if( std::fabs( a[i]-b[i] ) > 1e-7 ) return false; return true; }

int main(){
  std::vector<double> JF, J2;                                  // one-map (all differentiated) / two-map ({p} | {q})
  for( double s : { .5, 1. } ){ double const e = std::exp( -s*Pv*Qv ); JF.push_back( -s*Qv*e ); JF.push_back( -s*Pv*e ); }
  for( size_t k=0; k<2; ++k ){ J2.push_back( JF[2*k] ); J2.push_back( 0. ); }
  mc::FFODESLV OpODE;
  { // idep sharing
    mc::FFGraph G; Model M( G, 0 ); mc::FFVar P = G.add_var( "P" ), Q = G.add_var( "Q" ); std::vector<mc::FFVar> X{ P, Q };
    t_Map m{ { M.p, std::vector<mc::FFVar>{ P } }, { M.q, std::vector<mc::FFVar>{ Q } } };
    mc::FFVar const y0 = OpODE( 0, m, M.S ), y1 = OpODE( 1, m, M.S );
    std::vector<mc::FFVar> const ya = OpODE( m, M.S );
    check( same( { y0, y1 }, ya ), "idep: OpODE(0,m), OpODE(1,m), OpODE(m) share ONE operation" );
    check( close( jac( G, ya, X ), JF ), "  ... with the one-map derivatives" );
  }
  for( int order : { 0, 1 } ){ // one-map and two-map of the same solver, in both orders
    mc::FFGraph G; Model M( G, 0 ); mc::FFVar P = G.add_var( "P" ), Q = G.add_var( "Q" ); std::vector<mc::FFVar> X{ P, Q };
    t_Map m1{ { M.p, std::vector<mc::FFVar>{ P } }, { M.q, std::vector<mc::FFVar>{ Q } } };
    t_Map d2{ { M.p, std::vector<mc::FFVar>{ P } } }, r2{ { M.q, std::vector<mc::FFVar>{ Q } } };
    std::vector<mc::FFVar> y1, y2;
    if( order == 0 ){ y1 = OpODE( m1, M.S ); y2 = OpODE( d2, r2, M.S ); }
    else            { y2 = OpODE( d2, r2, M.S ); y1 = OpODE( m1, M.S ); }
    std::string const o = order? "two-map THEN one-map": "one-map THEN two-map";
    check( !same( y1, y2 ), o + ": two operations" );
    check( close( jac( G, y1, X ), JF ) && close( jac( G, y2, X ), J2 ), o + ": each with its own derivatives" );
    if( order ){ mc::FFVar const z0 = OpODE( 0, d2, r2, M.S ); check( z0.id() == y2[0].id(), "idep with the two-map form shares the two-map operation" ); }
  }
  { // policies
    mc::FFGraph G; Model M( G, 0 ); mc::FFVar P = G.add_var( "P" ), Q = G.add_var( "Q" );
    t_Map m{ { M.p, std::vector<mc::FFVar>{ P } }, { M.q, std::vector<mc::FFVar>{ Q } } };
    std::vector<mc::FFVar> const yc = OpODE( m, M.S, mc::FFODESLV::COPY ), ys = OpODE( m, M.S, mc::FFODESLV::SHALLOW );
    check( !same( yc, ys ), "COPY and SHALLOW of the same solver: two operations" );
  }
  { // different solvers
    mc::FFGraph G; Model M( G, 0 ), M2( G, 1 ); mc::FFVar P = G.add_var( "P" ), Q = G.add_var( "Q" );
    t_Map m { { M.p,  std::vector<mc::FFVar>{ P } }, { M.q,  std::vector<mc::FFVar>{ Q } } };
    t_Map m2{ { M2.p, std::vector<mc::FFVar>{ P } }, { M2.q, std::vector<mc::FFVar>{ Q } } };
    check( !same( OpODE( m, M.S ), OpODE( m2, M2.S ) ), "different solvers: two operations" );
  }
  std::printf( "\n  FFODE_lt: %d passed, %d failed -- %s\n", npass, nfail, nfail? "FAILURES": "ALL PASS" );
  return nfail? 1: 0;
}
