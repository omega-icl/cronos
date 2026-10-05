// OCFE_pencil.cpp -- what CRONOS says today about (1) complex characteristics and (2) a singular evolution
// matrix, the two cases Martinson & Barton's items 2-3 turn into verdicts.  MEASURING driver: it asserts only
// that setup returns SOMETHING and prints the classification, so the next revision can be scoped from evidence.
#include <iostream>
#include <iomanip>
#include <vector>
#ifndef OCFE_OCFESLV_HEADER
#define OCFE_OCFESLV_HEADER "ocfeslv.hpp"
#endif
#include OCFE_OCFESLV_HEADER
using namespace mc;

static int g_pass = 0, g_fail = 0;
static void check( char const* nm, bool ok )
{ std::cout << "  " << std::left << std::setw(58) << nm << std::right << ( ok? " PASS": " FAIL" ) << std::endl;
  ( ok? g_pass: g_fail )++; }

static void run( char const* tag, int kind )
{
  std::cout << "\n---- " << tag << " ----\n";
  FFGraph D; OCFESLV oc( &D );
  FFVar t = D.add_var("t"), x = D.add_var("x");
  FFVar u1 = D.add_var("u1(t,x)"), u2 = D.add_var("u2(t,x)");
  FFPartial OpP;
  int const T_INT = FFDom::ALL - FFDom::LB, X_INT = FFDom::ALL - FFDom::LB - FFDom::UB;
  oc.add_domain( t, FFDom( 0., 0.2, 2, FFDom::LGR, 4 ) );
  oc.add_domain( x, FFDom( 0., 1.0, 2, FFDom::LGL, 5 ) );
  oc.add_state( u1, {t,x} ); oc.add_state( u2, {t,x} );
  oc.set_evolution_domain( t );
  oc.update_ref( u1, 0. ); oc.update_ref( u2, 0. );
  OCFESLV::EqnOptions io( OCFESLV::EqnRole::INTERIOR, 0 ), ii( OCFESLV::EqnRole::INITIAL, 0 ),
                      ib( OCFESLV::EqnRole::BOUNDARY, 0 );
  if( kind == 0 ){            // COMPLEX characteristics: u1_t + u2_x = 0, u2_t - u1_x = 0  (eigenvalues +-i)
    oc.add_equation( OpP(u1,t) + OpP(u2,x), {t,x}, {T_INT, X_INT}, io );
    oc.add_equation( OpP(u2,t) - OpP(u1,x), {t,x}, {T_INT, X_INT}, io );
  }
  else{                       // SINGULAR evolution matrix (descriptor): 0*u1_t + u2_x = 0, u1_t + u1_x = 0
    oc.add_equation( OpP(u2,x) - 0.,             {t,x}, {T_INT, X_INT}, io );
    oc.add_equation( OpP(u1,t) + OpP(u1,x),      {t,x}, {T_INT, X_INT}, io );
  }
  oc.add_equation( u1 - 0., {t,x}, {FFDom::LB,FFDom::ALL}, ii );
  oc.add_equation( u2 - 0., {t,x}, {FFDom::LB,FFDom::ALL}, ii );
  oc.add_equation( u1 - 0., {t,x}, {T_INT,FFDom::LB},      ib );
  oc.add_equation( u2 - 0., {t,x}, {T_INT,FFDom::UB},      ib );
  oc.options.REDUCE.ORDER      = OCFESLV::Options::RED_FULL;
  oc.options.DISPLAY_LEVEL     = 1;
  oc.options.FATAL.REDUCED_DOF = false;
  bool ok = false;
  try { ok = oc.setup(); } catch( ... ) { std::cout << "  setup THREW\n"; }
  std::cout << "  setup=" << ok << "  balance=" << oc.dof_balance().str << "\n";
  for( auto const& [bid, cls] : oc.block_classification() )
    std::cout << "  block " << bid << ": " << OCFESLV::pde_type_name( cls.type )
              << "  hyperbolic=" << cls.evolution_hyperbolic
              << "  At_singular=" << cls.At_singular
              << "  descriptor="  << cls.descriptor
              << "  degenerate="  << cls.degenerate << "\n";
  for( auto const& fc : oc.face_conditions() )
    std::cout << "  face " << fc.direction << ( fc.face == FFDom::LB? " LB": " UB" )
              << ": leaving=" << fc.outgoing << " entering=" << fc.incoming
              << " rows_here=" << fc.rows_at_face << "\n";

  std::cout << "  well-posedness: " << oc.wellposedness().size() << " finding(s)\n";
  bool complex_found = false, singular_found = false, imbalance_found = false;
  for( auto const& f : oc.wellposedness() ){
    std::cout << "    - " << f.detail.substr( 0, 96 ) << "\n";
    if( f.kind == OCFESLV::t_Finding::COMPLEX_CHARACTERISTICS ) complex_found  = true;
    if( f.kind == OCFESLV::t_Finding::SINGULAR_INDEX )          singular_found = true;
    if( f.kind == OCFESLV::t_Finding::DOF_IMBALANCE )           imbalance_found = true;
  }
  if( kind == 0 ){
    check( "complex characteristics are reported as a FINDING", complex_found );
    check( "and the deficit in conditions is reported too",     imbalance_found );
  }
  else{
    check( "the descriptor block is reported as SINGULAR in t", singular_found );
    check( "and the deficit in conditions is reported too",     imbalance_found );
  }
}

int main()
{
  std::cout << "OCFESLV ** header: " << OCFESLV::HEADER_ID << "\n";
  run( "complex characteristics (eigenvalues +-i): the paper's item 2", 0 );
  run( "singular evolution matrix (descriptor): the paper's item 3",    1 );
  std::cout << "\n  OCFE_pencil: " << g_pass << " passed, " << g_fail << " failed\n";
  return g_fail? 1: 0;
}
