// OCFE_inpcolloc (sandbox test, 2026-09-11): the case rev193a's gates missed --
//  (a) an input with its own collocation order/type on its domain's grid,
//  (b) setup() run TWICE (which re-runs _reset(): declared data must survive it),
//  (c) a deep copy that must carry the declaration,
//  (d) a domain redeclared AFTER its inputs (semantics recorded in NOTES_20260911m Sec.3).
#include <iostream>
#include <iomanip>
#include "ffunc.hpp"
#include OCFE_OCFESLV_HEADER
using namespace mc;

static void report( char const* tag, OCFESLV& oc, FFVar const& u, FFVar const& x )
{
  // node_colloc() returns one entry per collocation node (each a coordinate tuple).  The input shares the
  // domain's ELEMENTS with the state and differs only in how many nodes sit on each, so the node COUNTS show
  // the declared orders: nelem*n_node_input against nelem*n_node_state.
  auto const un = oc.node_colloc( u );
  auto const xn = oc.node_colloc( x );
  std::cout << "  " << std::left << std::setw(30) << tag << std::right
            << " input nodes=" << un.size() << "  state nodes=" << xn.size()
            << "  n_colloc_inp=" << oc.n_colloc_inp();
  if( !un.empty() && !xn.empty() )
    std::cout << "  first/last input node=" << un.front().front() << "/" << un.back().front()
              << "  last state node=" << xn.back().front();
  std::cout << std::endl;
}

int main()
{
  std::cout << "OCFESLV ** header: " << OCFESLV::HEADER_ID << "\n";
  FFGraph D; OCFESLV oc( &D );
  FFVar t; t.set( &D ); FFVar x; x.set( &D ); FFVar u; u.set( &D );
  oc.add_domain( t, FFDom( 0., 1., 4, FFDom::LGL, 3 ) );
  oc.add_state( x, { t } );
  oc.add_input( u, { t }, FFDom::LGR, 1 );          // piecewise-constant control on the same grid
  oc.add_equation( x - u, std::vector<FFVar>{ t }, std::vector<int>{} );
  oc.add_equation( x - 1., std::vector<FFVar>{ t }, std::vector<int>{ FFDom::LB }, OCFESLV::EqnOptions( OCFESLV::EqnRole::INITIAL, 0 ) );
  std::cout << "setup #1: " << oc.setup() << "\n";  report( "after setup #1", oc, u, x );
  std::cout << "setup #2: " << oc.setup() << "\n";  report( "after setup #2", oc, u, x );
  OCFESLV cp( oc );
  std::cout << "copy setup: " << cp.setup() << "\n"; report( "copy after setup", cp, u, x );
  oc.add_domain( t, FFDom( 0., 1., 2, FFDom::LGL, 5 ) );   // redeclared coarser, higher order
  std::cout << "setup #3: " << oc.setup() << "\n";  report( "after domain redeclared", oc, u, x );
  return 0;
}
