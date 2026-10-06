// OCFE_zchain.cpp -- GATE for WORKPLAN 1.5 (2026-10-06): a chain HIGH-INDEX IN A SPATIAL DIRECTION sets up under every
// imposition.  p = g(z) at every node, p_z = q, q_z = r (OCFE_zindex's chain, coupled to a heat equation forced by r).
// Under IC_STRONG/IC_TRACE the exact continuity claims on p, q, r at the z-element boundary are redundant with the
// constraint; the plan's pass 0 RESTORES silenced claims (EXACT_NATURAL) and pass 1 flags the redundant ones for a
// DROP -- and the restore used to consume the single DROP retry, so setup failed at the exclusion it needed.  The
// restore now has its own retry.  Checked: setup, balance, and the exact solution p = g, q = g', r = g''.
#include <cstdio>
#include <cmath>
#include "ocfeslv.hpp"
using namespace mc;  typedef OCFESLV::EqnRole Role;  typedef OCFESLV::EqnOptions EO;
static int nfail = 0;
static void check( char const* what, bool ok ){ std::printf( "  %-62s %s\n", what, ok? "PASS": "FAIL" ); nfail += !ok; }
int main(){
  std::printf( "OCFE_zchain -- a z-high-index chain under the three impositions\n" );
  char const* NAME[] = { "IC_STRONG", "IC_WEAK", "IC_TRACE" };
  OCFESLV::Options::ImpositionType const IMP[] = { OCFESLV::Options::IC_STRONG, OCFESLV::Options::IC_WEAK, OCFESLV::Options::IC_TRACE };
  char buf[128];
  for( int i = 0; i < 3; ++i ){
    FFGraph D;  FFVar t = D.add_var( "t" ), z = D.add_var( "z" );
    FFVar u = D.add_var( "u" ), p = D.add_var( "p" ), q = D.add_var( "q" ), r = D.add_var( "r" );
    FFPartial OpP;  FFEval OpE;  int const T_INT = FFDom::ALL - FFDom::LB, Z_INT = FFDom::ALL - FFDom::LB - FFDom::UB;
    OCFESLV oc( &D );
    oc.add_domain( t, FFDom( 0., 0.2, 2, FFDom::LGR, 4 ) );  oc.add_domain( z, FFDom( 0., 1.0, 2, FFDom::LGL, 5 ) );
    for( auto const& s : { u, p, q, r } ) oc.add_state( s, {t,z} );
    oc.set_evolution_domain( t );
    oc.update_ref( u, 0. ); oc.update_ref( p, 1. ); oc.update_ref( q, 0.5 ); oc.update_ref( r, -0.6 );
    FFVar G = 1.0 + 0.5*z - 0.3*z*z;
    EO io( Role::INTERIOR, 0 ), ii( Role::INITIAL, 0 ), ib( Role::BOUNDARY, 0 );
    oc.add_equation( OpP(u,t) - OpP(u,{{z,2}}) - r, {t,z}, {T_INT, Z_INT}, io );
    oc.add_equation( p - G,        {t,z}, {(int)FFDom::ALL,(int)FFDom::ALL}, io );
    oc.add_equation( OpP(p,z) - q, {t,z}, {(int)FFDom::ALL,(int)FFDom::ALL}, io );
    oc.add_equation( OpP(q,z) - r, {t,z}, {(int)FFDom::ALL,(int)FFDom::ALL}, io );
    oc.add_equation( u, {t,z}, {(int)FFDom::LB,(int)FFDom::ALL}, ii );
    oc.add_equation( u, {t,z}, {T_INT,(int)FFDom::LB}, ib );  oc.add_equation( u, {t,z}, {T_INT,(int)FFDom::UB}, ib );
    // p, q, r at (t = 0.2, z = 0.7): exact 1.203, 0.08, -0.6
    oc.add_output( OpE( p, {{t,1},{z,1}}, {{t,0.2},{z,0.7}} ) );
    oc.add_output( OpE( q, {{t,1},{z,1}}, {{t,0.2},{z,0.7}} ) );
    oc.add_output( OpE( r, {{t,1},{z,1}}, {{t,0.2},{z,0.7}} ) );
    oc.options.REDUCE.ORDER = OCFESLV::Options::RED_FULL;  oc.options.DISPLAY_LEVEL = 0;  oc.options.FATAL.REDUCED_DOF = false;
    oc.options.INTERFACE.IMPOSITION = IMP[i];  oc.options.SOLVE.RES_TOL = 1e-12;
    bool const ok = oc.setup();
    std::snprintf( buf, sizeof buf, "%s: the chain sets up", NAME[i] );   check( buf, ok );
    std::snprintf( buf, sizeof buf, "%s: it is balanced", NAME[i] );      check( buf, ok && oc.dof_balance().balanced );
    double e = INFINITY;
    if( ok ){
      std::vector<double> var, inp;  oc.init( var, inp );  auto rep = oc.solve( var.data(), inp.data(), nullptr );
      auto const& f = oc.val_functions();
      if( rep.converged && f.size() >= 3 ) e = std::max( { std::fabs( f[0] - 1.203 ), std::fabs( f[1] - 0.08 ), std::fabs( f[2] + 0.6 ) } );
    }
    std::snprintf( buf, sizeof buf, "%s: p = g, q = g', r = g'' at (0.2, 0.7) (1e-8; %.1e)", NAME[i], e );   check( buf, e < 1e-8 );
  }
  std::printf( "  OCFE_zchain: %d failed\n", nfail );
  return nfail? 1: 0;
}
