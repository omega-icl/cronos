// DAESLV_molsetup -- MoL phase 3 (P3e): OCFEMOL::setup on the heated rod and BZ.  The peel, the working model without
// t, the shared spatial layout and plan, the DOF audit; the semi-discrete shape asserted.  No integration yet.
#include <iomanip>
#include <iostream>
#include "ocfemol.hpp"
using namespace mc;
static int n_pass = 0, n_fail = 0;
static void check( std::string const& w, bool ok )
{ ok? ++n_pass: ++n_fail; std::cout << "    " << std::left << std::setw(66) << w << ( ok? "PASS": "FAIL" ) << "\n"; }

int main()
{
  std::cout.setf( std::ios::unitbuf );   // an abort inside setup must not eat the lines before it
  std::cout << "\n--- DAESLV_molsetup: OCFEMOL::setup ------------------------------------------------\n";
  // ---- A. heated rod: u_t = 0.1 u_zz, u(t,0)=0, 0.1 u_z(t,1)=1, u(0,z)=0 ; z: 3 elements x 5 nodes
  {
    std::cout << "\n  A. heated rod\n";
    FFGraph G; OCFEMOL oc( &G ); FFPartial OpP;
    FFVar t = G.add_var("t"), z = G.add_var("z"), u = G.add_var("u");
    oc.add_domain( t, FFDom( 0., 1., 1, FFDom::LGR, 2 ) );       // the span; nodes irrelevant under MoL
    oc.add_domain( z, FFDom( 0., 1., 3, FFDom::LGL, 5 ) );
    oc.set_evolution_domain( t );
    oc.add_state( u, {t,z} ); oc.update_ref( u, 0. );
    int const TI = FFDom::ALL - FFDom::LB, ZI = FFDom::ALL - FFDom::LB - FFDom::UB;
    OCFEMOL::EqnOptions io( OCFEMOL::EqnRole::INTERIOR, 0 ), ib( OCFEMOL::EqnRole::BOUNDARY, 0 ), ii( OCFEMOL::EqnRole::INITIAL, 0 );
    oc.add_equation( OpP(u,t) - OpP( 0.1*OpP(u,z), z ), {t,z}, {TI, ZI},                io );
    oc.add_equation( u,                                 {t,z}, {TI, FFDom::LB},         ib );
    oc.add_equation( 0.1*OpP(u,z) - 1.,                 {t,z}, {TI, FFDom::UB},         ib );
    oc.add_equation( u,                                 {t,z}, {FFDom::LB, FFDom::ALL}, ii );
    oc.options.DISPLAY_LEVEL = 1;
    bool const ok = oc.setup();
    check( "rod: OCFEMOL::setup", ok );
    if( !ok ){ std::cout << "      status " << (int)oc.setup_status() << "\n"; }
    else {
      auto const& P = oc.peel(); auto const& S = oc.semi_discrete();
      check( "rod: 2 fields, 1 differential",            P.states.size() == 2 && P.n_differential() == 1 );
      check( "rod: 1 INITIAL row aside",                  S.n_ini == 1 );
      check( "rod: semi-discrete system is square",       S.n_rows == S.n_states );
      // 3 elements x 5 nodes, face nodes duplicated: 15 nodes per field, 2 fields = 30 nodal unknowns (+ tau)
      check( "rod: 30 nodal unknowns (2 fields x 3x5 nodes)", S.n_states - S.n_tau == 30 );
      check( "rod: t span [0,1]",                         S.t_lo == 0. && S.t_hi == 1. );
      check( "rod: PARABOLIC along t (the base classifies the spatial operator: ELLIPTIC)", S.evolution == OCFEMOL::t_SemiDiscrete::PARABOLIC );
      std::cout << "      shape: " << S.n_states << " unknowns (" << S.n_tau << " tau), " << S.n_rows << " rows\n";
    }
  }
  // ---- B. BZ, two fields, no-flux ends
  {
    std::cout << "\n  B. BZ\n";
    FFGraph G; OCFEMOL oc( &G ); FFPartial OpP;
    FFVar t = G.add_var("t"), z = G.add_var("z"), u1 = G.add_var("u1"), u2 = G.add_var("u2");
    oc.add_domain( t, FFDom( 0., 1., 1, FFDom::LGR, 2 ) );
    oc.add_domain( z, FFDom( -30., 20., 4, FFDom::LGL, 5 ) );
    oc.set_evolution_domain( t );
    oc.add_state( u1, {t,z} ); oc.update_ref( u1, 0.5 );
    oc.add_state( u2, {t,z} ); oc.update_ref( u2, 0.5 );
    int const TI = FFDom::ALL - FFDom::LB, ZI = FFDom::ALL - FFDom::LB - FFDom::UB;
    OCFEMOL::EqnOptions io( OCFEMOL::EqnRole::INTERIOR, 0 ), ib( OCFEMOL::EqnRole::BOUNDARY, 0 ), ii( OCFEMOL::EqnRole::INITIAL, 0 );
    oc.add_equation( OpP(u1,t) - OpP(OpP(u1,z),z) - u1*(1.-u1-10.*u2), {t,z}, {TI, ZI}, io );
    oc.add_equation( OpP(u2,t) - OpP(OpP(u2,z),z) + 1.25*u1*u2,        {t,z}, {TI, ZI}, io );
    oc.add_equation( OpP(u1,z), {t,z}, {TI, FFDom::LB}, ib );  oc.add_equation( OpP(u1,z), {t,z}, {TI, FFDom::UB}, ib );
    oc.add_equation( OpP(u2,z), {t,z}, {TI, FFDom::LB}, ib );  oc.add_equation( OpP(u2,z), {t,z}, {TI, FFDom::UB}, ib );
    oc.add_equation( u1 - 0.5*( 1. + tanh( z/2. ) ), {t,z}, {FFDom::LB, FFDom::ALL}, ii );
    oc.add_equation( u2 - 0.5*( 1. - tanh( z/2. ) ), {t,z}, {FFDom::LB, FFDom::ALL}, ii );
    oc.options.DISPLAY_LEVEL = 1;
    bool const ok = oc.setup();
    check( "BZ: OCFEMOL::setup", ok );
    if( ok ){
      auto const& P = oc.peel(); auto const& S = oc.semi_discrete();
      check( "BZ: 4 fields, 2 differential",             P.states.size() == 4 && P.n_differential() == 2 );
      check( "BZ: 2 INITIAL rows aside",                  S.n_ini == 2 );
      check( "BZ: square",                                S.n_rows == S.n_states );
      check( "BZ: 80 nodal unknowns (4 fields x 4x5)",    S.n_states - S.n_tau == 80 );
      check( "BZ: PARABOLIC along t",                      S.evolution == OCFEMOL::t_SemiDiscrete::PARABOLIC );
      std::cout << "      shape: " << S.n_states << " unknowns (" << S.n_tau << " tau), " << S.n_rows << " rows\n";
    }
  }
  std::cout << "\n==========  DAESLV_molsetup: " << n_pass << " passed, " << n_fail << " failed -- " << ( n_fail? "FAIL": "ALL PASS" ) << "  ==========\n\n";
  return n_fail? 1: 0;
}
