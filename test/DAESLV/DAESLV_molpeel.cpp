// DAESLV_molpeel -- MoL phase 2: FFModel::peel_evolution on PDE models (the heated rod, BZ in flux form),
// with M and R printed by row class and the counts asserted; then the refusals (a mixed derivative, a
// derivative in an INITIAL row) and one acceptance that matters (a dynamic boundary condition: a mass
// entry on a BOUNDARY row).  The models are declared through OCFESLV exactly as DAESLV_mol's leg 1 does,
// set up, and read along t.  Time-only models take the same code through _extract_from_model, which the
// ODESLV and DAESLV corpora gate separately.
#include <cmath>
#include <iomanip>
#include <iostream>
#include "ocfeslv.hpp"

using namespace mc;
static int n_pass = 0, n_fail = 0;
static void check( std::string const& w, bool ok )
{ ok? ++n_pass: ++n_fail; std::cout << "    " << std::left << std::setw(62) << w << ( ok? "PASS": "FAIL" ) << "\n"; }
static char const* role_name( FFModel::EqnRole r )
{ static char const* N[] = { "AUTO","INTERIOR","INITIAL","BOUNDARY","INTERFACE","LINK","SURFACE","DIAGNOSTIC" }; return N[(int)r]; }

static void print_peel( FFModel::t_Peel const& P )
{
  std::cout << "      states: ";
  for( size_t i = 0; i < P.states.size(); ++i ) std::cout << P.states[i] << ( P.differential[i]? "(diff) ": "(alg) " );
  std::cout << "\n      index along the direction: " << P.index << "\n";
  for( auto const& r : P.rows ){
    std::cout << "      row #" << r.eqn << " " << std::setw(9) << std::left << role_name( r.role ) << " M = {";
    for( auto const& [i, co] : r.mass ) std::cout << " d_t " << P.states[i] << " x " << co << ";";
    std::cout << " }  R = " << r.rest << "\n";
  }
  std::cout << "      INITIAL rows: " << P.initial.size() << "   LINK/INTERFACE rows: " << P.link.size() << "\n";
}

int main()
{
  std::cout << "\n--- DAESLV_molpeel: the PDE-level peel along t --------------------------------------\n";

  // ---- A. the heated rod:  u_t = a u_zz,  u(t,0)=0,  a u_z(t,1)=q,  u(0,z)=0
  {
    std::cout << "\n  A. heated rod\n";
    FFGraph G; OCFESLV oc( &G ); FFPartial OpP;
    FFVar t = G.add_var("t"), z = G.add_var("z"), u = G.add_var("u");
    oc.add_domain( t, FFDom( 0., 1., 2, FFDom::LGR, 3 ) );
    oc.add_domain( z, FFDom( 0., 1., 3, FFDom::LGL, 5 ) );
    oc.set_evolution_domain( t );
    oc.add_state( u, {t,z} ); oc.update_ref( u, 0. );
    int const TI = FFDom::ALL - FFDom::LB, ZI = FFDom::ALL - FFDom::LB - FFDom::UB;
    OCFESLV::EqnOptions io( OCFESLV::EqnRole::INTERIOR, 0 ), ib( OCFESLV::EqnRole::BOUNDARY, 0 ), ii( OCFESLV::EqnRole::INITIAL, 0 );
    oc.add_equation( OpP(u,t) - OpP( 0.1*OpP(u,z), z ), {t,z}, {TI, ZI},             io );
    oc.add_equation( u,                                 {t,z}, {TI, FFDom::LB},      ib );
    oc.add_equation( 0.1*OpP(u,z) - 1.,                 {t,z}, {TI, FFDom::UB},      ib );
    oc.add_equation( u,                                 {t,z}, {FFDom::LB, FFDom::ALL}, ii );
    oc.options.DISPLAY_LEVEL = 0;
    check( "rod: OCFESLV::setup", oc.setup() );

    FFModel::t_Peel P;
    bool const ok = oc.peel_evolution( t, P );
    check( "rod: peel_evolution accepts the model", ok );
    if( !ok ) std::cout << "      blocker: " << P.blocker << "\n";
    print_peel( P );
    check( "rod: 2 states after order reduction (u and its auxiliary)",   P.states.size() == 2 );
    check( "rod: exactly 1 differential state",                            P.n_differential() == 1 );
    check( "rod: 3 peeled rows (1 INTERIOR + 2 BOUNDARY)",                 P.rows.size() == 3 );
    check( "rod: 1 INITIAL row in the Ini stream",                         P.initial.size() == 1 );
    check( "rod: 1 LINK row (the order-reduction auxiliary)",              P.link.size() == 1 );
    size_t n_mass = 0; for( auto const& r : P.rows ) n_mass += r.mass.size();
    check( "rod: one mass entry in total, on the INTERIOR row",            n_mass == 1 && P.rows[0].mass.size() == 1 );
    check( "rod: the BOUNDARY rows are algebraic (no mass entry)",         P.rows[1].mass.empty() && P.rows[2].mass.empty() );
    check( "rod: structural index along t <= 1",                           P.index <= 1 );
  }

  // ---- B. BZ in flux form, two fields, no-flux at both ends
  {
    std::cout << "\n  B. Belousov-Zhabotinskii, flux form\n";
    FFGraph G; OCFESLV oc( &G ); FFPartial OpP;
    FFVar t = G.add_var("t"), z = G.add_var("z"), u1 = G.add_var("u1"), u2 = G.add_var("u2");
    oc.add_domain( t, FFDom( 0., 1., 2, FFDom::LGR, 3 ) );
    oc.add_domain( z, FFDom( -30., 20., 4, FFDom::LGL, 5 ) );
    oc.set_evolution_domain( t );
    oc.add_state( u1, {t,z} ); oc.update_ref( u1, 0.5 );
    oc.add_state( u2, {t,z} ); oc.update_ref( u2, 0.5 );
    int const TI = FFDom::ALL - FFDom::LB, ZI = FFDom::ALL - FFDom::LB - FFDom::UB;
    OCFESLV::EqnOptions io( OCFESLV::EqnRole::INTERIOR, 0 ), ib( OCFESLV::EqnRole::BOUNDARY, 0 ), ii( OCFESLV::EqnRole::INITIAL, 0 );
    oc.add_equation( OpP(u1,t) - OpP(OpP(u1,z),z) - u1*(1.-u1-10.*u2), {t,z}, {TI, ZI}, io );
    oc.add_equation( OpP(u2,t) - OpP(OpP(u2,z),z) + 1.25*u1*u2,        {t,z}, {TI, ZI}, io );
    oc.add_equation( OpP(u1,z), {t,z}, {TI, FFDom::LB}, ib );  oc.add_equation( OpP(u1,z), {t,z}, {TI, FFDom::UB}, ib );
    oc.add_equation( OpP(u2,z), {t,z}, {TI, FFDom::LB}, ib );  oc.add_equation( OpP(u2,z), {t,z}, {TI, FFDom::UB}, ib );
    oc.add_equation( u1 - 0.5*( 1. + tanh( z/2. ) ), {t,z}, {FFDom::LB, FFDom::ALL}, ii );
    oc.add_equation( u2 - 0.5*( 1. - tanh( z/2. ) ), {t,z}, {FFDom::LB, FFDom::ALL}, ii );
    oc.options.DISPLAY_LEVEL = 0;
    check( "BZ: OCFESLV::setup", oc.setup() );

    FFModel::t_Peel P;
    bool const ok = oc.peel_evolution( t, P );
    check( "BZ: peel_evolution accepts the model", ok );
    if( !ok ) std::cout << "      blocker: " << P.blocker << "\n";
    print_peel( P );
    check( "BZ: 4 states (2 fields + 2 auxiliaries)",                      P.states.size() == 4 );
    check( "BZ: 2 differential states",                                   P.n_differential() == 2 );
    check( "BZ: 6 peeled rows (2 INTERIOR + 4 BOUNDARY)",                 P.rows.size() == 6 );
    check( "BZ: 2 INITIAL rows",                                          P.initial.size() == 2 );
    check( "BZ: 2 LINK rows",                                             P.link.size() == 2 );
    size_t n_mass = 0; for( auto const& r : P.rows ) n_mass += r.mass.size();
    check( "BZ: two mass entries in total (one per INTERIOR row)",        n_mass == 2 );
    check( "BZ: structural index along t <= 1",                           P.index <= 1 );
  }

  // ---- C. refusals and the one acceptance that matters
  {
    std::cout << "\n  C. refusals\n";
    // C1: a mixed derivative d_t d_z u in the interior row
    {
      FFGraph G; OCFESLV oc( &G ); FFPartial OpP;
      FFVar t = G.add_var("t"), z = G.add_var("z"), u = G.add_var("u");
      oc.add_domain( t, FFDom( 0., 1., 2, FFDom::LGR, 3 ) ); oc.add_domain( z, FFDom( 0., 1., 2, FFDom::LGL, 4 ) );
      oc.set_evolution_domain( t ); oc.add_state( u, {t,z} ); oc.update_ref( u, 0. );
      int const TI = FFDom::ALL - FFDom::LB, ZI = FFDom::ALL - FFDom::LB - FFDom::UB;
      OCFESLV::EqnOptions io( OCFESLV::EqnRole::INTERIOR, 0 ), ib( OCFESLV::EqnRole::BOUNDARY, 0 ), ii( OCFESLV::EqnRole::INITIAL, 0 );
      oc.add_equation( OpP(u,t) - OpP( OpP(u,t), z ) - OpP( OpP(u,z), z ), {t,z}, {TI, ZI}, io );   // d_z(d_t u)
      oc.add_equation( u, {t,z}, {TI, FFDom::LB}, ib ); oc.add_equation( u - 1., {t,z}, {TI, FFDom::UB}, ib );
      oc.add_equation( u - z, {t,z}, {FFDom::LB, FFDom::ALL}, ii );
      oc.options.DISPLAY_LEVEL = 0;
      bool const su = oc.setup();
      FFModel::t_Peel P; bool const ok = su && oc.peel_evolution( t, P );
      check( "mixed d_z(d_t u): peel REFUSES (or setup already does)", !ok );
      if( su ) std::cout << "      blocker: " << P.blocker << "\n";
      // reduce_order runs before the peel and makes d_t u an AUXILIARY with a LINK row  d_t u - D = 0, so the
      // mixed derivative reaches the peel as a LINK row carrying d_t -- which the LINK rule refuses.  That is the
      // workplan's own wording ("an M entry on an auxiliary that a LINK row also defines").  The direct
      // mixed-derivative test (_mixed_derivative_in) guards a model read BEFORE reduction.
      check( "mixed derivative: refused as a d_t in a LINK row (post-reduction form)",
             !su || ( P.blocker.find( "LINK row" ) != std::string::npos && P.blocker.find( "derivative along t" ) != std::string::npos ) );
    }
    // C2: a derivative along t in an INITIAL row
    {
      FFGraph G; OCFESLV oc( &G ); FFPartial OpP;
      FFVar t = G.add_var("t"), z = G.add_var("z"), u = G.add_var("u");
      oc.add_domain( t, FFDom( 0., 1., 2, FFDom::LGR, 3 ) ); oc.add_domain( z, FFDom( 0., 1., 2, FFDom::LGL, 4 ) );
      oc.set_evolution_domain( t ); oc.add_state( u, {t,z} ); oc.update_ref( u, 0. );
      int const TI = FFDom::ALL - FFDom::LB, ZI = FFDom::ALL - FFDom::LB - FFDom::UB;
      OCFESLV::EqnOptions io( OCFESLV::EqnRole::INTERIOR, 0 ), ib( OCFESLV::EqnRole::BOUNDARY, 0 ), ii( OCFESLV::EqnRole::INITIAL, 0 );
      oc.add_equation( OpP(u,t) - OpP( OpP(u,z), z ), {t,z}, {TI, ZI}, io );
      oc.add_equation( u, {t,z}, {TI, FFDom::LB}, ib ); oc.add_equation( u - 1., {t,z}, {TI, FFDom::UB}, ib );
      oc.add_equation( OpP(u,t) - z, {t,z}, {FFDom::LB, FFDom::ALL}, ii );          // a RATE as initial condition
      oc.options.DISPLAY_LEVEL = 0;
      bool const su = oc.setup();
      FFModel::t_Peel P; bool const ok = su && oc.peel_evolution( t, P );
      check( "d_t in an INITIAL row: peel REFUSES (or setup already does)", !ok );
      if( su ) std::cout << "      blocker: " << P.blocker << "\n";
    }
    // C3: a dynamic boundary condition  u_t = -u  at z = 1  -- accepted, with a mass entry on a BOUNDARY row
    {
      FFGraph G; OCFESLV oc( &G ); FFPartial OpP;
      FFVar t = G.add_var("t"), z = G.add_var("z"), u = G.add_var("u");
      oc.add_domain( t, FFDom( 0., 1., 2, FFDom::LGR, 3 ) ); oc.add_domain( z, FFDom( 0., 1., 2, FFDom::LGL, 4 ) );
      oc.set_evolution_domain( t ); oc.add_state( u, {t,z} ); oc.update_ref( u, 0. );
      int const TI = FFDom::ALL - FFDom::LB, ZI = FFDom::ALL - FFDom::LB - FFDom::UB;
      OCFESLV::EqnOptions io( OCFESLV::EqnRole::INTERIOR, 0 ), ib( OCFESLV::EqnRole::BOUNDARY, 0 ), ii( OCFESLV::EqnRole::INITIAL, 0 );
      oc.add_equation( OpP(u,t) - OpP( OpP(u,z), z ), {t,z}, {TI, ZI}, io );
      oc.add_equation( u, {t,z}, {TI, FFDom::LB}, ib );
      oc.add_equation( OpP(u,t) + u, {t,z}, {TI, FFDom::UB}, ib );                 // dynamic BC
      oc.add_equation( u - z, {t,z}, {FFDom::LB, FFDom::ALL}, ii );
      oc.options.DISPLAY_LEVEL = 0;
      bool const su = oc.setup();
      FFModel::t_Peel P; bool const ok = su && oc.peel_evolution( t, P );
      check( "dynamic BC: setup and peel accept", ok );
      if( su && !ok ) std::cout << "      blocker: " << P.blocker << "\n";
      if( ok ){
        print_peel( P );
        bool bnd_mass = false;
        for( auto const& r : P.rows ) if( r.role == FFModel::EqnRole::BOUNDARY && !r.mass.empty() ) bnd_mass = true;
        check( "dynamic BC: a BOUNDARY row carries a mass entry", bnd_mass );
      }
    }
  }

  std::cout << "\n==========  DAESLV_molpeel: " << n_pass << " passed, " << n_fail << " failed -- "
            << ( n_fail? "FAIL": "ALL PASS" ) << "  ==========\n\n";
  return n_fail? 1: 0;
}
