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
      check( "rod: structural index along t <= 1 (checked on the full (t,z) model)", P.index <= 1 );
      check( "rod: semi-discrete system is square",       S.n_rows == S.n_states );
      // 3 elements x 5 nodes, face nodes duplicated: 15 nodes per field, 2 fields = 30 nodal unknowns (+ tau)
      check( "rod: 30 nodal unknowns (2 fields x 3x5 nodes)", S.n_states - S.n_tau == 30 );
      check( "rod: t span [0,1]",                         S.t_lo == 0. && S.t_hi == 1. );
      check( "rod: the full model classifies PARABOLIC (FFModel's classifier; the plan then sees ELLIPTIC)", S.model_class.size() == 1 && S.model_class.begin()->second == OCFEMOL::PARABOLIC );
      std::cout << "      Mc: " << oc.mass().size() << " entries; differential unknowns " << S.n_diff_cols << ", algebraic rows " << S.n_alg_rows
                << "; Ini stream " << S.n_ini_stream << ", dropped " << S.n_ini_dropped << "\n";
      check( "rod/ALGEBRAIC: 11 differential unknowns (15 u-nodes - Dirichlet - Neumann - 2 pivoted duplicates)", S.n_diff_cols == 11 );
      check( "rod/ALGEBRAIC: one Mc entry per differential unknown",   oc.mass().size() == S.n_diff_cols );
      check( "rod/ALGEBRAIC: Ini stream covers every differential unknown", S.n_ini_stream == S.n_diff_cols );
      check( "rod/ALGEBRAIC: the u-INITIAL row is dropped at the 4 algebraic u-nodes", S.n_ini_dropped == 4 );
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
      check( "BZ: structural index along t <= 1",            P.index <= 1 );
      check( "BZ: square",                                S.n_rows == S.n_states );
      check( "BZ: 80 nodal unknowns (4 fields x 4x5)",    S.n_states - S.n_tau == 80 );
      check( "BZ: the full model classifies PARABOLIC",         S.model_class.size() == 1 && S.model_class.begin()->second == OCFEMOL::PARABOLIC );
      std::cout << "      shape: " << S.n_states << " unknowns (" << S.n_tau << " tau), " << S.n_rows << " rows\n";
    }
  }
  // ---- C. the rod under IMPOSITION.NEUMANN_NODE = NEUMANN_DIFFERENTIAL: u keeps its d_t row at the flux node
  {
    std::cout << "\n  C. heated rod, NEUMANN_DIFFERENTIAL\n";
    FFGraph G; OCFEMOL oc( &G ); FFPartial OpP;
    FFVar t = G.add_var("t"), z = G.add_var("z"), u = G.add_var("u");
    oc.add_domain( t, FFDom( 0., 1., 1, FFDom::LGR, 2 ) );
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
    oc.options.IMPOSITION.NEUMANN_NODE = OCFEMOL::Options::NEUMANN_DIFFERENTIAL;
    bool const ok = oc.setup();
    check( "rod/DIFFERENTIAL: OCFEMOL::setup", ok );
    if( ok ){
      auto const& S = oc.semi_discrete();
      std::cout << "      shape: " << S.n_states << " unknowns, " << S.n_rows << " rows; differential " << S.n_diff_cols
                << "; Ini stream " << S.n_ini_stream << ", dropped " << S.n_ini_dropped << "\n";
      check( "rod/DIFFERENTIAL: still square",                       S.n_rows == S.n_states );
      check( "rod/DIFFERENTIAL: 12 differential unknowns (the flux node joins)", S.n_diff_cols == 12 );
      check( "rod/DIFFERENTIAL: Ini stream covers every differential unknown",  S.n_ini_stream == S.n_diff_cols );
      check( "rod/DIFFERENTIAL: 3 INITIAL evaluations dropped",             S.n_ini_dropped == 3 );
    }
  }

  // ---- D. index 1: an algebraic field pinned by its own row.  u_t = u_zz + v,  0 = v - u^2
  {
    std::cout << "\n  D. index-1 PDAE: algebraic field pinned by its own constraint\n";
    FFGraph G; OCFEMOL oc( &G ); FFPartial OpP;
    FFVar t = G.add_var("t"), z = G.add_var("z"), u = G.add_var("u"), v = G.add_var("v");
    oc.add_domain( t, FFDom( 0., 1., 1, FFDom::LGR, 2 ) ); oc.add_domain( z, FFDom( 0., 1., 2, FFDom::LGL, 4 ) );
    oc.set_evolution_domain( t );
    oc.add_state( u, {t,z} ); oc.update_ref( u, 0.5 ); oc.add_state( v, {t,z} ); oc.update_ref( v, 0.25 );
    int const TI = FFDom::ALL - FFDom::LB, ZI = FFDom::ALL - FFDom::LB - FFDom::UB;
    OCFEMOL::EqnOptions io( OCFEMOL::EqnRole::INTERIOR, 0 ), ib( OCFEMOL::EqnRole::BOUNDARY, 0 ), ii( OCFEMOL::EqnRole::INITIAL, 0 );
    oc.add_equation( OpP(u,t) - OpP(OpP(u,z),z) - v, {t,z}, {TI, ZI}, io );
    oc.add_equation( v - u*u,                        {t,z}, {TI, FFDom::ALL}, io );     // algebraic, pins v everywhere
    oc.add_equation( u - 0.5, {t,z}, {TI, FFDom::LB}, ib ); oc.add_equation( u - 0.5, {t,z}, {TI, FFDom::UB}, ib );
    oc.add_equation( u - 0.5, {t,z}, {FFDom::LB, FFDom::ALL}, ii );
    oc.options.DISPLAY_LEVEL = 1;
    bool const ok = oc.setup();
    check( "index-1: OCFEMOL::setup accepts", ok );
    if( !ok ) std::cout << "      blocker: " << oc.peel().blocker << "\n";
    else {
      auto const& P = oc.peel(); auto const& S = oc.semi_discrete();
      check( "index-1: structural index along t == 1",               P.index == 1 );
      check( "index-1: v algebraic, u differential",                 P.n_differential() == 1 && P.states.size() == 3 );
      check( "index-1: square",                                      S.n_rows == S.n_states );
      check( "index-1: no INITIAL evaluation for v (none declared)", S.n_ini_stream == S.n_diff_cols );
    }
  }
  // ---- E. index 2: a multiplier field pinned only through a constraint on the differential field.
  //         u_t = u_zz + lambda,  0 = u - g(z): lambda is unmatched; one differentiation of the constraint exposes it.
  {
    std::cout << "\n  E. index-2 PDAE: multiplier field (must be REFUSED)\n";
    FFGraph G; OCFEMOL oc( &G ); FFPartial OpP;
    FFVar t = G.add_var("t"), z = G.add_var("z"), u = G.add_var("u"), lam = G.add_var("lam");
    oc.add_domain( t, FFDom( 0., 1., 1, FFDom::LGR, 2 ) ); oc.add_domain( z, FFDom( 0., 1., 2, FFDom::LGL, 4 ) );
    oc.set_evolution_domain( t );
    oc.add_state( u, {t,z} ); oc.update_ref( u, 0.5 ); oc.add_state( lam, {t,z} ); oc.update_ref( lam, 0. );
    int const TI = FFDom::ALL - FFDom::LB, ZI = FFDom::ALL - FFDom::LB - FFDom::UB;
    OCFEMOL::EqnOptions io( OCFEMOL::EqnRole::INTERIOR, 0 ), ib( OCFEMOL::EqnRole::BOUNDARY, 0 ), ii( OCFEMOL::EqnRole::INITIAL, 0 );
    oc.add_equation( OpP(u,t) - OpP(OpP(u,z),z) - lam, {t,z}, {TI, ZI}, io );
    oc.add_equation( u - z*(1.-z),                      {t,z}, {TI, FFDom::ALL}, io );   // constraint on u: lam nowhere
    oc.add_equation( u, {t,z}, {TI, FFDom::LB}, ib ); oc.add_equation( u, {t,z}, {TI, FFDom::UB}, ib );
    oc.add_equation( u - z*(1.-z), {t,z}, {FFDom::LB, FFDom::ALL}, ii );
    oc.options.DISPLAY_LEVEL = 0;
    bool const ok = oc.setup();
    check( "index-2: OCFEMOL::setup REFUSES", !ok );
    std::cout << "      blocker: " << oc.peel().blocker << "\n";
    check( "index-2: the blocker names structural index 2 along t",
           oc.peel().blocker.find( "structural index 2" ) != std::string::npos && oc.peel().blocker.find( "phase 8" ) != std::string::npos );
    check( "index-2: setup status INCONSISTENT_MODEL", oc.setup_status() == OCFEMOL::SetupStatus::INCONSISTENT_MODEL );
  }
  // ---- F. the DOF audit's refusal: the rod with NO initial condition
  {
    std::cout << "\n  F. heated rod without an INITIAL row (must be REFUSED by the DOF audit)\n";
    FFGraph G; OCFEMOL oc( &G ); FFPartial OpP;
    FFVar t = G.add_var("t"), z = G.add_var("z"), u = G.add_var("u");
    oc.add_domain( t, FFDom( 0., 1., 1, FFDom::LGR, 2 ) ); oc.add_domain( z, FFDom( 0., 1., 3, FFDom::LGL, 5 ) );
    oc.set_evolution_domain( t ); oc.add_state( u, {t,z} ); oc.update_ref( u, 0. );
    int const TI = FFDom::ALL - FFDom::LB, ZI = FFDom::ALL - FFDom::LB - FFDom::UB;
    OCFEMOL::EqnOptions io( OCFEMOL::EqnRole::INTERIOR, 0 ), ib( OCFEMOL::EqnRole::BOUNDARY, 0 );
    oc.add_equation( OpP(u,t) - OpP( 0.1*OpP(u,z), z ), {t,z}, {TI, ZI}, io );
    oc.add_equation( u, {t,z}, {TI, FFDom::LB}, ib ); oc.add_equation( 0.1*OpP(u,z) - 1., {t,z}, {TI, FFDom::UB}, ib );
    oc.options.DISPLAY_LEVEL = 1;
    bool const ok = oc.setup();
    check( "no-IC: OCFEMOL::setup REFUSES (DOF audit: differential unknowns without an initial condition)", !ok );
    check( "no-IC: setup status INCONSISTENT_MODEL", oc.setup_status() == OCFEMOL::SetupStatus::INCONSISTENT_MODEL );
  }

  // ---- G. the index-2 PDE of E, REDUCED by hand at the PDE level (what phase 8 would produce): the constraint
  //         0 = u - g is differentiated once along t and u_t substituted from the PDE row -> the HIDDEN constraint
  //         0 = u_zz + lambda pins lambda; the original constraint survives only as the initial condition u(0,z) = g.
  {
    std::cout << "\n  G. index-2 PDE reduced by hand to index 1 (hidden constraint), then MoL\n";
    FFGraph G; OCFEMOL oc( &G ); FFPartial OpP;
    FFVar t = G.add_var("t"), z = G.add_var("z"), u = G.add_var("u"), lam = G.add_var("lam");
    oc.add_domain( t, FFDom( 0., 1., 1, FFDom::LGR, 2 ) ); oc.add_domain( z, FFDom( 0., 1., 2, FFDom::LGL, 4 ) );
    oc.set_evolution_domain( t );
    oc.add_state( u, {t,z} ); oc.update_ref( u, 0.25 ); oc.add_state( lam, {t,z} ); oc.update_ref( lam, 2. );
    int const TI = FFDom::ALL - FFDom::LB, ZI = FFDom::ALL - FFDom::LB - FFDom::UB;
    OCFEMOL::EqnOptions io( OCFEMOL::EqnRole::INTERIOR, 0 ), ib( OCFEMOL::EqnRole::BOUNDARY, 0 ), ii( OCFEMOL::EqnRole::INITIAL, 0 );
    oc.add_equation( OpP(u,t) - OpP(OpP(u,z),z) - lam, {t,z}, {TI, ZI}, io );
    oc.add_equation( OpP(OpP(u,z),z) + lam,            {t,z}, {TI, FFDom::ALL}, io );   // the hidden constraint: pins lambda
    oc.add_equation( u, {t,z}, {TI, FFDom::LB}, ib ); oc.add_equation( u, {t,z}, {TI, FFDom::UB}, ib );
    oc.add_equation( u - z*(1.-z), {t,z}, {FFDom::LB, FFDom::ALL}, ii );              // the original constraint, at t0 only
    oc.options.DISPLAY_LEVEL = 1;
    bool const ok = oc.setup();
    check( "reduced-by-hand: OCFEMOL::setup accepts", ok );
    if( !ok ) std::cout << "      blocker: " << oc.peel().blocker << "\n";
    else {
      auto const& P = oc.peel(); auto const& S = oc.semi_discrete();
      check( "reduced-by-hand: structural index along t == 1",        P.index == 1 );
      check( "reduced-by-hand: lambda algebraic, u differential",     P.n_differential() == 1 );
      check( "reduced-by-hand: square",                               S.n_rows == S.n_states );
      check( "reduced-by-hand: the Ini stream covers u's differential nodes", S.n_ini_stream == S.n_diff_cols );
    }
  }
  // ---- H. the same index-2 PDE, AUTOMATIC reduction (IMPOSITION.REDUCE_HIGH_INDEX): an experiment, reported, not asserted
  {
    std::cout << "\n  H. index-2 PDE with IMPOSITION.REDUCE_HIGH_INDEX (experiment)\n";
    FFGraph G; OCFEMOL oc( &G ); FFPartial OpP;
    FFVar t = G.add_var("t"), z = G.add_var("z"), u = G.add_var("u"), lam = G.add_var("lam");
    oc.add_domain( t, FFDom( 0., 1., 1, FFDom::LGR, 2 ) ); oc.add_domain( z, FFDom( 0., 1., 2, FFDom::LGL, 4 ) );
    oc.set_evolution_domain( t );
    oc.add_state( u, {t,z} ); oc.update_ref( u, 0.25 ); oc.add_state( lam, {t,z} ); oc.update_ref( lam, 2. );
    int const TI = FFDom::ALL - FFDom::LB, ZI = FFDom::ALL - FFDom::LB - FFDom::UB;
    OCFEMOL::EqnOptions io( OCFEMOL::EqnRole::INTERIOR, 0 ), ib( OCFEMOL::EqnRole::BOUNDARY, 0 ), ii( OCFEMOL::EqnRole::INITIAL, 0 );
    oc.add_equation( OpP(u,t) - OpP(OpP(u,z),z) - lam, {t,z}, {TI, ZI}, io );
    oc.add_equation( u - z*(1.-z),                      {t,z}, {TI, FFDom::ALL}, io );
    oc.add_equation( u, {t,z}, {TI, FFDom::LB}, ib ); oc.add_equation( u, {t,z}, {TI, FFDom::UB}, ib );
    oc.add_equation( u - z*(1.-z), {t,z}, {FFDom::LB, FFDom::ALL}, ii );
    oc.options.DISPLAY_LEVEL = 1;
    oc.options.IMPOSITION.REDUCE_HIGH_INDEX = true;
    bool ok = false;
    try{ ok = oc.setup(); } catch( std::exception const& e ){ std::cout << "      threw: " << e.what() << "\n"; }
    check( "automatic reduction: OCFEMOL::setup accepts (phase 8's reduction, on the full (t,z) model)", ok );
    if( !ok ) std::cout << "      blocker: " << oc.peel().blocker << "\n";
    else {
      auto const& P = oc.peel(); auto const& S = oc.semi_discrete();
      check( "automatic reduction: index 1 after the reduction",            P.index == 1 );
      check( "automatic reduction: one hidden level materialised at t0 (2 INITIAL rows aside)", S.n_ini == 2 );
      check( "automatic reduction: square",                                 S.n_rows == S.n_states );
      check( "automatic reduction: the same semi-discrete shape as the hand reduction (27 unknowns, 6 differential)",
             S.n_states == 27 && S.n_diff_cols == 6 );
    }
  }

  // ---- I. no spatial direction at all: an ODE.  MoL on an ODE must be the ODE: one unknown, one row, one IC.
  {
    std::cout << "\n  I. an ODE (no spatial direction): x' = -x, x(0) = 1\n";
    FFGraph G; OCFEMOL oc( &G ); FFPartial OpP;
    FFVar t = G.add_var("t"), x = G.add_var("x");
    oc.add_domain( t, FFDom( 0., 1., 1, FFDom::LGR, 2 ) );
    oc.set_evolution_domain( t );
    oc.add_state( x, {t} ); oc.update_ref( x, 1. );
    int const TI = FFDom::ALL - FFDom::LB;
    OCFEMOL::EqnOptions io( OCFEMOL::EqnRole::INTERIOR, 0 ), ii( OCFEMOL::EqnRole::INITIAL, 0 );
    oc.add_equation( OpP(x,t) + x, {t}, {TI},        io );
    oc.add_equation( x - 1.,       {t}, {FFDom::LB}, ii );
    oc.options.DISPLAY_LEVEL = 1;
    bool ok = false;
    try{ ok = oc.setup(); } catch( std::exception const& e ){ std::cout << "      threw: " << e.what() << "\n"; }
    check( "ODE: OCFEMOL::setup accepts", ok );
    if( !ok ) std::cout << "      status " << (int)oc.setup_status() << "  blocker: " << oc.peel().blocker << "\n";
    else {
      auto const& P = oc.peel(); auto const& S = oc.semi_discrete();
      std::cout << "      full-model class: " << OCFEMOL::pde_type_name( S.model_class.begin()->second ) << "\n";
      check( "ODE: the full model classifies DIFFERENTIAL_ORDINARY", S.model_class.size() == 1 && S.model_class.begin()->second == OCFEMOL::DIFFERENTIAL_ORDINARY );
      check( "ODE: 1 unknown, 1 row, 1 INITIAL row aside",    S.n_states == 1 && S.n_rows == 1 && S.n_ini == 1 );
      check( "ODE: the unknown is differential, Ini stream 1", S.n_diff_cols == 1 && S.n_ini_stream == 1 );
      check( "ODE: index 0",                                  P.index == 0 );
    }
  }
  // ---- J. a semi-explicit DAE, no spatial direction: x' = -x + y, 0 = y - x^2, x(0) = 1
  {
    std::cout << "\n  J. a DAE (no spatial direction): x' = -x + y, 0 = y - x^2\n";
    FFGraph G; OCFEMOL oc( &G ); FFPartial OpP;
    FFVar t = G.add_var("t"), x = G.add_var("x"), y = G.add_var("y");
    oc.add_domain( t, FFDom( 0., 1., 1, FFDom::LGR, 2 ) );
    oc.set_evolution_domain( t );
    oc.add_state( x, {t} ); oc.update_ref( x, 1. ); oc.add_state( y, {t} ); oc.update_ref( y, 1. );
    int const TI = FFDom::ALL - FFDom::LB;
    OCFEMOL::EqnOptions io( OCFEMOL::EqnRole::INTERIOR, 0 ), ii( OCFEMOL::EqnRole::INITIAL, 0 );
    oc.add_equation( OpP(x,t) + x - y, {t}, {TI},        io );
    oc.add_equation( y - x*x,          {t}, {TI},        io );
    oc.add_equation( x - 1.,           {t}, {FFDom::LB}, ii );
    oc.options.DISPLAY_LEVEL = 1;
    bool ok = false;
    try{ ok = oc.setup(); } catch( std::exception const& e ){ std::cout << "      threw: " << e.what() << "\n"; }
    check( "DAE: OCFEMOL::setup accepts", ok );
    if( !ok ) std::cout << "      status " << (int)oc.setup_status() << "  blocker: " << oc.peel().blocker << "\n";
    else {
      auto const& P = oc.peel(); auto const& S = oc.semi_discrete();
      std::cout << "      full-model class: " << OCFEMOL::pde_type_name( S.model_class.begin()->second ) << "\n";
      check( "DAE: the full model classifies DIFFERENTIAL_ALGEBRAIC, not an ODE", S.model_class.size() == 1 && S.model_class.begin()->second == OCFEMOL::DIFFERENTIAL_ALGEBRAIC );
      check( "DAE: 2 unknowns, 2 rows, 1 INITIAL row aside",   S.n_states == 2 && S.n_rows == 2 && S.n_ini == 1 );
      check( "DAE: x differential, y algebraic",               S.n_diff_cols == 1 && S.n_alg_rows == 1 );
      check( "DAE: index 1",                                   P.index == 1 );
      check( "DAE: Ini stream 1, nothing dropped",             S.n_ini_stream == 1 && S.n_ini_dropped == 0 );
    }
  }

  std::cout << "\n==========  DAESLV_molsetup: " << n_pass << " passed, " << n_fail << " failed -- " << ( n_fail? "FAIL": "ALL PASS" ) << "  ==========\n\n";
  return n_fail? 1: 0;
}
