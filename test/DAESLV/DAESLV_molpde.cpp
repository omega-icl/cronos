// DAESLV_molpde.cpp -- a model DISTRIBUTED IN SPACE declared on DAESLV itself (MoL workplan 7'b): DAESLV::setup
// semi-discretises it (OCFEMOL) and makes the semi-discrete system its own WORKING model, the declared model staying
// the user's PDE; solve, results and statistics are the solver's own.
// Gates: (A) the internal path reproduces the explicit one -- OCFEMOL + export_semidiscrete + DAESLV, on the same
// declarations -- BIT FOR BIT at every stage time, through solve(P) and through the by-name solve; (B) ODESLV refuses
// such a model, naming DAESLV; (C) what 7'b does not take yet -- inputs, outputs -- is refused by name.
#include <iostream>
#include <iomanip>
#include <cmath>
#include "daeslv_idas.hpp"
#include "odeslv_cvodes.hpp"
#include "ocfeslv.hpp"
using namespace mc;
static int npass = 0, nfail = 0;
static void check( bool ok, std::string const& what ){ std::cout << "  " << std::left << std::setw(70) << what << ( ok? "PASS": "FAIL" ) << "\n"; ( ok? npass: nfail )++; }

static constexpr double ALPHA = 0.1, QFLUX = 1.0, TEND = 1.0;
// the heated rod of DAESLV_molgate, declared on any FFModel-based solver
template <typename OC>
static void declare_rod( OC& oc, FFGraph& G, FFVar& t, FFVar& z, FFVar& u, size_t NEL, size_t NNOD )
{
  t = G.add_var("t"); z = G.add_var("z"); u = G.add_var("u");
  oc.add_domain( t, FFDom( 0., TEND, 2, FFDom::LGR, 3 ) );            // for DAESLV: the stages
  oc.add_domain( z, FFDom( 0., 1., NEL, FFDom::LGL, NNOD ) );
  oc.set_evolution_domain( t );
  oc.add_state( u, {t,z} ); oc.update_ref( u, 0. );
  FFPartial OpP;
  int const TI = FFDom::ALL - FFDom::LB, ZI = FFDom::ALL - FFDom::LB - FFDom::UB;
  typename OC::EqnOptions io( OC::EqnRole::INTERIOR, 0 ), ib( OC::EqnRole::BOUNDARY, 0 ), ii( OC::EqnRole::INITIAL, 0 );
  oc.add_equation( OpP(u,t) - OpP( ALPHA*OpP(u,z), z ), {t,z}, {TI, ZI},                io );
  oc.add_equation( u,                                   {t,z}, {TI, FFDom::LB},         ib );
  oc.add_equation( ALPHA*OpP(u,z) - QFLUX,              {t,z}, {TI, FFDom::UB},         ib );
  oc.add_equation( u,                                   {t,z}, {FFDom::LB, FFDom::ALL}, ii );
}
// the same rod with its boundary flux and initial profile either INPUTS (q lumped, u_ic(z) distributed) or written in
static double f_ic( double z ){ return 0.1*z*( 1. - z ); }
template <typename OC>
static void declare_rod_io( OC& oc, FFGraph& G, FFVar& t, FFVar& z, FFVar& u, FFVar* q, FFVar* uic,
                            int uic_type = -1, size_t uic_nnode = 0 )
{
  t = G.add_var("t"); z = G.add_var("z"); u = G.add_var("u");
  oc.add_domain( t, FFDom( 0., TEND, 2, FFDom::LGR, 3 ) );
  oc.add_domain( z, FFDom( 0., 1., 3, FFDom::LGL, 5 ) );
  oc.set_evolution_domain( t );
  oc.add_state( u, {t,z} ); oc.update_ref( u, 0. );
  if( q ){ *q = G.add_var("q");  oc.add_input( *q, {} ); }
  if( uic ){
    *uic = G.add_var("u_ic");
    if( uic_type < 0 ) oc.add_input( *uic, {z} );
    else               oc.add_input( *uic, {z}, static_cast<FFDom::TYPE>( uic_type ), uic_nnode );   // its own collocation
  }
  FFPartial OpP;
  int const TI = FFDom::ALL - FFDom::LB, ZI = FFDom::ALL - FFDom::LB - FFDom::UB;
  typename OC::EqnOptions io( OC::EqnRole::INTERIOR, 0 ), ib( OC::EqnRole::BOUNDARY, 0 ), ii( OC::EqnRole::INITIAL, 0 );
  oc.add_equation( OpP(u,t) - OpP( ALPHA*OpP(u,z), z ),       {t,z}, {TI, ZI},                io );
  oc.add_equation( u,                                         {t,z}, {TI, FFDom::LB},         ib );
  oc.add_equation( ALPHA*OpP(u,z) - ( q? *q: FFVar( QFLUX ) ), {t,z}, {TI, FFDom::UB},         ib );
  oc.add_equation( u - ( uic? *uic: 0.1*z*( 1. - z ) ),        {t,z}, {FFDom::LB, FFDom::ALL}, ii );
}
// the rod with two outputs -- its inventory at the end, its outlet value integrated over time -- on a given t grid
template <typename OC>
static void declare_rod_out( OC& oc, FFGraph& G, FFVar& t, FFVar& z, FFVar& u, size_t NELT, size_t NODT )
{
  t = G.add_var("t"); z = G.add_var("z"); u = G.add_var("u");
  oc.add_domain( t, FFDom( 0., TEND, NELT, FFDom::LGR, NODT ) );
  oc.add_domain( z, FFDom( 0., 1., 3, FFDom::LGL, 5 ) );
  oc.set_evolution_domain( t );
  oc.add_state( u, {t,z} ); oc.update_ref( u, 0. );
  FFPartial OpP;  FFIntegral OpI;
  int const TI = FFDom::ALL - FFDom::LB, ZI = FFDom::ALL - FFDom::LB - FFDom::UB;
  typename OC::EqnOptions io( OC::EqnRole::INTERIOR, 0 ), ib( OC::EqnRole::BOUNDARY, 0 ), ii( OC::EqnRole::INITIAL, 0 );
  oc.add_equation( OpP(u,t) - OpP( ALPHA*OpP(u,z), z ), {t,z}, {TI, ZI},                io );
  oc.add_equation( u,                                   {t,z}, {TI, FFDom::LB},         ib );
  oc.add_equation( ALPHA*OpP(u,z) - QFLUX,              {t,z}, {TI, FFDom::UB},         ib );
  // initial data COMPATIBLE with both boundary conditions (u = 0 at z = 0, alpha u_z = 1 at z = 1) and not a steady
  // state: the solution stays smooth, so the two time discretisations converge fast.  (With 0.1 z (1 - z), whose slope
  // at the outlet is -0.1 against the flux's 10, the outlet has a sqrt(t) layer that a fixed time grid does not
  // resolve -- measured: the outlet integral 8.7e-4 apart.)
  oc.add_equation( u - ( 10.*z + z*( 1. - z )*( 1. - z ) ), {t,z}, {FFDom::LB, FFDom::ALL}, ii );
  oc.add_output( OpI( u, z ), {t}, {TEND} );                         // 0: the inventory at the end
  oc.add_output( OpI( u, t ), {z}, {1.} );                           // 1: the outlet value, integrated over time
}
static double maxdiff( std::vector< std::vector<double> > const& a, std::vector< std::vector<double> > const& b )
{
  if( a.size() != b.size() || a.empty() ) return 1e30;
  double d = 0.;
  for( size_t k = 0; k < a.size(); ++k ){
    if( a[k].size() != b[k].size() ) return 1e30;
    for( size_t i = 0; i < a[k].size(); ++i ) d = std::max( d, std::fabs( a[k][i] - b[k][i] ) );
  }
  return d;
}

int main()
{
  double const RTOL = 1e-10, ATOL = 1e-12;
  for( auto const& [NEL, NNOD] : std::vector< std::pair<size_t,size_t> >{ {3,5}, {6,5} } ){
    std::cout << "ROD, " << NEL << " x " << NNOD << " LGL in z\n";
    // the explicit path: OCFEMOL, exported to a DAESLV by hand (phase 5)
    std::vector< std::vector<double> > ref;  size_t nref = 0;
    {
      FFGraph G;  OCFEMOL mol( &G );  FFVar t, z, u;
      declare_rod( mol, G, t, z, u, NEL, NNOD );
      mol.options.DISPLAY_LEVEL = 0;
      FFGraph H;  DAESLV_IDAS S;  std::vector<FFVar> xc;  FFVar tt;  std::string err;
      bool ok = mol.setup() && export_semidiscrete( mol, H, S, FFDom( 0., TEND, 2, FFDom::LGR, 3 ), xc, tt, &err );
      if( ok ){
        std::vector<double> x0( mol.semi_discrete().n_states );
        mol.init_seed( 0., x0.data() );
        for( size_t i = 0; i < x0.size(); ++i ) S.update_ref( xc[i], x0[i] );
        S.options.SOLVE.RTOL = RTOL;  S.options.SOLVE.ATOL = ATOL;  S.options.SOLVE.VERBOSE = 0;  S.options.DISPLAY_LEVEL = 0;
        ok = S.setup() && S.solve( std::vector<double>() ) == DAESLV_IDAS::NORMAL;
      }
      if( ok ){ ref = S.val_state();  nref = xc.size(); }
      check( ok, "explicit path (OCFEMOL + export + DAESLV) solves" );
    }
    // the internal path: the same declarations on DAESLV itself
    {
      FFGraph G;  DAESLV_IDAS S( &G );  FFVar t, z, u;
      declare_rod( S, G, t, z, u, NEL, NNOD );
      S.options.SOLVE.RTOL = RTOL;  S.options.SOLVE.ATOL = ATOL;  S.options.SOLVE.VERBOSE = 0;  S.options.DISPLAY_LEVEL = 0;
      bool const st = S.setup();
      check( st && S.semidiscrete_model(), "DAESLV::setup takes the model distributed in space" + std::string( st? "": " -- " + S.extract_error() ) );
      if( !st ) continue;
      check( S.semidiscrete_states().size() == nref && S.distributed_in_space(),
             "the same semi-discrete system (" + std::to_string( S.semidiscrete_states().size() ) + " unknowns); declared: the PDE" );
      bool const sv = S.solve( std::vector<double>() ) == DAESLV_IDAS::NORMAL;
      double const d1 = sv? maxdiff( S.val_state(), ref ): 1e30;
      std::cout << std::scientific << std::setprecision(2) << "    " << ref.size() << " stage times, |internal - explicit| = " << d1 << "\n";
      check( sv && d1 == 0., "(A) solve(P): the explicit path's states, bit for bit, every stage" );
      check( sv && S.stats_solve.numSteps > 0, "the solve and its statistics are the solver's own (" + std::to_string( S.stats_solve.numSteps ) + " steps)" );
      bool const sn = S.solve( std::vector<FFModel::InputVal>() ) == DAESLV_IDAS::NORMAL;
      double const d2 = sn? maxdiff( S.val_state(), ref ): 1e30;
      check( sn && d2 == 0., "(A) the by-name solve({}): the same, bit for bit" );
      bool const s2 = S.setup() && S.solve( std::vector<double>() ) == DAESLV_IDAS::NORMAL;
      check( s2 && maxdiff( S.val_state(), ref ) == 0., "a second setup() and solve: the same, bit for bit" );
    }
  }
  // (D) options.MOL reaches the semi-discretisation: NEUMANN_DIFFERENTIAL set on DAESLV = set on the OCFEMOL by hand,
  // and different from the default -- the option took effect
  {
    auto run = [&]( bool on_daeslv, OCFEMOL::Options::NeumannNode nn, std::vector< std::vector<double> >& xs ){
      FFGraph G;  FFVar t, z, u;
      if( on_daeslv ){
        DAESLV_IDAS S( &G );  declare_rod( S, G, t, z, u, 3, 5 );
        S.options.SOLVE.RTOL = RTOL;  S.options.SOLVE.ATOL = ATOL;  S.options.SOLVE.VERBOSE = 0;  S.options.DISPLAY_LEVEL = 0;
        S.options.MOL.IMPOSITION.NEUMANN_NODE = nn;
        bool const ok = S.setup() && S.solve( std::vector<double>() ) == DAESLV_IDAS::NORMAL;
        if( ok ) xs = S.val_state();
        return ok;
      }
      OCFEMOL mol( &G );  declare_rod( mol, G, t, z, u, 3, 5 );
      mol.options.DISPLAY_LEVEL = 0;  mol.options.IMPOSITION.NEUMANN_NODE = nn;
      FFGraph H;  DAESLV_IDAS S;  std::vector<FFVar> xc;  FFVar tt;  std::string err;
      if( !mol.setup() || !export_semidiscrete( mol, H, S, FFDom( 0., TEND, 2, FFDom::LGR, 3 ), xc, tt, &err ) ) return false;
      std::vector<double> x0( mol.semi_discrete().n_states );
      mol.init_seed( 0., x0.data() );
      for( size_t i = 0; i < x0.size(); ++i ) S.update_ref( xc[i], x0[i] );
      S.options.SOLVE.RTOL = RTOL;  S.options.SOLVE.ATOL = ATOL;  S.options.SOLVE.VERBOSE = 0;  S.options.DISPLAY_LEVEL = 0;
      bool const ok = S.setup() && S.solve( std::vector<double>() ) == DAESLV_IDAS::NORMAL;
      if( ok ) xs = S.val_state();
      return ok; };
    std::vector< std::vector<double> > dD, dE, dA;
    bool const ok = run( true, OCFEMOL::Options::NEUMANN_DIFFERENTIAL, dD ) && run( false, OCFEMOL::Options::NEUMANN_DIFFERENTIAL, dE )
                 && run( true, OCFEMOL::Options::NEUMANN_ALGEBRAIC, dA );
    double const same = ok? maxdiff( dD, dE ): 1e30, apart = ok? maxdiff( dD, dA ): 0.;
    std::cout << std::scientific << std::setprecision(2) << "    NEUMANN_DIFFERENTIAL: |DAESLV - by hand| = " << same
              << ", |differential - algebraic| = " << apart << "\n";
    check( ok && same == 0., "(D) options.MOL: as set on the OCFEMOL by hand, bit for bit" );
    check( ok && apart > 1e-8, "(D) options.MOL: and it takes effect (differs from the default)" );
  }
  // (B) ODESLV refuses a model distributed in space, naming DAESLV
  {
    FFGraph G;  ODESLV_CVODES O( &G );  FFVar t, z, u;
    declare_rod( O, G, t, z, u, 3, 5 );
    O.options.DISPLAY_LEVEL = 0;
    bool const st = O.setup();
    check( !st && O.extract_error().find( "distributed in space" ) != std::string::npos
               && O.extract_error().find( "DAESLV" ) != std::string::npos, "(B) ODESLV refuses it, naming DAESLV" );
  }
  // (E) inputs: the flux and the initial profile given BY NAME = the same rod with them written in
  {
    std::vector< std::vector<double> > ref;
    {
      FFGraph G;  DAESLV_IDAS S( &G );  FFVar t, z, u;
      declare_rod_io( S, G, t, z, u, nullptr, nullptr );
      S.options.SOLVE.RTOL = RTOL;  S.options.SOLVE.ATOL = ATOL;  S.options.SOLVE.VERBOSE = 0;  S.options.DISPLAY_LEVEL = 0;
      if( S.setup() && S.solve( std::vector<double>() ) == DAESLV_IDAS::NORMAL ) ref = S.val_state();
    }
    FFGraph G;  DAESLV_IDAS S( &G );  FFVar t, z, u, q, uic;
    declare_rod_io( S, G, t, z, u, &q, &uic );
    S.options.SOLVE.RTOL = RTOL;  S.options.SOLVE.ATOL = ATOL;  S.options.SOLVE.VERBOSE = 0;  S.options.DISPLAY_LEVEL = 0;
    bool const st = S.setup();
    check( st && S.input_dofs( uic ).size() == 15, "(E) inputs: q (lumped), u_ic(z) (15 DOFs: 3 elements x 5 LGL)" + std::string( st? "": " -- " + S.extract_error() ) );
    if( st ){
      bool const sv = S.solve({ { q, { QFLUX } }, { uic, FFModel::t_Fun( [&]( FFModel::t_Coord const& c ){ return f_ic( c.at( z ) ); } ) } }) == DAESLV_IDAS::NORMAL;
      double const d = sv? maxdiff( S.val_state(), ref ): 1e30;
      std::cout << std::scientific << std::setprecision(2) << "    by name vs written in: " << d << "  per stage time:";
      if( sv ) for( size_t k = 0; k < ref.size() && k < S.val_state().size(); ++k ){
        double dk = 0.;  for( size_t i = 0; i < ref[k].size(); ++i ) dk = std::max( dk, std::fabs( S.val_state()[k][i] - ref[k][i] ) );
        std::cout << " " << dk; }
      std::cout << "\n";
      // the initial state -- where the profile enters -- to round-off (the coordinates are computed two ways); the
      // trajectory within the integration tolerance (bit-level differences at t0 steer IDAS's steps: measured 3e-11 on
      // states of order 10, relative 3e-12 against RTOL = 1e-10).  A wrong flux or profile is off by 1e-2.
      double d0 = 1e30;
      if( sv && !ref.empty() && !S.val_state().empty() ){
        d0 = 0.;  for( size_t i = 0; i < ref[0].size(); ++i ) d0 = std::max( d0, std::fabs( S.val_state()[0][i] - ref[0][i] ) ); }
      check( sv && d0 <= 1e-15, "(E) by name = written in: the initial state, to round-off" );
      check( sv && d <= 10.*RTOL, "(E) by name = written in: the trajectory, within the tolerance" );
      auto refused = [&]( std::vector<FFModel::InputVal> const& v, std::string const& must ){
        try{ S.solve( v ); }catch( std::invalid_argument const& e ){ return std::string( e.what() ).find( must ) != std::string::npos; }
        return false; };
      check( refused( { { q, { QFLUX } } }, "missing: u_ic" ),                    "(E) a missing input is refused, named" );
      check( refused( { { q, { QFLUX } }, { uic, { 1., 2. } } }, "u_ic has 15 DOFs" ), "(E) a wrong count is refused, named" );
    }
  }
  // (F) an input at its OWN collocation, read by the INITIAL row: interpolated to the field's nodes -- an extrapolation
  // to the element's ends from Gauss nodes -- as OCFESLV's collocation arithmetic evaluates it.  The profile (quadratic)
  // is exact on either node set, so the initial state equals the rod with it written in, to round-off.
  {
    std::vector< std::vector<double> > ref;
    {
      FFGraph G;  DAESLV_IDAS S( &G );  FFVar t, z, u;
      declare_rod_io( S, G, t, z, u, nullptr, nullptr );
      S.options.SOLVE.RTOL = RTOL;  S.options.SOLVE.ATOL = ATOL;  S.options.SOLVE.VERBOSE = 0;  S.options.DISPLAY_LEVEL = 0;
      if( S.setup() && S.solve( std::vector<double>() ) == DAESLV_IDAS::NORMAL ) ref = S.val_state();
    }
    for( auto const& [ty, nn, what] : std::vector< std::tuple<int,size_t,std::string> >{
           { (int)FFDom::LGL, 3, "3 LGL (the field 5 LGL)" },
           { (int)FFDom::LG,  5, "5 LG, no end nodes (the field 5 LGL)" } } ){
      FFGraph G;  DAESLV_IDAS S( &G );  FFVar t, z, u, q, uic;
      declare_rod_io( S, G, t, z, u, &q, &uic, ty, nn );
      S.options.SOLVE.RTOL = RTOL;  S.options.SOLVE.ATOL = ATOL;  S.options.SOLVE.VERBOSE = 0;  S.options.DISPLAY_LEVEL = 0;
      bool const st = S.setup();
      bool const sv = st && S.solve({ { q, { QFLUX } }, { uic, FFModel::t_Fun( [&]( FFModel::t_Coord const& c ){ return f_ic( c.at( z ) ); } ) } })
                            == DAESLV_IDAS::NORMAL;
      double d0 = 1e30, d = 1e30;
      if( sv && !ref.empty() ){
        d0 = 0.;  for( size_t i = 0; i < ref[0].size(); ++i ) d0 = std::max( d0, std::fabs( S.val_state()[0][i] - ref[0][i] ) );
        d = maxdiff( S.val_state(), ref ); }
      std::cout << std::scientific << std::setprecision(2) << "    u_ic at " << what << ": " << S.input_dofs( uic ).size()
                << " DOFs; vs written in: t0 " << d0 << ", trajectory " << d << ( st? "": "  -- " + S.extract_error() ) << "\n";
      check( sv && d0 <= 1e-15 && d <= 10.*RTOL, "(F) u_ic at " + what + ": interpolated, = written in" );
    }
  }
  // (G) outputs: the spatial part of each, the time part DAESLV's own.  The WIRING: DAESLV's terminal output = OCFEMOL's own
  // evaluation of that spatial output from the state at T.  The MEANING: both outputs = OCFESLV's on the same declared
  // model at a high time order (the time discretisations differ: agreement to their error, not to round-off).
  {
    FFGraph G;  DAESLV_IDAS S( &G );  FFVar t, z, u;
    declare_rod_out( S, G, t, z, u, 2, 3 );
    S.options.SOLVE.RTOL = RTOL;  S.options.SOLVE.ATOL = ATOL;  S.options.SOLVE.VERBOSE = 0;  S.options.DISPLAY_LEVEL = 0;
    bool const ok = S.setup() && S.solve( std::vector<FFModel::InputVal>() ) == DAESLV_IDAS::NORMAL && S.val_function().size() == 2;
    check( ok, "(G) a model with outputs sets up and solves on DAESLV" + std::string( ok? "": " -- " + S.extract_error() ) );
    if( ok ){
      std::vector<double> f( S.semidiscrete_model()->n_outputs() );
      std::vector<double> const& xT = S.val_state().back();
      bool const ev = const_cast<OCFEMOL*>( S.semidiscrete_model() )->eval_outputs( TEND, xT.data(), nullptr, f.data() );
      double const dw = ev? std::fabs( f[0] - S.val_function()[0] ): 1e30;
      check( ev && dw <= 1e-15*std::max( 1., std::fabs( f[0] ) ), "(G) wiring: the terminal output = the spatial output from the state at T" );
      FFGraph G2;  OCFESLV oc( &G2 );  FFVar t2, z2, u2;
      declare_rod_out( oc, G2, t2, z2, u2, 4, 6 );
      oc.options.DISPLAY_LEVEL = 0;
      bool const oo = oc.setup() && oc.solve( std::vector<FFModel::InputVal>() ) == OCFESLV::NORMAL && oc.val_function().size() == 2;
      double e0 = 1e30, e1 = 1e30;
      if( oo ){ e0 = std::fabs( S.val_function()[0] - oc.val_function()[0] ) / std::fabs( oc.val_function()[0] );
                e1 = std::fabs( S.val_function()[1] - oc.val_function()[1] ) / std::fabs( oc.val_function()[1] ); }
      std::cout << std::scientific << std::setprecision(4) << "    inventory(T): DAESLV " << S.val_function()[0]
                << ( oo? ", OCFESLV " + std::to_string( oc.val_function()[0] ): std::string( "" ) ) << "  rel " << std::setprecision(2) << e0
                << "\n    outlet integral: DAESLV " << std::setprecision(4) << S.val_function()[1]
                << ( oo? ", OCFESLV " + std::to_string( oc.val_function()[1] ): std::string( "" ) ) << "  rel " << std::setprecision(2) << e1 << "\n";
      check( oo && e0 <= 1e-6 && e1 <= 1e-6, "(G) meaning: both outputs = OCFESLV's (time order 4 x 6), to 1e-6" );
    }
  }
  // (C) a time-varying input is taken (MoL 6c; its gate is DAESLV_moltime): declared over {t}, it is an input over t of the
  // semi-discrete DAE, its levels the integrator's parameters -- one declared DOF per level
  {
    FFGraph G;  DAESLV_IDAS S( &G );  FFVar t, z, u;
    declare_rod( S, G, t, z, u, 3, 5 );
    FFVar k = G.add_var( "k" );  S.add_input( k, {t}, FFDom::LGL, 1 );            // piecewise constant: 2 levels on 2 stages
    S.options.DISPLAY_LEVEL = 0;
    bool const ok = S.setup();
    check( ok && S.input_dofs( k ).size() == 2 && S.np() == 2,
           "(C) a time-varying input is taken: 2 levels, the integrator's 2 parameters" + ( ok? std::string(): " -- " + S.extract_error() ) );
  }
  {
    FFGraph G;  DAESLV_IDAS S( &G );  FFVar t, z, u;
    declare_rod( S, G, t, z, u, 3, 5 );
    S.add_output( u, {t,z}, std::vector<int>{ FFDom::ALL, FFDom::ALL } );   // a trajectory, distributed over t
    S.options.DISPLAY_LEVEL = 0;
    check( !S.setup() && S.extract_error().find( "distributed over the evolution direction" ) != std::string::npos,
           "(C) an output distributed over t: refused, named" );
  }
  std::cout << "\n==========  DAESLV_molpde: " << npass << " passed, " << nfail << " failed" << ( nfail? "": " -- ALL PASS" ) << "  ==========\n";
  return nfail? 1: 0;
}
