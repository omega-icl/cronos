// DAESLV_molsetup -- MoL phase 3 (P3e): OCFEMOL::setup on the heated rod and BZ.  The peel, the working model without
// t, the shared spatial layout and plan, the DOF audit; the semi-discrete shape asserted.  No integration yet.
#include <iomanip>
#include <iostream>
#include <random>
#include "ffocfemol.hpp"
using namespace mc;
static int n_pass = 0, n_fail = 0;
static void check( std::string const& w, bool ok )
{ ok? ++n_pass: ++n_fail; std::cout << "    " << std::left << std::setw(66) << w << ( ok? "PASS": "FAIL" ) << "\n"; }

// The Jacobian of Rc by compressed FAD against central differences of eval_residual at a random point: the largest
// absolute discrepancy over the pattern, and the largest entry OUTSIDE the pattern (a pattern that misses a
// dependency would show there).
static void check_jacobian( std::string const& tag, OCFEMOL& oc )
{
  auto const& S = oc.semi_discrete();
  std::mt19937 rng( 11 ); std::uniform_real_distribution<double> U( 0.2, 0.8 );
  std::vector<double> xc( S.n_states ), r0( S.n_rows ), rp( S.n_rows ), rm( S.n_rows ), rfad( S.n_rows );
  for( auto& v : xc ) v = U( rng );
  std::vector<size_t> ptr, col; std::vector<double> val;
  bool const ok = oc.deriv_residual( 0.37, xc.data(), rfad.data(), ptr, col, val, nullptr );
  check( tag + ": deriv_residual runs", ok );
  if( !ok ) return;
  oc.eval_residual( 0.37, xc.data(), r0.data() );
  double maxres = 0.; for( size_t r = 0; r < S.n_rows; ++r ) maxres = std::max( maxres, std::fabs( rfad[r] - r0[r] ) );
  check( tag + ": the FAD pass's residual equals eval_residual's", maxres < 1e-13 );
  double const h = 1e-6; double maxin = 0., maxout = 0.; size_t nnz = col.size();
  for( size_t c = 0; c < S.n_states; ++c ){
    std::vector<double> xp = xc, xm = xc; xp[c] += h; xm[c] -= h;
    oc.eval_residual( 0.37, xp.data(), rp.data() ); oc.eval_residual( 0.37, xm.data(), rm.data() );
    for( size_t r = 0; r < S.n_rows; ++r ){
      double const fd = ( rp[r] - rm[r] ) / ( 2.*h );
      double jac = 0.; bool inpat = false;
      for( size_t k = ptr[r]; k < ptr[r+1]; ++k ) if( col[k] == c ){ jac = val[k]; inpat = true; break; }
      if( inpat ) maxin = std::max( maxin, std::fabs( jac - fd ) ); else maxout = std::max( maxout, std::fabs( fd ) );
    }
  }
  std::cout << "      Jacobian: " << nnz << " nonzeros over " << S.n_rows << "x" << S.n_states << "; max |FAD-FD| in pattern "
            << std::scientific << std::setprecision(2) << maxin << ", max |FD| outside pattern " << maxout << std::defaultfloat << "\n";
  check( tag + ": FAD matches central differences on the pattern (< 1e-6)", maxin < 1e-6 );
  check( tag + ": nothing outside the pattern (< 1e-9)",                  maxout < 1e-9 );
}

// The INDEPENDENT check of Rc and Mc together: at a state the discretisation represents exactly, the semi-discrete
// right-hand side must equal the PDE's.  Fields set from @p fu (the primary field) and @p fd (every other field, the
// order-reduction auxiliary); for each differential column, the implied u_t = -sum_r Rc[r] / Mc[r,c] over its single
// mass row must equal @p expect.  Self-consistency checks (FAD vs FD of the same Rc) cannot see a wrong Mc -- this can.
template <typename FU, typename FD>
static void check_rhs( std::string const& tag, OCFEMOL& oc, FFVar const& u, FU fu, FD fd, double expect )
{
  size_t const n = oc.semi_discrete().n_states, m = oc.semi_discrete().n_rows;
  std::vector<double> xc( n, 0. ), rc( m );
  for( auto const& st : oc.peel().states ){
    auto const& cols = oc.field_columns( st ); auto const& nodes = oc.field_nodes( st );
    for( size_t k = 0; k < cols.size(); ++k ){ double const z = nodes[k].empty()? 0.: nodes[k][0].second;
      xc[ cols[k] ] = ( st.id() == u.id() )? fu( z ): fd( z ); }
  }
  oc.eval_residual( 0., xc.data(), rc.data() );
  std::vector<int> nm( n, 0 ); double worst = 0.; size_t ndiff = 0;
  for( auto const& e : oc.mass() ) ++nm[e.col];
  bool single = true;  for( size_t c = 0; c < n; ++c ) if( nm[c] > 1 ) single = false;
  check( tag + ": each differential unknown has exactly one mass row", single );
  for( auto const& e : oc.mass() ){
    double const ut = -rc[e.row] / e.coeff.num().val();
    worst = std::max( worst, std::fabs( ut - expect ) ); ++ndiff;
  }
  std::cout << "      implied u_t at " << ndiff << " differential unknowns: max |u_t - " << expect << "| = "
            << std::scientific << std::setprecision(2) << worst << std::defaultfloat << "\n";
  check( tag + ": the semi-discrete RHS equals the PDE's on an exactly represented state", worst < 1e-10 );
}

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
      check_jacobian( "rod", oc );
      // u = z^2, flux auxiliary D = 0.1 u_z = 0.2 z  ->  u_t = 0.1 u_zz = 0.2 everywhere
      check_rhs( "rod", oc, u, []( double z ){ return z*z; }, []( double z ){ return 0.2*z; }, 0.2 );
      // the Ini stream as data: 11 stream + 4 dropped entries of the row  u - 0  at u's nodes; evaluated, each equals u there
      {
        auto const& IS = oc.ini_stream(); size_t nd = 0; for( auto const& e : IS ) nd += e.dropped;
        check( "rod: Ini table has 15 entries (11 stream + 4 dropped)", IS.size() == 15 && nd == 4 );
        std::vector<double> xc( S.n_states ), ini( IS.size() ); std::mt19937 rng( 5 ); std::uniform_real_distribution<double> U( 0., 1. );
        for( auto& v : xc ) v = U( rng );
        check( "rod: eval_initial runs", oc.eval_initial( 0., xc.data(), ini.data() ) );
        double dm = 0.; for( size_t k = 0; k < IS.size(); ++k ) dm = std::max( dm, std::fabs( ini[k] - xc[ IS[k].target ] ) );
        check( "rod: each INITIAL evaluation (u - 0) equals u at its node",   dm < 1e-14 );
        std::vector<size_t> ptr, col; std::vector<double> val; bool okd = oc.deriv_initial( 0., xc.data(), nullptr, ptr, col, val ) && col.size() == IS.size();
        for( size_t k = 0; okd && k < IS.size(); ++k ) okd = ( col[k] == IS[k].target && std::fabs( val[k] - 1. ) < 1e-14 );
        check( "rod: deriv_initial is the identity on the targets",           okd );
        // the export's INITIAL rows: initial_expression(k) in a separate DAG must evaluate exactly as eval_initial
        FFGraph H;  std::vector<FFVar> xv( S.n_states );  for( auto& v : xv ) v = H.add_var();  FFVar th = H.add_var();
        std::vector<FFVar> ex;  for( size_t k = 0; k < IS.size(); ++k ) ex.push_back( oc.initial_expression( k, xv, th ) );
        std::vector<FFVar> lv = xv;  lv.push_back( th );  std::vector<double> lx = xc;  lx.push_back( 0. );
        std::vector<double> ev( ex.size() );  H.eval( ex, ev, lv, lx );
        double de = 0.;  for( size_t k = 0; k < IS.size(); ++k ) de = std::max( de, std::fabs( ev[k] - ini[k] ) );
        check( "rod: initial_expression(k) evaluates exactly as eval_initial",  de == 0. );
      }
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
      check_jacobian( "BZ", oc );
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
      check_rhs( "rod/DIFFERENTIAL", oc, u, []( double z ){ return z*z; }, []( double z ){ return 0.2*z; }, 0.2 );
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
  //         (On the plan's three passes here and in H: the PDE row enters MoL's plan as its remainder -u_zz - lambda,
  //          the hidden constraint IS u_zz + lambda, so the spatial rows receiving the interface claims are dependent
  //          (k=1); the trace projector promotes two Schur-eliminated taus to explicit and rebuilds, the exact-natural
  //          restore rebuilds once more, and 3 explicit taus survive.  OCFESLV's row carries u_t and is independent, so
  //          its plan takes one pass.  A property of reduced models under MoL, not a defect -- 1 Nov 2026.)
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
      check_jacobian( "DAE", oc );
    }
  }

  // ---- K. two spatial directions: PDE5's manufactured 2-D heat equation,  u_t - k (u_xx + u_yy) = f(t,x,y),  Dirichlet
  //         on all four sides, on a 2x2-element CGL mesh with 4 nodes per direction.  Exercises the tensor column
  //         placement, face resolution in two directions, the corner handling of the plan, the Ini coverage iteration.
  {
    std::cout << "\n  K. PDE5: 2-D heat equation (t, x, y)\n";
    FFGraph G; OCFEMOL oc( &G ); FFPartial OpP;
    FFVar t = G.add_var("t"), x = G.add_var("x"), y = G.add_var("y"), U = G.add_var("U");
    double const kappa = 0.1, x0 = 0.3, y0 = 0.6, width0 = 0.05, spread = 0.02, A = 1., bg = 0.2;
    FFVar DX = x - x0, DY = y - y0, S = width0 + spread*t, RR = DX*DX + DY*DY;
    FFVar BE = A*width0*exp( -RR/S )/S, UE = bg + BE;
    FFVar FE = BE*spread*( -1./S + RR/(S*S) ) - kappa*BE*( 4.*RR/(S*S) - 4./S );
    FFVar PDE = OpP(U,t) - kappa*( OpP(U,{x,2}) + OpP(U,{y,2}) ) - FE, BC = U - UE;
    size_t const NE = 2, NN = 4;
    oc.add_domain( t, FFDom( 0., 1., 1, FFDom::LGR, 2 ) );
    oc.add_domain( x, FFDom( 0., 1., NE, FFDom::CGL, NN ) );
    oc.add_domain( y, FFDom( 0., 1., NE, FFDom::CGL, NN ) );
    oc.set_evolution_domain( t );
    oc.add_state( U, {t,x,y} ); oc.update_ref( U, bg );
    int const TI = FFDom::ALL - FFDom::LB, XI = FFDom::ALL - FFDom::LB - FFDom::UB, YI = XI;
    OCFEMOL::EqnOptions io( OCFEMOL::EqnRole::INTERIOR, 0 ), ib( OCFEMOL::EqnRole::BOUNDARY, 0 ), ii( OCFEMOL::EqnRole::INITIAL, 0 );
    oc.add_equation( PDE, {t,x,y}, {TI, XI, YI}, io );
    oc.add_equation( BC,  {t,x,y}, {FFDom::LB, FFDom::ALL, FFDom::ALL}, ii );
    oc.add_equation( BC,  {t,x,y}, {TI, FFDom::LB, YI}, ib );  oc.add_equation( BC, {t,x,y}, {TI, FFDom::UB, YI}, ib );
    oc.add_equation( BC,  {t,x,y}, {TI, FFDom::ALL, FFDom::LB}, ib );  oc.add_equation( BC, {t,x,y}, {TI, FFDom::ALL, FFDom::UB}, ib );
    oc.options.DISPLAY_LEVEL = 1;
    bool ok = false;
    try{ ok = oc.setup(); } catch( std::exception const& e ){ std::cout << "      threw: " << e.what() << "\n"; }
    check( "2-D: OCFEMOL::setup accepts", ok );
    if( !ok ) std::cout << "      status " << (int)oc.setup_status() << "  blocker: " << oc.peel().blocker << "\n";
    else {
      auto const& P = oc.peel(); auto const& S = oc.semi_discrete();
      std::cout << "      shape: " << S.n_states << " unknowns (" << S.n_tau << " tau), " << S.n_rows << " rows; differential "
                << S.n_diff_cols << "; Ini stream " << S.n_ini_stream << ", dropped " << S.n_ini_dropped << "\n";
      // fields: U and its two order-reduction auxiliaries (d_x U, d_y U); nodes per field: (2x4)^2 = 64
      check( "2-D: 3 fields (U, dU/dx, dU/dy), 1 differential",            P.states.size() == 3 && P.n_differential() == 1 );
      check( "2-D: the full model classifies PARABOLIC",                    S.model_class.size() == 1 && S.model_class.begin()->second == OCFEMOL::PARABOLIC );
      check( "2-D: square",                                                 S.n_rows == S.n_states );
      check( "2-D: 192 nodal unknowns (3 fields x 64 nodes)",               S.n_states - S.n_tau == 192 );
      // U's interior nodes that keep their PDE row: 64 nodes - the 28 on the four sides (Dirichlet) - the pivoted
      // duplicates of the interior faces; the exact number depends on the plan's pivots, so assert the invariants
      check( "2-D: the Ini stream covers exactly the differential unknowns", S.n_ini_stream == S.n_diff_cols );
      check( "2-D: differential unknowns strictly fewer than U's 64 nodes",   S.n_diff_cols > 0 && S.n_diff_cols < 64 );
      check( "2-D: structural index along t <= 1",                           P.index <= 1 );
      check_jacobian( "2-D", oc );
    }
  }

  // ---- L. the operation: Rc as a DAG external operation in a SEPARATE graph (the integrator's), evaluated in double
  //         through the DAG and differentiated by MC++'s sparse forward mode -- the path _extract_from_model takes.
  {
    std::cout << "\n  L. FFMOLRES: Rc(t, xc) as an external operation in the integrator's DAG (rod)\n";
    FFGraph G; OCFEMOL oc( &G ); FFPartial OpP;
    FFVar t = G.add_var("t"), z = G.add_var("z"), u = G.add_var("u");
    oc.add_domain( t, FFDom( 0., 1., 1, FFDom::LGR, 2 ) ); oc.add_domain( z, FFDom( 0., 1., 3, FFDom::LGL, 5 ) );
    oc.set_evolution_domain( t ); oc.add_state( u, {t,z} ); oc.update_ref( u, 0. );
    int const TI = FFDom::ALL - FFDom::LB, ZI = FFDom::ALL - FFDom::LB - FFDom::UB;
    OCFEMOL::EqnOptions io( OCFEMOL::EqnRole::INTERIOR, 0 ), ib( OCFEMOL::EqnRole::BOUNDARY, 0 ), ii( OCFEMOL::EqnRole::INITIAL, 0 );
    oc.add_equation( OpP(u,t) - OpP( 0.1*OpP(u,z), z ), {t,z}, {TI, ZI}, io );
    oc.add_equation( u, {t,z}, {TI, FFDom::LB}, ib ); oc.add_equation( 0.1*OpP(u,z) - 1., {t,z}, {TI, FFDom::UB}, ib );
    oc.add_equation( u, {t,z}, {FFDom::LB, FFDom::ALL}, ii );
    oc.options.DISPLAY_LEVEL = 0;
    check( "op: OCFEMOL::setup", oc.setup() );
    size_t const n = oc.semi_discrete().n_states, m = oc.semi_discrete().n_rows;
    FFGraph H; std::vector<FFVar> xc( n ); for( auto& v : xc ) v = H.add_var(); FFVar tt = H.add_var("t");
    FFMOLRES RES; std::vector<FFVar> rc = RES( xc, tt, &oc );
    check( "op: n_rows outputs inserted", rc.size() == m );
    std::mt19937 rng( 3 ); std::uniform_real_distribution<double> U( 0., 1. );
    std::vector<double> x( n+1 ); for( auto& v : x ) v = U( rng );
    std::vector<FFVar> leaves = xc; leaves.push_back( tt );
    std::vector<double> r1( m ), r2( m ); H.eval( rc, r1, leaves, x ); oc.eval_residual( x[n], x.data(), r2.data() );
    double d = 0.; for( size_t i = 0; i < m; ++i ) d = std::max( d, std::fabs( r1[i] - r2[i] ) );
    check( "op: double evaluation through the DAG equals eval_residual", d == 0. );
    auto J = H.SFAD( rc, leaves ); auto const& Jr = std::get<0>(J); auto const& Jc = std::get<1>(J); auto const& Je = std::get<2>(J);
    std::vector<double> jv( Je.size() ); H.eval( Je, jv, leaves, x );
    std::vector<size_t> ptr, col; std::vector<double> val, drdt( m ); oc.deriv_residual( x[n], x.data(), nullptr, ptr, col, val, drdt.data() );
    size_t matched = 0; double dj = 0.;
    for( size_t k = 0; k < Je.size(); ++k ){
      unsigned const r = Jr[k], c = Jc[k]; double ref = 0.; bool found = false;
      if( c == n ){ ref = drdt[r]; found = true; } else for( size_t j = ptr[r]; j < ptr[r+1]; ++j ) if( col[j] == c ){ ref = val[j]; found = true; break; }
      if( found ){ ++matched; dj = std::max( dj, std::fabs( jv[k] - ref ) ); }
    }
    std::cout << "      SFAD: " << Je.size() << " nonzeros, " << matched << " matched to deriv_residual, max diff " << dj << "\n";
    check( "op: MC++'s sparse forward mode through the DAG reproduces deriv_residual exactly", matched == Je.size() && dj == 0. );
    // the Ini stream as pointwise expressions in the same DAG (initial_expression -- the vectorised FFMOLINI was retired,
    // the integrator peels INITIAL rows one at a time): each evaluates like eval_initial, and its derivative on the rod
    // (u - 0) is the identity on the target
    auto const& IS = oc.ini_stream();
    std::vector<FFVar> ini;  for( size_t k = 0; k < IS.size(); ++k ) ini.push_back( oc.initial_expression( k, xc, tt ) );
    std::vector<double> i1( ini.size() ), i2( ini.size() ); H.eval( ini, i1, leaves, x ); oc.eval_initial( x[n], x.data(), i2.data() );
    double di = 0.; for( size_t k = 0; k < ini.size(); ++k ) di = std::max( di, std::fabs( i1[k] - i2[k] ) );
    check( "expr: one initial_expression per Ini entry, evaluating like eval_initial", ini.size() == IS.size() && di == 0. );
    auto JI = H.SFAD( ini, leaves ); auto const& Ir = std::get<0>(JI); auto const& Ic = std::get<1>(JI); auto const& Ie = std::get<2>(JI);
    std::vector<double> iv( Ie.size() ); H.eval( Ie, iv, leaves, x );
    bool okI = Ie.size() == IS.size(); double dI = 0.;
    for( size_t k = 0; okI && k < Ie.size(); ++k ){ okI = ( Ic[k] == IS[ Ir[k] ].target ); dI = std::max( dI, std::fabs( iv[k] - 1. ) ); }
    check( "expr: SFAD of the initial expressions is the identity on the targets (rod: u - 0)", okI && dI == 0. );
  }

  std::cout << "\n==========  DAESLV_molsetup: " << n_pass << " passed, " << n_fail << " failed -- " << ( n_fail? "FAIL": "ALL PASS" ) << "  ==========\n\n";
  return n_fail? 1: 0;
}
