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
  // (B) ODESLV refuses a model distributed in space, naming DAESLV
  {
    FFGraph G;  ODESLV_CVODES O( &G );  FFVar t, z, u;
    declare_rod( O, G, t, z, u, 3, 5 );
    O.options.DISPLAY_LEVEL = 0;
    bool const st = O.setup();
    check( !st && O.extract_error().find( "distributed in space" ) != std::string::npos
               && O.extract_error().find( "DAESLV" ) != std::string::npos, "(B) ODESLV refuses it, naming DAESLV" );
  }
  // (C) not yet taken: inputs and outputs, refused by name
  {
    FFGraph G;  DAESLV_IDAS S( &G );  FFVar t, z, u;
    declare_rod( S, G, t, z, u, 3, 5 );
    FFVar k = G.add_var( "k" );  S.add_input( k, {} );
    S.options.DISPLAY_LEVEL = 0;
    check( !S.setup() && S.extract_error().find( "inputs are taken from MoL phase 6a" ) != std::string::npos, "(C) an input: refused, naming phase 6a" );
  }
  {
    FFGraph G;  DAESLV_IDAS S( &G );  FFVar t, z, u;
    declare_rod( S, G, t, z, u, 3, 5 );
    FFEval OpE;  S.add_output( OpE( u, t, TEND ) );
    S.options.DISPLAY_LEVEL = 0;
    check( !S.setup() && S.extract_error().find( "outputs are taken from MoL phase 6a" ) != std::string::npos, "(C) an output: refused, naming phase 6a" );
  }
  std::cout << "\n==========  DAESLV_molpde: " << npass << " passed, " << nfail << " failed" << ( nfail? "": " -- ALL PASS" ) << "  ==========\n";
  return nfail? 1: 0;
}
