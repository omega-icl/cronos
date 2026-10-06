// OCFE_distout.cpp -- GATE for distributed outputs (WORKPLAN 1.1, 2026-10-06).  Under marching (the default mode)
// each window is a one-element problem, so an output distributed over the EVOLUTION direction used to be sized for one
// element's nodes, evaluated with the window's masks and merged last-window-wins: the last element's values only.
// Now every window emits the rows of ITS time element (the evolution mask applied with the global element index) and
// they are placed in the monolithic layout, values and sensitivities alike.  Checked on a nonlinear heated rod
// against an explicit monolithic solve and point outputs -- under IC_STRONG, where marching and monolithic are the
// same discretisation; under IC_WEAK they now agree too (WORKPLAN 1.2/1.3), which is reported.
#include <cstdio>
#include <cmath>
#include "ocfeslv.hpp"
using namespace mc;  typedef OCFESLV::EqnRole Role;  typedef OCFESLV::EqnOptions EO;
static int nfail = 0;
static void check( char const* what, bool ok ){ std::printf( "  %-66s %s\n", what, ok? "PASS": "FAIL" ); nfail += !ok; }

struct Run { std::vector<double> f, J; std::vector<size_t> row0, nrow; bool marching = false; };
// field: 0 = u on [t,z] (ALL, UB); 1 = u on [t,z] (NLB, UB); 2 = u on [t,z] (LB, ALL); 3 = profile u on [z] at t=1
static Run run( int force_mode, bool with_field, bool sens, bool strong = false ){
  FFPartial OpP;  FFEval OpE;  int const NLB = FFDom::ALL - FFDom::LB, INN = FFDom::ALL - FFDom::LB - FFDom::UB;
  FFGraph G;  FFVar t = G.add_var( "t" ), z = G.add_var( "z" ), u = G.add_var( "u" ), q = G.add_var( "q" );
  OCFESLV S( &G );
  S.add_domain( t, FFDom( 0., 1., 4, FFDom::LGR, 3 ) );  S.add_domain( z, FFDom( 0., 1., 4, FFDom::LGL, 5 ) );
  S.set_evolution_domain( t );  S.add_state( u, {t, z} );  S.add_input( q );  S.update_ref( q, 1. );
  S.add_equation( OpP( u, t ) - 0.1 * ( 1. + 0.5 * u ) * OpP( u, {{z,2}} ), {t, z}, {NLB, INN}, EO( Role::INTERIOR ) );
  S.add_equation( u, {t, z}, {NLB, (int)FFDom::LB}, EO( Role::BOUNDARY ) );
  S.add_equation( 0.1 * OpP( u, z ) - q, {t, z}, {NLB, (int)FFDom::UB}, EO( Role::BOUNDARY ) );
  S.add_equation( u, {t, z}, {(int)FFDom::LB, (int)FFDom::ALL}, EO( Role::INITIAL ) );
  for( double tk : { 0., 0.25, 0.5, 0.75 } ) S.add_output( u, {t, z}, std::vector<double>{ tk, 1. } );   // 0-3: points
  S.add_output( OpE( u, {{t,1},{z,1}}, {{t,1.},{z,0.5}} ) );                                             // 4: u(1, 0.5)
  if( with_field ){
    S.add_output( u, {t, z}, {(int)FFDom::ALL, (int)FFDom::UB} );                                         // 5
    S.add_output( u, {t, z}, {NLB, (int)FFDom::UB} );                                                     // 6
    S.add_output( u, {t, z}, {(int)FFDom::LB, (int)FFDom::ALL} );                                         // 7
  }
  S.add_output( u, {z}, {(int)FFDom::ALL}, std::vector<FFVar>{ t }, std::vector<double>{ 1. } );          // profile at t=1
  S.register_control( q );
  if( force_mode >= 0 ) S.options.SOLVE.MARCHING = ( force_mode == 1 );
  // explicit either way (the default became IC_STRONG on 2026-10-06; the weak-gap checks need IC_WEAK)
  S.options.INTERFACE.IMPOSITION = strong? OCFESLV::Options::IC_STRONG: OCFESLV::Options::IC_WEAK;
  S.options.SOLVE.RES_TOL = 1e-12;  S.options.DISPLAY_LEVEL = 0;
  Run R;  if( !S.setup() ) return R;
  R.marching = S.is_marching();
  std::vector<double> var, inp;  S.init( var, inp );
  if( sens ){ S.solve_fsens( var.data(), inp.data(), nullptr );  R.J = S.sens_jacobian();
              S.init( var, inp ); S.solve_asens( var.data(), inp.data(), nullptr );
              auto Ja = S.sens_jacobian();  R.J.insert( R.J.end(), Ja.begin(), Ja.end() ); }
  else S.solve( var.data(), inp.data(), nullptr );
  R.f = sens? S.sens_functions(): S.val_functions();
  for( size_t k = 0; k < S.var_output().size(); ++k ){ auto b = S.blk_fct( k ); R.row0.push_back( b.first ); R.nrow.push_back( b.second ); }
  return R;
}

int main(){
  std::printf( "OCFE_distout -- distributed outputs over the evolution direction (default settings)\n" );
  Run D = run( -1, true, false, true ), M = run( 0, true, false, true );
  check( "the default solve has run (a field output over t is present)", !D.f.empty() );
  check( "it MARCHES (a field output over t does not change the mode)", D.marching );
  check( "u on (ALL, UB): 12 rows = 4 elements x 3 Radau nodes", D.nrow.size() > 5 && D.nrow[5] == 12 );
  check( "u on (NLB, UB): 11 rows (t = 0 excluded)", D.nrow.size() > 6 && D.nrow[6] == 11 );
  check( "u on (LB, ALL): 20 rows (t = 0, 4 elements x 5 Lobatto nodes)", D.nrow.size() > 7 && D.nrow[7] == 20 );
  bool starts = D.nrow.size() > 5;
  for( size_t k = 0; starts && k < 4; ++k ) starts = std::fabs( D.f[ D.row0[5] + 3*k ] - D.f[k] ) < 1e-10;
  check( "rows at the element starts = point outputs at t = 0, .25, .5, .75", starts );
  bool same = D.f.size() == M.f.size();
  for( size_t i = 0; same && i < D.f.size(); ++i ) same = std::fabs( D.f[i] - M.f[i] ) < 1e-12;
  check( "IC_STRONG: every row = an explicit monolithic solve (1e-12)", same );
  { Run Dw = run( -1, true, false ), Mw = run( 0, true, false );  double g = 0.;
    for( size_t i = 0; i < Dw.f.size() && i < Mw.f.size(); ++i ) g = std::max( g, std::fabs( Dw.f[i] - Mw.f[i] ) );
    std::printf( "    IC_WEAK: max |marching - monolithic| over all rows = %.1e (the weak scheme's level)\n", g );
    check( "IC_WEAK: same row counts, rows within the weak scheme's level (1e-4)", Dw.f.size() == Mw.f.size() && g < 1e-4 ); }
  // a profile at a FIXED time does not need the horizon: it still marches, with the monolithic values
  // (under IC_WEAK -- set explicitly -- the gap is reported; under IC_STRONG they are the same
  // discretisation, so the profile must agree to round-off)
  Run P = run( -1, false, false ), PM = run( 0, false, false );
  check( "without a field output over t, the default solve MARCHES", P.marching );
  double dw = 0.;  for( size_t i = 0; i < P.f.size() && i < PM.f.size(); ++i ) dw = std::max( dw, std::fabs( P.f[i] - PM.f[i] ) );
  std::printf( "    IC_WEAK: max |marching - monolithic| over the outputs = %.1e\n", dw );
  Run Ps = run( -1, false, false, true ), PMs = run( 0, false, false, true );
  bool prof = Ps.marching && Ps.f.size() == PMs.f.size() && !Ps.f.empty();
  double ds = 0.;  for( size_t i = 0; prof && i < Ps.f.size(); ++i ) ds = std::max( ds, std::fabs( Ps.f[i] - PMs.f[i] ) );
  std::printf( "    IC_STRONG: max |marching - monolithic| = %.1e\n", ds );
  check( "IC_STRONG: a profile at t = 1, marching = monolithic (1e-9)", prof && ds < 1e-9 );
  // sensitivities of the field outputs: forward = adjoint
  Run S = run( -1, true, true, true );
  size_t const n = S.J.size() / 2;  double d = 0.;
  for( size_t i = 0; i < n; ++i ) d = std::max( d, std::fabs( S.J[i] - S.J[n + i] ) );
  std::printf( "    %zu sensitivity entries, max |forward - adjoint| = %.1e\n", n, d );
  check( "IC_STRONG: field-output sensitivities, forward = adjoint (1e-10)", n > 0 && d < 1e-10 );
  Run SM = run( 0, true, true, true );  double dm = 0.;
  for( size_t i = 0; i < S.J.size() && i < SM.J.size(); ++i ) dm = std::max( dm, std::fabs( S.J[i] - SM.J[i] ) );
  { size_t const nJ = S.J.size() / 2;   // forward entries first; one control -> one entry per output row
    for( size_t o = 0; o < S.row0.size() && o < SM.row0.size(); ++o ){ double dd = 0.;
      for( size_t r = S.row0[o]; r < S.row0[o] + S.nrow[o] && r < nJ && r < SM.J.size(); ++r ) dd = std::max( dd, std::fabs( S.J[r] - SM.J[r] ) );
      std::printf( "    output %zu (%zu rows): max |marching - monolithic| forward sensitivity = %.1e\n", o, S.nrow[o], dd ); } }
  check( "IC_STRONG: field-output sensitivities, marching = monolithic (1e-9)", S.J.size() == SM.J.size() && dm < 1e-9 );
  std::printf( "  OCFE_distout: %d failed\n", nfail );
  return nfail? 1: 0;
}
