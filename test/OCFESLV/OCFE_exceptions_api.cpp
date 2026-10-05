// ============================================================================
//  OCFE_homotopy.cpp -- gate for OCFESLV's declarative homotopy/continuation (phase 1)
//
//  Validates add_homotopy() + solve_homotopy() against a hand-rolled ramp on a small
//  stiff MMS problem, and checks that the STAGE GRAMMAR reproduces the three schedules
//  the MBC drivers use:
//      all stage 0      -> SIMULTANEOUS  (one s drives every parameter)
//      stages 0,1,2     -> STAIRCASE     (each parameter ramped in turn)
//      stages 0,0,1     -> PARTIAL grouping
//
//  Manufactured problem (steady, 1-D, deliberately stiff in kap so a direct solve fails):
//      u''(z) = kap*Da0*u*w - lam*S(z) ,  u(0)=0, u(1)=1
//      w      = 1 + mu*z                  (an algebraic companion, ramped by mu)
//  The three homotopy parameters mirror MBC's roles: lam turns on a source, mu deforms a
//  profile, kap ramps a rate constant over decades (geometric map -> exercises the
//  s -> value override).
//
//  CHECKS
//    1. schedules all converge and reach the SAME root (path fidelity: ||dx||_inf < 1e-9)
//    2. solve_homotopy matches a hand-rolled staircase ramp to solver tolerance
//    3. a DIRECT solve (all parameters at 1, no continuation) is contrasted -- it may fail,
//       which is the point of having continuation at all
//    4. report bookkeeping is self-consistent (solves = accepts + backtracks per stage)
//
//  NOTE: there is no separate solve_homotopy() -- solve() honours a registered schedule, exactly
//  as it already honours SOLVE_MARCHING.  Per-stage detail comes from continuation_report().
//
//  Build: -DOCFE_OCFESLV_HEADER='"ocfeslv_homotopy.hpp"'
// ============================================================================
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <string>
#include <utility>   // std::pair -- per-stage cap overrides

#ifndef OCFE_OCFESLV_HEADER
#define OCFE_OCFESLV_HEADER "ocfeslv.hpp"
#endif
#include OCFE_OCFESLV_HEADER
using namespace mc;

static int g_pass=0, g_fail=0;
static void check( char const* nm, bool ok )
{
  (ok?g_pass:g_fail)++;
  std::cout<<"  "<<std::left<<std::setw(52)<<nm<<std::right<<(ok?"PASS":"FAIL")<<"\n";
}

static double const gDa0 = 2.0e2;      // base rate; kap ramps it geometrically to 2e5

struct Vars { FFVar z,u,w,lam,mu,kap; };

static void build( FFGraph& DAG, OCFESLV& oc, Vars& V )
{
  V.z  = DAG.add_var("z");
  V.u  = DAG.add_var("u(z)");
  V.w  = DAG.add_var("w(z)");
  V.lam= DAG.add_var("lam");
  V.mu = DAG.add_var("mu");
  V.kap= DAG.add_var("kap");
  FFPartial OpP;

  oc.add_domain( V.z, FFDom( 0., 1., 4, FFDom::CGL, 5 ) );
  oc.add_state ( V.u, { V.z } );
  oc.add_state ( V.w, { V.z } );
  oc.add_input ( V.lam, 0.0, true );
  oc.add_input ( V.mu , 0.0, true );
  oc.add_input ( V.kap, 0.0, true );

  // kap in [0,1] maps geometrically onto Da in [Da0, 1000*Da0]
  FFVar Da  = gDa0 * pow( 1000.0, V.kap );
  FFVar PDE = OpP(OpP(V.u,V.z),V.z) - Da*V.u*V.w + V.lam*( 1.0 + V.z );
  FFVar ALG = V.w - ( 1.0 + V.mu*V.z );

  int const Z_INT = FFDom::ALL - FFDom::LB - FFDom::UB;
  oc.add_equation( PDE,      { V.z }, { Z_INT     }, OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( ALG,      { V.z }, { FFDom::ALL}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( V.u,      { V.z }, { FFDom::LB }, OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );
  oc.add_equation( V.u-1.0,  { V.z }, { FFDom::UB }, OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );

  oc.update_ref( V.u, [](OCFESLV::t_Coord const&){ return 0.0; } );
  oc.update_ref( V.w, [](OCFESLV::t_Coord const&){ return 1.0; } );

  oc.options.REDUCE.ORDER  = OCFESLV::Options::RED_FULL;
  oc.options.CLASSIFY.MODE      = OCFESLV::Options::CLASS_AUTO;
  oc.options.DISPLAY_LEVEL = 0;
}

#include <type_traits>
#include <stdexcept>

// ============================================================================================
// OCFE_exceptions_api (sandbox test, 2026-09-11, rev192 gate): exceptions thrown by the model layer
// and by OCFESLV must be caught as OCFESLV::Exceptions, OCBase::Exceptions and FFModel::Exceptions with
// unchanged codes and what() strings; default domain masks and AUTO roles must resolve as before.
// Output must be identical under rev191 and rev192.
// ============================================================================================
static void tryit( char const* nm, std::function<void()> const& f )
{
  try{ f(); std::cout << "  " << nm << ": no throw\n"; }
  catch( OCFESLV::Exceptions& e ){ std::cout << "  " << nm << ": OCFESLV::Exceptions ierr=" << e.ierr() << " what=" << e.what() << "\n"; }
  catch( ... ){ std::cout << "  " << nm << ": OTHER exception type\n"; }
}

int main()
{
  std::cout << "OCFESLV ** header: " << OCFESLV::HEADER_ID << "\n";
  FFGraph D; OCFESLV oc(&D); Vars V; build(D,oc,V);
  FFVar bogus = D.add_var("bogus");
  tryit( "update_ref(unknown var)",       [&]{ oc.update_ref( bogus, 1.0 ); } );
  tryit( "update_ref(unknown var, fun)",  [&]{ oc.update_ref( bogus, [](OCFESLV::t_Coord const&){ return 0.; } ); } );
  tryit( "ref(unknown var)",              [&]{ oc.ref( bogus ); } );
  tryit( "add_output point size mismatch",[&]{ oc.add_output( V.u, std::vector<FFVar>{ V.z }, std::vector<double>{ 0.1, 0.2 } ); } );
  tryit( "add_output mask size mismatch", [&]{ oc.add_output( V.u, std::vector<FFVar>{ V.z }, std::vector<int>{ 0, 0 } ); } );
  OCFESLV notset( &D );
  tryit( "deep_copy_from(not set up)",    [&]{ oc.deep_copy_from( notset ); } );
  try{ oc.ref( bogus ); } catch( OCBase::Exceptions& e ){ std::cout << "  caught as OCBase::Exceptions  ierr=" << e.ierr() << "\n"; }
  try{ oc.ref( bogus ); } catch( FFModel::Exceptions& e ){ std::cout << "  caught as FFModel::Exceptions ierr=" << e.ierr() << "\n"; }
  std::cout << "  same type: " << std::is_same<OCFESLV::Exceptions, OCBase::Exceptions>::value
            << std::is_same<FFModel::Exceptions, OCBase::Exceptions>::value << "\n";

  // rev198: the same object must also be catchable as a standard exception, with the SAME text, and an
  // FFDom error (thrown as a bare enum before rev198) must be catchable at all.
  try{ oc.ref( bogus ); }
  catch( std::runtime_error& e ){ std::cout << "  caught as std::runtime_error  what=" << e.what() << "\n"; }
  try{ oc.ref( bogus ); }
  catch( std::exception& e ){ std::cout << "  caught as std::exception      what=" << e.what() << "\n"; }
  try{ FFDom bad( 1.0, 0.0, 1, FFDom::LGL, 3 ); (void)bad; std::cout << "  FFDom bad bounds: no throw\n"; }
  catch( FFDom::Exceptions& e ){ std::cout << "  FFDom bad bounds: FFDom::Exceptions ierr=" << e.ierr() << "\n"; }
  catch( ... ){ std::cout << "  FFDom bad bounds: NOT catchable as a class (pre-rev198 bare enum)\n"; }

  // masks and AUTO role resolution
  size_t const n0 = oc.var_equation().size();
  oc.add_equation( V.w - 1.0, { V.z }, { FFDom::LB } );                      // AUTO + pure LB -> BOUNDARY
  oc.add_equation( V.w - 2.0, { V.z }, {} );                                 // AUTO + default mask ALL -> INTERIOR
  oc.add_equation( V.w - 3.0, { V.z }, { FFDom::ALL - FFDom::LB } );         // AUTO + ALL-LB -> INTERIOR
  oc.add_equation( V.w - 4.0, { V.z }, { FFDom::UB }, OCFESLV::EqnOptions( OCFESLV::EqnRole::DIAGNOSTIC, 0 ) );
  auto const& eq = oc.var_equation();
  for( size_t i = n0; i < eq.size(); ++i ){
    std::cout << "  eqn " << i << ": role=" << static_cast<int>( eq[i].opt->role )
              << " classify=" << eq[i].opt->participate_in_classification << " sat=" << OCFESLV::_sopt( *eq[i].opt ).receive_sat
              << " donate=" << OCFESLV::_sopt( *eq[i].opt ).donate_for_state_continuity << " masks:";
    for( auto const& [d,m] : eq[i].dom ) std::cout << " " << m;
    std::cout << "\n";
  }
  oc.add_output( V.u, { V.z }, std::vector<int>{} );                         // distributed output, default mask
  auto const& fo = oc.var_output().back();
  std::cout << "  output grid masks:"; for( auto const& [d,m] : fo.grid ) std::cout << " " << m; std::cout << "\n";
  // rev198: both layers derive from std::runtime_error, so ONE standard handler catches them all.
  auto as_std = []( char const* nm, std::function<void()> const& f ){
    try{ f(); std::cout << "  " << nm << ": no throw\n"; }
    catch( std::runtime_error& e ){ std::cout << "  " << nm << ": std::runtime_error what=" << e.what() << "\n"; }
    catch( ... ){ std::cout << "  " << nm << ": OTHER\n"; }
  };
  as_std( "model layer as std::runtime_error", [&]{ oc.ref( bogus ); } );
  as_std( "FFDom  as std::runtime_error",      [&]{ FFDom bad( 1.0, 0.0, 1, FFDom::LGL, 3 ); (void)bad; } );
  try{ oc.ref( bogus ); }
  catch( std::exception& e ){ std::cout << "  caught as std::exception   what=" << e.what() << "\n"; }
  std::cout << "  derives from std::runtime_error: "
            << std::is_base_of<std::runtime_error, OCFESLV::Exceptions>::value
            << std::is_base_of<std::runtime_error, FFDom::Exceptions>::value << "\n";

  std::cout << "  done\n";
  return 0;
}
