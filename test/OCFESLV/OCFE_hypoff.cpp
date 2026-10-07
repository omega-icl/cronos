// OCFE_hypoff.cpp -- GATE for AUTO.HYP_CLOSURE (WORKPLAN 3.C, 2026-10-07).  The automatic outflow closure of hyperbolic
// blocks is an OPTION again (it was the environment-only CRONOS_AUTO_HYP_CLOSURE).  Measured before: with it off and the
// outflow face left open, setup succeeded with only a warning and the solve "converged" to an answer 10% wrong (an
// UNDER-determined system).  Now setup refuses (HYP_CLOSURE_MISSING) -- the rule of consistent initial data for a
// high-index DAE.  Scalar advection u_t + a u_z = 0, exact u = sin(2 pi (z - a t)), inflow condition at the inflow face:
//   on, either speed sign: correct;  off, outflow open: REFUSED, either sign;
//   off, the model closing the face itself (its PDE extended there): sets up, = on;  on + that closure: = on (rev317);
//   off, the face closed by a SEPARATE row written for the purpose: not refused, = on (a first version of the refusal
//   counted only the PDE-extension idiom and refused OCFE_PDE14/15's manual oracles).
#include <cstdio>
#include <cmath>
#include <sstream>
#include <iostream>
#include "ocfeslv.hpp"
using namespace mc;  typedef OCFESLV::EqnRole Role;  typedef OCFESLV::EqnOptions EO;
static int nfail = 0;
static void check( char const* what, bool ok ){ std::printf( "  %-66s %s\n", what, ok? "PASS": "FAIL" ); nfail += !ok; }
struct R { bool ok = false; OCFESLV::SetupStatus st = OCFESLV::SetupStatus::OK; double v = NAN, ex = NAN; std::string msg; };
static R run( double a0, bool hyp, int own_closure /*0 none, 1 PDE extended, 2 a separate closure row*/ ){
  double const T = 0.25, PI = 3.14159265358979;
  FFGraph G; FFPartial OpP; FFEval OpE;  int const NLB = FFDom::ALL - FFDom::LB, INN = FFDom::ALL - FFDom::LB - FFDom::UB;
  FFVar t = G.add_var( "t" ), z = G.add_var( "z" ), u = G.add_var( "u" ), a = G.add_var( "a" );
  OCFESLV S( &G );
  S.add_domain( t, FFDom( 0., T, 4, FFDom::LGR, 4 ) );  S.add_domain( z, FFDom( 0., 1., 8, FFDom::LGL, 6 ) );
  S.set_evolution_domain( t );  S.add_state( u, {t, z} );  S.add_input( a );  S.update_ref( a, a0 );
  int const out = a0 > 0? FFDom::UB: FFDom::LB, in = a0 > 0? FFDom::LB: FFDom::UB;
  // the PDE on the interior; with own_closure, extended to the OUTFLOW face too (the OCFE_PDE6/7/8/10 idiom)
  S.add_equation( OpP( u, t ) + a * OpP( u, z ), {t, z}, {NLB, own_closure == 1? ( FFDom::ALL - in ): INN}, EO( Role::INTERIOR ) );
  if( own_closure == 2 )   // a closure row WRITTEN FOR THE PURPOSE (the OCFE_PDE14/15 manual-oracle idiom), not the PDE's mask
    S.add_equation( OpP( u, t ) + a * OpP( u, z ), {t, z}, {NLB, out},
                    EO( Role::INTERIOR, 0, OCFESLV::Options::IC_AUTO, /*classify=*/false, /*sat=*/true ) );   // as OCFE_PDE14:
    // an INTERIOR row at the face, kept out of classification -- a BOUNDARY row there would be a condition at the
    // outflow end, which the incoming-BC validator refuses
  S.add_equation( u - sin( 2.*PI*z ), {t, z}, {(int)FFDom::LB, (int)FFDom::ALL}, EO( Role::INITIAL ) );
  S.add_equation( a0 > 0? u - sin( -2.*PI*a*t ): u - sin( 2.*PI*( 1. - a*t ) ), {t, z}, {NLB, in}, EO( Role::BOUNDARY ) );
  S.add_output( OpE( u, {{t,1},{z,1}}, {{t,T},{z,0.5}} ) );
  S.options.AUTO.HYP_CLOSURE = hyp;  S.options.SOLVE.MARCHING = false;  S.options.DISPLAY_LEVEL = 0;
  std::ostringstream os; auto* old = std::cerr.rdbuf( os.rdbuf() );
  R r;  r.ok = S.setup();  r.st = S.setup_status();  r.msg = os.str();
  if( r.ok ){ std::vector<double> var, inp; S.init( var, inp ); S.solve( var.data(), inp.data(), nullptr ); r.v = S.val_functions()[0]; }
  std::cerr.rdbuf( old );
  r.ex = std::sin( 2.*PI*( 0.5 - a0*T ) );
  (void)out;
  return r;
}
// The 2x2 system u_t + A u_z = 0, A = [[0,1],[1,0]] (speeds +1, -1; characteristics w+ = u1+u2 right-going,
// w- = u1-u2 left-going), exact w+ = sin(2 pi (z - t)), w- = 0.5 sin(2 pi (z + t)).  Inflow data on the ENTERING
// characteristic at each end (w+ at LB, w- at UB); the OUTGOING one closed by hand (w- at LB, w+ at UB) with rows that
// are characteristic COMBINATIONS of the PDEs (PDE1 -+ PDE2) -- the OCFE_PDE14/15 manual-oracle idiom, which a refusal
// counting only the PDE-extension idiom wrongly refuses.  manual: 0 none, 1 combination rows.
static R run2( bool hyp, int manual ){
  double const T = 0.25, PI = 3.14159265358979;
  FFGraph G; FFPartial OpP; FFEval OpE;  int const NLB = FFDom::ALL - FFDom::LB, INN = FFDom::ALL - FFDom::LB - FFDom::UB;
  FFVar t = G.add_var( "t" ), z = G.add_var( "z" ), u1 = G.add_var( "u1" ), u2 = G.add_var( "u2" );
  OCFESLV S( &G );
  S.add_domain( t, FFDom( 0., T, 4, FFDom::LGR, 4 ) );  S.add_domain( z, FFDom( 0., 1., 8, FFDom::LGL, 6 ) );
  S.set_evolution_domain( t );  S.add_state( u1, {t, z} );  S.add_state( u2, {t, z} );
  FFVar const P1 = OpP( u1, t ) + OpP( u2, z ), P2 = OpP( u2, t ) + OpP( u1, z );
  S.add_equation( P1, {t, z}, {NLB, INN}, EO( Role::INTERIOR ) );
  S.add_equation( P2, {t, z}, {NLB, INN}, EO( Role::INTERIOR ) );
  FFVar const wp0 = sin( 2.*PI*z ), wm0 = 0.5*sin( 2.*PI*z );
  S.add_equation( u1 - 0.5*( wp0 + wm0 ), {t, z}, {(int)FFDom::LB, (int)FFDom::ALL}, EO( Role::INITIAL ) );
  S.add_equation( u2 - 0.5*( wp0 - wm0 ), {t, z}, {(int)FFDom::LB, (int)FFDom::ALL}, EO( Role::INITIAL ) );
  S.add_equation( u1 + u2 - sin( -2.*PI*t ),            {t, z}, {NLB, (int)FFDom::LB}, EO( Role::BOUNDARY ) );   // w+ enters at LB
  S.add_equation( u1 - u2 - 0.5*sin( 2.*PI*( 1. + t ) ), {t, z}, {NLB, (int)FFDom::UB}, EO( Role::BOUNDARY ) );   // w- enters at UB
  if( manual ){
    EO const clo( Role::INTERIOR, 0, OCFESLV::Options::IC_AUTO, /*classify=*/false, /*sat=*/true );
    S.add_equation( P1 - P2, {t, z}, {NLB, (int)FFDom::LB}, clo );   // w- leaves at LB
    S.add_equation( P1 + P2, {t, z}, {NLB, (int)FFDom::UB}, clo );   // w+ leaves at UB
  }
  S.add_output( OpE( u1, {{t,1},{z,1}}, {{t,T},{z,0.5}} ) );
  S.options.AUTO.HYP_CLOSURE = hyp;  S.options.SOLVE.MARCHING = false;  S.options.DISPLAY_LEVEL = 0;
  std::ostringstream os; auto* old = std::cerr.rdbuf( os.rdbuf() );
  R r;  r.ok = S.setup();  r.st = S.setup_status();  r.msg = os.str();
  if( r.ok ){ std::vector<double> var, inp; S.init( var, inp ); S.solve( var.data(), inp.data(), nullptr ); r.v = S.val_functions()[0]; }
  std::cerr.rdbuf( old );
  r.ex = 0.5*( std::sin( 2.*PI*( 0.5 - T ) ) + 0.5*std::sin( 2.*PI*( 0.5 + T ) ) );
  return r;
}
int main(){
  std::printf( "OCFE_hypoff -- AUTO.HYP_CLOSURE on / off\n" );
  for( double a0 : { +1., -1. } ){
    char buf[128];  char const* sg = a0 > 0? "a=+1": "a=-1";
    R on = run( a0, true, 0 ), off = run( a0, false, 0 ), offc = run( a0, false, 1 ), onc = run( a0, true, 1 ), offs = run( a0, false, 2 );
    std::snprintf( buf, sizeof buf, "%s on: sets up, u(T,.5) = exact to the discretisation (1e-3; %.1e)", sg, std::fabs( on.v - on.ex ) );
    check( buf, on.ok && std::fabs( on.v - on.ex ) < 1e-3 );
    std::snprintf( buf, sizeof buf, "%s off, outflow face open: REFUSED, HYP_CLOSURE_MISSING", sg );
    check( buf, !off.ok && off.st == OCFESLV::SetupStatus::HYP_CLOSURE_MISSING );
    std::snprintf( buf, sizeof buf, "%s ...and the message names the outflow face (%s)", sg, a0 > 0? "UB": "LB" );
    check( buf, off.msg.find( a0 > 0? "z UB": "z LB" ) != std::string::npos );
    std::snprintf( buf, sizeof buf, "%s off, the model closes the face itself: sets up, = on (%.1e)", sg, std::fabs( offc.v - on.v ) );
    check( buf, offc.ok && std::fabs( offc.v - on.v ) < 1e-10 );
    std::snprintf( buf, sizeof buf, "%s on + the model's own closure: not double-closed, = on (%.1e)", sg, std::fabs( onc.v - on.v ) );
    check( buf, onc.ok && std::fabs( onc.v - on.v ) < 1e-10 );
    std::snprintf( buf, sizeof buf, "%s off, a SEPARATE closure row at the face: not refused, = on (%.1e)", sg, std::fabs( offs.v - on.v ) );
    check( buf, offs.ok && std::fabs( offs.v - on.v ) < 1e-10 );
  }
  { char buf[128];
    R on = run2( true, 0 ), off = run2( false, 0 ), man = run2( false, 1 );
    std::snprintf( buf, sizeof buf, "2x2 on: sets up, u1(T,.5) = exact to the discretisation (1e-3; %.1e)", std::fabs( on.v - on.ex ) );
    check( buf, on.ok && std::fabs( on.v - on.ex ) < 1e-3 );
    check( "2x2 off, both outflow faces open: REFUSED, HYP_CLOSURE_MISSING", !off.ok && off.st == OCFESLV::SetupStatus::HYP_CLOSURE_MISSING );
    std::snprintf( buf, sizeof buf, "2x2 off, closed by characteristic COMBINATIONS: not refused, = on (%.1e)", std::fabs( man.v - on.v ) );
    check( buf, man.ok && std::fabs( man.v - on.v ) < 1e-8 ); }
  std::printf( "  OCFE_hypoff: %d failed\n", nfail );
  return nfail? 1: 0;
}
