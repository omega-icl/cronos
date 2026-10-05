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


// ============================================================================================
// OCFE_homotopy_reset (sandbox test, 2026-09-11, rev189 gate): reset_homotopy() must remove the
// registered parameters AND the per-stage step caps.  No corpus driver calls reset_homotopy().
//   A  fresh environment, staircase schedule, global cap              -> reference (solves, x)
//   B  same, but first: register + tight per-stage caps + reset_homotopy(), then re-register
//      without caps                                                    -> must EQUAL A
//   C  caps set and NOT reset (positive control of the test)          -> must DIFFER from A
// ============================================================================================
static bool run_case( char which, std::vector<double>& xout, int& nsolve, int& nparam )
{
  FFGraph DAG; OCFESLV oc(&DAG); Vars V; build(DAG,oc,V);
  oc.options.HOMOTOPY.STEP_CAP = 0.1;
  if( !oc.setup() ) return false;
  std::vector<double> xv,inp;
  if( !oc.init(xv,inp,nullptr) ) return false;
  if( which == 'B' ){
    oc.add_homotopy( V.lam, 0.0, 1.0, 0 );
    oc.add_homotopy( V.mu , 0.0, 1.0, 1 );
    oc.add_homotopy( V.kap, 0.0, 1.0, 2 );
    for( int s=0; s<3; ++s ) oc.set_homotopy_cap( s, 0.02 );
    oc.reset_homotopy();
  }
  oc.add_homotopy( V.lam, 0.0, 1.0, 0 );
  oc.add_homotopy( V.mu , 0.0, 1.0, 1 );
  oc.add_homotopy( V.kap, 0.0, 1.0, 2 );
  if( which == 'C' )
    for( int s=0; s<3; ++s ) oc.set_homotopy_cap( s, 0.02 );
  OCFESLV::SolveReport sr = oc.solve( xv.data(), inp.data(), nullptr );
  OCFESLV::ContinuationReport const& r = oc.continuation_report();
  nsolve = r.solves; nparam = 0;
  for( auto const& st : r.stage ) nparam += st.nparam;
  xout = xv;
  return sr.converged;
}

int main()
{
  std::cout << "OCFESLV ** header: " << OCFESLV::HEADER_ID << "\n";
  std::vector<double> xA, xB, xC; int nA=0, nB=0, nC=0, pA=0, pB=0, pC=0;
  bool const okA = run_case( 'A', xA, nA, pA );
  bool const okB = run_case( 'B', xB, nB, pB );
  bool const okC = run_case( 'C', xC, nC, pC );
  double dAB = 0.; for( size_t i=0; i<xA.size() && i<xB.size(); ++i ) dAB = std::max( dAB, std::abs(xA[i]-xB[i]) );
  std::cout << std::setprecision(10)
            << "  A: converged=" << okA << " solves=" << nA << " params=" << pA << "\n"
            << "  B: converged=" << okB << " solves=" << nB << " params=" << pB << "  max|xB-xA|=" << dAB << "\n"
            << "  C: converged=" << okC << " solves=" << nC << " params=" << pC << "\n";
  check( "A converges",                                      okA );
  check( "B converges",                                      okB );
  check( "B == A: reset removed registrations (param count)", pB == pA );
  check( "B == A: reset removed caps (solve count)",          nB == nA );
  check( "B == A: same solution",                             xA.size()==xB.size() && dAB == 0. );
  check( "C != A: caps change the run (test can fail)",       nC != nA );
  std::cout << "\n  OCFE_homotopy_reset: " << ( g_fail ? "FAIL" : "PASS" ) << " (" << g_pass << " pass, " << g_fail << " fail)\n";
  return g_fail ? 1 : 0;
}
