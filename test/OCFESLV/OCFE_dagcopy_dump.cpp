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
// OCFE_dagcopy_dump (sandbox instrument, 2026-09-11, rev191 gate): dumps the working DAG of deep
// copies -- every variable (type, index) and every operation (type, info, operand ids, result ids) --
// for a reduced, classified model with point and distributed outputs and reference functions.
// rev191 moves the copy's single FFGraph::insert (model roots followed by the solver's roots) into
// FFModel::_deep_copy_model; the dump under rev190 and rev191 must be IDENTICAL.  (A first rev191 design
// used two inserts and id lookups; it aborted -- insert() re-creates auxiliary nodes with new ids.)
// ============================================================================================
static void vid( std::ostream& os, FFVar const* pv )
{
  auto const id = pv->id();
  os << id.first << ":";
  if( id.first == FFVar::VAR || id.first == FFVar::AUX ) os << id.second;
  else os << std::setprecision(17) << *pv;
}

static void dump( char const* tag, FFGraph const* dag )
{
  std::cout << "== " << tag << "  nvar=" << dag->Vars().size() << "  nop=" << dag->Ops().size() << "\n";
  size_t k = 0;
  for( auto const* pv : dag->Vars() ){ std::cout << "  V" << k++ << " "; vid( std::cout, pv ); std::cout << "\n"; }
  k = 0;
  for( auto const* po : dag->Ops() ){
    std::cout << "  O" << k++ << " t=" << po->type << " i=" << po->info << " in:";
    for( auto const* pv : po->varin ){ std::cout << " "; vid( std::cout, pv ); }
    std::cout << " out:";
    for( auto const* pv : po->varout ){ std::cout << " "; vid( std::cout, pv ); }
    std::cout << "\n";
  }
}

int main()
{
  std::cout << "OCFESLV ** header: " << OCFESLV::HEADER_ID << "\n";
  FFGraph D; OCFESLV src(&D); Vars V; build(D,src,V);
  src.add_output( V.u, { V.z }, { 0.5 } );                 // point output
  src.add_output( V.w, { V.z }, { FFDom::ALL } );          // distributed output
  if( !src.setup() ){ std::cout << "setup FAILED\n"; return 1; }
  std::cout << "source: nvar=" << src.dag()->Vars().size() << " nop=" << src.dag()->Ops().size()
            << "  eqns=" << src.var_equation().size() << " states=" << src.var_state().size()
            << " outputs=" << src.var_output().size() << "\n";
  OCFESLV c1( src );                 dump( "copy of source", c1.dag() );
  OCFESLV c2( c1 );                  dump( "copy of copy", c2.dag() );
  FFGraph D2; OCFESLV a( &D2 );  a = src;  dump( "assignment over fresh env", a.dag() );
  std::cout << "c1: eqns=" << c1.var_equation().size() << " states=" << c1.var_state().size()
            << " outputs=" << c1.var_output().size() << " ref(u)=" << std::setprecision(17) << c1.ref( V.u )
            << " ref(z)=" << c1.ref( V.z ) << "\n";
  return 0;
}
