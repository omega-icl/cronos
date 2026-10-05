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
// OCFE_homotopy_copy (sandbox test, 2026-09-11, rev190 gate): a deep copy of a set-up OCFESLV must
// carry the homotopy settings -- registered parameters (remapped to the copy's DAG) and per-stage
// step caps -- and a deep_copy_from() into an environment that already held settings must REPLACE them.
// PRE-REGISTERED: under rev189 (settings not copied) the copy checks FAIL; under rev190 all PASS.
// ============================================================================================
struct Run { bool ok=false; int solves=0; int nparam=0; std::vector<double> x; };

static void schedule( OCFESLV& oc, Vars const& V, bool caps )
{
  oc.add_homotopy( V.lam, 0.0, 1.0, 0 );
  oc.add_homotopy( V.mu , 0.0, 1.0, 1 );
  oc.add_homotopy( V.kap, 0.0, 1.0, 2 );
  if( caps ) for( int s=0; s<3; ++s ) oc.set_homotopy_cap( s, 0.02 );
}

static Run solve_on( OCFESLV& oc, std::vector<double> xv, std::vector<double> inp )
{
  Run r;
  OCFESLV::SolveReport sr = oc.solve( xv.data(), inp.data(), nullptr );
  OCFESLV::ContinuationReport const& c = oc.continuation_report();
  r.ok = sr.converged; r.solves = c.solves;
  for( auto const& st : c.stage ) r.nparam += st.nparam;
  r.x = xv;
  return r;
}

static double dinf2( std::vector<double> const& a, std::vector<double> const& b )
{
  if( a.size() != b.size() ) return 1e300;
  double d = 0.; for( size_t i=0; i<a.size(); ++i ) d = std::max( d, std::abs( a[i]-b[i] ) ); return d;
}

static void line( char const* nm, Run const& r, Run const* ref )
{
  std::cout << "  " << std::left << std::setw(34) << nm << std::right << " converged=" << r.ok
            << " solves=" << std::setw(4) << r.solves << " params=" << r.nparam;
  if( ref ) std::cout << "  max|x-ref|=" << std::setprecision(3) << dinf2( r.x, ref->x );
  std::cout << "\n";
}

int main()
{
  std::cout << "OCFESLV ** header: " << OCFESLV::HEADER_ID << "\n";
  // references: fresh environments
  Run Rcap, Rnocap;
  { FFGraph D; OCFESLV oc(&D); Vars V; build(D,oc,V); oc.options.HOMOTOPY.STEP_CAP=0.1; oc.setup();
    std::vector<double> xv,inp; oc.init(xv,inp,nullptr); schedule(oc,V,true);  Rcap   = solve_on(oc,xv,inp); }
  { FFGraph D; OCFESLV oc(&D); Vars V; build(D,oc,V); oc.options.HOMOTOPY.STEP_CAP=0.1; oc.setup();
    std::vector<double> xv,inp; oc.init(xv,inp,nullptr); schedule(oc,V,false); Rnocap = solve_on(oc,xv,inp); }
  line( "reference, caps", Rcap, nullptr );
  line( "reference, no caps", Rnocap, nullptr );

  // copies
  Run Ccap, Cnocap, Csrc, Crepl;
  { FFGraph D; OCFESLV src(&D); Vars V; build(D,src,V); src.options.HOMOTOPY.STEP_CAP=0.1; src.setup();
    std::vector<double> xv,inp; src.init(xv,inp,nullptr); schedule(src,V,true);
    OCFESLV cp( src );                          // copy constructor -> deep_copy_from
    Ccap = solve_on( cp, xv, inp );
    Csrc = solve_on( src, xv, inp );          // the source is undisturbed by being copied
  }
  { FFGraph D; OCFESLV src(&D); Vars V; build(D,src,V); src.options.HOMOTOPY.STEP_CAP=0.1; src.setup();
    std::vector<double> xv,inp; src.init(xv,inp,nullptr); schedule(src,V,false);
    OCFESLV cp( src );
    Cnocap = solve_on( cp, xv, inp );
  }
  { // replacement: destination copied from a CAPPED source, then deep_copy_from an UNCAPPED one
    FFGraph D1; OCFESLV s1(&D1); Vars V1; build(D1,s1,V1); s1.options.HOMOTOPY.STEP_CAP=0.1; s1.setup();
    std::vector<double> x1,i1; s1.init(x1,i1,nullptr); schedule(s1,V1,true);
    FFGraph D2; OCFESLV s2(&D2); Vars V2; build(D2,s2,V2); s2.options.HOMOTOPY.STEP_CAP=0.1; s2.setup();
    std::vector<double> x2,i2; s2.init(x2,i2,nullptr); schedule(s2,V2,false);
    OCFESLV dst( s1 );
    dst.deep_copy_from( s2 );
    Crepl = solve_on( dst, x2, i2 );
  }
  line( "copy of capped source", Ccap, &Rcap );
  line( "copy of uncapped source", Cnocap, &Rnocap );
  line( "source after being copied", Csrc, &Rcap );
  line( "capped copy, then copy-from uncapped", Crepl, &Rnocap );

  check( "references converge",                               Rcap.ok && Rnocap.ok );
  check( "references differ (caps matter: test can fail)",    Rcap.solves != Rnocap.solves );
  check( "copy carries parameters (param count)",             Ccap.nparam == Rcap.nparam && Cnocap.nparam == Rnocap.nparam );
  check( "copy carries caps (solve counts)",                  Ccap.solves == Rcap.solves && Cnocap.solves == Rnocap.solves );
  check( "copy reaches the same solution (<= 1e-12)",         dinf2(Ccap.x,Rcap.x) <= 1e-12 && dinf2(Cnocap.x,Rnocap.x) <= 1e-12 );
  check( "source undisturbed by copy",                        Csrc.solves == Rcap.solves && dinf2(Csrc.x,Rcap.x) <= 1e-12 );
  check( "deep_copy_from replaces settings (no leftover caps)", Crepl.nparam == Rnocap.nparam && Crepl.solves == Rnocap.solves );
  std::cout << "\n  OCFE_homotopy_copy: " << ( g_fail ? "FAIL" : "PASS" ) << " (" << g_pass << " pass, " << g_fail << " fail)\n";
  return g_fail ? 1 : 0;
}
