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

// ---- schedule via solve_homotopy: stages supplied by the caller --------------------------
static bool run_homotopy( int s_lam, int s_mu, int s_kap, std::vector<double>& xout,
                          OCFESLV::ContinuationReport& rep, double cap )
{
  FFGraph DAG; OCFESLV oc(&DAG); Vars V; build(DAG,oc,V);
  oc.options.HOMOTOPY.STEP_CAP = cap;
  if( !oc.setup() ) return false;
  std::vector<double> xv,inp;
  if( !oc.init(xv,inp,nullptr) ) return false;
  oc.add_homotopy( V.lam, 0.0, 1.0, s_lam );
  oc.add_homotopy( V.mu , 0.0, 1.0, s_mu  );
  oc.add_homotopy( V.kap, 0.0, 1.0, s_kap );      // linear in kap; Da is geometric in the model
  // Single entry point: solve() honours the registered schedule (options.HOMOTOPY.ENABLE).
  OCFESLV::SolveReport sr = oc.solve( xv.data(), inp.data(), nullptr );
  rep  = oc.continuation_report();     // per-stage detail of the run just performed
  xout = xv;
  return sr.converged;
}

// ---- per-stage cap + auto-destage probe ---------------------------------------------------
// Runs an all-stage-0 (simultaneous) schedule at a caller-chosen cap, optionally overriding the cap
// on individual stages.  Reports whether OCFESLV had to AUTO-DESTAGE, i.e. split the stage into
// sequential single-parameter sub-stages after failing to move the parameters together.
static bool run_capped( int s_lam, int s_mu, int s_kap, double cap,
                        std::vector<std::pair<int,double>> const& overrides,
                        std::vector<double>& xout, int& destaged, int& nsolve )
{
  FFGraph DAG; OCFESLV oc(&DAG); Vars V; build(DAG,oc,V);
  oc.options.HOMOTOPY.STEP_CAP = cap;
  if( !oc.setup() ) return false;
  std::vector<double> xv,inp;
  if( !oc.init(xv,inp,nullptr) ) return false;
  oc.add_homotopy( V.lam, 0.0, 1.0, s_lam );
  oc.add_homotopy( V.mu , 0.0, 1.0, s_mu  );
  oc.add_homotopy( V.kap, 0.0, 1.0, s_kap );
  for( auto const& o : overrides ) oc.set_homotopy_cap( o.first, o.second );
  OCFESLV::SolveReport sr = oc.solve( xv.data(), inp.data(), nullptr );
  OCFESLV::ContinuationReport const& r = oc.continuation_report();
  destaged = 0; nsolve = r.solves;
  for( auto const& st : r.stage ) destaged += st.destaged;
  xout = xv;
  return sr.converged;
}

// ---- hand-rolled staircase, the pattern solve_homotopy replaces ---------------------------
static bool run_manual( std::vector<double>& xout, int& nsolve )
{
  FFGraph DAG; OCFESLV oc(&DAG); Vars V; build(DAG,oc,V);
  if( !oc.setup() ) return false;
  std::vector<double> xv,inp;
  if( !oc.init(xv,inp,nullptr) ) return false;
  FFVar const* par[3] = { &V.lam, &V.mu, &V.kap };
  nsolve = 0;
  for( int k=0;k<3;++k ){
    std::vector<double> saved = xv;
    double cur=0., d=0.1;
    while( std::abs(cur-1.) > 1e-9 ){
      double const trial = cur + std::min(d, 1.-cur);
      oc.set_input_values( *par[k], {trial}, inp.data() );
      xv = saved;
      OCFESLV::SolveReport sr = oc.solve( xv.data(), inp.data(), nullptr );
      ++nsolve;
      if( sr.converged ){ cur=trial; saved=xv; d=std::min(0.1,d*1.4); }
      else{ xv=saved; oc.set_input_values(*par[k],{cur},inp.data()); d*=0.5; if(d<1e-5) return false; }
    }
  }
  xout = xv;
  return true;
}

// ---- direct solve at full physics, no continuation ---------------------------------------
static bool run_direct( double& resid )
{
  FFGraph DAG; OCFESLV oc(&DAG); Vars V; build(DAG,oc,V);
  if( !oc.setup() ) return false;
  std::vector<double> xv,inp;
  if( !oc.init(xv,inp,nullptr) ) return false;
  oc.set_input_values(V.lam,{1.0},inp.data());
  oc.set_input_values(V.mu ,{1.0},inp.data());
  oc.set_input_values(V.kap,{1.0},inp.data());
  OCFESLV::SolveReport sr = oc.solve( xv.data(), inp.data(), nullptr );   // no schedule registered
  resid = sr.final_residual;
  return sr.converged;
}

static double dinf( std::vector<double> const& a, std::vector<double> const& b )
{
  if( a.size()!=b.size() || a.empty() ) return 1e30;
  double m=0.; for(size_t i=0;i<a.size();++i) m=std::max(m,std::fabs(a[i]-b[i]));
  return m;
}

static void report( char const* nm, OCFESLV::ContinuationReport const& r )
{
  std::cout<<"    "<<std::left<<std::setw(14)<<nm<<std::right
           <<" stages="<<r.stage.size()<<"  solves="<<std::setw(3)<<r.solves
           <<"  acc="<<std::setw(3)<<r.accepts<<"  bt="<<std::setw(2)<<r.backtracks
           <<"  pred+="<<std::setw(3)<<r.pred_used<<"  pred-="<<std::setw(3)<<r.pred_rejected
           <<"  its="<<std::setw(4)<<r.iterations
           <<"  |r|="<<std::scientific<<std::setprecision(2)<<r.final_residual<<std::fixed<<"\n";
}

int main()
{
  std::cout<<"================================================================\n"
           <<"  OCFESLV declarative homotopy (phase 1): stages, adaptive step, report\n"
           <<"================================================================\n";

  // ---- contrast: no continuation ----
  double rdirect=1e30;
  bool const dok = run_direct( rdirect );
  std::cout<<"\n  direct solve (no continuation): converged="<<(dok?"yes":"NO")
           <<"  |r|="<<std::scientific<<std::setprecision(2)<<rdirect<<std::fixed<<"\n";

  // ---- the three schedules the stage grammar must express ----
  std::cout<<"\n  schedules via solve_homotopy:\n";
  std::vector<double> xSim, xStair, xPart, xMan;
  OCFESLV::ContinuationReport rSim, rStair, rPart;
  bool const okSim   = run_homotopy( 0,0,0, xSim,   rSim,   0.1 );  // simultaneous
  if(okSim)   report("simultaneous", rSim);
  bool const okStair = run_homotopy( 0,1,2, xStair, rStair, 0.1 );  // staircase
  if(okStair) report("staircase",    rStair);
  bool const okPart  = run_homotopy( 0,0,1, xPart,  rPart,  0.1 );  // partial grouping
  if(okPart)  report("partial",      rPart);

  int nman=0;
  bool const okMan = run_manual( xMan, nman );
  std::cout<<"    "<<std::left<<std::setw(14)<<"manual"<<std::right
           <<" (hand-rolled staircase)  solves="<<nman<<"\n";

  std::cout<<"\n  ---- verdict ----\n";
  check("simultaneous schedule converges",        okSim  );
  check("staircase schedule converges",           okStair);
  check("partial-grouping schedule converges",    okPart );
  check("hand-rolled staircase converges",        okMan  );

  if( okStair && okMan )
    check("solve_homotopy == hand-rolled staircase", dinf(xStair,xMan) < 1e-9 );
  if( okSim && okStair )
    check("simultaneous and staircase: SAME root",   dinf(xSim,xStair) < 1e-9 );
  if( okPart && okStair )
    check("partial and staircase: SAME root",        dinf(xPart,xStair) < 1e-9 );

  if( okStair ){
    check("staircase used 3 stages",  rStair.stage.size()==3 );
    bool book = true;
    for( auto const& st : rStair.stage ) book &= ( st.solves == st.accepts + st.backtracks );
    check("report bookkeeping consistent", book );
    bool one = true;
    for( auto const& st : rStair.stage ) one &= ( st.nparam == 1 );
    check("staircase stages drive 1 param each", one );
  }
  if( okSim ){
    check("simultaneous used 1 stage", rSim.stage.size()==1 );
    check("simultaneous stage drives 3 params",
          !rSim.stage.empty() && rSim.stage[0].nparam==3 );
    check("simultaneous costs fewer solves than staircase",
          okStair && rSim.solves < rStair.solves );
  }

  // ---- per-stage cap override ----
  std::cout<<"\n  per-stage cap override (staircase, stage 0 coarse / stage 2 fine):\n";
  std::vector<double> xCap; int dCap=0, nCap=0;
  bool const okCap = run_capped( 0,1,2, 0.1, {{0,0.5},{2,0.02}}, xCap, dCap, nCap );
  std::cout<<"    converged="<<(okCap?"yes":"NO")<<"  solves="<<nCap<<"  destaged="<<dCap<<"\n";
  check("per-stage cap: schedule converges", okCap );
  if( okCap && okStair )
    check("per-stage cap: SAME root as uniform cap", dinf(xCap,xStair) < 1e-9 );

  // ---- auto-destage: force a coarse simultaneous jump and see whether OCFESLV recovers ----
  std::cout<<"\n  auto-destage probe (all stage 0, deliberately coarse cap=1.0):\n";
  std::vector<double> xDes; int dDes=0, nDes=0;
  bool const okDes = run_capped( 0,0,0, 1.0, {}, xDes, dDes, nDes );
  std::cout<<"    converged="<<(okDes?"yes":"NO")<<"  solves="<<nDes
           <<"  destaged="<<dDes<<( dDes? "  (auto-destage FIRED)" : "  (not needed)" )<<"\n";
  check("coarse cap: converges (directly or via destage)", okDes );
  if( okDes && okStair )
    check("coarse cap: SAME root as the staircase",  dinf(xDes,xStair) < 1e-9 );

  std::cout<<"\n  RESULT: "<<g_pass<<" passed, "<<g_fail<<" failed\n";
  return g_fail?1:0;
}
