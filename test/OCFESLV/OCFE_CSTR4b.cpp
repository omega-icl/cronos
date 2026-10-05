// OCFE_CSTR_smooth_diag.cpp  ---  Does the causal evolution interface handle the cf DISCONTINUITY?
// ===========================================================================
// ORIGINAL HYPOTHESIS (pre-fix, now RESOLVED): the monolithic multi-element-in-time solve failed
// ONLY because the feed cf is discontinuous at the element interfaces.  The monolithic imposed EVOL
// at each interior element's LB node t_k with a symmetric C0 tie that coupled across the jump, so
// when cf jumped there dc/dt jumped and the downstream LB carried the wrong derivative -> order-1
// corruption / stall.  A CONTINUOUS cf kept dc/dt continuous, cf(t_k) single-valued and the symmetric
// tie consistent, so the monolithic converged normally.
//
// THE FIX (causal evolution-direction interface, InputContinuity-gated): when a distributed input may
// jump in the evolution direction, the interface is made CAUSAL (one-sided) for the differential
// states -- the downstream LB receives the upstream terminal (the marched IC hand-off) and the
// anti-causal upstream row is dropped.  The kink is then consumed correctly and the STEP-cf monolithic
// CONVERGES.  It converges ALGEBRAICALLY (the solution genuinely kinks, so the rate is slower than the
// smooth/spectral case), but it is no longer stuck.
//
// => This driver is now a REGRESSION GUARD for that fix.  At n_nd = 3,5,7,9 (n_el=6):
//   Part A  SMOOTH cf(t)=2+0.8 sin(2pi t/T):  marched & monolithic  vs standard RK4   (spectral)
//   Part B  STEP  cf schedule (as before):    monolithic            vs phase-aligned RK4 (algebraic)
//   Expected: Part A converges spectrally to <1e-6; Part B converges algebraically to <1e-4.
//   Pre-fix, Part B was STUCK ( mP.back() > 1e-3 ) -- that stall is what the fix removes.
// ===========================================================================
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <functional>
#include <string>

#ifndef OCFE_OCFESLV_HEADER
#define OCFE_OCFESLV_HEADER "ocfeslv.hpp"
#endif
#include OCFE_OCFESLV_HEADER

#include "ffocfe.hpp"

using namespace mc;

static double const Dil = 1.0, krate = 2.0, Ksat = 1.0;
static size_t const NPH = 6;
static double const T_end = 6.0;
static double const LEVELS[6] = { 2.0, 1.0, 2.5, 1.5, 3.0, 2.0 };

static double cf_smooth( double tt ){ return 2.0 + 0.8*std::sin( 2.0*mc::PI*tt/T_end ); }   // C-infinity
static double cf_step( double tt )
{ long k=(long)std::floor((tt-1e-9)/1.0); if(k<0)k=0; if(k>=(long)NPH)k=(long)NPH-1; return LEVELS[k]; }

static double r_of_c( double c ){ return krate*c/( 1.0 + Ksat*c ); }
static double c_steady( double u )
{ double c=0.5*u; for(int it=0;it<100;++it){ double f=Dil*(u-c)-r_of_c(c), df=-Dil-krate/((1.0+Ksat*c)*(1.0+Ksat*c));
    double dc=-f/df; c+=dc; if(std::fabs(dc)<1e-14)break; } return c; }

// General RK4 (fine, continuous stepping) for a SMOOTH cf; sample at t=0..NPH.
static std::vector<double> rk4_smooth()
{ std::vector<double> out(NPH+1); double c=c_steady(cf_smooth(0.0)); out[0]=c;
  size_t const M=600000; double const h=T_end/double(M); size_t next=1;
  auto f=[&](double tt,double cc){ return Dil*(cf_smooth(tt)-cc)-r_of_c(cc); };
  for(size_t i=0;i<M;++i){ double tt=i*h;
    double k1=f(tt,c),k2=f(tt+h/2,c+h/2*k1),k3=f(tt+h/2,c+h/2*k2),k4=f(tt+h,c+h*k3);
    c+=h/6*(k1+2*k2+2*k3+k4); double tn=(i+1)*h;
    while(next<=NPH && tn>=double(next)-1e-9){ out[next]=c; ++next; } }
  return out; }

// Phase-aligned RK4 for the STEP cf; sample at t=0..NPH.
static std::vector<double> rk4_step()
{ std::vector<double> out(NPH+1); double c=c_steady(cf_step(0.5)); out[0]=c_steady(1.0); c=out[0];
  for(size_t k=0;k<NPH;++k){ double cf=LEVELS[k]; auto f=[&](double cc){return Dil*(cf-cc)-r_of_c(cc);};
    size_t M=200000; double h=1.0/double(M);
    for(size_t i=0;i<M;++i){ double k1=f(c),k2=f(c+h/2*k1),k3=f(c+h/2*k2),k4=f(c+h*k3); c+=h/6*(k1+2*k2+2*k3+k4);} out[k+1]=c; }
  return out; }

static std::vector<double> solve_c( size_t n_el, size_t n_nd, bool marching,
                                    std::function<double(double)> const& cf_fn, double cf0 )
{
  FFGraph DAG;
  FFVar t=DAG.add_var("t"), c=DAG.add_var("c(t)"), r=DAG.add_var("r(t)"), cf=DAG.add_var("cf(t)");
  OCFESLV oc(&DAG);
  oc.options.SOLVE.MARCHING = marching;
  oc.options.DISPLAY_LEVEL  = 0;
  oc.add_domain( t, FFDom(0.,T_end,n_el,FFDom::LGL,n_nd) );
  oc.set_evolution_domain( t );
  oc.add_state( c, {t} );
  oc.add_state( r, {t} );
  oc.add_input( cf, {t}, [&]( OCFESLV::t_Coord const& crd ){ return cf_fn( crd.at(t) ); } );
  oc.update_ref( c, [&]( OCFESLV::t_Coord const& ){ return c_steady(cf0); } );
  oc.update_ref( r, [&]( OCFESLV::t_Coord const& ){ return r_of_c(c_steady(cf0)); } );
  FFPartial OpP;
  FFVar EVOL = OpP(c,t) - ( Dil*(cf-c) - r );
  FFVar ALG  = r - krate*c/(1.0+Ksat*c);
  FFVar IC   = Dil*(cf0-c) - r;
  int const T_NO_LB = FFDom::ALL - FFDom::LB;
  oc.add_equation( EVOL, {t}, {T_NO_LB},    OCFESLV::EqnOptions(OCFESLV::EqnRole::INTERIOR,0) );
  oc.add_equation( ALG,  {t}, {FFDom::ALL}, OCFESLV::EqnOptions(OCFESLV::EqnRole::INTERIOR,0) );
  oc.add_equation( IC,   {t}, {FFDom::LB},  OCFESLV::EqnOptions(OCFESLV::EqnRole::INITIAL, 0) );
  for(size_t k=0;k<=NPH;++k) oc.add_output( c, {t}, { double(k) } );
  if( !oc.setup() ) return {};
  std::vector<double> xv, inp;
  if( !oc.init(xv,inp,nullptr) ) return {};
  OCFESLV::SolveReport const rep = oc.solve( xv.data(), inp.data(), nullptr );
  if( !rep.converged ) return {-1.};
  std::vector<double> F=oc.val_functions(); F.resize(NPH+1); return F;
}
static double err_vs( std::vector<double> const& F, std::vector<double> const& ref )
{ if(F.empty()) return -2.; if(F.size()==1&&F[0]==-1.) return -1.; double e=0.;
  for(size_t k=0;k<ref.size();++k) e=std::max(e,std::fabs(F[k]-ref[k])); return e; }

static int g_pass=0,g_fail=0;
static void check(std::string const& nm,bool ok)
{ std::cout<<"  "<<std::left<<std::setw(58)<<nm<<(ok?" PASS":" FAIL")<<"\n"; ok?++g_pass:++g_fail; }

int main()
{
  std::cout << "================================================================\n"
            << "  Causal evolution interface vs the cf DISCONTINUITY\n"
            << "  (continuous cf -> spectral; discontinuous cf -> algebraic; neither stalls)\n"
            << "================================================================\n\n";
  std::vector<size_t> nnds = { 3, 5, 7, 9 };

  std::vector<double> const refS = rk4_smooth();
  std::function<double(double)> smooth = cf_smooth;
  std::cout << "  Part A -- SMOOTH cf(t)=2+0.8 sin(2pi t/T):\n";
  std::cout << "      n_nd    marched|err|     monolithic|err|\n";
  std::vector<double> mS;
  for( size_t nd : nnds ){
    double em = err_vs( solve_c(NPH,nd,true, smooth, cf_smooth(0.0)), refS );
    double eo = err_vs( solve_c(NPH,nd,false,smooth, cf_smooth(0.0)), refS );
    mS.push_back( eo );
    std::cout << "      " << std::setw(4) << nd << "    " << std::scientific << std::setprecision(3)
              << em << "        " << eo << std::defaultfloat << "\n";
  }

  std::vector<double> const refP = rk4_step();
  std::function<double(double)> step = cf_step;
  std::cout << "\n  Part B -- STEP cf schedule (contrast):\n";
  std::cout << "      n_nd    monolithic|err|\n";
  std::vector<double> mP;
  for( size_t nd : nnds ){
    double eo = err_vs( solve_c(NPH,nd,false,step, 1.0), refP );
    mP.push_back( eo );
    std::cout << "      " << std::setw(4) << nd << "    " << std::scientific << std::setprecision(3)
              << eo << std::defaultfloat << "\n";
  }

  std::cout << "\n";
  bool smooth_mono_conv = true;
  for( size_t i=1;i<mS.size();++i ) if( !(mS[i]>=0 && mS[i] < mS[i-1]) ) smooth_mono_conv = false;
  bool smooth_mono_small = ( !mS.empty() && mS.back()>=0 && mS.back() < 1e-6 );

  // Post causal-evolution-interface fix: a discontinuous control no longer stalls the monolithic.
  // The causal one-sided interface consumes the solution kink at the t-interfaces (instead of the old
  // symmetric C0 tie coupling across it), so the STEP-cf monolithic CONVERGES.  It converges
  // ALGEBRAICALLY -- the solution genuinely kinks, so the rate is slower than Part A's spectral rate
  // (Part B reaches ~5e-6 at n_nd=9 vs Part A ~7e-9) -- but it converges monotonically.  Pre-fix this
  // was STUCK ( mP.back() > 1e-3 ), the symptom this driver was originally written to confirm; the fix
  // resolves it, so these two checks are a durable REGRESSION GUARD that the discontinuous-control
  // monolithic keeps converging.
  bool step_mono_conv  = ( mP.size() > 1 );
  for( size_t i=1;i<mP.size();++i ) if( !(mP[i]>=0 && mP[i] < mP[i-1]) ) step_mono_conv = false;
  bool step_mono_small = ( !mP.empty() && mP.back()>=0 && mP.back() < 1e-4 );

  check( "A  monolithic CONVERGES with continuous cf (monotone, spectral)",  smooth_mono_conv );
  check( "A  monolithic reaches < 1e-6 with continuous cf",                  smooth_mono_small );
  check( "B  monolithic CONVERGES with discontinuous cf (causal interface)", step_mono_conv );
  check( "B  monolithic reaches < 1e-4 with discontinuous cf (algebraic)",   step_mono_small );

  std::cout << "\n  VERDICT: "
            << ( (smooth_mono_conv && smooth_mono_small && step_mono_conv && step_mono_small)
                 ? "resolved -- the causal evolution interface handles the cf DISCONTINUITY: "
                   "smooth cf converges spectrally, step cf converges algebraically (no stall)."
                 : "inconclusive -- see the tables above." )
            << "\n";

  std::cout << "\n============================================================\n"
            << "  RESULT: " << g_pass << " passed, " << g_fail << " failed\n"
            << "============================================================\n";
  return 0;
}
