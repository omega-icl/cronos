// OCFE_CSTR_stepfix.cpp  ---  Verify the interior-side input-sampling fix.
// ===========================================================================
// After the fix (fill_block_reference / sample_input sample a DISCONTINUOUS input on each element's
// INTERIOR side at interior interfaces), the MONOLITHIC multi-element-in-time solve with a STEP cf
// should CONVERGE under p-refinement -- the downstream element's LB now carries its own feed
// levels[k], not the upstream levels[k-1].  Regression guards: the SMOOTH-cf monolithic and the
// marched path must be unchanged (the fix is a no-op for continuous inputs and single-element
// windows).
//
//   Cstep   monolithic STEP-cf error is monotone-decreasing in n_nd and reaches < 1e-4  (FIXED)
//   Csmooth monolithic SMOOTH-cf error still converges to < 1e-6                          (regression)
//   Cmarch  marched STEP-cf error still converges to < 1e-4                               (regression)
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

static double cf_smooth( double tt ){ return 2.0 + 0.8*std::sin( 2.0*mc::PI*tt/T_end ); }
static double cf_step( double tt )
{ long k=(long)std::floor((tt-1e-9)/1.0); if(k<0)k=0; if(k>=(long)NPH)k=(long)NPH-1; return LEVELS[k]; }
static double r_of_c( double c ){ return krate*c/(1.0+Ksat*c); }
static double c_steady( double u )
{ double c=0.5*u; for(int it=0;it<100;++it){ double f=Dil*(u-c)-r_of_c(c), df=-Dil-krate/((1.0+Ksat*c)*(1.0+Ksat*c));
    double dc=-f/df; c+=dc; if(std::fabs(dc)<1e-14)break; } return c; }

static std::vector<double> rk4_smooth()
{ std::vector<double> out(NPH+1); double c=c_steady(cf_smooth(0.0)); out[0]=c;
  size_t const M=600000; double const h=T_end/double(M); size_t next=1;
  auto f=[&](double tt,double cc){ return Dil*(cf_smooth(tt)-cc)-r_of_c(cc); };
  for(size_t i=0;i<M;++i){ double tt=i*h; double k1=f(tt,c),k2=f(tt+h/2,c+h/2*k1),k3=f(tt+h/2,c+h/2*k2),k4=f(tt+h,c+h*k3);
    c+=h/6*(k1+2*k2+2*k3+k4); double tn=(i+1)*h; while(next<=NPH && tn>=double(next)-1e-9){ out[next]=c; ++next; } }
  return out; }
static std::vector<double> rk4_step()
{ std::vector<double> out(NPH+1); double c=c_steady(1.0); out[0]=c;
  for(size_t k=0;k<NPH;++k){ double cf=LEVELS[k]; auto f=[&](double cc){return Dil*(cf-cc)-r_of_c(cc);};
    size_t M=200000; double h=1.0/double(M);
    for(size_t i=0;i<M;++i){ double k1=f(c),k2=f(c+h/2*k1),k3=f(c+h/2*k2),k4=f(c+h*k3); c+=h/6*(k1+2*k2+2*k3+k4);} out[k+1]=c; }
  return out; }

static std::vector<double> solve_c( size_t n_nd, bool marching,
                                    std::function<double(double)> const& cf_fn, double cf0 )
{
  FFGraph DAG;
  FFVar t=DAG.add_var("t"), c=DAG.add_var("c(t)"), r=DAG.add_var("r(t)"), cf=DAG.add_var("cf(t)");
  OCFESLV oc(&DAG);
  oc.options.SOLVE.MARCHING=marching; oc.options.DISPLAY_LEVEL=0;
  oc.add_domain( t, FFDom(0.,T_end,NPH,FFDom::LGL,n_nd) );
  oc.set_evolution_domain( t );
  oc.add_state( c, {t} ); oc.add_state( r, {t} );
  oc.add_input( cf, {t}, [&]( OCFESLV::t_Coord const& crd ){ return cf_fn( crd.at(t) ); } );
  oc.update_ref( c, [&]( OCFESLV::t_Coord const& ){ return c_steady(cf0); } );
  oc.update_ref( r, [&]( OCFESLV::t_Coord const& ){ return r_of_c(c_steady(cf0)); } );
  FFPartial OpP;
  FFVar EVOL=OpP(c,t)-(Dil*(cf-c)-r), ALG=r-krate*c/(1.0+Ksat*c), IC=Dil*(cf0-c)-r;
  int const T_NO_LB=FFDom::ALL-FFDom::LB;
  oc.add_equation( EVOL, {t}, {T_NO_LB},    OCFESLV::EqnOptions(OCFESLV::EqnRole::INTERIOR,0) );
  oc.add_equation( ALG,  {t}, {FFDom::ALL}, OCFESLV::EqnOptions(OCFESLV::EqnRole::INTERIOR,0) );
  oc.add_equation( IC,   {t}, {FFDom::LB},  OCFESLV::EqnOptions(OCFESLV::EqnRole::INITIAL, 0) );
  for(size_t k=0;k<=NPH;++k) oc.add_output( c, {t}, { double(k) } );
  if( !oc.setup() ) return {};
  std::vector<double> xv, inp;
  if( !oc.init(xv,inp,nullptr) ) return {};
  OCFESLV::SolveReport const rep=oc.solve(xv.data(),inp.data(),nullptr);
  if( !rep.converged ) return {-1.};
  std::vector<double> F=oc.val_functions(); F.resize(NPH+1); return F;
}
static double emax( std::vector<double> const& F, std::vector<double> const& ref )
{ if(F.empty()||(F.size()==1&&F[0]==-1.)) return -1.; double e=0.; for(size_t k=0;k<ref.size();++k) e=std::max(e,std::fabs(F[k]-ref[k])); return e; }

static int g_pass=0,g_fail=0;
static void check(std::string const& nm,bool ok){ std::cout<<"  "<<std::left<<std::setw(50)<<nm<<(ok?" PASS":" FAIL")<<"\n"; ok?++g_pass:++g_fail; }

int main()
{
  std::cout << "================================================================\n"
            << "  Interior-side input sampling fix -- verification\n"
            << "================================================================\n\n";
  std::vector<size_t> nnds={3,5,7,9};
  std::vector<double> const refP=rk4_step(), refS=rk4_smooth();
  std::function<double(double)> step=cf_step, smooth=cf_smooth;

  std::cout << "  STEP cf:   n_nd   marched|err|   monolithic|err|\n";
  std::vector<double> mstep, mmono;
  for(size_t nd:nnds){
    double em=emax(solve_c(nd,true, step,1.0),refP);
    double eo=emax(solve_c(nd,false,step,1.0),refP);
    mstep.push_back(em); mmono.push_back(eo);
    std::cout<<"             "<<std::setw(4)<<nd<<"   "<<std::scientific<<std::setprecision(3)<<em<<"      "<<eo<<std::defaultfloat<<"\n";
  }
  std::cout << "\n  SMOOTH cf: n_nd   monolithic|err|\n";
  std::vector<double> msm;
  for(size_t nd:nnds){
    double eo=emax(solve_c(nd,false,smooth,cf_smooth(0.0)),refS);
    msm.push_back(eo);
    std::cout<<"             "<<std::setw(4)<<nd<<"   "<<std::scientific<<std::setprecision(3)<<eo<<std::defaultfloat<<"\n";
  }

  std::cout << "\n";
  bool step_mono_mono=!mmono.empty(); for(size_t i=1;i<mmono.size();++i) if(!(mmono[i]>=0&&mmono[i]<mmono[i-1])) step_mono_mono=false;
  check("Cstep  monolithic STEP-cf now converges (monotone)", step_mono_mono);
  check("Cstep  monolithic STEP-cf finest < 1e-4",            !mmono.empty()&&mmono.back()>=0&&mmono.back()<1e-4);
  check("Csmooth monolithic SMOOTH-cf still < 1e-6 (regress)",!msm.empty()&&msm.back()>=0&&msm.back()<1e-6);
  check("Cmarch marched STEP-cf still < 1e-4 (regression)",   !mstep.empty()&&mstep.back()>=0&&mstep.back()<1e-4);

  std::cout << "\n============================================================\n"
            << "  RESULT: "<<g_pass<<" passed, "<<g_fail<<" failed  -- "<<(g_fail==0?"ALL PASS":"SOME FAILED")<<"\n"
            << "============================================================\n";
  return g_fail?1:0;
}
