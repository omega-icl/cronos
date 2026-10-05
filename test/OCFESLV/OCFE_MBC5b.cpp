// OCFE_MBC5b.cpp  ---  Stage 3b, step 2 (STATE FORM): implicit wetting omega(z) as a solved {z}-state
// ============================================================================================
// Reinstates omega(z) as a SOLVED STATE on {z} (vs the prescribed-expression MBC5), closed locally by
// the Young-Laplace / pore-size-distribution relationship of Quek-Shah-Chachuat (2018):
//
//     critical pore radius (Young-Laplace):  delta_w(z) = 2*sigma_st*|cos(theta)| / dP(z)
//     wetting ratio (Eq.33 surrogate for Eq.9's cumulative pore-area integral):
//                    Zeta(delta_w) = 1/2 [ 1 + tanh(a0 + a1 u + a2 u^2 + a3 u^3) ],  u = delta_w/dbar
//                    (a0=7.966,a1=-14.08,a2=8.418,a3=-2.041; dbar=0.08um; Eq.32/Table 2)
//     closure (state):  OMG_CL = omega - [ (1-mu)*omega_fix + mu*Zeta(delta_w(z)) ] = 0  on {z}
//
// omega is add_state(.,{z}); its coefficients (rw,Drds,Drws; aMC,aPC,aMM,aPP; rDWf,rWLf; cgf,cgfV)
// feed the 2-D {rwc,z}/{rd,z} equations -- the exact cross-domain-state provisioning that previously
// threw "missing variable omega".  This is the MINIMAL test of whether a solved {z}-state in 2-D
// coefficients works; build the header with -DMC__OCFESLV_DEPEQN_PROBE to dump per-equation state
// dependencies and localize any drop.  mu ramps the closure (mu=0 pins omega=0.20, the MBC3 gate).
//
// dP(z) is PRESCRIBED (dP0 + dPslope*z/L).  The T-dependence of sigma_st (reaction->dT->wetting) and
// flow-coupled dP arrive with the energy balance; only then is the closure bidirectional.
//
// Gates: mu=0 reproduces MBC3 (Phi=4.369e-3); mu=1 solves omega -> Zeta(delta_w) ~ 0.040, and the
// SOLVED omega state must match the analytic surrogate; z-varying dP gives z-varying omega.
// ============================================================================================
#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <cmath>
#include <string>

#ifndef OCFE_OCFESLV_HEADER
#define OCFE_OCFESLV_HEADER "ocfeslv.hpp"
#endif
#include OCFE_OCFESLV_HEADER
using namespace mc;

static int g_pass=0,g_fail=0;
static void check_true(char const* n,bool ok){ (ok?g_pass:g_fail)++;
  std::cout<<"  "<<std::left<<std::setw(56)<<n<<"  "<<(ok?"PASS":"FAIL")<<"\n"; }

static std::vector<double> graded(int nel,double ratio){
  double w0=(ratio-1.0)/(std::pow(ratio,nel)-1.0);
  std::vector<double> w(nel); for(int i=0;i<nel;++i) w[i]=w0*std::pow(ratio,i); return w; }

// analytic surrogate value (for cross-check of the solved omega)
static double Zeta_ana(double dP,double sigma_st,double abscos,double dbar,
                       double a0,double a1,double a2,double a3){
  double u=(2.0*sigma_st*abscos/dP)/dbar;
  double arg=a0+a1*u+a2*u*u+a3*u*u*u;
  return 1.0/(1.0+std::exp(-2.0*arg)); }

struct Res { bool conv=false; double lam=0.0, mu=0.0, resid=1e30;
             double PhiA=0.0, PhiB=0.0, omgMin=0.0, omgMax=0.0, omgAna=0.0; };

static Res run(double dPslope)
{
  Res R;
  // ---- physical constants (Stage-2b/3) -----------------------------------------------------
  double const L=2.3,r1=225e-6,r2=550e-6,eta=0.40,epsm=0.41,tau=6.1; int const nf=209;
  double const r3=r2/std::sqrt(eta);
  double const Drl=r3-r2,Rref=r2;
  double const r1s=r1/Rref,r2s=r2/Rref,r3s=r3/Rref; (void)r3s;
  double const Drls=Drl/Rref;
  double const Dg=8.8e-8,DlC=1e-9,DlM=5e-10,DlP=5e-10;
  double const DwC=(epsm/tau)*DlC,DwM=(epsm/tau)*DlM,DwP=(epsm/tau)*DlP;
  double const DlCs=DlC/(Rref*Rref),DlMs=DlM/(Rref*Rref),DlPs=DlP/(Rref*Rref);
  double const H=3000.0;
  double const Z=0.9,Rg=8.314,Tg=293.0,P=54e5,yin=0.24;
  double const Cref=yin*P/(Z*Rg*Tg), Mref=0.37*1040/0.11916, Pref=0.06*1040/0.08614;
  double const MWg=yin*0.044+(1-yin)*0.028,rhog=P*MWg/(Z*Rg*Tg);
  double const vgIn=(3.0/3600.0)/(nf*M_PI*r1*r1*rhog), vlbar=(10e-3/3600.0)/(nf*M_PI*r3*r3*(1-eta));
  double const ClInS=0.10;
  double const partH=Z*Rg*Tg/H;
  double const kM=9.853e-2, kP=8.784e-1;
  double const bMC=kM*Mref, bPC=kP*Pref, bMM=kM*Cref, bPP=kP*Cref;

  // ---- wetting-closure constants (Quek-Shah-Chachuat 2018) --------------------------------
  double const dbar=0.08e-6;                                   // mean pore radius [m]  (Table 2)
  double const a0=7.966,a1=-14.08,a2=8.418,a3=-2.041;          // Eq.33 surrogate coefficients
  double const sigma_st=0.046;                                 // solvent surface tension [N/m]
  double const theta=92.5*M_PI/180.0, abscos=std::abs(std::cos(theta));  // contact angle
  double const dP0=30.0e3;                                     // transmembrane pressure [Pa]
  double const omega_fix=0.20;                                 // step-1 anchor (mu=0)
  R.omgAna = Zeta_ana(dP0,sigma_st,abscos,dbar,a0,a1,a2,a3);   // analytic omega at const dP0

  // ---- DAG / model -------------------------------------------------------------------------
  FFGraph DAG;
  FFVar z=DAG.add_var("z"),rd=DAG.add_var("rd"),rwc=DAG.add_var("rwc"),rl=DAG.add_var("rl");
  FFVar Cg=DAG.add_var("Cg"),Vg=DAG.add_var("Vg"),Cd=DAG.add_var("Cd");
  FFVar Cw=DAG.add_var("Cw"),Mw=DAG.add_var("Mw"),Pw=DAG.add_var("Pw");
  FFVar Cl=DAG.add_var("Cl"),Ml=DAG.add_var("Ml"),Pl=DAG.add_var("Pl");
  FFVar omg=DAG.add_var("omega");                 // wetting ratio -- STATE on {z}
  FFVar lam=DAG.add_var("lambda"), mu=DAG.add_var("mu");
  FFPartial OpP;

  // membrane coefficients as functions of the STATE omega (the cross-domain coupling under test:
  // a {z}-state feeding the coefficients of the 2-D {rwc,z}/{rd,z} equations).
  FFVar Drw_z = (r2-r1)*omg;
  FFVar Drd_z = (r2-r1)*(1.0-omg);
  FFVar rw_z  = r2 - (r2-r1)*omg;
  FFVar Drws_z= Drw_z/Rref, Drds_z=Drd_z/Rref, rws_z=rw_z/Rref;
  FFVar aMC_z =(Drw_z*Drw_z/DwC)*kM*Mref, aPC_z=(Drw_z*Drw_z/DwC)*kP*Pref;
  FFVar aMM_z =(Drw_z*Drw_z/DwM)*kM*Cref, aPP_z=(Drw_z*Drw_z/DwP)*kP*Cref;
  FFVar rDWf_z= DwC*Drd_z/(Dg*Drw_z);
  FFVar rWLf_z=(epsm/tau)*(Drl/Drw_z);
  FFVar cgf_z =(2.0*Dg/r1)/Drd_z;
  FFVar cgfV_z= cgf_z*yin/vgIn;

  // Young-Laplace + PSD surrogate closure: omega(z) is a STATE, pinned by OMG_CL at every z-node.
  // mu ramps the closure target from omega_fix (=step-1 anchor) to the full Zeta(delta_w(z)).
  FFVar dP_z    = dP0 + dPslope*(z/L);
  FFVar delta_w = (2.0*sigma_st*abscos)/dP_z;
  FFVar uu      = delta_w/dbar;
  FFVar arg     = a0 + a1*uu + a2*uu*uu + a3*uu*uu*uu;
  FFVar Zeta_z  = 1.0/(1.0+exp(-2.0*arg));                 // = (1+tanh(arg))/2
  FFVar OMG_CL  = 1e4*(omg - ( (1.0-mu)*omega_fix + mu*Zeta_z ));// wetting closure (algebraic, per z-node)

  auto cyl_d=[&](FFVar const& C,FFVar const& xi,double r0,FFVar const& dr){ return OpP(C,{xi,2})+(dr/(r0+xi*dr))*OpP(C,xi); };
  auto cyl_v=[&](FFVar const& C,FFVar const& xi,FFVar const& r0,FFVar const& dr){ return OpP(C,{xi,2})+(dr/(r0+xi*dr))*OpP(C,xi); };
  auto rad  =[&](FFVar const& C){ return OpP(C,{rl,2})/(Drls*Drls)+(1.0/(r2s+rl*Drls))*OpP(C,rl)/Drls; };

  FFVar GAS_C = vgIn*Vg*OpP(Cg,z) + cgf_z*(1.0-yin*Cg)*OpP(Cd,rd);
  FFVar GAS_V = OpP(Vg,z) + cgfV_z*OpP(Cd,rd);
  FFVar DRY_C = cyl_d(Cd,rd,r1s,Drds_z);
  FFVar WET_C = cyl_v(Cw,rwc,rws_z,Drws_z) - lam*( aMC_z*Cw*Mw + aPC_z*Cw*Pw );
  FFVar WET_M = cyl_v(Mw,rwc,rws_z,Drws_z) - lam*aMM_z*Cw*Mw;
  FFVar WET_P = cyl_v(Pw,rwc,rws_z,Drws_z) - lam*aPP_z*Cw*Pw;
  FFVar LIQ_C = vlbar*OpP(Cl,z) - DlCs*rad(Cl) + lam*( bMC*Cl*Ml + bPC*Cl*Pl );
  FFVar LIQ_M = vlbar*OpP(Ml,z) - DlMs*rad(Ml) + lam*bMM*Cl*Ml;
  FFVar LIQ_P = vlbar*OpP(Pl,z) - DlPs*rad(Pl) + lam*bPP*Cl*Pl;

  FFVar GAS_C_IN=Cg-1.0, GAS_V_IN=Vg-1.0;
  FFVar GD_VAL=Cd-Cg, DW_VAL=Cw-partH*Cd, DW_FLUX=OpP(Cd,rd)-rDWf_z*OpP(Cw,rwc);
  FFVar DW_MFLUX=OpP(Mw,rwc), DW_PFLUX=OpP(Pw,rwc);
  FFVar WL_CVAL=Cw-Cl, WL_CFLUX=rWLf_z*OpP(Cw,rwc)-OpP(Cl,rl);
  FFVar WL_MVAL=Mw-Ml, WL_MFLUX=rWLf_z*OpP(Mw,rwc)-OpP(Ml,rl);
  FFVar WL_PVAL=Pw-Pl, WL_PFLUX=rWLf_z*OpP(Pw,rwc)-OpP(Pl,rl);
  FFVar L_C_WALL=OpP(Cl,rl), L_M_WALL=OpP(Ml,rl), L_P_WALL=OpP(Pl,rl);
  FFVar L_C_IN=Cl-ClInS, L_M_IN=Ml-1.0, L_P_IN=Pl-1.0;

  OCFESLV oc(&DAG);
  oc.add_domain(z,  FFDom(0.,L,  3,FFDom::CGL,6));
  oc.add_domain(rd, FFDom(0.,1.0,3,FFDom::CGL,7));
  oc.add_domain(rwc,FFDom(0.0, graded(6,3.3), FFDom::CGL,8));
  oc.add_domain(rl, FFDom(0.0, graded(7,3.3), FFDom::CGL,8));
  oc.add_state(Cg,{z}); oc.add_state(Vg,{z}); oc.add_state(Cd,{rd,z});
  oc.add_state(Cw,{rwc,z}); oc.add_state(Mw,{rwc,z}); oc.add_state(Pw,{rwc,z});
  oc.add_state(Cl,{rl,z}); oc.add_state(Ml,{rl,z}); oc.add_state(Pl,{rl,z});
  oc.add_state(omg,{z});                              // wetting ratio state
  oc.add_input(lam,{}); oc.add_input(mu,{});
  oc.update_ref(Cg,[&](OCFESLV::t_Coord const&){return 1.0;});
  oc.update_ref(Vg,[&](OCFESLV::t_Coord const&){return 1.0;});
  oc.update_ref(Cd,[&](OCFESLV::t_Coord const&){return 1.0;});
  oc.update_ref(Cw,[&](OCFESLV::t_Coord const&){return partH;});
  oc.update_ref(Mw,[&](OCFESLV::t_Coord const&){return 1.0;});
  oc.update_ref(Pw,[&](OCFESLV::t_Coord const&){return 1.0;});
  oc.update_ref(Cl,[&](OCFESLV::t_Coord const&){return ClInS;});
  oc.update_ref(Ml,[&](OCFESLV::t_Coord const&){return 1.0;});
  oc.update_ref(Pl,[&](OCFESLV::t_Coord const&){return 1.0;});
  oc.update_ref(omg,[&](OCFESLV::t_Coord const&){return omega_fix;});

  int const Z_NO_UB=FFDom::ALL-FFDom::UB, Z_NO_LB=FFDom::ALL-FFDom::LB, R_INT=FFDom::ALL-FFDom::LB-FFDom::UB;
  OCFESLV::EqnOptions blk(OCFESLV::EqnRole::INTERIOR,0),bnd(OCFESLV::EqnRole::BOUNDARY,0),itf(OCFESLV::EqnRole::INTERFACE,0);

  oc.add_equation(GAS_V,   {z,rd},    {Z_NO_UB,FFDom::LB},                blk);
  oc.add_equation(GAS_C,   {z,rd},    {Z_NO_UB,FFDom::LB},                blk);
  oc.add_equation(GAS_V_IN,{z},       {FFDom::UB},                        bnd);
  oc.add_equation(GAS_C_IN,{z},       {FFDom::UB},                        bnd);
  oc.add_equation(DRY_C,   {rd,z},    {R_INT,FFDom::ALL},                 blk);
  oc.add_equation(WET_C,   {rwc,z},   {R_INT,FFDom::ALL},                 blk);
  oc.add_equation(WET_M,   {rwc,z},   {R_INT,FFDom::ALL},                 blk);
  oc.add_equation(WET_P,   {rwc,z},   {R_INT,FFDom::ALL},                 blk);
  oc.add_equation(LIQ_C,   {rl,z},    {R_INT,Z_NO_LB},                    blk);
  oc.add_equation(LIQ_M,   {rl,z},    {R_INT,Z_NO_LB},                    blk);
  oc.add_equation(LIQ_P,   {rl,z},    {R_INT,Z_NO_LB},                    blk);
  oc.add_equation(OMG_CL,  {z},       {FFDom::ALL},                       blk);   // wetting closure
  oc.add_equation(GD_VAL,  {rd,z},    {FFDom::LB,FFDom::ALL},             itf);
  oc.add_equation(DW_VAL,  {rd,rwc,z},{FFDom::UB,FFDom::LB,FFDom::ALL},   itf);
  oc.add_equation(DW_FLUX, {rd,rwc,z},{FFDom::UB,FFDom::LB,FFDom::ALL},   itf);
  oc.add_equation(DW_MFLUX,{rwc,z},   {FFDom::LB,FFDom::ALL},             itf);
  oc.add_equation(DW_PFLUX,{rwc,z},   {FFDom::LB,FFDom::ALL},             itf);
  oc.add_equation(WL_CVAL, {rwc,rl,z},{FFDom::UB,FFDom::LB,FFDom::ALL},   itf);
  oc.add_equation(WL_CFLUX,{rwc,rl,z},{FFDom::UB,FFDom::LB,FFDom::ALL},   itf);
  oc.add_equation(WL_MVAL, {rwc,rl,z},{FFDom::UB,FFDom::LB,FFDom::ALL},   itf);
  oc.add_equation(WL_MFLUX,{rwc,rl,z},{FFDom::UB,FFDom::LB,FFDom::ALL},   itf);
  oc.add_equation(WL_PVAL, {rwc,rl,z},{FFDom::UB,FFDom::LB,FFDom::ALL},   itf);
  oc.add_equation(WL_PFLUX,{rwc,rl,z},{FFDom::UB,FFDom::LB,FFDom::ALL},   itf);
  oc.add_equation(L_C_WALL,{rl,z},    {FFDom::UB,FFDom::ALL},             bnd);
  oc.add_equation(L_M_WALL,{rl,z},    {FFDom::UB,FFDom::ALL},             bnd);
  oc.add_equation(L_P_WALL,{rl,z},    {FFDom::UB,FFDom::ALL},             bnd);
  oc.add_equation(L_C_IN,  {rl,z},    {R_INT,FFDom::LB},                  bnd);
  oc.add_equation(L_M_IN,  {rl,z},    {R_INT,FFDom::LB},                  bnd);
  oc.add_equation(L_P_IN,  {rl,z},    {R_INT,FFDom::LB},                  bnd);

  oc.options.REDUCE.ORDER=OCFESLV::Options::RED_FULL;
  oc.options.CLASSIFY.MODE=OCFESLV::Options::CLASS_AUTO;
  oc.options.INTERFACE.TYPE=OCFESLV::Options::IC_AUTO;
  oc.options.INTERFACE.IMPOSITION=OCFESLV::Options::IC_STRONG;//WEAK;//TRACE;//STRONG;
  oc.options.INTERFACE.SAT_SIGMA0=10.0;
  oc.options.DISPLAY_LEVEL=0; oc.options.SOLVE.VERBOSE=true;//false;
#if defined(CRONOS__WITH_SPQR)
  oc.options.SOLVE.FACTORIZATION=OCFESLV::Options::SOLVE_SPQR;
#else
  oc.options.SOLVE.FACTORIZATION=OCFESLV::Options::SOLVE_SUPERLU;
#endif

  if(!oc.setup()){ std::cerr<<"  setup() FAILED\n"; return R; }
  std::vector<double> xv,inp;
  if(!oc.init(xv,inp,nullptr)){ std::cerr<<"  init() FAILED\n"; return R; }

  auto phiof=[&](double const* x){
    double const* ip=inp.empty()?nullptr:inp.data();
    auto ev=[&](FFVar const& V,OCFESLV::t_Coord const& pt){ return oc.eval_colloc<double>(V,pt,x,ip,nullptr); };
    OCFESLV::t_Coord gi;gi[z]=L; OCFESLV::t_Coord go;go[z]=0.0; OCFESLV::t_Coord vL;vL[z]=L; OCFESLV::t_Coord v0;v0[z]=0.0;
    return (r1/(2.0*L))*vgIn*Cref*( ev(Vg,vL)*ev(Cg,gi) - ev(Vg,v0)*ev(Cg,go) ); };

  OCFESLV::SolveReport rep;
  // ---- phase A: mu=0 (omega pinned to omega_fix), ramp lam 0->1  [reproduces step 1] --------
  oc.set_input_values(mu,{0.0},inp.data());
  double const lam_seq[]={0.0,0.01,0.03,0.06,0.1,0.2,0.35,0.5,0.65,0.8,1.0};
  for(double lv:lam_seq){
    oc.set_input_values(lam,{lv},inp.data());
    rep=oc.solve(xv.data(),inp.data(),nullptr); R.lam=lv; R.resid=rep.final_residual;
    if(!rep.converged){ R.conv=false; return R; }
  }
  R.PhiA=phiof(xv.data());

  // Expected omega(z) at the full closure (mu=1) = Zeta(delta_w(z)); the SOLVED state is read back
  // after phase B and must match this.  Kept here so it is reported even if phase B fails.
  { double omn=1e30,omx=-1e30;
    for(int i=0;i<=40;++i){ double zz=L*i/40.0;
      double o=Zeta_ana(dP0+dPslope*zz/L, sigma_st,abscos,dbar,a0,a1,a2,a3);
      omn=std::min(omn,o); omx=std::max(omx,o); }
    R.omgMin=omn; R.omgMax=omx; }

  // ---- phase B: adaptively ramp mu 0->1 (bisect on failure) ---------------------------------
  // A uniform low omega (const-dP: ~0.04 everywhere) is a stiffer continuation target than a
  // z-varying one -- the wet-liquid flux ratio rWLf ~ 1/omega quintuples across the whole column at
  // once -- so fixed mu steps stall near mu=1.  Bisection reaches it robustly.
  {
    std::vector<double> saved=xv; double cur=0.0, d=0.1; bool okB=true;
    while(cur < 1.0-1e-9){
      double trial=std::min(1.0, cur+d);
      oc.set_input_values(mu,{trial},inp.data());
      rep=oc.solve(xv.data(),inp.data(),nullptr); R.mu=trial; R.resid=rep.final_residual;
      if(rep.converged){ cur=trial; saved=xv; d=std::min(0.1,d*1.4); }
      else { xv=saved; d*=0.5; if(d<1e-4){ okB=false; break; } }
    }
    R.conv=okB;
    if(!okB) return R;
  }
  R.PhiB=phiof(xv.data());

  // read back the SOLVED omega(z) state (the cross-domain-state test): it must match the analytic
  // surrogate.  Overwrites the analytic min/max with the actually-solved values.
  { double const* ip=inp.empty()?nullptr:inp.data();
    double omn=1e30,omx=-1e30;
    for(int i=0;i<=40;++i){ OCFESLV::t_Coord pt; pt[z]=L*i/40.0;
      double o=oc.eval_colloc<double>(omg,pt,xv.data(),ip,nullptr);
      omn=std::min(omn,o); omx=std::max(omx,o); }
    R.omgMin=omn; R.omgMax=omx; }
  return R;
}

int main()
{
  std::cout<<"================================================================\n"
           <<"  Stage 3b step 2: implicit wetting omega(z) via Young-Laplace + PSD surrogate\n"
           <<"================================================================\n";
  double const PhiRef=4.369e-3;

  std::cout<<"  "<<std::left<<std::setw(22)<<"dP profile"<<std::setw(8)<<"conv"
           <<std::setw(12)<<"|r|"<<std::setw(13)<<"Phi(mu=0)"<<std::setw(13)<<"Phi(mu=1)"
           <<std::setw(18)<<"omega[min,max]"<<std::setw(11)<<"omega_ana"<<"\n";

  double const slopes[]={0.0, 6.0e3};   // const 30kPa ; +6kPa linear -> z-varying wetting
  std::vector<Res> out;
  for(double sl:slopes){
    Res R=run(sl); out.push_back(R);
    std::string convs=R.conv?std::string("y"):("n@lam"+std::to_string(R.lam)+"_mu"+std::to_string(R.mu));
    std::cout<<"  "<<std::left<<std::setw(22)<<(sl==0.0?"const 30kPa":"30+6*z/L kPa")<<std::setw(8)<<convs
             <<std::scientific<<std::setprecision(2)<<std::setw(12)<<R.resid
             <<std::setw(13)<<R.PhiA<<std::setw(13)<<R.PhiB
             <<std::fixed<<std::setprecision(4)<<"["<<R.omgMin<<","<<R.omgMax<<"]      "
             <<R.omgAna<<"\n";
  }

  std::cout<<"\n  ---- verdict ----\n";
  Res const& c=out[0];
  double gateErr=std::abs(c.PhiA-PhiRef)/PhiRef;
  std::cout<<"    gate: Phi(mu=0)="<<std::scientific<<std::setprecision(4)<<c.PhiA<<" vs MBC3 ref "<<PhiRef
           <<"  rel.err="<<gateErr<<"\n";
  check_true("mu=0 (omega pinned 0.20) reproduces MBC3 -- machinery correct", gateErr<2e-3);
  check_true("full Young-Laplace closure converges (mu=1, const & slope)", c.conv && out[1].conv);
  double omgErr = std::abs(c.omgMin-c.omgAna)/std::max(1e-12,c.omgAna);
  std::cout<<"    SOLVED omega(const dP)="<<std::fixed<<std::setprecision(4)<<c.omgMin
           <<" vs analytic surrogate "<<c.omgAna<<"  rel.err="<<std::scientific<<std::setprecision(2)<<omgErr<<"\n";
  check_true("solved omega state matches Young-Laplace surrogate (<1%)", c.conv && omgErr<1e-2);
  check_true("solved omega in the physical band (0.02-0.15)", c.omgMin>0.02 && c.omgMin<0.15);
  if(c.conv)
    std::cout<<"    physical wetting shifts operating point: Phi(0.20)="<<std::scientific<<std::setprecision(3)<<c.PhiA
             <<" -> Phi(omega="<<std::fixed<<std::setprecision(3)<<c.omgMin<<")="<<std::scientific<<std::setprecision(3)<<c.PhiB<<"\n";
  check_true("less wetting raises the flux (Phi(mu=1) > Phi(mu=0))", c.conv && c.PhiB>c.PhiA);
  check_true("z-varying dP gives z-varying omega (slope case)", out[1].conv && (out[1].omgMax-out[1].omgMin)>1e-4);

  std::cout<<"\n  RESULT: "<<g_pass<<" passed, "<<g_fail<<" failed\n";
  return g_fail?1:0;
}
