// OCFE_MBC5c.cpp  ---  Stage 3b step 2: state-form wetting with RENORMALIZED flux continuity
// ============================================================================================
// Fixes the omega-in-denominator singularity found in MBC5b: the dry-wet and wet-liquid flux rows
// carried coefficients rDWf,rWLf ~ 1/Drw ~ 1/omega, which blow up as a *solved* omega drifts toward
// 0 (the stiff reaction pulls it down), diverging Newton.  Multiplying those rows through by Drw
// gives an equivalent equation (zero at the solution) with BOUNDED coefficients across omega in
// (0,1):
//     DW_FLUX :  Drw*dCd/drd - (DwC*Drd/Dg)*dCw/drwc = 0        [was dCd/drd - rDWf*dCw/drwc]
//     WL_*FLUX:  (eps/tau)*Drl*dXw/drwc - Drw*dXl/drl = 0       [was rWLf*dXw/drwc - dXl/drl]
//
// Also provides a CROSS-CHECK between the two wetting formulations: solve the (working) EXPRESSION
// form (omega a prescribed z-expression), transfer its converged solution into the STATE form
// (omega a solved DOF, omega:=Zeta), and evaluate the state-form residual WITHOUT solving.  Because
// at omega=Zeta the two forms' physical equations are identical and OMG_CL=0, the residual must be
// tiny if the equations (and the renormalization) are coded consistently -- an equation-level
// correctness check independent of whether the state form converges.
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

// ---- wetting-closure constants (Quek-Shah-Chachuat 2018), file scope, reused by both build paths ----
static double const gL=2.3, gdP0=30.0e3, gsig=0.046, gdbar=0.08e-6;
static double const ga0=7.966, ga1=-14.08, ga2=8.418, ga3=-2.041;
static double const gabscos=std::abs(std::cos(92.5*M_PI/180.0));
static double const gomega_fix=0.20;

static double Zeta_ana(double dP){
  double u=(2.0*gsig*gabscos/dP)/gdbar;
  double arg=ga0+ga1*u+ga2*u*u+ga3*u*u*u;
  return 1.0/(1.0+std::exp(-2.0*arg)); }

// ---- shared model build (omega as solved STATE if as_state, else prescribed z-expression) --------
static void build_mbc( FFGraph& DAG, OCFESLV& oc, bool as_state, double dPslope,
  FFVar& z, FFVar& rd, FFVar& rwc, FFVar& rl,
  FFVar& Cg, FFVar& Vg, FFVar& Cd, FFVar& Cw, FFVar& Mw, FFVar& Pw,
  FFVar& Cl, FFVar& Ml, FFVar& Pl, FFVar& omg, FFVar& lam, FFVar& mu )
{
  double const L=gL,r1=225e-6,r2=550e-6,eta=0.40,epsm=0.41,tau=6.1; int const nf=209;
  double const r3=r2/std::sqrt(eta);
  double const Drl=r3-r2,Rref=r2;
  double const r1s=r1/Rref,r2s=r2/Rref;
  double const Drls=Drl/Rref;
  double const Dg=8.8e-8,DlC=1e-9,DlM=5e-10,DlP=5e-10;
  double const DwC=(epsm/tau)*DlC,DwM=(epsm/tau)*DlM,DwP=(epsm/tau)*DlP;
  double const DlCs=DlC/(Rref*Rref),DlMs=DlM/(Rref*Rref),DlPs=DlP/(Rref*Rref);
  double const H=3000.0;
  double const Z=0.9,Rg=8.314,Tg=293.0,P=54e5,yin=0.24;
  double const Cref=yin*P/(Z*Rg*Tg), Mref=0.37*1040/0.11916, Pref=0.06*1040/0.08614;
  double const MWg=yin*0.044+(1-yin)*0.028,rhog=P*MWg/(Z*Rg*Tg);
  double const vgIn=(3.0/3600.0)/(nf*M_PI*r1*r1*rhog), vlbar=(10e-3/3600.0)/(nf*M_PI*r3*r3*(1-eta));
  double const ClInS=0.10, partH=Z*Rg*Tg/H;
  double const kM=9.853e-2, kP=8.784e-1;
  double const bMC=kM*Mref, bPC=kP*Pref, bMM=kM*Cref, bPP=kP*Cref;

  z=DAG.add_var("z"); rd=DAG.add_var("rd"); rwc=DAG.add_var("rwc"); rl=DAG.add_var("rl");
  Cg=DAG.add_var("Cg"); Vg=DAG.add_var("Vg"); Cd=DAG.add_var("Cd");
  Cw=DAG.add_var("Cw"); Mw=DAG.add_var("Mw"); Pw=DAG.add_var("Pw");
  Cl=DAG.add_var("Cl"); Ml=DAG.add_var("Ml"); Pl=DAG.add_var("Pl");
  lam=DAG.add_var("lambda"); mu=DAG.add_var("mu");
  FFPartial OpP;

  // Young-Laplace + PSD surrogate wetting target Zeta(delta_w(z))
  FFVar dP_z    = gdP0 + dPslope*(z/L);
  FFVar delta_w = (2.0*gsig*gabscos)/dP_z;
  FFVar uu      = delta_w/gdbar;
  FFVar arg     = ga0 + ga1*uu + ga2*uu*uu + ga3*uu*uu*uu;
  FFVar Zeta_z  = 1.0/(1.0+exp(-2.0*arg));
  FFVar target  = (1.0-mu)*gomega_fix + mu*Zeta_z;

  if( as_state ) omg = DAG.add_var("omega");         // solved {z}-state
  else           omg = target;                       // prescribed z-expression (working form)

  // membrane coefficients as functions of omega
  FFVar Drw_z = (r2-r1)*omg;
  FFVar Drd_z = (r2-r1)*(1.0-omg);
  FFVar rw_z  = r2 - (r2-r1)*omg;
  FFVar Drws_z= Drw_z/Rref, Drds_z=Drd_z/Rref, rws_z=rw_z/Rref;
  FFVar aMC_z =(Drw_z*Drw_z/DwC)*kM*Mref, aPC_z=(Drw_z*Drw_z/DwC)*kP*Pref;
  FFVar aMM_z =(Drw_z*Drw_z/DwM)*kM*Cref, aPP_z=(Drw_z*Drw_z/DwP)*kP*Cref;
  FFVar cgf_z =(2.0*Dg/r1)/Drd_z;
  FFVar cgfV_z= cgf_z*yin/vgIn;
  // RENORMALIZED flux coefficients (bounded across omega): DWf=rDWf*Drw, WLf=rWLf*Drw (constant)
  FFVar DWf_z = DwC*Drd_z/Dg;             // = rDWf*Drw  -> (DwC/Dg)*(r2-r1)*(1-omega)
  double const WLf = (epsm/tau)*Drl;      // = rWLf*Drw  -> (eps/tau)*Drl  (omega-independent)

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
  FFVar GD_VAL=Cd-Cg, DW_VAL=Cw-partH*Cd;
  FFVar DW_FLUX = Drw_z*OpP(Cd,rd) - DWf_z*OpP(Cw,rwc);        // renormalized (x Drw)
  FFVar DW_MFLUX=OpP(Mw,rwc), DW_PFLUX=OpP(Pw,rwc);
  FFVar WL_CVAL=Cw-Cl, WL_CFLUX=WLf*OpP(Cw,rwc)-Drw_z*OpP(Cl,rl);   // renormalized (x Drw)
  FFVar WL_MVAL=Mw-Ml, WL_MFLUX=WLf*OpP(Mw,rwc)-Drw_z*OpP(Ml,rl);
  FFVar WL_PVAL=Pw-Pl, WL_PFLUX=WLf*OpP(Pw,rwc)-Drw_z*OpP(Pl,rl);
  FFVar L_C_WALL=OpP(Cl,rl), L_M_WALL=OpP(Ml,rl), L_P_WALL=OpP(Pl,rl);
  FFVar L_C_IN=Cl-ClInS, L_M_IN=Ml-1.0, L_P_IN=Pl-1.0;
  FFVar OMG_CL = omg - target;

  oc.add_domain(z,  FFDom(0.,L,  3,FFDom::CGL,6));
  oc.add_domain(rd, FFDom(0.,1.0,3,FFDom::CGL,7));
  oc.add_domain(rwc,FFDom(0.0, graded(6,3.3), FFDom::CGL,8));
  oc.add_domain(rl, FFDom(0.0, graded(7,3.3), FFDom::CGL,8));
  oc.add_state(Cg,{z}); oc.add_state(Vg,{z}); oc.add_state(Cd,{rd,z});
  oc.add_state(Cw,{rwc,z}); oc.add_state(Mw,{rwc,z}); oc.add_state(Pw,{rwc,z});
  oc.add_state(Cl,{rl,z}); oc.add_state(Ml,{rl,z}); oc.add_state(Pl,{rl,z});
  if( as_state ) oc.add_state(omg,{z});
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
  if( as_state ) oc.update_ref(omg,[&](OCFESLV::t_Coord const&){return gomega_fix;});

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
  if( as_state ) oc.add_equation(OMG_CL, {z}, {FFDom::ALL},              blk);
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
  oc.options.INTERFACE.IMPOSITION=OCFESLV::Options::IC_STRONG;
  oc.options.INTERFACE.SAT_SIGMA0=10.0;
  oc.options.DISPLAY_LEVEL=1; oc.options.SOLVE.VERBOSE=true;
#if defined(CRONOS__WITH_SPQR)
  oc.options.SOLVE.FACTORIZATION=OCFESLV::Options::SOLVE_SPQR;
#else
  oc.options.SOLVE.FACTORIZATION=OCFESLV::Options::SOLVE_SUPERLU;
#endif
  oc.setup();
}

static double phi_of( OCFESLV& oc, FFVar const& z, FFVar const& Cg, FFVar const& Vg,
                      double const* xv, double const* inp, double L )
{
  double const r1=225e-6,r2=550e-6,eta=0.40; int const nf=209;
  double const Z=0.9,Rg=8.314,Tg=293.0,P=54e5,yin=0.24;
  double const Cref=yin*P/(Z*Rg*Tg);
  double const rhog=P*(yin*0.044+(1-yin)*0.028)/(Z*Rg*Tg);
  double const vgIn=(3.0/3600.0)/(nf*M_PI*r1*r1*rhog); (void)r2;(void)eta;
  auto ev=[&](FFVar const& V,double zz){ OCFESLV::t_Coord pt; pt[z]=zz; return oc.eval_colloc<double>(V,pt,xv,inp,nullptr); };
  return (r1/(2.0*L))*vgIn*Cref*( ev(Vg,L)*ev(Cg,L) - ev(Vg,0.0)*ev(Cg,0.0) );
}

struct Res { bool conv=false; double lam=0.0,mu=0.0,resid=1e30,PhiA=0.0,PhiB=0.0,omgMin=0.0,omgMax=0.0; };

static Res run(bool as_state, double dPslope)
{
  Res R;
  FFGraph DAG; OCFESLV oc(&DAG);
  FFVar z,rd,rwc,rl,Cg,Vg,Cd,Cw,Mw,Pw,Cl,Ml,Pl,omg,lam,mu;
  build_mbc(DAG,oc,as_state,dPslope, z,rd,rwc,rl,Cg,Vg,Cd,Cw,Mw,Pw,Cl,Ml,Pl,omg,lam,mu);
  std::vector<double> xv,inp; if(!oc.init(xv,inp,nullptr)){ std::cerr<<"  init failed\n"; return R; }

  // phase A checkpoint (mu=0) for the gate
  OCFESLV::SolveReport rep;
  double const lam_seq[]={0.0,0.01,0.03,0.06,0.1,0.2,0.35,0.5,0.65,0.8,1.0};
  oc.set_input_values(mu,{0.0},inp.data());
  bool okA=true;
  for(double lv:lam_seq){ oc.set_input_values(lam,{lv},inp.data());
    rep=oc.solve(xv.data(),inp.data(),nullptr); R.lam=lv; R.resid=rep.final_residual;
    if(!rep.converged){ okA=false; break; } }
  if(okA) R.PhiA=phi_of(oc,z,Cg,Vg,xv.data(),inp.data(),gL);
  { double omn=1e30,omx=-1e30; for(int i=0;i<=40;++i){ double zz=gL*i/40.0;
      double o=Zeta_ana(gdP0+dPslope*zz/gL); omn=std::min(omn,o); omx=std::max(omx,o);} R.omgMin=omn;R.omgMax=omx; }
  if(!okA){ R.conv=false; return R; }

  // phase B: adaptive mu 0->1
  std::vector<double> saved=xv; double cur=0.0,d=0.1; bool okB=true;
  while(cur<1.0-1e-9){ double trial=std::min(1.0,cur+d); oc.set_input_values(mu,{trial},inp.data());
    rep=oc.solve(xv.data(),inp.data(),nullptr); R.mu=trial; R.resid=rep.final_residual;
    if(rep.converged){ cur=trial; saved=xv; d=std::min(0.1,d*1.4); }
    else { xv=saved; d*=0.5; if(d<1e-4){ okB=false; break; } } }
  R.conv=okB; if(!okB) return R;
  R.PhiB=phi_of(oc,z,Cg,Vg,xv.data(),inp.data(),gL);
  if(as_state){ double const* ip=inp.empty()?nullptr:inp.data(); double omn=1e30,omx=-1e30;
    for(int i=0;i<=40;++i){ OCFESLV::t_Coord pt; pt[z]=gL*i/40.0;
      double o=oc.eval_colloc<double>(omg,pt,xv.data(),ip,nullptr); omn=std::min(omn,o);omx=std::max(omx,o);} R.omgMin=omn;R.omgMax=omx; }
  return R;
}

int main()
{
  std::cout<<"================================================================\n"
           <<"  Stage 3b step 2c: renormalized flux (x Drw) -- state vs expression cross-check\n"
           <<"================================================================\n";
  double const PhiRef=4.369e-3;
  double const slopes[]={0.0, 6.0e3};

  // Solve BOTH formulations (solved {z}-state vs prescribed expression) and compare their converged
  // solutions.  With the algebraic-classification fix the state form now converges, so the two must
  // land on the same solution -- the direct, robust form of "are the residuals similarly small?".
  std::cout<<"  "<<std::left<<std::setw(16)<<"dP profile"<<std::setw(8)<<"form"<<std::setw(8)<<"conv"
           <<std::setw(12)<<"|r|"<<std::setw(13)<<"Phi(mu=0)"<<std::setw(13)<<"Phi(mu=1)"<<"omega[min,max]\n";
  std::vector<Res> St, Ex;
  for(double sl:slopes){
    Res rs=run(true, sl);  St.push_back(rs);
    Res re=run(false,sl);  Ex.push_back(re);
    for(int k=0;k<2;++k){ Res const& R=(k==0?rs:re);
      std::string cs=R.conv?std::string("y"):("n@l"+std::to_string(R.lam)+"_m"+std::to_string(R.mu));
      std::cout<<"  "<<std::left<<std::setw(16)<<(k?std::string(""):(sl==0.0?"const 30kPa":"30+6*z/L kPa"))
               <<std::setw(8)<<(k==0?"state":"expr")<<std::setw(8)<<cs
               <<std::scientific<<std::setprecision(2)<<std::setw(12)<<R.resid
               <<std::setw(13)<<R.PhiA<<std::setw(13)<<R.PhiB
               <<std::fixed<<std::setprecision(4)<<"["<<R.omgMin<<","<<R.omgMax<<"]\n"; }
  }

  std::cout<<"\n  ---- verdict ----\n";
  double gateErr=std::abs(St[0].PhiA-PhiRef)/PhiRef;
  std::cout<<"    gate: state Phi(mu=0)="<<std::scientific<<std::setprecision(4)<<St[0].PhiA<<" vs MBC3 ref "<<PhiRef<<"  rel.err="<<gateErr<<"\n";
  check_true("state mu=0 reproduces MBC3", gateErr<2e-3);
  check_true("state form converges (mu=1, const & slope)", St[0].conv && St[1].conv);
  check_true("expression form converges (mu=1, const & slope)", Ex[0].conv && Ex[1].conv);

  bool phi_match=true, omg_match=true;
  for(size_t i=0;i<2;++i){
    if(St[i].conv && Ex[i].conv){
      double dphi=std::abs(St[i].PhiB-Ex[i].PhiB)/std::max(1e-30,std::abs(Ex[i].PhiB));
      double domg=std::max(std::abs(St[i].omgMin-Ex[i].omgMin),std::abs(St[i].omgMax-Ex[i].omgMax));
      std::cout<<"    "<<(i==0?"const":"slope")<<": |dPhi|/Phi="<<std::scientific<<std::setprecision(2)<<dphi
               <<"  |domega|="<<domg<<"   (state |r|="<<St[i].resid<<", expr |r|="<<Ex[i].resid<<")\n";
      if(dphi>=1e-3) phi_match=false;
      if(domg>=1e-3) omg_match=false;
    } else { phi_match=false; omg_match=false; }
  }
  check_true("state & expression forms agree in Phi (<0.1%)", phi_match);
  check_true("state & expression forms agree in omega (<1e-3)", omg_match);
  if(St[0].conv) std::cout<<"    physical wetting: Phi(0.20)="<<std::scientific<<std::setprecision(3)<<St[0].PhiA
                          <<" -> Phi("<<std::fixed<<std::setprecision(3)<<St[0].omgMin<<")="<<std::scientific<<std::setprecision(3)<<St[0].PhiB<<"\n";
  check_true("less wetting raises the flux (Phi(mu=1) > Phi(mu=0))", St[0].conv && St[0].PhiB>St[0].PhiA);

  std::cout<<"\n  RESULT: "<<g_pass<<" passed, "<<g_fail<<" failed\n";
  return g_fail?1:0;
}
