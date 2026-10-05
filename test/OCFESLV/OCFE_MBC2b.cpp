// OCFE_MBC2b.cpp  ---  Stage 2b-i: three-species (CO2/MDEA/PZ) dimensionless MBC
// ===========================================================================
// Extends the validated 2a dimensionless scaffold to EXPLICIT MDEA + PZ, with the
// real reaction STRUCTURE  R_CO2 = -(k_M C_CO2 C_MDEA + k_P C_CO2 C_PZ),
// R_MDEA = -k_M C_CO2 C_MDEA, R_PZ = -k_P C_CO2 C_PZ  (paper Appendix A).
// 9 states: Cg,Vg (gas) | Cd (dry) | Cw,Mw,Pw (wet) | Cl,Ml,Pl (liquid).
// 2b-i uses MODERATE rates (wet Da_M~0.2, Da_P~0.4) to confirm the 9-state model
// and the k_P>k_M distinction solve; the REAL Arrhenius rates (Da_P~3.6e6, a film
// regime) are ramped via the lam continuation input in 2b-ii.  Non-dimensionalized
// throughout (concentrations by Cref/Mref/Pref, velocity by vgIn, geometry by Rref);
// cylindrical Laplacian; row-scaled interface fluxes; IC_STRONG.
// ===========================================================================
#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <cmath>

#ifndef OCFE_OCFESLV_HEADER
#define OCFE_OCFESLV_HEADER "ocfeslv.hpp"
#endif
#include OCFE_OCFESLV_HEADER

using namespace mc;

static int g_pass=0,g_fail=0;
static void check_true(char const* n,bool ok){ (ok?g_pass:g_fail)++;
  std::cout<<"  "<<std::left<<std::setw(48)<<n<<"  "<<(ok?"PASS":"FAIL")<<"\n"; }

int main()
{
  double const L=2.3,r1=225e-6,r2=550e-6,eta=0.40,epsm=0.41,tau=6.1; int const nf=209;
  double const omega=0.20,rw=r2-(r2-r1)*omega,r3=r2/std::sqrt(eta);
  double const Drd=rw-r1,Drw=r2-rw,Drl=r3-r2,Rref=r2;
  double const r1s=r1/Rref,r2s=r2/Rref,rws=rw/Rref;
  double const Drds=Drd/Rref,Drws=Drw/Rref,Drls=Drl/Rref;
  double const Dg=8.8e-8,DlC=1e-9,DlM=5e-10,DlP=5e-10;
  double const DwC=(epsm/tau)*DlC,DwM=(epsm/tau)*DlM,DwP=(epsm/tau)*DlP;
  double const DlCs=DlC/(Rref*Rref),DlMs=DlM/(Rref*Rref),DlPs=DlP/(Rref*Rref);
  double const H=3000.0;
  double const Z=0.9,Rg=8.314,Tg=293.0,P=54e5,yin=0.24;
  double const Cref=yin*P/(Z*Rg*Tg), Mref=0.37*1040/0.11916, Pref=0.06*1040/0.08614;
  double const MWg=yin*0.044+(1-yin)*0.028,rhog=P*MWg/(Z*Rg*Tg);
  double const vgIn=(3.0/3600.0)/(nf*M_PI*r1*r1*rhog), vlbar=(10e-3/3600.0)/(nf*M_PI*r3*r3*(1-eta));
  double const ClInS=0.10;                          // lean CO2 loading (scaled)
  double const partH=Z*Rg*Tg/H;
  // 2b-i MODERATE rates (real Arrhenius in 2b-ii)
  double const kM=9.853e-2, kP=8.784e-1;   // 1e5x moderate -> Da_P~40000 (plateau check)
  double const cgf=(2.0*Dg/r1)/Drd, cgfV=cgf*yin/vgIn, rDWf=DwC*Drd/(Dg*Drw), rWLf=(epsm/tau)*(Drl/Drw);
  // scaled reaction coefficients (lam multiplies)
  double const aMC=(Drw*Drw/DwC)*kM*Mref, aPC=(Drw*Drw/DwC)*kP*Pref;   // wet CO2 sinks
  double const aMM=(Drw*Drw/DwM)*kM*Cref, aPP=(Drw*Drw/DwP)*kP*Cref;   // wet amine sinks
  double const bMC=kM*Mref, bPC=kP*Pref, bMM=kM*Cref, bPP=kP*Cref;     // liquid sinks

  std::cout<<"================================================================\n"
           <<"  Stage 2b-i: three-species (CO2/MDEA/PZ) dimensionless MBC\n"
           <<"  Mref/Cref="<<std::fixed<<std::setprecision(2)<<Mref/Cref
           <<"  Pref/Cref="<<Pref/Cref<<"  (moderate Da_M~0.2 Da_P~0.4)\n"
           <<"================================================================\n";

  FFGraph DAG;
  FFVar z=DAG.add_var("z"),rd=DAG.add_var("rd"),rwc=DAG.add_var("rwc"),rl=DAG.add_var("rl");
  FFVar Cg=DAG.add_var("Cg"),Vg=DAG.add_var("Vg"),Cd=DAG.add_var("Cd");
  FFVar Cw=DAG.add_var("Cw"),Mw=DAG.add_var("Mw"),Pw=DAG.add_var("Pw");
  FFVar Cl=DAG.add_var("Cl"),Ml=DAG.add_var("Ml"),Pl=DAG.add_var("Pl");
  FFVar lam=DAG.add_var("lambda");
  FFPartial OpP;
  auto cyl=[&](FFVar const& C,FFVar const& xi,double r0,double dr){ return OpP(C,{xi,2})+(dr/(r0+xi*dr))*OpP(C,xi); };
  auto rad=[&](FFVar const& C){ return OpP(C,{rl,2})/(Drls*Drls)+(1.0/(r2s+rl*Drls))*OpP(C,rl)/Drls; };

  FFVar GAS_C = vgIn*Vg*OpP(Cg,z) + cgf*(1.0-yin*Cg)*OpP(Cd,rd);
  FFVar GAS_V = OpP(Vg,z) + cgfV*OpP(Cd,rd);
  FFVar DRY_C = cyl(Cd,rd,r1s,Drds);
  FFVar WET_C = cyl(Cw,rwc,rws,Drws) - lam*( aMC*Cw*Mw + aPC*Cw*Pw );
  FFVar WET_M = cyl(Mw,rwc,rws,Drws) - lam*aMM*Cw*Mw;
  FFVar WET_P = cyl(Pw,rwc,rws,Drws) - lam*aPP*Cw*Pw;
  FFVar LIQ_C = vlbar*OpP(Cl,z) - DlCs*rad(Cl) + lam*( bMC*Cl*Ml + bPC*Cl*Pl );
  FFVar LIQ_M = vlbar*OpP(Ml,z) - DlMs*rad(Ml) + lam*bMM*Cl*Ml;
  FFVar LIQ_P = vlbar*OpP(Pl,z) - DlPs*rad(Pl) + lam*bPP*Cl*Pl;

  FFVar GAS_C_IN=Cg-1.0, GAS_V_IN=Vg-1.0;
  FFVar GD_VAL=Cd-Cg, DW_VAL=Cw-partH*Cd, DW_FLUX=OpP(Cd,rd)-rDWf*OpP(Cw,rwc);
  FFVar DW_MFLUX=OpP(Mw,rwc), DW_PFLUX=OpP(Pw,rwc);                    // no amine into dry pores
  FFVar WL_CVAL=Cw-Cl, WL_CFLUX=rWLf*OpP(Cw,rwc)-OpP(Cl,rl);
  FFVar WL_MVAL=Mw-Ml, WL_MFLUX=rWLf*OpP(Mw,rwc)-OpP(Ml,rl);
  FFVar WL_PVAL=Pw-Pl, WL_PFLUX=rWLf*OpP(Pw,rwc)-OpP(Pl,rl);
  FFVar L_C_WALL=OpP(Cl,rl), L_M_WALL=OpP(Ml,rl), L_P_WALL=OpP(Pl,rl);
  FFVar L_C_IN=Cl-ClInS, L_M_IN=Ml-1.0, L_P_IN=Pl-1.0;

  // geometric grading: element widths (sum=1), smallest first -> clustered at the low end
  auto graded=[](int nel,double ratio){
    double w0=(ratio-1.0)/(std::pow(ratio,nel)-1.0);
    std::vector<double> w(nel); for(int i=0;i<nel;++i) w[i]=w0*std::pow(ratio,i); return w; };

  OCFESLV oc(&DAG);
  oc.add_domain(z,  FFDom(0.,L,  3,FFDom::CGL,6));
  oc.add_domain(rd, FFDom(0.,1.0,3,FFDom::CGL,7));
  oc.add_domain(rwc,FFDom(0.0, graded(6,3.3), FFDom::CGL,8));   // two-scale ratio 3.3: ~0.12um film elem (rwc=0)
  oc.add_domain(rl, FFDom(0.0, graded(7,3.3), FFDom::CGL,8));   // two-scale ratio 3.3: ~0.12um film elem (rl=0)
  oc.add_state(Cg,{z}); oc.add_state(Vg,{z}); oc.add_state(Cd,{rd,z});
  oc.add_state(Cw,{rwc,z}); oc.add_state(Mw,{rwc,z}); oc.add_state(Pw,{rwc,z});
  oc.add_state(Cl,{rl,z}); oc.add_state(Ml,{rl,z}); oc.add_state(Pl,{rl,z});
  oc.add_input(lam,{});
  oc.update_ref(Cg,[&](OCFESLV::t_Coord const&){return 1.0;});
  oc.update_ref(Vg,[&](OCFESLV::t_Coord const&){return 1.0;});
  oc.update_ref(Cd,[&](OCFESLV::t_Coord const&){return 1.0;});
  oc.update_ref(Cw,[&](OCFESLV::t_Coord const&){return partH;});
  oc.update_ref(Mw,[&](OCFESLV::t_Coord const&){return 1.0;});
  oc.update_ref(Pw,[&](OCFESLV::t_Coord const&){return 1.0;});
  oc.update_ref(Cl,[&](OCFESLV::t_Coord const&){return ClInS;});
  oc.update_ref(Ml,[&](OCFESLV::t_Coord const&){return 1.0;});
  oc.update_ref(Pl,[&](OCFESLV::t_Coord const&){return 1.0;});

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
  oc.options.DISPLAY_LEVEL=1;
#if defined(CRONOS__WITH_SPQR)
  oc.options.SOLVE.FACTORIZATION=OCFESLV::Options::SOLVE_SPQR;
#else
  oc.options.SOLVE.FACTORIZATION=OCFESLV::Options::SOLVE_SUPERLU;
#endif

  if(!oc.setup()){ std::cerr<<"  setup() FAILED\n"; return 1; }
  check_true("setup() completed (9-state three-species)",true);
  std::vector<double> xv,inp;
  if(!oc.init(xv,inp,nullptr)){ std::cerr<<"  init() FAILED\n"; return 1; }
  double const* ip=inp.empty()?nullptr:inp.data();

  OCFESLV::SolveReport rep;
  double const lam_seq[]={0.0,0.01,0.03,0.06,0.1,0.2,0.35,0.5,0.65,0.8,1.0};
  for(double lv:lam_seq){
    oc.set_input_values(lam,{lv},inp.data());
    rep=oc.solve(xv.data(),inp.data(),nullptr);
    std::cout<<"  [continuation] lambda="<<std::fixed<<std::setprecision(2)<<lv
             <<"  converged="<<(rep.converged?"y":"n")
             <<"  |r|="<<std::scientific<<std::setprecision(3)<<rep.final_residual<<"\n";
    if(!rep.converged) break;
  }
  check_true("three-species MBC solve converged",rep.converged);
  if(!rep.converged){ std::cout<<"\n  RESULT: "<<g_pass<<" passed, "<<g_fail<<" failed\n"; return 1; }

  auto ev=[&](FFVar const& V,OCFESLV::t_Coord const& pt){ return oc.eval_colloc<double>(V,pt,xv.data(),ip,nullptr); };
  OCFESLV::t_Coord gi;gi[z]=L;   double Cgi=ev(Cg,gi);
  OCFESLV::t_Coord go;go[z]=0.0; double Cgo=ev(Cg,go);
  OCFESLV::t_Coord mi;mi[rl]=0.5;mi[z]=0.0; double Mli=ev(Ml,mi), Pli=ev(Pl,mi);
  OCFESLV::t_Coord mo;mo[rl]=0.5;mo[z]=L;   double Mlo=ev(Ml,mo), Plo=ev(Pl,mo);
  OCFESLV::t_Coord vL;vL[z]=L; double VgL=ev(Vg,vL);
  OCFESLV::t_Coord v0;v0[z]=0.0; double Vg0=ev(Vg,v0);
  double const Phi=(r1/(2.0*L))*vgIn*Cref*( VgL*Cgi - Vg0*Cgo );   // Eq 31 CO2 absorption flux [mol/m2/s]
  double const removal=(Cgi>1e-12)?(1.0-Cgo/Cgi):0.0;
  std::cout<<std::scientific<<std::setprecision(3)
           <<"\n  CO2 absorption flux Phi = "<<Phi<<" mol/m2/s   (Eq 31; paper Fig 7A ~2-3e-3)\n"
           <<std::fixed<<std::setprecision(4)
           <<"  gas CO2 : inlet "<<Cgi<<"  outlet "<<Cgo<<"  (removal "<<std::setprecision(2)<<100.0*removal<<" %)\n"
           <<std::setprecision(4)
           <<"  MDEA    : inlet "<<Mli<<"  outlet "<<Mlo<<"\n"
           <<"  PZ      : inlet "<<Pli<<"  outlet "<<Plo<<"\n";
  check_true("gas CO2 decreases (absorption)",Cgo<Cgi-1e-4);
  check_true("MDEA consumed (Ml decreases)",Mlo<Mli-1e-9);
  check_true("PZ consumed (Pl decreases)",Plo<Pli-1e-9);
  check_true("PZ consumed faster than MDEA (kP>kM)",(Pli-Plo)>(Mli-Mlo));
  check_true("all fields finite",std::isfinite(Cgo)&&std::isfinite(Mlo)&&std::isfinite(Plo));

  std::cout<<"\n============================================================\n"
           <<"  RESULT: "<<g_pass<<" passed, "<<g_fail<<" failed  -- "<<(g_fail==0?"ALL PASS":"SOME FAILED")<<"\n"
           <<"============================================================\n";
  return g_fail==0?0:1;
}
