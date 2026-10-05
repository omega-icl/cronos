// OCFE_MBC4.cpp  ---  Stage 3b, step 1: PRESCRIBED membrane wetting omega(z)
// ============================================================================================
// Adds partial membrane wetting to the Stage-2b/3 three-species MBC by letting the wetting ratio
// vary along the fiber axis:  omega(z) = omega0 + dOmega*(z/L).  The physical dry/wet radial
// extents then become z-dependent,
//
//     Drw(z) = (r2-r1)*omega(z)      (wet, liquid-filled pores)
//     Drd(z) = (r2-r1)*(1-omega(z))  (dry, gas-filled pores)
//     rw(z)  = r2 - (r2-r1)*omega(z) (dry-wet interface radius)
//
// Because membrane transport is PURELY RADIAL at each z (no d/dz in DRY_C/WET_*), rescaling in the
// vertical introduces NO moving-mesh/ALE term: the scaled radial domains rd,rwc in [0,1] and their
// collocation stay fixed, and only the radial-Laplacian positions and the reaction/flux COEFFICIENTS
// become functions of z.  For a PRESCRIBED omega(z) these are ordinary z-varying coefficients
// (evaluated per z-node) -- no new state, no cross-domain coupling; that arrives in step 2 when
// omega(z) is closed implicitly by a radial integral.
//
// VALIDATION GATE: omega(z)=const=0.20 (dOmega=0) must reproduce the constant-wetting MBC3 result
// (Phi=4.369e-3).  Then two linear profiles (dOmega=+/-0.10) must converge.
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
  std::cout<<"  "<<std::left<<std::setw(52)<<n<<"  "<<(ok?"PASS":"FAIL")<<"\n"; }

static std::vector<double> graded(int nel,double ratio){
  double w0=(ratio-1.0)/(std::pow(ratio,nel)-1.0);
  std::vector<double> w(nel); for(int i=0;i<nel;++i) w[i]=w0*std::pow(ratio,i); return w; }

struct Res { bool conv=false; double lam=0.0, resid=1e30, Phi=0.0, omgMin=0.0, omgMax=0.0; };

static Res run(double omega0, double dOmega)
{
  Res R;
  // ---- physical constants (identical to Stage-2b/3) ----------------------------------------
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
  double const kM=9.853e-2, kP=8.784e-1;                            // Stage-2b moderate rates
  double const bMC=kM*Mref, bPC=kP*Pref, bMM=kM*Cref, bPP=kP*Cref;  // liquid sinks (no wetting)
  R.omgMin=std::min(omega0,omega0+dOmega); R.omgMax=std::max(omega0,omega0+dOmega);

  // ---- DAG / model with prescribed omega(z) ------------------------------------------------
  FFGraph DAG;
  FFVar z=DAG.add_var("z"),rd=DAG.add_var("rd"),rwc=DAG.add_var("rwc"),rl=DAG.add_var("rl");
  FFVar Cg=DAG.add_var("Cg"),Vg=DAG.add_var("Vg"),Cd=DAG.add_var("Cd");
  FFVar Cw=DAG.add_var("Cw"),Mw=DAG.add_var("Mw"),Pw=DAG.add_var("Pw");
  FFVar Cl=DAG.add_var("Cl"),Ml=DAG.add_var("Ml"),Pl=DAG.add_var("Pl");
  FFVar lam=DAG.add_var("lambda");
  FFPartial OpP;

  // prescribed wetting profile (fixed expression in z -> ordinary z-varying coefficients)
  FFVar omg   = omega0 + dOmega*(z/L);
  FFVar Drw_z = (r2-r1)*omg;                 // physical wet extent   (z)
  FFVar Drd_z = (r2-r1)*(1.0-omg);           // physical dry extent   (z)
  FFVar rw_z  = r2 - (r2-r1)*omg;            // dry-wet interface radius (z)
  FFVar Drws_z= Drw_z/Rref, Drds_z=Drd_z/Rref, rws_z=rw_z/Rref;
  FFVar aMC_z =(Drw_z*Drw_z/DwC)*kM*Mref, aPC_z=(Drw_z*Drw_z/DwC)*kP*Pref;   // wet CO2 sinks
  FFVar aMM_z =(Drw_z*Drw_z/DwM)*kM*Cref, aPP_z=(Drw_z*Drw_z/DwP)*kP*Cref;   // wet amine sinks
  FFVar rDWf_z= DwC*Drd_z/(Dg*Drw_z);        // dry-wet CO2 flux ratio
  FFVar rWLf_z=(epsm/tau)*(Drl/Drw_z);       // wet-liquid flux ratio
  FFVar cgf_z =(2.0*Dg/r1)/Drd_z;            // gas-dry coupling (~1/Drd(z)); for a fixed gas-film
                                             // coefficient use instead: cgf_z = (2.0*Dg/r1)/((r2-r1)*(1.0-omega0));
  FFVar cgfV_z= cgf_z*yin/vgIn;

  // cylindrical radial Laplacian (scaled): r0 constant (dry inner radius) or FFVar (wet position).
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
  oc.options.DISPLAY_LEVEL=0; oc.options.SOLVE.VERBOSE=false;
#if defined(CRONOS__WITH_SPQR)
  oc.options.SOLVE.FACTORIZATION=OCFESLV::Options::SOLVE_SPQR;
#else
  oc.options.SOLVE.FACTORIZATION=OCFESLV::Options::SOLVE_SUPERLU;
#endif

  if(!oc.setup()){ std::cerr<<"  setup() FAILED\n"; return R; }
  std::vector<double> xv,inp;
  if(!oc.init(xv,inp,nullptr)){ std::cerr<<"  init() FAILED\n"; return R; }

  OCFESLV::SolveReport rep;
  double const lam_seq[]={0.0,0.01,0.03,0.06,0.1,0.2,0.35,0.5,0.65,0.8,1.0};
  for(double lv:lam_seq){
    oc.set_input_values(lam,{lv},inp.data());
    rep=oc.solve(xv.data(),inp.data(),nullptr);
    R.lam=lv; R.resid=rep.final_residual;
    if(!rep.converged) break;
  }
  R.conv=rep.converged;
  if(!R.conv) return R;

  double const* ip=inp.empty()?nullptr:inp.data();
  auto ev=[&](FFVar const& V,OCFESLV::t_Coord const& pt){ return oc.eval_colloc<double>(V,pt,xv.data(),ip,nullptr); };
  OCFESLV::t_Coord gi;gi[z]=L;   double Cgi=ev(Cg,gi);
  OCFESLV::t_Coord go;go[z]=0.0; double Cgo=ev(Cg,go);
  OCFESLV::t_Coord vL;vL[z]=L;   double VgL=ev(Vg,vL);
  OCFESLV::t_Coord v0;v0[z]=0.0; double Vg0=ev(Vg,v0);
  R.Phi=(r1/(2.0*L))*vgIn*Cref*( VgL*Cgi - Vg0*Cgo );
  return R;
}

int main()
{
  std::cout<<"================================================================\n"
           <<"  Stage 3b step 1: prescribed membrane wetting omega(z)=omega0+dOmega*(z/L)\n"
           <<"================================================================\n";
  double const PhiRef=4.369e-3;   // constant-wetting MBC3 reference (omega=0.20)

  struct Case{ const char* tag; double o0,dO; };
  Case cases[]={ {"const  omega=0.20",0.20, 0.00},
                 {"rising 0.20->0.30",0.20, 0.10},
                 {"falling0.20->0.10",0.20,-0.10} };

  std::cout<<"  "<<std::left<<std::setw(20)<<"profile"<<std::setw(8)<<"conv"
           <<std::setw(12)<<"|r|"<<std::setw(13)<<"Phi"<<"omg[min,max]\n";
  std::vector<Res> out;
  for(auto const& c:cases){
    Res R=run(c.o0,c.dO); out.push_back(R);
    std::string convs=R.conv?std::string("y"):("n@"+std::to_string(R.lam));
    std::cout<<"  "<<std::left<<std::setw(20)<<c.tag<<std::setw(8)<<convs
             <<std::scientific<<std::setprecision(2)<<std::setw(12)<<R.resid
             <<std::setw(13)<<R.Phi<<std::fixed<<std::setprecision(2)
             <<"["<<R.omgMin<<","<<R.omgMax<<"]\n";
  }

  std::cout<<"\n  ---- verdict ----\n";
  Res const& c0=out[0];
  double relErr = std::abs(c0.Phi-PhiRef)/PhiRef;
  std::cout<<"    const-wetting Phi="<<std::scientific<<std::setprecision(4)<<c0.Phi
           <<"  vs MBC3 ref "<<PhiRef<<"  rel.err="<<relErr<<"\n";
  check_true("omega(z)=const=0.20 reproduces MBC3 (z-coeff machinery correct)", c0.conv && relErr<2e-3);
  check_true("rising  wetting profile converges", out[1].conv);
  check_true("falling wetting profile converges", out[2].conv);
  if(out[0].conv&&out[1].conv)
    std::cout<<"    wetting sensitivity: dPhi(rising)="<<std::scientific<<std::setprecision(3)
             <<(out[1].Phi-out[0].Phi)<<"  dPhi(falling)="<<(out[2].Phi-out[0].Phi)<<"\n";

  std::cout<<"\n  RESULT: "<<g_pass<<" passed, "<<g_fail<<" failed\n";
  return g_fail?1:0;
}
