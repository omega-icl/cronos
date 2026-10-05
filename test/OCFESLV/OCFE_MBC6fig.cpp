// OCFE_MBC6.cpp  ---  NG-pilot MBC (Quek thesis, https://doi.org/10.25560/95305) with UPDATED Chapter-4 characteristics:
//   * Table 4.3 pilot module (supersedes Table 3.2 / paper Table 2)
//   * Eq. (4.9)+Table 4.4  light-hydrocarbon Henry's constants  H_i,l(C_MDEA, Tl, Pl)
//   * Eq. (4.10)+Eq. (3.37) CO2 Henry in semi-lean amine (N2O analogy)  H_CO2 = H_N2O*0.330*exp(240/Tl)
//   * Eq. (4.3) Raoult solvent evaporation in treated gas (water vapour pressure, Antoine)
// Controls (add_input is_decision=true): yin, Mgin, Flin, dP_TMPD, Tg, Tl, fCO2  (Table 3.1)
// KPI outputs (add_output @ z=0): flux(3.39), removal, loading, outlet T(3.40), HC loss(4.1-4.2), evap(4.3)
// kPZ on adaptive homotopy 'kap' (Da_P ~2e4 -> ~1e6).
// ============================================================================================
#include <iostream>
#include <iomanip>
#include <fstream>
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
  std::cout<<"  "<<std::left<<std::setw(50)<<n<<"  "<<(ok?"PASS":"FAIL")<<"\n"; }

// ---- Table 4.3 (pilot module) ----
static double const gL=2.3,gr1=225e-6,gr2=550e-6,gphi=0.38,gRm=0.115,gRt=0.08;   // fibers, packing, module, tube insert
static double const gr3=gr2/std::sqrt(gphi);
static int    const gNf=8400;
static double const gAm=68.0;                                                    // membrane area [m2]
static double const gZ=0.9,gRg=8.314,gTg=308.0,gP=54e5,grhog=49.2,grhol=986.0,gGrav=9.81;
static double const gyin=0.052,gCref=gyin*gP/(gZ*gRg*gTg),gMref=0.39*grhol/0.11916,gPref=0.05*grhol/0.08614;
// wetting closure (Table 4.3 PSD + Eq. 3.33 surrogate)
static double const gsig=0.046,gdbar=0.06e-6,ga0=9.029,ga1=-17.209,ga2=11.222,ga3=-2.938;
static double const gabscos=std::abs(std::cos(92.4*M_PI/180.0)),gomega_fix=0.20;
static double const gkP_lo=0.878,gkP_hi=40.0;
static double const gdHr=60000.0,gCp=3600.0,gCmdea=0.39,gPliq=5.43e6;            // |dHr|; Cp; MDEA mass frac; liquid P

static std::vector<double> graded(int nel,double ratio){
  double w0=(ratio-1.0)/(std::pow(ratio,nel)-1.0);
  std::vector<double> w(nel); for(int i=0;i<nel;++i) w[i]=w0*std::pow(ratio,i); return w; }

struct Vars { FFVar z,rd,rwc,rl,Cg,Vg,Cd,Cw,Mw,Pw,Cl,Ml,Pl,omg,lam,mu,kap,
              yin,Mgin,Flin,dP0,Tg,Tl,fCO2,
              kFlux,kEta,kLoad,kTout,kHC,kEvap,kMDEA,kPZ; };

static void build_mbc( FFGraph& DAG, OCFESLV& oc, Vars& V, double hydro )
{
  double const L=gL,r1=gr1,r2=gr2,eta=gphi,Rm=gRm,Rt=gRt,epsm=0.41,tau=6.1; int const nf=gNf;
  double const r3=gr3,Drl=r3-r2,Rref=r2,r1s=r1/Rref,r2s=r2/Rref,Drls=Drl/Rref;
  double const DCO2g=3.23e-7,Dg=(epsm/tau)*DCO2g;
  double const DlC=1.46e-9,DlM=2.37e-10,DlP=5.31e-10;
  double const DwC=(epsm/tau)*DlC,DwM=(epsm/tau)*DlM,DwP=(epsm/tau)*DlP;
  double const DlCs=DlC/(Rref*Rref),DlMs=DlM/(Rref*Rref),DlPs=DlP/(Rref*Rref);
  double const Z=gZ,Rg=gRg,rhog=grhog;
  double const Cref=gCref,Mref=gMref,Pref=gPref,kM=9.83e-3;
  double const bMC=kM*Mref,bMM=kM*Cref;
  double const WLf=(epsm/tau)*Drl;
  double const Dhap=-(4*r2*r2*r3*r3-r2*r2*r2*r2-3*r3*r3*r3*r3-4*r3*r3*r3*r3*std::log(r2/r3)); // sign: vl_prof>0, avg=vlbar

  V.z=DAG.add_var("z"); V.rd=DAG.add_var("rd"); V.rwc=DAG.add_var("rwc"); V.rl=DAG.add_var("rl");
  V.Cg=DAG.add_var("Cg"); V.Vg=DAG.add_var("Vg"); V.Cd=DAG.add_var("Cd");
  V.Cw=DAG.add_var("Cw"); V.Mw=DAG.add_var("Mw"); V.Pw=DAG.add_var("Pw");
  V.Cl=DAG.add_var("Cl"); V.Ml=DAG.add_var("Ml"); V.Pl=DAG.add_var("Pl");
  V.omg=DAG.add_var("omega");
  V.lam=DAG.add_var("lambda"); V.mu=DAG.add_var("mu"); V.kap=DAG.add_var("kappa");
  V.yin=DAG.add_var("yin"); V.Mgin=DAG.add_var("Mgin"); V.Flin=DAG.add_var("Flin");
  V.dP0=DAG.add_var("dPtm"); V.Tg=DAG.add_var("Tg"); V.Tl=DAG.add_var("Tl"); V.fCO2=DAG.add_var("fCO2");
  FFVar &z=V.z,&rd=V.rd,&rwc=V.rwc,&rl=V.rl,&Cg=V.Cg,&Vg=V.Vg,&Cd=V.Cd,&Cw=V.Cw,&Mw=V.Mw,&Pw=V.Pw;
  FFVar &Cl=V.Cl,&Ml=V.Ml,&Pl=V.Pl,&omg=V.omg,&lam=V.lam,&mu=V.mu,&kap=V.kap;
  FFPartial OpP;

  FFVar kPZ=gkP_lo+kap*(gkP_hi-gkP_lo);
  FFVar vgIn=V.Mgin/(nf*M_PI*r1*r1*rhog);                        // Eq. 3.14
  FFVar vlbar=V.Flin/(M_PI*(Rm*Rm-Rt*Rt)*(1.0-eta));            // Eq. 3.6 with tube insert (Table 4.3)
  FFVar dP_z=V.dP0+hydro*(L-z);
  // CO2 Henry in semi-lean amine: Eq. (4.10) H_N2O, then N2O->CO2 analogy (Eq. 3.37, Versteeg)
  FFVar HN2O=-3.30e4-3.79e4*V.fCO2+1.70e2*V.Tl-4.37e3*V.fCO2*V.fCO2-1.34e-1*V.Tl*V.Tl+1.45e2*V.fCO2*V.Tl;
  FFVar HCO2=HN2O*0.3298*exp(240.0/V.Tl);
  FFVar partH=Z*Rg*V.Tl/HCO2;
  FFVar r_l=r2+rl*(r3-r2);
  FFVar vl_prof=2.0*vlbar*(r3*r3-r2*r2)/Dhap*(r_l*r_l-r2*r2+2.0*r3*r3*log(r2/r_l));   // Eq. 3.17

  FFVar delta_w=(2.0*gsig*gabscos)/dP_z, uu=delta_w/gdbar;
  FFVar arg=ga0+ga1*uu+ga2*uu*uu+ga3*uu*uu*uu, Zeta_z=1.0/(1.0+exp(-2.0*arg));
  FFVar target=(1.0-mu)*gomega_fix+mu*Zeta_z;

  FFVar Drw_z=(r2-r1)*omg, Drd_z=(r2-r1)*(1.0-omg), rw_z=r2-(r2-r1)*omg;
  FFVar Drws_z=Drw_z/Rref, Drds_z=Drd_z/Rref, rws_z=rw_z/Rref;
  FFVar aMC_z=(Drw_z*Drw_z/DwC)*bMC, aPC_z=(Drw_z*Drw_z/DwC)*kPZ*Pref;
  FFVar aMM_z=(Drw_z*Drw_z/DwM)*kM*Cref, aPP_z=(Drw_z*Drw_z/DwP)*kPZ*Cref;
  FFVar cgf_z=(2.0*Dg/r1)/Drd_z, cgfV_z=cgf_z*V.yin/vgIn;
  FFVar DWf_z=DwC*Drd_z/Dg;
  FFVar bPC=kPZ*Pref, bPP=kPZ*Cref;

  auto cyl_d=[&](FFVar const& C,FFVar const& xi,double r0,FFVar const& dr){ return OpP(C,{xi,2})+(dr/(r0+xi*dr))*OpP(C,xi); };
  auto cyl_v=[&](FFVar const& C,FFVar const& xi,FFVar const& r0,FFVar const& dr){ return OpP(C,{xi,2})+(dr/(r0+xi*dr))*OpP(C,xi); };
  auto rad  =[&](FFVar const& C){ return OpP(C,{rl,2})/(Drls*Drls)+(1.0/(r2s+rl*Drls))*OpP(C,rl)/Drls; };

  FFVar GAS_C=vgIn*Vg*OpP(Cg,z)+cgf_z*(1.0-V.yin*Cg)*OpP(Cd,rd);
  FFVar GAS_V=OpP(Vg,z)+cgfV_z*OpP(Cd,rd);
  FFVar DRY_C=cyl_d(Cd,rd,r1s,Drds_z);
  FFVar WET_C=cyl_v(Cw,rwc,rws_z,Drws_z)-lam*(aMC_z*Cw*Mw+aPC_z*Cw*Pw);
  FFVar WET_M=cyl_v(Mw,rwc,rws_z,Drws_z)-lam*aMM_z*Cw*Mw;
  FFVar WET_P=cyl_v(Pw,rwc,rws_z,Drws_z)-lam*aPP_z*Cw*Pw;
  FFVar LIQ_C=vl_prof*OpP(Cl,z)-DlCs*rad(Cl)+lam*(bMC*Cl*Ml+bPC*Cl*Pl);
  FFVar LIQ_M=vl_prof*OpP(Ml,z)-DlMs*rad(Ml)+lam*bMM*Cl*Ml;
  FFVar LIQ_P=vl_prof*OpP(Pl,z)-DlPs*rad(Pl)+lam*bPP*Cl*Pl;
  FFVar GAS_C_IN=Cg-1.0,GAS_V_IN=Vg-1.0;
  FFVar GD_VAL=Cd-Cg,DW_VAL=Cw-partH*Cd;
  FFVar DW_FLUX=Drw_z*OpP(Cd,rd)-DWf_z*OpP(Cw,rwc);
  FFVar DW_MFLUX=OpP(Mw,rwc),DW_PFLUX=OpP(Pw,rwc);
  FFVar WL_CVAL=Cw-Cl,WL_CFLUX=WLf*OpP(Cw,rwc)-Drw_z*OpP(Cl,rl);
  FFVar WL_MVAL=Mw-Ml,WL_MFLUX=WLf*OpP(Mw,rwc)-Drw_z*OpP(Ml,rl);
  FFVar WL_PVAL=Pw-Pl,WL_PFLUX=WLf*OpP(Pw,rwc)-Drw_z*OpP(Pl,rl);
  FFVar L_C_WALL=OpP(Cl,rl),L_M_WALL=OpP(Ml,rl),L_P_WALL=OpP(Pl,rl);
  FFVar L_C_IN=Cl-0.0,L_M_IN=Ml-1.0,L_P_IN=Pl-1.0;
  FFVar OMG_CL=omg-target;

  // ---- KPI output functions (evaluated at gas outlet z=0) ----
  FFVar co2rem=(nf*M_PI*r1*r1*Cref)*vgIn*(1.0-Vg*Cg);            // CO2 removal rate [mol/s]  (Eq. 3.39 * Am)
  V.kFlux=co2rem/gAm;                                           // (1) absorption flux [mol/m2/s]
  V.kEta =1.0-Vg*Cg;                                            // (2) removal efficiency [-]
  V.kLoad=V.fCO2+co2rem/((Mref+Pref)*V.Flin);                   // (3) CO2 loading, solvent outlet [mol/mol]
  V.kTout=V.Tl+co2rem*gdHr/(grhol*V.Flin*gCp);                  // (4) solvent outlet T [K]  (Eq. 3.40)
  // (5) light-HC (CH4) loss: Eq. (4.9)+Table 4.4 Henry, then Eqs. (4.1-4.2)
  FFVar H_CH4=-2.90e5+3.56e5*gCmdea+1.26e3*V.Tl+2.40e-2*gPliq
              -1.42e3*gCmdea*V.Tl+4.41e-3*gCmdea*gPliq-7.43e-5*V.Tl*gPliq;
  V.kHC=((1.0-V.yin)*gP/H_CH4)*V.Flin;
  // (6-8) solvent evaporation into treated gas: Raoult (Eq. 4.3)  y_i = x_i P_i^vap / Pg,  i in {H2O,MDEA,PZ}
  //   x_i from 39/5 wt% MDEA/PZ (rest water); vapour pressures via Antoine (T in Celsius).
  //   vapour-pressure Antoine coeffs fitted to the extracted Fig 4.13 data (P[Pa]=exp(A-B/(T-C)), T[K]).
  double const xH2O=0.8897, xMDEA=0.0937, xPZ=0.0166;
  FFVar Pv_H2O =exp(18.8981-1981.33/(V.Tl-114.50));           // water [Pa]  Antoine fit to Fig 4.13 (T in K)
  FFVar Pv_MDEA=exp(21.5532-4391.20/(V.Tl-95.00));            // MDEA  [Pa]  Antoine fit to Fig 4.13
  FFVar Pv_PZ  =exp(26.0538-6348.82/(V.Tl+23.50));            // PZ    [Pa]  Antoine fit to Fig 4.13
  FFVar molflow=(nf*M_PI*r1*r1*gP/(gZ*gRg*V.Tg))*vgIn*Vg;               // treated-gas molar flow [mol/s]
  V.kEvap=(xH2O *Pv_H2O /gP)*molflow;                          // (6) water evaporation [mol/s]
  V.kMDEA=(xMDEA*Pv_MDEA/gP)*molflow;                          // (7) MDEA makeup (evaporation) [mol/s]
  V.kPZ  =(xPZ  *Pv_PZ  /gP)*molflow;                          // (8) PZ makeup (evaporation) [mol/s]

  oc.add_domain(z,FFDom(0.,L,3,FFDom::CGL,6));
  oc.add_domain(rd,FFDom(0.,1.0,3,FFDom::CGL,7));
  oc.add_domain(rwc,FFDom(0.0,graded(6,3.3),FFDom::CGL,8));
  oc.add_domain(rl,FFDom(0.0,graded(7,3.3),FFDom::CGL,8));
  oc.add_state(Cg,{z}); oc.add_state(Vg,{z}); oc.add_state(Cd,{rd,z});
  oc.add_state(Cw,{rwc,z}); oc.add_state(Mw,{rwc,z}); oc.add_state(Pw,{rwc,z});
  oc.add_state(Cl,{rl,z}); oc.add_state(Ml,{rl,z}); oc.add_state(Pl,{rl,z});
  oc.add_state(omg,{z});
  oc.add_input(lam,{}); oc.add_input(mu,{}); oc.add_input(kap,{});
  oc.add_input(V.yin, gyin,          true);                     // Table 3.1 controls
  oc.add_input(V.Mgin,75.0/3600.0,   true);
  oc.add_input(V.Flin,220e-3/3600.0, true);
  oc.add_input(V.dP0, 30.0e3,        true);
  oc.add_input(V.Tg,  308.0,         true);
  oc.add_input(V.Tl,  308.0,         true);
  oc.add_input(V.fCO2,0.01,          true);
  oc.update_ref(Cg,[&](OCFESLV::t_Coord const&){return 1.0;});
  oc.update_ref(Vg,[&](OCFESLV::t_Coord const&){return 1.0;});
  oc.update_ref(Cd,[&](OCFESLV::t_Coord const&){return 1.0;});
  oc.update_ref(Cw,[&](OCFESLV::t_Coord const&){return 0.478;});   // nominal partH (Eq. 4.10 @ Tl=308,f=0.01)
  oc.update_ref(Mw,[&](OCFESLV::t_Coord const&){return 1.0;});
  oc.update_ref(Pw,[&](OCFESLV::t_Coord const&){return 1.0;});
  oc.update_ref(Cl,[&](OCFESLV::t_Coord const&){return 0.0;});
  oc.update_ref(Ml,[&](OCFESLV::t_Coord const&){return 1.0;});
  oc.update_ref(Pl,[&](OCFESLV::t_Coord const&){return 1.0;});
  oc.update_ref(omg,[&](OCFESLV::t_Coord const&){return gomega_fix;});

  int const Z_NO_UB=FFDom::ALL-FFDom::UB, Z_NO_LB=FFDom::ALL-FFDom::LB, R_INT=FFDom::ALL-FFDom::LB-FFDom::UB;
  OCFESLV::EqnOptions blk(OCFESLV::EqnRole::INTERIOR,0),bnd(OCFESLV::EqnRole::BOUNDARY,0),itf(OCFESLV::EqnRole::INTERFACE,0);
  oc.add_equation(GAS_V,{z,rd},{Z_NO_UB,FFDom::LB},blk);
  oc.add_equation(GAS_C,{z,rd},{Z_NO_UB,FFDom::LB},blk);
  oc.add_equation(GAS_V_IN,{z},{FFDom::UB},bnd);
  oc.add_equation(GAS_C_IN,{z},{FFDom::UB},bnd);
  oc.add_equation(DRY_C,{rd,z},{R_INT,FFDom::ALL},blk);
  oc.add_equation(WET_C,{rwc,z},{R_INT,FFDom::ALL},blk);
  oc.add_equation(WET_M,{rwc,z},{R_INT,FFDom::ALL},blk);
  oc.add_equation(WET_P,{rwc,z},{R_INT,FFDom::ALL},blk);
  oc.add_equation(LIQ_C,{rl,z},{R_INT,Z_NO_LB},blk);
  oc.add_equation(LIQ_M,{rl,z},{R_INT,Z_NO_LB},blk);
  oc.add_equation(LIQ_P,{rl,z},{R_INT,Z_NO_LB},blk);
  oc.add_equation(OMG_CL,{z},{FFDom::ALL},blk);
  oc.add_equation(GD_VAL,{rd,z},{FFDom::LB,FFDom::ALL},itf);
  oc.add_equation(DW_VAL,{rd,rwc,z},{FFDom::UB,FFDom::LB,FFDom::ALL},itf);
  oc.add_equation(DW_FLUX,{rd,rwc,z},{FFDom::UB,FFDom::LB,FFDom::ALL},itf);
  oc.add_equation(DW_MFLUX,{rwc,z},{FFDom::LB,FFDom::ALL},itf);
  oc.add_equation(DW_PFLUX,{rwc,z},{FFDom::LB,FFDom::ALL},itf);
  oc.add_equation(WL_CVAL,{rwc,rl,z},{FFDom::UB,FFDom::LB,FFDom::ALL},itf);
  oc.add_equation(WL_CFLUX,{rwc,rl,z},{FFDom::UB,FFDom::LB,FFDom::ALL},itf);
  oc.add_equation(WL_MVAL,{rwc,rl,z},{FFDom::UB,FFDom::LB,FFDom::ALL},itf);
  oc.add_equation(WL_MFLUX,{rwc,rl,z},{FFDom::UB,FFDom::LB,FFDom::ALL},itf);
  oc.add_equation(WL_PVAL,{rwc,rl,z},{FFDom::UB,FFDom::LB,FFDom::ALL},itf);
  oc.add_equation(WL_PFLUX,{rwc,rl,z},{FFDom::UB,FFDom::LB,FFDom::ALL},itf);
  oc.add_equation(L_C_WALL,{rl,z},{FFDom::UB,FFDom::ALL},bnd);
  oc.add_equation(L_M_WALL,{rl,z},{FFDom::UB,FFDom::ALL},bnd);
  oc.add_equation(L_P_WALL,{rl,z},{FFDom::UB,FFDom::ALL},bnd);
  oc.add_equation(L_C_IN,{rl,z},{R_INT,FFDom::LB},bnd);
  oc.add_equation(L_M_IN,{rl,z},{R_INT,FFDom::LB},bnd);
  oc.add_equation(L_P_IN,{rl,z},{R_INT,FFDom::LB},bnd);

  oc.add_output(V.kFlux,{z},{0.0});
  oc.add_output(V.kEta, {z},{0.0});
  oc.add_output(V.kLoad,{z},{0.0});
  oc.add_output(V.kTout,{z},{0.0});
  oc.add_output(V.kHC,  {z},{0.0});
  oc.add_output(V.kEvap,{z},{0.0});
  oc.add_output(V.kMDEA,{z},{0.0});
  oc.add_output(V.kPZ,  {z},{0.0});

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

static bool adaptive_ramp( OCFESLV& oc, std::vector<double>& xv, std::vector<double>& inp,
                           FFVar const& param, double start, double target, double& resid )
{
  OCFESLV::SolveReport rep; std::vector<double> saved=xv;
  double cur=start, d=0.1, dir=(target>=start?1.0:-1.0);
  while( std::abs(cur-target) > 1e-9 ){
    double trial = cur + dir*std::min(d, std::abs(target-cur));
    oc.set_input_values(param,{trial},inp.data());
    rep=oc.solve(xv.data(),inp.data(),nullptr); resid=rep.final_residual;
    if(rep.converged){ cur=trial; saved=xv; d=std::min(0.1,d*1.4); }
    else { xv=saved; d*=0.5; if(d<1e-5) return false; }
  }
  return true;
}

struct KPI { bool conv=false; double resid=1e30,flux=0,eta=0,load=0,Tout=0,hc=0,evap=0,mdea=0,pz=0; };

static void write_profiles( OCFESLV& oc, Vars const& V, double const* xv, double const* inp, std::string const& tag )
{
  std::ofstream fC("fig10_"+tag+"_CO2.dat"),fM("fig10_"+tag+"_MDEA.dat"),fP("fig10_"+tag+"_PZ.dat");
  fC<<std::scientific<<std::setprecision(6); fM<<std::scientific<<std::setprecision(6); fP<<std::scientific<<std::setprecision(6);
  auto ev=[&](FFVar const& X,OCFESLV::t_Coord const& pt){ return oc.eval_colloc<double>(X,pt,xv,inp,nullptr); };
  int const NR=80,NZ=80;
  for(int jz=0;jz<=NZ;++jz){ double zz=gL*jz/NZ; OCFESLV::t_Coord pz; pz[V.z]=zz;
    double rw=gr2-(gr2-gr1)*ev(V.omg,pz);
    for(int ir=0;ir<=NR;++ir){ double r=gr3*ir/NR; double co2=0,mdea=0,pz_=0;
      if(r<gr1){ co2=ev(V.Cg,pz)*gCref; }
      else if(r<rw){ OCFESLV::t_Coord pt=pz; pt[V.rd]=(r-gr1)/(rw-gr1); co2=ev(V.Cd,pt)*gCref; }
      else if(r<gr2){ OCFESLV::t_Coord pt=pz; pt[V.rwc]=(r-rw)/(gr2-rw);
                      co2=ev(V.Cw,pt)*gCref; mdea=ev(V.Mw,pt)*gMref; pz_=ev(V.Pw,pt)*gPref; }
      else { OCFESLV::t_Coord pt=pz; pt[V.rl]=(r-gr2)/(gr3-gr2);
             co2=ev(V.Cl,pt)*gCref; mdea=ev(V.Ml,pt)*gMref; pz_=ev(V.Pl,pt)*gPref; }
      fC<<r*1e6<<" "<<zz<<" "<<co2<<"\n"; fM<<r*1e6<<" "<<zz<<" "<<mdea<<"\n"; fP<<r*1e6<<" "<<zz<<" "<<pz_<<"\n"; }
    fC<<"\n"; fM<<"\n"; fP<<"\n"; }
  std::cout<<"    wrote fig10_"<<tag<<"_{CO2,MDEA,PZ}.dat\n";
}

static void write_gnuplot()
{
  std::ofstream g("fig10.gp");
  g<<"set terminal pngcairo size 1500,950 font 'Helvetica,11'\nset output 'fig10.png'\n"
   <<"set pm3d map interpolate 2,2\nset palette defined (0 'white', 0.5 '#4fa3d1', 1 '#08306b')\n"
   <<"r1="<<gr1*1e6<<"; r2="<<gr2*1e6<<"; r3="<<gr3*1e6<<"\n"
   <<"set xlabel 'radial position r [{/Symbol m}m]'\nset ylabel 'axial position z [m]'\n"
   <<"set xrange [0:"<<gr3*1e6<<"]\nset yrange [0:"<<gL<<"]\nset xtics 200\nunset key\n"
   <<"set arrow 1 from r1,0 to r1,"<<gL<<" nohead lc rgb 'red' lw 1 front\n"
   <<"set arrow 2 from r2,0 to r2,"<<gL<<" nohead lc rgb 'red' lw 1 front\n"
   <<"set arrow 3 from r3,0 to r3,"<<gL<<" nohead lc rgb 'red' lw 1 front\n"
   <<"set multiplot layout 2,3 title 'MBC single-fibre concentration fields (mol m^{-3}) -- corrected pilot'\n"
   <<"set cbrange [0:130];  set title 'CO_2 (vertical)';   splot 'fig10_vert_CO2.dat'  u 1:2:3\n"
   <<"set cbrange [0:3500]; set title 'MDEA (vertical)';   splot 'fig10_vert_MDEA.dat' u 1:2:3\n"
   <<"set cbrange [0:600];  set title 'PZ (vertical)';     splot 'fig10_vert_PZ.dat'   u 1:2:3\n"
   <<"set cbrange [0:130];  set title 'CO_2 (horizontal)'; splot 'fig10_horiz_CO2.dat' u 1:2:3\n"
   <<"set cbrange [0:3500]; set title 'MDEA (horizontal)'; splot 'fig10_horiz_MDEA.dat'u 1:2:3\n"
   <<"set cbrange [0:600];  set title 'PZ (horizontal)';   splot 'fig10_horiz_PZ.dat'  u 1:2:3\n"
   <<"unset multiplot\n";
  std::cout<<"  wrote fig10.gp\n";
}

static KPI run(double hydro,std::string const& tag)
{
  KPI K;
  FFGraph DAG; OCFESLV oc(&DAG); Vars V;
  build_mbc(DAG,oc,V,hydro);
  std::vector<double> xv,inp; if(!oc.init(xv,inp,nullptr)){ std::cerr<<"  ["<<tag<<"] init failed\n"; return K; }
  oc.set_input_values(V.lam,{0.0},inp.data());
  oc.set_input_values(V.mu,{0.0},inp.data());
  oc.set_input_values(V.kap,{0.0},inp.data());
  if(!adaptive_ramp(oc,xv,inp,V.lam,0.0,1.0,K.resid)){ std::cerr<<"  ["<<tag<<"] lam ramp failed\n"; return K; }
  if(!adaptive_ramp(oc,xv,inp,V.mu, 0.0,1.0,K.resid)){ std::cerr<<"  ["<<tag<<"] mu ramp failed\n"; return K; }
  if(!adaptive_ramp(oc,xv,inp,V.kap,0.0,1.0,K.resid)){ std::cerr<<"  ["<<tag<<"] kap ramp failed\n"; return K; }
  K.conv=true;
  double const* ip=inp.data();
  OCFESLV::t_Coord p0; p0[V.z]=0.0;
  // eval_colloc interpolates STATES (outputs are for FFOCFESLV; monolithic val_functions() is empty)
  double Vg0=oc.eval_colloc<double>(V.Vg,p0,xv.data(),ip,nullptr);
  double Cg0=oc.eval_colloc<double>(V.Cg,p0,xv.data(),ip,nullptr);
  // KPIs from outlet states + nominal controls (Eqs. 3.39/3.40, 4.9, 4.1-4.3)
  double const yinv=gyin,Mginv=75.0/3600.0,Flinv=220e-3/3600.0,Tlv=308.0,Tgv=308.0,fCO2v=0.01;
  double vgIn=Mginv/(gNf*M_PI*gr1*gr1*grhog);
  double co2rem=gCref*(Mginv/grhog)*(1.0-Vg0*Cg0);                     // CO2 removal rate [mol/s] (Eq. 3.39*Am)
  K.flux=co2rem/gAm;
  K.eta =1.0-Vg0*Cg0;
  K.load=fCO2v+co2rem/((gMref+gPref)*Flinv);
  K.Tout=Tlv+co2rem*gdHr/(grhol*Flinv*gCp);                            // Eq. 3.40
  double H_CH4=-2.90e5+3.56e5*gCmdea+1.26e3*Tlv+2.40e-2*gPliq-1.42e3*gCmdea*Tlv+4.41e-3*gCmdea*gPliq-7.43e-5*Tlv*gPliq;
  K.hc=((1.0-yinv)*gP/H_CH4)*Flinv;                                    // Eqs. 4.9 / 4.1-4.2
  double PvH2O =std::exp(18.8981-1981.33/(Tlv-114.50));               // Raoult evap (Eq. 4.3); Antoine fits to Fig 4.13 [Pa]
  double PvMDEA=std::exp(21.5532-4391.20/(Tlv-95.00));
  double PvPZ  =std::exp(26.0538-6348.82/(Tlv+23.50));
  double molflow=(gNf*M_PI*gr1*gr1*gP/(gZ*gRg*Tgv))*vgIn*Vg0;
  K.evap=(0.8897*PvH2O /gP)*molflow;
  K.mdea=(0.0937*PvMDEA/gP)*molflow;
  K.pz  =(0.0166*PvPZ  /gP)*molflow;
  write_profiles(oc,V,xv.data(),inp.data(),tag);
  return K;
}

int main()
{
  std::cout<<"================================================================\n"
           <<"  MBC6 test : Table 4.3 module, Eq 4.9/4.10 Henry, KPI outputs\n"
           <<"================================================================\n";
  KPI Vt=run(grhol*gGrav,"vert");
  KPI Hz=run(0.0,        "horiz");
  if(Vt.conv&&Hz.conv) write_gnuplot();

  std::cout<<"\n  "<<std::left<<std::setw(34)<<"KPI"<<std::setw(15)<<"vertical"<<std::setw(15)<<"horizontal"<<"\n";
  auto row=[&](char const* nm,double a,double b,char const* fmt){
    char ba[32],bb[32]; std::snprintf(ba,32,fmt,a); std::snprintf(bb,32,fmt,b);
    std::cout<<"  "<<std::left<<std::setw(34)<<nm<<std::setw(15)<<ba<<std::setw(15)<<bb<<"\n"; };
  if(Vt.conv&&Hz.conv){
    row("CO2 absorption flux [mol/m2/s]", Vt.flux,Hz.flux,"%.3e");
    row("CO2 removal efficiency [%]",     100*Vt.eta,100*Hz.eta,"%.1f");
    row("CO2 loading, outlet [mol/mol]",  Vt.load,Hz.load,"%.4f");
    row("solvent outlet T [K]",           Vt.Tout,Hz.Tout,"%.2f");
    row("light-HC (CH4) loss [mol/s]",    Vt.hc,Hz.hc,"%.3e");
    row("solvent (H2O) evap [mol/s]",     Vt.evap,Hz.evap,"%.3e");
    row("MDEA makeup (evap) [mol/s]",      Vt.mdea,Hz.mdea,"%.3e");
    row("PZ makeup (evap) [mol/s]",        Vt.pz,Hz.pz,"%.3e");
  }
  std::cout<<"\n  ---- verdict ----\n";
  check_true("vertical converges (full kPZ=40)",   Vt.conv);
  check_true("horizontal converges (full kPZ=40)", Hz.conv);
  check_true("removal efficiencies in (0,1]",      Vt.conv&&Hz.conv&&Vt.eta>0&&Vt.eta<=1.0&&Hz.eta>0&&Hz.eta<=1.0);
  check_true("outlet T rise positive",             Vt.conv&&Hz.conv&&Vt.Tout>=308.0&&Hz.Tout>=308.0);
  std::cout<<"\n  RESULT: "<<g_pass<<" passed, "<<g_fail<<" failed\n";
  return g_fail?1:0;
}
