// OCFE_MBC5.cpp  ---  NG-pilot MBC (Quek thesis, https://doi.org/10.25560/95305), CHAPTER-3 module + latest correlations:
//   * Table 3.2 Chapter-3 pilot module: r1=431, r2=846, L=2.0, phi=0.55, N=10007, Am=106.5, NO tube insert
//   * Chapter-3 wetting surrogate (Eq. 3.33): a0=7.966,a1=-14.08,a2=8.418,a3=-2.041; dbar=0.08um; theta=92.5deg
//     (this module + wetting is the basis for the vertical/horizontal orientation study, Fig. 3.11)
//   * Eq. (4.9)+Table 4.4  light-hydrocarbon Henry's constants  H_i,l(C_MDEA, Tl, Pl)
//   * CO2 Henry: CHAPTER-3 N2O analogy (Eqs 3.34-3.38, Table A.1 pilot NG -> H_CO2,l~3588, partH~0.64)
//     -- REPLACES the Chapter-4 Eq (4.10) polynomial that MBC5 previously carried (partH 0.478).
//   * Eq. (4.3) Raoult solvent evaporation in treated gas (water vapour pressure, Antoine)
//   * TEMPERATURE CORRECTION (Eq 3.40, Fig 3.11): liquid temperature Tlo(z) is a DISTRIBUTED-CONSTANT
//     state on {z} (flat transport dTlo/dz=0), closing the lumped adiabatic energy balance at the z=0
//     boundary:  rho_l*Fl*Cp*(Tlo-Tl_in) = Phi*Am*|dHr|  with Phi from the gas-outlet flux (1-Vg(0)Cg(0)).
//     (A true 0-dim scalar state has no collocation node -> rectangular symbol; Tlo(z) lives like Cg/Vg.)
//     Tlo feeds sigma(T) (0.046 -> 0.044 over dT~12 K) -> extra wetting.  Homotopy 'tauT' ramps the
//     sigma feedback 0->1 (tauT=0 reproduces the isothermal baseline; Henry stays at inlet-Tl, sigma-only).
// Controls (add_input is_decision=true): yin, Mgin, Flin, dP_TMPD, Tg, Tl, fCO2  (Table 3.1)
// KPI outputs (add_output @ z=0): flux(3.39), removal, loading, outlet T(3.40), HC loss(4.1-4.2), evap(4.3)
// kPZ on adaptive homotopy 'kap' (Da_P ~2e4 -> ~1e6); sigma(T) feedback on homotopy 'tauT'.
// ============================================================================================
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <string>
#include <chrono>

#ifndef OCFE_OCFESLV_HEADER
#define OCFE_OCFESLV_HEADER "ocfeslv.hpp"
#endif
#include OCFE_OCFESLV_HEADER

using namespace mc;

static int g_pass=0,g_fail=0;
static void check_true(char const* n,bool ok){ (ok?g_pass:g_fail)++;
  std::cout<<"  "<<std::left<<std::setw(50)<<n<<"  "<<(ok?"PASS":"FAIL")<<"\n"; }

// ---- Table 4.3 (pilot module) ----
static double const gL=2.0,gr1=431e-6,gr2=846e-6,gphi=0.55,gRm=0.115;   // Chapter-3 pilot module (Table 3.2), no tube insert
static double const gr3=gr2/std::sqrt(gphi);
static int    const gNf=10007;
static double const gAm=106.5;                                                    // membrane area [m2]
static double const gZ=0.9,gRg=8.314,gTg=308.0,gP=54e5,grhol=986.0,gGrav=9.81;
static double const gyin=0.052,gCref=gyin*gP/(gZ*gRg*gTg),gMref=0.39*grhol/0.11916,gPref=0.05*grhol/0.08614;
// Real-gas inlet density computed self-consistently as rho_g = Mavg*P/(Z R T), rather than the
// Table A.1 value 49.2 kg/m3 (which is inconsistent with Z=0.9 unless Mavg~21; PR gives Z~0.90 and
// rho_g~41 for the CH4-dominated pilot NG).  Set gMhc to the actual pilot NG molar mass.
static double const gMco2=44.01, gMhc=16.04;                          // CO2 and NG-hydrocarbon molar masses [g/mol] (16.04 = pure CH4)
static double const gMavg=gyin*gMco2+(1.0-gyin)*gMhc;                 // inlet mean molar mass [g/mol]
static double const grhog=1.0e-3*gMavg*gP/(gZ*gRg*gTg);              // -> ~41.0 kg/m3 (CH4-dominated); was 49.2
// wetting closure (Table 4.3 PSD + Eq. 3.33 surrogate)
static double const gsig=0.046,gdbar=0.08e-6,ga0=7.966,ga1=-14.08,ga2=8.418,ga3=-2.041;
static double const gTlref=308.0, gsig_dT=-1.6667e-4;   // sigma(T)=gsig+gsig_dT*(Tl-gTlref): 0.046->0.044 over dT~12 K (Sec 3.4.3)
static double const gabscos=std::abs(std::cos(92.5*M_PI/180.0)),gomega_fix=0.20;
static double const gkP_lo=0.878,gkP_hi=40.0;
static double const gdHr=60000.0,gCp=3600.0,gCmdea=0.39,gPliq=5.43e6;            // |dHr|; Cp; MDEA mass frac; liquid P

static std::vector<double> graded(int nel,double ratio){
  double w0=(ratio-1.0)/(std::pow(ratio,nel)-1.0);
  std::vector<double> w(nel); for(int i=0;i<nel;++i) w[i]=w0*std::pow(ratio,i); return w; }

struct Vars { FFVar z,rd,rwc,rl,Cg,Vg,Cd,Cw,Mw,Pw,Cl,Ml,Pl,omg,Tlo,lam,mu,kap,tauT,
              yin,Mgin,Flin,dP0,Tg,Tl,fCO2,
              kFlux,kEta,kLoad,kTout,kHC,kEvap,kMDEA,kPZ; };

static void build_mbc( FFGraph& DAG, OCFESLV& oc, Vars& V, double hydro )
{
  double const L=gL,r1=gr1,r2=gr2,eta=gphi,Rm=gRm,epsm=0.41,tau=6.1; int const nf=gNf;
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
  V.omg=DAG.add_var("omega"); V.Tlo=DAG.add_var("Tlo");
  V.lam=DAG.add_var("lambda"); V.mu=DAG.add_var("mu"); V.kap=DAG.add_var("kappa"); V.tauT=DAG.add_var("tauT");
  V.yin=DAG.add_var("yin"); V.Mgin=DAG.add_var("Mgin"); V.Flin=DAG.add_var("Flin");
  V.dP0=DAG.add_var("dPtm"); V.Tg=DAG.add_var("Tg"); V.Tl=DAG.add_var("Tl"); V.fCO2=DAG.add_var("fCO2");
  FFVar &z=V.z,&rd=V.rd,&rwc=V.rwc,&rl=V.rl,&Cg=V.Cg,&Vg=V.Vg,&Cd=V.Cd,&Cw=V.Cw,&Mw=V.Mw,&Pw=V.Pw;
  FFVar &Cl=V.Cl,&Ml=V.Ml,&Pl=V.Pl,&omg=V.omg,&Tlo=V.Tlo,&lam=V.lam,&mu=V.mu,&kap=V.kap,&tauT=V.tauT;
  FFPartial OpP;

  FFVar kPZ=gkP_lo+kap*(gkP_hi-gkP_lo);
  FFVar vgIn=V.Mgin/(nf*M_PI*r1*r1*rhog);                        // Eq. 3.14
  FFVar vlbar=V.Flin/(M_PI*Rm*Rm*(1.0-eta));                    // Eq. 3.6 (no tube insert, Chapter 3)
  FFVar dP_z=V.dP0+hydro*(L-z);
  // CO2 Henry in lean amine -- CHAPTER-3 N2O analogy (Eqs 3.34-3.38), evaluated at INLET Tl
  // (sigma-only temperature correction: solubility held at the inlet-T value, per Sec 3.4.3).
  // Volume fractions of MDEA/PZ/H2O in 39/5 wt% aq. solution (rho_MDEA=1042, rho_PZ=1100, rho_H2O=994):
  double const phiM=0.381, phiP=0.046, phiW=0.573;                       // sum=1; Eq 3.37 uses volume fractions
  FFVar Hco2w = 2.82e6*exp(-2044.0/V.Tl);                                // Eq 3.35  H_CO2,H2O
  FFVar Hn2ow = 8.55e6*exp(-2284.0/V.Tl);                                // Eq 3.36  H_N2O,H2O
  FFVar Hn2oa = 1.52e5*exp(-1312.7/V.Tl);                                // Eq 3.38  H_N2O,MDEA = H_N2O,PZ
  FFVar lnHn2ol = (phiM+phiP)*log(Hn2oa) + phiW*log(Hn2ow)
                + phiM*phiW*(-2.899 + 1405.43/V.Tl);                     // Eq 3.37  ln H_N2O,l
  FFVar HCO2 = exp(lnHn2ol)*(Hco2w/Hn2ow);                               // Eq 3.34  H_CO2,l (~3666 @308K; Table A.1: 3588)
  FFVar partH=Z*Rg*V.Tl/HCO2;                                            // partition C_wet/C_dry (~0.63; Table A.1: 0.642)
  FFVar r_l=r2+rl*(r3-r2);
  FFVar vl_prof=2.0*vlbar*(r3*r3-r2*r2)/Dhap*(r_l*r_l-r2*r2+2.0*r3*r3*log(r2/r_l));   // Eq. 3.17

  FFVar sigT=gsig+tauT*gsig_dT*(Tlo-gTlref);                    // Eq 3.40 feedback: surface tension drops with solved dT
  FFVar delta_w=(2.0*sigT*gabscos)/dP_z, uu=delta_w/gdbar;
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
  // ---- lumped energy balance (Eq 3.40), realized as a DISTRIBUTED-CONSTANT liquid temperature Tlo(z) ----
  //   Tlo lives on {z} (like Cg/Vg) so it has collocation nodes; a true 0-dim scalar has no node to attach
  //   to and makes the principal symbol rectangular.  Flat transport pins it constant in z; the single
  //   z=0 (LB) boundary row carries the energy balance, where Vg,Cg are already the gas-OUTLET values:
  //      total CO2 removed [mol/s] = (N*pi*r1^2*Cref)*vgIn*(1 - Vg(0)Cg(0)).
  FFVar TL_FLAT = OpP(Tlo,z);                                              // dTlo/dz = 0  -> Tlo constant along the fiber
  FFVar co2rem_eb = (nf*M_PI*r1*r1*Cref)*vgIn*(1.0-Vg*Cg);                 // at z=0 (LB): total CO2 removal [mol/s]
  FFVar TL_EB = Tlo - V.Tl - co2rem_eb*gdHr/(grhol*V.Flin*gCp);            // Eq 3.40: Tlo(0) = Tl_in + Phi*Am*|dHr|/(rho_l Fl Cp)

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
  oc.add_state(Tlo,{z});                                        // distributed-constant liquid temperature (lumped Eq 3.40)
  oc.add_input(lam,{}); oc.add_input(mu,{}); oc.add_input(kap,{}); oc.add_input(tauT,{});
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
  oc.update_ref(Cw,[&](OCFESLV::t_Coord const&){return 0.63;});    // nominal partH (Ch-3 N2O analogy @ Tl=308)
  oc.update_ref(Mw,[&](OCFESLV::t_Coord const&){return 1.0;});
  oc.update_ref(Pw,[&](OCFESLV::t_Coord const&){return 1.0;});
  oc.update_ref(Cl,[&](OCFESLV::t_Coord const&){return 0.0;});
  oc.update_ref(Ml,[&](OCFESLV::t_Coord const&){return 1.0;});
  oc.update_ref(Pl,[&](OCFESLV::t_Coord const&){return 1.0;});
  oc.update_ref(omg,[&](OCFESLV::t_Coord const&){return gomega_fix;});
  oc.update_ref(Tlo,[&](OCFESLV::t_Coord const&){return 308.0;});   // start at inlet Tl (dT=0)

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
  oc.add_equation(TL_FLAT,{z},{Z_NO_LB},blk);                    // dTlo/dz=0 on all z except LB (mirrors GAS_V/GAS_V_IN split)
  oc.add_equation(TL_EB,{z},{FFDom::LB},bnd);                    // energy balance as the single z=0 boundary condition
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

// =============================================================================================
//  CONTINUATION -- now DECLARED, not hand-rolled.
//
//  The parameter ramps this driver used to implement by hand are registered with OCFESLV via
//  add_homotopy() and executed inside the ordinary solve(), exactly as SOLVE_MARCHING is.  The
//  STAGE index is the whole schedule grammar:
//      all stage 0        -> SIMULTANEOUS (one s: 0->1 drives every parameter)
//      stages 0,1,2,3     -> STAIRCASE    (each parameter ramped in turn)
//      stages 0,0,1,1     -> PARTIAL grouping
//  Per-stage cost/backtrack/predictor detail comes from oc.continuation_report().
//
//  *** WHY THE CONTINUATION CANNOT BE REPLACED (measured 2026-07-23) ***
//  The collocated MBC system has MULTIPLE roots.  Jumping a parameter straight to 1 and ramping it
//  both converge to |r|~1e-11 but land on DIFFERENT solutions (||dx||_inf = 6.4).  The homotopy
//  path is a ROOT-SELECTION mechanism, not merely a convergence aid:
//    * a single full-physics solve diverges (more iterations make it worse),
//    * one-jump-per-parameter staging lands on a spurious root and is ~3x SLOWER,
//    * a reference calibrated to the converged profiles gives a 293x better initial residual
//      and NO speedup -- proximity does not select a root.
//  Consequence: "converged" is NOT evidence of correctness here, and HOM_STEP_CAP is a CORRECTNESS
//  knob, not just a speed knob.  Validate any schedule/cap change with mode 'simulcmp'.
// =============================================================================================

//! MEASURED (2026-07-23): MBC5 tracks the simultaneous path safely at 0.1 on BOTH vert and horiz
//! -- ||x_stair - x_simul||_inf ~ 5e-12, i.e. the same root -- for a 2.6x speedup (40 -> 10 solves).
static double gHomCap = 0.1;

//! Schedule: stage index per parameter, in the order (lam, mu, kap, tauT).
//! Default = all stage 0 = simultaneous, i.e. a SINGLE ramping parameter.
static int gStage[4] = { 0, 0, 0, 0 };

//! Register the schedule on an OCFESLV already through setup() (build_mbc calls it).
static void set_schedule( OCFESLV& oc, Vars const& V, int const stg[4], double cap )
{
  oc.options.HOMOTOPY.STEP_CAP = cap;
  oc.add_homotopy( V.lam , 0.0, 1.0, stg[0] );
  oc.add_homotopy( V.mu  , 0.0, 1.0, stg[1] );
  oc.add_homotopy( V.kap , 0.0, 1.0, stg[2] );
  oc.add_homotopy( V.tauT , 0.0, 1.0, stg[3] );
}

//! Print the per-stage breakdown of the continuation just performed.
static void report_continuation( OCFESLV const& oc, std::string const& tag )
{
  OCFESLV::ContinuationReport const& r = oc.continuation_report();
  for( auto const& st : r.stage )
    std::cout<<"    [stage "<<st.stage<<"] "<<st.nparam<<" param(s)  "
             <<std::setw(3)<<st.solves<<" solves  "<<std::setw(3)<<st.accepts<<" acc  "
             <<std::setw(2)<<st.backtracks<<" bt  "<<std::setw(3)<<st.pred_used<<" pred+ ("
             <<st.pred_ord1<<" lin/"<<st.pred_ord2<<" quad)  "
             <<std::setw(3)<<st.pred_rejected<<" pred-  "<<std::setw(4)<<st.iterations<<" its  "
             <<std::fixed<<std::setprecision(2)<<st.seconds<<" s\n";
  std::cout<<"    ["<<tag<<" TOTAL] "<<r.solves<<" solves, "<<r.iterations<<" its, "
           <<std::fixed<<std::setprecision(2)<<r.seconds<<" s\n";
}

struct KPI { bool conv=false; double resid=1e30,flux=0,eta=0,load=0,Tout=0,hc=0,evap=0,mdea=0,pz=0; };

static void write_profiles( OCFESLV& oc, Vars const& V, double const* xv, double const* inp, std::string const& tag )
{
  std::ofstream fC("fig10_"+tag+"_CO2.dat"),fM("fig10_"+tag+"_MDEA.dat"),fP("fig10_"+tag+"_PZ.dat");
  fC<<std::scientific<<std::setprecision(6); fM<<std::scientific<<std::setprecision(6); fP<<std::scientific<<std::setprecision(6);
  auto ev=[&](FFVar const& X,OCFESLV::t_Coord const& pt){ return oc.eval_colloc<double>(X,pt,xv,inp,nullptr); };
  int const NR=100,NZ=100;
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

static KPI run(double hydro,std::string const& tag, bool predict=true)
{
  KPI K;
  FFGraph DAG; OCFESLV oc(&DAG); Vars V;
  build_mbc(DAG,oc,V,hydro);
  std::vector<double> xv,inp; if(!oc.init(xv,inp,nullptr)){ std::cerr<<"  ["<<tag<<"] init failed\n"; return K; }
  oc.set_input_values(V.lam,{0.0},inp.data());
  oc.set_input_values(V.mu,{0.0},inp.data());
  oc.set_input_values(V.kap,{0.0},inp.data());
  oc.set_input_values(V.tauT,{0.0},inp.data());                        // sigma(T) feedback OFF: isothermal baseline
  set_schedule( oc, V, gStage, gHomCap );      // schedule declared; solve() runs it
  std::cout<<"  ["<<tag<<"] continuation (cap "<<gHomCap<<"):\n";
  OCFESLV::SolveReport sr = oc.solve( xv.data(), inp.data(), nullptr );
  K.resid = sr.final_residual;
  report_continuation( oc, tag );
  if( !sr.converged ){ std::cerr<<"  ["<<tag<<"] continuation failed\n"; return K; }
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

// ---- Schedule comparison: do two schedules reach the SAME root, and at what cost? -------------
// Because the system is multi-rooted, BOTH runs report converged regardless -- only the state
// difference reveals whether a schedule skipped a branch.
static bool run_schedule( double hydro, int const stg[4], double cap,
                          std::vector<double>& xout, int& nsolve, double& secs )
{
  FFGraph DAG; OCFESLV oc(&DAG); Vars V; build_mbc(DAG,oc,V,hydro);
  std::vector<double> xv,inp;
  if(!oc.init(xv,inp,nullptr)) return false;
  set_schedule( oc, V, stg, cap );
  OCFESLV::SolveReport sr = oc.solve( xv.data(), inp.data(), nullptr );
  OCFESLV::ContinuationReport const& r = oc.continuation_report();
  nsolve = r.solves; secs = r.seconds; xout = xv;
  return sr.converged;
}

static double dinf( std::vector<double> const& a, std::vector<double> const& b )
{
  if( a.size()!=b.size() || a.empty() ) return 1e30;
  double m=0.; for(size_t i=0;i<a.size();++i) m=std::max(m,std::fabs(a[i]-b[i]));
  return m;
}

//! STAIRCASE vs SIMULTANEOUS at the adopted cap -- the standing root-identity guard.
static void simul_compare( double hydro, std::string const& tag )
{
  std::cout<<"\n================ staircase vs simultaneous ["<<tag<<"]  cap="<<gHomCap<<" ================\n";
  int const stair[4] = { 0,1,2,3 }, simul[4] = { 0,0,0,0 };
  std::vector<double> xs,xm; int ns=0,nm=0; double ss=0.,sm=0.;
  if( !run_schedule(hydro,stair,gHomCap,xs,ns,ss) ){ std::cerr<<"  staircase FAILED\n"; return; }
  std::cout<<"  staircase   : "<<ns<<" solves, "<<std::fixed<<std::setprecision(2)<<ss<<" s\n";
  if( !run_schedule(hydro,simul,gHomCap,xm,nm,sm) ){ std::cerr<<"  simultaneous FAILED\n"; return; }
  std::cout<<"  simultaneous: "<<nm<<" solves, "<<sm<<" s   speedup x"<<(sm>0?ss/sm:0.)<<"\n";
  double const d = dinf(xs,xm);
  std::cout<<"  ||x_stair - x_simul||_inf = "<<std::scientific<<std::setprecision(3)<<d<<"\n";
  std::cout<<"  -> "<<( d<1e-6 ? "SAME root: the simultaneous schedule is SAFE at this cap"
                               : "DIFFERENT roots: this cap skips a branch -- REDUCE gHomCap" )<<"\n";
}

//! Cap sweep: smallest cap at which the simultaneous schedule still reproduces the staircase.
static void simul_sweep( double hydro, std::string const& tag )
{
  std::cout<<"\n================ cap sweep ["<<tag<<"] ================\n";
  int const stair[4] = { 0,1,2,3 }, simul[4] = { 0,0,0,0 };
  std::vector<double> xref; int nref=0; double sref=0.;
  if( !run_schedule(hydro,stair,gHomCap,xref,nref,sref) ){ std::cerr<<"  staircase FAILED\n"; return; }
  std::cout<<"  staircase reference: "<<nref<<" solves, "<<std::fixed<<std::setprecision(2)<<sref<<" s\n\n";
  std::cout<<"    "<<std::left<<std::setw(8)<<"cap"<<std::setw(9)<<"solves"<<std::setw(10)<<"time[s]"
           <<std::setw(11)<<"speedup"<<std::setw(13)<<"||dx||_inf"<<"verdict\n";
  for( double cap : {0.1,0.05,0.025,0.0125} ){
    std::vector<double> xs; int ns=0; double ss=0.;
    if( !run_schedule(hydro,simul,cap,xs,ns,ss) ){
      std::cout<<"    "<<std::left<<std::setw(8)<<std::fixed<<std::setprecision(4)<<cap<<"FAILED\n"; continue; }
    double const d = dinf(xs,xref);
    bool const same = (d<1e-6);
    std::cout<<"    "<<std::left<<std::setw(8)<<std::fixed<<std::setprecision(4)<<cap
             <<std::setw(9)<<ns<<std::setw(10)<<std::setprecision(2)<<ss
             <<std::setw(11)<<(ss>0?sref/ss:0.)
             <<std::setw(13)<<std::scientific<<std::setprecision(2)<<d
             <<( same ? "SAME root" : "DIFFERENT root" )<<"\n";
    if( same ){ std::cout<<"    -> smallest validated cap: "<<std::fixed<<std::setprecision(4)<<cap
                         <<"  (set gHomCap to this)\n"; break; }
  }
}

int main(int argc,char** argv)
{
  std::string const mode = (argc>1) ? argv[1] : "fig";
  std::cout<<"================================================================\n"
           <<"  MBC5 test : Chapter-3 module (phi=0.55, no insert), latest correlations\n"
           <<"================================================================\n";

  if( mode=="simulsweep" ){             // how fine must the diagonal ramp be to stay on-branch?
    simul_sweep( grhol*gGrav, "vert" );
    simul_sweep( 0.0,         "horiz" );
    return 0;
  }
  if( mode=="simulcmp" ){               // staircase vs simultaneous ramp: same root? how much faster?
    simul_compare( grhol*gGrav, "vert" );
    simul_compare( 0.0,         "horiz" );
    return 0;
  }

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
