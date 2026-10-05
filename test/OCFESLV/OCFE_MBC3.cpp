// OCFE_MBC3.cpp  ---  Stage 3: PZ Damköhler robustness sweep for the three-species MBC
// ============================================================================================
// Drives kP (PZ rate) alone from the Stage-2b plateau-check value up to the real Arrhenius
// regime (Da_P ~ 3.6e6), on two meshes, to separate PHYSICS from RESOLUTION:
//
//   (A) FIXED mesh  = the Stage-2b graded(6/7, 3.3) wet/liquid films  -> maps where it breaks.
//   (B) ADAPTIVE mesh = wet/liquid element counts sized so the finest element tracks the PZ
//       reaction-layer thickness  delta_P ~ 1/sqrt(aPC*Pw)  -> holds resolution as Da grows.
//
// The PZ film balance is  d2Cw/drwc^2 = aPC*Pw*Cw + aMC*Mw*Cw,  aPC = (Drw^2/DwC)*kP*Pref,
// so the reaction layer at the dry-wet interface (rwc=0) has scaled width delta_P = 1/sqrt(aPC).
// Because a geometric graded mesh needs finest element w0 = (r-1)/(r^nel-1) ~ delta_P, nel grows
// only ~0.5*log(Da)/log(r): Da_P 4e4 -> nel 6, 3.6e6 -> nel 8.  No ill-conditioning, which is why
// a coordinate-stretched film subdomain is not needed at this Da.
//
// Diagnostics per run:
//   * continuation convergence (does the lam homotopy reach 1) and final |r|;
//   * CO2 absorption flux Phi (Eq 31) -- DIFFUSION-LIMITED, so it must PLATEAU across the sweep;
//   * regime probe: gas-side amine Pw(0)/Mw(0) vs liquid-side Pw(1)/Mw(1).  The reaction is NOT a
//     gas-side pseudo-first-order film -- PZ/MDEA enter from the liquid (rwc=1) and are consumed
//     before reaching rwc=0, so it sits at an interior front; we locate that front and the CO2
//     1%-penetration depth, and count the mesh nodes bracketing the front.
//
// PASS: the ADAPTIVE mesh reaches Da_P=3.6e6 converged, with Phi on the low-Da plateau, and the
// amine-limited interior-front regime confirmed (Pw(0)<<Pw(1)).  The mesh-vs-front node count is
// reported (and flagged if the rwc=0 clustering under-resolves the interior front -- the metric that
// matters for the amine-consumption optimization objective, distinct from the transport-limited Phi).
// ============================================================================================
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
  std::cout<<"  "<<std::left<<std::setw(52)<<n<<"  "<<(ok?"PASS":"FAIL")<<"\n"; }

// geometric grading: element widths (sum=1), smallest first -> clustered at the low end (interface).
static std::vector<double> graded(int nel,double ratio){
  double w0=(ratio-1.0)/(std::pow(ratio,nel)-1.0);
  std::vector<double> w(nel); for(int i=0;i<nel;++i) w[i]=w0*std::pow(ratio,i); return w; }

struct Res {
  bool   conv=false;
  double lam=0.0, resid=1e30;
  double DaP=0.0, deltaP=0.0, Phi=0.0;
  int    nwet=0, nliq=0, nodesInLayer=0;
  double Cw_i=0.0, gradAct=0.0, gradFilm=0.0;
  double Pw0=0.0, Mw0=0.0, Pw1=0.0, Mw1=0.0;   // gas-side (rwc=0) vs liquid-side (rwc=1) amine
  double rwc_pen=1.0, rwc_front=0.0;           // CO2 1%-penetration depth; reaction-rate-max location
  int    nodesNearFront=0;                     // mesh nodes bracketing the interior front
};

static Res run(double kP_scale, bool adaptive, double ratio=3.3,
               int imp_sel = -1,                                   // -1: keep the MBC3_WEAK selector
               OCFESLV::Options::ReductionType red = OCFESLV::Options::RED_FULL)
{
  Res R;
  // ---- physical constants (identical to Stage-2b) ------------------------------------------
  double const L=2.3,r1=225e-6,r2=550e-6,eta=0.40,epsm=0.41,tau=6.1; int const nf=209;
  double const omega=0.20,rw=r2-(r2-r1)*omega,r3=r2/std::sqrt(eta);
  double const Drd=rw-r1,Drw=r2-rw,Drl=r3-r2,Rref=r2;
  double const r1s=r1/Rref,r2s=r2/Rref,rws=rw/Rref,r3s=r3/Rref; (void)r3s;
  double const Drds=Drd/Rref,Drws=Drw/Rref,Drls=Drl/Rref;
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
  // PZ rate swept alone; MDEA rate fixed.
  double const kM=9.853e-2, kP=8.784e-1*kP_scale;
  double const cgf=(2.0*Dg/r1)/Drd, cgfV=cgf*yin/vgIn, rDWf=DwC*Drd/(Dg*Drw), rWLf=(epsm/tau)*(Drl/Drw);
  double const aMC=(Drw*Drw/DwC)*kM*Mref, aPC=(Drw*Drw/DwC)*kP*Pref;   // wet CO2 sinks
  double const aMM=(Drw*Drw/DwM)*kM*Cref, aPP=(Drw*Drw/DwP)*kP*Cref;   // wet amine sinks
  double const bMC=kM*Mref, bPC=kP*Pref, bMM=kM*Cref, bPP=kP*Cref;     // liquid sinks

  R.DaP    = aPC;
  R.deltaP = 1.0/std::sqrt(aPC);   // scaled reaction-layer thickness at the interface (Pw~1)

  // ---- mesh: fixed (Stage-2b) or Da-adaptive (finest element ~ delta_P) --------------------
  int nwet, nliq;
  if(adaptive){
    // choose nel so the finest element w0=(r-1)/(r^nel-1) <= delta_P (~1 element per layer, 8 nodes)
    nwet = (int)std::ceil( std::log( (ratio-1.0)/R.deltaP + 1.0 )/std::log(ratio) );
    nwet = std::max(6, nwet);
    nliq = nwet + 1;
  } else { nwet = 6; nliq = 7; }
  if( char const* e = std::getenv("MBC3_NWET") ){ nwet = std::atoi(e); nliq = nwet + 1; }
  R.nwet=nwet; R.nliq=nliq;
  std::vector<double> wwet=graded(nwet,ratio), wliq=graded(nliq,ratio);

  // ---- DAG / model (identical structure to Stage-2b) ---------------------------------------
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
  FFVar DW_MFLUX=OpP(Mw,rwc), DW_PFLUX=OpP(Pw,rwc);
  FFVar WL_CVAL=Cw-Cl, WL_CFLUX=rWLf*OpP(Cw,rwc)-OpP(Cl,rl);
  FFVar WL_MVAL=Mw-Ml, WL_MFLUX=rWLf*OpP(Mw,rwc)-OpP(Ml,rl);
  FFVar WL_PVAL=Pw-Pl, WL_PFLUX=rWLf*OpP(Pw,rwc)-OpP(Pl,rl);
  FFVar L_C_WALL=OpP(Cl,rl), L_M_WALL=OpP(Ml,rl), L_P_WALL=OpP(Pl,rl);
  FFVar L_C_IN=Cl-ClInS, L_M_IN=Ml-1.0, L_P_IN=Pl-1.0;

  OCFESLV oc(&DAG);
  oc.add_domain(z,  FFDom(0.,L,  3,FFDom::CGL,6));
  oc.add_domain(rd, FFDom(0.,1.0,3,FFDom::CGL,7));
  oc.add_domain(rwc,FFDom(0.0, wwet, FFDom::CGL,8));
  oc.add_domain(rl, FFDom(0.0, wliq, FFDom::CGL,8));
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

  oc.options.REDUCE.ORDER=red;
  oc.options.CLASSIFY.MODE=OCFESLV::Options::CLASS_AUTO;
  oc.options.INTERFACE.TYPE=OCFESLV::Options::IC_AUTO;
  // 2026-09-19: the matrix passes the imposition explicitly (imp_sel 0/1/2); MBC3_WEAK stays the default.
  oc.options.INTERFACE.IMPOSITION = ( imp_sel == 0 ) ? OCFESLV::Options::IC_WEAK
                             : ( imp_sel == 1 ) ? OCFESLV::Options::IC_TRACE
                             : ( imp_sel == 2 ) ? OCFESLV::Options::IC_STRONG
                             : ( std::getenv("MBC3_WEAK") ? OCFESLV::Options::IC_WEAK : OCFESLV::Options::IC_STRONG );
  oc.options.INTERFACE.SAT_SIGMA0=10.0;
  if(char const* e=std::getenv("MBC3_SIG0")) oc.options.INTERFACE.SAT_SIGMA0=std::atof(e);
  oc.options.DISPLAY_LEVEL=0; oc.options.SOLVE.VERBOSE=false;
#if defined(CRONOS__WITH_SPQR)
  oc.options.SOLVE.FACTORIZATION=OCFESLV::Options::SOLVE_SPQR;
#else
  oc.options.SOLVE.FACTORIZATION=OCFESLV::Options::SOLVE_SUPERLU;
#endif

  if(!oc.setup()){ std::cerr<<"  setup() FAILED (DaP="<<R.DaP<<")\n"; return R; }
  std::vector<double> xv,inp;
  if(!oc.init(xv,inp,nullptr)){ std::cerr<<"  init() FAILED\n"; return R; }
  double const* ip=inp.empty()?nullptr:inp.data();

  // lam continuation (denser near 1 for the stiff high-Da film).
  OCFESLV::SolveReport rep;
  double const lam_seq[]={0.0,0.01,0.03,0.06,0.1,0.2,0.35,0.5,0.65,0.8,0.9,0.95,1.0};
  for(double lv:lam_seq){
    oc.set_input_values(lam,{lv},inp.data());
    rep=oc.solve(xv.data(),inp.data(),nullptr);
    R.lam=lv; R.resid=rep.final_residual;
    if(!rep.converged) break;
  }
  R.conv=rep.converged;
  if(!R.conv) return R;

  // ---- diagnostics ----
  auto ev=[&](FFVar const& V,OCFESLV::t_Coord const& pt){ return oc.eval_colloc<double>(V,pt,xv.data(),ip,nullptr); };
  OCFESLV::t_Coord gi;gi[z]=L;   double Cgi=ev(Cg,gi);
  OCFESLV::t_Coord go;go[z]=0.0; double Cgo=ev(Cg,go);
  OCFESLV::t_Coord vL;vL[z]=L;   double VgL=ev(Vg,vL);
  OCFESLV::t_Coord v0;v0[z]=0.0; double Vg0=ev(Vg,v0);
  R.Phi=(r1/(2.0*L))*vgIn*Cref*( VgL*Cgi - Vg0*Cgo );

  // reaction-layer resolution: count rwc collocation nodes within delta_P of the interface (rwc=0).
  FFDom drwc(0.0, wwet, FFDom::CGL, 8); drwc.set_nodes();
  int nin=0;
  for(size_t ie=0;ie<drwc.n_elem;++ie){ auto e=drwc.lgnodes(drwc.elem_lo(ie),drwc.elem_up(ie));
    for(double rn:e){ if(rn<=R.deltaP) ++nin; } }
  R.nodesInLayer=nin;

  // ---- reaction-front / regime characterization at a mid-column --------------------------
  // The film check assuming a gas-side pseudo-first-order layer does NOT apply here: PZ/MDEA
  // enter from the LIQUID side (rwc=1) and are consumed before reaching rwc=0, so the reaction
  // sits at an interior front, not the interface.  Probe the depletion and locate the front.
  double const zc=0.5*L, dd=std::max(1e-9, 0.1*R.deltaP);
  auto atrwc=[&](FFVar const& V,double rr){ OCFESLV::t_Coord pp; pp[rwc]=rr; pp[z]=zc; return ev(V,pp); };
  double const Cw0=atrwc(Cw,0.0), Cw1=atrwc(Cw,dd);
  R.Cw_i    = Cw0;
  R.Pw0 = atrwc(Pw,0.0);  R.Mw0 = atrwc(Mw,0.0);    // gas-side interface amine (depletion probe)
  R.Pw1 = atrwc(Pw,1.0);  R.Mw1 = atrwc(Mw,1.0);    // liquid-side amine (the source)
  R.gradAct = std::abs( (Cw1-Cw0)/dd );             // CO2 interface gradient (transport-limited)
  R.gradFilm= std::abs( Cw0 )*std::sqrt( std::max(0.0, aPC*R.Pw0 + aMC*R.Mw0) );  // guarded film ref

  // scan the wet membrane: CO2 1%-penetration depth and the reaction-rate-max (front) location.
  int const NS=400;
  double const Cwref=std::max(1e-30,std::abs(Cw0));
  double rate_max=0.0; bool penHit=false;
  for(int i=0;i<=NS;++i){
    double rr=(double)i/NS;
    double cw=atrwc(Cw,rr), pw=atrwc(Pw,rr), mw=atrwc(Mw,rr);
    double rate=std::abs(cw)*( aPC*std::max(0.0,pw) + aMC*std::max(0.0,mw) );
    if(rate>rate_max){ rate_max=rate; R.rwc_front=rr; }
    if(!penHit && std::abs(cw)<=0.01*Cwref){ R.rwc_pen=rr; penHit=true; }
  }
  if(!penHit) R.rwc_pen=1.0;

  // mesh nodes bracketing the front (front resolution -- what the amine-consumption objective needs).
  { FFDom d(0.0, wwet, FFDom::CGL, 8); d.set_nodes(); int cnt=0;
    for(size_t ie=0;ie<d.n_elem;++ie){ auto e=d.lgnodes(d.elem_lo(ie),d.elem_up(ie));
      for(double rn:e) if(std::abs(rn-R.rwc_front)<=0.05) ++cnt; }
    R.nodesNearFront=cnt; }
  return R;
}

int main()
{
  std::cout<<"================================================================\n"
           <<"  Stage 3: PZ Damkohler robustness sweep (kP alone -> Da_P=3.6e6)\n"
           <<"================================================================\n";

  double const scales[]={1.0,10.0,30.0,90.0};   // Da_P ~ 4e4, 4e5, 1.2e6, 3.6e6
  double const targetDa=3.6e6;

  auto sweep=[&](bool adaptive){
    std::cout<<"\n  ---- "<<(adaptive?"ADAPTIVE mesh (finest ~ delta_P)":"FIXED mesh (Stage-2b 6/7)")<<" ----\n";
    std::cout<<"  "<<std::left
             <<std::setw(11)<<"Da_P"<<std::setw(6)<<"nwet"<<std::setw(7)<<"conv"
             <<std::setw(11)<<"|r|"<<std::setw(11)<<"Phi"<<std::setw(11)<<"Pw(0)"
             <<std::setw(11)<<"Pw(1)"<<std::setw(10)<<"r_front"<<std::setw(10)<<"r_pen"<<std::setw(7)<<"nfront"<<"\n";
    std::vector<Res> out;
    for(double s:scales){
      Res R=run(s,adaptive);
      out.push_back(R);
      std::string convs = R.conv ? std::string("y") : ("n@"+std::to_string(R.lam));
      std::cout<<"  "<<std::left<<std::scientific<<std::setprecision(2)
               <<std::setw(11)<<R.DaP<<std::setw(6)<<R.nwet<<std::setw(7)<<convs
               <<std::setw(11)<<R.resid<<std::setw(11)<<R.Phi<<std::setw(11)<<R.Pw0
               <<std::setw(11)<<R.Pw1<<std::setw(10)<<R.rwc_front<<std::setw(10)<<R.rwc_pen
               <<std::setw(7)<<R.nodesNearFront<<"\n";
    }
    return out;
  };

  std::vector<Res> fixedR = sweep(false);
  std::vector<Res> adaptR = sweep(true);

  // ---- verdicts ----
  std::cout<<"\n  ---- verdict ----\n";
  Res const& lo = adaptR.front();              // reference plateau value (lowest Da, resolved)
  Res const& hi = adaptR.back();               // target Da_P ~ 3.6e6
  check_true("adaptive mesh converges at Da_P=3.6e6", hi.conv && hi.DaP>=0.9*targetDa);
  if(lo.conv && hi.conv){
    double const plateauErr = std::abs(hi.Phi-lo.Phi)/std::max(1e-30,std::abs(lo.Phi));
    std::cout<<"    Phi plateau: lowDa="<<std::scientific<<std::setprecision(3)<<lo.Phi
             <<"  highDa="<<hi.Phi<<"  rel.dev="<<plateauErr<<"\n";
    check_true("Phi stays on the diffusion-limited plateau (<5%)", plateauErr<5e-2);

    // regime characterization at Da_P=3.6e6: amine-limited INTERIOR front, not a gas-side film.
    std::cout<<std::scientific<<std::setprecision(2)
             <<"    regime @Da_P=3.6e6: gas-side amine Pw(0)="<<hi.Pw0<<" Mw(0)="<<hi.Mw0
             <<" ; liquid-side Pw(1)="<<hi.Pw1<<" Mw(1)="<<hi.Mw1<<"\n"
             <<std::fixed<<std::setprecision(3)
             <<"                       reaction front at rwc="<<hi.rwc_front
             <<" ; CO2 1%-penetration rwc="<<hi.rwc_pen
             <<" ; mesh nodes near front="<<hi.nodesNearFront<<"\n";
    bool const depleted = hi.Pw0 < 0.1*std::max(1e-30,hi.Pw1);
    check_true("amine-limited interior-front regime confirmed (Pw(0)<<Pw(1))", depleted);
    if(hi.nodesNearFront<5)
      std::cout<<"    note: front is under-resolved by the rwc=0 clustering (nfront="<<hi.nodesNearFront
               <<"); to resolve the amine-consumption objective the grading should target rwc~"
               <<std::setprecision(2)<<hi.rwc_front<<", not the gas-side film.\n";
    else
      check_true("reaction front resolved by mesh (>=5 nodes near front)", true);
  }
  std::cout<<std::scientific<<std::setprecision(3);
  std::cout<<"    (fixed mesh at Da_P=3.6e6: conv="<<(fixedR.back().conv?"y":"n")
           <<" nfront="<<fixedR.back().nodesNearFront
           <<" Phi="<<fixedR.back().Phi<<")\n";

  // -------------------------------------------------------------------------------------------------------
  // 2026-09-19 -- THE SYSTEMATIC MATRIX: 3 impositions x 2 reductions at the mildest Damkohler (scale 1, fixed
  // mesh).  MBC3 is the only MBC driver that already ran a second imposition, and the MBC family carries the
  // corpus's entire `non-constant` refusal population (47,586 edges), which M7b measured INERT to the rescue
  // fallback under STRONG.  This matrix says whether that holds under TRACE and under RED_MAIN as well.
  // A cell expected to differ gets a documented XFAIL naming the mechanism, never a loosened bar.
  // -------------------------------------------------------------------------------------------------------
  std::cout<<"\n---- SYSTEMATIC MATRIX (Da_P scale 1, fixed mesh): imposition x reduction ----\n";
  {
    struct MCell { char const* imp; int sel; char const* red; OCFESLV::Options::ReductionType rt; };
    MCell const cells[6] = {
      { "WEAK",   0, "RED_FULL", OCFESLV::Options::RED_FULL },
      { "TRACE",  1, "RED_FULL", OCFESLV::Options::RED_FULL },
      { "STRONG", 2, "RED_FULL", OCFESLV::Options::RED_FULL },
      { "WEAK",   0, "RED_MAIN", OCFESLV::Options::RED_MAIN },
      { "TRACE",  1, "RED_MAIN", OCFESLV::Options::RED_MAIN },
      { "STRONG", 2, "RED_MAIN", OCFESLV::Options::RED_MAIN } };
    bool matrix_ok = true;
    std::cout<<"  "<<std::left<<std::setw(10)<<"IMPOSITION"<<std::setw(11)<<"REDUCTION"
             <<std::setw(7)<<"conv"<<std::setw(12)<<"|r|"<<std::setw(12)<<"Phi"<<"VERDICT\n";
    for( auto const& mc : cells ){
      Res R = run( 1.0, false, 3.3, mc.sel, mc.rt );
      bool const xfail = false;            // none expected at the mildest Damkohler
      if( !R.conv && !xfail ) matrix_ok = false;
      std::cout<<"  "<<std::left<<std::setw(10)<<mc.imp<<std::setw(11)<<mc.red
               <<std::setw(7)<<( R.conv ? "y" : "n" )<<std::scientific<<std::setprecision(2)
               <<std::setw(12)<<R.resid<<std::setw(12)<<R.Phi
               <<( R.conv ? ( xfail ? "PASS (unexpected: the documented defect is gone?)" : "PASS" )
                          : ( xfail ? "XFAIL (documented)" : "FAIL" ) )<<"\n";
    }
    check_true( "matrix: all 6 imposition x reduction cells converge", matrix_ok );
    std::cout<<"  MATRIX: "<<( matrix_ok ? "PASS" : "FAIL" )<<"\n";
  }

  std::cout<<"\n  RESULT: "<<g_pass<<" passed, "<<g_fail<<" failed\n";
  return g_fail?1:0;
}
