// OCFE_MBC2a.cpp  ---  Stage 2a: DIMENSIONLESS multidomain MBC (simple reaction)
// ===========================================================================
// Dimensional geometry/transport from the paper (Table 1/2), but the model is
// solved in NON-DIMENSIONAL form: CO2 concentrations scaled by Cref=CgIn, amine by
// Sref=CsolIn, gas velocity by vgIn.  This is the PRINCIPLED conditioning fix:
// Stage 1 (dimensionless) solved this exact multidomain structure + reaction +
// interface corner-contention to machine precision; only the O(600) dimensional
// scale ill-conditioned the IC_STRONG geometric-C0 Schur-pivot fallback (evidenced
// by the mesh-dependent floor: nelr 3->4 worsened it 1e-12 -> 5e-8 relative).  With
// O(1) variables and O(1) coefficients the well-conditioned regime is restored and
// SOLVE_RES_TOL=1e-9 is again correct.  Radial subdomains rescaled to [0,1] with
// physical extents in coefficients; cylindrical Laplacian; real-gas Z; interface
// fluxes row-scaled.  Reaction fraction lambda in [0,1] is the continuation input.
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

static int g_pass = 0, g_fail = 0;
static void check_true( char const* n, bool ok ){ (ok?g_pass:g_fail)++;
  std::cout << "  " << std::left << std::setw(48) << n << "  " << (ok?"PASS":"FAIL") << "\n"; }

int main()
{
  // ---- geometry (Table 2) ----
  double const L=2.3, r1=225e-6, r2=550e-6, eta=0.40, epsm=0.41, tau=6.1;
  int    const nf=209;
  double const omega=0.20, rw=r2-(r2-r1)*omega, r3=r2/std::sqrt(eta);
  double const Drd=rw-r1, Drw=r2-rw, Drl=r3-r2;
  // radial length-scale rescaling -> O(1) dimensions in ALL FFVar intermediates
  double const Rref=r2;
  double const r1s=r1/Rref, r2s=r2/Rref, rws=rw/Rref;
  double const Drds=Drd/Rref, Drws=Drw/Rref, Drls=Drl/Rref;
  // ---- transport (representative SI) ----
  double const Dg=8.8e-8, DwC=(epsm/tau)*1.0e-9, DwS=(epsm/tau)*5.0e-10, DlC=1.0e-9, DlS=5.0e-10, H=3000.0;
  // ---- gas/liquid state (Table 1 lab point) ----
  double const Z=0.9, Rg=8.314, Tg=293.0, P=54e5, yin=0.24;
  double const Cref=yin*P/(Z*Rg*Tg), Sref=3000.0;                 // scaling references
  double const MWg=yin*0.044+(1-yin)*0.028, rhog=P*MWg/(Z*Rg*Tg);
  double const vgIn=(3.0/3600.0)/(nf*M_PI*r1*r1*rhog);
  double const vlbar=(10e-3/3600.0)/(nf*M_PI*r3*r3*(1-eta));
  double const ClInS=(0.22*Sref)/Cref;                            // scaled liquid CO2 feed
  double const partH=Z*Rg*Tg/H;                                   // Henry partition ~0.73
  double const krxn=2.0e-6, nuS=1.0;                              // Da~0.4 target (real k in 2b)
  // ---- dimensionless coefficients (all O(1) or smaller) ----
  double const DlCs=DlC/(Rref*Rref), DlSs=DlS/(Rref*Rref);        // effective radial diffusivities
  double const cgf   = (2.0*Dg/r1)/Drd;                           // ~3.0
  double const cgfV  = cgf*yin/vgIn;                              // ~2.26
  double const rDWf  = DwC*Drd/(Dg*Drw);                          // ~3.1e-3 (row-scaled flux ratio)
  double const rWLf  = (epsm/tau)*(Drl/Drw);                      // ~0.33
  double const kWC   = (Drw*Drw/DwC)*Sref*krxn;                   // WET_C reaction coeff ~0.38
  double const kWS   = (Drw*Drw/DwS)*Cref*krxn*nuS;               // WET_S reaction coeff ~0.15
  double const kLC   = Sref*krxn;                                 // LIQ_C reaction coeff ~6e-3
  double const kLS   = Cref*krxn*nuS;                             // LIQ_S reaction coeff ~1.2e-3

  std::cout << "================================================================\n"
            << "  Stage 2a: DIMENSIONLESS multidomain MBC (simple reaction)\n"
            << "  SI geometry; concentrations/velocity scaled to O(1)\n"
            << "  Cref=" << std::fixed << std::setprecision(1) << Cref
            << " mol/m3  Sref=" << std::setprecision(0) << Sref << " mol/m3  vgIn="
            << std::setprecision(3) << vgIn << " m/s\n"
            << "================================================================\n";

  FFGraph DAG;
  FFVar z=DAG.add_var("z"), rd=DAG.add_var("rd"), rwc=DAG.add_var("rwc"), rl=DAG.add_var("rl");
  FFVar Cg=DAG.add_var("Cg(z)"), Vg=DAG.add_var("Vg(z)"), Cd=DAG.add_var("Cd(rd,z)");
  FFVar Cw=DAG.add_var("Cw(rwc,z)"), Sw=DAG.add_var("Sw(rwc,z)");
  FFVar Cl=DAG.add_var("Cl(rl,z)"),  Sl=DAG.add_var("Sl(rl,z)");
  FFVar lam=DAG.add_var("lambda");                                // reaction fraction (continuation)
  FFPartial OpP;

  auto cyl=[&]( FFVar const& C, FFVar const& xi, double r0, double dr ){
    return OpP(C,{xi,2}) + ( dr/( r0 + xi*dr ) )*OpP(C,xi); };

  // gas plug-flow with real-gas Z (Eqs 2-3), corrected bracket (1 - y_CO2)=(1 - yin*Cg)
  FFVar GAS_C = vgIn*Vg*OpP(Cg,z) + cgf*( 1.0 - yin*Cg )*OpP(Cd,rd);
  FFVar GAS_V = OpP(Vg,z)         + cgfV*OpP(Cd,rd);
  FFVar DRY_C = cyl(Cd,rd,r1s,Drds);
  FFVar WET_C = cyl(Cw,rwc,rws,Drws) - kWC*lam*Cw*Sw;
  FFVar WET_S = cyl(Sw,rwc,rws,Drws) - kWS*lam*Cw*Sw;
  FFVar LIQ_C = vlbar*OpP(Cl,z) - DlCs*( OpP(Cl,{rl,2})/(Drls*Drls) + (1.0/(r2s+rl*Drls))*OpP(Cl,rl)/Drls ) + kLC*lam*Cl*Sl;
  FFVar LIQ_S = vlbar*OpP(Sl,z) - DlSs*( OpP(Sl,{rl,2})/(Drls*Drls) + (1.0/(r2s+rl*Drls))*OpP(Sl,rl)/Drls ) + kLS*lam*Cl*Sl;

  FFVar GAS_C_IN = Cg - 1.0,  GAS_V_IN = Vg - 1.0;                // scaled gas inlet (z=L)
  FFVar GD_VAL   = Cd - Cg;
  FFVar DW_VAL   = Cw - partH*Cd;
  FFVar DW_FLUX  = OpP(Cd,rd) - rDWf*OpP(Cw,rwc);                 // row-scaled
  FFVar DW_SFLUX = OpP(Sw,rwc);
  FFVar WL_CVAL  = Cw - Cl;
  FFVar WL_CFLUX = rWLf*OpP(Cw,rwc) - OpP(Cl,rl);                 // row-scaled
  FFVar WL_SVAL  = Sw - Sl;
  FFVar WL_SFLUX = rWLf*OpP(Sw,rwc) - OpP(Sl,rl);
  FFVar L_C_WALL = OpP(Cl,rl), L_S_WALL = OpP(Sl,rl);
  FFVar L_C_IN   = Cl - ClInS, L_S_IN = Sl - 1.0;                 // scaled liquid feed (z=0)

  OCFESLV oc( &DAG );
  oc.add_domain( z,   FFDom( 0., L,   3, FFDom::CGL, 6 ) );
  oc.add_domain( rd,  FFDom( 0., 1.0, 3, FFDom::CGL, 7 ) );
  oc.add_domain( rwc, FFDom( 0., 1.0, 3, FFDom::CGL, 8 ) );
  oc.add_domain( rl,  FFDom( 0., 1.0, 5, FFDom::CGL, 8 ) );   // 5 elements: resolve liquid depletion front
  oc.add_state( Cg,{z} );      oc.add_state( Vg,{z} );      oc.add_state( Cd,{rd,z} );
  oc.add_state( Cw,{rwc,z} );  oc.add_state( Sw,{rwc,z} );
  oc.add_state( Cl,{rl,z} );   oc.add_state( Sl,{rl,z} );
  oc.add_input( lam, {} );                                        // continuation parameter
  oc.update_ref( Cg,[&](OCFESLV::t_Coord const&){ return 1.0;    } );
  oc.update_ref( Vg,[&](OCFESLV::t_Coord const&){ return 1.0;    } );
  oc.update_ref( Cd,[&](OCFESLV::t_Coord const&){ return 1.0;    } );
  oc.update_ref( Cw,[&](OCFESLV::t_Coord const&){ return partH;  } );
  oc.update_ref( Sw,[&](OCFESLV::t_Coord const&){ return 1.0;    } );
  oc.update_ref( Cl,[&](OCFESLV::t_Coord const&){ return ClInS;  } );
  oc.update_ref( Sl,[&](OCFESLV::t_Coord const&){ return 1.0;    } );

  int const Z_NO_UB = FFDom::ALL - FFDom::UB;
  int const Z_NO_LB = FFDom::ALL - FFDom::LB;
  int const R_INT   = FFDom::ALL - FFDom::LB - FFDom::UB;
  OCFESLV::EqnOptions blk(OCFESLV::EqnRole::INTERIOR,0), bnd(OCFESLV::EqnRole::BOUNDARY,0), itf(OCFESLV::EqnRole::INTERFACE,0);

  oc.add_equation( GAS_V,    {z,rd},    {Z_NO_UB, FFDom::LB},                 blk );
  oc.add_equation( GAS_C,    {z,rd},    {Z_NO_UB, FFDom::LB},                 blk );
  oc.add_equation( GAS_V_IN, {z},       {FFDom::UB},                         bnd );
  oc.add_equation( GAS_C_IN, {z},       {FFDom::UB},                         bnd );
  oc.add_equation( DRY_C,    {rd,z},    {R_INT, FFDom::ALL},                 blk );
  oc.add_equation( WET_C,    {rwc,z},   {R_INT, FFDom::ALL},                 blk );
  oc.add_equation( WET_S,    {rwc,z},   {R_INT, FFDom::ALL},                 blk );
  oc.add_equation( LIQ_C,    {rl,z},    {R_INT, Z_NO_LB},                    blk );
  oc.add_equation( LIQ_S,    {rl,z},    {R_INT, Z_NO_LB},                    blk );
  oc.add_equation( GD_VAL,   {rd,z},    {FFDom::LB, FFDom::ALL},             itf );
  oc.add_equation( DW_VAL,   {rd,rwc,z},{FFDom::UB, FFDom::LB, FFDom::ALL},  itf );
  oc.add_equation( DW_FLUX,  {rd,rwc,z},{FFDom::UB, FFDom::LB, FFDom::ALL},  itf );
  oc.add_equation( DW_SFLUX, {rwc,z},   {FFDom::LB, FFDom::ALL},             itf );
  oc.add_equation( WL_CVAL,  {rwc,rl,z},{FFDom::UB, FFDom::LB, FFDom::ALL},  itf );
  oc.add_equation( WL_CFLUX, {rwc,rl,z},{FFDom::UB, FFDom::LB, FFDom::ALL},  itf );
  oc.add_equation( WL_SVAL,  {rwc,rl,z},{FFDom::UB, FFDom::LB, FFDom::ALL},  itf );
  oc.add_equation( WL_SFLUX, {rwc,rl,z},{FFDom::UB, FFDom::LB, FFDom::ALL},  itf );
  oc.add_equation( L_C_WALL, {rl,z},    {FFDom::UB, FFDom::ALL},             bnd );
  oc.add_equation( L_S_WALL, {rl,z},    {FFDom::UB, FFDom::ALL},             bnd );
  oc.add_equation( L_C_IN,   {rl,z},    {R_INT, FFDom::LB},                  bnd );
  oc.add_equation( L_S_IN,   {rl,z},    {R_INT, FFDom::LB},                  bnd );

  oc.options.REDUCE.ORDER     = OCFESLV::Options::RED_FULL;
  oc.options.CLASSIFY.MODE         = OCFESLV::Options::CLASS_AUTO;
  oc.options.INTERFACE.TYPE   = OCFESLV::Options::IC_AUTO;
  oc.options.INTERFACE.IMPOSITION  = OCFESLV::Options::IC_STRONG;
  oc.options.INTERFACE.SAT_SIGMA0       = 10.0;
  oc.options.DISPLAY_LEVEL    = 1;                                // SOLVE_RES_TOL left at default 1e-9
#if defined(CRONOS__WITH_SPQR)
  oc.options.SOLVE.FACTORIZATION = OCFESLV::Options::SOLVE_SPQR;
#else
  oc.options.SOLVE.FACTORIZATION = OCFESLV::Options::SOLVE_SUPERLU;
#endif

  if( !oc.setup() ){ std::cerr << "  setup() FAILED\n"; return 1; }
  check_true( "setup() completed (dimensionless multidomain)", true );
  std::vector<double> xv, inp;
  if( !oc.init( xv, inp, nullptr ) ){ std::cerr << "  init() FAILED\n"; return 1; }
  double const* ip = inp.empty()?nullptr:inp.data();

  // Full reaction solves DIRECTLY once the model is properly scaled -- no continuation
  // needed at this (mild) Da.  The lam input + set_input_values machinery is retained for
  // 2b, whose real Arrhenius kinetics (Da~1e3) will require reaction-rate continuation.
  oc.set_input_values( lam, { 1.0 }, inp.data() );
  OCFESLV::SolveReport const rep = oc.solve( xv.data(), inp.data(), nullptr );
  std::cout << "  solve: converged=" << (rep.converged?"yes":"no")
            << " iters=" << rep.iterations
            << " |r|=" << std::scientific << std::setprecision(3) << rep.final_residual << "\n";

  {
    std::vector<double> reqn( oc.n_colloc_eqn(), 0. ), rfct( oc.n_colloc_fct(), 0. );
    oc.eval( reqn.data(), rfct.data(), xv.data(), ip, nullptr );
    std::vector<size_t> ix( reqn.size() ); for( size_t i=0;i<ix.size();++i ) ix[i]=i;
    size_t const K = std::min<size_t>( 8, ix.size() );
    std::partial_sort( ix.begin(), ix.begin()+K, ix.end(), [&]( size_t a, size_t b ){ return std::fabs(reqn[a])>std::fabs(reqn[b]); } );
    std::cout << "  top residual rows: ";
    for( size_t i=0;i<K;++i ) std::cout << ix[i] << "=" << std::scientific << std::setprecision(2) << std::fabs(reqn[ix[i]]) << " ";
    std::cout << "\n";
  }
  check_true( "dimensionless MBC solve converged (conditioning)", rep.converged );
  if( !rep.converged ){ std::cout << "\n  RESULT: " << g_pass << " passed, " << g_fail << " failed\n"; return 1; }

  auto ev=[&](FFVar const& V, OCFESLV::t_Coord const& pt){ return oc.eval_colloc<double>(V,pt,xv.data(),ip,nullptr); };
  OCFESLV::t_Coord gi; gi[z]=L;   double Cg_in =ev(Cg,gi);
  OCFESLV::t_Coord go; go[z]=0.0; double Cg_out=ev(Cg,go);
  OCFESLV::t_Coord si; si[rl]=0.5; si[z]=0.0; double Sl_in =ev(Sl,si);
  OCFESLV::t_Coord so; so[rl]=0.5; so[z]=L;   double Sl_out=ev(Sl,so);
  double const removal = (Cg_in>1e-12)?(1.0-Cg_out/Cg_in):0.0;
  std::cout << std::fixed << std::setprecision(4)
            << "\n  gas CO2 (scaled): inlet " << Cg_in << "  outlet " << Cg_out << "\n"
            << "  amine   (scaled): inlet " << Sl_in << "  outlet " << Sl_out << "\n"
            << "  CO2 removal efficiency = " << std::setprecision(3) << 100.0*removal << " %\n";
  check_true( "gas CO2 decreases (absorption)",  Cg_out < Cg_in - 1e-4 );
  check_true( "amine consumed (Sl decreases)",   Sl_out < Sl_in - 1e-6 );
  check_true( "removal in (0,1)",                removal>1e-5 && removal<1.0 );
  check_true( "all fields finite",               std::isfinite(Cg_out)&&std::isfinite(Sl_out) );

  std::cout << "\n============================================================\n"
            << "  RESULT: " << g_pass << " passed, " << g_fail << " failed  -- "
            << (g_fail==0?"ALL PASS":"SOME FAILED") << "\n"
            << "============================================================\n";
  return g_fail==0?0:1;
}
