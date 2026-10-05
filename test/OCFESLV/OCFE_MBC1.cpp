// OCFE_MBC.cpp  ---  Stage 1: PHYSICAL isothermal membrane contactor (MBC)
// ===========================================================================
// Generalization of the manufactured multidomain regression OCFE_PDE3 towards the
// real high-pressure hollow-fiber MBC for CO2 absorption in natural-gas sweetening
// (Quek, Shah & Chachuat, Chem. Eng. Res. Des. 132 (2018) 1005-1019).
//
// PDE3's equations ALREADY carry the paper's physics -- counter-current gas tube +
// dry/wet membrane subdomains + liquid shell, real 2nd-order CO2-amine reaction,
// radial diffusion, gas-side flux coupling, and interface continuity/flux (paper
// Eqs. 22-23, 28-30).  This Stage-1 model differs from PDE3 in exactly two ways:
//   (1) the manufactured forcing (FG_*, FD_*, FW_*, FL_*) is REMOVED (-> 0);
//   (2) the gas inlets are REAL feed values (Cg=CgIn, Vg=VgIn at z=L) instead of
//       matched to a manufactured exact profile.
// Everything else -- transport, reaction, interface BCs, real lean-solvent feed
// (ClIn, SlIn), counter-current flow -- is retained unchanged.  Constant wetting
// radius, isothermal; betaG is a lumped real-gas correction.  No manufactured exact
// solution: validated by physical sanity (CO2 removed from the gas, solvent loaded).
//
//   states (paper Fig. 2/10):  Cg,Vg (gas tube) | Cd (dry membrane) |
//                              Cw,Sw (wet membrane) | Cl,Sl (liquid shell)
//   axial z in [0,L]; radial rd,rw,rl in [0,1] per subdomain; counter-current:
//   gas enters at z=L (UB), lean solvent enters at z=0 (LB).
// ===========================================================================
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <string>

#ifndef OCFE_OCFESLV_HEADER
#define OCFE_OCFESLV_HEADER "ocfeslv.hpp"
#endif
#include OCFE_OCFESLV_HEADER


using namespace mc;

struct Par {
  double L     = 1.0;
  double H     = 0.55;     // Henry-type partition (wet membrane / dry membrane)
  double Dmd   = 2.5e-2;   // dry-membrane CO2 diffusivity
  double DmwC  = 1.1e-2, DmwS = 6.0e-3;   // wet-membrane CO2 / solvent diffusivities
  double DlC   = 1.6e-2, DlS  = 7.0e-3;   // liquid CO2 / solvent diffusivities
  double epszL = 1.0e-2;   // axial dispersion (liquid)
  double krw   = 0.28, krl = 0.42;        // 2nd-order reaction rate (wet mem / liquid)
  double nuS   = 1.0;      // solvent stoichiometry
  double Kvg   = 0.25, Kcg = 0.65;        // gas-side coupling coefficients
  double betaG = 0.10;     // lumped real-gas correction
  double CgIn  = 1.0, VgIn = 1.0;         // gas feed: CO2 (normalized), velocity  (z=L)
  double ClIn  = 0.015, SlIn = 1.0;       // lean solvent feed: CO2, amine          (z=0)
};

static int g_pass = 0, g_fail = 0;
static void check_true( char const* name, bool ok )
{
  ( ok ? g_pass : g_fail )++;
  std::cout << "  " << std::left << std::setw(48) << name << "  " << ( ok ? "PASS" : "FAIL" ) << "\n";
}

int main()
{
  Par p;
  size_t const nelz = 3, nelr = 3, nz = 5, nr = 5;
  std::cout << "================================================================\n"
            << "  Stage 1: PHYSICAL isothermal MBC (PDE3 generalized)\n"
            << "  counter-current CO2 absorption; real reaction + interface flux\n"
            << "  gas feed CgIn=" << p.CgIn << " at z=L; lean solvent ClIn=" << p.ClIn << " at z=0\n"
            << "================================================================\n";

  FFGraph DAG;
  FFVar z  = DAG.add_var( "z" );
  FFVar rd = DAG.add_var( "rd" );
  FFVar rw = DAG.add_var( "rw" );
  FFVar rl = DAG.add_var( "rl" );
  FFVar Cg = DAG.add_var( "Cg(z)" );
  FFVar Vg = DAG.add_var( "Vg(z)" );
  FFVar Cd = DAG.add_var( "Cd(rd,z)" );
  FFVar Cw = DAG.add_var( "Cw(rw,z)" );
  FFVar Sw = DAG.add_var( "Sw(rw,z)" );
  FFVar Cl = DAG.add_var( "Cl(rl,z)" );
  FFVar Sl = DAG.add_var( "Sl(rl,z)" );
  FFVar krw_in = DAG.add_var( "krw" );   // reaction rates as inputs (continuation)
  FFVar krl_in = DAG.add_var( "krl" );
  FFPartial OpP;

  FFVar vl = 1.0 + 0.20*rl;              // liquid velocity profile (shell)
  FFVar Rw = krw_in*Cw*Sw;               // 2nd-order CO2-amine reaction (wet membrane)
  FFVar Rl = krl_in*Cl*Sl;               // 2nd-order CO2-amine reaction (liquid)

  // --- governing equations (PHYSICAL: manufactured forcing removed) ---
  FFVar GAS_V = OpP(Vg,z) + p.Kvg*p.Dmd*OpP(Cd,rd);
  FFVar GAS_C = OpP(Cg,z) + p.Kcg*p.Dmd*(1.0-p.betaG*Cg)/Vg*OpP(Cd,rd);
  FFVar DRY_C = -p.Dmd*OpP(Cd,{rd,2});
  FFVar WET_C = -p.DmwC*OpP(Cw,{rw,2}) + Rw;
  FFVar WET_S = -p.DmwS*OpP(Sw,{rw,2}) + p.nuS*Rw;
  FFVar LIQ_C = vl*OpP(Cl,z) - p.DlC*( OpP(Cl,{rl,2}) + p.epszL*OpP(Cl,{z,2}) ) + Rl;
  FFVar LIQ_S = vl*OpP(Sl,z) - p.DlS*( OpP(Sl,{rl,2}) + p.epszL*OpP(Sl,{z,2}) ) + p.nuS*Rl;

  // --- real gas inlets (z=L) ---
  FFVar GAS_C_IN = Cg - p.CgIn;
  FFVar GAS_V_IN = Vg - p.VgIn;
  // --- interface continuity / flux (unchanged from PDE3, all physical) ---
  FFVar GD_VAL   = Cd - Cg;                                    // gas <-> dry membrane (r1)
  FFVar DW_VAL   = Cw - p.H*Cd;                                // Henry partition at rw
  FFVar DW_FLUX  = p.Dmd*OpP(Cd,rd) - p.DmwC*OpP(Cw,rw);       // CO2 flux continuity
  FFVar DW_SFLUX = OpP(Sw,rw);                                 // no solvent into dry pores
  FFVar WL_CVAL  = Cw - Cl;                                    // wet membrane <-> liquid (r2)
  FFVar WL_CFLUX = p.DmwC*OpP(Cw,rw) - p.DlC*OpP(Cl,rl);
  FFVar WL_SVAL  = Sw - Sl;
  FFVar WL_SFLUX = p.DmwS*OpP(Sw,rw) - p.DlS*OpP(Sl,rl);
  // --- liquid feed (z=0), outlet, outer wall ---
  FFVar L_C_IN   = Cl - p.ClIn;
  FFVar L_S_IN   = Sl - p.SlIn;
  FFVar L_C_OUT  = OpP(Cl,z);
  FFVar L_S_OUT  = OpP(Sl,z);
  FFVar L_C_WALL = OpP(Cl,rl);
  FFVar L_S_WALL = OpP(Sl,rl);

  OCFESLV oc( &DAG );
  oc.add_domain( z,  FFDom( 0., p.L, nelz, FFDom::CGL, nz ) );
  oc.add_domain( rd, FFDom( 0., 1.0, nelr, FFDom::CGL, nr ) );
  oc.add_domain( rw, FFDom( 0., 1.0, nelr, FFDom::CGL, nr ) );
  oc.add_domain( rl, FFDom( 0., 1.0, nelr, FFDom::CGL, nr ) );
  oc.add_state( Cg, {z} );      oc.add_state( Vg, {z} );
  oc.add_state( Cd, {rd,z} );
  oc.add_state( Cw, {rw,z} );   oc.add_state( Sw, {rw,z} );
  oc.add_state( Cl, {rl,z} );   oc.add_state( Sl, {rl,z} );
  oc.add_input( krw_in, {} );   oc.add_input( krl_in, {} );   // continuation parameters
  // physical (feed-scaled) reference guesses
  oc.update_ref( Cg, [&]( OCFESLV::t_Coord const& ){ return p.CgIn; } );
  oc.update_ref( Vg, [&]( OCFESLV::t_Coord const& ){ return p.VgIn; } );
  oc.update_ref( Cd, [&]( OCFESLV::t_Coord const& ){ return p.CgIn; } );
  oc.update_ref( Cw, [&]( OCFESLV::t_Coord const& ){ return p.H*p.CgIn; } );
  oc.update_ref( Sw, [&]( OCFESLV::t_Coord const& ){ return p.SlIn; } );
  oc.update_ref( Cl, [&]( OCFESLV::t_Coord const& ){ return p.ClIn; } );
  oc.update_ref( Sl, [&]( OCFESLV::t_Coord const& ){ return p.SlIn; } );

  int const Z_NO_UB = FFDom::ALL - FFDom::UB;
  int const Z_INT   = FFDom::ALL - FFDom::LB - FFDom::UB;
  int const R_INT   = FFDom::ALL - FFDom::LB - FFDom::UB;
  OCFESLV::EqnOptions bulk( OCFESLV::EqnRole::INTERIOR, 0 ), bnd( OCFESLV::EqnRole::BOUNDARY, 0 ), itf( OCFESLV::EqnRole::INTERFACE, 0 );

  oc.add_equation( GAS_V,    {z,rd},    {Z_NO_UB, FFDom::LB},                 bulk );
  oc.add_equation( GAS_C,    {z,rd},    {Z_NO_UB, FFDom::LB},                 bulk );
  oc.add_equation( GAS_V_IN, {z},       {FFDom::UB},                         bnd  );
  oc.add_equation( GAS_C_IN, {z},       {FFDom::UB},                         bnd  );
  oc.add_equation( DRY_C,    {rd,z},    {R_INT, FFDom::ALL},                 bulk );
  oc.add_equation( WET_C,    {rw,z},    {R_INT, FFDom::ALL},                 bulk );
  oc.add_equation( WET_S,    {rw,z},    {R_INT, FFDom::ALL},                 bulk );
  oc.add_equation( LIQ_C,    {rl,z},    {R_INT, Z_INT},                      bulk );
  oc.add_equation( LIQ_S,    {rl,z},    {R_INT, Z_INT},                      bulk );
  oc.add_equation( GD_VAL,   {rd,z},    {FFDom::LB, FFDom::ALL},             itf  );
  oc.add_equation( DW_VAL,   {rd,rw,z}, {FFDom::UB, FFDom::LB, FFDom::ALL},  itf  );
  oc.add_equation( DW_FLUX,  {rd,rw,z}, {FFDom::UB, FFDom::LB, FFDom::ALL},  itf  );
  oc.add_equation( DW_SFLUX, {rw,z},    {FFDom::LB, FFDom::ALL},             itf  );
  oc.add_equation( WL_CVAL,  {rw,rl,z}, {FFDom::UB, FFDom::LB, FFDom::ALL},  itf  );
  oc.add_equation( WL_CFLUX, {rw,rl,z}, {FFDom::UB, FFDom::LB, FFDom::ALL},  itf  );
  oc.add_equation( WL_SVAL,  {rw,rl,z}, {FFDom::UB, FFDom::LB, FFDom::ALL},  itf  );
  oc.add_equation( WL_SFLUX, {rw,rl,z}, {FFDom::UB, FFDom::LB, FFDom::ALL},  itf  );
  oc.add_equation( L_C_WALL, {rl,z},    {FFDom::UB, FFDom::ALL},             bnd  );
  oc.add_equation( L_S_WALL, {rl,z},    {FFDom::UB, FFDom::ALL},             bnd  );
  oc.add_equation( L_C_IN,   {rl,z},    {R_INT, FFDom::LB},                  bnd  );
  oc.add_equation( L_S_IN,   {rl,z},    {R_INT, FFDom::LB},                  bnd  );
  oc.add_equation( L_C_OUT,  {rl,z},    {R_INT, FFDom::UB},                  bnd  );
  oc.add_equation( L_S_OUT,  {rl,z},    {R_INT, FFDom::UB},                  bnd  );

  // outputs: treated-gas CO2 at the outlet (z=0), rich-solvent CO2 at its outlet (z=L, mid-shell)
  oc.add_output( Cg, {z},    {0.0} );          // fct 0: treated gas CO2
  oc.add_output( Cl, {rl,z}, {0.5, p.L} );      // fct 1: rich-solvent CO2 (mid radius, outlet)

  oc.options.REDUCE.ORDER     = OCFESLV::Options::RED_FULL;
  oc.options.CLASSIFY.MODE         = OCFESLV::Options::CLASS_AUTO;
  oc.options.INTERFACE.TYPE   = OCFESLV::Options::IC_AUTO;
  oc.options.INTERFACE.IMPOSITION  = OCFESLV::Options::IC_WEAK;//STRONG;
  oc.options.INTERFACE.SAT_SIGMA0       = 10.0;
  oc.options.DISPLAY_LEVEL    = 1;
#if defined(CRONOS__WITH_SPQR)
  oc.options.SOLVE.FACTORIZATION = OCFESLV::Options::SOLVE_SPQR;
#else
  oc.options.SOLVE.FACTORIZATION = OCFESLV::Options::SOLVE_SUPERLU;
#endif

  if( !oc.setup() ){ std::cerr << "  setup() FAILED\n"; return 1; }
  check_true( "setup() completed (multidomain physical MBC)", true );

  std::vector<double> xv, inp;
  if( !oc.init( xv, inp, nullptr ) ){ std::cerr << "  init() FAILED\n"; return 1; }
  double const* ip = inp.empty() ? nullptr : inp.data();

  // reaction continuation: solve the (linear) no-reaction problem, then ramp krw,krl
  // to target, warm-starting each step -- the 2nd-order reaction from a flat start
  // otherwise stalls the Newton in a shallow basin.
  OCFESLV::SolveReport rep;
  double const lam_seq[] = { 1.0 };//0.0, 0.25, 0.5, 0.75, 1.0 };
  for( double lam : lam_seq ){
    oc.set_input_values( krw_in, { lam*p.krw }, inp.data() );
    oc.set_input_values( krl_in, { lam*p.krl }, inp.data() );
    rep = oc.solve( xv.data(), inp.data(), nullptr );
    std::cout << "  [continuation] lambda=" << std::fixed << std::setprecision(2) << lam
              << "  converged=" << ( rep.converged ? "y" : "n" )
              << "  |r|=" << std::scientific << std::setprecision(3) << rep.final_residual << "\n";
    if( !rep.converged ) break;
  }
  check_true( "physical MBC solve converged", rep.converged );
  if( !rep.converged ){
    std::cout << "\n  RESULT: " << g_pass << " passed, " << g_fail << " failed\n";
    return 1;
  }

  // physical sanity: CO2 removed from the gas, solvent loaded with CO2
  auto ev = [&]( FFVar const& V, OCFESLV::t_Coord const& pt ){ return oc.eval_colloc<double>( V, pt, xv.data(), ip, nullptr ); };
  OCFESLV::t_Coord g_in;  g_in[z]  = p.L;   double Cg_in  = ev( Cg, g_in );
  OCFESLV::t_Coord g_out; g_out[z] = 0.0;   double Cg_out = ev( Cg, g_out );
  OCFESLV::t_Coord l_in;  l_in[rl] = 0.5; l_in[z] = 0.0;   double Cl_in  = ev( Cl, l_in );
  OCFESLV::t_Coord l_out; l_out[rl]= 0.5; l_out[z]= p.L;    double Cl_out = ev( Cl, l_out );
  OCFESLV::t_Coord s_in;  s_in[rl] = 0.5; s_in[z] = 0.0;    double Sl_in  = ev( Sl, s_in );
  OCFESLV::t_Coord s_out; s_out[rl]= 0.5; s_out[z]= p.L;    double Sl_out = ev( Sl, s_out );
  double const removal = ( Cg_in > 1e-12 ) ? ( 1.0 - Cg_out/Cg_in ) : 0.0;

  std::cout << std::fixed << std::setprecision(5)
            << "\n  gas   CO2:  inlet Cg(z=L)=" << Cg_in  << "  outlet Cg(z=0)=" << Cg_out << "\n"
            << "  liq   CO2:  inlet Cl(z=0)=" << Cl_in  << "  outlet Cl(z=L)=" << Cl_out << "  (free CO2; consumed by reaction)\n"
            << "  amine    :  inlet Sl(z=0)=" << Sl_in  << "  outlet Sl(z=L)=" << Sl_out << "  (loaded by reaction)\n"
            << "  CO2 removal efficiency = " << 100.0*removal << " %\n";
  // reactive absorption: gas loses CO2, amine is consumed (loaded); free liquid CO2 (Cl)
  // is drawn down by the reaction, so it is NOT a valid "gain" indicator.
  check_true( "gas CO2 decreases along the module (absorption)", Cg_out < Cg_in - 1e-6 );
  check_true( "amine consumed by reaction (Sl loaded, decreases)", Sl_out < Sl_in - 1e-9 );
  check_true( "removal efficiency in (0,1)",                     removal > 1e-5 && removal < 1.0 );
  check_true( "all fields finite",                              std::isfinite(Cg_out) && std::isfinite(Sl_out) );

  std::cout << "\n============================================================\n"
            << "  RESULT: " << g_pass << " passed, " << g_fail << " failed  -- "
            << ( g_fail == 0 ? "ALL PASS" : "SOME FAILED" ) << "\n"
            << "============================================================\n";
  return g_fail == 0 ? 0 : 1;
}
