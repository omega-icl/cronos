// OCFE_PDE40_solve2.cpp  ---  PSA case study, RUNG 7 INDEX-2 / stage A1
//                             (MULTI-COMPONENT, constant-P total-continuity closure, MMS)
// ===========================================================================
// The multi-component analogue of the PDE39 index-2 reproducer, built as the
// first stage of the full physical index-2 bed that closes the PSA case study.
//
// It is EXACTLY rung 7b (PDE33) with the index-1 pressure-Darcy momentum closure
// replaced by the genuinely INDEX-2 total-continuity closure:
//
//   * DROP the pressure state P and the DARCY equation (u <- dP/dz).
//   * The EOS becomes a DERIVATIVE-FREE constraint on the TOTAL gas density:
//         EOS:  (c1 + c2) - rho(z,t) = 0        [ rho = P/(Rg T); prescribed ]
//     so ONE combination of the two differentiated components is algebraically
//     pinned while c1,c2 both appear differentiated in their balances.
//   * u is pinned by the total continuity that EMERGES from summing the two
//     component balances (velocity varies to accommodate the adsorption sink,
//     the constant-pressure PSA relation), with a single inlet BC_u.
//
// Structural Pantelides: EOS pins c1+c2 (derivative-free) but c1,c2 are
// differentiated in CONT1,CONT2 -> differentiate EOS once in t -> combined with
// CONT1+CONT2 pins d_z(u(c1+c2)) -> pins u.  INDEX 2, witness {u}.  Because the
// algebraic witness u carries NO explicit IC, the reducer must SYNTHESIZE u's
// consistent IC at the evolution-LB (t=0) -- the PDE39 mechanism, now on the
// multi-component structure.  This stage validates that de-indexing against a
// manufactured solution (machine-precision check) before the non-isothermal (A2)
// and physical (B) stages, which have no exact solution to check against.
//
// Everything OTHER than the momentum closure (competitive Langmuir, LDF kinetics,
// Danckwerts inlet, zero-gradient outlet, axial diffusion, MMS sources) is
// identical to the validated PDE33, so a failure here isolates the new
// multi-component EOS->u index-2 chain.
//
// Manufactured (identical fields/sources to PDE33; u kept as the PDE33 u_man so
// the CONT sources are byte-identical):
//   c1 = 1  + Ac1(1-z)^2 + Cc1 t,     c2 = C20 + Ac2(1-z)^2 + Cc2 t,
//   q1 = Q10 + Qz1 z + Qt1 t,         q2 = Q20 + Qz2 z + Qt2 t,
//   rho = c1 + c2,                    u  = 2 kappa_g (Ac1+Ac2)(1-z).
// Imposition sweep over {CGL,LGL} x {IC_WEAK,IC_STRONG}; PSA default IC_STRONG.
// ===========================================================================

#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <string>
#include <limits>

#ifndef OCFE_OCFESLV_HEADER
#define OCFE_OCFESLV_HEADER "ocfeslv.hpp"
#endif
#include OCFE_OCFESLV_HEADER

using namespace mc;

// transport / kinetics  (identical to PDE33)
static double const Kperm = 0.5, Rg = 1.0;
static double const D_ax = 0.1, F_ph = 0.5, k_ldf = 2.0;
static double const kappa_g = Kperm*Rg;                // u_man = 2 kappa_g (Ac1+Ac2)(1-z)
// competitive Langmuir
static double const qs1 = 1.0, qs2 = 1.0, b1_L = 1.0, b2_L = 0.5;
// manufactured coefficients -- component 1
static double const Ac1 = 0.3, Cc1 = 0.2;
static double const Q10 = 0.4, Qz1 = 0.1, Qt1 = 0.3;
// manufactured coefficients -- component 2
static double const C20 = 0.5, Ac2 = 0.2, Cc2 = 0.15;
static double const Q20 = 0.25, Qz2 = 0.08, Qt2 = 0.12;
static double const A_sum = Ac1 + Ac2;

static inline double c1M_exact( double z, double t ){ return 1.0 + Ac1*(1.0-z)*(1.0-z) + Cc1*t; }
static inline double c2M_exact( double z, double t ){ return C20 + Ac2*(1.0-z)*(1.0-z) + Cc2*t; }
static inline double q1M_exact( double z, double t ){ return Q10 + Qz1*z + Qt1*t; }
static inline double q2M_exact( double z, double t ){ return Q20 + Qz2*z + Qt2*t; }
static inline double uM_exact ( double z, double   ){ return 2.0*kappa_g*A_sum*(1.0-z); }

static char const* imp_name( OCFESLV::Options::ImpositionType it )
{
  switch( it ){
  case OCFESLV::Options::IC_WEAK:   return "IC_WEAK";
  case OCFESLV::Options::IC_TRACE:  return "IC_TRACE";
  case OCFESLV::Options::IC_STRONG: return "IC_STRONG";
  default:                        return "IC_?";
  }
}

struct Outcome
{
  bool   setup_ok=false, reduced=false, audit_ran=false, audit_ok=false, square=false;
  int    index=-99;
  size_t nVar=0, nEqn=0, nTrace=0;
  double A_resid = std::numeric_limits<double>::quiet_NaN();
  bool   converged=false;
  int    iters=0;
  double final_r=std::numeric_limits<double>::quiet_NaN();
  double err_c1=0, err_c2=0, err_q1=0, err_q2=0, err_u=0;
};

static Outcome run_index2mc( FFDom::TYPE coltype, OCFESLV::Options::ImpositionType imposition,
                             double eps, bool verbose )
{
  Outcome R;
  size_t const n_el = 3, n_nd = 5;   // n_nd=5 avoids the n_nd=6 CGL/IC_STRONG basin pocket (mapped by PDE40b); still MMS-exact

  FFGraph DAG;
  FFVar t  = DAG.add_var( "t" );
  FFVar z  = DAG.add_var( "z" );
  FFVar c1 = DAG.add_var( "c1(t,z)" );
  FFVar c2 = DAG.add_var( "c2(t,z)" );
  FFVar q1 = DAG.add_var( "q1(t,z)" );
  FFVar q2 = DAG.add_var( "q2(t,z)" );
  FFVar u  = DAG.add_var( "u(t,z)" );
  FFPartial OpP;

  // manufactured fields
  FFVar c1Man = 1.0 + Ac1*(1.0-z)*(1.0-z) + Cc1*t;
  FFVar c2Man = C20 + Ac2*(1.0-z)*(1.0-z) + Cc2*t;
  FFVar q1Man = Q10 + Qz1*z + Qt1*t;
  FFVar q2Man = Q20 + Qz2*z + Qt2*t;
  FFVar uMan  = 2.0*kappa_g*A_sum*(1.0-z);
  FFVar rhoMan= c1Man + c2Man;                          // prescribed total density (EOS RHS)

  FFVar denMan    = 1.0 + b1_L*c1Man + b2_L*c2Man;
  FFVar q1starMan = qs1*b1_L*c1Man/denMan;
  FFVar q2starMan = qs2*b2_L*c2Man/denMan;

  // MMS sources (identical to PDE33: u_man kept, so CONT sources are byte-identical)
  FFVar s_c1 = Cc1 - 2.0*kappa_g*A_sum*( 1.0 + 3.0*Ac1*(1.0-z)*(1.0-z) + Cc1*t ) - 2.0*D_ax*Ac1 + F_ph*Qt1;
  FFVar s_c2 = Cc2 - 2.0*kappa_g*A_sum*( C20 + 3.0*Ac2*(1.0-z)*(1.0-z) + Cc2*t ) - 2.0*D_ax*Ac2 + F_ph*Qt2;
  FFVar s_q1 = Qt1 - k_ldf*( q1starMan - q1Man );
  FFVar s_q2 = Qt2 - k_ldf*( q2starMan - q2Man );

  // Danckwerts inlet fluxes at z=0:  g_i = u_man(0,t) c_iMan(0,t) - D d_z c_iMan(0,t)
  //   c_iMan(0,t)=base_i+Ac_i+Cc_i t,  d_z c_iMan(0,t) = -2 Ac_i,  u_man(0,t)=2 kappa_g A_sum
  FFVar g1_in = 2.0*kappa_g*A_sum*( 1.0 + Ac1 + Cc1*t ) + 2.0*D_ax*Ac1;
  FFVar g2_in = 2.0*kappa_g*A_sum*( C20 + Ac2 + Cc2*t ) + 2.0*D_ax*Ac2;

  // live isotherm
  FFVar den    = 1.0 + b1_L*c1 + b2_L*c2;
  FFVar q1star = qs1*b1_L*c1/den;
  FFVar q2star = qs2*b2_L*c2/den;

  // equations
  FFVar CONT1 = OpP( c1, t ) + OpP( u*c1, z ) - D_ax*OpP( OpP( c1, z ), z ) + F_ph*OpP( q1, t ) - s_c1;
  FFVar CONT2 = OpP( c2, t ) + OpP( u*c2, z ) - D_ax*OpP( OpP( c2, z ), z ) + F_ph*OpP( q2, t ) - s_c2;
  FFVar LDF1  = OpP( q1, t ) - k_ldf*( q1star - q1 ) - s_q1;
  FFVar LDF2  = OpP( q2, t ) - k_ldf*( q2star - q2 ) - s_q2;
  FFVar EOS   = ( c1 + c2 ) - rhoMan;                   // derivative-free TOTAL-density constraint (index-2 trigger)
  FFVar IC_c1 = c1 - ( 1.0 + Ac1*(1.0-z)*(1.0-z) );
  FFVar IC_c2 = c2 - ( C20 + Ac2*(1.0-z)*(1.0-z) );
  FFVar IC_q1 = q1 - ( Q10 + Qz1*z );
  FFVar IC_q2 = q2 - ( Q20 + Qz2*z );
  FFVar BC_L1 = u*c1 - D_ax*OpP( c1, z ) - g1_in;       // Danckwerts inlet
  FFVar BC_L2 = u*c2 - D_ax*OpP( c2, z ) - g2_in;
  FFVar BC_U1 = OpP( c1, z );                           // zero-gradient outlet
  FFVar BC_U2 = OpP( c2, z );
  FFVar BC_u  = u - 2.0*kappa_g*A_sum;                  // inlet velocity (divergence integration constant)

  OCFESLV oc( &DAG );
  oc.add_domain( t, FFDom( 0., 1., n_el, coltype, n_nd ) );
  oc.add_domain( z, FFDom( 0., 1., n_el, coltype, n_nd ) );
  oc.add_state ( c1, {t,z} );
  oc.add_state ( c2, {t,z} );
  oc.add_state ( q1, {t,z} );
  oc.add_state ( q2, {t,z} );
  oc.add_state ( u,  {t,z} );
  oc.update_ref( c1, [&]( OCFESLV::t_Coord const& cr ){ return c1M_exact( cr.at(z), cr.at(t) ); } );
  oc.update_ref( c2, [&]( OCFESLV::t_Coord const& cr ){ return c2M_exact( cr.at(z), cr.at(t) ); } );
  oc.update_ref( q1, [&]( OCFESLV::t_Coord const& cr ){ return q1M_exact( cr.at(z), cr.at(t) ); } );
  oc.update_ref( q2, [&]( OCFESLV::t_Coord const& cr ){ return q2M_exact( cr.at(z), cr.at(t) ); } );
  oc.update_ref( u,  [&]( OCFESLV::t_Coord const& cr ){ return uM_exact ( cr.at(z), cr.at(t) ); } );
  oc.set_evolution_domain( t );

  int const T_NO_LB = FFDom::ALL - FFDom::LB;
  int const Z_NO_LB = FFDom::ALL - FFDom::LB;
  int const Z_INT   = FFDom::ALL - FFDom::LB - FFDom::UB;

  // CONT1/CONT2 pin the components + drive u across the interior; EOS pins the
  // total density off the inlet (the feed BCs pin the inlet); u's inlet BC spans
  // t=LB to pin the (t=0,z=0) corner while the reducer synthesizes u's IC on the
  // interior/UB of the t=0 face.
  oc.add_equation( CONT1, {t,z}, {T_NO_LB,   Z_INT},     OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( CONT2, {t,z}, {T_NO_LB,   Z_INT},     OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( LDF1,  {t,z}, {T_NO_LB,   FFDom::ALL},OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( LDF2,  {t,z}, {T_NO_LB,   FFDom::ALL},OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( EOS,   {t,z}, {FFDom::ALL,Z_NO_LB},   OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( IC_c1, {t,z}, {FFDom::LB, FFDom::ALL},OCFESLV::EqnOptions( OCFESLV::EqnRole::INITIAL,  0 ) );
  oc.add_equation( IC_c2, {t,z}, {FFDom::LB, FFDom::ALL},OCFESLV::EqnOptions( OCFESLV::EqnRole::INITIAL,  0 ) );
  oc.add_equation( IC_q1, {t,z}, {FFDom::LB, FFDom::ALL},OCFESLV::EqnOptions( OCFESLV::EqnRole::INITIAL,  0 ) );
  oc.add_equation( IC_q2, {t,z}, {FFDom::LB, FFDom::ALL},OCFESLV::EqnOptions( OCFESLV::EqnRole::INITIAL,  0 ) );
  oc.add_equation( BC_L1, {t,z}, {T_NO_LB,   FFDom::LB}, OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );
  oc.add_equation( BC_L2, {t,z}, {T_NO_LB,   FFDom::LB}, OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );
  oc.add_equation( BC_U1, {t,z}, {T_NO_LB,   FFDom::UB}, OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );
  oc.add_equation( BC_U2, {t,z}, {T_NO_LB,   FFDom::UB}, OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );
  oc.add_equation( BC_u,  {t,z}, {FFDom::ALL,FFDom::LB}, OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );

  oc.options.REDUCE.ORDER     = OCFESLV::Options::RED_FULL;
  oc.options.CLASSIFY.MODE         = OCFESLV::Options::CLASS_AUTO;
  oc.options.INTERFACE.TYPE   = OCFESLV::Options::IC_AUTO;
  oc.options.INTERFACE.IMPOSITION  = imposition;
  // This driver validates the index-2 reduction with FULL, CONSISTENT initial data (MMS-exact c1, c2, q1, q2): the
  // former convention.  Under REDUCE.HIDDEN_IC = true (the default since 2026-10-04) the EOS level at t = 0 -- imposed
  // off the inlet, N_z - 1 rows -- would be materialised on top of it: an over-determined initial point.  (The free
  // data in the new convention: c1 everywhere, c2 at the inlet node only, q1, q2.)
  oc.options.REDUCE.HIDDEN_IC = false;
  oc.options.INTERFACE.SAT_SIGMA0       = 10.0;
  oc.options.DISPLAY_LEVEL    = verbose ? 1 : 0;
  oc.options.SOLVE.VERBOSE    = false;
  oc.options.SOLVE.MAX_ITER   = 60;
  oc.options.SOLVE.RES_TOL    = 1e-9;
  oc.options.FATAL.REDUCED_DOF = false;           // advisory during bring-up of a NEW structure
#if defined(CRONOS__WITH_SPQR)
  oc.options.SOLVE.FACTORIZATION = OCFESLV::Options::SOLVE_SPQR;
#else
  oc.options.SOLVE.FACTORIZATION = OCFESLV::Options::SOLVE_SUPERLU;
#endif

  try { R.setup_ok = oc.setup(); }
  catch( ... ) { R.setup_ok = false; }

  R.index   = oc.pde_type().differential_index;
  R.reduced = !oc.reduction_plan().empty();
  { OCFESLV::t_DofAudit const& A = oc.reduced_dof_audit();
    R.audit_ran = A.ran; R.audit_ok = A.ok(); }
  if( !R.setup_ok ) return R;

  R.nVar = oc.n_colloc_sta(); R.nEqn = oc.n_colloc_eqn(); R.nTrace = oc.n_colloc_trace();
  R.square = ( R.nVar == R.nEqn );

  if( verbose ){
    std::cout << "  [struct] index=" << R.index
              << " character=" << OCFESLV::pde_type_name( oc.pde_type().type )
              << " reduced=" << (R.reduced?"yes":"no")
              << " | audit ran=" << (R.audit_ran?"yes":"no")
              << " ok=" << (R.audit_ok?"yes":"NO")
              << " | nVar=" << R.nVar << " nEqn=" << R.nEqn
              << " nTrace=" << R.nTrace << " square=" << (R.square?"yes":"NO") << "\n";
  }

  std::vector<double> varInit, inpInit;
  if( !oc.init( varInit, inpInit, nullptr ) ) return R;

  std::vector<double> res( R.nEqn, 0. );
  if( oc.eval( res.data(), nullptr, varInit.data(), nullptr, nullptr ) ){
    R.A_resid=0.; for( double v : res ) R.A_resid = std::max( R.A_resid, std::fabs(v) );
  }

  std::vector<double> xv = varInit;
  for( size_t i=0;i<xv.size();++i ) xv[i] += eps*std::sin( 0.7*double(i) + 0.2 );

  OCFESLV::SolveReport rep = oc.solve( xv.data(), nullptr, nullptr );
  R.converged = rep.converged; R.iters = (int)rep.iterations; R.final_r = rep.final_residual;

  double const sg[4] = { 0.15, 0.35, 0.65, 0.85 };
  for( double zs : sg ) for( double ts : sg ){
    OCFESLV::t_Coord pt; pt[z]=zs; pt[t]=ts;
    R.err_c1 = std::max( R.err_c1, std::fabs( oc.eval_colloc<double>( c1, pt, xv.data(), nullptr, nullptr ) - c1M_exact(zs,ts) ) );
    R.err_c2 = std::max( R.err_c2, std::fabs( oc.eval_colloc<double>( c2, pt, xv.data(), nullptr, nullptr ) - c2M_exact(zs,ts) ) );
    R.err_q1 = std::max( R.err_q1, std::fabs( oc.eval_colloc<double>( q1, pt, xv.data(), nullptr, nullptr ) - q1M_exact(zs,ts) ) );
    R.err_q2 = std::max( R.err_q2, std::fabs( oc.eval_colloc<double>( q2, pt, xv.data(), nullptr, nullptr ) - q2M_exact(zs,ts) ) );
    R.err_u  = std::max( R.err_u,  std::fabs( oc.eval_colloc<double>( u,  pt, xv.data(), nullptr, nullptr ) - uM_exact (zs,ts) ) );
  }
  return R;
}

// One mode: prints the row and returns PASS iff structure correct (index-2, reduced, square,
// MMS-exact) AND the root recovered from the perturbed seed.
static bool sweep_row( FFDom::TYPE coltype, std::string const& cname,
                       OCFESLV::Options::ImpositionType imposition )
{
  double const eps = 1e-4;   // modest perturbation; n_nd=5 keeps every mode inside the Newton basin
  Outcome R = run_index2mc( coltype, imposition, eps, /*verbose=*/false );
  double errmax = std::max( { R.err_c1, R.err_c2, R.err_q1, R.err_q2, R.err_u } );
  bool recovered = ( R.converged && errmax < 1e-7 );
  bool struct_ok = ( R.index==2 && R.reduced && R.square && R.A_resid < 1e-9 );
  bool pass = struct_ok && recovered;
  char const* verdict = !R.setup_ok  ? "SETUP-FAIL"
                      : !R.square     ? "NOT-SQUARE"
                      : !struct_ok    ? "STRUCT-FAIL"
                      : !R.converged  ? "no-converge"
                      : recovered     ? "PASS"
                                      : "wrong root";
  std::cout << "  " << std::left << std::setw(6) << cname
            << std::setw(11) << imp_name(imposition)
            << " idx=" << R.index << " red=" << (R.reduced?"y":"n")
            << " sq=" << (R.square?"y":"N")
            << " audit=" << (R.audit_ran?(R.audit_ok?"ok":"NO"):"-")
            << " | A|r|=" << std::scientific << std::setprecision(1) << R.A_resid
            << " conv=" << (R.converged?"y":"n") << " it=" << R.iters
            << " |r|=" << std::setprecision(1) << R.final_r
            << " errmax=" << std::setprecision(2) << errmax
            << "  " << verdict << "\n";
  return pass;
}

int main()
{
  std::cout << "================================================================\n";
  std::cout << "  PSA RUNG 7 INDEX-2 / stage A1 -- multi-component, constant-P\n";
  std::cout << "  total-continuity closure (MMS).  PDE33 with Darcy -> index-2.\n";
  std::cout << "================================================================\n";

  // verbose structural confirmation on the PSA default (CGL / IC_STRONG)
  std::cout << "\n---- structural confirmation (CGL / IC_STRONG) ----\n";
  Outcome S = run_index2mc( FFDom::CGL, OCFESLV::Options::IC_STRONG, /*eps=*/1e-4, /*verbose=*/true );
  std::cout << "  [A] residual at manufactured exact: " << std::scientific << std::setprecision(3)
            << S.A_resid << ( S.A_resid<1e-9 ? "  (consistent)" : "  (NOT exact!)" ) << "\n";
  std::cout << "  [B] solve: converged=" << (S.converged?"yes":"no") << " iters=" << S.iters
            << " final|r|=" << std::setprecision(3) << S.final_r << "\n";
  std::cout << "      recovery: err_c1=" << std::setprecision(2) << S.err_c1
            << " err_c2=" << S.err_c2 << " err_q1=" << S.err_q1
            << " err_q2=" << S.err_q2 << " err_u=" << S.err_u << "\n";

  // full mode sweep -- gate PASS on structure + root recovery in every mode
  std::cout << "\n---- mode sweep (structural + recovery) ----\n";
  bool pass = true;
  pass &= sweep_row( FFDom::CGL, "CGL", OCFESLV::Options::IC_WEAK   );
  pass &= sweep_row( FFDom::CGL, "CGL", OCFESLV::Options::IC_STRONG );
  pass &= sweep_row( FFDom::LGL, "LGL", OCFESLV::Options::IC_WEAK   );
  pass &= sweep_row( FFDom::LGL, "LGL", OCFESLV::Options::IC_STRONG );

  std::cout << "\n  Overall: " << ( pass ? "PASS" : "FAIL" )
            << "  (all modes: idx=2, red=y, sq=y, MMS-exact, root recovered)\n";
  return pass ? 0 : 1;
}
