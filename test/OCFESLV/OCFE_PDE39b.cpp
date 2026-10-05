// OCFE_PDE39.cpp  ---  PSA INDEX-2, consistent-IC + re-pivot formulation (regression)
// ===========================================================================
// The split-IC counterpart of PDE39: a separate INITIAL IC on c plus re-pivoted coverage.  This
// formulation DE-INDEXES the index-2 system -- the root is a nonsingular fixed point, confirming
// the _reduce_high_index diagnosis.
//
// The symbolic reduction is identical across all four {CGL,LGL}x{IC_WEAK,IC_STRONG} (same
// G_reduced, same masks), and the sweep shows every mode recovers the root cleanly and UNIFORMLY
// up to eps=1e-2 (errors ~1e-13..1e-11, 4-5 iters); the finite Newton basin's edge is ~0.05 for
// all modes alike.  (The single-shot 0.05 that once suggested a "near-null CGL" limitation simply
// sat at that uniform basin edge; there is no basis- or imposition-specific conditioning gap.
// PDE39 is the same problem via the IC-synthesis path.)
//
// Formulation:
//   IC_c  {t=LB,  z=ALL}    INITIAL    pins c at t=0
//   CONT  {t>LB,  z=ALL}    INTERIOR   pins c at t>0
//   BC_u  {t=ALL, z=LB}     BOUNDARY   pins u at z=0
//   EOS   {t=ALL, z>LB}     INTERIOR -> G_reduced pins u at z>0
//
// REGRESSION PASS: for every mode -- structure sound (index==2, reduced) and the root recovered
// (err<1e-7) from every perturbed seed with eps <= 1e-2 (5e-2 is the informational basin edge).
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


static double const Ac = 0.3, Cc = 0.2;
static double const U0 = 1.0, Uz = 0.2, Ut = 0.15;

static inline double cM_exact( double z, double t ){ return 1.0 + Ac*(1.0-z)*(1.0-z) + Cc*t; }
static inline double uM_exact( double z, double t ){ return U0 + Uz*z + Ut*t; }

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
  bool   setup_ok=false, reduced=false, audit_ok=false, converged=false;
  int    index=-99, iters=0;
  size_t nTrace=0;
  double A_resid=std::nan(""), final_r=std::nan(""), err_c=std::nan(""), err_u=std::nan("");
};

static Outcome run_split( FFDom::TYPE coltype, OCFESLV::Options::ImpositionType imposition,
                          double eps, bool verbose )
{
  // Quiet runs rely on DISPLAY_LEVEL=0 + SOLVE_VERBOSE=false (set below) to silence the
  // framework's setup/solve diagnostics; the opt-in [cond] conditioning probe still prints.
  Outcome R;
  size_t const n_el = 3, n_nd = 6;

  FFGraph DAG;
  FFVar t = DAG.add_var( "t" );
  FFVar z = DAG.add_var( "z" );
  FFVar c = DAG.add_var( "c(t,z)" );
  FFVar u = DAG.add_var( "u(t,z)" );
  FFPartial OpP;

  FFVar c_Man = 1.0 + Ac*(1.0-z)*(1.0-z) + Cc*t;
  FFVar u_Man = U0 + Uz*z + Ut*t;
  FFVar dt_cM = Cc;
  FFVar dz_cM = -2.0*Ac*(1.0-z);
  FFVar dz_uM = Uz;
  FFVar div_uc_M = dz_uM*c_Man + u_Man*dz_cM;
  FFVar s_c = dt_cM + div_uc_M;

  FFVar CONT = OpP( c, t ) + OpP( u*c, z ) - s_c;
  FFVar EOS  = c - c_Man;
  FFVar IC_c = c - c_Man;
  FFVar BC_u = u - ( U0 + Ut*t );

  OCFESLV oc( &DAG );
  oc.add_domain( t, FFDom( 0., 1., n_el, coltype, n_nd ) );
  oc.add_domain( z, FFDom( 0., 1., n_el, coltype, n_nd ) );
  oc.add_state ( c, {t,z} );
  oc.add_state ( u, {t,z} );
  oc.update_ref( c, [&]( OCFESLV::t_Coord const& cr ){ return cM_exact( cr.at(z), cr.at(t) ); } );
  oc.update_ref( u, [&]( OCFESLV::t_Coord const& cr ){ return uM_exact( cr.at(z), cr.at(t) ); } );
  oc.set_evolution_domain( t );

  int const T_NO_LB = FFDom::ALL - FFDom::LB;
  int const Z_NO_LB = FFDom::ALL - FFDom::LB;
  oc.add_equation( IC_c, {t,z}, {FFDom::LB, FFDom::ALL}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INITIAL,  0 ) );
  oc.add_equation( CONT, {t,z}, {T_NO_LB,   FFDom::ALL}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( BC_u, {t,z}, {FFDom::ALL, FFDom::LB},  OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );
  oc.add_equation( EOS,  {t,z}, {FFDom::ALL, Z_NO_LB},    OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );

  oc.options.REDUCE.ORDER     = OCFESLV::Options::RED_FULL;
  oc.options.CLASSIFY.MODE         = OCFESLV::Options::CLASS_AUTO;
  oc.options.INTERFACE.TYPE   = OCFESLV::Options::IC_AUTO;
  oc.options.INTERFACE.IMPOSITION  = imposition;
  oc.options.INTERFACE.SAT_SIGMA0       = 10.0;
  oc.options.DISPLAY_LEVEL    = verbose ? 1 : 0;
  oc.options.SOLVE.VERBOSE    = verbose;
  oc.options.SOLVE.MAX_ITER   = 60;
  oc.options.FATAL.REDUCED_DOF = false;
#if defined(CRONOS__WITH_SPQR)
  oc.options.SOLVE.FACTORIZATION = OCFESLV::Options::SOLVE_SPQR;
#else
  oc.options.SOLVE.FACTORIZATION = OCFESLV::Options::SOLVE_SUPERLU;
#endif

  try { R.setup_ok = oc.setup(); } catch( ... ) {}
  R.index   = oc.pde_type().differential_index;
  R.reduced = !oc.reduction_plan().empty();
  R.audit_ok = oc.reduced_dof_audit().ok();
  if( !R.setup_ok ) return R;
  R.nTrace = oc.n_colloc_trace();

  size_t const nEqn=oc.n_colloc_eqn();
  std::vector<double> varInit, inpInit;
  if( !oc.init( varInit, inpInit, nullptr ) ) return R;

  std::vector<double> res( nEqn, 0. );
  if( oc.eval( res.data(), nullptr, varInit.data(), nullptr, nullptr ) ){
    R.A_resid=0.; for( double v : res ) R.A_resid = std::max( R.A_resid, std::fabs(v) );
  }

  std::vector<double> xv = varInit;
  for( size_t i=0;i<xv.size();++i ) xv[i] += eps*std::sin( 0.7*double(i) + 0.2 );
  OCFESLV::SolveReport rep = oc.solve( xv.data(), nullptr, nullptr );
  R.converged = rep.converged; R.iters = (int)rep.iterations; R.final_r = rep.final_residual;

  double ce=0., ue=0.;
  double const sg[4] = { 0.15, 0.35, 0.65, 0.85 };
  for( double zs : sg ) for( double ts : sg ){
    OCFESLV::t_Coord pt; pt[z]=zs; pt[t]=ts;
    ce = std::max( ce, std::fabs( oc.eval_colloc<double>( c, pt, xv.data(), nullptr, nullptr ) - cM_exact(zs,ts) ) );
    ue = std::max( ue, std::fabs( oc.eval_colloc<double>( u, pt, xv.data(), nullptr, nullptr ) - uM_exact(zs,ts) ) );
  }
  R.err_c=ce; R.err_u=ue;
  return R;
}

struct SweepResult { bool struct_ok=false; bool recovered_gated=true; };

// One mode: prints the perturbation sweep, captures structure (eps=0 row) and whether the root is
// recovered for every seed with eps <= GATE_EPS (all four modes are gated).
static SweepResult sweep( FFDom::TYPE coltype, std::string const& cname,
                          OCFESLV::Options::ImpositionType imposition )
{
  std::cout << "\n  sweep  basis=" << cname << "  imposition=" << imp_name(imposition) << "   [GATED]\n";
  std::cout << "    " << std::left << std::setw(12) << "eps"
            << std::setw(11) << "converged" << std::setw(7) << "iters"
            << std::setw(13) << "final|r|" << std::setw(13) << "err_c"
            << std::setw(13) << "err_u" << "verdict\n";
  double const GATE_EPS = 1e-2;   // all modes recover uniformly to 1e-2; 5e-2 is the informational basin edge
  double const epslist[] = { 0.0, 1e-10, 1e-8, 1e-6, 1e-4, 1e-2, 5e-2 };
  SweepResult SR;
  for( double eps : epslist ){
    Outcome R = run_split( coltype, imposition, eps, false );
    if( eps == 0.0 ) SR.struct_ok = ( R.index==2 && R.reduced );
    bool recovered = ( R.converged && R.err_c < 1e-7 && R.err_u < 1e-7 );
    bool ineps     = ( eps <= GATE_EPS );
    if( ineps && !recovered ) SR.recovered_gated = false;
    char const* verdict = recovered ? "ROOT recovered"
                        : ineps      ? "FAIL (gated)"
                        : R.converged ? "basin edge (wrong root)"
                                      : "basin edge (no-converge)";
    std::cout << "    " << std::left << std::setw(12) << std::scientific << std::setprecision(0) << eps
              << std::setw(11) << (R.converged?"yes":"no") << std::setw(7) << R.iters
              << std::setw(13) << std::setprecision(2) << R.final_r
              << std::setw(13) << R.err_c << std::setw(13) << R.err_u << verdict << "\n";
  }
  std::cout << "    -> structure " << (SR.struct_ok?"ok":"BAD")
            << " | root recovery (eps<=1e-2) " << (SR.recovered_gated?"ok":"no") << "\n";
  return SR;
}

int main()
{
  std::cout << "================================================================\n";
  std::cout << "  PSA INDEX-2 consistent-IC + re-pivot formulation (regression)\n";
  std::cout << "  de-indexes the index-2 system; every {CGL,LGL}x{WEAK,STRONG} mode\n";
  std::cout << "  recovers the exact root uniformly to eps=1e-2 (basin edge ~0.05).\n";
  std::cout << "================================================================\n";

  std::cout << "\n---- structural confirmation (CGL / IC_WEAK, eps=0) ----\n";
  Outcome S = run_split( FFDom::CGL, OCFESLV::Options::IC_WEAK, 0.0, false );
  std::cout << "  index=" << S.index << " reduced=" << (S.reduced?"yes":"no")
            << " audit_ok=" << (S.audit_ok?"yes":"no")
            << " [A]=" << std::scientific << std::setprecision(2) << S.A_resid << "\n";

  bool struct_all = true, recov_all = true;
  auto acc = [&]( SweepResult r ){ struct_all &= r.struct_ok; recov_all &= r.recovered_gated; };
  acc( sweep( FFDom::CGL, "CGL", OCFESLV::Options::IC_WEAK   ) );
  acc( sweep( FFDom::CGL, "CGL", OCFESLV::Options::IC_STRONG ) );
  acc( sweep( FFDom::LGL, "LGL", OCFESLV::Options::IC_WEAK   ) );
  acc( sweep( FFDom::LGL, "LGL", OCFESLV::Options::IC_STRONG ) );

  bool pass = struct_all && recov_all;
  std::cout << "\n  structure sound (index-2, reduced) all modes : " << (struct_all?"yes":"NO")
            << "\n  root recovered all modes (eps<=1e-2)         : " << (recov_all?"yes":"NO")
            << "\n  Overall: " << ( pass ? "PASS" : "FAIL" )
            << "  (uniform basin across all four modes; 5e-2 is the informational basin edge)\n";
  return pass ? 0 : 1;
}
