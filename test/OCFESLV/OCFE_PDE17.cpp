// ===========================================================================
// OCFE_PDE17_solve2.cpp  --  reference-robustness probe oracle (sign flip)
// ===========================================================================
// Validates OCFESLV::_probe_reference_robustness (MC__OCFESLV_REFERENCE_ROBUSTNESS_
// PROBE).  All setup-time structural analysis linearizes ONCE at the single
// _classification_reference.  For a NONLINEAR/quasilinear A_z(u) that one point
// can get the STRUCTURE wrong, and -- unlike the solver's per-iteration Newton/LM
// re-linearization -- a wrong structure is baked into the assembled system at
// setup and cannot be recovered.  The probe re-evaluates the SAME principal symbol
// at extra state references and flags type / characteristic-split changes.
//
// This oracle builds a QUASILINEAR hyperbolic block whose wave speed is a STATE:
//
//     d_t c +     c d_z u = f_c           A_z(c) = c [[0,1],[g,0]] ,  g = 4
//     d_t u + g   c d_z c = f_u           eigenvalues  +- 2c ,  sg = sqrt(g) = 2
//
//   left eigvecs l_+- = (sg,+-1) ; characteristics w_+- = 2c +- u, speeds +- 2c.
//   The eigenVECTORS are c-INDEPENDENT (c only scales A_z), but the eigenVALUE
//   SIGNS follow sign(c): at c>0, w_+ (speed +2c) is incoming@LB; at c<0 every
//   speed reverses, so w_- becomes incoming@LB -- a PAIRED SWAP.  The in/out COUNT
//   per face is unchanged (1/1 either way); only WHICH characteristic is incoming
//   flips, rotating the incoming subspace Vin by ~53 deg.  This is the miniature
//   of a PSA flow reversal across cycle steps, and it is exactly the case a
//   count-only gate misses -- GATE 2 must test Vin invariance.
//
// Manufactured solution keeps c>0 (Pc(z) ~ 2), so the BASE classifies EVOL_
// HYPERBOLIC with w_+ incoming@LB and the square solve recovers it.  Two sampled
// states are then supplied to the probe:
//   * REVERSAL  {c = -2}  -> eigenvalue signs flip, Vin rotates ~53 deg
//                            -> GATE 2 FAIL  -> reference_robustness_ok()==false.
//   * BENIGN    {c = +2.4}-> same signs, eigvecs c-independent (0 deg rotation)
//                            -> all gates PASS -> reference_robustness_ok()==true.
// The oracle PASSES iff the base solve recovers AND each probe verdict matches its
// expectation -- i.e. the probe catches the genuine reversal and clears the benign
// perturbation.
//
// A second block exercises GATE 2's OTHER branch -- a COUNT change rather than a
// swap -- via SCALAR advection (Burgers)  d_t u + u d_z u = f,  A_z = u (1x1).  Here
// the single speed IS u: base u>0 gives split 1/0, a u<0 sample gives 0/1, so the
// in/out COUNT flips (the swap test above leaves the count at 1/1).  Same structure
// (base well-posed + recovering; only the sample reveals the danger), opposite
// detection path -- so both GATE 2 branches are covered in one driver.
//
// Build: the probe is opt-in in the header; this oracle self-enables it.  Auto-
// closure + base guard are default-on.  + -DCRONOS__WITH_SPQR / -DCRONOS__WITH_UMFPACK as
// usual.  Knobs reuse the TEST_HYP_* family.
// ===========================================================================

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <string>
#include <vector>

#include <armadillo>

// This oracle exists to exercise the reference-robustness probe: enable it.
#ifndef MC__OCFESLV_REFERENCE_ROBUSTNESS_PROBE
#define MC__OCFESLV_REFERENCE_ROBUSTNESS_PROBE
#endif

#ifndef OCFE_OCFESLV_HEADER
#define OCFE_OCFESLV_HEADER "ocfeslv.hpp"
#endif
#include OCFE_OCFESLV_HEADER

using namespace mc;

#ifndef TEST_HYP_NEL_T
#define TEST_HYP_NEL_T 3
#endif
#ifndef TEST_HYP_NEL_Z
#define TEST_HYP_NEL_Z 3
#endif
#ifndef TEST_HYP_NT
#define TEST_HYP_NT 6
#endif
#ifndef TEST_HYP_NZ
#define TEST_HYP_NZ 10
#endif
#ifndef TEST_HYP_TF
#define TEST_HYP_TF 0.5
#endif
#ifndef TEST_HYP_ZF
#define TEST_HYP_ZF 1.0
#endif
#ifndef TEST_HYP_MAXIT
#define TEST_HYP_MAXIT 40
#endif
#ifndef TEST_HYP_SOLVE_TOL
#define TEST_HYP_SOLVE_TOL 1e-9
#endif
#ifndef TEST_HYP_EXACT_TOL
#define TEST_HYP_EXACT_TOL 1e-6
#endif
// Sampled "reversed-flow" and "benign" wave-speed states fed to the probe.
#ifndef TEST_PDE17_C_REVERSAL
#define TEST_PDE17_C_REVERSAL (-2.0)
#endif
#ifndef TEST_PDE17_C_BENIGN
#define TEST_PDE17_C_BENIGN ( 2.4)
#endif

struct Par {
  double tf  = TEST_HYP_TF;
  double zf  = TEST_HYP_ZF;
  double g   = 4.0;            // A_z = c[[0,1],[g,0]]
  double sg  = 2.0;            // sqrt(g): left eigvecs (sg,+-1)
  // Pc(z) ~ 2 > 0 so the wave speed +-2c is sign-definite (forward) at the base
  // and the solve is well-posed; degree 3 <= N-1 so spectral d_z is exact.
  double cc[4] = { 2.0, 0.30, -0.20, 0.10 };  // Pc(z) > 0 on [0,zf]
  double uu[4] = { 0.4, 0.70,  0.30,-0.30 };  // Pu(z)
  double alp = 0.3;            // linear-in-time growth
};

static double polyval( double const c[4], double z ){ return c[0]+z*(c[1]+z*(c[2]+z*c[3])); }
static double polyder( double const c[4], double z ){ return c[1]+z*(2.0*c[2]+z*3.0*c[3]); }

// Manufactured solution: c = Pc(z)(1+alp t),  u = Pu(z)(1+alp t).
static double C_exact( double t, double z, Par const& p ){ return polyval(p.cc,z)*(1.0+p.alp*t); }
static double U_exact( double t, double z, Par const& p ){ return polyval(p.uu,z)*(1.0+p.alp*t); }

static char const* imp_name( OCFESLV::Options::ImpositionType t ){
  switch( t ){
    case OCFESLV::Options::IC_WEAK:   return "IC_WEAK";
    case OCFESLV::Options::IC_TRACE:  return "IC_TRACE";
    case OCFESLV::Options::IC_STRONG: return "IC_STRONG";
    default:                        return "IC_?";
  }
}
static double max_abs( std::vector<double> const& v ){
  double m=0.0; for( double x: v ) m=std::max(m,std::fabs(x)); return m;
}
static bool check_close( std::string const& label, double value, double tol ){
  bool const ok = ( value <= tol );
  std::cout << std::left << std::setw(40) << label
            << " value=" << std::scientific << std::setprecision(6) << value
            << " tol=" << tol << "  " << (ok?"PASS":"FAIL") << "\n";
  return ok;
}

struct CaseResult {
  std::string name;
  bool   solved=false, square=false, recovered=false;
  bool   probe_ok=true, probe_expected=true, probe_match=false;
  double final_res=0.0, eC=0.0, eU=0.0;
  bool   ok=false;
};

// reversal=true supplies a c<0 sample (expect probe FAIL); false a benign c>0
// sample (expect probe PASS).
static CaseResult run_case( OCFESLV::Options::ImpositionType imp, Par const& p, bool reversal )
{
  CaseResult R;
  R.name = std::string(reversal? "REVERSAL c<0" : "BENIGN c>0") + " [" + imp_name(imp) + "]";
  R.probe_expected = !reversal;            // reversal must trip the probe
  bool ok = true;

  std::cout << "\n====== quasilinear hyperbolic block, wave speed = state c ======\n";
  std::cout << "case: " << R.name
            << ",  sampled c = " << (reversal? TEST_PDE17_C_REVERSAL : TEST_PDE17_C_BENIGN)
            << "  (base c_ref ~ +" << p.cc[0] << ")\n";
  std::cout << "PDE: d_t c + c d_z u = f_c ; d_t u + g c d_z c = f_u,  g=" << p.g
            << " (eigenvalues +-2c, characteristics 2c+-u)\n";

  FFGraph DAG;
  FFVar t = DAG.add_var("t");
  FFVar z = DAG.add_var("z");
  FFVar c = DAG.add_var("c(t,z)");
  FFVar u = DAG.add_var("u(t,z)");
  FFPartial OpP;

  FFVar G   = 1.0 + p.alp*t;
  FFVar Pc  = p.cc[0] + p.cc[1]*z + p.cc[2]*z*z + p.cc[3]*z*z*z;
  FFVar Pu  = p.uu[0] + p.uu[1]*z + p.uu[2]*z*z + p.uu[3]*z*z*z;
  FFVar dPc = p.cc[1] + 2.0*p.cc[2]*z + 3.0*p.cc[3]*z*z;
  FFVar dPu = p.uu[1] + 2.0*p.uu[2]*z + 3.0*p.uu[3]*z*z;
  FFVar CE  = Pc*G;
  FFVar UE  = Pu*G;
  // Manufactured forcing (nonlinear: the c-coefficient makes the advection term
  // quadratic in the growth factor):
  //   f_c = d_t c +   c d_z u = alp Pc +   Pc Pu'(z) G^2
  //   f_u = d_t u + g c d_z c = alp Pu + g Pc Pc'(z) G^2
  FFVar FC = p.alp*Pc +       Pc*dPu*G*G;
  FFVar FU = p.alp*Pu + p.g*  Pc*dPc*G*G;

  FFVar PDEC = OpP(c,t) +        c*OpP(u,z) - FC;   // d_t c + c d_z u - f_c
  FFVar PDEU = OpP(u,t) + p.g*  c*OpP(c,z) - FU;   // d_t u + g c d_z c - f_u
  FFVar ICC  = c - CE;
  FFVar ICU  = u - UE;

  // Value inflow BCs on the manufactured characteristics (c>0 base: w_+ in@LB,
  // w_- in@UB).  Auto-closure transports the outgoing characteristics.
  FFVar grow  = G;
  FFVar wpVal = p.sg*c + u;                       // w_+ = 2c+u
  FFVar wmVal = p.sg*c - u;                       // w_- = 2c-u
  double const wpLB0 = p.sg*polyval(p.cc,0.0)  + polyval(p.uu,0.0);
  double const wmUB0 = p.sg*polyval(p.cc,p.zf) - polyval(p.uu,p.zf);
  FFVar BC_LB = wpVal - wpLB0*grow;               // prescribe w_+ @LB
  FFVar BC_UB = wmVal - wmUB0*grow;               // prescribe w_- @UB

  OCFESLV oc(&DAG);
  oc.add_domain( t, FFDom(0., p.tf, TEST_HYP_NEL_T, FFDom::LGR, TEST_HYP_NT) );
  oc.add_domain( z, FFDom(0., p.zf, TEST_HYP_NEL_Z, FFDom::LGL, TEST_HYP_NZ) );
  oc.add_state( c, {t,z} );
  oc.add_state( u, {t,z} );
  oc.update_ref( c, [&]( OCFESLV::t_Coord const& crd ){ return C_exact(crd.at(t),crd.at(z),p); } );
  oc.update_ref( u, [&]( OCFESLV::t_Coord const& crd ){ return U_exact(crd.at(t),crd.at(z),p); } );

  // Reference-robustness sample: only the wave-speed state c matters to A_z(c).
  double const c_sample = reversal ? double(TEST_PDE17_C_REVERSAL)
                                   : double(TEST_PDE17_C_BENIGN);
  oc.add_classification_reference_sample( { { c, c_sample } } );

  OCFESLV::EqnOptions int_opt( OCFESLV::EqnRole::INTERIOR, 0 );
  OCFESLV::EqnOptions ini_opt( OCFESLV::EqnRole::INITIAL,  0 );
  OCFESLV::EqnOptions bnd_opt( OCFESLV::EqnRole::BOUNDARY, 0 );

  int const T_INT = FFDom::ALL - FFDom::LB;
  int const Z_INT = FFDom::ALL - FFDom::LB - FFDom::UB;

  oc.add_equation( PDEC,  {t,z}, {T_INT, Z_INT},         int_opt );
  oc.add_equation( PDEU,  {t,z}, {T_INT, Z_INT},         int_opt );
  oc.add_equation( ICC,   {t,z}, {FFDom::LB, FFDom::ALL}, ini_opt );
  oc.add_equation( ICU,   {t,z}, {FFDom::LB, FFDom::ALL}, ini_opt );
  oc.add_equation( BC_LB, {t,z}, {T_INT, FFDom::LB},     bnd_opt );
  oc.add_equation( BC_UB, {t,z}, {T_INT, FFDom::UB},     bnd_opt );

  oc.set_evolution_domain( t );
  oc.options.REDUCE.ORDER    = OCFESLV::Options::RED_FULL;
  oc.options.CLASSIFY.MODE        = OCFESLV::Options::CLASS_AUTO;
  oc.options.INTERFACE.TYPE  = OCFESLV::Options::IC_AUTO;
  oc.options.INTERFACE.IMPOSITION = imp;
  oc.options.INTERFACE.SAT_SIGMA0      = 10.0;
  oc.options.DISPLAY_LEVEL   = 1;

  if( !oc.setup() ){
    std::cerr << "ERROR: OCFESLV::setup() failed for " << R.name
              << ": " << OCFESLV::setup_status_str( oc.setup_status() ) << "\n";
    R.ok = false; return R;
  }

  // Probe verdict (computed during classify, inside setup()).
  R.probe_ok    = oc.reference_robustness_ok();
  R.probe_match = ( R.probe_ok == R.probe_expected );

  size_t const nVar=oc.n_colloc_sta();
  size_t const nEqn=oc.n_colloc_eqn();
  R.square = (nVar==nEqn);
  std::cout << "nVar=" << nVar << " nEqn=" << nEqn
            << " square=" << (R.square?"yes":"no") << "\n";
  {
    auto const& cls = oc.pde_type();
    std::cout << "base PDE type: " << OCFESLV::pde_type_name(cls.type)
              << "  evolution_hyperbolic=" << (cls.evolution_hyperbolic?"yes":"no") << "\n";
  }

  // Nonlinear solve: initialise at the reference (manufactured) field so Newton/LM
  // starts in the c>0 basin (a zero start would make A_z(0)=0 singular).
  std::vector<double> var(nVar,0.0);
  {
    size_t off=0, sidx=0;
    for( auto const& st: oc.states_colloc() ){
      auto nodes = oc.node_colloc(st);
      for( size_t i=0;i<nodes.size();++i)
        var[off+i] = (sidx==0)? C_exact(nodes[i][0],nodes[i][1],p)
                              : U_exact(nodes[i][0],nodes[i][1],p);
      off += nodes.size(); ++sidx;
    }
  }

  oc.options.SOLVE.MAX_ITER = TEST_HYP_MAXIT;
  oc.options.SOLVE.RES_TOL  = TEST_HYP_SOLVE_TOL;
  OCFESLV::SolveReport const srep = oc.solve( var.data() );
  R.solved = srep.converged;
  if( !R.solved )
    std::cerr << "OCFESLV::solve did not converge: final|r|=" << srep.final_residual
              << " after " << srep.iterations << " it\n";

  std::vector<double> res(nEqn,0.0);
  oc.eval(res.data(),nullptr,var.data(),nullptr,nullptr);
  R.final_res = max_abs(res);

  {
    double eC=0.0, eU=0.0; size_t off=0, sidx=0;
    for( auto const& st: oc.states_colloc() ){
      auto nodes = oc.node_colloc(st);
      for( size_t i=0; i<nodes.size(); ++i ){
        double const ex = (sidx==0)? C_exact(nodes[i][0],nodes[i][1],p)
                                   : U_exact(nodes[i][0],nodes[i][1],p);
        double const e = std::fabs( var[off+i] - ex );
        if( sidx==0 ) eC = std::max(eC,e); else eU = std::max(eU,e);
      }
      off += nodes.size(); ++sidx;
    }
    R.eC=eC; R.eU=eU;
  }
  R.recovered = ( R.eC <= TEST_HYP_EXACT_TOL && R.eU <= TEST_HYP_EXACT_TOL );

  std::cout << "-- base solve --\n";
  ok &= check_close("square system (nVar==nEqn)", R.square?0.0:1.0, 0.0);
  ok &= R.solved;
  ok &= check_close("final max residual", R.final_res, TEST_HYP_SOLVE_TOL);
  ok &= check_close("max |c - c_exact|", R.eC, TEST_HYP_EXACT_TOL);
  ok &= check_close("max |u - u_exact|", R.eU, TEST_HYP_EXACT_TOL);

  std::cout << "-- reference-robustness probe --\n";
  std::cout << "  expected verdict: " << (R.probe_expected? "PASS (benign)":"FAIL (reversal)")
            << " ;  actual: " << (R.probe_ok? "PASS":"FAIL")
            << "  => " << (R.probe_match? "MATCH":"MISMATCH") << "\n";
  ok &= R.probe_match;

  R.ok = ok;
  std::cout << "case " << R.name << ": " << (ok?"PASS":"FAIL") << "\n";
  return R;
}

// ---------------------------------------------------------------------------
// COUNT-change variant: scalar advection (Burgers)  d_t u + u d_z u = f.
// A_z = u (1x1), so the single characteristic speed IS the reference value of u:
// u>0 is right-going (split 1/0); a NEGATIVE sample reverses it (split 0/1) -- the
// +/- COUNT changes.  This exercises GATE 2's count branch, which the symmetric
// +-lambda run_case above leaves untouched (there the count stays 1/1 and only the
// threshold-free Vin-swap test fires).  Same Tier-2 point: the base problem (u>0)
// is well-posed and recovers; only the c<0 sample reveals the reference-dependence.
// ---------------------------------------------------------------------------
static CaseResult run_scalar_case( OCFESLV::Options::ImpositionType imp, Par const& p,
                                   bool reversal )
{
  CaseResult R;
  R.name = std::string(reversal? "SCALAR REV u<0" : "SCALAR BEN u>0") + " [" + imp_name(imp) + "]";
  R.probe_expected = !reversal;
  bool ok = true;

  double const u_sample = reversal ? double(TEST_PDE17_C_REVERSAL)
                                   : double(TEST_PDE17_C_BENIGN);
  std::cout << "\n====== scalar advection (Burgers), wave speed = state u ======\n";
  std::cout << "case: " << R.name << ",  sampled u = " << u_sample
            << "  (base u_ref ~ +" << p.cc[0] << ")\n";
  std::cout << "PDE: d_t u + u d_z u = f,  A_z = u (eigenvalue = u; sign = flow direction)\n";

  FFGraph DAG;
  FFVar t = DAG.add_var("t");
  FFVar z = DAG.add_var("z");
  FFVar u = DAG.add_var("u(t,z)");
  FFPartial OpP;

  // Manufactured u = Pc(z)(1+alp t) (Pc>0 so base is forward); f = d_t u + u d_z u.
  FFVar G   = 1.0 + p.alp*t;
  FFVar Pc  = p.cc[0] + p.cc[1]*z + p.cc[2]*z*z + p.cc[3]*z*z*z;
  FFVar dPc = p.cc[1] + 2.0*p.cc[2]*z + 3.0*p.cc[3]*z*z;
  FFVar UE  = Pc*G;
  FFVar FU  = p.alp*Pc + Pc*dPc*G*G;            // alp Pc + Pc Pc'(z) G^2
  FFVar PDEU = OpP(u,t) + u*OpP(u,z) - FU;
  FFVar ICU  = u - UE;
  double const uLB0 = polyval(p.cc,0.0);
  FFVar BC_LB = u - uLB0*G;                      // inflow @LB (correct for u>0)

  OCFESLV oc(&DAG);
  oc.add_domain( t, FFDom(0., p.tf, TEST_HYP_NEL_T, FFDom::LGR, TEST_HYP_NT) );
  oc.add_domain( z, FFDom(0., p.zf, TEST_HYP_NEL_Z, FFDom::LGL, TEST_HYP_NZ) );
  oc.add_state( u, {t,z} );
  oc.update_ref( u, [&]( OCFESLV::t_Coord const& crd ){
                       return polyval(p.cc,crd.at(z))*(1.0+p.alp*crd.at(t)); } );
  oc.add_classification_reference_sample( { { u, u_sample } } );

  OCFESLV::EqnOptions int_opt( OCFESLV::EqnRole::INTERIOR, 0 );
  OCFESLV::EqnOptions ini_opt( OCFESLV::EqnRole::INITIAL,  0 );
  OCFESLV::EqnOptions bnd_opt( OCFESLV::EqnRole::BOUNDARY, 0 );
  int const T_INT = FFDom::ALL - FFDom::LB;
  int const Z_INT = FFDom::ALL - FFDom::LB - FFDom::UB;

  oc.add_equation( PDEU,  {t,z}, {T_INT, Z_INT},         int_opt );
  oc.add_equation( ICU,   {t,z}, {FFDom::LB, FFDom::ALL}, ini_opt );
  oc.add_equation( BC_LB, {t,z}, {T_INT, FFDom::LB},     bnd_opt );  // inflow @LB
  // outgoing @UB supplied by the framework auto-closure

  oc.set_evolution_domain( t );
  oc.options.REDUCE.ORDER    = OCFESLV::Options::RED_FULL;
  oc.options.CLASSIFY.MODE        = OCFESLV::Options::CLASS_AUTO;
  oc.options.INTERFACE.TYPE  = OCFESLV::Options::IC_AUTO;
  oc.options.INTERFACE.IMPOSITION = imp;
  oc.options.INTERFACE.SAT_SIGMA0      = 10.0;
  oc.options.DISPLAY_LEVEL   = 1;

  if( !oc.setup() ){
    std::cerr << "ERROR: OCFESLV::setup() failed for " << R.name
              << ": " << OCFESLV::setup_status_str( oc.setup_status() ) << "\n";
    R.ok = false; return R;
  }

  R.probe_ok    = oc.reference_robustness_ok();
  R.probe_match = ( R.probe_ok == R.probe_expected );

  size_t const nVar=oc.n_colloc_sta();
  size_t const nEqn=oc.n_colloc_eqn();
  R.square = (nVar==nEqn);
  std::cout << "nVar=" << nVar << " nEqn=" << nEqn << " square=" << (R.square?"yes":"no") << "\n";
  {
    auto const& cls = oc.pde_type();
    std::cout << "base PDE type: " << OCFESLV::pde_type_name(cls.type)
              << "  evolution_hyperbolic=" << (cls.evolution_hyperbolic?"yes":"no") << "\n";
  }

  std::vector<double> var(nVar,0.0);
  {
    size_t off=0;
    for( auto const& st: oc.states_colloc() ){
      auto nodes = oc.node_colloc(st);
      for( size_t i=0;i<nodes.size();++i)
        var[off+i] = polyval(p.cc,nodes[i][1])*(1.0+p.alp*nodes[i][0]);
      off += nodes.size();
    }
  }
  oc.options.SOLVE.MAX_ITER = TEST_HYP_MAXIT;
  oc.options.SOLVE.RES_TOL  = TEST_HYP_SOLVE_TOL;
  OCFESLV::SolveReport const srep = oc.solve( var.data() );
  R.solved = srep.converged;

  std::vector<double> res(nEqn,0.0);
  oc.eval(res.data(),nullptr,var.data(),nullptr,nullptr);
  R.final_res = max_abs(res);
  {
    double e=0.0; size_t off=0;
    for( auto const& st: oc.states_colloc() ){
      auto nodes = oc.node_colloc(st);
      for( size_t i=0;i<nodes.size();++i)
        e = std::max(e, std::fabs(var[off+i] - polyval(p.cc,nodes[i][1])*(1.0+p.alp*nodes[i][0])));
      off += nodes.size();
    }
    R.eC = e; R.eU = 0.0;
  }
  R.recovered = ( R.eC <= TEST_HYP_EXACT_TOL );

  std::cout << "-- base solve --\n";
  ok &= check_close("square system (nVar==nEqn)", R.square?0.0:1.0, 0.0);
  ok &= R.solved;
  ok &= check_close("final max residual", R.final_res, TEST_HYP_SOLVE_TOL);
  ok &= check_close("max |u - u_exact|", R.eC, TEST_HYP_EXACT_TOL);

  std::cout << "-- reference-robustness probe --\n";
  std::cout << "  expected verdict: " << (R.probe_expected? "PASS (benign)":"FAIL (reversal)")
            << " ;  actual: " << (R.probe_ok? "PASS":"FAIL")
            << "  => " << (R.probe_match? "MATCH":"MISMATCH") << "\n";
  ok &= R.probe_match;

  R.ok = ok;
  std::cout << "case " << R.name << ": " << (ok?"PASS":"FAIL") << "\n";
  return R;
}

int main()
{
  std::cout << "=============================================================\n";
  std::cout << " PDE17: reference-robustness probe oracle (characteristic sign flip)\n";
  std::cout << " SWAP path:  quasilinear A_z(c)=c[[0,1],[g,0]]; base c>0, sampled c<0\n";
  std::cout << " COUNT path: scalar Burgers A_z=u; base u>0, sampled u<0\n";
  std::cout << "=============================================================\n";

  Par p;

  // SWAP path (+-lambda system: count stays 1/1, incoming subspace swaps ~53 deg).
  CaseResult rev = run_case( OCFESLV::Options::IC_WEAK, p, /*reversal=*/true  );
  CaseResult ben = run_case( OCFESLV::Options::IC_WEAK, p, /*reversal=*/false );
  // COUNT path (scalar advection: split 1/0 -> 0/1).
  CaseResult srev = run_scalar_case( OCFESLV::Options::IC_WEAK, p, /*reversal=*/true  );
  CaseResult sben = run_scalar_case( OCFESLV::Options::IC_WEAK, p, /*reversal=*/false );

  std::cout << "\n====================================================================\n";
  std::cout << std::left << std::setw(22) << "case"
            << std::setw(10) << "solved" << std::setw(12) << "recovered"
            << std::setw(14) << "probe(exp)" << std::setw(14) << "probe(act)"
            << "result\n";
  auto row = []( CaseResult const& R ){
    std::cout << std::left << std::setw(22) << R.name
              << std::setw(10) << (R.solved?"yes":"no")
              << std::setw(12) << (R.recovered?"yes":"no")
              << std::setw(14) << (R.probe_expected?"PASS":"FAIL")
              << std::setw(14) << (R.probe_ok?"PASS":"FAIL")
              << (R.ok?"PASS":"FAIL") << "\n";
  };
  row(rev); row(ben); row(srev); row(sben);

  bool const all_ok = rev.ok && ben.ok && srev.ok && sben.ok;
  std::cout << "\nThe base problem is identical within each path and recovers the manufactured\n"
            << "solution; only the SAMPLED state differs.  GATE 2 fires on both reversals --\n"
            << "the +-lambda SWAP (count 1/1, incoming subspace rotates ~53 deg) and the scalar\n"
            << "COUNT change (1/0 -> 0/1) -- and clears both benign perturbations: a structural\n"
            << "danger invisible to the (well-posed, exactly-recovering) base solve.\n";
  std::cout << "\nreference-robustness sign-flip oracle (swap + count branches): "
            << (all_ok?"PASS":"FAIL") << "\n";
  return all_ok ? 0 : 1;
}
