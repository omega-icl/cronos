// ===========================================================================
// OCFE_PDE20_solve2.cpp
//
// STRUCTURAL-ANALYSIS VALIDATION CORPUS (Stage 1a).
//
// A range of DAE and IPDAE models of KNOWN (analytical) differentiation index,
// to validate the structural index/order analysis that Stage 1a will add to
// classify (symbolic, reference-free: FFDep incidence + symbolic dispersion +
// the setup-time matcher).  Each model documents its textbook index; the corpus
// reads what the framework currently reports and is the regression harness once
// Stage 1a exposes a structural-index accessor.
//
// The index ladder (pure DAE, time-only -- isolates the TEMPORAL index, no
// spatial order-raising to confound it):
//
//   M0  index 0  ODE          x' = -x + s(t)
//   M1  index 1  semi-explicit x' = y ;        0 = y - x - s(t)
//                              (constraint CONTAINS the algebraic var y:
//                               dG/dy = 1 != 0  -> one solve closes it)
//   M2  index 2  hidden var    x' = y ;        0 = x - a(t)
//                              (constraint on x only; y is HIDDEN -> differentiate
//                               0=x-a once -> y - a' = 0 exposes y)
//   M3  index 3  mechanical    x' = v ; v' = L ; 0 = x - a(t)
//                              (position constraint, no L; differentiate TWICE:
//                               v - a' = 0, then L - a'' = 0 exposes L)
//
// The IPDAE case (adds the SPATIAL order-raising on top of an index-1 link):
//
//   M4  index 1 (temporal), spatial order 2 (PARABOLIC) -- the PDE19 reduced
//       Darcy block:  d_t c + d_z(u c) = f ;  d_z c + R u = g
//       (u algebraic, pinned by momentum -> temporal index 1; eliminating u
//        raises the c-equation to 2nd order -> hidden parabolic)
//
// Non-DAE regression anchors (must keep classifying normally under Stage 1a):
//
//   M5  parabolic heat       d_t u = D d_zz u + f      -> PARABOLIC
//   M6  hyperbolic advection d_t u + a d_z u = f        -> HYPERBOLIC
//
// EXPECTED structural read (what Stage 1a must produce):
//   M0: index 0, regular evolution (ODE), square.
//   M1: index 1, reduces to an ODE in x, square.
//   M2: index 2 (1 differentiation), square after reduction.
//   M3: index 3 (2 differentiations), square after reduction.
//   M4: index 1, spatial order 2 -> PARABOLIC, square after reduction.
//
// BASELINE (current header, pre-Stage-1a) expectation: M0 classifies as a normal
// evolution/ODE and is square; M1-M4 hit the singular-evolution short-circuit ->
// DESCRIPTOR and are NON-square (no reduction/closure).  That gap is the Stage-1a
// target.  This corpus is DIAGNOSTIC-FIRST: it reports pde_type / square / counts
// per model and compares against the documented index; the index ASSERTIONS wire
// up once Stage 1a adds the structural-index accessor.
//
// Build: + -DCRONOS__WITH_SPQR.  Knobs: -DTEST_DAE_NT/_NZ, -DTEST_DAE_R.
// ===========================================================================

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include <armadillo>

#ifndef OCFE_OCFESLV_HEADER
#define OCFE_OCFESLV_HEADER "ocfeslv.hpp"
#endif
#include OCFE_OCFESLV_HEADER

using namespace mc;

#ifndef TEST_DAE_REDUCE
#define TEST_DAE_REDUCE RED_FULL
#endif

#ifndef TEST_DAE_NEL_T
#define TEST_DAE_NEL_T 3
#endif
#ifndef TEST_DAE_NT
#define TEST_DAE_NT 6
#endif
#ifndef TEST_DAE_NEL_Z
#define TEST_DAE_NEL_Z 3
#endif
#ifndef TEST_DAE_NZ
#define TEST_DAE_NZ 10
#endif
#ifndef TEST_DAE_R
#define TEST_DAE_R 1.5
#endif
#ifndef TEST_DAE_TF
#define TEST_DAE_TF 0.5
#endif
#ifndef TEST_DAE_ZF
#define TEST_DAE_ZF 0.5
#endif

static double max_abs( std::vector<double> const& v ){
  double m=0.0; for( double x: v ) m=std::max(m,std::fabs(x)); return m;
}

struct ModelResult {
  std::string name;
  int         known_index = -1;
  std::string expect_type;          // documented expected reduced type
  std::string got_type;             // current pde_type
  bool        got_evol_hyp = false;
  bool        setup_ok = false;
  bool        threw = false;        // setup() threw (e.g. unsupported model form)
  std::string status;
  size_t      nVar=0, nEqn=0, nTrace=0;
  bool        square=false;
};

// Common option/diagnostic tail: setup, read classification + counts, print.
static void finish( OCFESLV& oc, ModelResult& R )
{
  oc.options.REDUCE.ORDER    = OCFESLV::Options::TEST_DAE_REDUCE;
  oc.options.CLASSIFY.MODE        = OCFESLV::Options::CLASS_AUTO;
  oc.options.INTERFACE.TYPE  = OCFESLV::Options::IC_AUTO;
  oc.options.INTERFACE.IMPOSITION = OCFESLV::Options::IC_WEAK;
  // Algebraic-state boundary closure: load-bearing for M4 (and M7 once Phase 2
  // lands), inert on the index-0 / no-genuine-algebraic modes.  Was the global
  // -DMC__OCFESLV_AUTO_ALG_CLOSURE build flag; now a runtime Option set here so it
  // applies uniformly to every M0-M7 mode through this shared finish() helper.
  oc.options.DISPLAY_LEVEL   = 0;          // silence OCFESLV::setup [timing] + verbose stream
  oc.options.SOLVE.MAX_ITER  = 50;
  oc.options.SOLVE.RES_TOL   = 1e-9;
#ifdef CRONOS__WITH_SPQR
  oc.options.SOLVE.FACTORIZATION = OCFESLV::Options::SOLVE_SPQR;
#endif

  try {
    R.setup_ok = oc.setup();
  }
  catch( ... ) {
    R.threw    = true;
    R.setup_ok = false;
    R.got_type = "(setup threw)";
    R.status   = "setup() threw an exception";
    std::cout << "  setup=THREW  (model form not supported by setup; "
                 "see principal_symbol FAD)\n";
    return;
  }
  auto const pt = oc.pde_type();
  R.got_type     = OCFESLV::pde_type_name( pt.type );
  R.got_evol_hyp = pt.evolution_hyperbolic;
  R.status       = OCFESLV::setup_status_str( oc.setup_status() );

  std::cout << "  setup=" << (R.setup_ok?"OK":"FAIL")
            << "  status=" << R.status
            << "  pde_type=" << R.got_type
            << "  evolution_hyperbolic=" << (R.got_evol_hyp?"yes":"no") << "\n";

  if( R.setup_ok ){
    R.nVar=oc.n_colloc_sta(); R.nEqn=oc.n_colloc_eqn(); R.nTrace=oc.n_colloc_trace();
    R.square=(R.nVar==R.nEqn);
    std::cout << "  nVar=" << R.nVar << " nEqn=" << R.nEqn
              << " nTrace=" << R.nTrace << " square=" << (R.square?"yes":"no") << "\n";
  }
  std::cout << "  [known index " << R.known_index
            << ", expect " << R.expect_type << "]\n";
}

// ---- M0: index-0 ODE  x' = -x + s(t) --------------------------------------
static ModelResult run_M0()
{
  ModelResult R; R.name="M0 index-0 ODE  x'=-x+s"; R.known_index=0;
  R.expect_type="regular evolution (ODE), square";
  std::cout << "\n------ " << R.name << " ------\n";

  FFGraph DAG;
  FFVar t = DAG.add_var("t");
  FFVar x = DAG.add_var("x(t)");
  FFPartial OpP;
  // x_exact = 1 + 0.5 t - 0.3 t^2 ; s = x' + x
  FFVar Xe = 1.0 + 0.5*t - 0.3*t*t;
  FFVar Xt = 0.5 - 0.6*t;
  FFVar S  = Xt + Xe;
  FFVar ODE = OpP(x,t) + x - S;     // x' + x - s = 0
  FFVar ICx = x - Xe;

  OCFESLV oc(&DAG);
  oc.add_domain( t, FFDom(0.,TEST_DAE_TF,TEST_DAE_NEL_T,FFDom::LGR,TEST_DAE_NT) );
  oc.add_state( x, {t} );
  oc.update_ref( x, [&](OCFESLV::t_Coord const& c){ double tt=c.at(t); return 1.0+0.5*tt-0.3*tt*tt; } );
  OCFESLV::EqnOptions io(OCFESLV::EqnRole::INTERIOR,0), ii(OCFESLV::EqnRole::INITIAL,0);
  int const T_INT = FFDom::ALL - FFDom::LB;
  oc.add_equation( ODE, {t}, {T_INT},     io );
  oc.add_equation( ICx, {t}, {FFDom::LB}, ii );
  oc.set_evolution_domain( t );
  finish( oc, R );
  return R;
}

// ---- M1: index-1 semi-explicit  x'=y ; 0 = y - x - s(t) --------------------
static ModelResult run_M1()
{
  ModelResult R; R.name="M1 index-1 semi-explicit  x'=y; 0=y-x-s"; R.known_index=1;
  R.expect_type="index-1 -> ODE in x, square";
  std::cout << "\n------ " << R.name << " ------\n";

  FFGraph DAG;
  FFVar t = DAG.add_var("t");
  FFVar x = DAG.add_var("x(t)");
  FFVar y = DAG.add_var("y(t)");         // ALGEBRAIC
  FFPartial OpP;
  // x_exact poly; y = x + s ; x' = y  =>  s = x' - x ; y_exact = x'
  FFVar Xe = 1.0 + 0.5*t - 0.3*t*t;
  FFVar Xt = 0.5 - 0.6*t;
  FFVar S  = Xt - Xe;
  FFVar DIFF = OpP(x,t) - y;             // x' - y = 0
  FFVar ALG  = y - x - S;                // 0 = y - x - s  (CONTAINS y -> index 1)
  FFVar ICx  = x - Xe;                   // IC for x only (y algebraic)

  OCFESLV oc(&DAG);
  oc.add_domain( t, FFDom(0.,TEST_DAE_TF,TEST_DAE_NEL_T,FFDom::LGR,TEST_DAE_NT) );
  oc.add_state( x, {t} );
  oc.add_state( y, {t} );
  oc.update_ref( x, [&](OCFESLV::t_Coord const& c){ double tt=c.at(t); return 1.0+0.5*tt-0.3*tt*tt; } );
  oc.update_ref( y, [&](OCFESLV::t_Coord const& c){ double tt=c.at(t); return 0.5-0.6*tt; } );
  OCFESLV::EqnOptions io(OCFESLV::EqnRole::INTERIOR,0), ii(OCFESLV::EqnRole::INITIAL,0);
  int const T_INT = FFDom::ALL - FFDom::LB;
  oc.add_equation( DIFF, {t}, {T_INT},      io );
  oc.add_equation( ALG,  {t}, {FFDom::ALL}, io );   // algebraic: holds at all t
  oc.add_equation( ICx,  {t}, {FFDom::LB},  ii );
  oc.set_evolution_domain( t );
  finish( oc, R );
  return R;
}

// ---- M2: index-2  x'=y ; 0 = x - a(t)  (y hidden) -------------------------
static ModelResult run_M2()
{
  ModelResult R; R.name="M2 index-2 hidden  x'=y; 0=x-a"; R.known_index=2;
  R.expect_type="index-2 (1 differentiation), square after reduction";
  std::cout << "\n------ " << R.name << " ------\n";

  FFGraph DAG;
  FFVar t = DAG.add_var("t");
  FFVar x = DAG.add_var("x(t)");
  FFVar y = DAG.add_var("y(t)");         // ALGEBRAIC, HIDDEN (not in the constraint)
  FFPartial OpP;
  // x_exact = a(t) ; y_exact = a'(t)
  FFVar Ae = 1.0 + 0.5*t - 0.3*t*t;
  FFVar At = 0.5 - 0.6*t;
  FFVar DIFF = OpP(x,t) - y;             // x' - y = 0
  FFVar ALG  = x - Ae;                   // 0 = x - a  (NO y -> index 2)
  FFVar ICx  = x - Ae;

  OCFESLV oc(&DAG);
  oc.add_domain( t, FFDom(0.,TEST_DAE_TF,TEST_DAE_NEL_T,FFDom::LGR,TEST_DAE_NT) );
  oc.add_state( x, {t} );
  oc.add_state( y, {t} );
  oc.update_ref( x, [&](OCFESLV::t_Coord const& c){ double tt=c.at(t); return 1.0+0.5*tt-0.3*tt*tt; } );
  oc.update_ref( y, [&](OCFESLV::t_Coord const& c){ double tt=c.at(t); return 0.5-0.6*tt; } );
  OCFESLV::EqnOptions io(OCFESLV::EqnRole::INTERIOR,0), ii(OCFESLV::EqnRole::INITIAL,0);
  int const T_INT = FFDom::ALL - FFDom::LB;
  oc.add_equation( DIFF, {t}, {T_INT},      io );
  oc.add_equation( ALG,  {t}, {FFDom::ALL}, io );
  oc.add_equation( ICx,  {t}, {FFDom::LB},  ii );
  oc.set_evolution_domain( t );
  finish( oc, R );
  return R;
}

// ---- M3: index-3 mechanical  x'=v ; v'=L ; 0 = x - a(t) -------------------
static ModelResult run_M3()
{
  ModelResult R; R.name="M3 index-3 mechanical  x'=v; v'=L; 0=x-a"; R.known_index=3;
  R.expect_type="index-3 (2 differentiations), square after reduction";
  std::cout << "\n------ " << R.name << " ------\n";

  FFGraph DAG;
  FFVar t = DAG.add_var("t");
  FFVar x = DAG.add_var("x(t)");
  FFVar v = DAG.add_var("v(t)");
  FFVar L = DAG.add_var("L(t)");         // ALGEBRAIC multiplier, doubly hidden
  FFPartial OpP;
  // x_exact = a ; v_exact = a' ; L_exact = a''
  FFVar Ae = 1.0 + 0.5*t - 0.3*t*t;      // a
  FFVar At = 0.5 - 0.6*t;                // a'
  // a'' = -0.6
  FFVar DIFFx = OpP(x,t) - v;            // x' - v = 0
  FFVar DIFFv = OpP(v,t) - L;            // v' - L = 0
  FFVar ALG   = x - Ae;                  // 0 = x - a  (NO v, NO L -> index 3)
  FFVar ICx   = x - Ae;
  FFVar ICv   = v - At;

  OCFESLV oc(&DAG);
  oc.add_domain( t, FFDom(0.,TEST_DAE_TF,TEST_DAE_NEL_T,FFDom::LGR,TEST_DAE_NT) );
  oc.add_state( x, {t} );
  oc.add_state( v, {t} );
  oc.add_state( L, {t} );
  oc.update_ref( x, [&](OCFESLV::t_Coord const& c){ double tt=c.at(t); return 1.0+0.5*tt-0.3*tt*tt; } );
  oc.update_ref( v, [&](OCFESLV::t_Coord const& c){ double tt=c.at(t); return 0.5-0.6*tt; } );
  oc.update_ref( L, [&](OCFESLV::t_Coord const& c){ (void)c; return -0.6; } );
  OCFESLV::EqnOptions io(OCFESLV::EqnRole::INTERIOR,0), ii(OCFESLV::EqnRole::INITIAL,0);
  int const T_INT = FFDom::ALL - FFDom::LB;
  oc.add_equation( DIFFx, {t}, {T_INT},      io );
  oc.add_equation( DIFFv, {t}, {T_INT},      io );
  oc.add_equation( ALG,   {t}, {FFDom::ALL}, io );
  oc.add_equation( ICx,   {t}, {FFDom::LB},  ii );
  oc.add_equation( ICv,   {t}, {FFDom::LB},  ii );
  oc.set_evolution_domain( t );
  finish( oc, R );
  return R;
}

// ---- M4: index-1 IPDAE (PDE19 reduced Darcy)  d_t c + d_z(uc)=f ; d_z c+Ru=g
static ModelResult run_M4()
{
  ModelResult R; R.name="M4 index-1 IPDAE (Darcy)  d_t c+d_z(uc)=f; d_z c+Ru=g"; R.known_index=1;
  R.expect_type="index-1, spatial order 2 -> PARABOLIC, square after reduction";
  std::cout << "\n------ " << R.name << " ------\n";

  double const R_=TEST_DAE_R, alp=0.3;
  FFGraph DAG;
  FFVar t = DAG.add_var("t");
  FFVar z = DAG.add_var("z");
  FFVar c = DAG.add_var("c(t,z)");
  FFVar u = DAG.add_var("u(t,z)");       // ALGEBRAIC
  FFPartial OpP;
  FFVar Ce  = (1.0 + 0.4*z - 0.3*z*z)*(1.0+alp*t);
  FFVar Cze = (0.4 - 0.6*z)*(1.0+alp*t);
  FFVar Ue  = (0.8 + 0.3*z)*(1.0+alp*t);
  FFVar tfac= 1.0+alp*t;
  // f_c = d_t c + d_z(u c) ; g = d_z c + R u  (manufactured)
  FFVar FC = alp*(1.0+0.4*z-0.3*z*z)
           + tfac*tfac*( 0.3*(1.0+0.4*z-0.3*z*z) + (0.8+0.3*z)*(0.4-0.6*z) );
  FFVar G  = Cze + R_*Ue;
  FFVar PDEC = OpP(c,t) + u*OpP(c,z) + c*OpP(u,z) - FC;
  FFVar PDEU = OpP(c,z) + R_*u - G;
  FFVar ICC  = c - Ce;
  FFVar BCC  = c - Ce;

  OCFESLV oc(&DAG);
  oc.add_domain( t, FFDom(0.,TEST_DAE_TF,TEST_DAE_NEL_T,FFDom::LGR,TEST_DAE_NT) );
  oc.add_domain( z, FFDom(0.,TEST_DAE_ZF,TEST_DAE_NEL_Z,FFDom::LGL,TEST_DAE_NZ) );
  oc.add_state( c, {t,z} );
  oc.add_state( u, {t,z} );
  oc.update_ref( c, [&](OCFESLV::t_Coord const& cc){ double tt=cc.at(t),zz=cc.at(z); return (1.0+0.4*zz-0.3*zz*zz)*(1.0+alp*tt); } );
  oc.update_ref( u, [&](OCFESLV::t_Coord const& cc){ double tt=cc.at(t),zz=cc.at(z); return (0.8+0.3*zz)*(1.0+alp*tt); } );
  OCFESLV::EqnOptions io(OCFESLV::EqnRole::INTERIOR,0), ii(OCFESLV::EqnRole::INITIAL,0), bo(OCFESLV::EqnRole::BOUNDARY,0);
  int const T_INT = FFDom::ALL - FFDom::LB;
  int const Z_INT = FFDom::ALL - FFDom::LB - FFDom::UB;
  oc.add_equation( PDEC, {t,z}, {T_INT,      Z_INT},     io );
  oc.add_equation( PDEU, {t,z}, {FFDom::ALL, Z_INT},     io );
  oc.add_equation( ICC,  {t,z}, {FFDom::LB,  FFDom::ALL}, ii );
  oc.add_equation( BCC,  {t,z}, {T_INT,      FFDom::LB}, bo );
  oc.add_equation( BCC,  {t,z}, {T_INT,      FFDom::UB}, bo );
  oc.set_evolution_domain( t );
  finish( oc, R );
  return R;
}

// ---- M5: parabolic heat (non-DAE anchor)  d_t u = D d_zz u + f -------------
static ModelResult run_M5()
{
  ModelResult R; R.name="M5 parabolic heat  d_t u = D d_zz u + f"; R.known_index=0;
  R.expect_type="PARABOLIC, square (non-DAE regression anchor)";
  std::cout << "\n------ " << R.name << " ------\n";

  double const D=0.05, alp=0.3;
  FFGraph DAG;
  FFVar t = DAG.add_var("t");
  FFVar z = DAG.add_var("z");
  FFVar u = DAG.add_var("u(t,z)");
  FFPartial OpP;
  // u_exact = (1 + 0.4 z - 0.3 z^2)(1+alp t); u_zz = -0.6 (1+alp t); u_t = alp(1+0.4z-0.3z^2)
  FFVar Ue  = (1.0 + 0.4*z - 0.3*z*z)*(1.0+alp*t);
  FFVar F   = alp*(1.0+0.4*z-0.3*z*z) - D*(-0.6)*(1.0+alp*t);
  FFVar PDE = OpP(u,t) - D*OpP(u,{z,2}) - F;      // d_t u - D d_zz u - f
  FFVar ICU = u - Ue;
  FFVar BCU = u - Ue;                             // Dirichlet at both ends (2 BCs)

  OCFESLV oc(&DAG);
  oc.add_domain( t, FFDom(0.,TEST_DAE_TF,TEST_DAE_NEL_T,FFDom::LGR,TEST_DAE_NT) );
  oc.add_domain( z, FFDom(0.,TEST_DAE_ZF,TEST_DAE_NEL_Z,FFDom::LGL,TEST_DAE_NZ) );
  oc.add_state( u, {t,z} );
  oc.update_ref( u, [&](OCFESLV::t_Coord const& c){ double tt=c.at(t),zz=c.at(z); return (1.0+0.4*zz-0.3*zz*zz)*(1.0+alp*tt); } );
  OCFESLV::EqnOptions io(OCFESLV::EqnRole::INTERIOR,0), ii(OCFESLV::EqnRole::INITIAL,0), bo(OCFESLV::EqnRole::BOUNDARY,0);
  int const T_INT = FFDom::ALL - FFDom::LB;
  int const Z_INT = FFDom::ALL - FFDom::LB - FFDom::UB;
  oc.add_equation( PDE, {t,z}, {T_INT,     Z_INT},     io );
  oc.add_equation( ICU, {t,z}, {FFDom::LB, FFDom::ALL}, ii );
  oc.add_equation( BCU, {t,z}, {T_INT,     FFDom::LB}, bo );
  oc.add_equation( BCU, {t,z}, {T_INT,     FFDom::UB}, bo );
  oc.set_evolution_domain( t );
  finish( oc, R );
  return R;
}

// ---- M6: hyperbolic advection (non-DAE anchor)  d_t u + a d_z u = f --------
static ModelResult run_M6()
{
  ModelResult R; R.name="M6 hyperbolic advection  d_t u + a d_z u = f"; R.known_index=0;
  R.expect_type="HYPERBOLIC, square (non-DAE regression anchor)";
  std::cout << "\n------ " << R.name << " ------\n";

  double const a=1.0, alp=0.3;
  FFGraph DAG;
  FFVar t = DAG.add_var("t");
  FFVar z = DAG.add_var("z");
  FFVar u = DAG.add_var("u(t,z)");
  FFPartial OpP;
  FFVar Ue  = (1.0 + 0.4*z - 0.3*z*z)*(1.0+alp*t);
  FFVar Uz  = (0.4 - 0.6*z)*(1.0+alp*t);
  FFVar F   = alp*(1.0+0.4*z-0.3*z*z) + a*(0.4-0.6*z)*(1.0+alp*t);
  (void)Uz;
  FFVar PDE = OpP(u,t) + a*OpP(u,z) - F;          // d_t u + a d_z u - f
  FFVar ICU = u - Ue;
  FFVar BCU = u - Ue;                             // 1 inflow BC at z=LB (a>0)

  OCFESLV oc(&DAG);
  oc.add_domain( t, FFDom(0.,TEST_DAE_TF,TEST_DAE_NEL_T,FFDom::LGR,TEST_DAE_NT) );
  oc.add_domain( z, FFDom(0.,TEST_DAE_ZF,TEST_DAE_NEL_Z,FFDom::LGL,TEST_DAE_NZ) );
  oc.add_state( u, {t,z} );
  oc.update_ref( u, [&](OCFESLV::t_Coord const& c){ double tt=c.at(t),zz=c.at(z); return (1.0+0.4*zz-0.3*zz*zz)*(1.0+alp*tt); } );
  OCFESLV::EqnOptions io(OCFESLV::EqnRole::INTERIOR,0), ii(OCFESLV::EqnRole::INITIAL,0), bo(OCFESLV::EqnRole::BOUNDARY,0);
  int const T_INT = FFDom::ALL - FFDom::LB;
  int const Z_INT = FFDom::ALL - FFDom::LB - FFDom::UB;
  oc.add_equation( PDE, {t,z}, {T_INT,     Z_INT},     io );
  oc.add_equation( ICU, {t,z}, {FFDom::LB, FFDom::ALL}, ii );
  oc.add_equation( BCU, {t,z}, {T_INT,     FFDom::LB}, bo );   // inflow only
  oc.set_evolution_domain( t );
  finish( oc, R );
  return R;
}

// ---- M7: conservative-flux constraint  d_z(w c)=g  (w only under d_z) -------
//   VALIDATION of conservative-form support.  A modeller may legitimately write
//   a flux in conservative form d_z(w*c).  This once aborted: principal_symbol's
//   derivative-proxy pass only recognises OpP(state,dir), so the product-operand
//   partial reached FAD raw and FFPartial::deriv aborted.  _reduce_order now
//   chain-rule normalises any partial of a compound operand,
//       d_z(f) -> (df/dz)|_explicit + sum_s (df/ds) d_z(s) ,
//   so d_z(w*c) -> c*d_z(w) + w*d_z(c) BEFORE classification.  Consequences,
//   all confirmed: setup succeeds; the matcher sees w BARE in w*d_z(c) and pins
//   it -> INDEX 1; character is parabolic and pde_type DESCRIPTOR -- structurally
//   identical to the expanded twin M4 (index-1, parabolic, non-square), i.e. the
//   index is form-invariant.  The try/catch guard in finish() stays as defensive
//   infrastructure for any future throwing model but is no longer exercised here.
static ModelResult run_M7()
{
  ModelResult R; R.name="M7 conservative flux  d_z(w c)=g (w only under d_z)"; R.known_index=1;
  R.expect_type="conservative d_z(wc) chain-rule normalised in reduce_order -> w bare -> index-1 [parabolic], matches twin M4";
  std::cout << "\n------ " << R.name << " ------\n";

  double const alp=0.3;
  FFGraph DAG;
  FFVar t = DAG.add_var("t");
  FFVar z = DAG.add_var("z");
  FFVar c = DAG.add_var("c(t,z)");
  FFVar w = DAG.add_var("w(t,z)");                 // ALGEBRAIC, appears ONLY inside d_z(w*c)
  FFPartial OpP;
  FFVar Fc = alp*(1.0+0.4*z-0.3*z*z) + (0.4-0.6*z)*(1.0+alp*t);    // d_t c + d_z c
  FFVar G  = (0.62 - 0.24*z - 0.27*z*z)*(1.0+alp*t)*(1.0+alp*t);   // d_z(w*c) at exact
  FFVar Ce = (1.0 + 0.4*z - 0.3*z*z)*(1.0+alp*t);
  FFVar PDEC   = OpP(c,t) + OpP(c,z) - Fc;          // dynamic advection in c
  FFVar CONSTR = OpP(w*c, z) - G;                   // conservative flux: w ONLY inside d_z(product)
  FFVar ICC = c - Ce;
  FFVar BCC = c - Ce;

  OCFESLV oc(&DAG);
  oc.add_domain( t, FFDom(0.,TEST_DAE_TF,TEST_DAE_NEL_T,FFDom::LGR,TEST_DAE_NT) );
  oc.add_domain( z, FFDom(0.,TEST_DAE_ZF,TEST_DAE_NEL_Z,FFDom::LGL,TEST_DAE_NZ) );
  oc.add_state( c, {t,z} );
  oc.add_state( w, {t,z} );
  oc.update_ref( c, [&](OCFESLV::t_Coord const& cc){ double tt=cc.at(t),zz=cc.at(z); return (1.0+0.4*zz-0.3*zz*zz)*(1.0+alp*tt); } );
  oc.update_ref( w, [&](OCFESLV::t_Coord const& cc){ double tt=cc.at(t),zz=cc.at(z); return (0.8+0.3*zz)*(1.0+alp*tt); } );
  OCFESLV::EqnOptions io(OCFESLV::EqnRole::INTERIOR,0), ii(OCFESLV::EqnRole::INITIAL,0), bo(OCFESLV::EqnRole::BOUNDARY,0);
  int const T_INT = FFDom::ALL - FFDom::LB;
  int const Z_INT = FFDom::ALL - FFDom::LB - FFDom::UB;
  oc.add_equation( PDEC,   {t,z}, {T_INT,      Z_INT},      io );
  oc.add_equation( CONSTR, {t,z}, {FFDom::ALL, Z_INT},      io );
  oc.add_equation( ICC,    {t,z}, {FFDom::LB,  FFDom::ALL},  ii );
  oc.add_equation( BCC,    {t,z}, {T_INT,      FFDom::LB},  bo );
  oc.set_evolution_domain( t );
  finish( oc, R );
  return R;
}

int main()
{
  std::cout << "PDE20 structural-analysis corpus  (index ladder + IPDAE)\n";
  std::cout << "build: SPQR="
#ifdef CRONOS__WITH_SPQR
            << "ON"
#else
            << "off"
#endif
            << "\n";

  std::vector<ModelResult> all;
  all.push_back( run_M0() );
  all.push_back( run_M1() );
  all.push_back( run_M2() );
  all.push_back( run_M3() );
  all.push_back( run_M4() );
  all.push_back( run_M5() );
  all.push_back( run_M6() );
  all.push_back( run_M7() );

  std::cout << "\n==================== SUMMARY ====================\n";
  std::cout << std::left << std::setw(46) << "model"
            << std::setw(7) << "index"
            << std::setw(14) << "pde_type"
            << std::setw(8) << "square" << "  expected\n";
  for( auto const& R : all ){
    std::cout << std::left << std::setw(46) << R.name.substr(0,45)
              << std::setw(7) << R.known_index
              << std::setw(14) << (R.threw ? std::string("(setup THREW)")
                                  : R.setup_ok? R.got_type : std::string("(setup FAIL)"))
              << std::setw(8) << (R.setup_ok ? (R.square?"yes":"NO") : "-")
              << "  " << R.expect_type << "\n";
  }
  std::cout << "\nNote: index ASSERTIONS pending the Stage-1a structural-index accessor.\n"
            << "Baseline gap = any algebraic model (M1-M4) reported DESCRIPTOR / non-square.\n";
  return 0;
}
