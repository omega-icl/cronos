// ============================================================================
//  OCFE_DAE1.cpp  --  The classical INDEX-3 pendulum posed DIRECTLY (position
//                     constraint), to exercise the AUTOMATIC index reduction
//                     (Pantelides + dummy-derivative) and the post-reduction
//                     structural DoF audit.
//
//  OCFE_DAE0/DAE0b pose the pendulum in the hand-reduced INDEX-1 form.  Here we
//  instead hand the framework the raw INDEX-3 system and let setup() detect the
//  index, build and execute the reduction plan, and audit the reduced DOF:
//
//        x' = u
//        y' = v
//        u' = -lam x
//        v' = -lam y - g
//        0  = x^2 + y^2 - L^2          <-- POSITION constraint; lam is HIDDEN
//                                          (absent from the constraint) => index 3
//
//  Structural picture (cf. the minimal M3 gate x'=v, v'=L, 0=x-a): lam is
//  unpinned by the algebraic constraint and is exposed only after TWO
//  differentiations of the constraint along the flow,
//        g0 = x^2+y^2-L^2                        (position, level 0)
//        g1 = x u + y v                          (velocity, level 1)
//        g2 = u^2+v^2 - lam(x^2+y^2) - g y       (acceleration, level 2; pins lam)
//  so the differential index is 3, witness {lam}.  The reducer replaces the
//  interior constraint with g2 (pinning lam) and keeps g0 at the initial node;
//  g1 is the hidden velocity consistency.  IC convention (per the M3 gate):
//  ONE INITIAL-role IC per DYNAMIC state (x,y,u,v); the algebraic/hidden lam
//  gets none; the reducer resolves the LB redundancy.
//
//  WHAT THIS ASSERTS (monolithic, per imposition)
//    S1  pde_type().differential_index == 3
//    S2  reduction happened  (reduction_plan() non-empty, max_index==3, resolved)
//    S3  reduced_dof_audit().ok()  (ran, square, full structural rank)
//    C0  the de-indexed system solves (from an oracle+perturbation seed)
//    P*  constraint drift / energy / RK4-oracle / lam checks (as in DAE0b)
//    XM  IC_WEAK ~ IC_STRONG reduced solution (INFORMATIONAL delta; each monolithic
//        reduction is gated vs the RK4 oracle in gate(), above)
//
//  Marching is run as an INFORMATIONAL diagnostic only (index reduction is a
//  monolithic setup-time analysis; marching a reduced high-index DAE is not
//  asserted here) -- its lines print but do not gate the suite.
//
//  Build (SPQR strongly recommended for the reduced saddle system):
//    g++ -std=c++17 <suite flags> -DOCFE_OCFESLV_HEADER='"ocfeslv.hpp"' \
//        OCFE_DAE1.cpp -o OCFE_DAE1 <libs>
// ============================================================================

#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#ifndef OCFE_OCFESLV_HEADER
#define OCFE_OCFESLV_HEADER "ocfeslv.hpp"
#endif
#include OCFE_OCFESLV_HEADER

using namespace mc;

static double const L_rod  = 1.0;
static double const g_grav = 9.81;
static double const th0    = M_PI / 3.0;
static double const om0    = 0.0;
static double const T_end  = 0.5;
static size_t const N_EL   = 8;
static size_t const N_ND   = 5;

// ---- independent oracle: RK4 of  theta'' = -(g/L) sin theta ----
struct Oracle {
  size_t              N;
  double              h;
  std::vector<double> th, om;
  Oracle( size_t Ngrid ) : N( Ngrid ), h( T_end / double( Ngrid ) ),
                           th( Ngrid + 1 ), om( Ngrid + 1 )
  {
    double t = th0, w = om0;
    th[0] = t; om[0] = w;
    auto f_t = []( double, double w_ ){ return w_; };
    auto f_w = []( double th_, double ){ return -( g_grav / L_rod ) * std::sin( th_ ); };
    for( size_t i = 0; i < N; ++i ){
      double k1t = f_t( t, w ),                          k1w = f_w( t, w );
      double k2t = f_t( t + 0.5*h*k1t, w + 0.5*h*k1w ),  k2w = f_w( t + 0.5*h*k1t, w + 0.5*h*k1w );
      double k3t = f_t( t + 0.5*h*k2t, w + 0.5*h*k2w ),  k3w = f_w( t + 0.5*h*k2t, w + 0.5*h*k2w );
      double k4t = f_t( t + h*k3t, w + h*k3w ),          k4w = f_w( t + h*k3t, w + h*k3w );
      t += ( h / 6.0 ) * ( k1t + 2*k2t + 2*k3t + k4t );
      w += ( h / 6.0 ) * ( k1w + 2*k2w + 2*k3w + k4w );
      th[i+1] = t; om[i+1] = w;
    }
  }
  void at_interp( double t, double& theta, double& omega ) const
  {
    if( t <= 0.0 ){ theta = th.front(); omega = om.front(); return; }
    if( t >= T_end ){ theta = th.back();  omega = om.back();  return; }
    double s = t / h; size_t i = (size_t)s; double f = s - double( i );
    if( i >= N ){ theta = th.back(); omega = om.back(); return; }
    theta = ( 1.0 - f ) * th[i] + f * th[i+1];
    omega = ( 1.0 - f ) * om[i] + f * om[i+1];
  }
  void at_node( double t, double& theta, double& omega ) const
  {
    long idx = (long)std::llround( t / h );
    if( idx < 0 ) idx = 0;
    if( idx > (long)N ) idx = (long)N;
    theta = th[idx]; omega = om[idx];
  }
  static void cart( double theta, double omega,
                    double& x, double& y, double& u, double& v, double& lam )
  {
    x =  L_rod * std::sin( theta );
    y = -L_rod * std::cos( theta );
    u =  L_rod * omega * std::cos( theta );
    v =  L_rod * omega * std::sin( theta );
    lam = omega*omega + g_grav * std::cos( theta ) / L_rod;
  }
};

static Oracle const ORC( 200000 );
static std::vector<double> const T_SAMPLE = { 0.0, 0.25*T_end, 0.5*T_end, 0.75*T_end, T_end };

static double energy0()
{
  double x,y,u,v,lam; Oracle::cart( th0, om0, x,y,u,v,lam );
  return 0.5*( u*u + v*v ) + g_grav * y;
}

static int g_pass = 0, g_fail = 0;
static void check_true( char const* name, bool ok )
{
  ( ok ? g_pass : g_fail )++;
  std::cout << "  " << std::left << std::setw(42) << name << "  " << ( ok ? "PASS" : "FAIL" ) << "\n";
}
static void check_le( char const* name, double got, double bound )
{
  bool const ok = ( got <= bound ) && std::isfinite( got );
  ( ok ? g_pass : g_fail )++;
  std::cout << "  " << std::left << std::setw(42) << name
            << std::right << std::scientific << std::setprecision(4)
            << " val=" << std::setw(12) << got << " <= " << std::setw(10) << bound
            << "  " << ( ok ? "PASS" : "FAIL" ) << "\n";
}
static void check_close( char const* name, double got, double want, double tol )
{
  bool const ok = std::fabs( got - want ) <= tol;
  ( ok ? g_pass : g_fail )++;
  std::cout << "  " << std::left << std::setw(42) << name
            << std::right << std::scientific << std::setprecision(4)
            << " |d|=" << std::setw(12) << std::fabs( got - want )
            << " tol=" << std::setw(10) << tol
            << "  " << ( ok ? "PASS" : "FAIL" ) << "\n";
}

struct ModeOut {
  bool   setup_ok  = false, converged = false;
  int    index     = -99;
  bool   reduced   = false, resolved = false;
  int    max_index = -1;
  bool   audit_ran = false, audit_sq = false, audit_ok = false;
  size_t rows = 0, cols = 0, rank = 0, defc = 0, nAssign = 0;
  double maxg0 = 9e9, maxg1 = 9e9, maxE = 9e9, maxTraj = 9e9, maxLam = 9e9;
  std::vector<double> F;
};

static ModeOut run_mode( OCFESLV::Options::ImpositionType imp, char const* strimp, bool marching,
                         OCFESLV::Options::ReductionType red = OCFESLV::Options::RED_FULL )
{
  ModeOut R;
  char const* mstr = marching ? "march" : "mono";
  std::cout << "\n---------------- IC_" << strimp << " [" << mstr << "] (index-3 direct) ----------------\n";

  FFGraph DAG;
  FFVar t   = DAG.add_var( "t"     );
  FFVar x   = DAG.add_var( "x(t)"  );
  FFVar y   = DAG.add_var( "y(t)"  );
  FFVar u   = DAG.add_var( "u(t)"  );
  FFVar v   = DAG.add_var( "v(t)"  );
  FFVar lam = DAG.add_var( "lam(t)");
  FFPartial OpP;

  OCFESLV oc( &DAG );
  oc.add_domain( t, FFDom( 0., T_end, N_EL, FFDom::LGR, N_ND ) );
  oc.add_state( x,   {t} );
  oc.add_state( y,   {t} );
  oc.add_state( u,   {t} );
  oc.add_state( v,   {t} );
  oc.add_state( lam, {t} );
  oc.set_evolution_domain( t );

  // Seed from the oracle + a small perturbation, so the de-indexed solve is
  // non-vacuous (must actively re-pin lam and correct back to the trajectory).
  auto seed = [&]( int which ){
    return [&,which]( OCFESLV::t_Coord const& c ) -> double {
      double th,om,X,Y,U,V,Lm; ORC.at_interp( c.at(t), th, om ); Oracle::cart(th,om,X,Y,U,V,Lm);
      double base = ( which==0?X : which==1?Y : which==2?U : which==3?V : Lm );
      double pert = 1e-3 * std::sin( 3.1*c.at(t) + 0.5*which );
      return base + pert;
    };
  };
  oc.update_ref( x,   seed(0) );
  oc.update_ref( y,   seed(1) );
  oc.update_ref( u,   seed(2) );
  oc.update_ref( v,   seed(3) );
  oc.update_ref( lam, seed(4) );

  double x0,y0,u0,v0,l0; Oracle::cart( th0, om0, x0,y0,u0,v0,l0 );

  FFVar EVOL_X = OpP( x, t ) - u;
  FFVar EVOL_Y = OpP( y, t ) - v;
  FFVar EVOL_U = OpP( u, t ) + lam * x;
  FFVar EVOL_V = OpP( v, t ) + lam * y + g_grav;
  FFVar ALG    = x*x + y*y - L_rod*L_rod;      // INDEX-3 position constraint (lam hidden)
  FFVar IC_X = x - x0;
  FFVar IC_Y = y - y0;
  FFVar IC_U = u - u0;
  FFVar IC_V = v - v0;

  int const T_NO_LB = FFDom::ALL - FFDom::LB;
  OCFESLV::EqnOptions const io( OCFESLV::EqnRole::INTERIOR, 0 );
  OCFESLV::EqnOptions const ii( OCFESLV::EqnRole::INITIAL,  0 );
  oc.add_equation( EVOL_X, {t}, {T_NO_LB},    io );
  oc.add_equation( EVOL_Y, {t}, {T_NO_LB},    io );
  oc.add_equation( EVOL_U, {t}, {T_NO_LB},    io );
  oc.add_equation( EVOL_V, {t}, {T_NO_LB},    io );
  oc.add_equation( ALG,    {t}, {FFDom::ALL}, io );
  // 2026-09-23: the FREE initial data only.  This system is index 3, so of its four differential states just TWO
  // initial values may be chosen; y and v follow from the hidden constraints (position x^2+y^2=L^2 and velocity
  // xu+yv=0), which REDUCE.HIDDEN_IC materialises at t=0.  The reference values pick the branch (y < 0 here).
  // Declaring all four was accepted and square, but two of the rows had to be made consistent by hand.
  oc.add_equation( IC_X,   {t}, {FFDom::LB},  ii );   // free
  oc.add_equation( IC_U,   {t}, {FFDom::LB},  ii );   // free
  (void)IC_Y; (void)IC_V;                             // implied by the hidden constraints, not declared

  FFVar OUT_G0 = x*x + y*y - L_rod*L_rod;
  FFVar OUT_G1 = x*u + y*v;
  FFVar OUT_E  = 0.5*( u*u + v*v ) + g_grav * y;
  for( double ts : T_SAMPLE ){
    oc.add_output( x,      {t}, { ts } );
    oc.add_output( y,      {t}, { ts } );
    oc.add_output( u,      {t}, { ts } );
    oc.add_output( v,      {t}, { ts } );
    oc.add_output( lam,    {t}, { ts } );
    oc.add_output( OUT_G0, {t}, { ts } );
    oc.add_output( OUT_G1, {t}, { ts } );
    oc.add_output( OUT_E,  {t}, { ts } );
  }

  oc.options.REDUCE.ORDER     = red;
  oc.options.REDUCE.HIDDEN_IC = true;     // the model declares the free initial data; the hidden levels are rows
  oc.options.CLASSIFY.MODE         = OCFESLV::Options::CLASS_AUTO;
  oc.options.INTERFACE.TYPE   = OCFESLV::Options::IC_AUTO;
  oc.options.INTERFACE.IMPOSITION  = imp;
  oc.options.FATAL.REDUCED_DOF = false;   // inspect the audit rather than hard-fail
  oc.options.SOLVE.MARCHING   = marching;
  oc.options.SOLVE.MAX_ITER   = 100;
  oc.options.SOLVE.RES_TOL    = 1.0e-9;
  oc.options.DISPLAY_LEVEL    = 0;
#if defined(CRONOS__WITH_SPQR)
  oc.options.SOLVE.FACTORIZATION = OCFESLV::Options::SOLVE_SPQR;
#endif

  try { R.setup_ok = oc.setup(); }
  catch( ... ) { R.setup_ok = false; std::cout << "  setup() THREW\n"; }

  R.index    = oc.pde_type().differential_index;
  { auto const& plan = oc.reduction_plan();
    R.reduced = !plan.empty(); R.max_index = plan.max_index;
    R.resolved = plan.resolved; R.nAssign = plan.assigns.size(); }
  { OCFESLV::t_DofAudit const& A = oc.reduced_dof_audit();
    R.audit_ran = A.ran; R.audit_sq = A.square; R.audit_ok = A.ok();
    R.rows = A.rows; R.cols = A.cols; R.rank = A.rank; R.defc = A.deficiency; }

  std::cout << "  setup=" << ( R.setup_ok ? "ok" : "FAIL" )
            << "  differential_index=" << R.index
            << "  reduced=" << ( R.reduced ? "yes" : "no" )
            << " (assigns=" << R.nAssign << ", max_index=" << R.max_index
            << ", resolved=" << ( R.resolved ? "yes" : "no" ) << ")\n";
  std::cout << "  reduced_dof_audit: ran=" << ( R.audit_ran ? "yes" : "no" )
            << " rows=" << R.rows << " cols=" << R.cols << " rank=" << R.rank
            << " deficiency=" << R.defc << " square=" << ( R.audit_sq ? "yes" : "no" )
            << " ok=" << ( R.audit_ok ? "YES" : "no" ) << "\n";

  if( !R.setup_ok ) return R;

  std::vector<double> xv, inp;
  if( !oc.init( xv, inp, nullptr ) ){ std::cout << "  init() FAILED\n"; return R; }
  double* const pinp = inp.empty() ? nullptr : inp.data();

  OCFESLV::SolveReport const rep = oc.solve( xv.data(), pinp, nullptr );
  R.converged = rep.converged;
  std::cout << "  solve: converged=" << ( rep.converged ? "yes" : "no" )
            << " iters=" << rep.iterations
            << " final|r|=" << std::scientific << std::setprecision(3) << rep.final_residual << "\n";
  if( !rep.converged ) return R;

  R.F = oc.val_functions();
  size_t const per = 8;
  if( R.F.size() != per * T_SAMPLE.size() ) return R;

  double const E0 = energy0();
  R.maxg0 = R.maxg1 = R.maxE = R.maxTraj = R.maxLam = 0.;
  for( size_t kk = 0; kk < T_SAMPLE.size(); ++kk ){
    double const ts = T_SAMPLE[kk];
    double const xs = R.F[per*kk+0], ys = R.F[per*kk+1], us = R.F[per*kk+2],
                 vs = R.F[per*kk+3], ls = R.F[per*kk+4],
                 g0 = R.F[per*kk+5], g1 = R.F[per*kk+6], Es = R.F[per*kk+7];
    double th, om, xo, yo, uo, vo, lo; ORC.at_node( ts, th, om ); Oracle::cart( th, om, xo, yo, uo, vo, lo );
    R.maxg0   = std::max( R.maxg0,   std::fabs( g0 ) );
    R.maxg1   = std::max( R.maxg1,   std::fabs( g1 ) );
    R.maxE    = std::max( R.maxE,    std::fabs( Es - E0 ) );
    R.maxLam  = std::max( R.maxLam,  std::fabs( ls - lo ) );
    R.maxTraj = std::max( R.maxTraj, std::fabs( xs - xo ) );
    R.maxTraj = std::max( R.maxTraj, std::fabs( ys - yo ) );
    R.maxTraj = std::max( R.maxTraj, std::fabs( us - uo ) );
    R.maxTraj = std::max( R.maxTraj, std::fabs( vs - vo ) );
  }
  return R;
}

// gated evaluation (monolithic): asserts the structural + physical claims
static void gate( ModeOut const& R )
{
  check_true( "S1 differential_index == 3",           R.index == 3 );
  check_true( "S2 reduced (plan non-empty, idx3)",    R.reduced && R.max_index == 3 && R.resolved );
  check_true( "S3 reduced_dof_audit ok (square,rank)",R.audit_ran && R.audit_ok );
  check_true( "C0 de-indexed solve converged",        R.converged );
  if( R.converged ){
    // same pendulum, same resolution as DAE0 -> same realistic drift tolerances
    // (index-3 reduced form; drifts ~2.5e-5 position/velocity, ~2.9e-4 energy,
    //  ~6e-5 trajectory, ~1.4e-3 lam).
    check_le( "P0 position drift max|g0|",   R.maxg0,   1e-4 );
    check_le( "P1 velocity drift max|g1|",   R.maxg1,   1e-4 );
    check_le( "P2 energy drift max|E-E0|",   R.maxE,    1e-3 );
    check_le( "P3 trajectory vs RK4 oracle", R.maxTraj, 2e-4 );
    check_le( "P4 lam vs (u^2+v^2-g y)/L^2", R.maxLam,  5e-3 );
  }
}

int main()
{
  std::cout << "================================================================\n"
            << "  OCFE_DAE1 : INDEX-3 pendulum posed DIRECTLY\n"
            << "  exercises automatic Pantelides index reduction + DoF audit\n"
            << "  L=" << L_rod << " g=" << g_grav << " theta0=" << th0
            << " omega0=" << om0 << " T_end=" << T_end << "\n"
            << "  E0=" << std::scientific << std::setprecision(6) << energy0() << "\n"
            << "================================================================\n";

  std::cout << "\n#### MONOLITHIC (gated: reduction + DoF audit + solve) ####\n";
  ModeOut const wm = run_mode( OCFESLV::Options::IC_WEAK,   "WEAK",   false );
  std::cout << "\n  [gate IC_WEAK]\n";   gate( wm );
  ModeOut const sm = run_mode( OCFESLV::Options::IC_STRONG, "STRONG", false );
  std::cout << "\n  [gate IC_STRONG]\n"; gate( sm );

  std::cout << "\n  [cross-mode]\n";
  // Both monolithic reductions are gated against the RK4 oracle in gate() above; the WEAK
  // vs STRONG reduced solutions differ only at the discretisation scale -> print the delta.
  if( wm.converged && sm.converged && wm.F.size() == sm.F.size() && !wm.F.empty() ){
    double e = 0.; for( size_t i = 0; i < wm.F.size(); ++i ) e = std::max( e, std::fabs( wm.F[i] - sm.F[i] ) );
    std::cout << "  XM WEAK reduced ~ STRONG reduced  delta=" << std::scientific
              << std::setprecision(3) << e << "  (discretisation-scale; informational)\n";
  }
  else check_true( "XM WEAK/STRONG cross-mode available (both converged)", wm.converged && sm.converged );

  // ---- marching: INFORMATIONAL only (not gated) ----
  std::cout << "\n#### MARCHING (informational; index reduction is a monolithic analysis) ####\n";
  ModeOut const wc = run_mode( OCFESLV::Options::IC_WEAK, "WEAK", true );
  std::cout << "  [info] marching: setup=" << ( wc.setup_ok ? "ok" : "FAIL" )
            << " index=" << wc.index << " reduced=" << ( wc.reduced ? "yes" : "no" )
            << " audit_ok=" << ( wc.audit_ok ? "yes" : "no" )
            << " converged=" << ( wc.converged ? "yes" : "no" );
  if( wc.converged )
    std::cout << " maxTraj=" << std::scientific << std::setprecision(3) << wc.maxTraj;
  std::cout << "  (not gated)\n";

  // ---------------------------------------------------------------------------------------------------------
  // 2026-09-19 -- THE SYSTEMATIC MATRIX: 3 impositions x 2 reductions, monolithic.  DAE1 carries auto-closure,
  // a high index and marching, and its TRACE cells had never been run under RED_MAIN.  It is also one of the
  // three programs the rev273 multiplier experiment broke through the DOF audit, so its exact modes are worth
  // exercising on both reductions.  A cell expected to differ gets a documented XFAIL naming the mechanism.
  // ---------------------------------------------------------------------------------------------------------
  std::cout << "\n---- SYSTEMATIC MATRIX (monolithic): imposition x reduction ----\n";
  {
    struct MCell { char const* imp; OCFESLV::Options::ImpositionType it;
                   char const* red; OCFESLV::Options::ReductionType rt; };
    MCell const cells[6] = {
      { "WEAK",   OCFESLV::Options::IC_WEAK,   "RED_FULL", OCFESLV::Options::RED_FULL },
      { "TRACE",  OCFESLV::Options::IC_TRACE,  "RED_FULL", OCFESLV::Options::RED_FULL },
      { "STRONG", OCFESLV::Options::IC_STRONG, "RED_FULL", OCFESLV::Options::RED_FULL },
      { "WEAK",   OCFESLV::Options::IC_WEAK,   "RED_MAIN", OCFESLV::Options::RED_MAIN },
      { "TRACE",  OCFESLV::Options::IC_TRACE,  "RED_MAIN", OCFESLV::Options::RED_MAIN },
      { "STRONG", OCFESLV::Options::IC_STRONG, "RED_MAIN", OCFESLV::Options::RED_MAIN } };
    bool matrix_ok = true;
    std::cout << "  " << std::left << std::setw(10) << "IMPOSITION" << std::setw(11) << "REDUCTION"
              << std::setw(8) << "setup" << std::setw(7) << "conv" << std::setw(8) << "index" << "VERDICT\n";
    for( auto const& mc : cells ){
      ModeOut m = run_mode( mc.it, mc.imp, false, mc.rt );
      bool const cell  = m.setup_ok && m.converged;
      bool const xfail = false;
      if( !cell && !xfail ) matrix_ok = false;
      std::cout << "  " << std::left << std::setw(10) << mc.imp << std::setw(11) << mc.red
                << std::setw(8) << ( m.setup_ok ? "ok" : "FAILED" ) << std::setw(7) << ( m.converged ? "y" : "n" )
                << std::setw(8) << m.index
                << ( cell ? ( xfail ? "PASS (unexpected: the documented defect is gone?)" : "PASS" )
                          : ( xfail ? "XFAIL (documented)" : "FAIL" ) ) << "\n";
    }
    check_true( "matrix: all 6 imposition x reduction cells set up and converge", matrix_ok );
    std::cout << "  MATRIX: " << ( matrix_ok ? "PASS" : "FAIL" ) << "\n";
  }

  std::cout << "\n================================================================\n"
            << "  OCFE_DAE1: " << g_pass << " passed, " << g_fail << " failed -- "
            << ( g_fail == 0 ? "ALL PASS" : "SOME FAILED" ) << "\n"
            << "================================================================\n";
  return g_fail == 0 ? 0 : 1;
}
