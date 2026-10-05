// ===========================================================================
//  OCFE_nlpartial.cpp  --  reduce_order: FFPartial nested in a NONLINEAR
//                          expression (division by a derivative)
//
//  THE OPEN ISSUE
//  --------------
//  reduce_order already peels high-order partials, and already materialises the
//  partials carried INSIDE an OpI/OpEval operand (materialize_operand_partials),
//  because a nonlocal reduction LINK cannot expose an inner OpP.  It did NOT
//  touch a first-order OpP(state,dir) sitting in an ordinary residual.  That is
//  correct when the residual is AFFINE in the partial -- `D*u_z`, `u*u_z`,
//  `exp(-u)*u_z` all classify and assemble fine.  It is NOT correct when the
//  derivative is wrapped NONLINEARLY, which is exactly what an ALE / moving-mesh
//  formulation does:
//
//      u_z  = u_xi / x_xi                 <- DIVISION by a derivative
//      u_zz = u_xixi/x_xi^2 - x_xixi u_xi/x_xi^3
//
//  Every moving-mesh driver so far had to work around this by hand-rolling the
//  physical derivatives as extra states with hand-written defining equations in
//  MULTIPLY form (q*x_xi - u_xi = 0, xg - x_xi = 0), because a hand-rolled
//  auxiliary is invisible to the classifier AS an auxiliary -- which is how the
//  multiply-defined flux states and the rectangular symbols arose.
//
//  THE FIX UNDER TEST  (options.REDUCE.NONLINEAR_PARTIALS)
//  ------------------------------------------------------
//  A first-order partial of a bare state is materialised as a derivative
//  auxiliary Dp with a defining LINK
//      OpP(state,dir) - Dp = 0
//  iff the residual is not AFFINE in THAT partial alone (mc::FFDep type != L) --
//  the same mint / LINK / extra_subst_map-reuse protocol the OpI/OpEval operand
//  path already uses.  The nonlinearity then closes over a STATE, and the raw
//  OpP survives only inside its LINK, which is an ordinary differential equation
//  the principal-symbol proxy resolves.
//
//  The trigger is per-partial and therefore minimal.  In
//      (c - x_t) * u_xi / x_xi
//  ONLY x_xi is materialised; x_t and u_xi are each affine on their own and stay
//  ordinary variable-coefficient partials.  Any residual that is affine in each
//  of its first-order partials -- the whole existing corpus -- is untouched.
//
//  WHAT THIS DRIVER DOES
//  ---------------------
//  A STEADY 1-D MMS advection-diffusion on a SOLVED equidistribution mesh,
//  written in the NATURAL form -- divisions by derivatives left in place, no
//  hand-rolled flux states.  Two states only: u and x.
//
//      MESH:  W * x_xixi + W_xi * x_xi = 0            (equidistribution)
//      PDE :  c*(u_xi/x_xi)
//             - D*( u_xixi/x_xi^2 - x_xixi*u_xi/x_xi^3 ) - f = 0
//      x(0)=0, x(1)=L,  u(0), u(1) from the exact solution.
//
//  W = 1/x*_xi (prescribed), so the equidistribution solution is EXACTLY the
//  sinh oracle map x*(xi) and the check against it is sharp.  f is the MMS
//  source that makes u*(z) = 0.5(1-tanh((z-s)/delta)) exact.
//
//  It is run TWICE, and the two runs are the test:
//      REDUCE_NONLINEAR_PARTIALS = false  -> expected to FAIL (this is the bug)
//      REDUCE_NONLINEAR_PARTIALS = true   -> expected to PASS
//  A single build, one flag apart, so nothing else can explain a difference.
//
//  Build:
//    g++ -std=c++17 <suite flags> -DOCFE_OCFESLV_HEADER='"ocfeslv_nlpart.hpp"'
//        OCFE_nlpartial.cpp -o OCFE_nlpartial <libs>
//  Run:
//    ./OCFE_nlpartial              # Pe=100, nel = 4, 8, 16
//    ./OCFE_nlpartial -v 100 8     # DISPLAY_LEVEL 1 (shows what is materialised)
// ===========================================================================

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>

#ifndef OCFE_OCFESLV_HEADER
#define OCFE_OCFESLV_HEADER "ocfeslv_nlpart.hpp"
#endif
#include OCFE_OCFESLV_HEADER

using namespace mc;

static double const L_DOM = 1.0, C_ADV = 1.0, S_FRONT = 0.5;
static double const ERRTOL = 5.0e-4;    // u accuracy
static double const XTOL   = 1.0e-6;    // mesh must match the oracle map
static size_t const NND_XI = 8;
// Accuracy gates are DISCRETISATION statements and need enough resolution to mean
// anything.  With delta=1e-2 a 4-element mesh cannot resolve the front: errNode is
// ~2e-2 there BY CONSTRUCTION (the solve itself converges to max|r|~1e-12, square,
// 2 Newton iterations).  Gating |x-x*| and errNode at nel=4 therefore reports a
// framework failure that is really a resolution statement -- so below NEL_ACC the
// row is printed for the convergence table but not gated.  Measured order of
// convergence here is ~8 for |x-x*| and ~6 for errNode, i.e. spectral, as expected.
static size_t const NEL_ACC = 8;

struct Params
{
  double L = L_DOM, c = C_ADV, s = S_FRONT, Pe = 1.0e2;
  double D()     const { return c * L / Pe; }
  double delta() const { return D() / c; }
  double dm()    const { return delta(); }
};

// --- oracle equidistribution map x*(xi) (sinh, clustered at the front) -----
struct Oracle
{
  Params p;
  double a() const { return std::asinh( p.s / p.dm() )
                          + std::asinh( ( p.L - p.s ) / p.dm() ); }
  double b() const { return std::asinh( p.s / p.dm() ) / a(); }
  double x   ( double xi ) const { return p.s + p.dm() * std::sinh( a() * ( xi - b() ) ); }
  double x_xi( double xi ) const { return a() * p.dm() * std::cosh( a() * ( xi - b() ) ); }
};

static double u_exact( double z, Params const& p )
  { return 0.5 * ( 1.0 - std::tanh( ( z - p.s ) / p.delta() ) ); }

struct Result
{
  bool   setup_ok = false, square = false, conv = false, threw = false;
  std::string status = "(not reached)", exmsg;
  size_t nVar = 0, nEqn = 0, nAux = 0;
  double errNode = std::numeric_limits<double>::infinity();
  double errMesh = std::numeric_limits<double>::infinity();
};

// ---------------------------------------------------------------------------
static Result run( Params const& p, size_t nel, int display, bool nl_partials )
{
  Result R;
  Oracle O{ p };

  FFGraph DAG;
  FFVar xi = DAG.add_var( "xi" );
  FFVar u  = DAG.add_var( "u(xi)" );
  FFVar x  = DAG.add_var( "x(xi)" );
  FFPartial OpP;

  double const D = p.D(), del = p.delta(), dm = p.dm();

  // --- prescribed monitor W(xi) = 1/x*_xi(xi), and its xi-derivative --------
  auto ash = []( FFVar const& y ){ return log( y + sqrt( y * y + 1.0 ) ); };
  double const a1 = std::asinh( p.s / dm );
  double const a2 = std::asinh( ( p.L - p.s ) / dm );
  double const aa = a1 + a2, bb = a1 / aa;
  FFVar const yv        = aa * ( xi - bb );
  FFVar const xstar_xi  = aa * dm * cosh( yv );
  FFVar const xstar_xx  = aa * aa * dm * sinh( yv );
  FFVar const W         = 1.0 / xstar_xi;
  FFVar const W_xi      = -xstar_xx / ( xstar_xi * xstar_xi );
  (void)ash;

  // --- NATURAL ALE form: divisions by derivatives, no hand-rolled states ----
  FFVar const Xxi   = OpP( x, xi );
  FFVar const Xxixi = OpP( x, { xi, 2 } );
  FFVar const Uxi   = OpP( u, xi );
  FFVar const Uxixi = OpP( u, { xi, 2 } );

  FFVar const Uz  = Uxi / Xxi;                                       // u_z
  FFVar const Uzz = Uxixi / ( Xxi * Xxi )                            // u_zz
                  - Xxixi * Uxi / ( Xxi * Xxi * Xxi );

  // MMS source at the SOLVED physical point x, so u*(z) is exact
  FFVar const TT  = tanh( ( x - p.s ) / del );
  FFVar const uz_ = -0.5 / del * ( 1.0 - TT * TT );
  FFVar const uzz_ = ( 1.0 / ( del * del ) ) * TT * ( 1.0 - TT * TT );
  FFVar const FSR = p.c * uz_ - D * uzz_;

  FFVar PDE  = p.c * Uz - D * Uzz - FSR;
  FFVar MESH = W * Xxixi + W_xi * Xxi;                 // d/dxi( W x_xi ) = 0
  FFVar U_BC = u - 0.5 * ( 1.0 - tanh( ( x - p.s ) / del ) );
  FFVar X_LB = x - 0.0;
  FFVar X_UB = x - p.L;

  OCFESLV oc( &DAG );
  oc.add_domain( xi, FFDom( 0., 1.0, nel, FFDom::LGL, NND_XI ) );
  oc.add_state ( u, { xi } );
  oc.add_state ( x, { xi } );

  { Oracle const Ov = O; Params const pv = p; FFVar xv = xi;
    oc.update_ref( x, [Ov,xv]( OCFESLV::t_Coord const& q ){ return Ov.x( q.at( xv ) ); } );
    oc.update_ref( u, [Ov,pv,xv]( OCFESLV::t_Coord const& q ){
        return u_exact( Ov.x( q.at( xv ) ), pv ); } );
  }

  OCFESLV::EqnOptions int_opt( OCFESLV::EqnRole::INTERIOR, 0 );
  OCFESLV::EqnOptions bnd_opt( OCFESLV::EqnRole::BOUNDARY, 0 );
  int const XI_INT = FFDom::ALL - FFDom::LB - FFDom::UB;

  oc.add_equation( PDE,   { xi }, { XI_INT    }, int_opt );
  oc.add_equation( U_BC,  { xi }, { FFDom::LB }, bnd_opt );
  oc.add_equation( U_BC,  { xi }, { FFDom::UB }, bnd_opt );
  oc.add_equation( MESH,  { xi }, { XI_INT    }, int_opt );
  oc.add_equation( X_LB,  { xi }, { FFDom::LB }, bnd_opt );
  oc.add_equation( X_UB,  { xi }, { FFDom::UB }, bnd_opt );

  oc.options.REDUCE.ORDER              = OCFESLV::Options::RED_FULL;
  oc.options.REDUCE.NONLINEAR_PARTIALS = nl_partials;      // <-- the one variable
  oc.options.CLASSIFY.MODE                  = OCFESLV::Options::CLASS_AUTO;
  oc.options.SOLVE.MAX_ITER            = 60;
  oc.options.SOLVE.RES_TOL             = 1.0e-9;
  oc.options.DISPLAY_LEVEL             = display;

  try { R.setup_ok = oc.setup(); }
  catch( OCFESLV::Exceptions& e ){
    R.threw = true; R.exmsg = "OCFESLV ierr=" + std::to_string( e.ierr() ); }
  catch( FFBase::Exceptions& e ){
    R.threw = true; R.exmsg = std::string( "FFBase: " ) + e.what(); }
  catch( std::exception& e ){ R.threw = true; R.exmsg = e.what(); }
  catch( ... ){ R.threw = true; R.exmsg = "unknown exception"; }

  R.status = OCFESLV::setup_status_str( oc.setup_status() );
  try { R.nVar = oc.n_colloc_sta(); R.nEqn = oc.n_colloc_eqn(); } catch(...){}
  try { R.nAux = oc.auxiliary_states().size(); } catch(...){}
  R.square = ( R.nVar && R.nVar == R.nEqn );
  if( !R.setup_ok || !R.square ) return R;

  std::vector<double> varInit, inpInit;
  if( !oc.init( varInit, inpInit, nullptr ) ) return R;
  std::vector<double> xv = varInit;
  OCFESLV::SolveReport rep = oc.solve( xv.data(), nullptr, nullptr );
  R.conv = rep.converged;

  double emax = 0., xmax = 0.;
  int const PPE = 6;
  for( size_t e = 0; e < nel; ++e )
    for( int s = 0; s <= PPE; ++s ){
      double const qq = ( double( e ) + double( s ) / double( PPE ) ) / double( nel );
      OCFESLV::t_Coord pt; pt[xi] = qq;
      double un = 0., xn = 0.;
      try { un = oc.eval_colloc<double>( u, pt, xv.data(), nullptr, nullptr );
            xn = oc.eval_colloc<double>( x, pt, xv.data(), nullptr, nullptr ); }
      catch( ... ){ continue; }
      if( !std::isfinite( un ) || !std::isfinite( xn ) ){
        emax = std::numeric_limits<double>::infinity(); continue; }
      xmax = std::max( xmax, std::fabs( xn - O.x( qq ) ) );
      emax = std::max( emax, std::fabs( un - u_exact( O.x( qq ), p ) ) );
    }
  R.errNode = emax; R.errMesh = xmax;
  return R;
}

// ---------------------------------------------------------------------------
int main( int argc, char** argv )
{
  double Pe = 100.0; bool pe_set = false; int display = 0;
  std::vector<size_t> nelList;
  for( int i = 1; i < argc; ++i ){
    if( !std::strcmp( argv[i], "-v" ) ){ display = 1; continue; }
    if( !pe_set && nelList.empty() && std::atof( argv[i] ) >= 10.0 )
      { Pe = std::atof( argv[i] ); pe_set = true; continue; }
    nelList.push_back( (size_t)std::atol( argv[i] ) );
  }
  if( nelList.empty() ) nelList = { 4, 8, 16 };
  Params p; p.Pe = Pe;

  std::cout << "================================================================\n"
            << "  reduce_order: FFPartial nested in a NONLINEAR expression\n"
            << "  NATURAL ALE form -- divisions by derivatives, 2 states (u,x),\n"
            << "  NO hand-rolled flux states:\n"
            << "     PDE : c*(u_xi/x_xi) - D*( u_xixi/x_xi^2\n"
            << "                             - x_xixi*u_xi/x_xi^3 ) - f = 0\n"
            << "     MESH: W*x_xixi + W_xi*x_xi = 0,  W = 1/x*_xi (prescribed)\n"
            << "  Same build, run twice, ONE option apart:\n"
            << "     REDUCE_NONLINEAR_PARTIALS = false -> the open issue\n"
            << "     REDUCE_NONLINEAR_PARTIALS = true  -> the fix\n"
            << "  Pe=" << Pe << "  delta=" << std::scientific << std::setprecision(2)
            << p.delta() << "  xi order " << NND_XI << "\n"
            << "  gates: setup square, solve converges, |x-x*|<" << XTOL
            << ", errNode<" << ERRTOL << "\n"
            << "================================================================\n";
  std::cout << "  " << std::left
            << std::setw(6) << "nel" << std::setw(10) << "NLPART"
            << std::setw(8) << "nVar" << std::setw(6) << "nAux"
            << std::setw(8) << "square" << std::setw(7) << "conv"
            << std::setw(12) << "|x-x*|" << std::setw(12) << "errNode"
            << "  status\n";

  bool all_ok = true, any_effect = false;
  for( size_t ne : nelList ){
    Result Roff, Ron;
    for( int k = 0; k < 2; ++k ){
      bool const nl = ( k == 1 );
      Result R = run( p, ne, display, nl );
      bool const resolved = ( ne >= NEL_ACC );
      bool const ok = R.setup_ok && R.square && R.conv
                   && std::isfinite( R.errMesh ) && std::isfinite( R.errNode )
                   && ( !resolved || ( R.errMesh < XTOL && R.errNode < ERRTOL ) );
      // Only the NLPART=on runs are gated; the off runs document the defect.
      if( nl ) all_ok &= ok;
      if( nl ) Ron = R; else Roff = R;
      std::cout << "  " << std::left
                << std::setw(6) << ne << std::setw(10) << ( nl ? "on" : "off" )
                << std::setw(8) << R.nVar << std::setw(6) << R.nAux
                << std::setw(8) << ( R.square ? "yes" : "NO" )
                << std::setw(7) << ( R.conv ? "yes" : "no" )
                << std::scientific << std::setprecision(2)
                << std::setw(12) << R.errMesh << std::setw(12) << R.errNode
                << "  " << R.status;
      if( R.threw ) std::cout << "  THREW: " << R.exmsg;
      std::cout << ( ok ? "   [OK]" : "   [--]" )
                << ( resolved ? "" : "  (accuracy not gated: nel<NEL_ACC)" ) << "\n";
    }

    // Does the option actually CHANGE anything here?  It need not: when the model
    // also carries u_xixi / x_xixi, RED_FULL's high-order peel already mints the
    // same first-order auxiliary and its extra_subst_map reuse folds the divisions
    // onto it -- both routes converge on an identical reduced system.  Say so,
    // loudly, rather than let an [OK] imply the feature was exercised.
    // Compare the reduced SYSTEM, not the error norms.  nVar/nAux is the exact,
    // noise-free signal for "did the option change anything"; the error norms are
    // solver output and at fine resolution sit at residual level (errMesh 4e-10 vs
    // a 1.4e-10 residual), so no relative tolerance on them can be calibrated.
    bool const same = ( Roff.nVar == Ron.nVar ) && ( Roff.nAux == Ron.nAux );
    any_effect = any_effect || !same;
    if( same )
      std::cout << "  " << std::left << std::setw(6) << ne
                << "  [note] option had NO effect: identical nVar/nAux/errors."
                   "  RED_FULL's high-order peel\n         already mints this"
                   " auxiliary and reuses it in the divisions -- this case does"
                   " NOT isolate the feature.\n";
  }

  std::cout << "\n================================================================\n"
            << "  nlpartial: " << ( all_ok ? "PASS" : "FAIL" )
            << ( any_effect ? "   (option changed the reduced system)"
                            : "   -- BUT THE OPTION NEVER CHANGED ANYTHING:"
                              " this driver does not isolate the feature" ) << "\n"
            << "  PASS means the NATURAL ALE form -- derivatives divided by\n"
            << "  derivatives, only u and x declared -- sets up square, solves,\n"
            << "  and reproduces the oracle mesh and the MMS solution, with\n"
            << "  reduce_order minting the derivative auxiliaries itself.\n"
            << "  The NLPART=off rows are the defect they replace: they are NOT\n"
            << "  gated, they are printed so the difference is visible in one run.\n"
            << "  Run -v to see which partials are materialised.\n"
            << "================================================================\n";
  return all_ok ? 0 : 1;
}
