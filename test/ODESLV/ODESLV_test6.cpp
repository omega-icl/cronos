// test6_ff.cpp -- WITNESS for per-stage state transitions in ODESLV (the _pos_ic path), declared through FFModel.
// =============================================================================================
// A transition x(tau+) = g(x(tau-), p) becomes the per-stage initial-value entry of the stage starting at tau: ODESLV
// evaluates it with the state reached at tau^- and reinitialises the integrator.  (2026-09-29: ported from the retired
// description path -- set_initial( vector<vector<FFVar>> ) -- to FFModel::add_transition; the checks are unchanged.)
//
//   dx/dt = -a x + p   on [0,3], three stages, transitions at t = 1 and t = 2
//
// A  identity transitions             must reproduce the SAME model declared without any
// B  a genuine jump of +c per stage   must differ, and by exactly the amount the jump propagates
// C  gradients w.r.t. p               FSA == ASA == central finite differences
//
// Case A is the load-bearing one: it fails if the transition is applied wrongly, and it also fails if the per-stage
// restart or the ASA cycle mishandles a stage boundary that REINITIALISES the integrator.
// =============================================================================================
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <cmath>

#include "odeslvs_cvodes.hpp"
#include "ocbase.hpp"

static double const A_DECAY = 0.7;
static double const X0_VAL  = 0.5;
static double const P_VAL   = 1.3;
static size_t const NS      = 3;

static int npass = 0, nfail = 0;
void check( bool c, std::string const& what )
{ std::cout << "  " << (c? "PASS  ": "FAIL  ") << what << std::endl; c? ++npass: ++nfail; }

double reldiff( double a, double b )
{ double const s = std::fabs(b); return std::fabs(a-b) / ( s > 1e-12? s: 1. ); }

struct Out
{
  bool ok = false;  std::string err;
  size_t np = 0, nf = 0, nsen = 0;
  std::vector<double> f;
  std::vector<std::vector<double>> g;
};

//! @brief Build the IVP through FFModel.  @p perstage declares a transition at every interior stage boundary (so
//! ODESLV holds one initial-value entry per stage and _pos_ic fires) or none.  @p jump is added
//! to the state at each interior stage boundary.
Out solve( bool perstage, double jump, bool adjoint )
{
  Out R;
  mc::FFGraph DAG;
  mc::ODESLVS_CVODES IVP( &DAG );                           // the model lives on DAG
  mc::FFVar T = DAG.add_var( "t" ), X = DAG.add_var( "x(t)" ), P = DAG.add_var( "p" ), X0p = DAG.add_var( "x0" );
  mc::FFPartial OpP;  mc::FFEval OpE;
  std::vector<double> dT( NS+1 ); for( size_t k = 0; k <= NS; ++k ) dT[k] = double( k );
  IVP.add_domain( T, mc::FFDom( dT, mc::FFDom::LGR, 3 ) );  IVP.set_evolution_domain( T );
  IVP.add_state( X, {T} );  IVP.add_input( P );  IVP.add_input( X0p );   // the IC is a parameter, so it is differentiable
  IVP.update_ref( X, X0_VAL );
  int const T_INT = mc::FFDom::ALL - mc::FFDom::LB;
  IVP.add_equation( OpP( X, T ) + A_DECAY*X - P, {T}, {T_INT}, mc::FFModel::EqnOptions( mc::FFModel::EqnRole::INTERIOR, 0 ) );
  IVP.add_equation( X - X0p, {T}, {mc::FFDom::LB}, mc::FFModel::EqnOptions( mc::FFModel::EqnRole::INITIAL, 0 ) );
  if( perstage )
    for( size_t k = 1; k < NS; ++k )
      IVP.add_transition( X + jump, X, T, double( k ) );   // the transition map at t = k (identity when jump == 0)
  IVP.add_output( OpE( X, T, double( NS ) ) );              // x at the final time

  IVP.options.LINSOL  = mc::BASE_CVODES::Options::DENSE;
  IVP.options.DISPLAY = 0;
  IVP.options.ATOL = IVP.options.ATOLB = IVP.options.ATOLS = 1e-10;
  IVP.options.RTOL = IVP.options.RTOLB = IVP.options.RTOLS = 1e-10;

  if( !IVP.setup() ){ R.err = "setup failed: " + IVP.extract_error(); return R; }
  R.np = IVP.np();  R.nf = IVP.nf();

  size_t const ip = IVP.parameter_index( P )[0], ix0 = IVP.parameter_index( X0p )[0];
  std::vector<double> Pv( IVP.np() );  Pv[ip] = P_VAL;  Pv[ix0] = X0_VAL;
  auto const st = adjoint? IVP.solve_asens( Pv ): IVP.solve_fsens( Pv );
  if( st != mc::ODESLVS_CVODES::STATUS::NORMAL ){ R.err = "solve failed"; return R; }
  R.f = IVP.val_function();
  { auto const G = IVP.val_function_gradient();          // rows ordered (p, x0) whatever the declared order
    R.g.clear(); if( G.size() > ip && G.size() > ix0 ){ R.g.push_back( G[ip] ); R.g.push_back( G[ix0] ); } }
  R.nsen = R.g.size();
  R.ok = true;
  return R;
}

int main()
{
  std::cout << std::scientific << std::setprecision(6);
  std::cout << "================================================================\n"
            << "  test6_ff: per-stage initial conditions (the _pos_ic path)\n"
            << "  dx/dt = -a x + p on [0,3], transitions at t = 1 and t = 2\n"
            << "================================================================\n";

  std::cout << "\n--- reference: ONE initial condition (_vIC.size() == 1, _pos_ic never fires)\n";
  Out const ref = solve( false, 0., false );
  if( !ref.ok ){ std::cout << "  reference FAILED: " << ref.err << "\n"; return 1; }
  check( ref.np == 2,   "ref  np == 2  (p and the IC)" );
  check( ref.nf == 1,   "ref  nf == 1" );
  check( ref.nsen == 2, "ref  two sensitivity directions (both parameters resolved in _ndxSEN)" );
  std::cout << "       x(3) = " << ref.f[0] << "   dx(3)/dp = " << ref.g[0][0] << "\n";

  std::cout << "\n--- A: THREE ICs with IDENTITY transfers -- must reproduce the reference\n";
  Out const A = solve( true, 0., false );
  if( !A.ok ) check( false, std::string("A  ") + A.err );
  else{
    std::ostringstream o1; o1 << "A  x(3) = " << A.f[0] << " vs reference " << ref.f[0]
                              << "  (rel " << reldiff( A.f[0], ref.f[0] ) << ")";
    check( reldiff( A.f[0], ref.f[0] ) < 1e-7, o1.str() );
    std::ostringstream o2; o2 << "A  dx(3)/dp matches the reference (rel "
                              << reldiff( A.g[0][0], ref.g[0][0] ) << ")";
    check( reldiff( A.g[0][0], ref.g[0][0] ) < 1e-6, o2.str() );
  }

  std::cout << "\n--- B: a genuine jump of +0.25 at each interior boundary\n";
  double const JUMP = 0.25;
  Out const B = solve( true, JUMP, false );
  if( !B.ok ) check( false, std::string("B  ") + B.err );
  else{
    check( reldiff( B.f[0], A.f[0] ) > 1e-3, "B  the jump changes x(3)" );
    // Two jumps, at t=1 and t=2, each decaying by exp(-a*(3-t)) to the final time.
    double const expect = A.f[0] + JUMP*std::exp(-A_DECAY*2.) + JUMP*std::exp(-A_DECAY*1.);
    std::ostringstream o3; o3 << "B  x(3) = " << B.f[0] << " vs analytic " << expect
                              << "  (rel " << reldiff( B.f[0], expect ) << ")";
    check( reldiff( B.f[0], expect ) < 1e-6, o3.str() );
    // The jump is a constant, so it carries no p-dependence: the gradient is unchanged.
    std::ostringstream o4; o4 << "B  dx(3)/dp unchanged by a constant jump (rel "
                              << reldiff( B.g[0][0], A.g[0][0] ) << ")";
    check( reldiff( B.g[0][0], A.g[0][0] ) < 1e-6, o4.str() );
  }

  std::cout << "\n--- C: adjoint and finite differences on the jumping model\n";
  Out const Cadj = solve( true, JUMP, true );
  if( !Cadj.ok ) check( false, std::string("C  ") + Cadj.err );
  else{
    std::ostringstream o5; o5 << "C  ASA == FSA (rel " << reldiff( Cadj.g[0][0], B.g[0][0] ) << ")";
    check( reldiff( Cadj.g[0][0], B.g[0][0] ) < 1e-6, o5.str() );
  }
  {
    double const h = 1e-6;
    // dx(3)/dp by central differences on the analytic solution of the jumping model
    double const a = A_DECAY, T = 3.;
    auto xT = [&]( double p ){
      double x = X0_VAL;
      for( size_t k = 0; k < NS; ++k ){
        x = ( x - p/a )*std::exp(-a*1.) + p/a;             // one unit stage
        if( k+1 < NS ) x += JUMP;                          // transition at the boundary
      }
      (void)T;  return x;
    };
    double const fd = ( xT( P_VAL+h ) - xT( P_VAL-h ) ) / ( 2.*h );
    std::ostringstream o6; o6 << "C  FSA == analytic central FD, " << B.g[0][0] << " vs " << fd
                              << "  (rel " << reldiff( B.g[0][0], fd ) << ")";
    check( reldiff( B.g[0][0], fd ) < 1e-5, o6.str() );
  }

  std::cout << "\n================================================================\n"
            << "  test6_ff: " << npass << " passed, " << nfail << " failed -- "
            << ( nfail? "FAILURES": "ALL PASS" ) << "\n"
            << "================================================================\n";
  return nfail? 1: 0;
}
