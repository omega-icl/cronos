// OCFE_hypclosure.cpp -- rev317: the hyperbolic boundary closure appends only the per-face DEFICIT (2026-09-21).
//
// A model may close an outflow face itself by extending its own PDE's mask to that face -- the idiom of
// OCFE_PDE6/7/8/10.  Before rev317 the automatic closure then added its own characteristic row at the same face,
// double-closing it, so those drivers had to switch AUTO.HYP_CLOSURE off.  This driver checks the mechanism
// directly (PDE6/7/8/10 need a helper header that is not part of the project files):
//   K0  scalar advection u_t + a u_z = f (a>0): the closure row IS the PDE, so auto-only, manual-only and
//       manual+auto must give the SAME solution;
//   K1  the 2x2 system of OCFE_PDE15 (speeds +-a*sqrt(g)) closed the PDE7 way -- the first PDE extended to LB,
//       the second to UB.  Neither row is the characteristic combination, so coverage rests on the projection
//       test (their weights projected onto the outgoing rowspace);
//   K2  a decoupled pair with BOTH characteristics leaving at UB (n_out = 2), only one PDE extended: PARTIAL
//       coverage -- the closure must append exactly the missing direction.
// Built against rev316 the manual+auto variants are double-closed (not square); against rev317 they must match
// the manual-only reference.
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <string>
#include <algorithm>
#ifndef OCFE_OCFESLV_HEADER
#define OCFE_OCFESLV_HEADER "ocfeslv.hpp"
#endif
#include OCFE_OCFESLV_HEADER
using namespace mc;

static int g_pass = 0, g_fail = 0;
static void check( std::string const& nm, bool ok, std::string const& detail = "" )
{
  ( ok ? g_pass : g_fail )++;
  std::cout << "  " << std::left << std::setw(70) << nm << std::right << ( ok ? "PASS" : "FAIL" )
            << ( detail.empty() ? "" : "   " + detail ) << "\n";
}
static std::string sci( double v ){ char b[32]; std::snprintf( b, sizeof b, "%.3e", v ); return b; }

struct Run { bool setup = false, square = false, solved = false; size_t nVar = 0, nEqn = 0; std::vector<double> var; };

// kind 0: scalar advection; 1: 2x2 system (PDE15); 2: decoupled pair, both outgoing at UB.
// manual: bit0 extends the first PDE's mask to its outflow face, bit1 the second PDE's.
static Run run( int kind, bool auto_on, int manual )
{
  Run R;
  FFGraph DAG; OCFESLV oc( &DAG ); FFPartial OpP;
  FFVar t = DAG.add_var( "t" ), z = DAG.add_var( "z" );
  FFVar u = DAG.add_var( "u(t,z)" ), v = DAG.add_var( "v(t,z)" );
  oc.add_domain( t, FFDom( 0., 1., 2, FFDom::LGR, 3 ) );
  oc.add_domain( z, FFDom( 0., 1., 3, FFDom::LGL, 4 ) );
  oc.add_state( u, { t, z } );
  if( kind > 0 ) oc.add_state( v, { t, z } );
  typedef OCFESLV::EqnOptions EO; typedef OCFESLV::EqnRole ER;
  EO io( ER::INTERIOR, 0 ), ini( ER::INITIAL, 0 ), bnd( ER::BOUNDARY, 0 );
  int const T_INT = FFDom::ALL - FFDom::LB, Z_INT = FFDom::ALL - FFDom::LB - FFDom::UB;
  int const Z_UB = FFDom::ALL - FFDom::LB, Z_LB = FFDom::ALL - FFDom::UB;   // interior + one face
  double const a = 1.0, g = 4.0, sg = std::sqrt( g ), b = 0.5, k = 1.3, al = 0.4;
  // manufactured solutions and their partials, written out (as OCFE_PDE15 does) rather than via OpP of an expression
  FFVar U  = sin( k * z + 0.2 ) * ( 1.0 + al * t ),  Ut = al * sin( k * z + 0.2 ),  Uz =  k * cos( k * z + 0.2 ) * ( 1.0 + al * t );
  FFVar V  = cos( k * z + 0.7 ) * ( 1.0 + al * t ),  Vt = al * cos( k * z + 0.7 ),  Vz = -k * sin( k * z + 0.7 ) * ( 1.0 + al * t );
  if( kind == 0 ){
    FFVar PDE = OpP( u, t ) + a * OpP( u, z ) - ( Ut + a * Uz );
    oc.add_equation( PDE,   { t, z }, { T_INT, ( manual & 1 ) ? Z_UB : Z_INT }, io );
    oc.add_equation( u - U, { t, z }, { FFDom::LB, FFDom::ALL }, ini );
    oc.add_equation( u - U, { t, z }, { T_INT, FFDom::LB }, bnd );              // inflow
  }
  else if( kind == 1 ){
    FFVar PDEC = OpP( u, t ) + a * OpP( v, z ) - ( Ut + a * Vz );
    FFVar PDEU = OpP( v, t ) + g * a * OpP( u, z ) - ( Vt + g * a * Uz );
    oc.add_equation( PDEC, { t, z }, { T_INT, ( manual & 1 ) ? Z_LB : Z_INT }, io );   // PDE7 idiom: to LB
    oc.add_equation( PDEU, { t, z }, { T_INT, ( manual & 2 ) ? Z_UB : Z_INT }, io );   //             to UB
    oc.add_equation( u - U, { t, z }, { FFDom::LB, FFDom::ALL }, ini );
    oc.add_equation( v - V, { t, z }, { FFDom::LB, FFDom::ALL }, ini );
    oc.add_equation( sg * ( u - U ) + ( v - V ), { t, z }, { T_INT, FFDom::LB }, bnd ); // incoming w+ at LB
    oc.add_equation( sg * ( u - U ) - ( v - V ), { t, z }, { T_INT, FFDom::UB }, bnd ); // incoming w- at UB
  }
  else{
    FFVar PU = OpP( u, t ) + a * OpP( u, z ) - ( Ut + a * Uz );
    FFVar PW = OpP( v, t ) + b * OpP( v, z ) - ( Vt + b * Vz );
    oc.add_equation( PU, { t, z }, { T_INT, ( manual & 1 ) ? Z_UB : Z_INT }, io );
    oc.add_equation( PW, { t, z }, { T_INT, ( manual & 2 ) ? Z_UB : Z_INT }, io );
    oc.add_equation( u - U, { t, z }, { FFDom::LB, FFDom::ALL }, ini );
    oc.add_equation( v - V, { t, z }, { FFDom::LB, FFDom::ALL }, ini );
    oc.add_equation( u - U, { t, z }, { T_INT, FFDom::LB }, bnd );                // both inflow at LB
    oc.add_equation( v - V, { t, z }, { T_INT, FFDom::LB }, bnd );
  }
  oc.options.DISPLAY_LEVEL         = 0;
  oc.options.REDUCE.ORDER          = OCFESLV::Options::RED_FULL;
  oc.options.CLASSIFY.MODE         = OCFESLV::Options::CLASS_AUTO;
  oc.options.INTERFACE.TYPE        = OCFESLV::Options::IC_AUTO;
  oc.options.INTERFACE.IMPOSITION  = OCFESLV::Options::IC_WEAK;
  oc.options.INTERFACE.SAT_SIGMA0  = 10.0;
  oc.options.AUTO.HYP_CLOSURE     = auto_on;   // an option again since 2026-10-07 (was CRONOS_AUTO_HYP_CLOSURE)
  oc.options.SOLVE.MARCHING        = false;
  oc.options.SOLVE.MAX_ITER        = 10;
  oc.options.SOLVE.RES_TOL         = 1e-11;
#if defined(CRONOS__WITH_SPQR)
  oc.options.SOLVE.FACTORIZATION   = OCFESLV::Options::SOLVE_SPQR;   // available only with SuiteSparseQR (GPL build)
#endif
  R.setup = oc.setup();
  if( !R.setup ) return R;
  R.nVar = oc.n_colloc_sta(); R.nEqn = oc.n_colloc_eqn(); R.square = ( R.nVar == R.nEqn );
  if( !R.square ) return R;
  R.var.assign( R.nVar, 0.0 );
  R.solved = oc.solve( R.var.data() ).converged;
  return R;
}
static double maxdiff( Run const& A, Run const& B )
{
  if( A.var.size() != B.var.size() || A.var.empty() ) return INFINITY;
  double m = 0; for( size_t i = 0; i < A.var.size(); ++i ) m = std::max( m, std::fabs( A.var[i] - B.var[i] ) );
  return m;
}
static std::string cnt( Run const& R ){ return "nVar=" + std::to_string( R.nVar ) + " nEqn=" + std::to_string( R.nEqn ); }

int main()
{
  std::cout << "================================================================\n"
            << "  hyperbolic closure: append only the per-face deficit   header: " << OCFE_OCFESLV_HEADER << "\n"
            << "================================================================\n";
  char const* nm[3] = { "K0 scalar advection", "K1 2x2 system (PDE15), PDE7-style closure", "K2 decoupled pair, both out at UB" };
  for( int kind = 0; kind < 3; ++kind ){
    std::cout << "\n---- " << nm[kind] << " ----\n";
    int const full = ( kind == 0 ) ? 1 : 3;
    Run A  = run( kind, true,  0    );    // automatic closure only
    Run M  = run( kind, false, full );    // manual closure only: the reference
    Run MA = run( kind, true,  full );    // manual closure AND automatic closure on
    check( "automatic closure only: square and solved", A.square && A.solved, cnt( A ) );
    check( "manual closure only (reference): square and solved", M.square && M.solved, cnt( M ) );
    check( "manual + automatic: NOT double-closed (square) and solved", MA.square && MA.solved, cnt( MA ) );
    check( "manual + automatic == manual only (the same system)", maxdiff( MA, M ) <= 1e-13, "max|d|=" + sci( maxdiff( MA, M ) ) );
    if( kind != 1 )   // K0/K2: the automatic row is +-the PDE itself, so all closures give one solution
      check( "automatic only == manual only (the closure row is the PDE)", maxdiff( A, M ) <= 1e-10, "max|d|=" + sci( maxdiff( A, M ) ) );
    if( kind == 2 ){
      Run P = run( kind, true, 1 );       // PARTIAL: only the first PDE extended; the closure adds the other
      check( "PARTIAL (1 of 2 supplied): square and solved", P.square && P.solved, cnt( P ) );
      check( "PARTIAL == manual only (the missing direction appended)", maxdiff( P, M ) <= 1e-10, "max|d|=" + sci( maxdiff( P, M ) ) );
    }
  }
  std::cout << "\n  " << g_pass << " passed, " << g_fail << " failed\n"
            << "OCFE_hypclosure: " << ( g_fail ? "FAIL" : "PASS" ) << "\n";
  return g_fail ? 1 : 0;
}
