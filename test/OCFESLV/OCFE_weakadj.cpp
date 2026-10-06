// OCFE_weakadj.cpp -- GATE for the marching window transfer under weak imposition (WORKPLAN 1.2/1.3, 2026-10-06).
// terminal_profile() used to evaluate each spatial node at its COORDINATE, so at a spatial element boundary both copies
// of the node read one element's terminal.  Under IC_STRONG/IC_TRACE the copies are equal and nothing showed; under
// IC_WEAK they are distinct unknowns: marching differed from monolithic (whose causal continuity rows pin each copy to
// its own element), and the marching ADJOINT -- the transpose of a per-element map -- disagreed with the forward
// sensitivities and with finite differences (1.9e-5 here).  Now each copy reads its own element's terminal.
// The same flaw sat in transfer_terminal(), which writes the terminal into an IC INPUT (the canonical path) by sampling at
// coordinates -- there even the forward gradient missed its own finite differences.  Both paths are checked: the
// initial condition as a constant (override path) and through an input (canonical path), on a nonlinear heated rod
// with 4 spatial elements, for the three impositions: forward = adjoint, marching = monolithic (values and gradients),
// and the marching gradient = central finite differences.
#include <cstdio>
#include <cmath>
#include "ocfeslv.hpp"
using namespace mc;  typedef OCFESLV::EqnRole Role;  typedef OCFESLV::EqnOptions EO;
static int nfail = 0;
static void check( char const* what, bool ok ){ std::printf( "  %-70s %s\n", what, ok? "PASS": "FAIL" ); nfail += !ok; }

struct Res { double f = NAN; std::vector<double> Jf, Ja, fd; };
static Res run( OCFESLV::Options::ImpositionType imp, bool march, bool with_fd, bool ic_input ){
  FFPartial OpP;  FFEval OpE;  int const NLB = FFDom::ALL - FFDom::LB, INN = FFDom::ALL - FFDom::LB - FFDom::UB;
  FFGraph G;  FFVar t = G.add_var( "t" ), z = G.add_var( "z" ), u = G.add_var( "u" ), q = G.add_var( "q" ), k = G.add_var( "k" ), u0 = G.add_var( "u0" );
  OCFESLV S( &G );
  S.add_domain( t, FFDom( 0., 1., 4, FFDom::LGR, 3 ) );  S.add_domain( z, FFDom( 0., 1., 4, FFDom::LGL, 5 ) );
  S.set_evolution_domain( t );  S.add_state( u, {t, z} );
  S.add_input( q );  S.update_ref( q, 1. );  S.add_input( k );  S.update_ref( k, 0.1 );
  if( ic_input ){ S.add_input( u0, {z} );  S.update_ref( u0, 0. ); }
  S.add_equation( OpP( u, t ) - k * ( 1. + 0.5 * u ) * OpP( u, {{z,2}} ), {t, z}, {NLB, INN}, EO( Role::INTERIOR ) );
  S.add_equation( u, {t, z}, {NLB, (int)FFDom::LB}, EO( Role::BOUNDARY ) );
  S.add_equation( k * OpP( u, z ) - q, {t, z}, {NLB, (int)FFDom::UB}, EO( Role::BOUNDARY ) );
  S.add_equation( ic_input? u - u0: u, {t, z}, {(int)FFDom::LB, (int)FFDom::ALL}, EO( Role::INITIAL ) );
  S.add_output( OpE( u, {{t,1},{z,1}}, {{t,1.},{z,1.}} ) );
  S.register_control( q );  S.register_control( k );
  S.options.INTERFACE.IMPOSITION = imp;  S.options.SOLVE.MARCHING = march;
  S.options.SOLVE.RES_TOL = 1e-12;  S.options.DISPLAY_LEVEL = 0;
  Res R;  if( !S.setup() ) return R;
  std::vector<double> var, inp;  S.init( var, inp );  std::vector<double> const v0 = var, i0 = inp;
  S.solve_fsens( var.data(), inp.data(), nullptr );  R.Jf = S.sens_jacobian();  R.f = S.sens_functions()[0];
  var = v0; inp = i0;  S.solve_asens( var.data(), inp.data(), nullptr );  R.Ja = S.sens_jacobian();
  if( with_fd ){                                    // inputs: q, k first (in declaration order)
    double const h = 1e-6;
    for( size_t c = 0; c < 2; ++c ){
      std::vector<double> vp = v0, ip = i0, vm = v0, im = i0;  ip[c] += h;  im[c] -= h;
      S.solve( vp.data(), ip.data(), nullptr );  double const fp = S.val_functions()[0];
      S.solve( vm.data(), im.data(), nullptr );  double const fm = S.val_functions()[0];
      R.fd.push_back( ( fp - fm ) / ( 2. * h ) );
    }
  }
  return R;
}
static double maxdiff( std::vector<double> const& a, std::vector<double> const& b ){
  if( a.size() != b.size() || a.empty() ) return INFINITY;
  double d = 0.;  for( size_t i = 0; i < a.size(); ++i ) d = std::max( d, std::fabs( a[i] - b[i] ) );  return d;
}

int main(){
  std::printf( "OCFE_weakadj -- marching window transfer: forward, adjoint, monolithic, finite differences\n" );
  char const* NAME[] = { "IC_STRONG", "IC_WEAK", "IC_TRACE" };
  OCFESLV::Options::ImpositionType const IMP[] = { OCFESLV::Options::IC_STRONG, OCFESLV::Options::IC_WEAK, OCFESLV::Options::IC_TRACE };
  char buf[160];
  for( int path = 0; path < 2; ++path ) for( int i = 0; i < 3; ++i ){
    bool const ic_input = ( path == 1 );
    if( i == 0 ) std::printf( "  -- initial condition %s\n", ic_input? "through an INPUT (canonical path: transfer_terminal)": "as a constant (override path: terminal_profile)" );
    Res M = run( IMP[i], false, false, ic_input ), K = run( IMP[i], true, true, ic_input );
    double const dfa = maxdiff( K.Jf, K.Ja ), dmk = maxdiff( K.Jf, M.Jf ), dfd = maxdiff( K.Jf, K.fd ), dv = std::fabs( K.f - M.f );
    std::printf( "  %-9s | du/dk marching fwd %+.10f adj %+.10f fd %+.10f | monolithic %+.10f\n",
                 NAME[i], K.Jf.size() > 1? K.Jf[1]: NAN, K.Ja.size() > 1? K.Ja[1]: NAN, K.fd.size() > 1? K.fd[1]: NAN, M.Jf.size() > 1? M.Jf[1]: NAN );
    std::snprintf( buf, sizeof buf, "%s%s: marching forward = adjoint (1e-10; %.1e)", NAME[i], ic_input? " [input IC]": "", dfa );        check( buf, dfa < 1e-10 );
    std::snprintf( buf, sizeof buf, "%s%s: marching = monolithic, value and gradient (1e-9; %.1e, %.1e)", NAME[i], ic_input? " [input IC]": "", dv, dmk ); check( buf, dv < 1e-9 && dmk < 1e-9 );
    std::snprintf( buf, sizeof buf, "%s%s: marching gradient = finite differences (1e-6; %.1e)", NAME[i], ic_input? " [input IC]": "", dfd );              check( buf, dfd < 1e-6 );
  }
  std::printf( "  OCFE_weakadj: %d failed\n", nfail );
  return nfail? 1: 0;
}
