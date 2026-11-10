////////////////////////////////////////////////////////////////////////
// DAESLV_molsat -- the ENERGY-STABLE SAT form of IC_WEAK (8 Nov 2026).  The legacy SAT gives the semi-discrete operator a
// growing mode quadratic in SAT_SIGMA0 (+9.2e4 at 10 on the rod below), damped on OCFESLV by Radau's L-stability, fatal to
// an integrator; INTERFACE.SAT_FORM = SAT_ENERGY is the summation-by-parts form -- the alternating (LDG) lifting of the
// value jump into the auxiliary's LINK row and of the flux jump into the PDE row, an upwind term for a constant-
// coefficient convective part, optionally an interior penalty -- all strengths in OCFEBASE::_sat_strength.  Under the
// method of lines IC_WEAK takes it whatever SAT_FORM says; on OCFESLV SAT_LEGACY stays the default.
//   (A) u_t = 0.1 u_zz, u(t,0) = 0, u_z(t,1) = 0, u0 = sin(pi z/2), degree 4: IC_WEAK under MoL converges at p+1 = 5,
//       below IC_STRONG at 12 elements (measured 8.8e-9 vs 5.7e-8);
//   (B) the same rod with convection, u_t + u_z = 0.1 u_zz (Peclet 10), Dirichlet inflow: stable, within the measured 3x
//       of IC_STRONG at 6 elements (the convective accuracy gap is recorded, not closed);
//   (B2) Peclet 100 on 3 elements of degree 4, to t = 20: the column settles to its steady state u = 1.  The upwind term
//       of the convective part had been generated only when SAT_PENALTY != 0 -- switched off with the penalty's
//       default of 0 -- leaving a growing mode (+0.46 in the semi-discrete spectrum): u(20,1) = -1240 without the term,
//       1.000000 with it (9 Nov 2026).  IC_STRONG is printed for comparison.
//   (C) OCFESLV under SAT_ENERGY = DAESLV under MoL on (A): one spatial scheme; and SAT_LEGACY still OCFESLV's default.
//   (D) the central flux (SAT_FLUX_CENTRAL): refused on a field with a convective part (it diverges: 5e28 at 12
//       elements, Peclet 10, with or without the upwind term), taken on pure diffusion.
////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <iomanip>
#include <cmath>
#include "ocfeslv.hpp"
#include "daeslv_idas.hpp"
using namespace mc;

static int nfail = 0;
static void check( std::string const& w, bool ok ){ std::cout << "  " << std::left << std::setw(96) << w << ( ok? "PASS": "FAIL" ) << "\n"; if( !ok ) ++nfail; }
static double const D = 0.1, TEND = 1.;

template <typename S> static void declare( S& s, FFGraph& G, FFVar& t, FFVar& z, FFVar& u, int nel, int nnz, double U, int nelt, int nodt )
{
  t = G.add_var("t"); z = G.add_var("z"); u = G.add_var("u");
  s.add_domain( t, FFDom( 0., TEND, nelt, FFDom::LGR, nodt ) ); s.add_domain( z, FFDom( 0., 1., nel, FFDom::LGL, nnz ) );
  s.set_evolution_domain( t ); s.add_state( u, {t,z} ); s.update_ref( u, 0.5 );
  FFPartial OpP; FFEval OpE; FFIntegral OpI;
  int const TI = FFDom::ALL - FFDom::LB, ZI = FFDom::ALL - FFDom::LB - FFDom::UB;
  typename S::EqnOptions io( S::EqnRole::INTERIOR, 0 ), ib( S::EqnRole::BOUNDARY, 0 ), ii( S::EqnRole::INITIAL, 0 );
  if( U == 0. ){
    s.add_equation( OpP(u,t) - D*OpP(OpP(u,z),z), {t,z}, {TI, ZI}, io );
    s.add_equation( u, {t,z}, {TI, FFDom::LB}, ib );  s.add_equation( OpP(u,z), {t,z}, {TI, FFDom::UB}, ib );
    s.add_equation( u - sin( 0.5*M_PI*z ), {t,z}, {FFDom::LB, FFDom::ALL}, ii );
  }
  else{
    s.add_equation( OpP(u,t) + U*OpP(u,z) - D*OpP(OpP(u,z),z), {t,z}, {TI, ZI}, io );
    s.add_equation( u - 1., {t,z}, {TI, FFDom::LB}, ib );  s.add_equation( OpP(u,z), {t,z}, {TI, FFDom::UB}, ib );
    s.add_equation( u - 0.5 - 0.5*(1.-z)*(1.-z), {t,z}, {FFDom::LB, FFDom::ALL}, ii );
  }
  s.add_output( OpE( u, {{t,1},{z,1}}, {{t,TEND},{z,1.}} ) );
  s.add_output( OpI( u, z ), {t}, {TEND} );
}
static bool mol( int nel, int nnz, double U, bool strong, std::vector<double>& F, int* form = nullptr, bool central = false )
{
  FFGraph G; DAESLV_IDAS S( &G ); FFVar t, z, u;  declare( S, G, t, z, u, nel, nnz, U, 2, 3 );
  S.options.INTERFACE.IMPOSITION = strong? OCFEMOL::Options::IC_STRONG: OCFEMOL::Options::IC_WEAK;   // SAT_FORM left at its default
  if( central ) S.options.INTERFACE.SAT_FLUX = OCFEMOL::Options::SAT_FLUX_CENTRAL;
  S.options.SOLVE.RTOL = S.options.SOLVE.ATOL = 1e-12;  S.options.SOLVE.VERBOSE = 0;  S.options.SOLVE.NMAX = 100000;  S.options.DISPLAY_LEVEL = 0;
  if( !S.setup() || S.solve( std::vector<FFModel::InputVal>{} ) != DAESLV_IDAS::NORMAL ) return false;
  F = S.val_function();
  if( form ) *form = (int)S.semidiscrete_model()->options.INTERFACE.SAT_FORM;
  return true;
}
// (B2): u_t + u_z = Dd u_zz on [0,1], u(t,0) = 1, u_z(t,1) = 0, to t = TE; the outlet value at TE
static bool mol_pe( int nel, double Dd, double TE, bool strong, double& uT )
{
  FFGraph G; DAESLV_IDAS S( &G ); FFPartial OpP; FFEval OpE;
  FFVar t = G.add_var("t"), z = G.add_var("z"), u = G.add_var("u");
  S.add_domain( t, FFDom( 0., TE, 2, FFDom::LGR, 3 ) ); S.add_domain( z, FFDom( 0., 1., nel, FFDom::LGL, 5 ) );
  S.set_evolution_domain( t ); S.add_state( u, {t,z} ); S.update_ref( u, 0.5 );
  int const TI = FFDom::ALL - FFDom::LB, ZI = FFDom::ALL - FFDom::LB - FFDom::UB;
  DAESLV_IDAS::EqnOptions io( DAESLV_IDAS::EqnRole::INTERIOR, 0 ), ib( DAESLV_IDAS::EqnRole::BOUNDARY, 0 ), ii( DAESLV_IDAS::EqnRole::INITIAL, 0 );
  S.add_equation( OpP(u,t) + OpP(u,z) - Dd*OpP(OpP(u,z),z), {t,z}, {TI, ZI}, io );
  S.add_equation( u - 1., {t,z}, {TI, FFDom::LB}, ib );  S.add_equation( OpP(u,z), {t,z}, {TI, FFDom::UB}, ib );
  S.add_equation( u - 0.5 - 0.5*(1.-z)*(1.-z), {t,z}, {FFDom::LB, FFDom::ALL}, ii );
  S.add_output( OpE( u, {{t,1},{z,1}}, {{t,TE},{z,1.}} ) );
  S.options.INTERFACE.IMPOSITION = strong? OCFEMOL::Options::IC_STRONG: OCFEMOL::Options::IC_WEAK;
  S.options.SOLVE.RTOL = S.options.SOLVE.ATOL = 1e-9;  S.options.SOLVE.VERBOSE = 0;  S.options.SOLVE.NMAX = 200000;  S.options.DISPLAY_LEVEL = 0;
  if( !S.setup() || S.solve( std::vector<FFModel::InputVal>{} ) != DAESLV_IDAS::NORMAL ) return false;
  uT = S.val_function()[0];
  return std::isfinite( uT );
}
static double err( std::vector<double> const& F, double a, double b ){ return std::max( std::fabs( F[0] - a )/std::fabs( a ), std::fabs( F[1] - b )/std::fabs( b ) ); }

int main()
{
  std::cout << std::scientific << std::setprecision(3);
  double const e1 = std::exp( -D*0.25*M_PI*M_PI*TEND ), eI = 2./M_PI*e1;

  std::cout << "(A) pure diffusion, degree 4: IC_WEAK under MoL (the energy form) against the closed form\n";
  { std::vector<double> F;  double ew[3] = { 0, 0, 0 }, es = 0.;  int form = -1;  bool ok = true;  int const nels[3] = { 3, 6, 12 };
    for( int k = 0; k < 3; ++k ){ ok = ok && mol( nels[k], 5, 0., false, F, k == 0? &form: nullptr );  if( ok ) ew[k] = err( F, e1, eI ); }
    ok = ok && mol( 12, 5, 0., true, F );  if( ok ) es = err( F, e1, eI );
    std::cout << "    IC_WEAK: " << ew[0] << " " << ew[1] << " " << ew[2] << "  rate " << std::fixed << std::setprecision(2) << std::log2( ew[1]/ew[2] ) << std::scientific << std::setprecision(3)
              << ";  IC_STRONG at 12: " << es << "\n";
    check( "(A) IC_WEAK under MoL sets up with the energy-stable form and solves at 3, 6, 12 elements", ok && form == (int)OCFEMOL::Options::SAT_ENERGY );
    check( "(A) converges at p+1: rate from 6 to 12 elements >= 4.5", ok && std::log2( ew[1]/ew[2] ) >= 4.5 );
    check( "(A) at 12 elements below IC_STRONG's error", ok && ew[2] < es ); }

  std::cout << "(B) convection-diffusion, Peclet 10: stable, within the measured gap of IC_STRONG\n";
  { std::vector<double> Fr, Fw, Fs;
    bool const okr = mol( 24, 7, 1., true, Fr ), okw = mol( 6, 5, 1., false, Fw ), oks = mol( 6, 5, 1., true, Fs );
    double const ew = okr && okw? err( Fw, Fr[0], Fr[1] ): NAN, es = okr && oks? err( Fs, Fr[0], Fr[1] ): NAN;
    std::cout << "    at 6 elements against 24 x LGL(7): IC_WEAK " << ew << "  IC_STRONG " << es << "\n";
    check( "(B) solves under MoL with the upwind convective term", okr && okw && oks );
    check( "(B) error within 4x of IC_STRONG's (measured 2.6x: the convective gap, recorded)", okr && okw && oks && ew < 4.*es ); }

  std::cout << "(B2) Peclet 100, 3 elements of degree 4, to t = 20: settles to u = 1\n";
  { double uw = NAN, us = NAN;
    bool const okw = mol_pe( 3, 0.01, 20., false, uw ), oks = mol_pe( 3, 0.01, 20., true, us );
    std::cout << "    u(20, 1): IC_WEAK " << ( okw? std::to_string( uw ): std::string( "FAILED" ) ) << "   IC_STRONG " << ( oks? std::to_string( us ): std::string( "FAILED" ) ) << " (for comparison)\n";
    check( "(B2) IC_WEAK under MoL stays bounded and settles: |u(20,1) - 1| < 1e-3 (the upwind term active)", okw && std::fabs( uw - 1. ) < 1e-3 ); }

  std::cout << "(C) OCFESLV: SAT_ENERGY = DAESLV under MoL on (A); SAT_LEGACY still the default\n";
  { std::vector<double> Fm;  bool const okm = mol( 6, 5, 0., false, Fm );
    FFGraph G; OCFESLV O( &G ); FFVar t, z, u;  declare( O, G, t, z, u, 6, 5, 0., 4, 8 );
    bool const legacyDefault = O.options.INTERFACE.SAT_FORM == OCFESLV::Options::SAT_LEGACY;
    O.options.INTERFACE.IMPOSITION = OCFESLV::Options::IC_WEAK;  O.options.INTERFACE.SAT_FORM = OCFESLV::Options::SAT_ENERGY;
    O.options.SOLVE.MARCHING = false;  O.options.SOLVE.RES_TOL = 1e-12;  O.options.DISPLAY_LEVEL = 0;
    bool const oko = O.setup() && O.solve( std::vector<FFModel::InputVal>{} ) == OCFESLV::NORMAL;
    double const gap = okm && oko? err( O.val_function(), Fm[0], Fm[1] ): NAN;
    std::cout << "    OCFESLV vs DAESLV, 6 elements: " << gap << "\n";
    check( "(C) SAT_LEGACY is OCFESLV's default", legacyDefault );
    check( "(C) OCFESLV under SAT_ENERGY = DAESLV under MoL to 1e-8 (one spatial scheme, time errors apart)", okm && oko && gap < 1e-8 ); }

  std::cout << "(D) the central flux: refused with convection, taken on pure diffusion\n";
  { std::vector<double> F;
    bool const conv = mol( 3, 5, 1., false, F, nullptr, true ), diff = mol( 3, 5, 0., false, F, nullptr, true );
    check( "(D) SAT_FLUX_CENTRAL refused at setup on the Peclet 10 rod (a convective part)", !conv );
    check( "(D) SAT_FLUX_CENTRAL taken on the pure-diffusion rod", diff ); }

  std::cout << "==========  DAESLV_molsat: " << ( nfail? std::to_string( nfail ) + " FAILED": "ALL PASS" ) << "  ==========\n";
  return nfail? 1: 0;
}
