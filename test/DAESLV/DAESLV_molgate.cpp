// DAESLV_molgate -- MoL phase 5: the gate's third leg.  The heated rod, discretised in z by OCFEMOL, exported to
// DAESLV_IDAS through export_semidiscrete() and integrated to T, against
//   leg 1  OCFESLV, collocation in BOTH directions -- the same spatial plan, so under IMPOSITION.NEUMANN_ALGEBRAIC
//          legs 1 and 3 are the SAME discrete system in space and differ only by the time discretisation;
//   leg 2  the hand-expanded collocated rows of DAESLV_mol, integrated by the same IDAS.
// Legs 2 and 3 share the spatial discretisation and the integrator, so they must agree to the integration
// tolerance -- that is the exactness claim of the plan, and the first check here.  Leg 1 is compared on a time
// ladder: refining ITS time collocation must drive |1-3| down to the integration tolerance.
// Report-only timings (CRONOS_SNAPSHOT suppresses them, as in DAESLV_molbench).
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <chrono>
#include "daeslv_idas.hpp"
#include "ffocfemol.hpp"
#include "ocfeslv.hpp"

using namespace mc;
static int n_pass = 0, n_fail = 0;
static void check( std::string const& w, bool ok )
{ ok? ++n_pass: ++n_fail; std::cout << "    " << std::left << std::setw(62) << w << ( ok? "PASS": "FAIL" ) << "\n"; }
static void check_lt( std::string const& w, double got, double bound )
{
  bool const ok = ( got == got ) && got < bound;
  ok? ++n_pass: ++n_fail;
  std::cout << "    " << std::left << std::setw(50) << w << std::right << std::scientific << std::setprecision(2)
            << std::setw(10) << got << " < " << std::setw(8) << bound << "  " << ( ok? "PASS": "FAIL" ) << std::defaultfloat << "\n";
}

static constexpr double ALPHA = 0.1, QFLUX = 1.0, TEND = 1.0;

struct Mesh
{
  FFDom dom; size_t nel, nnod; std::vector<double> z;
  Mesh( size_t nel_, size_t nnod_, double lo = 0., double hi = 1. ) : dom( lo, hi, nel_, FFDom::LGL, nnod_ ), nel( nel_ ), nnod( nnod_ )
  {
    dom.set_nodes();
    for( size_t e = 0; e < nel; ++e ){ auto ze = dom.lgnodes( dom.elem_lo(e), dom.elem_up(e) ); for( size_t i = 0; i < nnod; ++i ) z.push_back( ze[i] ); }
  }
  size_t n() const { return nel*nnod; }
  size_t g( size_t e, size_t i ) const { return e*nnod + i; }
};

// The closed form: u = (q/a) z - (2q/a) sum_n (-1)^(n-1) lam_n^-2 exp(-a lam_n^2 t) sin(lam_n z), lam_n = (2n-1) pi/2.
// u(1,1) = 3.5682340045.  The one reference here that is not another discretisation.
static double u_exact( double t, double z, int NTERM = 400 )
{
  double s = 0.;
  for( int n = NTERM; n >= 1; --n ){
    double const lam = ( 2.*n - 1. ) * M_PI / 2.;
    s += ( ( n % 2 )? 1.: -1. ) / ( lam*lam ) * std::exp( -ALPHA*lam*lam*t ) * std::sin( lam*z );
  }
  return QFLUX/ALPHA*z - 2.*QFLUX/ALPHA*s;
}

static double maxdiff( std::vector<double> const& a, std::vector<double> const& b )
{
  if( a.size() != b.size() || a.empty() ) return 1e30;
  double d = 0.; for( size_t i = 0; i < a.size(); ++i ) d = std::max( d, std::fabs( a[i] - b[i] ) );
  return d;
}

// ---------------------------------------------------------------- the model, declared on any OCFE solver
template <typename OC>
static void declare_rod( OC& oc, FFGraph& G, FFVar& t, FFVar& z, FFVar& u, Mesh const& M, size_t NELT, size_t NODT )
{
  t = G.add_var("t"); z = G.add_var("z"); u = G.add_var("u");
  oc.add_domain( t, FFDom( 0., TEND, NELT, FFDom::LGR, NODT ) );
  oc.add_domain( z, FFDom( 0., 1., M.nel, FFDom::LGL, M.nnod ) );
  oc.set_evolution_domain( t );
  oc.add_state( u, {t,z} ); oc.update_ref( u, 0. );
  FFPartial OpP;
  int const TI = FFDom::ALL - FFDom::LB, ZI = FFDom::ALL - FFDom::LB - FFDom::UB;
  typename OC::EqnOptions io( OC::EqnRole::INTERIOR, 0 ), ib( OC::EqnRole::BOUNDARY, 0 ), ii( OC::EqnRole::INITIAL, 0 );
  oc.add_equation( OpP(u,t) - OpP( ALPHA*OpP(u,z), z ), {t,z}, {TI, ZI},                io );
  oc.add_equation( u,                                   {t,z}, {TI, FFDom::LB},         ib );
  oc.add_equation( ALPHA*OpP(u,z) - QFLUX,              {t,z}, {TI, FFDom::UB},         ib );
  oc.add_equation( u,                                   {t,z}, {FFDom::LB, FFDom::ALL}, ii );
}

// ---------------------------------------------------------------- leg 1: OCFESLV, collocation in t and z
static std::vector<double> leg1( Mesh const& M, size_t NELT, size_t NODT )
{
  FFGraph G; OCFESLV oc( &G ); FFVar t, z, u;
  declare_rod( oc, G, t, z, u, M, NELT, NODT );
  oc.options.INTERFACE.IMPOSITION = OCFESLV::Options::IC_STRONG;
  oc.options.SOLVE.MARCHING = false;  oc.options.SOLVE.RES_TOL = 1e-12;  oc.options.DISPLAY_LEVEL = 0;
  if( !oc.setup() ) return {};
  std::vector<double> xv, inp;
  if( !oc.init( xv, inp, nullptr ) ) return {};
  auto const rep = oc.solve( xv.data(), inp.empty()? nullptr: inp.data(), nullptr );
  if( !rep.converged ) return {};
  std::vector<double> out( M.n() );
  for( size_t k = 0; k < M.n(); ++k ){ std::map<FFVar,double,lt_FFVar> pt{ { t, TEND }, { z, M.z[k] } }; out[k] = oc.eval_solution( u, pt ); }
  return out;
}

// ---------------------------------------------------------------- leg 3: OCFEMOL exported to DAESLV
struct Leg3 { std::vector<double> unode; size_t n_states = 0, n_diff = 0, n_rows = 0; int index = -1; bool ok = false; double setup_us = 0., solve_us = 0.; };
static Leg3 leg3( Mesh const& M, double rtol, double atol, bool verbose, bool neumann_diff = false )
{
  Leg3 R;
  auto const t0 = std::chrono::steady_clock::now();
  FFGraph G; OCFEMOL mol( &G ); FFVar t, z, u;
  declare_rod( mol, G, t, z, u, M, 1, 2 );                 // the t domain here is only the span: MoL does not collocate it
  mol.options.INTERFACE.IMPOSITION = OCFEMOL::Options::IC_STRONG;
  mol.options.IMPOSITION.NEUMANN_NODE = neumann_diff? OCFEMOL::Options::NEUMANN_DIFFERENTIAL: OCFEMOL::Options::NEUMANN_ALGEBRAIC;
  mol.options.DISPLAY_LEVEL = verbose? 1: 0;
  if( !mol.setup() ){ std::cout << "      leg3: OCFEMOL::setup failed\n"; return R; }
  R.n_states = mol.semi_discrete().n_states; R.n_diff = mol.semi_discrete().n_diff_cols; R.n_rows = mol.semi_discrete().n_rows;

  FFGraph H; DAESLV_IDAS S; std::vector<FFVar> xc; FFVar tt; std::string err;
  if( !export_semidiscrete( mol, H, S, FFDom( 0., TEND, 2, FFDom::LGR, 3 ), xc, tt, &err ) ){
    std::cout << "      leg3: export failed: " << err << "\n"; return R;
  }
  // the start: OCFEMOL seeds (references, the Ini stream); DAESLV's sparse KINSOL corrector solves for consistency
  std::vector<double> x0( mol.semi_discrete().n_states );
  if( !mol.init_seed( 0., x0.data() ) ){ std::cout << "      leg3: init_seed failed\n"; return R; }
  for( size_t i = 0; i < x0.size(); ++i ) S.update_ref( xc[i], x0[i] );

  FFEval OpE;
  for( size_t k = 0; k < xc.size(); ++k ) S.add_output( OpE( xc[k], tt, TEND ) );
  S.options.SOLVE.RTOL = rtol;  S.options.SOLVE.ATOL = atol;  S.options.SOLVE.VERBOSE = 0;  S.options.DISPLAY_LEVEL = 0;
  if( !S.setup() ){ std::cout << "      leg3: DAESLV::setup failed: " << S.extract_error() << "\n"; return R; }
  R.index = S.structural_index( 0, tt ).index;
  std::cout << "      leg3: exported DAE -- " << mol.semi_discrete().n_states << " unknowns ("
            << mol.semi_discrete().n_tau << " explicit tau), " << mol.semi_discrete().n_rows << " rows, "
            << mol.semi_discrete().n_diff_cols << " differential; structural index along t = " << R.index << "\n";
  auto const t1 = std::chrono::steady_clock::now();
  std::vector<double> P( S.np(), 0. );
  if( S.solve( P ) != DAESLV_IDAS::NORMAL ){ std::cout << "      leg3: solve failed: " << S.extract_error() << "\n"; R.index = R.index; return R; }
  auto const t2 = std::chrono::steady_clock::now();
  auto const all = S.val_function();
  // the u-field values at the mesh nodes, in mesh order: the columns of the field u
  auto const& ucols = mol.field_columns( u );          // u's columns in mesh order
  R.unode.resize( ucols.size() );
  for( size_t k = 0; k < ucols.size(); ++k ) R.unode[k] = all[ ucols[k] ];
  R.ok = ( ucols.size() == M.n() );
  if( !R.ok ) std::cout << "      leg3: u has " << ucols.size() << " columns, the mesh has " << M.n() << " nodes\n";
  R.setup_us = std::chrono::duration<double,std::micro>( t1 - t0 ).count();
  R.solve_us = std::chrono::duration<double,std::micro>( t2 - t1 ).count();
  return R;
}


// ---------------------------------------------------------------- leg 2: the hand-expanded collocated rows
// The SAME spatial discretisation written out node by node as ordinary DAG rows and integrated by the same IDAS.
// Leg 2 treats the initial point exactly as the exported DAE must: the boundary rows hold at t0 too (a DAE's
// algebraic constraints cannot be switched off at the initial point), and only the differential nodes carry an
// initial condition.  So legs 2 and 3 must agree to the integration tolerance -- that is the export's own check.
static std::vector<double> leg2( Mesh const& M, double rtol, double atol )
{
  FFGraph G; DAESLV_IDAS S; S.set_dag( &G );
  FFPartial OpP; FFEval OpE;
  FFVar t = G.add_var("t");
  size_t const NG = M.n();
  std::vector<FFVar> u( NG ), D( NG );
  for( size_t k = 0; k < NG; ++k ){ u[k] = G.add_var(); D[k] = G.add_var(); }
  S.add_domain( t, FFDom( 0., TEND, 2, FFDom::LGR, 3 ) );
  S.set_evolution_domain( t );
  for( size_t k = 0; k < NG; ++k ){ S.add_state( u[k], {t} ); S.update_ref( u[k], 0. ); S.add_state( D[k], {t} ); S.update_ref( D[k], QFLUX/ALPHA ); }
  auto ddz = [&]( std::vector<FFVar> const& f ){
    std::vector<FFVar> out( NG );
    for( size_t e = 0; e < M.nel; ++e ){
      std::vector<FFVar> fe, dfe;
      for( size_t i = 0; i < M.nnod; ++i ) fe.push_back( f[ e*M.nnod + i ] );
      M.dom.eval_lagrange( dfe, 1, fe, 2./M.dom.elem_width(e) );
      for( size_t i = 0; i < M.nnod; ++i ) out[ e*M.nnod + i ] = dfe[i];
    }
    return out;
  };
  std::vector<FFVar> const dudz = ddz( u ), dDdz = ddz( D );
  auto const INTER = FFModel::EqnRole::INTERIOR, INIT = FFModel::EqnRole::INITIAL;
  int const TI = FFDom::ALL - FFDom::LB;
  for( size_t e = 0; e < M.nel; ++e )
  for( size_t i = 0; i < M.nnod; ++i ){
    size_t const k = e*M.nnod + i;
    bool const dup = ( i == 0 && e > 0 ), z0 = ( e == 0 && i == 0 ), z1 = ( e == M.nel-1 && i == M.nnod-1 );
    if( z0 )       S.add_equation( u[k],                                  {t}, {FFDom::ALL}, INTER );
    else if( dup ) S.add_equation( u[k] - u[ (e-1)*M.nnod + M.nnod-1 ],   {t}, {FFDom::ALL}, INTER );
    else         { S.add_equation( OpP(u[k],t) - ALPHA*dDdz[k],           {t}, {TI},         INTER );
                   S.add_equation( u[k],                                  {t}, {FFDom::LB},  INIT  ); }
    bool const ltrace = ( i == M.nnod-1 && e+1 < M.nel );
    if( z1 )          S.add_equation( ALPHA*D[k] - QFLUX,                 {t}, {FFDom::ALL}, INTER );
    else if( ltrace ) S.add_equation( D[k] - D[ (e+1)*M.nnod ],           {t}, {FFDom::ALL}, INTER );
    else              S.add_equation( dudz[k] - D[k],                     {t}, {FFDom::ALL}, INTER );
  }
  for( size_t k = 0; k < NG; ++k ) S.add_output( OpE( u[k], t, TEND ) );
  S.options.SOLVE.RTOL = rtol;  S.options.SOLVE.ATOL = atol;  S.options.SOLVE.VERBOSE = 0;  S.options.DISPLAY_LEVEL = 0;
  if( !S.setup() ){ std::cout << "      leg2 setup failed: " << S.extract_error() << "\n"; return {}; }
  std::vector<double> P( S.np(), 0. );
  if( S.solve( P ) != DAESLV_IDAS::NORMAL ){ std::cout << "      leg2 solve failed: " << S.extract_error() << "\n"; return {}; }
  return S.val_function();
}


// ======================================================================== case B: Belousov-Zhabotinskii, flux form
// Two fields, no-flux at both ends, tanh fronts of width 2 (phase 0's case B, verbatim).  Nonlinear, so the
// corrector's line search does real work at t0; no closed form, so the references are leg 1 (the same spatial
// system: legs 1 and 3 must converge onto each other as leg 1's time error is refined away) and leg 2 (a different
// spatial discretisation: spectral convergence between 2 and 3 under spatial refinement).
double const BZ_ZL = -30., BZ_ZR = 20., BZ_R = 10., BZ_B = 1.25, BZ_T = 1., BZ_W = 2.;
static double bz_f1( double z ){ return 0.5*( 1. + std::tanh( z/BZ_W ) ); }
static double bz_f2( double z ){ return 0.5*( 1. - std::tanh( z/BZ_W ) ); }

template <typename OC>
static void declare_bz( OC& oc, FFGraph& G, FFVar& t, FFVar& z, FFVar& u1, FFVar& u2, Mesh const& M, size_t NELT, size_t NODT )
{
  FFPartial OpP;
  t = G.add_var("t"); z = G.add_var("z"); u1 = G.add_var("u1"); u2 = G.add_var("u2");
  oc.add_domain( t, FFDom( 0., BZ_T, NELT, FFDom::LGR, NODT ) );
  oc.add_domain( z, FFDom( BZ_ZL, BZ_ZR, M.nel, FFDom::LGL, M.nnod ) );
  oc.set_evolution_domain( t );
  oc.add_state( u1, {t,z} ); oc.update_ref( u1, 0.5 );  oc.add_state( u2, {t,z} ); oc.update_ref( u2, 0.5 );
  int const TI = FFDom::ALL - FFDom::LB, ZI = FFDom::ALL - FFDom::LB - FFDom::UB;
  typename OC::EqnOptions io( OC::EqnRole::INTERIOR, 0 ), ib( OC::EqnRole::BOUNDARY, 0 ), ii( OC::EqnRole::INITIAL, 0 );
  oc.add_equation( OpP(u1,t) - OpP(OpP(u1,z),z) - u1*(1.-u1-BZ_R*u2), {t,z}, {TI, ZI}, io );
  oc.add_equation( OpP(u2,t) - OpP(OpP(u2,z),z) + BZ_B*u1*u2,         {t,z}, {TI, ZI}, io );
  oc.add_equation( OpP(u1,z), {t,z}, {TI, FFDom::LB}, ib );  oc.add_equation( OpP(u1,z), {t,z}, {TI, FFDom::UB}, ib );
  oc.add_equation( OpP(u2,z), {t,z}, {TI, FFDom::LB}, ib );  oc.add_equation( OpP(u2,z), {t,z}, {TI, FFDom::UB}, ib );
  FFVar const zz = z;
  oc.add_equation( u1 - 0.5*( 1. + tanh( zz/BZ_W ) ), {t,z}, {FFDom::LB, FFDom::ALL}, ii );
  oc.add_equation( u2 - 0.5*( 1. - tanh( zz/BZ_W ) ), {t,z}, {FFDom::LB, FFDom::ALL}, ii );
}

// both fields at T, u1 then u2, in mesh order
static std::vector<double> bz_leg1( Mesh const& M, size_t NELT, size_t NODT )
{
  FFGraph G; OCFESLV oc( &G ); FFVar t, z, u1, u2;
  declare_bz( oc, G, t, z, u1, u2, M, NELT, NODT );
  oc.options.SOLVE.MARCHING = false;  oc.options.SOLVE.RES_TOL = 1e-12;  oc.options.DISPLAY_LEVEL = 0;
  if( !oc.setup() ) return {};
  std::vector<double> xv, inp;
  if( !oc.init( xv, inp, nullptr ) ) return {};
  if( !oc.solve( xv.data(), inp.empty()? nullptr: inp.data(), nullptr ).converged ) return {};
  std::vector<double> out( 2*M.n() );
  for( size_t k = 0; k < M.n(); ++k ){
    std::map<FFVar,double,lt_FFVar> pt{ { t, BZ_T }, { z, M.z[k] } };
    out[k] = oc.eval_solution( u1, pt );  out[M.n()+k] = oc.eval_solution( u2, pt );
  }
  return out;
}

static std::vector<double> bz_leg2( Mesh const& M, double rtol, double atol )   // DAESLV_mol's bz_hand, verbatim in substance
{
  FFGraph G; DAESLV_IDAS S; S.set_dag( &G ); FFPartial OpP; FFEval OpE;
  FFVar t = G.add_var("t");
  size_t const NG = M.n();
  std::vector<FFVar> u1(NG), u2(NG), D1(NG), D2(NG);
  auto is_end = [&]( size_t k ){ return k == 0 || k == NG-1; };
  for( size_t k = 0; k < NG; ++k ){
    u1[k]=G.add_var(); u2[k]=G.add_var();
    if( is_end(k) ){ D1[k] = FFVar( 0. ); D2[k] = FFVar( 0. ); } else { D1[k]=G.add_var(); D2[k]=G.add_var(); }
  }
  S.add_domain( t, FFDom( 0., BZ_T, 2, FFDom::LGR, 3 ) );  S.set_evolution_domain( t );
  for( size_t k = 0; k < NG; ++k ){
    S.add_state( u1[k], {t} ); S.update_ref( u1[k], 0.5 );  S.add_state( u2[k], {t} ); S.update_ref( u2[k], 0.5 );
    if( is_end(k) ) continue;
    S.add_state( D1[k], {t} ); S.update_ref( D1[k], 0. );   S.add_state( D2[k], {t} ); S.update_ref( D2[k], 0. );
  }
  auto ddz = [&]( std::vector<FFVar> const& f ){
    std::vector<FFVar> out( NG );
    for( size_t e = 0; e < M.nel; ++e ){
      std::vector<FFVar> fe, dfe;
      for( size_t i = 0; i < M.nnod; ++i ) fe.push_back( f[M.g(e,i)] );
      M.dom.eval_lagrange( dfe, 1, fe, 2./M.dom.elem_width(e) );
      for( size_t i = 0; i < M.nnod; ++i ) out[M.g(e,i)] = dfe[i];
    }
    return out;
  };
  auto const du1 = ddz(u1), du2 = ddz(u2), dD1 = ddz(D1), dD2 = ddz(D2);
  int const TI = FFDom::ALL - FFDom::LB;
  auto const INTER = FFModel::EqnRole::INTERIOR, INIT = FFModel::EqnRole::INITIAL;
  for( size_t e = 0; e < M.nel; ++e )
  for( size_t i = 0; i < M.nnod; ++i ){
    size_t const k = M.g(e,i);
    bool const dup = ( i == 0 && e > 0 ), ltrace = ( i == M.nnod-1 && e+1 < M.nel );
    bool const end = ( e == 0 && i == 0 ) || ( e == M.nel-1 && i == M.nnod-1 );
    if( dup ){
      S.add_equation( u1[k] - u1[M.g(e-1,M.nnod-1)], {t}, {FFDom::ALL}, INTER );
      S.add_equation( u2[k] - u2[M.g(e-1,M.nnod-1)], {t}, {FFDom::ALL}, INTER );
    }
    else{
      S.add_equation( OpP(u1[k],t) - dD1[k] - u1[k]*(1.-u1[k]-BZ_R*u2[k]), {t}, {TI}, INTER );
      S.add_equation( OpP(u2[k],t) - dD2[k] + BZ_B*u1[k]*u2[k],           {t}, {TI}, INTER );
      S.add_equation( u1[k] - bz_f1(M.z[k]), {t}, {FFDom::LB}, INIT );
      S.add_equation( u2[k] - bz_f2(M.z[k]), {t}, {FFDom::LB}, INIT );
    }
    if( end ){ }
    else if( ltrace ){ S.add_equation( D1[k] - D1[M.g(e+1,0)], {t}, {FFDom::ALL}, INTER );
                       S.add_equation( D2[k] - D2[M.g(e+1,0)], {t}, {FFDom::ALL}, INTER ); }
    else{ S.add_equation( du1[k] - D1[k], {t}, {FFDom::ALL}, INTER );  S.add_equation( du2[k] - D2[k], {t}, {FFDom::ALL}, INTER ); }
  }
  for( size_t k = 0; k < NG; ++k ) S.add_output( OpE( u1[k], t, BZ_T ) );
  for( size_t k = 0; k < NG; ++k ) S.add_output( OpE( u2[k], t, BZ_T ) );
  S.options.SOLVE.RTOL = rtol;  S.options.SOLVE.ATOL = atol;  S.options.SOLVE.VERBOSE = 0;  S.options.DISPLAY_LEVEL = 0;
  if( !S.setup() ) return {};
  std::vector<double> P( S.np(), 0. );
  if( S.solve( P ) != DAESLV_IDAS::NORMAL ){ std::cout << "      bz leg2 solve failed: " << S.extract_error() << "\n"; return {}; }
  return S.val_function();
}

static std::vector<double> bz_leg3( Mesh const& M, double rtol, double atol )
{
  FFGraph G; OCFEMOL mol( &G ); FFVar t, z, u1, u2;
  declare_bz( mol, G, t, z, u1, u2, M, 1, 2 );
  mol.options.DISPLAY_LEVEL = 0;
  if( !mol.setup() ){ std::cout << "      bz leg3: OCFEMOL::setup failed\n"; return {}; }
  FFGraph H; DAESLV_IDAS S; std::vector<FFVar> xc; FFVar tt; std::string err;
  if( !export_semidiscrete( mol, H, S, FFDom( 0., BZ_T, 2, FFDom::LGR, 3 ), xc, tt, &err ) ){ std::cout << "      bz leg3: " << err << "\n"; return {}; }
  std::vector<double> x0( mol.semi_discrete().n_states );
  if( !mol.init_seed( 0., x0.data() ) ) return {};
  for( size_t i = 0; i < x0.size(); ++i ) S.update_ref( xc[i], x0[i] );
  FFEval OpE;
  for( size_t k = 0; k < xc.size(); ++k ) S.add_output( OpE( xc[k], tt, BZ_T ) );
  S.options.SOLVE.RTOL = rtol;  S.options.SOLVE.ATOL = atol;  S.options.SOLVE.VERBOSE = 0;  S.options.DISPLAY_LEVEL = 0;
  if( !S.setup() ){ std::cout << "      bz leg3: DAESLV::setup failed: " << S.extract_error() << "\n"; return {}; }
  std::vector<double> P( S.np(), 0. );
  if( S.solve( P ) != DAESLV_IDAS::NORMAL ){ std::cout << "      bz leg3: solve failed: " << S.extract_error() << "\n"; return {}; }
  auto const all = S.val_function();
  auto const& c1 = mol.field_columns( u1 ); auto const& c2 = mol.field_columns( u2 );
  if( c1.size() != M.n() || c2.size() != M.n() ) return {};
  std::vector<double> out( 2*M.n() );
  for( size_t k = 0; k < M.n(); ++k ){ out[k] = all[ c1[k] ]; out[M.n()+k] = all[ c2[k] ]; }
  return out;
}


// ======================================================================== case K: PDE5, the 2-D manufactured heat equation
// u_t = k (u_xx + u_yy) + f on the unit square, Dirichlet on all four sides, a Gaussian bump spreading in time; the exact
// solution is known.  Two directions exercise the export's column map, the face resolution and the fold in 2-D.
namespace pde5 {
  double const KAPPA = 0.1, X0 = 0.3, Y0 = 0.6, W0 = 0.05, SPREAD = 0.02, A = 1., BG = 0.2, T = 1.;
  double exact( double t, double x, double y ){
    double const S = W0 + SPREAD*t, RR = (x-X0)*(x-X0) + (y-Y0)*(y-Y0);
    return BG + A*W0*std::exp( -RR/S )/S;
  }
  template <typename OC>
  void declare( OC& oc, FFGraph& G, FFVar& t, FFVar& x, FFVar& y, FFVar& U, size_t NE, size_t NN, size_t NELT, size_t NODT )
  {
    FFPartial OpP;
    t = G.add_var("t"); x = G.add_var("x"); y = G.add_var("y"); U = G.add_var("U");
    FFVar DX = x - X0, DY = y - Y0, S = W0 + SPREAD*t, RR = DX*DX + DY*DY;
    FFVar BE = A*W0*exp( -RR/S )/S, UE = BG + BE;
    FFVar FE = BE*SPREAD*( -1./S + RR/(S*S) ) - KAPPA*BE*( 4.*RR/(S*S) - 4./S );
    FFVar PDE = OpP(U,t) - KAPPA*( OpP(U,{x,2}) + OpP(U,{y,2}) ) - FE, BC = U - UE;
    oc.add_domain( t, FFDom( 0., T, NELT, FFDom::LGR, NODT ) );
    oc.add_domain( x, FFDom( 0., 1., NE, FFDom::CGL, NN ) );
    oc.add_domain( y, FFDom( 0., 1., NE, FFDom::CGL, NN ) );
    oc.set_evolution_domain( t );
    oc.add_state( U, {t,x,y} ); oc.update_ref( U, BG );
    int const TI = FFDom::ALL - FFDom::LB, XI = FFDom::ALL - FFDom::LB - FFDom::UB, YI = XI;
    typename OC::EqnOptions io( OC::EqnRole::INTERIOR, 0 ), ib( OC::EqnRole::BOUNDARY, 0 ), ii( OC::EqnRole::INITIAL, 0 );
    oc.add_equation( PDE, {t,x,y}, {TI, XI, YI}, io );
    oc.add_equation( BC,  {t,x,y}, {FFDom::LB, FFDom::ALL, FFDom::ALL}, ii );
    oc.add_equation( BC,  {t,x,y}, {TI, FFDom::LB, YI}, ib );  oc.add_equation( BC, {t,x,y}, {TI, FFDom::UB, YI}, ib );
    oc.add_equation( BC,  {t,x,y}, {TI, FFDom::ALL, FFDom::LB}, ib );  oc.add_equation( BC, {t,x,y}, {TI, FFDom::ALL, FFDom::UB}, ib );
  }
  // leg 3: U at T at its nodes, with the nodes' coordinates
  std::vector<double> leg3( size_t NE, size_t NN, std::vector<std::pair<double,double>>& xy, double rtol = 1e-10, double atol = 1e-12 )
  {
    FFGraph G; OCFEMOL mol( &G ); FFVar t, x, y, U;
    declare( mol, G, t, x, y, U, NE, NN, 1, 2 );
    mol.options.DISPLAY_LEVEL = 0;
    if( !mol.setup() ){ std::cout << "      pde5 leg3: setup failed\n"; return {}; }
    FFGraph H; DAESLV_IDAS S; std::vector<FFVar> xc; FFVar tt; std::string err;
    if( !export_semidiscrete( mol, H, S, FFDom( 0., T, 2, FFDom::LGR, 3 ), xc, tt, &err ) ){ std::cout << "      pde5 leg3: " << err << "\n"; return {}; }
    std::vector<double> x0( mol.semi_discrete().n_states );
    if( !mol.init_seed( 0., x0.data() ) ) return {};
    for( size_t i = 0; i < x0.size(); ++i ) S.update_ref( xc[i], x0[i] );
    auto const& cols = mol.field_columns( U ); auto const& nodes = mol.field_nodes( U );
    FFEval OpE;
    for( size_t c : cols ) S.add_output( OpE( xc[c], tt, T ) );
    S.options.SOLVE.RTOL = rtol;  S.options.SOLVE.ATOL = atol;  S.options.SOLVE.VERBOSE = 0;  S.options.DISPLAY_LEVEL = 0;
    S.options.SOLVE.NMAX = 200000;   // a narrow bump at a tight tolerance: BDF takes many accuracy-limited steps
    if( !S.setup() ){ std::cout << "      pde5 leg3: DAESLV::setup failed: " << S.extract_error() << "\n"; return {}; }
    std::vector<double> P( S.np(), 0. );
    if( S.solve( P ) != DAESLV_IDAS::NORMAL ){ std::cout << "      pde5 leg3: solve failed: " << S.extract_error() << "\n"; return {}; }
    xy.clear();
    for( auto const& nd : nodes ){
      double xv = 0., yv = 0.;
      for( auto const& [d, v] : nd ){ if( d.id() == x.id() ) xv = v; else if( d.id() == y.id() ) yv = v; }
      xy.emplace_back( xv, yv );
    }
    return S.val_function();
  }
  // leg 1: OCFESLV's U at T at the given coordinates
  std::vector<double> leg1( size_t NE, size_t NN, size_t NELT, size_t NODT, std::vector<std::pair<double,double>> const& xy )
  {
    FFGraph G; OCFESLV oc( &G ); FFVar t, x, y, U;
    declare( oc, G, t, x, y, U, NE, NN, NELT, NODT );
    oc.options.SOLVE.MARCHING = false;  oc.options.SOLVE.RES_TOL = 1e-12;  oc.options.DISPLAY_LEVEL = 0;
    if( !oc.setup() ) return {};
    std::vector<double> xv, inp;
    if( !oc.init( xv, inp, nullptr ) ) return {};
    if( !oc.solve( xv.data(), inp.empty()? nullptr: inp.data(), nullptr ).converged ) return {};
    std::vector<double> out;
    for( auto const& [a, b] : xy ){ std::map<FFVar,double,lt_FFVar> pt{ { t, T }, { x, a }, { y, b } }; out.push_back( oc.eval_solution( U, pt ) ); }
    return out;
  }
}

int main()
{
  std::cout.setf( std::ios::unitbuf );
  bool const show_time = std::getenv( "CRONOS_SNAPSHOT" ) == nullptr;
  std::cout << "\n--- DAESLV_molgate: leg 3 (OCFEMOL -> DAESLV) against leg 1 (OCFESLV) ------------------\n";

  Mesh const M( 3, 5 );
  std::cout << "\n  mesh: " << M.nel << " elements x " << M.nnod << " nodes = " << M.n() << " nodes\n\n";

  auto const L3 = leg3( M, 1e-10, 1e-12, true );
  std::cout << "      (the exported DAE's structural index along t: " << L3.index << ")\n";
  check( "leg 3: the exported DAE has structural index <= 1", L3.index >= 0 && L3.index <= 1 );
  check( "leg 3: OCFEMOL set up, exported and integrated", L3.ok );
  if( !L3.ok ){ std::cout << "\n==========  DAESLV_molgate: " << n_pass << " passed, " << ++n_fail << " failed -- FAIL  ==========\n\n"; return 1; }
  std::cout << "      semi-discrete: " << L3.n_states << " unknowns (" << L3.n_diff << " differential), " << L3.n_rows << " rows\n";
  if( show_time ) std::cout << "      leg 3: setup " << std::fixed << std::setprecision(1) << L3.setup_us/1000. << " ms, solve "
                            << L3.solve_us/1000. << " ms\n" << std::defaultfloat;
  // leg 2: the same spatial rows, the same integrator, the same treatment of t0 -- the export's own check
  auto const u2 = leg2( M, 1e-10, 1e-12 );
  double const d23 = maxdiff( u2, L3.unode );
  std::cout << "\n  leg 2 (hand-expanded rows, same IDAS):  |leg2 - leg3| = " << std::scientific << std::setprecision(3) << d23 << std::defaultfloat << "\n";
  auto const u1ref = leg1( M, 8, 4 );
  std::cout << "  |leg1 - leg2| = " << std::scientific << std::setprecision(3) << maxdiff( u1ref, u2 ) << std::defaultfloat
            << "   (legs 1 and 2 of the phase-0 gate)\n";
  auto const L3d = leg3( M, 1e-10, 1e-12, false, true );      // the same export under NEUMANN_DIFFERENTIAL
  if( L3d.ok ){
    std::cout << "  leg 3 under NEUMANN_DIFFERENTIAL (" << L3d.n_diff << " differential): |leg2 - leg3| = "
              << std::scientific << std::setprecision(3) << maxdiff( u2, L3d.unode ) << "   |leg1 - leg3| = "
              << maxdiff( u1ref, L3d.unode ) << std::defaultfloat << "\n";
  }
  // legs 2 and 3 are DIFFERENT spatial discretisations of the same PDE (leg 2 sacrifices one duplicate's equation at each
  // interface, the plan's exact elimination folds both); the claim is spectral convergence under spatial refinement
  std::cout << "\n  spatial ladder, legs 2 vs 3 (nnod 5, elements doubled):\n";
  std::vector<double> d23s;
  for( size_t ne : { 3u, 6u, 12u } ){
    Mesh const Ms( ne, 5 );
    auto const a = leg2( Ms, 1e-11, 1e-13 ); auto const b = leg3( Ms, 1e-11, 1e-13, false );
    double const d = b.ok? maxdiff( a, b.unode ): 1e30;  d23s.push_back( d );
    std::cout << "      " << std::setw(2) << ne << " elements   |leg2 - leg3| = " << std::scientific << std::setprecision(3) << d << std::defaultfloat << "\n";
  }
  check_lt( "legs 2 and 3 converge spectrally: >= 3x per doubling",  3.*d23s[1], d23s[0] );
  check_lt( "   ... and again",                                       3.*d23s[2], d23s[1] );

  // where does the difference live?  node by node at T, the three legs side by side
  std::cout << "\n   node     z        leg1          leg2        leg3(ALG)     leg3(DIF)      1-3dif\n";
  for( size_t k = 0; k < M.n(); ++k )
    std::cout << "   " << std::setw(4) << k << std::fixed << std::setprecision(3) << std::setw(7) << M.z[k]
              << std::setprecision(6) << std::setw(14) << ( u1ref.size() > k? u1ref[k]: 0. )
              << std::setw(14) << ( u2.size() > k? u2[k]: 0. )
              << std::setw(14) << L3.unode[k] << std::setw(14) << ( L3d.ok? L3d.unode[k]: 0. )
              << std::scientific << std::setprecision(1) << std::setw(12)
              << ( ( u1ref.size() > k && L3d.ok )? std::fabs( u1ref[k] - L3d.unode[k] ): 0. ) << std::defaultfloat << "\n";

  // against the closed form: leg 3 under spatial refinement, at the nodes (the exact solution has no spatial
  // discretisation to share, so this is the independent check of the whole construction)
  std::cout << "\n  leg 3 against the closed form, at T = 1 (nnod 5, elements doubled):\n";
  std::vector<double> dex;
  for( size_t ne : { 3u, 6u, 12u } ){
    Mesh const Ms( ne, 5 );
    auto const b = leg3( Ms, 1e-11, 1e-13, false );
    double d = 1e30;
    if( b.ok ){ d = 0.; for( size_t k = 0; k < Ms.n(); ++k ) d = std::max( d, std::fabs( b.unode[k] - u_exact( TEND, Ms.z[k] ) ) ); }
    dex.push_back( d );
    std::cout << "      " << std::setw(2) << ne << " elements   |leg3 - exact| = " << std::scientific << std::setprecision(3) << d << std::defaultfloat << "\n";
  }
  std::cout << "      (u(1,1): exact " << std::setprecision(10) << u_exact( TEND, 1. ) << std::defaultfloat << ")\n";
  check_lt( "leg 3 converges to the closed form: >= 3x per doubling",  3.*dex[1], dex[0] );
  check_lt( "   ... and again",                                         3.*dex[2], dex[1] );
  check_lt( "leg 3 at 12x5 within 1e-5 of the closed form",            dex[2], 1e-5 );

  // the time ladder on leg 1: refining its time collocation must drive |1-3| to the integration tolerance
  std::cout << "\n  time ladder on leg 1 (its t-collocation refined; leg 3 fixed):\n";
  std::vector<double> d13;
  for( size_t nelt : { 2u, 4u, 8u, 16u, 32u } ){
    auto const u1 = leg1( M, nelt, 4 );
    double const d = maxdiff( u1, L3.unode );
    d13.push_back( d );
    std::cout << "      NELT = " << std::setw(2) << nelt << "   |leg1 - leg3| = " << std::scientific << std::setprecision(3) << d << std::defaultfloat << "\n";
  }
  // legs 1 and 3 are the SAME spatial system; they differ by leg 1's time error, which the ladder drives to zero
  // "to the integration tolerance" is relative to the solution's scale (u reaches 3.57 here): the ladder's last ratio
  // flattening is the floor set by the two solvers' own tolerances, not a disagreement
  double umax = 0.; for( double v : L3.unode ) umax = std::max( umax, std::fabs( v ) );
  check_lt( "legs 1 and 3 agree at the finest time mesh (relative)", d13.back() / std::max( 1., umax ), 1e-7 );
  check_lt( "  ... converging: >= 100x over the ladder",          d13.back()*100., d13.front() );

  // ---------------------------------------------------------------- case B
  std::cout << "\n  B. Belousov-Zhabotinskii, flux form, two fields, no-flux ends\n";
  {
    Mesh const Mb( 6, 6, BZ_ZL, BZ_ZR );
    auto const b3 = bz_leg3( Mb, 1e-10, 1e-12 );
    check( "BZ leg 3: OCFEMOL set up, exported and integrated (nonlinear initial block)", !b3.empty() );
    if( !b3.empty() ){
      // the same-system check, by refining leg 1's time collocation in ORDER: spectral, so it reaches the
      // integration tolerance in a few steps where halving elements crawls (x7 per halving)
      auto const b3t = bz_leg3( Mb, 1e-12, 1e-14 );
      std::cout << "\n  leg 1 refined in time ORDER (8 elements; leg 3 at rtol 1e-12):\n";
      std::vector<double> d;
      for( size_t nodt : { 3u, 4u, 5u, 6u } ){
        auto const b1 = bz_leg1( Mb, 8, nodt );  d.push_back( maxdiff( b1, b3t ) );
        std::cout << "      NODT = " << nodt << "   |leg1 - leg3| = " << std::scientific << std::setprecision(3) << d.back() << std::defaultfloat << "\n";
      }
      check_lt( "BZ: legs 1 and 3 converge spectrally in leg 1's time order: >= 10x per order", 10.*d[3], d[2] );
      check_lt( "BZ: legs 1 and 3 agree to 1e-7 at order 6 (fields of order 1)",              d[3], 1e-7 );
      std::cout << "\n  spatial ladder, legs 2 vs 3 (phase 0's meshes):\n";
      std::vector<double> e;
      for( auto const& [ne, nn] : std::vector<std::pair<size_t,size_t>>{ {6,6}, {8,7}, {12,8} } ){
        Mesh const Ms( ne, nn, BZ_ZL, BZ_ZR );
        e.push_back( maxdiff( bz_leg2( Ms, 1e-10, 1e-12 ), bz_leg3( Ms, 1e-10, 1e-12 ) ) );
        std::cout << "      " << std::setw(2) << ne << "x" << nn << "   |leg2 - leg3| = " << std::scientific << std::setprecision(3) << e.back() << std::defaultfloat << "\n";
      }
      check_lt( "BZ: legs 2 and 3 converge under spatial refinement: >= 3x per step", 3.*e[1], e[0] );
      check_lt( "   ... and again",                                                    3.*e[2], e[1] );
    }
  }

  // ---------------------------------------------------------------- case K
  std::cout << "\n  K. PDE5: 2-D manufactured heat equation on (t, x, y), Dirichlet on all sides\n";
  {
    std::vector<std::pair<double,double>> xy;
    auto const k3 = pde5::leg3( 2, 4, xy, 1e-12, 1e-14 );
    check( "2-D leg 3: OCFEMOL set up, exported and integrated", !k3.empty() && k3.size() == xy.size() );
    if( !k3.empty() ){
      // Until 2 Nov 2026 OCFESLV's plan marked an arbitrary subset of tangential-derivative claims redundant on this
      // three-direction model once a time element had more nodes than a spatial one (the aux-implied rule walked the
      // face code's TIME digit instead of the derivative's); the ladder then stopped at order 4.  Fixed (O1): the
      // full order ladder, as for the rod and BZ.
      std::cout << "\n  leg 1 refined in time ORDER (2x2 elements x 4 nodes in space; 4 time elements):\n";
      std::vector<double> d;
      for( size_t nodt : { 3u, 4u, 5u, 6u } ){
        auto const k1 = pde5::leg1( 2, 4, 4, nodt, xy );  d.push_back( maxdiff( k1, k3 ) );
        std::cout << "      NODT = " << nodt << "   |leg1 - leg3| = " << std::scientific << std::setprecision(3) << d.back() << std::defaultfloat << "\n";
      }
      check_lt( "2-D: legs 1 and 3 converge spectrally in leg 1's time order: >= 10x per order", 10.*d[3], d[2] );
      check_lt( "2-D: legs 1 and 3 agree to 1e-7 at order 6",                                   d[3], 1e-7 );
      std::cout << "\n  leg 3 against the exact solution at T, spatial refinement (NN 5):\n";
      std::vector<double> e;
      for( size_t ne : { 2u, 3u, 4u } ){
        std::vector<std::pair<double,double>> xyn;
        auto const u = pde5::leg3( ne, 5, xyn );
        double m = u.empty()? 1e30: 0.;
        for( size_t k = 0; k < u.size(); ++k ) m = std::max( m, std::fabs( u[k] - pde5::exact( pde5::T, xyn[k].first, xyn[k].second ) ) );
        e.push_back( m );
        std::cout << "      " << ne << "x" << ne << " elements (" << u.size() << " nodes)   |leg3 - exact| = " << std::scientific << std::setprecision(3) << m << std::defaultfloat << "\n";
      }
      check_lt( "2-D: leg 3 converges to the exact solution under spatial refinement", e[2], e[0] );
      // 4x4 elements: element lines bounded by interior interfaces at both ends.  At their corners the parent's claim
      // across the other face lifts into the tangential auxiliary's definition, so the tangential claim there is NOT
      // implied; until the rule read that from the plan's receivers it was dropped, and MoL sat 2x further from the
      // exact solution than OCFESLV.  The 2x2 case above cannot see it -- every line there ends at the domain boundary.
      {
        std::vector<std::pair<double,double>> xy4;
        auto const u3 = pde5::leg3( 4, 4, xy4, 1e-12, 1e-14 );
        auto const u1 = pde5::leg1( 4, 4, 4, 6, xy4 );
        double const d4 = maxdiff( u1, u3 );
        std::cout << "\n  4x4 elements x 4 nodes:  |OCFESLV (order 6) - MoL| = " << std::scientific << std::setprecision(3) << d4 << std::defaultfloat << "\n";
        check_lt( "2-D, 4x4: legs 1 and 3 agree to 1e-7 (interior element lines)", d4, 1e-7 );
      }
    }
  }

  std::cout << "\n==========  DAESLV_molgate: " << n_pass << " passed, " << n_fail << " failed -- "
            << ( n_fail? "FAIL": "ALL PASS" ) << "  ==========\n\n";
  return n_fail? 1: 0;
}
