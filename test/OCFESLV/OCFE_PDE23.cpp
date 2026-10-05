// OCFE_PDE23_solve2.cpp
// ---------------------------------------------------------------------------
// Positive-path validation of the distributed-input interface-continuity check
// (_check_input_continuity, run at the top of the numeric OCFESLV::eval()).
//
// The model is a fixed homogeneous Poisson BVP  u'' + w = 0,  u(0)=u(1)=0  on a
// 2-element z-grid (one interior interface at z=0.5).  w is a distributed input
// declared with a per-direction InputContinuity level; we then feed eval() a
// crafted w-profile and assert eval()'s BOOLEAN return:
//
//   * eval() computes a residual for ANY input data, so the ONLY thing that makes
//     it return false is the continuity assertion firing.  Hence
//        eval()==false  <=>  the declared continuity was violated by the data.
//
// Cases (each built on a fresh OCFESLV; only the declaration + data differ):
//   1. C0 declared, JUMP data            -> eval REJECTS  (value jump)
//   2. C0 declared, sloped-CONTINUOUS    -> eval ACCEPTS  (false-reject guard:
//                                           a value-only check would wrongly fail)
//   3. C0 declared, smooth (z^2)         -> eval ACCEPTS
//   4. C1 declared, KINK (C0 not C1)     -> eval REJECTS  (slope jump)
//   5. C1 declared, smooth (z^2)         -> eval ACCEPTS
//   6. no declaration, JUMP data         -> eval ACCEPTS  (control: eval is fine
//                                           with a discontinuous input when nothing
//                                           is declared)
//
// Run for CGL and LGL (endpoint-inclusive => coincident interface nodes).
// ---------------------------------------------------------------------------

#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <string>
#include <functional>

#ifndef OCFE_OCFESLV_HEADER
#define OCFE_OCFESLV_HEADER "ocfeslv.hpp"
#endif
#include OCFE_OCFESLV_HEADER

using namespace mc;

// ---- crafted w(z) profiles, keyed by (coordinate, element index) -----------
static double prof_jump  ( double  , size_t iel ){ return iel==0 ? 0.3 : 0.7; }       // discontinuous
static double prof_sloped( double z, size_t     ){ return z; }                        // continuous, slope 1
static double prof_smooth( double z, size_t     ){ return z*z; }                       // C-infinity
static double prof_kink  ( double z, size_t iel ){ return iel==0 ? z                   // C0, slope 1 left
                                                                 : 0.5 + 2.0*(z-0.5); }// slope 2 right (kink)

using ProfileFn = std::function<double(double,size_t)>;

// decl_level: -1 = no continuity declaration (DISCONTINUOUS default, no check);
//             otherwise an OCFESLV::InputContinuity value.
static bool run_case( FFDom::TYPE coltype, std::string const& family,
                      std::string const& name, int decl_level,
                      ProfileFn profile, bool expect_eval_ok )
{
  size_t const n_el_z = 2, n_nd_z = 4;

  FFGraph DAG;
  FFVar z = DAG.add_var( "z" );
  FFVar u = DAG.add_var( "u(z)" );
  FFVar w = DAG.add_var( "w(z)" );

  FFPartial OpP;
  FFVar PDE = OpP( OpP( u, z ), z ) + w;            // u'' = -w
  FFVar BC0 = u;                                    // u(0)=0
  FFVar BC1 = u;                                    // u(1)=0

  OCFESLV oc( &DAG );
  oc.add_domain( z, FFDom( 0., 1., n_el_z, coltype, n_nd_z ) );
  oc.add_state ( u, {z} );
  oc.update_ref( u, [&]( OCFESLV::t_Coord const& ){ return 0.; } );  // linear: ref value immaterial

  if( decl_level < 0 )
    oc.add_input ( w, {z}, std::optional<double>( 0. ) );                 // no declaration
  else
    oc.add_input ( w, {z},
      std::vector<OCFESLV::InputContinuity>{ static_cast<OCFESLV::InputContinuity>( decl_level ) },
      std::optional<double>( 0. ) );

  oc.reset_evolution_domain();

  int const Z_INT = FFDom::ALL - FFDom::LB - FFDom::UB;
  oc.add_equation( PDE, {z}, {Z_INT},     OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( BC0, {z}, {FFDom::LB}, OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );
  oc.add_equation( BC1, {z}, {FFDom::UB}, OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );

  oc.options.REDUCE.ORDER    = OCFESLV::Options::RED_FULL;
  oc.options.CLASSIFY.MODE        = OCFESLV::Options::CLASS_AUTO;
  oc.options.INTERFACE.TYPE  = OCFESLV::Options::IC_AUTO;
  oc.options.INTERFACE.IMPOSITION = OCFESLV::Options::IC_WEAK;
  oc.options.INTERFACE.SAT_SIGMA0      = 1.0;
  oc.options.DISPLAY_LEVEL   = 0;

  if( !oc.setup() ){
    std::cerr << "ERROR: setup failed (" << family << " / " << name << ")\n";
    return false;
  }

  size_t const nVar = oc.n_colloc_sta(), nEqn = oc.n_colloc_eqn(), nInp = oc.n_colloc_inp();

  // Crafted input data on the (matching) state z-grid.
  auto vn_w = oc.node_colloc( w );
  if( vn_w.size() != nInp || nInp != n_el_z*n_nd_z ){
    std::cerr << "ERROR: w input layout unexpected (nodes=" << vn_w.size()
              << " nInp=" << nInp << ")\n";
    return false;
  }
  std::vector<double> winp( nInp, 0. );
  for( size_t i=0; i<nInp; ++i ){
    size_t const iel = i / n_nd_z;
    winp[i] = profile( vn_w[i][0], iel );
  }

  std::vector<double> var( nVar, 0. );             // state values are irrelevant to the check
  std::vector<double> res( nEqn, 0. );

  std::cout << "  [" << family << "] " << std::left << std::setw(26) << name
            << "  eval=" << std::flush;
  bool const eval_ok = oc.eval( res.data(), nullptr, var.data(), winp.data(), nullptr );
  bool const pass    = ( eval_ok == expect_eval_ok );

  std::cout << ( eval_ok ? "accept" : "REJECT" )
            << "  expected=" << ( expect_eval_ok ? "accept" : "REJECT" )
            << "   -> " << ( pass ? "PASS" : "FAIL" ) << "\n";
  return pass;
}

int main()
{
  std::cout << "================================================================\n";
  std::cout << "  PDE23  input interface-continuity eval-check (positive path)\n";
  std::cout << "         model: u'' + w = 0, u(0)=u(1)=0, 2 z-elements\n";
  std::cout << "================================================================\n";

  bool all = true;
  for( auto const& [tag,typ] : std::vector<std::pair<std::string,FFDom::TYPE>>{
         { "CGL", FFDom::CGL }, { "LGL", FFDom::LGL } } )
  {
    std::cout << "\n--- " << tag << " ---\n";
    all &= run_case( typ, tag, "C0 + jump",          OCFESLV::CONTINUOUS_C0, prof_jump,   false );
    all &= run_case( typ, tag, "C0 + sloped-cont",   OCFESLV::CONTINUOUS_C0, prof_sloped, true  );
    all &= run_case( typ, tag, "C0 + smooth",        OCFESLV::CONTINUOUS_C0, prof_smooth, true  );
    all &= run_case( typ, tag, "C1 + kink",          OCFESLV::SMOOTH_C1,     prof_kink,   false );
    all &= run_case( typ, tag, "C1 + smooth",        OCFESLV::SMOOTH_C1,     prof_smooth, true  );
    all &= run_case( typ, tag, "no-decl + jump",     -1,                   prof_jump,   true  );
  }

  std::cout << "\n  Overall: " << ( all ? "ALL PASS" : "SOME FAILED" ) << "\n";
  return all ? 0 : 1;
}
