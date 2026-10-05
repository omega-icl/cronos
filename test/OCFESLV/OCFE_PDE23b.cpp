// OCFE_PDE23b.cpp
// ---------------------------------------------------------------------------
// Regression for the IC_FLUX C1/flux safety guard in _resolve_interface_type.
//
// The guard downgrades an explicit IC_FLUX to IC_VALUE, with a warning, in any
// direction that carries a DISTRIBUTED input which may jump there -- i.e. one
// collocated in that direction with an interior interface (n_elem >= 2) and NOT
// declared CONTINUOUS_C0 / SMOOTH_C1.  A jumping input forces the matched
// state's first derivative to jump, so imposing C1/flux continuity would
// contradict the true (kinked) solution; IC_VALUE (C0) is always safe.
//
// Model: first-order advection  d_t u + c d_z u = w  on (t,z), t evolution, z
// spatial with one interior interface (2 z-elements, LGL => coincident interface
// nodes).  First order in z => NO reduce_order aux link in z, so the aux-link
// downgrade does not pre-empt this guard; z is not the evolution face either.
// Hence, in the explicit-IC_FLUX branch, THIS guard is the only thing that can
// turn IC_FLUX into IC_VALUE for the z-direction -- a clean isolation.
//
// The block is set up ONCE with IC_AUTO (always square / well-posed).  The guard
// is then queried as a pure resolution: flip options.INTERFACE.TYPE to IC_FLUX
// and call the PUBLIC diagnostic accessor resolved_interface_type(), which routes
// through the guarded resolver WITHOUT re-running setup.  This isolates the
// guard's decision and avoids the (expected) row over-count of a genuine IC_FLUX
// setup on a first-order block -- we are testing the resolution verdict, which is
// exactly the thing the guard changes.
//
// Cases (fresh OCFESLV each; only the w input differs):
//   A  w distributed over {z}, UNDECLARED (DISCONTINUOUS default)  -> IC_VALUE + warning
//   B  w distributed over {z}, declared CONTINUOUS_C0             -> IC_FLUX  (no warning)
//   C  w distributed over {t} only, undeclared                    -> IC_FLUX  (no warning)
//        Q1 check: an input NOT distributed in z (here: distributed only in t;
//        the same code path also skips constants, which never enter _mInp) must
//        be ignored -- no C1 downgrade may be triggered by it.
//   D  w distributed over {z}, declared SMOOTH_C1                 -> IC_FLUX  (no warning)
//        (C1 >= C0: a smoother-than-C0 declaration also suppresses the guard.)
//
// Build: same flags as the rest of the OCFE suite; needs only the header.
// ---------------------------------------------------------------------------

#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <string>
#include <sstream>
#include <optional>

#ifndef OCFE_OCFESLV_HEADER
#define OCFE_OCFESLV_HEADER "ocfeslv.hpp"
#endif
#include OCFE_OCFESLV_HEADER

using namespace mc;

// RAII std::cerr capture -- the guard's warning is emitted on cerr (ungated).
struct CerrCapture
{
  std::ostringstream oss;
  std::streambuf*    old;
  CerrCapture()  : old( std::cerr.rdbuf( oss.rdbuf() ) ) {}
  ~CerrCapture() { std::cerr.rdbuf( old ); }
  std::string str() const { return oss.str(); }
};

enum class InpDom { Z, T };   // direction the input w is distributed over

static const char* iftype_name( OCFESLV::Options::InterfaceType it )
{
  switch( it ){
    case OCFESLV::Options::IC_VALUE:  return "IC_VALUE";
    case OCFESLV::Options::IC_FLUX:   return "IC_FLUX";
    case OCFESLV::Options::IC_UPWIND: return "IC_UPWIND";
    case OCFESLV::Options::IC_AUTO:   return "IC_AUTO";
  }
  return "?";
}

// decl_level: -1 => no continuity declaration (DISCONTINUOUS default);
//             else an OCFESLV::InputContinuity value (CONTINUOUS_C0 / SMOOTH_C1).
static bool run_case( std::string const& name, InpDom inpdom, int decl_level,
                      OCFESLV::Options::InterfaceType expect_type, bool expect_warn )
{
  size_t const n_el_t = 2, n_nd_t = 4;
  size_t const n_el_z = 2, n_nd_z = 4;      // one interior z-interface at z=0.5
  double const c = 0.6;                      // forward advection (inflow at z=LB)

  FFGraph DAG;
  FFVar t = DAG.add_var( "t" );
  FFVar z = DAG.add_var( "z" );
  FFVar u = DAG.add_var( "u(t,z)" );
  FFVar w = DAG.add_var( "w" );

  FFPartial OpP;
  FFVar PDE = OpP( u, t ) + c*OpP( u, z ) - w;   // first order in z (and t): no aux link in z
  FFVar REF = u;                                  // IC / inflow BC (ref value immaterial for setup)

  OCFESLV oc( &DAG );
  oc.add_domain( t, FFDom( 0., 1., n_el_t, FFDom::LGR, n_nd_t ) );
  oc.add_domain( z, FFDom( 0., 1., n_el_z, FFDom::LGL, n_nd_z ) );
  oc.add_state ( u, {t,z} );
  oc.update_ref( u, [&]( OCFESLV::t_Coord const& ){ return 0.; } );

  std::vector<FFVar> const wdom = ( inpdom == InpDom::Z )
                                ? std::vector<FFVar>{ z }
                                : std::vector<FFVar>{ t };
  if( decl_level < 0 )
    oc.add_input( w, wdom, std::optional<double>( 0. ) );
  else
    oc.add_input( w, wdom,
      std::vector<OCFESLV::InputContinuity>{ static_cast<OCFESLV::InputContinuity>( decl_level ) },
      std::optional<double>( 0. ) );

  oc.set_evolution_domain( t );

  int const T_INT = FFDom::ALL - FFDom::LB;
  oc.add_equation( PDE, {t,z}, {T_INT,     FFDom::ALL-FFDom::LB},
                   OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );  // t>0, z>LB (interior + outflow UB)
  oc.add_equation( REF, {t,z}, {FFDom::LB, FFDom::ALL},
                   OCFESLV::EqnOptions( OCFESLV::EqnRole::INITIAL, 0 ) );   // t=LB initial
  oc.add_equation( REF, {t,z}, {T_INT,     FFDom::LB},
                   OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );  // z=LB inflow

  oc.options.REDUCE.ORDER    = OCFESLV::Options::RED_FULL;
  oc.options.CLASSIFY.MODE        = OCFESLV::Options::CLASS_AUTO;
  oc.options.INTERFACE.TYPE  = OCFESLV::Options::IC_AUTO;   // set up with AUTO (always square)
  oc.options.INTERFACE.IMPOSITION = OCFESLV::Options::IC_WEAK;
  oc.options.INTERFACE.SAT_SIGMA0      = 1.0;
  oc.options.DISPLAY_LEVEL   = 0;

  if( !oc.setup() ){
    std::cout << "  " << std::left << std::setw(30) << name
              << "  setup FAILED  -> FAIL\n";
    return false;
  }

  // Query the guard as a pure resolution: flip to IC_FLUX, call the diagnostic
  // accessor (routes through _resolve_interface_type -> guard), then restore.
  oc.options.INTERFACE.TYPE = OCFESLV::Options::IC_FLUX;
  OCFESLV::Options::InterfaceType rz = OCFESLV::Options::IC_AUTO;
  std::string warn;
  {
    CerrCapture cap;
    rz   = oc.resolved_interface_type( 0, z, OCFESLV::EqnRole::INTERIOR, /*weak=*/true );
    warn = cap.str();
  }
  oc.options.INTERFACE.TYPE = OCFESLV::Options::IC_AUTO;

  bool const warned  = ( warn.find( "IC_FLUX" )     != std::string::npos )
                    && ( warn.find( "Downgrading" ) != std::string::npos );
  bool const type_ok = ( rz == expect_type );
  bool const warn_ok = ( warned == expect_warn );
  bool const pass    = type_ok && warn_ok;

  std::cout << "  " << std::left << std::setw(30) << name
            << "  z-resolved=" << std::setw(9) << iftype_name( rz )
            << " (expect " << std::setw(9) << iftype_name( expect_type ) << ")"
            << "  warn="   << ( warned ? "yes" : "no " )
            << " (expect " << ( expect_warn ? "yes" : "no " ) << ")"
            << "   -> "    << ( pass ? "PASS" : "FAIL" ) << "\n";
  return pass;
}

int main()
{
  std::cout << "================================================================\n";
  std::cout << "  PDE23b IC_FLUX C1/flux safety guard (_resolve_interface_type)\n";
  std::cout << "         model: d_t u + c d_z u = w, first order in z (no aux link),\n";
  std::cout << "         one interior z-interface; guard queried via resolved_interface_type\n";
  std::cout << "================================================================\n\n";

  bool all = true;
  all &= run_case( "A z-input undeclared (jumpy)", InpDom::Z, -1,
                   OCFESLV::Options::IC_VALUE, true  );
  all &= run_case( "B z-input CONTINUOUS_C0",      InpDom::Z, OCFESLV::CONTINUOUS_C0,
                   OCFESLV::Options::IC_FLUX,  false );
  all &= run_case( "C t-input (not in z) undecl",  InpDom::T, -1,
                   OCFESLV::Options::IC_FLUX,  false );
  all &= run_case( "D z-input SMOOTH_C1",          InpDom::Z, OCFESLV::SMOOTH_C1,
                   OCFESLV::Options::IC_FLUX,  false );

  std::cout << "\n  Overall: " << ( all ? "ALL PASS" : "SOME FAILED" ) << "\n";
  return all ? 0 : 1;
}
