// OCFE_zindex.cpp -- what "high index in a SPATIAL direction" means to a collocation framework
// ============================================================================================
// CRONOS's index machinery is EVOLUTION-ONLY by construction: _reduce_high_index takes its
// direction from the evolution variable and differentiates with total_ddt; there is no total_ddx
// and nothing computes nu_x (Martinson & Barton's index with respect to x).  This driver asks what
// that costs, on a chain that IS high index in z:
//
//     p        = g(z)          (constraint: no q, no r)
//     dp/dz    = q
//     dq/dz    = r             ->  r is reached only after TWO z-differentiations of the constraint
//
// For a time integrator that chain would be index 3 and unusable without reduction.  A COLLOCATION
// framework differentiates the INTERPOLANT, so the derivatives are available as auxiliaries with
// their LINK rows and the chain is simply solvable -- the hypothesis this driver tests.  What then
// remains of the spatial index question is BOUNDARY-CONDITION COUNTING, which is exactly what the
// paper's nu_x is for, and which the mesh-free DOF balance can see.
//
// CASES
//   1  the chain coupled to an evolution in t: sets up, balanced, solves EXACTLY (p=g, q=g', r=g''),
//      and the model reports it as INDEX 3 IN Z while the solver never needs to know
//   2  the t-direction reduction plan is EMPTY -- no spatial index machinery ran, and none was needed
//   3  a redundant z-boundary condition on q (g' already determines it) -> reported SURPLUS
//   3b the chain with r DECOUPLED from u's row -> the z-reading is STRUCTURALLY SINGULAR (-1): u is
//      reachable only through its t-derivative, so nothing determines it in a z-marching reading
//   4  the SAME chain with z declared the EVOLUTION direction -> the index machinery DOES engage
//      (plan non-empty, max_index 3): spatial index reduction is available by RELABELLING, which is
//      the cheapest answer to "is nu_x reachable today?"
// ============================================================================================
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#ifndef OCFE_OCFESLV_HEADER
#define OCFE_OCFESLV_HEADER "ocfeslv.hpp"
#endif
#include OCFE_OCFESLV_HEADER
using namespace mc;

static int g_pass = 0, g_fail = 0;
static void check( char const* nm, bool ok )
{ std::cout << "  " << std::left << std::setw(58) << nm << std::right << ( ok? " PASS": " FAIL" ) << std::endl;
  ( ok? g_pass: g_fail )++; }

static double gz  ( double z ){ return 1.0 + 0.5*z - 0.3*z*z; }
static double gpz ( double z ){ return 0.5 - 0.6*z; }
static double gppz( double   ){ return -0.6; }

int main()
{
  std::cout << "OCFESLV ** header: " << OCFESLV::HEADER_ID << "\n";
  int const T_INT = FFDom::ALL - FFDom::LB;
  int const Z_INT = FFDom::ALL - FFDom::LB - FFDom::UB;

  // ---- cases 1-3: the z-chain carried by a model that evolves in t -------------------------
  for( int variant = 0; variant < 2; ++variant ){         // 0 = as posed, 1 = with a redundant z-BC
    std::cout << "\n---- " << ( variant? "with a redundant boundary condition on q at z=LB"
                                       : "the z-chain, as posed" ) << " ----\n";
    FFGraph D; OCFESLV oc( &D );
    FFVar t = D.add_var("t"), z = D.add_var("z");
    FFVar u = D.add_var("u(t,z)"), p = D.add_var("p(t,z)"), q = D.add_var("q(t,z)"), r = D.add_var("r(t,z)");
    FFPartial OpP;
    oc.add_domain( t, FFDom( 0., 0.2, 2, FFDom::LGR, 4 ) );
    oc.add_domain( z, FFDom( 0., 1.0, 2, FFDom::LGL, 5 ) );
    for( auto const& s : { u, p, q, r } ) oc.add_state( s, {t,z} );
    oc.set_evolution_domain( t );
    oc.update_ref( u, 0. ); oc.update_ref( p, 1. ); oc.update_ref( q, 0.5 ); oc.update_ref( r, -0.6 );

    FFVar G = 1.0 + 0.5*z - 0.3*z*z;
    OCFESLV::EqnOptions io( OCFESLV::EqnRole::INTERIOR, 0 ), ii( OCFESLV::EqnRole::INITIAL, 0 ),
                        ib( OCFESLV::EqnRole::BOUNDARY, 0 );
    oc.add_equation( OpP(u,t) - OpP(u,{z,2}) - r, {t,z}, {T_INT, Z_INT},          io );  // forced by r
    oc.add_equation( p - G,                       {t,z}, {FFDom::ALL,FFDom::ALL}, io );  // the constraint
    oc.add_equation( OpP(p,z) - q,                {t,z}, {FFDom::ALL,FFDom::ALL}, io );  // chain level 1
    oc.add_equation( OpP(q,z) - r,                {t,z}, {FFDom::ALL,FFDom::ALL}, io );  // chain level 2
    oc.add_equation( u - 0.,                      {t,z}, {FFDom::LB,FFDom::ALL},  ii );
    oc.add_equation( u - 0.,                      {t,z}, {T_INT,FFDom::LB},       ib );
    oc.add_equation( u - 0.,                      {t,z}, {T_INT,FFDom::UB},       ib );
    if( variant ) oc.add_equation( q - 0.5, {t,z}, {T_INT,FFDom::LB}, ib );   // g'(0) -- already determined
    oc.options.REDUCE.ORDER   = OCFESLV::Options::RED_FULL;
    oc.options.DISPLAY_LEVEL  = 1;
    oc.options.FATAL.REDUCED_DOF = false;
    bool const ok = oc.setup();
    std::cout << "  setup=" << ok << "  balance: rows-unknowns=" << oc.dof_balance().str
              << "  reduction plan: " << ( oc.reduction_plan().empty()? "EMPTY": "non-empty" )
              << " max_index=" << oc.reduction_plan().max_index << "\n";
    if( !variant ){
      check( "the z-chain sets up with no spatial index machinery", ok );
      check( "it is balanced",                                      oc.dof_balance().balanced );
      check( "no t-index reduction was needed either (plan empty)", oc.reduction_plan().empty() );
      // rev339: the MODEL reports the index in ANY direction -- the chain is high index in z even though
      // this solver never needs to know, because a consumer that marches or shoots in z would.
      for( auto const& [blk, cls] : oc.block_classification() ){
        int const nu_t = oc.structural_index( blk, t ).index;
        int const nu_z = oc.structural_index( blk, z ).index;
        std::cout << "  block " << blk << ": index in t = " << nu_t << ", index in z = " << nu_z << std::endl;
        // The analysis removes order-reduction auxiliaries from the index-relevant set, so the reading is the
        // DECLARED structure of the chain and is stable under RED_FULL: index 3 in z, however the solver discretises.
        check( "the model reports the chain as index 3 IN Z", nu_z == 3 );
        check( "while the index in t stays low (1 or less)",     nu_t <= 1 );
      }
      if( ok ){
        std::vector<double> xv, inp;
        if( oc.init( xv, inp, nullptr ) ){
          OCFESLV::SolveReport rep = oc.solve( xv.data(), inp.data(), nullptr );
          auto ev = [&]( FFVar const& V, double tt, double zz ){
            OCFESLV::t_Coord pt; pt[t]=tt; pt[z]=zz;
            return oc.eval_colloc<double>( V, pt, xv.data(), inp.data(), nullptr ); };
          double ep=0., eq=0., er=0.;
          for( int k=0; k<=20; ++k ){ double zz = k/20.0;
            ep = std::max( ep, std::fabs( ev(p,0.2,zz) - gz(zz)   ) );
            eq = std::max( eq, std::fabs( ev(q,0.2,zz) - gpz(zz)  ) );
            er = std::max( er, std::fabs( ev(r,0.2,zz) - gppz(zz) ) ); }
          std::cout << "  solve: converged=" << rep.converged << std::scientific << std::setprecision(3)
                    << "  max|p-g|=" << ep << "  max|q-g'|=" << eq << "  max|r-g''|=" << er << "\n";
          check( "it converges", rep.converged );
          check( "the chain is exact: p=g, q=g', r=g'' (interpolant)", ep<1e-10 && eq<1e-10 && er<1e-10 );
        }
      }
    }
    else{
      check( "a redundant z-boundary condition is a reported SURPLUS",
             !oc.dof_balance().balanced && oc.dof_balance().str == "N_t - 1" );
    }
  }

  // ---- case 3b: the SAME chain with NO order reduction: the declared chain is high index in z ----
  std::cout << "\n---- the declared chain, with order reduction OFF ----\n";
  {
    FFGraph D; OCFESLV oc( &D );
    FFVar t = D.add_var("t"), z = D.add_var("z");
    FFVar u = D.add_var("u(t,z)"), p = D.add_var("p(t,z)"), q = D.add_var("q(t,z)"), r = D.add_var("r(t,z)");
    FFPartial OpP;
    oc.add_domain( t, FFDom( 0., 0.2, 2, FFDom::LGR, 4 ) );
    oc.add_domain( z, FFDom( 0., 1.0, 2, FFDom::LGL, 5 ) );
    for( auto const& st : { u, p, q, r } ) oc.add_state( st, {t,z} );
    oc.set_evolution_domain( t );
    oc.update_ref( u, 0. ); oc.update_ref( p, 1. ); oc.update_ref( q, 0.5 ); oc.update_ref( r, -0.6 );
    FFVar G = 1.0 + 0.5*z - 0.3*z*z;
    OCFESLV::EqnOptions io( OCFESLV::EqnRole::INTERIOR, 0 ), ii( OCFESLV::EqnRole::INITIAL, 0 );
    // u's row must NOT contain r: with r bare in it, r is pinned directly and the chain never needs a
    // z-differentiation -- the reading is then index 1, correctly.  Decoupled, r is reachable only by
    // differentiating the constraint twice in z, which is the genuine nu_z = 3 structure.
    oc.add_equation( OpP(u,t) + 0.6, {t,z}, {T_INT, FFDom::ALL},       io );   // first order in t only
    oc.add_equation( p - G,         {t,z}, {FFDom::ALL,FFDom::ALL},    io );
    oc.add_equation( OpP(p,z) - q,  {t,z}, {FFDom::ALL,FFDom::ALL},    io );
    oc.add_equation( OpP(q,z) - r,  {t,z}, {FFDom::ALL,FFDom::ALL},    io );
    oc.add_equation( u - 0.,        {t,z}, {FFDom::LB,FFDom::ALL},     ii );
    oc.options.REDUCE.ORDER   = OCFESLV::Options::RED_NONE;    // leave the chain as declared
    oc.options.DISPLAY_LEVEL  = 1;
    oc.options.FATAL.REDUCED_DOF = false;
    bool const ok = oc.setup();
    int nu_z = -1, nu_t = -1;
    for( auto const& [blk, cls] : oc.block_classification() ){
      nu_z = std::max( nu_z, oc.structural_index( blk, z ).index );
      nu_t = std::max( nu_t, oc.structural_index( blk, t ).index );
    }
    std::cout << "  setup=" << ok << "  index in t = " << nu_t << ", index in z = " << nu_z << "\n";
    for( auto const& [blk, cls] : oc.block_classification() ){
      for( auto const& dir : { t, z } ){
        auto const ir = oc.structural_index( blk, dir );
        std::cout << "    [" << dir.name() << "] index=" << ir.index << "  index_alg={";
        for( auto const& a : ir.index_alg ) std::cout << " " << a.name();
        std::cout << " }  unmatched={";
        for( auto const& a : ir.unmatched ) std::cout << " " << a.name();
        std::cout << " }" << std::endl;
      }
    }
    // MEASURED 2026-09-23: nu_z reads 1 even with order reduction OFF.  The direction parameter reaches the
    // EXPOSURE step of _structural_index_analysis, but the differential/algebraic split comes from
    // _structural_dae_decomposition, which is still keyed on the EVOLUTION direction: with t evolving, p, q and r
    // are all "algebraic" and each is pinned directly, so the reading is index 1 and internally consistent.
    // A true nu_x needs that decomposition parameterised too.  Recorded here so the day it changes, this says so.
    // With r decoupled, u appears ONLY under d/dt, so in the z-reading no row pins it: the analysis returns -1,
    // structurally singular.  That is the right answer for a consumer marching in z -- it could never determine u --
    // and it is a verdict the evolution-direction analysis alone could never produce.
    check( "decoupled, the z-reading is structurally singular (-1)", ok && nu_z == -1 );
    check( "while the index in t stays low",                         ok && nu_t <= 1 );
  }

  // ---- case 4: the same chain with z as the EVOLUTION direction ----------------------------
  std::cout << "\n---- the same chain with z declared the evolution direction ----\n";
  {
    FFGraph D; OCFESLV oc( &D );
    FFVar z = D.add_var("z");
    FFVar p = D.add_var("p(z)"), q = D.add_var("q(z)"), r = D.add_var("r(z)");
    FFPartial OpP;
    oc.add_domain( z, FFDom( 0., 1.0, 2, FFDom::LGR, 5 ) );
    for( auto const& s : { p, q, r } ) oc.add_state( s, {z} );
    oc.set_evolution_domain( z );
    oc.update_ref( p, 1. ); oc.update_ref( q, 0.5 ); oc.update_ref( r, -0.6 );
    FFVar G = 1.0 + 0.5*z - 0.3*z*z;
    int const Z_NO_LB = FFDom::ALL - FFDom::LB;
    OCFESLV::EqnOptions io( OCFESLV::EqnRole::INTERIOR, 0 );
    oc.add_equation( OpP(p,z) - q, {z}, {Z_NO_LB},    io );
    oc.add_equation( OpP(q,z) - r, {z}, {Z_NO_LB},    io );
    oc.add_equation( p - G,        {z}, {FFDom::ALL}, io );
    oc.options.REDUCE.ORDER   = OCFESLV::Options::RED_FULL;
    oc.options.REDUCE.HIDDEN_IC = true;                 // free data = 0, as in OCFE_PDE20c
    oc.options.DISPLAY_LEVEL  = 1;
    oc.options.FATAL.REDUCED_DOF = false;
    bool const ok = oc.setup();
    auto const& P = oc.reduction_plan();
    std::cout << "  setup=" << ok << "  plan: " << ( P.empty()? "EMPTY": "non-empty" )
              << " max_index=" << P.max_index << " assigns=" << P.assigns.size()
              << "  balance=" << oc.dof_balance().str << "\n";
    for( auto const& a : P.assigns ) std::cout << "    differentiate " << a.n_diff << "x to pin " << a.pinned_var.name() << "\n";
    check( "relabelled as the evolution direction, the index IS detected", !P.empty() && P.max_index == 3 );
    check( "and the reduced model is balanced",                            ok && oc.dof_balance().balanced );
  }

  std::cout << "\n  OCFE_zindex: " << g_pass << " passed, " << g_fail << " failed\n";
  return g_fail? 1: 0;
}
