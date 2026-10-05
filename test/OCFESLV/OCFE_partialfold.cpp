// OCFE_partialfold.cpp -- FFPartial's nested-derivative FOLD (rev293, CRONOS_FOLD_PARTIALS).
//   Without the knob a nested d/dy(d/dx u) is TWO nodes and every spelling of the same derivative is a distinct
//   node; with it, one node, hash-consed, so spelling stops mattering.  The driver asserts the DAG-level facts
//   in both states and, on a model, that the reduction mints ONE auxiliary per distinct derivative and the
//   answer is unchanged.
#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <string>
#include <cstdlib>
#include <cmath>
#include "ffunc.hpp"
#ifndef OCFE_OCFESLV_HEADER
#define OCFE_OCFESLV_HEADER "ocfeslv.hpp"
#endif
#include OCFE_OCFESLV_HEADER

using namespace mc;
using mc::FFGraph; using mc::FFVar; using mc::FFPartial;
typedef OCFESLV::Options O;
static int g_pass=0, g_fail=0;
static void expect( std::string const& what, bool ok, std::string const& detail="" )
{ std::cout << "  " << std::left << std::setw(66) << what << ( ok?"PASS":"FAIL" ) << ( detail.empty()?"":"   "+detail ) << "\n"; (ok?g_pass:g_fail)++; }
static std::string idof( FFVar const& v ){ std::ostringstream o; o << v.id().second; return o.str(); }

int main()
{
  bool const fold = FFPartial::FOLD_NESTED;   // rev297: a public option, defaulted from the environment
  std::cout << "================================================================\n"
            << "  FFPartial nested-derivative fold   CRONOS_FOLD_PARTIALS=" << (fold?1:0) << "\n"
            << "  header: " << OCFE_OCFESLV_HEADER << "\n"
            << "================================================================\n";

  // ---- DAG-level -------------------------------------------------------------------------------------------
  { FFGraph DAG;
    FFVar x = DAG.add_var("x"), y = DAG.add_var("y"), u = DAG.add_var("u(x,y)");
    FFPartial OpP;
    FFVar const uxx_nested = OpP( OpP( u, x ), x );
    FFVar const uxx_direct = OpP( u, {x,2} );
    FFVar const uxy        = OpP( OpP( u, x ), y );
    FFVar const uyx        = OpP( OpP( u, y ), x );
    FFVar const ux         = OpP( u, x );
    std::cout << "  ids: d2u/dx2 nested=" << idof(uxx_nested) << " direct=" << idof(uxx_direct)
              << " | dxdy=" << idof(uxy) << " dydx=" << idof(uyx) << " | du/dx=" << idof(ux) << "\n";
    expect( "nested d/dx(d/dx u) == direct d2u/dx2 (same node)", ( uxx_nested.id() == uxx_direct.id() ) == fold,
            fold?"folded":"not folded (expected without the knob)" );
    expect( "mixed d/dy(d/dx u) == d/dx(d/dy u) (Clairaut, same node)", ( uxy.id() == uyx.id() ) == fold );
    expect( "first-order node is untouched by the fold", ux.id() != uxx_direct.id() );
    auto sg = DAG.subgraph( 1, &uxx_nested );
    size_t np=0; for( auto const* op : sg.l_op ) if( op && op->sameid( typeid(FFPartial) ) ) ++np;
    expect( "nested expression holds ONE FFPartial when folded, two otherwise", np == (fold?1u:2u), "np="+std::to_string(np) );
  }

  // ---- model level: one auxiliary per distinct derivative, same answer -------------------------------------
  //   u_t = a u_xx, written with the second derivative spelled BOTH ways in two equations of the same block:
  //   the PDE uses the nested spelling, a diagnostic output uses the direct one.
  { double const a=0.1, xf=1.0, tf=0.5;
    auto run=[&]( bool nested_in_out, double& err, size_t& naux )->bool{
      FFGraph DAG;
      FFVar t=DAG.add_var("t"), x=DAG.add_var("x"), u=DAG.add_var("u(t,x)"), s=DAG.add_var("s(t,x)");
      FFPartial OpP;
      OCFESLV oc( &DAG );
      oc.add_domain( t, FFDom( 0., tf, 2, FFDom::LGR, 4 ) );
      oc.add_domain( x, FFDom( 0., xf, 3, FFDom::LGL, 5 ) );
      double const k = 3.14159265358979323846/2./xf;
      auto uex=[&]( double tt, double xx ){ return std::exp(-a*k*k*tt)*std::sin(k*xx); };
      oc.add_state( u, { t, x } );
      oc.update_ref( u, [&]( OCFESLV::t_Coord const& c ){ return uex( c.at(t), c.at(x) ); } );
      oc.add_state( s, { t, x } );
      oc.update_ref( s, [&]( OCFESLV::t_Coord const& c ){ return -k*k*uex( c.at(t), c.at(x) ); } );
      oc.set_evolution_domain( t );
      FFVar const uxx_nested = OpP( OpP( u, x ), x );
      FFVar const uxx_direct = OpP( u, {x,2} );
      FFVar PDE = OpP( u, t ) - a*uxx_nested;                       // one spelling here ...
      FFVar ALG = s - ( nested_in_out ? uxx_nested : uxx_direct );  // ... and the other (or the same) here
      FFVar INI = u - sin(k*x), BCL = u, BCU = OpP( u, x );
      typedef OCFESLV::EqnOptions EO; typedef OCFESLV::EqnRole ER;
      oc.add_equation( PDE, {t,x}, { FFDom::ALL-FFDom::LB, FFDom::ALL-FFDom::LB-FFDom::UB }, EO( ER::INTERIOR, 0 ) );
      oc.add_equation( ALG, {t,x}, { FFDom::ALL, FFDom::ALL }, EO( ER::INTERIOR, 0 ) );
      oc.add_equation( INI, {t,x}, { FFDom::LB, FFDom::ALL }, EO( ER::INITIAL, 0 ) );
      oc.add_equation( BCL, {t,x}, { FFDom::ALL-FFDom::LB, FFDom::LB }, EO( ER::BOUNDARY, 0 ) );
      oc.add_equation( BCU, {t,x}, { FFDom::ALL-FFDom::LB, FFDom::UB }, EO( ER::BOUNDARY, 0 ) );
      oc.options.REDUCE.ORDER = O::RED_FULL; oc.options.CLASSIFY.MODE = O::CLASS_AUTO;
      oc.options.INTERFACE.IMPOSITION = O::IC_WEAK; oc.options.DISPLAY_LEVEL = 0;
      if( !oc.setup() ) return false;
      naux = oc.auxiliary_states().size();
      std::vector<double> xv, inp;
      if( !oc.init( xv, inp, nullptr ) ) return false;
      auto rep = oc.solve( xv.data(), inp.data(), nullptr );
      if( !rep.converged ) return false;
      err = 0.; size_t off = 0;
      for( auto const& st : oc.states_colloc() ){
        auto nodes = oc.node_colloc( st );
        if( st.name() == std::string("u(t,x)") )
          for( size_t i=0; i<nodes.size(); ++i ) err = std::max( err, std::fabs( xv[off+i] - uex(nodes[i][0],nodes[i][1]) ) );
        off += nodes.size(); }
      return true; };
    double e_same=0., e_mixed=0.; size_t n_same=0, n_mixed=0;
    bool const ok1 = run( true,  e_same,  n_same  );    // both equations spell it the same way
    bool const ok2 = run( false, e_mixed, n_mixed );    // the two equations spell it differently
    std::cout << "  same spelling : setup " << (ok1?"ok":"FAILED") << "  auxiliaries=" << n_same  << "  max|u-u*|=" << std::scientific << std::setprecision(4) << e_same  << "\n";
    std::cout << "  mixed spelling: setup " << (ok2?"ok":"FAILED") << "  auxiliaries=" << n_mixed << "  max|u-u*|=" << e_mixed << "\n";
    expect( "both spellings set up and solve", ok1 && ok2 );
    // MEASURED 2026-09-18: one auxiliary per distinct derivative in BOTH knob states -- _reduce_order peels to
    // the first-order chain and keys its reuse on the INNER node, which both spellings already share.  So the
    // fold's value is not de-duplication here; it is node identity (asserted above) and the bare-state-operand
    // blind spot in the scans.  The invariant worth protecting is spelling- AND knob-independence:
    expect( "one auxiliary per distinct derivative, whatever the spelling or the knob", n_mixed == n_same && n_same == 1,
            "same="+std::to_string(n_same)+" mixed="+std::to_string(n_mixed) );
    expect( "the answer does not depend on the spelling", ok1 && ok2 && std::fabs(e_same-e_mixed) <= 1e-12*std::max(1.,e_same),
            "same="+std::to_string(e_same)+" mixed="+std::to_string(e_mixed) );
  }

  // ---- CELL 3: the option is settable IN PROCESS (rev297) --------------------------------------------------
  //   The old predicate read the environment once into a static local, so a single run could only exercise one
  //   state.  Each case gets its OWN FFGraph: nodes already inserted keep their identity, so the flag must be
  //   set before the expressions are built -- which is what this cell also documents.
  std::cout << "\n---- CELL 3: FOLD_NESTED settable in process ----\n";
  { bool const saved = FFPartial::FOLD_NESTED;
    auto nested_is_direct = []( bool flag ){
      FFPartial::FOLD_NESTED = flag;
      FFGraph G; FFVar x = G.add_var("x"), u = G.add_var("u(x)"); FFPartial OpP;
      return OpP( OpP( u, x ), x ).id() == OpP( u, {x,2} ).id(); };
    bool const on = nested_is_direct( true ), off = nested_is_direct( false );
    FFPartial::FOLD_NESTED = saved;
    expect( "C3 ON folds and OFF does not, in the same process", on && !off,
            std::string("on=") + (on?"folded":"not") + " off=" + (off?"folded":"not") );
    expect( "C3 the flag is restored", FFPartial::FOLD_NESTED == saved );
  }

  // ---- CELL 4: the UNIFORM vector fold (rev297) -------------------------------------------------------------
  //   One FFOp carries ONE _Indep for all its outputs, so a vector nesting DIFFERENT derivatives cannot fold
  //   into a single operation.  Uniform vectors can, and must; mixed ones must be left alone.
  std::cout << "\n---- CELL 4: uniform vector fold ----\n";
  { bool const saved = FFPartial::FOLD_NESTED;
    FFPartial::FOLD_NESTED = fold;
    FFGraph G;
    FFVar x = G.add_var("x"), y = G.add_var("y"), u = G.add_var("u(x,y)"), v = G.add_var("v(x,y)");
    FFPartial OpP;
    std::vector<FFVar> const ux_vx = { OpP( u, x ), OpP( v, x ) };          // uniform: both d/dx
    std::vector<FFVar> const uni   = OpP( ux_vx, x );                       // -> d2/dx2 of {u,v}
    std::vector<FFVar> const direct= OpP( std::vector<FFVar>{ u, v }, {x,2} );
    std::vector<FFVar> const mixed_in = { OpP( u, x ), OpP( v, y ) };       // NOT uniform
    std::vector<FFVar> const mixed = OpP( mixed_in, x );
    bool const uni_same = ( uni.size() == 2 && direct.size() == 2
                         && uni[0].id() == direct[0].id() && uni[1].id() == direct[1].id() );
    auto npartial = [&]( std::vector<FFVar> const& vs ){
      auto sg = G.subgraph( vs.size(), const_cast<FFVar*>( vs.data() ) ); size_t n = 0;
      for( auto const* op : sg.l_op ) if( op && op->sameid( typeid(FFPartial) ) ) ++n; return n; };
    std::cout << "    uniform: " << ( uni_same ? "== direct" : "!= direct" )
              << "   FFPartial ops: uniform=" << npartial( uni ) << " mixed=" << npartial( mixed ) << "\n";
    expect( "C4 a UNIFORM nested vector folds to the direct form", uni_same == fold );
    expect( "C4 a MIXED nested vector is left alone (two levels survive)", npartial( mixed ) >= 3 );
    FFPartial::FOLD_NESTED = saved;
  }

  std::cout << "\n  " << g_pass << " PASS, " << g_fail << " FAIL\nOCFE_partialfold: " << ( g_fail?"FAIL":"PASS" ) << "\n";
  return g_fail ? 1 : 0;
}
