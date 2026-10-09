// Copyright (C) Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the EPL-2.0 with GPL-2.0-or-later as a Secondary License; see the LICENSE file.
// SPDX-License-Identifier: EPL-2.0 OR GPL-2.0-or-later
//
// ============================================================================
//  OCFE_ffmodel.cpp -- gate for mc::FFModel used ON ITS OWN
//
//  Includes ffmodel.hpp ONLY: no solver header, no collocation environment.  That
//  is the point of the test -- the model layer must set up, reduce order and
//  classify with no mesh and no derived class present, which is what ODESLV and
//  the claim classifier (C0) rely on.  It is built in the sandbox with ONLY
//  ffmodel.hpp, ocbase.hpp and the MC++-level headers on the include path.
//
//  Every check prints PASS or FAIL, and any FAIL makes the program exit 1.  Each
//  case also prints a DIGEST of what the model layer produced (states, equations
//  by role, LINK-row masks, classification): the sweep's byte comparison is then a
//  golden-file check of the model layer itself, independent of any solver.
//
//  1  declare a small ODE directly on an FFModel and set it up.
//  2  copy that bare declared model into a SECOND FFModel with set_model(): it
//     must set up to the same digest.
//  3  set_model() over an FFModel that already holds a DIFFERENT model REPLACES
//     it: the digest must again match case 1.
//  4  set_model() with no user DAG returns false.
//  5  heat equation u_t = u_xx: order reduction mints ONE auxiliary and ONE LINK
//     row, and the block classifies PARABOLIC -- under RED_FULL and RED_MAIN.
//  6  biharmonic chain u_t + u_xxxx = 0 with a value and two derivative BCs at
//     x=LB (OCFE_PDE12's structure): three auxiliaries, and the boundary-LINK
//     displacement retires the traces of TWO LINK rows at that face -- under
//     BOTH reductions (rev322: RED_MAIN retired only one, leaving the system
//     over-determined by one row per time node).
//  7  linear advection u_t + u_x = 0: (A) the PDE on the interior only -- the
//     automatic outgoing-characteristic closure must append exactly one row at
//     the outflow face; (B) the PDE extended to the outflow face -- the face is
//     covered, and nothing is appended.
//
//  8  the claim classifier (C0): a state differentiated in a direction by a balance row is claimed
//     NATURAL C0 there (the ODE's x in t, the heat equation's u in t, advection's u in t and x); the heat
//     equation's full report is printed as a golden digest, and a set_model() copy must report the same.
//  9  automatic differential elimination (OCFE_autoelim's one-species model): off, the block symbol is
//     rectangular and the consumer stays on z alone; on, it is square and the consumer is relocated onto
//     the source's face (r=LB) -- with the equation count unchanged.
//  10 deferred values: an integral OpI(x,t) and a point value OpEval(x,t,0.25) in OUTPUTS are captured as
//     two inputs (..._Q, ..._L), and the model PUBLISHES the contract for them through var_deferred(): one
//     record per capture, saying which input holds the value, of what expression, and whether it is the
//     integral over the evolution domain or the value at tau -- readable with no solver present, which is
//     what lets a consumer that HOLDS a model (ODESLV, DAESLV) evaluate them with its own quadrature and
//     dense output.
//  12 the mesh-free DEGREE-OF-FREEDOM BALANCE: the heat equation balances; one boundary condition too many is a
//     surplus of N_t - 1 rows, one too few a deficit of the same -- counted with no mesh, as a polynomial in the
//     domains' node counts.  (On OCFE_PDE12's marching window, N_t = 6, so N_t - 1 is the +5 surplus that its
//     RED_MAIN systems carried before rev322.)
//  13 FFModel::report(): the model as the model sees it -- domains, states, inputs, constants, equations, outputs,
//     deferred values, the classification with the evolution direction and index character, and the DOF balance.
//  11 the model PREFLIGHT: four malformed models are REJECTED by FFModel::setup() with exactly the messages
//     OCFESLV prints for them (no DAG; a state on an undeclared domain; a variable not declared; a state from
//     another DAG).  Until rev328 a bare FFModel accepted three of them silently and threw on the fourth.
//
//  A companion driver, OCFE_setmodel.cpp, checks the same copy against a full
//  OCFESLV solve; here nothing but the model layer is linked.
// ============================================================================
#include <iostream>
#include <iomanip>
#include <sstream>
#include <istream>
#include <vector>
#include <map>
#include <cmath>
#include "ffunc.hpp"
#include "ffmodel.hpp"
using namespace mc;

static int g_fail = 0, g_pass = 0;
static void check( std::string const& label, bool ok )
{
  std::cout << "    " << std::left << std::setw(64) << label << std::right << " " << ( ok? "PASS": "FAIL" ) << std::endl;
  ( ok? g_pass: g_fail )++;
}

// ---------------------------------------------------------------- digests
static char const* role_name( FFModel::EqnRole r )
{
  switch( r ){
    case FFModel::EqnRole::AUTO:       return "AUTO";
    case FFModel::EqnRole::INTERIOR:   return "INTERIOR";
    case FFModel::EqnRole::INITIAL:    return "INITIAL";
    case FFModel::EqnRole::BOUNDARY:   return "BOUNDARY";
    case FFModel::EqnRole::INTERFACE:  return "INTERFACE";
    case FFModel::EqnRole::LINK:       return "LINK";
    case FFModel::EqnRole::SURFACE:    return "SURFACE";
    case FFModel::EqnRole::DIAGNOSTIC: return "DIAGNOSTIC";
  }
  return "?";
}
static size_t count_role( FFModel const& M, FFModel::EqnRole r )
{
  size_t n = 0;
  for( auto const& e : M.var_equation() ) if( e.opt && e.opt->role == r ) ++n;
  return n;
}
static int mask_of( FFModel::t_Eqn const& e, FFVar const& dom )
{
  auto it = e.dom.find( dom );
  return it == e.dom.end() ? FFDom::ALL : it->second;
}
//! LINK rows whose mask in @p dom equals @p mask (e.g. (a,b] = ALL-LB: the trace at LB retired)
static size_t count_link_mask( FFModel const& M, FFVar const& dom, int mask )
{
  size_t n = 0;
  for( auto const& e : M.var_equation() )
    if( e.opt && e.opt->role == FFModel::EqnRole::LINK && mask_of( e, dom ) == mask ) ++n;
  return n;
}
static std::string klass( FFModel const& M )
{
  return M.is_setup() ? FFModel::pde_type_name( M.pde_type().type ) : "n/a";
}
static std::string digest( FFModel const& M, std::vector<FFVar> const& doms )
{
  std::ostringstream os;
  os << "setup=" << M.is_setup() << " states=" << M.var_state().size() << " inputs=" << M.var_input().size()
     << " eqns=" << M.var_equation().size() << " domains=" << M.var_domain().size() << " class=" << klass( M );
  std::map<std::string,size_t> byrole;
  for( auto const& e : M.var_equation() ) if( e.opt ) ++byrole[ role_name( e.opt->role ) ];
  os << " |";
  for( auto const& kv : byrole ) os << " " << kv.first << "=" << kv.second;
  if( !doms.empty() ){
    os << " | LINK masks:";
    for( auto const& e : M.var_equation() ){
      if( !e.opt || e.opt->role != FFModel::EqnRole::LINK ) continue;
      os << " [";
      for( size_t k = 0; k < doms.size(); ++k ) os << ( k? ",": "" ) << mask_of( e, doms[k] );
      os << "]";
    }
  }
  return os.str();
}
static void report( char const* tag, FFModel const& M, std::vector<FFVar> const& doms = {} )
{
  std::cout << "  " << std::left << std::setw(30) << tag << std::right << " " << digest( M, doms ) << std::endl;
}

static bool natural_c0( FFModel const& M, FFVar const& s, FFVar const& d )
{
  auto const c = M.claim_info( s, d );
  return c.found && c.tag == 'A' && c.order == 0;
}
static std::string claims_of( FFModel const& M )
{
  std::ostringstream os; M.claim_report( os ); return os.str();
}

// ---------------------------------------------------------------- 1-4: the ODE
struct Model { FFVar t, x, u; };
static void declare( FFGraph& D, FFModel& M, Model& V )
{
  V.t = D.add_var("t");  V.x = D.add_var("x(t)");  V.u = D.add_var("u(t)");
  M.add_domain( V.t, FFDom( 0., 1., 4, FFDom::LGL, 4 ) );
  M.add_state ( V.x, { V.t } );
  M.add_input ( V.u, { V.t } );
  FFPartial OpP;
  M.add_equation( OpP( V.x, V.t ) - V.u, std::vector<FFVar>{ V.t }, std::vector<int>{} );
  M.add_equation( V.x - 1., std::vector<FFVar>{ V.t }, std::vector<int>{ FFDom::LB },
                  FFModel::EqnOptions( FFModel::EqnRole::INITIAL, 0 ) );
  M.update_ref( V.x, 1.0 );
  M.update_ref( V.u, 0.5 );
}

// ---------------------------------------------------------------- 5-7: PDE blocks
int const T_INT = FFDom::ALL - FFDom::LB;                 // (a,b]: every time node but the initial one
int const X_INT = FFDom::ALL - FFDom::LB - FFDom::UB;     // (a,b): the open interior in x
int const X_NLB = FFDom::ALL - FFDom::LB;                 // (a,b]: x without its lower face
struct PDEModel { FFGraph D; FFModel M; FFVar t, x, u; PDEModel(): M( &D ) {} };
static void pde_domains( PDEModel& P )
{
  P.t = P.D.add_var("t");  P.x = P.D.add_var("x");  P.u = P.D.add_var("u(t,x)");
  P.M.add_domain( P.t, FFDom( 0., 0.5, 3, FFDom::LGR, 6 ) );
  P.M.add_domain( P.x, FFDom( 0., 1.0, 3, FFDom::LGL, 8 ) );
  P.M.add_state ( P.u, { P.t, P.x } );
  P.M.set_evolution_domain( P.t );
  P.M.update_ref( P.u, 1.0 );
  P.M.options.CLASSIFY.MODE = FFModel::Options::CLASS_AUTO;
}
static FFModel::EqnOptions opt( FFModel::EqnRole r ){ return FFModel::EqnOptions( r, 0 ); }

static void heat( PDEModel& P, FFModel::Options::ReductionType red )
{
  pde_domains( P ); FFPartial OpP;
  FFVar const& t = P.t; FFVar const& x = P.x; FFVar const& u = P.u;
  P.M.add_equation( OpP(u,t) - OpP(u,{x,2}), {t,x}, {T_INT, X_INT},         opt( FFModel::EqnRole::INTERIOR ) );
  P.M.add_equation( u - 1.,                  {t,x}, {FFDom::LB, FFDom::ALL}, opt( FFModel::EqnRole::INITIAL  ) );
  P.M.add_equation( u - 0.,                  {t,x}, {T_INT, FFDom::LB},     opt( FFModel::EqnRole::BOUNDARY ) );
  P.M.add_equation( u - 0.,                  {t,x}, {T_INT, FFDom::UB},     opt( FFModel::EqnRole::BOUNDARY ) );
  P.M.options.REDUCE.ORDER = red;
}
static void biharmonic( PDEModel& P, FFModel::Options::ReductionType red )
{
  pde_domains( P ); FFPartial OpP;
  FFVar const& t = P.t; FFVar const& x = P.x; FFVar const& u = P.u;
  P.M.add_equation( OpP(u,t) + OpP(u,{x,4}), {t,x}, {T_INT, X_INT},         opt( FFModel::EqnRole::INTERIOR ) );
  P.M.add_equation( u - 0.,                  {t,x}, {FFDom::LB, FFDom::ALL}, opt( FFModel::EqnRole::INITIAL  ) );
  P.M.add_equation( u - 0.,                  {t,x}, {T_INT, FFDom::LB},     opt( FFModel::EqnRole::BOUNDARY ) );  // value
  P.M.add_equation( OpP(u,x) - 0.,           {t,x}, {T_INT, FFDom::LB},     opt( FFModel::EqnRole::BOUNDARY ) );  // u_x
  P.M.add_equation( OpP(u,{x,2}) - 0.,       {t,x}, {T_INT, FFDom::LB},     opt( FFModel::EqnRole::BOUNDARY ) );  // u_xx
  P.M.add_equation( u - 0.,                  {t,x}, {T_INT, FFDom::UB},     opt( FFModel::EqnRole::BOUNDARY ) );
  P.M.options.REDUCE.ORDER = red;
}
static void advection( PDEModel& P, bool pde_on_outflow_face )
{
  pde_domains( P ); FFPartial OpP;
  FFVar const& t = P.t; FFVar const& x = P.x; FFVar const& u = P.u;
  P.M.add_equation( OpP(u,t) + OpP(u,x), {t,x}, {T_INT, pde_on_outflow_face? X_NLB: X_INT}, opt( FFModel::EqnRole::INTERIOR ) );
  P.M.add_equation( u - 1.,              {t,x}, {FFDom::LB, FFDom::ALL}, opt( FFModel::EqnRole::INITIAL  ) );
  P.M.add_equation( u - 1.,              {t,x}, {T_INT, FFDom::LB},     opt( FFModel::EqnRole::BOUNDARY ) );   // inflow
}

int main()
{
  std::cout << "OCFE_ffmodel ** mc::FFModel standalone (ffmodel.hpp only)" << std::endl;
  std::cout << "FFModel ** header: " << FFModel::HEADER_ID << std::endl;

  // 1 -- the model layer alone
  FFGraph D1; FFModel M1( &D1 ); Model V1; declare( D1, M1, V1 );
  std::cout << "1: an ODE on FFModel alone" << std::endl;
  check( "M1.setup() succeeds", M1.setup() );
  report( "M1", M1 );
  check( "1 state, 1 input, 2 equations, 1 domain",
         M1.var_state().size()==1 && M1.var_input().size()==1 && M1.var_equation().size()==2 && M1.var_domain().size()==1 );
  check( "classified DIFFERENTIAL_ORDINARY", klass( M1 ) == "DIFFERENTIAL_ORDINARY" );
  std::string const d1 = digest( M1, {} );

  // 2 -- copied into a second FFModel
  std::cout << "2: set_model() into a second FFModel" << std::endl;
  FFGraph D2; FFModel M2( &D2 );
  check( "set_model( M2 <- M1 ) succeeds", M2.set_model( M1 ) );
  check( "M2.setup() succeeds", M2.setup() );
  report( "M2", M2 );
  check( "M2's digest equals M1's", digest( M2, {} ) == d1 );

  // 3 -- over a DIFFERENT model: set_model replaces
  std::cout << "3: set_model() over a different model replaces it" << std::endl;
  FFGraph D3; FFModel M3( &D3 ); Model V3; declare( D3, M3, V3 );
  FFVar extra = D3.add_var("extra(t)");
  M3.add_state( extra, { V3.t } );
  M3.add_equation( extra - 2., std::vector<FFVar>{ V3.t }, std::vector<int>{} );
  check( "set_model( M3 <- M1 ) succeeds", M3.set_model( M1 ) );
  check( "M3.setup() succeeds", M3.setup() );
  report( "M3", M3 );
  check( "M3's digest equals M1's (the extra state and equation are gone)", digest( M3, {} ) == d1 );

  // 4 -- no user DAG
  std::cout << "4: set_model() with no user DAG" << std::endl;
  FFModel M4;
  check( "set_model( M4 <- M1 ) returns false", !M4.set_model( M1 ) );

  // 5 -- heat equation: one auxiliary, one LINK row, PARABOLIC, in both reductions
  std::cout << "5: heat equation u_t = u_xx" << std::endl;
  for( auto red : { FFModel::Options::RED_FULL, FFModel::Options::RED_MAIN } ){
    char const* rn = ( red == FFModel::Options::RED_FULL ? "RED_FULL" : "RED_MAIN" );
    PDEModel P; heat( P, red );
    check( std::string( rn ) + ": setup() succeeds", P.M.setup() );
    report( rn, P.M, { P.t, P.x } );
    check( std::string( rn ) + ": one auxiliary (2 states)",  P.M.var_state().size() == 2 );
    check( std::string( rn ) + ": one LINK row",              count_role( P.M, FFModel::EqnRole::LINK ) == 1 );
    check( std::string( rn ) + ": classified PARABOLIC",      klass( P.M ) == "PARABOLIC" );
  }

  // 6 -- biharmonic chain: the boundary-LINK displacement, in both reductions
  std::cout << "6: biharmonic chain u_t + u_xxxx = 0, value + u_x + u_xx at x=LB" << std::endl;
  size_t restricted[2] = { 0, 0 }; int k = 0;
  for( auto red : { FFModel::Options::RED_FULL, FFModel::Options::RED_MAIN } ){
    char const* rn = ( red == FFModel::Options::RED_FULL ? "RED_FULL" : "RED_MAIN" );
    PDEModel P; biharmonic( P, red );
    check( std::string( rn ) + ": setup() succeeds", P.M.setup() );
    report( rn, P.M, { P.t, P.x } );
    check( std::string( rn ) + ": three auxiliaries (4 states)", P.M.var_state().size() == 4 );
    restricted[k] = count_link_mask( P.M, P.x, X_NLB );
    check( std::string( rn ) + ": TWO LINK rows with their trace at x=LB retired", restricted[k] == 2 );
    ++k;
  }
  check( "RED_FULL and RED_MAIN retire the same number of LINK traces", restricted[0] == restricted[1] );

  // 7 -- linear advection: the automatic outgoing-characteristic closure appends only a deficit
  std::cout << "7: linear advection u_t + u_x = 0" << std::endl;
  size_t appended[2] = { 0, 0 };
  for( int ext = 0; ext < 2; ++ext ){
    char const* tag = ext ? "(B) PDE on the outflow face" : "(A) PDE on the interior only";
    PDEModel P; advection( P, ext );
    size_t const declared = P.M.var_equation().size();
    check( std::string( tag ) + ": setup() succeeds", P.M.setup() );
    report( tag, P.M, { P.t, P.x } );
    check( std::string( tag ) + ": classified EVOL_HYPERBOLIC", klass( P.M ) == "EVOL_HYPERBOLIC" );
    appended[ext] = P.M.var_equation().size() - declared;
  }
  check( "(A) the automatic closure appends exactly one row", appended[0] == 1 );
  check( "(B) the face is covered: nothing is appended",      appended[1] == 0 );

  // 8 -- the claim classifier
  std::cout << "8: the claim classifier (C0)" << std::endl;
  check( "ODE: x in t is NATURAL C0",              natural_c0( M1, V1.x, V1.t ) );
  { PDEModel P; advection( P, false ); P.M.setup();
    check( "advection: u in t is NATURAL C0",      natural_c0( P.M, P.u, P.t ) );
    check( "advection: u in x is NATURAL C0",      natural_c0( P.M, P.u, P.x ) ); }
  { PDEModel P; heat( P, FFModel::Options::RED_FULL ); P.M.setup();
    check( "heat (RED_FULL): u in t is NATURAL C0", natural_c0( P.M, P.u, P.t ) );
    std::string const rep = claims_of( P.M );
    std::cout << "  claim report, heat (RED_FULL):\n" << rep;
    FFGraph DC; FFModel MC( &DC );
    bool const copied = MC.set_model( P.M ) && ( MC.options.REDUCE.ORDER = FFModel::Options::RED_FULL,
                                                 MC.options.CLASSIFY.MODE = FFModel::Options::CLASS_AUTO, MC.setup() );
    check( "a set_model() copy sets up and reports the same claims", copied && claims_of( MC ) == rep ); }

  // 9 -- automatic differential elimination
  std::cout << "9: automatic differential elimination (OCFE_autoelim, one species)" << std::endl;
  size_t eqns[2] = { 0, 0 };
  for( int on = 0; on < 2; ++on ){
    FFGraph D; FFModel M( &D );
    FFVar z = D.add_var("z"), r = D.add_var("r"), T = D.add_var("T(z)"), c1 = D.add_var("c1(z)"),
          d = D.add_var("d(r,z)"), w = D.add_var("w(z)");
    FFPartial OpP;
    M.add_domain( z, FFDom( 0., 1., 2, FFDom::LGL, 4 ) );  M.add_domain( r, FFDom( 0., 1., 2, FFDom::LGL, 4 ) );
    M.add_state( T, {z} ); M.add_state( c1, {z} ); M.add_state( d, {r,z} ); M.add_state( w, {z} );
    M.update_ref( T, 2. ); M.update_ref( c1, 1. ); M.update_ref( d, .5 ); M.update_ref( w, 1. );
    int const Z_NO_LB = FFDom::ALL - FFDom::LB, R_INT = FFDom::ALL - FFDom::LB - FFDom::UB;
    M.add_equation( OpP(T,z) + OpP(c1,z) - 4.,        {z},   {Z_NO_LB},               opt( FFModel::EqnRole::INTERIOR ) );  // E
    M.add_equation( T - 2.,                           {z},   {FFDom::LB},             opt( FFModel::EqnRole::BOUNDARY ) );
    M.add_equation( OpP(c1,z) + OpP(d,r) - w - 1.,    {z,r}, {Z_NO_LB, FFDom::LB},    opt( FFModel::EqnRole::INTERIOR ) );  // G1
    M.add_equation( c1 - 1.,                          {z},   {FFDom::LB},             opt( FFModel::EqnRole::BOUNDARY ) );
    M.add_equation( w - ( 1. + z ),                   {z},   {FFDom::ALL},            opt( FFModel::EqnRole::INTERIOR ) );  // W
    M.add_equation( OpP(OpP(d,r),r) - 2.*( 1. + z ), {r,z}, {R_INT, FFDom::ALL},     opt( FFModel::EqnRole::INTERIOR ) );
    M.add_equation( d - 0.,                           {r,z}, {FFDom::LB, FFDom::ALL}, opt( FFModel::EqnRole::BOUNDARY ) );
    M.add_equation( d - 2.*( 1. + z ),                {r,z}, {FFDom::UB, FFDom::ALL}, opt( FFModel::EqnRole::BOUNDARY ) );
    M.options.REDUCE.ORDER = FFModel::Options::RED_FULL;  M.options.CLASSIFY.MODE = FFModel::Options::CLASS_AUTO;
    M.options.AUTO.DIFF_ELIM = on;
    char const* tag = on ? "DIFF_ELIM on " : "DIFF_ELIM off";
    check( std::string( tag ) + ": setup() succeeds", M.setup() );
    report( tag, M, { z, r } );
    size_t interior_no_r = 0, on_g1_face = 0;
    for( auto const& e : M.var_equation() ){
      if( !e.opt || e.opt->role != FFModel::EqnRole::INTERIOR ) continue;
      if( !e.dom.count( r ) ) ++interior_no_r;
      if( mask_of( e, z ) == Z_NO_LB && e.dom.count( r ) && mask_of( e, r ) == FFDom::LB ) ++on_g1_face;
    }
    eqns[on] = M.var_equation().size();
    if( !on ){
      check( "off: the block symbol is rectangular",                        M.pde_type().symbol_rectangular );
      check( "off: E stays on z alone (E and W carry no r)",                interior_no_r == 2 && on_g1_face == 1 );
    }
    else{
      check( "on:  the block symbol is square",                             !M.pde_type().symbol_rectangular );
      check( "on:  E is relocated onto G1's face (z=(a,b], r=LB)",          interior_no_r == 1 && on_g1_face == 2 );
    }
  }
  check( "the elimination rewrites in place: same equation count", eqns[0] == eqns[1] );

  // 10 -- deferred values
  std::cout << "10: deferred values -- OpI and OpEval over t in outputs" << std::endl;
  { FFGraph D; FFModel M( &D ); FFVar t = D.add_var("t"), x = D.add_var("x(t)");
    FFPartial OpP; FFIntegral OpI; FFEval OpEval;
    M.add_domain( t, FFDom( 0., 1., 4, FFDom::LGR, 4 ) );  M.add_state( x, {t} );  M.set_evolution_domain( t );
    M.add_equation( OpP(x,t) - 1., {t}, {FFDom::ALL - FFDom::LB}, opt( FFModel::EqnRole::INTERIOR ) );
    M.add_equation( x - 0.,        {t}, {FFDom::LB},              opt( FFModel::EqnRole::INITIAL  ) );
    M.add_output( OpI( x, t ) );  M.add_output( OpEval( x, t, 0.25 ) );  M.update_ref( x, .5 );
    size_t const in0 = M.var_input().size();
    check( "setup() succeeds", M.setup() );
    size_t q = 0, l = 0; std::string names;
    for( auto const& kv : M.var_input() ){
      std::string const n = kv.first.name(); names += " " + n;
      if( n.size() > 2 && n.compare( n.size()-2, 2, "_Q" ) == 0 ) ++q;
      if( n.size() > 2 && n.compare( n.size()-2, 2, "_L" ) == 0 ) ++l;
    }
    std::cout << "  captured inputs:" << names << std::endl;
    check( "two captured inputs, where there were none", in0 == 0 && M.var_input().size() == 2 );
    check( "one integral capture (_Q) and one point capture (_L)", q == 1 && l == 1 );
    // the published contract, read with no solver present
    auto const& dv = M.var_deferred();
    check( "var_deferred() publishes one record per capture", dv.size() == 2 );
    size_t acc = 0, pt = 0; bool inputs_match = ( dv.size() == 2 );
    for( auto const& C : dv ){
      ( C.accumulate ? acc : pt )++;
      if( !C.input.dag() || !C.source.dag() || !M.var_input().count( C.input ) ) inputs_match = false;
      std::cout << "  deferred: input=" << C.input.name() << ( C.accumulate ? "  INTEGRAL over t" : "  VALUE at tau" )
                << ( C.accumulate ? std::string() : "=" + std::to_string( C.tau ).substr(0,4) )
                << "  of " << C.source.name() << "  block=" << C.block_id << std::endl;
    }
    check( "one INTEGRAL record and one VALUE-at-tau record", acc == 1 && pt == 1 );
    check( "each record names an input of the model", inputs_match );
    for( auto const& C : dv ) if( !C.accumulate )
      check( "the point record carries tau = 0.25", std::fabs( C.tau - 0.25 ) < 1e-12 );
  }

  // 10b -- a capture whose SOURCE IS AN EXPRESSION, not a bare state: the model materialises it
  std::cout << "10b: a capture of an expression -- INT (x*x) dt" << std::endl;
  { FFGraph D; FFModel M( &D ); FFVar t = D.add_var("t"), x = D.add_var("x(t)");
    FFPartial OpP; FFIntegral OpI;
    M.add_domain( t, FFDom( 0., 1., 4, FFDom::LGR, 4 ) );  M.add_state( x, {t} );  M.set_evolution_domain( t );
    M.add_equation( OpP(x,t) - 1., {t}, {FFDom::ALL - FFDom::LB}, opt( FFModel::EqnRole::INTERIOR ) );
    M.add_equation( x - 0.,        {t}, {FFDom::LB},              opt( FFModel::EqnRole::INITIAL  ) );
    M.add_output( OpI( x * x, t ) );  M.update_ref( x, .5 );
    check( "setup() succeeds", M.setup() );
    report( "INT (x*x) dt", M, { t } );
    auto const& dv = M.var_deferred();
    check( "one deferred value", dv.size() == 1 );
    bool src_is_state = false, src_is_x = false;
    if( dv.size() == 1 ){
      src_is_state = M.var_state().count( dv[0].source ) > 0;
      src_is_x     = ( dv[0].source.id() == x.id() );
      std::cout << "  source=" << dv[0].source.name() << ( src_is_state? "  (a state of the model)": "  (an expression)" ) << std::endl;
    }
    // rev333: the model records the EXPRESSION.  Materialising it as a state is OCFESLV's evaluation strategy and
    // happens in the solver, so a bare model keeps its one state, no defining row, and stays an ODE -- which is what
    // an integrator with a quadrature right-hand side needs.
    check( "the source stays an expression (not a state of the model)", !src_is_state && !src_is_x );
    check( "no state is minted for it (x alone)",                       M.var_state().size() == 1 );
    check( "no LINK row is added",                                      count_role( M, FFModel::EqnRole::LINK ) == 0 );
    check( "the model stays DIFFERENTIAL_ORDINARY",                     klass( M ) == "DIFFERENTIAL_ORDINARY" );
    check( "the record carries the source's directions",                dv.size() == 1 && dv[0].source_dom.size() == 1 );
  }

  // 12 -- the degree-of-freedom balance
  std::cout << "12: the degree-of-freedom balance (mesh-free)" << std::endl;
  { PDEModel P; heat( P, FFModel::Options::RED_FULL );
    check( "the heat equation sets up", P.M.setup() );
    std::cout << "  heat:            rows - unknowns = " << P.M.dof_balance().str << std::endl;
    check( "the heat equation is balanced", P.M.dof_balance().balanced && P.M.dof_balance().str == "0" ); }
  { PDEModel P; heat( P, FFModel::Options::RED_FULL );
    P.M.add_equation( P.u - 0., {P.t,P.x}, {T_INT, FFDom::UB}, opt( FFModel::EqnRole::BOUNDARY ) );  // one too many
    P.M.options.DISPLAY_LEVEL = 1;                                        // the report of an unbalanced model
    check( "with a duplicated BC it still sets up", P.M.setup() );
    std::cout << "  one BC too many: rows - unknowns = " << P.M.dof_balance().str << std::endl;
    check( "a duplicated boundary condition is a surplus of N_t - 1",
           !P.M.dof_balance().balanced && P.M.dof_balance().str == "N_t - 1" ); }
  { PDEModel P; pde_domains( P ); FFPartial OpP;
    FFVar const& t = P.t; FFVar const& x = P.x; FFVar const& u = P.u;
    P.M.add_equation( OpP(u,t) - OpP(u,{x,2}), {t,x}, {T_INT, X_INT},          opt( FFModel::EqnRole::INTERIOR ) );
    P.M.add_equation( u - 1.,                  {t,x}, {FFDom::LB, FFDom::ALL}, opt( FFModel::EqnRole::INITIAL  ) );
    P.M.add_equation( u - 0.,                  {t,x}, {T_INT, FFDom::LB},      opt( FFModel::EqnRole::BOUNDARY ) );
    check( "with a missing BC it still sets up", P.M.setup() );            // no condition at x=UB
    std::cout << "  one BC missing:  rows - unknowns = " << P.M.dof_balance().str << std::endl;
    check( "a missing boundary condition is a deficit of N_t - 1",
           !P.M.dof_balance().balanced && P.M.dof_balance().str == "-N_t + 1" ); }

  // 13 -- the model's own report
  std::cout << "13: FFModel::report()" << std::endl;
  { PDEModel P; biharmonic( P, FFModel::Options::RED_MAIN );
    check( "the biharmonic sets up", P.M.setup() );
    std::ostringstream rep; P.M.report( rep );
    std::cout << rep.str();
    auto has = [&]( char const* w ){ return rep.str().find( w ) != std::string::npos; };
    check( "the report has every section",
           has("DOMAINS") && has("STATES") && has("INPUTS") && has("CONSTANTS") && has("EQUATIONS")
        && has("OUTPUTS") && has("CLASSIFICATION") && has("DEGREES OF FREEDOM") );
    check( "it names the evolution direction as set",  has("<= EVOLUTION (set)") );
    check( "it reports the balance",                   has("rows - unknowns = 0   (balanced)") );
    check( "it shows the LINK rows with their regions", has("LINK") && has("on t in") ); }   // regions since 20261001a
  { FFGraph D; FFModel M( &D ); FFVar t = D.add_var("t"), x = D.add_var("x(t)");
    FFPartial OpP; FFIntegral OpI; FFEval OpEval;
    M.add_domain( t, FFDom( 0., 1., 4, FFDom::LGR, 4 ) );  M.add_state( x, {t} );  M.set_evolution_domain( t );
    M.add_equation( OpP(x,t) - 1., {t}, {FFDom::ALL - FFDom::LB}, opt( FFModel::EqnRole::INTERIOR ) );
    M.add_equation( x - 0.,        {t}, {FFDom::LB},              opt( FFModel::EqnRole::INITIAL  ) );
    M.add_output( OpI( x, t ) );  M.add_output( OpEval( x, t, 0.25 ) );  M.update_ref( x, .5 );
    check( "the capture model sets up", M.setup() );
    std::ostringstream rep; M.report( rep );
    auto has = [&]( char const* w ){ return rep.str().find( w ) != std::string::npos; };
    // 2026-10-04: the report shows the inputs and outputs the USER declared; the deferred-value plumbing (the
    // inputs that hold the values, their sources and LINK rows) is an implementation detail, published to a
    // consumer through var_deferred() and no longer printed
    std::cout << "  (its OUTPUTS section)" << std::endl;
    { std::istringstream is( rep.str() ); std::string l; bool in = false;
      while( std::getline( is, l ) ){
        if( l.rfind("OUTPUTS",0) == 0 ) in = true;
        else if( in && l.empty() ) break;
        if( in ) std::cout << "  " << l << std::endl; } }
    bool plumbing_hidden = !has("DEFERRED VALUES") && !has("holds a deferred value");
    for( auto const& C : M.var_deferred() ) plumbing_hidden = plumbing_hidden && !has( C.input.name().c_str() );
    check( "the report lists the two outputs as declared",   has("OUTPUTS (2)") );
    check( "it hides the deferred-value plumbing",           plumbing_hidden );
    check( "the contract stays available (var_deferred)",    M.var_deferred().size() == 2 );
    // 2026-10-05: every state is listed -- also one a deferred value is taken of (a point output of a state at the
    // end of the evolution direction once hid that state from the report)
    bool all_states = true;
    for( auto const& [var,dom] : M.var_state() )
      all_states = all_states && has( ( "\n  " + var.name() + " " ).c_str() );
    check( "it lists every state of the model",              all_states ); }

  // 11 -- the model preflight: each malformed model is rejected, with OCFESLV's own message
  std::cout << "11: the model preflight -- malformed models are rejected" << std::endl;
  auto rejects = [&]( char const* what, char const* msg, auto f ){
    std::ostringstream err; auto* old = std::cerr.rdbuf( err.rdbuf() );
    bool ok = false; std::string res;
    try{ bool const s = f(); ok = !s && err.str().find( msg ) != std::string::npos;
         res = s ? "setup() returned TRUE" : ( "setup() returned false: " + err.str().substr( 0, err.str().find('\n') ) ); }
    catch( std::exception const& e ){ res = std::string( "setup() threw: " ) + e.what(); }
    std::cerr.rdbuf( old );
    std::cout << "  " << what << " -> " << res << std::endl;
    check( std::string( what ) + ": rejected with \"" + msg + "\"", ok );
  };
  rejects( "a state on an undeclared domain", "**ERROR: DOMAIN t MISSING FOR STATE x(t)", [](){
       FFGraph D; FFModel M( &D ); FFVar t = D.add_var("t"), x = D.add_var("x(t)"); FFPartial P;
       M.add_state( x, {t} ); M.add_equation( P(x,t) - 1., std::vector<FFVar>{t}, std::vector<int>{} ); return M.setup(); } );
  rejects( "an equation with an undeclared variable", "**ERROR: VARIABLE y MISSING", [](){
       FFGraph D; FFModel M( &D ); FFVar t = D.add_var("t"), x = D.add_var("x(t)"), y = D.add_var("y"); FFPartial P;
       M.add_domain( t, FFDom( 0., 1., 2, FFDom::LGL, 3 ) ); M.add_state( x, {t} );
       M.add_equation( P(x,t) - y, std::vector<FFVar>{t}, std::vector<int>{} ); return M.setup(); } );
  rejects( "a state from another DAG", "**ERROR: STATE x(t) DOES NOT BELONG TO THE MODEL DAG", [](){
       FFGraph D, E; FFModel M( &D ); FFVar t = D.add_var("t"), x = E.add_var("x(t)");
       M.add_domain( t, FFDom( 0., 1., 2, FFDom::LGL, 3 ) ); M.add_state( x, {t} ); return M.setup(); } );
  rejects( "no DAG at all", "**ERROR: UNDEFINED DAG", [](){ FFModel M; return M.setup(); } );

  std::cout << "\n  OCFE_ffmodel: " << g_pass << " passed, " << g_fail << " failed" << std::endl;
  return g_fail ? 1 : 0;
}
