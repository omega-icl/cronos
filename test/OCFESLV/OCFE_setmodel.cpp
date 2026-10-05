// Copyright (C) Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// ============================================================================
//  OCFE_setmodel.cpp -- gate for FFModel::set_model() through a full OCFESLV solve
//
//  set_model() copies the BARE DECLARED model -- constants, domains, inputs with
//  their declared collocation and continuity, states, equations, outputs, the
//  evolution domain and the classification references -- and resets everything
//  derived.  This driver checks that a model reached that way is indistinguishable
//  from one declared directly, all the way through setup and the solve:
//
//    A  declare directly, set up, solve                          -- the reference
//    B  declare in one environment, set_model into a SECOND one
//       (its own user DAG), set up, solve                        -- must equal A
//    C  set_model into an environment that already held a
//       DIFFERENT model                                          -- must equal A
//    D  set_model with no user DAG set                           -- must return false
//    E  the SOURCE environment, after being copied from          -- must equal A
//
//  The companion driver OCFE_ffmodel.cpp exercises the same call with the model
//  layer alone (ffmodel.hpp only, no solver).
// ============================================================================
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <string>
#include <utility>   // std::pair -- per-stage cap overrides

#ifndef OCFE_OCFESLV_HEADER
#define OCFE_OCFESLV_HEADER "ocfeslv.hpp"
#endif
#include OCFE_OCFESLV_HEADER
using namespace mc;

static int g_pass=0, g_fail=0;
static void check( char const* nm, bool ok )
{
  (ok?g_pass:g_fail)++;
  std::cout<<"  "<<std::left<<std::setw(52)<<nm<<std::right<<" "<<(ok?"PASS":"FAIL")<<"\n";
}

static double const gDa0 = 2.0e2;      // base rate; kap ramps it geometrically to 2e5

struct Vars { FFVar z,u,w,lam,mu,kap; };

static void build( FFGraph& DAG, OCFESLV& oc, Vars& V )
{
  V.z  = DAG.add_var("z");
  V.u  = DAG.add_var("u(z)");
  V.w  = DAG.add_var("w(z)");
  V.lam= DAG.add_var("lam");
  V.mu = DAG.add_var("mu");
  V.kap= DAG.add_var("kap");
  FFPartial OpP;

  oc.add_domain( V.z, FFDom( 0., 1., 4, FFDom::CGL, 5 ) );
  oc.add_state ( V.u, { V.z } );
  oc.add_state ( V.w, { V.z } );
  oc.add_input ( V.lam, 0.0, true );
  oc.add_input ( V.mu , 0.0, true );
  oc.add_input ( V.kap, 0.0, true );

  // kap in [0,1] maps geometrically onto Da in [Da0, 1000*Da0]
  FFVar Da  = gDa0 * pow( 1000.0, V.kap );
  FFVar PDE = OpP(OpP(V.u,V.z),V.z) - Da*V.u*V.w + V.lam*( 1.0 + V.z );
  FFVar ALG = V.w - ( 1.0 + V.mu*V.z );

  int const Z_INT = FFDom::ALL - FFDom::LB - FFDom::UB;
  oc.add_equation( PDE,      { V.z }, { Z_INT     }, OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( ALG,      { V.z }, { FFDom::ALL}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( V.u,      { V.z }, { FFDom::LB }, OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );
  oc.add_equation( V.u-1.0,  { V.z }, { FFDom::UB }, OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );

  oc.update_ref( V.u, [](OCFESLV::t_Coord const&){ return 0.0; } );
  oc.update_ref( V.w, [](OCFESLV::t_Coord const&){ return 1.0; } );

  oc.options.REDUCE.ORDER  = OCFESLV::Options::RED_FULL;
  oc.options.CLASSIFY.MODE      = OCFESLV::Options::CLASS_AUTO;
  oc.options.DISPLAY_LEVEL = 0;
}


// ============================================================================================
// OCFE_setmodel (sandbox test, 2026-09-13): FFModel::set_model() copies the BARE DECLARED model.
//   A  declare directly, set up, solve                        -> reference
//   B  declare in one env, set_model into a SECOND env (own user DAG), set up, solve -> must equal A
//   C  set_model into an env that already held a DIFFERENT model  -> must equal A (it replaces)
//   D  set_model with no user DAG set                          -> must return false
//   E  the SOURCE is untouched by being copied                 -> still equals A
// ============================================================================================
struct Run { bool ok=false; size_t nvar=0, neqn=0; double r=0.; };

static Run solve_it( OCFESLV& oc )
{
  Run out;
  if( !oc.setup() ) return out;
  out.nvar = oc.var_state().size(); out.neqn = oc.var_equation().size();
  std::vector<double> xv, inp;
  if( !oc.init( xv, inp, nullptr ) ) return out;
  OCFESLV::SolveReport sr = oc.solve( xv.data(), inp.data(), nullptr );
  out.ok = sr.converged;
  for( double v : xv ) out.r = std::max( out.r, std::abs( v ) );
  return out;
}

static bool same( Run const& r, Run const& ref )
{
  return r.ok==ref.ok && r.nvar==ref.nvar && r.neqn==ref.neqn && r.r==ref.r;
}

static void line( char const* tag, Run const& r, Run const* ref )
{
  std::cout << "  " << std::left << std::setw(34) << tag << std::right
            << " conv=" << r.ok << " states=" << r.nvar << " eqns=" << r.neqn
            << " max|x|=" << std::setprecision(12) << r.r;
  if( ref ) std::cout << "  equals A: " << ( r.ok==ref->ok && r.nvar==ref->nvar && r.neqn==ref->neqn
                                             && r.r==ref->r );
  std::cout << std::endl;
}

int main()
{
  std::cout << "OCFESLV ** header: " << OCFESLV::HEADER_ID << std::endl;
  FFGraph DA; OCFESLV A( &DA ); Vars VA; build( DA, A, VA );
  Run rA = solve_it( A );
  line( "A: declared directly", rA, nullptr );
  check( "A: declared directly -- converges", rA.ok );

  FFGraph DS; OCFESLV S( &DS ); Vars VS; build( DS, S, VS );      // the source, never set up
  FFGraph DB; OCFESLV B( &DB );
  bool const okB = B.set_model( S );
  std::cout << "  set_model(B <- S): " << okB << std::endl;
  check( "set_model(B <- S) succeeds", okB );
  Run rB = solve_it( B );
  line( "B: copied before setup", rB, &rA );
  check( "B: copied before setup -- converges and equals A", same( rB, rA ) );

  FFGraph DC; OCFESLV C( &DC ); Vars VC; build( DC, C, VC );
  FFVar extra = DC.add_var("extra(z)");                          // a DIFFERENT model first
  C.add_state( extra, { VC.z } );
  C.add_equation( extra - 1., std::vector<FFVar>{ VC.z }, std::vector<int>{} );
  bool const okC = C.set_model( S );
  std::cout << "  set_model(C <- S) over a different model: " << okC << std::endl;
  check( "set_model(C <- S) over a different model succeeds", okC );
  Run rC = solve_it( C );
  line( "C: replaced an existing model", rC, &rA );
  check( "C: replaced an existing model -- equals A", same( rC, rA ) );

  OCFESLV Dn;                                                      // no user DAG
  bool const okD = Dn.set_model( S );
  std::cout << "  set_model(D <- S) with no user DAG: " << okD << " (expect 0)" << std::endl;
  check( "set_model(D <- S) with no user DAG returns false", !okD );

  Run rS = solve_it( S );
  line( "E: the source, after copying", rS, &rA );
  check( "E: the source, solved after being copied -- equals A", same( rS, rA ) );
  std::cout << "\n  OCFE_setmodel: " << g_pass << " passed, " << g_fail << " failed" << std::endl;
  return g_fail ? 1 : 0;
}
