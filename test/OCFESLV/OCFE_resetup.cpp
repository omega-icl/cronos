// OCFE_resetup (2026-09-11): lifecycle regression driver for the re-setup / deep_copy segfault (Arc B).
// Walks the REUSE MATRIX on one small model.  Every transition is fenced by a banner flushed before the
// call, so the crashing transition is identified by the LAST LINE PRINTED even with no debugger.
#include <iostream>
#include <iomanip>
#include <vector>
#include "ffunc.hpp"
#include OCFE_OCFESLV_HEADER
using namespace mc;

static int g_step = 0;
static void banner( char const* s )
{ std::cout << "\n[" << std::setw(2) << ++g_step << "] " << s << std::endl; }

struct Model { FFVar t, x, u; };

// small index-1 model on one domain: dx/dt = u, x(0) = 1, u an input
static void build( FFGraph& D, OCFESLV& oc, Model& M, size_t nelem = 3, size_t nnode = 3 )
{
  M.t.set( &D ); M.x.set( &D ); M.u.set( &D );
  oc.add_domain( M.t, FFDom( 0., 1., nelem, FFDom::LGL, nnode ) );
  oc.add_state ( M.x, { M.t } );
  oc.add_input ( M.u, { M.t }, FFDom::LGR, 1, 0.5 );
  FFPartial OpP;
  oc.add_equation( OpP( M.x, M.t ) - M.u, std::vector<FFVar>{ M.t }, std::vector<int>{} );
  oc.add_equation( M.x - 1., std::vector<FFVar>{ M.t }, std::vector<int>{ FFDom::LB },
                   OCFESLV::EqnOptions( OCFESLV::EqnRole::INITIAL, 0 ) );
}

static bool solve_once( OCFESLV& oc, char const* tag )
{
  std::vector<double> xv, inp;
  if( !oc.init( xv, inp, nullptr ) ){ std::cout << "    " << tag << ": init FAILED" << std::endl; return false; }
  OCFESLV::SolveReport sr = oc.solve( xv.data(), inp.data(), nullptr );
  std::cout << "    " << tag << ": converged=" << sr.converged << " n=" << xv.size() << std::endl;
  return sr.converged;
}

int main()
{
  std::cout << "OCFESLV ** header: " << OCFESLV::HEADER_ID << std::endl;
  std::cout << "reuse matrix; the last banner printed is the crashing transition." << std::endl;

  // 1. setup -> setup again (same object, no solve between)
  { banner( "setup -> setup (same object)" );
    FFGraph D; OCFESLV oc( &D ); Model M; build( D, oc, M );
    std::cout << "    setup#1=" << oc.setup() << std::endl;
    std::cout << "    setup#2=" << oc.setup() << std::endl; }

  // 2. setup -> solve -> setup (working state populated, then re-setup)
  { banner( "setup -> solve -> setup" );
    FFGraph D; OCFESLV oc( &D ); Model M; build( D, oc, M );
    std::cout << "    setup#1=" << oc.setup() << std::endl;
    solve_once( oc, "solve#1" );
    std::cout << "    setup#2=" << oc.setup() << std::endl;
    solve_once( oc, "solve#2" ); }

  // 3. setup -> deep copy -> both solve
  { banner( "setup -> deep copy -> both solve" );
    FFGraph D; OCFESLV oc( &D ); Model M; build( D, oc, M );
    std::cout << "    setup=" << oc.setup() << std::endl;
    try{
      OCFESLV cp( oc );
      solve_once( oc, "source" );
      solve_once( cp, "copy  " );
    }
    catch( OCFESLV::Exceptions& e ){ std::cout << "    copy ctor threw ierr=" << e.ierr() << std::endl; } }

  // 4. setup -> solve -> deep copy -> copy re-setup and solve
  { banner( "setup -> solve -> deep copy -> copy setup+solve" );
    FFGraph D; OCFESLV oc( &D ); Model M; build( D, oc, M );
    oc.setup(); solve_once( oc, "source" );
    try{
      OCFESLV cp( oc );
      std::cout << "    copy setup=" << cp.setup() << std::endl;
      solve_once( cp, "copy  " );
    }
    catch( OCFESLV::Exceptions& e ){ std::cout << "    copy ctor threw ierr=" << e.ierr() << std::endl; } }

  // 5. deep copy of an environment that was NEVER set up (the SETUP_ONLY shape)
  { banner( "deep copy of a NOT-set-up environment (expects SETUP throw)" );
    FFGraph D; OCFESLV src( &D ); Model M; build( D, src, M );
    FFGraph D2; OCFESLV dst( &D2 );
    try{ dst.deep_copy_from( src ); std::cout << "    no throw" << std::endl; }
    catch( OCFESLV::Exceptions& e ){ std::cout << "    threw ierr=" << e.ierr() << std::endl; } }

  // 6. imposition change between two setups on the same object
  { banner( "setup -> change IMPOSITION_TYPE -> setup" );
    FFGraph D; OCFESLV oc( &D ); Model M; build( D, oc, M );
    oc.options.INTERFACE.IMPOSITION = OCFESLV::Options::IC_WEAK;
    std::cout << "    setup weak=" << oc.setup() << std::endl;
    solve_once( oc, "weak  " );
    oc.options.INTERFACE.IMPOSITION = OCFESLV::Options::IC_STRONG;
    std::cout << "    setup strong=" << oc.setup() << std::endl;
    solve_once( oc, "strong" ); }

  // 7. model mutated between setups (declaration changed on a set-up object)
  { banner( "setup -> add_equation -> setup" );
    FFGraph D; OCFESLV oc( &D ); Model M; build( D, oc, M );
    oc.setup(); solve_once( oc, "before" );
    oc.add_equation( M.x - 0.5, std::vector<FFVar>{ M.t }, std::vector<int>{},
                     OCFESLV::EqnOptions( OCFESLV::EqnRole::DIAGNOSTIC, 0 ) );
    std::cout << "    setup#2=" << oc.setup() << std::endl;
    solve_once( oc, "after " ); }

  // 8. set() a new DAG on a used object, then rebuild and solve
  { banner( "setup -> solve -> set(new DAG) -> rebuild -> solve" );
    FFGraph D; OCFESLV oc( &D ); Model M; build( D, oc, M );
    oc.setup(); solve_once( oc, "first " );
    FFGraph D2; oc.set( &D2 ); Model M2; build( D2, oc, M2 );
    std::cout << "    setup#2=" << oc.setup() << std::endl;
    solve_once( oc, "second" ); }

  std::cout << "\nOCFE_resetup: completed all transitions" << std::endl;
  return 0;
}
