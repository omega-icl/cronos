// OCFE_faces.cpp -- conditions per face, counted against the characteristics (Martinson & Barton, items 5-7)
// =========================================================================================================
// A 2x2 wave system,  u1_t + c u2_x = 0 ,  u2_t + c u1_x = 0 , has characteristic speeds +-c: ONE characteristic
// enters at each end, so a well-posed problem takes ONE condition at x=0 and ONE at x=L.  Their telegrapher
// example fails exactly by putting both at the same end -- the substation -- and no count notices: the system
// stays SQUARE, so the DOF balance is silent, and the hyperbolic closure only ever FILLS a face that is short.
// rev340 keeps the per-face account, so the excess is visible: more condition rows at a face than characteristics
// leaving it means the extra one pins an INCOMING direction, which is data at the wrong end.
//   A  one condition at each end   -> sets up; the model-level face account records outgoing/covered/rows
//   B  both conditions at x=LB     -> REFUSED by OCFESLV's existing hyperbolic BC guard ("needs 1 incoming
//      condition but 2 BOUNDARY rows are collocated there"), while the DOF balance stays silent: the system is
//      square, so only a characteristic count can see it.  rev340 puts the same count in the MODEL, where a
//      consumer that is not OCFESLV can read it.
// =========================================================================================================
#include <iostream>
#include <iomanip>
#include <vector>
#ifndef OCFE_OCFESLV_HEADER
#define OCFE_OCFESLV_HEADER "ocfeslv.hpp"
#endif
#include OCFE_OCFESLV_HEADER
using namespace mc;

static int g_pass = 0, g_fail = 0;
static void check( char const* nm, bool ok )
{ std::cout << "  " << std::left << std::setw(56) << nm << std::right << ( ok? " PASS": " FAIL" ) << std::endl;
  ( ok? g_pass: g_fail )++; }

int main()
{
  std::cout << "OCFESLV ** header: " << OCFESLV::HEADER_ID << "\n";
  int const T_INT = FFDom::ALL - FFDom::LB;
  int const X_INT = FFDom::ALL - FFDom::LB - FFDom::UB;
  double const c = 1.0;

  for( int variant = 0; variant < 2; ++variant ){
    std::cout << "\n---- " << ( variant? "B: BOTH conditions at x=LB (the telegrapher mistake)"
                                       : "A: one condition at each end" ) << " ----\n";
    FFGraph D; OCFESLV oc( &D );
    FFVar t = D.add_var("t"), x = D.add_var("x");
    FFVar u1 = D.add_var("u1(t,x)"), u2 = D.add_var("u2(t,x)");
    FFPartial OpP;
    oc.add_domain( t, FFDom( 0., 0.2, 2, FFDom::LGR, 4 ) );
    oc.add_domain( x, FFDom( 0., 1.0, 2, FFDom::LGL, 5 ) );
    oc.add_state( u1, {t,x} ); oc.add_state( u2, {t,x} );
    oc.set_evolution_domain( t );
    oc.update_ref( u1, 0. ); oc.update_ref( u2, 0. );
    OCFESLV::EqnOptions io( OCFESLV::EqnRole::INTERIOR, 0 ), ii( OCFESLV::EqnRole::INITIAL, 0 ),
                        ib( OCFESLV::EqnRole::BOUNDARY, 0 );
    oc.add_equation( OpP(u1,t) + c*OpP(u2,x), {t,x}, {T_INT, X_INT},          io );
    oc.add_equation( OpP(u2,t) + c*OpP(u1,x), {t,x}, {T_INT, X_INT},          io );
    oc.add_equation( u1 - 0.,                 {t,x}, {FFDom::LB,FFDom::ALL},  ii );
    oc.add_equation( u2 - 0.,                 {t,x}, {FFDom::LB,FFDom::ALL},  ii );
    if( !variant ){
      oc.add_equation( u1 + u2, {t,x}, {T_INT, FFDom::LB}, ib );      // the characteristic entering at LB
      oc.add_equation( u1 - u2, {t,x}, {T_INT, FFDom::UB}, ib );      // the one entering at UB
    }
    else{
      oc.add_equation( u1 + u2, {t,x}, {T_INT, FFDom::LB}, ib );      // both at the same end
      oc.add_equation( u1 - u2, {t,x}, {T_INT, FFDom::LB}, ib );
    }
    oc.options.REDUCE.ORDER      = OCFESLV::Options::RED_FULL;
    oc.options.DISPLAY_LEVEL     = 1;
    oc.options.FATAL.REDUCED_DOF = false;
    bool const ok = oc.setup();
    std::cout << "  setup=" << ok << "  balance: rows-unknowns=" << oc.dof_balance().str << "\n";
    size_t excess = 0, faces = 0;
    for( auto const& fc : oc.face_conditions() ){
      ++faces;
      std::cout << "  block " << fc.block_id << "  " << fc.direction << " "
                << ( fc.face == FFDom::LB? "LB": "UB" ) << ": leaving=" << fc.outgoing
                << " (covered=" << fc.covered << " appended=" << fc.appended << ")"
                << " entering=" << fc.incoming << " rows_here=" << fc.rows_at_face
                << ( fc.rows_at_face > fc.incoming? "   <= EXCESS": "" ) << "\n";
      if( fc.rows_at_face > fc.incoming ) ++excess;
    }
    if( !variant ){
      check( "A: the model sets up and is balanced", ok && oc.dof_balance().balanced );
      check( "A: a face account was produced",       faces > 0 );
      check( "A: no face carries an excess",         excess == 0 );
      bool one_each = !oc.face_conditions().empty();
      for( auto const& fc : oc.face_conditions() )
        if( !( fc.incoming == 1 && fc.outgoing == 1 && fc.rows_at_face == 1 ) ) one_each = false;
      check( "A: each face has 1 entering, 1 leaving, 1 condition row", one_each );
    }
    else{
      // MEASURED: OCFESLV already carries a HYPERBOLIC BC GUARD that counts incoming characteristics per face and
      // REFUSES the model ("face z=LB needs 1 incoming condition(s) but 2 BOUNDARY row(s) are collocated there").
      // So the telegrapher mistake is caught today -- fatally, by the solver -- and the count stays silent because
      // the system is square.  This gates that behaviour.
      check( "B: the solver's hyperbolic BC guard REFUSES it", !ok );
      check( "B: while the DOF balance stays silent (square)", oc.dof_balance().str == "0" );
    }
  }
  std::cout << "\n  OCFE_faces: " << g_pass << " passed, " << g_fail << " failed\n";
  return g_fail? 1: 0;
}
