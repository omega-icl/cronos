// ============================================================================
//  OCFE_capture_nested_refuse.cpp  --  STEP-3 refusal for nested evolution-
//  direction reductions.
//
//  An evolution reduction nested inside another evolution reduction routes the
//  inner reduction's POST-solve captured value into the outer reduction's
//  IN-solve materialised operand -- which under marching reads a window-lagged
//  value (demonstrated silently wrong by the earlier diagnostic gate:
//  N1=INT(x+INT x)=0.593 vs 1, N2=INT(x*x(1/2))=0.139 vs 1/4).  This is the
//  causality-refused pattern, so setup() must REJECT it loudly rather than
//  return a wrong number.
//
//  Expected: setup() returns FALSE with status CAPTURE_NESTED_REFUSED.
//
//  (Spatial-inside-evolution nesting is NOT refused -- there the inner reduction
//  is an in-solve aux STATE, available to the outer; the refusal keys strictly on
//  captured evolution-reduction INPUTS.  See the passing OCFE_scalar_nest* gates.)
//
//  Build: -DOCFE_OCFESLV_HEADER='"ocfeslv.hpp"'
// ============================================================================
#include <iostream>
#include <vector>

#include OCFE_OCFESLV_HEADER

using namespace mc;

int main()
{
  std::cout << "===============================================================\n"
            << "  refusal: evolution reduction nested in evolution reduction\n"
            << "===============================================================\n";

  FFGraph DAG;
  FFVar t  = DAG.add_var( "t" );
  FFVar x  = DAG.add_var( "x(t)" );

  FFPartial  OpP;
  FFIntegral OpI;

  OCFESLV oc( &DAG );
  oc.set_evolution_domain( t );
  oc.add_domain( t, FFDom( 0., 1., 3, FFDom::LGR, 4 ) );
  oc.add_state ( x,  { t } );

  FFVar dxdt = OpP( x, t );
  oc.add_equation( dxdt - 1., { t }, { FFDom::ALL - FFDom::LB },
                   OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( x, { t }, { FFDom::LB },
                   OCFESLV::EqnOptions( OCFESLV::EqnRole::INITIAL, 0 ) );

  FFVar I_in = OpI( x, t );                 // inner evolution integral
  FFVar N1   = OpI( x + I_in, t );          // OUTER evolution integral over an operand that
                                            // depends on the inner capture -> must be refused
  oc.add_output( N1 );

  oc.options.REDUCE.ORDER   = OCFESLV::Options::RED_FULL;
  oc.options.CLASSIFY.MODE       = OCFESLV::Options::CLASS_AUTO;
  oc.options.SOLVE.MARCHING = true;
  oc.options.DISPLAY_LEVEL  = 1;

  bool const ok = oc.setup();
  std::cout << "  setup() returned " << ( ok ? "true" : "false" ) << "\n"
            << "  setup_status = " << OCFESLV::setup_status_str( oc.setup_status() ) << "\n\n";

  bool pass = true;
  if( ok ){
    std::cout << "  FAIL: setup() succeeded; expected a refusal of the nested evolution reduction\n";
    pass = false;
  }
  else if( oc.setup_status() != OCFESLV::SetupStatus::CAPTURE_NESTED_REFUSED ){
    std::cout << "  FAIL: setup() failed, but not with CAPTURE_NESTED_REFUSED\n";
    pass = false;
  }
  else{
    std::cout << "  PASS: nested evolution-in-evolution reduction refused (CAPTURE_NESTED_REFUSED)\n";
  }

  std::cout << "\n  RESULT: " << ( pass ? "ALL PASS" : "FAIL" ) << "\n";
  return pass ? 0 : 1;
}
