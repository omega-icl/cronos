// does the EqnRole overload resolve, and give the same model as EqnOptions?
#include <iostream>
#include "daeslv_idas.hpp"
using namespace mc;
int main(){
  FFPartial OpP; FFEval OpE; int const TI = FFDom::ALL - FFDom::LB;
  double f[2];
  for( int v = 0; v < 2; ++v ){
    FFGraph G; DAESLV_IDAS M; M.set_dag( &G );
    FFVar t = G.add_var("t"), x = G.add_var("x(t)"), k = G.add_var("k");
    M.add_domain( t, FFDom( 0., 1., 2, FFDom::LGR, 3 ) );
    M.add_state( x, {t} ); M.add_input( k ); M.set_evolution_domain( t ); M.update_ref( x, 1. );
    if( v == 0 ){   // the EqnOptions form
      M.add_equation( OpP(x,t) + k*x, {t}, {TI}, FFModel::EqnOptions( FFModel::EqnRole::INTERIOR, 0 ) );
      M.add_equation( x - 1.,         {t}, {FFDom::LB}, FFModel::EqnOptions( FFModel::EqnRole::INITIAL, 0 ) );
    } else {        // the EqnRole overload
      M.add_equation( OpP(x,t) + k*x, {t}, {TI}, FFModel::EqnRole::INTERIOR );
      M.add_equation( x - 1.,         {t}, {FFDom::LB}, FFModel::EqnRole::INITIAL );
    }
    M.add_output( OpE( x, t, 1. ) );
    M.options.SOLVE.RTOL = M.options.SOLVE.ATOL = 1e-10; M.options.DISPLAY_LEVEL = 0;
    if( !M.setup() ){ std::cout << "setup failed: " << M.extract_error() << "\n"; return 1; }
    std::vector<double> P{0.8}; M.solve( P ); f[v] = M.val_function()[0];
  }
  bool const same = ( f[0] == f[1] );
  std::cout << "    EqnOptions: " << f[0] << "   EqnRole: " << f[1] << "\n";
  std::cout << "    " << ( same? "PASS": "FAIL" ) << "  the EqnRole overload declares the SAME model\n";
  std::cout << "\n==========  DAESLV_eqnrole: " << ( same? 1: 0 ) << " passed, " << ( same? 0: 1 )
            << " failed -- " << ( same? "ALL PASS": "FAILURES" ) << "  ==========\n\n";
  return same? 0: 1; }
