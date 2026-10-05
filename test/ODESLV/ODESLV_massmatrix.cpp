// A nonsingular COUPLED mass matrix:  2 x0' + 1 x1' = f0 ;  1 x0' + 3 x1' = f1   (det M = 5)
// What does the classifier call it, what does dynamic_form() say, and what does the extraction do?
#include <iostream>
#include <iomanip>
#include "odeslv_base.hpp"
using namespace mc;
struct Probe : public ODESLV_BASE {
  Probe( FFGraph* dag ) : FFModel( dag ), ODESLV_BASE() {}
  bool extract(){ return _extract_from_model(); }
  std::string const& err() const { return _extractError; }
  size_t nx() const { return _nx; }
  std::vector<std::vector<FFVar>> const& de() const { return _mDE; }
  std::vector<FFVar> const& lhs() const { return _mLHS; }
  std::vector<FFVar> const& states() const { return _mX; }
  FFGraph* mdag() const { return FFModel::_dag; }
};
int main()
{
  FFGraph D; Probe M( &D );
  FFVar t = D.add_var("t"), x0 = D.add_var("x0(t)"), x1 = D.add_var("x1(t)");
  FFPartial OpP;
  int const T_INT = FFDom::ALL - FFDom::LB;
  M.add_domain( t, FFDom( 0., 1., 2, FFDom::LGR, 4 ) );
  M.add_state( x0, {t} ); M.add_state( x1, {t} );
  M.set_evolution_domain( t ); M.update_ref( x0, 1. ); M.update_ref( x1, 1. );
  auto io = FFModel::EqnOptions( FFModel::EqnRole::INTERIOR, 0 );
  auto ii = FFModel::EqnOptions( FFModel::EqnRole::INITIAL,  0 );
  M.add_equation( 2.*OpP(x0,t) + 1.*OpP(x1,t) - ( x0 + 2.*x1 ), {t}, {T_INT},     io );
  M.add_equation( 1.*OpP(x0,t) + 3.*OpP(x1,t) - ( 3.*x0 - x1 ), {t}, {T_INT},     io );
  M.add_equation( x0 - 1.,                                      {t}, {FFDom::LB}, ii );
  M.add_equation( x1 - 2.,                                      {t}, {FFDom::LB}, ii );
  M.options.DISPLAY_LEVEL = 1;
  std::cout << "setup: " << M.setup() << std::endl;
  for( auto const& [bid,cls] : M.block_classification() )
    std::cout << "  classified as " << FFModel::pde_type_name( cls.type )
              << "  At_singular=" << cls.At_singular << "  descriptor=" << cls.descriptor << std::endl;
  auto const& F = M.dynamic_form();
  std::cout << "  dynamic_form: lumped=" << F.lumped << " decoupled=" << F.decoupled
            << " differential=" << F.differential
            << " algebraic=" << F.algebraic << " index=" << F.declared_index
            << "  blocker='" << F.blocker << "'" << std::endl;
  std::cout << "  structural index in t: " << M.structural_index( 0, t ).index << std::endl;
  bool const ok = M.extract();
  std::cout << "  extraction: " << ok << ( ok? "": "  err: " + M.err() ) << std::endl;
  if( ok ){
    // does the "right-hand side" still depend on the OTHER state's derivative?  Evaluate it twice with that
    // derivative's row... simplest test: print whether the extracted RHS subgraph contains a PARTIAL operation.
    for( size_t i = 0; i < M.nx(); ++i ){
      auto sg = M.mdag()->subgraph( 1, &M.de()[0][i] );
      size_t npar = 0;
      for( auto const& op : sg.l_op ) if( op->sameid( typeid(FFPartial) ) ) ++npar;
      std::cout << "    rhs[" << i << "] holds " << npar << " derivative operation(s)"
                << ( npar? "   <= STILL CONTAINS A TIME DERIVATIVE": "" ) << std::endl;
    }
  }
  return 0;
}
