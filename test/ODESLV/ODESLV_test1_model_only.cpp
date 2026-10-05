// The MODEL half of test1_ff.cpp, compiled without SUNDIALS: does the description go in and come out right?
#include <iostream>
#include <iomanip>
#include "odeslv_base.hpp"
using namespace mc;
struct Probe : public ODESLV_BASE {
  Probe( FFGraph* dag ) : FFModel( dag ), ODESLV_BASE() {}
  bool solver_setup(){ if( !FFModel::is_setup() && !FFModel::setup() ) return false; return _SETUP(); }
  size_t nx() const { return _nx; } size_t np() const { return _np; }
  size_t nq() const { return _nq; } size_t nf() const { return _nf; } size_t ns() const { return _ns; }
  std::vector<double> const& dT() const { return _dT; }
  std::vector<std::map<size_t,FFVar>> const& fct() const { return _mFCT; }
};
int main()
{
  FFGraph D; Probe IVP( &D );
  double const t0 = 0., tf = 10.; size_t const NS = 4;
  FFVar t = D.add_var("t"), X0 = D.add_var("x0(t)"), X1 = D.add_var("x1(t)"), P0 = D.add_var("p0"), P1 = D.add_var("p1");
  FFPartial OpP; FFIntegral OpI; FFEval OpEval;
  IVP.add_domain( t, FFDom( t0, tf, NS, FFDom::LGR, 4 ) );
  IVP.add_state( X0, {t} ); IVP.add_state( X1, {t} ); IVP.add_input( P0 ); IVP.add_input( P1 );
  IVP.set_evolution_domain( t ); IVP.update_ref( X0, 1.2 ); IVP.update_ref( X1, 1.1 );
  int const T_INT = FFDom::ALL - FFDom::LB;
  auto const io = FFModel::EqnOptions( FFModel::EqnRole::INTERIOR, 0 );
  auto const ii = FFModel::EqnOptions( FFModel::EqnRole::INITIAL,  0 );
  IVP.add_equation( OpP(X0,t) - P0*X0*( 1. - X1 ), {t}, {T_INT},     io );
  IVP.add_equation( OpP(X1,t) - P0*X1*( X0 - 1. ), {t}, {T_INT},     io );
  IVP.add_equation( X0 - 1.2,                      {t}, {FFDom::LB}, ii );
  IVP.add_equation( X1 - ( 1.1 + 0.01*P1 ),        {t}, {FFDom::LB}, ii );
  IVP.add_output( OpI( X1, t ) );
  IVP.add_output( OpEval( X0 * X1, t, tf ) );
  IVP.add_output( OpEval( X0, t, 3.7 ) );
  bool const ok = IVP.solver_setup();
  std::cout << "setup (model + local copy): " << ok << ( ok? "": "  err: " + IVP.extract_error() ) << std::endl;
  std::cout << "  nx=" << IVP.nx() << " np=" << IVP.np() << " nq=" << IVP.nq()
            << " nf=" << IVP.nf() << " ns=" << IVP.ns() << std::endl;
  std::cout << "  stages:"; for( double d : IVP.dT() ) std::cout << " " << d; std::cout << std::endl;
  for( size_t k = 0; k < IVP.fct().size(); ++k )
    for( auto const& [i,f] : IVP.fct()[k] )
      std::cout << "  function " << i << " evaluated at t=" << IVP.dT()[k] << std::endl;
  return ok? 0: 1;
}
