// test2's model through FFModel, model half only: what does the extraction make of a PIECEWISE-CONSTANT input?
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
  std::vector<FFVar> const& pars() const { return _mP; }
  std::map<std::string,std::vector<FFVar>> const& levels() const { return _mInpLevels; }
  std::vector<std::vector<FFVar>> const& de() const { return _mDE; }
  FFGraph* mdag() const { return FFModel::_dag; }
};
int main()
{
  FFGraph D; Probe IVP( &D );
  size_t const NS = 3; double const t0 = 0., tf = 1.;
  FFVar t  = D.add_var("t");
  FFVar X0 = D.add_var("x0(t)"), X1 = D.add_var("x1(t)");
  FFVar TF = D.add_var("TF"), X10 = D.add_var("x10"), U = D.add_var("u(t)");
  FFPartial OpP; FFEval OpEval;
  int const T_INT = FFDom::ALL - FFDom::LB;

  IVP.add_domain( t, FFDom( t0, tf, NS, FFDom::LGR, 4 ) );      // NS stages = NS elements
  IVP.add_state( X0, {t} ); IVP.add_state( X1, {t} );
  IVP.add_input( TF );  IVP.add_input( X10 );                   // time-invariant
  IVP.add_input( U, { t }, FFDom::LGR, 1 );                     // piecewise constant: one level per stage
  IVP.set_evolution_domain( t );
  IVP.update_ref( X0, 0. ); IVP.update_ref( X1, 0.5 );
  auto const io = FFModel::EqnOptions( FFModel::EqnRole::INTERIOR, 0 );
  auto const ii = FFModel::EqnOptions( FFModel::EqnRole::INITIAL,  0 );
  IVP.add_equation( OpP(X0,t) - TF * X1,                       {t}, {T_INT},     io );
  IVP.add_equation( OpP(X1,t) - TF * ( U*X0 - 2.*X1 ) * t,     {t}, {T_INT},     io );
  IVP.add_equation( X0 - 0.,                                   {t}, {FFDom::LB}, ii );
  IVP.add_equation( X1 - X10,                                  {t}, {FFDom::LB}, ii );
  IVP.add_output( OpEval( X0 - 1., t, tf ) );
  IVP.add_output( OpEval( X1,      t, tf ) );

  std::cout << "model setup: " << IVP.setup() << std::endl;
  std::cout << "  var_input() after setup:";
  for( auto const& [var,dom] : IVP.var_input() )
    std::cout << " " << var.name() << ( IVP.is_deferred_input( var )? "[deferred]": dom.empty()? "(scalar)": "(on t)" );
  std::cout << "\n  var_declared_input():";
  for( auto const& [var,dom] : IVP.var_declared_input() )
    std::cout << " " << var.name() << ( dom.empty()? "(scalar)": "(on t)" );
  std::cout << std::endl;
  for( auto const& [dv,dd] : IVP.var_declared_input() )
    for( auto const& [wv,wd] : IVP.var_input() )
      if( wv.name() == dv.name() )
        std::cout << "    " << dv.name() << ": declared id=" << dv.id().second
                  << " working id=" << wv.id().second
                  << ( dv.id() == wv.id()? "  SAME NODE": "  DIFFERENT NODE" ) << std::endl;
  bool const ok = IVP.solver_setup();
  std::cout << "solver _SETUP: " << ok << ( ok? "": "  err: " + IVP.extract_error() ) << std::endl;
  if( ok ){
    std::cout << "  nx=" << IVP.nx() << " np=" << IVP.np() << " nq=" << IVP.nq()
              << " nf=" << IVP.nf() << " ns=" << IVP.ns() << std::endl;
    std::cout << "  parameters:"; for( auto const& p : IVP.pars() ) std::cout << " " << p.name();
    std::cout << std::endl;
    std::cout << "  stage-wise right-hand sides: " << IVP.de().size() << " entr(y/ies)" << std::endl;
    for( auto const& [n,lv] : IVP.levels() ) std::cout << "  " << n << " -> " << lv.size() << " level parameters\n";
    // does the extracted RHS still hold the time-varying input?
    for( size_t k = 0; k < IVP.de().size(); ++k )
      for( size_t i = 0; i < IVP.de()[k].size(); ++i ){
        auto sg = IVP.mdag()->subgraph( 1, &IVP.de()[k][i] );
        bool holds_u = false;
        for( auto const& op : sg.l_op )
          for( size_t jj = 0; jj < op->varin.size(); ++jj )
            if( op->varin[jj] && op->varin[jj]->name() == "u(t)" ) holds_u = true;
        if( holds_u ) std::cout << "    rhs[" << k << "][" << i << "] STILL HOLDS the input u(t)" << std::endl;
      }
  }
  return 0;
}
