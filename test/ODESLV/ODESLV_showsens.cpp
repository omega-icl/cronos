#include "odeslvs_cvodes.hpp"
#include <iostream>
int main(){
  mc::FFGraph DAG; mc::ODESLVS_CVODES IVP( &DAG );
  mc::FFVar t=DAG.add_var("t"), X0=DAG.add_var("x0(t)"), X1=DAG.add_var("x1(t)"), TF=DAG.add_var("TF"), X10=DAG.add_var("x10"), U=DAG.add_var("u(t)");
  mc::FFPartial OpP; mc::FFEval OpE;
  IVP.add_domain( t, mc::FFDom( 0., 1., 3, mc::FFDom::LGR, 4 ) );
  IVP.add_state( X0, {t} ); IVP.add_state( X1, {t} );
  IVP.add_input( TF ); IVP.add_input( X10 ); IVP.add_input( U, {t}, mc::FFDom::LGR, 1 );
  IVP.set_evolution_domain( t );
  int const T_INT = mc::FFDom::ALL - mc::FFDom::LB;
  auto io = mc::FFModel::EqnOptions( mc::FFModel::EqnRole::INTERIOR, 0 ), ii = mc::FFModel::EqnOptions( mc::FFModel::EqnRole::INITIAL, 0 );
  IVP.add_equation( OpP(X0,t) - TF*X1, {t}, {T_INT}, io );
  IVP.add_equation( OpP(X1,t) - TF*(U*X0 - 2.*X1)*t, {t}, {T_INT}, io );
  IVP.add_equation( X0 - 0., {t}, {mc::FFDom::LB}, ii );
  IVP.add_equation( X1 - X10, {t}, {mc::FFDom::LB}, ii );
  IVP.add_output( OpE( X0 - 1., t, 1. ) ); IVP.add_output( OpE( X1, t, 1. ) );
  mc::ODESLVS_CVODES SENS( &DAG ); mc::FFVar du; std::vector<mc::FFVar> vS; std::string err;
  IVP.fdiff( SENS, U, du, vS, err );
  auto show=[&](char const* tag, mc::FFVar const& f){ mc::FFSubgraph sg=DAG.subgraph(1,&f); std::cout<<"    "<<tag<<"  "<<mc::FFExpr::subgraph(&DAG,sg)[0]<<"\n"; };
  std::cout<<"  states : "; for( auto const& [x,d] : SENS.var_state() ) std::cout<<x.name()<<" "; std::cout<<"\n";
  std::cout<<"  inputs : "; for( auto const& [w,d] : SENS.var_declared_input() ) std::cout<<w.name()<<(d.empty()?"":"(t)")<<" "; std::cout<<"\n";
  std::cout<<"  equations:\n"; for( auto const& e : SENS.var_equation() ) show("E", e.var);
  std::cout<<"  outputs:\n";   for( auto const& f : SENS.var_output() ) show("F", f.var);
  SENS.setup(); std::cout<<"  parameters after extraction (np="<<SENS.np()<<"): "; for( auto const& p : SENS.var_parameter() ) std::cout<<p.name()<<" "; std::cout<<"\n";
  return 0; }
