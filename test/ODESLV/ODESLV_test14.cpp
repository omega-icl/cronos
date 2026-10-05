// test3_ff.cpp -- GATE for selective sensitivity analysis in ODESLV (step 2).
// =============================================================================================
// test2's model, solved repeatedly while REGISTERING DIFFERENT SUBSETS of its inputs as controls.
// ODESLV today differentiates w.r.t. every declared input at once; the step-2 change makes it honour
// FFModel's control registry, so the sensitivity runs over the registered controls only, in the
// CANONICAL CONTROL-VECTOR ORDER (FFVar order, offsets ascending).
//
// The model has np = 5 parameters, in _mP order:
//        0: TF          (time-invariant)
//        1: x10         (time-invariant)
//        2,3,4: u(t)[0..2]   (piecewise constant, one level per stage -- ONE control, 3 DOFs)
//
// WHY THIS DRIVER EXISTS.  test1 (np=2) and test2 (np=5, everything registered) cannot exercise
// _ndxSEN != identity, so every site mis-classified between "parameter vector" and "sensitivity
// direction" would still produce a plausible answer there.  Here a wrong classification produces a
// column count or a column VALUE that disagrees with the all-parameters run, which is checkable.
//
// BOTH CONTROLS.  Case A registers nothing and is the reference (and must stay byte-comparable to
// test2_ff).  Cases B/C/D register strict subsets.  BEFORE step 2 they must FAIL -- ODESLV ignores
// the registry, so each returns 5 directions where n_control_dof() says 3, 4 and 1.  AFTER step 2
// they must PASS with the selected columns equal to the reference columns.  Case C is the
// discriminating one: {TF, u} is NON-CONTIGUOUS in _mP (it skips x10), so a change that merely
// truncated the parameter list rather than indexing through _ndxSEN would pass B and D and fail C.
//
// The last block feeds the gradient through FFModel::control_columns() -- once _Dfp is
// nf x n_control_dof in control-vector order, the same call that slices OCFESLV::sens_jacobian()
// must slice ODESLV's function gradient.  That is the whole point of putting the layout in FFModel.

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <cmath>

#include "ffode.hpp"
#include <algorithm>
#include <functional>

size_t const NS = 3;
double const t0 = 0., tf = 1.;
std::vector<double> const P0( { 6., 0.5, 0.5, 0.5, 0.5 } );   // TF, x10, then one level per stage

enum Reg { REG_NONE = 0, REG_TF = 1, REG_X10 = 2, REG_U = 4 };

//! @brief Build test2's model into a fresh solver, register the requested controls, and return the
//! function gradients from BOTH sensitivity routes.  A fresh FFGraph per run, so no re-setup
//! semantics are involved and the runs cannot contaminate one another.
struct Run
{
  size_t ncd = 0, np = 0, nf = 0;
  std::vector<std::vector<double>> fsa, asa;   // [direction][function]
  bool ok = false;
};

Run build_and_solve( int reg, bool verbose = false )
{
  Run R;
  mc::FFGraph DAG;
  mc::ODESLVS_CVODES IVP( &DAG );

  mc::FFVar t   = DAG.add_var( "t" );
  mc::FFVar X0  = DAG.add_var( "x0(t)" ), X1 = DAG.add_var( "x1(t)" );
  mc::FFVar TF  = DAG.add_var( "TF" ), X10 = DAG.add_var( "x10" );
  mc::FFVar U   = DAG.add_var( "u(t)" );
  mc::FFPartial OpP;  mc::FFEval OpEval;

  IVP.add_domain( t, mc::FFDom( t0, tf, NS, mc::FFDom::LGR, 4 ) );
  IVP.add_state( X0, {t} );  IVP.add_state( X1, {t} );
  IVP.add_input( TF );  IVP.add_input( X10 );
  IVP.add_input( U, { t }, mc::FFDom::LGR, 1 );
  IVP.set_evolution_domain( t );
  IVP.update_ref( X0, 0. );  IVP.update_ref( X1, 0.5 );

  int const T_INT = mc::FFDom::ALL - mc::FFDom::LB;
  auto const io = mc::FFModel::EqnOptions( mc::FFModel::EqnRole::INTERIOR, 0 );
  auto const ii = mc::FFModel::EqnOptions( mc::FFModel::EqnRole::INITIAL,  0 );

  IVP.add_equation( OpP( X0, t ) - TF * X1,                       {t}, {T_INT},          io );
  IVP.add_equation( OpP( X1, t ) - TF * ( U * X0 - 2. * X1 ) * t, {t}, {T_INT},          io );
  IVP.add_equation( X0 - 0.,                                      {t}, {mc::FFDom::LB},  ii );
  IVP.add_equation( X1 - X10,                                     {t}, {mc::FFDom::LB},  ii );

  IVP.add_output( OpEval( X0 - 1., t, tf ) );
  IVP.add_output( OpEval( X1,      t, tf ) );

  IVP.options.INTMETH   = mc::BASE_CVODES::Options::MSBDF;
  IVP.options.NLINSOL   = mc::BASE_CVODES::Options::NEWTON;
  IVP.options.LINSOL    = mc::BASE_CVODES::Options::SPARSE;
  IVP.options.NMAX      = 0;
  IVP.options.DISPLAY   = verbose? 1: 0;
  IVP.options.ATOL      = IVP.options.ATOLB = IVP.options.ATOLS = 1e-9;
  IVP.options.RTOL      = IVP.options.RTOLB = IVP.options.RTOLS = 1e-9;

  // Registration BEFORE setup(): FFModel::control_ndof reads the declaration, and setup() reindexes,
  // so this is the path the OCFESLV drivers (which register after setup) never exercise.
  if( reg & REG_TF  ) IVP.register_control( TF  );
  if( reg & REG_X10 ) IVP.register_control( X10 );
  if( reg & REG_U   ) IVP.register_control( U   );

  if( !IVP.setup() ){
    std::cerr << "  setup failed: " << IVP.extract_error() << std::endl;
    return R;
  }
  R.np  = IVP.np();
  R.nf  = IVP.nf();
  R.ncd = IVP.n_control_dof();

  if( IVP.solve_fsens( P0 ) != mc::ODESLVS_CVODES::STATUS::NORMAL ) return R;
  R.fsa = IVP.val_function_gradient();
  if( IVP.solve_asens( P0 ) != mc::ODESLVS_CVODES::STATUS::NORMAL ) return R;
  R.asa = IVP.val_function_gradient();
  R.ok  = true;
  return R;
}

static int npass = 0, nfail = 0;
void check( bool c, std::string const& what )
{ std::cout << "  " << (c? "PASS  ": "FAIL  ") << what << std::endl; c? ++npass: ++nfail; }
double rel( double a, double b ){ double s=std::fabs(b); return std::fabs(a-b)/( s>1e-12? s: 1. ); }

//! @brief Build test3's model on @p DAG into @p IVP; returns the handles the gate needs.
struct H { mc::FFVar t, X0, X1, TF, X10, U; };
H build( mc::FFGraph& DAG, mc::ODESLVS_CVODES& IVP, size_t n_node = 1 )
{
  H h; h.t = DAG.add_var("t"); h.X0 = DAG.add_var("x0(t)"); h.X1 = DAG.add_var("x1(t)");
  h.TF = DAG.add_var("TF"); h.X10 = DAG.add_var("x10"); h.U = DAG.add_var("u(t)");
  mc::FFPartial OpP; mc::FFEval OpEval;
  IVP.add_domain( h.t, mc::FFDom( t0, tf, NS, mc::FFDom::LGR, 4 ) );
  IVP.add_state( h.X0, {h.t} ); IVP.add_state( h.X1, {h.t} );
  IVP.add_input( h.TF ); IVP.add_input( h.X10 ); IVP.add_input( h.U, { h.t }, mc::FFDom::LGL, n_node );
  IVP.set_evolution_domain( h.t ); IVP.update_ref( h.X0, 0. ); IVP.update_ref( h.X1, 0.5 );
  int const T_INT = mc::FFDom::ALL - mc::FFDom::LB;
  auto const io = mc::FFModel::EqnOptions( mc::FFModel::EqnRole::INTERIOR, 0 );
  auto const ii = mc::FFModel::EqnOptions( mc::FFModel::EqnRole::INITIAL,  0 );
  IVP.add_equation( OpP( h.X0, h.t ) - h.TF * h.X1,                           {h.t}, {T_INT},         io );
  IVP.add_equation( OpP( h.X1, h.t ) - h.TF * ( h.U * h.X0 - 2. * h.X1 ) * h.t, {h.t}, {T_INT},         io );
  IVP.add_equation( h.X0 - 0.,                                                {h.t}, {mc::FFDom::LB}, ii );
  IVP.add_equation( h.X1 - h.X10,                                             {h.t}, {mc::FFDom::LB}, ii );
  IVP.add_output( OpEval( h.X0 - 1., h.t, tf ) ); IVP.add_output( OpEval( h.X1, h.t, tf ) );
  IVP.options.LINSOL = mc::BASE_CVODES::Options::DENSE; IVP.options.DISPLAY = 0;
  IVP.options.ATOL = IVP.options.ATOLS = 1e-10; IVP.options.RTOL = IVP.options.RTOLS = 1e-10;
  return h;
}

//! @brief A fresh IVP on @p DAG, set up; @p regU registers u(t) as a control.
static mc::ODESLVS_CVODES* make_ivp( mc::FFGraph& DAG, H& h, bool regU, size_t n_node = 1 )
{
  auto* IVP = new mc::ODESLVS_CVODES( &DAG );  h = build( DAG, *IVP, n_node );
  if( regU ) IVP->register_control( h.U );
  if( !IVP->setup() ){ delete IVP; return nullptr; }
  return IVP;
}

//! @brief test3's model with the coefficient 2 as K: a declared CONSTANT (kconst) or a declared INPUT.
static mc::ODESLVS_CVODES* make_K( mc::FFGraph& DAG, H& h, mc::FFVar& K, bool kconst )
{
  auto* IVP = new mc::ODESLVS_CVODES( &DAG );
  h.t = DAG.add_var("t"); h.X0 = DAG.add_var("x0(t)"); h.X1 = DAG.add_var("x1(t)");
  h.TF = DAG.add_var("TF"); h.X10 = DAG.add_var("x10"); h.U = DAG.add_var("u(t)"); K = DAG.add_var("K");
  mc::FFPartial OpP; mc::FFEval OpEval;
  IVP->add_domain( h.t, mc::FFDom( t0, tf, NS, mc::FFDom::LGR, 4 ) );
  IVP->add_state( h.X0, {h.t} ); IVP->add_state( h.X1, {h.t} );
  IVP->add_input( h.TF ); IVP->add_input( h.X10 ); IVP->add_input( h.U, { h.t }, mc::FFDom::LGR, 1 );
  if( kconst ) IVP->FFModel::set_constant( { K }, { 2. } ); else IVP->add_input( K );
  IVP->set_evolution_domain( h.t ); IVP->update_ref( h.X0, 0. ); IVP->update_ref( h.X1, 0.5 );
  int const T_INT = mc::FFDom::ALL - mc::FFDom::LB;
  auto const io = mc::FFModel::EqnOptions( mc::FFModel::EqnRole::INTERIOR, 0 );
  auto const ii = mc::FFModel::EqnOptions( mc::FFModel::EqnRole::INITIAL,  0 );
  IVP->add_equation( OpP( h.X0, h.t ) - h.TF * h.X1,                               {h.t}, {T_INT},         io );
  IVP->add_equation( OpP( h.X1, h.t ) - h.TF * ( h.U * h.X0 - K * h.X1 ) * h.t,   {h.t}, {T_INT},         io );
  IVP->add_equation( h.X0 - 0.,                                                    {h.t}, {mc::FFDom::LB}, ii );
  IVP->add_equation( h.X1 - h.X10,                                                 {h.t}, {mc::FFDom::LB}, ii );
  IVP->add_output( OpEval( h.X0 - 1., h.t, tf ) ); IVP->add_output( OpEval( h.X1, h.t, tf ) );
  IVP->options.LINSOL = mc::BASE_CVODES::Options::DENSE; IVP->options.DISPLAY = 0;
  IVP->options.ATOL = IVP->options.ATOLS = 1e-10; IVP->options.RTOL = IVP->options.RTOLS = 1e-10;
  if( !IVP->setup() ){ delete IVP; return nullptr; }
  return IVP;
}

int main()
{
  std::cout << std::scientific << std::setprecision(6)
            << "================================================================\n"
            << "  test14_ff: ODESLV values BY NAME -- set/get_parameter_values, named\n"
            << "  solve_state/sensitivity/adjoint, val_function_gradient(var[,ndx_el])\n"
            << "================================================================\n";
  mc::FFGraph DAG; H h; mc::FFVar K; mc::ODESLVS_CVODES* IVP = make_K( DAG, h, K, true );
  if( !IVP ){ std::cout << "  setup failed\n"; return 1; }
  size_t const np = IVP->np(), nf = IVP->nf();
  double const Kv = 2.;  std::vector<double> const Kc{ Kv };

  std::cout << "\n--- set / get_parameter_values\n";
  std::vector<double> P( np, 0. );
  check( IVP->set_parameter_values( h.TF, { P0[0] }, P ) && IVP->set_parameter_values( h.X10, { P0[1] }, P )
      && IVP->set_parameter_values( h.U, { P0[2], P0[3], P0[4] }, P ), "set_parameter_values for TF, x10, u(t)" );
  bool byname = true;                                         // independent check: the minted NAMES
  for( size_t k=0; k<np; ++k ){ std::string const n = IVP->var_parameter()[k].name(); double want =
      n=="TF"? P0[0]: n=="x10"? P0[1]: P0[ 2 + std::stoi( n.substr( n.find('[')+1 ) ) ]; byname = byname && P[k] == want; }
  check( byname, "every slot of P holds the value of the parameter named there" );
  check( IVP->get_parameter_values( h.U, P ) == std::vector<double>{ P0[2], P0[3], P0[4] }, "get_parameter_values( u ) returns its 3 DOFs in order" );
  check( !IVP->set_parameter_values( h.U, { 1., 2. }, P ), "a wrong value count is refused" );
  check( !IVP->set_parameter_values( K, { 1. }, P ), "a constant is not a parameter: refused" );
  check( !IVP->fix_input( h.X10, { 0.5 } ) && !IVP->fixed_input( h.X10 ), "fix_input after setup(): REFUSED on ODESLV (it would silently do nothing)" );

  std::cout << "\n--- named solve_state == flat solve_state\n";
  IVP->solve( P, Kc ); auto const Fflat = IVP->val_function();
  IVP->solve( { { h.TF, { P0[0] } }, { h.X10, { P0[1] } }, { h.U, { P0[2], P0[3], P0[4] } }, { K, { Kv } } } );
  auto const Fn = IVP->val_function();
  double w = 0.; for( size_t j=0;j<nf;++j ) w = std::max( w, rel( Fn[j], Fflat[j] ) );
  { std::ostringstream o; o << "named == flat (worst rel " << w << ")"; check( w == 0., o.str() ); }
  IVP->solve( { { K, { Kv } }, { h.U, [&]( mc::FFModel::DofIndex const& d ){ return P0[ 2 + d.element.at( h.t ) ]; } },
                      { h.X10, { P0[1] } }, { h.TF, { P0[0] } } } );
  auto const Fg = IVP->val_function(); w = 0.; for( size_t j=0;j<nf;++j ) w = std::max( w, rel( Fg[j], Fflat[j] ) );
  { std::ostringstream o; o << "listing order reversed, u(t) by generator: identical (worst rel " << w << ")"; check( w == 0., o.str() ); }

  std::cout << "\n--- refusals, each naming the entry\n";
  auto refused = [&]( std::function<void()> f, std::string const& frag ){
    try{ f(); } catch( std::invalid_argument& e ){ return std::string( e.what() ).find( frag ) != std::string::npos; } return false; };
  check( refused( [&]{ IVP->solve( { { h.TF, {P0[0]} }, { h.U, {P0[2],P0[3],P0[4]} }, { K, {Kv} } } ); }, "missing: x10" ), "input missing: refused, naming x10" );
  check( refused( [&]{ IVP->solve( { { h.TF, {P0[0]} }, { h.X10, {P0[1]} }, { h.U, {P0[2],P0[3],P0[4]} } } ); }, "missing: K" ), "constant missing: refused, naming K" );
  check( refused( [&]{ IVP->solve( { { h.TF, {P0[0]} }, { h.TF, {P0[0]} }, { h.X10, {P0[1]} }, { h.U, {P0[2],P0[3],P0[4]} }, { K, {Kv} } } ); }, "more than once" ), "input given twice: refused" );
  check( refused( [&]{ IVP->solve( { { h.TF, {P0[0]} }, { h.X10, {P0[1]} }, { h.U, {P0[2],P0[3]} }, { K, {Kv} } } ); }, "u(t) has 3 DOFs" ), "wrong count: refused, both counts" );
  check( refused( [&]{ IVP->solve( { { h.TF, {P0[0]} }, { h.X10, {P0[1]} }, { h.U, {P0[2],P0[3],P0[4]} }, { K, {Kv, 1.} } } ); }, "exactly one value" ), "constant with two values: refused" );

  std::cout << "\n--- val_function_gradient( var [, ndx_el] ) after named solve_sensitivity / solve_adjoint\n";
  IVP->clear_controls();                                           // reference: every parameter a direction
  IVP->solve_fsens( P, Kc ); auto const Gall = IVP->val_function_gradient();
  IVP->register_control( h.U ); IVP->register_control( h.TF );     // now the directions are {TF, u}
  for( int adj : { 0, 1 } ){
    std::vector<mc::FFModel::InputVal> const vIn{ { h.TF, {P0[0]} }, { h.X10, {P0[1]} }, { h.U, {P0[2],P0[3],P0[4]} }, { K, {Kv} } };
    if( adj ) IVP->solve_asens( vIn ); else IVP->solve_fsens( vIn );
    std::string const tag = adj? "adjoint": "forward";
    auto const GU = IVP->val_function_gradient( h.U ), GT = IVP->val_function_gradient( h.TF );
    auto const iU = IVP->parameter_index( h.U ), iT = IVP->parameter_index( h.TF );
    double wg = 0.;
    check( GU.size() == 3 && GT.size() == 1, tag + ": u has 3 rows, TF 1" );
    for( size_t k=0;k<3 && GU.size()==3;++k ) for( size_t j=0;j<nf;++j ) wg = std::max( wg, rel( GU[k][j], Gall[ iU[k] ][j] ) );
    for( size_t j=0;j<nf && GT.size()==1;++j ) wg = std::max( wg, rel( GT[0][j], Gall[ iT[0] ][j] ) );
    { std::ostringstream o; o << tag << ": rows for u and TF == the full-gradient rows of their parameters (worst rel " << wg << ")"; check( wg < ( adj? 1e-5: 1e-7 ), o.str() ); }   // two FSA runs over DIFFERENT direction sets: integrator tolerance
    double we = 0.;
    for( size_t e=0;e<NS;++e ){ auto const Ge = IVP->val_function_gradient( h.U, {{h.t,e}} );
      if( Ge.size() != 1 ){ we = 1.; break; } for( size_t j=0;j<nf;++j ) we = std::max( we, rel( Ge[0][j], GU[e][j] ) ); }
    check( we == 0., tag + ": val_function_gradient( u, {t:e} ) == row e of u's block, every element" );
  }
  check( refused( [&]{ IVP->val_function_gradient( h.X10 ); }, "not a registered control" ), "x10 is not a control: refused" );
  check( refused( [&]{ IVP->val_function_gradient( h.U, {{h.t,NS}} ); }, "element out of range" ), "element out of range: refused" );
  delete IVP;

  std::cout << "\n--- n_node = 2: per-element blocks of 2 rows\n";
  { mc::FFGraph DAG2; H h2; mc::ODESLVS_CVODES* I2 = make_ivp( DAG2, h2, false, 2 );
    std::vector<double> uv( 6 ); for( size_t k=0;k<6;++k ) uv[k] = 0.35 + 0.1*k;
    std::vector<mc::FFModel::InputVal> const vIn{ { h2.TF, {P0[0]} }, { h2.X10, {P0[1]} },
      { h2.U, [&]( mc::FFModel::DofIndex const& d ){ return uv[ 2*d.element.at( h2.t ) + d.node.at( h2.t ) ]; } } };
    I2->solve_fsens( vIn ); auto const G2 = I2->val_function_gradient();
    I2->register_control( h2.U ); I2->solve_fsens( vIn );
    auto const iU = I2->parameter_index( h2.U ); double we = 0.; bool shape = true;
    for( size_t e=0;e<NS;++e ){ auto const Ge = I2->val_function_gradient( h2.U, {{h2.t,e}} ); shape = shape && Ge.size() == 2;
      for( size_t j=0;j<2 && Ge.size()==2;++j ) for( size_t f=0;f<nf;++f ) we = std::max( we, rel( Ge[j][f], G2[ iU[2*e+j] ][f] ) ); }
    std::ostringstream o; o << "u(t) with n_node=2: each element 2 rows == the full-gradient rows of u[e,0], u[e,1] (worst rel " << we << ")";
    check( shape && we < 1e-7, o.str() );   // 6 vs 8 directions: sensitivities are in CVODES' error control
    delete I2; }

  std::cout << "\n  test14_ff: " << npass << " passed, " << nfail << " failed -- " << (nfail? "FAILURES":"ALL PASS") << "\n";
  return nfail? 1: 0;
}
