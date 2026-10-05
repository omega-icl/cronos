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
            << "  test12_ff: FFODESLV TWO-MAP form -- map 1 sets the registry, map 2 carries\n"
            << "  every other input and constant through the DAG; SYMDIFF over constants\n"
            << "================================================================\n";
  mc::FFGraph DAG0; H h0; mc::ODESLVS_CVODES* REF = make_ivp( DAG0, h0, false );
  REF->solve_fsens( P0 ); auto const G0 = REF->val_function_gradient(); auto const F0 = REF->val_function();
  size_t const nf = REF->nf();
  std::vector<double> P1 = P0; P1[1] *= 1.3;  REF->solve( P1 ); auto const F1 = REF->val_function();  delete REF;

  auto dcol = [&]( mc::FFGraph& NLP, std::vector<mc::FFVar> const& F, std::vector<mc::FFVar> const& vX,
                   std::vector<double> const& vXv, std::vector<mc::FFVar> const& sym ){
    mc::FFODESLV::options.SYMDIFF = sym;
    auto const dF = NLP.FAD( F, vX ); std::vector<double> dFv( dF.size() ); NLP.eval( dF, dFv, vX, vXv );
    mc::FFODESLV::options.SYMDIFF.clear();  return dFv; };                                   // [k*nX + u]

  for( int policy : { (int)mc::FFODESLV::SHALLOW, (int)mc::FFODESLV::COPY } ){
    std::string const pn = policy? "COPY": "SHALLOW";
    std::cout << "\n--- " << pn << ": map 1 {TF flat, u(t) generator}, map 2 {x10}\n";
    mc::FFGraph DAG; H h; mc::ODESLVS_CVODES* IVP = make_ivp( DAG, h, true );             // caller registered {u}
    mc::FFGraph NLP;  mc::FFVar tf = NLP.add_var("tf"), x10 = NLP.add_var("x10");
    std::vector<mc::FFVar> uv( NS ); for( size_t e=0;e<NS;++e ) uv[e] = NLP.add_var( "u"+std::to_string(e) );
    mc::FFODESLV OpODE;  std::vector<mc::FFVar> F;
    try{
      F = OpODE( { { h.TF, std::vector<mc::FFVar>{ tf } },
                   { h.U,  [&]( mc::FFModel::DofIndex const& d ){ return uv.at( d.element.at( h.t ) ); } } },
                 { { h.X10, std::vector<mc::FFVar>{ x10 } } }, IVP, policy );
    } catch( std::exception& e ){ check( false, pn + " two-map embedding threw: " + e.what() ); delete IVP; continue; }
    size_t const nsen = IVP->sensitivity_index().size();
    check( policy? nsen == 3: nsen == 4, pn + ( policy? ": the caller's registry {u} was RESTORED"
                                                     : ": the shared registry is now map 1 {TF,u}" ) );
    std::vector<mc::FFVar> vX{ tf };  for( auto const& v : uv ) vX.push_back( v );  vX.push_back( x10 );
    std::vector<double> vXv{ P0[0], P0[2], P0[3], P0[4], P0[1] };
    size_t const nX = vX.size();  std::vector<size_t> const pcol{ 0, 2, 3, 4, 1 };
    std::vector<double> Fv( nf ); NLP.eval( F, Fv, vX, vXv );
    double w = 0.; for( size_t k=0;k<nf;++k ) w = std::max( w, rel( Fv[k], F0[k] ) );
    { std::ostringstream o; o << pn << ": value == solve_state (worst rel " << w << ")"; check( w < 1e-8, o.str() ); }
    // THE point of map 2: x10 is a DAG variable -- change its value, no rebuild
    std::vector<double> vXv1 = vXv; vXv1[4] = P1[1];  NLP.eval( F, Fv, vX, vXv1 );
    w = 0.; for( size_t k=0;k<nf;++k ) w = std::max( w, rel( Fv[k], F1[k] ) );
    { std::ostringstream o; o << pn << ": x10 driven through the DAG (x1.3, same op) == solve_state (worst rel " << w << ")"; check( w < 1e-8, o.str() ); }
    try{
      auto const dN = dcol( NLP, F, vX, vXv, {} );                          // numeric
      double wn = 0., z = 0.;
      for( size_t k=0;k<nf;++k ){ for( size_t u=0;u<4;++u ) wn = std::max( wn, rel( dN[k*nX+u], G0[pcol[u]][k] ) ); z = std::max( z, std::fabs( dN[k*nX+4] ) ); }
      { std::ostringstream o; o << pn << " numeric: dF/d(map 1) (worst rel " << wn << ")"; check( wn < 1e-6, o.str() ); }
      { std::ostringstream o; o << pn << " numeric: dF/dx10 (map 2) is ZERO by contract (max |.| " << z << ")"; check( z == 0., o.str() ); }
      auto const dS = dcol( NLP, F, vX, vXv, std::vector<mc::FFVar>( vX.begin(), vX.begin()+4 ) );
      double ws = 0.; for( size_t k=0;k<nf;++k ) for( size_t u=0;u<4;++u ) ws = std::max( ws, rel( dS[k*nX+u], G0[pcol[u]][k] ) );
      { std::ostringstream o; o << pn << " symbolic SYMDIFF=map 1: (worst rel " << ws << ")"; check( ws < 1e-6, o.str() ); }
      auto const dX = dcol( NLP, F, vX, vXv, { x10 } );                    // SYMDIFF reaches into map 2
      double wx = 0.; for( size_t k=0;k<nf;++k ) wx = std::max( wx, rel( dX[k*nX+4], G0[1][k] ) );
      { std::ostringstream o; o << pn << " symbolic SYMDIFF={x10}: a map-2 input differentiated symbolically (worst rel " << wx << ")"; check( wx < 1e-6, o.str() ); }
    } catch( std::exception& e ){ check( false, pn + " derivatives threw: " + e.what() ); }
    if( !policy ){
      IVP->clear_controls(); bool threw = false; std::string msg;
      try{ NLP.eval( F, Fv, vX, vXv ); } catch( std::exception& e ){ threw = true; msg = e.what(); }
      check( threw && msg.find("changed after this operation was built") != std::string::npos, "SHALLOW: evaluation refused after the caller changed the registry" );
    }
    delete IVP;
  }

  std::cout << "\n--- refusals at embedding\n";
  { mc::FFGraph DAG; H h; mc::ODESLVS_CVODES* IVP = make_ivp( DAG, h, false );
    mc::FFGraph NLP; mc::FFVar a = NLP.add_var("a"), b = NLP.add_var("b"); std::vector<mc::FFVar> u3{ a, b, a };
    auto refused = [&]( std::function<void()> f, std::string const& frag ){
      try{ f(); } catch( std::exception& e ){ return std::string( e.what() ).find( frag ) != std::string::npos; } return false; };
    mc::FFODESLV Op;
    check( refused( [&]{ Op( { { h.TF, {a} }, { h.U, u3 } }, std::vector<mc::FFModel::InputArg>{}, IVP, mc::FFODESLV::COPY ); }, "missing: x10" ),
           "an input in neither map: refused, naming it" );
    check( refused( [&]{ Op( { { h.TF, {a} }, { h.U, u3 } }, { { h.X10, {b} }, { h.TF, {b} } }, IVP, mc::FFODESLV::COPY ); }, "more than once" ),
           "an input in both maps: refused" );
    check( refused( [&]{ Op( { { h.TF, {a} }, { h.U, std::vector<mc::FFVar>{ a, b } } }, { { h.X10, {b} } }, IVP, mc::FFODESLV::COPY ); }, "u(t) has 3 DOFs" ),
           "a wrong DOF count: refused, naming the input and both counts" );
    check( IVP->sensitivity_index().size() == IVP->np(), "refusals left the solver's registry untouched" );
    delete IVP; }

  std::cout << "\n--- CONSTANTS: K declared as a constant; the reference is the same model with K an input\n";
  { mc::FFGraph DAGr; H hr; mc::FFVar Kr; mc::ODESLVS_CVODES* RK = make_K( DAGr, hr, Kr, false );
    std::vector<double> PK( RK->np() );
    PK[ RK->pos_input( hr.TF ) ] = P0[0]; PK[ RK->pos_input( hr.X10 ) ] = P0[1]; PK[ RK->pos_input( Kr ) ] = 2.;
    for( size_t e=0;e<NS;++e ) PK[ RK->pos_input( hr.U, {{hr.t,e}} ) ] = P0[2+e];
    RK->solve_fsens( PK ); auto const GK = RK->val_function_gradient(); auto const FK = RK->val_function();
    size_t const iK = RK->pos_input( Kr );  delete RK;

    mc::FFGraph DAG; H h; mc::FFVar K; mc::ODESLVS_CVODES* IVP = make_K( DAG, h, K, true );
    check( IVP && IVP->var_constant().size() == 1, "model with one declared constant K" );
    mc::FFGraph NLP;  mc::FFVar tf = NLP.add_var("tf"), x10 = NLP.add_var("x10"), k = NLP.add_var("k");
    std::vector<mc::FFVar> uv( NS ); for( size_t e=0;e<NS;++e ) uv[e] = NLP.add_var( "u" );
    mc::FFODESLV Op;
    check( [&]{ try{ Op( { { h.TF, {tf} }, { K, {k} } }, { { h.X10, {x10} }, { h.U, uv } }, IVP, mc::FFODESLV::COPY ); }
                catch( std::exception& e ){ return std::string( e.what() ).find("cannot be in map 1") != std::string::npos; } return false; }(),
           "a constant in map 1: refused (the solver's sensitivities run over inputs only)" );
    check( [&]{ try{ Op( { { h.TF, {tf} }, { h.U, uv }, { h.X10, {x10} } }, IVP, mc::FFODESLV::COPY ); }
                catch( std::exception& e ){ return std::string( e.what() ).find("missing: K") != std::string::npos; } return false; }(),
           "one-map form on a model with a constant: refused, naming the constant (use the two-map form)" );
    auto F = Op( { { h.TF, {tf} }, { h.U, uv } }, { { h.X10, {x10} }, { K, {k} } }, IVP, mc::FFODESLV::COPY );
    { mc::FFODESLV Op2; mc::FFVar& F0k = Op2( 0, { { h.TF, {tf} }, { h.U, uv } }, { { h.X10, {x10} }, { K, {k} } }, IVP, mc::FFODESLV::COPY );
      std::vector<mc::FFVar> one{ F0k }; std::vector<double> v1( 1 );
      std::vector<mc::FFVar> vXk{ tf }; for( auto const& v : uv ) vXk.push_back( v ); vXk.push_back( x10 ); vXk.push_back( k );
      NLP.eval( one, v1, vXk, std::vector<double>{ P0[0], P0[2], P0[3], P0[4], P0[1], 2. } );
      std::ostringstream o; o << "two-map, idep=0: the single output == F[0] (rel " << rel( v1[0], FK[0] ) << ")"; check( rel( v1[0], FK[0] ) < 1e-8, o.str() ); }
    std::vector<mc::FFVar> vX{ tf }; for( auto const& v : uv ) vX.push_back( v ); vX.push_back( x10 ); vX.push_back( k );
    std::vector<double> vXv{ P0[0], P0[2], P0[3], P0[4], P0[1], 2. };  size_t const nX = vX.size();
    std::vector<double> Fv( nf ); NLP.eval( F, Fv, vX, vXv );
    double w = 0.; for( size_t j=0;j<nf;++j ) w = std::max( w, rel( Fv[j], FK[j] ) );
    { std::ostringstream o; o << "value, K passed through the DAG as a constant (worst rel " << w << ")"; check( w < 1e-8, o.str() ); }
    try{
      auto const dN = dcol( NLP, F, vX, vXv, {} );
      double z = 0.; for( size_t j=0;j<nf;++j ) z = std::max( z, std::fabs( dN[j*nX+5] ) );
      { std::ostringstream o; o << "numeric: dF/dK is ZERO by contract (max |.| " << z << ")"; check( z == 0., o.str() ); }
      auto const dK = dcol( NLP, F, vX, vXv, { k } );
      double wk = 0.; for( size_t j=0;j<nf;++j ) wk = std::max( wk, rel( dK[j*nX+5], GK[iK][j] ) );
      { std::ostringstream o; o << "SYMDIFF={K}: symbolic derivative along a CONSTANT == FSA with K as an input (worst rel " << wk << ")"; check( wk < 1e-6, o.str() ); }
      auto const dB = dcol( NLP, F, vX, vXv, { tf, k } );
      double wb = 0.; for( size_t j=0;j<nf;++j ){ wb = std::max( wb, rel( dB[j*nX+0], GK[ 0 ][j] ) ); wb = std::max( wb, rel( dB[j*nX+5], GK[iK][j] ) ); }
      { std::ostringstream o; o << "SYMDIFF={TF,K}: an input and a constant together (worst rel " << wb << ")"; check( wb < 1e-6, o.str() ); }
    } catch( std::exception& e ){ check( false, std::string("constant derivatives threw: ") + e.what() ); }
    delete IVP; }

  std::cout << "\n--- n_node = 2: every input in map 1 (generator), map 2 empty\n";
  { mc::FFGraph DAG; H h; mc::ODESLVS_CVODES* IVP = make_ivp( DAG, h, false, 2 );
    std::vector<std::vector<double>> Uval( NS, std::vector<double>( 2 ) );
    for( size_t e=0;e<NS;++e ) for( size_t j=0;j<2;++j ) Uval[e][j] = 0.35 + 0.1*(2*e+j);
    std::vector<double> P2( IVP->np() );
    P2[ IVP->pos_input( h.TF ) ] = P0[0]; P2[ IVP->pos_input( h.X10 ) ] = P0[1];
    for( size_t e=0;e<NS;++e ) for( size_t j=0;j<2;++j ) P2[ IVP->pos_input( h.U, {{h.t,e}} ) + j ] = Uval[e][j];
    IVP->solve_fsens( P2 ); auto const G2 = IVP->val_function_gradient(); auto const F2 = IVP->val_function();
    mc::FFGraph NLP; mc::FFVar tf = NLP.add_var("tf"), x10 = NLP.add_var("x10");
    std::vector<std::vector<mc::FFVar>> U( NS, std::vector<mc::FFVar>( 2 ) );
    for( size_t e=0;e<NS;++e ) for( size_t j=0;j<2;++j ) U[e][j] = NLP.add_var( "U" );
    mc::FFODESLV Op;
    auto F = Op( { { h.U, [&]( mc::FFModel::DofIndex const& d ){ return U.at( d.element.at( h.t ) ).at( d.node.at( h.t ) ); } },
                   { h.X10, {x10} }, { h.TF, {tf} } }, IVP, mc::FFODESLV::COPY );                 // ONE-MAP form
    std::vector<mc::FFVar> vX{ tf, x10 }; std::vector<double> vXv{ P0[0], P0[1] };
    std::vector<size_t> pc{ IVP->pos_input( h.TF ), IVP->pos_input( h.X10 ) };
    for( size_t e=0;e<NS;++e ) for( size_t j=0;j<2;++j ){ vX.push_back( U[e][j] ); vXv.push_back( Uval[e][j] ); pc.push_back( IVP->pos_input( h.U, {{h.t,e}} ) + j ); }
    std::vector<double> Fv( nf ); NLP.eval( F, Fv, vX, vXv );
    double w = 0.; for( size_t j=0;j<nf;++j ) w = std::max( w, rel( Fv[j], F2[j] ) );
    { std::ostringstream o; o << "n_node=2: value (worst rel " << w << ")"; check( w < 1e-8, o.str() ); }
    { mc::FFODESLV Op1;                                                                     // ONE-MAP, single output
      mc::FFVar& F1 = Op1( 1, { { h.U, [&]( mc::FFModel::DofIndex const& d ){ return U.at( d.element.at( h.t ) ).at( d.node.at( h.t ) ); } },
                                { h.X10, {x10} }, { h.TF, {tf} } }, IVP, mc::FFODESLV::COPY );
      std::vector<mc::FFVar> one{ F1 }; std::vector<double> v1( 1 ); NLP.eval( one, v1, vX, vXv );
      std::ostringstream o; o << "one-map, idep=1: the single output == F[1] (rel " << rel( v1[0], F2[1] ) << ")"; check( rel( v1[0], F2[1] ) < 1e-8, o.str() ); }
    auto const dS = dcol( NLP, F, vX, vXv, vX ); size_t const nX = vX.size();
    double ws = 0.; for( size_t j=0;j<nf;++j ) for( size_t u=0;u<nX;++u ) ws = std::max( ws, rel( dS[j*nX+u], G2[pc[u]][j] ) );
    { std::ostringstream o; o << "n_node=2 symbolic: dF/d(8 inputs) (worst rel " << ws << ")"; check( ws < 1e-6, o.str() ); }
    delete IVP; }

  std::cout << "\n  test12_ff: " << npass << " passed, " << nfail << " failed -- " << (nfail? "FAILURES":"ALL PASS") << "\n";
  return nfail? 1: 0;
}
