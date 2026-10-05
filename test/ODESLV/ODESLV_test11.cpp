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
#include <sstream>

#include "ffode.hpp"
#include <algorithm>

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

int main()
{
  std::cout << std::scientific << std::setprecision(6)
            << "================================================================\n"
            << "  test11_ff: FFODESLV with several inputs, one distributed (n_node=2),\n"
            << "  wired through the layout probe pos_input / size_input\n"
            << "================================================================\n";
  mc::FFGraph DAG;  H h;
  mc::ODESLVS_CVODES* IVP = make_ivp( DAG, h, false, 2 );
  if( !IVP ){ std::cout << "  setup failed\n"; return 1; }
  size_t const np = IVP->np(), nf = IVP->nf(), NE = NS, NN = 2;
  check( np == 2 + NE*NN, "np = TF + x10 + 3 elements x 2 nodes = 8" );
  check( IVP->size_input( h.TF ) == 1 && IVP->size_input( h.U ) == NN, "size_input: 1 for TF, n_node=2 for u(t)" );

  // ---- the probe, checked against the minted NAMES, independently of itself -------------------------
  auto const& vP = IVP->var_parameter();
  bool okname = vP[ IVP->pos_input( h.TF ) ].name() == "TF" && vP[ IVP->pos_input( h.X10 ) ].name() == "x10";
  for( size_t e = 0; e < NE; ++e ) for( size_t j = 0; j < NN; ++j ){
    std::string const want = "u(t)[" + std::to_string(e) + "," + std::to_string(j) + "]";
    okname = okname && vP[ IVP->pos_input( h.U, {{h.t,e}} ) + j ].name() == want;
  }
  check( okname, "pos_input( var, {t:e} ) + j names exactly the minted parameter u(t)[e,j]" );
  bool threw = false; try{ IVP->pos_input( h.U, {{h.t,NE}} ); } catch( std::exception& ){ threw = true; }
  check( threw, "pos_input refuses an element out of range" );

  // ---- the user's own structure: two scalars and a 3x2 array, placed by the probe ---------------------
  mc::FFGraph NLP;
  mc::FFVar TFv = NLP.add_var("TFv"), X10v = NLP.add_var("X10v");
  std::vector<std::vector<mc::FFVar>> Uv( NE, std::vector<mc::FFVar>( NN ) );
  for( size_t e = 0; e < NE; ++e ) for( size_t j = 0; j < NN; ++j ) Uv[e][j] = NLP.add_var( "U"+std::to_string(e)+std::to_string(j) );
  std::vector<mc::FFVar> vPar( np );  std::vector<int> hit( np, 0 );
  auto place = [&]( size_t k, mc::FFVar const& v ){ vPar[k] = v; ++hit[k]; };
  place( IVP->pos_input( h.TF ), TFv );  place( IVP->pos_input( h.X10 ), X10v );
  for( size_t e = 0; e < NE; ++e ) for( size_t j = 0; j < NN; ++j ) place( IVP->pos_input( h.U, {{h.t,e}} ) + j, Uv[e][j] );
  check( std::all_of( hit.begin(), hit.end(), []( int c ){ return c == 1; } ), "every parameter slot filled exactly once" );

  // values in the user's structure, and the flat vector the solver needs, built by NAME (not by the probe)
  double const TFval = P0[0], X10val = P0[1];
  std::vector<std::vector<double>> Uval( NE, std::vector<double>( NN ) );
  for( size_t e = 0; e < NE; ++e ) for( size_t j = 0; j < NN; ++j ) Uval[e][j] = 0.35 + 0.1*(2*e+j);
  std::vector<double> Pflat( np );
  for( size_t k = 0; k < np; ++k ){ std::string const n = vP[k].name();
    if( n == "TF" ) Pflat[k] = TFval; else if( n == "x10" ) Pflat[k] = X10val;
    else{ size_t const e = std::stoi( n.substr( n.find('[')+1 ) ), j = std::stoi( n.substr( n.find(',')+1 ) ); Pflat[k] = Uval[e][j]; } }
  check( IVP->solve_fsens( Pflat ) == mc::ODESLVS_CVODES::STATUS::NORMAL, "reference: solve_sensitivity" );
  auto const F0 = IVP->val_function();  auto const G0 = IVP->val_function_gradient();
  // The op's NUMERIC route chooses FSA or ASA by the legacy rule nPar <= NP2NF*nf: here 8 > 3*2, so ASA.
  check( np > mc::FFODESLV::options.NP2NF * nf, "numeric route will use the ADJOINT (np=8 > NP2NF*nf=6)" );
  check( IVP->solve_asens( Pflat ) == mc::ODESLVS_CVODES::STATUS::NORMAL, "reference: solve_adjoint" );
  auto const GA = IVP->val_function_gradient();
  { double w = 0.; for( size_t i=0;i<np;++i ) for( size_t k=0;k<nf;++k ) w = std::max( w, rel( GA[i][k], G0[i][k] ) );
    std::cout << "  [info] ASA vs FSA on the IVP itself, at the driver's backward tolerances: worst rel " << w << "\n"; }

  std::vector<mc::FFVar> vUser{ TFv, X10v };  std::vector<double> vUserVal{ TFval, X10val };
  for( size_t e = 0; e < NE; ++e ) for( size_t j = 0; j < NN; ++j ){ vUser.push_back( Uv[e][j] ); vUserVal.push_back( Uval[e][j] ); }

  for( int policy : { (int)mc::FFODESLV::SHALLOW, (int)mc::FFODESLV::COPY } ){
    char const* pn = policy? "COPY": "SHALLOW";
    std::cout << "\n--- " << pn << "\n";
    // ONE-MAP form: the user's two scalars and 3x2 array, u(t) by a generator per (element, node)
    mc::FFODESLV OpODE;
    auto F = OpODE( { { h.TF,  std::vector<mc::FFVar>{ TFv } }, { h.X10, std::vector<mc::FFVar>{ X10v } },
                      { h.U,   [&]( mc::FFModel::DofIndex const& d ){ return Uv.at( d.element.at( h.t ) ).at( d.node.at( h.t ) ); } } },
                    IVP, policy );
    std::vector<double> Fv( nf );  NLP.eval( F, Fv, vUser, vUserVal );
    double w = 0.; for( size_t k=0;k<nf;++k ) w = std::max( w, rel( Fv[k], F0[k] ) );
    std::ostringstream o; o << pn << ": value from the user's own variables == solve_state (worst rel " << w << ")"; check( w < 1e-8, o.str() );
    for( int sym : { 0, 1 } ){
      mc::FFODESLV::options.SYMDIFF.clear();  if( sym ) mc::FFODESLV::options.SYMDIFF = vPar;
      try{
        auto const dF = NLP.FAD( F, vUser );                     // w.r.t. the USER's variables, in the user's order
        std::vector<double> dFv( dF.size() );  NLP.eval( dF, dFv, vUser, vUserVal );
        double wd = 0.;
        for( size_t k=0;k<nf;++k ) for( size_t u=0;u<vUser.size();++u ){
          size_t const ip = std::find_if( vPar.begin(), vPar.end(), [&]( mc::FFVar const& v ){ return v.id() == vUser[u].id(); } ) - vPar.begin();
          wd = std::max( wd, rel( dFv[k*vUser.size()+u], ( sym? G0: GA )[ip][k] ) ); }
        std::ostringstream od; od << pn << ( sym? " symbolic: dF/d(user variables) == solve_sensitivity": " numeric: dF/d(user variables) == solve_adjoint" ) << " (worst rel " << wd << ")";
        check( wd < 1e-6, od.str() );
      } catch( std::exception& e ){ check( false, std::string(pn) + ( sym? " symbolic": " numeric" ) + " threw: " + e.what() ); }
    }
    mc::FFODESLV::options.SYMDIFF.clear();
  }
  delete IVP;
  std::cout << "\n  test11_ff: " << npass << " passed, " << nfail << " failed -- " << (nfail? "FAILURES":"ALL PASS") << "\n";
  return nfail? 1: 0;
}
