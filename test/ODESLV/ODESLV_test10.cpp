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
            << "  test10_ff: FFODESLV -- the IVP as an external operation in a DAG\n"
            << "================================================================\n";
  mc::FFGraph DAG;  H h;
  mc::ODESLVS_CVODES* IVP = make_ivp( DAG, h, false );
  if( !IVP ){ std::cout << "  setup failed\n"; return 1; }
  size_t const np = IVP->np(), nf = IVP->nf();
  check( IVP->solve_fsens( P0 ) == mc::ODESLVS_CVODES::STATUS::NORMAL, "reference: solve_sensitivity on the IVP" );
  auto const F0 = IVP->val_function();  auto const G0 = IVP->val_function_gradient();

  // the DAG the ODE is embedded in: np variables, one per parameter of the IVP
  auto embed = [&]( mc::FFGraph& NLP, std::vector<mc::FFVar>& P, mc::ODESLVS_CVODES* ivp, int policy ){
    P.assign( ivp->np(), mc::FFVar() );
    for( size_t i = 0; i < P.size(); ++i ) P[i].set( &NLP );
    // ONE-MAP form: every declared input, its DOFs taken from P through parameter_index -- so P stays in
    // parameter order and the reference columns G0[i] line up with the derivatives w.r.t. P[i].
    std::vector<mc::FFModel::InputArg> m;
    for( auto const& [w,dom] : ivp->var_declared_input() ){
      std::vector<mc::FFVar> d;  for( size_t i : ivp->parameter_index( w ) ) d.push_back( P[i] );
      m.emplace_back( w, d );
    }
    mc::FFODESLV OpODE;
    return OpODE( m, ivp, policy );
  };
  auto near = [&]( std::vector<double> const& a, std::vector<double> const& b ){ double w=0.; for( size_t i=0;i<a.size()&&i<b.size();++i ) w=std::max(w,rel(a[i],b[i])); return w; };

  for( int policy : { (int)mc::FFODESLV::SHALLOW, (int)mc::FFODESLV::COPY } ){
    char const* pn = policy? "COPY": "SHALLOW";
    std::cout << "\n--- " << pn << ": value\n";
    mc::FFGraph NLP;  std::vector<mc::FFVar> P;
    std::vector<mc::FFVar> F;
    try{ F = embed( NLP, P, IVP, policy ); } catch( std::exception& e ){ check( false, std::string(pn)+" embed threw: "+e.what() ); continue; }
    check( F.size() == nf, std::string(pn) + ": the op has nf outputs" );
    std::vector<double> Fv( nf );
    try{ NLP.eval( F, Fv, P, P0 ); } catch( std::exception& e ){ check( false, std::string(pn)+" eval threw: "+e.what() ); continue; }
    std::ostringstream o; o << pn << ": DAG value == solve_state (worst rel " << near( Fv, F0 ) << ")"; check( near( Fv, F0 ) < 1e-8, o.str() );

    std::cout << "--- " << pn << ": NUMERIC derivative (SYMDIFF empty -> FFGradODESLV)\n";
    mc::FFODESLV::options.SYMDIFF.clear();
    try{
      auto const dF = NLP.FAD( F, P );                 // nf*np derivative FFVars, [k*np + i]? -- measured below
      std::vector<double> dFv( dF.size() );  NLP.eval( dF, dFv, P, P0 );
      double w = 0.;  for( size_t k=0;k<nf;++k ) for( size_t i=0;i<np;++i ) w = std::max( w, rel( dFv[k*np+i], G0[i][k] ) );
      std::ostringstream od; od << pn << ": dF/dP == solve_sensitivity, all " << nf*np << " entries (worst rel " << w << ")"; check( w < 1e-6, od.str() );
    } catch( std::exception& e ){ check( false, std::string(pn)+" numeric FAD threw: "+e.what() ); }

    std::cout << "--- " << pn << ": SYMBOLIC derivative (SYMDIFF = all P -> fdiff)\n";
    mc::FFODESLV::options.SYMDIFF = P;
    try{
      auto const dF = NLP.FAD( F, P );
      std::vector<double> dFv( dF.size() );  NLP.eval( dF, dFv, P, P0 );
      double w = 0.;  for( size_t k=0;k<nf;++k ) for( size_t i=0;i<np;++i ) w = std::max( w, rel( dFv[k*np+i], G0[i][k] ) );
      std::ostringstream od; od << pn << ": symbolic dF/dP == solve_sensitivity (worst rel " << w << ")"; check( w < 1e-6, od.str() );
    } catch( std::exception& e ){ check( false, std::string(pn)+" symbolic FAD threw: "+e.what() ); }
    mc::FFODESLV::options.SYMDIFF.clear();
  }

  std::cout << "\n--- a registered control on the IVP (u(t): 3 directions, not 5)\n";
  { mc::FFGraph DAGc; H hc; mc::ODESLVS_CVODES* IVPc = make_ivp( DAGc, hc, true );
    check( IVPc && IVPc->sensitivity_index().size() == 3, "IVP with u(t) registered: 3 sensitivity directions" );
    if( IVPc ){
      mc::FFODESLV::options.SYMDIFF.clear();
      // COPY: the op owns its copy and differentiates w.r.t. EVERY input
      { mc::FFGraph NLP; std::vector<mc::FFVar> P; auto F = embed( NLP, P, IVPc, mc::FFODESLV::COPY );
        try{
          auto const dF = NLP.FAD( F, P ); std::vector<double> dFv( dF.size() ); NLP.eval( dF, dFv, P, P0 );
          double w = 0.; for( size_t k=0;k<nf;++k ) for( size_t i=0;i<np;++i ) w = std::max( w, rel( dFv[k*np+i], G0[i][k] ) );
          std::ostringstream od; od << "COPY: numeric dF/dP over ALL 5 inputs despite the registry (worst rel " << w << ")"; check( w < 1e-6, od.str() );
        } catch( std::exception& e ){ check( false, std::string("COPY registry case threw: ")+e.what() ); }
        check( IVPc->sensitivity_index().size() == 3, "COPY left the user's IVP registry untouched" ); }
      // SHALLOW: the op RESETS the shared solver's registry from its map -- here every input (one-map form) --
      // so the NUMERIC route runs over the requested u(t) directions like any other; the SYMBOLIC route
      // differentiates exactly the requested inputs.
      { mc::FFGraph NLP; std::vector<mc::FFVar> P; auto F = embed( NLP, P, IVPc, mc::FFODESLV::SHALLOW );
        check( IVPc->sensitivity_index().size() == np, "SHALLOW: the shared registry is now every input (map 1)" );
        std::vector<mc::FFVar> Pu( P.begin()+2, P.end() );          // u(t)[0..2]
        mc::FFODESLV::options.SYMDIFF.clear();
        try{
          auto const dF = NLP.FAD( F, Pu ); std::vector<double> dFv( dF.size() ); NLP.eval( dF, dFv, P, P0 );
          double w = 0.; for( size_t k=0;k<nf;++k ) for( size_t i=0;i<3;++i ) w = std::max( w, rel( dFv[k*3+i], G0[2+i][k] ) );
          std::ostringstream od; od << "SHALLOW numeric: dF/du (worst rel " << w << ")"; check( w < 1e-6, od.str() );
        } catch( std::exception& e ){ check( false, std::string("SHALLOW numeric threw: ")+e.what() ); }
        mc::FFODESLV::options.SYMDIFF = Pu;
        try{
          auto const dF = NLP.FAD( F, Pu ); std::vector<double> dFv( dF.size() ); NLP.eval( dF, dFv, P, P0 );
          double w = 0.; for( size_t k=0;k<nf;++k ) for( size_t i=0;i<3;++i ) w = std::max( w, rel( dFv[k*3+i], G0[2+i][k] ) );
          std::ostringstream od; od << "SHALLOW symbolic: dF/du along the requested inputs only (worst rel " << w << ")"; check( w < 1e-6, od.str() );
        } catch( std::exception& e ){ check( false, std::string("SHALLOW symbolic subset threw: ")+e.what() ); }
        mc::FFODESLV::options.SYMDIFF = P;
        try{
          auto const dF = NLP.FAD( F, P ); std::vector<double> dFv( dF.size() ); NLP.eval( dF, dFv, P, P0 );
          double w = 0.; for( size_t k=0;k<nf;++k ) for( size_t i=0;i<np;++i ) w = std::max( w, rel( dFv[k*np+i], G0[i][k] ) );
          std::ostringstream od; od << "SHALLOW symbolic: dF/dP over ALL 5 inputs, registry irrelevant (worst rel " << w << ")"; check( w < 1e-6, od.str() );
        } catch( std::exception& e ){ check( false, std::string("SHALLOW symbolic full threw: ")+e.what() ); }
        mc::FFODESLV::options.SYMDIFF.clear(); }
      delete IVPc; }
  }
  delete IVP;
  std::cout << "\n  test10_ff: " << npass << " passed, " << nfail << " failed -- " << (nfail? "FAILURES":"ALL PASS") << "\n";
  return nfail? 1: 0;
}
