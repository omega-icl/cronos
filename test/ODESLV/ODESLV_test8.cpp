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

#include "odeslvs_cvodes.hpp"

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

int main()
{
  std::cout << std::scientific << std::setprecision(6)
            << "================================================================\n"
            << "  test8_ff: FFModel::fdiff (Gateaux) vs CVODES forward sensitivity\n"
            << "================================================================\n";

  // ---- reference: FSA on the original -------------------------------------------------------------
  mc::FFGraph DAG; mc::ODESLVS_CVODES IVP( &DAG ); H h = build( DAG, IVP );
  if( !IVP.setup() ){ std::cout << "  setup failed\n"; return 1; }
  check( IVP.solve_fsens( P0 ) == mc::ODESLVS_CVODES::STATUS::NORMAL, "reference FSA" );
  auto const G = IVP.val_function_gradient();     // [direction][function]; directions = _mP order TF,x10,u0,u1,u2
  size_t const nf = IVP.nf();

  auto run_case = [&]( mc::FFVar const& u, char const* name, std::vector<size_t> const& refcols ){
    std::cout << "\n--- fdiff along " << name << "  (" << refcols.size() << " DOF)\n";
    mc::ODESLVS_CVODES SENS( &DAG ); SENS.options = IVP.options;
    mc::FFVar du; std::vector<mc::FFVar> vS; std::string err;
    bool const okf = IVP.fdiff( SENS, u, du, vS, err );
    check( okf, std::string("fdiff declares the sensitivity model") + ( err.empty()? "": "  -- " + err ) );
    if( !okf ) return;
    check( vS.size() == 2, "one sensitivity state per state" );
    if( !SENS.setup() ){ check( false, "setup of the sensitivity model" ); return; }
    check( SENS.nx() == 2*IVP.nx(), "sensitivity model has 2*nx states" );
    check( SENS.nf() == 2*nf,       "sensitivity model has 2*nf outputs (F then dF)" );
    // parameter vector of SENS: place the originals by name, du's DOFs as e_k
    auto const& vP = SENS.var_parameter();
    std::vector<size_t> ndx_du;
    std::vector<double> P( vP.size(), 0. );
    std::string const dun = du.name();
    for( size_t i = 0; i < vP.size(); ++i ){
      std::string const n = vP[i].name();
      if      ( n == "TF" )  P[i] = P0[0];
      else if ( n == "x10" ) P[i] = P0[1];
      else if ( n.rfind("u(t)",0) == 0 ) P[i] = P0[2 + std::stoi( n.substr( n.find('[')+1 ) )];
      else if ( n.rfind(dun,0) == 0 )    ndx_du.push_back( i );
    }
    std::ostringstream o; o << "du has " << ndx_du.size() << " DOF(s) in the sensitivity model (expected " << refcols.size() << ")";
    check( ndx_du.size() == refcols.size(), o.str() );
    double worst = 0.;
    for( size_t k = 0; k < ndx_du.size() && k < refcols.size(); ++k ){
      std::vector<double> Pk = P; Pk[ ndx_du[k] ] = 1.;           // unit perturbation of DOF k
      if( SENS.solve( Pk ) != mc::ODESLVS_CVODES::STATUS::NORMAL ){ check( false, "solve_state on the sensitivity model" ); return; }
      auto const& F = SENS.val_function();                          // [F0, F1, dF0, dF1]
      for( size_t j = 0; j < nf; ++j ) worst = std::max( worst, rel( F[nf+j], G[refcols[k]][j] ) );
    }
    std::ostringstream w; w << "dF along every DOF of " << name << " == FSA columns (worst rel " << worst << ")";
    check( worst < 1e-6, w.str() );
  };
  { std::cout << "\n--- probe: fdiff on a model that has NOT been set up\n";
    mc::FFGraph DAG2; mc::ODESLVS_CVODES IVP2( &DAG2 ); H h2 = build( DAG2, IVP2 );
    mc::ODESLVS_CVODES S2( &DAG2 ); mc::FFVar du2; std::vector<mc::FFVar> vS2; std::string e2;
    bool const ok2 = IVP2.fdiff( S2, h2.TF, du2, vS2, e2 );
    check( ok2, std::string("fdiff before setup()") + ( e2.empty()? "": "  -- " + e2 ) ); }
  run_case( h.TF, "TF   (time-invariant)",    { 0 } );
  run_case( h.U,  "u(t) (piecewise constant)", { 2, 3, 4 } );

  // ---- n_node = 2: the case the Gateaux form exists for -- 6 nodal DOFs, no DOF exists in the model DAG
  std::cout << "\n=== n_node = 2 (piecewise LINEAR control, 6 DOFs) ===\n";
  mc::FFGraph DAG3; mc::ODESLVS_CVODES IVP3( &DAG3 ); H h3 = build( DAG3, IVP3, 2 );
  if( !IVP3.setup() ){ check( false, "setup n_node=2" ); }
  else{
    std::vector<double> P3( IVP3.np(), 0.5 ); P3[0] = P0[0]; P3[1] = P0[1];
    for( size_t k = 0; k < 6; ++k ) P3[2+k] = 0.35 + 0.1*k;                 // a genuinely varying profile
    check( IVP3.solve_fsens( P3 ) == mc::ODESLVS_CVODES::STATUS::NORMAL, "reference FSA, n_node=2" );
    auto const G3 = IVP3.val_function_gradient();
    check( G3.size() == 8, "reference has 8 directions (TF, x10, 6 nodal values)" );
    mc::ODESLVS_CVODES S3( &DAG3 ); S3.options = IVP3.options;
    mc::FFVar du3; std::vector<mc::FFVar> vS3; std::string e3;
    check( IVP3.fdiff( S3, h3.U, du3, vS3, e3 ), "fdiff along u(t), n_node=2" );
    if( S3.setup() ){
      auto const& vP = S3.var_parameter(); std::vector<size_t> ndx; std::vector<double> P( vP.size(), 0. );
      for( size_t i = 0; i < vP.size(); ++i ){
        std::string const n = vP[i].name();
        if( n == "TF" ) P[i] = P3[0]; else if( n == "x10" ) P[i] = P3[1];
        else if( n.rfind("u(t)",0) == 0 ){
          size_t const e = std::stoi( n.substr( n.find('[')+1 ) ), j = std::stoi( n.substr( n.find(',')+1 ) );
          P[i] = P3[2 + 2*e + j]; }
        else if( n.rfind( du3.name(),0 ) == 0 ) ndx.push_back( i );
      }
      check( ndx.size() == 6, "du has 6 DOFs" );
      double worst = 0.;
      for( size_t k = 0; k < ndx.size(); ++k ){
        std::vector<double> Pk = P; Pk[ndx[k]] = 1.;
        if( S3.solve( Pk ) != mc::ODESLVS_CVODES::STATUS::NORMAL ){ check( false, "solve" ); break; }
        auto const& F = S3.val_function();
        for( size_t j = 0; j < nf; ++j ) worst = std::max( worst, rel( F[nf+j], G3[2+k][j] ) );
      }
      std::ostringstream w; w << "dF along all 6 nodal DOFs == FSA columns (worst rel " << worst << ")";
      check( worst < 1e-6, w.str() );
    } else check( false, "setup of the n_node=2 sensitivity model" );
  }

  std::cout << "\n  test8_ff: " << npass << " passed, " << nfail << " failed -- " << (nfail? "FAILURES":"ALL PASS") << "\n";
  return nfail? 1: 0;
}
