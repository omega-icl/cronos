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

void check( bool cond, std::string const& what )
{
  std::cout << "  " << ( cond? "PASS  ": "FAIL  " ) << what << std::endl;
  cond? ++npass: ++nfail;
}

//! @brief Compare one selected column against the reference column it should reproduce.
void check_col( std::vector<std::vector<double>> const& sel, size_t j,
                std::vector<std::vector<double>> const& ref, size_t jref,
                std::string const& tag )
{
  if( j >= sel.size() || jref >= ref.size() ){
    check( false, tag + "  (column out of range: sel " + std::to_string(sel.size())
                      + ", ref " + std::to_string(ref.size()) + ")" );
    return;
  }
  double worst = 0.;
  for( size_t f = 0; f < sel[j].size() && f < ref[jref].size(); ++f ){
    double const d = std::fabs( sel[j][f] - ref[jref][f] );
    double const s = std::fabs( ref[jref][f] );
    worst = std::max( worst, d / ( s > 1e-12? s: 1. ) );
  }
  std::ostringstream os;
  os << tag << "  col " << j << " == reference col " << jref
     << "   (worst rel " << std::scientific << std::setprecision(2) << worst << ")";
  check( worst < 1e-6, os.str() );
}

void run_case( int reg, std::string const& name, std::vector<size_t> const& expect, Run const& ref )
{
  std::cout << "\n--- " << name << " : expected directions {";
  for( size_t i = 0; i < expect.size(); ++i ) std::cout << ( i? ",": "" ) << expect[i];
  std::cout << "}\n";
  Run const R = build_and_solve( reg );
  if( !R.ok ){ check( false, name + "  solve" ); return; }

  std::ostringstream on;
  on << name << "  n_control_dof()=" << R.ncd << " (expected " << expect.size() << ")";
  check( R.ncd == expect.size(), on.str() );

  std::ostringstream of;
  of << name << "  FSA returned " << R.fsa.size() << " directions (expected " << expect.size()
     << ", np=" << R.np << ")";
  check( R.fsa.size() == expect.size(), of.str() );

  std::ostringstream oa;
  oa << name << "  ASA returned " << R.asa.size() << " directions (expected " << expect.size() << ")";
  check( R.asa.size() == expect.size(), oa.str() );

  for( size_t j = 0; j < expect.size(); ++j ){
    check_col( R.fsa, j, ref.fsa, expect[j], name + " FSA" );
    check_col( R.asa, j, ref.asa, expect[j], name + " ASA" );
  }
}

int main()
{
  std::cout << std::scientific << std::setprecision(6);
  std::cout << "================================================================\n"
            << "  test3_ff: selective sensitivity analysis in ODESLV\n"
            << "  _mP order: 0=TF  1=x10  2,3,4=u(t)[0..2]\n"
            << "================================================================\n";

  std::cout << "\n--- A: no controls registered (reference; must be today's behaviour)\n";
  Run const ref = build_and_solve( REG_NONE );
  if( !ref.ok ){ std::cerr << "reference run FAILED\n"; return 1; }
  check( ref.np == 5, "A  np == 5" );
  check( ref.fsa.size() == 5, "A  FSA returns all 5 directions with an empty registry" );
  check( ref.asa.size() == 5, "A  ASA returns all 5 directions with an empty registry" );
  for( size_t j = 0; j < ref.fsa.size() && j < ref.asa.size(); ++j )
    check_col( ref.asa, j, ref.fsa, j, "A  ASA==FSA" );

  run_case( REG_U,          "B: u(t) only",        { 2, 3, 4 },    ref );
  run_case( REG_TF|REG_U,   "C: TF + u(t)",        { 0, 2, 3, 4 }, ref );   // NON-CONTIGUOUS
  run_case( REG_X10,        "D: x10 only",         { 1 },          ref );

  // ---- extraction: control_columns() on ODESLV's gradient, as on OCFESLV's sens_jacobian --------
  std::cout << "\n--- E: FFModel::control_columns() applied to ODESLV's function gradient\n";
  {
    mc::FFGraph DAG;
    mc::ODESLVS_CVODES IVP( &DAG );
    mc::FFVar t  = DAG.add_var( "t" );
    mc::FFVar X0 = DAG.add_var( "x0(t)" ), X1 = DAG.add_var( "x1(t)" );
    mc::FFVar TF = DAG.add_var( "TF" ), X10 = DAG.add_var( "x10" );
    mc::FFVar U  = DAG.add_var( "u(t)" );
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
    IVP.options.LINSOL  = mc::BASE_CVODES::Options::SPARSE;
    IVP.options.DISPLAY = 0;
    IVP.options.ATOL = IVP.options.ATOLB = IVP.options.ATOLS = 1e-9;
    IVP.options.RTOL = IVP.options.RTOLB = IVP.options.RTOLS = 1e-9;
    IVP.register_control( TF );
    IVP.register_control( U );
    if( !IVP.setup() ){ std::cerr << "  setup failed\n"; return 1; }
    if( IVP.solve_fsens( P0 ) == mc::ODESLVS_CVODES::STATUS::NORMAL ){
      auto const& g = IVP.val_function_gradient();          // [direction][function]
      size_t const ncd = IVP.n_control_dof(), nf = IVP.nf();
      std::ostringstream os;  os << "E  n_control_dof()=" << ncd << " and " << g.size() << " directions";
      check( ncd == g.size() && ncd > 0, os.str() );
      if( ncd == g.size() && ncd ){
        std::vector<double> Jrm( nf * ncd, 0. );            // row-major nf x ncd, as sens_jacobian()
        for( size_t j = 0; j < ncd; ++j )
          for( size_t f = 0; f < nf && f < g[j].size(); ++f ) Jrm[f*ncd + j] = g[j][f];
        auto const bTF = IVP.control_columns( TF, Jrm, nf );
        auto const bU  = IVP.control_columns( U,  Jrm, nf );
        check( bTF.size() == nf * 1, "E  control_columns(TF) has nf x 1 entries" );
        check( bU.size()  == nf * 3, "E  control_columns(u)  has nf x 3 entries" );
        bool okTF = true, okU = true;
        for( size_t f = 0; f < nf; ++f ){
          if( bTF.size() > f && std::fabs( bTF[f] - Jrm[f*ncd + 0] ) > 0. ) okTF = false;
          for( size_t k = 0; k < 3; ++k )
            if( bU.size() > f*3+k && std::fabs( bU[f*3+k] - Jrm[f*ncd + 1 + k] ) > 0. ) okU = false;
        }
        check( okTF, "E  control_columns(TF) == the TF block of the gradient" );
        check( okU,  "E  control_columns(u)  == the u block of the gradient" );
      }
    }
    else check( false, "E  solve_sensitivity" );
  }

  std::cout << "\n================================================================\n"
            << "  test3_ff: " << npass << " passed, " << nfail << " failed -- "
            << ( nfail? "FAILURES": "ALL PASS" ) << "\n"
            << "================================================================\n";
  return nfail? 1: 0;
}
