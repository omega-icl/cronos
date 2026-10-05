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

int main()
{
  std::cout << std::scientific << std::setprecision(6);
  std::cout << "================================================================\n"
            << "  test7_ff: fdiff (symbolic augmentation) vs CVODES forward sensitivity\n"
            << "================================================================\n";
  mc::FFGraph DAG;
  mc::ODESLVS_CVODES IVP( &DAG );
  mc::FFVar t = DAG.add_var( "t" ), X0 = DAG.add_var( "x0(t)" ), X1 = DAG.add_var( "x1(t)" );
  mc::FFVar TF = DAG.add_var( "TF" ), X10 = DAG.add_var( "x10" ), U = DAG.add_var( "u(t)" );
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
  IVP.add_equation( OpP( X0, t ) - TF * X1,                       {t}, {T_INT},         io );
  IVP.add_equation( OpP( X1, t ) - TF * ( U * X0 - 2. * X1 ) * t, {t}, {T_INT},         io );
  IVP.add_equation( X0 - 0.,                                      {t}, {mc::FFDom::LB}, ii );
  IVP.add_equation( X1 - X10,                                     {t}, {mc::FFDom::LB}, ii );
  IVP.add_output( OpEval( X0 - 1., t, tf ) );
  IVP.add_output( OpEval( X1,      t, tf ) );
  IVP.options.LINSOL = mc::BASE_CVODES::Options::DENSE;  IVP.options.DISPLAY = 0;
  IVP.options.ATOL = IVP.options.ATOLS = 1e-10;  IVP.options.RTOL = IVP.options.RTOLS = 1e-10;
  if( !IVP.setup() ){ std::cout << "  setup failed\n"; return 1; }

  std::cout << "\n--- CVODES forward sensitivity on the original problem\n";
  check( IVP.solve_fsens( P0 ) == mc::ODESLVS_CVODES::STATUS::NORMAL, "FSA solve" );
  auto const G = IVP.val_function_gradient();          // [direction][function]
  size_t const np = IVP.np(), nf = IVP.nf();
  check( G.size() == np, "FSA returns np directions" );

  std::cout << "\n--- fdiff: the augmented system, differentiated w.r.t. all " << np << " parameters\n";
  check( IVP.var_parameter().size() == np, "var_parameter() exposes np parameters, minted levels included" );
  check( IVP.sensitivity_index().size() == np, "no controls registered -> every parameter is a direction" );
  std::vector<std::vector<mc::FFVar>> vDU, vS; std::string err;
  std::vector<mc::FFVar> vAll = IVP.var_input().size()? std::vector<mc::FFVar>{} : std::vector<mc::FFVar>{};
  for( auto const& [w,d] : IVP.var_declared_input() ) vAll.push_back( w );
  mc::ODESLVS_CVODES* SEN = IVP.fdiff( vAll, 0, vDU, vS, err );    // one direction per DOF of every input
  check( SEN != nullptr, std::string("fdiff returned a product") + ( err.empty()? "": "  -- " + err ) );
  if( !SEN ) return 1;
  bool const ok_setup = SEN->setup();
  check( ok_setup, "setup() on the fdiff product (description path -- failed with UNDEFINED DAG before today)" );
  if( ok_setup ){
    check( SEN->nx() == (1+np)*IVP.nx(), "product has nx*(1+np) states" );
    check( SEN->nf() == (1+np)*nf,       "product has (1+np)*nf functions: F then the np gradient blocks" );
    // the product's parameters: originals by name, then every direction seeded as its own DOF unit vector
    std::vector<double> P( SEN->np(), 0. );
    { auto const& vP = SEN->var_parameter();
      for( size_t i = 0; i < vP.size(); ++i ){ std::string const n = vP[i].name();
        if( n=="TF" ) P[i]=P0[0]; else if( n=="x10" ) P[i]=P0[1];
        else if( n.rfind("u(t)",0)==0 ) P[i]=P0[2+std::stoi(n.substr(n.find('[')+1))]; } }
    bool okseed = true;
    for( size_t k = 0; k < vDU.size(); ++k ) okseed = okseed && mc::FFModel::fdiff_seed( *SEN, vDU, k, P );
    // fdiff_seed zeroes all directions each call: seed all of them for the ONE-SOLVE Jacobian instead
    for( size_t k = 0; k < vDU.size(); ++k ){ std::vector<double> Pk( SEN->np(), 0. ); mc::FFModel::fdiff_seed( *SEN, vDU, k, Pk );
      for( size_t i = 0; i < Pk.size(); ++i ) if( Pk[i] != 0. ) P[i] = Pk[i]; }
    check( okseed, "every direction has a DOF to seed" );
    check( SEN->solve( P ) == mc::ODESLVS_CVODES::STATUS::NORMAL, "ONE solve_state on the product" );
    auto const& F = SEN->val_function();
    double worst = 0.;  bool ok = ( F.size() == (1+np)*nf );
    for( size_t i = 0; ok && i < np; ++i )
      for( size_t j = 0; j < nf; ++j ){
        double const a = F[nf + i*nf + j], b = G[i][j];
        worst = std::max( worst, std::fabs(a-b)/( std::fabs(b)>1e-12? std::fabs(b): 1. ) );
      }
    std::ostringstream os;  os << "fdiff functions == FSA gradients, all " << np*nf
                               << " entries (worst rel " << worst << ")";
    check( ok && worst < 1e-6, os.str() );
  }
  delete SEN;

  std::cout << "\n--- subset: register u(t) only -- fdiff() must follow the control registry\n";
  {
    mc::FFGraph DAG2;  mc::ODESLVS_CVODES IVP2( &DAG2 );
    mc::FFVar t2 = DAG2.add_var("t"), Y0 = DAG2.add_var("x0(t)"), Y1 = DAG2.add_var("x1(t)");
    mc::FFVar TF2 = DAG2.add_var("TF"), X102 = DAG2.add_var("x10"), U2 = DAG2.add_var("u(t)");
    IVP2.add_domain( t2, mc::FFDom( t0, tf, NS, mc::FFDom::LGR, 4 ) );
    IVP2.add_state( Y0, {t2} );  IVP2.add_state( Y1, {t2} );
    IVP2.add_input( TF2 );  IVP2.add_input( X102 );  IVP2.add_input( U2, { t2 }, mc::FFDom::LGR, 1 );
    IVP2.set_evolution_domain( t2 );  IVP2.update_ref( Y0, 0. );  IVP2.update_ref( Y1, 0.5 );
    IVP2.add_equation( OpP( Y0, t2 ) - TF2 * Y1,                         {t2}, {T_INT},         io );
    IVP2.add_equation( OpP( Y1, t2 ) - TF2 * ( U2 * Y0 - 2. * Y1 ) * t2, {t2}, {T_INT},         io );
    IVP2.add_equation( Y0 - 0.,                                          {t2}, {mc::FFDom::LB}, ii );
    IVP2.add_equation( Y1 - X102,                                        {t2}, {mc::FFDom::LB}, ii );
    IVP2.add_output( OpEval( Y0 - 1., t2, tf ) );  IVP2.add_output( OpEval( Y1, t2, tf ) );
    IVP2.options = IVP.options;
    IVP2.register_control( U2 );
    if( !IVP2.setup() ){ check( false, "subset setup" ); }
    else{
      check( IVP2.sensitivity_index().size() == 3, "u(t) only -> 3 directions" );
      std::vector<std::vector<mc::FFVar>> vDU2, vS2; std::string e2;
      mc::ODESLVS_CVODES* S2 = IVP2.fdiff( vDU2, vS2, e2 );          // registered controls, one direction per DOF
      check( S2 && S2->setup() && S2->nf() == 4*nf, "product has (1+3)*nf functions: fdiff followed the registry" );
      if( S2 ){
        std::vector<double> P( S2->np(), 0. );
        { auto const& vP = S2->var_parameter();
          for( size_t i = 0; i < vP.size(); ++i ){ std::string const n = vP[i].name();
            if( n=="TF" ) P[i]=P0[0]; else if( n=="x10" ) P[i]=P0[1];
            else if( n.rfind("u(t)",0)==0 ) P[i]=P0[2+std::stoi(n.substr(n.find('[')+1))]; } }
        for( size_t k = 0; k < vDU2.size(); ++k ){ std::vector<double> Pk( S2->np(), 0. ); mc::FFModel::fdiff_seed( *S2, vDU2, k, Pk );
          for( size_t i = 0; i < Pk.size(); ++i ) if( Pk[i] != 0. ) P[i] = Pk[i]; }
        check( S2->solve( P ) == mc::ODESLVS_CVODES::STATUS::NORMAL, "subset product solves (one solve)" );
        auto const& F2 = S2->val_function();
        double worst = 0.;
        for( size_t i = 0; i < 3 && F2.size() == 4*nf; ++i )
          for( size_t j = 0; j < nf; ++j ){
            double const a = F2[nf + i*nf+j], b = G[2+i][j];      // reference columns 2,3,4 are u[0..2]
            worst = std::max( worst, std::fabs(a-b)/( std::fabs(b)>1e-12? std::fabs(b): 1. ) );
          }
        std::ostringstream os;  os << "subset fdiff == reference FSA columns {2,3,4} (worst rel " << worst << ")";
        check( F2.size() == 4*nf && worst < 1e-6, os.str() );
        delete S2;
      }
    }
  }
  std::cout << "\n  test7_ff: " << npass << " passed, " << nfail << " failed -- "
            << ( nfail? "FAILURES": "ALL PASS" ) << "\n";
  return nfail? 1: 0;
}
