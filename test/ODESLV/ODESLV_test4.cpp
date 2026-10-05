// test4_ff.cpp -- GATE for n_node > 1 controls in ODESLV (step 3).
// =============================================================================================
// ODESLV refuses a time-varying input declared with more than one value per element:
//
//    "an integrator needs its interpolant evaluated at t, which is not implemented --
//     declare it piecewise constant, or hold it fixed"                 (odeslv_base_ff.hpp)
//
// Step 3 lifts that.  On element e the control becomes u_e(t) = sum_j L_j(s) U[e,j] with
// s = (t - t_e)/h_e the element-local coordinate and L_j the Lagrange basis on the declared
// node family -- so U[e,j] are n_elem * n_node parameters, which is exactly the DOF count
// FFModel::control_ndof() already reports.
//
// THE GATE IS AN EXACT EQUIVALENCE, NOT A CONVERGENCE CHECK.  A Lagrange interpolant through
// n_node EQUAL nodal values is the constant, whatever the node family.  So an n_node=2 control
// with U[e,0] == U[e,1] == c_e must reproduce, to integrator tolerance, the n_node=1 control
// with level c_e -- same states, same outputs, same gradients w.r.t. c_e once the two nodal
// derivatives are summed (dJ/dc = dJ/dU[e,0] + dJ/dU[e,1] by the chain rule).
//
// Case F is the sensitivity check that a wrong interpolant would survive: a NON-constant
// control, with the gradient w.r.t. each nodal value checked against central finite
// differences on that value.  FD shares the RHS, so it cannot catch a wrong interpolant --
// that is case B's job -- but it does catch a wrong slot, a wrong seed or a wrong stride.
//
// BEFORE step 3, B..F must all FAIL: setup() refuses the model outright.

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <cmath>

#include "odeslvs_cvodes.hpp"

size_t const NS = 3;
double const t0 = 0., tf = 1.;
double const TF_VAL = 6., X10_VAL = 0.5;

static int npass = 0, nfail = 0;
void check( bool c, std::string const& what )
{ std::cout << "  " << (c? "PASS  ": "FAIL  ") << what << std::endl; c? ++npass: ++nfail; }

struct Out
{
  bool   ok = false;
  size_t np = 0, ncd = 0, nf = 0;
  std::vector<double> f;                        // output values
  std::vector<std::vector<double>> g;           // [direction][function]
  std::string err;
};

//! @brief test2's model with the control declared at @p n_node values per element, carrying the
//! nodal values @p U (n_elem * n_node of them, element-major).  Registers u as the only control.
Out solve( size_t n_node, std::vector<double> const& U, bool adjoint = false )
{
  Out R;
  mc::FFGraph DAG;
  mc::ODESLVS_CVODES IVP( &DAG );
  mc::FFVar t  = DAG.add_var( "t" );
  mc::FFVar X0 = DAG.add_var( "x0(t)" ), X1 = DAG.add_var( "x1(t)" );
  mc::FFVar TF = DAG.add_var( "TF" ), X10 = DAG.add_var( "x10" );
  mc::FFVar Uv = DAG.add_var( "u(t)" );
  mc::FFPartial OpP;  mc::FFEval OpEval;

  IVP.add_domain( t, mc::FFDom( t0, tf, NS, mc::FFDom::LGR, 4 ) );
  IVP.add_state( X0, {t} );  IVP.add_state( X1, {t} );
  IVP.add_input( TF );  IVP.add_input( X10 );
  IVP.add_input( Uv, { t }, mc::FFDom::LGL, n_node );      // LGL: nodes include both element ends
  IVP.set_evolution_domain( t );
  IVP.update_ref( X0, 0. );  IVP.update_ref( X1, 0.5 );

  int const T_INT = mc::FFDom::ALL - mc::FFDom::LB;
  auto const io = mc::FFModel::EqnOptions( mc::FFModel::EqnRole::INTERIOR, 0 );
  auto const ii = mc::FFModel::EqnOptions( mc::FFModel::EqnRole::INITIAL,  0 );
  IVP.add_equation( OpP( X0, t ) - TF * X1,                        {t}, {T_INT},         io );
  IVP.add_equation( OpP( X1, t ) - TF * ( Uv * X0 - 2. * X1 ) * t, {t}, {T_INT},         io );
  IVP.add_equation( X0 - 0.,                                       {t}, {mc::FFDom::LB}, ii );
  IVP.add_equation( X1 - X10,                                      {t}, {mc::FFDom::LB}, ii );
  IVP.add_output( OpEval( X0 - 1., t, tf ) );
  IVP.add_output( OpEval( X1,      t, tf ) );

  IVP.options.LINSOL  = mc::BASE_CVODES::Options::SPARSE;
  IVP.options.DISPLAY = 0;
  IVP.options.ATOL = IVP.options.ATOLB = IVP.options.ATOLS = 1e-10;
  IVP.options.RTOL = IVP.options.RTOLB = IVP.options.RTOLS = 1e-10;

  IVP.register_control( Uv );
  if( !IVP.setup() ){ R.err = IVP.extract_error(); return R; }
  R.np = IVP.np();  R.nf = IVP.nf();  R.ncd = IVP.n_control_dof();

  std::vector<double> P( R.np, 0. );                 // _mP order: TF, x10, then the nodal values
  P[0] = TF_VAL;  P[1] = X10_VAL;
  if( U.size() + 2 != R.np ){
    R.err = "driver: supplied " + std::to_string(U.size()) + " nodal values for np="
          + std::to_string(R.np);
    return R;
  }
  for( size_t i = 0; i < U.size(); ++i ) P[2+i] = U[i];

  auto const st = adjoint? IVP.solve_asens( P ): IVP.solve_fsens( P );
  if( st != mc::ODESLVS_CVODES::STATUS::NORMAL ){ R.err = "solve failed"; return R; }
  R.g = IVP.val_function_gradient();
  R.f = IVP.val_function();
  R.ok = true;
  return R;
}

double reldiff( double a, double b )
{ double const s = std::fabs(b); return std::fabs(a-b) / ( s > 1e-12? s: 1. ); }

int main()
{
  std::cout << std::scientific << std::setprecision(6);
  std::cout << "================================================================\n"
            << "  test4_ff: n_node > 1 controls in ODESLV\n"
            << "================================================================\n";

  std::vector<double> const C = { 0.4, 0.7, 0.2 };            // one level per element

  std::cout << "\n--- A: n_node = 1 (piecewise constant) -- reference\n";
  Out const A = solve( 1, C );
  if( !A.ok ){ std::cout << "  reference FAILED: " << A.err << "\n"; return 1; }
  check( A.np == 2 + NS,     "A  np == 2 + NS" );
  check( A.ncd == NS,        "A  n_control_dof == NS" );
  check( A.g.size() == NS,   "A  NS sensitivity directions" );

  std::cout << "\n--- B: n_node = 2, both nodal values equal -- must reproduce A exactly\n";
  std::vector<double> U2;                                     // element-major: U[e,0], U[e,1]
  for( size_t e = 0; e < NS; ++e ){ U2.push_back( C[e] ); U2.push_back( C[e] ); }
  Out const B = solve( 2, U2 );
  if( !B.ok ){
    check( false, std::string("B  setup/solve: ") + B.err );
    std::cout << "\n  (B..F all depend on n_node>1 being accepted)\n";
  }
  else{
    check( B.np  == 2 + 2*NS, "B  np == 2 + 2*NS  (one parameter per element NODE)" );
    check( B.ncd == 2*NS,     "B  n_control_dof == 2*NS" );
    check( B.g.size() == 2*NS,"B  2*NS sensitivity directions" );
    bool okf = ( B.f.size() == A.f.size() );
    for( size_t i = 0; okf && i < A.f.size(); ++i ) okf = reldiff( B.f[i], A.f[i] ) < 1e-6;
    check( okf, "B  output values equal the piecewise-constant reference" );
    // chain rule: dJ/dc_e = dJ/dU[e,0] + dJ/dU[e,1]
    bool okg = ( B.g.size() == 2*NS && A.g.size() == NS );
    double worst = 0.;
    for( size_t e = 0; okg && e < NS; ++e )
      for( size_t k = 0; k < A.g[e].size(); ++k ){
        double const s = B.g[2*e][k] + B.g[2*e+1][k];
        worst = std::max( worst, reldiff( s, A.g[e][k] ) );
      }
    std::ostringstream os;
    os << "B  dJ/dU[e,0] + dJ/dU[e,1] == dJ/dc_e  (worst rel " << worst << ")";
    check( okg && worst < 1e-5, os.str() );
  }

  std::cout << "\n--- C: n_node = 2, UNEQUAL nodal values -- must NOT match A (test sensitivity)\n";
  std::vector<double> U2b;
  for( size_t e = 0; e < NS; ++e ){ U2b.push_back( C[e] - 0.15 ); U2b.push_back( C[e] + 0.15 ); }
  Out const Cc = solve( 2, U2b );
  if( !Cc.ok ) check( false, std::string("C  setup/solve: ") + Cc.err );
  else{
    bool differs = false;
    for( size_t i = 0; i < Cc.f.size() && i < A.f.size(); ++i )
      if( reldiff( Cc.f[i], A.f[i] ) > 1e-4 ) differs = true;
    check( differs, "C  a genuinely varying control changes the outputs" );
  }

  std::cout << "\n--- D: n_node = 2 -- adjoint agrees with forward sensitivity\n";
  Out const D = solve( 2, U2b, true );
  if( !D.ok ) check( false, std::string("D  setup/solve: ") + D.err );
  else{
    bool ok = ( D.g.size() == Cc.g.size() );
    double worst = 0.;
    for( size_t j = 0; ok && j < D.g.size(); ++j )
      for( size_t k = 0; k < D.g[j].size(); ++k )
        worst = std::max( worst, reldiff( D.g[j][k], Cc.g[j][k] ) );
    std::ostringstream os;  os << "D  ASA == FSA on every nodal direction (worst rel " << worst << ")";
    check( ok && worst < 1e-5, os.str() );
  }

  std::cout << "\n--- E: n_node = 3 -- DOF count follows the declaration\n";
  std::vector<double> U3;
  for( size_t e = 0; e < NS; ++e ){ U3.push_back( C[e] ); U3.push_back( C[e] ); U3.push_back( C[e] ); }
  Out const E = solve( 3, U3 );
  if( !E.ok ) check( false, std::string("E  setup/solve: ") + E.err );
  else{
    check( E.np == 2 + 3*NS,  "E  np == 2 + 3*NS" );
    check( E.ncd == 3*NS,     "E  n_control_dof == 3*NS" );
    bool okf = ( E.f.size() == A.f.size() );
    for( size_t i = 0; okf && i < A.f.size(); ++i ) okf = reldiff( E.f[i], A.f[i] ) < 1e-6;
    check( okf, "E  three equal nodal values still reproduce the constant control" );
  }

  std::cout << "\n--- F: n_node = 2 -- nodal gradients against central finite differences\n";
  if( !Cc.ok ) check( false, "F  needs case C" );
  else{
    double const h = 1e-5;
    double worst = 0.;  bool ok = true;
    for( size_t d = 0; d < U2b.size() && ok; ++d ){
      std::vector<double> Up = U2b, Um = U2b;
      Up[d] += h;  Um[d] -= h;
      Out const P = solve( 2, Up ), M = solve( 2, Um );
      if( !P.ok || !M.ok ){ ok = false; break; }
      for( size_t k = 0; k < P.f.size() && k < Cc.g[d].size(); ++k ){
        double const fd = ( P.f[k] - M.f[k] ) / ( 2.*h );
        worst = std::max( worst, std::fabs( fd - Cc.g[d][k] )
                                 / ( std::fabs(fd) > 1e-8? std::fabs(fd): 1. ) );
      }
    }
    std::ostringstream os;  os << "F  FSA == central FD on all " << U2b.size()
                               << " nodal directions (worst rel " << worst << ")";
    check( ok && worst < 1e-4, os.str() );
  }

  std::cout << "\n================================================================\n"
            << "  test4_ff: " << npass << " passed, " << nfail << " failed -- "
            << ( nfail? "FAILURES": "ALL PASS" ) << "\n"
            << "================================================================\n";
  return nfail? 1: 0;
}
