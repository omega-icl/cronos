// test5_ff.cpp -- OCFE_DAE3q's model and ORACLE, ported to ODESLV.
// =============================================================================================
// Recovered from the OCEnv-era driver OCFE_DAE3q.cpp (thread "Test examples for ODE and DAE
// models", 2026-08-02).  That driver is written against ocenv.hpp / OCEnv::Options and does not
// compile here; what is ported is its MODEL and its ANALYTIC RECURRENCE, which is the part worth
// having -- an oracle derived independently of the interpolant under test.
//
//     dc/dt = -a c + u(t),     c(0) = c_ic,      u piecewise QUADRATIC (n_node = 3, LGR)
//
// u is DISCONTINUOUS at every interior element interface by construction, which is exactly the
// case step 3 must handle: the control is declared per element, so the interpolant is rebuilt per
// element and nothing is shared across the interface.
//
// WHY THIS IS A STRONGER GATE THAN test4.  test4 checks the interpolant against its own degenerate
// case -- equal nodal values must give the constant.  That catches a mis-normalised basis but not a
// basis evaluated at the WRONG NODES, because with equal values every node gives the same answer.
// Here u is exactly quadratic on each element, so an n_node=3 interpolant reproduces it EXACTLY,
// and the closed-form recurrence below is then the exact solution of the DISCRETISED problem.  If
// the node positions or the basis are wrong, the reconstructed u is a different quadratic and c(t_m)
// misses the recurrence.  The driver supplies NODAL VALUES u(s_j), not the coefficients A,B,C --
// the solver has to rebuild the polynomial from them.
//
//   Particular quadratic for dc/dt = -a c + (A + B tau + C tau^2):
//     p2 = C/a,   p1 = B/a - 2C/a^2,   p0 = A/a - B/a^2 + 2C/a^3
//     c(h) = (c_k - p0) e^{-a h} + p0 + p1 h + p2 h^2

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <cmath>

#include "odeslvs_cvodes.hpp"

// --- OCFE_DAE3q's parameters, verbatim -------------------------------------------------------
static double const a_decay = 1.0;
static double const c0_ic   = 0.5;
static size_t const Nel     = 5;                 // elements on [0,Nel], width 1
static double const T_end   = double( Nel );
static size_t const u_nd    = 3;                 // INPUT nodes (LGR) -- piecewise quadratic

static double const Acoef[5] = { 1.0,  1.4,  0.8,  1.3,  0.7 };
static double const Bcoef[5] = { 0.6, -0.4,  0.7, -0.3,  0.5 };
static double const Ccoef[5] = { 0.3,  0.25,-0.20, 0.35,-0.15 };

//! @brief The intended control: exactly quadratic on each element, jumping at every interface.
static double u_fun( size_t k, double tau )
{ return Acoef[k] + Bcoef[k]*tau + Ccoef[k]*tau*tau; }

//! @brief Closed-form c at each element end -- the exact solution of the discretised problem.
static std::vector<double> recurrence()
{
  std::vector<double> c( Nel, 0. );
  double ck = c0_ic;
  double const h = 1., a = a_decay;
  for( size_t k = 0; k < Nel; ++k ){
    double const p2 = Ccoef[k]/a;
    double const p1 = Bcoef[k]/a - 2.*Ccoef[k]/(a*a);
    double const p0 = Acoef[k]/a - Bcoef[k]/(a*a) + 2.*Ccoef[k]/(a*a*a);
    ck = ( ck - p0 )*std::exp( -a*h ) + p0 + p1*h + p2*h*h;
    c[k] = ck;
  }
  return c;
}

static int npass = 0, nfail = 0;
void check( bool c, std::string const& what )
{ std::cout << "  " << (c? "PASS  ": "FAIL  ") << what << std::endl; c? ++npass: ++nfail; }

struct Out
{
  bool ok = false;  std::string err;
  size_t np = 0, ncd = 0, nf = 0;
  std::vector<double> f;
  std::vector<std::vector<double>> g;
};

//! @brief Build the model, carrying @p U as the nodal values of the control (element-major).
Out solve( std::vector<double> const& U, bool adjoint = false )
{
  Out R;
  mc::FFGraph DAG;
  mc::ODESLVS_CVODES IVP( &DAG );
  mc::FFVar t  = DAG.add_var( "t" );
  mc::FFVar C  = DAG.add_var( "c(t)" );
  mc::FFVar CIC = DAG.add_var( "c_ic" );
  mc::FFVar Uv = DAG.add_var( "u(t)" );
  mc::FFPartial OpP;  mc::FFEval OpEval;

  IVP.add_domain( t, mc::FFDom( 0., T_end, Nel, mc::FFDom::LGL, 7 ) );
  IVP.add_state( C, {t} );
  IVP.add_input( CIC );
  IVP.add_input( Uv, { t }, mc::FFDom::LGR, u_nd );      // the declaration under test
  IVP.set_evolution_domain( t );
  IVP.update_ref( C, c0_ic );

  int const T_INT = mc::FFDom::ALL - mc::FFDom::LB;
  auto const io = mc::FFModel::EqnOptions( mc::FFModel::EqnRole::INTERIOR, 0 );
  auto const ii = mc::FFModel::EqnOptions( mc::FFModel::EqnRole::INITIAL,  0 );
  IVP.add_equation( OpP( C, t ) + a_decay * C - Uv, {t}, {T_INT},         io );
  IVP.add_equation( C - CIC,                        {t}, {mc::FFDom::LB}, ii );

  for( size_t k = 1; k <= Nel; ++k )                     // c at every element end
    IVP.add_output( OpEval( C, t, double(k) ) );

  IVP.options.LINSOL  = mc::BASE_CVODES::Options::SPARSE;
  IVP.options.DISPLAY = 0;
  IVP.options.ATOL = IVP.options.ATOLB = IVP.options.ATOLS = 1e-11;
  IVP.options.RTOL = IVP.options.RTOLB = IVP.options.RTOLS = 1e-11;

  IVP.register_control( Uv );
  if( !IVP.setup() ){ R.err = IVP.extract_error(); return R; }
  R.np = IVP.np();  R.nf = IVP.nf();  R.ncd = IVP.n_control_dof();

  std::vector<double> P( R.np, 0. );
  P[0] = c0_ic;                                          // _mP order: c_ic, then the nodal values
  if( U.size() + 1 != R.np ){
    R.err = "driver: " + std::to_string(U.size()) + " nodal values for np=" + std::to_string(R.np);
    return R;
  }
  for( size_t i = 0; i < U.size(); ++i ) P[1+i] = U[i];

  auto const st = adjoint? IVP.solve_asens( P ): IVP.solve_fsens( P );
  if( st != mc::ODESLVS_CVODES::STATUS::NORMAL ){ R.err = "solve failed"; return R; }
  R.f = IVP.val_function();
  R.g = IVP.val_function_gradient();
  R.ok = true;
  return R;
}

double reldiff( double a, double b )
{ double const s = std::fabs(b); return std::fabs(a-b) / ( s > 1e-12? s: 1. ); }

int main()
{
  std::cout << std::scientific << std::setprecision(6);
  std::cout << "================================================================\n"
            << "  test5_ff: OCFE_DAE3q's piecewise-QUADRATIC control, on ODESLV\n"
            << "  dc/dt = -a c + u(t),  u exactly quadratic per element (LGR, n_node=3)\n"
            << "================================================================\n";

  // The nodal values the driver hands over: u evaluated at the element-local LGR nodes.  The node
  // POSITIONS are the declaration (same BASE_OC family the solver reads); the BASIS that rebuilds
  // the quadratic from them is what is under test.
  mc::BASE_OC quad;
  if( !quad.set_lgnodes( (int)mc::FFDom::LGR, u_nd, -1., 1. ) ){
    std::cout << "  cannot generate LGR nodes\n"; return 1; }
  std::vector<double> const s = quad.lgnodes( 0., 1. );   // element width is 1, so tau == s
  std::cout << "\n  element-local LGR nodes:";
  for( auto const& v : s ) std::cout << " " << v;
  std::cout << "\n";

  std::vector<double> U;
  for( size_t e = 0; e < Nel; ++e )
    for( size_t j = 0; j < u_nd; ++j ) U.push_back( u_fun( e, s[j] ) );

  std::cout << "\n--- A: the model is accepted and sized from the declaration\n";
  Out const A = solve( U );
  if( !A.ok ){ check( false, std::string("A  setup/solve: ") + A.err );
               std::cout << "\n  test5_ff: " << npass << " passed, " << ++nfail << " failed\n"; return 1; }
  check( A.np  == 1 + Nel*u_nd, "A  np == 1 + Nel*u_nd" );
  check( A.ncd == Nel*u_nd,     "A  n_control_dof == Nel*u_nd == 15" );
  check( A.nf  == Nel,          "A  one output per element end" );

  std::cout << "\n--- B: c(t_m) against the closed-form recurrence  [the independent oracle]\n";
  std::vector<double> const cref = recurrence();
  double worst = 0.;
  for( size_t k = 0; k < Nel && k < A.f.size(); ++k ){
    double const d = reldiff( A.f[k], cref[k] );
    worst = std::max( worst, d );
    std::ostringstream os;
    os << "B  c(" << (k+1) << ") = " << A.f[k] << " vs recurrence " << cref[k]
       << "   (rel " << d << ")";
    check( d < 1e-6, os.str() );
  }
  std::ostringstream ow;  ow << "B  worst relative deviation over all element ends = " << worst;
  check( worst < 1e-6, ow.str() );

  std::cout << "\n--- C: adjoint agrees with forward sensitivity on all 15 nodal directions\n";
  Out const Cc = solve( U, true );
  if( !Cc.ok ) check( false, std::string("C  ") + Cc.err );
  else{
    double w = 0.;  bool ok = ( Cc.g.size() == A.g.size() );
    for( size_t j = 0; ok && j < Cc.g.size(); ++j )
      for( size_t k = 0; k < Cc.g[j].size(); ++k )
        w = std::max( w, reldiff( Cc.g[j][k], A.g[j][k] ) );
    std::ostringstream os;  os << "C  ASA == FSA (worst rel " << w << ")";
    check( ok && w < 1e-5, os.str() );
  }

  std::cout << "\n--- D: nodal gradients against central finite differences\n";
  {
    double const h = 1e-6;  double w = 0.;  bool ok = true;
    for( size_t d = 0; d < U.size() && ok; ++d ){
      std::vector<double> Up = U, Um = U;
      Up[d] += h;  Um[d] -= h;
      Out const P = solve( Up ), M = solve( Um );
      if( !P.ok || !M.ok ){ ok = false; break; }
      for( size_t k = 0; k < P.f.size() && k < A.g[d].size(); ++k ){
        double const fd = ( P.f[k] - M.f[k] ) / ( 2.*h );
        w = std::max( w, std::fabs( fd - A.g[d][k] ) / ( std::fabs(fd) > 1e-8? std::fabs(fd): 1. ) );
      }
    }
    std::ostringstream os;  os << "D  FSA == central FD on all " << U.size()
                               << " nodal directions (worst rel " << w << ")";
    check( ok && w < 1e-4, os.str() );
  }

  std::cout << "\n--- E: causality -- a nodal value in element e cannot affect c before element e\n";
  {
    bool ok = true;  double leak = 0.;
    for( size_t e = 0; e < Nel; ++e )
      for( size_t j = 0; j < u_nd; ++j ){
        size_t const d = e*u_nd + j;
        for( size_t k = 0; k + 1 <= e && k < A.g[d].size(); ++k )   // outputs c(1..e), i.e. before e ends
          if( std::fabs( A.g[d][k] ) > 1e-8 ){ ok = false; leak = std::max( leak, std::fabs(A.g[d][k]) ); }
      }
    std::ostringstream os;  os << "E  no upstream leakage (worst |dc/dU| before the element = " << leak << ")";
    check( ok, os.str() );
  }

  std::cout << "\n================================================================\n"
            << "  test5_ff: " << npass << " passed, " << nfail << " failed -- "
            << ( nfail? "FAILURES": "ALL PASS" ) << "\n"
            << "================================================================\n";
  return nfail? 1: 0;
}
