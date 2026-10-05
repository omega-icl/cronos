// OCFE_PDE5b.cpp  ---  OUTPUT FUNCTIONALS + FORWARD/ADJOINT SENSITIVITY on the PDE5 model
// =============================================================================================
// Extends the manufactured dynamic 2D Laplace/heat model of OCFE_PDE5 with two output
// functionals and validates their VALUES and their DERIVATIVES w.r.t. two lumped controls that
// parameterise the initial Gaussian:
//
//   control A   = initial Gaussian magnitude   (peak amplitude at t=0)
//   control W0  = initial Gaussian spread       (S(0) = W0, the initial width)
//
//   output f0 = U(tf, x0, y0)               (bump peak at the final time)
//   output f1 = double integral over x,y of U(tf, . , .)   (spatial inventory at final time)
//
// Because the forcing, IC and boundary data are all manufactured from U_exact(.;A,W0), the solve
// reproduces U == U_exact(.;A,W0) for EVERY control value.  Hence both functionals are known in
// closed form, giving analytic reference derivatives:
//
//   S(tf) = W0 + spread*tf
//   f0 = bg + A*W0/S(tf)
//        df0/dA  = W0/S(tf)
//        df0/dW0 = A*spread*tf / S(tf)^2
//   f1 = bg + pi*A*W0*erf(0.5/sqrt(S(tf)))^2                     (unit square, centred bump)
//        df1/dA  = pi*W0*E^2                              E = erf(0.5/sqrt(S(tf)))
//        df1/dW0 = pi*A*E^2 + 2*pi*A*W0*E*dE/dW0
//                  dE/dW0 = -(0.5/(sqrt(pi)*S^1.5))*exp(-0.25/S)
//
// Validation ladder, monolithic AND marching:
//   (1) solve_fsens  vs  solve_asens   -- identical reduced Jacobian (machinery cross-check)
//   (2) sens_functions()               -- vs analytic f0,f1 (values)
//   (3) sens_jacobian()                -- vs analytic dF/d(control) AND central finite differences
//   (4) marching sens_jacobian()       -- vs monolithic sens_jacobian() (the reinit test)
//
// Under marching, (4) exercises the override-path terminal-sensitivity transfer between windows:
// window k's LB state sensitivity must equal window k-1's terminal sensitivity, propagated through
// the SAME _marchICRows override the primal uses.  Until that extension lands, the marching path is
// expected to fail (4) -- it is included so the fix can be validated by flipping a switch.
// =============================================================================================

#include <iostream>
#include <iomanip>
#include <vector>
#include <array>
#include <cmath>
#include <string>

#ifndef OCFE_OCFESLV_HEADER
#define OCFE_OCFESLV_HEADER "ocfeslv.hpp"
#endif
#include OCFE_OCFESLV_HEADER

using namespace mc;

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifndef PDE5B_NEL_T
#define PDE5B_NEL_T 2
#endif
#ifndef PDE5B_NEL_X
#define PDE5B_NEL_X 1
#endif
#ifndef PDE5B_NEL_Y
#define PDE5B_NEL_Y 1
#endif
#ifndef PDE5B_NT
#define PDE5B_NT 8
#endif
#ifndef PDE5B_NX
#define PDE5B_NX 12
#endif
#ifndef PDE5B_NY
#define PDE5B_NY 12
#endif

namespace {

struct Par
{
  double tf     = 1.0;
  double kappa  = 3.0e-2;
  double bg     = 0.0;
  double spread = 8.0e-2;
  double x0     = 0.5;
  double y0     = 0.5;
  double A0     = 4.0e-1;   // nominal magnitude control
  double W0_0   = 1.0e-1;   // nominal spread control
};

// ---- exact manufactured profile, parameterised by the two controls (A,W0) ----
static double Soft( double t, double W0, Par const& p ){ return W0 + p.spread*t; }
static double R2( double x, double y, Par const& p )
{ double const dx=x-p.x0, dy=y-p.y0; return dx*dx+dy*dy; }
static double bump( double t, double x, double y, double A, double W0, Par const& p )
{ double const S=Soft(t,W0,p); return A*W0/S*std::exp(-R2(x,y,p)/S); }
static double U_exact( double t, double x, double y, double A, double W0, Par const& p )
{ return p.bg + bump(t,x,y,A,W0,p); }

// ---- closed-form output functionals and their control derivatives at t=tf ----
struct Fana { double f0, f1, df0dA, df0dW0, df1dA, df1dW0; };
static Fana analytic( double A, double W0, Par const& p )
{
  double const S = Soft(p.tf,W0,p);
  double const E = std::erf( 0.5/std::sqrt(S) );
  double const dEdW0 = -( 0.5/(std::sqrt(M_PI)*std::pow(S,1.5)) )*std::exp(-0.25/S);
  Fana o;
  o.f0     = p.bg + A*W0/S;
  o.df0dA  = W0/S;
  o.df0dW0 = A*p.spread*p.tf/(S*S);
  o.f1     = p.bg + M_PI*A*W0*E*E;
  o.df1dA  = M_PI*W0*E*E;
  o.df1dW0 = M_PI*A*E*E + 2.0*M_PI*A*W0*E*dEdW0;
  return o;
}

static void check_close( char const* name, double got, double want, double tol, bool& ok )
{
  double const err = std::abs(got-want);
  bool const pass = std::isfinite(err) && err <= tol;
  ok = ok && pass;
  std::cout << std::left << std::setw(42) << name
            << " got=" << std::right << std::scientific << std::setprecision(6) << std::setw(14) << got
            << " want=" << std::setw(14) << want
            << " |err|=" << std::setw(11) << err
            << " tol=" << tol << "  " << (pass?"PASS":"FAIL") << "\n";
}

// -------------------------------------------------------------------------------------------
// Build + solve the model at the given control values; return the two output-function values.
// Used both for the primal/sensitivity run and as the finite-difference reference.
// -------------------------------------------------------------------------------------------
static bool solve_values( bool march, double Aval, double W0val, Par const& p,
                          std::array<double,2>& fout )
{
  FFGraph DAG;
  FFVar t = DAG.add_var("t");
  FFVar x = DAG.add_var("x");
  FFVar y = DAG.add_var("y");
  FFVar U = DAG.add_var("U(t,x,y)");
  FFVar A = DAG.add_var("A");
  FFVar W0= DAG.add_var("W0");
  FFPartial OpP;
  FFIntegral OpI;

  FFVar DX = x-p.x0, DY = y-p.y0;
  FFVar S  = W0 + p.spread*t;
  FFVar RR = DX*DX + DY*DY;
  FFVar BE = A*W0*exp(-RR/S)/S;
  FFVar UE = p.bg + BE;
  FFVar FE = BE*p.spread*( -1.0/S + RR/(S*S) )
           - p.kappa*BE*( 4.0*RR/(S*S) - 4.0/S );
  FFVar PDE = OpP(U,t) - p.kappa*( OpP(U,{x,2}) + OpP(U,{y,2}) ) - FE;
  FFVar BC  = U - UE;

  OCFESLV oc(&DAG);
  oc.options.DISPLAY_LEVEL   = 0;
  oc.options.SOLVE.VERBOSE   = false;
  oc.options.SOLVE.MARCHING  = march;
  oc.options.INTERFACE.IMPOSITION = OCFESLV::Options::IC_STRONG;//WEAK;   // distributed-IC override path
  oc.options.INTERFACE.SAT_SIGMA0      = 100.0;                     // stiff interface continuity (PDE5 fix)
#if defined(CRONOS__WITH_SPQR)
  oc.options.SOLVE.FACTORIZATION = OCFESLV::Options::SOLVE_SPQR;
#endif

  oc.add_domain( t, FFDom(0., p.tf, PDE5B_NEL_T, FFDom::LGR, PDE5B_NT) );
  oc.add_domain( x, FFDom(0., 1.0,  PDE5B_NEL_X, FFDom::CGL, PDE5B_NX) );
  oc.add_domain( y, FFDom(0., 1.0,  PDE5B_NEL_Y, FFDom::CGL, PDE5B_NY) );
  oc.add_state( U, {t,x,y} );
  oc.add_input( A,  {} );          // lumped magnitude control
  oc.add_input( W0, {} );          // lumped spread control
  oc.update_ref( U, [&]( OCFESLV::t_Coord const& c ){ return U_exact(c.at(t),c.at(x),c.at(y),Aval,W0val,p); } );

  OCFESLV::EqnOptions int_opt( OCFESLV::EqnRole::INTERIOR, 0 );
  OCFESLV::EqnOptions ini_opt( OCFESLV::EqnRole::INITIAL,  0 );
  OCFESLV::EqnOptions bnd_opt( OCFESLV::EqnRole::BOUNDARY, 0 );
  int const T_INT = FFDom::ALL - FFDom::LB;
  int const X_INT = FFDom::ALL - FFDom::LB - FFDom::UB;
  int const Y_INT = FFDom::ALL - FFDom::LB - FFDom::UB;

  oc.add_equation( PDE, {t,x,y}, {T_INT, X_INT, Y_INT}, int_opt );
  oc.add_equation( BC,  {t,x,y}, {FFDom::LB, FFDom::ALL, FFDom::ALL}, ini_opt );
  oc.add_equation( BC,  {t,x,y}, {T_INT, FFDom::LB, Y_INT}, bnd_opt );
  oc.add_equation( BC,  {t,x,y}, {T_INT, FFDom::UB, Y_INT}, bnd_opt );
  oc.add_equation( BC,  {t,x,y}, {T_INT, FFDom::ALL, FFDom::LB}, bnd_opt );
  oc.add_equation( BC,  {t,x,y}, {T_INT, FFDom::ALL, FFDom::UB}, bnd_opt );
  oc.set_evolution_domain( t );

  // f0: peak value at (tf,x0,y0).  f1: double x,y integral of U at tf.
  oc.add_output( U, {t,x,y}, std::vector<double>{ p.tf, p.x0, p.y0 } );
  FFVar Iu = OpI( OpI( U, x ), y );
  oc.add_output( Iu, {t}, std::vector<double>{ p.tf } );

  if( !oc.setup() ){ std::cerr << "  setup FAILED\n"; return false; }

  std::vector<double> var, inp;
  if( !oc.init( var, inp, nullptr ) ){ std::cerr << "  init FAILED\n"; return false; }
  oc.set_input_values( A,  std::vector<double>{ Aval  }, inp.data() );
  oc.set_input_values( W0, std::vector<double>{ W0val }, inp.data() );

  OCFESLV::SolveReport const rep = oc.solve( var.data(), inp.data(), nullptr );
  if( !rep.converged ){ std::cerr << "  primal solve did NOT converge\n"; return false; }

  auto const& fv = oc.val_functions();
  if( fv.size() < 2 ){ std::cerr << "  val_functions() missing outputs\n"; return false; }
  fout = { fv[0], fv[1] };
  return true;
}

// -------------------------------------------------------------------------------------------
static bool run_mode( bool march, Par const& p )
{
  std::string const tag = march ? "MARCHING" : "MONOLITHIC";
  std::cout << "\n================ PDE5b output sensitivity : " << tag << " ================\n";

  FFGraph DAG;
  FFVar t = DAG.add_var("t");
  FFVar x = DAG.add_var("x");
  FFVar y = DAG.add_var("y");
  FFVar U = DAG.add_var("U(t,x,y)");
  FFVar A = DAG.add_var("A");
  FFVar W0= DAG.add_var("W0");
  FFPartial OpP;
  FFIntegral OpI;

  FFVar DX = x-p.x0, DY = y-p.y0;
  FFVar S  = W0 + p.spread*t;
  FFVar RR = DX*DX + DY*DY;
  FFVar BE = A*W0*exp(-RR/S)/S;
  FFVar UE = p.bg + BE;
  FFVar FE = BE*p.spread*( -1.0/S + RR/(S*S) )
           - p.kappa*BE*( 4.0*RR/(S*S) - 4.0/S );
  FFVar PDE = OpP(U,t) - p.kappa*( OpP(U,{x,2}) + OpP(U,{y,2}) ) - FE;
  FFVar BC  = U - UE;

  OCFESLV oc(&DAG);
  oc.options.DISPLAY_LEVEL   = 1;
  oc.options.SOLVE.VERBOSE   = true;//false;
  oc.options.SOLVE.MARCHING  = march;
  oc.options.INTERFACE.IMPOSITION = OCFESLV::Options::IC_WEAK;   // distributed-IC override path (nTrace=0)
  oc.options.INTERFACE.SAT_SIGMA0      = 100.0;                     // stiff interface continuity (PDE5 fix)
  oc.options.SOLVE.MAX_ITER  = 20;
  oc.options.SOLVE.RES_TOL   = 1e-9;

  oc.add_domain( t, FFDom(0., p.tf, PDE5B_NEL_T, FFDom::LGR, PDE5B_NT) );
  oc.add_domain( x, FFDom(0., 1.0,  PDE5B_NEL_X, FFDom::CGL, PDE5B_NX) );
  oc.add_domain( y, FFDom(0., 1.0,  PDE5B_NEL_Y, FFDom::CGL, PDE5B_NY) );
  oc.add_state( U, {t,x,y} );
  oc.add_input( A,  {} );
  oc.add_input( W0, {} );
  oc.update_ref( U, [&]( OCFESLV::t_Coord const& c ){ return U_exact(c.at(t),c.at(x),c.at(y),p.A0,p.W0_0,p); } );

  OCFESLV::EqnOptions int_opt( OCFESLV::EqnRole::INTERIOR, 0 );
  OCFESLV::EqnOptions ini_opt( OCFESLV::EqnRole::INITIAL,  0 );
  OCFESLV::EqnOptions bnd_opt( OCFESLV::EqnRole::BOUNDARY, 0 );
  int const T_INT = FFDom::ALL - FFDom::LB;
  int const X_INT = FFDom::ALL - FFDom::LB - FFDom::UB;
  int const Y_INT = FFDom::ALL - FFDom::LB - FFDom::UB;

  oc.add_equation( PDE, {t,x,y}, {T_INT, X_INT, Y_INT}, int_opt );
  oc.add_equation( BC,  {t,x,y}, {FFDom::LB, FFDom::ALL, FFDom::ALL}, ini_opt );
  oc.add_equation( BC,  {t,x,y}, {T_INT, FFDom::LB, Y_INT}, bnd_opt );
  oc.add_equation( BC,  {t,x,y}, {T_INT, FFDom::UB, Y_INT}, bnd_opt );
  oc.add_equation( BC,  {t,x,y}, {T_INT, FFDom::ALL, FFDom::LB}, bnd_opt );
  oc.add_equation( BC,  {t,x,y}, {T_INT, FFDom::ALL, FFDom::UB}, bnd_opt );
  oc.set_evolution_domain( t );

  oc.add_output( U, {t,x,y}, std::vector<double>{ p.tf, p.x0, p.y0 } );
  FFVar Iu = OpI( OpI( U, x ), y );
  oc.add_output( Iu, {t}, std::vector<double>{ p.tf } );

  if( !oc.setup() ){ std::cerr << "setup FAILED\n"; return false; }

  oc.register_control( A );
  oc.register_control( W0 );
  size_t const nf  = oc.n_colloc_fct();
  size_t const ncd = oc.n_control_dof();
  std::cout << "outputs nf=" << nf << "  controls ncd=" << ncd
            << " (canonical order: A, W0)\n";
  if( nf != 2 || ncd != 2 ){ std::cerr << "unexpected nf/ncd\n"; return false; }

  std::vector<double> var, inp;
  if( !oc.init( var, inp, nullptr ) ){ std::cerr << "init FAILED\n"; return false; }
  oc.set_input_values( A,  std::vector<double>{ p.A0   }, inp.data() );
  oc.set_input_values( W0, std::vector<double>{ p.W0_0 }, inp.data() );

  // control order in the reduced vector: canonical FFVar order.  Identify each control's offset.
  auto const& C = oc.controls();
  size_t const offA = C.at(A).offset, offW = C.at(W0).offset;

  bool ok = true;
  Fana const ref = analytic( p.A0, p.W0_0, p );

  // ---- solve ----
  {
    std::vector<double> xv( var );
    OCFESLV::SolveReport const srep = oc.solve( xv.data(), inp.data(), nullptr );
    if( !srep.converged )
      std::cerr << "OCFESLV::solve did not converge: final|r|=" << srep.final_residual
                << " after " << srep.iterations << " it\n";
  }

  // ---- forward sensitivity ----
  std::vector<double> Jf, Ff;
  {
    std::vector<double> xv( var );
    if( !oc.solve_fsens( xv.data(), inp.data(), nullptr ) ){ std::cerr << "solve_fsens FAILED\n"; return false; }
    Ff = oc.sens_functions();
    Jf = oc.sens_jacobian();          // row-major nf x ncd
  }
  // ---- adjoint sensitivity ----
  std::vector<double> Ja, Fa;
  {
    std::vector<double> xv( var );
    if( !oc.solve_asens( xv.data(), inp.data(), nullptr ) ){ std::cerr << "solve_asens FAILED\n"; return false; }
    Fa = oc.sens_functions();
    Ja = oc.sens_jacobian();
  }
  if( Jf.size() != nf*ncd || Ja.size() != nf*ncd ){ std::cerr << "jacobian size mismatch\n"; return false; }

  auto Jf_ = [&]( size_t f, size_t off ){ return Jf[ f*ncd + off ]; };
  auto Ja_ = [&]( size_t f, size_t off ){ return Ja[ f*ncd + off ]; };

  // (1) fsens vs asens : identical functions and reduced Jacobian
  std::cout << "\n[1] forward vs adjoint (machinery cross-check)\n";
  check_close( "f0 value  (fsens vs asens)", Ff[0], Fa[0], 1e-9, ok );
  check_close( "f1 value  (fsens vs asens)", Ff[1], Fa[1], 1e-9, ok );
  check_close( "df0/dA    (fsens vs asens)", Jf_(0,offA), Ja_(0,offA), 1e-8, ok );
  check_close( "df0/dW0   (fsens vs asens)", Jf_(0,offW), Ja_(0,offW), 1e-8, ok );
  check_close( "df1/dA    (fsens vs asens)", Jf_(1,offA), Ja_(1,offA), 1e-8, ok );
  check_close( "df1/dW0   (fsens vs asens)", Jf_(1,offW), Ja_(1,offW), 1e-8, ok );

  // (2) output values vs analytic
  std::cout << "\n[2] output values vs analytic (bump peak; unit-square Gaussian inventory)\n";
  check_close( "f0 = A*W0/S(tf)",           Ff[0], ref.f0, 5e-3, ok );
  check_close( "f1 = pi*A*W0*erf^2",        Ff[1], ref.f1, 5e-3, ok );

  // (3) reduced Jacobian vs analytic
  std::cout << "\n[3] reduced Jacobian vs analytic dF/d(control)\n";
  check_close( "df0/dA",  Jf_(0,offA), ref.df0dA,  5e-3, ok );
  check_close( "df0/dW0", Jf_(0,offW), ref.df0dW0, 5e-3, ok );
  check_close( "df1/dA",  Jf_(1,offA), ref.df1dA,  5e-3, ok );
  check_close( "df1/dW0", Jf_(1,offW), ref.df1dW0, 5e-3, ok );

  // (3b) independent central-finite-difference reference (re-solves the primal at A+-h, W0+-h)
  std::cout << "\n[3b] reduced Jacobian vs central finite differences\n";
  double const hA = 1e-6*std::max(1.0,std::abs(p.A0));
  double const hW = 1e-6*std::max(1.0,std::abs(p.W0_0));
  std::array<double,2> fAp, fAm, fWp, fWm;
  bool fd_ok = solve_values(march, p.A0+hA, p.W0_0, p, fAp)
            && solve_values(march, p.A0-hA, p.W0_0, p, fAm)
            && solve_values(march, p.A0, p.W0_0+hW, p, fWp)
            && solve_values(march, p.A0, p.W0_0-hW, p, fWm);
  if( fd_ok ){
    double const fd_f0_A = (fAp[0]-fAm[0])/(2*hA), fd_f0_W = (fWp[0]-fWm[0])/(2*hW);
    double const fd_f1_A = (fAp[1]-fAm[1])/(2*hA), fd_f1_W = (fWp[1]-fWm[1])/(2*hW);
    check_close( "df0/dA  (vs FD)",  Jf_(0,offA), fd_f0_A, 1e-5, ok );
    check_close( "df0/dW0 (vs FD)",  Jf_(0,offW), fd_f0_W, 1e-5, ok );
    check_close( "df1/dA  (vs FD)",  Jf_(1,offA), fd_f1_A, 1e-5, ok );
    check_close( "df1/dW0 (vs FD)",  Jf_(1,offW), fd_f1_W, 1e-5, ok );
  }
  else { std::cerr << "  finite-difference reference solve FAILED\n"; ok = false; }

  std::cout << "\nPDE5b " << tag << ": " << (ok?"PASS":"FAIL") << "\n";
  return ok;
}

} // namespace

int main()
{
  Par p;
  std::cout << "PDE5b: output functionals (peak, xy-inventory) + fwd/adj sensitivity\n"
            << "controls: A (magnitude), W0 (spread); nominal A=" << p.A0 << " W0=" << p.W0_0 << "\n";

  bool const okm = run_mode(false, p);   // monolithic (works with current header)
  bool const okc = run_mode(true,  p);   // marching   (requires override-path sensitivity reinit)

  std::cout << "\n==================== PDE5b summary ====================\n"
            << "  monolithic : " << (okm?"PASS":"FAIL") << "\n"
            << "  marching   : " << (okc?"PASS":"FAIL")
            << (okc?"":"   (expected until override-path sens transfer lands)") << "\n"
            << "======================================================\n";
  return (okm && okc) ? 0 : 1;
}
