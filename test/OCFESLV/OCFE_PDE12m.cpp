// ===========================================================================
// OCFE_PDE12m_solve2.cpp  --  MANUAL reduced-first-order verification, ORDER 4
//
// Hand-built reduced form of a 4th-order (dissipative) evolution test
//
//   d_t u + d_xxxx u = f(t,x)          (stable: mode e^{ikx} decays as e^{-k^4 t})
//
// reduced by hand (RED_NONE) into the explicit depth-3 first-order chain
//
//   states:  u , D1 , D2 , D3
//   PDE  :   d_t u + d_x D3 = f         (D3 = u_xxx, so d_x D3 = u_xxxx)
//   LINK1:   D1 - d_x u  = 0            (D1 = u_x)
//   LINK2:   D2 - d_x D1 = 0            (D2 = u_xx)
//   LINK3:   D3 - d_x D2 = 0            (D3 = u_xxx)
//
// PURPOSE.  Oracle for the *multi-displacement* boundary closure that a 4th-order
// operator needs (deficit = order-2 = 2).  An asymmetric, well-posed BC split is
// used so that ONE face is over-determined by TWO -- the hardest case, which the
// 3rd-order test (deficit 1) never reaches:
//
//   x=LB :  value(u) + deriv1(D1) + deriv2(D2)        [3 user conditions]
//   x=UB :  value(u)                                   [1 user condition]
//
// Closure under test: at x=LB the two derivative BCs pin D1 and D2, so LINK1 AND
// LINK2 are displaced there (LINK3 kept); at x=UB nothing is displaced.  Per-node
// row count is then exactly 4 = #states everywhere:
//   (t=LB , any x ):  IC(u)              + LINK1 + LINK2 + LINK3
//   (t>LB , x int ):  PDE                + LINK1 + LINK2 + LINK3
//   (t>LB , x=LB  ):  value+deriv1+deriv2                + LINK3   [LINK1,LINK2 dropped]
//   (t>LB , x=UB  ):  value              + LINK1 + LINK2 + LINK3
//
// Dropping the WRONG LINKs (e.g. LINK2,LINK3 instead of LINK1,LINK2) leaves an
// aux unpinned -> singular; this driver hard-wires the correct drops so the run
// confirms which drops a rank-based auto resolver (avenue B) must reproduce.
//
// SINGLE element in each direction isolates the boundary closure from interfaces.
// ===========================================================================

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include <armadillo>

#ifndef OCFE_OCFESLV_HEADER
#define OCFE_OCFESLV_HEADER "ocfeslv.hpp"
#endif
#include OCFE_OCFESLV_HEADER

using namespace mc;

#ifndef TEST_KDV_NEL_T
#define TEST_KDV_NEL_T 1
#endif
#ifndef TEST_KDV_NEL_X
#define TEST_KDV_NEL_X 1
#endif
#ifndef TEST_KDV_NT
#define TEST_KDV_NT 6
#endif
#ifndef TEST_KDV_NX
#define TEST_KDV_NX 20
#endif
#ifndef TEST_KDV_TF
#define TEST_KDV_TF 0.5
#endif
#ifndef TEST_KDV_XF
#define TEST_KDV_XF 1.0
#endif
#ifndef TEST_KDV_MAXIT
#define TEST_KDV_MAXIT 40
#endif
#ifndef TEST_KDV_SOLVE_TOL
#define TEST_KDV_SOLVE_TOL 1e-9
#endif
// Order-aware accuracy gates: each spectral differentiation amplifies truncation
// error by ~N^2, so D_j is inherently coarser than u; gates tighten as NX grows.
#ifndef TEST_KDV_U_TOL
#define TEST_KDV_U_TOL  1e-6
#endif
#ifndef TEST_KDV_D1_TOL
#define TEST_KDV_D1_TOL 1e-6
#endif
#ifndef TEST_KDV_D2_TOL
#define TEST_KDV_D2_TOL 1e-4
#endif
#ifndef TEST_KDV_D3_TOL
#define TEST_KDV_D3_TOL 1e-3
#endif

static double const kPi = 3.14159265358979323846;

struct Par {
  double tf  = TEST_KDV_TF;
  double xf  = TEST_KDV_XF;
  double k   = 2.0*kPi;
  double phi = 0.7;
  double a   = 0.3;
  double u0  = 1.0;
};

// u = u0 sin(kx+phi)(1+at) and its x-derivatives (d/dx: sin->cos->-sin->-cos->sin).
static double U_exact  ( double t, double x, Par const& p ){ return  p.u0*std::sin(p.k*x+p.phi)*(1.0+p.a*t); }
static double Ux_exact ( double t, double x, Par const& p ){ return  p.u0*std::pow(p.k,1)*std::cos(p.k*x+p.phi)*(1.0+p.a*t); }
static double Uxx_exact( double t, double x, Par const& p ){ return -p.u0*std::pow(p.k,2)*std::sin(p.k*x+p.phi)*(1.0+p.a*t); }
static double Uxxx_exact(double t, double x, Par const& p ){ return -p.u0*std::pow(p.k,3)*std::cos(p.k*x+p.phi)*(1.0+p.a*t); }

static char const* imp_name( OCFESLV::Options::ImpositionType t ){
  switch( t ){
    case OCFESLV::Options::IC_WEAK:   return "IC_WEAK";
    case OCFESLV::Options::IC_TRACE:  return "IC_TRACE";
    case OCFESLV::Options::IC_STRONG: return "IC_STRONG";
    default:                        return "IC_?";
  }
}
static double max_abs( std::vector<double> const& v ){
  double m=0.0; for( double x: v ) m = std::max(m,std::fabs(x)); return m;
}
static bool check_close( std::string const& label, double value, double tol ){
  bool ok = ( value <= tol );
  std::cout << std::left << std::setw(34) << label
            << " value=" << std::scientific << std::setprecision(6) << value
            << " tol=" << tol << "  " << (ok?"PASS":"FAIL") << "\n";
  return ok;
}

struct ModeResult { std::string name; bool ok=false; size_t nVar=0,nEqn=0,nTrace=0;
                    bool solved=false; double final_res=0., eU=0., eD1=0., eD2=0., eD3=0.; };

static ModeResult run_mode( OCFESLV::Options::ImpositionType imp, Par const& p )
{
  ModeResult R; R.name = imp_name(imp); bool ok=true;

  std::cout << "\n===== MANUAL reduced 4th-order test ("<< R.name <<") =====\n";
  std::cout << "elements: t=" << TEST_KDV_NEL_T << " x=" << TEST_KDV_NEL_X
            << "  nodes/elem: t=" << TEST_KDV_NT << " x=" << TEST_KDV_NX
            << "  k=" << p.k << "\n";

  FFGraph DAG;
  FFVar t  = DAG.add_var("t");
  FFVar x  = DAG.add_var("x");
  FFVar u  = DAG.add_var("u(t,x)");
  FFVar D1 = DAG.add_var("D1(t,x)");
  FFVar D2 = DAG.add_var("D2(t,x)");
  FFVar D3 = DAG.add_var("D3(t,x)");
  FFPartial OpP;

  FFVar sinx = sin( p.k*x + p.phi );
  FFVar cosx = cos( p.k*x + p.phi );
  double const k4 = p.k*p.k*p.k*p.k;
  FFVar UE   = p.u0*sinx*( 1.0 + p.a*t );
  FFVar UEx  = p.u0*p.k*cosx*( 1.0 + p.a*t );
  FFVar UExx = -p.u0*(p.k*p.k)*sinx*( 1.0 + p.a*t );
  // f = u_t + u_xxxx = u0 a sin + u0 k^4 sin (1+at).
  FFVar FE   = p.u0*p.a*sinx + p.u0*k4*sinx*( 1.0 + p.a*t );

  // Hand-reduced depth-3 first-order chain (RED_NONE; only first-order OpP).
  FFVar PDE   = OpP(u,t) + OpP(D3,x) - FE;   // d_t u + d_x D3 = f
  FFVar LINK1 = D1 - OpP(u,x);               // D1 = u_x
  FFVar LINK2 = D2 - OpP(D1,x);              // D2 = D1_x
  FFVar LINK3 = D3 - OpP(D2,x);              // D3 = D2_x
  FFVar BCVAL = u  - UE;                      // u  = UE
  FFVar BCDX  = D1 - UEx;                     // D1 = UE_x
  FFVar BCDXX = D2 - UExx;                    // D2 = UE_xx

  OCFESLV oc(&DAG);
  oc.add_domain( t, FFDom(0., p.tf, TEST_KDV_NEL_T, FFDom::LGR, TEST_KDV_NT) );
  oc.add_domain( x, FFDom(0., p.xf, TEST_KDV_NEL_X, FFDom::LGL, TEST_KDV_NX) );
  oc.add_state( u,  {t,x} );
  oc.add_state( D1, {t,x} );
  oc.add_state( D2, {t,x} );
  oc.add_state( D3, {t,x} );
  oc.update_ref( u,  [&]( OCFESLV::t_Coord const& c ){ return U_exact  (c.at(t),c.at(x),p); } );
  oc.update_ref( D1, [&]( OCFESLV::t_Coord const& c ){ return Ux_exact (c.at(t),c.at(x),p); } );
  oc.update_ref( D2, [&]( OCFESLV::t_Coord const& c ){ return Uxx_exact(c.at(t),c.at(x),p); } );
  oc.update_ref( D3, [&]( OCFESLV::t_Coord const& c ){ return Uxxx_exact(c.at(t),c.at(x),p); } );

  OCFESLV::EqnOptions int_opt( OCFESLV::EqnRole::INTERIOR, 0 );
  OCFESLV::EqnOptions ini_opt( OCFESLV::EqnRole::INITIAL,  0 );
  OCFESLV::EqnOptions bnd_opt( OCFESLV::EqnRole::BOUNDARY, 0 );

  int const T_INT   = FFDom::ALL - FFDom::LB;
  int const X_INT   = FFDom::ALL - FFDom::LB - FFDom::UB;
  int const X_NO_LB = FFDom::ALL - FFDom::LB;

  // PDE on x-interior, t>0.
  oc.add_equation( PDE,   {t,x}, {T_INT, X_INT},          int_opt );
  // Initial data, whole spatial grid.
  oc.add_equation( BCVAL, {t,x}, {FFDom::LB, FFDom::ALL},  ini_opt );
  // LINK1, LINK2 : everywhere EXCEPT (t>LB, x=LB) -- displaced there by deriv1, deriv2.
  oc.add_equation( LINK1, {t,x}, {T_INT, X_NO_LB},         int_opt );
  oc.add_equation( LINK1, {t,x}, {FFDom::LB, FFDom::ALL},  int_opt );
  oc.add_equation( LINK2, {t,x}, {T_INT, X_NO_LB},         int_opt );
  oc.add_equation( LINK2, {t,x}, {FFDom::LB, FFDom::ALL},  int_opt );
  // LINK3 : everywhere.
  oc.add_equation( LINK3, {t,x}, {FFDom::ALL, FFDom::ALL}, int_opt );
  // Asymmetric boundary conditions (t>0): value+deriv1+deriv2 at LB, value at UB.
  oc.add_equation( BCVAL, {t,x}, {T_INT, FFDom::LB},       bnd_opt );  // u  = UE    at x=LB
  oc.add_equation( BCDX,  {t,x}, {T_INT, FFDom::LB},       bnd_opt );  // D1 = UE_x  at x=LB  (takes LINK1's slot)
  oc.add_equation( BCDXX, {t,x}, {T_INT, FFDom::LB},       bnd_opt );  // D2 = UE_xx at x=LB  (takes LINK2's slot)
  oc.add_equation( BCVAL, {t,x}, {T_INT, FFDom::UB},       bnd_opt );  // u  = UE    at x=UB

  oc.set_evolution_domain( t );
  oc.options.REDUCE.ORDER    = OCFESLV::Options::RED_NONE;
  oc.options.CLASSIFY.MODE        = OCFESLV::Options::CLASS_AUTO;
  oc.options.INTERFACE.TYPE  = OCFESLV::Options::IC_AUTO;
  oc.options.INTERFACE.IMPOSITION = imp;
  oc.options.INTERFACE.SAT_SIGMA0      = 10.0;
  oc.options.DISPLAY_LEVEL   = 1;

  if( !oc.setup() ){
    std::cerr << "ERROR: OCFESLV::setup() failed for " << R.name << "\n";
    R.ok=false; return R;
  }

  size_t const nVar=oc.n_colloc_sta();
  size_t const nEqn=oc.n_colloc_eqn();
  size_t const nTrace=oc.n_colloc_trace();
  R.nVar=nVar; R.nEqn=nEqn; R.nTrace=nTrace;
  std::cout << "nVar=" << nVar << " nEqn=" << nEqn << " nTrace=" << nTrace
            << " square=" << (nVar==nEqn?"yes":"no") << "\n";
  ok &= check_close("square (nVar==nEqn)", nVar==nEqn?0.0:1.0, 0.0);

  std::cout << "States:";
  for( auto const& st: oc.states_colloc() ) std::cout << ' ' << st.name();
  std::cout << "\nPDE type: " << OCFESLV::pde_type_name(oc.pde_type().type) << "\n";

  // Initial guess: exact primitive, auxes at zero -> the chain + BCs must recover
  // u_x, u_xx, u_xxx if the multi-displacement closure is non-singular.
  std::vector<double> var(nVar,0.0);
  {
    size_t off=0, sidx=0;
    for( auto const& st: oc.states_colloc() ){
      auto nodes = oc.node_colloc(st);
      for( size_t i=0; i<nodes.size(); ++i )
        var[off+i] = (sidx==0) ? U_exact( nodes[i][0], nodes[i][1], p ) : 0.0;
      off += nodes.size(); ++sidx;
    }
  }

  oc.options.SOLVE.MAX_ITER = TEST_KDV_MAXIT;
  oc.options.SOLVE.RES_TOL  = TEST_KDV_SOLVE_TOL;
#if defined(CRONOS__WITH_SPQR)
  oc.options.SOLVE.FACTORIZATION = OCFESLV::Options::SOLVE_SPQR;
#endif
  OCFESLV::SolveReport const srep = oc.solve( var.data() );
  R.solved = srep.converged;
  if( !R.solved )
    std::cerr << "OCFESLV::solve did not converge: final|r|=" << srep.final_residual
              << " after " << srep.iterations << " it\n";
  ok &= R.solved;

  std::vector<double> res(nEqn,0.0);
  oc.eval(res.data(),nullptr,var.data(),nullptr,nullptr);
  R.final_res = max_abs(res);
  ok &= check_close("final max residual", R.final_res, TEST_KDV_SOLVE_TOL);

  {
    double eU=0.,eD1=0.,eD2=0.,eD3=0.; size_t off=0, sidx=0;
    for( auto const& st: oc.states_colloc() ){
      auto nodes = oc.node_colloc(st);
      for( size_t i=0; i<nodes.size(); ++i ){
        double tt=nodes[i][0], xx=nodes[i][1], val=var[off+i];
        if     ( sidx==0 ) eU  = std::max( eU,  std::fabs(val - U_exact  (tt,xx,p)) );
        else if( sidx==1 ) eD1 = std::max( eD1, std::fabs(val - Ux_exact (tt,xx,p)) );
        else if( sidx==2 ) eD2 = std::max( eD2, std::fabs(val - Uxx_exact(tt,xx,p)) );
        else if( sidx==3 ) eD3 = std::max( eD3, std::fabs(val - Uxxx_exact(tt,xx,p)) );
      }
      off += nodes.size(); ++sidx;
    }
    R.eU=eU; R.eD1=eD1; R.eD2=eD2; R.eD3=eD3;
    ok &= check_close("max |u  - u_exact|",  eU,  TEST_KDV_U_TOL);
    ok &= check_close("max |D1 - u_x|",      eD1, TEST_KDV_D1_TOL);
    ok &= check_close("max |D2 - u_xx|",     eD2, TEST_KDV_D2_TOL);
    ok &= check_close("max |D3 - u_xxx|",    eD3, TEST_KDV_D3_TOL);
  }

  std::cout << "MANUAL reduced 4th-order (" << R.name << "): " << (ok?"PASS":"FAIL") << "\n";
  R.ok = ok;
  return R;
}

int main()
{
  Par p;
  std::vector<ModeResult> results;
  results.push_back( run_mode(OCFESLV::Options::IC_WEAK,   p) );
  results.push_back( run_mode(OCFESLV::Options::IC_TRACE,  p) );
  results.push_back( run_mode(OCFESLV::Options::IC_STRONG, p) );

  std::cout << "\n========== MANUAL reduced 4th-order summary ==========\n";
  std::cout << std::left << std::setw(11) << "mode" << std::setw(7) << "square"
            << std::setw(13) << "final|r|" << std::setw(13) << "|u-ex|"
            << std::setw(13) << "|D1-ex|" << std::setw(13) << "|D2-ex|"
            << std::setw(13) << "|D3-ex|" << "result\n";
  bool all_ok=true;
  for( auto const& r: results ){
    std::cout << std::left << std::setw(11) << r.name
              << std::setw(7) << (r.nVar==r.nEqn?"yes":"no")
              << std::scientific << std::setprecision(3)
              << std::setw(13) << r.final_res << std::setw(13) << r.eU
              << std::setw(13) << r.eD1 << std::setw(13) << r.eD2
              << std::setw(13) << r.eD3 << (r.ok?"PASS":"FAIL") << "\n";
    all_ok &= r.ok;
  }
  std::cout << "======================================================\n";
  std::cout << "MANUAL reduced 4th-order (double-displacement) closure: " << (all_ok?"PASS":"FAIL") << "\n";
  return all_ok ? 0 : 1;
}
