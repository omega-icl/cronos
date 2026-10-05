// ============================================================================
//  OCFE_capture_mixed.cpp  --  MIXED captured outputs f(cap, state) with a nonzero
//  direct state dependence (df/dvar != 0), probing the Stage-1b adjoint over-seed
//  caveat.  Control a: dx/dt = a, x(0)=0  =>  x(t) = a t.
//    I1 = INT_0^1 x dt = a/2 ,   x(1) = a  (terminal state)
//  Outputs (presented at t=1 so the direct-state term is the terminal value):
//    F0 = I1        (pure cap, linear)      dF0/da = 1/2
//    F1 = I1 + x    (cap + terminal state)  dF1/da = 3/2
//    F2 = I1 * x    (cap * terminal state)  dF2/da = a
//  At a=1.2:  F0=0.6, F1=1.8, F2=0.72 ; dF0=0.5, dF1=1.5, dF2=1.2.
//
//  GRIDS: run on BOTH uniform and NON-UNIFORM evolution grids.  The graded grid puts
//  the terminal state x(1) on an irregular LAST window, stressing the capW window-gating
//  of the direct-state adjoint (the over-seed guard) with a non-constant jac.
//
//  Build: -DOCFE_OCFESLV_HEADER='"ocfeslv.hpp"'
// ============================================================================
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

#ifndef OCFE_OCFESLV_HEADER
#define OCFE_OCFESLV_HEADER "ocfeslv.hpp"
#endif
#include OCFE_OCFESLV_HEADER

using namespace mc;

static int g_fail = 0;
static void ck( char const* nm, double got, double want, double tol )
{
  bool const p = std::fabs( got - want ) < tol;
  if( !p ) ++g_fail;
  std::cout << "    " << std::left << std::setw(28) << nm << std::right
            << std::scientific << std::setprecision(6)
            << " got=" << std::setw(14) << got << " want=" << std::setw(14) << want
            << "  " << ( p ? "PASS" : "FAIL" ) << "\n";
}

static FFDom t_dom( bool nu ){ return nu ? FFDom( std::vector<double>{ 0.0, 0.3, 0.6, 1.0 }, FFDom::LGR, 4 )
                                         : FFDom( 0., 1., 3, FFDom::LGR, 4 ); }

static void build( FFGraph& DAG, OCFESLV& oc, FFVar& t, FFVar& x, FFVar& a, bool march, bool nu )
{
  t = DAG.add_var("t"); x = DAG.add_var("x(t)"); a = DAG.add_var("a");
  FFPartial OpP; FFIntegral OpI;
  oc.set_evolution_domain( t );
  oc.add_domain( t, t_dom( nu ) );
  oc.add_state ( x, { t } );
  oc.add_input ( a, {} );
  oc.add_equation( OpP(x,t) - a, { t }, { FFDom::ALL - FFDom::LB }, OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( x, { t }, { FFDom::LB }, OCFESLV::EqnOptions( OCFESLV::EqnRole::INITIAL, 0 ) );
  FFVar I1 = OpI( x, t );
  oc.add_output( I1,       { t }, { 1.0 } );   // F0 pure cap (presented at t=1; I1 is t-free)
  oc.add_output( I1 + x,   { t }, { 1.0 } );   // F1 cap + terminal state x(1)
  oc.add_output( I1 * x,   { t }, { 1.0 } );   // F2 cap * terminal state x(1)
  oc.options.REDUCE.ORDER   = OCFESLV::Options::RED_FULL;
  oc.options.SOLVE.MARCHING = march;
  oc.options.DISPLAY_LEVEL  = 0;
}

static std::vector<double> value_at( bool march, double aval, bool nu )
{
  FFGraph DAG; OCFESLV oc( &DAG ); FFVar t,x,a; build( DAG, oc, t, x, a, march, nu );
  if( !oc.setup() ) return {};
  std::vector<double> vi, ii; oc.init( vi, ii, nullptr );
  oc.set_input_values( a, { aval }, ii.data() );
  std::vector<double> xv = vi;
  if( !oc.solve( xv.data(), ii.data(), nullptr ).converged ) return {};
  return oc.val_functions();
}

static void run( bool march, double aval, bool nu )
{
  std::cout << "\n[" << ( march ? "MARCH" : "MONO" ) << "]  a=" << aval << "\n";
  FFGraph DAG; OCFESLV oc( &DAG ); FFVar t,x,a; build( DAG, oc, t, x, a, march, nu );
  if( !oc.setup() ){ std::cout << "    setup FAILED\n"; ++g_fail; return; }
  std::vector<double> vi, ii; oc.init( vi, ii, nullptr );
  oc.set_input_values( a, { aval }, ii.data() );
  oc.register_control( a );

  std::vector<double> xv = vi;
  if( !oc.solve_fsens( xv.data(), ii.data(), nullptr ) ){ std::cout << "    solve_fsens FAILED\n"; ++g_fail; return; }
  std::vector<double> const Jf = oc.sens_jacobian();
  std::vector<double> const Fv = oc.sens_functions();
  std::vector<double> xa = vi;
  if( !oc.solve_asens( xa.data(), ii.data(), nullptr ) ){ std::cout << "    solve_asens FAILED\n"; ++g_fail; return; }
  std::vector<double> const Ja = oc.sens_jacobian();
  if( Jf.size() < 3 || Ja.size() < 3 || Fv.size() < 3 ){ std::cout << "    jacobian/functions too small\n"; ++g_fail; return; }

  double const h = 1e-6;
  std::vector<double> Fp = value_at( march, aval + h, nu ), Fm = value_at( march, aval - h, nu );
  std::vector<double> gFD( 3, 0. );
  for( size_t r = 0; r < 3 && r < Fp.size() && r < Fm.size(); ++r ) gFD[r] = ( Fp[r]-Fm[r] )/( 2*h );

  double const I1 = aval/2., xT = aval;
  double const vF0 = I1, vF1 = I1 + xT, vF2 = I1 * xT;
  double const dF0 = 0.5, dF1 = 1.5, dF2 = xT/2. + I1;   // = a
  ck( "F0 value (cap)",        Fv[0], vF0, 1e-6 );
  ck( "F1 value (cap+state)",  Fv[1], vF1, 1e-6 );
  ck( "F2 value (cap*state)",  Fv[2], vF2, 1e-6 );
  ck( "fwd dF0/da",            Jf[0], dF0, 1e-5 );
  ck( "fwd dF1/da (mixed)",    Jf[1], dF1, 1e-5 );
  ck( "fwd dF2/da (mixed)",    Jf[2], dF2, 1e-5 );
  ck( "adj dF1 == fwd (mixed)", Ja[1], Jf[1], 1e-8 );
  ck( "adj dF2 == fwd (mixed)", Ja[2], Jf[2], 1e-8 );
  ck( "fwd dF1 == FD",         Jf[1], gFD[1], 1e-4 );
  ck( "fwd dF2 == FD",         Jf[2], gFD[2], 1e-4 );
}

int main()
{
  std::cout << "===============================================================\n"
            << "  MIXED captured outputs f(cap, state): adjoint over-seed probe\n"
            << "  (validated on UNIFORM and NON-UNIFORM evolution grids)\n"
            << "===============================================================\n";
  for( bool nu : { false, true } ){
    std::cout << "\n========== " << ( nu ? "NON-UNIFORM" : "UNIFORM" ) << " grid ==========\n";
    run( true,  1.2, nu );
    run( false, 1.2, nu );
  }
  std::cout << "\n  RESULT: " << ( g_fail==0 ? "ALL PASS" : "FAIL" ) << "  (" << g_fail << " failing)\n";
  return g_fail==0 ? 0 : 1;
}
