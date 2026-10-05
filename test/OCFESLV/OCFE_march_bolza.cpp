// ============================================================================
//  OCFE_march_bolza.cpp  --  sensitivities of CAPTURED (nonlinear) evolution-
//  integral outputs (ACCUM / OpI), harmonised with the linear (fctsum) path.
//
//  Control a (LUMPED input).  dx/dt = a,  x(0)=0   =>   x(t) = a t.
//    I1 = INT_0^1 x  dt   = a/2 ,   I2 = INT_0^1 x^2 dt = a^2/3
//  Outputs:  F0 = I1 (linear) ;  F1 = I1^2 (captured) ;  F2 = I1*I2 (captured cross)
//  Analytic:  dF0/da = 1/2 ;  dF1/da = a/2 ;  dF2/da = a^2/2.
//  At a=1.2:  dF0=0.5, dF1=0.6, dF2=0.72.
//
//  GRIDS: run on BOTH uniform and NON-UNIFORM evolution grids.  The integrals are
//  exact for the (degree<=2) MMS on any grid, so analytic answers are grid-independent;
//  a graded grid exercises the ACCUM quadrature transpose (jac*w_ref) on VARIABLE jac.
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
  std::cout << "    " << std::left << std::setw(30) << nm << std::right
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
  FFVar I1 = OpI( x, t ), I2 = OpI( x*x, t );
  oc.add_output( I1      );                                // F0 linear
  oc.add_output( I1 * I1 );                                // F1 captured (square)
  oc.add_output( I1 * I2 );                                // F2 captured (product cross term)
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
  oc.solve( xv.data(), ii.data(), nullptr );
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
  std::vector<double> xa = vi;
  if( !oc.solve_asens( xa.data(), ii.data(), nullptr ) ){ std::cout << "    solve_asens FAILED\n"; ++g_fail; return; }
  std::vector<double> const Ja = oc.sens_jacobian();
  if( Jf.size() < 3 || Ja.size() < 3 ){ std::cout << "    jacobian too small\n"; ++g_fail; return; }

  double const h = 1e-6;
  std::vector<double> Fp = value_at( march, aval + h, nu ), Fm = value_at( march, aval - h, nu );
  std::vector<double> gFD( 3, 0. );
  for( size_t r = 0; r < 3 && r < Fp.size() && r < Fm.size(); ++r ) gFD[r] = ( Fp[r] - Fm[r] ) / ( 2*h );

  double const dF0 = 0.5, dF1 = aval/2., dF2 = aval*aval/2.;
  ck( "fwd dF0/da (linear)",   Jf[0], dF0, 1e-6 );
  ck( "fwd dF1/da (captured)", Jf[1], dF1, 1e-6 );
  ck( "fwd dF2/da (captured)", Jf[2], dF2, 1e-6 );
  ck( "adj dF0/da == fwd",     Ja[0], Jf[0], 1e-9 );
  ck( "adj dF1/da == fwd",     Ja[1], Jf[1], 1e-9 );
  ck( "adj dF2/da == fwd",     Ja[2], Jf[2], 1e-9 );
  ck( "fwd dF1/da == FD",      Jf[1], gFD[1], 1e-4 );
  ck( "fwd dF2/da == FD",      Jf[2], gFD[2], 1e-4 );
}

int main()
{
  std::cout << "===============================================================\n"
            << "  captured f(INT R) output sensitivities (ACCUM)\n"
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
