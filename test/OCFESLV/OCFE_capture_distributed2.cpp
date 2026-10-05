// ============================================================================
//  OCFE_capture_distributed2.cpp  --  coverage gate for DISTRIBUTED LATCH
//  captured reductions: cap = R(tau) (OpEval point-eval at a fixed evolution
//  coordinate tau), distributed over the spatial dimension z.  Mirrors
//  OCFE_capture_distributed.cpp (ACCUM / OpI flavor) so BOTH capture flavors are
//  exercised over a spatial dimension, value AND sensitivity, mono AND marched,
//  linear AND nonlinear.
//
//  Manufactured {t,z} problem, control a (LUMPED input):
//    d_t c - D d_zz c = a z ,  c(0,z)=z ,  c(t,0)=0 ,  d_z c|_{z=1}=1+a t
//    => c(t,z) = z (1 + a t)
//  LATCH at tau=0.25:  L(z) = c(tau,z) = z (1 + a tau)      (distributed over z)
//    Llin = L|_{z=1}   = 1 + a tau          (LINEAR   LATCH)
//    Lsq  = L^2|_{z=1} = (1 + a tau)^2      (NONLINEAR LATCH)
//  With tau=0.25, a=1.2:  Llin=1.3, Lsq=1.69, dLlin/da=tau=0.25, dLsq/da=0.65.
//
//  GRIDS: run on BOTH uniform and NON-UNIFORM grids, in t (evolution) AND z (spatial).
//  The MMS is degree-1 in t and z, so the analytic answers are grid-independent, but a
//  graded grid makes jac vary per element, lands tau at an arbitrary interior xi of an
//  irregular EVOLUTION element (stressing the latch element-find + L_m(tau) transpose),
//  and grades the SPATIAL dofs (stressing the latch scatter's spatial indexing) --
//  exactly the structure a uniform grid masked in the latch-adjoint bug.
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
  std::cout << "    " << std::left << std::setw(34) << nm << std::right
            << std::scientific << std::setprecision(6)
            << " got=" << std::setw(14) << got << " want=" << std::setw(14) << want
            << "  " << ( p ? "PASS" : "FAIL" ) << "\n";
}

static double const Dz  = 0.1;
static double const TAU = 0.25;

// Uniform vs NON-UNIFORM element boundaries.  tau=0.25 lands in the FIRST t-element for both
// grids (interior, non-terminal), at xi = +0.667 on the non-uniform [0,0.3,...] element.
static FFDom t_dom( bool nu ){ return nu ? FFDom( std::vector<double>{ 0.0, 0.3, 0.6, 1.0 }, FFDom::LGR, 4 )
                                         : FFDom( 0., 1., 3, FFDom::LGR, 4 ); }
static FFDom z_dom( bool nu ){ return nu ? FFDom( std::vector<double>{ 0.0, 0.3, 1.0 },      FFDom::LGL, 4 )
                                         : FFDom( 0., 1., 2, FFDom::LGL, 4 ); }

static void build( FFGraph& DAG, OCFESLV& oc, FFVar& t, FFVar& z, FFVar& c, FFVar& a,
                   bool march, bool nu )
{
  t = DAG.add_var("t"); z = DAG.add_var("z"); c = DAG.add_var("c(t,z)"); a = DAG.add_var("a");
  FFPartial OpP; FFEval OpEval;

  FFVar PDE    = OpP(c,t) - Dz*OpP(c,{z,2}) - a*z;      // d_t c - D d_zz c = a z
  FFVar BCINIT = c - z;                                 // c(0,z) = z
  FFVar BCLO   = c;                                     // c(t,0) = 0
  FFVar BCHI   = OpP(c,z) - ( 1.0 + a*t );              // d_z c|_{z=1} = 1 + a t

  oc.add_domain( t, t_dom( nu ) );
  oc.add_domain( z, z_dom( nu ) );
  oc.add_state ( c, { t, z } );
  oc.add_input ( a, {} );

  OCFESLV::EqnOptions io( OCFESLV::EqnRole::INTERIOR, 0 ), ii( OCFESLV::EqnRole::INITIAL, 0 ), ib( OCFESLV::EqnRole::BOUNDARY, 0 );
  int const T_INT = FFDom::ALL - FFDom::LB;
  int const Z_INT = FFDom::ALL - FFDom::LB - FFDom::UB;

  oc.add_equation( PDE,    { t, z }, { T_INT,     Z_INT     }, io );
  oc.add_equation( BCINIT, { t, z }, { FFDom::LB, FFDom::ALL }, ii );
  oc.add_equation( BCLO,   { t, z }, { T_INT,     FFDom::LB  }, ib );
  oc.add_equation( BCHI,   { t, z }, { T_INT,     FFDom::UB  }, ib );

  FFVar L = OpEval( c, t, TAU );                        // c(tau,z) -- distributed LATCH over z
  oc.add_output( L,     { z }, { 1.0 } );               // Llin = L|z=1     (LINEAR LATCH)
  oc.add_output( L * L, { z }, { 1.0 } );               // Lsq  = L^2|z=1   (NONLINEAR LATCH)

  oc.set_evolution_domain( t );
  oc.options.REDUCE.ORDER   = OCFESLV::Options::RED_FULL;
  oc.options.SOLVE.MARCHING = march;
  oc.options.DISPLAY_LEVEL  = 0;
}

static std::vector<double> value_at( bool march, double aval, bool nu )
{
  FFGraph DAG; OCFESLV oc( &DAG ); FFVar t,z,c,a; build( DAG, oc, t, z, c, a, march, nu );
  if( !oc.setup() ) return {};
  std::vector<double> vi, ii; oc.init( vi, ii, nullptr );
  oc.set_input_values( a, { aval }, ii.data() );
  std::vector<double> xv = vi;
  if( !oc.solve( xv.data(), ii.data(), nullptr ).converged ) return {};
  return oc.val_functions();
}

static void run_value( bool march, double aval, bool nu )
{
  std::cout << "\n[" << ( march ? "MARCH" : "MONO" ) << " value]  a=" << aval << "\n";
  std::vector<double> F = value_at( march, aval, nu );
  if( F.size() < 2 ){ std::cout << "    solve/value FAILED\n"; ++g_fail; return; }
  double const Llin = 1.0 + aval*TAU, Lsq = Llin*Llin;
  ck( "Llin = L|z=1  (linear LATCH)",    F[0], Llin, 1e-6 );
  ck( "Lsq  = L^2|z=1 (nonlinear LATCH)", F[1], Lsq,  1e-6 );
}

static void run_sens( bool march, double aval, bool nu )
{
  std::cout << "\n[" << ( march ? "MARCH" : "MONO" ) << " sens]   a=" << aval << "\n";
  FFGraph DAG; OCFESLV oc( &DAG ); FFVar t,z,c,a; build( DAG, oc, t, z, c, a, march, nu );
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
  if( Jf.size() < 2 || Ja.size() < 2 ){ std::cout << "    jacobian too small\n"; ++g_fail; return; }

  double const h = 1e-6;
  std::vector<double> Fp = value_at( march, aval + h, nu ), Fm = value_at( march, aval - h, nu );
  double const fdlin = ( Fp.size()>=2 && Fm.size()>=2 ) ? ( Fp[0]-Fm[0] )/( 2*h ) : 0.;
  double const fdsq  = ( Fp.size()>=2 && Fm.size()>=2 ) ? ( Fp[1]-Fm[1] )/( 2*h ) : 0.;

  double const Llin = 1.0 + aval*TAU;
  double const dLlin = TAU, dLsq = 2.0 * Llin * TAU;
  ck( "fwd dLlin/da",  Jf[0], dLlin, 1e-5 );
  ck( "fwd dLsq/da",   Jf[1], dLsq,  1e-5 );
  ck( "adj dLlin == fwd", Ja[0], Jf[0], 1e-8 );
  ck( "adj dLsq  == fwd", Ja[1], Jf[1], 1e-8 );
  ck( "fwd dLlin == FD", Jf[0], fdlin, 1e-4 );
  ck( "fwd dLsq  == FD", Jf[1], fdsq,  1e-4 );
}

int main()
{
  std::cout << "===============================================================\n"
            << "  DISTRIBUTED LATCH captured reductions R(tau): value + sensitivity\n"
            << "  (validated on UNIFORM and NON-UNIFORM grids, in t AND z)\n"
            << "===============================================================\n";
  for( bool nu : { false, true } ){
    std::cout << "\n========== " << ( nu ? "NON-UNIFORM" : "UNIFORM" ) << " grid ==========\n";
    run_value( false, 1.2, nu );
    run_value( true,  1.2, nu );
    run_sens ( false, 1.2, nu );
    run_sens ( true,  1.2, nu );
  }
  std::cout << "\n  RESULT: " << ( g_fail==0 ? "ALL PASS" : "FAIL" ) << "  (" << g_fail << " failing)\n";
  return g_fail==0 ? 0 : 1;
}
