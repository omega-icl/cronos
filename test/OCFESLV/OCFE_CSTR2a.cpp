// OCFE_CSTR_march.cpp  ---  MARCHING the isothermal CSTR DAE vs monolithic (stage 2a)
// ===========================================================================
// A well-mixed (0-D space, time-only) isothermal CSTR modelled as an index-1 DAE:
//
//   differential   dc/dt = D*(cf - c) - r          (reactant balance, feed cf)
//   algebraic      0     = r - k*c/(1 + K*c)        (quasi-steady reaction rate)
//
// The reactor starts AT STEADY STATE for the pre-step feed cf0, then the feed
// steps to cf immediately after start.  The distinction between the t=0 input
// (cf0) and the t>0 input (cf) is expressed at the equation level:
//
//   INITIAL (t=0, LB) :  0 = D*(cf0 - c0) - r0     (steady state at cf0; dc/dt=0)
//   EVOLUTION (t>0)   :  dc/dt = D*(cf - c) - r     (dynamics at cf)
//   ALGEBRAIC (all t) :  r = k*c/(1+K*c)
//
// The monolithic solve therefore CONSISTENTLY INITIALISES the DAE at the LB:
// c0 and the algebraic r0 are solved from the steady-state condition at cf0
// (not prescribed), then the transient relaxes toward the cf steady state.
//
//   C0  solve converged
//   C1  LB state c(0) == analytic steady state at cf0            [consistent init]
//   C2  algebraic r(0) == k*c0/(1+K*c0)                          [algebraic consistency at LB]
//   C3  terminal c(T_end) ~ analytic steady state at cf          [relaxation]
//   C4  reduced Jacobian forward == adjoint                      [tight]
//   C5  b/cf columns vs central-difference FD                    [loose]
//   C6  d c(0)/d cf == 0 and d r(0)/d cf == 0                    [IC is causally upstream of cf]
//   Marching mechanics for the reactor DAE: per-window algebraic re-solve of r, value
//   continuity of c across windows, vs the monolithic reference.  The cf0 steady state is
//   the (precomputed) initial value here; option (b) will later compute it in-framework.
//   M0 marched == monolithic ; M1 c(0)==cf0 ss ; M2 r(0) algebraic ; M3 c(T)~cf ss ;
//   M4 fwd==adj ; M5 FD ; M6 dc(0)/dcf==0
// ===========================================================================
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <string>
#include <sstream>

#ifndef OCFE_OCFESLV_HEADER
#define OCFE_OCFESLV_HEADER "ocfeslv.hpp"
#endif
#include OCFE_OCFESLV_HEADER
#include "ffocfe.hpp"

using namespace mc;

// --- CSTR parameters ---
static double const Dil = 1.0;    // dilution rate q/V
static double const krate = 2.0;  // rate constant
static double const Ksat = 1.0;   // saturation constant
static double const cf0_nom = 1.0; // pre-step feed  (t=0 steady state)
static double const cf_nom  = 2.0; // post-step feed (t>0 dynamics)
static double const T_end = 5.0;

// analytic steady state c_ss(u): solve  Dil*(u - c) - krate*c/(1+Ksat*c) = 0  (Newton)
static double c_steady( double u )
{
  double c = 0.5*u;
  for( int it=0; it<100; ++it ){
    double f  = Dil*( u - c ) - krate*c/( 1.0 + Ksat*c );
    double df = -Dil - krate*( 1.0 )/( ( 1.0 + Ksat*c )*( 1.0 + Ksat*c ) );
    double dc = -f/df; c += dc;
    if( std::fabs( dc ) < 1e-14 ) break;
  }
  return c;
}
static double r_of_c( double c ){ return krate*c/( 1.0 + Ksat*c ); }

static int g_pass = 0, g_fail = 0;
static void check_true( char const* nm, bool ok )
{ std::cout << "  " << std::left << std::setw(46) << nm << ( ok ? " PASS" : " FAIL" ) << "\n"; ok ? ++g_pass : ++g_fail; }
static void check_close( char const* nm, double got, double want, double tol )
{ double e = std::fabs( got - want );
  std::cout << "  " << std::left << std::setw(46) << nm << " |got-want|=" << std::scientific << std::setprecision(3) << e
            << " tol=" << tol << ( e <= tol ? "  PASS" : "  FAIL" ) << std::defaultfloat << "\n";
  ( e <= tol ) ? ++g_pass : ++g_fail; }

int main()
{
  std::cout << "================================================================\n"
            << "  MARCHING isothermal CSTR DAE vs monolithic (stage 2a)\n"
            << "  steady state at cf0=" << cf0_nom << " -> step to cf=" << cf_nom << "\n"
            << "================================================================\n\n";

  size_t const n_el = 6, n_nd = 5;
  FFDom::TYPE const coltype = FFDom::LGL;

  FFGraph DAG;
  FFVar t   = DAG.add_var( "t" );
  FFVar c   = DAG.add_var( "c(t)" );    // differential state
  FFVar r   = DAG.add_var( "r(t)" );    // algebraic state
  FFVar c_ic = DAG.add_var( "c_ic" );   // initial value of c (= cf0 steady state) -- transfer/decision
  FFVar cf   = DAG.add_var( "cf"  );    // dynamics feed                            -- decision

  OCFESLV oc( &DAG );
  oc.add_domain( t, FFDom( 0., T_end, n_el, coltype, n_nd ) );
  oc.set_evolution_domain( t );   // REQUIRED for marching: the collapse runs before auto-inference

  oc.add_state( c, {t} );
  oc.add_state( r, {t} );
  double const css0 = c_steady( cf0_nom );              // consistent initial state (precomputed here)
  oc.add_input( c_ic, css0,   /*is_decision=*/true );   // IC value -> value continuity across windows
  oc.add_input( cf,   cf_nom, /*is_decision=*/true );   // dynamics feed

  // reference profiles (initial guess): start everywhere near the cf0 steady state
  oc.update_ref( c, [&]( OCFESLV::t_Coord const& ){ return css0; } );
  oc.update_ref( r, [&]( OCFESLV::t_Coord const& ){ return r_of_c( css0 ); } );

  FFPartial OpP;   // partial-derivative operator (d/dt); builds the DAG evolution derivative

  FFVar EVOL = OpP( c, t ) - ( Dil*( cf  - c ) - r );          // t>0 dynamics (cf)
  FFVar ALG  = r - krate*c/( 1.0 + Ksat*c );                   // algebraic, all t
  FFVar IC   = c - c_ic;                                       // value continuity (cf0 ss at k=0; terminal at k>0)

  int const T_NO_LB = FFDom::ALL - FFDom::LB;
  oc.add_equation( EVOL, {t}, {T_NO_LB},   OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( ALG,  {t}, {FFDom::ALL},OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( IC,   {t}, {FFDom::LB}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INITIAL,  0 ) );

  oc.add_output( c, {t}, {0.0}   );   // F0 = c(0)      (initial steady state)
  oc.add_output( c, {t}, {T_end} );   // F1 = c(T_end)  (relaxed state)
  oc.add_output( r, {t}, {0.0}   );   // F2 = r(0)      (initial algebraic)

  if( !oc.setup() ){ std::cerr << "  setup() FAILED\n"; return 1; }

  std::vector<double> xv, inp;
  if( !oc.init( xv, inp, nullptr ) ){ std::cerr << "  init() FAILED\n"; return 1; }
  size_t const ncd = oc.n_control_dof();
  size_t const ncf = oc.n_colloc_fct();
  auto const& C = oc.controls();
  size_t const icic = C.at( c_ic ).offset, icf = C.at( cf ).offset;
  std::vector<double> p0; oc.encode_controls( inp.data(), p0 );

  std::cout << "  is_marching()=" << ( oc.is_marching() ? "true" : "false" )
            << " ;  n_march_steps=" << oc.n_march_steps()
            << " ;  ncd=" << ncd << " (c_ic@" << icic << ", cf@" << icf << ")  ncf=" << ncf << "\n\n";

  check_true( "M-pre is_marching()", oc.is_marching() );

  // ---- MARCHED solve (value continuity of c; r re-solved per window) ----
  std::vector<double> xm( xv ), im( inp );
  OCFESLV::SolveReport const repM = oc.solve( xm.data(), im.data(), nullptr );
  std::vector<double> Fm = oc.val_functions();
  check_true( "M0a marched DAE solve converged", repM.converged );
  if( !repM.converged ){ std::cerr << "  abort\n"; return 1; }

  // ---- MONOLITHIC reference (same model, SOLVE_MARCHING off) ----
  oc.options.SOLVE.MARCHING = false;
  std::vector<double> xo( xv ), io( inp );
  OCFESLV::SolveReport const repO = oc.solve( xo.data(), io.data(), nullptr );
  std::vector<double> Fo = oc.val_functions();
  oc.options.SOLVE.MARCHING = true;
  check_true( "M0b monolithic reference converged", repO.converged );

  // ---- M0: marched == monolithic ----
  { double e=0.; for( size_t a=0;a<ncf && a<Fm.size() && a<Fo.size();++a ) e=std::max(e,std::fabs(Fm[a]-Fo[a]));
    check_close( "M0 marched functions == monolithic", e, 0., 1e-8 ); }

  double const c0 = Fm[0], cT = Fm[1], r0 = Fm[2];
  // ---- M1: LB is the cf0 steady state (value continuity IC) ----
  check_close( "M1 c(0) == cf0 steady state", c0, css0, 1e-8 );
  // ---- M2: algebraic variable consistent at the LB (per-window re-solve) ----
  check_close( "M2 r(0) == k c0/(1+K c0)", r0, r_of_c( c0 ), 1e-10 );
  // ---- M3: transient relaxes toward the cf steady state ----
  check_close( "M3 c(T_end) ~ cf steady state", cT, c_steady( cf_nom ), 1e-3 );

  // ---- reduced Jacobian through the march: forward and adjoint ----
  std::vector<std::vector<double>> Jf( ncf, std::vector<double>( ncd, 0. ) ),
                                   Ja( ncf, std::vector<double>( ncd, 0. ) );
  { std::vector<double> x( xv ), i( inp );
    if( !oc.solve_fsens( x.data(), i.data(), nullptr ) ){ std::cerr << "  solve_fsens FAILED\n"; return 1; }
    auto const& J = oc.sens_jacobian(); for( size_t a=0;a<ncf;++a) for(size_t b=0;b<ncd;++b) Jf[a][b]=J[a*ncd+b]; }
  { std::vector<double> x( xv ), i( inp );
    if( !oc.solve_asens( x.data(), i.data(), nullptr ) ){ std::cerr << "  solve_asens FAILED\n"; return 1; }
    auto const& J = oc.sens_jacobian(); for( size_t a=0;a<ncf;++a) for(size_t b=0;b<ncd;++b) Ja[a][b]=J[a*ncd+b]; }
  { double e=0.; for(size_t a=0;a<ncf;++a) for(size_t b=0;b<ncd;++b) e=std::max(e,std::fabs(Jf[a][b]-Ja[a][b]));
    check_close( "M4 marched Jacobian forward == adjoint", e, 0., 1e-9 ); }

  // ---- M5: FD check of both control columns (through the march) ----
  {
    double const h = 1e-6;
    for( size_t col=0; col<ncd; ++col ){
      std::vector<double> pp( p0 ), pm2( p0 ); pp[col]+=h; pm2[col]-=h;
      std::vector<double> ip( inp ), im2( inp ), xp( xv ), xm2( xv );
      oc.decode_controls( pp,  ip.data()  ); oc.solve( xp.data(),  ip.data(),  nullptr ); std::vector<double> Fp = oc.val_functions();
      oc.decode_controls( pm2, im2.data() ); oc.solve( xm2.data(), im2.data(), nullptr ); std::vector<double> Fmn = oc.val_functions();
      for( size_t a=0; a<ncf; ++a ){
        double fd = ( Fp[a]-Fmn[a] )/( 2.0*h );
        double tol = 1e-2*std::max( 1e-4, std::fabs( Jf[a][col] ) ) + 1e-5;
        std::ostringstream nm; nm << "M5 dF" << a << "/dp" << col << " ~ FD";
        check_close( nm.str().c_str(), Jf[a][col], fd, tol );
      }
    }
  }

  // ---- M6: c(0) depends on c_ic, not on the dynamics feed cf ----
  check_close( "M6 d c(0)/d cf == 0", Jf[0][icf], 0., 1e-9 );
  std::cout << "\n============================================================\n"
            << "  RESULT: " << g_pass << " passed, " << g_fail << " failed  -- "
            << ( g_fail==0 ? "ALL PASS" : "SOME FAILED" ) << "\n"
            << "============================================================\n";
  return g_fail ? 1 : 0;
}
