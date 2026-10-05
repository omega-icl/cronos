// OCFE_PDE31b_chain.cpp  ---  reproducer step 2: the two-level algebraic chain
// ===========================================================================
// PDE31 (c + value-slaved P=Rg c) PASSED under STRONG/TRACE -> a value-slaved
// value-state is NOT the cause.  Rung 6's extra ingredient is u = -Kperm dP/dz
// with P = Rg c: a SECOND algebraic value-state slaved to the DERIVATIVE of the
// first.  This couples three redundancies into one cluster:
//     P-value == c-value      (P = Rg c)
//     P-flux  == c-flux        (d/dz of the above)
//     u-value == P-flux == c-flux   (u = -Kperm Rg dc/dz),
// and c-flux is independently claimed by the parabolic Dz_c machinery.
//
//   heat/conv:  dc/dt [+ d(u c)/dz] - D d2c/dz2 = s_c     (parabolic)
//   Darcy:      u + Kperm dP/dz = 0                        (algebraic value-state)
//   EOS:        P - Rg c = 0                               (algebraic value-state)
//   manufactured  c=1+a(1-z)^2+b t,  P=Rg c,  u=2 Kperm Rg a (1-z)
//
// Two bisection axes:
//   USE_CONV=false : u defined by Darcy but NOT used in continuity (chain present,
//                    decoupled) -> tests whether the CLAIM cluster alone is singular.
//   USE_CONV=true  : full rung-6-minus-q (u in d(u c)/dz) -> tests whether the
//                    convection coupling is required to tip it singular.
// Sweep STRONG / TRACE / WEAK for each.  WEAK is the always-pass control.
// ===========================================================================

#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <string>

#define MC__OCFESLV_STRONG_DUMMY_DIAG
#define MC__OCFESLV_SYMBOL_TAG_PROBE

#ifndef OCFE_OCFESLV_HEADER
#define OCFE_OCFESLV_HEADER "ocfeslv.hpp"
#endif
#include OCFE_OCFESLV_HEADER
using namespace mc;

static double const Kperm = 0.5, Rg = 1.0, D_ax = 0.1, a_c = 0.3, b_c = 0.2;
static double const kappa_g = Kperm*Rg;

static inline double cM( double z, double t ){ return 1.0 + a_c*(1.0-z)*(1.0-z) + b_c*t; }
static inline double PM( double z, double t ){ return Rg*cM(z,t); }
static inline double uM( double z, double   ){ return 2.0*kappa_g*a_c*(1.0-z); }

static void run( char const* tag, bool use_conv, int imposition )
{
  std::cout << "\n---- " << tag << " ----\n";
  size_t const n_el = 2, n_nd = 4;

  FFGraph DAG;
  FFVar t=DAG.add_var("t"), z=DAG.add_var("z");
  FFVar c=DAG.add_var("c(t,z)"), P=DAG.add_var("P(t,z)"), u=DAG.add_var("u(t,z)");
  FFPartial OpP;

  FFVar s_c = use_conv
            ? ( b_c - 2.0*kappa_g*a_c*( 1.0 + 3.0*a_c*(1.0-z)*(1.0-z) + b_c*t ) - D_ax*( 2.0*a_c ) )
            : ( b_c - D_ax*( 2.0*a_c ) );
  FFVar conv = use_conv ? ( u*OpP(c,z) + c*OpP(u,z) ) : FFVar( 0.0 );
  FFVar HEAT = OpP(c,t) + conv - D_ax*OpP(OpP(c,z),z) - s_c;
  FFVar DARCY = u + Kperm*OpP(P,z);
  FFVar EOS   = P - Rg*c;
  FFVar IC_c  = c - ( 1.0 + a_c*(1.0-z)*(1.0-z) );
  FFVar BC_Lc = c - ( 1.0 + a_c + b_c*t );
  FFVar BC_Uc = OpP(c,z);

  OCFESLV oc( &DAG );
  oc.add_domain( t, FFDom( 0., 1., n_el, FFDom::LGL, n_nd ) );
  oc.add_domain( z, FFDom( 0., 1., n_el, FFDom::LGL, n_nd ) );
  oc.add_state( c, {t,z} ); oc.add_state( P, {t,z} ); oc.add_state( u, {t,z} );
  oc.update_ref( c, [&](OCFESLV::t_Coord const& cr){ return cM(cr.at(z),cr.at(t)); } );
  oc.update_ref( P, [&](OCFESLV::t_Coord const& cr){ return PM(cr.at(z),cr.at(t)); } );
  oc.update_ref( u, [&](OCFESLV::t_Coord const& cr){ return uM(cr.at(z),cr.at(t)); } );
  oc.set_evolution_domain( t );

  int const T_NO_LB = FFDom::ALL - FFDom::LB;
  int const Z_INT   = FFDom::ALL - FFDom::LB - FFDom::UB;
  oc.add_equation( HEAT,  {t,z}, {T_NO_LB,   Z_INT},      OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( DARCY, {t,z}, {FFDom::ALL,FFDom::ALL}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( EOS,   {t,z}, {FFDom::ALL,FFDom::ALL}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( IC_c,  {t,z}, {FFDom::LB, FFDom::ALL}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INITIAL,  0 ) );
  oc.add_equation( BC_Lc, {t,z}, {T_NO_LB,   FFDom::LB},  OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );
  oc.add_equation( BC_Uc, {t,z}, {T_NO_LB,   FFDom::UB},  OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );

  oc.options.REDUCE.ORDER     = OCFESLV::Options::RED_FULL;
  oc.options.CLASSIFY.MODE         = OCFESLV::Options::CLASS_AUTO;
  oc.options.INTERFACE.TYPE   = OCFESLV::Options::IC_AUTO;
  oc.options.INTERFACE.IMPOSITION  = static_cast<OCFESLV::Options::ImpositionType>( imposition );
  oc.options.INTERFACE.SAT_SIGMA0       = 10.0;
  oc.options.DISPLAY_LEVEL    = 1;

  if( !oc.setup() ){ std::cerr << "  setup failed\n"; return; }
  size_t const nVar = oc.n_colloc_sta(), nEqn = oc.n_colloc_eqn();
  std::cout << "  nVar=" << nVar << " nEqn=" << nEqn << " nTrace=" << oc.n_colloc_trace()
            << " square=" << (nVar==nEqn?"yes":"NO") << "\n";
  if( nVar != nEqn ) return;

  std::vector<double> varInit, inpInit;
  if( !oc.init( varInit, inpInit, nullptr ) ){ std::cerr << "  init failed\n"; return; }
  std::vector<double> res( nEqn, 0. );
  double Aerr=0.;
  if( oc.eval( res.data(), nullptr, varInit.data(), nullptr, nullptr ) )
    for( double v : res ) Aerr = std::max( Aerr, std::fabs(v) );

  std::vector<double> xv = varInit;
  for( size_t i=0;i<xv.size();++i ) xv[i] += 0.05*std::sin( 0.7*double(i)+0.2 );
  OCFESLV::SolveReport rep = oc.solve( xv.data(), nullptr, nullptr );

  double cerr=0., uerr=0., Perr=0.;
  double const sg[3]={0.2,0.5,0.8};
  for( double zs : sg ) for( double ts : sg ){
    OCFESLV::t_Coord pt; pt[z]=zs; pt[t]=ts;
    cerr=std::max(cerr,std::fabs(oc.eval_colloc<double>(c,pt,xv.data(),nullptr,nullptr)-cM(zs,ts)));
    uerr=std::max(uerr,std::fabs(oc.eval_colloc<double>(u,pt,xv.data(),nullptr,nullptr)-uM(zs,ts)));
    Perr=std::max(Perr,std::fabs(oc.eval_colloc<double>(P,pt,xv.data(),nullptr,nullptr)-PM(zs,ts)));
  }
  std::cout << std::scientific << std::setprecision(3)
            << "  [A]=" << Aerr << "  conv=" << (rep.converged?"y":"n")
            << " it=" << rep.iterations << " |r|=" << rep.final_residual
            << "  [C] c=" << cerr << " u=" << uerr << " P=" << Perr
            << "  " << ( cerr<1e-9 && uerr<1e-9 && Perr<1e-9 ? "PASS" : "fail" ) << "\n";
}

int main()
{
  std::cout << "================================================================\n";
  std::cout << "  Two-level chain reproducer: c -> P=Rg c -> u=-Kperm dP/dz\n";
  std::cout << "  (compile with -DMC__OCFESLV_STRONG_DUMMY_DIAG for the tau partition)\n";
  std::cout << "================================================================\n";

  std::cout << "\n=== AXIS A: chain present, u DECOUPLED from continuity (no convection) ===\n";
  run( "(1) no-conv  STRONG", false, OCFESLV::Options::IC_STRONG );
  run( "(2) no-conv  WEAK",   false, OCFESLV::Options::IC_WEAK   );

  std::cout << "\n=== AXIS B: full rung-6-minus-q, u in d(u c)/dz (with convection) ===\n";
  run( "(3) conv     STRONG", true,  OCFESLV::Options::IC_STRONG );
  run( "(4) conv     TRACE",  true,  OCFESLV::Options::IC_TRACE  );
  run( "(5) conv     WEAK",   true,  OCFESLV::Options::IC_WEAK   );

  std::cout << "\n  Read: if (1) already fails -> the chain's CLAIM cluster is singular by itself.\n";
  std::cout << "        if (1) passes but (3) fails -> the convection coupling tips it.\n";
  return 0;
}
