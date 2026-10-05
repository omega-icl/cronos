// OCFE_PDE31_eosmin.cpp  ---  minimal reproducer for the IC_STRONG / IC_TRACE
//                             value-slaved-state continuity singularity
// ===========================================================================
// Strips rung 6 to its essence: a single DYNAMIC state c (heat equation) and a
// single ALGEBRAIC state P value-slaved to c by a derivative-free EOS, P=Rg c.
// No q, no u, no Darcy.  Three configs isolate the cause:
//
//   (1) c only,  IC_STRONG          -> baseline; must be non-singular (PASS).
//   (2) c + P,   IC_STRONG          -> add the value-slaved P; singular if the
//                                       implied P-continuity claim is the cause.
//   (3) c + P,   IC_WEAK            -> no tau block; must recover (PASS).
//
//   heat:  dc/dt - D d2c/dz2 = s_c           (parabolic; gets 2 BCs + IC)
//   EOS:   P - Rg c = 0                       (algebraic; value-slaved, no IC/BC)
//   manufactured  c = 1 + a(1-z)^2 + b t ,  P = Rg c   (degree-2, exact-rep)
//
// Hypothesis: P's value-continuity (P_L=P_R) is exactly Rg*(c_L=c_R) THROUGH the
// EOS, but as a RAW claim it distributes to different receiver rows than c's, so
// the receiver-coupling redundancy proxy never flags it.  Run also with
//   -DMC__OCFESLV_STRONG_DUMMY_DIAG
// to print the tau column-dependence partition: the prediction is rank(all rows)
// == ntrace (NO genuine column-dependence seen) while the solve is still singular
// -- i.e. the dependence is NOT in the receiver coupling the detector inspects.
// ===========================================================================

#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <string>

#ifndef OCFE_OCFESLV_HEADER
#define OCFE_OCFESLV_HEADER "ocfeslv.hpp"
#endif
#include OCFE_OCFESLV_HEADER
using namespace mc;

static double const D_ax = 0.1, Rg = 1.0, a_c = 0.3, b_c = 0.2;

static inline double cM( double z, double t ){ return 1.0 + a_c*(1.0-z)*(1.0-z) + b_c*t; }
static inline double PM( double z, double t ){ return Rg*cM(z,t); }

static void run( char const* tag, bool with_P, int imposition )
{
  std::cout << "\n---- " << tag << " ----\n";

  size_t const n_el = 2, n_nd = 4;
  FFGraph DAG;
  FFVar t=DAG.add_var("t"), z=DAG.add_var("z");
  FFVar c=DAG.add_var("c(t,z)");
  FFVar P = with_P ? DAG.add_var("P(t,z)") : FFVar();
  FFPartial OpP;

  FFVar s_c  = b_c - D_ax*( 2.0*a_c );
  FFVar HEAT = OpP(c,t) - D_ax*OpP(OpP(c,z),z) - s_c;
  FFVar EOS  = with_P ? ( P - Rg*c ) : FFVar();
  FFVar IC_c = c - ( 1.0 + a_c*(1.0-z)*(1.0-z) );
  FFVar BC_Lc = c - ( 1.0 + a_c + b_c*t );      // Dirichlet inlet
  FFVar BC_Uc = OpP(c,z);                         // zero-gradient outlet

  OCFESLV oc( &DAG );
  oc.add_domain( t, FFDom( 0., 1., n_el, FFDom::LGL, n_nd ) );
  oc.add_domain( z, FFDom( 0., 1., n_el, FFDom::LGL, n_nd ) );
  oc.add_state( c, {t,z} );
  oc.update_ref( c, [&](OCFESLV::t_Coord const& cr){ return cM(cr.at(z),cr.at(t)); } );
  if( with_P ){
    oc.add_state( P, {t,z} );
    oc.update_ref( P, [&](OCFESLV::t_Coord const& cr){ return PM(cr.at(z),cr.at(t)); } );
  }
  oc.set_evolution_domain( t );

  int const T_NO_LB = FFDom::ALL - FFDom::LB;
  int const Z_INT   = FFDom::ALL - FFDom::LB - FFDom::UB;
  oc.add_equation( HEAT,  {t,z}, {T_NO_LB,   Z_INT},      OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  if( with_P )
  oc.add_equation( EOS,   {t,z}, {FFDom::ALL,FFDom::ALL}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( IC_c,  {t,z}, {FFDom::LB, FFDom::ALL}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INITIAL,  0 ) );
  oc.add_equation( BC_Lc, {t,z}, {T_NO_LB,   FFDom::LB},  OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );
  oc.add_equation( BC_Uc, {t,z}, {T_NO_LB,   FFDom::UB},  OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );

  oc.options.REDUCE.ORDER     = OCFESLV::Options::RED_FULL;
  oc.options.CLASSIFY.MODE         = OCFESLV::Options::CLASS_AUTO;
  oc.options.INTERFACE.TYPE   = OCFESLV::Options::IC_AUTO;
  oc.options.INTERFACE.IMPOSITION  = static_cast<OCFESLV::Options::ImpositionType>( imposition );
  oc.options.INTERFACE.SAT_SIGMA0       = 10.0;
  oc.options.DISPLAY_LEVEL    = 1;   // show redundancy pre-pass + B_S rcond

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

  double cerr=0., Perr=0.;
  double const sg[3]={0.2,0.5,0.8};
  for( double zs : sg ) for( double ts : sg ){
    OCFESLV::t_Coord pt; pt[z]=zs; pt[t]=ts;
    cerr=std::max(cerr,std::fabs(oc.eval_colloc<double>(c,pt,xv.data(),nullptr,nullptr)-cM(zs,ts)));
    if( with_P )
      Perr=std::max(Perr,std::fabs(oc.eval_colloc<double>(P,pt,xv.data(),nullptr,nullptr)-PM(zs,ts)));
  }
  std::cout << std::scientific << std::setprecision(3)
            << "  [A]=" << Aerr << "  conv=" << (rep.converged?"y":"n")
            << " it=" << rep.iterations << " |r|=" << rep.final_residual
            << "  [C] c=" << cerr << ( with_P ? "" : "" );
  if( with_P ) std::cout << " P=" << Perr;
  std::cout << "  " << ( cerr<1e-9 && (!with_P||Perr<1e-9) ? "PASS" : "fail" ) << "\n";
}

int main()
{
  std::cout << "================================================================\n";
  std::cout << "  EOS-redundancy minimal reproducer: dynamic c + value-slaved P=Rg c\n";
  std::cout << "  (run also with -DMC__OCFESLV_STRONG_DUMMY_DIAG for the tau partition)\n";
  std::cout << "================================================================\n";
  run( "(1) c only        IC_STRONG  -- baseline, expect PASS",  false, OCFESLV::Options::IC_STRONG );
  run( "(2) c + P=Rg c    IC_STRONG  -- expect singular/fail",   true,  OCFESLV::Options::IC_STRONG );
  run( "(3) c + P=Rg c    IC_TRACE   -- expect singular/fail",   true,  OCFESLV::Options::IC_TRACE  );
  run( "(4) c + P=Rg c    IC_WEAK    -- expect PASS",            true,  OCFESLV::Options::IC_WEAK   );
  std::cout << "\n  If (1)PASS (2,3)fail (4)PASS: value-slaved continuity claim is the sole cause.\n";
  return 0;
}
