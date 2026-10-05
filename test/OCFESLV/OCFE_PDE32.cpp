// OCFE_PDE32_solve2.cpp   ---  PSA case study, RUNG 7a (MULTI-COMPONENT, const velocity)
// ===========================================================================
// First MULTI-COMPONENT rung.  A binary gas mixture (c1,c2) on a packed bed at
// PRESCRIBED CONSTANT velocity, competitive (multi-component) Langmuir uptake,
// and the ideal-gas equation of state on the TOTAL pressure  P = Rg (c1 + c2).
//
//   continuity_i: dci/dt + U dci/dz - D d2ci/dz2 + F dqi/dt = s_ci   (i=1,2)
//   LDF_i:        dqi/dt - k ( qi*(c1,c2) - qi )                = s_qi
//   EOS:          P - Rg ( c1 + c2 ) = 0          (algebraic; TOTAL pressure)
//   qi*(c1,c2) = qs_i b_i c_i / ( 1 + b1 c1 + b2 c2 )   (competitive Langmuir)
//
// WHY THIS RUNG.  Rung 6 (PDE30) carried P = Rg c, a value-slaved algebraic
// state that is an exact SCALAR MULTIPLE of a single state.  Here P is slaved to
// the SUM of two independently evolving components -- a "genuinely independent
// pressure" that is no longer a multiple of any one state.  Structurally it is
// still the same value-slaved case (derivative-free defining constraint, P the
// lone ALGEBRAIC bare state in the EOS while c1,c2 are dynamic), so the
// value-slaved eligibility-gate drop should flag P alone and IC_STRONG should
// recover the exact polynomial root to ~1e-13, exactly as rung 6.  This is the
// first rung where purity/recovery OUTPUTS are physically meaningful, so it also
// exercises the solve() output-overflow fix (_nCollFct > 0) for real.
//
// Staging: 7a here keeps velocity CONSTANT (no Darcy) to isolate the
// multi-component + value-slaved-on-a-sum + outputs machinery; 7b (PDE33) adds
// the pressure-Darcy momentum so P genuinely drives u.
//
// Manufactured (both ci use a (1-z)^2 profile so the zero-gradient outlet holds
// exactly; qi linear):
//   c1 = 1   + Ac1 (1-z)^2 + Cc1 t,   c2 = C20 + Ac2 (1-z)^2 + Cc2 t,
//   q1 = Q10 + Qz1 z + Qt1 t,         q2 = Q20 + Qz2 z + Qt2 t,
//   P  = Rg ( c1 + c2 ).
// Imposition: IC_STRONG (the PSA default), SPQR factorization, SIGMA0=10.
// ===========================================================================

#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <string>
#include <limits>

#ifndef OCFE_OCFESLV_HEADER
#define OCFE_OCFESLV_HEADER "ocfeslv.hpp"
#endif
#include OCFE_OCFESLV_HEADER

using namespace mc;

// transport / kinetics (single shared bed)
static double const U_VEL = 1.0, D_ax = 0.1, F_ph = 0.5, k_ldf = 2.0, Rg = 1.0;
// competitive Langmuir (distinct affinities so the two components compete)
static double const qs1 = 1.0, qs2 = 1.0, b1_L = 1.0, b2_L = 0.5;
// manufactured coefficients -- component 1
static double const Ac1 = 0.3, Cc1 = 0.2;
static double const Q10 = 0.4, Qz1 = 0.1, Qt1 = 0.3;
// manufactured coefficients -- component 2
static double const C20 = 0.5, Ac2 = 0.2, Cc2 = 0.15;
static double const Q20 = 0.25, Qz2 = 0.08, Qt2 = 0.12;

static inline double c1M_exact( double z, double t ){ return 1.0 + Ac1*(1.0-z)*(1.0-z) + Cc1*t; }
static inline double c2M_exact( double z, double t ){ return C20 + Ac2*(1.0-z)*(1.0-z) + Cc2*t; }
static inline double q1M_exact( double z, double t ){ return Q10 + Qz1*z + Qt1*t; }
static inline double q2M_exact( double z, double t ){ return Q20 + Qz2*z + Qt2*t; }
static inline double PM_exact ( double z, double t ){ return Rg*( c1M_exact(z,t) + c2M_exact(z,t) ); }

// analytic functional values for the manufactured solution
static double const c1_out_exact = c1M_exact( 1.0, 1.0 );                         // 1.20
static double const c2_out_exact = c2M_exact( 1.0, 1.0 );                         // 0.65
static double const P_out_exact  = PM_exact ( 1.0, 1.0 );                         // 1.85
static double const y1_out_exact = c1_out_exact / ( c1_out_exact + c2_out_exact );// 0.6486...
static double const F1_in_exact  = U_VEL*( ( 1.0 + Ac1 ) + 0.5*Cc1 );             // 1.40
static double const Q1_T_exact   = ( Q10 + Qt1 ) + 0.5*Qz1;                        // 0.75
static double const C1_T_exact   = ( 1.0 + Cc1 ) + Ac1/3.0;                        // 1.30

struct OutSpec { std::string name; double exact; };

static bool run_psa7a( FFDom::TYPE coltype, std::string const& cname )
{
  std::cout << "\n================================================================\n";
  std::cout << "  PSA rung 7a  multi-component, const velocity (EOS P=Rg(c1+c2))  " << cname << "\n";
  std::cout << "================================================================\n";

  size_t const n_el = 3, n_nd = 6;

  FFGraph DAG;
  FFVar t  = DAG.add_var( "t" );
  FFVar z  = DAG.add_var( "z" );
  FFVar c1 = DAG.add_var( "c1(t,z)" );
  FFVar c2 = DAG.add_var( "c2(t,z)" );
  FFVar q1 = DAG.add_var( "q1(t,z)" );
  FFVar q2 = DAG.add_var( "q2(t,z)" );
  FFVar P  = DAG.add_var( "P(t,z)" );

  FFPartial  OpP;
  FFIntegral OpI;

  // manufactured fields (for the sources)
  FFVar c1Man = 1.0 + Ac1*(1.0-z)*(1.0-z) + Cc1*t;
  FFVar c2Man = C20 + Ac2*(1.0-z)*(1.0-z) + Cc2*t;
  FFVar q1Man = Q10 + Qz1*z + Qt1*t;
  FFVar q2Man = Q20 + Qz2*z + Qt2*t;
  FFVar denMan = 1.0 + b1_L*c1Man + b2_L*c2Man;
  FFVar q1starMan = qs1*b1_L*c1Man/denMan;
  FFVar q2starMan = qs2*b2_L*c2Man/denMan;

  // sources (constant U: d(U ci)/dz = U dci/dz = U*(-2 Aci (1-z)))
  FFVar s_c1 = Cc1 - 2.0*U_VEL*Ac1*(1.0-z) - 2.0*D_ax*Ac1 + F_ph*Qt1;
  FFVar s_c2 = Cc2 - 2.0*U_VEL*Ac2*(1.0-z) - 2.0*D_ax*Ac2 + F_ph*Qt2;
  FFVar s_q1 = Qt1 - k_ldf*( q1starMan - q1Man );
  FFVar s_q2 = Qt2 - k_ldf*( q2starMan - q2Man );
  // Danckwerts inlet fluxes  g_i = U ci(0,t) - D dci/dz|_{z=0},  dci/dz|_0 = -2 Aci
  FFVar g1_in = U_VEL*( 1.0 + Ac1 + Cc1*t ) + 2.0*D_ax*Ac1;
  FFVar g2_in = U_VEL*( C20 + Ac2 + Cc2*t ) + 2.0*D_ax*Ac2;

  // governing residuals
  FFVar den   = 1.0 + b1_L*c1 + b2_L*c2;
  FFVar q1star = qs1*b1_L*c1/den;
  FFVar q2star = qs2*b2_L*c2/den;
  FFVar CONT1 = OpP( c1, t ) + U_VEL*OpP( c1, z ) - D_ax*OpP( OpP( c1, z ), z ) + F_ph*OpP( q1, t ) - s_c1;
  FFVar CONT2 = OpP( c2, t ) + U_VEL*OpP( c2, z ) - D_ax*OpP( OpP( c2, z ), z ) + F_ph*OpP( q2, t ) - s_c2;
  FFVar LDF1  = OpP( q1, t ) - k_ldf*( q1star - q1 ) - s_q1;
  FFVar LDF2  = OpP( q2, t ) - k_ldf*( q2star - q2 ) - s_q2;
  FFVar EOS   = P - Rg*( c1 + c2 );                    // value-slaved (P is the lone alg. bare state)
  FFVar IC_c1 = c1 - ( 1.0 + Ac1*(1.0-z)*(1.0-z) );
  FFVar IC_c2 = c2 - ( C20 + Ac2*(1.0-z)*(1.0-z) );
  FFVar IC_q1 = q1 - ( Q10 + Qz1*z );
  FFVar IC_q2 = q2 - ( Q20 + Qz2*z );
  FFVar BC_L1 = U_VEL*c1 - D_ax*OpP( c1, z ) - g1_in;  // Danckwerts inlet
  FFVar BC_L2 = U_VEL*c2 - D_ax*OpP( c2, z ) - g2_in;
  FFVar BC_U1 = OpP( c1, z );                          // zero-gradient outlet
  FFVar BC_U2 = OpP( c2, z );

  // output functionals
  FFVar PUR1 = c1/( c1 + c2 );          // outlet mole fraction (purity) of component 1
  FFVar F1in = OpI( U_VEL*c1, t );       // cumulative inlet feed of comp 1, integrated over t
  FFVar Q1T  = OpI( q1, z );             // bed loading of comp 1, integrated over z
  FFVar C1T  = OpI( c1, z );             // gas inventory of comp 1, integrated over z

  OCFESLV oc( &DAG );
  oc.add_domain( t, FFDom( 0., 1., n_el, coltype, n_nd ) );
  oc.add_domain( z, FFDom( 0., 1., n_el, coltype, n_nd ) );
  oc.add_state ( c1, {t,z} );
  oc.add_state ( c2, {t,z} );
  oc.add_state ( q1, {t,z} );
  oc.add_state ( q2, {t,z} );
  oc.add_state ( P,  {t,z} );
  oc.update_ref( c1, [&]( OCFESLV::t_Coord const& cr ){ return c1M_exact( cr.at(z), cr.at(t) ); } );
  oc.update_ref( c2, [&]( OCFESLV::t_Coord const& cr ){ return c2M_exact( cr.at(z), cr.at(t) ); } );
  oc.update_ref( q1, [&]( OCFESLV::t_Coord const& cr ){ return q1M_exact( cr.at(z), cr.at(t) ); } );
  oc.update_ref( q2, [&]( OCFESLV::t_Coord const& cr ){ return q2M_exact( cr.at(z), cr.at(t) ); } );
  oc.update_ref( P,  [&]( OCFESLV::t_Coord const& cr ){ return PM_exact ( cr.at(z), cr.at(t) ); } );
  oc.set_evolution_domain( t );

  int const T_NO_LB = FFDom::ALL - FFDom::LB;
  int const Z_INT   = FFDom::ALL - FFDom::LB - FFDom::UB;

  oc.add_equation( CONT1, {t,z}, {T_NO_LB,   Z_INT},      OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( CONT2, {t,z}, {T_NO_LB,   Z_INT},      OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( LDF1,  {t,z}, {T_NO_LB,   FFDom::ALL}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( LDF2,  {t,z}, {T_NO_LB,   FFDom::ALL}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( EOS,   {t,z}, {FFDom::ALL,FFDom::ALL}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( IC_c1, {t,z}, {FFDom::LB, FFDom::ALL}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INITIAL,  0 ) );
  oc.add_equation( IC_c2, {t,z}, {FFDom::LB, FFDom::ALL}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INITIAL,  0 ) );
  oc.add_equation( IC_q1, {t,z}, {FFDom::LB, FFDom::ALL}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INITIAL,  0 ) );
  oc.add_equation( IC_q2, {t,z}, {FFDom::LB, FFDom::ALL}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INITIAL,  0 ) );
  oc.add_equation( BC_L1, {t,z}, {T_NO_LB,   FFDom::LB},  OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );
  oc.add_equation( BC_L2, {t,z}, {T_NO_LB,   FFDom::LB},  OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );
  oc.add_equation( BC_U1, {t,z}, {T_NO_LB,   FFDom::UB},  OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );
  oc.add_equation( BC_U2, {t,z}, {T_NO_LB,   FFDom::UB},  OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );

  // scalar/integral outputs (order defines fct-row order)
  std::vector<OutSpec> ospec;
  oc.add_output( c1,   {z,t}, {1.0, 1.0} );  ospec.push_back( { "c1_out=c1(1,1)",      c1_out_exact } );
  oc.add_output( c2,   {z,t}, {1.0, 1.0} );  ospec.push_back( { "c2_out=c2(1,1)",      c2_out_exact } );
  oc.add_output( P,    {z,t}, {1.0, 1.0} );  ospec.push_back( { "P_out=P(1,1)",        P_out_exact  } );
  oc.add_output( PUR1, {z,t}, {1.0, 1.0} );  ospec.push_back( { "y1_out=c1/(c1+c2)|out",y1_out_exact } );
  oc.add_output( F1in, {z},   {0.0} );       ospec.push_back( { "F1_in=int(U c1)dt|z=0",F1_in_exact  } );
  oc.add_output( Q1T,  {t},   {1.0} );       ospec.push_back( { "Q1_T=int q1 dz|t=1",   Q1_T_exact   } );
  oc.add_output( C1T,  {t},   {1.0} );       ospec.push_back( { "C1_T=int c1 dz|t=1",   C1_T_exact   } );

  oc.options.REDUCE.ORDER     = OCFESLV::Options::RED_FULL;
  oc.options.CLASSIFY.MODE         = OCFESLV::Options::CLASS_AUTO;
  oc.options.INTERFACE.TYPE   = OCFESLV::Options::IC_AUTO;
  oc.options.INTERFACE.IMPOSITION  = OCFESLV::Options::IC_STRONG;   // PSA default; value-slaved P drop de-singularizes
  oc.options.INTERFACE.SAT_SIGMA0       = 10.0;                        // basin robustness
  oc.options.DISPLAY_LEVEL    = 1;                           // show classification + value-slaved drop on a NEW rung
#if defined(CRONOS__WITH_SPQR)
  oc.options.SOLVE.FACTORIZATION = OCFESLV::Options::SOLVE_SPQR;
#else
  oc.options.SOLVE.FACTORIZATION = OCFESLV::Options::SOLVE_SUPERLU;
#endif

  if( !oc.setup() ){ std::cerr << "ERROR: setup failed\n"; return false; }

  size_t const nVar = oc.n_colloc_sta();
  size_t const nEqn = oc.n_colloc_eqn();
  size_t const nFct = oc.n_colloc_fct();
  std::cout << "  type=" << OCFESLV::pde_type_name( oc.pde_type().type )
            << " nVar=" << nVar << " nEqn=" << nEqn << " nTrace=" << oc.n_colloc_trace()
            << " nFct=" << nFct << " square=" << (nVar==nEqn?"yes":"NO") << "\n";
  if( nVar != nEqn ){ std::cerr << "ERROR: not square\n"; return false; }
  if( nFct != ospec.size() ){
    std::cerr << "ERROR: nFct=" << nFct << " != " << ospec.size() << " scalar outputs\n";
    return false;
  }

  std::vector<double> varInit, inpInit;
  if( !oc.init( varInit, inpInit, nullptr ) ){ std::cerr << "ERROR: init failed\n"; return false; }

  std::vector<double> res( nEqn, 0. ), fct( nFct, 0. );
  if( oc.eval( res.data(), fct.data(), varInit.data(), nullptr, nullptr ) ){
    double rm=0.; for( double v : res ) rm = std::max( rm, std::fabs(v) );
    std::cout << "  [A] residual at manufactured exact: max|r|="
              << std::scientific << std::setprecision(4) << rm
              << "  " << ( rm<1e-9 ? "PASS" : "(not exact)" ) << "\n";
  }

  std::vector<double> xv = varInit;
  for( size_t i=0;i<xv.size();++i ) xv[i] += 0.05*std::sin( 0.7*double(i) + 0.2 );
  OCFESLV::SolveReport rep = oc.solve( xv.data(), nullptr, nullptr );
  std::cout << "  [B] solve: converged=" << (rep.converged?"yes":"no")
            << " iters=" << rep.iterations
            << " final|r|=" << std::scientific << std::setprecision(4) << rep.final_residual << "\n";
  if( !rep.converged ){ std::cerr << "ERROR: solve did not converge\n"; return false; }

  // state recovery
  double c1e=0., c2e=0., q1e=0., q2e=0., Pe=0.;
  double const sg[4] = { 0.15, 0.35, 0.65, 0.85 };
  for( double zs : sg ) for( double ts : sg ){
    OCFESLV::t_Coord pt; pt[z]=zs; pt[t]=ts;
    c1e = std::max( c1e, std::fabs( oc.eval_colloc<double>( c1, pt, xv.data(), nullptr, nullptr ) - c1M_exact(zs,ts) ) );
    c2e = std::max( c2e, std::fabs( oc.eval_colloc<double>( c2, pt, xv.data(), nullptr, nullptr ) - c2M_exact(zs,ts) ) );
    q1e = std::max( q1e, std::fabs( oc.eval_colloc<double>( q1, pt, xv.data(), nullptr, nullptr ) - q1M_exact(zs,ts) ) );
    q2e = std::max( q2e, std::fabs( oc.eval_colloc<double>( q2, pt, xv.data(), nullptr, nullptr ) - q2M_exact(zs,ts) ) );
    Pe  = std::max( Pe,  std::fabs( oc.eval_colloc<double>( P,  pt, xv.data(), nullptr, nullptr ) - PM_exact (zs,ts) ) );
  }
  std::cout << std::scientific << std::setprecision(4)
            << "  [C] state err c1=" << c1e << " c2=" << c2e
            << " q1=" << q1e << " q2=" << q2e << " P=" << Pe << "\n";

  // outputs at the converged solution
  std::fill( fct.begin(), fct.end(), 0. );
  // Output functionals from val_functions() (window-summed for evolution integrals) -- correct for
  // both monolithic and marching; a re-eval on xv would see only the last window.
  if( oc.val_functions().size() < nFct ){
    std::cerr << "ERROR: output functionals unavailable\n"; return false;
  }
  fct = oc.val_functions();

  std::cout << "  [D] purity/recovery outputs:\n";
  std::cout << "      " << std::left << std::setw(26) << "functional"
            << std::right << std::setw(16) << "computed"
            << std::setw(16) << "exact"
            << std::setw(13) << "abs err" << "  result\n";
  bool ok = ( c1e<1e-7 && c2e<1e-7 && q1e<1e-7 && q2e<1e-7 && Pe<1e-7 );
  for( size_t kk=0; kk<ospec.size(); ++kk ){
    size_t const r = oc.row_fct( kk );
    double const val = ( r==std::numeric_limits<size_t>::max() ) ? NAN : fct[ r - nEqn ];
    double const err = std::fabs( val - ospec[kk].exact );
    bool const pass = ( err < 1e-7 );
    ok &= pass;
    std::cout << "      " << std::left << std::setw(26) << ospec[kk].name
              << std::right << std::scientific << std::setprecision(8)
              << std::setw(16) << val << std::setw(16) << ospec[kk].exact
              << std::setprecision(3) << std::setw(13) << err
              << "  " << (pass?"PASS":"FAIL") << "\n";
  }

  std::cout << "  PSA7a " << cname << ": " << (ok?"PASS":"FAIL") << "\n";
  return ok;
}

int main()
{
  std::cout << "================================================================\n";
  std::cout << "  PSA rung 7a: multi-component, constant velocity\n";
  std::cout << "  binary mixture (c1,c2) + competitive Langmuir + EOS P=Rg(c1+c2)\n";
  std::cout << "  + purity/recovery outputs   (IC_STRONG, SPQR, SIGMA0=10)\n";
  std::cout << "================================================================\n";
  bool ok = true;
  ok &= run_psa7a( FFDom::CGL, "CGL" );
  ok &= run_psa7a( FFDom::LGL, "LGL" );
  std::cout << "\n  Overall: " << ( ok ? "ALL PASS" : "SOME FAILED" ) << "\n";
  return ok ? 0 : 1;
}
