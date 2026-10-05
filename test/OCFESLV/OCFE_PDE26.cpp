// OCFE_PDE26_solve2.cpp   ---  PSA case study, RUNG 3
// ===========================================================================
// Same isothermal single-component Darcy-velocity bed as rung 2, now exercising
// SCALAR OUTPUT FUNCTIONALS (the rung that uses solve()'s output path):
//
//   point  :  c_mid = c(z=0.5, t=1)        c_out = c(z=1, t=1)
//   integral: F_in  = [ \int_0^1 (u c) dt ]_{z=0}   (inlet cumulative feed)
//             Q_T   = [ \int_0^1 q dz ]_{t=1}       (bed loading at final time)
//             C_T   = [ \int_0^1 c dz ]_{t=1}       (gas inventory at final time)
//
// Integrals are FFIntegral (OpI) quadratures over a collocation direction; each
// output is evaluated at a fixed point in the remaining direction.  Outputs land
// in the eval() `fct` buffer (NOT the equation set) and are validated against the
// analytic values of the manufactured polynomials.
//
// PSA standard from here on: IMPOSITION_TYPE = IC_STRONG (structural continuity,
// robust to initial guess on the quadratic-gradient nonlinearity) and a bumped
// SAT_SIGMA0 = 10 for margin.
//
//   gas:  dc/dt + u dc/dz + c du/dz - D d2c/dz2 + F dq/dt = s_c
//   LDF:  dq/dt - k ( qs b c/(1+b c) - q )                = s_q
//   mom:  u + K dc/dz = 0
// Manufactured: c = 1 + (Az/2)(1-z)^2 + Ct t,  q = Q0+Qz z+Qt t,  u = K Az (1-z).
// ===========================================================================

#include <iostream>
#include <cstdlib>
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

static double const D_ax = 0.1, F_ph = 0.5, k_ldf = 2.0, qs_L = 1.0, b_L = 1.0;
static double const K_dar = 2.0;
static double const Az = 0.5, Ct = 0.2;
static double const Q0 = 0.4, Qz = 0.1, Qt = 0.3;

static inline double cM_exact( double z, double t ){ return 1.0 + 0.5*Az*(1.0-z)*(1.0-z) + Ct*t; }
static inline double qM_exact( double z, double t ){ return Q0 + Qz*z + Qt*t; }
static inline double uM_exact( double z, double   ){ return K_dar*Az*(1.0-z); }

// Analytic functional values for the manufactured solution.
//   F_in = int_0^1 u(0,t) c(0,t) dt,  u(0,t)=K Az,  c(0,t)=1+Az/2+Ct t
static double const F_in_exact = K_dar*Az*( 1.0 + 0.5*Az ) + K_dar*Az*Ct/2.0;   // 1.35
static double const Q_T_exact  = ( Q0 + Qt ) + Qz/2.0;                          // 0.75
static double const C_T_exact  = ( 1.0 + Ct ) + 0.5*Az/3.0;                     // 1.283333...

struct OutSpec { std::string name; double exact; };

static bool run_psa3( FFDom::TYPE coltype, std::string const& cname , OCFESLV::Options::ImpositionType imp = OCFESLV::Options::IC_STRONG,
                          OCFESLV::Options::ReductionType red = OCFESLV::Options::RED_FULL,
                          char const* tag = nullptr )
{
  std::cout << "\n================================================================\n";
  std::cout << "  PSA rung 3  scalar outputs (point + integral)  " << cname
            << ( tag ? std::string( "  [" ) + tag + "]" : std::string() ) << "\n";
  std::cout << "================================================================\n";

  size_t const n_el = 3, n_nd = 6;

  FFGraph DAG;
  FFVar t = DAG.add_var( "t" );
  FFVar z = DAG.add_var( "z" );
  FFVar c = DAG.add_var( "c(t,z)" );
  FFVar q = DAG.add_var( "q(t,z)" );
  FFVar u = DAG.add_var( "u(t,z)" );

  FFPartial  OpP;
  FFIntegral OpI;

  FFVar cMan  = 1.0 + 0.5*Az*(1.0-z)*(1.0-z) + Ct*t;
  FFVar qMan  = Q0 + Qz*z + Qt*t;
  FFVar uMan  = K_dar*Az*(1.0-z);
  FFVar dzc   = -Az*(1.0-z);
  FFVar qstarM = qs_L*b_L*cMan / ( 1.0 + b_L*cMan );
  FFVar s_c   = Ct + uMan*dzc + cMan*( -K_dar*Az ) - D_ax*Az + F_ph*Qt;
  FFVar s_q   = Qt - k_ldf*( qstarM - qMan );
  FFVar g_in  = K_dar*Az*( 1.0 + 0.5*Az + Ct*t ) + D_ax*Az;

  FFVar PDE_c = OpP( c, t ) + u*OpP( c, z ) + c*OpP( u, z ) - D_ax*OpP( OpP( c, z ), z ) + F_ph*OpP( q, t ) - s_c;
  FFVar LDF_q = OpP( q, t ) - k_ldf*( qs_L*b_L*c/( 1.0 + b_L*c ) - q ) - s_q;
  FFVar MOM   = u + K_dar*OpP( c, z );
  FFVar IC_c  = c - ( 1.0 + 0.5*Az*(1.0-z)*(1.0-z) );
  FFVar IC_q  = q - ( Q0 + Qz*z );
  FFVar BC_L  = u*c - D_ax*OpP( c, z ) - g_in;
  FFVar BC_U  = OpP( c, z );

  // Output functionals.
  FFVar F_in = OpI( u*c, t );   // integrate inlet flux over t -> function of z
  FFVar Q_T  = OpI( q, z );     // integrate loading over z   -> function of t
  FFVar C_T  = OpI( c, z );     // integrate concentration over z -> function of t

  OCFESLV oc( &DAG );
  oc.add_domain( t, FFDom( 0., 1., n_el, coltype, n_nd ) );
  oc.add_domain( z, FFDom( 0., 1., n_el, coltype, n_nd ) );
  oc.add_state ( c, {t,z} );
  oc.add_state ( q, {t,z} );
  oc.add_state ( u, {t,z} );
  oc.update_ref( c, [&]( OCFESLV::t_Coord const& cr ){ return cM_exact( cr.at(z), cr.at(t) ); } );
  oc.update_ref( q, [&]( OCFESLV::t_Coord const& cr ){ return qM_exact( cr.at(z), cr.at(t) ); } );
  oc.update_ref( u, [&]( OCFESLV::t_Coord const& cr ){ return uM_exact( cr.at(z), cr.at(t) ); } );
  oc.set_evolution_domain( t );

  int const T_NO_LB = FFDom::ALL - FFDom::LB;
  int const Z_INT   = FFDom::ALL - FFDom::LB - FFDom::UB;

  oc.add_equation( PDE_c, {t,z}, {T_NO_LB,   Z_INT},      OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( LDF_q, {t,z}, {T_NO_LB,   FFDom::ALL}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( MOM,   {t,z}, {FFDom::ALL,FFDom::ALL}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( IC_c,  {t,z}, {FFDom::LB, FFDom::ALL}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INITIAL,  0 ) );
  oc.add_equation( IC_q,  {t,z}, {FFDom::LB, FFDom::ALL}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INITIAL,  0 ) );
  oc.add_equation( BC_L,  {t,z}, {T_NO_LB,   FFDom::LB},  OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );
  oc.add_equation( BC_U,  {t,z}, {T_NO_LB,   FFDom::UB},  OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );

  // Scalar outputs (order defines fct-row order); all evaluate to a single value.
  std::vector<OutSpec> ospec;
  oc.add_output( c,    {z,t}, {0.5, 1.0} );  ospec.push_back( { "c_mid=c(0.5,1)", cM_exact(0.5,1.0) } );
  oc.add_output( c,    {z,t}, {1.0, 1.0} );  ospec.push_back( { "c_out=c(1,1)",   cM_exact(1.0,1.0) } );
  oc.add_output( F_in, {z},   {0.0} );       ospec.push_back( { "F_in=int(u c)dt|z=0", F_in_exact } );
  oc.add_output( Q_T,  {t},   {1.0} );       ospec.push_back( { "Q_T=int q dz|t=1",     Q_T_exact } );
  oc.add_output( C_T,  {t},   {1.0} );       ospec.push_back( { "C_T=int c dz|t=1",     C_T_exact } );

  oc.options.REDUCE.ORDER     = red;
  oc.options.CLASSIFY.MODE         = OCFESLV::Options::CLASS_AUTO;
  oc.options.INTERFACE.TYPE   = OCFESLV::Options::IC_AUTO;
  oc.options.INTERFACE.IMPOSITION  = imp;
  // 2026-09-09: PSA_IMP overrides the imposition so the gated solve can be run in the
  // modes this driver does not otherwise test.  MEASURED: OCFE_PDE30 gates IC_STRONG
  // only and fails in IC_TRACE with the claim drop unlocked -- a real defect the corpus
  // could not see.  Env-only; unset leaves the driver's own choice untouched.
  if( char const* v = std::getenv( "PSA_IMP" ) ){
    std::string const m_( v );
    if     ( m_ == "WEAK"   ) oc.options.INTERFACE.IMPOSITION = OCFESLV::Options::IC_WEAK;
    else if( m_ == "TRACE"  ) oc.options.INTERFACE.IMPOSITION = OCFESLV::Options::IC_TRACE;
    else if( m_ == "STRONG" ) oc.options.INTERFACE.IMPOSITION = OCFESLV::Options::IC_STRONG;
  }   // PSA standard
  oc.options.INTERFACE.SAT_SIGMA0       = 10.0;                        // bumped for margin
  oc.options.DISPLAY_LEVEL    = 0;
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
            << " nVar=" << nVar << " nEqn=" << nEqn
            << " nTrace=" << oc.n_colloc_trace()
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

  // State recovery (sanity that the bed itself is right before judging outputs).
  double ce=0., qe=0., ue=0.;
  double const sg[4] = { 0.15, 0.35, 0.65, 0.85 };
  for( double zs : sg ) for( double ts : sg ){
    OCFESLV::t_Coord pt; pt[z]=zs; pt[t]=ts;
    ce = std::max( ce, std::fabs( oc.eval_colloc<double>( c, pt, xv.data(), nullptr, nullptr ) - cM_exact(zs,ts) ) );
    qe = std::max( qe, std::fabs( oc.eval_colloc<double>( q, pt, xv.data(), nullptr, nullptr ) - qM_exact(zs,ts) ) );
    ue = std::max( ue, std::fabs( oc.eval_colloc<double>( u, pt, xv.data(), nullptr, nullptr ) - uM_exact(zs,ts) ) );
  }
  std::cout << std::scientific << std::setprecision(4)
            << "  [C] state err c=" << ce << " q=" << qe << " u=" << ue << "\n";

  // Output functionals at the converged solution.
  std::fill( fct.begin(), fct.end(), 0. );
  // Output functionals from val_functions() (window-summed for evolution integrals) -- correct for
  // both monolithic and marching; a re-eval on xv would see only the last window.
  if( oc.val_functions().size() < nFct ){
    std::cerr << "ERROR: output functionals unavailable\n"; return false;
  }
  fct = oc.val_functions();

  std::cout << "  [D] scalar outputs:\n";
  std::cout << "      " << std::left << std::setw(24) << "functional"
            << std::right << std::setw(16) << "computed"
            << std::setw(16) << "exact"
            << std::setw(13) << "abs err" << "  result\n";
  bool ok = ( ce<1e-7 && qe<1e-7 && ue<1e-7 );
  for( size_t kk=0; kk<ospec.size(); ++kk ){
    size_t const r = oc.row_fct( kk );
    double const val = ( r==std::numeric_limits<size_t>::max() ) ? NAN : fct[ r - nEqn ];
    double const err = std::fabs( val - ospec[kk].exact );
    bool const pass = ( err < 1e-7 );
    ok &= pass;
    std::cout << "      " << std::left << std::setw(24) << ospec[kk].name
              << std::right << std::scientific << std::setprecision(8)
              << std::setw(16) << val << std::setw(16) << ospec[kk].exact
              << std::setprecision(3) << std::setw(13) << err
              << "  " << (pass?"PASS":"FAIL") << "\n";
  }

  std::cout << "  PSA3 " << cname << ( tag ? std::string( " [" ) + tag + "]" : std::string() )
            << ": " << (ok?"PASS":"FAIL") << "\n";
  return ok;
}

int main()
{
  std::cout << "================================================================\n";
  std::cout << "  PSA rung 3: scalar recovery/loading outputs (IC_STRONG, SIGMA0=10)\n";
  std::cout << "================================================================\n";
  bool ok = true;
  ok &= run_psa3( FFDom::CGL, "CGL" );
  ok &= run_psa3( FFDom::LGL, "LGL" );

  // -------------------------------------------------------------------------------------------------------
  // 2026-09-19 -- THE SYSTEMATIC MATRIX: 3 impositions x 2 reductions on the CGL mesh.  PDE26 is the corpus's
  // only model combining a NESTED derivative, an INTEGRAL (the one path _strip_nonlocal_terms_for_symbol acts
  // on) and auto-closure, and it ran IC_STRONG / RED_FULL only.  A cell expected to differ gets a documented
  // XFAIL naming the mechanism, never a loosened bar; an unexpected PASS is flagged.
  // -------------------------------------------------------------------------------------------------------
  std::cout << "\n---- SYSTEMATIC MATRIX (CGL): imposition x reduction ----\n";
  {
    struct MCell { char const* imp; OCFESLV::Options::ImpositionType it;
                   char const* red; OCFESLV::Options::ReductionType rt; };
    MCell const cells[6] = {
      { "WEAK",   OCFESLV::Options::IC_WEAK,   "RED_FULL", OCFESLV::Options::RED_FULL },
      { "TRACE",  OCFESLV::Options::IC_TRACE,  "RED_FULL", OCFESLV::Options::RED_FULL },
      { "STRONG", OCFESLV::Options::IC_STRONG, "RED_FULL", OCFESLV::Options::RED_FULL },
      { "WEAK",   OCFESLV::Options::IC_WEAK,   "RED_MAIN", OCFESLV::Options::RED_MAIN },
      { "TRACE",  OCFESLV::Options::IC_TRACE,  "RED_MAIN", OCFESLV::Options::RED_MAIN },
      { "STRONG", OCFESLV::Options::IC_STRONG, "RED_MAIN", OCFESLV::Options::RED_MAIN } };
    bool matrix_ok = true;
    for( auto const& mc : cells ){
      std::string const tag = std::string( mc.imp ) + "/" + mc.red;
      bool const cell = run_psa3( FFDom::CGL, "CGL", mc.it, mc.rt, tag.c_str() );
      bool const xfail = false;   // none expected: PDE26's balance rows keep a derivative under RED_MAIN
      std::cout << "  MATRIX CELL " << std::left << std::setw(18) << tag
                << ( cell ? ( xfail ? "PASS (unexpected: the documented defect is gone?)" : "PASS" )
                          : ( xfail ? "XFAIL (documented)" : "FAIL" ) ) << "\n";
      if( !cell && !xfail ) matrix_ok = false;
    }
    std::cout << "  MATRIX: " << ( matrix_ok ? "PASS" : "FAIL" ) << "\n";
    ok &= matrix_ok;
  }
  std::cout << "\n  Overall: " << ( ok ? "ALL PASS" : "SOME FAILED" ) << "\n";
  return ok ? 0 : 1;
}
