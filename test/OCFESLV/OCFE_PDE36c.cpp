// OCFE_PDE36grad_march.cpp  ---  FUNCTION GRADIENTS under SOLVE_MARCHING (FD reference)
// ===========================================================================
// Validates that output-function GRADIENTS with respect to a control are correct
// under marching -- the prerequisite for using marching in optimization.
//
// Per the agreed design, EVERY differentiable control is an INPUT (so dF/dcontrol
// is dF/dinp, already in deriv(); no constant-sensitivity path needed).  Here the
// control is the isotherm affinity b01 of species 1, declared as a LUMPED (0-D)
// scalar input -- a coefficient that genuinely shifts the physical breakthrough
// (unlike the manufactured cases, whose forcing is built to cancel coefficients,
// giving a degenerate zero gradient).
//
// The reduced gradient dF/db01 = direct + (dF/dstate)*(dstate/db01), and the
// second term propagates through the causal marching chain.  Since the marched
// VALUES F(b01) are already correct, a finite-difference gradient is correct by
// construction and serves as the REFERENCE for a future analytic adjoint.  The
// test computes central-difference dEff1/db01 and dInv1/db01 under BOTH modes and
// checks that the MARCHED gradient matches the MONOLITHIC gradient.
// ===========================================================================
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

static double const U_VEL = 1.0, D_ax = 0.1, F_ph = 0.5, k_ldf = 2.0;
static double const qs1 = 1.0, qs2 = 1.0, b02_L = 1.0;   // b01 is now a control input
static double const beta1 = 2.0, beta2 = 1.0, T0_ref = 1.0;
static double const Cp_e = 1.0, G_cv = 1.0, lam = 0.1, dH1 = 1.0, dH2 = 0.5, hw = 0.5, Tw = 1.0;
static double const c0_1 = 0.5, c0_2 = 0.5, tau_in = 0.15, T_end = 5.0;
static double const b01_nom = 3.0;

static inline double c1feed_d( double t ){ return c0_1*( 1.0 - std::exp( -t/tau_in ) ); }
static inline double c2feed_d( double t ){ return c0_2*( 1.0 - std::exp( -t/tau_in ) ); }
static inline double q1star_g( double c1g, double c2g, double b01 )
{ return qs1*b01*c1g/( 1.0 + b01*c1g + b02_L*c2g ); }
static inline double q2star_g( double c1g, double c2g, double b01 )
{ return qs2*b02_L*c2g/( 1.0 + b01*c1g + b02_L*c2g ); }

static int g_pass = 0, g_fail = 0;
static void check_close( char const* name, double got, double want, double tol )
{
  bool ok = std::fabs( got - want ) <= tol;
  ( ok ? g_pass : g_fail )++;
  std::cout << "  " << std::left << std::setw(46) << name
            << std::scientific << std::setprecision(3)
            << " |got-want|=" << std::fabs( got - want ) << " tol=" << tol
            << "  " << ( ok ? "PASS" : "FAIL" ) << "\n";
}
static void check_true( char const* name, bool ok )
{
  ( ok ? g_pass : g_fail )++;
  std::cout << "  " << std::left << std::setw(46) << name << "  " << ( ok ? "PASS" : "FAIL" ) << "\n";
}

// Solve the physical model at b01=@p b01val and return { Eff1, Inv1 }.
// (march: clean-bed ICs auto-detected + transferred; b01 lumped input, no reference,
//  so it is neither re-seeded nor transferred -- it stays the control value.)
static std::array<double,2> run_at_b01( FFDom::TYPE coltype, size_t n_el, size_t n_nd,
                                        bool march, double b01val, double Tic_off, bool& ok,
                                        int agrad_kind = 0, std::array<double,2>* agrad_out = nullptr,
                                        std::array<double,4>* reuse_out = nullptr )
{
  ok = false;
  FFGraph DAG;
  FFVar t  = DAG.add_var( "t" );
  FFVar z  = DAG.add_var( "z" );
  FFVar c1 = DAG.add_var( "c1(t,z)" );
  FFVar c2 = DAG.add_var( "c2(t,z)" );
  FFVar q1 = DAG.add_var( "q1(t,z)" );
  FFVar q2 = DAG.add_var( "q2(t,z)" );
  FFVar T  = DAG.add_var( "T(t,z)"  );
  FFVar c1_ic = DAG.add_var( "c1_ic(z)" );
  FFVar c2_ic = DAG.add_var( "c2_ic(z)" );
  FFVar q1_ic = DAG.add_var( "q1_ic(z)" );
  FFVar q2_ic = DAG.add_var( "q2_ic(z)" );
  FFVar T_ic  = DAG.add_var( "T_ic(z)"  );
  FFVar b01   = DAG.add_var( "b01" );          // isotherm affinity -- LUMPED input CONTROL

  FFPartial  OpP;
  FFIntegral OpI;

  FFVar c1feed = c0_1*( 1.0 - exp( -t/tau_in ) );
  FFVar c2feed = c0_2*( 1.0 - exp( -t/tau_in ) );
  FFVar b1 = b01  *exp( beta1*( 1.0/T - 1.0/T0_ref ) );   // b01 (input) enters the isotherm
  FFVar b2 = b02_L*exp( beta2*( 1.0/T - 1.0/T0_ref ) );
  FFVar den = 1.0 + b1*c1 + b2*c2;
  FFVar q1star = qs1*b1*c1/den;
  FFVar q2star = qs2*b2*c2/den;

  FFVar CONT1 = OpP( c1, t ) + U_VEL*OpP( c1, z ) - D_ax*OpP( OpP( c1, z ), z ) + F_ph*OpP( q1, t );
  FFVar CONT2 = OpP( c2, t ) + U_VEL*OpP( c2, z ) - D_ax*OpP( OpP( c2, z ), z ) + F_ph*OpP( q2, t );
  FFVar LDF1  = OpP( q1, t ) - k_ldf*( q1star - q1 );
  FFVar LDF2  = OpP( q2, t ) - k_ldf*( q2star - q2 );
  FFVar ENE_T = Cp_e*OpP( T, t ) + G_cv*OpP( T, z ) - lam*OpP( OpP( T, z ), z )
              - F_ph*( dH1*OpP( q1, t ) + dH2*OpP( q2, t ) ) + hw*( T - Tw );

  FFVar IC_c1 = march ? ( c1 - c1_ic ) : ( c1 );
  FFVar IC_c2 = march ? ( c2 - c2_ic ) : ( c2 );
  FFVar IC_q1 = march ? ( q1 - q1_ic ) : ( q1 );
  FFVar IC_q2 = march ? ( q2 - q2_ic ) : ( q2 );
  FFVar IC_T  = march ? ( T  - T_ic  ) : ( T - T0_ref );

  FFVar BC_L1 = U_VEL*c1 - D_ax*OpP( c1, z ) - U_VEL*c1feed;
  FFVar BC_L2 = U_VEL*c2 - D_ax*OpP( c2, z ) - U_VEL*c2feed;
  FFVar BC_LT = G_cv*T  - lam*OpP( T,  z ) - G_cv*T0_ref;
  FFVar BC_U1 = OpP( c1, z );
  FFVar BC_U2 = OpP( c2, z );
  FFVar BC_UT = OpP( T,  z );

  FFVar Inv1 = OpI( c1 + F_ph*q1, z );   // bed inventory of sp1 (terminal)
  FFVar Eff1 = OpI( U_VEL*c1, t );        // cumulative effluent of sp1 (accumulated)

  OCFESLV oc( &DAG );
  oc.add_domain( t, FFDom( 0., T_end, n_el, coltype, n_nd ) );
  oc.add_domain( z, FFDom( 0., 1.0,   n_el, coltype, n_nd ) );
  oc.add_state ( c1, {t,z} );  oc.add_state ( c2, {t,z} );
  oc.add_state ( q1, {t,z} );  oc.add_state ( q2, {t,z} );
  oc.add_state ( T,  {t,z} );
  oc.add_input ( b01, {} );              // LUMPED scalar input control (no reference)

  auto c1g = [&]( OCFESLV::t_Coord const& cr ){ return c1feed_d( cr.at(t) )*( 1.0 - 0.5*cr.at(z) ); };
  auto c2g = [&]( OCFESLV::t_Coord const& cr ){ return c2feed_d( cr.at(t) )*( 1.0 - 0.5*cr.at(z) ); };
  oc.update_ref( c1, c1g );
  oc.update_ref( c2, c2g );
  oc.update_ref( q1, [&]( OCFESLV::t_Coord const& cr ){ return q1star_g( c1g(cr), c2g(cr), b01val ); } );
  oc.update_ref( q2, [&]( OCFESLV::t_Coord const& cr ){ return q2star_g( c1g(cr), c2g(cr), b01val ); } );
  oc.update_ref( T,  [&]( OCFESLV::t_Coord const& cr ){
    return T0_ref + F_ph*( dH1*q1star_g(c1g(cr),c2g(cr),b01val) + dH2*q2star_g(c1g(cr),c2g(cr),b01val) )/Cp_e; } );

  if( march ){
    oc.add_input ( c1_ic, {z} );  oc.add_input ( c2_ic, {z} );
    oc.add_input ( q1_ic, {z} );  oc.add_input ( q2_ic, {z} );  oc.add_input ( T_ic, {z} );
    oc.update_ref( c1_ic, [&]( OCFESLV::t_Coord const& cr ){ return 0.0;    } );
    oc.update_ref( c2_ic, [&]( OCFESLV::t_Coord const& cr ){ return 0.0;    } );
    oc.update_ref( q1_ic, [&]( OCFESLV::t_Coord const& cr ){ return 0.0;    } );
    oc.update_ref( q2_ic, [&]( OCFESLV::t_Coord const& cr ){ return 0.0;    } );
    oc.update_ref( T_ic,  [&]( OCFESLV::t_Coord const& cr ){ return T0_ref; } );
  }

  oc.set_evolution_domain( t );

  int const T_NO_LB = FFDom::ALL - FFDom::LB;
  int const Z_INT   = FFDom::ALL - FFDom::LB - FFDom::UB;
  oc.add_equation( CONT1, {t,z}, {T_NO_LB,   Z_INT},      OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( CONT2, {t,z}, {T_NO_LB,   Z_INT},      OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( LDF1,  {t,z}, {T_NO_LB,   FFDom::ALL}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( LDF2,  {t,z}, {T_NO_LB,   FFDom::ALL}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( ENE_T, {t,z}, {T_NO_LB,   Z_INT},      OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( IC_c1, {t,z}, {FFDom::LB, FFDom::ALL}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INITIAL,  0 ) );
  oc.add_equation( IC_c2, {t,z}, {FFDom::LB, FFDom::ALL}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INITIAL,  0 ) );
  oc.add_equation( IC_q1, {t,z}, {FFDom::LB, FFDom::ALL}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INITIAL,  0 ) );
  oc.add_equation( IC_q2, {t,z}, {FFDom::LB, FFDom::ALL}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INITIAL,  0 ) );
  oc.add_equation( IC_T,  {t,z}, {FFDom::LB, FFDom::ALL}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INITIAL,  0 ) );
  oc.add_equation( BC_L1, {t,z}, {T_NO_LB,   FFDom::LB},  OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );
  oc.add_equation( BC_L2, {t,z}, {T_NO_LB,   FFDom::LB},  OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );
  oc.add_equation( BC_LT, {t,z}, {T_NO_LB,   FFDom::LB},  OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );
  oc.add_equation( BC_U1, {t,z}, {T_NO_LB,   FFDom::UB},  OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );
  oc.add_equation( BC_U2, {t,z}, {T_NO_LB,   FFDom::UB},  OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );
  oc.add_equation( BC_UT, {t,z}, {T_NO_LB,   FFDom::UB},  OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );
  oc.add_output( Inv1, {t}, {T_end} );   // fct 0
  oc.add_output( Eff1, {z}, {1.0}   );   // fct 1

  oc.options.REDUCE.ORDER     = OCFESLV::Options::RED_FULL;
  oc.options.CLASSIFY.MODE         = OCFESLV::Options::CLASS_AUTO;
  oc.options.INTERFACE.TYPE   = OCFESLV::Options::IC_AUTO;
  oc.options.INTERFACE.IMPOSITION  = OCFESLV::Options::IC_STRONG;
  oc.options.INTERFACE.SAT_SIGMA0       = 10.0;
  oc.options.DISPLAY_LEVEL    = 0;
#if defined(CRONOS__WITH_SPQR)
  oc.options.SOLVE.FACTORIZATION = OCFESLV::Options::SOLVE_SPQR;
#else
  oc.options.SOLVE.FACTORIZATION = OCFESLV::Options::SOLVE_SUPERLU;
#endif
  if( march ) oc.options.SOLVE.WARMSTART = OCFESLV::Options::BROADCAST_IC;
  else        oc.options.SOLVE.MARCHING  = false;

  if( !oc.setup() ){ std::cerr << "  setup() FAILED\n"; return {0.,0.}; }
  std::vector<double> xv, inp;
  if( !oc.init( xv, inp, nullptr ) ){ std::cerr << "  init() FAILED\n"; return {0.,0.}; }
  oc.set_input_values( b01, { b01val }, inp.data() );          // coefficient control value
  if( march ){
    size_t const ntic = oc.get_input_values( T_ic, inp.data() ).size();
    oc.set_input_values( T_ic, std::vector<double>( ntic, T0_ref + Tic_off ), inp.data() );  // IC control value
  }
  double const* ip = inp.empty() ? nullptr : inp.data();

  if( agrad_kind == 4 && march && reuse_out ){
    // reduced-space CONTROL REGISTRY: register controls, encode/decode round-trip, and read the
    // reduced Jacobian from the public solve_fsens() -> sens_jacobian() (nf x ncd, row-major).
    oc.register_control( b01 );          // coefficient (lumped)
    oc.register_control( T_ic );         // IC-input (distributed)
    std::vector<double> p0; oc.encode_controls( inp.data(), p0 );
    std::vector<double> inp2( inp ); oc.decode_controls( p0, inp2.data() );
    std::vector<double> p1; oc.encode_controls( inp2.data(), p1 );
    double rt = 0.; for( size_t i = 0; i < p0.size() && i < p1.size(); ++i ) rt = std::max( rt, std::fabs( p0[i]-p1[i] ) );
    std::vector<double> xvS( xv ), inpS( inp );
    if( !oc.solve_fsens( xvS.data(), inpS.data(), nullptr ) ){ std::cerr << "  fsens FAILED\n"; return {0.,0.}; }
    size_t const ncd = oc.n_control_dof(), ncf = oc.n_colloc_fct();
    auto const& C = oc.controls();
    std::vector<double> const& J = oc.sens_jacobian();
    size_t const ib01 = C.at(b01).offset, itic0 = C.at(T_ic).offset, iticN = C.at(T_ic).ndof;
    double dEff1_db01 = ( ncf > 1 ) ? J[ 1*ncd + ib01 ] : 0.;                 // Eff1 = fct row 1
    double dEff1_dtic = 0.; if( ncf > 1 ) for( size_t j = 0; j < iticN; ++j ) dEff1_dtic += J[ 1*ncd + itic0 + j ];
    (*reuse_out)[0] = dEff1_db01;    // dEff1/db01 via registry
    (*reuse_out)[1] = dEff1_dtic;    // dEff1/dTic (uniform shift) via registry
    (*reuse_out)[2] = rt;            // encode/decode round-trip error
    (*reuse_out)[3] = double( ncd ); // control-vector dimension
    ok = true;
    return {0.,0.};
  }
  if( agrad_kind == 3 && march && reuse_out ){
    // dEff1/dInv1 w.r.t. b01 and a uniform T_ic shift via the public reduced Jacobian (solve_fsens):
    // lumped b01 -> its single control column; distributed T_ic -> the all-ones (uniform-shift) sum.
    oc.register_control( b01 );          // coefficient (lumped)
    oc.register_control( T_ic );         // IC-input (distributed)
    std::vector<double> xvS( xv ), inpS( inp );
    if( !oc.solve_fsens( xvS.data(), inpS.data(), nullptr ) ){ std::cerr << "  fsens FAILED\n"; return {0.,0.}; }
    size_t const ncd = oc.n_control_dof(), ncf = oc.n_colloc_fct();
    auto const& C = oc.controls();
    std::vector<double> const& J = oc.sens_jacobian();
    if( ncf < 2 || J.size() < 2*ncd ){ std::cerr << "  apply FAILED\n"; return {0.,0.}; }
    size_t const ib01 = C.at(b01).offset, itic0 = C.at(T_ic).offset, iticN = C.at(T_ic).ndof;
    auto colB = [&]( size_t r ){ return J[ r*ncd + ib01 ]; };
    auto colT = [&]( size_t r ){ double s = 0.; for( size_t j = 0; j < iticN; ++j ) s += J[ r*ncd + itic0 + j ]; return s; };
    (*reuse_out)[0] = colB(1);   // dEff1/db01
    (*reuse_out)[1] = colB(0);   // dInv1/db01
    (*reuse_out)[2] = colT(1);   // dEff1/dTic
    (*reuse_out)[3] = colT(0);   // dInv1/dTic
    ok = true;
    return {0.,0.};
  }
  if( agrad_kind && march ){
    // analytic dF/d(control) via the public reduced Jacobian (solve_fsens); lumped b01 -> its single
    // control column, distributed T_ic -> the all-ones (uniform-shift) column sum.
    oc.register_control( b01 );          // coefficient (lumped)
    oc.register_control( T_ic );         // IC-input (distributed)
    std::vector<double> xvS( xv ), inpS( inp );
    if( !oc.solve_fsens( xvS.data(), inpS.data(), nullptr ) ){ std::cerr << "  solve_fsens FAILED\n"; return {0.,0.}; }
    size_t const ncd = oc.n_control_dof(), ncf = oc.n_colloc_fct();
    auto const& C = oc.controls();
    std::vector<double> const& J = oc.sens_jacobian();
    if( ncf < 2 || J.size() < 2*ncd ){ std::cerr << "  solve_fsens FAILED\n"; return {0.,0.}; }
    std::vector<double> dFdp( ncf, 0. );
    if( agrad_kind == 1 ){                          // b01: lumped coefficient
      size_t const ib01 = C.at(b01).offset;
      for( size_t r = 0; r < ncf; ++r ) dFdp[r] = J[ r*ncd + ib01 ];
    }
    else {                                          // T_ic: distributed IC input, uniform-shift (all-ones)
      size_t const itic0 = C.at(T_ic).offset, iticN = C.at(T_ic).ndof;
      for( size_t r = 0; r < ncf; ++r ){ double s = 0.; for( size_t j = 0; j < iticN; ++j ) s += J[ r*ncd + itic0 + j ]; dFdp[r] = s; }
    }
    (*agrad_out)[0] = dFdp[0];   // dInv1/d(control)  (fct row 0)
    (*agrad_out)[1] = dFdp[1];   // dEff1/d(control)  (fct row 1)
    ok = true;
    return {0.,0.};
  }

  OCFESLV::SolveReport const rep = oc.solve( xv.data(), inp.data(), nullptr );
  if( !rep.converged ){ std::cerr << "  solve FAILED at b01=" << b01val << "\n"; return {0.,0.}; }

  std::vector<double> fct;
  if( march ) fct = oc.val_functions();
  else{
    fct.assign( oc.n_colloc_fct(), 0. );
    std::vector<double> eqn( oc.n_colloc_eqn(), 0. );
    oc.eval( eqn.data(), fct.data(), xv.data(), ip, nullptr );
  }
  if( fct.size() < 2 ) return {0.,0.};
  ok = true;
  return { fct[1], fct[0] };   // { Eff1, Inv1 }
}

int main()
{
  std::cout << "================================================================\n"
            << "  FUNCTION GRADIENTS under SOLVE_MARCHING -- FD reference\n"
            << "  control: isotherm affinity b01 (LUMPED input); outputs Eff1, Inv1\n"
            << "  central-difference dF/db01, marched vs monolithic\n"
            << "================================================================\n";
  size_t const n_el = 5, n_nd = 6;
  double const eps  = 1e-4;

  auto grad = [&]( bool march, double& dEff, double& dInv ) -> bool {
    bool o1=false, o2=false;
    std::array<double,2> fp = run_at_b01( FFDom::LGL, n_el, n_nd, march, b01_nom+eps, 0.0, o1 );
    std::array<double,2> fm = run_at_b01( FFDom::LGL, n_el, n_nd, march, b01_nom-eps, 0.0, o2 );
    if( !( o1 && o2 ) ) return false;
    dEff = ( fp[0] - fm[0] )/( 2.0*eps );
    dInv = ( fp[1] - fm[1] )/( 2.0*eps );
    return true;
  };

  double dEff_mono=0., dInv_mono=0., dEff_march=0., dInv_march=0.;
  bool const okM = grad( /*march=*/false, dEff_mono,  dInv_mono  );
  bool const okm = grad( /*march=*/true,  dEff_march, dInv_march );

  if( okM && okm ){
    std::cout << std::scientific << std::setprecision(4)
              << "\n  dEff1/db01:  monolithic=" << dEff_mono  << "  marched=" << dEff_march << "\n"
              <<   "  dInv1/db01:  monolithic=" << dInv_mono  << "  marched=" << dInv_march << "\n";
    // marched gradient must match the monolithic gradient (both FD, same model) to the
    // discretization/FD level; a scale-relative tolerance keeps it meaningful.
    double const teff = 1e-3*std::max( 1e-6, std::fabs(dEff_mono) ) + 1e-6;
    double const tinv = 1e-3*std::max( 1e-6, std::fabs(dInv_mono) ) + 1e-6;
    check_close( "dEff1/db01 marched ~ monolithic", dEff_march, dEff_mono, teff );
    check_close( "dInv1/db01 marched ~ monolithic", dInv_march, dInv_mono, tinv );
    std::cout << "  (gradients are non-trivial: |dEff/db01|~" << std::fabs(dEff_mono)
              << ", so this exercises a real coefficient sensitivity)\n";

    // --- analytic forward-sensitivity march vs the FD reference ---
    std::array<double,2> agrad = {0.,0.}; bool oka=false;
    run_at_b01( FFDom::LGL, n_el, n_nd, /*march=*/true, b01_nom, 0.0, oka, 1, &agrad );
    if( oka ){
      double const dInv_an = agrad[0], dEff_an = agrad[1];
      std::cout << std::scientific << std::setprecision(4)
                << "\n  analytic (forward-sensitivity march):\n"
                << "  dEff1/db01:  analytic=" << dEff_an << "  FD(march)=" << dEff_march << "\n"
                << "  dInv1/db01:  analytic=" << dInv_an << "  FD(march)=" << dInv_march << "\n";
      check_close( "dEff1/db01 analytic ~ FD", dEff_an, dEff_march, teff );
      check_close( "dInv1/db01 analytic ~ FD", dInv_an, dInv_march, tinv );
    }
    else{ ++g_fail; std::cout << "  analytic gradient run FAILED\n"; }

    // ===== IC-input control (uniform initial-T shift) : marched FD vs analytic =====
    double const epsT = 1e-4;
    bool ot1=false, ot2=false;
    std::array<double,2> tp = run_at_b01( FFDom::LGL, n_el, n_nd, true, b01_nom, +epsT, ot1 );
    std::array<double,2> tm = run_at_b01( FFDom::LGL, n_el, n_nd, true, b01_nom, -epsT, ot2 );
    std::array<double,2> agT = {0.,0.}; bool okT=false;
    run_at_b01( FFDom::LGL, n_el, n_nd, true, b01_nom, 0.0, okT, 2, &agT );
    if( ot1 && ot2 && okT ){
      double const dEffT_fd = ( tp[0]-tm[0] )/( 2.0*epsT ), dInvT_fd = ( tp[1]-tm[1] )/( 2.0*epsT );
      double const dInvT_an = agT[0], dEffT_an = agT[1];
      std::cout << std::scientific << std::setprecision(4)
                << "\n  IC-input control (uniform T(0,z) shift):\n"
                << "  dEff1/dTic:  analytic=" << dEffT_an << "  FD(march)=" << dEffT_fd << "\n"
                << "  dInv1/dTic:  analytic=" << dInvT_an << "  FD(march)=" << dInvT_fd << "\n";
      double const tET = 1e-3*std::max( 1e-6, std::fabs(dEffT_fd) ) + 1e-6;
      double const tIT = 1e-3*std::max( 1e-6, std::fabs(dInvT_fd) ) + 1e-6;
      check_close( "dEff1/dTic analytic ~ FD (IC control)", dEffT_an, dEffT_fd, tET );
      check_close( "dInv1/dTic analytic ~ FD (IC control)", dInvT_an, dInvT_fd, tIT );
      std::cout << "  (confirms IC-input control + a second control type via the same forward march)\n";

      // ===== reuse: one primal setup, two sensitivity applies (b01 + T_ic) =====
      std::array<double,4> ru = {0.,0.,0.,0.}; bool okr=false;
      run_at_b01( FFDom::LGL, n_el, n_nd, true, b01_nom, 0.0, okr, 3, nullptr, &ru );
      if( okr ){
        std::cout << std::scientific << std::setprecision(4)
                  << "\n  reuse (1 primal setup + 2 sensitivity applies):\n"
                  << "  dEff1/db01=" << ru[0] << "  dEff1/dTic=" << ru[2] << "\n";
        check_close( "reuse dEff1/db01 == FD", ru[0], dEff_march, teff );
        check_close( "reuse dEff1/dTic == FD", ru[2], dEffT_fd,   tET );
        check_close( "reuse dInv1/db01 == FD", ru[1], dInv_march, tinv );
        check_close( "reuse dInv1/dTic == FD", ru[3], dInvT_fd,   tIT );
        std::cout << "  (proves 1 primal march + N cheap sensitivity marches: shared factorizations)\n";
      }
      else{ ++g_fail; std::cout << "  reuse run FAILED\n"; }

      // ===== reduced-space control registry: encode/decode + sensitivity via decode =====
      std::array<double,4> rg = {0.,0.,0.,0.}; bool okg=false;
      run_at_b01( FFDom::LGL, n_el, n_nd, true, b01_nom, 0.0, okg, 4, nullptr, &rg );
      if( okg ){
        std::cout << std::scientific << std::setprecision(4)
                  << "\n  control registry (reduced input space, dim=" << int(rg[3]) << "):\n"
                  << "  encode/decode round-trip err=" << rg[2] << "\n"
                  << "  dEff1/db01(via decode)=" << rg[0] << "  dEff1/dTic(via decode)=" << rg[1] << "\n";
        check_close( "registry encode/decode round-trip", rg[2], 0., 1e-12 );
        check_close( "registry dEff1/db01 (via decode) == FD", rg[0], dEff_march, teff );
        check_close( "registry dEff1/dTic (via decode) == FD", rg[1], dEffT_fd,   tET );
        check_true ( "registry n_control_dof > 1 (lumped + distributed)", rg[3] > 1.5 );
        std::cout << "  (encode/decode + fsens compose into the reduced input->output map)\n";
      }
      else{ ++g_fail; std::cout << "  control-registry run FAILED\n"; }
    }
    else{ ++g_fail; std::cout << "  IC-control gradient run FAILED\n"; }
  }
  else{ ++g_fail; std::cout << "  FD gradient runs FAILED\n"; }

  std::cout << "\n============================================================\n"
            << "  RESULT: " << g_pass << " passed, " << g_fail << " failed  -- "
            << ( g_fail == 0 ? "ALL PASS" : "SOME FAILED" ) << "\n"
            << "============================================================\n";
  return g_fail == 0 ? 0 : 1;
}
