// OCFE_autoelim.cpp  ---  MMS de-risking probe for AUTO_DIFF_ELIM (1- and 2-species)
// ============================================================================
// Reproduces MBC's ELIMINABLE trigger in the smallest steady systems, OUTSIDE the
// 30-equation MBC model, with manufactured exact solutions so the flag-ON run
// asserts SOLUTION-EXACTNESS (not just that the block squared).
//
// The MBC pattern (E_TL references d/dz Cg; GAS_C -- an algebraic advection in
// vAlgEqn -- determines d/dz Cg via the wall flux d/drd Cd) reduced to essentials:
//
//   states:  T(z)      -- differential ODE state          (mirrors Tl)
//            c_i(z)     -- advected concentration(s)       (mirror Cg [, Vg])
//            d(r,z)     -- 2-D membrane state              (mirrors Cd; cross-dir d/dr d)
//            w(z)       -- VALUE-SLAVED algebraic state    (pushes SOURCE(s) -> vAlgEqn)
//
//   E   (consumer, differential):  d/dz T + SUM_i d/dz c_i - fE       = 0
//   G_i (source  -> vAlgEqn):       d/dz c_i + d/dr d - w - fG_i       = 0   @ r=LB
//   W   (value-slaved def):         w - (1 + z)                        = 0   (deriv-free -> vAlgEqn)
//   D   (2-D diffusion):            d2/dr2 d - 2(1+z)                   = 0
//   + Dirichlet BCs on T(0), c_i(0), d(r=0), d(r=1).
//
// NSPECIES=1 : one derivative into the consumer (single elimination).
// NSPECIES=2 : TWO derivatives (d/dz c1, d/dz c2) into the SAME consumer E -- the
//              exact MBC multiplicity (Cg & Vg into E_TL).  Exercises the fixpoint:
//              eliminate c1, re-classify, re-FIND the (rewritten) consumer, eliminate c2.
//
// Manufactured exact:  c1=1+z, c2=1+2z, T=2+3z, d=(1+z)(r+r^2), w=1+z.
//
// VALIDATION (DM diagnostic is the stage-1 oracle):
//   * flag OFF: expect [DM-global] ... c_i[... ELIMINABLE], block RECTANGULAR (defect=NSPECIES).
//   * flag ON : expect ** auto_diff_elim: eliminated d/dz c_i ... (NSPECIES times),
//               block SQUARE, MMS state errors ~1e-12.
// ============================================================================

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

static inline double c1M( double z ){ return 1.0 + z; }
static inline double c2M( double z ){ return 1.0 + 2.0*z; }
static inline double TM ( double z ){ return 2.0 + 3.0*z; }
static inline double dM ( double r, double z ){ return (1.0+z)*(r + r*r); }
static inline double wM ( double z ){ return 1.0 + z; }

static void run( char const* tag, bool auto_elim, int nspecies )
{
  std::cout << "\n---- " << tag << " ----\n";
  size_t const n_el = 2, n_nd = 4;

  FFGraph DAG;
  FFVar z = DAG.add_var("z"), r = DAG.add_var("r");
  FFVar T = DAG.add_var("T(z)"), c1 = DAG.add_var("c1(z)"), c2 = DAG.add_var("c2(z)");
  FFVar d = DAG.add_var("d(r,z)"), w = DAG.add_var("w(z)");
  FFPartial OpP;

  double const fE = ( nspecies == 2 ) ? 6.0 : 4.0;      // dz T + dz c1 [+ dz c2] = 3+1[+2]
  FFVar E = ( nspecies == 2 )
          ? ( OpP(T,z) + OpP(c1,z) + OpP(c2,z) - fE )
          : ( OpP(T,z) + OpP(c1,z)             - fE );

  FFVar G1 = OpP(c1,z) + OpP(d,r) - w - 1.0;            // dz c1*=1 -> fG1=1
  FFVar G2 = OpP(c2,z) + OpP(d,r) - w - 2.0;            // dz c2*=2 -> fG2=2
  FFVar W  = w - ( 1.0 + z );                           // value-slaved
  FFVar D  = OpP(OpP(d,r),r) - 2.0*( 1.0 + z );         // 2-D diffusion

  FFVar T_BC  = T  - 2.0;
  FFVar c1_BC = c1 - 1.0;
  FFVar c2_BC = c2 - 1.0;
  FFVar d_BC0 = d  - 0.0;
  FFVar d_BC1 = d  - 2.0*( 1.0 + z );

  OCFESLV oc( &DAG );
  oc.add_domain( z, FFDom( 0., 1., n_el, FFDom::LGL, n_nd ) );
  oc.add_domain( r, FFDom( 0., 1., n_el, FFDom::LGL, n_nd ) );
  oc.add_state( T,  {z} );
  oc.add_state( c1, {z} );
  if( nspecies == 2 ) oc.add_state( c2, {z} );
  oc.add_state( d,  {r,z} );
  oc.add_state( w,  {z} );
  oc.update_ref( T,  [&](OCFESLV::t_Coord const& cr){ return TM (cr.at(z)); } );
  oc.update_ref( c1, [&](OCFESLV::t_Coord const& cr){ return c1M(cr.at(z)); } );
  if( nspecies == 2 )
    oc.update_ref( c2, [&](OCFESLV::t_Coord const& cr){ return c2M(cr.at(z)); } );
  oc.update_ref( d,  [&](OCFESLV::t_Coord const& cr){ return dM(cr.at(r),cr.at(z)); } );
  oc.update_ref( w,  [&](OCFESLV::t_Coord const& cr){ return wM (cr.at(z)); } );

  int const Z_NO_LB = FFDom::ALL - FFDom::LB;
  int const R_INT   = FFDom::ALL - FFDom::LB - FFDom::UB;

  oc.add_equation( E,     {z},   {Z_NO_LB},              OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( T_BC,  {z},   {FFDom::LB},            OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );
  oc.add_equation( G1,    {z,r}, {Z_NO_LB, FFDom::LB},   OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( c1_BC, {z},   {FFDom::LB},            OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );
  if( nspecies == 2 ){
    oc.add_equation( G2,    {z,r}, {Z_NO_LB, FFDom::LB}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
    oc.add_equation( c2_BC, {z},   {FFDom::LB},          OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );
  }
  oc.add_equation( W,     {z},   {FFDom::ALL},           OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( D,     {r,z}, {R_INT, FFDom::ALL},    OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( d_BC0, {r,z}, {FFDom::LB, FFDom::ALL},OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );
  oc.add_equation( d_BC1, {r,z}, {FFDom::UB, FFDom::ALL},OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );

  oc.options.REDUCE.ORDER     = OCFESLV::Options::RED_FULL;
  oc.options.CLASSIFY.MODE         = OCFESLV::Options::CLASS_AUTO;
  oc.options.INTERFACE.TYPE   = OCFESLV::Options::IC_AUTO;
  oc.options.INTERFACE.IMPOSITION  = OCFESLV::Options::IC_STRONG;
  oc.options.AUTO.DIFF_ELIM   = auto_elim;
  oc.options.DISPLAY_LEVEL    = 1;
  oc.options.SOLVE.VERBOSE    = true;
  
  if( !oc.setup() ){ std::cerr << "  setup failed\n"; return; }
  size_t const nVar = oc.n_colloc_sta(), nEqn = oc.n_colloc_eqn();
  std::cout << "  nVar=" << nVar << " nEqn=" << nEqn
            << " square=" << (nVar==nEqn?"yes":"NO") << "\n";
  if( nVar != nEqn ){ std::cout << "  (rectangular -- flag-OFF baseline; skip solve)\n"; return; }

  std::vector<double> varInit, inpInit;
  if( !oc.init( varInit, inpInit, nullptr ) ){ std::cerr << "  init failed\n"; return; }
  std::vector<double> res( nEqn, 0. );
  double Aerr = 0.;
  if( oc.eval( res.data(), nullptr, varInit.data(), nullptr, nullptr ) )
    for( double v : res ) Aerr = std::max( Aerr, std::fabs(v) );

  std::vector<double> xv = varInit;
  for( size_t i=0;i<xv.size();++i ) xv[i] += 0.05*std::sin( 0.7*double(i)+0.2 );
  OCFESLV::SolveReport rep = oc.solve( xv.data(), nullptr, nullptr );

  double Terr=0., c1err=0., c2err=0., derr=0., werr=0.;
  double const sg[3] = { 0.2, 0.5, 0.8 };
  for( double zs : sg ){
    OCFESLV::t_Coord pz; pz[z]=zs;
    Terr  = std::max( Terr,  std::fabs(oc.eval_colloc<double>(T ,pz,xv.data(),nullptr,nullptr)-TM (zs)) );
    c1err = std::max( c1err, std::fabs(oc.eval_colloc<double>(c1,pz,xv.data(),nullptr,nullptr)-c1M(zs)) );
    werr  = std::max( werr,  std::fabs(oc.eval_colloc<double>(w ,pz,xv.data(),nullptr,nullptr)-wM (zs)) );
    if( nspecies == 2 )
      c2err = std::max( c2err, std::fabs(oc.eval_colloc<double>(c2,pz,xv.data(),nullptr,nullptr)-c2M(zs)) );
    for( double rs : sg ){
      OCFESLV::t_Coord pr; pr[z]=zs; pr[r]=rs;
      derr = std::max( derr, std::fabs(oc.eval_colloc<double>(d,pr,xv.data(),nullptr,nullptr)-dM(rs,zs)) );
    }
  }
  bool const ok = rep.converged && Terr<1e-9 && c1err<1e-9 && c2err<1e-9 && derr<1e-9 && werr<1e-9;
  std::cout << std::scientific << std::setprecision(3)
            << "  [A]=" << Aerr << "  conv=" << (rep.converged?"y":"n")
            << " it=" << rep.iterations << " |r|=" << rep.final_residual
            << "  [C] T=" << Terr << " c1=" << c1err << " c2=" << c2err
            << " d=" << derr << " w=" << werr
            << "  " << ( ok ? "PASS" : "fail" ) << "\n";
}

int main()
{
  std::cout << "================================================================\n";
  std::cout << "  AUTO_DIFF_ELIM de-risking probe (MBC ELIMINABLE pattern, MMS)\n";
  std::cout << "================================================================\n";

  std::cout << "\n=== 1-SPECIES (single elimination) ===\n";
  run( "1sp flag OFF (expect [DM-global] c1 ELIMINABLE, rectangular)", false, 1 );
  run( "1sp flag ON  (expect elimination -> square -> MMS-exact)",     true,  1 );

  std::cout << "\n=== 2-SPECIES (two eliminations into one consumer + fixpoint = MBC multiplicity) ===\n";
  run( "2sp flag OFF (expect c1 & c2 ELIMINABLE, defect=2)",           false, 2 );
  run( "2sp flag ON  (expect 2 eliminations -> square -> MMS-exact)",  true,  2 );

  std::cout << "\n  Read: 2sp flag-ON must print TWO 'eliminated d/dz c*' lines and square the block.\n";
  return 0;
}
