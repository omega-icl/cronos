// OCFE_STRONGCORNER1.cpp  ---  minimal IC_STRONG tensor-corner gauge probe
// =========================================================================
// PURPOSE (2026-09-07).  Two facts are measured but not yet isolated:
//   (i)  in IC_STRONG, k' == the number of EXPLICIT tau columns (PDE31b 3 and 4,
//        PDE3 12, blk0 0, PDE1 0);
//   (ii) explicit tau columns are protected by the block-eligibility gate, which
//        re-derives the plan (and falls back to legacy) rather than lose them.
// Both were read off models where other structure is mixed in.  This driver
// isolates the effect on the SMALLEST system that can carry it, and bisects on
// the three candidate ingredients, so that "what condition does the explicit
// column carry" can be answered from a system small enough to print whole.
//
//   heat:   dc/dt - D d2c/dz2 = s_c          (parabolic; reduce_order -> Dz_c)
//   EOS:    P - Rg c = 0                     (optional value-slaved state)
//   Darcy:  u + K dP/dz = 0                  (optional derivative-slaved state)
//   manufactured  c = 1 + a(1-z)^2 + b t,  P = Rg c,  u = 2 K Rg a (1-z)
//
// BISECTION AXES
//   LADDER  0 = c alone, 1 = c+P, 2 = c+P+u, 3 = c+P+u WITH convection d(uc)/dz
//           (rung 3 is PDE31b's AXIS B, the only rung it reports as singular;
//            rungs 0-2 measured k'=0 here, so convection is the tipping term)
//   NELT    1 or 2 elements in t.  With NELT=1 there is NO t-seam, hence NO
//           tensor corner, while every other structure is unchanged.  This is
//           the discriminating axis: if k' collapses at NELT=1 the deficiency is
//           corner-borne; if it survives, it is not.
//   NELZ    1 or 2 elements in z (same logic in the other direction).
//
// READ: the [determinacy] line's n / rank / k' for each cell.  Run with
// CRONOS_PLAN_TRACE_DOFS=1 to add "realisation: ... (explicit N)" and check
// N == k' on every cell; with CRONOS_STRONG_NO_EXPLICIT_TAU=1 to watch the
// block-eligibility gate refuse to give the columns up.
// =========================================================================

#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <string>

#define MC__OCFESLV_STRONG_DUMMY_DIAG

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

struct Cell { int ladder; size_t nelt, nelz; int imposition; };

static void run( Cell const& C, char const* impname )
{
  size_t const n_nd = 4;
  std::cout << "\n---- ladder=" << C.ladder << " NELT=" << C.nelt << " NELZ=" << C.nelz
            << " " << impname << " ----\n";

  FFGraph DAG;
  FFVar t=DAG.add_var("t"), z=DAG.add_var("z");
  FFVar c=DAG.add_var("c(t,z)"), P=DAG.add_var("P(t,z)"), u=DAG.add_var("u(t,z)");
  FFPartial OpP;

  bool  const use_conv = ( C.ladder >= 3 );
  FFVar s_c  = use_conv
             ? ( b_c - 2.0*kappa_g*a_c*( 1.0 + 3.0*a_c*(1.0-z)*(1.0-z) + b_c*t ) - D_ax*( 2.0*a_c ) )
             : ( b_c - D_ax*( 2.0*a_c ) );
  FFVar conv = use_conv ? ( u*OpP(c,z) + c*OpP(u,z) ) : FFVar( 0.0 );
  FFVar HEAT = OpP(c,t) + conv - D_ax*OpP(OpP(c,z),z) - s_c;
  FFVar EOS  = P - Rg*c;
  FFVar DARCY= u + Kperm*OpP(P,z);
  FFVar IC_c = c - ( 1.0 + a_c*(1.0-z)*(1.0-z) );
  FFVar BC_Lc= c - ( 1.0 + a_c + b_c*t );
  FFVar BC_Uc= OpP(c,z);

  OCFESLV oc( &DAG );
  oc.add_domain( t, FFDom( 0., 1., C.nelt, FFDom::LGL, n_nd ) );
  oc.add_domain( z, FFDom( 0., 1., C.nelz, FFDom::LGL, n_nd ) );
  oc.add_state( c, {t,z} );
  oc.update_ref( c, [&](OCFESLV::t_Coord const& cr){ return cM(cr.at(z),cr.at(t)); } );
  if( C.ladder >= 1 ){
    oc.add_state( P, {t,z} );
    oc.update_ref( P, [&](OCFESLV::t_Coord const& cr){ return PM(cr.at(z),cr.at(t)); } );
  }
  if( C.ladder >= 2 ){
    oc.add_state( u, {t,z} );
    oc.update_ref( u, [&](OCFESLV::t_Coord const& cr){ return uM(cr.at(z),cr.at(t)); } );
  }
  oc.set_evolution_domain( t );

  int const T_NO_LB = FFDom::ALL - FFDom::LB;
  int const Z_INT   = FFDom::ALL - FFDom::LB - FFDom::UB;
  oc.add_equation( HEAT,  {t,z}, {T_NO_LB,   Z_INT},      OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  if( C.ladder >= 1 )
    oc.add_equation( EOS,   {t,z}, {FFDom::ALL,FFDom::ALL}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  if( C.ladder >= 2 )
    oc.add_equation( DARCY, {t,z}, {FFDom::ALL,FFDom::ALL}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( IC_c,  {t,z}, {FFDom::LB, FFDom::ALL}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INITIAL,  0 ) );
  oc.add_equation( BC_Lc, {t,z}, {T_NO_LB,   FFDom::LB},  OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );
  oc.add_equation( BC_Uc, {t,z}, {T_NO_LB,   FFDom::UB},  OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );

  oc.options.REDUCE.ORDER     = OCFESLV::Options::RED_FULL;
  oc.options.CLASSIFY.MODE         = OCFESLV::Options::CLASS_AUTO;
  oc.options.INTERFACE.TYPE   = OCFESLV::Options::IC_AUTO;
  oc.options.INTERFACE.IMPOSITION  = static_cast<OCFESLV::Options::ImpositionType>( C.imposition );
  oc.options.INTERFACE.SAT_SIGMA0       = 10.0;
  oc.options.DISPLAY_LEVEL    = 1;

  if( !oc.setup() ){ std::cout << "  setup failed\n"; return; }
  size_t const nVar = oc.n_colloc_sta(), nEqn = oc.n_colloc_eqn();
  std::cout << "  nVar=" << nVar << " nEqn=" << nEqn << " nTrace=" << oc.n_colloc_trace()
            << " square=" << (nVar==nEqn?"yes":"NO") << "\n";
  if( nVar != nEqn ) return;

  std::vector<double> varInit, inpInit;
  if( !oc.init( varInit, inpInit, nullptr ) ){ std::cout << "  init failed\n"; return; }
  std::vector<double> xv = varInit;
  for( size_t i=0;i<xv.size();++i ) xv[i] += 0.05*std::sin( 0.7*double(i)+0.2 );
  OCFESLV::SolveReport rep = oc.solve( xv.data(), nullptr, nullptr );

  double cerr=0.;
  double const sg[3]={0.2,0.5,0.8};
  for( double zs : sg ) for( double ts : sg ){
    OCFESLV::t_Coord pt; pt[z]=zs; pt[t]=ts;
    cerr=std::max(cerr,std::fabs(oc.eval_colloc<double>(c,pt,xv.data(),nullptr,nullptr)-cM(zs,ts)));
  }
  std::cout << std::scientific << std::setprecision(3)
            << "  conv=" << (rep.converged?"y":"n") << " it=" << rep.iterations
            << " |r|=" << rep.final_residual << "  c_err=" << cerr
            << "  " << ( rep.converged && cerr<1e-9 ? "PASS" : "fail" ) << "\n";
}

int main( int argc, char** argv )
{
  bool only_strong = false;
  for( int i=1;i<argc;++i ) if( std::string(argv[i]) == "--strong" ) only_strong = true;

  std::cout << "==================================================================\n";
  std::cout << "  Minimal IC_STRONG tensor-corner gauge probe\n";
  std::cout << "  ladder 0=c  1=c+P(value-slave)  2=c+P+u(derivative-slave)\n";
  std::cout << "  NELT=1 removes the t-seam (hence every tensor corner)\n";
  std::cout << "  READ the [determinacy] n / rank / k' of each cell.\n";
  std::cout << "==================================================================\n";

  for( int ladder = 0; ladder <= 3; ++ladder )
    for( size_t nelt : { size_t(1), size_t(2) } )
      for( size_t nelz : { size_t(1), size_t(2) } ){
        run( { ladder, nelt, nelz, OCFESLV::Options::IC_STRONG }, "STRONG" );
        if( !only_strong ){
          run( { ladder, nelt, nelz, OCFESLV::Options::IC_TRACE }, "TRACE" );
          run( { ladder, nelt, nelz, OCFESLV::Options::IC_WEAK  }, "WEAK"  );
        }
      }
  return 0;
}
