// OCFE_PDE20h.cpp  ---  a declared initial condition IMPLIED by a high-index reduction's hidden constraints (8 Nov 2026)
// ===========================================================================
// The high-index reduction materialises its hidden constraint levels as INITIAL rows at the initial point.  A declared
// INITIAL row those levels already imply over-determined the initial point (OCFE_PDE20b: 13 rows for 12 unknowns,
// refused by the DOF audit).  FFModel now SETS IT ASIDE (a maximum matching of the INITIAL rows to the differential
// states, the materialised levels first) and the solver checks it against the solution
// (OCFESLV::initial_consistency_defect).  Corpus M2, x' = y, 0 = x - a(t), a = 1 + 0.5 t - 0.3 t^2 (exact x = a, y = a'):
//   (A) the declared x(0) = a(0), CONSISTENT: set aside, the solution exact, the defect 0;
//   (B) a declared x(0) = a(0) + 0.5, INCONSISTENT: set aside, the hidden constraint holds (x = a), the defect 0.5;
//   (C) plus a free state z' = -z, z(0) = 2: x(0) set aside, z(0) KEPT (it is free data), z(TF) = 2 exp(-TF).
#include <iostream>
#include <iomanip>
#include <cmath>
#include "ocfeslv.hpp"
using namespace mc;

static double const TF = 0.5;
static size_t const NEL_T = 3, NT = 6;
static int nfail = 0;
static void check( std::string const& w, bool ok ){ std::cout << "  " << std::left << std::setw(92) << w << ( ok? "PASS": "FAIL" ) << "\n"; if( !ok ) ++nfail; }
static double a_ex( double t ){ return 1.0 + 0.5*t - 0.3*t*t; }
static double ap_ex( double t ){ return 0.5 - 0.6*t; }

struct Result { bool ok = false; size_t set_aside = 0, nchecked = 0, declared = 0, free = 0; double defect = -1., ex = 1., ey = 1., z = 0.; };

static Result run( double offset, bool with_z )
{
  Result R;
  FFGraph DAG;  FFPartial OpP;
  FFVar t = DAG.add_var( "t" ), x = DAG.add_var( "x(t)" ), y = DAG.add_var( "y(t)" ), z = DAG.add_var( "z(t)" );
  FFVar Ae = 1.0 + 0.5*t - 0.3*t*t;
  OCFESLV oc( &DAG );
  oc.add_domain( t, FFDom( 0., TF, NEL_T, FFDom::LGR, NT ) );
  oc.add_state( x, {t} );  oc.add_state( y, {t} );
  oc.update_ref( x, []( OCFESLV::t_Coord const& ){ return 1.0; } );  oc.update_ref( y, []( OCFESLV::t_Coord const& ){ return 0.0; } );
  OCFESLV::EqnOptions io( OCFESLV::EqnRole::INTERIOR, 0 ), ii( OCFESLV::EqnRole::INITIAL, 0 );
  int const T_INT = FFDom::ALL - FFDom::LB;
  oc.add_equation( OpP( x, t ) - y, {t}, {T_INT},      io );
  oc.add_equation( x - Ae,          {t}, {FFDom::ALL}, io );
  oc.add_equation( x - Ae - offset, {t}, {FFDom::LB},  ii );         // the declared x(0)
  if( with_z ){
    oc.add_state( z, {t} );  oc.update_ref( z, []( OCFESLV::t_Coord const& ){ return 1.0; } );
    oc.add_equation( OpP( z, t ) + z, {t}, {T_INT},     io );
    oc.add_equation( z - 2.,          {t}, {FFDom::LB}, ii );         // free data: kept
  }
  oc.set_evolution_domain( t );
  oc.options.REDUCE.ORDER = OCFESLV::Options::RED_MAIN;
  oc.options.INTERFACE.IMPOSITION = OCFESLV::Options::IC_WEAK;
  oc.options.DISPLAY_LEVEL = 1;  oc.options.SOLVE.MAX_ITER = 50;  oc.options.SOLVE.RES_TOL = 1e-11;
#if defined(CRONOS__WITH_SPQR)
  oc.options.SOLVE.FACTORIZATION = OCFESLV::Options::SOLVE_SPQR;
#endif
  try{ if( !oc.setup() ) return R; } catch( ... ){ return R; }
  R.set_aside = oc.initial_set_aside().size();  R.declared = oc.initial_data().declared;  R.free = oc.initial_data().free;
  if( oc.solve( std::vector<FFModel::InputVal>{} ) != OCFESLV::NORMAL ) return R;
  R.ok = true;
  R.defect = oc.initial_consistency_defect( &R.nchecked );
  R.ex = R.ey = 0.;
  for( int k = 0; k <= 20; ++k ){
    double const tk = TF*k/20.;
    R.ex = std::max( R.ex, std::fabs( oc.eval_solution( x, {{t,tk}} ) - a_ex( tk ) ) );
    R.ey = std::max( R.ey, std::fabs( oc.eval_solution( y, {{t,tk}} ) - ap_ex( tk ) ) );
  }
  if( with_z ) R.z = oc.eval_solution( z, {{t,TF}} );
  return R;
}

int main()
{
  std::cout << std::scientific << std::setprecision(3);
  std::cout << "(A) x(0) = a(0) declared, consistent with the materialised hidden constraint\n";
  Result A = run( 0., false );
  std::cout << "    set aside " << A.set_aside << ", declared " << A.declared << ", free " << A.free << ";  max|x-a| " << A.ex << "  max|y-a'| " << A.ey << "  defect " << A.defect << " (" << A.nchecked << " checked)\n";
  check( "(A) sets up and solves: the declared x(0) set aside, the model balanced (declared = free = 0)", A.ok && A.set_aside == 1 && A.declared == 0 && A.free == 0 );
  check( "(A) x = a(t), y = a'(t) to 1e-10, y included at t = 0", A.ok && A.ex < 1e-10 && A.ey < 1e-10 );
  check( "(A) the set-aside row is consistent with the solution: defect < 1e-12", A.ok && A.nchecked == 1 && A.defect < 1e-12 );

  std::cout << "(B) x(0) = a(0) + 0.5 declared, INCONSISTENT with the hidden constraint\n";
  Result B = run( 0.5, false );
  std::cout << "    set aside " << B.set_aside << ";  max|x-a| " << B.ex << "  max|y-a'| " << B.ey << "  defect " << B.defect << "\n";
  check( "(B) set aside as in (A): the hidden constraint holds, x = a(t) and y = a'(t) to 1e-10", B.ok && B.set_aside == 1 && B.ex < 1e-10 && B.ey < 1e-10 );
  check( "(B) the inconsistency is reported, not hidden: defect = 0.5 to 1e-10", B.ok && B.nchecked == 1 && std::fabs( B.defect - 0.5 ) < 1e-10 );

  std::cout << "(C) plus a free state z' = -z, z(0) = 2\n";
  Result C = run( 0., true );
  double const zex = 2.*std::exp( -TF );
  std::cout << "    set aside " << C.set_aside << ", declared " << C.declared << ", free " << C.free << ";  z(TF) " << C.z << " (exact " << zex << ")  max|x-a| " << C.ex << "\n";
  check( "(C) only the implied row set aside: x(0) aside, z(0) kept (declared = free = 1)", C.ok && C.set_aside == 1 && C.declared == 1 && C.free == 1 );
  check( "(C) z(TF) = 2 exp(-TF) to 1e-8, x = a(t) to 1e-10", C.ok && std::fabs( C.z - zex ) < 1e-8 && C.ex < 1e-10 );

  std::cout << "==========  OCFE_PDE20h: " << ( nfail? std::to_string( nfail ) + " FAILED": "ALL PASS" ) << "  ==========\n";
  return nfail? 1: 0;
}
