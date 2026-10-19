// DAESLV_icseed -- the initial-condition system at its SEED: rank, and where the seed comes from.
// =====================================================================================================
// ic_probe() drives the very callbacks KINSOL drives, so the matrix below is the one the solver sees.
// Gaussian elimination with full pivoting then gives the rank and, when it is deficient, the row and
// column that run out -- which is the thing worth knowing: WHICH row is dependent on the others.
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include "daeslv_idas.hpp"

using namespace mc;

int main()
{
  FFPartial OpP;
  int const TI = FFDom::ALL - FFDom::LB;
  auto io = FFModel::EqnOptions( FFModel::EqnRole::INTERIOR, 0 );
  auto ii = FFModel::EqnOptions( FFModel::EqnRole::INITIAL,  0 );

  FFGraph G; DAESLV_IDAS M; M.set_dag( &G );
  FFVar t = G.add_var("t"), x = G.add_var("x(t)"), y = G.add_var("y(t)"),
        u = G.add_var("u(t)"), v = G.add_var("v(t)"), lam = G.add_var("lam(t)");
  std::vector<std::string> sname{ "x", "y", "u", "v", "lam" };
  double const GR = 1.;

  M.add_domain( t, FFDom( 0., 1., 1, FFDom::LGR, 4 ) );
  for( auto s : { x, y, u, v, lam } ) M.add_state( s, {t} );
  M.set_evolution_domain( t );
  M.update_ref( x, 0.6 ); M.update_ref( y, 0.8 );
  M.update_ref( u, 0. );  M.update_ref( v, 0. ); M.update_ref( lam, -GR*0.8 );
  M.add_equation( OpP(x,t) - u,            {t}, {TI}, io );
  M.add_equation( OpP(y,t) - v,            {t}, {TI}, io );
  M.add_equation( OpP(u,t) + x*lam,        {t}, {TI}, io );
  M.add_equation( OpP(v,t) + y*lam + GR,   {t}, {TI}, io );
  M.add_equation( x*x + y*y - 1.,          {t}, {FFDom::ALL}, io );     // INDEX 3
  M.add_equation( x - 0.6, {t}, {FFDom::LB}, ii );                      // free data {x, u}
  M.add_equation( u - 0.,  {t}, {FFDom::LB}, ii );
  M.options.DISPLAY_LEVEL = 0;

  if( !M.setup() ){ std::cout << "setup failed: " << M.extract_error() << "\n"; return 1; }

  auto const& R = M.reduction_plan();
  std::cout << "\n--- index-3 pendulum, reduced: max_index=" << R.max_index
            << " resolved=" << R.resolved << "  differential=" << M.n_differential()
            << " algebraic=" << M.n_algebraic() << "\n";

  std::vector<double> z0, F, J;
  if( !M.ic_probe( z0, F, J ) ){ std::cout << "no IC system was built\n"; return 1; }
  size_t const n = M.ic_unknowns().size();
  auto const& dotOf = M.ic_dot_of();

  std::vector<std::string> ucol;
  for( size_t i = 0; i < M.nx(); ++i ) ucol.push_back( sname[i] );
  for( size_t d = 0; d < dotOf.size(); ++d ) ucol.push_back( "d" + sname[ dotOf[d] ] );

  std::cout << "\n    seed and residual at it:\n";
  for( size_t i = 0; i < n; ++i )
    std::cout << "      " << std::setw(5) << ucol[i] << " = " << std::setw(12) << std::fixed
              << std::setprecision(6) << z0[i] << "      F[" << i << "] = "
              << std::scientific << std::setprecision(3) << F[i] << "\n";

  std::cout << "\n    Jacobian (rows = IC rows, cols = unknowns):\n          ";
  for( size_t c = 0; c < n; ++c ) std::cout << std::setw(9) << ucol[c];
  std::cout << "\n";
  for( size_t r = 0; r < n; ++r ){
    std::cout << "    r" << std::setw(2) << r << "   ";
    for( size_t c = 0; c < n; ++c )
      std::cout << std::setw(9) << std::fixed << std::setprecision(3) << J[ r + n*c ];
    std::cout << "\n";
  }

  // --- rank by Gaussian elimination with FULL pivoting; report what runs out
  std::vector<double> A = J;
  std::vector<size_t> rowOf( n ), colOf( n );
  for( size_t i = 0; i < n; ++i ){ rowOf[i] = i; colOf[i] = i; }
  size_t rank = 0;
  double const tol = 1e-10;
  for( size_t k = 0; k < n; ++k ){
    double best = 0.; size_t br = k, bc = k;
    for( size_t r = k; r < n; ++r )
      for( size_t c = k; c < n; ++c ){
        double const a = std::fabs( A[ rowOf[r] + n*colOf[c] ] );
        if( a > best ){ best = a; br = r; bc = c; }
      }
    if( best < tol ){
      std::cout << "\n    RANK DEFICIENT: pivot " << k << " is " << std::scientific << best
                << " (tol " << tol << ")\n    rank = " << rank << " of " << n
                << "   nullity = " << ( n - rank ) << "\n";
      std::cout << "    rows never pivoted:";
      for( size_t r = k; r < n; ++r ) std::cout << " r" << rowOf[r];
      std::cout << "\n    unknowns never pinned:";
      for( size_t c = k; c < n; ++c ) std::cout << " " << ucol[ colOf[c] ];
      std::cout << "\n";
      break;
    }
    std::swap( rowOf[k], rowOf[br] ); std::swap( colOf[k], colOf[bc] );
    double const p = A[ rowOf[k] + n*colOf[k] ];
    for( size_t r = k+1; r < n; ++r ){
      double const f = A[ rowOf[r] + n*colOf[k] ] / p;
      if( f == 0. ) continue;
      for( size_t c = k; c < n; ++c ) A[ rowOf[r] + n*colOf[c] ] -= f * A[ rowOf[k] + n*colOf[c] ];
    }
    ++rank;
  }
  if( rank == n ) std::cout << "\n    FULL RANK (" << rank << " of " << n
                            << ") -- the Jacobian at the seed is NOT singular\n";

  // ---- is it the REDUCTION that loses update_ref, or does the seed ignore it generally?
  // Same states and same refs, but NO algebraic constraint, so nothing is reduced.
  { std::cout << "\n--- control: the same 4 ODEs with NO constraint (nothing reduced)\n";
    FFGraph G2; DAESLV_IDAS M2; M2.set_dag( &G2 );
    FFVar t2 = G2.add_var("t"), x2 = G2.add_var("x(t)"), y2 = G2.add_var("y(t)"),
          u2 = G2.add_var("u(t)"), v2 = G2.add_var("v(t)");
    M2.add_domain( t2, FFDom( 0., 1., 1, FFDom::LGR, 4 ) );
    for( auto sv : { x2, y2, u2, v2 } ) M2.add_state( sv, {t2} );
    M2.set_evolution_domain( t2 );
    M2.update_ref( x2, 0.6 ); M2.update_ref( y2, 0.8 );
    M2.update_ref( u2, 0.3 ); M2.update_ref( v2, 0.4 );
    M2.add_equation( OpP(x2,t2) - u2,  {t2}, {TI}, io );
    M2.add_equation( OpP(y2,t2) - v2,  {t2}, {TI}, io );
    M2.add_equation( OpP(u2,t2) + x2,  {t2}, {TI}, io );
    M2.add_equation( OpP(v2,t2) + y2,  {t2}, {TI}, io );
    // a NON-explicit initial row, so an IC system is actually built
    M2.add_equation( x2*x2 + y2*y2 - 1., {t2}, {FFDom::LB}, ii );
    M2.add_equation( y2 - 0.8, {t2}, {FFDom::LB}, ii );
    M2.add_equation( u2 - 0.3, {t2}, {FFDom::LB}, ii );
    M2.add_equation( v2 - 0.4, {t2}, {FFDom::LB}, ii );
    M2.options.DISPLAY_LEVEL = 0;
    if( !M2.setup() ){ std::cout << "    setup failed: " << M2.extract_error() << "\n"; return 1; }
    std::cout << "    reduction: " << M2.reduction_plan().assigns.size() << " assign(s)\n";
    std::vector<double> z2, F2, J2;
    if( !M2.ic_probe( z2, F2, J2 ) ){ std::cout << "    no IC system\n"; return 1; }
    std::cout << "    model's own ref():";
    for( auto const& vv : { x2, y2, u2, v2 } ) std::cout << " " << std::fixed << std::setprecision(3) << M2.ref( vv );
    std::cout << "    nx=" << M2.nx() << " var_state=" << M2.var_state().size() << "\n";
    std::cout << "    seed =";
    for( double const q : z2 ) std::cout << " " << std::fixed << std::setprecision(3) << q;
    std::cout << "    (refs were 0.6 0.8 0.3 0.4)\n"; }

  std::cout << "\n";
  return 0;
}
