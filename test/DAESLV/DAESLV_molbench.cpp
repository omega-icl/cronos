// DAESLV_molbench -- MoL phase 4, the baseline: the cost of one evaluation of the collocated residual Rc(t, xc)
// on the heated rod, by two routes, per mesh size:
//   A. OCFEMOL::eval_residual -- ONE pass of the OCVar environment over the layout (the interior expression
//      evaluated as a vector over all nodes, the plan's rows emitted alongside);
//   B. the hand-expanded DAG rows of DAESLV_mol's leg 2 -- the same spatial residual written out node by node
//      as ordinary DAG rows and evaluated by FFGraph::eval, which is what an integrator sees when the model is
//      discretised by hand (or by DAESLV_scaling's finite differences, which have fewer operations still).
// Reported, not asserted: microseconds per evaluation, the ratio B/A, and the sizes.  The design claim of the
// MoL plan (section 1) is that A scales with the expression, not with the node count; this is where that is
// measured before any external operation exists.
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <random>
#include "ocfemol.hpp"
using namespace mc;

static constexpr double ALPHA = 0.1, QFLUX = 1.0;

struct Mesh
{
  FFDom dom; size_t nel, nnod; std::vector<double> z;
  Mesh( size_t nel_, size_t nnod_ ) : dom( 0., 1., nel_, FFDom::LGL, nnod_ ), nel( nel_ ), nnod( nnod_ )
  {
    dom.set_nodes();
    for( size_t e = 0; e < nel; ++e ){
      auto ze = dom.lgnodes( dom.elem_lo(e), dom.elem_up(e) );
      for( size_t i = 0; i < nnod; ++i ) z.push_back( ze[i] );
    }
  }
  size_t n() const { return nel*nnod; }
  size_t g( size_t e, size_t i ) const { return e*nnod + i; }
};

template <typename F> static double time_us( size_t reps, F&& f )
{
  f();                                                   // warm
  auto const t0 = std::chrono::steady_clock::now();
  for( size_t k = 0; k < reps; ++k ) f();
  auto const t1 = std::chrono::steady_clock::now();
  return std::chrono::duration<double,std::micro>( t1 - t0 ).count() / reps;
}

int main( int argc, char* argv[] )
{
  std::cout.setf( std::ios::unitbuf );
  std::cout << "\n--- DAESLV_molbench: one evaluation of the collocated residual, OCFEMOL vs hand-expanded DAG rows ------\n\n";
  std::cout << "  " << std::setw(6) << "nel" << std::setw(6) << "nnod" << std::setw(8) << "nodes" << std::setw(9) << "unknowns"
            << std::setw(14) << "A: OCFEMOL" << std::setw(14) << "B: DAG rows" << std::setw(9) << "B/A" << std::setw(12) << "DAG ops"
            << std::setw(14) << "A: Jacobian" << std::setw(9) << "nnz\n";
  std::cout << "  " << std::setw(43) << "" << std::setw(14) << "[us/eval]" << std::setw(14) << "[us/eval]" << std::setw(21) << "" << std::setw(14) << "[us/eval]\n";

  size_t const NNOD = 5;
  std::vector<size_t> meshes = { 5, 10, 20, 40, 80, 160 };
  if( argc > 1 ) meshes = { (size_t)std::atoi( argv[1] ) };
  std::mt19937 rng( 7 ); std::uniform_real_distribution<double> U( 0., 1. );

  for( size_t nel : meshes ){
    Mesh const M( nel, NNOD );
    size_t const NG = M.n();
    size_t const reps = std::max<size_t>( 20, 20000 / NG );

    // ---- A. OCFEMOL
    double usA = 0., usJ = 0.; size_t nA = 0, ncolour = 0;
    {
      FFGraph G; OCFEMOL oc( &G ); FFPartial OpP;
      FFVar t = G.add_var("t"), z = G.add_var("z"), u = G.add_var("u");
      oc.add_domain( t, FFDom( 0., 1., 1, FFDom::LGR, 2 ) );
      oc.add_domain( z, FFDom( 0., 1., nel, FFDom::LGL, NNOD ) );
      oc.set_evolution_domain( t );
      oc.add_state( u, {t,z} ); oc.update_ref( u, 0. );
      int const TI = FFDom::ALL - FFDom::LB, ZI = FFDom::ALL - FFDom::LB - FFDom::UB;
      OCFEMOL::EqnOptions io( OCFEMOL::EqnRole::INTERIOR, 0 ), ib( OCFEMOL::EqnRole::BOUNDARY, 0 ), ii( OCFEMOL::EqnRole::INITIAL, 0 );
      oc.add_equation( OpP(u,t) - OpP( ALPHA*OpP(u,z), z ), {t,z}, {TI, ZI},                io );
      oc.add_equation( u,                                   {t,z}, {TI, FFDom::LB},         ib );
      oc.add_equation( ALPHA*OpP(u,z) - QFLUX,              {t,z}, {TI, FFDom::UB},         ib );
      oc.add_equation( u,                                   {t,z}, {FFDom::LB, FFDom::ALL}, ii );
      oc.options.DISPLAY_LEVEL = 0;
      if( !oc.setup() ){ std::cout << "  OCFEMOL::setup failed at nel=" << nel << "\n"; continue; }
      auto const& S = oc.semi_discrete();
      nA = S.n_states;
      std::vector<double> xc( S.n_states ), rc( S.n_rows );
      for( auto& v : xc ) v = U( rng );
      usA = time_us( reps, [&]{ oc.eval_residual( 0.3, xc.data(), rc.data() ); } );
      std::vector<size_t> ptr, col; std::vector<double> val;
      usJ = time_us( std::max<size_t>( 5, reps/4 ), [&]{ oc.deriv_residual( 0.3, xc.data(), rc.data(), ptr, col, val, nullptr ); } );
      std::vector<size_t> p2, c2; oc.jacobian_pattern( p2, c2 ); ncolour = 0;
      // the colour count is in the setup report; recover it from the cost ratio only approximately, so print nnz instead
      ncolour = c2.size();
    }

    // ---- B. the hand-expanded rows (leg 2's construction), the time derivative as a leaf vector xp
    double usB = 0.; size_t nops = 0;
    {
      FFGraph G;
      std::vector<FFVar> u( NG ), D( NG ), xp( NG );
      for( size_t k = 0; k < NG; ++k ){ u[k] = G.add_var(); D[k] = G.add_var(); xp[k] = G.add_var(); }
      auto ddz = [&]( std::vector<FFVar> const& f ){
        std::vector<FFVar> out( NG );
        for( size_t e = 0; e < M.nel; ++e ){
          std::vector<FFVar> fe, dfe;
          for( size_t i = 0; i < M.nnod; ++i ) fe.push_back( f[ M.g(e,i) ] );
          M.dom.eval_lagrange( dfe, 1, fe, 2./M.dom.elem_width(e) );
          for( size_t i = 0; i < M.nnod; ++i ) out[ M.g(e,i) ] = dfe[i];
        }
        return out;
      };
      std::vector<FFVar> const dudz = ddz( u ), dDdz = ddz( D );
      std::vector<FFVar> rows;
      for( size_t e = 0; e < M.nel; ++e )
      for( size_t i = 0; i < M.nnod; ++i ){
        size_t const k = M.g(e,i);
        bool const dup = ( i == 0 && e > 0 ), z0 = ( e == 0 && i == 0 ), z1 = ( e == M.nel-1 && i == M.nnod-1 );
        if( z0 )       rows.push_back( u[k] );
        else if( dup ) rows.push_back( u[k] - u[ M.g(e-1,M.nnod-1) ] );
        else           rows.push_back( xp[k] - ALPHA*dDdz[k] );
        bool const ltrace = ( i == M.nnod-1 && e+1 < M.nel );
        if( z1 )          rows.push_back( ALPHA*D[k] - QFLUX );
        else if( ltrace ) rows.push_back( D[k] - D[ M.g(e+1,0) ] );
        else              rows.push_back( dudz[k] - D[k] );
      }
      std::vector<FFVar> leaves; leaves.insert( leaves.end(), u.begin(), u.end() ); leaves.insert( leaves.end(), D.begin(), D.end() ); leaves.insert( leaves.end(), xp.begin(), xp.end() );
      std::vector<double> x( leaves.size() ), r( rows.size() );
      for( auto& v : x ) v = U( rng );
      FFSubgraph sg = G.subgraph( rows );
      nops = sg.l_op.size();
      std::vector<double> wk;
      usB = time_us( reps, [&]{ G.eval( sg, wk, rows, r, leaves, x ); } );
    }

    std::cout << "  " << std::setw(6) << nel << std::setw(6) << NNOD << std::setw(8) << NG << std::setw(9) << nA
              << std::setw(14) << std::fixed << std::setprecision(1) << usA << std::setw(14) << usB
              << std::setw(9) << std::setprecision(2) << ( usA > 0.? usB/usA: 0. ) << std::setw(12) << nops
              << std::setw(14) << std::setprecision(1) << usJ << std::setw(9) << ncolour << "\n";
  }
  std::cout << "\n  (A scales with the expression and the node count through the OCVar arithmetic; B with the expanded graph.)\n\n";
  return 0;
}
