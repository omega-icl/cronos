// Copyright (C) Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// OCFE_ffipdae.cpp
// ---------------------
// Regression test for FFOCFERES as an external DAG operation.
//
// The test builds one OCFESLV for a first-order advection PDE with two lumped
// inputs, then compares:
//
//   1) direct OCFESLV::eval / OCFESLV::deriv;
//   2) DAG evaluation of FFOCFERES with SHALLOW policy;
//   3) DAG sparse automatic differentiation of FFOCFERES with SHALLOW policy;
//   4) DAG evaluation of FFOCFERES with COPY policy;
//   5) DAG sparse automatic differentiation of FFOCFERES with COPY policy.
//
// Derivative testing deliberately uses FFGraph::SFAD(row_ext,xp), followed by
// FFGraph::eval() on the derivative DAG returned by SFAD.  The test does not
// instantiate FFGradOCFERES directly; FFGradOCFERES is created only through
// FFOCFERES::deriv() during SFAD.
//
// No constants are used deliberately.  The FFOCFERES::eval<FFDep> path in the
// current header can only use numeric constants, whereas this test focuses on
// state/input coupling and the external-operation derivative hook.

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

#ifndef OCFE_OCFESLV_HEADER
#define OCFE_OCFESLV_HEADER "ocfeslv.hpp"
#endif
#include OCFE_OCFESLV_HEADER

#include "ffocfe.hpp"

using namespace mc;

static double u_exact
( double t, double x, double c, double u0, double xf )
{
  return u0 * std::sin( 2. * PI / xf * ( x - c * t ) );
}

static double max_abs_diff
( std::vector<double> const& a, std::vector<double> const& b )
{
  if( a.size() != b.size() ) return std::numeric_limits<double>::infinity();
  double out = 0.;
  for( size_t i=0; i<a.size(); ++i )
    out = std::max( out, std::abs( a[i] - b[i] ) );
  return out;
}

static void print_check
( std::string const& label, double const err, double const tol, bool& ok )
{
  bool const pass = std::isfinite(err) && err <= tol;
  ok &= pass;
  std::cout << std::left << std::setw(66) << label
            << " maxerr=" << std::scientific << std::setprecision(6) << err
            << " tol=" << tol << "  " << ( pass? "PASS": "FAIL" ) << "\n";
}

static void print_bool_check
( std::string const& label, bool const pass, bool& ok )
{
  ok &= pass;
  std::cout << std::left << std::setw(66) << label
            << ( pass? " PASS": " FAIL" ) << "\n";
}

static std::set<size_t> dep_cols
( FFDep const& d )
{
  std::set<size_t> out;
  for( auto const& kv : d.dep() ) out.insert( static_cast<size_t>( kv.first ) );
  return out;
}

static bool run_external_policy
( std::string const& label, int const policy, OCFESLV& oc,
  size_t const nVar, size_t const nInp, size_t const nRows, size_t const nCol,
  size_t const nNnz,
  std::vector<double> const& xAll,
  std::vector<double> const& row_ref,
  std::vector<double> const& jac_ref,
  std::vector< std::vector<size_t> > const& col_store,
  std::vector<size_t> const& nnz )
{
  bool ok = true;
  std::cout << "\n-- FFOCFERES policy: " << label << " --\n";

  FFGraph DAGext;
  std::vector<FFVar> z( nVar );
  std::vector<FFVar> p( nInp );
  for( size_t i=0; i<nVar; ++i ){
    std::ostringstream os; os << label << "_z" << i;
    z[i] = DAGext.add_var( os.str() );
  }
  for( size_t i=0; i<nInp; ++i ){
    std::ostringstream os; os << label << "_p" << i;
    p[i] = DAGext.add_var( os.str() );
  }

  std::vector<FFVar> xp;
  xp.reserve( nCol );
  xp.insert( xp.end(), z.begin(), z.end() );
  xp.insert( xp.end(), p.begin(), p.end() );

  FFOCFERES OpColl;
  std::vector<FFVar> row_ext = OpColl( z, p, {}, &oc, policy, label );
  if( row_ext.size() != nRows ){
    std::cerr << "ERROR: FFOCFERES(" << label << ") returned "
              << row_ext.size() << " rows, expected " << nRows << "\n";
    return false;
  }

  std::vector<double> row_val( nRows, 0.0 );
  DAGext.eval( row_ext, row_val, xp, xAll );
  print_check( "FFOCFERES value [" + label + "] vs OCFESLV::eval",
               max_abs_diff( row_val, row_ref ), 1e-10, ok );

  // Differentiate the external operation sparsely, then evaluate the derivative
  // DAG.  This is the path that should create FFGradOCFERES internally via
  // FFOCFERES::deriv(); the test does not instantiate FFGradOCFERES.
  auto sJac = DAGext.SFAD( row_ext, xp );
  std::vector<unsigned> const& irow = std::get<0>( sJac );
  std::vector<unsigned> const& jcol = std::get<1>( sJac );
  std::vector<FFVar>    const& drow = std::get<2>( sJac );

  print_bool_check( "SFAD nonzero count [" + label + "]",
                    drow.size() == nNnz, ok );
  if( irow.size() != drow.size() || jcol.size() != drow.size() ){
    std::cerr << "ERROR: inconsistent SFAD tuple sizes for " << label << "\n";
    return false;
  }

  std::vector<double> dval( drow.size(), 0.0 );
  DAGext.eval( drow, dval, xp, xAll );

  std::map< std::pair<unsigned,unsigned>, double > sfad_map;
  bool indices_ok = true;
  for( size_t k=0; k<drow.size(); ++k ){
    if( irow[k] >= nRows || jcol[k] >= nCol ){
      indices_ok = false;
      continue;
    }
    auto const key = std::make_pair( irow[k], jcol[k] );
    if( sfad_map.count( key ) ) indices_ok = false;
    sfad_map[key] = dval[k];
  }
  print_bool_check( "SFAD indices unique/in range [" + label + "]", indices_ok, ok );

  std::vector<double> jac_val( nRows*nCol, 0.0 );
  for( auto const& kv : sfad_map )
    jac_val[ size_t(kv.first.first)*nCol + size_t(kv.first.second) ] = kv.second;

  print_check( "SFAD derivative DAG eval [" + label + "] vs OCFESLV::deriv",
               max_abs_diff( jac_val, jac_ref ), 1e-10, ok );

  bool pattern_ok = true;
  for( size_t i=0; i<nRows; ++i ){
    std::set<unsigned> got;
    for( auto const& kv : sfad_map )
      if( kv.first.first == i ) got.insert( kv.first.second );

    std::set<unsigned> ref;
    for( size_t j=0; j<nnz[i]; ++j ) ref.insert( static_cast<unsigned>( col_store[i][j] ) );

    if( got != ref ){
      pattern_ok = false;
      std::cerr << "SFAD sparsity mismatch [" << label << "] at row " << i
                << ": got nnz=" << got.size()
                << " ref nnz=" << ref.size() << "\n";
      break;
    }
  }
  print_bool_check( "SFAD sparsity pattern [" + label + "] vs OCFESLV cached pattern",
                    pattern_ok, ok );

  // Check FFDep propagation through the value operation as a separate sparsity
  // route; this does not differentiate the DAG, but catches value-op dependency
  // regressions for both OCFESLV ownership policies.
  std::vector<FFDep> xdep( nCol );
  for( size_t i=0; i<nCol; ++i ) xdep[i].indep( static_cast<int>( i ) );
  std::vector<FFDep> row_dep( nRows );
  DAGext.eval( row_ext, row_dep, xp, xdep );

  bool dep_pattern_ok = true;
  for( size_t i=0; i<nRows; ++i ){
    std::set<size_t> got = dep_cols( row_dep[i] );
    std::set<size_t> ref( col_store[i].begin(), col_store[i].end() );
    if( got != ref ){
      dep_pattern_ok = false;
      std::cerr << "FFDep sparsity mismatch [" << label << "] at row " << i
                << ": got nnz=" << got.size()
                << " ref nnz=" << ref.size() << "\n";
      break;
    }
  }
  print_bool_check( "FFOCFERES FFDep sparsity [" + label + "] vs OCFESLV cached pattern",
                    dep_pattern_ok, ok );

  return ok;
}

int main()
{
  bool ok = true;

  std::cout << "\n========== OCFE_ffipdae: external DAG PDE collocation =========="
            << "\n";

  // ---------------------------------------------------------------------
  // 1. Build the reference OCFESLV.
  // ---------------------------------------------------------------------
  FFGraph DAG;
  FFVar c   = DAG.add_var( "c" );
  FFVar u0v = DAG.add_var( "u0" );
  FFVar t   = DAG.add_var( "t" );
  FFVar x   = DAG.add_var( "x" );
  FFVar u   = DAG.add_var( "u(t,x)" );

  double const tf = 1.;
  double const xf = 1.;
  double const c_ref  = 0.5;
  double const u0_ref = 1.25;

  FFPartial OpP;
  FFVar PDE  = OpP( u, t ) + c * OpP( u, x );
  FFVar INI  = u - u0v * sin( 2. * PI / xf * x );
  FFVar BCIN = u + u0v * sin( 2. * PI * c / xf * t );

  size_t const n_el_t = 1;
  size_t const n_nd_t = 4;
  size_t const n_el_x = 1;
  size_t const n_nd_x = 4;

  FFDom dom_t( 0., tf, n_el_t, FFDom::LGR, n_nd_t );
  FFDom dom_x( 0., xf, n_el_x, FFDom::LGL, n_nd_x );
  if( !dom_t.set_nodes() || !dom_x.set_nodes() ){
    std::cerr << "ERROR: node construction failed\n";
    return 1;
  }

  auto phys_node = []( FFDom const& dom, size_t const iel, size_t const inode ){
    return dom.lo_dom + dom.w_elem * ( static_cast<double>( iel )
         + 0.5 * ( dom.lgnodes()[inode] + 1.0 ) );
  };
  double const t_obs = phys_node( dom_t, 0, 2 );
  double const x_obs = phys_node( dom_x, 0, 3 );

  OCFESLV oc( &DAG );
  oc.options.DISPLAY_LEVEL   = 0;
  oc.options.REDUCE.ORDER    = OCFESLV::Options::RED_NONE;
  oc.options.CLASSIFY.MODE        = OCFESLV::Options::CLASS_AUTO;
  oc.options.INTERFACE.IMPOSITION = OCFESLV::Options::IC_STRONG;
  oc.options.INTERFACE.TYPE  = OCFESLV::Options::IC_AUTO;

  oc.add_domain( t, dom_t );
  oc.add_domain( x, dom_x );
  oc.add_state ( u, {t,x} );
  oc.add_input ( {c,u0v} );
  oc.update_ref( c,   c_ref );
  oc.update_ref( u0v, u0_ref );
  oc.update_ref( u, [&]( OCFESLV::t_Coord const& coord ){
    return u_exact( coord.at( t ), coord.at( x ), c_ref, u0_ref, xf );
  } );
  oc.set_evolution_domain( t );

  oc.add_equation( PDE,  {t,x}, {FFDom::ALL-FFDom::LB,FFDom::ALL-FFDom::LB},
                   OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc.add_equation( INI,  {t,x}, {FFDom::LB,FFDom::ALL},
                   OCFESLV::EqnOptions( OCFESLV::EqnRole::INITIAL, 0 ) );
  oc.add_equation( BCIN, {t,x}, {FFDom::ALL-FFDom::LB,FFDom::LB},
                   OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );

  // Output rows are included to test that FFOCFERES stacks equation rows
  // followed by scalar output rows exactly like OCFESLV::n_colloc_rows().
  oc.add_output( u,   {t,x}, {t_obs,x_obs} );
  oc.add_output( u*u, {t,x}, {t_obs,x_obs} );

  if( !oc.setup() ){
    std::cerr << "ERROR: OCFESLV::setup failed\n";
    return 1;
  }

  size_t const nVar  = oc.n_colloc_sta();
  size_t const nInp  = oc.n_colloc_inp();
  size_t const nEqn  = oc.n_colloc_eqn();
  size_t const nFct  = oc.n_colloc_fct();
  size_t const nRows = oc.n_colloc_rows();
  size_t const nCol  = nVar + nInp;
  size_t const nNnz  = oc.n_colloc_eqn_nnz() + oc.n_colloc_fct_nnz();

  std::cout << "OCFESLV sizes: nVar=" << nVar
            << " nInp=" << nInp
            << " nEqn=" << nEqn
            << " nFct=" << nFct
            << " nRows=" << nRows
            << " nnz=" << nNnz << "\n";

  if( nInp != 2 || nRows != nEqn+nFct ){
    std::cerr << "ERROR: unexpected OCFESLV dimensions\n";
    return 1;
  }

  // ---------------------------------------------------------------------
  // 2. Create a nontrivial test point: initialise from references, then
  //    perturb states and inputs so residuals and input sensitivities are nonzero.
  // ---------------------------------------------------------------------
  std::vector<double> var, inp;
  if( !oc.init( var, inp, nullptr ) ){
    std::cerr << "ERROR: OCFESLV::init failed\n";
    return 1;
  }
  if( var.size() != nVar || inp.size() != nInp ){
    std::cerr << "ERROR: OCFESLV::init returned inconsistent vector sizes\n";
    return 1;
  }

  for( size_t i=0; i<var.size(); ++i )
    var[i] += 1e-3 * std::cos( 0.31 * double(i) + 0.7 );
  inp[0] += 1.25e-2;
  inp[1] -= 8.0e-2;

  std::vector<double> xAll;
  xAll.reserve( nCol );
  xAll.insert( xAll.end(), var.begin(), var.end() );
  xAll.insert( xAll.end(), inp.begin(), inp.end() );

  // ---------------------------------------------------------------------
  // 3. Direct OCFESLV evaluation and derivatives.
  // ---------------------------------------------------------------------
  std::vector<double> eqn_ref( nEqn, 0.0 );
  std::vector<double> fct_ref( nFct, 0.0 );
  if( !oc.eval( eqn_ref.data(), fct_ref.data(), var.data(), inp.data(), nullptr ) ){
    std::cerr << "ERROR: OCFESLV::eval failed\n";
    return 1;
  }

  std::vector<double> row_ref;
  row_ref.reserve( nRows );
  row_ref.insert( row_ref.end(), eqn_ref.begin(), eqn_ref.end() );
  row_ref.insert( row_ref.end(), fct_ref.begin(), fct_ref.end() );

  std::vector<size_t> nnz( nRows, 0 );
  if( !oc.deriv( nnz.data(), nullptr ) ){
    std::cerr << "ERROR: OCFESLV::deriv sparsity count failed\n";
    return 1;
  }
  std::vector< std::vector<size_t> > col_store( nRows );
  std::vector<size_t*> col_ptr( nRows, nullptr );
  for( size_t i=0; i<nRows; ++i ){
    col_store[i].resize( nnz[i] );
    col_ptr[i] = col_store[i].data();
  }
  if( !oc.deriv( nnz.data(), col_ptr.data() ) ){
    std::cerr << "ERROR: OCFESLV::deriv sparsity columns failed\n";
    return 1;
  }

  size_t nnz_sum = std::accumulate( nnz.begin(), nnz.end(), size_t(0) );
  if( nnz_sum != nNnz ){
    std::cerr << "ERROR: n_colloc_*_nnz mismatch: rows sum " << nnz_sum
              << " vs reported " << nNnz << "\n";
    return 1;
  }

  std::vector<double> grad_eqn( oc.n_colloc_eqn_nnz(), 0.0 );
  std::vector<double> grad_fct( oc.n_colloc_fct_nnz(), 0.0 );
  if( !oc.deriv( grad_eqn.data(), grad_fct.data(), var.data(), inp.data(), nullptr ) ){
    std::cerr << "ERROR: OCFESLV::deriv values failed\n";
    return 1;
  }

  std::vector<double> grad_ref;
  grad_ref.reserve( nNnz );
  grad_ref.insert( grad_ref.end(), grad_eqn.begin(), grad_eqn.end() );
  grad_ref.insert( grad_ref.end(), grad_fct.begin(), grad_fct.end() );

  std::vector<double> jac_ref( nRows * nCol, 0.0 );
  size_t kflat = 0;
  for( size_t i=0; i<nRows; ++i ){
    for( size_t j=0; j<nnz[i]; ++j, ++kflat ){
      if( col_store[i][j] >= nCol ){
        std::cerr << "ERROR: sparsity column out of range at row " << i
                  << ": " << col_store[i][j] << " >= " << nCol << "\n";
        return 1;
      }
      jac_ref[i*nCol + col_store[i][j]] = grad_ref[kflat];
    }
  }

  // ---------------------------------------------------------------------
  // 4. External operation graphs for SHALLOW and COPY policies.
  // ---------------------------------------------------------------------
  ok &= run_external_policy( "SHALLOW", FFOCFERES::SHALLOW, oc,
                             nVar, nInp, nRows, nCol, nNnz,
                             xAll, row_ref, jac_ref, col_store, nnz );

  ok &= run_external_policy( "COPY", FFOCFERES::COPY, oc,
                             nVar, nInp, nRows, nCol, nNnz,
                             xAll, row_ref, jac_ref, col_store, nnz );

  std::cout << "\nExternal collocation operation checks: "
            << ( ok? "PASS": "FAIL" ) << "\n";
  return ok? 0: 1;
}
