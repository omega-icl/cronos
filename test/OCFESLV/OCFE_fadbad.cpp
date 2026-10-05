// Copyright (C) Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// OCFE_fadbad.cpp
// ----------------
// Regression test for FFOCFERES FAD & BAD overloads:
//   * FADType<double>  : values via OCFESLV::eval<double> and forward sensitivities via OCFESLV::deriv<double>
//   * BADType<double>  : values via OCFESLV::eval<double> and reverse sensitivities via OCFESLV::deriv<double>
//   * FADType<FFVar>   : symbolic forward AD through FFGradOCFERES
//
// All derivative paths are compared against OCFESLV::deriv<double>().

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#ifndef OCFE_OCFESLV_HEADER
#define OCFE_OCFESLV_HEADER "ocfeslv.hpp"
#endif
#include OCFE_OCFESLV_HEADER

#include "ffocfe.hpp"

using namespace mc;

namespace {

static double u_exact( double t, double U0 )
{
  return U0 * std::exp( -t );
}

static double max_abs_diff( std::vector<double> const& a,
                            std::vector<double> const& b,
                            size_t* imax = nullptr )
{
  if( a.size() != b.size() )
    throw std::runtime_error( "max_abs_diff: size mismatch" );
  double err = 0.;
  size_t idx = 0;
  for( size_t i=0; i<a.size(); ++i ){
    double const ei = std::abs( a[i] - b[i] );
    if( ei > err ){
      err = ei;
      idx = i;
    }
  }
  if( imax ) *imax = idx;
  return err;
}

static bool report_check( std::string const& label, double err, double tol )
{
  bool const pass = err <= tol;
  std::cout << std::left << std::setw(48) << label
            << " maxerr=" << std::scientific << std::setprecision(3) << err
            << " tol=" << tol << "  " << ( pass? "PASS": "FAIL" ) << "\n";
  return pass;
}

struct CaseData
{
  std::unique_ptr<FFGraph> model_dag;
  std::unique_ptr<OCFESLV>   oc;
  std::unique_ptr<FFGraph> ext_dag;

  std::vector<FFVar> x;       // states, inputs, constants for the external op
  std::vector<FFVar> y;       // residual/output rows from FFOCFERES
  std::vector<double> xval;
  std::vector<double> yref;
  std::vector<double> jref_dense; // row-major nRow x nArg; constants are zero
};

static CaseData make_case()
{
  CaseData c;
  double const U0 = 2.0;

  c.model_dag = std::make_unique<FFGraph>();
  FFVar t   = c.model_dag->add_var( "t" );
  FFVar u   = c.model_dag->add_var( "u(t)" );
  FFVar U0v = c.model_dag->add_var( "U0" );

  FFPartial OpP;
  FFVar ODE = OpP( u, t ) + u;
  FFVar IC  = u - U0v;

  c.oc = std::make_unique<OCFESLV>( c.model_dag.get() );
  c.oc->options.DISPLAY_LEVEL = 0;
  c.oc->options.SOLVE.MARCHING = false;   // direct eval/deriv reference test: keep the full
                                          // (uncollapsed) t-partition -- do not auto-march.
  c.oc->add_domain  ( t, FFDom( 0., 1., 3, FFDom::LGR, 5 ) );
  c.oc->add_state   ( u, {t} );
  c.oc->update_ref  ( u, [&]( OCFESLV::t_Coord const& coord ){
    return u_exact( coord.at( t ), U0 );
  } );
  c.oc->set_constant( {U0v}, {U0} );
  c.oc->add_equation( ODE, {t}, {FFDom::ALL-FFDom::LB},
                      OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  c.oc->add_equation( IC,  {t}, {FFDom::LB},
                      OCFESLV::EqnOptions( OCFESLV::EqnRole::INITIAL, 0 ) );
  c.oc->add_output  ( u, {t}, {1.0} );
  c.oc->set_evolution_domain( t );

  if( !c.oc->setup() )
    throw std::runtime_error( "OCFESLV::setup() failed" );

  size_t const nSta = c.oc->n_colloc_sta();
  size_t const nInp = c.oc->n_colloc_inp();
  size_t const nCst = c.oc->var_constant().size();
  size_t const nEqn = c.oc->n_colloc_eqn();
  size_t const nFct = c.oc->n_colloc_fct();
  size_t const nRow = c.oc->n_colloc_rows();
  size_t const nArg = nSta + nInp + nCst;

  if( nInp != 0 || nCst != 1 || nFct != 1 )
    throw std::runtime_error( "Unexpected test dimensions" );

  std::vector<double> var, inp;
  double const cst[] = { U0 };
  if( !c.oc->init( var, inp, cst ) )
    throw std::runtime_error( "OCFESLV::init() failed" );

  c.ext_dag = std::make_unique<FFGraph>();
  c.x.reserve( nArg );
  for( size_t i=0; i<nSta; ++i ){
    std::ostringstream os;
    os << "z" << i;
    c.x.push_back( c.ext_dag->add_var( os.str() ) );
  }
  for( size_t i=0; i<nInp; ++i ){
    std::ostringstream os;
    os << "q" << i;
    c.x.push_back( c.ext_dag->add_var( os.str() ) );
  }
  c.x.push_back( c.ext_dag->add_var( "U0ext" ) );

  c.xval = var;
  c.xval.insert( c.xval.end(), inp.begin(), inp.end() );
  c.xval.push_back( U0 );

  FFOCFERES OpColl;
  FFVar** ppY = OpColl( nSta, c.x.data(),
                        nInp, nInp? c.x.data()+nSta: nullptr,
                        nCst, c.x.data()+nSta+nInp,
                        c.oc.get(), FFOCFERES::COPY, "fadbad_ode" );
  c.y.resize( nRow );
  for( size_t i=0; i<nRow; ++i ) c.y[i] = *ppY[i];

  // Reference values from OCFESLV::eval<double>().
  std::vector<double> eqn( nEqn, 0. ), fct( nFct, 0. );
  if( !c.oc->eval( eqn.data(), fct.data(), var.data(), nullptr, cst ) )
    throw std::runtime_error( "OCFESLV::eval<double>() failed" );
  c.yref = eqn;
  c.yref.insert( c.yref.end(), fct.begin(), fct.end() );

  // Reference sparse derivatives from OCFESLV::deriv<double>(), densified with
  // columns [states, inputs, constants]; constant columns are expected zero.
  std::vector<size_t> nnz( nRow, 0 );
  if( !c.oc->deriv( nnz.data(), nullptr ) )
    throw std::runtime_error( "OCFESLV::deriv pattern counts failed" );

  std::vector< std::vector<size_t> > col_store( nRow );
  std::vector<size_t*> col_ptr( nRow, nullptr );
  for( size_t i=0; i<nRow; ++i ){
    col_store[i].resize( nnz[i] );
    col_ptr[i] = col_store[i].empty()? nullptr: col_store[i].data();
  }
  if( !c.oc->deriv( nnz.data(), col_ptr.data() ) )
    throw std::runtime_error( "OCFESLV::deriv pattern columns failed" );

  std::vector<double> grad_eqn( c.oc->n_colloc_eqn_nnz(), 0. );
  std::vector<double> grad_fct( c.oc->n_colloc_fct_nnz(), 0. );
  if( !c.oc->deriv( grad_eqn.empty()? nullptr: grad_eqn.data(),
                    grad_fct.empty()? nullptr: grad_fct.data(),
                    var.data(), nullptr, cst ) )
    throw std::runtime_error( "OCFESLV::deriv<double>() failed" );

  std::vector<double> grad_stack = grad_eqn;
  grad_stack.insert( grad_stack.end(), grad_fct.begin(), grad_fct.end() );

  c.jref_dense.assign( nRow * nArg, 0. );
  size_t iel = 0;
  for( size_t irow=0; irow<nRow; ++irow ){
    for( size_t inz=0; inz<nnz[irow]; ++inz, ++iel ){
      size_t const jcol = col_store[irow][inz];
      if( jcol < nArg ) c.jref_dense[irow*nArg+jcol] = grad_stack[iel];
    }
  }
  if( iel != grad_stack.size() )
    throw std::runtime_error( "Unexpected stacked derivative length" );

  return c;
}

static FFOCFERES const* external_op( CaseData const& c )
{
  if( c.y.empty() || !c.y.front().opdef().first ) return nullptr;
  return dynamic_cast<FFOCFERES const*>( c.y.front().opdef().first );
}

} // namespace

int main()
{
  double const tol_val = 1e-11;
  double const tol_der = 1e-10;
  bool ok = true;

  std::cout << "========== OCFE_fadbad: FFOCFERES FADBAD overloads ==========\n";

  CaseData c = make_case();
  size_t const nRow = c.y.size();
  size_t const nArg = c.x.size();

  std::cout << "dimensions: nArg=" << nArg
            << " nRow=" << nRow
            << " nEqn=" << c.oc->n_colloc_eqn()
            << " nFct=" << c.oc->n_colloc_fct()
            << " nnz=" << c.oc->n_colloc_nnz() << "\n\n";

  FFOCFERES const* op = external_op( c );
  if( !op ){
    std::cerr << "ERROR: dependent is not an FFOCFERES operation\n";
    return 1;
  }

  // Baseline external-op double evaluation.
  std::vector<double> y_double( nRow, 0. );
  op->eval( static_cast<unsigned>( nRow ), y_double.data(),
            static_cast<unsigned>( nArg ), c.xval.data(), nullptr );
  ok &= report_check( "FFOCFERES::eval<double> vs OCFESLV::eval", 
                      max_abs_diff( y_double, c.yref ), tol_val );

  // FADType<double>: seed all arguments and compare the generated
  // dense forward Jacobian against OCFESLV::deriv<double>().
  std::vector<FADType<double>> x_Fd( nArg ), y_Fd( nRow );
  for( size_t i=0; i<nArg; ++i ){
    x_Fd[i] = c.xval[i];
    x_Fd[i].diff( static_cast<unsigned>( i ), static_cast<unsigned>( nArg ) );
  }
  op->eval( static_cast<unsigned>( nRow ), y_Fd.data(),
            static_cast<unsigned>( nArg ), x_Fd.data(), nullptr );
  std::vector<double> y_Fd_val( nRow, 0. );
  std::vector<double> jFd_dense( nRow * nArg, 0. );
  for( size_t i=0; i<nRow; ++i ){
    y_Fd_val[i] = y_Fd[i].val();
    for( size_t j=0; j<nArg; ++j )
      jFd_dense[i*nArg+j] = y_Fd[i].deriv( static_cast<unsigned>( j ) );
  }
  ok &= report_check( "FADType<double> values vs OCFESLV::eval", 
                      max_abs_diff( y_Fd_val, c.yref ), tol_val );
  size_t imax_Fd = 0;
  double const err_Fd = max_abs_diff( jFd_dense, c.jref_dense, &imax_Fd );
  ok &= report_check( "FADType<double> gradients vs OCFESLV::deriv", 
                      err_Fd, tol_der );
  if( err_Fd > tol_der ){
    std::cerr << "  worst F<double> derivative entry: row=" << imax_Fd/nArg
              << " col=" << imax_Fd%nArg
              << " fadbad=" << jFd_dense[imax_Fd]
              << " ref=" << c.jref_dense[imax_Fd] << "\n";
  }

  // BADType<double>: evaluate values once, then recompute one scalar reverse
  // sweep per residual/output row to assemble a dense Jacobian.
  std::vector<BADType<double>> x_Bd_val( nArg ), y_Bd_val_ad( nRow );
  for( size_t i=0; i<nArg; ++i ) x_Bd_val[i] = c.xval[i];
  op->eval( static_cast<unsigned>( nRow ), y_Bd_val_ad.data(),
            static_cast<unsigned>( nArg ), x_Bd_val.data(), nullptr );
  std::vector<double> y_Bd_val( nRow, 0. );
  for( size_t i=0; i<nRow; ++i ) y_Bd_val[i] = y_Bd_val_ad[i].val();
  ok &= report_check( "BADType<double> values vs OCFESLV::eval", 
                      max_abs_diff( y_Bd_val, c.yref ), tol_val );

  std::vector<double> jBd_dense( nRow * nArg, 0. );
  for( size_t irow=0; irow<nRow; ++irow ){
    std::vector<BADType<double>> x_Bd( nArg ), y_Bd( nRow );
    for( size_t j=0; j<nArg; ++j ) x_Bd[j] = c.xval[j];
    op->eval( static_cast<unsigned>( nRow ), y_Bd.data(),
              static_cast<unsigned>( nArg ), x_Bd.data(), nullptr );
    y_Bd[irow].diff( 0, 1 );
    for( size_t j=0; j<nArg; ++j )
      jBd_dense[irow*nArg+j] = x_Bd[j].d(0);
  }
  size_t imax_Bd = 0;
  double const err_Bd = max_abs_diff( jBd_dense, c.jref_dense, &imax_Bd );
  ok &= report_check( "BADType<double> gradients vs OCFESLV::deriv", 
                      err_Bd, tol_der );
  if( err_Bd > tol_der ){
    std::cerr << "  worst B<double> derivative entry: row=" << imax_Bd/nArg
              << " col=" << imax_Bd%nArg
              << " fadbad=" << jBd_dense[imax_Bd]
              << " ref=" << c.jref_dense[imax_Bd] << "\n";
  }

  // FADType<FFVar>: FFGraph::SFAD triggers the overload and creates an
  // FFGradOCFERES derivative DAG. Evaluate that DAG in double arithmetic.
  auto sder = c.ext_dag->SFAD( c.y, c.x );
  auto const& irow = std::get<0>( sder );
  auto const& jcol = std::get<1>( sder );
  auto const& dvar = std::get<2>( sder );

  std::vector<double> dval( dvar.size(), 0. );
  if( !dvar.empty() ){
    c.ext_dag->eval( static_cast<unsigned>( dvar.size() ), dvar.data(), dval.data(),
                     static_cast<unsigned>( nArg ), c.x.data(), c.xval.data() );
  }

  std::vector<double> jsfad_dense( nRow * nArg, 0. );
  for( size_t k=0; k<dvar.size(); ++k ){
    if( irow[k] >= nRow || jcol[k] >= nArg ){
      std::cerr << "ERROR: SFAD returned out-of-range entry ("
                << irow[k] << "," << jcol[k] << ")\n";
      return 1;
    }
    jsfad_dense[irow[k]*nArg+jcol[k]] += dval[k];
  }

  size_t imax = 0;
  double const err_sfad = max_abs_diff( jsfad_dense, c.jref_dense, &imax );
  ok &= report_check( "FADType<FFVar>/SFAD values vs OCFESLV::deriv", 
                      err_sfad, tol_der );
  if( err_sfad > tol_der ){
    std::cerr << "  worst derivative entry: row=" << imax/nArg
              << " col=" << imax%nArg
              << " sfad=" << jsfad_dense[imax]
              << " ref=" << c.jref_dense[imax] << "\n";
  }

  // Dense FAD is an additional check that the same F<FFVar> external overload
  // agrees with the sparse derivative reference.
  std::vector<FFVar> dfad = c.ext_dag->FAD( c.y, c.x );
  std::vector<double> dfad_val( dfad.size(), 0. );
  if( !dfad.empty() ){
    c.ext_dag->eval( static_cast<unsigned>( dfad.size() ), dfad.data(), dfad_val.data(),
                     static_cast<unsigned>( nArg ), c.x.data(), c.xval.data() );
  }
  ok &= report_check( "FADType<FFVar>/FAD values vs OCFESLV::deriv", 
                      max_abs_diff( dfad_val, c.jref_dense ), tol_der );

  std::cout << "\nOverall: " << ( ok? "PASS": "FAIL" ) << "\n";
  return ok? 0: 1;
}
