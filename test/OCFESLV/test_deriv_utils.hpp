#ifndef MC__TEST_DERIV_UTILS_HPP
#define MC__TEST_DERIV_UTILS_HPP

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <utility>
#include <vector>

namespace mc_test {

struct DerivCheckOptions
{
  double fd_abs_step  = 1e-7;
  double fd_rel_step  = 1e-6;
  double abs_tol      = 2e-5;
  double rel_tol      = 5e-4;
  double missing_tol  = 1e-6;
  double res_tol      = 1e-9;
  size_t max_columns  = 0;
  size_t max_print    = 8;
};

inline std::vector<size_t>
_select_columns( std::vector<size_t> cols, size_t const max_columns )
{
  std::sort( cols.begin(), cols.end() );
  cols.erase( std::unique( cols.begin(), cols.end() ), cols.end() );
  if( !max_columns || cols.size() <= max_columns ) return cols;

  std::vector<size_t> out;
  out.reserve( max_columns );
  for( size_t i = 0; i < max_columns; ++i ){
    size_t const pos = ( max_columns == 1 )
                     ? 0
                     : ( i * ( cols.size() - 1 ) ) / ( max_columns - 1 );
    if( out.empty() || out.back() != cols[pos] ) out.push_back( cols[pos] );
  }
  return out;
}

inline bool
check_oc_derivatives
( mc::OCFESLV& oc, std::string const& label,
  std::vector<double> const& var, double const* inp, double const* cst,
  DerivCheckOptions const& opt = DerivCheckOptions() )
{
  size_t const nVar = oc.n_colloc_sta();
  size_t const nInp = oc.n_colloc_inp();
  size_t const nDeriv = nVar + nInp;
  size_t const nEqn = oc.n_colloc_eqn();

  if( var.size() != nVar ){
    std::cerr << "Derivative check [" << label << "] ERROR: var.size()="
              << var.size() << " but n_colloc_sta()=" << nVar << "\n";
    return false;
  }

  std::vector<size_t> nnz( nEqn, 0 );
  if( !oc.deriv( nullptr, nullptr, nnz.data(), nullptr,
                 nullptr, nullptr, nullptr, nullptr,
                 nullptr, nullptr, cst ) ){
    std::cerr << "Derivative check [" << label
              << "] ERROR: deriv() failed during equation sparsity query.\n";
    return false;
  }

  std::vector< std::vector<size_t> > col_store( nEqn );
  std::vector< std::vector<double> > grad_store( nEqn );
  std::vector<size_t*> col_ptr( nEqn, nullptr );
  std::vector<double*> grad_ptr( nEqn, nullptr );

  size_t total_nnz = 0, max_row_nnz = 0;
  for( size_t k = 0; k < nEqn; ++k ){
    total_nnz += nnz[k];
    max_row_nnz = std::max( max_row_nnz, nnz[k] );
    col_store[k].resize( nnz[k] );
    grad_store[k].resize( nnz[k] );
    col_ptr[k]  = col_store[k].data();
    grad_ptr[k] = grad_store[k].data();
  }

  std::vector<double> res_deriv( nEqn, 0.0 );
  if( !oc.deriv( res_deriv.data(), grad_ptr.data(), nnz.data(), col_ptr.data(),
                 nullptr, nullptr, nullptr, nullptr,
                 var.data(), inp, cst ) ){
    std::cerr << "Derivative check [" << label
              << "] ERROR: deriv() failed during equation AD evaluation.\n";
    return false;
  }

  std::vector<double> res_eval( nEqn, 0.0 );
  if( !oc.eval( res_eval.data(), nullptr, var.data(), inp, cst ) ){
    std::cerr << "Derivative check [" << label
              << "] ERROR: eval() failed for residual comparison.\n";
    return false;
  }

  double max_res_diff = 0.0;
  for( size_t k = 0; k < nEqn; ++k )
    max_res_diff = std::max( max_res_diff, std::abs( res_deriv[k] - res_eval[k] ) );

  std::vector< std::vector< std::pair<size_t,size_t> > > col_entries( nDeriv );
  std::vector<size_t> structural_cols;
  structural_cols.reserve( total_nnz );

  bool pattern_sorted = true;
  bool pattern_range = true;
  for( size_t k = 0; k < nEqn; ++k ){
    std::vector<size_t> row_sorted = col_store[k];
    std::sort( row_sorted.begin(), row_sorted.end() );
    if( row_sorted != col_store[k] ) pattern_sorted = false;

    for( size_t j = 0; j < col_store[k].size(); ++j ){
      size_t const col = col_store[k][j];
      if( col >= nDeriv ){
        std::cerr << "Derivative check [" << label << "] ERROR: row " << k
                  << " contains out-of-range stacked column " << col
                  << " but nVar+nInp=" << nDeriv << "\n";
        pattern_range = false;
        continue;
      }
      structural_cols.push_back( col );
      col_entries[col].push_back( {k,j} );
    }
  }
  if( !pattern_range ) return false;

  std::vector<size_t> fd_cols = _select_columns( structural_cols, opt.max_columns );

  double max_abs_err = 0.0, max_rel_err = 0.0, max_missing = 0.0;
  size_t n_bad = 0;
  std::vector<double> vp( var ), vm( var );
  std::vector<double> rp( nEqn ), rm( nEqn );
  std::vector<unsigned char> in_pattern( nEqn, 0 );

  for( size_t const col : fd_cols ){
    if( col >= nVar ) continue; // state-only finite differences in this helper
    double const h = opt.fd_abs_step + opt.fd_rel_step * std::max( 1.0, std::abs( var[col] ) );
    vp[col] += h;
    vm[col] -= h;

    if( !oc.eval( rp.data(), nullptr, vp.data(), inp, cst ) ||
        !oc.eval( rm.data(), nullptr, vm.data(), inp, cst ) ){
      std::cerr << "Derivative check [" << label
                << "] ERROR: eval() failed during finite differences.\n";
      return false;
    }

    std::fill( in_pattern.begin(), in_pattern.end(), 0 );
    for( auto const& e : col_entries[col] ) in_pattern[e.first] = 1;

    for( auto const& e : col_entries[col] ){
      size_t const k = e.first, j = e.second;
      double const fd = ( rp[k] - rm[k] ) / ( 2.0 * h );
      double const ad = grad_store[k][j];
      double const err = std::abs( ad - fd );
      double const scale = std::max( 1.0, std::max( std::abs(ad), std::abs(fd) ) );
      double const rel = err / scale;
      max_abs_err = std::max( max_abs_err, err );
      max_rel_err = std::max( max_rel_err, rel );
      if( err > opt.abs_tol + opt.rel_tol * scale ){
        if( n_bad < opt.max_print ){
          std::cerr << "Derivative check [" << label << "] mismatch: row=" << k
                    << " col=" << col << " AD=" << std::scientific << ad
                    << " FD=" << fd << " abs_err=" << err
                    << " rel_err=" << rel << "\n";
        }
        ++n_bad;
      }
    }

    for( size_t k = 0; k < nEqn; ++k ){
      if( in_pattern[k] ) continue;
      double const fd = ( rp[k] - rm[k] ) / ( 2.0 * h );
      max_missing = std::max( max_missing, std::abs(fd) );
      if( std::abs(fd) > opt.missing_tol ){
        if( n_bad < opt.max_print ){
          std::cerr << "Derivative check [" << label << "] missing pattern entry: row="
                    << k << " col=" << col << " FD=" << std::scientific << fd << "\n";
        }
        ++n_bad;
      }
    }

    vp[col] = var[col];
    vm[col] = var[col];
  }

  bool const ok = pattern_sorted && max_res_diff <= opt.res_tol && n_bad == 0;
  std::cout << "Derivative check [" << label << "]: "
            << "nEqn=" << nEqn
            << " nVar=" << nVar
            << " nInp=" << nInp
            << " nnz=" << total_nnz
            << " max_row_nnz=" << max_row_nnz
            << " fd_cols=" << fd_cols.size() << "/" << _select_columns( structural_cols, 0 ).size()
            << " max|res_deriv-res_eval|=" << std::scientific << max_res_diff
            << " max|AD-FD|=" << max_abs_err
            << " max_rel=" << max_rel_err
            << " max_missing=" << max_missing
            << " sorted_rows=" << ( pattern_sorted ? "yes" : "NO" )
            << "  " << ( ok ? "PASS" : "FAIL" ) << "\n";
  return ok;
}

inline bool
check_oc_output_derivatives
( mc::OCFESLV& oc, std::string const& label,
  std::vector<double> const& var, double const* inp, double const* cst,
  DerivCheckOptions const& opt = DerivCheckOptions() )
{
  size_t const nVar = oc.n_colloc_sta();
  size_t const nInp = oc.n_colloc_inp();
  size_t const nDeriv = nVar + nInp;
  size_t const nOut = oc.n_colloc_fct();

  if( !nOut ){
    std::cout << "Output derivative check [" << label << "]: no outputs - SKIP\n";
    return true;
  }
  if( var.size() != nVar ){
    std::cerr << "Output derivative check [" << label << "] ERROR: var.size()="
              << var.size() << " but n_colloc_sta()=" << nVar << "\n";
    return false;
  }

  std::vector<size_t> nnz( nOut, 0 );
  if( !oc.deriv( nullptr, nullptr, nullptr, nullptr,
                 nullptr, nullptr, nnz.data(), nullptr,
                 nullptr, nullptr, cst ) ){
    std::cerr << "Output derivative check [" << label
              << "] ERROR: deriv() failed during output sparsity query.\n";
    return false;
  }

  std::vector< std::vector<size_t> > col_store( nOut );
  std::vector< std::vector<double> > grad_store( nOut );
  std::vector<size_t*> col_ptr( nOut, nullptr );
  std::vector<double*> grad_ptr( nOut, nullptr );

  size_t total_nnz = 0, max_row_nnz = 0;
  for( size_t k = 0; k < nOut; ++k ){
    total_nnz += nnz[k];
    max_row_nnz = std::max( max_row_nnz, nnz[k] );
    col_store[k].resize( nnz[k] );
    grad_store[k].resize( nnz[k] );
    col_ptr[k]  = col_store[k].data();
    grad_ptr[k] = grad_store[k].data();
  }

  std::vector<double> fct_deriv( nOut, 0.0 );
  if( !oc.deriv( nullptr, nullptr, nullptr, nullptr,
                 fct_deriv.data(), grad_ptr.data(), nnz.data(), col_ptr.data(),
                 var.data(), inp, cst ) ){
    std::cerr << "Output derivative check [" << label
              << "] ERROR: deriv() failed during output AD evaluation.\n";
    return false;
  }

  std::vector<double> fct_eval( nOut, 0.0 );
  if( !oc.eval( nullptr, fct_eval.data(), var.data(), inp, cst ) ){
    std::cerr << "Output derivative check [" << label
              << "] ERROR: eval() failed for output comparison.\n";
    return false;
  }

  double max_fct_diff = 0.0;
  for( size_t k = 0; k < nOut; ++k )
    max_fct_diff = std::max( max_fct_diff, std::abs( fct_deriv[k] - fct_eval[k] ) );

  std::vector< std::vector< std::pair<size_t,size_t> > > col_entries( nDeriv );
  std::vector<size_t> structural_cols;
  structural_cols.reserve( total_nnz );
  bool pattern_sorted = true, pattern_range = true;

  for( size_t k = 0; k < nOut; ++k ){
    std::vector<size_t> row_sorted = col_store[k];
    std::sort( row_sorted.begin(), row_sorted.end() );
    if( row_sorted != col_store[k] ) pattern_sorted = false;
    for( size_t j = 0; j < col_store[k].size(); ++j ){
      size_t const col = col_store[k][j];
      if( col >= nDeriv ){
        std::cerr << "Output derivative check [" << label << "] ERROR: row " << k
                  << " contains out-of-range stacked column " << col
                  << " but nVar+nInp=" << nDeriv << "\n";
        pattern_range = false;
        continue;
      }
      structural_cols.push_back( col );
      col_entries[col].push_back( {k,j} );
    }
  }
  if( !pattern_range ) return false;

  std::vector<size_t> fd_cols = _select_columns( structural_cols, opt.max_columns );
  double max_abs_err = 0.0, max_rel_err = 0.0, max_missing = 0.0;
  size_t n_bad = 0;
  std::vector<double> vp( var ), vm( var ), fp( nOut ), fm( nOut );
  std::vector<unsigned char> in_pattern( nOut, 0 );

  for( size_t const col : fd_cols ){
    if( col >= nVar ) continue; // state-only finite differences in this helper
    double const h = opt.fd_abs_step + opt.fd_rel_step * std::max( 1.0, std::abs( var[col] ) );
    vp[col] += h;
    vm[col] -= h;
    if( !oc.eval( nullptr, fp.data(), vp.data(), inp, cst ) ||
        !oc.eval( nullptr, fm.data(), vm.data(), inp, cst ) ){
      std::cerr << "Output derivative check [" << label
                << "] ERROR: eval() failed during output finite differences.\n";
      return false;
    }

    std::fill( in_pattern.begin(), in_pattern.end(), 0 );
    for( auto const& e : col_entries[col] ) in_pattern[e.first] = 1;

    for( auto const& e : col_entries[col] ){
      size_t const k = e.first, j = e.second;
      double const fd = ( fp[k] - fm[k] ) / ( 2.0 * h );
      double const ad = grad_store[k][j];
      double const err = std::abs( ad - fd );
      double const scale = std::max( 1.0, std::max( std::abs(ad), std::abs(fd) ) );
      double const rel = err / scale;
      max_abs_err = std::max( max_abs_err, err );
      max_rel_err = std::max( max_rel_err, rel );
      if( err > opt.abs_tol + opt.rel_tol * scale ){
        if( n_bad < opt.max_print ){
          std::cerr << "Output derivative check [" << label << "] mismatch: out=" << k
                    << " col=" << col << " AD=" << std::scientific << ad
                    << " FD=" << fd << " abs_err=" << err
                    << " rel_err=" << rel << "\n";
        }
        ++n_bad;
      }
    }

    for( size_t k = 0; k < nOut; ++k ){
      if( in_pattern[k] ) continue;
      double const fd = ( fp[k] - fm[k] ) / ( 2.0 * h );
      max_missing = std::max( max_missing, std::abs(fd) );
      if( std::abs(fd) > opt.missing_tol ){
        if( n_bad < opt.max_print ){
          std::cerr << "Output derivative check [" << label << "] missing pattern entry: out="
                    << k << " col=" << col << " FD=" << std::scientific << fd << "\n";
        }
        ++n_bad;
      }
    }
    vp[col] = var[col];
    vm[col] = var[col];
  }

  bool const ok = pattern_sorted && max_fct_diff <= opt.res_tol && n_bad == 0;
  std::cout << "Output derivative check [" << label << "]: "
            << "nOut=" << nOut
            << " nVar=" << nVar
            << " nInp=" << nInp
            << " nnz=" << total_nnz
            << " max_row_nnz=" << max_row_nnz
            << " fd_cols=" << fd_cols.size() << "/" << _select_columns( structural_cols, 0 ).size()
            << " max|fct_deriv-fct_eval|=" << std::scientific << max_fct_diff
            << " max|AD-FD|=" << max_abs_err
            << " max_rel=" << max_rel_err
            << " max_missing=" << max_missing
            << " sorted_rows=" << ( pattern_sorted ? "yes" : "NO" )
            << "  " << ( ok ? "PASS" : "FAIL" ) << "\n";
  return ok;
}

inline bool
check_oc_derivatives_varinp
( mc::OCFESLV& oc, std::string const& label,
  std::vector<double> const& var,
  std::vector<double> const& inp,
  double const* cst,
  DerivCheckOptions const& opt = DerivCheckOptions() )
{
  size_t const nVar   = oc.n_colloc_sta();
  size_t const nInp   = inp.size();
  size_t const nDeriv = nVar + nInp;
  size_t const nEqn   = oc.n_colloc_eqn();

  if( var.size() != nVar ){
    std::cerr << "Derivative check [" << label << "] ERROR: var.size()="
              << var.size() << " but n_colloc_sta()=" << nVar << "\n";
    return false;
  }

  std::vector<size_t> nnz( nEqn, 0 );
  if( !oc.deriv( nullptr, nullptr, nnz.data(), nullptr,
                 nullptr, nullptr, nullptr, nullptr,
                 nullptr, nullptr, cst ) ){
    std::cerr << "Derivative check [" << label
              << "] ERROR: deriv() failed during equation sparsity query.\n";
    return false;
  }

  std::vector< std::vector<size_t> > col_store( nEqn );
  std::vector< std::vector<double> > grad_store( nEqn );
  std::vector<size_t*> col_ptr( nEqn, nullptr );
  std::vector<double*> grad_ptr( nEqn, nullptr );

  size_t total_nnz = 0, max_row_nnz = 0;
  for( size_t k = 0; k < nEqn; ++k ){
    total_nnz += nnz[k];
    max_row_nnz = std::max( max_row_nnz, nnz[k] );
    col_store[k].resize( nnz[k] );
    grad_store[k].resize( nnz[k] );
    col_ptr[k]  = col_store[k].data();
    grad_ptr[k] = grad_store[k].data();
  }

  std::vector<double> res_deriv( nEqn, 0.0 );
  if( !oc.deriv( res_deriv.data(), grad_ptr.data(), nnz.data(), col_ptr.data(),
                 nullptr, nullptr, nullptr, nullptr,
                 var.data(), inp.empty()? nullptr: inp.data(), cst ) ){
    std::cerr << "Derivative check [" << label
              << "] ERROR: deriv() failed during AD evaluation.\n";
    return false;
  }

  std::vector<double> res_eval( nEqn, 0.0 );
  if( !oc.eval( res_eval.data(), nullptr, var.data(), inp.empty()? nullptr: inp.data(), cst ) ){
    std::cerr << "Derivative check [" << label
              << "] ERROR: eval() failed for residual comparison.\n";
    return false;
  }

  double max_res_diff = 0.0;
  for( size_t k = 0; k < nEqn; ++k )
    max_res_diff = std::max( max_res_diff, std::abs( res_deriv[k] - res_eval[k] ) );

  std::vector< std::vector< std::pair<size_t,size_t> > > col_entries( nDeriv );
  std::vector<size_t> structural_cols;
  structural_cols.reserve( total_nnz );

  bool pattern_sorted = true, pattern_range = true, saw_input_col = false;
  for( size_t k = 0; k < nEqn; ++k ){
    std::vector<size_t> row_sorted = col_store[k];
    std::sort( row_sorted.begin(), row_sorted.end() );
    if( row_sorted != col_store[k] ) pattern_sorted = false;
    for( size_t j = 0; j < col_store[k].size(); ++j ){
      size_t const col = col_store[k][j];
      if( col >= nDeriv ){
        std::cerr << "Derivative check [" << label << "] ERROR: row " << k
                  << " contains out-of-range stacked column " << col
                  << " but nVar+nInp=" << nDeriv << "\n";
        pattern_range = false;
        continue;
      }
      if( col >= nVar ) saw_input_col = true;
      structural_cols.push_back( col );
      col_entries[col].push_back( {k,j} );
    }
  }
  if( !pattern_range ) return false;

  std::vector<size_t> fd_cols = _select_columns( structural_cols, opt.max_columns );
  double max_abs_err = 0.0, max_rel_err = 0.0, max_missing = 0.0;
  size_t n_bad = 0;
  std::vector<double> vp( var ), vm( var ), ip( inp ), im( inp );
  std::vector<double> rp( nEqn ), rm( nEqn );
  std::vector<unsigned char> in_pattern( nEqn, 0 );

  for( size_t const col : fd_cols ){
    double const xj = ( col < nVar ) ? var[col] : inp[col-nVar];
    double const h = opt.fd_abs_step + opt.fd_rel_step * std::max( 1.0, std::abs( xj ) );
    if( col < nVar ){
      vp[col] += h; vm[col] -= h;
    } else {
      ip[col-nVar] += h; im[col-nVar] -= h;
    }

    if( !oc.eval( rp.data(), nullptr, vp.data(), ip.empty()? nullptr: ip.data(), cst ) ||
        !oc.eval( rm.data(), nullptr, vm.data(), im.empty()? nullptr: im.data(), cst ) ){
      std::cerr << "Derivative check [" << label
                << "] ERROR: eval() failed during finite differences.\n";
      return false;
    }

    std::fill( in_pattern.begin(), in_pattern.end(), 0 );
    for( auto const& e : col_entries[col] ) in_pattern[e.first] = 1;

    for( auto const& e : col_entries[col] ){
      size_t const k = e.first, j = e.second;
      double const fd = ( rp[k] - rm[k] ) / ( 2.0 * h );
      double const ad = grad_store[k][j];
      double const err = std::abs( ad - fd );
      double const scale = std::max( 1.0, std::max( std::abs(ad), std::abs(fd) ) );
      double const rel = err / scale;
      max_abs_err = std::max( max_abs_err, err );
      max_rel_err = std::max( max_rel_err, rel );
      if( err > opt.abs_tol + opt.rel_tol * scale ){
        if( n_bad < opt.max_print ){
          std::cerr << "Derivative check [" << label << "] mismatch: row=" << k
                    << " col=" << col << ( col < nVar ? " var" : " inp" )
                    << " AD=" << std::scientific << ad
                    << " FD=" << fd << " abs_err=" << err
                    << " rel_err=" << rel << "\n";
        }
        ++n_bad;
      }
    }

    for( size_t k = 0; k < nEqn; ++k ){
      if( in_pattern[k] ) continue;
      double const fd = ( rp[k] - rm[k] ) / ( 2.0 * h );
      max_missing = std::max( max_missing, std::abs(fd) );
      if( std::abs(fd) > opt.missing_tol ){
        if( n_bad < opt.max_print ){
          std::cerr << "Derivative check [" << label << "] missing pattern entry: row="
                    << k << " col=" << col << ( col < nVar ? " var" : " inp" )
                    << " FD=" << std::scientific << fd << "\n";
        }
        ++n_bad;
      }
    }

    if( col < nVar ){
      vp[col] = var[col]; vm[col] = var[col];
    } else {
      ip[col-nVar] = inp[col-nVar]; im[col-nVar] = inp[col-nVar];
    }
  }

  bool const ok = pattern_sorted && saw_input_col && max_res_diff <= opt.res_tol && n_bad == 0;
  std::cout << "Derivative check [" << label << "]: "
            << "nEqn=" << nEqn
            << " nVar=" << nVar
            << " nInp=" << nInp
            << " nDeriv=" << nDeriv
            << " nnz=" << total_nnz
            << " max_row_nnz=" << max_row_nnz
            << " fd_cols=" << fd_cols.size() << "/" << _select_columns( structural_cols, 0 ).size()
            << " max|res_deriv-res_eval|=" << std::scientific << max_res_diff
            << " max|AD-FD|=" << max_abs_err
            << " max_rel=" << max_rel_err
            << " max_missing=" << max_missing
            << " sorted_rows=" << ( pattern_sorted ? "yes" : "NO" )
            << " input_cols=" << ( saw_input_col ? "yes" : "NO" )
            << "  " << ( ok ? "PASS" : "FAIL" ) << "\n";
  return ok;
}

inline bool
check_oc_derivatives_stacked
( mc::OCFESLV& oc, std::string const& label,
  std::vector<double> const& var,
  std::vector<double> const& inp,
  double const* cst,
  DerivCheckOptions const& opt = DerivCheckOptions() )
{
  return check_oc_derivatives_varinp( oc, label, var, inp, cst, opt );
}

inline bool
check_deriv_pattern_matches_ffdep
( mc::OCFESLV& oc, std::string const& label,
  size_t const nInp, double const* cst )
{
  size_t const nVar   = oc.n_colloc_sta();
  size_t const nDeriv = nVar + nInp;
  size_t const nEqn   = oc.n_colloc_eqn();

  std::vector<mc::FFDep> var_dep( nVar );
  for( size_t i = 0; i < nVar; ++i )
    var_dep[i].indep( static_cast<int>( i ) );

  std::vector<mc::FFDep> inp_dep( nInp );
  for( size_t i = 0; i < nInp; ++i )
    inp_dep[i].indep( static_cast<int>( nVar + i ) );

  std::vector<mc::FFDep> res_dep( nEqn );
  if( !oc.eval( res_dep.data(), nullptr,
                var_dep.empty()? nullptr: var_dep.data(),
                inp_dep.empty()? nullptr: inp_dep.data(), cst ) ){
    std::cerr << "FFDep sparsity [" << label << "] ERROR: eval<FFDep>() failed.\n";
    return false;
  }

  std::vector<size_t> nnz( nEqn, 0 );
  if( !oc.deriv( nullptr, nullptr, nnz.data(), nullptr,
                 nullptr, nullptr, nullptr, nullptr,
                 nullptr, nullptr, nullptr ) ){
    std::cerr << "FFDep sparsity [" << label
              << "] ERROR: deriv() failed during sparsity query.\n";
    return false;
  }

  std::vector< std::vector<size_t> > col_store( nEqn );
  std::vector<size_t*> col_ptr( nEqn, nullptr );
  for( size_t k = 0; k < nEqn; ++k ){
    col_store[k].resize( nnz[k] );
    col_ptr[k] = col_store[k].data();
  }

  if( !oc.deriv( nullptr, nullptr, nnz.data(), col_ptr.data(),
                 nullptr, nullptr, nullptr, nullptr,
                 nullptr, nullptr, nullptr ) ){
    std::cerr << "FFDep sparsity [" << label
              << "] ERROR: deriv() failed while retrieving columns.\n";
    return false;
  }

  size_t total_ffdep = 0, total_deriv = 0, n_mismatch = 0;
  size_t input_deps = 0, nonlinear_deps = 0;

  for( size_t k = 0; k < nEqn; ++k ){
    std::vector<size_t> pat_ffdep;
    pat_ffdep.reserve( res_dep[k].dep().size() );
    for( auto const& it : res_dep[k].dep() ){
      int const icol = it.first;
      auto const typ = it.second;
      if( icol < 0 ) continue;
      size_t const col = static_cast<size_t>( icol );
      if( col >= nDeriv ) continue;
      pat_ffdep.push_back( col );
      if( col >= nVar ) ++input_deps;
      if( typ != mc::FFDep::L ) ++nonlinear_deps;
    }
    std::sort( pat_ffdep.begin(), pat_ffdep.end() );
    pat_ffdep.erase( std::unique( pat_ffdep.begin(), pat_ffdep.end() ), pat_ffdep.end() );

    std::vector<size_t> pat_deriv = col_store[k];
    std::sort( pat_deriv.begin(), pat_deriv.end() );
    pat_deriv.erase( std::unique( pat_deriv.begin(), pat_deriv.end() ), pat_deriv.end() );

    total_ffdep += pat_ffdep.size();
    total_deriv += pat_deriv.size();

    if( pat_ffdep != pat_deriv ){
      if( n_mismatch < 8 ){
        std::cerr << "FFDep sparsity [" << label << "] mismatch row " << k << "\n  FFDep:";
        for( auto c : pat_ffdep ) std::cerr << " " << c;
        std::cerr << "\n  deriv:";
        for( auto c : pat_deriv ) std::cerr << " " << c;
        std::cerr << "\n";
      }
      ++n_mismatch;
    }
  }

  bool const ok = n_mismatch == 0 && ( nInp == 0 || input_deps > 0 );
  std::cout << "FFDep sparsity [" << label << "]: "
            << "nEqn=" << nEqn
            << " nVar=" << nVar
            << " nInp=" << nInp
            << " nDeriv=" << nDeriv
            << " nnz_ffdep=" << total_ffdep
            << " nnz_deriv=" << total_deriv
            << " input_deps=" << input_deps
            << " nonlinear_dep_tags=" << nonlinear_deps
            << " mismatched_rows=" << n_mismatch
            << "  " << ( ok ? "PASS" : "FAIL" ) << "\n";
  return ok;
}

} // namespace mc_test

#endif // MC_TEST_DERIV_UTILS_HPP
