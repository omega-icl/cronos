// Copyright (C) Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// OCFE_parallel.cpp
// ----------------------------------
// Regression/diagnostic test for deep-copying DAGs that contain FFOCFERES
// and FFGradOCFERES external operations and then evaluating the copies on
// separate std::threads.
//
// The test builds a tiny ODE collocation model
//
//     u_t + u = 0,   u(0) = U0,   t in [0,1],
//
// wraps it as an FFOCFERES external operation in an outer DAG, differentiates
// that outer DAG using sparse forward AD (SFAD), then deep-copies both
//
//   1) the residual DAG, and
//   2) the sparse derivative DAG containing FFGradOCFERES,
//
// and evaluates all copies concurrently.
//
// Important ownership convention checked here:
//   - FFOp::data remains the canonical source identity pointer used by the
//     generic external-operation comparator/deduplication logic;
//   - FFOCFERES/FFGradOCFERES::_pOCFESLV is the executable environment and
//     must be distinct in every deep-copied DAG, under COPY and under SHALLOW
//     alike: since 2026-10-03 an operation re-inserted into ANOTHER DAG owns a
//     deep copy of its solver whatever its policy (FFBaseOCFE::_into_other_dag),
//     so that FFGraph::veval's per-thread DAG copies never share a solver.
//     SHALLOW only governs the ORIGINAL operation's relation to the user's
//     solver.  Both policies are therefore validated serially and in parallel.
//
// Build example:
//   g++ -std=c++17 -O0 -g -pthread -DMC__USE_THREAD -I/path/to/cronos/src \
//       OCFE_parallel.cpp -o OCFE_parallel \
//       -larmadillo
//
// Useful diagnostics:
//   compile with -DCRONOS__FFIPDAE_TRACE or -DMC__FFUNC_EXTERN_DEBUG
//   valgrind --tool=memcheck ./OCFE_parallel 8 100

#include <algorithm>
#include <atomic>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
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

static double max_abs_diff
( std::vector<double> const& a, std::vector<double> const& b )
{
  if( a.size() != b.size() )
    throw std::runtime_error( "max_abs_diff size mismatch" );

  double out = 0.;
  for( size_t i = 0; i < a.size(); ++i )
    out = std::max( out, std::abs( a[i] - b[i] ) );
  return out;
}

template <typename Op>
static Op const* find_external_op( std::vector<FFVar> const& dep )
{
  for( auto const& var : dep ){
    FFOp const* op = var.opdef().first;
    if( auto const* ptr = dynamic_cast<Op const*>( op ) )
      return ptr;
  }
  return nullptr;
}

struct SourceDAG
{
  std::unique_ptr<FFGraph> dag;
  std::vector<FFVar> x;        // external-operation arguments: states, inputs, constants
  std::vector<double> xval;    // numerical values for x

  std::vector<FFVar> r;        // residual/output rows: FFOCFERES
  std::vector<double> rval;    // serial reference r(x)

  std::vector<unsigned> drow;  // sparse derivative row indices from SFAD(r,x)
  std::vector<unsigned> dcol;  // sparse derivative col indices from SFAD(r,x)
  std::vector<FFVar> dr;       // sparse derivative values: FFGradOCFERES
  std::vector<double> drval;   // serial reference dr(x)

  FFOCFERES const* residual_op() const
  { return find_external_op<FFOCFERES>( r ); }

  FFGradOCFERES const* derivative_op() const
  { return find_external_op<FFGradOCFERES>( dr ); }
};

struct EvalDAG
{
  std::unique_ptr<FFGraph> dag;
  std::string label;
  std::vector<FFVar> x;
  std::vector<FFVar> y;
  std::vector<double> xval;
  std::vector<double> yref;

  FFOCFERES const* residual_op() const
  { return find_external_op<FFOCFERES>( y ); }

  FFGradOCFERES const* derivative_op() const
  { return find_external_op<FFGradOCFERES>( y ); }
};

static EvalDAG copy_case
( SourceDAG const& src, std::vector<FFVar> const& dep,
  std::vector<double> const& ref, std::string const& label )
{
  EvalDAG dst;
  dst.label = label;
  dst.dag = std::make_unique<FFGraph>();
  dst.x.resize( src.x.size() );
  dst.y.resize( dep.size() );

  // Copy variables first so that the external operation copied below binds to
  // the destination DAG's variables rather than creating duplicate arguments.
  dst.dag->insert( src.dag.get(), static_cast<unsigned>( src.x.size() ),
                   src.x.data(), dst.x.data() );
  dst.dag->insert( src.dag.get(), static_cast<unsigned>( dep.size() ),
                   dep.data(), dst.y.data() );

  dst.xval = src.xval;
  dst.yref = ref;
  return dst;
}

static void eval_case
( EvalDAG const& item, std::vector<double>& y )
{
  y.assign( item.y.size(), std::numeric_limits<double>::quiet_NaN() );
  item.dag->eval( static_cast<unsigned>( item.y.size() ), item.y.data(), y.data(),
                  static_cast<unsigned>( item.x.size() ), item.x.data(), item.xval.data() );
}

static SourceDAG make_source_dag( int collPolicy, std::string const& opname )
{
  double const U0 = 2.0;

  // Inner/user DAG defining the collocation model.
  auto model_dag = std::make_unique<FFGraph>();
  FFVar t   = model_dag->add_var( "t" );
  FFVar u   = model_dag->add_var( "u(t)" );
  FFVar U0v = model_dag->add_var( "U0" );

  FFPartial OpP;
  FFVar ODE = OpP( u, t ) + u;
  FFVar IC  = u - U0v;

  auto oc = std::make_unique<OCFESLV>( model_dag.get() );
  oc->options.DISPLAY_LEVEL = 0;
  oc->add_domain  ( t, FFDom( 0., 1., 3, FFDom::LGR, 5 ) );
  oc->add_state   ( u, {t} );
  oc->update_ref  ( u, [&]( OCFESLV::t_Coord const& coord ){
    return u_exact( coord.at( t ), U0 );
  } );
  oc->set_constant( {U0v}, {U0} );
  oc->add_equation( ODE, {t}, {FFDom::ALL-FFDom::LB},
                    OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
  oc->add_equation( IC,  {t}, {FFDom::LB},
                    OCFESLV::EqnOptions( OCFESLV::EqnRole::INITIAL, 0 ) );
  oc->set_evolution_domain( t );

  if( !oc->setup() )
    throw std::runtime_error( "OCFESLV::setup() failed" );

  size_t const nSta = oc->n_colloc_sta();
  size_t const nInp = oc->n_colloc_inp();
  size_t const nCst = oc->var_constant().size();
  size_t const nRow = oc->n_colloc_rows();
  if( nInp != 0 || nCst != 1 )
    throw std::runtime_error( "Unexpected test dimensions" );

  std::vector<double> var, inp;
  if( !oc->init( var, inp, &U0 ) )
    throw std::runtime_error( "OCFESLV::init() failed" );
  if( var.size() != nSta || inp.size() != nInp )
    throw std::runtime_error( "Unexpected OCFESLV::init() vector size" );

  // Outer DAG containing the external collocation operation.
  SourceDAG src;
  src.dag = std::make_unique<FFGraph>();
  src.x.reserve( nSta + nInp + nCst );
  for( size_t i = 0; i < nSta; ++i ){
    std::ostringstream nm;
    nm << "z" << i;
    src.x.push_back( src.dag->add_var( nm.str() ) );
  }
  src.x.push_back( src.dag->add_var( "U0ext" ) );

  src.xval = var;
  src.xval.push_back( U0 );

  FFOCFERES OpColl;
  FFVar** ppR = OpColl( nSta, src.x.data(),
                        nInp, nullptr,
                        nCst, src.x.data() + nSta + nInp,
                        oc.get(), collPolicy, opname );
  src.r.resize( nRow );
  for( size_t i = 0; i < nRow; ++i ) src.r[i] = *ppR[i];

  src.rval.assign( src.r.size(), 0. );
  src.dag->eval( static_cast<unsigned>( src.r.size() ), src.r.data(), src.rval.data(),
                 static_cast<unsigned>( src.x.size() ), src.x.data(), src.xval.data() );

  // Build the sparse derivative DAG. This is the path that inserts an
  // FFGradOCFERES external operation through FFOCFERES::deriv().
  std::tie( src.drow, src.dcol, src.dr ) = src.dag->SFAD( src.r, src.x );
  if( src.dr.empty() )
    throw std::runtime_error( "SFAD returned an empty derivative DAG" );

  src.drval.assign( src.dr.size(), 0. );
  src.dag->eval( static_cast<unsigned>( src.dr.size() ), src.dr.data(), src.drval.data(),
                 static_cast<unsigned>( src.x.size() ), src.x.data(), src.xval.data() );

  // Keep the user model alive deliberately.  The executable environment used
  // by the inserted FFOCFERES/FFGradOCFERES operations should nevertheless
  // be their owned/deep-copied _pOCFESLV when COPY semantics are used.
  static std::vector<std::unique_ptr<FFGraph>> keep_model_dags;
  static std::vector<std::unique_ptr<OCFESLV>>   keep_ocfeslvs;
  keep_model_dags.push_back( std::move( model_dag ) );
  keep_ocfeslvs.push_back( std::move( oc ) );

  return src;
}

template <typename Op>
static bool print_identity_diagnostics
( char const* title, Op const* src_op, std::vector<EvalDAG> const& copies,
  Op const* (EvalDAG::*get_op)() const )
{
  if( !src_op ){
    std::cerr << "ERROR: source " << title << " dependent has no expected external operation\n";
    return false;
  }

  std::cout << "\n" << title << " external-operation identity diagnostics:\n";
  std::cout << "  source op=" << src_op
            << "  pOCFESLV=" << src_op->pOCFESLV()
            << "  FFOp::data=" << src_op->data << "\n";

  bool ok = true;
  for( size_t i = 0; i < copies.size(); ++i ){
    Op const* op = (copies[i].*get_op)();
    if( !op ){
      std::cerr << "ERROR: copy " << i << " has no expected " << title << " external operation\n";
      return false;
    }

    // every copy must execute on its OWN solver (a deep copy), whatever the source's policy
    bool const unique_exec_oc = op->pOCFESLV() && op->pOCFESLV() != src_op->pOCFESLV();
    bool const shared_exec_oc = op->pOCFESLV() && op->pOCFESLV() == src_op->pOCFESLV();
    bool const exec_policy_ok = unique_exec_oc;
    bool const same_identity  = op->data == src_op->data;
    bool const identity_not_exec = op->data != static_cast<void const*>( op->pOCFESLV() );

    std::cout << "  copy[" << i << "] op=" << op
              << "  pOCFESLV=" << op->pOCFESLV()
              << "  FFOp::data=" << op->data
              << "  unique_exec_oc=" << ( unique_exec_oc ? "yes" : "NO" )
              << "  shared_exec_oc=" << ( shared_exec_oc ? "yes" : "NO" )
              << "  same_identity=" << ( same_identity ? "yes" : "NO" )
              << "  identity_not_exec=" << ( identity_not_exec ? "yes" : "NO" )
              << "\n";

    ok = ok && exec_policy_ok && same_identity;
  }
  return ok;
}

static bool serial_check
( std::vector<EvalDAG> const& copies, double tol )
{
  std::cout << "\nSerial copy checks:\n";
  bool ok = true;
  for( size_t i = 0; i < copies.size(); ++i ){
    std::vector<double> y;
    eval_case( copies[i], y );
    double const err = max_abs_diff( y, copies[i].yref );
    std::cout << "  " << copies[i].label << " copy[" << i << "] max|y-yref|="
              << std::scientific << err << ( err <= tol ? " PASS" : " FAIL" ) << "\n";
    ok = ok && err <= tol;
  }
  return ok;
}

static bool parallel_check
( std::vector<EvalDAG> const& copies, size_t nRepeat, double tol )
{
  std::cout << "\nParallel evaluation check for " << copies.front().label << ": ";

  std::atomic<bool> go{ false };
  std::atomic<bool> fail{ false };
  std::atomic<size_t> first_bad_thread{ static_cast<size_t>( -1 ) };
  std::atomic<double> first_bad_err{ 0.0 };

  std::vector<std::thread> pool;
  pool.reserve( copies.size() );
  for( size_t tid = 0; tid < copies.size(); ++tid ){
    pool.emplace_back( [&, tid]() {
      while( !go.load( std::memory_order_acquire ) )
        std::this_thread::yield();

      std::vector<double> y;
      for( size_t rep = 0; rep < nRepeat && !fail.load( std::memory_order_relaxed ); ++rep ){
        eval_case( copies[tid], y );
        double const err = max_abs_diff( y, copies[tid].yref );
        if( !( err <= tol ) ){
          first_bad_thread.store( tid, std::memory_order_relaxed );
          first_bad_err.store( err, std::memory_order_relaxed );
          fail.store( true, std::memory_order_release );
          break;
        }
      }
    } );
  }

  go.store( true, std::memory_order_release );
  for( auto& th : pool ) th.join();

  bool const ok = !fail.load();
  std::cout << ( ok ? "PASS" : "FAIL" );
  if( !ok )
    std::cout << "  first_bad_thread=" << first_bad_thread.load()
              << "  err=" << std::scientific << first_bad_err.load();
  std::cout << "\n";
  return ok;
}

} // namespace

int main( int argc, char** argv )
{
  size_t const nThread = argc > 1 ? static_cast<size_t>( std::stoul( argv[1] ) ) : 4;
  size_t const nRepeat = argc > 2 ? static_cast<size_t>( std::stoul( argv[2] ) ) : 20;
  double const tol = 1e-11;

  std::cout << "========== FFOCFERES / FFGradOCFERES parallel deep-copy diagnostic =========="
            << "\nthreads=" << nThread << " repeats/thread=" << nRepeat << "\n";

  SourceDAG source = make_source_dag( FFOCFERES::COPY, "parallel_ode_copy" );
  std::cout << "source residual rows        : " << source.r.size() << "\n";
  std::cout << "source sparse derivative nnz: " << source.dr.size() << "\n";

  std::vector<EvalDAG> residual_copies;
  std::vector<EvalDAG> derivative_copies;
  residual_copies.reserve( nThread );
  derivative_copies.reserve( nThread );
  for( size_t i = 0; i < nThread; ++i ){
    residual_copies.push_back( copy_case( source, source.r,  source.rval,  "residual" ) );
    derivative_copies.push_back( copy_case( source, source.dr, source.drval, "derivative" ) );
  }

  bool const residual_identity_ok = print_identity_diagnostics<FFOCFERES>(
    "FFOCFERES", source.residual_op(), residual_copies, &EvalDAG::residual_op );
  bool const derivative_identity_ok = print_identity_diagnostics<FFGradOCFERES>(
    "FFGradOCFERES", source.derivative_op(), derivative_copies, &EvalDAG::derivative_op );

  bool const residual_serial_ok   = serial_check( residual_copies, tol );
  bool const derivative_serial_ok = serial_check( derivative_copies, tol );

  bool const residual_parallel_ok   = parallel_check( residual_copies,   nRepeat, tol );
  bool const derivative_parallel_ok = parallel_check( derivative_copies, nRepeat, tol );

  // SHALLOW primal: the source operation refers to the user's solver, but every copy into another DAG owns a
  // deep copy (2026-10-03) -- validated exactly as COPY, serially and in parallel, for the residual and the
  // derivative operations
  std::cout << "\nSHALLOW policy check (copies own their solver, as under COPY):\n";
  SourceDAG shallow_source = make_source_dag( FFOCFERES::SHALLOW, "parallel_ode_shallow" );
  std::vector<EvalDAG> shallow_residual_copies, shallow_derivative_copies;
  shallow_residual_copies.reserve( nThread );
  shallow_derivative_copies.reserve( nThread );
  for( size_t i = 0; i < nThread; ++i ){
    shallow_residual_copies.push_back(
      copy_case( shallow_source, shallow_source.r,  shallow_source.rval,  "shallow residual" ) );
    shallow_derivative_copies.push_back(
      copy_case( shallow_source, shallow_source.dr, shallow_source.drval, "shallow derivative" ) );
  }
  bool const shallow_residual_identity_ok = print_identity_diagnostics<FFOCFERES>(
    "FFOCFERES(SHALLOW)", shallow_source.residual_op(), shallow_residual_copies, &EvalDAG::residual_op );
  bool const shallow_derivative_identity_ok = print_identity_diagnostics<FFGradOCFERES>(
    "FFGradOCFERES(SHALLOW primal)", shallow_source.derivative_op(),
    shallow_derivative_copies, &EvalDAG::derivative_op );
  bool const shallow_residual_serial_ok     = serial_check( shallow_residual_copies, tol );
  bool const shallow_derivative_serial_ok   = serial_check( shallow_derivative_copies, tol );
  bool const shallow_residual_parallel_ok   = parallel_check( shallow_residual_copies,   nRepeat, tol );
  bool const shallow_derivative_parallel_ok = parallel_check( shallow_derivative_copies, nRepeat, tol );

  std::cout << "\nInterpretation:\n";
  std::cout << "  residual identity check     : " << ( residual_identity_ok ? "PASS" : "FAIL" ) << "\n";
  std::cout << "  derivative identity check   : " << ( derivative_identity_ok ? "PASS" : "FAIL" ) << "\n";
  std::cout << "  residual serial values      : " << ( residual_serial_ok ? "PASS" : "FAIL" ) << "\n";
  std::cout << "  derivative serial values    : " << ( derivative_serial_ok ? "PASS" : "FAIL" ) << "\n";
  std::cout << "  residual parallel values    : " << ( residual_parallel_ok ? "PASS" : "FAIL" ) << "\n";
  std::cout << "  derivative parallel values  : " << ( derivative_parallel_ok ? "PASS" : "FAIL" ) << "\n";
  std::cout << "  shallow residual identity   : " << ( shallow_residual_identity_ok ? "PASS" : "FAIL" ) << "\n";
  std::cout << "  shallow derivative identity : " << ( shallow_derivative_identity_ok ? "PASS" : "FAIL" ) << "\n";
  std::cout << "  shallow residual serial     : " << ( shallow_residual_serial_ok ? "PASS" : "FAIL" ) << "\n";
  std::cout << "  shallow derivative serial   : " << ( shallow_derivative_serial_ok ? "PASS" : "FAIL" ) << "\n";
  std::cout << "  shallow residual parallel   : " << ( shallow_residual_parallel_ok ? "PASS" : "FAIL" ) << "\n";
  std::cout << "  shallow derivative parallel : " << ( shallow_derivative_parallel_ok ? "PASS" : "FAIL" ) << "\n";

  if( !derivative_identity_ok ){
    std::cout << "\nDIAGNOSIS: the derivative DAG contains an FFGradOCFERES operation whose "
                 "deep-copied executable _pOCFESLV is not unique, or whose canonical "
                 "FFOp::data identity was rebound.  For parallel derivative-DAG "
                 "evaluation, FFGradOCFERES must own/deep-copy _pOCFESLV while "
                 "leaving FFOp::data as the comparator identity.\n";
  }

  return ( residual_identity_ok && derivative_identity_ok
        && residual_serial_ok && derivative_serial_ok
        && residual_parallel_ok && derivative_parallel_ok
        && shallow_residual_identity_ok && shallow_derivative_identity_ok
        && shallow_residual_serial_ok && shallow_derivative_serial_ok
        && shallow_residual_parallel_ok && shallow_derivative_parallel_ok ) ? 0 : 1;
}
