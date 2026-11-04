// Copyright (C) Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the EPL-2.0 with GPL-2.0-or-later as a Secondary License; see the LICENSE file.
// SPDX-License-Identifier: EPL-2.0 OR GPL-2.0-or-later

#ifndef CRONOS__FFOCFEMOL_HPP
#define CRONOS__FFOCFEMOL_HPP

#include <sstream>
#include <stdexcept>
#include "ocfemol.hpp"

namespace mc
{

//! @brief The method-of-lines residual Rc(t, xc) as a DAG external operation
////////////////////////////////////////////////////////////////////////
//! mc::FFMOLRES inserts into a DAG the n_rows outputs rc = Rc(t, xc) of a set-up mc::OCFEMOL, from the
//! n_states + 1 inputs (xc, then t).  Evaluation goes through OCFEMOL::eval_residual; derivatives through
//! OCFEMOL::deriv_residual, sparse, with FFGradMOLRES as the companion operation that produces the Jacobian's
//! entries.  Mirrors FFOCFERES/FFGradOCFERES (ffocfe.hpp), reduced to what an integrator needs: double, forward
//! AD in double, dependence, symbolic insertion, and the sparse derivative.  The mass matrix is not part of the
//! operation: the integrator's residual is  Mc(xc) d_t xc - Rc(t, xc)  with Mc from OCFEMOL::mass().
//! MoL workplan 2026-10-31, phase 4.
////////////////////////////////////////////////////////////////////////
class FFGradMOLRES;
class FFMOLRES
: public FFOp
{
public:
  //! @brief Copy policy for the solver the operation holds
  enum Policy { SHALLOW = 0, COPY = 1, TRANSFER = -1 };

protected:
  OCFEMOL*     _pMOL = nullptr;   //!< the set-up solver
  bool         _own  = false;     //!< whether this operation owns (and deletes) it
  size_t       _nSta = 0;         //!< n_states
  size_t       _nRow = 0;         //!< n_rows
  std::string  _name;
  //! @brief The gradient operation's outputs, built once per operation instance and reused by every deriv call:
  //! the derivative propagation asks for this operation's derivative once per row it differentiates through.
  mutable std::vector<FFVar>  _gradOut;
  mutable std::vector<long>   _gradIn;      //!< the input ids the cached gradient was built for

  FFVar const* _gradient
    ( FFVar const* vVar )
    const;

  FFVar** _set
    ( size_t const nSta, FFVar const* pSta, FFVar const& t, OCFEMOL* pMOL, int policy, std::string const& name )
    {
      if( !pMOL || !pMOL->semi_discrete().ran || nSta != pMOL->semi_discrete().n_states )
        throw std::runtime_error( "FFMOLRES ** the solver is not set up or the state count does not match\n" );
      if( _own && _pMOL ) delete _pMOL;
      _own  = ( policy > 0 );
      _pMOL = pMOL;
      _nSta = nSta;
      _nRow = pMOL->semi_discrete().n_rows;
      _name = name;
      data    = pMOL;
      owndata = false;
      std::vector<FFVar> vIn( pSta, pSta + nSta );  vIn.push_back( t );
      FFVar** ppRes = insert_external_operation( *this, _nRow, vIn.size(), vIn.data() );
      _own = false;
      FFOp* pOp = (*ppRes)->opdef().first;
      if( policy > 0 )      _pMOL = dynamic_cast<FFMOLRES*>( pOp )->_pMOL;
      else if( policy < 0 ) dynamic_cast<FFMOLRES*>( pOp )->_own = true;
      return ppRes;
    }

public:
  FFMOLRES
    ()
    : FFOp( EXTERN )
    { this->sparse = true; }   // MC++ calls the SPARSE deriv only when this is set; unset, it allocates a dense
                                // n_rows x (n_states+1) array per call and calls the dense deriv (FFOCFERES sets it too)

  virtual ~FFMOLRES
    ()
    { if( _own && _pMOL ) delete _pMOL; }

  //! @brief Insert the operation: the n_rows outputs for inputs @p vSta (n_states) and @p t.
  std::vector<FFVar> operator()
    ( std::vector<FFVar> const& vSta, FFVar const& t, OCFEMOL* pMOL, int policy = SHALLOW, std::string const& name = "" )
    {
      FFVar** ppRes = _set( vSta.size(), vSta.data(), t, pMOL, policy, name );
      std::vector<FFVar> vDep( _nRow );
      for( size_t i = 0; i < _nRow; ++i ) vDep[i] = *ppRes[i];
      return vDep;
    }

  OCFEMOL* solver() const { return _pMOL; }
  size_t n_states() const { return _nSta; }
  size_t n_rows()   const { return _nRow; }

  // ---- the DAG's type dispatch
  virtual void feval
    ( std::type_info const& idU, unsigned const nRes, void* vRes, unsigned const nVar, void const* vVar, unsigned const* mVar )
    const
    {
      if( mc::same_type( idU, typeid( FFVar ) ) )
        return eval( nRes, static_cast<FFVar*>(vRes), nVar, static_cast<FFVar const*>(vVar), mVar );
      else if( mc::same_type( idU, typeid( FFDep ) ) )
        return eval( nRes, static_cast<FFDep*>(vRes), nVar, static_cast<FFDep const*>(vVar), mVar );
      else if( mc::same_type( idU, typeid( double ) ) )
        return eval( nRes, static_cast<double*>(vRes), nVar, static_cast<double const*>(vVar), mVar );
      else if( mc::same_type( idU, typeid( FADType<double> ) ) )
        return eval( nRes, static_cast<FADType<double>*>(vRes), nVar, static_cast<FADType<double> const*>(vVar), mVar );
      throw std::runtime_error( "FFMOLRES::feval ** no evaluation method for type " + std::string( idU.name() ) + "\n" );
    }

  // ---- evaluation
  void eval
    ( unsigned const nRes, double* vRes, unsigned const nVar, double const* vVar, unsigned const* mVar )
    const;
  void eval
    ( unsigned const nRes, FADType<double>* vRes, unsigned const nVar, FADType<double> const* vVar, unsigned const* mVar )
    const;
  void eval
    ( unsigned const nRes, FFDep* vRes, unsigned const nVar, FFDep const* vVar, unsigned const* mVar )
    const;
  void eval
    ( unsigned const nRes, FFVar* vRes, unsigned const nVar, FFVar const* vVar, unsigned const* mVar )
    const;
  template <typename U>
  void eval
    ( unsigned const nRes, U* vRes, unsigned const nVar, U const* vVar, unsigned const* mVar )
    const
    { throw std::runtime_error( "FFMOLRES::eval ** evaluation type not supported\n" ); }

  // ---- derivatives: sparse, through FFGradMOLRES
  void deriv
    ( unsigned const nRes, FFVar const* vRes, unsigned const nVar, FFVar const* vVar, FFVar** vDer )
    const;
  void deriv
    ( unsigned const nRes, FFVar const* vRes, unsigned const nVar, FFVar const* vVar, FFVar** vDer,
      size_t* nnz, size_t** colnz )
    const;

  // ---- identity in the DAG
  bool lt
    ( FFOp const* op )
    const
    {
      if( data != op->data ) return data < op->data;
      return false;
    }
  std::string name
    ()
    const
    { std::ostringstream oss; oss << ( _name.empty()? "MOLRES": _name ) << "[" << _nRow << "x" << _nSta + 1 << "]"; return oss.str(); }
  bool commutative
    ()
    const
    { return false; }
};

//! @brief The Jacobian entries of FFMOLRES as a DAG external operation: the nnz values of dRc/dxc in CSR order,
//! then the n_rows values of dRc/dt.
////////////////////////////////////////////////////////////////////////
class FFGradMOLRES
: public FFOp
{
protected:
  OCFEMOL*     _pMOL = nullptr;
  bool         _own  = false;
  size_t       _nSta = 0, _nRow = 0, _nnz = 0;
  std::string  _name;

public:
  FFGradMOLRES
    ()
    : FFOp( EXTERN )
    {}
  virtual ~FFGradMOLRES
    ()
    { if( _own && _pMOL ) delete _pMOL; }

  //! @brief Insert: outputs are the nnz Jacobian values (CSR order of OCFEMOL::jacobian_pattern) followed by dRc/dt.
  FFVar** operator()
    ( size_t const nSta, FFVar const* pSta, FFVar const& t, OCFEMOL* pMOL, bool own, std::string const& name )
    {
      _pMOL = pMOL; _own = false; _nSta = nSta; _nRow = pMOL->semi_discrete().n_rows; _name = name;
      std::vector<size_t> ptr, col;  pMOL->jacobian_pattern( ptr, col );  _nnz = col.size();
      data = pMOL;  owndata = false;
      std::vector<FFVar> vIn( pSta, pSta + nSta );  vIn.push_back( t );
      FFVar** ppRes = insert_external_operation( *this, _nnz + _nRow, vIn.size(), vIn.data() );
      if( own ) dynamic_cast<FFGradMOLRES*>( (*ppRes)->opdef().first )->_own = true;
      return ppRes;
    }

  virtual void feval
    ( std::type_info const& idU, unsigned const nRes, void* vRes, unsigned const nVar, void const* vVar, unsigned const* mVar )
    const
    {
      if( mc::same_type( idU, typeid( FFVar ) ) )
        return eval( nRes, static_cast<FFVar*>(vRes), nVar, static_cast<FFVar const*>(vVar), mVar );
      else if( mc::same_type( idU, typeid( double ) ) )
        return eval( nRes, static_cast<double*>(vRes), nVar, static_cast<double const*>(vVar), mVar );
      throw std::runtime_error( "FFGradMOLRES::feval ** no evaluation method for type " + std::string( idU.name() ) + "\n" );
    }

  void eval
    ( unsigned const nRes, double* vRes, unsigned const nVar, double const* vVar, unsigned const* mVar )
    const
    {
      (void)mVar;
      if( !_pMOL || nVar != _nSta + 1 || nRes != _nnz + _nRow )
        throw std::runtime_error( "FFGradMOLRES::eval double ** size mismatch\n" );
      std::vector<size_t> ptr, col; std::vector<double> val;
      std::vector<double> drdt( _nRow );
      if( !_pMOL->deriv_residual( vVar[_nSta], vVar, nullptr, ptr, col, val, drdt.data() ) )
        throw std::runtime_error( "FFGradMOLRES::eval double ** OCFEMOL::deriv_residual failed\n" );
      if( val.size() != _nnz ) throw std::runtime_error( "FFGradMOLRES::eval double ** pattern changed\n" );
      for( size_t k = 0; k < _nnz; ++k ) vRes[k] = val[k];
      for( size_t r = 0; r < _nRow; ++r ) vRes[_nnz + r] = drdt[r];
    }
  void eval
    ( unsigned const nRes, FFVar* vRes, unsigned const nVar, FFVar const* vVar, unsigned const* mVar )
    const
    {
      (void)mVar;
      FFGradMOLRES op;
      FFVar** ppRes = op( _nSta, vVar, vVar[_nSta], _pMOL, false, _name );
      for( unsigned k = 0; k < nRes; ++k ) vRes[k] = *ppRes[k];
    }
  template <typename U>
  void eval
    ( unsigned const nRes, U* vRes, unsigned const nVar, U const* vVar, unsigned const* mVar )
    const
    { throw std::runtime_error( "FFGradMOLRES::eval ** evaluation type not supported\n" ); }

  bool lt( FFOp const* op ) const { if( data != op->data ) return data < op->data; return false; }
  std::string name() const { std::ostringstream oss; oss << ( _name.empty()? "GradMOLRES": _name ) << "[" << _nnz + _nRow << "]"; return oss.str(); }
  bool commutative() const { return false; }
};


//! @brief The Ini stream as a DAG external operation: every entry of OCFEMOL::ini_stream() (the stream first use, the
//! dropped pairs for the consistency check at t0) evaluated from (xc, t).  Derivative through FFGradMOLINI.
////////////////////////////////////////////////////////////////////////
class FFGradMOLINI
: public FFOp
{
protected:
  OCFEMOL* _pMOL = nullptr; size_t _nSta = 0, _nIni = 0, _nnz = 0; std::string _name;
public:
  FFGradMOLINI() : FFOp( EXTERN ) {}
  FFVar** operator()
    ( size_t const nSta, FFVar const* pSta, FFVar const& t, OCFEMOL* pMOL, std::string const& name )
    {
      _pMOL = pMOL; _nSta = nSta; _nIni = pMOL->ini_stream().size(); _name = name;
      _nnz = 0; for( auto const& e : pMOL->ini_stream() ) _nnz += e.cols.size();
      data = pMOL; owndata = false;
      std::vector<FFVar> vIn( pSta, pSta + nSta ); vIn.push_back( t );
      return insert_external_operation( *this, _nnz, vIn.size(), vIn.data() );
    }
  virtual void feval
    ( std::type_info const& idU, unsigned const nRes, void* vRes, unsigned const nVar, void const* vVar, unsigned const* mVar )
    const
    {
      if( mc::same_type( idU, typeid( FFVar ) ) ) return eval( nRes, static_cast<FFVar*>(vRes), nVar, static_cast<FFVar const*>(vVar), mVar );
      else if( mc::same_type( idU, typeid( double ) ) ) return eval( nRes, static_cast<double*>(vRes), nVar, static_cast<double const*>(vVar), mVar );
      throw std::runtime_error( "FFGradMOLINI::feval ** no evaluation method for type " + std::string( idU.name() ) + "\n" );
    }
  void eval( unsigned const nRes, double* vRes, unsigned const nVar, double const* vVar, unsigned const* mVar ) const
    {
      (void)mVar; if( !_pMOL || nVar != _nSta + 1 || nRes != _nnz ) throw std::runtime_error( "FFGradMOLINI::eval ** size mismatch\n" );
      std::vector<size_t> ptr, col; std::vector<double> val;
      if( !_pMOL->deriv_initial( vVar[_nSta], vVar, nullptr, ptr, col, val ) || val.size() != _nnz ) throw std::runtime_error( "FFGradMOLINI::eval ** deriv_initial failed\n" );
      for( size_t k = 0; k < _nnz; ++k ) vRes[k] = val[k];
    }
  void eval( unsigned const nRes, FFVar* vRes, unsigned const nVar, FFVar const* vVar, unsigned const* mVar ) const
    { (void)mVar; FFGradMOLINI op; FFVar** pp = op( _nSta, vVar, vVar[_nSta], _pMOL, _name ); for( unsigned k = 0; k < nRes; ++k ) vRes[k] = *pp[k]; }
  template <typename U> void eval( unsigned const, U*, unsigned const, U const*, unsigned const* ) const
    { throw std::runtime_error( "FFGradMOLINI::eval ** evaluation type not supported\n" ); }
  bool lt( FFOp const* op ) const { if( data != op->data ) return data < op->data; return false; }
  std::string name() const { std::ostringstream oss; oss << ( _name.empty()? "GradMOLINI": _name ) << "[" << _nnz << "]"; return oss.str(); }
  bool commutative() const { return false; }
};

class FFMOLINI
: public FFOp
{
protected:
  OCFEMOL* _pMOL = nullptr; size_t _nSta = 0, _nIni = 0; std::string _name;
public:
  FFMOLINI() : FFOp( EXTERN ) { this->sparse = true; }   // the sparse deriv, as for FFMOLRES
  //! @brief Insert: one output per entry of OCFEMOL::ini_stream(), in that order.
  std::vector<FFVar> operator()
    ( std::vector<FFVar> const& vSta, FFVar const& t, OCFEMOL* pMOL, std::string const& name = "" )
    {
      if( !pMOL || !pMOL->semi_discrete().ran || vSta.size() != pMOL->semi_discrete().n_states )
        throw std::runtime_error( "FFMOLINI ** the solver is not set up or the state count does not match\n" );
      _pMOL = pMOL; _nSta = vSta.size(); _nIni = pMOL->ini_stream().size(); _name = name;
      data = pMOL; owndata = false;
      std::vector<FFVar> vIn( vSta ); vIn.push_back( t );
      FFVar** pp = insert_external_operation( *this, _nIni, vIn.size(), vIn.data() );
      std::vector<FFVar> out( _nIni ); for( size_t k = 0; k < _nIni; ++k ) out[k] = *pp[k]; return out;
    }
  virtual void feval
    ( std::type_info const& idU, unsigned const nRes, void* vRes, unsigned const nVar, void const* vVar, unsigned const* mVar )
    const
    {
      if( mc::same_type( idU, typeid( FFVar ) ) ) return eval( nRes, static_cast<FFVar*>(vRes), nVar, static_cast<FFVar const*>(vVar), mVar );
      else if( mc::same_type( idU, typeid( FFDep ) ) ) return eval( nRes, static_cast<FFDep*>(vRes), nVar, static_cast<FFDep const*>(vVar), mVar );
      else if( mc::same_type( idU, typeid( double ) ) ) return eval( nRes, static_cast<double*>(vRes), nVar, static_cast<double const*>(vVar), mVar );
      else if( mc::same_type( idU, typeid( FADType<double> ) ) ) return eval( nRes, static_cast<FADType<double>*>(vRes), nVar, static_cast<FADType<double> const*>(vVar), mVar );
      throw std::runtime_error( "FFMOLINI::feval ** no evaluation method for type " + std::string( idU.name() ) + "\n" );
    }
  void eval( unsigned const nRes, double* vRes, unsigned const nVar, double const* vVar, unsigned const* mVar ) const
    { (void)mVar; if( !_pMOL || nVar != _nSta + 1 || nRes != _nIni || !_pMOL->eval_initial( vVar[_nSta], vVar, vRes ) ) throw std::runtime_error( "FFMOLINI::eval double ** failure\n" ); }
  void eval( unsigned const nRes, FADType<double>* vRes, unsigned const nVar, FADType<double> const* vVar, unsigned const* mVar ) const
    {
      (void)mVar; std::vector<double> x( nVar ); for( unsigned i = 0; i < nVar; ++i ) x[i] = vVar[i].val();
      unsigned nDir = 0; for( unsigned i = 0; i < nVar; ++i ) if( vVar[i].depend() ){ nDir = vVar[i].size(); break; }
      std::vector<size_t> ptr, col; std::vector<double> val, v( nRes );
      if( !_pMOL->deriv_initial( x[_nSta], x.data(), v.data(), ptr, col, val ) ) throw std::runtime_error( "FFMOLINI::eval FAD ** failure\n" );
      for( unsigned k = 0; k < nRes; ++k ){
        vRes[k] = v[k]; if( !nDir ) continue;
        vRes[k].diff( 0, nDir ); vRes[k].d( 0 ) = 0.;
        for( unsigned d = 0; d < nDir; ++d ){ double s = 0.; for( size_t j = ptr[k]; j < ptr[k+1]; ++j ) if( vVar[col[j]].depend() ) s += val[j]*vVar[col[j]].deriv( d ); vRes[k].d( d ) = s; }
      }
    }
  void eval( unsigned const nRes, FFDep* vRes, unsigned const nVar, FFDep const* vVar, unsigned const* mVar ) const
    {
      (void)mVar; (void)nVar; auto const& IS = _pMOL->ini_stream();
      for( unsigned k = 0; k < nRes; ++k ){ vRes[k] = FFDep(); for( auto const& [y, c] : IS[k].cols ) vRes[k] += vVar[c]; vRes[k] += vVar[_nSta]; vRes[k].update( FFDep::TYPE::N ); }
    }
  void eval( unsigned const nRes, FFVar* vRes, unsigned const nVar, FFVar const* vVar, unsigned const* mVar ) const
    { (void)mVar; FFMOLINI op; std::vector<FFVar> vs( vVar, vVar + _nSta ); auto out = op( vs, vVar[_nSta], _pMOL, _name ); for( unsigned k = 0; k < nRes; ++k ) vRes[k] = out[k]; }
  template <typename U> void eval( unsigned const, U*, unsigned const, U const*, unsigned const* ) const
    { throw std::runtime_error( "FFMOLINI::eval ** evaluation type not supported\n" ); }
  void deriv( unsigned const nRes, FFVar const* vRes, unsigned const nVar, FFVar const* vVar, FFVar** vDer ) const
    {
      auto const& IS = _pMOL->ini_stream(); FFGradMOLINI gop; FFVar** pg = gop( _nSta, vVar, vVar[_nSta], _pMOL, _name ); size_t j = 0;
      for( unsigned k = 0; k < nRes; ++k ){ for( unsigned i = 0; i < nVar; ++i ) vDer[k][i] = 0.; for( auto const& [y, c] : IS[k].cols ) vDer[k][c] = *pg[j++]; }
    }
  void deriv( unsigned const nRes, FFVar const* vRes, unsigned const nVar, FFVar const* vVar, FFVar** vDer, size_t* nnz, size_t** colnz ) const
    {
      auto const& IS = _pMOL->ini_stream();
      for( unsigned k = 0; k < nRes; ++k ) nnz[k] = IS[k].cols.size();
      if( colnz ) for( unsigned k = 0; k < nRes; ++k ){ if( !colnz[k] ) colnz[k] = new size_t[ std::max<size_t>( nnz[k], 1 ) ]; size_t m = 0; for( auto const& [y, c] : IS[k].cols ) colnz[k][m++] = c; }
      if( !vDer ) return;
      FFGradMOLINI gop; FFVar** pg = gop( _nSta, vVar, vVar[_nSta], _pMOL, _name ); size_t j = 0;
      for( unsigned k = 0; k < nRes; ++k ){ if( !vDer[k] ) vDer[k] = new FFVar[ std::max<size_t>( nnz[k], 1 ) ]; for( size_t m = 0; m < nnz[k]; ++m ) vDer[k][m] = *pg[j++]; }
    }
  bool lt( FFOp const* op ) const { if( data != op->data ) return data < op->data; return false; }
  std::string name() const { std::ostringstream oss; oss << ( _name.empty()? "MOLINI": _name ) << "[" << _nIni << "x" << _nSta + 1 << "]"; return oss.str(); }
  bool commutative() const { return false; }
};

inline FFVar const*
FFMOLRES::_gradient
( FFVar const* vVar )
const
{
  std::vector<long> ids( _nSta + 1 );  for( size_t i = 0; i <= _nSta; ++i ) ids[i] = vVar[i].id().second;
  if( _gradOut.empty() || ids != _gradIn ){
    FFGradMOLRES gop;
    FFVar** ppG = gop( _nSta, vVar, vVar[_nSta], _pMOL, false, _name );
    std::vector<size_t> ptr, col;  _pMOL->jacobian_pattern( ptr, col );
    _gradOut.assign( col.size() + _nRow, FFVar() );
    for( size_t k = 0; k < _gradOut.size(); ++k ) _gradOut[k] = *ppG[k];
    _gradIn = ids;
  }
  return _gradOut.data();
}

// -------------------------------------------------------------------------------------------------------- FFMOLRES
inline void
FFMOLRES::eval
( unsigned const nRes, double* vRes, unsigned const nVar, double const* vVar, unsigned const* mVar )
const
{
  (void)mVar;
  if( !_pMOL || nVar != _nSta + 1 || nRes != _nRow ) throw std::runtime_error( "FFMOLRES::eval double ** size mismatch\n" );
  if( !_pMOL->eval_residual( vVar[_nSta], vVar, vRes ) ) throw std::runtime_error( "FFMOLRES::eval double ** OCFEMOL::eval_residual failed\n" );
}

inline void
FFMOLRES::eval
( unsigned const nRes, FADType<double>* vRes, unsigned const nVar, FADType<double> const* vVar, unsigned const* mVar )
const
{
  (void)mVar;
  if( !_pMOL || nVar != _nSta + 1 || nRes != _nRow ) throw std::runtime_error( "FFMOLRES::eval FAD ** size mismatch\n" );
  std::vector<double> x( nVar );
  for( unsigned i = 0; i < nVar; ++i ) x[i] = vVar[i].val();
  bool dependent = false;  unsigned nDir = 0;
  for( unsigned i = 0; i < nVar; ++i ){
    if( !vVar[i].depend() ) continue;
    if( !dependent ) nDir = vVar[i].size();
    else if( vVar[i].size() != nDir ) throw std::runtime_error( "FFMOLRES::eval FAD ** inconsistent derivative vector sizes\n" );
    dependent = true;
  }
  std::vector<double> r( nRes );
  if( !dependent ){
    if( !_pMOL->eval_residual( x[_nSta], x.data(), r.data() ) ) throw std::runtime_error( "FFMOLRES::eval FAD ** eval failed\n" );
    for( unsigned k = 0; k < nRes; ++k ) vRes[k] = r[k];
    return;
  }
  // one compressed Jacobian, then each output's directions as a sparse combination of the inputs' directions
  std::vector<size_t> ptr, col; std::vector<double> val, drdt( nRes );
  if( !_pMOL->deriv_residual( x[_nSta], x.data(), r.data(), ptr, col, val, drdt.data() ) )
    throw std::runtime_error( "FFMOLRES::eval FAD ** OCFEMOL::deriv_residual failed\n" );
  for( unsigned k = 0; k < nRes; ++k ){
    vRes[k] = r[k];
    vRes[k].diff( 0, nDir );  vRes[k].d( 0 ) = 0.;      // allocate nDir directions
    for( unsigned d = 0; d < nDir; ++d ){
      double s = 0.;
      for( size_t j = ptr[k]; j < ptr[k+1]; ++j ) if( vVar[col[j]].depend() ) s += val[j] * vVar[col[j]].deriv( d );
      if( vVar[_nSta].depend() ) s += drdt[k] * vVar[_nSta].deriv( d );
      vRes[k].d( d ) = s;
    }
  }
}

inline void
FFMOLRES::eval
( unsigned const nRes, FFDep* vRes, unsigned const nVar, FFDep const* vVar, unsigned const* mVar )
const
{
  (void)mVar;
  if( !_pMOL || nVar != _nSta + 1 || nRes != _nRow ) throw std::runtime_error( "FFMOLRES::eval FFDep ** size mismatch\n" );
  std::vector<size_t> ptr, col;  _pMOL->jacobian_pattern( ptr, col );
  for( unsigned k = 0; k < nRes; ++k ){
    vRes[k] = FFDep();
    for( size_t j = ptr[k]; j < ptr[k+1]; ++j ) vRes[k] += vVar[col[j]];
    vRes[k] += vVar[_nSta];                               // every row may depend on t (conservative)
    vRes[k].update( FFDep::TYPE::N );                     // nonlinear in general
  }
}

inline void
FFMOLRES::eval
( unsigned const nRes, FFVar* vRes, unsigned const nVar, FFVar const* vVar, unsigned const* mVar )
const
{
  (void)mVar;
  FFMOLRES op;
  FFVar** ppRes = op._set( _nSta, vVar, vVar[_nSta], _pMOL, SHALLOW, _name );
  for( unsigned k = 0; k < nRes; ++k ) vRes[k] = *ppRes[k];
}

inline void
FFMOLRES::deriv
( unsigned const nRes, FFVar const* vRes, unsigned const nVar, FFVar const* vVar, FFVar** vDer )
const
{
  // dense: every derivative, from the sparse entries
  std::vector<size_t> ptr, col;  _pMOL->jacobian_pattern( ptr, col );
  FFVar const* ppG = _gradient( vVar );                 // one gradient operation per instance, reused
  for( unsigned k = 0; k < nRes; ++k ){
    for( unsigned i = 0; i < nVar; ++i ) vDer[k][i] = 0.;
    for( size_t j = ptr[k]; j < ptr[k+1]; ++j ) vDer[k][ col[j] ] = ppG[j];
    vDer[k][_nSta] = ppG[ col.size() + k ];
  }
}

inline void
FFMOLRES::deriv
( unsigned const nRes, FFVar const* vRes, unsigned const nVar, FFVar const* vVar, FFVar** vDer,
  size_t* nnz, size_t** colnz )
const
{
  std::vector<size_t> ptr, col;  _pMOL->jacobian_pattern( ptr, col );
  for( unsigned k = 0; k < nRes; ++k ) nnz[k] = ( ptr[k+1] - ptr[k] ) + 1;      // + the t column
  if( colnz ){
    for( unsigned k = 0; k < nRes; ++k ){
      if( !colnz[k] ) colnz[k] = new size_t[ nnz[k] ];
      size_t m = 0;
      for( size_t j = ptr[k]; j < ptr[k+1]; ++j ) colnz[k][m++] = col[j];
      colnz[k][m] = _nSta;
    }
  }
  if( !vDer ) return;
  FFVar const* ppG = _gradient( vVar );                 // one gradient operation per instance, reused
  for( unsigned k = 0; k < nRes; ++k ){
    if( !vDer[k] ) vDer[k] = new FFVar[ nnz[k] ];
    size_t m = 0;
    for( size_t j = ptr[k]; j < ptr[k+1]; ++j ) vDer[k][m++] = ppG[j];
    vDer[k][m] = ppG[ col.size() + k ];
  }
}


//! @brief Write the semi-discrete system of a set-up OCFEMOL into a DAE solver's model.
////////////////////////////////////////////////////////////////////////
//! Declares into @p S, on the DAG @p dag: the unknowns xc (one per column of the semi-discrete system, with the
//! owning field's reference as their reference) on the evolution direction @p t with the discretisation @p tdom;
//! one INTERIOR row per row of Rc -- `sum_c Mc[r,c] d_t xc_c + Rc_r` for a row with mass entries, collocated
//! everywhere but the initial point, and `Rc_r` elsewhere, collocated everywhere; and one INITIAL row per entry of
//! the Ini stream, as the ordinary scalar expression it is (OCFEMOL::initial_expression).  The peel's remainder IS Rc, so the sign is a PLUS: the model row was  M d_t x + R.
//! The Ini entries at algebraic nodes are NOT declared -- they are the consistency check at t0 (OCFEMOL::ini_stream).
//! Mass coefficients must be CONSTANT in this first version; a state-dependent one is refused by name.
//! The references of xc are the owning field's; for a good start call OCFEMOL::init_seed() and pass its values to
//! update_ref() afterwards -- DAESLV's corrector then solves for the consistent state.
//! @return false with @p err set on a state-dependent mass coefficient or a solver that is not set up.
//! MoL workplan, phase 5.
////////////////////////////////////////////////////////////////////////
//! @brief The integrator options a semi-discrete DAE needs, set on @p S (export_semidiscrete's solver half).
template <typename SOLVER>
inline void
configure_semidiscrete
( SOLVER& S )
{
  // The consistent initialisation is DAESLV's own KINSOL corrector (INITIAL rows with the interior rows at t0).  A
  // semi-discrete system is large and sparse: the corrector's dense default would be O(n^3).  The user may override
  // after the call.  (The step bound needs nothing here: its automatic default is right for any seed.)
  S.options.INITIAL.KIN.SPARSE = true;
  // ... and the integrator's own linear solver: with the DENSE default DAESLV builds the model's full Jacobian by dense
  // forward mode, which through ONE operation with n inputs and n outputs is n^3 symbolic nodes (476 MB at 144
  // unknowns, 1.8 GB at 224 -- measured).  The sparse solver asks the operation for its sparse derivative.
  S.options.SOLVE.LINSOL = SOLVER::Options::SPARSE;
  // ... and algebraic components out of IDAS's local error test.  The corrector gives a consistent x(t0) and the
  // DIFFERENTIAL x_t(t0), but leaves an algebraic component's derivative at zero (it appears in no row).  A semi-
  // discrete system has many algebraic unknowns slaved to boundary data, and where that data moves at t0 -- a ramped
  // feed -- the first step's error estimate on them cannot pass (measured: the PSA column, h -> 4e-9 at t = 0).
  S.options.SOLVE.SUPPRESSALG = true;
}

//! @brief The semi-discrete system of a set-up @p mol, DECLARED into the model @p M on @p dag (export_semidiscrete's
//! model half): the domain @p tdom of t, one state per nodal unknown xc_i on {t} with its field's reference, the rows
//! Mc dxc/dt + Rc(xc, t) = 0 (a row with mass is differential and holds after t0; an algebraic row holds everywhere)
//! and the INITIAL rows at t0.  @p M may be any FFModel -- a solver, or a plain record whose declarations a solver
//! then sets up from (DAESLV, for a model distributed in space).
inline bool
declare_semidiscrete
( OCFEMOL& mol, FFGraph& dag, FFModel& M, FFDom const& tdom, std::vector<FFVar>& xc, FFVar& t, std::string* err = nullptr )
{
  auto fail = [&]( std::string const& m ){ if( err ) *err = m; return false; };
  auto const& D = mol.semi_discrete();
  if( !D.ran ) return fail( "the OCFEMOL is not set up" );
  for( auto const& m : mol.mass() )
    if( !m.coeff.cst() )
      return fail( "the mass coefficient of row " + std::to_string( m.row ) + " is state-dependent; this hand-over takes"
                   " constant mass coefficients only (a node-wise coefficient operation is the next step)" );

  M.set_dag( &dag );
  t = dag.add_var( "t" );
  xc.assign( D.n_states, FFVar() );
  for( size_t i = 0; i < D.n_states; ++i ) xc[i] = dag.add_var( "xc_" + std::to_string( i ) );
  M.add_domain( t, tdom );
  M.set_evolution_domain( t );
  FFVar fld;  double ref = 0.;
  for( size_t i = 0; i < D.n_states; ++i ){
    M.add_state( xc[i], {t} );
    M.update_ref( xc[i], mol.column_field( i, fld, ref )? ref: 0. );
  }

  // Rc and the Ini stream as operations over (xc, t)
  FFMOLRES RES;  std::vector<FFVar> rc  = RES( xc, t, &mol );

  // the rows: mass terms + the remainder.  A row with mass is differential and does not hold at t0 (its INITIAL
  // rows do); an algebraic row holds everywhere, t0 included.
  FFPartial OpP;
  std::vector< std::vector< std::pair<size_t,double> > > mass_of_row( D.n_rows );
  for( auto const& m : mol.mass() ) mass_of_row[ m.row ].emplace_back( m.col, m.coeff.num().val() );
  auto const INTER = FFModel::EqnRole::INTERIOR, INIT = FFModel::EqnRole::INITIAL;
  int const TI = FFDom::ALL - FFDom::LB;
  for( size_t r = 0; r < D.n_rows; ++r ){
    if( mass_of_row[r].empty() ){ M.add_equation( rc[r], {t}, {FFDom::ALL}, INTER ); continue; }
    FFVar row = rc[r];
    for( auto const& [c, a] : mass_of_row[r] ) row += a * OpP( xc[c], t );
    M.add_equation( row, {t}, {TI}, INTER );
  }
  // the INITIAL rows as ordinary pointwise expressions, not as outputs of FFMOLINI (see OCFEMOL::initial_expression)
  auto const& IS = mol.ini_stream();
  for( size_t k = 0; k < IS.size(); ++k )
    if( !IS[k].dropped ) M.add_equation( mol.initial_expression( k, xc, t ), {t}, {FFDom::LB}, INIT );
  return true;
}

//! @brief declare_semidiscrete into the solver @p S, then configure_semidiscrete( @p S ) -- the explicit two-object
//! hand-over of MoL phase 5.  DAESLV takes a model distributed in space itself (DAESLV_IDAS::setup).
template <typename SOLVER>
inline bool
export_semidiscrete
( OCFEMOL& mol, FFGraph& dag, SOLVER& S, FFDom const& tdom, std::vector<FFVar>& xc, FFVar& t, std::string* err = nullptr )
{
  if( !declare_semidiscrete( mol, dag, S, tdom, xc, t, err ) ) return false;
  configure_semidiscrete( S );
  return true;
}

} // namespace mc

#endif
