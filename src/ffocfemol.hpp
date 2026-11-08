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

//! @brief What FFMOLRES / FFGradMOLRES ask an OCFEMOL, by kind: the residual rows or the outputs' spatial parts
namespace molop
{
inline size_t rows( OCFEMOL const* m, int kind ){ return kind? m->n_outputs(): m->semi_discrete().n_rows; }
inline void pattern_x( OCFEMOL const* m, int kind, std::vector<size_t>& p, std::vector<size_t>& c )
{ if( kind ) m->output_state_pattern( p, c ); else m->jacobian_pattern( p, c ); }
inline void pattern_u( OCFEMOL const* m, int kind, std::vector<size_t>& p, std::vector<size_t>& c )
{ if( kind ) m->output_input_pattern( p, c ); else m->input_pattern( p, c ); }
inline bool eval( OCFEMOL* m, int kind, double t, double const* xc, double const* uc, double* out )
{ return kind? m->eval_outputs( t, xc, uc, out ): m->eval_residual( t, xc, uc, out ); }
inline bool deriv_full( OCFEMOL* m, int kind, double t, double const* xc, double const* uc, double* out,
                        std::vector<size_t>& p, std::vector<size_t>& c, std::vector<double>& v )
{ return kind? m->deriv_outputs_full( t, xc, uc, out, p, c, v ): m->deriv_residual_full( t, xc, uc, out, p, c, v ); }
}

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
  //! @brief What the operation's outputs are: the residual rows Rc(xc, uc), or the outputs' spatial parts F(xc, uc)
  //! (OCFEMOL::eval_outputs) -- one implementation, OCFEMOL asked by kind (see molop)
  enum Kind { RESIDUAL = 0, OUTPUTS = 1 };

protected:
  OCFEMOL*     _pMOL = nullptr;   //!< the set-up solver
  bool         _own  = false;     //!< whether this operation owns (and deletes) it
  int          _kind = RESIDUAL;  //!< see Kind
  size_t       _nSta = 0;         //!< n_states
  size_t       _nInp = 1;         //!< n_input_slots: the collocated input vector uc, t in its own slot
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
    ( size_t const nSta, FFVar const* pSta, size_t const nInp, FFVar const* pInp, OCFEMOL* pMOL, int policy, std::string const& name )
    {
      if( !pMOL || !pMOL->semi_discrete().ran || nSta != pMOL->semi_discrete().n_states )
        throw std::runtime_error( "FFMOLRES ** the solver is not set up or the state count does not match\n" );
      if( nInp != std::max<size_t>( 1, pMOL->n_input_slots() ) )
        throw std::runtime_error( "FFMOLRES ** the input vector must have n_input_slots() entries, t in time_slot()\n" );
      if( _own && _pMOL ) delete _pMOL;
      _own  = ( policy > 0 );
      _pMOL = pMOL;
      _nSta = nSta;
      _nInp = nInp;
      _nRow = molop::rows( pMOL, _kind );
      _name = name;
      data    = pMOL;
      owndata = false;
      std::vector<FFVar> vIn( pSta, pSta + nSta );  vIn.insert( vIn.end(), pInp, pInp + nInp );
      FFVar** ppRes = insert_external_operation( *this, _nRow, vIn.size(), vIn.data() );
      _own = false;
      FFOp* pOp = (*ppRes)->opdef().first;
      if( policy > 0 )      _pMOL = dynamic_cast<FFMOLRES*>( pOp )->_pMOL;
      else if( policy < 0 ) dynamic_cast<FFMOLRES*>( pOp )->_own = true;
      return ppRes;
    }

public:
  FFMOLRES
    ( Kind kind = RESIDUAL )
    : FFOp( EXTERN ), _kind( kind )
    { this->sparse = true; }   // MC++ calls the SPARSE deriv only when this is set; unset, it allocates a dense
                                // n_rows x (n_states+1) array per call and calls the dense deriv (FFOCFERES sets it too)

  virtual ~FFMOLRES
    ()
    { if( _own && _pMOL ) delete _pMOL; }

  //! @brief Insert the operation: the n_rows outputs Rc(xc, uc) for the states @p vSta (n_states) and the collocated
  //! input vector @p vInp (n_input_slots(): t in time_slot(), each other slot a variable of the integrator).
  std::vector<FFVar> operator()
    ( std::vector<FFVar> const& vSta, std::vector<FFVar> const& vInp, OCFEMOL* pMOL, int policy = SHALLOW,
      std::string const& name = "" )
    {
      FFVar** ppRes = _set( vSta.size(), vSta.data(), vInp.size(), vInp.data(), pMOL, policy, name );
      std::vector<FFVar> vDep( _nRow );
      for( size_t i = 0; i < _nRow; ++i ) vDep[i] = *ppRes[i];
      return vDep;
    }
  //! @brief For a model whose only input is t: inputs @p vSta and @p t.
  std::vector<FFVar> operator()
    ( std::vector<FFVar> const& vSta, FFVar const& t, OCFEMOL* pMOL, int policy = SHALLOW, std::string const& name = "" )
    {
      if( pMOL && pMOL->n_input_slots() > 1 )
        throw std::runtime_error( "FFMOLRES ** the model has inputs: give the collocated input vector, not t alone\n" );
      return operator()( vSta, std::vector<FFVar>{ t }, pMOL, policy, name );
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
      return _kind < static_cast<FFMOLRES const*>( op )->_kind;      // the residual and the outputs of one OCFEMOL
    }
  std::string name
    ()
    const
    { std::ostringstream oss; oss << ( _name.empty()? ( _kind? "MOLFCT": "MOLRES" ): _name ) << "[" << _nRow << "x" << _nSta + _nInp << "]"; return oss.str(); }
  bool commutative
    ()
    const
    { return false; }
};

//! @brief The Jacobian entries of FFMOLRES as a DAG external operation: the values of dRc/dxc in CSR order
//! (OCFEMOL::jacobian_pattern), then those of dRc/d(the other input slots) in CSR order (OCFEMOL::input_pattern),
//! then the n_rows values of dRc/dt.
////////////////////////////////////////////////////////////////////////
class FFGradMOLRES
: public FFOp
{
protected:
  OCFEMOL*     _pMOL = nullptr;
  bool         _own  = false;
  int          _kind = FFMOLRES::RESIDUAL;
  size_t       _nSta = 0, _nInp = 1, _nRow = 0, _nnz = 0, _nnu = 0;
  std::string  _name;

public:
  FFGradMOLRES
    ( int kind = FFMOLRES::RESIDUAL )
    : FFOp( EXTERN ), _kind( kind )
    {}
  virtual ~FFGradMOLRES
    ()
    { if( _own && _pMOL ) delete _pMOL; }

  //! @brief Insert over the inputs [xc, uc]: see the class.
  FFVar** operator()
    ( size_t const nSta, FFVar const* pSta, size_t const nInp, FFVar const* pInp, OCFEMOL* pMOL, bool own, std::string const& name )
    {
      _pMOL = pMOL; _own = false; _nSta = nSta; _nInp = nInp; _nRow = molop::rows( pMOL, _kind ); _name = name;
      std::vector<size_t> ptr, col;  molop::pattern_x( pMOL, _kind, ptr, col );  _nnz = col.size();
      molop::pattern_u( pMOL, _kind, ptr, col );  _nnu = col.size();
      data = pMOL;  owndata = false;
      std::vector<FFVar> vIn( pSta, pSta + nSta );  vIn.insert( vIn.end(), pInp, pInp + nInp );
      FFVar** ppRes = insert_external_operation( *this, _nnz + _nnu + _nRow, vIn.size(), vIn.data() );
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
      if( !_pMOL || nVar != _nSta + _nInp || nRes != _nnz + _nnu + _nRow )
        throw std::runtime_error( "FFGradMOLRES::eval double ** size mismatch\n" );
      size_t const it = _pMOL->time_slot(), tc = _nSta + it;
      std::vector<size_t> ptr, col; std::vector<double> val;
      if( !molop::deriv_full( _pMOL, _kind,  vVar[tc], vVar, vVar + _nSta, nullptr, ptr, col, val ) )
        throw std::runtime_error( "FFGradMOLRES::eval double ** OCFEMOL::deriv_residual_full failed\n" );
      // split each row's entries: xc first (jacobian_pattern's order), the other input slots next, t last (dense)
      size_t kx = 0, ku = _nnz;
      for( size_t r = 0; r + 1 < ptr.size(); ++r ){
        vRes[_nnz + _nnu + r] = 0.;
        for( size_t j = ptr[r]; j < ptr[r+1]; ++j ){
          if( col[j] < _nSta )  vRes[kx++] = val[j];
          else if( col[j] == tc ) vRes[_nnz + _nnu + r] = val[j];
          else                  vRes[ku++] = val[j];
        }
      }
      if( kx != _nnz || ku != _nnz + _nnu ) throw std::runtime_error( "FFGradMOLRES::eval double ** pattern changed\n" );
    }
  void eval
    ( unsigned const nRes, FFVar* vRes, unsigned const nVar, FFVar const* vVar, unsigned const* mVar )
    const
    {
      (void)mVar;
      FFGradMOLRES op( _kind );
      FFVar** ppRes = op( _nSta, vVar, _nInp, vVar + _nSta, _pMOL, false, _name );
      for( unsigned k = 0; k < nRes; ++k ) vRes[k] = *ppRes[k];
    }
  template <typename U>
  void eval
    ( unsigned const nRes, U* vRes, unsigned const nVar, U const* vVar, unsigned const* mVar )
    const
    { throw std::runtime_error( "FFGradMOLRES::eval ** evaluation type not supported\n" ); }

  bool lt( FFOp const* op ) const { if( data != op->data ) return data < op->data; return _kind < static_cast<FFGradMOLRES const*>( op )->_kind; }
  std::string name() const { std::ostringstream oss; oss << ( _name.empty()? ( _kind? "GradMOLFCT": "GradMOLRES" ): _name ) << "[" << _nnz + _nnu + _nRow << "]"; return oss.str(); }
  bool commutative() const { return false; }
};


// (FFMOLINI, the Ini stream as one vectorised operation, and FFGradMOLINI were retired on 5 Nov 2026: the integrator
//  PEELS INITIAL rows one at a time, composing each with its target state, and every such compose rebuilt the whole
//  operation -- the INITIAL rows are declared as pointwise expressions (OCFEMOL::initial_expression), and so are the
//  transitions (transition_expressions).  FFMOLRES keeps the kinds whose vectors the integrator consumes whole.)

inline FFVar const*
FFMOLRES::_gradient
( FFVar const* vVar )
const
{
  std::vector<long> ids( _nSta + _nInp );  for( size_t i = 0; i < _nSta + _nInp; ++i ) ids[i] = vVar[i].id().second;
  if( _gradOut.empty() || ids != _gradIn ){
    FFGradMOLRES gop( _kind );
    FFVar** ppG = gop( _nSta, vVar, _nInp, vVar + _nSta, _pMOL, false, _name );
    std::vector<size_t> ptr, col, pu, cu;  molop::pattern_x( _pMOL, _kind, ptr, col );  molop::pattern_u( _pMOL, _kind, pu, cu );
    _gradOut.assign( col.size() + cu.size() + _nRow, FFVar() );
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
  if( !_pMOL || nVar != _nSta + _nInp || nRes != _nRow ) throw std::runtime_error( "FFMOLRES::eval double ** size mismatch\n" );
  if( !molop::eval( _pMOL, _kind,  vVar[_nSta + _pMOL->time_slot()], vVar, vVar + _nSta, vRes ) )
    throw std::runtime_error( "FFMOLRES::eval double ** OCFEMOL::eval_residual failed\n" );
}

inline void
FFMOLRES::eval
( unsigned const nRes, FADType<double>* vRes, unsigned const nVar, FADType<double> const* vVar, unsigned const* mVar )
const
{
  (void)mVar;
  if( !_pMOL || nVar != _nSta + _nInp || nRes != _nRow ) throw std::runtime_error( "FFMOLRES::eval FAD ** size mismatch\n" );
  size_t const tc = _nSta + _pMOL->time_slot();
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
    if( !molop::eval( _pMOL, _kind,  x[tc], x.data(), x.data() + _nSta, r.data() ) ) throw std::runtime_error( "FFMOLRES::eval FAD ** eval failed\n" );
    for( unsigned k = 0; k < nRes; ++k ) vRes[k] = r[k];
    return;
  }
  // one compressed Jacobian over [xc | uc], then each output's directions as a sparse combination of the inputs' --
  // xc's entries, the other input slots', and t's last (its column dense, as the operation declares it)
  std::vector<size_t> ptr, col; std::vector<double> val;
  if( !molop::deriv_full( _pMOL, _kind,  x[tc], x.data(), x.data() + _nSta, r.data(), ptr, col, val ) )
    throw std::runtime_error( "FFMOLRES::eval FAD ** OCFEMOL::deriv_residual_full failed\n" );
  for( unsigned k = 0; k < nRes; ++k ){
    vRes[k] = r[k];
    vRes[k].diff( 0, nDir );  vRes[k].d( 0 ) = 0.;      // allocate nDir directions
    double drdt = 0.;
    for( size_t j = ptr[k]; j < ptr[k+1]; ++j ) if( col[j] == tc ) drdt = val[j];
    for( unsigned d = 0; d < nDir; ++d ){
      double s = 0.;
      for( size_t j = ptr[k]; j < ptr[k+1]; ++j ) if( col[j] != tc && vVar[col[j]].depend() ) s += val[j] * vVar[col[j]].deriv( d );
      if( vVar[tc].depend() ) s += drdt * vVar[tc].deriv( d );
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
  if( !_pMOL || nVar != _nSta + _nInp || nRes != _nRow ) throw std::runtime_error( "FFMOLRES::eval FFDep ** size mismatch\n" );
  std::vector<size_t> ptr, col, pu, cu;  molop::pattern_x( _pMOL, _kind, ptr, col );  molop::pattern_u( _pMOL, _kind, pu, cu );
  for( unsigned k = 0; k < nRes; ++k ){
    vRes[k] = FFDep();
    for( size_t j = ptr[k]; j < ptr[k+1]; ++j ) vRes[k] += vVar[col[j]];
    for( size_t j = pu[k]; j < pu[k+1]; ++j ) vRes[k] += vVar[ _nSta + cu[j] ];
    vRes[k] += vVar[ _nSta + _pMOL->time_slot() ];        // every row may depend on t (conservative)
    vRes[k].update( FFDep::TYPE::N );                     // nonlinear in general
  }
}

inline void
FFMOLRES::eval
( unsigned const nRes, FFVar* vRes, unsigned const nVar, FFVar const* vVar, unsigned const* mVar )
const
{
  (void)mVar;
  FFMOLRES op( static_cast<Kind>( _kind ) );
  FFVar** ppRes = op._set( _nSta, vVar, _nInp, vVar + _nSta, _pMOL, SHALLOW, _name );
  for( unsigned k = 0; k < nRes; ++k ) vRes[k] = *ppRes[k];
}

inline void
FFMOLRES::deriv
( unsigned const nRes, FFVar const* vRes, unsigned const nVar, FFVar const* vVar, FFVar** vDer )
const
{
  // dense: every derivative, from the sparse entries
  std::vector<size_t> ptr, col, pu, cu;  molop::pattern_x( _pMOL, _kind, ptr, col );  molop::pattern_u( _pMOL, _kind, pu, cu );
  size_t const tc = _nSta + _pMOL->time_slot();
  FFVar const* ppG = _gradient( vVar );                 // one gradient operation per instance, reused
  for( unsigned k = 0; k < nRes; ++k ){
    for( unsigned i = 0; i < nVar; ++i ) vDer[k][i] = 0.;
    for( size_t j = ptr[k]; j < ptr[k+1]; ++j ) vDer[k][ col[j] ] = ppG[j];
    for( size_t j = pu[k]; j < pu[k+1]; ++j ) vDer[k][ _nSta + cu[j] ] = ppG[ col.size() + j ];
    vDer[k][tc] = ppG[ col.size() + cu.size() + k ];
  }
}

inline void
FFMOLRES::deriv
( unsigned const nRes, FFVar const* vRes, unsigned const nVar, FFVar const* vVar, FFVar** vDer,
  size_t* nnz, size_t** colnz )
const
{
  std::vector<size_t> ptr, col, pu, cu;  molop::pattern_x( _pMOL, _kind, ptr, col );  molop::pattern_u( _pMOL, _kind, pu, cu );
  size_t const tc = _nSta + _pMOL->time_slot();
  // each row: its xc columns, its other input slots, and t's column -- the last always (declared conservatively)
  for( unsigned k = 0; k < nRes; ++k ) nnz[k] = ( ptr[k+1] - ptr[k] ) + ( pu[k+1] - pu[k] ) + 1;
  if( colnz ){
    for( unsigned k = 0; k < nRes; ++k ){
      if( !colnz[k] ) colnz[k] = new size_t[ nnz[k] ];
      size_t m = 0;
      for( size_t j = ptr[k]; j < ptr[k+1]; ++j ) colnz[k][m++] = col[j];
      for( size_t j = pu[k]; j < pu[k+1]; ++j ) colnz[k][m++] = _nSta + cu[j];
      colnz[k][m] = tc;
    }
  }
  if( !vDer ) return;
  FFVar const* ppG = _gradient( vVar );                 // one gradient operation per instance, reused
  for( unsigned k = 0; k < nRes; ++k ){
    if( !vDer[k] ) vDer[k] = new FFVar[ nnz[k] ];
    size_t m = 0;
    for( size_t j = ptr[k]; j < ptr[k+1]; ++j ) vDer[k][m++] = ppG[j];
    for( size_t j = pu[k]; j < pu[k+1]; ++j ) vDer[k][m++] = ppG[ col.size() + j ];
    vDer[k][m] = ppG[ col.size() + cu.size() + k ];
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
//! The model's inputs: each slot of a time-invariant input becomes an input of @p M with no domain -- a parameter of the
//! integrator -- in input_dofs() order; @p uc, if given, receives the collocated input vector (t in time_slot(), the
//! parameter of each other slot; a time-varying input's slots are inputs over t, with its declared collocation).
inline bool
declare_semidiscrete
( OCFEMOL& mol, FFGraph& dag, FFModel& M, FFDom const& tdom, std::vector<FFVar>& xc, FFVar& t, std::string* err = nullptr,
  std::vector<FFVar>* uc = nullptr )
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

  // the collocated input vector: t in its slot; each slot of a time-invariant input a PARAMETER of the integrator,
  // each slot of a time-varying one an INPUT OVER t with the declared collocation in t (its levels per element are
  // the integrator's parameters; within a stage it is the interpolant of them -- BASE_DE's own machinery)
  std::vector<FFVar> U( std::max<size_t>( 1, mol.n_input_slots() ), FFVar( 0. ) );
  U[ mol.time_slot() ] = t;
  for( auto const& L : mol.inputs_layout() ){
    FFDom::TYPE ttype = FFDom::LGR;  size_t tnode = 0;
    bool const tcoll = L.time_varying && mol.input_time_collocation( L.input, ttype, tnode );
    for( size_t k = 0; k < L.cols.size(); ++k ){
      FFVar const p = dag.add_var( L.input.name() + "_" + std::to_string( k ) );
      if( !L.time_varying ) M.add_input( p, {} );
      else if( tcoll )      M.add_input( p, {t}, ttype, tnode );
      else                  M.add_input( p, {t} );
      U.at( L.cols[k] ) = p;
    }
  }
  if( uc ) *uc = U;

  // Rc as one operation over (xc, uc)
  FFMOLRES RES;  std::vector<FFVar> rc  = RES( xc, U, &mol );

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
    if( !IS[k].dropped ) M.add_equation( mol.initial_expression( k, xc, t, U ), {t}, {FFDom::LB}, INIT );
  // the transitions (MoL 6c): each declared one applied at its state's differential nodes, both sides composed there;
  // the integrator's stop machinery does the rest (the algebraic nodes re-solved, the sensitivities carried across)
  for( size_t k = 0; k < mol.n_transitions(); ++k ){
    std::vector<FFVar> lft, rgt;  std::string why;
    if( !mol.transition_expressions( k, xc, t, lft, rgt, why, U ) ) return fail( "transition " + std::to_string( k ) + ": " + why );
    M.add_transition( lft, rgt, t, mol.transition_time( k ) );
  }
  // the outputs: their spatial parts as ONE operation over (xc, uc), each value given its time part back -- a point in
  // t, or an integral over t -- in the declared order, so the integrator's outputs are the model's (val_function())
  if( mol.n_outputs() ){
    FFMOLRES FCT( FFMOLRES::OUTPUTS );
    std::vector<FFVar> const fc = FCT( xc, U, &mol );
    FFEval OpE;  FFIntegral OpI;
    auto const& T = mol.outputs_time();
    auto spatial = [&]( size_t k )->FFVar{                 // the scalar spatial part of output k, at every t
      if( T[k].sparts.empty() ) return fc.at( mol.output_rows( k ).first );
      std::vector<FFVar> vars( T[k].sph ), vals;         // hs over the nested reductions' rows, the constants by value
      for( size_t i : T[k].sparts ) vals.push_back( fc.at( mol.output_rows( i ).first ) );
      for( size_t q = 0; q < mol.constant_values().size(); ++q ){ vars.push_back( mol.constant_vars()[q] );  vals.push_back( FFVar( mol.constant_values()[q] ) ); }
      std::vector<FFVar> out( 1 );
      mol.dag()->eval( std::vector<FFVar>{ T[k].hs }, out, vars, vals );
      return out[0]; };
    auto reduction = [&]( size_t k )->FFVar{               // the reduction in t of the scalar spatial part of output k
      return T[k].kind == OCFEMOL::t_OutputTime::INTEGRAL? OpI( spatial( k ), t ): OpE( spatial( k ), t, T[k].tau, T[k].side ); };
    size_t const ndecl = mol.n_declared_outputs();          // the parts of composed outputs follow the declared ones
    for( size_t k = 0; k < ndecl && k < T.size(); ++k ){
      if( T[k].kind == OCFEMOL::t_OutputTime::COMPOSED ){
        // h over the parts' reductions: the deferred inputs replaced by them, the constants by their values
        std::vector<FFVar> vars( T[k].dvars ), vals;
        for( size_t i : T[k].parts ){
          if( mol.output_rows( i ).second != 1 ) return fail( "a composed output's part is not a scalar in space" );
          vals.push_back( reduction( i ) );
        }
        for( size_t q = 0; q < mol.constant_values().size(); ++q ){ vars.push_back( mol.constant_vars()[q] );  vals.push_back( FFVar( mol.constant_values()[q] ) ); }
        std::vector<FFVar> out( 1 );
        mol.dag()->eval( std::vector<FFVar>{ T[k].h }, out, vars, vals );
        M.add_output( out[0] );
        continue;
      }
      if( !T[k].sparts.empty() ){ M.add_output( reduction( k ) );  continue; }
      if( T[k].kind == OCFEMOL::t_OutputTime::TRAJECTORY ){          // one trajectory per spatial row: space-major
        auto const [r0, nr] = mol.output_rows( k );
        for( size_t r = r0; r < r0 + nr; ++r ) M.add_output( fc.at( r ), { t }, std::vector<int>{ T[k].tmask } );
        continue;
      }
      auto const [r0, nr] = mol.output_rows( k );
      for( size_t r = r0; r < r0 + nr; ++r )
        M.add_output( T[k].kind == OCFEMOL::t_OutputTime::INTEGRAL? OpI( fc.at( r ), t ): OpE( fc.at( r ), t, T[k].tau, T[k].side ) );
    }
  }
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
