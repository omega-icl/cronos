// Copyright (C) Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

#ifndef CRONOS__FFODE_HPP
#define CRONOS__FFODE_HPP

#include "ffexpr.hpp"
#include "ffdep.hpp"
#include "slift.hpp"
#include "odeslvs_cvodes.hpp"

namespace mc
{

//! @brief C++ class defining IVP in parametric ODEs as external DAG operations in MC++.
////////////////////////////////////////////////////////////////////////
//! mc::FFBASEODE is a C++ base class for defining the options in an
//! IVP in parametric ODEs as external DAG operations in MC++.
////////////////////////////////////////////////////////////////////////
class FFBaseODE
: public FFOp
{

protected:

  // Pointer to ODE solver
  ODESLVS_CVODES*     _pODESLV;
  // Whether this class owns _pODESLV
  bool                _ownODESLV;
  // Number of parameters
  size_t              _nPar;
  // Number of constants
  size_t              _nCst;

public:

  // Default constructor
  FFBaseODE
    ()
    : FFOp( EXTERN ),
      _pODESLV( nullptr ),
      _ownODESLV( false )
    {}

  // Destructor
  virtual ~FFBaseODE
    ()
    {
#ifdef CRONOS__FFODE_TRACE
      std::cout << "FFBaseODE::destructor\n";
#endif
      if( _ownODESLV && _pODESLV )
        delete _pODESLV;
    }

  // Copy constructor
  FFBaseODE
    ( FFBaseODE const& Op )
    : FFOp( Op ),
      _nPar( Op._nPar ),
      _nCst( Op._nCst )
    {
#ifdef CRONOS__FFODE_TRACE
      std::cout << "FFBaseODE::copy constructor\n";
#endif
      if( !Op._pODESLV )
        throw std::runtime_error( "FFBaseODE::copy constructor ** Undefined ODE solver\n" );

      _ownODESLV = Op._ownODESLV;      
      if( _ownODESLV ){
        _pODESLV = new ODESLVS_CVODES;
        _pODESLV->set( *Op._pODESLV );
        _pODESLV->options = Op._pODESLV->options;
        _pODESLV->setup( *Op._pODESLV ); // Create IVP copy from local copy in Op._pODESLV
#ifdef CRONOS__FFODE_TRACE
        std::cerr << "ODE address copied: " << _pODESLV << std::endl;
#endif
      }
      else
        _pODESLV   = Op._pODESLV;
    }

  //! @brief MCODE options
  static struct Options
  {
    //! @brief Constructor
    Options():
      DIFF(NUM_P), NP2NF(3)
      {}
    //! @brief Assignment operator
    Options& operator= ( Options const& options ){
        DIFF   = options.DIFF;
        return *this;
      }
    //! @brief Enumeration type for ODE differentiation
    enum DERIV_TYPE{
      NUM_P=0,  //!< Derivatives w.r.t. parameters only through forward or adjoint sensitivity integration
      SYM_P,    //!< Derivatives w.r.t. parameters only through symbolic differentiation of ODEs
      SYM_C,    //!< Derivatives w.r.t. constants only through symbolic differentiation of ODEs
      SYM_PC,   //!< Derivatives w.r.t. parameters and constants jointly through symbolic differentiation of ODEs
    };
    //! @brief Selected ODE differentiation
    DERIV_TYPE                DIFF;
    //! @brief parameter-to-function-size ratio above which adjoint sensitivity is applied instead of forward sensitivity
    double                    NP2NF;
  } options;
};

inline FFBaseODE::Options FFBaseODE::options;

//! @brief C++ class defining IVP in parametric ODEs as external DAG operations in MC++.
////////////////////////////////////////////////////////////////////////
//! mc::FFODE is a C++ class for defining an IVP in parametric ODEs as
//! external DAG operations in MC++.
////////////////////////////////////////////////////////////////////////
class FFODE
: public FFBaseODE
{

public:

  //! @brief Enumeration type for ODE copy policy
  enum POLICY_TYPE{
    SHALLOW=0,  //!< Shallow copy of ODE system in FFGraph (without ownership)
    COPY=1,     //!< Deep copy of ODE system in FFGraph (with ownership)
    TRANSFER=-1 //!< Shallow copy of ODE system in FFGraph (with ownership transfer)
  };

  // Default constructor
  FFODE
    ()
    : FFBaseODE()
    {}

  // Destructor
  virtual ~FFODE
    ()
    {}

  // Define operation
  FFVar& operator()
    ( unsigned const idep, unsigned const nPar, FFVar const* pPar, ODESLVS_CVODES* pODESLV, int policy=COPY )
//    const
    {
#ifdef CRONOS__FFODE_CHECK
      assert( idep < pODESLV->nf()*nPar );
#endif
      return *(_set( nPar, pPar, 0, nullptr, pODESLV, policy )[idep]);
    }

  FFVar** operator()
    ( unsigned const nPar, FFVar const* pPar, ODESLVS_CVODES* pODESLV, int policy=COPY )
//    const
    {
      return _set( nPar, pPar, 0, nullptr, pODESLV, policy );
    }

  FFVar& operator()
    ( unsigned const idep, unsigned const nPar, FFVar const* pPar, unsigned const nCst, FFVar const* pCst,
      ODESLVS_CVODES* pODESLV, int policy=COPY )
//    const
    {
#ifdef CRONOS__FFODE_CHECK
      assert( idep < pODESLV->nf()*nPar );
#endif
      return *(_set( nPar, pPar, nCst, pCst, pODESLV, policy )[idep]);
    }

  FFVar** operator()
    ( unsigned const nPar, FFVar const* pPar, unsigned const nCst, FFVar const* pCst,
      ODESLVS_CVODES* pODESLV, int policy=COPY )
//    const
    {
      return _set( nPar, pPar, nCst, pCst, pODESLV, policy );
    }

  ODESLVS_CVODES* pODESLV
    ()
    const
    { //std::cerr << "ODE address retreived: " << _pODESLV << std::endl;
      return _pODESLV; }

  // Evaluation overloads
  virtual void feval
    ( std::type_info const& idU, unsigned const nRes, void* vRes, unsigned const nVar,
      void const* vVar, unsigned const* mVar )
    const
    {
      if( idU == typeid( FFVar ) )
        return eval( nRes, static_cast<FFVar*>(vRes), nVar, static_cast<FFVar const*>(vVar), mVar );
      else if( idU == typeid( fadbad::F<FFVar> ) )
        return eval( nRes, static_cast<fadbad::F<FFVar>*>(vRes), nVar, static_cast<fadbad::F<FFVar> const*>(vVar), mVar );
      else if( idU == typeid( FFDep ) )
        return eval( nRes, static_cast<FFDep*>(vRes), nVar, static_cast<FFDep const*>(vVar), mVar );
      else if( idU == typeid( double ) )
        return eval( nRes, static_cast<double*>(vRes), nVar, static_cast<double const*>(vVar), mVar );
      else if( idU == typeid( fadbad::F<double> ) )
        return eval( nRes, static_cast<fadbad::F<double>*>(vRes), nVar, static_cast<fadbad::F<double> const*>(vVar), mVar );
      else if( idU == typeid( SLiftVar ) )
        return eval( nRes, static_cast<SLiftVar*>(vRes), nVar, static_cast<SLiftVar const*>(vVar), mVar );
//      else if( idU == typeid( FFExpr ) )
//        return eval( nRes, static_cast<FFExpr*>(vRes), nVar, static_cast<FFExpr const*>(vVar), mVar );

      throw std::runtime_error( "FFODE::feval ** No evaluation method for type"+std::string(idU.name())+"\n" );
    }

  void eval
    ( unsigned const nRes, double* vRes, unsigned const nVar, double const* vVar, unsigned const* mVar )
    const;

  void eval
    ( unsigned const nRes, fadbad::F<double>* vRes, unsigned const nVar, fadbad::F<double> const* vVar,
      unsigned const* mVar )
    const;

  void eval
    ( unsigned const nRes, FFVar* vRes, unsigned const nVar, FFVar const* vVar, unsigned const* mVar )
    const;

  void eval
    ( unsigned const nRes, FFDep* vRes, unsigned const nVar, FFDep const* vVar, unsigned const* mVar )
    const;

  void eval
    ( unsigned const nRes, fadbad::F<FFVar>* vRes, unsigned const nVar, fadbad::F<FFVar> const* vVar,
      unsigned const* mVar )
    const;

  void eval
    ( unsigned const nRes, SLiftVar* vRes, unsigned const nVar, SLiftVar const* vVar, unsigned const* mVar )
    const;

  // Derivatives
  void deriv
    ( unsigned const nRes, FFVar const* vRes, unsigned const nVar, FFVar const* vVar, FFVar** vDer )
    const;

//  // Ordering
//  bool lt
//    ( FFOp const* op )
//    const;

  // Properties
  std::string name
    ()
    const
    //{ std::ostringstream oss; oss << data; return "ODE[" + oss.str() + "]"; }
    { std::ostringstream oss; oss << _pODESLV; return "ODE[" + oss.str() + "]"; }

  //! @brief Return whether or not operation is commutative
  bool commutative
    ()
    const
    { return false; }

protected:

  FFVar** _set
    ( unsigned const nPar, FFVar const* pPar, unsigned const nCst, FFVar const* pCst,
      ODESLVS_CVODES* pODESLV, int policy )
//    const
    {
#ifdef CRONOS__FFODE_CHECK
      assert( pODESLV && nPar == _pODESLV->np() && ( !pCst || nCst == _pODESLV->nc() ) );
#endif
      if( _ownODESLV && _pODESLV )
        delete _pODESLV;
      _ownODESLV = ( policy>0? true: false ); //copy;
      _pODESLV = pODESLV;
      _nPar = nPar;
      _nCst = nCst;

      data = pODESLV;
      owndata = false;
      size_t const nFun = pODESLV->nf();
      FFVar** ppRes = ( pCst? insert_external_operation( *this, nFun, nPar, pPar, nCst, pCst ):
                              insert_external_operation( *this, nFun, nPar, pPar ) );

      _ownODESLV = false;
      FFOp* pOp = (*ppRes)->opdef().first;
      if( policy > 0 )
        _pODESLV = dynamic_cast<FFODE*>(pOp)->_pODESLV; // set pointer to DAG copy
      else if( policy < 0 )
        dynamic_cast<FFODE*>(pOp)->_ownODESLV = true;   // transfer ownership
#ifdef CRONOS__FFODE_TRACE
      std::cerr << "ODE operation address: " << this << std::endl;
      std::cerr << "ODE address in DAG: " << _pODESLV << std::endl;
#endif
      return ppRes;
    }
};

//! @brief C++ class defining gradient of IVP in parametric ODEs as external DAG operations in MC++.
////////////////////////////////////////////////////////////////////////
//! mc::FFGRADODE is a C++ class for defining the gradient of an IVP
//! in parametric ODEs as external DAG operations in MC++.
////////////////////////////////////////////////////////////////////////
class FFGRADODE
: public FFBaseODE
{
  friend class FFODE;

public:

  //! @brief Enumeration type for ODE copy policy
  enum POLICY_TYPE{
    SHALLOW=0,  //!< Shallow copy of ODE system in FFGraph (without ownership)
    COPY=1      //!< Deep copy of ODE system in FFGraph (with ownership)
  };

  // Constructors
  FFGRADODE
    ()
    : FFBaseODE()
    {}

  // Destructor
  virtual ~FFGRADODE
    ()
    {}

  // Define operation
  FFVar& operator()
    ( unsigned const idep, unsigned const nPar, FFVar const* pPar, ODESLVS_CVODES* pODESLV, int policy=COPY )
//    const
    {
#ifdef CRONOS__FFODE_CHECK
      assert( idep < pODESLV->nf()*nPar );
#endif
      return *(_set( nPar, pPar, 0, nullptr, pODESLV, policy )[idep]);
    }

  FFVar** operator()
    ( unsigned const nPar, FFVar const* pPar, ODESLVS_CVODES* pODESLV, int policy=COPY )
//    const
    {
      return _set( nPar, pPar, 0, nullptr, pODESLV, policy );
    }

  FFVar& operator()
    ( unsigned const idep, unsigned const nPar, FFVar const* pPar, unsigned const nCst, FFVar const* pCst,
      ODESLVS_CVODES* pODESLV, int policy=COPY )
//    const
    {
#ifdef CRONOS__FFODE_CHECK
      assert( idep < pODESLV->nf()*nPar );
#endif
      return *(_set( nPar, pPar, nCst, pCst, pODESLV, policy )[idep]);
    }

  FFVar** operator()
    ( unsigned const nPar, FFVar const* pPar, unsigned const nCst, FFVar const* pCst,
      ODESLVS_CVODES* pODESLV, int policy=COPY )
//    const
    {
      return _set( nPar, pPar, nCst, pCst, pODESLV, policy );
    }

  // Evaluation overloads
  virtual void feval
    ( std::type_info const& idU, unsigned const nRes, void* vRes, unsigned const nVar,
      void const* vVar, unsigned const* mVar )
    const
    {
      if( idU == typeid( FFVar ) )
        return eval( nRes, static_cast<FFVar*>(vRes), nVar, static_cast<FFVar const*>(vVar), mVar );
      else if( idU == typeid( FFDep ) )
        return eval( nRes, static_cast<FFDep*>(vRes), nVar, static_cast<FFDep const*>(vVar), mVar );
      else if( idU == typeid( double ) )
        return eval( nRes, static_cast<double*>(vRes), nVar, static_cast<double const*>(vVar), mVar );
      else if( idU == typeid( SLiftVar ) )
        return eval( nRes, static_cast<SLiftVar*>(vRes), nVar, static_cast<SLiftVar const*>(vVar), mVar );
//      else if( idU == typeid( FFExpr ) )
//        return eval( nRes, static_cast<FFExpr*>(vRes), nVar, static_cast<FFExpr const*>(vVar), mVar );

      throw std::runtime_error( "FFGRADODE::feval ** No evaluation method for type"+std::string(idU.name())+"\n" );
    }

  void eval
    ( unsigned const nRes, double* vRes, unsigned const nVar, double const* vVar, unsigned const* mVar )
    const;

  void eval
    ( unsigned const nRes, FFVar* vRes, unsigned const nVar, FFVar const* vVar, unsigned const* mVar )
    const;

  void eval
    ( unsigned const nRes, FFDep* vRes, unsigned const nVar, FFDep const* vVar, unsigned const* mVar )
    const;

  void eval
    ( unsigned const nRes, SLiftVar* vRes, unsigned const nVar, SLiftVar const* vVar, unsigned const* mVar )
    const;

//  // Ordering
//  bool lt
//    ( FFOp const* op )
//    const;

  // Properties
  std::string name
    ()
    const
    //{ std::ostringstream oss; oss << data; return "GRADODE[" + oss.str() + "]"; }
    { std::ostringstream oss; oss << _pODESLV; return "GRADODE[" + oss.str() + "]"; }

  //! @brief Return whether or not operation is commutative
  bool commutative
    ()
    const
    { return false; }

protected:

  FFVar** _set
    ( unsigned const nPar, FFVar const* pPar, unsigned const nCst, FFVar const* pCst,
      ODESLVS_CVODES* pODESLV, int policy )
//    const
    {
#ifdef CRONOS__FFODE_CHECK
  assert( pODESLV && nPar == _pODESLV->np() && ( !pCst || nCst == _pODESLV->nc() ) );
#endif
      if( _ownODESLV && _pODESLV )
        delete _pODESLV;
      _ownODESLV = ( policy!=0? true: false );
      _pODESLV = pODESLV;
      _nPar = nPar;
      _nCst = nCst;

      data = pODESLV;
      owndata = false;
      size_t const nFun = pODESLV->nf();
      FFVar** ppRes = ( pCst? insert_external_operation( *this, nFun*nPar, nPar, pPar, nCst, pCst ):
                              insert_external_operation( *this, nFun*nPar, nPar, pPar ) );

      _ownODESLV = false;
      FFOp* pOp = (*ppRes)->opdef().first;
      if( policy > 0 )
        _pODESLV = dynamic_cast<FFGRADODE*>(pOp)->_pODESLV; // set pointer to DAG copy
#ifdef CRONOS__FFODE_TRACE
      std::cerr << "GRADODE address in DAG: " << _pODESLV << std::endl;
#endif
      return ppRes;
    }

};

inline void
FFODE::eval
( unsigned const nRes, FFVar* vRes, unsigned const nVar, FFVar const* vVar,
  unsigned const* mVar )
const
{
#ifdef CRONOS__FFODE_TRACE
  std::cout << "FFODE::eval: FFVar\n";
  std::cerr << "ODE operation address: " << this << std::endl;
  std::cerr << "ODE address in DAG: " << _pODESLV << std::endl;
#endif
#ifdef CRONOS__FFODE_CHECK
  assert( _pODESLV && nRes == _pODESLV->nf() && nVar == _nPar+_nCst );
#endif

  FFVar** ppRes = ( _nCst? insert_external_operation( *this, nRes, _nPar, vVar, _nCst, vVar+_nPar ):
                           insert_external_operation( *this, nRes, _nPar, vVar ) );
  //FFVar** ppRes = _set( _nPar, vVar, _nCst, vVar+_nPar, static_cast<ODESLVS_CVODES*>(data), _ownODESLV );
  for( unsigned j=0; j<nRes; ++j )
    vRes[j] = *(ppRes[j]);
}

inline void
FFODE::eval
( unsigned const nRes, FFDep* vRes, unsigned const nVar, FFDep const* vVar,
  unsigned const* mVar )
const
{
#ifdef CRONOS__FFODE_TRACE
  std::cout << "FFODE::eval: FFDep\n";
#endif
#ifdef CRONOS__FFODE_CHECK
  assert( _pODESLV && nRes == _pODESLV->nf() && nVar == _nPar+_nCst );
#endif

  vRes[0] = 0;
  for( unsigned i=0; i<nVar; ++i ) vRes[0] += vVar[i];
  vRes[0].update( FFDep::TYPE::N );
  for( unsigned j=1; j<nRes; ++j ) vRes[j] = vRes[0];
}

inline void
FFODE::eval
( unsigned const nRes, double* vRes, unsigned const nVar, double const* vVar,
  unsigned const* mVar )
const
{
#ifdef CRONOS__FFODE_TRACE
  std::cerr << "FFODE::eval: double w/ ODE address " << _pODESLV << std::endl;
  for( unsigned i=0; i<nVar; ++i ) std::cout << "vVar[" << i << "] = " << vVar[i] << std::endl;
#endif
#ifdef CRONOS__FFODE_CHECK
  assert( _pODESLV && nRes == _pODESLV->nf() && nVar == _nPar+_nCst );
#endif

  if( _pODESLV->solve_state( vVar, _nCst? vVar+_nPar: nullptr ) != ODESLVS_CVODES::NORMAL ){
    //for( unsigned i=0; i<nVar; ++i ) std::cout << "vVar[" << i << "] = " << vVar[i] << std::endl;
    throw std::runtime_error( "FFODE::eval double ** State integration failure\n" );
  }

  for( unsigned i=0; i<nRes; ++i )
    vRes[i] = _pODESLV->val_function()[i];
}

inline void
FFODE::eval
( unsigned const nRes, fadbad::F<double>* vRes, unsigned const nVar, fadbad::F<double> const* vVar,
  unsigned const* mVar )
const
{
#ifdef CRONOS__FFODE_TRACE
  std::cout << "FFODE::eval: fadbad::F<double>\n";
#endif
#ifdef CRONOS__FFODE_CHECK
  assert( _pODESLV && nRes == _pODESLV->nf() && nVar == _nPar+_nCst );
#endif

  std::vector<double> vVarVal( nVar );
  for( unsigned i=0; i<nVar; ++i ) vVarVal[i] = vVar[i].val();
  if( _nPar <= options.NP2NF * nRes ){
    if( _pODESLV->solve_sensitivity( vVarVal.data(), _nCst? vVarVal.data()+_nPar: nullptr ) != ODESLVS_CVODES::NORMAL )
      throw std::runtime_error( "FFODE::eval fadbad::F<double> ** Forward sensitivity integration failure\n" );
  }
  else{
    if( _pODESLV->solve_adjoint( vVarVal.data(), _nCst? vVarVal.data()+_nPar: nullptr ) != ODESLVS_CVODES::NORMAL )
      throw std::runtime_error( "FFODE::eval fadbad::F<double> ** Adjoint sensitivity integration failure\n" );
  }

  for( unsigned k=0; k<nRes; ++k ){
    vRes[k] = _pODESLV->val_function()[k];
    for( unsigned i=0; i<_nPar; ++i )
      vRes[k].setDepend( vVar[i] );
    for( unsigned j=0; j<vRes[k].size(); ++j ){
      vRes[k][j] = 0.;
      for( unsigned i=0; i<_nPar; ++i ){
        if( vVar[i][j] == 0. ) continue;
        vRes[k][j] += _pODESLV->val_function_gradient()[i][k] * vVar[i][j];
      }
    }
  }
}

inline void
FFODE::eval
( unsigned const nRes, fadbad::F<FFVar>* vRes, unsigned const nVar, fadbad::F<FFVar> const* vVar,
  unsigned const* mVar )
const
{
#ifdef CRONOS__FFODE_TRACE
  std::cout << "FFODE::eval: fadbad::F<FFVar>\n";
#endif
#ifdef CRONOS__FFODE_CHECK
  assert( _pODESLV && nRes == _pODESLV->nf() && nVar >= _pODESLV->np() );
#endif

  std::vector<FFVar> vVarVal( nVar );
  for( unsigned i=0; i<nVar; ++i )
    vVarVal[i] = vVar[i].val();
  FFVar const*const* vResVal = ( _nCst? insert_external_operation( *this, nRes, _nPar, vVarVal.data(), _nCst, vVarVal.data()+_nPar ):
                                        insert_external_operation( *this, nRes, _nPar, vVarVal.data() ) );
  //FFVar const*const* vResVal = _set( _nPar, vVarVal.data(), _nCst, vVarVal.data()+_nPar, static_cast<ODESLVS_CVODES*>(data), _ownODESLV );

  if( options.DIFF == Options::NUM_P ){
    FFGRADODE ResDer;
    // No DAG copy of ODE - reuse FFODE DAG copy
    // Caveat is that passing a pointer to the orginal ODESLV object will create a separate object
    // FFVar const*const* vResDer = ResDer( nVar, vVarVal.data(), _pODESLV, false );
    // DAG copy of ODE - no resuse of FFODE DAG copy
    // Caveat is external data pointer may change
    FFVar const*const* vResDer = ResDer._set( _nPar, vVarVal.data(), _nCst, vVarVal.data()+_nPar, static_cast<ODESLVS_CVODES*>(data), _ownODESLV );
    for( unsigned k=0; k<nRes; ++k ){
      vRes[k] = *vResVal[k];
      for( unsigned i=0; i<_nPar; ++i )
        vRes[k].setDepend( vVar[i] );
      for( unsigned j=0; j<vRes[k].size(); ++j ){
        vRes[k][j] = 0.;
        for( unsigned i=0; i<_nPar; ++i ){
          if( vVar[i][j].cst() && vVar[i][j].num().val() == 0. ) continue;
          vRes[k][j] += *vResDer[k+nRes*i] * vVar[i][j];
        }
      }
    }
  }

  else if( options.DIFF == Options::SYM_P ){
    FFODE ResDer;
    auto pODESLVSEN = _pODESLV->fdiff( _nPar, _pODESLV->var_parameter().data() );
    pODESLVSEN->setup();
    FFVar const*const* vResDer = ResDer._set( _nPar, vVarVal.data(), _nCst, vVarVal.data()+_nPar, pODESLVSEN, -1 ); // ask to transfer ownership of ODE data
    for( unsigned k=0; k<nRes; ++k ){
      vRes[k] = *vResVal[k];
      for( unsigned i=0; i<_nPar; ++i )
        vRes[k].setDepend( vVar[i] );
      for( unsigned j=0; j<vRes[k].size(); ++j ){
        vRes[k][j] = 0.;
        for( unsigned i=0; i<_nPar; ++i ){
          if( vVar[i][j].cst() && vVar[i][j].num().val() == 0. ) continue;
          vRes[k][j] += *vResDer[k+nRes*i] * vVar[i][j];
        }
      }
    }
  }

  else if( options.DIFF == Options::SYM_C ){
    FFODE ResDer;
    auto pODESLVSEN = _pODESLV->fdiff( _nCst, _pODESLV->var_constant().data() );
    pODESLVSEN->setup();
    FFVar const*const* vResDer = ResDer._set( _nPar, vVarVal.data(), _nCst, vVarVal.data()+_nPar, pODESLVSEN, -1 ); // ask to transfer ownership of ODE data
    for( unsigned k=0; k<nRes; ++k ){
      vRes[k] = *vResVal[k];
      for( unsigned i=0; i<_nCst; ++i )
        vRes[k].setDepend( vVar[_nPar+i] );
      for( unsigned j=0; j<vRes[k].size(); ++j ){
        vRes[k][j] = 0.;
        for( unsigned i=0; i<_nCst; ++i ){
          if( vVar[_nPar+i][j].cst() && vVar[_nPar+i][j].num().val() == 0. ) continue;
          vRes[k][j] += *vResDer[k+nRes*i] * vVar[_nPar+i][j];
        }
      }
    }
  }

  else if( options.DIFF == Options::SYM_PC ){
    FFODE ResDer;
    std::vector<FFVar> var_cst = _pODESLV->var_parameter();
    var_cst.insert( var_cst.end(), _pODESLV->var_constant().begin(), _pODESLV->var_constant().end() );
    auto pODESLVSEN = _pODESLV->fdiff( nVar, var_cst.data() );
    pODESLVSEN->setup();
    FFVar const*const* vResDer = ResDer._set( _nPar, vVarVal.data(), _nCst, vVarVal.data()+_nPar, pODESLVSEN, -1 ); // ask to transfer ownership of ODE data
    for( unsigned k=0; k<nRes; ++k ){
      vRes[k] = *vResVal[k];
      for( unsigned i=0; i<nVar; ++i )
        vRes[k].setDepend( vVar[i] );
      for( unsigned j=0; j<vRes[k].size(); ++j ){
        vRes[k][j] = 0.;
        for( unsigned i=0; i<nVar; ++i ){
          if( vVar[i][j].cst() && vVar[i][j].num().val() == 0. ) continue;
          vRes[k][j] += *vResDer[k+nRes*i] * vVar[i][j];
        }
      }
    }
  }
}

inline void
FFODE::deriv
( unsigned const nRes, FFVar const* vRes, unsigned const nVar, FFVar const* vVar, FFVar** vDer )
const
{
#ifdef CRONOS__FFODE_TRACE
  std::cout << "FFODE::deriv\n";
#endif
#ifdef CRONOS__FFODE_CHECK
  assert( _pODESLV && nRes == _pODESLV->nf() && nVar == _nPar+_nCst );
#endif

  if( options.DIFF == Options::NUM_P ){
    FFGRADODE ResDer;
    // No DAG copy of ODE - reuse FFODE DAG copy
    // Caveat is that passing a pointer to the orginal ODESLV object will create a separate object
    // FFVar const*const* vResDer = ResDer( nVar, vVar, _pODESLV, false ); // no copy of ODE data
    // DAG copy of ODE - no resuse of FFODE DAG copy
    // Caveat is external data pointer may change
    FFVar const*const* vResDer = ResDer._set( _nPar, vVar, _nCst, vVar+_nPar, static_cast<ODESLVS_CVODES*>(data), _ownODESLV );
    for( unsigned k=0; k<nRes; ++k )
      for( unsigned i=0; i<nVar; ++i )
        vDer[k][i] = i<_nPar? *vResDer[k+nRes*i]: 0;
  }

  else if( options.DIFF == Options::SYM_P ){
    FFODE ResDer;
    auto pODESLVSEN = _pODESLV->fdiff( _nPar, _pODESLV->var_parameter().data() );
    pODESLVSEN->setup();
    FFVar const*const* vResDer = ResDer._set( _nPar, vVar, _nCst, vVar+_nPar, pODESLVSEN, -1 ); // ask to transfer ownership of ODE data
    for( unsigned k=0; k<nRes; ++k )
      for( unsigned i=0; i<nVar; ++i )
        vDer[k][i] = i<_nPar? *vResDer[k+nRes*i]: 0;
  }

  else if( options.DIFF == Options::SYM_C ){
    FFODE ResDer;
    auto pODESLVSEN = _pODESLV->fdiff( _nCst, _pODESLV->var_constant().data() );
    pODESLVSEN->setup();
    FFVar const*const* vResDer = ResDer._set( _nPar, vVar, _nCst, vVar+_nPar, pODESLVSEN, -1 ); // ask to transfer ownership of ODE data
    for( unsigned k=0; k<nRes; ++k )
      for( unsigned i=0; i<nVar; ++i )
        vDer[k][i] = i>=_nPar? *vResDer[k+nRes*(i-_nPar)]: 0;
  }

  else if( options.DIFF == Options::SYM_PC ){
    FFODE ResDer;
    std::vector<FFVar> var_cst = _pODESLV->var_parameter();
    var_cst.insert( var_cst.end(), _pODESLV->var_constant().begin(), _pODESLV->var_constant().end() );
    auto pODESLVSEN = _pODESLV->fdiff( nVar, var_cst.data() );
    pODESLVSEN->setup();
    FFVar const*const* vResDer = ResDer._set( _nPar, vVar, _nCst, vVar+_nPar, pODESLVSEN, -1 ); // ask to transfer ownership of ODE data
    for( unsigned k=0; k<nRes; ++k )
      for( unsigned i=0; i<nVar; ++i )
        vDer[k][i] = *vResDer[k+nRes*i];
  }
}

inline void
FFODE::eval
( unsigned const nRes, SLiftVar* vRes, unsigned const nVar, SLiftVar const* vVar,
  unsigned const* mVar )
const
{
#ifdef CRONOS__FFODE_TRACE
  std::cout << "FFODE::eval: SLiftVar\n";
#endif
#ifdef CRONOS__FFODE_CHECK
  assert( _pODESLV && nRes == _pODESLV->nf() && nVar == _pODESLV->np() );
#endif

  vVar->env()->lift( nRes, vRes, nVar, vVar );
}

inline void
FFGRADODE::eval
( unsigned const nRes, FFVar* vRes, unsigned const nVar, FFVar const* vVar,
  unsigned const* mVar )
const
{
#ifdef CRONOS__FFODE_TRACE
  std::cout << "FFGRADODE::eval: FFVar\n";
#endif
#ifdef CRONOS__FFODE_CHECK
  assert( _pODESLV && nVar == _nPar+_nCst && nRes == _pODESLV->nf()*_nPar );
#endif

  FFVar** ppRes = ( _nCst? insert_external_operation( *this, nRes, _nPar, vVar, _nCst, vVar+_nPar ):
                           insert_external_operation( *this, nRes, _nPar, vVar ) );
  //FFVar** ppRes = _set( _nPar, vVar, _nCst, vVar+_nPar, static_cast<ODESLVS_CVODES*>(data), _ownODESLV );
  for( unsigned j=0; j<nRes; ++j ) vRes[j] = *(ppRes[j]);
}

inline void
FFGRADODE::eval
( unsigned const nRes, FFDep* vRes, unsigned const nVar, FFDep const* vVar,
  unsigned const* mVar )
const
{
#ifdef CRONOS__FFODE_TRACE
  std::cout << "FFGRADODE::eval: FFDep\n";
#endif
#ifdef CRONOS__FFODE_CHECK
  assert( _pODESLV && nVar == _nPar+_nCst && nRes == _pODESLV->nf()*_nPar );
#endif

  vRes[0] = 0;
  for( unsigned i=0; i<nVar; ++i ) vRes[0] += vVar[i];
  vRes[0].update( FFDep::TYPE::N );
  for( unsigned j=1; j<nRes; ++j ) vRes[j] = vRes[0];
}

inline void
FFGRADODE::eval
( unsigned const nRes, double* vRes, unsigned const nVar, double const* vVar,
  unsigned const* mVar )
const
{
#ifdef CRONOS__FFODE_TRACE
  std::cout << "FFGRADODE::eval: double\n";
#endif
#ifdef CRONOS__FFODE_CHECK
  assert( _pODESLV && nVar == _nPar+_nCst && nRes == _pODESLV->nf()*_nPar );
#endif

  if( _nPar <= options.NP2NF * _pODESLV->nf() ){
    if( _pODESLV->solve_sensitivity( vVar, _nCst? vVar+_nPar: nullptr ) != ODESLVS_CVODES::NORMAL )
      throw std::runtime_error( "FFGRADODE::eval double ** Forward sensitivity integration failure\n" );
  }
  else{
    if( _pODESLV->solve_adjoint( vVar, _nCst? vVar+_nPar: nullptr ) != ODESLVS_CVODES::NORMAL )
      throw std::runtime_error( "FFGRADODE::eval double ** Adjoint sensitivity integration failure\n" );
  }

  for( unsigned i=0, k=0; i<_nPar; ++i )
    for( unsigned j=0; j<_pODESLV->nf(); ++j, ++k )
      vRes[k] = _pODESLV->val_function_gradient()[i][j];
}

inline void
FFGRADODE::eval
( unsigned const nRes, SLiftVar* vRes, unsigned const nVar, SLiftVar const* vVar,
  unsigned const* mVar )
const
{
#ifdef CRONOS__FFODE_TRACE
  std::cout << "FFGRADODE::eval: SLiftVar\n";
#endif
#ifdef CRONOS__FFODE_CHECK
  assert( _pODESLV && nVar == _nPar+_nCst && nRes == _pODESLV->nf()*_nPar );
#endif

  vVar->env()->lift( nRes, vRes, nVar, vVar );
}
/*
inline bool
FFODE::lt
( FFOp const* op )
const
{
#ifdef MC__FFODE_TRACE
  std::cout << "FFODE::lt\n";
#endif

  // Compare data fields first
  if( data < op->data ) return true;
  if( data > op->data ) return false;

  // Compare constant vector size and values next
  auto const* extop = dynamic_cast<FFODE const*>(op);
  if( _constants.size() < extop->_constants.size() ) return true;
  if( _constants.size() > extop->_constants.size() ) return false;
  for( auto it1=_constants.cbegin(), it2=extop->_constants.cbegin();
       it1!=_constants.cend();
       ++it1, ++it2 ){
    if( *it1 < *it2 ) return true;
    if( *it1 > *it2 ) return false;
  }
  return false;    
}

inline bool
FFGRADODE::lt
( FFOp const* op )
const
{
#ifdef MC__FFODE_TRACE
  std::cout << "FFODE::lt\n";
#endif

  // Compare data fields first
  if( data < op->data ) return true;
  if( data > op->data ) return false;

  // Compare constant vector size and values next
  auto const* extop = dynamic_cast<FFGRADODE const*>(op);
  if( _constants.size() < extop->_constants.size() ) return true;
  if( _constants.size() > extop->_constants.size() ) return false;
  for( auto it1=_constants.cbegin(), it2=extop->_constants.cbegin();
       it1!=_constants.cend();
       ++it1, ++it2 ){
    if( *it1 < *it2 ) return true;
    if( *it1 > *it2 ) return false;
  }
  return false;    
}
*/
} // end namescape mc

#endif
