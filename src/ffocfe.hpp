// Copyright (C) Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

#ifndef CRONOS__FFOCFE_HPP
#define CRONOS__FFOCFE_HPP

#include "ffexpr.hpp"
#include "ffdep.hpp"
#include "slift.hpp"
#include "ocfeslv.hpp"

namespace mc
{

//! @brief C++ base class defining IPDAE collocation as external DAG operations in MC++.
////////////////////////////////////////////////////////////////////////
//! mc::FFBaseOCFE is a C++ base class for defining IPDAE
//! collocation as external DAG operation in MC++.
////////////////////////////////////////////////////////////////////////
class FFBaseOCFE
: public FFOp
{

protected:
  //! @brief Input layout of the reduced operations (two-map contract).  Op inputs are
  //! [ map-1 DOFs in control order | map-2 input DOFs in listing order | map-2 constants in listing order ].
  struct Layout
  {
  size_t                               nCtrl = 0;  //!< map-1 DOFs (= n_control_dof())
  std::vector<std::pair<FFVar,size_t>> rest;       //!< map-2 inputs: declared FFVar and DOF count, listing order
  std::vector<size_t>                  cstNdx;     //!< map-2 constants: index into var_constant(), listing order
  std::vector<FFVar>                   ctrlRec;    //!< the control registry when the op was built (SHALLOW guard)
  size_t nRest() const { size_t n = 0; for( auto const& r : rest ) n += r.second; return n; }
  size_t nIn  () const { return nCtrl + nRest() + cstNdx.size(); }
  };

  //! @brief Refuse if the solver's control registry is no longer the one the op was built with.
  static void _check_controls( OCFESLV const* oc, Layout const& L );
  //! @brief From the op inputs (layout L): state starting values from the references (init), the input vector from
  //! map 1 (decode_controls) and map 2 (set_input_values), and the constants.
  static void _assemble( OCFESLV* oc, Layout const& L, double const* v, std::vector<double>& xv, std::vector<double>& inp, std::vector<double>& cst );
  //! @brief Resolve the two maps, RESET the registry to map 1, and return the op inputs in layout order.
  static std::vector<FFVar> _bind( std::vector<FFModel::InputArg> const& vDiff, std::vector<FFModel::InputArg> const& vRest, OCFESLV* oc, std::vector<FFVar>& saved, Layout& L );
  //! @brief Reduced Jacobian of the collocated outputs w.r.t. the op inputs (see its definition).
  static bool _reduced_jacobian( OCFESLV* pOCFESLV, Layout const& L, double const* v, size_t const nFct, int mode, double np2nf, std::vector<double>& Jflat /* size nFct*L.nCtrl, row-major (irow*nCtrl+jcol) */ );


  //! @brief Pointer to collocation environment
  OCFESLV*              _pOCFESLV;
  //! @brief Whether this class owns _pOCFESLV
  bool                _ownOCFESLV;
  //! @brief Number of constants
  size_t              _nCst;
  //! @brief Number of collocated states
  size_t              _nStaColl;
  //! @brief Number of collocated inputs
  size_t              _nInpColl;
  //! @brief Number of collocated states
  size_t              _nEqnColl;
  //! @brief Number of collocated inputs
  size_t              _nFctColl;
  //! @brief Number of nonzero entries in collocated states
  size_t              _nnzEqnColl;
  //! @brief Number of nonzero entries in collocated inputs
  size_t              _nnzFctColl;
  //! @brief Name of collocated IPDAE model
  std::string         _name;

  //! @brief Is an FFVar evaluation re-inserting this operation into ANOTHER DAG (FFGraph::insert -- e.g. the
  //! per-thread copies of FFGraph::veval)?  Then the copy must own a deep copy of the solver: a SHALLOW operation
  //! would otherwise share one OCFESLV between threads, which solve concurrently (2026-10-03: wrong values, failed
  //! scenarios and crashes in the OCFESLV tutorial's veval with MAXTHREAD = 0).
  bool _into_other_dag
    ( unsigned const nVar, FFVar const* vVar )
    const
    {
      if( _ownOCFESLV || !nVar || varin.empty() ) return false;
      auto const* tgt = vVar[0].dag();
      for( unsigned i = 1; !tgt && i < nVar; ++i ) tgt = vVar[i].dag();
      decltype( tgt ) src = nullptr;
      for( auto const* v : varin ) if( v && v->dag() ){ src = v->dag(); break; }
      return tgt && src && tgt != src;
    }

  //! @brief The layout part of the operation's identity (empty: none -- FFOCFERES).  FFOCFESLV and FFGradOCFESLV
  //! return their input layout, so that two embeddings of the same solver over the same variables but with a
  //! different map split are different operations.
  virtual std::vector<size_t> _layout_key
    ()
    const
    { return std::vector<size_t>(); }

  //! @brief An input layout as a comparable key: map-1 DOF count and control registry, map-2 inputs with their DOF
  //! counts, map-2 constants.
  static std::vector<size_t> _layout_key_of
    ( Layout const& L )
    {
      std::vector<size_t> k{ L.nCtrl, L.ctrlRec.size() };
      for( auto const& u : L.ctrlRec ) k.push_back( static_cast<size_t>( u.id().second ) );
      k.push_back( L.rest.size() );
      for( auto const& [w, n] : L.rest ){ k.push_back( static_cast<size_t>( w.id().second ) ); k.push_back( n ); }
      k.push_back( L.cstNdx.size() );
      k.insert( k.end(), L.cstNdx.begin(), L.cstNdx.end() );
      return k;
    }

public:
  //! @brief Ordering of the operation in the DAG (2026-10-01).  lt_FFOp has already compared the type, info and
  //! operands; then: (1) the solver (data), (2) whether the op owns a copy of it (COPY) or refers to it (SHALLOW),
  //! (3) the input layout (_layout_key), so that a one-map and a two-map embedding of the same solver over the same
  //! variables are two operations (their derivatives differ) -- as FFBaseODESLV::lt.
  bool lt
    ( FFOp const* op )
    const
    {
      if( data != op->data ) return data < op->data;
      auto const* o = dynamic_cast<FFBaseOCFE const*>( op );
      if( !o ) return false;
      if( _ownOCFESLV != o->_ownOCFESLV ) return _ownOCFESLV < o->_ownOCFESLV;
      return _layout_key() < o->_layout_key();
    }


  //! @brief Default constructor
  FFBaseOCFE
    ()
    : FFOp        ( EXTERN ),
      _pOCFESLV     ( nullptr ),
      _ownOCFESLV   ( false ),
      _nCst       ( 0 ),
      _nStaColl   ( 0 ),
      _nInpColl   ( 0 ),
      _nEqnColl   ( 0 ),
      _nFctColl   ( 0 ),
      _nnzEqnColl ( 0 ),
      _nnzFctColl ( 0 ),
      _name       ( "" )
    {}

  //! @brief Destructor
  virtual ~FFBaseOCFE
    ()
    {
#ifdef CRONOS__FFOCFE_TRACE
      std::cout << "FFBaseOCFE::destructor\n";
#endif
      if( _ownOCFESLV && _pOCFESLV )
        delete _pOCFESLV;
    }

  //! @brief Copy constructor
  FFBaseOCFE
    ( FFBaseOCFE const& Op )
    : FFOp        ( Op ),
      _nCst       ( Op._nCst ),
      _nStaColl   ( Op._nStaColl ),
      _nInpColl   ( Op._nInpColl ),
      _nEqnColl   ( Op._nEqnColl ),
      _nFctColl   ( Op._nFctColl ),
      _nnzEqnColl ( Op._nnzEqnColl ),
      _nnzFctColl ( Op._nnzFctColl ),
      _name      ( Op._name )
    {
#ifdef CRONOS__FFOCFE_TRACE
      std::cout << "FFBaseOCFE::copy constructor\n";
#endif
      if( !Op._pOCFESLV )
        throw std::runtime_error( "FFBaseOCFE::copy constructor ** Null pointer to collocation environment\n" );

      _ownOCFESLV = Op._ownOCFESLV;      
      if( _ownOCFESLV ){
        _pOCFESLV = new OCFESLV;
        _pOCFESLV->deep_copy_from( *Op._pOCFESLV );
#ifdef CRONOS__FFOCFE_TRACE
        std::cerr << "OCFESLV address copied: " << _pOCFESLV << std::endl;
#endif
      }
      else
        _pOCFESLV = Op._pOCFESLV;
    }

  //! @brief Reduced-space gradient mode for FFOCFESLV / FFGradOCFESLV: FORWARD uses one
  //! forward-sensitivity march per control DOF (n_control_dof applies); ADJOINT uses one
  //! reverse march per output function (n_colloc_fct sweeps); AUTO picks adjoint when
  //! n_control_dof > NP2NF * n_colloc_fct (mirrors FFGradODE's NP2NF heuristic).
  enum GRADIENT_TYPE{ FORWARD=0, ADJOINT=1, AUTO=2 };

  //! @brief options - static so it remains accessible.modifiable after setup
  static struct Options
  {
    //! @brief Constructor
    Options():
      SYMDIFF (),
      NP2NF   (3.),
      GRADIENT( AUTO )
      {}
    //! @brief Assignment operator
    Options& operator= ( Options const& options ){
        SYMDIFF   = options.SYMDIFF;
        NP2NF     = options.NP2NF;
        GRADIENT  = options.GRADIENT;
        return *this;
      }
    //! @brief Variable selection for symbolic differentiation - needs to be participating inputs or constants. Applies numerical differentiation w.r.t. states/inputs if empty
    std::vector<FFVar>        SYMDIFF;
    //! @brief control-to-function ratio above which the reduced-space AUTO gradient uses adjoint instead of forward sensitivity
    double                    NP2NF;
    //! @brief reduced-space gradient mode (FORWARD / ADJOINT / AUTO) for FFOCFESLV / FFGradOCFESLV
    int                       GRADIENT;
  } options;
};

inline FFBaseOCFE::Options FFBaseOCFE::options;

//! @brief C++ class defining IPDAE collocation evaluation as external DAG operations in MC++.
////////////////////////////////////////////////////////////////////////
//! mc::FFOCFERES is a C++ class for defining IPDAE collocation
//! evaluation as external DAG operation in MC++.
////////////////////////////////////////////////////////////////////////
class FFOCFERES
: public FFBaseOCFE
{

protected:

  FFVar** _set
    ( size_t const nStaColl, FFVar const* pStaColl,
      size_t const nInpColl, FFVar const* pInpColl, 
      size_t const nCst,     FFVar const* pCst,
      OCFESLV* pOCFESLV, int policy, std::string const& name )
    {
#ifdef CRONOS__FFOCFE_CHECK
      assert( pOCFESLV
           && nStaColl == pOCFESLV->n_colloc_sta() 
           && nInpColl == pOCFESLV->n_colloc_inp() 
           && ( !pCst || nCst == pOCFESLV->var_constant().size() ) );
#endif
      if( _ownOCFESLV && _pOCFESLV )
        delete _pOCFESLV;
      _ownOCFESLV   = ( policy>0? true: false ); //copy;
      _pOCFESLV     = pOCFESLV;
      _nCst       = nCst;
      _nStaColl   = nStaColl;
      _nInpColl   = nInpColl;
      _nEqnColl   = pOCFESLV->n_colloc_eqn();
      _nFctColl   = pOCFESLV->n_colloc_fct();
      _nnzEqnColl = pOCFESLV->n_colloc_eqn_nnz();
      _nnzFctColl = pOCFESLV->n_colloc_fct_nnz();
      _name     = name;

      data    = pOCFESLV;
      owndata = false;
      size_t const nDep = _nEqnColl + _nFctColl;
      FFVar** ppRes = insert_external_operation( *this, nDep,
                                                 nStaColl, pStaColl,
                                                 nInpColl, pInpColl,
                                                 nCst,     pCst );

      _ownOCFESLV = false;
      FFOp* pOp = (*ppRes)->opdef().first;
      if( policy > 0 )
        _pOCFESLV = dynamic_cast<FFOCFERES*>(pOp)->_pOCFESLV; // set pointer to DAG copy
      else if( policy < 0 )
        dynamic_cast<FFOCFERES*>(pOp)->_ownOCFESLV = true;   // transfer ownership
#ifdef CRONOS__FFOCFE_TRACE
      std::cerr << "FFOCFERES operation address: " << this << std::endl;
      std::cerr << "OCFESLV address in DAG: " << _pOCFESLV << std::endl;
#endif
      return ppRes;
    }

public:

  //! @brief Enumeration type for OCFESLV copy policy
  enum POLICY_TYPE{
    SHALLOW=0,  //!< Shallow copy of OCFESLV in FFGraph (without ownership)
    COPY=1,     //!< Deep copy of OCFESLV in FFGraph (with ownership)
    TRANSFER=-1 //!< Shallow copy of OCFESLV in FFGraph (with ownership transfer)
  };

  // Default constructor
  FFOCFERES
    ()// bool const sparse=true )
    : FFBaseOCFE()
    {
      this->sparse = true;//sparse;
    }

  // Copy constructor
  FFOCFERES
    (   FFOCFERES const& Other )
    :   FFBaseOCFE( Other )
    {}

  // Destructor
  virtual ~FFOCFERES
    ()
    {}

  // Define operation
  std::vector<FFVar> operator()
    ( std::vector<FFVar> const& vStaColl, std::vector<FFVar> const& vInpColl,
      std::vector<FFVar> const& vCst, OCFESLV* pOCFESLV, int policy=COPY, std::string const& name="" )
    {
      FFVar** ppRes = _set( vStaColl.size(), vStaColl.data(), vInpColl.size(), vInpColl.data(),
                            vCst.size(), vCst.data(), pOCFESLV, policy, name );
      std::vector<FFVar> vDep( _nEqnColl + _nFctColl );
      for( size_t i=0; i<vDep.size(); ++i ) vDep[i] = *ppRes[i];
      return vDep;
    }

  FFVar& operator()
    ( size_t const iDep, std::vector<FFVar> const& vStaColl, std::vector<FFVar> const& vInpColl,
      std::vector<FFVar> const& vCst, OCFESLV* pOCFESLV, int policy=COPY, std::string const& name="" )
    {
#ifdef CRONOS__FFOCFE_CHECK
      assert( iDep < _nEqnColl + _nFctColl );
#endif
      return *(_set( vStaColl.size(), vStaColl.data(), vInpColl.size(), vInpColl.data(),
                     vCst.size(), vCst.data(), pOCFESLV, policy, name )[iDep]);
    }

  FFVar** operator()
    ( size_t const nStaColl, FFVar const* pStaColl, size_t const nInpColl, FFVar const* pInpColl,
      size_t const nCst, FFVar const* pCst, OCFESLV* pOCFESLV, int policy=COPY, std::string const& name="" )
    {
      return _set( nStaColl, pStaColl, nInpColl, pInpColl, nCst, pCst, pOCFESLV, policy, name );
    }

  FFVar& operator()
    ( size_t const iDep, size_t const nStaColl, FFVar const* pStaColl, size_t const nInpColl, FFVar const* pInpColl,
      size_t const nCst, FFVar const* pCst, OCFESLV* pOCFESLV, int policy=COPY, std::string const& name="" )
    {
#ifdef CRONOS__FFOCFE_CHECK
      assert( iDep < _nEqnColl + _nFctColl );
#endif
      return *(_set( nStaColl, pStaColl, nInpColl, pInpColl, nCst, pCst, pOCFESLV, policy, name )[iDep]);
    }

  OCFESLV* pOCFESLV
    ()
    const
    { //std::cerr << "OCFESLV address retreived: " << _pOCFESLV << std::endl;
      return _pOCFESLV; }

  // Evaluation overloads
  virtual void feval
    ( std::type_info const& idU, unsigned const nRes, void* vRes, unsigned const nVar,
      void const* vVar, unsigned const* mVar )
    const
    {
      if( mc::same_type( idU, typeid( FFVar ) ) )
        return eval( nRes, static_cast<FFVar*>(vRes), nVar, static_cast<FFVar const*>(vVar), mVar );
      else if( mc::same_type( idU, typeid( FADType<FFVar> ) ) )
        return eval( nRes, static_cast<FADType<FFVar>*>(vRes), nVar, static_cast<FADType<FFVar> const*>(vVar), mVar );
      else if( mc::same_type( idU, typeid( FFDep ) ) )
        return eval( nRes, static_cast<FFDep*>(vRes), nVar, static_cast<FFDep const*>(vVar), mVar );
      else if( mc::same_type( idU, typeid( double ) ) )
        return eval( nRes, static_cast<double*>(vRes), nVar, static_cast<double const*>(vVar), mVar );
      else if( mc::same_type( idU, typeid( FADType<double> ) ) )
        return eval( nRes, static_cast<FADType<double>*>(vRes), nVar, static_cast<FADType<double> const*>(vVar), mVar );
      else if( mc::same_type( idU, typeid( BADType<double> ) ) )
        return eval( nRes, static_cast<BADType<double>*>(vRes), nVar, static_cast<BADType<double> const*>(vVar), mVar );
      else if( mc::same_type( idU, typeid( SLiftVar ) ) )
        return eval( nRes, static_cast<SLiftVar*>(vRes), nVar, static_cast<SLiftVar const*>(vVar), mVar );
      else if( mc::same_type( idU, typeid( FFExpr ) ) )
        return eval( nRes, static_cast<FFExpr*>(vRes), nVar, static_cast<FFExpr const*>(vVar), mVar );

      throw std::runtime_error( "FFOCFERES::feval ** No evaluation method for type"+std::string(idU.name())+"\n" );
    }

  void eval
    ( unsigned const nRes, double* vRes, unsigned const nVar, double const* vVar, unsigned const* mVar )
    const;

  void eval
    ( unsigned const nRes, FADType<double>* vRes, unsigned const nVar, FADType<double> const* vVar,
      unsigned const* mVar )
    const;

  void eval
    ( unsigned const nRes, BADType<double>* vRes, unsigned const nVar, BADType<double> const* vVar,
      unsigned const* mVar )
    const;

  void eval
    ( unsigned const nRes, FFVar* vRes, unsigned const nVar, FFVar const* vVar, unsigned const* mVar )
    const;

  void eval
    ( unsigned const nRes, FFDep* vRes, unsigned const nVar, FFDep const* vVar, unsigned const* mVar )
    const;

  void eval
    ( unsigned const nRes, FADType<FFVar>* vRes, unsigned const nVar, FADType<FFVar> const* vVar,
      unsigned const* mVar )
    const;

  void eval
    ( unsigned const nRes, SLiftVar* vRes, unsigned const nVar, SLiftVar const* vVar, unsigned const* mVar )
    const;

  void eval
    ( unsigned const nRes, FFExpr* vRes, unsigned const nVar, FFExpr const* vVar, unsigned const* mVar )
    const;

  // Derivatives
  void deriv
    ( unsigned const nRes, FFVar const* vRes, unsigned const nVar, FFVar const* vVar, FFVar** vDer )
    const;

  void deriv
    ( unsigned const nRes, FFVar const* vRes, unsigned const nVar, FFVar const* vVar, FFVar** vDer,
      size_t* nnz, size_t** colnz )
    const;

//  // Ordering
  // Ordering in the DAG: FFBaseOCFE::lt (solver, policy, layout)
//  bool lt
//    ( FFOp const* op )
//    const;

  // Properties
  std::string name
    ()
    const
    { std::ostringstream oss;
      if( !_name.empty() )  oss << _name;
      else                  oss << _pOCFESLV;
      return "OCFE_RES[" + oss.str() + "]"; }

  //! @brief Return whether or not operation is commutative
  bool commutative
    ()
    const
    { return false; }
};

//! @brief C++ class defining IPDAE collocation gradient evaluation as external DAG operations in MC++.
////////////////////////////////////////////////////////////////////////
//! mc::FFGradOCFERES is a C++ class for defining IPDAE collocation
//! gradient evaluation as external DAG operation in MC++.
////////////////////////////////////////////////////////////////////////
class FFGradOCFERES
: public FFBaseOCFE
{

protected:

  FFVar** _set
    ( size_t const nStaColl, FFVar const* pStaColl,
      size_t const nInpColl, FFVar const* pInpColl, 
      size_t const nCst,     FFVar const* pCst,
      OCFESLV* pOCFESLV, int policy, std::string const& name )
    {
#ifdef CRONOS__FFOCFE_CHECK
      assert( pOCFESLV
           && nStaColl == pOCFESLV->n_colloc_sta() 
           && nInpColl == pOCFESLV->n_colloc_inp() 
           && ( !pCst || nCst == pOCFESLV->var_constant().size() ) );
#endif
      if( _ownOCFESLV && _pOCFESLV )
        delete _pOCFESLV;
      _ownOCFESLV = ( policy>0? true: false ); //copy;
      _pOCFESLV   = pOCFESLV;
      _nCst       = nCst;
      _nStaColl   = nStaColl;
      _nInpColl   = nInpColl;
      _nEqnColl   = pOCFESLV->n_colloc_eqn();
      _nFctColl   = pOCFESLV->n_colloc_fct();
      _nnzEqnColl = pOCFESLV->n_colloc_eqn_nnz();
      _nnzFctColl = pOCFESLV->n_colloc_fct_nnz();
      _name     = name;

      data    = pOCFESLV;
      owndata = false;
      size_t const nDep = _nnzEqnColl + _nnzFctColl;
      FFVar** ppRes = insert_external_operation( *this, nDep,
                                                 nStaColl, pStaColl,
                                                 nInpColl, pInpColl,
                                                 nCst,     pCst );

      _ownOCFESLV = false;
      FFOp* pOp = (*ppRes)->opdef().first;
      if( policy > 0 )
        _pOCFESLV = dynamic_cast<FFGradOCFERES*>(pOp)->_pOCFESLV; // set pointer to DAG copy
      //else if( policy < 0 )
      //  dynamic_cast<FFGradOCFERES*>(pOp)->_ownOCFESLV = true;   // transfer ownership
#ifdef CRONOS__FFOCFE_TRACE
      std::cerr << "FFGradOCFERES operation address: " << this << std::endl;
      std::cerr << "OCFESLV address in DAG: " << _pOCFESLV << std::endl;
#endif
      return ppRes;
    }

public:

  //! @brief Enumeration type for OCFESLV copy policy
  enum POLICY_TYPE{
    SHALLOW=0,  //!< Shallow copy of OCFESLV in FFGraph (without ownership)
    COPY=1     //!< Deep copy of OCFESLV in FFGraph (with ownership)
    //TRANSFER=-1 //!< Shallow copy of OCFESLV in FFGraph (with ownership transfer)
  };

  // Default constructor
  FFGradOCFERES
    ()// bool const sparse=true )
    : FFBaseOCFE()
    {
      this->sparse = true;//sparse;
    }

  // Copy constructor
  FFGradOCFERES
    (   FFGradOCFERES const& Other )
    :   FFBaseOCFE( Other )
    {}

  // Destructor
  virtual ~FFGradOCFERES
    ()
    {}

  // Define operation
  std::vector<FFVar> operator()
    ( std::vector<FFVar> const& vStaColl, std::vector<FFVar> const& vInpColl,
      std::vector<FFVar> const& vCst, OCFESLV* pOCFESLV, int policy=COPY, std::string const& name="" )
    {
      FFVar** ppRes = _set( vStaColl.size(), vStaColl.data(), vInpColl.size(), vInpColl.data(),
                            vCst.size(), vCst.data(), pOCFESLV, policy, name );
      std::vector<FFVar> vDep( _nnzEqnColl + _nnzFctColl );
      for( size_t i=0; i<vDep.size(); ++i ) vDep[i] = *ppRes[i];
      return vDep;
    }

  FFVar& operator()
    ( size_t const iDep, std::vector<FFVar> const& vStaColl, std::vector<FFVar> const& vInpColl,
      std::vector<FFVar> const& vCst, OCFESLV* pOCFESLV, int policy=COPY, std::string const& name="" )
    {
#ifdef CRONOS__FFOCFE_CHECK
      assert( iDep < pOCFESLV->n_colloc_rows() );
#endif
      return *(_set( vStaColl.size(), vStaColl.data(), vInpColl.size(), vInpColl.data(),
                     vCst.size(), vCst.data(), pOCFESLV, policy, name )[iDep]);
    }

  FFVar** operator()
    ( size_t const nStaColl, FFVar const* pStaColl, size_t const nInpColl, FFVar const* pInpColl,
      size_t const nCst, FFVar const* pCst, OCFESLV* pOCFESLV, int policy=COPY, std::string const& name="" )
    {
      return _set( nStaColl, pStaColl, nInpColl, pInpColl, nCst, pCst, pOCFESLV, policy, name );
    }

  FFVar& operator()
    ( size_t const iDep, size_t const nStaColl, FFVar const* pStaColl, size_t const nInpColl, FFVar const* pInpColl,
      size_t const nCst, FFVar const* pCst, OCFESLV* pOCFESLV, int policy=COPY, std::string const& name="" )
    {
#ifdef CRONOS__FFOCFE_CHECK
      assert( iDep < pOCFESLV->n_colloc_rows() );
#endif
      return *(_set( nStaColl, pStaColl, nInpColl, pInpColl, nCst, pCst, pOCFESLV, policy, name )[iDep]);
    }

  OCFESLV* pOCFESLV
    ()
    const
    { //std::cerr << "OCFESLV address retreived: " << _pOCFESLV << std::endl;
      return _pOCFESLV; }

  // Evaluation overloads
  virtual void feval
    ( std::type_info const& idU, unsigned const nRes, void* vRes, unsigned const nVar,
      void const* vVar, unsigned const* mVar )
    const
    {
      if( mc::same_type( idU, typeid( FFVar ) ) )
        return eval( nRes, static_cast<FFVar*>(vRes), nVar, static_cast<FFVar const*>(vVar), mVar );
      else if( mc::same_type( idU, typeid( FFDep ) ) )
        return eval( nRes, static_cast<FFDep*>(vRes), nVar, static_cast<FFDep const*>(vVar), mVar );
      else if( mc::same_type( idU, typeid( double ) ) )
        return eval( nRes, static_cast<double*>(vRes), nVar, static_cast<double const*>(vVar), mVar );
      else if( mc::same_type( idU, typeid( SLiftVar ) ) )
        return eval( nRes, static_cast<SLiftVar*>(vRes), nVar, static_cast<SLiftVar const*>(vVar), mVar );
      else if( mc::same_type( idU, typeid( FFExpr ) ) )
        return eval( nRes, static_cast<FFExpr*>(vRes), nVar, static_cast<FFExpr const*>(vVar), mVar );

      throw std::runtime_error( "FFGradOCFERES::feval ** No evaluation method for type"+std::string(idU.name())+"\n" );
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

  void eval
    ( unsigned const nRes, FFExpr* vRes, unsigned const nVar, FFExpr const* vVar, unsigned const* mVar )
    const;

//  // Ordering
  // Ordering in the DAG: FFBaseOCFE::lt (solver, policy, layout)
//  bool lt
//    ( FFOp const* op )
//    const;

  // Properties
  std::string name
    ()
    const
    { std::ostringstream oss;
      if( !_name.empty() )  oss << _name;
      else                  oss << _pOCFESLV;
      return "Grad_OCFE_RES[" + oss.str() + "]"; }

  //! @brief Return whether or not operation is commutative
  bool commutative
    ()
    const
    { return false; }
};


inline void
FFOCFERES::eval
( unsigned const nRes, FFVar* vRes, unsigned const nVar, FFVar const* vVar,
  unsigned const* mVar )
const
{
#ifdef CRONOS__FFOCFE_TRACE
  std::cout << "FFOCFERES::eval: FFVar\n";
  std::cerr << "IPDAEColl operation address: " << this << std::endl;
  std::cerr << "OCFESLV address in DAG: " << _pOCFESLV << std::endl;
#endif
#ifdef CRONOS__FFOCFE_CHECK
  assert( _pOCFESLV && nRes == _nEqnColl+_nFctColl && nVar == _nStaColl+_nInpColl+_nCst );
#endif

  FFVar** ppRes = nullptr;
  if( _into_other_dag( nVar, vVar ) ){
    auto op = *this;                  // shares the solver ...
    op._ownOCFESLV = true;            // ... so that the INSERTED copy deep-copies it (FFBaseOCFE copy constructor)
    ppRes = insert_external_operation( op, nRes, _nStaColl, vVar, _nInpColl, vVar+_nStaColl, _nCst, vVar+_nStaColl+_nInpColl );
    op._ownOCFESLV = false;           // the temporary must not delete the shared solver
  }
  else
    ppRes = insert_external_operation( *this, nRes, _nStaColl, vVar, _nInpColl, vVar+_nStaColl, _nCst, vVar+_nStaColl+_nInpColl );
  for( unsigned j=0; j<nRes; ++j )
    vRes[j] = *(ppRes[j]);
}

inline void
FFOCFERES::eval
( unsigned const nRes, FFDep* vRes, unsigned const nVar, FFDep const* vVar,
  unsigned const* mVar )
const
{
#ifdef CRONOS__FFOCFERES_TRACE
  std::cout << "FFOCFERES::eval: FFDep\n";
#endif
#ifdef CRONOS__FFOCFERES_CHECK
  assert( _pOCFESLV && nRes == _nEqnColl+_nFctColl && nVar == _nStaColl+_nInpColl+_nCst );
#endif

  std::vector<double> vCst( _nCst, 0. );
  if( !_pOCFESLV->eval( vRes, vRes+_nEqnColl, vVar, vVar+_nStaColl,
                      _nCst ? vCst.data() : nullptr ) )
    throw std::runtime_error( "FFOCFERES::eval FFDep ** Failure\n" );
}

inline void
FFOCFERES::eval
( unsigned const nRes, SLiftVar* vRes, unsigned const nVar, SLiftVar const* vVar,
  unsigned const* mVar )
const
{
#ifdef CRONOS__FFOCFERES_TRACE
  std::cout << "FFOCFERES::eval: SLiftVar\n";
#endif
#ifdef CRONOS__FFOCFERES_CHECK
  assert( _pOCFESLV && nRes == _nEqnColl+_nFctColl && nVar == _nStaColl+_nInpColl+_nCst );
#endif

  vVar->env()->lift( nRes, vRes, nVar, vVar );
}

inline void
FFOCFERES::eval
( unsigned const nRes, FFExpr* vRes, unsigned const nVar, FFExpr const* vVar,
  unsigned const* mVar )
const
{
#ifdef CRONOS__FFOCFERES_TRACE
  std::cout << "FFOCFERES::eval: FFExpr\n";
#endif
#ifdef CRONOS__FFOCFERES_CHECK
  assert( _pOCFESLV && nRes == _nEqnColl+_nFctColl && nVar == _nStaColl+_nInpColl+_nCst );
#endif

  switch( FFExpr::options.LANG ){
   case FFExpr::Options::DAG:
    for( unsigned j=0; j<nRes; ++j ){
      std::ostringstream os; os << name() << "[" << j << "]";
      vRes[j] = FFExpr::compose( os.str(), nVar, vVar );
    }
    break;
   case FFExpr::Options::GAMS:
   default:
    throw typename FFExpr::Exceptions( FFExpr::Exceptions::UNDEF );
  }
}

inline void
FFOCFERES::eval
( unsigned const nRes, double* vRes, unsigned const nVar, double const* vVar,
  unsigned const* mVar )
const
{
#ifdef CRONOS__FFOCFERES_TRACE
  std::cerr << "FFOCFERES::eval: double w/ OCFESLV address " << _pOCFESLV << std::endl;
  for( unsigned i=0; i<nVar; ++i ) std::cout << "vVar[" << i << "] = " << vVar[i] << std::endl;
#endif
#ifdef CRONOS__FFOCFERES_CHECK
  assert( _pOCFESLV && nRes == _nEqnColl+_nFctColl && nVar == _nStaColl+_nInpColl+_nCst );
#endif

  if( !_pOCFESLV->eval( vRes, vRes+_nEqnColl, vVar, vVar+_nStaColl, vVar+_nStaColl+_nInpColl ) )
    throw std::runtime_error( "FFOCFERES::eval double ** Failure\n" );
}


inline void
FFOCFERES::eval
( unsigned const nRes, FADType<double>* vRes, unsigned const nVar, FADType<double> const* vVar,
  unsigned const* mVar )
const
{
#ifdef CRONOS__FFOCFERES_TRACE
  std::cout << "FFOCFERES::eval: FADType<double> via OCFESLV::eval/deriv<double>\n";
#endif
#ifdef CRONOS__FFOCFERES_CHECK
  assert( _pOCFESLV && nRes == _nEqnColl+_nFctColl && nVar == _nStaColl+_nInpColl+_nCst );
#endif

  (void)mVar;
  std::vector<double> vVarVal( nVar );
  std::vector<double> vResVal( nRes );
  for( unsigned i=0; i<nVar; ++i ) vVarVal[i] = vVar[i].val();

  eval( nRes, vResVal.data(), nVar, vVarVal.data(), nullptr );
  for( unsigned k=0; k<nRes; ++k ) vRes[k] = vResVal[k];

  bool dependent = false;
  unsigned nDir = 0;
  for( unsigned i=0; i<nVar; ++i ){
    if( !vVar[i].depend() ) continue;
    if( !dependent ) nDir = vVar[i].size();
    else if( vVar[i].size() != nDir )
      throw std::runtime_error( "FFOCFERES::eval FADType<double> ** Inconsistent derivative vector sizes\n" );
    dependent = true;
  }
  if( !dependent ) return;

  std::vector<size_t> nnz( nRes, 0 );
  if( !_pOCFESLV->deriv( nnz.data(), nullptr ) )
    throw std::runtime_error( "FFOCFERES::eval FADType<double> ** Failure retrieving sparsity counts\n" );

  std::vector< std::vector<size_t> > col_store( nRes );
  std::vector<size_t*> col_ptr( nRes, nullptr );
  for( unsigned irow=0; irow<nRes; ++irow ){
    col_store[irow].resize( nnz[irow] );
    col_ptr[irow] = col_store[irow].empty()? nullptr: col_store[irow].data();
  }
  if( !_pOCFESLV->deriv( nnz.data(), col_ptr.data() ) )
    throw std::runtime_error( "FFOCFERES::eval FADType<double> ** Failure retrieving sparsity columns\n" );

  std::vector<double> gradEqn( _nnzEqnColl, 0. );
  std::vector<double> gradFct( _nnzFctColl, 0. );
  if( !_pOCFESLV->deriv( gradEqn.empty()? nullptr: gradEqn.data(),
                       gradFct.empty()? nullptr: gradFct.data(),
                       vVarVal.data(),
                       _nInpColl? vVarVal.data()+_nStaColl: nullptr,
                       _nCst?     vVarVal.data()+_nStaColl+_nInpColl: nullptr ) )
    throw std::runtime_error( "FFOCFERES::eval FADType<double> ** Failure evaluating sparse derivative\n" );

  std::vector<double> grad( gradEqn );
  grad.insert( grad.end(), gradFct.begin(), gradFct.end() );

  size_t iel = 0;
  for( unsigned irow=0; irow<nRes; ++irow ){
    bool rowDependent = false;
    for( size_t inz=0; inz<nnz[irow]; ++inz ){
      size_t const jcol = col_store[irow][inz];
      if( jcol < nVar && vVar[jcol].depend() ){
        vRes[irow].setDepend( vVar[jcol] );
        rowDependent = true;
      }
    }
    if( rowDependent ){
      for( unsigned idir=0; idir<nDir; ++idir ){
        vRes[irow][idir] = 0.;
        for( size_t inz=0; inz<nnz[irow]; ++inz ){
          size_t const jcol = col_store[irow][inz];
          if( jcol >= nVar || !vVar[jcol].depend() ) continue;
          vRes[irow][idir] += grad[iel+inz] * vVar[jcol][idir];
        }
      }
    }
    iel += nnz[irow];
  }
  if( iel != grad.size() )
    throw std::runtime_error( "FFOCFERES::eval FADType<double> ** Sparse derivative size mismatch\n" );
}

inline void
FFOCFERES::eval
( unsigned const nRes, BADType<double>* vRes, unsigned const nVar, BADType<double> const* vVar,
  unsigned const* mVar )
const
{
#ifdef CRONOS__FFOCFERES_TRACE
  std::cout << "FFOCFERES::eval: BADType<double> via OCFESLV::eval/deriv<double>\n";
#endif
#ifdef CRONOS__FFOCFERES_CHECK
  assert( _pOCFESLV && nRes == _nEqnColl+_nFctColl && nVar == _nStaColl+_nInpColl+_nCst );
#endif

  (void)mVar;
  std::vector<double> vVarVal( nVar );
  std::vector<double> vResVal( nRes );
  for( unsigned i=0; i<nVar; ++i ) vVarVal[i] = vVar[i].val();

  eval( nRes, vResVal.data(), nVar, vVarVal.data(), nullptr );

  std::vector<size_t> nnz( nRes, 0 );
  if( !_pOCFESLV->deriv( nnz.data(), nullptr ) )
    throw std::runtime_error( "FFOCFERES::eval BADType<double> ** Failure retrieving sparsity counts\n" );

  std::vector< std::vector<size_t> > col_store( nRes );
  std::vector<size_t*> col_ptr( nRes, nullptr );
  for( unsigned irow=0; irow<nRes; ++irow ){
    col_store[irow].resize( nnz[irow] );
    col_ptr[irow] = col_store[irow].empty()? nullptr: col_store[irow].data();
  }
  if( !_pOCFESLV->deriv( nnz.data(), col_ptr.data() ) )
    throw std::runtime_error( "FFOCFERES::eval BADType<double> ** Failure retrieving sparsity columns\n" );

  std::vector<double> gradEqn( _nnzEqnColl, 0. );
  std::vector<double> gradFct( _nnzFctColl, 0. );
  if( !_pOCFESLV->deriv( gradEqn.empty()? nullptr: gradEqn.data(),
                       gradFct.empty()? nullptr: gradFct.data(),
                       vVarVal.data(),
                       _nInpColl? vVarVal.data()+_nStaColl: nullptr,
                       _nCst?     vVarVal.data()+_nStaColl+_nInpColl: nullptr ) )
    throw std::runtime_error( "FFOCFERES::eval BADType<double> ** Failure evaluating sparse derivative\n" );

  std::vector<double> grad( gradEqn );
  grad.insert( grad.end(), gradFct.begin(), gradFct.end() );

  size_t iel = 0;
  for( unsigned irow=0; irow<nRes; ++irow ){
    vRes[irow] = vResVal[irow];
    for( size_t inz=0; inz<nnz[irow]; ++inz ){
      size_t const jcol = col_store[irow][inz];
      if( jcol >= nVar || grad[iel+inz] == 0. ) continue;
      vRes[irow] += grad[iel+inz] * ( vVar[jcol] - vVarVal[jcol] );
    }
    iel += nnz[irow];
  }
  if( iel != grad.size() )
    throw std::runtime_error( "FFOCFERES::eval BADType<double> ** Sparse derivative size mismatch\n" );
}

inline void
FFOCFERES::eval
( unsigned const nRes, FADType<FFVar>* vRes, unsigned const nVar, FADType<FFVar> const* vVar,
  unsigned const* mVar )
const
{
#ifdef CRONOS__FFOCFERES_TRACE
  std::cout << "FFOCFERES::eval: FADType<FFVar> via FFGradOCFERES\n";
#endif
#ifdef CRONOS__FFOCFERES_CHECK
  assert( _pOCFESLV && nRes == _nEqnColl+_nFctColl && nVar == _nStaColl+_nInpColl+_nCst );
#endif

  (void)mVar;
  std::vector<FFVar> vVarVal( nVar );
  for( unsigned i=0; i<nVar; ++i ) vVarVal[i] = vVar[i].val();

  FFVar** vResVal = insert_external_operation( *this, nRes,
                                                _nStaColl, vVarVal.data(),
                                                _nInpColl, vVarVal.data()+_nStaColl,
                                                _nCst,     vVarVal.data()+_nStaColl+_nInpColl );

  std::vector<size_t> nnz( nRes, 0 );
  if( !_pOCFESLV->deriv( nnz.data(), nullptr ) )
    throw std::runtime_error( "FFOCFERES::eval FADType<FFVar> ** Failure retrieving sparsity counts\n" );

  std::vector< std::vector<size_t> > col_store( nRes );
  std::vector<size_t*> col_ptr( nRes, nullptr );
  for( unsigned irow=0; irow<nRes; ++irow ){
    col_store[irow].resize( nnz[irow] );
    col_ptr[irow] = col_store[irow].empty()? nullptr: col_store[irow].data();
  }
  if( !_pOCFESLV->deriv( nnz.data(), col_ptr.data() ) )
    throw std::runtime_error( "FFOCFERES::eval FADType<FFVar> ** Failure retrieving sparsity columns\n" );

  FFGradOCFERES ResDer;
  int const gradPolicy = _ownOCFESLV ? FFGradOCFERES::COPY
                                  : FFGradOCFERES::SHALLOW;
  FFVar** vResDer = ResDer( _nStaColl, vVarVal.data(),
                            _nInpColl, vVarVal.data()+_nStaColl,
                            _nCst,     vVarVal.data()+_nStaColl+_nInpColl,
                            _pOCFESLV, gradPolicy, _name );

  size_t iel = 0;
  for( unsigned irow=0; irow<nRes; ++irow ){
    vRes[irow] = *vResVal[irow];

    bool dependent = false;
    for( size_t inz=0; inz<nnz[irow]; ++inz ){
      size_t const jcol = col_store[irow][inz];
      if( jcol >= nVar || !vVar[jcol].depend() ) continue;
      vRes[irow].setDepend( vVar[jcol] );
      dependent = true;
    }

    if( dependent ){
      for( unsigned ider=0; ider<vRes[irow].size(); ++ider ){
        vRes[irow][ider] = 0.;
        for( size_t inz=0; inz<nnz[irow]; ++inz ){
          size_t const jcol = col_store[irow][inz];
          if( jcol >= nVar || !vVar[jcol].depend() ) continue;
          vRes[irow][ider] += *vResDer[iel+inz] * vVar[jcol][ider];
        }
      }
    }

    iel += nnz[irow];
  }
  assert( iel == _nnzEqnColl + _nnzFctColl );
}

inline void
FFOCFERES::deriv
( unsigned const nRes, FFVar const* vRes, unsigned const nVar, FFVar const* vVar,
  FFVar** vDer )
const
{
#ifdef MC__FFDAGEXT_CHECK
  assert( _pOCFESLV && nRes == _nEqnColl+_nFctColl && nVar == _nStaColl+_nInpColl+_nCst );
#endif
  std::vector<size_t> nnz( nRes, 0 );
  std::vector< std::vector<size_t> > col_store( nRes );
  std::vector<size_t*> col_ptr( nRes, nullptr );

  deriv( nRes, vRes, nVar, vVar, nullptr, nnz.data(), nullptr );
  for( unsigned i=0; i<nRes; ++i ){
    col_store[i].resize( nnz[i] );
    col_ptr[i] = col_store[i].data();
  }

  std::vector< std::vector<FFVar> > grad_store( nRes );
  std::vector<FFVar*> grad_ptr( nRes, nullptr );
  for( unsigned i=0; i<nRes; ++i ){
    grad_store[i].resize( nnz[i] );
    grad_ptr[i] = grad_store[i].data();
  }

  deriv( nRes, vRes, nVar, vVar, grad_ptr.data(), nnz.data(), col_ptr.data() );

  for( unsigned i=0; i<nRes; ++i ){
    for( unsigned j=0; j<nVar; ++j ) vDer[i][j] = 0.;
    for( size_t k=0; k<nnz[i]; ++k )
      if( col_store[i][k] < nVar ) vDer[i][ col_store[i][k] ] = grad_store[i][k];
  }
}

inline void
FFOCFERES::deriv
( unsigned const nRes, FFVar const* vRes, unsigned const nVar, FFVar const* vVar,
  FFVar** vDer, size_t* nnz, size_t** colnz )
const
{
#ifdef MC__FFDAGEXT_TRACE
  std::cout << "FFDAGEXT::deriv (sparse)\n";
#endif
#ifdef MC__FFDAGEXT_CHECK
  assert( _pOCFESLV && nRes == _nEqnColl+_nFctColl && nVar == _nStaColl+_nInpColl+_nCst );
  assert( this->sparse && nnz && colnz );
#endif

  if( options.SYMDIFF.empty() ){
    // FFGraph sparse external differentiation expects this method to size the
    // per-output derivative arrays.  OCFESLV::deriv(nnz,colnz) only fills user-
    // allocated column arrays, so allocate missing row storage here.
    if( !_pOCFESLV->deriv( nnz, nullptr ) )
      throw std::runtime_error( "FFOCFERES::deriv ** Failure retrieving sparsity counts\n" );

    if( colnz ){
      for( size_t irow=0; irow<_nEqnColl+_nFctColl; ++irow )
        if( !colnz[irow] ) colnz[irow] = new size_t[nnz[irow]];
      if( !_pOCFESLV->deriv( nnz, colnz ) )
        throw std::runtime_error( "FFOCFERES::deriv ** Failure retrieving sparsity columns\n" );
    }

    if( !vDer ) return;

    for( size_t irow=0; irow<_nEqnColl+_nFctColl; ++irow )
      if( !vDer[irow] ) vDer[irow] = new FFVar[nnz[irow]];

    FFGradOCFERES ResDer;
    // Preserve the ownership policy of the primal collocation operation when
    // creating the derivative external operation.  If the primal operation
    // owns a private OCFESLV (COPY/TRANSFER), the FFGradOCFERES node must own
    // its own private executable OCFESLV too; otherwise derivative-DAG copies
    // remain shallow aliases of the source OCFESLV and parallel evaluation can
    // race through OCFESLV's mutable workspaces.
    int const gradPolicy = _ownOCFESLV ? FFGradOCFERES::COPY
                                    : FFGradOCFERES::SHALLOW;
    FFVar** vResDer = ResDer( _nStaColl, vVar, _nInpColl, vVar+_nStaColl,
                              _nCst, vVar+_nStaColl+_nInpColl,
                              _pOCFESLV, gradPolicy, _name );

    size_t iel = 0;
    for( size_t irow=0; irow<_nEqnColl+_nFctColl; ++irow )
      for( size_t jcol=0; jcol<nnz[irow]; ++jcol, ++iel )
        vDer[irow][jcol] = *vResDer[iel];
    assert( iel == _nnzEqnColl + _nnzFctColl );
    return;
  }

/*
    // Match DAG variables to IVP parameters
    std::vector<FFVar> vPar; 
    vPar.reserve( options.SYMDIFF.size() );
    std::vector<unsigned> ndxPar; 
    ndxPar.reserve( options.SYMDIFF.size() );
    for( auto const& dVar : options.SYMDIFF ){
      for( unsigned i=0; i<nVar; ++i ){
        if( dVar.id() != vVar[i].id() ) continue;
#ifdef CRONOS__FFOCFERES_CHECK
        std::cout << "Sensitivity parameter #" << i << ": " << dVar << std::endl;
#endif
        vPar.push_back( i<_nPar? _pODESLV->var_parameter()[i]: _pODESLV->var_constant()[i-_nPar] );
        ndxPar.push_back( i );
        break;
      }
    }

    FFOCFERES ResDer;
    auto pODESLVSEN = _pODESLV->fdiff( vPar );
    pODESLVSEN->setup();
    FFVar const*const* vResDer = ResDer._set( _nPar, vVar, _nCst, vVar+_nPar, pODESLVSEN, -1, _name ); // transfer ownership of ODE data to DAG
    for( unsigned k=0; k<nRes; ++k ){
      for( unsigned i=0; i<nVar; ++i )
        vDer[k][i] = 0;
      for( unsigned ii=0; ii<ndxPar.size(); ++ii )
        vDer[k][ndxPar[ii]] = *vResDer[k+nRes*ii];
    }
*/

  throw std::runtime_error( "FFOCFERES::deriv ** Symbolic differentiation not yet implementedFailure\n" );
}

inline void
FFGradOCFERES::eval
( unsigned const nRes, FFVar* vRes, unsigned const nVar, FFVar const* vVar,
  unsigned const* mVar )
const
{
#ifdef CRONOS__FFOCFERES_TRACE
  std::cout << "FFGradOCFERES::eval: FFVar\n";
#endif
#ifdef CRONOS__FFOCFERES_CHECK
  assert( _pOCFESLV && nRes == _nnzEqnColl+ _nnzFctColl && nVar == _nStaColl+_nInpColl+_nCst );
#endif

  FFVar** ppRes = nullptr;
  if( _into_other_dag( nVar, vVar ) ){
    auto op = *this;                  // shares the solver ...
    op._ownOCFESLV = true;            // ... so that the INSERTED copy deep-copies it (FFBaseOCFE copy constructor)
    ppRes = insert_external_operation( op, nRes, _nStaColl, vVar, _nInpColl, vVar+_nStaColl, _nCst, vVar+_nStaColl+_nInpColl );
    op._ownOCFESLV = false;           // the temporary must not delete the shared solver
  }
  else
    ppRes = insert_external_operation( *this, nRes, _nStaColl, vVar, _nInpColl, vVar+_nStaColl, _nCst, vVar+_nStaColl+_nInpColl );
  for( unsigned j=0; j<nRes; ++j )
    vRes[j] = *(ppRes[j]);
}
/*
inline void
FFGradOCFERES::eval
( unsigned const nRes, FFDep* vRes, unsigned const nVar, FFDep const* vVar,
  unsigned const* mVar )
const
{
#ifdef CRONOS__FFOCFERES_TRACE
  std::cout << "FFGradOCFERES::eval: FFDep\n";
#endif
#ifdef CRONOS__FFOCFERES_CHECK
  assert( _pOCFESLV && nRes == _nnzEqnColl+ _nnzFctColl && nVar == _nStaColl+_nInpColl+_nCst );
#endif

  vRes[0] = 0;
  for( unsigned i=0; i<nVar; ++i ) vRes[0] += vVar[i];
  vRes[0].update( FFDep::TYPE::N );
  for( unsigned j=1; j<nRes; ++j ) vRes[j] = vRes[0];
}
*/

inline void
FFGradOCFERES::eval
( unsigned const nRes, FFDep* vRes, unsigned const nVar, FFDep const* vVar,
  unsigned const* mVar )
const
{
#ifdef CRONOS__FFOCFE_CHECK
  assert( _pOCFESLV && nRes == _nnzEqnColl+_nnzFctColl && nVar == _nStaColl+_nInpColl+_nCst );
#endif
  FFDep dep;
  for( unsigned i=0; i<nVar; ++i ) dep += vVar[i];
  dep.update( FFDep::TYPE::N );
  for( unsigned j=0; j<nRes; ++j ) vRes[j] = dep;
}

inline void
FFGradOCFERES::eval
( unsigned const nRes, SLiftVar* vRes, unsigned const nVar, SLiftVar const* vVar,
  unsigned const* mVar )
const
{
#ifdef CRONOS__FFOCFERES_TRACE
  std::cout << "FFGradOCFERES::eval: SLiftVar\n";
#endif
#ifdef CRONOS__FFOCFERES_CHECK
  assert( _pOCFESLV && nRes == _nnzEqnColl+ _nnzFctColl && nVar == _nStaColl+_nInpColl+_nCst );
#endif

  vVar->env()->lift( nRes, vRes, nVar, vVar );
}

inline void
FFGradOCFERES::eval
( unsigned const nRes, FFExpr* vRes, unsigned const nVar, FFExpr const* vVar,
  unsigned const* mVar )
const
{
#ifdef CRONOS__FFOCFERES_TRACE
  std::cout << "FFGradOCFERES::eval: FFExpr\n";
#endif
#ifdef CRONOS__FFOCFERES_CHECK
  assert( _pOCFESLV && nRes == _nnzEqnColl+ _nnzFctColl && nVar == _nStaColl+_nInpColl+_nCst );
#endif

  switch( FFExpr::options.LANG ){
   case FFExpr::Options::DAG:
    for( unsigned j=0; j<nRes; ++j ){
      std::ostringstream os; os << name() << "[" << j << "]";
      vRes[j] = FFExpr::compose( os.str(), nVar, vVar );
    }
    break;
   case FFExpr::Options::GAMS:
   default:
    throw typename FFExpr::Exceptions( FFExpr::Exceptions::UNDEF );
  }
}

inline void
FFGradOCFERES::eval
( unsigned const nRes, double* vRes, unsigned const nVar, double const* vVar,
  unsigned const* mVar )
const
{
#ifdef CRONOS__FFOCFERES_TRACE
  std::cout << "FFGradOCFERES::eval: double\n";
#endif
#ifdef CRONOS__FFOCFERES_CHECK
  assert( _pOCFESLV && nRes == _nnzEqnColl+ _nnzFctColl && nVar == _nStaColl+_nInpColl+_nCst );
#endif

  if( !_pOCFESLV->deriv( vRes, vRes+_nnzEqnColl, vVar, vVar+_nStaColl, vVar+_nStaColl+_nInpColl ) )
    throw std::runtime_error( "FFGradOCFERES::eval double ** Failure\n" );
}


//! @brief C++ class defining REDUCED-SPACE IPDAE collocation as an external DAG operation in MC++.
////////////////////////////////////////////////////////////////////////
//! mc::FFOCFESLV is the reduced-space sibling of mc::FFOCFERES.  Where
//! FFOCFERES exposes the FULL collocation space (states + inputs) and
//! dispatches feval over the full residual + function stack, FFOCFESLV
//! solves the discretized state system INTERNALLY (marching backend) and
//! exposes ONLY the reduced control inputs -> output functions map:
//!
//!     F : R^{n_control_dof}  ->  R^{n_colloc_fct}
//!
//! so an outer NLP shrinks from O(state DOFs) to O(controls).  feval dispatch
//! mirrors FFOCFERES but with the reduced contract:
//!   double            -> decode_controls(p) -> solve (marching) -> functions
//!   FADType<double> -> ONE march_fsens_setup (primal + factorization),
//!                        then ONE march_fsens_apply per seed direction:
//!                        (dF/dp).dp  -- no Jacobian is ever formed
//!   BADType<double> -> reverse (adjoint) march  [NOT YET BUILT: throws]
//!   FADType<FFVar>  -> symbolic forward AD via FFGradOCFESLV (mirrors
//!                        FFOCFERES -> FFGradOCFERES)
//!   FFVar/FFDep/SLiftVar/FFExpr -> symbolic / dependency passthrough
//! The symbolic Jacobian dF/dp is DENSE (n_colloc_fct x n_control_dof): every
//! output function generically couples to every control through the internal
//! state solve.
//!
//! The state is never a DAG variable of this operation; the nominal input
//! collocation vector (_inp0, with the control-linked entries overwritten by
//! decode_controls at eval time) and the state initial guess (_xv0) are held
//! as numeric data on the operation.
////////////////////////////////////////////////////////////////////////

//! @brief Refuse if the solver's control registry is no longer the one the op was built with.
inline void
FFBaseOCFE::_check_controls
( OCFESLV const* oc, Layout const& L )
{
  bool same = ( oc->n_control_dof() == L.nCtrl && oc->controls().size() == L.ctrlRec.size() );
  if( same ){ size_t i = 0; for( auto const& [u,spec] : oc->controls() ) if( u.id() != L.ctrlRec[i++].id() ){ same = false; break; } }
  if( !same )
    throw std::runtime_error( "FFOCFESLV ** the OCFESLV control registry changed after this operation was built"
      " (SHALLOW policy): its inputs no longer match the solver's controls. Rebuild the operation, or embed with COPY.\n" );
}

//! @brief From the op inputs @p v (layout @p L): state starting values from the references (init), the input
//! vector from map 1 (decode_controls) and map 2 (set_input_values), and the constants.
inline void
FFBaseOCFE::_assemble
( OCFESLV* oc, Layout const& L, double const* v,
  std::vector<double>& xv, std::vector<double>& inp, std::vector<double>& cst )
{
  cst.assign( oc->var_constant().size(), 0. );
  size_t off = L.nCtrl + L.nRest();
  for( size_t i=0; i<L.cstNdx.size(); ++i ) cst[ L.cstNdx[i] ] = v[ off+i ];
  if( !oc->init( xv, inp, cst.empty()? nullptr: cst.data() ) )
    throw std::runtime_error( "FFOCFESLV ** OCFESLV::init failed (state references)\n" );
  std::vector<double> const pc( v, v+L.nCtrl );
  oc->decode_controls( pc, inp.data() );
  off = L.nCtrl;
  for( auto const& [w,nd] : L.rest ){
    std::vector<double> const blk( v+off, v+off+nd );
    oc->set_input_values( w, blk, inp.data() );
    off += nd;
  }
}

//! @brief Resolve the two maps, RESET the registry to map 1, and return the op inputs in layout order.
//! @p saved receives the previous registry (for COPY to restore).  Refuses an undeclared entry, a wrong DOF
//! count, an entry listed twice, a constant in map 1, and any declared input or constant missing from both.
inline std::vector<FFVar>
FFBaseOCFE::_bind
( std::vector<FFModel::InputArg> const& vDiff, std::vector<FFModel::InputArg> const& vRest,
  OCFESLV* oc, std::vector<FFVar>& saved, Layout& L )
{
  if( !oc ) throw std::runtime_error( "FFOCFESLV ** Null pointer to OCFESLV\n" );
  auto const& vC = oc->var_constant();
  auto cst_index = [&]( FFVar const& v )->int {
    for( size_t c=0; c<vC.size(); ++c ) if( vC[c].id().second == v.id().second ) return (int)c;
    return -1; };
  auto is_input = [&]( FFVar const& v ){
    for( auto const& [w,dom] : oc->var_declared_input() ) if( w.id().second == v.id().second ) return true;
    return false; };
  std::vector<FFVar> seen;
  auto once = [&]( FFVar const& v ){
    for( auto const& u : seen ) if( u.id().second == v.id().second )
      throw std::runtime_error( "FFOCFESLV ** " + v.name() + " appears more than once in the two maps\n" );
    seen.push_back( v ); };
  std::string err;
  std::vector<std::vector<FFVar>> dDiff( vDiff.size() ), dRest( vRest.size() );
  std::vector<int> cRest( vRest.size(), -1 );
  for( size_t i=0; i<vDiff.size(); ++i ){
    once( vDiff[i].input );
    if( cst_index( vDiff[i].input ) >= 0 )
      throw std::runtime_error( "FFOCFESLV ** constant " + vDiff[i].input.name() + " cannot be in map 1:"
        " the reduced sensitivities run over inputs only\n" );
    if( !oc->resolve_input_arg( vDiff[i], dDiff[i], err ) ) throw std::runtime_error( "FFOCFESLV ** " + err + "\n" );
  }
  for( size_t i=0; i<vRest.size(); ++i ){
    once( vRest[i].input );
    cRest[i] = cst_index( vRest[i].input );
    if( cRest[i] >= 0 ){
      if( vRest[i].gen ) dRest[i].assign( 1, vRest[i].gen( FFModel::DofIndex() ) );
      else               dRest[i] = vRest[i].dofs;
      if( dRest[i].size() != 1 )
        throw std::runtime_error( "FFOCFESLV ** constant " + vRest[i].input.name() + " takes exactly one variable\n" );
    }
    else if( !oc->resolve_input_arg( vRest[i], dRest[i], err ) ) throw std::runtime_error( "FFOCFESLV ** " + err + "\n" );
  }
  for( auto const& u : seen )
    if( !is_input( u ) && cst_index( u ) < 0 )
      throw std::runtime_error( "FFOCFESLV ** " + u.name() + " is neither a declared input nor a constant of the model\n" );
  std::string missing;
  for( auto const& [w,dom] : oc->var_declared_input() ){
    if( oc->fixed_input( w ) ) continue;                  // the model supplies its values (fix_input)
    bool in = false; for( auto const& u : seen ) if( u.id().second == w.id().second ){ in = true; break; }
    if( !in ) missing += " " + w.name();
  }
  for( auto const& c : vC ){
    bool in = false; for( auto const& u : seen ) if( u.id().second == c.id().second ){ in = true; break; }
    if( !in ) missing += " " + c.name();
  }
  if( !missing.empty() )
    throw std::runtime_error( "FFOCFESLV ** every declared input and constant must be mapped (map 1 or map 2) so its"
      " value can propagate through the DAG; missing:" + missing + "\n" );
  // the registry is map 1
  saved.clear();
  for( auto const& [u,spec] : oc->controls() ) saved.push_back( u );
  oc->clear_controls();
  for( auto const& a : vDiff ) oc->register_control( a.input );
  size_t nd1 = 0; for( auto const& d : dDiff ) nd1 += d.size();
  if( oc->n_control_dof() != nd1 ){
    oc->clear_controls(); for( auto const& u : saved ) oc->register_control( u );
    throw std::runtime_error( "FFOCFESLV ** map 1 does not match the solver's control layout\n" );
  }
  L = Layout();
  L.nCtrl = nd1;
  for( auto const& [u,spec] : oc->controls() ) L.ctrlRec.push_back( u );
  std::vector<FFVar> vIn( nd1 );
  for( size_t i=0; i<vDiff.size(); ++i ){
    auto const b = oc->control_block( vDiff[i].input );
    for( size_t k=0; k<b.ndof; ++k ) vIn[ b.offset+k ] = dDiff[i][k];
  }
  for( size_t i=0; i<vRest.size(); ++i ) if( cRest[i] < 0 ){
    L.rest.emplace_back( vRest[i].input, dRest[i].size() );
    vIn.insert( vIn.end(), dRest[i].begin(), dRest[i].end() );
  }
  for( size_t i=0; i<vRest.size(); ++i ) if( cRest[i] >= 0 ){
    L.cstNdx.push_back( (size_t)cRest[i] );
    vIn.push_back( dRest[i][0] );
  }
  return vIn;
}

class FFOCFESLV
: public FFBaseOCFE
{

protected:

  //! @brief Number of map-1 DOFs (the reduced controls; the Jacobian's columns)
  size_t              _nCtrl;
  //! @brief The op-input layout (two-map contract)
  Layout    _L;
  std::vector<size_t> _layout_key() const override { return _layout_key_of( _L ); }

  FFVar** _set
    ( size_t const nIn, FFVar const* pIn, Layout const& L,
      OCFESLV* pOCFESLV, int policy, std::string const& name )
    {
      if( !pOCFESLV || nIn != L.nIn() || L.nCtrl != pOCFESLV->n_control_dof() )
        throw std::runtime_error( "FFOCFESLV ** inconsistent operation layout\n" );
      if( _ownOCFESLV && _pOCFESLV )
        delete _pOCFESLV;
      _ownOCFESLV   = ( policy>0? true: false ); //copy;
      _pOCFESLV     = pOCFESLV;
      _nCst       = 0;
      _nStaColl   = 0;
      _nInpColl   = pOCFESLV->n_colloc_inp();
      _nEqnColl   = pOCFESLV->n_colloc_eqn();
      _nFctColl   = pOCFESLV->n_colloc_fct();
      _nnzEqnColl = 0;
      _nnzFctColl = 0;
      _nCtrl      = L.nCtrl;
      _L          = L;
      _name       = name;

      data    = pOCFESLV;
      owndata = false;
      FFVar** ppRes = insert_external_operation( *this, _nFctColl, nIn, pIn );

      _ownOCFESLV = false;
      FFOp* pOp = (*ppRes)->opdef().first;
      if( policy > 0 )
        _pOCFESLV = dynamic_cast<FFOCFESLV*>(pOp)->_pOCFESLV; // set pointer to DAG copy
      else if( policy < 0 )
        dynamic_cast<FFOCFESLV*>(pOp)->_ownOCFESLV = true;   // transfer ownership
      return ppRes;
    }

  //! @brief Bind the two maps, insert the op, and -- for COPY -- restore the caller's registry.
  FFVar** _embed
    ( std::vector<FFModel::InputArg> const& vDiff, std::vector<FFModel::InputArg> const& vRest,
      OCFESLV* pOCFESLV, int policy, std::string const& name )
    {
      std::vector<FFVar> saved;  Layout L;
      std::vector<FFVar> const vIn = _bind( vDiff, vRest, pOCFESLV, saved, L );
      FFVar** ppRes = _set( vIn.size(), vIn.data(), L, pOCFESLV, policy, name );
      if( policy > 0 ){                                 // the op's copy keeps the registry
        pOCFESLV->clear_controls();
        for( auto const& u : saved ) pOCFESLV->register_control( u );
      }
      return ppRes;
    }

public:

  //! @brief Enumeration type for OCFESLV copy policy
  enum POLICY_TYPE{
    SHALLOW=0,  //!< Shallow copy of OCFESLV in FFGraph (without ownership)
    COPY=1,     //!< Deep copy of OCFESLV in FFGraph (with ownership)
    TRANSFER=-1 //!< Shallow copy of OCFESLV in FFGraph (with ownership transfer)
  };

  // Default constructor
  FFOCFESLV
    ()
    : FFBaseOCFE(),
      _nCtrl( 0 )
    {
      this->sparse = true;
    }

  // Copy constructor
  FFOCFESLV
    ( FFOCFESLV const& Other )
    : FFBaseOCFE( Other ),
      _nCtrl( Other._nCtrl ),
      _L    ( Other._L )
    {}

  // Destructor
  virtual ~FFOCFESLV
    ()
    {}

  //! @brief TWO-MAP form.  @p vDiff (map 1) lists the inputs the NUMERICAL gradient is taken w.r.t.: it RESETS the
  //! OCFESLV control registry to exactly these.  @p vRest (map 2) lists every other declared input and every
  //! constant: their values pass through the DAG, with ZERO numerical-derivative columns by contract.  Each entry
  //! is a flat vector in control_dofs() order or a generator called per DofIndex; a constant takes one variable.
  //! Every declared input and every constant must appear in exactly one map.  The state starting values come
  //! from the references (OCFESLV::init) at every evaluation.  COPY leaves the caller's registry as it was;
  //! SHALLOW and TRANSFER keep the new one.  Op inputs: [ map 1, control order | map-2 inputs | map-2 constants ],
  //! map 2 in listing order.  Returns the n_colloc_fct outputs.
  std::vector<FFVar> operator()
    ( std::vector<FFModel::InputArg> const& vDiff, std::vector<FFModel::InputArg> const& vRest,
      OCFESLV* pOCFESLV, int policy=SHALLOW, std::string const& name="" )
    {
      FFVar** ppRes = _embed( vDiff, vRest, pOCFESLV, policy, name );
      std::vector<FFVar> vDep( _nFctColl );
      for( size_t i=0; i<vDep.size(); ++i ) vDep[i] = *ppRes[i];
      return vDep;
    }

  //! @brief TWO-MAP form, output @p iDep only.
  FFVar& operator()
    ( size_t const iDep, std::vector<FFModel::InputArg> const& vDiff, std::vector<FFModel::InputArg> const& vRest,
      OCFESLV* pOCFESLV, int policy=SHALLOW, std::string const& name="" )
    {
      FFVar** ppRes = _embed( vDiff, vRest, pOCFESLV, policy, name );
      if( iDep >= _nFctColl ) throw std::runtime_error( "FFOCFESLV ** output index out of range\n" );
      return *ppRes[iDep];
    }

  //! @brief ONE-MAP form: every declared input is differentiated and the model has no constant -- the two-map
  //! form with an empty map 2, refused (naming what is missing) otherwise.
  std::vector<FFVar> operator()
    ( std::vector<FFModel::InputArg> const& vDiff, OCFESLV* pOCFESLV, int policy=SHALLOW, std::string const& name="" )
    { return (*this)( vDiff, std::vector<FFModel::InputArg>(), pOCFESLV, policy, name ); }

  //! @brief ONE-MAP form, output @p iDep only.
  FFVar& operator()
    ( size_t const iDep, std::vector<FFModel::InputArg> const& vDiff, OCFESLV* pOCFESLV, int policy=SHALLOW,
      std::string const& name="" )
    { return (*this)( iDep, vDiff, std::vector<FFModel::InputArg>(), pOCFESLV, policy, name ); }

  //! @brief Map 1 built from the solver's CURRENT control registry, with @p pCtrl in control order
  //! (pCtrl[ control_block(u).offset + k ] for DOF k of control u) -- for callers that keep control vectors.
  static std::vector<FFModel::InputArg> controls_map
    ( OCFESLV const& oc, std::vector<FFVar> const& pCtrl )
    {
      std::vector<FFModel::InputArg> m;
      for( auto const& [u,spec] : oc.controls() )
        m.emplace_back( u, std::vector<FFVar>( pCtrl.begin()+spec.offset, pCtrl.begin()+spec.offset+spec.ndof ) );
      return m;
    }

  //! @brief The direction of op input @p i: a (control, DOF), a (map-2 input, DOF), or a (constant, 0).
  std::pair<FFVar,size_t> _op_input_direction
    ( size_t const i )
    const
    {
      if( i < _nCtrl ){
        for( auto const& [u,spec] : _pOCFESLV->controls() )
          if( i >= spec.offset && i < spec.offset + spec.ndof ) return { u, i - spec.offset };
      }
      else if( i < _nCtrl + _L.nRest() ){
        size_t off = _nCtrl;
        for( auto const& [w,nd] : _L.rest ){ if( i < off + nd ) return { w, i - off }; off += nd; }
      }
      else{
        size_t const j = i - _nCtrl - _L.nRest();
        if( j < _L.cstNdx.size() ) return { _pOCFESLV->var_constant()[ _L.cstNdx[j] ], 0 };
      }
      throw std::runtime_error( "FFOCFESLV ** op input index out of range\n" );
    }

  //! @brief Number of op inputs (map 1 + map-2 inputs + constants)
  size_t n_input
    ()
    const
    { return _L.nIn(); }

  OCFESLV* pOCFESLV
    ()
    const
    { return _pOCFESLV; }

  //! @brief Number of reduced control DOFs (external-operation inputs)
  size_t n_control
    ()
    const
    { return _nCtrl; }

  // Evaluation overloads
  virtual void feval
    ( std::type_info const& idU, unsigned const nRes, void* vRes, unsigned const nVar,
      void const* vVar, unsigned const* mVar )
    const
    {
      if( mc::same_type( idU, typeid( FFVar ) ) )
        return eval( nRes, static_cast<FFVar*>(vRes), nVar, static_cast<FFVar const*>(vVar), mVar );
      else if( mc::same_type( idU, typeid( FADType<FFVar> ) ) )
        return eval( nRes, static_cast<FADType<FFVar>*>(vRes), nVar, static_cast<FADType<FFVar> const*>(vVar), mVar );
      else if( mc::same_type( idU, typeid( FADType<double> ) ) )
        return eval( nRes, static_cast<FADType<double>*>(vRes), nVar, static_cast<FADType<double> const*>(vVar), mVar );
      else if( mc::same_type( idU, typeid( BADType<double> ) ) )
        return eval( nRes, static_cast<BADType<double>*>(vRes), nVar, static_cast<BADType<double> const*>(vVar), mVar );
      else if( mc::same_type( idU, typeid( FFDep ) ) )
        return eval( nRes, static_cast<FFDep*>(vRes), nVar, static_cast<FFDep const*>(vVar), mVar );
      else if( mc::same_type( idU, typeid( double ) ) )
        return eval( nRes, static_cast<double*>(vRes), nVar, static_cast<double const*>(vVar), mVar );
      else if( mc::same_type( idU, typeid( SLiftVar ) ) )
        return eval( nRes, static_cast<SLiftVar*>(vRes), nVar, static_cast<SLiftVar const*>(vVar), mVar );
      else if( mc::same_type( idU, typeid( FFExpr ) ) )
        return eval( nRes, static_cast<FFExpr*>(vRes), nVar, static_cast<FFExpr const*>(vVar), mVar );

      throw std::runtime_error( "FFOCFESLV::feval ** No evaluation method for type"+std::string(idU.name())+"\n" );
    }

  void eval
    ( unsigned const nRes, double* vRes, unsigned const nVar, double const* vVar, unsigned const* mVar )
    const;

  void eval
    ( unsigned const nRes, FADType<double>* vRes, unsigned const nVar, FADType<double> const* vVar,
      unsigned const* mVar )
    const;

  void eval
    ( unsigned const nRes, BADType<double>* vRes, unsigned const nVar, BADType<double> const* vVar,
      unsigned const* mVar )
    const;

  void eval
    ( unsigned const nRes, FADType<FFVar>* vRes, unsigned const nVar, FADType<FFVar> const* vVar,
      unsigned const* mVar )
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

  void eval
    ( unsigned const nRes, FFExpr* vRes, unsigned const nVar, FFExpr const* vVar, unsigned const* mVar )
    const;

  // Derivatives (symbolic): dense reduced Jacobian dF/dp via FFGradOCFESLV
  void deriv
    ( unsigned const nRes, FFVar const* vRes, unsigned const nVar, FFVar const* vVar, FFVar** vDer )
    const;

  void deriv
    ( unsigned const nRes, FFVar const* vRes, unsigned const nVar, FFVar const* vVar, FFVar** vDer,
      size_t* nnz, size_t** colnz )
    const;

  // Properties
  std::string name
    ()
    const
    { std::ostringstream oss;
      if( !_name.empty() )  oss << _name;
      else                  oss << _pOCFESLV;
      return "OCFE_SLV[" + oss.str() + "]"; }

  //! @brief Return whether or not operation is commutative
  bool commutative
    ()
    const
    { return false; }
};

//! @brief C++ class defining REDUCED-SPACE IPDAE collocation gradient evaluation as an external DAG operation.
////////////////////////////////////////////////////////////////////////
//! mc::FFGradOCFESLV is the reduced-space sibling of mc::FFGradOCFERES.
//! It exposes the DENSE reduced Jacobian dF/dp of the FFOCFESLV map,
//! flattened ROW-MAJOR: output at index irow*n_control_dof + jcol holds
//! d F_irow / d p_jcol,  with irow over 0..n_colloc_fct-1 and jcol over 0..n_control_dof-1.
//! eval<double> assembles it from ONE march_fsens_setup followed by ONE
//! march_fsens_apply per control column (one-hot direction).
////////////////////////////////////////////////////////////////////////
class FFGradOCFESLV
: public FFBaseOCFE
{
  friend class FFOCFESLV;               // built only by FFOCFESLV::deriv / eval( FADType<FFVar> ), via _set


protected:

  //! @brief Number of reduced control DOFs (= number of external-operation inputs)
  size_t              _nCtrl;
  //! @brief The op-input layout (as the FFOCFESLV it differentiates)
  Layout    _L;
  std::vector<size_t> _layout_key() const override { return _layout_key_of( _L ); }

  FFVar** _set
    ( size_t const nIn, FFVar const* pIn, Layout const& L,
      OCFESLV* pOCFESLV, int policy, std::string const& name )
    {
      if( !pOCFESLV || nIn != L.nIn() ) throw std::runtime_error( "FFGradOCFESLV ** inconsistent operation layout\n" );
      if( _ownOCFESLV && _pOCFESLV )
        delete _pOCFESLV;
      _ownOCFESLV   = ( policy>0? true: false ); //copy;
      _pOCFESLV     = pOCFESLV;
      _nCst       = 0;
      _nStaColl   = 0;
      _nInpColl   = pOCFESLV->n_colloc_inp();
      _nEqnColl   = pOCFESLV->n_colloc_eqn();
      _nFctColl   = pOCFESLV->n_colloc_fct();
      _nnzEqnColl = 0;
      _nnzFctColl = 0;
      _nCtrl      = L.nCtrl;
      _L          = L;
      _name       = name;

      data    = pOCFESLV;
      owndata = false;
      size_t const nDep = _nFctColl * _nCtrl;   // dense reduced Jacobian, row-major
      FFVar** ppRes = insert_external_operation( *this, nDep, nIn, pIn );

      _ownOCFESLV = false;
      FFOp* pOp = (*ppRes)->opdef().first;
      if( policy > 0 )
        _pOCFESLV = dynamic_cast<FFGradOCFESLV*>(pOp)->_pOCFESLV; // set pointer to DAG copy
#ifdef CRONOS__FFOCFE_TRACE
      std::cerr << "FFGradOCFESLV operation address: " << this << std::endl;
      std::cerr << "OCFESLV address in DAG: " << _pOCFESLV << std::endl;
#endif
      return ppRes;
    }

public:

  //! @brief Enumeration type for OCFESLV copy policy
  enum POLICY_TYPE{
    SHALLOW=0,  //!< Shallow copy of OCFESLV in FFGraph (without ownership)
    COPY=1      //!< Deep copy of OCFESLV in FFGraph (with ownership)
  };

  // Default constructor
  FFGradOCFESLV
    ()
    : FFBaseOCFE(),
      _nCtrl( 0 )
    {
      this->sparse = true;
    }

  // Copy constructor
  FFGradOCFESLV
    ( FFGradOCFESLV const& Other )
    : FFBaseOCFE( Other ),
      _nCtrl( Other._nCtrl ),
      _L    ( Other._L )
    {}

  // Destructor
  virtual ~FFGradOCFESLV
    ()
    {}

  // Define operation: controls -> flattened dense Jacobian dF/dp
  // No public call operator: built by FFOCFESLV::deriv / eval( FADType<FFVar> ) through _set(), with the
  // same op-input layout as the FFOCFESLV it differentiates.

  // Evaluation overloads
  virtual void feval
    ( std::type_info const& idU, unsigned const nRes, void* vRes, unsigned const nVar,
      void const* vVar, unsigned const* mVar )
    const
    {
      if( mc::same_type( idU, typeid( FFVar ) ) )
        return eval( nRes, static_cast<FFVar*>(vRes), nVar, static_cast<FFVar const*>(vVar), mVar );
      else if( mc::same_type( idU, typeid( FFDep ) ) )
        return eval( nRes, static_cast<FFDep*>(vRes), nVar, static_cast<FFDep const*>(vVar), mVar );
      else if( mc::same_type( idU, typeid( double ) ) )
        return eval( nRes, static_cast<double*>(vRes), nVar, static_cast<double const*>(vVar), mVar );
      else if( mc::same_type( idU, typeid( SLiftVar ) ) )
        return eval( nRes, static_cast<SLiftVar*>(vRes), nVar, static_cast<SLiftVar const*>(vVar), mVar );
      else if( mc::same_type( idU, typeid( FFExpr ) ) )
        return eval( nRes, static_cast<FFExpr*>(vRes), nVar, static_cast<FFExpr const*>(vVar), mVar );

      throw std::runtime_error( "FFGradOCFESLV::feval ** No evaluation method for type"+std::string(idU.name())+"\n" );
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

  void eval
    ( unsigned const nRes, FFExpr* vRes, unsigned const nVar, FFExpr const* vVar, unsigned const* mVar )
    const;

  // Properties
  std::string name
    ()
    const
    { std::ostringstream oss;
      if( !_name.empty() )  oss << _name;
      else                  oss << _pOCFESLV;
      return "Grad_OCFE_SLV[" + oss.str() + "]"; }

  //! @brief Return whether or not operation is commutative
  bool commutative
    ()
    const
    { return false; }
};

//======================= reduced-space Jacobian assembly (mode-aware) =======================

//! @brief Assemble the DENSE reduced Jacobian dF/dp (row-major, n_colloc_fct x n_control_dof),
//! choosing forward vs adjoint per @p mode (FFBaseOCFE::FORWARD / ADJOINT / AUTO with the
//! NP2NF ratio).  The complete OCFESLV::solve_fsens()/solve_asens() do the primal solve, the per-DOF
//! forward apply / per-function reverse sweep, and the control-space decode/encode internally, so
//! this is now just a mode dispatch + read-back of sens_jacobian().  Returns false on march failure.
inline bool
FFBaseOCFE::_reduced_jacobian
( OCFESLV* pOCFESLV, Layout const& L, double const* v, size_t const nFct, int mode, double np2nf,
  std::vector<double>& Jflat /* size nFct*L.nCtrl, row-major (irow*nCtrl+jcol) */ )
{
  size_t const nCtrl = L.nCtrl;
  bool const adjoint = ( mode == FFBaseOCFE::ADJOINT )
                    || ( mode == FFBaseOCFE::AUTO && (double)nCtrl > np2nf * (double)nFct );
  _check_controls( pOCFESLV, L );
  std::vector<double> xv, inp, cst;
  _assemble( pOCFESLV, L, v, xv, inp, cst );
  double const* pc = cst.empty()? nullptr: cst.data();
  bool const ok = adjoint ? pOCFESLV->solve_asens( xv.data(), inp.data(), pc )
                          : pOCFESLV->solve_fsens( xv.data(), inp.data(), pc );
  if( !ok ) return false;
  Jflat = pOCFESLV->sens_jacobian();                // n_colloc_fct x n_control_dof, row-major, control space
  return Jflat.size() >= nFct*nCtrl;
}

//======================= FFOCFESLV inline implementations =======================

inline void
FFOCFESLV::eval
( unsigned const nRes, FFVar* vRes, unsigned const nVar, FFVar const* vVar,
  unsigned const* mVar )
const
{
#ifdef CRONOS__FFOCFE_TRACE
  std::cout << "FFOCFESLV::eval: FFVar\n";
#endif
#ifdef CRONOS__FFOCFE_CHECK
  assert( _pOCFESLV && nRes == _nFctColl && nVar == _L.nIn() );
#endif
  (void)mVar;

  FFVar** ppRes = nullptr;
  if( _into_other_dag( nVar, vVar ) ){
    auto op = *this;                  // shares the solver ...
    op._ownOCFESLV = true;            // ... so that the INSERTED copy deep-copies it (FFBaseOCFE copy constructor)
    ppRes = insert_external_operation( op, nRes, nVar, vVar );
    op._ownOCFESLV = false;           // the temporary must not delete the shared solver
  }
  else
    ppRes = insert_external_operation( *this, nRes, nVar, vVar );
  for( unsigned j=0; j<nRes; ++j )
    vRes[j] = *(ppRes[j]);
}

inline void
FFOCFESLV::eval
( unsigned const nRes, FFDep* vRes, unsigned const nVar, FFDep const* vVar,
  unsigned const* mVar )
const
{
#ifdef CRONOS__FFOCFE_CHECK
  assert( _pOCFESLV && nRes == _nFctColl && nVar == _L.nIn() );
#endif
  (void)mVar;

  // Declare dependency from the structural reduced pattern: output f depends only on the control
  // DOFs P[f][j] (dense for a coupled monolithic block / marching; sparse when decomposable).
  std::vector<std::vector<bool>> const& P = _pOCFESLV->reduced_dependency_pattern();
  for( unsigned j=0; j<nRes; ++j ){
    FFDep dep;
    bool any = false;
    for( unsigned i=0; i<nVar; ++i )
      if( j < P.size() && i < P[j].size() && P[j][i] ){ dep += vVar[i]; any = true; }
    if( any ) dep.update( FFDep::TYPE::N );
    vRes[j] = dep;
  }
}

inline void
FFOCFESLV::eval
( unsigned const nRes, SLiftVar* vRes, unsigned const nVar, SLiftVar const* vVar,
  unsigned const* mVar )
const
{
#ifdef CRONOS__FFOCFE_CHECK
  assert( _pOCFESLV && nRes == _nFctColl && nVar == _L.nIn() );
#endif
  (void)mVar;

  vVar->env()->lift( nRes, vRes, nVar, vVar );
}

inline void
FFOCFESLV::eval
( unsigned const nRes, FFExpr* vRes, unsigned const nVar, FFExpr const* vVar,
  unsigned const* mVar )
const
{
#ifdef CRONOS__FFOCFE_CHECK
  assert( _pOCFESLV && nRes == _nFctColl && nVar == _L.nIn() );
#endif
  (void)mVar;

  switch( FFExpr::options.LANG ){
   case FFExpr::Options::DAG:
    for( unsigned j=0; j<nRes; ++j ){
      std::ostringstream os; os << name() << "[" << j << "]";
      vRes[j] = FFExpr::compose( os.str(), nVar, vVar );
    }
    break;
   case FFExpr::Options::GAMS:
   default:
    throw typename FFExpr::Exceptions( FFExpr::Exceptions::UNDEF );
  }
}

inline void
FFOCFESLV::eval
( unsigned const nRes, double* vRes, unsigned const nVar, double const* vVar,
  unsigned const* mVar )
const
{
#ifdef CRONOS__FFOCFE_TRACE
  std::cerr << "FFOCFESLV::eval: double w/ OCFESLV address " << _pOCFESLV << std::endl;
#endif
#ifdef CRONOS__FFOCFE_CHECK
  assert( _pOCFESLV && nRes == _nFctColl && nVar == _L.nIn() );
#endif
  (void)mVar;

  // Reduced input -> output map, value only: decode controls into the nominal input
  // vector, solve the state system internally (marching or monolithic), read functions.
  _check_controls( _pOCFESLV, _L );
  std::vector<double> xv, inp, cst;
  _assemble( _pOCFESLV, _L, vVar, xv, inp, cst );   // states from init; inputs from the two maps

  OCFESLV::SolveReport const rep = _pOCFESLV->solve( xv.data(), inp.data(), cst.empty()? nullptr: cst.data() );
  if( !rep.converged )
    throw std::runtime_error( "FFOCFESLV::eval double ** internal solve did not converge\n" );

  // solve() populates val_functions() for both marching (per-window accumulation) and
  // monolithic (function eval at the converged state).
  std::vector<double> const& F = _pOCFESLV->val_functions();
  if( F.size() < nRes )
    throw std::runtime_error( "FFOCFESLV::eval double ** function count mismatch\n" );
  for( unsigned j=0; j<nRes; ++j ) vRes[j] = F[j];
}

inline void
FFOCFESLV::eval
( unsigned const nRes, FADType<double>* vRes, unsigned const nVar, FADType<double> const* vVar,
  unsigned const* mVar )
const
{
#ifdef CRONOS__FFOCFE_TRACE
  std::cout << "FFOCFESLV::eval: FADType<double> via forward sensitivity (march or monolithic)\n";
#endif
#ifdef CRONOS__FFOCFE_CHECK
  assert( _pOCFESLV && nRes == _nFctColl && nVar == _L.nIn() );
#endif
  (void)mVar;

  // Mode-aware reduced Jacobian dF/dp (forward or adjoint per options.GRADIENT), computed by the
  // complete OCFESLV solve; sens_functions() then holds the primal values at the current controls.
  std::vector<double> p( nVar );
  for( unsigned i=0; i<nVar; ++i ) p[i] = vVar[i].val();
  std::vector<double> Jflat;
  if( !_reduced_jacobian( _pOCFESLV, _L, p.data(), nRes,
                                     options.GRADIENT, options.NP2NF, Jflat ) )
    throw std::runtime_error( "FFOCFESLV::eval FADType<double> ** reduced-Jacobian solve failed\n" );

  std::vector<double> F = _pOCFESLV->sens_functions();
  if( F.size() < nRes )
    throw std::runtime_error( "FFOCFESLV::eval FADType<double> ** function count mismatch\n" );
  for( unsigned k=0; k<nRes; ++k ) vRes[k] = F[k];   // sets value, resets derivatives

  // Number of seed directions (with consistency check across dependent inputs).
  bool dependent = false;
  unsigned nDir = 0;
  for( unsigned i=0; i<nVar; ++i ){
    if( !vVar[i].depend() ) continue;
    if( !dependent ) nDir = vVar[i].size();
    else if( vVar[i].size() != nDir )
      throw std::runtime_error( "FFOCFESLV::eval FADType<double> ** Inconsistent derivative vector sizes\n" );
    dependent = true;
  }
  if( !dependent ) return;

  // Mark the (dense) dependency: every output depends on every dependent control.
  for( unsigned irow=0; irow<nRes; ++irow )
    for( unsigned jcol=0; jcol<_nCtrl; ++jcol )           // map 2 / constants: zero by contract
      if( vVar[jcol].depend() ) vRes[irow].setDepend( vVar[jcol] );

  // Propagate each seed direction through the precomputed Jacobian: dF_i = sum_j (dF_i/dp_j) dp_j.
  for( unsigned idir=0; idir<nDir; ++idir ){
    for( unsigned irow=0; irow<nRes; ++irow ){
      size_t const base = size_t(irow) * _nCtrl;   // Jflat: nFct x nCtrl
      double s = 0.;
      for( unsigned jcol=0; jcol<_nCtrl; ++jcol )
        if( vVar[jcol].depend() ) s += Jflat[ base + jcol ] * vVar[jcol][idir];
      vRes[irow][idir] = s;
    }
  }
}

inline void
FFOCFESLV::eval
( unsigned const nRes, BADType<double>* vRes, unsigned const nVar, BADType<double> const* vVar,
  unsigned const* mVar )
const
{
#ifdef CRONOS__FFOCFE_TRACE
  std::cout << "FFOCFESLV::eval: BADType<double> via adjoint (reverse) march\n";
#endif
#ifdef CRONOS__FFOCFE_CHECK
  assert( _pOCFESLV && nRes == _nFctColl && nVar == _L.nIn() );
#endif
  (void)mVar;

  // Reduced Jacobian by the ADJOINT solve (reverse-mode's natural fit: one sweep per output);
  // the same solve leaves sens_functions() holding the primal values, so no separate value solve.
  std::vector<double> p( nVar );
  for( unsigned i=0; i<nVar; ++i ) p[i] = vVar[i].val();
  std::vector<double> Jflat;
  if( !_reduced_jacobian( _pOCFESLV, _L, p.data(), nRes,
                                     FFBaseOCFE::ADJOINT, options.NP2NF, Jflat ) )
    throw std::runtime_error( "FFOCFESLV::eval BADType<double> ** adjoint solve failed\n" );
  std::vector<double> Fval = _pOCFESLV->sens_functions();
  if( Fval.size() < nRes )
    throw std::runtime_error( "FFOCFESLV::eval BADType<double> ** function count mismatch\n" );

  // Attach the linear sensitivity so fadbad's reverse sweep propagates output adjoints
  // to the control inputs: vRes[i] = F_i + sum_j (dF_i/dp_j) ( p_j - p_j^value ).
  for( unsigned irow=0; irow<nRes; ++irow ){
    vRes[irow] = Fval[irow];
    size_t const base = size_t(irow) * _nCtrl;   // Jflat: nFct x nCtrl
    for( unsigned jcol=0; jcol<_nCtrl; ++jcol ){
      double const jv = Jflat[ base + jcol ];
      if( jv != 0. ) vRes[irow] += jv * ( vVar[jcol] - p[jcol] );
    }
  }
}

inline void
FFOCFESLV::eval
( unsigned const nRes, FADType<FFVar>* vRes, unsigned const nVar, FADType<FFVar> const* vVar,
  unsigned const* mVar )
const
{
#ifdef CRONOS__FFOCFE_TRACE
  std::cout << "FFOCFESLV::eval: FADType<FFVar> via FFGradOCFESLV\n";
#endif
#ifdef CRONOS__FFOCFE_CHECK
  assert( _pOCFESLV && nRes == _nFctColl && nVar == _L.nIn() );
#endif
  (void)mVar;

  std::vector<FFVar> vVarVal( nVar );
  for( unsigned i=0; i<nVar; ++i ) vVarVal[i] = vVar[i].val();

  // Primal (value) subgraph.
  FFVar** vResVal = insert_external_operation( *this, nRes, nVar, vVarVal.data() );

  // Dense reduced Jacobian subgraph (row-major, n_colloc_fct x n_control_dof).
  FFGradOCFESLV ResDer;
  int const gradPolicy = _ownOCFESLV ? FFGradOCFESLV::COPY
                                   : FFGradOCFESLV::SHALLOW;
  FFVar** vResDer = ResDer._set( nVar, vVarVal.data(), _L, _pOCFESLV, gradPolicy, _name );

  for( unsigned irow=0; irow<nRes; ++irow ){
    vRes[irow] = *vResVal[irow];
    size_t const base = size_t(irow) * _nCtrl;   // Jflat: nFct x nCtrl

    bool dependent = false;
    for( unsigned jcol=0; jcol<_nCtrl; ++jcol ){
      if( !vVar[jcol].depend() ) continue;
      vRes[irow].setDepend( vVar[jcol] );
      dependent = true;
    }

    if( dependent ){
      for( unsigned ider=0; ider<vRes[irow].size(); ++ider ){
        vRes[irow][ider] = 0.;
        for( unsigned jcol=0; jcol<_nCtrl; ++jcol ){
          if( !vVar[jcol].depend() ) continue;
          vRes[irow][ider] += *vResDer[base + jcol] * vVar[jcol][ider];
        }
      }
    }
  }
}

inline void
FFOCFESLV::deriv
( unsigned const nRes, FFVar const* vRes, unsigned const nVar, FFVar const* vVar,
  FFVar** vDer )
const
{
#ifdef CRONOS__FFOCFE_CHECK
  assert( _pOCFESLV && nRes == _nFctColl && nVar == _L.nIn() );
#endif
  (void)vRes;

  // Dense scatter of the (here fully dense) reduced Jacobian.
  std::vector<size_t> nnz( nRes, 0 );
  std::vector< std::vector<size_t> > col_store( nRes );
  std::vector<size_t*> col_ptr( nRes, nullptr );

  deriv( nRes, vRes, nVar, vVar, nullptr, nnz.data(), nullptr );
  for( unsigned i=0; i<nRes; ++i ){
    col_store[i].resize( nnz[i] );
    col_ptr[i] = col_store[i].empty()? nullptr: col_store[i].data();
  }

  std::vector< std::vector<FFVar> > grad_store( nRes );
  std::vector<FFVar*> grad_ptr( nRes, nullptr );
  for( unsigned i=0; i<nRes; ++i ){
    grad_store[i].resize( nnz[i] );
    grad_ptr[i] = grad_store[i].empty()? nullptr: grad_store[i].data();
  }

  deriv( nRes, vRes, nVar, vVar, grad_ptr.data(), nnz.data(), col_ptr.data() );

  for( unsigned i=0; i<nRes; ++i ){
    for( unsigned j=0; j<nVar; ++j ) vDer[i][j] = 0.;
    for( size_t k=0; k<nnz[i]; ++k )
      if( col_store[i][k] < nVar ) vDer[i][ col_store[i][k] ] = grad_store[i][k];
  }
}

inline void
FFOCFESLV::deriv
( unsigned const nRes, FFVar const* vRes, unsigned const nVar, FFVar const* vVar,
  FFVar** vDer, size_t* nnz, size_t** colnz )
const
{
#ifdef CRONOS__FFOCFE_CHECK
  assert( _pOCFESLV && nRes == _nFctColl && nVar == _L.nIn() && this->sparse && nnz );
#endif
  (void)vRes;

  if( !options.SYMDIFF.empty() ){
    // SYMDIFF: a NEW FFOCFESLV on the forward-differentiated model (OCFESLV::fdiff) along the listed op inputs --
    // control DOFs, map-2 DOFs or constants alike -- with the same op-input layout; unlisted columns stay ZERO
    // (the contract downstream libraries rely on).  Its outputs are f[i*nf+j] = dF_j/d(listed input i).
    std::vector<unsigned> cols;
    for( unsigned i=0; i<nVar; ++i )
      for( auto const& v : options.SYMDIFF ) if( v.id() == vVar[i].id() ){ cols.push_back( i ); break; }
    size_t const nd = cols.size();
    for( unsigned irow=0; irow<nRes; ++irow ) nnz[irow] = nd;
    if( colnz )
      for( unsigned irow=0; irow<nRes; ++irow ){
        if( !colnz[irow] ) colnz[irow] = new size_t[ nd? nd: 1 ];
        for( size_t kk=0; kk<nd; ++kk ) colnz[irow][kk] = cols[kk];
      }
    if( !vDer || !nd ) return;
    std::vector<std::pair<FFVar,size_t>> dirs;
    for( unsigned i : cols ) dirs.push_back( _op_input_direction( i ) );
    std::string err;
    OCFESLV* P = _pOCFESLV->fdiff( dirs, err );
    if( !P ) throw std::runtime_error( "FFOCFESLV::deriv ** SYMDIFF: " + err + "\n" );
    if( !P->setup() ){ delete P; throw std::runtime_error( "FFOCFESLV::deriv ** SYMDIFF: the differentiated model failed to set up\n" ); }
    for( auto const& u : _L.ctrlRec ) P->register_control( u );
    Layout LP = _L;  LP.ctrlRec.clear();
    for( auto const& [u,spec] : P->controls() ) LP.ctrlRec.push_back( u );
    if( P->n_control_dof() != _L.nCtrl ){ delete P; throw std::runtime_error( "FFOCFESLV::deriv ** SYMDIFF: control layout mismatch\n" ); }
    FFOCFESLV ResDer;
    FFVar** ppRes = ResDer._set( nVar, vVar, LP, P, TRANSFER, _name + "_fdiff" );   // the product owns P
    for( unsigned irow=0; irow<nRes; ++irow ){
      if( !vDer[irow] ) vDer[irow] = new FFVar[ nd ];
      for( size_t kk=0; kk<nd; ++kk ) vDer[irow][kk] = *ppRes[ kk*nRes + irow ];
    }
    return;
  }

  // Dense in the map-1 DOFs: every output depends on every control DOF; map 2 and constants have no column.
  for( unsigned irow=0; irow<nRes; ++irow )
    nnz[irow] = _nCtrl;

  if( colnz ){
    for( unsigned irow=0; irow<nRes; ++irow ){
      if( !colnz[irow] ) colnz[irow] = new size_t[_nCtrl];
      for( unsigned j=0; j<_nCtrl; ++j ) colnz[irow][j] = j;
    }
  }

  if( !vDer ) return;

  for( unsigned irow=0; irow<nRes; ++irow )
    if( !vDer[irow] ) vDer[irow] = new FFVar[_nCtrl];

  // Preserve the ownership policy of the primal operation (see FFGradOCFERES).
  FFGradOCFESLV ResDer;
  int const gradPolicy = _ownOCFESLV ? FFGradOCFESLV::COPY
                                   : FFGradOCFESLV::SHALLOW;
  FFVar** vResDer = ResDer._set( nVar, vVar, _L, _pOCFESLV, gradPolicy, _name );

  size_t iel = 0;
  for( unsigned irow=0; irow<nRes; ++irow )
    for( unsigned jcol=0; jcol<_nCtrl; ++jcol, ++iel )
      vDer[irow][jcol] = *vResDer[iel];
}

//======================= FFGradOCFESLV inline implementations =======================

inline void
FFGradOCFESLV::eval
( unsigned const nRes, double* vRes, unsigned const nVar, double const* vVar,
  unsigned const* mVar )
const
{
#ifdef CRONOS__FFOCFE_TRACE
  std::cout << "FFGradOCFESLV::eval: double (dense reduced Jacobian, mode=" << options.GRADIENT << ")\n";
#endif
#ifdef CRONOS__FFOCFE_CHECK
  assert( _pOCFESLV && nRes == _nFctColl*_nCtrl && nVar == _L.nIn() );
#endif
  (void)nRes; (void)mVar;

  // Dense reduced Jacobian dF/dp (row-major), forward or adjoint per options.GRADIENT.
  std::vector<double> Jflat;
  if( !_reduced_jacobian( _pOCFESLV, _L, vVar, _nFctColl,
                                     options.GRADIENT, options.NP2NF, Jflat ) )
    throw std::runtime_error( "FFGradOCFESLV::eval double ** reduced-Jacobian march failed\n" );
  for( unsigned e=0; e<_nFctColl*_nCtrl; ++e ) vRes[e] = Jflat[e];
}

inline void
FFGradOCFESLV::eval
( unsigned const nRes, FFVar* vRes, unsigned const nVar, FFVar const* vVar,
  unsigned const* mVar )
const
{
#ifdef CRONOS__FFOCFE_CHECK
  assert( _pOCFESLV && nRes == _nFctColl*_nCtrl && nVar == _L.nIn() );
#endif
  (void)mVar;

  FFVar** ppRes = nullptr;
  if( _into_other_dag( nVar, vVar ) ){
    auto op = *this;                  // shares the solver ...
    op._ownOCFESLV = true;            // ... so that the INSERTED copy deep-copies it (FFBaseOCFE copy constructor)
    ppRes = insert_external_operation( op, nRes, nVar, vVar );
    op._ownOCFESLV = false;           // the temporary must not delete the shared solver
  }
  else
    ppRes = insert_external_operation( *this, nRes, nVar, vVar );
  for( unsigned j=0; j<nRes; ++j )
    vRes[j] = *(ppRes[j]);
}

inline void
FFGradOCFESLV::eval
( unsigned const nRes, FFDep* vRes, unsigned const nVar, FFDep const* vVar,
  unsigned const* mVar )
const
{
#ifdef CRONOS__FFOCFE_CHECK
  assert( _pOCFESLV && nRes == _nFctColl*_nCtrl && nVar == _L.nIn() );
#endif
  (void)mVar;

  FFDep dep;
  for( unsigned i=0; i<nVar; ++i ) dep += vVar[i];
  dep.update( FFDep::TYPE::N );
  for( unsigned j=0; j<nRes; ++j ) vRes[j] = dep;
}

inline void
FFGradOCFESLV::eval
( unsigned const nRes, SLiftVar* vRes, unsigned const nVar, SLiftVar const* vVar,
  unsigned const* mVar )
const
{
#ifdef CRONOS__FFOCFE_CHECK
  assert( _pOCFESLV && nRes == _nFctColl*_nCtrl && nVar == _L.nIn() );
#endif
  (void)mVar;

  vVar->env()->lift( nRes, vRes, nVar, vVar );
}

inline void
FFGradOCFESLV::eval
( unsigned const nRes, FFExpr* vRes, unsigned const nVar, FFExpr const* vVar,
  unsigned const* mVar )
const
{
#ifdef CRONOS__FFOCFE_CHECK
  assert( _pOCFESLV && nRes == _nFctColl*_nCtrl && nVar == _L.nIn() );
#endif
  (void)mVar;

  switch( FFExpr::options.LANG ){
   case FFExpr::Options::DAG:
    for( unsigned j=0; j<nRes; ++j ){
      std::ostringstream os; os << name() << "[" << j << "]";
      vRes[j] = FFExpr::compose( os.str(), nVar, vVar );
    }
    break;
   case FFExpr::Options::GAMS:
   default:
    throw typename FFExpr::Exceptions( FFExpr::Exceptions::UNDEF );
  }
}

} // end namescape mc

#endif
