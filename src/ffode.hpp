// Copyright (C) Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the EPL-2.0 with GPL-2.0-or-later as a Secondary License; see the LICENSE file.
// SPDX-License-Identifier: EPL-2.0 OR GPL-2.0-or-later

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
class FFBaseODESLV
: public FFOp
{

protected:

  //! @brief Pointer to ODE solver
  ODESLVS_CVODES*     _pODESLV;
  //! @brief Whether this class owns _pODESLV
  bool                _ownODESLV;
  //! @brief Number of parameters
  size_t              _nPar;
  //! @brief Number of constants
  size_t              _nCst;
  //! @brief Name of ODE model
  std::string         _name;
  //! @brief The solver's sensitivity directions when the op was built: op input i < nDir is parameter
  //! _ndxRec[i].  Checked at every evaluation, since a SHALLOW solver's registry could change under the op.
  std::vector<size_t> _ndxRec;
  //! @brief The other parameters, ascending: op input nDir+j is parameter _ndxRest[j] (map 2 of the call).
  std::vector<size_t> _ndxRest;

public:

  //! @brief Ordering of the operation in the DAG.  lt_FFOp has already compared the type, info and operands; then:
  //! (1) the solver (data), (2) whether the op owns a copy of it (COPY) or refers to it (SHALLOW), (3) the
  //! differentiated set _ndxRec (map 1, in gradient-column order).  The OUTPUT INDEX is NOT compared, so calls that
  //! differ only by idep share one operation; a one-map and a two-map embedding of the same solver over the same
  //! variables are two operations (their derivatives differ).
  //! @brief Is an FFVar evaluation re-inserting this operation into ANOTHER DAG (FFGraph::insert -- e.g. the
  //! per-thread copies of FFGraph::veval)?  Then the copy must own a deep copy of the solver, or threads would share
  //! one ODESLVS_CVODES and integrate concurrently (2026-10-03, as FFBaseOCFE::_into_other_dag).
  bool _into_other_dag
    ( unsigned const nVar, FFVar const* vVar )
    const
    {
      if( _ownODESLV || !nVar || varin.empty() ) return false;
      auto const* tgt = vVar[0].dag();
      for( unsigned i = 1; !tgt && i < nVar; ++i ) tgt = vVar[i].dag();
      decltype( tgt ) src = nullptr;
      for( auto const* v : varin ) if( v && v->dag() ){ src = v->dag(); break; }
      return tgt && src && tgt != src;
    }

  bool lt
    ( FFOp const* op )
    const
    {
      if( data != op->data ) return data < op->data;
      auto const* o = dynamic_cast<FFBaseODESLV const*>( op );
      if( !o ) return false;
      if( _ownODESLV != o->_ownODESLV ) return _ownODESLV < o->_ownODESLV;
      return _ndxRec < o->_ndxRec;
    }

protected:

public:

  //! @brief Default constructor
  FFBaseODESLV
    ()
    : FFOp( EXTERN ),
      _pODESLV( nullptr ),
      _ownODESLV( false ),
      _nPar( 0 ),
      _nCst( 0 ),
      _name( "" )
    {}

  //! @brief Destructor
  virtual ~FFBaseODESLV
    ()
    {
#ifdef CRONOS__FFODESLV_TRACE
      std::cout << "FFBaseODESLV::destructor\n";
#endif
      if( _ownODESLV && _pODESLV )
        delete _pODESLV;
    }

  //! @brief Copy constructor
  FFBaseODESLV
    ( FFBaseODESLV const& Op )
    : FFOp( Op ),
      _nPar( Op._nPar ),
      _nCst( Op._nCst ),
      _name( Op._name ),
      _ndxRec( Op._ndxRec ),
      _ndxRest( Op._ndxRest )
    {
#ifdef CRONOS__FFODESLV_TRACE
      std::cout << "FFBaseODESLV::copy constructor\n";
#endif
      if( !Op._pODESLV )
        throw std::runtime_error( "FFBaseODESLV::copy constructor ** Null pointer to ODE solver\n" );

      _ownODESLV = Op._ownODESLV;      
      if( _ownODESLV ){
        _pODESLV = new ODESLVS_CVODES;
        if( !_pODESLV->setup( *Op._pODESLV ) ){   // an independent copy: same description, full setup
          delete _pODESLV;  _pODESLV = nullptr;
          throw std::runtime_error( "FFBaseODESLV::copy constructor ** Setup of the ODE solver copy failed\n" );
        }
#ifdef CRONOS__FFODESLV_TRACE
        std::cerr << "ODE address copied: " << _pODESLV << std::endl;
#endif
      }
      else
        _pODESLV   = Op._pODESLV;
    }

  //! @brief d f_k / d p_i from the solver's last sensitivity solve.  Since the control registry (2026-09-24) the
  //! solver's gradient runs over its sensitivity DIRECTIONS, in control order -- row r is parameter
  //! sensitivity_index()[r] -- not over its parameters.  A COPY clears its registry so every parameter is a
  //! direction; a SHALLOW op on a solver with registered controls has no derivative for the others.
  double _dfdp
    ( size_t const i, size_t const k )
    const
    {
      _check_directions();
      return _pODESLV->val_function_gradient()[i][k];     // row i IS op input i: the registry is the op's inputs
    }

  //! @brief Record the op's contract on the solver it is built on: its directions and the rest.
  void _record
    ( ODESLVS_CVODES const* pODESLV )
    {
      _ndxRec = pODESLV->sensitivity_index();
      _ndxRest.clear();
      std::vector<bool> isdir( pODESLV->np(), false );
      for( size_t i : _ndxRec ) isdir[i] = true;
      for( size_t p=0; p<pODESLV->np(); ++p ) if( !isdir[p] ) _ndxRest.push_back( p );
    }

  //! @brief The number of op inputs the numerical gradient is taken w.r.t. (map 1).
  size_t _nDir
    () const
    { return _ndxRec.size(); }

  //! @brief The parameter index of op input @p i (i < _nPar).
  size_t _param_of_input
    ( size_t const i ) const
    { return i < _nDir()? _ndxRec[i]: _ndxRest[i-_nDir()]; }

  //! @brief Refuse if the solver's sensitivity directions are no longer those the op was built with.
  void _check_directions
    ()
    const
    {
      if( _pODESLV->sensitivity_index() != _ndxRec )
        throw std::runtime_error( "FFODESLV ** the ODE solver's registered controls changed after this operation was built"
          " (SHALLOW policy): its inputs no longer match the solver's sensitivity directions."
          " Rebuild the operation, or embed with COPY.\n" );
    }

  //! @brief The FULL parameter vector for the solver, from the op inputs [ directions | the other parameters ].
  void _params
    ( double const* vVar, std::vector<double>& P )
    const
    {
      P.assign( _nPar, 0. );
      for( size_t i=0; i<_nPar; ++i ) P[ _param_of_input( i ) ] = vVar[i];
    }

  //! @brief Resolve the TWO MAPS -- @p vDiff (the numerical gradient's inputs; sets the registry) and @p vRest
  //! (every other declared input, and every constant) -- into the op inputs @p vPar = [ map-1 DOFs in control
  //! order | the other parameters in parameter order ] and @p vCst (constants, declared order).  Refuses an
  //! undeclared entry, a wrong DOF count, an entry listed twice, a constant in map 1, and any declared input
  //! (other than a fixed one) or constant missing from both maps.  @p saved receives the previous registry.
  static void _bind_maps
    ( std::vector<FFModel::InputArg> const& vDiff, std::vector<FFModel::InputArg> const& vRest,
      ODESLVS_CVODES* pODESLV, std::vector<FFVar>& saved, std::vector<FFVar>& vPar, std::vector<FFVar>& vCst )
    {
      if( !pODESLV ) throw std::runtime_error( "FFODESLV ** Null pointer to ODE solver\n" );
      auto const& vC = pODESLV->var_constant();
      auto cst_index = [&]( FFVar const& v )->int {
        for( size_t c=0; c<vC.size(); ++c ) if( vC[c].id().second == v.id().second ) return (int)c;
        return -1; };
      auto is_input = [&]( FFVar const& v ){
        for( auto const& [w,dom] : pODESLV->var_declared_input() ) if( w.id().second == v.id().second ) return true;
        return false; };
      std::vector<FFVar> seen;
      auto once = [&]( FFVar const& v ){
        for( auto const& u : seen ) if( u.id().second == v.id().second )
          throw std::runtime_error( "FFODESLV ** " + v.name() + " appears more than once in the two maps\n" );
        seen.push_back( v ); };
      std::string err;
      std::vector<std::vector<FFVar>> dDiff( vDiff.size() ), dRest( vRest.size() );
      std::vector<int> cRest( vRest.size(), -1 );
      for( size_t i=0; i<vDiff.size(); ++i ){
        once( vDiff[i].input );
        if( cst_index( vDiff[i].input ) >= 0 )
          throw std::runtime_error( "FFODESLV ** constant " + vDiff[i].input.name() + " cannot be in map 1: the solver's"
            " sensitivities run over inputs only (name it in SYMDIFF for a symbolic derivative)\n" );
        if( !pODESLV->resolve_input_arg( vDiff[i], dDiff[i], err ) ) throw std::runtime_error( "FFODESLV ** " + err + "\n" );
      }
      for( size_t i=0; i<vRest.size(); ++i ){
        once( vRest[i].input );
        cRest[i] = cst_index( vRest[i].input );
        if( cRest[i] >= 0 ){                                   // a constant: exactly one variable
          if( vRest[i].gen ) dRest[i].assign( 1, vRest[i].gen( FFModel::DofIndex() ) );
          else               dRest[i] = vRest[i].dofs;
          if( dRest[i].size() != 1 )
            throw std::runtime_error( "FFODESLV ** constant " + vRest[i].input.name() + " takes exactly one variable\n" );
        }
        else if( !pODESLV->resolve_input_arg( vRest[i], dRest[i], err ) ) throw std::runtime_error( "FFODESLV ** " + err + "\n" );
      }
      // completeness: every declared input (not fixed) and every constant, in exactly one map
      std::string missing;
      for( auto const& [w,dom] : pODESLV->var_declared_input() ){
        if( pODESLV->fixed_input( w ) ) continue;
        bool in = false; for( auto const& u : seen ) if( u.id().second == w.id().second ){ in = true; break; }
        if( !in ) missing += " " + w.name();
      }
      for( auto const& c : vC ){
        bool in = false; for( auto const& u : seen ) if( u.id().second == c.id().second ){ in = true; break; }
        if( !in ) missing += " " + c.name();
      }
      for( auto const& u : seen )
        if( !is_input( u ) && cst_index( u ) < 0 )
          throw std::runtime_error( "FFODESLV ** " + u.name() + " is neither a declared input nor a constant of the model\n" );
      if( !missing.empty() )
        throw std::runtime_error( "FFODESLV ** every declared input and constant must be mapped (map 1 or map 2) so its"
          " value can propagate through the DAG; missing:" + missing + "\n" );
      // the registry is map 1
      saved.clear();
      for( auto const& [u,spec] : pODESLV->controls() ) saved.push_back( u );
      pODESLV->clear_controls();
      for( auto const& a : vDiff ) pODESLV->register_control( a.input );
      size_t const nsen = pODESLV->sensitivity_index().size();
      if( nsen != pODESLV->n_control_dof() ){
        _restore_controls( pODESLV, saved );
        throw std::runtime_error( "FFODESLV ** map 1 could not be resolved to parameters of the extracted model: "
          + pODESLV->extract_error() + "\n" );
      }
      // op inputs: [ map-1 DOFs, control order | the other parameters, parameter order ]
      vPar.assign( pODESLV->np(), FFVar() );
      for( size_t i=0; i<vDiff.size(); ++i ){
        auto const b = pODESLV->control_block( vDiff[i].input );
        for( size_t k=0; k<b.ndof; ++k ) vPar[ b.offset+k ] = dDiff[i][k];
      }
      std::vector<bool> isdir( pODESLV->np(), false );
      for( size_t p : pODESLV->sensitivity_index() ) isdir[p] = true;
      std::vector<size_t> posRest( pODESLV->np(), 0 );
      for( size_t p=0, j=0; p<pODESLV->np(); ++p ) if( !isdir[p] ) posRest[p] = j++;
      for( size_t i=0; i<vRest.size(); ++i ){
        if( cRest[i] >= 0 ) continue;
        auto const pidx = pODESLV->parameter_index( vRest[i].input );
        for( size_t k=0; k<pidx.size() && k<dRest[i].size(); ++k ) vPar[ nsen + posRest[ pidx[k] ] ] = dRest[i][k];
      }
      vCst.assign( vC.size(), FFVar() );
      for( size_t i=0; i<vRest.size(); ++i ) if( cRest[i] >= 0 ) vCst[ cRest[i] ] = dRest[i][0];
    }

  static void _restore_controls
    ( ODESLVS_CVODES* pODESLV, std::vector<FFVar> const& saved )
    { pODESLV->clear_controls();  for( auto const& u : saved ) pODESLV->register_control( u ); }

  //! @brief Sensitivity mode for FFODESLV / FFGradODESLV (as FFBaseOCFE::GRADIENT_TYPE): FORWARD uses forward
  //! sensitivities (one direction per parameter); ADJOINT uses adjoint sensitivities (one backward sweep per output
  //! function); AUTO picks adjoint when np > NP2NF * nf -- the rule FFODESLV always applied.
  enum GRADIENT_TYPE{ FORWARD=0, ADJOINT=1, AUTO=2 };
  //! @brief FFODESLV options - static so that can still be accessed/modified after setup
  static struct Options
  {
    //! @brief Constructor
    Options():
      SYMDIFF (),
      NP2NF   (3),
      GRADIENT( AUTO )
      {}
    //! @brief Assignment operator
    Options& operator= ( Options const& options ){
        SYMDIFF   = options.SYMDIFF;
        NP2NF     = options.NP2NF;
        GRADIENT  = options.GRADIENT;
        return *this;
      }
    //! @brief Variable selection for symbolic differentiation - needs to be participating parameters or constants. Applies numerical differentiation w.r.t. parameters if empty
    std::vector<FFVar>        SYMDIFF;
    //! @brief parameter-to-function-size ratio above which the AUTO gradient uses adjoint instead of forward sensitivity
    double                    NP2NF;
    //! @brief sensitivity mode (FORWARD / ADJOINT / AUTO) for FFODESLV / FFGradODESLV
    int                       GRADIENT;
  } options;
  //! @brief Whether a derivative with @p np parameters and @p nf functions uses ADJOINT sensitivities (options.GRADIENT;
  //! AUTO: np > NP2NF * nf).  The single rule of both FFODESLV and FFGradODESLV.
  static bool _use_adjoint
    ( size_t const np, size_t const nf )
    { return options.GRADIENT == ADJOINT || ( options.GRADIENT == AUTO && (double)np > options.NP2NF * (double)nf ); }
};

inline FFBaseODESLV::Options FFBaseODESLV::options;

//! @brief C++ class defining IVP in parametric ODEs as external DAG operations in MC++.
////////////////////////////////////////////////////////////////////////
//! mc::FFODESLV is a C++ class for defining an IVP in parametric ODEs as
//! external DAG operations in MC++.
////////////////////////////////////////////////////////////////////////
class FFODESLV
: public FFBaseODESLV
{

public:

  //! @brief Enumeration type for ODE copy policy
  enum POLICY_TYPE{
    SHALLOW=0,  //!< Shallow copy of ODE system in FFGraph (without ownership)
    COPY=1,     //!< Deep copy of ODE system in FFGraph (with ownership)
    TRANSFER=-1 //!< Shallow copy of ODE system in FFGraph (with ownership transfer)
  };

  // Default constructor
  FFODESLV
    ()
    : FFBaseODESLV()
    {}

  // Destructor
  virtual ~FFODESLV
    ()
    {}

  //! @brief Bind the two maps, insert the op (via _set), and -- for COPY -- restore the caller's registry.
  FFVar** _embed
    ( std::vector<FFModel::InputArg> const& vDiff, std::vector<FFModel::InputArg> const& vRest,
      ODESLVS_CVODES* pODESLV, int policy, std::string const& name );

  //! @brief TWO-MAP form.  @p vDiff (map 1) lists the inputs the NUMERICAL gradient is taken w.r.t.: it RESETS the
  //! solver's control registry to exactly these.  @p vRest (map 2) lists every other declared input and every
  //! constant: their values pass through the DAG, with ZERO numerical-derivative columns by contract.  Each entry
  //! is a flat vector in control_dofs() order or a generator called per DofIndex; a constant takes one variable.
  //! Every declared input (except a fixed one) and every constant must appear in exactly one map.  SYMDIFF may
  //! name any op input or constant for a SYMBOLIC derivative.  COPY leaves the caller's registry as it was;
  //! SHALLOW and TRANSFER keep the new one.  Returns the nf outputs.
  std::vector<FFVar> operator()
    ( std::vector<FFModel::InputArg> const& vDiff, std::vector<FFModel::InputArg> const& vRest,
      ODESLVS_CVODES* pODESLV, int policy=COPY, std::string const& name="" )
    {
      FFVar** ppRes = _embed( vDiff, vRest, pODESLV, policy, name );
      std::vector<FFVar> vRes( pODESLV->nf() );
      for( size_t i=0; i<vRes.size(); ++i ) vRes[i] = *ppRes[i];
      return vRes;
    }

  //! @brief TWO-MAP form, output @p idep only.
  FFVar& operator()
    ( unsigned const idep, std::vector<FFModel::InputArg> const& vDiff, std::vector<FFModel::InputArg> const& vRest,
      ODESLVS_CVODES* pODESLV, int policy=COPY, std::string const& name="" )
    {
      FFVar** ppRes = _embed( vDiff, vRest, pODESLV, policy, name );
      if( idep >= pODESLV->nf() ) throw std::runtime_error( "FFODESLV ** output index out of range\n" );
      return *ppRes[idep];
    }

  //! @brief ONE-MAP form: every declared input is differentiated and the model has no constant -- the two-map
  //! form with an empty map 2, refused (naming what is missing) otherwise.
  std::vector<FFVar> operator()
    ( std::vector<FFModel::InputArg> const& vDiff, ODESLVS_CVODES* pODESLV, int policy=COPY,
      std::string const& name="" )
    { return (*this)( vDiff, std::vector<FFModel::InputArg>(), pODESLV, policy, name ); }

  //! @brief ONE-MAP form, output @p idep only.
  FFVar& operator()
    ( unsigned const idep, std::vector<FFModel::InputArg> const& vDiff, ODESLVS_CVODES* pODESLV,
      int policy=COPY, std::string const& name="" )
    { return (*this)( idep, vDiff, std::vector<FFModel::InputArg>(), pODESLV, policy, name ); }

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
      else if( mc::same_type( idU, typeid( SLiftVar ) ) )
        return eval( nRes, static_cast<SLiftVar*>(vRes), nVar, static_cast<SLiftVar const*>(vVar), mVar );
      else if( mc::same_type( idU, typeid( FFExpr ) ) )
        return eval( nRes, static_cast<FFExpr*>(vRes), nVar, static_cast<FFExpr const*>(vVar), mVar );

      throw std::runtime_error( "FFODESLV::feval ** No evaluation method for type"+std::string(idU.name())+"\n" );
    }

  void eval
    ( unsigned const nRes, double* vRes, unsigned const nVar, double const* vVar, unsigned const* mVar )
    const;

  void eval
    ( unsigned const nRes, FADType<double>* vRes, unsigned const nVar, FADType<double> const* vVar,
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

  // Ordering in the DAG: FFBaseODESLV::lt (solver, policy, differentiated set -- not the output index)

  // Properties
  std::string name
    ()
    const
    { std::ostringstream oss;
      if( !_name.empty() ) oss << _name;
      else                  oss << _pODESLV;
      return "ODE[" + oss.str() + "]"; }

  //! @brief Return whether or not operation is commutative
  bool commutative
    ()
    const
    { return false; }

protected:

  FFVar** _set
    ( unsigned const nPar, FFVar const* pPar, unsigned const nCst, FFVar const* pCst,
      ODESLVS_CVODES* pODESLV, int policy, std::string const& name,
      bool const named=false )
//    const
    {
#ifdef CRONOS__FFODESLV_CHECK
      assert( pODESLV && ( !pCst || nCst == pODESLV->var_constant().size() ) );
#endif
      // The op's inputs are every parameter, laid out [ directions (the registry) | the other parameters ];
      // the positional form (named=false) means parameter order, so it needs NO control registered.
      if( !pODESLV ) throw std::runtime_error( "FFODESLV ** Null pointer to ODE solver\n" );
      if( nPar != pODESLV->np() )
        throw std::runtime_error( "FFODESLV ** the operation must receive all " + std::to_string( pODESLV->np() )
          + " parameters of the solver, " + std::to_string( nPar ) + " given\n" );
      if( !named && pODESLV->sensitivity_index().size() != pODESLV->np() )
        throw std::runtime_error( "FFODESLV ** the solver has registered controls, so its inputs are not in parameter order:"
          " use the two-map form, which sets the registered controls from map 1\n" );
      if( _ownODESLV && _pODESLV )
        delete _pODESLV;
      _ownODESLV = ( policy>0? true: false ); //copy;
      _pODESLV = pODESLV;
      _nPar = nPar;
      _nCst = nCst;
      _name = name;
      _record( pODESLV );

      data = pODESLV;
      owndata = false;
      size_t const nFun = pODESLV->nf();
      FFVar** ppRes = ( pCst? insert_external_operation( *this, nFun, nPar, pPar, nCst, pCst ):
                              insert_external_operation( *this, nFun, nPar, pPar ) );

      _ownODESLV = false;
      FFOp* pOp = (*ppRes)->opdef().first;
      if( policy > 0 )
        _pODESLV = dynamic_cast<FFODESLV*>(pOp)->_pODESLV; // set pointer to DAG copy
      else if( policy < 0 )
        dynamic_cast<FFODESLV*>(pOp)->_ownODESLV = true;   // transfer ownership
#ifdef CRONOS__FFODESLV_TRACE
      std::cerr << "ODE operation address: " << this << std::endl;
      std::cerr << "ODE address in DAG: " << _pODESLV << std::endl;
#endif
      return ppRes;
    }
};

//! @brief C++ class defining gradient of IVP in parametric ODEs as external DAG operations in MC++.
////////////////////////////////////////////////////////////////////////
//! mc::FFGradODESLV is a C++ class for defining the gradient of an IVP
//! in parametric ODEs as external DAG operations in MC++.
////////////////////////////////////////////////////////////////////////
class FFGradODESLV
: public FFBaseODESLV
{
  friend class FFODESLV;

public:

  //! @brief Enumeration type for ODE copy policy
  enum POLICY_TYPE{
    SHALLOW=0,  //!< Shallow copy of ODE system in FFGraph (without ownership)
    COPY=1      //!< Deep copy of ODE system in FFGraph (with ownership)
  };

  // Constructors
  FFGradODESLV
    ()
    : FFBaseODESLV()
    {}

  // Destructor
  virtual ~FFGradODESLV
    ()
    {}

  // No public call operator: FFGradODESLV is built by FFODESLV::deriv / eval( FADType<FFVar> ) through _set(),
  // with the same op-input layout as the FFODESLV it differentiates.

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

      throw std::runtime_error( "FFGradODESLV::feval ** No evaluation method for type"+std::string(idU.name())+"\n" );
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

  // Ordering in the DAG: FFBaseODESLV::lt (solver, policy, differentiated set -- not the output index)

  // Properties
  std::string name
    ()
    const
    { std::ostringstream oss;
      if( !_name.empty() ) oss << _name;
      else                  oss << _pODESLV;
      return "GRAD_ODE[" + oss.str() + "]"; }

  //! @brief Return whether or not operation is commutative
  bool commutative
    ()
    const
    { return false; }

protected:

  FFVar** _set
    ( unsigned const nPar, FFVar const* pPar, unsigned const nCst, FFVar const* pCst,
      ODESLVS_CVODES* pODESLV, int policy, std::string const& name="",
      bool const named=false )
    {
#ifdef CRONOS__FFODESLV_CHECK
  assert( pODESLV && ( !pCst || nCst == pODESLV->var_constant().size() ) );
#endif
      // The op's inputs are every parameter, laid out [ directions (the registry) | the other parameters ];
      // the positional form (named=false) means parameter order, so it needs NO control registered.
      if( !pODESLV ) throw std::runtime_error( "FFODESLV ** Null pointer to ODE solver\n" );
      if( nPar != pODESLV->np() )
        throw std::runtime_error( "FFODESLV ** the operation must receive all " + std::to_string( pODESLV->np() )
          + " parameters of the solver, " + std::to_string( nPar ) + " given\n" );
      if( !named && pODESLV->sensitivity_index().size() != pODESLV->np() )
        throw std::runtime_error( "FFODESLV ** the solver has registered controls, so its inputs are not in parameter order:"
          " use the two-map form, which sets the registered controls from map 1\n" );
      if( _ownODESLV && _pODESLV )
        delete _pODESLV;
      _ownODESLV = ( policy!=0? true: false );
      _pODESLV = pODESLV;
      _nPar = nPar;
      _nCst = nCst;
      _name = name;
      _record( pODESLV );

      data = pODESLV;
      owndata = false;
      size_t const nFun = pODESLV->nf(), nDir = _nDir();   // the gradient w.r.t. map 1 only
      FFVar** ppRes = ( pCst? insert_external_operation( *this, nFun*nDir, nPar, pPar, nCst, pCst ):
                              insert_external_operation( *this, nFun*nDir, nPar, pPar ) );

      _ownODESLV = false;
      FFOp* pOp = (*ppRes)->opdef().first;
      if( policy > 0 )
        _pODESLV = dynamic_cast<FFGradODESLV*>(pOp)->_pODESLV; // set pointer to DAG copy
#ifdef CRONOS__FFODESLV_TRACE
      std::cerr << "FFGradODESLV address in DAG: " << _pODESLV << std::endl;
#endif
      return ppRes;
    }
};

inline FFVar**
FFODESLV::_embed
( std::vector<FFModel::InputArg> const& vDiff, std::vector<FFModel::InputArg> const& vRest,
  ODESLVS_CVODES* pODESLV, int policy, std::string const& name )
{
  std::vector<FFVar> saved, vPar, vCst;
  _bind_maps( vDiff, vRest, pODESLV, saved, vPar, vCst );
  FFVar** ppRes = vCst.empty()? _set( vPar.size(), vPar.data(), 0, nullptr, pODESLV, policy, name, true )
                              : _set( vPar.size(), vPar.data(), vCst.size(), vCst.data(), pODESLV, policy, name, true );
  if( policy > 0 ) _restore_controls( pODESLV, saved );     // the op's copy keeps the registry
  return ppRes;
}

inline void
FFODESLV::eval
( unsigned const nRes, FFVar* vRes, unsigned const nVar, FFVar const* vVar,
  unsigned const* mVar )
const
{
#ifdef CRONOS__FFODESLV_TRACE
  std::cout << "FFODESLV::eval: FFVar\n";
  std::cerr << "ODE operation address: " << this << std::endl;
  std::cerr << "ODE address in DAG: " << _pODESLV << std::endl;
#endif
#ifdef CRONOS__FFODESLV_CHECK
  assert( _pODESLV && nRes == _pODESLV->nf() && nVar == _nPar+_nCst );
#endif

  FFVar** ppRes = nullptr;
  if( _into_other_dag( nVar, vVar ) ){
    auto op = *this;                  // shares the solver ...
    op._ownODESLV = true;             // ... so that the INSERTED copy deep-copies it (FFBaseODESLV copy constructor)
    ppRes = ( _nCst? insert_external_operation( op, nRes, _nPar, vVar, _nCst, vVar+_nPar ):
                     insert_external_operation( op, nRes, _nPar, vVar ) );
    op._ownODESLV = false;            // the temporary must not delete the shared solver
  }
  else
    ppRes = ( _nCst? insert_external_operation( *this, nRes, _nPar, vVar, _nCst, vVar+_nPar ):
                     insert_external_operation( *this, nRes, _nPar, vVar ) );
  //FFVar** ppRes = _set( _nPar, vVar, _nCst, vVar+_nPar, static_cast<ODESLVS_CVODES*>(data), _ownODESLV );
  for( unsigned j=0; j<nRes; ++j )
    vRes[j] = *(ppRes[j]);
}

inline void
FFODESLV::eval
( unsigned const nRes, FFDep* vRes, unsigned const nVar, FFDep const* vVar,
  unsigned const* mVar )
const
{
#ifdef CRONOS__FFODESLV_TRACE
  std::cout << "FFODESLV::eval: FFDep\n";
#endif
#ifdef CRONOS__FFODESLV_CHECK
  assert( _pODESLV && nRes == _pODESLV->nf() && nVar == _nPar+_nCst );
#endif

  vRes[0] = 0;
  for( unsigned i=0; i<nVar; ++i ) vRes[0] += vVar[i];
  vRes[0].update( FFDep::TYPE::N );
  for( unsigned j=1; j<nRes; ++j ) vRes[j] = vRes[0];
}

inline void
FFODESLV::eval
( unsigned const nRes, FFExpr* vRes, unsigned const nVar, FFExpr const* vVar,
  unsigned const* mVar )
const
{
#ifdef CRONOS__FFODESLV_TRACE
  std::cout << "FFODESLV::eval: FFExpr\n";
#endif
#ifdef CRONOS__FFODESLV_CHECK
  assert( _pODESLV && nVar == _nPar+_nCst && nRes == _pODESLV->nf() );
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
FFODESLV::eval
( unsigned const nRes, double* vRes, unsigned const nVar, double const* vVar,
  unsigned const* mVar )
const
{
#ifdef CRONOS__FFODESLV_TRACE
  std::cerr << "FFODESLV::eval: double w/ ODE address " << _pODESLV << std::endl;
  for( unsigned i=0; i<nVar; ++i ) std::cout << "vVar[" << i << "] = " << vVar[i] << std::endl;
#endif
#ifdef CRONOS__FFODESLV_CHECK
  assert( _pODESLV && nRes == _pODESLV->nf() && nVar == _nPar+_nCst );
#endif

  _check_directions();
  std::vector<double> P;  _params( vVar, P );
  if( _pODESLV->solve( P.data(), _nCst? vVar+_nPar: nullptr ) != ODESLVS_CVODES::NORMAL ){
    //for( unsigned i=0; i<nVar; ++i ) std::cout << "vVar[" << i << "] = " << vVar[i] << std::endl;
    throw std::runtime_error( "FFODESLV::eval double ** State integration failure\n" );
  }

  for( unsigned i=0; i<nRes; ++i )
    vRes[i] = _pODESLV->val_function()[i];
}

inline void
FFODESLV::eval
( unsigned const nRes, FADType<double>* vRes, unsigned const nVar, FADType<double> const* vVar,
  unsigned const* mVar )
const
{
#ifdef CRONOS__FFODESLV_TRACE
  std::cout << "FFODESLV::eval: FADType<double>\n";
#endif
#ifdef CRONOS__FFODESLV_CHECK
  assert( _pODESLV && nRes == _pODESLV->nf() && nVar == _nPar+_nCst );
#endif

  std::vector<double> vVarVal( nVar );
  for( unsigned i=0; i<nVar; ++i ) vVarVal[i] = vVar[i].val();
  _check_directions();
  std::vector<double> P;  _params( vVarVal.data(), P );
  if( !_use_adjoint( _nPar, nRes ) ){
    if( _pODESLV->solve_fsens( P.data(), _nCst? vVarVal.data()+_nPar: nullptr ) != ODESLVS_CVODES::NORMAL )
      throw std::runtime_error( "FFODESLV::eval FADType<double> ** Forward sensitivity integration failure\n" );
  }
  else{
    if( _pODESLV->solve_asens( P.data(), _nCst? vVarVal.data()+_nPar: nullptr ) != ODESLVS_CVODES::NORMAL )
      throw std::runtime_error( "FFODESLV::eval FADType<double> ** Adjoint sensitivity integration failure\n" );
  }

  for( unsigned k=0; k<nRes; ++k ){
    vRes[k] = _pODESLV->val_function()[k];
    for( unsigned i=0; i<_nPar; ++i )
      vRes[k].setDepend( vVar[i] );
    for( unsigned j=0; j<vRes[k].size(); ++j ){
      vRes[k][j] = 0.;
      for( unsigned i=0; i<_nDir(); ++i ){              // map 2 and constants: zero by contract
        if( vVar[i][j] == 0. ) continue;
        vRes[k][j] += _dfdp( i, k ) * vVar[i][j];
      }
    }
  }
}

inline void
FFODESLV::eval
( unsigned const nRes, FADType<FFVar>* vRes, unsigned const nVar, FADType<FFVar> const* vVar,
  unsigned const* mVar )
const
{
#ifdef CRONOS__FFODESLV_TRACE
  std::cout << "FFODESLV::eval: FADType<FFVar>\n";
#endif
#ifdef CRONOS__FFODESLV_CHECK
  assert( _pODESLV && nRes == _pODESLV->nf() && nVar >= _pODESLV->np() );
#endif

  std::vector<FFVar> vVarVal( nVar );
  for( unsigned i=0; i<nVar; ++i )
    vVarVal[i] = vVar[i].val();
  FFVar const*const* vResVal = ( _nCst? insert_external_operation( *this, nRes, _nPar, vVarVal.data(), _nCst, vVarVal.data()+_nPar ):
                                        insert_external_operation( *this, nRes, _nPar, vVarVal.data() ) );

  if( options.SYMDIFF.empty() ){
    FFGradODESLV ResDer;
    // No DAG copy of ODE - reuse FFODESLV DAG copy
    // Caveat is that passing a pointer to the orginal ODESLV object will create a separate object
    // FFVar const*const* vResDer = ResDer( nVar, vVarVal.data(), _pODESLV, false );
    // DAG copy of ODE - no resuse of FFODESLV DAG copy
    // Caveat is external data pointer may change
    FFVar const*const* vResDer = ResDer._set( _nPar, vVarVal.data(), _nCst, vVarVal.data()+_nPar, _pODESLV, _ownODESLV, _name, true );   // the op's OWN solver: a COPY's registry is not the caller's
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

  else{
    // Match DAG variables to IVP parameters
    std::vector<FFVar> vPar;
    vPar.reserve( options.SYMDIFF.size() );
    std::vector<unsigned> ndxPar;
    ndxPar.reserve( options.SYMDIFF.size() );
    for( auto const& dVar : options.SYMDIFF ){
      for( unsigned i=0; i<nVar; ++i ){
        if( dVar.id() != vVarVal[i].id() ) continue;
#ifdef CRONOS__FFODESLV_CHECK
        std::cout << "Sensitivity parameter #" << i << ": " << dVar << std::endl;
#endif
        vPar.push_back( i<_nPar? _pODESLV->var_parameter()[ _param_of_input( i ) ]: _pODESLV->var_constant()[i-_nPar] );
        ndxPar.push_back( i );
        break;
      }
    }

    FFODESLV ResDer;
    auto pODESLVSEN = _pODESLV->fdiff( vPar );
    if( !pODESLVSEN->setup() ){
      std::string const why = pODESLVSEN->extract_error();  delete pODESLVSEN;
      throw std::runtime_error( "FFODESLV ** SYMDIFF: the forward-differentiated model failed to set up: " + why + "\n" );
    }
    // the product has exactly the original parameters; give it the original's registry so its op takes the
    // same inputs, in the same order, with the same baseline
    for( auto const& [u,spec] : _pODESLV->controls() ) pODESLVSEN->register_control( u );
    FFVar const*const* vResDer = ResDer._set( _nPar, vVarVal.data(), _nCst, vVarVal.data()+_nPar, pODESLVSEN, -1, _name, true ); // transfer ownership of ODE data to DAG
    for( unsigned k=0; k<nRes; ++k ){
      vRes[k] = *vResVal[k];
      for( unsigned ii=0; ii<ndxPar.size(); ++ii )
        vRes[k].setDepend( vVar[ndxPar[ii]] );
      for( unsigned j=0; j<vRes[k].size(); ++j ){
        vRes[k][j] = 0.;
        for( unsigned ii=0; ii<ndxPar.size(); ++ii ){
          if( vVar[ndxPar[ii]][j].cst() && vVar[ndxPar[ii]][j].num().val() == 0. ) continue;
          vRes[k][j] += *vResDer[k+nRes*ii] * vVar[ndxPar[ii]][j];
        }
      }
    }
  }
}

inline void
FFODESLV::deriv
( unsigned const nRes, FFVar const* vRes, unsigned const nVar, FFVar const* vVar, FFVar** vDer )
const
{
#ifdef CRONOS__FFODESLV_TRACE
  std::cout << "FFODESLV::deriv\n";
#endif
#ifdef CRONOS__FFODESLV_CHECK
  assert( _pODESLV && nRes == _pODESLV->nf() && nVar == _nPar+_nCst );
#endif

  if( options.SYMDIFF.empty() ){
    FFGradODESLV ResDer;
    // No DAG copy of ODE - reuse FFODESLV DAG copy
    // Caveat is that passing a pointer to the orginal ODESLV object will create a separate object
    // FFVar const*const* vResDer = ResDer( nVar, vVar, _pODESLV, false ); // no copy of ODE data
    // DAG copy of ODE - no resuse of FFODESLV DAG copy
    // Caveat is external data pointer may change
    FFVar const*const* vResDer = ResDer._set( _nPar, vVar, _nCst, vVar+_nPar, _pODESLV, _ownODESLV, _name, true );   // the op's OWN solver: a COPY's registry is not the caller's
    for( unsigned k=0; k<nRes; ++k )
      for( unsigned i=0; i<nVar; ++i )
        vDer[k][i] = i<_nDir()? *vResDer[k+nRes*i]: 0;    // map 2 and constants: zero by contract
  }

  else{
    // Match DAG variables to IVP parameters
    std::vector<FFVar> vPar; 
    vPar.reserve( options.SYMDIFF.size() );
    std::vector<unsigned> ndxPar; 
    ndxPar.reserve( options.SYMDIFF.size() );
    for( auto const& dVar : options.SYMDIFF ){
      for( unsigned i=0; i<nVar; ++i ){
        if( dVar.id() != vVar[i].id() ) continue;
#ifdef CRONOS__FFODESLV_CHECK
        std::cout << "Sensitivity parameter #" << i << ": " << dVar << std::endl;
#endif
        vPar.push_back( i<_nPar? _pODESLV->var_parameter()[ _param_of_input( i ) ]: _pODESLV->var_constant()[i-_nPar] );
        ndxPar.push_back( i );
        break;
      }
    }

    FFODESLV ResDer;
    auto pODESLVSEN = _pODESLV->fdiff( vPar );
    if( !pODESLVSEN->setup() ){
      std::string const why = pODESLVSEN->extract_error();  delete pODESLVSEN;
      throw std::runtime_error( "FFODESLV ** SYMDIFF: the forward-differentiated model failed to set up: " + why + "\n" );
    }
    // the product has exactly the original parameters; give it the original's registry so its op takes the
    // same inputs, in the same order, with the same baseline
    for( auto const& [u,spec] : _pODESLV->controls() ) pODESLVSEN->register_control( u );
    FFVar const*const* vResDer = ResDer._set( _nPar, vVar, _nCst, vVar+_nPar, pODESLVSEN, -1, _name, true ); // transfer ownership of ODE data to DAG
    for( unsigned k=0; k<nRes; ++k ){
      for( unsigned i=0; i<nVar; ++i )
        vDer[k][i] = 0;
      for( unsigned ii=0; ii<ndxPar.size(); ++ii )
        vDer[k][ndxPar[ii]] = *vResDer[k+nRes*ii];
    }
  }
}

inline void
FFODESLV::eval
( unsigned const nRes, SLiftVar* vRes, unsigned const nVar, SLiftVar const* vVar,
  unsigned const* mVar )
const
{
#ifdef CRONOS__FFODESLV_TRACE
  std::cout << "FFODESLV::eval: SLiftVar\n";
#endif
#ifdef CRONOS__FFODESLV_CHECK
  assert( _pODESLV && nRes == _pODESLV->nf() && nVar == _nPar+_nCst );
#endif

  vVar->env()->lift( nRes, vRes, nVar, vVar );
}

inline void
FFGradODESLV::eval
( unsigned const nRes, FFVar* vRes, unsigned const nVar, FFVar const* vVar,
  unsigned const* mVar )
const
{
#ifdef CRONOS__FFODESLV_TRACE
  std::cout << "FFGradODESLV::eval: FFVar\n";
#endif
#ifdef CRONOS__FFODESLV_CHECK
  assert( _pODESLV && nVar == _nPar+_nCst && nRes == _pODESLV->nf()*_nDir() );
#endif

  FFVar** ppRes = nullptr;
  if( _into_other_dag( nVar, vVar ) ){
    auto op = *this;                  // shares the solver ...
    op._ownODESLV = true;             // ... so that the INSERTED copy deep-copies it (FFBaseODESLV copy constructor)
    ppRes = ( _nCst? insert_external_operation( op, nRes, _nPar, vVar, _nCst, vVar+_nPar ):
                     insert_external_operation( op, nRes, _nPar, vVar ) );
    op._ownODESLV = false;            // the temporary must not delete the shared solver
  }
  else
    ppRes = ( _nCst? insert_external_operation( *this, nRes, _nPar, vVar, _nCst, vVar+_nPar ):
                     insert_external_operation( *this, nRes, _nPar, vVar ) );
  //FFVar** ppRes = _set( _nPar, vVar, _nCst, vVar+_nPar, static_cast<ODESLVS_CVODES*>(data), _ownODESLV );
  for( unsigned j=0; j<nRes; ++j ) vRes[j] = *(ppRes[j]);
}

inline void
FFGradODESLV::eval
( unsigned const nRes, FFDep* vRes, unsigned const nVar, FFDep const* vVar,
  unsigned const* mVar )
const
{
#ifdef CRONOS__FFODESLV_TRACE
  std::cout << "FFGradODESLV::eval: FFDep\n";
#endif
#ifdef CRONOS__FFODESLV_CHECK
  assert( _pODESLV && nVar == _nPar+_nCst && nRes == _pODESLV->nf()*_nDir() );
#endif

  vRes[0] = 0;
  for( unsigned i=0; i<nVar; ++i ) vRes[0] += vVar[i];
  vRes[0].update( FFDep::TYPE::N );
  for( unsigned j=1; j<nRes; ++j ) vRes[j] = vRes[0];
}

inline void
FFGradODESLV::eval
( unsigned const nRes, FFExpr* vRes, unsigned const nVar, FFExpr const* vVar,
  unsigned const* mVar )
const
{
#ifdef CRONOS__FFODESLV_TRACE
  std::cout << "FFGradODESLV::eval: FFExpr\n";
#endif
#ifdef CRONOS__FFODESLV_CHECK
  assert( _pODESLV && nVar == _nPar+_nCst && nRes == _pODESLV->nf()*_nDir() );
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
FFGradODESLV::eval
( unsigned const nRes, double* vRes, unsigned const nVar, double const* vVar,
  unsigned const* mVar )
const
{
#ifdef CRONOS__FFODESLV_TRACE
  std::cout << "FFGradODESLV::eval: double\n";
#endif
#ifdef CRONOS__FFODESLV_CHECK
  assert( _pODESLV && nVar == _nPar+_nCst && nRes == _pODESLV->nf()*_nDir() );
#endif

  _check_directions();
  std::vector<double> P;  _params( vVar, P );
  if( !_use_adjoint( _nPar, _pODESLV->nf() ) ){
    if( _pODESLV->solve_fsens( P.data(), _nCst? vVar+_nPar: nullptr ) != ODESLVS_CVODES::NORMAL )
      throw std::runtime_error( "FFGradODESLV::eval double ** Forward sensitivity integration failure\n" );
  }
  else{
    if( _pODESLV->solve_asens( P.data(), _nCst? vVar+_nPar: nullptr ) != ODESLVS_CVODES::NORMAL )
      throw std::runtime_error( "FFGradODESLV::eval double ** Adjoint sensitivity integration failure\n" );
  }

  for( unsigned i=0, k=0; i<_nDir(); ++i )
    for( unsigned j=0; j<_pODESLV->nf(); ++j, ++k )
      vRes[k] = _dfdp( i, j );
}

inline void
FFGradODESLV::eval
( unsigned const nRes, SLiftVar* vRes, unsigned const nVar, SLiftVar const* vVar,
  unsigned const* mVar )
const
{
#ifdef CRONOS__FFODESLV_TRACE
  std::cout << "FFGradODESLV::eval: SLiftVar\n";
#endif
#ifdef CRONOS__FFODESLV_CHECK
  assert( _pODESLV && nVar == _nPar+_nCst && nRes == _pODESLV->nf()*_nDir() );
#endif

  vVar->env()->lift( nRes, vRes, nVar, vVar );
}

} // end namescape mc

#endif
