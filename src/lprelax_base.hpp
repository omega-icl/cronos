// Copyright (C) 2015 Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

//TODO: 
//- [TO DO] Split LPRELAX_BASE into LPRELAX_CPLEX and LPRELAX_GUROBI

#ifndef MC__LPRELAX_BASE_HPP
#define MC__LPRELAX_BASE_HPP

#include <stdexcept>
#include <cassert>

#include "base_opt.hpp"
#include "polimage.hpp"
#include "cmodel.hpp"
#include "mctime.hpp"
#ifdef MC__USE_CPLEX
  #include "ilcplex/ilocplex.h"
#else
  #include "gurobi_c++.h"
#endif

#undef MC__LPRELAX_BASE_DEBUG
#undef MC__LPRELAX_BASE_TRACE
#undef MC__LPRELAX_BASE_SHOW_BREAKPTS

extern "C"{
  #include <fenv.h>
  int fedisableexcept( int );
}

namespace mc
{
//! @brief C++ base class for LP relaxation in complete search optimization
////////////////////////////////////////////////////////////////////////
//! mc::LPRELAX_BASE is a C++ base class for constructing and solving 
//! LP relaxations for use in complete search optimization. Relaxations
//! for the nonlinear or nonconvex participating terms are generated
//! using MC++ and the LP solvers are CPLEX and GUROBI.
////////////////////////////////////////////////////////////////////////
template< typename T>
class LPRELAX_BASE:
  public virtual BASE_OPT
{
  // Typedef's
#ifdef MC__USE_CPLEX
  typedef std::map< const PolVar<T>*, IloNumVar, lt_PolVar<T> > t_LPVar;
  typedef std::map< const PolCut<T>*, IloRange, lt_PolCut<T> > t_LPCut;
#else
  typedef std::map< const PolVar<T>*, GRBVar, lt_PolVar<T> > t_LPVar;
  typedef std::map< const PolCut<T>*, GRBConstr, lt_PolCut<T> > t_LPCut;
#endif

private:
#ifdef MC__USE_CPLEX
  //! @brief CPLEX environment for piecewise relaxation
  IloEnv* _ILOenv;
  //! @brief CPLEX model for piecewise relaxation
  IloModel* _ILOmodel;
  //! @brief CPLEX object for piecewise relaxation
  IloCplex* _ILOcplex;
#else
  //! @brief GUROBI environment for piecewise relaxation
  GRBEnv* _GRBenv;
  //! @brief GUROBI model for piecewise relaxation
  GRBModel* _GRBmodel;
#endif

  //! @brief map of LP variables in polyhedral image
  t_LPVar _LPvar;
  //! @brief map of LP cuts in polyhedral image
  t_LPCut _LPcut;
#ifdef MC__USE_CPLEX
  //! @brief LP objective in polyhedral image
  std::pair<bool,IloObjective> _LPobj;
  //! @brief LP constraint for incumbent cut in polyhedral image
  std::pair<bool,IloRange> _LPinc;
#else
  //! @brief GUROBI constraint for incumbent cut in polyhedral image
  std::pair<bool,GRBConstr> _LPinc;
#endif

protected:
  //! @brief Integrator status
  enum LP_STATUS{
     LP_OPTIMAL=0,	//!< Optimal solution found within tolerances
     LP_INFEASIBLE,	//!< Infeasible, but not unbounded
     LP_INFORUNBND,	//!< Infeasible or unbounded
     LP_OTHER		//!< Other status
  };

  //! @brief Polyhedral image environment
  PolImg<T> _POLenv;

  //! @brief Default option for LP solver
#ifdef MC__USE_CPLEX
  static const int LPALGO_DEFAULT = 0;
#else
  static const int LPALGO_DEFAULT = -1;
#endif

public:
  /** @defgroup LPRELAX_BASE Nonlinear Global Optimization using GUROBI
   *  @{
   */
  //! @brief Constructor
  LPRELAX_BASE()
    {
#ifdef MC__USE_CPLEX
      _ILOenv = new IloEnv;
      _ILOmodel = new IloModel(*_ILOenv);
      _ILOcplex = new IloCplex(*_ILOenv);
      _LPobj.first = false;
#else
      _GRBenv = new GRBEnv();
      _GRBmodel = new GRBModel(*_GRBenv);
#endif
      _LPinc.first = false;
    }

  //! @brief Destructor
  virtual ~LPRELAX_BASE()
    {
#ifdef MC__USE_CPLEX
      delete _ILOmodel;
      delete _ILOcplex;
      _ILOenv->end();
      delete _ILOenv;
#else
      delete _GRBmodel;
      delete _GRBenv;
#endif
    }

  //! @brief Value of DAG variable <a>X</a> after last LP optimization
  double get_variable
    ( const FFVar&X ) const
    {
      auto itp = _POLenv.Vars().find( const_cast<FFVar*>(&X) );
      auto itv = _LPvar.find( itp->second );
#ifdef MC__USE_CPLEX
      return _ILOcplex->getValue( itv->second );
#else
      //GRBVar XGRB = itv->second;
      //return XGRB.get(GRB_DoubleAttr_X);
      return itv->second.get(GRB_DoubleAttr_X);
#endif
    }

  //! @brief Optimal cost value after last LP optimization
  double get_objective() const
    {
#ifdef MC__USE_CPLEX
      return _ILOcplex->getObjValue();
#else
      return _GRBmodel->get( GRB_DoubleAttr_ObjVal );
#endif
    }

  //! @brief Status after last LP optimization
  LP_STATUS get_status() const
    {
#ifdef MC__USE_CPLEX
      switch( _ILOcplex->getStatus() ){
       case IloAlgorithm::Optimal:               return LP_OPTIMAL;
       case IloAlgorithm::Infeasible:            return LP_INFEASIBLE;
       case IloAlgorithm::InfeasibleOrUnbounded: return LP_INFORUNBND;
       default:                                  return LP_OTHER;
      }
#else
      switch( _GRBmodel->get( GRB_IntAttr_Status ) ){
       case GRB_OPTIMAL:      return LP_OPTIMAL;
       case GRB_INFEASIBLE:   return LP_INFEASIBLE;
       case GRB_INF_OR_UNBD:  return LP_INFORUNBND;
       default:               return LP_OTHER;
      }
#endif
    }

  //! @brief Pointer to last LP optimization model and solution
#ifdef MC__USE_CPLEX
  const IloModel* get_relaxed_model() const
    { return _ILOmodel; }
#else
  const GRBModel* get_relaxed_model() const
    { return _GRBmodel; }
#endif
  /** @} */

protected:
  //! @brief Solve LP model
  template <typename OPT, typename STAT> void _solve_LPmodel
    ( const OPT&options, STAT&stats, const std::vector<FFVar>&var );
  //! @brief Set objective in LP relaxation
  void _set_LPrelax
    ( const PolVar<T>*pObj, const t_OBJ&tObj, const bool feastest=false );
  //! @brief Set objective in LP relaxation
  void _set_LPrelax
    ( const unsigned nObj, const PolVar<T>*pObj, const double*cObj,
      const t_OBJ&tObj, const bool feastest=false );
  //! @brief Set parameter bound contracting objective in LP model
  void _set_LPcontract
    ( const PolVar<T>*pVar, const bool uplo, const PolVar<T>*pObj,
     const double*inc, const t_OBJ&tObj, const bool feastest=false );
  //! @brief Set variables and cuts in LP model
  void _set_LPcuts
    ();
 
  //! @brief Function computing Hausdorff distance between intervals
  template <typename U> static double _dH
    ( const U&X, const U&Y );
  //! @brief Function computing Hausdorff distance between interval vectors
  template <typename U> static double _reducrel
    ( const unsigned n, const U*Xred, const U*X );

  //! @brief Value of PolImg variable <a>X</a> after last LP optimization
  double _get_variable
    ( const PolVar<T>&X ) const
    {
      auto itv = _LPvar.find( const_cast<PolVar<T>*>(&X) );
#ifdef MC__USE_CPLEX
      return _ILOcplex->getValue( itv->second );
#else
      //GRBVar XGRB = itv->second;
      //return XGRB.get(GRB_DoubleAttr_X);
      return itv->second.get(GRB_DoubleAttr_X);
#endif
    }

private:
  //! @brief Set options of LP model
  template <typename OPT> void _set_LPoptions
    ( const OPT&options );
  //! @brief Reset LP model
  void _reset_LPmodel
    ();
  //! @brief Append polyhedral image variable/auxiliary into LP model
  std::pair<typename t_LPVar::iterator,bool> _set_LPvar
    ( const PolVar<T>*pVar );
  //! @brief Append polyhedral image cut into LP model
  void _set_LPcut
    ( const PolCut<T>*pCut );

  //! @brief Private methods to block default compiler methods
  LPRELAX_BASE
    (const LPRELAX_BASE&);
  LPRELAX_BASE& operator=
    (const LPRELAX_BASE&);
};

template <typename T> template <typename OPT, typename STAT> inline void
LPRELAX_BASE<T>::_solve_LPmodel
( const OPT&options, STAT&stats, const std::vector<FFVar>&var )
{
  stats.tLPSOL -= cpuclock();
  _set_LPoptions( options );
#ifdef MC__USE_CPLEX
  if( options.MIPFILE != "" ) _ILOcplex->exportModel( options.MIPFILE.c_str() );
  //{ int dum; std::cout << "PAUSED"; std::cin >> dum; }
  //_time = -_cplex->getCplexTime();
  _ILOcplex->solve();
  //_time += _cplex->getCplexTime();
#else
  _GRBmodel->update();
  if( options.MIPFILE != "" ) _GRBmodel->write( options.MIPFILE );
  fedisableexcept(FE_ALL_EXCEPT);
  _GRBmodel->optimize();
#endif
  stats.tLPSOL += cpuclock();
  stats.nLPSOL++;
#ifdef MC__LPRELAX_BASE_DEBUG
  std::cout << "LP solution complete\n";
  std::cout << std::scientific << std::setprecision(4);
  std::cout << "  fopt = " << get_objective() << std::endl;
  unsigned ivar = 0;
  for( auto itv=_var.begin(); itv!=_var.end(); ++itv, ivar++ )
    std::cout << "popt" << ivar << " = " << get_variable( *itv ) << std::endl;
  { int dum; std::cin >> dum; }
#endif
}

template <typename T> template <typename OPT> inline void
LPRELAX_BASE<T>::_set_LPoptions
( const OPT&options )
{
#ifdef MC__USE_CPLEX
  // CPLEX options
  _ILOcplex->extract(*_ILOmodel);
  _ILOcplex->setWarning( options.MIPDISPLAY? std::cout: _ILOenv->getNullStream() );
  _ILOcplex->setOut( options.MIPDISPLAY? std::cout: _ILOenv->getNullStream() );
  _ILOcplex->setParam( IloCplex::RootAlg, options.LPALGO );
  _ILOcplex->setParam( IloCplex::EpOpt,   options.LPOPTIMTOL );
  _ILOcplex->setParam( IloCplex::EpRHS,   options.LPFEASTOL );
#else
  // Gurobi options
  _GRBmodel->getEnv().set( GRB_IntParam_Method,            options.LPALGO );
  _GRBmodel->getEnv().set( GRB_DoubleParam_FeasibilityTol, options.LPFEASTOL );
  _GRBmodel->getEnv().set( GRB_DoubleParam_OptimalityTol,  options.LPOPTIMTOL );
  _GRBmodel->getEnv().set( GRB_DoubleParam_MIPGap,         options.MIPRELGAP );
  _GRBmodel->getEnv().set( GRB_DoubleParam_MIPGapAbs,      options.MIPABSGAP );
  _GRBmodel->getEnv().set( GRB_IntParam_Presolve,          options.LPPRESOLVE  );
  _GRBmodel->getEnv().set( GRB_IntParam_OutputFlag,        options.MIPDISPLAY );
  _GRBmodel->getEnv().set( GRB_DoubleParam_PreSOS2BigM,    options.PRESOS2BIGM );
  _GRBmodel->getEnv().set( GRB_IntParam_DualReductions,    0 ); // In order to avoid INF_OR_UNBD status
#endif
}

template <typename T> inline void
LPRELAX_BASE<T>::_reset_LPmodel
()
{
#ifdef MC__USE_CPLEX
  delete _ILOmodel; delete _ILOcplex;
  _ILOenv->end();
  delete _ILOenv;
  _ILOenv = new IloEnv;
  _ILOmodel = new IloModel( *_ILOenv );
  _ILOcplex = new IloCplex( *_ILOenv );
  _LPobj.first = false;
#else
  delete _GRBmodel;
  _GRBmodel = new GRBModel( *_GRBenv );
#endif
  _LPvar.clear(); _LPcut.clear(); _LPinc.first = false;
}

template <typename T> inline void
LPRELAX_BASE<T>::_set_LPcuts
()
{
  _reset_LPmodel();
  for( auto itv=_POLenv.Vars().begin(); itv!=_POLenv.Vars().end(); ++itv )
    _set_LPvar( itv->second );
    //if( itv->second->cuts() ) _set_LPvar( itv->second );
  for( auto itv=_POLenv.Aux().begin(); itv!=_POLenv.Aux().end(); ++itv )
    _set_LPvar( *itv );
#ifndef MC__USE_CPLEX
  _GRBmodel->update();
#endif
  for( auto itc=_POLenv.Cuts().begin(); itc!=_POLenv.Cuts().end(); ++itc )
    _set_LPcut( *itc );
}

template <typename T> inline void
LPRELAX_BASE<T>::_set_LPrelax
( const PolVar<T>*pObj, const t_OBJ&tObj, const bool feastest )
{
  // Set objective
  auto jtobj = _LPvar.find( pObj );
#ifdef MC__USE_CPLEX
  if( _LPobj.first ) _ILOmodel->remove( _LPobj.second );
  if( !feastest )
    switch( tObj ){
      case MIN: _LPobj.second = IloMinimize( *_ILOenv, jtobj->second ); break;
      case MAX: _LPobj.second = IloMaximize( *_ILOenv, jtobj->second ); break;
    }
  else
    _LPobj.second = IloMinimize( *_ILOenv, jtobj->second );
  _ILOmodel->add( _LPobj.second ); _LPobj.first = true;
#else
  _GRBmodel->setObjective( GRBLinExpr( jtobj->second,  1. ) );
  if( !feastest )
    switch( tObj ){
      case MIN: _GRBmodel->set( GRB_IntAttr_ModelSense,  1 ); break;
      case MAX: _GRBmodel->set( GRB_IntAttr_ModelSense, -1 ); break;
    }
  else
    _GRBmodel->set( GRB_IntAttr_ModelSense,  1 );
#endif

  // Remove incumbent constraint (if any)
  if( _LPinc.first ){
#ifdef MC__USE_CPLEX
    _ILOmodel->remove( _LPinc.second );
#else
    _GRBmodel->remove( _LPinc.second );
#endif
    _LPinc.first = false;
  }
}

template <typename T> inline void
LPRELAX_BASE<T>::_set_LPrelax
( const unsigned nObj, const PolVar<T>*pObj, const double*cObj,
  const t_OBJ&tObj, const bool feastest )
{
  // Set objective
#ifdef MC__USE_CPLEX
  auto ObjExpr = IloNumExpr( *_ILOenv );
  for( unsigned i=0; i<nObj; i++ ){
    auto it = _LPvar.find( &pObj[i] );
    assert( it != _LPvar.end() );
    ObjExpr += it->second * cObj[i];
  }
  if( _LPobj.first ) _ILOmodel->remove( _LPobj.second );
  if( !feastest )
    switch( tObj ){
      case MIN: _LPobj.second = IloMinimize( *_ILOenv, ObjExpr ); break;
      case MAX: _LPobj.second = IloMaximize( *_ILOenv, ObjExpr ); break;
    }
  else
    _LPobj.second = IloMinimize( *_ILOenv, ObjExpr );
  _ILOmodel->add( _LPobj.second ); _LPobj.first = true;
#else
  auto ObjExpr = GRBLinExpr( 0. );
  for( unsigned i=0; i<nObj; i++ ){
    auto it = _LPvar.find( &pObj[i] );
    assert( it != _LPvar.end() );
    ObjExpr += it->second * cObj[i];
  }
  _GRBmodel->setObjective( ObjExpr );
  if( !feastest )
    switch( tObj ){
      case MIN: _GRBmodel->set( GRB_IntAttr_ModelSense,  1 ); break;
      case MAX: _GRBmodel->set( GRB_IntAttr_ModelSense, -1 ); break;
    }
  else
    _GRBmodel->set( GRB_IntAttr_ModelSense,  1 );
#endif

  // Remove incumbent constraint (if any)
  if( _LPinc.first ){
#ifdef MC__USE_CPLEX
    _ILOmodel->remove( _LPinc.second );
#else
    _GRBmodel->remove( _LPinc.second );
#endif
    _LPinc.first = false;
  }
}

template <typename T> inline void
LPRELAX_BASE<T>::_set_LPcontract
( const PolVar<T>*pVar, const bool uplo, const PolVar<T>*pObj, 
  const double*inc, const t_OBJ&tObj, const bool feastest )
{
  // Add lower/upper parameter bound objective
  auto jtobj = _LPvar.find( pVar );
#ifdef MC__USE_CPLEX
  if( _LPobj.first ) _ILOmodel->remove( _LPobj.second );
  switch( (int)uplo ){
    case false: _LPobj.second = IloMinimize( *_ILOenv, jtobj->second ); break;
    case true:  _LPobj.second = IloMaximize( *_ILOenv, jtobj->second ); break;
  }
  _ILOmodel->add( _LPobj.second ); _LPobj.first = true;
#else
  _GRBmodel->setObjective( GRBLinExpr( jtobj->second,  1. ) );
  switch( uplo ){
    case false: _GRBmodel->set( GRB_IntAttr_ModelSense,  1 ); break;
    case true:  _GRBmodel->set( GRB_IntAttr_ModelSense, -1 ); break;
  }
#endif

  // Add/update incumbent or backoff constraint
  if( _LPinc.first ){
#ifdef MC__USE_CPLEX
    _ILOmodel->remove( _LPinc.second );
#else
    _GRBmodel->remove( _LPinc.second );
#endif
    _LPinc.first = false;
  }
  if( !inc ) return;
  jtobj = _LPvar.find( pObj );
#ifdef MC__USE_CPLEX
  IloExpr lhsinc( *_ILOenv ); lhsinc = jtobj->second;
#endif
  if( !feastest ){
    switch( tObj ){
#ifdef MC__USE_CPLEX
      case MIN: _LPinc.second = (lhsinc <= *inc); _ILOmodel->add( _LPinc.second ); break;
      case MAX: _LPinc.second = (lhsinc >= *inc); _ILOmodel->add( _LPinc.second ); break;
#else
      case MIN: _LPinc.second = _GRBmodel->addConstr( GRBLinExpr( jtobj->second,  1. ),
        GRB_LESS_EQUAL, *inc ); break;
      case MAX: _LPinc.second = _GRBmodel->addConstr( GRBLinExpr( jtobj->second,  1. ),
        GRB_GREATER_EQUAL, *inc ); break;
#endif
    }
  }
  else
#ifdef MC__USE_CPLEX
    _LPinc.second = (lhsinc <= *inc); _ILOmodel->add( _LPinc.second );
#else
    _LPinc.second = _GRBmodel->addConstr( GRBLinExpr( jtobj->second,  1. ),
        GRB_LESS_EQUAL, *inc );
#endif
  _LPinc.first = true;
}

template <typename T> inline std::pair<typename LPRELAX_BASE<T>::t_LPVar::iterator,bool>
LPRELAX_BASE<T>::_set_LPvar
( const PolVar<T>*pVar )
{
#ifdef MC__USE_CPLEX
  IloNumVar var;
#else
  GRBVar var;
#endif
  switch( pVar->id().first ){
    case PolVar<T>::VARCONT:
    case PolVar<T>::AUXCONT:
#ifdef MC__USE_CPLEX
      var = IloNumVar( *_ILOenv, Op<T>::l(pVar->range()), Op<T>::u(pVar->range()),
       ILOFLOAT, pVar->name().c_str() );
      _ILOmodel->add( var );
#else
      var = _GRBmodel->addVar( Op<T>::l(pVar->range()), Op<T>::u(pVar->range()),
        0., GRB_CONTINUOUS, pVar->name() );
#endif
      break;
    case PolVar<T>::VARINT:
    case PolVar<T>::AUXINT:
      if( isequal( Op<T>::l(pVar->range()), 0. ) && isequal( Op<T>::u(pVar->range()), 1. ) ){
#ifdef MC__USE_CPLEX
        var = IloNumVar( *_ILOenv, 0., 1., ILOBOOL, pVar->name().c_str() );
        _ILOmodel->add( var );
#else
        var = _GRBmodel->addVar( 0., 1., 0., GRB_BINARY, pVar->name() );
#endif
      }
      else{
#ifdef MC__USE_CPLEX
        var = IloNumVar( *_ILOenv, Op<T>::l(pVar->range()), Op<T>::u(pVar->range()),
         ILOINT, pVar->name().c_str() );
        _ILOmodel->add( var );
#else
        var = _GRBmodel->addVar( Op<T>::l(pVar->range()), Op<T>::u(pVar->range()),
          0., GRB_INTEGER, pVar->name() );
#endif
      }
      break;

    default:
      throw std::runtime_error("Invalid auxiliary variable type");
  }

  return _LPvar.insert( std::make_pair( pVar, var ) );
}

template <typename T> inline void
LPRELAX_BASE<T>::_set_LPcut
( const PolCut<T>*pCut )
{
#ifdef MC__USE_CPLEX
  bool isSOS = ( pCut->type() == PolCut<T>::SOS1
              || pCut->type() == PolCut<T>::SOS2 );
  IloNumVarArray VarSOS( *_ILOenv, (isSOS? pCut->nvar(): 0) );
  IloNumArray WeiSOS( *_ILOenv, (isSOS? pCut->nvar(): 0) );

  IloExpr lhs( *_ILOenv );
  for( unsigned k=0; k<pCut->nvar(); k++ ){
    auto ivar = _LPvar.find( pCut->var()+k );
    if( ivar==_LPvar.end() ) throw std::runtime_error("variable not found");
    if( pCut->type() == PolCut<T>::SOS1 || pCut->type() == PolCut<T>::SOS2 ){
      VarSOS[k] = ivar->second; WeiSOS[k] = pCut->coef()[k];
    }
    else
      lhs += pCut->coef()[k] * ivar->second;
  }

  IloRange ctr;
  try{
    switch( pCut->type() ){
      case PolCut<T>::SOS1:
        _ILOmodel->add( IloSOS1( *_ILOenv, VarSOS, WeiSOS ) );
        ctr = ( lhs == pCut->rhs() ); break;
      case PolCut<T>::SOS2:
        _ILOmodel->add( IloSOS2( *_ILOenv, VarSOS, WeiSOS ) );
        ctr = ( lhs == pCut->rhs() ); break;
      case PolCut<T>::EQ:
        ctr = ( lhs == pCut->rhs() ); break;
      case PolCut<T>::LE:
        ctr = ( lhs <= pCut->rhs() ); break;
      case PolCut<T>::GE:
        ctr = ( lhs >= pCut->rhs() ); break;
    }
    _ILOmodel->add( ctr );
    _LPcut.insert( std::make_pair( pCut, ctr ) );
  }
  catch(IloException& e){
#ifdef MC__LPRELAX_BASE_DEBUG
    std::cout << "IloException Caught - Error code = " << e.getMessage() << std::endl;
#endif
  }

#else
  GRBVar* VarSOS = 0;
  double* WeiSOS = 0;
  int TypSOS = 0;
  if( pCut->type() == PolCut<T>::SOS1 || pCut->type() == PolCut<T>::SOS2 ){
    VarSOS = new GRBVar[pCut->nvar()];
    WeiSOS = new double[pCut->nvar()];
    TypSOS = ( pCut->type() == PolCut<T>::SOS1? GRB_SOS_TYPE1: GRB_SOS_TYPE2 );
  }

  GRBLinExpr lhs;
  for( unsigned k=0; k<pCut->nvar(); k++ ){
    auto ivar = _LPvar.find( pCut->var()+k );
    if( ivar==_LPvar.end() ) throw std::runtime_error("variable not found");
    GRBVar&Var = ivar->second;
    if( pCut->type() == PolCut<T>::SOS1 || pCut->type() == PolCut<T>::SOS2 ){
      VarSOS[k] = Var; WeiSOS[k] = pCut->coef()[k];//(double)k/(double)pCut->nvar();//
    }
    else
      lhs += GRBLinExpr( Var, pCut->coef()[k] );
  }

  GRBConstr ctr;
  try{
    switch( pCut->type() ){
      case PolCut<T>::SOS1:
      case PolCut<T>::SOS2:
        _GRBmodel->addSOS( VarSOS, WeiSOS, pCut->nvar(), TypSOS );
        delete [] VarSOS;
        delete [] WeiSOS;
        break;
      case PolCut<T>::EQ:
        ctr = _GRBmodel->addConstr( lhs, GRB_EQUAL, pCut->rhs() );
        break;
      case PolCut<T>::LE:
        ctr = _GRBmodel->addConstr( lhs, GRB_LESS_EQUAL, pCut->rhs() );
        break;
      case PolCut<T>::GE:
        ctr = _GRBmodel->addConstr( lhs, GRB_GREATER_EQUAL, pCut->rhs() );
        break;
    }
    _LPcut.insert( std::make_pair( pCut, ctr ) );
  }
  catch(GRBException& e){
#ifdef MC__LPRELAX_BASE_DEBUG
    std::cout << "GRBException Caught - Error code = " << e.getErrorCode() << std::endl;
    std::cout << e.getMessage() << std::endl;
#endif
  }
#endif
}

template <typename T> template <typename U> inline double
LPRELAX_BASE<T>::_dH
( const U&X, const U&Y )
{
  return std::max( std::fabs(Op<U>::l(X)-Op<U>::l(Y)),
                   std::fabs(Op<U>::u(X)-Op<U>::u(Y)) );
}

template <typename T> template <typename U> inline double
LPRELAX_BASE<T>::_reducrel
( const unsigned n, const U*Xred, const U*X )
{
  double drel = 0.;
  for( unsigned ip=0; ip<n; ip++ )
    drel = std::max( drel, _dH( Xred[ip], X[ip] ) / Op<T>::diam( X[ip] ) );
  return drel;
}

} // end namescape mc

#endif
