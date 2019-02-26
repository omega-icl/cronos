// Copyright (C) 2014 Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

#ifndef MC__BASE_DO_HPP
#define MC__BASE_DO_HPP

#include "base_opt.hpp"
#include "base_de.hpp"

#undef MC__BASE_DO__DEBUG

namespace mc
{
//! @brief C++ structure for comparing variables in a dependency map
////////////////////////////////////////////////////////////////////////
//! mc::lt_DepVar is a C++ structure for ordering variables in a
//! dependency map.
////////////////////////////////////////////////////////////////////////
struct lt_DepVar
////////////////////////////////////////////////////////////////////////
{
  typedef std::pair< unsigned, int > key_DepVar;

  bool operator()
    ( const key_DepVar& Var1,
      const key_DepVar& Var2 ) const
    {
      // Order dependent variables w.r.t. their stages and variable indices
      if( Var1.first < Var2.first ) return true;
      if( Var1.first > Var2.first ) return false;
      // Order RLT variables w.r.t. to their index next
      return( Var1.second < Var2.second );
    }
};

//! @brief C++ base class for definition of dynamic optimization problems
////////////////////////////////////////////////////////////////////////
//! mc::BASE_DO is a C++ base class for definition of the variables,
//! objective and constraints participating in dynamic optimization
//! problems.
////////////////////////////////////////////////////////////////////////
class BASE_DO:
  public virtual BASE_OPT,
  public virtual BASE_DE
{
public:
  typedef std::pair< unsigned, unsigned > key_DepVar;
  typedef std::map< const key_DepVar, FFVar, lt_DepVar > t_DepVar;

  //! @brief Class constructor
  BASE_DO()
    : BASE_OPT(), BASE_DE()
    {}

  //! @brief Class destructor
  virtual ~BASE_DO()
    {}
  //! @brief Get constraint functions
  const std::tuple< std::vector<t_CTR>, std::vector< std::map<unsigned,FFVar> >,
                    std::vector<FFVar>, std::vector<bool> >& ctr
    ()
    const
    { return _ctr; }

  //! @brief Reset constraint functions
  void reset_ctr
    ()
    { std::get<0>(_ctr).clear(); std::get<1>(_ctr).clear(); std::get<2>(_ctr).clear();
      std::get<3>(_ctr).clear(); }

  //! @brief Add constraint function
  void add_ctr
    ( const t_CTR type, const std::map<unsigned,FFVar>&ctrmap, const bool is_redundant=false )
    { std::get<0>(_ctr).push_back( type );
      std::get<1>(_ctr).push_back( ctrmap );
      std::get<2>(_ctr).push_back( FFVar( ctrmap.begin()->second.dag() ) );
      std::get<3>(_ctr).push_back( is_redundant ); }
  void add_ctr
    ( const t_CTR type, const std::pair<unsigned,FFVar>&ctr, const bool is_redundant=false )
    { std::get<0>(_ctr).push_back( type );
      std::map<unsigned,FFVar> ctrmap; ctrmap.insert( ctr );
      std::get<1>(_ctr).push_back( ctrmap );
      std::get<2>(_ctr).push_back( FFVar( ctr.second.dag() ) );
      std::get<3>(_ctr).push_back( is_redundant ); }
  void add_ctr
    ( const t_CTR type, const FFVar&ctr, const bool is_redundant=false )
    { std::get<0>(_ctr).push_back( type );
      std::map<unsigned,FFVar> ctrmap; ctrmap.insert( std::make_pair(0,ctr) );
      std::get<1>(_ctr).push_back( ctrmap );
      std::get<2>(_ctr).push_back( FFVar( ctr.dag() ) );
      std::get<3>(_ctr).push_back( is_redundant ); }

  //! @brief Get objective function
  const std::tuple< std::vector<t_OBJ>, std::vector< std::map<unsigned,FFVar> >,
                    std::vector<FFVar> >& obj
    ()
    const
    { return _obj; }

  //! @brief Set objective function
  void set_obj
    ( const t_OBJ type, const std::map<unsigned,FFVar>&objmap )
    { std::get<0>(_obj).clear(); std::get<0>(_obj).push_back( type );
      std::get<1>(_obj).clear(); std::get<1>(_obj).push_back( objmap );
      std::get<2>(_obj).clear(); std::get<2>(_obj).push_back( FFVar( objmap.begin()->second.dag() ) ); }
  void set_obj
    ( const t_OBJ type, const std::pair<unsigned,FFVar>&obj )
    { std::get<0>(_obj).clear(); std::get<0>(_obj).push_back( type );
      std::map<unsigned,FFVar> objmap; objmap.insert( obj );
      std::get<1>(_obj).clear(); std::get<1>(_obj).push_back( objmap );
      std::get<2>(_obj).clear(); std::get<2>(_obj).push_back( FFVar( obj.second.dag() ) ); }
  void set_obj
    ( const t_OBJ type, const FFVar&obj )
    { std::get<0>(_obj).clear(); std::get<0>(_obj).push_back( type );
      std::map<unsigned,FFVar> objmap; objmap.insert( std::make_pair(0,obj) );
      std::get<1>(_obj).clear(); std::get<1>(_obj).push_back( objmap );
      std::get<2>(_obj).clear(); std::get<2>(_obj).push_back( FFVar( obj.dag() ) ); }

  //! @brief Copy equations
  void set
    ( const BASE_DO&pb )
    { BASE_DE::set(pb); _ctr = pb._ctr; _obj = pb._obj; }

  //! @brief Generate dependency information in IVP and functions
  bool set_depend
    ();

protected:
  //! @brief constraints (types, constraint variables, constraint multipliers, redundancy status)
  std::tuple< std::vector<t_CTR>, std::vector< std::map<unsigned,FFVar> >, std::vector<FFVar>,
              std::vector<bool> > _ctr;

  //! @brief objective (type, cost variable, cost multiplier)
  std::tuple< std::vector<t_OBJ>, std::vector< std::map<unsigned,FFVar> >,
              std::vector<FFVar> > _obj;

  //! @brief constraints (types, NCO variables), including dependent equations
  std::tuple< std::vector<t_CTR>, std::vector<FFVar> > _nco;

  //! @brief Get 1st-order necessary conditions for optimality (NCO)
  const std::tuple< std::vector<t_CTR>, std::vector<FFVar> >& nco
    ()
    const
    { return _nco; }

  //! @brief Reset NCO
  void reset_nco
    ()
    { std::get<0>(_nco).clear(); std::get<1>(_nco).clear(); }

  //! @brief Define NCO based on cost/constraint functions
  template<typename T>
  void set_nco
    ( const T*P, const unsigned*tvar );

  //! @brief Define state functions based on cost/constraint functions
  bool set_function
    ();

  //! @brief Test dependency w.r.t. state variables
  bool state_depend
    ( const unsigned ifct );

  //! @brief number of time stages
  unsigned _ns;

  //! @brief vector of function
  std::vector<FFVar> _vec_FCT_AGGR;

  //! @brief vector of parameter-dependent functions (no state dependency)
  std::vector<FFVar> _vec_FCT_PAR;

  //! @brief indices of state-dependent functions
  std::set<unsigned> _set_FCT_IVP;

  //! @brief vector of state-dependent functions (assumes additive contributions) in each stage of IVP-DAE system
  std::vector<FFVar> _vec_FCT_IVP_STG;

  //! @brief vector of variables mapping state-dependent functions
  std::vector<FFVar> _vec_VAR_FCT_IVP;

  //! @brief vector of variables mapping state-dependent function derivatives
  std::vector<FFVar> _vec_VAR_DFCT_IVP;

  //! @brief indices of parameters participating in state-dependent functions
  std::set<unsigned> _set_PAR_FCT_IVP;

  //! @brief parameters participating in state-dependent functions
  std::vector<FFVar> _vec_PAR_FCT_IVP;

  //! @brief indices of parameters participating in parametric IVP
  std::set<unsigned> _set_PAR_STA;

  //! @brief parameters participating in parametric IVP
  std::vector<FFVar> _vec_PAR_STA;

  //! @brief depency map for state-dependent functions
  t_DepVar _map_STA_FCT_IVP;

private:
  //! @brief Return number of time stages in cost/constraint functions
  bool _set_ns
    ();

  //! @brief vector holding function contributions
  std::vector<FFVar> _pFCT;

  //! @brief Build-up state function
  void _build_function
    ( const unsigned is, std::map<unsigned,FFVar>::iterator&jt,
      const unsigned ifct );

  //! @brief Private using directive to block method
  using BASE_DE::set_function;

  //! @brief Private methods to block default compiler methods
  BASE_DO(const BASE_DO&);
  BASE_DO& operator=(const BASE_DO&);
};

inline bool
BASE_DO::_set_ns
()
{
  // #stages in cost
  _ns = ( std::get<1>(_obj).empty()? 0: std::get<1>(_obj)[0].rbegin()->first+1 );
  // #stages in constraints
  auto it = std::get<1>(_ctr).begin();
  for( ; it!=std::get<1>(_ctr).end(); ++it ){
    const unsigned is = it->rbegin()->first+1;
    if( _ns < is ) _ns = is;
  }
  return( _ns <= _nsmax );
}

inline bool
BASE_DO::state_depend
( const unsigned ifct )
{
  assert( ifct < _nf );
  for( auto it=_map_STA_FCT_IVP.begin(); it!=_map_STA_FCT_IVP.end(); ++it )
    if( _vec_FCT_AGGR[ifct].dep().dep( it->second.id().second ).first )
      return true;
  return false;
}

inline bool
BASE_DO::set_depend
()
{
  if( !set_function() || !BASE_DE::set_depend( _ns ) ) return false;
  _set_FCT_IVP.clear();
  for( unsigned ifct=0; ifct<_nf; ifct++ ){
#ifdef MC__BASE_DO__DEBUG
    std::cout << "FCT[" << ifct << "]: " << state_depend(ifct) << std::endl;
#endif
    if( state_depend( ifct ) ) _set_FCT_IVP.insert( ifct );
  }

  // Parametric dependencies in IVP
  _set_PAR_STA.clear();
  for( auto it=_map_STA_FCT_IVP.begin(); it!=_map_STA_FCT_IVP.end(); ++it )
    for( unsigned ip=0; ip<_np; ip++ ){
      bool dep = it->first.second < _nx?
                 _depX[_nx*it->first.first+it->first.second].dep( _pP[ip].id().second ).first:
                 _depQ[_nq*it->first.first+it->first.second-_nx].dep( _pP[ip].id().second ).first;
      if( dep ) _set_PAR_STA.insert( ip );
    }

  // Parameters participating in IVP
  _vec_PAR_STA.clear();
  for( auto itp=_set_PAR_STA.begin(); itp!=_set_PAR_STA.end(); ++itp )
    _vec_PAR_STA.push_back( _pP[*itp] );

  // Create vectors of functions with function "variables"
  _vec_VAR_FCT_IVP = _vec_FCT_AGGR;
  unsigned ifct = 0;
  for( auto it=_set_FCT_IVP.begin(); it!=_set_FCT_IVP.end(); ++it, ifct++ )
    _vec_VAR_FCT_IVP[*it] = FFVar( _pDAG );

  // Function dependencies
  _vec_FCT_IVP_STG.resize( _ns*_set_FCT_IVP.size() );
  _vec_FCT_PAR.clear();
  _set_PAR_FCT_IVP = _set_PAR_STA;//.clear();
  for( unsigned ifct=0, ista=0; ifct<_nf; ifct++ ){
    auto itc = _set_FCT_IVP.find( ifct );
    // State-dependent functions
    if( itc != _set_FCT_IVP.end() ){
      for( unsigned is=0; is<_ns; is++ )
        _vec_FCT_IVP_STG[is*_set_FCT_IVP.size()+ista] = _vFCT[is][ifct];
      ista++;
      // Append parametric dependencies
      for( unsigned ip=0; ip<_np; ip++ ){
        bool dep = _depF[ifct].dep( _pP[ip].id().second ).first;
        if( dep ) _set_PAR_FCT_IVP.insert( ip );
      }
    }
    // Parameter-dependent functions
    else{
      _vec_FCT_PAR.push_back( _vec_FCT_AGGR[ifct] );
    }
  }

  // Parameters participating in state-dependent functions
  _vec_PAR_FCT_IVP.clear();
  for( auto itp=_set_PAR_FCT_IVP.begin(); itp!=_set_PAR_FCT_IVP.end(); ++itp )
    _vec_PAR_FCT_IVP.push_back( _pP[*itp] );

  return true;
}

inline void
BASE_DO::_build_function
( const unsigned is, std::map<unsigned,FFVar>::iterator&jt,
  const unsigned ifct )
{
  // Detect state/quadrature dependencies
  std::vector<FFVar> termsta, termdep;
  auto jd = jt->second.dep().dep().begin();
  for( ; jd!=jt->second.dep().dep().end(); ++jd ){
    auto jx = _ndxX.find(jd->first);
    auto jq = _ndxQ.find(jd->first);
    if( jx == _ndxX.end() && jq == _ndxQ.end() )
      continue; // not a state or quadrature variable
    const unsigned ixq = ( jx!=_ndxX.end()? jx->second: _nx+jq->second );
    const FFVar&pXQ = ( jx!=_ndxX.end()? _pX[jx->second]: _pQ[jq->second] );
    auto ndxq = std::make_pair(is,ixq); 
    if( _map_STA_FCT_IVP.find( ndxq ) == _map_STA_FCT_IVP.end() )
      _map_STA_FCT_IVP[ndxq] = FFVar( _pDAG );
    termsta.push_back( pXQ );
    termdep.push_back( _map_STA_FCT_IVP[ndxq] );
  }
  // Append objective term to _vec_FCT_AGGR[ifct]
  if( termdep.size() ){
    const mc::FFVar* term = _pDAG->compose( 1, &jt->second, termdep.size(),
                                            termsta.data(), termdep.data() );
    _vec_FCT_AGGR[ifct] += term[0];
    delete[] term;
  }
  else{
    _vec_FCT_AGGR[ifct] += jt->second;
  }
}

inline bool
BASE_DO::set_function
()
{
  // Resize function vector
  if( !_set_ns() ) return false;
  _nf = std::get<0>(_obj).size() + std::get<0>(_ctr).size();
  _pFCT.resize( _ns*_nf );
  _vec_FCT_AGGR.assign( _nf, FFVar(0) );
  _map_STA_FCT_IVP.clear();

  // Append cost contributions
  unsigned ifct = 0;
  auto itobj = std::get<1>(_obj).begin();
  for( ; itobj!=std::get<1>(_obj).end(); ++itobj, ++ifct ){
    unsigned is = 0;
    auto jt = itobj->begin();
    for( ; jt!=itobj->end(); ++jt ){
      for( ; is<jt->first; ++is ) _pFCT[_nf*is+ifct] = 0;
      _pFCT[_nf*jt->first+ifct] = jt->second;
      _build_function( is, jt, ifct );
      ++is;
    }
    for( ; is<_ns; ++is ) _pFCT[_nf*is+ifct] = 0;
#ifdef MC__BASE_DO__DEBUG
    std::cout << "_vec_FCT_AGGR[" << ifct << "] <-- " << _vec_FCT_AGGR[ifct].dep() << std::endl;
    //_pDAG->output( _pDAG->subgraph( 1, &_vec_FCT_AGGR[ifct] ) );
    //{ int dum; std::cout << "PAUSED "; std::cin >> dum; }
#endif
  }

  // Append constraint contributions
  auto itctr = std::get<1>(_ctr).begin();
  for( ; itctr!=std::get<1>(_ctr).end(); ++itctr, ++ifct ){
    unsigned is = 0;
    auto jt = itctr->begin();
    for( ; jt!=itctr->end(); ++jt ){
      for( ; is<jt->first; ++is ) _pFCT[_nf*is+ifct] = 0;
      _pFCT[_nf*jt->first+ifct] = jt->second;
      _build_function( is, jt, ifct );
      ++is;
    }
    for( ; is<_ns; ++is ) _pFCT[_nf*is+ifct] = 0;
#ifdef MC__BASE_DO__DEBUG
    std::cout << "_vec_FCT_AGGR[" << ifct << "] <-- " << _vec_FCT_AGGR[ifct].dep() << std::endl;
    //_pDAG->output( _pDAG->subgraph( 1, &_vec_FCT_AGGR[ifct] ) );
    //{ int dum; std::cout << "PAUSED "; std::cin >> dum; }
#endif
  }

#ifdef MC__BASE_DO__DEBUG
  for( auto itv=_map_STA_FCT_IVP.begin(); itv!=_map_STA_FCT_IVP.end(); ++itv )
    std::cout << "_map_STA_FCT_IVP[" << itv->first.first << "," << itv->first.second << "] = "
              << itv->second << std::endl;
  for( unsigned is=0; is<_ns; is++ )
    for( unsigned ifct=0; ifct<_nf; ifct++ )
      std::cout << "FCT[" << is << "][" << ifct << "] = " << _pFCT[_nf*is+ifct] << std::endl;
  int dum;  std::cin >> dum; 
#endif

  // vector of pointers
  _vFCT.clear();
  for( unsigned i=0; i<_ns; i++ ) _vFCT.push_back( _pFCT.data()+i*_nf );

  return true;
}

template<typename T>
inline void
BASE_DO::set_nco
( const T*P, const unsigned*tvar )
{
  reset_nco();
  BASE_DE::set_boundmultiplier();
  const unsigned nc = std::get<0>(_ctr).size();

  // Create vectors of functions and state-dependent function derivatives
  _vec_VAR_DFCT_IVP.resize( _set_FCT_IVP.size()*_vec_PAR_FCT_IVP.size() );
  unsigned ifct = 0;
  for( auto it=_set_FCT_IVP.begin(); it!=_set_FCT_IVP.end(); ++it, ifct++ )
    for( unsigned ip=0; ip<_vec_PAR_FCT_IVP.size(); ++ip )
      _vec_VAR_DFCT_IVP[ifct+ip*_set_FCT_IVP.size()] = FFVar( _pDAG );

  // Expressions of Lagrangian function and multiplier scaling
  FFVar lagr = 0., scal = -1.;
  switch( std::get<0>(_obj)[0] ){
   case BASE_OPT::MIN:
    lagr += std::get<2>(_obj)[0] * _vec_VAR_FCT_IVP[0];
    scal += std::get<2>(_obj)[0];
    break;
   case BASE_OPT::MAX:
    lagr -= std::get<2>(_obj)[0] * _vec_VAR_FCT_IVP[0];
    scal += std::get<2>(_obj)[0];
    break;
  }
  for( unsigned ip=0; ip<_np; ip++ ){
    if( tvar && tvar[ip] ) continue;
    lagr += ( _pMU[ip] - _pML[ip] ) * _pP[ip];
    scal += _pMU[ip] + _pML[ip];
  }
  for( unsigned ic=0; ic<nc; ++ic ){
    switch( std::get<0>(_ctr)[ic] ){
     case BASE_OPT::EQ:
      lagr += std::get<2>(_ctr)[ic] * _vec_VAR_FCT_IVP[ic+1];
      scal += sqr( std::get<2>(_ctr)[ic] );
      break;
     case BASE_OPT::LE:
      lagr += std::get<2>(_ctr)[ic] * _vec_VAR_FCT_IVP[ic+1];
      scal += std::get<2>(_ctr)[ic];
      break;
     case BASE_OPT::GE:
      lagr -= std::get<2>(_ctr)[ic] * _vec_VAR_FCT_IVP[ic+1];
      scal += std::get<2>(_ctr)[ic];
      break;
    }
  }
  std::get<0>(_nco).push_back( BASE_OPT::EQ );
  std::get<1>(_nco).push_back( scal );
  //std::cout << "scaling:";
  //_pDAG->output( _pDAG->subgraph( 1, &scal ) );

  // Lagrangian stationarity conditions
  std::vector<FFVar> vPF, dPF;
  //auto itsen = _mapFCTSEN.begin();
  for( unsigned ip=0; ip<_np; ip++ ){
    if( tvar && tvar[ip] ) continue;
    vPF.clear(); vPF.push_back( _pP[ip] );
    dPF.clear(); dPF.push_back( 1. );
    unsigned ifct = 0;
    for( auto it=_set_FCT_IVP.begin(); it!=_set_FCT_IVP.end(); ++it, ifct++ ){
      vPF.push_back( _vec_VAR_FCT_IVP[*it] ); 
      dPF.push_back( _vec_VAR_DFCT_IVP[ifct+ip*_set_FCT_IVP.size()] );
    }
    const FFVar* dlagr = _pDAG->FAD( 1, &lagr, vPF.size(), vPF.data(), dPF.data() );
    std::get<0>(_nco).push_back( BASE_OPT::EQ );
    std::get<1>(_nco).push_back( *dlagr );
#ifdef MC__BASE_DO__DEBUG
    std::cout << "dLagr/d" << _pP[ip] << ":";
    _pDAG->output( _pDAG->subgraph( 1, dlagr ) );
#endif
    delete[] dlagr;
  }

  for( unsigned ip=0; ip<_np; ip++ ){
    if( tvar && tvar[ip] ) continue;
    std::get<0>(_nco).push_back( BASE_OPT::EQ );
    std::get<1>(_nco).push_back( ( _pP[ip] - Op<T>::u(P[ip]) ) * _pMU[ip] );
    std::get<0>(_nco).push_back( BASE_OPT::EQ );
    std::get<1>(_nco).push_back( ( _pP[ip] - Op<T>::l(P[ip]) ) * _pML[ip] );
  }

  for( unsigned ic=0; ic<nc; ++ic ){
    switch( std::get<0>(_ctr)[ic] ){
     case BASE_OPT::EQ:
      break;
     case BASE_OPT::LE:
     case BASE_OPT::GE:
      std::get<0>(_nco).push_back( BASE_OPT::EQ );
      std::get<1>(_nco).push_back( _vec_VAR_FCT_IVP[ic+1] * std::get<2>(_ctr)[ic] );
      break;
    }
  }
}

} // end namescape mc

#endif

