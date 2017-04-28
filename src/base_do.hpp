// Copyright (C) 2014 Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

#ifndef MC__BASE_NLP_HPP
#define MC__BASE_NLP_HPP

#include "base_opt.hpp"
#include "base_de.hpp"

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
  //! @brief Get constraints
  const std::tuple< std::vector<t_CTR>, std::vector< std::map<unsigned,FFVar> >, std::vector<FFVar> >& ctr
    () const
    { return _ctr; }

  //! @brief Reset constraints
  void reset_ctr()
    { std::get<0>(_ctr).clear(); std::get<1>(_ctr).clear(); std::get<2>(_ctr).clear(); }

  //! @brief Add constraint
  void add_ctr
    ( const t_CTR type, const std::map<unsigned,FFVar>&ctrmap )
    { std::get<0>(_ctr).push_back( type );
      std::get<1>(_ctr).push_back( ctrmap );
      std::get<2>(_ctr).push_back( FFVar( ctrmap.begin()->second.dag() ) ); }
  void add_ctr
    ( const t_CTR type, const std::pair<unsigned,FFVar>&ctr )
    { std::get<0>(_ctr).push_back( type );
      std::map<unsigned,FFVar> ctrmap; ctrmap.insert( ctr );
      std::get<1>(_ctr).push_back( ctrmap );
      std::get<2>(_ctr).push_back( FFVar( ctr.second.dag() ) ); }
  void add_ctr
    ( const t_CTR type, const FFVar&ctr )
    { std::get<0>(_ctr).push_back( type );
      std::map<unsigned,FFVar> ctrmap; ctrmap.insert( std::make_pair(0,ctr) );
      std::get<1>(_ctr).push_back( ctrmap );
      std::get<2>(_ctr).push_back( FFVar( ctr.dag() ) ); }

  //! @brief Get objective
  const std::tuple< std::vector<t_OBJ>, std::vector< std::map<unsigned,FFVar> >, std::vector<FFVar> >& obj
    () const
    { return _obj; }

  //! @brief Set objective
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
  //! @brief constraints (type, constraint variable, constraint multiplier)
  std::tuple< std::vector<t_CTR>, std::vector< std::map<unsigned,FFVar> >, std::vector<FFVar> > _ctr;

  //! @brief objective (type, cost variable, cost multiplier)
  std::tuple< std::vector<t_OBJ>, std::vector< std::map<unsigned,FFVar> >, std::vector<FFVar> > _obj;

  //! @brief Define state functions based on cost/constraint functions
  bool set_function
    ();

  //! @brief Test dependency w.r.t. state variables
  bool state_depend
    ( const unsigned ifct );

  //! @brief number of time stages
  unsigned _ns;

  //! @brief vector of functionals
  std::vector<FFVar> _fct;

  //! @brief depency map for state-dependent functions
  t_DepVar _mapFCTSTA;

  //! @brief indices of state-dependent functions
  std::set<unsigned> _ndxFCTSTA;

  //! @brief vector of const pointers to state-dependent functions (assumes additive contributions) in each stage of IVP-DAE system
  std::vector<FFVar> _vFCTSTA;

  //! @brief vector of parameter-dependent functions (no state dependency)
  std::vector<FFVar> _vFCTPAR;

  //! @brief indices of parameters participating in state-dependent functions
  std::set<unsigned> _ndxFCTDEP;

  //! @brief parameters participating in state-dependent functions
  std::vector<FFVar> _vFCTDEP;

  //! @brief indices of parameters participating in parametric IVP
  std::set<unsigned> _ndxSTADEP;

  //! @brief parameters participating in parametric IVP
  std::vector<FFVar> _vSTADEP;

private:
  //! @brief Return number of time stages in cost/constraint functions
  bool _set_ns
    ();

  //! @brief vector holding function contributions
  std::vector<FFVar> _pFCT;

  //! @brief Build-up state function
  void _build_function
    ( const unsigned is, std::map<unsigned,FFVar>::iterator&jt,
      const unsigned ic );

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
  _ns = std::get<1>(_obj)[0].rbegin()->first+1;
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
  for( auto it=_mapFCTSTA.begin(); it!=_mapFCTSTA.end(); ++it )
    if( _fct[ifct].dep().dep( it->second.id().second ).first )
      return true;
  return false;
}

inline bool
BASE_DO::set_depend
()
{
  if( !set_function() || !BASE_DE::set_depend( _ns ) ) return false;
  _ndxFCTSTA.clear();
  for( unsigned ic=0; ic<_nf; ic++ ){
#ifdef MC__BASE_DO__DEBUG
    std::cout << "FCT[" << ic << "]: " << state_depend(ic) << std::endl;
#endif
    if( state_depend( ic ) ) _ndxFCTSTA.insert( ic );
  }

  // Function dependencies
  _vFCTSTA.resize( _ns*_ndxFCTSTA.size() );
  _vFCTPAR.clear();
  _ndxFCTDEP.clear();
  for( unsigned ic=0, icsta=0; ic<_nf; ic++ ){
    auto itc = _ndxFCTSTA.find( ic );
    // State-dependent functions
    if( itc != _ndxFCTSTA.end() ){
      for( unsigned is=0; is<_ns; is++ )
        _vFCTSTA[is*_ndxFCTSTA.size()+icsta] = _vFCT[is][ic];
      icsta++;
      // Append parametric dependencies
      for( unsigned ip=0; ip<_np; ip++ ){
        bool dep = _depF[ic].dep( _pP[ip].id().second ).first;
        if( dep ) _ndxFCTDEP.insert( ip );
      }
    }
    // Parameter-dependent functions
    else{
      _vFCTPAR.push_back( _fct[ic] );
    }
  }

  // Parameters participating in state-dependent functions
  _vFCTDEP.clear();
  for( auto itp=_ndxFCTDEP.begin(); itp!=_ndxFCTDEP.end(); ++itp )
    _vFCTDEP.push_back( _pP[*itp] );

  // Parametric dependencies in IVP
  _ndxSTADEP.clear();
  for( auto it=_mapFCTSTA.begin(); it!=_mapFCTSTA.end(); ++it )
    for( unsigned ip=0; ip<_np; ip++ ){
      bool dep = it->first.second < _nx?
                 _depX[_nx*it->first.first+it->first.second].dep( _pP[ip].id().second ).first:
                 _depQ[_nq*it->first.first+it->first.second-_nx].dep( _pP[ip].id().second ).first;
      if( dep ) _ndxSTADEP.insert( ip );
    }

  // Parameters participating in IVP
  _vSTADEP.clear();
  for( auto itp=_ndxSTADEP.begin(); itp!=_ndxSTADEP.end(); ++itp )
    _vSTADEP.push_back( _pP[*itp] );

  return true;
}

inline void
BASE_DO::_build_function
( const unsigned is, std::map<unsigned,FFVar>::iterator&jt,
  const unsigned ic )
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
    if( _mapFCTSTA.find( ndxq ) == _mapFCTSTA.end() )
      _mapFCTSTA[ndxq] = FFVar( _pDAG );
    termsta.push_back( pXQ );
    termdep.push_back( _mapFCTSTA[ndxq] );
  }
  // Append objective term to _fct[0]
  const mc::FFVar* term = _pDAG->compose( 1, &jt->second, termdep.size(),
                                          termsta.data(), termdep.data() );
  _fct[ic] += term[0];
  if( termdep.size() ) delete[] term;
}

inline bool
BASE_DO::set_function
()
{
  // Resize function vector
  if( !_set_ns() ) return false;
  _nf = 1 + std::get<0>(_ctr).size();
  _pFCT.resize( _ns*_nf );
  _fct.resize( _nf, FFVar(0) );
  _mapFCTSTA.clear();

  // Append cost contributions
  unsigned is = 0;
  auto jt = std::get<1>(_obj)[0].begin();
  for( ; jt!=std::get<1>(_obj)[0].end(); ++jt ){
    for( ; is<jt->first; ++is ) _pFCT[_nf*is] = 0;
    _pFCT[_nf*jt->first] = jt->second;
    _build_function( is, jt, 0 );
    ++is;
  }
  for( ; is<_ns; ++is ) _pFCT[_nf*is] = 0;
#ifdef MC__BASE_DO__DEBUG
  std::cout << "_fct[0] <-- " << _fct[0].dep() << std::endl;
  //_pDAG->output( _pDAG->subgraph( 1, &_fct[0] ) );
  //{ int dum; std::cout << "PAUSED "; std::cin >> dum; }
#endif

  // Append constraint contributions
  auto it = std::get<1>(_ctr).begin();
  for( unsigned ic=1; it!=std::get<1>(_ctr).end(); ++it, ++ic ){
    unsigned is = 0;
    auto jt = it->begin();
    for( ; jt!=it->end(); ++jt ){
      for( ; is<jt->first; ++is ) _pFCT[_nf*is+ic] = 0;
      _pFCT[_nf*jt->first+ic] = jt->second;
      _build_function( is, jt, ic );
      ++is;
    }
    for( ; is<_ns; ++is ) _pFCT[_nf*is+ic] = 0;
#ifdef MC__BASE_DO__DEBUG
    std::cout << "_fct[" << ic << "] <-- " << _fct[ic].dep() << std::endl;
    //_pDAG->output( _pDAG->subgraph( 1, &_fct[ic] ) );
    //{ int dum; std::cout << "PAUSED "; std::cin >> dum; }
#endif
  }

#ifdef MC__BASE_DO__DEBUG
  for( auto itv=_mapFCTSTA.begin(); itv!=_mapFCTSTA.end(); ++itv )
    std::cout << "_mapFCTSTA[" << itv->first.first << "," << itv->first.second << "] = "
              << itv->second << std::endl;
  for( unsigned is=0; is<_ns; is++ )
    for( unsigned ic=0; ic<_nf; ic++ )
      std::cout << "FCT[" << is << "][" << ic << "] = " << _pFCT[_nf*is+ic] << std::endl;
  int dum;  std::cin >> dum; 
#endif

  // vector of pointers
  _vFCT.clear();
  for( unsigned i=0; i<_ns; i++ ) _vFCT.push_back( _pFCT.data()+i*_nf );

  return true;
}

} // end namescape mc

#endif

