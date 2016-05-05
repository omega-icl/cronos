// Copyright (C) 2014 Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

#ifndef MC__BASE_NLP_HPP
#define MC__BASE_NLP_HPP

#include "base_opt.hpp"
#include "base_de.hpp"

namespace mc
{
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
  //! @brief Class constructor
  BASE_DO()
    : BASE_OPT(), BASE_DE()
    {}

  //! @brief Class destructor
  virtual ~BASE_DO()
    {}
  //! @brief Get constraints
  const std::tuple< std::vector<t_CTR>, std::vector< std::map<unsigned,FFVar> >, std::vector<FFVar> >& ctr() const
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
  const std::tuple< std::vector<t_OBJ>, std::vector< std::map<unsigned,FFVar> >, std::vector<FFVar> >& obj() const
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

protected:
  using BASE_DE::set_function;

  //! @brief constraints (types, constraint variables, constraint multipliers)
  std::tuple< std::vector<t_CTR>, std::vector< std::map<unsigned,FFVar> >, std::vector<FFVar> > _ctr;

  //! @brief objective (type, cost variable, cost multiplier)
  std::tuple< std::vector<t_OBJ>, std::vector< std::map<unsigned,FFVar> >, std::vector<FFVar> > _obj;

  //! @brief Private methods to block default compiler methods
  BASE_DO(const BASE_DO&);
  BASE_DO& operator=(const BASE_DO&);
};

} // end namescape mc

#endif

