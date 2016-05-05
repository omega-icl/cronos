// Copyright (C) 2014 Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

#ifndef MC__BASE_NLP_HPP
#define MC__BASE_NLP_HPP

#include "base_opt.hpp"
#include "base_ae.hpp"

namespace mc
{
//! @brief C++ base class for definition of nonlinear programs
////////////////////////////////////////////////////////////////////////
//! mc::BASE_NLP is a C++ base class for definition of the variables,
//! objective and constraints participating in nonlinear programs.
////////////////////////////////////////////////////////////////////////
class BASE_NLP:
  public virtual BASE_OPT,
  public virtual BASE_AE
{
public:
  //! @brief Class constructor
  BASE_NLP()
    : BASE_OPT(), BASE_AE()
    {}

  //! @brief Class destructor
  virtual ~BASE_NLP()
    {}

  //! @brief Get constraints
  const std::tuple< std::vector<t_CTR>, std::vector<FFVar>, std::vector<FFVar> >& ctr() const
    { return _ctr; }

  //! @brief Reset constraints
  void reset_ctr()
    { std::get<0>(_ctr).clear(); std::get<1>(_ctr).clear(); std::get<2>(_ctr).clear(); }

  //! @brief Add constraint
  void add_ctr
    ( const t_CTR type, const FFVar&ctr )
    { std::get<0>(_ctr).push_back( type );
      std::get<1>(_ctr).push_back( ctr );
      std::get<2>(_ctr).push_back( FFVar( ctr.dag() ) ); }

  //! @brief Get objective
  const std::tuple< std::vector<t_OBJ>, std::vector<FFVar>, std::vector<FFVar> >& obj() const
    { return _obj; }

  //! @brief Set objective
  void set_obj
    ( const t_OBJ type, const FFVar&obj )
    { std::get<0>(_obj).clear(); std::get<0>(_obj).push_back( type );
      std::get<1>(_obj).clear(); std::get<1>(_obj).push_back( obj );
      std::get<2>(_obj).clear(); std::get<2>(_obj).push_back( FFVar( obj.dag() ) ); }

  //! @brief Copy equations
  void set
    ( const BASE_NLP&nlp )
    { BASE_AE::set(nlp); _ctr = nlp._ctr; _obj = nlp._obj; }

protected:
  //! @brief constraints (types, constraint variables, constraint multipliers)
  std::tuple< std::vector<t_CTR>, std::vector<FFVar>, std::vector<FFVar> > _ctr;

  //! @brief objective (type, cost variable, cost multiplier)
  std::tuple< std::vector<t_OBJ>, std::vector<FFVar>, std::vector<FFVar> > _obj;

  //! @brief Private methods to block default compiler methods
  BASE_NLP(const BASE_NLP&);
  BASE_NLP& operator=(const BASE_NLP&);
};

} // end namescape mc

#endif

