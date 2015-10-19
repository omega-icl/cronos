// Copyright (C) 2015 Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

#ifndef MC__BASE_OPT_HPP
#define MC__BASE_OPT_HPP

#include "ffunc.hpp"

namespace mc
{
//! @brief C++ base class for definition of optimization problems
////////////////////////////////////////////////////////////////////////
//! mc::BASE_OPT is a C++ base class for definition of the objective 
//! and constraint functions participating in optimization problems.
////////////////////////////////////////////////////////////////////////
class BASE_OPT
{
public:
  //! @brief Infinity
  static double INF;
  //! @brief Enumeration type for objective function
  enum t_OBJ{
    MIN=0,	//!< Minimization
    MAX		//!< Maximization
  };
  //! @brief Enumeration type for constraints
  enum t_CTR{
    EQ=0,	//!< Equality constraint
    LE,		//!< Inequality constraint
    GE		//!< Inequality constraint
  };

protected:
  //! @brief constraints (types, constraint variables, constraint multipliers)
  std::tuple< std::vector<t_CTR>, std::vector<FFVar>, std::vector<FFVar> > _ctr;

  //! @brief objective (type, cost variable, cost multiplier)
  std::tuple< std::vector<t_OBJ>, std::vector<FFVar>, std::vector<FFVar> > _obj;

public:
  //! @brief Class constructor
  BASE_OPT()
    {}

  //! @brief Class destructor
  virtual ~BASE_OPT()
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
    ( const BASE_OPT&op )
    { _ctr = op._ctr; _obj = op._obj; }

protected:
  //! @brief Private methods to block default compiler methods
  BASE_OPT(const BASE_OPT&);
  BASE_OPT& operator=(const BASE_OPT&);
};

double BASE_OPT::INF = 1e20;

} // end namescape mc

#endif

