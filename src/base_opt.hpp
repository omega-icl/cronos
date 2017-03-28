// Copyright (C) 2015 Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

#ifndef MC__BASE_OPT_HPP
#define MC__BASE_OPT_HPP

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

  //! @brief Class constructor
  BASE_OPT()
    {}

  //! @brief Class destructor
  virtual ~BASE_OPT()
    {}

protected:
  //! @brief Private methods to block default compiler methods
  BASE_OPT(const BASE_OPT&);
  BASE_OPT& operator=(const BASE_OPT&);
};

double BASE_OPT::INF = 1e20;

} // end namescape mc

#endif

