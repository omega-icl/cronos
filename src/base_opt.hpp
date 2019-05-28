// Copyright (C) 2015 Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

#ifndef MC__BASE_OPT_HPP
#define MC__BASE_OPT_HPP

#include <vector>
#include <cmath>

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
  double INF;

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
  BASE_OPT(): INF(1e20)
    {}

  //! @brief Class destructor
  virtual ~BASE_OPT()
    {}

protected:
  //! @brief Private methods to block default compiler methods
  BASE_OPT(const BASE_OPT&);
  BASE_OPT& operator=(const BASE_OPT&);
};

//inline double BASE_OPT::INF = 1e20;

//! @brief C++ structure for holding the solution of optimization models
////////////////////////////////////////////////////////////////////////
//! mc::SOLUTION_OPT is a C++ structure for holding the solution of 
//! optimization models, including variables, cost and constraint
//! functions and multiplers.
////////////////////////////////////////////////////////////////////////
struct SOLUTION_OPT
{
  SOLUTION_OPT
    ()
    {}

  ~SOLUTION_OPT
    ()
    {}

  SOLUTION_OPT
    ( const SOLUTION_OPT &sol )
    : status( sol.status ), p( sol.p ), upL( sol.upL ), upU( sol.upU ),
      g( sol.g ), ug( sol.ug ), f( sol.f )
    {}

  void reset()
    { p.clear(); upL.clear(); upU.clear(); g.clear(); ug.clear(); f = NAN; }

  int status;
  std::vector<double> p;
  std::vector<double> upL;
  std::vector<double> upU;
  std::vector<double> g;
  std::vector<double> ug;
  double f;
};

} // end namescape mc

#endif

