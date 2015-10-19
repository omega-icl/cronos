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
class BASE_NLP: public virtual BASE_OPT, public virtual BASE_AE
{
public:
  //! @brief Class constructor
  BASE_NLP()
    : BASE_OPT(), BASE_AE()
    {}

  //! @brief Class destructor
  virtual ~BASE_NLP()
    {}

  //! @brief Copy equations
  void set
    ( const BASE_NLP&nlp )
    { BASE_OPT::set(nlp); BASE_AE::set(nlp); }

protected:
  //! @brief Private methods to block default compiler methods
  BASE_NLP(const BASE_NLP&);
  BASE_NLP& operator=(const BASE_NLP&);
};

//double BASE_NLP::INF = 1e20;

} // end namescape mc

#endif

