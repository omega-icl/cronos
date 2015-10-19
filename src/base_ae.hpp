// Copyright (C) 2014 Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

#ifndef MC__BASE_AE_HPP
#define MC__BASE_AE_HPP

#include <assert.h>
#include "ffunc.hpp"

namespace mc
{
//! @brief C++ base class for defining of parametric nonlinear equations
////////////////////////////////////////////////////////////////////////
//! mc::BASE_AE is a C++ base class for defining parametric algebraic
//! equations, distinguishing between independent and dependent
//! variables
////////////////////////////////////////////////////////////////////////
class BASE_AE
{
protected:
  //! @brief pointer to DAG of equation
  FFGraph* _dag;

  //! @brief decision variables
  std::vector<FFVar> _var;

  //! @brief dependent variables
  std::vector<FFVar> _dep;

  //! @brief equation system
  std::vector<FFVar> _sys;

public:
  /** @ingroup AEBND
   *  @{
   */

  //! @brief Class constructor
  BASE_AE()
    : _dag(0)
    {}

  //! @brief Class destructor
  virtual ~BASE_AE()
    {}

  //! @brief Get pointer to DAG
  FFGraph* dag() const
    { return _dag; }

  //! @brief Set pointer to DAG
  void set_dag
    ( FFGraph*pDAG )
    { _dag = pDAG; }

  //! @brief Get decision variables
  const std::vector<FFVar>& var() const
    { return _var; }

  //! @brief Set decision variables
  void set_var
    ( const std::vector<FFVar>&var )
    { _var = var; }

  //! @brief Set decision variables
  void set_var
    ( const unsigned nvar, const FFVar*var )
    { _var.assign( var, var+nvar ); }

  //! @brief Get dependent variables
  const std::vector<FFVar>& dep() const
    { return _dep; }

  //! @brief Get equation system
  const std::vector<FFVar>& sys() const
    { return _sys; }

  //! @brief Set dependent variables
  void set_dep
    ( const std::vector<FFVar>&dep, const std::vector<FFVar>&sys )
    { assert( dep.size()==sys.size() ); _dep = dep; _sys = sys; }

  //! @brief Set dependent variables
  void set_dep
    ( const unsigned ndep, const FFVar*dep, const FFVar*sys )
    { _dep.assign( dep, dep+ndep ); _sys.assign( sys, sys+ndep ); }

  //! @brief Copy equations
  void set
    ( const BASE_AE&aes )
    { _dag = aes._dag; _var = aes._var; _dep = aes._dep; _sys = aes._sys; };
  /** @} */

protected:
  //! @brief Private methods to block default compiler methods
  BASE_AE(const BASE_AE&);
  BASE_AE& operator=(const BASE_AE&);
};

} // end namescape mc

#endif

