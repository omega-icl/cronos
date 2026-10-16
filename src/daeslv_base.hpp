// Copyright (C) Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the EPL-2.0 with GPL-2.0-or-later as a Secondary License; see the LICENSE file.
// SPDX-License-Identifier: EPL-2.0 OR GPL-2.0-or-later

#ifndef CRONOS__DAESLV_BASE_HPP
#define CRONOS__DAESLV_BASE_HPP

#include <string>
#include <vector>

#include "base_dynamic.hpp"

namespace mc
{
//! @brief C++ base class for the DAE form of a model, as an IDAS integrator reads it
////////////////////////////////////////////////////////////////////////
//! mc::DAESLV_BASE is the DAE counterpart of mc::ODESLV_BASE, and is deliberately THIN: the model
//! description, the extraction and its shared machinery all live in mc::BASE_DYNAMIC.  What is here is
//! only what the two consumers genuinely disagree about.
//!
//! ODESLV reads the rows PER STATE -- row i supplies d/dt of state i, divided through so the integrator
//! gets x' = f.  IDAS reads them PER ROW: the residual is the stacked equation list, and nothing has to
//! belong to a particular state.  So this class takes BASE_DYNAMIC's general mass-matrix form,
//!
//!     sum over _mMASS[s][k] of ( coefficient * d _mX[index] / dt )  +  _mRES[s][k]  =  0
//!
//! and leaves _mLHS / _mDE -- ODESLV's diagonal projection -- untouched and empty.  That is what lets a
//! COUPLED mass matrix through: it has no row that supplies one state's derivative on its own, which is
//! exactly why ODESLV refuses it and points at IDAS.
////////////////////////////////////////////////////////////////////////
class DAESLV_BASE
: public virtual BASE_DYNAMIC
{
public:

  //! The revision of the DAESLV headers, so a run can say which it used (mirrors ODESLV_BASE::HEADER_ID,
  //! which ODESLV_revision prints).
  static constexpr char const* HEADER_ID
    = "daeslv  rev1  2026-10-16";

  //! @brief The revision of the DAESLV headers (HEADER_ID), e.g. for a bug report
  static char const* revision() { return HEADER_ID; }

  DAESLV_BASE()
    {}
  virtual ~DAESLV_BASE()
    {}

  //! @brief The id vector IDASetId() wants: 1 for a state whose derivative appears in some row
  //! (DIFFERENTIAL), 0 for one whose derivative appears nowhere (ALGEBRAIC).  Empty before extraction.
  std::vector<double> const& id
    ()
    const
    { return _mID; }

  //! @brief Number of differential states (entries of id() equal to 1).
  size_t n_differential
    ()
    const
    { size_t n = 0; for( double const v : _mID ) if( v != 0. ) ++n; return n; }

  //! @brief Number of algebraic states (entries of id() equal to 0).
  size_t n_algebraic
    ()
    const
    { return _mID.size() - n_differential(); }

  //! @brief Number of evolution rows of the residual, i.e. the size of the system IDAS integrates.
  size_t n_residual
    ()
    const
    { return _mMASS.empty()? 0: _mMASS.front().size(); }

protected:

  //! @brief 1 (differential) / 0 (algebraic) per state, for IDASetId().
  std::vector<double> _mID;

  //! @brief Collect the general mass-matrix form -- that IS the DAE residual.
  bool _wants_mass_rows
    ()
    const
    override
    { return true; }

  //! @brief An algebraic state has no row giving its derivative, by definition, so the explicit-row
  //! requirement that ODESLV imposes would refuse every DAE.
  bool _require_explicit_rows
    ()
    const
    override
    { return false; }

  //! @brief A DAE consumer is stopped by FEWER forms than an explicit integrator, and decides from the
  //! STRUCTURED fields rather than the blocker message.  It accepts the two forms that exist precisely
  //! for it -- a coupled mass matrix (decoupled == false) and algebraic states -- and keeps the rest.
  //! NOTE the last test is deliberately conservative: when the form is decoupled AND has no algebraic
  //! states yet still carries a blocker, that blocker is about something else and is honoured.
  bool _blocked_by_dynamic_form
    ()
    const
    override
    { auto const& F = dynamic_form();
      if( !F.resolved )      return true;   // a high index the reduction could not resolve
      if( !F.initial_value ) return true;   // something imposed at the FAR END: not an initial-value problem
      if( F.blocker.empty() ) return false;
      return F.decoupled && !F.algebraic; }

  //! @brief The per-state projection ODESLV uses.  Deliberately empty here: a DAE consumer reads the rows
  //! through _mMASS / _mRES, so returning false leaves _mLHS / _mDE untouched and lets no row claim a state.
  bool _read_evolution_row
    ( FFVar const&, FFVar const&, size_t, std::vector<FFVar>&, std::vector<FFVar>& )
    override
    { return false; }

  //! @brief Extract the model as a DAE: the shared extraction, then the id vector and the squareness check.
  bool _extract_dae
    ();

  //! @brief Fill _mID from the extracted rows: a state is DIFFERENTIAL exactly when its derivative appears
  //! in some row.  Cross-checked against dynamic_form().differential, which reaches the same answer by a
  //! different route -- a mismatch means the extraction and the classification disagree, and that is a bug
  //! worth failing on rather than integrating through.
  bool _build_id
    ();
};

inline
bool
DAESLV_BASE::_build_id
()
{
  _mID.assign( _nx, 0. );
  for( auto const& stage : _mMASS )
    for( auto const& row : stage )
      for( auto const& [ ix, co ] : row ){
        if( ix >= _nx ){ _extractError = "mass-matrix row names a state out of range"; return false; }
        _mID[ ix ] = 1.;
      }
  size_t const ndiff = n_differential();
  if( ndiff != dynamic_form().differential ){
    _extractError = "the extracted rows give " + std::to_string( ndiff )
                  + " differential state(s) but the classification says "
                  + std::to_string( dynamic_form().differential );
    return false;
  }
  return true;
}

inline
bool
DAESLV_BASE::_extract_dae
()
{
  if( !_extract_from_model() ) return false;
  if( _mMASS.empty() ){ _extractError = "no evolution rows were collected"; return false; }
  if( !_build_id() ) return false;
  // IDAS integrates a SQUARE residual: one row per unknown.  The unknowns are the states, so this is the
  // place a mis-posed model shows up, before any IDACalcIC is asked to make sense of it.
  if( n_residual() != _nx ){
    _extractError = std::to_string( n_residual() ) + " evolution row(s) for " + std::to_string( _nx )
                  + " state(s): the residual IDAS integrates must be square";
    return false;
  }
  return true;
}

} // end namespace mc

#endif
