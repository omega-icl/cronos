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

  //! @brief number of AE block
  unsigned _noblk;

  //! @brief starting position of AE blocks
  std::vector<unsigned> _pblk;

  //! @brief size of AE blocks
  std::vector<unsigned> _nblk;

  //! @brief Current block
  unsigned _iblk;

  //! @brief Structural singularity of system
  bool _singsys;

  //! @brief Linearity of problem w.r.t *all of* the dependents
  bool _linsys;

  //! @brief Linearity of blocks w.r.t. the block variables
  std::vector<bool> _linblk;

  //! @brief Linearity of dependents in blocks
  std::vector<bool> _lindep;

  //! @brief Lower and upper band width of system
  std::pair<long,long> _bwsys;

  //! @brief Lower and upper band width of system blocks
  std::vector< std::pair<long,long> > _bwblk;

  //! @brief variable indices after possible permutation (forward)
  std::vector<unsigned> _fpdep;

  //! @brief variable indices after possible permutation (reverse)
  std::vector<unsigned> _rpdep;

  //! @brief Perform block decomposition of system
  bool set_block
    ( const bool disp=false, std::ostream&os=std::cout );

  //! @brief Reset block decomposition
  bool reset_block
    ();

public:
  /** @ingroup AEBND
   *  @{
   */

  //! @brief Class constructor
  BASE_AE()
    : _dag(0), _noblk(0), _singsys(false), _linsys(false)
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

  //! @brief Reset decision variables
  void reset_var
    ()
    { _var.clear(); }

  //! @brief Add decision variable
  void add_var
    ( const FFVar&var )
    { _var.push_back( var ); }

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
    ( const unsigned ndep, const FFVar*dep, const FFVar*eq )
    { _dep.assign( dep, dep+ndep ); _sys.assign( eq, eq+ndep ); }

  //! @brief Reset dependent variables
  void reset_dep
    ()
    { _dep.clear(); }

  //! @brief Add decision variable
  void add_dep
    ( const FFVar&dep )
    { _dep.push_back( dep ); }

  //! @brief Reset dependent variables
  void reset_sys
    ()
    { _sys.clear(); }

  //! @brief Add decision variable
  void add_sys
    ( const FFVar&eq )
    { _sys.push_back( eq ); }

  //! @brief Copy equations
  void set
    ( const BASE_AE&aes )
    { _dag = aes._dag; _var = aes._var; _dep = aes._dep; _sys = aes._sys; };

  //! @brief Number of blocks
  unsigned int noblk
    () const
    { return _noblk; }

  //! @brief Size of block #ib
  unsigned int nblk
    ( const unsigned ib ) const
    { return ib<_noblk? _nblk[ib]: 0; }

  //! @brief Current block
  unsigned int iblk
    () const
    { return _iblk; }

  //! @brief Linearity of block #ib
  bool linblk
    ( const unsigned ib ) const
    { return ib<_noblk? _linblk[ib]: false; }

  //! @brief Equations in block #ib
  const FFVar* eqblk
    ( const unsigned ib ) const
    { return ib<_noblk? _sys.data()+_pblk[ib]: 0; }

  //! @brief Variables in block #ib
  const FFVar* depblk
    ( const unsigned ib ) const
    { return ib<_noblk? _dep.data()+_pblk[ib]: 0; }

  //! @brief Linearity of block #ib
  bool lindepblk
    ( const unsigned ib, const unsigned j ) const
    { return ib<_noblk && j<_nblk[ib]? _lindep[_pblk[ib]+j]: false; }

  //! @brief Forward permutation of dependent variables
  unsigned int pblk
    ( const unsigned ib ) const
    { return ib<_noblk? _pblk[ib]: 0; }

  //! @brief Forward permutation of dependent variables
  unsigned int fpdep
    ( const unsigned i ) const
    { return i<_dep.size()? _fpdep[i]: 0; }

  //! @brief Reverse permutation of dependent variables
  unsigned int rpdep
    ( const unsigned i ) const
    { return i<_dep.size()? _rpdep[i]: 0; }

  //! @brief Reverse permutation of dependent variables
  unsigned int rpdep
    ( const unsigned ib, const unsigned j ) const
    { return ib<_noblk && j<_nblk[ib]? _rpdep[_pblk[ib]+j]: 0; }

  //! @brief Infinity
  static double INF;
  /** @} */

protected:
  //! @brief Private methods to block default compiler methods
  BASE_AE(const BASE_AE&);
  BASE_AE& operator=(const BASE_AE&);
};

double BASE_AE::INF = 1e20;

inline bool
BASE_AE::reset_block
()
{
  const unsigned int ndep = _dep.size();
  _noblk = 1;
  _nblk.resize(1); _nblk[0] = ndep;
  _pblk.resize(1); _pblk[0] = 0;
  _fpdep.resize(ndep); _rpdep.resize(ndep);
  for( unsigned int i=0; i<ndep; i++ )
    _fpdep[i] = _rpdep[i] = i;
  return true;
}

inline bool
BASE_AE::set_block
( const bool disp, std::ostream&os )
{
  const unsigned int ndep = _dep.size();
  if( !ndep || _sys.size() != ndep ) return false;

  // Perform block lower-triangular decomposition using MC21A/MC13D
  int NB = 1;
  std::vector<int> IPERM(ndep), IOR(ndep), IB(ndep);
  _singsys = !_dag->MC13( ndep, _sys.data(), _dep.data(), IPERM.data(),
    IOR.data(), IB.data(), NB, disp?true:false, os );
  if( _singsys ) return reset_block();

  // Permute order of equation system AND variables in vector _pAE and _pVAR,
  // now arranged in upper-triangular block form
  // Keep track of forward and reverse permutations in _bVAR and _bVARrev
  std::vector<FFVar> sys(ndep), var(ndep);
  _fpdep.resize(ndep); _rpdep.resize(ndep);
  for( unsigned int i=0; i<ndep; i++ ){
    sys[i] = _sys[IPERM[IOR[ndep-i-1]-1]-1];
    var[i] = _dep[IOR[ndep-i-1]-1];
    _fpdep[IOR[i]-1] = ndep-i-1;
    _rpdep[ndep-i-1] = IOR[i]-1;
  }
  _sys = sys;
  _dep = var;

  // Keep track of first row and size of each block in permuted matrix in _pblk and _nblk
  _noblk = NB;
  _nblk.resize(_noblk);
  _pblk.resize(_noblk);
  for( int i=0; i<NB; i++ ){
    _nblk[i] = ( i==NB-1? ndep+1: IB[i+1] ) - IB[i]; 
    _pblk[i] = ndep+1 - IB[i] - _nblk[i];   
  }

  // Systam & block properties (linearity, Jacobian bandwidth)
  var.insert( var.end(), _var.begin(), _var.end() );
  _linblk.resize(_noblk);
  _lindep.resize(ndep);
  _bwblk.resize(_noblk);
  _linsys = true;

  std::vector<FFDep> depsys(ndep), depvar(var.size());
  for( unsigned i=0; i<ndep; i++ ) depvar[i].indep(i);
  //for( unsigned i=0; i<var.size(); i++ ) std::cout << var[i] << ": " << depvar[i] << std::endl;
  _dag->eval( ndep, sys.data(), depsys.data(), var.size(), var.data(),
              depvar.data() );
  _bwsys.first = _bwsys.second = 0;
  for( unsigned i=0; i<ndep; i++ ){
    auto cit = depsys[i].dep().begin();
    for( ; cit != depsys[i].dep().end(); ++cit ){
      if( _bwsys.first  < (int)i-(*cit).first ) // Updating lower band width
        _bwsys.first  = (int)i-(*cit).first;
      if( _bwsys.second < (*cit).first-(int)i ) // updating upper band width
        _bwsys.second = (*cit).first-(int)i;
    }
  }

  for( unsigned ib=0; ib<_noblk; ib++ ){
    std::vector<FFDep> depblk(_nblk[ib]), varblk(var.size()-_pblk[ib]);
    for( unsigned i=0; i<_nblk[ib]; i++ ) varblk[i].indep(i);
    _dag->eval( _nblk[ib], sys.data()+_pblk[ib], depblk.data(),
                var.size()-_pblk[ib], var.data()+_pblk[ib], varblk.data() );
    _linblk[ib] = true;
    for( unsigned i=0; i<_nblk[ib]; i++ ) _lindep[_pblk[ib]+i] = true;
    _bwblk[ib].first = _bwblk[ib].second = 0;

    for( unsigned i=0; i<_nblk[ib]; i++ ){
      auto cit = depblk[i].dep().begin();
      for( ; cit != depblk[i].dep().end(); ++cit ){
        if( !(*cit).second ) // Detecting block/dependent linearity
          _linblk[ib] = _lindep[_pblk[ib]+(*cit).first] = false;
        if( _bwblk[ib].first  < (int)i-(*cit).first ) // Updating lower band width
          _bwblk[ib].first  = (int)i-(*cit).first;
        if( _bwblk[ib].second < (*cit).first-(int)i ) // updating upper band width
          _bwblk[ib].second = (*cit).first-(int)i;
      }
      if( !_linblk[ib] ) _linsys = false; // Overall linearity
    }
#ifdef MC__BASE_AE_DEBUG
    std::cout << "BLOCK #" << ib << ": "
              << (_linblk[ib]?"L":"NL") << ", BW "
              << _bwblk[ib].first << "," << _bwblk[ib].second << std::endl;
#endif
  }

  return true;
}

} // end namescape mc

#endif

