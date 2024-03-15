// Copyright (C) Benoit Chachuat, Imperial College London.
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
template <typename... ExtOps>
class BASE_AE
{
public:
  //! @brief Infinity
  static double INF;

protected:
  //! @brief pointer to DAG of equation
  FFGraph<ExtOps...>* _dag;

  //! @brief parameters
  std::vector<FFVar> _par;

  //! @brief decision variables
  std::vector<FFVar> _var;

  //! @brief dependent variables
  std::vector<FFVar> _dep;

  //! @brief equation system
  std::vector<FFVar> _sys;

  //! @brief Whether the system equations have changed
  bool _newsys;

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

  //! @brief Perform block decomposition of system using Harwell's MC13
  bool set_block
    ( const bool disp=false, std::ostream&os=std::cout );

  //! @brief Set block decomposition of system using user arrays IOR, IB and IPERM (per Harwell's MC13 input/output)
  bool set_block
    ( int const NB, int const* IOR=nullptr, int const* IB=nullptr, int const* IPERM=nullptr,
      const bool disp=false, std::ostream&os=std::cout );

  //! @brief Reset block decomposition
  bool reset_block
    ();

public:
  /** @ingroup AEBND
   *  @{
   */

  //! @brief Class constructor
  BASE_AE()
    : _dag(nullptr), _newsys(true), _noblk(0), _singsys(false), _linsys(false)
    {}

  //! @brief Class destructor
  virtual ~BASE_AE()
    {}

  //! @brief Get pointer to DAG
  FFGraph<ExtOps...>* dag()
    const
    { return _dag; }

  //! @brief Set pointer to DAG
  void set_dag
    ( FFGraph<ExtOps...>* dag )
    { _dag = dag; }

  //! @brief Get parameters
  std::vector<FFVar> const& par()
    const
    { return _par; }

  //! @brief Set parameters
  void set_par
    ( std::vector<FFVar> const& par, std::vector<double> const& val=std::vector<double>() )
    { _par = par;
      for( unsigned i=0; i<val.size() && i<_par.size(); i++ ){
        _par[i].set( val[i] );
      }
    }

  //! @brief Add parameters
  void add_par
    ( std::vector<FFVar> const& par, std::vector<double> const& val=std::vector<double>() )
    { _par.insert( _par.end(), par.begin(), par.end() );
      for( unsigned i=0; i<val.size() && i<par.size(); i++ ){
        _par[_par.size()-par.size()+i].set( val[i] );
      }
    }

  //! @brief Set parameters
  void set_par
    ( unsigned const npar, FFVar const* par, double const* val=0 )
    { _par.assign( par, par+npar );
      for( unsigned i=0; val && i<_par.size(); i++ ){
        _par[i].set( val[i] );
      }
    }

  //! @brief Add parameters
  void add_par
    ( unsigned const npar, FFVar const* par, double const* val=0 )
    { _par.insert( _par.end(), par, par+npar );
      for( unsigned i=0; val && i<npar; i++ ){
        _par[_par.size()-npar+i].set( val[i] );
      }
    }

  //! @brief Set parameters
  void set_par
    ( FFVar const& par )
    { _par.assign( &par, &par+1 );
    }

  //! @brief Set parameters
  void set_par
    ( FFVar const& par, double const val )
    { _par.assign( &par, &par+1 );
      _par[0].set( val );
    }

  //! @brief Add parameter
  void add_par
    ( const FFVar&par )
    { _par.push_back( par ); }

  //! @brief Add parameter
  void add_par
    ( const FFVar&par, const double val )
    { _par.push_back( par );
      _par.back().set( val ); }

  //! @brief Reset parameters
  void reset_par
    ()
    { _par.clear(); }

  //! @brief Get decision variables
  std::vector<FFVar> const& var() const
    { return _var; }

  //! @brief Set decision variables
  void set_var
    ( std::vector<FFVar> const& var )
    { _var = var; }

  //! @brief Add decision variables
  void add_var
    ( std::vector<FFVar> const& var )
    { _var.insert( _var.end(), var.begin(), var.end() ); }

  //! @brief Set decision variables
  void set_var
    ( unsigned const nvar, FFVar const* var )
    { _var.assign( var, var+nvar ); }

  //! @brief Add decision variables
  void add_var
    ( unsigned const nvar, FFVar const* var )
    { _var.insert( _var.end(), var, var+nvar ); }

  //! @brief Add decision variable
  void add_var
    ( FFVar const& var )
    { _var.push_back( var ); }

  //! @brief Reset decision variables
  void reset_var
    ()
    { _var.clear(); }

  //! @brief Get dependent variables
  std::vector<FFVar> const& dep
    ()
    const
    { return _dep; }

  //! @brief Get equation system
  std::vector<FFVar> const& sys
    ()
    const
    { return _sys; }

  //! @brief Set dependent variables
  void set_dep
    ( std::vector<FFVar> const& dep, std::vector<FFVar> const& sys )
    { assert( dep.size() == sys.size() );
      _dep = dep;
      _sys = sys;
      _newsys = true; }

  //! @brief Set dependent variables
  void set_dep
    ( unsigned const ndep, FFVar const* dep, FFVar const* sys )
    { _dep.assign( dep, dep+ndep );
      _sys.assign( sys, sys+ndep );
      _newsys = true; }

  //! @brief Add dependent variable
  void add_dep
    ( FFVar const& dep )
    { _dep.push_back( dep ); }
    
  //! @brief Reset dependent variables
  void reset_dep
    ()
    { _dep.clear(); _newsys = true; }

  //! @brief Add algebraic equation
  void add_sys
    ( FFVar const& sys )
    { _sys.push_back( sys );
      _newsys = true; }

  //! @brief Reset algebraic equations
  void reset_sys
    ()
    { _sys.clear(); _newsys = true; }

  //! @brief Copy algebraic system and structure
  void set
    ( BASE_AE<ExtOps...> const& aes )
    { _dag = aes._dag; //std::cout << "DAG: " << aes._dag << std::endl;
      _var = aes._var; _dep = aes._dep; _sys = aes._sys; _newsys = aes._newsys;
      _noblk = aes._noblk; _pblk  = aes._pblk;  _nblk   = aes._nblk;
      _singsys = aes._singsys; _linsys = aes._linsys; _linblk = aes._linblk;
      _lindep = aes._lindep; _bwsys = aes._bwsys; _bwblk = aes._bwblk;
      _fpdep = aes._fpdep; _rpdep = aes._rpdep; }

  //! @brief Number of blocks
  unsigned int noblk
    ()
    const
    { return _noblk; }

  //! @brief Size of block ib
  unsigned int nblk
    ( const unsigned ib )
    const
    { return ib<_noblk? _nblk[ib]: 0; }

  //! @brief Current block
  unsigned int iblk
    ()
    const
    { return _iblk; }

  //! @brief Linearity of block ib
  bool linblk
    ( const unsigned ib )
    const
    { return ib<_noblk? _linblk[ib]: false; }

  //! @brief Equations in block ib
  FFVar const* eqblk
    ( const unsigned ib )
    const
    { return ib<_noblk? _sys.data()+_pblk[ib]: 0; }

  //! @brief Variables in block ib
  FFVar const* depblk
    ( const unsigned ib )
    const
    { return ib<_noblk? _dep.data()+_pblk[ib]: 0; }

  //! @brief Linearity of block ib
  bool lindepblk
    ( const unsigned ib, const unsigned j )
    const
    { return ib<_noblk && j<_nblk[ib]? _lindep[_pblk[ib]+j]: false; }

  //! @brief Forward permutation of dependent variables
  unsigned int pblk
    ( const unsigned ib )
    const
    { return ib<_noblk? _pblk[ib]: 0; }

  //! @brief Forward permutation of dependent variables
  unsigned int fpdep
    ( const unsigned i )
    const
    { return i<_dep.size()? _fpdep[i]: 0; }

  //! @brief Reverse permutation of dependent variables
  unsigned int rpdep
    ( const unsigned i )
    const
    { return i<_dep.size()? _rpdep[i]: 0; }

  //! @brief Reverse permutation of dependent variables
  unsigned int rpdep
    ( const unsigned ib, const unsigned j )
    const
    { return ib<_noblk && j<_nblk[ib]? _rpdep[_pblk[ib]+j]: 0; }
  /** @} */

protected:

  //! @brief Private methods to block default compiler methods
  BASE_AE( BASE_AE<ExtOps...> const& );
  BASE_AE<ExtOps...>& operator=( BASE_AE<ExtOps...> const& );
};

template <typename... ExtOps>
inline double BASE_AE<ExtOps...>::INF = 1E30;

template <typename... ExtOps>
inline
bool
BASE_AE<ExtOps...>::reset_block
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

template <typename... ExtOps>
inline
bool
BASE_AE<ExtOps...>::set_block
( const bool disp, std::ostream&os )
{
  const unsigned int ndep = _dep.size();
  if( !ndep || _sys.size() != ndep ) return false;
  if( !_newsys ) return true;
  //_newsys = false;

  // Perform block lower-triangular decomposition using MC21A/MC13D
  int NB = 1;
  std::vector<int> IPERM(ndep), IOR(ndep), IB(ndep);
  _singsys = !_dag->MC13( ndep, _sys.data(), _dep.data(), IPERM.data(),
                          IOR.data(), IB.data(), NB, disp?true:false, os );
  if( _singsys ) return reset_block();

  return set_block( NB, IOR.data(), IB.data(), IPERM.data() );
}

template <typename... ExtOps>
inline
bool
BASE_AE<ExtOps...>::set_block
( int const NB, int const* IOR, int const* IB, int const* IPERM,
  const bool disp, std::ostream&os )
{
  const unsigned int ndep = _dep.size();
  if( !ndep || _sys.size() != ndep ) return false;
  if( !_newsys ) return true;
  _newsys  = false;
  _singsys = false;

//  // Display permuted system structure
//  if( disp ){
//    std::cout << std::endl << "Number of Blocks: " << NB << std::endl;
//    os << "Lower-triangular block structure:" << std::endl
//       << std::right << "     ";
//    for( unsigned j=0; j<nDep; j++ )
//    //  os << " " << std::setw(3) << IOR[j]-1;
//      os << " " << std::setw(4) << _dep[IOR[j]-1];
//    os << std::endl;
//    for( unsigned i=0; i<nDep; i++ ){
//      //os << std::setw(3) << IPERM[IOR[i]-1]-1 << " ";
//      os << std::setw(4) << _sys[IPERM[IOR[i]-1]-1] << " ";
//      for( unsigned j=0; j<nDep; j++ )
//        os << std::setw(3) << " "
//           << (vDep[IPERM[IOR[i]-1]-1].dep(IOR[j]-1).first?"X ":"  ");
//      os << std::endl;
//    }
//    os << std::endl;
//  }

  // Permute order of equation system AND variables in vectors sys and var,
  // now arranged in upper-triangular block form
  // Keep track of forward and reverse permutations in _fpdep and _rpdep
  std::vector<FFVar> sys(ndep), var(ndep);
  _fpdep.resize(ndep); _rpdep.resize(ndep);
  for( unsigned int i=0; i<ndep; i++ ){
    if( IOR ){
      sys[i] = IPERM? _sys[IPERM[IOR[ndep-i-1]-1]-1]: _sys[IOR[ndep-i-1]-1];
      var[i] = _dep[IOR[ndep-i-1]-1];
      _fpdep[IOR[i]-1] = ndep-i-1;
      _rpdep[ndep-i-1] = IOR[i]-1;
    }
    else{
      sys[i] = IPERM? _sys[IPERM[ndep-i-1]-1]: _sys[ndep-i-1];
      var[i] = _dep[ndep-i-1];
      _fpdep[i] = ndep-i-1;
      _rpdep[ndep-i-1] = i;    
    }
  }
  _sys = sys;
  _dep = var;

  // Keep track of first row and size of each block in permuted matrix in _pblk and _nblk
  _noblk = NB;
  _nblk.resize(_noblk);
  _pblk.resize(_noblk);
  for( int i=0; i<NB; i++ ){
    if( IB ){
      _nblk[i] = ( i==NB-1? ndep+1: IB[i+1] ) - IB[i];
      _pblk[i] = ndep - (IB[i]-1) - _nblk[i];
    }
    else{
      assert( NB == (int)ndep );
      _nblk[i] = 1;
      _pblk[i] = ndep-i-1;
    }
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
  _dag->eval( ndep, sys.data(), depsys.data(), var.size(), var.data(), depvar.data() );
  _bwsys.first = _bwsys.second = 0;
  for( unsigned i=0; i<ndep; i++ ){
    auto cit = depsys[i].dep().begin();
    for( ; cit != depsys[i].dep().end(); ++cit ){
      if( (*cit).second ) // Detecting overall system linearity
          _linsys = false;
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
        if( (*cit).second ) // Detecting block/dependent linearity
          _linblk[ib] = _lindep[_pblk[ib]+(*cit).first] = false;
        if( _bwblk[ib].first  < (int)i-(*cit).first ) // Updating lower band width
          _bwblk[ib].first  = (int)i-(*cit).first;
        if( _bwblk[ib].second < (*cit).first-(int)i ) // updating upper band width
          _bwblk[ib].second = (*cit).first-(int)i;
      }
      //if( !_linblk[ib] ) _linsys = false; // Overall linearity
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

