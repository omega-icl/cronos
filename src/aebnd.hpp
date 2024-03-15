// Copyright (C) Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

#ifndef MC__AEBND_HPP
#define MC__AEBND_HPP

//#undef  MC__AEBND_DEBUG

#include <stdexcept>
#include <cassert>
#include <cmath>
#include <fstream>
#include <vector>

#include "base_ae.hpp"
#include "mctime.hpp"
#include "tmodel.hpp"
#include "ffunc.hpp"

namespace mc
{
//! @brief C++ class computing enclosures of the solution set of parametric AEs.
////////////////////////////////////////////////////////////////////////
//! mc::AEBND is a C++ class that computes enclosures of the solution
//! set of parametric algebraic equations (AEs). It implements the 
//! methods of Gaussian elimination and both regular and simplified
//! Gauss-Siedel and Krawczyk methods depending on linearity in the
//! dependent variables. Besides interval bounds, polynomial models
//! with interval or ellipsoidal remainders are used to enable higher-
//! order convergence. 
////////////////////////////////////////////////////////////////////////
template < typename T,
           typename PMT = TModel<T>,
           typename PVT = TVar<T>,
           typename... ExtOps >
class AEBND
: public virtual BASE_AE<ExtOps...>
{
public:

  using BASE_AE<ExtOps...>::set;
  using BASE_AE<ExtOps...>::dag;
  using BASE_AE<ExtOps...>::set_dag;
  using BASE_AE<ExtOps...>::par;
  using BASE_AE<ExtOps...>::set_par;
  using BASE_AE<ExtOps...>::add_par;
  using BASE_AE<ExtOps...>::reset_par;
  using BASE_AE<ExtOps...>::var;
  using BASE_AE<ExtOps...>::set_var;
  using BASE_AE<ExtOps...>::add_var;
  using BASE_AE<ExtOps...>::reset_var;
  using BASE_AE<ExtOps...>::dep;
  using BASE_AE<ExtOps...>::set_dep;
  using BASE_AE<ExtOps...>::add_dep;
  using BASE_AE<ExtOps...>::reset_dep;
  using BASE_AE<ExtOps...>::sys;
  using BASE_AE<ExtOps...>::add_sys;
  using BASE_AE<ExtOps...>::reset_sys;

private:

  using BASE_AE<ExtOps...>::_dag;
  using BASE_AE<ExtOps...>::_sys;
  using BASE_AE<ExtOps...>::_iblk;
  using BASE_AE<ExtOps...>::INF;

  typedef std::tuple< unsigned, unsigned const*, unsigned const*, FFVar const*> FFSDer;

  //! @brief number of all (dependent and independent) variables
  unsigned                      _nvar;

  //! @brief number of dependent variables
  unsigned                      _ndep;

  //! @brief subgraph of each AE block
  std::vector< FFSubgraph >     _sgsys;

  //! @brief subgraph of each AE block Jacobian
  std::vector< FFSubgraph >     _sgjac;

  //! @brief maximum number of operations among any AE block evaluation
  unsigned                      _maxop;

  //! @brief intermediate operations during AE evaluation in T arithmetic
  std::vector<T>                _Iwk;

  //! @brief intermediate operations during AE evaluation in PM arithmetic
  std::vector<PVT>              _PMwk;

  //! @brief number of AE block
  unsigned                      _noblk;

  //! @brief starting position of AE blocks
  std::vector<unsigned>         _pblk;

  //! @brief size of AE blocks
  std::vector<unsigned>         _nblk;

  //! @brief maximum block size among all AE blocks
  unsigned                      _nblkmax;

  //! @brief leading dimension of AE blocks (accounts for block recursivity)
  std::vector<unsigned>         _ldblk;

  //! @brief maximum leading dimension among all AE blocks
  unsigned                      _ldblkmax;

  //! @brief variable indices after possible permutation (forward)
  std::vector<unsigned>         _fpdep;

  //! @brief variable indices after possible permutation (reverse)
  std::vector<unsigned>         _rpdep;

  //! @brief AE system Jacobian entries
  std::vector<FFSDer>           _jac;

  //! @brief dependent and independent variables for DAG evaluation
  std::vector<FFVar>            _var;

  //! @brief variable bounds for DAG evaluation
  std::vector<T>                _Ivar;

  //! @brief reference variable bounds for DAG evaluation
  std::vector<T>                _Iref;

  //! @brief pointer to dependent interval bounds **DO NOT FREE**
  T*                            _Ix;

  //! @brief pointer to parameter interval bounds **DO NOT FREE**
  T*                            _Ip;

  //! @brief function interval bounds
  std::vector<T>                _If;

  //! @brief function Jacobian interval bounds
  std::vector<std::vector<T>>   _Idfdx;

  //! @brief variable bounds for DAG evaluation
  std::vector<PVT>              _PMvar;

  //! @brief variable bounds/values for DAG evaluation
  std::vector<PVT>              _PMref;

  //! @brief pointer to state polynomial models **DO NOT FREE**
  PVT*                          _PMx;

  //! @brief pointer to parameter polynomial models **DO NOT FREE**
  PVT*                          _PMp;

  //! @brief function polynomial models
  std::vector<PVT>              _PMf;

  //! @brief function Jacobian interval bounds
  std::vector<std::vector<PVT>> _PMdfdx;

  //! @brief Preconditioning matrix
  CPPL::dgematrix               _Y;

public:
  /** @defgroup AEBND Set-valued techniques for solution of parametric AEs
   *  @{
   */
  //! @brief Default constructor
  AEBND();

  //! @brief Virtual destructor
  virtual ~AEBND();

  //! @brief Solver status
  enum STATUS{
     NORMAL=0,	//!< Normal execution
     EMPTY,     //!< Empty solution set
     SINGULAR,	//!< Structural or numerical singularity encountered
     FAILURE 	//!< Error encountered in the underlying arithmetic
  };

  //! @brief Integrator options
  struct Options
  {
    //! @brief Constructor
    Options( int iDISP=1 ):
      BOUNDER(ALGORITHM::AUTO), PRECOND(PRECONDITIONING::INVMD),
      BLKDEC(DECOMPOSITION::RECUR), AUTODIFF(DIFFERENTIATION::FSA),
      MAXREC(0), INTERBND(true), MAXIT(10), RTOL(1e-7), ATOL(machprec()),
      PMNOREM(false), DISPLEVEL(iDISP)
      {}
    //! @brief Enumeration of bounding methods (for polynomial models only)
    enum class ALGORITHM{
      AUTO=0,  //!< Automatic selection
      KRAW=1,  //!< Krawczyk method
      KRAWS=2, //!< Simplified Krawczyk method
      GS=3,    //!< Gauss-Seidel
      GSS=4,   //!< Simplified Gauss-Seidel
      GE=5     //!< Gauss Elimination
    };
    //! @brief Bounding method
    ALGORITHM BOUNDER;
    //! @brief Enumeration of preconditioning methods
    enum class PRECONDITIONING{
      NONE=0,   //!< No preconditioning used
      INVMD=1,  //!< Inverse of Jacobian mid-point (dense)
      INVMB=2,  //!< Inverse of Jacobian mid-point (banded)
      QRMD=3    //!< QR decomposition of Jacobian mid-point (dense)
          };
    //! @brief Bounding method
    PRECONDITIONING PRECOND;
    //! @brief Enumeration of preconditioning methods
    enum class DECOMPOSITION{
      NONE=0,   //!< No decomposition used
      DIAG=1,   //!< Full decomposition with independent diagonal blocks
      RECUR=2   //!< Full decomposition with recursive triangular blocks
    };
    //! @brief Level of block decomposition (default, recursive: 2, diagonal: 1, no decomposition: 0)
    DECOMPOSITION BLKDEC;
    //! @brief Differentiation strategy
    enum DIFFERENTIATION{
      FSA=0,      //!< Forward sensitivity analysis
      ASA         //!< Adjoint sensitivity analysis
    };
    //! @brief Differentiation strategy for Jacobian matrix
    DIFFERENTIATION AUTODIFF;
    //! @brief Maximum number of recursive blocks (default: 0 - no limit)
    unsigned int MAXREC;
    //! @brief Whether or not to intersect bounds with previous iteration in iterative methods (default: true)
    bool INTERBND;
    //! @brief Maximum number of iterations (default: 10, no limit: 0)
    unsigned int MAXIT;
    //! @brief Relative stopping tolerance (default: 1e-7)
    double RTOL;
    //! @brief Absolute stopping tolerance (default: MACHPREC)
    double ATOL;
    //! @brief Whether or not to cancel the remainder term in polynomial models in order to determine a polynomial approximant only (Default: false)
    bool PMNOREM;
    //! @brief Display level (default: 1)
    int DISPLEVEL;
  } options;

  //! @brief Structure for AE bounding exceptions
  class Exceptions
  {
  public:
    //! @brief Enumeration type for AEBND exception handling
    enum TYPE{
      APRIORI=1,	//!< A priori enclosure required for selected method
      PRECOND,		//!< System preconditionning failed, e.g. due to a singular matrix
      GAUSSEL,		//!< Gauss elimination may not be applied to nonlinear implicit systems
      INTERN=-33	//!< Internal error
    };
    //! @brief Constructor for error <a>ierr</a>
    Exceptions( TYPE ierr ) : _ierr( ierr ){}
    //! @brief Inline function returning the error flag
    int ierr(){ return _ierr; }
    //! @brief Inline function returning the error description
    std::string what(){
      switch( _ierr ){
      case APRIORI:
        return "AEBND::Exceptions  A priori enclosure required for selected method";
      case PRECOND:
        return "AEBND::Exceptions  System preconditionning failed, e.g. due to a singular matrix";
      case GAUSSEL:
        return "AEBND::Exceptions  Gauss elimination  may not be applied to nonlinear equation systems";
      case INTERN: default:
        return "AEBND::Exceptions  Internal error";
      }
    }
  private:
    TYPE _ierr;
  };

  //! @brief Setup block decomposition of parametric AEs, subgraphs, and required differentiation for each block
  STATUS setup
    ( bool const updtdec = true, std::ostream& os = std::cout );

  //! @brief Setup block decomposition of parametric AEs, subgraphs, and required differentiation for each block using user arrays IOR, IB and IPERM (per Harwell's MC13 input/output)
  STATUS setup
    ( int const NB, int const* IOR, int const* IB, int const* IPERM = nullptr,
      std::ostream& os = std::cout );

  //! @brief Compute interval enclosure of solution set of parametric AEs
  STATUS solve
    ( T const* Ip, T* Ix, T const* Ix0 = nullptr, std::ostream& os = std::cout );

  //! @brief Compute polynomial model of solution set of parametric AEs
  template <typename INI>
  STATUS solve
    ( PVT const* PMp, PVT* PMx, INI const* PMx0 = nullptr, std::ostream& os = std::cout );

  //! @brief Compute symbolic solution (DAG) of parametric AEs
  //STATUS solve
  //  ( FFVar* X, std::ostream& os = std::cout );

  //! @brief Whether the blocks contain a unique solution branch
  const std::vector<bool>& uniblk
    ()
    const
    { return _uniblk; }

  //! @brief Whether the blocks contain a unique solution branch
  bool uniblk
    ( unsigned const ib )
    const
    { return _uniblk.at( ib ); }
  /** @} */

private:

  //! @brief Flag for setup function
  bool _issetup;

  //! @brief Lower-bound crossing by dependent branch
  std::vector<int> _xlodep;

  //! @brief Upper-bound crossing by dependent branch
  std::vector<int> _xupdep;

  //! @brief Solution uniqueness in AE blocks
  std::vector<bool> _uniblk;

  //! @brief Linearity of problem blocks w.r.t. the block variables
  std::vector<bool> _linblk;

  //! @brief Lower and upper band width of problem blocks
  std::vector< std::pair<long,long> > _bwblk;


  //! @brief Setup block decomposition of parametric AEs, subgraphs, and required differentiation for each block
  STATUS _setup
    ( std::ostream& os );

  //! @brief Compute preconditionned equation system LHS and RHS for given block
  template <typename U, typename V>
  void _precondlin
    ( unsigned const ib, std::vector<U>& A, std::vector<U>& b,
      std::vector<U>& wkf, U* f, U* ref, std::vector<V>& wkjacf,
      std::vector<V>* jacf, V* jacvar, std::ostream& os );

  //! @brief Initialize bounding
  bool _init
    ();

  //! @brief Initialize bounding for interval case
  bool _init
    ( const unsigned np, T const* Ip );

  //! @brief Initialize bounding for polynomial model case
  bool _init
    ( const unsigned np, PVT const* PMp );

  //! @brief Set reference value (mid-point) for interval bounds
  void _reference
    ( unsigned const ib, T* var,  T* ref, T* jacvar, bool zeroref=false );

  //! @brief Set reference model (centered polynomial) for polynomial models
  void _reference
    ( unsigned const ib, PVT* var, PVT* ref, PVT* jacvar, bool zeroref=false );

  //! @brief Set reference model (centered polynomial) for polynomial models
  void _reference
    ( unsigned const ib, PVT* var, PVT* ref, T* jacvar, bool zeroref=false );

  //! @brief Test convergence for interval bounds
  bool _cvgtest
    ( T const& x )
    const;

  //! @brief Test convergence for polynomial models
  bool _cvgtest
    ( PVT const& x )
    const;

  //! @brief Test convergence for interval bounds
  bool _cvgtest
    ( unsigned const nx, T const* x, T const* x0 )
    const;

  //! @brief Test convergence for polynomial models
  bool _cvgtest
    ( unsigned const nx, PVT const* x, PVT const* x0 )
    const;

  //! @brief Test for a slution branch crossing interval bounds - return value is whether a solution branch may cross any dependent lower or upper bound (true) or not (false)
  bool _crosstest
    ( unsigned const nx, T const* x, T const* x0, int* xlodep, int* xupdep )
    const;

  //! @brief Test for a solution branch crossing polynomial models - return value is whether a solution branch may cross any dependent lower or upper bound (true) or not (false)
  bool _crosstest
    ( unsigned const nx, PVT const* x, PVT const* x0, int* xlodep, int* xupdep )
    const;

  //! @brief Cancel remainder term in polynomial models
  T _cancelrem
    ( T& x )
    const;

  //! @brief Cancel remainder term in polynomial models
  PVT _cancelrem
    ( PVT& x )
    const;

  //! @brief Apply Gauss-Seidel method to given block
  template <typename U, typename V>
  STATUS _gs
    ( unsigned const ib, U* var, std::vector<U>& wkf, U* f, U* ref,
      std::vector<V>& wkjacf, std::vector<V>* jacf, V* jacvar, bool usekraw,
      std::ostream& os );

  //! @brief Apply Gauss elimination method to given block (linear system only)
  template <typename U>
  STATUS _ge
    ( unsigned const ib, U* var, std::vector<U>& wk, U* f, U* ref,
      std::vector<U>* jacf, std::ostream& os );

  //! @brief Apply Gauss elimination method for symbolic solution (linear system only)
  STATUS _ge
    ( std::vector<FFVar>&X, std::ostream&os );

  //! @brief Index position in 2d-array
  static unsigned _ndx
    ( const unsigned i, const unsigned j, const unsigned n )
    { return i*n+j; } // Jacobian in DAG is transpose

  //! @brief Structure storing implicit solver statistics
  struct Stats
  {
    //! @brief Constructor
    Stats():
      cputime(0.), precond(0.), evalrhs(0.), evaljac(0.), maxIter(0)
      {}
    //! @brief Constructor
    void reset()
      { cputime = 0.; precond = 0.; evalrhs = 0.; evaljac = 0.; maxIter = 0; }
    //! @brief CPU time
    double cputime;
    //! @brief preconditioning time
    double precond;
    //! @brief RHS evaluation time
    double evalrhs;
    //! @brief Jacobian evaluation time
    double evaljac;
    //! @brief Iteration count
    unsigned long maxIter;
  } _stats_ae;

  //! @brief Function to initialize implicit solver statistics
  static void _init_stats
    ( Stats& stats );

  //! @brief Function to finalize implicit solver statistics
  static void _final_stats
    ( Stats& stats );

  //! @brief Function to display implicit solver statistics
  static void _output_stats
    ( Stats const& stats, std::ostream& os = std::cout );

  //! @brief Private methods to block default compiler methods
  AEBND( AEBND<T,PMT,PVT,ExtOps...> const& );
  AEBND<T,PMT,PVT,ExtOps...>& operator=( AEBND<T,PMT,PVT,ExtOps...> const& );
};

template <typename T, typename PMT, typename PVT, typename... ExtOps>
inline
AEBND<T,PMT,PVT,ExtOps...>::AEBND
()
: BASE_AE<ExtOps...>()
{
  // Initalize state/parameter arrays
  _Ix = _Ip = nullptr;
  _PMp = _PMx = nullptr;
}

template <typename T, typename PMT, typename PVT, typename... ExtOps>
inline
AEBND<T,PMT,PVT,ExtOps...>::~AEBND
()
{
  auto it = _jac.begin();
  for( ; it != _jac.end(); ++it ){
    if( !std::get<0>(*it) ) continue;
    delete[] std::get<1>(*it);
    delete[] std::get<2>(*it);
    delete[] std::get<3>(*it);
  }
  /* DO NOT DELETE _Ix, _Ip */
  /* DO NOT DELETE _PMx, _PMp */
}

template <typename T, typename PMT, typename PVT, typename... ExtOps>
inline
void
AEBND<T,PMT,PVT,ExtOps...>::_init_stats
( Stats& stats )
{
  // Initialize statistics
  stats.reset();
  stats.cputime = -cpuclock();
}

template <typename T, typename PMT, typename PVT, typename... ExtOps>
inline
void
AEBND<T,PMT,PVT,ExtOps...>::_final_stats
( Stats& stats )
{
  // Get final CPU time
  stats.cputime += cpuclock();
}

template <typename T, typename PMT, typename PVT, typename... ExtOps>
inline 
void
AEBND<T,PMT,PVT,ExtOps...>::_output_stats
( Stats const& stats, std::ostream& os )
{
  // Statistics
  os << " No MAX ITERATIONS  " << stats.maxIter
     << std::endl
     << " CPU TIME (SEC)     " << std::fixed << std::left
                               << std::setprecision(5) << stats.cputime
     << std::endl
     << " PRECOND: " << std::setprecision(1) << stats.precond/stats.cputime*1e2
     << "%, EVALRHS: " << std::setprecision(1) << stats.evalrhs/stats.cputime*1e2
     << "%, EVALJAC: " << std::setprecision(1) << stats.evaljac/stats.cputime*1e2
     << "%" << std::endl;
  return;
}

template <typename T, typename PMT, typename PVT, typename... ExtOps>
inline
typename AEBND<T,PMT,PVT,ExtOps...>::STATUS
AEBND<T,PMT,PVT,ExtOps...>::setup
( int const NB, int const* IOR, int const* IB, int const* IPERM, std::ostream& os )
{
  _issetup = false;

  // Perform block decomposition
  if( !BASE_AE<ExtOps...>::set_block( NB, IOR, IB, IPERM, options.DISPLEVEL, os ) ) return FAILURE;
  if( BASE_AE<ExtOps...>::_singsys ) return SINGULAR;

  return _setup( os );
}

template <typename T, typename PMT, typename PVT, typename... ExtOps>
inline
typename AEBND<T,PMT,PVT,ExtOps...>::STATUS
AEBND<T,PMT,PVT,ExtOps...>::setup
( bool const updtdec, std::ostream& os )
{
  _issetup = false;

  // Perform block decomposition
  if( updtdec && !BASE_AE<ExtOps...>::set_block( options.DISPLEVEL, os ) ) return FAILURE;
  if( BASE_AE<ExtOps...>::_singsys ) return SINGULAR;

  return _setup( os );
}

template <typename T, typename PMT, typename PVT, typename... ExtOps>
inline
typename AEBND<T,PMT,PVT,ExtOps...>::STATUS
AEBND<T,PMT,PVT,ExtOps...>::_setup
( std::ostream& os )
{
  // (Re)size variable arrays
  _ndep = BASE_AE<ExtOps...>::_dep.size();
  _nvar = _ndep + BASE_AE<ExtOps...>::_var.size();
  if( _var.size() < _nvar ){
    _var.resize( _nvar );
    _Ivar.clear();  _Iref.clear();
    _PMvar.clear(); _PMref.clear(); 
  }

  // Populate variable arrays: dependents first
  unsigned ivar=0;
  for( auto id=BASE_AE<ExtOps...>::_dep.begin(); id!=BASE_AE<ExtOps...>::_dep.end(); ++id, ivar++ )
    _var[ivar] = *id;
  for( auto iv=BASE_AE<ExtOps...>::_var.begin(); iv!=BASE_AE<ExtOps...>::_var.end(); ++iv, ivar++ )
    _var[ivar] = *iv;

  // Initialize block computations according to decomposition option
  switch( options.BLKDEC ){
   case Options::DECOMPOSITION::NONE:
    _noblk = 1;
    _nblk.resize(1); _nblk[0] = _ndep;
    _pblk.resize(1); _pblk[0] = 0;
    //_fpdep.resize(_ndep); _rpdep.resize(_ndep);
    //for( unsigned i=0; i<_ndep; i++ ) _fpdep[i] = _rpdep[i] = i;
    _linblk.resize(1); _linblk[0] = BASE_AE<ExtOps...>::_linsys;
    _bwblk.resize(1);  _bwblk[0]  = BASE_AE<ExtOps...>::_bwsys;
    break;

   default:
    _noblk = BASE_AE<ExtOps...>::_noblk;
    _nblk  = BASE_AE<ExtOps...>::_nblk;
    _pblk  = BASE_AE<ExtOps...>::_pblk;
    _linblk = BASE_AE<ExtOps...>::_linblk;
    _bwblk  = BASE_AE<ExtOps...>::_bwblk;
  }
  _fpdep = BASE_AE<ExtOps...>::_fpdep;
  _rpdep = BASE_AE<ExtOps...>::_rpdep;
//#ifdef MC__AEBND_SHOW_PERMUTATION
  std::cout << "VARIABLE PERMUTATION:";
  for( unsigned const& i : _fpdep ) std::cout << " " << i;
  std::cout << std::endl;
//#endif

  _sgsys.resize(_noblk);
  auto it = _jac.begin();
  for( ; it != _jac.end(); ++it ){
    if( !std::get<0>(*it) ) continue;
    delete[] std::get<1>(*it);
    delete[] std::get<2>(*it);
    delete[] std::get<3>(*it);
  }
  _jac.resize(_noblk);
  _sgjac.resize(_noblk);

  _ldblk.resize(_noblk);
  _nblkmax = _ldblkmax = _maxop = 0;

  for( unsigned ib=0; ib<_noblk; ib++ ){
    // Set maximal block dimensions
    _ldblk[ib] = _nblk[ib];
    if( options.BLKDEC == Options::DECOMPOSITION::RECUR ){
      const unsigned int jbmin = ((options.MAXREC && ib>options.MAXREC)? ib-options.MAXREC: 0);
      for( unsigned jb=ib; jb>jbmin; jb-- ) _ldblk[ib] += _nblk[jb-1];
    }
    if( _ldblkmax < _ldblk[ib] ) _ldblkmax = _ldblk[ib];
    if( _nblkmax < _nblk[ib] ) _nblkmax = _nblk[ib];

    // Initialize block operations and Jacobian
    _sgsys[ib] = _dag->subgraph( _nblk[ib], _sys.data()+_pblk[ib] );
    switch( options.AUTODIFF ){
      case Options::DIFFERENTIATION::FSA:
        _jac[ib] = _dag->SFAD( _nblk[ib], _sys.data()+_pblk[ib],
                               _ldblk[ib], _var.data()+_pblk[ib] );
        break;
      case Options::DIFFERENTIATION::ASA:
        _jac[ib] = _dag->SBAD( _nblk[ib], _sys.data()+_pblk[ib],
                               _ldblk[ib], _var.data()+_pblk[ib] ); 
        break;
    }
    _sgjac[ib] = _dag->subgraph( std::get<0>(_jac[ib]), std::get<3>(_jac[ib]) );
    if( _maxop < _sgsys[ib].l_op.size() ) _maxop = _sgsys[ib].l_op.size();
    if( _maxop < _sgjac[ib].l_op.size() ) _maxop = _sgjac[ib].l_op.size();
  }

  // (Re)size whole-system preconditioning matrix
  _Y.resize( _ndep, _ndep );

  _issetup = true;
  return NORMAL;
}

template <typename T, typename PMT, typename PVT, typename... ExtOps>
inline
bool
AEBND<T,PMT,PVT,ExtOps...>::_init
()
{
  if( !_issetup ) return false;

  // Initializing bound-crossing by dependent branch to true
  _xlodep.assign( _ndep, true );
  _xupdep.assign( _ndep, true );
  _uniblk.assign( _noblk, false );

  // Initialize preconditionning matrix to zero
  _Y.zero();

  return true;
}

template <typename T, typename PMT, typename PVT, typename... ExtOps>
inline
bool
AEBND<T,PMT,PVT,ExtOps...>::_init
( unsigned const np, T const* Ip )
{
  if( !_init() ) return false;

  // (Re)size variable arrays
  _Ivar.resize(_nvar);
  _Iref.resize(_nvar);
  _Ix = _Ivar.data(); _Ip = _Ivar.data() + _ndep;
  _If.resize(_ndep);
  _Idfdx.resize(_noblk);
  for( unsigned ib=0; ib<_noblk; ib++ )
    _Idfdx[ib].resize( std::get<0>(_jac[ib]) );
  _Iwk.reserve( _maxop );

  _PMx = _PMp = nullptr;
  _PMf.clear();
  _PMdfdx.clear();
  _PMwk.clear();

  // Populate parameters in variable arrays
  for( unsigned i=0; i<np; i++ )
    _Ivar[_ndep+i] = _Iref[_ndep+i] = Ip[i];

  return _issetup;
}

template <typename T, typename PMT, typename PVT, typename... ExtOps>
inline
bool
AEBND<T,PMT,PVT,ExtOps...>::_init
( unsigned const np, PVT const* PMp )
{
  if( !_init() ) return false;

  // (Re)size variable arrays
  _Ivar.resize(_nvar);
  _PMvar.resize(_nvar);
  _PMref.resize(_nvar);
  _Ix = _Ivar.data(); _Ip = _Ivar.data() + _ndep;
  _Idfdx.resize(_noblk);
  for( unsigned ib=0; ib<_noblk; ib++ )
    _Idfdx[ib].resize( std::get<0>(_jac[ib]) );
  _Iwk.reserve( _maxop );

  _PMx = _PMvar.data(); _PMp = _PMvar.data() + _ndep;
  _PMf.resize(_ndep);
  _PMdfdx.resize(_noblk);
  for( unsigned ib=0; ib<_noblk; ib++ )
    _PMdfdx[ib].resize( std::get<0>(_jac[ib]) );
  _PMwk.reserve( _maxop );

  // Populate parameters in variable arrays
  for( unsigned i=0; i<np; i++ ){
    _Ivar[_ndep+i] = PMp[i].bound();
    _PMvar[_ndep+i] = _PMref[_ndep+i] = PMp[i];
  }

  return _issetup;
}

template <typename T, typename PMT, typename PVT, typename... ExtOps>
inline
bool
AEBND<T,PMT,PVT,ExtOps...>::_cvgtest
( T const& x )
const
{
  return Op<T>::diam(x) <= options.ATOL? true: false;
}

template <typename T, typename PMT, typename PVT, typename... ExtOps>
inline
bool
AEBND<T,PMT,PVT,ExtOps...>::_cvgtest
( PVT const& x )
const
{
  return Op<T>::diam(x.R()) <= options.ATOL? true: false;
}

template <typename T, typename PMT, typename PVT, typename... ExtOps>
inline
bool
AEBND<T,PMT,PVT,ExtOps...>::_crosstest
( unsigned const nx, T const* x, T const* x0, int* xlodep, int* xupdep )
const
{
  // when applying a parametric interval Newton-type method, if any improvement is observed
  // on any of the bounds, it is known that a solution curve does not cross that bound
  // see: Stuber, M. (2013), PhD Thesis, Corollary 3.4.3.
  bool cross = false;
  for( unsigned i=0; i<nx; i++ ){	
    if( xlodep[i] && Op<T>::l(x[i]) > Op<T>::l(x0[i]) ) xlodep[i] = false;
    if( xupdep[i] && Op<T>::u(x[i]) < Op<T>::u(x0[i]) ) xupdep[i] = false;
#ifdef MC__AEBND_DEBUG
    std::cout << "(" << (xlodep[i]?"T":"F") << "," << (xupdep[i]?"T":"F") << ") ";
#endif
    cross = xlodep[i] || xupdep[i];
  }
#ifdef MC__AEBND_DEBUG
  std::cout << std::endl;
#endif
  return cross;
}

template <typename T, typename PMT, typename PVT, typename... ExtOps>
inline
bool
AEBND<T,PMT,PVT,ExtOps...>::_cvgtest
( unsigned const nx, T const* x, T const* x0 )
const
{
  // Relative tolerance uses the largest improvement divided by interval width in any direction
  // Absolute tolerance uses the largest improvement in any direction
  double rtol = 0., atol = 0.;
  for( unsigned i=0; i<nx; i++ ){	
    double diamxi = Op<T>::diam(x[i]), diamx0i = Op<T>::diam(x0[i]);
    if( diamxi > 0. ){
      rtol = std::max( rtol, std::fabs( diamx0i - diamxi ) / diamxi );
      atol = std::max( atol, std::fabs( diamx0i - diamxi ) );
    }
  }
#ifdef MC__AEBND_DEBUG
  std::cout << "RTOL =" << rtol  << "  ATOL =" << atol << std::endl;
#endif
  return rtol <= options.RTOL || atol <= options.ATOL? true: false;
}

template <typename T, typename PMT, typename PVT, typename... ExtOps>
inline
bool
AEBND<T,PMT,PVT,ExtOps...>::_crosstest
( unsigned const nx, PVT const* x, PVT const* x0, int* xlodep, int* xupdep )
const
{
  // when applying a parametric intevval Newton-type method, if any improvement is observed
  // on any of the bounds, it is known that a solution curve does not cross that bound
  // see: Stuber, M. (2013), PhD Thesis, Corollary 3.4.3.
  bool cross = false;
  for( unsigned i=0; i<nx; i++ ){	
    if( xlodep[i] && Op<T>::l(x[i].R()) > Op<T>::l(x0[i].R()) ) xlodep[i] = false;
    if( xupdep[i] && Op<T>::u(x[i].R()) < Op<T>::u(x0[i].R()) ) xupdep[i] = false;
#ifdef MC__AEBND_DEBUG
    std::cout << "(" << (xlodep[i]?"T":"F") << "," << (xupdep[i]?"T":"F") << ") ";
#endif
    cross = xlodep[i] || xupdep[i];
  }
  return cross;
}

template <typename T, typename PMT, typename PVT, typename... ExtOps>
inline
bool
AEBND<T,PMT,PVT,ExtOps...>::_cvgtest
( unsigned const nx, PVT const* x, PVT const* x0 )
const
{
  if( options.PMNOREM ) return false;
  // Relative tolerance uses the largest improvement divided by interval width in any direction
  // Absolute tolerance uses the largest improvement in any direction
  double rtol = 0., atol = 0.;
#ifdef MC__AEBND_DEBUG
  std::cout << "diam(R)" << std::scientific << std::setprecision(3);
#endif
  for( unsigned i=0; i<nx; i++ ){	
    double diamxi = Op<T>::diam(x[i].R()), diamx0i = Op<T>::diam(x0[i].R());
    if( diamxi > 0. ){
      rtol = std::max( rtol, std::fabs( diamx0i - diamxi ) / diamxi );
      atol = std::max( atol, std::fabs( diamx0i - diamxi ) );
    }
#ifdef MC__AEBND_DEBUG
    std::cout << "  " << diamxi;
#endif
  }
#ifdef MC__AEBND_DEBUG
  std::cout << std::endl;
  std::cout << "RTOL =" << rtol  << "  ATOL =" << atol << std::endl;
#endif
  return rtol <= options.RTOL || atol <= options.ATOL? true: false;
}

template <typename T, typename PMT, typename PVT, typename... ExtOps>
inline
T
AEBND<T,PMT,PVT,ExtOps...>::_cancelrem
( T& x )
const
{
  return x;
}

template <typename T, typename PMT, typename PVT, typename... ExtOps>
inline
PVT
AEBND<T,PMT,PVT,ExtOps...>::_cancelrem
( PVT& x )
const
{
  return x.C().P();
}

template <typename T, typename PMT, typename PVT, typename... ExtOps>
inline
void
AEBND<T,PMT,PVT,ExtOps...>::_reference
( unsigned const ib, T* var, T* ref, T* jacvar, bool zeroref )
{
  _stats_ae.precond -= cpuclock();
  unsigned ivar = _pblk[ib];
  for( ; ivar<_pblk[ib]+_nblk[ib]; ivar++ )
    ref[ivar] = zeroref? 0.: Op<T>::mid( var[ivar] );
  for( ; ivar<_pblk[ib]+_ldblk[ib]; ivar++ )
    ref[ivar] = Op<T>::mid( var[ivar] );
  for( ; ivar<_ndep; ivar++ )
    ref[ivar] = var[ivar];
  if( var != jacvar ) throw Exceptions( Exceptions::INTERN );
  _stats_ae.precond += cpuclock();
}

template <typename T, typename PMT, typename PVT, typename... ExtOps>
inline
void
AEBND<T,PMT,PVT,ExtOps...>::_reference
( unsigned const ib, PVT* var, PVT* ref, PVT* jacvar, bool zeroref )
{
  _stats_ae.precond -= cpuclock();
  unsigned ivar = _pblk[ib];
  for( ; ivar<_pblk[ib]+_nblk[ib]; ivar++ )
    ref[ivar] = zeroref? 0.: var[ivar].C().P();
  for( ; ivar<_pblk[ib]+_ldblk[ib]; ivar++ )
    ref[ivar] = var[ivar].C().P();
  for( ; ivar<_ndep; ivar++ )
    ref[ivar] = var[ivar];
  if( var != jacvar ) throw Exceptions( Exceptions::INTERN );
  _stats_ae.precond += cpuclock();
}

template <typename T, typename PMT, typename PVT, typename... ExtOps>
inline
void
AEBND<T,PMT,PVT,ExtOps...>::_reference
( unsigned const ib, PVT* var, PVT* ref, T* jacvar, bool zeroref )
{
  _stats_ae.precond -= cpuclock();
  unsigned ivar = _pblk[ib];
  for( ; ivar<_pblk[ib]+_nblk[ib]; ivar++ ){
    ref[ivar] = zeroref? 0.: var[ivar].C().P();
    jacvar[ivar] = var[ivar].B();
  }
  for( ; ivar<_pblk[ib]+_ldblk[ib]; ivar++ ){
    ref[ivar] = var[ivar].C().P();
    jacvar[ivar] = var[ivar].B();
  }
  for( ; ivar<_ndep; ivar++ ){
    ref[ivar] = var[ivar];
    jacvar[ivar] = var[ivar].B();
  }
  _stats_ae.precond += cpuclock();
}

template <typename T, typename PMT, typename PVT, typename... ExtOps> 
template <typename U, typename V>
inline
void
AEBND<T,PMT,PVT,ExtOps...>::_precondlin
( unsigned const ib, std::vector<U>& G, std::vector<U>& b,
  std::vector<U>& wkf, U* f, U* ref, std::vector<V>& wkjacf, std::vector<V>* jacf,
  V* jacvar, std::ostream& os )
{
  try{
    if( jacf ){
      // Compute Jacobian 
      const unsigned neJAC = std::get<0>(_jac[ib]);
      const unsigned *rJAC = std::get<1>(_jac[ib]);
      const unsigned *cJAC = std::get<2>(_jac[ib]);
      const FFVar *pJAC = std::get<3>(_jac[ib]);
      V *vJAC = jacf[ib].data();
#ifdef MC__AEBND_DEBUG
      std::ofstream o_jacf( "jacf.dot", std::ios_base::out );
      _dag->dot_script( neJAC, pJAC, o_jacf );
      o_jacf.close();
      std::ofstream o_f( "f.dot", std::ios_base::out );
      _dag->dot_script( _ldblk[ib], _sys.data()+_pblk[ib], o_f );
      o_f.close();
#endif
      _stats_ae.evaljac -= cpuclock();
      _dag->eval( _sgjac[ib], wkjacf, neJAC, pJAC, vJAC, _nvar-_pblk[ib],
                  _var.data()+_pblk[ib], jacvar+_pblk[ib] );
      _stats_ae.evaljac += cpuclock();
#ifdef MC__AEBND_DEBUG
      std::cout << "Non-Preconditioned LHS Block #" << ib+1 << ":\n";
      for( unsigned ie=0; ie<neJAC; ie++ )
        std::cout << "(" << rJAC[ie] << "," << cJAC[ie] << ") : "
                  << pJAC[ie] << "=" << vJAC[ie] << std::endl;
      std::cout << std::endl;
#endif

      // Compute diagonal preconditioning block in _Y
      _stats_ae.precond -= cpuclock();
      CPPL::dgematrix mA, Y, Q, R;
      CPPL::dgbmatrix mB;
      switch( options.PRECOND ){

      case Options::PRECONDITIONING::INVMD:
        mA.resize(_nblk[ib],_nblk[ib]);
        mA.zero();
        for( unsigned ie=0; ie<neJAC; ie++ )
          if( cJAC[ie] < _nblk[ib] )
            mA(rJAC[ie],cJAC[ie]) = Op<U>::mid( vJAC[ie] );
#ifdef MC__AEBND_DEBUG
        std::cout << "Mid-point of LHS Block #" << ib+1 << ":\n" << mA << std::endl;
#endif
        if( _nblk[ib]==1 ){
          if( isequal(mA(0,0),0.) ) goto nocond;//throw Exceptions( Exceptions::PRECOND );
          Y.resize(1,1); Y(0,0) = 1./mA(0,0); break;
        }
        if( dgesv( mA, Y ) ) goto nocond;//throw Exceptions( Exceptions::PRECOND );
        break;

      case Options::PRECONDITIONING::INVMB:
        mB.resize(_nblk[ib],_nblk[ib],_bwblk[ib].first,_bwblk[ib].second);
        mB.zero();
        for( unsigned ie=0; ie<neJAC; ie++ )
          if( cJAC[ie] < _nblk[ib] )
            mB(rJAC[ie],cJAC[ie]) = Op<U>::mid( vJAC[ie] );
#ifdef MC__AEBND_DEBUG
        std::cout << "Mid-point of LHS Block #" << ib+1 << ":\n" << mB << std::endl;
#endif
        if( _nblk[ib]==1 ){
          if( isequal(mB(0,0),0.) ) goto nocond;//throw Exceptions( Exceptions::PRECOND );
          Y.resize(1,1); Y(0,0) = 1./mB(0,0); break;
        }
        if( dgbsv( mB, Y ) ) goto nocond;//throw Exceptions( Exceptions::PRECOND );
        break;

      case Options::PRECONDITIONING::QRMD:
        mA.resize(_nblk[ib],_nblk[ib]);
        mA.zero();
        for( unsigned ie=0; ie<neJAC; ie++ )
          if( cJAC[ie] < _nblk[ib] )
            mA(rJAC[ie],cJAC[ie]) = Op<U>::mid( vJAC[ie] );
#ifdef MC__AEBND_DEBUG
        std::cout << "Mid-point of LHS Block #" << ib+1 << ":\n" << mA << std::endl;
#endif
        if( _nblk[ib]==1 ){
          if( isequal(mA(0,0),0.) ) goto nocond;//throw Exceptions( Exceptions::PRECOND );
          Y.resize(1,1); Y(0,0) = 1.; break;
        }
        if( dgeqrf( mA, Q, R ) ) goto nocond;//throw Exceptions( Exceptions::PRECOND );
#ifdef MC__AEBND_DEBUG
        std::cout << "QR decomposition of LHS Block #" << ib+1 << ":\n" << Q << R << Q*R << std::endl;
#endif
        Y = t(Q);
        break;

      case Options::PRECONDITIONING::NONE:
      nocond:
        Y.resize(_nblk[ib],_nblk[ib]); Y.identity();
        break;
      }
#ifdef MC__AEBND_DEBUG
      std::cout << "Preconditioning matrix Block #" << ib+1 << ":\n" << Y << std::endl;
#endif
      for( unsigned i=0; i<_nblk[ib]; i++ )
        for( unsigned j=0; j<_nblk[ib]; j++ )
          _Y(_pblk[ib]+i,_pblk[ib]+j) = Y(i,j);        

      // Compute off-diagonal preconditioning blocks in _Y (Schur complement)
      if( _nblk[ib] < _ldblk[ib] ){
        CPPL::dgematrix mC;
        mC.resize(_nblk[ib],_ldblk[ib]-_nblk[ib]);
        mC.zero();
        std::set<unsigned> colY;
        for( unsigned ie=0; ie<neJAC; ie++ )
          if( cJAC[ie] >= _nblk[ib] ){
            auto&& it = colY.insert( cJAC[ie]-_nblk[ib] );
            mC(rJAC[ie],*it.first) = Op<U>::mid( vJAC[ie] );
            //mC(rJAC[ie],cJAC[ie]-_nblk[ib]) = Op<U>::mid( vJAC[ie] );
          }
#ifdef MC__AEBND_DEBUG
        std::cout << "Intermediate matrix Yii:\n" << Y << std::endl;
#endif
        Y *= mC;
#ifdef MC__AEBND_DEBUG
        std::cout << "Intermediate matrix Yii·mC:\n" << Y << std::endl << " { ";
        for( auto&& kt=colY.begin(); kt!=colY.end(); ++kt ) std::cout << *kt << " ";
        std::cout << "}" << std::endl;
#endif
        for( unsigned i=0; i<_nblk[ib]; i++ )
          for( unsigned j=0; j<_ldblk[ib]-_nblk[ib]; j++ ){
            _Y(_pblk[ib]+i,_pblk[ib]+_nblk[ib]+j) = 0.;
            for( auto&& kt=colY.begin(); kt!=colY.end(); ++kt )
              _Y(_pblk[ib]+i,_pblk[ib]+_nblk[ib]+j)
                -= Y(i,*kt) * _Y(_pblk[ib]+_nblk[ib]+*kt,_pblk[ib]+_nblk[ib]+j);
          }
      }
#ifdef MC__AEBND_SHOW_PRECONDITIONING
        std::cout << std::scientific << std::setprecision(5)
                  << "Full preconditioning matrix:\n" << _Y << std::endl;
#endif
      //_stats_ae.precond += cpuclock();

      // Set G = Y*dfdx
      for( unsigned i=0; i<_nblk[ib]; i++ ){
        for( unsigned j=0; j<_ldblk[ib]; j++ )
          G[_ndx(i,j,_ldblk[ib])] = 0.;
        const int jbmin = ((options.MAXREC && ib>options.MAXREC)? ib-options.MAXREC: 0);
        //std::cout << "jbmin: " << jbmin << std::endl;
        for( int jb=ib; jb>=jbmin; jb-- ){
          const unsigned neJAC = std::get<0>(_jac[jb]);
          const unsigned *rJAC = std::get<1>(_jac[jb]);
          const unsigned *cJAC = std::get<2>(_jac[jb]);
          const V *vJAC = jacf[jb].data();
          for( unsigned ie=0; ie<neJAC; ie++ ){
            const unsigned j = _pblk[jb]-_pblk[ib]+cJAC[ie];
            const unsigned k = _pblk[jb]-_pblk[ib]+rJAC[ie];
            //std::cout << "(i,j,k): " << i << " " << j << " " << k << std::endl;
	    if( j < _ldblk[ib] )
              G[_ndx(i,j,_ldblk[ib])] += _Y(_pblk[ib]+i,_pblk[ib]+k) * vJAC[ie];
	    //else
	    //  std::cout << " -- not used" << std::endl;
          }
          if( _nblk[ib] == _ldblk[ib] ) break; // diagonal preconditioning
        }
      }
#ifdef MC__AEBND_DEBUG
      std::cout << "Preconditioned LHS Block #" << ib+1 << ":\n";
      for( unsigned i=0; i<_nblk[ib]; i++ ){
        for( unsigned j=0; j<_ldblk[ib]; j++ )
          std::cout << G[_ndx(i,j,_ldblk[ib])] << "  ";
        std::cout << std::endl;
      }
      std::cout << std::endl;
#endif
      _stats_ae.precond += cpuclock();
    }

    if( f ){
      //setting b = -Y*f(x)
      _stats_ae.evalrhs -= cpuclock();
      _dag->eval( _sgsys[ib], wkf, _nblk[ib], _sys.data()+_pblk[ib], f+_pblk[ib],
                  _nvar-_pblk[ib], _var.data()+_pblk[ib], ref+_pblk[ib] );
      // Update function value for previous block (if any)
//      if( ib && _nblk[ib] < _ldblk[ib] )
//        _dag->eval( _sgsys[ib-1], wkf, _nblk[ib-1], _sys.data()+_pblk[ib-1], f+_pblk[ib-1],
//                    _nvar-_pblk[ib-1], _var.data()+_pblk[ib-1], ref+_pblk[ib-1] );  
      const int jbmin = ((options.MAXREC && ib>options.MAXREC)? ib-options.MAXREC: 0);
      for( int jb=ib-1; jb>=jbmin; jb-- ){
        _dag->eval( _sgsys[jb], wkf, _nblk[jb], _sys.data()+_pblk[jb], f+_pblk[jb],
                    _nvar-_pblk[jb], _var.data()+_pblk[jb], ref+_pblk[jb] );
        if( !jbmin ) break; // no need to recompute ealier blocks since identical
        //std::cout << "NO BREAK!" << std::endl; int dum; std::cin >> dum;
      }
      _stats_ae.evalrhs += cpuclock();
#ifdef MC__AEBND_DEBUG
      //_dag->output( _sgsys[ib] );
      //std::cout << "Intermediates in Block #" << ib+1 << ":\n";
      //for (unsigned i=0; i<_sgsys[ib].size(); i++ ) std::cout << opf[i] << std::endl;
      std::cout << "reference in Block #" << ib+1 << ":\n";
      for (unsigned i=0; i<_nvar-_pblk[ib]; i++ )
        std::cout << _var[_pblk[ib]+i] << " = " << ref[_pblk[ib]+i] << std::endl;
      std::cout << "Non-Preconditioned RHS Block #" << ib+1 << ":\n";
      for (unsigned i=0; i<_ldblk[ib]; i++ )
        std::cout << _sys[_pblk[ib]+i] << ": " << f[_pblk[ib]+i] << std::endl;
      std::cout << std::endl;
      { int dum; std::cin >> dum; }
#endif
      _stats_ae.precond -= cpuclock();
      for (unsigned i=0; i<_nblk[ib]; i++ ){
        b[i] = 0.;
        for( unsigned j=0; j<_ldblk[ib]; j++)
          b[i] -= _Y(_pblk[ib]+i,_pblk[ib]+j) * f[_pblk[ib]+j];
      }
      _stats_ae.precond += cpuclock();
#ifdef MC__AEBND_DEBUG
      std::cout << "Preconditioned RHS Block #" << ib+1 << ":\n";
      for (unsigned i=0; i<_nblk[ib]; i++ ) std::cout << b[i] << std::endl;
      std::cout << std::endl;
#endif
    }
  }
  catch(...){
    throw Exceptions( Exceptions::PRECOND );
  }
}

template <typename T, typename PMT, typename PVT, typename... ExtOps>
template <typename U, typename V>
inline
typename AEBND<T,PMT,PVT,ExtOps...>::STATUS
AEBND<T,PMT,PVT,ExtOps...>::_gs
( unsigned const ib, U* var, std::vector<U>& wkf, U* f, U* ref,
  std::vector<V>& wkjacf, std::vector<V>* jacf, V* jacvar, bool usekraw,
  std::ostream& os )
{
  U* varblk = var + _pblk[ib];
  U* refblk = ref + _pblk[ib];

  // Display
  if( options.DISPLEVEL >= 2  ){
    os << std::endl << "Block #" << ib+1 << "  GS Initial:\n";
    for( unsigned i=0; i<_nblk[ib]; i++ )
     os << _var[_pblk[ib]+i] << " = " << varblk[i] << std::endl;
  }

  std::vector<U> G(_nblk[ib]*_ldblk[ib]), b(_nblk[ib]), varblk0; 
  for( unsigned iter=0; iter<options.MAXIT; iter++ ){
    if( _stats_ae.maxIter <= iter ) _stats_ae.maxIter = iter+1;

    // Cancel remainder
    for( unsigned i=0; options.PMNOREM && i<_nblk[ib]; i++ ){
       varblk[i] = _cancelrem( varblk[i] );
#ifdef  MC__AEBND_DEBUG
      std::cout << "X[" << i << "] = " << varblk[i];
#endif
    }

    try{
      // Update reference
      _reference( ib, var, ref, jacvar );

      // (Re)compute preconditionned LHS matrix (only if needed) and RHS vector
      _precondlin( ib, G, b, wkf, f, ref, wkjacf,
                   (!iter||!_linblk[ib]?jacf:0), jacvar, os );

      // Keep track of current iterate
      varblk0.assign( varblk, varblk+_nblk[ib] );

      for( unsigned i=0; i<_nblk[ib]; i++ ){
        //if( _cvgtest( varblk[i] ) ){
        //  std::cout << "converged: " << varblk[i] << std::endl;
        //  continue;
        //}
        U Xk, temp(0.); 
        // Apply componentwise Gauss-Seidel step
        if( !usekraw
         && ( Op<U>::l(G[_ndx(i,i,_ldblk[ib])]) > 0.
           || Op<U>::u(G[_ndx(i,i,_ldblk[ib])]) < 0. ) ){
          for( unsigned j=0; j<_ldblk[ib]; j++ )
            if( j != i ) temp += G[_ndx(i,j,_ldblk[ib])] * ( varblk[j] - refblk[j] );
          Xk = refblk[i] + ( b[i] - temp ) / G[_ndx(i,i,_ldblk[ib])];
        }
        // Apply componentwise Krawczyk step
        else{
          for( unsigned j=0; j<_ldblk[ib]; j++ ){
            if( j != i )  temp -= G[_ndx(i,j,_ldblk[ib])] * ( varblk[j] - refblk[j] );
            else temp += ( 1. - G[_ndx(i,i,_ldblk[ib])] ) * ( varblk[j] - refblk[j] );
          }
          Xk = refblk[i] + b[i] + temp;
        }
#ifdef  MC__AEBND_DEBUG
        std::cout << "varblk[" << i << "] = " << varblk[i];
        std::cout << "refblk[" << i << "] = " << refblk[i];
        std::cout << "b[" << i << "] = " << b[i];
        std::cout << "temp = " << temp;
        std::cout << "X[" << i << "] = " << Xk;
//        { int dum; std::cout << "paused"; std::cin >> dum; }
#endif
        // Remainder operations
        if( options.PMNOREM ) varblk[i] = _cancelrem( Xk );
        else if( !options.INTERBND ) varblk[i] = Xk;
        else if( !Op<U>::inter( varblk[i], Xk, varblk[i] ) ){
#ifdef  MC__AEBND_SHOW_INTER_FAIL
          os << _var[_pblk[ib]+i] << " = " << Xk << " && " << varblk[i] << std::endl;
#endif
#ifdef  MC__AEBND_IGNORE_INTER_FAIL
          varblk[i] = Xk;
#else
          return EMPTY;
#endif
        }
      }

#ifdef MC__AEBND_DEBUG
      std::cout << "Iter #" << iter << " bounds:\n";
      for( unsigned i=0; i<_nblk[ib]; i++ ) std::cout << varblk[i] << std::endl;
#endif

      // Display
      if( options.DISPLEVEL >= 2  ){
        os << std::endl << "Block #" << ib+1 << "  GS Iteration #" << iter+1 << ":\n";
        for( unsigned i=0; i<_nblk[ib]; i++ )
          os << _var[_pblk[ib]+i] << " = " << varblk[i] << std::endl;
      }

      // Update solution branch uniqueness
      if( ( _nblk[ib] == 1 && ( Op<U>::l(G[_ndx(0,0,_ldblk[ib])]) > 0
                             || Op<U>::u(G[_ndx(0,0,_ldblk[ib])]) < 0 ) )
        || !_crosstest( _nblk[ib], varblk, varblk0.data(), _xlodep.data(), _xupdep.data() ) )
        _uniblk[ib] = true;

      // Check convergence
      if( _cvgtest( _nblk[ib], varblk, varblk0.data() ) ) break;
    }
    catch(...){
      if( options.DISPLEVEL >= 2  )
        os << std::endl << "Block #" << ib+1 << "  GS Failed\n";
      return FAILURE;
    }
  }

  if( options.DISPLEVEL >= 2  )
    std::cout << "Unique Solution: " << (_uniblk[ib]?"T\n":"F\n");

  return NORMAL;
}
template <typename T, typename PMT, typename PVT, typename... ExtOps>
template <typename U>
inline
typename AEBND<T,PMT,PVT,ExtOps...>::STATUS
AEBND<T,PMT,PVT,ExtOps...>::_ge
( unsigned const ib, U* var, std::vector<U>& wk, U* f, U* ref,
  std::vector<U>* jacf, std::ostream& os )
{
  if( !_stats_ae.maxIter ) _stats_ae.maxIter = 1;
  U* varblk = var + _pblk[ib];
  U* refblk = ref + _pblk[ib];

  try{
    // Compute preconditionned LHS matrix and RHS vector
    std::vector<U> G(_nblk[ib]*_ldblk[ib]), b(_nblk[ib]); 
    //std::cout << "Block: " << ib << "  G is " << _nblk[ib] << "x" << _ldblk[ib] << std::endl;
    _reference( ib, var, ref, var, true ); // zero reference to compute right-hand-side vector
    _precondlin( ib, G, b, wk, f, ref, wk, jacf, var, os );

    // Perform RHS vector correction for block recursivity
    for( unsigned k=0; k<_nblk[ib]; k++ )
      for( unsigned j=_nblk[ib]; j<_ldblk[ib]; j++ )
        b[k] -= G[_ndx(k,j,_ldblk[ib])] * ( varblk[j] - refblk[j] );

    // Perform forward elimination
    for( unsigned k=0; k<_nblk[ib]-1; k++ ){
      if( Op<U>::l( G[_ndx(k,k,_ldblk[ib])] ) <= 0. 
       && Op<U>::u( G[_ndx(k,k,_ldblk[ib])] ) >= 0. ) return SINGULAR;
      for( unsigned i=k+1; i<_nblk[ib]; i++ ){
        U factor = G[_ndx(i,k,_ldblk[ib])] / G[_ndx(k,k,_ldblk[ib])];
        for( unsigned j=k+1; j<_nblk[ib]; j++ )
          G[_ndx(i,j,_ldblk[ib])] -= factor * G[_ndx(k,j,_ldblk[ib])];
        b[i] -= factor * b[k];
      }
    }
#ifdef MC__AEBND_DEBUG
    std::cout << "LHS Block #" << ib+1 << " After Forward Elimination:\n";
    for( unsigned i=0; i<_nblk[ib]; i++ ){
      for( unsigned j=0; j<_nblk[ib]; j++ )
        std::cout << G[_ndx(i,j,_ldblk[ib])] << "  ";
      std::cout << std::endl;
    }
    std::cout << std::endl;
    std::cout << "RHS Block #" << ib+1 << " After Forward Elimination:\n";
    for (unsigned i=0; i<_nblk[ib]; i++ ) std::cout << b[i] << std::endl;
    std::cout << std::endl;
#endif

    // Check for singularity
    if( Op<U>::l( G[_ndx(_nblk[ib]-1,_nblk[ib]-1,_ldblk[ib])] ) <= 0. 
     && Op<U>::u( G[_ndx(_nblk[ib]-1,_nblk[ib]-1,_ldblk[ib])] ) >= 0. ){
      if( options.DISPLEVEL >= 2  )
        os << std::endl << "Block #" << ib+1 << "  GE Singular\n";
      return SINGULAR;
    }

    // Perform backward substitution
    varblk[_nblk[ib]-1] = b[_nblk[ib]-1] / G[_ndx(_nblk[ib]-1,_nblk[ib]-1,_ldblk[ib])];
    for( unsigned i=_nblk[ib]-1; i>0; i-- ){
      U sum = b[i-1];
      for( unsigned j=i; j<_nblk[ib]; j++ )
        sum -= G[_ndx(i-1,j,_ldblk[ib])] * varblk[j];
      varblk[i-1] = sum / G[_ndx(i-1,i-1,_ldblk[ib])];
    }

    // Display
    if( options.DISPLEVEL >= 2  ){
      os << std::endl << "Block #" << ib+1 << "  GE Final:\n";
      for( unsigned i=0; i<_nblk[ib]; i++ )
        os << _var[_pblk[ib]+i] << " = " << varblk[i] << std::endl;
    }
  }
  catch(...){
    if( options.DISPLEVEL >= 2  )
      os << std::endl << "Block #" << ib+1 << "  GE Failed\n";
    return FAILURE;
  }
  
  return NORMAL;
}

//! @fn template <typename T, typename PMT, typename PVT, typename... ExtOps> inline typename AEBND<T,PMT,PVT,ExtOps...>::STATUS AEBND<T,PMT,PVT,ExtOps...>::solve( const T*Ip, T*Ix, const T*Ix0=0, std::ostream&os=std::cout )
//!
//! This function computes an interval enclosure of the solution set of
//! the parametric AEs:
//!   - <a>Ip</a> [input] interval parameter set
//!   - <a>Ix</a> [output] interval state enclosure
//!   - <a>Ix0</a> [input] a priori interval state bounds (default: NULL)
//!   - <a>os</a> [input] output stream (default: std::cout)
//! .
//! The return value is the status.
template <typename T, typename PMT, typename PVT, typename... ExtOps>
inline
typename AEBND<T,PMT,PVT,ExtOps...>::STATUS
AEBND<T,PMT,PVT,ExtOps...>::solve
( T const* Ip, T* Ix, T const* Ix0, std::ostream& os )
{
  _init_stats( _stats_ae );
  _iblk=0;
  if( BASE_AE<ExtOps...>::_singsys )  return SINGULAR;
  if( !_init( BASE_AE<ExtOps...>::_var.size(), Ip ) ) return FAILURE;

  // Loop over each block
  STATUS flag = NORMAL;
  for( ; _iblk<_noblk; _iblk++ ){
    switch( options.BOUNDER ){

    // Automatic selection of bounder
    case Options::ALGORITHM::AUTO:
      // No a priori box supplied
      if( !Ix0 ){
        if( !_linblk[_iblk] ) throw Exceptions( Exceptions::APRIORI );
        flag = _ge( _iblk, _Ivar.data(), _Iwk, _If.data(), _Iref.data(),
                    _Idfdx.data(), os );
      }
  
      // A priori box supplied and linear case
      else if( _linblk[_iblk] ){
        if( options.DISPLEVEL >= 2  ){
          os << std::endl << "Block #" << _iblk+1 << "  AUTO Initial:\n";
          for( unsigned i=0; i<_nblk[_iblk]; i++ )
            os << _var[_pblk[_iblk]+i] << " = " << Ix0[_rpdep[_pblk[_iblk]+i]] << std::endl;
        }
        flag = _ge( _iblk, _Ivar.data(), _Iwk, _If.data(), _Iref.data(),
                    _Idfdx.data(), os );

        // abnormal termination
        if( flag == SINGULAR || flag == FAILURE ){
          for( unsigned i=0; i<_nblk[_iblk]; i++ )
            _Ivar[_pblk[_iblk]+i] = Ix0[_rpdep[_pblk[_iblk]+i]];
          break;
        }
        else if( flag != NORMAL )
          break;

        // check whether supplied bounds are better than GE bounds in parts
        bool improvecheck = false;
        for( unsigned i=0; i<_nblk[_iblk]; i++ ){
          if( Op<T>::l(_Ivar[_pblk[_iblk]+i]) < Op<T>::l(Ix0[_rpdep[_pblk[_iblk]+i]])
           || Op<T>::u(_Ivar[_pblk[_iblk]+i]) > Op<T>::u(Ix0[_rpdep[_pblk[_iblk]+i]]) ){
            if( !Op<T>::inter( _Ivar[_pblk[_iblk]+i], _Ivar[_pblk[_iblk]+i],
                               Ix0[_rpdep[_pblk[_iblk]+i]] ) ) return EMPTY;
            improvecheck = true;
          }
        }
        if( improvecheck )
          flag = _gs( _iblk, _Ivar.data(), _Iwk, _If.data(), _Iref.data(),
                      _Iwk, _Idfdx.data(), _Ivar.data(), false, os );
      }

      // A priori box supplied and nonlinear case
      else{
        for( unsigned i=0; i<_nblk[_iblk]; i++ )
          _Ivar[_pblk[_iblk]+i] = Ix0[_rpdep[_pblk[_iblk]+i]];
        flag = _gs( _iblk, _Ivar.data(), _Iwk, _If.data(), _Iref.data(),
                    _Iwk, _Idfdx.data(), _Ivar.data(), false, os );
      }
      break;

    // Gauss Elimination method
    case Options::ALGORITHM::GE:
      // Check system is linear in dependents
      if( !_linblk[_iblk] ) throw Exceptions( Exceptions::GAUSSEL );
      //for( unsigned irep=0; irep<1000; ++irep )
      flag = _ge( _iblk, _Ivar.data(), _Iwk, _If.data(), _Iref.data(),
                  _Idfdx.data(), os );
      break;

    // Gauss-Seidel method
    case Options::ALGORITHM::GS: case Options::ALGORITHM::GSS:
      // Check initial box is given
      if( !Ix0 ) throw Exceptions( Exceptions::APRIORI );
      for( unsigned i=0; i<_nblk[_iblk]; i++ )
        _Ivar[_pblk[_iblk]+i] = Ix0[_rpdep[_pblk[_iblk]+i]];
      //for( unsigned irep=0; irep<1000; ++irep )
      flag = _gs( _iblk, _Ivar.data(), _Iwk, _If.data(), _Iref.data(),
                  _Iwk, _Idfdx.data(), _Ivar.data(), false, os );
      break;

    // Krawczyk method
    case Options::ALGORITHM::KRAW: case Options::ALGORITHM::KRAWS:
      // Check initial box is given
      if( !Ix0 ) throw Exceptions( Exceptions::APRIORI );
      for( unsigned i=0; i<_nblk[_iblk]; i++ )
        _Ivar[_pblk[_iblk]+i] = Ix0[_rpdep[_pblk[_iblk]+i]];
      flag = _gs( _iblk, _Ivar.data(), _Iwk, _If.data(), _Iref.data(),
                  _Iwk, _Idfdx.data(), _Ivar.data(), true, os );
      break;
    } // end switch

    // Copy result bounds for current block
    if( flag == NORMAL )
      for( unsigned i=0; i<_nblk[_iblk]; i++ )
        Ix[_rpdep[_pblk[_iblk]+i]] = _Ivar[_pblk[_iblk]+i];
    else
      for( unsigned i=0; i<_nblk[_iblk]; i++ )
        if( Ix0 ) Ix[_rpdep[_pblk[_iblk]+i]] = Ix0[_rpdep[_pblk[_iblk]+i]];
        else      Ix[_rpdep[_pblk[_iblk]+i]] = T(-INF,INF);
  }

  // Copy result bounds for unsucessful blocks
  //for( ; _iblk<_noblk; _iblk++ )
  //  for( unsigned i=0; i<_nblk[_iblk]; i++ )
  //    Ix[_rpdep[_pblk[_iblk]+i]] = Ix0? Ix0[_rpdep[_pblk[_iblk]+i]]: T(-INF,INF);

  _final_stats( _stats_ae );
  if( flag != NORMAL && flag != FAILURE ) return flag;

  for( unsigned i=0; i<_ndep; i++ ) Ix[i] = _Ivar[_fpdep[i]];
  if( options.DISPLEVEL >= 1 ){
    std::cout << std::endl << "Final Bounds:\n";
    for( unsigned i=0; i<_ndep; i++ )
      os << _var[_fpdep[i]] << " = " << Ix[i] << std::endl;
    _output_stats( _stats_ae, os );
  }
  return NORMAL;
}

//! @fn template <typename T, typename PMT, typename PVT, typename... ExtOps> inline typename AEBND<T,PMT,PVT,ExtOps...>::STATUS AEBND<T,PMT,PVT,ExtOps...>::solve( const PVT*PMp, PVT*PMx, const PVT*PMx0=0, std::ostream&os=std::cout )
//!
//! This function computes an interval enclosure of the solution set of
//! the parametric AEs:
//!   - <a>PMp</a> [input] polynomial model of parameter set
//!   - <a>PMx</a> [output] polynomial model of state enclosure
//!   - <a>PMx0</a> [input] a priori polynomial model of state bounds (default: NULL)
//!   - <a>os</a> [input] output stream (default: std::cout)
//! .
//! The return value is the status.
template <typename T, typename PMT, typename PVT, typename... ExtOps>
template <typename INI>
inline
typename AEBND<T,PMT,PVT,ExtOps...>::STATUS
AEBND<T,PMT,PVT,ExtOps...>::solve
( PVT const* PMp, PVT* PMx, INI const* PMx0, std::ostream& os )
{
  _init_stats( _stats_ae );
  _iblk=0;
  if( BASE_AE<ExtOps...>::_singsys )   return SINGULAR;
  if( !_init( BASE_AE<ExtOps...>::_var.size(), PMp ) ) return FAILURE;

  // Loop over each block
  STATUS flag = NORMAL;
  for( ; _iblk<_noblk; _iblk++ ){
    switch( options.BOUNDER ){

    // Automatic selection of bounder
    case Options::ALGORITHM::AUTO:
      // No a priori box supplied
      if( !PMx0 ){
        if( !_linblk[_iblk] ) throw Exceptions( Exceptions::APRIORI );
        flag = _ge( _iblk, _PMvar.data(), _PMwk, _PMf.data(), _PMref.data(),
                    _PMdfdx.data(), os );
      }
  
      // A priori box supplied and linear case
      else if( _linblk[_iblk] ){
        flag = _ge( _iblk, _PMvar.data(), _PMwk, _PMf.data(), _PMref.data(),
                    _PMdfdx.data(), os );
        if( flag != NORMAL ) break;

        // check whether supplied bounds are better than GE bounds in parts
        bool improvecheck = false;
        for( unsigned i=0; i<_nblk[_iblk]; i++ ){
          if( Op<PVT>::l(_PMvar[_pblk[_iblk]+i]) < Op<PVT>::l(PMx0[_rpdep[_pblk[_iblk]+i]])
           || Op<PVT>::u(_PMvar[_pblk[_iblk]+i]) > Op<PVT>::u(PMx0[_rpdep[_pblk[_iblk]+i]]) ){
            if( !Op<PVT>::inter( _PMvar[_pblk[_iblk]+i], _PMvar[_pblk[_iblk]+i],
                                 PMx0[_rpdep[_pblk[_iblk]+i]] ) ) return EMPTY;
            improvecheck = true;
          }
        }
        if( improvecheck )
          flag = _gs( _iblk, _PMvar.data(), _PMwk, _PMf.data(), _PMref.data(),
                      _Iwk, _Idfdx.data(), _Ivar.data(), false, os );
      }

      // A priori box supplied and nonlinear case
      else{
        for( unsigned i=0; i<_nblk[_iblk]; i++ )
          _PMvar[_pblk[_iblk]+i] = PMx0[_rpdep[_pblk[_iblk]+i]];
        flag = _gs( _iblk, _PMvar.data(), _PMwk, _PMf.data(), _PMref.data(),
                    _PMwk, _PMdfdx.data(), _PMvar.data(), false, os );
      }
      break;

    // Gauss Elimination method
    case Options::ALGORITHM::GE:
      // Check system is linear in dependents
      if( !_linblk[_iblk] ) throw Exceptions( Exceptions::GAUSSEL );
      flag = _ge( _iblk, _PMvar.data(), _PMwk, _PMf.data(), _PMref.data(),
                  _PMdfdx.data(), os );
      break;

    // Gauss-Seidel method
    case Options::ALGORITHM::GS:
      // Check initial box is given
      if( !PMx0 ) throw Exceptions( Exceptions::APRIORI );
      for( unsigned i=0; i<_nblk[_iblk]; i++ )
        _PMvar[_pblk[_iblk]+i] = PMx0[_rpdep[_pblk[_iblk]+i]];
      flag = _gs( _iblk, _PMvar.data(), _PMwk, _PMf.data(), _PMref.data(),
                  _PMwk, _PMdfdx.data(), _PMvar.data(), false, os );
      break;

    // Simplified Gauss-Seidel method
    case Options::ALGORITHM::GSS:
      // Check initial box is given
      if( !PMx0 ) throw Exceptions( Exceptions::APRIORI );
      for( unsigned i=0; i<_nblk[_iblk]; i++ )
        _PMvar[_pblk[_iblk]+i] = PMx0[_rpdep[_pblk[_iblk]+i]];
      flag = _gs( _iblk, _PMvar.data(), _PMwk, _PMf.data(), _PMref.data(),
                  _Iwk, _Idfdx.data(), _Ivar.data(), false, os );
      break;

    // Krawczyk method
    case Options::ALGORITHM::KRAW:
      // Check initial box is given
      if( !PMx0 ) throw Exceptions( Exceptions::APRIORI );
      for( unsigned i=0; i<_nblk[_iblk]; i++ )
        _PMvar[_pblk[_iblk]+i] = PMx0[_rpdep[_pblk[_iblk]+i]];
      flag = _gs( _iblk, _PMvar.data(), _PMwk, _PMf.data(), _PMref.data(),
                  _PMwk, _PMdfdx.data(), _PMvar.data(), true, os );
      break;

    // Simplified Krawczyk method
    case Options::ALGORITHM::KRAWS:
      // Check initial box is given
      if( !PMx0 ) throw Exceptions( Exceptions::APRIORI );
      for( unsigned i=0; i<_nblk[_iblk]; i++ )
        _PMvar[_pblk[_iblk]+i] = PMx0[_rpdep[_pblk[_iblk]+i]];
      flag = _gs( _iblk, _PMvar.data(), _PMwk, _PMf.data(), _PMref.data(),
                  _Iwk, _Idfdx.data(), _Ivar.data(), true, os );
      break;
    } // end switch

    // Copy result bounds for current block
    if( flag == NORMAL )
      for( unsigned i=0; i<_nblk[_iblk]; i++ )
        PMx[_rpdep[_pblk[_iblk]+i]] = _PMvar[_pblk[_iblk]+i];
    else
      for( unsigned i=0; i<_nblk[_iblk]; i++ )
        if( PMx0 ) PMx[_rpdep[_pblk[_iblk]+i]] = PMx0[_rpdep[_pblk[_iblk]+i]];
        else       PMx[_rpdep[_pblk[_iblk]+i]] = T(-INF,INF);
  }

  // Copy result bounds for unsucessful blocks
  //for( ; _iblk<_noblk; _iblk++ )
  //  for( unsigned i=0; i<_nblk[_iblk]; i++ )
  //    if( PMx0 ) PMx[_rpdep[_pblk[_iblk]+i]] = PMx0[_rpdep[_pblk[_iblk]+i]];
  //    else       PMx[_rpdep[_pblk[_iblk]+i]] = T(-INF,INF);

  _final_stats( _stats_ae );
  if( flag != NORMAL && flag != FAILURE ) return flag;

  for( unsigned i=0; i<_ndep; i++ ) PMx[i] = _PMvar[_fpdep[i]];
  if( options.DISPLEVEL >= 1 ){
    std::cout << std::endl << "Final Bounds:\n";
    for( unsigned i=0; i<_ndep; i++ )
      os << _var[_fpdep[i]] << " = " << PMx[i] << std::endl;
    _output_stats( _stats_ae, os );
  }
  return NORMAL;
}
#if 0
template <typename T, typename PMT, typename PVT, typename... ExtOps>
inline typename AEBND<T,PMT,PVT,ExtOps...>::STATUS
AEBND<T,PMT,PVT,ExtOps...>::_ge
( std::vector<FFVar>&X, std::ostream&os )
{
 try{
  _stats_ae.maxIter = 1;
  for( unsigned ib=0; ib<_noblk; ib++ ){
    // Construct LHS matrix and RHS vector
    const unsigned neJAC = std::get<0>( _jac[ib] );
    const unsigned *rJAC = std::get<1>( _jac[ib] );
    const unsigned *cJAC = std::get<2>( _jac[ib] );
    const FFVar *pJAC = std::get<3>( _jac[ib] );
    FFVar*Xblk = X.data() + _pblk[ib];
    const FFVar*Fblk = _sys.data() + _pblk[ib];
    std::vector<FFVar> A( _nblk[ib]*_nblk[ib], 0. ), b( _nblk[ib] );
#ifdef MC__AEBND_DEBUG
    std::cout << "\nBefore forward substitution:\n";
#endif
    //for( unsigned i=0; i<_nblk[ib]*_nblk[ib]; i++ ){
    //  A[i] = _jac[ib][i];
#ifdef MC__AEBND_DEBUG
    //  std::cout << "A[" << i << "]: " << A[i] << std::endl;
#endif
    //}
    for( unsigned ie=0; ie<neJAC; ie++ ){
      if( cJAC[ie] < _nblk[ib] )
        A[_ndx(rJAC[ie],cJAC[ie],_nblk[ib])] = pJAC[ie];
#ifdef MC__AEBND_DEBUG
      std::cout << "A[" << _ndx(rJAC[ie],cJAC[ie],_nblk[ib]) << "]: "
                << A[_ndx(rJAC[ie],cJAC[ie],_nblk[ib])] << std::endl;
#endif
    }
    std::vector<FFVar> zeros( _nblk[ib], 0. );
    const FFVar*Fblk0 = _dag->compose( _nblk[ib], Fblk, _nblk[ib],
                                       _var.data()+_pblk[ib], zeros.data() );
    for( unsigned i=0; i<_nblk[ib]; i++ ){
      b[i] = -Fblk0[i];
#ifdef MC__AEBND_DEBUG
      std::cout << "b[" << i << "]: " << b[i] << std::endl;
#endif
    }
    delete[] Fblk0;

    // Proceed with forward substitution
    for( unsigned k=0; k<_nblk[ib]-1; k++ ){
      for( unsigned i=k+1; i<_nblk[ib]; i++ ){
        FFVar factor = A[_ndx(i,k,_nblk[ib])] / A[_ndx(k,k,_nblk[ib])];
        for( unsigned j=k+1; j<_nblk[ib]; j++ )
          A[_ndx(i,j,_nblk[ib])] -= factor * A[_ndx(k,j,_nblk[ib])];
        b[i] -= factor * b[k];
      }
    }
#ifdef MC__AEBND_DEBUG
    std::cout << "\nAfter forward substitution:\n";
    for( unsigned i=0; i<_nblk[ib]*_nblk[ib]; i++ )
      std::cout << "A[" << i << "]: " << A[i] << std::endl;
    for( unsigned i=0; i<_nblk[ib]; i++ )
      std::cout << "b[" << i << "]: " << b[i] << std::endl;
#endif

    // Proceed with backward elimination
    Xblk[_nblk[ib]-1] = b[_nblk[ib]-1] / A[_ndx(_nblk[ib]-1,_nblk[ib]-1,_nblk[ib])];
    for( unsigned i=_nblk[ib]; i>0; i-- ){
      FFVar sum = b[i-1];
      for( unsigned j=i; j<_nblk[ib]; j++ )
        sum -= A[_ndx(i-1,j,_nblk[ib])] * Xblk[j];
      Xblk[i-1] = sum / A[_ndx(i-1,i-1,_nblk[ib])];
    }

    // Compose with Solutions of previous blocks
    if( ib ){
#ifdef MC__AEBND_DEBUG
      std::cout << "\nComposition with previous solutions:\n";
      for( unsigned i=0; i<_ndep-_pblk[ib]-_nblk[ib]; i++ )
        std::cout << (_var.data()+_pblk[ib]+_nblk[ib])[i] << " -> " << (Xblk+_nblk[ib])[i] << std::endl;
#endif
      const FFVar*Xblkcomp = _dag->compose( _nblk[ib], Xblk, _ndep-_pblk[ib]-_nblk[ib],
        _var.data()+_pblk[ib]+_nblk[ib], Xblk+_nblk[ib] );
      for( unsigned i=0; i<_nblk[ib]; i++ ) Xblk[i] = Xblkcomp[i];
      delete[] Xblkcomp;
    }
#ifdef MC__AEBND_DEBUG
    std::cout << "\nAfter backward elimination:\n";
    for( unsigned i=0; i<_nblk[ib]*_nblk[ib]; i++ )
      std::cout << "A[" << i << "] :" << A[i] << std::endl;
    for( unsigned i=0; i<_nblk[ib]; i++ )
      std::cout << "b[" << i << "] :" << b[i] << std::endl;
    for( unsigned i=0; i<_nblk[ib]; i++ )
      std::cout << "X[" << i << "] :" << Xblk[i] << std::endl;
#endif
   }
  }
  catch(...){
    return FAILURE;
  }
#ifdef MC__AEBND_DEBUG
    std::cout << *_dag << std::endl;
#endif
  return NORMAL;
}

//! @fn template <typename T, typename PMT, typename PVT, typename... ExtOps> inline typename AEBND<T,PMT,PVT,ExtOps...>::STATUS AEBND<T,PMT,PVT,ExtOps...>::bounds(
//! FFVar*X, std::ostream&os=std::cout )
//!
//! This function computes a symbolic expression of the solution set of
//! the parametric AEs:
//!   - <a>X</a> [output] symbolic solution
//!   - <a>os</a> [input] output stream (default: std::cout)
//! .
//! The return value is the status.
template <typename T, typename PMT, typename PVT, typename... ExtOps>
inline typename AEBND<T,PMT,PVT,ExtOps...>::STATUS
AEBND<T,PMT,PVT,ExtOps...>::solve
( FFVar*X, std::ostream&os )
{
  _init_stats( _stats_ae );
  if( !_issetup )   return FAILURE;
  if( BASE_AE<ExtOps...>::_singsys ) return SINGULAR;
  std::vector<FFVar> var(_ndep);
  STATUS flag = NORMAL;

  switch( options.BOUNDER ){

  // Gauss Elimination method
  case Options::ALGORITHM::AUTO:
  case Options::ALGORITHM::GE:
    if( !BASE_AE<ExtOps...>::_linsys ) throw Exceptions( Exceptions::DAG );
    flag = _ge( var, os );
    break;

  // Gauss-Seidel method
  case Options::ALGORITHM::GS: case Options::ALGORITHM::GSS:
  case Options::ALGORITHM::KRAW: case Options::ALGORITHM::KRAWS:
    throw Exceptions( Exceptions::DAG );
  }

  _final_stats( _stats_ae );
  if( flag != NORMAL ) return flag;

  for( unsigned i=0; i<_ndep; i++ ) X[i] = var[_fpdep[i]];
  if( options.DISPLEVEL >= 1 ){
    std::cout << std::endl << "Symbolic Solution:\n";
    for( unsigned i=0; i<_ndep; i++ )
      os << std::setw(5) << _var[i] << " = " << X[i] << std::endl;
    _output_stats( _stats_ae, os );
  }
  return NORMAL;
}
#endif
} // end namescape mc

#endif

