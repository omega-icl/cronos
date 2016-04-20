// Copyright (C) 2014-2015 Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

#ifndef MC__AEBND_HPP
#define MC__AEBND_HPP

#undef  MC__AEBND_DEBUG
#undef  MC__AEBND_DISABLE_MC21A

#include <stdexcept>
#include <cassert>
#include <cmath>
#include <fstream>
#include <vector>

#include "base_ae.hpp"
#include "mcfunc.hpp"
#include "ffunc.hpp"

#include "ellipsoid.hpp"
#include "tmodel.hpp"
#include "cmodel.hpp"

namespace mc
{
// For block decomposition
extern "C" void mc13d_
  ( const int*, const int*, const int*, const int*, const int*, int*, int*, int*, int* );
extern "C" void mc21a_
  ( const int*, const int*, const int*, const int*, const int*, int*, int*, int* );

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
template <typename T, typename PMT=mc::TModel<T>, typename PVT=mc::TVar<T> >
class AEBND: public virtual BASE_AE
{
private:
  //! @brief number of independent variables
  unsigned _nvar;

  //! @brief number of dependent variables
  unsigned _ndep;

  //! @brief list of operations in each AE block
  std::vector< std::list<const FFOp*> > _opAE;

  //! @brief list of operations in each AE block Jacobian
  std::vector< std::list<const FFOp*> > _opAEJAC;

  //! @brief maximum number of operations among any AE block evaluation
  unsigned _maxopAE;

  //! @brief intermediate operations during AE evaluation in T arithmetic
  std::vector<T> _Iop;

  //! @brief intermediate operations during AE evaluation in PM arithmetic
  std::vector<PVT> _PMop;

  //! @brief AE system entries
  std::vector<FFVar> _pAE;

  //! @brief maximum size among any AE block
  unsigned _maxpAE;

  //! @brief equation indices after possible permutation
  std::vector<unsigned> _bAE;

  //! @brief AE system Jacobian entries
  std::vector<const FFVar*> _pAEJAC;

  //! @brief variables for DAG evaluation
  std::vector<FFVar> _pVAR;

  //! @brief variable indices after possible permutation (forward)
  std::vector<unsigned> _bVAR;

  //! @brief variable indices after possible permutation (reverse)
  std::vector<unsigned> _bVARrev;

  //! @brief variable bounds for DAG evaluation
  std::vector<T> _IVAR;

  //! @brief reference variable bounds for DAG evaluation
  std::vector<T> _IREF;

  //! @brief pointer to dependent interval bounds **DO NOT FREE**
  T *_Ix;

  //! @brief pointer to parameter interval bounds **DO NOT FREE**
  T *_Ip;

  //! @brief function interval bounds
  std::vector<T> _If;

  //! @brief function Jacobian interval bounds
  std::vector<T> _Idfdx;

  //! @brief polynomial model environment
  PMT *_PMenv;

  //! @brief variable bounds for DAG evaluation
  std::vector<PVT> _PMVAR;

  //! @brief variable bounds/values for DAG evaluation
  std::vector<PVT> _PMREF;

  //! @brief pointer to state polynomial models **DO NOT FREE**
  PVT *_PMx;

  //! @brief pointer to parameter polynomial models **DO NOT FREE**
  PVT *_PMp;

  //! @brief function polynomial models
  std::vector<PVT> _PMf;

  //! @brief function Jacobian interval bounds
  std::vector<PVT> _PMdfdx;

  //! @brief Preconditioning matrix
  CPPL::dgematrix _Y;

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
     NORMAL=0,	 //!< Normal execution
     EMPTY,	 //!< Empty solution set
     SINGULAR,	 //!< Structural or numerical singularity encountered
     FAILURE 	 //!< Error encountered in the underlying arithmetic
  };

  //! @brief Integrator options
  struct Options
  {
    //! @brief Constructor
    Options( int iDISP=1 ):
      BOUNDER(AUTO), PRECOND(INVMD), BLKDEC(false), INTERBND(true),
      MAXIT(10), RTOL(1e-7), ATOL(machprec()), DISPLAY(iDISP)
      {}
    //! @brief Enumeration of bounding methods (for polynomial models only)
    enum BOUNDING_METHOD{
      AUTO=0,  //!< Automatic selection
      KRAW=1,  //!< Krawczyk method
      KRAWS=2, //!< Simplified Krawczyk method
      GS=3,    //!< Gauss-Seidel
      GSS=4,   //!< Simplified Gauss-Seidel
      GE=5     //!< Gauss Elimination
    };
    //! @brief Bounding method
    BOUNDING_METHOD BOUNDER;
    //! @brief Enumeration of preconditioning methods
    enum PRECOND_METHOD{
      NONE=0,   //!< No preconditioning used
      INVMD=1,  //!< Inverse of Jacobian mid-point (dense)
      INVMB=2,  //!< Inverse of Jacobian mid-point (banded)
      QRM=3     //!< Unitary Q matrix in QR decomposition of Jacobian mid-point
    };
    //! @brief Bounding method
    PRECOND_METHOD PRECOND;
    //! @brief Whether to apply block decomposition (default: true)
    bool BLKDEC;
    //! @brief Whether to intersect bounds from one iteration to the next (default: true)
    bool INTERBND;
    //! @brief Maximum number of iterations (default: 10, no limit: 0)
    unsigned int MAXIT;
    //! @brief Relative stopping tolerance (default: 1e-7)
    double RTOL;
    //! @brief Absolute stopping tolerance (default: MACHPREC)
    double ATOL;
    //! @brief Display level (default: 1)
    int DISPLAY;
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
      DAG,		//!< DAG may not be obtained for the solution of nonlinear implicit systems or using iterative methods
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
        return "AEBND::Exceptions  Gauss elimination  may not be applied to nonlinear implicit systems";
      case DAG:
        return "AEBND::Exceptions  DAG may not be obtained for the solution of nonlinear implicit systems or using iterative methods";
      case INTERN: default:
        return "AEBND::Exceptions  Internal error";
      }
    }
  private:
    TYPE _ierr;
  };

  //! @brief Perform block decomposition of AE system
  STATUS setup
    ( std::ostream&os=std::cout );
  //! @brief Compute interval enclosure of solution set of parametric AEs
  STATUS solve
    ( const T*Ip, T*Ix, const T*Ix0=0, std::ostream&os=std::cout );
  //! @brief Compute polynomial model of solution set of parametric AEs
  STATUS solve
    ( const PVT*PMp, PVT*PMx, const T*Ix0, std::ostream&os=std::cout );
  //! @brief Compute polynomial model of solution set of parametric AEs
  STATUS solve
    ( const PVT*PMp, PVT*PMx, const PVT*PMx0=0, std::ostream&os=std::cout );
  //! @brief Compute symbolic solution (DAG) of parametric AEs
  STATUS solve
    ( FFVar*X, std::ostream&os=std::cout );
  /** @} */

private:
  //! @brief Flag for setup function
  bool _issetup;
  //! @brief Structural singularity of equation system
  bool _singstruct;
  //! @brief Linearity of problem w.r.t the dependents
  bool _islindep;
  //! @brief Linearity of problem blocks w.r.t. the block variables
  std::vector<bool> _islinblk;
  //! @brief Lower and upper band width of problem blocks
  std::vector< std::pair<long,long> > _bandblk;
  //! @brief Perform block decomposition of AE system
  void _blockdec
    ( std::ostream&os );
  //! @brief Compute preconditionned equation system LHS and RHS for given block
  template <typename U, typename V>
  void _precondlin
    ( const unsigned iblk, const unsigned ndepblk, const unsigned posblk,
      std::vector<U>&A, std::vector<U>&b, U*opf, U*f, U*ref, V*opjacf,
      V*jacf, V*jacvar, std::ostream&os );
  //! @brief Initialize bounding for interval case
  bool _init
    ( const T*Ip );
  //! @brief Initialize bounding for polynomial model case
  bool _init
    ( const PVT*PMp );
  //! @brief Set reference value (mid-point) for interval bounds
  void _reference
    ( const unsigned ndepblk, const unsigned posblk, T*var,
      T*ref, T*jacvar );
  //! @brief Set reference model (centered polynomial) for polynomial models
  void _reference
    ( const unsigned ndepblk, const unsigned posblk, PVT*var,
      PVT*ref, PVT*jacvar );
  //! @brief Set reference model (centered polynomial) for polynomial models
  void _reference
    ( const unsigned ndepblk, const unsigned posblk, PVT*var,
      PVT*ref, T*jacvar );
  //! @brief Test convergence for interval bounds
  bool _cvgtest
    ( const unsigned nx, const T*x, const T*x0 ) const;
  //! @brief Test convergence for polynomial models
  bool _cvgtest
    ( const unsigned nx, const PVT*x, const PVT*x0 ) const;
  //! @brief Apply Gauss-Seidel method (linear system)
  template <typename U, typename V>
  STATUS _gs
    ( U*var, U*opf, U*f, U*ref, V*opjacf, V*jacf, V*jacvar, bool usekraw,
      std::ostream&os );
  //! @brief Apply Gauss-Seidel method to given block (linear system)
  template <typename U, typename V>
  STATUS _gs
    ( const unsigned iblk, const unsigned ndepblk, U*var, U*opf, U*f,
      U*ref, V*opjacf, V*jacf, V*jacvar, bool usekraw, std::ostream&os );
  //! @brief Apply Gauss elimination method (linear system)
  template <typename U>
  STATUS _ge
    ( U*var, U*op, U*f, U*jacf, std::ostream&os );
  //! @brief Apply Gauss elimination method to given block (linear system)
  template <typename U>
  STATUS _ge
    ( const unsigned iblk, const unsigned ndepblk, U*var, U*op, U*f, U*jacf,
      std::ostream&os );
  //! @brief Apply Gauss elimination method for symbolic solution (linear system)
  STATUS _ge
    ( std::vector<FFVar>&X, std::ostream&os );
  //! @brief Index position in 2d-array
  static unsigned _ndx
    ( const unsigned i, const unsigned j, const unsigned n )
    { return i*n+j; } // b/c Jacobian in DAG is transpose
    //{ return i+j*n; }

  //! @brief Structure storing implicit solver statistics
  struct Stats
  {
    //! @brief Constructor
    Stats():
      cputime(0.), maxIter(0)
      {}
    //! @brief Constructor
    void reset()
      { cputime = 0.; maxIter = 0; }

    //! @brief CPU time
    double cputime;
    //! @brief Iteration count
    unsigned long maxIter;
  } _stats_ae;

  //! @brief Function to initialize implicit solver statistics
  static void _init_stats
    ( Stats&stats );
  //! @brief Function to finalize implicit solver statistics
  static void _final_stats
    ( Stats&stats );
  //! @brief Function to display implicit solver statistics
  static void _output_stats
    ( const Stats&stats, std::ostream&os=std::cout );

  //! @brief Private methods to block default compiler methods
  AEBND(const AEBND&);
  AEBND& operator=(const AEBND&);
};

template <typename T, typename PMT, typename PVT> inline
AEBND<T,PMT,PVT>::AEBND
()
: BASE_AE(), _singstruct(false), _islindep(false)
{
  // Initalize state/parameter arrays
  _Ix = _Ip = 0;
  _PMp = _PMx = 0;
  _PMenv = 0;
}

template <typename T, typename PMT, typename PVT>
inline
AEBND<T,PMT,PVT>::~AEBND
()
{
  std::vector<const FFVar*>::iterator it = _pAEJAC.begin();
  for( ; it != _pAEJAC.end(); ++it ) delete[] *it;
  /* DO NOT DELETE _Ix, _Ip */
  /* DO NOT DELETE _PMx, _PMp, _PMenv */
}

template <typename T, typename PMT, typename PVT>
inline void
AEBND<T,PMT,PVT>::_init_stats
( Stats&stats )
{
  // Initialize statistics
  stats.reset();
  timeval time;
  gettimeofday(&time, 0) ;
  stats.cputime = - time.tv_sec - time.tv_usec*1e-6;
}

template <typename T, typename PMT, typename PVT>
inline void
AEBND<T,PMT,PVT>::_final_stats
( Stats&stats )
{
  // Get final CPU time
  timeval time;
  gettimeofday(&time, 0);
  stats.cputime += time.tv_sec + time.tv_usec*1e-6;
}

template <typename T, typename PMT, typename PVT>
inline void
AEBND<T,PMT,PVT>::_output_stats
( const Stats&stats, std::ostream&os )
{
  // Statistics
  os << " No MAX ITERATIONS  " << stats.maxIter
     << std::endl
     << " CPU TIME (SEC)     " << std::fixed << std::left
                               << std::setprecision(5) << stats.cputime
     << std::endl;
  return;
}

template <typename T, typename PMT, typename PVT>
inline typename AEBND<T,PMT,PVT>::STATUS
AEBND<T,PMT,PVT>::setup
( std::ostream&os )
{
  _issetup = false;

  // Check problem size
  _nvar = _var.size();
  _ndep = _dep.size();
  if( !_ndep || _sys.size()!=_ndep ) return FAILURE;

  // (Re)size variable arrays
  const unsigned nVAR = _nvar+_ndep;
  if( _pVAR.size() < nVAR ){
    _pVAR.resize(nVAR);
    _IVAR.clear();  _IREF.clear();
    _PMVAR.clear(); _PMREF.clear(); 
  }

  // Populate variable arrays
  unsigned ivar=0;
  for( unsigned idep=0; idep<_ndep; idep++, ivar++ ){
    _pVAR[ivar] = _dep[idep];
#ifdef MC__AEBND_DEBUG
    os << "dependent #" << idep << ": " << _pVAR[ivar] << std::endl;
#endif
  }
  for( unsigned ipar=0; ipar<_nvar; ipar++, ivar++ ){
    _pVAR[ivar] = _var[ipar];
#ifdef MC__AEBND_DEBUG
    os << "independent #" << ipar << ": " << _pVAR[ivar] << std::endl;
#endif
  }

  // Initialize DAG operations w/ block decomposition
  _pAE.resize(_ndep);
  for( unsigned i=0; i<_ndep; i++ ) _pAE[i] = _sys[i];
  _blockdec( os );
  if( _singstruct ) return SINGULAR;

  _opAE.resize(_bAE.size());
  std::vector<const FFVar*>::iterator it = _pAEJAC.begin();
  for( ; it != _pAEJAC.end(); ++it ) delete[] *it;
  _pAEJAC.clear(); _opAEJAC.clear();
  _islinblk.resize(_bAE.size());
  _bandblk.resize(_bAE.size());

  _maxpAE = _maxopAE = 0;
  for( unsigned iblk=0; iblk<_bAE.size(); iblk++ ){
    // Initialize block operations and Jacobian
#ifdef MC__AEBND_DEBUG
    double blkoper_cpu = -time();
#endif
    const unsigned ndepblk = ( iblk==_bAE.size()-1? _ndep: _bAE[iblk+1] ) - _bAE[iblk]; 
    const unsigned posblk = _ndep - _bAE[iblk] - ndepblk;   
    _opAE[iblk] = _dag->subgraph( ndepblk, _pAE.data()+posblk );

    const FFVar* pJACi = _dag->FAD( ndepblk, _pAE.data()+posblk, ndepblk, _pVAR.data()+posblk ); 
    _pAEJAC.push_back( pJACi );
    std::list<const FFOp*> opJACi = _dag->subgraph( ndepblk*ndepblk, pJACi );  
    _opAEJAC.push_back( opJACi );

    if( _maxpAE < ndepblk ) _maxpAE = ndepblk;
    if( _maxopAE < _opAE[iblk].size() ) _maxopAE = _opAE[iblk].size();
    if( _maxopAE < opJACi.size() ) _maxopAE = opJACi.size();
#ifdef MC__AEBND_DEBUG
    blkoper_cpu += time();
#endif

    // Detect block linearity
#ifdef MC__AEBND_DEBUG
    double blkprop_cpu = -time();
#endif
    std::vector<FFDep> depblk(ndepblk), varblk(_pVAR.size()-posblk);
    for( unsigned i(0); i<ndepblk; i++ ) varblk[i].indep(i);
    _dag->eval( _opAE[iblk], ndepblk, _pAE.data()+posblk, depblk.data(),
                 _pVAR.size()-posblk, _pVAR.data()+posblk, varblk.data() );
    _islinblk[iblk] = true;
    _bandblk[iblk].first = _bandblk[iblk].second = 0;
    for( unsigned i=0; i<ndepblk; i++ ){
      std::map<int,bool>::const_iterator cit = depblk[i].dep().begin();
      for( ; cit != depblk[i].dep().end(); ++cit ){
        if( !(*cit).second ) _islinblk[iblk] = false; // Detecting linearity
        if( _bandblk[iblk].first  < (int)i-(*cit).first ) // Updating lower band width
          _bandblk[iblk].first  = (int)i-(*cit).first;
        if( _bandblk[iblk].second < (*cit).first-(int)i ) // updating upper band width
          _bandblk[iblk].second = (*cit).first-(int)i;
      }
    }
#ifdef MC__AEBND_DEBUG
    blkprop_cpu += time();
    std::cout << "Block #" << iblk << ": "
              << (_islinblk[iblk]?"linear":"nonlinear") << ", bandwidth "
              << _bandblk[iblk].first << "," << _bandblk[iblk].second << std::endl
              << "CPU: " << blkoper_cpu << " (oper)  " << blkprop_cpu << " (prop)\n";
#endif
  }

  _maxopAE*=2;
  _issetup = true;
  return NORMAL;
}

template <typename T, typename PMT, typename PVT>
inline void
AEBND<T,PMT,PVT>::_blockdec
( std::ostream&os )
{
  // Compute structure with mc::FFDep
  const unsigned nVAR = _pVAR.size();
  std::vector<FFDep> VAR(nVAR), SYS(_ndep);
  for( unsigned i(0); i<_ndep; i++ ) VAR[i].indep(i);
  _dag->eval( _ndep, _pAE.data(), SYS.data(), nVAR, _pVAR.data(), VAR.data() );

  // Detect linearity and populate sparse arrays
  std::vector<int> IP(_ndep), LENR(_ndep), ICN;
  _islindep = true;
  for( unsigned i=0; i<_ndep; i++ ){
    IP[i] = ICN.size()+1;
    LENR[i] = SYS[i].dep().size();
    std::map<int,bool>::const_iterator cit = SYS[i].dep().begin();
    for( ; cit != SYS[i].dep().end(); ++cit ){
      ICN.push_back( (*cit).first+1 );
      if( !(*cit).second ) _islindep = false; // Detecting linearity
    }
  }
#ifdef MC__AEBND_DEBUG
  std::cout << "Linearity: " << (_islindep?'Y':'N') << std::endl;
#endif

  // Make a row permutation to remove nonzeros on diagonal: MC21A
  int N = IP.size(), LICN = ICN.size();
  std::vector<int> IPERM(_ndep), IW(4*_ndep);
#ifdef MC__AEBND_DISABLE_MC21A
  for( unsigned i=0; i<_ndep; i++ ) IPERM[i] = i+1;
  _singstruct = false;
#else
  int NUMNZ; 
  mc21a_( &N, ICN.data(), &LICN, IP.data(), LENR.data(), IPERM.data(),
          &NUMNZ, IW.data() );
  _singstruct = NUMNZ<N? true: false;
#ifdef MC__AEBND_DEBUG
  std::cout << "Structural singularity: " << (_singstruct?'Y':'N') << std::endl;
#endif
  // return if structurally singular matrix or block decomposition not requested
  if( _singstruct || !options.BLKDEC ){
    _bAE.resize(1); _bAE[0] = 0;
    _bVAR.resize(_ndep); _bVARrev.resize(_ndep);
    for( unsigned i=0; i<_ndep; i++ ) _bVAR[i] = _bVARrev[i] = i;
    std::vector<FFVar> pAEblk(_ndep);
    for( unsigned i=0; i<_ndep; i++ ) pAEblk[i] = _pAE[IPERM[i]-1];
    _pAE.swap(pAEblk);

    // Display permuted system structure
    if( options.DISPLAY >= 2 ){
      os << std::endl << "Number of Blocks: 1" << std::endl;
      os << "   ";
      for( unsigned j=0; j<_ndep; j++ )
        os << " " << std::setw(3) << j;
      os << std::endl;
      for( unsigned i=0; i<_ndep; i++ ){
        os << std::setw(3) << IPERM[i]-1;
        for( unsigned j=0; j<_ndep; j++ )
          os << std::setw(3) << " " << (SYS[IPERM[i]-1].dep(j).first?"X":" ");
        os << std::endl;
      }
    }

    return;
  }

  // Permute order of equation system using IPERM (!!Fortran style indices!!)
  ICN.clear();
  for( unsigned i=0; i<_ndep; i++ ){
    IP[i] = ICN.size()+1;
    LENR[i] = SYS[IPERM[i]-1].dep().size();
    std::map<int,bool>::const_iterator cit = SYS[IPERM[i]-1].dep().begin();
    for( ; cit != SYS[IPERM[i]-1].dep().end(); ++cit )
      ICN.push_back( (*cit).first+1 );
  }
#ifdef MC__AEBND_DEBUG
  std::cout << "Row reordering: ";
  for( unsigned i=0; i<_ndep; i++) std::cout << " " << IPERM[i];
  std::cout << std::endl;
#endif
#endif

  // Make a block lower-triangular decomposition: MC13D
  int NUM; 
  std::vector<int> IOR(_ndep), IB(_ndep);
  _bAE.resize(_ndep);
  mc13d_( &N, ICN.data(), &LICN, IP.data(), LENR.data(), IOR.data(),
          IB.data(), &NUM, IW.data() );
  _bAE.resize(NUM);
  for( int i=0; i<NUM; i++ ) _bAE[i] = IB[i]-1;
#ifdef MC__AEBND_DEBUG
  std::cout << "Number of blocks in permuted matrix: " << NUM << std::endl;
  std::cout << "Row/Column reordering: ";
  for( unsigned i=0; i<_ndep; i++) std::cout << " " << IOR[i];
  std::cout << std::endl;
#endif

  // Permute order of equation system AND variables using IOR
  // and keep track of inverse permutation
  std::vector<FFVar> pAEblk(_ndep), pVARblk(_pVAR);
  _bVAR.resize(_ndep); _bVARrev.resize(_ndep);
  for( unsigned i=0; i<_ndep; i++ ){
    pAEblk[i]  = _pAE[IPERM[IOR[_ndep-i-1]-1]-1];
    pVARblk[i] = _pVAR[IOR[_ndep-i-1]-1];
    _bVAR[IOR[i]-1] = _ndep-i-1;
    _bVARrev[_ndep-i-1] = IOR[i]-1;
  }
  //for( int k=0; k<NUM; k++ ){
  //  const unsigned iL = IB[0]-1, iU = ( k==NUM-1? _ndep: IB[k+1]-1 );
  //  for( unsigned i=iL; i<iU; i++ ){
  //    pVARblk[i] = _pVAR[IOR[_ndep-iU+(i-iL)]-1];
  //    _bVAR[IOR[i]-1] = _ndep-iU+(i-iL);
  //  }
  //}
  _pAE.swap(pAEblk);
  _pVAR.swap(pVARblk);

  // Display permuted system structure
  if( options.DISPLAY >= 2 ){
    std::cout << std::endl << "Number of Blocks: " << NUM << std::endl;
    os << "   ";
    for( unsigned j=0; j<_ndep; j++ )
      os << " " << std::setw(3) << IOR[_ndep-j-1]-1;
    os << std::endl;
    for( unsigned i=0; i<_ndep; i++ ){
      os << std::setw(3) << IPERM[IOR[_ndep-i-1]-1]-1;
      for( unsigned j=0; j<_ndep; j++ )
        os << std::setw(3) << " "
           << (SYS[IPERM[IOR[_ndep-i-1]-1]-1].dep(IOR[_ndep-j-1]-1).first?"X":" ");
      os << std::endl;
    }
  }
}

template <typename T, typename PMT, typename PVT>
inline bool
AEBND<T,PMT,PVT>::_init
( const T*Ip )
{
  if( !_issetup ) return false;

  // (Re)size variable arrays
  const unsigned nVAR = _nvar+_ndep;
  if( !_IVAR.size() ) _IVAR.resize(nVAR);
  if( !_IREF.size() ) _IREF.resize(nVAR);
  _Ix = _IVAR.data(); _Ip = _IVAR.data() + _ndep;
  _PMx = _PMp = 0;
  _If.resize(_maxpAE);
  _Idfdx.resize(_maxpAE*_maxpAE);
  _PMf.clear();
  _PMdfdx.clear();
  _PMenv = 0;

  // Populate parameters in variable arrays
  for( unsigned ipar=0; ipar<_nvar; ipar++ )
    _IVAR[_ndep+ipar] = _IREF[_ndep+ipar] = Ip[ipar];

  // (Re)size DAG evaluation arrays
  _Iop.resize( _maxopAE );
  _PMop.clear();

  return _issetup;
}

template <typename T, typename PMT, typename PVT>
inline bool
AEBND<T,PMT,PVT>::_init
( const PVT*PMp )
{
  if( !_issetup ) return false;

  // (Re)size variable arrays
  const unsigned nVAR = _nvar+_ndep;
  if( !_IVAR.size() )  _IVAR.resize(nVAR);
  if( !_PMVAR.size() ) _PMVAR.resize(nVAR);
  if( !_PMREF.size() ) _PMREF.resize(nVAR);
  _Ix = _IVAR.data(); _Ip = _IVAR.data() + _ndep;
  _Idfdx.resize(_maxpAE*_maxpAE);
  _PMx = _PMVAR.data(); _PMp = _PMVAR.data() + _ndep;
  _PMf.resize(_maxpAE);
  _PMdfdx.resize(_maxpAE*_maxpAE);
  _PMenv = 0;
  for( unsigned i=0; i<_nvar; i++ )
    if( PMp[i].env() ) _PMenv = PMp[i].env();

  // Populate parameters in variable arrays
  for( unsigned ipar=0; ipar<_nvar; ipar++ ){
    _IVAR[_ndep+ipar] = PMp[ipar].bound();
    _PMVAR[_ndep+ipar] = _PMREF[_ndep+ipar] = PMp[ipar];
  }

  // (Re)size DAG evaluation arrays
  _Iop.resize( _maxopAE );
  _PMop.resize( _maxopAE );

  return _issetup;
}

template <typename T, typename PMT, typename PVT>
inline bool
AEBND<T,PMT,PVT>::_cvgtest
( const unsigned nx, const T*x, const T*x0 ) const
{
  // Relative tolerance uses the largest improvement divided by interval width in any direction
  // Absolute tolerance uses the largest improvement in any direction
  double rtol = 0., atol = 0.;
  for( unsigned i=0; i<nx; i++ ){	
    double diamxi = Op<T>::diam(x[i]);
    if( diamxi > 0. ){
        rtol = std::max(rtol, (Op<T>::l(x[i]) - Op<T>::l(x0[i]))/diamxi);
        rtol = std::max(rtol, (Op<T>::u(x0[i]) - Op<T>::u(x[i]))/diamxi);
    atol = std::max(atol, (Op<T>::l(x[i]) - Op<T>::l(x0[i])));
    atol = std::max(atol, (Op<T>::u(x0[i]) - Op<T>::u(x[i])));
    }
  }
#ifdef MC__AEBND_DEBUG
  std::cout << "RTOL =" << rtol  << "  ATOL =" << atol << std::endl;
#endif
  return rtol <= options.RTOL || atol <= options.ATOL? true: false;
}

template <typename T, typename PMT, typename PVT>
inline bool
AEBND<T,PMT,PVT>::_cvgtest
( const unsigned nx, const PVT*x, const PVT*x0 ) const
{
  // Relative tolerance uses the largest improvement divided by interval width in any direction
  // Absolute tolerance uses the largest improvement in any direction
  double rtol = 0., atol = 0.;
  for( unsigned i=0; i<nx; i++ ){	
    double diamxi = Op<T>::diam(x[i].R());
    if( diamxi > 0. ){
        rtol = std::max(rtol, (Op<T>::l(x[i].R()) - Op<T>::l(x0[i].R()))/diamxi);
        rtol = std::max(rtol, (Op<T>::u(x0[i].R()) - Op<T>::u(x[i].R()))/diamxi);
    atol = std::max(atol, (Op<T>::l(x[i].R()) - Op<T>::l(x0[i].R())));
    atol = std::max(atol, (Op<T>::u(x0[i].R()) - Op<T>::u(x[i].R())));
    }
  }
#ifdef MC__AEBND_DEBUG
  std::cout << "RTOL =" << rtol  << "  ATOL =" << atol << std::endl;
#endif
  return rtol <= options.RTOL || atol <= options.ATOL? true: false;
}

template <typename T, typename PMT, typename PVT>
inline void
AEBND<T,PMT,PVT>::_reference
( const unsigned ndepblk, const unsigned posblk, T*var,
  T*ref, T*jacvar )
{
  unsigned ivar = posblk;
  for( ; ivar<posblk+ndepblk; ivar++ )
    ref[ivar] = Op<T>::mid( var[ivar] );
  for( ; ivar<_ndep; ivar++ )
    ref[ivar] = var[ivar];
  if( var != jacvar ) throw Exceptions( Exceptions::INTERN );
}

template <typename T, typename PMT, typename PVT>
inline void
AEBND<T,PMT,PVT>::_reference
( const unsigned ndepblk, const unsigned posblk, PVT*var,
  PVT*ref, PVT*jacvar )
{
  unsigned ivar = posblk;
  for( ; ivar<posblk+ndepblk; ivar++ )
    ref[ivar] = var[ivar].C().P();
  for( ; ivar<_ndep; ivar++ )
    ref[ivar] = var[ivar];
  if( var != jacvar ) throw Exceptions( Exceptions::INTERN );
}

template <typename T, typename PMT, typename PVT>
inline void
AEBND<T,PMT,PVT>::_reference
( const unsigned ndepblk, const unsigned posblk, PVT*var,
  PVT*ref, T*jacvar )
{
  unsigned ivar = posblk;
  for( ; ivar<posblk+ndepblk; ivar++ ){
    ref[ivar] = var[ivar].center().polynomial();
    jacvar[ivar] = var[ivar].B();
  }
  for( ; ivar<_ndep; ivar++ ){
    ref[ivar] = var[ivar];
    jacvar[ivar] = var[ivar].B();
  }
}

template <typename T, typename PMT, typename PVT> 
template <typename U, typename V>
inline void
AEBND<T,PMT,PVT>::_precondlin
( const unsigned iblk, const unsigned ndepblk, const unsigned posblk,
  std::vector<U>&A, std::vector<U>&b, U*opf, U*f, U*ref, V*opjacf,
  V*jacf, V*jacvar, std::ostream&os )
{
  try{
    if( jacf ){
      // Get Jacobian and preconditioning matrix
      _dag->eval( _opAEJAC[iblk], opjacf, ndepblk*ndepblk, _pAEJAC[iblk], jacf,
                   _pVAR.size()-posblk, _pVAR.data()+posblk, jacvar+posblk );
#ifdef MC__AEBND_DEBUG
      std::cout << "Non-Preconditioned LHS Block #" << iblk+1 << ":\n";
      for( unsigned i=0; i<ndepblk; i++ ){
        for( unsigned j=0; j<ndepblk; j++ )
          std::cout << jacf[_ndx(i,j,ndepblk)] << "  ";
        std::cout << std::endl;
      }
      std::cout << std::endl;
#endif
      CPPL::dgematrix mA, Q, R;
      CPPL::dgbmatrix mB;
      switch( options.PRECOND ){
      case Options::INVMD:
        if( ndepblk==1 ){ _Y.resize(ndepblk,ndepblk); _Y.identity(); break; }
        mA.resize(ndepblk,ndepblk);
        for( unsigned i=0; i<ndepblk; i++ )
          for( unsigned j=0; j<ndepblk; j++ )
            mA(i,j) = Op<U>::mid( jacf[_ndx(i,j,ndepblk)] );
#ifdef MC__AEBND_DEBUG
        std::cout << "Preconditioning matrix Block #" << iblk+1 << ":\n" << mA << std::endl;
#endif
        if( dgesv( mA, _Y ) ) throw Exceptions( Exceptions::PRECOND );
        break;
      case Options::INVMB:
        if( ndepblk==1 ){ _Y.resize(ndepblk,ndepblk); _Y.identity(); break; }
        mB.resize(ndepblk,ndepblk,_bandblk[iblk].first,_bandblk[iblk].second);
        for( unsigned i=0; i<ndepblk; i++ )
          for( unsigned j=0; j<ndepblk; j++ ){
            if( _bandblk[iblk].first  < (int)i-(int)j
             || _bandblk[iblk].second < (int)j-(int)i ) continue;
            mB(i,j) = Op<U>::mid( jacf[_ndx(i,j,ndepblk)] );
          }
#ifdef MC__AEBND_DEBUG
        std::cout << "Preconditioning matrix Block #" << iblk+1 << ":\n" << mA << std::endl;
#endif
        if( dgbsv( mB, _Y ) ) throw Exceptions( Exceptions::PRECOND );
        break;
      case Options::QRM:
        if( ndepblk==1 ){ _Y.resize(ndepblk,ndepblk); _Y.identity(); break; }
        mA.resize(ndepblk,ndepblk);  
        for( unsigned i=0; i<ndepblk; i++ )
          for( unsigned j=0; j<ndepblk; j++ )
            mA(j,i) = Op<U>::mid( jacf[_ndx(i,j,ndepblk)] );
#ifdef MC__AEBND_DEBUG
        std::cout << "Preconditioning matrix Block #" << iblk+1 << ":\n" << mA << std::endl;
#endif
        if( dgeqrf( mA, Q, R ) ) throw Exceptions( Exceptions::PRECOND );
        _Y = Q;
        break;
      case Options::NONE:
        _Y.resize(ndepblk,ndepblk); _Y.identity();
        break;
      }
#ifdef MC__AEBND_DEBUG
      std::cout << "Preconditioning matrix Block #" << iblk+1 << ":\n" << _Y << std::endl;
#endif

      //setting A = Y*dfdx
      for( unsigned i=0; i<ndepblk; i++ )
        for( unsigned j=0; j<ndepblk; j++ ){
          A[_ndx(i,j,ndepblk)] = 0.;
          for( unsigned k=0; k<ndepblk; k++ )
            A[_ndx(i,j,ndepblk)] += _Y(i,k) * jacf[_ndx(k,j,ndepblk)];        
        }
#ifdef MC__AEBND_DEBUG
      std::cout << "Preconditioned LHS Block #" << iblk+1 << ":\n";
      for( unsigned i=0; i<ndepblk; i++ ){
        for( unsigned j=0; j<ndepblk; j++ )
          std::cout << A[_ndx(i,j,ndepblk)] << "  ";
        std::cout << std::endl;
      }
      std::cout << std::endl;
#endif
    }

    if( f ){
      //setting b = Y*f(x)
      _dag->eval( _opAE[iblk], opf, ndepblk, _pAE.data()+posblk, f,
                   _pVAR.size()-posblk, _pVAR.data()+posblk, ref+posblk );
#ifdef MC__AEBND_DEBUG
      _dag->output( _opAE[iblk] );
      std::cout << "Intermediates in Block #" << iblk+1 << ":\n";
      for (unsigned i=0; i<_opAE[iblk].size(); i++ ) std::cout << opf[i] << std::endl;
      std::cout << "Reference in Block #" << iblk+1 << ":\n";
      for (unsigned i=0; i<ndepblk; i++ ) std::cout << ref[posblk+i] << std::endl;
      std::cout << "Non-Preconditioned RHS Block #" << iblk+1 << ":\n";
      for (unsigned i=0; i<ndepblk; i++ ) std::cout << -f[i] << std::endl;
      std::cout << std::endl;
      { int dum; std::cin >> dum; }
#endif
      for (unsigned i=0; i<ndepblk; i++ ){
        b[i] = 0.;
        for( unsigned j=0; j<ndepblk; j++)
          b[i] -= _Y(i,j) * f[j];
      }
#ifdef MC__AEBND_DEBUG
      std::cout << "Preconditioned RHS Block #" << iblk+1 << ":\n";
      for (unsigned i=0; i<ndepblk; i++ ) std::cout << b[i] << std::endl;
      std::cout << std::endl;
#endif
    }
  }
  catch(...){
    throw Exceptions( Exceptions::PRECOND );
  }
}

template <typename T, typename PMT, typename PVT>
template <typename U, typename V>
inline typename AEBND<T,PMT,PVT>::STATUS
AEBND<T,PMT,PVT>::_gs
( U*var, U*opf, U*f, U*ref, V*opjacf, V*jacf, V*jacvar, bool usekraw, std::ostream&os )
{
  for( unsigned iblk=0; iblk<_bAE.size(); iblk++ ){
    const unsigned ndepblk = ( iblk==_bAE.size()-1? _ndep: _bAE[iblk+1] ) - _bAE[iblk];
    STATUS flag = _gs( iblk, ndepblk, var, opf, f, ref, opjacf, jacf, jacvar, usekraw, os );
    if( flag != NORMAL ) return flag;
  }
  return NORMAL;
}

template <typename T, typename PMT, typename PVT>
template <typename U, typename V>
inline typename AEBND<T,PMT,PVT>::STATUS
AEBND<T,PMT,PVT>::_gs
( const unsigned iblk, const unsigned ndepblk, U*var, U*opf, U*f, U*ref,
  V*opjacf, V*jacf, V*jacvar, bool usekraw, std::ostream&os )
{
  const unsigned posblk = _ndep - _bAE[iblk] - ndepblk;
  U*varblk  = var + posblk;
  U*refblk = ref + posblk;

#ifdef MC__AEBND_DEBUG
  std::cout << "Initial bounds:\n";
  for( unsigned i=0; i<ndepblk; i++ ) std::cout << varblk[i] << std::endl;
#endif

  std::vector<U> G(ndepblk*ndepblk), b(ndepblk), varblk0; 
  for( unsigned iter=0; iter<options.MAXIT; iter++ ){
    if( _stats_ae.maxIter <= iter ) _stats_ae.maxIter = iter+1;

    // Update reference
    _reference( ndepblk, posblk, var, ref, jacvar );

    // (Re)compute preconditionned LHS matrix (only if needed) and RHS vector
    _precondlin( iblk, ndepblk, posblk, G, b, opf, f, ref,
                 opjacf, (!iter||!_islinblk[iblk]?jacf:0), jacvar, os );

    // Keep track of current iterate
    varblk0.assign( varblk, varblk+ndepblk );

    try{
      // Apply componentwise Krawczyk step
      if( usekraw ){
        for( unsigned i=0; i<ndepblk; i++ ){
          U Xk, temp(0.); 
          for( unsigned j=0; j<ndepblk; j++ ){
            if( j != i )  temp -= G[_ndx(i,j,ndepblk)] * ( varblk[j] - refblk[j] );
            else temp += ( 1. - G[_ndx(i,i,ndepblk)] ) * ( varblk[j] - refblk[j] );
          }
          Xk = refblk[i] + b[i] + temp;
          if( options.INTERBND ) varblk[i] = Xk;
          else if( !Op<U>::inter( varblk[i], Xk, varblk[i] ) ) return EMPTY;
          //if( !Op<U>::inter( varblk[i], Xk, varblk[i] ) ) return EMPTY;
        }
      }

      // Apply componentwise Gauss-Seidel step
      else{
        for( unsigned i=0; i<ndepblk; i++ ){
          U Xk, temp(0.); 
          if( Op<U>::l(G[_ndx(i,i,ndepblk)]) > 0.
           || Op<U>::u(G[_ndx(i,i,ndepblk)]) < 0. ){
            for( unsigned j=0; j<ndepblk; j++ )
              if( j != i ) temp += G[_ndx(i,j,ndepblk)] * ( varblk[j] - refblk[j] );
            Xk = refblk[i] + ( b[i] - temp ) / G[_ndx(i,i,ndepblk)];
          }
          else{
            for( unsigned j=0; j<ndepblk; j++ ){
              if( j != i )  temp -= G[_ndx(i,j,ndepblk)] * ( varblk[j] - refblk[j] );
              else temp += ( 1. - G[_ndx(i,i,ndepblk)] ) * ( varblk[j] - refblk[j] );
            }
            Xk = refblk[i] + b[i] + temp;
          }
          if( options.INTERBND ) varblk[i] = Xk;
          else if( !Op<U>::inter( varblk[i], Xk, varblk[i] ) ) return EMPTY;
          //if( !Op<U>::inter( varblk[i], Xk, varblk[i] ) ) return EMPTY;
          //else if( options.DISPLAY >= 2 )
          //  os << " Skipped Iteration #" << iter+1
          //     << " for Variable X" << _bVARrev[posblk+i] << std::endl;
        }
      }
#ifdef MC__AEBND_DEBUG
      std::cout << "Iter #" << iter << " bounds:\n";
      for( unsigned i=0; i<ndepblk; i++ ) std::cout << varblk[i] << std::endl;
#endif

      // Display
      if( options.DISPLAY >= 2  ){
        os << std::endl << "Block #" << iblk+1 << "  Iteration #" << iter+1 << ":\n";
        for( unsigned i=0; i<ndepblk; i++ )
          //os << " X" << _bVAR[_ndep-posblk-i-1] << ": " << varblk[i] << std::endl;
          os << " X" << _bVARrev[posblk+i] << ": " << varblk[i] << std::endl;
      }

      // Check convergence
      if( _cvgtest( ndepblk, varblk, varblk0.data() ) ) break;

    }
    catch(...){ return FAILURE; }
  }

  return NORMAL;
}

template <typename T, typename PMT, typename PVT>
template <typename U>
inline typename AEBND<T,PMT,PVT>::STATUS
AEBND<T,PMT,PVT>::_ge
( U*var, U*op, U*f, U*jacf, std::ostream&os )
{
  for( unsigned iblk=0; iblk<_bAE.size(); iblk++ ){
    const unsigned ndepblk = ( iblk==_bAE.size()-1? _ndep: _bAE[iblk+1] ) - _bAE[iblk];
    STATUS flag = _ge( iblk, ndepblk, var, op, f, jacf, os );
    if( flag != NORMAL ) return flag;
  }
  return NORMAL;
}

template <typename T, typename PMT, typename PVT>
template <typename U>
inline typename AEBND<T,PMT,PVT>::STATUS
AEBND<T,PMT,PVT>::_ge
( const unsigned iblk, const unsigned ndepblk, U*var, U*op, U*f, U*jacf,
  std::ostream&os )
{
  if( !_stats_ae.maxIter ) _stats_ae.maxIter = 1;
  const unsigned posblk = _ndep - _bAE[iblk] - ndepblk;
  U*varblk = var + posblk;

  // Compute preconditionned LHS matrix and RHS vector
  std::vector<U> G(ndepblk*ndepblk), b(ndepblk); 
  for( unsigned i=0; i<ndepblk; i++ ) varblk[i] = 0.;
  _precondlin( iblk, ndepblk, posblk, G, b, op, f, var,
               op, jacf, var, os );

  try{
    // Perform forward substitution
    for( unsigned k=0; k<ndepblk-1; k++ ){
      if( Op<U>::l( G[_ndx(k,k,ndepblk)] ) <= 0. 
       && Op<U>::u( G[_ndx(k,k,ndepblk)] ) >= 0. ) return SINGULAR;

      for( unsigned i=k+1; i<ndepblk; i++ ){
        U factor = G[_ndx(i,k,ndepblk)] / G[_ndx(k,k,ndepblk)];
        for( unsigned j=k+1; j<ndepblk; j++ )
          G[_ndx(i,j,ndepblk)] -= factor * G[_ndx(k,j,ndepblk)];
        b[i] -= factor * b[k];
      }
    }

    // Perform backward elimination
//std::cout << "G[" << _ndx(ndepblk-1,ndepblk-1,ndepblk) << "] = " << G[_ndx(ndepblk-1,ndepblk-1,ndepblk)] << std::endl;
    if( Op<U>::l( G[_ndx(ndepblk-1,ndepblk-1,ndepblk)] ) <= 0. 
     && Op<U>::u( G[_ndx(ndepblk-1,ndepblk-1,ndepblk)] ) >= 0. )
      return SINGULAR;
//std::cout << "NON-SINGULAR!" << std::endl;
    varblk[ndepblk-1] = b[ndepblk-1] / G[_ndx(ndepblk-1,ndepblk-1,ndepblk)];
    for( unsigned i=ndepblk; i>0; i-- ){
      U sum = b[i-1];
      for( unsigned j=i; j<ndepblk; j++ )
        sum -= G[_ndx(i-1,j,ndepblk)] * varblk[j];
      varblk[i-1] = sum / G[_ndx(i-1,i-1,ndepblk)];
    }
#ifdef MC__AEBND_DEBUG
    std::cout << "GE bounds:\n";
    for( unsigned i=0; i<ndepblk; i++ ) std::cout << varblk[i] << std::endl;
#endif
  }
  catch(...){ return FAILURE; }
  return NORMAL;
}

//! @fn template <typename T, typename PMT, typename PVT> inline typename AEBND<T,PMT,PVT>::STATUS AEBND<T,PMT,PVT>::solve(
//! const T*Ip, T*Ix, const T*Ix0=0, std::ostream&os=std::cout )
//!
//! This function computes an interval enclosure of the solution set of
//! the parametric AEs:
//!   - <a>Ip</a> [input] interval parameter set
//!   - <a>Ix</a> [output] interval state enclosure
//!   - <a>Ix0</a> [input] a priori interval state bounds (default: NULL)
//!   - <a>os</a> [input] output stream (default: std::cout)
//! .
//! The return value is the status.
template <typename T, typename PMT, typename PVT>
inline typename AEBND<T,PMT,PVT>::STATUS
AEBND<T,PMT,PVT>::solve
( const T*Ip, T*Ix, const T*Ix0, std::ostream&os )
{
  _init_stats( _stats_ae );
  if( _singstruct )  return SINGULAR;
  if( !_init( Ip ) ) return FAILURE;

  // Loop over each block
  STATUS flag = NORMAL;
  for( unsigned iblk=0; iblk<_bAE.size() && flag == NORMAL; iblk++ ){

    const unsigned ndepblk = ( iblk==_bAE.size()-1? _ndep: _bAE[iblk+1] ) - _bAE[iblk];
    const unsigned posblk = _ndep - _bAE[iblk] - ndepblk;   

    switch( options.BOUNDER ){

    // Automatic selection of bounder
    case Options::AUTO:
      // No a priori box supplied
      if( !Ix0 ){
        if( !_islindep ) throw Exceptions( Exceptions::APRIORI );
        flag = _ge( iblk, ndepblk, _IVAR.data(), _Iop.data(), _If.data(), _Idfdx.data(), os );
      }
  
      // A priori box supplied and linear case
      else if( _islindep ){
        flag = _ge( iblk, ndepblk, _IVAR.data(), _Iop.data(), _If.data(), _Idfdx.data(), os );
        if( flag != NORMAL ) return flag;

        // check whether supplied bounds are better than GE bounds in parts
        bool improvecheck = false;
        for( unsigned i=0; i<ndepblk; i++ ){
          if( Op<T>::l(_IVAR[posblk+i]) < Op<T>::l(Ix0[_bVARrev[posblk+i]])
           || Op<T>::u(_IVAR[posblk+i]) > Op<T>::u(Ix0[_bVARrev[posblk+i]]) ){
            if( !Op<T>::inter( _IVAR[posblk+i], _IVAR[posblk+i], Ix0[_bVARrev[posblk+i]] ) )
              return EMPTY;
            improvecheck = true;
          }
        }
        if( improvecheck )
          flag = _gs( iblk, ndepblk, _IVAR.data(), _Iop.data(), _If.data(), _IREF.data(),
                      _Iop.data(), _Idfdx.data(), _IVAR.data(), false, os );
      }

      // A priori box supplied and nonlinear case
      else{
        for( unsigned i=0; i<ndepblk; i++ ) _IVAR[posblk+i] = Ix0[_bVARrev[posblk+i]];
        flag = _gs( iblk, ndepblk, _IVAR.data(), _Iop.data(), _If.data(), _IREF.data(),
                    _Iop.data(), _Idfdx.data(), _IVAR.data(), false, os );
      }
      break;

    // Gauss Elimination method
    case Options::GE:
      // Check system is linear in dependents
      if( !_islindep ) throw Exceptions( Exceptions::GAUSSEL );
      flag = _ge( iblk, ndepblk, _IVAR.data(), _Iop.data(), _If.data(), _Idfdx.data(), os );
      break;

    // Gauss-Seidel method
    case Options::GS: case Options::GSS:
      // Check initial box is given
      if( !Ix0 ) throw Exceptions( Exceptions::APRIORI );
      for( unsigned i=0; i<ndepblk; i++ ) _IVAR[posblk+i] = Ix0[_bVARrev[posblk+i]];
      flag = _gs( iblk, ndepblk, _IVAR.data(), _Iop.data(), _If.data(), _IREF.data(),
                  _Iop.data(), _Idfdx.data(), _IVAR.data(), false, os );
      break;

    // Krawczyk method
    case Options::KRAW: case Options::KRAWS:
      // Check initial box is given
      if( !Ix0 ) throw Exceptions( Exceptions::APRIORI );
      for( unsigned i=0; i<ndepblk; i++ ) _IVAR[posblk+i] = Ix0[_bVARrev[posblk+i]];
      flag = _gs( iblk, ndepblk, _IVAR.data(), _Iop.data(), _If.data(), _IREF.data(),
                  _Iop.data(), _Idfdx.data(), _IVAR.data(), true, os );
      break;
    }
  }

  _final_stats( _stats_ae );
  if( flag != NORMAL ) return flag;

  for( unsigned i=0; i<_ndep; i++ ) Ix[i] = _IVAR[_bVAR[i]];
  if( options.DISPLAY >= 1 ){
    std::cout << std::endl << "Computed Bounds:\n";
    //for( unsigned i=0; i<_ndep; i++ )
    //  os << " X" << i << ": " << _IVAR[i] << std::endl;
    //for( unsigned i=0; i<_ndep; i++ )
    //  os << " b" << i << ": " << _bVAR[i] << std::endl;
    for( unsigned i=0; i<_ndep; i++ )
      os << " X" << i << ": " << Ix[i] << std::endl;
    _output_stats( _stats_ae, os );
  }
  return NORMAL;
}

//! @fn template <typename T, typename PMT, typename PVT> inline typename AEBND<T,PMT,PVT>::STATUS AEBND<T,PMT,PVT>::solve(
//! const PVT*PMp, PVT*PMx, const T*Ix0=0, std::ostream&os=std::cout )
//!
//! This function computes an interval enclosure of the solution set of
//! the parametric AEs:
//!   - <a>PMp</a> [input] polynomial model of parameter set
//!   - <a>PMx</a> [output] polynomial model of state enclosure
//!   - <a>Ix0</a> [input] a priori interval state bounds (default: NULL)
//!   - <a>os</a> [input] output stream (default: std::cout)
//! .
//! The return value is the status.
template <typename T, typename PMT, typename PVT>
inline typename AEBND<T,PMT,PVT>::STATUS
AEBND<T,PMT,PVT>::solve
( const PVT*PMp, PVT*PMx, const T*Ix0, std::ostream&os )
{
  if( Ix0 ){
    std::vector<PVT> PMx0(_ndep);
    for( unsigned i=0; i<_ndep; i++ ) PMx0[i] = Ix0[i];
    return solve( PMp, PMx, PMx0.data(), os );
  }
  return solve( PMp, PMx, (const PVT*)0, os );
}

//! @fn template <typename T, typename PMT, typename PVT> inline typename AEBND<T,PMT,PVT>::STATUS AEBND<T,PMT,PVT>::solve(
//! const PVT*PMp, PVT*PMx, const PVT*PMx0=0, std::ostream&os=std::cout )
//!
//! This function computes an interval enclosure of the solution set of
//! the parametric AEs:
//!   - <a>PMp</a> [input] polynomial model of parameter set
//!   - <a>PMx</a> [output] polynomial model of state enclosure
//!   - <a>PMx0</a> [input] a priori polynomial model of state bounds (default: NULL)
//!   - <a>os</a> [input] output stream (default: std::cout)
//! .
//! The return value is the status.
template <typename T, typename PMT, typename PVT>
inline typename AEBND<T,PMT,PVT>::STATUS
AEBND<T,PMT,PVT>::solve
( const PVT*PMp, PVT*PMx, const PVT*PMx0, std::ostream&os )
{
  _init_stats( _stats_ae );
  if( _singstruct )   return SINGULAR;
  if( !_init( PMp ) ) return FAILURE;

  // Loop over each block
  STATUS flag = NORMAL;
  for( unsigned iblk=0; iblk<_bAE.size() && flag == NORMAL; iblk++ ){

    const unsigned ndepblk = ( iblk==_bAE.size()-1? _ndep: _bAE[iblk+1] ) - _bAE[iblk];
    const unsigned posblk = _ndep - _bAE[iblk] - ndepblk;   

    switch( options.BOUNDER ){

    // Automatic selection of bounder
    case Options::AUTO:
      // No a priori box supplied
      if( !PMx0 ){
        if( !_islindep ) throw Exceptions( Exceptions::APRIORI );
        flag = _ge( iblk, ndepblk, _PMVAR.data(), _PMop.data(), _PMf.data(), _PMdfdx.data(), os );
      }
  
      // A priori box supplied and linear case
      else if( _islindep ){
        flag = _ge( iblk, ndepblk, _PMVAR.data(), _PMop.data(), _PMf.data(), _PMdfdx.data(), os );
        if( flag != NORMAL ) return flag;

        // check whether supplied bounds are better than GE bounds in parts
        bool improvecheck = false;
        for( unsigned i=0; i<ndepblk; i++ ){
          if( Op<PVT>::l(_PMVAR[posblk+i]) < Op<PVT>::l(PMx0[_bVARrev[posblk+i]])
           || Op<PVT>::u(_PMVAR[posblk+i]) > Op<PVT>::u(PMx0[_bVARrev[posblk+i]]) ){
            if( !Op<PVT>::inter( _PMVAR[posblk+i], _PMVAR[posblk+i], PMx0[_bVARrev[posblk+i]] ) )
              return EMPTY;
            improvecheck = true;
          }
        }
        if( improvecheck )
          flag = _gs( iblk, ndepblk, _PMVAR.data(), _PMop.data(), _PMf.data(), _PMREF.data(),
                      _Iop.data(), _Idfdx.data(), _IVAR.data(), false, os );
      }

      // A priori box supplied and nonlinear case
      else{
        for( unsigned i=0; i<ndepblk; i++ ) _PMVAR[posblk+i] = PMx0[_bVARrev[posblk+i]];
        flag = _gs( iblk, ndepblk, _PMVAR.data(), _PMop.data(), _PMf.data(), _PMREF.data(),
                    _Iop.data(), _Idfdx.data(), _IVAR.data(), false, os );
      }
      break;

    // Gauss Elimination method
    case Options::GE:
      // Check system is linear in dependents
      if( !_islindep ) throw Exceptions( Exceptions::GAUSSEL );
      flag = _ge( iblk, ndepblk, _PMVAR.data(), _PMop.data(), _PMf.data(), _PMdfdx.data(), os );
      break;

    // Gauss-Seidel method
    case Options::GS:
      // Check initial box is given
      if( !PMx0 ) throw Exceptions( Exceptions::APRIORI );
      for( unsigned i=0; i<ndepblk; i++ ) _PMVAR[posblk+i] = PMx0[_bVARrev[posblk+i]];
      flag = _gs( iblk, ndepblk, _PMVAR.data(), _PMop.data(), _PMf.data(), _PMREF.data(),
                  _PMop.data(), _PMdfdx.data(), _PMVAR.data(), false, os );
      break;

    // Simplified Gauss-Seidel method
    case Options::GSS:
      // Check initial box is given
      if( !PMx0 ) throw Exceptions( Exceptions::APRIORI );
      for( unsigned i=0; i<ndepblk; i++ ) _PMVAR[posblk+i] = PMx0[_bVARrev[posblk+i]];
      flag = _gs( iblk, ndepblk, _PMVAR.data(), _PMop.data(), _PMf.data(), _PMREF.data(),
                  _Iop.data(), _Idfdx.data(), _IVAR.data(), false, os );
      break;

    // Krawczyk method
    case Options::KRAW:
      // Check initial box is given
      if( !PMx0 ) throw Exceptions( Exceptions::APRIORI );
      for( unsigned i=0; i<ndepblk; i++ ) _PMVAR[posblk+i] = PMx0[_bVARrev[posblk+i]];
      flag = _gs( iblk, ndepblk, _PMVAR.data(), _PMop.data(), _PMf.data(), _PMREF.data(),
                  _PMop.data(), _PMdfdx.data(), _PMVAR.data(), true, os );
      break;

    // Simplified Krawczyk method
    case Options::KRAWS:
      // Check initial box is given
      if( !PMx0 ) throw Exceptions( Exceptions::APRIORI );
      for( unsigned i=0; i<ndepblk; i++ ) _PMVAR[posblk+i] = PMx0[_bVARrev[posblk+i]];
      flag = _gs( iblk, ndepblk, _PMVAR.data(), _PMop.data(), _PMf.data(), _PMREF.data(),
                  _Iop.data(), _Idfdx.data(), _IVAR.data(), true, os );
      break;
    }
  }

  _final_stats( _stats_ae );
  if( flag != NORMAL) return flag;

  for( unsigned i=0; i<_ndep; i++ ) PMx[i] = _PMVAR[_bVAR[i]];
  if( options.DISPLAY >= 1 ){
    std::cout << std::endl << "Computed Bounds:\n";
    //for( unsigned i=0; i<_ndep; i++ )
    //  os << " X" << i << ": " << _PMVAR[i] << std::endl;
    //for( unsigned i=0; i<_ndep; i++ )
    //  os << " b" << i << ": " << _bVAR[i] << std::endl;
    for( unsigned i=0; i<_ndep; i++ )
      os << " X" << i << ": " << PMx[i] << std::endl;
    _output_stats( _stats_ae, os );
  }
  return NORMAL;

}

template <typename T, typename PMT, typename PVT>
inline typename AEBND<T,PMT,PVT>::STATUS
AEBND<T,PMT,PVT>::_ge
( std::vector<FFVar>&X, std::ostream&os )
{
 try{
  _stats_ae.maxIter = 1;
  for( unsigned iblk=0; iblk<_bAE.size(); iblk++ ){
    // Construct LHS matrix and RHS vector
    const unsigned ndepblk = ( iblk==_bAE.size()-1? _ndep: _bAE[iblk+1] ) - _bAE[iblk];
    const unsigned posblk = _ndep - _bAE[iblk] - ndepblk;
    FFVar*Xblk = X.data() + posblk;
    const FFVar*F = _pAE.data() + posblk;
    std::vector<FFVar> A(ndepblk*ndepblk), b(ndepblk);
#ifdef MC__AEBND_DEBUG
    std::cout << "\nBefore forward substitution:\n";
#endif
    for( unsigned i=0; i<ndepblk*ndepblk; i++ ){
      A[i] = _pAEJAC[iblk][i];
#ifdef MC__AEBND_DEBUG
      std::cout << "A[" << i << "]: " << A[i] << std::endl;
#endif
    }
    std::vector<FFVar> zeros( ndepblk, 0. );
    const FFVar*F0 = _dag->compose( ndepblk, F, ndepblk, _pVAR.data()+posblk, zeros.data() );
    for( unsigned i=0; i<ndepblk; i++ ){
      b[i] = -F0[i];
#ifdef MC__AEBND_DEBUG
      std::cout << "b[" << i << "]: " << b[i] << std::endl;
#endif
    }
    delete[] F0;

    // Proceed with forward substitution
    for( unsigned k=0; k<ndepblk-1; k++ ){
      for( unsigned i=k+1; i<ndepblk; i++ ){
        FFVar factor = A[_ndx(i,k,ndepblk)] / A[_ndx(k,k,ndepblk)];
        for( unsigned j=k+1; j<ndepblk; j++ )
          A[_ndx(i,j,ndepblk)] -= factor * A[_ndx(k,j,ndepblk)];
        b[i] -= factor * b[k];
      }
    }
#ifdef MC__AEBND_DEBUG
    std::cout << "\nAfter forward substitution:\n";
    for( unsigned i=0; i<ndepblk*ndepblk; i++ )
      std::cout << "A[" << i << "]: " << A[i] << std::endl;
    for( unsigned i=0; i<ndepblk; i++ )
      std::cout << "b[" << i << "]: " << b[i] << std::endl;
#endif

    // Proceed with backward elimination
    Xblk[ndepblk-1] = b[ndepblk-1] / A[_ndx(ndepblk-1,ndepblk-1,ndepblk)];
    for( unsigned i=ndepblk; i>0; i-- ){
      FFVar sum = b[i-1];
      for( unsigned j=i; j<ndepblk; j++ )
        sum -= A[_ndx(i-1,j,ndepblk)] * Xblk[j];
      Xblk[i-1] = sum / A[_ndx(i-1,i-1,ndepblk)];
    }

    // Compose with Solutions of previous blocks
    if( iblk ){
#ifdef MC__AEBND_DEBUG
      std::cout << "\nComposition with previous solutions:\n";
      for( unsigned i=0; i<_ndep-posblk-ndepblk; i++ )
        std::cout << (_pVAR.data()+posblk+ndepblk)[i] << " -> " << (Xblk+ndepblk)[i] << std::endl;
#endif
      const FFVar*Xblkcomp = _dag->compose( ndepblk, Xblk, _ndep-posblk-ndepblk,
        _pVAR.data()+posblk+ndepblk, Xblk+ndepblk );
      for( unsigned i=0; i<ndepblk; i++ ) Xblk[i] = Xblkcomp[i];
      delete[] Xblkcomp;
    }
#ifdef MC__AEBND_DEBUG
    std::cout << "\nAfter backward elimination:\n";
    for( unsigned i=0; i<ndepblk*ndepblk; i++ )
      std::cout << "A[" << i << "] :" << A[i] << std::endl;
    for( unsigned i=0; i<ndepblk; i++ )
      std::cout << "b[" << i << "] :" << b[i] << std::endl;
    for( unsigned i=0; i<ndepblk; i++ )
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

//! @fn template <typename T, typename PMT, typename PVT> inline typename AEBND<T,PMT,PVT>::STATUS AEBND<T,PMT,PVT>::bounds(
//! FFVar*X, std::ostream&os=std::cout )
//!
//! This function computes a symbolic expression of the solution set of
//! the parametric AEs:
//!   - <a>X</a> [output] symbolic solution
//!   - <a>os</a> [input] output stream (default: std::cout)
//! .
//! The return value is the status.
template <typename T, typename PMT, typename PVT>
inline typename AEBND<T,PMT,PVT>::STATUS
AEBND<T,PMT,PVT>::solve
( FFVar*X, std::ostream&os )
{
  _init_stats( _stats_ae );
  if( !_issetup )   return FAILURE;
  if( _singstruct ) return SINGULAR;
  std::vector<FFVar> VAR(_ndep);
  STATUS flag = NORMAL;

  switch( options.BOUNDER ){

  // Gauss Elimination method
  case Options::AUTO:
  case Options::GE:
    if( !_islindep ) throw Exceptions( Exceptions::DAG );
    flag = _ge( VAR, os );
    break;

  // Gauss-Seidel method
  case Options::GS: case Options::GSS:
  case Options::KRAW: case Options::KRAWS:
    throw Exceptions( Exceptions::DAG );
  }

  _final_stats( _stats_ae );
  if( flag != NORMAL ) return flag;

  for( unsigned i=0; i<_ndep; i++ ) X[i] = VAR[_bVAR[i]];
  if( options.DISPLAY >= 1 ){
    std::cout << std::endl << "Symbolic Solution:\n";
    for( unsigned i=0; i<_ndep; i++ )
      os << " X" << i << ": " << X[i] << std::endl;
    _output_stats( _stats_ae, os );
  }
  return NORMAL;
}

} // end namescape mc

#endif

