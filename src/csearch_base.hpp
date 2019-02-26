// Copyright (C) 2015 Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

//TODO: 
//- [OK]    Implement SBB method
//- [OK]    Enable use of polynomial models in relaxations
//- [OK]    Enable multiple rounds of PWR refinement
//- [NO]    Exclude linear variables from bound contraction
//- [TO DO] Enable KKT cuts and reduction constraints
//- [TO DO] Make it possible to add/remove a constraint from the model?
//- [TO DO] Exploit the degree of separability to introduce auxiliary variables

#ifndef MC__CSEARCH_BASE_HPP
#define MC__CSEARCH_BASE_HPP

#include <stdexcept>

#include "sbb.hpp"
#include "sbp.hpp"
#include "lprelax_base.hpp"
#include "base_nlp.hpp"
#include "scmodel.hpp"

#undef MC__NLGO_DEBUG
//#undef MC__NLGO_TRACE
#undef MC__NLGO_SHOW_BREAKPTS

namespace mc
{

//! @brief C++ base class for complete search algorithms
////////////////////////////////////////////////////////////////////////
//! mc::CSEARCH_BASE is a C++ base class (pure virtual) for complete
//! search algorithms. Relaxations for nonlinear / nonconvex
//! participating functions are generated using MC++.
////////////////////////////////////////////////////////////////////////
template <typename T>
class CSEARCH_BASE:
  public virtual BASE_OPT,
  public virtual LPRELAX_BASE<T>,
  protected SBB<T>,
  protected SBP<T>
{
  // Typedef's
  // Monomial representation: <total order, <variable index, order>>
  typedef std::pair< unsigned, std::map< unsigned, unsigned > > t_expmon;
  typedef std::map< t_expmon, double > t_coefmon;
  typedef std::set< FFVar*, lt_FFVar > t_FFVar;
  typedef typename LPRELAX_BASE<T>::LP_STATUS LP_STATUS;

protected:
  //! @brief pointer to problem DAG
  FFGraph* _dag;
  //! @brief number of reduced-space variables in problem
  unsigned _nrvar;
  //! @brief number of reduced-space dependents in problem
  unsigned _nrdep;
  //! @brief number of variables in problem (both independent & dependent)
  unsigned _nvar;
  //! @brief decision variables (both independent & dependent)
  std::vector<FFVar> _var;
  //! @brief variable types in problem (both independent & dependent)
  std::vector<unsigned> _tvar;
  //! @brief reordering of decision variables into independents and dependents (forward)
  std::vector<unsigned> _var_fperm;
  //! @brief reordering of decision variables into independents and dependents (reverse)
  std::vector<unsigned> _var_rperm;
  //! @brief objective (type, cost variable, cost multiplier)
  std::tuple< std::vector<t_OBJ>, std::vector<FFVar>, std::vector<FFVar> > _obj;
  //! @brief number of constraints in problem
  unsigned _nctr;
  //! @brief constraints (types, constraint variables, constraint multipliers), including dependent equations
  std::tuple< std::vector<t_CTR>, std::vector<FFVar>, std::vector<FFVar>, std::vector<bool> > _ctr;
  //! @brief number of equality constraints
  unsigned _neq;
  //! @brief equality constraints
  std::vector<FFVar> _eq;

  //! @brief pointer to subset of dependent variables showing improvement using implicit contractors
  std::set<unsigned> _ndxdep;
  //! @brief pointer to subset of dependent variables showing improvement using implicit interval contractor
  std::set<unsigned> _Indxdep;
  //! @brief pointer to subset of dependent variables showing improvement using implicit polynomial model contractor
  std::set<unsigned> _CMndxdep;
  //! @brief boolean flag to ignore the dependents in the relaxations
  bool _ignore_deps;

  //! @brief list of operations for objective evaluation
  std::list<const FFOp*> _op_f;
  //! @brief vector of lists of operations for constraint evaluation
  std::vector< std::list<const FFOp*> > _op_g;
  //! @brief list of operations for Lagragian gradient evaluation
  std::list<const FFOp*> _op_dL;
  //! @brief Storage vector for function evaluation in double arithmetic
  //std::vector< double > _wk_D;

  //! @brief set of variables excluded from branching
  std::set<unsigned> _var_excl;
  //! @brief set of linear particpating variables in problem
  std::set<unsigned> _var_lin;
  //! @brief set of linear (objective or constraints) functions in problem
  t_FFVar _fct_lin;

  //! @brief Polyhedral image environment
  using LPRELAX_BASE<T>::_POLenv;
  //! @brief Storage vector for function evaluation in polyhedral relaxation arithmetic
  std::vector< PolVar<T> > _wk_POL;
  //! @brief Polyhedral image decision variables
  std::vector< PolVar<T> > _POLvar;
  //! @brief Polyhedral image objective auxiliary
  PolVar<T> _POLobjaux;
  //! @brief Polyhedral image scaled decision variables
  std::vector< PolVar<T> > _POLscalvar;
  //! @brief Polyhedral image of Chebyshev basis map
  std::map< t_expmon, PolVar<T> > _POLbasis;
  //! @brief Polyhedral image objective variable
  PolVar<T> _POLobj;
  //! @brief Polyhedral image constraint variables
  std::vector< PolVar<T> > _POLctr;

  //! @brief Chebyshev model environment
  SCModel<T>* _CMenv;
  //! @brief Chebyshev variables
  std::vector< SCVar<T> > _CMvar;
  //! @brief Chebyshev objective variable
  SCVar<T> _CMobj;
  //! @brief Chebyshev constraint variables
  std::vector< SCVar<T> > _CMctr;
  //! @brief Storage vector for function evaluation in Chebyshev arithmetic
  std::vector< SCVar<T> > _wk_CM;

  //! @brief Chebyshev reduced-space [-1,1] scaled model environment
  SCModel<T>* _CMrenv;
  //! @brief Chebyshev reduced-space [-1,1] scaled variables
  std::vector< SCVar<T> > _CMrbas;
  //! @brief Chebyshev reduced-space variables
  std::vector< SCVar<T> > _CMrvar;
  //! @brief Chebyshev reduced-space dependents
  std::vector< SCVar<T> > _CMrdep;
  //! @brief Interval reduced-space variables
  std::vector<T> _Irvar;
  //! @brief Interval reduced-space dependents
  std::vector<T> _Irdep;

  //! @brief Current incumbent value
  double _f_inc;
  //! @brief Variable values at current incumbent
  std::vector<double> _p_inc;

  //! @brief maximum number of values displayed in a row
  static const unsigned int _LDISP = 4;
  //! @brief reserved space for integer variable display
  static const unsigned int _IPREC = 9;
  //! @brief reserved space for double variable display
  static const unsigned int _DPREC = 6;
  //! @brief stringstream for displaying results
  std::ostringstream _odisp;

private:
  //! @brief Chebyshev basis map
  std::map< t_expmon, FFVar > _basis;
  //! @brief Internal DAG variable for [-1,1] scaled variables
  std::vector< FFVar > _scalvar;
  //! @brief Internal DAG variable for objective function
  FFVar _objvar;
  //! @brief Internal DAG variable for Lagrangian function
  FFVar _lagr;
  //! @brief Lagrangian gradient
  const FFVar* _lagr_grad;

public:
  //! @brief Constructor
  CSEARCH_BASE
    ()
    : _dag(0), _nvar(0), _nctr(0), _CMenv(0), _CMrenv(0), _lagr_grad(0)
    {}

  //! @brief Destructor
  virtual ~CSEARCH_BASE
    ()
    {
      _cleanup();
      delete _CMenv;
      delete _CMrenv;
      // No need to delete _NLPSLV
    }

protected:
  //! @brief Setup internal variables, including structure detection
  template <typename OPT>
  void _setup
    ( const OPT&options, FFGraph*dag, std::ostream&os );

  //! @brief Initialize incumbent
  template <typename SLVLOC, typename OPT, typename STAT>
  void _solve_init_inc
    ( SLVLOC&locopt, const OPT&options, STAT&stats, const T*P,
      const unsigned*tvar, const double*p0, std::ostream&os );
  //! @brief Preprocess variable domain
  template <typename OPT, typename STAT>
  void _solve_init_rel
    ( const OPT&options, STAT&stats, const T*P0,
      const unsigned*tvar, std::vector<T>&P, std::ostream&os );
  //! @brief Solve optimization model to global optimality within variable range <a>P</a> and for variable types <a>tvar</a>
  template <typename SLVLOC, typename OPT, typename STAT>
  int _solve_pwl
    ( SLVLOC&locopt, const OPT&options, STAT&stats, const T*P,
      const unsigned*tvar, const double*p0, std::ostream&os );
  //! @brief Solve optimization model to global optimality within variable range <a>P</a> and for variable types <a>tvar</a> -- return value is mc::SBB::STATUS
  template <typename SLVLOC, typename OPT, typename STAT>
  int _solve_sbb
    ( SLVLOC&locopt, const OPT&options, STAT&stats, const T*P,
      const unsigned*tvar, const double*p0, std::ostream&os );
  //! @brief Solve constraint system using complete search within variable range <a>P</a> and for variable types <a>tvar</a> -- return value is mc::SBB::STATUS
  template <typename OPT, typename STAT>
  double _solve_sbp
    ( const OPT&options, STAT&stats, const T*P, const unsigned*tvar,
      std::ostream&os );

  //! @brief Solve (continuous) optimization model to local optimality using IPOPT in variable range <a>P</a> and from initial point <a>p0</a> -- return value is IPOPT status
  template <typename SLVLOC, typename STAT>
  int _local
    ( SLVLOC&locopt, STAT&stats, const T*P, const double*p0=0, const bool reset=true );
  //! @brief Solve relaxed optimization model (or relaxed feasibility model if <a>feastest</a> is true) within variable range <a>P</a> and for variable types <a>tvar</a>
  template <typename OPT, typename STAT>
  LP_STATUS _relax
    ( const OPT&options, STAT&stats, const T*P, const unsigned*tvar=0,
      const unsigned refine=0, const bool reset=true, const bool feastest=false );
  //! @brief Solve bound contraction problems from relaxed model, starting with variable range <a>P</a>, for variable types <a>tvar</a> and incumbent value <a>inc</a> (or constraint backoff if <a>feastest</a> is true), and using the options specified in <a>NLGO::OPT::DOMREDMAX</a> and <a>NLGO::OPT::DOMREDTHRES</a> -- returns updated variable bounds <a>P</a>, and number of iterative refinements <a>nred</a> 
  template <typename OPT, typename STAT>
  LP_STATUS _contract
    ( const OPT&options, STAT&stats, T*P, unsigned&nred, const unsigned*tvar=0,
      const double*inc=0, const bool reset=true, const bool feastest=false );
 //! @brief Tighten parameter bounds <a>P</a> for NLP model relaxation and incumbent value <a>inc</a> (or constraint backoff is <a>feastest</a> is true)
  template <typename OPT, typename STAT>
  LP_STATUS _contract
    ( const OPT&options, STAT&stats, T*P, const unsigned*tvar=0, const double*inc=0,
      const bool feastest=false );
  //! @brief Solve bound contraction problem for lower/upper component <a>ip</a> from relaxed model, starting with variable range <a>P</a> and incumbent value <a>inc</a> (or constraint backoff if <a>feastest</a> is true)
  template <typename OPT, typename STAT>
  LP_STATUS _contract
    ( const OPT&options, STAT&stats, const T*P, const unsigned ip, const bool uplo,
      const unsigned*tvar=0, const double*inc=0, const bool feastest=false );

  //! @brief Set model polyhedral relaxation
  template <typename OPT, typename STAT>
  void _set_polrelax
    ( const OPT&options, STAT&stats, const T*P, const unsigned*tvar,
      const bool feastest=false );
  //! @brief Update model polyhedral relaxation (bounds and types)
  template <typename OPT, typename STAT>
  void _update_polrelax
    ( const OPT&options, STAT&stats, const T*P, const unsigned*tvar,
      const bool feastest=false );
  //! @brief Refine model polyhedral relaxation by adding breakpoints
  template <typename OPT, typename STAT>
  void _refine_polrelax
    ( const OPT&options, STAT&stats, const T*P, const unsigned*tvar,
      const bool feastest=false );

  //! @brief Set dependent Chebyshev variables
  template <typename OPT, typename STAT>
  bool _set_depcheb
    ( const OPT&options, STAT&stats, const T*P=0 );
  //! @brief Update dependent Chebyshev variables
  template <typename OPT, typename STAT>
  bool _update_depcheb
    ( const OPT&options, STAT&stats, const T*P=0 );
  //! @brief Compute dependent Chebyshev variables
  virtual bool _get_depcheb
    ( const bool reset=true ) = 0;
  //! @brief Rescale dependent Chebyshev variables on reduced domain
  template <typename OPT, typename STAT>
  void _rescale_depcheb
    ( const OPT&options, STAT&stats, const T*P );
  //! @brief Set dependent interval variables
  template <typename STAT>
  bool _set_depbnd
    ( STAT&stats, const T*P );
  //! @brief Update dependent interval variables
  template <typename STAT>
  bool _update_depbnd
    ( STAT&stats, const T*P );
  //! @brief Compute dependent interval variables
  virtual bool _get_depbnd
    ( const bool reset=true ) = 0;
  //! @brief Tighten dependent variable bounds
  template <typename STAT>
  void _tighten_depbnd
    ( STAT&stats, T*P );

  //! @brief Set model linear cuts
  template <typename OPT>
  void _set_lincuts
    ( const OPT&options, const bool feastest );
  //! @brief Set model polyhedral outer-approximation cuts
  template <typename OPT>
  void _set_poacuts
    ( const OPT&options, const bool feastest );
  //! @brief Set model polyhedral Chebyshev-derived cuts
  template <typename OPT>
  void _set_chebcuts
    ( const OPT&options, const bool feastest );

  //! @brief Set options of SBB solver
  template <typename OPT>
  void _set_SBBoptions
    ( const OPT&options );
  //! @brief Set options of SBP solver
  template <typename OPT>
  void _set_SBPoptions
    ( const OPT&options ) const;
  //! @brief Set local optimizer
  virtual void _set_SLVLOC
    () = 0;
  //! @brief Get local optimum
  virtual const double* _get_SLVLOC
    ( const SOLUTION_NLP&locopt ) = 0;

  //! @brief Set branching scores for current node
  virtual void _set_branching_scores
    ( std::map<unsigned,double>&scores ) = 0;
//  //! @brief Append branching scores based on Chebyshev model enclosures
//  void _add_chebscores
//    ( std::map<unsigned,double>&scores, const FFVar&FV,
//      const SCVar<T>&SCV, const double weight, const bool linexcl=false ) const;
//  //! @brief Set up branching scores based on Chebyshev model enclosures
//  template <typename OPT>
//  void _set_chebscores
//    ( const OPT&options, std::map<unsigned,double>&scores, const bool feastest=false );
//  //! @brief Set up branching scores based on block structure
//  template <typename OPT>
//  void _set_blkscores
//    ( const OPT&options, std::map<unsigned,double>&scores );
//  //! @brief Set up branching scores based on block structure
//  virtual void _get_blkscores
//    ( std::map<unsigned,double>&scores ) = 0;

  //! @brief Solve LP model
  using LPRELAX_BASE<T>::_solve_LPmodel;
  //! @brief Set relaxed objective in LP model
  using LPRELAX_BASE<T>::_set_LPrelax;
  //! @brief Set parameter bound contracting objective in LP model
  using LPRELAX_BASE<T>::_set_LPcontract;
  //! @brief Set variables and cuts in LP model
  using LPRELAX_BASE<T>::_set_LPcuts;
  //! @brief Value of PolImg variable <a>X</a> after last LP optimization
  using LPRELAX_BASE<T>::_get_variable;

  //! @brief Value of DAG variable <a>X</a> after last LP optimization
  using LPRELAX_BASE<T>::get_variable;
  //! @brief Optimal cost value after last LP optimization
  using LPRELAX_BASE<T>::get_objective;
  //! @brief Status after last LP optimization
  using LPRELAX_BASE<T>::get_status;
  //! @brief Pointer to last LP optimization model and solution
  using LPRELAX_BASE<T>::get_relaxed_model;

  //! @brief Function computing Hausdorff distance between intervals
  using LPRELAX_BASE<T>::_dH;
  //! @brief Function computing Hausdorff distance between interval vectors
  using LPRELAX_BASE<T>::_reducrel;

  //! @brief Cleanup gradient/hessian storage
  void _cleanup()
  { delete[] _lagr_grad; _lagr_grad = 0; }

  //! @brief User-function to subproblems in SBB
  template <typename SLVLOC, typename OPT, typename STAT>
  typename SBB<T>::STATUS _subproblems
    ( SLVLOC&locopt, OPT&options, STAT&stats, const typename SBB<T>::TASK task,
      SBBNode<T>*node, std::vector<double>&p, double&f, const double INC );
  //! @brief Subproblem for local optimization
  template <typename SLVLOC, typename OPT, typename STAT>
  typename SBB<T>::STATUS _subproblem_local
    ( SLVLOC&locopt, const OPT&options, STAT&stats, const std::vector<T>&P,
      std::vector<double>&p, double&f );
  //! @brief Subproblem for feasibility test
  template <typename SLVLOC, typename OPT>
  typename SBB<T>::STATUS _subproblem_feasibility
    ( SLVLOC&locopt, const OPT&options, const double*p, double&f );
  //! @brief Subproblem for relaxation
  template <typename OPT, typename STAT>
  typename SBB<T>::STATUS _subproblem_relaxed
    ( const OPT&options, STAT&stats, std::vector<T>&P, std::vector<double>&p,
      double&f, const double INC );
  //! @brief Update incumbent
  void _update_incumbent
    ( const double*p, const double&f );
//  //! @brief Update incumbent
//  bool _check_feasibility
//    ( const double*p, const double FEASTOL );

  //! @brief User-function to subproblem assessment in SBP
  template <typename OPT, typename STAT>
  typename SBP<T>::STATUS _assess
    ( OPT&options, STAT&stats, SBPNode<T>*node );
  //! @brief Subproblem for node assessment
  template <typename OPT, typename STAT>
  typename SBP<T>::STATUS _subproblem_assess
    ( const OPT&options, STAT&stats, std::vector<T>&P );
  //! @brief Perform inclusion test based on constraint bounds
  template <typename OPT, typename BND>
  typename SBP<T>::STATUS _check_inclusion
    ( const OPT&options, const std::vector<BND>&C ) const;

  //! @brief Initialize display
  virtual void _display_pwl_init
    () = 0;
  //! @brief Final display
  virtual void _display_pwl_final
    ( const unsigned iter ) = 0;
  //! @brief Add double to display
  template <typename OPT>
  void _display_add
    ( const OPT&options, const double dval );
  //! @brief Add unsigned int to display
  template <typename OPT>
  void _display_add
    ( const OPT&options, const int ival );
  //! @brief Add string to display
  template <typename OPT>
  void _display_add
    ( const OPT&options, const std::string &sval );
  //! @brief Display current buffer stream and reset it
  void _display_flush
    ( std::ostream&os );

  //! @brief Private methods to block default compiler methods
  CSEARCH_BASE
    (const CSEARCH_BASE&);
  CSEARCH_BASE& operator=
    (const CSEARCH_BASE&);
};

template <typename T>
template <typename OPT>
inline void
CSEARCH_BASE<T>::_setup
( const OPT&options, FFGraph*dag, std::ostream&os )
{
  // setup DAG pointer
  _dag = dag;

  // setup [-1,1] scaled variables for Chebyshev model environment
  // do not downsize to avoid adding more variables into DAG
  const unsigned nscalvar = _scalvar.size();
  if( nscalvar < _nvar ){
    _scalvar.resize( _nvar );
    for( unsigned i=nscalvar; i<_nvar; i++ )
      _scalvar[i].set( _dag );
  }

  // setup for objective and constraints evaluation
  _objvar.set( _dag );
  if( std::get<0>(_obj).size() )
    _op_f = _dag->subgraph( 1, std::get<1>(_obj).data() );
  _op_g.resize( _nctr );
  for( unsigned i=0; i<_nctr; i++ ){
    _op_g[i] = _dag->subgraph( 1, std::get<1>(_ctr).data()+i );
  }

  // identify linear objective/constraint functions
  _fct_lin.clear();
  if( std::get<0>(_obj).size() ){
    FFDep fdep = std::get<1>(_obj)[0].dep();
    auto it = fdep.dep().cbegin();
    bool islin = true;
    for( ; islin && it != fdep.dep().cend(); ++it )
      if( it->second != FFDep::L ) islin = false;
    if( islin ) _fct_lin.insert( &std::get<1>(_obj)[0] );
  }
  for( unsigned j=0; j<_nctr; j++ ){
    FFDep gdep = std::get<1>(_ctr)[j].dep();
    auto it = gdep.dep().cbegin();
    bool islin = true;
    for( ; islin && it != gdep.dep().cend(); ++it )
      if( it->second != FFDep::L ) islin = false;
    if( islin ) _fct_lin.insert( &std::get<1>(_ctr)[j] );
  }

  if( options.DISPLAY ){
    os << "LINEAR OBJECTIVE/CONSTRAINT FUNCTIONS:    " << _fct_lin.size() << std::endl
       << "NONLINEAR OBJECTIVE/CONSTRAINT FUNCTIONS: " << _nctr-_fct_lin.size()+std::get<0>(_obj).size() << std::endl;
  }

  // Initialize local solution
  _set_SLVLOC();
}

template <typename T>
template <typename OPT, typename STAT>
inline void
CSEARCH_BASE<T>::_solve_init_rel
( const OPT&options, STAT&stats, const T*P0,
  const unsigned*tvar, std::vector<T>&P, std::ostream&os )
{
  tvar? _tvar.assign( tvar, tvar+_nvar ): _tvar.clear();
  P.assign( P0, P0+_nvar );
#ifdef MC__CSEARCH_SHOW_BOXES
  std::cout << "\nInitial Box:\n";
  for( unsigned i=0; i<_nvar; i++ )
    //if( _exclude_vars.find(i) == _exclude_vars.end() )
      std::cout << "P[" << i << "] = " << P[i] << std::endl;
#endif

  // Preprocessing using domain contraction
  if( options.PREPROC ){
    _ignore_deps = true;
    unsigned nred = 0;
    _contract( options, stats, P.data(), nred, tvar, !_p_inc.empty()? &_f_inc: 0, true );
#ifdef MC__CSEARCH_SHOW_BOXES
  std::cout << "\nReduced Box:\n";
  for( unsigned i=0; i<_nvar; i++ )
    //if( _exclud_vars.find(i) == _exclude_vars.end() )
      std::cout << "P[" << i << "] = " << P[i] << std::endl;
#endif
  }

  // Set-up NLP relaxed solver
  _ignore_deps = false;
  _set_depbnd( stats, P.data() );
  _set_depcheb( options, stats );
  _tighten_depbnd( stats, P.data() );
  _set_polrelax( options, stats, P.data(), tvar );
  //std::cout << *_dag;
  //std::cout << _POLenv;
}

template <typename T>
template <typename SLVLOC, typename OPT, typename STAT>
inline void
CSEARCH_BASE<T>::_solve_init_inc
( SLVLOC&locopt, const OPT&options, STAT&stats, const T*P,
  const unsigned*tvar, const double*p0, std::ostream&os )
{
  bool isMIP = false;
  for( unsigned ivar=0; tvar && !isMIP && ivar<_nvar; ivar++ )
    if( tvar[ivar] ) isMIP = true;
  _p_inc.clear();

  // Check feasibility of user-supplied point
  double f_loc;
  if( p0 ){
    if( _subproblem_feasibility( locopt, options, p0, f_loc ) == SBB<T>::NORMAL )
      _update_incumbent( p0, f_loc );
  }

  // Find local solution
  _local( locopt, stats, P, p0 );
  const double* p_loc = _get_SLVLOC( locopt->solution() );

  // Rounding solution point to nearest integer in case of MIP
  if( isMIP
   && _subproblem_feasibility( locopt, options, p_loc, f_loc ) == SBB<T>::NORMAL ){
    std::vector<T> Pround(_nvar);
    for( unsigned ivar=0; ivar<_nvar; ivar++ )
      Pround[ivar] = tvar[ivar]? std::round( p_loc[ivar] ): P[ivar];
    _local( locopt, stats, Pround.data(), p_loc );
    p_loc = _get_SLVLOC( locopt->solution() );
  }

  // Update incumbent
  if( _subproblem_feasibility( locopt, options, p_loc, f_loc ) == SBB<T>::NORMAL )
    _update_incumbent( p_loc, f_loc );
}

template <typename T>
template <typename OPT, typename STAT>
inline double
CSEARCH_BASE<T>::_solve_sbp
( const OPT&options, STAT&stats, const T*P0, const unsigned*tvar,
  std::ostream&os )
{
  // Start B&P from scratch
  if( P0 ){
    stats.reset();
    stats.tALL -= cpuclock();

    // Initialize solve and preprocess
    std::vector<T> P;
    _p_inc.clear();
    _solve_init_rel( options, stats, P0, tvar, P, os );

    // Set variable types and exclusions from B&P
    //_var_excl.insert( _var_lin.begin(), _var_lin.end() ); // excluded branching variables
    for( unsigned i=0; i<_nvar; i++ )
      if( Op<T>::diam( P[i] ) <= options.DOMREDMIG ) _var_excl.insert( i );
    SBP<T>::variables( _nvar, P.data(), _var_excl );
  }

  // Resume B&P after interruption
  else{
    stats.tALL -= cpuclock();
  }

  // Call B&P solver
  _set_SBPoptions( options );
  double open = SBP<T>::solve( os, P0? false: true );

  stats.tALL += cpuclock();
  return open;
}

template <typename T>
template <typename SLVLOC, typename OPT, typename STAT>
inline int
CSEARCH_BASE<T>::_solve_sbb
( SLVLOC&locopt, const OPT&options, STAT&stats, const T*P0,
  const unsigned*tvar, const double*p0, std::ostream&os )
{
  stats.reset();
  stats.tALL -= cpuclock();

  // Initialize solve and preprocess
  std::vector<T> P;
  _solve_init_inc( locopt, options, stats, P0, tvar, p0, os );
  _solve_init_rel( options, stats, P0, tvar, P, os );

  // Set variable types and exclusions from B&B
  _set_SBBoptions( options );
  _var_excl.insert( _var_lin.begin(), _var_lin.end() ); // excluded branching variables
  for( unsigned i=0; i<_nvar; i++ )
    if( Op<T>::diam( P[i] ) <= options.DOMREDMIG ) _var_excl.insert( i );
  SBB<T>::variables( _nvar, P.data(), _var_excl );

  // Run SBB solver
  const double* fINC = (!_p_inc.empty()? &_f_inc: 0);
  const double* pINC = (!_p_inc.empty()? _p_inc.data(): p0);
  switch( std::get<0>(_obj)[0] ){
    case MIN: SBB<T>::solve( SBB<T>::MIN, pINC, fINC, os ); break;
    case MAX: SBB<T>::solve( SBB<T>::MAX, pINC, fINC, os ); break;
  }
  _p_inc = SBB<T>::_p_inc;
  _f_inc = SBB<T>::_f_inc;

  stats.tALL += cpuclock();
  return SBB<T>::_status;
}

template <typename T>
template <typename SLVLOC, typename OPT, typename STAT>
inline int
CSEARCH_BASE<T>::_solve_pwl
( SLVLOC&locopt, const OPT&options, STAT&stats, const T*P0, 
  const unsigned*tvar, const double*p0, std::ostream&os )
{
  stats.reset();
  stats.tALL -= cpuclock();

  // Initialize solve and preprocess
  std::vector<T> P;
  _solve_init_inc( locopt, options, stats, P0, tvar, p0, os );
  _solve_init_rel( options, stats, P0, tvar, P, os );
  _display_pwl_init();
  std::vector<T> P_loc(_nvar);
  std::vector<double> p_rel(_nvar);
  const double* p_loc;
  double f_loc;

  // Iterative relaxation solution and refinement
  unsigned iter = 1, nred = 0;
  for( ; ; ++iter ){
    // iteration and preprocessing display
    std::ostringstream odisp, onred;
    odisp << std::right << std::setw(_IPREC) << iter;
    if( nred ) onred << "  R" << nred;
    else       onred << "-";
    odisp << std::right << std::setw(_IPREC) << onred.str();
    _display_add( options, odisp.str() );

    // Set-up relaxation and solve relaxed model
    _relax( options, stats, P.data(), tvar, 0, false );
    switch( get_status() ){
     case LPRELAX_BASE<T>::LP_OPTIMAL:
      _display_add( options, get_objective() );
      break;
     case LPRELAX_BASE<T>::LP_INFEASIBLE:
      _display_add( options, "-" );
      _display_add( options, "-" );
      _display_add( options, cpuclock()-stats.tALL );
      _display_add( options, "INFEASIBLE" );
      _display_flush( os );
      _display_pwl_final( iter );
      _display_flush( os );
      return get_status();
     default:
      _display_add( options, "-" );
      _display_add( options, "-" );
      _display_add( options, cpuclock()-stats.tALL );
      _display_add( options, "FAILURE" );
      _display_flush( os );
      _display_pwl_final( iter );
      _display_flush( os );
      return get_status();
    }

    // Solve local NLP model (integer variable bounds fixed to relaxed solution if MIP)
    unsigned ivar = 0;
    for( auto itv=_var.begin(); itv!=_var.end(); ++itv, ivar++ ){
      p_rel[ivar] = get_variable( *itv );
      P_loc[ivar] = (tvar && tvar[ivar])? p_rel[ivar]: P[ivar];
    }
    locopt->solve( P_loc.data(), p_rel.data() );
    p_loc = _get_SLVLOC( locopt->solution() );

    // Update incumbent
    if( _subproblem_feasibility( locopt, options, p_loc, f_loc ) == SBB<T>::NORMAL )
      _update_incumbent( p_loc, f_loc );

    // Test termination criteria
    if( !_p_inc.empty() ){
      _display_add( options, _f_inc );
      _display_add( options, cpuclock()-stats.tALL );
      if( std::fabs( _f_inc - get_objective() ) <= options.CVATOL 
       || std::fabs( _f_inc - get_objective() ) <= 0.5 * options.CVRTOL
          * std::fabs( _f_inc + get_objective() ) ){
        _display_add( options, "OPTIMAL" );
        _display_flush( os );
        break;
      }
    }
    else{
      _display_add( options, "-" );
      _display_add( options, cpuclock()-stats.tALL );
    }
    _display_add( options, "REFINE" );
    _display_flush( os );

    if( options.MAXITER && iter >= options.MAXITER ){
      _display_add( options, "INTERRUPT" );
      _display_flush( os );
      break;
    }

    // Refine relaxation via additional breakpoints
    _refine_polrelax( options, stats, P.data(), tvar );

    // Apply domain contraction
    // Do NOT test for infeasibility here, because contraction problem may
    // become infeasible due to round-off in LP solver
    _contract( options, stats, P.data(), nred, tvar, !_p_inc.empty()? &_f_inc: 0, false );

    //return 0;
  }

  // Final display
  _display_pwl_final( iter );
  _display_flush( os );

  return get_status();
}

template <typename T>
template <typename SLVLOC, typename STAT>
inline int
CSEARCH_BASE<T>::_local
( SLVLOC&locopt, STAT&stats, const T*P, const double*p0,
  const bool init )
{
  stats.tLOCSOL -= cpuclock();
  stats.nLOCSOL++;

  // Initialize local solver
  if( init ) _set_SLVLOC();

  // Solve local NLP model
  std::vector<double> p(_nvar);
  for( unsigned ip=0; ip<_nvar; ip++ )
    p[ip] = ( p0? p0[ip]: Op<T>::mid(P[ip]) );
  int status = locopt->solve( P, p.data() );

  stats.tLOCSOL += cpuclock();
  //std::cout << "tlocopt: " << stats.tlocopt << std::endl
  //          << "nlocopt: " << stats.nlocopt << std::endl;
  return status;
}

template <typename T>
template <typename OPT, typename STAT>
inline typename LPRELAX_BASE<T>::LP_STATUS
CSEARCH_BASE<T>::_relax
( const OPT&options, STAT&stats, const T*P, const unsigned*tvar,
  const unsigned nref, const bool reset, const bool feastest )
{
  // Reset/update polyhedral image, LP variables and cuts on request
  if( reset ){
    _set_depbnd( stats, P );
    _set_depcheb( options, stats );
    //_tighten_depbnd( stats, P );
    _set_polrelax( options, stats, P, tvar, feastest );
  }
  //else{
  //  _update_depbnd( stats, P );
  //  _update_depcheb( options, stats, P );
  //  _update_polrelax( options, stats, P, tvar, feastest );
  //}
#ifdef MC__NLGO_DEBUG
  std::cout << _POLenv;
#endif

  for( unsigned iref=0; ; iref++ ){
    // Set-up relaxed objective, options, and solve relaxed model
    _set_LPrelax( &_POLobjaux, std::get<0>(_obj)[0], feastest );
    _solve_LPmodel( options, stats, _var );

    // Break if relaxation unsuccessful or refinement iteration exceeded
    if( get_status() != LPRELAX_BASE<T>::LP_OPTIMAL || iref >= nref ) break;

    // Refine relaxation via additional breakpoints
    _refine_polrelax( options, stats, P, tvar );
  }

  return get_status();
}

template <typename T>
template <typename OPT, typename STAT>
inline typename LPRELAX_BASE<T>::LP_STATUS
CSEARCH_BASE<T>::_contract
( const OPT&options, STAT&stats, const T*P, const unsigned ip,
  const bool uplo, const unsigned*tvar, const double*inc,
  const bool feastest )
{
#ifdef MC__CSEARCH_DEBUG
  std::cout << "\nTIGHTENING OF VARIABLE " << ip << (uplo?"U":"L") << ":\n";
  std::cout << _POLenv;
#endif
  // Set-up relaxed objective, options, and solve relaxed model
  _set_LPcontract( &_POLvar[ip], uplo, &_POLobjaux, inc, std::get<0>(_obj)[0], feastest );
  _solve_LPmodel( options, stats, _var );
  //{ int dum; std::cin >> dum; }
  return get_status();
}

template <typename T>
template <typename OPT, typename STAT>
inline typename LPRELAX_BASE<T>::LP_STATUS
CSEARCH_BASE<T>::_contract
( const OPT&options, STAT&stats, T*P, const unsigned*tvar,
  const double*inc, const bool feastest )
{
  // solve reduction subproblems from closest to farthest from bounds
  std::multimap<double,std::pair<unsigned,bool>> vardomred, vardomredupd;
  std::pair<unsigned,bool> varini;
  for( unsigned ip=0; ip<_nvar; ip++ ){
    // do not contract variables whose domain is less than DOMREDMIG
    if( Op<T>::diam( P[ip] ) < options.DOMREDMIG ) continue;

    varini.first = ip;
    varini.second = false; // lower bound
    double dist = 1.;
    vardomred.insert( std::pair<double,std::pair<unsigned,bool>>(dist,varini) );
    varini.second = true;  // upper bound
    vardomred.insert( std::pair<double,std::pair<unsigned,bool>>(dist,varini) );
  }

  unsigned nred=0;
  for( ; !vardomred.empty(); nred++ ){
    // upper/lower range reduction for current subproblem
    auto itv = vardomred.begin();
    const unsigned ip = (*itv).second.first;
    const bool uplo = (*itv).second.second;
    double pL = Op<T>::l( P[ip] ), pU = Op<T>::u( P[ip] );
    _contract( options, stats, P, ip, uplo, tvar, inc, feastest );
    if( get_status() != LPRELAX_BASE<T>::LP_OPTIMAL ){
#ifdef MC__CSEARCH_PAUSE_INFEASIBLE
      int dum; std::cout << "Infeasible problem - PAUSED"; std::cin >> dum;
#endif
      return get_status();
    }
    switch( (int)uplo ){
     case false: // lower bound
      pL = get_objective();
      if( options.DOMREDBKOFF > 0. ) pL -= options.DOMREDBKOFF;
      if( !Op<T>::inter(  P[ip], P[ip], T(pL,pU+1.) ) ) P[ip] = pU;
      break;
     case true: // upper bound
      pU = get_objective(); 
      if( options.DOMREDBKOFF > 0. ) pU += options.DOMREDBKOFF;
      if( !Op<T>::inter(  P[ip], P[ip], T(pL-1.,pU) ) ) P[ip] = pL;
      break;
    }
#ifdef MC__CSEARCH_DEBUG
    std::cout << "  UPDATED RANGE OF VARIABLE #" << ip << ": " << P[ip] << std::endl;
#endif
    // update map of candidate reduction subproblems
    vardomredupd.clear();
    for( ++itv; itv!=vardomred.end(); ++itv ){
      const unsigned ip = (*itv).second.first;
      const bool uplo = (*itv).second.second;
      double dist = ( uplo? std::fabs( get_variable( _var[ip] ) - Op<T>::u( P[ip] ) ):
                            std::fabs( get_variable( _var[ip] ) - Op<T>::l( P[ip] ) ) )
                    / Op<T>::diam( P[ip] ); // consider relative distance to bound
      if( dist <= options.DOMREDTHRES ) continue;
      if( dist > (*itv).first ) dist = (*itv).first;
      vardomredupd.insert( std::pair<double,std::pair<unsigned,bool>>(dist,std::make_pair(ip,uplo)) );
    }
    vardomred.swap( vardomredupd );
  }
#ifdef MC__CSEARCH_DEBUG
    std::cout << "SOLVED " << nred << " RANGE REDUCTION LPs OUT OF " << 2*_nvar << std::endl;
#endif
  return get_status();
}

template <typename T>
template <typename OPT, typename STAT>
inline typename LPRELAX_BASE<T>::LP_STATUS
CSEARCH_BASE<T>::_contract
( const OPT&options, STAT&stats, T*P, unsigned&nred,
  const unsigned*tvar, const double*inc, const bool reset,
  const bool feastest )
{
  // Dependent bounds
  if( reset ){
    _set_depbnd( stats, P );
    _set_depcheb( options, stats );
    _tighten_depbnd( stats, P );
#ifdef MC__CSEARCH_SHOW_REDUC
    //std::cout << "Reduction #0: " << std::fixed << std::setprecision(1)
    //          << _reducrel( _nvar, P, P )*1e2 << "%\n";
#endif
    _set_polrelax( options, stats, P, tvar, feastest );
  }
  //else{
  //  _update_depbnd( stats, P );
  //  _update_depcheb( options, stats, P );
  //  _update_polrelax( options, stats, P, tvar, feastest );
  //}

  // Main loop for relaxation and domain reduction
  std::vector<T> P0( P, P+_nvar);
  std::vector<T> P1;
  for( nred=0; nred < options.DOMREDMAX; nred++ ){
    P1.assign( P, P+_nvar );
    if( nred ){
      _update_depcheb( options, stats, P );
      //_rescale_depcheb( options, stats, P );
      _update_polrelax( options, stats, P, tvar, feastest );
    }
    _contract( options, stats, P, tvar, inc, feastest );
    if( get_status() != LPRELAX_BASE<T>::LP_OPTIMAL ) return get_status();
    double vred = _reducrel( _nvar, P, P1.data(), P0.data() );
#ifdef MC__CSEARCH_SHOW_REDUC
    std::cout << "Reduction #" << nred+1 << ": "
              << std::fixed << std::setprecision(1) << vred*1e2 << "%\n";
#endif
    if( vred < options.DOMREDTHRES ) break;
  }
  return get_status();
}

template <typename T>
template <typename SLVLOC, typename OPT, typename STAT>
inline typename SBB<T>::STATUS
CSEARCH_BASE<T>::_subproblems
( SLVLOC&locopt, OPT&options, STAT&stats,
  const typename SBB<T>::TASK task, SBBNode<T>*node,
  std::vector<double>&p, double&f, const double INC )
{
  // Strong branching case
  unsigned DOMREDMAX = options.DOMREDMAX;
  if( node->strongbranch() && options.STGBCHDRMAX
   && options.STGBCHDRMAX < options.DOMREDMAX ){
    options.DOMREDMAX = options.STGBCHDRMAX; 
#if defined (MC__CSEARCH_DEBUG)
    std::cout << "on-going strong branching\n";
#endif
  }

  // Dependents in parent node
  _ndxdep = node->depend();

  typename SBB<T>::STATUS status = SBB<T>::FATAL;
  switch( task ){

  // UPPER BOUND
  case SBB<T>::UPPERBD:
   try{ 
    switch( std::get<0>(_obj)[0] ){
      case MIN: status = _subproblem_local( locopt, options, stats, node->P(), p, f ); break;
      case MAX: status = _subproblem_relaxed( options, stats, node->P(), p, f, INC );
        if( !node->strongbranch() ){
          _set_branching_scores( node->scores() );
          //_set_chebscores( options, node->scores(), false );
          node->depend() = _ndxdep;
        } break;
    }
   }
   catch(...){
     status = SBB<T>::FAILURE;
   } break;

  // LOWER BOUND
  case SBB<T>::LOWERBD:
   try{ 
    switch( std::get<0>(_obj)[0] ){
      case MAX: status = _subproblem_local( locopt, options, stats, node->P(), p, f ); break;
      case MIN: status = _subproblem_relaxed( options, stats, node->P(), p, f, INC );
        if( !node->strongbranch() ){
          _set_branching_scores( node->scores() );
          //_set_chebscores( options, node->scores(), false );
          node->depend() = _ndxdep;
        } break;
    }
   }
   catch(...){
     status = SBB<T>::FAILURE;
   } break;

  // FEASIBILITY TEST
  case SBB<T>::FEASTEST:
    status = _subproblem_feasibility( locopt, options, p.data(), f ); break;

  // PRE/POST-PROCESSING
  case SBB<T>::PREPROC:
  case SBB<T>::POSTPROC:
    status = SBB<T>::NORMAL; break;

  default:
    status = SBB<T>::FAILURE; break;
  }

  // Strong branching case
  if( node->strongbranch() ){
    options.DOMREDMAX = DOMREDMAX;
  }

  return status;
}

template <typename T>
template <typename OPT, typename STAT>
inline typename SBB<T>::STATUS
CSEARCH_BASE<T>::_subproblem_relaxed
( const OPT&options, STAT&stats, std::vector<T>&P,
  std::vector<double>&p, double&f, const double INC )
{
#ifdef MC__CSEARCH_SHOW_BOXES
  std::cout << "\nInitial Box:\n";
  for( unsigned i=0; i<_nvar; i++ )
    //if( _exclude_vars.find(i) == _exclude_vars.end() )
      std::cout << "P[" << i << "] = " << P[i] << std::endl;
#endif
  unsigned *tvar = !_tvar.empty()? _tvar.data(): 0;
  std::vector<T> P0 = P;

  // Update bounds for dependents
  if( _update_depbnd( stats, P.data() ) ){
    _ndxdep.insert( _Indxdep.begin(), _Indxdep.end() );
    if( _update_depcheb( options, stats ) )
      _ndxdep.insert( _CMndxdep.begin(), _CMndxdep.end() );
    _tighten_depbnd( stats, P.data() );
  }
  else{
    _CMndxdep.clear();
  }
#ifdef MC__CSEARCH_SHOW_REDUC
  std::cout << "Reduction #0: " << std::fixed << std::setprecision(1)
            << _reducrel( _nvar, P.data(), P0.data() )*1e2 << "%\n";
#endif
/*
  // Update bounds for dependents
  if( !_update_depbnd( stats, P.data() )
   || !_update_depcheb( options, stats ) )
    return SBB<T>::INFEASIBLE;
  _tighten_depbnd( stats, P.data() );
  _ndxdep.insert( _Indxdep.begin(), _Indxdep.end() );
  _ndxdep.insert( _CMndxdep.begin(), _CMndxdep.end() );
#ifdef MC__CSEARCH_SHOW_REDUC
  std::cout << "Reduction #0: " << std::fixed << std::setprecision(1)
            << _reducrel( _nvar, P.data(), P0.data() )*1e2 << "%\n";
#endif
*/
  // Main loop for relaxation and domain reduction
  std::vector<T> P1;
  for( unsigned nred=0; nred < options.DOMREDMAX; nred++ ){
    P1 = P;
    if( nred )
      _update_depcheb( options, stats, P.data() );
      //_rescale_depcheb( options, stats, P.data() );
    _update_polrelax( options, stats, P.data(), tvar );
    // Solve relaxation
    //_relax( options, stats, P.data(), tvar, 0, false );
    //if( get_status() == LPRELAX_BASE<T>::LP_INFEASIBLE ) return SBB<T>::INFEASIBLE;
    // Apply domain reduction
    _contract( options, stats, P.data(), tvar, &INC );
    //_tighten_depbnd( stats, P.data() );
    if( get_status() != LPRELAX_BASE<T>::LP_OPTIMAL ) break;
    // Test improvement
    double vred = _reducrel( _nvar, P.data(), P1.data(), P0.data() );
    if( vred < options.DOMREDTHRES ) break;
#ifdef MC__CSEARCH_SHOW_REDUC
    std::cout << "Reduction #" << nred+1 << ": "
              << std::fixed << std::setprecision(1) << vred*1e2 << "%\n";
#endif
  }

#ifdef MC__CSEARCH_SHOW_BOXES
  std::cout << "\nReduced Box:\n";
  for( unsigned i=0; i<_nvar; i++ )
    //if( _exclude_vars.find(i) == _exclude_vars.end() )
      std::cout << "P[" << i << "] = " << P[i] << std::endl;
#endif

  // Update bounds for dependents
  //_update_depbnd( stats, P.data() );
  _update_depcheb( options, stats, P.data() );
  _ndxdep.insert( _CMndxdep.begin(), _CMndxdep.end() );
  _tighten_depbnd( stats, P.data() );

  // Solve relaxed problem
  _relax( options, stats, P.data(), tvar, 0, false );
  switch( get_status() ){
   case LPRELAX_BASE<T>::LP_OPTIMAL:{
    f = get_objective();
#ifdef MC__CSEARCH_DEBUG
    std::cout << std::scientific << std::setprecision(4)
              << "  _CMobj = " << _CMobj
              << "  f_rel = " << f << std::endl;
#endif
    unsigned ivar = 0;
    for( auto itv=_var.begin(); itv!=_var.end(); ++itv, ivar++ ){
      p[ivar] = get_variable( *itv );
#ifdef MC__CSEARCH_DEBUG
      std::cout << std::scientific << std::setprecision(4)
                << "  p_rel(" << ivar << ") = " << p[ivar] << std::endl;
#endif
    }
    break;
   }

   case LPRELAX_BASE<T>::LP_INFEASIBLE:
    return SBB<T>::INFEASIBLE;

   default:
    return SBB<T>::FAILURE;
  }

  return SBB<T>::NORMAL;
}

template <typename T>
template <typename SLVLOC, typename OPT, typename STAT>
inline typename SBB<T>::STATUS
CSEARCH_BASE<T>::_subproblem_local
( SLVLOC&locopt, const OPT&options, STAT&stats,
  const std::vector<T>&P, std::vector<double>&p, double&f )
{
  stats.tLOCSOL -= cpuclock();
  stats.nLOCSOL++;

#if defined (MC__CSEARCH_SHOW_BOXES)
  std::cout << "\nLocal Box:\n";
  for( unsigned i=0; i<_nvar; i++ )
    //if( _exclude_vars.find(i) == _exclude_vars.end() )
      std::cout << "P[" << i << "] = " << P[i] << std::endl;
#endif

  // Solve local optimization model
  // Set integer variable bounds fixed to relaxed solution if MIP
  std::vector<T> P_loc;
  unsigned *tvar = !_tvar.empty()? _tvar.data(): 0;
  for( unsigned ivar=0; ivar<_nvar; ivar++ )
    P_loc.push_back( (tvar && tvar[ivar])? p[ivar]: P[ivar] );
  locopt->solve( P_loc.data(), p.data() );
  const double* p_loc = _get_SLVLOC( locopt->solution() );

  // Update incumbent on return
  for( unsigned ivar=0; ivar<_nvar; ivar++ )
    p[ivar] = p_loc[ivar];
  auto flag = _subproblem_feasibility( locopt, options, p.data(), f );

  stats.tLOCSOL += cpuclock(); 
  // returning INFEASIBLE will cause B&B to fathom...
  return flag == SBB<T>::NORMAL? SBB<T>::NORMAL: SBB<T>::FAILURE;
}

template <typename T>
template <typename SLVLOC,typename OPT>
inline typename SBB<T>::STATUS
CSEARCH_BASE<T>::_subproblem_feasibility
( SLVLOC&locopt, const OPT&options, const double*p, double&f )
{
  if( !locopt->test_feasibility( p,options.FEASTOL ) )
    return SBB<T>::INFEASIBLE;
  f = locopt->solution().f;
  return SBB<T>::NORMAL;

//  // Compute objective function value
//  _dag->eval( _op_f, _wk_D, 1, std::get<1>(_obj).data(), &f, _nvar, _var.data(), p );

//  // Check current point feasibility
//  return _check_feasibility( p, options.FEASTOL )? SBB<T>::NORMAL: SBB<T>::INFEASIBLE;
}

//template <typename T>
//inline bool
//CSEARCH_BASE<T>::_check_feasibility
//( const double*p, const double FEASTOL )
//{
//  //return false;

//  // Check current point feasibility
//  auto itc=std::get<0>(_ctr).begin();
//  auto itr=std::get<3>(_ctr).begin();
//  double gj, maxinfeas=0.;
//  for( unsigned j=0; itc!=std::get<0>(_ctr).end(); ++itc, ++itr, j++ ){
//    if( *itr ) continue; // Do not check redundant constraints for feasibility
//    _dag->eval( _op_g[j], _wk_D, 1, std::get<1>(_ctr).data()+j, &gj, _nvar, _var.data(), p );
//    switch( (*itc) ){
//      case EQ: maxinfeas = std::fabs(gj); break;
//      case LE: maxinfeas = gj;            break;
//      case GE: maxinfeas = -gj;           break;
//    }
//    //std::cout << "ctr#" << j << ": " << gj << "  (" << maxinfeas << ")\n";
//    if( maxinfeas > FEASTOL ) return false;
//  }
//  return true;
//}

template <typename T>
inline void
CSEARCH_BASE<T>::_update_incumbent
( const double*p, const double&f )
{
  if( _p_inc.empty()
   || ( std::get<0>(_obj)[0]==MIN && f<_f_inc )
   || ( std::get<0>(_obj)[0]==MAX && f>_f_inc ) ){
    _p_inc.resize(_nvar);
    for( unsigned ivar=0; ivar<_nvar; ivar++ )
      _p_inc[ivar] = p[ivar];
    _f_inc = f;
  }      
}

template <typename T>
template <typename OPT, typename STAT>
inline typename SBP<T>::STATUS
CSEARCH_BASE<T>::_assess
( OPT&options, STAT&stats, SBPNode<T>*node )
{
  // Strong branching case
  unsigned DOMREDMAX = options.DOMREDMAX;
  if( node->strongbranch() && options.STGBCHDRMAX
   && options.STGBCHDRMAX < options.DOMREDMAX ){
    options.DOMREDMAX = options.STGBCHDRMAX; 
#if defined (MC__CSEARCH_DEBUG)
    std::cout << "on-going strong branching\n";
#endif
  }

  // Dependents in parent node
  _ndxdep = node->depend();

  // Node assessment and domain contraction
  auto status = _subproblem_assess( options, stats, node->P() );

  // Strong branching case
  if( node->strongbranch() ){
    options.DOMREDMAX = DOMREDMAX;
    return status;
  }

  // Update scores and dependents
  //_set_chebscores( options, node->scores(), false );
  //_set_blkscores( options, node->scores() );
  _set_branching_scores( node->scores() );
  node->depend() = _ndxdep;

#if defined (MC__CSEARCH_SHOW_OUTER)
  if( status == SBP<T>::OUTER ){
    std::cout << "\nOuter Box:\n";
    for( unsigned i=0; i<_nvar; i++ )
      std::cout << node->P()[i] << std::endl;    
    {int dum; std::cout << "PAUSED"; std::cin >> dum;}
  }
#endif
  return status;
}

template <typename T>
template <typename OPT, typename STAT>
inline typename SBP<T>::STATUS
CSEARCH_BASE<T>::_subproblem_assess
( const OPT&options, STAT&stats, std::vector<T>&P )
{
#if defined (MC__CSEARCH_SHOW_BOXES)
  std::cout << "\nInitial Box:\n";
  for( unsigned i=0; i<_nvar; i++ )
    //if( _exclude_vars.find(i) == _exclude_vars.end() ) continue;
    std::cout << P[i] << std::endl;
#endif
  unsigned *tvar = !_tvar.empty()? _tvar.data(): 0;
  std::vector<T> P0 = P;

  // Update bounds for dependents
  if( _update_depbnd( stats, P.data() ) ){
    _ndxdep.insert( _Indxdep.begin(), _Indxdep.end() );
    if( _update_depcheb( options, stats ) )
      _ndxdep.insert( _CMndxdep.begin(), _CMndxdep.end() );
    _tighten_depbnd( stats, P.data() );
  }
  else{
    _CMndxdep.clear();
  }
#ifdef MC__CSEARCH_SHOW_REDUC
  std::cout << "Reduction #0: " << std::fixed << std::setprecision(1)
            << _reducrel( _nvar, P.data(), P0.data() )*1e2 << "%\n";
#endif

  // Perform inclusion test
  _update_polrelax( options, stats, P.data(), tvar );
  auto status = _check_inclusion( options, _CMctr );
  if( status != SBP<T>::UNDETERMINED ) return status;

  // Main loop for relaxation and domain reduction
  std::vector<T> P1;
  for( unsigned nred=0; nred < options.DOMREDMAX; nred++ ){
    P1 = P;
    if( nred ){
      _update_depbnd( stats, P.data() );
      _update_depcheb( options, stats, P.data() );
      //_rescale_depcheb( options, stats, P.data() );
      _update_polrelax( options, stats, P.data(), tvar );
    }
    // Apply domain reduction
    _contract( options, stats, P.data(), tvar, &options.CTRBACKOFF, true );
    //_tighten_depbnd( stats, P.data() );
    if( get_status() != LPRELAX_BASE<T>::LP_OPTIMAL ) break;
    // Perform inclusion test
    //std::cout << "Solving inclusion test:\n"; {int dum; std::cin >> dum; }
    auto status = _check_inclusion( options, _CMctr );
    if( status != SBP<T>::UNDETERMINED ) return status;
    // Test improvement
    double vred = _reducrel( _nvar, P.data(), P1.data(), P0.data() );
#ifdef MC__CSEARCH_SHOW_REDUC
    std::cout << "Reduction #" << nred+1 << ": "
              << std::fixed << std::setprecision(1) << vred*1e2 << "%\n";
    for( unsigned i=0; i<_nvar; i++ )
      std::cout << P[i] << std::endl;
#endif
    if( vred < options.DOMREDTHRES ) break;
  }

#if defined (MC__CSEARCH_SHOW_BOXES)
  std::cout << "\nReduced Box:\n";
  for( unsigned i=0; i<_nvar; i++ )
    //if( _exclude_vars.find(i) == _exclude_vars.end() ) continue;
    std::cout << P[i] << std::endl;
#endif

  // Update bounds for dependents
  //_update_depbnd( stats, P.data() );
  if( _update_depcheb( options, stats, P.data() ) ){
    _tighten_depbnd( stats, P.data() );
    _ndxdep.insert( _CMndxdep.begin(), _CMndxdep.end() );
  }

  // Solve relaxed feasibility problem
  _update_polrelax( options, stats, P.data(), tvar, true );
#if defined (MC__CSEARCH_DEBUG)
  std::cout << "Solving relaxed feasibility problem:\n"; //{int dum; std::cin >> dum; }
#endif
  _relax( options, stats, P.data(), tvar, 0, false, true );
#if defined (MC__CSEARCH_DEBUG)
  std::cout << "Status:" << get_status() << "\n"; {int dum; std::cin >> dum; }
#endif
  switch( get_status() ){
   // objective contains minimum backoffs
   case LPRELAX_BASE<T>::LP_OPTIMAL:  
#if defined (MC__CSEARCH_SHOW_BOXES)
    std::cout << "\nFinal Box (OBJ=" << std::scientific << std::setprecision(14) << get_objective() << "):\n";
    for( unsigned i=0; i<_nvar; i++ )
      std::cout << P[i] << std::endl;
#endif
    if( get_objective() > options.FEASTOL ){//LPOPTIMTOL ){
      return SBP<T>::OUTER;
    }
    break; // need more test to find out if INNER or UNDETERMINED
   // infeasibility may not happen w/ backoffs
   case LPRELAX_BASE<T>::LP_INFEASIBLE:
   default:
     return SBP<T>::FAILURE;
  }

  return SBP<T>::UNDETERMINED;
}

template <typename T>
template <typename OPT, typename BND>
inline typename SBP<T>::STATUS 
CSEARCH_BASE<T>::_check_inclusion
( const OPT&options, const std::vector<BND>&C ) const
{
  if( C.empty() ) return SBP<T>::UNDETERMINED;

  // Function enclosure and inclusion tests
  bool flag1 = true, flag2 = true;
  auto itc = std::get<0>(_ctr).begin();
  for( unsigned j=0; flag2 && itc!=std::get<0>(_ctr).end(); ++itc, j++ ){
#if defined (MC__CSEARCH_SHOW_INCLUSION)
    std::cout << "C" << j << ":" << C[j];
    std::cout << "C" << j << " L:" << Op<BND>::l(C[j])
                          << " U:" << Op<BND>::u(C[j]) << std::endl;
#endif
    switch( (*itc) ){
      case EQ: flag1 = false;
               if( Op<BND>::u(C[j]) < -options.FEASTOL
                || Op<BND>::l(C[j]) >  options.FEASTOL )
                 flag2 = false;
               break;
      case LE: if( Op<BND>::l(C[j]) > options.FEASTOL )
                 flag2 = false;
               else if( Op<BND>::u(C[j]) > -options.FEASTOL )
                 flag1 = false;
               break;
      case GE: if( Op<BND>::u(C[j]) < -options.FEASTOL )
                 flag2 = false;
               else if( Op<BND>::l(C[j]) < options.FEASTOL )
                 flag1 = false;
               break;
    }
  }
  if (flag1 && flag2)
    return SBP<T>::INNER;
  if (!flag2)
    return SBP<T>::OUTER;
  return SBP<T>::UNDETERMINED; 
}

template <typename T>
template <typename OPT, typename STAT>
inline void
CSEARCH_BASE<T>::_refine_polrelax
( const OPT&options, STAT&stats, const T*P, const unsigned*tvar,
  const bool feastest )
{

  // Update polyhedral main variables
  stats.tPOLIMG -= cpuclock();
  for( auto itv = _POLenv.Vars().begin(); itv!=_POLenv.Vars().end(); ++itv ){
    double Xval = _get_variable( *itv->second );
    itv->second->add_breakpt( Xval );
#ifdef MC__NLGO_SHOW_BREAKPTS
    std::cout << itv->second->name();
    for( auto it = itv->second->breakpts().begin(); it!=itv->second->breakpts().end(); ++it )
      std::cout << "  " << *it;
    std::cout << std::endl;
#endif
  }
  unsigned ip = 0;
  for( auto itv = _POLvar.begin(); itv!=_POLvar.end(); ++itv, ip++ ){
    double Xval = _get_variable( *itv );
    //double Xval = _LPvar.find( &(*itv) )->second.get( GRB_DoubleAttr_X );
    itv->add_breakpt( Xval );
    itv->update( P[ip], tvar?(tvar[ip]?false:true):true );
  }

  // Update polyhedral cuts
  _POLenv.reset_cuts();
  _set_lincuts( options, feastest );
  if( options.RELMETH==OPT::DRL || options.RELMETH==OPT::HYBRID )
    _set_poacuts( options, feastest );
  if( options.RELMETH==OPT::CHEB || options.RELMETH==OPT::HYBRID )
    _set_chebcuts( options, feastest );
  stats.tPOLIMG += cpuclock();

  stats.tLPSET -= cpuclock();
  _set_LPcuts();
  stats.tLPSET += cpuclock();
}

template <typename T>
template <typename OPT, typename STAT>
inline void
CSEARCH_BASE<T>::_set_polrelax
( const OPT&options, STAT&stats, const T*P, const unsigned*tvar,
  const bool feastest )
{
#ifdef MC__NLGO_TRACE
   std::cerr << "ENTERING: _set_porelax\n";
#endif
  // Reset polynomial image
  stats.tPOLIMG -= cpuclock();
  _POLenv.reset();
  _POLvar.clear();
  _POLenv.options = options.POLIMG;

  // Set polyhedral main variables and auxiliary objective
  auto itv = _var.begin();
  for( unsigned ip=0; itv!=_var.end(); ++itv, ip++ )
    _POLvar.push_back(
      PolVar<T>( &_POLenv, *itv, P[ip], tvar?(tvar[ip]?false:true):true )
    );
  _POLobjaux.set( &_POLenv, _objvar, T(-SBB<T>::INF,SBB<T>::INF), true );

  // Add linear cuts
  _set_lincuts( options, feastest );

  if( options.RELMETH==OPT::DRL || options.RELMETH==OPT::HYBRID )
    // Add polyhedral cuts
    _set_poacuts( options, feastest );

  if( options.RELMETH==OPT::CHEB || options.RELMETH==OPT::HYBRID ){
    // Chebyshev model environment reset
    if( _CMenv && (_CMenv->nvar() != _nvar || _CMenv->maxord() != options.CMODPROP) ){
      _CMvar.clear();
      delete _CMenv; _CMenv = 0;   
    }
    if( !_CMenv ){
      // Set Chebyshev model
      _CMenv = new SCModel<T>( options.CMODPROP, _nvar );
      _CMenv->options = options.CMODEL;
      _CMvar.resize( _nvar );
    }

    // Set/update Chebyshev variables
    for( unsigned ip=0; ip<_nvar; ip++ )
      _CMvar[ip].set( _CMenv, ip, P[ip] );

    // Set polyhedral main scaled variables and size Chebyshev basis
    _POLscalvar.clear();
    itv = _scalvar.begin();
    for( unsigned ip=0; itv!=_scalvar.end(); ++itv, ip++ )
      _POLscalvar.push_back( PolVar<T>( &_POLenv, *itv, Op<T>::zeroone()*2.-1., true ) );

    // Add polyhedral cuts
    _set_chebcuts( options, feastest );
  }
  stats.tPOLIMG += cpuclock();

  // Update cuts in relaxed model
  stats.tLPSET -= cpuclock();
  _set_LPcuts();
  stats.tLPSET += cpuclock();
#ifdef MC__NLGO_TRACE
   std::cerr << "EXITING: _set_porelax\n";
#endif
}

template <typename T>
template <typename OPT, typename STAT>
inline void
CSEARCH_BASE<T>::_update_polrelax
( const OPT&options, STAT&stats, const T*P, const unsigned*tvar,
  const bool feastest )
{
#ifdef MC__NLGO_TRACE
   std::cerr << "ENTERING: _update_porelax\n";
#endif
  // Reset polyhedral cuts
  stats.tPOLIMG -= cpuclock();
  _POLenv.reset_cuts();

  // Update variable bounds
  unsigned ip = 0;
  for( auto itv = _POLvar.begin(); itv!=_POLvar.end(); ++itv, ip++ )
    itv->update( P[ip], tvar?(tvar[ip]?false:true):true );

  // Add linear cuts
  _set_lincuts( options, feastest );

  if( options.RELMETH==OPT::DRL || options.RELMETH==OPT::HYBRID )
    // Add polyhedral cuts
    _set_poacuts( options, feastest );

  if( options.RELMETH==OPT::CHEB || options.RELMETH==OPT::HYBRID ){
    // Update Chebyshev variables
    for( unsigned ip=0; ip<_nvar; ip++ )
      _CMvar[ip].set( _CMenv, ip, P[ip] );
#ifdef MC__NLGO_DEBUG_CHEBDEPS
    std::cout << "CMenv:" << *_CMenv << std::endl;
    std::cout << "CMrenv:" << *_CMrenv << std::endl;
    { int dum; std::cout << "PAUSED"; std::cin >> dum; }
#endif

    // Add Chebyshev-derived polyhedral cuts
    _set_chebcuts( options, feastest );
  }
  stats.tPOLIMG += cpuclock();

  // Update cuts in relaxed model
  stats.tLPSET -= cpuclock();
  _set_LPcuts();
  stats.tLPSET += cpuclock();
#ifdef MC__NLGO_TRACE
   std::cerr << "EXITING: _update_porelax\n";
#endif
}

template <typename T>
template <typename OPT>
inline void
CSEARCH_BASE<T>::_set_lincuts
( const OPT&options, const bool feastest )
{
  // Add polyhedral cuts for objective - by-pass if feasibility test or no objective function defined
  auto ito = std::get<0>(_obj).begin();
  for( unsigned j=0; ito!=std::get<0>(_obj).end() && j<1; ++ito, j++ ){
    if( feastest || _fct_lin.find(std::get<1>(_obj).data()) == _fct_lin.end() ) continue;
    try{
      _dag->eval( _op_f, _wk_POL, 1, std::get<1>(_obj).data(), &_POLobj, _nvar, _var.data(), _POLvar.data() );
      _POLenv.generate_cuts( 1, &_POLobj, false );
      switch( std::get<0>(_obj)[0] ){
        case MIN: _POLenv.add_cut( PolCut<T>::GE, 0., _POLobjaux, 1., _POLobj, -1. ); break;
        case MAX: _POLenv.add_cut( PolCut<T>::LE, 0., _POLobjaux, 1., _POLobj, -1. ); break;
      }
    }
    catch(...){
      // No cut added for objective in case evaluation failed
    }
  }

  // Add polyhedral cuts for constraints - add slack to inequality constraints if feasibility test
  _POLctr.resize( _nctr );
  auto itc = std::get<0>(_ctr).begin();
  for( unsigned j=0; itc!=std::get<0>(_ctr).end(); ++itc, j++ ){
    if( _fct_lin.find(std::get<1>(_ctr).data()+j) == _fct_lin.end() ) continue;
    try{
      _dag->eval( _op_g[j], _wk_POL, 1, std::get<1>(_ctr).data()+j, _POLctr.data()+j, _nvar, _var.data(), _POLvar.data() );
      _POLenv.generate_cuts( 1, _POLctr.data()+j, false );
      switch( (*itc) ){
//        case EQ: _POLenv.add_cut( PolCut<T>::EQ, 0., _POLctr[j], 1. ); break;
        case EQ: if( !feastest ){ _POLenv.add_cut( PolCut<T>::EQ, 0., _POLctr[j], 1. ); break; }
                 else             _POLenv.add_cut( PolCut<T>::GE, 0., _POLctr[j], 1., _POLobjaux,  1. ); // no break
        case LE: if( !feastest ){ _POLenv.add_cut( PolCut<T>::LE, 0., _POLctr[j], 1. ); break; }
                 else           { _POLenv.add_cut( PolCut<T>::LE, 0., _POLctr[j], 1., _POLobjaux, -1. ); break; }
        case GE: if( !feastest ){ _POLenv.add_cut( PolCut<T>::GE, 0., _POLctr[j], 1. ); break; }
                 else           { _POLenv.add_cut( PolCut<T>::GE, 0., _POLctr[j], 1., _POLobjaux,  1. ); break; }
      }
    }
    catch(...){
      // No cut added for constraint #j in case evaluation failed
      continue;
    }
  }
#ifdef MC__CSEARCH_DEBUG
  std::cout << _POLenv;
  { int dum; std::cout << "PAUSED --"; std::cin >> dum; } 
#endif
}

template <typename T>
template <typename OPT>
inline void
CSEARCH_BASE<T>::_set_poacuts
( const OPT&options, const bool feastest )
{
  // Add polyhedral cuts for objective - by-pass if feasibility test or no objective function defined
  //if( !feastest && std::get<1>(_obj).size() ){
  auto ito = std::get<0>(_obj).begin();
  for( unsigned j=0; ito!=std::get<0>(_obj).end() && j<1; ++ito, j++ ){
    if( feastest || _fct_lin.find(std::get<1>(_obj).data()) != _fct_lin.end() ) continue;
    try{
      //std::cout << "objective" << std::endl;
      _dag->eval( _op_f, _wk_POL, 1, std::get<1>(_obj).data(), &_POLobj, _nvar, _var.data(), _POLvar.data() );
      _POLenv.generate_cuts( 1, &_POLobj, false );
      switch( std::get<0>(_obj)[0] ){
        case MIN: _POLenv.add_cut( PolCut<T>::GE, 0., _POLobjaux, 1., _POLobj, -1. ); break;
        case MAX: _POLenv.add_cut( PolCut<T>::LE, 0., _POLobjaux, 1., _POLobj, -1. ); break;
      }
    }
    catch(...){
      // No cut added for objective in case evaluation failed
    }
  }

  // Add polyhedral cuts for constraints - add slack to inequality constraints if feasibility test
  _POLctr.resize( _nctr );
  auto itc = std::get<0>(_ctr).begin();
  for( unsigned j=0; itc!=std::get<0>(_ctr).end(); ++itc, j++ ){
    if( _fct_lin.find(std::get<1>(_ctr).data()+j) != _fct_lin.end() ) continue;
    try{
      //std::cout << "constraint #" << j << std::endl;
      //_dag->output( _op_g[j] );
      _dag->eval( _op_g[j], _wk_POL, 1, std::get<1>(_ctr).data()+j, _POLctr.data()+j, _nvar, _var.data(), _POLvar.data() );
      _POLenv.generate_cuts( 1, _POLctr.data()+j, false );
      switch( (*itc) ){
        case EQ: if( !feastest ){ _POLenv.add_cut( PolCut<T>::EQ, 0., _POLctr[j], 1. ); break; }
                 else             _POLenv.add_cut( PolCut<T>::GE, 0., _POLctr[j], 1., _POLobjaux,  1. ); // no break
        case LE: if( !feastest ){ _POLenv.add_cut( PolCut<T>::LE, 0., _POLctr[j], 1. ); break; }
                 else           { _POLenv.add_cut( PolCut<T>::LE, 0., _POLctr[j], 1., _POLobjaux, -1. ); break; }
        case GE: if( !feastest ){ _POLenv.add_cut( PolCut<T>::GE, 0., _POLctr[j], 1. ); break; }
                 else           { _POLenv.add_cut( PolCut<T>::GE, 0., _POLctr[j], 1., _POLobjaux,  1. ); break; }
      }
    }
    catch(...){
      // No cut added for constraint #j in case evaluation failed
      continue;
    }
  }
#ifdef MC__CSEARCH_DEBUG
  std::cout << _POLenv;
#endif
}

template <typename T>
template <typename OPT>
inline void
CSEARCH_BASE<T>::_set_chebcuts
( const OPT&options, const bool feastest )
{
  // reset bases maps
  const unsigned basisord = (!options.CMODCUTS || options.CMODCUTS>options.CMODPROP)?
                            options.CMODPROP: options.CMODCUTS;
  _basis.clear();
  _POLbasis.clear();

  // Add cuts for variable scaling and basis polynomials (up to order options.CMODCUT)
  // TEST FOR PARTICIPATING VARIABLES FIRST?!?
  // if( _var_lin.find( ip ) != _var_lin.end() ) continue; // only nonlinearly participating variables
  auto itv = _POLscalvar.begin();
  for( unsigned i=0; itv!=_POLscalvar.end(); ++itv, i++ ){
    if( _CMenv->scalvar()[i] > options.CMODEL.MIN_FACTOR )
      _POLenv.add_cut( PolCut<T>::EQ, _CMenv->refvar()[i], _POLvar[i], 1., *itv, -_CMenv->scalvar()[i] );
    //else{
    //  _POLenv.add_cut( PolCut<T>::GE, _CMenv->refvar()[i]-_CMenv->scalvar()[i], _POLvar[i], 1. );
    //  _POLenv.add_cut( PolCut<T>::LE, _CMenv->refvar()[i]+_CMenv->scalvar()[i], _POLvar[i], 1. );
    //}
#ifdef MC__NLGO_CHEBCUTS_DEBUG
      std::cout << "Variable #" << i << ": " << _CMvar[i];
#endif
  }

  // Evaluate Chebyshev model for (nonlinear) objective and keep track of participating monomials
  auto ito = std::get<0>(_obj).begin();
  for( unsigned j=0; ito!=std::get<0>(_obj).end() && j<1; ++ito, j++ ){
    if( feastest || _fct_lin.find(std::get<1>(_obj).data()) != _fct_lin.end() ) continue; // only nonlinear objectives
    T Iobj(-SBB<T>::INF,SBB<T>::INF);
    try{
      // Polynomial model evaluation
      _dag->eval( _op_f, _wk_CM, 1, std::get<1>(_obj).data(), &_CMobj,
                  _nvar, _var.data(), _CMvar.data() );
#ifdef MC__NLGO_CHEBCUTS_DEBUG
      std::cout << "Objective remainder: " << _CMobj.R() << std::endl;
      std::cout << _CMobj;
      //if(_p_inc.data()) std::cout << _CMobj.P(_p_inc.data())+_CMobj.R() << std::endl;
#endif
      // Test for too large Chebyshev bounds or NaN
      if( !(Op<SCVar<T>>::diam(_CMobj) <= options.CMODDMAX) ) throw(0);

      // Keep track of participating Chebyshev monomials into '_basis'
      auto first_CMobj = _CMobj.coefmon().lower_bound( std::make_pair( 1, std::map<unsigned,unsigned>() ) );
      auto last_CMobj  = _CMobj.coefmon().lower_bound( std::make_pair( basisord+1, std::map<unsigned,unsigned>() ) );
      for( auto it=first_CMobj; it!=last_CMobj; ++it )
        _basis.insert( std::make_pair( it->first, FFVar() ) );
    }
    catch( int ecode ){
#ifdef MC__NLGO_CHEBCUTS_DEBUG
      std::cout << "Objective: " << "cut too weak!\n";
#endif
      Op<T>::inter( Iobj, Iobj, _CMobj.B() );
      _CMobj = Iobj;
    }
    catch(...){
#ifdef MC__NLGO_CHEBCUTS_DEBUG
      std::cout << "Objective: " << "cut failed!\n";
#endif
      //Op<T>::inter( Iobj, Iobj, _CMobj.B() );
      _CMobj = Iobj;
      // No cut added for objective in case evaluation failed
    }
  }

  // Evaluate Chebyshev model for (nonlinear) constraints and keep track of participating monomials
  _CMctr.resize( _nctr );
  auto itc = std::get<0>(_ctr).begin();
  for( unsigned j=0; itc!=std::get<0>(_ctr).end(); ++itc, j++ ){
    if( _fct_lin.find(std::get<1>(_ctr).data()+j) != _fct_lin.end() ) continue;
    T Ictr(-SBB<T>::INF,SBB<T>::INF);
    try{
      // Polynomial model evaluation
#ifdef MC__NLGO_CHEBCUTS_DEBUG
      _dag->output( _op_g[j] );
#endif
      _dag->eval( _op_g[j], _wk_CM, 1, std::get<1>(_ctr).data()+j, _CMctr.data()+j,
                  _nvar, _var.data(), _CMvar.data() );
#ifdef MC__NLGO_CHEBCUTS_DEBUG
      std::cout << "\nConstraint #" << j << ": " << _CMctr[j];
      //if(_p_inc.data()) std::cout << "optim: " << _CMctr[j].P(_p_inc.data())+_CMctr[j].R() << std::endl;
      //{ int dum; std::cout << "PAUSED"; std::cin >> dum; }
#endif
      // Test for too large Chebyshev bounds or NaN
      if( !(Op<SCVar<T>>::diam(_CMctr[j]) <= options.CMODDMAX) ) throw(0);

      // Keep track of participating Chebyshev monomials into '_basis'
      auto first_CMctr = _CMctr[j].coefmon().lower_bound( std::make_pair( 1, std::map<unsigned,unsigned>() ) );
      auto last_CMctr  = _CMctr[j].coefmon().lower_bound( std::make_pair( basisord+1, std::map<unsigned,unsigned>() ) );
      for( auto it=first_CMctr; it!=last_CMctr; ++it )
        _basis.insert( std::make_pair( it->first, FFVar() ) );
   }
    catch( int ecode ){
#ifdef MC__NLGO_CHEBCUTS_DEBUG
      std::cout << "Constraint #" << j << ": " << "cut too weak!\n";
#endif
      Op<T>::inter( Ictr, Ictr, _CMctr[j].B() );
      _CMctr[j] = Ictr;
    }
    catch(...){
#ifdef MC__NLGO_CHEBCUTS_DEBUG
      std::cout << "Constraint #" << j << ": " << "cut failed!\n";
#endif
      //Op<T>::inter( Ictr, Ictr, _CMctr[j].B() );
      _CMctr[j] = Ictr;
      // No cut added for constraint #j in case evaluation failed
    }
  }

  // Keep track of participating monomials in Chebyshev model for dependents
  for( auto it=_CMndxdep.begin(); options.CMODDEPS && options.CMODRED == OPT::APPEND && it!=_CMndxdep.end(); ++it ){
    const unsigned irdep = _var_fperm[*it]-_nrvar;
    T Irdep(-SBB<T>::INF,SBB<T>::INF);
    try{
      // Test for too large Chebyshev bounds or NaN
      if( !(Op<SCVar<T>>::diam(_CMrdep[irdep]) <= options.CMODDMAX) ) throw(0);
      // Keep track of participating Chebyshev monomials into '_basis'
      auto first_CMrdep = _CMrdep[irdep].coefmon().lower_bound( std::make_pair( 1, std::map<unsigned,unsigned>() ) );
      auto last_CMrdep  = _CMrdep[irdep].coefmon().lower_bound( std::make_pair( basisord+1, std::map<unsigned,unsigned>() ) );
      for( auto jt=first_CMrdep; jt!=last_CMrdep; ++jt ){
        t_expmon expmon_CMrdep;
        expmon_CMrdep.first = jt->first.first;
        // Match indices between the full and reduced Chebyshev models
        for( auto ie = jt->first.second.begin(); ie!=jt->first.second.end(); ++ie )
          expmon_CMrdep.second.insert( std::make_pair( _var_rperm[ie->first], ie->second ) );
        _basis.insert( std::make_pair( expmon_CMrdep, FFVar() ) );
      }
    }
    catch(...){
#ifdef MC__NLGO_CHEBCUTS_DEBUG
      std::cout << "Dependent #" << irdep << ": " << "cut failed!\n";
#endif
      Op<T>::inter( Irdep, Irdep, _CMrdep[irdep].B() );
      _CMrdep[irdep] = Irdep; // <- more efficient to interesect with dependent bounds?
      // No cut added for constraint #j in case evaluation failed
    }

  }

  // Populate '_basis' with references to the Chebyshev monomials in the DAG
  _CMenv->get_bndmon( _basis, _scalvar.data(), true );
#ifdef MC__NLGO_CHEBCUTS_DEBUG
  std::cout << "\nBASIS:\n";
  for( auto it=_basis.begin(); it!=_basis.end(); ++it ){
    std::cout << it->second << " = ";
    for( auto ie=it->first.second.begin(); ie!=it->first.second.end(); ++ie )
      std::cout << "T" << ie->second << "[" << ie->first << "] ";
    std::cout << std::endl;
  }
  std::list<const mc::FFOp*> op_basis  = _dag->subgraph( _basis );
  _dag->output( op_basis );
  { int dum; std::cin >> dum; }
#endif

  // Populate '_POLbasis' with references to the Chebyshev monomials in the polyhedral relaxation
  _dag->eval( _wk_POL, _basis, _POLbasis, _scalvar.size(), _scalvar.data(), _POLscalvar.data() );

  // Append cuts for the Chebyshev monomials in the polyhedral relaxation
  _POLenv.generate_cuts( _POLbasis, false );

  // Add Chebyshev-derived cuts for objective
  ito = std::get<0>(_obj).begin();
  for( unsigned j=0; ito!=std::get<0>(_obj).end() && j<1; ++ito, j++ ){
    if( feastest
     || _fct_lin.find(std::get<1>(_obj).data()) != _fct_lin.end()
     || !(Op<SCVar<T>>::diam(_CMobj) < SBB<T>::INF) ) continue; // only nonlinear / finite objectives

    // Constant, variable and bound on objective model
    T Robj = _CMobj.bndord( basisord+1 ) + _CMobj.remainder();
    const double a0 = (_CMobj.coefmon().empty() || _CMobj.coefmon().begin()->first.first)?
                      0.: _CMobj.coefmon().begin()->second;
    auto first_CMobj = _CMobj.coefmon().lower_bound( std::make_pair( 1, std::map<unsigned,unsigned>() ) );
    auto last_CMobj  = _CMobj.coefmon().lower_bound( std::make_pair( basisord+1, std::map<unsigned,unsigned>() ) );
    t_coefmon Cobj; Cobj.insert( first_CMobj, last_CMobj );

    // Linear objective sparse cut
    switch( std::get<0>(_obj)[0] ){
      case MIN: _POLenv.add_cut( PolCut<T>::LE, -Op<T>::l(Robj)-a0, _POLbasis, Cobj, _POLobjaux, -1. ); break;
      case MAX: _POLenv.add_cut( PolCut<T>::GE, -Op<T>::u(Robj)-a0, _POLbasis, Cobj, _POLobjaux, -1. ); break;
    }
    _POLobjaux.update( _CMobj.bound(), true );

  }
#ifdef MC__NLGO_DEBUG
  std::cout << _POLenv;
#endif

  // Add Chebyshev-derived cuts for constraints
  itc = std::get<0>(_ctr).begin();
  for( unsigned j=0; itc!=std::get<0>(_ctr).end(); ++itc, j++ ){
    if( _fct_lin.find(std::get<1>(_ctr).data()+j) != _fct_lin.end()
     || !(Op<SCVar<T>>::diam(_CMctr[j]) < SBB<T>::INF) ) continue; // only nonlinear / finite constraints

    // Constant, variable and bound on constraint model
    T Rctr = _CMctr[j].bndord( basisord+1 ) + _CMctr[j].remainder();
    const double a0 = (_CMctr[j].coefmon().empty() || _CMctr[j].coefmon().begin()->first.first)?
                      0.: _CMctr[j].coefmon().begin()->second;
    auto first_CMctr = _CMctr[j].coefmon().lower_bound( std::make_pair( 1, std::map<unsigned,unsigned>() ) );
    auto last_CMctr  = _CMctr[j].coefmon().lower_bound( std::make_pair( basisord+1, std::map<unsigned,unsigned>() ) );
    t_coefmon Cctr; Cctr.insert( first_CMctr, last_CMctr );

    // Nonlinear constraint sparse cut
    switch( (*itc) ){
      case EQ: if( !feastest ){ if( !Cctr.empty() ) _POLenv.add_cut( PolCut<T>::GE, -Op<T>::u(Rctr)-a0, _POLbasis, Cctr ); }// no break
               else { _POLenv.add_cut( PolCut<T>::GE, -Op<T>::u(Rctr)-a0, _POLbasis, Cctr, _POLobjaux,  1. ); }// no break
      case LE: if( !feastest ){ if( !Cctr.empty() ) _POLenv.add_cut( PolCut<T>::LE, -Op<T>::l(Rctr)-a0, _POLbasis, Cctr ); break; }
               else { _POLenv.add_cut( PolCut<T>::LE, -Op<T>::l(Rctr)-a0, _POLbasis, Cctr, _POLobjaux, -1. ); break; }
      case GE: if( !feastest ){ if( !Cctr.empty() ) _POLenv.add_cut( PolCut<T>::GE, -Op<T>::u(Rctr)-a0, _POLbasis, Cctr ); break; }
               else { _POLenv.add_cut( PolCut<T>::GE, -Op<T>::u(Rctr)-a0, _POLbasis, Cctr, _POLobjaux,  1. ); break; }
    }
  }

  // Add Chebyshev-derived cuts for dependents
  for( auto it=_CMndxdep.begin(); options.CMODDEPS && options.CMODRED == OPT::APPEND && it!=_CMndxdep.end(); ++it ){
    // Constant, variable and bound on constraint model
    const unsigned irdep = _var_fperm[*it]-_nrvar;
#ifdef MC__NLGO_CHEBCUTS_DEBUG
    std::cout << "\nDependent #" << irdep << ": " << _CMrdep[irdep];
#endif
    T Rrdep = _CMrdep[irdep].bndord( basisord+1 ) + _CMrdep[irdep].remainder();
    const double a0 = (_CMrdep[irdep].coefmon().empty() || _CMrdep[irdep].coefmon().begin()->first.first)?
                      0.: _CMrdep[irdep].coefmon().begin()->second;
    auto first_CMrdep = _CMrdep[irdep].coefmon().lower_bound( std::make_pair( 1, std::map<unsigned,unsigned>() ) );
    auto last_CMrdep  = _CMrdep[irdep].coefmon().lower_bound( std::make_pair( basisord+1, std::map<unsigned,unsigned>() ) );
    t_coefmon Crdep; //Crdep.insert( first_CMrdep, last_CMrdep );
    for( auto jt=first_CMrdep; jt!=last_CMrdep; ++jt ){
      t_expmon expmon_CMrdep;
      expmon_CMrdep.first = jt->first.first;
      // Match indices between the full and reduced Chebyshev models
      for( auto ie = jt->first.second.begin(); ie!=jt->first.second.end(); ++ie )
        expmon_CMrdep.second.insert( std::make_pair( _var_rperm[ie->first], ie->second ) );
      Crdep.insert( std::make_pair( expmon_CMrdep, jt->second ) );
    }

    // Dependent sparse cut
    if( !feastest ){
      _POLenv.add_cut( PolCut<T>::GE, -Op<T>::u(Rrdep)-a0, _POLbasis, Crdep, _POLvar[*it],  -1. );
      _POLenv.add_cut( PolCut<T>::LE, -Op<T>::l(Rrdep)-a0, _POLbasis, Crdep, _POLvar[*it],  -1. );
    }
    else{
      _POLenv.add_cut( PolCut<T>::GE, -Op<T>::u(Rrdep)-a0, _POLbasis, Crdep, _POLvar[*it],  -1., _POLobjaux,  1. );
      _POLenv.add_cut( PolCut<T>::LE, -Op<T>::l(Rrdep)-a0, _POLbasis, Crdep, _POLvar[*it],  -1., _POLobjaux, -1. );
    }
  }

#ifdef MC__NLGO_CHEBCUTS_DEBUG
  std::cout << _POLenv;
  { int dum; std::cin >> dum; }
#endif
}

//template <typename T>
//inline void
//CSEARCH_BASE<T>::_add_chebscores
//( std::map<unsigned,double>&scores, const FFVar&FV,
//  const SCVar<T>&SCV, const double weight, const bool linexcl )
//const
//{
//  if( weight <= 0. ) return;

//  // Contribution of polynomial model remainder
//  double rem = Op<T>::diam(SCV.R()), bnd = Op<T>::diam(SCV.B());
//  //double contrib = rem * std::fabs( weight ) / FV.dep().dep().size();
//  double contrib = rem < SBB<T>::INF? rem / bnd: 1.;
//  contrib *= std::fabs( weight ) / FV.dep().dep().size();
//  for( unsigned i=0; i<_nvar; i++ ){
//    auto ft = FV.dep().dep().find( _var[i].id().second );
//    if( ft == FV.dep().dep().end() || (linexcl && ft->second) ) continue;
//    auto its = scores.insert( std::make_pair( i, contrib ) );
//    if( !its.second ) its.first->second += contrib;
//    //if( !its.second && its.first->second < contrib )
//    //  its.first->second = contrib;

//  }
//  //if( rem >= SBB<T>::INF ) return;

//  // Polynomial model nonlinear dependencies contribution
//  auto kt = SCV.coefmon().lower_bound( std::make_pair( 2 ,std::map<unsigned,unsigned>() ) );
//  for( ; kt != SCV.coefmon().cend(); ++kt ){
//    //contrib = std::fabs( weight * kt->second );
//    contrib = std::fabs( weight * kt->second ) / bnd;
//    for( auto ie=kt->first.second.begin(); ie!=kt->first.second.end(); ++ie ){
//      auto its = scores.insert( std::make_pair( ie->first, contrib ) );
//      if( !its.second ) its.first->second += contrib;
//      //if( !its.second && its.first->second < contrib )
//      //  its.first->second = contrib;
//    }
//  }
//}

//template <typename T>
//template <typename OPT>
//inline void
//CSEARCH_BASE<T>::_set_blkscores
//( const OPT&options, std::map<unsigned,double>&scores )
//{
//  scores.clear();
//  if( !options.SCOBCHUSE ) return;
//  _get_blkscores( scores );
//}

//template <typename T>
//template <typename OPT>
//inline void
//CSEARCH_BASE<T>::_set_chebscores
//( const OPT&options, std::map<unsigned,double>&scores,
//  const bool feastest )
//{
//  scores.clear();
//  if( !options.SCOBCHUSE ) return;
//  if( options.RELMETH != OPT::CHEB && options.RELMETH != OPT::HYBRID ) return;

//  // Append scores based on Chebyshev model for (nonlinear) objective
//  auto ito = std::get<0>(_obj).begin();
//  for( unsigned j=0; ito!=std::get<0>(_obj).end() && j<1; ++ito, j++ ){
//    if( feastest || _fct_lin.find(std::get<1>(_obj).data()) != _fct_lin.end() ) continue; // only nonlinear objectives
//    _add_chebscores( scores, std::get<1>(_obj)[0], _CMobj, 1. );
//  }

//  // Append scores based on Chebyshev model for (nonlinear) constraints
//  auto itc = std::get<0>(_ctr).begin();
//  for( unsigned j=0; itc!=std::get<0>(_ctr).end(); ++itc, j++ ){
//    if( _fct_lin.find(std::get<1>(_ctr).data()+j) != _fct_lin.end()
//     || std::get<3>(_ctr)[j] ) continue;
//    switch( (*itc) ){ // only score equality and active inequality constraints
//      case LE: if( Op<SCVar<T>>::u(_CMctr[j]) < 0. ) continue;
//               break;
//      case GE: if( Op<SCVar<T>>::l(_CMctr[j]) > 0. ) continue;
//               break;
//      default: break;
//    }
//    _add_chebscores( scores, std::get<1>(_ctr)[j], _CMctr[j], 1. );
//  }

//#ifdef MC__CSEARCH_SHOW_SCORES
//  for( auto it=scores.cbegin(); it!=scores.cend(); ++it )
//    std::cout << std::right << std::setw(4) << it->first << " " << it->second << std::endl; 
//  { int dum; std::cout << "PAUSED"; std::cin >> dum; }
//#endif
//}

template <typename T>
template <typename OPT, typename STAT>
inline bool
CSEARCH_BASE<T>::_set_depcheb
( const OPT&options, STAT&stats, const T*P )
{
  if( !_nrdep || !options.CMODDEPS ){
    _CMrvar.clear();
    _CMrdep.clear();
    _CMndxdep.clear();
    delete _CMrenv; _CMrenv = 0;
    return true;
  }

  // Set-up reduced-space Chebyshev model
  stats.tDEPBND -= cpuclock();
  if( _CMrenv && (_CMrenv->nvar() != _nrvar || _CMrenv->maxord() != options.CMODDEPS) ){
    _CMrvar.clear();
    _CMrdep.clear();
    delete _CMrenv; _CMrenv = 0;
  }

  if( !_CMrenv ){
    // Create reduce-space Chebyshev model
    _CMrenv = new SCModel<T>( options.CMODDEPS, _nrvar );
    _CMrenv->options = options.CMODEL;
    _CMrvar.resize( _nrvar );
    _CMrdep.resize( _nrdep );
  }

  // Update reduced Chebyshev variables
  for( unsigned ip=0; ip<_nrvar; ip++ ){
    _CMrvar[ip].set( _CMrenv, ip, P? P[_var_rperm[ip]]: _Irvar[ip] );
#ifdef MC__CSEARCH_DEBUG
    std::cout << "SCVar #" << ip << ": " << _CMrenv->bndvar()[ip] << " "
              << _CMrenv->refvar()[ip] << " " << _CMrenv->scalvar()[ip] << " " << std::endl;
#endif
  }
  for( unsigned ip=0; ip<_nrdep; ip++ )
    _CMrdep[ip] = P? P[_var_rperm[_nrvar+ip]]: _Irdep[ip];

  // Get dependent Chebyshev variables
  const bool feas = _get_depcheb();
  stats.nDEPBND++;
  stats.tDEPBND += cpuclock();
  return feas;
}

template <typename T>
template <typename OPT, typename STAT>
inline bool
CSEARCH_BASE<T>::_update_depcheb
( const OPT&options, STAT&stats, const T*P )
{
  if( _ignore_deps || !_nrdep || !options.CMODDEPS ){
    _CMndxdep.clear();
    return true;
  }

  // Update reduced Chebyshev variables
  stats.tDEPBND -= cpuclock();
  for( unsigned ip=0; ip<_nrvar; ip++ ){
    _CMrvar[ip].set( _CMrenv, ip, P? P[_var_rperm[ip]]: _Irvar[ip] );
#ifdef MC__CSEARCH_DEBUG
    std::cout << "Var #" << ip << ": " << _CMrenv->bndvar()[ip] << " "
              << _CMrenv->refvar()[ip] << " " << _CMrenv->scalvar()[ip] << " " << std::endl;
#endif
  }
  for( unsigned ip=0; ip<_nrdep; ip++ )
    _CMrdep[ip] = P? P[_var_rperm[_nrvar+ip]]: _Irdep[ip];

  // Get dependent Chebyshev variables
  const bool feas = _get_depcheb();
  stats.nDEPBND++;
  stats.tDEPBND += cpuclock();
  return feas;
}

template <typename T>
template <typename OPT, typename STAT>
inline void
CSEARCH_BASE<T>::_rescale_depcheb
( const OPT&options, STAT&stats, const T*P )
{
  if( _ignore_deps || !_nrdep || !P || !options.CMODDEPS ){
    _CMndxdep.clear();
    return;
  }

  // Rescale dependent Chebyshev variables
  stats.tDEPBND -= cpuclock();
  for( unsigned iv=0; iv<_nrvar; ++iv ){
    for( auto it=_CMndxdep.begin(); it!=_CMndxdep.end(); ++it ){
      const unsigned ip = _var_fperm[*it]-_nrvar;
      _CMrdep[ip].scale( iv, P[_var_rperm[iv]] ); //.simplify();
#ifdef MC__CSEARCH_DEBUG_DEPCHEB
      if( iv+1 == _nrvar )
        std::cout << "CMrdep[" << ip << "] =" << _CMrdep[ip];
#endif
    }
  }
  for( unsigned iv=0; iv<_nrvar; ++iv )
    _CMrvar[iv].set( _CMrenv, iv, P[_var_rperm[iv]] );
  stats.tDEPBND += cpuclock();
}

template <typename T>
template <typename STAT>
inline bool
CSEARCH_BASE<T>::_set_depbnd
( STAT&stats, const T*P )
{
  if( !_nrdep ){
    _Indxdep.clear();
    return true;
  }

  // Set-up reduced-space interval bounds
  stats.tDEPBND -= cpuclock();
  _Irvar.resize( _nrvar );
  _Irdep.resize( _nrdep );

  // Update reduced Chebyshev variables
  for( unsigned ip=0; ip<_nrvar; ip++ )
    _Irvar[ip] = P[_var_rperm[ip]];
  for( unsigned ip=0; ip<_nrdep; ip++ )
    _Irdep[ip] = P[_var_rperm[_nrvar+ip]];

  // Get dependent interval bounds
  const bool feas = _get_depbnd();
  stats.tDEPBND += cpuclock();
  return feas;
}

template <typename T>
template <typename STAT>
inline bool
CSEARCH_BASE<T>::_update_depbnd
( STAT&stats, const T*P )
{
  if( !_nrdep ){
    _Indxdep.clear();
    return true;
  }

  // Update reduced Chebyshev variables
  stats.tDEPBND -= cpuclock();
  for( unsigned ip=0; ip<_nrvar; ip++ )
    _Irvar[ip] = P[_var_rperm[ip]];
  for( unsigned ip=0; ip<_nrdep; ip++ )   // <- necessary?
    _Irdep[ip] = P[_var_rperm[_nrvar+ip]];

  // Get dependent interval bounds
  const bool feas = _get_depbnd();
  stats.tDEPBND += cpuclock();
  return feas;
}

template <typename T>
template <typename STAT>
inline void
CSEARCH_BASE<T>::_tighten_depbnd
( STAT&stats, T*P )
{
  stats.tDEPBND -= cpuclock();

  // Update dependent bounds based on implicit interval contractor
  for( auto it=_Indxdep.begin(); it!=_Indxdep.end(); ++it ){
    const unsigned idep = _var_fperm[*it]-_nrvar;
    if( !Op<T>::inter( P[*it], P[*it], _Irdep[idep] ) ){
  //for( unsigned i=0; _Irdep.size() && i<_nrdep; i++ ){   // <- maybe not necessary?
    //if( !Op<T>::inter( P[_var_rperm[_nrvar+i]], P[_var_rperm[_nrvar+i]], _Irdep[i] ) ){
      //P[_nrvar+i] = _Irdep[i]; // we may need to revisit this (e.g., AEBND failure?)
//#ifdef MC__NLGO_DEBUG_CHEBDEPS
      std::cout << "IA intersection with dependent #" << idep << " failed!\n";
//#endif
      continue;
    }
  }

//  // Update dependent bounds based on implicit Chebyshev model contractor
//  for( auto it=_CMndxdep.begin(); it!=_CMndxdep.end(); ++it ){
//    const unsigned idep = _var_fperm[*it]-_nrvar;
//    if( !Op<T>::inter( P[*it], P[*it], _CMrdep[idep].B() ) ){
//  //for( unsigned i=0; _CMrdep.size() && i<_nrdep; i++ ){   // <- maybe not necessary?
//    //if( !Op<T>::inter( P[_var_rperm[_nrvar+i]], P[_var_rperm[_nrvar+i]], _CMrdep[i].B() ) ){
//      //P[_nrvar+i] = _CMrdep[i].B(); // we may need to revisit this (e.g., AEBND failure?)
//#ifdef MC__NLGO_DEBUG_CHEBDEPS
//      std::cout << "PM intersection with dependent #" << idep << " failed!\n";
//#endif
//      continue;
//    }
//  }

  stats.tDEPBND += cpuclock();
}

template <typename T>
template <typename OPT>
inline void
CSEARCH_BASE<T>::_set_SBBoptions
( const OPT&options )
{
  // SBB options
  SBB<T>::options.STOPPING_ABSTOL              = options.CVATOL;
  SBB<T>::options.STOPPING_RELTOL              = options.CVRTOL;
  SBB<T>::options.BRANCHING_STRATEGY           = options.BRANCHPT;
  //SBB<T>::options.BRANCHING_BOUND_THRESHOLD    = options.BRANCHDMIN;
  SBB<T>::options.BRANCHING_VARIABLE_CRITERION = options.BRANCHVAR;
  SBB<T>::options.BRANCHING_USERFUNCTION       = options.BRANCHSEL;
  SBB<T>::options.SCORE_BRANCHING_USE          = options.SCOBCHMETH;
  SBB<T>::options.SCORE_BRANCHING_MAXSIZE      = options.SCOBCHVMAX;
  SBB<T>::options.SCORE_BRANCHING_RELTOL       = options.SCOBCHRTOL;
  SBB<T>::options.SCORE_BRANCHING_ABSTOL       = options.SCOBCHATOL;
  SBB<T>::options.STRONG_BRANCHING_MAXDEPTH    = options.STGBCHDEPTH;
  SBB<T>::options.STRONG_BRANCHING_WEIGHT      = options.STGBCHWEIGHT;
  SBB<T>::options.MAX_NODES                    = options.MAXITER;
  SBB<T>::options.MAX_CPUTIME                  = options.MAXCPU;
  SBB<T>::options.DISPLAY                      = options.DISPLAY;
  // TREE_REREDUCE(true)
}

template <typename T>
template <typename OPT>
inline void
CSEARCH_BASE<T>::_set_SBPoptions
( const OPT&options ) const
{
  // SBP options
  SBP<T>::options.MEASURE                      = options.NODEMEAS;
  SBP<T>::options.STOPPING_TOLERANCE           = options.CVTOL;
  SBP<T>::options.BRANCHING_VARIABLE_CRITERION = options.BRANCHVAR;
  SBP<T>::options.BRANCHING_USERFUNCTION       = options.BRANCHSEL;
  SBP<T>::options.SCORE_BRANCHING_USE          = options.SCOBCHMETH;
  SBP<T>::options.SCORE_BRANCHING_MAXSIZE      = options.SCOBCHVMAX;
  SBP<T>::options.SCORE_BRANCHING_RELTOL       = options.SCOBCHRTOL;
  SBP<T>::options.SCORE_BRANCHING_ABSTOL       = options.SCOBCHATOL;
  SBP<T>::options.STRONG_BRANCHING_MAXDEPTH    = options.STGBCHDEPTH;
  SBP<T>::options.STRONG_BRANCHING_RELTOL      = options.STGBCHRTOL;
  SBP<T>::options.STRONG_BRANCHING_ABSTOL      = options.STGBCHATOL;
  SBP<T>::options.MAX_NODES                    = options.MAXITER;
  SBP<T>::options.MAX_CPUTIME                  = options.MAXCPU;
  SBP<T>::options.DISPLAY                      = options.DISPLAY;
}

template <typename T>
template <typename OPT>
inline void
CSEARCH_BASE<T>::_display_add
( const OPT&options, const double dval )
{
  if( options.DISPLAY <= 1 ) return;
  _odisp << std::right << std::scientific << std::setprecision(_DPREC)
         << std::setw(_DPREC+8) << dval;
}

template <typename T>
template <typename OPT>
inline void
CSEARCH_BASE<T>::_display_add
( const OPT&options, const int ival )
{
  if( options.DISPLAY <= 1 ) return;
  _odisp << std::right << std::setw(_IPREC) << ival;
}

template <typename T>
template <typename OPT>
inline void
CSEARCH_BASE<T>::_display_add
( const OPT&options, const std::string &sval )
{
  if( options.DISPLAY <= 1 ) return;
  _odisp << std::right << std::setw(_DPREC+8) << sval;
}

template <typename T>
inline void
CSEARCH_BASE<T>::_display_flush
( std::ostream &os )
{
  if( _odisp.str() == "" ) return;
  //if( options.DISPLAY > 0 ){
  os << _odisp.str() << std::endl;
  //}
  _odisp.str("");
  return;
}

} // end namescape mc

#endif
