// Copyright (C) 2012-2015 Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

#ifndef MC__SBB_HPP
#define MC__SBB_HPP

#include <utility>
#include <set>
#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <cmath>

#include "mcop.hpp"
#include "mcfunc.hpp"

#undef DEBUG__SBB_PSEUDOCOST
#undef DEBUG__SBB_STRONGBRANCHING
#undef DEBUG__SBB_BRANCHINGVAR
#undef USE__SBB_PSEUDOCOST_FILTER
#undef SAVE__SBB_NODES_TO_FILE

// TO DO:
// - [DONE] ACCOMODATE MIN AND MAX PROBLEMS
// - [DONE] INTRODUCE MAXIMUM EXECUTION TIME
// - [DONE] IMPLEMENT STRONG AND RELIABLE BRANCHING IN SBB
// - [DONE] IMPLEMENT DISPLAY IN SBB
// - DOCUMENTATION!!

namespace mc
{

template <typename T> class SBBNode;
template <typename T> struct lt_SBBNode;

//! @brief Pure virtual base class for spatial branch-and-bound search
////////////////////////////////////////////////////////////////////////
//! mc::SBB<T> is a pure virutal base C++ class implementing the
//! spatial branch-and-bound (a.k.a. branch-and-reduce) algorithm
//! for nonconvex, continuous optimization problems.
//!
//! REFERENCES:
//! - Horst, R., and H. Tuy, <A href="http://books.google.com/books?id=usFjGFvuBDEC&lpg=PA39&ots=DmNqBtTNOQ&dq=global%20optimization%20horst%20tuy&pg=PP1#v=onepage&q&f=false"><i>Global Optimization</i></A>, Springer-Verlag, Berlin, 1996.
//! - Ryoo, H.S., and N.V. Sahinidis, <A href="http://dx.doi.org/10.1016/0098-1354(94)00097-8">Global optimization of nonconvex NLPs and MINLPs with application in process design</A>, <i>Computers & Chemical Engineering</i>, <b>19</b>(5):551--566, 1995.
//! - Tawarmalani, M., and N. V. Sahinidis, <A HREF="http://books.google.com/books?id=MjueCVdGZfoC&lpg=PP1&ots=bgaQ0uqGU_&dq=tawarmalani%20sahinidis&pg=PP1#v=onepage&q&f=false">Convexification and Global Optimization in Continuous and Mixed-Integer Nonlinear Programming: Theory, Algorithms, Software, and Applications</A>, Kluwer Academic Publishers, Dordrecht, Vol. 65 in ``Nonconvex Optimization And Its Applications'' series, 2002.
//! .
////////////////////////////////////////////////////////////////////////
template <typename T>
class SBB
////////////////////////////////////////////////////////////////////////
{
  template <typename U> friend class SBBNode;

public:  
  typedef std::multiset< SBBNode<T>*, lt_SBBNode<T> > t_Nodes;
  typedef typename t_Nodes::iterator it_Nodes;
  template <typename U> friend std::ostream& operator<<
    ( std::ostream&, const SBBNode<U>& );

  //! @brief Status of subproblems
  enum STATUS{
    NORMAL=0,	//!< Normal termination
    INFEASIBLE,	//!< Termination w/ indication of infeasibility
    FAILURE,	//!< Termination after failure
    FATAL	//!< Termination after fatal error
  };

  //! @brief Task for subproblems
  enum TASK{
    LOWERBD=0,	//!< Solve lower bounding subproblem
    POSTPROC,	//!< Postprocess node (e.g. domain reduction)
    PREPROC,	//!< Preprocess node (e.g. constraint propagation)
    UPPERBD,	//!< Solve upper bounding subproblem
    FEASTEST	//!< Test feasibility of upper bounding subproblem
  };
  //! @brief Problem type
  enum TYPE{
    MINIMIZE=0,	//!< Minimize problem
    MAXIMIZE	//!< Maximize problem
  };

  //! @brief Prototype for user-supplied function
  virtual STATUS subproblems
    ( const TASK,		// task
      const unsigned, 	// number of variables
      T*,			// variable bounds
      double*,			// optimal solution point - initial guess on entry
      double&, 			// optimal solution value
      const double		// current incumbent
    )=0;

  //! @brief Public constructor
  SBB()
    : options(), _p_inc(0), _f_inc(options.INF), _np(0), _P_root(0),
      _P_tmp(0), _node_inc(0)
    {}

  //! @brief Class destructor
  virtual ~SBB()
    {
      _clean_stack();
      delete[] _P_root;
      delete[] _P_tmp;
      delete[] _p_inc;
    }

  //! @brief SBB options
  struct Options
  {
    //! @brief Constructor
    Options():
      ABSOLUTE_TOLERANCE(1e-3), RELATIVE_TOLERANCE(1e-3), TREE_REREDUCE(true),
      BRANCHING_STRATEGY(OMEGA), BRANCHING_VARIABLE_CRITERION(RGREL),
      BRANCHING_RELIABILITY_THRESHOLD(1), BRANCHING_SCORE_FACTOR(1./6.),
      BRANCHING_BOUND_THRESHOLD(1e-3), DISPLAY(2), MAX_CPU_TIME(1e6),
      MAX_NODES(0), INF(1e20)
      {}
    //! @brief Display
    void display
      ( std::ostream&out ) const;
    //! @brief Branching strategy
    enum STRATEGY{
      MIDPOINT=0,	//!< Bisection at mid-point
      OMEGA		//!< Bisection at (i) incumbent location if inside current bounds; otherwise (ii) at lower-bounding-problem solution if not at bound; otherwise (iii) at mid-point
    };
    //! @brief Branching variable criterion
    enum CRITERION{
      RGREL=0,	//!< Relative variable range diameter
      RGABS,	//!< Absolute variable range diameter
      RELIAB	//!< Reliability branching based on pseudocosts and strong branching
    };
    //! @brief Absolute stopping tolerance
    double ABSOLUTE_TOLERANCE;
    //! @brief Relative stopping tolerance
    double RELATIVE_TOLERANCE;
    //! @brief Whether the entire B&B tree is to be updated after every incumbent update
    bool TREE_REREDUCE;
    //! @brief Branching strategy
    STRATEGY BRANCHING_STRATEGY;
    //! @brief Variable selection criterion
    CRITERION BRANCHING_VARIABLE_CRITERION;
    //! @brief Minimum number of pseudocost evaluations after which a pseudocost is considered reliable
    unsigned BRANCHING_RELIABILITY_THRESHOLD;
    //! @brief Score factor (between 0 and 1) used to account for the left and right nodes in setting the pseudocost
    double BRANCHING_SCORE_FACTOR;
    //! @brief Relative tolerance within which a variable is considered to be at one of its bounds, i.e. excluded from branching variable selection
    double BRANCHING_BOUND_THRESHOLD;
    //! @brief Display option
    int DISPLAY;
    //! @brief Maximum CPU time limit
    double MAX_CPU_TIME;
    //! @brief Maximum number of nodes (0 for no limit)
    unsigned MAX_NODES;
    //! @brief Infinite
    double INF;
  } options;

  //! @brief SBB exceptions
  class Exceptions
  {
  public:
    //! @brief Enumeration type for SBB exception handling
    enum TYPE{
      BRANCH=0,		//!< Error due to an empty set of branching variables 
      INTERN=-3,	//!< SBB internal error
      UNDEF=-33	//!< Error due to calling a function/feature not yet implemented in SBB
    };
    //! @brief Constructor for error <a>ierr</a>
    Exceptions( TYPE ierr ) : _ierr( ierr ){}
    //! @brief Inline function returning the error flag
    int ierr()
      { return _ierr; }
    std::string what()
      {
        switch( _ierr ){
        case BRANCH:
          return "SBB Error: empty set of branching variables";
        case INTERN:
          return "SBB Internal Error";
        case UNDEF: default:
          return "SBB Error: calling a feature not yet implemented";
        }
      }
  private:
    TYPE _ierr;
  };

  //! @brief Set variables
  void variables
    ( const unsigned np, const T*P,
      const std::set<unsigned>&exclude=std::set<unsigned>() );

  //! @brief Retreive number of variables
  unsigned np() const
    { return _np; }

  //! @brief Apply branch-and-bound search
  std::pair<double,const double*> solve
    ( const TYPE pb, const double*p0=0, const double*f0=0,
      std::ostream&os=std::cout );

protected:  
  //! @brief Problem type 
  TYPE _pb;
  //! @brief Exclusion set for branching variables (e.g. linear variables)
  std::set<unsigned> _exclude_vars;
  //! @brief SBB status
  STATUS _status;
  //! @brief Variable values at current incumbent
  double *_p_inc;
  //! @brief Current incumbent value
  double _f_inc;

private:  
  //! @brief Status of bounding problems
  enum NODESTAT{
    FATHOM=0,		//!< Node should be fathomed
    NOTUPDATED,		//!< Incumbent for node has NOT been updated
    UPDATED,		//!< Incumbent for node HAS been updated
    ABORT		//!< Execution should be aborted
  };
  //! @brief Node type
  enum NODETYPE{
    ROOT=0,	//!< Root node
    LEFT=-1,	//!< Left child node
    RIGHT=1	//!< Right child node
  };

  //! @brief Number of variables in optimization problem
  unsigned _np;
  //! @brief Variable bounds at root node
  T *_P_root;
  //! @brief Default branching variable set
  std::set<unsigned> _branch_set;

  //! @brief Variable bounds (temporary)
  T *_P_tmp;
  //! @brief Current node index
  unsigned _node_index;
  //! @brief Index of node where current incumbent was found
  unsigned _node_inc;
  //! @brief Number of nodes introduced so far
  unsigned _node_cnt;
  //! @brief Max. number of nodes during the SBB search so far
  unsigned _node_max;
  //! @brief Set of nodes
  t_Nodes _Nodes;

  //! @brief Status after relaxed problem solution
  NODESTAT _REL_stat;
  //! @brief Status after original problem solution
  NODESTAT _ORI_stat;

  //! @brief Starting time
  double _tstart;
  //! @brief Current time
  double _tcur;
  //! @brief UDB cumulated time
  double _tUBD;
  //! @brief LDB cumulated time
  double _tLBD;

  //! @brief maximum number of values displayed in a row
  static const unsigned _LDISP = 4;
  //! @brief reserved space for integer variable display
  static const unsigned _IPREC = 6;
  //! @brief reserved space for double variable display
  static const unsigned _DPREC = 6;
  //! @brief stringstream for displaying results
  std::ostringstream _odisp;

  //! @brief Erase stored nodes
  void _clean_stack();
  //! @brief Reinitialize variables, incumbent, counters, time, etc.
  void _restart
    ( const double*p0, const double*f0 );
  //! @brief Default set of branching variables
  void _branching_variable_set();

  //! @brief Lower bound given node
  NODESTAT _lower_bound
    ( SBBNode<T>*pNode, const bool relaxed, const bool strongbranching=false );
  //! @brief Upper bound given node
  NODESTAT _upper_bound
    ( SBBNode<T>*pNode, const bool relaxed, const bool strongbranching=false );

  //! @brief Determines whether the relaxation solution point is feasible
  STATUS _feasible_relax
    ( SBBNode<T>*pNode );
  //! @brief Determines whether given node can be fathomed by value dominance
  bool _fathom_by_dominance
    ( SBBNode<T>*pNode );
  //! @brief Try and update incumbent
  bool _update_incumbent
    ( const double fINC, const double*pINC );

  //! @brief Apply preprocessing to given node
  NODESTAT _preprocess
    ( SBBNode<T>*pNode );
  //! @brief Apply postprocessing to given node
  NODESTAT _postprocess
    ( SBBNode<T>*pNode );
  //! @brief Branch current node and create subdomains
  NODESTAT _branch_node
    ( SBBNode<T>*pNode );
  //! @brief Selects the branching variable for given node
  std::pair<unsigned, double> _select_branching_variable
    ( SBBNode<T>*pNode, const std::set<unsigned>&set_branch );
  //! @brief Apply strong branching for branching variable selection
  std::pair<unsigned, double> _strong_branching
    ( SBBNode<T>*pNode, const std::set<unsigned>&set_branch );
  //! @brief Score branching variable
  double _score_branching
    ( const double fLEFT, const double fRIGHT ) const;
  //! @brief Partition branching variable domain for given node
  std::pair<const T, const T> _partition_variable_domain
    ( SBBNode<T>*pNode, const unsigned ip ) const;
  //! @brief Determine whether using pseudocost for branching is reliable
  bool _reliable_pseudocost
    ( SBBNode<T>*pNode, const std::set<unsigned>&set_branch );

  //! @brief Test entire SBB tree for fathoming by value dominance
  void _update_tree();
  //! @brief Apply postprocessing to entire SBB tree
  NODESTAT _postprocess_tree();

  //! @brief Add information for display
  void _display_add
    ( const double dval );
  void _display_add
    ( const unsigned ival );
  void _display_add
    ( const std::string &sval );

  //! @brief Initialize display
  void _display_init();
  //! @brief Final display
  void _display_final();

  //! @brief Display current buffer stream and reset it
  void _display
    ( std::ostream&os );
  //! @brief Display current nodes in tree
  void _display_tree( std::ostream&os );
  //! @brief Display current box (only 2 and 3 parameters)
  void _display_box
    ( const SBBNode<T>*pNode, std::ostream&os ) const;
};

//! @brief C++ base class for branch-and-bound nodes
////////////////////////////////////////////////////////////////////////
//! mc::SBBNode<T> is a C++ base class for defining nodes in the spatial
//! branch-and-bound algorithm.
////////////////////////////////////////////////////////////////////////
template <typename T>
class SBBNode
////////////////////////////////////////////////////////////////////////
{
  template <typename U> friend class SBB;
  template <typename U> friend struct lt_SBBNode;
  template <typename U> friend std::ostream& operator<<
    ( std::ostream&, const SBBNode<U>& );

public:  
  //! @brief Retreive node strength (based on parent bound)
  double strength() const
    { return _strength; }
  //! @brief Retreive node index
  unsigned index() const
    { return _index; }
  //! @brief Retreive node depth
  unsigned depth() const
    { return _depth; }
  //! @brief Retreive node iteration
  unsigned iter() const
    { return _iter; }
  //! @brief Retreive current lower bound point
  double pLB
    ( const unsigned ip ) const
    { assert( ip < _pSBB->_np ); return _LB_var[ip]; }
  //! @brief Retreive current upper bound point
  const double* pLB() const
    { return _LB_var; }
  //! @brief Retreive current lower bound value
  double& fLB()
    { return _LB_obj; }
  double fLB() const
    { return _LB_obj; }
  //! @brief Retreive current upper bound point
  double pUB
    ( const unsigned ip ) const
    { assert( ip < _pSBB->_np ); return _UB_var[ip]; }
  //! @brief Retreive current upper bound point
  const double* pUB() const
    { return _UB_var; }
  //! @brief Retreive current upper bound value
  double& fUB()
    { return _UB_obj; }
  double fUB() const
    { return _UB_obj; }

  //! @brief Retreive current variable bounds
  const T& P
    ( const unsigned ip ) const
    { assert( ip < _pSBB->_np ); return _P[ip]; }
  //! @brief Retreive/set bound for variable <a>ip</a> 
  T& P
    ( const unsigned ip )
    { assert( ip < _pSBB->_np ); return _P[ip]; }
  //! @brief Retreive pointer to variable bounds
  const T* P() const
    { return _P; }
  //! @brief Retreive/set pseudocost for variable <a>ip</a>
  std::pair<unsigned,double>& pc_var
    ( const unsigned ip )
    { assert( ip < 2*_pSBB->_np ); return _pc_var[ip]; }
  //! @brief Retreive variable pseudocosts
  const std::pair<unsigned,double>* pc_var() const
    { return _pc_var; }
  //! @brief Retreive/set node type (ROOT/LEFT/RIGHT)
  typename SBB<T>::NODETYPE& type()
    { return _type; }
  //! @brief Retreive/set parent branching variable and range
  std::pair<unsigned,T>& parent()
    { return _parent; }
  //! @brief Set strong branch status to true
  void strongbranch()
    { _strongbranch = true; }

  //! @brief Public constructor (for root node)
  SBBNode
    ( SBB<T>*pSBB, const T*P, const double*p0=0,
      const unsigned index=1, const unsigned iter=0 );

  //! @brief Public constructor (for regular node)
  SBBNode
    ( SBB<T>*pSBB, const T*P, const double strength,
      const unsigned index, const unsigned depth,
      const unsigned iter, const typename SBB<T>::NODETYPE type,
      const std::pair<unsigned,T> parent,
      const std::pair<unsigned,double>* pc_var );

  //! @brief Destructor
  ~SBBNode();

  //! @brief Lower bounding
  typename SBB<T>::STATUS lower_bound
    ( const double*var, const double inc );

  //! @brief Preprocessing
  typename SBB<T>::STATUS preprocess
    ( double*var, double&obj, const double inc );

  //! @brief Postprocessing
  typename SBB<T>::STATUS postprocess
    ( double*var, double&obj, const double inc );

  //! @brief Upper bounding
  typename SBB<T>::STATUS upper_bound
    ( const double*var, const double inc );

  //! @brief Feasibility test at <a>var</a>
  typename SBB<T>::STATUS test_feasibility
    ( double*var, double&obj, const double inc );

  //! @brief Check if the point <a>val</a> is at a bound
  bool at_bound
    ( const double val, const unsigned ip, const double rtol ) const;

  //! @brief Update the pseudocost of the current branching variable
  void update_pseudocost
    ( const typename SBB<T>::TYPE pb );

private:
  //! @brief Private default constructor
  SBBNode<T>(){};

  //! @brief Pointer to underlying branch-and-bound tree
  SBB<T> *_pSBB;
  //! @breif Strength of parent node
  double _strength;
  //! @brief Depth in SBB tree
  unsigned _depth;
  //! @brief Index in SBB tree
  unsigned _index;
  //! @brief Iteration when created in SBB tree
  unsigned _iter;
  //! @brief Node type (ROOT/LEFT/RIGHT)
  typename SBB<T>::NODETYPE _type;

  //! @brief Parent node branching variable and range
  std::pair<unsigned,T> _parent;
  //! @brief Whether or not strong branching has been applied
  bool _strongbranch;
  //! @brief Variable pseudocosts
  std::pair<unsigned,double> *_pc_var;

  //! @brief Variable bounds
  T *_P;
  //! @brief Backup variable bounds (for domain reduction)
  T *_P0;

  //! @brief Upper bound value
  double _UB_obj;
  //! @brief upper bound point
  double *_UB_var;

  //! @brief Lower bound value
  double _LB_obj;
  //! @brief Lower bound point
  double *_LB_var;
};

//! @brief C++ structure for comparing SBB Nodes
////////////////////////////////////////////////////////////////////////
//! mc::lt_SBBNode is a C++ structure for comparing nodes in branch-and-
//! bound tree based on their lower bound values.
////////////////////////////////////////////////////////////////////////
template <typename T>
struct lt_SBBNode
{
  bool operator()
    ( const SBBNode<T>*Node1, const SBBNode<T>*Node2 ) const
    { return( Node1->strength() < Node2->strength() ); }
};

///////////////////////////////   SBB   ////////////////////////////////

template <typename T>
inline void
SBB<T>::variables
( const unsigned np, const T*P, const std::set<unsigned>&exclude )
{
  if( !np || !P ){
    _np = 0;
    delete[] _P_root; _P_root = 0;
    delete[] _P_tmp;  _P_tmp = 0;
    delete[] _p_inc;  _p_inc = 0;
    return;
  }
  
  if( _np != np ){
    delete[] _P_root;
    delete[] _P_tmp;
    delete[] _p_inc;
    _np = np;
    _P_root = new T[_np];
    _P_tmp = new T[_np];
    _p_inc = 0;
  } 
  for( unsigned i=0; i<_np; i++ ) _P_root[i] = P[i];

  if( &exclude == &_exclude_vars ) return;
  _exclude_vars = exclude;

  return;
}

template <typename T>
inline std::pair< double, const double* >
SBB<T>::solve
( const TYPE pb, const double*p0, const double*f0, std::ostream&os )
{
  // Create and add root node to set _Nodes
  _pb = pb;
  _restart( p0, f0 );
  _display_init();
  _display( os );
  _Nodes.insert( new SBBNode<T>( this, _P_root, p0 ) );

#ifdef SAVE__SBB_NODES_TO_FILE
  std::ofstream osbbtree( "sbb_tree.out", std::ios_base::out );
#endif
  
  // keep branch-and-bound going until set _Nodes is empty
  for( _node_index = 1; !_Nodes.empty() && _tcur-_tstart < options.MAX_CPU_TIME
       && ( !options.MAX_NODES || _node_index <= options.MAX_NODES );
       _node_index++, _tcur=time() ){

    // intermediate display
    _display_add( _node_index );
    _display_add( (unsigned)_Nodes.size() );
    _display_add( time()-_tstart );
    switch( _pb ){
    case MINIMIZE:
      _display_add( (*_Nodes.begin())->strength() ); break;
    case MAXIMIZE:
      _display_add( -(*_Nodes.begin())->strength() ); break;
    }
    _display_add( _f_inc );

    // get pointer to next node and remove it for set _Nodes
    _display_tree( os );
    SBBNode<T>* pNode = *_Nodes.begin();
    _Nodes.erase( _Nodes.begin() );
    _display_add( pNode->iter() );

#ifdef SAVE__SBB_NODES_TO_FILE
    _display_box( pNode, osbbtree );
#endif

    // pre-processing
    typename SBB<T>::NODESTAT PRE_stat = _preprocess( pNode );
    if( PRE_stat == FATHOM ){
      delete pNode;
      _display_add( "SKIPPED" );
      _display_add( "SKIPPED" );
      _display_add( "FATHOM" );
      _display( os );
      continue;
    }
    else if( PRE_stat == ABORT ){
      delete pNode;
      _display_add( "SKIPPED" );
      _display_add( "SKIPPED" );
      _display_add( "ABORT" );
      _display( os );
      break;
    }

    // relaxation
    switch( _pb ){
    case MINIMIZE:
      _REL_stat = _lower_bound( pNode, true ); break;
    case MAXIMIZE:
      _REL_stat = _upper_bound( pNode, true ); break;
    }
    if( _REL_stat == FATHOM ){
      delete pNode;
      _display_add( "SKIPPED" );
      _display_add( "FATHOM" );
      _display( os );
      continue;
    }
    else if( _REL_stat == ABORT ){
      delete pNode;
      _display_add( "SKIPPED" );
      _display_add( "ABORT" );
      _display( os );
      break;
    }

    // tightening
    switch( _pb ){
    case MINIMIZE:
      _ORI_stat = _upper_bound( pNode, false ); break;
    case MAXIMIZE:
      _ORI_stat = _lower_bound( pNode, false ); break;
    }
    if( _ORI_stat == FATHOM ){
      delete pNode;
      _display_add( "FATHOM" );
      _display( os );
      continue;
    }
    else if( _ORI_stat == ABORT ){
      delete pNode;
      _display_add( "ABORT" );
      _display( os );
      break;
    }

    // SBB tree reduction
    if( _REL_stat==UPDATED || _ORI_stat==UPDATED ){
      _update_tree();
      if( _postprocess_tree() == ABORT ){
        delete pNode;
        _display_add( "ABORT" );
        _display( os );
        break;
      }
    }

    // post-processing
    typename SBB<T>::NODESTAT POST_stat = _postprocess( pNode );
    if( POST_stat == FATHOM ){
      delete pNode;
      _display_add( "FATHOM" );
      _display( os );
      continue;
    }
    else if( POST_stat == ABORT ){
      delete pNode;
      _display_add( "ABORT" );
      _display( os );
      break;
    }

    // domain branching
    _REL_stat = _branch_node( pNode );
    delete pNode;
    if( _REL_stat == ABORT ) _display_add( "ABORT" );
    _display( os );

    // keep track of the maximum nodes in memory
    if( _Nodes.size() > _node_max ) _node_max = _Nodes.size();
  }

#ifdef SAVE__SBB_NODES_TO_FILE
  osbbtree.close();
#endif

  _display_final();
  _display( os );
  _status = (!_Nodes.empty()? FAILURE: (_p_inc? NORMAL: INFEASIBLE));
  return std::make_pair( _f_inc, _p_inc );
}

template <typename T>
inline void
SBB<T>::_restart
( const double*p0, const double*f0 )
{
  _clean_stack();
  _branching_variable_set();
  
  delete[] _p_inc; _p_inc = 0;
  if( p0 && f0 ){
    _p_inc = new double[_np];
    for( unsigned ip=0; ip<_np; ip++ )
      _p_inc[ip] = p0[ip]; 
    _f_inc = *f0;
  }
  else
    _f_inc = ( _pb==MINIMIZE? options.INF: -options.INF );

  _node_index = _node_inc = _node_cnt = _node_max = 0;
  _tstart = _tcur = time();
  _tUBD = _tLBD = 0;
}

template <typename T>
inline void
SBB<T>::_clean_stack()
{
  it_Nodes it = _Nodes.begin();
  for( ; it != _Nodes.end(); it++ ) delete *it;
  _Nodes.clear();
}

template <typename T>
inline void
SBB<T>::_branching_variable_set()
{
  _branch_set.clear();
  for( unsigned ip=0; ip<_np; ip++ ){
    if( _exclude_vars.empty() || _exclude_vars.find(ip)==_exclude_vars.end() )
      _branch_set.insert( ip );
  }
  if( _branch_set.empty() ) throw Exceptions( Exceptions::BRANCH );
}

template <typename T>
inline typename SBB<T>::NODESTAT
SBB<T>::_lower_bound
( SBBNode<T>*pNode, const bool relaxed, const bool strongbranching )
{
  _tLBD -= time();
  NODESTAT stat;
  switch( pNode->lower_bound( pNode->pUB(), _f_inc ) ){

    case INFEASIBLE:
      if( !strongbranching ) _display_add( "INFEASIBLE" );
      if( relaxed ){
        pNode->fLB() = options.INF;
        pNode->update_pseudocost( _pb );
      }
      stat = FATHOM; break;
      
    case NORMAL:{
      if( !strongbranching ) _display_add( pNode->fLB() );
      if( relaxed ){
        // update variable pseudo-cost
        pNode->update_pseudocost( _pb );
        // fathom by value dominance?
        if( _fathom_by_dominance( pNode ) )
          { stat = FATHOM; break; }
        // try and update incumbent w/ relaxation solution point if applicable
        if( _feasible_relax( pNode )==NORMAL
         && _update_incumbent( pNode->fUB(), pNode->pLB() ) )
          { stat = UPDATED; break; }
      }
      else{
        // try and update incumbent
        if( _update_incumbent( pNode->fLB(), pNode->pLB() ) )
          { stat = UPDATED; break; }
      }
      stat = NOTUPDATED; break;
    }

    case FAILURE:
      if( !strongbranching ) _display_add( "FAILURE" );
      // retain node (to be sure...)
      pNode->fLB() = -options.INF;
      stat = NOTUPDATED; break;

    case FATAL: default:
      if( !strongbranching ) _display_add( "FATAL" );
      stat = ABORT; break;
  }
  _tLBD += time();
  return stat;
}

template <typename T>
inline typename SBB<T>::NODESTAT
SBB<T>::_upper_bound
( SBBNode<T>*pNode, const bool relaxed, const bool strongbranching )
{
  _tUBD -= time();
  NODESTAT stat;
  switch( pNode->upper_bound( pNode->pLB(), _f_inc ) ){

    case INFEASIBLE:
      if( !strongbranching ) _display_add( "INFEASIBLE" );
      if( relaxed ){
        pNode->fUB() = -options.INF;
        pNode->update_pseudocost( _pb );
      }
      stat = FATHOM; break;
      
    case NORMAL:{
      if( !strongbranching ) _display_add( pNode->fUB() );
      if( relaxed ){
        // update variable pseudo-cost
        pNode->update_pseudocost( _pb );
        // fathom by value dominance?
        if( _fathom_by_dominance( pNode ) )
          { stat = FATHOM; break; }
        // try and update incumbent w/ relaxation solution point if applicable
        if( _feasible_relax( pNode )==NORMAL
         && _update_incumbent( pNode->fLB(), pNode->pUB() ) )
          { stat = UPDATED; break; }
      }
      else{
        // try and update incumbent
        if( _update_incumbent( pNode->fUB(), pNode->pUB() ) ) stat = UPDATED;
      }
      stat = NOTUPDATED; break;
    }

    case FAILURE:
      if( !strongbranching ) _display_add( "FAILURE" );
      // retain node (to be sure...)
      pNode->fUB() = options.INF;
      stat = NOTUPDATED; break;

    case FATAL: default:
      if( !strongbranching ) _display_add( "FATAL" );
      stat = ABORT; break;
  }
  _tUBD += time();
  return stat;
}

template <typename T>
inline typename SBB<T>::STATUS
SBB<T>::_feasible_relax
( SBBNode<T>*pNode )
{
  switch( _pb ){
  case MINIMIZE:
    return pNode->test_feasibility( pNode->_LB_var, pNode->_UB_obj, _f_inc );
  case MAXIMIZE: default:
    return pNode->test_feasibility( pNode->_UB_var, pNode->_LB_obj, _f_inc );
  }
}

template <typename T>
inline bool
SBB<T>::_fathom_by_dominance
( SBBNode<T>*pNode )
{
  double f_eff = _f_inc;
  switch( _pb ){

  case MINIMIZE:
    f_eff -= std::max( options.ABSOLUTE_TOLERANCE,
      std::fabs(_f_inc)*options.RELATIVE_TOLERANCE );
    // try and fathom current node by value dominance
    return( pNode->fLB() > f_eff ? true: false );
    
  case MAXIMIZE: default:
    f_eff += std::max( options.ABSOLUTE_TOLERANCE,
      std::fabs(_f_inc)*options.RELATIVE_TOLERANCE );
    // try and fathom current node by value dominance
    return( pNode->fUB() < f_eff ? true: false );
  }
}

template <typename T>
inline bool
SBB<T>::_update_incumbent
( const double fINC, const double*pINC )
{
  // test incumbent w.r.t current solution point
  switch( _pb ){
  case MINIMIZE: if( fINC >= _f_inc ) return false; break;
  case MAXIMIZE: if( fINC <= _f_inc ) return false; break;
  }

  _f_inc = fINC;
  if( !_p_inc ) _p_inc = new double[_np];
  for( unsigned ip=0; ip<_np; ip++ ) _p_inc[ip] = pINC[ip];//Node->pLB(ip);
  _node_inc = _node_index;
  return true;
}

template <typename T>
inline std::pair<unsigned, double>
SBB<T>::_select_branching_variable
( SBBNode<T>*pNode, const std::set<unsigned>&set_branch )
{
  std::pair<unsigned, double> branchsel( _np, -options.INF );
  
  switch( options.BRANCHING_VARIABLE_CRITERION ){
  // Branching based on pseudocosts
  case Options::RELIAB:
    if( _reliable_pseudocost( pNode, set_branch ) ){
      std::set<unsigned>::iterator it = set_branch.begin();
#ifdef DEBUG__SBB_BRANCHINGVAR
      std::cout << "*** SCORES";
#endif
      for( ; it != set_branch.end(); ++it ){
        std::pair<const T, const T> partition
          = _partition_variable_domain( pNode, *it );
#ifdef USE__SBB_PSEUDOCOST_FILTER
        double score = _score_branching( Op<T>::diam( partition.second )
          * pNode->pc_var()[*it].second, Op<T>::diam( partition.first )
          * pNode->pc_var()[*it+_np].second );
#else
        double score = _score_branching( Op<T>::diam( partition.second )
          * pNode->pc_var()[*it].second / pNode->pc_var()[*it].first,
                                         Op<T>::diam( partition.first )
          * pNode->pc_var()[*it+_np].second / pNode->pc_var()[*it+_np].first );
#endif
#ifdef DEBUG__SBB_BRANCHINGVAR
        std::cout << "   " << *it << ": " << score;
#endif
        if( score > branchsel.second ) branchsel = std::make_pair( *it, score );
      }
#ifdef DEBUG__SBB_BRANCHINGVAR
      std::cout << std::endl;
#endif
      break;
    }
    // Unreliable pseudocosts - Use strong branching
    branchsel = _strong_branching( pNode, set_branch );
    if( _REL_stat == ABORT || branchsel.first < _np ) break;
    // No break - use RGREL criterion in case strong branching failed

  // Branching based on relative range diameter
  case Options::RGREL:{
    std::set<unsigned>::iterator it = set_branch.begin();
    for( ; it != set_branch.end(); ++it ){
      double score = Op<T>::diam( pNode->P(*it) )
                   / Op<T>::diam( _P_root[*it] );
      if( score > branchsel.second )
        branchsel = std::make_pair( *it, score );
    }
    break;
   }
  // Branching based on absolute range diameter
  case Options::RGABS:{
    std::set<unsigned>::iterator it = set_branch.begin();
    for( ; it != set_branch.end(); ++it ){
      double score = Op<T>::diam( pNode->P(*it) );
      if( score > branchsel.second )
        branchsel = std::make_pair( *it, score );
    }
    break;
   }
  }

  return branchsel;
}

template <typename T>
inline bool
SBB<T>::_reliable_pseudocost
( SBBNode<T>*pNode, const std::set<unsigned>&set_branch )
{
  std::set<unsigned>::iterator it = set_branch.begin();
  for( ; it != set_branch.end(); ++it ){
    if( !pNode->pc_var()[*it].first || pNode->pc_var()[*it].first
      < options.BRANCHING_RELIABILITY_THRESHOLD ) return false;
    if( !pNode->pc_var()[_np+*it].first || pNode->pc_var()[_np+*it].first
      < options.BRANCHING_RELIABILITY_THRESHOLD ) return false;
  }
  return true;
}

template <typename T>
inline std::pair<const T, const T>
SBB<T>::_partition_variable_domain
( SBBNode<T>*pNode, const unsigned ip_branch )
const
{
  std::pair<T,T> partition;

  switch( options.BRANCHING_STRATEGY ){
  case Options::OMEGA:
    // Branch at incumbent if interior to current variable range
    if( _p_inc && !pNode->at_bound(_p_inc[ip_branch], ip_branch,
      options.BRANCHING_BOUND_THRESHOLD) ){
      partition.first = mc::Op<T>::l(pNode->P(ip_branch)) + Op<T>::zeroone()
        *( _p_inc[ip_branch] - mc::Op<T>::l(pNode->P(ip_branch)) );
      partition.second = mc::Op<T>::u(pNode->P(ip_branch)) + Op<T>::zeroone()
        *( _p_inc[ip_branch] - mc::Op<T>::u(pNode->P(ip_branch)) );
      break;
    }
    // Otherwise branch at relaxation solution point if available +
    // interior to current variable range
    if( _pb == MINIMIZE && pNode->fLB() > -options.INF
     && !pNode->at_bound(pNode->pLB(ip_branch), ip_branch,
     options.BRANCHING_BOUND_THRESHOLD) ){ 
      partition.first = mc::Op<T>::l(pNode->P(ip_branch)) + Op<T>::zeroone()
        *( pNode->pLB(ip_branch) - mc::Op<T>::l(pNode->P(ip_branch)) );
      partition.second = mc::Op<T>::u(pNode->P(ip_branch)) + Op<T>::zeroone()
        *( pNode->pLB(ip_branch) - mc::Op<T>::u(pNode->P(ip_branch)) );
      break;
    }
    else if( _pb == MAXIMIZE && pNode->fUB() < options.INF
     && !pNode->at_bound(pNode->pUB(ip_branch), ip_branch,
     options.BRANCHING_BOUND_THRESHOLD) ){ 
      partition.first = mc::Op<T>::l(pNode->P(ip_branch)) + Op<T>::zeroone()
        *( pNode->pUB(ip_branch) - mc::Op<T>::l(pNode->P(ip_branch)) );
      partition.second = mc::Op<T>::u(pNode->P(ip_branch)) + Op<T>::zeroone()
        *( pNode->pUB(ip_branch) - mc::Op<T>::u(pNode->P(ip_branch)) );
      break;
    }
    
  case Options::MIDPOINT: default:
    // Branch at mid-point of current variable range
    partition.first = mc::Op<T>::l(pNode->P(ip_branch)) + Op<T>::zeroone()
      *Op<T>::diam(pNode->P(ip_branch))/2.;
    partition.second = mc::Op<T>::u(pNode->P(ip_branch)) - Op<T>::zeroone()
      *Op<T>::diam(pNode->P(ip_branch))/2.;
    break;
  }

  return partition;
}

template <typename T>
inline typename SBB<T>::NODESTAT
SBB<T>::_branch_node
( SBBNode<T>*pNode )
{
  // Branching variable selection
  const unsigned ip_branch = _select_branching_variable( pNode,
    _branch_set ).first;
  if( _REL_stat == ABORT ) return ABORT;

  // Partitionning
  std::pair<const T, const T> partition = _partition_variable_domain(
    pNode, ip_branch );
  for( unsigned ip=0; ip<_np; ip++ ) _P_tmp[ip] = pNode->P(ip);
  _P_tmp[ip_branch] = partition.first;
  switch( _pb ){
  case MINIMIZE:
    _Nodes.insert( new SBBNode<T>( this, _P_tmp, pNode->fLB(), ++_node_cnt,
      pNode->depth()+1, _node_index, LEFT, std::make_pair(ip_branch,
      pNode->P(ip_branch)), pNode->pc_var() ) ); break;
  case MAXIMIZE:
    _Nodes.insert( new SBBNode<T>( this, _P_tmp, -pNode->fUB(), ++_node_cnt,
      pNode->depth()+1, _node_index, LEFT, std::make_pair(ip_branch,
      pNode->P(ip_branch)), pNode->pc_var() ) ); break;
  }
  _P_tmp[ip_branch] = partition.second;
  switch( _pb ){
  case MINIMIZE:
    _Nodes.insert( new SBBNode<T>( this, _P_tmp, pNode->fLB(), ++_node_cnt,
      pNode->depth()+1, _node_index, RIGHT, std::make_pair(ip_branch,
      pNode->P(ip_branch)), pNode->pc_var() ) ); break;
  case MAXIMIZE:
    _Nodes.insert( new SBBNode<T>( this, _P_tmp, -pNode->fUB(), ++_node_cnt,
      pNode->depth()+1, _node_index, RIGHT, std::make_pair(ip_branch,
      pNode->P(ip_branch)), pNode->pc_var() ) ); break;
  }

  // Display
  std::ostringstream omsg;
  omsg << "BRANCH" << ip_branch;
  _display_add( omsg.str() );

  return UPDATED;
}

template <typename T>
inline double
SBB<T>::_score_branching
( const double fLEFT, const double fRIGHT ) const
{
  return (1-options.BRANCHING_SCORE_FACTOR) * fmin(fLEFT,fRIGHT)
          + options.BRANCHING_SCORE_FACTOR  * fmax(fLEFT,fRIGHT);
}

template <typename T>
inline std::pair<unsigned, double>
SBB<T>::_strong_branching
( SBBNode<T>*pNode, const std::set<unsigned>&set_branch )
{
  std::pair<unsigned, double> branchsel( _np, -options.INF );
  if( !options.BRANCHING_RELIABILITY_THRESHOLD ) return branchsel;

  // Create child subnode
  //for( unsigned ip=0; ip<_np; ip++ ) _P_tmp[ip] = pNode->P(ip);
  SBBNode<T>* pSubnode = 0;
  switch( _pb ){
  case MINIMIZE:
    if( pNode->fLB() == -options.INF ) return branchsel;
    pSubnode = new SBBNode<T>( this, pNode->P(), pNode->fLB(), pNode->index(),
      pNode->depth()+1, _node_index, LEFT, std::make_pair(0,pNode->P(0)), pNode->pc_var() ); break;
  case MAXIMIZE:
    if( pNode->fUB() == options.INF ) return branchsel;
    pSubnode = new SBBNode<T>( this, pNode->P(), -pNode->fUB(), pNode->index(),
      pNode->depth()+1, _node_index, LEFT, std::make_pair(0,pNode->P(0)), pNode->pc_var() ); break;
  }
  
  // Repeat for all variables in branching set
#ifdef DEBUG__SBB_STRONGBRANCHING
  std::cout << "*** SCORES";
#endif
  double fLEFT, fRIGHT, fSCORE;
  std::set<unsigned>::iterator it = set_branch.begin();
  for( ; it != set_branch.end(); ++it ){

    // Simulate partitionning
    pSubnode->parent() = std::make_pair( *it, pNode->P(*it) );
    const std::pair<const T, const T> partition
      = _partition_variable_domain( pNode, *it );

    // Relax left partition
    pSubnode->P(*it) = partition.first;
    pSubnode->type() = LEFT;
    switch( _pb ){
    case MINIMIZE:
      _REL_stat = _lower_bound( pSubnode, true, true );
      fLEFT = pSubnode->fLB(); break;
    case MAXIMIZE: default:
      _REL_stat = _upper_bound( pSubnode, true, true );
      fLEFT = -pSubnode->fUB(); break;
    }
    if( _REL_stat == ABORT ) break;
    // Left partition failed
    //if( fLEFT == -options.INF ) continue;
 
    // Relax right partition
    pSubnode->P(*it) = partition.second;
    pSubnode->type() = RIGHT;
    switch( _pb ){
    case MINIMIZE:
      _REL_stat = _lower_bound( pSubnode, true, true );
      fRIGHT = pSubnode->fLB(); break;
    case MAXIMIZE: default:
      _REL_stat = _upper_bound( pSubnode, true, true );
      fRIGHT = -pSubnode->fUB(); break;
    }
    if( _REL_stat == ABORT ) break;
    // Right partition failed
    //if( fRIGHT == -options.INF ) continue;

    // Score current variable and compare
    fSCORE = _score_branching( fLEFT, fRIGHT );
#ifdef DEBUG__SBB_STRONGBRANCHING
    std::cout << "   " << *it << ": " << fSCORE
              << " (" << fLEFT << "," << fRIGHT << ")";
#endif
    if( fSCORE > branchsel.second )
      branchsel = std::make_pair( *it, fSCORE );

    // Reset variable bounds
    for( unsigned ip=0; ip<_np; ip++ )
      pSubnode->P(ip) = pNode->P(ip);
  }
#ifdef DEBUG__SBB_STRONGBRANCHING
  std::cout << std::endl;
#endif

  // Update variable pseudocosts
  for( unsigned ip=0; ip<2*_np; ip++ )
    pNode->pc_var(ip) = pSubnode->pc_var()[ip];
  delete pSubnode;
  pNode->strongbranch();

  return branchsel;
}

template <typename T>
inline typename SBB<T>::NODESTAT
SBB<T>::_preprocess
( SBBNode<T>*pNode )
{
  switch( _pb==MINIMIZE?
          pNode->preprocess( pNode->_LB_var, pNode->_LB_obj, _f_inc ):
          pNode->preprocess( pNode->_UB_var, pNode->_UB_obj, _f_inc ) ){
  case INFEASIBLE:
    return FATHOM;
  case NORMAL: case FAILURE: default:
    return NOTUPDATED;
  }
}

template <typename T>
inline typename SBB<T>::NODESTAT
SBB<T>::_postprocess
( SBBNode<T>*pNode )
{
  if( _fathom_by_dominance( pNode ) ) return FATHOM;
  
  switch( _pb==MINIMIZE?
          pNode->postprocess( pNode->_LB_var, pNode->_LB_obj, _f_inc ):
          pNode->postprocess( pNode->_UB_var, pNode->_UB_obj, _f_inc ) ){
  case INFEASIBLE:
    return FATHOM;
  case NORMAL: case FAILURE:
    return NOTUPDATED;
  case FATAL: default:
    return ABORT;
  }
}

template <typename T>
inline void
SBB<T>::_update_tree()
{
  // try and reduce all nodes in B&B tree
  it_Nodes it = _Nodes.begin();
  for( ; it != _Nodes.end(); ++it )
    if( _fathom_by_dominance( *it ) ){ delete *it; _Nodes.erase( it ); }
  return;
}

template <typename T>
inline typename SBB<T>::NODESTAT
SBB<T>::_postprocess_tree()
{
  if( !options.TREE_REREDUCE ) return NOTUPDATED;

  // try and reduce all nodes in B&B tree
  it_Nodes it = _Nodes.begin();
  for( ; it != _Nodes.end(); ++it ){
    typename SBB<T>::NODESTAT POST_stat = _postprocess( *it );
    if( POST_stat == ABORT ) return ABORT;
    else if( POST_stat == FATHOM ){ delete *it; _Nodes.erase( it ); }
  }
  return UPDATED;
}

template <typename T>
inline void
SBB<T>::_display_tree
( std::ostream&os )
{
  if( options.DISPLAY <= 2 ) return;

  // Show current node
  if( options.DISPLAY == 3 ){
    os << *(*_Nodes.begin());
    return;
  }

  // Show all nodes in stack
  it_Nodes it = _Nodes.begin();
  for( unsigned id=1; it != _Nodes.end(); ++it, id++ ){
    os << "Node " << id << ":" << std::scientific << std::setprecision(_DPREC)
       << std::setw(_DPREC+8) << (*it)->strength() << *(*it);
  }
  return;
}

template <typename T>
inline void
SBB<T>::_display_init()
{
  _odisp.str("");
  if( options.DISPLAY <= 1 ) return;
  _odisp << std::right
  	 << std::setw(_IPREC) << "INDEX"
  	 << std::setw(_IPREC) << "STACK"
  	 << std::setw(_DPREC+8) << "CUMUL TIME"
  	 << std::setw(_DPREC+8) << "RELAX "
  	 << std::setw(_DPREC+8) << "INC   "
  	 << std::setw(_IPREC) << "PARENT"
  	 << std::setw(_DPREC+8) << ( _pb==MINIMIZE? "LBD   ": "UBD   " )
  	 << std::setw(_DPREC+8) << ( _pb==MINIMIZE? "UBD   ": "LBD   " )
  	 << std::setw(_DPREC+8) << "ACTION  "
  	 << std::endl;
  
}

template <typename T>
inline void
SBB<T>::_display_add
( const double dval )
{
  if( options.DISPLAY <= 1 ) return;
  _odisp << std::right << std::scientific << std::setprecision(_DPREC)
         << std::setw(_DPREC+8) << dval;
}

template <typename T>
inline void
SBB<T>::_display_add
( const unsigned ival )
{
  if( options.DISPLAY <= 1 ) return;
  _odisp << std::right << std::setw(_IPREC) << ival;
}

template <typename T>
inline void
SBB<T>::_display_add
( const std::string &sval )
{
  if( options.DISPLAY <= 1 ) return;
  _odisp << std::right << std::setw(_DPREC+8) << sval;
}

template <typename T>
inline void
SBB<T>::_display_final()
{
  if( options.DISPLAY <= 0 ) return;
  // Solution found within allowed time?
  if( _tcur-_tstart < options.MAX_CPU_TIME
    && ( !options.MAX_NODES || _node_index <= options.MAX_NODES ) )
    _odisp << std::endl << "#  NORMAL TERMINATION: ";
  else
    _odisp << std::endl << "#  EXECUTION STOPPED: ";
  _odisp << std::fixed << std::setprecision(6) << _tcur-_tstart << " CPU SEC"
         << "  (LBD:" << std::fixed << std::setprecision(1)
         << _tLBD/(_tcur-_tstart)*1e2 << "%  UBD:" << std::fixed
         << std::setprecision(1) << _tUBD/(_tcur-_tstart)*1e2 << "%)"
         << std::endl;

  // No feasible solution found
  if( !_p_inc ){
    _odisp << "#  NO FEASIBLE SOLUTION FOUND" << std::endl;
  }

  // Feasible solution found
  else{
    // Incumbent
    _odisp << "#  INCUMBENT VALUE:" << std::scientific
           << std::setprecision(_DPREC) << std::setw(_DPREC+8) << _f_inc
           << std::endl;
    _odisp << "#  INCUMBENT POINT:";
    for( unsigned ip=0, id=0; ip<_np; ip++, id++ ){
      if( id == _LDISP ){
        _odisp << std::endl << std::left << std::setw(19) << "#";
        id = 0;
      }
      _odisp << std::right << std::setw(_DPREC+8) << _p_inc[ip];
    }
    _odisp << std::endl
           << "#  INCUMBENT FOUND AT NODE: " << _node_inc << std::endl;
  }

  _odisp << "#  TOTAL NUMBER OF NODES:   " << _node_index-1 << std::endl
         << "#  MAXIMUM NODES IN STACK:  " << _node_max << std::endl;
}

template <typename T>
inline void
SBB<T>::_display
( std::ostream &os )
{
  if( _odisp.str() == "" ) return;
  //if( options.DISPLAY > 0 ){
  os << _odisp.str() << std::endl;
  //}
  _odisp.str("");
  return;
}

template <typename T>
inline void
SBB<T>::Options::display
( std::ostream&out ) const
{
  // Display SBB Options
  out << std::setw(60) << "  ABSOLUTE CONVERGENCE TOLERANCE"
      << std::scientific << std::setprecision(1)
      << ABSOLUTE_TOLERANCE << std::endl;
  out << std::setw(60) << "  RELATIVE CONVERGENCE TOLERANCE"
      << std::scientific << std::setprecision(1)
      << RELATIVE_TOLERANCE << std::endl;
  out << std::setw(60) << "  MAXIMUM NUMBER OF ITERATIONS";
  switch( MAX_NODES ){
  case 0:  out << "NO LIMIT\n"; break;
  default: out << MAX_NODES << std::endl; break;
  }
  out << std::setw(60) << "  MAXIMUM CPU TIME (SEC)"
      << std::scientific << std::setprecision(1)
      << MAX_CPU_TIME << std::endl;
  out << std::setw(60) << "  DISPLAY LEVEL"
      << DISPLAY << std::endl;
  out << std::setw(60) << "  INTERNAL VALUE FOR INFINITY"
      << std::scientific << std::setprecision(1)
      << INF << std::endl;
  out << std::setw(60) << "  UPDATE ENTIRE TREE AFTER AN INCUMBENT IMPROVEMENT?"
      << (TREE_REREDUCE?"Y\n":"N\n");
  out << std::setw(60) << "  BRANCHING STRATEGY FOR NODE PARTITIONING";
  switch( BRANCHING_STRATEGY ){
  case MIDPOINT: out << "MIDPOINT\n"; break;
  case OMEGA:    out << "OMEGA\n";    break;
  }
  out << std::setw(60) << "  BRANCHING TOLERANCE FOR VARIABLE AT BOUND"
      << std::scientific << std::setprecision(2)
      << BRANCHING_BOUND_THRESHOLD << std::endl;
  out << std::setw(60) << "  BRANCHING STRATEGY FOR VARIABLE SELECTION";
  switch( BRANCHING_VARIABLE_CRITERION ){
  case RGREL:  out << "RGREL\n";  break;
  case RGABS:  out << "RGABS\n";  break;
  case RELIAB: out << "RELIAB\n"; break;
  }
  out << std::setw(60) << "  INITIALIZATION THRESHOLD IN RELIABILITY BRANCHING"
      << BRANCHING_RELIABILITY_THRESHOLD << std::endl;
  out << std::setw(60) << "  SCORE PARAMETER IN RELIABILITY BRANCHING"
      << std::fixed << std::setprecision(2)
      << BRANCHING_SCORE_FACTOR << std::endl;
}

template <typename T>
inline void
SBB<T>::_display_box
( const SBBNode<T>*pNode, std::ostream&os ) const
{
  switch( _np ){
  case 2:
    os << Op<T>::l(pNode->P(0)) << "  " << Op<T>::l(pNode->P(1)) << std::endl
       << Op<T>::l(pNode->P(0)) << "  " << Op<T>::u(pNode->P(1)) << std::endl
       << std::endl;
    os << Op<T>::l(pNode->P(0)) << "  " << Op<T>::l(pNode->P(1)) << std::endl
       << Op<T>::u(pNode->P(0)) << "  " << Op<T>::l(pNode->P(1)) << std::endl
       << std::endl;
    os << Op<T>::u(pNode->P(0)) << "  " << Op<T>::u(pNode->P(1)) << std::endl
       << Op<T>::l(pNode->P(0)) << "  " << Op<T>::u(pNode->P(1)) << std::endl
       << std::endl;
    os << Op<T>::u(pNode->P(0)) << "  " << Op<T>::u(pNode->P(1)) << std::endl
       << Op<T>::u(pNode->P(0)) << "  " << Op<T>::l(pNode->P(1)) << std::endl
       << std::endl;
    break;

  case 3:
    os << Op<T>::l(pNode->P(0)) << "  " << Op<T>::l(pNode->P(1))
       << "  " << Op<T>::l(pNode->P(2)) << std::endl
       << Op<T>::l(pNode->P(0)) << "  " << Op<T>::u(pNode->P(1))
       << "  " << Op<T>::l(pNode->P(2)) << std::endl
       << std::endl << std::endl;
    os << Op<T>::l(pNode->P(0)) << "  " << Op<T>::l(pNode->P(1))
       << "  " << Op<T>::l(pNode->P(2)) << std::endl
       << Op<T>::u(pNode->P(0)) << "  " << Op<T>::l(pNode->P(1))
       << "  " << Op<T>::l(pNode->P(2)) << std::endl
       << std::endl << std::endl;
    os << Op<T>::u(pNode->P(0)) << "  " << Op<T>::u(pNode->P(1))
       << "  " << Op<T>::l(pNode->P(2)) << std::endl
       << Op<T>::l(pNode->P(0)) << "  " << Op<T>::u(pNode->P(1))
       << "  " << Op<T>::l(pNode->P(2)) << std::endl
       << std::endl << std::endl;
    os << Op<T>::u(pNode->P(0)) << "  " << Op<T>::u(pNode->P(1))
       << "  " << Op<T>::l(pNode->P(2)) << std::endl
       << Op<T>::u(pNode->P(0)) << "  " << Op<T>::l(pNode->P(1))
       << "  " << Op<T>::l(pNode->P(2)) << std::endl
       << std::endl << std::endl;

    os << Op<T>::l(pNode->P(0)) << "  " << Op<T>::l(pNode->P(1))
       << "  " << Op<T>::u(pNode->P(2)) << std::endl
       << Op<T>::l(pNode->P(0)) << "  " << Op<T>::u(pNode->P(1))
       << "  " << Op<T>::u(pNode->P(2)) << std::endl
       << std::endl << std::endl;
    os << Op<T>::l(pNode->P(0)) << "  " << Op<T>::l(pNode->P(1))
       << "  " << Op<T>::u(pNode->P(2)) << std::endl
       << Op<T>::u(pNode->P(0)) << "  " << Op<T>::l(pNode->P(1))
       << "  " << Op<T>::u(pNode->P(2)) << std::endl
       << std::endl << std::endl;
    os << Op<T>::u(pNode->P(0)) << "  " << Op<T>::u(pNode->P(1))
       << "  " << Op<T>::u(pNode->P(2)) << std::endl
       << Op<T>::l(pNode->P(0)) << "  " << Op<T>::u(pNode->P(1))
       << "  " << Op<T>::u(pNode->P(2)) << std::endl
       << std::endl << std::endl;
    os << Op<T>::u(pNode->P(0)) << "  " << Op<T>::u(pNode->P(1))
       << "  " << Op<T>::u(pNode->P(2)) << std::endl
       << Op<T>::u(pNode->P(0)) << "  " << Op<T>::l(pNode->P(1))
       << "  " << Op<T>::u(pNode->P(2)) << std::endl
       << std::endl << std::endl;

    os << Op<T>::l(pNode->P(0)) << "  " << Op<T>::l(pNode->P(1))
       << "  " << Op<T>::l(pNode->P(2)) << std::endl
       << Op<T>::l(pNode->P(0)) << "  " << Op<T>::l(pNode->P(1))
       << "  " << Op<T>::u(pNode->P(2)) << std::endl
       << std::endl << std::endl;
    os << Op<T>::u(pNode->P(0)) << "  " << Op<T>::l(pNode->P(1))
       << "  " << Op<T>::l(pNode->P(2)) << std::endl
       << Op<T>::u(pNode->P(0)) << "  " << Op<T>::l(pNode->P(1))
       << "  " << Op<T>::u(pNode->P(2)) << std::endl
       << std::endl << std::endl;
    os << Op<T>::l(pNode->P(0)) << "  " << Op<T>::u(pNode->P(1))
       << "  " << Op<T>::l(pNode->P(2)) << std::endl
       << Op<T>::l(pNode->P(0)) << "  " << Op<T>::u(pNode->P(1))
       << "  " << Op<T>::u(pNode->P(2)) << std::endl
       << std::endl << std::endl;
    os << Op<T>::u(pNode->P(0)) << "  " << Op<T>::u(pNode->P(1))
       << "  " << Op<T>::l(pNode->P(2)) << std::endl
       << Op<T>::u(pNode->P(0)) << "  " << Op<T>::u(pNode->P(1))
       << "  " << Op<T>::u(pNode->P(2)) << std::endl
       << std::endl << std::endl;
    break;

  default:
    break;
  }
}

/////////////////////////////// SBBNode ///////////////////////////////

template <typename T>
inline
SBBNode<T>::SBBNode
( SBB<T>*pSBB, const T*P, const double*p0,
  const unsigned index, const unsigned iter ):
_pSBB( pSBB ), _strength( -_pSBB->options.INF ), _depth( 0 ),
_index( index ), _iter( iter ), _type( SBB<T>::ROOT ),
_parent( 0, T(0.) ), _strongbranch( false )
{
  _P0 = new T[_pSBB->_np];
  _P  = new T[_pSBB->_np];
  for( unsigned i=0; i<_pSBB->_np; i++ ) _P[i] = P[i];
  
  _UB_obj = _pSBB->options.INF;
  _UB_var = new double[_pSBB->_np];
  _LB_obj = -_pSBB->options.INF;
  _LB_var = new double[_pSBB->_np];
  _pc_var = new std::pair<unsigned,double>[2*_pSBB->_np];
  
  // initial points and pseudocosts
  for( unsigned i=0; i<_pSBB->_np; i++ ){
    _LB_var[i] = _UB_var[i] = ( p0? p0[i]: mc::Op<T>::mid( P[i] ) );
    _pc_var[i] = _pc_var[_pSBB->_np+i] = std::make_pair(0,0.);
  }
}

template <typename T>
inline
SBBNode<T>::SBBNode
( SBB<T>*pSBB, const T*P, const double strength,
  const unsigned index, const unsigned depth,
  const unsigned iter, const typename SBB<T>::NODETYPE type,
  const std::pair<unsigned,T> parent,
  const std::pair<unsigned,double>* pc_var ):
_pSBB( pSBB ), _strength( strength ), _depth( depth ), _index( index ),
_iter( iter ), _type( type ), _parent( parent ), _strongbranch( false )
{
  _P0 = new T[_pSBB->_np];
  _P  = new T[_pSBB->_np];
  for( unsigned i=0; i<_pSBB->_np; i++ ) _P[i] = P[i];
  
  _UB_obj = _pSBB->options.INF;
  _UB_var = new double[_pSBB->_np];
  _LB_obj = -_pSBB->options.INF;
  _LB_var = new double[_pSBB->_np];
  _pc_var = new std::pair<unsigned,double>[2*_pSBB->_np];
  
  // initial points
  for( unsigned i=0; i<_pSBB->_np; i++ ){
    _LB_var[i] = _UB_var[i] = Op<T>::mid( P[i] );
    _pc_var[i] = pc_var[i];
    _pc_var[_pSBB->_np+i] = pc_var[_pSBB->_np+i];
  }
}

template <typename T>
inline
SBBNode<T>::~SBBNode()
{
  delete[] _P0;
  delete[] _P;
  delete[] _UB_var;
  delete[] _LB_var;
  delete[] _pc_var;
}

template <typename T>
inline typename SBB<T>::STATUS
SBBNode<T>::lower_bound
( const double*var, const double inc )
{
  for( unsigned i=0; i<_pSBB->_np; i++ ){
    _P0[i] = _P[i];
    if( var ) _LB_var[i] = var[i]; // initial guess
  }
  return _pSBB->subproblems( SBB<T>::LOWERBD, _pSBB->_np, _P, _LB_var,
    _LB_obj, inc );
}

template <typename T>
inline typename SBB<T>::STATUS
SBBNode<T>::preprocess
( double*var, double&obj, const double inc )
{
  for( unsigned i=0; i<_pSBB->_np; i++ ) _P0[i] = _P[i];
  return _pSBB->subproblems( SBB<T>::PREPROC, _pSBB->_np, _P, var,
    obj, inc );
}

template <typename T>
inline typename SBB<T>::STATUS
SBBNode<T>::postprocess
( double*var, double&obj, const double inc )
{
  for( unsigned i=0; i<_pSBB->_np; i++ ) _P0[i] = _P[i];
  return _pSBB->subproblems( SBB<T>::POSTPROC, _pSBB->_np, _P, var,
    obj, inc );
}

template <typename T>
inline typename SBB<T>::STATUS
SBBNode<T>::upper_bound
( const double*var, const double inc )
{
  for( unsigned i=0; i<_pSBB->_np; i++ ){
    _P0[i] = _P[i];
    if( var ) _UB_var[i] = var[i]; // initial guess
  }
  return _pSBB->subproblems( SBB<T>::UPPERBD, _pSBB->_np, _P, _UB_var,
    _UB_obj, inc );
}

template <typename T>
inline typename SBB<T>::STATUS
SBBNode<T>::test_feasibility
( double*var, double&obj, const double inc )
{
  return _pSBB->subproblems( SBB<T>::FEASTEST, _pSBB->_np, _P, var,
    obj, inc );
}

template <typename T>
inline bool
SBBNode<T>::at_bound
( const double val, const unsigned ip, const double rtol )
const
{
  assert( ip < _pSBB->_np );
  const double margin = Op<T>::diam(_P[ip]) * rtol;
  return( val < Op<T>::l(_P[ip]) + margin
       || val > Op<T>::u(_P[ip]) - margin );
}

template <typename T>
inline void
SBBNode<T>::update_pseudocost
( const typename SBB<T>::TYPE pb )
{
  if( _type == SBB<T>::ROOT || _strongbranch
   || _strength == -_pSBB->options.INF ) return;

  unsigned ipc = _parent.first;
  if( _type == SBB<T>::RIGHT ) ipc += _pSBB->_np;

  _pc_var[ipc].first++;
  switch( pb ){
  case SBB<T>::MINIMIZE:
    if( _LB_obj > -_pSBB->options.INF )
#ifdef USE__SBB_PSEUDOCOST_FILTER
#define PCFILTER 0.5
      _pc_var[ipc].second = ( 1 - PCFILTER ) * _pc_var[ipc].second
        + PCFILTER * ( _LB_obj - _strength )
        / ( Op<T>::diam( _parent.second ) - Op<T>::diam( _P[_parent.first] ) );
#else
      _pc_var[ipc].second += ( _LB_obj - _strength )
        / ( Op<T>::diam( _parent.second ) - Op<T>::diam( _P0[_parent.first] ) );
#endif
#ifdef DEBUG__SBB_PSEUDOCOST
    std::cout << "*** LB0 " << _strength
              << "  LB" << (ipc<_pSBB->_np?"-: ":"+: ") << _LB_obj
              << "  w(P): " << Op<T>::diam( _parent.second )
                             - Op<T>::diam( _P[_parent.first] ) << std::endl;
#endif
    break;
  case SBB<T>::MAXIMIZE:
    if( _UB_obj < _pSBB->options.INF )
#ifdef USE__SBB_PSEUDOCOST_FILTER
#define PCFILTER 0.5
      _pc_var[ipc].second = ( 1 - PCFILTER ) * _pc_var[ipc].second
        - PCFILTER * ( _UB_obj + _strength )
        / ( Op<T>::diam( _parent.second ) - Op<T>::diam( _P[_parent.first] ) );
#else
      _pc_var[ipc].second -= ( _UB_obj + _strength )
        / ( Op<T>::diam( _parent.second ) - Op<T>::diam( _P0[_parent.first] ) );
//        / ( Op<T>::diam( _parent.second ) - Op<T>::diam( _P[_parent.first] ) );
#endif
#ifdef DEBUG__SBB_PSEUDOCOST
    std::cout << "*** UB0 " << -_strength
              << "  UB" << (ipc<_pSBB->_np?"-: ":"+: ") << _UB_obj
              << "  w(P): " << Op<T>::diam( _parent.second )
                             - Op<T>::diam( _P[_parent.first] ) << std::endl;
#endif
    break;
  }
#ifdef DEBUG__SBB_PSEUDOCOST
  std::cout << "*** P " << _parent.first
            << (ipc<_pSBB->_np?"-: ":"+: ")
            << _pc_var[ipc].second/_pc_var[ipc].first
            << "  #" << _pc_var[ipc].first << std::endl;
#endif
}

template <typename T> inline std::ostream&
operator <<
( std::ostream&out, const SBBNode<T>&CV )
{
  std::set<unsigned>::iterator ivar = CV._pSBB->_branch_set.begin();
  for( unsigned i=0; ivar != CV._pSBB->_branch_set.end(); ++ivar, i++ )
    out << "  " << CV._P[i];
  out << std::endl;
  return out;
}

} // namespace mc

#endif
