// Copyright (C) 2012-2014 Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

#ifndef MC__SETINV_HPP
#define MC__SETINV_HPP

#include <utility>
#include <set>
#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <cassert>

#include "mcop.hpp"
#include "mcfunc.hpp"
#include "mclapack.hpp"

#undef MC__SETINV_RGREL // Branching based on absolute variable range only

namespace mc
{

template <typename T> class SetInvNode;
template <typename T> struct lt_SetInvNode;

//! @brief Pure virtual base class for set inversion
////////////////////////////////////////////////////////////////////////
//! mc::SetInv is a pure virtual base C++ class implementing a
//! rigorous set inversion algorithm for nonlinear, implicit functions.
//!
//! REFERENCES:
//! - Jaulin, L., and E. Walter, <A href="http://dx.doi.org/10.1016/0005-1098(93)90106-4">Set inversion via interval analysis for nonlinear bounded-error estimation</A>, <i>Automatica</i>, <b>29</b>(4):1053--1064, 1993.
//! - Lin, Y., and M. A. Stadtherr, <A href="http://dx.doi.org/10.1021/ie0707725"> Guaranteed state and parameter estimation for nonlinear continuous-time systems with bounded-error measurements,</A> <i>Industrial & Engineering Chemistry Research</I>, <B>46</B>:7198--7207, 2007.
//! - Kieffer, M., and E. Walter, <A href="http://dx.doi.org/10.1002/acs.1194">Guaranteed estimation of the parameters of nonlinear continuous-time models: Contributions of interval analysis</A>, <i>Intermation Journal of Adaptive Control & Signal Processing</i>, <b>25</b>:191--207, 2011.
//! .
////////////////////////////////////////////////////////////////////////
template < typename T, typename NODE=SetInvNode<T>, typename LT_NODE=lt_SetInvNode<T> >
class SetInv
////////////////////////////////////////////////////////////////////////
{
  //template <typename U> friend class SetInvNode;

public:
  typedef std::multiset< NODE*, LT_NODE > t_Nodes;
  typedef typename t_Nodes::iterator it_Nodes;
  typedef typename t_Nodes::const_iterator cit_Nodes;

  //! @brief Status of subproblems
  enum STATUS{
    INNER=0,		//!< Current node is inside the inversed set
    OUTER,		//!< Current node is outside of the inversed set
    UNDETERMINED,	//!< Current node is currently undetermined
    FAILURE,		//!< Current node remains undetermined due to failure
    ABORT		//!< Current node terminated after fatal error
  };

  //! @brief Node type
  enum NODETYPE{
    ROOT=0,	//!< Root node
    LEFT=-1,	//!< Left child node
    RIGHT=1	//!< Right child node
  };

  //! @brief Prototype for user-supplied function
  virtual STATUS assess
    ( NODE*	// pointer to current node - possibly modified on return due to bound contraction
    )=0;

  //! @brief Public constructor
  SetInv()
    : options(), _np(0)
    {}

  //! @brief Class destructor
  virtual ~SetInv()
    {
      _clean_stacks();
    }

  //! @brief mc::SetInv options
  struct Options
  {
    //! @brief Constructor
    Options():
      ABSOLUTE_TOLERANCE(1e-5), RELATIVE_TOLERANCE(1e-5),
      BRANCHING_VARIABLE_CRITERION(RGREL), DISPLAY(2), MAX_CPU_TIME(1e6),
      MAX_NODES(0)
      {}
    //! @brief Display
    void display
      ( std::ostream&out ) const;
    //! @brief Branching variable criterion
    enum CRITERION{
      RGREL=0,	//!< Relative variable range diameter
      RGABS	//!< Absolute variable range diameter
    };
    //! @brief Absolute stopping tolerance
    double ABSOLUTE_TOLERANCE;
    //! @brief Relative stopping tolerance
    double RELATIVE_TOLERANCE;
    //! @brief Variable selection criterion
    CRITERION BRANCHING_VARIABLE_CRITERION;
    //! @brief Display option
    int DISPLAY;
    //! @brief Maximum CPU time limit
    double MAX_CPU_TIME;
    //! @brief Maximum number of nodes (0 for no limit)
    unsigned int MAX_NODES;
  } options;

  //! @brief mc::SetInv exceptions
  class Exceptions
  {
  public:
    //! @brief Enumeration type for SetInv exception handling
    enum TYPE{
      BRANCH=0,		//!< Error due to an empty set of branching variables 
      INTERN=-3,	//!< SetInv internal error
      UNDEF=-33		//!< Error due to calling a function/feature not yet implemented in SetInv
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
          return "SetInv Error: empty set of branching variables";
        case INTERN:
          return "SetInv Internal Error";
        case UNDEF: default:
          return "SetInv Error: calling a feature not yet implemented";
        }
      }
  private:
    TYPE _ierr;
  };

  //! @brief Set variables
  void variables
    ( const unsigned int np, const T*P,
      const std::set<unsigned int>*exclude=0 );

  //! @brief Retreive number of variables
  unsigned int np() const
    { return _np; }

  //! @brief Apply set-inversion algorithm - on return, volumes of inner- and boundary-approximations 
  double solve
    ( std::ostream&os=std::cout );

  //! @brief Append all open nodes to <a>os_open</a> - Number of significant digits is set via <a>DPREC</a> (default=6)
  void output_stacks
    ( std::ostream&os_open, const bool PROJ=false, const unsigned int DPREC=6 ) const;

protected:  
  //! @brief Exclusion set for branching variables (e.g. non participating)
  std::set<unsigned int> _exclude_vars;

  //! @brief Compute partition volume
  static double volume
    ( const std::vector<T>&P )
    { double V=1.;
      for( unsigned int i=0; i<P.size(); i++ ) V *= Op<T>::diam(P[i]);
      return V; }

private:  
  //! @brief Number of variables in optimization problem
  unsigned int _np;
  //! @brief Variable bounds at root node
  std::vector<T> _P_root;
  //! @brief Partition volume at root node
  double _volume_root;
  //! @brief Default branching variable set
  std::set<unsigned int> _branch_set;

  //! @brief Variable bounds (temporary)
  std::vector<T> _P_tmp;
  //! @brief Current node index
  unsigned int _node_index;
  //! @brief Total number of nodes introduced so far
  unsigned int _node_cnt;
  //! @brief Max. number of open nodes so far
  unsigned int _node_max;
  //! @brief Set of open nodes
  t_Nodes _open_nodes;
  //! @brief Volume of open nodes
  double _open_volume;

  //! @brief Status after subproblem assessment
  STATUS _status;

  //! @brief Starting time
  double _tstart;
  //! @brief Current time
  double _tcur;

  //! @brief maximum number of values displayed in a row
  static const unsigned int _LDISP = 4;
  //! @brief reserved space for integer variable display
  static const unsigned int _IPREC = 9;
  //! @brief reserved space for double variable display
  static const unsigned int _DPREC = 6;
  //! @brief stringstream for displaying results
  std::ostringstream _odisp;

  //! @brief Erase stored nodes in stacks
  void _clean_stacks();
  //! @brief Reinitialize variables, stacks, counters, time, etc.
  void _restart();
  //! @brief Default set of branching variables
  void _branching_variable_set();

  //! @brief Branch node and create subdomains
  STATUS _branch_node
    ( NODE*pNode );
  //! @brief Selects the branching variable for given node
  std::pair<unsigned int, double> _select_branching_variable
    ( NODE*pNode, const std::set<unsigned int>&set_branch );
  //! @brief Partition branching variable domain for given node
  std::pair<const T, const T> _partition_variable_domain
    ( NODE*pNode, const unsigned int ip ) const;
  //! @brief Determine whether a termination criterion is met
  bool _terminate() const;

  //! @brief Add information for display
  void _display_add
    ( const double dval );
  void _display_add
    ( const unsigned int ival );
  void _display_add
    ( const std::string &sval );

  //! @brief Initialize display
  void _display_init();
  //! @brief Final display
  void _display_final();

  //! @brief Display current buffer stream and reset it
  void _display
    ( std::ostream&os );
  //! @brief Display current open nodes in tree
  void _display_open
    ( std::ostream&os ) const;
  //! @brief Rotate object
  template <typename U>
  U _rotate
    ( const unsigned i, const U*y, const CPPL::dgematrix&A ) const;
  //! @brief Rotate object
  template <typename U>
  U _rotate
    ( const unsigned i, const std::vector<U>&y, const CPPL::dgematrix&A ) const;
};

//! @brief C++ base class for set-inversion nodes
////////////////////////////////////////////////////////////////////////
//! mc::SetInvNode is a C++ base class for defining nodes in the
//! set-inversion algorithm.
////////////////////////////////////////////////////////////////////////
template <typename T>
class SetInvNode
////////////////////////////////////////////////////////////////////////
{
  template <typename U, typename NODE, typename LT_NODE> friend class SetInv;
  template <typename U> friend struct lt_SetInvNode;

public:
  //! @brief Retreive (parent) node strength
  double strength() const
    { return _strength; }
  //! @brief Retreive node index
  unsigned int index() const
    { return _index; }
  //! @brief Retreive node depth
  unsigned int depth() const
    { return _depth; }
  //! @brief Retreive node iteration
  unsigned int iter() const
    { return _iter; }
  //! @brief Retreive/set pointer to user data
  void*& data()
    { return _data; }

  //! @brief Compute partition volume
  double volume() const
    { double V=1.;
      for( unsigned int i=0; i<_pSetInv->np(); i++ ) V *= Op<T>::diam(_P[i]);
      return V; }
  //! @brief Compute partition width
  double width() const
    { double W=0.;
      for( unsigned int i=0; i<_pSetInv->np(); i++ )
        if( W < Op<T>::diam(_P[i]) ) W = Op<T>::diam(_P[i]);
      return W; }

  //! @brief Retreive current variable bounds
  const T& P
    ( const unsigned int ip ) const
    { assert( ip < _pSetInv->np() ); return _P[ip]; }
  //! @brief Retreive/set bound for variable <a>ip</a> 
  T& P
    ( const unsigned int ip )
    { assert( ip < _pSetInv->np() ); return _P[ip]; }
  //! @brief Retreive pointer to variable bounds
  const std::vector<T>& P() const
    { return _P; }
  //! @brief Retreive pointer to variable bounds
  std::vector<T>& P()
    { return _P; }

  //! @brief Retreive/set node type (ROOT/LEFT/RIGHT)
  typename SetInv< T,SetInvNode<T>,lt_SetInvNode<T> >::NODETYPE& type()
    { return _type; }
  //! @brief Typedef for parent nodes
  typedef std::pair<std::vector<T>,CPPL::dgematrix> PARENTNODE;
  //! @brief Retreive/set parent branching variable and range
  std::list<PARENTNODE>& parents()
    { return _parents; }

  //! @brief Retreive const pointer to variable basis
  const CPPL::dgematrix& A() const
    { return _A; }
  //! @brief Retreive/set pointer to variable basis
  CPPL::dgematrix& A()
    { return _A; }

  //! @brief Public constructor (root node)
  SetInvNode
    ( SetInv< T,SetInvNode<T>,lt_SetInvNode<T> >*pSetInv, const std::vector<T>&P,
      const unsigned int index=1, const unsigned int iter=0 );

  //! @brief Public constructor (child node)
  template <typename U> SetInvNode
    ( SetInv< T,SetInvNode<T>,lt_SetInvNode<T> >*pSetInv, const std::vector<T>&P,
      const CPPL::dgematrix&A, const double strength, const unsigned int index,
      const unsigned int depth, const unsigned int iter,
      const typename SetInv< T,SetInvNode<T>,lt_SetInvNode<T> >::NODETYPE type,
      const std::list<PARENTNODE>&parents, U*data );

  //! @brief Destructor
  virtual ~SetInvNode();

  //! @brief Assess
  typename SetInv< T,SetInvNode<T>,lt_SetInvNode<T> >::STATUS assess();

private:
  //! @brief Private default constructor
  SetInvNode<T>()
    {}

  //! @brief Pointer to underlying set-inversion problem
  SetInv< T,SetInvNode<T>,lt_SetInvNode<T> > *_pSetInv;
  //! @brief Stength of parent node
  double _strength;
  //! @brief Depth in SetInv tree
  unsigned int _depth;
  //! @brief Index in SetInv tree
  unsigned int _index;
  //! @brief Iteration when created in SetInv tree
  unsigned int _iter;
  //! @brief Node type (ROOT/LEFT/RIGHT)
  typename SetInv< T,SetInvNode<T>,lt_SetInvNode<T> >::NODETYPE _type;

  //! @brief Variable range and basis of parent nodes
  std::list<PARENTNODE> _parents;
  //! @brief Pointer to user data
  void* _data;

  //! @brief Variable basis
  CPPL::dgematrix _A;
  //! @brief Variable bounds
  std::vector<T> _P;
  //! @brief Backup variable bounds (for domain reduction)
  std::vector<T> _P0;
};

//! @brief C++ structure for comparing SetInv Nodes
////////////////////////////////////////////////////////////////////////
//! mc::lt_SetInvNode is a C++ structure for comparing nodes in branch-and-
//! bound tree based on their lower bound values.
////////////////////////////////////////////////////////////////////////
template <typename T>
struct lt_SetInvNode
{
  bool operator()
    ( const SetInvNode<T>*Node1, const SetInvNode<T>*Node2 ) const
    { return( Node1->strength() > Node2->strength() ); }
};

///////////////////////////////   SetInv   ////////////////////////////////

template <typename T, typename NODE, typename LT_NODE>
inline void
SetInv<T,NODE,LT_NODE>::variables
( const unsigned int np, const T*P, const std::set<unsigned int>*exclude )
{
  if( !np || !P ){
    _np = 0;
    _P_root.clear();
    _P_tmp.clear();
    return;
  }
  
  if( _np != np ){
    _np = np;
    _P_root.resize(_np);
    _P_tmp.resize(_np);
  } 
  for( unsigned int i=0; i<_np; i++ ) _P_root[i] = P[i];

  if( exclude == &_exclude_vars ) return;
  _exclude_vars.clear();
  if( exclude ) _exclude_vars = *exclude;

  return;
}

template <typename T, typename NODE, typename LT_NODE>
template <typename U>
inline U
SetInv<T,NODE,LT_NODE>::_rotate
( const unsigned i, const U*y, const CPPL::dgematrix&A ) const
{
  U xi = 0.;
  for( unsigned j=0; j<_np; j++ ) xi += A(i,j) * y[j];
  return xi;
}

template <typename T, typename NODE, typename LT_NODE>
template <typename U>
inline U
SetInv<T,NODE,LT_NODE>::_rotate
( const unsigned i, const std::vector<U>&y, const CPPL::dgematrix&A ) const
{
  U xi = 0.;
  for( unsigned j=0; j<_np; j++ ) xi += A(i,j) * y[j];
  return xi;
}

template <typename T, typename NODE, typename LT_NODE>
inline void
SetInv<T,NODE,LT_NODE>::output_stacks
( std::ostream&os_open, const bool PROJ, const unsigned DPREC ) const
{
  // append all open sets from _open_nodes to os_open
  if( os_open.good() ){
    os_open << std::scientific << std::setprecision(DPREC); 
    cit_Nodes cit = _open_nodes.begin();
    for( ; cit!=_open_nodes.end(); ++cit ){
      if( PROJ ){
        std::vector<double> p(_np);
        for( int i=0; i<mc::pow(2,_np); i++ ){
          for( unsigned j=0; j<_np; j++ ){
            //std::cout << "  p" << j << (i/mc::pow(2,j)%2?"U":"L");
            p[j] = i/mc::pow(2,j)%2? Op<T>::u((*cit)->P(j)): Op<T>::l((*cit)->P(j));
          }
          //std::cout << std::endl;
          for( unsigned j=0; j<_np; j++ )
            os_open << std::setw(DPREC+9) << _rotate(j,p,CPPL::t((*cit)->A()));
          os_open << std::endl;
        }
        //std::cout << std::endl;
        for( unsigned j=0; j<_np; j++ )
          p[j] = Op<T>::l((*cit)->P(j));
        for( unsigned j=0; j<_np; j++ )
          os_open << std::setw(DPREC+9) << _rotate(j,p,CPPL::t((*cit)->A()));
        os_open << std::endl << std::endl;
      }
      else{
        for( unsigned int ip=0; ip<_np; ip++ )
          os_open << std::setw(DPREC+9) << Op<T>::l( (*cit)->P(ip) )
                  << std::setw(DPREC+9) << Op<T>::u( (*cit)->P(ip) );
        for( unsigned int ip=0; ip<_np; ip++ )
          for( unsigned int jp=0; jp<_np; jp++ )
            os_open << std::setw(DPREC+9) << (*cit)->A()(ip,jp);
        for( unsigned int ip=0; ip<_np; ip++ )
          os_open << std::setw(DPREC+9) << Op<T>::l( _rotate(ip,(*cit)->P(),CPPL::t((*cit)->A())) )
                  << std::setw(DPREC+9) << Op<T>::u( _rotate(ip,(*cit)->P(),CPPL::t((*cit)->A())) );
        os_open << std::endl;
      }
    }
  }
}

template <typename T, typename NODE, typename LT_NODE>
inline double
SetInv<T,NODE,LT_NODE>::solve
( std::ostream&os )
{
  // Create and add root node to set _open_nodes
  _restart();
  _display_init();
  _display( os );
  NODE* rootNode = new NODE( this, _P_root );
  _open_nodes.insert( rootNode );
  _open_volume = _volume_root = rootNode->volume();
  
  // Run iterative set-inversion algorithm
  for( _node_index = 1; !_open_nodes.empty() && _tcur-_tstart < options.MAX_CPU_TIME
       && ( !options.MAX_NODES || _node_index <= options.MAX_NODES );
       _node_index++, _tcur=time() ){

    // Display
    _display_open( os );
    _display_add( _node_index );
    _display_add( (*_open_nodes.begin())->strength() );
    _display_add( (unsigned int)_open_nodes.size() );
    _display_add( _open_volume );

    // Select node and remove it for the set _open_nodes
    NODE* pNode = *_open_nodes.begin();
    _open_nodes.erase( _open_nodes.begin() );
    _open_volume -= pNode->volume();

    // Display
    _display_add( pNode->iter() );

    // Assess node
    _status = pNode->assess();
    switch( _status ){
    case OUTER:
      delete pNode;
      _display_add( "OUTER" );
      break;//continue;

    case INNER:
      delete pNode;
      _display_add( "INNER" );
      break;//continue;

    case UNDETERMINED:
    case FAILURE:
      // Domain branching
      if( _branch_node( pNode ) == ABORT )
        _display_add( "ABORT" );
      delete pNode;
      break;

    case ABORT: default:
      delete pNode;
      _display_add( "ABORT" );
      break;
    }

    _display_add( time()-_tstart );
    _display( os );

    // Keep track of the maximum nodes in memory
    if( _open_nodes.size() > _node_max ) _node_max = _open_nodes.size();

    // Check status and termination criteria
    if( _status == ABORT || _terminate() ) break;
  }

  _display_final();
  _display( os );
  
  return _open_volume;
}

template <typename T, typename NODE, typename LT_NODE>
inline bool
SetInv<T,NODE,LT_NODE>::_terminate() const
{
  if( std::pow(_open_volume,1./_np) < options.ABSOLUTE_TOLERANCE 
   || _open_volume < options.RELATIVE_TOLERANCE * _volume_root )
    return true;
  return false;
}

template <typename T, typename NODE, typename LT_NODE>
inline void
SetInv<T,NODE,LT_NODE>::_restart()
{
  _clean_stacks();
  _branching_variable_set();
  
  _node_index = _node_cnt = _node_max = 0;
  _tstart = _tcur = time();  
}

template <typename T, typename NODE, typename LT_NODE>
inline void
SetInv<T,NODE,LT_NODE>::_clean_stacks()
{
  // Clean open nodes stack
  it_Nodes it = _open_nodes.begin();
  for( ; it != _open_nodes.end(); it++ ) delete *it;
  _open_nodes.clear();
  _open_volume = 0.;
}

template <typename T, typename NODE, typename LT_NODE>
inline void
SetInv<T,NODE,LT_NODE>::_branching_variable_set()
{
  _branch_set.clear();
  for( unsigned int ip=0; ip<_np; ip++ ){
    if( _exclude_vars.empty() || _exclude_vars.find(ip) == _exclude_vars.end() )
      _branch_set.insert( ip );
  }
  if( _branch_set.empty() ) throw Exceptions( Exceptions::BRANCH );
}

template <typename T, typename NODE, typename LT_NODE>
inline std::pair<unsigned int, double>
SetInv<T,NODE,LT_NODE>::_select_branching_variable
( NODE*pNode, const std::set<unsigned int>&set_branch )
{
  std::pair<unsigned int, double> branchsel( _np, -1. ); // <- Can be any negative number
  std::set<unsigned int>::iterator it;

  switch( options.BRANCHING_VARIABLE_CRITERION ){

  // Branching based on relative range diameter
  case Options::RGREL:
    for( it = set_branch.begin(); it != set_branch.end(); ++it ){
      // Compare w.r.t. rotated root node bounds
      double score = Op<T>::diam( pNode->P(*it) )
                   / Op<T>::diam( _rotate( *it, _P_root, pNode->A() ) );
      if( score > branchsel.second ) branchsel = std::make_pair( *it, score );
    }
    break;

  // Branching based on absolute range diameter
  case Options::RGABS:
    for( it = set_branch.begin(); it != set_branch.end(); ++it ){
      double score = Op<T>::diam( pNode->P(*it) );
      if( score > branchsel.second ) branchsel = std::make_pair( *it, score );
    }
    break;
  }

  return branchsel;
}

template <typename T, typename NODE, typename LT_NODE>
inline std::pair<const T, const T>
SetInv<T,NODE,LT_NODE>::_partition_variable_domain
( NODE*pNode, const unsigned int ip_branch )
const
{
  // Branch at mid-point of current variable range
  std::pair<T,T> partition;
  partition.first = mc::Op<T>::l(pNode->P(ip_branch)) + Op<T>::zeroone()
    *Op<T>::diam(pNode->P(ip_branch))/2.;
  partition.second = mc::Op<T>::u(pNode->P(ip_branch)) - Op<T>::zeroone()
    *Op<T>::diam(pNode->P(ip_branch))/2.;
  return partition;
}

template <typename T, typename NODE, typename LT_NODE>
inline typename SetInv<T,NODE,LT_NODE>::STATUS
SetInv<T,NODE,LT_NODE>::_branch_node
( NODE*pNode )
{
  // Branching variable selection
  const unsigned int ip_branch
    = _select_branching_variable( pNode, _branch_set ).first;
  if( _status == ABORT ) return ABORT;

  // Partitionning
  std::pair<const T, const T> partition = _partition_variable_domain( pNode, ip_branch );
  _P_tmp = pNode->P();

  // Left child node
  _P_tmp[ip_branch] = partition.first;
  NODE*pNodeL = new NODE( this, _P_tmp, pNode->A(), pNode->volume(),
    ++_node_cnt, pNode->depth()+1, _node_index, LEFT,
    pNode->parents(), pNode->data() );
  _open_nodes.insert( pNodeL );
  _open_volume += pNodeL->volume();
  
  // Right chold node
  _P_tmp[ip_branch] = partition.second;
  NODE*pNodeR = new NODE( this, _P_tmp, pNode->A(), pNode->volume(),
    ++_node_cnt, pNode->depth()+1, _node_index, RIGHT,
    pNode->parents(), pNode->data() );
  _open_nodes.insert( pNodeR );
  _open_volume += pNodeR->volume();

  // Display
  std::ostringstream omsg;
  omsg << "BRANCH" << ip_branch;
  _display_add( omsg.str() );

  return UNDETERMINED;
}

template <typename T, typename NODE, typename LT_NODE>
inline void
SetInv<T,NODE,LT_NODE>::_display_open
( std::ostream&os ) const
{
  if( options.DISPLAY <= 2 ) return;
  cit_Nodes cit = _open_nodes.begin();
  for( unsigned int id=1; cit != _open_nodes.end(); ++cit, id++ ){
    std::set<unsigned int>::iterator ivar = _branch_set.begin();
    os << "Node " << id << ":" << std::scientific << std::setprecision(_DPREC)
       << std::setw(_DPREC+8) << (*cit)->strength();
    for( ; ivar != _branch_set.end(); ++ivar )
      os << (*cit)->P(*ivar) << "  ";
    os << std::endl;
  }
  return;
}

template <typename T, typename NODE, typename LT_NODE>
inline void
SetInv<T,NODE,LT_NODE>::_display_init()
{
  _odisp.str("");
  if( options.DISPLAY <= 1 ) return;
  _odisp << std::right
  	 << std::setw(_IPREC) << "INDEX"
  	 << std::setw(_DPREC+8) << "MEASURE  "
  	 << std::setw(_IPREC) << "STACK"
  	 << std::setw(_DPREC+8) << "OPEN    "
  	 << std::setw(_IPREC) << "PARENT"
  	 << std::setw(_DPREC+8) << "ACTION"
  	 << std::setw(_DPREC+8) << "CPU TIME "
  	 << std::endl;
  
}

template <typename T, typename NODE, typename LT_NODE>
inline void
SetInv<T,NODE,LT_NODE>::_display_add
( const double dval )
{
  if( options.DISPLAY <= 1 ) return;
  _odisp << std::right << std::scientific << std::setprecision(_DPREC)
         << std::setw(_DPREC+8) << dval;
}

template <typename T, typename NODE, typename LT_NODE>
inline void
SetInv<T,NODE,LT_NODE>::_display_add
( const unsigned int ival )
{
  if( options.DISPLAY <= 1 ) return;
  _odisp << std::right << std::setw(_IPREC) << ival;
}

template <typename T, typename NODE, typename LT_NODE>
inline void
SetInv<T,NODE,LT_NODE>::_display_add
( const std::string &sval )
{
  if( options.DISPLAY <= 1 ) return;
  _odisp << std::right << std::setw(_DPREC+8) << sval;
}

template <typename T, typename NODE, typename LT_NODE>
inline void
SetInv<T,NODE,LT_NODE>::_display_final()
{
  if( options.DISPLAY <= 0 ) return;

  // Solution found within allowed time?
  if( _tcur-_tstart < options.MAX_CPU_TIME
    && ( !options.MAX_NODES || _node_index <= options.MAX_NODES ) )
    _odisp << std::endl << "#  NORMAL TERMINATION:      ";
  else
    _odisp << std::endl << "#  EXECUTION STOPPED:       ";
  _odisp << std::fixed << std::setprecision(6) << _tcur-_tstart << " CPU SEC"
         << std::endl;

  // Set-inversion results
  _odisp << "#  BOUNDARY APPROXIMATION:  "
  	 << std::scientific << std::setprecision(_DPREC)
         << "MEASURE = " << std::setw(_DPREC+6) << _open_volume
         << ",  NODES = "  << _open_nodes.size() 
  	 << std::endl;
  _odisp << "#  TOTAL NUMBER OF NODES:   " << _node_index-1 << std::endl
         << "#  MAX OPEN NODES IN STACK: " << _node_max << std::endl;
}

template <typename T, typename NODE, typename LT_NODE>
inline void
SetInv<T,NODE,LT_NODE>::_display
( std::ostream &os )
{
  if( _odisp.str() == "" ) return;
  //if( options.DISPLAY > 0 ){
  os << _odisp.str() << std::endl;
  //}
  _odisp.str("");
  return;
}

template <typename T, typename NODE, typename LT_NODE>
inline void
SetInv<T,NODE,LT_NODE>::Options::display
( std::ostream&out ) const
{
  // Display SetInv Options
  out << std::setw(60) << "  ABSOLUTE CONVERGENCE TOLERANCE"
      << std::scientific << std::setprecision(1)
      << ABSOLUTE_TOLERANCE << std::endl;
  out << std::setw(60) << "  RELATIVE CONVERGENCE TOLERANCE"
      << std::scientific << std::setprecision(1)
      << RELATIVE_TOLERANCE << std::endl;
  out << std::setw(60) << "  MAXIMUM NUMBER OF ITERATIONS (PARTITIONS)";
  switch( MAX_NODES ){
    case 0:  out << "N/A\n"; break;
    default: out << MAX_NODES << std::endl; break;
  }
  out << std::setw(60) << "  MAXIMUM CPU TIME (SEC)"
      << std::scientific << std::setprecision(1)
      << MAX_CPU_TIME << std::endl;
  out << std::setw(60) << "  DISPLAY LEVEL"
      << DISPLAY << std::endl;
  out << std::setw(60) << "  BRANCHING STRATEGY FOR VARIABLE SELECTION";
#ifdef MC__SETINV_RGREL
  switch( BRANCHING_VARIABLE_CRITERION ){
    case RGREL:  out << "RGREL\n";  break;
    case RGABS:  out << "RGABS\n";  break;
  }
#else
  out << "RGABS\n";
#endif
}

/////////////////////////////// SetInvNode ///////////////////////////////

template <typename T>
inline
SetInvNode<T>::SetInvNode
( SetInv< T,SetInvNode<T>,lt_SetInvNode<T> >*pSetInv, const std::vector<T>&P,
  const unsigned int index, const unsigned int iter ):
_pSetInv( pSetInv ), _depth( 0 ), _index( index ), _iter( iter ),
_type( SetInv< T,SetInvNode<T>,lt_SetInvNode<T> >::ROOT ), _data(0)
{
  _P = P;
  _A.resize(_pSetInv->np(),_pSetInv->np()).identity();
  _strength = volume();
}

template <typename T> template <typename U>
inline
SetInvNode<T>::SetInvNode
( SetInv< T,SetInvNode<T>,lt_SetInvNode<T> >*pSetInv, const std::vector<T>&P,
  const CPPL::dgematrix&A, const double strength, const unsigned int index,
  const unsigned int depth, const unsigned int iter,
  const typename SetInv< T,SetInvNode<T>,lt_SetInvNode<T> >::NODETYPE type,
  const std::list<PARENTNODE>&parents, U*data ):
_pSetInv( pSetInv ), _strength( strength ), _depth( depth ), _index( index ),
_iter( iter ), _type( type ), _parents( parents ), _data( data ), _A( A ), _P( P )
{}

template <typename T>
inline
SetInvNode<T>::~SetInvNode()
{}

template <typename T>
inline typename SetInv< T,SetInvNode<T>,lt_SetInvNode<T> >::STATUS
SetInvNode<T>::assess()
{
  _P0 = _P;
  return _pSetInv->assess( this );
}

} // namespace mc

#endif
