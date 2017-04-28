// Copyright (C) 2012-2016 Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

#ifndef MC__SBP_HPP
#define MC__SBP_HPP

#include <utility>
#include <list>
#include <set>
#include <map>
#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <cassert>

#include "mcop.hpp"
#include "mctime.hpp"
#include "mclapack.hpp"

#undef  DEBUG__SBP_STRONGBRANCHING
//#undef  DEBUG__SBP_SCOREBRANCHING

namespace mc
{

template <typename T> class SBPNode;
template <typename T> struct lt_SBP;

// For block decomposition
extern "C" void mc13d_
  ( const int*, const int*, const int*, const int*, const int*, int*, int*, int*, int* );

//! @brief Pure virtual base class for branch-and-prune
////////////////////////////////////////////////////////////////////////
//! mc::SBP is a virtual base C++ class implementing a branch-and-prune
//! algorithm, also known as set-inversion or constraint-projection, in
//! order to find the pre-image of a nonlinear mapping.
//!
//! REFERENCES:
//! - Jaulin, L., and E. Walter, <A href="http://dx.doi.org/10.1016/0005-1098(93)90106-4">Set inversion via interval analysis for nonlinear bounded-error estimation</A>, <i>Automatica</i>, <b>29</b>(4):1053--1064, 1993.
//! - Lin, Y., and M. A. Stadtherr, <A href="http://dx.doi.org/10.1021/ie0707725"> Guaranteed state and parameter estimation for nonlinear continuous-time systems with bounded-error measurements,</A> <i>Industrial & Engineering Chemistry Research</I>, <B>46</B>:7198--7207, 2007.
//! - Baharev, A., and E. Rév, <A href="https://www.researchgate.net/profile/Ali_Baharev/publication/237677929_A_complete_nonlinear_system_solver_using_affine/links/0a85e5382708b1dfab000000.pdf#page=24">A complete nonlinear system solver using affine arithmetic</A>, International Workshop on Interval Analysis and Constraint Propagation for Applications (IntCP), 2009.
//! - Kieffer, M., and E. Walter, <A href="http://dx.doi.org/10.1002/acs.1194">Guaranteed estimation of the parameters of nonlinear continuous-time models: Contributions of interval analysis</A>, <i>Intermation Journal of Adaptive Control & Signal Processing</i>, <b>25</b>:191--207, 2011.
//! .
////////////////////////////////////////////////////////////////////////
template < typename T, typename NODE=SBPNode<T>, typename LT_NODE=lt_SBP<T> >
class SBP
////////////////////////////////////////////////////////////////////////
{
  //template <typename U> friend class SBPNode;

public:
  typedef std::multiset< NODE*, LT_NODE > t_Nodes;
  typedef typename t_Nodes::iterator it_Nodes;
  typedef typename t_Nodes::const_iterator cit_Nodes;

  //! @brief Status of subproblems
  enum STATUS{
    INNER=0,		//!< Current node is inside the inversed set
    OUTER,		    //!< Current node is outside of the inversed set
    UNDETERMINED,	//!< Current node is currently undetermined
    FAILURE,		//!< Current node remains undetermined due to failure
    ABORT		    //!< Current node terminated after fatal error
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
  SBP()
    : options(), _np(0)
    {}

  //! @brief Class destructor
  virtual ~SBP()
    {
      _clean_stacks();
    }

  //! @brief mc::SBP options
  struct Options
  {
    //! @brief Constructor
    Options():
      MEASURE(RELMAXLEN), STOPPING_TOLERANCE(1e-5),
      BRANCHING_VARIABLE_CRITERION(RGREL), BRANCHING_USERFUNCTION(0),
      STRONG_BRANCHING_MAXDEPTH(0), STRONG_BRANCHING_RELTOL(1e-3),
      STRONG_BRANCHING_ABSTOL(1e-3), SCORE_BRANCHING_USE(false),
      SCORE_BRANCHING_RELTOL(1e-1), SCORE_BRANCHING_ABSTOL(1e-2),
      SCORE_BRANCHING_MAXSIZE(0), DISPLAY(2), OUTPUT_ROTATED(false),
      MAX_CPUTIME(1e6), MAX_NODES(0)
      {}
    //! @brief Display
    void display
      ( std::ostream&out ) const;
    //! @brief Branching variable criterion
    enum CRITERION{
      RGREL=0,	//!< Relative variable range diameter
      RGABS	//!< Absolute variable range diameter
    };
    //! @brief Branching variable criterion - Relative is w.r.t. the root node 
    enum METRIC{
      ABSMAXLEN=0,	//!< Absolute max length of partition elements
      RELMAXLEN,	//!< Relative max length of partition elements
      ABSMEANLEN,	//!< Absolute mean length of partition elements (corresponding to width of a cube of equal volume)
      RELMEANLEN,	//!< Relative mean length of partition elements (corresponding to width of a cube of equal volume)
      ABSVOLUME,	//!< Absolute volume of partition elements
      RELVOLUME		//!< Relative volume of partition elements
    };
    //! @brief Branching selection user-function
    typedef std::set<unsigned> (*SELECTION)( const NODE* );
    //! @brief Partition measure
    int MEASURE;
    //! @brief Stopping tolerance in the chosen metric
    double STOPPING_TOLERANCE;
    //! @brief Branching variable selection criterion
    int BRANCHING_VARIABLE_CRITERION;
    //! @brief Branching variable selection user-function
    SELECTION BRANCHING_USERFUNCTION;
    //! @brief Maximum depth for strong branching interruption
    unsigned STRONG_BRANCHING_MAXDEPTH;
    //! @brief Relative tolerance for strong branching interruption
    double STRONG_BRANCHING_RELTOL;
    //! @brief Absolute tolerance for strong branching interruption
    double STRONG_BRANCHING_ABSTOL;
    //! @brief Whether to base branching variable selection on scores
    bool SCORE_BRANCHING_USE;
    //! @brief Relative tolerance for branching variable selection based on scores
    double SCORE_BRANCHING_RELTOL;
    //! @brief Absolute tolerance for branching variable selection based on scores
    double SCORE_BRANCHING_ABSTOL;
    //! @brief Maximal size of branching variable selection based on scores
    unsigned SCORE_BRANCHING_MAXSIZE;
    //! @brief Display option
    int DISPLAY;
    //! @brief Whether to display rotated boxes
    bool OUTPUT_ROTATED;
    //! @brief Maximum CPU time limit
    double MAX_CPUTIME;
    //! @brief Maximum number of nodes (0 for no limit)
    unsigned MAX_NODES;
  };
  mutable Options options;

  //! @brief mc::SBP exceptions
  class Exceptions
  {
  public:
    //! @brief Enumeration type for SBP exception handling
    enum TYPE{
      BRANCH=0,		//!< Error due to an empty set of branching variables 
      INTERN=-3,	//!< SBP internal error
      UNDEF=-33		//!< Error due to calling a function/feature not yet implemented in SBP
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
          return "SBP Error: empty set of branching variables";
        case INTERN:
          return "SBP Internal Error";
        case UNDEF: default:
          return "SBP Error: calling a feature not yet implemented";
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

  //! @brief Retreive root variable bounds
  const T& P_root
    ( const unsigned ip ) const
    { assert( ip < _np ); return _P_root[ip]; }

  //! @brief Apply set-inversion algorithm - on return, measure of set-boundary approximations 
  double solve
    ( std::ostream&os=std::cout, const bool resume=false );

  //! @brief Cluster boxes in set-boundary approximation - on return, list of box enclosures for clusters
  std::list<T*> clusters
    ( const double rtol=1e-4, const double atol=1e-8, std::ostream&os=std::cout ) const;

  //! @brief Append all open nodes to <a>os_open</a> - Number of significant digits is set via <a>DPREC</a> (default=6)
  void output_nodes
    ( std::ostream&os_open, const bool PROJ=false, const unsigned DPREC=6 ) const
    {
      _output_nodes( os_open,  _open_nodes,  PROJ, DPREC );
    }

  //! @brief Append all open nodes to <a>os_open</a> and inner nodes to <a>os_inner</a> - Number of significant digits is set via <a>DPREC</a> (default=6)
  void output_nodes
    ( std::ostream&os_open, std::ostream&os_inner, const bool PROJ=false,
      const unsigned DPREC=6 ) const
    {
      _output_nodes( os_open,  _open_nodes,  PROJ, DPREC );
      _output_nodes( os_inner, _inner_nodes, PROJ, DPREC );
    }

  //! @brief Return the current stack of boundary (open) nodes
  const t_Nodes& open_nodes() const
    { return _open_nodes; };

  //! @brief Return the current stack of inner nodes
  const t_Nodes& inner_nodes() const
    { return _inner_nodes; };

  //! @brief Append all clusters to <a>os_clus</a> - Number of significant digits is set via <a>DPREC</a> (default=6)
  void output_clusters
    ( std::ostream&os_clus, const double RTOL=1e-4, const double ATOL=1e-8,
      const unsigned DPREC=6 ) const;

  //! @brief Measure current partition
  double measure
    ( const NODE*pNode ) const
    { unsigned dim(0); 
      switch( options.MEASURE ){
        case Options::ABSVOLUME:          return volume( pNode, dim, false );
        case Options::RELVOLUME:          return volume( pNode, dim, true );
        case Options::ABSMEANLEN:         return meanwidth( pNode, false );
        case Options::RELMEANLEN:         return meanwidth( pNode, true );
        case Options::ABSMAXLEN:          return maxwidth( pNode, false );
        case Options::RELMAXLEN: default: return maxwidth( pNode, true );
      } }

protected:  
  //! @brief Exclusion set for branching variables (e.g. non participating)
  std::set<unsigned> _exclude_vars;
  //! @brief Set of open nodes
  t_Nodes _open_nodes;
  //! @brief Set of inner nodes
  t_Nodes _inner_nodes;
  //! @brief Curent measure of open nodes
  double _open_measure;
  //! @brief Variable bounds at root node
  std::vector<T> _P_root;

  //! @brief Compute element volume
  virtual double volume
    ( const NODE*pNode, unsigned &dim, const bool rel=true ) const
    { double V=1.;
      dim = _np;
      for( unsigned i=0; i<dim; i++ ){
        double scal = rel? Op<T>::diam(pNode->P0(i)): 1;
        V *= Op<T>::diam(pNode->P(i)) / scal;
      }
      return V; }
  //! @brief Compute element mean width (equivalent cube edge)
  virtual double meanwidth
    ( const NODE*pNode, const bool rel=true ) const
    { unsigned dim(0);
      double V = volume( pNode, dim, rel );
      return std::pow( V, 1./dim ); }
  //! @brief Compute element max width
  virtual double maxwidth
    ( const NODE*pNode, const bool rel=true ) const
    { double W=0.;
      for( unsigned i=0; i<_np; i++ ){
        double scal = rel? Op<T>::diam(pNode->P0(i)): 1;
        if( Op<T>::diam(pNode->P(i)) > W * scal )
          W = Op<T>::diam(pNode->P(i)) / scal;
      }
      return W; }

  //! @brief Append all nodes in <a>nodes</a> to <a>os</a> - Number of significant digits is set via <a>DPREC</a> (default=6)
  void _output_nodes
    ( std::ostream&os, const t_Nodes&nodes, const bool PROJ, const unsigned DPREC ) const;

private:  
  //! @brief Number of variables in optimization problem
  unsigned _np;
  //! @brief Partition measure at root node
  //double _root_measure;
  //! @brief Default branching variable set
  std::set<unsigned> _branch_set;
  //! @brief Current branching variable
  unsigned _branch_var;

  //! @brief Variable bounds (temporary)
  std::vector<T> _P_tmp;
  //! @brief Current node index
  unsigned _node_index;
  //! @brief Total number of nodes introduced so far
  unsigned _node_count;
  //! @brief Max. number of open nodes so far
  unsigned _node_max;

  //! @brief Status after subproblem assessment
  STATUS _status;

  //! @brief Starting time
  double _tstart;

  //! @brief maximum number of values displayed in a row
  static const unsigned _LDISP = 4;
  //! @brief reserved space for integer variable display
  static const unsigned _IPREC = 9;
  //! @brief reserved space for double variable display
  static const unsigned _DPREC = 6;
  //! @brief stringstream for displaying results
  std::ostringstream _odisp;

  //! @brief Erase stored nodes in stacks
  void _clean_stacks();
  //! @brief Reinitialize variables, stacks, counters, time, etc.
  void _restart();
  //! @brief Determine whether a termination criterion is met
  bool _terminate() const;

  //! @brief Default set of branching variables
  void _branching_variable_set
    ( const NODE*pNode );
  //! @brief Subset of branching variables based on scores
  void _branching_score_subset
    ( const NODE*pNode );
  //! @brief Branch node by selecting branching variable and creating subdomains
  std::pair<NODE*,NODE*> _branch_node
    ( NODE*pNode );
  //! @brief Branch node by bisecting along variable ivar and creating subdomains
  std::pair<NODE*,NODE*> _branch_node
    ( NODE*pNode, const unsigned ivar, unsigned&node_count );
  //! @brief Partition node domain along variable ivar
  std::pair<const T, const T> _partition_variable_domain
    ( NODE*pNode, const unsigned ivar ) const;
  //! @brief Selects the branching variable for given node
  std::pair<unsigned, double> _select_branching_variable
    ( NODE*pNode );
  //! @brief Apply strong branching for branching variable selection
  std::pair<unsigned, double> _strong_branching
    ( NODE*pNode, unsigned&var_branch, unsigned&node_count );

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
  void _display_final
    ( const double tcur );

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
  //! @brief Test range adjacency
  bool _adjacent
    ( const NODE*Node1, const NODE*Node2, const double rtol, const double atol ) const;
};

//! @brief C++ base class for set-inversion nodes
////////////////////////////////////////////////////////////////////////
//! mc::SBPNode is a C++ base class for defining nodes in the
//! set-inversion algorithm.
////////////////////////////////////////////////////////////////////////
template <typename T>
class SBPNode
////////////////////////////////////////////////////////////////////////
{
  template <typename U, typename NODE, typename LT_NODE> friend class SBP;
  template <typename U> friend struct lt_SBP;

public:
  //! @brief Retreive node measure
  double measure() const
    { return _measure; }
  //! @brief Retreive node index
  unsigned index() const
    { return _index; }
  //! @brief Retreive node depth
  unsigned depth() const
    { return _depth; }
  //! @brief Retreive node iteration
  unsigned iter() const
    { return _iter; }
  //! @brief Retreive/set pointer to user data
  void*& data()
    { return _data; }
  //! @brief Retreive/set node dependent variables
  std::set<unsigned>& depend()
    { return _depend; }
  //! @brief Retreive node dependent variables
  const std::set<unsigned>& depend() const
    { return _depend; }
  //! @brief Retreive/set node branching scores
  std::map<unsigned,double>& scores()
    { return _scores; }
  //! @brief Retreive node branching scores
  const std::map<unsigned,double>& scores() const
    { return _scores; }

  //! @brief Retreive root variable bounds
  const T& P0
    ( const unsigned ip ) const
    { assert( ip < _pSBP->np() ); return _pSBP->P_root(ip); }
  //! @brief Retreive current variable bounds
  const T& P
    ( const unsigned ip ) const
    { assert( ip < _pSBP->np() ); return _P[ip]; }
  //! @brief Retreive/set bound for variable <a>ip</a> 
  T& P
    ( const unsigned ip )
    { assert( ip < _pSBP->np() ); return _P[ip]; }
  //! @brief Retreive pointer to variable bounds
  const std::vector<T>& P() const
    { return _P; }
  //! @brief Retreive pointer to variable bounds
  std::vector<T>& P()
    { return _P; }
  //! @brief Output current node bounds in <a>os_nodes</a>    
  void output
    ( std::ostream& os_nodes, bool use=false ) const;

  //! @brief Retreive/set node type (ROOT/LEFT/RIGHT)
  typename SBP< T,SBPNode<T>,lt_SBP<T> >::NODETYPE& type()
    { return _type; }
  //! @brief Retreive status
  typename SBP< T,SBPNode<T>,lt_SBP<T> >::STATUS status() const
    { return _status; }
  //! @brief Retreive/set parent branching variable and range
  std::pair<unsigned,T>& parent()
    { return _parent; }
  //! @brief Retreive/set strong branching status
  bool& strongbranch()
    { return _strongbranch; }
  //! @brief Retreive/set converged status
  bool& converged()
    { return _converged; }

  //! @brief Retreive const pointer to variable basis
  const CPPL::dgematrix& A() const
    { return _A; }
  //! @brief Retreive/set pointer to variable basis
  CPPL::dgematrix& A()
    { return _A; }

  //! @brief Public constructor (root node)
  SBPNode
    ( SBP< T,SBPNode<T>,lt_SBP<T> >*pSBP, const std::vector<T>&P,
      const unsigned index=1, const unsigned iter=0 );

  //! @brief Public constructor (child node)
  template <typename U> SBPNode
    ( SBP< T,SBPNode<T>,lt_SBP<T> >*pSBP, const std::vector<T>&P,
      const CPPL::dgematrix&A, const unsigned index, const unsigned depth,
      const unsigned iter, const std::set<unsigned>&depend,
      const typename SBP<T,SBPNode<T>,lt_SBP<T>>::NODETYPE type,
      const std::pair<unsigned,T>&parent, const bool converged, U*data );

  //! @brief Destructor
  virtual ~SBPNode();

  //! @brief Assess
  typename SBP< T,SBPNode<T>,lt_SBP<T> >::STATUS assess();

private:
  //! @brief Private default constructor
  SBPNode<T>()
    {}

  //! @brief Pointer to underlying set-inversion problem
  SBP< T,SBPNode<T>,lt_SBP<T> > *_pSBP;
  //! @brief Node measure
  double _measure;
  //! @brief Depth in SBP tree
  unsigned _depth;
  //! @brief Index in SBP tree
  unsigned _index;
  //! @brief Iteration when created in SBP tree
  unsigned _iter;
  //! @brief Node type (ROOT/LEFT/RIGHT)
  typename SBP< T,SBPNode<T>,lt_SBP<T> >::NODETYPE _type;
  //! @brief Node status
  typename SBP< T,SBPNode<T>,lt_SBP<T> >::STATUS _status;

  //! @brief Subset of dependent variables in node
  std::set<unsigned> _depend;
  //! @brief Map of variables with scores for branching variable selection
  std::map<unsigned,double> _scores;
  //! @brief Parent branching variable and range
  std::pair<unsigned,T> _parent;
  //! @brief Strong branching status
  bool _strongbranch;
  //! @brief Converged status
  bool _converged;
  //! @brief Pointer to user data
  void* _data;

  //! @brief Variable basis
  CPPL::dgematrix _A;
  //! @brief Variable bounds
  std::vector<T> _P;
};

//! @brief C++ structure for comparing SBP Nodes
////////////////////////////////////////////////////////////////////////
//! mc::lt_SBP is a C++ structure for comparing nodes in branch-and-
//! bound tree based on their lower bound values.
////////////////////////////////////////////////////////////////////////
template <typename T>
struct lt_SBP
{
  bool operator()
    ( const SBPNode<T>*Node1, const SBPNode<T>*Node2 ) const
    //{ return( Node1->measure() > Node2->measure() ); }
    { if( Node1->measure() > Node2->measure() ) return true;
      if( Node1->measure() < Node2->measure() ) return false;
      return( Node1->depth() < Node2->depth() ); }
};

///////////////////////////////   SBP   ////////////////////////////////

template <typename T, typename NODE, typename LT_NODE>
inline void
SBP<T,NODE,LT_NODE>::variables
( const unsigned np, const T*P, const std::set<unsigned>&exclude )
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
  for( unsigned i=0; i<_np; i++ ) _P_root[i] = P[i];

  if( &exclude == &_exclude_vars ) return;
  _exclude_vars = exclude;
//#ifdef DEBUG__SBP_SHOW_EXCLUDED
  std::cout << "*** EXCLUDED VARIABLES: {";
  for( auto it=_exclude_vars.begin(); it!=_exclude_vars.end(); ++it )
    std::cout << " " << *it;
  std::cout << " }" << std::endl;
//#endif

  return;
}

template <typename T, typename NODE, typename LT_NODE>
template <typename U>
inline U
SBP<T,NODE,LT_NODE>::_rotate
( const unsigned i, const U*y, const CPPL::dgematrix&A ) const
{
  U xi = 0.;
  for( unsigned j=0; j<_np; j++ ) xi += A(i,j) * y[j];
  return xi;
}

template <typename T, typename NODE, typename LT_NODE>
template <typename U>
inline U
SBP<T,NODE,LT_NODE>::_rotate
( const unsigned i, const std::vector<U>&y, const CPPL::dgematrix&A ) const
{
  U xi = 0.;
  for( unsigned j=0; j<_np; j++ ) xi += A(i,j) * y[j];
  return xi;
}

template <typename T, typename NODE, typename LT_NODE>
inline bool
SBP<T,NODE,LT_NODE>::_adjacent
( const NODE*Node1, const NODE*Node2, const double rtol, const double atol ) const
{
  for( unsigned k=0; k<_np; k++ ){
    T P1k = _rotate( k, Node1->P(), Node1->A() );
    T P2k = _rotate( k, Node2->P(), Node2->A() );
    if( 2.*std::fabs(Op<T>::mid(P1k)-Op<T>::mid(P2k))
        > (Op<T>::diam(P1k)+Op<T>::diam(P2k))*(1.+rtol)+atol ) return false;
  }
  return true;
}

template <typename T, typename NODE, typename LT_NODE>
inline std::list<T*>
SBP<T,NODE,LT_NODE>::clusters
( const double RTOL, const double ATOL, std::ostream&os ) const
{
  // Setup indidence matrix
  const int N = _open_nodes.size();
  std::vector<int> ICN, IP(N), LENR(N);
  std::vector<const T*> range_nodes(N);
  auto it1 = _open_nodes.begin();
  for( int i=0; it1!=_open_nodes.end(); ++it1, i++ ){
    IP[i] = ICN.size()+1;
    auto it2 = _open_nodes.begin();
    for( int j=0; it2!=_open_nodes.end(); ++it2, j++ )
      if( i==j || _adjacent( *it1, *it2, RTOL, ATOL ) ) ICN.push_back(j+1);
    LENR[i] = ICN.size() - IP[i]+1;
    range_nodes[i] = (*it1)->P().data();
  }
  const int LICN = ICN.size();

  // Permute incidence matrix to lower triangular form using MC13
  std::vector<int> IOR(N), IB(N), IW(3*N);
  int NUM;
  mc13d_( &N, ICN.data(), &LICN, IP.data(), LENR.data(), IOR.data(),
          IB.data(), &NUM, IW.data() );

  // Determine clusters
  std::list<T*> L_clus;
  for( int iclus=0; iclus<NUM; iclus++ ){
    T* I_clus = new T[_np];
    int pclus = IB[iclus]-1;
    for( unsigned ip=0; ip<_np; ip++ ) // initialize cluster
      I_clus[ip] = range_nodes[IOR[pclus]-1][ip];
    int uclus = (iclus+1<NUM? IB[iclus+1]-1: N);
    for( ; pclus<uclus; pclus++ )
      for( unsigned ip=0; ip<_np; ip++ ) // expand cluster
        I_clus[ip] = Op<T>::hull( I_clus[ip], range_nodes[IOR[pclus]-1][ip] );
    L_clus.push_back( I_clus );
  }
  return L_clus;
}

template <typename T, typename NODE, typename LT_NODE>
inline void
SBP<T,NODE,LT_NODE>::output_clusters
( std::ostream&os_clus, const double RTOL, const double ATOL,
  const unsigned DPREC ) const
{
  // Determine clusters
  const std::list<T*> L_clus = clusters( RTOL, ATOL );

  // append all clusters from L_clus to os_clus
  if( os_clus.good() ){
    os_clus << std::scientific << std::setprecision(DPREC); 
    auto cit = L_clus.begin();
    for( ; cit!=L_clus.end(); ++cit ){
      for( unsigned ip=0; ip<_np; ip++ )
        os_clus << std::setw(DPREC+9) << Op<T>::l( (*cit)[ip] )
                << std::setw(DPREC+9) << Op<T>::u( (*cit)[ip] );
      os_clus << std::endl;
    }
  }

  // Clean-up
  for( auto it=L_clus.begin(); it!=L_clus.end(); ++it ) delete[] *it;

}

template <typename T, typename NODE, typename LT_NODE>
inline void
SBP<T,NODE,LT_NODE>::_output_nodes
( std::ostream&os, const t_Nodes&nodes, const bool PROJ, const unsigned DPREC ) const
{
  if( !os.good() ) return;

  os << std::scientific << std::setprecision(DPREC); 
  cit_Nodes cit = nodes.begin();
  for( ; cit!=nodes.end(); ++cit ){
    if( PROJ ){
      std::vector<double> p(_np);
      for( int i=0; i<mc::pow(2,_np); i++ ){
        for( unsigned j=0; j<_np; j++ ){
          //std::cout << "  p" << j << (i/mc::pow(2,j)%2?"U":"L");
          p[j] = i/mc::pow(2,j)%2? Op<T>::u((*cit)->P(j)): Op<T>::l((*cit)->P(j));
        }
        //std::cout << std::endl;
        for( unsigned j=0; j<_np; j++ )
          os << std::setw(DPREC+9) << _rotate(j,p,CPPL::t((*cit)->A()));
        os << std::endl;
      }
      //std::cout << std::endl;
      for( unsigned j=0; j<_np; j++ )
        p[j] = Op<T>::l((*cit)->P(j));
      for( unsigned j=0; j<_np; j++ )
        os << std::setw(DPREC+9) << _rotate(j,p,CPPL::t((*cit)->A()));
      os << std::endl << std::endl;
    }
    else{
      for( unsigned ip=0; ip<_np; ip++ )
        os << std::setw(DPREC+9) << Op<T>::l( (*cit)->P(ip) )
                << std::setw(DPREC+9) << Op<T>::u( (*cit)->P(ip) );
      for( unsigned ip=0; options.OUTPUT_ROTATED && ip<_np; ip++ )
        for( unsigned jp=0; jp<_np; jp++ )
          os << std::setw(DPREC+9) << (*cit)->A()(ip,jp);
      for( unsigned ip=0; options.OUTPUT_ROTATED && ip<_np; ip++ )
        os << std::setw(DPREC+9) << Op<T>::l( _rotate(ip,(*cit)->P(),CPPL::t((*cit)->A())) )
                << std::setw(DPREC+9) << Op<T>::u( _rotate(ip,(*cit)->P(),CPPL::t((*cit)->A())) );
      os << std::endl;
    }
  }
}

template <typename T, typename NODE, typename LT_NODE>
inline double
SBP<T,NODE,LT_NODE>::solve
( std::ostream&os, const bool resume )
{
  // Display header
  _display_init();
  _display( os );

  if( !resume ){
    // Create and assess root node
    _restart();
    NODE* rootNode = new NODE( this, _P_root );
     _status = rootNode->assess();
    //for( unsigned i=0; i<_np; i++ ) _P_root[i] = rootNode->P()[i];
    //_P_root = rootNode->P();
    //_root_measure = measure( _P_root );
    _open_measure = measure( rootNode );
    _open_nodes.insert( rootNode );
  }

  // Run iterative set-inversion algorithm
  for( ; !_open_nodes.empty()
         && cpuclock()-_tstart < options.MAX_CPUTIME
         && ( !options.MAX_NODES || _node_index <= options.MAX_NODES );
       _node_index++ ){

    // Select node, remove it for the set _open_nodes, and test termination criterion
    NODE* pNode = *_open_nodes.begin();
    _open_measure = pNode->measure();
    if( _terminate() ) break;
    _open_nodes.erase( _open_nodes.begin() );

    // Display node
    _display_add( _node_index );
    _display_add( pNode->measure() );
    _display_add( (unsigned)_open_nodes.size() );
    _display_add( pNode->iter() );
/*
    switch( pNode->status() ){
     case OUTER:
      _open_nodes.erase( _open_nodes.begin() );
      delete pNode;
      _display_add( cpuclock()-_tstart );
      _display_add( "OUTER" );
      _display( os );
      continue;
     case INNER:
      _open_nodes.erase( _open_nodes.begin() );
      _inner_nodes.insert( pNode );
      //delete pNode;
      _display_add( cpuclock()-_tstart );
      _display_add( "INNER" );
      _display( os );
      continue;
     case ABORT: default:
      _open_nodes.erase( _open_nodes.begin() );
      delete pNode;
      _display_add( cpuclock()-_tstart );
      _display_add( "ABORT" );
      _display( os );
      continue;
     case UNDETERMINED:
     case FAILURE:
      _open_nodes.erase( _open_nodes.begin() );
      break;
    }
*/
    // Branch node domain
    std::pair<NODE*,NODE*> pNewNodes = _branch_node( pNode );
    if( !pNewNodes.first || !pNewNodes.second || _status == ABORT ){
      _display_add( cpuclock()-_tstart );
      _display_add( "ABORT" );
      _display( os );
      break; // terminate if branching unsuccessful
    }
    std::ostringstream obranch;
    obranch << "BRANCH:" << std::setw(std::ceil(std::log10(_np))) << _branch_var;
    delete pNode;

    // Assess left child node
    _status = pNewNodes.first->assess();
    switch( _status ){
     case OUTER:
      delete pNewNodes.first;
      obranch << " L:OUT";
      break;

     case INNER:
      _inner_nodes.insert( pNewNodes.first );
      //delete pNewNodes.first;
      obranch << " L:INN";
      break;

     case UNDETERMINED:
     case FAILURE:
      _open_nodes.insert( pNewNodes.first );
      if( pNewNodes.first->converged() ){ obranch << " L:CVG"; break; }     
      obranch << " L:BND (" << std::scientific << std::setprecision(3) << pNewNodes.first->measure() << ")";
      break;

     case ABORT: default:
      delete pNewNodes.first;
      delete pNewNodes.second;
      obranch << " L:ABORT";
      break;
    }

    // Check whether to continue
    if( _status == ABORT ){
      _display_add( cpuclock()-_tstart );
      _display_add( obranch.str() );
      _display( os );
      break;
    }

    // Assess right child node
    _status = pNewNodes.second->assess();
    switch( _status ){
     case OUTER:
      delete pNewNodes.second;
      obranch << " R:OUT";
      break;

     case INNER:
      _inner_nodes.insert( pNewNodes.second );
      //delete pNewNodes.second;
      obranch << " R:INN";
      break;

     case UNDETERMINED:
     case FAILURE:
      _open_nodes.insert( pNewNodes.second );
      if( pNewNodes.second->converged() ){ obranch << " R:CVG"; break; }     
      obranch << " R:BND (" << std::scientific << std::setprecision(3) << pNewNodes.second->measure() << ")";
      break;

     case ABORT: default:
      delete pNewNodes.first;
      delete pNewNodes.second;
      obranch << " R:ABORT";
      break;
    }
    _display_add( cpuclock()-_tstart );
    _display_add( obranch.str() );
    _display( os );

    // Check whether to continue
    if( _status == ABORT ) break;

    // Keep track of maximum number of open nodes
    if( _open_nodes.size() > _node_max ) _node_max = _open_nodes.size();
  }

  _display_final( cpuclock() );
  _display( os );
  
  return _open_measure;
}

template <typename T, typename NODE, typename LT_NODE>
inline bool
SBP<T,NODE,LT_NODE>::_terminate() const
{
  return( _open_measure < options.STOPPING_TOLERANCE? true : false );
}

template <typename T, typename NODE, typename LT_NODE>
inline void
SBP<T,NODE,LT_NODE>::_restart()
{
  _clean_stacks();
  _open_measure = 0.;
  _node_index = 1;
  _node_count = _node_max = 0;
  _tstart = cpuclock();
}

template <typename T, typename NODE, typename LT_NODE>
inline void
SBP<T,NODE,LT_NODE>::_clean_stacks()
{
  // Clean open nodes stack
  it_Nodes it = _open_nodes.begin();
  for( ; it != _open_nodes.end(); it++ ) delete *it;
  _open_nodes.clear();
  it = _inner_nodes.begin();
  for( ; it != _inner_nodes.end(); it++ ) delete *it;
  _inner_nodes.clear();
  _open_measure = 0.;
}

template <typename T, typename NODE, typename LT_NODE>
inline typename std::pair<NODE*,NODE*>
SBP<T,NODE,LT_NODE>::_branch_node
( NODE*pNode )
{
  // Branching set update
  _branching_variable_set( pNode );

  // Branching variable selection
  _branch_var = _select_branching_variable( pNode ).first;
  if( _status == ABORT ) return std::make_pair( (NODE*)0, (NODE*)0 );

  // Subnode creation
  return _branch_node( pNode, _branch_var, _node_count );
}

template <typename T, typename NODE, typename LT_NODE>
inline typename std::pair<NODE*,NODE*>
SBP<T,NODE,LT_NODE>::_branch_node
( NODE*pNode, const unsigned ivar, unsigned&node_count )
{
  // Partitionning
  std::pair<const T, const T> partition = _partition_variable_domain( pNode, ivar );
  std::vector<T> P_tmp = pNode->P();

  // Left child node
  P_tmp[ivar] = partition.first;
  NODE*pNodeL = new NODE( this, P_tmp, pNode->A(), ++node_count,
    pNode->depth()+1, _node_index, pNode->depend(), LEFT,
    std::make_pair(ivar, pNode->P(ivar)), pNode->converged(),
    pNode->data() );
  
  // Right child node
  P_tmp[ivar] = partition.second;
  NODE*pNodeR = new NODE( this, P_tmp, pNode->A(), ++node_count,
    pNode->depth()+1, _node_index, pNode->depend(), RIGHT, 
    std::make_pair(ivar, pNode->P(ivar)), pNode->converged(),
    pNode->data() );

  return std::make_pair( pNodeL, pNodeR );
}

template <typename T, typename NODE, typename LT_NODE>
inline std::pair<const T, const T>
SBP<T,NODE,LT_NODE>::_partition_variable_domain
( NODE*pNode, const unsigned ivar )
const
{
  // Branch at mid-point of current variable range
  std::pair<T,T> partition;
  partition.first = mc::Op<T>::l(pNode->P(ivar)) + Op<T>::zeroone()
    *Op<T>::diam(pNode->P(ivar))/2.;
  partition.second = mc::Op<T>::u(pNode->P(ivar)) - Op<T>::zeroone()
    *Op<T>::diam(pNode->P(ivar))/2.;
  return partition;
}

template <typename T, typename NODE, typename LT_NODE>
inline void
SBP<T,NODE,LT_NODE>::_branching_variable_set
( const NODE*pNode )
{
  _branch_set.clear();
  typename Options::SELECTION psel = options.BRANCHING_USERFUNCTION;
  std::set<unsigned> ssel;
  if( psel ) ssel = psel( pNode );
  for( unsigned ip=0; ip<_np; ip++ ){
    // Allow branching on var #ip if part of the user selection (if any)
    if( psel && ssel.find(ip) == ssel.end() ) continue;
    // Allow branching on var #ip if not excluded (_exclude_vars)
    // or not a dependent in current node (_pNode->depend())
    if( _exclude_vars.find(ip) != _exclude_vars.end()
     || pNode->depend().find(ip) != pNode->depend().end() )
      continue;
    // Add variable index to branching set
    _branch_set.insert( ip );
  }

  // Preselect variables based on scores
  _branching_score_subset( pNode );

  // Interrupt if branch set is empty - internal error...
  if( _branch_set.empty() ) throw Exceptions( Exceptions::BRANCH );
}

template <typename T, typename NODE, typename LT_NODE>
inline void
SBP<T,NODE,LT_NODE>::_branching_score_subset
( const NODE*pNode )
{
  if( !options.SCORE_BRANCHING_USE || pNode->scores().empty() ) return;

  // Create map of scores
  struct lt_scores{
    bool operator()
      ( const std::pair<unsigned,double>&el1, const std::pair<unsigned,double>el2 )
      { return el1.second < el2.second; }
  };
  std::multiset< std::pair<unsigned,double>, lt_scores > allscores;
  for( auto it = _branch_set.begin(); it != _branch_set.end(); ++it ){
    auto its = pNode->scores().find( *it );
    if( its == pNode->scores().end() ) continue;
    allscores.insert( *its );
#ifdef DEBUG__SBP_SCOREBRANCHING
    std::cout << "Score variable #" << its->first << ": " << its->second << std::endl;
    { int dum; std::cout << "PAUSED"; std::cin >> dum; }
#endif
  }

  // Keep best candidates based on map of scores
  for( auto it=allscores.begin(); it!=allscores.end(); ++it ){
    if( it->second >= allscores.rbegin()->second*(1.-options.SCORE_BRANCHING_RELTOL)
                     -options.SCORE_BRANCHING_ABSTOL
     && ( !options.SCORE_BRANCHING_MAXSIZE
       || _branch_set.size() < options.SCORE_BRANCHING_MAXSIZE ) ) break;
    _branch_set.erase( it->first );
  }
#ifdef MC__SBP_SHOW_SCOREBRANCHING
  std::cout << "***Branching variable score subset: {";
  for( auto it = _branch_set.begin(); it != _branch_set.end(); ++it )
    std::cout << " " << *it << " ";
  std:: cout << "}\n";
  { int dum; std::cout << "PAUSED"; std::cin >> dum; }
#endif
}

template <typename T, typename NODE, typename LT_NODE>
inline std::pair<unsigned, double>
SBP<T,NODE,LT_NODE>::_select_branching_variable
( NODE*pNode )
{
  // Strong branching strategy
  if( pNode->depth() < options.STRONG_BRANCHING_MAXDEPTH && !pNode->converged()
      && _branch_set.size() > 1 )
    return _strong_branching( pNode, _branch_var, _node_count );
  
  std::pair<unsigned, double> branchsel( _np, -1. ); // <- Can be any negative number
  switch( options.BRANCHING_VARIABLE_CRITERION ){

   // Branching based on relative range diameter
   case Options::RGREL:
#ifdef MC__SBP_SHOW_SCOREBRANCHING
    std::cout << "***Branching scores:";
#endif
    for( auto it = _branch_set.begin(); it != _branch_set.end(); ++it ){
      // Compare w.r.t. rotated root node bounds
      double score = Op<T>::diam( pNode->P(*it) )
                   / Op<T>::diam( _rotate( *it, _P_root, pNode->A() ) );
#ifdef MC__SBP_SHOW_SCOREBRANCHING
      std::cout << " " << *it << ": " << score;
#endif
      if( score > branchsel.second ) branchsel = std::make_pair( *it, score );
    }
#ifdef MC__SBP_SHOW_SCOREBRANCHING
    std::cout << std::endl;
#endif
    break;

   // Branching based on absolute range diameter
   case Options::RGABS:
    for( auto it = _branch_set.begin(); it != _branch_set.end(); ++it ){
      double score = Op<T>::diam( pNode->P(*it) );
      if( score > branchsel.second ) branchsel = std::make_pair( *it, score );
    }
    break;

   default:
    _status = ABORT;
    break;
  }

  return branchsel;
}

template <typename T, typename NODE, typename LT_NODE>
inline std::pair<unsigned, double>
SBP<T,NODE,LT_NODE>::_strong_branching
( NODE*pNode, unsigned&var_branch, unsigned&node_count )
{
  var_branch = _np;
  double min_branch = 0./0., max_branch = 0./0.; //NaN to initialize;
  std::vector<T> P_upd = pNode->P();
#ifdef DEBUG__SBP_STRONGBRANCHING
  std::cout << "*** SCORES";
#endif
  for( auto it = _branch_set.begin(); it != _branch_set.end(); ++it ){

    // Branch node domain
    std::pair<NODE*,NODE*> pNewNodes = _branch_node( pNode, *it, node_count );
    if( !pNewNodes.first || !pNewNodes.second || _status == ABORT )
      return std::make_pair( _np, 0./0. );
    pNewNodes.first->strongbranch()  = true; // In order for assessment methods to know
    pNewNodes.second->strongbranch() = true; // In order for assessment methods to know
    double meas = 0.;

    // Assess left child node
    _status = pNewNodes.first->assess();
    switch( _status ){
     case OUTER:
     case INNER:
     default:
      break;

     case UNDETERMINED:
     case FAILURE:
      meas += pNewNodes.first->measure();
      break;
    }

    // Check whether to continue
    if( _status == ABORT ){
      delete pNewNodes.first;
      delete pNewNodes.second;
      return std::make_pair( _np, 0./0. );
    }

    // Assess right child node
    _status = pNewNodes.second->assess();
    switch( _status ){
     case OUTER:
     case INNER:
     default:
      break;

     case UNDETERMINED:
     case FAILURE:
      meas += pNewNodes.second->measure();
      break;
    }

    // Check whether to continue
    if( _status == ABORT ){
      delete pNewNodes.first;
      delete pNewNodes.second;
      return std::make_pair( _np, 0./0. );
    }

#ifdef DEBUG__SBP_STRONGBRANCHING
    std::cout << "  " << *it << ": " << meas << " (" << min_branch << ")";
#endif

    // Update domain hull (excl. dependent variables)
    for( unsigned ip=0; ip<_np; ip++ ){
      if( pNewNodes.first->depend().find(ip)  != pNewNodes.first->depend().end()
       || pNewNodes.second->depend().find(ip) != pNewNodes.second->depend().end() )
        continue;
      T hull = Op<T>::hull( pNewNodes.first->P(ip), pNewNodes.second->P(ip) ), inter;
      if( !Op<T>::inter( inter, P_upd[ip], hull ) ) continue;
      P_upd[ip] = inter;
    }

    // Update most promising branching var
    if( !(meas >= min_branch) ){
      var_branch = *it;
      min_branch = meas;
    }
    if( !(meas <= max_branch) )
      max_branch = meas;

    // Cleanup
    delete pNewNodes.first;
    delete pNewNodes.second;

    // Interrupt if most promising branching less than absolute threshold
    if( !(min_branch > options.STRONG_BRANCHING_ABSTOL) ){
#ifdef DEBUG__SBP_STRONGBRANCHING
      std::cout << "\nStrong branching interrupted" << std::endl;
      { int dum; std::cin >> dum; }
#endif
      return std::make_pair( var_branch, min_branch );
    }
  }

  // Update node domain with new domain hull (excl. dependent variables)
  for( unsigned ip=0; ip<_np; ip++ ){
    if( pNode->depend().find(ip) != pNode->depend().end() ) continue;
    T inter;
    if( !Op<T>::inter( inter, pNode->P(ip), P_upd[ip] ) ) continue;
    pNode->P(ip) = inter;
  }

  // Has strong branching converged?
#ifdef DEBUG__SBP_STRONGBRANCHING
    std::cout << std::scientific << std::setprecision(4)
              << "\n  MIN: " << min_branch << "  MAX: " << max_branch
              << "  CUR: " << pNode->measure() << std::endl;
#endif
  if( max_branch < pNode->measure()
   && max_branch < min_branch*(1.+options.STRONG_BRANCHING_RELTOL)
                   +options.STRONG_BRANCHING_ABSTOL ){ //*pNode->measure() ){
#ifdef DEBUG__SBP_STRONGBRANCHING
    std::cout << "Strong branching has converged" << std::endl;
    { int dum; std::cin >> dum; }
#endif
    pNode->converged() = true;
  }

  return std::make_pair( var_branch, min_branch );
}

template <typename T, typename NODE, typename LT_NODE>
inline void
SBP<T,NODE,LT_NODE>::_display_open
( std::ostream&os ) const
{
  if( options.DISPLAY <= 2 ) return;
  cit_Nodes cit = _open_nodes.begin();
  for( unsigned id=1; cit != _open_nodes.end(); ++cit, id++ ){
    std::set<unsigned>::iterator ivar = _branch_set.begin();
    os << "Node " << id << ":" << std::scientific << std::setprecision(_DPREC)
       << std::setw(_DPREC+8) << (*cit)->measure();
    for( ; ivar != _branch_set.end(); ++ivar )
      os << (*cit)->P(*ivar) << "  ";
    os << std::endl;
  }
  return;
}

template <typename T, typename NODE, typename LT_NODE>
inline void
SBP<T,NODE,LT_NODE>::_display_init()
{
  _odisp.str("");
  if( options.DISPLAY <= 1 ) return;
  _odisp << std::right
  	 << std::setw(_IPREC) << "#   INDEX"
  	 << std::setw(_DPREC+8) << "MEASURE  "
  	 << std::setw(_IPREC) << "STACK"
  	 << std::setw(_IPREC) << "PARENT"
  	 << std::setw(_DPREC+8) << "CPU TIME "
  	 << std::setw(_DPREC+8) << "STATUS"
  	 << std::endl;  
}

template <typename T, typename NODE, typename LT_NODE>
inline void
SBP<T,NODE,LT_NODE>::_display_add
( const double dval )
{
  if( options.DISPLAY <= 1 ) return;
  _odisp << std::right << std::scientific << std::setprecision(_DPREC)
         << std::setw(_DPREC+8) << dval;
}

template <typename T, typename NODE, typename LT_NODE>
inline void
SBP<T,NODE,LT_NODE>::_display_add
( const unsigned ival )
{
  if( options.DISPLAY <= 1 ) return;
  _odisp << std::right << std::setw(_IPREC) << ival;
}

template <typename T, typename NODE, typename LT_NODE>
inline void
SBP<T,NODE,LT_NODE>::_display_add
( const std::string &sval )
{
  if( options.DISPLAY <= 1 ) return;
  _odisp << "  " << std::left << sval;
}

template <typename T, typename NODE, typename LT_NODE>
inline void
SBP<T,NODE,LT_NODE>::_display_final
( const double tcur )
{
  if( options.DISPLAY <= 0 ) return;

  // Solution found within allowed time?
  if( tcur-_tstart < options.MAX_CPUTIME
    && ( !options.MAX_NODES || _node_index <= options.MAX_NODES ) )
    _odisp << std::endl << "#  NORMAL TERMINATION:      ";
  else
    _odisp << std::endl << "#  EXECUTION STOPPED:       ";
  _odisp << std::fixed << std::setprecision(6) << tcur-_tstart << " CPU SEC"
         << std::endl;

  // Set-inversion results
  _odisp << "#  PARTITION MEASURE:    "
  	 << std::scientific << std::setprecision(_DPREC)
         << std::setw(_DPREC+6) << _open_measure
         << " (" << _open_nodes.size() << " NODES)"
  	 << std::endl;
  _odisp << "#  TOTAL NUMBER OF NODES:   " << _node_index-1 << std::endl
         << "#  MAX OPEN NODES IN STACK: " << _node_max << std::endl;
}

template <typename T, typename NODE, typename LT_NODE>
inline void
SBP<T,NODE,LT_NODE>::_display
( std::ostream &os )
{
  if( _odisp.str() == "" ) return;
  os << _odisp.str() << std::endl;
  _odisp.str("");
  return;
}

template <typename T, typename NODE, typename LT_NODE>
inline void
SBP<T,NODE,LT_NODE>::Options::display
( std::ostream&out ) const
{
  // Display SBP Options
  out << std::left;
  out << std::setw(60) << "  STOPPING TOLERANCE"
      << std::scientific << std::setprecision(1)
      << STOPPING_TOLERANCE << std::endl;
  out << std::setw(60) << "  PARTITION MEASURE CRITERION";
  switch( MEASURE ){
    case ABSMAXLEN:   out << "ABSMAXLEN\n";  break;
    case RELMAXLEN:   out << "RELMAXLEN\n";  break;
    case ABSMEANLEN:  out << "ABSMEANLEN\n"; break;
    case RELMEANLEN:  out << "RELMEANLEN\n"; break;
    case ABSVOLUME:   out << "ABSVOLUME\n";  break;
    case RELVOLUME:   out << "RELVOLUME\n";  break;
  }
  out << std::setw(60) << "  BRANCHING STRATEGY FOR VARIABLE SELECTION";
  switch( BRANCHING_VARIABLE_CRITERION ){
    case RGREL:  out << "RGREL\n";   break;
    case RGABS:  out << "RGABS\n";   break;
  }
  out << std::setw(60) << "  USE SCORE BRANCHING";
  switch( (int)SCORE_BRANCHING_USE ){
   case 0:  out << "N\n"; break;
   default: out << "Y\n";
            out << std::setw(60) << "  RELATIVE TOLERANCE FOR SCORE BRANCHING"
                << std::scientific << std::setprecision(1)
                << SCORE_BRANCHING_RELTOL << std::endl;
            out << std::setw(60) << "  ABSOLUTE TOLERANCE FOR SCORE BRANCHING"
                << std::scientific << std::setprecision(1)
                << SCORE_BRANCHING_ABSTOL << std::endl; break;
  }
  out << std::setw(60) << "  MAXIMUM DEPTH FOR STRONG BRANCHING";
  switch( STRONG_BRANCHING_MAXDEPTH ){
   case 0:  out << "-\n"; break;
   default: out << STRONG_BRANCHING_MAXDEPTH << std::endl;
            out << std::setw(60) << "  RELATIVE TOLERANCE FOR STRONG BRANCHING"
                << std::scientific << std::setprecision(1)
                << STRONG_BRANCHING_RELTOL << std::endl;
            out << std::setw(60) << "  ABSOLUTE TOLERANCE FOR STRONG BRANCHING"
                << std::scientific << std::setprecision(1)
                << STRONG_BRANCHING_ABSTOL << std::endl; break;
  }
  out << std::setw(60) << "  MAXIMUM ITERATION COUNT";
  switch( MAX_NODES ){
   case 0:  out << "-\n"; break;
   default: out << MAX_NODES << std::endl; break;
  }
  out << std::setw(60) << "  MAXIMUM CPU TIME (SEC)"
      << std::scientific << std::setprecision(1)
      << MAX_CPUTIME << std::endl;
  out << std::setw(60) << "  DISPLAY LEVEL"
      << DISPLAY << std::endl;
}

/////////////////////////////// SBPNode ///////////////////////////////

template <typename T>
inline
SBPNode<T>::SBPNode
( SBP< T,SBPNode<T>,lt_SBP<T> >*pSBP, const std::vector<T>&P,
  const unsigned index, const unsigned iter ):
_pSBP( pSBP ), _depth( 0 ), _index( index ), _iter( iter ),
_type( SBP< T,SBPNode<T>,lt_SBP<T> >::ROOT ),
_strongbranch( false ), _converged( false ), _data(0)
{
  _P = P;
  _A.resize(_pSBP->np(),_pSBP->np()).identity();
  _measure = _pSBP->measure( this );
  _status = SBP< T,SBPNode<T>,lt_SBP<T> >::UNDETERMINED;
}

template <typename T> template <typename U>
inline
SBPNode<T>::SBPNode
( SBP< T,SBPNode<T>,lt_SBP<T> >*pSBP, const std::vector<T>&P,
  const CPPL::dgematrix&A, const unsigned index, const unsigned depth,
  const unsigned iter, const std::set<unsigned>&depend,

  const typename SBP< T,SBPNode<T>,lt_SBP<T> >::NODETYPE type,
  const std::pair<unsigned,T>&parent, const bool converged, U*data ):
  _pSBP( pSBP ), _depth( depth ), _index( index ),
  _iter( iter ), _type( type ), _depend( depend ), _scores(),
  _parent( parent ), _strongbranch( false ), _converged( false ),
  _data( data ), _A( A ), _P( P )
{
  _measure = _pSBP->measure( this );
  _status = SBP< T,SBPNode<T>,lt_SBP<T> >::UNDETERMINED;
}

template <typename T>
inline
SBPNode<T>::~SBPNode()
{}

template <typename T>
inline typename SBP< T,SBPNode<T>,lt_SBP<T> >::STATUS
SBPNode<T>::assess()
{
  _status = _pSBP->assess( this );
  std::set<unsigned> varsel;
  if( _pSBP->options.BRANCHING_USERFUNCTION )
    varsel = _pSBP->options.BRANCHING_USERFUNCTION( this );
  _measure = _pSBP->measure( this );
  if( _measure < _pSBP->options.STOPPING_TOLERANCE )
    _converged = true;

  return _status;
}

template <typename T>
inline void
SBPNode<T>::output
( std::ostream& os_nodes, bool use ) const
{
  if( use )
    for( unsigned i=0; i<_pSBP->np(); ++i )
      os_nodes << "I("<<mc::Op<T>::l(_P[i]) <<","<< mc::Op<T>::u(_P[i]) << "), ";
  else
    for( unsigned i=0; i<_pSBP->np(); ++i )
      os_nodes << mc::Op<T>::l(_P[i]) <<"\t"<< mc::Op<T>::u(_P[i]) << "\t";
  os_nodes <<"Node #" << _index <<"\n";
}

} // namespace mc

#endif
