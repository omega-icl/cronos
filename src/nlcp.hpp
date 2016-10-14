// Copyright (C) 2015-2016 Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

/*!
\page page_NLCP Nonlinear Constraint Projection using complete search
\author Benoit C. Chachuat <tt>(b.chachuat@imperial.ac.uk)</tt> and OMEGA Research Group (http://www3.imperial.ac.uk/environmentenergyoptimisation)
\version 1.0
\date 2015
\bug No known bugs.

Consider a set of nonlinear constraints in the form:
\f{align*}
\mathcal{C}:\quad & g_j(x_1,\ldots,x_n)\ \leq,=,\geq\ 0,\ \ j=1,\ldots,m\\
& \quad x_i^L\leq x_i\leq x_i^U,\ \ i=1,\ldots,n\,,
\f}
where \f$g_1, \ldots, g_m\f$ are factorable, potentially nonlinear, real-valued functions; and \f$x_1, \ldots, x_n\f$ can be either continuous or integer decision variables. The class mc::NLCP computes an enclosure of all the solutions of \f$(\mathcal{C})\f$ using complete search. The main method implemented in mc::NLCP is set inversion, as first proposed by Walter and coworkers. The bounds and relaxations for the nonlinear or nonconvex constraints in \f$(\mathcal{C}\f$ are computed using various arithmetics in <A href="https://projects.coin-or.org/MCpp">MC++</A>.

\section sec_NLCP_setup How do I setup my constraint projection model?

Consider the system of nonlinear constraints:
\f{align*}
(1.-R) * (D/10/(1+b1)-x) * exp(10*x/(1+10*x/g)) - x & = 0\\
x - (1+b2)*y + (1.-R)*(D/10-b1*x-(1+b2)*y)*exp(10*y/(1+10*y/g)) & = 0\,,
\f}
with \f$(x,y) \in [0,1]^2\f$, \f$R \in [0.93,0.99]\f$, and \f$g=1000\f$, \f$D=22\f$, \f$b1=2\f$, \f$b2=2\f$.

First, we define an mc::NLCP class as below:

\code
  mc::NLCP CP;
\endcode

Next, we set the variables and objective/constraint functions by creating a direct acyclic graph (DAG) of the problem: 

\code
  #include "NLCP.hpp"

  mc::FFGraph DAG;
  const unsigned NP = 3, NC = 2;
  mc::FFVar P[NP], C[NC];
  for( unsigned int i=0; i<NP; i++ ) P[i].set( &DAG );

  double g=1e3, D=22., b1=2., b2=2.;
  mc::FFVar R = P[2];
  Y[0] = (1.-R)*(D/10/(1+b1)-P[0])*mc::exp(10*P[0]/(1+10*P[0]/g))-P[0];
  Y[1] = P[0]-(1+b2)*P[1]+(1.-R)*(D/10-b1*P[0]-(1+b2)*P[1])*mc::exp(10*P[1]/(1+10*P[1]/g));

  mc::NLCP<I> CP;
  CP.set_dag( &DAG );
  CP.set_par( NP, P );
  CP.set_dep( NY, Y );

  mc::FFGraph DAG;
  const unsigned NP = 4; mc::FFVar p[NP];
  for( unsigned i=0; i<NP; i++ ) p[i].set( &DAG );

  NLP.set_dag( &DAG );  // DAG
  NLP.set_var( NP, p ); // decision variables
  NLP.set_obj( mc::BASE_NLP::MIN, (p[0]*p[3])*(p[0]+p[1]+p[2])+p[2] ); // ojective
  NLP.add_ctr( mc::BASE_NLP::GE,  (p[0]*p[3])*p[1]*p[2]-25 );          // constraints
  NLP.add_ctr( mc::BASE_NLP::EQ,  sqr(p[0])+sqr(p[1])+sqr(p[2])+sqr(p[3])-40 );
  NLP.setup();
\endcode

The variable bounds and types are passed to mc::NLGO in invoking the various methods, as described below.


\section sec_NLGO_methods What are the methods available?

Other options can be modified to tailor the search, including output level, maximum number of iterations, tolerances, maximum CPU time, etc. These options can be modified through the public member mc::NLCP::options. 
*/

#ifndef MC__NLCP_H
#define MC__NLCP_H

#include <iostream>
#include <list>
#include "setinv.hpp"
#include "nlgo.hpp"
#include "aebnd.hpp"

#undef  MC__NLCP_DEBUG

namespace mc
{
//! @brief C++ class for nonlinear constraint projection using complete search
////////////////////////////////////////////////////////////////////////
//! mc::NLCP is a C++ class for constraint projection using
//! complete search. Relaxations for the nonlinear or nonconvex
//! participating terms are generated using MC++. Further details can be
//! found at: \ref page_NLCP
////////////////////////////////////////////////////////////////////////

template <typename T>
class NLCP:
  public virtual BASE_NLP,
  protected NLGO<T>,
  protected AEBND<T,CModel<T>,CVar<T>>,
  protected SetInv<CVar<T>>
{
  // Overloading stdout operator
  template <typename U> friend std::ostream& operator<<
    ( std::ostream&os, const NLCP<U>& );

  using NLGO<T>::_nvar;
  using NLGO<T>::_tvar;
  using NLGO<T>::_nctr;
  using NLGO<T>::_nrvar;
  using NLGO<T>::_nrdep;
  using NLGO<T>::_CMctr;
  //using NLGO<T>::_CMexcl;
  using NLGO<T>::_POLvar;
  using NLGO<T>::_set_chebscores;
  using SetInv<CVar<T>>::_open_nodes;
  using SetInv<CVar<T>>::_exclude_vars;
  using SetInv<CVar<T>>::_P_root;

 public:

  using NLGO<T>::setup;
  using NLGO<T>::stats;
  using SetInv<CVar<T>>::open_nodes;
  using SetInv<CVar<T>>::clusters;
  using SetInv<CVar<T>>::output_clusters;
  typedef SetInvNode<CVar<T>> NODE;

  //Constructor
  NLCP()
    : NLGO<T>(), AEBND<T,CModel<T>,CVar<T>>(), SetInv<CVar<T>>(), _CMrenv(0)
    {};

  // Default Destructor
  ~NLCP()
    { delete _CMrenv; }

  //! @brief NLCP options
  struct Options: public NLGO<T>::Options
  {
    //! @brief Constructor
    Options():
      NLGO<T>::Options(), BRANCHSEL(0), STGBCHRTOL(1e-2), STGBCHATOL(0e0),
      NODEMEAS(SetInv<CVar<T>>::Options::MEANWIDTH), CTRBACKOFF(0e0), CMREDORD(0),
      CMREDTHRES(1e-3), CMREDWARMS(false), CMREDALL(true),
      CMREDOPT(typename AEBND<T,CModel<T>,CVar<T> >::Options())
      //, INNREDMAX(10), INNREDTHRES(0.2)
      { CMREDOPT.DISPLAY = 0; }
    //! @brief Assignment operator
    Options& operator= ( Options&options ){
        NLGO<T>::Options::operator=( options );
        BRANCHSEL   = options.BRANCHSEL;
        STGBCHRTOL  = options.STGBCHRTOL;
        STGBCHATOL  = options.STGBCHATOL;
        NODEMEAS    = options.NODEMEAS;
        CTRBACKOFF  = options.CTRBACKOFF;
        CMREDORD    = options.CMREDORD;
        CMREDOPT    = options.CMREDOPT;
        CMREDTHRES  = options.CMREDTHRES;
        CMREDWARMS  = options.CMREDWARMS;
        CMREDALL    = options.CMREDALL;
        //INNREDMAX   = options.INNREDMAX;
        //INNREDTHRES = options.INNREDTHRES;
        return *this;
      }
    //! @brief Display options
    void display
      ( std::ostream&out ) const;
    //! @brief Branching-variable selection user-function
    typename SetInv<CVar<T>>::Options::SELECTION BRANCHSEL;
    //! @brief Relative tolerance for strong branching interruption
    double STGBCHRTOL;
    //! @brief Absolute tolerance for strong branching interruption
    double STGBCHATOL;
    //! @brief Partition measure
    int NODEMEAS;
    //! @brief Constraint back-off (avoids cutting-off part of the solution set after LP domain reduction)
    double CTRBACKOFF;
    //! @brief Reduced-space Chebyhev model order (0: no reduction)
    unsigned CMREDORD;
    //! @brief Reduced-space Chebyhev model threshold
    double CMREDTHRES;
    //! @brief Whether or not to warm start the reduction with linear estimators
    bool CMREDWARMS;
    //! @brief Whether or not to wait for all dependent variables to improve simultaneously before reduction
    bool CMREDALL;
    //! @brief Reduced-space Chebyhev model options (AEBND solver)
    typename AEBND<T,CModel<T>,CVar<T> >::Options CMREDOPT;
    //! @brief Maximum number of domain inner-reduction rounds
    //unsigned INNREDMAX;
    //! @brief Threshold for repeating inner reduction (minimum domain reduction ratio)
    //double INNREDTHRES;
  } options;

  //! @brief Solve constraint projection problem -- return value is a measure of boundary volume
  double solve
    ( const T*P, std::ostream&os=std::cout );

  //! @brief Solve relaxed optimization model using GUROBI within variable range <a>P</a> and for variable types <a>tvar</a> -- return value is SetInv<CVar<T>>::STATUS
  typename SetInv<CVar<T>>::STATUS assess
    ( const T*P, const unsigned*tvar=0, const unsigned refine=0,
      const bool reset=true );

  //! @brief Append all open (set-boundary) nodes to <a>os_open</a> - Number of significant digits is set via <a>DPREC</a> (default=6)
  void output_nodes
    ( std::ostream&os_open, const unsigned int DPREC=6 ) const;
  //! @brief Append all open (set-boundary) nodes to <a>os_open</a> - Number of significant digits is set via <a>DPREC</a> (default=6)
  void output_nodes
    ( std::ostream&os_open, const double NSAMP, const unsigned DPREC ) const;

 protected:
  //! @brief Set options of SetInv solver
  void _set_SetInvoptions
    () const;
  //! @brief User-function to subproblem assessment in SetInv
  typename SetInv<CVar<T>>::STATUS assess
    ( SetInvNode<CVar<T>>*node );
  //! @brief Subproblem for node assessment
  typename SetInv<CVar<T>>::STATUS _subproblem_assess
    ( T*P, std::map<unsigned,double>&scores );
  //! @brief Subproblem for node projection (reduced space)
  void _subproblem_project
    ( CVar<T>*CVP, std::set<unsigned>&depend, bool&converged );

  //! @brief Perform iterative bound contraction from relaxed model with inclusion test at each iteration
  typename SetInv<CVar<T>>::STATUS _assess
    ( T*P, unsigned&nred, const unsigned*tvar=0, const double*inc=0,
      const bool reset=true, const bool feastest=false );
  //! @brief Perform inclusion test based on constraint Chebysev models (if available)
  template <typename U> typename SetInv<CVar<T>>::STATUS _inclusion_test
    ( const std::vector<U>&C ) const;

  //! @brief Compute element volume
  double volume
    ( const std::vector<CVar<T>>&CVP ) const
    { double V=1.;
      for( unsigned i=_nrdep?_nrvar:0; i<CVP.size(); i++ )
        //V *= Op<T>::diam(CVP[i].R());
        V *= Op<T>::diam(CVP[i].R()) / Op<T>::diam(_P_root[i].B());
      return V; }      
  //! @brief Compute element mean width
  double meanwidth
    ( const std::vector<CVar<T>>&CVP ) const
    { return std::pow( volume(CVP), 1./(_nrdep?CVP.size()-_nrvar:CVP.size()) ); }
  //! @brief Compute element max width
  double maxwidth
    ( const std::vector<CVar<T>>&CVP ) const
    { double W=0.;
      for( unsigned i=_nrdep?_nrvar:0; i<CVP.size(); i++ )
        //if( W < Op<T>::diam(CVP[i].R()) ) W = Op<T>::diam(CVP[i].R());
        if( W*Op<T>::diam(_P_root[i].B()) < Op<T>::diam(CVP[i].R()) )
          W = Op<T>::diam(CVP[i].R()) / Op<T>::diam(_P_root[i].B());
      return W; }

 private:
  //! @brief Chebyshev reduced-space [-1,1] scaled model environment
  CModel<T>* _CMrenv;
  //! @brief Chebyshev reduced-space [-1,1] scaled variables
  std::vector< CVar<T> > _CMrbas;
  //! @brief Chebyshev reduced-space variables
  std::vector< CVar<T> > _CMrvar;
  //! @brief Chebyshev reduced-space dependents
  std::vector< CVar<T> > _CMrdep;

  //! @brief Private methods to block default compiler methods
  NLCP(const NLCP&);
  NLCP& operator=(const NLCP&);
};

template <typename T>
inline void
NLCP<T>::output_nodes
( std::ostream&os_open, const unsigned DPREC ) const
{
  // append all open sets from _open_nodes to os_open
  if( os_open.good() ){
    os_open << std::scientific << std::setprecision(DPREC); 
    auto cit = _open_nodes.begin();
    for( ; cit!=_open_nodes.end(); ++cit ){
      for( unsigned i=0; i<_nvar; i++ )
        os_open << std::setw(DPREC+9) << Op<CVar<T>>::l( (*cit)->P(i) )
                << std::setw(DPREC+9) << Op<CVar<T>>::u( (*cit)->P(i) );
      for( unsigned i=0; options.CMREDORD && i<_nrdep; i++ ){
        if( !(*cit)->P(_nrvar+i).env() ){
          os_open << std::setw(DPREC+9) << (*cit)->P(_nrvar+i).coefmon().second[0];
          for( unsigned j=1; j<_CMrenv->nmon(); j++ )
            os_open << std::setw(DPREC+9) << 0.;
        }
        else
          for( unsigned j=0; j<_CMrenv->nmon(); j++ )
            os_open << std::setw(DPREC+9) << (*cit)->P(_nrvar+i).coefmon().second[j];
        os_open << std::setw(DPREC+9) << Op<T>::l( (*cit)->P(_nrvar+i).R() )
                << std::setw(DPREC+9) << Op<T>::u( (*cit)->P(_nrvar+i).R() );
      }
      os_open << std::endl;
    }
  }
  //return SetInv<CVar<T>>::output_nodes( os_open, false, DPREC );
}

template <typename T>
inline void
NLCP<T>::output_nodes
( std::ostream&os_open, const double NSAMP, const unsigned DPREC ) const
{
  // Only for 1 reduced variables
  if( _nrvar != 1 || !options.CMREDORD || !_nrdep ) return;

  // append all open sets from _open_nodes to os_open
  if( os_open.good() ){
    os_open << std::scientific << std::setprecision(DPREC); 
    auto cit = _open_nodes.begin();
    for( ; cit!=_open_nodes.end(); ++cit ){
      for( unsigned k=0; k<NSAMP; k++ ){
        double rvar = Op<CVar<T>>::l( (*cit)->P(0) )
                    + Op<CVar<T>>::diam( (*cit)->P(0) ) * k / (NSAMP-1.);
        //double svar = (*cit)->type()<0? -1.+k/(NSAMP-1.): k/(NSAMP-1.);
        double svar = -1. + 2. * k / (NSAMP-1.);
        os_open << std::setw(DPREC+9) << rvar;
        for( unsigned i=0; i<_nrdep; i++ )
          os_open << std::setw(DPREC+9) << (*cit)->P(_nrvar+i).P(&svar)+Op<T>::l( (*cit)->P(_nrvar+i).R() )
                  << std::setw(DPREC+9) << (*cit)->P(_nrvar+i).P(&svar)+Op<T>::u( (*cit)->P(_nrvar+i).R() );
        os_open << std::endl;
      }
      os_open << std::endl;
    }
  }
}

template <typename T>
inline double
NLCP<T>::solve
( const T*P0, std::ostream&os )
{
  // Keep track of execution times
  stats.reset();

  // Set variable types and exclusions from B&B
  NLGO<T>::options = options;
  _set_SetInvoptions();
  std::vector<CVar<T>> CVP( P0, P0+_nvar );
  SetInv<CVar<T>>::variables( _nvar, CVP.data() );//, _CMexcl );

  // Set-up constraitn relaxations
  _tvar.clear();
  std::vector<T> P( P0, P0+_nvar );
  NLGO<T>::_set_polrelax( P.data(), 0, true );
#if defined (MC__NLCP_DEBUG)
  std::cout << *_dag;
  std::cout << NLGO<T>::_POLenv;
  //return 0;
#endif

  // Set-up reduced-space Chebyshev model (if applicable)
  _nrvar = BASE_NLP::_var.size();
  _nrdep = BASE_NLP::_dep.size();
  if( _nrdep && options.CMREDORD ){
    if( _CMrenv && (_CMrenv->nvar() != _nrvar || _CMrenv->nord() != options.CMREDORD) ){
      delete _CMrenv; _CMrenv = 0;
      _CMrbas.clear();
    }
    if( !_CMrenv ){
      // Create reduce-space Chebyshev model and populate with [-1,1] scaled variables
      _CMrenv = new CModel<T>( _nrvar, options.CMREDORD );
      for( unsigned i=0; i<_nrvar; i++ )
        _CMrbas.push_back( CVar<T>( _CMrenv, i, T(-1,1) ) );
    }
    _CMrvar.resize( _nrvar );
    _CMrdep.resize( _nrdep );
    AEBND<T,CModel<T>,CVar<T> >::options = options.CMREDOPT;
    AEBND<T,CModel<T>,CVar<T> >::setup();
  }

  // Run set-inversion algorithm
  double res = SetInv<CVar<T>>::solve( os );
  return res;
}

template <typename T>
inline typename SetInv<CVar<T>>::STATUS
NLCP<T>::assess
( SetInvNode<CVar<T>>*node )
{
#if defined (MC__NLCP_SHOW_BOXES)
  std::cout << "\nInitial Box:\n";
  for( unsigned i=0; i<_nvar; i++ )
    if( _exclude_vars.find(i) == _exclude_vars.end() )
      std::cout << node->P()[i] << std::endl;
#endif

  // Strong branching case
  unsigned DOMREDMAX = options.DOMREDMAX;
  if( node->strongbranch() && options.STGBCHDRMAX
   && options.STGBCHDRMAX < options.DOMREDMAX ){
    options.DOMREDMAX = options.STGBCHDRMAX; 
#if defined (MC__NLCP_DEBUG)
    std::cout << "on-going strong branching\n";
#endif
  }

  // Node assessment and domain contraction
  std::vector<T> IP;
  for( auto it=node->P().begin(); it!=node->P().end(); ++it )
    IP.push_back( (*it).B() );
  auto status = _subproblem_assess( IP.data(), node->scores() );
  node->P().assign( IP.begin(), IP.end() ); // Could be buggy!!!
#if defined (MC__NLCP_SHOW_BOXES)
  std::cout << "\nReduced Box:\n";
  for( unsigned i=0; i<_nvar; i++ )
    if( _exclude_vars.find(i) == _exclude_vars.end() )
      std::cout << node->P()[i] << std::endl;
#endif

  // Strong branching case
  if( node->strongbranch() ){
    options.DOMREDMAX = DOMREDMAX;
    return status;
  }

  // Subspace reduction
  if( status != SetInv<CVar<T>>::UNDETERMINED || !_nrdep || !options.CMREDORD )
    return status;
  _subproblem_project( node->P().data(), node->depend(), node->converged() );
#if defined (MC__NLCP_SHOW_BOXES)
  std::cout << "\nProjected Box:\n";
  for( unsigned i=0; i<_nvar; i++ )
    if( _exclude_vars.find(i) == _exclude_vars.end() )
      std::cout << node->P()[i] << std::endl;
#endif

  return status;
}

template <typename T>
inline typename SetInv<CVar<T>>::STATUS
NLCP<T>::_subproblem_assess
( T*P, std::map<unsigned,double>&scores )
{
  // Apply domain contraction
  // Do NOT test for infeasibility here, because contraction problem may
  // become infeasible due to round-off in LP solver
  unsigned *tvar = !_tvar.empty()? _tvar.data(): 0;
  unsigned nred = 0;
  auto status = _assess( P, nred, tvar, &options.CTRBACKOFF, false, true );
  if( status != SetInv<CVar<T>>::UNDETERMINED ) return status;

  // Solve relaxed feasibility problem
  //NLGO<T>::_update_polrelax( P, tvar, true );
  //std::cout << "Solving relaxed feasibility problem:\n"; {int dum; std::cin >> dum; }
  NLGO<T>::relax( P, tvar, 0, false, true );
  switch( NLGO<T>::get_status() ){
   // objective contains minimum backoffs
   case LPRELAX_BASE<T>::LP_OPTIMAL:  
    //std::cout << "Feasibility:" << NLGO<T>::get_objective() << std::endl;
    if( NLGO<T>::get_objective() > options.LPOPTIMTOL )
      return SetInv<CVar<T>>::OUTER;
    break; // need more test to find out if INNER or UNDETERMINED
   // infeasibility may not happen w/ backoffs
   case LPRELAX_BASE<T>::LP_INFEASIBLE:
   default:
     return SetInv<CVar<T>>::FAILURE;
  }

  // Compute scores
  _set_chebscores( scores, false );

  return SetInv<CVar<T>>::UNDETERMINED;
}

template <typename T> inline typename SetInv<CVar<T>>::STATUS
NLCP<T>::_assess
( T*P, unsigned&nred, const unsigned*tvar, const double*inc,
  const bool reset, const bool feastest )
{
  // Main loop for relaxation and domain reduction
  std::vector<T> P0;
  for( nred=0; nred < options.DOMREDMAX; nred++ ){
    // Update relaxation
    if( reset ) NLGO<T>::_set_polrelax( P, tvar, feastest );
    else        NLGO<T>::_update_polrelax( P, tvar, feastest );
    // Perform inclusion test
    //std::cout << "Solving inclusion test:\n"; {int dum; std::cin >> dum; }
    auto status = _inclusion_test( _CMctr );
    if( status != SetInv<CVar<T>>::UNDETERMINED ) return status;
    // Perform domain reduction
    P0.assign( P, P+_nvar );
    NLGO<T>::_contract( P, tvar, inc, false, feastest );
    if( NLGO<T>::get_status() != LPRELAX_BASE<T>::LP_OPTIMAL ) break;
    double vred = NLGO<T>::_reducrel( _nvar, P, P0.data() );
    if( vred < options.DOMREDTHRES ) break;
#ifdef MC__NLCP_DEBUG
    std::cout << "Reduction #" << nred << ": " << vred*1e2 << "%\n";
#endif
  }
  return SetInv<CVar<T>>::UNDETERMINED;
}

template <typename T>
template <typename U>
inline typename SetInv<CVar<T>>::STATUS 
NLCP<T>::_inclusion_test
( const std::vector<U>&C ) const
{
  if( C.empty() ) return SetInv<CVar<T>>::UNDETERMINED;

  // Function enclosure and inclusion tests
  bool flag1 = true, flag2 = true;
  auto itc = std::get<0>(NLGO<T>::_ctr).begin();
  for( unsigned j=0; flag2 && itc!=std::get<0>(NLGO<T>::_ctr).end(); ++itc, j++ ){
#if defined (MC__NLCP_SHOW_INCLUSIONS)
    std::cout << "C" << j << ":" << C[j];
    std::cout << "C" << j << " L:" << Op<U>::l(C[j]) << " U:" << Op<U>::u(C[j]) << std::endl;
#endif
    switch( (*itc) ){
      case EQ: flag1 = false;
               //if( Op<U>::u(C[j]) < 0. )
               //  flag2 = false;
               //else if( Op<U>::l(C[j]) < 0. )
               //  flag1 = false;
               // no break;
      case LE: if( Op<U>::l(C[j]) > 0. )
                 flag2 = false;
               else if( Op<U>::u(C[j]) > 0. )
                 flag1 = false;
               break;
      case GE: if( Op<U>::u(C[j]) < 0. )
                 flag2 = false;
               else if( Op<U>::l(C[j]) < 0. )
                 flag1 = false;
               break;
    }
  }
  if (flag1 && flag2)
    return SetInv<CVar<T>>::INNER;
  if (!flag2)
    return SetInv<CVar<T>>::OUTER;
  return SetInv<CVar<T>>::UNDETERMINED; 
}

template <typename T>
inline void
NLCP<T>::_subproblem_project
( CVar<T>*CVP, std::set<unsigned>&depend, bool&converged )
{
  // Apply implicit algebraic (if applicable)
  for( unsigned i=0; i<_nrvar; i++ ){
    _CMrvar[i] = _CMrbas[i] * (Op<CVar<T>>::diam(CVP[i])/2.) + Op<CVar<T>>::mid(CVP[i]);
#if defined (MC__NLCP_DEBUG)
    std::cout << "CMrvar[" << i << "] =" << _CMrvar[i];
#endif
  }
  try{
    if( !options.CMREDWARMS ){
      // Compute polynomial model on the dependents
      AEBND<T,CModel<T>,CVar<T> >::options.PMNOREM = false;
      auto status = AEBND<T,CModel<T>,CVar<T>>::solve( _CMrvar.data(), _CMrdep.data(), CVP+_nrvar );
      if( status != AEBND<T,CModel<T>,CVar<T>>::NORMAL ) return;
    }
    else{
      // Compute polynomial approximant
      AEBND<T,CModel<T>,CVar<T> >::options.PMNOREM = true;
      auto status0 = AEBND<T,CModel<T>,CVar<T>>::solve( _CMrvar.data(), _CMrdep.data(), CVP+_nrvar );
      if( status0 != AEBND<T,CModel<T>,CVar<T>>::NORMAL ) return;

      // Compute remainder bounds for the dependents
      std::vector<double> cobj( _nrvar+1 );
      std::vector<PolVar<T>> vobj( _nrvar+1 );

      for( unsigned i=0; i<_nrdep; i++ ){
#if defined (MC__NLCP_DEBUG)
        std::cout << "_CMrdep[" << i << "] =" << _CMrdep[i];
#endif
        // Set-up contaction objective function
        for( unsigned j=0; j<_nrvar; j++ ){
          cobj[j] = -_CMrdep[i].linear(j);
          vobj[j] = _POLvar[j];
        }
        cobj[_nrvar] = 1.;
        vobj[_nrvar] = _POLvar[_nrvar+i];

        // Solve for lower bound
        NLGO<T>::_set_LPrelax( _nrvar+1, vobj.data(), cobj.data(), MIN );
        NLGO<T>::_solve_LPmodel( options, stats, _var );
        double RxL = NLGO<T>::get_objective();
#if defined (MC__NLCP_DEBUG)
        std::cout << "  RxL = " << RxL << std::endl;
#endif

        // Solve for upper bound
        NLGO<T>::_set_LPrelax( _nrvar+1, vobj.data(), cobj.data(), MAX );
        NLGO<T>::_solve_LPmodel( options, stats, _var );
        double RxU = NLGO<T>::get_objective();
#if defined (MC__NLCP_DEBUG)
        std::cout << "  RxU = " << RxU << std::endl;
#endif

        // Setup remainder bounds in first-order polynomial model
#if defined (MC__NLCP_DEBUG)
        std::cout << "  diam = " << Op<CVar<T>>::diam(CVP[_nrvar+i]) << std::endl;
#endif
        if( !( RxU - RxL < 2.*Op<CVar<T>>::diam(CVP[_nrvar+i]) ) ) return;
        _CMrdep[i] = T( RxL, RxU );
        for( unsigned j=0; j<_nrvar; j++ )
          _CMrdep[i] += _CMrdep[i].linear(j) * _CMrvar[j];
      }
      AEBND<T,CModel<T>,CVar<T> >::options.PMNOREM = false;
      auto status = AEBND<T,CModel<T>,CVar<T>>::solve( _CMrvar.data(), _CMrdep.data(), _CMrdep.data() );
      if( status != AEBND<T,CModel<T>,CVar<T>>::NORMAL ) return;
    }
  }
  catch(...){ return;}
#ifdef MC__NLCP_DEBUG
  for( unsigned i=0; i<_nrdep; i++ ) std::cout << _CMrdep[i];
  {int dum; std::cout << "PAUSED"; std::cin >> dum;}
#endif

  // Update bounds and list of dependents if improvement
  depend.clear();
  if( options.CMREDALL ){
    bool allred = true;
    for( unsigned i=0; allred && i<_nrdep; i++ ){
#if defined (MC__NLCP_DEBUG)
      std::cout << "CMrdep[" << i << "] =" << _CMrdep[i];
#endif
      if( !Op<T>::lt( _CMrdep[i].R(), CVP[_nrvar+i].R() )
       || Op<T>::diam( _CMrdep[i].R() ) > options.CMREDTHRES ) allred = false;
    }
    if( !allred ) return;
#ifdef MC__NLCP_DEBUG
    std::cout << "converged node\n";
    { int dum; std::cout << "PAUSED"; std::cin >> dum; }
#endif
    for( unsigned i=0; i<_nrdep; i++ ){
      depend.insert( _nrvar+i );
      CVP[_nrvar+i] = _CMrdep[i];
    }
    converged = true;
  }
  else{
    for( unsigned i=0; i<_nrdep; i++ )
      if( Op<T>::lt( _CMrdep[i].R(), CVP[_nrvar+i].R() ) 
       && Op<T>::diam( _CMrdep[i].R() ) < options.CMREDTHRES ){
        depend.insert( _nrvar+i );
        CVP[_nrvar+i] = _CMrdep[i];
        converged = true;
      }
  }
}

template <typename T> inline void
NLCP<T>::_set_SetInvoptions
() const
{
  // SetInv options
  SetInv<CVar<T>>::options.MEASURE                      = options.NODEMEAS;
  SetInv<CVar<T>>::options.STOPPING_ABSTOL              = options.CVATOL;
  SetInv<CVar<T>>::options.STOPPING_RELTOL              = options.CVRTOL;
  SetInv<CVar<T>>::options.BRANCHING_VARIABLE_CRITERION = options.BRANCHVAR;
  SetInv<CVar<T>>::options.BRANCHING_USERFUNCTION       = options.BRANCHSEL;
  SetInv<CVar<T>>::options.SCORE_BRANCHING_USE          = options.SCOBCHUSE;
  SetInv<CVar<T>>::options.SCORE_BRANCHING_MAXSIZE      = options.SCOBCHVMAX;
  SetInv<CVar<T>>::options.SCORE_BRANCHING_RELTOL       = options.SCOBCHRTOL;
  SetInv<CVar<T>>::options.SCORE_BRANCHING_ABSTOL       = options.SCOBCHATOL;
  SetInv<CVar<T>>::options.STRONG_BRANCHING_MAXDEPTH    = options.STGBCHDEPTH;
  SetInv<CVar<T>>::options.STRONG_BRANCHING_RELTOL      = options.STGBCHRTOL;
  SetInv<CVar<T>>::options.STRONG_BRANCHING_ABSTOL      = options.STGBCHATOL;
  SetInv<CVar<T>>::options.MAX_NODES                    = options.MAXITER;
  SetInv<CVar<T>>::options.MAX_CPUTIME                  = options.MAXCPU;
  SetInv<CVar<T>>::options.DISPLAY                      = options.DISPLAY;
}

template <typename T>
inline void
NLCP<T>::Options::display
( std::ostream&out ) const
{
  // Display NLGO Options
  out << std::left;
  out << std::setw(60) << "  MAXIMUM OPTIMIZATION-BASED REDUCTION LOOPS"
      << NLGO<T>::Options::DOMREDMAX << std::endl;
  out << std::setw(60) << "  THRESHOLD FOR OPTIMIZATION-BASED REDUCTION LOOP"
      << std::fixed << std::setprecision(0)
      << NLGO<T>::Options::DOMREDTHRES*1e2 << "%\n";
  out << std::setw(60) << "  BACKOFF FOR OPTIMIZATION-BASED REDUCTION"
      << std::scientific << std::setprecision(1)
      << NLGO<T>::Options::DOMREDBKOFF << std::endl;
  out << std::setw(60) << "  CONSTRAINT BACKOFF";
  if( !(CTRBACKOFF > 0. ) )
    out << "-\n";
  else
    out << std::scientific << std::setprecision(1) << CTRBACKOFF << std::endl;
  out << std::setw(60) << "  POLYHEDRAL RELAXATION APPROACH";
  switch( NLGO<T>::Options::RELMETH ){
   case NLGO<T>::Options::DRL:    out << "DRL"    << std::endl; break;
   case NLGO<T>::Options::CHEB:   out << "CHEB"   << std::endl; break;
   case NLGO<T>::Options::HYBRID: out << "HYBRID" << std::endl; break;
  }
  out << std::setw(60) << "  ORDER OF CHEBYSHEV MODEL PROPAGATION";
  switch( NLGO<T>::Options::CMODPROP ){
   case 0:  out << "-\n"; break;
   default: out << NLGO<T>::Options::CMODPROP << std::endl; break;
  }
  out << std::setw(60) << "  ORDER OF CHEBYSHEV-DERIVED CUTS";
  if( !NLGO<T>::Options::CMODCUTS)
    switch( NLGO<T>::Options::CMODPROP ){
     case 0:  out << "-\n"; break;
     default: out << NLGO<T>::Options::CMODPROP << std::endl; break;
    }
  else 
    switch( NLGO<T>::Options::CMODCUTS ){
     case 0:  out << "-\n"; break;
     default: out << std::min(NLGO<T>::Options::CMODPROP,NLGO<T>::Options::CMODCUTS) << std::endl; break;
    }
  out << std::setw(60) << "  ORDER OF REDUCED-SPACE CHEBYSHEV MODEL";
  switch( CMREDORD ){
   case 0:  out << "-\n"; break;
   default: out << CMREDORD << std::endl;
            out << std::setw(60) << "  THRESHOLD FOR REDUCED-SPACE CHEBYSHEV MODEL"
                << std::scientific << std::setprecision(1)
                << CMREDTHRES << std::endl;
            out << std::setw(60) << "  SIMULTANEOUS REDUCED-SPACE CHEBYSHEV MODEL"
                << (CMREDALL?"Y\n":"N\n"); break;
  }
}

template <typename T>
inline std::ostream&
operator <<
( std::ostream&out, const NLCP<T>&CP )
{
  out << std::endl
      << std::setfill('_') << std::setw(72) << "\n\n" << std::setfill(' ')
      << std::setw(55) << "NONLINEAR CONSTRAINT PROJECTION IN CRONOS\n"
      << std::setfill('_') << std::setw(72) << "\n\n" << std::setfill(' ');

  // Display SetInv and NLCP Options
  CP._set_SetInvoptions();
  CP.SetInv<CVar<T>>::options.display( out );
  CP.NLCP<T>::options.display( out );

  out << std::setfill('_') << std::setw(72) << " " << std::setfill(' ')
      << std::endl << std::endl;
  return out;
}

/*
template <typename T>
inline typename SetInv<T>::STATUS 
NLCP<T>::_reduction_inner
( const std::vector<CUT>&CP, std::vector<T>&P, double&maxred ) const
{
  typename SetInv<T>::STATUS status = SetInv<T>::UNDETERMINED;
  maxred = 0.;
  const std::vector<T> P0 = P;
#if defined (MC__NLCP_DEBUG)
  std::cout << "\n Box before inner domain reduction:\n";
  for( unsigned i=0; i<_np; i++ ) cout << P[i] << endl;
#endif

  // Intersect linear inner relaxation constraints
  for( unsigned i=0; i<_np && status==SetInv<T>::UNDETERMINED; i++ ){
    typename std::vector<CUT>::const_iterator it = CP.begin();
    bool unbnd = false, PinLset = false, PinUset = false; 
    double PinL(0.), PinU(0.);
    for( unsigned k=0; it != CP.end(); ++it, k++ ){
      if( isequal( (*it).first(i), 0. ) ){ unbnd = true; break; }
      T Pbnd = (*it).second;
      for( unsigned j=0; j<_np; j++ )
        if( i != j ) Pbnd -= (*it).first(j)*P[j];
      Pbnd /= (*it).first(i);
      if( (*it).first(i)>0 ){
        PinU = PinUset? std::max( PinU, Op<T>::u(Pbnd) ): Op<T>::u(Pbnd), PinUset = true;
#if defined (MC__NLCP_DEBUG)
        std::cout << "Cut #" << k << "  PinU = " << PinU << std::endl;
#endif
      }
      else{
        PinL = PinLset? std::min( PinL, Op<T>::l(Pbnd) ): Op<T>::l(Pbnd), PinLset = true;
#if defined (MC__NLCP_DEBUG)
        std::cout << "Cut #" << k << "  PinL = " << PinL << std::endl;
#endif
      }
    }
    // Variable range is unbounded or outer parts intersect - cannot conclude
    if( unbnd || PinL < PinU )
      continue;
    // Outer parts on both sides w/o intersecting - no inner reduction possible
    else if( PinL < Op<T>::u(P[i]) && PinU > Op<T>::l(P[i]) )
      continue;
    // Outer parts on neither sides - no update
    else if( PinL >= Op<T>::u(P[i]) && PinU <= Op<T>::l(P[i]) )
      continue;
    // Outer part on the left-side only
    else if( PinL >= Op<T>::u(P[i]) && PinU < Op<T>::u(P[i]) )
      P[i] = T( Op<T>::l(P[i]), PinU );
    // Outer part on the right-side only
    else if( PinU <= Op<T>::l(P[i]) && PinL > Op<T>::l(P[i]) )
      P[i] = T( PinL, Op<T>::u(P[i]) );
    maxred = std::max( 1.-Op<T>::diam(P[i])/Op<T>::diam(P0[i]), maxred );
  }

#if defined (MC__NLCP_DEBUG)
  std::cout << "\n Box after inner domain reduction:\n";
  for( unsigned i=0; i<_np; i++ ) cout << P[i] << endl;
#endif
  return status;
}

template <typename T>
template < typename U >
inline std::vector<typename NLCP<T>::CUT>
NLCP<T>::_cuttingplanes_inner
( const std::vector<U>&Y, const std::vector<T>&P ) const
{
  std::vector<CUT> CP;
  CPPL::dcovector a(_np);

  // add the two hyperplanes that is generated from each time measurement
  typename std::list<Data>::const_iterator ity = _Ym->begin();
  for( unsigned i=0; ity != _Ym->end(); ++ity, i++ ){
    if( !_excl_output.empty()
      && _excl_output.find(i) != _excl_output.end() ) continue;
    U Yi = Y[i];
    for( unsigned k=0; k<_np; k++ )
      a(k) = Yi.linear( k, true );

    double bU = (*ity).yl - Op<U>::l( Yi );
    for(unsigned k=0;k<_np;k++) bU += a(k) * Op<T>::mid( P[k] );
    CP.push_back( std::make_pair( a, bU ) );

    double bL = -(*ity).yu + Op<U>::u( Yi );  
    for(unsigned k=0;k<_np;k++) bL -= a(k) * Op<T>::mid( P[k] );
    CP.push_back( std::make_pair( -a, bL ) );
  }

#if defined (MC__NLCP_DEBUG)
  std::cout << " Inner Cutting Planes a^T.x <= b: " << std::endl;
  std::vector<CUT>::const_iterator it = CP.begin();
  for( ; it!=CP.end(); ++it ){
    std::cout << (*it).first;  
    std::cout <<" <=  " << (*it).second << std::endl;
  }
#endif

  return CP;
}
*/
} // namespace mc
#endif

