// Copyright (C) 2015-2018 Benoit Chachuat, Imperial College London.
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
#include "nlgo.hpp"

#undef  MC__NLCP_DEBUG

namespace mc
{
//! @brief C++ class for nonlinear constraint projection using complete search
////////////////////////////////////////////////////////////////////////
//! mc::NLCP is a C++ class for constraint projection using complete
//! search. Relaxations for nonlinear / nonconvex participating
//! functions are generated using MC++. Further details can be
//! found at: \ref page_NLCP
////////////////////////////////////////////////////////////////////////

template <typename T>
class NLCP:
  //public virtual CSEARCH_BASE<T>,
  public virtual BASE_NLP,
  protected NLGO<T>
{
public:

  //! @brief NLCP options
  struct Options: public NLGO<T>::Options
  {
    //! @brief Constructor
    Options():
      NLGO<T>::Options(), CVTOL(1e-3), BRANCHSEL(0), STGBCHRTOL(1e-2), STGBCHATOL(0e0),
      NODEMEAS(SBP<T>::Options::RELMEANLEN), CTRBACKOFF(0e0), VARMEAS(ALL)
      //, INNREDMAX(10), INNREDTHRES(0.2)
      {}
    //! @brief Assignment operator
    Options& operator= ( Options&options ){
        NLGO<T>::Options::operator=( options );
        CVTOL       = options.CVTOL;
        BRANCHSEL   = options.BRANCHSEL;
        STGBCHRTOL  = options.STGBCHRTOL;
        STGBCHATOL  = options.STGBCHATOL;
        NODEMEAS    = options.NODEMEAS;
        CTRBACKOFF  = options.CTRBACKOFF;
        VARMEAS     = options.VARMEAS;
        //INNREDMAX   = options.INNREDMAX;
        //INNREDTHRES = options.INNREDTHRES;
        return *this;
      }
    //! @brief Node measure strategy
    enum OPTMEAS{
      ALL=0,	//!< Measure all variables
      DEP,   	//!< Measure dependent variables only
      INDEP,   	//!< Measure independent variables only
      USER   	//!< Measure user-defined branching set
    };
    //! @brief Display options
    void display
      ( std::ostream&out ) const;
    //! @brief Convergence tolerance
    double CVTOL;
    //! @brief Branching-variable selection user-function
    typename SBP<T>::Options::SELECTION BRANCHSEL;
    //! @brief Relative tolerance for strong branching interruption
    double STGBCHRTOL;
    //! @brief Absolute tolerance for strong branching interruption
    double STGBCHATOL;
    //! @brief Partition measure strategy
    int NODEMEAS;
    //! @brief Constraint back-off (avoids cutting-off part of the solution set after LP domain reduction)
    double CTRBACKOFF;
    //! @brief Variable measure strategy
    OPTMEAS VARMEAS;
    //! @brief Maximum number of domain inner-reduction rounds
    //unsigned INNREDMAX;
    //! @brief Threshold for repeating inner reduction (minimum domain reduction ratio)
    //double INNREDTHRES;
  } options;

  // Overloading stdout operator
  template <typename U> friend std::ostream& operator<<
    ( std::ostream&os, const NLCP<U>& );

  //Constructor
  NLCP()
    : NLGO<T>()
    {};

  // Default Destructor
  ~NLCP()
    {}

  //! @brief Setup DAG for cost and constraint evaluation
  void setup
    ( const T*P, const unsigned*tvar=0, std::ostream&os=std::cout );

  //! @brief Solve bound contraction problems from relaxed constraints, starting with variable range <a>P</a>, for variable types <a>tvar</a> and using the options specified in <a>NLGO::Options::DOMREDMAX</a> and <a>NLGO::Options::DOMREDTHRES</a> -- returns updated variable bounds <a>P</a>, and number of iterative refinements <a>nred</a> 
  typename LPRELAX_BASE<T>::LP_STATUS contract
    ( T*P, unsigned&nred, const unsigned*tvar=0,
      const bool reset=true, const bool feastest=false );

  //! @brief Run constraint projection algorithm -- return value is a measure of boundary volume
  double solve
    ( const T*P, const unsigned*tvar=0, std::ostream&os=std::cout );

  //! @brief Resume constraint projection algorithm -- return value is a measure of boundary volume
  double solve
    ( std::ostream&os=std::cout );

  //! @brief Public members of SBP
  typedef SBPNode<T> NODE;
  using SBP<T>::open_nodes;
  using SBP<T>::clusters;
  using SBP<T>::output_clusters;
  using SBP<T>::output_nodes;
  using NLGO<T>::stats;

protected:

  //! @brief Protected members of CSEARCH_BASE
  using CSEARCH_BASE<T>::_nvar;
  using CSEARCH_BASE<T>::_nrvar;
  using CSEARCH_BASE<T>::_nrdep;
  using CSEARCH_BASE<T>::_var_fperm;
  using CSEARCH_BASE<T>::_CMndxdep;
  using CSEARCH_BASE<T>::_CMrdep;
  using NLGO<T>::_var_excl;
  using NLGO<T>::_Ivar;
  using NLGO<T>::_tvar;

  //! @brief Set local optimizer
  virtual void _set_SLVLOC
    ()
    {}

  //! @brief User-function to subproblems in SBB
  virtual typename SBP<T>::STATUS assess
    ( SBPNode<T>*node )
    { return CSEARCH_BASE<T>::_assess( options, stats, node ); }

  //! @brief Compute element volume
  virtual double volume
    ( const NODE*pNode, unsigned &dim, const bool rel=true ) const;

  //! @brief Compute element mean width (equivalent cube edge)
  virtual double meanwidth
    ( const NODE*pNode, const bool rel=true ) const;

  //! @brief Compute element max width
  virtual double maxwidth
    ( const NODE*pNode, const bool rel=true ) const;

 private:

  //! @brief Private methods to block default compiler methods
  NLCP(const NLCP&);
  NLCP& operator=(const NLCP&);
};

template <typename T>
inline void
NLCP<T>::setup
( const T*P, const unsigned*tvar, std::ostream&os )
{
  NLGO<T>::options = options;
  return NLGO<T>::setup( P, tvar, os );
}

template <typename T>
inline typename LPRELAX_BASE<T>::LP_STATUS 
NLCP<T>::contract
( T*P, unsigned&nred, const unsigned*tvar,
  const bool reset, const bool feastest )
{
  NLGO<T>::options = options;
  return NLGO<T>::contract( P, nred, tvar, 0, reset, feastest );
}

template <typename T>
inline double
NLCP<T>::solve
( const T*P, const unsigned*tvar, std::ostream&os )
{
  if( !NLGO<T>::_issetup || !NLGO<T>::_init( P, tvar ) )
    throw typename NLGO<T>::Exceptions( NLGO<T>::Exceptions::SETUP );
  return CSEARCH_BASE<T>::_solve_sbp( options, stats, _Ivar.data(), tvar?_tvar.data():0, os );
}

template <typename T>
inline double
NLCP<T>::solve
( std::ostream&os )
{
  if( !NLGO<T>::_issetup )
    throw typename NLGO<T>::Exceptions( NLGO<T>::Exceptions::SETUP );
  return CSEARCH_BASE<T>::_solve_sbp( options, stats, 0, 0, os );
}

template <typename T>
inline double
NLCP<T>::volume
( const NODE*pNode, unsigned &dim, const bool rel ) const
{
  dim = 0;
  double V=1.;
  auto sset = options.VARMEAS == Options::USER && options.BRANCHSEL?
              options.BRANCHSEL( pNode ): std::set<unsigned>();
  for( unsigned i=0; i<_nvar; i++ ){
    if( Op<T>::diam(pNode->P0(i)) < options.DOMREDMIG )
      continue; // Do not account for variables whose domain is too narrow
    switch( options.VARMEAS ){
    case Options::ALL: default:
      break;
    case Options::INDEP:
      if( _nrvar && _var_fperm[i] >= _nrvar ) continue;
      break;
    case Options::DEP:
      if( _nrdep && _var_fperm[i] < _nrvar ) continue;
      break;
    case Options::USER:
      if( !sset.empty() && sset.find(i) == sset.end() ) continue;
      break;
    }
    ++dim;
    double scal = rel? Op<T>::diam(pNode->P0(i)): 1;
    if( _var_fperm[i] < _nrvar || _CMndxdep.find( i ) == _CMndxdep.end() )
      V *= Op<T>::diam(pNode->P(i)) / scal;
    else
      V *= Op<T>::diam(_CMrdep[_var_fperm[i]-_nrvar].R()) / scal;
  }
  return V;
}

template <typename T>
inline double
NLCP<T>::meanwidth
( const NODE*pNode, const bool rel ) const
{
  unsigned dim(0);
  double V = volume( pNode, dim, rel );
  return std::pow( V, 1./dim );
}

template <typename T>
inline double
NLCP<T>::maxwidth
( const NODE*pNode, const bool rel ) const
{
  double W=0.;
  auto sset = options.VARMEAS == Options::USER && options.BRANCHSEL?
              options.BRANCHSEL( pNode ): std::set<unsigned>();
  for( unsigned i=0; i<_nvar; i++ ){
    if( Op<T>::diam(pNode->P0(i)) < options.DOMREDMIG )
      continue; // Do not account for variables whose domain is too narrow
    switch( options.VARMEAS ){
    case Options::ALL: default:
      break;
    case Options::INDEP:
      if( _nrvar && _var_fperm[i] >= _nrvar ) continue;
      break;
    case Options::DEP:
      if( _nrdep && _var_fperm[i] < _nrvar ) continue;
      break;
    case Options::USER:
      if( !sset.empty() && sset.find(i) == sset.end() ) continue;
      break;
    }
    double scal = rel? Op<T>::diam(pNode->P0(i)): 1;
    if( _var_fperm[i] < _nrvar || _CMndxdep.find( i ) == _CMndxdep.end() ){
      if( W < Op<T>::diam(pNode->P(i)) / scal )
        W = Op<T>::diam(pNode->P(i)) / scal;
    }
    else{
      if( W < Op<T>::diam(_CMrdep[_var_fperm[i]-_nrvar].R()) / scal )
        W = Op<T>::diam(_CMrdep[_var_fperm[i]-_nrvar].R()) / scal;
    }
  }
  return W;
}

template <typename T>
inline void
NLCP<T>::Options::display
( std::ostream&out ) const
{
  // Display NLCP Options
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
   case NLGO<T>::Options::ISM:    out << "ISM"    << std::endl; break;
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
  switch( NLGO<T>::Options::CMODDEPS ){
   case 0:  out << "-\n"; break;
   default: out << NLGO<T>::Options::CMODDEPS << std::endl;
            out << std::setw(60) << "  CONVERGENCE THRESHOLD FOR REDUCED-SPACE CHEBYSHEV MODEL"
                << std::scientific << std::setprecision(1)
                << NLGO<T>::Options::CMODRTOL << std::endl;
            out << std::setw(60) << "  USE OF CHEBYSHEV-REDUCTION CONSTRAINTS";
            switch( NLGO<T>::Options::CMODRED ){
             case NLGO<T>::Options::NOREDUC: out << "NOREDUC" << std::endl; break;
             case NLGO<T>::Options::APPEND:  out << "APPEND"  << std::endl; break;
            }
//            out << std::setw(60) << "  SIMULTANEOUS REDUCED-SPACE CHEBYSHEV MODEL"
//                << (REDALL?"Y\n":"N\n"); break;
  }
  out << std::setw(60) << "  APPEND NCO CUTS"
      << (NLGO<T>::Options::NCOCUTS?"Y\n":"N\n");
  if( NLGO<T>::Options::CMODCUTS ){
    out << std::setw(60) << "  METHOD FOR NCO CUTS";
    switch( NLGO<T>::Options::NCOMETH ){
     case NLGO<T>::Options::FSA: out << "FORWARD SENS.\n";
     case NLGO<T>::Options::ASA: out << "ADJOINT SENS.\n";
    }
  }
}

template <typename T>
inline std::ostream&
operator <<
( std::ostream&out, const NLCP<T>&CP )
{
  out << std::right << std::endl
      << std::setfill('_') << std::setw(72) << " " << std::endl << std::endl << std::setfill(' ')
      << std::setw(55) << "NONLINEAR CONSTRAINT PROJECTION IN CRONOS\n"
      << std::setfill('_') << std::setw(72) << " " << std::endl << std::endl << std::setfill(' ');

  // Display SBP and NLCP Options
  CP.CSEARCH_BASE<T>::_set_SBPoptions( CP.NLCP<T>::options );
  CP.SBP<T>::options.display( out );
  CP.NLCP<T>::options.display( out );

  out << std::setfill('_') << std::setw(72) << " " << std::endl << std::endl << std::setfill(' ');
  return out;
}
/*
template <typename T>
inline typename SBP<T>::STATUS 
NLCP<T>::_reduction_inner
( const std::vector<CUT>&CP, std::vector<T>&P, double&maxred ) const
{
  typename SBP<T>::STATUS status = SBP<T>::UNDETERMINED;
  maxred = 0.;
  const std::vector<T> P0 = P;
#if defined (MC__NLCP_DEBUG)
  std::cout << "\n Box before inner domain reduction:\n";
  for( unsigned i=0; i<_np; i++ ) cout << P[i] << endl;
#endif

  // Intersect linear inner relaxation constraints
  for( unsigned i=0; i<_np && status==SBP<T>::UNDETERMINED; i++ ){
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
  //return SBP<T>::output_nodes( os_open, false, DPREC );
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
*/
} // namespace mc
#endif


