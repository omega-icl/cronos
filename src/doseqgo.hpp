// Copyright (C) 2017 Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

/*!
\page page_DOSEQGO Global Dynamic Optimization using MC++
\author Benoit C. Chachuat <tt>(b.chachuat@imperial.ac.uk)</tt> and OMEGA Research Group (http://www3.imperial.ac.uk/environmentenergyoptimisation)
\version 1.0
\date 2017
\bug No known bugs.

Consider a dynamic optimization (DO) problem of the form
\f{align*}
\mathcal{P}:\quad & \min_{{\bf p}\in P}\ & \phi_0({\bf p}, {\bf x}(t_0),\ldots,{\bf x}(t_N)) + \int_{t_0}^{t_N} \psi_0({\bf p}, {\bf x}(t))\, dt\\
  & {\rm s.t.} \ \ \phi_k({\bf p}, {\bf x}(t_0),\ldots,{\bf x}(t_N)) + \int_{t_0}^{t_N} \psi_k({\bf p}, {\bf x}(t))\, dt \{\geq,=,\leq\} 0, \quad k=1,\ldots,n_c\\
  & \qquad \dot{{\bf x}}(t)={\bf f}({\bf p},{\bf x}(t)),\ t\in(t_0,t_{N}]; \quad {\bf x}(t_0)={\bf h}({\bf p})
\f}
where \f${\bf f}\f$, \f${\bf h}\f$, \f$\phi_k\f$ and \f$\psi_k\f$ are factorable real-valued functions and twice continuously differentiable in all their arguments; \f$p_1, \ldots, p_n\f$ are either continuous or integer decision variables in \f$P\subset\mathbb{R}^{n_p}\f$; and \f${\bf x}(t)\in\mathbb{R}^{n_x}\f$ are the parametric state trajectories,

. The class mc::DOSEQGO solves such DO or MIDO problems to global optimality using complete search (branch-and-bound). Bounds on the nonconvex participating terms and the nonlinear ODE solutions are generated using various arithmetics in <A href="https://projects.coin-or.org/MCpp">MC++</A> and the class mc::ODEBND (\ref page_ODEBND). 

\section sec_DOSEQGO_solve How to Solve a Dynamic Optimization Model using mc::DOSEQSLV_IPOPT?

Consider the following dynamic optimization problem:
\f{align*}
  \max_{t_{\rm f},x_1^0,u_1,\ldots,u_{N_{\rm s}}}\ & t_{\rm f} \\
  \text{s.t.} \ & x_1(1) = 1\\
  & x_2(1) = 0\\
  & \left.\begin{array}{l}
    \dot{x_1}(\tau) = x_2(\tau)\\
    \dot{x_2}(\tau) = (u_k-x_1(\tau)-2*x_2(\tau))
    \end{array}\right\},\ \tau\in({\textstyle\frac{k-1}{N},\frac{k}{N}}],\ k=1,\ldots,N_{\rm s}\\
  & x_1(0) = x_1^0,\ x_2(0) = 0
\f}
with \f$N_{\rm s}=10\f$ stages.

First, we define an mc::DOSEQGO class as below:

\code
  mc::DOSEQGO NLP;
\endcode

Next, we set the variables and objective/constraint functions by creating a direct acyclic graph (DAG) of the problem: 

\code
  #include "doseqgo.hpp"
  
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

The variable bounds and types are passed to mc::DOSEQGO in invoking the various methods, as described below.


\section sec_DOSEQGO_methods What are the methods available?


Given initial bounds \f$P\f$ and initial guesses \f$p_0\f$ on the decision variables, the NLP model is solved using branch-and-bound search (default) as follows:

\code
  #include "interval.hpp"

  typedef mc::Interval I;
  I Ip[NP] = { I(1,5), I(1,5), I(1,5), I(1,5) };

  std::cout << NLP;
  int status = NLP.solve( Ip );
\endcode

The following result is produced:

\verbatim
______________________________________________________________________

               NONLINEAR GLOBAL OPTIMIZATION IN CRONOS
______________________________________________________________________

  COMPLETE SEARCH METHOD                                    SBB
  CONVERGENCE ABSOLUTE TOLERANCE                            1.0e-03
  CONVERGENCE RELATIVE TOLERANCE                            1.0e-03
  ROOT NODE PROPROCESSING                                   Y
  OPTIMIZATION-BASED REDUCTION MAX LOOPS                    10
  OPTIMIZATION-BASED REDUCTION THRESHOLD LOOP               20%
  MAXIMUM ITERATION COUNT                                   -
  MAXIMUM CPU TIME (SEC)                                    7.2e+03
  DISPLAY LEVEL                                             2
 _______________________________________________________________________

 INDEX STACK    CUMUL TIME        RELAX         INC   PARENT        LBD           UBD         ACTION  

     1     1  0.000000e+00 -1.000000e+20  1.701402e+01     0  1.511037e+01  1.701402e+01       BRANCH2
     2     2  2.446400e-02  1.511037e+01  1.701402e+01     1  1.659980e+01  1.701402e+01       BRANCH3
     3     3  3.611400e-02  1.511037e+01  1.701402e+01     1  1.614164e+01  1.701402e+01       BRANCH2
     4     4  5.255800e-02  1.614164e+01  1.701402e+01     3  1.669923e+01  1.701402e+01       BRANCH3
     5     5  6.833200e-02  1.614164e+01  1.701402e+01     3  1.713529e+01       SKIPPED        FATHOM
     6     4  8.233300e-02  1.659980e+01  1.701402e+01     2  1.664991e+01  1.701402e+01       BRANCH2
     7     5  9.279500e-02  1.659980e+01  1.701402e+01     2  1.701402e+01       SKIPPED        FATHOM
     8     4  1.237580e-01  1.664991e+01  1.701402e+01     6  1.703133e+01       SKIPPED        FATHOM
     9     3  1.286230e-01  1.664991e+01  1.701402e+01     6  1.701402e+01       SKIPPED        FATHOM
    10     2  1.553660e-01  1.669923e+01  1.701402e+01     4  1.672776e+01  1.701402e+01       BRANCH3
    11     3  1.638520e-01  1.669923e+01  1.701402e+01     4  1.701406e+01       SKIPPED        FATHOM
    12     2  1.661390e-01  1.672776e+01  1.701402e+01    10  1.676720e+01  1.704213e+01       BRANCH2
    13     3  1.748900e-01  1.672776e+01  1.701402e+01    10  1.701402e+01       SKIPPED        FATHOM
    14     2  2.044000e-01  1.676720e+01  1.701402e+01    12    INFEASIBLE       SKIPPED        FATHOM
    15     1  2.149610e-01  1.676720e+01  1.701402e+01    12  1.703627e+01       SKIPPED        FATHOM

#  NORMAL TERMINATION: 0.222576 CPU SEC  (LBD:80.6%  UBD:19.3%)
#  INCUMBENT VALUE:  1.701402e+01
#  INCUMBENT POINT:  1.000000e+00  4.742955e+00  3.821208e+00  1.379400e+00
#  INCUMBENT FOUND AT NODE: 10
#  TOTAL NUMBER OF NODES:   15
#  MAXIMUM NODES IN STACK:  5
\endverbatim

Integrality of certain variables can be enforced by passing an additional argument as follows, here enforcing integrality of the second and third variables:

\code
  unsigned Tp[NP] = { 0, 1, 1, 0 };

  int status = NLP.solve( Ip, Tp );
\endcode

The following result is now produced:

\verbatim
 INDEX STACK    CUMUL TIME        RELAX         INC   PARENT        LBD           UBD         ACTION  

     1     1  0.000000e+00 -1.000000e+20  1.000000e+20     0  1.224000e+01  2.312461e+01       BRANCH0
     2     2  2.387200e-02  1.224000e+01  2.312461e+01     1  1.224000e+01  2.312461e+01       BRANCH1
     3     3  3.766500e-02  1.224000e+01  2.312461e+01     1  3.200000e+01       SKIPPED        FATHOM
     4     2  4.145500e-02  1.224000e+01  2.312461e+01     2  2.501911e+01       SKIPPED        FATHOM
     5     1  7.219900e-02  1.224000e+01  2.312461e+01     2  1.488304e+01  2.312461e+01       BRANCH2
     6     2  1.082510e-01  1.488304e+01  2.312461e+01     5  2.312460e+01       SKIPPED        FATHOM
     7     1  1.395070e-01  1.488304e+01  2.312461e+01     5  2.012448e+01  2.312461e+01       BRANCH1
     8     2  1.729230e-01  2.012448e+01  2.312461e+01     7  2.511437e+01       SKIPPED        FATHOM
     9     1  1.866310e-01  2.012448e+01  2.312461e+01     7  2.312458e+01       SKIPPED        FATHOM

#  NORMAL TERMINATION: 0.215300 CPU SEC  (LBD:94.3%  UBD:5.7%)
#  INCUMBENT VALUE:  2.312461e+01
#  INCUMBENT POINT:  1.000000e+00  5.000000e+00  3.000000e+00  2.236068e+00
#  INCUMBENT FOUND AT NODE: 7
#  TOTAL NUMBER OF NODES:   9
#  MAXIMUM NODES IN STACK:  3
\endverbatim

The complete search algorithm based on piecewise linearization can be selected by modifying the options as follows:

\code
  NLP.options.CSALGO = mc::DOSEQGO<I>::Options::PWL;
  NLP.options.POLIMG.BREAKPOINT_TYPE = mc::PolImg<I>::Options::BIN;
  int status = NLP.solve( Ip, Tp );
\endcode

The following result is now produced:

\verbatim
     ITER  PREPROC      RELAX          BEST        CPU TOT          STATUS
        1        -  1.224000e+01  2.312461e+01  5.066300e-02        REFINE
        2        -  1.410526e+01  2.312461e+01  1.267340e-01        REFINE
        3       R3  1.425993e+01  2.312461e+01  2.750010e-01        REFINE
        4       R1  1.488303e+01  2.312461e+01  3.852940e-01        REFINE
        5       R3  2.312456e+01  2.312461e+01  5.360410e-01       OPTIMAL

#  TERMINATION: 0.536041 CPU SEC
#  INCUMBENT VALUE:  2.312461e+01
#  INCUMBENT POINT:  1.000000e+00  5.000000e+00  3.000000e+00  2.236068e+00
#  NUMBER OF RELAXATION REFINEMENTS:   5
#  OPTIMALITY GAP:   4.746061e-05 (ABS)
                     2.052387e-06 (REL)
\endverbatim

Other options can be modified to tailor the search, including output level, maximum number of iterations, tolerances, maximum CPU time, etc. These options can be modified through the public member mc::DOSEQGO::options. 
*/

//TODO: 
//- [OK]    Implement SBB method
//- [OK]    Enable use of polynomial models in relaxations
//- [OK]    Enable multiple rounds of PWR refinement
//- [OK]    Support MINLP
//- [OK]    Exclude linear variables from bound contraction
//- [TO DO] Enable KKT cuts and reduction constraints
//- [TO DO] Make it possible to add/remove a constraint from the model?
//- [TO DO] Exploit the degree of separability to introduce auxiliary variables

#ifndef MC__DOSEQGO_HPP
#define MC__DOSEQGO_HPP

#include <stdexcept>

#include "csearch_base.hpp"
#include "doseqslv_ipopt.hpp"
#include "odebnds_sundials.hpp"

#undef MC__DOSEQGO_DEBUG
//#undef MC__DOSEQGO_TRACE
#undef MC__DOSEQGO_SHOW_BREAKPTS

namespace mc
{

//! @brief C++ class for global optimization of MIDO using complete search
////////////////////////////////////////////////////////////////////////
//! mc::DOSEQGO is a C++ class for global optimization of nonlinear /
//! mixed-integer dynamic optimization using complete search. Bounds for
//! the nonlinear ODEs are obtained using mc::ODEBND, whereas
//! relaxations for the nonconvex participating terms are generated
//! using MC++. Further details can be found at: \ref page_DOSEQGO
////////////////////////////////////////////////////////////////////////
template < typename T, typename PMT=CModel<T>, typename PVT=CVar<T> >
class DOSEQGO:
  public virtual CSEARCH_BASE<T>,
  public virtual BASE_DO
{
public:

  // Typedef's
  typedef typename LPRELAX_BASE<T>::LP_STATUS LP_STATUS;
  typedef ODEBNDS_SUNDIALS< T, PMT, PVT >  t_ODEBNDS;
  typedef Ipopt::SmartPtr<mc::DOSEQSLV_IPOPT> t_DOSEQSLV;

  //! @brief DOSEQGO options
  struct Options
  {
    //! @brief Constructor
    Options():
      CSALGO(BB), CVATOL(1e-3), CVRTOL(1e-3), FEASTOL(1e-5),
      BLKDECUSE(true), BRANCHPT(SBB<T>::Options::OMEGA), BRANCHDMIN(5e-2),
      BRANCHVAR(SBB<T>::Options::RGREL), BRANCHSEL(0), SCOBCHUSE(false),
      SCOBCHVMAX(0), SCOBCHRTOL(1e-1), SCOBCHATOL(0e0), STGBCHDEPTH(0),
      STGBCHDRMAX(1), STGBCHWEIGHT(1./6.), PREPROC(true), DOMREDMIG(1e-10),
      DOMREDMAX(5), DOMREDTHRES(0.1), DOMREDBKOFF(1e-8), RELMETH(CHEB),
      LPALGO( LPRELAX_BASE<T>::LPALGO_DEFAULT ), LPPRESOLVE(-1),
      LPFEASTOL(1e-9), LPOPTIMTOL(1e-9), MIPRELGAP(1e-7), MIPABSGAP(1e-7),
      PRESOS2BIGM(-1.), CMODEL(), CMODMIG(1e-10), CMODPROP(2), CMODCUTS(0),
      CMODDMAX(1e20), CMODDEPS(0), CMODRED(APPEND), CMODWARMS(false),
      CMODRTOL(2e-1), CMODATOL(1e-8), CMODJOINT(false), MAXITER(0),
      MAXCPU(7.2e3), DISPLAY(2), MIPDISPLAY(0), MIPFILE(""), POLIMG(),
      DOSEQSLV(), ODEBNDS()
      { ODEBNDS.DISPLAY = 0;
        CMODEL.MIXED_IA = true; }
    //! @brief Assignment operator
    Options& operator= ( Options&options ){
        CSALGO       = options.CSALGO;
        RELMETH      = options.RELMETH;
        PREPROC      = options.PREPROC;
        DOMREDMIG    = options.DOMREDMIG;
        DOMREDMAX    = options.DOMREDMAX;
        DOMREDTHRES  = options.DOMREDTHRES;
        DOMREDBKOFF  = options.DOMREDBKOFF;
        CVATOL       = options.CVATOL;
        CVRTOL       = options.CVRTOL;
        FEASTOL      = options.FEASTOL;
        BLKDECUSE    = options.BLKDECUSE;
        BRANCHPT     = options.BRANCHPT;
        BRANCHDMIN   = options.BRANCHDMIN;
        BRANCHVAR    = options.BRANCHVAR;
        BRANCHSEL    = options.BRANCHSEL;
        SCOBCHUSE    = options.SCOBCHUSE;
        SCOBCHVMAX   = options.SCOBCHVMAX;
        SCOBCHRTOL   = options.SCOBCHRTOL;
        SCOBCHATOL   = options.SCOBCHATOL;
        STGBCHDEPTH  = options.STGBCHDEPTH;
        STGBCHDRMAX  = options.STGBCHDRMAX;
        STGBCHWEIGHT = options.STGBCHWEIGHT;
        LPALGO       = options.LPALGO;
        LPPRESOLVE   = options.LPPRESOLVE;
        LPFEASTOL    = options.LPFEASTOL;
        LPOPTIMTOL   = options.LPOPTIMTOL;
        MIPRELGAP    = options.MIPRELGAP;
        MIPABSGAP    = options.MIPABSGAP;
        PRESOS2BIGM  = options.PRESOS2BIGM;
        CMODEL       = options.CMODEL;
        CMODMIG      = options.CMODMIG;
        CMODPROP     = options.CMODPROP;
        CMODCUTS     = options.CMODCUTS;
        CMODDMAX     = options.CMODDMAX;
        CMODDEPS     = options.CMODDEPS;
        CMODRED      = options.CMODRED;
        CMODWARMS    = options.CMODWARMS;
        CMODRTOL     = options.CMODRTOL;
        CMODATOL     = options.CMODATOL;
        CMODJOINT    = options.CMODJOINT;
        MAXITER      = options.MAXITER;
        MAXCPU       = options.MAXCPU;
        DISPLAY      = options.DISPLAY;
        MIPDISPLAY   = options.MIPDISPLAY;
        MIPFILE      = options.MIPFILE;
        POLIMG       = options.POLIMG;
        DOSEQSLV     = options.DOSEQSLV;
        ODEBNDS      = options.ODEBNDS;
        return *this;
      }
    //! @brief Complete search method
    enum CS{
      BB=0	//!< Spatial branch-and-bound
    };
    //! @brief Relaxation strategy
    enum RELAX{
      DRL=0,	//!< Decomposition-relaxation-linearization (Tawarmalani & Sahinidis)
      CHEB, 	//!< Chebyshev-derived relaxations, controlled by parameters CMODPROP and CMODCUT
      HYBRID    //!< Combination of DRL and CHEB
    };
    //! @brief Reduction strategy
    enum REDUC{
      NONE=0,	    //!< Do not use Chebyshev-reduction constraints
      SUBSTITUTE,	//!< Substitute Chebyshev-reduction constraints into the other constraints
      APPEND	    //!< Append Chebyshev-reduction constraints to the other constraints
    };
    //! @brief Complete search algorithm
    CS CSALGO;
    //! @brief Convergence absolute tolerance
    double CVATOL;
    //! @brief Convergence relative tolerance
    double CVRTOL;
    //! @brief Feasibility tolerance 
    double FEASTOL;
    //! @brief Whether to use block-triangular decomposition for branching variable preselection
    bool BLKDECUSE;
    //! @brief Branching strategy for bisection point
    int BRANCHPT;
    //typename SBB<T>::Options::STRATEGY BRANCHPT;
    //! @brief Relative tolerance within which a variable is considered to be at one of its bounds, i.e. excluded from branching variable selection
    double BRANCHDMIN;
    //! @brief Branching-variable selection criterion
    int BRANCHVAR;
    //typename SBB<T>::Options::CRITERION BRANCHVAR;
    //! @brief Branching-variable selection user-function
    //int BRANCHSEL;
    typename SBB<T>::Options::SELECTION BRANCHSEL;
    //! @brief Whether or not to preselect branching variables based on scores
    bool SCOBCHUSE;
    //! @brief Maximum variable subset from score branching selection
    unsigned SCOBCHVMAX;
    //! @brief Relative tolerance for score branching selection
    double SCOBCHRTOL;
    //! @brief Absolute tolerance for score branching selection
    double SCOBCHATOL;
    //! @brief Maximum depth for strong branching interruption (0: no strong branching)
    unsigned STGBCHDEPTH;
    //! @brief Maximum domain reduction loops during strong branching (0: no restriction)
    unsigned STGBCHDRMAX;
    //! @brief Weighting (between 0 and 1) used to account for the left and right nodes in strong branching
    double STGBCHWEIGHT;
    //! @brief Whether or not to preprocess the root node (local optimization, domain reduction)
    bool PREPROC;
    //! @brief Minimal diameter of variable domain for contraction
    double DOMREDMIG;
    //! @brief Maximum number of domain reduction rounds
    unsigned DOMREDMAX;
    //! @brief Threshold for repeating domain reduction (minimum domain reduction ratio)
    double DOMREDTHRES;
    //! @brief Backoff of reduced variable bounds to compensate for numerical errors
    double DOMREDBKOFF;
    //! @brief Relaxation method
    RELAX RELMETH;
    //! @brief LP algorithm
    int LPALGO;
    //! @brief LP presolve
    int LPPRESOLVE;
    //! @brief Tolerance on LP feasibility
    double LPFEASTOL;
     //! @brief Tolerance on LP optimality
    double LPOPTIMTOL;
    //! @brief Tolerance on relative MIP gap
    double MIPRELGAP;
    //! @brief Tolerance on absolute MIP gap
    double MIPABSGAP;
    //! @brief Parameter controlling SOS2 reformulations
    double PRESOS2BIGM;
    //! @brief CModel options
    typename SCModel<T>::Options CMODEL;
    //! @brief Minimal magnitude of coefficient in Chebyhev model
    double CMODMIG;
    //! @brief Chebyhev model propagation order (0: no propag.)
    unsigned CMODPROP;
    //! @brief Chebyhev model cut order (0: same as propag.)
    unsigned CMODCUTS;
    //! @brief Chebyhev model maximum diameter for cut generation
    double CMODDMAX;
     //! @brief Dependent Chebyhev model order (0: no CM)
    unsigned CMODDEPS;
     //! @brief Reduction method
    REDUC CMODRED;
    //! @brief Whether or not to warm start the reduced-space estimators with linear estimators
    bool CMODWARMS;
    //! @brief Reduced-space Chebyshev model relative tolerance
    double CMODRTOL;
    //! @brief Reduced-space Chebyshev model absolute tolerance
    double CMODATOL;
    //! @brief Whether or not to enforce convergence for all reduced-space Chebyshev models jointly
    double CMODJOINT;
    //! @brief Maximum number of iterations (0: no limit )
    unsigned MAXITER;
    //! @brief Maximum run time (seconds)
    double MAXCPU;
    //! @brief Display level for solver
    int DISPLAY;
    //! @brief Display level for MIP
    int MIPDISPLAY;
    //! @brief Name of output file for optimization model
    std::string MIPFILE;
    //! @brief PolImg (polyhedral relaxation) options
    typename PolImg<T>::Options POLIMG;
    //! @brief DOSEQSLV_IPOPT (local optinmization) options
    typename DOSEQSLV_IPOPT::Options DOSEQSLV;
    //! @brief AEBND (algebraic equation bounder) options
    typename t_ODEBNDS::Options ODEBNDS;
    //! @brief Display
    void display
      ( std::ostream&out ) const;
  };
  //! @brief DOSEQGO options
  Options options;

  //! @brief Class managing exceptions for DOSEQGO
  class Exceptions
  {
  public:
    //! @brief Enumeration type for DOSEQGO exception handling
    enum TYPE{
      MOBJ=1,		//!< Optimization problem may not have more than one objective functions
      SETUP,		//!< Incomplete setup before a solve
      INTERN=-33	//!< Internal error
    };
    //! @brief Constructor for error <a>ierr</a>
    Exceptions( TYPE ierr ) : _ierr( ierr ){}
    //! @brief Inline function returning the error flag
    int ierr(){ return _ierr; }
    //! @brief Inline function returning the error description
    std::string what(){
      switch( _ierr ){
      case MOBJ:
        return "DOSEQGO::Exceptions  Optimization problem may not have more than one objective functions";
      case SETUP:
        return "DOSEQGO::Exceptions  Incomplete setup before a solve";
      case INTERN: default:
        return "DOSEQGO::Exceptions  Internal error";
      }
    }
  private:
    TYPE _ierr;
  };

  //! @brief Structure holding current statistics
  struct Stats{
    void reset()
      { tDEPBND = tPOLIMG = tLPSOL = tLPSET = tLOCSOL = tALL = 0.;
        nDEPBND = nLPSOL = nLOCSOL = 0; }
    void display
      ( std::ostream&os=std::cout )
      { os << std::fixed << std::setprecision(2)
           << "#  TOTAL:  " << tALL    << " CPU SEC" << std::endl
           << "#  POLIMG: " << tPOLIMG << " CPU SEC" << std::endl
           << "#  LPSET:  " << tLPSET  << " CPU SEC" << std::endl;
        if( nLPSOL )  os << "#  LPSOL:  " << tLPSOL  << " CPU SEC  (" << nLPSOL << ")" << std::endl; 
        if( nLOCSOL ) os << "#  LOCSOL: " << tLOCSOL << " CPU SEC  (" << nLOCSOL << ")" << std::endl;
        if( nDEPBND ) os << "#  DEPBND: " << tDEPBND << " CPU SEC  (" << nDEPBND << ")" << std::endl;
        os << std::endl; }
    double tALL;
    double tDEPBND;
    unsigned nDEPBND;
    double tPOLIMG;
    double tLPSOL;
    unsigned nLPSOL;
    double tLPSET;
    double tLOCSOL;
    unsigned nLOCSOL;
  } stats;

  // Overloading stdout operator
  template <typename U> friend std::ostream& operator<<
    ( std::ostream&os, const DOSEQGO<U>& );

protected:
  // Typedef's
  typedef std::pair< unsigned, std::map< unsigned, unsigned > > t_expmon;
  typedef std::map< t_expmon, double > t_coefmon;
  typedef std::set< FFVar*, lt_FFVar > t_FFVar;

  //! @brief Members of CSEARCH_BASE
  using CSEARCH_BASE<T>::_nrdep;
  using CSEARCH_BASE<T>::_nrvar;
  using CSEARCH_BASE<T>::_nvar;
  using CSEARCH_BASE<T>::_var;
  using CSEARCH_BASE<T>::_var_fperm;
  using CSEARCH_BASE<T>::_var_rperm;
  using CSEARCH_BASE<T>::_var_lin;
  using CSEARCH_BASE<T>::_var_excl;
  using CSEARCH_BASE<T>::_nctr;
  using CSEARCH_BASE<T>::_ctr;
  using CSEARCH_BASE<T>::_neq;
  using CSEARCH_BASE<T>::_eq;
  using CSEARCH_BASE<T>::_obj;
  using CSEARCH_BASE<T>::_Indxdep;
  using CSEARCH_BASE<T>::_CMndxdep;
  using CSEARCH_BASE<T>::_p_inc;
  using CSEARCH_BASE<T>::_f_inc;
  using CSEARCH_BASE<T>::_ignore_deps;
  using CSEARCH_BASE<T>::_CMrenv;
  using CSEARCH_BASE<T>::_CMrvar;
  using CSEARCH_BASE<T>::_CMrdep;
  using CSEARCH_BASE<T>::_Irvar;
  using CSEARCH_BASE<T>::_Irdep;
  using CSEARCH_BASE<T>::_POLvar;
  using CSEARCH_BASE<T>::_odisp;
  using CSEARCH_BASE<T>::_LDISP;
  using CSEARCH_BASE<T>::_IPREC;
  using CSEARCH_BASE<T>::_DPREC;

  //! @brief Members of LPRELAX_BASE
  using LPRELAX_BASE<T>::_solve_LPmodel;
  using LPRELAX_BASE<T>::_set_LPrelax;

  //! @brief Members of BASE_DO
  //using BASE_DO::_ns;

  //! @brief Local NLP search
  t_DOSEQSLV _DOSEQSLV;

  //! @brief Differential equation bounder
  t_ODEBNDS _ODEBNDS;

public:
  //! @brief Constructor
  DOSEQGO()
    : _issetup(false), _PM_ENV(0)
    { _DOSEQSLV = new DOSEQSLV_IPOPT; }

  //! @brief Destructor
  virtual ~DOSEQGO()
    { delete _PM_ENV;
      for( auto it=_PM_STA.begin(); it!=_PM_STA.end(); ++it ) delete[] *it; }

  //! @brief Setup DAG for cost and constraint evaluation
  bool setup
    ( std::ostream&os=std::cout );

  //! @brief Solve optimization model to global optimality within variable range <a>P</a> and for variable types <a>tvar</a>
  int solve
    ( const T*P, const T*X=0, const unsigned*tvar=0, const double*p0=0, 
      std::ostream&os=std::cout );

  //! @brief Solve (continuous) optimization model to local optimality using IPOPT in variable range <a>P</a> and from initial point <a>p0</a> -- return value is IPOPT status
  int local
    ( const T*P, const double*p0=0, const bool reset=true )
    { return CSEARCH_BASE<T>::_local( _DOSEQSLV, stats, P, p0, reset ); }

  //! @brief Local solution information after last NLP optimization
  const DOSEQSLV_IPOPT::SOLUTION& get_local_solution() const
    { return _DOSEQSLV->solution(); }

  //! @brief Local solver
  const t_DOSEQSLV& local_solver() const
    { return _DOSEQSLV; }

  //! @brief Solve relaxed optimization model (or relaxed feasibility model if <a>feastest</a> is true) within variable range <a>P</a> and for variable types <a>tvar</a>
  LP_STATUS relax
    ( const T*P, const T*X=0, const unsigned*tvar=0, const unsigned refine=0,
      const bool reset=true, const bool feastest=false );

  //! @brief Solve bound contraction problems from relaxed model, starting with variable range <a>P</a>, for variable types <a>tvar</a> and incumbent value <a>inc</a> (or constraint backoff if <a>feastest</a> is true), and using the options specified in <a>DOSEQGO::Options::DOMREDMAX</a> and <a>DOSEQGO::Options::DOMREDTHRES</a> -- returns updated variable bounds <a>P</a>, and number of iterative refinements <a>nred</a> 
  LP_STATUS contract
    ( T*P, unsigned&nred, const T*X=0, const unsigned*tvar=0, const double*inc=0,
      const bool reset=true, const bool feastest=false );

  //! @brief Value of DAG variable <a>X</a> after last LP optimization
  using LPRELAX_BASE<T>::get_variable;
  //! @brief Optimal cost value after last LP optimization
  using LPRELAX_BASE<T>::get_objective;
  //! @brief Status after last LP optimization
  using LPRELAX_BASE<T>::get_status;
  //! @brief Pointer to last LP optimization model and solution
  using LPRELAX_BASE<T>::get_relaxed_model;

protected:
  //! @brief Flag for setup function
  bool _issetup;

  //! @brief Polynomial model environment for parametric IVP bounding
  PMT* _PM_ENV;

  //! @brief Polynomial model of states in parametric IVP
  std::vector<PVT*> _PM_STA;

  //! @brief Polynomial model of parameters participating in parametric IVP
  std::vector<PVT> _PM_PAR;

  //! @brief Set local optimizer
  virtual void _set_SLVLOC
    ();

  //! @brief Get local optimum
  const double* _get_SLVLOC
    ( const double*p );

  //! @brief Compute dependent Chebyshev variables
  bool _get_depcheb
    ( const bool reset=true );
  //! @brief Compute dependent interval variables
  bool _get_depbnd
    ( const bool reset=true );

  //! @brief User-function to subproblems in SBB
  typename SBB<T>::STATUS subproblems
    ( const typename SBB<T>::TASK task, SBBNode<T>*node,
      std::vector<double>&p, double&f, const double INC )
    { return CSEARCH_BASE<T>::_subproblems( _DOSEQSLV, options, stats, task, node, p, f, INC ); }

  //! @brief User-function to subproblems in SBP
  virtual typename SBP<T>::STATUS assess
    ( SBPNode<T>*node )
    { throw Exceptions( Exceptions::INTERN ); }

  //! @brief Initialize display
  void _display_pwl_init
    ()
    { throw Exceptions( Exceptions::INTERN ); }

  //! @brief Final display
  void _display_pwl_final
    ( const unsigned iter )
    { throw Exceptions( Exceptions::INTERN ); }

  //! @brief Conversion from dense Chebyshev model into sparse Chebyshev model
  SCVar<T> _dense2sparse
    ( const CVar<T>&CMx );

  //! @brief Conversion from dense Taylor model into sparse Chebyshev model
  SCVar<T> _dense2sparse
    ( const TVar<T>&TMx );

private:
  //! @brief Bounds on parameters and state dependents
  std::vector<T> _bPX;
  //! @brief Types of parameters and state dependents
  std::vector<unsigned> _tPX;
  //! @brief Values of parameters and state dependents
  std::vector<double> _vPX;

  //! @brief Conversion from dense Taylor model into sparse Chebyshev model
  bool _init
    ( const T*P, const T*X, const unsigned*tvar, const double*p0=0 );

  //! @brief Private methods to block default compiler methods
  DOSEQGO
    (const DOSEQGO&);
  DOSEQGO& operator=
    (const DOSEQGO&);
};

template <typename T, typename PMT, typename PVT>
inline bool
DOSEQGO<T,PMT,PVT>::setup
( std::ostream&os )
{
  _issetup = false;

  // Generate dependency information
  if( !set_depend() ) return false;
  _pDAG->output( _pDAG->subgraph( 1, BASE_DO::_fct.data() ) );
#ifdef MC__DOSEQGO__DEBUG
  { int dum; std::cout << "PAUSED--"; std::cin >> dum; }
#endif

  //_valFCTSTA.resize( nFCTSTA );
  //_valDFCTSTA.resize( nDEPSTA*nFCTSTA );
  //_valFCTPAR.resize( _nf-nFCTSTA );
  //_valDFCTPAR.resize( _np*(_nf-nFCTSTA) );

  // Setup parametric IVP bounder
  _ODEBNDS.set( *this );
  _ODEBNDS.set_time( _ns, _dT.data(), _pT );
  _ODEBNDS.set_parameter( _vSTADEP.size(), _vSTADEP.data() );
  _ODEBNDS.reset_function();
  //_ODEBNDS.set_function( _ns, nFCTSTA, _vFCTSTA.data() );
  _ODEBNDS.options = options.ODEBNDS;

  _PM_ENV = new PMT( _vSTADEP.size(), options.CMODDEPS );
  _PM_PAR.resize( _vSTADEP.size() );
  for( auto it=_PM_STA.begin(); it!=_PM_STA.end(); ++it )
    delete[] *it;
  _PM_STA.assign( _ns+1, 0 );

  // Independent decision variables: parameters
  _nrvar = _np;
  _var.assign( _pP, _pP+_nrvar );

  // Dependent decision variables: states
  _nrdep = _mapFCTSTA.size();
  for( auto it=_mapFCTSTA.begin(); it!=_mapFCTSTA.end(); ++it )
    _var.push_back( it->second );
  _nvar = _var.size();
  _var_fperm.resize(_nvar);
  _var_rperm.resize(_nvar);
  for( unsigned i=0; i<_nvar; i++ )
    _var_fperm[i] = _var_rperm[i] = i;

  // Cost function
  std::get<0>(_obj).clear();
  std::get<1>(_obj).clear();
  std::get<2>(_obj).clear();
  if( std::get<0>(BASE_DO::_obj).size() > 1 )
    throw Exceptions( Exceptions::MOBJ );
  else if( std::get<0>(BASE_DO::_obj).size() ){
    std::get<0>(_obj).push_back( std::get<0>(BASE_DO::_obj)[0] );
    std::get<1>(_obj).push_back( BASE_DO::_fct[0] );
    std::get<2>(_obj).push_back( std::get<2>(BASE_DO::_obj)[0] );
  }

  // Constraint functions
  std::get<0>(_ctr).clear();
  std::get<1>(_ctr).clear();
  std::get<2>(_ctr).clear();
  std::get<3>(_ctr).clear();
  for( unsigned ic=0; ic<std::get<0>(BASE_DO::_ctr).size(); ++ic ){
    std::get<0>(_ctr).push_back( std::get<0>(BASE_DO::_ctr)[ic] );
    std::get<1>(_ctr).push_back( BASE_DO::_fct[ic+1] );
    std::get<2>(_ctr).push_back( std::get<2>(BASE_DO::_ctr)[ic] );
    std::get<3>(_ctr).push_back( false );
  }
  _nctr = std::get<0>(_ctr).size();
  _eq.clear();
  _neq = _eq.size();

  // internal setup, including structure detection
  CSEARCH_BASE<T>::_setup( options, _pDAG, os );

  // Exclude dependent-state variables from branching
  _var_excl.clear();
  for( unsigned iv=0; iv<_nrdep; ++iv )
    _var_excl.insert( _nrvar+iv );

  // Identify linearly-participating (independent) variables
  FFDep fgdep;
  for( unsigned ic=0; ic<_nf; ic++ ){
    fgdep += _depF[ic];
    if( state_depend( ic ) ) fgdep += _fct[ic].dep();
  }
#ifdef MC__DOSEQGO__DEBUG
  std::cout << "DO <- " << fgdep << std::endl;
  //int dum; std::cin >> dum;
#endif
  _var_lin.clear();
  for( unsigned i=0; i<_nrvar; i++ ){ // only consider independents
    auto it = fgdep.dep().find( _var[i].id().second );
    if( it == fgdep.dep().end() || it->second )
      _var_lin.insert(i);
  }

  if( options.DISPLAY ){
    os << "LINEAR VARIABLES IN FUNCTIONS:            " << _var_lin.size() << std::endl
       << "NONLINEAR VARIABLES IN FUNCTIONS:         " << _nrvar-_var_lin.size() << std::endl;
  }

  _ignore_deps = false;
  _issetup = true;
  return _issetup;
}

template <typename T, typename PMT, typename PVT>
inline void
DOSEQGO<T,PMT,PVT>::_set_SLVLOC
()
{
  _DOSEQSLV->options = options.DOSEQSLV;
  _DOSEQSLV->set( *this );
  _DOSEQSLV->setup();
}

template <typename T, typename PMT, typename PVT>
inline const double*
DOSEQGO<T,PMT,PVT>::_get_SLVLOC
( const double*p )
{
  //_vPX.assign( _DOSEQSLV->solution().p, _DOSEQSLV->solution().p+_np );
  _vPX.assign( p, p+_np );
  std::vector<double*> xk( _ns+1, 0 );
  if( _DOSEQSLV->states( _vPX.data(), xk.data() ) != BASE_DE::NORMAL ){
    for( auto itx=xk.begin(); itx!=xk.end(); ++itx ) delete[] *itx;
    return 0;
  }

  auto it=_mapFCTSTA.begin();
  for( unsigned i=0; it!=_mapFCTSTA.end(); ++it, i++ )
    _vPX.push_back( xk[it->first.first+1][it->first.second] );
  for( auto itx=xk.begin(); itx!=xk.end(); ++itx ) delete[] *itx;
  return _vPX.data();
}

template <typename T, typename PMT, typename PVT>
inline bool
DOSEQGO<T,PMT,PVT>::_init
( const T*P, const T*X, const unsigned*tvar, const double*p0 )
{
  // Form vector of parameter and state bounds
  if( P ){
    _bPX.assign( P, P+_np );
    auto it=_mapFCTSTA.begin();
    for( unsigned i=0; it!=_mapFCTSTA.end(); ++it, i++ )
      _bPX.push_back( X? X[it->first.second]: T(-SBB<T>::INF,SBB<T>::INF) );
  }

  // Form vector of parameter and state types (continuous vs discrete)
  if( tvar ){
    _tPX.assign( tvar, tvar+_np );
    auto it=_mapFCTSTA.begin();
    for( unsigned i=0; it!=_mapFCTSTA.end(); ++it, i++ )
      _tPX.push_back( false ); // states are continuous variables
  }

  // Form vector of parameter and state values
  if( p0 ){
    _vPX.assign( p0, p0+_np );
    std::vector<double*> xk( _ns+1, 0 );
    if( _DOSEQSLV->states( p0, xk.data() ) != BASE_DE::NORMAL ){
      for( auto itx=xk.begin(); itx!=xk.end(); ++itx ) delete[] *itx;
      return false;
    }
    auto it=_mapFCTSTA.begin();
    for( unsigned i=0; it!=_mapFCTSTA.end(); ++it, i++ )
      _vPX.push_back( xk[it->first.first+1][it->first.second] );
    for( auto itx=xk.begin(); itx!=xk.end(); ++itx ) delete[] *itx;
  }

  return true;
}

template <typename T, typename PMT, typename PVT>
inline typename LPRELAX_BASE<T>::LP_STATUS
DOSEQGO<T,PMT,PVT>::relax
( const T*P, const T*X, const unsigned*tvar, const unsigned refine,
  const bool reset, const bool feastest )
{
  if( !_issetup || !_init( P, X, tvar ) ) throw Exceptions( Exceptions::SETUP );

  _ignore_deps = false;
  return CSEARCH_BASE<T>::_relax( options, stats, _bPX.data(), tvar?_tPX.data():0,
                                  refine, reset, feastest );
}

template <typename T, typename PMT, typename PVT>
inline typename LPRELAX_BASE<T>::LP_STATUS
DOSEQGO<T,PMT,PVT>::contract
( T*P, unsigned&nred, const T*X, const unsigned*tvar, const double*inc,
  const bool reset, const bool feastest )
{
  if( !_issetup || !_init( P, X, tvar ) ) throw Exceptions( Exceptions::SETUP );

  _ignore_deps = false;
  return CSEARCH_BASE<T>::_contract( options, stats, _bPX.data(), nred,
                                     tvar?_tPX.data():0, inc, reset, feastest );
}

template <typename T, typename PMT, typename PVT>
inline int
DOSEQGO<T,PMT,PVT>::solve
( const T*P, const T*X, const unsigned*tvar, const double*p0,
  std::ostream&os )
{
  if( !_issetup || !_init( P, X, tvar, p0 ) ) throw Exceptions( Exceptions::SETUP );

  switch( options.CSALGO ){
   case Options::BB: default:
     return CSEARCH_BASE<T>::_solve_sbb( _DOSEQSLV, options, stats, _bPX.data(),
                                         tvar?_tPX.data():0, p0?_vPX.data():0, os );
  }
}

template <typename T, typename PMT, typename PVT>
inline SCVar<T>
DOSEQGO<T,PMT,PVT>::_dense2sparse
( const CVar<T>&CMx )
{
  if( !CMx.env() )
    return SCVar<T>( CMx.coefmon().second[0] ) + CMx.remainder();

  typename SCVar<T>::t_coefmon SCmon;
  T SCrem = CMx.remainder();
  auto Cmon = CMx.coefmon();

  for( auto it=CMx.ndxmon().begin(); it!=CMx.ndxmon().end(); ++it ){
    // Dump into reminder term if too small
    if( std::fabs(Cmon.second[*it]) < options.CMODMIG ){
      SCrem += Cmon.second[*it] * (Op<T>::zeroone()*2-1);
      continue;
    }
    // Build monomial term
    typename SCVar<T>::t_expmon expmon( 0, std::map<unsigned,unsigned>() );
    for( unsigned k=0; k<CMx.nvar(); k++ ){
      const unsigned kexp = CMx.expmon(*it)[k];
      if( !kexp ) continue;
      expmon.first += kexp;
      expmon.second[k] = kexp;
    }
    SCmon[expmon] = Cmon.second[*it];
  }

  for( unsigned i=0; CMx.ndxmon().empty() && i<Cmon.first; i++ ){
    // Dump into reminder term if too small
    if( std::fabs(Cmon.second[i]) < options.CMODMIG ){
      SCrem += Cmon.second[i] * (Op<T>::zeroone()*2-1);
      continue;
    }
    // Build monomial term
    typename SCVar<T>::t_expmon expmon( 0, std::map<unsigned,unsigned>() );
    for( unsigned k=0; k<CMx.nvar(); k++ ){
      const unsigned kexp = CMx.expmon(i)[k];
      if( !kexp ) continue;
      expmon.first += kexp;
      expmon.second[k] = kexp;
    }
    SCmon[expmon] = Cmon.second[i];
  }

  SCVar<T> SCMx;
  SCMx.set( _CMrenv );
  SCMx.set( SCmon );
  return SCMx + SCrem;
}

template <typename T, typename PMT, typename PVT>
inline SCVar<T>
DOSEQGO<T,PMT,PVT>::_dense2sparse
( const TVar<T>&TMx )
{
  if( !TMx.env() )
    return SCVar<T>( TMx.coefmon().second[0] ) + TMx.remainder();

  SCVar<T> SCMx;
  auto Tmon = TMx.coefmon();

  for( unsigned i=0; i<Tmon.first; i++ ){
    // Dump into reminder term if too small
    if( Tmon.second[i] < options.CMODMIG ){
      SCMx += Tmon.second[i] * (Op<T>::zeroone()*2-1);
      continue;
    }
    // Build monomial term
    SCVar<T> term = Tmon.second[i];
    for( unsigned k=0; k<TMx.nvar(); k++ ){
      const unsigned kexp = TMx.expmon(i)[k];
      if( !kexp ) continue;
      term *= pow( _CMrvar[k], kexp );
    }
    SCMx += term;
  }

  return SCMx + TMx.remainder();
}

template <typename T, typename PMT, typename PVT>
inline bool
DOSEQGO<T,PMT,PVT>::_get_depcheb
( const bool reset )
{
#if defined (MC__DOSEQGO__DEBUG_CHEBDEPS)
  for( unsigned i=0; i<_nrvar; i++ )
    std::cout << "SCMrvar[" << i << "] =" << _CMrvar[i];
  for( unsigned i=0; i<_nrdep; i++ )
    std::cout << "SCMrdep[" << i << "] =" << _CMrdep[i];
#endif
  // Reset list of dependents
  if( reset ) _CMndxdep.clear();

  std::vector<SCVar<T>> CMrdep0 = _CMrdep; // local copy for convergence test later on

  int status = NORMAL;

  // Convert IVP particpating parameters into PVT format
  for( unsigned i=0; i<_nrvar; i++ ){
    _PM_PAR[i].set( _PM_ENV, i, _CMrvar[i].B() );
#if defined (MC__DOSEQGO__DEBUG_CHEBDEPS)
    std::cout << "CMrvar[" << i << "] =" << _PM_PAR[i].B();
    std::cout << "CVar #" << i << ": " << _PM_ENV->bndvar(i) << " "
              << _PM_ENV->refvar(i) << " " << _PM_ENV->scalvar(i) << " " << std::endl;
#endif
  }

  // Bound parametric ODE system
  try{
    status = _ODEBNDS.bounds( _PM_PAR.data(), _PM_STA.data() );
  }

  // ODE bounding unsuccessful?
  catch(...){
#ifdef MC__CSEARCH_SHOW_DEPS
    std::cout << "Dependents: Failed" << std::endl;
#endif
    status = FAILURE;
  }

  if( status != NORMAL && status != FAILURE )
    return true;

  // Convert IVP state enclosures from PVT format back into SCVar<T> format
  unsigned nsf = _ODEBNDS.final_stage();
  auto it=_mapFCTSTA.begin();
#ifdef MC__CSEARCH_SHOW_DEPS
  std::cout << "_CMndxdep = {";
#endif
  for( unsigned i=0; it!=_mapFCTSTA.end(); ++it, i++ ){
    if( it->first.first >= nsf ) continue;
    _CMrdep[i] = _dense2sparse( _PM_STA[it->first.first+1][it->first.second] );
    _CMndxdep.insert( _var_rperm[_nrvar+i] );
#ifdef MC__CSEARCH_SHOW_DEPS
    std::cout << " " << _var[_var_rperm[_nrvar+i]];
#endif
#ifdef MC__CSEARCH_DEBUG_BNDDEPS
    std::cout << "CMrdep[" << i << "] =" << _CMrdep[i];
    //std::cout << "CMrdep[" << i << "]@p* =" << _CMrdep[i].P(_p_inc.data())+_CMrdep[i].R()
    //          << " in " << _CMrdep[i].B() << std::endl;
    //std::cout << "CMrdep[" << i << "]@p* =" << _CMrdep[i].P(_p_inc.data()) << std::endl;
    //std::cout << "CMrdep[" << i << "] =" << _PM_STA[it->first.first+1][it->first.second];
    //std::cout << "CMrdep[" << i << "]@p* =" << _PM_STA[it->first.first+1][it->first.second].P(_p_inc.data()) << std::endl;
    //{int dum; std::cout << "PAUSED"; std::cin >> dum;}
#endif
  }
#ifdef MC__CSEARCH_SHOW_DEPS
  std::cout << " }" << std::endl;
#endif

/*
  // Append dependent showing improvement to _CMndxdep.
#ifdef MC__CSEARCH_SHOW_DEPS
  std::cout << "_CMndxdep = {";
#endif
  for( unsigned ib=0; ib<_AEBND.iblk(); ++ib )
    for( unsigned j=0; j<_AEBND.nblk(ib); j++ ){
#ifdef MC__CSEARCH_DEBUG_BNDDEPS
      std::cout << "diam CMrdep(" << _AEBND.rpdep(ib,j) << ") = "
                << Op<T>::diam( _CMrdep[_AEBND.rpdep(ib,j)].R() ) << " + "
                << Op<T>::diam( _CMrdep[_AEBND.rpdep(ib,j)].bndord( options.CMODPROP+1 ) )
                << std::endl;
      std::cout << "diam CMrdep0(" << _AEBND.rpdep(ib,j) << ") = "
                << Op<T>::diam( CMrdep0[_AEBND.rpdep(ib,j)].B() )
                << std::endl;
#endif
      if( Op<T>::diam( _CMrdep[_AEBND.rpdep(ib,j)].R() + _CMrdep[_AEBND.rpdep(ib,j)].bndord( options.CMODPROP+1 ) )
          < options.CMODATOL + options.CMODRTOL * Op<T>::diam( CMrdep0[_AEBND.rpdep(ib,j)].B() ) ){
        _CMrdep[_AEBND.rpdep(ib,j)].simplify();
        _CMndxdep.insert( _var_rperm[_nrvar+_AEBND.rpdep(ib,j)] );
#ifdef MC__CSEARCH_SHOW_DEPS
        std::cout << " " << _var[_var_rperm[_nrvar+_AEBND.rpdep(ib,j)]];
               // << "(" << _AEBND.rpdep(ib,j) << ") " << _Irdep[_AEBND.rpdep(ib,j)];
#endif
      }
    }
#ifdef MC__CSEARCH_SHOW_DEPS
  std::cout << " }" << std::endl;
#endif
*/
  return true;
}

template <typename T, typename PMT, typename PVT>
inline bool
DOSEQGO<T,PMT,PVT>::_get_depbnd
( const bool reset )
{
#if defined (MC__DOSEQGO_DEBUG_BNDDEPS)
  for( unsigned i=0; i<_nrdep; i++ )
    std::cout << "Irdep[" << i << "] =" << _Irdep[i];
#endif

  // Reset list of dependents
  if( reset) _Indxdep.clear();

  return true;
}

template <typename T, typename PMT, typename PVT>
inline void
DOSEQGO<T,PMT,PVT>::Options::display
( std::ostream&out ) const
{
  // Display DOSEQGO Options
  out << std::left;
  out << std::setw(60) << "  COMPLETE SEARCH METHOD";
  switch( CSALGO ){
   case BB: out << "BB" << std::endl; break;
  }
  out << std::setw(60) << "  CONVERGENCE ABSOLUTE TOLERANCE"
      << std::scientific << std::setprecision(1)
      << CVATOL << std::endl;
  out << std::setw(60) << "  CONVERGENCE RELATIVE TOLERANCE"
      << std::scientific << std::setprecision(1)
      << CVRTOL << std::endl;
  out << std::setw(60) << "  POLYHEDRAL RELAXATION APPROACH";
  switch( RELMETH ){
   case DRL:    out << "DRL"    << std::endl; break;
   case CHEB:   out << "CHEB"   << std::endl; break;
   case HYBRID: out << "HYBRID" << std::endl; break;
  }
  out << std::setw(60) << "  ORDER OF CHEBYSHEV MODEL PROPAGATION";
  switch( CMODPROP ){
   case 0:  out << "-\n"; break;
   default: out << CMODPROP << std::endl; break;
  }
  out << std::setw(60) << "  ORDER OF CHEBYSHEV-DERIVED CUTS";
  if( !CMODCUTS)
    switch( CMODPROP ){
     case 0:  out << "-\n"; break;
     default: out << CMODPROP << std::endl; break;
    }
  else 
    switch( CMODCUTS ){
     case 0:  out << "-\n"; break;
     default: out << std::min(CMODPROP,CMODCUTS) << std::endl; break;
    }
  out << std::setw(60) << "  USE OF CHEBYSHEV-REDUCTION CONSTRAINTS";
  switch( CMODRED ){
   case NONE:       out << "NONE" << std::endl; break;
   case SUBSTITUTE: out << "SUBSTITUTE" << std::endl; break;
   case APPEND:     out << "APPEND"     << std::endl; break;
  }
  out << std::setw(60) << "  ORDER OF CHEBYSHEV-REDUCTION CONSTRAINTS";
  switch( CMODDEPS ){
   case 0:  out << "-\n"; break;
   default: out << CMODDEPS << std::endl;
            out << std::setw(60) << "  CONVERGENCE THRESHOLD FOR REDUCED-SPACE CHEBYSHEV MODEL"
                << std::scientific << std::setprecision(1)
                << CMODRTOL << std::endl;
  }
  out << std::setw(60) << "  ROOT NODE PROPROCESSING"
      << (PREPROC?"Y\n":"N\n");
  out << std::setw(60) << "  MAXIMUM OPTIMIZATION-BASED REDUCTION LOOPS"
      << DOMREDMAX << std::endl;
  out << std::setw(60) << "  THRESHOLD FOR OPTIMIZATION-BASED REDUCTION LOOP"
      << std::fixed << std::setprecision(0)
      << DOMREDTHRES*1e2 << "%\n";
  out << std::setw(60) << "  BACKOFF FOR OPTIMIZATION-BASED REDUCTION"
      << std::scientific << std::setprecision(1)
      << DOMREDBKOFF << std::endl;
  out << std::setw(60) << "  MAXIMUM ITERATION COUNT";
  switch( MAXITER ){
   case 0:  out << "-\n"; break;
   default: out << MAXITER << std::endl; break;
  }
  out << std::setw(60) << "  MAXIMUM CPU TIME (SEC)"
      << std::scientific << std::setprecision(1)
      << MAXCPU << std::endl;
  out << std::setw(60) << "  DISPLAY LEVEL"
      << DISPLAY << std::endl;
}

template <typename T, typename PMT, typename PVT>
inline std::ostream&
operator <<
( std::ostream&out, const DOSEQGO<T,PMT,PVT>&NLP )
{
  out << std::right << std::endl
      << std::setfill('_') << std::setw(72) << " " << std::endl << std::endl << std::setfill(' ')
      << std::setw(55) << "DYNAMIC GLOBAL OPTIMIZATION IN CRONOS\n"
      << std::setfill('_') << std::setw(72) << " " << std::endl << std::endl << std::setfill(' ');

  // Display DOSEQGO Options
  //out << std::left << "SPATIAL BRANCH-AND-BOUND OPTIONS:\n\n";
  //NLP.SBB<T>::options.display( out );
  NLP.DOSEQGO<T,PMT,PVT>::options.display( out );

  out << std::setfill('_') << std::setw(72) << " " << std::endl << std::endl << std::setfill(' ');
  return out;
}

} // end namescape mc

#endif
