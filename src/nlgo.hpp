// Copyright (C) 2015-2016 Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

/*!
\page page_NLGO Nonlinear Global Optimization using MC++
\author Benoit C. Chachuat <tt>(b.chachuat@imperial.ac.uk)</tt> and OMEGA Research Group (http://www3.imperial.ac.uk/environmentenergyoptimisation)
\version 1.0
\date 2015
\bug No known bugs.

Consider a nonlinear optimization problem in the form:
\f{align*}
\mathcal{P}:\quad & \min_{x_1,\ldots,x_n}\ f(x_1,\ldots,x_n)\\
& {\rm s.t.}\ \ g_j(x_1,\ldots,x_n)\ \leq,=,\geq\ 0,\ \ j=1,\ldots,m\\
& \qquad x_i^L\leq x_i\leq x_i^U,\ \ i=1,\ldots,n\,,
\f}
where \f$f, g_1, \ldots, g_m\f$ are factorable, potentially nonlinear, real-valued functions; and \f$x_1, \ldots, x_n\f$ can be either continuous or integer decision variables. The class mc::NLGO solves such NLP or MINLP problems to global optimality using complete search. Two main methods are implemented in mc::NLGO:
- spatial branch-and-bound search
- hierarchy of semi-linear relaxations
.
In both methods, relaxations are generated for the nonlinear or nonconvex participating terms using various arithmetics in <A href="https://projects.coin-or.org/MCpp">MC++</A>.

\section sec_NLGO_setup How do I setup my optimization model?

Consider the following NLP:
\f{align*}
  \max_{\bf p}\ & p_1\,p_4\,(p_1+p_2+p_3)+p_3 \\
  \text{s.t.} \ & p_1\,p_2\,p_3\,p_4 \geq 25 \\
  & p_1^2+p_2^2+p_3^2+p_4^2 = 40 \\
  & 1 \leq p_1,p_2,p_3,p_4 \leq 5\,.
\f}

First, we define an mc::NLGO class as below:

\code
  mc::NLGO NLP;
\endcode

Next, we set the variables and objective/constraint functions by creating a direct acyclic graph (DAG) of the problem: 

\code
  #include "NLGO.hpp"
  
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
  NLP.options.CSALGO = mc::NLGO<I>::Options::PWL;
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

Other options can be modified to tailor the search, including output level, maximum number of iterations, tolerances, maximum CPU time, etc. These options can be modified through the public member mc::NLGO::options. 
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

#ifndef MC__NLGO_HPP
#define MC__NLGO_HPP

#include <stdexcept>

#include "csearch_base.hpp"
#include "nlpslv_ipopt.hpp"
#include "aebnd.hpp"

#undef MC__NLGO_DEBUG
//#undef MC__NLGO_TRACE
#undef MC__NLGO_SHOW_BREAKPTS

namespace mc
{

//! @brief C++ class for global optimization of MINLP using complete search
////////////////////////////////////////////////////////////////////////
//! mc::NLGO is a C++ class for global optimization of NLP and
//! MINLP using complete search. Relaxations for the nonlinea or
//! nonconvex participating terms are generated using MC++. Further
//! details can be found at: \ref page_NLGO
////////////////////////////////////////////////////////////////////////
template < typename T >
class NLGO:
  public virtual CSEARCH_BASE<T>,
  public virtual BASE_NLP
{
public:

  // Typedef's
  typedef typename LPRELAX_BASE<T>::LP_STATUS LP_STATUS;
  typedef AEBND< T, SCModel<T>, SCVar<T> >  t_AEBND;
  typedef Ipopt::SmartPtr<mc::NLPSLV_IPOPT> t_NLPSLV;

  //! @brief NLGO options
  struct Options
  {
    //! @brief Constructor
    Options():
      CSALGO(BB), CVATOL(1e-3), CVRTOL(1e-3), FEASTOL(1e-7),
      BRANCHPT(SBB<T>::Options::OMEGA), BRANCHDMIN(1e-10),
      BRANCHVAR(SBB<T>::Options::RGREL), BRANCHSEL(0), SCOBCHMETH(NOSCO),
      SCOBCHVMAX(0), SCOBCHRTOL(0e0), SCOBCHATOL(0e0), STGBCHDEPTH(0),
      STGBCHDRMAX(1), STGBCHWEIGHT(1./6.), PREPROC(true), DOMREDMIG(1e-7),
      DOMREDMAX(5), DOMREDTHRES(0.1), DOMREDBKOFF(1e-8), RELMETH(CHEB),
      LPALGO( LPRELAX_BASE<T>::LPALGO_DEFAULT ), LPPRESOLVE(-1),
      LPFEASTOL(1e-9), LPOPTIMTOL(1e-9), MIPRELGAP(1e-7), MIPABSGAP(1e-7),
      PRESOS2BIGM(-1.), CMODEL(), CMODPROP(2), CMODCUTS(0),
      CMODDMAX(1e20), DEPSUSE(true), CMODDEPS(0), CMODRED(APPEND),
      CMODWARMS(false), CMODRTOL(2e-1), CMODATOL(1e-8), CMODJOINT(false),
      NCOCUTS(false), NCOMETH(ASA), MAXITER(0), MAXCPU(7.2e3), DISPLAY(2),
      MIPDISPLAY(0), MIPFILE(""), POLIMG(), NLPSLV(), AEBND(0)
      { //AEBND.BOUNDER     = t_AEBND::Options::ALGORITHM::GS;
        //AEBND.BLKDEC      = t_AEBND::Options::DECOMPOSITION::DIAG;
        CMODEL.MIXED_IA   = true;
        CMODEL.MIN_FACTOR = 1e-8; }
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
        BRANCHPT     = options.BRANCHPT;
        BRANCHDMIN   = options.BRANCHDMIN;
        BRANCHVAR    = options.BRANCHVAR;
        BRANCHSEL    = options.BRANCHSEL;
        SCOBCHMETH   = options.SCOBCHMETH;
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
        CMODPROP     = options.CMODPROP;
        CMODCUTS     = options.CMODCUTS;
        CMODDMAX     = options.CMODDMAX;
        DEPSUSE      = options.DEPSUSE;
        CMODDEPS     = options.CMODDEPS;
        CMODRED      = options.CMODRED;
        CMODWARMS    = options.CMODWARMS;
        CMODRTOL     = options.CMODRTOL;
        CMODATOL     = options.CMODATOL;
        CMODJOINT    = options.CMODJOINT;
        NCOCUTS      = options.NCOCUTS;
        NCOMETH      = options.NCOMETH;
        MAXITER      = options.MAXITER;
        MAXCPU       = options.MAXCPU;
        DISPLAY      = options.DISPLAY;
        MIPDISPLAY   = options.MIPDISPLAY;
        MIPFILE      = options.MIPFILE;
        POLIMG       = options.POLIMG;
        NLPSLV       = options.NLPSLV;
        AEBND        = options.AEBND;
        return *this;
      }
    //! @brief Complete search method
    enum CS{
      BB=0,	//!< Spatial branch-and-bound
      PWL	//!< Hierarchy of piecewise-linear relaxations
    };
    //! @brief Relaxation strategy
    enum RELAX{
      DRL=0, //!< Decomposition-relaxation-linearization (Tawarmalani & Sahinidis)
      CHEB,  //!< Chebyshev-derived relaxations, controlled by parameters CMODPROP and CMODCUT
      HYBRID //!< Combination of DRL and CHEB
    };
    //! @brief Reduced-space strategy
    enum REDUC{
      NOREDUC=0, //!< Do not use Chebyshev-reduction constraints
      APPEND	 //!< Append Chebyshev-reduction constraints to the other constraints
    };
    //! @brief Scoring strategy
    enum SCORE{
      NOSCO=0, //!< No branching scores
      BLKSCO,  //!< Branching scores based on block-triangular decomposition of embedded equation systems
      CHEBSCO  //!< Branching scores based on nonlinear and remainder terms in chenyshev model inclusions
    };
    //! @brief Sensitivity strategy
    enum SENS{
      FSA=0,      //!< Forward sensitivity analysis
      ASA         //!< Adjoint sensitivity analysis
    };
    //! @brief Complete search algorithm
    CS CSALGO;
    //! @brief Convergence absolute tolerance
    double CVATOL;
    //! @brief Convergence relative tolerance
    double CVRTOL;
    //! @brief Feasibility tolerance 
    double FEASTOL;
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
    SCORE SCOBCHMETH;
    //! @brief Maximum variable subset from score branching selection (0: no restriction)
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
    //! @brief Minimal diameter of variable domain for branching or contraction
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
    //! @brief Chebyhev model propagation order (0: no propag.)
    unsigned CMODPROP;
    //! @brief Chebyhev model cut order (0: same as propag.)
    unsigned CMODCUTS;
    //! @brief Chebyhev model maximum diameter for cut generation
    double CMODDMAX;
    //! @brief Whether to exploit dependent variables in equation subsystems
    bool DEPSUSE;
     //! @brief Reduced-space Chebyhev model order (0: no reduction)
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
    //! @brief Whether or not to add NCO cuts
    bool NCOCUTS;
    //! @brief NCO method
    SENS NCOMETH;
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
    //! @brief NLPSLV_IPOPT (local optinmization) options
    typename NLPSLV_IPOPT::Options NLPSLV;
    //! @brief AEBND (algebraic equaiton bounder) options
    typename t_AEBND::Options AEBND;
    //! @brief Display
    void display
      ( std::ostream&out ) const;
  };
  //! @brief NLGO options
  Options options;

  //! @brief Class managing exceptions for NLGO
  class Exceptions
  {
  public:
    //! @brief Enumeration type for NLGO exception handling
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
        return "NLGO::Exceptions  Optimization problem may not have more than one objective functions";
      case SETUP:
        return "NLGO::Exceptions  Incomplete setup before a solve";
      case INTERN: default:
        return "NLGO::Exceptions  Internal error";
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
    ( std::ostream&os, const NLGO<U>& );

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
  using CSEARCH_BASE<T>::_fct_lin;
  using CSEARCH_BASE<T>::_nctr;
  using CSEARCH_BASE<T>::_ctr;
  using CSEARCH_BASE<T>::_neq;
  using CSEARCH_BASE<T>::_eq;
  using CSEARCH_BASE<T>::_obj;
  using CSEARCH_BASE<T>::_ndxdep;
  using CSEARCH_BASE<T>::_Indxdep;
  using CSEARCH_BASE<T>::_CMndxdep;
  using CSEARCH_BASE<T>::_p_inc;
  using CSEARCH_BASE<T>::_f_inc;
  using CSEARCH_BASE<T>::_ignore_deps;
  using CSEARCH_BASE<T>::_CMobj;
  using CSEARCH_BASE<T>::_CMctr;
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

  //! @brief Local NLP search
  t_NLPSLV _NLPSLV;

  //! @brief Implicit equation bounder
  t_AEBND _AEBND;

public:
  //! @brief Constructor
  NLGO()
    : CSEARCH_BASE<T>(), _issetup(false)
    { _NLPSLV = new NLPSLV_IPOPT; }

  //! @brief Destructor
  virtual ~NLGO()
    {}

  //! @brief Setup DAG for cost and constraint evaluation
  void setup
    ( const T*P, const unsigned*tvar=0, std::ostream&os=std::cout );

  //! @brief Solve optimization model to global optimality within variable range <a>P</a> and for variable types <a>tvar</a>
  int solve
    ( const T*P, const unsigned*tvar=0, const double*p0=0,
      std::ostream&os=std::cout );

  //! @brief Local solution information after last NLP optimization
  const SOLUTION_NLP& get_local_solution() const
    { return _NLPSLV->solution(); }

  //! @brief Solve relaxed optimization model (or relaxed feasibility model if <a>feastest</a> is true) within variable range <a>P</a> and for variable types <a>tvar</a>
  LP_STATUS relax
    ( const T*P, const unsigned*tvar=0, const unsigned refine=0,
      const bool reset=true, const bool feastest=false );

  //! @brief Solve bound contraction problems from relaxed model, starting with variable range <a>P</a>, for variable types <a>tvar</a> and incumbent value <a>inc</a> (or constraint backoff if <a>feastest</a> is true), and using the options specified in <a>NLGO::Options::DOMREDMAX</a> and <a>NLGO::Options::DOMREDTHRES</a> -- returns updated variable bounds <a>P</a>, and number of iterative refinements <a>nred</a> 
  LP_STATUS contract
    ( T*P, unsigned&nred, const unsigned*tvar=0, const double*inc=0,
      const bool reset=true, const bool feastest=false );

  //! @brief Global solution information after last NLP optimization
  std::pair<double,const double*> get_global_solution() const
    { return std::make_pair(_f_inc,_p_inc.data()); }

  //! @brief Solve (continuous) optimization model to local optimality using IPOPT in variable range <a>P</a> and from initial point <a>p0</a> -- return value is IPOPT status
  int local
    ( const T*P, const double*p0=0, const bool reset=true )
    { return CSEARCH_BASE<T>::_local( _NLPSLV, stats, P, p0, reset ); }

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

  //! @brief Bounds on all variables
  std::vector<T> _Ivar;
  //! @brief Types of all variables
  std::vector<unsigned> _tvar;
  //! @brief Values of all variables
  std::vector<double> _Dvar;

  //! @brief Initialize values, bounds and types for all variables
  bool _init
    ( const T*P, const unsigned*tvar, const double*p0=0 );

  //! @brief Set local optimizer
  virtual void _set_SLVLOC
    ();

  //! @brief Get local optimum
  const double* _get_SLVLOC
    ( const SOLUTION_NLP&locopt );

  //! @brief Compute dependent Chebyshev variables
  bool _get_depcheb
    ( const bool reset=true );
  //! @brief Compute dependent interval variables
  bool _get_depbnd
    ( const bool reset=true );

  //! @brief Set branching scores for current node
  void _set_branching_scores
    ( std::map<unsigned,double>&scores );
  //! @brief Attribute scores based on block decomposition
  void _set_block_scores
    ( std::map<unsigned,double>&scores );
  //! @brief Set up branching scores based on Chebyshev model enclosures
  void _set_cheb_scores
    ( std::map<unsigned,double>&scores );
  //! @brief Append branching scores based on Chebyshev model enclosures
  void _update_cheb_scores
    ( std::map<unsigned,double>&scores, const FFVar&FV,
      const SCVar<T>&SCV, const double weight, const bool linexcl=false ) const;

  //! @brief User-function to subproblems in SBB
  typename SBB<T>::STATUS subproblems
    ( const typename SBB<T>::TASK task, SBBNode<T>*node,
      std::vector<double>&p, double&f, const double INC )
    { return CSEARCH_BASE<T>::_subproblems( _NLPSLV, options, stats, task, node, p, f, INC ); }

  //! @brief User-function to subproblems in SBP
  virtual typename SBP<T>::STATUS assess
    ( SBPNode<T>*node )
    { throw Exceptions( Exceptions::INTERN ); }
    //{ return CSEARCH_BASE<T>::_assess( options, stats, node ); }

  //! @brief Initialize display
  void _display_pwl_init
    ();
  //! @brief Final display
  void _display_pwl_final
    ( const unsigned iter );

private:
  //! @brief Private methods to block default compiler methods
  NLGO
    (const NLGO&);
  NLGO& operator=
    (const NLGO&);
};

template <typename T> inline void
NLGO<T>::setup
( const T*P, const unsigned*tvar, std::ostream&os )
{
  _issetup = false;

  // full set of decision variables (independent & dependent)
  _nrvar = BASE_NLP::_var.size();
  _nrdep = BASE_NLP::_dep.size();
  _var = BASE_NLP::_var;
  _var.insert( _var.end(), BASE_NLP::_dep.begin(), BASE_NLP::_dep.end() );
  _nvar = _var.size();

  // objective function
  _obj = BASE_NLP::_obj;
  if( std::get<0>(_obj).size() > 1 ) throw Exceptions( Exceptions::MOBJ );

  // full set of constraints (main & equation system)
  _ctr = BASE_NLP::_ctr;
  _eq.clear();
  for( unsigned ie=0; ie<_sys.size(); ie++ ){
    _eq.push_back( ( _sys[ie] ) );
    std::get<0>(_ctr).push_back( EQ );
    std::get<1>(_ctr).push_back( _sys[ie] );
    std::get<2>(_ctr).push_back( _sysm[ie] );
    std::get<3>(_ctr).push_back( false );
  }

  // NCO cuts as redundant constraints
  if( options.NCOCUTS ){
    set_nco( options, P, tvar );
    _var.push_back( std::get<2>(BASE_NLP::_obj)[0] );
    _var.insert( _var.end(), std::get<2>(BASE_NLP::_ctr).begin(), std::get<2>(BASE_NLP::_ctr).end() );
    _var.insert( _var.end(), _sysm.begin(), _sysm.end() );
    for( unsigned ivar=0; ivar<_nrvar; ivar++ ){
      if( tvar && tvar[ivar] ) continue;
      _var.push_back( _varlm[ivar] );
      _var.push_back( _varum[ivar] );
    }
    for( unsigned idep=0; idep<_nrdep; idep++ ){
      if( tvar && tvar[_nrvar+idep] ) continue;
      _var.push_back( _deplm[idep] );
      _var.push_back( _depum[idep] );
    }
    for( unsigned ic=0; ic<std::get<0>(BASE_NLP::_nco).size(); ++ic ){
      std::get<0>(_ctr).push_back( std::get<0>(BASE_NLP::_nco)[ic] );
      std::get<1>(_ctr).push_back( std::get<1>(BASE_NLP::_nco)[ic] );
      std::get<2>(_ctr).push_back( FFVar() );
      std::get<3>(_ctr).push_back( true );
    }
  }

  _nctr = std::get<0>(_ctr).size();
  _neq = _eq.size();
  _nvar = _var.size();
  _var_fperm.resize(_nvar);
  _var_rperm.resize(_nvar);
  for( unsigned i=0; i<_nvar; i++ )
    _var_fperm[i] = _var_rperm[i] = i;

  // internal setup, including structure detection
  CSEARCH_BASE<T>::_setup( options, _dag, os );

  // Identify linear variables
  FFDep fgdep = std::get<0>(_obj).size()? std::get<1>(_obj)[0].dep(): 0.;
  for( unsigned j=0; j<_nctr; j++ )
    fgdep += std::get<1>(_ctr)[j].dep();
#ifdef MC__NLGO_DEBUG
  std::cout << "DEPS <- " << fgdep << std::endl;
  //int dum; std::cin >> dum;
#endif
  _var_lin.clear();
  for( unsigned i=0; i<_nvar; i++ ){
    auto it = fgdep.dep().find( _var[i].id().second );
    if( it == fgdep.dep().end() || it->second == FFDep::L )
      _var_lin.insert(i);
  }

  if( options.DISPLAY ){
    os << "LINEAR VARIABLES IN FUNCTIONS:            " << _var_lin.size() << std::endl
       << "NONLINEAR VARIABLES IN FUNCTIONS:         " << _nvar-_var_lin.size() << std::endl;
  }

  _var_excl.clear();
  if( options.DEPSUSE ){
    // User-defined equation subsystem
    if( _neq ){
      _nrvar = _nvar - _nrdep; // Multipliers may not be dependents in this case
      _AEBND.set( *this );
      _AEBND.options = options.AEBND;
      _AEBND.setup();
    }
    // Bordered-block decomposition of equation subsystem
    else{
      BASE_OPT::t_CTR* Ctype = std::get<0>(_ctr).data();
      FFVar* Cvar = std::get<1>(_ctr).data();
      for( unsigned i=0; i<_nctr; i++ )
        if( Ctype[i] == BASE_OPT::EQ ) _eq.push_back( Cvar[i] );
      _neq = _eq.size();
      if( _neq ){
        std::vector<int> IP(_neq), IQ(_nvar), IPROF(_nvar), IFLAG(3);
        _dag->MC33( _neq, _eq.data(), _nvar, _var.data(), IP.data(),
                    IQ.data(), IPROF.data(), IFLAG.data(), options.DISPLAY>1?true:false );
        // Identify linear variables in dependent blocks and exclude from branching
        _nrvar = IFLAG[2];
        _nrdep = IFLAG[1]-_nrvar;
        _AEBND.set_dag( _dag );
        _AEBND.reset_dep();
        _AEBND.reset_sys();
        for( unsigned i=0; i<_nrdep; ++i ){
          _AEBND.add_dep( _var[IQ[i]-1] );
          _AEBND.add_sys( _eq[IP[i]-1] );
          //std::cout << _eq[IP[i]-1] << std::endl;
          _var_rperm[_nrvar+i] = IQ[i]-1;
          _var_fperm[IQ[i]-1] = _nrvar+i;
        }
        _AEBND.reset_var();
        for( int i=_nrdep; i<IFLAG[1]; ++i ){
          _AEBND.add_var( _var[IQ[i]-1] );
          _var_rperm[i-_nrdep] = IQ[i]-1;
          _var_fperm[IQ[i]-1] = i-_nrdep;
        }
      }
    }
    _AEBND.options = options.AEBND;
    _AEBND.setup();

    // Count linear variables in dependent blocks
    unsigned ndeplin = 0; 
    for( unsigned ib=0; ib<_AEBND.noblk(); ib++ ){
      if( !_AEBND.linblk(ib) ) continue;
      for( unsigned j=0; j<_AEBND.nblk(ib); j++ ){ // variables in block
        auto jdep = _AEBND.depblk(ib)[j].dep().dep();
        assert( jdep.size() == 1 );
        if( !_AEBND.lindepblk(ib,j) ) continue;
        //_var_excl.insert( jdep.begin()->first );
        ndeplin++;
      }
      //if( !_AEBND.linblk(ib) ) continue;
      //for( unsigned j=0; j<_AEBND.nblk(ib); j++ ){ // variables in block
      //  auto jdep = _AEBND.depblk(ib)[j].dep().dep();
      //  assert( jdep.size() == 1 );
      //  _var_excl.insert( jdep.begin()->first );
      //  ndeplin++;
      //}
    }
    if( options.DISPLAY ){
      os << "LINEAR VARIABLES IN DEPENDENT BLOCKS:     " << ndeplin << std::endl
         << "TOTAL VARIABLES IN DEPENDENT BLOCKS:      " << _AEBND.dep().size() << std::endl;
#ifdef MC__NLGO_DEBUG
      { int dum; std::cout << "PAUSED--"; std::cin >> dum; }
#endif
    }
  }
  // No equation subsystem considered
  else{
    _nrvar = _nvar;
    _nrdep = 0;
    _AEBND.reset_dep();
    _AEBND.reset_sys();
  }

  _issetup = true;
  return;
}

template <typename T>
inline void
NLGO<T>::_set_SLVLOC
()
{
  _NLPSLV->options = options.NLPSLV;
  _NLPSLV->set( *this );
  _NLPSLV->setup();
}

template <typename T>
inline const double*
NLGO<T>::_get_SLVLOC
( const SOLUTION_NLP&locopt )
{
  //std::cout << locopt.n << " =?= " << BASE_NLP::_var.size() << " + " << BASE_NLP::_dep.size() << std::endl; 
  assert( locopt.n == (int)BASE_NLP::_var.size() + (int)BASE_NLP::_dep.size() );
  assert( locopt.m == (int)std::get<0>(BASE_NLP::_ctr).size() + (int)BASE_NLP::_sys.size() );
  _Dvar.assign( locopt.p, locopt.p+locopt.n );
  if( options.NCOCUTS ){
    //_Dvar.insert( _Dvar.end(), _nvar-locopt.n, 0. ); // <- MODIFY TO USE ACTUAL FJ MULTIPLIER VALUES (AFTER RESCALING)?
    _Dvar.push_back( 1. );
    double mult_ineq = 1., mult_eq = 0.;
    unsigned ig = 0;
    for( auto&& t_ctr : std::get<0>(BASE_NLP::_ctr) ){
      switch( t_ctr ){
        case BASE_OPT::EQ: 
          _Dvar.push_back( locopt.ug[ig] );
          mult_eq += mc::sqr( locopt.ug[ig] );
          break;
        case BASE_OPT::LE:
        case BASE_OPT::GE:
          _Dvar.push_back( std::fabs( locopt.ug[ig] ) );
          mult_ineq += std::fabs( locopt.ug[ig] );
          break;
      }
      ig++;
    }
    for( auto it=_sys.begin(); it!=_sys.end(); ++it ){
      _Dvar.push_back( locopt.ug[ig] );
      mult_eq += sqr( locopt.ug[ig] );
      ig++;
    }
    for( int ivar=0; ivar<locopt.n; ivar++ ){
      if( !_tvar.empty() && _tvar[ivar] ) continue;
      _Dvar.push_back( std::fabs( locopt.upL[ivar] ) );
      _Dvar.push_back( std::fabs( locopt.upU[ivar] ) );
      mult_ineq += std::fabs( locopt.upL[ivar] ) + std::fabs( locopt.upU[ivar] );
    }
    double mult_scal = ( - mult_ineq + std::sqrt( mult_ineq * mult_ineq + 4. * mult_eq ) )
                      / ( 2. * mult_eq );
    for( unsigned i=locopt.n; i<_nvar; i++ ){
      _Dvar[i] *= mult_scal;
      //std::cout << "Dvar[" << i << "] = " << _Dvar[i] << std::endl;
    }
    //for( int i=0; i<locopt.m; i++ )
    //  std::cout << "ug[" << i << "] = " << locopt.ug[i] << std::endl;
    //for( int i=0; i<locopt.n; i++ ){
    //  std::cout << "upL[" << i << "] = " << locopt.upL[i] << std::endl;
    //  std::cout << "upU[" << i << "] = " << locopt.upU[i] << std::endl;
    //}
    //std::cout << "mult_eq = " << mult_eq << std::endl;
    //std::cout << "mult_ineq = " << mult_ineq << std::endl;
    //std::cout << "mult_scal = " << mult_scal << std::endl;
  }
  return _Dvar.data();
}

template <typename T>
inline bool
NLGO<T>::_init
( const T*P, const unsigned*tvar, const double*p0 )
{
  const unsigned n0var = BASE_NLP::_var.size(), n0dep = BASE_NLP::_dep.size();
  // Form vector of parameter and state bounds
  _Ivar.assign( P, P+n0var+n0dep );
  // Extra multiplier variables for NCO cuts
  if( options.NCOCUTS ){
    //_Ivar.insert( _Ivar.end(), _nvar-n0, T(0.,1.) ); // multipliers are bounded in [0,1]
    _Ivar.push_back( T(0.,1.) );
    for( auto&& t_ctr : std::get<0>(BASE_NLP::_ctr) )
      switch( t_ctr ){
        case BASE_OPT::EQ: _Ivar.push_back( T(-1.,1.) ); break;
        case BASE_OPT::LE:
        case BASE_OPT::GE: _Ivar.push_back( T(0.,1.) ); break;
      }
    _Ivar.insert( _Ivar.end(), _sysm.size(), T(-1.,1.) );
    for( unsigned ivar=0; ivar<n0var; ivar++ ){
      if( tvar && tvar[ivar] ) continue;
      _Ivar.push_back( T(0.,1.) );
      _Ivar.push_back( T(0.,1.) );
    }
    for( unsigned idep=0; idep<n0dep; idep++ ){
      if( tvar && tvar[_nrvar+idep] ) continue;
      _Ivar.push_back( T(0.,1.) );
      _Ivar.push_back( T(0.,1.) );
    }
  }
  // Form vector of parameter and state types (continuous vs discrete)
  _tvar.clear();
  if( tvar ){
    _tvar.assign( tvar, tvar+n0var+n0dep );
    if( options.NCOCUTS )
      _tvar.insert( _tvar.end(), _nvar-n0var-n0dep, false ); // multipliers are continuous variables
  }
  // Form vector of parameter and state values
  _Dvar.clear();
  if( p0 ){
    _Dvar.assign( p0, p0+n0var+n0dep );
    if( options.NCOCUTS )
      _Dvar.insert( _Dvar.end(), _nvar-n0var-n0dep, 0. );
  }
  return true;
}

template <typename T>
inline typename LPRELAX_BASE<T>::LP_STATUS
NLGO<T>::relax
( const T*P, const unsigned*tvar, const unsigned refine,
  const bool reset, const bool feastest )
{
  if( !_issetup || !_init( P, tvar ) ) throw Exceptions( Exceptions::SETUP );
  _ignore_deps = false;
  return CSEARCH_BASE<T>::_relax( options, stats, _Ivar.data(), tvar?_tvar.data():0,
                                  refine, reset, feastest );
}

template <typename T>
inline typename LPRELAX_BASE<T>::LP_STATUS
NLGO<T>::contract
( T*P, unsigned&nred, const unsigned*tvar, const double*inc,
  const bool reset, const bool feastest )
{
  if( !_issetup || !_init( P, tvar ) ) throw Exceptions( Exceptions::SETUP );
  _ignore_deps = false;
  auto status = CSEARCH_BASE<T>::_contract( options, stats, _Ivar.data(), nred,
                                            tvar?_tvar.data():0, inc, reset, feastest );
  auto it=_Ivar.begin();
  for( unsigned i=0; i<_nrvar+_nrdep && it!=_Ivar.end(); ++it, i++ )
    P[i] = *it;
#if defined (MC__NLGO_SHOW_CONTRACT)
  std::cout << "Contracted bounds:" << std::endl;
  for( auto it=_Ivar.begin(); it!=_Ivar.end(); ++it )
    std::cout << " " << *it << std::endl;
#endif
  return status;
}

template <typename T>
inline int
NLGO<T>::solve
( const T*P, const unsigned*tvar, const double*p0, std::ostream&os )
{
  if( !_issetup || !_init( P, tvar, p0 ) ) throw Exceptions( Exceptions::SETUP );

  switch( options.CSALGO ){
   case Options::BB: default:
     return CSEARCH_BASE<T>::_solve_sbb( _NLPSLV, options, stats, _Ivar.data(),
                                         tvar?_tvar.data():0, p0?_Dvar.data():0, os );
   case Options::PWL:
     return CSEARCH_BASE<T>::_solve_pwl( _NLPSLV, options, stats, _Ivar.data(),
                                         tvar?_tvar.data():0, p0?_Dvar.data():0, os );
  }
}

template <typename T>
inline bool
NLGO<T>::_get_depcheb
( const bool reset )
{
#if defined (MC__NLGO_DEBUG_CHEBDEPS)
  for( unsigned i=0; i<_nrvar; i++ )
    std::cout << "CMrvar[" << i << "] =" << _CMrvar[i];
  for( unsigned i=0; i<_nrdep; i++ )
    std::cout << "CMrdep[" << i << "] =" << _CMrdep[i];
#endif
  // Reset list of dependents
  if( reset ) _CMndxdep.clear();

  std::vector<SCVar<T>> CMrdep0 = _CMrdep; // local copy for convergence test later on
  auto status = t_AEBND::NORMAL;

  // Warmstart dependent subsystem solutions
  if( options.CMODWARMS ){
    try{
      // Compute polynomial approximant
      _AEBND.options.PMNOREM = true;
      std::vector<double> ddep( _nrdep );
      for( unsigned i=0; i<_nrdep; i++ )
        ddep[i] = Op<SCVar<T>>::diam(_CMrdep[i]);

      std::vector<SCVar<T>> CMrdep1 = _CMrdep;
      if( _AEBND.solve( _CMrvar.data(), CMrdep1.data(), CMrdep0.data() )
       == t_AEBND::NORMAL ){
        // Compute remainder bounds for the dependents
        std::vector<double> cobj( _nrvar+1 );
        std::vector<PolVar<T>> vobj( _nrvar+1 );

        for( unsigned i=0; i<_nrdep; i++ ){
#if defined (MC__NLGO_DEBUG_CHEBDEPS)
          std::cout << "_CMrdep1[" << i << "] =" << CMrdep1[i];
#endif
          // Set-up contaction objective function
          for( unsigned j=0; j<_nrvar; j++ ){
            cobj[j] = -CMrdep1[i].linear(j);
            vobj[j] = _POLvar[_var_rperm[j]];
          }
          cobj[_nrvar] = 1.;
          vobj[_nrvar] = _POLvar[_var_rperm[_nrvar+i]];

          // Solve for lower bound
          _set_LPrelax( _nrvar+1, vobj.data(), cobj.data(), MIN );
          _solve_LPmodel( options, stats, _var );
          double RxL = get_objective();
#if defined (MC__NLGO_DEBUG_CHEBDEPS)
          std::cout << "  RxL = " << RxL << std::endl;
#endif

          // Solve for upper bound
          _set_LPrelax( _nrvar+1, vobj.data(), cobj.data(), MAX );
          _solve_LPmodel( options, stats, _var );
          double RxU = get_objective();
#if defined (MC__NLGO_DEBUG_CHEBDEPS)
          std::cout << "  RxU = " << RxU << std::endl;
#endif

          // Setup first-order polynomial model
#if defined (MC__NLGO_DEBUG_CHEBDEPS)
          std::cout << "  diam = " << ddep[i] << std::endl;
#endif
          //if( !( RxU - RxL < 2.*ddep[i] ) ) throw(0);
          if( !( RxU - RxL < ddep[i] ) ) continue;
          CMrdep1[i] = T( RxL, RxU );
          for( unsigned j=0; j<_nrvar; j++ )
            CMrdep1[i] += CMrdep1[i].linear(j) * _CMrvar[j];
        }
        // Warm-started solution
        _CMrdep = CMrdep1;
      }
    }

    // Warm-start unsuccessful?
    catch(...){
#ifdef MC__CSEARCH_SHOW_DEPS
    std::cout << "Warm-start: Failed" << std::endl;
#endif
      status = t_AEBND::FAILURE;
    }
  }

  // Bound dependent subsystem solutions
  try{
    status = _AEBND.solve( _CMrvar.data(), _CMrdep.data(), _CMrdep.data() );
    //if( status != t_AEBND::NORMAL ) throw(0);
  }

  // Bounding unsuccessful?
  catch(...){
#ifdef MC__CSEARCH_SHOW_DEPS
    std::cout << "Dependents: Failed" << std::endl;
#endif
    //for( unsigned ib=_AEBND.iblk(); ib<_AEBND.noblk(); ++ib )
    //  for( unsigned j=0; j<_AEBND.nblk(ib); j++ )
    //    _CMrdep[_AEBND.rpdep(ib,j)] = T(-_AEBND.INF,_AEBND.INF);
    status = t_AEBND::FAILURE;
  }

  // Dependent system is either structurally singular or infeasible
  //if( status == t_AEBND::EMPTY )
  //  return false;
  //else
  if( status != t_AEBND::NORMAL && status != t_AEBND::FAILURE )
    return true;

#ifdef MC__CSEARCH_DEBUG_BNDDEPS
  for( unsigned i=0; i<_nrdep; i++ )
    std::cout << "CMrdep[" << i << "] =" << _CMrdep[i];
  {int dum; std::cout << "PAUSED"; std::cin >> dum;}
#endif

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
        //_CMrdep[_AEBND.rpdep(ib,j)].simplify();
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
  return true;
}
/*
  // Append to set of dependents if enough improvement
#ifdef MC__NLGO_SHOW_DEPS
  std::cout << "Dependents: {";
#endif
  if( options.CMODJOINT ){
    bool allred = true;
    for( unsigned i=0; allred && i<_nrdep; i++ ){
      if( Op<T>::diam( _CMrdep[i].R() + _CMrdep[i].bndord( options.CMODPROP+1 ) )
          < options.CMODATOL + options.CMODRTOL * Op<T>::diam( _CMrdep[i].B() ) )
        continue;
      allred = false;
    }
    for( unsigned i=0; allred && i<_nrdep; i++ ){
      _ndxdep.insert( _var_rperm[_nrvar+i] );
#ifdef MC__NLGO_SHOW_DEPS
      std::cout << " " << _var_rperm[_nrvar+i];
#endif
    }
  }
  else{
    for( unsigned i=0; i<_nrdep; i++ ){
      if( Op<T>::diam( _CMrdep[i].R() + _CMrdep[i].bndord( options.CMODPROP+1 ) )
          < options.CMODATOL + options.CMODRTOL * Op<T>::diam( _CMrdep[i].B() ) ){
        _ndxdep.insert( _var_rperm[_nrvar+i] );
#ifdef MC__NLGO_SHOW_DEPS
        std::cout << " " << _var_rperm[_nrvar+i];
#endif
      }
    }
  }
#ifdef MC__NLGO_SHOW_DEPS
  std::cout << " }" << std::endl;
#endif
}
*/
template <typename T>
inline bool
NLGO<T>::_get_depbnd
( const bool reset )
{
#ifdef MC__CSEARCH_SHOW_DEPS
  for( unsigned i=0; i<_nrdep; i++ )
    std::cout << "Irdep[" << i << "] =" << _Irdep[i];
#endif

  // Reset list of dependents
  if( reset) _Indxdep.clear();

  std::vector<T> Irdep0 = _Irdep; // local copy for convergence test later on
  auto status = t_AEBND::NORMAL;

  // Bound dependent subsystem solutions
  try{
    status = _AEBND.solve( _Irvar.data(), _Irdep.data(), Irdep0.data() );
    //if( status != t_AEBND::NORMAL ) throw(0);
  }

  // Bounding unsuccessful?
  catch(...){
#ifdef MC__CSEARCH_SHOW_DEPS
    std::cout << "Dependents: Failed" << std::endl;
#endif
    //for( unsigned ib=_AEBND.iblk(); ib<_AEBND.noblk(); ++ib )
    //  for( unsigned j=0; j<_AEBND.nblk(ib); j++ )
    //    _Irdep[_AEBND.rpdep(ib,j)] = T(-_AEBND.INF,_AEBND.INF);
    status = t_AEBND::FAILURE;
  }

  // Dependent system is either structurally singular or infeasible
  if( status == t_AEBND::EMPTY )
    return false;
  else if( status != t_AEBND::NORMAL && status != t_AEBND::FAILURE )
    return true;

#ifdef MC__CSEARCH_DEBUG_BNDDEPS
  for( unsigned i=0; i<_nrdep; i++ )
    std::cout << "Irdep[" << i << "] =" << _Irdep[i] << std::endl;
  {int dum; std::cout << "PAUSED"; std::cin >> dum;}
#endif

  // Append dependent showing improvement to _Indxdep.
  for( unsigned ib=0; ib<_AEBND.iblk(); ++ib )
    if( _AEBND.uniblk(ib) ) 
      for( unsigned j=0; j<_AEBND.nblk(ib); j++ )
        _Indxdep.insert( _var_rperm[_nrvar+_AEBND.rpdep(ib,j)] );
#ifdef MC__CSEARCH_SHOW_DEPS
  std::cout << "_Indxdep = { ";
  for( auto&&ndx : _Indxdep ) std::cout << ndx << " ";
  std::cout << " }" << std::endl;
#endif

//#ifdef MC__CSEARCH_SHOW_DEPS
//  std::cout << "_Indxdep = {";
//#endif
//  for( unsigned ib=0; ib<_AEBND.iblk(); ++ib )
//    for( unsigned j=0; j<_AEBND.nblk(ib); j++ )
//      if( Op<T>::diam( _Irdep[_AEBND.rpdep(ib,j)] ) < Op<T>::diam( Irdep0[_AEBND.rpdep(ib,j)] ) ){
//        _Indxdep.insert( _var_rperm[_nrvar+_AEBND.rpdep(ib,j)] );
//#ifdef MC__CSEARCH_SHOW_DEPS
//        std::cout << " " << _var[_var_rperm[_nrvar+_AEBND.rpdep(ib,j)]];
//               // << "(" << _AEBND.rpdep(ib,j) << ") " << _Irdep[_AEBND.rpdep(ib,j)];
//#endif
//      }
//#ifdef MC__CSEARCH_SHOW_DEPS
//  std::cout << " }" << std::endl;
//#endif
  return true;
}

template <typename T>
inline void
NLGO<T>::_set_branching_scores
( std::map<unsigned,double>&scores )
{
  switch( options.SCOBCHMETH ){
    case Options::NOSCO:
      return;
    case Options::BLKSCO:
      return _set_block_scores( scores );
    case Options::CHEBSCO:
      return _set_cheb_scores( scores );
  }
}

template <typename T>
inline void
NLGO<T>::_set_block_scores
( std::map<unsigned,double>&scores )
{
  scores.clear();
  if( !options.DEPSUSE ) return;
#ifdef MC__CSEARCH_SHOW_SCORES
  std::cout << "_ndxdep = { ";
  for( auto&& ndx: _ndxdep ) std::cout << ndx << " ";
  std::cout << "}" << std::endl;
#endif
  std::set<unsigned> highdep;

  // Append scores based on block structure
  for( unsigned ib=0, isys=0;  ib<_AEBND.noblk(); ++ib ){
    const double score = 1.-(double)ib/_AEBND.noblk();
    // Check uniqueness of solution branch in current block
    bool uniblk = true;
    for( unsigned j=0; uniblk && j<_AEBND.nblk(ib); j++ ){
      const unsigned ivar = _var_rperm[_nrvar+_AEBND.rpdep(ib,j)];
      if( _ndxdep.find( ivar ) == _ndxdep.end() ) uniblk = false;
    }
    // Assign score to block independent when uniqueness is not guaranteed, 
    // or zero score if uniqueness is guaranteed
    for( unsigned j=0; j<_AEBND.nblk(ib); j++, isys++ ){
      auto sys = _AEBND.sys()[_nrdep-isys-1]; // current equation
      for( unsigned i=0; i<_nrvar; i++ ){
        const unsigned ivar = _var_rperm[i];
        // no score if independent ivar was previously scored;
        if( scores.find( ivar ) != scores.end() ) continue;
        // no score if independent ivar does not participate in equation
        if( !sys.dep().dep(ivar).first ) continue;
        // score if uniqueness is not guaranteed
        if( !uniblk ) scores.insert( std::make_pair( ivar, score ) );
        // otherwise, add ivar to higher dependency list 
        else highdep.insert( ivar );
      }
    }
    // Assign scores to independents from previous block 
    if( !uniblk ){
      for( auto&& ivar: highdep )
        scores.insert( std::make_pair( ivar, score ) );
      highdep.clear();
    }
    // Assign scores to block dependents when uniqueness is not guaranteed
    // and block is nonlinear
    for( unsigned j=0; !uniblk && !_AEBND.linblk(ib) && j<_AEBND.nblk(ib); j++ ){
      const unsigned ivar = _var_rperm[_nrvar+_AEBND.rpdep(ib,j)];
      scores.insert( std::make_pair( ivar, score ) );
    }
  }

#ifdef MC__CSEARCH_SHOW_SCORES
  for( auto it=scores.cbegin(); it!=scores.cend(); ++it )
    std::cout << std::right << std::setw(4) << it->first << " " << it->second << std::endl; 
  { int dum; std::cout << "PAUSED"; std::cin >> dum; }
#endif

  // Assign zero scores to all remaining independent variables
  unsigned ivar = 0;
  for( auto&& scovar: scores )
    for( ; ivar<scovar.first; ivar++ )
      scores.insert( std::make_pair( ivar, 0. ) );
  for( ivar++; ivar<_nvar; ivar++ )
  //for( auto&& ivar: highdep )
    scores.insert( std::make_pair( ivar, 0. ) );
}

//template <typename T>
//inline void
//NLGO<T>::_set_block_scores
//( std::map<unsigned,double>&scores )
//{
//  scores.clear();
//  if( !options.DEPSUSE ) return;
//#ifdef MC__CSEARCH_SHOW_SCORES
//  std::cout << "_ndxdep = { ";
//  for( auto&& ndx : _ndxdep ) std::cout << ndx << " ";
//  std::cout << "}" << std::endl;
//#endif

//  // Append scores based on block structure
//  for( unsigned ib=0, isys=0;  ib<_AEBND.noblk(); ++ib ){
//    const double score = 1.-(double)ib/_AEBND.noblk();
//    bool is_contracting = true;
//    // Assign scores to block dependents
//    for( unsigned j=0; j<_AEBND.nblk(ib); j++ ){
//      const unsigned ivar = _var_rperm[_nrvar+_AEBND.rpdep(ib,j)];
//      scores.insert( std::make_pair( ivar, score ) );
//      if( _ndxdep.find( ivar ) == _ndxdep.end() ) is_contracting = false;
//    }
//    if( is_contracting ) continue;
//    // Assign scores to block independent *only if* block is not yet contracting
//    for( unsigned j=0; j<_AEBND.nblk(ib); j++, isys++ ){
//      auto sys = _AEBND.sys()[_nrdep-isys-1];
//      for( unsigned i=0; i<_nrvar; i++ ){
//        const unsigned ivar = _var_rperm[i];
//        if( scores.find( ivar ) != scores.end()
//         || !sys.dep().dep(ivar).first ) continue;
//        scores.insert( std::make_pair( ivar, score ) );
//      }
//    }
//  }
//  // Assign zero scores to all remaining independent variables
//  for( unsigned j=0; j<_nrvar; j++ ){
//    const unsigned ivar = _var_rperm[j];
//    if( scores.find( ivar ) != scores.end() ) continue;
//    scores.insert( std::make_pair( ivar, 0. ) );
//  }

//#ifdef MC__CSEARCH_SHOW_SCORES
//  for( auto it=scores.cbegin(); it!=scores.cend(); ++it )
//    std::cout << std::right << std::setw(4) << it->first << " " << it->second << std::endl; 
//  { int dum; std::cout << "PAUSED"; std::cin >> dum; }
//#endif
//}

//template <typename T>
//inline void
//NLGO<T>::_set_block_scores
//( std::map<unsigned,double>&scores )
//{
//  scores.clear();
//  if( !options.DEPSUSE ) return;
//#ifdef MC__CSEARCH_SHOW_SCORES
//  std::cout << "_ndxdep = { ";
//  for( auto&& ndx : _ndxdep ) std::cout << ndx << " ";
//  std::cout << "}" << std::endl;
//#endif

//  // Append scores based on block structure
//  for( unsigned ib=0, isys=0;  ib<_AEBND.noblk(); ++ib ){
//    const double score = 1.-(double)ib/_AEBND.noblk();
//    // Assign scores to block independent *only if* block is not yet contracting
//    if( _AEBND.uniblk(ib) ) continue;
//    // Assign scores to block dependents
//    for( unsigned j=0; j<_AEBND.nblk(ib); j++ ){
//      const unsigned ivar = _var_rperm[_nrvar+_AEBND.rpdep(ib,j)];
//      scores.insert( std::make_pair( ivar, score ) );
//    }
//    // Assign scores to block independents
//    for( unsigned j=0; j<_AEBND.nblk(ib); j++, isys++ ){
//      auto sys = _AEBND.sys()[_nrdep-isys-1];
//      for( unsigned i=0; i<_nrvar; i++ ){
//        const unsigned ivar = _var_rperm[i];
//        if( scores.find( ivar ) != scores.end()
//         || !sys.dep().dep(ivar).first ) continue;
//        scores.insert( std::make_pair( ivar, score ) );
//      }
//    }
//  }
//  // Assign zero scores to all remaining independent variables
//  for( unsigned j=0; j<_nrvar; j++ ){
//    const unsigned ivar = _var_rperm[j];
//    if( scores.find( ivar ) != scores.end() ) continue;
//    scores.insert( std::make_pair( ivar, 0. ) );
//  }

//#ifdef MC__CSEARCH_SHOW_SCORES
//  for( auto it=scores.cbegin(); it!=scores.cend(); ++it )
//    std::cout << std::right << std::setw(4) << it->first << " " << it->second << std::endl; 
//  { int dum; std::cout << "PAUSED"; std::cin >> dum; }
//#endif
//}

template <typename T>
inline void
NLGO<T>::_set_cheb_scores
( std::map<unsigned,double>&scores )
{
  scores.clear();
  if( options.RELMETH != Options::CHEB && options.RELMETH != Options::HYBRID ) return;

  // Append scores based on Chebyshev model for (nonlinear) objective
  auto ito = std::get<0>(_obj).begin();
  for( unsigned j=0; ito!=std::get<0>(_obj).end() && j<1; ++ito, j++ ){
    if( _fct_lin.find(std::get<1>(_obj).data()) != _fct_lin.end() ) continue; // only nonlinear objectives
    _update_cheb_scores( scores, std::get<1>(_obj)[0], _CMobj, 1. );
  }

  // Append scores based on Chebyshev model for (nonlinear) constraints
  auto itc = std::get<0>(_ctr).begin();
  for( unsigned j=0; itc!=std::get<0>(_ctr).end(); ++itc, j++ ){
    if( _fct_lin.find(std::get<1>(_ctr).data()+j) != _fct_lin.end()
     || std::get<3>(_ctr)[j] ) continue;
    switch( (*itc) ){ // only score equality and active inequality constraints
      case LE: if( Op<SCVar<T>>::u(_CMctr[j]) < 0. ) continue;
               break;
      case GE: if( Op<SCVar<T>>::l(_CMctr[j]) > 0. ) continue;
               break;
      default: break;
    }
    _update_cheb_scores( scores, std::get<1>(_ctr)[j], _CMctr[j], 1. );
  }

#ifdef MC__CSEARCH_SHOW_SCORES
  for( auto it=scores.cbegin(); it!=scores.cend(); ++it )
    std::cout << std::right << std::setw(4) << it->first << " " << it->second << std::endl; 
  { int dum; std::cout << "PAUSED"; std::cin >> dum; }
#endif
}

template <typename T>
inline void
NLGO<T>::_update_cheb_scores
( std::map<unsigned,double>&scores, const FFVar&FV,
  const SCVar<T>&SCV, const double weight, const bool linexcl )
const
{
  // Case Chebyshev model bounding was not successful
  if( Op<T>::diam(SCV.R()) >= SBB<T>::INF ){
    double contrib = Op<T>::diam(SCV.R()) * std::fabs( weight );
    for( unsigned i=0; i<_nvar; i++ ){
      auto ft = FV.dep().dep().find( _var[i].id().second );
      // Discard variables that participate linearly or polynomially (no bound explosion possible)
      if( ft == FV.dep().dep().end() || ft->second == FFDep::L || ft->second == FFDep::P ) continue;
      auto its = scores.insert( std::make_pair( i, contrib ) );
      if( !its.second && its.first->second < contrib ) its.first->second = contrib;
    }
    return;
  }

  // Case Chebyshev model bounding was successful
  std::map<unsigned,double> newscores;
  // Scores through remainder
  double contrib = Op<T>::diam(SCV.R()) * std::fabs( weight );
  unsigned remdim = 0;
  for( unsigned i=0; i<_nvar; i++ ){
    auto ft = FV.dep().dep().find( _var[i].id().second );
    if( ft == FV.dep().dep().end() || ft->second == FFDep::L ) continue;
    auto its = newscores.insert( std::make_pair( i, contrib ) );
    assert( its.second );
    remdim++;
  }
  for( auto&& score : newscores ) score.second /= remdim;

  // Scores through order 1+ (linexcl=false) or 2+ (linexcl=true) coefficients
  auto kt = SCV.coefmon().lower_bound( std::make_pair( 2 ,std::map<unsigned,unsigned>() ) );
  //auto kt = SCV.coefmon().lower_bound( std::make_pair( linexcl?2:1 ,std::map<unsigned,unsigned>() ) );
  for( ; kt != SCV.coefmon().cend(); ++kt ){
    // Scale with the number of participating variables
    contrib = std::fabs( weight * kt->second ) / kt->first.first;
    for( auto ie=kt->first.second.begin(); ie!=kt->first.second.end(); ++ie ){
      auto its = newscores.insert( std::make_pair( ie->first, contrib ) );
      if( !its.second ) its.first->second += contrib;
    }
  }

  // Scale scores?
  //double bnd = Op<T>::diam(SCV.B());
  //for( auto&& score : newscores ) score.second /= bnd;

  // Update scores with the new scores (in newscores)
  for( auto&& score : newscores ){
    auto its = scores.insert( score );
    if( !its.second && its.first->second < score.second ) its.first->second = score.second;
  }
}

template <typename T> inline void
NLGO<T>::_display_pwl_init
()
{
  _odisp.str("");
  if( options.DISPLAY <= 1 ) return;
  _odisp << std::right
  	 << std::setw(_IPREC) << "ITER"
  	 << std::setw(_IPREC) << "PREPROC"
  	 << std::setw(_DPREC+8) << "RELAX   "
  	 << std::setw(_DPREC+8) << "BEST   "
  	 << std::setw(_DPREC+8) << "CPU TOT  "
  	 << std::setw(_DPREC+8) << "STATUS"
  	 << std::endl;  
}

template <typename T>
inline void
NLGO<T>::_display_pwl_final
( const unsigned iter )
{
  stats.tALL += cpuclock();
  if( options.DISPLAY <= 0 ) return;
  _odisp << std::endl << "#  TERMINATION: ";
  _odisp << std::fixed << std::setprecision(6) << stats.tALL << " CPU SEC"
         //<< "  (LBD:" << std::fixed << std::setprecision(1)
         //<< _tLBD/(_tcur-stats.tstart)*1e2 << "%  UBD:" << std::fixed
         //<< std::setprecision(1) << _tUBD/(_tcur-stats.tstart)*1e2 << "%)"
         << std::endl;

  // No feasible solution found
  if( _p_inc.empty() )
    _odisp << "#  NO FEASIBLE SOLUTION FOUND" << std::endl;

  // Feasible solution found
  else{
    // Incumbent
    _odisp << "#  INCUMBENT VALUE:" << std::scientific
           << std::setprecision(_DPREC) << std::setw(_DPREC+8) << _f_inc
           << std::endl;
    _odisp << "#  INCUMBENT POINT:";
    for( unsigned ip=0, id=0; ip<_nvar; ip++, id++ ){
      if( id == _LDISP ){
        _odisp << std::endl << std::left << std::setw(19) << "#";
        id = 0;
      }
      _odisp << std::right << std::setw(_DPREC+8) << _p_inc[ip];
    }
    _odisp << std::endl;
  }

  _odisp << "#  NUMBER OF RELAXATION REFINEMENTS:   " << iter << std::endl;
  if( get_status() != LPRELAX_BASE<T>::LP_OPTIMAL )
    _odisp << "#  SOLUTION INTERRUPTED AFTER RELAXATION FAILURE" << std::endl;
  else{
    _odisp << "#  OPTIMALITY GAP:   "
           << std::fabs( _f_inc - get_objective() ) << " (ABS)" << std::endl
           << "                     "
           << 2. * std::fabs( _f_inc - get_objective() )
                 / std::fabs( _NLPSLV->solution().f + get_objective() ) << " (REL)" << std::endl;
  }
}

template <typename T>
inline void
NLGO<T>::Options::display
( std::ostream&out ) const
{
  // Display NLGO Options
  out << std::left;
  out << std::setw(60) << "  COMPLETE SEARCH METHOD";
  switch( CSALGO ){
   case BB: out << "BB" << std::endl; break;
   case PWL: out << "PWL" << std::endl; break;
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
  //out << std::setw(60) << "  USE OF REDUCED-SPACE STRATEGY";
  //switch( CMODRED ){
  // case NOREDUC:    out << "NOREDUC" << std::endl; break;
  // case APPEND:     out << "APPEND"  << std::endl; break;
  //}
  if( CMODRED == APPEND ){
    out << std::setw(60) << "  ORDER OF CHEBYSHEV-REDUCTION CONSTRAINTS";
    switch( CMODDEPS ){
     case 0:  out << "-\n"; break;
     default: out << CMODDEPS << std::endl;
              out << std::setw(60) << "  CONVERGENCE THRESHOLD FOR REDUCED-SPACE CHEBYSHEV MODEL"
                  << std::scientific << std::setprecision(1)
                  << CMODRTOL << std::endl;
    }
  }
  out << std::setw(60) << "  APPEND NCO CUTS"
      << (NCOCUTS?"Y\n":"N\n");
  if( CMODCUTS ){
    out << std::setw(60) << "  METHOD FOR NCO CUTS";
    switch( NCOMETH ){
     case FSA: out << "FSA\n";
     case ASA: out << "ASA\n";
    }
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

template <typename T>
inline std::ostream&
operator <<
( std::ostream&out, const NLGO<T>&NLP )
{
  out << std::right << std::endl
      << std::setfill('_') << std::setw(72) << " " << std::endl << std::endl << std::setfill(' ')
      << std::setw(55) << "NONLINEAR GLOBAL OPTIMIZATION IN CRONOS\n"
      << std::setfill('_') << std::setw(72) << " " << std::endl << std::endl << std::setfill(' ');

  // Display NLGO Options
  //out << std::left << "SPATIAL BRANCH-AND-BOUND OPTIONS:\n\n";
  //NLP.SBB<T>::options.display( out );
  NLP.NLGO<T>::options.display( out );

  out << std::setfill('_') << std::setw(72) << " " << std::endl << std::endl << std::setfill(' ');
  return out;
}

} // end namescape mc

#endif
