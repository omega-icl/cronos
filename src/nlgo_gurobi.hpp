// Copyright (C) 2015 Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

/*!
\page page_NLGO_GUROBI Nonlinear Global Optimization using MC++ and GUROBI
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
where \f$f, g_1, \ldots, g_m\f$ are factorable, potentially nonlinear, real-valued functions; and \f$x_1, \ldots, x_n\f$ can be either continuous or integer decision variables. The class mc::NLGO_GUROBI solves such NLP or MINLP problems to global optimality using complete search. Two main methods are implemented in mc::NLGO_GUROBI:
- spatial branch-and-bound search
- hierarchy of semi-linear relaxations
.
In both methods, relaxations are generated for the nonlinear or nonconvex participating terms using various arithmetics in <A href="https://projects.coin-or.org/MCpp">MC++</A>.

\section sec_NLGO_GUROBI_setup How do I setup my optimization model?

Consider the following NLP:
\f{align*}
  \max_{\bf p}\ & p_1\,p_4\,(p_1+p_2+p_3)+p_3 \\
  \text{s.t.} \ & p_1\,p_2\,p_3\,p_4 \geq 25 \\
  & p_1^2+p_2^2+p_3^2+p_4^2 = 40 \\
  & 1 \leq p_1,p_2,p_3,p_4 \leq 5\,.
\f}

First, we define an mc::NLGO_GUROBI class as below:

\code
  mc::NLGO_GUROBI NLP;
\endcode

Next, we set the variables and objective/constraint functions by creating a direct acyclic graph (DAG) of the problem: 

\code
  #include "NLGO_GUROBI.hpp"
  
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

The variable bounds and types are passed to mc::NLGO_GUROBI in invoking the various methods, as described below.


\section sec_NLGO_GUROBI_methods What are the methods available?


Given initial bounds \f$P\f$ and initial guesses \f$p_0\f$ on the decision variables, the NLP model is solved using branch-and-bound search (default) as follows:

\code
  #include "interval.hpp"

  typedef mc::Interval I;
  I Ip[NP] = { I(1,5), I(1,5), I(1,5), I(1,5) };

  std::cout << NLP;
  int status = NLP.solve( Ip );
\endcode

The return value is the <A href="http://www.gurobi.com/documentation/6.0/refman/optimization_status_codes.html">GUROBI status</A>, which is of integer type. The following result is produced:

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
  NLP.options.CSALGO = mc::NLGO_GUROBI<I>::Options::PWL;
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

Other options can be modified to tailor the search, including output level, maximum number of iterations, tolerances, maximum CPU time, etc. These options can be modified through the public member mc::NLGO_GUROBI::options. 
*/

//TODO: 
//- [TO DO] Implement similar solver in CPLEX (maybe with a base class NLGO_base)
//- [TO DO] Enable KKT cuts and reduction constraints
//- [OK]    Implement SBB method
//- [OK]    Enable use of polynomial models in relaxations
//- [OK]    Enable multiple rounds of PWR refinement
//- [OK]    Support MINLP
//- [OK]    Exclude linear variables from bound contraction
//- [TO DO] Make it possible to add/remove a constraint from the model?

#ifndef MC__NLGO_GUROBI_HPP
#define MC__NLGO_GUROBI_HPP

#include <stdexcept>
#include <assert.h>

#include "sbb.hpp"
#include "nlpslv_ipopt.hpp"
#include "polimage.hpp"
#include "cmodel.hpp"
#ifdef MC__USE_CPLEX
  #include "ilcplex/ilocplex.h"
#else
  #include "gurobi_c++.h"
#endif

#undef MC__NLGO_GUROBI_DEBUG
#undef MC__NLGO_GUROBI_TRACE
#undef MC__NLGO_GUROBI_SHOW_BREAKPTS
//#undef MC__USE_CPLEX

extern "C"{
  #include <fenv.h>
  int fedisableexcept( int );
}

namespace mc
{

//! @brief C++ class for global optimization of NLP and MINLP using complete search
////////////////////////////////////////////////////////////////////////
//! mc::NLGO_GUROBI is a C++ class for global optimization of NLP and
//! MINLP using complete search. Relaxations for the nonlinea or
//! nonconvex participating terms are generated using MC++. Further
//! details can be found at: \ref page_NLGO_GUROBI
////////////////////////////////////////////////////////////////////////
template< typename T>
class NLGO_GUROBI:
  public virtual BASE_NLP,
  protected SBB<T>
{
  // Overloading stdout operator
  template <typename U> friend std::ostream& operator<<
    ( std::ostream&os, const NLGO_GUROBI<U>& );

  // Typedef's
#ifdef MC__USE_CPLEX
  typedef std::map< const PolVar<T>*, IloNumVar, lt_PolVar<T> > t_LPVar;
  typedef std::map< const PolCut<T>*, IloRange, lt_PolCut<T> > t_LPCut;
#else
  typedef std::map< const PolVar<T>*, GRBVar, lt_PolVar<T> > t_LPVar;
  typedef std::map< const PolCut<T>*, GRBConstr, lt_PolCut<T> > t_LPCut;
#endif
  typedef std::set< FFVar*, lt_FFVar > t_FFVar;
  using SBB<T>::_exclude_vars;


protected:
  //! @brief number of variables in problem (both independent & dependent)
  unsigned _nvar;
  //! @brief decision variables (both independent & dependent)
  std::vector<FFVar> _var;
  //! @brief variable types in problem (both independent & dependent)
  std::vector<unsigned> _tvar;
  //! @brief number of constraints in problem
  unsigned _nctr;
  //! @brief constraints (types, constraint variables, constraint multipliers), including dependent equations
  std::tuple< std::vector<t_CTR>, std::vector<FFVar>, std::vector<FFVar> > _ctr;

  //! @brief list of operations for objective evaluation
  std::list<const FFOp*> _op_f;
  //! @brief vector of lists of operations for constraint evaluation
  std::vector< std::list<const FFOp*> > _op_g;
  //! @brief list of operations for Lagragian gradient evaluation
  std::list<const FFOp*> _op_dL;

  //! @brief set of linear functions (objective or constraints) in problem
  t_FFVar _fct_lin;

  //! @brief Polyhedral image environment
  PolImg<T> _POLenv;
  //! @brief Storage vector for function evaluation in polyhedral relaxation arithmetic
  std::vector< PolVar<T> > _op_POLfg;

  //! @brief Chebyshev model dimension
  unsigned _CMdim;
  //! @brief Chebyshev excluded variables
  std::set<unsigned> _CMexcl;
  //! @brief Chebyshev model environment
  CModel<T>* _CMenv;
  //! @brief Chebyshev variables
  std::vector< CVar<T> > _CMvar;
  //! @brief Chebyshev objective variable
  CVar<T> _CMobj;
  //! @brief Chebyshev constraint variables
  std::vector< CVar<T> > _CMctr;
  //! @brief Storage vector for function evaluation in Chebyshev arithmetic
  std::vector< CVar<T> > _op_CMfg;

private:
  //! @brief Chebyshev full basis in DAG
  std::vector< FFVar > _basis;
  //! @brief Chebyshev full basis order
  unsigned _basisord;
  //! @brief Chebyshev full basis size
  unsigned _basisdim;
  //! @brief Internal DAG variable for [-1,1] scaled variables
  std::vector< FFVar > _scalvar;
  //! @brief Internal DAG variable for objective function
  FFVar _objvar;
  //! @brief Internal DAG variable for Lagrangian function
  FFVar _lagr;
  //! @brief Lagrangian gradient
  const FFVar* _lagr_grad;

#ifdef MC__USE_CPLEX
  //! @brief CPLEX environment for piecewise relaxation
  IloEnv* _ILOenv;
  //! @brief CPLEX model for piecewise relaxation
  IloModel* _ILOmodel;
  //! @brief CPLEX object for piecewise relaxation
  IloCplex* _ILOcplex;
#else
  //! @brief GUROBI environment for piecewise relaxation
  GRBEnv* _GRBenv;
  //! @brief GUROBI model for piecewise relaxation
  GRBModel* _GRBmodel;
#endif

  //! @brief IPOPT model for local search
  Ipopt::SmartPtr<mc::NLPSLV_IPOPT> _NLPSLV;

  //! @brief Polyhedral image decision variables
  std::vector< PolVar<T> > _POLvar;
  //! @brief Polyhedral image objective auxiliary
  PolVar<T> _POLobjaux;
  //! @brief Polyhedral image scaled decision variables
  std::vector< PolVar<T> > _POLscalvar;
  //! @brief Polyhedral image basis variables
  std::vector< PolVar<T> > _POLbasis;
  //! @brief Polyhedral image objective variable
  PolVar<T> _POLobj;
  //! @brief Polyhedral image variable of objective Chebyshev polynomial
  //PolVar<T> _POLobjcpol;
  //! @brief Polyhedral image constraint variables
  std::vector< PolVar<T> > _POLctr;
  //! @brief Polyhedral image variables of constraint Chebyshev polynomials
  //std::vector< PolVar<T> > _POLctrcpol;

  //! @brief map of LP variables in polyhedral image
  t_LPVar _LPvar;
  //! @brief map of LP cuts in polyhedral image
  t_LPCut _LPcut;
#ifdef MC__USE_CPLEX
  //! @brief LP objective in polyhedral image
  std::pair<bool,IloObjective> _LPobj;
  //! @brief LP constraint for incumbent cut in polyhedral image
  std::pair<bool,IloRange> _LPinc;
#else
  //! @brief GUROBI constraint for incumbent cut in polyhedral image
  std::pair<bool,GRBConstr> _LPinc;
#endif

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

  //! @brief Cleanup gradient/hessian storage
  void _cleanup()
  {
    delete[] _lagr_grad; _lagr_grad = 0;
  }

public:
  /** @defgroup NLGO_GUROBI Nonlinear Global Optimization using GUROBI
   *  @{
   */
  //! @brief Constructor
  NLGO_GUROBI()
    : _nvar(0), _nctr(0), _CMenv(0), _basisord(0), _basisdim(0), _lagr_grad(0)
    {
#ifdef MC__USE_CPLEX
      _ILOenv = new IloEnv;
      _ILOmodel = new IloModel(*_ILOenv);
      _ILOcplex = new IloCplex(*_ILOenv);
      _LPobj.first = false;
#else
      _GRBenv = new GRBEnv();
      _GRBmodel = new GRBModel(*_GRBenv);
#endif
      _LPinc.first = false;
      _NLPSLV = new NLPSLV_IPOPT;
    }

  //! @brief Destructor
  virtual ~NLGO_GUROBI()
    {
      _cleanup();
#ifdef MC__USE_CPLEX
      delete _ILOmodel;
      delete _ILOcplex;
      _ILOenv->end();
      delete _ILOenv;
#else
      delete _GRBmodel;
      delete _GRBenv;
#endif
      delete _CMenv;
      // No need to delete _NLPSLV
    }

  //! @brief Structure holding options for NLGO_GUROBI
  struct Options
  {
    //! @brief Constructor
    Options():
      CSALGO(SBB), RELMETH(DRL), PREPROC(true), DOMREDMAX(10), DOMREDTHRES(0.05),
      DOMREDBKOFF(1e-8), CVATOL(1e-3), CVRTOL(1e-3), LPPRESOLVE(-1),
      LPFEASTOL(1e-9), LPOPTIMTOL(1e-9), MIPRELGAP(1e-7), MIPABSGAP(1e-7),
      PRESOS2BIGM(-1.), CMODSPAR(false), CMODPROP(2), CMODCUTS(0), CMODDMAX(1e20),
      MAXITER(0), MAXCPU(7.2e3), DISPLAY(2), MIPDISPLAY(0), MIPFILE(""),
      POLIMG(), NLPSLV(), CMODEL()
      {
#ifdef MC__USE_CPLEX
        LPALGO = 0;
#else
        LPALGO = -1;
#endif
      }
    //! @brief Assignment operator
    Options& operator= ( Options&options ){
        CSALGO      = options.CSALGO;
        RELMETH     = options.RELMETH;
        PREPROC     = options.PREPROC;
        DOMREDMAX   = options.DOMREDMAX;
        DOMREDTHRES = options.DOMREDTHRES;
        DOMREDBKOFF = options.DOMREDBKOFF;
        CVATOL      = options.CVATOL;
        CVRTOL      = options.CVRTOL;
        LPALGO      = options.LPALGO;
        LPPRESOLVE  = options.LPPRESOLVE;
        LPFEASTOL   = options.LPFEASTOL;
        LPOPTIMTOL  = options.LPOPTIMTOL;
        MIPRELGAP   = options.MIPRELGAP;
        MIPABSGAP   = options.MIPABSGAP;
        PRESOS2BIGM = options.PRESOS2BIGM;
        CMODSPAR    = options.CMODSPAR;
        CMODPROP    = options.CMODPROP;
        CMODCUTS    = options.CMODCUTS;
        CMODDMAX    = options.CMODDMAX;
        MAXITER     = options.MAXITER;
        MAXCPU      = options.MAXCPU;
        DISPLAY     = options.DISPLAY;
        MIPDISPLAY  = options.MIPDISPLAY;
        MIPFILE     = options.MIPFILE;
        POLIMG      = options.POLIMG;
        NLPSLV      = options.NLPSLV;
        CMODEL      = options.CMODEL;
        return *this;
      }
    //! @brief Complete search methods
    enum CS{
      SBB=0,	//!< Spatial branch-and-bound
      PWL	//!< Hierarchy of piecewise-linear relaxations
    };
    //! @brief Relaxation methods
    enum RELAX{
      DRL=0,	//!< Decomposition-relaxation-linearization (Tawarmalani & Sahinidis)
      CHEB,	//!< Chebyshev-derived relaxations, controlled by parameters CMODPROP and CMODCUT
      HYBRID    //!< Combination of DRL and CHEB
    };
    //! @brief Complete search algorithm
    CS CSALGO;
    //! @brief Relaxation method
    RELAX RELMETH;
    //! @brief Display
    void display
      ( std::ostream&out ) const;
    //! @brief Preprocessing (local optimization, domain reduction)
    bool PREPROC;
    //! @brief Maximum number of domain reduction rounds
    unsigned DOMREDMAX;
    //! @brief Threshold for repeating domain reduction (minimum domain reduction ratio)
    double DOMREDTHRES;
    //! @brief Backoff of reduced variable bounds to compensate for numerical errors
    double DOMREDBKOFF;
    //! @brief Convergence absolute tolerance
    double CVATOL;
    //! @brief Convergence relative tolerance
    double CVRTOL;
    //! @brief LP algorithm
    int LPALGO;
    //! @brief LP presolve\f
    int LPPRESOLVE;
    //! @brief Tolerance on LP feasibility
    double LPFEASTOL;
     //! @brief Tolerance on LP optimality
    double LPOPTIMTOL;
    //! @brief Tolerance on relative MIP gap
    double MIPRELGAP;
    //! @brief Tolerance on absolute MIP gap
    double MIPABSGAP;
    //! @brief Parameter controlling SOS2 reformulations in Gurobi
    double PRESOS2BIGM;
     //! @brief Sparse Chebyshev model propagation
    bool CMODSPAR;
    //! @brief Chebyhev model propagation order (0: no propag.)
    unsigned CMODPROP;
    //! @brief Chebyhev model cut order (0: same as propag.)
    unsigned CMODCUTS;
    //! @brief Chebyhev model maximum diameter for cut generation
    double CMODDMAX;
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
    //! @brief PolImg options
    typename PolImg<T>::Options POLIMG;
    //! @brief NLPSLV_IPOPT options
    typename NLPSLV_IPOPT::Options NLPSLV;
    //! @brief CModel options
    typename CModel<T>::Options CMODEL;
  };
  //! @brief NLGO_GUROBI options
  Options options;

  //! @brief Class managing exceptions for NLGO_GUROBI
  class Exceptions
  {
  public:
    //! @brief Enumeration type for NLGO_GUROBI exception handling
    enum TYPE{
      NLPOBJ=1,		//!< Undefined objective function in optimization problem
      INTERN=-33	//!< Internal error
    };
    //! @brief Constructor for error <a>ierr</a>
    Exceptions( TYPE ierr ) : _ierr( ierr ){}
    //! @brief Inline function returning the error flag
    int ierr(){ return _ierr; }
    //! @brief Inline function returning the error description
    std::string what(){
      switch( _ierr ){
      case NLPOBJ:
        return "NLGO_GUROBI::Exceptions  Undefined objective function in optimization problem";
      case INTERN: default:
        return "NLGO_GUROBI::Exceptions  Internal error";
      }
    }
  private:
    TYPE _ierr;
  };

  //! @brief Structure holding current statistics
  struct Stats{
    void reset(){ tPOLIMG = tLPSOL = tLPSET = tNLPSOL = 0.; nLPSOL = 0; tstart = time(); }
    double tstart;
    double tPOLIMG;
    double tLPSOL;
    unsigned nLPSOL;
    double tLPSET;
    double tNLPSOL;
  } stats;

  //! @brief Setup DAG for cost and constraint evaluation
  void setup
    //();
    ( std::set<unsigned> CMexcl=std::set<unsigned>() );

  //! @brief Solve (continuous) optimization model to local optimality using IPOPT in variable range <a>P</a> and from initial point <a>p0</a> -- return value is IPOPT status
  Ipopt::ApplicationReturnStatus local
    ( const T*P, const double*p0=0, const bool reset=true );

  //! @brief Solve relaxed optimization model (or relaxed feasibility model if <a>feastest</a> is true) using GUROBI within variable range <a>P</a> and for variable types <a>tvar</a> -- return value is <A href="http://www.gurobi.com/documentation/6.0/refman/optimization_status_codes.html">GUROBI status</A>
  int relax
    ( const T*P, const unsigned*tvar=0, const unsigned refine=0,
      const bool reset=true, const bool feastest=false );

  //! @brief Solve bound contraction problems from relaxed model using GUROBI, starting with variable range <a>P</a>, for variable types <a>tvar</a> and incumbent value <a>inc</a> (or constraint backoff if <a>feastest</a> is true), and using the options specified in <a>NLGO_GUROBI::Options::DOMREDMAX</a> and <a>NLGO_GUROBI::Options::DOMREDTHRES</a> -- return is <A href="http://www.gurobi.com/documentation/6.0/refman/optimization_status_codes.html">GUROBI status</A>, updated variable bounds <a>P</a>, and number of iterative refinements <a>nred</a> 
  int contract
    ( T*P, unsigned&nred, const unsigned*tvar=0, const double*inc=0,
      const bool reset=true, const bool feastest=false );

  //! @brief Solve bound contraction problem for lower/upper component <a>ip</a> from relaxed model using GUROBI, starting with variable range <a>P</a> and incumbent value <a>inc</a> (or constraint backoff if <a>feastest</a> is true) -- return is <A href="http://www.gurobi.com/documentation/6.0/refman/optimization_status_codes.html">GUROBI status</A>
  int contract
    ( const T*P, const unsigned ip, const bool uplo, const unsigned*tvar=0,
      const double*inc=0, const bool reset=true, const bool feastest=false );

  //! @brief Solve optimization model to global optimality within variable range <a>P</a> and for variable types <a>tvar</a> -- return value is <A href="http://www.gurobi.com/documentation/6.0/refman/optimization_status_codes.html">GUROBI status</A>
  int solve
    ( const T*P, const unsigned*tvar=0, const double*p0=0,
      std::ostream&os=std::cout );

  //! @brief Value of DAG variable <a>X</a> after last LP optimization
  double get_variable
    ( const FFVar&X ) const
    {
      auto itp = _POLenv.Vars().find( const_cast<FFVar*>(&X) );
      auto itv = _LPvar.find( itp->second );
#ifdef MC__USE_CPLEX
      return _ILOcplex->getValue( itv->second );
#else
      //GRBVar XGRB = itv->second;
      //return XGRB.get(GRB_DoubleAttr_X);
      return itv->second.get(GRB_DoubleAttr_X);
#endif
    }

  //! @brief Optimal cost value after last LP optimization
  double get_objective() const
    {
#ifdef MC__USE_CPLEX
      return _ILOcplex->getObjValue();
#else
      return _GRBmodel->get( GRB_DoubleAttr_ObjVal );
#endif
    }

  //! @brief Status after last LP optimization
  int get_status() const
    {
#ifdef MC__USE_CPLEX
      return _ILOcplex->getStatus();
#else
      return _GRBmodel->get( GRB_IntAttr_Status );
#endif
    }

  //! @brief Check optimal status after last LP optimization
  bool optimal_status() const
    {
#ifdef MC__USE_CPLEX
      return( _ILOcplex->getStatus() == IloAlgorithm::Optimal );
#else
      return( _GRBmodel->get( GRB_IntAttr_Status ) == GRB_OPTIMAL );
#endif
    }

  //! @brief Pointer to last LP optimization model and solution
#ifdef MC__USE_CPLEX
  const IloModel* get_relaxed_model() const
    { return _ILOmodel; }
#else
  const GRBModel* get_relaxed_model() const
    { return _GRBmodel; }
#endif

  //! @brief Local solution information after last NLP optimization
  const NLPSLV_IPOPT::SOLUTION& get_local_solution() const
    {
      return _NLPSLV->solution();
    }
  /** @} */
protected:
  //! @brief Set model polyhedral relaxation
  void _set_polrelax
    ( const T*P, const unsigned*tvar, const bool feastest=false );
  //! @brief Refine model polyhedral relaxation by adding breakpoints
  void _refine_polrelax
    ( const T*P, const unsigned*tvar, const bool feastest=false );
  //! @brief Update model polyhedral relaxation (bounds and types)
  void _update_polrelax
    ( const T*P, const unsigned*tvar, const bool feastest=false );

  //! @brief Set relaxed objective in LP model
  void _set_LPrelax
    ( const bool feastest=false );
  //! @brief Set parameter bound contracting objective in LP model
  void _set_LPcontract
    ( const unsigned ip, const bool uplo, const double*inc,
      const bool feastest=false );
  //! @brief Solve LP model
  void _solve_LPmodel
    ();
 
 //! @brief Tighten parameter bounds <a>P</a> for NLP model relaxation and incumbent value <a>inc</a> (or constraint backoff is <a>feastest</a> is true) -- return value is <A href="http://www.gurobi.com/documentation/6.0/refman/optimization_status_codes.html">GUROBI status</A>
  int _contract
    ( T*P, const unsigned*tvar=0, const double*inc=0, const bool reset=true,
      const bool feastest=false );
  //! @brief Function computing Hausdorff distance between intervals
  template <typename U> static double _dH
    ( const U&X, const U&Y );
  //! @brief Function computing Hausdorff distance between interval vectors
  template <typename U> static double _reducrel
    ( const unsigned n, const U*Xred, const U*X );

  //! @brief Value of PolImg variable <a>X</a> after last LP optimization
  double _get_variable
    ( const PolVar<T>&X ) const
    {
      auto itv = _LPvar.find( const_cast<PolVar<T>*>(&X) );
#ifdef MC__USE_CPLEX
      return _ILOcplex->getValue( itv->second );
#else
      //GRBVar XGRB = itv->second;
      //return XGRB.get(GRB_DoubleAttr_X);
      return itv->second.get(GRB_DoubleAttr_X);
#endif
    }

private:
  //! @brief Set model linear cuts
  void _set_lincuts
    ( const bool feastest );
  //! @brief Set model polyhedral outer-approximation cuts
  void _set_poacuts
    ( const bool feastest );
  //! @brief Set model polyhedral Chebyshev-derived cuts
  void _set_chebcuts
    ( const bool feastest );

  //! @brief Set options of LP model
  void _set_LPoptions
    ();
  //! @brief Reset LP model
  void _reset_LPmodel
    ();
  //! @brief Set variables and cuts in LP model
  void _set_LPcuts
    ();
  //! @brief Append polyhedral image variable/auxiliary into LP model
  std::pair<typename t_LPVar::iterator,bool> _set_LPvar
    ( const PolVar<T>*pVar );
  //! @brief Append polyhedral image cut into LP model
  void _set_LPcut
    ( const PolCut<T>*pCut );

  //! @brief Set options of SBB solver
  void _set_SBBoptions
    ();
  //! @brief User-function to subproblems in SBB
  typename SBB<T>::STATUS subproblems
    ( const typename SBB<T>::TASK task, const unsigned int np, T*P,
      double*p, double&f, const double INC );
  //! @brief Subproblem for local optimization
  typename SBB<T>::STATUS _subproblem_local
    ( const T*P, double*p, double&f );
  //! @brief Subproblem for feasibility test
  typename SBB<T>::STATUS _subproblem_feasibility
    ( const double*p, double&f );
  //! @brief Subproblem for relaxation
  typename SBB<T>::STATUS _subproblem_relaxed
    ( T*P, double*p, double&f, const double INC );

  //! @brief Solve optimization model to global optimality within variable range <a>P</a> and for variable types <a>tvar</a> -- return value is <A href="http://www.gurobi.com/documentation/6.0/refman/optimization_status_codes.html">GUROBI status</A>
  int _solve_pwl
    ( const T*P, const unsigned*tvar=0, const double*p0=0,
      std::ostream&os=std::cout );
  //! @brief Solve optimization model to global optimality within variable range <a>P</a> and for variable types <a>tvar</a> -- return value is mc::SBB::STATUS
  int _solve_sbb
    ( const T*P, const unsigned*tvar=0, const double*p0=0,
      std::ostream&os=std::cout );

  //! @brief Initialize display
  void _display_init
    ();
  //! @brief Final display
  void _display_final
    ( const unsigned iter );
  //! @brief Add double to display
  void _display_add
    ( const double dval );
  //! @brief Add unsigned int to display
  void _display_add
    ( const int ival );
  //! @brief Add string to display
  void _display_add
    ( const std::string &sval );
  //! @brief Display current buffer stream and reset it
  void _display_flush
    ( std::ostream&os );

  //! @brief Private methods to block default compiler methods
  NLGO_GUROBI
    (const NLGO_GUROBI&);
  NLGO_GUROBI& operator=
    (const NLGO_GUROBI&);
};

template <typename T> inline void
NLGO_GUROBI<T>::setup
//()
( std::set<unsigned> CMexcl )
{
  // full set of decision variables (independent & dependent)
  _var = BASE_NLP::_var;
  _var.insert( _var.end(), _dep.begin(), _dep.end() );
  _nvar = _var.size();

  // full set of constraints (main & equation system)
  _ctr = BASE_NLP::_ctr;
  for( auto its=_sys.begin(); its!=_sys.end(); ++its ){
    std::get<0>(_ctr).push_back( EQ );
    std::get<1>(_ctr).push_back( (*its) );
    std::get<2>(_ctr).push_back( FFVar( _dag ) );
  }
  _nctr = std::get<0>(_ctr).size();

  // setup [-1,1] scaled variables for Chebyshev model environment
  // do not downsize to avoid adding more variables into DAG
  _CMexcl = CMexcl;
  //_CMexcl.clear();
  _CMdim = 0;
  for( unsigned i=0; i<_nvar; i++ )
    if( _CMexcl.find(i) == _CMexcl.end() ) ++_CMdim;
  const unsigned nscalvar = _scalvar.size();
  if( nscalvar < _CMdim ){
    _scalvar.resize( _CMdim );
    for( unsigned i=nscalvar; i<_CMdim; i++ )
      _scalvar[i].set( _dag );
  }

  // setup for objective and constraints evaluation
  _objvar.set( _dag );
  if( std::get<0>(_obj).size() > 1 ) throw Exceptions( Exceptions::NLPOBJ );
  if( std::get<0>(_obj).size() ) _op_f = _dag->subgraph( 1, std::get<1>(_obj).data() );
  _op_g.resize( _nctr );
  unsigned nfgop = _op_f.size();
  for( unsigned i=0; i<_nctr; i++ ){
    _op_g[i] = _dag->subgraph( 1, std::get<1>(_ctr).data()+i );
    if( nfgop <  _op_g[i].size() ) nfgop = _op_g[i].size();
  }
  _op_CMfg.resize(nfgop);
  _op_POLfg.resize(nfgop);
  //_op_g = _dag->subgraph( std::get<1>(_ctr).size(), std::get<1>(_ctr).data() );

  // identify linear objective/constraint functions in NLP
  _fct_lin.clear();
  if( std::get<0>(_obj).size() ){
    FFDep fdep = std::get<1>(_obj)[0].dep();
    auto it = fdep.dep().cbegin();
    bool islin = true;
    for( ; islin && it != fdep.dep().cend(); ++it )
      if( !it->second ) islin = false;
    if( islin ) _fct_lin.insert( &std::get<1>(_obj)[0] );
  }
  for( unsigned j=0; j<_nctr; j++ ){
    FFDep gdep = std::get<1>(_ctr)[j].dep();
    auto it = gdep.dep().cbegin();
    bool islin = true;
    for( ; islin && it != gdep.dep().cend(); ++it )
      if( !it->second ) islin = false;
    if( islin ) _fct_lin.insert( &std::get<1>(_ctr)[j] );
  }
  std::cout << "LINEAR OBJECTIVE/CONSTRAINT FUNCTIONS:    " << _fct_lin.size() << std::endl
            << "NONLINEAR OBJECTIVE/CONSTRAINT FUNCTIONS: " << _nctr-_fct_lin.size()+1 << std::endl;

  // identify linear variables in NLP and exclude from branching
  FFDep fgdep = std::get<0>(_obj).size()? std::get<1>(_obj)[0].dep(): 0.;
  for( unsigned j=0; j<_nctr; j++ )
    fgdep += std::get<1>(_ctr)[j].dep();
#ifdef MC__NLGO_GUROBI_DEBUG
  std::cout << "NLP <- " << fgdep << std::endl;
  //int dum; std::cin >> dum;
#endif
  _exclude_vars.clear();
  unsigned nvar_lin = 0;
  for( unsigned i=0; i<_nvar; i++ ){
    auto it = fgdep.dep().find(i);
    if( it == fgdep.dep().end() || it->second ){
      _exclude_vars.insert(i);
      nvar_lin++;
    }
  }
  std::cout << "LINEARLY PARTICIPATING VARIABLES:         " << nvar_lin << std::endl
            << "NONLINEARLY PARTICIPATING VARIABLES:      " << _nvar-nvar_lin << std::endl
            << std::endl;

/*
  // setup Lagrangian gradient evaluation
  _cleanup();
  FFVar lagr = std::get<1>(_obj)[0] * std::get<2>(_obj)[0];
  for( unsigned ic=0; ic<_nctr; ic++ )
    lagr += std::get<1>(_ctr)[ic] * std::get<2>(_ctr)[ic];
  //_lagr_grad = _dag->FAD( 1, &lagr, _nvar, _var.data() );
  _lagr_grad = _dag->BAD( 1, &lagr, _nvar, _var.data() );
  _op_dL = _dag->subgraph( _nvar, _lagr_grad );
*/
}

template <typename T> inline Ipopt::ApplicationReturnStatus
NLGO_GUROBI<T>::local
( const T*P, const double*p0, const bool reset )
{
  stats.tNLPSOL -= time();

  // Set-up local NLP solver
  if( reset ){
    _NLPSLV->options = options.NLPSLV;
    _NLPSLV->set( *this );
    _NLPSLV->setup();
  }

  // Solve local NLP model
  std::vector<double> p(_nvar);
  for( unsigned ip=0; ip<_nvar; ip++ )
    p[ip] = ( p0? p0[ip]: Op<T>::mid(P[ip]) );
  Ipopt::ApplicationReturnStatus status = _NLPSLV->solve( P, p.data() );

  stats.tNLPSOL += time();
  return status;
}

template <typename T> inline int
NLGO_GUROBI<T>::relax
( const T*P, const unsigned*tvar, const unsigned nref, const bool reset,
  const bool feastest )
{
  // Reset polyhedral image, LP variables and cuts on request
  if( reset ) _set_polrelax( P, tvar, feastest );
#ifdef MC__NLGO_GUROBI_DEBUG
  std::cout << _POLenv;
#endif

  for( unsigned iref=0; ; iref++ ){
    // Set-up relaxed objective, options, and solve GUROBI model
    _set_LPrelax( feastest );
    _solve_LPmodel();

    // Break if relaxation unsuccessful or refinement iteration exceeded
    if( !optimal_status() || iref >= nref ) break;

    // Refine relaxation via additional breakpoints
    _refine_polrelax( P, tvar );
  }

  return get_status();
}

template <typename T> inline int
NLGO_GUROBI<T>::contract
( const T*P, const unsigned ip, const bool uplo, const unsigned*tvar,
  const double*inc, const bool reset, const bool feastest )
{
  // Reset polyhedral image, LP variables and cuts on request
  if( reset ) _set_polrelax( P, tvar, feastest );
#ifdef MC__NLGO_GUROBI_DEBUG
  std::cout << "\nTIGHTENING OF VARIABLE " << ip << (uplo?"U":"L") << ":\n";
  std::cout << _POLenv;
#endif

  // Set-up relaxed objective, options, and solve GUROBI model
  _set_LPcontract( ip, uplo, inc, feastest );
  _solve_LPmodel();

  return get_status();
}

template <typename T> inline int
NLGO_GUROBI<T>::_contract
( T*P, const unsigned*tvar, const double*inc, const bool reset,
  const bool feastest )
{
  std::multimap<double,std::pair<unsigned,bool>> vardomred, vardomredupd;
  std::pair<unsigned,bool> varini;
  for( unsigned ip=0; ip<_nvar; ip++ ){
    varini.first = ip;
    varini.second = false; // lower bound
    vardomred.insert( std::pair<double,std::pair<unsigned,bool>>(1.,varini) );
    varini.second = true; // upper bound
    vardomred.insert( std::pair<double,std::pair<unsigned,bool>>(1.,varini) );
  }
  if( reset ) _set_polrelax( P, tvar, feastest );

  // solve reduction subproblems from closest to farthest from bounds
  unsigned nred=0;
  for( ; !vardomred.empty(); nred++ ){
    // upper/lower range reduction for current subproblem
    auto itv = vardomred.begin();
    const unsigned ip = (*itv).second.first;
    const bool uplo = (*itv).second.second;
    double pL = Op<T>::l( P[ip] ), pU = Op<T>::u( P[ip] );
    contract( P, ip, uplo, tvar, inc, false, feastest );
    if( !optimal_status() ) return get_status();
    switch( uplo ){
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
#ifdef MC__NLGO_GUROBI_DEBUG
    std::cout << "  UPDATED RANGE OF VARIABLE #" << ip << ": " << P[ip] << std::endl;
#endif
    // update map of candidate reduction subproblems
    vardomredupd.clear();
    for( ++itv; itv!=vardomred.end(); ++itv ){
      const unsigned ip = (*itv).second.first;
      const bool uplo = (*itv).second.second;
      double dist = uplo? std::fabs( get_variable( _var[ip] ) - Op<T>::u( P[ip] ) ):
                          std::fabs( get_variable( _var[ip] ) - Op<T>::l( P[ip] ) );
      if( dist <= options.DOMREDTHRES ) continue;
      if( dist > (*itv).first ) dist = (*itv).first;
      vardomredupd.insert( std::pair<double,std::pair<unsigned,bool>>(dist,std::make_pair(ip,uplo)) );
    }
    vardomred.swap( vardomredupd );
  }
#ifdef MC__NLGO_GUROBI_DEBUG
    std::cout << "SOLVED " << nred << " RANGE REDUCTION LPs OUT OF " << 2*_nvar << std::endl;
#endif
  return get_status();
}

template <typename T> inline int
NLGO_GUROBI<T>::contract
( T*P, unsigned&nred, const unsigned*tvar, const double*inc,
  const bool reset, const bool feastest )
{
  // Main loop for relaxation and domain reduction
  std::vector<T> P0;
  for( nred=0; nred < options.DOMREDMAX; nred++ ){
    P0.assign( P, P+_nvar );
    _contract( P, tvar, inc, nred? true: reset, feastest );
    if( !optimal_status() ) return get_status();
    double vred = _reducrel( _nvar, P, P0.data() );
    if( vred < options.DOMREDTHRES ) break;
#ifdef MC__NLGO_GUROBI_DEBUG
    std::cout << "Reduction: " << vred*1e2 << "%\n";
#endif
  }
  return get_status();
}

template <typename T> template <typename U> inline double
NLGO_GUROBI<T>::_dH
( const U&X, const U&Y )
{
  return std::max( std::fabs(Op<U>::l(X)-Op<U>::l(Y)),
                   std::fabs(Op<U>::u(X)-Op<U>::u(Y)) );
}

template <typename T> template <typename U> inline double
NLGO_GUROBI<T>::_reducrel
( const unsigned n, const U*Xred, const U*X )
{
  double drel = 0.;
  for( unsigned ip=0; ip<n; ip++ )
    drel = std::max( drel, _dH( Xred[ip], X[ip] ) / Op<T>::diam( X[ip] ) );
  return drel;
}

template <typename T>
inline typename SBB<T>::STATUS
NLGO_GUROBI<T>::_subproblem_local
( const T*P, double*p, double&f )
{
  // Solve local NLP model (integer variable bounds fixed to relaxed solution if MINLP)
  std::vector<T> P_loc;
  unsigned *tvar = !_tvar.empty()? _tvar.data(): 0;
  for( unsigned ivar=0; ivar<_nvar; ivar++ )
    P_loc.push_back( (tvar && tvar[ivar])? p[ivar]: P[ivar] );
  _NLPSLV->solve( P_loc.data(), p );

  // Update incumbent
  if( _NLPSLV->solution().status == Ipopt::Solve_Succeeded
   || _NLPSLV->solution().status == Ipopt::Solved_To_Acceptable_Level
   || _NLPSLV->solution().status == Ipopt::Feasible_Point_Found ){
    // Copy solution value and solution point
    f = _NLPSLV->solution().f;
#ifdef MC__NLGO_GUROBI_DEBUG
    std::cout << std::scientific << std::setprecision(4)
              << "  f_loc = " << f << std::endl;
#endif
    for( unsigned ivar=0; ivar<_nvar; ivar++ ){
      p[ivar] = _NLPSLV->solution().p[ivar];
#ifdef MC__NLGO_GUROBI_DEBUG
      std::cout << std::scientific << std::setprecision(4)
                << "  p_loc(" << ivar << ") = " << p[ivar] << std::endl;
#endif
    }
    return SBB<T>::NORMAL;
  }
  else 
    return SBB<T>::FAILURE;
}

template <typename T>
inline typename SBB<T>::STATUS
NLGO_GUROBI<T>::_subproblem_relaxed
( T*P, double*p, double&f, const double INC )
{
  // Apply domain contraction
  // Do NOT test for infeasibility here, because contraction problem may
  // become infeasible due to round-off in LP solver
  unsigned *tvar = !_tvar.empty()? _tvar.data(): 0;
  _update_polrelax( P, tvar );
  unsigned nred = 0;
  contract( P, nred, tvar, &INC, false );

  // Solve relaxed problem
  _update_polrelax( P, tvar );
  relax( P, tvar, 0, false );
  switch( get_status() ){
#ifdef MC__USE_CPLEX
   case IloAlgorithm::Optimal:{
#else
   case GRB_OPTIMAL:{
#endif
    f = get_objective();
#ifdef MC__NLGO_GUROBI_DEBUG
    std::cout << std::scientific << std::setprecision(4)
              << "  f_rel = " << f << std::endl;
#endif
    unsigned ivar = 0;
    for( auto itv=_var.begin(); itv!=_var.end(); ++itv, ivar++ ){
      p[ivar] = get_variable( *itv );
#ifdef MC__NLGO_GUROBI_DEBUG
      std::cout << std::scientific << std::setprecision(4)
                << "  p_rel(" << ivar << ") = " << p[ivar] << std::endl;
#endif
    }
    return SBB<T>::NORMAL;
   }

#ifdef MC__USE_CPLEX
   case IloAlgorithm::Infeasible:
   //case IloAlgorithm::InfeasibleOrUnbounded:
#else
   case GRB_INFEASIBLE:
   //case GRB_INF_OR_UNBD:
#endif
    return SBB<T>::INFEASIBLE;

   default:
    return SBB<T>::FAILURE;
  }
}

template <typename T>
inline typename SBB<T>::STATUS
NLGO_GUROBI<T>::_subproblem_feasibility
( const double*p, double&f )
{
  // Compute objective function value
  _dag->eval( _op_f, 1, std::get<1>(_obj).data(), &f, _nvar, _var.data(), p );

  // Check current point feasibility
  auto itc=std::get<0>(_ctr).begin();
  double gj, maxinfeas=0.;
  for( unsigned j=0; itc!=std::get<0>(_ctr).end(); ++itc, j++ ){
    _dag->eval( _op_g[j], 1, std::get<1>(_ctr).data()+j, &gj, _nvar, _var.data(), p );
    switch( (*itc) ){
      case EQ: maxinfeas = std::fabs(gj); break;
      case LE: maxinfeas = gj;            break;
      case GE: maxinfeas = -gj;           break;
    }
    if( maxinfeas > _NLPSLV->options.PRIMALTOL ) return SBB<T>::INFEASIBLE;
  }
  return SBB<T>::NORMAL;  
}

template <typename T>
inline typename SBB<T>::STATUS
NLGO_GUROBI<T>::subproblems
( const typename SBB<T>::TASK task, const unsigned int np, T*P,
  double*p, double&f, const double INC )
{
  switch( task ){

  // UPPER BOUND
  case SBB<T>::UPPERBD:
   try{ 
    switch( std::get<0>(_obj)[0] ){
      case MIN: return _subproblem_local( P, p, f );
      case MAX: return _subproblem_relaxed( P, p, f, INC );
    }
   }
   catch(...){
     return SBB<T>::FAILURE;
   }

  // LOWER BOUND
  case SBB<T>::LOWERBD:
   try{ 
    switch( std::get<0>(_obj)[0] ){
      case MAX: return _subproblem_local( P, p, f );
      case MIN: return _subproblem_relaxed( P, p, f, INC );
    }
   }
   catch(...){
     return SBB<T>::FAILURE;
   }

  // FEASIBILITY TEST
  case SBB<T>::FEASTEST:
    return _subproblem_feasibility( p, f );

  // PRE/POST-PROCESSING
  case SBB<T>::PREPROC:
  case SBB<T>::POSTPROC:
    return SBB<T>::NORMAL;

  default:
    return SBB<T>::FAILURE;
  }
}

template <typename T> inline int
NLGO_GUROBI<T>::solve
( const T*P0, const unsigned*tvar, const double*p0, std::ostream&os )
{
  switch( options.CSALGO ){
   case Options::SBB: default:
     return _solve_sbb( P0, tvar, p0, os );
   case Options::PWL:
     return _solve_pwl( P0, tvar, p0, os );
  }
}

template <typename T> inline int
NLGO_GUROBI<T>::_solve_sbb
( const T*P0, const unsigned*tvar, const double*p0, std::ostream&os )
{
  stats.reset();
  bool isMINLP = false;
  for( unsigned ivar=0; tvar && !isMINLP && ivar<_nvar; ivar++ )
    if( tvar[ivar] ) isMINLP = true;

  // Set-up local NLP solver
  _NLPSLV->options = options.NLPSLV;
  _NLPSLV->set( *this );
  _NLPSLV->setup();
  std::vector<T> P_loc(_nvar);
  std::vector<double> p_loc(_nvar);
  _p_inc.clear();

  // Set-up NLP relaxed solver
  std::vector<T> P( P0, P0+_nvar );
  _set_polrelax( P.data(), tvar );
  //std::cout << *_dag;
  //std::cout << _POLenv;
  //return 0;

  // Preprocessing
  if( options.PREPROC ){
    local( P.data(), p0 );

    // Rounding solution point to nearest integer in case of MINLP
    if( _NLPSLV->solution().status == Ipopt::Solve_Succeeded && isMINLP ){
      for( unsigned ivar=0; ivar<_nvar; ivar++ )
        P_loc[ivar] = tvar[ivar]? std::round( _NLPSLV->solution().p[ivar] ): P[ivar];
      local( P_loc.data(), _NLPSLV->solution().p );
    }

    // Update incumbent
    if( _NLPSLV->solution().status == Ipopt::Solve_Succeeded
     || _NLPSLV->solution().status == Ipopt::Solved_To_Acceptable_Level
     || _NLPSLV->solution().status == Ipopt::Feasible_Point_Found ){
      if( _p_inc.empty()
       || ( std::get<0>(_obj)[0]==MIN && _NLPSLV->solution().f<_f_inc )
       || ( std::get<0>(_obj)[0]==MAX && _NLPSLV->solution().f>_f_inc ) ){
        _p_inc.resize(_nvar);
        for( unsigned ivar=0; ivar<_nvar; ivar++ )
          _p_inc[ivar] = _NLPSLV->solution().p[ivar];
        _f_inc = _NLPSLV->solution().f;
      }
    }
  }

  // Set variable types and exclusions from B&B
  tvar? _tvar.assign( tvar, tvar+_nvar ): _tvar.clear();
  std::set<unsigned int>& exclude = _exclude_vars;
  _set_SBBoptions();
  SBB<T>::variables( _nvar, P0, exclude );

  // Run SBB solver
  const double* fINC = (!_p_inc.empty()? &_f_inc: 0);
  const double* pINC = (!_p_inc.empty()? _p_inc.data(): p0);
  switch( std::get<0>(_obj)[0] ){
    case MIN: SBB<T>::solve( SBB<T>::MINIMIZE, pINC, fINC, os ); break;
    case MAX: SBB<T>::solve( SBB<T>::MAXIMIZE, pINC, fINC, os ); break;
  }
  if( SBB<T>::_p_inc ) _p_inc.assign( SBB<T>::_p_inc, SBB<T>::_p_inc+_nvar );
  _f_inc = SBB<T>::_f_inc;
  return SBB<T>::_status;
}

template <typename T> inline int
NLGO_GUROBI<T>::_solve_pwl
( const T*P0, const unsigned*tvar, const double*p0, std::ostream&os )
{
  stats.reset();
  _display_init();

  unsigned nred = 0;
  bool isMINLP = false;
  for( unsigned ivar=0; tvar && !isMINLP && ivar<_nvar; ivar++ )
    if( tvar[ivar] ) isMINLP = true;

  // Set-up local NLP solver
  _NLPSLV->options = options.NLPSLV;
  _NLPSLV->set( *this );
  _NLPSLV->setup();
  std::vector<T> P_loc(_nvar);
  std::vector<double> p_loc(_nvar);
  _p_inc.clear();

  // Set-up NLP relaxed solver
  std::vector<T> P( P0, P0+_nvar );
  _set_polrelax( P.data(), tvar );
  //std::cout << *_dag;
  //std::cout << _POLenv;
  //return 0;

  // Preprocessing
  if( options.PREPROC ){
    local( P.data(), p0 );

    // Rounding solution point to nearest integer in case of MINLP
    if( _NLPSLV->solution().status == Ipopt::Solve_Succeeded && isMINLP ){
      for( unsigned ivar=0; ivar<_nvar; ivar++ )
        P_loc[ivar] = tvar[ivar]? std::round( _NLPSLV->solution().p[ivar] ): P[ivar];
      local( P_loc.data(), _NLPSLV->solution().p );
    }

    // Update incumbent
    if( _NLPSLV->solution().status == Ipopt::Solve_Succeeded
     || _NLPSLV->solution().status == Ipopt::Solved_To_Acceptable_Level
     || _NLPSLV->solution().status == Ipopt::Feasible_Point_Found ){
      if( _p_inc.empty()
       || ( std::get<0>(_obj)[0]==MIN && _NLPSLV->solution().f<_f_inc )
       || ( std::get<0>(_obj)[0]==MAX && _NLPSLV->solution().f>_f_inc ) ){
        _p_inc.resize(_nvar);
        for( unsigned ivar=0; ivar<_nvar; ivar++ )
          _p_inc[ivar] = _NLPSLV->solution().p[ivar];
        _f_inc = _NLPSLV->solution().f;
      }
    }

    // Apply domain contraction
    // Do NOT test for infeasibility here, because contraction problem may
    // become infeasible due to round-off in LP solver
    contract( P.data(), nred, tvar, !_p_inc.empty()? &_f_inc: 0, false );
  }

  // Iterative relaxation solution and refinement
  unsigned iter = 1;
  for( ; ; ++iter ){
    // iteration and preprocessing display
    std::ostringstream odisp, onred;
    odisp << std::right << std::setw(_IPREC) << iter;
    if( nred ) onred << "  R" << nred;
    else       onred << "-";
    odisp << std::right << std::setw(_IPREC) << onred.str();
    _display_add( odisp.str() );

    // Set-up relaxation and solve with GUROBI
    //_update_polrelax( P.data(), tvar );
    relax( P.data(), tvar, 0, false );
    switch( get_status() ){
#ifdef MC__USE_CPLEX
     case IloAlgorithm::Optimal:
#else
     case GRB_OPTIMAL:
#endif
      _display_add( get_objective() );
      break;
#ifdef MC__USE_CPLEX
     case IloAlgorithm::Infeasible:
     //case IloAlgorithm::InfeasibleOrUnbounded:
#else
     case GRB_INFEASIBLE:
     //case GRB_INF_OR_UNBD:
#endif
      _display_add( "-" );
      _display_add( "-" );
      _display_add( time()-stats.tstart );
      _display_add( "INFEASIBLE" );
      _display_flush( os );
      _display_final( iter );
      _display_flush( os );
      return get_status();
     default:
      _display_add( "-" );
      _display_add( "-" );
      _display_add( time()-stats.tstart );
      _display_add( "FAILURE" );
      _display_flush( os );
      _display_final( iter );
      _display_flush( os );
      return get_status();
    }

    // Solve local NLP model (integer variable bounds fixed to relaxed solution if MINLP)
    unsigned ivar = 0;
    for( auto itv=_var.begin(); itv!=_var.end(); ++itv, ivar++ ){
      p_loc[ivar] = get_variable( *itv );
      P_loc[ivar] = (tvar && tvar[ivar])? p_loc[ivar]: P[ivar];
    }
    _NLPSLV->solve( P_loc.data(), p_loc.data() );

    // Update incumbent
    if( _NLPSLV->solution().status == Ipopt::Solve_Succeeded
     || _NLPSLV->solution().status == Ipopt::Solved_To_Acceptable_Level
     || _NLPSLV->solution().status == Ipopt::Feasible_Point_Found ){
      if( _p_inc.empty()
       || ( std::get<0>(_obj)[0]==MIN && _NLPSLV->solution().f<_f_inc )
       || ( std::get<0>(_obj)[0]==MAX && _NLPSLV->solution().f>_f_inc ) ){
        _p_inc.resize(_nvar);
        for( unsigned ivar=0; ivar<_nvar; ivar++ )
          _p_inc[ivar] = _NLPSLV->solution().p[ivar];
        _f_inc = _NLPSLV->solution().f;
      }
    }

    // Test termination criteria
    if( !_p_inc.empty() ){
      _display_add( _f_inc );
      _display_add( time()-stats.tstart );
      if( std::fabs( _f_inc - get_objective() ) <= options.CVATOL 
       || std::fabs( _f_inc - get_objective() ) <= 0.5 * options.CVRTOL
          * std::fabs( _NLPSLV->solution().f + get_objective() ) ){
        _display_add( "OPTIMAL" );
        _display_flush( os );
        break;
      }
    }
    else{
      _display_add( "-" );
      _display_add( time()-stats.tstart );
    }
    _display_add( "REFINE" );
    _display_flush( os );

    if( options.MAXITER && iter >= options.MAXITER ){
      _display_add( "INTERRUPT" );
      _display_flush( os );
      break;
    }

    // Refine relaxation via additional breakpoints
    _refine_polrelax( P.data(), tvar );

    // Apply domain contraction
    // Do NOT test for infeasibility here, because contraction problem may
    // become infeasible due to round-off in LP solver
    contract( P.data(), nred, tvar, !_p_inc.empty()? &_f_inc: 0, false );

    //return 0;
  }

  // Final display
  _display_final( iter );
  _display_flush( os );

  return get_status();
}

template <typename T> inline void
NLGO_GUROBI<T>::_refine_polrelax
( const T*P, const unsigned*tvar, const bool feastest )
{

  // Update polyhedral main variables
  stats.tPOLIMG -= time();
  for( auto itv = _POLenv.Vars().begin(); itv!=_POLenv.Vars().end(); ++itv ){
    double Xval = _get_variable( *itv->second );
    itv->second->add_breakpt( Xval );
#ifdef MC__NLGO_GUROBI_SHOW_BREAKPTS
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
  _set_lincuts( feastest );
  if( options.RELMETH==Options::DRL || options.RELMETH==Options::HYBRID )
    _set_poacuts( feastest );
  if( options.RELMETH==Options::CHEB || options.RELMETH==Options::HYBRID )
    _set_chebcuts( feastest );
  stats.tPOLIMG += time();

  stats.tLPSET -= time();
  _set_LPcuts();
  stats.tLPSET += time();
}

template <typename T> inline void
NLGO_GUROBI<T>::_set_polrelax
( const T*P, const unsigned*tvar, const bool feastest )
{
  // Reset polynomial image
  stats.tPOLIMG -= time();
  _POLenv.reset();
  _POLvar.clear();
  _POLenv.options = options.POLIMG;

  // Set polyhedral main variables and auxiliary objective
  auto itv = _var.begin();
  for( unsigned ip=0; itv!=_var.end(); ++itv, ip++ )
    _POLvar.push_back(
      PolVar<T>( &_POLenv, *itv, P[ip], tvar?(tvar[ip]?false:true):true )
    );
  _POLobjaux.set( &_POLenv, _objvar, T(-SBB<T>::options.INF,SBB<T>::options.INF), true );

  _set_lincuts( feastest );

  if( options.RELMETH==Options::DRL || options.RELMETH==Options::HYBRID )
    // Add polyhedral cuts
    _set_poacuts( feastest );

  if( options.RELMETH==Options::CHEB || options.RELMETH==Options::HYBRID ){
    // Chebyshev model environment reset
    if( _CMenv && (_CMenv->nvar() != _CMdim || _CMenv->nord() != options.CMODPROP) ){
      _CMvar.clear();
      delete _CMenv; _CMenv = 0;   
    }

    if( !_CMenv ){
      // Set Chebyshev model
      _CMenv = new CModel<T>( _CMdim, options.CMODPROP, options.CMODSPAR );
      _CMenv->options = options.CMODEL;
      // Set Chebyshev variables
      for( unsigned ip=0, iCM=0; ip<_nvar; ip++ )
        if( _CMexcl.find(ip) == _CMexcl.end() )
          _CMvar.push_back( CVar<T>( _CMenv, iCM++, P[ip] ) );
        else
          _CMvar.push_back( CVar<T>( P[ip] ) );
      //std::cout << "reset CM\n";
    }
    else if( !_CMvar.empty() ){
      // Update Chebyshev variables
      for( unsigned ip=0, iCM=0; ip<_nvar; ip++ )
        if( _CMexcl.find(ip) == _CMexcl.end() )
          _CMvar[ip].set( _CMenv, iCM++, P[ip] );
        else
          _CMvar[ip] = P[ip];
    }

    // Set polyhedral main scaled variables and size Chebyshev basis
    _POLscalvar.clear();
    itv = _scalvar.begin();
    for( unsigned ip=0; itv!=_scalvar.end(); ++itv, ip++ )
      _POLscalvar.push_back(
        PolVar<T>( &_POLenv, *itv, Op<T>::zeroone()*2.-1., true )
      );
    _basisord = (!options.CMODCUTS || options.CMODCUTS>options.CMODPROP)?
                options.CMODPROP: options.CMODCUTS;
    _basisdim = _CMenv->posord()[_basisord+1];
    _basis.resize( _basisdim );
    _POLbasis.resize( _basisdim );

    // Add polyhedral cuts
    _set_chebcuts( feastest );
  }
  stats.tPOLIMG += time();

  // Update cuts in GUROBI
  stats.tLPSET -= time();
  _set_LPcuts();
  stats.tLPSET += time();
}

template <typename T> inline void
NLGO_GUROBI<T>::_update_polrelax
( const T*P, const unsigned*tvar, const bool feastest )
{
  // Reset polyhedral cuts
  stats.tPOLIMG -= time();
  _POLenv.reset_cuts();

  // Update variable bounds
  unsigned ip = 0;
  for( auto itv = _POLvar.begin(); itv!=_POLvar.end(); ++itv, ip++ )
    itv->update( P[ip], tvar?(tvar[ip]?false:true):true );

  _set_lincuts( feastest );

  if( options.RELMETH==Options::DRL || options.RELMETH==Options::HYBRID )
    // Add polyhedral cuts
    _set_poacuts( feastest );

  if( options.RELMETH==Options::CHEB || options.RELMETH==Options::HYBRID ){
    // Update Chebyshev variables
    for( unsigned ip=0, iCM=0; ip<_nvar; ip++ )
      if( _CMexcl.find(ip) == _CMexcl.end() )
        _CMvar[ip].set( _CMenv, iCM++, P[ip] );
      else
        _CMvar[ip] = P[ip];

    // Add Chebyshev-derived polyhedral cuts
    _set_chebcuts( feastest );
  }
  stats.tPOLIMG += time();

  // Update cuts in GUROBI
  stats.tLPSET -= time();
  _set_LPcuts();
  stats.tLPSET += time();
}

template <typename T> inline void
NLGO_GUROBI<T>::_set_lincuts
( const bool feastest )
{
  // Add polyhedral cuts for objective - by-pass if feasibility test or no objective function defined
  auto ito = std::get<0>(_obj).begin();
  for( unsigned j=0; ito!=std::get<0>(_obj).end() && j<1; ++ito, j++ ){
    if( feastest || _fct_lin.find(std::get<1>(_obj).data()) == _fct_lin.end() ) continue;
    try{
      _dag->eval( _op_f, _op_POLfg.data(), 1, std::get<1>(_obj).data(), &_POLobj, _nvar, _var.data(), _POLvar.data() );
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
      _dag->eval( _op_g[j], _op_POLfg.data(), 1, std::get<1>(_ctr).data()+j, _POLctr.data()+j, _nvar, _var.data(), _POLvar.data() );
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
#ifdef MC__NLGO_GUROBI_DEBUG
  std::cout << _POLenv;
#endif
}

template <typename T> inline void
NLGO_GUROBI<T>::_set_poacuts
( const bool feastest )
{
  // Add polyhedral cuts for objective - by-pass if feasibility test or no objective function defined
  //if( !feastest && std::get<1>(_obj).size() ){
  auto ito = std::get<0>(_obj).begin();
  for( unsigned j=0; ito!=std::get<0>(_obj).end() && j<1; ++ito, j++ ){
    if( feastest || _fct_lin.find(std::get<1>(_obj).data()) != _fct_lin.end() ) continue;
    try{
      //std::cout << "objective" << std::endl;
      _dag->eval( _op_f, _op_POLfg.data(), 1, std::get<1>(_obj).data(), &_POLobj, _nvar, _var.data(), _POLvar.data() );
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
      _dag->eval( _op_g[j], _op_POLfg.data(), 1, std::get<1>(_ctr).data()+j, _POLctr.data()+j, _nvar, _var.data(), _POLvar.data() );
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
#ifdef MC__NLGO_GUROBI_DEBUG
  std::cout << _POLenv;
#endif
}

template <typename T> inline void
NLGO_GUROBI<T>::_set_chebcuts
( const bool feastest )
{
  // Add cuts for variable scaling and basis polynomials (up to order options.CMODCUT)
  auto itv = _POLscalvar.begin();
  for( unsigned i=0; itv!=_POLscalvar.end(); ++itv, i++ ){
    _POLenv.add_cut( PolCut<T>::EQ, _CMenv->refvar()[i], _POLvar[i], 1., *itv, -_CMenv->scalvar()[i] );
#ifdef MC__NLGO_GUROBI_DEBUG
      std::cout << "Variable #" << i << ": " << _CMvar[i];
#endif
  }
  bool basis_dense = false;

  // Add Chebyshev-derived cuts for objective - by-pass if feasibility test or no objective function defined
  //if( !feastest && std::get<1>(_obj).size() ){
  auto ito = std::get<0>(_obj).begin();
  for( unsigned j=0; ito!=std::get<0>(_obj).end() && j<1; ++ito, j++ ){
    if( feastest || _fct_lin.find(std::get<1>(_obj).data()) != _fct_lin.end() ) continue;
    try{
       // Polynomial model evaluation
      _dag->eval( _op_f, _op_CMfg.data(), 1, std::get<1>(_obj).data(), &_CMobj,
                  _nvar, _var.data(), _CMvar.data() );
#ifdef MC__NLGO_GUROBI_DEBUG
       std::cout << _CMobj;
       if(_p_inc.data()) std::cout << _CMobj.P(_p_inc.data())+_CMobj.R() << std::endl;
#endif

      // Polynomial model remainder
      T Robj = _CMobj.remainder();
      for( unsigned i=_basisord+1; i<=_CMenv->nord(); i++ )
        Robj += _CMobj.bndord(i);

      // Dense polynomial relaxation
      if( _CMobj.ndxmon().empty() ){
        if( !basis_dense ){
          _CMenv->get_bndmon( _basisord, _basis.data(), _scalvar.data(), true );
          _dag->eval( _basis.size()-1, _basis.data()+1, _POLbasis.data()+1, _scalvar.size(),
                      _scalvar.data(), _POLscalvar.data() );
          _POLenv.generate_cuts( _basis.size()-1, _POLbasis.data()+1, false );
          basis_dense = true;
        }
        switch( std::get<0>(_obj)[0] ){
          case MIN: _POLenv.add_cut( PolCut<T>::LE, -Op<T>::l(Robj)-_CMobj.coefmon().second[0],
            _POLbasis.size()-1, _POLbasis.data()+1, _CMobj.coefmon().second+1, _POLobjaux, -1. ); break;
          case MAX: _POLenv.add_cut( PolCut<T>::GE, -Op<T>::u(Robj)-_CMobj.coefmon().second[0],
            _POLbasis.size()-1, _POLbasis.data()+1, _CMobj.coefmon().second+1, _POLobjaux, -1. ); break;
        }
      }

      // Sparse polynomial relaxation
      else{
        auto first_CMobj = _CMobj.ndxmon().cbegin();
        const double a0 = *first_CMobj? 0.: (++first_CMobj, _CMobj.coefmon().second[0]);
        auto last_CMobj = _CMobj.ndxmon().lower_bound(_basisdim);
        const std::set<unsigned> ndx_CMobj( first_CMobj, last_CMobj );
        if( !basis_dense ){
          _CMenv->get_bndmon( _basisord, _basis.data(), ndx_CMobj, _scalvar.data(), true );
          _dag->eval( ndx_CMobj, _basis.data(), _POLbasis.data(), _scalvar.size(),
                      _scalvar.data(), _POLscalvar.data() );
          _POLenv.generate_cuts( ndx_CMobj, _POLbasis.data(), false );
        }
        switch( std::get<0>(_obj)[0] ){
          case MIN: _POLenv.add_cut( PolCut<T>::LE, -Op<T>::l(Robj)-a0, ndx_CMobj,
            _POLbasis.data(), _CMobj.coefmon().second, _POLobjaux, -1. ); break;
          case MAX: _POLenv.add_cut( PolCut<T>::GE, -Op<T>::u(Robj)-a0, ndx_CMobj,
            _POLbasis.data(), _CMobj.coefmon().second, _POLobjaux, -1. ); break;
        }
      }
    }
    catch(...){
      _CMobj = T(-SBB<T>::options.INF,SBB<T>::options.INF);
      // No cut added for objective in case evaluation failed
    }
  }
#ifdef MC__NLGO_GUROBI_DEBUG
  std::cout << _POLenv;
#endif

  // Add Chebyshev-derived cuts for constraints - add slack to inequality constraints if feasibility test
  _CMctr.resize( _nctr );
  auto itc = std::get<0>(_ctr).begin();
  for( unsigned j=0; itc!=std::get<0>(_ctr).end(); ++itc, j++ ){
    if( _fct_lin.find(std::get<1>(_ctr).data()+j) != _fct_lin.end() ) continue;
    try{
      // Polynomial model evaluation
      _dag->eval( _op_g[j], _op_CMfg.data(), 1, std::get<1>(_ctr).data()+j, _CMctr.data()+j,
                  _nvar, _var.data(), _CMvar.data() );
#ifdef MC__NLGO_GUROBI_DEBUG
      std::cout << "Constraint #" << j << ": " << _CMctr[j];
      if(_p_inc.data()) std::cout << "optim: " << _CMctr[j].P(_p_inc.data())+_CMctr[j].R() << std::endl;
#endif
      // Test for too large Chebyshev bounds or NaN
      if( !(Op<CVar<T>>::diam(_CMctr[j]) <= options.CMODDMAX) ) throw(0);

      // Polynomial model remainder
      T Rctr = _CMctr[j].remainder();
      for( unsigned i=_basisord+1; i<=_CMenv->nord(); i++ )
        Rctr += _CMctr[j].bndord(i);

      // Dense polynomial relaxation
      if( _CMctr[j].ndxmon().empty() ){
        if( !basis_dense ){
          _CMenv->get_bndmon( _basisord, _basis.data(), _scalvar.data(), true );
          _dag->eval( _basis.size()-1, _basis.data()+1, _POLbasis.data()+1, _scalvar.size(),
                      _scalvar.data(), _POLscalvar.data() );
          _POLenv.generate_cuts( _basis.size()-1, _POLbasis.data()+1, false );
          basis_dense = true;
        }
        switch( (*itc) ){
          case EQ: if( !feastest ) _POLenv.add_cut( PolCut<T>::GE,
            -Op<T>::u(Rctr)-_CMctr[j].coefmon().second[0],
            _POLbasis.size()-1, _POLbasis.data()+1, _CMctr[j].coefmon().second+1 ); // no break
                   else            _POLenv.add_cut( PolCut<T>::GE,
            -Op<T>::u(Rctr)-_CMctr[j].coefmon().second[0],
            _POLbasis.size()-1, _POLbasis.data()+1, _CMctr[j].coefmon().second+1, _POLobjaux,  1. ); // no break
          case LE: if( !feastest ){ _POLenv.add_cut( PolCut<T>::LE,
            -Op<T>::l(Rctr)-_CMctr[j].coefmon().second[0],
            _POLbasis.size()-1, _POLbasis.data()+1, _CMctr[j].coefmon().second+1 ); break; }
                   else           { _POLenv.add_cut( PolCut<T>::LE,
            -Op<T>::l(Rctr)-_CMctr[j].coefmon().second[0],
            _POLbasis.size()-1, _POLbasis.data()+1, _CMctr[j].coefmon().second+1, _POLobjaux, -1. ); break; }
          case GE: if( !feastest ){ _POLenv.add_cut( PolCut<T>::GE,
            -Op<T>::u(Rctr)-_CMctr[j].coefmon().second[0],
            _POLbasis.size()-1, _POLbasis.data()+1, _CMctr[j].coefmon().second+1 ); break; }
                   else           { _POLenv.add_cut( PolCut<T>::GE,
            -Op<T>::u(Rctr)-_CMctr[j].coefmon().second[0],
            _POLbasis.size()-1, _POLbasis.data()+1, _CMctr[j].coefmon().second+1, _POLobjaux,  1. ); break; }
        }
      }

      // Sparse polynomial relaxation
      else{
        auto first_CMctr = _CMctr[j].ndxmon().cbegin();
        const double a0 = *first_CMctr? 0.: (++first_CMctr, _CMctr[j].coefmon().second[0]);
        auto last_CMctr = _CMctr[j].ndxmon().lower_bound(_basisdim);
        const std::set<unsigned> ndx_CMctr( first_CMctr, last_CMctr );
        if( !basis_dense ){
          _CMenv->get_bndmon( _basisord, _basis.data(), ndx_CMctr, _scalvar.data(), true );
          _dag->eval( ndx_CMctr, _basis.data(), _POLbasis.data(), _scalvar.size(),
                      _scalvar.data(), _POLscalvar.data() );
          _POLenv.generate_cuts( ndx_CMctr, _POLbasis.data(), false );
        }
        switch( (*itc) ){
          case EQ: if( !feastest ) _POLenv.add_cut( PolCut<T>::GE, -Op<T>::u(Rctr)-a0,
            ndx_CMctr, _POLbasis.data(), _CMctr[j].coefmon().second ); // no break
                   else            _POLenv.add_cut( PolCut<T>::GE, -Op<T>::u(Rctr)-a0,
            ndx_CMctr, _POLbasis.data(), _CMctr[j].coefmon().second, _POLobjaux,  1. ); // no break
          case LE: if( !feastest ){ _POLenv.add_cut( PolCut<T>::LE, -Op<T>::l(Rctr)-a0,
            ndx_CMctr, _POLbasis.data(), _CMctr[j].coefmon().second ); break; }
                   else           { _POLenv.add_cut( PolCut<T>::LE, -Op<T>::l(Rctr)-a0,
            ndx_CMctr, _POLbasis.data(), _CMctr[j].coefmon().second, _POLobjaux, -1. ); break; }
          case GE: if( !feastest ){ _POLenv.add_cut( PolCut<T>::GE, -Op<T>::u(Rctr)-a0,
            ndx_CMctr, _POLbasis.data(), _CMctr[j].coefmon().second ); break; }
                   else           { _POLenv.add_cut( PolCut<T>::GE, -Op<T>::u(Rctr)-a0,
            ndx_CMctr, _POLbasis.data(), _CMctr[j].coefmon().second, _POLobjaux,  1. ); break; }
        }
      }
    }
    catch(...){
#ifdef MC__NLGO_GUROBI_DEBUG
      std::cout << "Constraint #" << j << ": " << "bad cut!\n";
#endif
      _CMctr[j] = T(-SBB<T>::options.INF,SBB<T>::options.INF);
      // No cut added for constraint #j in case evaluation failed
      continue;
    }
  }
#ifdef MC__NLGO_GUROBI_DEBUG
  { int dum; std::cin >> dum; }
  //std::cout << _POLenv;
#endif
}

template <typename T> inline void
NLGO_GUROBI<T>::_set_SBBoptions
()
{
  // SBB options
  SBB<T>::options.ABSOLUTE_TOLERANCE = options.CVATOL;
  SBB<T>::options.RELATIVE_TOLERANCE = options.CVRTOL;
  SBB<T>::options.MAX_NODES          = options.MAXITER;
  SBB<T>::options.MAX_CPU_TIME       = options.MAXCPU;
  SBB<T>::options.DISPLAY            = options.DISPLAY;

  // TREE_REREDUCE(true)
  // BRANCHING_STRATEGY(OMEGA)
  // BRANCHING_VARIABLE_CRITERION(RGREL)
  // BRANCHING_RELIABILITY_THRESHOLD(1)
  // BRANCHING_SCORE_FACTOR(1./6.)
  // BRANCHING_BOUND_THRESHOLD(1e-3)
  // INF(1e20)
}

template <typename T> inline void
NLGO_GUROBI<T>::_solve_LPmodel
()
{
  stats.tLPSOL -= time();
  _set_LPoptions();
#ifdef MC__USE_CPLEX
  if( options.MIPFILE != "" ) _ILOcplex->exportModel( options.MIPFILE.c_str() );
  //_time = -_cplex->getCplexTime();
  _ILOcplex->solve();
  //_time += _cplex->getCplexTime();
#else
  _GRBmodel->update();
  if( options.MIPFILE != "" ) _GRBmodel->write( options.MIPFILE );
  fedisableexcept(FE_ALL_EXCEPT);
  _GRBmodel->optimize();
#endif
  stats.tLPSOL += time();
  stats.nLPSOL++;
#ifdef MC__NLGO_GUROBI_DEBUG
  std::cout << "LP solution complete\n";
  std::cout << std::scientific << std::setprecision(4);
  std::cout << "  fopt = " << get_objective() << std::endl;
  unsigned ivar = 0;
  for( auto itv=_var.begin(); itv!=_var.end(); ++itv, ivar++ )
    std::cout << "popt" << ivar << " = " << get_variable( *itv ) << std::endl;
  { int dum; std::cin >> dum; }
#endif
}

template <typename T> inline void
NLGO_GUROBI<T>::_set_LPoptions
()
{
#ifdef MC__USE_CPLEX
  // CPLEX options
  _ILOcplex->extract(*_ILOmodel);
  _ILOcplex->setWarning( options.MIPDISPLAY? std::cout: _ILOenv->getNullStream() );
  _ILOcplex->setOut( options.MIPDISPLAY? std::cout: _ILOenv->getNullStream() );
  _ILOcplex->setParam( IloCplex::RootAlg, options.LPALGO );
  _ILOcplex->setParam( IloCplex::EpOpt,   options.LPOPTIMTOL );
  _ILOcplex->setParam( IloCplex::EpRHS,   options.LPFEASTOL );
#else
  // Gurobi options
  _GRBmodel->getEnv().set( GRB_IntParam_Method,            options.LPALGO );
  _GRBmodel->getEnv().set( GRB_DoubleParam_FeasibilityTol, options.LPFEASTOL );
  _GRBmodel->getEnv().set( GRB_DoubleParam_OptimalityTol,  options.LPOPTIMTOL );
  _GRBmodel->getEnv().set( GRB_DoubleParam_MIPGap,         options.MIPRELGAP );
  _GRBmodel->getEnv().set( GRB_DoubleParam_MIPGapAbs,      options.MIPABSGAP );
  _GRBmodel->getEnv().set( GRB_IntParam_Presolve,          options.LPPRESOLVE  );
  _GRBmodel->getEnv().set( GRB_IntParam_OutputFlag,        options.MIPDISPLAY );
  _GRBmodel->getEnv().set( GRB_DoubleParam_PreSOS2BigM,    options.PRESOS2BIGM );
  _GRBmodel->getEnv().set( GRB_IntParam_DualReductions,    0 ); // In order to avoid INF_OR_UNBD status
#endif
}

template <typename T> inline void
NLGO_GUROBI<T>::_set_LPcuts
()
{
  _reset_LPmodel();
  for( auto itv=_POLenv.Vars().begin(); itv!=_POLenv.Vars().end(); ++itv )
    _set_LPvar( itv->second );
    //if( itv->second->cuts() ) _set_LPvar( itv->second );
  for( auto itv=_POLenv.Aux().begin(); itv!=_POLenv.Aux().end(); ++itv )
    _set_LPvar( *itv );
#ifndef MC__USE_CPLEX
  _GRBmodel->update();
#endif
  for( auto itc=_POLenv.Cuts().begin(); itc!=_POLenv.Cuts().end(); ++itc )
    _set_LPcut( *itc );
}

template <typename T> inline void
NLGO_GUROBI<T>::_reset_LPmodel
()
{
#ifdef MC__USE_CPLEX
  delete _ILOmodel; delete _ILOcplex;
  _ILOenv->end();
  delete _ILOenv;
  _ILOenv = new IloEnv;
  _ILOmodel = new IloModel( *_ILOenv );
  _ILOcplex = new IloCplex( *_ILOenv );
  _LPobj.first = false;
#else
  delete _GRBmodel;
  _GRBmodel = new GRBModel( *_GRBenv );
#endif
  _LPvar.clear(); _LPcut.clear(); _LPinc.first = false;
}

template <typename T> inline void
NLGO_GUROBI<T>::_set_LPrelax
( const bool feastest )
{
  // Set objective
  auto jtobj = _LPvar.find( &_POLobjaux );
#ifdef MC__USE_CPLEX
  if( _LPobj.first ) _ILOmodel->remove( _LPobj.second );
  if( !feastest )
    switch( std::get<0>(_obj)[0] ){
      case MIN: _LPobj.second = IloMinimize( *_ILOenv, jtobj->second ); break;
      case MAX: _LPobj.second = IloMaximize( *_ILOenv, jtobj->second ); break;
    }
  else
    _LPobj.second = IloMinimize( *_ILOenv, jtobj->second );
  _ILOmodel->add( _LPobj.second ); _LPobj.first = true;
#else
  _GRBmodel->setObjective( GRBLinExpr( jtobj->second,  1. ) );
  if( !feastest )
    switch( std::get<0>(_obj)[0] ){
      case MIN: _GRBmodel->set( GRB_IntAttr_ModelSense,  1 ); break;
      case MAX: _GRBmodel->set( GRB_IntAttr_ModelSense, -1 ); break;
    }
  else
    _GRBmodel->set( GRB_IntAttr_ModelSense,  1 );
#endif

  // Remove incumbent constraint (if any)
  if( _LPinc.first ){
#ifdef MC__USE_CPLEX
    _ILOmodel->remove( _LPinc.second );
#else
    _GRBmodel->remove( _LPinc.second );
#endif
    _LPinc.first = false;
  }
}

template <typename T> inline void
NLGO_GUROBI<T>::_set_LPcontract
( const unsigned ip, const bool uplo, const double*inc,
  const bool feastest )
{
  // Add lower/upper parameter bound objective
  auto jtobj = _LPvar.find( &_POLvar[ip] );
#ifdef MC__USE_CPLEX
  if( _LPobj.first ) _ILOmodel->remove( _LPobj.second );
  switch( uplo ){
    case false: _LPobj.second = IloMinimize( *_ILOenv, jtobj->second ); break;
    case true:  _LPobj.second = IloMaximize( *_ILOenv, jtobj->second ); break;
  }
  _ILOmodel->add( _LPobj.second ); _LPobj.first = true;
#else
  _GRBmodel->setObjective( GRBLinExpr( jtobj->second,  1. ) );
  switch( uplo ){
    case false: _GRBmodel->set( GRB_IntAttr_ModelSense,  1 ); break;
    case true:  _GRBmodel->set( GRB_IntAttr_ModelSense, -1 ); break;
  }
#endif

  // Add/update incumbent or backoff constraint
  if( _LPinc.first ){
#ifdef MC__USE_CPLEX
    _ILOmodel->remove( _LPinc.second );
#else
    _GRBmodel->remove( _LPinc.second );
#endif
    _LPinc.first = false;
  }
  if( !inc ) return;
  jtobj = _LPvar.find( &_POLobjaux );
#ifdef MC__USE_CPLEX
  IloExpr lhsinc( *_ILOenv ); lhsinc = jtobj->second;
#endif
  if( !feastest ){
    switch( std::get<0>(_obj)[0] ){
#ifdef MC__USE_CPLEX
      case MIN: _LPinc.second = (lhsinc <= *inc); _ILOmodel->add( _LPinc.second ); break;
      case MAX: _LPinc.second = (lhsinc >= *inc); _ILOmodel->add( _LPinc.second ); break;
#else
      case MIN: _LPinc.second = _GRBmodel->addConstr( GRBLinExpr( jtobj->second,  1. ),
        GRB_LESS_EQUAL, *inc ); break;
      case MAX: _LPinc.second = _GRBmodel->addConstr( GRBLinExpr( jtobj->second,  1. ),
        GRB_GREATER_EQUAL, *inc ); break;
#endif
    }
  }
  else
#ifdef MC__USE_CPLEX
    _LPinc.second = (lhsinc <= *inc); _ILOmodel->add( _LPinc.second );
#else
    _LPinc.second = _GRBmodel->addConstr( GRBLinExpr( jtobj->second,  1. ),
        GRB_LESS_EQUAL, *inc );
#endif
  _LPinc.first = true;
}

template <typename T> inline std::pair<typename NLGO_GUROBI<T>::t_LPVar::iterator,bool>
NLGO_GUROBI<T>::_set_LPvar
( const PolVar<T>*pVar )
{
#ifdef MC__USE_CPLEX
  IloNumVar var;
#else
  GRBVar var;
#endif
  switch( pVar->id().first ){
    case PolVar<T>::VARCONT:
    case PolVar<T>::AUXCONT:
#ifdef MC__USE_CPLEX
      var = IloNumVar( *_ILOenv, Op<T>::l(pVar->range()), Op<T>::u(pVar->range()),
       ILOFLOAT, pVar->name().c_str() );
      _ILOmodel->add( var );
#else
      var = _GRBmodel->addVar( Op<T>::l(pVar->range()), Op<T>::u(pVar->range()),
        0., GRB_CONTINUOUS, pVar->name() );
#endif
      break;
    case PolVar<T>::VARINT:
    case PolVar<T>::AUXINT:
      if( isequal( Op<T>::l(pVar->range()), 0. ) && isequal( Op<T>::u(pVar->range()), 1. ) ){
#ifdef MC__USE_CPLEX
        var = IloNumVar( *_ILOenv, 0., 1., ILOBOOL, pVar->name().c_str() );
        _ILOmodel->add( var );
#else
        var = _GRBmodel->addVar( 0., 1., 0., GRB_BINARY, pVar->name() );
#endif
      }
      else{
#ifdef MC__USE_CPLEX
        var = IloNumVar( *_ILOenv, Op<T>::l(pVar->range()), Op<T>::u(pVar->range()),
         ILOINT, pVar->name().c_str() );
        _ILOmodel->add( var );
#else
        var = _GRBmodel->addVar( Op<T>::l(pVar->range()), Op<T>::u(pVar->range()),
          0., GRB_INTEGER, pVar->name() );
#endif
      }
      break;

    default:
      throw std::runtime_error("Invalid auxiliary variable type");
  }

  return _LPvar.insert( std::make_pair( pVar, var ) );
}

template <typename T> inline void
NLGO_GUROBI<T>::_set_LPcut
( const PolCut<T>*pCut )
{
#ifdef MC__USE_CPLEX
  bool isSOS = ( pCut->type() == PolCut<T>::SOS1
              || pCut->type() == PolCut<T>::SOS2 );
  IloNumVarArray VarSOS( *_ILOenv, (isSOS? pCut->nvar(): 0) );
  IloNumArray WeiSOS( *_ILOenv, (isSOS? pCut->nvar(): 0) );

  IloExpr lhs( *_ILOenv );
  for( unsigned k=0; k<pCut->nvar(); k++ ){
    auto ivar = _LPvar.find( pCut->var()+k );
    if( ivar==_LPvar.end() ) throw std::runtime_error("variable not found");
    if( pCut->type() == PolCut<T>::SOS1 || pCut->type() == PolCut<T>::SOS2 ){
      VarSOS[k] = ivar->second; WeiSOS[k] = pCut->coef()[k];
    }
    else
      lhs += pCut->coef()[k] * ivar->second;
  }

  IloRange ctr;
  try{
    switch( pCut->type() ){
      case PolCut<T>::SOS1:
        _ILOmodel->add( IloSOS1( *_ILOenv, VarSOS, WeiSOS ) );
        ctr = ( lhs == pCut->rhs() ); break;
      case PolCut<T>::SOS2:
        _ILOmodel->add( IloSOS2( *_ILOenv, VarSOS, WeiSOS ) );
        ctr = ( lhs == pCut->rhs() ); break;
      case PolCut<T>::EQ:
        ctr = ( lhs == pCut->rhs() ); break;
      case PolCut<T>::LE:
        ctr = ( lhs <= pCut->rhs() ); break;
      case PolCut<T>::GE:
        ctr = ( lhs >= pCut->rhs() ); break;
    }
    _ILOmodel->add( ctr );
    _LPcut.insert( std::make_pair( pCut, ctr ) );
  }
  catch(IloException& e){
#ifdef MC__NLGO_GUROBI_DEBUG
    std::cout << "IloException Caught - Error code = " << e.getMessage() << std::endl;
#endif
  }

#else
  GRBVar* VarSOS = 0;
  double* WeiSOS = 0;
  int TypSOS = 0;
  if( pCut->type() == PolCut<T>::SOS1 || pCut->type() == PolCut<T>::SOS2 ){
    VarSOS = new GRBVar[pCut->nvar()];
    WeiSOS = new double[pCut->nvar()];
    TypSOS = ( pCut->type() == PolCut<T>::SOS1? GRB_SOS_TYPE1: GRB_SOS_TYPE2 );
  }

  GRBLinExpr lhs;
  for( unsigned k=0; k<pCut->nvar(); k++ ){
    auto ivar = _LPvar.find( pCut->var()+k );
    if( ivar==_LPvar.end() ) throw std::runtime_error("variable not found");
    GRBVar&Var = ivar->second;
    if( pCut->type() == PolCut<T>::SOS1 || pCut->type() == PolCut<T>::SOS2 ){
      VarSOS[k] = Var; WeiSOS[k] = pCut->coef()[k];//(double)k/(double)pCut->nvar();//
    }
    else
      lhs += GRBLinExpr( Var, pCut->coef()[k] );
  }

  GRBConstr ctr;
  try{
    switch( pCut->type() ){
      case PolCut<T>::SOS1:
      case PolCut<T>::SOS2:
        _GRBmodel->addSOS( VarSOS, WeiSOS, pCut->nvar(), TypSOS );
        delete [] VarSOS;
        delete [] WeiSOS;
        break;
      case PolCut<T>::EQ:
        ctr = _GRBmodel->addConstr( lhs, GRB_EQUAL, pCut->rhs() );
        break;
      case PolCut<T>::LE:
        ctr = _GRBmodel->addConstr( lhs, GRB_LESS_EQUAL, pCut->rhs() );
        break;
      case PolCut<T>::GE:
        ctr = _GRBmodel->addConstr( lhs, GRB_GREATER_EQUAL, pCut->rhs() );
        break;
    }
    _LPcut.insert( std::make_pair( pCut, ctr ) );
  }
  catch(GRBException& e){
#ifdef MC__NLGO_GUROBI_DEBUG
    std::cout << "GRBException Caught - Error code = " << e.getErrorCode() << std::endl;
    std::cout << e.getMessage() << std::endl;
#endif
  }
#endif
}

template <typename T> inline void
NLGO_GUROBI<T>::_display_init
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

template <typename T> inline void
NLGO_GUROBI<T>::_display_add
( const double dval )
{
  if( options.DISPLAY <= 1 ) return;
  _odisp << std::right << std::scientific << std::setprecision(_DPREC)
         << std::setw(_DPREC+8) << dval;
}

template <typename T> inline void
NLGO_GUROBI<T>::_display_add
( const int ival )
{
  if( options.DISPLAY <= 1 ) return;
  _odisp << std::right << std::setw(_IPREC) << ival;
}

template <typename T> inline void
NLGO_GUROBI<T>::_display_add
( const std::string &sval )
{
  if( options.DISPLAY <= 1 ) return;
  _odisp << std::right << std::setw(_DPREC+8) << sval;
}

template <typename T>
inline void
NLGO_GUROBI<T>::_display_flush
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
NLGO_GUROBI<T>::_display_final
( const unsigned iter )
{
  if( options.DISPLAY <= 0 ) return;
  _odisp << std::endl << "#  TERMINATION: ";
  _odisp << std::fixed << std::setprecision(6) << time()-stats.tstart << " CPU SEC"
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
  if( !optimal_status() )
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
NLGO_GUROBI<T>::Options::display
( std::ostream&out ) const
{
  // Display NLGO Options
  out << std::left;
  out << std::setw(60) << "  COMPLETE SEARCH METHOD";
  switch( CSALGO ){
   case SBB: out << "SBB" << std::endl; break;
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
( std::ostream&out, const NLGO_GUROBI<T>&NLP )
{
  out << std::endl
      << std::setfill('_') << std::setw(72) << "\n\n" << std::setfill(' ')
      << std::setw(55) << "NONLINEAR GLOBAL OPTIMIZATION IN CRONOS\n"
      << std::setfill('_') << std::setw(72) << "\n\n" << std::setfill(' ');

  // Display NLGO_GUROBI Options
  //out << std::left << "SPATIAL BRANCH-AND-BOUND OPTIONS:\n\n";
  //NLP.SBB<T>::options.display( out );
  NLP.NLGO_GUROBI<T>::options.display( out );

  out << std::setfill('_') << std::setw(72) << " " << std::setfill(' ')
      << std::endl << std::endl;
  return out;
}

} // end namescape mc

#endif
