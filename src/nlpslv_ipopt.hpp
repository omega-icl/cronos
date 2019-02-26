// Copyright (C) 2014-2016 Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

/*!
\page page_NLPSLV_IPOPT Local (Continuous) Optimization using IPOPT
\author Benoit C. Chachuat
\version 1.0
\date 2014
\bug No known bugs.

Consider a nonlinear optimization problem in the form:
\f{align*}
\mathcal{P}:\quad & \min_{x_1,\ldots,x_n}\ f(x_1,\ldots,x_n)\\
& {\rm s.t.}\ \ g_j(x_1,\ldots,x_n)\ \leq,=,\geq\ 0,\ \ j=1,\ldots,m\\
& \qquad x_i^L\leq x_i\leq x_i^U,\ \ i=1,\ldots,n\,,
\f}
where \f$f, g_1, \ldots, g_m\f$ are factorable, potentially nonlinear, real-valued functions; and \f$x_1, \ldots, x_n\f$ are continuous decision variables. The class mc::NLPSLV_IPOPT solves such NLP problems using the software package <A href="https://projects.coin-or.org/Ipopt">IPOPT</A>, which implements a local solution method (interior point). IPOPT requires the first and second derivatives as well as the sparsity pattern of the objective and constraint functions in the NLP model. This information is generated using direct acyclic graphs (DAG) in <A href="https://projects.coin-or.org/MCpp">MC++</A>.

\section sec_NLPSLV_solve How to Solve an NLP Model using mc::NLPSLV_IPOPT?

Consider the following NLP:
\f{align*}
  \max_{\bf p}\ & p_1+p_2 \\
  \text{s.t.} \ & p_1\,p_2 \leq 4 \\
  & 0 \leq p_1 \leq 6\\
  & 0 \leq p_2 \leq 4\,.
\f}

We start by defining an mc::NLPSLV_IPOPT class as below, whereby the class <A href="http://www.coin-or.org/Doxygen/Ipopt/class_ipopt_1_1_smart_ptr.html">Ipopt::SmartPtr</A> is used here to comply with the IPOPT C++ interface:

\code
  #include "nlpslv_ipopt.hpp"
  
  Ipopt::SmartPtr<mc::NLPSLV_IPOPT> NLP = new mc::NLPSLV_IPOPT;
\endcode

Next, we set the variables and objective/constraint functions by creating a DAG of the problem: 

\code
  mc::FFGraph DAG;
  const unsigned NP = 2; mc::FFVar P[NP];
  for( unsigned i=0; i<NP; i++ ) P[i].set( &DAG );
  NLP->set_dag( &DAG );                       // DAG
  NLP->set_var( NP, P );                      // decision variables
  NLP->set_obj( BASE_NLP::MAX, P[0]+P[1] );   // objective
  NLP->add_ctr( BASE_NLP::LE, P[0]*P[1]-4. ); // constraints
  NLP->setup();                               // DAG
\endcode

Given initial bounds \f$P\f$ and initial guesses \f$p_0\f$ on the decision variables, the NLP model is solved as:

\code
  #include "interval.hpp"

  typedef mc::interval I;
  I P[NP] = { I(0.,6.), I(0.,4.) };
  double p0[NP] = { 5., 1. };

  Ipopt::ApplicationReturnStatus status = NLP->solve( P, p0 );
\endcode

The return value is of the enumeration type <A href="http://www.coin-or.org/Doxygen/Ipopt/namespace_ipopt.html#efa0497854479cde8b0994cdf132c982">Ipopt::ApplicationReturnStatus</A>. Moreover, the optimization results can be retrieved as follows:

\code
  #include <iostream>
  #include <iomanip>

  if( status == Ipopt::Solve_Succeeded ){
    std::cout << "NLP (LOCAL) SOLUTION: " << std::endl;
    std::cout << "  f* = " << pNLP->solution().f << std::endl;
    for( unsigned int ip=0; ip<NP; ip++ )
      std::cout << "  p*(" << ip << ") = " << pNLP->solution().p[ip] << std::endl;
  }
\endcode

The following result is displayed here:

\verbatim
NLP (LOCAL) SOLUTION: 
  f* = 6.66667
  p*(0) = 6
  p*(1) = 0.666667
\endverbatim

Regarding options, the output level, maximum number of iterations, tolerance, maximum CPU time, etc can all be modified through the public member mc::NLPSLV_IPOPT::options. 
*/

#ifndef MC__NLPSLV_IPOPT_HPP
#define MC__NLPSLV_IPOPT_HPP

#include <stdexcept>
#include <cassert>
#include "coin/IpTNLP.hpp"
#include "coin/IpIpoptApplication.hpp"

#include "mclapack.hpp"
#include "base_nlp.hpp"

#ifdef MC__USE_SOBOL
  #include "sobol.hpp"
#endif

#undef MC__NLPSLV_IPOPT_DEBUG
#undef MC__NLPSLV_IPOPT_TRACE

/* TO DO:
- overload << operator for NLP solution
*/

namespace mc
{

//! @brief C++ derived class for NLP solution using IPOPT and FADBAD++
////////////////////////////////////////////////////////////////////////
//! mc::NLPSLV_IPOPT is a C++ derived class for solving NLP problems
//! using IPOPT and MC++
////////////////////////////////////////////////////////////////////////
class NLPSLV_IPOPT:
  public Ipopt::TNLP,
  public virtual BASE_NLP
{
  // Overloading stdout operator
  friend std::ostream& operator<<
    ( std::ostream&os, const NLPSLV_IPOPT& );

private:
  //! @brief total count of variables (independent and dependent) in problem
  unsigned _nvar;
  //! @brief vector of decision variables (independent and dependent)
  std::vector<FFVar> _var;

  //! @brief number of constraints in problem
  unsigned _nctr;
  //! @brief vector of constraints (main and equation system)
  std::tuple< std::vector<t_CTR>, std::vector<FFVar>, std::vector<FFVar> > _ctr;

  //! @brief list of operations for objective evaluation
  FFSubgraph _op_f;
  //! @brief list of operations for objective gradient evaluation
  FFSubgraph _op_df;
  //! @brief list of operations for constraint evaluation
  FFSubgraph _op_g;
  //! @brief list of operations for constraint gradient evaluation
  FFSubgraph _op_dg;
  //! @brief list of operations for Lagragian Hessian evaluation
  FFSubgraph _op_d2L;
  //! @brief Storage vector for function evaluation in double arithmetic
  std::vector<Ipopt::Number> _dwk;

  //! @brief objective gradient
  const FFVar* _obj_grad;
  //! @brief constraint gradient (sparse format)
  std::tuple< unsigned, const unsigned*, const unsigned*, const FFVar* > _ctr_grad;
  //! @brief Lagrangian Hessian (sparse format)
  std::tuple< unsigned, const unsigned*, const unsigned*, const FFVar* > _lagr_hess;

  //! @brief Internal scaling factor for the objective function (+1: minimize; -1: maximize)
  double _scaling;
  //! @brief Internal DAG variable for Lagrangian structure
  FFVar _lagr;

  //! @brief Structure holding optimization problem data
  struct DATA{
    const double*p0;
    std::pair<double,double>*P;
    DATA(): p0(0), P(0) {}
    ~DATA() { delete[] P; }
    template <typename T> void set
      ( const unsigned np, const double*p0_, const T*P_ )
      {
        p0 = p0_;
        delete[] P;
        P = new std::pair<double,double>[np];
        for( unsigned i=0; i<np; i++ )
          P[i] = std::make_pair( Op<T>::l(P_[i]), Op<T>::u(P_[i]) );
      }
  } _data;

  //! @brief Structure holding solution information
  SOLUTION_OPT _solution;

  //! @brief Cleanup gradient/hessian storage
  void _cleanup()
  {
    delete[] _obj_grad;               _obj_grad = 0;
    delete[] std::get<1>(_ctr_grad);  std::get<1>(_ctr_grad) = 0;
    delete[] std::get<2>(_ctr_grad);  std::get<2>(_ctr_grad) = 0;
    delete[] std::get<3>(_ctr_grad);  std::get<3>(_ctr_grad) = 0;
    delete[] std::get<1>(_lagr_hess); std::get<1>(_lagr_hess) = 0;
    delete[] std::get<2>(_lagr_hess); std::get<2>(_lagr_hess) = 0;
    delete[] std::get<3>(_lagr_hess); std::get<3>(_lagr_hess) = 0;
  }

public:
  /** @defgroup NLPSLV_IPOPT Local (Continuous) Optimization using IPOPT and MC++
   *  @{
   */
  //! @brief Constructor
  NLPSLV_IPOPT()
    : _nvar(0), _nctr(0), _obj_grad(0), _ctr_grad(0,0,0,0), _lagr_hess(0,0,0,0)
    {}

  //! @brief Destructor
  virtual ~NLPSLV_IPOPT()
    { _cleanup(); }

  //! @brief NLP solver options
  struct Options
  {
    //! @brief Constructor
    Options():
      CVTOL(1e-7), PRIMALTOL(1e-8), DUALTOL(1e-5), COMPLTOL(1e-7),
      MAXITER(100), MAXCPU(1e6), GRADIENT(FORWARD), HESSIAN(LBFGS),
      LINSLV(MA57), TESTDER(false), DISPLAY(0)
      {} 
    //! @brief Assignment operator
    Options& operator= ( Options&options ){
        CVTOL     = options.CVTOL;
        PRIMALTOL = options.PRIMALTOL;
        DUALTOL   = options.DUALTOL;
        COMPLTOL  = options.COMPLTOL;
        MAXITER   = options.MAXITER;
        MAXCPU    = options.MAXCPU;
        GRADIENT  = options.GRADIENT;
        HESSIAN   = options.HESSIAN;
        LINSLV    = options.LINSLV;
        TESTDER   = options.TESTDER;
        DISPLAY   = options.DISPLAY;
        return *this;
      }
    //! @brief Enumeration type for Hessian strategy
    enum HESSIAN_STRATEGY{
      EXACT=0, 	//!< Use exact second derivatives
      LBFGS,	//!< Perform a limited-memory quasi-Newton approximation
    };
    //! @brief Enumeration type for gradient strategy
    enum GRADIENT_STRATEGY{
      FORWARD=0,	//!< Forward AD
      BACKWARD		//!< Backward AD
    };
    //! @brief Enumeration type for gradient strategy
    enum LINEAR_SOLVER{
      MA27=0,       //!< use the Harwell routine MA27
      MA57,         //!< use the Harwell routine MA57
      MA77,         //!< use the Harwell routine HSL_MA77
      MA86,         //!< use the Harwell routine HSL_MA86
      MA97,         //!< use the Harwell routine HSL_MA97
      PARDISO,      //!< use the Pardiso package
      WSMP,         //!< use WSMP package
      MUMPS         //!< use MUMPS package
    };
    //! @brief Convergence tolerance
    double CVTOL;
    //! @brief Tolerance on primal feasibility
    double PRIMALTOL;
    //! @brief Tolerance on dual feasibility
    double DUALTOL;
     //! @brief Tolerance on complementarity conditions
    double COMPLTOL;
   //! @brief Maximum number of iterations
    int MAXITER;
    //! @brief Maximum run time (seconds)
    double MAXCPU;
    //! @brief Strategy for gradient computation
    GRADIENT_STRATEGY GRADIENT;
    //! @brief Strategy for Hessian computation during linesearch
    HESSIAN_STRATEGY HESSIAN;
    //! @brief Linear solver used for step computations
    LINEAR_SOLVER LINSLV;
    //! @brief Whether to test first- and second-derivatives
    bool TESTDER;
    //! @brief Display level in IPOPT
    int DISPLAY;
  } options;

  //! @brief NLP solver exceptions
  class Exceptions
  {
  public:
    //! @brief Enumeration type for NLPSLV_IPOPT exception handling
    enum TYPE{
      NOOBJ=1,		//!< Undefined objective function in optimization problem
      INTERN=-33	//!< Internal error
    };
    //! @brief Constructor for error <a>ierr</a>
    Exceptions( TYPE ierr ) : _ierr( ierr ){}
    //! @brief Inline function returning the error flag
    int ierr(){ return _ierr; }
    //! @brief Inline function returning the error description
    std::string what(){
      switch( _ierr ){
      case NOOBJ:
        return "NLPSLV_IPOPT::Exceptions  Undefined objective function in optimization problem";
      case INTERN: default:
        return "NLPSLV_IPOPT::Exceptions  Internal error";
      }
    }
  private:
    TYPE _ierr;
  };

  //! @brief Setup DAG for function/gradient/Hessian evaluation
  bool setup();

  //! @brief Solve NLP model -- return value is IPOPT status
  template <typename T>
  int solve
    ( const T*P, const double*p0=0 );

#ifdef MC__USE_SOBOL
  //! @brief Solve NLP model using multistart search -- return value is IPOPT status
  template <typename T>
  int solve
    ( const unsigned NSAM, const T*P, const bool*logscal=0, const bool DISP=true );
#endif

  //! @brief Test primal feasibility
  bool is_feasible
    ( const double*p, const double CTRTOL );

  //! @brief Test primal feasibility of current solution point
  bool is_feasible
    ( const double CTRTOL );

  //! @brief Test dual feasibility
  bool is_stationary
    ( const double*p, const double*ug, const double*upL, const double*upU,
      const double GRADTOL );

  //! @brief Test dual feasibility of current solution point
  bool is_stationary
    ( const double GRADTOL );//, const double NUMTOL=1e-8 );

  //! @brief Get IPOPT internal scaling value
  double get_scaling()
    {
      if( std::get<0>(_obj).size() != 1 ) throw Exceptions( Exceptions::NOOBJ );
      // set scaling factor. 1: minimize; -1: maximize
      switch( std::get<0>(_obj)[0] ){
        case MIN: _scaling =  1.; break;
        case MAX: _scaling = -1.; break;
      }
      return _scaling;
    }

  //! @brief Get IPOPT solution info
  const SOLUTION_OPT& solution() const
    {
      return _solution;
    }
  /** @} */

protected:
  //! @brief Method to return some info about the NLP
  virtual bool get_nlp_info
    ( Ipopt::Index& n, Ipopt::Index& m, Ipopt::Index& nnz_jac_g,
      Ipopt::Index& nnz_h_lag, IndexStyleEnum& index_style );

  //! @brief Method to return the variable and constraint bounds
  virtual bool get_bounds_info
    ( Ipopt::Index n, Ipopt::Number* x_l, Ipopt::Number* x_u,
      Ipopt::Index m, Ipopt::Number* g_l, Ipopt::Number* g_u );

  //! @brief Method to return the initial point for the NLP solver
  virtual bool get_starting_point
    ( Ipopt::Index n, bool init_x, Ipopt::Number* x, bool init_z,
      Ipopt::Number* z_L, Ipopt::Number* z_U, Ipopt::Index m,
      bool init_lambda, Ipopt::Number* lambda );

  //! @brief Method to return the objective function value
  virtual bool eval_f
    ( Ipopt::Index n, const Ipopt::Number* x, bool new_x,
      Ipopt::Number& obj_value );

  //! @brief Method to return the objective function gradient
  virtual bool eval_grad_f
    ( Ipopt::Index n, const Ipopt::Number* x, bool new_x,
      Ipopt::Number* grad_f );

  //! @brief Method to return the constraint residuals
  virtual bool eval_g
    ( Ipopt::Index n, const Ipopt::Number* x, bool new_x,
      Ipopt::Index m, Ipopt::Number* g );

  //! @brief Method to return the structure of the jacobian (if "values" is NULL) and the values of the jacobian (otherwise)
  virtual bool eval_jac_g
    ( Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Index m,
      Ipopt::Index nele_jac, Ipopt::Index* iRow, Ipopt::Index *jCol,
      Ipopt::Number* values );

  //! @brief Method to return the structure of the hessian of the lagrangian (if "values" is NULL) and the values of the hessian of the lagrangian (otherwise)
  virtual bool eval_h
    ( Ipopt::Index n, const Ipopt::Number* x, bool new_x,
      Ipopt::Number obj_factor, Ipopt::Index m, const Ipopt::Number* lambda,
      bool new_lambda, Ipopt::Index nele_hess, Ipopt::Index* iRow,
      Ipopt::Index* jCol, Ipopt::Number* values );

  //! @brief Method called when the algorithm is complete
  virtual void finalize_solution
    ( Ipopt::SolverReturn status, Ipopt::Index n, const Ipopt::Number* x,
      const Ipopt::Number* z_L, const Ipopt::Number* z_U, Ipopt::Index m,
      const Ipopt::Number* g, const Ipopt::Number* lambda,
      Ipopt::Number obj_value, const Ipopt::IpoptData* ip_data,
      Ipopt::IpoptCalculatedQuantities* ip_cq );

  //! @brief Method called at the end of each iteration
  virtual bool intermediate_callback
    ( Ipopt::AlgorithmMode mode, Ipopt::Index iter, Ipopt::Number obj_value,
      Ipopt::Number inf_pr, Ipopt::Number inf_du, Ipopt::Number mu,
      Ipopt::Number d_norm, Ipopt::Number regularization_size,
      Ipopt::Number alpha_du, Ipopt::Number alpha_pr,
      Ipopt::Index ls_trials, const Ipopt::IpoptData* ip_data,
      Ipopt::IpoptCalculatedQuantities* ip_cq );

private:
  //! @brief Set IPOPT options
  void set_options
    ( Ipopt::SmartPtr<Ipopt::IpoptApplication>&IpoptApp );

#ifdef MC__USE_SOBOL
  //! @brief Get next Sobol sampling point
  void _get_sobol
    ( const unsigned nr, double*r, long long int*pseed, const bool disp=false );
#endif

  //! @brief Private methods to block default compiler methods
  NLPSLV_IPOPT(const NLPSLV_IPOPT&);
  NLPSLV_IPOPT& operator=(const NLPSLV_IPOPT&);
};

inline
void
NLPSLV_IPOPT::set_options
( Ipopt::SmartPtr<Ipopt::IpoptApplication>&IpoptApp )
{
  IpoptApp->Options()->SetNumericValue( "tol", options.CVTOL<0.? 1e-12: options.CVTOL );
  IpoptApp->Options()->SetNumericValue( "dual_inf_tol", options.DUALTOL<=0.? 1e-12: options.DUALTOL );
  IpoptApp->Options()->SetNumericValue( "constr_viol_tol", options.PRIMALTOL<=0.? 1e-12: options.PRIMALTOL );
  IpoptApp->Options()->SetNumericValue( "compl_inf_tol", options.COMPLTOL<=0.? 1e-12: options.COMPLTOL );
  IpoptApp->Options()->SetIntegerValue( "max_iter", options.MAXITER );
  IpoptApp->Options()->SetNumericValue( "max_cpu_time", options.MAXCPU);
  IpoptApp->Options()->SetNumericValue( "obj_scaling_factor", get_scaling() );
  switch( options.HESSIAN ){
   case Options::EXACT:
    IpoptApp->Options()->SetStringValue( "hessian_approximation", "exact" );
    break;
   case Options::LBFGS:
    IpoptApp->Options()->SetStringValue( "hessian_approximation", "limited-memory" );
    break;
  }
  switch( options.LINSLV ){
   case Options::MA27:
    IpoptApp->Options()->SetStringValue( "linear_solver", "ma27" );
    break;
   case Options::MA57:
    IpoptApp->Options()->SetStringValue( "linear_solver", "ma57" );
    break;
   case Options::MA77:
    IpoptApp->Options()->SetStringValue( "linear_solver", "ma77" );
    break;
   case Options::MA86:
    IpoptApp->Options()->SetStringValue( "linear_solver", "ma86" );
    break;
   case Options::MA97:
    IpoptApp->Options()->SetStringValue( "linear_solver", "ma97" );
    break;
   case Options::PARDISO:
    IpoptApp->Options()->SetStringValue( "linear_solver", "pardiso" );
    break;
   case Options::WSMP:
    IpoptApp->Options()->SetStringValue( "linear_solver", "wsmp" );
    break;
   case Options::MUMPS:
    IpoptApp->Options()->SetStringValue( "linear_solver", "mumps" );
    break;
  }
  IpoptApp->Options()->SetStringValue( "derivative_test", options.TESTDER? "second-order": "none" );
  IpoptApp->Options()->SetIntegerValue( "print_level", options.DISPLAY<0? 0: (options.DISPLAY>12? 12: options.DISPLAY ) );
}

inline
bool
NLPSLV_IPOPT::setup
()
{
  // full set of decision variables (independent & dependent)
  _var = BASE_NLP::_var;
  _var.insert( _var.end(), _dep.begin(), _dep.end() );
  _nvar = _var.size();

  // full set of constraints (main & equation system)
  //_ctr = BASE_NLP::_ctr;
  std::get<0>(_ctr).clear();
  std::get<1>(_ctr).clear();
  std::get<2>(_ctr).clear();
  for( unsigned ic=0; ic<std::get<0>(BASE_NLP::_ctr).size(); ic++ ){
    if( std::get<3>(BASE_NLP::_ctr)[ic] ) continue; // redundant constraint - ignore
    std::get<0>(_ctr).push_back( std::get<0>(BASE_NLP::_ctr)[ic] );
    std::get<1>(_ctr).push_back( std::get<1>(BASE_NLP::_ctr)[ic] );
    std::get<2>(_ctr).push_back( std::get<2>(BASE_NLP::_ctr)[ic] );
  }
  for( auto its=_sys.begin(); its!=_sys.end(); ++its ){
    std::get<0>(_ctr).push_back( EQ );
    std::get<1>(_ctr).push_back( (*its) );
    std::get<2>(_ctr).push_back( FFVar( _dag ) );
  }
  _nctr = std::get<0>(_ctr).size();

  // setup objective and constraint evaluation
  _cleanup();
  if( std::get<0>(_obj).size() != 1 ) throw Exceptions( Exceptions::NOOBJ );
  _op_f = _dag->subgraph( std::get<1>(_obj).size(), std::get<1>(_obj).data() );
  _op_g = _dag->subgraph( std::get<1>(_ctr).size(), std::get<1>(_ctr).data() );

  // setup objective and constraint gradient evaluation
  switch( options.GRADIENT ){
    case Options::FORWARD:  // Forward AD
      _obj_grad = _dag->FAD( std::get<1>(_obj).size(), std::get<1>(_obj).data(),
                             _nvar, _var.data() ); 
      _ctr_grad = _dag->SFAD( std::get<1>(_ctr).size(), std::get<1>(_ctr).data(),
                              _nvar, _var.data() ); 
      break;
    case Options::BACKWARD: // Backward AD
      _obj_grad = _dag->BAD( std::get<1>(_obj).size(), std::get<1>(_obj).data(),
                             _nvar, _var.data() ); 
      _ctr_grad = _dag->SBAD( std::get<1>(_ctr).size(), std::get<1>(_ctr).data(),
                              _nvar, _var.data() ); 
      break;
  }
#ifdef MC__NLPSLV_IPOPT_DEBUG
  for( unsigned ie=0; ie<std::get<0>(_ctr_grad); ++ie )
     std::cout << "  dg[" << std::get<1>(_ctr_grad)[ie]
               << "," << std::get<2>(_ctr_grad)[ie] << "]";
  std::cout << std::endl;
#endif
  _op_df = _dag->subgraph( _nvar, _obj_grad );
  _op_dg = _dag->subgraph( std::get<0>(_ctr_grad), std::get<3>(_ctr_grad) );

  // setup Lagrangian Hessian evaluation
  //if( options.HESSIAN == Options::EXACT ){
    // form Lagrangian
    FFVar lagr = std::get<1>(_obj)[0] * std::get<2>(_obj)[0];
    for( unsigned ic=0; ic<_nctr; ic++ )
      lagr += std::get<1>(_ctr)[ic] * std::get<2>(_ctr)[ic];
    // differentiate twice
    const FFVar* lagr_grad = _dag->BAD( 1, &lagr, _nvar, _var.data() );
    _lagr_hess = _dag->SFAD( _nvar, lagr_grad, _nvar, _var.data(), true );
    delete[] lagr_grad;
#ifdef MC__NLPSLV_IPOPT_DEBUG
    _dag->output( _dag->subgraph( std::get<0>(_lagr_hess), std::get<3>(_lagr_hess) ), std::cout );
    for( unsigned ie=0; ie<std::get<0>(_lagr_hess); ++ie )
       std::cout << "  d2L[" << std::get<1>(_lagr_hess)[ie]
                 << "," << std::get<2>(_lagr_hess)[ie] << "]: "
                 << std::get<3>(_lagr_hess)[ie];
    std::cout << std::endl;
#endif
    _op_d2L = _dag->subgraph( std::get<0>(_lagr_hess), std::get<3>(_lagr_hess) );
  //}

  return true;
}

template <typename T>
inline
int
NLPSLV_IPOPT::solve
( const T*P, const double*p0 )
{
  Ipopt::SmartPtr<Ipopt::IpoptApplication> IpoptApp = new Ipopt::IpoptApplication();

  // Set (a few) IPOPT options
  set_options( IpoptApp );

  // Keep track of bounds and initial guess
  _data.set( _nvar, p0, P );

  // Run NLP solver
  _solution.status = IpoptApp->Initialize();
  if( _solution.status == Ipopt::Solve_Succeeded )
    _solution.status = IpoptApp->OptimizeTNLP( this );

  return _solution.status;
}

#ifdef MC__USE_SOBOL
inline
void
NLPSLV_IPOPT::_get_sobol
( const unsigned nr, double*r, long long int*pseed, const bool disp )
{
  //if( disp ) std::cout << std::setw(6) << *pseed << "  ";
  i8_sobol( nr, pseed, r );
  if( disp ){
    //std::cout << std::setw(6) << *pseed << "  ";
    for( unsigned j=0; j<nr; j++ ){
      std::cout << std::setw(6) << r[j] << "  ";
    }
    std::cout << std::endl;
  }
}

template <typename T>
inline
int
NLPSLV_IPOPT::solve
( const unsigned NSAM, const T*P, const bool*logscal, const bool DISP )
{
  long long SEED = 1;
  std::vector<double> vSAM(_nvar), p0(_nvar);

  bool found = false;
  SOLUTION_OPT best;

  Ipopt::SmartPtr<Ipopt::IpoptApplication> IpoptApp = new Ipopt::IpoptApplication();
  set_options( IpoptApp );
  if( DISP ) std::cout << "Multistart: ";

  for(unsigned k=0; k<NSAM; k++){
    // Sample variable domain
    _get_sobol( _nvar, vSAM.data(), &SEED );
    for( unsigned i=0; i<_nvar; i++ ){
      if( !logscal || !logscal[i] || Op<T>::l(P[i]) <= 0. )
        p0[i] = Op<T>::l(P[i]) + Op<T>::diam(P[i]) * vSAM[i];
      else
        p0[i] = std::exp( Op<T>::l(Op<T>::log(P[i])) + Op<T>::diam(Op<T>::log(P[i])) * vSAM[i] );
      //std::cout << "  " << p0[i];
    }
    //std::cout << std::endl;

    // Solve NLP
    _data.set( _nvar, p0.data(), P );
    _solution.status = IpoptApp->Initialize();
    if( _solution.status != Ipopt::Solve_Succeeded ) continue;
    _solution.status = IpoptApp->OptimizeTNLP( this );

    // Test feasibility and improvement
    if( !is_feasible( 1.1*options.PRIMALTOL ) ){
      if( DISP ) std::cout << "·";
      continue;
    }
    if( DISP ) std::cout << "*";
    if( !found ){
      best = _solution;
      found = true;
      continue;
    }
    switch( std::get<0>(_obj)[0] ){
      case MIN:
        if( _solution.f < best.f ) best = _solution;
        break;
      case MAX:
        if( _solution.f > best.f ) best = _solution;
        break;
    }
  }
  if( DISP ) std::cout << std::endl;

  _solution = best;
  return _solution.status;
}
#endif

inline
bool
NLPSLV_IPOPT::is_feasible
( const double*p, const double CTRTOL )
{
  _solution.p.assign( p, p+_nvar );
  return is_feasible( CTRTOL );
}

inline
bool
NLPSLV_IPOPT::is_feasible
( const double CTRTOL )
{
  if( !eval_f( _solution.p.size(), _solution.p.data(), true, _solution.f )
   || !eval_g( _solution.p.size(), _solution.p.data(), false, _solution.g.size(), _solution.g.data() ) )
  return false;

  double maxinfeas=0.;
  auto itc = std::get<0>(_ctr).begin();
  for( unsigned ic=0; itc!=std::get<0>(_ctr).end(); ++itc, ic++ ){
    switch( (*itc) ){
      case EQ: maxinfeas = std::fabs(_solution.g[ic]); break;
      case LE: maxinfeas = _solution.g[ic];            break;
      case GE: maxinfeas = -_solution.g[ic];           break;
    }
#ifdef MC__DOSEQSLV_IPOPT_DEBUG
    std::cout << "g[" << ic << "] = " << _solution.g[ic] << "  (" << maxinfeas << ")\n";
#endif
    if( maxinfeas > CTRTOL ) return false;
  }

  return true;
}

inline
bool
NLPSLV_IPOPT::is_stationary
( const double*p, const double*ug, const double*upL, const double*upU,
  const double GRADTOL )
{
  _solution.p.assign( p, p+_nvar );
  _solution.ug.assign( ug, ug+_nctr );
  _solution.upL.assign( upL, upL+_nvar );
  _solution.upU.assign( upU, upU+_nvar );
  return is_stationary( GRADTOL );
}

inline
bool
NLPSLV_IPOPT::is_stationary
( const double GRADTOL )
{
  CPPL::dcovector grad_f( _solution.p.size() ), jac_g( std::get<0>(_ctr_grad) );
  if( !eval_grad_f( _solution.p.size(), _solution.p.data(), true, grad_f.array )
   || !eval_jac_g( _solution.p.size(), _solution.p.data(), false, _solution.g.size(),
                   std::get<0>(_ctr_grad), 0, 0, jac_g.array ) )
  return false;

#ifdef MC__NLPSLV_IPOPT_DEBUG
  for( unsigned ip=0; ip<_solution.p.size(); ip++ )
    std::cout << ip << ": " << grad_f(ip) << std::endl;
#endif
  CPPL::dcovector mul_g( _solution.g.size() );
  CPPL::dgsmatrix sjac_g( _solution.p.size(), _solution.g.size(), std::get<0>(_ctr_grad) );
  for( unsigned ie=0; ie<std::get<0>(_ctr_grad); ie++ )
    sjac_g.put( std::get<2>(_ctr_grad)[ie], std::get<1>(_ctr_grad)[ie], jac_g(ie) );
#ifdef MC__NLPSLV_IPOPT_DEBUG
  std::cout << std::endl << sjac_g.to_dgematrix();
#endif

  CPPL::dcovector grad_L = grad_f; //_scaling * grad_f;
  for( unsigned ig=0; ig<_solution.g.size(); ig++ ){
#ifdef MC__NLPSLV_IPOPT_DEBUG
    std::cout << ig << ": " << _solution.ug[ ig ] << std::endl;
#endif
    mul_g(ig) = _solution.ug[ig];
  }
  grad_L += sjac_g * mul_g;

  for( unsigned ip=0; ip<_solution.p.size(); ip++ ){
#ifdef MC__NLPSLV_IPOPT_DEBUG
    std::cout << ip << ": " << _solution.upL[ ip ] << "  " << _solution.upU[ ip ] << "  "
              << grad_L(ip) << std::endl;
#endif
    grad_L(ip) += _solution.upU[ ip ] - _solution.upL[ ip ];
    if( std::fabs( grad_L(ip) ) > GRADTOL ) return false;
  }
  return true;
}

//inline
//bool
//NLPSLV_IPOPT::is_stationary
//( const double GRADTOL, const double NUMTOL )
//{
//  CPPL::dcovector grad_f( _solution.p.size() ), jac_g( std::get<0>(_ctr_grad) );
//  if( !eval_grad_f( _solution.p.size(), _solution.p.data(), true, grad_f.array )
//   || !eval_jac_g( _solution.p.size(), _solution.p.data(), false, _solution.g.size(),
//                   std::get<0>(_ctr_grad), 0, 0, jac_g.array ) )
//  return false;

//#ifdef MC__NLPSLV_IPOPT_DEBUG
//  for( unsigned ip=0; ip<_solution.p.size(); ip++ )
//    std::cout << ip << ": " << grad_f(ip) << std::endl;
//#endif
//  CPPL::dgsmatrix sjac_g( _solution.g.size()+_solution.p.size(), _solution.p.size(), std::get<0>(_ctr_grad) );
//  for( unsigned ie=0; ie<std::get<0>(_ctr_grad); ie++ ){
//#ifdef MC__NLPSLV_IPOPT_DEBUG
//    std::cout << std::get<1>(_ctr_grad)[ie] << ": " << _solution.ug[ std::get<1>(_ctr_grad)[ie] ] << std::endl;
//#endif
//    if( std::fabs( _solution.ug[ std::get<1>(_ctr_grad)[ie] ] ) < NUMTOL ) continue;
//    sjac_g.put( std::get<1>(_ctr_grad)[ie], std::get<2>(_ctr_grad)[ie], jac_g(ie) );
//  }
//  for( unsigned ip=0; ip<_solution.p.size(); ip++ ){
//#ifdef MC__NLPSLV_IPOPT_DEBUG
//    std::cout << ip << ": " << _solution.upL[ ip ] << "  " << _solution.upU[ ip ] << std::endl;
//#endif
//    if( std::fabs( _solution.upL[ ip ] ) < NUMTOL
//     && std::fabs( _solution.upU[ ip ] ) < NUMTOL ) continue;
//    sjac_g.put( _solution.g.size()+ip, ip, 1 );
//  }
//#ifdef MC__NLPSLV_IPOPT_DEBUG
//  std::cout << std::endl << sjac_g.to_dgematrix();
//#endif

//  CPPL::dcovector S;
//  CPPL::dgematrix U, VT, G = sjac_g.to_dgematrix();
//  G.dgesvd( S, U, VT );
//  unsigned nulldim = 0;
//  for( ; nulldim<_solution.p.size(); nulldim++ )
//    if( std::fabs( S(_solution.p.size()-nulldim-1) ) > NUMTOL ) break;
//  // Point is fully detemrined by constraints
//  if( !nulldim ) return true;
//  // Project cost gradient on constraint null space
//  CPPL::dgematrix proj( nulldim, _solution.p.size() );
//  for( unsigned i=0; i<nulldim; i++ )
//    for( unsigned j=0; j<_solution.p.size(); j++ )
//      proj(j,i) = VT(_solution.p.size()-i-1,j);
//#ifdef MC__NLPSLV_IPOPT_DEBUG
//  std::cout << proj * grad_f;
//#endif
//  return( CPPL::nrm2( proj * grad_f ) > GRADTOL ? false: true );
//}

inline
bool
NLPSLV_IPOPT::get_nlp_info
( Ipopt::Index& n, Ipopt::Index& m, Ipopt::Index& nnz_jac_g,
  Ipopt::Index& nnz_h_lag, IndexStyleEnum& index_style )
{
#ifdef MC__NLPSLV_IPOPT_TRACE
    std::cout << "  NLPSLV_IPOPT::get_nlp_info\n";
#endif
  // set size
  n = _nvar; m = _nctr;
  nnz_jac_g = std::get<0>(_ctr_grad);
  nnz_h_lag = std::get<0>(_lagr_hess);

  // use the C style indexing (0-based)
  index_style = Ipopt::TNLP::C_STYLE;

#ifdef MC__NLPSLV_IPOPT_DEBUG
  std::cout << "n:" << n << std::endl;
  std::cout << "m:" << m << std::endl;
  std::cout << "nnz_jac_g:" << nnz_jac_g << std::endl;
  std::cout << "nnz_h_lag:" << nnz_h_lag << std::endl;
#endif

  return true;
}

inline
bool
NLPSLV_IPOPT::get_bounds_info
( Ipopt::Index n, Ipopt::Number* x_l, Ipopt::Number* x_u,
  Ipopt::Index m, Ipopt::Number* g_l, Ipopt::Number* g_u )
{
#ifdef MC__NLPSLV_IPOPT_TRACE
    std::cout << "  NLPSLV_IPOPT::get_bounds_info\n";
#endif
  // set variable bounds
  unsigned ip = 0;
  for( ; ip<_nvar; ip++ ){   
    x_l[ip] = ( _data.P? _data.P[ip].first: -BASE_OPT::INF );
    x_u[ip] = ( _data.P? _data.P[ip].second: BASE_OPT::INF );
#ifdef MC__NLPSLV_IPOPT_DEBUG
    std::cout << "  x_l[" << ip << "] = " << x_l[ip]
              << "  x_u[" << ip << "] = " << x_u[ip] << std::endl;
#endif
  }

  // set constraint bounds
  auto itc = std::get<0>(_ctr).begin();
  for( unsigned ic=0; itc!=std::get<0>(_ctr).end(); ++itc, ic++ ){
    switch( (*itc) ){
      case EQ: g_l[ic] = g_u[ic] = 0.; break;
      case LE: g_l[ic] = -BASE_OPT::INF; g_u[ic] = 0.; break;
      case GE: g_l[ic] = 0.; g_u[ic] = BASE_OPT::INF; break;
    }
#ifdef MC__NLPSLV_IPOPT_DEBUG
    std::cout << "  g_l[" << ic << "] = " << g_l[ic]
              << "  g_u[" << ic << "] = " << g_u[ic] << std::endl;
#endif
  }
  return true;
}

inline
bool
NLPSLV_IPOPT::get_starting_point
( Ipopt::Index n, bool init_x, Ipopt::Number* x, bool init_z,
  Ipopt::Number* z_L, Ipopt::Number* z_U, Ipopt::Index m,
  bool init_lambda, Ipopt::Number* lambda )
{
#ifdef MC__NLPSLV_IPOPT_TRACE
    std::cout << "  NLPSLV_IPOPT::get_starting_point  "
              << init_x << init_z << init_lambda << std::endl;
#endif
  // Here, we assume we only have starting values for x, if you code
  // your own NLP, you can provide starting values for the dual variables
  // if you wish
  if( !init_x || init_z || init_lambda ) return false;

  // initialize to the given starting point
  //if( !_data.p0 ) return false;
  for( unsigned ip=0; ip<_nvar; ip++ ){   
    x[ip] = _data.p0? _data.p0[ip]: 0.;
#ifdef MC__NLPSLV_IPOPT_DEBUG
    std::cout << "  x_0[" << ip << "] = " << x[ip] << std::endl;
#endif
  }
  return true;
}

inline
bool
NLPSLV_IPOPT::eval_f
( Ipopt::Index n, const Ipopt::Number* x, bool new_x,
  Ipopt::Number& f )
{
  assert( (unsigned)n == _nvar );
#ifdef MC__NLPSLV_IPOPT_TRACE
  std::cout << "  NLPSLV_IPOPT::eval_f  " << new_x << std::endl;
  for( Ipopt::Index ip=0; ip<n; ip++ )
    std::cout << "  x[" << ip << "] = " << x[ip] << std::endl;
#endif
  // evaluate objective
  try{
    _dag->eval( _op_f, _dwk, 1, std::get<1>(_obj).data(), &f, n, _var.data(), x );
  }
  catch(...){
    return false;
  }
#ifdef MC__NLPSLV_IPOPT_DEBUG
  std::cout << "  f = " << f << std::endl;
#endif
  return true;
}

inline
bool
NLPSLV_IPOPT::eval_grad_f
( Ipopt::Index n, const Ipopt::Number* x, bool new_x,
  Ipopt::Number* df )
{
  assert( (unsigned)n == _nvar );
#ifdef MC__NLPSLV_IPOPT_TRACE
  std::cout << "  NLPSLV_IPOPT::eval_grad_f  " << new_x << std::endl;
  for( Ipopt::Index ip=0; ip<n; ip++ )
    std::cout << "  x[" << ip << "] = " << x[ip] << std::endl;
#endif
  // evaluate objective gradient
  try{
    _dag->eval( _op_df, _dwk, n, _obj_grad, df, n, _var.data(), x );
#ifdef MC__NLPSLV_IPOPT_DEBUG
    for( Ipopt::Index i=0; i<n; i++ )
      std::cout << "  df[" << i << "] = " << df[i] << std::endl;
#endif
  }
  catch(...){
    return false;
  }
  return true;
}

inline
bool
NLPSLV_IPOPT::eval_g
( Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Index m,
  Ipopt::Number* g )
{
  assert( (unsigned)n == _nvar && (unsigned)m == _nctr );
#ifdef MC__NLPSLV_IPOPT_TRACE
  std::cout << "  NLPSLV_IPOPT::eval_g  " << new_x << std::endl;
#endif
  // evaluate constraints
  try{
    _dag->eval( _op_g, _dwk, m, std::get<1>(_ctr).data(), g, n, _var.data(), x );
#ifdef MC__NLPSLV_IPOPT_DEBUG
    for( Ipopt::Index ic=0; ic<m; ic++ )
      std::cout << "  g[" << ic << "] = " << g[ic] << std::endl;
#endif
  }
  catch(...){
    return false;
  }
  return true;
}

inline
bool
NLPSLV_IPOPT::eval_jac_g
( Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Index m,
  Ipopt::Index nele_jac, Ipopt::Index* iRow, Ipopt::Index *jCol,
  Ipopt::Number* dg )
{
  assert( (unsigned)n == _nvar && (unsigned)m == _nctr );
#ifdef MC__NLPSLV_IPOPT_TRACE
  std::cout << "  NLPSLV_IPOPT::eval_jac_g  " << new_x << std::endl;
#endif
  if( (unsigned)nele_jac != std::get<0>(_ctr_grad) ) return false;

  // return the constraint Jacobian structure
  if( !dg ){
    Ipopt::Index ie = 0;
    for( ; ie<nele_jac; ++ie ){
      iRow[ie] = std::get<1>(_ctr_grad)[ie];
      jCol[ie] = std::get<2>(_ctr_grad)[ie];
#ifdef MC__NLPSLV_IPOPT_DEBUG
      std::cout << "  dg[" << iRow[ie] << ", " << jCol[ie] << "]" << std::endl;
#endif
    }
    return true;
  }

  // evaluate constraint gradient
  try{
    _dag->eval( _op_dg, _dwk, nele_jac, std::get<3>(_ctr_grad), dg, n, _var.data(), x );
#ifdef MC__NLPSLV_IPOPT_DEBUG
    for( Ipopt::Index ie=0; ie<nele_jac; ++ie ){
       std::cout << "  dg[" << std::get<1>(_ctr_grad)[ie]
                 << ", " << std::get<2>(_ctr_grad)[ie] << "] = " << dg[ie]
                 << std::endl;
    }
#endif
  }
  catch(...){
    return false;
  }

#ifdef MC__NLPSLV_IPOPT_DEBUG
  int dum; std::cin >> dum;
#endif
  return true;
}

inline
bool
NLPSLV_IPOPT::eval_h
( Ipopt::Index n, const Ipopt::Number* x, bool new_x,
  Ipopt::Number obj_factor, Ipopt::Index m, const Ipopt::Number* lambda,
  bool new_lambda, Ipopt::Index nele_hess, Ipopt::Index* iRow,
  Ipopt::Index* jCol, Ipopt::Number* d2L )
{
  assert( (unsigned)n == _nvar && (unsigned)m == _nctr );
#ifdef MC__NLPSLV_IPOPT_TRACE
  std::cout << "  NLPSLV_IPOPT::eval_h  " << new_x  << new_lambda << std::endl;
#endif
  if( (unsigned)nele_hess != std::get<0>(_lagr_hess) ) return false;

  // return the Lagrangian Hessian structure
  if( !d2L ){
    Ipopt::Index ie = 0;
    for( ; ie<nele_hess; ++ie ){
      iRow[ie] = std::get<1>(_lagr_hess)[ie];
      jCol[ie] = std::get<2>(_lagr_hess)[ie];
#ifdef MC__NLPSLV_IPOPT_DEBUG
      std::cout << "  d2L[" << iRow[ie] << ", " << jCol[ie] << "]" << std::endl;
#endif
    }
    return true;
  }

  // evaluate Lagrangian Hessian
  try{
    std::list<unsigned> nvar; std::list<const FFVar*> pvar; std::list<const double*> vvar;
    nvar.push_back(n); pvar.push_back(_var.data());              vvar.push_back(x);
    nvar.push_back(m); pvar.push_back(std::get<2>(_ctr).data()); vvar.push_back(lambda);
    nvar.push_back(1); pvar.push_back(std::get<2>(_obj).data()); vvar.push_back(&obj_factor);
    _dag->eval( _op_d2L, _dwk, nele_hess, std::get<3>(_lagr_hess), d2L, nvar, pvar, vvar );
#ifdef MC__NLPSLV_IPOPT_DEBUG
    for( Ipopt::Index ie=0; ie<nele_hess; ++ie ){
       std::cout << "  d2L[" << std::get<1>(_lagr_hess)[ie]
                 << ", " << std::get<2>(_lagr_hess)[ie] << "] = " << d2L[ie]
                 << std::endl;
    }
#endif
  }
  catch(...){
    return false;
  }
  return true;
}

inline
void
NLPSLV_IPOPT::finalize_solution
( Ipopt::SolverReturn status, Ipopt::Index n, const Ipopt::Number* p,
  const Ipopt::Number* upL, const Ipopt::Number* upU, Ipopt::Index m,
  const Ipopt::Number* g, const Ipopt::Number* ug, Ipopt::Number f,
  const Ipopt::IpoptData* ip_data, Ipopt::IpoptCalculatedQuantities* ip_cq )
{
#ifdef MC__NLPSLV_IPOPT_TRACE
    std::cout << "  DOSEQSLV_IPOPT::finalize_solution\n";
#endif
  _solution.status = status;

  // Successful (or near-successful) completion
  //if( status == Ipopt::SUCCESS || status == Ipopt::STOP_AT_ACCEPTABLE_POINT
  // || status == Ipopt::STOP_AT_TINY_STEP ){
    _solution.f = f;
    _solution.p.assign( p, p+n );
    _solution.upL.assign( upL, upL+n );
    _solution.upU.assign( upU, upU+n );
    _solution.g.assign( g, g+m );
    _solution.ug.assign( ug, ug+m );
  //}

  // Failure
  //else{
  //  _solution.p.clear();
  //  _solution.upL.clear();
  //  _solution.upU.clear();
  //  _solution.g.clear();
  //  _solution.ug.clear();
  //}
  return;
}

inline
bool
NLPSLV_IPOPT::intermediate_callback
( Ipopt::AlgorithmMode mode, Ipopt::Index iter, Ipopt::Number obj_value,
  Ipopt::Number inf_pr, Ipopt::Number inf_du, Ipopt::Number mu,
  Ipopt::Number d_norm, Ipopt::Number regularization_size,
  Ipopt::Number alpha_du, Ipopt::Number alpha_pr, Ipopt::Index ls_trials,
  const Ipopt::IpoptData* ip_data, Ipopt::IpoptCalculatedQuantities* ip_cq )
{
#ifdef MC__NLPSLV_IPOPT_TRACE
    std::cout << "  NLPSLV_IPOPT::intermediate_callback\n";
#endif
  return true;
}

} // end namescape mc

#endif
