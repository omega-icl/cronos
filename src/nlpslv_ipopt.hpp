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
    std::cout << "  f* = " << pNLP->get_objective() << std::endl;
    for( unsigned ip=0; ip<NP; ip++ )
      std::cout << "  p*(" << ip << ") = " << pNLP->get_variable(ip) << std::endl;
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

#include "base_nlp.hpp"

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
  std::list<const FFOp*> _op_f;
  //! @brief list of operations for objective gradient evaluation
  std::list<const FFOp*> _op_df;
  //! @brief list of operations for constraint evaluation
  std::list<const FFOp*> _op_g;
  //! @brief list of operations for constraint gradient evaluation
  std::list<const FFOp*> _op_dg;
  //! @brief list of operations for Lagragian Hessian evaluation
  std::list<const FFOp*> _op_d2L;

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
  struct SOLUTION{
    SOLUTION(): n(0), p(0), upL(0), upU(0), m(0), g(0), ug(0) {}
    ~SOLUTION() { delete[] p; delete[] upL; delete[] upU; delete[] g; delete[] ug; }
    Ipopt::ApplicationReturnStatus status;
    Ipopt::SolverReturn solverflag;
    Ipopt::Index n;
    Ipopt::Number*p;
    Ipopt::Number*upL;
    Ipopt::Number*upU;
    Ipopt::Index m;
    Ipopt::Number*g;
    Ipopt::Number*ug;
    Ipopt::Number f;
  } _solution;

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
      CVTOL(1e-8), PRIMALTOL(1e-8), DUALTOL(1e-4), COMPLTOL(1e-7),
      MAXITER(100), MAXCPU(1e6), GRADIENT(FORWARD), HESSIAN(LBFGS),
      TESTDER(false), DISPLAY(0)
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
  Ipopt::ApplicationReturnStatus solve
    ( const T*P, const double*p0=0 );

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
  const SOLUTION& solution() const
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
  _ctr = BASE_NLP::_ctr;
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
Ipopt::ApplicationReturnStatus
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
    x_l[ip] = ( _data.P? _data.P[ip].first: -INF );
    x_u[ip] = ( _data.P? _data.P[ip].second: INF );
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
      case LE: g_l[ic] = -INF; g_u[ic] = 0.; break;
      case GE: g_l[ic] = 0.; g_u[ic] = INF; break;
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
    _dag->eval( _op_f, 1, std::get<1>(_obj).data(), &f, n, _var.data(), x );
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
    _dag->eval( _op_df, n, _obj_grad, df, n, _var.data(), x );
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
    _dag->eval( _op_g, m, std::get<1>(_ctr).data(), g, n, _var.data(), x );
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
    _dag->eval( _op_dg, nele_jac, std::get<3>(_ctr_grad), dg, n, _var.data(), x );
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
    _dag->eval( _op_d2L, nele_hess, std::get<3>(_lagr_hess), d2L, nvar, pvar, vvar );
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
    std::cout << "  NLPSLV_IPOPT::finalize_solution\n";
#endif
  _solution.solverflag = status;

  // Successful (or near-successful) completion
  //if( status == Ipopt::SUCCESS || status == Ipopt::STOP_AT_ACCEPTABLE_POINT
  // || status == Ipopt::STOP_AT_TINY_STEP ){
    // resize solution arrays
    if( _solution.n != n ){
      if( _solution.n ){
        delete[] _solution.p; delete[] _solution.upL; delete[] _solution.upU;
      }
      _solution.n = n;
      _solution.p = new Ipopt::Number[n];
      _solution.upL = new Ipopt::Number[n];
      _solution.upU = new Ipopt::Number[n];
    }
    if( _solution.m != m ){
      if( _solution.m ){
        delete[] _solution.g; delete[] _solution.ug;
      }
      _solution.m = m;
      _solution.g = new Ipopt::Number[m];
      _solution.ug = new Ipopt::Number[m];
    }
    // copy solution values into _solution
    _solution.f = f;
    for( Ipopt::Index i=0; i<n; i++ ){
      _solution.p[i] = p[i]; _solution.upL[i] = upL[i]; _solution.upU[i] = upU[i];
    }
    for( Ipopt::Index j=0; j<m; j++ ){
      _solution.g[j] = g[j]; _solution.ug[j] = ug[j];
    }
  //}

  // Failure
  //else{
  //  if( _solution.n ){
  //    delete[] _solution.p; delete[] _solution.upL; delete[] _solution.upU;
  //    _solution.n = 0;
  //  }
  //  if( _solution.m ){
  //    delete[] _solution.g; delete[] _solution.ug;
  //    _solution.m = 0;
  //  }
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
