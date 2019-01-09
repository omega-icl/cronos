// Copyright (C) 2012-2017 Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

/*!
\page page_DOSEQSLV Local Dynamic Optimization using IPOPT, CVODES and MC++ (Sequential Approach)
\author Benoit C. Chachuat
\version 0.1
\date 2016
\bug No known bugs.

Consider a dynamic optimization (DO) problem of the form
\f{align*}
  \min_{{\bf p}\in P}\ & \phi_0({\bf p}, {\bf x}(t_0),\ldots,{\bf x}(t_N)) + \int_{t_0}^{t_N} \psi_0({\bf p}, {\bf x}(t))\, dt\\
  \text{s.t.} \ & \phi_k({\bf p}, {\bf x}(t_0),\ldots,{\bf x}(t_N)) + \int_{t_0}^{t_N} \psi_k({\bf p}, {\bf x}(t))\, dt \{\geq,=,\leq\} 0, \quad k=1,\ldots,n_c\\
  & \dot{{\bf x}}(t)={\bf f}({\bf p},{\bf x}(t)),\ t\in(t_0,t_{N}]; \quad {\bf x}(t_0)={\bf h}({\bf p})
\f}
where \f${\bf p}\in P\subset\mathbb{R}^{n_p}\f$, \f${\bf x}\in\mathbb{R}^{n_x}\f$, and the real-valued functions \f${\bf f}\f$, \f${\bf h}\f$, \f$\phi_k\f$ and \f$\psi_k\f$ are twice continuously differentiable in all their arguments and factorable. The class mc::DOSEQSLV_IPOPT solves such DO problems to local optimality, based on the sequential approach. It relies on the software packages <A href="https://projects.coin-or.org/Ipopt">IPOPT</A> for local solution of NLP problems, and <A href="https://computation.llnl.gov/casc/sundials/description/description.html">SUNDIALS/CVODES</A> for numerical integration and sensitivity analysis of IVPs in ODEs. The dynamic system, objective function and constraints are defined, manipulated, and analyzed using direct acyclic graphs (DAG) in <A href="https://projects.coin-or.org/MCpp">MC++</A>.

\section sec_DOSEQSLV_solve How to Solve a Dynamic Optimization Model using mc::DOSEQSLV_IPOPT?

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

We start by defining an mc::DOSEQSLV_IPOPT class as below, whereby the class <A href="http://www.coin-or.org/Doxygen/Ipopt/class_ipopt_1_1_smart_ptr.html">Ipopt::SmartPtr</A> is used here to comply with the IPOPT C++ interface:

\code
  #include "doseqslv_ipopt.hpp"

  Ipopt::SmartPtr<mc::DOSEQSLV_IPOPT> OC = new mc::DOSEQSLV_IPOPT;
\endcode

Next, we set the variables and objective/constraint functions by creating a DAG of the problem: 

\code
  mc::FFGraph IVP;             // DAG of the IVP
  OC->set_dag( &IVP );

  const unsigned int NS = 10;  // Time stages
  double tk[NS+1]; tk[0] = 0.;
  for( unsigned k=0; k<NS; k++ ) tk[k+1] = tk[k] + 1./(double)NS;
  OC->set_time( NS, tk );

  const unsigned NP = 2+NS;    // Decision variables
  mc::FFVar P[NP];
  for( unsigned int i=0; i<NP; i++ ) P[i].set( &IVP );
  mc::FFVar TF = P[0], X10 = P[1], *U = P+2;
  OC->set_parameter( NP, P );

  const unsigned NX = 2;       // State variables
  mc::FFVar X[NX];
  for( unsigned int i=0; i<NX; i++ ) X[i].set( &IVP );
  OC->set_state( NX, X );

  mc::FFVar RHS[NX*NS];        // Right-hand side function
  for( unsigned k=0; k<NS; k++ ){
    RHS[NX*k+0] = TF * X[1];
    RHS[NX*k+1] = TF * ( U[k]*X[0] - 2.*X[1] );
  }
  OC->set_differential( NS, NX, RHS );
 
  mc::FFVar IC[NX];            // Initial value function
  IC[0] = 0.;
  IC[1] = X10;
  OC->set_initial( NX, IC );

  OC->set_obj( mc::DOSEQSLV_IPOPT::MIN, TF );         // objective
  OC->add_ctr( mc::DOSEQSLV_IPOPT::EQ, std::make_pair( NS-1, X[0]-1. ) );        // constraint #1
  OC->add_ctr( mc::DOSEQSLV_IPOPT::EQ, std::make_pair( NS-1, X[1] ) );        // constraint #2
\endcode

Options, such as the output level, maximum number of iterations, tolerance, maximum CPU time, etc can all be modified through the public member mc::DOSEQSLV_IPOPT::Options:

\code
  OC->options.DISPLAY   = 5;
  OC->options.MAXITER   = 100;
  OC->options.TESTDER   = false;
  OC->options.CVTOL     = 1e-7;
  OC->options.GRADIENT  = mc::DOSEQSLV_IPOPT::Options::BACKWARD;
  OC->options.INTMETH   = mc::DOSEQSLV_IPOPT::Options::MSBDF;
  OC->options.JACAPPROX = mc::DOSEQSLV_IPOPT::Options::CV_DENSE;
  OC->options.NMAX      = 20000;
  OC->options.ATOL      = OC->options.ATOLB      = OC->options.ATOLS  = 1e-9;
  OC->options.RTOL      = OC->options.RTOLB      = OC->options.RTOLS  = 1e-9;
\endcode

The dynamic optimization model is solved by passing bounds and an initial guess on the decision variables as follows:

\code
  #include "interval.hpp"

  typedef mc::Interval I;
  I Ip[NP];
  double p0[NP];
  p0[0] = 6.;   Ip[0] = I( 0.1, 10. );
  p0[1] = 0.5;  Ip[1] = I( -1., 1. );
  for( unsigned int is=0; is<NS; is++ ){
    p0[2+is] = 0.5;  Ip[2+is]  = I( -0.5, 1.5 );
  }

  OC->setup();
  Ipopt::ApplicationReturnStatus status = OC->solve( Ip, p0 );
\endcode

The return value is of the enumeration type <A href="http://www.coin-or.org/Doxygen/Ipopt/namespace_ipopt.html#efa0497854479cde8b0994cdf132c982">Ipopt::ApplicationReturnStatus</A>. Moreover, the optimization results can be retrieved as follows:

\code
  #include <iostream>
  #include <iomanip>

  if( status == Ipopt::Solve_Succeeded ){
    std::cout << "OC (LOCAL) SOLUTION: " << std::endl;
    std::cout << "  f* = " << OC->solution().f << std::endl;
    for( unsigned int ip=0; ip<NP; ip++ )
      std::cout << "  p*(" << ip << ") = " << OC->solution().p[ip] << std::endl;
  }
\endcode

The following result is displayed here:

\verbatim
EXIT: Optimal Solution Found.
OC (LOCAL) SOLUTION: 
  f* = 2.35755
  p*(0) = 2.35755
  p*(1) = 1
  p*(2) = 1.5
  p*(3) = 1.5
  p*(4) = 1.5
  p*(5) = 1.5
  p*(6) = 1.5
  p*(7) = 1.5
  p*(8) = 1.5
  p*(9) = 0.475093
  p*(10) = -0.5
  p*(11) = -0.5
\endverbatim
*/

#ifndef MC__DOSEQSLV_IPOPT_HPP
#define MC__DOSEQSLV_IPOPT_HPP

#include <stdexcept>
#include <cassert>
#include "coin/IpTNLP.hpp"
#include "coin/IpIpoptApplication.hpp"

#include "base_do.hpp"
#include "odeslvs_sundials.hpp"
#include "interval.hpp"

#undef  MC__DOSEQSLV_IPOPT_DEBUG
#undef  MC__DOSEQSLV_IPOPT_TRACE

/* TO DO:
- Feasibility and optimality tests (both 1st-order and 2nd-order conditions)
*/

namespace mc
{

//! @brief C++ class for local solution of dynamic optimization using a direct sequential approach with IPOPT, CVODES and MC++
////////////////////////////////////////////////////////////////////////
//! mc::DOSEQSLV_IPOPT is a C++ class for solving dynamic optimization
//! problems to local optimality using a direct sequential approach
//! with IPOPT, CVODES and MC++
////////////////////////////////////////////////////////////////////////
class DOSEQSLV_IPOPT:
  public Ipopt::TNLP,
  public virtual BASE_DO
{
  // Overloading stdout operator
  friend std::ostream& operator<<
    ( std::ostream&os, const DOSEQSLV_IPOPT& );

private:
  //! @brief Internal scaling factor for the objective function (+1: minimize; -1: maximize)
  double _scaling;

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
        delete[] P; P = new std::pair<double,double>[np];
        for( unsigned i=0; i<np; i++ ) P[i] = std::make_pair( Op<T>::l(P_[i]), Op<T>::u(P_[i]) );
      }
  } _data;

  //! @brief Structure holding function evaluation results
  struct EVAL{
    EVAL(): ns(0), f(0), fp(0), sens(false), fstruct(0) {}
    //EVAL(): ns(0), f(0), xk(0), fp(0), xpk(0), sens(false), fstruct(0) {}
    ~EVAL() { cleanup(); }
    void cleanup
      ()
      {
        delete[] f; delete[] fp; delete[] fstruct;
        //for( unsigned is=0; xk && is<=ns; is++ )  delete[] xk[is];
        //for( unsigned is=0; xpk && is<=ns; is++ ) delete[] xpk[is];
        //delete[] xk; delete[] xpk;
      }
    void resize
      ( const unsigned ns_, const unsigned nf, const unsigned np )//, const unsigned nx )
      {
        cleanup();
        ns = ns_;
        //xk = new double*[ ns+1 ];
        //for( unsigned is=0; is<=ns; is++ ) xk[is] = 0;//new double[ nx ];
        f = new double[ nf ];
        //xpk = new double*[ ns+1 ];
        //for( unsigned is=0; is<=ns; is++ ) xpk[is] = 0;//new double[ nx*np ];
        fp = new double[ nf*np ];
        fstruct = new std::map<int,bool>[ nf ];
      }
    unsigned ns;
    double *f;
    //double **xk;
    double *fp;
    //double **xpk;
    bool sens;
    std::map<int,bool>* fstruct;
  } _eval;


  //! @brief Structure holding solution information
  SOLUTION_DO _solution;

  //! @brief IVP solution and sensitivity
  ODESLVS_SUNDIALS _ODESLVS;

public:
  /** @defgroup DOSEQSLV Local Dynamic Optimization using a direct sequential approach with IPOPT, CVODES and MC++
   *  @{
   */
  //! @brief Constructor
  DOSEQSLV_IPOPT()
    { _pDFCTPAR = 0; }

  //! @brief Destructor
  virtual ~DOSEQSLV_IPOPT()
    { delete[] _pDFCTPAR; }

  //! @brief Dynamic optimization options
  struct Options
  {
    //! @brief Constructor
    Options():
      CVTOL(1e-8), PRIMALTOL(1e-8), DUALTOL(1e-4), COMPLTOL(1e-7),
      MAXITER(100), MAXCPU(1e6), GRADIENT(FORWARD), HESSIAN(LBFGS),
      LINSLV(MA57), TESTDER(false), DISPLAY(0), ODESLVS()
      { ODESLVS.DISPLAY = 0; ODESLVS.RESRECORD = 0; }
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
        ODESLVS   = options.ODESLVS;
        return *this;
      }
    //! @brief Enumeration type for Hessian strategy
    enum HESSIAN_STRATEGY{
      LBFGS=0	//!< Perform a limited-memory quasi-Newton approximation
      // Exact Hessian currently unavailable
    };
    //! @brief Enumeration type for gradient strategy
    enum GRADIENT_STRATEGY{
      FORWARD=0,	//!< Forward Sensitivity
      BACKWARD		//!< Backward (Adjoint) Sensitivity
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
    //! @brief IVP solver options
    ODESLVS_SUNDIALS::Options ODESLVS;
  } options;

  //! @brief Dynamic optimization exceptions
  class Exceptions
  {
  public:
    //! @brief Enumeration type for DOSEQSLV_IPOPT exception handling
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
        return "DOSEQSLV_IPOPT::Exceptions  Undefined objective function in optimization problem";
      case INTERN: default:
        return "DOSEQSLV_IPOPT::Exceptions  Internal error";
      }
    }
  private:
    TYPE _ierr;
  };

  //! @brief Setup DAG for function/gradient/Hessian evaluation
  bool setup();

  //! @brief Solve DO model -- return value is IPOPT status
  template <typename T>
  int solve
    ( const T*P, const double*p0=0 );

  //! @brief Test primal feasibility of parameters <a>p</a>
  bool primal_feasible
    ( const double*p, const double FEASTOL );

  //! @brief Test primal feasibility of current solution point
  bool primal_feasible
    ( const double FEASTOL );

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
  const SOLUTION_DO& solution() const
    {
      return _solution;
    }

  //! @brief Access to ODE integrator
  //ODESLVS_SUNDIALS& ODESLVS
  //  ()
  //  { return _ODESLVS; }

  //! @brief Numerical solution of parametric ODEs
  BASE_DE::STATUS states
    ( const double*p, double**xk=0, double*f=0, std::ostream&os=std::cout )
    { _valSTADEP.clear();
      for( auto it=_ndxFCTDEP.begin(); it!=_ndxFCTDEP.end(); ++it )
        _valSTADEP.push_back( p[*it] );
      _ODESLVS.options = options.ODESLVS;
      return _ODESLVS.states( _valSTADEP.data(), xk, f, os ); }

  //! @brief Forward sensitivity of parametric ODEs
  BASE_DE::STATUS states_FSA
    ( const double*p, double**xk=0, double*f=0, double**xpk=0, double*fp=0, std::ostream&os=std::cout )
    { _valSTADEP.clear();
      for( auto it=_ndxFCTDEP.begin(); it!=_ndxFCTDEP.end(); ++it )
        _valSTADEP.push_back( p[*it] );
      _ODESLVS.options = options.ODESLVS;
      return _ODESLVS.states_FSA( _valSTADEP.data(), xk, f, xpk, fp, os ); }

  //! @brief Forward sensitivity of parametric ODEs
  BASE_DE::STATUS states_ASA
    ( const double*p, double**xk=0, double*f=0, double**lk=0, double*fp=0, std::ostream&os=std::cout )
    { _valSTADEP.clear();
      for( auto it=_ndxFCTDEP.begin(); it!=_ndxFCTDEP.end(); ++it )
        _valSTADEP.push_back( p[*it] );
      _ODESLVS.options = options.ODESLVS;
      return _ODESLVS.states_ASA( _valSTADEP.data(), xk, f, lk, fp, os ); }

  //! @brief Record states in file <a>ores</a>, with accuracy of <a>iprec</a> digits
  void record
    ( std::ofstream&ores, const unsigned iprec=5 ) const
    { return _ODESLVS.record( ores, iprec ); }
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

  //! @brief Method to update function values and gradients
  bool _update_functions
    ( const Ipopt::Number *p, const bool sens );

  //! @brief subgraph of parameter-dependent functions
  FFSubgraph _opFCTPAR;

  //! @brief const pointers to parameter-dependent function derivatives
  const FFVar* _pDFCTPAR;

  //! @brief subgraph of parameter-dependent function derivatives
  FFSubgraph _opDFCTPAR;

  //! @brief values of parameters participating in parametric IVP
  std::vector<double> _valSTADEP;

  //! @brief values of state-dependent functions
  std::vector<double> _valFCTSTA;

  //! @brief derivatives of state-dependent functions
  std::vector<double> _valDFCTSTA;

  //! @brief values of parameter-dependent functions
  std::vector<double> _valFCTPAR;

  //! @brief derivatives of parameter-dependent functions
  std::vector<double> _valDFCTPAR;

  //! @brief workspace for parameter-dependent functions
  std::vector<double> _wkFCTPAR;

  //! @brief Private methods to block default compiler methods
  DOSEQSLV_IPOPT(const DOSEQSLV_IPOPT&);
  DOSEQSLV_IPOPT& operator=(const DOSEQSLV_IPOPT&);
};

inline
void
DOSEQSLV_IPOPT::set_options
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
  //case Options::EXACT:
  //  IpoptApp->Options()->SetStringValue( "hessian_approximation", "exact" );
  //  break;
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

  _ODESLVS.options = options.ODESLVS;
}

inline
bool
DOSEQSLV_IPOPT::setup
()
{
  // Generate dependency inormation
  if( !set_depend() ) return false;
  _eval.resize( _ns, _nf, _np );//, _nx );
  for( unsigned ic=0; ic<_nf; ic++ ){
    _eval.fstruct[ic].clear();
    for( unsigned ip=0; ip<_np; ip++ ){
      auto idep = _depF[ic].dep( _pP[ip].id().second );
      if( idep.first ){
        _eval.fstruct[ic][ip] = idep.second;
      }
    }
  }

  // Setup parametric IVP solver
  _ODESLVS.set( *this );
  _ODESLVS.set_time( _ns, _dT.data(), _pT );
  _ODESLVS.set_parameter( _vFCTDEP.size(), _vFCTDEP.data() );
  _ODESLVS.set_function( _ns, _ndxFCTSTA.size(), _vFCTSTA.data() );

  // setup derivatives of parameter-dependent functions
  delete[] _pDFCTPAR;
  switch( options.GRADIENT ){
    case Options::FORWARD:  // Forward AD
      _pDFCTPAR = _pDAG->FAD( _vFCTPAR.size(), _vFCTPAR.data(), _np, _pP ); 
      break;
    case Options::BACKWARD: // Backward AD
      _pDFCTPAR = _pDAG->BAD( _vFCTPAR.size(), _vFCTPAR.data(), _np, _pP ); 
      break;
  }
  _opFCTPAR  = _pDAG->subgraph( _vFCTPAR.size(), _vFCTPAR.data() );
  _opDFCTPAR = _pDAG->subgraph( _vFCTPAR.size()*_np, _pDFCTPAR );
  _wkFCTPAR.resize( _opFCTPAR.l_op.size()>_opDFCTPAR.l_op.size()? _opFCTPAR.l_op.size(): _opDFCTPAR.l_op.size() );
  _valFCTSTA.resize( _ndxFCTSTA.size() );
  _valDFCTSTA.resize( _vFCTDEP.size()*_ndxFCTSTA.size() );
  _valFCTPAR.resize( _vFCTPAR.size() );
  _valDFCTPAR.resize( _np*_vFCTPAR.size() );

#ifdef MC__DOSEQSLV_IPOPT_DEBUG
  {int dum; std::cin >> dum;}
#endif

  return true;
}

template <typename T>
inline
int
DOSEQSLV_IPOPT::solve
( const T*P, const double*p0 )
{
  Ipopt::SmartPtr<Ipopt::IpoptApplication> IpoptApp = new Ipopt::IpoptApplication();

  // Set (a few) IPOPT options
  set_options( IpoptApp );

  // Keep track of bounds and initial guess
  _data.set( _np, p0, P );

  // Run NLP solver
  _solution.status = IpoptApp->Initialize();
  if( _solution.status == Ipopt::Solve_Succeeded )
    _solution.status = IpoptApp->OptimizeTNLP( this );

  return _solution.status;
}

inline
bool
DOSEQSLV_IPOPT::primal_feasible
( const double*p, const double FEASTOL )
{
  _solution.p.assign( p, p+_np );
  return primal_feasible( FEASTOL );
}

inline
bool
DOSEQSLV_IPOPT::primal_feasible
( const double FEASTOL )
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
    std::cout << "g[" << j << "] = " << _solution.g[ic] << "  (" << maxinfeas << ")\n";
#endif
    if( maxinfeas > FEASTOL ) return false;
  }

  return true;
}

inline
bool
DOSEQSLV_IPOPT::get_nlp_info
( Ipopt::Index& n, Ipopt::Index& m, Ipopt::Index& nnz_jac_g,
  Ipopt::Index& nnz_h_lag, IndexStyleEnum& index_style )
{
#ifdef MC__DOSEQSLV_IPOPT_TRACE
    std::cout << "  DOSEQSLV_IPOPT::get_nlp_info\n";
#endif
  // set size
  n = _np; m = _nf-1;
  nnz_jac_g = 0;
  for( unsigned int ic=1; ic<_nf; ic++ )
    nnz_jac_g += _eval.fstruct[ic].size();
  nnz_h_lag = 0;

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
DOSEQSLV_IPOPT::get_bounds_info
( Ipopt::Index n, Ipopt::Number* p_l, Ipopt::Number* p_u,
  Ipopt::Index m, Ipopt::Number* g_l, Ipopt::Number* g_u )
{
#ifdef MC__DOSEQSLV_IPOPT_TRACE
    std::cout << "  DOSEQSLV_IPOPT::get_bounds_info\n";
#endif
  // set variable bounds
  unsigned ip = 0;
  for( ; ip<_np; ip++ ){   
    p_l[ip] = ( _data.P? _data.P[ip].first: -INF );
    p_u[ip] = ( _data.P? _data.P[ip].second: INF );
#ifdef MC__DOSEQSLV_IPOPT_DEBUG
    std::cout << "  p_l[" << ip << "] = " << p_l[ip]
              << "  p_u[" << ip << "] = " << p_u[ip] << std::endl;
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
#ifdef MC__DOSEQSLV_IPOPT_DEBUG
    std::cout << "  g_l[" << ic << "] = " << g_l[ic]
              << "  g_u[" << ic << "] = " << g_u[ic] << std::endl;
#endif
  }
  return true;
}
inline
bool
DOSEQSLV_IPOPT::get_starting_point
( Ipopt::Index n, bool init_p, Ipopt::Number* p, bool init_up,
  Ipopt::Number* up_L, Ipopt::Number* up_U, Ipopt::Index m,
  bool init_ug, Ipopt::Number* ug )
{
#ifdef MC__DOSEQSLV_IPOPT_TRACE
    std::cout << "  DOSEQSLV_IPOPT::get_starting_point  "
              << init_p << init_up << init_ug << std::endl;
#endif
  // Here, we assume we only have starting values for x, if you code
  // your own NLP, you can provide starting values for the dual variables
  // if you wish
  if( !init_p || init_up || init_ug ) return false;

  // initialize to the given starting point
  //if( !_data.p0 ) return false;
  for( unsigned ip=0; ip<_np; ip++ ){
    p[ip] = _data.p0? _data.p0[ip]: 0.;
#ifdef MC__DOSEQSLV_IPOPT_DEBUG
    std::cout << "  p_0[" << ip << "] = " << p[ip] << std::endl;
#endif
  }
  return true;
}

inline
bool
DOSEQSLV_IPOPT::_update_functions
( const Ipopt::Number*p, const bool sens )
{
#ifdef MC__DOSEQSLV_IPOPT_TRACE
  std::cout << "  DOSEQSLV_IPOPT::update_functions " << (sens?"Y":"N") << std::endl;
#endif
#ifdef MC__DOSEQSLV_SHOW_VARITER
  std::cout << "\n" << std::scientific << std::setprecision(15);
  for( unsigned int i=0; i<_np; i++ )
     std::cout << "P[" << i << "] = " << p[i] << std::endl;
  std::cout << "\n";
#endif

  // Values and derivatives of parameter-dependent functions
  try{
    _pDAG->eval( _opFCTPAR, _wkFCTPAR, _vFCTPAR.size(), _vFCTPAR.data(), _valFCTPAR.data(), _np, _pP, p );
    if( sens )
      _pDAG->eval( _opDFCTPAR, _wkFCTPAR, _vFCTPAR.size()*_np, _pDFCTPAR, _valDFCTPAR.data(), _np, _pP, p );
  }
  catch(...){
    return false;
  }

  // Values and derivatives of state-dependent functions
  _valSTADEP.clear();
  for( auto it=_ndxFCTDEP.begin(); it!=_ndxFCTDEP.end(); ++it )
    _valSTADEP.push_back( p[*it] );

  STATUS flag = FATAL;
  if( !sens )
    flag = _ODESLVS.states( _valSTADEP.data(), 0, _valFCTSTA.data() );
  else{
    switch( options.GRADIENT ){
    case Options::FORWARD:
      flag = _ODESLVS.states_FSA( _valSTADEP.data(), 0, _valFCTSTA.data(), 0, _valDFCTSTA.data() );
      break;
    case Options::BACKWARD: default:
      flag = _ODESLVS.states_ASA( _valSTADEP.data(), 0, _valFCTSTA.data(), 0, _valDFCTSTA.data() );
      break;
    }
  }
  if( flag != NORMAL ) return false;

  // Keep track of evaluation results in _eval
  _eval.sens = sens;
  for( unsigned ic=0, icsta=0, icpar=0; ic<_nf; ic++ ){
    // Parameter-dependent functions
    if( _ndxFCTSTA.find( ic ) == _ndxFCTSTA.end() ){
      _eval.f[ic] = _valFCTPAR[icpar];
      for( unsigned ip=0; sens && ip<_np; ip++ )
        _eval.fp[ic+ip*_nf] = _valDFCTPAR[icpar+ip*_vFCTPAR.size()];
      icpar++;
    }
    // State-dependent functions
    else{
      _eval.f[ic] = _valFCTSTA[icsta];
      for( unsigned ip=0, ipsta=0; sens && ip<_np; ip++ ){
        if( _ndxFCTDEP.find( ip ) != _ndxFCTDEP.end() ){
          _eval.fp[ic+ip*_nf] = _valDFCTSTA[icsta+ipsta*_ndxFCTSTA.size()];
          ipsta++;
          continue;
        }
        _eval.fp[ic+ip*_nf] = 0.;
      }
      icsta++;
    }
  }

  return true;
}

inline
bool
DOSEQSLV_IPOPT::eval_f
( Ipopt::Index n, const Ipopt::Number* p, bool new_p,
  Ipopt::Number& f )
{
#ifdef MC__DOSEQSLV_IPOPT_TRACE
    std::cout << "  DOSEQSLV_IPOPT::eval_f  " << new_p << std::endl;
#endif
  if( new_p && !_update_functions( p, false ) ) return false;
  f = _eval.f[0];
#ifdef MC__DOSEQSLV_IPOPT_DEBUG
  std::cout << "  f = " << f << std::endl;
#endif
  return true;
}

inline
bool
DOSEQSLV_IPOPT::eval_grad_f
( Ipopt::Index n, const Ipopt::Number* p, bool new_p,
  Ipopt::Number* fp )
{
#ifdef MC__DOSEQSLV_IPOPT_TRACE
  std::cout << "  DOSEQSLV_IPOPT::eval_grad_f  " << new_p << std::endl;
#endif
  if( (new_p || !_eval.sens) && !_update_functions( p, true ) ) return false;
  for( unsigned ip=0; ip<_np; ip++ ){
    fp[ip] = _eval.fp[ip*_nf];
#ifdef MC__DOSEQSLV_IPOPT_DEBUG
    std::cout << "  fp[" << ip << "] = " << fp[ip] << std::endl;
#endif
  }
  return true;
}

inline
bool
DOSEQSLV_IPOPT::eval_g
( Ipopt::Index n, const Ipopt::Number* p, bool new_p, Ipopt::Index m,
  Ipopt::Number* g )
{
#ifdef MC__DOSEQSLV_IPOPT_TRACE
  std::cout << "  DOSEQSLV_IPOPT::eval_g  " << new_p << std::endl;
#endif
  if( new_p && !_update_functions( p, false ) ) return false;
  for( unsigned ic=1; ic<_nf; ic++ ){
    g[ic-1] = _eval.f[ic];
#ifdef MC__DOSEQSLV_IPOPT_DEBUG
    std::cout << "  g[" << ic-1 << "] = " << g[ic-1] << std::endl;
#endif
  }
  return true;
}

inline
bool
DOSEQSLV_IPOPT::eval_jac_g
( Ipopt::Index n, const Ipopt::Number* p, bool new_p, Ipopt::Index m,
  Ipopt::Index nele_jac, Ipopt::Index* iRow, Ipopt::Index *jCol,
  Ipopt::Number* gp )
{
#ifdef MC__DOSEQSLV_IPOPT_TRACE
  std::cout << "  DOSEQSLV_IPOPT::eval_jac_g  " << new_p << std::endl;
#endif
  // return the constraint Jacobian structure
  if( !gp ){
    unsigned ie = 0;
    for( unsigned ic=1; ic<_nf; ic++ ){
      std::map<int,bool>::const_iterator it = _eval.fstruct[ic].begin();
      for( ; it != _eval.fstruct[ic].end(); ++it, ie++ ){
        iRow[ie] = ic-1;
        jCol[ie] = it->first;
#ifdef MC__DOSEQSLV_IPOPT_DEBUG
        std::cout << "  gp[" << ic-1 << "][" << it->first << "]" << std::endl;
#endif
      }
    }
    assert( (int)ie == nele_jac );
    return true;
  }

  // evaluate constraint Jacobian
  if( (new_p || !_eval.sens) && !_update_functions( p, true ) ) return false;
  unsigned ie = 0;
  for( unsigned ic=1; ic<_nf; ic++ ){
    std::map<int,bool>::const_iterator it = _eval.fstruct[ic].begin();
    for( ; it != _eval.fstruct[ic].end(); ++it, ie++ ){
      gp[ie] = _eval.fp[ic+it->first*_nf]; 
#ifdef MC__DOSEQSLV_IPOPT_DEBUG
      std::cout << "  gp[" << (ic-1)*_np+it->first << "] = " << gp[ie]
  		<< std::endl;
#endif
    }
  }

#ifdef MC__DOSEQSLV_IPOPT_DEBUG
  int dum; std::cin >> dum;
#endif
  return true;
}

inline
bool
DOSEQSLV_IPOPT::eval_h
( Ipopt::Index n, const Ipopt::Number* p, bool new_p,
  Ipopt::Number obj_factor, Ipopt::Index m, const Ipopt::Number* ug,
  bool new_ug, Ipopt::Index nele_hess, Ipopt::Index* iRow,
  Ipopt::Index* jCol, Ipopt::Number* Lpp )
{
#ifdef MC__DOSEQSLV_IPOPT_TRACE
  std::cout << "  DOSEQSLV_IPOPT::eval_h  " << new_p  << new_ug << std::endl;
#endif
  return false;
}

inline
void
DOSEQSLV_IPOPT::finalize_solution
( Ipopt::SolverReturn status, Ipopt::Index n, const Ipopt::Number* p,
  const Ipopt::Number* upL, const Ipopt::Number* upU, Ipopt::Index m,
  const Ipopt::Number* g, const Ipopt::Number* ug, Ipopt::Number f,
  const Ipopt::IpoptData* ip_data, Ipopt::IpoptCalculatedQuantities* ip_cq )
{
#ifdef MC__DOSEQSLV_IPOPT_TRACE
    std::cout << "  DOSEQSLV_IPOPT::finalize_solution\n";
#endif
  _solution.status = status;

  // Successful (or near-successful) completion
  //if( status == Ipopt::SUCCESS || status == Ipopt::STOP_AT_ACCEPTABLE_POINT
  // || status == Ipopt::STOP_AT_TINY_STEP ){
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
DOSEQSLV_IPOPT::intermediate_callback
( Ipopt::AlgorithmMode mode, Ipopt::Index iter, Ipopt::Number obj_value,
  Ipopt::Number inf_pr, Ipopt::Number inf_du, Ipopt::Number mu,
  Ipopt::Number d_norm, Ipopt::Number regularization_size,
  Ipopt::Number alpha_du, Ipopt::Number alpha_pr, Ipopt::Index ls_trials,
  const Ipopt::IpoptData* ip_data, Ipopt::IpoptCalculatedQuantities* ip_cq )
{
#ifdef MC__DOSEQSLV_IPOPT_TRACE
    std::cout << "  DOSEQSLV_IPOPT::intermediate_callback\n";
#endif
  return true;
}

} // end namescape mc

#endif
