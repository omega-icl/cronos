// Copyright (C) 2012, 2014 Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

/*!
\page page_DOSEQ Local Dynamic Optimization using IPOPT, CVODES and FADBAD++ (Sequential Approach)
\author Benoit C. Chachuat
\version 0.1
\date 2011
\bug No known bugs.

Consider a dynamic optimization (DO) problem of the form
\f{align*}
  \min_{{\bf p}\in P}\ & \phi_0({\bf p}, {\bf x}(t_0),\ldots,{\bf x}(t_N)) + \int_{t_0}^{t_N} \psi_0({\bf p}, {\bf x}(t))\, dt\\
  \text{s.t.} \ & \phi_k({\bf p}, {\bf x}(t_0),\ldots,{\bf x}(t_N)) + \int_{t_0}^{t_N} \psi_k({\bf p}, {\bf x}(t))\, dt \{\geq,=,\leq\} 0, \quad k=1,\ldots,n_c\\
  & \dot{{\bf x}}(t)={\bf f}({\bf p},{\bf x}(t)),\ t\in(t_0,t_{N}]; \quad {\bf x}(t_0)={\bf h}({\bf p})
\f}
where \f${\bf p}\in P\subset\mathbb{R}^{n_p}\f$, \f${\bf x}\in\mathbb{R}^{n_x}\f$, and the real-valued functions \f${\bf f}\f$, \f${\bf h}\f$, \f$\phi_k\f$ and \f$\psi_k\f$ are twice continuously differentiable in all their arguments and factorable. The class mc::DOSEQ in MC++ solves such DO problems to local optimality, based on the sequential approach. It relies on the software packages <A href="https://projects.coin-or.org/Ipopt">IPOPT</A> for local solution of NLP problems, <A href="https://computation.llnl.gov/casc/sundials/description/description.html">SUNDIALS/CVODES</A> for numerical integration and sensitivity analysis of IVPs in ODEs, and <A href="http://www.fadbad.com/fadbad.html">FADBAD++</A> for automatic differentiation (AD).

mc::DOSEQ is templated in the DO problem to be solved, which has to be a class derived from mc::DOSTRUCT. For example, suppose we want to solve the following dynamic optimization problem:
\f{align*}
  \max_{T,x_1^0,u_1,\ldots,u_N}\ & T \\
  \text{s.t.} \ & x_1(1) = 1\\
  & x_2(1) = 0\\
  & \left.\begin{array}{l}
    \dot{x_1}(t) = x_2(t)\\
    \dot{x_2}(t) = (u_k-x_1(t)-2*x_2(t))\theta(t)\\
    \dot{\theta}(t) = T
    \end{array}\right\},\ t\in({\textstyle\frac{k-1}{N},\frac{k}{N}}],\ k=1,\ldots,N\\
  & x_1(0) = x_1^0,\ x_2(0) = 0,\ \theta(0) = 0
\f}
with \f$N=10\f$ stages. The following class is defined:
\code
      #include "doseq.h"

      // Number of time stages
      const unsigned int NS = 5;
      // Number of decision variables
      const unsigned int NP = 2+NS;
      // Number of state variables
      const unsigned int NX = 3;
      // Number of constraints
      const unsigned int NC = 2;

      class DO : public mc::DOSTRUCT
      {
      public:
        DYNOPT()
          : mc::DOSTRUCT( NP, NX, NC )
          {}

        // ODEs right-hand side
	template <typename TX, typename TP>
        TX RHS
          ( const unsigned int ix, const TP*p, const TX*x, const unsigned int is )
          {
            assert( ix < nx() );
            switch( ix ){
              case 0: return p[0]*x[1];
              case 1: return p[0]*(p[2+is]-x[0]-2*x[1])*x[2];
              case 2: return p[0];
              default: throw std::runtime_error("invalid size");
            }
          }
      
        // Initial conditions
        template <typename T>
        T IC
          ( const unsigned int ix, const T*p )
          {
            assert( ix < nx() );
            switch( ix ){
              case 0: return p[1];
              case 1: return 0.;
              case 2: return 0.;
              default: throw std::runtime_error("invalid size");
            }
          }

        // Objective functional
        template <typename T>
        std::pair<T,t_OBJ> OBJ
          ( const T*p, T* const*x, const unsigned int ns )
          {
            return std::make_pair( p[0], MIN );
          }

        // Constraint functionals
        template <typename T>
        std::pair<T,t_CTR> CTR
          ( const unsigned int ic, const T*p, T* const*x, const unsigned int ns )
          {
            switch( ic ){
              case 0: return std::make_pair( x[ns][0]-1., EQ );
              case 1: return std::make_pair( x[ns][1]-0., EQ );
              default: throw std::runtime_error("invalid size");
            }
          }
      };
\endcode

\section sec_DOSEQ_solve How to Solve a DO Model using mc::DOSEQ?

Start by defining the time stages \f$t_k\f$, parameter set \f$P\f$ and initial guess \f$p_0\f$:

\code
  double p0[NP], tk[NS+1];
  std::pair<double,double> P[NP];
  p0[0] = 1.;   P[0] = std::make_pair( 0.1, 2. );
  p0[1] = 0.5;  P[1] = std::make_pair( -1., 1. );
  tk[0] = 0.;
  for( unsigned int is=0; is<NS; is++ ){
    P[2+is]  = std::make_pair( -0.5, 1.5 );
    p0[2+is] = 0.;
    tk[1+is] = tk[is]+1./(double)NS;
  }
\endcode

Then, define an mc::DOSEQ object:

\code
  Ipopt::SmartPtr< mc::DOSEQ<DO> > pDO = new mc::DOSEQ<DO>;
\endcode

Note that the class <A href="http://www.coin-or.org/Doxygen/Ipopt/class_ipopt_1_1_smart_ptr.html">Ipopt::SmartPtr</A> is used here to comply with the IPOPT guidelines. 

The DO model is optimized as follows:

\code
  Ipopt::ApplicationReturnStatus status = pDO->solve( NS, tk, P, p0 );
\endcode

where the first abd second arguments are the number and value of the time stages, the third one is the variable bounds, and the last one is the initial guess used by the optimizer (IPOPT). The return value is of the enumeration type <A href="http://www.coin-or.org/Doxygen/Ipopt/namespace_ipopt.html#efa0497854479cde8b0994cdf132c982">Ipopt::ApplicationReturnStatus</A>.

Finally, the optimization results can be retrieved as follows:

\code
  #include <iostream>
  #include <iomanip>

  if( status == Ipopt::Solve_Succeeded ){
    std::cout << "DO (LOCAL) Solution: " << std::endl;
    std::cout << "  f* = " << pDO->solution().f << std::endl;
    for( unsigned int ip=0; ip<NP; ip++ )
      std::cout << "  p*(" << ip << ") = " << pDO->solution().p[ip] << std::endl;
  }
\endcode

The following result is displayed:

\verbatim
EXIT: Optimal Solution Found.
DO (LOCAL) SOLUTION: 
  f* = 1.53979
  p*(0) = 1.53979
  p*(1) = 1
  p*(2) = 1.5
  p*(3) = 1.47932
  p*(4) = -0.5
\endverbatim
*/

#ifndef MC__DOSEQ_HPP
#define MC__DOSEQ_HPP

#include <stdexcept>
#include <cassert>
#include "dostruct.hpp"
#include "mccvodes.hpp"
#include "fadiff.h"
#include "badiff.h"
#include "coin/IpTNLP.hpp"
#include "coin/IpIpoptApplication.hpp"

#undef  MC__DOSEQ_DEBUG
#undef  MC__DOSEQ_TRACE

/* TO DO:
- Documentation
- Feasibility and optimality tests (both 1st-order and 2nd-order conditions)
*/

namespace mc
{

//! @brief C++ class for local solution of dynamic optimization using IPOPT, CVODES and FADBAD++
////////////////////////////////////////////////////////////////////////
//! mc::DOSEQ is a C++ class for solving dynamic optimization problems
//! to local optimality using IPOPT, CVODES and FADBAD++
////////////////////////////////////////////////////////////////////////
template <typename DO>
class DOSEQ: public Ipopt::TNLP, public CVODES<DO>, public virtual DOSTRUCT
{
  // Overloading stdout operator
  template <typename T> friend std::ostream& operator<<
    ( std::ostream&os, const DOSEQ<T>& );

private:
  // AD types in FADBAD++
  typedef fadbad::F< double > t_F;
  typedef fadbad::B< double > t_B;
  typedef fadbad::B< fadbad::F< double > > t_BF;

  //! @brief Structure holding current optimization data
  struct Data{
    const double*p0;
    const std::pair<double,double>*P;
    unsigned int ns;
    const double *tk;  
  } _data;

  //! @brief Structure holding function evaluation results
  struct Evaluation{
    double f;  
    double *g;
    double *fp;  
    double **gp;
    double *Lpp;
    std::map<int,bool>* gstruct;
    std::map<int,bool>* Lstruct;
    double scaling;
  } _evaluation;

  //! @brief Structure holding solution information
  struct Solution{
    Ipopt::SolverReturn status;
    Ipopt::Index n;
    Ipopt::Number*p;
    Ipopt::Number*upL;
    Ipopt::Number*upU;
    Ipopt::Index m;
    Ipopt::Number*g;
    Ipopt::Number*ug;
    Ipopt::Number f;
  } _solution;

public:
  /** @defgroup DOSEQ Local Dynamic Optimization using IPOPT, CVODES and FADBAD++
   *  @{
   */
  //! @brief Constructor
  DOSEQ()
    : DOSTRUCT(DO())
    {
      _evaluation.gstruct = _evaluation.Lstruct = 0;
      _evaluation.g   = new double[_nc];
      _evaluation.fp  = new double[_np];
      _evaluation.gp  = new double*[_nc];
      for( unsigned int ic=0; ic<_nc; ic++ )
        _evaluation.gp[ic]  = new double[_np];
      _evaluation.Lpp = new double[_np*_np];
      _solution.n = _solution.m = 0;
    }

  //! @brief Destructor
  virtual ~DOSEQ()
    {
      delete[] _evaluation.gstruct; delete _evaluation.Lstruct;
      delete[] _evaluation.g; delete[] _evaluation.fp;
        for( unsigned int ic=0; ic<_nc; ic++ )
        delete[] _evaluation.gp[ic];   
      delete[] _evaluation.gp; delete[] _evaluation.Lpp;
      if( _solution.n ){
        delete[] _solution.p; delete[] _solution.upL; delete[] _solution.upU;
      }
      if( _solution.m ){
        delete[] _solution.g; delete[] _solution.ug;
      }
    }

  //! @brief Class for setting up/storing the solver options
  struct Options
  {
    //! @brief Constructor
    Options():
      CVTOL(1e-8), PRIMALTOL(1e-8), DUALTOL(1e-4), COMPLTOL(1e-7),
      MAXITER(100), MAXCPU(1e6), SPARSE(true), GRADIENT(FORWARD),
      HESSIAN(LBFGS), DISPLAY(0)
      {}
    //! @brief Enumeration type for Hessian strategy
    enum HESSIAN_STRATEGY{
      EXACT=0, 	//!< Use exact second derivatives
      LBFGS,	//!< Perform a limited-memory quasi-Newton approximation
    };
    //! @brief Enumeration type for gradient strategy
    enum GRADIENT_STRATEGY{
      FORWARD=0,	//!< Forward Sensitivity
      BACKWARD		//!< Backward (Adjoint) Sensitivity
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
    //! @brief Whether or not to account for the sparsity of the NLP model?
    bool SPARSE;
    //! @brief Strategy for gradient computation
    GRADIENT_STRATEGY GRADIENT;
    //! @brief Strategy for Hessian computation
    HESSIAN_STRATEGY HESSIAN;
    //! @brief Display level in IPOPT
    int DISPLAY;
  } options;

  //! @brief Solve DO model -- return value is status
  Ipopt::ApplicationReturnStatus solve
    ( const unsigned int ns, const double*tk,
      const std::pair<double,double>*P, const double*p0)
    {
      Ipopt::SmartPtr<Ipopt::IpoptApplication> IpoptApp
        = new Ipopt::IpoptApplication();

      // Keep track of bounds, initial guess and time stages
      _data.p0 = p0;
      _data.P  = P;
      _data.tk = tk;
      _data.ns = ns;

      // Set (a few) IPOPT options
      set_options( IpoptApp );

      // Run NLP solver
      Ipopt::ApplicationReturnStatus status = IpoptApp->Initialize();
      if( status == Ipopt::Solve_Succeeded ){
        status = IpoptApp->OptimizeTNLP( this );
      }
      return status;
    }

  //! @brief Get IPOPT internal scaling value
  double get_scaling()
    {
      Structure *pS  = new Structure[_np];
      Structure **xkS = new Structure*[_data.ns+1];
      for( unsigned int is=0; is<=_data.ns; is++ )
        xkS[is] = new Structure[_nx];
      
      // set scaling factor. 1: minimize; -1: maximize
      std::pair<Structure,t_OBJ> f = DO().OBJ( pS, xkS, _data.ns );
      switch( f.second ){
        case MIN: _evaluation.scaling =  1.; break;
        case MAX: _evaluation.scaling = -1.; break;
      }

      delete[] pS;
      for( unsigned int is=0; is<=_data.ns; is++ )
        delete[] xkS[is];
      delete[] xkS;

      return _evaluation.scaling;
    }

  //! @brief Get IPOPT solution info
  const Solution& solution() const
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
    ( Ipopt::SmartPtr<Ipopt::IpoptApplication>&IpoptApp )
    {
      IpoptApp->Options()->SetNumericValue( "tol", options.CVTOL<0.?
        1e-12: options.CVTOL );
      IpoptApp->Options()->SetNumericValue( "dual_inf_tol", options.DUALTOL<=0.?
        1e-12: options.DUALTOL );
      IpoptApp->Options()->SetNumericValue( "constr_viol_tol", options.PRIMALTOL<=0.?
        1e-12: options.PRIMALTOL );
      IpoptApp->Options()->SetNumericValue( "compl_inf_tol", options.COMPLTOL<=0.?
        1e-12: options.COMPLTOL );
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
      IpoptApp->Options()->SetIntegerValue( "print_level",
        (options.DISPLAY<0? 0: (options.DISPLAY>12? 12: options.DISPLAY )) );
    }

  //! @brief Method to update function values and gradients
  bool update_functions
    ( const Ipopt::Number *p );

  //! @brief Method to update lagragian DSO derivatives
  bool update_lagrangian
    ( const Ipopt::Number *p, const Ipopt::Number *ug );

  //! @brief Private methods to block default compiler methods
  DOSEQ(const DOSEQ&);
  DOSEQ& operator=(const DOSEQ&);
};


template <typename DO> inline
bool
DOSEQ<DO>::get_nlp_info
( Ipopt::Index& n, Ipopt::Index& m, Ipopt::Index& nnz_jac_g,
  Ipopt::Index& nnz_h_lag, IndexStyleEnum& index_style )
{
#ifdef MC__DOSEQ_TRACE
    std::cout << "  DOSEQ<DO>::get_nlp_info\n";
#endif
  // set size
  n = _np; m = _nc;

  // determine IVP sparsity
  double t = 0.; // dummy time
  Structure *pS = new Structure[_np];
  for( unsigned int ip=0; ip<_np; ip++ )
    pS[ip].indep(ip);   
  Structure **xkS = new Structure*[_data.ns+1];
  xkS[0] = new Structure[_nx];
  for( unsigned int ix=0; ix<_nx; ix++ ){
    xkS[0][ix] = DO().IC( ix, pS );
#ifdef MC__DOSEQ_DEBUG
    std::cout << "xkS[0][" << ix << "] <- " << xkS[0][ix]
              << std::endl;
#endif
  }
  for( unsigned int is=0; is<_data.ns; is++ ){
    xkS[is+1] = new Structure[_nx];
    for( unsigned int ix=0; ix<_nx; ix++ )
      xkS[is+1][ix] = DO().RHS( ix, pS, xkS[is], t, is ) + xkS[is][ix];
    bool iterate = true;
    while( iterate ){
      iterate = false;
      for( unsigned int ix=0; ix<_nx; ix++ ){
        Structure xiSnew = DO().RHS( ix, pS, xkS[is+1], t, is ) + xkS[is+1][ix];
        if( xkS[is+1][ix] != xiSnew ){
	  xkS[is+1][ix] = xiSnew; iterate = true;
	}
      }
    }
#ifdef MC__DOSEQ_DEBUG
    for( unsigned int ix=0; ix<_nx; ix++ )
      std::cout << "xkS[" << is+1 << "][" << ix << "] <- " << xkS[is+1][ix]
                << std::endl;
#endif
  }

  // determine constraint sparsity
  delete[] _evaluation.gstruct;
  _evaluation.gstruct = new std::map<int,bool>[_nc];
  nnz_jac_g = 0;
  for( unsigned int ic=0; ic<_nc; ic++ ){
    _evaluation.gstruct[ic] = DO().CTR( ic, pS, xkS, _data.ns ).first.dep();
    nnz_jac_g += _evaluation.gstruct[ic].size();
  }

  // determine Lagrangian sparsity
  delete _evaluation.Lstruct;
  Structure L = DO().OBJ( pS, xkS, _data.ns ).first;
#ifdef MC__DOSEQ_DEBUG
  std::cout << "h_lag:" << L << std::endl;
#endif
  for( unsigned int ic=0; ic<_nc; ic++ ){
    L += DO().CTR( ic, pS, xkS, _data.ns ).first;
#ifdef MC__DOSEQ_DEBUG
    std::cout << "h_lag:" << L << std::endl;
#endif
  }
  _evaluation.Lstruct = new std::map<int,bool>( L.dep() ); 
  if( !options.SPARSE )
    nnz_h_lag = (n*(n+1))/2;
  else{
    int nnz_jac_lag = 0;
    std::map<int,bool>::const_iterator cit = _evaluation.Lstruct->begin();
    for( ; cit != _evaluation.Lstruct->end(); ++cit )
      if( !cit->second ) nnz_jac_lag++;
    nnz_h_lag = (nnz_jac_lag*(nnz_jac_lag+1))/2;
  }

  // use the C style indexing (0-based)
  index_style = Ipopt::TNLP::C_STYLE;

#ifdef MC__DOSEQ_DEBUG
  std::cout << "n:" << n << std::endl;
  std::cout << "m:" << m << std::endl;
  std::cout << "nnz_jac_g:" << nnz_jac_g << std::endl;
  std::cout << "nnz_h_lag:" << nnz_h_lag << std::endl;
  std::cout << "h_lag:" << L << std::endl;
#endif

  // clean up
  delete[] pS;
  for( unsigned int is=0; is<=_data.ns; is++ )
    delete[] xkS[is];
  delete[] xkS;

  return true;
}

template <typename DO> inline
bool
DOSEQ<DO>::get_bounds_info
( Ipopt::Index n, Ipopt::Number* p_l, Ipopt::Number* p_u,
  Ipopt::Index m, Ipopt::Number* g_l, Ipopt::Number* g_u )
{
#ifdef MC__DOSEQ_TRACE
    std::cout << "  DOSEQ<DO>::get_bounds_info\n";
#endif
  // set variable bounds
  for( unsigned int ip=0; ip<_np; ip++ ){   
    p_l[ip] = ( _data.P? _data.P[ip].first: -INF );
    p_u[ip] = ( _data.P? _data.P[ip].second: INF );
#ifdef MC__DOSEQ_DEBUG
    std::cout << "  p_l[" << ip << "] = " << p_l[ip]
              << "  p_u[" << ip << "] = " << p_u[ip] << std::endl;
#endif
  }

  // set constraint bounds
  Structure *pS  = new Structure[_np];
  Structure **xkS = new Structure*[_data.ns+1];
  for( unsigned int is=0; is<=_data.ns; is++ )
    xkS[is] = new Structure[_nx];

  for( unsigned int ic=0; ic<_nc; ic++ ){
    switch( DO().CTR( ic, pS, xkS, _data.ns ).second ){
      case EQ: g_l[ic] = g_u[ic] = 0.; break;
      case LE: g_l[ic] = -INF; g_u[ic] = 0.; break;
      case GE: g_l[ic] = 0.; g_u[ic] = INF; break;
    }
#ifdef MC__DOSEQ_DEBUG
    std::cout << "  g_l[" << ic << "] = " << g_l[ic]
              << "  g_u[" << ic << "] = " << g_u[ic] << std::endl;
#endif
  }

  delete[] pS;
  for( unsigned int is=0; is<=_data.ns; is++ )
    delete[] xkS[is];
  delete[] xkS;

  return true;
}

template <typename DO> inline
bool
DOSEQ<DO>::get_starting_point
( Ipopt::Index n, bool init_p, Ipopt::Number* p, bool init_up,
  Ipopt::Number* up_L, Ipopt::Number* up_U, Ipopt::Index m,
  bool init_ug, Ipopt::Number* ug )
{
#ifdef MC__DOSEQ_TRACE
    std::cout << "  DOSEQ<DO>::get_starting_point  "
              << init_p << init_up << init_ug << std::endl;
#endif
  // Here, we assume we only have starting values for x, if you code
  // your own NLP, you can provide starting values for the dual variables
  // if you wish
  if( !init_p || init_up || init_ug ) return false;

  // initialize to the given starting point
  if( !_data.p0 ) return false;
  for( unsigned int ip=0; ip<_np; ip++ ){
    p[ip] = _data.p0[ip];
#ifdef MC__DOSEQ_DEBUG
    std::cout << "  p_0[" << ip << "] = " << p[ip] << std::endl;
#endif
  }
  return true;
}

template <typename DO> inline
bool
DOSEQ<DO>::update_functions
( const Ipopt::Number *p )
{
#ifdef MC__DOSEQ_TRACE
    std::cout << "  DOSEQ<DO>::update_functions" << std::endl;
#endif
  switch( options.GRADIENT ){
  case Options::FORWARD:
    return( this->functions_FSA( _data.ns, _data.tk, p, _evaluation.f,
      _evaluation.fp, _evaluation.g, _evaluation.gp )
      == mc::CVODES<DO>::NORMAL? true: false );
  case Options::BACKWARD: default:
    return( this->functions_ASA( _data.ns, _data.tk, p, _evaluation.f,
      _evaluation.fp, _evaluation.g, _evaluation.gp )
      == mc::CVODES<DO>::NORMAL? true: false );
  }
}

template <typename DO> inline
bool
DOSEQ<DO>::update_lagrangian
( const Ipopt::Number *p, const Ipopt::Number *ug )
{
#ifdef MC__DOSEQ_TRACE
    std::cout << "  DOSEQ<DO>::update_lagrangian" << std::endl;
#endif
  return( this->functions_DSOASA( _data.ns, _data.tk, p, ug, _evaluation.f,
    _evaluation.fp, _evaluation.g, _evaluation.gp, _evaluation.Lpp )
    == mc::CVODES<DO>::NORMAL? true: false );
}

template <typename DO> inline
bool
DOSEQ<DO>::eval_f
( Ipopt::Index n, const Ipopt::Number* p, bool new_p,
  Ipopt::Number& f )
{
#ifdef MC__DOSEQ_TRACE
    std::cout << "  DOSEQ<DO>::eval_f  " << new_p << std::endl;
#endif
  if( new_p && !update_functions( p ) ) return false;
  f = _evaluation.f;
#ifdef MC__DOSEQ_DEBUG
  std::cout << "  f = " << f << std::endl;
#endif
  return true;
}

template <typename DO> inline
bool
DOSEQ<DO>::eval_grad_f
( Ipopt::Index n, const Ipopt::Number* p, bool new_p,
  Ipopt::Number* fp )
{
#ifdef MC__DOSEQ_TRACE
    std::cout << "  DOSEQ<DO>::eval_grad_f  " << new_p << std::endl;
#endif
  if( new_p && !update_functions( p ) ) return false;
  for( unsigned int ip=0; ip<_np; ip++ ){
    fp[ip] = _evaluation.fp[ip];  
#ifdef MC__DOSEQ_DEBUG
    std::cout << "  fp[" << ip << "] = " << fp[ip] << std::endl;
#endif
  }
  return true;
}

template <typename DO> inline
bool
DOSEQ<DO>::eval_g
( Ipopt::Index n, const Ipopt::Number* p, bool new_p, Ipopt::Index m,
  Ipopt::Number* g )
{
#ifdef MC__DOSEQ_TRACE
    std::cout << "  DOSEQ<DO>::eval_g  " << new_p << std::endl;
#endif
  if( new_p && !update_functions( p ) ) return false;
  for( unsigned int ic=0; ic<_nc; ic++ ){
    g[ic] = _evaluation.g[ic];
#ifdef MC__DOSEQ_DEBUG
    std::cout << "  g[" << ic << "] = " << g[ic] << std::endl;
#endif
  }
  return true;
}

template <typename DO> inline
bool
DOSEQ<DO>::eval_jac_g
( Ipopt::Index n, const Ipopt::Number* p, bool new_p, Ipopt::Index m,
  Ipopt::Index nele_jac, Ipopt::Index* iRow, Ipopt::Index *jCol,
  Ipopt::Number* gp )
{
#ifdef MC__DOSEQ_TRACE
    std::cout << "  DOSEQ<DO>::eval_jac_g  " << new_p << std::endl;
#endif
  // return the constraint Jacobian structure
  if( !gp ){
    unsigned int ie = 0;
    for( unsigned int ic=0; ic<_nc; ic++ ){
      std::map<int,bool>::const_iterator it = _evaluation.gstruct[ic].begin();
      for( ; it != _evaluation.gstruct[ic].end(); ++it, ie++ ){
        iRow[ie] = ic;
        jCol[ie] = it->first;
#ifdef MC__DOSEQ_DEBUG
        std::cout << "  gp[" << ic << "][" << it->first << "]" << std::endl;
#endif
      }
    }
    assert( (int)ie == nele_jac );
    return true;
  }

  // evaluate constraint Jacobian
  if( new_p && !update_functions( p ) ) return false;
  unsigned int ie = 0;
  for( unsigned int ic=0; ic<_nc; ic++ ){
    std::map<int,bool>::const_iterator it = _evaluation.gstruct[ic].begin();
    for( ; it != _evaluation.gstruct[ic].end(); ++it, ie++ ){
      gp[ie] = _evaluation.gp[ic][it->first]; 
#ifdef MC__DOSEQ_DEBUG
      std::cout << "  gp[" << ic << "][" << it->first << "] = " << gp[ie]
  		<< std::endl;
#endif
    }
  }

#ifdef MC__DOSEQ_DEBUG
  int dum; std::cin >> dum;
#endif
  return true;
}

template <typename DO> inline
bool
DOSEQ<DO>::eval_h
( Ipopt::Index n, const Ipopt::Number* p, bool new_p,
  Ipopt::Number obj_factor, Ipopt::Index m, const Ipopt::Number* ug,
  bool new_ug, Ipopt::Index nele_hess, Ipopt::Index* iRow,
  Ipopt::Index* jCol, Ipopt::Number* Lpp )
{
#ifdef MC__DOSEQ_TRACE
    std::cout << "  DOSEQ<DO>::eval_h  " << new_p  << new_ug << std::endl;
#endif
  // return the Lagangian Hessian structure (dense)
  if( !Lpp ){
    unsigned int ie = 0;
    // Dense Hessian structure
    if( !options.SPARSE ){
      for( unsigned int ip=0; ip<_np; ip++ ){
        for( unsigned int jp=ip; jp<_np; jp++, ie++ ){
          iRow[ie] = ip; jCol[ie] = jp;
#ifdef MC__DOSEQ_DEBUG
          std::cout << "  Lpp[" << ip << "][" << jp << "]" << std::endl;
#endif
        }
      }
    }
    // Sparse Hessian structure
    else{
      std::map<int,bool>::const_iterator cit = _evaluation.Lstruct->begin();
      for( ; cit != _evaluation.Lstruct->end(); ++cit ){
        if( cit->second ) continue;
        std::map<int,bool>::const_iterator cjt = cit;
        for( ; cjt != _evaluation.Lstruct->end(); ++cjt ){
          if( cjt->second ) continue;
          iRow[ie] = cit->first; jCol[ie] = cjt->first;
#ifdef MC__DOSEQ_DEBUG
          std::cout << "  Lpp[" << cit->first << "][" << cjt->first << "]" << std::endl;
#endif
          ie++;
        }
      }
    }
    assert( (int)ie == nele_hess );
    return true;
  }
  
  // evaluate Lagragian Hessian-vector product
  if( ( new_p || new_ug ) && !update_lagrangian( p, ug ) ) return false;
#ifdef MC__DOSEQ_DEBUG
  switch( CVODES<DO>::options.HESSFORMAT ){
  case CVODES<DO>::Options::HESSFULL:
    for( unsigned int ip=0; ip<_np; ip++ ){
      for( unsigned int jp=0; jp<_np; jp++ )
  	std::cout << "  Lpp[" << ip << "][" << jp << "] = "
                  << _evaluation.Lpp[ip+jp*_np];
  	std::cout << std::endl;
    }
    break;
  case CVODES<DO>::Options::HESSLOWER:
    for( unsigned int ip=0; ip<_np; ip++ ){
      for( unsigned int jp=0; jp<=ip; jp++ )
        std::cout << "  Lpp[" << ip << "][" << jp << "] = "
                  << _evaluation.Lpp[jp*_np-(jp*(jp-1))/2+ip-jp];
        std::cout << std::endl;
    }
    break;
  case CVODES<DO>::Options::HESSUPPER:
    for( unsigned int ip=0; ip<_np; ip++ ){
      for( unsigned int jp=ip; jp<_np; jp++ )
        std::cout << "  Lpp[" << ip << "][" << jp << "] = "
                  << _evaluation.Lpp[jp*(jp+1)/2+ip];
        std::cout << std::endl;
    }
    break;
  }
#endif

  unsigned int ie = 0;
  // Dense Hessian structure
  if( !options.SPARSE ){
    for( unsigned int ip=0; ip<_np; ip++ ){
      for( unsigned int jp=ip; jp<_np; jp++, ie++ ){
        switch( CVODES<DO>::options.HESSFORMAT ){
        case CVODES<DO>::Options::HESSFULL:
          Lpp[ie] = (_evaluation.Lpp[ip+jp*_np]+_evaluation.Lpp[jp+ip*_np])/2.;
          break;
        case CVODES<DO>::Options::HESSLOWER:
          Lpp[ie] = _evaluation.Lpp[ip*_np-(ip*(ip+1))/2+jp];
          break;
        case CVODES<DO>::Options::HESSUPPER:
          Lpp[ie] = _evaluation.Lpp[jp*(jp+1)/2+ip];
          break;
        }
#ifdef MC__DOSEQ_DEBUG
        std::cout << "  Lpp[" << ip << "][" << jp << "] = " << Lpp[ie]
                  << std::endl;
#endif
      }
    }
  }
  // Sparse Hessian structure
  else{
    std::map<int,bool>::const_iterator cit = _evaluation.Lstruct->begin();
    for( ; cit != _evaluation.Lstruct->end(); ++cit ){
      if( cit->second ) continue;
      const unsigned int ip = cit->first;
      std::map<int,bool>::const_iterator cjt = cit;
      for( ; cjt != _evaluation.Lstruct->end(); ++cjt ){
        if( cjt->second ) continue;
        const unsigned int jp = cjt->first;
        switch( CVODES<DO>::options.HESSFORMAT ){
        case CVODES<DO>::Options::HESSFULL:
          Lpp[ie] = (_evaluation.Lpp[ip+jp*_np]+_evaluation.Lpp[jp+ip*_np])/2.;
          break;
        case CVODES<DO>::Options::HESSLOWER:
          Lpp[ie] = _evaluation.Lpp[ip*_np-(ip*(ip+1))/2+jp];
          break;
        case CVODES<DO>::Options::HESSUPPER:
          Lpp[ie] = _evaluation.Lpp[(jp*(jp+1))/2+ip];
          break;
        }

#ifdef MC__DOSEQ_DEBUG
        std::cout << "  Lpp[" << ip << "][" << jp << "] = " << Lpp[ie]
                  << std::endl;
        int dum; std::cin >> dum;
#endif
        ie++;
      }
    }
  }

  return true;
}

template <typename DO> inline
void
DOSEQ<DO>::finalize_solution
( Ipopt::SolverReturn status, Ipopt::Index n, const Ipopt::Number* p,
  const Ipopt::Number* upL, const Ipopt::Number* upU, Ipopt::Index m,
  const Ipopt::Number* g, const Ipopt::Number* ug, Ipopt::Number f,
  const Ipopt::IpoptData* ip_data, Ipopt::IpoptCalculatedQuantities* ip_cq )
{
#ifdef MC__DOSEQ_TRACE
    std::cout << "  DOSEQ<DO>::finalize_solution\n";
#endif
  _solution.status = status;

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

template <typename DO> inline
bool
DOSEQ<DO>::intermediate_callback
( Ipopt::AlgorithmMode mode, Ipopt::Index iter, Ipopt::Number obj_value,
  Ipopt::Number inf_pr, Ipopt::Number inf_du, Ipopt::Number mu,
  Ipopt::Number d_norm, Ipopt::Number regularization_size,
  Ipopt::Number alpha_du, Ipopt::Number alpha_pr, Ipopt::Index ls_trials,
  const Ipopt::IpoptData* ip_data, Ipopt::IpoptCalculatedQuantities* ip_cq )
{
  return true;
}

} // end namescape mc

#endif
