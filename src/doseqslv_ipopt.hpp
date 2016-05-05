// Copyright (C) 2012, 2014 Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

/*!
\page page_DOSEQ Local Dynamic Optimization using IPOPT, CVODES and MC++ (Sequential Approach)
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
      const unsigned NS = 5;
      // Number of decision variables
      const unsigned NP = 2+NS;
      // Number of state variables
      const unsigned NX = 3;
      // Number of constraints
      const unsigned NC = 2;

      class DO : public mc::DOSTRUCT
      {
      public:
        DYNOPT()
          : mc::DOSTRUCT( NP, NX, NC )
          {}

        // ODEs right-hand side
	template <typename TX, typename TP>
        TX RHS
          ( const unsigned ix, const TP*p, const TX*x, const unsigned is )
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
          ( const unsigned ix, const T*p )
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
          ( const T*p, T* const*x, const unsigned ns )
          {
            return std::make_pair( p[0], MIN );
          }

        // Constraint functionals
        template <typename T>
        std::pair<T,t_CTR> CTR
          ( const unsigned ic, const T*p, T* const*x, const unsigned ns )
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
  for( unsigned is=0; is<NS; is++ ){
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
    for( unsigned ip=0; ip<NP; ip++ )
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

#ifndef MC__DOSEQSLV_IPOPT_HPP
#define MC__DOSEQSLV_IPOPT_HPP

#include <stdexcept>
#include <cassert>
#include "coin/IpTNLP.hpp"
#include "coin/IpIpoptApplication.hpp"

#include "base_do.hpp"
#include "odeslvs_sundials.hpp"

#undef  MC__DOSEQSLV_IPOPT_DEBUG
#undef  MC__DOSEQSLV_IPOPT_TRACE

/* TO DO:
- Documentation
- Feasibility and optimality tests (both 1st-order and 2nd-order conditions)
*/

namespace mc
{

//! @brief C++ class for local solution of dynamic optimization using a direct sequential approach with IPOPT, CVODES and MC++
////////////////////////////////////////////////////////////////////////
//! mc::DOSEQSLV is a C++ class for solving dynamic optimization
//! problems to local optimality using a direct sequential approach
//! with IPOPT, CVODES and MC++
////////////////////////////////////////////////////////////////////////
class DOSEQSLV_IPOPT:
  public Ipopt::TNLP,
  public virtual BASE_DO,
  public virtual ODESLVS_SUNDIALS
{
  // Overloading stdout operator
  friend std::ostream& operator<<
    ( std::ostream&os, const DOSEQSLV_IPOPT& );

private:
  //! @brief Internal scaling factor for the objective function (+1: minimize; -1: maximize)
  double _scaling;

  //! @brief vector holding functions
  std::vector<FFVar> _fct;

  //! @brief Structure holding optimization problem data
  struct DATA{
    unsigned ns;
    const double*tk;
    const double*p0;
    std::pair<double,double>*P;
    DATA(): ns(0), tk(0), p0(0), P(0) {}
    ~DATA() { delete[] P; }
    template <typename T> void set
      ( const unsigned ns_, const double*tk_, const unsigned np, const double*p0_,
        const T*P_ )
      {
        ns = ns_;
        tk = tk_;
        p0 = p0_;
        delete[] P; P = new std::pair<double,double>[np];
        for( unsigned i=0; i<np; i++ ) P[i] = std::make_pair( Op<T>::l(P_[i]), Op<T>::u(P_[i]) );
      }
  } _data;

  //! @brief Structure holding function evaluation results
  struct EVAL{
    EVAL(): ns(0), f(0), xk(0), fp(0), xpk(0), sens(false), fstruct(0) {}
    ~EVAL() { cleanup(); }
    void cleanup
      ()
      {
        delete[] f; delete[] fp; delete[] fstruct;
        for( unsigned is=0; xk && is<=ns; is++ )  delete[] xk[is];
        for( unsigned is=0; xpk && is<=ns; is++ ) delete[] xpk[is];
        delete[] xk; delete[] xpk;
      }
    void resize
      ( const unsigned ns_, const unsigned nf, const unsigned np, const unsigned nx )
      {
        cleanup();
        ns = ns_;
        xk = new double*[ ns+1 ];
        for( unsigned is=0; is<=ns; is++ ) xk[is] = new double[ nx ];
        f = new double[ nf ];
        xpk = new double*[ ns+1 ];
        for( unsigned is=0; is<=ns; is++ ) xpk[is] = new double[ nx*np ];
        fp = new double[ nf*np ];
        fstruct = new std::map<int,bool>[ nf ];
      }
    unsigned ns;
    double *f;
    double **xk;
    double *fp;
    double **xpk;
    bool sens;
    std::map<int,bool>* fstruct;
  } _eval;

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

public:
  /** @defgroup DOSEQSLV Local Dynamic Optimization using a direct sequential approach with IPOPT, CVODES and MC++
   *  @{
   */
  //! @brief Constructor
  DOSEQSLV_IPOPT()
    : BASE_DE(), ODESLVS_SUNDIALS()
    {}

  //! @brief Destructor
  virtual ~DOSEQSLV_IPOPT()
    {}

  //! @brief Dynamic optimization options
  struct Options:
    public ODESLVS_SUNDIALS::Options
  {
    //! @brief Constructor
    Options():
      ODESLVS_SUNDIALS::Options(), CVTOL(1e-8), PRIMALTOL(1e-8), DUALTOL(1e-4),
      COMPLTOL(1e-7), MAXITER(100), MAXCPU(1e6), GRADIENT(FORWARD),
      HESSIAN(LBFGS), TESTDER(false), DISPLAY(0)
      {}
    //! @brief Assignment operator
    Options& operator= ( Options&options ){
        ODESLVS_SUNDIALS::Options::operator=(options);
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
      LBFGS=0	//!< Perform a limited-memory quasi-Newton approximation
      // Exact Hessian currently unavailable
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
    //! @brief Strategy for gradient computation
    GRADIENT_STRATEGY GRADIENT;
    //! @brief Strategy for Hessian computation during linesearch
    HESSIAN_STRATEGY HESSIAN;
    //! @brief Whether to test first- and second-derivatives
    bool TESTDER;
    //! @brief Display level in IPOPT
    int DISPLAY;
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
  Ipopt::ApplicationReturnStatus solve
    ( const unsigned ns, const double*tk, const T*P, const double*p0=0 );

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

  //! @brief Method to update function values and gradients
  bool _update_functions
    ( const Ipopt::Number *p, const bool sens );

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
  IpoptApp->Options()->SetStringValue( "derivative_test", options.TESTDER? "second-order": "none" );
  IpoptApp->Options()->SetIntegerValue( "print_level", options.DISPLAY<0? 0: (options.DISPLAY>12? 12: options.DISPLAY ) );

  ODESLVS_SUNDIALS::options = options;
  ODESLVS_SUNDIALS::options.DISPLAY = 0;
  ODESLVS_SUNDIALS::options.RESRECORD = 0;
}

inline
bool
DOSEQSLV_IPOPT::setup
()
{
  if( std::get<0>(_obj).size() != 1 ) return false;
  
  // Set number of functions and stages
  const unsigned nf = 1 + std::get<0>(_ctr).size();
  unsigned ns = std::get<1>(_obj)[0].rbegin()->first+1;
  auto it = std::get<1>(_ctr).begin();
  for( ; it!=std::get<1>(_ctr).end(); ++it ){
    const unsigned nsk = it->rbegin()->first+1;
    if( ns < nsk ) ns = nsk;
  }
  _fct.resize( ns*nf );

  // Set objective stage contributions
  unsigned is = 0;
  auto jt = std::get<1>(_obj)[0].begin();
  for( ; jt!=std::get<1>(_obj)[0].end(); ++jt ){
    for( ; is<jt->first; ++is ) _fct[nf*is] = 0;
    _fct[nf*jt->first] = jt->second; ++is;
  }
  for( ; is<ns; ++is ) _fct[nf*is] = 0;

  // Set constraint stage contributions
  it = std::get<1>(_ctr).begin();
  for( unsigned ic=1; it!=std::get<1>(_ctr).end(); ++it, ++ic ){
    unsigned is = 0;
    auto jt = it->begin();
    for( ; jt!=it->end(); ++jt ){
      for( ; is<jt->first; ++is ) _fct[nf*is+ic] = 0;
      _fct[nf*jt->first+ic] = jt->second; ++is;
    }
    for( ; is<ns; ++is ) _fct[nf*is+ic] = 0;
  }

#ifdef MC__DOSEQSLV_IPOPT_DEBUG
  for( unsigned is=0; is<ns; is++ )
    for( unsigned ic=0; ic<nf; ic++ )
      std::cout << "FCT[" << is << "][" << ic << "] = " << _fct[nf*is+ic] << std::endl;
  int dum;  std::cin >> dum; 
#endif
  set_function( ns, nf, _fct.data() );

  // Resize evaluation arrays
  _eval.resize( ns, nf, _np, _nx );

  // Generate parametric dependence of state
  std::vector<FFVar> Xk( _nx ), G(nf);
  if( !_IC_D_SET() ) return false;
  for( unsigned ix=0; ix<_nx; ix++ ){
    Xk[ix] = _pIC[ix];
#ifdef MC__DOSEQSLV_IPOPT_DEBUG
    std::cout << "X[0][" << ix << "]: " << Xk[ix].dep() << std::endl;
#endif
  }
  for( unsigned is=0; is<ns; ++is ){
    const unsigned pos_rhs  = ( _vRHS.size()<=1?  0: is );
    const unsigned pos_quad = ( _vQUAD.size()<=1? 0: is );
    if( (!is || pos_rhs || pos_quad)
     && !ODESLV_SUNDIALS::_RHS_D_SET( pos_rhs, pos_quad ) )
      return false;
    const FFVar* Fk = _pDAG->compose( _nx, ODESLV_SUNDIALS::_pRHS, _nx, _pX, Xk.data() );
    for( unsigned ix=0; ix<_nx; ix++ ){
      Xk[ix] += Fk[ix];
#ifdef MC__DOSEQSLV_IPOPT_DEBUG
      std::cout << "X[" << is+1 << "][" << ix << "]: " << Xk[ix].dep() << std::endl;
#endif
    }
    delete[] Fk;
    bool iterate = true;
    while( iterate ){
      iterate = false;
      const FFVar* Fk = _pDAG->compose( _nx, ODESLV_SUNDIALS::_pRHS, _nx, _pX, Xk.data() );
      for( unsigned ix=0; ix<_nx; ix++ ){
        FFVar Xki = Xk[ix];
        Xk[ix] += Fk[ix];
        if( Xk[ix].dep() != Xki.dep() ) iterate = true;
#ifdef MC__DOSEQSLV_IPOPT_DEBUG
        std::cout << "X[" << is+1 << "][" << ix << "]: " << Xk[ix].dep() << std::endl;
#endif
      }
      delete[] Fk;
    }
    // Need to repeat multiple times!
    const FFVar* pFCT = _vFCT.at( is );
#ifdef MC__DOSEQSLV_IPOPT_DEBUG
    for( unsigned ic=0; ic<nf; ic++ )
      std::cout << "FCT[" << is << "][" << ic << "] = " << pFCT[ic] << std::endl;
#endif
    const FFVar* Gk = _pDAG->compose( nf, pFCT, _nx, _pX, Xk.data() );
    for( unsigned ic=0; ic<nf; ic++ ){
      if( !is ) G[ic]  = Gk[ic];
      else      G[ic] += Gk[ic];
#ifdef MC__DOSEQSLV_IPOPT_DEBUG
      std::cout << "G[" << is+1 << "][" << ic << "]: " << Gk[ic].dep() << std::endl;
#endif
    }
    delete[] Gk;
  }
  for( unsigned ic=0; ic<nf; ic++ ){
    //_eval.fstruct[ic] = G[ic].dep().dep();
    _eval.fstruct[ic].clear();
    for( unsigned ip=0; ip<_np; ip++ ){
      auto dep = G[ic].dep().dep( _pP[ip].id().second );
      if( dep.first ) _eval.fstruct[ic].insert( std::make_pair( ip, dep.second ) );
    }
  }
#ifdef MC__DOSEQSLV_IPOPT_DEBUG
  {int dum; std::cin >> dum;}
#endif
/*
  // dense functions
  for( unsigned ic=0; ic<nf; ic++ ){
    _eval.fstruct[ic].clear();
    for( unsigned ip=0; ip<_np; ip++ )
      _eval.fstruct[ic].insert( std::make_pair( ip, false ) );
  }
*/
  return true;
}

template <typename T>
inline
Ipopt::ApplicationReturnStatus
DOSEQSLV_IPOPT::solve
( const unsigned ns, const double*tk, const T*P, const double*p0 )
{
  Ipopt::SmartPtr<Ipopt::IpoptApplication> IpoptApp = new Ipopt::IpoptApplication();

  // Set (a few) IPOPT options
  set_options( IpoptApp );

  // Keep track of bounds and initial guess
  _data.set( ns, tk, _np, p0, P );

  // Run NLP solver
  _solution.status = IpoptApp->Initialize();
  if( _solution.status == Ipopt::Solve_Succeeded )
    _solution.status = IpoptApp->OptimizeTNLP( this );

  return _solution.status;
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

//#ifdef MC__NLPSLV_IPOPT_DEBUG
  std::cout << "n:" << n << std::endl;
  std::cout << "m:" << m << std::endl;
  std::cout << "nnz_jac_g:" << nnz_jac_g << std::endl;
  std::cout << "nnz_h_lag:" << nnz_h_lag << std::endl;
//#endif

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
  if( !_data.p0 ) return false;
  for( unsigned ip=0; ip<_np; ip++ ){
    p[ip] = _data.p0[ip];
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
  _eval.sens = sens;

  if( !sens )
    return( states( _data.ns, _data.tk, p, _eval.xk, _eval.f ) == NORMAL? true: false );

  switch( options.GRADIENT ){
  case Options::FORWARD:
    return( states_FSA( _data.ns, _data.tk, p, _eval.xk, _eval.f, _eval.xpk, _eval.fp ) == NORMAL? true: false );
  case Options::BACKWARD: default:
    return( states_ASA( _data.ns, _data.tk, p, _eval.xk, _eval.f, _eval.xpk, _eval.fp ) == NORMAL? true: false );
  }
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
  std::cout << "  f = " << f[0] << std::endl;
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
