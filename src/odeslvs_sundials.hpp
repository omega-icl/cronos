// Copyright (C) 2012-2014 Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

/*!
\page page_MCCVODES Numerical Solution and Sensitivity Analysis of IVPs in ODEs using CVODES and FADBAD++
\author Benoit C. Chachuat
\version 0.1
\date 2011
\bug No known bugs.

mc::CVODES is a template C++ class for the numerical integration and sensitivity analysis of parametric IVPs in ODEs,
\f[\dot{{\bf x}}(t)={\bf f}({\bf p},{\bf x}(t)),\ t\in(t_0,t_{N}]; \quad {\bf x}(t_0)={\bf h}({\bf p}), \f]
where \f${\bf p}\in P\subset\mathbb{R}^{n_p}\f$, \f${\bf x}\in\mathbb{R}^{n_x}\f$, and \f${\bf f}: \mathbb{R}^{n_p}\times\mathbb{R}^{n_x}\to\mathbb{R}^{n_x}\f$ and \f${\bf h}:\mathbb{R}^{n_p}\to\mathbb{R}^{n_x}\f$ are twice continuously differentiable.

mc::CVODES also supports the evaluation and sensitivity analysis of functionals depending on the IVP solutions at a finite number of points \f$t_0,\ldots,t_N\f$, with \f$N\geq 1\f$, for instance the cost and constraints in a dynamic optimization problem such as
\f{align*}
  \min_{\bf p}\ & \phi_0({\bf p}, {\bf x}(t_0),\ldots,{\bf x}(t_N)) + \int_{t_0}^{t_N} \psi_0({\bf p}, {\bf x}(t))\, dt\\
  \text{s.t.} \ & \phi_k({\bf p}, {\bf x}(t_0),\ldots,{\bf x}(t_N)) + \int_{t_0}^{t_N} \psi_k({\bf p}, {\bf x}(t))\, dt \{\geq,=,\leq\} 0, \quad k=1,\ldots,n_c\\
  & \dot{{\bf x}}(t)={\bf f}({\bf p},{\bf x}(t)),\ t\in(t_0,t_{N}]; \quad {\bf x}(t_0)={\bf h}({\bf p})
\f}
where the functions \f$\phi_k\f$ and \f$\psi_k\f$ are twice continuously differentiable in all their arguments.

The class mc::CVODES is essentially a wrapper to the solver CVODES in the software package <A href="https://computation.llnl.gov/casc/sundials/description/description.html">SUNDIALS</A>, which implement both forward and backward sensitivity analysis for IVPs in ODEs. To perform optimally CVODES requires the first and second derivatives of the RHS and IC of the ODE model as well as those of the objective and constraint functions in the dynamic optimization model. This information is generated automatically in mc::CVODES using <A href="http://www.fadbad.com/fadbad.html">FADBAD++</A>, which implements both the forward and reverse mode of automatic differentiation (AD).

mc::CVODES is templated in the DO problem to be solved, which has to be a class derived from mc::DOSTRUCT. For example, suppose we want to solve the following dynamic optimization problem:
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
      #include "mccvodes.h"

      // Number of time stages
      const unsigned int NS = 10;
      // Number of decision variables
      const unsigned int NP = 2+NS;
      // Number of state variables
      const unsigned int NX = 3;
      // Number of constraints
      const unsigned int NC = 2;

      class DO : public virtual mc::DOSTRUCT
      {
      public:
        DO()
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

\section sec_CVODES_NOM How to Calculate a Numerical Solution to the IVP and Function Values?



*/

#ifndef MC__ODESLVS_SUNDIALS_HPP
#define MC__ODESLVS_SUNDIALS_HPP

#include <stdexcept>
#include <cassert>
#include <cmath>
#include <fstream>
#include "cvodes.h"
#include "nvector_serial.h"
#include "cvodes_dense.h"
#include "sundials_dense.h"
#include "sundials_types.h"

#include "base_de.hpp"

#undef  MC__ODESLVS_SUNDIALS_DEBUG

/* TO DO:
- Finish writing up documentation
- Create function to check DSOA via finite differences
- Remove static members in initialization functions
*/

namespace mc
{
//! @brief C++ class computing solutions of parametric ODEs using SUNDIALS-CVODES
////////////////////////////////////////////////////////////////////////
//! mc::ODESLVS_SUNDIALS is a C++ class for solution of IVPs in ODEs
//! using SUNDIALS-CVODES. Besides computing a numerical solution of
//! the IVP, it can also perform first- and directional second-order
//! sensitivity analysis of state functionals in the Mayer form.
////////////////////////////////////////////////////////////////////////
class ODESLVS_SUNDIALS: public virtual BASE_DE
{
private:

  //! @brief Pointer to the CVODE memory block (adjoint case)
  void *_cvode_mem;

  //! @brief Pointer to the CVODE memory block (adjoint sensitivity case)
  void *_cvode_memS;

  //! @brief N_Vector object holding current states
  N_Vector _Nx;

  //! @brief N_Vector object holding current state quadratures
  N_Vector _Nxq;

  //! @brief array of N_Vector objects holding current state sensitivities
  N_Vector *_Nxp;

  //! @brief array of N_Vector objects holding current adjoints
  N_Vector *_Nl;

  //! @brief array of N_Vector objects holding current adjoint quadratures
  N_Vector *_Nlq;

  //! @brief N_Vector object holding current state sensitivities for parameter _iparBS
  N_Vector *_Nxpi;

  //! @brief array of N_Vector objects holding current adjoint sensitivities for parameter _iparBS
  N_Vector _Nlpi;

  //! @brief array of N_Vector objects holding current adjoint sensitivity quadratures for parameter _iparBS
  N_Vector _Nqpi;

  //! @brief pointer to array holding current states
  double *_x;

  //! @brief pointer to t_F array holding current states in forward AD
  //t_F *_Fx;

  //! @brief pointer to t_B array holding current states in backward AD
  //t_B *_Bx;

  //! @brief pointer to t_BF array holding current states in forward-backward AD
  //t_BF *_BFx;

  //! @brief pointer to array holding current adjoints
  double *_l;

  //! @brief pointer to array holding current adjoint sensitivities
  double *_lp;

  //! @brief pointer to array holding current parameters -- do NOT free!
  const double *_p;

  //! @brief pointer to t_F array holding current parameters in forward AD
  //t_F *_Fp;

  //! @brief pointer to t_B array holding current parameters in backward AD
  //t_B *_Bp;

  //! @brief pointer to t_BF array holding current parameters in forward-backward AD
  //t_BF *_BFp;

  //! @brief current time
  realtype _t;

  //! @brief current stage
  unsigned int _istg;

  //! @brief pointer to array holding identifiers of the backward problems
  int* _indexB;

  //! @brief identifier of the backward problems with sensitivities
  int _indexBS;

  //! @brief parameter index for current BS problem
  unsigned int _iparBS;

  //! @brief Return flag for SUNDIALS methods
  int _flag;

  //! @brief number of variables for DAG evaluation
  unsigned _nVAR;

  //! @brief pointer to variables for DAG evaluation
  FFVar* _pVAR;

  //! @brief pointer to variable values for DAG evaluation
  realtype* _dVAR;

  //! @brief list of operations in RHS evaluation
  std::list<const FFOp*> _opRHS;

  //! @brief list of operations in RHS quadrature
  std::list<const FFOp*> _opRHSQ;

  //! @brief list of operations in RHS Jacobian
  std::list<const FFOp*> _opJAC;

  //! @brief list of operations in IC evaluation
  std::list<const FFOp*> _opIC;

  //! @brief pointer to RHS function in current stage of ODE system
  const FFVar* _pRHS;

  //! @brief pointer to RHS quadrature in current stage of ODE system
  const FFVar* _pRHSQ;

  //! @brief preallocated array for RHS evaluation
  realtype* _dRHS;

  //! @brief RHS Jacobian in current stage of ODE system (sparse format)
  std::tuple< unsigned, const unsigned*, const unsigned*, const FFVar* > _pJAC;

  //! @brief preallocated array for RHS Jacobian evaluation
  realtype* _dJAC;

  //! @brief const pointer to IC function in current stage of ODE system
  const FFVar* _pIC;


public:
  /** @defgroup ODESLV_SUNDIALS Real-valued (non-validated) integration of parametric ODEs using CVODES-SUNDIALS
   *  @{
   */
  //! @brief Constructor
  ODESLVS_SUNDIALS()
    : BASE_DE(), _nVAR(0), _pVAR(0), _dVAR(0), _pRHS(0), _pRHSQ(0), _dRHS(0),
      _pJAC(0,0,0,0), _dJAC(0), _pIC(0)
    {
      // Size state arrays
      _x  = new double[ _nx ];
      _l  = new double[ _nx ];
      _lp  = new double[ 2*_nx ];
      _Nx  = 0;//N_VNew_Serial( _nx );
      _Nxq = 0;//N_VNew_Serial( _nq );
      _Nxp = 0;//N_VCloneVectorArray_Serial( _np, _Nx );
      _Nl  = 0;//N_VCloneVectorArray_Serial( _nf, _Nx );
      _Nlq = 0;//N_VCloneVectorArrayEmpty_Serial( _nf, _Nx );
      //for( unsigned int i=0; i<_nf; i++ ){
        //N_VDestroy_Serial( _Nlq[i] );
	//_Nlq[i] = N_VNew_Serial( _np );
      //}
      _Nxpi = 0;//N_VCloneVectorArray_Serial( 1, _Nx );
      _Nlpi = 0;//N_VNew_Serial( 2*_nx );
      _Nqpi = 0;//N_VNew_Serial( _np );

      // Create the solver memory and specify the BDF method and the use of a
      // Newton iteration
      _cvode_mem = CVodeCreate( CV_BDF, CV_NEWTON );
      if( _check_flag((void *)_cvode_mem, "CVodeCreate", 0) ) return;
      _cvode_memS = CVodeCreate( CV_BDF, CV_NEWTON );
      if( _check_flag((void *)_cvode_memS, "CVodeCreate", 0) ) return;

      // Create the array of identifiers for the backward solver
      _indexB  = new int[ _nf ];
      _indexBS = 0;
    }

  //! @brief Default destructor
  virtual ~ODESLVS_SUNDIALS()
    {
      // Free integrators memory
      CVodeFree( &_cvode_mem );
      CVodeFree( &_cvode_memS );

      // Free state arrays
      if( _Nx )   N_VDestroy_Serial( _Nx );
      if( _Nxq )  N_VDestroy_Serial( _Nxq );
      if( _Nxp )  N_VDestroyVectorArray_Serial( _Nxp, _np );
      if( _Nl )   N_VDestroyVectorArray_Serial( _Nl, _nf );
      if( _Nlq )  N_VDestroyVectorArray_Serial( _Nlq, _nf );
      if( _Nxpi ) N_VDestroyVectorArray_Serial( _Nxpi, 1 );
      if( _Nlpi ) N_VDestroy_Serial( _Nlpi );
      if( _Nqpi ) N_VDestroy_Serial( _Nqpi );
      delete [] _x;
      //delete [] _Fx;
      //delete [] _Bx;
      //delete [] _BFx;
      delete [] _l;
      delete [] _lp;
      //delete [] _Fp;
      //delete [] _Bp;
      //delete [] _BFp;
      delete [] _indexB;

      delete[] _dVAR;
      delete[] _pVAR;
      /* DO NOT FREE _pIC, _pRHS, _pRHSQ */
      delete[] _dRHS;
      delete[] _dJAC;
      delete[] std::get<1>(_pJAC);  std::get<1>(_pJAC) = 0;
      delete[] std::get<2>(_pJAC);  std::get<2>(_pJAC) = 0;
      delete[] std::get<3>(_pJAC);  std::get<3>(_pJAC) = 0;
    }

  //! @brief Structure for setting up storing the solver options
  struct Options
  {
    //! @brief Constructor
    Options():
      RELTOL(1e-7), ABSTOL(1e-8), QERR(true), MAXFAIL(10), AUTOTOLS(true),
      RELTOLS(1e-7), ABSTOLS(1e-8), FSACORR(STAGGERED), FSAERR(true),
      RELTOLB(1e-7), ABSTOLB(1e-8), QERRB(true), ASAINTERP(HERMITE),
      ASACHKPT(100), HESSFORMAT(HESSLOWER), FDRTOL(1e-3), FDATOL(1e-3),
      FDCEN(true), DISPLAY(0)
      {}
    //! @brief Enumeration type for FSA method
    enum FSA_STRATEGY{
      SIMULTANEOUS=CV_SIMULTANEOUS,//!< Simultaneous state/sensitivity correction
      STAGGERED=CV_STAGGERED,	//!< Simultaneous sensitivity corrections after state corrections
      STAGGERED1=CV_STAGGERED1	//!< Sequential sensitivity corrections after state corrections
    };
    //! @brief Enumeration type for ASA method
    enum ASA_STRATEGY{
      HERMITE=CV_HERMITE,	//!< Cubic Hermite interpolation
      POLYNOMIAL=CV_POLYNOMIAL	//!< Variable degree polynomial interpolation
    };
    //! @brief Enumeration type for DSOA result format
    enum DSOA_FORMAT{
      HESSFULL=0,		//!< Full np x np Hessian matrix
      HESSUPPER,		//!< Upper triangular Hessian matrix (columnwise storage)
      HESSLOWER			//!< Lower triangular Hessian matrix (columnwise storage)
    };
    //! @brief Relative (scalar) integration tolerance
    realtype RELTOL;
    //! @brief Absolute (scalar) integration tolerance
    realtype ABSTOL;
    //! @brief Whether or not (state) quadrature error control is performed?
    bool QERR;
    //! @brief Maximum number of error test failures
    int MAXFAIL;
    //! @brief Whether integration tolerances for FS are to be set automatically?
    bool AUTOTOLS;
    //! @brief Relative (scalar) integration tolerance for FS
    realtype RELTOLS;
    //! @brief Absolute (scalar) integration tolerance for FS
    realtype ABSTOLS;
    //! @brief FSA correction strategy
    FSA_STRATEGY FSACORR;
    //! @brief Whether or not FSA error control is performed?
    bool FSAERR;
    //! @brief Relative (scalar) integration tolerance for AS
    realtype RELTOLB;
    //! @brief Absolute (scalar) integration tolerance for AS
    realtype ABSTOLB;
    //! @brief Whether or not ASA quadrature error control is performed?
    bool QERRB;
    //! @brief ASA interpolation strategy
    ASA_STRATEGY ASAINTERP;
    //! @brief Number of steps between each check point for AS
    int ASACHKPT;
    //! @brief Format of Lagrangian hessian matrix
    DSOA_FORMAT HESSFORMAT;
    //! @brief Relative tolerance for finite differences
    double FDRTOL;
    //! @brief Absolute tolerance for finite differences
    double FDATOL;
    //! @brief Whether or not centered finite differences are used?
    bool FDCEN;
    //! @brief Display level
    int DISPLAY;
  } options;

  //! @brief static pointer to this class
  static ODESLVS_SUNDIALS *pODESLVS;

  //! @brief Numerical solution of IVP in ODEs over <a>ns</a> stages -- return value is status
  STATUS states
    ( const unsigned int ns, const double*tk, const double*p, double**xk,
      double*q, double*f, std::ostream&os=std::cout );
/*
  //! @brief Compute function values using FSA -- return value is status
  STATUS functions
    ( const unsigned int ns, const double*tk, const double*p, double&f,
      double*g, std::ostream&os=std::cout );

  //! @brief Compute function values using FSA -- return value is status
  STATUS feasibility_test
    ( const unsigned int ns, const double*tk, const double*p, double&f,
      double&maxinfeas, std::ostream&os=std::cout );

  //! @brief Numerical FS analysis of IVP in ODEs over <a>ns</a> stages -- return value is status
  STATUS states_FSA
    ( const unsigned int ns, const double*tk, const double*p, double**xk,
      double**xpk, std::ostream&os=std::cout );

  //! @brief Compute function derivatives using FSA -- return value is status
  STATUS functions_FSA
    ( const unsigned int ns, const double*tk, const double*p, double&f,
      double*fp, double*g, double**gp, std::ostream&os=std::cout );

  //! @brief Numerical AS analysis of IVP in ODEs over <a>ns</a> stages -- return value is status
  STATUS states_ASA
  ( const unsigned int ns, const double*tk, const double*p, double**xk,
    double**lk, double*qk, std::ostream&os=std::cout );

  //! @brief Compute function derivatives using ASA -- return value is status
  STATUS functions_ASA
    ( const unsigned int ns, const double*tk, const double*p, double&f,
      double*fp, double*g, double**gp, std::ostream&os=std::cout );

  //! @brief Numerical DSO AS analysis of IVP in ODEs over <a>ns</a> stages -- return value is status
  STATUS states_DSOASA
  ( const unsigned int ns, const double*tk, const double*p, const double*mu,
    double**xk, double**xkp, double**lpk, double*qpk, std::ostream&os=std::cout );

  //! @brief Compute function and Lagragian derivatives using DSOASA -- return value is status
  STATUS functions_DSOASA
    ( const unsigned int ns, const double*tk, const double*p, const double*mu,
      double&f, double*fp, double*g, double**gp, double*Lpp,
      std::ostream&os=std::cout );

  //! @brief Numerical finite differences of IVP in ODEs over <a>ns</a> stages -- return value is status
  STATUS states_FDSA
    ( const unsigned int ns, const double*tk, const double*p, double**xk,
      double**xpk, std::ostream&os=std::cout );

  //! @brief Compute function derivatives using finite differences -- return value is status
  STATUS functions_FDSA
    ( const unsigned int ns, const double*tk, const double*p, double&f,
      double*fp, double*g, double**gp, std::ostream&os=std::cout );

  //! @brief Compute function and Lagragian derivatives using finite differences -- return value is status
  STATUS functions_FDSOASA
    ( const unsigned int ns, const double*tk, const double*p, const double*mu,
      double&f, double*fp, double*g, double**gp, double*Lpp,
      std::ostream&os=std::cout );
*/
private:

  //! @brief Function to initialize CVode for forward integration
  bool _CVode_init();
/*
  //! @brief Function to initialize CVode for ODEs FSA
  STATUS _CVode_FSA_init();

  //! @brief Function to initialize CVode for forward-backward integration
  STATUS _CVode_initB();

  //! @brief Function to initialize CVode for forward-backward integration with sensitivies
  STATUS _CVode_initBS();

  //! @brief Function to initialize CVode for ODEs ASA
  STATUS _CVode_ASA_init
    ( const double tf );

  //! @brief Function to initialize CVode for ODEs DSOFSA
  STATUS _CVode_DSOFSA_init
    ( const double t0 );

  //! @brief Function to initialize CVode for ODEs DSOASA
  STATUS _CVode_DSOASA_init
    ( const double tf );
*/
  //! @brief Function to initialize ODEs solution
  void _INI
    ( const double*p );

  //! @brief Evaluate state initial conditions
  bool _IC();

  //! @brief Evaluate state continuity conditions
  bool _CC();

  //! @brief Set state RHS pointer and corresponding Jacobian
  bool _SET
    ( const unsigned iRHS );

  //! @brief Evaluate the functions at intermediate/end point
  bool _FCT
    ( const unsigned iFCT, double*f );
/*
  //! @brief Function to initialize ODEs FSA
  void _states_FSA_init
    ( const double*p );

  //! @brief Function to initialize ODEs ASA
  void _states_ASA_init
    ( const unsigned int ns, const double* const*xk );

  //! @brief Function to initialize ODEs DSOFSA
  void _states_DSOFSA_init
    ( const double*p );

  //! @brief Function to initialize ODEs DSOASA
  void _states_DSOASA_init
    ( const unsigned int ns, const double*mu, const double* const*xk,
      const double* const*xkp );
*/
  //! @brief Static wrapper to function to calculate the ODEs RHS values
  static int MC_CVRHS__
    ( realtype t, N_Vector Nx, N_Vector Nxdot, void *user_data );

  //! @brief Function to calculate the ODEs RHS values
  int _RHS
    ( realtype t, N_Vector Nx, N_Vector Nxdot, void *user_data );

  //! @brief Static wrapper to function to calculate the ODEs RHS quadratures
  static int MC_CVRHSQ__
    ( realtype t, N_Vector Nx, N_Vector Nqdot, void *user_data );

  //! @brief Function to calculate the ODEs RHS quadratures
  int _RHSQ
    ( realtype t, N_Vector Nx, N_Vector Nqdot, void *user_data );

  //! @brief Static wrapper to function to calculate the ODEs RHS derivatives
  static int MC_CVJAC__
    ( long int N, realtype t, N_Vector y, N_Vector fy, DlsMat Jac, void *user_data,
      N_Vector tmp1, N_Vector tmp2, N_Vector tmp3 );

  //! @brief Function to calculate the ODEs RHS derivatives
  int _JAC
    ( long int N, realtype t, N_Vector y, N_Vector fy, DlsMat Jac, void *user_data );
/*
  //! @brief Static wrapper to function to calculate the ODEs RHS derivatives
  static int MC_CVFSA__
    ( int Ns, realtype t, N_Vector y, N_Vector ydot, int iS, N_Vector yS,
      N_Vector ySdot, void *user_data, N_Vector tmp1, N_Vector tmp2 );

  //! @brief Function to calculate the ODEs RHS derivatives
  int _FSA
    ( int Ns, realtype t, N_Vector y, N_Vector ydot, int iS, N_Vector yS,
      N_Vector ySdot, void *user_data );

  //! @brief Static wrapper to function to calculate the ODEs RHSB values
  static int MC_CVRHSB__
    ( realtype t, N_Vector Nx, N_Vector NxB, N_Vector NxBdot, void *user_data );

  //! @brief Function to calculate the ODEs RHSB values
  int _RHSB
    ( realtype t, N_Vector Nx, N_Vector NxB, N_Vector NxBdot, void *user_data );

  //! @brief Static wrapper to function to calculate the ODEs RHSB derivatives
  static int MC_CVJACB__
    ( long int NB, realtype t, N_Vector y, N_Vector yB, N_Vector fyB, DlsMat JacB,
      void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3 );

  //! @brief Function to calculate the ODEs RHSB derivatives
  int _JACB
    ( long int NB, realtype t, N_Vector y, N_Vector yB, N_Vector fyB, DlsMat JacB,
      void *user_data );

  //! @brief Static wrapper to function to calculate the ODEs RHSB quadratures
  static int MC_CVRHSBQ__
    ( realtype t, N_Vector Nx, N_Vector NxB, N_Vector NqBdot, void *user_data );

  //! @brief Function to calculate the ODEs RHSB quadratures
  int _RHSBQ
    ( realtype t, N_Vector Nx, N_Vector NxB, N_Vector NqBdot, void *user_data );

  //! @brief Static wrapper to function to calculate the ODEs RHS derivatives
  static int MC_CVDSOFSA__
    ( int Ns, realtype t, N_Vector y, N_Vector ydot, N_Vector *yS,
      N_Vector *ySdot, void *user_data, N_Vector tmp1, N_Vector tmp2 );

  //! @brief Function to calculate the ODEs RHS derivatives
  int _DSOFSA
    ( int Ns, realtype t, N_Vector y, N_Vector ydot, N_Vector *yS,
      N_Vector *ySdot, void *user_data );

  //! @brief Static wrapper to function to calculate the ODEs RHSBS values
  static int MC_CVRHSBS__
    ( realtype t, N_Vector Nx, N_Vector *NxS, N_Vector NxB, N_Vector NxBdot,
      void *user_data );

  //! @brief Function to calculate the ODEs RHSBS values
  int _RHSBS
    ( realtype t, N_Vector Nx, N_Vector *NxS, N_Vector NxB, N_Vector NxBdot,
      void *user_data );

  //! @brief Static wrapper to function to calculate the ODEs RHSBS quadratures
  static int MC_CVRHSBSQ__
    ( realtype t, N_Vector Nx, N_Vector *NxS, N_Vector NxB, N_Vector NqBdot,
      void *user_data );

  //! @brief Function to calculate the ODEs RHSBS quadratures
  int _RHSBSQ
    ( realtype t, N_Vector Nx, N_Vector *NxS, N_Vector NxB, N_Vector NqBdot,
      void *user_data );
*/
  //! @brief Function to check function return values in SUNDIALS
  static bool _check_flag
    ( void *flagvalue, std::string funcname, int opt );

  //! @brief Function to display final statistics
  static void _print_final_stats
    ( void *cvode_mem, bool sens=false, std::ostream&os=std::cout );

  //! @brief Function to display intermediate results
  template<typename U> static void _print_interm
    ( const double t, const unsigned nx, const U*x, const std::string&var,
      std::ostream&os=std::cout );

  //! @brief Function to display intermediate results
  template<typename U> static void _print_interm
    ( const unsigned nx, const U*x, const std::string&var,
      std::ostream&os=std::cout );

  //! @brief Private methods to block default compiler methods
  ODESLVS_SUNDIALS(const ODESLVS_SUNDIALS&);
  ODESLVS_SUNDIALS& operator=(const ODESLVS_SUNDIALS&);
};

ODESLVS_SUNDIALS* ODESLVS_SUNDIALS::pODESLVS = 0;

inline bool
ODESLVS_SUNDIALS::_CVode_init
()
{
  // Initialize the integrator memory and specify the user's right hand
  // side function in y'=rhs(t,y), the inital time T0, and the initial
  // dependent variable vector y
  static unsigned nx = 0, nq = 0;
  if( nx != _nx ){
    _flag = CVodeInit( _cvode_mem, MC_CVRHS__, _t, _Nx );
    if( _check_flag(&_flag, "CVodeInit", 1) ) return false;
    nx = _nx;
  }
  else{
    _flag = CVodeReInit( _cvode_mem, _t, _Nx );
    if( _check_flag(&_flag, "CVodeReInit", 1) ) return false;
  }

  // Specify the relative and absolute tolerances for states
  _flag = CVodeSStolerances( _cvode_mem, options.RELTOL, options.ABSTOL );
  if( _check_flag(&_flag, "CVodeSStolerances", 1) ) return false;

  // Specify the CVDENSE dense linear solver
  //_flag = CVLapackDense( _cvode_mem, _nx );
  _flag = CVDense( _cvode_mem, _nx );
  if( _check_flag(&_flag, "CVDense", 1)) return false;

  // Set the Jacobian routine
  _flag = CVDlsSetDenseJacFn( _cvode_mem, MC_CVJAC__ );
  if ( _check_flag(&_flag, "CVDlsSetDenseJacFn", 1) ) return false;

  // Set maximum number of error test failures
  _flag = CVodeSetMaxErrTestFails( _cvode_mem, options.MAXFAIL );
  if ( _check_flag(&_flag, "CVodeSetMaxErrTestFails", 1) ) return false;

  // Initialize the integrator memory for the quadrature variables
  if( !_nq ){
    CVodeQuadFree( _cvode_mem );
    nq = 0;
    return true;
  }
  else if( nq != _nq ){
    _flag = CVodeQuadInit( _cvode_mem, MC_CVRHSQ__, _Nxq );
    if( _check_flag(&_flag, "CVodeQuadInit", 1) ) return false;
    nq = _nq;
  }
  else{
    _flag = CVodeQuadReInit( _cvode_mem, _Nxq );
    if( _check_flag(&_flag, "CVodeQuadReInit", 1) ) return false;
  }

  // Specify whether or not to perform error control on quadrature
  _flag = CVodeSetQuadErrCon( _cvode_mem, options.QERR );
  if( _check_flag(&_flag, "CVodeSetQuadErrCon", 1) ) return false;

  // Specify the relative and absolute tolerances for quadratures
  _flag = CVodeQuadSStolerances( _cvode_mem, options.RELTOL, options.ABSTOL );
  if( _check_flag(&_flag, "CVodeQuadSStolerances", 1) ) return false;

  return true;
}

inline void
ODESLVS_SUNDIALS::_INI
( const double*p )
{
  // Size and set DAG evaluation arrays
  if( _nVAR < _nx+_np+1+_nq ){
    _nVAR = _nx+_np+1+_nq;
    delete[] _pVAR; _pVAR = new FFVar[_nVAR];
    delete[] _dVAR; _dVAR = new realtype[_nVAR];
  }
  for( unsigned ix=0; ix<_nx; ix++ ) _pVAR[ix] = _pX[ix];
  for( unsigned ip=0; ip<_np; ip++ ) _pVAR[_nx+ip] = _pP[ip];
  _pVAR[_nx+_np] = (_pT? *_pT: 0. );
  for( unsigned iq=0; iq<_nq; iq++ ) _pVAR[_nx+_np+1+iq] = _pQ?_pQ[iq]:0.;
  for( unsigned ip=0; ip<_np; ip++ ) _dVAR[_nx+ip] = p[ip];

  // Size SUNDIALS arrays
  if( _Nx )  N_VDestroy_Serial( _Nx );
  if( _Nxq ) N_VDestroy_Serial( _Nxq );
  _Nx  = N_VNew_Serial( _nx );
  _Nxq = N_VNew_Serial( _nq );
}

inline bool
ODESLVS_SUNDIALS::_IC
()
{
  if( !_vIC.size() || _nx0 != _nx ) return false;

  _pIC = _vIC.at(0);
  _opIC = _pDAG->subgraph( _nx, _pIC );
  _dVAR[_nx+_np] = _t; // set initial time
  realtype* _Nx_dat = NV_DATA_S( _Nx ); // pointer to states
  _pDAG->eval( _opIC, _nx, _pIC, _Nx_dat, _np+1, _pVAR+_nx, _dVAR+_nx );
  realtype* _Nxq_dat = NV_DATA_S( _Nxq ); // pointer to state quadratures
  for( unsigned iq=0; iq<_nq; iq++ ) _Nxq_dat[iq] = 0.; // initial quadratures at 0

  // Reset ODE solver
  return _CVode_init();
}

inline bool
ODESLVS_SUNDIALS::_CC
()
{
  if( !_istg || _vIC.size() == 1 ) return true;
  if( _nx0 != _nx ) return false;

  _pIC = _vIC.at(_istg);
  _opIC = _pDAG->subgraph( _nx, _pIC );
  _dVAR[_nx+_np] = _t; // set current time
  realtype* _Nx_dat = NV_DATA_S( _Nx ); // pointer to state
  for( unsigned i=0; i<_nx; i++ ) _dVAR[i] = _Nx_dat[i]; // set current state values
  _pDAG->eval( _opIC, _nx, _pIC, _Nx_dat, _nVAR-_nq, _pVAR, _dVAR );
  // quadratures are continuous (_Nxq unchanged)

  // Reset ODE solver - needed to account for discontinuity
  return _CVode_init();
}

inline bool
ODESLVS_SUNDIALS::_SET
( const unsigned iRHS )
{
  if( _vRHS.size() <= iRHS ) return false;
  if( _nq && _vQUAD.size() <= iRHS ) return false;

  _pRHS = _vRHS.at( iRHS );
  _opRHS = _pDAG->subgraph( _nx, _pRHS );
  unsigned nopRHS = _opRHS.size();
  
  _pRHSQ = _nq? _vQUAD.at( iRHS ): 0;
  _opRHSQ = _pDAG->subgraph( _nq, _pRHSQ );
  if( nopRHS < _opRHSQ.size() ) nopRHS = _opRHSQ.size();
  delete[] _dRHS; _dRHS = new realtype[ nopRHS ];

  delete[] std::get<1>(_pJAC); delete[] std::get<2>(_pJAC); delete[] std::get<3>(_pJAC);
  _pJAC = _pDAG->SFAD( _nx, _pRHS, _nx, _pX ); // Jacobian in sparse format
  _opJAC = _pDAG->subgraph( std::get<0>(_pJAC), std::get<3>(_pJAC) );
  delete[] _dJAC; _dJAC  = new double[ std::get<0>(_pJAC)+_opJAC.size() ];

  return true;
}

inline bool
ODESLVS_SUNDIALS::_FCT
( const unsigned iFCT, double*f )
{
  if( !_nf || !f ) return true;

  _dVAR[_nx+_np] = _t; // set current time
  realtype* _Nx_dat = NV_DATA_S( _Nx ); // pointer to state
  for( unsigned i=0; i<_nx; i++ ) _dVAR[i] = _Nx_dat[i]; // set current state values
  if( _nq ){
    _flag = CVodeGetQuad( _cvode_mem, &_t, _Nxq );
    if( _check_flag(&_flag, "CVodeGetQuad", 1) ) return FAILURE;
    realtype* _Nxq_dat = NV_DATA_S( _Nxq ); // pointer to state
    for( unsigned i=0; i<_nq; i++ ) _dVAR[i+_np+1+i] = _Nxq_dat[i]; // set current state values
  }
  const FFVar* pFCT = _vFCT.at( iFCT );
  _pDAG->eval( _nf, pFCT, f, _nVAR, _pVAR, _dVAR, iFCT?true:false );

  return true;
}

inline typename ODESLVS_SUNDIALS::STATUS
ODESLVS_SUNDIALS::states
( const unsigned int ns, const double*tk, const double*p, double**xk,
  double*q, double*f, std::ostream&os )
{
  try{
    // Initialize trajectory integration with CVODES-SUNDIALS
    _INI( p );

    // Initialize state integration
    _t = tk[0];
    if( !_IC() ) return FATAL;
    realtype* _Nx_dat = NV_DATA_S( _Nx ); // pointer to states
    realtype* _Nxq_dat = NV_DATA_S( _Nxq ); // pointer to quadratures
    if( options.DISPLAY >= 1 ){
      _print_interm( _t, _nx, _Nx_dat, "x", os );
      _print_interm( _nq, _Nxq_dat, "q", os );
    }
    for( unsigned int ix=0; ix<_nx; ix++ ) xk[0][ix] = _Nx_dat[ix];

    // Integrate ODEs through each stage using CVODES
    pODESLVS = this;
    for( _istg=0; _istg<ns; _istg++ ){

      // Reinitialize states at intermediate times
      if( !_CC() ) return FATAL; 

      // Update list of operations in RHS and JAC
      const unsigned iRHS = ( _vRHS.size()<2? 0: _istg );
      if( (!_istg || iRHS) && !_SET( iRHS ) ) return FATAL;

      // Integrate until end of time stage
      _flag = CVodeSetStopTime( _cvode_mem, tk[_istg+1] );
      if( _check_flag(&_flag, "CVodeSetStopTime", 1) ) return FATAL;
      _flag = CVode( _cvode_mem, tk[_istg+1], _Nx, &_t, CV_NORMAL );
      if ( _check_flag(&_flag, "CVode", 1) ) return FAILURE;

      // Display and store states at intermediate time
      _Nx_dat = NV_DATA_S( _Nx ); // pointer to state
      _Nxq_dat = NV_DATA_S( _Nxq ); // pointer to quadratures
      if( options.DISPLAY >= 1 ){
        _print_interm( _t, _nx, _Nx_dat, "x", os );
        _print_interm( _nq, _Nxq_dat, "q", os );
      }
      for( unsigned int ix=0; ix<_nx; ix++ ) xk[_istg+1][ix] = _Nx_dat[ix];

      // Add intermediate function terms
      const unsigned iFCT = ( _vFCT.size()>=ns? _istg:0 );
      if( (_vFCT.size()>=ns || _istg==ns-1) && !_FCT( iFCT, f ) ) return FATAL;
    }

    // Final quadrature and function values
    if( _nq ){
      _flag = CVodeGetQuad( _cvode_mem, &_t, _Nxq );
      if( _check_flag(&_flag, "CVodeGetQuad", 1) ) return FAILURE;
    }
    if( options.DISPLAY >= 1 ) _print_interm( _nf, f, "f", os );
  }
  catch(...){
    if( options.DISPLAY >= 1 ) _print_final_stats( _cvode_mem, false, os );
    return FAILURE;
  }

  // Display some final statistics
  if( options.DISPLAY >= 1 ) _print_final_stats( _cvode_mem, false, os );
  return NORMAL;
}

inline int
ODESLVS_SUNDIALS::_RHS
( realtype t, N_Vector Nx, N_Vector Nxdot, void *user_data )
{
  _dVAR[_nx+_np] = t; // set current time
  realtype* Nx_dat = NV_DATA_S( Nx ); // pointer to state
  for( unsigned i=0; i<_nx; i++ ) _dVAR[i] = Nx_dat[i]; // set current state values
  realtype* Nxdot_dat = NV_DATA_S( Nxdot ); // pointer to state derivatives
  _pDAG->eval( _opRHS, _dRHS, _nx, _pRHS, Nxdot_dat, _nVAR-_nq, _pVAR, _dVAR );
  return 0;  
}

inline int
ODESLVS_SUNDIALS::MC_CVRHS__
( realtype t, N_Vector Nx, N_Vector Nxdot, void *user_data )
{
  ODESLVS_SUNDIALS *pODESLVS = ODESLVS_SUNDIALS::pODESLVS;
  int flag = pODESLVS->_RHS( t, Nx, Nxdot, user_data );
  ODESLVS_SUNDIALS::pODESLVS = pODESLVS;
  return flag;
}

inline int
ODESLVS_SUNDIALS::_JAC
( long int N, realtype t, N_Vector Nx, N_Vector Nxdot, DlsMat Jac, void *user_data )
{
  _dVAR[_nx+_np] = t; // set current time
  realtype* Nx_dat = NV_DATA_S( Nx ); // pointer to state
  for( unsigned i=0; i<_nx; i++ ) _dVAR[i] = Nx_dat[i]; // set current state values
  _pDAG->eval( _opJAC, _dJAC+std::get<0>(_pJAC), std::get<0>(_pJAC),
               std::get<3>(_pJAC), _dJAC, _nVAR-_nq, _pVAR, _dVAR );
  for( unsigned ie=0; ie<std::get<0>(_pJAC); ++ie ){
    DENSE_ELEM( Jac, std::get<1>(_pJAC)[ie], std::get<2>(_pJAC)[ie] ) = _dJAC[ie];
#ifdef MC__ODESLVS_SUNDIALS_DEBUG
    std::cout << "  jac[" << std::get<1>(_pJAC)[ie] << ", "
              << std::get<2>(_pJAC)[ie] << "] = " << _dJAC[ie] << std::endl;
#endif
  }
  return 0;
}

inline int
ODESLVS_SUNDIALS::MC_CVJAC__
( long int N, realtype t, N_Vector Nx, N_Vector Nxdot, DlsMat Jac, void *user_data,
  N_Vector tmp1, N_Vector tmp2, N_Vector tmp3 )
{
  ODESLVS_SUNDIALS *pODESLVS = ODESLVS_SUNDIALS::pODESLVS;
  int flag = pODESLVS->_JAC( N, t, Nx, Nxdot, Jac, user_data );
  ODESLVS_SUNDIALS::pODESLVS = pODESLVS;
  return flag;
}

inline int
ODESLVS_SUNDIALS::_RHSQ
( realtype t, N_Vector Nx, N_Vector Nqdot, void *user_data )
{
  _dVAR[_nx+_np] = t; // set current time
  realtype* Nx_dat = NV_DATA_S( Nx ); // pointer to state
  for( unsigned i=0; i<_nx; i++ ) _dVAR[i] = Nx_dat[i]; // set current state values
  realtype* Nqdot_dat = NV_DATA_S( Nqdot ); // pointer to quadrature derivatives
  _pDAG->eval( _opRHSQ, _dRHS, _nq, _pRHSQ, Nqdot_dat, _nVAR-_nq, _pVAR, _dVAR );
  return 0;  
}

inline int
ODESLVS_SUNDIALS::MC_CVRHSQ__
( realtype t, N_Vector Nx, N_Vector Nqdot, void *user_data )
{
  ODESLVS_SUNDIALS *pODESLVS = ODESLVS_SUNDIALS::pODESLVS;
  int flag = pODESLVS->_RHSQ( t, Nx, Nqdot, user_data );
  ODESLVS_SUNDIALS::pODESLVS = pODESLVS;
  return flag;
}
/*
inline typename ODESLVS_SUNDIALS::STATUS
ODESLVS_SUNDIALS::functions
( const unsigned int ns, const double*tk, const double*p, double&f,
  double*g, std::ostream&os )
{
  // Get states at intermediate times
  double **xk  = new double*[ns+1];
  for( unsigned int is=0; is<=ns; is++ ) xk[is]  = new double[_nx];
  STATUS flag = states( ns, tk, p, xk, os );
  if( flag != NORMAL ){
    for( unsigned int is=0; is<=ns; is++ ) delete[] xk[is];
    delete[] xk;
    return flag;
  }

  // Compute objective and constraint function values
  f = _pb.OBJ( p, xk, ns ).first;
  for( unsigned int ic=0; ic<_nc; ic++ )
    g[ic] = _pb.CTR( ic, p, xk, ns ).first;

  for( unsigned int is=0; is<=ns; is++ ) delete[] xk[is];
  delete[] xk;
  return NORMAL;
}  

inline typename ODESLVS_SUNDIALS::STATUS
ODESLVS_SUNDIALS::feasibility_test
( const unsigned int ns, const double*tk, const double*p, double&f,
  double&maxinfeas, std::ostream&os )
{
  // Get states at intermediate times
  double **xk  = new double*[ns+1];
  for( unsigned int is=0; is<=ns; is++ ) xk[is]  = new double[_nx];
  STATUS flag = states( ns, tk, p, xk, os );
  if( flag != NORMAL ){
    for( unsigned int is=0; is<=ns; is++ ) delete[] xk[is];
    delete[] xk;
    return flag;
  }

  // Test feasibility
  std::pair<double,DOSTRUCT::t_CTR> g;
  maxinfeas = 0.;
  for( unsigned int ic=0; ic<_nc; ic++ ){
    g = _pb.CTR( ic, p, xk, ns );
    switch( g.second ){
    case DOSTRUCT::EQ:
      maxinfeas = std::max(maxinfeas,std::fabs(g.first)); break;
    case DOSTRUCT::LE:
      maxinfeas = std::max(maxinfeas,g.first); break;
    case DOSTRUCT::GE:
      maxinfeas = std::max(maxinfeas,-g.first); break;
    }
  }

  // Compute objective function value
  f = _pb.OBJ( p, xk, ns ).first;

  // Clean-up
  for( unsigned int is=0; is<=ns; is++ ) delete[] xk[is];
  delete[] xk;
  return NORMAL;
}  

inline void
ODESLVS_SUNDIALS::_states_FSA_init
( const double*p )
{
  // Initial conditions and parameters
  _p = p;
  for( unsigned int ip=0; ip<_np; ip++ ){
    _Fp[ip] = _p[ip];
    _Fp[ip].diff(ip,_np);
  }
  for( unsigned int ix=0; ix<_nx; ix++ ){
    t_F IC = _pb.IC( ix, _Fp );
    NV_Ith_S( _Nx, ix ) = IC.x();
    for( unsigned int jp=0; jp<_np; jp++ )
      NV_Ith_S( _Nxp[jp], ix ) = IC.d(jp);
  }
}

inline typename ODESLVS_SUNDIALS::STATUS
ODESLVS_SUNDIALS::_CVode_FSA_init()
{
  // Activate forward sensitivity computations and allocate internal memory
  // related to sensitivity calculations
  static bool first_call = true;
  if( first_call ){
    _flag = CVodeSensInit1( _cvode_mem, _np, options.FSACORR, MC_CVFSA__, _Nxp );
    if( _check_flag(&_flag, "CVodeSensInit", 1) ) return FATAL;
    first_call = false;
  }
  else{
    _flag = CVodeSensReInit( _cvode_mem, options.FSACORR, _Nxp );
    if( _check_flag(&_flag, "CVodeSensReInit", 1) ) return FATAL;
  }

  // Specify the relative and absolute tolerances for sensitivity variables
  if( options.AUTOTOLS ){
    _flag = CVodeSensEEtolerances( _cvode_mem );
    if( _check_flag(&_flag, "CVodeSensEEtolerances", 1)) return FATAL;
  }
  else{
    realtype ATOLS[_np];
    for( unsigned int ip=0; ip<_np; ip++ ) ATOLS[ip] = options.ABSTOLS;
    _flag = CVodeSensSStolerances( _cvode_mem, options.RELTOLS, ATOLS );
    if( _check_flag(&_flag, "CVodeSensSStolerances", 1)) return FATAL;
  }

  // Specify the error control strategy for sensitivity variables
  _flag = CVodeSetSensErrCon( _cvode_mem, options.FSAERR );
  if( _check_flag(&_flag, "CVodeSetSensErrCon", 1) ) return FATAL;

  // Specify problem parameter information for sensitivity calculations
  _flag = CVodeSetSensParams( _cvode_mem, 0, 0, 0 );
  if( _check_flag(&_flag, "CVodeSetSensParams", 1) ) return FATAL;
  
  return NORMAL;
}

inline typename ODESLVS_SUNDIALS::STATUS
ODESLVS_SUNDIALS::states_FSA
( const unsigned int ns, const double*tk, const double*p, double**xk,
  double**xpk, std::ostream&os )
{
  // Initialize state integration
  _states_FSA_init( p );
  if( _CVode_init( tk[0] ) != NORMAL ) return FATAL;
  if( _CVode_FSA_init() != NORMAL ) return FATAL;
  for( unsigned int ix=0; ix<_nx; ix++ ){
    xk[0][ix] = NV_Ith_S( _Nx, ix );
    for( unsigned int ip=0; ip<_np; ip++ )
      xpk[0][ix+ip*_nx] = NV_Ith_S( _Nxp[ip], ix );
  }

  // Integrate ODEs through each stage using CVODES
  pODESLVS = this;
  for( _istg=0; _istg<ns; _istg++ ){
    _flag = CVodeSetStopTime( _cvode_mem, tk[_istg+1] );
    if( _check_flag(&_flag, "CVodeSetStopTime", 1) ) return FATAL;
    //_flag = CVodeSetInitStep( _cvode_mem, 1e-9 );
    //if( _check_flag(&_flag, "CVodeSetInitStep", 1) ) return FATAL;
    //_flag = CVodeSetMinStep( _cvode_mem, 1e-9 );
    //if( _check_flag(&_flag, "CVodeSetMinStep", 1) ) return FATAL;
    _flag = CVode( _cvode_mem, tk[_istg+1], _Nx, &_t, CV_NORMAL );
    if ( _check_flag(&_flag, "CVode", 1) ) return FAILURE;
    // Retreive state sensitivities at intermediate time
    _flag = CVodeGetSens( _cvode_mem, &_t, _Nxp );
    if ( _check_flag( &_flag, "CVodeGetSens", 1) ) return FAILURE;
    // Store states and sensitivities at intermediate time in xk and xpk
    for( unsigned int ix=0; ix<_nx; ix++ ){
      xk[_istg+1][ix] = NV_Ith_S( _Nx, ix );
      for( unsigned int ip=0; ip<_np; ip++ )
        xpk[_istg+1][ix+ip*_nx] = NV_Ith_S( _Nxp[ip], ix );
    }
  }

  // Display some final statistics
  if( options.DISPLAY >= 1 ) _print_final_stats( _cvode_mem, true, os );

  return NORMAL;
}

inline typename ODESLVS_SUNDIALS::STATUS
ODESLVS_SUNDIALS::functions_FSA
( const unsigned int ns, const double*tk, const double*p, double&f,
  double*fp, double*g, double**gp, std::ostream&os )
{
  // Get state sensitivities at intermediate times
  double **xk  = new double*[ns+1];
  double **xkp = new double*[ns+1];
  for( unsigned int is=0; is<=ns; is++ ){
    xk[is]  = new double[_nx];
    xkp[is] = new double[_nx*_np];
  }  
  STATUS flag = states_FSA( ns, tk, p, xk, xkp, os );
  if( flag != NORMAL ){
    for( unsigned int is=0; is<=ns; is++ ){ delete[] xk[is]; delete[] xkp[is]; }
    delete[] xk; delete[] xkp;
    return flag;
  }
  
  // Compute objective and constraint function derivatives
  t_F **Fxk  = new t_F*[ns+1];
  for( unsigned int is=0; is<=ns; is++ ) Fxk[is] = new t_F[_nx];

  for( unsigned int jp=0; jp<_np; jp++ ){
    for( unsigned int is=0; is<=ns; is++ ){
      for( unsigned int ix=0; ix<_nx; ix++ ){
        Fxk[is][ix] = xk[is][ix];
        Fxk[is][ix].diff(0,1) = xkp[is][ix+jp*_nx];
      }
    }
    for( unsigned int ip=0; ip<_np; ip++ ){
      _Fp[ip] = _p[ip];
      _Fp[ip].diff(0,1) = ( ip==jp? 1: 0 );
    }

    t_F OBJ = _pb.OBJ( _Fp, Fxk, ns ).first;
    f = OBJ.x();
    fp[jp] = OBJ.d(0);

    for( unsigned int ic=0; ic<_nc; ic++ ){
      t_F CTR = _pb.CTR( ic, _Fp, Fxk, ns ).first;
      g[ic] = CTR.x();
      gp[ic][jp] = CTR.d(0);
    }
  }

  for( unsigned int is=0; is<=ns; is++ )
    { delete[] xk[is]; delete[] xkp[is]; delete[] Fxk[is]; }
  delete[] xk; delete[] xkp; delete[] Fxk;
  return NORMAL;
}

inline typename ODESLVS_SUNDIALS::STATUS
ODESLVS_SUNDIALS::_CVode_initB()
{
  // Activate adjoint sensitivity computations and allocate internal memory
  // related to the forward-backward problem
  static bool first_call = true;
  if( first_call ){
    _flag = CVodeAdjInit( _cvode_mem, options.ASACHKPT, options.ASAINTERP );
    if( _check_flag(&_flag, "CVodeAdjInit", 1) ) return FATAL;
    first_call = false;
  }

  for( unsigned int ic=0; ic<=_nc; ic++ ){
    N_VConst( 0., _Nl[ic] );
    N_VConst( 0., _Nq[ic] );
  }

  return NORMAL;
}

inline void
ODESLVS_SUNDIALS::_states_ASA_init
( const unsigned int ns, const double* const*xk )
{
  // Transition adjoint conditions
  for( unsigned int ip=0; ip<_np; ip++ )
    _Fp[ip] = _p[ip];

  t_F **Fxk = new t_F*[ns+1];
  for( unsigned int is=0; is<=ns; is++ ){
    Fxk[is] = new t_F[_nx];
    for( unsigned int ix=0; ix<_nx; ix++ ){
      Fxk[is][ix] = xk[is][ix];
      if( is == _istg+1 ) Fxk[is][ix].diff(ix,_nx);
    }
  }

  t_F OBJ = _pb.OBJ( _Fp, Fxk, ns ).first;
  for( unsigned int ix=0; ix<_nx; ix++ )
    NV_Ith_S( _Nl[0], ix ) += OBJ.d(ix);

  for( unsigned int ic=0; ic<_nc; ic++ ){
    t_F CTR = _pb.CTR( ic, _Fp, Fxk, ns ).first;
    for( unsigned int ix=0; ix<_nx; ix++ )
      NV_Ith_S( _Nl[1+ic], ix ) += CTR.d(ix);
  }

  for( unsigned int is=0; is<=ns; is++ ) delete[] Fxk[is];
  delete[] Fxk;
}

inline typename ODESLVS_SUNDIALS::STATUS
ODESLVS_SUNDIALS::_CVode_ASA_init
( const double tf )
{
  // Create and allocate CVODES memory for backward run -> MOVE TO CONSTRUCTOR???
  static bool first_call = true;
  if( first_call ){
    for( unsigned int ic=0; ic<1+_nc; ic++ ){
      _flag = CVodeCreateB( _cvode_mem, CV_BDF, CV_NEWTON, &_indexB[ic] );
      if( _check_flag(&_flag, "CVodeCreateB", 1) ) return FATAL;
    }
  }

  for( unsigned int ic=0; ic<1+_nc; ic++ ){
    if( first_call ){
      // Initialize backward problem ic
      _flag = CVodeInitB( _cvode_mem, _indexB[ic], MC_CVRHSB__, tf, _Nl[ic] );
      if( _check_flag(&_flag, "CVodeInitB", 1)) return FATAL;

      // Specify the CVDENSE dense linear solver for backward problem ic
      _flag = CVDenseB( _cvode_mem, _indexB[ic], _nx );
      if( _check_flag(&_flag, "CVDenseB", 1)) return FATAL;

      // Set the Jacobian routine for backward problem ic
      _flag = CVDlsSetDenseJacFnB( _cvode_mem, _indexB[ic], MC_CVJACB__ );
      if( _check_flag(&_flag, "CVDlsSetDenseJacFnB", 1) ) return FATAL;

      // Initialize backward quadrature problem ic
      _flag = CVodeQuadInitB( _cvode_mem, _indexB[ic], MC_CVRHSBQ__, _Nq[ic] );
      if( _check_flag(&_flag, "CVodeQuadInitB", 1) ) return FATAL;
    }

    else{
      // Reinitialize backward problem ic
      _flag = CVodeReInitB( _cvode_mem, _indexB[ic], tf, _Nl[ic] );
      if( _check_flag(&_flag, "CVodeReInitB", 1) ) return FATAL;

      // Reinitialize backward quadrature problem ic
      _flag = CVodeQuadReInitB( _cvode_mem, _indexB[ic], _Nq[ic] );
      if( _check_flag(&_flag, "CVodeQuadReInitB", 1)) return FATAL;
    }

    // Toggle off forward sensitivities
    _flag = CVodeSensToggleOff( _cvode_mem );
    if( _check_flag(&_flag, "CVodeSensToggleOff", 1) ) return FATAL;

    // No checkpointing data for forward sensitivities
    _flag = CVodeSetAdjNoSensi( _cvode_mem );
    if( _check_flag(&_flag, "CVodeSetAdjNoSensi", 1) ) return FATAL;

    // Specify the relative and absolute tolerances for backward problem ic
    _flag = CVodeSStolerancesB( _cvode_mem, _indexB[ic], options.RELTOLB,
      options.ABSTOLB );
    if( _check_flag(&_flag, "CVodeSStolerancesB", 1) ) return FATAL;

    // Specify the relative and absolute tolerances for backward quadrature ic
    _flag = CVodeQuadSStolerancesB( _cvode_mem, _indexB[ic], options.RELTOLB,
      options.ABSTOLB );
    if( _check_flag(&_flag, "CVodeQuadSStolerancesB", 1) ) return FATAL;

    // Specify whether or not to perform error control on backward quadrature ic
    _flag = CVodeSetQuadErrConB( _cvode_mem, _indexB[ic], options.ASAQERR );
    if( _check_flag(&_flag, "CVodeSetQuadErrConB", 1) ) return FATAL;
  }

  first_call = false;
  return NORMAL;
}

inline typename ODESLVS_SUNDIALS::STATUS
ODESLVS_SUNDIALS::states_ASA
( const unsigned int ns, const double*tk, const double*p, double**xk,
  double**lk, double*q, std::ostream&os )
{
  // Initialize forward-backward integration
  _states_init( p );
  if( _CVode_init( tk[0] ) != NORMAL || _CVode_initB() != NORMAL) return FATAL;
  for( unsigned int ix=0; ix<_nx; ix++ )
    xk[0][ix] = NV_Ith_S( _Nx, ix );

  // Forward ODE integration through each stage using CVODES
  int nchk;
  pODESLVS = this;
  for( _istg=0; _istg<ns; _istg++ ){
    _flag = CVodeSetStopTime( _cvode_mem, tk[_istg+1] );
    if( _check_flag(&_flag, "CVodeSetStopTime", 1) ) return FATAL;
    _flag = CVodeF( _cvode_mem, tk[_istg+1], _Nx, &_t, CV_NORMAL, &nchk );
    if ( _check_flag(&_flag, "CVodeF", 1) ) return FAILURE;
    // Store states at intermediate time in xk
    for( unsigned int ix=0; ix<_nx; ix++ )
      xk[_istg+1][ix] = NV_Ith_S( _Nx, ix );
  }

  // Display some final statistics
  if( options.DISPLAY >= 1 ) _print_final_stats( _cvode_mem, false, os );

  // Bacward ODE integration through each stage using CVODES
  pODESLVS = this;
  _istg = ns-1;
  for( unsigned is=ns; is>0; is--, _istg-- ){
    // Update adjoints at intermediate time
    _states_ASA_init( ns, xk );
    if( _CVode_ASA_init( tk[_istg+1] ) != NORMAL ) return FATAL;
    // Store adjoints at intermediate time in lk
    for( unsigned int ic=0; ic<1+_nc; ic++ )
      for( unsigned int ix=0; ix<_nx; ix++ )
        lk[_istg+1][ix+ic*_nx] = NV_Ith_S( _Nl[ic], ix );
    // Integrate adjoint system backward
    //_flag = CVodeSetStopTimeB( _cvode_mem, tk[_istg] );
    //if( _check_flag(&_flag, "CVodeSetStopTimeB", 1) ) return FATAL;
    _flag = CVodeB( _cvode_mem, tk[_istg], CV_NORMAL );
    if ( _check_flag(&_flag, "CVodeB", 1) ) return FAILURE;
    // Retreive adjoints at intermediate time
    for( unsigned int ic=0; ic<1+_nc; ic++ ){ 
      _flag = CVodeGetB( _cvode_mem, _indexB[ic], &_t, _Nl[ic] );
      if ( _check_flag(&_flag, "CVodeGetB", 1) ) return FAILURE;
      _flag = CVodeGetQuadB( _cvode_mem, _indexB[ic], &_t, _Nq[ic] );
      if ( _check_flag(&_flag, "CVodeGetQuadB", 1) ) return FAILURE;
    }
    if( options.DISPLAY >= 2 )
      std::cout << "t(" << _istg << ")=" << _t << std::endl;
  }

  // Store adjoints and quadratures at initial time in lk and qk
  for( unsigned int ic=0; ic<1+_nc; ic++ ){
    for( unsigned int ix=0; ix<_nx; ix++ )
      lk[0][ix+ic*_nx] = NV_Ith_S( _Nl[ic], ix );
    for( unsigned int ip=0; ip<_np; ip++ ){
      q[ip+ic*_np] = NV_Ith_S( _Nq[ic], ip );
      //std::cout << "q[" << ic << "]["  << ip << "] =" << q[ip+ic*_np]  << std::endl;
    }
  }

  return NORMAL;
}

inline typename ODESLVS_SUNDIALS::STATUS
ODESLVS_SUNDIALS::functions_ASA
( const unsigned int ns, const double*tk, const double*p, double&f,
  double*fp, double*g, double**gp, std::ostream&os )
{
  // Get adjoints and quadratures at initial time
  double **xk  = new double*[ns+1];
  double **lk = new double*[ns+1];
  for( unsigned int is=0; is<=ns; is++ ){
    xk[is]  = new double[_nx];
    lk[is] = new double[_nx*(_nc+1)];
  }
  double *q = new double[_np*(_nc+1)];
  STATUS flag = states_ASA( ns, tk, p, xk, lk, q, os );
  if( flag != NORMAL ){
    for( unsigned int is=0; is<=ns; is++ ){ delete[] xk[is]; delete[] lk[is]; }
    delete[] xk; delete[] lk; delete[] q;
    return flag;
  }
  
  // Compute objective and constraint function derivatives
  double *xp0 = new double[_nx*_np];
  for( unsigned int ip=0; ip<_np; ip++ ){
    _Fp[ip] = p[ip];
    _Fp[ip].diff(ip,_np);
  }
  for( unsigned int ix=0; ix<_nx; ix++ ){
    t_F IC = _pb.IC( ix, _Fp );
    for( unsigned int ip=0; ip<_np; ip++ )
      xp0[ix+ip*_nx] = IC.d(ip);
  }

  t_F **Fxk = new t_F*[ns+1];
  for( unsigned int is=0; is<=ns; is++ )
    Fxk[is] = new t_F[_nx];
  for( unsigned int ip=0; ip<_np; ip++ ){
    for( unsigned int is=0; is<=ns; is++ ){
      for( unsigned int ix=0; ix<_nx; ix++ ){
        Fxk[is][ix] = xk[is][ix];
        if( is==0 ) Fxk[is][ix].diff(0,1) = xp0[ix+ip*_nx];
      }
    }
    for( unsigned int jp=0; jp<_np; jp++ ){
      _Fp[jp] = p[jp];
      _Fp[jp].diff(0,1) = ( jp==ip? 1: 0 );
    }

    t_F OBJ = _pb.OBJ( _Fp, Fxk, ns ).first;
    f = OBJ.x();
    fp[ip] = q[ip] + OBJ.d(0);
    for( unsigned int ix=0; ix<_nx; ix++ )
      fp[ip] += lk[0][ix]*xp0[ix+ip*_nx];  

    for( unsigned int ic=0; ic<_nc; ic++ ){
      t_F CTR = _pb.CTR( ic, _Fp, Fxk, ns ).first;
      g[ic] = CTR.x();
      gp[ic][ip] = q[ip+(ic+1)*_np] + CTR.d(0);
      for( unsigned int ix=0; ix<_nx; ix++ )
        gp[ic][ip] += lk[0][ix+(ic+1)*_nx]*xp0[ix+ip*_nx]; 
    }
  }

  for( unsigned int is=0; is<=ns; is++ )
    { delete[] xk[is]; delete[] lk[is]; delete[] Fxk[is]; }
  delete[] xk; delete[] lk; delete[] q; delete[] xp0; delete[] Fxk;
  return NORMAL;
}

inline void
ODESLVS_SUNDIALS::_states_DSOFSA_init
( const double*p )
{
  // Initial conditions for forward sensitivity problem _iparBS
  _p = p;
  for( unsigned int ip=0; ip<_np; ip++ ){
    _Fp[ip] = _p[ip];
    _Fp[ip].diff(0,1) = ( ip==_iparBS? 1: 0 );
  }
  for( unsigned int ix=0; ix<_nx; ix++ ){
    t_F IC = _pb.IC( ix, _Fp );
    NV_Ith_S( _Nx, ix ) = IC.x();
    NV_Ith_S( _Nxpi[0], ix ) = IC.d(0);
  }
}

inline typename ODESLVS_SUNDIALS::STATUS
ODESLVS_SUNDIALS::_CVode_DSOFSA_init
( const double t0 )
{
  // Activate forward sensitivity computations and allocate internal memory
  // related to sensitivity calculations
  static bool first_call = true;
  if( first_call ){
    _flag = CVodeInit( _cvode_memS, MC_CVRHS__, t0, _Nx );
    if( _check_flag(&_flag, "CVodeInit", 1) ) return FATAL;
    _flag = CVodeSensInit( _cvode_memS, 1, options.FSACORR, MC_CVDSOFSA__, _Nxpi );
    if( _check_flag(&_flag, "CVodeSensInit", 1) ) return FATAL;
  }
  else{
    _flag = CVodeReInit( _cvode_memS, t0, _Nx );
    if( _check_flag(&_flag, "CVodeReInit", 1) ) return FATAL;
    _flag = CVodeSensReInit( _cvode_memS, options.FSACORR, _Nxpi );
    if( _check_flag(&_flag, "CVodeSensReInit", 1) ) return FATAL;
  }

  // Specify the relative and absolute tolerances
  _flag = CVodeSStolerances( _cvode_memS, options.RELTOL, options.ABSTOL );
  if( _check_flag(&_flag, "CVodeSVtolerances", 1) ) return FATAL;

  // Specify the CVDENSE dense linear solver
  _flag = CVDense( _cvode_memS, _nx );
  if( _check_flag(&_flag, "CVDense", 1)) return FATAL;

  // Set the Jacobian routine
  _flag = CVDlsSetDenseJacFn( _cvode_memS, MC_CVJAC__ );
  if ( _check_flag(&_flag, "CVDlsSetDenseJacFn", 1) ) return FATAL;

  // Specify the relative and absolute tolerances for sensitivity variables
  if( options.AUTOTOLS ){
    _flag = CVodeSensEEtolerances( _cvode_memS );
    if( _check_flag(&_flag, "CVodeSensEEtolerances", 1)) return FATAL;
  }
  else{
    realtype ATOLS[_np];
    for( unsigned int ip=0; ip<_np; ip++ ) ATOLS[ip] = options.ABSTOLS;
    _flag = CVodeSensSStolerances( _cvode_memS, options.RELTOLS, ATOLS );
    if( _check_flag(&_flag, "CVodeSensSStolerances", 1)) return FATAL;
  }

  // Specify the error control strategy for sensitivity variables
  _flag = CVodeSetSensErrCon( _cvode_memS, options.FSAERR );
  if( _check_flag(&_flag, "CVodeSetSensErrCon", 1) ) return FATAL;

  // Specify problem parameter information for sensitivity calculations
  _flag = CVodeSetSensParams( _cvode_memS, 0, 0, 0 );
  if( _check_flag(&_flag, "CVodeSetSensParams", 1) ) return FATAL;

  // Set maximum number of error test failures
  _flag = CVodeSetMaxErrTestFails( _cvode_memS, options.MAXFAIL );
  if ( _check_flag(&_flag, "CVodeSetMaxErrTestFails", 1) ) return FATAL;

  first_call = false;
  return NORMAL;
}

inline typename ODESLVS_SUNDIALS::STATUS
ODESLVS_SUNDIALS::_CVode_initBS()
{
  // Activate adjoint sensitivity computations and allocate internal memory
  // related to the forward-backward sensitivity problem
  static bool first_call = true;
  if( first_call ){
    _flag = CVodeAdjInit( _cvode_memS, options.ASACHKPT, options.ASAINTERP );
    if( _check_flag(&_flag, "CVodeAdjInit", 1) ) return FATAL;
    first_call = false;
  }

  N_VConst( 0., _Nlpi );
  N_VConst( 0., _Nqpi );

  return NORMAL;
}

inline void
ODESLVS_SUNDIALS::_states_DSOASA_init
( const unsigned int ns, const double*mu, const double* const*xk,
  const double* const*xpk )
{
  // Temporary AD variables
  t_F **Fxk  = new t_F*[ns+1];
  t_B **Bxk  = new t_B*[ns+1];
  t_BF **BFxk = new t_BF*[ns+1];
  for( unsigned int is=0; is<=ns; is++ ){
    Fxk[is]  = new t_F[_nx];
    Bxk[is]  = new t_B[_nx];
    BFxk[is] = new t_BF[_nx];
  }
    
  // Transition adjoint sensitivity conditions for backward sensitivity problem _iparBS
  for( unsigned int ip=0; ip<_np; ip++ ){
    _Fp[ip] = _p[ip];
    _Fp[ip].diff(0,1) = ( ip==_iparBS? 1: 0 );
    _Bp[ip] = _p[ip];
    _BFp[ip] = _Fp[ip];
  }
  for( unsigned int is=0; is<=ns; is++ ){
    for( unsigned int ix=0; ix<_nx; ix++ ){
      Fxk[is][ix] = xk[is][ix];
      Fxk[is][ix].diff(0,1) = xpk[is][ix+_iparBS*_nx];
      Bxk[is][ix] = xk[is][ix];
      BFxk[is][ix] = Fxk[is][ix];
    }
  }

  t_B LAGR   = _pb.OBJ( _Bp,  Bxk,  ns ).first;
  t_BF DLAGR = _pb.OBJ( _BFp, BFxk, ns ).first;
  for( unsigned int ic=0; ic<_nc; ic++ ){
    LAGR  += mu[ic] * _pb.CTR( ic, _Bp,  Bxk,  ns ).first;
    DLAGR += mu[ic] * _pb.CTR( ic, _BFp, BFxk, ns ).first;
  }
  LAGR.diff(0,1);
  DLAGR.diff(0,1);
  for( unsigned int ix=0; ix<_nx; ix++ ){
    NV_Ith_S( _Nlpi, ix ) += Bxk[_istg+1][ix].d(0);
    NV_Ith_S( _Nlpi, _nx+ix ) += BFxk[_istg+1][ix].d(0).d(0);
  }  

  // Clean up
  for( unsigned int is=0; is<=ns; is++ )
  { delete[] Fxk[is]; delete[] Bxk[is]; delete[] BFxk[is]; }
  delete[] Fxk; delete[] Bxk; delete[] BFxk;
}

inline typename ODESLVS_SUNDIALS::STATUS
ODESLVS_SUNDIALS::_CVode_DSOASA_init
( const double tf )
{
  static bool first_call = true;
  if( first_call ){
    // Create and allocate CVODES memory for backward run -> MOVE TO CONSTRUCTOR???
    _flag = CVodeCreateB( _cvode_memS, CV_BDF, CV_NEWTON, &_indexBS );
    if( _check_flag(&_flag, "CVodeCreateB", 1) ) return FATAL;

    // Initialize BS problem
    _flag = CVodeInitBS( _cvode_memS, _indexBS, MC_CVRHSBS__, tf, _Nlpi );
    if( _check_flag(&_flag, "CVodeInitBS", 1)) return FATAL;

    // Specify the CVDENSE dense linear solver for BS problem
    _flag = CVDenseB( _cvode_memS, _indexBS, 2*_nx );
    if( _check_flag(&_flag, "CVDenseB", 1)) return FATAL;

    // Set the Jacobian routine for BS problem
   //_flag = CVDlsSetDenseJacFnB( _cvode_memS, _indexBS, MC_CVJACB__ );
    //if( _check_flag(&_flag, "CVDlsSetDenseJacFnB", 1) ) return FATAL;

    // Initialize backward quadrature BS problem
    _flag = CVodeQuadInitBS( _cvode_memS, _indexBS, MC_CVRHSBSQ__, _Nqpi );
    if( _check_flag(&_flag, "CVodeQuadInitBS", 1) ) return FATAL;
  }

  else{
    // Reinitialize backward BS problem
    _flag = CVodeReInitB( _cvode_memS, _indexBS, tf, _Nlpi );
    if( _check_flag(&_flag, "CVodeReInitBS", 1) ) return FATAL;

    // Reinitialize backward quadrature BS problem
    _flag = CVodeQuadReInitB( _cvode_memS, _indexBS, _Nqpi );
    if( _check_flag(&_flag, "CVodeQuadReInitBS", 1)) return FATAL;
  }

  // Specify the relative and absolute tolerances for BS problem
  _flag = CVodeSStolerancesB( _cvode_memS, _indexBS, options.RELTOLB,
    options.ABSTOLB );
  if( _check_flag(&_flag, "CVodeSStolerancesB", 1) ) return FATAL;

  // Specify the relative and absolute tolerances for BS quadrature
  _flag = CVodeQuadSStolerancesB( _cvode_memS, _indexBS, options.RELTOLB,
    options.ABSTOLB );
  if( _check_flag(&_flag, "CVodeQuadSStolerancesB", 1) ) return FATAL;

  // Specify whether or not to perform error control on BS quadrature ip
  _flag = CVodeSetQuadErrConB( _cvode_memS, _indexBS, options.ASAQERR );
  if( _check_flag(&_flag, "CVodeSetQuadErrConB", 1) ) return FATAL;
  

  first_call = false;
  return NORMAL;
}

inline typename ODESLVS_SUNDIALS::STATUS
ODESLVS_SUNDIALS::states_DSOASA
( const unsigned int ns, const double*tk, const double*p, const double*mu,
  double**xk, double**xpk, double**lpk, double*qp, std::ostream&os )
{
  for( _iparBS=0; _iparBS<_np; _iparBS++ ){

    // Initialize state integration
    _states_DSOFSA_init( p );
    if( _CVode_initBS() != NORMAL
     || _CVode_DSOFSA_init( tk[0] ) != NORMAL ) return FATAL;
    for( unsigned int ix=0; ix<_nx; ix++ ){
      xk[0][ix]  = NV_Ith_S( _Nx, ix );
      xpk[0][ix+_iparBS*_nx] = NV_Ith_S( _Nxpi[0], ix );
    }

    // Integrate ODEs through each stage using CVODES
    int nchk;
    pODESLVS = this;
    for( _istg=0; _istg<ns; _istg++ ){
      _flag = CVodeSetStopTime( _cvode_memS, tk[_istg+1] );
      if( _check_flag(&_flag, "CVodeSetStopTime", 1) ) return FATAL;
      _flag = CVodeF( _cvode_memS, tk[_istg+1], _Nx, &_t, CV_NORMAL, &nchk );
      if ( _check_flag(&_flag, "CVodeF", 1) ) return FAILURE;
      // Retreive state sensitivities at intermediate time
      _flag = CVodeGetSens( _cvode_memS, &_t, _Nxpi );
      if ( _check_flag( &_flag, "CVodeGetSens", 1) ) return FAILURE;
      // Store states and sensitivities at intermediate time in xk and xpk
      for( unsigned int ix=0; ix<_nx; ix++ ){
        xk[_istg+1][ix]  = NV_Ith_S( _Nx, ix );
        xpk[_istg+1][ix+_iparBS*_nx] = NV_Ith_S( _Nxpi[0], ix );
      }
    }

    // Display some final statistics
    if( options.DISPLAY >= 1 ) _print_final_stats( _cvode_memS, false, os );

    // Bacward ODE integration through each stage using CVODES
    pODESLVS = this;
    _istg = ns-1;
    for( unsigned is=ns; is>0; is--, _istg-- ){
      // Update adjoint sensitivities at intermediate time
      _states_DSOASA_init( ns, mu, xk, xpk );
      if( _CVode_DSOASA_init( tk[_istg+1] ) != NORMAL ) return FATAL;
      // Store adjoint sensitivities at intermediate time in lpk
      for( unsigned int ix=0; ix<2*_nx; ix++ ){
        lpk[_istg+1][ix+_iparBS*2*_nx] = NV_Ith_S( _Nlpi, ix );
      }
      // Integrate adjoint sensitivity system backward
      _flag = CVodeB( _cvode_memS, tk[_istg], CV_NORMAL );
      if ( _check_flag(&_flag, "CVodeB", 1) ) return FAILURE;
      // Retreive adjoints at intermediate time
      _flag = CVodeGetB( _cvode_memS, _indexBS, &_t, _Nlpi );
      if ( _check_flag(&_flag, "CVodeGetB", 1) ) return FAILURE;
      _flag = CVodeGetQuadB( _cvode_memS, _indexBS, &_t, _Nqpi );
      if ( _check_flag(&_flag, "CVodeGetQuadB", 1) ) return FAILURE;
      if( options.DISPLAY >= 2 )
        std::cout << "t(" << _istg << ")=" << _t << std::endl;
    }

    // Store adjoints and quadratures at initial time in lpk and qp
    for( unsigned int jp=0; jp<_np; jp++ ){
      for( unsigned int ix=0; ix<2*_nx; ix++ )
        lpk[0][ix+_iparBS*2*_nx] = NV_Ith_S( _Nlpi, ix );
      for( unsigned int ip=0; ip<_np; ip++ )
        qp[ip+_iparBS*_np] = NV_Ith_S( _Nqpi, ip );
    }
  }
  
  return NORMAL;
}

inline typename ODESLVS_SUNDIALS::STATUS
ODESLVS_SUNDIALS::functions_DSOASA
( const unsigned int ns, const double*tk, const double*p, const double*mu,
  double&f, double*fp, double*g, double**gp, double*Lpp, std::ostream&os )
{
  // Get state, adjoint and quadrature sensitivities at intermediate times
  double **xk  = new double*[ns+1];
  double **xkp = new double*[ns+1];
  double **lkp = new double*[ns+1];
  for( unsigned int is=0; is<=ns; is++ ){
    xk[is]  = new double[_nx];
    xkp[is] = new double[_nx*_np];
    lkp[is] = new double[2*_nx*_np];
  }  
  double *qp = new double[_np*_np];
  STATUS flag = states_DSOASA( ns, tk, p, mu, xk, xkp, lkp, qp, os );
  if( flag != NORMAL ){
    for( unsigned int is=0; is<=ns; is++ )
    { delete[] xk[is]; delete[] xkp[is]; delete[] lkp[is]; }
    delete[] xk; delete[] xkp; delete[] lkp; delete[] qp;
    return flag;
  }

  // Intermediate variables
  t_F **Fxk = new t_F*[ns+1];
  t_B **Bxk = new t_B*[ns+1];
  t_BF **BFxk = new t_BF*[ns+1];
  for( unsigned int is=0; is<=ns; is++ ){
    Fxk[is] = new t_F[_nx];
    Bxk[is] = new t_B[_nx];
    BFxk[is] = new t_BF[_nx];
  }

  for( unsigned int jp=0; jp<_np; jp++ ){

    // Initialize variables for forward-backward AD
    for( unsigned int ip=0; ip<_np; ip++ ){
      _Fp[ip] = p[ip];
      _Fp[ip].diff(0,1) = ( ip==jp? 1.:0. );
      _Bp[ip] = p[ip];
      _BFp[ip] = _Fp[ip];
    }
    for( unsigned int is=0; is<=ns; is++ ){
      for( unsigned int ix=0; ix<_nx; ix++ ){
        Fxk[is][ix] = xk[is][ix];
        Fxk[is][ix].diff(0,1) = xkp[is][ix+jp*_nx];
        Bxk[is][ix] = xk[is][ix];
        BFxk[is][ix] = Fxk[is][ix];
      }
    }

    // Compute objective and constraint function values and first derivatives
    t_F OBJ = _pb.OBJ( _Fp, Fxk, ns ).first;
    f = OBJ.x();
    fp[jp] = OBJ.d(0);
    for( unsigned int ic=0; ic<_nc; ic++ ){
      t_F CTR = _pb.CTR( ic, _Fp, Fxk, ns ).first;
      g[ic] = CTR.x();
      gp[ic][jp] = CTR.d(0);
    }

    // Compute Lagrangian DSO derivatives
    t_B  LAGR  = _pb.OBJ( _Bp,  Bxk,  ns ).first;
    t_BF DLAGR = _pb.OBJ( _BFp, BFxk, ns ).first;
    for( unsigned int ic=0; ic<_nc; ic++ ){
      LAGR  += mu[ic] * _pb.CTR( ic, _Bp,  Bxk,  ns ).first;
      DLAGR += mu[ic] * _pb.CTR( ic, _BFp, BFxk, ns ).first;
    }
    LAGR.diff(0,1);

    t_BF DIC = 0.;
    for( unsigned int ix=0; ix<_nx; ix++ )
      DIC += ( Bxk[0][ix].d(0) + lkp[0][ix+jp*2*_nx] )
        * _pb.IC( ix, _BFp );
    DLAGR.diff(0,2);
    DIC.diff(1,2);    

    switch( options.HESSFORMAT ){
    case Options::HESSFULL:
      for( unsigned int ip=0; ip<_np; ip++ ){
        Lpp[ip+jp*_np] = _BFp[ip].d(0).d(0) + _BFp[ip].d(1).d(0) + qp[ip+jp*_np];
        for( unsigned int ix=0; ix<_nx; ix++ )
          Lpp[ip+jp*_np] += ( BFxk[0][ix].d(0).d(0) + lkp[0][_nx+ix+jp*2*_nx] )
	    * xkp[0][ix+ip*_nx];
      }
      break;
    case Options::HESSLOWER:
      for( unsigned int ip=jp; ip<_np; ip++ ){
        unsigned int ie = jp*_np-(jp*(jp-1))/2+ip-jp;
        Lpp[ie] = ( _BFp[ip].d(0).d(0) + _BFp[ip].d(1).d(0)
	  + qp[ip+jp*_np] + _BFp[jp].d(0).d(0) + _BFp[jp].d(1).d(0)
	  + qp[jp+ip*_np] ) / 2.;
        for( unsigned int ix=0; ix<_nx; ix++ )
          Lpp[ie] += (( BFxk[0][ix].d(0).d(0) + lkp[0][_nx+ix+jp*2*_nx] )
	    * xkp[0][ix+ip*_nx] + ( BFxk[0][ix].d(0).d(0) + lkp[0][_nx+ix+ip*2*_nx] )
	    * xkp[0][ix+jp*_nx] ) / 2.;
      }
      break;
    case Options::HESSUPPER:
      for( unsigned int ip=0; ip<=jp; ip++ ){
        unsigned int ie = (jp*(jp+1))/2+ip;
        Lpp[ie] = ( _BFp[ip].d(0).d(0) + _BFp[ip].d(1).d(0)
	  + qp[ip+jp*_np] + _BFp[jp].d(0).d(0) + _BFp[jp].d(1).d(0)
	  + qp[jp+ip*_np] ) / 2.;
        for( unsigned int ix=0; ix<_nx; ix++ )
          Lpp[ie] += (( BFxk[0][ix].d(0).d(0) + lkp[0][_nx+ix+jp*2*_nx] )
	    * xkp[0][ix+ip*_nx] + ( BFxk[0][ix].d(0).d(0) + lkp[0][_nx+ix+ip*2*_nx] )
	    * xkp[0][ix+jp*_nx] ) / 2.;
      }
      break;
    }
  }

  // Clean-up
  for( unsigned int is=0; is<=ns; is++ )
    { delete[] xk[is];  delete[] xkp[is]; delete[] lkp[is];
      delete[] Fxk[is]; delete[] Bxk[is]; delete[] BFxk[is]; }
  delete[] xk;  delete[] xkp; delete[] lkp; delete[] qp;
  delete[] Fxk; delete[] Bxk; delete[] BFxk; 
  return NORMAL;
}

inline typename ODESLVS_SUNDIALS::STATUS
ODESLVS_SUNDIALS::states_FDSA
( const unsigned int ns, const double*tk, const double*p, double**xk,
  double**xpk, std::ostream&os )
{
  STATUS flag = NORMAL;

  // Get nominal states at intermediate times
  flag = states( ns, tk, p, xk, os );
  if( flag != NORMAL ) return flag;
 
   // Get perturbed states at intermediate times
  double*pp = new double[_np];
  double **xkp = new double*[ns];
  double*pm = ( options.FDCEN? new double[_np]: const_cast<double*>(p) );
  double **xkm = ( options.FDCEN? new double*[ns]: xk );
  for( unsigned int is=0; is<ns; is++ ){
    xkp[is] = new double[_nx];
    if( options.FDCEN ) xkm[is] = new double[_nx];
  }

  for( unsigned int ip=0; ip<_np; ip++ ){
    for( unsigned int jp=0; jp<_np; jp++ )
      pp[jp] = ( ip!=jp? p[jp]: p[jp]+std::fabs(p[jp])*options.FDRTOL
        +options.FDATOL );
    flag = states( ns, tk, pp, xkp, os );
    if( flag != NORMAL ){
      for( unsigned int is=0; is<ns; is++ ) delete[] xkp[is];
      delete[] pp; delete[] xkp;
      if( options.FDCEN ){
        for( unsigned int is=0; is<ns; is++ ) delete[] xkm[is];
        delete[] pm; delete[] xkm;
      }
      return flag;
    }
    if( options.FDCEN ){
      for( unsigned int jp=0; jp<_np; jp++ )
        pm[jp] = ( ip!=jp? p[jp]: p[jp]-std::fabs(p[jp])*options.FDRTOL
          -options.FDATOL );
      flag = states( ns, tk, pm, xkm, os );
      if( flag != NORMAL ){
        for( unsigned int is=0; is<ns; is++ ) delete[] xkp[is];
        delete[] pp; delete[] xkp;
        if( options.FDCEN ){
          for( unsigned int is=0; is<ns; is++ ) delete[] xkm[is];
          delete[] pm; delete[] xkm;
        }
        return flag;
      }
    }
    for( unsigned int is=0; is<ns; is++ ){
      for( unsigned int ix=0; ix<_nx; ix++ ){
        xpk[is][ix+ip*_nx] = (xkp[is][ix]-xkm[is][ix])/(pp[ip]-pm[ip]);
      }
    }
  }

  for( unsigned int is=0; is<ns; is++ ) delete[] xkp[is];
  delete[] pp; delete[] xkp;
  if( options.FDCEN ){
    for( unsigned int is=0; is<ns; is++ ) delete[] xkm[is];
    delete[] pm; delete[] xkm;
  }
  return NORMAL;
}

inline typename ODESLVS_SUNDIALS::STATUS
ODESLVS_SUNDIALS::functions_FDSA
( const unsigned int ns, const double*tk, const double*p, double&f,
  double*dfdp, double*g, double**dgdp, std::ostream&os )
{
  STATUS flag = NORMAL;

  // Get nominal states at intermediate times
  flag = functions( ns, tk, p, f, g, os );
  if( flag != NORMAL ) return flag;
 
   // Get perturbed states at intermediate times
  double*pp = new double[_np];
  double fp, *gp = new double[_nc];
  double*pm = ( options.FDCEN? new double[_np]: const_cast<double*>(p) );
  double fm = f, *gm = ( options.FDCEN? new double[_nc]: g );

  for( unsigned int ip=0; ip<_np; ip++ ){
    for( unsigned int jp=0; jp<_np; jp++ )
      pp[jp] = ( ip!=jp? p[jp]: p[jp]+std::fabs(p[jp])*options.FDRTOL
        +options.FDATOL );
    flag = functions( ns, tk, pp, fp, gp, os );
    if( flag != NORMAL ){
      delete[] pp; delete[] gp;
      if( options.FDCEN ){ delete[] pm; delete[] gm; }
      return flag;
    }
    if( options.FDCEN ){
      for( unsigned int jp=0; jp<_np; jp++ )
        pm[jp] = ( ip!=jp? p[jp]: p[jp]-std::fabs(p[jp])*options.FDRTOL
          -options.FDATOL );
      flag = functions( ns, tk, pm, fm, gm, os );
      if( flag != NORMAL ){
        delete[] pp; delete[] gp;
        if( options.FDCEN ){ delete[] pm; delete[] gm; }
        return flag;
      }
    }
    dfdp[ip] = (fp-fm)/(pp[ip]-pm[ip]);
    for( unsigned int ic=0; ic<_nc; ic++ )
      dgdp[ic][ip] = (gp[ic]-gm[ic])/(pp[ip]-pm[ip]);
  }

  delete[] pp; delete[] gp;
  if( options.FDCEN ){ delete[] pm; delete[] gm; }
  return NORMAL;
}

inline typename ODESLVS_SUNDIALS::STATUS
ODESLVS_SUNDIALS::functions_FDSOASA
( const unsigned int ns, const double*tk, const double*p, const double*mu,
  double&f, double*dfdp, double*g, double**dgdp, double*Lpp, std::ostream&os )
{
  STATUS flag = NORMAL;

  // Get nominal states at intermediate times
  flag = functions( ns, tk, p, f, g, os );
  if( flag != NORMAL ) return flag;
 
   // Get perturbed states at intermediate times
  double*pp = new double[_np];
  double fp, *gp = new double[_nc];
  double*ppp = new double[_np], *ppm = new double[_np];
  double fpp, fpm, *gpp = new double[_nc], *gpm = new double[_nc];
  double*pm = new double[_np];
  double fm, *gm = new double[_nc];
  double*pmp = new double[_np], *pmm = new double[_np];
  double fmp, fmm, *gmp = new double[_nc], *gmm = new double[_nc];

  for( unsigned int ip=0; ip<_np; ip++ ){
    for( unsigned int kp=0; kp<_np; kp++ )
      pp[kp] = ( kp!=ip? p[kp]: p[kp]+std::fabs(p[kp])*options.FDRTOL
        +options.FDATOL );
    flag = functions( ns, tk, pp, fp, gp, os );
    if( flag != NORMAL ){
      delete[] pp; delete[] gp; delete[] ppp; delete[] ppm; delete[] gpp; delete[] gpm;
      delete[] pm; delete[] gm; delete[] pmp; delete[] pmm; delete[] gmp; delete[] gmm;
      return flag;
    }
    for( unsigned int kp=0; kp<_np; kp++ )
      pm[kp] = ( kp!=ip? p[kp]: p[kp]-std::fabs(p[kp])*options.FDRTOL
        -options.FDATOL );
    flag = functions( ns, tk, pm, fm, gm, os );
    if( flag != NORMAL ){
      delete[] pp; delete[] gp; delete[] ppp; delete[] ppm; delete[] gpp; delete[] gpm;
      delete[] pm; delete[] gm; delete[] pmp; delete[] pmm; delete[] gmp; delete[] gmm;
      return flag;
    }
    dfdp[ip] = (fp-fm)/(pp[ip]-pm[ip]);
    for( unsigned int ic=0; ic<_nc; ic++ )
      dgdp[ic][ip] = (gp[ic]-gm[ic])/(pp[ip]-pm[ip]);

    switch( options.HESSFORMAT ){
    case Options::HESSFULL:
      Lpp[ip+ip*_np] = (fp-2.*f+fm)/((pp[ip]-pm[ip])*(pp[ip]-pm[ip]))*4.;
      for( unsigned int ic=0; ic<_nc; ic++ )
        Lpp[ip+ip*_np] += mu[ic]*(gp[ic]-2.*g[ic]+gm[ic])/((pp[ip]-pm[ip])*(pp[ip]-pm[ip]))*4.;
      break;
    case Options::HESSLOWER:{
      unsigned int ie = ip*_np-(ip*(ip-1))/2;
      Lpp[ie] = (fp-2.*f+fm)/((pp[ip]-pm[ip])*(pp[ip]-pm[ip]))*4.;
      for( unsigned int ic=0; ic<_nc; ic++ )
        Lpp[ie] += mu[ic]*(gp[ic]-2.*g[ic]+gm[ic])/((pp[ip]-pm[ip])*(pp[ip]-pm[ip]))*4.;
      }
      break;
    case Options::HESSUPPER:{
      unsigned int ie = (ip*(ip+1))/2+ip;
      Lpp[ie] = (fp-2.*f+fm)/((pp[ip]-pm[ip])*(pp[ip]-pm[ip]))*4.;
      for( unsigned int ic=0; ic<_nc; ic++ )
        Lpp[ie] += mu[ic]*(gp[ic]-2.*g[ic]+gm[ic])/((pp[ip]-pm[ip])*(pp[ip]-pm[ip]))*4.;
      }
      break;
    }
    
    for( unsigned int jp=0; jp<ip; jp++ ){
      for( unsigned int kp=0; kp<_np; kp++ )
        ppp[kp] = ( kp!=jp? pp[kp]: pp[kp]+std::fabs(pp[kp])*options.FDRTOL
          +options.FDATOL );
      flag = functions( ns, tk, ppp, fpp, gpp, os );
      if( flag != NORMAL ){
        delete[] pp; delete[] gp; delete[] ppp; delete[] ppm; delete[] gpp; delete[] gpm;
        delete[] pm; delete[] gm; delete[] pmp; delete[] pmm; delete[] gmp; delete[] gmm;
        return flag;
      }      
      for( unsigned int kp=0; kp<_np; kp++ )
        ppm[kp] = ( kp!=jp? pp[kp]: pp[kp]-std::fabs(pp[kp])*options.FDRTOL
          -options.FDATOL );
      flag = functions( ns, tk, ppm, fpm, gpm, os );
      if( flag != NORMAL ){
        delete[] pp; delete[] gp; delete[] ppp; delete[] ppm; delete[] gpp; delete[] gpm;
        delete[] pm; delete[] gm; delete[] pmp; delete[] pmm; delete[] gmp; delete[] gmm;
        return flag;
      }
      for( unsigned int kp=0; kp<_np; kp++ )
        pmp[kp] = ( kp!=jp? pm[kp]: pm[kp]+std::fabs(pm[kp])*options.FDRTOL
          +options.FDATOL );
      flag = functions( ns, tk, pmp, fmp, gmp, os );
      if( flag != NORMAL ){
        delete[] pp; delete[] gp; delete[] ppp; delete[] ppm; delete[] gpp; delete[] gpm;
        delete[] pm; delete[] gm; delete[] pmp; delete[] pmm; delete[] gmp; delete[] gmm;
        return flag;
      }      
      for( unsigned int kp=0; kp<_np; kp++ )
        pmm[kp] = ( kp!=jp? pm[kp]: pm[kp]-std::fabs(pm[kp])*options.FDRTOL
          -options.FDATOL );
      flag = functions( ns, tk, pmm, fmm, gmm, os );
      if( flag != NORMAL ){
        delete[] pp; delete[] gp; delete[] ppp; delete[] ppm; delete[] gpp; delete[] gpm;
        delete[] pm; delete[] gm; delete[] pmp; delete[] pmm; delete[] gmp; delete[] gmm;
        return flag;
      }

      switch( options.HESSFORMAT ){
      case Options::HESSFULL:
        Lpp[ip+jp*_np] = (fpp-fpm-fmp+fmm)/((ppp[ip]-pmp[ip])*(ppp[jp]-ppm[jp]));
        for( unsigned int ic=0; ic<_nc; ic++ )
          Lpp[ip+jp*_np] += mu[ic]*(gpp[ic]-gpm[ic]-gmp[ic]+gmm[ic])/((ppp[ip]-pmp[ip])*(ppp[jp]-ppm[jp]));
        Lpp[jp+ip*_np] = Lpp[ip+jp*_np];
        break;
      case Options::HESSLOWER:{
        unsigned int ie = jp*_np-(jp*(jp-1))/2+ip-jp;
        Lpp[ie] = (fpp-fpm-fmp+fmm)/((ppp[ip]-pmp[ip])*(ppp[jp]-ppm[jp]));
        for( unsigned int ic=0; ic<_nc; ic++ )
          Lpp[ie] += mu[ic]*(gpp[ic]-gpm[ic]-gmp[ic]+gmm[ic])/((ppp[ip]-pmp[ip])*(ppp[jp]-ppm[jp]));
        }
        break;
      case Options::HESSUPPER:{
        unsigned int ie = (jp*(jp+1))/2+ip;
        Lpp[ie] = (fpp-fpm-fmp+fmm)/((ppp[ip]-ppm[ip])*(pmp[jp]-pmm[jp]));
        for( unsigned int ic=0; ic<_nc; ic++ )
          Lpp[ie] += mu[ic]*(gpp[ic]-gpm[ic]-gmp[ic]+gmm[ic])/((ppp[ip]-pmp[ip])*(ppp[jp]-ppm[jp]));
        }
        break;
      }
    }
    
  }

  delete[] pp; delete[] gp; delete[] ppp; delete[] ppm; delete[] gpp; delete[] gpm;
  delete[] pm; delete[] gm; delete[] pmp; delete[] pmm; delete[] gmp; delete[] gmm;
  return NORMAL;

}

inline int
ODESLVS_SUNDIALS::_FSA
( int Np, realtype t, N_Vector Nx, N_Vector Nxdot, int ipcur, N_Vector Nxp,
  N_Vector Nxpdot, void *user_data )
{
  for( unsigned int ix=0; ix<_nx; ix++ ){
    _Fx[ix] = NV_Ith_S( Nx, ix );
    _Fx[ix].diff(0,1) = NV_Ith_S( Nxp, ix );
  }
  for( unsigned int ip=0; ip<_np; ip++ ){
    _Fp[ip] = _p[ip];
    _Fp[ip].diff(0,1) = ( ip==(unsigned int)ipcur? 1.: 0. );
  }

  for( unsigned int ix=0; ix<_nx; ix++ ){
    t_F RHS = _pb.RHS( ix, _Fp, _Fx, t, _istg );
    NV_Ith_S( Nxpdot, ix ) = RHS.d(0);
  }
  return 0;
}

inline int
ODESLVS_SUNDIALS::MC_CVFSA__
( int Np, realtype t, N_Vector Nx, N_Vector Nxdot, int ip, N_Vector Nxp,
  N_Vector Nxpdot, void *user_data, N_Vector tmp1, N_Vector tmp2 )
{
  ODESLVS_SUNDIALS *pODESLVS = ODESLVS_SUNDIALS::pODESLVS;
  int flag = pODESLVS_SUNDIALS->_FSA( Np, t, Nx, Nxdot, ip, Nxp, Nxpdot,
    user_data );
  ODESLVS_SUNDIALS::pODESLVS = pODESLVS;
  return flag;
}

inline int
ODESLVS_SUNDIALS::_DSOFSA
( int Np, realtype t, N_Vector Nx, N_Vector Nxdot, N_Vector *Nxp,
  N_Vector *Nxpdot, void *user_data )
{
  for( unsigned int ip=0; ip<_np; ip++ ){
    _Fp[ip] = _p[ip];
    _Fp[ip].diff(0,1) = ( ip==_iparBS? 1.: 0. );
  }
  for( unsigned int ix=0; ix<_nx; ix++ ){
    _Fx[ix] = NV_Ith_S( Nx, ix );
    _Fx[ix].diff(0,1) = NV_Ith_S( Nxp[0], ix );
  }

  for( unsigned int ix=0; ix<_nx; ix++ ){
    t_F RHSFS = _pb.RHS( ix, _Fp, _Fx, t, _istg );
    NV_Ith_S( Nxpdot[0], ix ) = RHSFS.d(0);
  }
  return 0;
}

inline int
ODESLVS_SUNDIALS::MC_CVDSOFSA__
( int Np, realtype t, N_Vector Nx, N_Vector Nxdot, N_Vector *Nxp,
  N_Vector *Nxpdot, void *user_data, N_Vector tmp1, N_Vector tmp2 )
{
  ODESLVS_SUNDIALS *pODESLVS = ODESLVS_SUNDIALS::pODESLVS;
  int flag = pODESLVS->_DSOFSA( Np, t, Nx, Nxdot, Nxp, Nxpdot,
    user_data );
  ODESLVS_SUNDIALS::pODESLVS = pODESLVS;
  return flag;
}

inline int
ODESLVS_SUNDIALS::_RHSB
( realtype t, N_Vector Nx, N_Vector NxB, N_Vector NxBdot, void *user_data )
{
  for( unsigned int ix=0; ix<_nx; ix++ ) _Bx[ix] = NV_Ith_S( Nx, ix );
  for( unsigned int ip=0; ip<_np; ip++ ) _Bp[ip] = _p[ip];

  t_B RHSB = 0.;
  for( unsigned int ix=0; ix<_nx; ix++ )
    RHSB -= NV_Ith_S( NxB, ix ) * _pb.RHS( ix, _Bp, _Bx, t, _istg );
  RHSB.diff(0,1); 
  for( unsigned int ix=0; ix<_nx; ix++ )
    NV_Ith_S( NxBdot, ix ) = _Bx[ix].d(0);

  return 0;
}

inline int
ODESLVS_SUNDIALS::MC_CVRHSB__
( realtype t, N_Vector Nx, N_Vector NxB, N_Vector NxBdot, void *user_data )
{
  ODESLVS_SUNDIALS *pODESLVS = ODESLVS_SUNDIALS::pODESLVS;
  int flag = pODESLVS->_RHSB( t, Nx, NxB, NxBdot, user_data );
  ODESLVS_SUNDIALS::pODESLVS = pODESLVS;
  return flag;
}

inline int
ODESLVS_SUNDIALS::_JACB
( long int NB, realtype t, N_Vector Nx, N_Vector NxB, N_Vector NxBdot, DlsMat JacB,
  void *user_data )
{
  for( unsigned int ix=0; ix<_nx; ix++ ){
    _Fx[ix] = NV_Ith_S( Nx, ix );
    _Fx[ix].diff(ix,_nx);
  }
  for( unsigned int ix=0; ix<_nx; ix++ ){
    t_F RHSB = -_pb.RHS( ix, _p, _Fx, t, _istg );
    for( unsigned int jx=0; jx<_nx; jx++ )
      DENSE_ELEM( JacB, jx, ix ) = RHSB.d(jx);
  }
  return 0;
}

inline int
ODESLVS_SUNDIALS::MC_CVJACB__
( long int NB, realtype t, N_Vector Nx, N_Vector NxB, N_Vector NxBdot, DlsMat JacB,
  void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3 )
{
  ODESLVS_SUNDIALS *pODESLVS = ODESLVS_SUNDIALS::pODESLVS;
  int flag = pODESLVS->_JACB( NB, t, Nx, NxB, NxBdot, JacB, user_data );
  ODESLVS_SUNDIALS::pODESLVS = pODESLVS;
  return flag;
}

inline int
ODESLVS_SUNDIALS::_RHSBQ
( realtype t, N_Vector Nx, N_Vector NxB, N_Vector NqBdot, void *user_data )
{
  for( unsigned int ix=0; ix<_nx; ix++ ) _Bx[ix] = NV_Ith_S( Nx, ix );
  for( unsigned int ip=0; ip<_np; ip++ ) _Bp[ip] = _p[ip];

  t_B RHSBQ = 0.;
  for( unsigned int ix=0; ix<_nx; ix++ )
    RHSBQ -= NV_Ith_S( NxB, ix ) * _pb.RHS( ix, _Bp, _Bx, t, _istg );
  RHSBQ.diff(0,1); 
  for( unsigned int ip=0; ip<_np; ip++ )
    NV_Ith_S( NqBdot, ip ) = _Bp[ip].d(0);

  return 0;
}

inline int
ODESLVS_SUNDIALS::MC_CVRHSBQ__
( realtype t, N_Vector Nx, N_Vector NxB, N_Vector NqBdot, void *user_data )
{
  ODESLVS_SUNDIALS *pODESLVS = ODESLVS_SUNDIALS::pODESLVS;
  int flag = pODESLVS->_RHSBQ( t, Nx, NxB, NqBdot, user_data );
  ODESLVS_SUNDIALS::pODESLVS = pODESLVS;
  return flag;
}

inline int
ODESLVS_SUNDIALS::_RHSBS
( realtype t, N_Vector Nx, N_Vector *NxS, N_Vector NxB, N_Vector NxBdot,
  void *user_data )
{
  for( unsigned int ip=0; ip<_np; ip++ ){
    _Fp[ip] = _p[ip];
    _Fp[ip].diff(0,1) = ( ip==_iparBS? 1: 0 );
    _Bp[ip] = _p[ip];
    _BFp[ip] = _Fp[ip];
  }
  for( unsigned int ix=0; ix<_nx; ix++ ){
    _Fx[ix] = NV_Ith_S( Nx, ix );
    _Fx[ix].diff(0,1) = NV_Ith_S( NxS[0], ix );
    _Bx[ix] = NV_Ith_S( Nx, ix );
    _BFx[ix] = _Fx[ix];
  }

  t_B RHSB    = 0.;
  t_B RHSBp   = 0.;
  t_BF DRHSBp = 0.;
  for( unsigned int ix=0; ix<_nx; ix++ ){
    RHSB   -= NV_Ith_S( NxB, ix ) * _pb.RHS( ix, _Bp, _Bx, t, _istg );
    RHSBp  -= NV_Ith_S( NxB, _nx+ix ) * _pb.RHS( ix, _Bp, _Bx, t, _istg );
    DRHSBp -= NV_Ith_S( NxB, ix ) * _pb.RHS( ix, _BFp, _BFx, t, _istg );
  }
  RHSB.diff(0,2); 
  RHSBp.diff(1,2); 
  DRHSBp.diff(0,1); 
  for( unsigned int ix=0; ix<_nx; ix++ ){
    NV_Ith_S( NxBdot, ix ) = _Bx[ix].d(0);
    NV_Ith_S( NxBdot, _nx+ix ) = _Bx[ix].d(1) + _BFx[ix].d(0).d(0);
  }
  return 0;
}

inline int
ODESLVS_SUNDIALS::MC_CVRHSBS__
( realtype t, N_Vector Nx, N_Vector *NxS, N_Vector NxB, N_Vector NxBdot,
  void *user_data )
{
  ODESLVS_SUNDIALS *pODESLVS = ODESLVS_SUNDIALS::pODESLVS;
  int flag = pODESLVS->_RHSBS( t, Nx, NxS, NxB, NxBdot, user_data );
  ODESLVS_SUNDIALS::pODESLVS = pODESLVS;
  return flag;
}

inline int
ODESLVS_SUNDIALS::_RHSBSQ
( realtype t, N_Vector Nx, N_Vector *NxS, N_Vector NxB, N_Vector NqBdot,
  void *user_data )
{
  for( unsigned int ip=0; ip<_np; ip++ ){
    _Fp[ip] = _p[ip];
    _Fp[ip].diff(0,1) = ( ip==_iparBS? 1: 0 );
    _Bp[ip] = _p[ip];
    _BFp[ip] = _Fp[ip];
  }
  for( unsigned int ix=0; ix<_nx; ix++ ){
    _Fx[ix] = NV_Ith_S( Nx, ix );
    _Fx[ix].diff(0,1) = NV_Ith_S( NxS[0], ix );
    _Bx[ix] = NV_Ith_S( Nx, ix );
    _BFx[ix] = _Fx[ix];
  }

  t_B RHSBQp  = 0.;
  t_BF DRHSBQp = 0.;
  for( unsigned int ix=0; ix<_nx; ix++ ){
    RHSBQp  -= NV_Ith_S( NxB, _nx+ix ) * _pb.RHS( ix, _Bp, _Bx, t, _istg );
    DRHSBQp -= NV_Ith_S( NxB, ix ) * _pb.RHS( ix, _BFp, _BFx, t, _istg );
  }
  RHSBQp.diff(0,1); 
  DRHSBQp.diff(0,1); 
  for( unsigned int ip=0; ip<_np; ip++ ){
    NV_Ith_S( NqBdot, ip ) = _Bp[ip].d(0) + _BFp[ip].d(0).d(0);
  }
  return 0;
}

inline int
ODESLVS_SUNDIALS::MC_CVRHSBSQ__
( realtype t, N_Vector Nx, N_Vector *NxS, N_Vector NxB, N_Vector NqBdot,
  void *user_data )
{
  ODESLVS_SUNDIALS *pODESLVS = ODESLVS_SUNDIALS::pODESLVS;
  int flag = pODESLVS->_RHSBSQ( t, Nx, NxS, NxB, NqBdot, user_data );
  ODESLVS_SUNDIALS::pODESLVS = pODESLVS;
  return flag;
}
*/
inline bool
ODESLVS_SUNDIALS::_check_flag
( void *flagvalue, std::string funcname, int opt )
{
  // Check if SUNDIALS function returned NULL pointer - no memory allocated
  if (opt == 0 && flagvalue == NULL) {
    std::cerr << "\nSUNDIALS_ERROR: " << funcname
              << "() failed - returned NULL pointer\n\n";
    return(true); }

  // Check if flag < 0
  else if (opt == 1) {
    int *errflag = (int *) flagvalue;
    if (*errflag < 0) {
      std::cerr << "\nSUNDIALS_ERROR: " << funcname
                << "() failed with flag = " << *errflag << "\n\n";
      return(true); }}

  // Check if function returned NULL pointer - no memory allocated
  else if (opt == 2 && flagvalue == NULL) {
    std::cerr << "\nMEMORY_ERROR: " << funcname
              << "() failed - returned NULL pointer\n\n";
    return(true); }

  return(false);
}

inline void
ODESLVS_SUNDIALS::_print_final_stats
( void *cvode_mem, bool sens, std::ostream&os )
{
  long int nst, nfe, nsetups, nje, nfeLS, nni, ncfn, netf;

  int flag = CVodeGetNumSteps(cvode_mem, &nst);
  if( _check_flag(&flag, "CVodeGetNumSteps", 1) ) return;
  flag = CVodeGetNumRhsEvals(cvode_mem, &nfe);
  if( _check_flag(&flag, "CVodeGetNumRhsEvals", 1) ) return;
  flag = CVodeGetNumLinSolvSetups(cvode_mem, &nsetups);
  if( _check_flag(&flag, "CVodeGetNumLinSolvSetups", 1) ) return;
  flag = CVodeGetNumErrTestFails(cvode_mem, &netf);
  if( _check_flag(&flag, "CVodeGetNumErrTestFails", 1) ) return;
  flag = CVodeGetNumNonlinSolvIters(cvode_mem, &nni);
  if( _check_flag(&flag, "CVodeGetNumNonlinSolvIters", 1) ) return;
  flag = CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);
  if( _check_flag(&flag, "CVodeGetNumNonlinSolvConvFails", 1) ) return;

  flag = CVDlsGetNumJacEvals(cvode_mem, &nje);
  if( _check_flag(&flag, "CVDlsGetNumJacEvals", 1) ) return;
  flag = CVDlsGetNumRhsEvals(cvode_mem, &nfeLS);
  if( _check_flag(&flag, "CVDlsGetNumRhsEvals", 1) ) return;

  os << "\nFinal Statistics:\n"
     << "   nst = " << nst << "   nfe  = " << nfe << "   nsetups = " << nsetups
     << "   nfeLS = " << nfeLS << "   nje = " << nje << std::endl
     << "   nni = " << nni << "   ncfn = " << ncfn << "   netf = " << netf
     << std::endl;

  if( sens ){
    long int nfSe, nfeS, nsetupsS, nniS, ncfnS, netfS;

    flag = CVodeGetSensNumRhsEvals(cvode_mem, &nfSe);
    if( _check_flag(&flag, "CVodeGetSensNumRhsEvals", 1) ) return;
    flag = CVodeGetNumRhsEvalsSens(cvode_mem, &nfeS);
    if( _check_flag(&flag, "CVodeGetNumRhsEvalsSens", 1) ) return;
    flag = CVodeGetSensNumLinSolvSetups(cvode_mem, &nsetupsS);
    if( _check_flag(&flag, "CVodeGetSensNumLinSolvSetups", 1) ) return;
    flag = CVodeGetSensNumErrTestFails(cvode_mem, &netfS);
    if( _check_flag(&flag, "CVodeGetSensNumErrTestFails", 1) ) return;
    flag = CVodeGetSensNumNonlinSolvIters(cvode_mem, &nniS);
    if( _check_flag(&flag, "CVodeGetSensNumNonlinSolvIters", 1) ) return;
    flag = CVodeGetSensNumNonlinSolvConvFails(cvode_mem, &ncfnS);
    if( _check_flag(&flag, "CVodeGetSensNumNonlinSolvConvFails", 1) ) return;

    os << "   nfSe = " << nfSe << "   nfeS  = " << nfeS
       << "   nsetupsS = " << nsetupsS << "   netfS = " << netfS
       << "   nniS = " << nniS << "   ncfnS = " << ncfnS << std::endl
       << std::endl;
  }
}

template<typename U> inline void
ODESLVS_SUNDIALS::_print_interm
( const double t, const unsigned nx, const U*x, const std::string&var,
  std::ostream&os )
{
  os << " @t = " << std::scientific << std::setprecision(6)
                 << std::left << t << " :" << std::endl;
  for( unsigned ix=0; ix<nx; ix++ )
    os << " " << var.c_str() << "[" << ix << "] = " << x[ix] << std::endl;
  return;
}

template<typename U> inline void
ODESLVS_SUNDIALS::_print_interm
( const unsigned nx, const U*x, const std::string&var, std::ostream&os )
{
  if( !nx || !x ) return;
  os << std::scientific << std::setprecision(6) << std::left;
  for( unsigned ix=0; ix<nx; ix++ )
    os << " " << var.c_str() << "[" << ix << "] = " << x[ix] << std::endl;
  return;
}

} // end namescape mc

#endif
