// Copyright (C) Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

#ifndef MC__ODESLVS_CVODES_HPP
#define MC__ODESLVS_CVODES_HPP

#undef  MC__ODESLVS_CVODES_DEBUG

#include <sstream>
#include "odeslvs_base.hpp"
#include "odeslv_cvodes.hpp"

#define MC__ODESLVS_CVODES_USE_BAD
#undef  MC__ODESLVS_CVODES_DEBUG

namespace mc
{
//! @brief C++ class computing solutions of parametric ODEs with forward/adjoint sensitivity analysis capability using SUNDIALS and MC++.
////////////////////////////////////////////////////////////////////////
//! mc::ODESLV_CVODES is a C++ class for solution of IVPs in ODEs
//! with forward/adjoint sensitivity analysis capability using the code
//! CVODES in SUNDIALS and MC++.
////////////////////////////////////////////////////////////////////////
template <typename... ExtOps>
class ODESLVS_CVODES
: public virtual BASE_CVODES
, public virtual ODESLV_CVODES<ExtOps...>
, public virtual ODESLVS_BASE<ExtOps...>
, public virtual BASE_DE<ExtOps...>
{
 protected:
 
  using ODESLV_BASE<ExtOps...>::_print_interm;
  using ODESLV_BASE<ExtOps...>::_record;
  using ODESLV_BASE<ExtOps...>::_IC_D_QUAD;
  using ODESLV_BASE<ExtOps...>::_D2vec;
  using ODESLV_BASE<ExtOps...>::_Dx;
  using ODESLV_BASE<ExtOps...>::_Dq;
  using ODESLV_BASE<ExtOps...>::_Df;
      
  using ODESLVS_BASE<ExtOps...>::_nsmax;
  using ODESLVS_BASE<ExtOps...>::_istg;
  using ODESLVS_BASE<ExtOps...>::_t;
  using ODESLVS_BASE<ExtOps...>::_dT;
  using ODESLVS_BASE<ExtOps...>::_nx;
  using ODESLVS_BASE<ExtOps...>::_nq;
  using ODESLVS_BASE<ExtOps...>::_nf;
  using ODESLVS_BASE<ExtOps...>::_np;
  using ODESLVS_BASE<ExtOps...>::_Dfp;
  using ODESLVS_BASE<ExtOps...>::_ny;
  using ODESLVS_BASE<ExtOps...>::_Dy;
  using ODESLVS_BASE<ExtOps...>::_Dyq;
  using ODESLVS_BASE<ExtOps...>::_vIC;
  using ODESLVS_BASE<ExtOps...>::_vRHS;
  using ODESLVS_BASE<ExtOps...>::_vQUAD;
  using ODESLVS_BASE<ExtOps...>::_vFCT;

  using ODESLVS_BASE<ExtOps...>::_IC_SET_ASA;
  using ODESLVS_BASE<ExtOps...>::_CC_SET_ASA;
  using ODESLVS_BASE<ExtOps...>::_TC_SET_ASA;
  using ODESLVS_BASE<ExtOps...>::_RHS_SET_ASA;
  using ODESLVS_BASE<ExtOps...>::_IC_D_QUAD_ASA;
  using ODESLVS_BASE<ExtOps...>::_CC_D_QUAD_ASA;
  using ODESLVS_BASE<ExtOps...>::_TC_D_QUAD_ASA;
  using ODESLVS_BASE<ExtOps...>::_IC_SET_FSA;
  using ODESLVS_BASE<ExtOps...>::_CC_SET_FSA;
  using ODESLVS_BASE<ExtOps...>::_RHS_SET_FSA;
  using ODESLVS_BASE<ExtOps...>::_RHS_D_SET;
  using ODESLVS_BASE<ExtOps...>::_INI_D_SEN;
  using ODESLVS_BASE<ExtOps...>::_GET_D_SEN;
  using ODESLVS_BASE<ExtOps...>::_IC_D_SEN;
  using ODESLVS_BASE<ExtOps...>::_CC_D_SEN;
  using ODESLVS_BASE<ExtOps...>::_TC_D_SEN;
  using ODESLVS_BASE<ExtOps...>::_RHS_D_SEN;
  using ODESLVS_BASE<ExtOps...>::_RHS_D_QUAD;
  using ODESLVS_BASE<ExtOps...>::_JAC_D_SEN;
  using ODESLVS_BASE<ExtOps...>::_FCT_D_SEN;

  using ODESLV_CVODES<ExtOps...>::_cv_mem;
  using ODESLV_CVODES<ExtOps...>::_cv_flag;
  using ODESLV_CVODES<ExtOps...>::_Nx;
  using ODESLV_CVODES<ExtOps...>::_Nq;
  using ODESLV_CVODES<ExtOps...>::_pos_ic;
  using ODESLV_CVODES<ExtOps...>::_pos_rhs;
  using ODESLV_CVODES<ExtOps...>::_pos_quad;
  using ODESLV_CVODES<ExtOps...>::_pos_fct;
  using ODESLV_CVODES<ExtOps...>::_END_STA;
  using ODESLV_CVODES<ExtOps...>::_states;
  using ODESLV_CVODES<ExtOps...>::_states_stage;
  using ODESLV_CVODES<ExtOps...>::MC_CVRHS__;
  using ODESLV_CVODES<ExtOps...>::MC_CVQUAD__;
  using ODESLV_CVODES<ExtOps...>::MC_CVJAC__;
  using ODESLV_CVODES<ExtOps...>::_xk;
  using ODESLV_CVODES<ExtOps...>::_qk;
  using ODESLV_CVODES<ExtOps...>::_f;

  //! @brief Dense SUNMatrix for use in linear solves
  SUNMatrix _sun_matB;

  //! @brief dense SUNLinearSolver object for use by CVodeS
  SUNLinearSolver _sun_lsB;

  //! @brief SUNNonlinearSolver object for use by CVodeS adjoint
  SUNNonlinearSolver _sun_nlsB;

  //! @brief SUNNonlinearSolver object for use by CVodeS forward
  SUNNonlinearSolver _sun_nlsF;

  //! @brief current function index
  unsigned _ifct;

  //! @brief current parameter sensitivity index
  unsigned _isen;

  //! @brief size of N_vector arrays
  unsigned _nvec;

  //! @brief N_Vector object holding current adjoints
  N_Vector* _Ny;

  //! @brief N_Vector object holding current quadratures
  N_Vector* _Nyq;

  //! @brief pointer to array holding identifiers of the backward problems
  int* _indexB;
  
  //! @brief pointer to array holding identifiers of the backward problems
  unsigned* _iusrB;

  //! @brief position in _Ny
  unsigned _pos_adj;

  //! @brief position in _Nq
  unsigned _pos_adjquad; 

  //! @brief state sensitivity values at stage times
  std::vector< std::vector< std::vector< double > > > _xpk;

  //! @brief state adjoint values at stage times
  std::vector< std::vector< std::vector< double > > > _lk;

  //! @brief quadrature sensitivity values at stage times
  std::vector< std::vector< std::vector< double > > > _qpk;

  //! @brief function derivatives
  std::vector< std::vector< double > > _fp;

 public:

  using BASE_DE<ExtOps...>::set;

  /** @ingroup ODESLV
   *  @{
   */
  typedef typename BASE_DE<ExtOps...>::STATUS STATUS;
  typedef typename ODESLV_BASE<ExtOps...>::Results Results;
  typedef typename ODESLV_CVODES<ExtOps...>::Options Options;
  typedef typename ODESLV_CVODES<ExtOps...>::Exceptions Exceptions;
  using BASE_DE<ExtOps...>::np;
  using BASE_DE<ExtOps...>::nf;
  using ODESLV_CVODES<ExtOps...>::options;
  using ODESLV_CVODES<ExtOps...>::solve_state;
  using ODESLV_CVODES<ExtOps...>::val_state;
  using ODESLV_CVODES<ExtOps...>::val_quadrature;
  using ODESLV_CVODES<ExtOps...>::val_function;
  using ODESLV_CVODES<ExtOps...>::results_state;
  using ODESLV_CVODES<ExtOps...>::stats_state;
  
  //! @brief Default constructor
  ODESLVS_CVODES
    ();

  //! @brief Virtual destructor
  virtual ~ODESLVS_CVODES
    ();

  //! @brief Statistics for sensitivity/adjoint integration
  Stats stats_sensitivity;

  //! @brief Vector storing adjoint/sensitivity trajectyories (see Options::RESRECORD)
  std::vector< std::vector< Results > > results_sensitivity;

 //! @brief Propagate states and state-sensitivities forward in time through every time stages
  STATUS solve_sensitivity
    ( std::vector<double> const& p, std::ostream& os=std::cout );

 //! @brief Propagate states and state-sensitivities forward in time through every time stages
  STATUS solve_sensitivity
    ( double const* p, std::ostream& os=std::cout );

  //! @brief Propagate states and adjoints forward and backward in time through every time stages
  STATUS solve_adjoint
    ( std::vector<double> const& p, std::ostream& os=std::cout );

  //! @brief Propagate states and adjoints forward and backward in time through every time stages
  STATUS solve_adjoint
    ( double const* p, std::ostream& os=std::cout );

  //! @brief Setup local copy of parametric ODEs
  bool setup
    ()
    { return ODESLVS_BASE<ExtOps...>::_SETUP(); }

  //! @brief Setup local copy of parametric ODEs based on IVP
  bool setup
    ( ODESLVS_CVODES<ExtOps...> const& IVP )
    { return ODESLVS_BASE<ExtOps...>::_SETUP( IVP ); }

  //! @brief Record state and sensitivity trajectories in files <a>obndsta</a> and <a>obndsa</a>, with accuracy of <a>iprec</a> digits
  void record
    ( std::ofstream& obndsta, std::ofstream* obndsen, unsigned const iprec=5 )
    const
    { this->ODESLV_CVODES<ExtOps...>::record( obndsta, iprec );
      for( unsigned isen=0; isen<results_sensitivity.size(); ++isen )
        this->ODESLV_BASE<ExtOps...>::_record( obndsen[isen], results_sensitivity[isen], iprec ); }

  //! @brief Record state trajectories in files <a>obndsta</a>, with accuracy of <a>iprec</a> digits
  void record
    ( std::ofstream& obndsta, unsigned const iprec=5 )
    const
    { ODESLV_CVODES<ExtOps...>::record( obndsta, iprec ); }

  //! @brief Retreive state sensitivity values at stage times
  std::vector< std::vector< std::vector< double > > > const& val_state_sensitivity
    ()
    const
    { return _xpk; }

  //! @brief Retreive adjoint sensitivity values at stage times
  std::vector< std::vector< std::vector< double > > > const& val_state_adjoint
    ()
    const
    { return _lk; }

  //! @brief Retreive quadrature sensitivity values at stage times
  std::vector< std::vector< std::vector< double > > > const& val_quadrature_sensitivity
    ()
    const
    { return _qpk; }

  //! @brief Retreive function values at stage times
  std::vector< std::vector< double > > const& val_function_gradient
    ()
    const
    { return _fp; }
  /** @} */

 protected:
  //! @brief Propagate states and state-sensitivities forward in time through every time stages
  STATUS _states_FSA
    ( double const* p, std::ostream& os=std::cout );

  //! @brief Propagate states and adjoints forward and backward in time through every time stages
  STATUS _states_ASA
    ( double const* p, std::ostream& os=std::cout );

 private:
  //! @brief Function to initialize CVodes memory block (virtual)
  virtual bool _INI_CVODE
    ();

  //! @brief Function to initialize CVodeS memory block for forward sensitivity
  bool _INI_CVODES_FSA
    ();

  //! @brief Function to initialize CVodeS memory block for adjoint sensitivity
  bool _INI_CVODES_ASA
    ( unsigned const ifct, int& indexB, unsigned& iusrB );

  //! @brief Function to reinitialize CVodeS memory block for forward sensitivity
  bool _CC_CVODES_FSA
    ();

  //! @brief Function to reinitialize CVodeS memory block for forward quadrature sensitivity
  bool _CC_CVODES_QUAD
    ();

  //! @brief Function to reinitialize CVodeS memory block for adjoint sensitivity
  bool _CC_CVODES_ASA
    ( unsigned const ifct, int const indexB );

  //! @brief Function to finalize sensitivity/adjoint bounding
  void _END_SEN
    ();

  //! @brief Function to reinitialize sensitivity analysis
  bool _REINI_SEN
    ();

  //! @brief Function to initialize adjoint sensitivity analysis
  bool _INI_ASA
    ( double const* p );

  //! @brief Static wrapper to calculate the adjoint ODEs RHS derivatives
  static int MC_CVRHSB__
    ( sunrealtype t, N_Vector x, N_Vector y, N_Vector ydot, void* user_data );

  //! @brief Virtual function to calculate the adjoint ODEs RHS derivatives
  virtual int CVRHSB__
    ( sunrealtype t, N_Vector x, N_Vector y, N_Vector ydot, void* user_data );

  //! @brief Static wrapper to calculate the adjoint quadrature RHS derivatives
  static int MC_CVQUADB__
    ( sunrealtype t, N_Vector x, N_Vector y, N_Vector qdot, void* user_data );

  //! @brief Virtual function to calculate the adjoint quadrature RHS derivatives
  virtual int CVQUADB__
    ( sunrealtype t, N_Vector x, N_Vector y, N_Vector qdot, void* user_data );

  //! @brief Static wrapper to calculate the adjoint ODEs RHS Jacobian
  static int MC_CVJACB__
    ( sunrealtype t, N_Vector y, N_Vector yB, N_Vector fyB, SUNMatrix JacB,
      void* user_dataB, N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B );

  //! @brief Virtual function to calculate the adjoint ODEs RHS Jacobian
  virtual int CVJACB__
    ( sunrealtype t, N_Vector y, N_Vector yB, N_Vector fyB, SUNMatrix JacB,
      void* user_dataB, N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B );

  //! @brief Function to initialize forward sensitivity analysis
  bool _INI_FSA
    ( double const* p );

  //! @brief Static wrapper to calculate the sensitivity ODEs RHS derivatives
  static int MC_CVRHSF__
    ( int Ns, sunrealtype t, N_Vector x, N_Vector xdot, int is, N_Vector y,
      N_Vector ydot, void* user_data, N_Vector tmp1, N_Vector tmp2 );

  //! @brief Pure virtual function to calculate the sensitivity ODEs RHS derivatives
  virtual int CVRHSF__
    ( int Ns, sunrealtype t, N_Vector x, N_Vector xdot, int is, N_Vector y,
      N_Vector ydot, void* user_data, N_Vector tmp1, N_Vector tmp2 );

  //! @brief Static wrapper to calculate the sensitivity quadrature RHS derivatives
  static int MC_CVQUADF__
    ( int Ns, sunrealtype t, N_Vector x, N_Vector* y, N_Vector qdot, N_Vector* qSdot, 
      void *user_data, N_Vector tmp1, N_Vector tmp2 );

  //! @brief Pure virtual function to calculate the sensitivity quadrature RHS derivatives
  virtual int CVQUADF__
    ( int Ns, sunrealtype t, N_Vector x, N_Vector* y, N_Vector qdot, N_Vector* qSdot, 
      void *user_data, N_Vector tmp1, N_Vector tmp2 );

  //! @brief Private methods to block default compiler methods
  ODESLVS_CVODES( ODESLVS_CVODES<ExtOps...> const& ) = delete;
  ODESLVS_CVODES<ExtOps...>& operator=( ODESLVS_CVODES<ExtOps...> const& ) = delete;
};

template <typename... ExtOps>
ODESLVS_CVODES<ExtOps...>::ODESLVS_CVODES
()
: _sun_matB(nullptr), _sun_lsB(nullptr), _sun_nlsB(nullptr), _sun_nlsF(nullptr),
  _ifct(0), _isen(0), _nvec(0), _Ny(nullptr), _Nyq(nullptr),
  _indexB(nullptr), _iusrB(nullptr)
{}

template <typename... ExtOps>
ODESLVS_CVODES<ExtOps...>::~ODESLVS_CVODES
()
{
  if( _Ny )  N_VDestroyVectorArray( _Ny,  _nvec );
  if( _Nyq ) N_VDestroyVectorArray( _Nyq, _nvec );
  delete[] _indexB;
  delete[] _iusrB;
  if( _sun_nlsF ) SUNNonlinSolFree( _sun_nlsF ); /* Free the nonlinear solver memory */
  if( _sun_nlsB ) SUNNonlinSolFree( _sun_nlsB ); /* Free the nonlinear solver memory */
  if( _sun_lsB )  SUNLinSolFree( _sun_lsB );     /* Free the linear solver memory */
  if( _sun_matB ) SUNMatDestroy( _sun_matB );    /* Free the matrix memory */
}

template <typename... ExtOps>
bool
ODESLVS_CVODES<ExtOps...>::_INI_CVODE
()
{
  // Call _INI_CVODE in ODEBND_CVODES
  this->ODESLV_CVODES<ExtOps...>::_INI_CVODE();

  // Allocate memory for adjoint integration
  _cv_flag = CVodeAdjInit( _cv_mem, options.ASACHKPT, options.ASAINTERP );
  if( _check_cv_flag( &_cv_flag, "CVodeAdjInit", 1 ) ) return false;

  // Reinitialize adjoint holding vectors
  delete[] _indexB; _indexB = new int[_nf];
  delete[] _iusrB;  _iusrB  = new unsigned[_nf];

  return true;
}

template <typename... ExtOps>
bool
ODESLVS_CVODES<ExtOps...>::_INI_CVODES_ASA
( unsigned const ifct, int& indexB, unsigned& iusrB )
{
  // Create CVodeS memory block for the ADAMS or BDF method
  _cv_flag = CVodeCreateB( _cv_mem, options.INTMETH, &indexB );
  if( _check_cv_flag( &_cv_flag, "CVodeCreateB", 1 ) ) return false;

  // Specify error output
//  if( options.DISPLAY < 0 )
//    _cv_flag = CVodeSetErrFile( _cv_mem, NULL );
//  else
//    _cv_flag = CVodeSetErrFile( _cv_mem, stderr );
//  if( _check_cv_flag( &_cv_flag, "CVodeSetErrFile", 1 ) ) return false;

  // Initialize CVodeS memory and specify the adjoint RHS function,
  // terminal time _t, and terminal adjoint _Ny
  _cv_flag = CVodeInitB( _cv_mem, indexB, MC_CVRHSB__, _t, _Ny[ifct] );
  if( _check_cv_flag( &_cv_flag, "CVodeInitB", 1 ) ) return false;

  // Specify the user_data to pass the function index corresponding to
  // the RHS function
  iusrB = _ifct;
  _cv_flag = CVodeSetUserDataB( _cv_mem, indexB, &iusrB );
  if( _check_cv_flag( &_cv_flag, "CVodeSetUserDataB", 1 ) ) return false;

  // Specify the nonlinear solver
  if( !ifct ){
    if( _sun_nlsB ){ SUNNonlinSolFree( _sun_nlsB );  _sun_nlsB = nullptr; } /* Free the nonlinear solver memory */
    if( _sun_lsB ) { SUNLinSolFree( _sun_lsB );      _sun_lsB  = nullptr; } /* Free the linear solver memory */
    if( _sun_matB ){ SUNMatDestroy( _sun_matB );     _sun_matB = nullptr; } /* Free the matrix memory */
  }
  switch( options.NLINSOL ){
   // Fixed point nonlinear solver
   case Options::FIXEDPOINT:
    if( !ifct ){
      _sun_nlsB = SUNNonlinSol_FixedPoint( _Ny[ifct], 0, sunctx );
      if( _check_cv_flag( (void *)_sun_nlsB, "SUNNonlinSol_FixedPoint", 0 ) ) return false;
    }
    break;
   
   // Newton nonlinear solver
   case Options::NEWTON:
    // Specify the linear solver and Jacobian approximation
    switch( options.LINSOL ){
     case Options::DIAG: default:
       _cv_flag = CVDiagB( _cv_mem, indexB );
       if( _check_cv_flag( &_cv_flag, "CVDiag", 1) ) return false;
       break;

     // Dense Jacobian
     case Options::DENSE:
     case Options::DENSEDQ:
       if( !ifct ){
         // Create dense SUNMatrix for use in linear solves
         _sun_matB = SUNDenseMatrix( NV_LENGTH_S( _Ny[ifct] ), NV_LENGTH_S( _Ny[ifct] ), sunctx );
         if( _check_cv_flag( (void*)_sun_matB, "SUNDenseMatrix", 0 ) ) return false;
         // Create dense SUNLinearSolver object for use by CVode
         _sun_lsB = SUNLinSol_Dense( _Ny[ifct], _sun_matB, sunctx );
         if( _check_cv_flag( (void *)_sun_lsB, "SUNLinSol_Dense", 0 ) ) return false;
       }
       // Attach the matrix and linear solver
       _cv_flag = CVodeSetLinearSolverB( _cv_mem, indexB, _sun_lsB, _sun_matB );
       if( _check_cv_flag( &_cv_flag, "CVodeSetLinearSolverB", 1 ) ) return false;
       // Set the user-supplied Jacobian routine Jac
       _cv_flag = CVodeSetJacFnB( _cv_mem, indexB, options.LINSOL==Options::DENSE? MC_CVJACB__: nullptr );
       if ( _check_cv_flag( &_cv_flag, "CVodeSetJacFnB", 1 ) ) return false;
       break;
    }
    if( !ifct ){
      _sun_nlsB = SUNNonlinSol_Newton( _Ny[ifct], sunctx );
      if( _check_cv_flag( (void *)_sun_nlsB, "SUNNonlinSol_Newton", 0 ) ) return false;
    }
    break;
  }
  _cv_flag = CVodeSetNonlinearSolverB( _cv_mem, indexB, _sun_nlsB );
  if( _check_cv_flag( &_cv_flag, "CVodeSetNonlinearSolverB", 1 ) ) return false;

  // Specify the relative and absolute tolerances for states
  _cv_flag = CVodeSStolerancesB( _cv_mem, indexB, options.RTOLB, options.ATOLB );
  if( _check_cv_flag( &_cv_flag, "CVodeSStolerancesB", 1 ) ) return false;

  // Specify minimum stepsize
  _cv_flag = CVodeSetMinStepB( _cv_mem, indexB, options.HMIN>0.? options.HMIN:0. );
  if( _check_cv_flag( &_cv_flag, "CVodeSetMinStepB", 1 ) ) return false;

  // Specify maximum stepsize
  _cv_flag = CVodeSetMaxStepB( _cv_mem, indexB, options.HMAX>0.? options.HMAX: 0. );
  if( _check_cv_flag( &_cv_flag, "CVodeSetMaxStepB", 1 ) ) return false;

  // Specify maximum number of steps between two stage times
  _cv_flag = CVodeSetMaxNumStepsB( _cv_mem, indexB, options.NMAX );
  if( _check_cv_flag( &_cv_flag, "CVodeSetMaxNumStepsB", 1 ) ) return false;

  // Initialize the integrator memory for the quadrature variables
  if( !_Nyq ) return true;
  _cv_flag = CVodeQuadInitB( _cv_mem, indexB, MC_CVQUADB__, _Nyq[ifct] );
  if( _check_cv_flag( &_cv_flag, "CVodeQuadInitB", 1 ) ) return false;

  // Specify whether or not to perform error control on quadrature
  _cv_flag = CVodeSetQuadErrConB( _cv_mem, indexB, options.QERRB );
  if( _check_cv_flag( &_cv_flag, "CVodeSetQuadErrConB", 1 ) ) return false;

  // Specify the relative and absolute tolerances for quadratures
  _cv_flag = CVodeQuadSStolerancesB( _cv_mem, indexB, options.RTOLB, options.ATOLB );
  if( _check_cv_flag( &_cv_flag, "CVodeQuadSStolerancesB", 1 ) ) return false;

  return true;
}

template <typename... ExtOps>
bool
ODESLVS_CVODES<ExtOps...>::_INI_CVODES_FSA
()
{
  // Allocate memory for sensitivity integration
  _cv_flag = CVodeSensInit1( _cv_mem, _np, options.FSACORR, MC_CVRHSF__, _Ny );
  if( _check_cv_flag( &_cv_flag, "CVodeSensInit1", 1 ) ) return false;

  // Specify error output
//  if( options.DISPLAY < 0 )
//    _cv_flag = CVodeSetErrFile( _cv_mem, NULL );
//  else
//    _cv_flag = CVodeSetErrFile( _cv_mem, stderr );
//  if( _check_cv_flag( &_cv_flag, "CVodeSetErrFile", 1 ) ) return false;

  // Specify absolute and relative tolerances for sensitivities
  if( options.AUTOTOLS ){
    _cv_flag = CVodeSensEEtolerances( _cv_mem );
    if( _check_cv_flag( &_cv_flag, "CVodeSensEEtolerances", 1) ) return false;
  }
  else{
    std::vector<sunrealtype> ATOLS( _np, options.ATOLS );
    _cv_flag = CVodeSensSStolerances( _cv_mem, options.RTOLS, ATOLS.data() );
    if( _check_cv_flag( &_cv_flag, "CVodeSensSStolerances", 1) ) return false;
  }

  // Specify the error control strategy for sensitivity variables
  _cv_flag = CVodeSetSensErrCon( _cv_mem, options.FSAERR );
  if( _check_cv_flag( &_cv_flag, "CVodeSetSensErrCon", 1 ) ) return false;

  // Specify problem parameter information for sensitivity calculations
  _cv_flag = CVodeSetSensParams( _cv_mem, 0, 0, 0 );
  if( _check_cv_flag( &_cv_flag, "CVodeSetSensParams", 1 ) ) return false;

  // Specify the nonlinear solver for sensitivity calculations
  if( _sun_nlsF ){ SUNNonlinSolFree( _sun_nlsF );  _sun_nlsF = nullptr; } /* Free the nonlinear solver memory */
  switch( options.NLINSOL ){
   // Fixed point nonlinear solver
   case Options::FIXEDPOINT:
    switch( options.FSACORR ){
     case Options::SIMULTANEOUS:
      _sun_nlsF = SUNNonlinSol_FixedPointSens( _np+1, _Nx, 0, sunctx);
      break;
     case Options::STAGGERED:
      _sun_nlsF = SUNNonlinSol_FixedPointSens( _np, _Nx, 0, sunctx);
      break;
     case Options::STAGGERED1:
      _sun_nlsF = SUNNonlinSol_FixedPoint( _Nx, 0, sunctx);
      break;
    }
    break;
   
   // Newton nonlinear solver
   case Options::NEWTON:
    switch( options.FSACORR ){
     case Options::SIMULTANEOUS:
      _sun_nlsF = SUNNonlinSol_NewtonSens( _np+1, _Nx, sunctx);
      break;
     case Options::STAGGERED:
      _sun_nlsF = SUNNonlinSol_NewtonSens( _np, _Nx, sunctx);
      break;
     case Options::STAGGERED1:
      _sun_nlsF = SUNNonlinSol_Newton( _Nx, sunctx);
      break;
    }
    if( _check_cv_flag( (void *)_sun_nlsF, "SUNNonlinSol_FixedPointSens", 0 ) ) return false;
  }
  
  // Attach sensitivity nonlinear solver to CVodeS
  switch( options.FSACORR ){
   case Options::SIMULTANEOUS:
    _cv_flag = CVodeSetNonlinearSolverSensSim( _cv_mem, _sun_nlsF );
    if( _check_cv_flag( &_cv_flag, "CVodeSetNonlinearSolverSensSim", 1 ) ) return false;
    break;
   case Options::STAGGERED:
    _cv_flag = CVodeSetNonlinearSolverSensStg( _cv_mem, _sun_nlsF );
    if( _check_cv_flag( &_cv_flag, "CVodeSetNonlinearSolverSensStg", 1 ) ) return false;
    break;
   case Options::STAGGERED1:
    _cv_flag = CVodeSetNonlinearSolverSensStg1( _cv_mem, _sun_nlsF );
    if( _check_cv_flag( &_cv_flag, "CVodeSetNonlinearSolverSensStg1", 1 ) ) return false;
    break;
  }
  
  // Initialize integrator memory for quadratures
  if( !_nq ) return true;
  _cv_flag = CVodeQuadSensInit( _cv_mem, MC_CVQUADF__, _Nyq );
  if( _check_cv_flag( &_cv_flag, "CVodeQuadSensInit", 1 ) ) return false;
  
  // Specify whether or not to perform error control on quadrature
  _cv_flag = CVodeSetQuadSensErrCon( _cv_mem, options.QERRS );
  if( _check_cv_flag( &_cv_flag, "CVodeSetQuadSensErrCon", 1 ) ) return false;

  // Specify absolute and relative tolerances for quadratures
  if( options.AUTOTOLS ){
    _cv_flag = CVodeQuadSensEEtolerances( _cv_mem );
    if( _check_cv_flag( &_cv_flag, "CVodeQuadSensEEtolerances", 1 ) ) return false;
  }
  else{
    std::vector<sunrealtype> ATOLS( _np, options.ATOLS );
    _cv_flag = CVodeQuadSensSStolerances( _cv_mem, options.RTOLS, ATOLS.data() );
    if( _check_cv_flag( &_cv_flag, "CVodeQuadSensSStolerances", 1 ) ) return false;
  }

  return true;
}

template <typename... ExtOps>
bool
ODESLVS_CVODES<ExtOps...>::_CC_CVODES_ASA
( unsigned const ifct, int const indexB )
{

  // std::cout << "Calling CVodeReInitB @t=" << _t << std::endl;
  // Reinitialize CVodeS memory block for current time _t and adjoint _Ny
  _cv_flag = CVodeReInitB( _cv_mem, indexB, _t, _Ny[ifct] );
  if( _check_cv_flag( &_cv_flag, "CVodeReInitB", 1 ) ) return false;

  // Reinitialize CVodeS memory block for current adjoint quarature _Nyq
  if( !_np ) return true;
  _cv_flag = CVodeQuadReInitB( _cv_mem, indexB, _Nyq[ifct] );
  if( _check_cv_flag( &_cv_flag, "CVodeQuadReInitB", 1 ) ) return false;

  return true;
}

template <typename... ExtOps>
bool
ODESLVS_CVODES<ExtOps...>::_CC_CVODES_FSA
()
{
  // Reinitialize CVodeS memory block for current sensitivity _Ny
  _cv_flag = CVodeSensReInit( _cv_mem, options.FSACORR, _Ny );
  if( _check_cv_flag( &_cv_flag, "CVodeSensReInit", 1 ) ) return false;

  return true;
}

template <typename... ExtOps>
bool
ODESLVS_CVODES<ExtOps...>::_CC_CVODES_QUAD
()
{
  // Reinitialize CVode memory block for current sensitivity quarature _Nyq
  if( !_nq ) return true;
  _cv_flag = CVodeQuadSensReInit( _cv_mem, _Nyq );
  if( _check_cv_flag( &_cv_flag, "CVodeQuadSensReInit", 1 ) ) return false;

  return true;
}

template <typename... ExtOps>
void
ODESLVS_CVODES<ExtOps...>::_END_SEN
()
{
  // Get final CPU time
  _final_stats( stats_sensitivity );
}

template <typename... ExtOps>
bool
ODESLVS_CVODES<ExtOps...>::_REINI_SEN
()
{
  // reset at time stages
  _xpk.clear(); _xpk.reserve(_nsmax);
  _lk.clear();  _lk.reserve(_nsmax);
  _qpk.clear(); _qpk.reserve(_nsmax);
  _fp.clear();  _fp.reserve(_nf);

  return true;
}

template <typename... ExtOps>
inline
bool
ODESLVS_CVODES<ExtOps...>::_INI_ASA
( double const* p )
{
  // Initialize bound propagation
  if( !_INI_D_SEN( p, _nf, _np ) || !_REINI_SEN() )
    return false;

  // Set SUNDIALS adjoint/quadrature arrays
  if( _nvec != _nf ){
    if( _Ny )   N_VDestroyVectorArray( _Ny,  _nvec );
    if( _Nyq )  N_VDestroyVectorArray( _Nyq, _nvec );
    _nvec = _nf;
    _Ny  = N_VCloneVectorArray( _nvec, _Nx );
    _Nyq = N_VCloneVectorArray( _nvec, _Nx );
  }
  for( unsigned i=0; i<_nf; i++ ){
    if( !_Ny[i] || NV_LENGTH_S( _Ny[i] ) != _ny ){
      if( _Ny[i] ) N_VDestroy( _Ny[i] );
      _Ny[i] = N_VNew_Serial( _ny, sunctx );
    }
    if( !_Nyq[i] || NV_LENGTH_S( _Nyq[i] ) != (sunindextype)_np ){
      if( _Nyq[i] ) N_VDestroy( _Nyq[i] );
      _Nyq[i] = N_VNew_Serial( _np, sunctx );
    }
  }

  // Reset result record and statistics
  results_sensitivity.clear();
  results_sensitivity.resize( _nf );
  _init_stats( stats_sensitivity );

  return true;
}

template <typename... ExtOps>
inline
int
ODESLVS_CVODES<ExtOps...>::MC_CVRHSB__
( sunrealtype t, N_Vector x, N_Vector y, N_Vector ydot, void* user_data )
{
#ifdef MC__BASE_CVODES_CHECK
  //std::cout << "BASE_CVODES::PTR_BASE_CVODES: " << BASE_CVODES::PTR_BASE_CVODES
  //          << "  PTR_CVRHSB: " << PTR_CVRHSB << std::endl;
  assert( BASE_CVODES::PTR_BASE_CVODES != nullptr && PTR_CVRHSB != nullptr);
#endif
  return (BASE_CVODES::PTR_BASE_CVODES->*PTR_CVRHSB)( t, x, y, ydot, user_data );
}

template <typename... ExtOps>
inline
int
ODESLVS_CVODES<ExtOps...>::CVRHSB__
( sunrealtype t, N_Vector x, N_Vector y, N_Vector ydot, void* user_data )
{
  _ifct = *static_cast<unsigned*>( user_data );
#ifdef MC__ODESLVS_CVODES_DEBUG
  std::cout << std::scientific << std::setprecision(6) << t;
  for( unsigned i=0; i<NV_LENGTH_S( x ); i++ ) std::cout << "  " << NV_Ith_S( x, i );
  for( unsigned i=0; i<NV_LENGTH_S( y ); i++ ) std::cout << "  " << NV_Ith_S( y, i );
#endif
  bool flag = _RHS_D_SEN( t, NV_DATA_S( x ), NV_DATA_S( y ), NV_DATA_S( ydot ), _ifct );
#ifdef MC__ODESLVS_CVODES_DEBUG
  std::cout << "  " << _ifct;
  for( unsigned i=0; i<NV_LENGTH_S( ydot ); i++ ) std::cout << "  " << NV_Ith_S( ydot, i );
  std::cout << std::endl;
  { int dum; std::cin >> dum; }
#endif
  stats_sensitivity.numRHS++;
  return( flag? 0: -1 );
}

template <typename... ExtOps>
inline
int
ODESLVS_CVODES<ExtOps...>::MC_CVJACB__
( sunrealtype t, N_Vector x, N_Vector y, N_Vector fy, SUNMatrix Jac,
  void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3 )
{
#ifdef MC__BASE_CVODES_CHECK
  //std::cout << "BASE_CVODES::PTR_BASE_CVODES: " << BASE_CVODES::PTR_BASE_CVODES
  //          << "  PTR_CVJACB: " << PTR_CVJACB << std::endl;
  assert( BASE_CVODES::PTR_BASE_CVODES != nullptr && PTR_CVJACB != nullptr);
#endif
  return (BASE_CVODES::PTR_BASE_CVODES->*PTR_CVJACB)( t, x, y, fy, Jac, user_data, tmp1, tmp2, tmp3 );
}

template <typename... ExtOps>
inline
int
ODESLVS_CVODES<ExtOps...>::CVJACB__
( sunrealtype t, N_Vector x, N_Vector y, N_Vector fy, SUNMatrix Jac,
  void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3 )
{
  //std::cerr << "Entering: ODESLV_CVODES<ExtOps...>::MC_CVJACB__\n";
  bool flag = false;
  switch( options.LINSOL ){
   case Options::DIAG:
   case Options::DENSEDQ:
   //case Options::SPARSE:
    flag = false;
    break;
   case Options::DENSE:
    flag = _JAC_D_SEN( t, NV_DATA_S( x ), NV_DATA_S( y ), SM_COLS_D(Jac) );
    break;
  }
  stats_sensitivity.numJAC++; // increment JAC counter
  return( flag? 0: -1 );
}

template <typename... ExtOps>
inline
int
ODESLVS_CVODES<ExtOps...>::MC_CVQUADB__
( sunrealtype t, N_Vector x, N_Vector y, N_Vector qdot, void* user_data )
{
#ifdef MC__BASE_CVODES_CHECK
  //std::cout << "BASE_CVODES::PTR_BASE_CVODES: " << BASE_CVODES::PTR_BASE_CVODES
  //          << "  PTR_CVQUADB: " << PTR_CVQUADB << std::endl;
  assert( BASE_CVODES::PTR_BASE_CVODES != nullptr && PTR_CVQUADB != nullptr);
#endif
  return (BASE_CVODES::PTR_BASE_CVODES->*PTR_CVQUADB)( t, x, y, qdot, user_data );
}

template <typename... ExtOps>
inline
int
ODESLVS_CVODES<ExtOps...>::CVQUADB__
( sunrealtype t, N_Vector x, N_Vector y, N_Vector qdot, void* user_data )
{
  _ifct = *static_cast<unsigned*>( user_data );
  bool flag = _RHS_D_QUAD( _np, NV_DATA_S( qdot ), _ifct );
  return( flag? 0: -1 );
}

//! @fn template <typename... ExtOps> inline typename ODESLVS_CVODES<ExtOps...>::STATUS ODESLVS_CVODES<ExtOps...>::solve_adjoint
//!( std::vector<double> const& p, std::ostream& os=std::cout )
//!
//! This function computes a solution to the parametric ODEs with adjoint
//! sensitivity analysis:
//!  - <a>p</a>  [input]  parameter values
//!  - <a>os</a> [input]  output stream [default: std::cout]
//! .
//! The return value is the status.
template <typename... ExtOps>
typename ODESLVS_CVODES<ExtOps...>::STATUS
ODESLVS_CVODES<ExtOps...>::solve_adjoint
( std::vector<double> const& p, std::ostream& os )
{
  registration();
  STATUS flag = _states_ASA( p.data(), os );
  unregistration();
  return flag;
}

//! @fn template <typename... ExtOps> inline typename ODESLVS_CVODES<ExtOps...>::STATUS ODESLVS_CVODES<ExtOps...>::solve_adjoint
//!( double const* p, std::ostream& os=std::cout )
//!
//! This function computes a solution to the parametric ODEs with adjoint
//! sensitivity analysis:
//!  - <a>p</a>  [input]  parameter values
//!  - <a>os</a> [input]  output stream [default: std::cout]
//! .
//! The return value is the status.
template <typename... ExtOps>
typename ODESLVS_CVODES<ExtOps...>::STATUS
ODESLVS_CVODES<ExtOps...>::solve_adjoint
( double const* p, std::ostream& os )
{
  registration();
  STATUS flag = _states_ASA( p, os );
  unregistration();
  return flag;
}

template <typename... ExtOps>
typename ODESLVS_CVODES<ExtOps...>::STATUS
ODESLVS_CVODES<ExtOps...>::_states_ASA
( double const* p, std::ostream& os )
{
  // Compute state bounds and store intermediate results
  STATUS flag = STATUS::NORMAL;
  flag = ODESLV_CVODES<ExtOps...>::_states( p, true, os );
  if( flag != STATUS::NORMAL ) return flag;

  // Nothing to do if no functions or parameters are defined
  if( !_nf || !_np ) return STATUS::NORMAL;

  try{
    // Initialize adjoint integration
    if( !_INI_ASA( p ) ) return STATUS::FATAL;
    _t = _dT[_nsmax];
    const unsigned NSTEP = options.RESRECORD? options.RESRECORD: 1;

    // Terminal adjoint & quadrature values
    _lk.resize( _nsmax+1, std::vector<std::vector<double>>( _nf ) );
    _qpk.resize( _nsmax+1, std::vector<std::vector<double>>( _nf ) );
//    if( lk && !lk[_nsmax] ) lk[_nsmax] = new double[(_nx+_np)*_nf];
    _pos_fct = ( _vFCT.size()>=_nsmax? _nsmax-1:0 );
    for( _ifct=0; _ifct < _nf; _ifct++ ){
      if( !_TC_SET_ASA( _pos_fct, _ifct )
       || !_TC_D_SEN( _t, _xk[_nsmax].data(), NV_DATA_S(_Ny[_ifct]) )
       || ( _Nyq && _Nyq[_ifct] && !_TC_D_QUAD_ASA( NV_DATA_S(_Nyq[_ifct]) ) ) )
        { _END_SEN(); return STATUS::FATAL; }
      _GET_D_SEN( NV_DATA_S(_Ny[_ifct]), _np, _Nyq? NV_DATA_S(_Nyq[_ifct]): 0 );
      for( unsigned iq=0; iq<_np; iq++ )
        _Dfp[iq*_nf+_ifct] = _Dyq[iq];

      // Display / record / return adjoint terminal values
      _lk[_nsmax].push_back( std::vector<double>( _Dy, _Dy+_ny ) );
      _qpk[_nsmax].push_back( std::vector<double>( _Dyq, _Dyq+_np ) );
      if( options.DISPLAY >= 1 ){
        std::ostringstream ol; ol << " l[" << _ifct << "]";
        if( !_ifct ) _print_interm( _dT[_nsmax], _nx, _Dy, ol.str(), os );
        else         _print_interm( _nx, _Dy, ol.str(), os );
        std::ostringstream oq; oq << " qp[" << _ifct << "]";
        _print_interm( _np, _Dyq, oq.str(), os );
      }
      if( options.RESRECORD )
        results_sensitivity[_ifct].push_back( Results( _t, _nx, _Dy, _np, _Dyq ) );
//      for( unsigned iy=0; lk && iy<_ny+_np; iy++ )
//        lk[_nsmax][_ifct*(_nx+_np)+iy] = iy<_ny? _Dy[iy]: _Dyq[iy-_ny];
    }

    // Initialization of adjoint integration
    for( _ifct=0; _ifct < _nf; _ifct++ )
      if( !_INI_CVODES_ASA( _ifct, _indexB[_ifct], _iusrB[_ifct] ) )
        { _END_SEN(); return STATUS::FATAL;}

    // Integrate adjoint ODEs through each stage using SUNDIALS
    for( _istg=_nsmax; _istg>0; _istg-- ){

      // Is forward state evaluation needed after discontinuity?
      if( _istg<_nsmax && _vIC.size()>=_nsmax ){
        _t = _dT[_istg-1];
        _D2vec( _xk[_istg-1].data(), _nx, NV_DATA_S( _Nx ) );
        _IC_D_QUAD( NV_DATA_S( _Nq ) );
        std::cout << "RESTARTING FORWARD INTEGRATION at t=" << _t << std::endl;
        _states_stage( _istg-1, _t, _Nx, _Nq, true, true, false, os );
      }

      // Update list of operations in RHSADJ and QUADADJ
      _pos_rhs  = ( _vRHS.size() <=1? 0: _istg-1 );
      _pos_quad = ( _vQUAD.size()<=1? 0: _istg-1 );
      _pos_fct  = ( _vFCT.size()>=_nsmax? _istg: ( _vFCT.size()==1 && _istg==_nsmax? 1: 0 ) );   
      //_vFCT.size()>=_nsmax? _istg: 0 );
#ifdef MC__ODESLVS_CVODES_DEBUG
      std::cout << "pos_fct: " << _pos_fct << std::endl;
#endif
      if( !ODESLV_BASE<ExtOps...>::_RHS_D_SET( _pos_rhs, _pos_quad )
       || !_RHS_SET_ASA( _pos_rhs, _pos_quad, _pos_fct )
       || !_RHS_D_SET( _nf, _np ) )
        { _END_SEN(); return STATUS::FATAL; }

      // Propagate bounds backward to previous stage time
      _t = _dT[_istg];
      const double TSTEP = ( _t - _dT[_istg-1] ) / NSTEP;
      double TSTOP = _t-TSTEP;
      for( unsigned k=0; k<NSTEP; k++, TSTOP-=TSTEP ){
        if( k+1 == NSTEP ) TSTOP = _dT[_istg-1];
        _cv_flag = CVodeB( _cv_mem, TSTOP, CV_NORMAL );
        if( _check_cv_flag( &_cv_flag, "CVodeB", 1 ) )
          { _END_SEN(); return STATUS::FATAL; }

        // intermediate record
        if( options.RESRECORD ){
          for( _ifct=0; _ifct < _nf; _ifct++ ){
            _cv_flag = CVodeGetB( _cv_mem, _indexB[_ifct], &_t, _Ny[_ifct]);
            if( _check_cv_flag( &_cv_flag, "CVodeGetB", 1) )
              { _END_SEN(); return STATUS::FATAL; }
#ifdef MC__ODESLVS_CVODES_DEBUG
            std::cout << "Adjoint #" << _ifct << ": " << _t << std::endl;
            for( unsigned iy=0; iy<_ny; iy++ )
              std::cout << "_Ny" << _ifct << "[" << iy << "] = " << NV_Ith_S(_Ny[_ifct],iy) << std::endl;
#endif
            _cv_flag = CVodeGetQuadB( _cv_mem, _indexB[_ifct], &_t, _Nyq[_ifct]);
            if( _check_cv_flag( &_cv_flag, "CVodeGetQuadB", 1) )
              { _END_SEN(); return STATUS::FATAL; }
#ifdef MC__ODESLVS_CVODES_DEBUG
            for( unsigned ip=0; ip<_np; ip++ )
              std::cout << "_Nyq" << _ifct << "[" << ip << "] = " << NV_Ith_S(_Nyq[_ifct],ip) << std::endl;
#endif
            results_sensitivity[_ifct].push_back( Results( _t, _nx, NV_DATA_S(_Ny[_ifct]), _np, _Nyq && _Nyq[_ifct]? NV_DATA_S(_Nyq[_ifct]): 0 ) );
          }
        }
      }
      for( _ifct=0; _ifct < _nf; _ifct++ ){
        void *cv_memB = CVodeGetAdjCVodeBmem( _cv_mem, _indexB[_ifct] );
        long int nstpB;
        _cv_flag = CVodeGetNumSteps( cv_memB, &nstpB );
        stats_sensitivity.numSteps += nstpB;
#ifdef MC__ODESLVS_CVODES_DEBUG
        std::cout << "Number of steps for adjoint #" << _ifct << ": " 
                  << nstpB << std::endl;
#endif
      }

      // states/adjoints/quadratures at stage time
#ifdef MC__ODESLVS_CVODES_DEBUG
      for( unsigned ix=0; ix<_nx; ix++ )
        std::cout << "xk[" << _istg-1 << "][" << ix << "] = " << _xk[_istg-1][ix] << std::endl;
#endif
//      if( lk && !lk[_istg-1] ) lk[_istg-1] = new double[(_nx+_np)*_nf];
      for( _ifct=0; _ifct < _nf; _ifct++ ){
        _cv_flag = CVodeGetB( _cv_mem, _indexB[_ifct], &_t, _Ny[_ifct]);
        if( _check_cv_flag( &_cv_flag, "CVodeGetB", 1) )
          { _END_SEN(); return STATUS::FATAL; }
#ifdef MC__ODESLVS_CVODES_DEBUG
        for( unsigned iy=0; iy<_ny; iy++ )
          std::cout << "_Ny" << _ifct << "[" << iy << "] = " << NV_Ith_S(_Ny[_ifct],iy) << std::endl;
#endif
        _cv_flag = CVodeGetQuadB( _cv_mem, _indexB[_ifct], &_t, _Nyq[_ifct]);
        if( _check_cv_flag( &_cv_flag, "CVodeGetQuadB", 1) )
          { _END_SEN(); return STATUS::FATAL; }
#ifdef MC__ODESLVS_CVODES_DEBUG
        for( unsigned ip=0; ip<_np; ip++ )
          std::cout << "_Nyq" << _ifct << "[" << ip << "] = " << NV_Ith_S(_Nyq[_ifct],ip) << std::endl;
#endif
        // Add function contribution to adjoint values (discontinuities)
        if( _istg > 1  ){
          _pos_ic  = ( _vIC.size() >=_nsmax? _istg-1: 0 );
          _pos_fct = ( _vFCT.size()>=_nsmax? _istg-1: 0 );
          if( ( _pos_fct || _pos_ic )
           && ( !_CC_SET_ASA( _pos_ic, _pos_fct, _ifct )
             || !_CC_D_SEN( _t, _xk[_istg-1].data(), NV_DATA_S(_Ny[_ifct]) )
             || ( _Nyq && _Nyq[_ifct] && !_CC_D_QUAD_ASA( NV_DATA_S(_Nyq[_ifct]) ) ) ) )
            { _END_SEN(); return STATUS::FATAL; }
          //else if( !_pos_fct )
            _GET_D_SEN( NV_DATA_S(_Ny[_ifct]), _np, _Nyq && _Nyq[_ifct]? NV_DATA_S(_Nyq[_ifct]): 0 );
#ifdef MC__ODESLVS_CVODES_DEBUG
          for( unsigned iy=0; iy<_ny; iy++ )
            std::cout << "_Ny" << _ifct << "[" << iy << "] = " << NV_Ith_S(_Ny[_ifct],iy) << std::endl;
          for( unsigned ip=0; ip<_np; ip++ )
            std::cout << "_Nyq" << _ifct << "[" << ip << "] = " << NV_Ith_S(_Nyq[_ifct],ip) << std::endl;
#endif

          // Reset ODE solver - needed in case of discontinuity
          if( !_CC_CVODES_ASA( _ifct, _indexB[_ifct] ) )
            { _END_SEN(); return STATUS::FATAL; }
        }
        
        // Add initial state contribution to function derivatives 
        else if( !_IC_SET_ASA()
              || !_IC_D_SEN( _t, _xk[_istg-1].data(), NV_DATA_S(_Ny[_ifct]) )
              || ( _Nyq && _Nyq[_ifct] && !_IC_D_QUAD_ASA( NV_DATA_S(_Nyq[_ifct]) ) ) )
          { _END_SEN(); return STATUS::FATAL; }

        // Display / record / return adjoint terminal values
        _lk[_istg-1].push_back( std::vector<double>( _Dy, _Dy+_ny ) );
        _qpk[_istg-1].push_back( std::vector<double>( _Dyq, _Dyq+_np ) );
        if( options.DISPLAY >= 1 ){
          std::ostringstream ol; ol << " l[" << _ifct << "]";
          if( !_ifct ) _print_interm( _dT[_istg-1], _nx, _Dy, ol.str(), os );
          else         _print_interm( _nx, _Dy, ol.str(), os );
          std::ostringstream oq; oq << " qp[" << _ifct << "]";
          _print_interm( _np, _Dyq, oq.str(), os );
        }
        if( options.RESRECORD )
          results_sensitivity[_ifct].push_back( Results( _t, _nx, _Dy, _np, _Dyq ) );
//        for( unsigned iy=0; lk && iy<_ny+_np; iy++ )
//          lk[_istg-1][_ifct*(_nx+_np)+iy] = iy<_ny? _Dy[iy]: _Dyq[iy-_ny];

        // Keep track of function derivatives
        for( unsigned iq=0; iq<_np; iq++ ) _Dfp[iq*_nf+_ifct] = _Dyq[iq];
      }
    }

    // Display / return function derivatives
    for( unsigned ip=0; ip<_np; ++ip )
      _fp.push_back( std::vector<double>( _Dfp.data()+ip*_nf, _Dfp.data()+(ip+1)*_nf ) );
//    for( unsigned i=0; fp && i<_nf*_np; i++ ) fp[i] = _Dfp[i];
    if( options.DISPLAY >= 1 ){
      for( unsigned iq=0; iq<_np; iq++ ){
        std::ostringstream ofp; ofp << " fp[" << iq << "]";
        _print_interm( _nf, _Dfp.data()+iq*_nf, ofp.str(), os );
      }
    }
  }
  catch(...){
    _END_SEN();
    if( options.DISPLAY >= 1 ) _print_stats( stats_sensitivity, os );
    return STATUS::FAILURE;
  }

  _END_SEN();
  if( options.DISPLAY >= 1 ) _print_stats( stats_sensitivity, os );
  return STATUS::NORMAL;
}

template <typename... ExtOps>
bool
ODESLVS_CVODES<ExtOps...>::_INI_FSA
( double const* p )
{
  // Initialize bound propagation
  if( !_INI_D_SEN( p, _np, _nq ) || !_REINI_SEN() )
    return false;

  // Set SUNDIALS sensitivity/quadrature arrays
  if( _nvec != _np ){
    if( _Ny )   N_VDestroyVectorArray( _Ny,  _nvec );
    _Ny = nullptr;
    if( _Nyq )  N_VDestroyVectorArray( _Nyq, _nvec );
    _Nyq = nullptr;
    _nvec = _np;
    _Ny  = N_VCloneVectorArray( _np, _Nx );
    if( _nq ) _Nyq = N_VCloneVectorArray( _np, _Nq );
  }
  for( unsigned i=0; i<_np; i++ ){
    if( !_Ny[i] || NV_LENGTH_S( _Ny[i] ) != (sunindextype)_nx ){
      if( _Ny[i] ) N_VDestroy( _Ny[i] );
      _Ny[i] = N_VNew_Serial( _nx, sunctx );
    }
    if( _Nyq && (!_Nyq[i] || NV_LENGTH_S( _Nyq[i] ) != (sunindextype)_nq) ){
      if( _Nyq[i] ) N_VDestroy( _Nyq[i] );
      _Nyq[i] = _nq? N_VNew_Serial( _nq, sunctx ): nullptr;
    }
  }

  // Reset result record and statistics
  results_sensitivity.clear();
  results_sensitivity.resize( _np );
  _init_stats( stats_sensitivity );

  return true;
}

template <typename... ExtOps>
inline
int
ODESLVS_CVODES<ExtOps...>::MC_CVRHSF__
( int Ns, sunrealtype t, N_Vector x, N_Vector xdot, int is, N_Vector y,
  N_Vector ydot, void* user_data, N_Vector tmp1, N_Vector tmp2 )
{
#ifdef MC__BASE_CVODES_CHECK
  //std::cout << "BASE_CVODES::PTR_BASE_CVODES: " << BASE_CVODES::PTR_BASE_CVODES
  //          << "  PTR_CVRHSF: " << PTR_CVRHSF << std::endl;
  assert( BASE_CVODES::PTR_BASE_CVODES != nullptr && PTR_CVRHSF != nullptr);
#endif
  return (BASE_CVODES::PTR_BASE_CVODES->*PTR_CVRHSF)( Ns, t, x, xdot, is, y, ydot, user_data, tmp1, tmp2 );
}

template <typename... ExtOps>
inline
int
ODESLVS_CVODES<ExtOps...>::CVRHSF__
( int Ns, sunrealtype t, N_Vector x, N_Vector xdot, int is, N_Vector y,
  N_Vector ydot, void* user_data, N_Vector tmp1, N_Vector tmp2 )
{
#ifdef MC__ODESLVS_CVODES_DEBUG
  std::cout << "@t=" << t << "\nx:\n";
  for( unsigned i=0; i<NV_LENGTH_S( x ); i++ ) std::cout << NV_Ith_S( x, i ) << std::endl;
  std::cout << "y[" << is << "]:\n";
  for( unsigned i=0; i<NV_LENGTH_S( y ); i++ ) std::cout << NV_Ith_S( y, i ) << std::endl;
#endif
  bool flag = _RHS_D_SEN( t, NV_DATA_S( x ), NV_DATA_S( y ), NV_DATA_S( ydot ), is );
#ifdef MC__ODESLVS_CVODES_DEBUG
  std::cout << "ydot[" << is << "]:\n";
  for( unsigned i=0; i<NV_LENGTH_S( ydot ); i++ ) std::cout << NV_Ith_S( ydot, i ) << std::endl;
  { int dum; std::cin >> dum; }
#endif
  stats_sensitivity.numRHS++;
  return( flag? 0: -1 );
}

template <typename... ExtOps>
inline
int
ODESLVS_CVODES<ExtOps...>::MC_CVQUADF__
( int Ns, sunrealtype t, N_Vector x, N_Vector* y, N_Vector qdot, N_Vector* qSdot, 
  void* user_data, N_Vector tmp1, N_Vector tmp2 )
{
#ifdef MC__BASE_CVODES_CHECK
  //std::cout << "BASE_CVODES::PTR_BASE_CVODES: " << BASE_CVODES::PTR_BASE_CVODES
  //          << "  PTR_CVQUADF: " << PTR_CVQUADF << std::endl;
  assert( BASE_CVODES::PTR_BASE_CVODES != nullptr && PTR_CVQUADF != nullptr);
#endif
  return (BASE_CVODES::PTR_BASE_CVODES->*PTR_CVQUADF)( Ns, t, x, y, qdot, qSdot, user_data, tmp1, tmp2 );
}

template <typename... ExtOps>
inline
int
ODESLVS_CVODES<ExtOps...>::CVQUADF__
( int Ns, sunrealtype t, N_Vector x, N_Vector* y, N_Vector qdot, N_Vector* qSdot, 
  void* user_data, N_Vector tmp1, N_Vector tmp2 )
{
#ifdef MC__ODESLVS_CVODES_DEBUG
  std::cout << "@t=" << t << "\nx:\n";
  for( unsigned i=0; i<NV_LENGTH_S( x ); i++ ) std::cout << NV_Ith_S( x, i ) << std::endl;
  std::cout << "qdot:\n";
  for( unsigned i=0; i<NV_LENGTH_S( qdot ); i++ ) std::cout << NV_Ith_S( qdot, i ) << std::endl;
#endif
  bool flag = true;
  for( int is=0; is<Ns && flag; is++ ){
    _GET_D_SEN( NV_DATA_S(x), NV_DATA_S(y[is]), (sunrealtype*)nullptr, 0, (sunrealtype*)nullptr );
#ifdef MC__ODESLVS_CVODES_DEBUG
    std::cout << "y:\n";
    for( unsigned i=0; i<NV_LENGTH_S( y[is] ); i++ ) std::cout << NV_Ith_S( y[is], i ) << std::endl;
#endif
    flag = _RHS_D_QUAD( _nq, NV_DATA_S( qSdot[is] ), is );
#ifdef MC__ODESLVS_CVODES_DEBUG
    std::cout << "qSdot:\n";
    for( unsigned i=0; i<NV_LENGTH_S( qSdot[is] ); i++ ) std::cout << NV_Ith_S( qSdot[is], i ) << std::endl;
    { int dum; std::cin >> dum; }
#endif
  }
  return( flag? 0: -1 );
}

//! @fn template <typename... ExtOps> inline typename ODESLVS_CVODES<ExtOps...>::STATUS ODESLVS_CVODES<ExtOps...>::solve_sensitivity(
//! std::vector<double> const& p, std::ostream& os=std::cout )
//!
//! This function computes a solution to the parametric ODEs with forward
//! sensitivity analysis:
//!  - <a>p</a>   [input]  parameter values
//!  - <a>os</a>  [input]  output stream [default: std::cout]
//! .
//! The return value is the status.
template <typename... ExtOps>
typename ODESLVS_CVODES<ExtOps...>::STATUS
inline
ODESLVS_CVODES<ExtOps...>::solve_sensitivity
( std::vector<double> const& p, std::ostream& os )
{
  registration();
  STATUS flag = _states_FSA( p.data(), os );
  unregistration();
  return flag;
}

//! @fn template <typename... ExtOps> inline typename ODESLVS_CVODES<ExtOps...>::STATUS ODESLVS_CVODES<ExtOps...>::solve_sensitivity(
//! std::vector<double> const& p, std::ostream& os=std::cout )
//!
//! This function computes a solution to the parametric ODEs with forward
//! sensitivity analysis:
//!  - <a>p</a>   [input]  parameter values
//!  - <a>os</a>  [input]  output stream [default: std::cout]
//! .
//! The return value is the status.
template <typename... ExtOps>
typename ODESLVS_CVODES<ExtOps...>::STATUS
inline
ODESLVS_CVODES<ExtOps...>::solve_sensitivity
( double const* p, std::ostream& os )
{
  registration();
  STATUS flag = _states_FSA( p, os );
  unregistration();
  return flag;
}

template <typename... ExtOps>
typename ODESLVS_CVODES<ExtOps...>::STATUS
inline
ODESLVS_CVODES<ExtOps...>::_states_FSA
( double const* p, std::ostream& os )
{
  // Check arguments
  if( !_np )
    return _states( nullptr, false, os );
  else if( !p )
    return STATUS::FATAL;

  try{
    // Initialize trajectory integration
    if( !ODESLV_CVODES<ExtOps...>::_INI_STA( p ) 
     || !_INI_FSA( p ) ) return STATUS::FATAL;
    _t = _dT[0];
    const unsigned NSTEP = options.RESRECORD? options.RESRECORD: 1;

    // Initial state/quadrature values
    if( !ODESLV_BASE<ExtOps...>::_IC_D_SET()
     || !ODESLV_BASE<ExtOps...>::_IC_D_STA( _t, NV_DATA_S( _Nx ) )
     || ( _Nq && !ODESLV_BASE<ExtOps...>::_IC_D_QUAD( NV_DATA_S( _Nq ) ) ) )
      { _END_STA(); _END_SEN(); return STATUS::FATAL; }
    ODESLV_BASE<ExtOps...>::_GET_D_STA( NV_DATA_S(_Nx), _nq && _Nq? NV_DATA_S(_Nq): nullptr );

    // Display / record / return initial results
    _xk.push_back( std::vector<double>( _Dx, _Dx+_nx ) );
    if( _nq ) _qk.push_back( std::vector<double>( _Dq, _Dq+_nq ) );
    if( options.DISPLAY >= 1 ){
      _print_interm( _t, _nx, _Dx, " x", os );
      _print_interm( _nq, _Dq, " q", os );
    }
//    if( options.RESRECORD )
//      results_state.push_back( Results( _t, _nx, NV_DATA_S(_Nx), _nq, _nq? NV_DATA_S(_Nq): nullptr ) );

    // Initial state/quadrature sensitivities
    _xpk.push_back( std::vector<std::vector<double>>( _np ) );
    if( _nq ) _qpk.push_back( std::vector<std::vector<double>>( _np ) );
    for( _isen=0; _isen<_np; _isen++ ){
      if( !_IC_SET_FSA( _isen )
       || !ODESLV_BASE<ExtOps...>::_IC_D_STA( _t, NV_DATA_S(_Ny[_isen]) )
       || ( _Nyq && _Nyq[_isen] && !ODESLV_BASE<ExtOps...>::_IC_D_QUAD( NV_DATA_S(_Nyq[_isen]) ) ) ) 
        { _END_STA(); _END_SEN(); return STATUS::FATAL; }
      _GET_D_SEN( NV_DATA_S(_Ny[_isen]), _nq, _nq && _Nyq? NV_DATA_S(_Nyq[_isen]): nullptr );

      // Display / record / return initial results
      _xpk[0].push_back( std::vector<double>( _Dy, _Dy+_nx ) );
      if( _nq ) _qpk[0].push_back( std::vector<double>( _Dyq, _Dyq+_nq ) );
      if( options.DISPLAY >= 1 ){
        std::ostringstream oxp; oxp << " xp[" << _isen << "]";
        _print_interm( _nx, _Dy, oxp.str(), os );
        std::ostringstream oqp; oqp << " qp[" << _isen << "]";
        _print_interm( _nq, _Dyq, oqp.str(), os );
      }
//      if( options.RESRECORD )
//        results_sensitivity[_isen].push_back( Results( _t, _nx, NV_DATA_S(_Ny[_isen]), _nq, _nq? NV_DATA_S(_Nyq[_isen]):nullptr ) );
//      for( unsigned ix=0; xpk && ix<_nx+_nq; ix++ )
//        xpk[0][(_nx+_nq)*_isen+ix] = ix<_nx? _Dy[ix]: _Dyq[ix-_nx];
    }

    // Integrate ODEs through each stage using SUNDIALS
    if( !ODESLV_CVODES<ExtOps...>::_INI_CVODE()
     || !_INI_CVODES_FSA() )
      { _END_STA(); _END_SEN(); return STATUS::FATAL; }

    for( _istg=0; _istg<_nsmax; _istg++ ){
    
      // Account for state/sensitivity discontinuities (if any) at stage times
      // and solver reinitialization (if applicable)
      _pos_ic = ( _vIC.size()>=_nsmax? _istg:0 );
      if( _pos_ic
       && ( !ODESLV_BASE<ExtOps...>::_CC_D_SET( _pos_ic )
         || !ODESLV_BASE<ExtOps...>::_CC_D_STA( _t, NV_DATA_S( _Nx ) )
         || !ODESLV_CVODES<ExtOps...>::_CC_CVODE_STA() ) )
        { _END_STA(); _END_SEN(); return STATUS::FAILURE; }
      if( _istg 
       //&& !ODESLV_CVODES<ExtOps...>::_CC_CVODE_QUAD() )
       && ( ( _Nq && !ODESLV_BASE<ExtOps...>::_IC_D_QUAD( NV_DATA_S( _Nq ) ) ) // quadrature reinitialization
         || !ODESLV_CVODES<ExtOps...>::_CC_CVODE_QUAD() ) )
        { _END_STA(); _END_SEN(); return STATUS::FAILURE; }
      if( options.RESRECORD )
        results_state.push_back( Results( _t, _nx, NV_DATA_S(_Nx), _nq, _nq? NV_DATA_S(_Nq): nullptr ) );

      for( _isen=0; _isen<_np; _isen++ ){
        if( _pos_ic
         && ( !_CC_SET_FSA( _pos_ic, _isen )
           || !_CC_D_SEN( _t, NV_DATA_S( _Nx ), NV_DATA_S(_Ny[_isen]) )
           || ( _isen==_np-1 && !_CC_CVODES_FSA() ) ) )
            { _END_STA(); _END_SEN(); return STATUS::FATAL; }
        if( _istg
         //&& ( _isen==_np-1 && !_CC_CVODES_QUAD() ) )
         && ( ( _Nyq && _Nyq[_isen] && !ODESLV_BASE<ExtOps...>::_IC_D_QUAD( NV_DATA_S(_Nyq[_isen]) ) ) //quadrature sensitivity reinitialization
           || ( _isen==_np-1 && !_CC_CVODES_QUAD() ) ) )
            { _END_STA(); _END_SEN(); return STATUS::FATAL; }
#ifdef MC__ODESLVS_CVODES_DEBUG
        for( unsigned iy=0; iy<NV_LENGTH_S(_Ny[_isen]); iy++ )
          std::cout << "_Ny" << _isen << "[" << iy << "] = " << NV_Ith_S(_Ny[_isen],iy) << std::endl;
        for( unsigned iy=0; _nq && iy<NV_LENGTH_S(_Nyq[_isen]); iy++ )
          std::cout << "_Nyq" << _isen << "[" << iy << "] = " << NV_Ith_S(_Nyq[_isen],iy) << std::endl;
#endif
        if( options.RESRECORD )
          results_sensitivity[_isen].push_back( Results( _t, _nx, NV_DATA_S(_Ny[_isen]), _nq, _nq? NV_DATA_S(_Nyq[_isen]): nullptr ) );
      }

      // update list of operations in RHS, JAC, QUAD, RHSFSA and QUADFSA
      _pos_rhs  = ( _vRHS.size()<=1?  0: _istg );
      _pos_quad = ( _vQUAD.size()<=1? 0: _istg );
      if( (!_istg || _pos_rhs || _pos_quad)
        && ( !ODESLV_BASE<ExtOps...>::_RHS_D_SET( _pos_rhs, _pos_quad )
          || !_RHS_SET_FSA( _pos_rhs, _pos_quad )
          || !_RHS_D_SET( _np, _nq ) ) )
        { _END_STA(); _END_SEN(); return STATUS::FATAL; }

      // integrate till end of time stage
      _cv_flag = CVodeSetStopTime( _cv_mem, _dT[_istg+1] );
      if( _check_cv_flag(&_cv_flag, "CVodeSetStopTime", 1) )
        { _END_STA(); _END_SEN(); return STATUS::FATAL; }

      const double TSTEP = ( _dT[_istg+1] - _t ) / NSTEP;
      double TSTOP = _t+TSTEP;
      for( unsigned k=0; k<NSTEP; k++, TSTOP+=TSTEP ){
        if( k+1 == NSTEP ) TSTOP = _dT[_istg+1];
        _cv_flag = CVode( _cv_mem, TSTOP, _Nx, &_t, CV_NORMAL );
        if( _check_cv_flag( &_cv_flag, "CVode", 1 ) )
         //|| (options.NMAX && stats_sensitivity.numSteps > options.NMAX) )
          throw Exceptions( Exceptions::INTERN );

        // intermediate record
        if( options.RESRECORD ){
          if( _nq ){
            _cv_flag = CVodeGetQuad( _cv_mem, &_t, _Nq );
            if( _check_cv_flag(&_cv_flag, "CVodeGetQuad", 1) )
              { _END_STA(); return STATUS::FATAL; }
          }
          results_state.push_back( Results( _t, _nx, NV_DATA_S(_Nx), _nq, _nq? NV_DATA_S(_Nq): nullptr ) );
          for( _isen=0; _isen<_np; _isen++ ){
            _cv_flag = CVodeGetSens1(_cv_mem, &_t, _isen, _Ny[_isen] );
            if( _check_cv_flag( &_cv_flag, "CVodeGetSens", 1) )
             { _END_STA(); _END_SEN(); return STATUS::FATAL; }
            if( _nq ){
              _cv_flag = CVodeGetQuadSens1(_cv_mem, &_t, _isen, _Nyq[_isen]);
              if( _check_cv_flag( &_cv_flag, "CVodeGetQuadSens", 1) )
                { _END_STA(); _END_SEN(); return STATUS::FATAL; }
              //for( unsigned iq=0; iq<_nq; ++iq )
              //  std::cerr << "_Nyq[" << _isen << "][" << iq << "] = " << NV_DATA_S(_Nyq[_isen])[iq] << std::endl;
            }
            results_sensitivity[_isen].push_back( Results( _t, _nx, NV_DATA_S(_Ny[_isen]), _nq, _nq? NV_DATA_S(_Nyq[_isen]): nullptr ) );
          }
        }
      }

      // Intermediate states and quadratures
      if( _nq ){
        _cv_flag = CVodeGetQuad( _cv_mem, &_t, _Nq );
        if( _check_cv_flag(&_cv_flag, "CVodeGetQuad", 1) )
          { _END_STA(); return STATUS::FATAL; }
      }
      ODESLV_BASE<ExtOps...>::_GET_D_STA( NV_DATA_S(_Nx), _nq && _Nq? NV_DATA_S(_Nq): nullptr );

      // Display / return stage results
      _xk.push_back( std::vector<double>( _Dx, _Dx+_nx ) );
      if( _nq ) _qk.push_back( std::vector<double>( _Dq, _Dq+_nq ) );
      if( options.DISPLAY >= 1 ){
        _print_interm( _t, _nx, _Dx, " x", os );
        _print_interm( _nq, _Dq, " q", os );
      }

      // Add intermediate function terms
      _pos_fct = ( _vFCT.size()>=_nsmax? _istg:0 );
      if( (_vFCT.size()>=_nsmax || _istg==_nsmax-1)
       && !ODESLV_BASE<ExtOps...>::_FCT_D_STA( _pos_fct, _t ) )
        { _END_STA(); _END_SEN(); return STATUS::FATAL; }

      // Intermediate state and quadrature sensitivities
      _xpk.push_back( std::vector<std::vector<double>>( _np ) );
      if( _nq ) _qpk.push_back( std::vector<std::vector<double>>( _np ) );
      for( _isen=0; _isen<_np; _isen++ ){
        _cv_flag = CVodeGetSens1(_cv_mem, &_t, _isen, _Ny[_isen] );
        if( _check_cv_flag( &_cv_flag, "CVodeGetSens", 1) )
         { _END_STA(); _END_SEN(); return STATUS::FATAL; }
        if( _nq ){
          _cv_flag = CVodeGetQuadSens1(_cv_mem, &_t, _isen, _Nyq[_isen]);
          if( _check_cv_flag( &_cv_flag, "CVodeGetQuadSens", 1) )
            { _END_STA(); _END_SEN(); return STATUS::FATAL; }
          //for( unsigned iq=0; iq<_nq; ++iq )
          //  std::cerr << "_Nyq[" << _isen << "][" << iq << "] = " << NV_DATA_S(_Nyq[_isen])[iq] << std::endl;
        }
        _GET_D_SEN( NV_DATA_S(_Nx), NV_DATA_S(_Ny[_isen]), _nq && _Nq? NV_DATA_S(_Nq): nullptr,
                    _nq, _nq && _Nyq[_isen]? NV_DATA_S(_Nyq[_isen]): nullptr );

        // Display / return stage results
        _xpk[_istg+1].push_back( std::vector<double>( _Dy, _Dy+_nx ) );
        if( _nq ) _qpk[_istg+1].push_back( std::vector<double>( _Dyq, _Dyq+_nq ) );
        if( options.DISPLAY >= 1 ){
          std::ostringstream oxp; oxp << " xp[" << _isen << "]";
          _print_interm( _nx, _Dy, oxp.str(), os );
          std::ostringstream oqp; oqp << " qp[" << _isen << "]";
          _print_interm( _nq, _Dyq, oqp.str(), os );
        }
//        for( unsigned ix=0; xpk && ix<_nx+_nq; ix++ )
//          xpk[_istg+1][(_nx+_nq)*_isen+ix] = ix<_nx? _Dy[ix]: _Dyq[ix-_nx];

        // Add intermediate function derivative terms
        if( (_vFCT.size()>=_nsmax || _istg==_nsmax-1)
         && !_FCT_D_SEN( _pos_fct, _isen, _t ) )
          { _END_STA(); _END_SEN(); return STATUS::FATAL; }
      }
    }

    // Display / return function values and derivatives
    _f = _Df;
    for( unsigned ip=0; ip<_np; ++ip )
      _fp.push_back( std::vector<double>( _Dfp.data()+ip*_nf, _Dfp.data()+(ip+1)*_nf ) );
//    for( unsigned i=0; f && i<_nf; i++ ) f[i] = _Df[i];
//    for( unsigned i=0; fp && i<_nf*_np; i++ ) fp[i] = _Dfp[i];
    if( options.DISPLAY >= 1 ){
      _print_interm( _nf, _Df.data(), " f", os );
      for( unsigned iq=0; iq<_np; iq++ ){
        std::ostringstream ofp; ofp << " fp[" << iq << "]";
        _print_interm( _nf, _Dfp.data()+iq*_nf, ofp.str(), os );
      }
    }
  }
  catch(...){
    _END_STA(); _END_SEN();
    long int nstp;
    _cv_flag = CVodeGetNumSteps( _cv_mem, &nstp );
    stats_state.numSteps += nstp;
    stats_sensitivity.numSteps += nstp;
    if( options.DISPLAY >= 1 ) _print_stats( stats_sensitivity, os );
    return STATUS::FAILURE;
  }

  long int nstp;
  _cv_flag = CVodeGetNumSteps( _cv_mem, &nstp );
  stats_state.numSteps += nstp;
  stats_sensitivity.numSteps += nstp;
#ifdef MC__ODESLVS_CVODES_DEBUG
  std::cout << "number of steps: " << nstp << std::endl;
#endif

  _END_STA(); _END_SEN();
  if( options.DISPLAY >= 1 ) _print_stats( stats_sensitivity, os );
  return STATUS::NORMAL;
}

} // end namescape mc


#include "slift.hpp"

namespace mc
{

//! @brief C++ class defining IVP in parametric ODEs as external DAG operations in MC++.
////////////////////////////////////////////////////////////////////////
//! mc::FFODE is a C++ class for defining an IVP in parametric ODEs as
//! external DAG operations in MC++.
////////////////////////////////////////////////////////////////////////
template<unsigned int ID, typename... ExtOps>
class FFODE
: public FFOp
{
public:
  // Constructors
  FFODE
    ()
    : FFOp( (int)EXTERN )
    {}

  // Functor
  FFVar& operator()
    ( unsigned const idep, unsigned const nVar, FFVar const* pVar, ODESLVS_CVODES<ExtOps...>*& pODE,
      bool copy=true )
    const
    {
      data = pODE;
      info = ID;
      unsigned const nFun = pODE->nf();
      FFVar** pRes = insert_external_operation( *this, nFun, nVar, pVar );
      FFOp* pOp = (*pRes)->opdef().first;
      if( copy && !pOp->owndata ){//pOp->data != pODE ){
        ODESLVS_CVODES<ExtOps...>* pODEloc = new ODESLVS_CVODES<ExtOps...>;
        pODEloc->set( *pODE );
        pODEloc->options = pODE->options;
        pODEloc->setup( *pODE );
	auto opins = update_data( *this, pOp, pODEloc, true );
	assert( opins.second );
	pOp = opins.first;
        data = pODE = static_cast<ODESLVS_CVODES<ExtOps...>*>( pOp->data );
        data = pOp->data;
      }
      assert( idep < nFun*nVar );
      return *(pRes[idep]);
    }
  FFVar** operator()
    ( unsigned const nVar, FFVar const* pVar, ODESLVS_CVODES<ExtOps...>*& pODE,
      bool copy=true )
    const
    {
      data = pODE;
      info = ID;
      unsigned const nFun = pODE->nf();
      FFVar** pRes = insert_external_operation( *this, nFun, nVar, pVar );
      FFOp* pOp = (*pRes)->opdef().first;
      if( copy && !pOp->owndata ){//pOp->data != pODE ){
        ODESLVS_CVODES<ExtOps...>* pODEloc = new ODESLVS_CVODES<ExtOps...>;
        pODEloc->set( *pODE );
        pODEloc->options = pODE->options;
        pODEloc->setup( *pODE );
	auto opins = update_data( *this, pOp, pODEloc, true );
	assert( opins.second );
	pOp = opins.first;
        pODE = static_cast<ODESLVS_CVODES<ExtOps...>*>( pOp->data );
        data = pOp->data;
      }
      return pRes;
    }

  // Evaluation overloads
  template< typename T >
  void eval
    ( unsigned const nRes, T* vRes, unsigned const nVar, T const* vVar, unsigned const* mVar )
    const
    {
      throw std::runtime_error("Error: No generic implementation for FFODE\n");
    }

  void eval
    ( unsigned const nRes, double* vRes, unsigned const nVar, double const* vVar, unsigned const* mVar )
    const;

  void eval
    ( unsigned const nRes, fadbad::F<double>* vRes, unsigned const nVar, fadbad::F<double> const* vVar,
      unsigned const* mVar )
    const;

  void eval
    ( unsigned const nRes, FFVar* vRes, unsigned const nVar, FFVar const* vVar, unsigned const* mVar )
    const;

  void eval
    ( unsigned const nRes, FFDep* vRes, unsigned const nVar, FFDep const* vVar, unsigned const* mVar )
    const;

  void eval
    ( unsigned const nRes, fadbad::F<FFVar>* vRes, unsigned const nVar, fadbad::F<FFVar> const* vVar,
      unsigned const* mVar )
    const;

  void eval
    ( unsigned const nRes, SLiftVar* vRes, unsigned const nVar, SLiftVar const* vVar, unsigned const* mVar )
    const;

  // Derivatives
  void deriv
    ( unsigned const nRes, FFVar const* vRes, unsigned const nVar, FFVar const* vVar, FFVar** vDer )
    const;

  // Properties
  std::string name
    ()
    const
    //{ return "ODE"; }
    { std::ostringstream oss; oss << data; return "ODE[" + oss.str() + "]"; }

  //! @brief Return whether or not operation is commutative
  bool commutative
    ()
    const
    { return false; }

  // Data cleanup
  bool cleanup
    ()
    const
    { 
#ifdef MC__FFODE_TRACE
      std::cout << "FFODE::cleanup\n"; 
#endif
      delete static_cast<ODESLVS_CVODES<ExtOps...>*>( data );
      data = nullptr;
      return true;
    }
};

//! @brief C++ class defining gradient of IVP in parametric ODEs as external DAG operations in MC++.
////////////////////////////////////////////////////////////////////////
//! mc::FFGRADODE is a C++ class for defining the gradient of an IVP
//! in parametric ODEs as external DAG operations in MC++.
////////////////////////////////////////////////////////////////////////
template<unsigned int ID, typename... ExtOps>
class FFGRADODE
: public FFOp
{
public:
  // Constructors
  FFGRADODE
    ()
    : FFOp( (int)EXTERN )
    {}

  // Functor
  FFVar& operator()
    ( unsigned const idep, unsigned const nVar, FFVar const* pVar, ODESLVS_CVODES<ExtOps...>*& pODE,
      bool copy=true )
    const
    {
      data = pODE;
      info = ID+1;
      unsigned const nFun = pODE->nf();
      FFVar** pRes = insert_external_operation( *this, nFun*nVar, nVar, pVar );
      FFOp* pOp = (*pRes)->opdef().first;
      if( copy && !pOp->owndata ){//pOp->data != pODE ){
        ODESLVS_CVODES<ExtOps...>* pODEloc = new ODESLVS_CVODES<ExtOps...>;
        pODEloc->set( *pODE );
        pODEloc->options = pODE->options;
        pODEloc->setup( *pODE );
	auto opins = update_data( *this, pOp, pODEloc, true );
	assert( opins.second );
	pOp = opins.first;
        pODE = static_cast<ODESLVS_CVODES<ExtOps...>*>( pOp->data );
        data = pOp->data;
      }
#ifdef MC__FFODE_CHECK
      assert( idep < nFun*nVar );
#endif
      return *(pRes[idep]);
    }
  FFVar** operator()
    ( unsigned const nVar, FFVar const* pVar, ODESLVS_CVODES<ExtOps...>*& pODE,
      bool copy=true )
    const
    {
      data = pODE;
      info = ID+1;
      unsigned const nFun = pODE->nf();
      FFVar** pRes = insert_external_operation( *this, nFun*nVar, nVar, pVar );
      FFOp* pOp = (*pRes)->opdef().first;
      if( copy && !pOp->owndata ){//pOp->data != pODE ){
        ODESLVS_CVODES<ExtOps...>* pODEloc = new ODESLVS_CVODES<ExtOps...>;
        pODEloc->set( *pODE );
        pODEloc->options = pODE->options;
        pODEloc->setup( *pODE );
	auto opins = update_data( *this, pOp, pODEloc, true );
	assert( opins.second );
	pOp = opins.first;
        pODE = static_cast<ODESLVS_CVODES<ExtOps...>*>( pOp->data );
        data = pOp->data;
      }
      return pRes;
    }

  // Evaluation overloads
  template< typename T >
  void eval
    ( unsigned const nRes, T* vRes, unsigned const nVar, T const* vVar, unsigned const* mVar )
    const
    {
      throw std::runtime_error("Error: No generic implementation for FFGRADODE\n");
    }

  void eval
    ( unsigned const nRes, double* vRes, unsigned const nVar, double const* vVar, unsigned const* mVar )
    const;

  void eval
    ( unsigned const nRes, FFVar* vRes, unsigned const nVar, FFVar const* vVar, unsigned const* mVar )
    const;

  void eval
    ( unsigned const nRes, FFDep* vRes, unsigned const nVar, FFDep const* vVar, unsigned const* mVar )
    const;

  // Properties
  std::string name
    ()
    const
    //{ return "GRADODE"; }
    { std::ostringstream oss; oss << data; return "GRADODE[" + oss.str() + "]"; }

  //! @brief Return whether or not operation is commutative
  bool commutative
    ()
    const
    { return false; }

  // Data cleanup
  bool cleanup
    ()
    const
    {
#ifdef MC__FFODE_TRACE
      std::cout << "FFGRADODE::cleanup\n"; 
#endif
      delete static_cast<ODESLVS_CVODES<ExtOps...>*>( data );
      data = nullptr;
      return true;
    }
};

template<unsigned int ID, typename... ExtOps>
inline void
FFODE<ID,ExtOps...>::eval
( unsigned const nRes, double* vRes, unsigned const nVar, double const* vVar,
  unsigned const* mVar )
const
{
#ifdef MC__FFODE_TRACE
  std::cout << "FFODE::eval: double\n";
#endif
  ODESLVS_CVODES<ExtOps...>* pODE = static_cast<ODESLVS_CVODES<ExtOps...>*>( data );
#ifdef MC__FFODE_CHECK
  assert( pODE && nRes == pODE->nf() && nVar == pODE->np() );
#endif

  if( pODE->solve_state( vVar ) != ODESLVS_CVODES<ExtOps...>::NORMAL ){
    for( unsigned i=0; i<nVar; ++i ) std::cout << "vVar[" << i << "] = " << vVar[i] << std::endl;
    std::cerr << "FFODE::eval double ** Integration failure\n";
    throw typename ODESLVS_CVODES<ExtOps...>::Exceptions( ODESLVS_CVODES<ExtOps...>::Exceptions::INTERN );
  }
  for( unsigned i=0; i<nRes; ++i ) vRes[i] = pODE->val_function()[i];
}

template<unsigned int ID, typename... ExtOps>
inline void
FFODE<ID,ExtOps...>::eval
( unsigned const nRes, FFVar* vRes, unsigned const nVar, FFVar const* vVar,
  unsigned const* mVar )
const
{
#ifdef MC__FFODE_TRACE
  std::cout << "FFODE::eval: FFVar\n";
#endif
  ODESLVS_CVODES<ExtOps...>* pODE = static_cast<ODESLVS_CVODES<ExtOps...>*>( data );
#ifdef MC__FFODE_CHECK
  assert( pODE && nRes == pODE->nf() && nVar == pODE->np() );
#endif

  FFVar** pRes = operator()( nVar, vVar, pODE );
  for( unsigned j=0; j<nRes; ++j ) vRes[j] = *(pRes[j]);
}

template<unsigned int ID, typename... ExtOps>
inline void
FFODE<ID,ExtOps...>::eval
( unsigned const nRes, FFDep* vRes, unsigned const nVar, FFDep const* vVar,
  unsigned const* mVar )
const
{
#ifdef MC__FFODE_TRACE
  std::cout << "FFODE::eval: FFDep\n";
#endif
#ifdef MC__FFODE_CHECK
  ODESLVS_CVODES<ExtOps...>* pODE = static_cast<ODESLVS_CVODES<ExtOps...>*>( data );
  assert( pODE && nRes == pODE->nf() && nVar == pODE->np() );
#endif

  vRes[0] = 0;
  for( unsigned i=0; i<nVar; ++i ) vRes[0] += vVar[i];
  vRes[0].update( FFDep::TYPE::N );
  for( unsigned j=1; j<nRes; ++j ) vRes[j] = vRes[0];
}

template<unsigned int ID, typename... ExtOps>
inline void
FFODE<ID,ExtOps...>::eval
( unsigned const nRes, fadbad::F<double>* vRes, unsigned const nVar, fadbad::F<double> const* vVar,
  unsigned const* mVar )
const
{
#ifdef MC__FFODE_TRACE
  std::cout << "FFODE::eval: fadbad::F<double>\n";
#endif
  ODESLVS_CVODES<ExtOps...>* pODE = static_cast<ODESLVS_CVODES<ExtOps...>*>( data );
#ifdef MC__FFODE_CHECK
  assert( pODE && nRes == pODE->nf() && nVar == pODE->np() );
#endif

  std::vector<double> vVarVal( nVar );
  for( unsigned i=0; i<nVar; ++i ) vVarVal[i] = vVar[i].val();
  if( nVar <= 3*nRes ){
    if( pODE->solve_sensitivity( vVarVal ) != ODESLVS_CVODES<ExtOps...>::NORMAL ){
      std::cerr << "FFODE::eval fadbad::F<double> ** Integration failure\n";
      throw typename ODESLVS_CVODES<ExtOps...>::Exceptions( ODESLVS_CVODES<ExtOps...>::Exceptions::INTERN );
    }
  }
  else{
    if( pODE->solve_adjoint( vVarVal ) != ODESLVS_CVODES<ExtOps...>::NORMAL ){
      std::cerr << "FFODE::eval ** Integration failure\n";
      throw typename ODESLVS_CVODES<ExtOps...>::Exceptions( ODESLVS_CVODES<ExtOps...>::Exceptions::INTERN );
    }
  }

  for( unsigned k=0; k<nRes; ++k ){
    vRes[k] = pODE->val_function()[k];
    for( unsigned i=0; i<nVar; ++i )
      vRes[k].setDepend( vVar[i] );
    for( unsigned j=0; j<vRes[k].size(); ++j ){
      vRes[k][j] = 0.;
      for( unsigned i=0; i<nVar; ++i ){
        if( vVar[i][j] == 0. ) continue;
        vRes[k][j] += pODE->val_function_gradient()[i][k] * vVar[i][j];
      }
    }
  }
}

template<unsigned int ID, typename... ExtOps>
inline void
FFODE<ID,ExtOps...>::eval
( unsigned const nRes, fadbad::F<FFVar>* vRes, unsigned const nVar, fadbad::F<FFVar> const* vVar,
  unsigned const* mVar )
const
{
#ifdef MC__FFODE_TRACE
  std::cout << "FFODE::eval: fadbad::F<FFVar>\n";
#endif
  ODESLVS_CVODES<ExtOps...>* pODE = static_cast<ODESLVS_CVODES<ExtOps...>*>( data );
#ifdef MC__FFODE_CHECK
  assert( pODE && nRes == pODE->nf() && nVar == pODE->np() );
#endif

  std::vector<FFVar> vVarVal( nVar );
  for( unsigned i=0; i<nVar; ++i )
    vVarVal[i] = vVar[i].val();
  FFVar const*const* vResVal = operator()( nVar, vVarVal.data(), pODE );
  FFGRADODE<ID,ExtOps...> ResDer;
  FFVar const*const* vResDer = ResDer( nVar, vVarVal.data(), pODE, false ); // no copy of ODE data

  for( unsigned k=0; k<nRes; ++k ){
    vRes[k] = *vResVal[k];
    for( unsigned i=0; i<nVar; ++i )
      vRes[k].setDepend( vVar[i] );
    for( unsigned j=0; j<vRes[k].size(); ++j ){
      vRes[k][j] = 0.;
      for( unsigned i=0; i<nVar; ++i ){
        if( vVar[i][j].cst() && vVar[i][j].num().val() == 0. ) continue;
        vRes[k][j] += *vResDer[k+nRes*i] * vVar[i][j];
      }
    }
  }
}

template<unsigned int ID, typename... ExtOps>
inline void
FFODE<ID,ExtOps...>::deriv
( unsigned const nRes, FFVar const* vRes, unsigned const nVar, FFVar const* vVar, FFVar** vDer )
const
{
#ifdef MC__FFODE_TRACE
  std::cout << "FFODE::eval: fadbad::F<FFVar>\n";
#endif
  ODESLVS_CVODES<ExtOps...>* pODE = static_cast<ODESLVS_CVODES<ExtOps...>*>( data );
#ifdef MC__FFODE_CHECK
  assert( pODE && nRes == pODE->nf() && nVar == pODE->np() );
#endif

  FFGRADODE<ID,ExtOps...> ResDer;
  FFVar const*const* vResDer = ResDer( nVar, vVar, pODE, false ); // no copy of ODE data
  for( unsigned k=0; k<nRes; ++k )
    for( unsigned i=0; i<nVar; ++i )
      vDer[k][i] = *vResDer[k+nRes*i];
}

template<unsigned int ID, typename... ExtOps>
inline void
FFODE<ID,ExtOps...>::eval
( unsigned const nRes, SLiftVar* vRes, unsigned const nVar, SLiftVar const* vVar,
  unsigned const* mVar )
const
{
#ifdef MC__FFODE_TRACE
  std::cout << "FFODE::eval: fadbad::F<FFVar>\n";
#endif
  ODESLVS_CVODES<ExtOps...>* pODE = static_cast<ODESLVS_CVODES<ExtOps...>*>( data );
#ifdef MC__FFODE_CHECK
  assert( pODE && nRes == pODE->nf() && nVar == pODE->np() );
#endif

  vVar->env()->lift( nRes, vRes, nVar, vVar );
}

template<unsigned int ID, typename... ExtOps>
inline void
FFGRADODE<ID,ExtOps...>::eval
( unsigned const nRes, double* vRes, unsigned const nVar, double const* vVar,
  unsigned const* mVar )
const
{
#ifdef MC__FFODE_TRACE
  std::cout << "FFGRADODE::eval: double\n";
#endif
  ODESLVS_CVODES<ExtOps...>* pODE = static_cast<ODESLVS_CVODES<ExtOps...>*>( data );
#ifdef MC__FFODE_CHECK
  assert( pODE && nRes == pODE->nf()*pODE->np() && nVar == pODE->np() );
#endif

  if( nVar <= 3*nRes ){
    if( pODE->solve_sensitivity( vVar ) != ODESLVS_CVODES<ExtOps...>::NORMAL ){
      std::cerr << "FFGRADODE::eval ** Integration failure\n";
      throw typename ODESLVS_CVODES<ExtOps...>::Exceptions( ODESLVS_CVODES<ExtOps...>::Exceptions::INTERN );
    }
  }
  else{
    if( pODE->solve_adjoint( vVar ) != ODESLVS_CVODES<ExtOps...>::NORMAL ){
      std::cerr << "FFGRADODE::eval ** Integration failure\n";
      throw typename ODESLVS_CVODES<ExtOps...>::Exceptions( ODESLVS_CVODES<ExtOps...>::Exceptions::INTERN );
    }
  }
  for( unsigned i=0, k=0; i<pODE->np(); ++i )
    for( unsigned j=0; j<pODE->nf(); ++j, ++k )
      vRes[k] = pODE->val_function_gradient()[i][j]; // DOUBLE-CHECK ENTRY ORDERING!!!
}

template<unsigned int ID, typename... ExtOps>
inline void
FFGRADODE<ID,ExtOps...>::eval
( unsigned const nRes, FFVar* vRes, unsigned const nVar, FFVar const* vVar,
  unsigned const* mVar )
const
{
#ifdef MC__FFODE_TRACE
  std::cout << "FFGRADODE::eval: FFVar\n";
#endif
  ODESLVS_CVODES<ExtOps...>* pODE = static_cast<ODESLVS_CVODES<ExtOps...>*>( data );
#ifdef MC__FFODE_CHECK
  assert( pODE && nRes == pODE->nf()*pODE->np() && nVar == pODE->np() );
#endif

  FFVar** pRes = operator()( nVar, vVar, pODE );
  for( unsigned j=0; j<nRes; ++j ) vRes[j] = *(pRes[j]);
}

template<unsigned int ID, typename... ExtOps>
inline void
FFGRADODE<ID,ExtOps...>::eval
( unsigned const nRes, FFDep* vRes, unsigned const nVar, FFDep const* vVar,
  unsigned const* mVar )
const
{
#ifdef MC__FFODE_TRACE
  std::cout << "FFGRADODE::eval: FFDep\n";
#endif
#ifdef MC__FFODE_CHECK
  ODESLVS_CVODES<ExtOps...>* pODE = static_cast<ODESLVS_CVODES<ExtOps...>*>( data );
  assert( pODE && nRes == pODE->nf()*pODE->np() && nVar == pODE->np() );
#endif

  vRes[0] = 0;
  for( unsigned i=0; i<nVar; ++i ) vRes[0] += vVar[i];
  vRes[0].update( FFDep::TYPE::N );
  for( unsigned j=1; j<nRes; ++j ) vRes[j] = vRes[0];
}

} // end namescape mc

#endif

