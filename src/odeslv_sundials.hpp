// Copyright (C) 2019 Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

#ifndef MC__ODESLV_SUNDIALS_HPP
#define MC__ODESLV_SUNDIALS_HPP

#undef  MC__ODESLV_SUNDIALS_DEBUG

#include "base_sundials.hpp"
#include "odeslv_base.hpp"

namespace mc
{
//! @brief C++ class computing solutions of parametric ODEs using SUNDIALS and MC++.
////////////////////////////////////////////////////////////////////////
//! mc::ODESLV_SUNDIALS is a C++ class for solution of IVPs in ODEs
//! using the code CVODES in SUNDIALS and MC++.
////////////////////////////////////////////////////////////////////////
class ODESLV_SUNDIALS:
  public virtual BASE_SUNDIALS,
  public virtual ODESLV_BASE,
  public virtual BASE_DE
{
  typedef int (*CVRhsFn)( realtype t, N_Vector y, N_Vector ydot, void *user_data );
  typedef int (*CVDlsDenseJacFn)( long int N, realtype t, N_Vector y, N_Vector fy,
     DlsMat Jac, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3 );

 protected:
  //! @brief Pointer to the CVODE memory block
  void *_cv_mem;

  //! @brief Return flag for SUNDIALS methods
  int _cv_flag;

  //! @brief N_Vector object holding current states
  N_Vector _Nx;

  //! @brief N_Vector object holding current state quadratures
  N_Vector _Nq;

  //! @brief stepsize
  double _h;

  //! @brief position in _vIC
  unsigned _pos_ic;

  //! @brief position in _vRHS
  unsigned _pos_rhs;

  //! @brief position in _vQUAD
  unsigned _pos_quad;
  
  //! @brief position in _vFCT
  unsigned _pos_fct;

  //! @brief checkpoints for adjoint integration
  int _nchk;

  //! @brief state parameterizations at time stages for adjoint integration
  std::vector< std::vector<realtype> > _vec_sta;

  //! @brief static pointer to class
  static ODESLV_SUNDIALS *pODESLV;

public:
  /** @ingroup ODESLV
   *  @{
   */
  //! @brief Default constructor
  ODESLV_SUNDIALS
    ();

  //! @brief Virtual destructor
  virtual ~ODESLV_SUNDIALS
    ();

  //! @brief Integrator options
  struct Options:
    public BASE_SUNDIALS::Options
  {
    //! @brief Constructor
    Options():
      BASE_SUNDIALS::Options(), DISPLAY(1), RESRECORD(false)
      { JACAPPROX = CV_DENSE; }
    //! @brief Assignment operator
    template <typename OPT> Options& operator=
      ( OPT&options ){
        BASE_SUNDIALS::Options::operator=(options);
        DISPLAY   = options.DISPLAY;
        RESRECORD = options.RESRECORD;
        return *this;
      }
    //! @brief Display level (default: 1)
    int DISPLAY;
    //! @brief Whether or not to record results (default: false)
    unsigned RESRECORD;
  } options;

  //! @brief Structure for setting up storing the solver exceptions
  class Exceptions
  {
  public:
    //! @brief Enumeration type for ODESLV_SUNDIALS exception handling
    enum TYPE{
      INTERN=-33	//!< Internal error
    };
    //! @brief Constructor for error <a>ierr</a>
    Exceptions( TYPE ierr ) : _ierr( ierr ){}
    //! @brief Inline function returning the error flag
    int ierr(){ return _ierr; }
    //! @brief Inline function returning the error description
    std::string what(){
      switch( _ierr ){
      case INTERN: default:
        return "ODESLV_SUNDIALS::Exceptions  Internal error";
       }
    }
  private:
    TYPE _ierr;
  };

  //! @brief Vector storing results (upon request only)
  std::vector< Results > results_sta;

  //! @brief Statistics for state integration
  Stats stats_sta;

  //! @brief Computes solution of parametric ODEs
  STATUS states
    ( const double*p, double**xk=0, double*f=0, std::ostream&os=std::cout );

  //! @brief Record results in file <a>ores</a>, with accuracy of <a>iprec</a> digits
  void record
    ( std::ofstream&ores, const unsigned iprec=5 ) const
    { return _record( ores, results_sta, iprec ); }
  /** @} */

protected:

  //! @brief Function to initialize CVode memory block
  virtual bool _INI_CVODE
    ( CVRhsFn MC_CVRHS, CVRhsFn MC_CVQUAD, CVDlsDenseJacFn MC_CVJAC );

  //! @brief Function to reinitialize CVode memory block
  bool _CC_CVODE_STA
    ();

  //! @brief Function to reinitialize CVode memory block
  bool _CC_CVODE_QUAD
    ();

  //! @brief Function to finalize state integration
  void _END_STA
    ();

  //! @brief Function to initialize state integration
  bool _INI_STA
    ( const double*p );

  //! @brief Function to initialize state integration
  bool _INI_STA
    ();

  //! @brief Static wrapper to function to calculate the ODEs RHS derivatives
  static int MC_CVRHSD__
    ( realtype t, N_Vector Nx, N_Vector Nxdot, void *user_data );

  //! @brief Static wrapper to function to calculate the quadrature RHS derivatives
  static int MC_CVQUADD__
    ( realtype t, N_Vector Nx, N_Vector Nqdot, void *user_data );

  //! @brief Static wrapper to function to calculate the ODEs RHS Jacobian
  static int MC_CVJACD__
    ( long int N, realtype t, N_Vector y, N_Vector fy, DlsMat Jac, void *user_data,
      N_Vector tmp1, N_Vector tmp2, N_Vector tmp3 );

  //! @brief Solve parametric ODEs forward in time through every time stages
  virtual STATUS _states
    ( const double*p, double**xk, double*f, const bool store, std::ostream&os );

  //! @brief Private methods to block default compiler methods
  ODESLV_SUNDIALS(const ODESLV_SUNDIALS&);
  ODESLV_SUNDIALS& operator=(const ODESLV_SUNDIALS&);
};

} // end namescape mc

#endif

