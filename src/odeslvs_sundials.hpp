// Copyright (C) 2019 Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

#ifndef MC__ODESLVS_SUNDIALS_HPP
#define MC__ODESLVS_SUNDIALS_HPP

#undef  MC__ODESLVS_SUNDIALS_DEBUG

#include <sstream>
#include "odeslvs_base.hpp"
#include "odeslv_sundials.hpp"

#define MC__ODESLVS_SUNDIALS_USE_BAD
#undef  MC__ODESLVS_SUNDIALS_DEBUG

namespace mc
{
//! @brief C++ class computing solutions of parametric ODEs with forward/adjoint sensitivity analysis capability using SUNDIALS and MC++.
////////////////////////////////////////////////////////////////////////
//! mc::ODESLV_SUNDIALS is a C++ class for solution of IVPs in ODEs
//! with forward/adjoint sensitivity analysis capability using the code
//! CVODES in SUNDIALS and MC++.
////////////////////////////////////////////////////////////////////////
class ODESLVS_SUNDIALS:
  public virtual ODESLV_SUNDIALS,
  public virtual BASE_DE,
  public virtual ODESLVS_BASE
{
 protected:
  typedef int (*CVRhsFn)( realtype t, N_Vector y, N_Vector ydot, void *user_data );
  typedef int (*CVADJRhsFn)( realtype t, N_Vector x, N_Vector y, N_Vector ydot, void *user_data );
  typedef int (*CVSENRhs1Fn)( int Ns, realtype t, N_Vector x, N_Vector xdot, int is,
    N_Vector y, N_Vector ydot, void *user_data, N_Vector tmp1, N_Vector tmp2 );
  typedef int (*CVSENQuadFn)( int Ns, realtype t, N_Vector x, N_Vector *y, N_Vector qdot,
    N_Vector *qSdot, void *user_data, N_Vector tmp1, N_Vector tmp2 );
  typedef int (*CVDlsDenseJacFn)( long int N, realtype t, N_Vector y, N_Vector fy,
     DlsMat Jac, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3 );
  typedef int (*CVDlsDenseJacFnB)( long int NeqB, realtype t, N_Vector y, N_Vector yB,
     N_Vector fyB, DlsMat JacB, void *user_dataB, N_Vector tmp1B, N_Vector tmp2B,
     N_Vector tmp3B );

 protected:
  //! @brief current function index
  unsigned _ifct;

  //! @brief current parameter sensitivity index
  unsigned _isen;

  //! @brief size of N_vector arrays
  unsigned _nvec;

  //! @brief N_Vector object holding current adjoints
  N_Vector *_Ny;

  //! @brief N_Vector object holding current quadratures
  N_Vector *_Nyq;
  
  //! @brief pointer to array holding identifiers of the backward problems
  int* _indexB;
  
  //! @brief pointer to array holding identifiers of the backward problems
  unsigned* _iusrB;

  //! @brief position in _Ny
  unsigned _pos_adj;

  //! @brief position in _Nq
  unsigned _pos_adjquad; 

  //! @brief static pointer to class
  static ODESLVS_SUNDIALS *pODESLVS;

 public:
  /** @ingroup ODESLV
   *  @{
   */
  //! @brief Default constructor
  ODESLVS_SUNDIALS
    ();

  //! @brief Virtual destructor
  virtual ~ODESLVS_SUNDIALS
    ();

  //! @brief Statistics for sensitivity/adjoint integration
  Stats stats_sen;

  //! @brief Vector storing adjoint/sensitivity trajectyories (see Options::RESRECORD)
  std::vector< std::vector< Results > > results_sen;

 //! @brief Propagate states and state-sensitivities forward in time through every time stages
  STATUS states_FSA
    ( const double*p, double**xk=0, double*f=0, double**xpk=0, double*fp=0, std::ostream&os=std::cout );

  //! @brief Propagate states and adjoints forward and backward in time through every time stages
  STATUS states_ASA
    ( const double*p, double**xk=0, double*f=0, double**lk=0, double*fp=0, std::ostream&os=std::cout );

  //! @brief Record state and sensitivity trajectories in files <a>obndsta</a> and <a>obndsa</a>, with accuracy of <a>iprec</a> digits
  void record
    ( std::ofstream&obndsta, std::ofstream*obndsen, const unsigned iprec=5 ) const
    { this->ODESLV_SUNDIALS::record( obndsta, iprec );
      for( unsigned isen=0; isen<results_sen.size(); ++isen )
        this->ODESLV_BASE::_record( obndsen[isen], results_sen[isen], iprec ); }

  //! @brief Record state trajectories in files <a>obndsta</a>, with accuracy of <a>iprec</a> digits
  void record
    ( std::ofstream&obndsta, const unsigned iprec=5 ) const
    { ODESLV_SUNDIALS::record( obndsta, iprec ); }
  /** @} */

 private:
  //! @brief Function to initialize CVodes memory block (virtual)
  virtual bool _INI_CVODE
    ( CVRhsFn MC_CVRHS, CVRhsFn MC_CVQUAD, CVDlsDenseJacFn MC_CVJAC );

  //! @brief Function to initialize CVodeS memory block for forward sensitivity
  bool _INI_CVODES
    ( CVSENRhs1Fn MC_CVSENRHS, CVSENQuadFn MC_CVSENQUAD );

  //! @brief Function to initialize CVodeS memory block for adjoint sensitivity
  bool _INI_CVODES
    ( CVADJRhsFn MC_CVADJRHS, CVADJRhsFn MC_CVADJQUAD, CVDlsDenseJacFnB MC_CVADJJAC,
      const unsigned ifct, int&indexB, unsigned&iusrB );

  //! @brief Function to reinitialize CVodeS memory block for forward sensitivity
  bool _CC_CVODES_FSA
    ();

  //! @brief Function to reinitialize CVodeS memory block for forward quadrature sensitivity
  bool _CC_CVODES_QUAD
    ();

  //! @brief Function to reinitialize CVodeS memory block for adjoint sensitivity
  bool _CC_CVODES_ASA
    ( const unsigned ifct, const int indexB );

  //! @brief Function to finalize sensitivity/adjoint bounding
  void _END_SEN
    ();

  //! @brief Function to initialize adjoint sensitivity analysis
  bool _INI_ASA
    ( const double *p );

  //! @brief Static wrapper to function computing the adjoint RHS
  static int MC_CVASARHSD__
    ( realtype t, N_Vector x, N_Vector y, N_Vector ydot, void *user_data );

  //! @brief Static wrapper to function computing the adjoint quadrature RHS
  static int MC_CVASAQUADD__
    ( realtype t, N_Vector x, N_Vector y, N_Vector qdot, void *user_data );

  //! @brief Static wrapper to function to calculate the adjoint RHS Jacobian
  static int MC_CVASAJACD__
    ( long int NeqB, realtype t, N_Vector y, N_Vector yB, N_Vector fyB, DlsMat JacB,
      void *user_dataB, N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B );

  //! @brief Function to initialize forward sensitivity analysis
  bool _INI_FSA
    ( const double *p );

  //! @brief Static wrapper to function to calculate the sensitivity RHS
  static int MC_CVFSARHSD__
    ( int Ns, realtype t, N_Vector x, N_Vector xdot, int is, N_Vector y,
      N_Vector ydot, void *user_data, N_Vector tmp1, N_Vector tmp2 );

  //! @brief Static wrapper to function computing the sensitivity quadrature RHS
  static int MC_CVFSAQUADD__
    ( int Ns, realtype t, N_Vector x, N_Vector *y, N_Vector qdot, N_Vector *qSdot, 
      void *user_data, N_Vector tmp1, N_Vector tmp2 );

  //! @brief Static wrapper to function to calculate the sensitivity RHS Jacobian
  static int MC_CVFSAJACD__
    ( long int NeqB, realtype t, N_Vector y, N_Vector yB, N_Vector fyB, DlsMat JacB,
      void *user_dataB, N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B );

  //! @brief Private methods to block default compiler methods
  ODESLVS_SUNDIALS(const ODESLVS_SUNDIALS&);
  ODESLVS_SUNDIALS& operator=(const ODESLVS_SUNDIALS&);
};

} // end namescape mc

#endif

