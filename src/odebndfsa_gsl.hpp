// Copyright (C) 2015 Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

#ifndef MC__ODEBNDS_GSL_HPP
#define MC__ODEBNDS_GSL_HPP

#undef  MC__ODEBNDS_GSL_DINEQI_DEBUG
#undef  MC__ODEBNDS_GSL_DINEQPM_DEBUG

#include "odebnds_base.hpp"
#include "base_gsl.hpp"
#include "mesh_gsl.hpp"

#include "tmodel.hpp"
#include "cmodel.hpp"

namespace mc
{
//! @brief C++ class computing enclosures of the reachable set of parametric ODEs using continuous-time set-valued integration with sensitivity analysis.
////////////////////////////////////////////////////////////////////////
//! mc::ODEBNDS_GSL is a C++ class that computes enclosures of the
//! reachable set of parametric ordinary differential equations
//! (ODEs) using continuous-time set-valued integration as well as
//! enclosures for the first-order derivative information (both forward
//! and adjoint sensitivity analysis).
//! It implements the methods of differential inequalities, whereby
//! polynomial models with interval or ellipsoidal remainders are used
//! to enable high-order convergence. The use of ellipsoidal remainders
//! enables stability of the enclosures for asymptotically stable ODE
//! systems when the parameter host is sufficiently small. The numerical
//! integrator is gsl_odeiv2 in GSL.
////////////////////////////////////////////////////////////////////////
template <typename T, typename PMT=mc::TModel<T>, typename PVT=mc::TVar<T> >
class ODEBNDS_GSL:
  public virtual ODEBNDS_BASE<T,PMT,PVT>,
  public virtual ODEBND_GSL<T,PMT,PVT>
{
 typedef Ellipsoid E;
 typedef BASE_DE::STATUS STATUS;
 using ODEBND_BASE<T,PMT,PVT>::NORMAL;
 using ODEBND_BASE<T,PMT,PVT>::FAILURE;
 using ODEBND_BASE<T,PMT,PVT>::FATAL;

 using ODEBND_BASE<T,PMT,PVT>::_nx;
 using ODEBND_BASE<T,PMT,PVT>::_np;
 using ODEBND_BASE<T,PMT,PVT>::_nq;
 using ODEBND_BASE<T,PMT,PVT>::_nf;
 using ODEBND_BASE<T,PMT,PVT>::_vRHS;
 using ODEBND_BASE<T,PMT,PVT>::_vQUAD;
 using ODEBND_BASE<T,PMT,PVT>::_vIC;
 using ODEBND_BASE<T,PMT,PVT>::_vFCT;

 using ODEBND_BASE<T,PMT,PVT>::_Q;
 using ODEBND_BASE<T,PMT,PVT>::_Er;
 using ODEBND_BASE<T,PMT,PVT>::_Ir;
 using ODEBND_BASE<T,PMT,PVT>::_pref;
 using ODEBND_BASE<T,PMT,PVT>::_Ip;
 using ODEBND_BASE<T,PMT,PVT>::_B;
 using ODEBND_BASE<T,PMT,PVT>::_xref;
 using ODEBND_BASE<T,PMT,PVT>::_Ix;
 using ODEBND_BASE<T,PMT,PVT>::_Iq;
 using ODEBND_BASE<T,PMT,PVT>::_vec2I;
 using ODEBND_BASE<T,PMT,PVT>::_vec2E;
 using ODEBND_BASE<T,PMT,PVT>::_IC_I_STA;
 using ODEBND_BASE<T,PMT,PVT>::_IC_I_QUAD;
 using ODEBND_BASE<T,PMT,PVT>::_CC_I_STA;
 using ODEBND_BASE<T,PMT,PVT>::_SET_I_STA;
 using ODEBND_BASE<T,PMT,PVT>::_RHS_I_STA;
 using ODEBND_BASE<T,PMT,PVT>::_RHS_I_QUAD;
 using ODEBND_BASE<T,PMT,PVT>::_JAC_I_STA;
 using ODEBND_BASE<T,PMT,PVT>::_FCT_I_STA;

 using ODEBND_BASE<T,PMT,PVT>::_PMenv;
 using ODEBND_BASE<T,PMT,PVT>::_PMx;
 using ODEBND_BASE<T,PMT,PVT>::_PMp;
 using ODEBND_BASE<T,PMT,PVT>::_PMq;
 using ODEBND_BASE<T,PMT,PVT>::_vec2PMI;
 using ODEBND_BASE<T,PMT,PVT>::_vec2PME;
 using ODEBND_BASE<T,PMT,PVT>::_IC_PM_STA;
 using ODEBND_BASE<T,PMT,PVT>::_IC_PM_QUAD;
 using ODEBND_BASE<T,PMT,PVT>::_CC_PM_STA;
 using ODEBND_BASE<T,PMT,PVT>::_SET_PM_STA;
 using ODEBND_BASE<T,PMT,PVT>::_RHS_PM_STA;
 using ODEBND_BASE<T,PMT,PVT>::_RHS_PM_QUAD;
 using ODEBND_BASE<T,PMT,PVT>::_JAC_PM_STA;
 using ODEBND_BASE<T,PMT,PVT>::_FCT_PM_STA;

 using ODEBND_BASE<T,PMT,PVT>::_diam;
 using ODEBND_BASE<T,PMT,PVT>::_dH;
 using ODEBND_BASE<T,PMT,PVT>::_hausdorff;
 using ODEBND_BASE<T,PMT,PVT>::_remainders;
 using ODEBND_BASE<T,PMT,PVT>::_print_interm;

 using ODEBND_GSL<T,PMT,PVT>::_bounds;

 template <typename U, typename PMU, typename PVU>
 friend int MC_GSLRHSFSAPM__
   ( double t, const double* x, double* xdot, void* user_data );
 template <typename U, typename PMU, typename PVU>
 friend int MC_GSLJACFSAPM__
   ( double t, const double* x, double* jac, double* xdot, void* user_data );

 private:
  //! @brief GSL data type for bounding ODE sensitivity system
  gsl_odeiv2_system _sys_sa;

  //! @brief GSL driver for bounding ODE sensitivity system
  gsl_odeiv2_driver *_driver_sa;

 protected:
  //! @brief full GSL array (states/quadratures + sensitivities)
  double *_vec_sa;

  //! @brief pointer to quadratures in full GSL array
  double *_vec_quads;

  //! @brief offset to quadratures in full GSL array
  unsigned _offset_quads;

  //! @brief stepsize
  double _h;

  //! @brief Mesh storing state bound parameterizations 
  MESH_GSL _mesh_fsa;

  //! @brief position in _vIC
  unsigned _pos_ic;

  //! @brief position in _vRHS
  unsigned _pos_rhs;

  //! @brief position in _vQUAD
  unsigned _pos_quad;
  
  //! @brief position in _vFCT
  unsigned _pos_fct;

  //! @brief Static const member for mesh type
  static const int PMOFFSET = 10;

  //! @brief static pointer to class
  static ODEBND_GSL<T,PMT,PVT> *_pODEBND;

public:
  /** @defgroup ODEBND_GSL Continuous-time set-valued integration of parametric ODEs
   *  @{
   */
  //! @brief Default constructor
  ODEBNDS_GSL();

  //! @brief Virtual destructor
  virtual ~ODEBNDS_GSL();

  //! @brief Integrator options
  struct Options: public ODEBND_GSL<T,PMT,PVT>::Options
  {
    //! @brief Constructor
    Options():
      ODEBND_GSL<T,PMT,PVT>::Options()
      {}
    //! @brief Assignment operator
    template <typename U> Options& operator=
      ( U&options ){
        ODEBND_GSL<T,PMT,PVT>::Options::operator=(options);
        return *this;
      }
  } options;

  //! @brief Structure for setting up storing the solver exceptions
  class Exceptions
  {
  public:
    //! @brief Enumeration type for ODEBND_GSL exception handling
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
        return "ODEBNDS_GSL::Exceptions  Internal error";
       }
    }
  private:
    TYPE _ierr;
  };

  //! @brief Vector storing interval sensitivity/adjoint bounds (see Options::RESRECORD)
  std::vector< Results > results_sa;

  //! @brief Statistics for state sensitivity/adjoint bounds integration
  Stats stats_sa;

  //! @brief Computes polynomial model enclosure of reachable set of parametric ODEs
  STATUS bounds
    ( const unsigned ns, const double*tk, const PVT*PMp, PVT**PMxk=0,
      PVT*PMq=0, PVT*PMf=0, std::ostream&os=std::cout );
    { ODEBND_GSL<T,PMT,PVT>::options = options;
      return ODEBND_GSL<T,PMT,PVT>::_bounds( ns, tk, PMp, PMxk, PMq, PMf, false, os); }

  //! @brief Propagate state and adjoint interval bounds forward and backward in time through every time stages
  STATUS bounds_FSA
    ( const unsigned ns, const double*tk, const PVT*PMp, PVT**PMxk,
      PVT*PMq, PVT*PMf, PVT**PMxpk, PVT*PMqp, PVT*PMfp, std::ostream&os=std::cout );

  //! @brief Record state and sensitivity bounds in files <a>obndsta</a> and <a>obndsa</a>, with accuracy of <a>iprec</a> digits
  void record
    ( std::ofstream&obndsta, std::ofstream&obndsa, const unsigned iprec=5 ) const;
  /** @} */

protected:
  //! @brief Function to finalize state bounding
  void _END_SA();

  //! @brief Function to initialize GSL for foward sensitivity polynomial models
  bool _INI_PM_FSA
    ( const unsigned np, const PVT*PMp, const unsigned ns );

  //! @brief Static wrapper to function to calculate the DINEQ-PMs FSA RHS values
  static int MC_GSLRHSFSAPM__
    ( double t, const double*x, double*xdot, void*user_data );

  //! @brief Static wrapper to function to calculate the DINEQ-PMs FSA RHS derivatives
  static int MC_GSLJACFSAPM__
    ( double t, const double*x, double*jac, double*xdot, void*user_data );

  //! @brief Propagate state/quadrature polynomial models forward in time through every time stages
  STATUS _bounds_FSA
    ( const unsigned ns, const double*tk, const PVT*PMp, PVT**PMxk,
      PVT*PMq, PVT*PMf, PVT**PMxpk, PVT*PMqp, PVT*PMfp, const bool store,
      std::ostream&os );

  //! @brief Private methods to block default compiler methods
  ODEBNDS_GSL(const ODEBNDS_GSL&);
  ODEBNDS_GSL& operator=(const ODEBNDS_GSL&);
};

template <typename T, typename PMT, typename PVT>
 ODEBNDS_GSL<T,PMT,PVT>* ODEBNDS_GSL<T,PMT,PVT>::_pODEBNDS = 0;

template <typename T, typename PMT, typename PVT> inline
ODEBNDS_GSL<T,PMT,PVT>::ODEBNDS_GSL
()
: ODEBND_GSL<T,PMT,PVT>(), ODEBNDS_BASE<T,PMT,PVT>(), _driver_sa(0),
  _vec_sa(0), _vec_quads(0)
{}

template <typename T, typename PMT, typename PVT> inline
ODEBNDS_GSL<T,PMT,PVT>::~ODEBNDS_GSL
()
{
  delete[] _vec_sa;
  if( _driver_sa ) gsl_odeiv2_driver_free( _driver_sa );
}

template <typename T, typename PMT, typename PVT> inline void
ODEBNDS_GSL<T,PMT,PVT>::_END_SA()
{
  // Get final CPU time
  _final_stats( stats_sa );
}

template <typename T, typename PMT, typename PVT> inline int
ODEBNDS_GSL<T,PMT,PVT>::MC_GSLRHSFSAPM__
( double t, const double* x, double* xdot, void* user_data )
{
  ODEBNDS_GSL<T,PMT,PVT> *pODEBNDS = ODEBNDS_GSL<T,PMT,PVT>::_pODEBNDS;
  bool flag = pODEBNDS->_RHS_PM_FSA( pODEBNDS->options, t, x, xdot );
  if( flag && pODEBNDS->_nq ){
    double* qdot = xdot + pODEBNDS->_offset_quads;
    flag = pODEBNDS->_RHS_PM_QUADS( pODEBNDS->options, t, x, qdot );
  }
  pODEBNDS->stats_sa.numRHS++;
  ODEBNDS_GSL<T,PMT,PVT>::_pODEBNDS = pODEBNDS;
  return( flag? GSL_SUCCESS: GSL_EBADFUNC );
}

template <typename T, typename PMT, typename PVT> inline int
ODEBNDS_GSL<T,PMT,PVT>::MC_GSLJACFSAPM__
( double t, const double* x, double* jac, double* xdot, void* user_data )
{
  ODEBNDS_GSL<T,PMT,PVT> *pODEBNDS = ODEBNDS_GSL<T,PMT,PVT>::_pODEBNDS;
  bool flag = pODEBNDS->_JAC_PM_FSA( pODEBNDS->options, t, x, jac, xdot );
  pODEBNDS->stats_sa.numJAC++;
  ODEBNDS_GSL<T,PMT,PVT>::_pODEBNDS = pODEBNDS;
  return( flag? GSL_SUCCESS: GSL_EBADFUNC );
}

template <typename T, typename PMT, typename PVT> inline bool
ODEBNDS_GSL<T,PMT,PVT>::_INI_PM_FSA
( const unsigned np, const PVT* PMp, const unsigned ns )
{

  // Initialize bound propagation
  if( !ODEBND_BASE<T,PMT,PVT>::_INI_PM_FSA( options, np, PMp, ns ) )
    return false;

  // Define ODE system in GSL format
  _sys_sa.function = MC_GSLRHSFSAPM__;
  _sys_sa.jacobian = MC_GSLJACFSAPM__;
  _sys_sa.params = 0;
  _sys_sa.dimension = _PMenv->nmon()*(_nx+_nq)+_nq*(np+1);
  switch( options.WRAPMIT){
  case Options::NONE:
    _sys_sa.dimension += _nx*(np+1);
    break;
  case Options::DINEQ:
    _sys_sa.dimension += 2*_nx*(np+1);
    break;
  case Options::ELLIPS:
  default:
    _sys_sa.dimension += _nx*(_nx+1)/2*(np+1);
    break;
  }
  ODEBND_BASE<T,PMT,PVT>::_INI_GSL( _sys_sa, _driver_sa );
  delete [] _vec_sa;
  _vec_sa = new double[ _sys_sa.dimension ];
  _offset_quadsa = _sys_sa.dimension - _PMenv->nmon()*_nq-_nq*(np+1);
  _vec_quadsa = _vec_sa + _offset_quads;

  // Initialize mesh
  if( !_mesh_fsa.set( ns, _offset_quadsa,
    options.MESHPREALLOC, PMOFFSET+(int)options.WRAPMIT ) ) return false;

  // Reset result record and statistics
  results_sa.clear();  
  _init_stats( stats_sa );

  return true;
}

template <typename T, typename PMT, typename PVT>
inline typename ODEBND_GSL<T,PMT,PVT>::STATUS
ODEBND_GSL<T,PMT,PVT>::_bounds_FSA
( const unsigned ns, const double*tk, const PVT*PMp, PVT**PMxk,
  PVT*PMq, PVT*PMf, PVT**PMxpk, PVT*PMqp, PVT*PMfp, const bool store,
  std::ostream&os )
{
  // Check arguments
  if( !tk || !PMxk || !PMp || (_nf && !PMf) ) return FATAL;

  try{
    // Initialize trajectory integration with GSL
    if( !_INI_PM_FSA( _np, PMp, ns ) ) return FATAL;

    // Bounds on initial states/quadratures and sensitivities
    _t = tk[0];
    if( !_IC_PM_FSA_STA( options, _vec_sa )
     || !_IC_PM_FSA_QUAD( options, _vec_quadsa ) )
      { _END_SA(); return FATAL; }
    if( options.DISPLAY >= 1 ){
      _print_interm( _t, _nx, _PMx, "x", os );
      _print_interm( _nq, _PMq, "q", os );
      for( unsigned ip=0; ip<_np; ip++ ){
        _print_interm( _nx, _PMxp+ip*_nx, "x_p"+ip, os );
        _print_interm( _nq, _PMqp+ip*_nq, "q_p"+ip, os );
      }
    }
    if( PMxk && !PMxk[0] ) PMxk[0] = new PVT[_nx];
    for( unsigned ix=0; PMxk[0] && ix<_nx; ix++ ) PMxk[0][ix] = _PMx[ix];

    // Record initial results
    if( options.RESRECORD )
      results_sta.push_back( Results( tk[0], _nx, PMxk[0] ) );

    // Integrate ODEs through each stage using GSL
    _h = options.H0;
    _pODEBND = this;

    for( _istg=0; _istg<ns; _istg++ ){
      // Bounds on state discontinuities (if any) at stage times
      // and integrator reinitialization (if applicable)
      _pos_ic = ( _vIC.size()>=ns? _istg:0 );
      if( _pos_ic ){
        if( !_CC_PM_STA( options, _pos_ic, _t, _vec_sta ) )
          { _END_STA(); return FATAL; }
        gsl_odeiv2_driver_reset( _driver_sta );
      }

      // update list of operations in RHS, JAC and QUAD
      _pos_rhs = ( _vRHS.size()<=1? 0: _istg );
      _pos_quad = ( _vQUAD.size()<=1? 0: _istg );
      if( (!_istg || _pos_rhs || _pos_quad)
        && !_SET_PM_STA( options, _pos_rhs, _pos_quad ) )
        { _END_STA(); return FATAL; }

      // Store mesh point at time stage
      if( store && !_mesh_sta.add( _istg, _t, _vec_sta, true ) )
        { _END_STA(); return FATAL; }

      // integrate till end of time stage
      while( _t < tk[_istg+1] ){
        if( gsl_odeiv2_evolve_apply( _driver_sta->e, _driver_sta->c, _driver_sta->s,
            &_sys_sta, &_t, tk[_istg+1], &_h, _vec_sta ) != GSL_SUCCESS
         || _h < options.HMIN
         || (options.NMAX && stats_sta.numSteps > options.NMAX)
         || _diam( _nx, _PMx ) > options.DMAX
         || _diam( _nx, _Ir  ) > options.DMAX )
          throw Exceptions( Exceptions::INTERN );
        stats_sta.numSteps++;
        if( options.HMAX > 0 && _h > options.HMAX ) _h = options.HMAX;

        // Store mesh point at intermediate time
        if( store && !_mesh_sta.add( _istg, _t, _vec_sta, false ) )
          { _END_STA(); return FATAL; }
      }

      // Bounds on intermediate states/quadratures
      switch( options.WRAPMIT){
      case Options::NONE:
      case Options::DINEQ:
        _vec2PMI( _vec_sta,  _PMenv, _nx, _PMx, false );
        break;
      case Options::ELLIPS:
      default:
        _vec2PME( _vec_sta, _PMenv, _nx, _PMx, _Q, _Er, _Ir );
        //std::cout << _Er.eigQ().first;
        break;
      }
      _vec2PMI( _vec_quad, _PMenv, _nq, _PMq, true );

      // Keep track/display/record stage results
      if( options.DISPLAY >= 1 ){
        _print_interm( _t, _nx, _PMx, "x", os );
        _print_interm( _nq, _PMq, "q", os );
      }
      if( PMxk && !PMxk[_istg+1] ) PMxk[_istg+1] = new PVT[_nx];
      for( unsigned ix=0; PMxk[_istg+1] && ix<_nx; ix++ ) PMxk[_istg+1][ix] = _PMx[ix];
      if( options.RESRECORD )
        results_sta.push_back( Results( tk[_istg+1], _nx, PMxk[_istg+1] ) );

      // Add intermediate function terms
      _pos_fct = ( _vFCT.size()>=ns? _istg:0 );
      if( (_vFCT.size()>=ns || _istg==ns-1) && !_FCT_PM_STA( _pos_fct, _t, PMf ) )
        { _END_STA(); return FATAL; }
    }

    // Bounds on final quadratures and functions
    for( unsigned iq=0; PMq && iq<_nq; iq++ ) PMq[iq] = _PMq[iq];
    if( options.DISPLAY >= 1 ) _print_interm( _nf, PMf, "f", os );
  }
  catch(...){
    _END_STA();
    if( options.DISPLAY >= 1 ) _print_stats( stats_sta, os );
    return FAILURE;
  }

  _END_STA(); 
  if( options.DISPLAY >= 1 ) _print_stats( stats_sta, os );
  return NORMAL;
}

//! @fn template <typename T, typename PMT, typename PVT> inline typename ODEBND_GSL<T,PMT,PVT>::STATUS ODEBND_GSL<T,PMT,PVT>::bounds(
//! const unsigned ns, const double*tk, const PVT*PMp, PVT**PMxk, PVT*PMq,
//! PVT*PMf, PVT**PMxpk, PVT*PMqp, PVT*PMfp, std::ostream&os=std::cout )
//!
//! This function computes an enclosure of the reachable set of the parametric ODEs
//! using propagation of polynomial models with convex remainders (intervals, ellipsoids):
//!   - <a>ns</a> [input] number of time stages
//!   - <a>tk</a> [input] stage times, including the initial time
//!   - <a>PMp</a> [input] polynomial model of parameter set
//!   - <a>PMxk</a> [output] polynomial model of state variables at stage times
//!   - <a>PMq</a> [output] polynomial model of quadrature variables
//!   - <a>PMf</a> [output] polynomial model of functionals
//!   - <a>PMxpk</a> [output] polynomial model of state sensitivities at stage times
//!   - <a>PMqp</a> [output] polynomial model of quadrature sensitivities
//!   - <a>PMfp</a> [output] polynomial model of functional derivitives
//!   - <a>os</a> [input] output stream (default: std::cout)
//! .
//! The return value is the status.
template <typename T, typename PMT, typename PVT>
inline typename ODEBND_GSL<T,PMT,PVT>::STATUS
ODEBND_GSL<T,PMT,PVT>::bounds_FSA
( const unsigned ns, const double*tk, const PVT*PMp, PVT**PMxk,
  PVT*PMq, PVT*PMf, PVT**PMxpk, PVT*PMqp, PVT*PMfp, std::ostream&os );
{
  return _bounds_FSA( ns, tk, PMp, PMxk, PMq, PMf, PMxpk, PMqp, PMfp, false, os );
}

template <typename T, typename PMT, typename PVT> inline void
ODEBNDS_GSL<T,PMT,PVT>::record
( std::ofstream&obndsta, std::ofstream&obndsa, const unsigned iprec ) const
{
  ODEBND_GSL<T,PMT,PVT>::record( obndsta, iprec );
  ODEBND_GSL<T,PMT,PVT>::_record( obndsa, results_sa, iprec ); }
}

} // end namescape mc

#endif

