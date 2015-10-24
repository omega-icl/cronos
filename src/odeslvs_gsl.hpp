// Copyright (C) 2012-2014 Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

#ifndef MC__ODESLVS_GSL_HPP
#define MC__ODESLVS_GSL_HPP

#undef  MC__ODESLVS_GSL_SAMPLE_DEBUG
#undef  MC__ODESLVS_GSL_DEBUG_INTERP
#undef  MC__ODESLVS_GSL_DEBUG_ADJOINT
#define MC__ODESLVS_GSL_CHECK

#include "odeslv_gsl.hpp"

#define MC__ODESLVS_GSL_USE_BAD

// *** TO DO
// - Implement state quadrature in adjoint sensitivity

namespace mc
{
//! @brief C++ class computing solutions of parametric ODEs and adjoint sensitivity analysis using non-validated integration.
////////////////////////////////////////////////////////////////////////
//! mc::ODESLVS_GSL is a C++ class that computes solutions of parametric
//! ordinary differential equations (ODEs) and performs adjoint
//! sensitivity analysis using GSL.
////////////////////////////////////////////////////////////////////////
template <typename T>
class ODESLVS_GSL: public ODESLV_GSL<T>
{
  using ODESLV_GSL<T>::_nx;
  using ODESLV_GSL<T>::_np;
  using ODESLV_GSL<T>::_nq;
  using ODESLV_GSL<T>::_nf;
  using ODESLV_GSL<T>::_pDAG;
  using ODESLV_GSL<T>::_pX;
  using ODESLV_GSL<T>::_pY;
  using ODESLV_GSL<T>::_pP;
  using ODESLV_GSL<T>::_pQ;
  using ODESLV_GSL<T>::_pT;
  using ODESLV_GSL<T>::_nVAR;
  using ODESLV_GSL<T>::_pVAR;
  using ODESLV_GSL<T>::_dVAR;
  using ODESLV_GSL<T>::_vRHS;
  using ODESLV_GSL<T>::_pRHS;
  using ODESLV_GSL<T>::_vQUAD;
  using ODESLV_GSL<T>::_pJAC;
  using ODESLV_GSL<T>::_dJAC;
  using ODESLV_GSL<T>::_opJAC;
  using ODESLV_GSL<T>::_vIC;
  using ODESLV_GSL<T>::_pIC;
  using ODESLV_GSL<T>::_vFCT;
  using ODESLV_GSL<T>::_mesh_sta;
  using ODESLV_GSL<T>::_vec_sta;
  using ODESLV_GSL<T>::_t;
  using ODESLV_GSL<T>::_istg;
  using ODESLV_GSL<T>::_print_interm;
  using ODESLV_GSL<T>::_init_stats;
  using ODESLV_GSL<T>::_print_stats;
  using ODESLV_GSL<T>::_END_STA;
  using ODESLV_GSL<T>::_results_sta;
  using ODESLV_GSL<T>::_states;
  using ODESLV_GSL<T>::NORMAL;
  using ODESLV_GSL<T>::FAILURE;
  using ODESLV_GSL<T>::FATAL;

private:
  //! @brief GSL data type for adjoint ODE integration
  gsl_odeiv2_system _sys_adj;

  //! @brief GSL drivers for adjoint ODE integration
  std::vector<gsl_odeiv2_driver*> _driver_adj;

  //! @brief list of operations in adjoint RHS evaluation
  std::list<const FFOp*> _opADJ;

  //! @brief list of operations in adjoint TC evaluation
  std::list<const FFOp*> _opTC;

  //! @brief const pointer to adjoint RHS function in current stage of ODE system
  const FFVar* _pADJ;

  //! @brief preallocated array for evaluation of adjoint RHS function
  double* _dADJ;

  //! @brief const pointer to adjoint TC function in current stage of ODE system
  const FFVar* _pTC;

  //! @brief array for storing TC Jacobian in current stage
  double* _dTC;

public:
  typedef BASE_GSL::STATUS STATUS;
  typedef BASE_GSL::Stats Stats;
  typedef typename ODESLV_GSL<T>::Results Results;

  /** @defgroup ODESLVS_GSL Real-valued (non-validated) integration of parametric ODEs
   *  @{
   */
  //! @brief Default class constructor
  ODESLVS_GSL()
    : ODESLV_GSL<T>(), _pADJ(0), _dADJ(0), _pTC(0), _dTC(0), _ifct(0)
    {}

  //! @brief Default destructor
  virtual ~ODESLVS_GSL()
    {    
      for( unsigned i=0; i<_driver_adj.size(); i++ )
        if( _driver_adj[i] )  gsl_odeiv2_driver_free( _driver_adj[i] );
      delete[] _dADJ;
      delete[] _dTC;
      delete[] _pADJ;
      delete[] _pTC;
    }

  //! @brief Integrator options
  struct Options: public ODESLV_GSL<T>::Options
  {
    //! @brief Constructor
    Options():
      ODESLV_GSL<T>::Options(), INTERPMETH(MESH_GSL::CSPLINE)
      {}
    //! @brief Assignment operator
    template <typename U> Options& operator=
      ( U&options ){
        ODESLV_GSL<T>::Options::operator=(options);
        INTERPMETH   = options.INTERPMETH;
        return *this;
      }
    //! @brief Numerical interpolation method
    MESH_GSL::INTERPOLATION_METHOD INTERPMETH;
  } options;

  //! @brief Statistics for adjoint integration
  Stats stats_adj;

  //! @brief Integrate trajectory of parametric ODEs
  STATUS states
    ( const unsigned ns, const double*tk, const double*p, double**xk,
      double*q, double*f, std::ostream&os=std::cout )
    { ODESLV_GSL<T>::options = options;
      return ODESLV_GSL<T>::states( ns, tk, p, xk, q, f, os ); }

  //! @brief Compute approximate interval enclosure of reachable set of parametric ODEs using parameter sampling
  STATUS bounds
    ( const unsigned ns, const double*tk, const T*Ip, T**Ixk,
      T*Iq, T*If, const unsigned nsamp, std::ostream&os=std::cout )
    { ODESLV_GSL<T>::options = options;
      return ODESLV_GSL<T>::bounds( ns, tk, Ip, Ixk, Iq, If, nsamp, os ); }

  //! @brief Integrate trajectory of parametric ODEs and adjoint ODEs
  STATUS states_ASA
    ( const unsigned ns, const double*tk, const double*p, double**xk,
      double*q, double*f, double**lk, double*df, std::ostream&os=std::cout );

  //! @brief Compute approximate interval enclosure of reachable set of parametric ODEs and adjoint ODEs using parameter sampling
  STATUS bounds_ASA
    ( const unsigned ns, const double*tk, const T*Ip, T**Ixk,
      T*Iq, T*If, T**Ilk, T*Idf, const unsigned nsamp, std::ostream&os=std::cout );

  //! @brief Record state and adjoint bounds in files <a>bndsta</a> and <a>bndadj</a>, respectively, with accuracy of <a>iprec</a> digits
  void record
    ( std::ofstream&bndsta, std::ofstream&bndadj, const unsigned iprec=5 ) const;

  //! @brief Record state bounds in file <a>bndsta</a>, with accuracy of <a>iprec</a> digits
  void record
    ( std::ofstream&bndsta, const unsigned iprec=5 ) const
    { ODESLV_GSL<T>::record( bndsta, iprec ); }

  //! @brief static pointer to class
  static ODESLVS_GSL<T> *pODESLVS_GSL;
  /** @} */

protected:
  //! @brief Vector storing interval adjoint bounds (see Options::RESRECORD)
  std::vector< Results > _results_adj;

  //! @brief current function
  unsigned _ifct;

  //! @brief vector storing stepsize during adjoint integration
  std::vector<double> _h_adj;

  //! @brief Function to initialize GSL numerical integration drivers
  void _INI_GSL
    ( gsl_odeiv2_system &sys, std::vector<gsl_odeiv2_driver*>&driver );

  //! @brief Initialize GSL for adjoint integration
  void _INI_ADJ
    ( const double* p );

  //! @brief Set adjoint RHS pointer and corresponding Jacobian
  bool _RHS_ADJ_SET
    ( const unsigned iRHS, const unsigned iQUAD, const unsigned iFCT );

  //! @brief Integrate adjoint trajectories on current stages
  STATUS _adjoints_traj
    ( const double tf, double&h );

  //! @brief Integrate adjoint trajectories on every time stages
  STATUS _adjoints
    ( const unsigned ns, const double*tk, const double*p,
      const double*const*xk, double**lk, double*df, std::ostream&os );

  //! @brief Static wrapper to function to calculate the adjoint ODEs RHS values
  static int MC_GSLADJRHS__
    ( double t, const double* l, double* ldot, void* user_data );

  //! @brief Function to calculate the adjoint ODEs RHS values
  int _RHS_ADJ
    ( double t, const double* l, double* ldot, void* user_data );

  //! @brief Static wrapper to function to calculate the adjoint ODEs RHS derivatives
  static int MC_GSLADJJAC__
    ( double t, const double* l, double* jac, double* xdot, void* user_data );

  //! @brief Function to calculate the adjoint ODEs RHS derivatives
  int _JAC_ADJ
    ( double t, const double* l, double* jac, double* xdot, void* user_data );

  //! @brief Recursive function computing bounds on solutions of IVP in ODEs using sampling
  STATUS _states_ASA
    ( const unsigned ns, const double*tk, const T*Ip, T**Ixk, T*Iq,
      T*If, T**Ilk, T*Idf, const unsigned nsamp, unsigned int* vsamp,
      const unsigned ip, double*p, double**xk, double*q, double*f,
      double**lk, double*df, std::ostream&os );

private:
  //! @brief Private methods to block default compiler methods
  ODESLVS_GSL(const ODESLVS_GSL&);
  ODESLVS_GSL& operator=(const ODESLVS_GSL&);
};

template <typename T>
ODESLVS_GSL<T>* ODESLVS_GSL<T>::pODESLVS_GSL = 0;

template <typename T> inline void
ODESLVS_GSL<T>::_INI_ADJ
( const double*p )
{
  // Define adjoint ODE system in GSL format
  _sys_adj.function = MC_GSLADJRHS__;
  _sys_adj.jacobian = MC_GSLADJJAC__;
  _sys_adj.dimension = _nx+_np;
  _sys_adj.params = 0;
  
  // Set GSL drivers for adjoint ODE integration
  _INI_GSL( _sys_adj, _driver_adj );

  // Set intermediate storage
  delete [] _vec_sta; _vec_sta = new double[ _sys_adj.dimension ];
  delete [] _dTC; _dTC = new double[ _sys_adj.dimension*_nf ];

  // Size and set DAG evaluation arrays
  if( _nVAR != 2*_nx+_np+1 ){
    _nVAR = 2*_nx+_np+1;
    delete[] _pVAR; _pVAR = new FFVar[_nVAR];
    delete[] _dVAR; _dVAR = new double[_nVAR];
  }

  BASE_DE::set_adjoint();
  for( unsigned ix=0; ix<_nx; ix++ ) _pVAR[ix] = _pY[ix];
  for( unsigned ix=0; ix<_nx; ix++ ) _pVAR[_nx+ix] = _pX[ix];
  for( unsigned ip=0; ip<_np; ip++ ) _pVAR[2*_nx+ip] = _pP[ip];
  _pVAR[2*_nx+_np] = (_pT? *_pT: 0. );

  for( unsigned ip=0; ip<_np; ip++ ) _dVAR[2*_nx+ip] = p[ip];

  // Initialize statistics
  _init_stats( stats_adj );

  return;
}

template <typename T> inline void
ODESLVS_GSL<T>::_INI_GSL
( gsl_odeiv2_system &sys, std::vector<gsl_odeiv2_driver*>&driver )
{
  // Reset GSL numerical integration drivers
  for( unsigned i=0; i<driver.size(); i++ )
    gsl_odeiv2_driver_free( driver[i] );
  driver.clear();

  for( unsigned i=0; i<_nf; i++ ){
    switch( options.INTMETH ){
    case Options::RK8PD:
      driver.push_back( gsl_odeiv2_driver_alloc_y_new( &sys,
        gsl_odeiv2_step_rk8pd, options.H0, options.ATOL, options.RTOL ) );
      break;
    case Options::MSADAMS:
      driver.push_back( gsl_odeiv2_driver_alloc_y_new( &sys,
        gsl_odeiv2_step_msadams, options.H0, options.ATOL, options.RTOL ) );
      break;
    case Options::MSBDF:
      driver.push_back( gsl_odeiv2_driver_alloc_y_new( &sys,
        gsl_odeiv2_step_msbdf, options.H0, options.ATOL, options.RTOL ) );
      break;
    case Options::RKF45: default:
      driver.push_back( gsl_odeiv2_driver_alloc_y_new( &sys,
        gsl_odeiv2_step_rkf45, options.H0, options.ATOL, options.RTOL ) );
      break;
    }
    gsl_odeiv2_driver_set_hmin( driver[i], options.HMIN );  
    gsl_odeiv2_driver_set_nmax( driver[i], options.NMAX );  
  }

  return;
}

template <typename T> inline int
ODESLVS_GSL<T>::_RHS_ADJ
( double t, const double* y, double* ydot, void* user_data )
{
  stats_adj.numRHS++; // increment RHS counter
  if( !_pADJ ) return GSL_EBADFUNC; // **error** RHS not defined
  if( !_mesh_sta.eval( _istg, -t, _dVAR+_nx ) ) return GSL_EBADFUNC; // set interpolated state values
  for( unsigned i=0; i<_nx; i++ ) _dVAR[i] = y[i]; // set current adjoint values
  _dVAR[2*_nx+_np] = -t; // set current time
  _pDAG->eval( _opADJ, _dADJ, _nx+_np, _pADJ, ydot, _nVAR, _pVAR, _dVAR );

  return GSL_SUCCESS;  
}

template <typename T> inline int
ODESLVS_GSL<T>::MC_GSLADJRHS__
( double t, const double* y, double* ydot, void* user_data )
{
  ODESLVS_GSL<T> *pODESLVS_GSL = ODESLVS_GSL<T>::pODESLVS_GSL;
  int flag = pODESLVS_GSL->_RHS_ADJ( t, y, ydot, user_data );
  ODESLVS_GSL<T>::pODESLVS_GSL = pODESLVS_GSL;
  return flag;
}

template <typename T> inline int
ODESLVS_GSL<T>::_JAC_ADJ
( double t, const double* y, double* jac, double* ydot, void* user_data )
{
  stats_adj.numJAC++;
  if( !_pADJ || !_pJAC ) return GSL_EBADFUNC; // **error** RHS or JAC not defined
  if( !_mesh_sta.eval( _istg, -t, _dVAR+_nx ) ) return GSL_EBADFUNC; // set interpolated state values
  for( unsigned i=0; i<_nx; i++ ) _dVAR[i] = y[i]; // set current adjoint/quadrature values
  _dVAR[2*_nx+_np] = -t; // set current time
  _pDAG->eval( _opADJ, _dADJ, _nx+_np, _pADJ, ydot, _nVAR, _pVAR, _dVAR );
  _pDAG->eval( _opJAC, _dJAC, (_nx+_np)*(_nx+_np), _pJAC, jac, _nVAR, _pVAR, _dVAR );

  return GSL_SUCCESS;
}

template <typename T> inline int
ODESLVS_GSL<T>::MC_GSLADJJAC__
( double t, const double* y, double* jac, double* ydot, void* user_data )
{
  ODESLVS_GSL<T> *pODESLVS_GSL = ODESLVS_GSL<T>::pODESLVS_GSL;
  int flag = pODESLVS_GSL->_JAC_ADJ( t, y, jac, ydot, user_data );
  ODESLVS_GSL<T>::pODESLVS_GSL = pODESLVS_GSL;
  return flag;
}

template <typename T> inline bool
ODESLVS_GSL<T>::_RHS_ADJ_SET
( const unsigned iRHS, const unsigned iQUAD, const unsigned iFCT )
{
  if( _vRHS.size() <= iRHS ) return false; 
  if( _nq && _vQUAD.size() <= iQUAD ) return false;

  FFVar pHAM( 0. );
  delete[] _pRHS; _pRHS = new FFVar[_nx+_nq];
  for( unsigned i=0; i<_nx; i++ ) _pRHS[i] = _vRHS.at( iRHS )[i];
  for( unsigned ix=0; ix<_nx; ix++ ) pHAM += _pVAR[ix] * _pRHS[ix];
  for( unsigned i=0; i<_nq; i++ ) _pRHS[_nx+i] = _vQUAD.at( iQUAD )[i];
  const FFVar* pFCT = _vFCT.at(iFCT)+_ifct;
#ifndef MC__ODESLVS_GSL_USE_BAD
  delete[] _pTC; _pTC = _nq? _pDAG->FAD( 1, pFCT, _nq, _pQ ): 0;
#else
  delete[] _pTC; _pTC = _nq? _pDAG->BAD( 1, pFCT, _nq, _pQ ): 0;
#endif
  for( unsigned iq=0; iq<_nq; iq++ ){
    if( !_pTC[iq].cst() ) return false; // quadrature appears nonlinearly in function
    pHAM += _pRHS[_nx+iq] * _pTC[iq];
  }

#ifndef MC__ODESLVS_GSL_USE_BAD
  delete[] _pADJ; _pADJ = _pDAG->FAD( 1, &pHAM, _nx+_np, _pVAR+_nx );
#else
  delete[] _pADJ; _pADJ = _pDAG->BAD( 1, &pHAM, _nx+_np, _pVAR+_nx );
#endif
  _opADJ = _pDAG->subgraph( _nx+_np, _pADJ );
  delete[] _dADJ; _dADJ = new double[ _opADJ.size() ];
  
  switch( options.INTMETH ){
  case Options::MSADAMS:
  case Options::MSBDF:
    delete[] _pJAC; _pJAC  = _pDAG->FAD( _nx+_np, _pADJ, _nx+_np, _pVAR+_nx );
    _opJAC = _pDAG->subgraph( (_nx+_np)*(_nx+_np), _pJAC );
    delete[] _dJAC; _dJAC  = new double[ _opJAC.size() ];
    break;
  case Options::RK8PD:
  case Options::RKF45:
  default:
    delete[] _pJAC; _pJAC = 0;
    _opJAC.clear();
    delete[] _dJAC; _dJAC = 0;
    break;
  }

  return true;
}

template <typename T> inline typename ODESLVS_GSL<T>::STATUS
ODESLVS_GSL<T>::_adjoints_traj
( const double tf, double&h )
{
  // integrate till end of time stage
#ifdef MC__ODESLVS_GSL_DEBUG_ADJOINT
    std::cout << "adjoint #" << _ifct << "  " << "stage #" << _istg << "  times:\n";
#endif
  while( _t < tf ){
    if( gsl_odeiv2_evolve_apply( _driver_adj[_ifct]->e, _driver_adj[_ifct]->c,
        _driver_adj[_ifct]->s, &_sys_adj, &_t, tf, &h, _vec_sta ) != GSL_SUCCESS
     || h < options.HMIN
     || (options.NMAX && stats_adj.numSteps > options.NMAX) ) return FAILURE;
    stats_adj.numSteps++;
    if( options.HMAX > 0 && h > options.HMAX ) h = options.HMAX;
#ifdef MC__ODESLVS_GSL_DEBUG_ADJOINT
    std::cout << -_t << std::endl;
#endif
  }
  return NORMAL;
}

template <typename T> inline typename ODESLVS_GSL<T>::STATUS
ODESLVS_GSL<T>::_adjoints
( const unsigned ns, const double*tk, const double*p, const double*const*xk,
  double**lk, double*q, std::ostream&os )
{
  // Initialize trajectory integration with GSL
  _INI_ADJ( p );

  // Adjoint terminal conditions
  if( _vFCT.size() != 1 && _vFCT.size() < ns )
    { _END_STA( stats_adj ); return FATAL; }
  unsigned iFCT = ( _vFCT.size()>1? ns-1:0 );
  const FFVar* pFCT = _vFCT.at(iFCT);
#ifndef MC__ODESLVS_GSL_USE_BAD
  delete[] _pTC; _pTC = _pDAG->FAD( _nf, pFCT, _nx+_np, _pVAR+_nx );
#else
  delete[] _pTC; _pTC = _pDAG->BAD( _nf, pFCT, _nx+_np, _pVAR+_nx );
#endif
  _opTC = _pDAG->subgraph( (_nx+_np)*_nf, _pTC );

  for( unsigned i=0; i<_nx; i++ ) _dVAR[_nx+i] = xk[ns][i]; // set current state values
  _dVAR[2*_nx+_np] = tk[ns]; // set current time
  _pDAG->eval( _opTC, (_nx+_np)*_nf, _pTC, _dTC, _nx+_np+1, _pVAR+_nx, _dVAR+_nx );

  for( unsigned k=0, ik=0; k<_nf; k++ ){
    for( unsigned i=0; i<_nx; i++, ik++ ) lk[ns][k*_nx+i] = _dTC[ik];
    for( unsigned i=0; i<_np; i++, ik++ ) q[k*_np+i] = _dTC[ik];
  }
  if( options.DISPLAY >= 1 ){
    _print_interm( tk[ns], _nf*_nx, lk[ns], "l", os );
    //_print_interm( tk[ns], _nf*_np, q, "q", os );
  }

  // Integrate adjoint ODEs backward through each stage using GSL
  _h_adj.assign( _nf, options.H0 );
  pODESLVS_GSL = this;

  for( _istg=ns; _istg>0; _istg-- ){

    // Adjoint discontinuity at stage times
    iFCT = ( _vFCT.size()>1? _istg-1:0 );
    if( _istg<ns ){

      if( _vFCT.size() > 1 ){
        pFCT = _vFCT.at(iFCT);
#ifndef MC__ODESLVS_GSL_USE_BAD
        delete[] _pTC; _pTC = _pDAG->FAD( _nf, pFCT, _nx+_np, _pVAR+_nx );
#else
        delete[] _pTC; _pTC = _pDAG->BAD( _nf, pFCT, _nx+_np, _pVAR+_nx );
#endif
        _opTC = _pDAG->subgraph( (_nx+_np)*_nf, _pTC );

        _dVAR[2*_nx+_np] = tk[_istg]; // set current time
        for( unsigned i=0; i<_nx; i++ ) _dVAR[_nx+i] = xk[_istg][i]; // set current state values
        _pDAG->eval( _opTC, (_nx+_np)*_nf, _pTC, _dTC, _nx+_np+1, _pVAR+_nx, _dVAR+_nx );
        // Reset adjoint ODE solver
        for( _ifct=0; _ifct<_nf; _ifct++ ) gsl_odeiv2_driver_reset( _driver_adj[_ifct] );
      }
      else
        for( unsigned i=0; i<(_nx+_np)*_nf; i++ ) _dTC[i] = 0.;

      // State discontinuity contribution
      if( _vIC.size() > 1 ){
        _pIC = _vIC.at(_istg);
        FFVar pHAM( 0. );
        for( unsigned ix=0; ix<_nx; ix++ ) pHAM += _pVAR[ix] * _pIC[ix];
#ifndef MC__ODESLVS_GSL_USE_BAD
        delete[] _pTC; _pTC = _pDAG->FAD( 1, &pHAM, _nx+_np, _pVAR+_nx );
#else
        delete[] _pTC; _pTC = _pDAG->BAD( 1, &pHAM, _nx+_np, _pVAR+_nx );
#endif
        _opTC = _pDAG->subgraph( _nx+_np, _pTC );

        _dVAR[2*_nx+_np] = tk[_istg]; // set current time
        for( _ifct=0; _ifct<_nf; _ifct++ ){
          for( unsigned i=0; i<_nx; i++ ){
            _dVAR[i] = lk[_istg][_ifct*_nx+i] + _dTC[_ifct*(_nx+_np)+i]; // l(ti+)+dphi/dx(ti)
            _dTC[_ifct*(_nx+_np)+i] = 0.;
            _dVAR[_nx+i] = xk[_istg][i]; // set current state values
          }
          for( unsigned i=0; i<_np; i++ )
            _dTC[_ifct*(_nx+_np)+_nx+i] += q[_ifct*_np+i];
          _pDAG->eval( _opTC, _nx+_np, _pTC, _dTC+_ifct*(_nx+_np), _nVAR, _pVAR, _dVAR, true ); // add result to _dTC
          // Reset adjoint ODE solver
          gsl_odeiv2_driver_reset( _driver_adj[_ifct] );
        }
      }
      else
        for( _ifct=0; _ifct<_nf; _ifct++ ){
          for( unsigned i=0; i<_nx; i++ ) _dTC[_ifct*(_nx+_np)+i] += lk[_istg][_ifct*_nx+i];
          for( unsigned i=0; i<_np; i++ ) _dTC[_ifct*(_nx+_np)+_nx+i] += q[_ifct*_np+i];
        }
    }

    // Interpolate state mesh
    if( !_mesh_sta.interp( _istg, options.INTERPMETH ) )
      { _END_STA( stats_adj ); return FAILURE; }

    // Integrate backward along time stage for each function
    for( _ifct=0; _ifct<_nf; _ifct++ ){
      // Update list of operations in adjoint RHS and JAC
      const unsigned iRHS  = ( _vRHS.size()<2?  0: _istg-1 );
      const unsigned iQUAD = ( _vQUAD.size()<2? 0: _istg-1 );
      //if( (_istg==ns || iRHS) && !_RHS_ADJ_SET( iRHS ) )
      if( !_RHS_ADJ_SET( iRHS, iQUAD, iFCT ) )
        { _END_STA( stats_adj ); return FATAL; }

      _t = -tk[_istg]; // need to reinitialize stepsize too???
      for( unsigned i=0; i<_nx+_np; i++ ) _vec_sta[i] = _dTC[_ifct*(_nx+_np)+i];
      if( _adjoints_traj( -tk[_istg-1], _h_adj[_ifct] ) == FAILURE )
        { _END_STA( stats_adj ); return FAILURE; }
      for( unsigned i=0; i<_nx; i++ ) lk[_istg-1][_ifct*_nx+i] = _vec_sta[i];
      for( unsigned i=0; i<_np; i++ ) q[_ifct*_np+i] = _vec_sta[_nx+i];
    }

    // Initial state contribution
    if( _istg==1 ){
      _pIC = _vIC.at(0);
      FFVar pHAM( 0. );
      for( unsigned ix=0; ix<_nx; ix++ ) pHAM += _pVAR[ix] * _pIC[ix];
#ifndef MC__ODESLVS_GSL_USE_BAD
      delete[] _pTC; _pTC = _pDAG->FAD( 1, &pHAM, _np, _pVAR+2*_nx );
#else
      delete[] _pTC; _pTC = _pDAG->BAD( 1, &pHAM, _np, _pVAR+2*_nx );
#endif
      _opTC = _pDAG->subgraph( _np, _pTC );

      _dVAR[2*_nx+_np] = tk[0]; // set initial time
      for( _ifct=0; _ifct<_nf; _ifct++ ){
        for( unsigned i=0; i<_nx; i++ ) _dVAR[i] = lk[0][_ifct*_nx+i]; // set initial adjoint
        _pDAG->eval( _opTC, _np, _pTC, q+_ifct*_np, _nVAR, _pVAR, _dVAR, true ); // add result to q
      }
    }

    if( options.DISPLAY >= 1 )
      _print_interm( tk[_istg-1], _nf*_nx, lk[_istg-1], "l", os );
  }

  if( options.DISPLAY >= 1 )
    _print_interm( _nf*_np, q, "df", os );

  _END_STA( stats_adj );
  if( options.DISPLAY >= 1 ) _print_stats( stats_adj, os );
  return NORMAL;
}

template <typename T> inline typename ODESLVS_GSL<T>::STATUS
ODESLVS_GSL<T>::_states_ASA
( const unsigned ns, const double*tk, const T*Ip, T**Ixk, T*Iq,
  T*If, T**Ilk, T*Idf, const unsigned nsamp, unsigned int* vsamp,
  const unsigned ip, double*p, double**xk, double*q, double*f,
  double**lk, double*df, std::ostream&os )
{
  STATUS flag = NORMAL;

  // Update bounds for all sampling points
  for( unsigned isamp=0; isamp<nsamp; isamp++ ){
    vsamp[ip] = isamp;

    // Continue recursive call
    if( ip+1 < _np ){
      flag = _states_ASA( ns, tk, Ip, Ixk, Iq, If, Ilk, Idf, nsamp, vsamp,
                          ip+1, p, xk, q, f, lk, df, os );
      if( flag != NORMAL ) return flag;
      continue;
    }

    // Update bounds for current point
#ifdef MC__ODESLVS_GSL_SAMPLE_DEBUG
    std::cout << "Sample: ";
#endif
    for( unsigned ip=0; ip<_np; ip++ ){
      p[ip] = Op<T>::l( Ip[ip] ) + vsamp[ip]/(nsamp-1.) * Op<T>::diam( Ip[ip] );
#ifdef MC__ODESLVS_GSL_SAMPLE_DEBUG
      std::cout << p[ip] << "  ";
#endif
    }
#ifdef MC__ODESLVS_GSL_SAMPLE_DEBUG
    std::cout << std::endl;
#endif
    flag = states_ASA( ns, tk, p, xk, q, f, lk, df, os );
    if( flag != NORMAL ) return flag;
    for( unsigned is=0; is<=ns; is++ ){
      for( unsigned ix=0; ix<_nx; ix++ )
        Ixk[is][ix] = Op<T>::hull( xk[is][ix], Ixk[is][ix] );
      for( unsigned ix=0; ix<_nf*_nx; ix++ )
        Ilk[is][ix] = Op<T>::hull( lk[is][ix], Ilk[is][ix] );
    }
    for( unsigned iq=0; Iq && iq<_nq; iq++ )
      Iq[iq] = Op<T>::hull( q[iq], Iq[iq] );
    for( unsigned ifn=0; If && ifn<_nf; ifn++ ){
      If[ifn] = Op<T>::hull( f[ifn], If[ifn] );
      for( unsigned ip=0; ip<_np; ip++ )
        Idf[ifn*_np+ip] = Op<T>::hull( df[ifn*_np+ip], Idf[ifn*_np+ip] );
    }
  }

  return flag;
}

//! @fn template <typename T> inline typename ODESLVS_GSL<T>::STATUS ODESLVS_GSL<T>::states_ASA(
//! const unsigned ns, const double*tk, const double*p, double**xk,
//! double*q, double*f, double**lk, double*df, std::ostream&os )
//!
//! This function computes the solution of the parametric ODEs defined in IVP:
//!   - <a>ns</a> [input] number of time stages
//!   - <a>tk</a> [input] stage times, including the initial time
//!   - <a>p</a> [input] parameter values
//!   - <a>xk</a> [output] state values at stage times
//!   - <a>q</a>  [output] quadrature values at final time (only if q != 0)
//!   - <a>f</a>  [output] function values (only if f != 0)
//!   - <a>lk</a> [output] adjoint values at stage times
//!   - <a>f</a>  [output] function derivatives (stored row-wise)
//!   - <a>os</a> [input] output stream
//!   .
//! The return value is the status.
template <typename T> inline typename ODESLVS_GSL<T>::STATUS
ODESLVS_GSL<T>::states_ASA
( const unsigned ns, const double*tk, const double*p, double**xk,
  double*q, double*f, double**lk, double*df, std::ostream&os )
{
  ODESLV_GSL<T>::options = options;
  STATUS flag = NORMAL;

  flag = _states( ns, tk, p, xk, q, f, true, os );
  if( flag != NORMAL ) return flag;

  flag = _adjoints( ns, tk, p, xk, lk, df, os );
  return flag;
}

//! @fn template <typename T> inline typename ODESLVS_GSL<T>::STATUS ODESLVS_GSL<T>::bounds_ASA(
//! const unsigned ns, const double*tk, const T*Ip, T**Ixk, T*Iq, T*If,
//! T**Ilk, T*df, const unsigned nsamp, std::ostream&os )
//!
//! This function computes an approximate interval enclosure of the
//! reachable set of the parametric ODEs defined in IVP using equally
//! spaced samples:
//!   - <a>ns</a> [input] number of time stages
//!   - <a>tk</a> [input] stage times, including the initial time
//!   - <a>Ip</a> [input] interval parameter set
//!   - <a>Ixk</a> [output] approximate interval state enclosures at stage times
//!   - <a>If</a>  [output] approximate quadrature enclosure (only if Iq != 0)
//!   - <a>If</a>  [output] approximate function enclosure (only if If != 0)
//!   - <a>Ilk</a> [output] approximate interval adjoint enclosures at stage times
//!   - <a>Idf</a>  [output] approximate function derivative enclosures (stored row-wise)
//!   - <a>nsamp</a> [input] number of samples for each parameter
//!   - <a>os</a> [input] output stream
//! .
//! The return value is the status.
template <typename T> inline typename ODESLVS_GSL<T>::STATUS
ODESLVS_GSL<T>::bounds_ASA
( const unsigned ns, const double*tk, const T*Ip, T**Ixk,
  T*Iq, T*If, T**Ilk, T*Idf, const unsigned nsamp, std::ostream&os )
{
  int DISPLAY_SAVE = options.DISPLAY;
  options.DISPLAY = 0;
  ODESLV_GSL<T>::options = options;
  STATUS flag = NORMAL;
  
  // Initialization of sampled bounds at parameter lower bound
  double *p = new double[_np];
  for( unsigned ip=0; ip<_np; ip++ )
    p[ip] = Op<T>::l(Ip[ip]);
  double **xk = new double*[ns+1];
  double **lk = new double*[ns+1];
  for( unsigned is=0; is<=ns; is++ ){
    xk[is] = new double[_nx];
    lk[is] = new double[_nf*_nx];
  }
  double *q = Iq&&_nq? new double[_nq]: 0;
  double *f = If&&_nf? new double[_nf]: 0;
  double *df = new double[_nf*_nx];
  flag = states_ASA( ns, tk, p, xk, q, f, lk, df, os );
  if( flag != NORMAL || nsamp <= 1 ){
    delete[] p;
    for( unsigned is=0; is<=ns; is++ ) delete[] xk[is]; delete[] xk;
    for( unsigned is=0; is<=ns; is++ ) delete[] lk[is]; delete[] lk;
    delete[] f; delete[] df;
    return flag;
  }   
  for( unsigned is=0; is<=ns; is++ ){
    for( unsigned ix=0; ix<_nx; ix++ )
      Ixk[is][ix] = xk[is][ix];
    for( unsigned ix=0; ix<_nf*_nx; ix++ )
      Ilk[is][ix] = lk[is][ix];
  }
  for( unsigned iq=0; Iq && iq<_nq; iq++ )
    Iq[iq] = q[iq];
  for( unsigned ifn=0; If && ifn<_nf; ifn++ ){
    If[ifn] = f[ifn];
    for( unsigned ip=0; ip<_np; ip++ )
      Idf[ifn*_np+ip] = df[ifn*_np+ip];
  }
  
  // Start sampling process
  unsigned int* vsamp = new unsigned int[_np];
  flag = _states_ASA( ns, tk, Ip, Ixk, Iq, If, Ilk, Idf, nsamp, vsamp,
                      0, p, xk, q, f, lk, df, os );

  // Display results
  options.DISPLAY = DISPLAY_SAVE;
  if( options.DISPLAY >= 1 ){
    for( unsigned is=0; is<=ns; is++ )
      _print_interm( tk[is], _nx, Ixk[is], "x", os );
    if( Iq ) _print_interm( _nq, Iq, "q", os );
    if( If ) _print_interm( _nf, If, "f", os );
    for( unsigned is=0; is<=ns; is++ )
      _print_interm( tk[ns-is], _nx*_nf, Ilk[ns-is], "l", os );
    _print_interm( _nf*_np, Idf, "df", os );
  }

  // Record intermediate results
  _results_sta.clear();
  _results_adj.clear();
  if( options.RESRECORD )
    for( unsigned is=0; is<=ns; is++ ){
      _results_sta.push_back( Results( tk[is], _nx, Ixk[is] ) );
      _results_adj.push_back( Results( tk[is], _nf*_nx, Ilk[is] ) );
    }
  // Clean-up
  delete[] p; delete[] q; delete[] f; delete[] df;
  for( unsigned is=0; is<=ns; is++ ) delete[] xk[is]; delete[] xk;
  for( unsigned is=0; is<=ns; is++ ) delete[] lk[is]; delete[] lk;
  delete[] vsamp;
  
  return flag;
}

template <typename T> inline void
ODESLVS_GSL<T>::record
( std::ofstream&bndsta, std::ofstream&bndadj, const unsigned iprec ) const
{
  ODESLV_GSL<T>::record( bndsta, iprec );
  if( !bndadj ) return;

  // Specify format
  bndadj << std::right << std::scientific << std::setprecision(iprec);

  // Record computed adjoint interval bounds at stage times
  typename std::vector< Results >::const_iterator it = _results_adj.begin();
  for( ; it != _results_adj.end(); ++it ){
    bndadj << std::setw(iprec+9) << (*it).t;
    for( unsigned ix=0; ix<(*it).nx; ix++ )
      bndadj << std::setw(iprec+9) << mc::Op<T>::l( (*it).X[ix] )
             << std::setw(iprec+9) << mc::Op<T>::u( (*it).X[ix] );
    bndadj << std::endl;
  }
}

} // end namescape mc

#endif

