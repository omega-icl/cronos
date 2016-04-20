// Copyright (C) 2012-2016 Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

#ifndef MC__BASE_DE_HPP
#define MC__BASE_DE_HPP

#include "ffunc.hpp"

namespace mc
{
//! @brief C++ base class for defining of parametric IVP-DAEs
////////////////////////////////////////////////////////////////////////
//! mc::BASE_DE is a C++ base class for defining the variables,
//! parameters and functions participating in parametric differential-
//! algebraic equations (DAEs)
////////////////////////////////////////////////////////////////////////
class BASE_DE
{
protected:
  //! @brief pointer to DAG of IVP-DAE system
  FFGraph* _pDAG;

  //! @brief pointer to time/independent variable
  const FFVar* _pT;

  //! @brief number of differential states
  unsigned _nd;

  //! @brief number of algebraic states
  unsigned _na;

  //! @brief pointer to state time derivatives
  FFVar* _pDX;

  //! @brief number of states
  unsigned _nx;

  //! @brief pointer to state variables - first _nd are differential states, next _na algebraic states
  const FFVar* _pX;

  //! @brief number of initial states
  unsigned _nx0;

  //! @brief number of sensitivity/adjoint variables
  unsigned _ny;

  //! @brief pointer to sensitivity/adjoint variables
  FFVar* _pY;

  //! @brief number of parameters
  unsigned _np;

  //! @brief pointer to parameters
  const FFVar* _pP;

  //! @brief number of quadratures
  unsigned _nq;

  //! @brief pointer to quadrature variables
  const FFVar* _pQ;

  //! @brief number of sensitivity/adjoint quadratures
  unsigned _nyq;

  //! @brief pointer to sensitivity/adjoint quadratures
  FFVar* _pYQ;

  //! @brief vector of const pointers to initial value in each stage of IVP-DAE system
  std::vector<const FFVar*> _vIC;

  //! @brief vector of const pointers to RHS of differential equations in each stage of IVP-DAE system
  std::vector<const FFVar*> _vRHS;

  //! @brief vector of const pointers to algebraic equations in each stage of IVP-DAE system
  std::vector<const FFVar*> _vAE;

  //! @brief vector of const pointers to integrand of quadratures in each stage of IVP-DAE system
  std::vector<const FFVar*> _vQUAD;

  //! @brief number of invariants
  unsigned _ni;

  //! @brief vector of const pointers to invariant equations in each stage of IVP-DAE system
  std::vector<const FFVar*> _vINV;

  //! @brief number of state functionals
  unsigned _nf;

  //! @brief vector of const pointers to functionals (assumes additive contributions) in each stage of IVP-DAE system
  std::vector<const FFVar*> _vFCT;

  //! @brief current time
  double _t;

  //! @brief current stage
  unsigned _istg;

public:
  /** @ingroup ODESLV_GSL
   *  @ingroup ODEBND_GSL
   *  @ingroup ODEBND_VAL
   *  @{
   */

  //! @brief Class constructor
  BASE_DE()
    : _pDAG(0), _pT(0), _nd(0), _na(0), _pDX(0), _nx(0), _pX(0),
      _nx0(0), _ny(0), _pY(0), _np(0), _pP(0), _nq(0), _pQ(0),
      _nyq(0), _pYQ(0), _ni(0), _nf(0)
    {}

  //! @brief Class destructor
  virtual ~BASE_DE()
    { delete[] _pY; delete[] _pYQ; delete[] _pDX; }

  //! @brief Integrator status
  enum STATUS{
     NORMAL=0,	//!< Normal execution
     FAILURE,	//!< Integration breakdown (bounds explosion)
     FATAL	//!< Interruption due to errors in third-party libraries
  };

  //! @brief Get pointer to DAG
  FFGraph* dag() const
    { return _pDAG; }

  //! @brief Set pointer to DAG
  void set_dag
    ( FFGraph*pDAG )
    { _pDAG = pDAG; }

  //! @brief Get size of state
  unsigned nx() const
    { return _nx; }

  //! @brief Get size of initial condition
  unsigned nx0() const
    { return _nx0; }

  //! @brief Get size of differential state
  unsigned nd() const
    { return _nd; }

  //! @brief Get size of algebraic state
  unsigned na() const
    { return _na; }

  //! @brief Get pointer to state
  const FFVar* state() const
    { return _pX; }

  //! @brief Get pointer to state
  const FFVar* state_deriv() const
    { return _pDX; }

  //! @brief Set pointer to state
  void set_state
    ( const unsigned nx, const FFVar*pX, const FFVar*pDX=0 )
    { if( nx != _nx){ delete[] _pDX; _pDX = (nx?new FFVar[nx]:0); }
      _nx = nx; _pX = pX;
      for( unsigned ix=0; ix<_nx; ix++ ){
        if( pDX ) _pDX[ix] = pDX[ix];
        else if( _pDX[ix].dag() != _pDAG ) _pDX[ix].set( _pDAG );
      } }

  //! @brief Return size of parameter
  unsigned np() const
    { return _np; }

  //! @brief Get pointer to parameter
  const FFVar* parameter() const
    { return _pP; }

  //! @brief Set pointer to parameter
  void set_parameter
    ( const unsigned np, const FFVar*pP )
    { _np = np; _pP = pP; }

  //! @brief Get pointer to time
  const FFVar* time() const
    { return _pT; }

  //! @brief Set pointer to time
  void set_time
    ( const FFVar*pT )
    { _pT = pT; }

  //! @brief Define RHS of differential equations in single-stage IVP-DAE
  void set_differential
    ( const unsigned nd, const FFVar*const RHS )
    { _nd = nd; _vRHS.clear(); _vRHS.push_back( RHS ); }

  //! @brief Define RHS of differential equations in multi-stage IVP-DAE with <a>ns</a> stages
  void set_differential
    ( const unsigned ns, const unsigned nd, const FFVar*const RHS )
    { _nd = nd; _vRHS.clear(); for( unsigned i=0; i<ns; i++ ) _vRHS.push_back( RHS+i*_nd ); }

  //! @brief Define algebraic equations in single-stage IVP-DAE
  void set_algebraic
    ( const unsigned na, const FFVar*const AE )
    { _na = na; _vAE.clear(); _vAE.push_back( AE ); }

  //! @brief Define algebraic equations in multi-stage IVP-DAE with <a>ns</a> stages
  void set_algebraic
    ( const unsigned ns, const unsigned na, const FFVar*const AE )
    { _na = na; _vAE.clear(); for( unsigned i=0; i<ns; i++ ) _vAE.push_back( AE+i*_na ); }

  //! @brief Define quadrature equations in single-stage IVP-DAE
  void set_quadrature
    ( const unsigned nq, const FFVar*const QUAD, const FFVar*pQ )
    { _nq = nq; _pQ = pQ; _vQUAD.clear(); _vQUAD.push_back( QUAD ); }

  //! @brief Define quadrature equations in multi-stage IVP-DAE with <a>ns</a> stages
  void set_quadrature
    ( const unsigned ns, const unsigned nq, const FFVar*const QUAD, const FFVar*pQ )
    { _nq = nq; _pQ = pQ; _vQUAD.clear(); for( unsigned i=0; i<ns; i++ ) _vQUAD.push_back( QUAD+i*_nq ); }

  //! @brief Define invariant equations in single-stage IVP-DAE
  void set_invariant
    ( const unsigned ni, const FFVar*const INV )
    { _ni = ni; _vINV.clear(); _vINV.push_back( INV ); }

  //! @brief Define invariant equations in multi-stage IVP-DAE with <a>ns</a> stages
  void set_invariant
    ( const unsigned ns, const unsigned ni, const FFVar*const INV )
    { _ni = ni; _vINV.clear(); for( unsigned i=0; i<ns; i++ ) _vINV.push_back( INV+i*_ni ); }

  //! @brief Define IC of differential variables in single-stage IVP-DAE
  void set_initial
    ( const unsigned nx0, const FFVar*const IC )
    { _nx0 = nx0; _vIC.clear(); _vIC.push_back( IC ); }

  //! @brief Define IC of differential variables in multi-stage IVP-DAE with <a>ns</a> stages
  void set_initial
    ( const unsigned ns, const unsigned nx0, const FFVar*const IC )
    { _nx0 = nx0; _vIC.clear(); for( unsigned i=0; i<ns; i++ ) _vIC.push_back( IC+i*_nx0 ); }

  //! @brief Define state function in single-stage IVP-DAE
  void set_function
    ( const unsigned nf, const FFVar*const FCT )
    { _nf = nf; _vFCT.clear(); _vFCT.push_back( FCT ); }

  //! @brief Define state function in multi-stage IVP-DAE with <a>ns</a> stages
  void set_function
    ( const unsigned ns, const unsigned nf, const FFVar*const FCT )
    { _nf = nf; _vFCT.clear(); for( unsigned i=0; i<ns; i++ ) _vFCT.push_back( FCT+i*_nf ); }

  //! @brief Copy DAE-IVP
  void set
    ( const BASE_DE&de )
    { _pDAG = de._pDAG;
      _nx = de._nx; _nx0 = de._nx0; _nq = de._nq; _np = de._np; _ni = de._ni; _nf = de._nf;
      _pT = de._pT; _pX = de._pX; _pP = de._pP; _pQ = de._pQ;
      delete[] _pDX; _pDX = (de._pDX && _nx? new FFVar[_nx]:0);
      delete[] _pY; delete[] _pYQ; _pY = _pYQ = 0; _ny = _nyq = 0;
      _vIC = de._vIC; _vRHS = de._vRHS; _vAE = de._vAE; _vQUAD = de._vQUAD;
      _vINV = de._vINV; _vFCT = de._vFCT; };

  //! @brief Return last successful integration time
  double final_time() const
    { return _t; }

  //! @brief Return last successful integration stage
  double final_stage() const
    { return _istg; }

  /** @} */

protected:
  //! @brief Get pointer to sensitivity/adjoint variables
  const FFVar* sensitivity() const
    { return _pY; }

  //! @brief Get pointer to quadratures
  const FFVar* quadrature() const
    { return _pQ; }

  //! @brief Get pointer to sensitivity/adjoint quadratures
  const FFVar* sensquadrature() const
    { return _pYQ; }

  //! @brief Set sensitivity/adjoint arrays
  void set_sensitivity
    ( const unsigned ny, const unsigned nyq )
    { if( _ny != ny ){ delete[] _pY; _ny = ny; _pY = new FFVar[_ny]; }
      for( unsigned iy=0; iy<_ny; iy++ )
        if( _pY[iy].dag() != _pDAG ) _pY[iy].set( _pDAG );
      if( _nyq != nyq ){ delete[] _pYQ; _nyq = nyq; _pYQ = new FFVar[_nyq]; }
      for( unsigned iyq=0; iyq<_nyq; iyq++ )
        if( _pYQ[iyq].dag() != _pDAG ) _pYQ[iyq].set( _pDAG ); }

  //! @brief Function to display intermediate results
  template<typename U> static void _print_interm
    ( const unsigned nx, const U*x, const std::string&var, std::ostream&os=std::cout );

  //! @brief Function to display intermediate results
  template<typename U> static void _print_interm
    ( const double t, const unsigned nx, const U*x, const std::string&var,
      std::ostream&os=std::cout );

  //! @brief Function to display intermediate results
  template<typename U, typename V> static void _print_interm
    ( const double t, const unsigned nx, const U*x, const V&r,
      const std::string&var, std::ostream&os=std::cout );

  //! @brief Private methods to block default compiler methods
  BASE_DE(const BASE_DE&);
  BASE_DE& operator=(const BASE_DE&);
};

template <typename U> inline void
BASE_DE::_print_interm
( const double t, const unsigned nx, const U*x, const std::string&var,
  std::ostream&os )
{
  os << " @t = " << std::scientific << std::setprecision(4)
                 << std::left << t << " :" << std::endl;
  _print_interm( nx, x, var, os );
  return;
}

template <typename U, typename V> inline void
BASE_DE::_print_interm
( const double t, const unsigned nx, const U*x, const V&r,
  const std::string&var, std::ostream&os )
{
  os << " @t = " << std::scientific << std::setprecision(4)
                 << std::left << t << " :" << std::endl;
  _print_interm( nx, x, var, os );
  os << " " << "R" << var.c_str() << " =" << r << std::endl;
  return;
}

template <typename U> inline void
BASE_DE::_print_interm
( const unsigned nx, const U*x, const std::string&var, std::ostream&os )
{
  if( !x ) return;
  for( unsigned ix=0; ix<nx; ix++ )
    os << " " << var.c_str() << "[" << ix << "] = " << x[ix] << std::endl;
  return;
}

} // end namescape mc

#endif

