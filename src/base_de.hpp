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

  //! @brief max number of time stages
  unsigned _nsmax;

  //! @brief pointer to stage times (size _nsmax+1)
  std::vector<double> _dT;

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

  //! @brief map between DAG index and state index
  std::map<int,unsigned> _ndxX;

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

  //! @brief map between DAG index and quadrature index
  std::map<int,unsigned> _ndxQ;

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

  //! @brief state parametric dependency 
  std::vector<FFDep> _depX;

  //! @brief quadrature parametric dependency 
  std::vector<FFDep> _depQ;

  //! @brief function parametric dependency 
  std::vector<FFDep> _depF;

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
    : _pDAG(0), _nsmax(0), _dT(), _pT(0), _nd(0), _na(0), _pDX(0),
      _nx(0), _pX(0), _nx0(0), _ny(0), _pY(0), _np(0), _pP(0),
      _nq(0), _pQ(0), _nyq(0), _pYQ(0), _ni(0), _nf(0)
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

  //! @brief Set state
  void set_state
    ( const unsigned nx, const FFVar*pX, const FFVar*pDX=0 )
    { if( nx != _nx){ delete[] _pDX; _pDX = (nx?new FFVar[nx]:0); }
      _nx = nx; _pX = pX; _ndxX.clear();
      for( unsigned ix=0; ix<_nx; ix++ ){
        //std::cout << _pX[ix].id().second << std::endl;
        _ndxX.insert( std::make_pair( _pX[ix].id().second, ix ) );
        //_ndxX[_pX[ix].id().second] = ix;
        if( pDX ) _pDX[ix] = pDX[ix];
        else if( _pDX[ix].dag() != _pDAG ) _pDX[ix].set( _pDAG );
      } }

  //! @brief Return size of parameter
  unsigned np() const
    { return _np; }

  //! @brief Get pointer to parameter
  const FFVar* parameter() const
    { return _pP; }

  //! @brief Set parameter
  void set_parameter
    ( const unsigned np, const FFVar*pP )
    { _np = np; _pP = pP; }

  //! @brief Return max number of time stages
  unsigned nsmax() const
    { return _nsmax; }

  //! @brief Get pointer to stage times
  const double* stage() const
    { return _dT.data(); }

  //! @brief Get pointer to time
  const FFVar* time() const
    { return _pT; }

  //! @brief Set time
  void set_time
    ( const unsigned ns, const double*dT, const FFVar*pT=0 )
    { _nsmax = ns; _dT.assign( dT, dT+ns+1 ); _pT = pT; }

  //! @brief Set time
  void set_time
    ( const double t0, const double tf, const FFVar*pT=0 )
    { _nsmax = 1;
      _dT.clear();
      _dT.push_back( t0 );
      _dT.push_back( tf );
      //_dT.resize(2);
      //_dT[0] = t0;
      //_dT[1] = tf;
      _pT = pT; }

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
    { _nq = nq; _pQ = pQ; _ndxQ.clear();
      for( unsigned iq=0; iq<_nq; iq++ )
        _ndxQ.insert( std::make_pair( _pQ[iq].id().second, iq ) );
        //_ndxQ[_pQ[iq].id().second] = iq;
      _vQUAD.clear(); _vQUAD.push_back( QUAD ); }

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

  //! @brief Reset state functions
  void reset_function
    ()
    { _nf = 0; _vFCT.clear(); }

  //! @brief Copy DAE-IVP
  void set
    ( const BASE_DE&de )
    { _pDAG = de._pDAG;
      _nsmax = de._nsmax; _nx = de._nx; _nx0 = de._nx0; _nq = de._nq; _np = de._np; _ni = de._ni; _nf = de._nf;
      _dT = de._dT; _pT = de._pT; _pX = de._pX; _pP = de._pP; _pQ = de._pQ;
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

  //! @brief Get pointer to initial/transition value function
  const FFVar* _pIC
    ( const unsigned is=0 ) const;

  //! @brief Get pointer to right-hand-side function
  const FFVar* _pRHS
    ( const unsigned is ) const;

  //! @brief Get pointer to quadrature function
  const FFVar* _pQUAD
    ( const unsigned is ) const;

  //! @brief Set state/quadrature dependencies w.r.t. parameters
  bool set_depend
    ( const unsigned ns );

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

inline const FFVar*
BASE_DE::_pIC
( const unsigned is ) const
{
  return( is<_vIC.size()? _vIC.at( is ): 0 );
}

inline const FFVar*
BASE_DE::_pRHS
( const unsigned is ) const
{
  return( is<_vRHS.size()? _vRHS.at( is ): 0 );
}

inline const FFVar*
BASE_DE::_pQUAD
( const unsigned is ) const
{
  return( is<_vQUAD.size()? _vQUAD.at( is ): 0 );
}

inline bool
BASE_DE::set_depend
( const unsigned ns )
{
  if( _nd != _nx || _nx0 != _nx ) return false; // ODE systems only

  // Generate parametric dependence of state
  _depX.resize( _nx*ns );
  _depQ.resize( _nq*ns );
  _depF.resize( _nf );

  // Intermediates
  std::vector<FFVar> vVAR( _pP, _pP+_np );
  vVAR.push_back( _pT? *_pT: 0. );
  vVAR.insert( vVAR.end(), _pX, _pX+_nx );
  vVAR.insert( vVAR.end(), _pQ, _pQ+_nq );
  const unsigned nVAR = vVAR.size();
  std::vector<FFDep> depVAR( nVAR, 0 );
  for( unsigned ip=0; ip<_np; ++ip ) depVAR[ip] = _pP[ip].dep();

  // Initial/transition condition
  const FFVar* pIC = _pIC();
  if( !pIC ) return false;
  _pDAG->eval( _nx, pIC, _depX.data(), _np+1, vVAR.data(), depVAR.data() ); 
  for( unsigned ix=0; ix<_nx; ix++ ){
    depVAR[_np+1+ix] = _depX[ix];
#ifdef MC__DEBUG__BASE_DE
    std::cout << "X[0][" << ix << "]: " << _depX[ix] << std::endl;
#endif
  }

  for( unsigned is=0; is<ns; ++is ){

    const unsigned _pos_ic = ( _vIC.size()>=ns? is:0 );
    if( _pos_ic ){
      const FFVar* pIC = _pIC(is);
      _pDAG->eval( _nx, pIC, _depX.data()+_nx*is, _np+1+_nx, vVAR.data(), depVAR.data() ); 
      for( unsigned ix=0; ix<_nx; ix++ ){
        depVAR[_np+1+ix] = _depX[ix];
#ifdef MC__DEBUG__BASE_DE
        std::cout << "X[" << is+1 << "][" << ix << "]: " << _depX[_nx*is+ix] << std::endl;
#endif
      }
    }

    const unsigned pos_rhs  = ( _vRHS.size()<=1?  0: is );
    const FFVar* pRHS  = _pRHS(pos_rhs);
    if( !pRHS ) return false;
    _pDAG->eval( _nx, pRHS, _depX.data()+_nx*is, _np+1+_nx, vVAR.data(), depVAR.data() ); 
    for( unsigned ix=0; ix<_nx; ix++ ){
#ifdef MC__DEBUG__BASE_DE
      std::cout << "X[" << is+1 << "][" << ix << "]: " << _depX[_nx*is+ix] << std::endl;
#endif
      depVAR[_np+1+ix] = _depX[_nx*is+ix];
    }
    bool iterate = true;
    while( iterate ){
      iterate = false;
      _pDAG->eval( _nx, pRHS, _depX.data()+_nx*is, _np+1+_nx, vVAR.data(), depVAR.data() ); 
      for( unsigned ix=0; ix<_nx; ix++ ){
#ifdef MC__DEBUG__BASE_DE
        std::cout << "X[" << is+1 << "][" << ix << "]: " << _depX[_nx*is+ix] << std::endl;
#endif
        if( depVAR[_np+1+ix] == _depX[_nx*is+ix] ) continue;
        iterate = true;
        depVAR[_np+1+ix] = _depX[_nx*is+ix];
      }
    }

    // Quadrature variables
    if( _nq ){
      const unsigned pos_quad = ( _vQUAD.size()<=1? 0: is );
      const FFVar* pQUAD = _pQUAD(pos_quad);
      if( !pQUAD ) return false;
      _pDAG->eval( _nq, pQUAD, _depQ.data()+_nq*is, _np+1+_nx, vVAR.data(), depVAR.data() ); 
      for( unsigned iq=0; iq<_nq; iq++ ){
#ifdef MC__DEBUG__BASE_DE
        std::cout << "Q[" << is+1 << "][" << iq << "]: " << _depQ[_nq*is+iq] << std::endl;
#endif
        depVAR[_np+1+_nx+iq] = _depQ[_nq*is+iq];
      }
    }

    // Function stage contribution
    if( _nf ){
      const FFVar* pFCT = _vFCT.at( is );
      _pDAG->eval( _nf, pFCT, _depF.data(), nVAR, vVAR.data(), depVAR.data(), true );
#ifdef MC__DEBUG__BASE_DE
      for( unsigned ic=0; ic<_nf; ic++ )
        std::cout << "G[" << is+1 << "][" << ic << "]: " << _depF[ic] << std::endl;
#endif
    }
  }

  return true;
}

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

