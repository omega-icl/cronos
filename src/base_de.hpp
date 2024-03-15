// Copyright (C) Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

#ifndef MC__BASE_DE_HPP
#define MC__BASE_DE_HPP

#undef  MC__DEBUG__BASE_DE

#include <assert.h>
#include "ffunc.hpp"

namespace mc
{
//! @brief C++ base class for defining of parametric differential-algebraic equations
////////////////////////////////////////////////////////////////////////
//! mc::BASE_DE is a C++ base class for defining the variables,
//! parameters and functions participating in parametric differential-
//! algebraic equations (DAEs)
////////////////////////////////////////////////////////////////////////
template <typename... ExtOps>
class BASE_DE
{
protected:
  //! @brief pointer to DAG of equation
  FFGraph<ExtOps...>* _dag;

  //! @brief max size of time stages
  unsigned _nsmax;

  //! @brief pointer to stage times (size _nsmax+1)
  std::vector<double> _dT;

  //! @brief pointer to time/independent variable
  const FFVar* _pT;

  //! @brief size of differential states
  unsigned _nd;

  //! @brief size of algebraic states
  unsigned _na;

  //! @brief pointer to state time derivatives
  FFVar* _pDX;

  //! @brief size of states
  unsigned _nx;

  //! @brief pointer to state variables - first _nd are differential states, next _na algebraic states
  const FFVar* _pX;

  //! @brief map between DAG index and state index
  std::map<int,unsigned> _ndxX;

  //! @brief size of initial states
  unsigned _nx0;

//  //! @brief size of sensitivity/adjoint variables
//  unsigned _ny;

//  //! @brief pointer to sensitivity/adjoint variables
//  FFVar* _pY;

  //! @brief size of parameters
  unsigned _np;

  //! @brief pointer to parameters
  const FFVar* _pP;

  //! @brief size of quadratures
  unsigned _nq;

  //! @brief pointer to quadrature variables
  const FFVar* _pQ;

  //! @brief map between DAG index and quadrature index
  std::map<int,unsigned> _ndxQ;

  //! @brief size of bound multipliers
  unsigned _nbm;

  //! @brief pointer to lower bound multipliers
  FFVar* _pML;

  //! @brief pointer to upper bound multipliers
  FFVar* _pMU;

  //! @brief vector of const pointers to initial value in each stage of IVP-DAE system
  std::vector<const FFVar*> _vIC;

  //! @brief vector of const pointers to RHS of differential equations in each stage of IVP-DAE system
  std::vector<const FFVar*> _vRHS;

  //! @brief vector of const pointers to algebraic equations in each stage of IVP-DAE system
  std::vector<const FFVar*> _vAE;

  //! @brief vector of const pointers to integrand of quadratures in each stage of IVP-DAE system
  std::vector<const FFVar*> _vQUAD;

  //! @brief size of invariants
  unsigned _ni;

  //! @brief vector of const pointers to invariant equations in each stage of IVP-DAE system
  std::vector<const FFVar*> _vINV;

  //! @brief size of state functionals
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
  //! @brief Class constructor
  BASE_DE()
    : _dag(nullptr), _nsmax(0), _dT(), _pT(nullptr),
      _nd(0), _na(0), _pDX(nullptr), _nx(0), _pX(nullptr), _nx0(0),
      _np(0), _pP(nullptr), _nq(0), _pQ(nullptr),
      _nbm(0), _pML(nullptr), _pMU(nullptr), _ni(0), _nf(0)
    {}

  //! @brief Class destructor
  virtual ~BASE_DE()
    { delete[] _pDX; delete[] _pML; delete[] _pMU; }

  //! @brief Integrator status
  enum STATUS{
     NORMAL=0,	//!< Normal execution
     FAILURE,	//!< Integration breakdown (bounds explosion)
     FATAL	//!< Interruption due to errors in third-party libraries
  };

  //! @brief Get pointer to DAG
  FFGraph<ExtOps...>* dag()
    const
    { return _dag; }

  //! @brief Set pointer to DAG
  void set_dag
    ( FFGraph<ExtOps...>* dag )
    { _dag = dag; }

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
    ( const unsigned nx, const FFVar*pX, const FFVar*pDX=nullptr )
    { if( nx != _nx){ delete[] _pDX; _pDX = (nx?new FFVar[nx]:nullptr); }
      _nx = nx; _pX = pX; _ndxX.clear();
      for( unsigned ix=0; ix<_nx; ix++ ){
        //std::cout << _pX[ix].id().second << std::endl;
        _ndxX.insert( std::make_pair( _pX[ix].id().second, ix ) );
        if( pDX ) _pDX[ix] = pDX[ix];
        else if( _pDX[ix].dag() != _dag ) _pDX[ix].set( _dag );
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

  //! @brief Return max size of time stages
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
    ( const unsigned ns, const double*dT, const FFVar*pT=nullptr )
    { _nsmax = ns; _dT.assign( dT, dT+ns+1 ); _pT = pT; }

  //! @brief Set time
  void set_time
    ( const double t0, const double tf, const FFVar*pT=nullptr )
    { _nsmax = 1;
      _dT.clear();
      _dT.push_back( t0 );
      _dT.push_back( tf );
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

  //! @brief Return size of quadrature
  unsigned nq() const
    { return _nq; }

  //! @brief Get pointer to quadratures
  const FFVar* quadrature() const
    { return _pQ; }

  //! @brief Define quadrature equations in single-stage IVP-DAE
  void set_quadrature
    ( const unsigned nq, const FFVar*const QUAD, const FFVar*pQ )
    { _nq = nq; _pQ = pQ; _ndxQ.clear();
      for( unsigned iq=0; iq<_nq; iq++ )
        _ndxQ.insert( std::make_pair( _pQ[iq].id().second, iq ) );
      _vQUAD.clear(); _vQUAD.push_back( QUAD ); }

  //! @brief Define quadrature equations in multi-stage IVP-DAE with <a>ns</a> stages
  void set_quadrature
    ( const unsigned ns, const unsigned nq, const FFVar*const QUAD, const FFVar*pQ )
    { _nq = nq; _pQ = pQ; _ndxQ.clear();
      for( unsigned iq=0; iq<_nq; iq++ )
        _ndxQ.insert( std::make_pair( _pQ[iq].id().second, iq ) );
      _vQUAD.clear(); for( unsigned i=0; i<ns; i++ ) _vQUAD.push_back( QUAD+i*_nq ); }

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

  //! @brief Return size of function
  unsigned nf() const
    { return _nf; }

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
    ( BASE_DE<ExtOps...> const& de )
    { _dag = de._dag;
      _nsmax = de._nsmax; _nx = de._nx; _nx0 = de._nx0; _nd = de._nd; _na = de._na;
      _nq = de._nq; _np = de._np; _ni = de._ni; _nf = de._nf;
      _dT = de._dT; _pT = de._pT; _pX = de._pX; _pP = de._pP; _pQ = de._pQ;
      _ndxX = de._ndxX; _ndxQ = de._ndxQ;
      delete[] _pDX; _pDX = (de._pDX && _nx? new FFVar[_nx]: nullptr);
      delete[] _pML; delete[] _pMU; _pML = _pMU = nullptr;
      _nbm = 0;
      _vIC = de._vIC; _vRHS = de._vRHS; _vAE = de._vAE; _vQUAD = de._vQUAD;
      _vINV = de._vINV; _vFCT = de._vFCT; };

  //! @brief Return last successful integration time
  double final_time() const
    { return _t; }

  //! @brief Return last successful integration stage
  double final_stage() const
    { return _istg; }

protected:

  //! @brief Get pointer to lower bound multipliers
  const FFVar* lowerboundmultiplier() const
    { return _pML; }

  //! @brief Get pointer to upper bound multipliers
  const FFVar* upperboundmultiplier() const
    { return _pMU; }

  //! @brief Set bound multiplier array
  void set_boundmultiplier
    ()
    { if( _nbm != _np ){
        delete[] _pML; delete[] _pMU;
        _nbm = _np; _pML = new FFVar[_nbm]; _pMU = new FFVar[_nbm];
      }
      for( unsigned ibm=0; ibm<_nbm; ibm++ ){
        if( _pML[ibm].dag() != _dag ) _pML[ibm].set( _dag );
        if( _pMU[ibm].dag() != _dag ) _pMU[ibm].set( _dag );
      } }

  //! @brief Get pointer to initial/transition value function
  const FFVar* _pIC
    ( const unsigned is=0 ) const
    { return( is<_vIC.size()? _vIC.at( is ): nullptr ); }

  //! @brief Get pointer to right-hand-side function
  const FFVar* _pRHS
    ( const unsigned is ) const
    { return( is<_vRHS.size()? _vRHS.at( is ): nullptr ); }

  //! @brief Get pointer to quadrature function
  const FFVar* _pQUAD
    ( const unsigned is ) const
    { return( is<_vQUAD.size()? _vQUAD.at( is ): nullptr ); }

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
  BASE_DE( BASE_DE<ExtOps...> const& ) = delete;
  BASE_DE<ExtOps...>& operator=( BASE_DE<ExtOps...> const& ) = delete;
};

template <typename... ExtOps>
template <typename U>
inline
void
BASE_DE<ExtOps...>::_print_interm
( const double t, const unsigned nx, const U*x, const std::string&var,
  std::ostream&os )
{
  os << " @t = " << std::scientific << std::setprecision(4)
                 << std::left << t << " :" << std::endl;
  _print_interm( nx, x, var, os );
  return;
}

template <typename... ExtOps>
template <typename U, typename V>
inline
void
BASE_DE<ExtOps...>::_print_interm
( const double t, const unsigned nx, const U*x, const V&r,
  const std::string&var, std::ostream&os )
{
  os << " @t = " << std::scientific << std::setprecision(4)
                 << std::left << t << " :" << std::endl;
  _print_interm( nx, x, var, os );
  os << " " << "R" << var << " =" << r << std::endl;
  return;
}

template <typename... ExtOps>
template <typename U>
inline
void
BASE_DE<ExtOps...>::_print_interm
( const unsigned nx, const U*x, const std::string&var, std::ostream&os )
{
  if( !x ) return;
  for( unsigned ix=0; ix<nx; ix++ )
    os << " " << var << "[" << ix << "] = " << x[ix] << std::endl;
  return;
}

template <typename... ExtOps>
inline
bool
BASE_DE<ExtOps...>::set_depend
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
  for( unsigned ip=0; ip<_np; ++ip ) depVAR[ip].indep( ip );//_pP[ip].id().second );

  // Initial/transition condition
  const FFVar* pIC = _pIC();
  if( !pIC ) return false;
  _dag->eval( _nx, pIC, _depX.data(), _np+1, vVAR.data(), depVAR.data() ); 
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
      _dag->eval( _nx, pIC, _depX.data()+_nx*is, _np+1+_nx, vVAR.data(), depVAR.data() ); 
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
#ifdef MC__DEBUG__BASE_DE
    for( unsigned ix=0; ix<_nx; ix++ )
      std::cout << "X[" << is+1 << "][" << ix << "]: " << depVAR[_np+1+ix] << std::endl;
#endif
    _dag->eval( _nx, pRHS, _depX.data()+_nx*is, _np+1+_nx, vVAR.data(), depVAR.data() ); 
    for( unsigned ix=0; ix<_nx; ix++ ){
      depVAR[_np+1+ix] += _depX[_nx*is+ix];
#ifdef MC__DEBUG__BASE_DE
      std::cout << "X[" << is+1 << "][" << ix << "]: " << depVAR[_np+1+ix] << std::endl;
#endif
    }
    bool iterate = true;
    while( iterate ){
      iterate = false;
      _dag->eval( _nx, pRHS, _depX.data()+_nx*is, _np+1+_nx, vVAR.data(), depVAR.data() ); 
      for( unsigned ix=0; ix<_nx; ix++ ){
        if( depVAR[_np+1+ix] == depVAR[_np+1+ix]+_depX[_nx*is+ix] ) continue;
        iterate = true;
        depVAR[_np+1+ix] += _depX[_nx*is+ix];
#ifdef MC__DEBUG__BASE_DE
        std::cout << "X[" << is+1 << "][" << ix << "]: " << depVAR[_np+1+ix] << std::endl;
#endif
      }
    }

    // Quadrature variables
    if( _nq ){
      const unsigned pos_quad = ( _vQUAD.size()<=1? 0: is );
      const FFVar* pQUAD = _pQUAD(pos_quad);
      if( !pQUAD ) return false;
#ifdef MC__DEBUG__BASE_DE
      _dag->output( _dag->subgraph( _nq, pQUAD ) );
#endif
      _dag->eval( _nq, pQUAD, _depQ.data()+_nq*is, _np+1+_nx, vVAR.data(), depVAR.data() ); 
      for( unsigned iq=0; iq<_nq; iq++ ){
#ifdef MC__DEBUG__BASE_DE
        std::cout << "Q[" << is+1 << "][" << iq << "]: " << _depQ[_nq*is+iq] << std::endl;
#endif
        depVAR[_np+1+_nx+iq] += _depQ[_nq*is+iq];
      }
    }

    // Function stage contribution
    if( _nf ){
      const FFVar* pFCT = _vFCT.at( is );
      _dag->eval( _nf, pFCT, _depF.data(), nVAR, vVAR.data(), depVAR.data(), true );
#ifdef MC__DEBUG__BASE_DE
      for( unsigned ic=0; ic<_nf; ic++ )
        std::cout << "G[" << is+1 << "][" << ic << "]: " << _depF[ic] << std::endl;
#endif
    }
  }

  return true;
}

} // end namescape mc

#endif

