// Copyright (C) Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

#ifndef CRONOS__BASE_DE_HPP
#define CRONOS__BASE_DE_HPP

#undef  CRONOS__DEBUG__BASE_DE

#include <assert.h>
#include "ffdep.hpp"
#include "ffunc.hpp"

namespace mc
{
//! @brief C++ base class for defining of parametric differential-algebraic equations
////////////////////////////////////////////////////////////////////////
//! mc::BASE_DE is a C++ base class for defining the variables,
//! parameters and functions participating in parametric differential-
//! algebraic equations (DAEs)
////////////////////////////////////////////////////////////////////////
class BASE_DE
{
protected:
  //! @brief pointer to DAG of equation
  FFGraph* _dag;

  //! @brief max size of time stages
  size_t _nsmax;

  //! @brief pointer to stage times (size _nsmax+1)
  std::vector<double> _dT;

  //! @brief pointer to time/independent variable
  std::vector<FFVar> _vT;

  //! @brief size of differential states
  size_t _nd;

  //! @brief size of algebraic states
  size_t _na;

  //! @brief pointer to state time derivatives
  std::vector<FFVar> _vDX;

  //! @brief size of states
  size_t _nx;

  //! @brief vector of state variables - first _nd are differential states, next _na algebraic states
  std::vector<FFVar> _vX;

  //! @brief map between DAG index and state index
  std::map<int,size_t> _ndxX;

  //! @brief size of initial states
  size_t _nx0;

//  //! @brief size of sensitivity/adjoint variables
//  size_t _ny;

//  //! @brief pointer to sensitivity/adjoint variables
//  std::vector<FFVar> _vY;

  //! @brief size of constants
  size_t _nc;

  //! @brief pointer to parameters
  std::vector<FFVar> _vC;

  //! @brief size of parameters
  size_t _np;

  //! @brief pointer to parameters
  std::vector<FFVar> _vP;

  //! @brief size of quadratures
  size_t _nq;

  //! @brief pointer to quadrature variables
  std::vector<FFVar> _vQ;

  //! @brief map between DAG index and quadrature index
  std::map<int,size_t> _ndxQ;

  //! @brief size of bound multipliers
  size_t _nbm;

  //! @brief pointer to lower bound multipliers
  std::vector<FFVar>  _vML;

  //! @brief pointer to upper bound multipliers
  std::vector<FFVar>  _vMU;

  //! @brief vector of const pointers to initial value in each stage of IVP-DAE system
  std::vector<std::vector<FFVar>> _vIC;

  //! @brief vector of const pointers to RHS of differential equations in each stage of IVP-DAE system
  std::vector<std::vector<FFVar>> _vDE;

  //! @brief vector of const pointers to algebraic equations in each stage of IVP-DAE system
  std::vector<std::vector<FFVar>> _vAE;

  //! @brief vector of const pointers to integrand of quadratures in each stage of IVP-DAE system
  std::vector<std::vector<FFVar>> _vQUAD;

  //! @brief size of invariants
  size_t _ni;

  //! @brief vector of const pointers to invariant equations in each stage of IVP-DAE system
  std::vector<std::vector<FFVar>> _vINV;

  //! @brief size of state functionals
  size_t _nf;

  //! @brief vector of const pointers to functionals (assumes additive contributions) in each stage of IVP-DAE system
  std::vector<std::vector<FFVar>> _vFCT;

  //! @brief state parametric dependency 
  std::vector<FFDep> _depX;

  //! @brief quadrature parametric dependency 
  std::vector<FFDep> _depQ;

  //! @brief function parametric dependency 
  std::vector<FFDep> _depF;

  //! @brief state parametric dependency 
  size_t _nnzjac;

  //! @brief current time
  double _t;

  //! @brief current stage
  size_t _istg;

public:
  //! @brief Class constructor
  BASE_DE
    ()
    : _dag(nullptr), _nsmax(0), _nd(0), _na(0), _nx(0), _nx0(0),
      _nc(0), _np(0), _nq(0), _nbm(0), _ni(0), _nf(0), _nnzjac(0)
    {}

  //! @brief Class destructor
  virtual ~BASE_DE
    ()
    {}

  //! @brief Integrator status
  enum STATUS{
     NORMAL=0,	//!< Normal execution
     FAILURE,	//!< Integration breakdown
     FATAL	//!< Interruption due to errors in third-party libraries
  };

  //! @brief Get pointer to DAG
  FFGraph* dag
    ()
    const
    { return _dag; }

  //! @brief Set pointer to DAG
  void set_dag
    ( FFGraph* dag )
    { _dag = dag; }

  //! @brief Get size of state
  size_t nx
    ()
    const
    { return _nx; }

  //! @brief Get size of differential state
  size_t nd
    ()
    const
    { return _nd; }

  //! @brief Get size of algebraic state
  size_t na
    ()
    const
    { return _na; }

  //! @brief Get pointer to state
  std::vector<FFVar> const& var_state
    ()
    const
    { return _vX; }

  //! @brief Get pointer to state
  std::vector<FFVar> const& var_state_deriv
    ()
    const
    { return _vDX; }

  //! @brief Set state
  void set_state
    ( std::vector<FFVar> const& vX, std::vector<FFVar> const& vDX=std::vector<FFVar>() )
    { _nx = vX.size();
      _vX = vX;
      _vDX.resize( _nx );
      _ndxX.clear();
      for( size_t ix=0; ix<_nx; ix++ ){
        _ndxX.insert( std::make_pair( _vX[ix].id().second, ix ) );
        if( vDX.size() > ix ) _vDX[ix] = vDX[ix];
        else                  _vDX[ix].set( _dag );
      } }

  //! @brief Return size of constant
  size_t nc
    ()
    const
    { return _nc; }

  //! @brief Get pointer to constant
  std::vector<FFVar> const& var_constant
    ()
    const
    { return _vC; }

  //! @brief Set constant
  void set_constant
    ( std::vector<FFVar> const& vC )
    { _nc = vC.size(); _vC = vC; }

  //! @brief Return size of parameter
  size_t np
    ()
    const
    { return _np; }

  //! @brief Get pointer to parameter
  std::vector<FFVar> const& var_parameter
    ()
    const
    { return _vP; }

  //! @brief Set parameter
  void set_parameter
    ( std::vector<FFVar> const& vP )
    { _np = vP.size(); _vP = vP; }

  //! @brief Return max size of time stages
  size_t nsmax
    ()
    const
    { return _nsmax; }

  //! @brief Get pointer to stage times
  std::vector<double> const& val_stage
    ()
    const
    { return _dT; }

  //! @brief Get pointer to time
  FFVar const* var_time
    ()
    const
    { return _vT.size()? _vT.data(): nullptr; }

  //! @brief Set time
  void set_time
    ( double const& t0, double const& tf, FFVar const* pT=nullptr )
    { set_time( std::vector<double>({t0,tf}), pT ); }

  //! @brief Set time
  void set_time
    ( std::vector<double> const& dT, FFVar const* pT=nullptr )
    { _nsmax = dT.size()-1; _dT = dT;
      if( pT ) _vT.assign( pT, pT+1 );
      else     _vT.clear(); }

  //! @brief Get pointer to differential expressions
  std::vector<std::vector<FFVar>> const& eqn_differential
    ()
    const
    { return _vDE; }

  //! @brief Define differential equations in single-stage problem
  void set_differential
    ( std::vector<FFVar> const& vDE )
    { _nd = vDE.size(); _vDE.assign( { vDE } ); }

  //! @brief Define differential equations in multi-stage IVP-DAE
  void set_differential
    ( std::vector<std::vector<FFVar>> const& vDE )
    { _nd = vDE.front().size(); _vDE = vDE; }

  //! @brief Get pointer to algebraic expressions
  std::vector<std::vector<FFVar>> const& eqn_algebraic
    ()
    const
    { return _vAE; }

  //! @brief Define algebraic equations in single-stage IVP-DAE
  void set_algebraic
    ( std::vector<FFVar> const& vAE )
    { _na = vAE.size(); _vAE.assign( { vAE } ); }

  //! @brief Define algebraic equations in multi-stage IVP-DAE
  void set_algebraic
    ( std::vector<std::vector<FFVar>> const& vAE )
    { _na = vAE.front().size(); _vAE = vAE; }

  //! @brief Return size of quadrature
  size_t nq
    ()
    const
    { return _nq; }

  //! @brief Get pointer to quadratures
  std::vector<FFVar> const& var_quadrature
    ()
    const
    { return _vQ; }

  //! @brief Get quadrature expressions
  std::vector<std::vector<FFVar>> const& eqn_quadrature
    ()
    const
    { return _vQUAD; }

  //! @brief Define quadrature equations in single-stage IVP-DAE
  void set_quadrature
    ( std::vector<FFVar> const& vQUAD, std::vector<FFVar> const& vQ )
    { set_quadrature( std::vector<std::vector<FFVar>>({vQUAD}), vQ ); }

  //! @brief Define quadrature equations in multi-stage IVP-DAE with <a>ns</a> stages
  void set_quadrature
    ( std::vector<std::vector<FFVar>> const& vQUAD, std::vector<FFVar> const& vQ )
    { _nq = vQ.size();
      _vQ = vQ;
      _ndxQ.clear();
      for( size_t iq=0; iq<_nq; iq++ )
        _ndxQ.insert( std::make_pair( _vQ[iq].id().second, iq ) );
      _vQUAD = vQUAD; }

  //! @brief Return size of invariant
  size_t ni
    ()
    const
    { return _ni; }

  //! @brief Get invariant equations
  std::vector<std::vector<FFVar>> const& eqn_invariant
    ()
    const
    { return _vINV; }

  //! @brief Define invariant equations in single-stage IVP-DAE
  void set_invariant
    ( std::vector<FFVar> const& vINV )
    { _ni = vINV.size(); _vINV.assign( { vINV } ); }

  //! @brief Define initial conditions in multi-stage IVP-DAE
  void set_invariant
    ( std::vector<std::vector<FFVar>> const& vINV )
    { _ni = vINV.front().size(); _vINV = vINV; }

  //! @brief Return size of function
  size_t nx0
    ()
    const
    { return _nx0; }

  //! @brief Get initial conditions
  std::vector<std::vector<FFVar>> const& eqn_initial
    ()
    const
    { return _vIC; }

  //! @brief Define initial conditions in single-stage IVP-DAE
  void set_initial
    ( std::vector<FFVar> const& vIC )
    { _nx0 = vIC.size(); _vIC.assign( { vIC } ); }

  //! @brief Define initial conditions in multi-stage IVP-DAE
  void set_initial
    ( std::vector<std::vector<FFVar>> const& vIC )
    { _nx0 = vIC.front().size(); _vIC = vIC; }

  //! @brief Return size of function
  size_t nf
    ()
    const
    { return _nf; }

  //! @brief Get initial conditions
  std::vector<std::vector<FFVar>> const& eqn_function
    ()
    const
    { return _vFCT; }

  //! @brief Define state function in single-stage IVP-DAE
  void set_function
    ( std::vector<FFVar> const& vFCT )
    { _nf = vFCT.size(); _vFCT.assign( { vFCT } ); }

  //! @brief Define state functions in multi-stage IVP-DAE
  void set_function
    ( std::vector<std::vector<FFVar>> const& vFCT )
    { _nf = vFCT.front().size(); _vFCT = vFCT; }

  //! @brief Reset state functions
  void reset_function
    ()
    { _nf = 0; _vFCT.clear(); }

  //! @brief Copy DAE-IVP
  void set
    ( BASE_DE const& de )
    { _dag = de._dag;
      _nsmax = de._nsmax; _nx = de._nx; _nx0 = de._nx0; _nd = de._nd; _na = de._na;
      _nq = de._nq; _nc = de._nc; _np = de._np; _ni = de._ni; _nf = de._nf; _nnzjac = de._nnzjac;
      _dT = de._dT; _vT = de._vT; _vX = de._vX; _vC = de._vC; _vP = de._vP; _vQ = de._vQ; _vDX = de._vDX;
      _ndxX = de._ndxX; _ndxQ = de._ndxQ;
      _nbm = de._nbm; _vML = de._vML; _vMU = de._vMU;
      _vIC = de._vIC; _vDE = de._vDE; _vAE = de._vAE; _vQUAD = de._vQUAD;
      _vINV = de._vINV; _vFCT = de._vFCT; };

  //! @brief Return last successful integration time
  double final_time
    ()
    const
    { return _t; }

  //! @brief Return last successful integration stage
  size_t final_stage
    ()
    const
    { return _istg; }

protected:

  //! @brief Get pointer to lower bound multipliers
  std::vector<FFVar> lowerboundmultiplier
    ()
    const
    { return _vML; }

  //! @brief Get pointer to upper bound multipliers
  std::vector<FFVar> upperboundmultiplier
    ()
    const
    { return _vMU; }

  //! @brief Set bound multiplier array
  void set_boundmultiplier
    ()
    { _nbm = _np;
      _vML.resize( _np );
      for( auto& ML : _vML ) ML.set( _dag );
      _vMU.resize( _np );
      for( auto& MU : _vMU ) MU.set( _dag ); }

  //! @brief Get pointer to initial/transition value function
  FFVar const* _pIC
    ( size_t const is=0 )
    const
    { return( is<_vIC.size()? _vIC.at( is ).data(): nullptr ); }

  //! @brief Get pointer to right-hand-side function
  const FFVar* _pRHS
    ( size_t const is )
    const
    { return( is<_vDE.size()? _vDE.at( is ).data(): nullptr ); }

  //! @brief Get pointer to quadrature function
  const FFVar* _pQUAD
    ( size_t const is )
    const
    { return( is<_vQUAD.size()? _vQUAD.at( is ).data(): nullptr ); }

  //! @brief Get pointer to state function
  const FFVar* _pFCT
    ( size_t const is )
    const
    { return( is<_vFCT.size()? _vFCT.at( is ).data(): nullptr ); }
/*
  //! @brief Set state/quadrature dependencies w.r.t. parameters
  bool set_depend
    ( size_t const ns );
*/
  //! @brief Set RHS Jacobain sparsity size
  bool set_sparse
    ();

  //! @brief Function to display intermediate results
  template<typename U>
  static void _print_interm
    ( size_t const nx, U const* x, std::string const& var, std::ostream& os=std::cout );

  //! @brief Function to display intermediate results
  template<typename U>
  static void _print_interm
    ( double const& t, size_t const nx, U const* x, std::string const& var,
      std::ostream& os=std::cout );

  //! @brief Function to display intermediate results
  template<typename U, typename V>
  static void _print_interm
    ( double const& t, size_t const nx, U const* x, V const& r,
      std::string const& var, std::ostream& os=std::cout );

  //! @brief Private methods to block default compiler methods
  BASE_DE( BASE_DE const& ) = delete;
  BASE_DE& operator=( BASE_DE const& ) = delete;
};

//! @brief Display IVP-DAE
inline std::ostream&
operator <<
( std::ostream & out, BASE_DE const& IVP )
{
//  out << std::left << std::endl
//      << std::setfill('_') << std::setw(72) << "#" << std::endl << "#" << std::endl << std::setfill(' ')
//      << "#  LOCAL MIXED-INTEGER NONLINEAR OPTIMIZATION IN CANON\n"
//      << std::setfill('_') << std::setw(72) << "#" << std::endl << "#" << std::endl << std::setfill(' ');

//  // Display MINLPSLV Options
//  MINLP.options.display( out );

//  out << std::left
//      << std::setfill('_') << std::setw(72) << "#" << std::endl << std::endl << std::setfill(' ');
  return out;
}

template <typename U>
inline
void
BASE_DE::_print_interm
( double const& t, size_t const nx, U const* x, std::string const& var,
  std::ostream& os )
{
  os << " @t = " << std::scientific << std::setprecision(6)
                 << std::left << t << " :" << std::endl;
  _print_interm( nx, x, var, os );
  return;
}

template <typename U, typename V>
inline
void
BASE_DE::_print_interm
( double const& t, size_t const nx, U const* x, V const& r,
  std::string const& var, std::ostream& os )
{
  os << " @t = " << std::scientific << std::setprecision(6)
                 << std::left << t << " :" << std::endl;
  _print_interm( nx, x, var, os );
  os << " " << "R" << var << " =" << r << std::endl;
  return;
}

template <typename U>
inline
void
BASE_DE::_print_interm
( size_t const nx, U const* x, std::string const& var, std::ostream& os )
{
  if( !x ) return;
  for( size_t ix=0; ix<nx; ix++ )
    os << " " << var << "[" << ix << "] = " << x[ix] << std::endl;
  return;
}

inline
bool
BASE_DE::set_sparse
()
{
  _nnzjac = 0;
  if( _nd != _nx || _nx0 != _nx ) return false; // ODE systems only

  // Intermediates
  std::vector<FFVar> vVAR = _vX;
  vVAR.insert( vVAR.end(), _vP.cbegin(), _vP.cend() );
  vVAR.push_back( _vT.size()? _vT.front(): 0. );
  size_t const nVAR = vVAR.size();
  std::vector<FFDep> depSTA( nVAR, 0 ), depRHS( _nx );
  for( size_t ix=0; ix<_nx; ++ix ) depSTA[ix].indep( ix );

  // RHS dependencies
  for( size_t is=0; is<_nsmax; ++is ){
    size_t const pos_rhs  = ( _vDE.size()<=1?  0: is );
    FFVar const* pRHS  = _pRHS( pos_rhs );
    if( !pRHS ) return false;
    if( !_nc ) _dag->eval( _nx, pRHS, depRHS.data(), nVAR, vVAR.data(), depSTA.data() ); 
    else       _dag->eval( _nx, pRHS, depRHS.data(), nVAR, vVAR.data(), depSTA.data(),
                           _nc, _vC.data(), std::vector<FFDep>(_nc, 0 ).data() ); 
    size_t nnz = 0;
    for( size_t ix=0; ix<_nx; ix++ ){
#ifdef CRONOS__DEBUG__BASE_DE
      std::cout << "RHS[" << is << "][" << ix << "]: " << depRHS[ix] << std::endl;
#endif
      nnz += depRHS[ix].dep().size();
    }

    if( _nnzjac < nnz ) _nnzjac = nnz;
#ifdef CRONOS__DEBUG__BASE_DE
    std::cout << "NNZ: " << nnz << "   MAX: " << _nnzjac << std::endl;
#endif
  }

  return true;
}
/*
inline
bool
BASE_DE::set_depend
( size_t const ns )
{
  if( _nd != _nx || _nx0 != _nx ) return false; // ODE systems only

  // Generate parametric dependence of state
  _depX.resize( _nx*ns );
  _depQ.resize( _nq*ns );
  _depF.resize( _nf );

  // Intermediates
  std::vector<FFVar> vVAR = _vP;
  vVAR.push_back( _vT.size()? _vT.front(): 0. );
  vVAR.insert( vVAR.end(), _vX.cbegin(), _vX.cend() );
  vVAR.insert( vVAR.end(), _vQ.cbegin(), _vQ.cend() );
  size_t const nVAR = vVAR.size();
  std::vector<FFDep> depVAR( nVAR, 0 );
  for( size_t ip=0; ip<_np; ++ip ) depVAR[ip].indep( ip );

  // Initial/transition condition
  FFVar const* pIC = _pIC();
  if( !pIC ) return false;
  _dag->eval( _nx, pIC, _depX.data(), _np+1, vVAR.data(), depVAR.data() ); 
  for( size_t ix=0; ix<_nx; ix++ ){
    depVAR[_np+1+ix] = _depX[ix];
#ifdef CRONOS__DEBUG__BASE_DE
    std::cout << "X[0][" << ix << "]: " << _depX[ix] << std::endl;
#endif
  }

  for( size_t is=0; is<ns; ++is ){

    size_t const _pos_ic = ( _vIC.size()>=ns? is:0 );
    if( _pos_ic ){
      FFVar const* pIC = _pIC( is );
      _dag->eval( _nx, pIC, _depX.data()+_nx*is, _np+1+_nx, vVAR.data(), depVAR.data() ); 
      for( size_t ix=0; ix<_nx; ix++ ){
        depVAR[_np+1+ix] = _depX[ix];
#ifdef CRONOS__DEBUG__BASE_DE
        std::cout << "X[" << is+1 << "][" << ix << "]: " << _depX[_nx*is+ix] << std::endl;
#endif
      }
    }

    size_t const pos_rhs  = ( _vDE.size()<=1?  0: is );
    FFVar const* pRHS  = _pRHS( pos_rhs );
    if( !pRHS ) return false;
#ifdef CRONOS__DEBUG__BASE_DE
    for( size_t ix=0; ix<_nx; ix++ )
      std::cout << "X[" << is+1 << "][" << ix << "]: " << depVAR[_np+1+ix] << std::endl;
#endif
    _dag->eval( _nx, pRHS, _depX.data()+_nx*is, _np+1+_nx, vVAR.data(), depVAR.data() ); 
    for( size_t ix=0; ix<_nx; ix++ ){
      depVAR[_np+1+ix] += _depX[_nx*is+ix];
#ifdef CRONOS__DEBUG__BASE_DE
      std::cout << "X[" << is+1 << "][" << ix << "]: " << depVAR[_np+1+ix] << std::endl;
#endif
    }

    bool iterate = true;
    while( iterate ){
      iterate = false;
      _dag->eval( _nx, pRHS, _depX.data()+_nx*is, _np+1+_nx, vVAR.data(), depVAR.data() ); 
      for( size_t ix=0; ix<_nx; ix++ ){
        if( depVAR[_np+1+ix] == depVAR[_np+1+ix]+_depX[_nx*is+ix] ) continue;
        iterate = true;
        depVAR[_np+1+ix] += _depX[_nx*is+ix];
#ifdef CRONOS__DEBUG__BASE_DE
        std::cout << "X[" << is+1 << "][" << ix << "]: " << depVAR[_np+1+ix] << std::endl;
#endif
      }
    }

    // Quadrature variables
    if( _nq ){
      size_t const pos_quad = ( _vQUAD.size()<=1? 0: is );
      FFVar const* pQUAD = _pQUAD( pos_quad );
      if( !pQUAD ) return false;
#ifdef CRONOS__DEBUG__BASE_DE
      _dag->output( _dag->subgraph( _nq, pQUAD ) );
#endif
      _dag->eval( _nq, pQUAD, _depQ.data()+_nq*is, _np+1+_nx, vVAR.data(), depVAR.data() ); 
      for( size_t iq=0; iq<_nq; iq++ ){
#ifdef CRONOS__DEBUG__BASE_DE
        std::cout << "Q[" << is+1 << "][" << iq << "]: " << _depQ[_nq*is+iq] << std::endl;
#endif
        depVAR[_np+1+_nx+iq] += _depQ[_nq*is+iq];
      }
    }

    // Function stage contribution
    size_t const pos_fct = ( _vFCT.size()<=1? 0: is );
    FFVar const* pFCT = _pFCT( pos_fct );
    if( _nf && ((is && pos_fct == is) || (!is && ){

      _dag->eval( _nf, pFCT, _depF.data(), nVAR, vVAR.data(), depVAR.data(), true );
#ifdef CRONOS__DEBUG__BASE_DE
      for( size_t ic=0; ic<_nf; ic++ )
        std::cout << "G[" << is+1 << "][" << ic << "]: " << _depF[ic] << std::endl;
#endif
    }
  }

  return true;
}
*/
} // end namescape mc

#endif

