// Copyright (C) Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

#ifndef CRONOS__ODESLV_BASE_HPP
#define CRONOS__ODESLV_BASE_HPP

#undef  CRONOS__ODESLV_BASE_DEBUG

#include <stdexcept>
#include <cassert>
#include <cmath>
#include <fstream>
#include <vector>
#include <list>
#include <sys/time.h>

#include "base_de.hpp"

//#ifdef CRONOS__MBDOE_SETUP_DEBUG
// #include "ffexpr.hpp"
//#endif

namespace mc
{
//! @brief C++ base class for computing solutions of parametric ODEs using continuous-time real-valued integration.
////////////////////////////////////////////////////////////////////////
//! mc::ODESLV_BASE is a C++ base class for computing solutions of
//! parametric ordinary differential equation (ODEs) using
//! continuous-time real-valued integration.
////////////////////////////////////////////////////////////////////////
class ODESLV_BASE
: public virtual BASE_DE
{
 public:
  /** @defgroup ODESLV Continuous-time real-valued integration of parametric ODEs
   *  @{
   */
  //! @brief Default constructor
  ODESLV_BASE();

  //! @brief Virtual destructor
  virtual ~ODESLV_BASE();

  //! @brief Integration results at a given time instant
  struct Results
  {
    //! @brief Constructors
    Results
      ( double const& tk, unsigned const nx1, double const* x1,
        unsigned const nx2=0, double const* x2=nullptr ):
      t( tk )
      { x.resize( nx1+nx2 );
        unsigned int ix=0;
        for( ; ix<nx1; ix++ )      x[ix] = x1[ix];
        for( ; ix<nx1+nx2; ix++ )  x[ix] = x2[ix-nx1]; }
    Results
      ( Results const& res ):
      t( res.t ), x( res.x )
      {}
    //! @brief Time point
    double t;
    //! @brief Solution point
    std::vector<double> x;
  };
  /** @} */

protected:
  using BASE_DE::_nsmax;
  using BASE_DE::_nx;
  using BASE_DE::_nx0;
  using BASE_DE::_nc;
  using BASE_DE::_np;
  using BASE_DE::_nq;
  using BASE_DE::_nf;
  using BASE_DE::_t;
  using BASE_DE::_istg;
  using BASE_DE::_nnzjac;

  //! @brief local copy of DAG
  FFGraph* _dag;

  //! @brief local copy of initial value functions
  std::vector<FFVar*> _vIC;

  //! @brief local copy of right-hand-side functions
  std::vector<FFVar*> _vRHS;

  //! @brief local copy of quadrature functions
  std::vector<FFVar*> _vQUAD;

  //! @brief local copy of output functions
  std::vector<FFVar*> _vFCT;

  //! @brief local copy of constants
  FFVar* _pC;

  //! @brief local copy of parameters
  FFVar* _pP;

  //! @brief local copy of time/independent variable
  FFVar* _pT;

  //! @brief local copy of state variables
  FFVar* _pX;

  //! @brief local copy of quadrature variables
  FFVar* _pQ;

  //! @brief Sugraph of ODE RHS function
  FFSubgraph _opRHS;

  //! @brief Subgraph of ODE RHS Jacobian
  FFSubgraph _opJAC;

  //! @brief Subgraph of quadrature RHS function
  FFSubgraph _opQUAD;

  //! @brief const pointer to RHS function in current stage of ODE system
  FFVar const* _pRHS;

  //! @brief const pointer to quadrature integrand in current stage of ODE system
  FFVar const* _pQUAD;

  //! @brief const pointer to IC function in current stage of ODE system
  FFVar const* _pIC;

  //! @brief sparse representation of RHS Jacobian in current stage of ODE system
  std::tuple< unsigned, unsigned const*, unsigned const*, FFVar const* > _pJAC;

  //! @brief sparse representation of RHS Jacobian in current stage of ODE system: i-th entry is the index in data where the first non-zero matrix entry of the i-th column is stored (length NEQ + 1), last entry is number of non-zeros
  std::vector< size_t > _pJACCOLNDX;

#if 0
  //! @brief sparse representation of RHS Jacobian in current stage of ODE system
  int* _NDXPJAC;
#endif

  //! @brief number of variables for DAG evaluation
  unsigned _nVAR;

  //! @brief number of variables for DAG evaluation (without quadratures)
  unsigned _nVAR0;

  //! @brief array of variables for DAG evaluation
  FFVar* _pVAR;

  //! @brief array of variable values for DAG evaluation
  double* _DVAR;

  //! @brief pointer to time value **DO NOT FREE**
  double* _Dt;

  //! @brief pointer to state values **DO NOT FREE**
  double* _Dx;

  //! @brief pointer to parameter values **DO NOT FREE**
  double* _Dp;

  //! @brief pointer to quadrature values **DO NOT FREE**
  double* _Dq;

  //! @brief vector to hold function values
  std::vector<double> _Df;

  //! @brief vector to hold Jacobian evaluation results
  std::vector<double> _DJAC;

  //! @brief storage vector for DAG evaluation
  std::vector<double> _DWRK;

  //! @brief Function setting up local DAG
  bool _SETUP
    ();

  //! @brief Function setting up DAG of IVP
  bool _SETUP
    ( ODESLV_BASE const& IVP );

  //! @brief Function converting integrator array to internal format
  template <typename REALTYPE>
  static void _vec2D
    ( REALTYPE const* vec, unsigned const n, double* d );

  //! @brief Function converting integrator array to internal format
  template <typename REALTYPE>
  static void _D2vec
    ( double const* d, unsigned const n, REALTYPE* vec );

  //! @brief Function to initialize state integration
  bool _INI_D_STA
    ( double const* p, double const* c );

  //! @brief Function to finalize state integration
  bool _END_D_STA
    ();

  //! @brief Function to retreive state bounds
  template <typename REALTYPE>
  void _GET_D_STA
    ( REALTYPE const* x, REALTYPE const* q );

  //! @brief Function to set state/quadrature at initial time
  bool _IC_D_SET
    ();

  //! @brief Function to initialize quaratures
  template <typename REALTYPE>
  bool _IC_D_QUAD
    ( REALTYPE* vec );

  //! @brief Function to initialize states
  template <typename REALTYPE>
  bool _IC_D_STA
    ( double const& t, REALTYPE* vec );

  //! @brief Function to reset state at intermediate time
  bool _CC_D_SET
    ( unsigned const iIC );

  //! @brief Function to reinitialize state at intermediate time
  template <typename REALTYPE>
  bool _CC_D_STA
    ( double const& t, REALTYPE* vec );

  //! @brief Function to set RHS and QUAD pointers
  bool _RHS_D_SET
    ( unsigned const iRHS, unsigned const iQUAD );

  //! @brief Function to set RHS and QUAD pointers
  bool _RHS_D_SET
    ();

  //! @brief Function to calculate the ODE RHS
  template <typename REALTYPE>
  bool _RHS_D_STA
    ( double const& t, REALTYPE const* x, REALTYPE* xdot );

  //! @brief Function to calculate the ODE quadrature
  template <typename REALTYPE>
  bool _RHS_D_QUAD
    ( double const& t, REALTYPE const* x, REALTYPE* qdot );

  //! @brief Function to calculate the ODE Jacobian
  template <typename REALTYPE, typename INDEXTYPE>
  bool _JAC_D_STA
    ( double const& t, REALTYPE const* x, REALTYPE* jac, INDEXTYPE* ptrs, INDEXTYPE* vals );

  //! @brief Function to calculate the ODE Jacobian
  template <typename REALTYPE>
  bool _JAC_D_STA
    ( double const& t, REALTYPE const* x, REALTYPE** jac );

  //! @brief Function to calculate the functions at intermediate/end point
  bool _FCT_D_STA
    ( unsigned const iFCT, double const& t );

  //! @brief Record results in file <a>bndrec</a>, with accuracy of <a>iprec</a> digits
  template <typename VRES>
  static void _record
    ( std::ofstream& ofile, VRES const& bnd, unsigned const iprec=5 );

  //! @brief Block default compiler methods
  ODESLV_BASE( ODESLV_BASE const& ) = delete;
  ODESLV_BASE& operator=( ODESLV_BASE const& ) = delete;
};

inline 
ODESLV_BASE::ODESLV_BASE
()
: BASE_DE(),
  _dag(nullptr),
  _pC(nullptr), _pP(nullptr), _pT(nullptr), _pX(nullptr), _pQ(nullptr),
  _pRHS(nullptr), _pQUAD(nullptr), _pIC(nullptr),
  _pJAC(0,nullptr,nullptr,nullptr), 
  _nVAR(0), _nVAR0(0), _pVAR(nullptr),
  _DVAR(nullptr), _Dt(nullptr), _Dx(nullptr), _Dp(nullptr), _Dq(nullptr)
{}

inline
ODESLV_BASE::~ODESLV_BASE
()
{
  /* DO NOT FREE _pRHS, _pQUAD, _pIC */
  for( auto& ic   : _vIC )   delete[] ic;
  for( auto& rhs  : _vRHS )  delete[] rhs;
  for( auto& quad : _vQUAD ) delete[] quad;
  for( auto& fct  : _vFCT )  delete[] fct;

  delete[] std::get<1>(_pJAC);  std::get<1>(_pJAC) = nullptr;
  delete[] std::get<2>(_pJAC);  std::get<2>(_pJAC) = nullptr;
  delete[] std::get<3>(_pJAC);  std::get<3>(_pJAC) = nullptr;
  delete[] _pVAR;
  delete[] _DVAR;
  delete[] _pX;
  delete[] _pQ;
  delete[] _pP;
  delete[] _pC;
  delete   _pT;
  delete   _dag;
}

inline
bool
ODESLV_BASE::_SETUP
()
{
  delete _dag; _dag = new FFGraph;
#ifdef CRONOS__ODESLV_BASE_DEBUG
  std::cout << "ODESLV_BASE:: Original DAG: " << BASE_DE::_dag << std::endl;
  std::cout << "ODESLV_BASE:: Copied DAG:   " << _dag << std::endl;
#endif

  delete[] _pC; _pC = nullptr;
  if( _nc ){
    _pC  = new FFVar[_nc];
    _dag->insert( BASE_DE::_dag, _nc, BASE_DE::_vC.data(), _pC );
  }

  delete[] _pP; _pP = nullptr;
  if( _np ){
    _pP  = new FFVar[_np];
    _dag->insert( BASE_DE::_dag, _np, BASE_DE::_vP.data(), _pP );
  }

  delete _pT; _pT = nullptr;
  if( BASE_DE::_vT.size() ){
    _pT  = new FFVar;
    _dag->insert( BASE_DE::_dag, 1, BASE_DE::_vT.data(), _pT );
  }

  delete[] _pX; _pX = nullptr;
  if( _nx ){
    _pX  = new FFVar[_nx];
    _dag->insert( BASE_DE::_dag, _nx, BASE_DE::_vX.data(), _pX );
  }

  delete[] _pQ; _pQ = nullptr;
  if( _nq ){
    _pQ  = new FFVar[_nq];
    _dag->insert( BASE_DE::_dag, _nq, BASE_DE::_vQ.data(), _pQ );
  }
  
  for( auto& ic : _vIC ) delete[] ic;
  _vIC.clear();
  _vIC.reserve( BASE_DE::_vIC.size() );
  for( auto const& ic0 : BASE_DE::_vIC ){
    FFVar* ic = new FFVar[_nx0];
    _dag->insert( BASE_DE::_dag, _nx0, ic0.data(), ic );
    _vIC.push_back( ic );
  }

  for( auto& rhs  : _vRHS )  delete[] rhs;
  _vRHS.clear();
  _vRHS.reserve( BASE_DE::_vDE.size() );
  for( auto const& rhs0 : BASE_DE::_vDE ){
#ifdef CRONOS__ODESLV_BASE_DEBUG
    //BASE_DE::_dag->output( BASE_DE::_dag->subgraph( _nx, rhs0.data() ), " - Before insert" );
    FFSubgraph sgrhs0 = BASE_DE::_dag->subgraph( _nx, rhs0.data() );
    std::vector<FFExpr> exprrhs0 = FFExpr::subgraph( BASE_DE::_dag, sgrhs0 ); 
    for( unsigned j=0; j<_nx; ++j )
        std::cout << "RHS0[" << j << "] = " << exprrhs0[j] << std::endl;
#endif
    FFVar* rhs = new FFVar[_nx];
    //for( unsigned j=0; j<_nx; ++j ){
    //  BASE_DE::_dag->output( BASE_DE::_dag->subgraph( 1, rhs0.data()+j ), " - Before insert" );
    //  _dag->insert( BASE_DE::_dag, 1, rhs0.data()+j, rhs+j );
    //}
    _dag->insert( BASE_DE::_dag, _nx, rhs0.data(), rhs );
    _vRHS.push_back( rhs );
#ifdef CRONOS__ODESLV_BASE_DEBUG
    //_dag->output( _dag->subgraph( _nx, rhs ), " - After insert" );
    FFSubgraph sgrhs = _dag->subgraph( _nx, rhs );
    std::vector<FFExpr> exprrhs = FFExpr::subgraph( _dag, sgrhs ); 
    for( unsigned j=0; j<_nx; ++j )
        std::cout << "RHS[" << j << "] = " << exprrhs[j] << std::endl;
#endif
  }
  //std::cout << *_dag;

  for( auto& quad : _vQUAD ) delete[] quad;
  _vQUAD.clear();
  _vQUAD.reserve( BASE_DE::_vQUAD.size() );
  for( auto const& quad0 : BASE_DE::_vQUAD ){
    FFVar* quad = new FFVar[_nq];
    _dag->insert( BASE_DE::_dag, _nq, quad0.data(), quad );
    _vQUAD.push_back( quad );
  }

  for( auto& fct  : _vFCT )  delete[] fct;
  _vFCT.clear();
  _vFCT.reserve( BASE_DE::_vFCT.size() );
  for( auto const& fct0 : BASE_DE::_vFCT ){
    FFVar* fct = new FFVar[_nf];
    _dag->insert( BASE_DE::_dag, _nf, fct0.data(), fct );
    _vFCT.push_back( fct );
  }

  return true;
}

inline
bool
ODESLV_BASE::_SETUP
( ODESLV_BASE const& IVP )
{
  delete _dag; _dag = new FFGraph;

  delete _pT; _pT = nullptr;
  if( IVP._pT ){
    _pT  = new FFVar;
    _dag->insert( IVP._dag, 1, IVP._pT, _pT );
  }

  delete[] _pX; _pX = nullptr;
  if( _nx ){
    _pX  = new FFVar[_nx];
    _dag->insert( IVP._dag, _nx, IVP._pX, _pX );
  }

  delete[] _pQ; _pQ = nullptr;
  if( _nq ){
    _pQ  = new FFVar[_nq];
    _dag->insert( IVP._dag, _nq, IVP._pQ, _pQ );
  }

  delete[] _pC; _pC = nullptr;
  if( _nc ){
    _pC  = new FFVar[_nc];
    _dag->insert( IVP._dag, _nc, IVP._pC, _pC );
  }

  delete[] _pP; _pP = nullptr;
  if( _np ){
    _pP  = new FFVar[_np];
    _dag->insert( IVP._dag, _np, IVP._pP, _pP );
  }
  
  for( auto& ic : _vIC ) delete[] ic;
  _vIC.clear();
  _vIC.reserve( IVP._vIC.size() );
  for( auto const& ic0 : IVP._vIC ){
    FFVar* ic = new FFVar[_nx0];
    _dag->insert( IVP._dag, _nx0, ic0, ic );
    _vIC.push_back( ic );
  }
    
  for( auto& rhs : _vRHS )  delete[] rhs;
  _vRHS.clear();
  _vRHS.reserve( IVP._vRHS.size() );
  for( auto const& rhs0 : IVP._vRHS ){
    //IVP._dag->output( IVP._dag->subgraph( 1, rhs0 ), " - Before insert" );
    FFVar* rhs = new FFVar[_nx];
    _dag->insert( IVP._dag, _nx, rhs0, rhs );
    _vRHS.push_back( rhs );
    //_dag->output( _dag->subgraph( 1, rhs ), " - After insert" );
  }
  //std::cout << *_dag;

  for( auto& quad : _vQUAD ) delete[] quad;
  _vQUAD.clear();
  _vQUAD.reserve( IVP._vQUAD.size() );
  for( auto const& quad0 : IVP._vQUAD ){
    FFVar* quad = new FFVar[_nq];
    _dag->insert( IVP._dag, _nq, quad0, quad );
    _vQUAD.push_back( quad );
  }

  for( auto& fct : _vFCT )  delete[] fct;
  _vFCT.clear();
  _vFCT.reserve( IVP._vFCT.size() );
  for( auto const& fct0 : IVP._vFCT ){
    FFVar* fct = new FFVar[_nf];
    _dag->insert( IVP._dag, _nf, fct0, fct );
    _vFCT.push_back( fct );
  }

  return true;
}

template <typename VRES> 
inline
void
ODESLV_BASE::_record
( std::ofstream& ofile, VRES const& res, unsigned const iprec )
{
  if( !ofile ) return;

  // Specify format
  ofile << std::right << std::scientific << std::setprecision(iprec);

  // Record computed states at stage times
  auto it = res.begin(), it0 = it;
  for( ; it != res.end(); ++it ){
    if( it != res.begin() && it->t == it0->t )
      ofile << std::endl;
    ofile << std::setw(iprec+9) << it->t;
    for( auto const& xi : it->x )
      ofile << std::setw(iprec+9) << xi;
    ofile << std::endl;
    it0 = it;
  }
}

template <typename REALTYPE>
inline
void
ODESLV_BASE::_vec2D
( REALTYPE const* vec, unsigned const n, double* d )
{
  for( unsigned i=0; i<n; i++  ) d[i] = vec[i];
  return;
}

template <typename REALTYPE>
inline
void
ODESLV_BASE::_D2vec
( double const* d, unsigned const n, REALTYPE* vec )
{
  for( unsigned i=0; i<n; i++  ) vec[i] = d[i];
  return;
}

inline
bool
ODESLV_BASE::_INI_D_STA
( double const* p, double const* c )
{
  // Set constants
  for( unsigned ic=0; c && ic<_nc; ++ic ) _pC[ic].set( c[ic] );

  // Size and set DAG evaluation arrays
  _nVAR0 = _nx+_np+1;
  _nVAR  = _nVAR0+_nq;
  delete[] _pVAR; _pVAR = new FFVar[_nVAR];
  delete[] _DVAR; _DVAR = new double[_nVAR];
  for( unsigned ix=0; ix<_nx; ix++ ) _pVAR[ix] = _pX[ix];
  for( unsigned ip=0; ip<_np; ip++ ) _pVAR[_nx+ip] = _pP[ip];
  _pVAR[_nx+_np] = (_pT? *_pT: 0. );
  for( unsigned iq=0; iq<_nq; iq++ ) _pVAR[_nx+_np+1+iq] = _pQ?_pQ[iq]:0.;

  _Dx = _DVAR;
  _Dp = _Dx + _nx;
  _Dt = _Dp + _np;
  _Dq = _Dt + 1;
  for( unsigned ip=0; ip<_np; ip++ ) _Dp[ip] = p[ip];

  _Df.resize( _nf );

  return true;
}

inline
bool
ODESLV_BASE::_END_D_STA
()
{
  // Set constants
  for( unsigned ic=0; ic<_nc; ++ic ) _pC[ic].unset();

  return true;
}

template <typename REALTYPE>
inline
void
ODESLV_BASE::_GET_D_STA
( REALTYPE const* x, REALTYPE const* q )
{
  _vec2D( x, _nx, _Dx );
  if( q ) _vec2D( q, _nq, _Dq );
}

inline
bool
ODESLV_BASE::_IC_D_SET
()
{
  if( !_vIC.size() || _nx0 != _nx ) return false;
  _pIC = _vIC.at(0);
  return true;
}

template <typename REALTYPE>
inline
bool
ODESLV_BASE::_IC_D_STA
( double const& t, REALTYPE* x )
{
  *_Dt = t; // current time
  _dag->eval( _nx, _pIC, (double*)x, _np+1, _pVAR+_nx, _DVAR+_nx );
  return true;
}

template <typename REALTYPE>
inline
bool
ODESLV_BASE::_IC_D_QUAD
( REALTYPE* q )
{
  for( unsigned iq=0; iq<_nq; iq++ ) q[iq] = 0.;
  return true;
}

inline
bool
ODESLV_BASE::_CC_D_SET
( unsigned const iIC )
{
  if( _vIC.size() <= iIC || _nx0 != _nx ) return false;
  _pIC = _vIC.at( iIC );
  return true;
}

template <typename REALTYPE>
inline
bool
ODESLV_BASE::_CC_D_STA
( double const& t, REALTYPE* x )
{
  *_Dt = t; // current time
  _vec2D( x, _nx, _Dx ); // current state
  _dag->eval( _nx, _pIC, (double*)x, _nVAR0, _pVAR, _DVAR );
  return true;
}

inline
bool
ODESLV_BASE::_RHS_D_SET
( unsigned const iRHS, unsigned const iQUAD )
{
  if( _vRHS.size() <= iRHS ) return false;
  _pRHS = _vRHS.at( iRHS );

  if( _nq && _vQUAD.size() <= iQUAD ) return false;
  _pQUAD = _nq? _vQUAD.at( iQUAD ): 0;

  // Generate Jacobian using sparse forward AD
  delete[] std::get<1>(_pJAC); delete[] std::get<2>(_pJAC); delete[] std::get<3>(_pJAC);
  _pJAC = _dag->SFAD( _nx, _pRHS, _nx, _pX ); // Jacobian in sparse format, ordered columnwise
  _pJACCOLNDX.resize( _nx+1 );
  for( unsigned ie=0, ic=0; ie<std::get<0>(_pJAC); ++ie ){
#ifdef CRONOS__ODESLV_BASE_DEBUG
    std::cout << "  JAC[" << std::get<1>(_pJAC)[ie] << ", " << std::get<2>(_pJAC)[ie] << "]" << std::endl;
#endif
    for( ; std::get<2>(_pJAC)[ie] >= ic; ++ic ){
      _pJACCOLNDX[ic] = ie;
#ifdef CRONOS__ODESLV_BASE_DEBUG
      std::cout << "  JACCOLNDX[" << ic << "] = " << ie << std::endl;
#endif
    }
  }
  _pJACCOLNDX[_nx] = std::get<0>(_pJAC);
#ifdef CRONOS__ODESLV_BASE_DEBUG
  std::cout << "  JACCOLNDX[" << _nx << "] = " << std::get<0>(_pJAC) << std::endl;
  std::cout << "PAUSED - <1> TO CONTINUE"; int dum; std::cin >> dum;
#endif

  return _RHS_D_SET();
}

inline
bool
ODESLV_BASE::_RHS_D_SET
()
{
  _opRHS  = _dag->subgraph( _nx, _pRHS );

  if( _pQUAD ) _opQUAD = _dag->subgraph( _nq, _pQUAD );

  _opJAC = _dag->subgraph( std::get<0>(_pJAC), std::get<3>(_pJAC) );
  _DJAC.resize( std::get<0>(_pJAC) );

  return true;
}

template <typename REALTYPE>
inline
bool
ODESLV_BASE::_RHS_D_STA
( double const& t, REALTYPE const* x, REALTYPE* xdot )
{
  if( !_pRHS ) return false;
  *_Dt = t; // current time
  _vec2D( x, _nx, _Dx ); // current state
  _dag->eval( _opRHS, _DWRK, _nx, _pRHS, (double*)xdot, _nVAR0, _pVAR, _DVAR );
  return true;
}

template <typename REALTYPE>
inline
bool
ODESLV_BASE::_RHS_D_QUAD
( double const& t, REALTYPE const* x, REALTYPE* qdot )
{
  if( !_pQUAD ) return false;
  // No need to update _DVAR
  _dag->eval( _opQUAD, _DWRK, _nq, _pQUAD, (double*)qdot, _nVAR0, _pVAR, _DVAR );
  return true;
}

template <typename REALTYPE>
inline
bool
ODESLV_BASE::_JAC_D_STA
( double const& t, REALTYPE const* x, REALTYPE** jac )
{
  // No need to update _DVAR
  _dag->eval( _opJAC, _DWRK, std::get<0>(_pJAC), std::get<3>(_pJAC),
              _DJAC.data(), _nVAR0, _pVAR, _DVAR );
  for( unsigned ie=0; ie<std::get<0>(_pJAC); ++ie ){
    jac[std::get<2>(_pJAC)[ie]][std::get<1>(_pJAC)[ie]] = _DJAC[ie];
#ifdef CRONOS__ODESLV_BASE_DEBUG
    std::cout << "  jac[" << std::get<1>(_pJAC)[ie] << ", "
              << std::get<2>(_pJAC)[ie] << "] = " << _DJAC[ie] << std::endl;
#endif
  }
  return true;
}

#if defined( CRONOS__WITH_KLU )
template <typename REALTYPE, typename INDEXTYPE>
inline
bool
ODESLV_BASE::_JAC_D_STA
( double const& t, REALTYPE const* x, REALTYPE* jac, INDEXTYPE* ptrs, INDEXTYPE* vals )
{
  // No need to update _DVAR
  _dag->eval( _opJAC, _DWRK, std::get<0>(_pJAC), std::get<3>(_pJAC),
              (double*)jac, _nVAR0, _pVAR, _DVAR );
  for( unsigned ie=0; ie<std::get<0>(_pJAC); ++ie ){
    vals[ie] = (INDEXTYPE)std::get<1>(_pJAC)[ie];
#ifdef CRONOS__ODESLV_BASE_DEBUG
    std::cout << "  jac[" << ie << "] = " << jac[ie] << std::endl;
    std::cout << "  vals[" << ie << "] = " << vals[ie] << std::endl;
#endif
  }
  for( unsigned ic=0; ic<=_nx; ++ic ){
    ptrs[ic] = (INDEXTYPE) _pJACCOLNDX[ic];
#ifdef CRONOS__ODESLV_BASE_DEBUG
    std::cout << "  ptrs[" << ic << "] = " << ptrs[ic] << std::endl;
#endif
  }
  return true;
}
#endif

inline
bool
ODESLV_BASE::_FCT_D_STA
( unsigned const iFCT, double const& t )
{
  if( !_nf ) return true;
  *_Dt = t; // current time
  FFVar const* pFCT = _vFCT.at( iFCT );
  //std::cout << "evaluating functions @" << t << std::endl;
  _dag->eval( _nf, pFCT, _Df.data(), _nVAR, _pVAR, _DVAR, iFCT?true:false );
  return true;
}

} // end namescape mc

#endif

