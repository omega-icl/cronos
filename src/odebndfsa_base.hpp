// Copyright (C) 2015 Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

#ifndef MC__ODEBNDS_BASE_HPP
#define MC__ODEBNDS_BASE_HPP

#include "odebnd_base.hpp"

namespace mc
{
//! @brief C++ base class for computing enclosures of parametric ODEs forward/adjoint sensitivity using continuous-time set-valued integration.
////////////////////////////////////////////////////////////////////////
//! mc::ODEBNDS_BASE is a C++ base class for computing enclosures of the
//! reachable set of parametric ODEs forward/adjoint sensitivity
//! using continuous-time set-valued integration.
////////////////////////////////////////////////////////////////////////
template <typename T, typename PMT, typename PVT>
class ODEBNDS_BASE: public virtual BASE_DE, public virtual ODEBND_BASE
{
 typedef Ellipsoid E;

 public:
  /** @defgroup ODEBNDS_BASE Continuous-time set-valued integration of parametric ODEs with sensitivity analysis
   *  @{
   */
  //! @brief Default constructor
  ODEBNDS_BASE();

  //! @brief Virtual destructor
  virtual ~ODEBNDS_BASE();

 protected:
  //! @brief list of operations in RHS sensitivity
  std::list<const FFOp*> _opFSARHS;

  //! @brief array of list of operations in individual RHS sensitivity
  std::list<const FFOp*> *_opFSARHSi;

  //! @brief list of operations in quadrature Jacobian
  std::list<const FFOp*> _opDQUAD;

  //! @brief const pointer to RHS sensitivity function in current stage of ODE system
  const FFVar* _pFSARHS;

  //! @brief const pointer to quadrature Jacobian in current stage of ODE system
  const FFVar* _pDQUAD;

  //! @brief state sensitivity interval enclosures
  T *_Ixp;

  //! @brief state sensitivity derivative interval enclosures
  T *_Ixpdot;

  //! @brief quadrature sensitivity interval enclosures
  T *_PMqp;

  //! @brief quadrature sensitivity derivative interval enclosures
  T *_Iqpdot;

  //! @brief state sensitivity polynomial models
  PVT *_PMxp;

  //! @brief state sensitivity remainder bounds
  T *_Rxp;

  //! @brief state sensitivity derivative polynomial models
  PVT *_PMxpdot;

  //! @brief state sensitivity derivative remainder bounds
  T *_Rxpdot;

  //! @brief quadrature sensitivity polynomial models
  PVT *_PMqp;

  //! @brief quadrature sensitivity remainder bounds
  T *_Rqp;

  //! @brief quadrature sensitivity derivative polynomial models
  PVT *_PMqpdot;

  //! @brief quadrature sensitivity derivative remainder bounds
  T *_Rqpdot;

  //! @brief Function to initialize state/sensitivity interval bounding
  template <typename OPT> bool _INI_I_FSA
    ( const OPT&options, const unsigned np, const T*Ip, const unsigned ns );

  //! @brief Function to initialize quarature/sensitivity interval bounds
  template <typename REALTYPE, typename OPT> bool _IC_I_FSA_QUAD
    ( const OPT&options, REALTYPE*vec );

  //! @brief Function to initialize state/sensitivity interval bounds
  template <typename REALTYPE, typename OPT> bool _IC_I_FSA_STA
    ( const OPT&options, REALTYPE*vec );

  //! @brief Function to reinitialize state/sensitivity interval bounds
  template <typename REALTYPE, typename OPT> bool _CC_I_FSA_STA
    ( const OPT&options, const unsigned iIC, const double t, REALTYPE*vec );

  //! @brief Function to set RHS and QUAD pointers
  template <typename OPT> bool _SET_I_FSA_STA
    ( const OPT&options, const unsigned iRHS, const unsigned iQUAD );

  //! @brief Function to calculate the RHS of auxiliary ODEs in interval arithmetic
  template <typename REALTYPE, typename OPT> bool _RHS_I_FSA_STA
    ( const OPT&options, double t, const REALTYPE*x, REALTYPE*xdot );

  //! @brief Function to calculate the RHS of auxiliary ODEs in interval arithmetic
  template <typename REALTYPE, typename OPT> bool _RHS_I_FSA_QUAD
    ( const OPT&options, double t, const REALTYPE*x, REALTYPE*qdot,
      const bool bndinit=true );

  //! @brief Function to calculate the Jacobian of auxiliary ODEs in interval arithmetic
  template <typename REALTYPE, typename OPT> bool _JAC_I_FSA_STA
    ( const OPT&options, double t, const REALTYPE*x, REALTYPE*jac, REALTYPE*xdot );

  //! @brief Function to calculate the functional derivatives at intermediate/end point
  bool _FCT_I_FSA_STA
    ( const unsigned iFCT, const double t, T*If );

  //! @brief Function to initialize GSL for state sensitivity polynomial models
  template <typename OPT> bool _INI_PM_FSA_STA
    ( const OPT&options, const unsigned np, const PVT*PMp, const unsigned ns );

  //! @brief Function to initialize quarature sensitivity polynomial models
  template <typename REALTYPE, typename OPT> bool _IC_PM_FSA_QUAD
    ( const OPT&options, REALTYPE*vec );

  //! @brief Function to initialize state sensitivity polynomial models
  template <typename REALTYPE, typename OPT> bool _IC_PM_FSA_STA
    ( const OPT&options, REALTYPE*vec );

  //! @brief Function to reinitialize state sensitivity polynomial bounds
  template <typename REALTYPE, typename OPT> bool _CC_PM_STA
    ( const OPT&options, const unsigned iIC, const double t, REALTYPE*vec );

  //! @brief Function to set RHS and QUAD pointers
  template <typename OPT> bool _SET_PM_FSA_STA
    ( const OPT&options, const unsigned iRHS, const unsigned iQUAD );

  //! @brief Function to calculate the RHS of auxiliary ODEs in polynomial mode arithmetic
  template <typename REALTYPE, typename OPT> bool _RHS_PM_FSA_STA
    ( const OPT&options, double t, const REALTYPE*x, REALTYPE*xdot );

  //! @brief Function to calculate the RHS of auxiliary ODEs in polynomial mode arithmetic
  template <typename REALTYPE, typename OPT> bool _RHS_PM_FSA_QUAD
    ( const OPT&options, double t, const REALTYPE*x, REALTYPE*qdot,
      const bool bndinit=true );

  //! @brief Function to calculate the Jacobian of auxiliary ODEs in polynomial mode arithmetic
  template <typename REALTYPE, typename OPT> bool _JAC_PM_FSA_STA
    ( const OPT&options, double t, const REALTYPE*x, REALTYPE*jac, REALTYPE*xdot );
iprec=5 );

  //! @brief Private methods to block default compiler methods
  ODEBNDS_BASE(const ODEBNDS_BASE&);
  ODEBNDS_BASE& operator=(const ODEBNDS_BASE&);
};

template <typename T, typename PMT, typename PVT>
inline
ODEBNDS_BASE<T,PMT,PVT>::ODEBNDS_BASE
()
: ODEBND_BASE(), _pFSARHS(0), _pDQUAD(0)
{
  // Initalize state/parameter arrays
  _opFSARHSi = 0;
  _Ixp = _Ixpdot = _Iqp = _Iqpdot = _Rxp = _Rxpdot = _Rqp = _Rqpdot = 0;
  _PMxp = _PMxpdot = _PMqp = _PMqpdot = 0;
}

template <typename T, typename PMT, typename PVT>
inline
ODEBND_BASE<T,PMT,PVT>::~ODEBND_BASE
()
{
  delete[] _opRHSi;
  delete[] _pFSARHS;
  delete[] _pDQUAD;

  // Free state/quadrature arrays -- Do *NOT* delete _Ip _PMp _PMenv
  delete[] _Ixp;
  delete[] _Ixpdot;
  delete[] _Iqp;
  delete[] _Iqpdot;
  delete[] _Rxp;
  delete[] _Rxpdot;
  delete[] _Rqp;
  delete[] _Rqpdot;
  delete[] _PMxp;
  delete[] _PMxpdot;
  delete[] _PMqp;
  delete[] _PMqpdot;
}

template <typename T, typename PMT, typename PVT>
template <typename REALTYPE, typename OPT> inline bool
ODEBND_BASE<T,PMT,PVT>::_IC_I_FSA_QUAD
( const OPT&options, REALTYPE*vec )
{
  if( !_vQUAD.size() || !_nq ) return true;
  for( unsigned iq=0; iq<_nq; iq++ ) _Iq[iq] = 0.;
  _I2vec( _nq, _Iq, vec );
  return true;
}

template <typename T, typename PMT, typename PVT>
template <typename REALTYPE, typename OPT> inline bool
ODEBND_BASE<T,PMT,PVT>::_IC_I_FSA_STA
( const OPT&options, REALTYPE*vec )
{
  if( !_vIC.size() || _nx0 != _nx ) return false;
  _pIC = _vIC.at(0);
  _opIC = _pDAG->subgraph( _nx, _pIC );

  switch( options.WRAPMIT){

  case OPT::NONE:
  case OPT::DINEQ:
    _pDAG->eval( _opIC, _nx, _pIC, _Ix, _npar, _pVAR+_nx, _Ip );
    _I2vec( _nx, _Ix, vec );
    break;

  case OPT::ELLIPS:
  default:
    for( unsigned ip=0; ip<_npar; ip++ ){
      _pref[ip] = Op<T>::mid( _Ip[ip] );
      _MVPp[ip].set( _MVPenv, ip, _Ip[ip] );
    }
    _IC_I_ELL( _pDAG, _opIC, _nx, _pIC, _npar, _pVAR+_nx, _MVPp, _MVPx,
               _xref, _Q, _npar, _B, _Ix );
    _E2vec( _nx, _npar, _xref, _Q, _B, vec );
    break;
  }

  return true;
}

template <typename T, typename PMT, typename PVT>
template <typename REALTYPE, typename OPT> inline bool
ODEBND_BASE<T,PMT,PVT>::_CC_I_FSA_STA
( const OPT&options, const unsigned iIC, const double t, REALTYPE*vec )
{
  if( _vIC.size() <= iIC || _nx0 != _nx ) return false;
  _pIC = _vIC.at( iIC );
  _opIC = _pDAG->subgraph( _nx, _pIC );
  delete[] _IIC;   _IIC = 0;
  delete[] _PMIC;  _PMIC = 0;

  switch( options.WRAPMIT){
  case OPT::NONE:
  case OPT::DINEQ:{
    _vec2I( vec, _nx, _Ix ); // current state bounds
    _IVAR[_nx+_npar] = t; // current time
    _IIC = new T[_opIC.size()];
    _pDAG->eval( _opIC, _IIC, _nx, _pIC, _Ixdot, _nVAR-_nq, _pVAR, _IVAR );
    _I2vec( _nx, _Ixdot, vec );
    break;
   }
  case OPT::ELLIPS:
  default:{
    _vec2E( vec, _nx, _npar, _Q, _Er, _Ir, _pref, _Ip, _B, _xref, _Ix ); // current state enclosure
    for( unsigned jx=0; jx<_nx; jx++ )
      _MVXPd[jx].set( _MVXPenv, jx, _Ir[jx] );
    for( unsigned jp=0; jp<_npar; jp++ )
      _MVXPp[jp].set( _MVXPenv, _nx+jp, _Ip[jp] );
    _ep2x( _nx, _npar, _MVXPd, _pref, _MVXPp, _B, _xref, _MVXPx );
    _MVXPVAR[_nx+_npar] = t; // current time
    _PMIC = new PVT[_opIC.size()];
    _CC_I_ELL( _pDAG, _opIC, _PMIC, _nx, _pIC, _nVAR-_nq, _pVAR, _MVXPVAR,
       _MVXPf, _Er, _A, _npar, _xrefdot, _Bdot, _Irdot, _Qdot, options.QTOL,
       machprec() );
    _E2vec( _nx, _npar, _xrefdot, _Qdot, _Bdot, vec );
    break;
   }
  }

  return true;
}

template <typename T, typename PMT, typename PVT>
template <typename REALTYPE, typename OPT> inline bool
ODEBND_BASE<T,PMT,PVT>::_RHS_I_FSA_STA
( const OPT&options, double t, const REALTYPE*x, REALTYPE*xdot )
{
  if( !_pRHS ) return false;

  switch( options.WRAPMIT){
  case OPT::NONE:
    _vec2I( x, _nx, _Ix );   // set current state bounds
    _IVAR[_nx+_npar] = t; // set current time
    _RHS_I_NONE( _pDAG, _opRHS, _IRHS, _nx, _pRHS, _nVAR-_nq, _pVAR, _IVAR, _Ixdot );
    _I2vec( _nx, _Ixdot, xdot );
    return true;
   
  case OPT::DINEQ:
    if( !_opRHSi ) return GSL_EBADFUNC;
    _vec2I( x, _nx, _Ix );   // set current state bounds
    _IVAR[_nx+_npar] = t; // set current time
    _RHS_I_DI( _pDAG, _opRHSi, _IRHS, _nx, _pRHS, _nVAR-_nq, _pVAR, _IVAR, _Ixdot,
               _xLdot, _xUdot );
    _I2vec( _nx, _xLdot, _xUdot, xdot );
    return true;
   
  case OPT::ELLIPS:
  default:{
    _vec2E( x, _nx, _npar, _Q, _Er, _Ir, _pref, _Ip, _B, _xref, _Ix ); // set current state enclosure
#ifdef MC__ODEBND_BASE_DINEQI_DEBUG
    std::cout << "@t=" << t << std::endl;
    //E::options.PSDCHK = true;
    for( unsigned ix=0; ix<_nx; ix++ )
      std::cout << "xref[" << ix << "] = " << _xref[ix] << std::endl;
    for( unsigned ix=0; ix<_nx; ix++ ){
      std::cout << "B[" << ix << ",#] = ";
      for( unsigned jp=0; jp<_npar; jp++ )
        std::cout << _B[jp*_nx+ix] << "  ";
      std::cout << std::endl;
    }
    std::cout << "Er = " << _Er << std::endl;
    //{ int dum; std::cin >> dum; }
#endif

    // Compute polynomial expansion of ODE's RHS
    for( unsigned jx=0; jx<_nx; jx++ )
      _MVXPd[jx].set( _MVXPenv, jx, _Ir[jx] );
    for( unsigned jp=0; jp<_npar; jp++ )
      _MVXPp[jp].set( _MVXPenv, _nx+jp, _Ip[jp] );
    _ep2x( _nx, _npar, _MVXPd, _pref, _MVXPp, _B, _xref, _MVXPx );
#ifdef MC__ODEBND_BASE_DINEQI_DEBUG
    for( unsigned ix=0; ix<_nx; ix++ )
      std::cout << "MVXPx[ " << ix << "] = " << _MVXPx[ix] << std::endl;
    { int dum; std::cin >> dum; }
#endif
    _MVXPVAR[_nx+_npar] = t; // set current time
    _RHS_I_ELL( _pDAG, _opRHS, _PMRHS, _nx, _pRHS, _nVAR-_nq, _pVAR, _MVXPVAR,
       _MVXPf, _Q, _A, _npar, _xrefdot, _Bdot, _Irdot, _Qdot, options.QTOL,
       machprec() );
    _E2vec( _nx, _npar, _xrefdot, _Qdot, _Bdot, xdot );
#ifdef MC__ODEBND_BASE_DINEQI_DEBUG
    //E::options.PSDCHK = true;
    for( unsigned ix=0; ix<_nx; ix++ )
      std::cout << "xrefdot[" << ix << "] = " << _xrefdot[ix] << std::endl;
    for( unsigned ix=0; ix<_nx; ix++ ){
      std::cout << "Bdot[" << ix << ",#] = ";
      for( unsigned jp=0; jp<_npar; jp++ )
        std::cout << _Bdot[jp*_nx+ix] << "  ";
      std::cout << std::endl;
    }
    std::cout << "Er = " << E(_nx,_Qdot) << std::endl;
    { int dum; std::cin >> dum; }
#endif
    
    return true;
   }
  }
}

template <typename T, typename PMT, typename PVT>
template <typename REALTYPE, typename OPT> inline bool
ODEBND_BASE<T,PMT,PVT>::_RHS_I_FSA_QUAD
( const OPT&options, double t, const REALTYPE*x, REALTYPE*qdot,
  const bool bndinit )
{
  if( !_pQUAD ) return false;

  switch( options.WRAPMIT){
  case OPT::NONE:
  case OPT::DINEQ:
    if( !bndinit ) break;
    _vec2I( x, _nx, _Ix );   // set current state bounds
    break;

  case OPT::ELLIPS:
  default:
    if( !bndinit ) break;
    _vec2E( x, _nx, _npar, _Q, _Er, _Ir, _pref, _Ip, _B, _xref, _Ix ); // set current state enclosure
    break;
  }

  _IVAR[_nx+_npar] = t; // set current time
  _QUAD_I( _pDAG, _opQUAD, _IRHS, _nq, _pQUAD, _nVAR-_nq, _pVAR, _IVAR, _Iqdot );
  _I2vec( _nq, _Iqdot, qdot );
    
   return true;
}

template <typename T, typename PMT, typename PVT>
template <typename REALTYPE, typename OPT> inline bool
ODEBND_BASE<T,PMT,PVT>::_JAC_I_FSA_STA
( const OPT&options, double t, const REALTYPE*x, REALTYPE*jac, REALTYPE*xdot )
{
  // Jacobian not (yet) implemented
  return false;
}

template <typename T, typename PMT, typename PVT>
inline bool
ODEBND_BASE<T,PMT,PVT>::_FCT_I_FSA_STA
( const unsigned iFCT, const double t, T*If )
{
  if( !_nf || !If ) return true;

  _IVAR[_nx+_npar] = t; // set current time
  const FFVar* pFCT = _vFCT.at( iFCT );
  _pDAG->eval( _nf, pFCT, If, _nVAR, _pVAR, _IVAR, iFCT?true:false );
  return true;
}

template <typename T, typename PMT, typename PVT>
template <typename OPT> inline bool
ODEBND_BASE<T,PMT,PVT>::_INI_I_FSA_STA
( const OPT&options, const unsigned np, const T*Ip, const unsigned ns )
{
  // Update effective number of parameters
  // (possibly larger than _np if lifting is used)
  _npar = np;

  // Size and set DAG evaluation arrays
  _nVAR = _nx+_npar+1+_nq;
  delete[] _pVAR; _pVAR = new FFVar[_nVAR];
  delete[] _IVAR; _IVAR = new T[_nVAR];
  for( unsigned ix=0; ix<_nx; ix++ ) _pVAR[ix] = _pX[ix];
  for( unsigned ip=0; ip<_npar; ip++ ) _pVAR[_nx+ip] = _pP[ip];
  _pVAR[_nx+_npar] = (_pT? *_pT: 0. );
  for( unsigned iq=0; iq<_nq; iq++ ) _pVAR[_nx+_npar+1+iq] = _pQ?_pQ[iq]:0.;
  for( unsigned ip=0; ip<_npar; ip++ ) _IVAR[_nx+ip] = Ip[ip];
  _Ix = _IVAR;
  _Ip = _IVAR + _nx;
  _Iq = _IVAR + _nx+_npar+1;

  // Reset _MVXPVenv and related variables
  if( _MVXPenv && ( _MVXPenv->nord() != options.ORDMIT 
                 || _MVXPenv->nvar() != _nx+_npar ) ){
    delete[] _MVXPf;   _MVXPf = 0;
    delete[] _MVXPd;   _MVXPd = 0;
    delete[] _MVXPVAR; _MVXPVAR = 0;
    delete   _MVXPenv; _MVXPenv = 0;
  }

  // Reset _MVPVenv and related variables
  if( _MVPenv && ( _MVPenv->nord() != 1 
                || _MVPenv->nvar() != _npar ) ){
    delete[] _MVXPx;   _MVPx = 0;
    delete[] _MVXPp;   _MVPp = 0;
    delete   _MVPenv; _MVPenv = 0;
  }

  // Set parameterization variables
  delete[] _Iqdot; _Iqdot = new T[_nq];

  switch( options.WRAPMIT){
  case OPT::DINEQ:
    delete[] _xLdot;   _xLdot = new double[_nx];
    delete[] _xUdot;   _xUdot = new double[_nx];
  case OPT::NONE:
    delete[] _Ixdot;   _Ixdot = new T[_nx];
    break;

  case OPT::ELLIPS:
  default:
    delete[] _pref;    _pref    = new double[_npar];
    delete[] _xref;    _xref    = new double[_nx];
    delete[] _A;       _A       = new double[_nx*_nx];
    delete[] _B;       _B       = new double[_nx*_npar];
    delete[] _Q;       _Q       = new double[_nx*(_nx+1)/2];
    delete[] _Ir;      _Ir      = new T[_nx];
    delete[] _xrefdot; _xrefdot = new double[_nx];
    delete[] _Bdot;    _Bdot    = new double[_nx*_npar];
    delete[] _Qdot;    _Qdot    = new double[_nx*(_nx+1)/2];
    delete[] _Irdot;   _Irdot   = new T[_nx];
    delete   _MVXPenv; _MVXPenv = new PMT( _nx+_npar, options.ORDMIT );
    _MVXPenv->options = options.PMOPT;
    delete[] _MVXPd;   _MVXPd   = new PVT[_nx];
    delete[] _MVXPf;   _MVXPf   = new PVT[_nx];
    delete[] _MVXPVAR; _MVXPVAR = new PVT[_nVAR-_nq];
    _MVXPp = _MVXPVAR + _nx;
    _MVXPx = _MVXPVAR;
    delete   _MVPenv;  _MVPenv  = new PMT( _npar, 1 );
    _MVPenv->options = options.PMOPT;
    delete[] _MVPp;    _MVPp    = new PVT[_npar];
    delete[] _MVPx;    _MVPx    = new PVT[_nx];
    break;
  }

  return true;
}

template <typename T, typename PMT, typename PVT>
template <typename OPT> inline bool
ODEBND_BASE<T,PMT,PVT>::_SET_I_FSA_STA
( const OPT&options, const unsigned iRHS, const unsigned iQUAD )
{
  if( _vRHS.size() <= iRHS ) return false;
  if( _nq && _vQUAD.size() <= iQUAD ) return false;

  _pRHS   = _vRHS.at( iRHS );
  _opRHS  = _pDAG->subgraph( _nx, _pRHS );
  _pQUAD  = _nq? _vQUAD.at( iQUAD ): 0;
  _opQUAD = _pDAG->subgraph( _nq, _pQUAD );
  const unsigned opmax = _opRHS.size()>_opQUAD.size()?_opRHS.size():_opQUAD.size();

  delete[] _IRHS;   _IRHS = 0;
  delete[] _PMRHS;  _PMRHS = 0;
  delete[] _opRHSi; _opRHSi = 0;
  switch( options.WRAPMIT){
  case OPT::NONE:
    _IRHS = new T[ opmax ];
    break;
  case OPT::DINEQ:
    _opRHSi = new std::list<const FFOp*>[_nx];
    for( unsigned ix=0; ix<_nx; ix++ )
      _opRHSi[ix] = _pDAG->subgraph( 1, _pRHS+ix );
    _IRHS = new T[ opmax ];
    break;   
  case OPT::ELLIPS:
  default:
    _IRHS = new T[ _opQUAD.size() ];
    _PMRHS = new PVT[ _opRHS.size() ];
    break;
  }

  delete[] _pJAC;
  _pJAC = 0;
  _opJAC.clear();

  return true;
}

template <typename T, typename PMT, typename PVT>
template <typename REALTYPE, typename OPT> inline bool
ODEBND_BASE<T,PMT,PVT>::_IC_PM_FSA_STA
( const OPT&options, REALTYPE*vec )
{
  if( !_vIC.size() || _nx0 != _nx ) return false;
  _pIC = _vIC.at(0);
  _opIC = _pDAG->subgraph( _nx, _pIC );

  switch( options.WRAPMIT){

  case OPT::NONE:
  case OPT::DINEQ:
    _pDAG->eval( _opIC, _nx, _pIC, _PMx, _npar, _pVAR+_nx, _PMp );
    // Whether or not to ignore the remainder
    if( !options.PMNOREM )
      _PMI2vec( _PMenv, _nx, _PMx, vec, false );
    else
      _PMI2vec( _PMenv, _nx, _PMx, 0, vec );
    break;

  case OPT::ELLIPS:
  default:
    _IC_PM_ELL( _pDAG, _opIC, _nx, _pIC, _npar, _pVAR+_nx, _PMp, _PMx, _Q, _Er, _Ir );
    // Whether or not to ignore the remainder
    if( !options.PMNOREM )
      _PME2vec( _PMenv, _nx, _PMx, _Q, vec );
    else
      _PME2vec( _PMenv, _nx, _PMx, 0, vec );
    break;
  }
  return true;
}

template <typename T, typename PMT, typename PVT>
template <typename REALTYPE, typename OPT> inline bool
ODEBND_BASE<T,PMT,PVT>::_CC_PM_FSA_STA
( const OPT&options, const unsigned iIC, const double t, REALTYPE*vec )
{
  if( _vIC.size() <= iIC || _nx0 != _nx ) return false;
  _pIC = _vIC.at( iIC );
  _opIC = _pDAG->subgraph( _nx, _pIC );
  delete[] _IIC;   _IIC = 0;
  delete[] _PMIC;  _PMIC = 0;
  delete[] _IDIC;  _IDIC = 0;
  delete[] _PMDIC; _PMDIC = 0;
  delete[] _pDIC;  _pDIC = 0; _opDIC.clear();

  switch( options.WRAPMIT){
  case OPT::NONE:
  case OPT::DINEQ:
    _vec2PMI( vec, _PMenv, _nx, _PMx ); // current state/quadrature polynomial model
    _PMVAR[_nx+_npar] = t; // current time
    _PMIC = new PVT[_opIC.size()];
    _pDAG->eval( _opIC, _PMIC, _nx, _pIC, _PMxdot, _nVAR-_nq, _pVAR, _PMVAR );
    // Whether or not to ignore the remainder
    if( !options.PMNOREM )
      _PMI2vec( _PMenv, _nx, _PMxdot, vec, false );
    else
      _PMI2vec( _PMenv, _nx, _PMxdot, 0, vec );
    break;
   
  case OPT::ELLIPS:
  default:{
    _vec2PME( vec, _PMenv, _nx, _PMx, _Q, _Er, _Ir ); // current state/quadrature polynomial model
    _PMVAR[_nx+_npar] = t; // current time   
    _PMIC = new PVT[_opIC.size()];

    // In this variant a bound on the Jacobian matrix is computed and the
    // linear part is taken as the mid-point of this matrix
    if( !options.ORDMIT ){
      for( unsigned ix=0; ix<_nx; ix++ ){
        _IVAR[ix] = _PMx[ix].bound(); // set current state bounds
        _PMVAR[ix].center().set( T(0.) ); // cancel remainder term
      }
      _IVAR[_nx+_npar] = t; // current time
      _pDIC = _pDAG->FAD( _nx, _pIC, _nx, _pVAR );
      _opDIC = _pDAG->subgraph( _nx*_nx, _pDIC );
      _IDIC = new T[ _opDIC.size() ];
      _RHS_PM_ELL0( _pDAG, _opIC, _PMIC, _opDIC, _IDIC, _nx, _pIC,
                    _pDIC, _nVAR-_nq, _pVAR, _PMVAR, _IVAR, _PMxdot, _Idfdx,
                    _Ir, _A, _Irdot );
    }

    // In this variant a polynomial model of the Jacobian matrix is computed and the
    // linear part is taken as the mid-point of this matrix
    else if( _PMenv->nord() >= _MVXPenv->nord() ){
#ifdef MC__ODEBND_BASE_MVXP_USE
      for( unsigned jx=0; jx<_nx; jx++ ){
        _MVXPd[jx].set( _MVXPenv, _npar+jx, _Ir[jx] );
        _MVXPx[jx].set( _MVXPenv ).set( _PMx[jx].center().set( T(0.) ), true );
      }
      _e2x( _nx, _MVXPd, _MVXPx, false );
#else
      for( unsigned jx=0; jx<_nx; jx++ ){
        _MVXPx[jx].set( _MVXPenv ).set( _PMx[jx].center(), true );
        _PMx[jx].set( T(0.) );
      }
#endif
      _MVXPVAR[_nx+_npar] = t; // current time
      _pDIC = _pDAG->FAD( _nx, _pIC, _nx, _pVAR );
      _opDIC = _pDAG->subgraph( _nx*_nx, _pDIC );
      _PMDIC = new PVT[ _opDIC.size() ];
      _RHS_PM_ELL1( _pDAG, _opIC, _PMIC, _opDIC, _PMDIC, _nx, _pIC,
                    _pDIC, _nVAR-_nq, _pVAR, _PMVAR, _MVXPVAR, _PMxdot,
                    _MVXPdfdx, _Ir, _A, _Irdot );
    }

    // In this variant a polynomial model in the joint state-parameter and
    // of the same order as the parameter polynomial model is computed
    else{
      for( unsigned jx=0; jx<_nx; jx++ ){
        _MVXPd[jx].set( _MVXPenv, _npar+jx, _Ir[jx] );
        _MVXPx[jx].set( _MVXPenv ).set( _PMx[jx].center().set( T(0.) ), true );
      }
      _e2x( _nx, _MVXPd, _MVXPx, false );
      _MVXPVAR[_nx+_npar] = t; // current time   
      _RHS_PM_ELL2( _pDAG, _opIC, _PMIC, _nx, _pRHS, _nVAR-_nq, _pVAR,
                    _MVXPVAR, _PMenv, _PMxdot, _MVXPf, _npar, _Ir, _A, _Irdot );
    }

    // Whether or not to ignore the remainder
    if( !options.PMNOREM ){
      _CC_PM_ELL( _nx, _Er, _A, _Irdot, _Qdot, options.QTOL, machprec() );
      _PME2vec( _PMenv, _nx, _PMxdot, _Qdot, vec );
    }
    else
      _PME2vec( _PMenv, _nx, _PMxdot, 0, vec );
    break;
   }
  }

  return true;
}

template <typename T, typename PMT, typename PVT>
template <typename REALTYPE, typename OPT> inline bool
ODEBND_BASE<T,PMT,PVT>::_RHS_PM_FSA_STA
( const OPT&options, double t, const REALTYPE*x, REALTYPE*xdot )
{
  if( !_pRHS ) return false;

  switch( options.WRAPMIT){
  case OPT::NONE:
    _vec2PMI( x, _PMenv, _nx, _PMx );   // set current state polynomial model
#ifdef MC__ODEBND_BASE_DINEQPM_DEBUG
    _print_interm( t, _nx, _PMx, "PMx Intermediate", std::cerr );
#endif
    _PMVAR[_nx+_npar] = t; // set current time
    _QUAD_PM( _pDAG, _opQUAD, _PMRHS, _nq, _pQUAD, _nVAR-_nq, _pVAR, _PMVAR,
              _PMqdot );
    _RHS_PM_NONE( _pDAG, _opRHS, _PMRHS, _nx, _pRHS, _nVAR-_nq, _pVAR, _PMVAR,
                  _PMxdot );
    // Whether or not to ignore the remainder
    if( !options.PMNOREM )
      _PMI2vec( _PMenv, _nx, _PMxdot, xdot, false );
    else
      _PMI2vec( _PMenv, _nx, _PMxdot, 0, xdot );
#ifdef MC__ODEBND_BASE_DINEQPM_DEBUG
    _print_interm( t, _nx, _PMxdot, "PMxdot Intermediate", std::cerr );
    {int dum; std::cin >> dum;}
#endif
    return true;  

  case OPT::DINEQ:
    _vec2PMI( x, _PMenv, _nx, _PMx );   // set current state polynomial model
    for( unsigned ix=0; ix<_nx; ix++ ) _PMx[ix].center();
#ifdef MC__ODEBND_BASE_DINEQPM_DEBUG
    _print_interm( t, _nx, _PMx, "PMx Intermediate", std::cerr );
#endif
    _PMVAR[_nx+_npar] = t; // set current time
    _QUAD_PM( _pDAG, _opQUAD, _PMRHS, _nq, _pQUAD, _nVAR-_nq, _pVAR, _PMVAR,
              _PMqdot );
    _RHS_PM_DI( _pDAG, _opRHSi, _PMRHS, _nx, _pRHS, _nVAR-_nq, _pVAR, _PMVAR,
                _PMxdot, _RxLdot, _RxUdot );
    // Whether or not to ignore the remainder
    if( !options.PMNOREM )
      _PMI2vec( _PMenv, _nx, _PMxdot, _RxLdot, _RxUdot, xdot );
    else
      _PMI2vec( _PMenv, _nx, _PMxdot, 0, xdot );
#ifdef MC__ODEBND_BASE_DINEQPM_DEBUG
    _print_interm( t, _nx, _PMxdot, "PMxdot Intermediate", std::cerr );
    {int dum; std::cin >> dum;}
#endif
    return true;  
   
  case OPT::ELLIPS:
  default:{
    _vec2PME( x, _PMenv, _nx, _PMx, _Q, _Er, _Ir ); // set current state polynomial model
    _PMVAR[_nx+_npar] = t; // set current time   
#ifdef MC__ODEBND_BASE_DINEQPM_DEBUG
    _print_interm( t, _nx, _PMx, _Er, "PMx Intermediate", std::cerr );
#endif
    _QUAD_PM( _pDAG, _opQUAD, _PMRHS, _nq, _pQUAD, _nVAR-_nq, _pVAR, _PMVAR,
              _PMqdot );

    // In this variant a bound on the Jacobian matrix is computed and the
    // linear part is taken as the mid-point of this matrix
    if( !options.ORDMIT ){
      for( unsigned ix=0; ix<_nx; ix++ ){
        _IVAR[ix] = _PMx[ix].bound(); // set current state bounds
        _PMVAR[ix].center().set( T(0.) ); // cancel remainder term
      }
      _IVAR[_nx+_npar] = t; // set current time
      _RHS_PM_ELL0( _pDAG, _opRHS, _PMRHS, _opJAC, _IJAC, _nx, _pRHS,
                    _pJAC, _nVAR-_nq, _pVAR, _PMVAR, _IVAR, _PMxdot, _Idfdx,
                    _Ir, _A, _Irdot );
    }

    // In this variant a polynomial model of the Jacobian matrix is computed and the
    // linear part is taken as the mid-point of this matrix
    else if( _PMenv->nord() >= _MVXPenv->nord() ){
#ifdef MC__ODEBND_BASE_MVXP_USE
      for( unsigned jx=0; jx<_nx; jx++ ){
        _MVXPd[jx].set( _MVXPenv, _npar+jx, _Ir[jx] );
        _MVXPx[jx].set( _MVXPenv ).set( _PMx[jx].center().set( T(0.) ), true );
      }
      _e2x( _nx, _MVXPd, _MVXPx, false );
#else
      for( unsigned jx=0; jx<_nx; jx++ ){
        _MVXPx[jx].set( _MVXPenv ).set( _PMx[jx].center(), true );
        _PMx[jx].set( T(0.) );
      }
#endif
#ifdef MC__ODEBND_BASE_DINEQPM_DEBUG
      _print_interm( t, _nx, _MVXPx, "MVXPx Intermediate", std::cerr );
#endif
      _MVXPVAR[_nx+_npar] = t; // set current time
      _RHS_PM_ELL1( _pDAG, _opRHS, _PMRHS, _opJAC, _PMJAC, _nx, _pRHS,
                    _pJAC, _nVAR-_nq, _pVAR, _PMVAR, _MVXPVAR, _PMxdot,
                    _MVXPdfdx, _Ir, _A, _Irdot );
    }

    // In this variant a polynomial model in the joint state-parameter and
    // of the same order as the parameter polynomial model is computed
    else{
      for( unsigned jx=0; jx<_nx; jx++ ){
        _MVXPd[jx].set( _MVXPenv, _npar+jx, _Ir[jx] );
        _MVXPx[jx].set( _MVXPenv ).set( _PMx[jx].center().set( T(0.) ), true );
      }
      _e2x( _nx, _MVXPd, _MVXPx, false );
#ifdef MC__ODEBND_BASE_DINEQPM_DEBUG
      _print_interm( t, _nx, _MVXPx, "MVXPx Intermediate", std::cerr );
#endif
      _MVXPVAR[_nx+_npar] = t; // set current time   
      _RHS_PM_ELL2( _pDAG, _opRHS, _PMRHS, _nx, _pRHS, _nVAR-_nq, _pVAR,
                    _MVXPVAR, _PMenv, _PMxdot, _MVXPf, _npar, _Ir, _A, _Irdot );
    }

    _RHS_PM_ELL( _nx, _Q, _A, _Irdot, _Qdot, options.QTOL, machprec() );
    // Whether or not to ignore the remainder
    if( !options.PMNOREM )
      _PME2vec( _PMenv, _nx, _PMxdot, _Qdot, xdot );
    else
      _PME2vec( _PMenv, _nx, _PMxdot, 0, xdot );
#ifdef MC__ODEBND_BASE_DINEQPM_DEBUG
    _print_interm( t, _nx, _PMxdot, "PMxdot Intermediate", std::cerr );
    {int dum; std::cin >> dum;}
#endif
    return true;
   }
  }
}

template <typename T, typename PMT, typename PVT>
template <typename REALTYPE, typename OPT> inline bool
ODEBND_BASE<T,PMT,PVT>::_RHS_PM_FSA_QUAD
( const OPT&options, double t, const REALTYPE*x, REALTYPE*qdot,
  const bool reinit )
{
  if( !_pQUAD ) return false;

  switch( options.WRAPMIT){
  case OPT::NONE:
  case OPT::DINEQ:
    if( !reinit ) break;
    _vec2PMI( x, _PMenv, _nx, _PMx );   // set current state polynomial model
#ifdef MC__ODEBND_BASE_DINEQPM_DEBUG
    _print_interm( t, _nx, _PMx, "PMx Intermediate", std::cerr );
#endif
    break;
   
  case OPT::ELLIPS:
  default:
    if( !reinit ) break;
    _vec2PME( x, _PMenv, _nx, _PMx, _Q, _Er, _Ir ); // set current state polynomial model
#ifdef MC__ODEBND_BASE_DINEQPM_DEBUG
    _print_interm( t, _nx, _PMx, _Er, "PMx Intermediate", std::cerr );
#endif
    break;
  }

  _PMVAR[_nx+_npar] = t; // set current time
  _QUAD_PM( _pDAG, _opQUAD, _PMRHS, _nq, _pQUAD, _nVAR-_nq, _pVAR, _PMVAR,
            _PMqdot );
#ifdef MC__ODEBND_BASE_DINEQPM_DEBUG
    _print_interm( t, _nq, _PMqdot, "PMqdot Intermediate", std::cerr );
    {int dum; std::cin >> dum;}
#endif

  // Whether or not to ignore the remainder
  if( !options.PMNOREM )
    _PMI2vec( _PMenv, _nq, _PMqdot, qdot, true );
  else
    _PMI2vec( _PMenv, _nq, _PMqdot, 0, qdot );

  return true;  
}

template <typename T, typename PMT, typename PVT>
template <typename REALTYPE, typename OPT> inline bool
ODEBND_BASE<T,PMT,PVT>::_JAC_PM_FSA_STA
( const OPT&options, double t, const REALTYPE*x, REALTYPE*jac, REALTYPE*xdot )
{
  return false;
}

template <typename T, typename PMT, typename PVT> inline bool
ODEBND_BASE<T,PMT,PVT>::_FCT_PM_FSA_STA
( const unsigned iFCT, const double t, PVT*PMf )
{
  if( !_nf || !PMf ) return true;

  _PMVAR[_nx+_npar] = t; // set current time
  const FFVar* pFCT = _vFCT.at( iFCT );
  _pDAG->eval( _nf, pFCT, PMf, _nVAR, _pVAR, _PMVAR, iFCT?true:false );
  return true;
}

template <typename T, typename PMT, typename PVT>
template <typename OPT> inline bool
ODEBND_BASE<T,PMT,PVT>::_SET_PM_FSA_STA
( const OPT&options, const unsigned iRHS, const unsigned iQUAD )
{
  if( _vRHS.size() <= iRHS ) return false;
  if( _nq && _vQUAD.size() <= iQUAD ) return false;

  _pRHS   = _vRHS.at( iRHS );
  _opRHS  = _pDAG->subgraph( _nx, _pRHS );
  _pQUAD  = _nq? _vQUAD.at( iQUAD ): 0;
  _opQUAD = _pDAG->subgraph( _nq, _pQUAD );
  const unsigned opmax = _opRHS.size()>_opQUAD.size()?_opRHS.size():_opQUAD.size();

  delete[] _IRHS;   _IRHS = 0;
  delete[] _PMRHS;  _PMRHS = 0;
  delete[] _opRHSi; _opRHSi = 0;
  delete[] _IJAC;   _IJAC = 0;
  delete[] _PMJAC;  _PMJAC = 0;
  delete[] _pJAC;   _pJAC = 0; _opJAC.clear();

  _PMRHS = new PVT[ opmax ];

  switch( options.WRAPMIT){
  case OPT::NONE: default:
    break;

  case OPT::DINEQ:
    _opRHSi = new std::list<const FFOp*>[_nx];
    for( unsigned ix=0; ix<_nx; ix++ )
      _opRHSi[ix] = _pDAG->subgraph( 1, _pRHS+ix );
    break;

  case OPT::ELLIPS:
    _pJAC = _pDAG->FAD( _nx, _pRHS, _nx, _pVAR );
    _opJAC = _pDAG->subgraph( _nx*_nx, _pJAC );
    if( !options.ORDMIT )
      _IJAC = new T[ _opJAC.size() ];
    else if( _PMenv->nord() >= _MVXPenv->nord() )
      _PMJAC = new PVT[ _opJAC.size() ];
    break;
  }

  return true;
}

template <typename T, typename PMT, typename PVT>
template <typename OPT> inline bool
ODEBND_BASE<T,PMT,PVT>::_INI_PM_FSA
( const OPT&options, const unsigned np, const PVT* PMp, const unsigned ns )
{
  if( !_INI_PM_FSA( options, np, PMp, ns ) ) return false;

  // Size state/quadrature sensitivity arrays
  delete[] _PMxp;    _PMxp    = new PVT[_npar*_nx];
  delete[] _Rxp;     _Rxp     = new T[_npar*_nx];
  delete[] _PMqp;    _PMqp    = _nq? new PVT[_npar*_nq];
  delete[] _Rqp;     _Rqp     = _nq? new T[_npar*_nq];
  delete[] _PMxpdot; _PMxpdot = new PVT[_npar*_nx];
  delete[] _Rxpdot;  _Rxpdot  = new PVT[_npar*_nx];
  delete[] _PMqpdot; _PMqpdot = _nq? new PVT[_npar*_nq]: 0;
  delete[] _Rqpdot;  _Rqpdot  = _nq? new double[_npar*_nq]: 0;

  return true;
}

template <typename T, typename PMT, typename PVT>
template <typename REALTYPE, typename OPT> inline bool
ODEBND_BASE<T,PMT,PVT>::_IC_PM_FSA_QUAD
( const OPT&options, REALTYPE*vec )
{
  if( !_vQUAD.size() || !_nq ) return true;
  for( unsigned iq=0; iq<_nq; iq++ )
    _PMq[iq] = 0.;
  _PMI2vec( _PMenv, _nq, _PMq, vec, true );
  for( unsigned iqp=0; iqp<_nq*_npar; iqp++ ){
    _PMqp[iq] = 0.;
    _Rqp[iqp] = _PMqp[iq].R();
  }
  _I2vec( _nq, _Rqp, vec+_PMenv->nmon()+_nq, true );
  return true;
}

} // end namescape mc

#endif

