// Copyright (C) 2015 Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

#ifndef MC__IODEBND_BASE_HPP
#define MC__IODEBND_BASE_HPP

#undef  MC__IODEBND_BASE_DINEQI_DEBUG
#undef  MC__IODEBND_BASE_DINEQPM_DEBUG

#include "odebnd_base.hpp"
#include "aebnd.hpp"

namespace mc
{
//! @brief C++ base class for computing enclosures of the reachable set of parametric IODEs or DAEs using continuous-time set-valued integration.
////////////////////////////////////////////////////////////////////////
//! mc::IODEBND_GSL is a C++ base class that computes enclosures of the
//! reachable set of parametric implicit ordinary differential equations
//! (IODEs) or differential-differential equations (DAEs) using
//! continuous-time set-valued integration. The methods for DAEs relies
//! on the underlying ODEs.
////////////////////////////////////////////////////////////////////////
template <typename T, typename PMT=mc::TModel<T>, typename PVT=mc::TVar<T> >
class IODEBND_BASE:
  public virtual BASE_DE,
  public virtual ODEBND_BASE<T,PMT,PVT>
{
  typedef Ellipsoid E;
  typedef BASE_DE::STATUS STATUS;

protected:
  using ODEBND_BASE<T,PMT,PVT>::_npar;
  using ODEBND_BASE<T,PMT,PVT>::_nVAR;
  using ODEBND_BASE<T,PMT,PVT>::_pVAR;
  using ODEBND_BASE<T,PMT,PVT>::_pRHS;
  using ODEBND_BASE<T,PMT,PVT>::_pJAC;
  using ODEBND_BASE<T,PMT,PVT>::_pQUAD;

  using ODEBND_BASE<T,PMT,PVT>::_opRHS;
  using ODEBND_BASE<T,PMT,PVT>::_opRHSi;
  using ODEBND_BASE<T,PMT,PVT>::_IRHS;
  using ODEBND_BASE<T,PMT,PVT>::_PMRHS;
  using ODEBND_BASE<T,PMT,PVT>::_opJAC;
  using ODEBND_BASE<T,PMT,PVT>::_IJAC;
  using ODEBND_BASE<T,PMT,PVT>::_PMJAC;

  using ODEBND_BASE<T,PMT,PVT>::_IVAR;
  using ODEBND_BASE<T,PMT,PVT>::_It;
  using ODEBND_BASE<T,PMT,PVT>::_A;
  using ODEBND_BASE<T,PMT,PVT>::_Q;
  using ODEBND_BASE<T,PMT,PVT>::_Qdot;
  using ODEBND_BASE<T,PMT,PVT>::_Ir;
  using ODEBND_BASE<T,PMT,PVT>::_Irdot;
  using ODEBND_BASE<T,PMT,PVT>::_Er;
  using ODEBND_BASE<T,PMT,PVT>::_pref;
  using ODEBND_BASE<T,PMT,PVT>::_Ip;
  using ODEBND_BASE<T,PMT,PVT>::_B;
  using ODEBND_BASE<T,PMT,PVT>::_Bdot;
  using ODEBND_BASE<T,PMT,PVT>::_xref;
  using ODEBND_BASE<T,PMT,PVT>::_xrefdot;
  using ODEBND_BASE<T,PMT,PVT>::_Ix;
  using ODEBND_BASE<T,PMT,PVT>::_Ixdot;
  using ODEBND_BASE<T,PMT,PVT>::_xLdot;
  using ODEBND_BASE<T,PMT,PVT>::_xUdot;
  using ODEBND_BASE<T,PMT,PVT>::_MVXPenv;
  using ODEBND_BASE<T,PMT,PVT>::_MVXPVAR;
  using ODEBND_BASE<T,PMT,PVT>::_MVXPt;
  using ODEBND_BASE<T,PMT,PVT>::_MVXPx;
  using ODEBND_BASE<T,PMT,PVT>::_MVXPp;
  using ODEBND_BASE<T,PMT,PVT>::_MVXPd;
  using ODEBND_BASE<T,PMT,PVT>::_MVXPf;
  using ODEBND_BASE<T,PMT,PVT>::_MVPt;
  using ODEBND_BASE<T,PMT,PVT>::_MVPx;
  using ODEBND_BASE<T,PMT,PVT>::_MVPp;
  using ODEBND_BASE<T,PMT,PVT>::_vec2I;
  using ODEBND_BASE<T,PMT,PVT>::_vec2E;
  using ODEBND_BASE<T,PMT,PVT>::_ep2x;
  using ODEBND_BASE<T,PMT,PVT>::_E2vec;
  using ODEBND_BASE<T,PMT,PVT>::_I2vec;
  using ODEBND_BASE<T,PMT,PVT>::_IC_I_ELL;
  using ODEBND_BASE<T,PMT,PVT>::_CC_I_ELL;
  using ODEBND_BASE<T,PMT,PVT>::_RHS_I_ELL;

  using ODEBND_BASE<T,PMT,PVT>::_PMenv;
  using ODEBND_BASE<T,PMT,PVT>::_PMVAR;
  using ODEBND_BASE<T,PMT,PVT>::_PMt;
  using ODEBND_BASE<T,PMT,PVT>::_PMx;
  using ODEBND_BASE<T,PMT,PVT>::_PMxdot;
  using ODEBND_BASE<T,PMT,PVT>::_RxLdot;
  using ODEBND_BASE<T,PMT,PVT>::_RxUdot;
  using ODEBND_BASE<T,PMT,PVT>::_PMp;
  using ODEBND_BASE<T,PMT,PVT>::_Idfdx;
  using ODEBND_BASE<T,PMT,PVT>::_MVXPdfdx;
  using ODEBND_BASE<T,PMT,PVT>::_vec2PMI;
  using ODEBND_BASE<T,PMT,PVT>::_vec2PME;
  using ODEBND_BASE<T,PMT,PVT>::_PMI2vec;
  using ODEBND_BASE<T,PMT,PVT>::_PME2vec;
  using ODEBND_BASE<T,PMT,PVT>::_IC_PM_ELL;
  using ODEBND_BASE<T,PMT,PVT>::_CC_PM_ELL;
  using ODEBND_BASE<T,PMT,PVT>::_RHS_PM_ELL0;
  using ODEBND_BASE<T,PMT,PVT>::_RHS_PM_ELL1;
  using ODEBND_BASE<T,PMT,PVT>::_RHS_PM_ELL2;
  using ODEBND_BASE<T,PMT,PVT>::_RHS_PM_ELL;
  using ODEBND_BASE<T,PMT,PVT>::_e2x;

  //! @brief AE bounder for nonlinear algebraic equations
  AEBND<T,PMT,PVT>* _ICBND;

  //! @brief AE bounder for time derivatives of underlying ODEs
  AEBND<T,PMT,PVT>* _RHSBND;

  //! @brief const pointer to AE function in current stage of IODE or DAE **DO NOT FREE**
  const FFVar* _pAE;

  //! @brief const pointer to underlying IC function in current stage of IODE or DAE
  FFVar* _pUIC;

  //! @brief const pointer to underlying RHS function in current stage of IODE or DAE
  FFVar* _pURHS;

  //! @brief const pointer to underlying residualfunction in current stage of IODE or DAE
  FFVar* _pURES;

  //! @brief vector of a priori state bounds at initial/stage times
  std::vector<const T*> _vIx0;

  //! @brief vector of a priori state bounds at initial/stage times
  std::vector<const PVT*> _vPMx0;

public:
  //! @brief Default constructor
  IODEBND_BASE();

  //! @brief Virtual destructor
  virtual ~IODEBND_BASE();

  //! @brief Define a priori interval bounds for the state variables at initial time
  void set_apriori
    ( const T*const Ix0 )
    { _vIx0.clear(); _vIx0.push_back( Ix0 ); }

  //! @brief Define a priori interval bounds for the state variables at initial/stage times
  void set_apriori
    ( const unsigned ns, const T*const Ix0 )
    { _vIx0.clear(); for( unsigned i=0; i<ns; i++ ) _vIx0.push_back( Ix0+i*_nx ); }

  //! @brief Define a priori polynomial bounds for the state variables at initial times
  void set_apriori
    ( const PVT*const PMx0 )
    { _vPMx0.clear(); _vPMx0.push_back( PMx0 ); }

  //! @brief Define a priori polynomial bounds for the state variables at initial/stage times
  void set_apriori
    ( const unsigned ns, const PVT*const PMx0 )
    { _vPMx0.clear(); for( unsigned i=0; i<ns; i++ ) _vPMx0.push_back( PMx0+i*_nx ); }

protected:

  //! @brief Function to set RHS and QUAD pointers to underlying ODEs
  template <typename OPT> bool _URHS_STA
    ( const OPT&options, const unsigned iRHS, const unsigned iQUAD );

  //! @brief Function to initialize state interval bounding
  template <typename OPT> bool _INI_I_STA
    ( const OPT&options, const unsigned np, const T*Ip, const unsigned ns );

  //! @brief Function to initialize state interval bounds
  template <typename REALTYPE, typename OPT> bool _IC_I_STA
    ( const OPT&options, const double t, REALTYPE*vec );

  //! @brief Function to reinitialize state interval bounds
  template <typename REALTYPE, typename OPT> bool _CC_I_STA
    ( const OPT&options, const unsigned iIC, const double t, REALTYPE*vec );

  //! @brief Function to calculate the RHS of auxiliary ODEs in interval arithmetic
  template <typename REALTYPE, typename OPT> bool _RHS_I_STA
    ( const OPT&options, double t, const REALTYPE*x, REALTYPE*xdot );

  //! @brief Function to initialize GSL for state polynomial models
  template <typename OPT> bool _INI_PM_STA
    ( const OPT&options, const unsigned np, const PVT*PMp, const unsigned ns );

  //! @brief Function to initialize state polynomial models
  template <typename REALTYPE, typename OPT> bool _IC_PM_STA
    ( const OPT&options, const double t, REALTYPE*vec );

  //! @brief Function to reinitialize state polynomial bounds
  template <typename REALTYPE, typename OPT> bool _CC_PM_STA
    ( const OPT&options, const unsigned iIC, const double t, REALTYPE*vec );

  //! @brief Function to calculate the RHS of auxiliary ODEs in polynomial mode arithmetic
  template <typename REALTYPE, typename OPT> bool _RHS_PM_STA
    ( const OPT&options, double t, const REALTYPE*x, REALTYPE*xdot );

  //! @brief Private methods to block default compiler methods
  IODEBND_BASE(const IODEBND_BASE&);
  IODEBND_BASE& operator=(const IODEBND_BASE&);
};

template <typename T, typename PMT, typename PVT> inline
IODEBND_BASE<T,PMT,PVT>::IODEBND_BASE
()
: BASE_DE(), ODEBND_BASE<T,PMT,PVT>(), _pAE(0), _pUIC(0), _pURHS(0), _pURES(0)
{
  // Initalize implicit equation bounders
  _ICBND = new AEBND<T,PMT,PVT>();
  _RHSBND = new AEBND<T,PMT,PVT>();
}

template <typename T, typename PMT, typename PVT> inline
IODEBND_BASE<T,PMT,PVT>::~IODEBND_BASE
()
{
  delete[] _pUIC;
  delete[] _pURHS;
  delete[] _pURES;
  delete _ICBND;
  delete _RHSBND;
}

template <typename T, typename PMT, typename PVT>
template <typename OPT> inline bool
IODEBND_BASE<T,PMT,PVT>::_URHS_STA
( const OPT&options, const unsigned iRHS, const unsigned iQUAD )
{
  if( _vRHS.size() <= iRHS || (_na && _vAE.size() <= iRHS) )
    return false;

  // Create residuals for underlying ODEs
  delete[] _pURHS; _pURHS = new FFVar[_nx];
  delete[] _pURES; _pURES = new FFVar[_nx];
  for( unsigned id=0; id<_nd; id++ )
    _pURES[id] = _vRHS.at( iRHS )[id]; // _pDX[id] already part of IODE residual
  if( _na ){
    delete[] _pJAC; _pJAC = 0;
    _pAE = _vAE.at( iRHS );
    try{ _pJAC = _pDAG->FAD( _na, _pAE, _nx, _pVAR ); }
    catch(...){ return false; } 
    // ordering in pJac: [ dAE1/dx1 | ... | dAE1/dxNX | dAE2/dx1 | ... | dAE2/dxNX | ... ] 
    for( unsigned ia=0; ia<_na; ia++ ){
      _pURES[_nd+ia] = 0.;
      for( unsigned jd=0; jd<_nd; jd++ )
        _pURES[_nd+ia] += _pJAC[ia*_nx+jd] * _pDX[jd];
      for( unsigned ja=0; ja<_na; ja++ )
        _pURES[_nd+ia] += _pJAC[ia*_nx+(_nd+ja)] * _pDX[_nd+ja];
    }
  }

  // Set implicit equation bounder for underlying ODEs
  _RHSBND->set_dag( _pDAG );
  _RHSBND->set_var( _nx+_npar+1, _pVAR );
  _RHSBND->set_dep( _nx, _pDX, _pURES );
  _RHSBND->options = options.RHSBNDOPT;
  if( _RHSBND->setup()!=AEBND<T,PMT,PVT>::NORMAL
   || _RHSBND->solve( _pURHS )!=AEBND<T,PMT,PVT>::NORMAL )
    return false;
  _pRHS = _pURHS;
#ifdef MC__IODEBND_BASE_DEBUG
  std::ofstream o_sol( "UODE_RHS.dot", std::ios_base::out );
  _pDAG->dot_script( _nx, _pURHS, o_sol );
  o_sol.close();
  throw(1);
#endif

  // Quadrature equations
  if( _nq && _vQUAD.size() <= iQUAD ) return false;
  _pQUAD  = _nq? _vQUAD.at( iQUAD ): 0;

  return true;
}

template <typename T, typename PMT, typename PVT>
template <typename OPT> inline bool
IODEBND_BASE<T,PMT,PVT>::_INI_I_STA
( const OPT&options, const unsigned np, const T*Ip, const unsigned ns )
{
  if( _na+_nd != _nx
   || !ODEBND_BASE<T,PMT,PVT>::_INI_I_STA( options, np, Ip, ns ) )
    return false;
  return true;
}

template <typename T, typename PMT, typename PVT>
template <typename REALTYPE, typename OPT> inline bool
IODEBND_BASE<T,PMT,PVT>::_IC_I_STA
( const OPT&options, const double t, REALTYPE*vec )
{
  // Setup AE bounder
  if( !_vIC.size() || (_na && !_vAE.size()) || _nx0 != _nd )
    return false;
  delete[] _pUIC; _pUIC = new FFVar[_nx];
  for( unsigned i=0; i<_nx0; i++ ) 
    _pUIC[i] = _vIC.at(0)[i];
  for( unsigned i=0; i<_na; i++ )
    _pUIC[_nx0+i] = _vAE.at(0)[i];
  _ICBND->set_dag( _pDAG );
  _ICBND->set_var( _npar+1, _pVAR+_nx );
  _ICBND->set_dep( _nx, _pVAR, _pUIC );
  _ICBND->options = options.ICBNDOPT;
  _ICBND->setup();
  const T*Ix0 = _vIx0.size()? _vIx0.at(0): 0;

  switch( options.WRAPMIT){

  case OPT::NONE:
  case OPT::DINEQ:
    *_It = t; // current time
    if( _ICBND->solve( _Ip, _Ix, Ix0 ) != AEBND<T,PMT,PVT>::NORMAL )
      return false;
    _I2vec( _nx, _Ix, vec );
    break;

  case OPT::ELLIPS:
  default:
    *_MVPt = t; // current time
    if( _ICBND->solve( _MVPp, _MVPx, Ix0 ) != AEBND<T,PMT,PVT>::NORMAL )
      return false;
    _IC_I_ELL( _nx, _MVPx, _xref, _Q, _npar, _B, _Ix );
    _E2vec( _nx, _npar, _xref, _Q, _B, vec );
    break;
  }

  return true;
}

template <typename T, typename PMT, typename PVT>
template <typename REALTYPE, typename OPT> inline bool
IODEBND_BASE<T,PMT,PVT>::_CC_I_STA
( const OPT&options, const unsigned iRHS, const double t, REALTYPE*vec )
{
  // Only handles discontinuity in algebraic states
  if( (_na && _vAE.size() <= iRHS) || _nx0 != _nd )
    return false;
  _pAE = _vAE.at(iRHS);
  std::vector<FFVar> vVAR( _nd+_npar+1 );
  for( unsigned i=0; i<_nd; i++ ) vVAR[i] = _pX[i];
  for( unsigned i=0; i<_npar+1; i++ ) vVAR[_nd+i] = _pVAR[_nx+i];
  _ICBND->set_dag( _pDAG );
  _ICBND->set_var( _nd+_npar+1, vVAR.data() );
  _ICBND->set_dep( _na, _pVAR+_nd, _pAE );
  _ICBND->options = options.ICBNDOPT;
  _ICBND->setup();
  const T*Ixa0 = _vIx0.size() >= iRHS? _vIx0.at(iRHS)+_nd:
                (_vIx0.size()? _vIx0.at(0)+_nd: 0 );

  switch( options.WRAPMIT){

  case OPT::NONE:
  case OPT::DINEQ:{
    _vec2I( vec, _nx, _Ix ); // current state bounds
    *_It = t; // current time
    std::vector<T> IVAR( _nd+_npar+1 );
    for( unsigned i=0; i<_nd; i++ ) IVAR[i] = _Ix[i];
    for( unsigned i=0; i<_npar+1; i++ ) IVAR[_nd+i] = _IVAR[_nx+i];
    if( _ICBND->solve( IVAR.data(), _Ix+_nd, Ixa0 ) != AEBND<T,PMT,PVT>::NORMAL )
      return false;
    _I2vec( _nx, _Ix, vec );
    break;
   }
  case OPT::ELLIPS:
  default:{
    _vec2E( vec, _nx, _npar, _Q, _Er, _Ir, _pref, _Ip, _B, _xref, _Ix ); // current state enclosure
    for( unsigned jx=0; jx<_nx; jx++ )
      _MVXPd[jx].set( _MVXPenv, jx, _Ir[jx] );
    _ep2x( _nx, _npar, _MVXPd, _pref, _MVXPp, _B, _xref, _MVXPx );
    *_MVXPt = t; // current time
    std::vector<PVT> PVAR( _nd+_npar+1 );
    for( unsigned i=0; i<_nd; i++ ) PVAR[i] = _MVXPx[i];
    for( unsigned i=0; i<_npar+1; i++ ) PVAR[_nd+i] = _MVXPVAR[_nx+i];
    if( _ICBND->solve( PVAR.data(), _MVXPx+_nd, Ixa0 ) != AEBND<T,PMT,PVT>::NORMAL )
      return false;
    _CC_I_ELL( _nx, _MVXPx, _Er, _A, _npar, _xrefdot, _Bdot, _Irdot, _Qdot,
               options.QTOL, machprec() );
    _E2vec( _nx, _npar, _xrefdot, _Qdot, _Bdot, vec );
    break;
   }
  }

  return true;
}

template <typename T, typename PMT, typename PVT>
template <typename REALTYPE, typename OPT> inline bool
IODEBND_BASE<T,PMT,PVT>::_RHS_I_STA
( const OPT&options, double t, const REALTYPE*x, REALTYPE*xdot )
{
  switch( options.WRAPMIT){
  case OPT::NONE:
    _vec2I( x, _nx, _Ix );   // set current state bounds
    *_It = t; // set current time
    // Solve residuals for underlying ODEs
    if( options.RHSNUMER ){
      if( _RHSBND->solve( _IVAR, _Ixdot ) != AEBND<T,PMT,PVT>::NORMAL )
        return false;
    }
    else
      _pDAG->eval( _opRHS, _IRHS, _nx, _pRHS, _Ixdot, _nVAR-_nq, _pVAR, _IVAR );
    _I2vec( _nx, _Ixdot, xdot );
    return true;
   
  case OPT::DINEQ:
    _vec2I( x, _nx, _Ix );   // set current state bounds
    *_It = t; // set current time
    for( unsigned ix=0; ix<_nx; ix++ ){
      T Ixi = _IVAR[ix];
      for( unsigned up=0; up<2; up++ ){ // separate lower/upper bounding subproblems
        _IVAR[ix] = up? Op<T>::u( Ixi ): Op<T>::l( Ixi );
        if( options.RHSNUMER ){
          if( _RHSBND->solve( _IVAR, _Ixdot ) != AEBND<T,PMT,PVT>::NORMAL )
            { _IVAR[ix] = Ixi; return false; }
        }
        else
          _pDAG->eval( _opRHSi[ix], _IRHS, 1, _pRHS+ix, _Ixdot+ix, _nVAR-_nq,
                       _pVAR, _IVAR );
        if( up ) _xUdot[ix] = Op<T>::u( _Ixdot[ix] );
        else     _xLdot[ix] = Op<T>::l( _Ixdot[ix] );
      }
      _IVAR[ix] = Ixi;
    }
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

    // Setup polynomial model expansion of RHS
    *_MVXPt = t; // set current time
    for( unsigned jx=0; jx<_nx; jx++ )
      _MVXPd[jx].set( _MVXPenv, jx, _Ir[jx] );
    _ep2x( _nx, _npar, _MVXPd, _pref, _MVXPp, _B, _xref, _MVXPx );
#ifdef MC__ODEBND_BASE_DINEQI_DEBUG
    for( unsigned ix=0; ix<_nx; ix++ )
      std::cout << "MVXPx[ " << ix << "] = " << _MVXPx[ix] << std::endl;
    //{ int dum; std::cin >> dum; }
#endif

    // Construct the ellipsoidal remainder derivatives
    if( options.RHSNUMER ){
      if( _RHSBND->solve( _MVXPVAR, _MVXPf ) != AEBND<T,PMT,PVT>::NORMAL )
        return false;
    }
    else
      _pDAG->eval( _opRHS, _PMRHS, _nx, _pRHS, _MVXPf, _nVAR-_nq, _pVAR, _MVXPVAR );
    if( options.QSCALE )
      _RHS_I_ELL( _nx, _MVXPf, _Q, _A, _npar, _xrefdot, _Bdot, _Irdot,
         _Qdot, options.QTOL, machprec(), _Ix );
    else
      _RHS_I_ELL( _nx, _MVXPf, _Q, _A, _npar, _xrefdot, _Bdot, _Irdot,
         _Qdot, options.QTOL, machprec() );
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
template <typename OPT> inline bool
IODEBND_BASE<T,PMT,PVT>::_INI_PM_STA
( const OPT&options, const unsigned np, const PVT*PMp, const unsigned ns )
{
  if( _na+_nd != _nx
   || !ODEBND_BASE<T,PMT,PVT>::_INI_PM_STA( options, np, PMp, ns ) )
    return false;
  return true;
}

template <typename T, typename PMT, typename PVT>
template <typename REALTYPE, typename OPT> inline bool
IODEBND_BASE<T,PMT,PVT>::_IC_PM_STA
( const OPT&options, const double t, REALTYPE*vec )
{
  // Setup AE bounder
  if( !_vIC.size() || (_na && !_vAE.size()) || _nx0 != _nd )
    return false;
  delete[] _pUIC; _pUIC = new FFVar[_nx];
  for( unsigned i=0; i<_nx0; i++ ) 
    _pUIC[i] = _vIC.at(0)[i];
  for( unsigned i=0; i<_na; i++ )
    _pUIC[_nx0+i] = _vAE.at(0)[i];
  _ICBND->set_dag( _pDAG );
  _ICBND->set_var( _npar+1, _pVAR+_nx );
  _ICBND->set_dep( _nx, _pVAR, _pUIC );
  _ICBND->options = options.ICBNDOPT;
  _ICBND->setup();
  const T*Ix0 = _vIx0.size()? _vIx0.at(0): 0;

  *_PMt = t; // current time
  if( _ICBND->solve( _PMp, _PMx, Ix0 ) != AEBND<T,PMT,PVT>::NORMAL )
    return false;

  switch( options.WRAPMIT){
  case OPT::NONE:
    // Whether or not to ignore the remainder
    if( !options.PMNOREM )
      _PMI2vec( _PMenv, _nx, _PMx, vec, true );
    else
      _PMI2vec( _PMenv, _nx, _PMx, 0, vec );
    break;

  case OPT::DINEQ:
    // Whether or not to ignore the remainder
    if( !options.PMNOREM )
      _PMI2vec( _PMenv, _nx, _PMx, vec, false );
    else
      _PMI2vec( _PMenv, _nx, _PMx, 0, vec );
    break;

  case OPT::ELLIPS:
  default:
    _IC_PM_ELL( _nx, _PMx, _Q, _Er, _Ir );
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
IODEBND_BASE<T,PMT,PVT>::_CC_PM_STA
( const OPT&options, const unsigned iRHS, const double t, REALTYPE*vec )
{
  // Only handles discontinuity in algebraic states
  if( (_na && _vAE.size() <= iRHS) || _nx0 != _nd )
    return false;
  _pAE = _vAE.at(iRHS);
  std::vector<FFVar> vVAR( _nd+_npar+1 );
  for( unsigned i=0; i<_nd; i++ ) vVAR[i] = _pX[i];
  for( unsigned i=0; i<_npar+1; i++ ) vVAR[_nd+i] = _pVAR[_nx+i];
  _ICBND->set_dag( _pDAG );
  _ICBND->set_var( _nd+_npar+1, vVAR.data() );
  _ICBND->set_dep( _na, _pVAR+_nd, _pAE );
  _ICBND->options = options.ICBNDOPT;
  _ICBND->setup();
  //const T*Ixa0 = _vIx0.size() >= iRHS? _vIx0.at(iRHS)+_nd:
  //              (_vIx0.size()? _vIx0.at(0)+_nd: 0 );

  return false;
}

template <typename T, typename PMT, typename PVT>
template <typename REALTYPE, typename OPT> inline bool
IODEBND_BASE<T,PMT,PVT>::_RHS_PM_STA
( const OPT&options, double t, const REALTYPE*x, REALTYPE*xdot )
{
  switch( options.WRAPMIT ){
  case OPT::NONE:
    _vec2PMI( x, _PMenv, _nx, _PMx, true );   // set current state polynomial model
#ifdef MC__ODEBND_BASE_DINEQPM_DEBUG
    _print_interm( t, _nx, _PMx, "PMx Intermediate", std::cerr );
#endif
    *_PMt = t; // set current time
    if( options.RHSNUMER ){
      if( _RHSBND->solve( _PMVAR, _PMxdot ) != AEBND<T,PMT,PVT>::NORMAL )
        return false;
    }
    else
      _pDAG->eval( _opRHS, _PMRHS, _nx, _pRHS, _PMxdot, _nVAR-_nq, _pVAR, _PMVAR );
    // Whether or not to ignore the remainder
    if( !options.PMNOREM )
      _PMI2vec( _PMenv, _nx, _PMxdot, xdot, true );
    else
      _PMI2vec( _PMenv, _nx, _PMxdot, 0, xdot );
#ifdef MC__ODEBND_BASE_DINEQPM_DEBUG
    _print_interm( t, _nx, _PMxdot, "PMxdot Intermediate", std::cerr );
    {int dum; std::cin >> dum;}
#endif
    return true;  

  case OPT::DINEQ:
    _vec2PMI( x, _PMenv, _nx, _PMx, false );   // set current state polynomial model
    for( unsigned ix=0; ix<_nx; ix++ ) _PMx[ix].center();
#ifdef MC__ODEBND_BASE_DINEQPM_DEBUG
    _print_interm( t, _nx, _PMx, "PMx Intermediate", std::cerr );
#endif
    *_PMt = t; // set current time
    for( unsigned ix=0; ix<_nx; ix++ ){
      T Rxi = _PMVAR[ix].remainder();
      for( unsigned up=0; up<2; up++ ){ // separate lower/upper bounding subproblems
        if( up ) _PMVAR[ix].set( Op<T>::u( Rxi ) );
        else     _PMVAR[ix].set( Op<T>::l( Rxi ) );
        if( options.RHSNUMER ){
          if( _RHSBND->solve( _PMVAR, _PMxdot ) != AEBND<T,PMT,PVT>::NORMAL )
            { _PMVAR[ix].set( Rxi ); return false; }
        }
        else
          _pDAG->eval( _opRHSi[ix], _PMRHS, 1, _pRHS+ix, _PMxdot+ix, _nVAR-_nq,
                       _pVAR, _PMVAR );
        if( up ) _RxUdot[ix] = Op<T>::u( _PMxdot[ix].remainder() );
        else     _RxLdot[ix] = Op<T>::l( _PMxdot[ix].remainder() );
      }
      _PMVAR[ix].set( Rxi );
    }
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
    *_PMt = t; // set current time   
#ifdef MC__ODEBND_BASE_DINEQPM_DEBUG
    _print_interm( t, _nx, _PMx, _Er, "PMx Intermediate", std::cerr );
#endif

    // In this variant a bound on the Jacobian matrix is computed and the
    // linear part is taken as the mid-point of this matrix
    if( !options.ORDMIT ){
      *_It = t; // set current time
      for( unsigned ix=0; ix<_nx; ix++ ){
        _IVAR[ix] = _PMx[ix].bound(); // set current state bounds
        _PMVAR[ix].center().set( T(0.) ); // cancel remainder term
      }
      if( options.RHSNUMER ){
        if( _RHSBND->solve( _PMVAR, _PMxdot ) != AEBND<T,PMT,PVT>::NORMAL )
          return false;
      }
      else
        _pDAG->eval( _opRHS, _PMRHS, _nx, _pRHS, _PMxdot, _nVAR-_nq, _pVAR, _PMVAR );
      _pDAG->eval( _opJAC, _IJAC, _nx*_nx, _pJAC, _Idfdx, _nVAR-_nq, _pVAR, _IVAR );
      _RHS_PM_ELL0( _nx, _PMxdot, _Idfdx, _Ir, _A, _Irdot );
    }

    // In this variant a polynomial model of the Jacobian matrix is computed and the
    // linear part is taken as the mid-point of this matrix
    else if( _PMenv->nord() > _MVXPenv->nord() ){
    //else if( _PMenv->nord() >= _MVXPenv->nord() ){
      *_MVXPt = t; // set current time
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
      if( options.RHSNUMER ){
        if( _RHSBND->solve( _PMVAR, _PMxdot ) != AEBND<T,PMT,PVT>::NORMAL )
          return false;
      }
      else
        _pDAG->eval( _opRHS, _PMRHS, _nx, _pRHS, _PMxdot, _nVAR-_nq, _pVAR, _PMVAR );
      _pDAG->eval( _opJAC, _PMJAC, _nx*_nx, _pJAC, _MVXPdfdx, _nVAR-_nq, _pVAR, _MVXPVAR );
      _RHS_PM_ELL1( _nx, _PMxdot, _MVXPdfdx, _Ir, _A, _Irdot );
    }

    // In this variant a polynomial model in the joint state-parameter and
    // of the same order as the parameter polynomial model is computed
    else{
      *_MVXPt = t; // set current time   
      for( unsigned jx=0; jx<_nx; jx++ ){
        _MVXPd[jx].set( _MVXPenv, _npar+jx, _Ir[jx] );
        _MVXPx[jx].set( _MVXPenv ).set( _PMx[jx].center().set( T(0.) ), true );
      }
      _e2x( _nx, _MVXPd, _MVXPx, false );
#ifdef MC__ODEBND_BASE_DINEQPM_DEBUG
      _print_interm( t, _nx, _MVXPx, "MVXPx Intermediate", std::cerr );
#endif
      if( options.RHSNUMER ){
        if( _RHSBND->solve( _MVXPVAR, _MVXPf ) != AEBND<T,PMT,PVT>::NORMAL )
          return false;
      }
      else
        _pDAG->eval( _opRHS, _PMRHS, _nx, _pRHS, _MVXPf, _nVAR-_nq, _pVAR, _MVXPVAR );
      _RHS_PM_ELL2( _nx, _PMenv, _PMxdot, _MVXPf, _npar, _Ir, _A, _Irdot );
    }
    // Construct the ellipsoidal remainder derivatives
    if( options.QSCALE )
      _RHS_PM_ELL( _nx, _Q, _A, _Irdot, _Qdot, options.QTOL, machprec(), _PMx );
    else
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

} // end namescape mc

#endif

