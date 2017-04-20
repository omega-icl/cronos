// Copyright (C) 2017-... Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

#ifndef MC__ODEBND_EXPAND_HPP
#define MC__ODEBND_EXPAND_HPP

#undef  MC__ODEBND_EXPAND_DEBUG
#undef  MC__ODEBND_EXPAND_NO_TRIVIAL_RESIDUAL

#include "odebnd_base.hpp"
#include "base_expand.hpp"
#include "aebnd.hpp"
#include "odeslv_sundials.hpp"

namespace mc
{
//! @brief C++ class computing enclosures of the reachable set of parametric ODEs using discretization and algebraic equation bounding in AEBND.
////////////////////////////////////////////////////////////////////////
//! mc::ODEBND_EXPAND is a C++ class that computes enclosures of the
//! reachable set of parametric ordinary differential equations
//! (ODEs) using discretization. The available discretization schemes
//! include Taylor expansion, explicit and implicit Runge-Kutta, and
//! Radau orthogonal collocation.
////////////////////////////////////////////////////////////////////////
template <typename T, typename PMT, typename PVT>
class ODEBND_EXPAND:
  public virtual BASE_EXPAND,
  public virtual ODEBND_BASE<T,PMT,PVT>,
  public virtual BASE_DE
{
  typedef BASE_DE::STATUS STATUS;
  typedef AEBND<T,PMT,PVT> t_AEBND;

 protected:
  //using ODEBND_BASE<T,PMT,PVT>::_diam;
  using ODEBND_BASE<T,PMT,PVT>::_hausdorff;
  using ODEBND_BASE<T,PMT,PVT>::_bounds;
  using ODEBND_BASE<T,PMT,PVT>::_print_interm;

  //using ODEBND_BASE<T,PMT,PVT>::_PMenv;
  //using ODEBND_BASE<T,PMT,PVT>::_PMx;

  using ODEBND_BASE<T,PMT,PVT>::_pIC;
  using ODEBND_BASE<T,PMT,PVT>::_pRHS;
  using ODEBND_BASE<T,PMT,PVT>::_pQUAD;

  //using ODEBND_BASE<T,PMT,PVT>::_INI_PM_STA;
  using ODEBND_BASE<T,PMT,PVT>::_IC_SET;
  //using ODEBND_BASE<T,PMT,PVT>::_IC_PM_QUAD;
  //using ODEBND_BASE<T,PMT,PVT>::_IC_PM_STA;
  //using ODEBND_BASE<T,PMT,PVT>::_CC_PM_STA;
  using ODEBND_BASE<T,PMT,PVT>::_RHS_SET;
  using ODEBND_BASE<T,PMT,PVT>::_FCT_I_STA;
  using ODEBND_BASE<T,PMT,PVT>::_FCT_PM_STA;

 protected:
  //! @brief Implicit equation bounder
  t_AEBND _AEBND;

  //! @brief stepsize
  double _h;

  //! @brief position in _vIC
  unsigned _pos_ic;

  //! @brief position in _vRHS
  unsigned _pos_rhs;

  //! @brief position in _vQUAD
  unsigned _pos_quad;
  
  //! @brief position in _vFCT
  unsigned _pos_fct;

  //! @brief pointer to discretized ODE residuals
  std::vector<std::vector<FFVar>> _vRES;

  //! @brief pointer to discretized ODE dependents
  std::vector<FFVar> _vDEP;

  //! @brief pointer to discretized ODE independents
  std::vector<FFVar> _vVAR;

  //! @brief Interval bounds on discretized ODE dependents
  std::vector<T> _IDEP;

  //! @brief Interval bounds on discretized ODE independents
  std::vector<T> _IVAR;

  //! @brief Polynomial model bounds on discretized ODE dependents
  std::vector<PVT> _PMDEP;

  //! @brief Polynomial model bounds on discretized ODE independents
  std::vector<PVT> _PMVAR;

  //! @brief static pointer to local integrator class
  ODESLV_SUNDIALS pODESLV;

public:
  typedef typename ODEBND_BASE<T,PMT,PVT>::Results Results;

  //! @brief Default constructor
  ODEBND_EXPAND
    ();

  //! @brief Virtual destructor
  virtual ~ODEBND_EXPAND
    ();

  //! @brief Integrator options
  struct Options: public BASE_EXPAND::Options
  {
    //! @brief Constructor
    Options():
      BASE_EXPAND::Options(), TORD(4), LBLK(3), DBLK(2), DMAX(1e20), DISPLAY(1), 
      RESRECORD(false), ODESLV(typename ODESLV_SUNDIALS::Options()),
      AEBND(typename AEBND<T,PMT,PVT>::Options())
      { AEBND.DISPLAY = 0; }
    //! @brief Assignment operator
    template <typename OPT> Options& operator=
      ( OPT&options ){
        BASE_SUNDIALS::Options::operator=(options);
        TORD      = options.TORD;
        LBLK      = options.LBLK;
        DBLK      = options.DBLK;
        DMAX      = options.DMAX;
        DISPLAY   = options.DISPLAY;
        RESRECORD = options.RESRECORD;
        ODESLV    = options.ODESLV;
        AEBND     = options.AEBND;
        return *this;
      }
    //! @brief Expansion order (Default: 4)
    unsigned int TORD;
    //! @brief Block length for simultaneous bounding (Default: 3)
    unsigned int LBLK;
    //! @brief Block shift at each step (Default: 2)
    unsigned int DBLK;
    //! @brief Maximum enclosure diameter, \f$D_{\rm max}\f$ (Default: 1e20)
    double DMAX;
    //! @brief Display level (default: 1)
    int DISPLAY;
    //! @brief Whether or not to record results (default: false)
    bool RESRECORD;
    //! @brief Options of real-valued integrator for bounds sampling
    typename ODESLV_SUNDIALS::Options ODESLV;
    //! @brief Options of implicit equation bounder
    typename AEBND<T,PMT,PVT>::Options AEBND;
  } options;

  //! @brief Structure for setting up storing the solver exceptions
  class Exceptions
  {
  public:
    //! @brief Enumeration type for ODEBND_EXPAND exception handling
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
        return "ODEBND_EXPAND::Exceptions  Internal error";
       }
    }
  private:
    TYPE _ierr;
  };

  //! @brief Vector storing interval bound results (upon request only)
  std::vector< Results > results_sta;

  //! @brief Statistics for state bounds integration
  Stats stats_sta;

  //! @brief Set-up set-valued integration of parametric ODEs
  STATUS setup
    ();

  //! @brief Computes interval enclosure of reachable set of parametric ODEs
  STATUS bounds
    ( const T*Ip, T**Ixk, T*If, std::ostream&os=std::cout );

  //! @brief Computes Hausdorff distance between interval enclosure and actual reachable set of parametric ODEs, using parameter sampling
  STATUS hausdorff
    ( const T*Ip, double**Hxk, double*Hf, const unsigned nsamp, std::ostream&os=std::cout );

  //! @brief Computes polynomial model enclosure of reachable set of parametric ODEs
  STATUS bounds
    ( const PVT*PMp, PVT**PMxk, PVT*PMf, std::ostream&os=std::cout );

  //! @brief Computes Hausdorff distance between polynomial model remainder enclosure and actual remainder function range, using parameter sampling
  STATUS hausdorff
    ( const PVT*PMp, double**Hxk, double*Hf, const unsigned nsamp, std::ostream&os=std::cout );

 //! @brief Compute approximate interval enclosure of reachable set of parametric ODEs using parameter sampling
  STATUS bounds
    ( const T*Ip, T**Ixk, T*If, const unsigned nsamp, std::ostream&os=std::cout );

  //! @brief Record results in file <a>bndrec</a>, with accuracy of <a>iprec</a> digits
  void record
    ( std::ofstream&bndrec, const unsigned iprec=5 ) const
    { return ODEBND_BASE<T,PMT,PVT>::_record( bndrec, results_sta, iprec ); }

protected:
  //! @brief Function to initialize state interval bounding
  void _INI_STA
    ();

  //! @brief Function to finalize state bounding
  void _END_STA
    ();

  //! @brief Function to construct IC residuals
  void _AE_SET_ICRES
    ( const FFVar*x, FFVar*res );

  //! @brief Function to construct (trivial) CC residuals
  void _AE_SET_CCRES
    ( const FFVar*x0, const FFVar*x, FFVar*res );

  //! @brief Function to construct ODE residuals using Taylor series expansion
  void _AE_SET_RHSRES
    ( const unsigned ORD, const FFVar*x1, const FFVar*x2, const FFVar*t1,
      const FFVar*t2,  const FFVar&h, FFVar*res, const bool ITS=false );

  //! @brief Function to construct nonlinear equation system from ODE residuals
  void _AE_SET
    ( const bool init );

  //! @brief Function to bound ODE residual solutions using interval arithmetic
  typename t_AEBND::STATUS _AE_SOLVE
    ( const bool init, const T*Ip, const T*Ix0, const T*Ixap );

  //! @brief Function to bound ODE residual solutions using polynomial model arithmetic
  typename t_AEBND::STATUS _AE_SOLVE
    ( const bool init, const PVT*PMp, const PVT*PMx0, const PVT*PMxap );

  //! @brief Private methods to block default compiler methods
  ODEBND_EXPAND(const ODEBND_EXPAND&);
  ODEBND_EXPAND& operator=(const ODEBND_EXPAND&);
};

template <typename T, typename PMT, typename PVT> inline
ODEBND_EXPAND<T,PMT,PVT>::ODEBND_EXPAND
()
: //BASE_DE(), BASE_EXPAND(), ODEBND_BASE<T,PMT,PVT>(),
  _AEBND(), pODESLV()
{}

template <typename T, typename PMT, typename PVT> inline
ODEBND_EXPAND<T,PMT,PVT>::~ODEBND_EXPAND
()
{}

template <typename T, typename PMT, typename PVT> inline void
ODEBND_EXPAND<T,PMT,PVT>::_INI_STA
()
{
  // Reset result record and statistics
  results_sta.clear();
  _init_stats( stats_sta );
}

template <typename T, typename PMT, typename PVT> inline void
ODEBND_EXPAND<T,PMT,PVT>::_END_STA
()
{
  // Get final CPU time
  _final_stats( stats_sta );
}

template <typename T, typename PMT, typename PVT>
inline typename ODEBND_EXPAND<T,PMT,PVT>::STATUS
ODEBND_EXPAND<T,PMT,PVT>::setup
()
{
  const unsigned LBLK = options.LBLK;
  const unsigned TORD = options.TORD;

  _vVAR.resize( _np+_nx+_nq+2+LBLK+1 );  // independents
  for( unsigned i=0; i<_np; i++)
    _vVAR[i] = _pP[i];                          // parameters
  if( _pT ) _vVAR[_np+_nx+_nq] = *_pT;          // initial time
  _vVAR[_np+_nx+_nq+1].set( _pDAG );            // time step
  for( unsigned k=0; _pT && k<LBLK+1; k++)
    _vVAR[_np+_nx+_nq+2+k].set( _pDAG );        // stage times

  _vDEP.resize( (_nx+_nq)*(LBLK+1) );    // discretized states
  for( unsigned i=0; i<(_nx+_nq)*(LBLK+1); i++ )
    _vDEP[i].set( _pDAG );

  _vRES.resize( _nsmax );                    // residuals
  
  // Define AE systems for the discretized ODEs in each stage
  for( _istg=0; _istg<_nsmax; _istg++ ){
    _vRES[_istg].resize( (_nx+_nq)*(LBLK+2) );

    _pos_ic = ( _vIC.size()>=_nsmax? _istg:0 );
    if( !_IC_SET( _pos_ic ) ) return FATAL;
    _AE_SET_ICRES( _vDEP.data(), _vRES[_istg].data() );

    _pos_rhs  = ( _vRHS.size()<=1? 0: _istg );
    _pos_quad = ( _vQUAD.size()<=1? 0: _istg );
    if( !_RHS_SET( _pos_rhs, _pos_quad ) ) return FATAL;
    for( unsigned k=0; k<LBLK; k++ )
      _AE_SET_RHSRES( TORD, _vDEP.data()+k*(_nx+_nq),
        _vDEP.data()+(k+1)*(_nx+_nq), _pT?_vVAR.data()+_np+_nx+_nq+2+k:0,
        _pT?_vVAR.data()+_np+_nx+_nq+2+k+1:0, _vVAR[_np+_nx+_nq+1],
        _vRES[_istg].data()+(k+1)*(_nx+_nq), false );

    _AE_SET_CCRES( _pX, _vDEP.data(), _vRES[_istg].data()+(LBLK+1)*(_nx+_nq) );
  }

#ifdef MC__ODEBND_EXPAND_DEBUG
  std::ofstream oRES( "RES_TS.dot", std::ios_base::out );
  //_pDAG->dot_script( (_nx+_nq)*(LBLK+1), _vRES[0].data(), oRES );
  _pDAG->dot_script( (_nx+_nq)*LBLK, _vRES[0].data()+_nx+_nq, oRES );
  _pDAG->output( _pDAG->subgraph( (_nx+_nq)*LBLK, _vRES[0].data()+_nx+_nq ) );
  oRES.close();
  { int dum; std::cin >> dum; std::cout << "--PAUSED"; }
#endif

  // Size variable bound vectors
  _IVAR.resize( _np+_nx+_nq+2+LBLK+1 );
  _PMVAR.resize( _np+_nx+_nq+2+LBLK+1 );
  _IDEP.resize( (_nx+_nq)*(LBLK+1) );
  _PMDEP.resize( (_nx+_nq)*(LBLK+1) );

  // Reset AE bounder
  _AEBND.set_dag( _pDAG );
  _AEBND.options = options.AEBND;

  return NORMAL;
}

template <typename T, typename PMT, typename PVT> inline void
ODEBND_EXPAND<T,PMT,PVT>::_AE_SET_ICRES
( const FFVar*x, FFVar*res )
{
  for( unsigned i=0; i<_nx+_nq; i++ ){
    res[i] = x[i];
    if( i<_nx ) res[i] -= _pIC[i];
  }
}

template <typename T, typename PMT, typename PVT> inline void
ODEBND_EXPAND<T,PMT,PVT>::_AE_SET_CCRES
( const FFVar*x0, const FFVar*x, FFVar*res )
{
  for( unsigned i=0; i<_nx+_nq; i++ ){
    res[i] = x[i] - x0[i];
  }
}

template <typename T, typename PMT, typename PVT> inline void
ODEBND_EXPAND<T,PMT,PVT>::_AE_SET_RHSRES
( const unsigned TORD, const FFVar*x1, const FFVar*x2, const FFVar*t1,
  const FFVar*t2, const FFVar&h, FFVar*res, const bool ITS )
{
  std::vector<FFVar> fct( _pRHS, _pRHS+_nx );
  if( _nq ) fct.insert( fct.end(), _pQUAD, _pQUAD+_nq );

  const mc::FFVar* tmp1 = _pDAG->compose( _nx+_nq, fct.data(), _pT?1:0, _pT, ITS?t2:t1 );
#ifdef MC__ODEBND_EXPAND_DEBUG
  _pDAG->output( _pDAG->subgraph( _nx+_nq, tmp1 ) );
#endif
  const mc::FFVar* tmp2 = _pDAG->compose( _nx+_nq, tmp1, _nx, _pX, ITS?x2:x1 );
#ifdef MC__ODEBND_EXPAND_DEBUG
  _pDAG->output( _pDAG->subgraph( _nx+_nq, tmp2 ) );
  { int dum; std::cout << "0 TO CONTINUE "; std::cin >> dum; }
#endif
  const mc::FFVar* derx = _pDAG->TAD( TORD, _nx+_nq, tmp2, _nx+_nq, ITS?x2:x1, ITS?t2:t1 );
  for( unsigned k=0; k<TORD; k++ ){
    for( unsigned i=0; i<_nx+_nq; i++ ){
      if( !k ) res[i]  = derx[TORD*(_nx+_nq)+i];
      else     res[i] += derx[(TORD-k)*(_nx+_nq)+i];
      res[i] *= h;
    }
  }
  if( _pT ) delete[] tmp1;
  delete[] tmp2;
  delete[] derx;
  for( unsigned i=0; i<_nx+_nq; i++ )
    res[i] += x1[i] - x2[i];
}

template <typename T, typename PMT, typename PVT> inline void
ODEBND_EXPAND<T,PMT,PVT>::_AE_SET
( const bool init )
{
  const unsigned LBLK = options.LBLK;
  if( init ){ // independents to include state at end of previous stage
    for( unsigned i=0; i<_nx; i++)
      _vVAR[_np+i] = _pX[i];
    _AEBND.set_var( _np+_nx+2+LBLK+1, _vVAR.data() );
    _AEBND.set_dep( (_nx+_nq)*(LBLK+1), _vDEP.data(),
                    _vRES[_istg].data() );
  }
  else{       // independents to include state at initial time
    for( unsigned i=0; i<_nx+_nq; i++)
#ifdef MC__ODEBND_EXPAND_NO_TRIVIAL_RESIDUAL
      _vVAR[_np+i] = _vDEP[i];
    _AEBND.set_var( _np+_nx+2+LBLK+1, _vVAR.data() );
    _AEBND.set_dep( (_nx+_nq)*LBLK, _vDEP.data()+_nx+_nq,
                    _vRES[_istg].data()+_nx+_nq );
#else
      _vVAR[_np+i] = _pX[i];
    _AEBND.set_var( _np+_nx+2+LBLK+1, _vVAR.data() );
    _AEBND.set_dep( (_nx+_nq)*(LBLK+1), _vDEP.data(),
                    _vRES[_istg].data()+_nx+_nq );
#endif
  }
#ifdef MC__ODEBND_EXPAND_DEBUG
  _AEBND.options.DISPLAY = 2;
  _AEBND.options.BLKDEC = false;//true;
#endif
  _AEBND.setup();
}

template <typename T, typename PMT, typename PVT>
inline typename AEBND<T,PMT,PVT>::STATUS
ODEBND_EXPAND<T,PMT,PVT>::_AE_SOLVE
( const bool init, const T*Ip, const T*Ix0, const T*Ixap )
{
  const unsigned LBLK = options.LBLK;
  assert( _nq == 0 ); // Quadrature enclosures not implemented

  for( unsigned i=0; i<_np; i++)
    _IVAR[i] = Ip[i];                     // parameters
  for( unsigned i=0; i<_nx; i++)
    _IVAR[_np+i] = Ix0[i];                // initial state
  if( _pT ) _IVAR[_np+_nx+_nq] = _t;       // initial time
  _IVAR[_np+_nx+_nq+1] = _h;               // time step
  for( unsigned k=0; _pT && k<LBLK+1; k++)
    _IVAR[_np+_nx+_nq+2+k] = _t+k*_h;      // stage times

  for( unsigned k=0, ik=0; k<LBLK+1; k++, ik+=_nx+_nq)
    for( unsigned i=0; i<_nx; i++ )
      _IDEP[ik+i] = Ixap[i];              // state enclosure

  return init? _AEBND.solve( _IVAR.data(), _IDEP.data(), _IDEP.data() ):
#ifdef MC__ODEBND_EXPAND_NO_TRIVIAL_RESIDUAL
               _AEBND.solve( _IVAR.data(), _IDEP.data()+_nx+_nq, _IDEP.data()+_nx+_nq );
#else
               _AEBND.solve( _IVAR.data(), _IDEP.data(), _IDEP.data() );
#endif
}

//! @fn template <typename T, typename PMT, typename PVT> inline typename ODEBND_EXPAND<T,PMT,PVT>::STATUS ODEBND_EXPAND<T,PMT,PVT>::bounds(
//! const T*Ip, T**Ixk, T*If, std::ostream&os=std::cout )
//!
//! This function computes an enclosure of the reachable set of the parametric ODEs
//! using propagation of polynomial models with convex remainders (intervals, ellipsoids):
//!   - <a>Ip</a>   [input]  interval bound on parameter set
//!   - <a>Ixk</a>  [output] interval bound on state varibles at stage times
//!   - <a>If</a>   [output] interval bound on state-dependent functions
//!   - <a>os</a>   [input]  output stream [default: std::cout]
//! .
//! The return value is the status.
template <typename T, typename PMT, typename PVT>
inline typename ODEBND_EXPAND<T,PMT,PVT>::STATUS
ODEBND_EXPAND<T,PMT,PVT>::bounds
( const T*Ip, T**Ixk, T*If, std::ostream&os )
{
  // Check arguments
  if( !Ixk || !Ip || (_nf && !If) ) return FATAL;

  const unsigned LBLK = options.LBLK;
  const unsigned DBLK = options.DBLK;
  const double H0 = options.H0, TOL = 1e-8;
  _INI_STA();

  try{
    // Integrate ODEs through each stage
    for( _istg=0; _istg<_nsmax; _istg++ ){
      _t = _dT[_istg];
      _h = ( _t + LBLK*options.H0 > _dT[_istg+1]-TOL? (_dT[_istg+1]-_t)/LBLK: H0 ); 

      // First step of stage w/ initial conditions
      _AE_SET( true );
      auto flag = _AE_SOLVE( true, Ip, Ixk[_istg], Ixk[_istg+1] );
      if( flag != t_AEBND::NORMAL ){ _END_STA(); return FAILURE; }

      // Initial state bounds
      if( !_istg ){
        for( unsigned ix=0; ix<_nx; ix++ )
          Ixk[_istg][ix] = _IDEP[ix];
        if( options.DISPLAY >= 1 )
          _print_interm( _t, _nx, Ixk[0], "x", os );
        for( unsigned k=0; options.RESRECORD && k<=DBLK; k++ )
          results_sta.push_back( Results( _t+k*_h, _nx, _IDEP.data()+(_nx+_nq)*k ) );
      }

      // Next steps w/o initial conditions
      bool reinit = true;
      while( _t+LBLK*_h < _dT[_istg+1]-TOL ){
        _t += DBLK*_h;
        _h = ( _t + LBLK*options.H0 > _dT[_istg+1]-TOL? (_dT[_istg+1]-_t)/LBLK: H0 ); 

        if( reinit ){ _AE_SET( false ); reinit = false; }
        auto flag = _AE_SOLVE( false, Ip, _IDEP.data()+(_nx+_nq)*DBLK, Ixk[_istg+1] );
        if( flag != t_AEBND::NORMAL ){ _END_STA(); return FAILURE; }

        for( unsigned k=1; options.RESRECORD && k<=DBLK; k++ )
          results_sta.push_back( Results( _t+k*_h, _nx, _IDEP.data()+(_nx+_nq)*k ) );
      }
      _t += LBLK*_h;

      // End-stage state bounds
      for( unsigned ix=0; ix<_nx; ix++ )
        Ixk[_istg+1][ix] = _IDEP[(_nx+_nq)*LBLK+ix];
      if( options.DISPLAY >= 1 )
        _print_interm( _t, _nx, Ixk[_istg+1], "x", os );

      // Add intermediate function terms
      _pos_fct = ( _vFCT.size()>=_nsmax? _istg:0 );
      if( _nf
       && (_vFCT.size()>=_nsmax || _istg==_nsmax-1)
       && !_FCT_I_STA( _pos_fct, _t, If ) )
        { _END_STA(); return FATAL; }
    }

    // Bounds on final quadratures and functions
    if( options.DISPLAY >= 1 ) _print_interm( _nf, If, "f", os );
  }


  catch( typename t_AEBND::Exceptions &expt ){
    if( options.DISPLAY >= 1 ){
      os << expt.what() << std::endl;
      _END_STA();
      os << " ABORT TIME  " << std::scientific << std::left
                            << std::setprecision(5) << _t << std::endl;
      _print_stats( stats_sta, os );
    }
    return FAILURE;
  }

  catch(...){
    _END_STA();
    if( options.DISPLAY >= 1 ){
      os << " ABORT TIME  " << std::scientific << std::left
                            << std::setprecision(5) << _t << std::endl;
      _print_stats( stats_sta, os );
    }
    return FAILURE;
  }

  _END_STA();
  if( options.DISPLAY >= 1 ) _print_stats( stats_sta, os );
  return NORMAL;
}

template <typename T, typename PMT, typename PVT>
inline typename AEBND<T,PMT,PVT>::STATUS
ODEBND_EXPAND<T,PMT,PVT>::_AE_SOLVE
( const bool init, const PVT*PMp, const PVT*PMx0, const PVT*PMxap )
{
  const unsigned LBLK = options.LBLK;
  assert( _nq == 0 ); // Quadrature enclosures not implemented

  for( unsigned i=0; i<_np; i++)
    _PMVAR[i] = PMp[i];                     // parameters
  for( unsigned i=0; i<_nx; i++)
    _PMVAR[_np+i] = PMx0[i];                // initial state
  if( _pT ) _PMVAR[_np+_nx+_nq] = _t;       // initial time
  _PMVAR[_np+_nx+_nq+1] = _h;               // time step
  for( unsigned k=0; _pT && k<LBLK+1; k++)
    _PMVAR[_np+_nx+_nq+2+k] = _t+k*_h;      // stage times

  for( unsigned k=0, ik=0; k<LBLK+1; k++, ik+=_nx+_nq)
    for( unsigned i=0; i<_nx; i++ )
      _PMDEP[ik+i] = PMxap[i];              // state enclosure

  return init? _AEBND.solve( _PMVAR.data(), _PMDEP.data(), _PMDEP.data() ):
#ifdef MC__ODEBND_EXPAND_NO_TRIVIAL_RESIDUAL
               _AEBND.solve( _PMVAR.data(), _PMDEP.data()+_nx+_nq, _PMDEP.data()+_nx+_nq );
#else
               _AEBND.solve( _PMVAR.data(), _PMDEP.data(), _PMDEP.data() );
#endif
}

//! @fn template <typename T, typename PMT, typename PVT> inline typename ODEBND_EXPAND<T,PMT,PVT>::STATUS ODEBND_EXPAND<T,PMT,PVT>::bounds(
//! const PVT*PMp, PVT**PMxk, PVT*PMf, std::ostream&os=std::cout )
//!
//! This function computes an enclosure of the reachable set of the parametric ODEs
//! using propagation of polynomial models with convex remainders (intervals, ellipsoids):
//!   - <a>PMp</a>  [input]  polynomial model of parameter set
//!   - <a>PMxk</a> [output] polynomial model of state varibles at stage times
//!   - <a>PMf</a>  [output] polynomial model of state-dependent functions
//!   - <a>os</a>   [input]  output stream [default: std::cout]
//! .
//! The return value is the status.
template <typename T, typename PMT, typename PVT>
inline typename ODEBND_EXPAND<T,PMT,PVT>::STATUS
ODEBND_EXPAND<T,PMT,PVT>::bounds
( const PVT*PMp, PVT**PMxk, PVT*PMf, std::ostream&os )
{
  // Check arguments
  if( !PMxk || !PMp || (_nf && !PMf) ) return FATAL;

  const unsigned LBLK = options.LBLK;
  const unsigned DBLK = options.DBLK;
  const double H0 = options.H0, TOL = 1e-8;
  _INI_STA();

  try{
    // Integrate ODEs through each stage
    for( _istg=0; _istg<_nsmax; _istg++ ){
      _t = _dT[_istg];
      _h = ( _t + LBLK*options.H0 > _dT[_istg+1]-TOL? (_dT[_istg+1]-_t)/LBLK: H0 ); 

      // First step of stage w/ initial conditions
      _AE_SET( true );
      auto flag = _AE_SOLVE( true, PMp, PMxk[_istg], PMxk[_istg+1] );
      if( flag != t_AEBND::NORMAL ){ _END_STA(); return FAILURE; }

      // Initial state bounds
      if( !_istg ){
        for( unsigned ix=0; ix<_nx; ix++ )
          PMxk[_istg][ix] = _PMDEP[ix];
        if( options.DISPLAY >= 1 )
          _print_interm( _t, _nx, PMxk[0], "x", os );
        for( unsigned k=0; options.RESRECORD && k<=DBLK; k++ )
          results_sta.push_back( Results( _t+k*_h, _nx, _PMDEP.data()+(_nx+_nq)*k ) );
      }

      // Next steps w/o initial conditions
      bool reinit = true;
      while( _t+LBLK*_h < _dT[_istg+1]-TOL ){
        _t += DBLK*_h;
        _h = ( _t + LBLK*options.H0 > _dT[_istg+1]-TOL? (_dT[_istg+1]-_t)/LBLK: H0 ); 

        if( reinit ){ _AE_SET( false ); reinit = false; }
        auto flag = _AE_SOLVE( false, PMp, _PMDEP.data()+(_nx+_nq)*DBLK, PMxk[_istg+1] );
        if( flag != t_AEBND::NORMAL ){ _END_STA(); return FAILURE; }

        for( unsigned k=1; options.RESRECORD && k<=DBLK; k++ )
          results_sta.push_back( Results( _t+k*_h, _nx, _PMDEP.data()+(_nx+_nq)*k ) );
      }
      _t += LBLK*_h;

      // End-stage state bounds
      for( unsigned ix=0; ix<_nx; ix++ )
        PMxk[_istg+1][ix] = _PMDEP[(_nx+_nq)*LBLK+ix];
      if( options.DISPLAY >= 1 )
        _print_interm( _t, _nx, PMxk[_istg+1], "x", os );

      // Add intermediate function terms
      _pos_fct = ( _vFCT.size()>=_nsmax? _istg:0 );
      if( _nf
       && (_vFCT.size()>=_nsmax || _istg==_nsmax-1)
       && !_FCT_PM_STA( _pos_fct, _t, PMf ) )
        { _END_STA(); return FATAL; }
    }

    // Bounds on final quadratures and functions
    if( options.DISPLAY >= 1 ) _print_interm( _nf, PMf, "f", os );
  }

  catch( typename t_AEBND::Exceptions &expt ){
    if( options.DISPLAY >= 1 ){
      os << expt.what() << std::endl;
      _END_STA();
      os << " ABORT TIME  " << std::scientific << std::left
                            << std::setprecision(5) << _t << std::endl;
      _print_stats( stats_sta, os );
    }
    return FAILURE;
  }

  catch(...){
    _END_STA();
    if( options.DISPLAY >= 1 ){
      os << " ABORT TIME  " << std::scientific << std::left
                            << std::setprecision(5) << _t << std::endl;
      _print_stats( stats_sta, os );
    }
    return FAILURE;
  }

  _END_STA();
  if( options.DISPLAY >= 1 ) _print_stats( stats_sta, os );
  return NORMAL;
}

//! @fn template <typename T, typename PMT, typename PVT> inline typename ODEBND_EXPAND<T,PMT,PVT>::STATUS ODEBND_EXPAND<T,PMT,PVT>::hausdorff(
//! const T*Ip, double**Hxk, double**Hf, const unsigned nsamp, std::ostream&os )
//!
//! This function computes the Hausdorff distance between the interval enclosure
//! and the exact reachable set projected onto each variable:
//!   - <a>Ip</a>    [input]  interval parameter set
//!   - <a>Hxk</a>   [output] Hausdorff distance between the interval enclosure
//!                           and the exact reachable set projected onto each
//!                           variable for states at stage times
//!   - <a>Hf</a>    [output] Hausdorff distance between the interval enclosure
//!                           and the exact reachable set projected onto each
//!                           variable for state-dependent functions
//!   - <a>nsamp</a> [input]  number of samples for each parameter
//!   - <a>os</a>    [input]  output stream
//! .
//! The return value is the status.
template <typename T, typename PMT, typename PVT>
inline typename ODEBND_EXPAND<T,PMT,PVT>::STATUS
ODEBND_EXPAND<T,PMT,PVT>::hausdorff
( const T*Ip, double**Hxk, double*Hf, const unsigned nsamp, std::ostream&os )
{
  pODESLV.set( *this );
  pODESLV.options = options.ODESLV;
  return _hausdorff( Ip, Hxk, Hf, *this, pODESLV, nsamp, os )? NORMAL: FAILURE;
}

//! @fn template <typename T, typename PMT, typename PVT> inline typename ODEBND_EXPAND<T,PMT,PVT>::STATUS ODEBND_EXPAND<T,PMT,PVT>::hausdorff(
//! const PVT*PMp, double**Hxk, const unsigned nsamp, std::ostream&os=std::cout )
//!
//! This function computes the Hausdorff distance between the polynomial model
//! remainder and the actual (sampled) range of the remainder function
//! in projection onto each variable and for each stage time
//! remainder and the actual range of the remainder function:
//!   - <a>ns</a>    [input]  number of time stages
//!   - <a>tk</a>    [input]  stage times, including the initial time
//!   - <a>PMp</a>   [input]  polynomial model of parameter set
//!   - <a>Hxk</a>   [output] Hausdorff distance between the polynomial model
//!                           remainder and the actual (sampled) range of the
//!                           remainder term for states at stage times
//!   - <a>Hf</a>    [output] Hausdorff distance between the polynomial model
//!                           remainder and the actual (sampled) range of the
//!                           remainder term for state-dependent functions
//!   - <a>nsamp</a> [input]  number of samples for each parameter
//!   - <a>os</a>    [input]  output stream [default: std::cout]

template <typename T, typename PMT, typename PVT>
inline typename ODEBND_EXPAND<T,PMT,PVT>::STATUS
ODEBND_EXPAND<T,PMT,PVT>::hausdorff
( const PVT*PMp, double**Hxk, double*Hf, const unsigned nsamp, std::ostream&os )
{
  pODESLV.set( *this );
  pODESLV.options = options.ODESLV;
  return _hausdorff( PMp, Hxk, Hf, *this, pODESLV, nsamp, os )? NORMAL: FAILURE;
}

//! @fn template <typename T, typename PMT, typename PVT> inline typename ODEBND_EXPAND<T,PMT,PVT>::STATUS ODEBND_EXPAND<T,PMT,PVT>::bounds(
//! const T*Ip, T**Ixk, T*q, T*f, const unsigned nsamp, std::ostream&os=std::cout )
//!
//! This function computes projections of an inner-approximation enclosure of
//! the reachable set of the parametric ODEs using sampling and continuous-time
//! integration:
//!   - <a>Ip</a>    [input] interval parameter set
//!   - <a>Ixk</a>   [output] approximate interval state enclosures at stage times
//!   - <a>If</a>    [output] approximate function enclosures
//!   - <a>nsamp</a> [input] number of samples for each parameter
//!   - <a>os</a>    [input] output stream
//! .
//! The return value is the status.
template <typename T, typename PMT, typename PVT>
inline typename ODEBND_EXPAND<T,PMT,PVT>::STATUS
ODEBND_EXPAND<T,PMT,PVT>::bounds
( const T*Ip, T**Ixk, T*If, const unsigned nsamp, std::ostream&os )
{
  // Sample inner approximation
  STATUS flag = NORMAL;
  pODESLV.set( *this );
  pODESLV.options = options.ODESLV;
  if( !_bounds( Ip, Ixk, If, pODESLV, nsamp, os ) )
    flag = FAILURE;

  // Display results
  if( options.DISPLAY >= 1 ){
    for( unsigned is=0; Ixk && is<=_nsmax; is++ )
      _print_interm( _dT[is], _nx, Ixk[is], "x", os );
    if( If ) _print_interm( _nf, If, "f", os );
  }

  // Record intermediate results
  results_sta.clear();
  if( options.RESRECORD )
    for( unsigned is=0; Ixk && is<=_nsmax; is++ )
      results_sta.push_back( Results( _dT[is], _nx, Ixk[is] ) );
  
  return flag;
}

} // end namescape mc

#endif

