// Copyright (C) Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the EPL-2.0 with GPL-2.0-or-later as a Secondary License; see the LICENSE file.
// SPDX-License-Identifier: EPL-2.0 OR GPL-2.0-or-later

#ifndef CRONOS__DAESLV_IDAS_HPP
#define CRONOS__DAESLV_IDAS_HPP

#include <map>
#include <set>
#include <string>
#include <tuple>
#include <vector>

#include "daeslv_base.hpp"
#include "base_idas.hpp"

namespace mc
{
//! @brief C++ class integrating a model's DAE form with SUNDIALS IDAS
////////////////////////////////////////////////////////////////////////
//! mc::DAESLV_IDAS is the IDAS counterpart of mc::ODESLV_CVODES.  It takes the general mass-matrix form
//! that DAESLV_BASE extracts,
//!
//!     sum over _mMASS[0][k] of ( coefficient * dy[index] )  +  _mRES[0][k]  =  0 ,
//!
//! and hands it to IDAS as the residual, with the id vector marking which unknowns are differential.
//!
//! WHAT v0 DOES NOT DO YET, each deliberate rather than overlooked:
//!   * ONE SET OF ROWS -- a model with TIME-VARYING INPUTS is refused by setup(), since its rows differ
//!     stage by stage; without them the _ns stages hold identical rows and one set covers the horizon;
//!   * DIFFERENCE-QUOTIENT Jacobian (IDAS supplies it when no IDALsJacFn is set).  The analytic one is a
//!     later step; the residual is symbolic, so it is available, just not wired;
//!   * no quadratures, no point outputs, no sensitivities.  The eight BASE_IDAS sensitivity callbacks are
//!     present and REFUSE, so phases 5/6 fill them in rather than finding a silent wrong answer;
//!   * consistent initialisation is IDACalcIC (option A of the plan).  Option C -- a one-element OCFESLV
//!     solve -- is the fallback for multi-rooted ICs and is not wired here.
//!
//! The residual is evaluated in the MODEL's dag rather than a solver-local copy, so nothing is inserted
//! or renamed; the coefficients and remainders come out of one subgraph in a single eval, and the
//! derivative values are folded in afterwards.  That keeps v0 small; a solver-local dag is what
//! ODESLV_CVODES does, and is the thing to revisit when threads or a second concurrent solve appear.
////////////////////////////////////////////////////////////////////////
class DAESLV_IDAS
: public virtual DAESLV_BASE
, public virtual BASE_IDAS
{
public:
  DAESLV_IDAS()
    : _nres(0), _ncoef(0)
    {}
  virtual ~DAESLV_IDAS()
    { unregistration(); delete _dag; }

  //! @brief The INTEGRATOR's options.  Both FFModel and BASE_IDAS declare `options`, so an unqualified
  //! `solver.options` would be ambiguous; this names the IDAS one, and the model's stay reachable as
  //! `solver.FFModel::options`.  (ODESLV_CVODES instead declares its own Options member that hides both --
  //! worth revisiting here once DAESLV_IDAS needs options of its own.)
  using BASE_IDAS::options;

  //! @brief Statistics of the last solve
  Stats stats;

  //! @brief Extract the model as a DAE and build the evaluation plan.  False with err() set.
  bool setup
    ();

  //! @brief How many of the np() parameters are per-stage INPUT LEVELS rather than model parameters.  They
  //! are passed in solve()'s parameter vector like any other.  Non-zero exactly when the model has
  //! time-varying inputs, which is also when the rows differ stage by stage and the integrator restarts.
  size_t n_input_levels
    ()
    const
    { return _vInpLvl.size(); }

  //! @brief Integrate from the initial time to each time in @p tout, with parameter values @p p.
  //! The state at the last time is left in val_state().  False with err() set.
  bool solve
    ( std::vector<double> const& p, std::vector<double> const& tout );

  //! @brief Nonzeros of the Jacobian IDAS factorizes -- the UNION of the dF/dy and dF/dy' patterns.
  size_t nnz_jac
    ()
    const
    { return _nnzjac; }

  //! @brief The state vector at the end of the last solve.
  std::vector<double> const& val_state
    ()
    const
    { return _dX; }

  //! @brief The state vector right after consistent initialisation, before any step.
  std::vector<double> const& val_state_init
    ()
    const
    { return _dX0; }

  //! @brief The derivative vector right after consistent initialisation.
  std::vector<double> const& val_deriv_init
    ()
    const
    { return _dXP0; }

protected:

  size_t _nres;                              //!< rows of the residual (== _nx)
  size_t _ncoef;                             //!< mass-matrix entries actually present
  size_t _offJy, _offJyp;                    //!< where the dF/dy and dF/dy' entries start in _dEVAL
  //! SOLVER-LOCAL dag, as ODESLV_BASE::_SETUP builds one: the model's description is INSERTED here and
  //! everything downstream -- residual, Jacobian blocks, evaluation -- happens in this graph.  Note the
  //! model's dag is never mutated: the mass coefficients and remainders are inserted FIRST and the
  //! derivative variables created HERE, so dy/dt exists only in the solver's own graph.
  FFGraph* _dag = nullptr;
  std::vector<FFVar> _pT, _pX, _pP, _pYP;    //!< the inserted variables, and dy/dt, all local
  struct StagePlan {
    std::vector<FFVar> F;                    //!< the residual rows of this stage
    std::vector<FFVar> eval;                 //!< [ F, dF/dy entries, dF/dy' entries ]
    std::vector<int>   ixJy, ixJyp;          //!< where each entry lands in the CSC value array
    std::vector<sunindextype> colptr, rowidx;
    size_t nJy = 0, nJyp = 0, nnz = 0, offJy = 0, offJyp = 0;
  };
  std::vector<StagePlan> _plan;              //!< one per stage (one entry when the rows do not change)
  std::vector<FFVar> _vYP;                   //!< unused; kept so older call sites still compile
  size_t _nJy = 0, _nJyp = 0;                //!< entry counts of the two blocks, for the active stage
  std::vector<int> _ixJy, _ixJyp;            //!< where each block entry lands in the CSC value array
  std::vector<sunindextype> _colptr, _rowidx;//!< the union pattern, CSC
  size_t _nnzjac = 0;                        //!< nonzeros of dF/dy + c_j dF/dy' for the CURRENT stage
  size_t _nnzjacmax = 0;                     //!< the worst stage's count, which the matrix is sized on
  SUNMatrix _sun_mat = nullptr;              //!< kept so a restart can reinitialize KLU on it
  SUNLinearSolver _sun_ls = nullptr;
  std::vector<FFVar>  _vInpLvl;              //!< the per-stage input LEVEL variables, in declaration order
  std::vector< std::pair<size_t,size_t> > _coefAt;  //!< (row, state) of each coefficient, in eval order
  std::vector<FFVar>  _vEVAL;                //!< [ every coefficient, then every remainder ]
  std::vector<FFVar>  _vVAR;                 //!< [ t, states..., parameters... ]
  std::vector<double> _dVAR;                 //!< values of _vVAR
  std::vector<double> _dEVAL;                //!< values of _vEVAL
  std::vector<double> _dWRK;                 //!< eval work array
  std::vector<FFVar>  _vIC;                  //!< initial values, as expressions
  std::vector<double> _dX, _dX0, _dXP0;      //!< state now, and (state, derivative) after IDACalcIC

  //! @brief The DAE residual: F = M(t,y,p) * yp + R(t,y,p).
  int IDRES__
    ( sunrealtype tt, N_Vector yy, N_Vector yp, N_Vector rr, void* user_data )
    override;

  //! @brief v0 has no quadratures.
  int IDQUAD__
    ( sunrealtype, N_Vector, N_Vector, N_Vector, void* )
    override
    { return -1; }

  //! @brief v0 uses the IDAS difference-quotient Jacobian; this is never installed.
  //! @brief The analytic Jacobian  dF/dy + c_j dF/dy' .  With F = M(t,y,p) y' + R(t,y,p),
  //!     dF_k/dy_j  =  sum_i ( dM_ki/dy_j ) y'_i  +  dR_k/dy_j        and   dF_k/dy'_i = M_ki ,
  //! so every piece comes from the same symbolic rows the residual uses -- nothing is re-derived.
  int IDJAC__
    ( sunrealtype t, sunrealtype c_j, N_Vector yy, N_Vector yp, N_Vector rr, SUNMatrix Jac,
      void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3 )
    override;

  //! @brief Adjoint sensitivity is phase 6; refusing beats a silent wrong gradient.
  int IDRESB__
    ( sunrealtype, N_Vector, N_Vector, N_Vector, N_Vector, N_Vector, void* )
    override
    { return -1; }
  int IDQUADB__
    ( sunrealtype, N_Vector, N_Vector, N_Vector, N_Vector, N_Vector, void* )
    override
    { return -1; }
  //! @brief Adjoint Jacobian; phase 6, refuses for now.
  int IDJACB__
    ( sunrealtype, sunrealtype, N_Vector, N_Vector, N_Vector, N_Vector, N_Vector, SUNMatrix, void*,
      N_Vector, N_Vector, N_Vector )
    override
    { return -1; }

  //! @brief Forward sensitivity is phase 5.
  int IDRESF__
    ( int, sunrealtype, N_Vector, N_Vector, N_Vector, N_Vector*, N_Vector*, N_Vector*, void*, N_Vector,
      N_Vector, N_Vector )
    override
    { return -1; }
  int IDQUADF__
    ( int, sunrealtype, N_Vector, N_Vector, N_Vector*, N_Vector*, N_Vector, N_Vector*, void*, N_Vector,
      N_Vector, N_Vector )
    override
    { return -1; }

  //! @brief Evaluate the rows and their y-derivatives at (t, y, p) into _dEVAL.
  void _eval_rows
    ( double const t, double const* y, double const* yp );

  //! @brief Build the evaluation plan for stage @p s: the residual as ONE expression in (t, y, y', p), and
  //! its two sparse Jacobian blocks by SFAD, mirroring ODESLV_BASE's _pJAC = _dag->SFAD( nx, RHS, nx, X ).
  bool _build_rows
    ( size_t const s );

  //! @brief Build the solver-local dag and every stage's residual + Jacobian plan.
  bool _build_plans
    ();

  //! @brief Mass-matrix entries of stage @p s.
  size_t _ncoefOf
    ( size_t const s )
    const
    { size_t n = 0; for( auto const& row : _mMASS[s] ) n += row.size(); return n; }
};

//! @brief C trampoline for the residual: IDAS calls this, which reaches the registered solver's IDRES__.
//! Named as ODESLV_CVODES names its own (MC_CVRHS__, MC_CVJAC__).
inline
int
MC_IDRES__
( sunrealtype tt, N_Vector yy, N_Vector yp, N_Vector rr, void* user_data )
{
  return ( BASE_IDAS::PTR_BASE_IDAS->*PTR_IDRES )( tt, yy, yp, rr, user_data );
}

//! @brief C trampoline for the analytic Jacobian.
inline
int
MC_IDJAC__
( sunrealtype tt, sunrealtype c_j, N_Vector yy, N_Vector yp, N_Vector rr, SUNMatrix Jac, void* user_data,
  N_Vector tmp1, N_Vector tmp2, N_Vector tmp3 )
{
  return ( BASE_IDAS::PTR_BASE_IDAS->*PTR_IDJAC )( tt, c_j, yy, yp, rr, Jac, user_data, tmp1, tmp2, tmp3 );
}

inline
bool
DAESLV_IDAS::_build_plans
()
{
  delete _dag; _dag = new FFGraph;
  _pT.assign( 1, FFVar() );   _dag->insert( FFModel::_dag, 1,   _mT.data(), _pT.data() );
  _pX.assign( _nx, FFVar() ); _dag->insert( FFModel::_dag, _nx, _mX.data(), _pX.data() );
  _pP.assign( _np, FFVar() );
  if( _np ) _dag->insert( FFModel::_dag, _np, _mP.data(), _pP.data() );
  _pYP.clear();
  for( size_t i = 0; i < _nx; ++i ) _pYP.push_back( _dag->add_var( "d" + _mX[i].name() + "/dt" ) );

  _vVAR.clear();
  _vVAR.push_back( _pT[0] );
  for( auto const& v : _pX  ) _vVAR.push_back( v );
  for( auto const& v : _pYP ) _vVAR.push_back( v );
  for( auto const& v : _pP  ) _vVAR.push_back( v );
  _dVAR.assign( _vVAR.size(), 0. );

  _plan.assign( _mMASS.size(), StagePlan() );
  _nnzjacmax = 0;
  for( size_t sg = 0; sg < _mMASS.size(); ++sg ){
    StagePlan& P = _plan[sg];
    // insert this stage's coefficients and remainders, then build F HERE
    std::vector<FFVar> M( _ncoefOf( sg ) ), R( _mRES[sg].size() );
    std::vector<FFVar> Msrc;
    std::vector< std::pair<size_t,size_t> > at;
    for( size_t k = 0; k < _mMASS[sg].size(); ++k )
      for( auto const& [ ix, co ] : _mMASS[sg][k] ){ at.emplace_back( k, ix ); Msrc.push_back( co ); }
    if( !Msrc.empty() ) _dag->insert( FFModel::_dag, Msrc.size(), Msrc.data(), M.data() );
    _dag->insert( FFModel::_dag, _mRES[sg].size(), const_cast<FFVar*>( _mRES[sg].data() ), R.data() );
    P.F = R;
    for( size_t c = 0; c < at.size(); ++c ) P.F[ at[c].first ] += M[c] * _pYP[ at[c].second ];

    std::tuple<unsigned,unsigned*,unsigned*,FFVar*> Jy, Jyp;
    try{
      Jy  = _dag->SFAD( P.F.size(), P.F.data(), _nx, _pX.data()  );
      Jyp = _dag->SFAD( P.F.size(), P.F.data(), _nx, _pYP.data() );
    }
    catch( ... ){ _extractError = "SFAD failed on the DAE residual"; return false; }

    std::set< std::pair<size_t,size_t> > pat;                       // (col, row)
    for( unsigned e = 0; e < std::get<0>(Jy);  ++e ) pat.emplace( std::get<2>(Jy)[e],  std::get<1>(Jy)[e]  );
    for( unsigned e = 0; e < std::get<0>(Jyp); ++e ) pat.emplace( std::get<2>(Jyp)[e], std::get<1>(Jyp)[e] );
    P.nnz = pat.size();
    P.colptr.assign( _nx + 1, 0 );
    P.rowidx.clear();
    std::map< std::pair<size_t,size_t>, int > where;
    for( auto const& [ c, r ] : pat ){ where[{c,r}] = (int)P.rowidx.size(); P.rowidx.push_back( (sunindextype)r ); }
    for( auto const& [ c, r ] : pat ) ++P.colptr[ c + 1 ];
    for( size_t c = 0; c < _nx; ++c ) P.colptr[c+1] += P.colptr[c];
    P.nJy = std::get<0>(Jy);   P.ixJy.resize( P.nJy );
    for( unsigned e = 0; e < P.nJy;  ++e ) P.ixJy[e]  = where[{ std::get<2>(Jy)[e],  std::get<1>(Jy)[e]  }];
    P.nJyp = std::get<0>(Jyp); P.ixJyp.resize( P.nJyp );
    for( unsigned e = 0; e < P.nJyp; ++e ) P.ixJyp[e] = where[{ std::get<2>(Jyp)[e], std::get<1>(Jyp)[e] }];

    P.eval = P.F;
    P.offJy  = P.eval.size();  P.eval.insert( P.eval.end(), std::get<3>(Jy),  std::get<3>(Jy)  + P.nJy  );
    P.offJyp = P.eval.size();  P.eval.insert( P.eval.end(), std::get<3>(Jyp), std::get<3>(Jyp) + P.nJyp );
    for( auto* t : { &Jy, &Jyp } ){ delete[] std::get<1>(*t); delete[] std::get<2>(*t); delete[] std::get<3>(*t); }
    if( P.nnz > _nnzjacmax ) _nnzjacmax = P.nnz;
  }
  return _build_rows( 0 );
}

inline
bool
DAESLV_IDAS::_build_rows
( size_t const s )
{
  size_t const ss = ( s < _plan.size() )? s: _plan.size() - 1;
  StagePlan const& P = _plan[ss];
  _vEVAL = P.eval;
  _offJy = P.offJy;  _offJyp = P.offJyp;
  _nJy = P.nJy;      _nJyp = P.nJyp;
  _ixJy = P.ixJy;    _ixJyp = P.ixJyp;
  _colptr = P.colptr; _rowidx = P.rowidx;
  _nnzjac = P.nnz;
  _dEVAL.assign( _vEVAL.size(), 0. );
  return true;
}

inline
int
DAESLV_IDAS::IDJAC__
( sunrealtype t, sunrealtype c_j, N_Vector yy, N_Vector yp, N_Vector, SUNMatrix Jac, void*, N_Vector,
  N_Vector, N_Vector )
{
  _eval_rows( t, N_VGetArrayPointer( yy ), N_VGetArrayPointer( yp ) );
  std::vector<double> v( _nnzjac, 0. );
  for( size_t e = 0; e < _nJy;  ++e ) v[ _ixJy[e]  ] += _dEVAL[ _offJy  + e ];
  for( size_t e = 0; e < _nJyp; ++e ) v[ _ixJyp[e] ] += c_j * _dEVAL[ _offJyp + e ];
  if( SUNMatGetID( Jac ) == SUNMATRIX_SPARSE ){
    sunindextype* cp = SUNSparseMatrix_IndexPointers( Jac );
    sunindextype* ri = SUNSparseMatrix_IndexValues( Jac );
    sunrealtype*  d  = SUNSparseMatrix_Data( Jac );
    for( size_t c = 0; c <= _nx; ++c ) cp[c] = _colptr[c];
    for( size_t e = 0; e < _nnzjac; ++e ){ ri[e] = _rowidx[e]; d[e] = v[e]; }
  }
  else{
    SUNMatZero( Jac );
    for( size_t c = 0; c < _nx; ++c )
      for( sunindextype e = _colptr[c]; e < _colptr[c+1]; ++e )
        SM_ELEMENT_D( Jac, _rowidx[e], c ) = v[e];
  }
  ++stats.numJAC;
  return 0;
}

inline
bool
DAESLV_IDAS::setup
()
{
  if( !FFModel::setup() ){ _extractError = "FFModel::setup() failed"; return false; }
  if( !_extract_dae() ) return false;
  _nres = n_residual();

  if( _nc ){
    _extractError = "this model declares constants; DAESLV_IDAS v1 does not insert them into its local dag yet";
    return false;
  }
  if( !_build_plans() ) return false;
  _vInpLvl.clear();
  for( auto const& [ name, levels ] : _mInpLevels )
    for( auto const& lv : levels ) _vInpLvl.push_back( lv );
  _vIC.assign( _mIC[0].size(), FFVar() );
  _dag->insert( FFModel::_dag, _mIC[0].size(), const_cast<FFVar*>( _mIC[0].data() ), _vIC.data() );
  _dX.assign( _nx, 0. );
  _dX0.assign( _nx, 0. );
  _dXP0.assign( _nx, 0. );
  registration();
  return true;
}

inline
void
DAESLV_IDAS::_eval_rows
( double const t, double const* y, double const* yp )
{
  _dVAR[0] = t;
  for( size_t i = 0; i < _nx; ++i ) _dVAR[1+i] = y[i];
  for( size_t i = 0; i < _nx; ++i ) _dVAR[1+_nx+i] = yp[i];
  _dag->eval( _dWRK, _vEVAL.size(), _vEVAL.data(), _dEVAL.data(),
              _vVAR.size(), _vVAR.data(), _dVAR.data() );
}

inline
int
DAESLV_IDAS::IDRES__
( sunrealtype tt, N_Vector yy, N_Vector yp, N_Vector rr, void* )
{
  sunrealtype const* y  = N_VGetArrayPointer( yy );
  sunrealtype const* dy = N_VGetArrayPointer( yp );
  sunrealtype*       r  = N_VGetArrayPointer( rr );
  _eval_rows( tt, y, dy );
  for( size_t k = 0; k < _nres; ++k ) r[k] = _dEVAL[k];
  ++stats.numRES;
  return 0;
}

inline
bool
DAESLV_IDAS::solve
( std::vector<double> const& p, std::vector<double> const& tout )
{
  if( _vEVAL.empty() ){ _extractError = "solve() before a successful setup()"; return false; }
  if( p.size() != _np ){
    _extractError = std::to_string( p.size() ) + " parameter value(s) for " + std::to_string( _np )
                  + " parameter(s)";
    return false;
  }
  for( size_t j = 0; j < _np; ++j ) _dVAR[ 1 + 2*_nx + j ] = p[j];

  // initial values: the INITIAL rows, which may depend on the parameters
  std::vector<double> y0( _nx, 0. );
  _dVAR[0] = _dT[0];
  _dag->eval( _dWRK, _vIC.size(), _vIC.data(), y0.data(),
              _vVAR.size(), _vVAR.data(), _dVAR.data() );

  _init_stats( stats );
  SUNContext ctx = nullptr;
  if( SUNContext_Create( SUN_COMM_NULL, &ctx ) < 0 ){ _extractError = "SUNContext_Create failed"; return false; }
  N_Vector yy = N_VNew_Serial( (sunindextype)_nx, ctx );
  N_Vector yp = N_VNew_Serial( (sunindextype)_nx, ctx );
  N_Vector id = N_VNew_Serial( (sunindextype)_nx, ctx );
  for( size_t i = 0; i < _nx; ++i ){
    N_VGetArrayPointer( yy )[i] = y0[i];
    N_VGetArrayPointer( yp )[i] = 0.;
    N_VGetArrayPointer( id )[i] = _mID[i];
  }

  bool ok = true;
  std::string why;
  void* mem = IDACreate( ctx );
  // mirrors ODESLV_CVODES's switch, minus DIAG, which IDAS does not have
#if !defined( CRONOS__WITH_KLU )
  if( BASE_IDAS::options.LINSOL == BASE_IDAS::Options::SPARSE ){
    std::cerr << "DAESLV_IDAS ** LINSOL = SPARSE requires KLU, which this build does not include"
              << "; choose DENSE or DENSEDQ, or use a build with KLU" << std::endl;
    N_VDestroy( yy ); N_VDestroy( yp ); N_VDestroy( id ); IDAFree( &mem ); SUNContext_Free( &ctx );
    _extractError = "LINSOL = SPARSE without KLU";
    return false;
  }
#endif
  SUNMatrix A = nullptr;
  SUNLinearSolver LS = nullptr;
#if defined( CRONOS__WITH_KLU )
  if( BASE_IDAS::options.LINSOL == BASE_IDAS::Options::SPARSE ){
    A  = SUNSparseMatrix( (sunindextype)_nx, (sunindextype)_nx, (sunindextype)_nnzjacmax, CSC_MAT, ctx );
    LS = SUNLinSol_KLU( yy, A, ctx );
  }
  else
#endif
  {
    A  = SUNDenseMatrix( (sunindextype)_nx, (sunindextype)_nx, ctx );
    LS = SUNLinSol_Dense( yy, A, ctx );
  }
  _sun_mat = A; _sun_ls = LS;
  reregistration();
  do{
    if( IDAInit( mem, MC_IDRES__, _dT[0], yy, yp ) != IDA_SUCCESS ){ ok = false; why = "IDAInit"; break; }
    IDASStolerances( mem, BASE_IDAS::options.RTOL, BASE_IDAS::options.ATOL );
    if( IDASetLinearSolver( mem, LS, A ) != IDA_SUCCESS ){ ok = false; why = "IDASetLinearSolver"; break; }
    IDASetId( mem, id );
    if( BASE_IDAS::options.LINSOL != BASE_IDAS::Options::DENSEDQ )
      IDASetJacFn( mem, MC_IDJAC__ );         // analytic; DENSEDQ leaves IDAS its difference quotient
    IDASetSuppressAlg( mem, BASE_IDAS::options.SUPPRESSALG );
    IDASetMaxOrd( mem, BASE_IDAS::options.MAXORD );
    if( BASE_IDAS::options.NMAX ) IDASetMaxNumSteps( mem, (long)BASE_IDAS::options.NMAX );
    if( BASE_IDAS::options.HMAX > 0. ) IDASetMaxStep( mem, BASE_IDAS::options.HMAX );
    if( BASE_IDAS::options.H0 > 0. ) IDASetInitStep( mem, BASE_IDAS::options.H0 );
    IDASetMaxNumStepsIC( mem, BASE_IDAS::options.NMAXIC );
    IDASetMaxNumJacsIC( mem, BASE_IDAS::options.MAXJACIC );
    IDASetMaxNumItersIC( mem, BASE_IDAS::options.MAXNLIIC );
    IDASetLineSearchOffIC( mem, BASE_IDAS::options.LSOFFIC );
    if( BASE_IDAS::options.STEPTOLIC > 0. ) IDASetStepToleranceIC( mem, BASE_IDAS::options.STEPTOLIC );

    if( BASE_IDAS::options.CALCIC != BASE_IDAS::Options::NONE ){
      double t1 = BASE_IDAS::options.TOUT1;
      if( t1 == 0. ) t1 = tout.empty()? _dT[0] + 1.: tout.front();
      int const icopt = ( BASE_IDAS::options.CALCIC == BASE_IDAS::Options::YA_YDP )? IDA_YA_YDP_INIT: IDA_Y_INIT;
      if( IDACalcIC( mem, icopt, t1 ) != IDA_SUCCESS ){ ok = false; why = "IDACalcIC"; break; }
      IDAGetConsistentIC( mem, yy, yp );
      ++stats.numIC;
    }
    for( size_t i = 0; i < _nx; ++i ){
      _dX0[i]  = N_VGetArrayPointer( yy )[i];
      _dXP0[i] = N_VGetArrayPointer( yp )[i];
    }

    // The rows differ stage by stage exactly when an input is time-varying -- _mMASS then holds _ns
    // entries rather than one.  At such a boundary the BDF history describes a problem that no longer
    // exists (ODESLV's _stage_restart rule), AND, unlike an ODE, the state itself is no longer
    // consistent: the DIFFERENTIAL components carry over, but the ALGEBRAIC ones are whatever the new
    // constraints say and generally JUMP.  So the restart is IDAReInit followed by IDACalcIC in
    // IDA_YA_YDP_INIT, which holds the differential components and re-solves the algebraic ones and the
    // derivatives -- the split the id vector already carries.  Merely re-initialising would resume from
    // an inconsistent point, and with SUPPRESSALG on that first step can look plausible rather than fail.
    bool const restarts = ( _mMASS.size() > 1 );
    for( size_t sg = 0; sg < _ns && ok; ++sg ){
      double const ta = _dT[sg], tb = _dT[sg+1];
      if( sg && restarts ){
        if( !_build_rows( sg ) ){ ok = false; why = "_build_rows"; break; }
        if( IDAReInit( mem, ta, yy, yp ) != IDA_SUCCESS ){ ok = false; why = "IDAReInit"; break; }
#if defined( CRONOS__WITH_KLU )
        // the stage's rows changed, so the Jacobian's nonzero count or structure may have too: force a new
        // SYMBOLIC and numeric factorization, exactly as ODESLV_CVODES does at its own restarts
        if( BASE_IDAS::options.LINSOL == BASE_IDAS::Options::SPARSE
         && SUNLinSol_KLUReInit( _sun_ls, _sun_mat, (sunindextype)_nnzjacmax, 2 ) ){
          ok = false; why = "SUNLinSol_KLUReInit"; break;
        }
#endif
        if( BASE_IDAS::options.CALCIC != BASE_IDAS::Options::NONE ){
          int const icopt = ( BASE_IDAS::options.CALCIC == BASE_IDAS::Options::YA_YDP )? IDA_YA_YDP_INIT: IDA_Y_INIT;
          if( IDACalcIC( mem, icopt, tb ) != IDA_SUCCESS ){ ok = false; why = "IDACalcIC at a stage boundary"; break; }
          IDAGetConsistentIC( mem, yy, yp );
          ++stats.numIC;
        }
      }
      sunrealtype tret = ta;
      for( double const t : tout ){
        if( t <= ta || t > tb ) continue;
        if( IDASolve( mem, t, &tret, yy, yp, IDA_NORMAL ) != IDA_SUCCESS ){ ok = false; why = "IDASolve"; break; }
      }
      if( ok && sg + 1 < _ns && tret < tb )
        if( IDASolve( mem, tb, &tret, yy, yp, IDA_NORMAL ) != IDA_SUCCESS ){ ok = false; why = "IDASolve to the stage boundary"; }
    }
    if( !ok ) break;
    for( size_t i = 0; i < _nx; ++i ) _dX[i] = N_VGetArrayPointer( yy )[i];
    long n = 0;
    IDAGetNumSteps( mem, &n );    stats.numSteps = (unsigned)n;
    IDAGetNumJacEvals( mem, &n ); stats.numJAC   = (unsigned)n;
  } while( false );

  SUNLinSolFree( LS );
  SUNMatDestroy( A );
  N_VDestroy( yy ); N_VDestroy( yp ); N_VDestroy( id );
  IDAFree( &mem );
  SUNContext_Free( &ctx );
  _final_stats( stats );
  if( !ok ) _extractError = why + " failed";
  return ok;
}

} // end namespace mc

#endif
