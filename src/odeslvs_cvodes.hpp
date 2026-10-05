// Copyright (C) Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

#ifndef CRONOS__ODESLVS_CVODES_HPP
#define CRONOS__ODESLVS_CVODES_HPP

#undef  CRONOS__ODESLVS_CVODES_DEBUG

#include <sstream>
#include "odeslvs_base.hpp"
#include "odeslv_cvodes.hpp"

#define CRONOS__ODESLVS_CVODES_USE_BAD
#undef  CRONOS__ODESLVS_CVODES_DEBUG

namespace mc
{
//! @brief C++ class computing solutions of parametric ODEs with forward/adjoint sensitivity analysis capability using SUNDIALS and MC++.
////////////////////////////////////////////////////////////////////////
//! mc::ODESLV_CVODES is a C++ class for solution of IVPs in ODEs
//! with forward/adjoint sensitivity analysis capability using the code
//! CVODES in SUNDIALS and MC++.
////////////////////////////////////////////////////////////////////////
class ODESLVS_CVODES
: public virtual BASE_CVODES
, public virtual ODESLV_CVODES
, public virtual ODESLVS_BASE
, public virtual FFModel
{
 protected:
 
  using FFModel::_dag;

  using ODESLV_BASE::_print_interm;
  using ODESLV_BASE::_record;
  using ODESLV_BASE::_IC_D_QUAD;
  using ODESLV_BASE::_D2vec;
  using ODESLV_BASE::_Dx;
  using ODESLV_BASE::_Dq;
  using ODESLV_BASE::_Df;
      
  using ODESLVS_BASE::_ns;
  using ODESLVS_BASE::_istg;
  using ODESLVS_BASE::_t;
  using ODESLVS_BASE::_dT;
  using ODESLVS_BASE::_nx;
  using ODESLVS_BASE::_nq;
  using ODESLVS_BASE::_nf;
  using ODESLVS_BASE::_np;
  using ODESLVS_BASE::_ndxSEN;
  using ODESLVS_BASE::_nsen;
  using ODESLVS_BASE::_Dfp;
  using ODESLVS_BASE::_ny;
  using ODESLVS_BASE::_Dy;
  using ODESLVS_BASE::_Dyq;
  using ODESLVS_BASE::_vIC;
  using ODESLVS_BASE::_vRHS;
  using ODESLVS_BASE::_vQUAD;
  using ODESLVS_BASE::_vFCT;
  
  using ODESLV_BASE::_END_D_STA;
  
  using ODESLVS_BASE::_IC_SET_ASA;
  using ODESLVS_BASE::_CC_SET_ASA;
  using ODESLVS_BASE::_TC_SET_ASA;
  using ODESLVS_BASE::_RHS_SET_ASA;
  using ODESLVS_BASE::_IC_D_QUAD_ASA;
  using ODESLVS_BASE::_CC_D_QUAD_ASA;
  using ODESLVS_BASE::_TC_D_QUAD_ASA;
  using ODESLVS_BASE::_IC_SET_FSA;
  using ODESLVS_BASE::_CC_SET_FSA;
  using ODESLVS_BASE::_RHS_SET_FSA;
  using ODESLVS_BASE::_RHS_D_SET;
  using ODESLVS_BASE::_INI_D_SEN;
  using ODESLVS_BASE::_GET_D_SEN;
  using ODESLVS_BASE::_IC_D_SEN;
  using ODESLVS_BASE::_CC_D_SEN;
  using ODESLVS_BASE::_TC_D_SEN;
  using ODESLVS_BASE::_RHS_D_SEN;
  using ODESLVS_BASE::_RHS_D_QUAD;
  using ODESLVS_BASE::_JAC_D_SEN;
  using ODESLVS_BASE::_FCT_D_SEN;

  using ODESLV_CVODES::_cv_mem;
  using ODESLV_CVODES::_cv_flag;
  using ODESLV_CVODES::_Nx;
  using ODESLV_CVODES::_Nq;
  using ODESLV_CVODES::_pos_ic;
  using ODESLV_CVODES::_pos_rhs;
  using ODESLV_CVODES::_pos_quad;
  using ODESLV_CVODES::_pos_fct;
  using ODESLV_CVODES::_END_STA;
  using ODESLV_CVODES::_states;
  using ODESLV_CVODES::_states_stage;
  using ODESLV_CVODES::MC_CVRHS__;
  using ODESLV_CVODES::MC_CVQUAD__;
  using ODESLV_CVODES::MC_CVJAC__;
  using ODESLV_CVODES::_xk;
  using ODESLV_CVODES::_qk;
  using ODESLV_CVODES::_f;
  using ODESLV_CVODES::_nnzjac;

  //! @brief SUNMatrix, SUNLinearSolver and SUNNonlinearSolver of EACH backward problem (index ifct).  SUNDIALS
  //! requires one of each per backward problem: sharing them across the _nf backward problems (as before
  //! 2026-09-27) corrupted the KLU factorisation and crashed CVodeB (ODESLV_NIFTE, SPARSE, nf = 2).
  std::vector<SUNMatrix>          _sun_matB;
  std::vector<SUNLinearSolver>    _sun_lsB;
  std::vector<SUNNonlinearSolver> _sun_nlsB;

  //! @brief Free the backward solvers of every backward problem
  void _free_backward_solvers
    ()
    {
      for( auto& n : _sun_nlsB ) if( n ){ SUNNonlinSolFree( n ); n = nullptr; }
      for( auto& l : _sun_lsB )  if( l ){ SUNLinSolFree( l );    l = nullptr; }
      for( auto& m : _sun_matB ) if( m ){ SUNMatDestroy( m );    m = nullptr; }
      _sun_nlsB.clear();  _sun_lsB.clear();  _sun_matB.clear();
    }

  //! @brief SUNNonlinearSolver object for use by CVodeS forward
  SUNNonlinearSolver _sun_nlsF;

  //! @brief current function index
  unsigned _ifct;

  //! @brief current parameter sensitivity index
  unsigned _isen;

  //! @brief size of N_vector arrays
  unsigned _nvec;

  //! @brief N_Vector object holding current adjoints
  N_Vector* _Ny;

  //! @brief N_Vector object holding current quadratures
  N_Vector* _Nyq;

  //! @brief pointer to array holding identifiers of the backward problems
  int* _indexB;
  
  //! @brief pointer to array holding identifiers of the backward problems
  unsigned* _iusrB;

  //! @brief position in _Ny
  unsigned _pos_adj;

  //! @brief position in _Nq
  unsigned _pos_adjquad; 

  //! @brief state sensitivity values at stage times
  std::vector< std::vector< std::vector< double > > > _xpk;

  //! @brief state adjoint values at stage times
  std::vector< std::vector< std::vector< double > > > _lk;

  //! @brief quadrature sensitivity values at stage times
  std::vector< std::vector< std::vector< double > > > _qpk;

  //! @brief function derivatives
  std::vector< std::vector< double > > _fp;

 public:

  using FFModel::set_model;

  /** @ingroup ODESLV
   *  @{
   */
  typedef typename ODESLV_BASE::STATUS STATUS;
  typedef typename ODESLV_BASE::Results Results;
  typedef typename ODESLV_CVODES::Options Options;
  typedef typename ODESLV_CVODES::Exceptions Exceptions;
  using ODESLV_BASE::np;
  using ODESLV_BASE::nf;
  using ODESLV_CVODES::options;
  using ODESLV_CVODES::solve;
  using ODESLV_CVODES::val_state;
  using ODESLV_CVODES::val_quadrature;
  using ODESLV_CVODES::val_function;
  using ODESLV_CVODES::results_solve;
  using ODESLV_CVODES::stats_solve;
  
  //! @brief Constructor taking the DAG.  FFModel is a VIRTUAL base, so the MOST-DERIVED class must initialise
  //! it -- delegating and then calling FFModel::set( dag ) is what FFModel's own constructor does.
  ODESLVS_CVODES
    ( FFGraph* dag )
    : ODESLVS_CVODES()
    { FFModel::set( dag ); }

  //! @brief Default constructor
  ODESLVS_CVODES
    ();

  //! @brief Virtual destructor
  virtual ~ODESLVS_CVODES
    ();

  //! @brief Statistics of the last forward sensitivity integration (solve_fsens)
  Stats stats_fsens;

  //! @brief Statistics of the last adjoint integration (solve_asens)
  Stats stats_asens;

  //! @brief Forward sensitivity trajectories of the last solve_fsens (see Options::RESRECORD): one list of records
  //! per sensitivity direction, each record holding the state then the quadrature sensitivities
  std::vector< std::vector< Results > > results_fsens;

  //! @brief Adjoint trajectories of the last solve_asens (see Options::RESRECORD), in BACKWARD time: one list of
  //! records per output function, each record holding the adjoint states then the running gradient
  std::vector< std::vector< Results > > results_asens;

 //! @brief Propagate states and state-sensitivities forward in time through every time stages
  STATUS solve_fsens
    ( std::vector<double> const& p, std::vector<double> const& c=std::vector<double>(), std::ostream& os=std::cout );

 //! @brief Propagate states and state-sensitivities forward in time through every time stages
  STATUS solve_fsens
    ( double const* p, double const* c=nullptr, std::ostream& os=std::cout );

  //! @brief solve_fsens with the input and constant values given BY NAME (see ODESLV_CVODES::solve).
  STATUS solve_fsens
    ( std::vector<FFModel::InputVal> const& vIn, std::ostream& os=std::cout )
    {
      std::vector<double> P, c;  std::string err;
      if( !assemble_values( vIn, P, c, err ) ) throw std::invalid_argument( "ODESLVS::solve_fsens ** " + err + "\n" );
      return solve_fsens( P, c, os );
    }

  //! @brief Propagate states and adjoints forward and backward in time through every time stages
  STATUS solve_asens
    ( std::vector<double> const& p, std::vector<double> const& c=std::vector<double>(), std::ostream& os=std::cout );

  //! @brief Propagate states and adjoints forward and backward in time through every time stages
  //! KNOWN LIMITATION (2026-09-27): the adjoint is WRONG for models whose initial conditions are given PER STAGE
  //! (per-stage initial values: add_transition) -- even an identity transition x -> x.
  //! Single-stage models and multi-stage models with ONE initial condition (state continuous across stages) are
  //! correct.  Use solve_fsens for per-stage initial conditions.  To be revisited with FFModel::add_interface.
  //! See KNOWN_ISSUE_ASA_per_stage_IC.md; reproducer probe_asa_stage2.cpp.
  STATUS solve_asens
    ( double const* p, double const* c=nullptr, std::ostream& os=std::cout );

  //! @brief solve_asens with the input and constant values given BY NAME (see ODESLV_CVODES::solve).
  STATUS solve_asens
    ( std::vector<FFModel::InputVal> const& vIn, std::ostream& os=std::cout )
    {
      std::vector<double> P, c;  std::string err;
      if( !assemble_values( vIn, P, c, err ) ) throw std::invalid_argument( "ODESLVS::solve_asens ** " + err + "\n" );
      return solve_asens( P, c, os );
    }

  //! @brief Setup local copy of parametric ODEs
  bool setup
    ()
    { std::lock_guard<std::recursive_mutex> dag_lock_( dag_mutex() );   // see FFModel::dag_mutex()
      // as ODESLV_CVODES::setup(): split `options` (inherited from ODESLV_CVODES) into the two parents, set the
      // model up if it has not been, then take the local copy of the description.
      FFModel::options     = static_cast<FFModel::Options const&>( options );
      BASE_CVODES::options = static_cast<BASE_CVODES::Options const&>( options );
      // See ODESLV_CVODES::setup(): the DESCRIPTION path populates no FFModel model, so its validation
      // must be skipped -- _SETUP() takes the staging vectors as given.
      if( !FFModel::is_setup() && !FFModel::setup() ){
        _extractError = std::string( "model setup refused: " ) + setup_status_str( setup_status() );
        return false;
      }
      return ODESLVS_BASE::_SETUP(); }

  //! @brief Setup as an INDEPENDENT COPY of @p IVP: the same description (see _copy_description), the same
  //! options, then a full setup().  The user DAG @p IVP was declared on must outlive this solver.
  bool setup
    ( ODESLVS_CVODES const& IVP )
    { std::lock_guard<std::recursive_mutex> dag_lock_( dag_mutex() );  _copy_description( IVP ); options = IVP.options; return setup(); }

  using FFModel::fdiff;

  //! @brief The forward SENSITIVITY MODEL of this problem along the controls @p vU, as a new solver on the same
  //! user DAG with this solver's options -- FFModel::fdiff in returning form.  @p nDir directions (0: one per
  //! DOF of @p vU, the full Jacobian in one solve); the product's outputs are [F | dF^(1) | ... | dF^(nDir)],
  //! its direction inputs are registered controls, and FFModel::fdiff_seed() seeds direction k.  Returns
  //! nullptr, with the reason in @p err, if the model cannot be differentiated.  The caller owns the result.
  ODESLVS_CVODES* fdiff
    ( std::vector<FFVar> const& vU, size_t const nDir,
      std::vector<std::vector<FFVar>>& vDU, std::vector<std::vector<FFVar>>& vS, std::string& err )
    const
    {
      ODESLVS_CVODES* sens = new ODESLVS_CVODES( _usr._dagUsr );
      sens->options = options;
      if( FFModel::fdiff( *sens, vU, nDir, vDU, vS, err ) ) return sens;
      delete sens;  return nullptr;
    }

  //! @brief LEGACY INTERFACE ON THE NEW ENGINE (2026-09-27).  Symbolic forward differentiation of this problem
  //! w.r.t. the PARAMETERS @p pPar -- entries of var_parameter(), minted levels included.  The product is a new
  //! solver on the same user DAG with EXACTLY this problem's parameters; its states are the originals plus one
  //! sensitivity set per parameter, and its nPar*nf functions are the gradient components
  //! f[i*nf+j] = d f_j / d p_i.  Built on FFModel::fdiff: each parameter is mapped back to its declared input
  //! and DOF, one direction is declared per parameter, and every direction input is FIXED to its unit seed so
  //! none survives extraction.  This is the form FFODESLV::deriv plugs into SYMDIFF.  Returns nullptr with the
  //! reason in @p err on failure; the caller owns the result.
  ODESLVS_CVODES* fdiff
    ( size_t const nPar, FFVar const* pPar, std::string& err )
    const
    {
      err.clear();
      // 1. each entry -> a parameter (declared input, DOF) or a CONSTANT; distinct inputs / constants in order
      std::vector<FFVar> vIn( nPar );  std::vector<size_t> vDof( nPar );  std::vector<int> vCi( nPar, -1 );
      std::vector<FFVar> vU, vC;
      auto const& vCst = var_constant();
      for( size_t i = 0; i < nPar; ++i ){
        if( input_of_parameter( pPar[i], vIn[i], vDof[i] ) ){
          bool seen = false; for( auto const& u : vU ) if( u.id().second == vIn[i].id().second ){ seen = true; break; }
          if( !seen ) vU.push_back( vIn[i] );
          continue;
        }
        int ic = -1;
        for( size_t c = 0; c < vCst.size(); ++c ) if( vCst[c].id().second == pPar[i].id().second ){ ic = (int)c; break; }
        if( ic < 0 ){ err = "fdiff: " + pPar[i].name() + " is neither a parameter nor a constant of the extracted model"; return nullptr; }
        FFVar const& cdecl = _usr._vCstUsr[ic];            // constants are extracted in declared order
        int pos = -1; for( size_t c = 0; c < vC.size(); ++c ) if( vC[c].id().second == cdecl.id().second ){ pos = (int)c; break; }
        if( pos < 0 ){ pos = (int)vC.size(); vC.push_back( cdecl ); }
        vCi[i] = pos;
      }
      // 2. one direction per entry, legacy output layout; constant directions by literal seed
      std::vector<std::vector<double>> seedC( nPar, std::vector<double>( vC.size(), 0. ) );
      for( size_t i = 0; i < nPar; ++i ) if( vCi[i] >= 0 ) seedC[i][ vCi[i] ] = 1.;
      std::vector<std::vector<FFVar>> vDU, vS;
      ODESLVS_CVODES* sens = new ODESLVS_CVODES( _usr._dagUsr );
      sens->options = options;
      if( !FFModel::fdiff( *sens, vU, nPar, vDU, vS, err, /*keep_originals=*/false, vC, seedC ) ){ delete sens; return nullptr; }
      // 3. direction i: unit seed at its DOF if it is along an input, zero in every input otherwise
      for( size_t i = 0; i < nPar; ++i )
        for( size_t c = 0; c < vU.size(); ++c ){
          std::vector<double> e( sens->control_ndof( vDU[i][c] ), 0. );
          if( vCi[i] < 0 && vU[c].id().second == vIn[i].id().second && vDof[i] < e.size() ) e[ vDof[i] ] = 1.;
          if( !sens->fix_input( vDU[i][c], e ) ){ err = "fdiff: cannot fix the direction input " + vDU[i][c].name(); delete sens; return nullptr; }
        }
      return sens;
    }
  ODESLVS_CVODES* fdiff
    ( size_t const nPar, FFVar const* pPar )
    const
    { std::string err; ODESLVS_CVODES* p = fdiff( nPar, pPar, err ); if( !p ) std::cerr << "  **ERROR: " << err << std::endl; return p; }
  ODESLVS_CVODES* fdiff
    ( std::vector<FFVar> const& vPar )
    const
    { return fdiff( vPar.size(), vPar.data() ); }

  //! @brief The same, addressed by DECLARED INPUTS: every DOF of each input in @p vU, in order, resolved through
  //! parameter_index().
  ODESLVS_CVODES* fdiff
    ( std::vector<FFVar> const& vU, std::string& err )
    const
    {
      std::vector<FFVar> vPar;
      for( auto const& u : vU ){
        auto const ndx = parameter_index( u );
        if( ndx.empty() ){ err = "fdiff: " + u.name() + " is not an input of the extracted model (call setup() first)"; return nullptr; }
        for( size_t i : ndx ) vPar.push_back( _mP[i] );
      }
      return fdiff( vPar.size(), vPar.data(), err );
    }

  //! @brief Same, along the REGISTERED CONTROLS, one direction per DOF: the full Jacobian in one solve.
  ODESLVS_CVODES* fdiff
    ( std::vector<std::vector<FFVar>>& vDU, std::vector<std::vector<FFVar>>& vS, std::string& err )
    const
    {
      std::vector<FFVar> vU;
      for( auto const& [u,spec] : controls() ) vU.push_back( u );
      return fdiff( vU, 0, vDU, vS, err );
    }

  //! @brief Legacy layout over the REGISTERED CONTROLS -- the same directions, in the same order, as
  //! solve_fsens() and solve_asens(), so f[i*nf+j] lines up with val_function_gradient()[i][j].
  ODESLVS_CVODES* fdiff
    () const
    {
      std::vector<FFVar> vPar( _nsen );
      for( size_t i=0; i<_nsen; ++i ) vPar[i] = _mP[_ndxSEN[i]];
      return fdiff( vPar );
    }
    
  //! @brief Record state and sensitivity trajectories in files <a>obndsta</a> and <a>obndsa</a>, with accuracy of <a>iprec</a> digits
  void record
    ( std::ofstream& obndsta, std::ofstream* obndsen, unsigned const iprec=5 )
    const
    { this->ODESLV_CVODES::record( obndsta, iprec );
      for( unsigned isen=0; isen<results_fsens.size(); ++isen )
        this->ODESLV_BASE::_record( obndsen[isen], results_fsens[isen], iprec ); }

  //! @brief Record state trajectories in files <a>obndsta</a>, with accuracy of <a>iprec</a> digits
  void record
    ( std::ofstream& obndsta, unsigned const iprec=5 )
    const
    { ODESLV_CVODES::record( obndsta, iprec ); }

  //! @brief Retreive state sensitivity values at stage times
  std::vector< std::vector< std::vector< double > > > const& val_state_sensitivity
    ()
    const
    { return _xpk; }

  //! @brief Retreive adjoint sensitivity values at stage times
  std::vector< std::vector< std::vector< double > > > const& val_state_adjoint
    ()
    const
    { return _lk; }

  //! @brief Retreive quadrature sensitivity values at stage times
  std::vector< std::vector< std::vector< double > > > const& val_quadrature_sensitivity
    ()
    const
    { return _qpk; }

  //! @brief Retreive function values at stage times
  std::vector< std::vector< double > > const& val_function_gradient
    ()
    const
    { return _fp; }
  /** @} */

  //! @brief The gradient rows of control @p var: [dof][function], its DOFs in control_dofs() order, from the last
  //! solve_fsens() or solve_asens().  Throws if @p var is not a registered control (it then has no
  //! sensitivity direction) or no sensitivity solve has been run.
  std::vector< std::vector< double > > val_function_gradient
    ( FFVar const& var )
    const
    {
      auto const b = control_block( var );
      auto const& G = val_function_gradient();
      if( !b.ndof )
        throw std::invalid_argument( "ODESLVS::val_function_gradient ** " + var.name() + " is not a registered control\n" );
      if( G.size() < b.offset + b.ndof )
        throw std::invalid_argument( "ODESLVS::val_function_gradient ** no sensitivity solve over the current controls\n" );
      return std::vector< std::vector< double > >( G.begin()+b.offset, G.begin()+b.offset+b.ndof );
    }

  //! @brief The gradient rows of control @p var on ONE element: [node][function], the n_node DOFs of element
  //! @p ndx_el (keyed by the evolution domain, as pos_input takes it).  A time-invariant control ignores @p ndx_el.
  std::vector< std::vector< double > > val_function_gradient
    ( FFVar const& var, std::map<FFVar,size_t,lt_FFVar> const& ndx_el )
    const
    {
      auto const all = val_function_gradient( var );
      size_t const nn = size_input( var );
      if( all.size() == 1 ) return all;                                   // time-invariant
      size_t e = 0;
      if( !_mT.empty() ){ auto const it = ndx_el.find( _mT[0] ); if( it != ndx_el.cend() ) e = it->second; }
      if( (e+1)*nn > all.size() )
        throw std::invalid_argument( "ODESLVS::val_function_gradient ** element out of range for " + var.name() + "\n" );
      return std::vector< std::vector< double > >( all.begin()+e*nn, all.begin()+(e+1)*nn );
    }

 protected:
  //! @brief Propagate states and state-sensitivities forward in time through every time stages
  STATUS _states_FSA
    ( double const* p, double const* c, std::ostream& os );

  //! @brief Propagate states and adjoints forward and backward in time through every time stages
  STATUS _states_ASA
    ( double const* p, double const* c, std::ostream& os );

 private:
  //! @brief Function to initialize CVodes memory block (virtual)
  virtual bool _INI_CVODE
    ();

  //! @brief Function to initialize CVodeS memory block for forward sensitivity
  bool _INI_CVODES_FSA
    ();

  //! @brief Function to initialize CVodeS memory block for adjoint sensitivity
  bool _INI_CVODES_ASA
    ( unsigned const ifct, int& indexB, unsigned& iusrB );

  //! @brief Function to reinitialize CVodeS memory block for forward sensitivity
  bool _CC_CVODES_FSA
    ();

  //! @brief Function to reinitialize CVodeS memory block for forward quadrature sensitivity
  bool _CC_CVODES_QUAD
    ();

  //! @brief Function to reinitialize CVodeS memory block for adjoint sensitivity
  bool _CC_CVODES_ASA
    ( unsigned const ifct, int const indexB );

  //! @brief Function to finalize sensitivity/adjoint bounding, closing the statistics @p stats
  void _END_SEN
    ( Stats& stats );

  //! @brief Function to reinitialize sensitivity analysis
  bool _REINI_SEN
    ();

  //! @brief Function to initialize adjoint sensitivity analysis
  //! @brief Arm the adjoint machinery, after the forward sweep
  bool _ARM_ASA
    ();

  //! @brief Arm the adjoint inside _INI_CVODE (one-go CVodeF sweep) rather than after the sweep.  Set by
  //! _states_ASA from _sweep_restarts(): false when a stage boundary restarts the integrator.
  bool _armOnInit = false;

  bool _INI_ASA
    ( double const* p );

  //! @brief Static wrapper to calculate the adjoint ODEs RHS derivatives
  static int MC_CVRHSB__
    ( sunrealtype t, N_Vector x, N_Vector y, N_Vector ydot, void* user_data );

  //! @brief Virtual function to calculate the adjoint ODEs RHS derivatives
  virtual int CVRHSB__
    ( sunrealtype t, N_Vector x, N_Vector y, N_Vector ydot, void* user_data );

  //! @brief Static wrapper to calculate the adjoint quadrature RHS derivatives
  static int MC_CVQUADB__
    ( sunrealtype t, N_Vector x, N_Vector y, N_Vector qdot, void* user_data );

  //! @brief Virtual function to calculate the adjoint quadrature RHS derivatives
  virtual int CVQUADB__
    ( sunrealtype t, N_Vector x, N_Vector y, N_Vector qdot, void* user_data );

  //! @brief Static wrapper to calculate the adjoint ODEs RHS Jacobian
  static int MC_CVJACB__
    ( sunrealtype t, N_Vector y, N_Vector yB, N_Vector fyB, SUNMatrix JacB,
      void* user_dataB, N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B );

  //! @brief Virtual function to calculate the adjoint ODEs RHS Jacobian
  virtual int CVJACB__
    ( sunrealtype t, N_Vector y, N_Vector yB, N_Vector fyB, SUNMatrix JacB,
      void* user_dataB, N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B );

  //! @brief Function to initialize forward sensitivity analysis
  bool _INI_FSA
    ( double const* p );

  //! @brief Static wrapper to calculate the sensitivity ODEs RHS derivatives
  static int MC_CVRHSF__
    ( int Ns, sunrealtype t, N_Vector x, N_Vector xdot, int is, N_Vector y,
      N_Vector ydot, void* user_data, N_Vector tmp1, N_Vector tmp2 );

  //! @brief Pure virtual function to calculate the sensitivity ODEs RHS derivatives
  virtual int CVRHSF__
    ( int Ns, sunrealtype t, N_Vector x, N_Vector xdot, int is, N_Vector y,
      N_Vector ydot, void* user_data, N_Vector tmp1, N_Vector tmp2 );

  //! @brief Static wrapper to calculate the sensitivity quadrature RHS derivatives
  static int MC_CVQUADF__
    ( int Ns, sunrealtype t, N_Vector x, N_Vector* y, N_Vector qdot, N_Vector* qSdot, 
      void *user_data, N_Vector tmp1, N_Vector tmp2 );

  //! @brief Pure virtual function to calculate the sensitivity quadrature RHS derivatives
  virtual int CVQUADF__
    ( int Ns, sunrealtype t, N_Vector x, N_Vector* y, N_Vector qdot, N_Vector* qSdot, 
      void *user_data, N_Vector tmp1, N_Vector tmp2 );

  //! @brief Private methods to block default compiler methods
  ODESLVS_CVODES( ODESLVS_CVODES const& ) = delete;
  ODESLVS_CVODES& operator=( ODESLVS_CVODES const& ) = delete;
};

inline
ODESLVS_CVODES::ODESLVS_CVODES
()
: _sun_nlsF(nullptr),
  _ifct(0), _isen(0), _nvec(0), _Ny(nullptr), _Nyq(nullptr),
  _indexB(nullptr), _iusrB(nullptr)
{}

inline
ODESLVS_CVODES::~ODESLVS_CVODES
()
{
  if( _Ny )  N_VDestroyVectorArray( _Ny,  _nvec );
  if( _Nyq ) N_VDestroyVectorArray( _Nyq, _nvec );
  delete[] _indexB;
  delete[] _iusrB;
  if( _sun_nlsF ) SUNNonlinSolFree( _sun_nlsF ); /* Free the nonlinear solver memory */
  _free_backward_solvers();
}

//! @brief Arm the adjoint machinery.  Must be called AFTER the forward sweep: CVodeAdjInit'd memory is
//! invalidated by any plain CVode step, and neither CVodeAdjReInit nor a fresh CVodeAdjInit recovers it.
inline
bool
ODESLVS_CVODES::_ARM_ASA
()
{
  _cv_flag = CVodeAdjInit( _cv_mem, options.ASACHKPT, options.ASAINTERP );
  if( _check_cv_flag( &_cv_flag, "CVodeAdjInit", 1 ) ) return false;
  return true;
}

inline
bool
ODESLVS_CVODES::_INI_CVODE
()
{
  // Call _INI_CVODE in ODEBND_CVODES
  this->ODESLV_CVODES::_INI_CVODE();

  // A solver armed by CVodeAdjInit must never be stepped by plain CVode -- that poisons it unrecoverably
  // (see _ARM_ASA).  So the adjoint is armed here, right after CVodeCreate, ONLY when the sweep will be
  // one uninterrupted CVodeF (_armOnInit, set by _states_ASA when nothing restarts); a restarting model
  // runs its sweep unarmed on plain CVode and _states_ASA arms afterwards.
  if( _armOnInit && !_ARM_ASA() ) return false;

  // Reinitialize adjoint holding vectors
  delete[] _indexB; _indexB = new int[_nf];
  delete[] _iusrB;  _iusrB  = new unsigned[_nf];

  return true;
}

inline
bool
ODESLVS_CVODES::_INI_CVODES_ASA
( unsigned const ifct, int& indexB, unsigned& iusrB )
{
  // Create CVodeS memory block for the ADAMS or BDF method
  _cv_flag = CVodeCreateB( _cv_mem, options.INTMETH, &indexB );
  if( _check_cv_flag( &_cv_flag, "CVodeCreateB", 1 ) ) return false;

  // Specify error output
//  if( options.DISPLAY < 0 )
//    _cv_flag = CVodeSetErrFile( _cv_mem, NULL );
//  else
//    _cv_flag = CVodeSetErrFile( _cv_mem, stderr );
//  if( _check_cv_flag( &_cv_flag, "CVodeSetErrFile", 1 ) ) return false;

  // Initialize CVodeS memory and specify the adjoint RHS function,
  // terminal time _t, and terminal adjoint _Ny
  _cv_flag = CVodeInitB( _cv_mem, indexB, MC_CVRHSB__, _t, _Ny[ifct] );
  if( _check_cv_flag( &_cv_flag, "CVodeInitB", 1 ) ) return false;

  // Specify the user_data to pass the function index corresponding to
  // the RHS function
  iusrB = _ifct;
  _cv_flag = CVodeSetUserDataB( _cv_mem, indexB, &iusrB );
  if( _check_cv_flag( &_cv_flag, "CVodeSetUserDataB", 1 ) ) return false;

  // Specify the nonlinear solver -- ONE matrix, linear solver and nonlinear solver PER backward problem
  if( !ifct ) _free_backward_solvers();
  if( _sun_matB.size() <= ifct ){ _sun_matB.resize( ifct+1, nullptr ); _sun_lsB.resize( ifct+1, nullptr ); _sun_nlsB.resize( ifct+1, nullptr ); }
  switch( options.NLINSOL ){
   // Fixed point nonlinear solver
   case Options::FIXEDPOINT:
    _sun_nlsB[ifct] = SUNNonlinSol_FixedPoint( _Ny[ifct], 0, sunctx );
    if( _check_cv_flag( (void *)_sun_nlsB[ifct], "SUNNonlinSol_FixedPoint", 0 ) ) return false;
    break;

   // Newton nonlinear solver
   case Options::NEWTON:
    // Specify the linear solver and Jacobian approximation
    switch( options.LINSOL ){
     case Options::DIAG: default:
       _cv_flag = CVDiagB( _cv_mem, indexB );
       if( _check_cv_flag( &_cv_flag, "CVDiag", 1) ) return false;
       break;

     // Dense Jacobian
     case Options::DENSE:
     case Options::DENSEDQ:
       _sun_matB[ifct] = SUNDenseMatrix( _ny, _ny, sunctx );
       if( _check_cv_flag( (void*)_sun_matB[ifct], "SUNDenseMatrix", 0 ) ) return false;
       _sun_lsB[ifct] = SUNLinSol_Dense( _Ny[ifct], _sun_matB[ifct], sunctx );
       if( _check_cv_flag( (void *)_sun_lsB[ifct], "SUNLinSol_Dense", 0 ) ) return false;
       _cv_flag = CVodeSetLinearSolverB( _cv_mem, indexB, _sun_lsB[ifct], _sun_matB[ifct] );
       if( _check_cv_flag( &_cv_flag, "CVodeSetLinearSolverB", 1 ) ) return false;
       _cv_flag = CVodeSetJacFnB( _cv_mem, indexB, options.LINSOL==Options::DENSE? MC_CVJACB__: nullptr );
       if ( _check_cv_flag( &_cv_flag, "CVodeSetJacFnB", 1 ) ) return false;
       break;

#if defined( CRONOS__WITH_KLU )
     // Sparse Jacobian
     case Options::SPARSE:
       _sun_matB[ifct] = SUNSparseMatrix( _ny, _ny, _nnzjac, CSR_MAT, sunctx );
       if( _check_cv_flag( (void*)_sun_matB[ifct], "SUNSparseMatrix", 0 ) ) return false;
       _sun_lsB[ifct] = SUNLinSol_KLU( _Ny[ifct], _sun_matB[ifct], sunctx );
       if( _check_cv_flag( (void *)_sun_lsB[ifct], "SUNLinSol_KLU", 0 ) ) return false;
       _cv_flag = CVodeSetLinearSolverB( _cv_mem, indexB, _sun_lsB[ifct], _sun_matB[ifct] );
       if( _check_cv_flag( &_cv_flag, "CVodeSetLinearSolverB", 1 ) ) return false;
       _cv_flag = CVodeSetJacFnB( _cv_mem, indexB, MC_CVJACB__ );
       if ( _check_cv_flag( &_cv_flag, "CVodeSetJacFnB", 1 ) ) return false;
       break;
#endif
    }
    _sun_nlsB[ifct] = SUNNonlinSol_Newton( _Ny[ifct], sunctx );
    if( _check_cv_flag( (void *)_sun_nlsB[ifct], "SUNNonlinSol_Newton", 0 ) ) return false;
    break;
  }
  _cv_flag = CVodeSetNonlinearSolverB( _cv_mem, indexB, _sun_nlsB[ifct] );
  if( _check_cv_flag( &_cv_flag, "CVodeSetNonlinearSolverB", 1 ) ) return false;

  // Specify the relative and absolute tolerances for states
  _cv_flag = CVodeSStolerancesB( _cv_mem, indexB, options.RTOLB, options.ATOLB );
  if( _check_cv_flag( &_cv_flag, "CVodeSStolerancesB", 1 ) ) return false;

  // Specify minimum stepsize
  _cv_flag = CVodeSetMinStepB( _cv_mem, indexB, options.HMIN>0.? options.HMIN:0. );
  if( _check_cv_flag( &_cv_flag, "CVodeSetMinStepB", 1 ) ) return false;

  // Specify maximum stepsize
  _cv_flag = CVodeSetMaxStepB( _cv_mem, indexB, options.HMAX>0.? options.HMAX: 0. );
  if( _check_cv_flag( &_cv_flag, "CVodeSetMaxStepB", 1 ) ) return false;

  // Specify maximum number of steps between two stage times
  _cv_flag = CVodeSetMaxNumStepsB( _cv_mem, indexB, options.NMAX );
  if( _check_cv_flag( &_cv_flag, "CVodeSetMaxNumStepsB", 1 ) ) return false;

  // Initialize the integrator memory for the quadrature variables
  if( !_Nyq ) return true;
  _cv_flag = CVodeQuadInitB( _cv_mem, indexB, MC_CVQUADB__, _Nyq[ifct] );
  if( _check_cv_flag( &_cv_flag, "CVodeQuadInitB", 1 ) ) return false;

  // Specify whether or not to perform error control on quadrature
  _cv_flag = CVodeSetQuadErrConB( _cv_mem, indexB, options.QERRB );
  if( _check_cv_flag( &_cv_flag, "CVodeSetQuadErrConB", 1 ) ) return false;

  // Specify the relative and absolute tolerances for quadratures
  _cv_flag = CVodeQuadSStolerancesB( _cv_mem, indexB, options.RTOLB, options.ATOLB );
  if( _check_cv_flag( &_cv_flag, "CVodeQuadSStolerancesB", 1 ) ) return false;

  return true;
}

inline
bool
ODESLVS_CVODES::_INI_CVODES_FSA
()
{
  // Allocate memory for sensitivity integration
  _cv_flag = CVodeSensInit1( _cv_mem, _nsen, options.FSACORR, MC_CVRHSF__, _Ny );
  if( _check_cv_flag( &_cv_flag, "CVodeSensInit1", 1 ) ) return false;

  // Specify error output
//  if( options.DISPLAY < 0 )
//    _cv_flag = CVodeSetErrFile( _cv_mem, NULL );
//  else
//    _cv_flag = CVodeSetErrFile( _cv_mem, stderr );
//  if( _check_cv_flag( &_cv_flag, "CVodeSetErrFile", 1 ) ) return false;

  // Specify absolute and relative tolerances for sensitivities
  if( options.AUTOTOLS ){
    _cv_flag = CVodeSensEEtolerances( _cv_mem );
    if( _check_cv_flag( &_cv_flag, "CVodeSensEEtolerances", 1) ) return false;
  }
  else{
    std::vector<sunrealtype> ATOLS( _nsen, options.ATOLS );
    _cv_flag = CVodeSensSStolerances( _cv_mem, options.RTOLS, ATOLS.data() );
    if( _check_cv_flag( &_cv_flag, "CVodeSensSStolerances", 1) ) return false;
  }

  // Specify the error control strategy for sensitivity variables
  _cv_flag = CVodeSetSensErrCon( _cv_mem, options.FSAERR );
  if( _check_cv_flag( &_cv_flag, "CVodeSetSensErrCon", 1 ) ) return false;

  // Specify problem parameter information for sensitivity calculations
  _cv_flag = CVodeSetSensParams( _cv_mem, 0, 0, 0 );
  if( _check_cv_flag( &_cv_flag, "CVodeSetSensParams", 1 ) ) return false;

  // Specify the nonlinear solver for sensitivity calculations
  if( _sun_nlsF ){ SUNNonlinSolFree( _sun_nlsF );  _sun_nlsF = nullptr; } /* Free the nonlinear solver memory */
  switch( options.NLINSOL ){
   // Fixed point nonlinear solver
   case Options::FIXEDPOINT:
    switch( options.FSACORR ){
     case Options::SIMULTANEOUS:
      _sun_nlsF = SUNNonlinSol_FixedPointSens( _nsen+1, _Nx, 0, sunctx);
      break;
     case Options::STAGGERED:
      _sun_nlsF = SUNNonlinSol_FixedPointSens( _nsen, _Nx, 0, sunctx);
      break;
     case Options::STAGGERED1:
      _sun_nlsF = SUNNonlinSol_FixedPoint( _Nx, 0, sunctx);
      break;
    }
    break;
   
   // Newton nonlinear solver
   case Options::NEWTON:
    switch( options.FSACORR ){
     case Options::SIMULTANEOUS:
      _sun_nlsF = SUNNonlinSol_NewtonSens( _nsen+1, _Nx, sunctx);
      break;
     case Options::STAGGERED:
      _sun_nlsF = SUNNonlinSol_NewtonSens( _nsen, _Nx, sunctx);
      break;
     case Options::STAGGERED1:
      _sun_nlsF = SUNNonlinSol_Newton( _Nx, sunctx);
      break;
    }
    if( _check_cv_flag( (void *)_sun_nlsF, "SUNNonlinSol_FixedPointSens", 0 ) ) return false;
  }
  
  // Attach sensitivity nonlinear solver to CVodeS
  switch( options.FSACORR ){
   case Options::SIMULTANEOUS:
    _cv_flag = CVodeSetNonlinearSolverSensSim( _cv_mem, _sun_nlsF );
    if( _check_cv_flag( &_cv_flag, "CVodeSetNonlinearSolverSensSim", 1 ) ) return false;
    break;
   case Options::STAGGERED:
    _cv_flag = CVodeSetNonlinearSolverSensStg( _cv_mem, _sun_nlsF );
    if( _check_cv_flag( &_cv_flag, "CVodeSetNonlinearSolverSensStg", 1 ) ) return false;
    break;
   case Options::STAGGERED1:
    _cv_flag = CVodeSetNonlinearSolverSensStg1( _cv_mem, _sun_nlsF );
    if( _check_cv_flag( &_cv_flag, "CVodeSetNonlinearSolverSensStg1", 1 ) ) return false;
    break;
  }
  
  // Initialize integrator memory for quadratures
  if( !_nq ) return true;
  _cv_flag = CVodeQuadSensInit( _cv_mem, MC_CVQUADF__, _Nyq );
  if( _check_cv_flag( &_cv_flag, "CVodeQuadSensInit", 1 ) ) return false;
  
  // Specify whether or not to perform error control on quadrature
  _cv_flag = CVodeSetQuadSensErrCon( _cv_mem, options.QERRS );
  if( _check_cv_flag( &_cv_flag, "CVodeSetQuadSensErrCon", 1 ) ) return false;

  // Specify absolute and relative tolerances for quadratures
  if( options.AUTOTOLS ){
    _cv_flag = CVodeQuadSensEEtolerances( _cv_mem );
    if( _check_cv_flag( &_cv_flag, "CVodeQuadSensEEtolerances", 1 ) ) return false;
  }
  else{
    std::vector<sunrealtype> ATOLS( _nsen, options.ATOLS );
    _cv_flag = CVodeQuadSensSStolerances( _cv_mem, options.RTOLS, ATOLS.data() );
    if( _check_cv_flag( &_cv_flag, "CVodeQuadSensSStolerances", 1 ) ) return false;
  }

  return true;
}

inline
bool
ODESLVS_CVODES::_CC_CVODES_ASA
( unsigned const ifct, int const indexB )
{

  // std::cout << "Calling CVodeReInitB @t=" << _t << std::endl;
  // Reinitialize CVodeS memory block for current time _t and adjoint _Ny
  _cv_flag = CVodeReInitB( _cv_mem, indexB, _t, _Ny[ifct] );
  if( _check_cv_flag( &_cv_flag, "CVodeReInitB", 1 ) ) return false;

#if defined( CRONOS__WITH_KLU )
  switch( options.LINSOL ){
    case Options::SPARSE:
      // Function SUNLinSol_KLUReInit(SUNLinearSolver S, SUNMatrix A, sunindextype nnz, int reinit_type)
      // needed to reinitialize memory and flag for a new factorization (symbolic and numeric) to be conducted at
      // the next solver setup call. This routine is useful in the cases where the number of nonzeroes has changed or if
      // the structure of the linear system has changed which would require a new symbolic (and numeric factorization).
      // std::cout << "Calling SUNLinSol_KLUReInit" << std::endl;
      _cv_flag = SUNLinSol_KLUReInit( _sun_lsB[ifct], _sun_matB[ifct], _nnzjac, 2 );
      if( _cv_flag ) return false;
      break;
    default:
      break;
  }
#endif

  // Reinitialize CVodeS memory block for current adjoint quarature _Nyq
  if( !_nsen ) return true;
  _cv_flag = CVodeQuadReInitB( _cv_mem, indexB, _Nyq[ifct] );
  if( _check_cv_flag( &_cv_flag, "CVodeQuadReInitB", 1 ) ) return false;

  return true;
}

inline
bool
ODESLVS_CVODES::_CC_CVODES_FSA
()
{
  // Reinitialize CVodeS memory block for current sensitivity _Ny
  _cv_flag = CVodeSensReInit( _cv_mem, options.FSACORR, _Ny );
  if( _check_cv_flag( &_cv_flag, "CVodeSensReInit", 1 ) ) return false;

  return true;
}

inline
bool
ODESLVS_CVODES::_CC_CVODES_QUAD
()
{
  // Reinitialize CVode memory block for current sensitivity quarature _Nyq
  if( !_nq ) return true;
  _cv_flag = CVodeQuadSensReInit( _cv_mem, _Nyq );
  if( _check_cv_flag( &_cv_flag, "CVodeQuadSensReInit", 1 ) ) return false;

  return true;
}

inline
void
ODESLVS_CVODES::_END_SEN
( Stats& stats )
{
  // Unset constants - only if states are not stored for adjoints
  _END_D_STA();
  
  // Get final CPU time
  _final_stats( stats );
}

inline
bool
ODESLVS_CVODES::_REINI_SEN
()
{
  // reset at time stages
  _xpk.clear(); _xpk.reserve(_ns);
  _lk.clear();  _lk.reserve(_ns);
  _qpk.clear(); _qpk.reserve(_ns);
  _fp.clear();  _fp.reserve(_nf);

  return true;
}

inline
bool
ODESLVS_CVODES::_INI_ASA
( double const* p )
{
  // Initialize bound propagation
  if( !_INI_D_SEN( p, _nf, _nsen ) || !_REINI_SEN() )
    return false;

  // Set SUNDIALS adjoint/quadrature arrays
  if( _Ny )  N_VDestroyVectorArray( _Ny,  _nvec );
  if( _Nyq ) N_VDestroyVectorArray( _Nyq, _nvec );
  _nvec = _nf;
  _Ny = N_VCloneVectorArray( _nvec, _Nx );
  //_Nyq = N_VCloneVectorArray( _nvec, _Nx );
  _Nyq = N_VNewVectorArray( _nvec, sunctx );
  for( unsigned i=0; i<_nvec; i++ )
    _Nyq[i] = N_VNew_Serial( _nsen, sunctx );

  // Reset result record and statistics
  results_asens.clear();
  results_asens.resize( _nf );
  _init_stats( stats_asens );

  return true;
}

inline
int
ODESLVS_CVODES::MC_CVRHSB__
( sunrealtype t, N_Vector x, N_Vector y, N_Vector ydot, void* user_data )
{
#ifdef CRONOS__BASE_CVODES_CHECK
  //std::cout << "BASE_CVODES::PTR_BASE_CVODES: " << BASE_CVODES::PTR_BASE_CVODES
  //          << "  PTR_CVRHSB: " << PTR_CVRHSB << std::endl;
  assert( BASE_CVODES::PTR_BASE_CVODES != nullptr && PTR_CVRHSB != nullptr);
#endif
  //return (BASE_CVODES::PTR_BASE_CVODES->*PTR_CVRHSB)( t, x, y, ydot, user_data );
  BASE_CVODES::PTR_BASE_CVODES->reregistration();
  auto PTR_BASE_CVODES_ = BASE_CVODES::PTR_BASE_CVODES;
  auto flag = (BASE_CVODES::PTR_BASE_CVODES->*PTR_CVRHSB)( t, x, y, ydot, user_data );
  BASE_CVODES::PTR_BASE_CVODES = PTR_BASE_CVODES_;
  //std::cout << "BASE_CVODES::PTR_BASE_CVODES: RHS Exit " << BASE_CVODES::PTR_BASE_CVODES << std::endl;
  return flag;
}

inline
int
ODESLVS_CVODES::CVRHSB__
( sunrealtype t, N_Vector x, N_Vector y, N_Vector ydot, void* user_data )
{
  _ifct = *static_cast<unsigned*>( user_data );
#ifdef CRONOS__ODESLVS_CVODES_DEBUG
  std::cout << std::scientific << std::setprecision(6) << t;
  for( unsigned i=0; i<NV_LENGTH_S( x ); i++ ) std::cout << "  " << NV_Ith_S( x, i );
  for( unsigned i=0; i<NV_LENGTH_S( y ); i++ ) std::cout << "  " << NV_Ith_S( y, i );
#endif
  bool flag = _RHS_D_SEN( t, NV_DATA_S( x ), NV_DATA_S( y ), NV_DATA_S( ydot ), _ifct );
#ifdef CRONOS__ODESLVS_CVODES_DEBUG
  std::cout << "  " << _ifct;
  for( unsigned i=0; i<NV_LENGTH_S( ydot ); i++ ) std::cout << "  " << NV_Ith_S( ydot, i );
  std::cout << std::endl;
  { int dum; std::cin >> dum; }
#endif
  stats_asens.numRHS++;
  return( flag? 0: -1 );
}

inline
int
ODESLVS_CVODES::MC_CVJACB__
( sunrealtype t, N_Vector x, N_Vector y, N_Vector fy, SUNMatrix Jac,
  void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3 )
{
#ifdef CRONOS__BASE_CVODES_CHECK
  //std::cout << "BASE_CVODES::PTR_BASE_CVODES: " << BASE_CVODES::PTR_BASE_CVODES
  //          << "  PTR_CVJACB: " << PTR_CVJACB << std::endl;
  assert( BASE_CVODES::PTR_BASE_CVODES != nullptr && PTR_CVJACB != nullptr);
#endif
  //return (BASE_CVODES::PTR_BASE_CVODES->*PTR_CVJACB)( t, x, y, fy, Jac, user_data, tmp1, tmp2, tmp3 );
  BASE_CVODES::PTR_BASE_CVODES->reregistration();
  auto PTR_BASE_CVODES_ = BASE_CVODES::PTR_BASE_CVODES;
  auto flag = (BASE_CVODES::PTR_BASE_CVODES->*PTR_CVJACB)( t, x, y, fy, Jac, user_data, tmp1, tmp2, tmp3 );
  BASE_CVODES::PTR_BASE_CVODES = PTR_BASE_CVODES_;
  //std::cout << "BASE_CVODES::PTR_BASE_CVODES: RHS Exit " << BASE_CVODES::PTR_BASE_CVODES << std::endl;
  return flag;
}

inline
int
ODESLVS_CVODES::CVJACB__
( sunrealtype t, N_Vector x, N_Vector y, N_Vector fy, SUNMatrix Jac,
  void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3 )
{
  //std::cerr << "Entering: ODESLV_CVODES::MC_CVJACB__\n";
  bool flag = false;
  switch( options.LINSOL ){
   case Options::DIAG:
   case Options::DENSEDQ:
    flag = false;
    break;
    
   case Options::DENSE:
    flag = _JAC_D_SEN( t, NV_DATA_S( x ), NV_DATA_S( y ), SM_COLS_D(Jac) );
    break;

#if defined( CRONOS__WITH_KLU )
   case Options::SPARSE:
    flag = _JAC_D_SEN( t, NV_DATA_S( x ), NV_DATA_S( y ), SUNSparseMatrix_Data(Jac),
                       SUNSparseMatrix_IndexPointers(Jac), SUNSparseMatrix_IndexValues(Jac) );
    break;
#endif
  }
  stats_asens.numJAC++; // increment JAC counter
  return( flag? 0: -1 );
}

inline
int
ODESLVS_CVODES::MC_CVQUADB__
( sunrealtype t, N_Vector x, N_Vector y, N_Vector qdot, void* user_data )
{
#ifdef CRONOS__BASE_CVODES_CHECK
  //std::cout << "BASE_CVODES::PTR_BASE_CVODES: " << BASE_CVODES::PTR_BASE_CVODES
  //          << "  PTR_CVQUADB: " << PTR_CVQUADB << std::endl;
  assert( BASE_CVODES::PTR_BASE_CVODES != nullptr && PTR_CVQUADB != nullptr);
#endif
  //return (BASE_CVODES::PTR_BASE_CVODES->*PTR_CVQUADB)( t, x, y, qdot, user_data );
  BASE_CVODES::PTR_BASE_CVODES->reregistration();
  auto PTR_BASE_CVODES_ = BASE_CVODES::PTR_BASE_CVODES;
  auto flag = (BASE_CVODES::PTR_BASE_CVODES->*PTR_CVQUADB)( t, x, y, qdot, user_data );
  BASE_CVODES::PTR_BASE_CVODES = PTR_BASE_CVODES_;
  //std::cout << "BASE_CVODES::PTR_BASE_CVODES: RHS Exit " << BASE_CVODES::PTR_BASE_CVODES << std::endl;
  return flag;
}

inline
int
ODESLVS_CVODES::CVQUADB__
( sunrealtype t, N_Vector x, N_Vector y, N_Vector qdot, void* user_data )
{
  _ifct = *static_cast<unsigned*>( user_data );
  // 2026-09-28: load THIS call's (t, x, y) -- it used to integrate the adjoint quadrature (the gradient) at whatever
  // _DVAR held from the last RHS call, e.g. a forward-REPLAY state at another time
  *_Dt = t;  _vec2D( NV_DATA_S( x ), _nx, _Dx );  _vec2D( NV_DATA_S( y ), _ny, _Dy );
  bool flag = _RHS_D_QUAD( _nsen, NV_DATA_S( qdot ), _ifct );
  return( flag? 0: -1 );
}

//! @fn inline typename ODESLVS_CVODES::STATUS ODESLVS_CVODES::solve_asens
//!( std::vector<double> const& p, std::vector<double> const& c=std::vector<double>(), std::ostream& os=std::cout )
//!
//! This function computes a solution to the parametric ODEs with adjoint
//! sensitivity analysis:
//!  - <a>p</a>  [input]  parameter values
//!  - <a>c</a>  [input]  constant values
//!  - <a>os</a> [input/output]  output stream [default: std::cout]
//! .
//! The return value is the status.
inline
typename ODESLVS_CVODES::STATUS
ODESLVS_CVODES::solve_asens
( std::vector<double> const& p, std::vector<double> const& c, std::ostream& os )
{
  registration();
  STATUS flag = _states_ASA( p.data(), c.data(), os );
  unregistration();
  return flag;
}

//! @fn inline typename ODESLVS_CVODES::STATUS ODESLVS_CVODES::solve_asens
//!( double const* p, double const* c=nullptr, std::ostream& os=std::cout )
//!
//! This function computes a solution to the parametric ODEs with adjoint
//! sensitivity analysis:
//!  - <a>p</a>  [input]  parameter values
//!  - <a>c</a>  [input]  constant values
//!  - <a>os</a> [input/output]  output stream [default: std::cout]
//! .
//! The return value is the status.
inline
typename ODESLVS_CVODES::STATUS
ODESLVS_CVODES::solve_asens
( double const* p, double const* c, std::ostream& os )
{
  registration();
  STATUS flag = _states_ASA( p, c, os );
  unregistration();
  return flag;
}

inline
typename ODESLVS_CVODES::STATUS
ODESLVS_CVODES::_states_ASA
( double const* p, double const* c, std::ostream& os )
{
  //std::cerr << "&c: " << c << std::endl;
  //if( c ) std::cerr << "c[0]: " << c[0] << std::endl;

  // Compute state bounds and store intermediate results.  Only the LAST stage is checkpointed: the
  // backward loop re-integrates every earlier stage itself, one clean CVodeF at a time, because a
  // checkpoint set must come from a single uninterrupted forward integration.
  // Per-stage checkpointing is needed only where the sweep would restart.  With no restart the sweep is
  // already one uninterrupted integration, so a single horizon-wide CVodeF is valid and nothing is
  // re-integrated -- the pre-restructure path, at the pre-restructure cost.
  // Two shapes, chosen by the same predicate as the stage restart:
  //  - nothing restarts: arm inside _INI_CVODE and run ONE CVodeF over the horizon (valid because the sweep
  //    is one uninterrupted integration); the backward loop re-integrates nothing.
  //  - some boundary restarts: run the sweep UNARMED on plain CVode (recording _xk only), arm afterwards,
  //    and re-integrate every stage below on checkpoints of its own.  A CVodeF sweep with a mid-way
  //    CVodeReInit would corrupt its checkpoints, and plain CVode on an armed solver poisons it.
  bool const perstage = _sweep_restarts();
  _armOnInit          = !perstage;
  _storeLastStageOnly = false;
  STATUS flag = STATUS::NORMAL;
  flag = ODESLV_CVODES::_states( p, c, !perstage, os );
  if( perstage && flag == STATUS::NORMAL && !_ARM_ASA() ) return STATUS::FATAL;

  // The forward sweep ends with _END_D_STA, which UNSETS the constants' values on _pC.  The backward phase evaluates
  // the same DAG -- terminal conditions, adjoint right-hand sides and quadratures, and the per-stage forward
  // re-integration -- so the constants are set again here, and unset on EVERY exit from this function.
  // (Forward sensitivity never hit this: it integrates during the forward sweep, while the constants are set.)
  struct ConstantsGuard {
    FFVar* pC; unsigned nc;
    ConstantsGuard( FFVar* p_, unsigned n_, double const* c_ ) : pC( p_ ), nc( c_? n_: 0 )
      { for( unsigned i=0; i<nc; ++i ) pC[i].set( c_[i] ); }
    ~ConstantsGuard()
      { for( unsigned i=0; i<nc; ++i ) pC[i].unset(); }
  } constants_guard( _pC, _nc, c );
  if( flag != STATUS::NORMAL ) return flag;

  // Nothing to do if no functions or parameters are defined
  if( !_nf || !_nsen ) return STATUS::NORMAL;

  //std::cerr << "&c: " << c << std::endl;
  //if( c ) std::cerr << "c[0]: " << c[0] << std::endl;

  try{
    // Initialize adjoint integration
    if( !_INI_ASA( p ) ) return STATUS::FATAL;
    _t = _dT[_ns];
    const unsigned NSTEP = options.RESRECORD? options.RESRECORD: 1;

    // Terminal adjoint & quadrature values
    _lk.resize( _ns+1, std::vector<std::vector<double>>( _nf ) );
    _qpk.resize( _ns+1, std::vector<std::vector<double>>( _nf ) );
//    if( lk && !lk[_ns] ) lk[_ns] = new double[(_nx+_np)*_nf];
    _pos_fct = _ns;// ( _vFCT.size()>=_ns? _ns-1:0 );
    for( _ifct=0; _ifct < _nf; _ifct++ ){
      if( !_TC_SET_ASA( _pos_fct, _ifct )
       || !_TC_D_SEN( _t, _xk[_ns].data(), NV_DATA_S(_Ny[_ifct]) )
       || ( _Nyq && _Nyq[_ifct] && !_TC_D_QUAD_ASA( NV_DATA_S(_Nyq[_ifct]) ) ) )
        { _END_SEN( stats_asens ); return STATUS::FATAL; }
      _GET_D_SEN( NV_DATA_S(_Ny[_ifct]), _nsen, _Nyq? NV_DATA_S(_Nyq[_ifct]): nullptr );
      for( unsigned iq=0; iq<_nsen; iq++ )
        _Dfp[iq*_nf+_ifct] = _Dyq[iq];

      // Display / record / return adjoint terminal values
      _lk[_ns].push_back( std::vector<double>( _Dy, _Dy+_ny ) );
      _qpk[_ns].push_back( std::vector<double>( _Dyq, _Dyq+_nsen ) );
      if( options.DISPLAY >= 1 ){
        std::ostringstream ol; ol << " l[" << _ifct << "]";
        if( !_ifct ) _print_interm( _dT[_ns], _nx, _Dy, ol.str(), os );
        else         _print_interm( _nx, _Dy, ol.str(), os );
        std::ostringstream oq; oq << " qp[" << _ifct << "]";
        _print_interm( _nsen, _Dyq, oq.str(), os );
      }
      if( options.RESRECORD )
        results_asens[_ifct].push_back( Results( _t, _nx, _Dy, _nsen, _Dyq ) );
//      for( unsigned iy=0; lk && iy<_ny+_np; iy++ )
//        lk[_ns][_ifct*(_nx+_np)+iy] = iy<_ny? _Dy[iy]: _Dyq[iy-_ny];
    }

    // Initialization of adjoint integration
    for( _ifct=0; _ifct < _nf; _ifct++ )
      if( !_INI_CVODES_ASA( _ifct, _indexB[_ifct], _iusrB[_ifct] ) )
        { _END_SEN( stats_asens ); return STATUS::FATAL;}

    // Integrate adjoint ODEs through each stage using SUNDIALS
    for( _istg=_ns; _istg>0; _istg-- ){

      // Re-integrate this stage forward, on checkpoints of its own.  The sweep left checkpoints for the
      // LAST stage only, so every earlier one is rebuilt here: drop the previous stage's checkpoint list,
      // restart the forward problem from the stored stage-entry state, and run one clean CVodeF over this
      // stage.  CVodeAdjReInit keeps the backward problems (cvodea.c:303 never touches cvB_mem), so
      // _indexB[] stays valid and _CC_CVODES_ASA below still re-anchors lambda per stage.
      if( perstage ){
        if( _istg < _ns ){            // the first stage runs on the freshly-armed memory
          _cv_flag = CVodeAdjReInit( _cv_mem );
          if( _check_cv_flag( &_cv_flag, "CVodeAdjReInit", 1 ) )
            { _END_SEN( stats_asens ); return STATUS::FATAL; }
        }
        _t = _dT[_istg-1];
        _D2vec( _xk[_istg-1].data(), _nx, NV_DATA_S( _Nx ) );
        if(_Nq ) _IC_D_QUAD( NV_DATA_S( _Nq ) );
#ifdef CRONOS__ODESLVS_CVODES_DEBUG
        std::cout << "RESTARTING FORWARD INTEGRATION at t=" << _t << std::endl;
#endif
        _states_stage( _istg-1, _t, _Nx, _Nq, true, true, false, os );
      }

      // Update list of operations in RHSADJ and QUADADJ
      _pos_rhs  = ( _vRHS.size() <=1? 0: _istg-1 );
      _pos_quad = ( _vQUAD.size()<=1? 0: _istg-1 );
      _pos_fct  = _istg;// ( _vFCT.size()>=_ns? _istg: ( _vFCT.size()==1 && _istg==_ns? 1: 0 ) );   
      //_vFCT.size()>=_ns? _istg: 0 );
#ifdef CRONOS__ODESLVS_CVODES_DEBUG
      std::cout << "pos_fct: " << _pos_fct << std::endl;
#endif
      if( !ODESLV_BASE::_RHS_D_SET( _pos_rhs, _pos_quad )
       || !_RHS_SET_ASA( _pos_rhs, _pos_quad, _pos_fct )
       || !_RHS_D_SET( _nf, _nsen ) )
        { _END_SEN( stats_asens ); return STATUS::FATAL; }

      // Propagate adjoints backward to previous stage time
      _t = _dT[_istg];
      const double TSTEP = ( _t - _dT[_istg-1] ) / NSTEP;
      double TSTOP = _t-TSTEP;
      for( unsigned k=0; k<NSTEP; k++, TSTOP-=TSTEP ){
        if( k+1 == NSTEP ) TSTOP = _dT[_istg-1];
        _cv_flag = CVodeB( _cv_mem, TSTOP, CV_NORMAL );
        if( _check_cv_flag( &_cv_flag, "CVodeB", 1 ) )
          { _END_SEN( stats_asens ); return STATUS::FATAL; }

        // intermediate record
        if( options.RESRECORD ){
          for( _ifct=0; _ifct < _nf; _ifct++ ){
            _cv_flag = CVodeGetB( _cv_mem, _indexB[_ifct], &_t, _Ny[_ifct]);
            if( _check_cv_flag( &_cv_flag, "CVodeGetB", 1) )
              { _END_SEN( stats_asens ); return STATUS::FATAL; }
#ifdef CRONOS__ODESLVS_CVODES_DEBUG
            std::cout << "Adjoint #" << _ifct << ": " << _t << std::endl;
            for( unsigned iy=0; iy<_ny; iy++ )
              std::cout << "_Ny" << _ifct << "[" << iy << "] = " << NV_Ith_S(_Ny[_ifct],iy) << std::endl;
#endif
            _cv_flag = CVodeGetQuadB( _cv_mem, _indexB[_ifct], &_t, _Nyq[_ifct]);
            if( _check_cv_flag( &_cv_flag, "CVodeGetQuadB", 1) )
              { _END_SEN( stats_asens ); return STATUS::FATAL; }
#ifdef CRONOS__ODESLVS_CVODES_DEBUG
            for( unsigned ip=0; ip<_nsen; ip++ )
              std::cout << "_Nyq" << _ifct << "[" << ip << "] = " << NV_Ith_S(_Nyq[_ifct],ip) << std::endl;
#endif
            results_asens[_ifct].push_back( Results( _t, _nx, NV_DATA_S(_Ny[_ifct]), _nsen, _Nyq && _Nyq[_ifct]? NV_DATA_S(_Nyq[_ifct]): 0 ) );
          }
        }
      }
      for( _ifct=0; _ifct < _nf; _ifct++ ){
        void *cv_memB = CVodeGetAdjCVodeBmem( _cv_mem, _indexB[_ifct] );
        long int nstpB;
        _cv_flag = CVodeGetNumSteps( cv_memB, &nstpB );
        stats_asens.numSteps += nstpB;
#ifdef CRONOS__ODESLVS_CVODES_DEBUG
        std::cout << "Number of steps for adjoint #" << _ifct << ": " 
                  << nstpB << std::endl;
#endif
      }

      // states/adjoints/quadratures at stage time
#ifdef CRONOS__ODESLVS_CVODES_DEBUG
      for( unsigned ix=0; ix<_nx; ix++ )
        std::cout << "xk[" << _istg-1 << "][" << ix << "] = " << _xk[_istg-1][ix] << std::endl;
#endif
//      if( lk && !lk[_istg-1] ) lk[_istg-1] = new double[(_nx+_np)*_nf];
      for( _ifct=0; _ifct < _nf; _ifct++ ){
        _cv_flag = CVodeGetB( _cv_mem, _indexB[_ifct], &_t, _Ny[_ifct]);
        if( _check_cv_flag( &_cv_flag, "CVodeGetB", 1) )
          { _END_SEN( stats_asens ); return STATUS::FATAL; }
#ifdef CRONOS__ODESLVS_CVODES_DEBUG
        for( unsigned iy=0; iy<_ny; iy++ )
          std::cout << "_Ny" << _ifct << "[" << iy << "] = " << NV_Ith_S(_Ny[_ifct],iy) << std::endl;
#endif
        _cv_flag = CVodeGetQuadB( _cv_mem, _indexB[_ifct], &_t, _Nyq[_ifct]);
        if( _check_cv_flag( &_cv_flag, "CVodeGetQuadB", 1) )
          { _END_SEN( stats_asens ); return STATUS::FATAL; }
#ifdef CRONOS__ODESLVS_CVODES_DEBUG
        for( unsigned ip=0; ip<_nsen; ip++ )
          std::cout << "_Nyq" << _ifct << "[" << ip << "] = " << NV_Ith_S(_Nyq[_ifct],ip) << std::endl;
#endif
        // Add function contribution to adjoint values (discontinuities)
        if( _istg > 1 ){
          _pos_ic  = ( _vIC.size() >=_ns? _istg-1: 0 );
          _pos_fct = _istg-1;//( _vFCT.size()>=_ns? _istg-1: 0 );
          if( ( _pos_fct || _pos_ic )
           && ( !_CC_SET_ASA( _pos_ic, _pos_fct, _ifct )
             || !_CC_D_SEN( _t, _xk[_istg-1].data(), NV_DATA_S(_Ny[_ifct]) )
             || ( _Nyq && _Nyq[_ifct] && !_CC_D_QUAD_ASA( NV_DATA_S(_Nyq[_ifct]) ) ) ) )
            { _END_SEN( stats_asens ); return STATUS::FATAL; }

#ifdef CRONOS__ODESLVS_CVODES_DEBUG
          for( unsigned iy=0; iy<_ny; iy++ )
            std::cout << "_Ny" << _ifct << "[" << iy << "] = " << NV_Ith_S(_Ny[_ifct],iy) << std::endl;
          for( unsigned ip=0; ip<_nsen; ip++ )
            std::cout << "_Nyq" << _ifct << "[" << ip << "] = " << NV_Ith_S(_Nyq[_ifct],ip) << std::endl;
#endif
          _GET_D_SEN( NV_DATA_S(_Ny[_ifct]), _nsen, _Nyq && _Nyq[_ifct]? NV_DATA_S(_Nyq[_ifct]): 0 );
          
          // Reset ODE solver - needed in case of discontinuity
          if( !_CC_CVODES_ASA( _ifct, _indexB[_ifct] ) )
            { _END_SEN( stats_asens ); return STATUS::FATAL; }
        }
        
        // Add initial state contribution to function derivatives 
        else{
          if( !_IC_SET_ASA( _ifct )
           || !_IC_D_SEN( _t, _xk[_istg-1].data(), NV_DATA_S(_Ny[_ifct]) )
           || ( _Nyq && _Nyq[_ifct] && !_IC_D_QUAD_ASA( NV_DATA_S(_Nyq[_ifct]) ) ) )
          { _END_SEN( stats_asens ); return STATUS::FATAL; }
          
#ifdef CRONOS__ODESLVS_CVODES_DEBUG
          for( unsigned iy=0; iy<_ny; iy++ )
            std::cout << "_Ny" << _ifct << "[" << iy << "] = " << NV_Ith_S(_Ny[_ifct],iy) << std::endl;
          for( unsigned ip=0; ip<_nsen; ip++ )
            std::cout << "_Nyq" << _ifct << "[" << ip << "] = " << NV_Ith_S(_Nyq[_ifct],ip) << std::endl;
#endif
          _GET_D_SEN( NV_DATA_S(_Ny[_ifct]), _nsen, _Nyq && _Nyq[_ifct]? NV_DATA_S(_Nyq[_ifct]): 0 );
        }

        // Display / record / return adjoint terminal values
        _lk[_istg-1].push_back( std::vector<double>( _Dy, _Dy+_ny ) );
        _qpk[_istg-1].push_back( std::vector<double>( _Dyq, _Dyq+_nsen ) );
        if( options.DISPLAY >= 1 ){
          std::ostringstream ol; ol << " l[" << _ifct << "]";
          if( !_ifct ) _print_interm( _dT[_istg-1], _nx, _Dy, ol.str(), os );
          else         _print_interm( _nx, _Dy, ol.str(), os );
          std::ostringstream oq; oq << " qp[" << _ifct << "]";
          _print_interm( _nsen, _Dyq, oq.str(), os );
        }
        if( options.RESRECORD )
          results_asens[_ifct].push_back( Results( _t, _nx, _Dy, _nsen, _Dyq ) );
//        for( unsigned iy=0; lk && iy<_ny+_np; iy++ )
//          lk[_istg-1][_ifct*(_nx+_np)+iy] = iy<_ny? _Dy[iy]: _Dyq[iy-_ny];

        // Keep track of function derivatives
        for( unsigned iq=0; iq<_nsen; iq++ ) _Dfp[iq*_nf+_ifct] = _Dyq[iq];
      }
    }

    // Display / return function derivatives
    for( unsigned ip=0; ip<_nsen; ++ip )
      _fp.push_back( std::vector<double>( _Dfp.data()+ip*_nf, _Dfp.data()+(ip+1)*_nf ) );
//    for( unsigned i=0; fp && i<_nf*_np; i++ ) fp[i] = _Dfp[i];
    if( options.DISPLAY >= 1 ){
      for( unsigned iq=0; iq<_nsen; iq++ ){
        std::ostringstream ofp; ofp << " fp[" << iq << "]";
        _print_interm( _nf, _Dfp.data()+iq*_nf, ofp.str(), os );
      }
    }
  }
  catch(...){
    _END_SEN( stats_asens );
    if( options.DISPLAY >= 1 ) _print_stats( stats_asens, os );
    return STATUS::FAILURE;
  }

  _END_SEN( stats_asens );
  if( options.DISPLAY >= 1 ) _print_stats( stats_asens, os );
  return STATUS::NORMAL;
}

inline
bool
ODESLVS_CVODES::_INI_FSA
( double const* p )
{
  // Initialize bound propagation
  if( !_INI_D_SEN( p, _nsen, _nq ) || !_REINI_SEN() )
    return false;

  // Set SUNDIALS sensitivity/quadrature arrays
  if( _Ny )   N_VDestroyVectorArray( _Ny,  _nvec );
  if( _Nyq )  N_VDestroyVectorArray( _Nyq, _nvec );
  _nvec = _nsen;
  _Ny  = N_VCloneVectorArray( _nvec, _Nx );
  _Nyq = _nq? N_VCloneVectorArray( _nvec, _Nq ): nullptr;

  // Reset result record and statistics
  results_fsens.clear();
  results_fsens.resize( _nsen );
  _init_stats( stats_fsens );

  return true;
}

inline
int
ODESLVS_CVODES::MC_CVRHSF__
( int Ns, sunrealtype t, N_Vector x, N_Vector xdot, int is, N_Vector y,
  N_Vector ydot, void* user_data, N_Vector tmp1, N_Vector tmp2 )
{
#ifdef CRONOS__BASE_CVODES_CHECK
  //std::cout << "BASE_CVODES::PTR_BASE_CVODES: " << BASE_CVODES::PTR_BASE_CVODES
  //          << "  PTR_CVRHSF: " << PTR_CVRHSF << std::endl;
  assert( BASE_CVODES::PTR_BASE_CVODES != nullptr && PTR_CVRHSF != nullptr);
#endif
  //return (BASE_CVODES::PTR_BASE_CVODES->*PTR_CVRHSF)( Ns, t, x, xdot, is, y, ydot, user_data, tmp1, tmp2 );
  BASE_CVODES::PTR_BASE_CVODES->reregistration();
  auto PTR_BASE_CVODES_ = BASE_CVODES::PTR_BASE_CVODES;
  auto flag = (BASE_CVODES::PTR_BASE_CVODES->*PTR_CVRHSF)( Ns, t, x, xdot, is, y, ydot, user_data, tmp1, tmp2 );
  BASE_CVODES::PTR_BASE_CVODES = PTR_BASE_CVODES_;
  //std::cout << "BASE_CVODES::PTR_BASE_CVODES: RHS Exit " << BASE_CVODES::PTR_BASE_CVODES << std::endl;
  return flag;
}

inline
int
ODESLVS_CVODES::CVRHSF__
( int Ns, sunrealtype t, N_Vector x, N_Vector xdot, int is, N_Vector y,
  N_Vector ydot, void* user_data, N_Vector tmp1, N_Vector tmp2 )
{
#ifdef CRONOS__ODESLVS_CVODES_DEBUG
  std::cout << "@t=" << t << "\nx:\n";
  for( unsigned i=0; i<NV_LENGTH_S( x ); i++ ) std::cout << NV_Ith_S( x, i ) << std::endl;
  std::cout << "y[" << is << "]:\n";
  for( unsigned i=0; i<NV_LENGTH_S( y ); i++ ) std::cout << NV_Ith_S( y, i ) << std::endl;
#endif
  bool flag = _RHS_D_SEN( t, NV_DATA_S( x ), NV_DATA_S( y ), NV_DATA_S( ydot ), is );
#ifdef CRONOS__ODESLVS_CVODES_DEBUG
  std::cout << "ydot[" << is << "]:\n";
  for( unsigned i=0; i<NV_LENGTH_S( ydot ); i++ ) std::cout << NV_Ith_S( ydot, i ) << std::endl;
  { int dum; std::cin >> dum; }
#endif
  stats_fsens.numRHS++;
  return( flag? 0: -1 );
}

inline
int
ODESLVS_CVODES::MC_CVQUADF__
( int Ns, sunrealtype t, N_Vector x, N_Vector* y, N_Vector qdot, N_Vector* qSdot, 
  void* user_data, N_Vector tmp1, N_Vector tmp2 )
{
#ifdef CRONOS__BASE_CVODES_CHECK
  //std::cout << "BASE_CVODES::PTR_BASE_CVODES: " << BASE_CVODES::PTR_BASE_CVODES
  //          << "  PTR_CVQUADF: " << PTR_CVQUADF << std::endl;
  assert( BASE_CVODES::PTR_BASE_CVODES != nullptr && PTR_CVQUADF != nullptr);
#endif
  //return (BASE_CVODES::PTR_BASE_CVODES->*PTR_CVQUADF)( Ns, t, x, y, qdot, qSdot, user_data, tmp1, tmp2 );
  BASE_CVODES::PTR_BASE_CVODES->reregistration();
  auto PTR_BASE_CVODES_ = BASE_CVODES::PTR_BASE_CVODES;
  auto flag = (BASE_CVODES::PTR_BASE_CVODES->*PTR_CVQUADF)( Ns, t, x, y, qdot, qSdot, user_data, tmp1, tmp2 );
  BASE_CVODES::PTR_BASE_CVODES = PTR_BASE_CVODES_;
  //std::cout << "BASE_CVODES::PTR_BASE_CVODES: RHS Exit " << BASE_CVODES::PTR_BASE_CVODES << std::endl;
  return flag;
}

inline
int
ODESLVS_CVODES::CVQUADF__
( int Ns, sunrealtype t, N_Vector x, N_Vector* y, N_Vector qdot, N_Vector* qSdot, 
  void* user_data, N_Vector tmp1, N_Vector tmp2 )
{
#ifdef CRONOS__ODESLVS_CVODES_DEBUG
  std::cout << "@t=" << t << "\nx:\n";
  for( unsigned i=0; i<NV_LENGTH_S( x ); i++ ) std::cout << NV_Ith_S( x, i ) << std::endl;
  std::cout << "qdot:\n";
  for( unsigned i=0; i<NV_LENGTH_S( qdot ); i++ ) std::cout << NV_Ith_S( qdot, i ) << std::endl;
#endif
  bool flag = true;
  *_Dt = t;   // 2026-09-28: _GET_D_SEN loads x and y[is] but not t
  for( int is=0; is<Ns && flag; is++ ){
    _GET_D_SEN( NV_DATA_S(x), NV_DATA_S(y[is]), (sunrealtype*)nullptr, 0, (sunrealtype*)nullptr );
#ifdef CRONOS__ODESLVS_CVODES_DEBUG
    std::cout << "y:\n";
    for( unsigned i=0; i<NV_LENGTH_S( y[is] ); i++ ) std::cout << NV_Ith_S( y[is], i ) << std::endl;
#endif
    flag = _RHS_D_QUAD( _nq, NV_DATA_S( qSdot[is] ), is );
#ifdef CRONOS__ODESLVS_CVODES_DEBUG
    std::cout << "qSdot:\n";
    for( unsigned i=0; i<NV_LENGTH_S( qSdot[is] ); i++ ) std::cout << NV_Ith_S( qSdot[is], i ) << std::endl;
    { int dum; std::cin >> dum; }
#endif
  }
  return( flag? 0: -1 );
}

//! @fn inline typename ODESLVS_CVODES::STATUS ODESLVS_CVODES::solve_fsens(
//! std::vector<double> const& p, std::vector<double> const& c=std::vector<double>(), std::ostream& os=std::cout )
//!
//! This function computes a solution to the parametric ODEs with forward
//! sensitivity analysis:
//!  - <a>p</a>   [input]  parameter values
//!  - <a>c</a>   [input]  constant values
//!  - <a>os</a>  [input/output]  output stream [default: std::cout]
//! .
//! The return value is the status.
inline
typename ODESLVS_CVODES::STATUS
ODESLVS_CVODES::solve_fsens
( std::vector<double> const& p, std::vector<double> const& c, std::ostream& os )
{
  registration();
  STATUS flag = _states_FSA( p.data(), c.data(), os );
  unregistration();
  return flag;
}

//! @fn inline typename ODESLVS_CVODES::STATUS ODESLVS_CVODES::solve_fsens(
//! double const* p, double const* c=nullptr, std::ostream& os=std::cout )
//!
//! This function computes a solution to the parametric ODEs with forward
//! sensitivity analysis:
//!  - <a>p</a>   [input]  parameter values
//!  - <a>c</a>   [input]  constant values
//!  - <a>os</a>  [input/output]  output stream [default: std::cout]
//! .
//! The return value is the status.
inline
typename ODESLVS_CVODES::STATUS
ODESLVS_CVODES::solve_fsens
( double const* p, double const* c, std::ostream& os )
{
  registration();
  STATUS flag = _states_FSA( p, c, os );
  unregistration();
  return flag;
}

inline
typename ODESLVS_CVODES::STATUS
ODESLVS_CVODES::_states_FSA
( double const* p, double const* c, std::ostream& os )
{
  // Check arguments
  if( !_nsen )
    return _states( nullptr, c, false, os );
  else if( !p )
    return STATUS::FATAL;

  try{
    // Initialize trajectory integration
    if( !ODESLV_CVODES::_INI_STA( p, c ) 
     || !_INI_FSA( p ) ) return STATUS::FATAL;
    _t = _dT[0];
    const unsigned NSTEP = options.RESRECORD? options.RESRECORD: 1;

    // Initial state/quadrature values
    if( !ODESLV_BASE::_IC_D_SET()
     || !ODESLV_BASE::_IC_D_STA( _t, NV_DATA_S( _Nx ) )
     || ( _Nq && !ODESLV_BASE::_IC_D_QUAD( NV_DATA_S( _Nq ) ) ) )
      { _END_STA(); _END_SEN( stats_fsens ); return STATUS::FATAL; }
    ODESLV_BASE::_GET_D_STA( NV_DATA_S(_Nx), _nq && _Nq? NV_DATA_S(_Nq): nullptr );

    // Add initial function terms
    _pos_fct = 0;
    if( !ODESLV_BASE::_FCT_D_STA( _pos_fct, _t ) )
      { _END_STA(); _END_SEN( stats_fsens ); return STATUS::FATAL; }

    // Display / record / return initial results
    _xk.push_back( std::vector<double>( _Dx, _Dx+_nx ) );
    if( _nq ) _qk.push_back( std::vector<double>( _Dq, _Dq+_nq ) );
    if( options.DISPLAY >= 1 ){
      _print_interm( _t, _nx, _Dx, " x", os );
      _print_interm( _nq, _Dq, " q", os );
    }
//    if( options.RESRECORD )
//      results_solve.push_back( Results( _t, _nx, NV_DATA_S(_Nx), _nq, _nq? NV_DATA_S(_Nq): nullptr ) );

    // Initial state/quadrature sensitivities
    _xpk.push_back( std::vector<std::vector<double>>( _nsen ) );
    if( _nq ) _qpk.push_back( std::vector<std::vector<double>>( _nsen ) );
    for( _isen=0; _isen<_nsen; _isen++ ){
      if( !_IC_SET_FSA( _isen )
       || !_IC_D_SEN( _t, NV_DATA_S(_Ny[_isen]) )
       || ( _Nyq && _Nyq[_isen] && !ODESLV_BASE::_IC_D_QUAD( NV_DATA_S(_Nyq[_isen]) ) ) ) 
        { _END_STA(); _END_SEN( stats_fsens ); return STATUS::FATAL; }
      _GET_D_SEN( NV_DATA_S(_Ny[_isen]), _nq, _nq && _Nyq? NV_DATA_S(_Nyq[_isen]): nullptr );

      // Display / record / return initial results
      _xpk[0].push_back( std::vector<double>( _Dy, _Dy+_nx ) );
      if( _nq ) _qpk[0].push_back( std::vector<double>( _Dyq, _Dyq+_nq ) );
      if( options.DISPLAY >= 1 ){
        std::ostringstream oxp; oxp << " xp[" << _isen << "]";
        _print_interm( _nx, _Dy, oxp.str(), os );
        std::ostringstream oqp; oqp << " qp[" << _isen << "]";
        _print_interm( _nq, _Dyq, oqp.str(), os );
      }
//      if( options.RESRECORD )
//        results_fsens[_isen].push_back( Results( _t, _nx, NV_DATA_S(_Ny[_isen]), _nq, _nq? NV_DATA_S(_Nyq[_isen]):nullptr ) );
//      for( unsigned ix=0; xpk && ix<_nx+_nq; ix++ )
//        xpk[0][(_nx+_nq)*_isen+ix] = ix<_nx? _Dy[ix]: _Dyq[ix-_nx];

      // Add initial function derivative terms
      if( !_FCT_D_SEN( _pos_fct, _isen, _t ) )
          { _END_STA(); _END_SEN( stats_fsens ); return STATUS::FATAL; }
    }

    // Integrate ODEs through each stage using SUNDIALS
    if( !ODESLV_CVODES::_INI_CVODE()
     || !_INI_CVODES_FSA() )
      { _END_STA(); _END_SEN( stats_fsens ); return STATUS::FATAL; }

    for( _istg=0; _istg<_ns; _istg++ ){
    
      // Account for state/sensitivity discontinuities (if any) at stage times
      // and solver reinitialization (if applicable)
      _pos_ic = ( _vIC.size()>=_ns? _istg:0 );
      if( _pos_ic
       && ( !ODESLV_BASE::_CC_D_SET( _pos_ic )
         || !ODESLV_BASE::_CC_D_STA( _t, NV_DATA_S( _Nx ) )
         || !ODESLV_CVODES::_CC_CVODE_STA() ) )
        { _END_STA(); _END_SEN( stats_fsens ); return STATUS::FAILURE; }
      if( _istg 
       //&& !ODESLV_CVODES::_CC_CVODE_QUAD() )
       && ( ( _Nq && !ODESLV_BASE::_IC_D_QUAD( NV_DATA_S( _Nq ) ) ) // quadrature reinitialization
         || !ODESLV_CVODES::_CC_CVODE_QUAD() ) )
        { _END_STA(); _END_SEN( stats_fsens ); return STATUS::FAILURE; }
      if( options.RESRECORD )
        results_solve.push_back( Results( _t, _nx, NV_DATA_S(_Nx), _nq, _nq? NV_DATA_S(_Nq): nullptr ) );

      for( _isen=0; _isen<_nsen; _isen++ ){
        if( _pos_ic
         && ( !_CC_SET_FSA( _pos_ic, _isen )
           || !_CC_D_SEN( _t, NV_DATA_S( _Nx ), NV_DATA_S(_Ny[_isen]) )
           || ( _isen==_nsen-1 && !_CC_CVODES_FSA() ) ) )
            { _END_STA(); _END_SEN( stats_fsens ); return STATUS::FATAL; }
        if( _istg
         //&& ( _isen==_np-1 && !_CC_CVODES_QUAD() ) )
         && ( ( _Nyq && _Nyq[_isen] && !ODESLV_BASE::_IC_D_QUAD( NV_DATA_S(_Nyq[_isen]) ) ) //quadrature sensitivity reinitialization
           || ( _isen==_nsen-1 && !_CC_CVODES_QUAD() ) ) )
            { _END_STA(); _END_SEN( stats_fsens ); return STATUS::FATAL; }
#ifdef CRONOS__ODESLVS_CVODES_DEBUG
        for( unsigned iy=0; iy<NV_LENGTH_S(_Ny[_isen]); iy++ )
          std::cout << "_Ny" << _isen << "[" << iy << "] = " << NV_Ith_S(_Ny[_isen],iy) << std::endl;
        for( unsigned iy=0; _nq && iy<NV_LENGTH_S(_Nyq[_isen]); iy++ )
          std::cout << "_Nyq" << _isen << "[" << iy << "] = " << NV_Ith_S(_Nyq[_isen],iy) << std::endl;
#endif
        if( options.RESRECORD )
          results_fsens[_isen].push_back( Results( _t, _nx, NV_DATA_S(_Ny[_isen]), _nq, _nq? NV_DATA_S(_Nyq[_isen]): nullptr ) );
      }

      // update list of operations in RHS, JAC, QUAD, RHSFSA and QUADFSA
      _pos_rhs  = ( _vRHS.size()<=1?  0: _istg );
      _pos_quad = ( _vQUAD.size()<=1? 0: _istg );
      if( (!_istg || _pos_rhs || _pos_quad)
        && ( !ODESLV_BASE::_RHS_D_SET( _pos_rhs, _pos_quad )
          || !_RHS_SET_FSA( _pos_rhs, _pos_quad )
          || !_RHS_D_SET( _nsen, _nq ) ) )
        { _END_STA(); _END_SEN( stats_fsens ); return STATUS::FATAL; }

      // integrate till end of time stage
      _cv_flag = CVodeSetStopTime( _cv_mem, _dT[_istg+1] );
      if( _check_cv_flag(&_cv_flag, "CVodeSetStopTime", 1) )
        { _END_STA(); _END_SEN( stats_fsens ); return STATUS::FATAL; }

      const double TSTEP = ( _dT[_istg+1] - _t ) / NSTEP;
      double TSTOP = _t+TSTEP;
      for( unsigned k=0; k<NSTEP; k++, TSTOP+=TSTEP ){
        if( k+1 == NSTEP ) TSTOP = _dT[_istg+1];
        _cv_flag = CVode( _cv_mem, TSTOP, _Nx, &_t, CV_NORMAL );
        if( _check_cv_flag( &_cv_flag, "CVode", 1 ) )
         //|| (options.NMAX && stats_fsens.numSteps > options.NMAX) )
          throw Exceptions( Exceptions::INTERN );

        // intermediate record
        if( options.RESRECORD ){
          if( _nq ){
            _cv_flag = CVodeGetQuad( _cv_mem, &_t, _Nq );
            if( _check_cv_flag(&_cv_flag, "CVodeGetQuad", 1) )
              { _END_STA(); _END_SEN( stats_fsens ); return STATUS::FATAL; }
          }
          results_solve.push_back( Results( _t, _nx, NV_DATA_S(_Nx), _nq, _nq? NV_DATA_S(_Nq): nullptr ) );
          for( _isen=0; _isen<_nsen; _isen++ ){
            _cv_flag = CVodeGetSens1(_cv_mem, &_t, _isen, _Ny[_isen] );
            if( _check_cv_flag( &_cv_flag, "CVodeGetSens", 1) )
             { _END_STA(); _END_SEN( stats_fsens ); return STATUS::FATAL; }
            if( _nq ){
              _cv_flag = CVodeGetQuadSens1(_cv_mem, &_t, _isen, _Nyq[_isen]);
              if( _check_cv_flag( &_cv_flag, "CVodeGetQuadSens", 1) )
                { _END_STA(); _END_SEN( stats_fsens ); return STATUS::FATAL; }
              //for( unsigned iq=0; iq<_nq; ++iq )
              //  std::cerr << "_Nyq[" << _isen << "][" << iq << "] = " << NV_DATA_S(_Nyq[_isen])[iq] << std::endl;
            }
            results_fsens[_isen].push_back( Results( _t, _nx, NV_DATA_S(_Ny[_isen]), _nq, _nq? NV_DATA_S(_Nyq[_isen]): nullptr ) );
          }
        }
      }

      // Intermediate states and quadratures
      if( _nq ){
        _cv_flag = CVodeGetQuad( _cv_mem, &_t, _Nq );
        if( _check_cv_flag(&_cv_flag, "CVodeGetQuad", 1) )
          { _END_STA(); _END_SEN( stats_fsens ); return STATUS::FATAL; }
      }
      ODESLV_BASE::_GET_D_STA( NV_DATA_S(_Nx), _nq && _Nq? NV_DATA_S(_Nq): nullptr );

      // Display / return stage results
      _xk.push_back( std::vector<double>( _Dx, _Dx+_nx ) );
      if( _nq ) _qk.push_back( std::vector<double>( _Dq, _Dq+_nq ) );
      if( options.DISPLAY >= 1 ){
        _print_interm( _t, _nx, _Dx, " x", os );
        _print_interm( _nq, _Dq, " q", os );
      }

      // Add intermediate function terms
      _pos_fct = _istg+1;//( _vFCT.size()>=_ns? _istg:0 );
//      if( (_vFCT.size()>=_ns || _istg==_ns-1)
//       && !ODESLV_BASE::_FCT_D_STA( _pos_fct, _t ) )
      if( !ODESLV_BASE::_FCT_D_STA( _pos_fct, _t ) )
        { _END_STA(); _END_SEN( stats_fsens ); return STATUS::FATAL; }

      // Intermediate state and quadrature sensitivities
      _xpk.push_back( std::vector<std::vector<double>>( _nsen ) );
      if( _nq ) _qpk.push_back( std::vector<std::vector<double>>( _nsen ) );
      for( _isen=0; _isen<_nsen; _isen++ ){
        _cv_flag = CVodeGetSens1(_cv_mem, &_t, _isen, _Ny[_isen] );
        if( _check_cv_flag( &_cv_flag, "CVodeGetSens", 1) )
         { _END_STA(); _END_SEN( stats_fsens ); return STATUS::FATAL; }
        if( _nq ){
          _cv_flag = CVodeGetQuadSens1(_cv_mem, &_t, _isen, _Nyq[_isen]);
          if( _check_cv_flag( &_cv_flag, "CVodeGetQuadSens", 1) )
            { _END_STA(); _END_SEN( stats_fsens ); return STATUS::FATAL; }
          //for( unsigned iq=0; iq<_nq; ++iq )
          //  std::cerr << "_Nyq[" << _isen << "][" << iq << "] = " << NV_DATA_S(_Nyq[_isen])[iq] << std::endl;
        }
        _GET_D_SEN( NV_DATA_S(_Nx), NV_DATA_S(_Ny[_isen]), _nq && _Nq? NV_DATA_S(_Nq): nullptr,
                    _nq, _nq && _Nyq[_isen]? NV_DATA_S(_Nyq[_isen]): nullptr );

        // Display / return stage results
        _xpk[_istg+1].push_back( std::vector<double>( _Dy, _Dy+_nx ) );
        if( _nq ) _qpk[_istg+1].push_back( std::vector<double>( _Dyq, _Dyq+_nq ) );
        if( options.DISPLAY >= 1 ){
          std::ostringstream oxp; oxp << " xp[" << _isen << "]";
          _print_interm( _nx, _Dy, oxp.str(), os );
          std::ostringstream oqp; oqp << " qp[" << _isen << "]";
          _print_interm( _nq, _Dyq, oqp.str(), os );
        }
//        for( unsigned ix=0; xpk && ix<_nx+_nq; ix++ )
//          xpk[_istg+1][(_nx+_nq)*_isen+ix] = ix<_nx? _Dy[ix]: _Dyq[ix-_nx];

        // Add intermediate function derivative terms
//        if( (_vFCT.size()>=_ns || _istg==_ns-1)
//         && !_FCT_D_SEN( _pos_fct, _isen, _t ) )
        if( !_FCT_D_SEN( _pos_fct, _isen, _t ) )
          { _END_STA(); _END_SEN( stats_fsens ); return STATUS::FATAL; }
      }
    }

    // Display / return function values and derivatives
    _f = _Df;
    for( unsigned ip=0; ip<_nsen; ++ip )
      _fp.push_back( std::vector<double>( _Dfp.data()+ip*_nf, _Dfp.data()+(ip+1)*_nf ) );
//    for( unsigned i=0; f && i<_nf; i++ ) f[i] = _Df[i];
//    for( unsigned i=0; fp && i<_nf*_np; i++ ) fp[i] = _Dfp[i];
    if( options.DISPLAY >= 1 ){
      _print_interm( _nf, _Df.data(), " f", os );
      for( unsigned iq=0; iq<_nsen; iq++ ){
        std::ostringstream ofp; ofp << " fp[" << iq << "]";
        _print_interm( _nf, _Dfp.data()+iq*_nf, ofp.str(), os );
      }
    }
  }
  catch(...){
    _END_STA(); _END_SEN( stats_fsens );
    long int nstp;
    _cv_flag = CVodeGetNumSteps( _cv_mem, &nstp );
    stats_solve.numSteps += nstp;
    stats_fsens.numSteps += nstp;
    if( options.DISPLAY >= 1 ) _print_stats( stats_fsens, os );
    return STATUS::FAILURE;
  }

  long int nstp;
  _cv_flag = CVodeGetNumSteps( _cv_mem, &nstp );
  stats_solve.numSteps += nstp;
  stats_fsens.numSteps += nstp;
#ifdef CRONOS__ODESLVS_CVODES_DEBUG
  std::cout << "number of steps: " << nstp << std::endl;
#endif

  _END_STA(); _END_SEN( stats_fsens );
  if( options.DISPLAY >= 1 ) _print_stats( stats_fsens, os );
  return STATUS::NORMAL;
}

} // end namescape mc

#endif

