// Copyright (C) Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

#ifndef MC__FFMODEL_HPP
#define MC__FFMODEL_HPP

// The model layer of the collocation stack (rev199: the principal symbol moved in, F2a step 1): the declared model and its declaration API.
// ocbase.hpp is included for the DAG operators the equations are written with (FFPartial,
// FFIntegral, FFEval) and the OCVar arithmetic they dispatch on; FFModel itself holds no
// discretisation state -- no mesh, no basis, no FFDom member, no collocation call.
#include <mutex>
#include "ocbase.hpp"
#if defined(__unix__) || defined(__APPLE__)
# include <dlfcn.h>     // dlsym: the run-time thread setters (_thread_hooks); none on Windows (the cap then reports it)
#endif
#include <typeinfo>   // rev205: t_ThreadCap resolves the OpenMP/BLAS thread setters at run time
#include <cassert>
#include <sstream>
#include <cmath>
#include <cfloat>
#include <ostream>
#include <functional>
#include <map>
#include <optional>
#include <set>
#include <string>
#include <memory>
#include <vector>

namespace mc
{

//! @brief Model-level base class: the declared model and its declaration API.
////////////////////////////////////////////////////////////////////////
//! mc::FFModel holds a factorable model as declared by the user -- the DAG, constants,
//! domains, inputs, states, equations with their options, output functions, the evolution
//! domain, and classification reference values and functions -- and,
//! once a derived solver has run setup(), the working-DAG copies of these collections.
//! It performs no discretisation and caches nothing solver-side: model mutations call
//! _on_model_changed(), through which a derived solver drops what the change makes stale.
//!
//! FFModel owns its working DAG: it creates it when a derived solver sets up or deep-copies the
//! model, and deletes it on release, reset of the user DAG, or destruction.  It is not copyable by
//! itself (a copy would share the working DAG); derived solvers provide deep copies.
//!
//! WHAT IS SETTLED, and what a reader can rely on:
//!  - the model sets up STANDALONE -- setup() runs, reduces order and classifies with no solver, no mesh and
//!    no collocation present.  OCFE_ffmodel.cpp links ffmodel.hpp alone and proves it;
//!  - setup() is NOT virtual: it is FFModel's skeleton and a derived solver contributes through named hooks,
//!    so a solver cannot skip the model work;
//!  - set_model() copies another model's BARE DECLARATIONS, resetting everything derived.
//!
//! STILL TRANSITIONAL (see the current WORKPLAN):
//!  - domains are declared as FFDom objects, which carry the mesh, and FFModel derives from OCBase.  This is
//!    DELIBERATE and is not scheduled for replacement: a separate mesh-free declaration would duplicate FFDom
//!    for no gain.  Mesh-freeness is a property of the model's ANALYSES -- setup(), order reduction and the
//!    classification run without collocation (see WHAT IS SETTLED above) -- not of the declaration type;
//!  - the model has no preflight of its own: _check_consistency stayed solver-side because it also validates
//!    the marching grid, so a standalone FFModel is not yet checked before setup().
//!
//! PROVENANCE: extracted from OCFESLV -- the lifecycle, the analyses (symbol, classification, order reduction),
//! the setup() skeleton, the mesh-free classification reference, set_model() and add_input(), and the claim
//! classifier.
////////////////////////////////////////////////////////////////////////
class FFModel
////////////////////////////////////////////////////////////////////////
{
public:
  //! @brief Identifies this ffmodel.hpp.  The solver header carries its own OCFESLV::HEADER_ID, and a binary can
  //! mix the two (a sweep has run one revision's model layer under another's solver), so both are printed.
  static constexpr char const* HEADER_ID
    = "ffmodel  rev362  2026-10-09";

  //! @brief The revision of this ffmodel.hpp (HEADER_ID), e.g. for a bug report; the model report does not print it.
  static char const* revision() { return HEADER_ID; }

  //! @brief Kind of model change notified to a derived solver by _on_model_changed().
  enum class ModelChange
  {
    DERIVATIVES = 0  //!< declared model or reference values changed: derivative/Jacobian caches are stale
  };

  //! @brief Empty model with no user DAG.
  FFModel()
    : _pClassification( new t_Classify ), _classification( *_pClassification )
    {}

  //! @brief Empty model declared on the user DAG @p dag (see set()).
  explicit FFModel
    ( FFGraph* dag )
    : _pClassification( new t_Classify ), _classification( *_pClassification )
    { set( dag ); }

  //! @brief Not copyable: a copy would share the working DAG.  Derived solvers provide deep copies.
  FFModel( FFModel const& ) = delete;
  FFModel& operator=( FFModel const& ) = delete;

  //! @brief Destructor: clears the declared model and deletes the working DAG, if owned.
  virtual ~FFModel()
    {
      _usr.clear();
      _release_working_dag();
    }

  //! @brief Attach the user DAG @p dag, clearing the declared model, the working model and the working DAG.
  //! Derived solvers override it to clear their own state as well, and call this version.
  virtual void set
    ( FFGraph* dag )
    {
      _reset_model();
      _usr.clear();
      _usr._dagUsr = dag;
      _release_working_dag();
    }

  //! @brief Compatibility wrapper for set().
  void set_dag
    ( FFGraph* dag )
    { set( dag ); }

protected:
public:
  //! @brief Replace this model by a copy of @p src's DECLARED model, before any transformation: constants and
  //! values, domains, inputs with their declared collocation and continuity, states, equations, outputs, the
  //! evolution domain and the classification references.  Nothing derived is copied -- reduction, classification
  //! and the working DAG are reset, because none of them is part of the model as declared.
  //! @p src's declared roots are imported into THIS model's user DAG, so the two models are independent
  //! afterwards; set() the user DAG first.
  //! @return false if no user DAG has been set here, or if the import fails
  bool set_model
    ( FFModel const& src );


protected:
  //! @brief C0 -- the continuity-CLAIM classification of a state in a direction (structural, no reference point).
  //! PROTECTED, not public.  The claim classification is an internal analysis that C1/C2 consume; a
  //! CALLER has no use for it, and exposing it invited the drivers to depend on a shape still under design.
  //! claim_report() stays public: printing a report is a legitimate thing for a driver to ask for.
  enum class ClaimTag { NATURAL, IMPLIED, CLOSURE, NONE };

  //! @brief One (state, direction) entry of the claim classification.
  struct t_Claim
  {
    ClaimTag           tag   = ClaimTag::NONE;
    int                order = -1;      //!< continuity order: 1 = C1 (flux-claimed), 0 = C0 (value-claimed), -1 = none
    bool               algebraic = false;
    size_t             via   = size_t(-1);  //!< row index in the working equation list justifying the tag
    char               coeff = '-';     //!< IMPLIED/CLOSURE: FFInv type of the relation in s (L, S, N, U)
    std::vector<FFVar> rests_on;        //!< IMPLIED/CLOSURE: the states the determining relation rests on
  };

  //! @brief state -> direction -> claim
  typedef std::map< FFVar, std::map<FFVar,t_Claim,lt_FFVar>, lt_FFVar >  t_ClaimMap;

protected:
  //! @brief The states a working row contains BARE, i.e. outside every derivative/integral/evaluation node.
  //! Only such a state can be determined POINTWISE by that row: one appearing solely under a derivative is
  //! constrained in the differentiated direction, not at a point.
  std::set<FFVar,lt_FFVar> _row_bare_states
    ( size_t const row )
    const;

  //! @brief FFInv invertibility type of one working row in each STATE it depends on (inputs and constants
  //! unseeded): L constant-coefficient linear, S linear with a coefficient depending on other states, N
  //! separably nonlinear, U undetermined (in particular: entering only through a partial).
  std::map<FFVar::pt_idVar,char> _row_invertibility
    ( size_t const row )
    const;

  //! @brief Compute the claim classification for every block, or one block (block_id >= 0).  Report only:
  //! nothing in the model or a derived solver's plan is changed.  Requires setup().
  t_ClaimMap claim_classification
    ( int const block_id = -1 )
    const;

public:
  //! @brief What the C0 classification says about ONE (state, direction) -- the PUBLIC view of an analysis
  //! whose internals (ClaimTag, t_Claim, t_ClaimMap, claim_classification) are protected.
  //! added because making those protected broke OCFE_claim.cpp, the driver that VALIDATES the
  //! classification.  A test asserting on the analysis is a legitimate consumer; exposing the whole map to it
  //! is not.  This is the minimum a consumer needs, and it is stable even while the internals are in design.
  struct ClaimInfo
  {
    bool  found = false;   //!< false when the state does not depend on that direction, or is not classified
    char  tag   = '?';     //!< 'A' natural, 'I' implied, 'C' closure, 'N' none
    char  coeff = '-';     //!< FFInv letter of the determining relation in the state (L/S/N/U), '-' if n/a
    int   order = -1;      //!< continuity order: 1 = C1 (flux-claimed), 0 = C0 (value-claimed), -1 = none
  };

  //! @brief The claim classification for one (state, direction); @p block_id < 0 searches every block.
  ClaimInfo claim_info
    ( FFVar const& state, FFVar const& dir, int const block_id = -1 )
    const
    {
      ClaimInfo out;
      auto const cm = claim_classification( block_id );
      auto is = cm.find( state );
      if( is == cm.end() ) return out;
      auto id = is->second.find( dir );
      if( id == is->second.end() ) return out;
      out.found = true;
      out.tag   = ( id->second.tag == ClaimTag::NATURAL )? 'A':
                  ( id->second.tag == ClaimTag::IMPLIED )? 'I':
                  ( id->second.tag == ClaimTag::CLOSURE )? 'C': 'N';
      out.coeff = id->second.coeff;
      out.order = id->second.order;
      return out;
    }

  //! @brief Print the claim classification, one line per (state, direction), with its justification.
  void claim_report
    ( std::ostream& os, int const block_id = -1 )
    const;

  //! @brief Set the model up: check it, import it onto a private working DAG, reduce its order and run the
  //! auto-elimination.  NOT VIRTUAL: a derived solver contributes through the hooks below, so that it can add
  //! work at defined moments but can never skip the model's own.
  //! @return false if any phase fails
  bool setup
    ();

protected:

  //! @brief Wall-clock marks for the setup phase timing (DISPLAY_LEVEL >= 2).  Members rather than locals
  //! because setup() is a skeleton and the phases run partly in a derived solver's hooks.
  std::chrono::steady_clock::time_point _t_phase{};
  std::chrono::steady_clock::time_point _t_setup0{};

  //! @brief Report the time since the previous mark under DISPLAY_LEVEL >= 2, and re-mark.
  void _phase
    ( char const* name )
    {
      if( options.DISPLAY_LEVEL >= 2 ){
        auto const now = std::chrono::steady_clock::now();
        double const ms = std::chrono::duration<double,std::milli>( now - _t_phase ).count();
        _disp( 2 ) << "OCFESLV::setup [timing] " << std::left << std::setw(30) << name
                  << std::right << std::fixed << std::setprecision(2)
                  << std::setw(11) << ms << " ms" << std::defaultfloat << "\n";
        _t_phase = now;
      }
    }

  //! @brief Reference values at which the classification evaluates the symbol's coefficients.  FFModel uses the
  //! MODEL's own reference point: a declared reference function evaluated at the domains' reference coordinates,
  //! else the declared scalar reference, else zero; constants take their values and domains their reference
  //! coordinate (the midpoint unless declared).  No mesh is involved.
  //! A derived solver may override this -- OCFESLV samples the collocated field instead, which it does by
  //! default.
  virtual void _classification_reference
    ( std::vector<double>& state_ref, std::vector<double>& input_ref,
      std::vector<double>& cst_ref,   std::vector<double>& dom_ref );


  //! @brief Compute and cache the PDE classification at the reference point _classification_reference() gives.
  bool _classify_pde
    ();

  //! @brief Test whether a set of variables is a subset of another set
  static bool _subset
    ( std::set<FFVar,lt_FFVar> const& set1, std::set<FFVar,lt_FFVar> const& set2 );

  //! @brief Calculate the set difference between a set and a subset
  static std::set<FFVar,lt_FFVar> _setminus
    ( std::set<FFVar,lt_FFVar> set1, std::set<FFVar,lt_FFVar> const& set2 );

  //! @brief Create a string of the elements of a set
  static std::string _strset
    ( std::set<FFVar,lt_FFVar> const& set1, std::string const& sep );

  //! @brief Test whether a set of variables is a subset of another set
  template <typename U>
  static bool _subset
    ( std::set<FFVar,lt_FFVar> const& set1, std::map<FFVar,U, lt_FFVar> const& set2 );

  //! @brief Calculate the set difference between a set and a subset
  template <typename U>
  static std::map<FFVar,U,lt_FFVar> _setminus
    ( std::map<FFVar,U,lt_FFVar> set1, std::set<FFVar,lt_FFVar> const& set2 );

  //! @brief Create a string of the elements of a set
  template <typename U>
  static std::string _strset
    ( std::map<FFVar,U,lt_FFVar> const& set1, std::string const& sep );

  //! @brief Preflight the declared model before it is imported: _check_model() on the declared collections.  A
  //! derived solver may add its own checks, calling this version first.
  //! @return false to fail the setup with SetupStatus::INCONSISTENT_MODEL
  //! @brief The bounds [lo,up] an output point on domain @p var must lie in; _check_model() passes the domain's
  //! own.  A derived solver may widen them (OCFESLV: the march span, whose working domain is one element).
  virtual void _output_point_bounds
    ( FFVar const& var, double& lo, double& up )
    const
    { (void)var; (void)lo; (void)up; }

  virtual bool _validate_model
    ()
    { return _check_model( _usr._dagUsr, _usr._vCstUsr, &_usr._vCstValUsr, _usr._mDomUsr, _usr._mVarUsr,
                           _usr._mInpUsr, _usr._mEqnUsr, _usr._mFctUsr ); }

  //! @brief The declared model has been replaced by set_model(); does nothing in FFModel.  A derived solver
  //! drops its own DECLARED data here (OCFESLV: the moving-mesh specifications).
  virtual void _on_model_replaced
    ()
    {}

  //! @brief A setup is starting; does nothing in FFModel.  A derived solver clears what a re-setup invalidates.
  //! @return false to fail the setup
  virtual bool _on_setup_begin
    ()
    { return true; }

  //! @brief The model is imported and about to be reduced; does nothing in FFModel.
  //! @return false to fail the setup
  virtual bool _on_before_reduction
    ()
    { return true; }

  //! @brief The model is reduced and eliminated, and is about to be classified; does nothing in FFModel.  A
  //! derived solver discretises it here.
  //! @return false to fail the setup
  //! @brief apply the caller's options before setup() reads any of them.  FFModel's own options ARE
  //! the applied configuration; a derived solver overrides this to slice its own (richer) options in, which
  //! is the single, greppable moment of transfer -- and the only place a solver-side option may influence
  //! what the model is asked to do.  Default: nothing to apply.
  virtual void _apply_options
    ()
    {}

  virtual bool _on_model_discretise
    ()
    { return true; }

  //! @brief The reduction plan has run and the boundary closures are about to be generated; does nothing in
  //! FFModel.
  //! @return false to fail the setup
  virtual bool _on_before_closures
    ()
    { return true; }

  //! @brief The model is classified, reduced to index one and closed; does nothing in FFModel.  A derived solver
  //! completes its own setup here.
  //! @return false to fail the setup
  virtual bool _on_model_ready
    ()
    { return true; }

  //! @brief The declared model has just been imported onto a new working DAG; does nothing in FFModel.  The
  //! positional index the import produced is passed on: a derived class has no other way to remap its own
  //! members onto that DAG.
  virtual void _on_model_imported
    ( std::vector<FFVar> const&, std::map<FFVar::pt_idVar,size_t> const& )
    {}

  //! @brief The classification is about to be recomputed; does nothing in FFModel.
  virtual void _on_classify_begin
    ()
    {}

  //! @brief The classification is determined but not yet published; does nothing in FFModel.  The reference
  //! values it was computed at are passed on, since a derived solver has no other way to obtain them.
  virtual void _on_classify_ready
    ( std::vector<double> const&, std::vector<double> const&,
      std::vector<double> const&, std::vector<double> const& )
    {}

  //! @brief The classification has been published (_classified set, serial bumped); does nothing in FFModel.
  //! @return false to fail the classification
  virtual bool _on_classify_published
    ()
    { return true; }

  //! @brief Notification that the declared model changed; does nothing in FFModel.
  //! Called after each mutation by the add_*(), reset_*() and set_*() methods and update_ref().  A
  //! derived solver overrides it to invalidate its caches (OCFESLV clears _jacReady).
  virtual void _on_model_changed( ModelChange const ) {}

public:

  //! @brief Exception type of the declaration API (codes e.g. INDEX, CSTVAL).  The class is defined in OCBase and
  //! named here: the same type, so `catch( OCFESLV::Exceptions& )` also catches model-layer throws.  Transitional
  //! until the header split (NOTES_20260911k_rev192 Sec.3).
  typedef OCBase::Exceptions Exceptions;


  //! @brief Map from domain variable to its mc::FFDom: the extent, and the finite-element mesh a solver will
  //! use on it.  Owned here, as EqnOptions is: FFModel stores the mesh and passes it on, and interprets none of
  //! it.  The type is defined in ocbase.hpp (one definition, so that OCBase::var_domain() can be overridden
  //! with this map).
  typedef OCBase::t_Dom  t_Dom;


  //! @brief Role of an equation in the model.
  //!
  //! INTERIOR equations hold on the interior of their domains (volume PDEs and algebraic residuals)
  //! and are the rows of the block principal symbol.  INITIAL, BOUNDARY and INTERFACE equations are
  //! lower-dimensional trace constraints -- initial conditions, exterior boundary conditions, and
  //! user-supplied transmission conditions between subdomains -- and never enter the volume
  //! principal-symbol classification.  LINK equations are the pointwise auxiliary definitions
  //! generated by order reduction; they may enter classification.  SURFACE equations are PDEs posed
  //! on a lower-dimensional manifold and may form their own classification block.  DIAGNOSTIC
  //! equations are evaluated only.  AUTO is resolved when the equation is added
  //! (_normalise_equation_options()).  How a solver imposes each role -- SAT reception, trace
  //! donation -- is documented with its interface plan (OCFESLV::_resolve_interface_type()).
  enum class EqnRole
  {
    AUTO = 0,  //!< resolved on add: BOUNDARY if some domain mask is LB or UB alone, else INTERIOR
    INTERIOR,  //!< volume/interior PDE or algebraic residual
    INITIAL,  //!< initial-condition trace residual
    BOUNDARY,  //!< exterior boundary-condition trace residual
    INTERFACE,  //!< user-supplied transmission residual between physical subdomains
    LINK,  //!< auxiliary first-order link equation generated by reduce_order()
    SURFACE,  //!< PDE or algebraic residual posed on a lower-dimensional boundary or interface
    DIAGNOSTIC  //!< evaluated only: not classified, receives no SAT, never donated
  };

  //! @brief Options
  struct Options
  {
    //! @brief Read a boolean option default from the environment.
    //!
    //! The AUDIT_* diagnostics are wanted across the whole driver corpus, and editing
    //! every driver to flip a flag is neither practical nor reversible.  reset() is
    //! called from the Options constructor, so a variable set for the process supplies
    //! the default for every OCFESLV built in it -- no driver edit, no API change, and a
    //! driver that sets the field explicitly still wins.
    //! False spellings: unset, "", "0", "false", "off", "no" (case-insensitive).
    //!
    //! Companion for numeric options.  Without it INTERFACE.TRACE_PROJ_EPS would have no run-time hook, and the
    //! one knob whose value the projector's own note says must be swept -- "Sweep
    //! INTERFACE.TRACE_PROJ_EPS to confirm the product is eps-independent before reading it" -- could
    //! not be swept without editing a driver.  A diagnostic that names a sweep and offers no
    //! way to perform it is not a diagnostic.
    static double _env_dbl
      ( char const* name, double dflt )
      {
        char const* v = std::getenv( name );
        if( !v || !*v ) return dflt;
        char* end = nullptr;
        double const d = std::strtod( v, &end );
        return ( end && end != v ) ? d : dflt;
      }

  //! @brief ONE knob for the PAIR (R1 + R2), CRONOS_REUSE_AWARE.
  //! DEFAULT ON.  Evidence: knob-ON corpus sweep 20260910_200924 (only PDE20f moves --
  //! n=m 565->605, the R2 fingerprint; both XFAILs -> PASS; R1 recouples 4,312 edges on 12
  //! drivers with zero outcome change; R2 touches one other driver harmlessly), the
  //! perturbed-start matrix on the PSA_IMP drivers GREEN, and OCFE_PDE20g's sigma-reuse family
  //! EQUIVALENT (RED_FULL == RED_MAIN) on all five models, 30 cells.  CRONOS_REUSE_AWARE=0
  //! restores the pre-pair behaviour for bisection.
  //! MEASURED on OCFE_PDE20f model C / RED_FULL (NOTES_20260910m): each half alone is partial
  //! or harmful; together, every exact-imposition cell reproduces RED_MAIN to the digit.
  //!   R2  the value-slaved drop must see through reduction aliases: a defining constraint
  //!       whose bare terms include an order-reduction auxiliary (a lifted derivative) is
  //!       NOT derivative-free, so the state it defines is NOT value-slaved.  Under sigma
  //!       reuse  w - cf*Dz_u = 0  fooled the literal test and w's 20 seam claims were dropped
  //!       as implied -- by a carrier (Dz_u's rescued edges) that was itself degenerate.
  //!   R1  S1c on the EXACT path too: coupling = c0 on algebraic receivers under IC_TRACE /
  //!       IC_STRONG, not only IC_WEAK.  Un-degenerates the 20 "tau-block 2/3" clusters.
  //!       Without R2 the tau it introduces contaminates the two-term identity (4.49e-04,
  //!       measured 2026-09-10j); with R2 the identity's state is claimed and carried on the
  //!       PDE row, and the pair is complete.
  //! Plan-changing under exact imposition by construction (cluster rank).  Gate: the
  //! perturbed-start matrix, not a sweep alone.
  static int _reuse_aware()
  {
    static int const v = 1;   // CRONOS_REUSE_AWARE (retired 2026-10-07, WORKPLAN 3.C); rev183: default ON   // rev183: default ON
    return v;
  }

  //! @brief Force solver verbosity ON from the environment (never off).
    //!
    //! FORCE-ON ONLY, deliberately.  A driver that sets SOLVE.VERBOSE=false does so to keep
    //! its own table readable; a diagnostic run overrides that, and the reverse never
    //! happens silently.
    //!
    //! Applied at BOTH the constructor and the option-ADOPTION site.  Applied at construction
    //! alone it is overwritten by any driver that sets the field afterwards (OCFE_MBC5.cpp does
    //! `oc.options.SOLVE.VERBOSE = false;`), and CRONOS_SOLVE_VERBOSE=1 then changes nothing
    //! while the output still looks plausible -- a diagnostic switch that only appears to work.
    //! The adoption site is the last
    //! point before the solver reads the flag, so a force-on there beats any driver
    //! assignment.
    static void _force_verbose_from_env( bool& flag )
    {
      char const* v = std::getenv( "CRONOS_SOLVE_VERBOSE" );
      if( !v ) return;
      std::string const w( v );
      if( w == "1" || w == "true" || w == "TRUE" ) flag = true;
      else if( w != "0" && w != "false" && w != "FALSE" )
        std::cerr << "OCFESLV ** CRONOS_SOLVE_VERBOSE='" << w << "' not recognised (1|0);"
                     " leaving SOLVE.VERBOSE as the driver set it." << std::endl;
    }

    //! @brief Integer environment knob, clamped to [lo,hi]; @p def when unset or unparseable.
    static int _env_int
      ( char const* name, int def, int lo, int hi )
      {
        char const* e = std::getenv( name );
        if( !e || !*e ) return def;
        int v = def;
        try{ v = std::stoi( e ); } catch(...){ return def; }
        return v < lo? lo: ( v > hi? hi: v );
      }

    static bool _env_flag
      ( char const* name, bool dflt = false )
      {
        char const* v = std::getenv( name );
        if( !v || !*v ) return dflt;
        std::string s( v );
        for( auto& c : s ) c = static_cast<char>( std::tolower( static_cast<unsigned char>( c ) ) );
        if( s == "0" || s == "false" || s == "off" || s == "no" ) return false;
        return true;
      }

    //! @brief Reset options
    //! @brief apply every CRONOS_* override on top of the current values.  Called by reset()
    //! so today's behaviour is unchanged; separated so the model/solver ENVIRONMENT BOUNDARY is one place.
    //! Each read uses the CURRENT value as its default, so a caller may set fields programmatically and then
    //! apply the environment on top -- and D4b can move this CALL to OCFESLV, leaving FFModel plain settings.
    void apply_environment_defaults
      ()
      {
        // rev142b: 0 = no cap imposed (the runtime keeps whatever the environment set).  rev321: back HERE -- rev308
        // had moved this line into OCFESLV::Options, leaving FFModel::options.MAXTHREAD uninitialised until setup().
        MAXTHREAD = (size_t)_env_dbl( "CRONOS_MAXTHREAD", 0. );
        CLASSIFY.ROBUST     = ROBUST_OFF;   // CRONOS_CLASSIFY_ROBUST retired 2026-10-07 (WORKPLAN 3.A): the option decides
      }

    void reset
      ()
      {
        DISPLAY_LEVEL       = 0;
        REDUCE.ORDER        = RED_FULL;
        REDUCE.NONLINEAR_PARTIALS = false;
        REDUCE.HIDDEN_IC    = true;
        TTOL                = 1e-9;
        CLASSIFY.MODE            = CLASS_AUTO;
        CLASSIFY.NSAMPLE    = 32;    // rev284: default on (sweep H)
        CLASSIFY.ROBUST     = ROBUST_OFF;
        CLASSIFY.IMAG_TOL   = 1e-8;
        AUTO.HYP_CLOSURE    = true;    //!< auto outflow closure of hyperbolic blocks (was CRONOS_AUTO_HYP_CLOSURE)
        AUTO.DIFF_ELIM      = false;   //!< auto differential-elimination of algebraically-determined derivatives (off by default)
        apply_environment_defaults();   // rev307 (D4a): the boundary, in one place
      }
      
    //! @brief Constructor
    Options()
      {
        reset(); 
      }
      
    //! @brief Assignment of mc::OCFESLV::Options
    Options& operator=
      ( Options const& opt )
      {
        DISPLAY_LEVEL       = opt.DISPLAY_LEVEL;
        REDUCE.ORDER        = opt.REDUCE.ORDER;
        REDUCE.NONLINEAR_PARTIALS = opt.REDUCE.NONLINEAR_PARTIALS;
        REDUCE.HIDDEN_IC    = opt.REDUCE.HIDDEN_IC;
        TTOL                = opt.TTOL;
        CLASSIFY.MODE            = opt.CLASSIFY.MODE;
        CLASSIFY.NSAMPLE    = opt.CLASSIFY.NSAMPLE;
        CLASSIFY.ROBUST     = opt.CLASSIFY.ROBUST;
        CLASSIFY.IMAG_TOL   = opt.CLASSIFY.IMAG_TOL;
        MAXTHREAD           = opt.MAXTHREAD;
        AUTO.HYP_CLOSURE    = opt.AUTO.HYP_CLOSURE;
        AUTO.DIFF_ELIM      = opt.AUTO.DIFF_ELIM;
        return *this;
      }

    //! @brief Console verbosity for setup()/eval() (default: 0).  Tiered:
    //!   0 = silent: only errors, validation failures, and warnings print.
    //!   1 = key structural decisions: closures, interface-plan DROP / KEEP-
    //!       EXPLICIT, continuity redundancy / value-slaved drops, order and
    //!       tau reduction summaries.
    //!   2 = detailed diagnostics (includes everything at level 1): phase and
    //!       [iplan] timings, trace conditioning (rcond / SVD / sparsity),
    //!       per-predicate timing accumulators, verbose detection and edge
    //!       breakdowns.
    //! Errors and warnings are NEVER gated by this level.
    int DISPLAY_LEVEL;
    enum ReductionType
    {
      RED_NONE = 0,  //!< Do not perform order reduction during setup
      //! Auxiliaries are SHARED in both modes -- one per distinct derivative, always; two names for
      //! one continuous quantity with nothing constraining them to agree is not a policy choice (PDE12's
      //! biharmonic minted 288 duplicate unknowns and its interface claims pinned nothing).  What the mode
      //! selects is whether rows processed BEFORE an auxiliary existed are rewritten to use it.
      RED_MAIN = 1,  //!< Reduce high-order derivatives, leaving earlier rows their explicit derivative
      RED_FULL = 2   //!< Reduce high-order derivatives and rewrite earlier rows to use the auxiliaries
    };


    enum RobustnessMode
    {
      ROBUST_OFF    = 0,  //!< Do not sample: classify at the reference point only (default)
      ROBUST_REPORT,      //!< Sample and report disagreement; the verdict still comes from the base point
      ROBUST_STRICT       //!< Sample and treat disagreement as the classification being reference-sensitive
    };


    //! @brief How setup() classifies the PDE blocks -- the values of CLASSIFY.MODE.
    enum ClassifyType
    {
      CLASS_NONE   = 0,  //!< Do not perform PDE classification during setup
      CLASS_AUTO   = 1,  //!< Attempt PDE classification during setup and proceed even if it fails
      CLASS_STRICT = 2   //!< Attempt PDE classification during setup and interrupt when it fails
    };

    //! @brief Order reduction of high-order derivatives.
    struct t_Reduce
    {
      //! @brief Order-reduction mode applied automatically by setup.
      //! RED_NONE leaves user equations unchanged.  RED_MAIN performs first-order
      //! reduction without substituting matching first-order derivative nodes in
      //! pre-existing equations.  RED_FULL also reuses the new auxiliary states
      //! in matching pre-existing equations; this is the default.
      ReductionType ORDER;
      //! @brief Materialise first-order state-partials that appear NON-AFFINELY in a residual.
      //!
      //! reduce_order already peels high-order partials, and already materialises the partials
      //! carried INSIDE an OpI/OpEval operand (materialize_operand_partials), because a nonlocal
      //! reduction LINK cannot expose an inner OpP.  It does NOT touch a first-order OpP(state,dir)
      //! sitting in an ordinary residual, which is right when the residual is AFFINE in it -- the
      //! principal-symbol proxy and the collocation Jacobian both handle a variable-coefficient
      //! derivative.  It is NOT right when the derivative is wrapped nonlinearly: a DIVISION by a
      //! derivative (the ALE/moving-mesh 1/x_xi, u_xi/x_xi, u_xi/x_xi^2) or a transcendental of one
      //! (an arclength monitor sqrt(1+u_xi^2)).  Those force a driver to hand-roll the derivative as
      //! a state with its own defining equation -- and hand-rolled auxiliaries are invisible to the
      //! classifier AS auxiliaries, which is how the multiply-defined flux states and the
      //! rectangular symbols in the moving-mesh drivers arose.
      //!
      //! When true, such a partial is minted as a derivative auxiliary Dp with a defining LINK
      //!     OpP(state,dir) - Dp = 0
      //! and substituted everywhere -- exactly the mint/LINK/reuse protocol the OpI/OpEval operand
      //! path already uses.  The residual then contains only the STATE Dp inside the nonlinearity,
      //! and the raw OpP survives solely in a LINK, which is an ordinary differential equation the
      //! proxy pass resolves.
      //!
      //! The trigger is per-partial and minimal: a partial is materialised iff the residual is not
      //! AFFINE in THAT partial alone (mc::FFDep type != L).  So `D*u_z`, `u*u_z`, `exp(-u)*u_z` are
      //! untouched, and in `(c-x_t)*u_xi/x_xi` ONLY x_xi is materialised -- x_t and u_xi remain
      //! ordinary variable-coefficient partials.  Every residual that is affine in each of its
      //! first-order partials -- the entire existing corpus -- is therefore bit-identical.
      //!
      //! Default false: opt in, sweep, then consider promoting.  Requires REDUCE.ORDER != RED_NONE.
      bool NONLINEAR_PARTIALS;

      //! @brief Materialise the HIDDEN constraints of a high-index reduction as INITIAL rows at the evolution
      //! boundary, so the model declares only the FREE initial data (default true since 2026-10-04; false was the
      //! former convention of one declared initial condition per differential state, made consistent by hand).  A reduction that differentiates a
      //! constraint k times implies k lower levels which hold at the initial point but are not rows: with this off,
      //! the model must declare one initial condition per evolution-differential state and make the redundant ones
      //! consistent by hand (initial_data() says how many they are); with it on, those levels ARE rows, the free data
      //! alone makes the system square, and a full declaration is reported as a surplus by the DOF audit.
      //! Which data is free is a modelling choice -- for the index-3 pendulum {x,u} works, {x,y} does not.
      bool HIDDEN_IC;
    } REDUCE;

    //! @brief Classification of each block's principal symbol.
    struct t_Classify
    {
      ClassifyType MODE;          //!< PDE classification during setup(): CLASS_AUTO (default), CLASS_NONE or CLASS_STRICT
      //! @brief Reference-robustness sampling (see RobustnessMode).
      //! (0 = off, 1 = report, 2 = strict).
      RobustnessMode ROBUST;
      //! @brief Number of angular samples for no-evolution-domain classification.
      unsigned NSAMPLE;
      //! @brief Imaginary-part/rank tolerance used by automatic classification.
      double IMAG_TOL;
    } CLASSIFY;

    //! @brief Automatic model transformations applied at setup(): boundary closures and differential elimination.
    struct t_Auto
    {
      bool   HYP_CLOSURE;     //!< auto closure of the OUTFLOW faces of hyperbolic blocks: appends the missing
                              //!< outgoing-characteristic rows (default on).  Off: the model must close every outflow
                              //!< face itself (e.g. its PDE extended to the face) -- setup REFUSES (HYP_CLOSURE_MISSING)
                              //!< when rows are missing, as an under-determined system would be solved silently wrong
      bool   DIFF_ELIM;       //!< auto differential-elimination of algebraically-determined derivatives (substitute + relocate consumer onto source face); OFF by default
    } AUTO;




    //! @brief Two time points closer than this are the SAME point (1e-9).  A consumer that stages a model -- an
    //! integrator, a marching solve, anything that stops at boundaries -- must decide whether a function's time
    //! point coincides with an element boundary, and an exact comparison is the wrong test for a value that was
    //! typed or computed.  Absolute, in the units of the evolution domain.
    double TTOL;

    size_t MAXTHREAD;            //!< cap on the threads the numerical backends may
                                 //!< use inside setup() and solve().  SPQR, SuperLU, UMFPACK,
                                 //!< KLU and every dense LAPACK call take their parallelism from
                                 //!< the SAME BLAS/OpenMP layer, so ONE cap covers them all --
                                 //!< SPQR has no thread pool of its own (its TBB path was removed
                                 //!< upstream, and cc.SPQR_nthreads is now inert).
                                 //!< 0 (default) imposes NO cap: the runtime is left exactly as
                                 //!< the environment set it.  n>0 caps to n for the duration of
                                 //!< the call and RESTORES the previous setting on every exit.
                                 //!< NOTE the deliberate difference from FFGraph::Options::
                                 //!< MAXTHREAD, which reads 0 as "hardware_concurrency": here 0
                                 //!< means "do not touch", so an outer OMP_NUM_THREADS or a
                                 //!< cluster allocation is respected rather than overridden.
                                 //!< MEASURED (PDE5, 9 audit windows): cap=1 is 18% FASTER than
                                 //!< the uncapped run -- 152.4 s against 186.0 s of rank time.
                                 //!< Environment: CRONOS_MAXTHREAD.
    //!@}
  } options;

  //! @brief Per-equation options: model role and classification block, plus solver annotations.
  //! role, block_id and participate_in_classification are model attributes.  interface_type,
  //! receive_sat and donate_for_state_continuity annotate how a solver may impose interface
  //! conditions on the equation; FFModel stores them without interpreting them.  add_equation()
  //! passes the options through _normalise_equation_options(), so the stored options may differ from
  //! those supplied.
  //! POLYMORPHIC: a solver derives from this to carry its own per-equation fields, and t_Eqn
  //! holds it by shared_ptr so that derived object is never sliced on storage.
  struct EqnOptions
  {
    EqnRole role;                         //!< equation role (AUTO resolved on add)
    bool    participate_in_classification;  //!< enters block principal-symbol classification
    int     block_id;                     //!< classification block
    //! @brief True when the role was resolved from AUTO (by _normalise_equation_options()): such a role is
    //! re-resolved once the evolution domain is known (_reresolve_auto_roles()) -- an equation on LB alone of the
    //! EVOLUTION domain is INITIAL, not BOUNDARY.
    bool    role_auto = false;

    //! @brief Role-derived defaults: INITIAL, BOUNDARY, INTERFACE and DIAGNOSTIC equations do not enter
    //! classification.  (how a SOLVER imposes each role -- SAT reception, trace donation, the local
    //! interface type -- is no longer carried here; a solver derives from this struct to add it.)
    EqnOptions
      ( EqnRole r = EqnRole::AUTO, int b = 0 )
      : role( r ), participate_in_classification( true ), block_id( b )
      {
        if( role == EqnRole::INITIAL || role == EqnRole::BOUNDARY ||
            role == EqnRole::INTERFACE || role == EqnRole::DIAGNOSTIC )
          participate_in_classification = false;
      }

    EqnOptions( EqnOptions const& other ) = default;
    EqnOptions& operator=( EqnOptions const& other ) = default;

    //! @brief EqnOptions is a base.  Without this, storing a derived object slices it away silently.
    virtual ~EqnOptions() = default;

    //! @brief an independent copy of the SAME dynamic type.  A derived options type MUST override this,
    //! or a copied environment would silently lose its derived fields.
    virtual std::shared_ptr<EqnOptions> clone
      ()
      const
      { return std::make_shared<EqnOptions>( *this ); }
  };


  //! @brief automatic algebraic-constraint boundary closure for parabolic/DAE blocks (default ON).  The
  //! closure is COVERAGE-AWARE -- it adds a closure at a face only where no existing equation already pins the
  //! target state -- so leaving it on is harmless when the user writes the closure too, and switching it off
  //! only exposes the model's UNclosed system (a diagnostic, e.g. when chasing a determinacy failure).  No
  //! driver ever set it false.  (The CRONOS_AUTO_ALG_CLOSURE variable is retired, 2026-10-07.)
  static bool _knob_AUTO_ALG_CLOSURE
    ()
    { return true; }   // CRONOS_AUTO_ALG_CLOSURE (retired 2026-10-07, WORKPLAN 3.B batch 2b)

  //! @brief Chain-aware principal symbol.  A participating INTERIOR row with no differentiated state but a bare
  //! reduction auxiliary whose definition is a derivative (RED_FULL's fully reduced balance row) is read, for the symbol
  //! only, with that auxiliary inlined one level, so it keeps its derivative column and is the natural receiver of the
  //! parent's continuity claim (as under RED_MAIN).  Default on; CRONOS_CHAIN_SYMBOL=0 switches it off.
  static bool _knob_CHAIN_SYMBOL
    ()
    { return true; }   // CRONOS_CHAIN_SYMBOL (retired 2026-10-07, WORKPLAN 3.C)

  //! @brief Map from state or input to the domain variables it depends on (empty set: lumped).
  typedef std::map< FFVar, std::set< FFVar, lt_FFVar >, lt_FFVar >                                    t_Var;

  //! @brief Map from domain variable to a domain mask (FFDom::ALL, LB, UB, or combinations such as ALL-LB).
  typedef std::map< FFVar, int, lt_FFVar >                                                            t_EqnDom;

  //! @brief Map from domain variable to a fixed coordinate.
  typedef std::map< FFVar, double, lt_FFVar >                                                         t_FctDom;

  //! @brief Equation record, in declaration order.
  struct t_Eqn
  {
    FFVar       var;  //!< residual expression (equation: var = 0)
    t_EqnDom    dom;  //!< domain masks per domain variable
    //! Held by pointer so a solver's derived options survive storage.  Never null after
    //! add_equation(); the model reads only the base fields through it.
    std::shared_ptr<EqnOptions> opt;  //!< normalised equation options (polymorphic)
  };

  //! @brief Output-function kind.
  enum class FctKind
  {
    POINT,  //!< one scalar value, evaluated at a fixed point
    DISTRIBUTED  //!< one value per discretisation node selected by the grid masks
  };

  //! @brief Output-function record, in declaration order.
  struct t_Fct
  {
    FFVar       var;                         //!< Output expression.
    FctKind     kind = FctKind::POINT;       //!< Point or distributed output.
    t_FctDom    point;  //!< fixed coordinates
    std::map<FFVar,int,lt_FFVar> side;  //!< per point direction: FFDom::PLUS where given (absent: FFDom::MINUS)
    t_EqnDom    grid;  //!< domain masks of the distributed dimensions
    mutable size_t row0 = 0;  //!< first row in the output array (set by the solver)
    mutable size_t nrow = 1;  //!< number of output rows (set by the solver for DISTRIBUTED)
  };

  //! @brief Declared collocation of an input on the elements of one domain: the input shares the domain's
  //! grid and differs only in its local representation -- n_node values per element (0: as the domain's
  //! states, i.e. the same order) at nodes of the given type.  With the declared continuity (t_Cont) this
  //! fixes the control's function space: n_node = 1 with DISCONTINUOUS is a piecewise-constant control.
  struct InpColloc
  {
    FFDom::TYPE  type;         //!< node family on each element
    size_t       n_node = 0;   //!< values per element; 0: as the domain
  };

  //! @brief Map from input to its declared collocation per domain variable.
  typedef std::map< FFVar, std::map<FFVar,InpColloc,lt_FFVar>, lt_FFVar >  t_InpColloc;

  //! @brief Map from input to its declared continuity level per domain variable (see
  //! OCFESLV::InputContinuity); a missing entry means DISCONTINUOUS.
  typedef std::map< FFVar, std::map< FFVar, int, lt_FFVar >, lt_FFVar >                                 t_Cont;

  //! @brief Equations in declaration order.
  typedef std::vector< t_Eqn >                                                                         t_Eqns;

  //! @brief Output functions in declaration order.
  typedef std::vector< t_Fct >                                                                         t_Fcts;

  //! @brief A TRANSITION on the evolution direction (add_transition): left(tau^-) - right(tau^+) = 0, componentwise.
  //! States mentioned by no transition are continuous at tau.
  struct t_Transition
  {
    std::vector<FFVar> left;   //!< evaluated at tau^- (the limit from below)
    std::vector<FFVar> right;  //!< evaluated at tau^+ (the limit from above)
    FFVar              dom;    //!< the evolution direction
    double             tau = 0.;
    bool               evaluated = false;  //!< declared with EVALUATIONS (add_transition( left, right )): dom and tau
                                           //!< are read from them at setup, and each is replaced by its operand
  };
  typedef std::vector< t_Transition > t_Transitions;

  //! @brief The states a group of transitions at one tau makes JUMP: the post-jump set S+ (the states the right sides,
  //! read at tau^+, depend on), with each transition component matched to a distinct state of S+ it depends on.  A
  //! solver breaks exactly these states' continuity at tau, and imposes the transition rows in its place.
  struct t_TransitionJump
  {
    double                          tau = 0.;
    std::vector<size_t>             comp;     //!< (transition index, component) packed as index*65536 + component
    std::vector<FFVar>              jump;     //!< S+, in matching order: comp[k] is matched to jump[k]
  };
  typedef std::vector< t_TransitionJump > t_TransitionJumps;

  //! @brief The model's DEGREE-OF-FREEDOM BALANCE, counted symbolically: one symbol N_d per domain, so no mesh is
  //! needed.  Each equation contributes, for every domain it mentions, N (the whole domain), N-1 (one face excluded),
  //! N-2 (both faces) or 1 (a single face), multiplied together; each state contributes the product of its domains'
  //! symbols.  @a diff is rows - unknowns as a multilinear polynomial (monomial -> coefficient, the empty monomial
  //! being the constant), @a str its readable form, and the model is @a balanced iff every coefficient is zero.
  //! Evaluated on any grid, a non-zero @a diff is exactly the surplus (or, negative, the deficit) of collocated rows.
  //! The balance is over the WHOLE model: the model has no state-to-block map, so a multi-block result says the global
  //! system is square, not that each block is.
  struct t_DofBalance
  {
    bool                                     balanced = true;
    std::map< std::set<std::string>, long >  diff;
    std::string                              str = "0";
  };

  //! @brief What the model's INITIAL rows amount to once the index reduction is taken into account.  A high-index
  //! system carries HIDDEN constraints -- the lower differentiations of the constraint the reduction differentiated
  //! (for the index-3 pendulum, the position and velocity levels) -- which hold at the initial point but are NOT rows
  //! of this model.  So only @a free of the initial values may be CHOSEN (@a differential evolution-differential
  //! states, less @a hidden levels) -- a property of the model, whatever @a declared says -- and @a redundant declared
  //! rows must be CONSISTENT with constraints the model does not enforce; nothing here checks that they are.  For an
  //! index-0 or index-1 model @a hidden is zero and every declared row is free.
  struct t_InitialData
  {
    size_t declared     = 0;   //!< INITIAL rows in the working model
    size_t differential = 0;   //!< states differentiated in the evolution direction
    size_t hidden       = 0;   //!< hidden constraint levels implied by the reduction
    size_t free         = 0;   //!< differential - hidden: what may be chosen
    size_t redundant    = 0;   //!< declared - free, when more is declared than may be chosen
    bool   materialised = false; //!< whether the hidden levels were materialised as rows (REDUCE.HIDDEN_IC)
    bool   reduced      = false;
  };

  //! @brief A DEFERRED VALUE: a functional of the solution that cannot be evaluated during the solve, because it
  //! reduces the evolution direction away -- the integral of an expression over that direction, or its value at one
  //! point of it.  FFModel::setup() decides these (see _is_deferred_reduction()), rewrites the system so the reduction
  //! node becomes the INPUT @a input (initially 0, no state, no equation, no initial condition), and records one of
  //! these per capture; the value feeds outputs and objectives only, and an inner deferred value may not feed another
  //! evolution reduction (the causality refusal).  THE MODEL NEVER COMPUTES IT: a consumer evaluates each record after
  //! solving and writes it into @a input.  The contract per record is:
  //!  - @a accumulate == true:  value = INTEGRAL of @a source over the whole evolution domain.  How that integral is
  //!    formed, and whether the solve is windowed, is the consumer's business.
  //!  - @a accumulate == false: value = @a source evaluated at @a tau on the evolution direction.
  //! In BOTH kinds @a source is the expression the functional is taken OF -- the integrand, or the operand evaluated
  //! at tau -- never the reduction node itself, which is not evaluable pointwise.  (The comment this contract
  //! replaces said the integral kind recorded the OpI node; the records say otherwise, and so does the arithmetic in
  //! OCFESLV::_fill_captures, which evaluates @a source at quadrature points.)
  //! @a block_id names the model block the capture belongs to; a capture distributed over the remaining directions has
  //! @a input declared on them (a profile capture), otherwise it is scalar.
  //! The two kinds are what an ODE/DAE integrator already offers -- a quadrature variable and dense output at a point
  //! -- so a consumer that holds a model (rather than deriving from one) can implement this contract directly.
  //! @a source is an EXPRESSION of the model's states and inputs: a consumer evaluates it as it sees fit -- an
  //! integrator passes it to a quadrature right-hand side; OCFESLV, which reads values from its collocated solution,
  //! materialises it as a state with a defining row when it is not one already (in _on_model_discretise, before the
  //! classification), which is why such a model classifies DIFFERENTIAL_ALGEBRAIC under OCFESLV and not on its own.
  //! OCFESLV realises the contract as: LATCH -- @a source at tau, written when the window or element containing tau
  //! completes (static schedule, left-window tie-break at a boundary tau); ACCUM -- an element quadrature of
  //! @a source over the collocation, accumulated per window under marching, one full-domain evaluation otherwise.
  struct t_DeferredValue
  {
    FFVar  input;              //!< the model input that holds the value; written by the consumer after the solve
    FFVar  source;             //!< the expression the functional is taken of: the integrand, or the operand at @a tau
    std::set<FFVar,lt_FFVar> source_dom; //!< the directions @a source varies over (the evolution one included)
    double tau = 0.;           //!< LATCH only: the point on the evolution direction
    int    side = 0;           //!< LATCH only: FFDom::MINUS (tau^-, default) or FFDom::PLUS (tau^+)
    bool   accumulate = false; //!< true: the integral over the evolution domain; false: the value at @a tau
    int    block_id = 0;       //!< the model block the capture belongs to
  };

protected:
  //! @brief Check a model's collections for consistency -- the DAG exists; constants, domains, states, inputs,
  //! equations and outputs belong to it; every domain a variable, equation or output uses is declared; equation
  //! and output masks are valid; every variable an equation or output uses is declared, with domains consistent
  //! with it; equation domains are tight; output points lie inside their domains.  Prints the first violation to
  //! std::cerr and returns false.  Validation only -- nothing is counted or written -- except that each
  //! equation's and output's subgraph is stored in @p sgEqn / @p sgFct when given.
  bool _check_model
    ( FFGraph* dag, std::vector<FFVar> const& vCst, std::vector<double> const* vCstVal,
      t_Dom const& mDom, t_Var const& mVar, t_Var const& mInp,
      t_Eqns const& mEqn, t_Fcts const& mFct,
      std::vector<FFSubgraph>* sgEqn = nullptr, std::vector<FFSubgraph>* sgFct = nullptr )
    const;

public:

  //! @brief Coordinates keyed by domain variable, argument of reference functions.
  typedef std::map<FFVar,double,lt_FFVar>                                                              t_Coord;
  typedef std::map<FFVar,int,lt_FFVar>                                                                 t_Side;   //!< FFDom::MINUS / FFDom::PLUS per direction (absent: MINUS)

  //! @brief Reference function of the domain coordinates, for states and inputs.
  //! The argument is keyed by domain variable.  setup() and deep copies wrap user functions so that
  //! both the user-DAG and the working-DAG domain variables are accepted as keys.
  typedef std::function<double(t_Coord const&)>                                                        t_Fun;

protected:
  //! @brief The model's own reference value for @p var at @p coord: a declared reference function or scalar, a
  //! constant's value, or -- for an auxiliary minted by order reduction -- the derivative of its parent's value
  //! in the direction it differentiates, by central difference.  Zero when nothing is declared, which is the
  //! right answer for a constant reference and the wrong one wherever the value divides: see the non-finite
  //! warning in _classification_reference().
  double _model_reference_value
    ( FFVar const& var, t_Coord const& coord, unsigned const depth = 0 )
    const;

public:

protected:

  //! @brief The model as declared, on the user DAG.
  //! Filled by the declaration API.  A solver's setup() imports the roots into its private working DAG
  //! with one FFGraph::insert call and rebuilds the working collections (_vCst, _vCstVal, _mDom, _mInp,
  //! _mVar, _mEqn, _mFct and the classification references) on it.
  struct UsrData
  {
    FFGraph*                           _dagUsr = nullptr;  //!< user DAG owning the declared variables
    std::vector<FFVar>                 _vCstUsr;  //!< constants
    std::vector<double>                _vCstValUsr;  //!< constant values (empty: not given)
    t_Dom                              _mDomUsr;  //!< declared domains (extent and mesh)
    t_Var                              _mInpUsr;  //!< inputs and their domains
    std::map<FFVar,std::vector<double>,lt_FFVar> _mInpFixUsr;   //!< inputs held at known values (see fix_input)
    t_InpColloc                        _mInpDiscUsr;  //!< solver: input discretisation (OCFESLV::add_input; moves in F1c)
    t_Cont                             _mInpContUsr;  //!< declared input continuity
    t_Var                              _mVarUsr;  //!< states and their domains
    t_Eqns                             _mEqnUsr;  //!< equations
    t_Fcts                             _mFctUsr;  //!< output functions
    t_Transitions                      _mTrnUsr;  //!< transitions (add_transition)
    FFVar                              _evolution_dom_varUsr;  //!< evolution direction (meaningful when _evolution_dom_setUsr)
    bool                               _evolution_dom_userUsr = false;  //!< evolution direction given by set_evolution_domain()
    bool                               _evolution_dom_setUsr  = false;  //!< _evolution_dom_varUsr holds a direction
    std::map<FFVar,double,lt_FFVar>    _classVarRefUsr;  //!< scalar references of states
    std::map<FFVar,double,lt_FFVar>    _classInpRefUsr;  //!< scalar references of inputs
    std::map<FFVar,double,lt_FFVar>    _classDomRefUsr;  //!< reference coordinates of domains
    std::map<FFVar,t_Fun,lt_FFVar>     _classVarRefFunUsr;  //!< reference functions of states
    std::map<FFVar,t_Fun,lt_FFVar>     _classInpRefFunUsr;  //!< reference functions of inputs


    //! @brief Clear the declared model; the user-DAG pointer is kept.
    void clear_model()
    {
      _vCstUsr.clear();
      _vCstValUsr.clear();
      _mDomUsr.clear();
      _mInpUsr.clear();
      _mInpContUsr.clear();
      _mVarUsr.clear();
      _mEqnUsr.clear();
      _mFctUsr.clear();
      _mTrnUsr.clear();
      _classVarRefUsr.clear();
      _classInpRefUsr.clear();
      _classDomRefUsr.clear();
      _classVarRefFunUsr.clear();
      _classInpRefFunUsr.clear();
      _evolution_dom_userUsr = false;
      _evolution_dom_setUsr  = false;
    }

    //! @brief Clear the declared model and the user-DAG pointer.
    void clear()
    {
      clear_model();
      _dagUsr = nullptr;
    }
  };

  //! @brief Private working DAG created by setup(); null before setup() and after its release.
  //! Never the user DAG.
  FFGraph*                           _dag =  nullptr;

  //! @brief True when _dag is a private working DAG owned, and deleted, by this object.
  bool                               _dagOwned = false;

  //! @brief The model as declared (see UsrData).
  UsrData                            _usr;

  //! @brief True after a successful setup(); cleared by model mutations.  Selects between the
  //! working and the declared collections in the var_*() getters.
  bool                               _issetup = false;

  //! @brief Domains (working copy, valid after setup()): the declarations remapped to the working DAG.  A
  //! solver may rewrite an entry in place -- OCFESLV collapses the evolution domain to the current window when
  //! marching -- so this is the extent in force, which is what ref() must read.
  t_Dom                              _mDom;

  //! @brief Constants (working copy, valid after setup()).
  std::vector<FFVar>                 _vCst;

  //! @brief Constant values (working copy).
  std::vector<double>                _vCstVal;

  //! @brief Inputs and their domains (working copy).
  t_Var                              _mInp;

  //! @brief States and their domains (working copy).
  t_Var                              _mVar;

  //! @brief Equations in declaration order (working copy; setup may append generated equations).
  t_Eqns                             _mEqn;

  //! @brief Output functions in declaration order (working copy).
  t_Fcts                             _mFct;
  t_Transitions                      _mTrn;   //!< transitions, on the working DAG (see add_transition)

public:

  //! @brief Private working DAG (null before setup()).
  FFGraph* dag
    ()
    const
    {
      return _dag;
    };

  //! @brief Declare the model constants, replacing previous ones, with optional values.
  //! @param C     constants
  //! @param valC  their values, one per constant, or empty (see also update_ref())
  void set_constant
    ( std::vector<FFVar> const& C, std::vector<double> const& valC=std::vector<double>() )
    {
      assert( valC.empty() || valC.size() == C.size() );
      auto& vCst    = _usr._vCstUsr;
      auto& vCstVal = _usr._vCstValUsr;
      vCst = C;
      vCstVal = valC;
      _issetup = false;
      _classified = false;
      _on_model_changed( ModelChange::DERIVATIVES );
    }

  //! @brief Declare the evolution (marching) direction used by PDE classification.
  //! Supersedes set_time_domain(), which remains as a deprecated wrapper.
  void set_evolution_domain
    ( FFVar const& dom )
    {
      auto& evolution_dom_var  = _usr._evolution_dom_varUsr;
      auto& evolution_dom_user = _usr._evolution_dom_userUsr;
      auto& evolution_dom_set  = _usr._evolution_dom_setUsr;
      evolution_dom_var  = dom;
      evolution_dom_user = true;
      evolution_dom_set  = true;
      _reresolve_auto_roles( _usr._mEqnUsr, &evolution_dom_var );   // AUTO: LB alone of t is INITIAL
      _issetup           = false;
      _classified        = false;
      _on_model_changed( ModelChange::DERIVATIVES );
    }

  //! @brief Remove the declared evolution direction (setup may then infer one).
  void reset_evolution_domain
    ()
    {
      auto& evolution_dom_user = _usr._evolution_dom_userUsr;
      auto& evolution_dom_set  = _usr._evolution_dom_setUsr;
      evolution_dom_user = false;
      evolution_dom_set  = false;
      _issetup           = false;
      _classified        = false;
      _on_model_changed( ModelChange::DERIVATIVES );
    }

  //! @brief Deprecated: use set_evolution_domain().
  void set_time_domain
    ( FFVar const& dom )
    { set_evolution_domain( dom ); }

  //! @brief Remove all constants and their values.
  void reset_constant
    ()
    {
      auto& vCst    = _usr._vCstUsr;
      auto& vCstVal = _usr._vCstValUsr;
      vCst.clear();
      vCstVal.clear();
      _issetup    = false;
      _classified = false;
      _on_model_changed( ModelChange::DERIVATIVES );
    }

  //! @brief Constants: the working copy after setup(), the declared ones before.
  std::vector<FFVar> const& var_constant
    ()
    const
    { return _issetup ? _vCst : _usr._vCstUsr; }

  //! @brief Constant values: the working copy after setup(), the declared ones before.
  std::vector<double> const& val_constant
    ()
    const
    { return _issetup ? _vCstVal : _usr._vCstValUsr; }

  //! @brief Declare several domains sharing one FFDom, with optional reference coordinates.
  //! @param vVar    domain variables
  //! @param optDom  extent and finite-element mesh, applied to every variable
  //! @param ref     reference coordinates for classification, one per variable, or empty
  void add_domain
    ( std::vector<FFVar> const& vVar, FFDom const& optDom,
      std::vector<double> const& ref=std::vector<double>() )
    {
      assert( ref.empty() || ref.size() == vVar.size() );
      for( size_t i=0; i<vVar.size(); ++i ){
        std::optional<double> optref;
        if( !ref.empty() ) optref = ref[i];
        add_domain( vVar[i], optDom, optref );
      }
    }

  //! @brief Declare a domain with its FFDom (extent and mesh), replacing any previous declaration.
  //! @param ref  reference coordinate for classification (default: the domain midpoint, see ref())
  void add_domain
    ( FFVar const& Var, FFDom const& optDom,
      std::optional<double> ref=std::nullopt )
    {
      auto& mDom = _usr._mDomUsr;
      auto& classDomRef = _usr._classDomRefUsr;
      mDom[Var] = optDom;
      if( ref ) classDomRef[Var] = *ref;
      _issetup    = false;
      _classified = false;
      _on_model_changed( ModelChange::DERIVATIVES );
    }

  //! @brief Remove all domains and their reference coordinates (states and equations are kept).
  void reset_domain
    ()
    {
      auto& mDom = _usr._mDomUsr;
      auto& classDomRef = _usr._classDomRefUsr;
      mDom.clear();
      classDomRef.clear();
      _issetup    = false;
      _classified = false;
      _on_model_changed( ModelChange::DERIVATIVES );
    }

  //! @brief True after a successful setup(), false before it and after any model change.  Also the value a
  //! copy inherits: a copy of a not-set-up environment is itself not set up (see OCFESLV's copy constructor).
  bool is_setup
    ()
    const
    { return _issetup; }

  //! @brief Domains: the working copy after setup(), the declared ones before.  Overrides OCBase's accessor
  //! through OCFESLV, so the collocation layer reads this very map.
  t_Dom const& var_domain
    ()
    const
    { return _issetup ? _mDom : _usr._mDomUsr; }

  //! @brief Inputs: the working copy after setup(), the declared ones before.  AFTER SETUP THE WORKING COPY ALSO
  //! HOLDS THE DEFERRED-VALUE INPUTS setup() created -- one per captured output, through which a consumer reads
  //! that value back (var_deferred() lists exactly those, and is_deferred_input() tests one).  A consumer that
  //! means "what the modeller declared" wants var_declared_input(), or it will take a captured value for a
  //! parameter of the model.
  t_Var const& var_input
    ()
    const
    { return _issetup ? _mInp : _usr._mInpUsr; }

  //! @brief The inputs the MODELLER declared, before and after setup() alike: the deferred-value inputs are not
  //! among them.  Use this where "input" means something given from outside -- a parameter, a control profile.
  t_Var const& var_declared_input
    ()
    const
    { return _usr._mInpUsr; }

  //! @brief Whether @p var is a deferred-value input setup() created, rather than one the modeller declared.
  bool is_deferred_input
    ( FFVar const& var )
    const
    { for( auto const& C : _deferredValue )
        if( C.input.dag() && var.dag() && C.input.id() == var.id() ) return true;
      return false; }

  //! @brief Declare a state over the domains @p vDom, with an optional scalar reference value.
  //! Re-declaring a state replaces its domains.  A scalar reference replaces any reference
  //! function set for the state.
  void add_state
    ( FFVar const& Var, std::vector<FFVar> const& vDom,
      std::optional<double> ref=std::nullopt )
    {
      auto& mVar = _usr._mVarUsr;
      auto& classStateRef = _usr._classVarRefUsr;
      auto& classStateRefFun = _usr._classVarRefFunUsr;
      auto [it,ins] = mVar.insert( {Var,{}} );
      if( !ins ) it->second.clear();
      for( auto const& d : vDom )
        it->second.insert( d );
      if( ref ){
        classStateRef[Var] = *ref;
        classStateRefFun.erase( Var );
      }
      _issetup    = false;
      _classified = false;
      _on_model_changed( ModelChange::DERIVATIVES );
    }

  //! @brief Declare a state with a reference function of the domain coordinates (see update_ref()).
  void add_state
    ( FFVar const& Var, std::vector<FFVar> const& vDom, t_Fun const& ref )
    { add_state( Var, vDom ); update_ref( Var, ref ); }

  //! @brief Declare several states over the same domains, with optional scalar references (one per
  //! state, or empty).
  void add_state
    ( std::vector<FFVar> const& vVar, std::vector<FFVar> const& vDom,
      std::vector<double> const& ref=std::vector<double>() )
    {
      assert( ref.empty() || ref.size() == vVar.size() );
      for( size_t i=0; i<vVar.size(); ++i ){
        std::optional<double> optref;
        if( !ref.empty() ) optref = ref[i];
        add_state( vVar[i], vDom, optref );
      }
    }

  //! @brief Remove all states and their reference values and functions.
  void reset_state
    ()
    {
      auto& mVar = _usr._mVarUsr;
      auto& classStateRef = _usr._classVarRefUsr;
      auto& classStateRefFun = _usr._classVarRefFunUsr;
      mVar.clear();
      classStateRef.clear();
      classStateRefFun.clear();
      _issetup    = false;
      _classified = false;
      _on_model_changed( ModelChange::DERIVATIVES );
    }

  //! @brief States: the working copy after setup(), the declared ones before.
  t_Var const& var_state
    ()
    const
    { return _issetup ? _mVar : _usr._mVarUsr; }

  //! @brief Declare an equation (residual = 0) over the domains @p vDom.
  //! @param Eqn     residual expression
  //! @param vDom    domain variables
  //! @param vLim    domain mask per variable (FFDom::ALL, LB, UB, or combinations such as ALL-LB);
  //!                missing masks default to ALL.  The sizes are checked only when
  //!                CRONOS__FFCOLLOC_CHECK is defined.
  //! @param EqnOpt  options, normalised by _normalise_equation_options()
  virtual void add_equation
    ( FFVar const& Eqn, std::vector<FFVar> const& vDom,
      std::vector<int> const& vLim, EqnOptions const& EqnOpt=EqnOptions() )
    {
#ifdef CRONOS__FFCOLLOC_CHECK
      assert( vLim.empty() || vDom.size() == vLim.size() );
#endif
      auto& mEqn = _usr._mEqnUsr;
      t_EqnDom md;
      for( size_t i=0; i<vDom.size(); ++i )
        md.insert( { vDom[i], i < vLim.size() ? vLim[i] : FFDom::ALL } );
      mEqn.push_back( { Eqn, md, _normalised_options( EqnOpt, md ) } );
      _issetup    = false; 
      _classified = false;
      _on_model_changed( ModelChange::DERIVATIVES );
    }

  //! @brief Declare several equations over the same domains and masks, with AUTO options.
  virtual void add_equation
    ( std::vector<FFVar> const& vEqn, std::vector<FFVar> const& vDom={},
      std::vector<int> const& vLim={} )
    {
      for( auto const& e : vEqn ) add_equation( e, vDom, vLim );
    }

  //! @brief Declare several equations over the same domains and masks, with the same options.
  virtual void add_equation
    ( std::vector<FFVar> const& vEqn, std::vector<FFVar> const& vDom,
      std::vector<int> const& vLim, EqnOptions const& opt )
    {
      for( auto const& e : vEqn ) add_equation( e, vDom, vLim, opt );
    }


  //! @brief Remove all equations.
  void reset_equation
    ()
    {
      auto& mEqn = _usr._mEqnUsr;
      mEqn.clear();
      _issetup    = false;
      _classified = false;
      _on_model_changed( ModelChange::DERIVATIVES );
    }

  //! @brief Equations: the working copy after setup(), the declared ones before.
  t_Eqns const& var_equation
    ()
    const
    { return _issetup ? _mEqn : _usr._mEqnUsr; }

  //! @brief Declare a scalar output evaluated at one point.
  //! @param Fct   output expression
  //! @param vDom  domain variables fixed by the point
  //! @return the INDEX of the new output (the first one added, for the list form): outputs are numbered in declaration
  //!         order from 0, and the index is the output's position in the solvers' blk_fct(); it stays valid after setup()
  //! @param vVal  coordinates, one per domain variable, or empty: each coordinate then defaults
  //!              to the domain's reference coordinate (ref(): the midpoint unless set)
  //! @throw Exceptions::INDEX if vVal is neither empty nor of the size of vDom
  size_t add_output
    ( FFVar const& Fct, std::vector<FFVar> const& vDom={},
      std::vector<double> const& vVal=std::vector<double>() )
    {
#ifdef CRONOS__FFCOLLOC_CHECK
      assert( vVal.empty() || vDom.size() == vVal.size() );
#endif
      if( !vVal.empty() && vDom.size() != vVal.size() )
        throw Exceptions( Exceptions::INDEX );
      t_FctDom md;
      for( size_t i=0; i<vDom.size(); ++i )
        md.insert( { vDom[i], i < vVal.size()? vVal[i]: ref( vDom[i] ) } );
      return _append_output_point( Fct, md );
    }

  //! @brief Point output with an explicit SIDE per direction (FFDom::MINUS / FFDom::PLUS), for a coordinate where
  //! the value may differ on either side (an element boundary, a transition): e.g. add_output(f,{t},{T1},{FFDom::PLUS})
  size_t add_output
    ( FFVar const& Fct, std::vector<FFVar> const& vDom, std::vector<double> const& vVal,
      std::vector<int> const& vSide )
    {
      if( vVal.size() != vDom.size() || vSide.size() != vDom.size() ) throw Exceptions( Exceptions::INDEX );
      t_FctDom md;  std::map<FFVar,int,lt_FFVar> sd;
      for( size_t i=0; i<vDom.size(); ++i ){
        md.insert( { vDom[i], vVal[i] } );
        if( vSide[i] == FFDom::PLUS ) sd[ vDom[i] ] = FFDom::PLUS;
        else if( vSide[i] != FFDom::MINUS ) throw Exceptions( Exceptions::INDEX );
      }
      return _append_output_point( Fct, md, sd );
    }

  //! @brief Point output with brace-list coordinates, e.g. add_output(f,{t},{1.0}); disambiguates
  //! from the distributed overload taking masks.
  size_t add_output
    ( FFVar const& Fct, std::vector<FFVar> const& vDom,
      std::initializer_list<double> vVal )
    {
      return add_output( Fct, vDom, std::vector<double>( vVal ) );
    }

  //! @brief Declare a distributed output, evaluated at the discretisation nodes selected by masks.
  //! vDom and vLim are as in add_equation() (missing masks default to ALL).  The rows appear in
  //! the output array after evaluation; they are not equations.
  //! @throw Exceptions::INDEX if vLim is neither empty nor of the size of vDom
  size_t add_output
    ( FFVar const& Fct, std::vector<FFVar> const& vDom,
      std::vector<int> const& vLim )
    {
#ifdef CRONOS__FFCOLLOC_CHECK
      assert( vLim.empty() || vDom.size() == vLim.size() );
#endif
      if( !vLim.empty() && vDom.size() != vLim.size() )
        throw Exceptions( Exceptions::INDEX );
      t_EqnDom grid;
      for( size_t i=0; i<vDom.size(); ++i )
        grid.insert( { vDom[i], i < vLim.size()? vLim[i]: FFDom::ALL } );
      return _append_output_distributed( Fct, grid, t_FctDom() );
    }

  //! @brief Distributed output with brace-list masks, e.g. add_output(g,{t},{FFDom::ALL-FFDom::LB});
  //! disambiguates from the point overload.
  size_t add_output
    ( FFVar const& Fct, std::vector<FFVar> const& vDom,
      std::initializer_list<int> vLim )
    {
      return add_output( Fct, vDom, std::vector<int>( vLim ) );
    }

  //! @brief Declare several distributed outputs over the same domains and masks.
  size_t add_output
    ( std::vector<FFVar> const& vFct, std::vector<FFVar> const& vDom,
      std::vector<int> const& vLim )
    {
      size_t const first = _usr._mFctUsr.size();     // the index of the first output added
      for( auto const& f : vFct ) add_output( f, vDom, vLim );
      return first;
    }

  //! @brief Declare a distributed output over @p vGridDom (masks @p vLim) at the fixed coordinates
  //! @p vPointVal of the domains @p vPointDom.
  //! @throw Exceptions::INDEX on inconsistent sizes
  size_t add_output
    ( FFVar const& Fct, std::vector<FFVar> const& vGridDom,
      std::vector<int> const& vLim, std::vector<FFVar> const& vPointDom,
      std::vector<double> const& vPointVal )
    {
#ifdef CRONOS__FFCOLLOC_CHECK
      assert( vLim.empty() || vGridDom.size() == vLim.size() );
      assert( vPointDom.size() == vPointVal.size() );
#endif
      if( ( !vLim.empty() && vGridDom.size() != vLim.size() )
       || vPointDom.size() != vPointVal.size() )
        throw Exceptions( Exceptions::INDEX );

      t_EqnDom grid;
      for( size_t i=0; i<vGridDom.size(); ++i )
        grid.insert( { vGridDom[i], i < vLim.size()? vLim[i]: FFDom::ALL } );
      t_FctDom point;
      for( size_t i=0; i<vPointDom.size(); ++i )
        point.insert( { vPointDom[i], vPointVal[i] } );
      return _append_output_distributed( Fct, grid, point );
    }

  //! @brief Remove all outputs (classification is not invalidated: outputs do not enter it).
  void reset_output
    ()
    {
      auto& mFct   = _usr._mFctUsr;
      mFct.clear();
      _issetup  = false;
      _on_model_changed( ModelChange::DERIVATIVES );
    }

  //! @brief Declare a TRANSITION at @p tau on the evolution direction @p dom: left(tau^-) = right(tau^+),
  //! componentwise.  @p left is evaluated at tau^-, @p right at tau^+; both are pointwise expressions of states,
  //! inputs and constants.  States no transition mentions are continuous at tau.  Checked at setup().
  void add_transition
    ( std::vector<FFVar> const& left, std::vector<FFVar> const& right, FFVar const& dom, double const tau )
    {
      if( left.empty() || left.size() != right.size() ) throw Exceptions( Exceptions::INDEX );
      if( _usr._mDomUsr.find( dom ) == _usr._mDomUsr.end() ) throw Exceptions( Exceptions::INDEX );
      _usr._mTrnUsr.push_back( t_Transition{ left, right, dom, tau } );
      _on_model_changed( ModelChange::DERIVATIVES );
    }
  //! @brief Scalar form of add_transition
  void add_transition
    ( FFVar const& left, FFVar const& right, FFVar const& dom, double const tau )
    { add_transition( std::vector<FFVar>{ left }, std::vector<FFVar>{ right }, dom, tau ); }

  //! @brief Declare a TRANSITION through EVALUATIONS, with no domain or time given: every state in @p left and @p right
  //! appears inside an FFEval on the evolution direction, all at one point tau -- at tau^- (FFDom::MINUS, the default)
  //! on the left, at tau^+ (FFDom::PLUS) on the right.  E.g. add_transition( OpE( x, t, 0.4 ) + d, OpE( x, t, 0.4,
  //! FFDom::PLUS ) ).  tau and the direction are read from the evaluations at setup(), which checks all of the above.
  void add_transition
    ( std::vector<FFVar> const& left, std::vector<FFVar> const& right )
    {
      if( left.empty() || left.size() != right.size() ) throw Exceptions( Exceptions::INDEX );
      t_Transition tr;  tr.left = left;  tr.right = right;  tr.evaluated = true;
      _usr._mTrnUsr.push_back( tr );
      _on_model_changed( ModelChange::DERIVATIVES );
    }
  //! @brief Scalar form of add_transition through evaluations
  void add_transition
    ( FFVar const& left, FFVar const& right )
    { add_transition( std::vector<FFVar>{ left }, std::vector<FFVar>{ right } ); }

  //! @brief Validate the declared transitions (setup): the evolution direction, tau strictly inside it, pointwise
  //! expressions of declared variables.  Then REFUSE them if this model does not take transitions
  //! (_transitions_supported(): true for a bare FFModel, which only analyses the model, and for ODESLV and OCFESLV;
  //! false for any other derived class) -- a transition must never be silently ignored.
  //! @brief CAUSALITY REFUSAL for equations (2026-09-30).  A reduction along the EVOLUTION direction -- a value at a point
  //! of it (FFEval) or an integral over it (FFIntegral), alone or together with spatial directions -- may not appear in
  //! a user EQUATION, in any solver or solve mode:
  //!  - if its operand depends on a STATE, the value is a functional of the solution, which cannot feed the solve that
  //!    produces it (monolithic returned a silent 0 from the post-solve latch; marching failed or integrated one window);
  //!  - if its operand depends on INPUTS / constants only, the value is known before the solve, but evaluating it there
  //!    is not implemented yet (open item) -- today it read 0 or failed in the same ways.
  //! Outputs are unaffected: they defer such reductions correctly in both modes.  Transitions' internal rows are added
  //! after this check and are not user equations.
  bool _check_equation_reductions
    ()
    const
    {
      FFVar const& t = _usr._evolution_dom_varUsr;
      if( !_usr._evolution_dom_setUsr || !t.dag() ) return true;
      size_t ie = 0;
      for( auto const& e : _usr._mEqnUsr ){
        ++ie;
        if( !e.var.dag() ) continue;
        FFSubgraph sg = e.var.dag()->subgraph( 1, &e.var );
        for( auto const* op : sg.l_op ){
          if( !op ) continue;
          bool evo = false;  char const* kind = nullptr;
          if( auto const* pe = mc::type_cast<FFEval const>( op ) ){
            for( auto const& [dv, z0] : pe->Coord() ){ (void)z0; if( dv.id() == t.id() ) evo = true; }
            kind = "an evaluation at a point of";
          }
          else if( auto const* pi = mc::type_cast<FFIntegral const>( op ) ){
            for( auto const& [dv, ord] : pi->Indep().expr ){ (void)ord; if( dv.id() == t.id() ) evo = true; }
            kind = "an integral over";
          }
          if( !evo ) continue;
          bool on_state = false;
          for( auto const* vi : op->varin ){
            if( !vi || !vi->dag() ) continue;
            FFSubgraph s2 = vi->dag()->subgraph( 1, vi );
            for( auto const* o2 : s2.l_op )
              if( o2 && o2->type == FFOp::VAR && o2->varout[0] && _usr._mVarUsr.count( *o2->varout[0] ) ) on_state = true;
          }
          std::cerr << "FFModel::setup ** EQUATION REFUSED: equation #" << ie << " contains " << kind << " the evolution direction "
                    << t.name() << ( on_state
                       ? " of an expression that depends on a STATE -- the value is a functional of the solution and cannot feed"
                         " the solve that produces it (causality); compute it as an OUTPUT"
                       : " of an expression of INPUTS only -- evaluating it before the solve is not supported yet; compute it as"
                         " an OUTPUT, or precompute it" ) << std::endl;
          return false;
        }
      }
      return true;
    }

  bool _check_transitions
    ()
    {
      for( size_t k = 0; k < _usr._mTrnUsr.size(); ++k ){
        auto& trw = _usr._mTrnUsr[k];
        auto const& tr = trw;
        auto fail = [&]( std::string const& why ){
          std::cerr << "FFModel::setup ** TRANSITION " << k << " (tau = " << tr.tau << ") REFUSED: " << why << std::endl;
          return false; };
        if( tr.evaluated ){                                         // add_transition( left, right ): read tau and the
          if( !_usr._evolution_dom_setUsr ) return fail( "no evolution direction is set" );   // direction off the evaluations
          FFVar const evo = _usr._evolution_dom_varUsr;
          bool have_tau = false;  double tau = 0.;  std::string why;
          std::set<FFVar const*> seen;
          // walk an expression, stopping at evaluations: a state reached OUTSIDE one has no time in a transition
          std::function<bool( FFVar const&, int )> walk = [&]( FFVar const& v, int const want ) -> bool {
            if( !v.dag() ) return true;
            auto const* op = v.opdef().first;
            if( !op || !seen.insert( &v ).second ) return true;
            if( op->type == FFOp::VAR ){
              if( _usr._mVarUsr.count( v ) ){ why = "state " + v.name() + " appears OUTSIDE an evaluation -- in a transition declared"
                                                " without a time point every state is evaluated, at tau^- (left) or tau^+ (right)"; return false; }
              return true;
            }
            if( mc::type_cast<FFPartial const>( op ) || mc::type_cast<FFIntegral const>( op ) ){
              why = "a transition expression must be POINTWISE -- no derivative or integral operator";  return false; }
            if( auto const* e = mc::type_cast<FFEval const>( op ) ){
              auto const& C = e->Coord();
              if( C.size() != 1 || C.begin()->first.id() != evo.id() ){ why = "an evaluation must consume the evolution direction only"; return false; }
              double const c = C.begin()->second;
              if( have_tau && std::fabs( c - tau ) > 64.*std::numeric_limits<double>::epsilon()*std::max( 1., std::fabs( tau ) ) ){
                why = "its evaluations are at DIFFERENT points (" + std::to_string( tau ) + " and " + std::to_string( c ) + ")";  return false; }
              have_tau = true;  tau = c;
              int const sd = e->side( C.begin()->first );
              if( sd != want ){ why = want == FFDom::MINUS? "the LEFT side is evaluated at tau^- -- an evaluation there is FFDom::PLUS"
                                                           : "the RIGHT side is evaluated at tau^+ -- give its evaluations FFDom::PLUS";  return false; }
              for( auto const* vi : op->varin ){                     // inside: pointwise, no nested evaluation
                if( !vi || !vi->dag() ) continue;
                FFSubgraph sg = vi->dag()->subgraph( 1, vi );
                for( auto const* o2 : sg.l_op )
                  if( o2 && ( mc::type_cast<FFPartial const>( o2 ) || mc::type_cast<FFIntegral const>( o2 ) || mc::type_cast<FFEval const>( o2 ) ) ){
                    why = "inside an evaluation a transition expression must be POINTWISE -- no derivative, integral or nested evaluation";  return false; }
              }
              return true;
            }
            for( auto const* vi : op->varin ) if( vi && !walk( *vi, want ) ) return false;
            return true;
          };
          for( auto const& ex : tr.left )  if( !walk( ex, FFDom::MINUS ) ) return fail( why );
          seen.clear();
          for( auto const& ex : tr.right ) if( !walk( ex, FFDom::PLUS ) ) return fail( why );
          if( !have_tau ) return fail( "declared without a time point, it must contain at least one evaluation" );
          trw.dom = evo;  trw.tau = tau;                            // what the evaluations say
        }
        if( !_usr._evolution_dom_setUsr || tr.dom.id() != _usr._evolution_dom_varUsr.id() )
          return fail( "its domain is not the evolution direction" );
        auto const itd = _usr._mDomUsr.find( tr.dom );
        if( itd == _usr._mDomUsr.end() ) return fail( "its domain is not declared" );
        double const lo = itd->second.lo_dom, up = itd->second.up_dom;
        double const tol = 64. * std::numeric_limits<double>::epsilon() * std::max( { 1., std::fabs( lo ), std::fabs( up ) } );
        if( !( tr.tau > lo + tol && tr.tau < up - tol ) )
          return fail( "tau must lie strictly inside the evolution domain" );
        for( auto const* side : { &tr.left, &tr.right } )
          for( auto const& ex : *side ){
            if( !ex.dag() ) continue;                                  // a constant
            FFSubgraph sg = ex.dag()->subgraph( 1, &ex );
            for( auto const* op : sg.l_op ){
              if( !op ) continue;
              if( !tr.evaluated && ( mc::type_cast<FFPartial const>( op ) || mc::type_cast<FFIntegral const>( op ) || mc::type_cast<FFEval const>( op ) ) )
                return fail( "a transition expression must be POINTWISE -- no derivative, integral or evaluation operator" );
              if( op->type != FFOp::VAR || !op->varout[0] ) continue;
              FFVar const& v = *op->varout[0];
              bool const known = _usr._mVarUsr.count( v ) || _usr._mInpUsr.count( v ) || v.id() == tr.dom.id() || _usr._mDomUsr.count( v )
                || std::any_of( _usr._vCstUsr.begin(), _usr._vCstUsr.end(), [&]( FFVar const& c ){ return c.id() == v.id(); } );
              if( !known ) return fail( "variable " + v.name() + " is not a declared state, input or constant" );
            }
          }
      }
      if( !_usr._mTrnUsr.empty() ){                          // square + structurally determined jumps, per tau
        t_TransitionJumps J;  std::string why;
        auto is_state = [&]( FFVar const& v ){ return _usr._mVarUsr.count( v ) > 0; };
        if( !_transition_jumps( _usr._mTrnUsr, is_state, J, why ) ){
          std::cerr << "FFModel::setup ** TRANSITION REFUSED: " << why << std::endl;
          return false;
        }
        // A post-jump state must be DIFFERENTIATED along the evolution direction in some equation: an algebraic state
        // is re-determined after tau by its algebraic equations -- it cannot be jumped, and a transition on it would
        // be silently ineffective.  (2026-09-30)
        std::set<long> differentiated;
        for( auto const& e : _usr._mEqnUsr ){
          if( !e.var.dag() ) continue;
          FFSubgraph sg = e.var.dag()->subgraph( 1, &e.var );
          for( auto const* op : sg.l_op ){
            auto const* pp = op? mc::type_cast<FFPartial const>( op ): nullptr;
            if( !pp ) continue;
            bool along = false;  for( auto const& [dv, ord] : pp->Indep().expr ){ (void)ord; if( dv.id() == _usr._evolution_dom_varUsr.id() ) along = true; }
            if( !along ) continue;
            for( auto const* vi : op->varin ) if( vi && vi->dag() ){
              FFSubgraph s2 = vi->dag()->subgraph( 1, vi );
              for( auto const* o2 : s2.l_op ) if( o2 && o2->type == FFOp::VAR && o2->varout[0] && _usr._mVarUsr.count( *o2->varout[0] ) )
                differentiated.insert( o2->varout[0]->id().second );
            }
          }
        }
        for( auto const& g : J )
          for( auto const& x : g.jump )
            if( !differentiated.count( x.id().second ) ){
              std::cerr << "FFModel::setup ** TRANSITION REFUSED at tau = " << g.tau << ": " << x.name() << " is not differentiated along"
                           " the evolution direction -- an ALGEBRAIC state is re-determined after tau by its equations and cannot"
                           " be jumped (map the differential states it depends on instead)" << std::endl;
              return false;
            }
      }
      if( !_usr._mTrnUsr.empty() && !_transitions_supported() ){
        std::cerr << "FFModel::setup ** add_transition: " << _usr._mTrnUsr.size() << " transition(s) declared and valid, but"
                     " this solver does not implement transitions -- refused rather than ignored" << std::endl;
        return false;
      }
      return true;
    }

  //! @brief Group @p trn by tau and determine, per group, the post-jump states S+ and a matching of components to them
  //! (see t_TransitionJump).  Refuses -- false, with @p why -- a group that is not SQUARE (components != |S+|) or not
  //! STRUCTURALLY determined (no matching saturates the components): its jump would be under- or over-determined.
  //! @p is_state tells a state from any other variable.
  static bool _transition_jumps
    ( t_Transitions const& trn, std::function<bool( FFVar const& )> const& is_state,
      t_TransitionJumps& out, std::string& why )
    {
      out.clear();
      std::vector<size_t> order( trn.size() );  for( size_t k = 0; k < trn.size(); ++k ) order[k] = k;
      std::sort( order.begin(), order.end(), [&]( size_t a, size_t b ){ return trn[a].tau < trn[b].tau; } );
      auto states_of = [&]( FFVar const& ex ){                // the states an expression depends on, by id
        std::vector<FFVar> st;
        if( !ex.dag() ) return st;
        FFSubgraph sg = ex.dag()->subgraph( 1, &ex );
        for( auto const* op : sg.l_op )
          if( op && op->type == FFOp::VAR && op->varout[0] && is_state( *op->varout[0] ) ) st.push_back( *op->varout[0] );
        return st; };
      for( size_t i = 0; i < order.size(); ){
        double const tau = trn[order[i]].tau;
        double const tol = 64.*std::numeric_limits<double>::epsilon()*std::max( 1., std::fabs( tau ) );
        t_TransitionJump J;  J.tau = tau;
        std::vector< std::vector<size_t> > adj;               // component -> indices into S+
        std::vector<FFVar> Splus;
        auto index_of = [&]( FFVar const& v ){
          for( size_t k = 0; k < Splus.size(); ++k ) if( Splus[k].id() == v.id() ) return k;
          Splus.push_back( v );  return Splus.size()-1; };
        for( ; i < order.size() && std::fabs( trn[order[i]].tau - tau ) <= tol; ++i ){
          auto const& tr = trn[order[i]];
          for( size_t c = 0; c < tr.right.size(); ++c ){
            J.comp.push_back( order[i]*65536 + c );
            std::vector<size_t> a;  for( auto const& v : states_of( tr.right[c] ) ) a.push_back( index_of( v ) );
            adj.push_back( a );
          }
        }
        std::ostringstream os;  os << "transition(s) at tau = " << tau << ": ";
        if( J.comp.size() != Splus.size() ){
          os << J.comp.size() << " component(s) but " << Splus.size() << " post-jump state(s) on the right-hand sides"
                " -- a jump must determine exactly the states it breaks (one component per state read at tau^+)";
          why = os.str();  return false;
        }
        // maximum bipartite matching, components -> S+ (augmenting paths)
        std::vector<long> matchS( Splus.size(), -1 );
        std::function<bool( size_t, std::vector<bool>& )> augment = [&]( size_t c, std::vector<bool>& vis ) -> bool {
          for( size_t s : adj[c] ){ if( vis[s] ) continue; vis[s] = true;
            if( matchS[s] < 0 || augment( (size_t)matchS[s], vis ) ){ matchS[s] = (long)c; return true; } }
          return false; };
        for( size_t c = 0; c < adj.size(); ++c ){
          std::vector<bool> vis( Splus.size(), false );
          if( !augment( c, vis ) ){
            os << "the components cannot each be matched to a distinct post-jump state they depend on -- the jump is"
                  " STRUCTURALLY singular (e.g. two components involving only the same state)";
            why = os.str();  return false;
          }
        }
        std::vector<size_t> compOrdered( J.comp.size() );
        J.jump.resize( Splus.size() );
        for( size_t s = 0; s < Splus.size(); ++s ){ compOrdered[s] = J.comp[ matchS[s] ]; J.jump[s] = Splus[s]; }
        J.comp = compOrdered;
        out.push_back( J );
      }
      return true;
    }

  //! @brief Whether this solver imposes transitions through LIFTED AUXILIARIES (OCFESLV): each post-jump state x_i at
  //! tau gets an auxiliary w_i standing for x_i(tau^+), defined by the transition rows (LINK) as a reduce_order()
  //! auxiliary is, and the solver joins x_i to w_i across tau instead of to x_i(tau^-).  False: the solver handles
  //! transitions its own way (ODESLV: per-stage initial values), and _lower_transitions() does nothing.
  virtual bool _transitions_lifted
    ()
    const
    { return false; }

  //! @brief Whether a lifting solver can place a transition at @p tau (default: yes); false with the reason in @p why
  virtual bool _accept_transition_tau
    ( double const tau, std::string& why )
    const
    { (void)tau; (void)why; return true; }
  //! @brief Does the consumer validate the transition times even when it does NOT lift the transitions (marching: a
  //! transfer map applies them at the window seams, so a tau between two seams is never met and was silently
  //! ignored -- WORKPLAN 1.6, 2026-10-06)?  Default: no (ODESLV, a bare FFModel).
  virtual bool _validate_transition_tau_unlifted
    ()
    const
    { return false; }

  //! @brief One lifted post-jump state: x_i(tau^+) is the auxiliary aux (see _transitions_lifted)
  struct t_TransitionLift
  {
    FFVar  state;
    FFVar  aux;
    double tau = 0.;
  };
  //! @brief The lifted post-jump states of the working model (filled by _lower_transitions)
  std::vector<t_TransitionLift> _trnLift;

  //! @brief Lower the working transitions for a solver that lifts them (see _transitions_lifted): per tau, per post-jump
  //! state, mint the auxiliary w_i, register it with parent x_i, and append each transition component as a LINK row
  //!   left( x(tau^-) ) - right( w for the post-jump states, x(tau^+) for the others ) = 0.
  //! Runs AFTER order reduction: the rows read x(tau^-/+) in-solve through pinned evaluations, which the reduction
  //! would capture as post-solve latches.
  bool _lower_transitions
    ()
    {
      _trnLift.clear();
      if( _mTrn.empty() ) return true;
      if( !_transitions_lifted() ){
        // not lifted (marching): the transfer map applies the jumps at the window seams -- still validate the times
        if( !_validate_transition_tau_unlifted() ) return true;
        t_TransitionJumps Ju;  std::string whyu;
        auto is_state_u = [&]( FFVar const& v ){ return _mVar.count( v ) > 0; };
        if( !_transition_jumps( _mTrn, is_state_u, Ju, whyu ) ){
          std::cerr << "FFModel::setup ** TRANSITION REFUSED: " << whyu << std::endl;  return false; }
        for( auto const& G : Ju ){ std::string w2;
          if( !_accept_transition_tau( G.tau, w2 ) ){
            std::cerr << "FFModel::setup ** TRANSITION REFUSED at tau = " << G.tau << ": " << w2 << std::endl;  return false; } }
        return true;
      }
      t_TransitionJumps J;  std::string why;
      auto is_state = [&]( FFVar const& v ){ return _mVar.count( v ) > 0; };
      if( !_transition_jumps( _mTrn, is_state, J, why ) ){
        std::cerr << "FFModel::setup ** TRANSITION REFUSED: " << why << std::endl;  return false; }
      FFVar const t = _evolution_dom_var;
      for( auto const& G : J ){ std::string w2;
        if( !_accept_transition_tau( G.tau, w2 ) ){
          std::cerr << "FFModel::setup ** TRANSITION REFUSED at tau = " << G.tau << ": " << w2 << std::endl;  return false; } }
      FFEval OpE;
      std::vector<FFVar> evo_states;                         // the states that evolve (the only ones a transition reads)
      for( auto const& [x, dm] : _mVar ) if( dm.count( t ) ) evo_states.push_back( x );
      for( auto const& G : J ){
        std::map<long,FFVar> w_of;                             // state id -> its auxiliary at this tau
        for( auto const& x : G.jump ){
          std::vector<FFVar> wdom;  for( auto const& d : _mVar.at( x ) ) if( d.id() != t.id() ) wdom.push_back( d );
          // v1.5: a state that also spans spatial directions jumps as a PROFILE: w lives on those directions
          std::ostringstream nm;  nm << "Wtrn_" << x.name() << "@" << G.tau;
          FFVar const w = _dag->add_var( nm.str() );
          _add_working_state( w, wdom );
          _register_aux_def( w, x, wdom, std::set<FFVar,lt_FFVar>(), 0, x );
          w_of[ x.id().second ] = w;
          _trnLift.push_back( t_TransitionLift{ x, w, G.tau } );
        }
        std::vector<FFVar> tL, rL, tR, rR;
        for( auto const& x : evo_states ){
          tL.push_back( x );  rL.push_back( OpE( x, t, G.tau, FFDom::MINUS ) );
          tR.push_back( x );  auto const iw = w_of.find( x.id().second );
          rR.push_back( iw != w_of.end()? iw->second: OpE( x, t, G.tau, FFDom::PLUS ) );
        }
        for( size_t const k : G.comp ){
          auto const& tr = _mTrn[ k / 65536 ];  size_t const c = k % 65536;
          FFVar const L = _dag->compose( std::vector<FFVar>{ tr.left[c]  }, tL, rL )[0];
          FFVar const R = _dag->compose( std::vector<FFVar>{ tr.right[c] }, tR, rR )[0];
          // the row lives on the spatial directions of the state matched to this component (pointwise there; a scalar
          // state gives a domain-less row)
          std::map<FFVar,int,lt_FFVar> rdom;
          for( size_t q = 0; q < G.comp.size(); ++q ) if( G.comp[q] == k )
            for( auto const& d : _mVar.at( G.jump[q] ) ) if( d.id() != t.id() ) rdom[d] = FFDom::ALL;
          _mEqn.push_back( { L - R, rdom, _make_equation_options( _link_eqn_options( 0, true ) ) } );
        }
      }
      return _derive_transition_lifts();
    }

  //! @brief DERIVED transitions: an evolving AUXILIARY (reduce_order's, or a solver's materialised capture source) whose
  //! definition reads a quantity that jumps at tau jumps with it -- a(tau^+) = expr( x(tau^+) ).  Lift it the same way
  //! (its own w, its claim at tau joined to w, the LINK row w - expr(tau^+)), so that no claim demands continuity of a
  //! discontinuous quantity.  Repeated until nothing new is lifted (an auxiliary defined from another follows the
  //! chain); idempotent, so a solver calls it again after minting auxiliaries of its own.
  bool _derive_transition_lifts
    ()
    {
      if( _trnLift.empty() ) return true;
      FFVar const t = _evolution_dom_var;
      FFEval OpE;
      std::vector<double> taus;
      for( auto const& L : _trnLift ){
        bool seen = false;  for( double const u : taus ) if( std::fabs( u - L.tau ) <= 64.*DBL_EPSILON*std::max( 1., std::fabs( u ) ) ) seen = true;
        if( !seen ) taus.push_back( L.tau );
      }
      for( double const tau : taus ){
        auto same = [&]( double u ){ return std::fabs( u - tau ) <= 64.*DBL_EPSILON*std::max( 1., std::fabs( tau ) ); };
        for( bool grew = true; grew; ){
          grew = false;
          std::map<long,FFVar> w_of;
          for( auto const& L : _trnLift ) if( same( L.tau ) ) w_of[ L.state.id().second ] = L.aux;
          for( auto const& def : _auxDef ){
            auto const itv = _mVar.find( def.aux );
            if( itv == _mVar.end() || !itv->second.count( t ) || w_of.count( def.aux.id().second ) || !def.expr.dag() ) continue;
            if( def.aux.id() == def.expr.id() ) continue;
            bool reads_jump = false, pointwise = true;
            FFSubgraph sg = _dag->subgraph( 1, &def.expr );
            for( auto const* op : sg.l_op ){
              if( !op ) continue;
              // a derivative, integral or evaluation ACROSS the jump (along spatial directions only) carries over to the
              // post-jump profile -- the substitution x -> w rebuilds d x/dz as d w/dz, OpE(x,z,z0) as OpE(w,z,z0) -- but
              // one ALONG the evolution direction does not (the x-dot(tau+) family)
              if( auto const* pp = mc::type_cast<FFPartial const>( op ) )
                for( auto const& [dv, ord] : pp->Indep().expr ){ (void)ord; if( dv.id() == t.id() ) pointwise = false; }
              if( auto const* pi = mc::type_cast<FFIntegral const>( op ) )
                for( auto const& [dv, ord] : pi->Indep().expr ){ (void)ord; if( dv.id() == t.id() ) pointwise = false; }
              if( auto const* pe = mc::type_cast<FFEval const>( op ) )
                for( auto const& [dv, z0] : pe->Coord() ){ (void)z0; if( dv.id() == t.id() ) pointwise = false; }
              if( op->type == FFOp::VAR && op->varout[0] && w_of.count( op->varout[0]->id().second ) ) reads_jump = true;
            }
            if( !reads_jump ) continue;
            if( !pointwise ){
              std::cerr << "FFModel::setup ** TRANSITION REFUSED at tau = " << tau << ": auxiliary " << def.aux.name() << " is defined"
                           " by a derivative along the evolution direction, an integral or an evaluation of a jumping state -- not"
                           " derivable at tau+ yet" << std::endl;
              return false; }
            std::vector<FFVar> wdom;  for( auto const& d : itv->second ) if( d.id() != t.id() ) wdom.push_back( d );
            std::ostringstream nm;  nm << "Wtrn_" << def.aux.name() << "@" << tau;
            FFVar const w = _dag->add_var( nm.str() );
            FFVar const aux = def.aux, expr = def.expr;  int const blk = def.block_id;   // copies: _auxDef grows below
            std::vector<FFVar> tA, rA;
            for( auto const& [x, dm] : _mVar ){
              if( !dm.count( t ) || x.id() == aux.id() ) continue;
              auto const iw = w_of.find( x.id().second );
              tA.push_back( x );  rA.push_back( iw != w_of.end()? iw->second: OpE( x, t, tau, FFDom::PLUS ) );
            }
            FFVar const E = _dag->compose( std::vector<FFVar>{ expr }, tA, rA )[0];
            _add_working_state( w, wdom );
            _register_aux_def( w, aux, wdom, std::set<FFVar,lt_FFVar>(), blk, aux );
            _trnLift.push_back( t_TransitionLift{ aux, w, tau } );
            std::map<FFVar,int,lt_FFVar> rdom;  for( auto const& d : wdom ) rdom[d] = FFDom::ALL;
            _mEqn.push_back( { w - E, rdom, _make_equation_options( _link_eqn_options( blk, true ) ) } );
            grew = true;  break;                               // _auxDef and _mVar grew: restart the scan
          }
        }
      }
      return true;
    }

  //! @brief The stream for an informational message of display level @p n (DISPLAY_REVIEW_20261001): std::cout when
  //! options.DISPLAY_LEVEL >= n, a null stream otherwise.  Refusals and failures go to std::cerr at every level.
  std::ostream& _disp
    ( int const n )
    const
    { static std::ostream null( nullptr ); return options.DISPLAY_LEVEL >= n? std::cout: null; }

  //! @brief Whether this model takes transitions (add_transition); false: setup() refuses them, never ignores them.
  //! A bare FFModel only ANALYSES the model, so it takes them (validated, and listed in the report); a class derived
  //! from FFModel must override this to take them (ODESLV and OCFESLV do), so a solver that does not implement
  //! transitions still refuses them.
  virtual bool _transitions_supported
    ()
    const
    { return typeid( *this ) == typeid( FFModel ); }

  //! @brief The declared transitions
  t_Transitions const& var_transition
    ()
    const
    { return _issetup ? _mTrn : _usr._mTrnUsr; }   // like var_output: the working copy once set up

  //! @brief Outputs: the working copy after setup(), the declared ones before.
  t_Fcts const& var_output
    ()
    const
    { return _issetup ? _mFct : _usr._mFctUsr; }

  //! @brief Write a report of the model to @p os: the domains; the states and inputs with their domains (an input
  //! holding a deferred value is marked); the constants; the equations with role, block and masks; the outputs; the
  //! deferred values; the classification per block, with the evolution direction (set or found) and the index
  //! character; and the degree-of-freedom balance.  Before setup() this is the DECLARED model, after it the working
  //! one -- the auxiliaries, LINK rows and closures setup() added included.  A derived solver has its own display
  //! (OCFESLV's operator<<, which adds the collocated counts); this one needs no discretisation.
  void report
    ( std::ostream& os )
    const;

  //! @brief Write report() to @p os.
  friend std::ostream& operator<<
    ( std::ostream& os, FFModel const& mod )
    { mod.report( os ); return os; }

  //! @brief The account of the model's initial data after the last setup() (see t_InitialData).
  t_InitialData const& initial_data
    ()
    const
    { return _initialData; }

  //! @brief The model's degree-of-freedom balance, as computed by the last setup() (see t_DofBalance).
  t_DofBalance const& dof_balance
    ()
    const
    { return _dofBalance; }

  //! @brief The deferred values setup() created -- the functionals of the solution a consumer must evaluate after
  //! solving and write into their inputs (see t_DeferredValue for the contract).  Empty before setup(), and for a
  //! model whose reductions are all spatial.
  std::vector<t_DeferredValue> const& var_deferred
    ()
    const
    { return _deferredValue; }


public:

  //! @brief Set the scalar reference of a declared state, input or domain, or the value of a constant.
  //! Looked up in that order.  For a state or input it replaces any reference function; for a
  //! constant it sets the value (values default to 0 when none were given); for a domain it sets
  //! the reference coordinate.  Invalidates setup and classification.
  //! @throw Exceptions::INDEX if @p var is not declared; Exceptions::CSTVAL on inconsistent
  //! constant values
  void update_ref
    ( FFVar const& var, double const& val );

  //! @brief Set the reference function of a declared state or input; it replaces any scalar reference.
  //! @throw Exceptions::INDEX if @p var is not a declared state or input
  void update_ref
    ( FFVar const& var, t_Fun const& fun );

  //! @brief Effective reference value of a state, input, constant or domain.
  //! State or input: its reference function evaluated at the domains' reference coordinates, else
  //! its scalar reference, else 0.  Constant: its value (0 if none).  Domain: its reference
  //! coordinate, else its midpoint.  Looked up in the working model then the declared one after
  //! setup(), in the reverse order before.
  //! @throw Exceptions::INDEX if @p var is unknown; Exceptions::CSTVAL on inconsistent constant values
  double ref
    ( FFVar const& var )
    const;

protected:

  //! @brief Append a point-output record to the declared outputs (sizes checked by the caller).
  size_t _append_output_point
    ( FFVar const& Fct, t_FctDom const& point,
      std::map<FFVar,int,lt_FFVar> const& side = std::map<FFVar,int,lt_FFVar>() )
    {
      t_Fct rec;
      rec.var   = Fct;
      rec.kind  = FctKind::POINT;
      rec.point = point;
      rec.side  = side;
      rec.grid.clear();
      rec.row0 = 0;
      rec.nrow = 1;
      size_t const ndx = _usr._mFctUsr.size();       // outputs are numbered in declaration order, from 0
      _usr._mFctUsr.push_back( std::move( rec ) );
      _issetup  = false;
      _on_model_changed( ModelChange::DERIVATIVES );
      return ndx;
    }

  //! @brief Append a distributed-output record to the declared outputs (sizes checked by the caller).
  size_t _append_output_distributed
    ( FFVar const& Fct, t_EqnDom const& grid, t_FctDom const& point=t_FctDom() )
    {
      t_Fct rec;
      rec.var   = Fct;
      rec.kind  = FctKind::DISTRIBUTED;
      rec.point = point;
      rec.grid  = grid;
      rec.row0 = 0;
      rec.nrow = 0;
      size_t const ndx = _usr._mFctUsr.size();       // outputs are numbered in declaration order, from 0
      _usr._mFctUsr.push_back( std::move( rec ) );
      _issetup  = false;
      _on_model_changed( ModelChange::DERIVATIVES );
      return ndx;
    }

  //! @brief True once classify_pde() has succeeded on the current model; cleared by model mutations
  //! other than output changes.
  bool                                      _classified = false;

  //! @brief Scalar references of states (working copy; declared values are in _usr).
  std::map<FFVar,double,lt_FFVar>           _classVarRef;

  //! @brief Scalar references of inputs (working copy).
  std::map<FFVar,double,lt_FFVar>           _classInpRef;

  //! @brief Reference coordinates of domains (working copy).
  std::map<FFVar,double,lt_FFVar>           _classDomRef;

  //! @brief Reference functions of states (working copy).
  std::map<FFVar,t_Fun,lt_FFVar>         _classVarRefFun;

  //! @brief Reference functions of inputs (working copy).
  std::map<FFVar,t_Fun,lt_FFVar>         _classInpRefFun;

  //! @brief Resolve and constrain equation options given the equation's domain masks.
  //! AUTO becomes BOUNDARY if some mask is LB or UB alone, INTERIOR otherwise.  INITIAL, BOUNDARY,
  //! INTERFACE and DIAGNOSTIC equations do not enter classification; DIAGNOSTIC equations receive
  //! no SAT and are never donated; LINK equations are never donated.  Other flags are kept.
  //! IN PLACE and VIRTUAL.  The base applies the MODEL rules (resolve AUTO, classification); a solver
  //! override calls it and then applies its own rules to its own fields.  In place is safe because every
  //! environment owns its options (EqnOptions::clone).
  virtual void _normalise_equation_options
    ( EqnOptions& opt, std::map<FFVar,int,lt_FFVar> const& eqndom )
    const;

  //! @brief the FACTORY for every stored options object -- user-added and internally minted alike.
  //! The base clones the seed (preserving its dynamic type); a solver overrides it to promote a base seed to
  //! its own derived type, so that no stored equation ever lacks the solver's fields.
  virtual std::shared_ptr<EqnOptions> _make_equation_options
    ( EqnOptions const& seed )
    const
    { return seed.clone(); }

  //! @brief build through the factory, then normalise in place.
  std::shared_ptr<EqnOptions> _normalised_options
    ( EqnOptions const& seed, std::map<FFVar,int,lt_FFVar> const& eqndom )
    const
    { auto p = _make_equation_options( seed ); _normalise_equation_options( *p, eqndom ); return p; }

  //! @brief the options of a LINK row minted by reduce_order().  Classification is the model's
  //! decision; SAT reception and donation are the solver's, supplied by its factory's role defaults.
  static EqnOptions _link_eqn_options
    ( int block, bool classify )
    { EqnOptions o( EqnRole::LINK, block ); o.participate_in_classification = classify; return o; }

  //! @brief True for the trace roles INITIAL, BOUNDARY and INTERFACE.
  //! @brief The role AUTO resolves to for an equation on @p eqndom (2026-10-01): INITIAL if its mask on the
  //! evolution domain @p evol is LB alone; else BOUNDARY if any mask is LB or UB alone; else INTERIOR.
  static EqnRole _auto_role
    ( std::map<FFVar,int,lt_FFVar> const& eqndom, FFVar const* evol )
    {
      if( evol )
        for( auto const& [v,lim] : eqndom )
          if( v.id() == evol->id() && lim == FFDom::LB ) return EqnRole::INITIAL;
      for( auto const& [v,lim] : eqndom )
        if( lim == FFDom::LB || lim == FFDom::UB ) return EqnRole::BOUNDARY;
      return EqnRole::INTERIOR;
    }

  //! @brief Re-resolve the AUTO-derived roles of @p eqns for the evolution domain @p evol (nullptr: none): only
  //! between BOUNDARY and INITIAL, which classification treats alike (neither participates).
  static void _reresolve_auto_roles
    ( t_Eqns& eqns, FFVar const* evol )
    {
      for( auto& eqn : eqns ){
        if( !eqn.opt || !eqn.opt->role_auto ) continue;
        EqnRole const r = eqn.opt->role;
        if( r != EqnRole::BOUNDARY && r != EqnRole::INITIAL ) continue;
        EqnRole const n = _auto_role( eqn.dom, evol );
        if( n == EqnRole::BOUNDARY || n == EqnRole::INITIAL ) eqn.opt->role = n;
      }
    }

  static bool _ordinary_trace_role
    ( EqnRole role )
    { return role == EqnRole::INITIAL || role == EqnRole::BOUNDARY
          || role == EqnRole::INTERFACE; }

  //! @brief Clear everything a re-setup must discard.  FFModel clears the working model; a derived solver
  //! overrides this to clear its own state as well and calls this version.
  virtual void _reset
    ()
    { _reset_model(); }

  //! @brief Delete the working DAG if owned; the user DAG is never touched.
  void _release_working_dag
    ()
    {
      if( _dagOwned ){
        delete _dag;
        _dag = nullptr;
      }
      _dagOwned = false;
    }

  //! @brief Clear the working model (working collections, classification references) and the setup and
  //! classification flags.  The declared model and the working DAG are kept.
  void _reset_model
    ();

  //! @brief Replace this model by a copy of @p src: options, declared model, a new working DAG, and the
  //! working collections remapped to it.  The DAG is imported with ONE FFGraph::insert of the model roots
  //! (constants, domains, states, inputs, equations, outputs) followed by @p extraRoots, a derived
  //! solver's own roots.  insert() keeps the ids of independent variables but re-creates auxiliary nodes by
  //! evaluation (new ids), so members are remapped through the positional index it returns:
  //! @p ndxAllVars maps a source id to a position in @p vAllVarsLoc, returned for the solver's own remap.
  //! Setup and classification flags are left to the caller.  Returns false if the import fails.
  bool _deep_copy_model
    ( FFModel const& src, std::vector<FFVar> const& extraRoots,
      std::vector<FFVar>& vAllVarsLoc, std::map<FFVar::pt_idVar,size_t>& ndxAllVars );

public:

  //! @brief Equation types based on principal symbol analysis.
  //!
  //! Use the flags in t_Classify, especially evolution_hyperbolic
  //! and spatially_characteristic, to distinguish causal/evolution
  //! hyperbolicity from a merely characteristic spatial symbol.
  enum EqnType
  {
    DIFFERENTIAL_ORDINARY = 0,    //!< Pure initial/marching problem with one effective domain
    DIFFERENTIAL_ALGEBRAIC,       //!< Singular evolution matrix and no transverse domain
    ALGEBRAIC_FIELD,              //!< DISTRIBUTED states, but NO derivative of any state in any
                                  //!< equation: the principal symbol is EMPTY, not inconclusive.
                                  //!< Distinct from UNDETERMINED, which means "could not tell".
                                  //!< Consequences differ: an algebraic field needs NO initial
                                  //!< condition, has nothing to march, and has no characteristic
                                  //!< treatment.  Continuity across element seams is IMPLIED by
                                  //!< the pointwise closure, so the interface plan correctly
                                  //!< emits no claims (measured on OCFE_ALGFIELD1: nTau=0 in all
                                  //!< three imposition modes, results bit-identical, |y-y*|
                                  //!< converging at order ~7).
    DIFFERENTIAL_IMPLICIT,        //!< M * x' = f with M NONSINGULAR: an ODE of index 0, and collocation solves it
                                  //!< as it stands, but the derivatives are coupled through that matrix, so an
                                  //!< explicit integrator needs it inverted while IDAS takes the residual as is.
                                  //!< Distinguished from DIFFERENTIAL_ORDINARY by a STRUCTURAL test on the
                                  //!< evolution coefficient matrix -- one nonzero per row AND per column -- not by
                                  //!< diagonality, which would misjudge a scaled or permuted assignment.
    RECTANGULAR_RETIRED_,         //!< sentinel for the retired DIFFERENTIAL_RECTANGULAR.  It
                                  //!< was a SHAPE fact written into the CHARACTER axis, and it
                                  //!< is now carried by t_Classify::symbol_rectangular, with
                                  //!< symbol_analysed (D3) distinguishing "no square test
                                  //!< applies" from "analysis inconclusive".  Kept as a named
                                  //!< sentinel rather than deleted outright so a stale integer
                                  //!< comparison fails loudly.  ORIGINAL DOC: the principal
                                  //!< symbol is RECTANGULAR because the retained
                                  //!< differential rows differentiate MORE states than there are
                                  //!< rows -- one row carrying several derivatives (an ALE
                                  //!< equation with q'xi, u't and x't, say) is enough.  The FULL
                                  //!< system is square; only the differential core is not.
                                  //!< Distinct from UNDETERMINED, which means "analysis
                                  //!< inconclusive": here the shape is KNOWN and it is known to
                                  //!< be uncharacterisable by a square-symbol test.
                                  //!< MEASURED on OCFE_MMPDE27 (--evolve --gdirect): sym=3x4,
                                  //!< retained rows differentiate u,q,x,g while Z26 alone carries
                                  //!< three of them.  Routing definitional rows to vAlgEqn is
                                  //!< CORRECT here and still cannot yield a square core.
    ALGEBRAIC_LUMPED,             //!< ZERO collocation domains: every state is a lumped scalar and
                                  //!< the system is a plain square F(x)=0 (OCFE_AE0/AE1/AE2).
                                  //!< Named rather than skipped so the verdict is reportable.
    EVOL_HYPERBOLIC,              //!< Evolution direction known/detected; real characteristic speeds
    PARABOLIC,                    //!< Descriptor system with recognised diffusion/link structure
    ELLIPTIC,                     //!< No real characteristics detected in the sampled principal symbol
    COMPLEX_CHARACTERISTIC,       //!< The sampled principal symbol has COMPLEX characteristic
                                  //!< speeds (max |Im(lambda)| > imag_tol) and the evolution
                                  //!< sweep did not succeed: the block is not hyperbolic in the
                                  //!< sampled directions.
                                  //!< This value was once called MIXED and
                                  //!< documented as "classification changes with sampled
                                  //!< direction" -- a DIFFERENT condition, tested nowhere.  The
                                  //!< assignment has always been `max_imag_eig > imag_tol`, so
                                  //!< the name now matches the test.  MEASURED: 0 occurrences in
                                  //!< 1410 corpus classification lines; OCFE_SYMBOLCASES1's
                                  //!< COMPLEX case (A = [[0,-1],[1,0]], eigenvalues +-i) reaches
                                  //!< it and reports evolution_hyperbolic=0, max_imag_eig=1.0.
    UNDETERMINED,                 //!< Non-square system or analysis inconclusive
    SPATIALLY_CHARACTERISTIC,     //!< Real characteristic directions, but no evolution direction
    DESCRIPTOR,                   //!< Singular evolution matrix, not recognised as parabolic
    DEGENERATE,                   //!< Principal symbol singular in all sampled directions
    DEGENERATE_LAST_              //!< sentinel for the retired WEAK_HYPERBOLIC (it was a
                                  //!< qualifier written into the character axis).  Weak
                                  //!< hyperbolicity now lives ONLY in t_Classify::weak_hyperbolic,
                                  //!< which every consumer already read.  Kept as a named
                                  //!< sentinel rather than deleted outright so that any stale
                                  //!< integer comparison fails loudly at compile time.
  };


protected:

  //! @brief vector of domain variables
  std::vector<FFVar>                 _vDom;

  //! @brief Fast identity lookup for auxiliary state variables.
  std::set<FFVar::pt_idVar>          _auxVarID;

  //! @brief vector of inputs
  std::vector<FFVar>                 _vInp;

  //! @brief vector of states
  std::vector<FFVar>                 _vVar;

  //! @brief Symbolic principal symbol of a first-order PDE system.
  //! Populated by _principal_symbol(); all FFVar fields are live DAG nodes.
  struct t_Symbol
  {
    //! Interior equation root nodes (rows of each coefficient matrix)
    std::vector<FFVar>                vEqn;
    //! Participating equations with NO differentiated state (algebraic closures/constraints, e.g.
    //! omega = g(...)).  Excluded from the differential symbol (they would be zero rows); their
    //! bare states form the algebraic block handled by the index analysis.
    std::vector<FFVar>                vAlgEqn;
    //! State variables appearing in at least one derivative term (columns)
    std::vector<FFVar>                vState;
    //! Domain variables in map order (one coefficient matrix each)
    std::vector<FFVar>                vDom;
    //! Symbolic coefficient matrices, row-major, one entry per domain variable:
    //! vCoeff[i][ k * vState.size() + j ] = dF_k / d(du_j/dx_i)  as FFVar
    std::vector< std::vector<FFVar> > vCoeff;
    //! Zeroth-order state-Jacobian coefficients from the SAME proxied system
    //! (derivative terms held as independent proxies), row-major:
    //! vCoeff0[ k * vState.size() + j ] = dF_k / d(u_j-value)  as FFVar.
    //! This captures the lower-order algebraic coupling (e.g. advection) that
    //! the principal symbol omits.  Populated by _principal_symbol() and consumed
    //! by _compute_interface_reads() (item 11 Part-2).
    std::vector<FFVar>                vCoeff0;
    //! the SAME prepared zeroth-order Jacobian, taken over the ALGEBRAIC rows
    //! (vAlgEqn) and over ALL states (vAlgState, not just vState -- the states that need it
    //! are precisely the ones absent from vState).  vAlgCoeff0[ g*vAlgState.size() + j ]
    //! = d(vAlgEqn[g]) / d(state_j-value), derivatives held as proxies, nonlocal terms
    //! stripped, exactly as for vCoeff0.  Empty if the computation was not possible.
    std::vector<FFVar>                vAlgState;
    std::vector<FFVar>                vAlgCoeff0;

    //! @brief C2's OWN zeroth-order coefficients -- all participating rows (differential and algebraic)
    //! over all states, with every FFPartial output proxied and nothing else substituted and no stripping, so a
    //! state that appears BARE in a row is differentiated as itself.  The symbol's own tables are unchanged;
    //! this one exists because C2 needs the row's coefficient, not the symbol's view of it.
    std::vector<FFVar>                vC2Eqn;
    std::vector<FFVar>                vC2State;
    std::vector<FFVar>                vC2Coef0;
    //! vCoeff0 over ALL states (vAllState) for the DIFFERENTIAL rows, so a
    //! differential receiver's zeroth-order coefficient of an auxiliary outside vState is
    //! available too.  vCoeff0All[ k*vAllState.size() + j ].  Same prepared FAD as vCoeff0.
    std::vector<FFVar>                vAllState;
    std::vector<FFVar>                vCoeff0All;
    //! Derivative proxy VARs created by _principal_symbol (flat order
    //! proxy[i*nState+j] = d(state j)/d(dom i)).  vCoeff0 entries reference these
    //! when a derivative's coefficient is state-dependent; they must be supplied
    //! to any eval of vCoeff0.  Their value is immaterial to the magnitude read,
    //! which consumes only aux columns of A0 (proxy-free).  (item 11 Part-2)
    std::vector<FFVar>                vDerivProxy;
  };

  //! @brief Return the cached principal symbol associated with a block id.
  t_Symbol const& _symbol_for_block
    ( int block_id ) const;

  //! @brief Cached principal symbol (populated by classify_pde)
  t_Symbol                                  _symbol;

  //! @brief Block-level principal symbols, classifications and face data.
  std::map<int,t_Symbol>                    _blockSymbol;

  //! @brief Strip lower-order nonlocal terms before principal-symbol FAD.
  std::vector<FFVar> _strip_nonlocal_terms_for_symbol
    ( std::vector<FFVar> const& eqns )
    const;


public:

  //! @brief Human-readable PDE type name.
  static const char* pde_type_name
    ( EqnType type );

  //! @brief Return cached principal symbol (valid after setup())
  t_Symbol const& symbol_cached
    ()
    const
    { return _symbol; }


protected:

  //! @brief Extract the symbolic principal symbol of the first-order system.
  //! Must be called after reduce_order() so that all equations are first-order.
  t_Symbol _principal_symbol
    ()
    const;

  //! @brief Build principal symbol for a specific equation block.
  t_Symbol _principal_symbol
    ( int const block_id )
    const;

  //! @brief Evaluate the coefficient matrices of the principal symbol at a
  //! given reference point.
  //! @param sym        Output of _principal_symbol()
  //! @param state_vals Values for sym.vState  (length = vState.size())
  //! @param cst_vals   Values for _vCst constants (length = _vCst.size())
  //! @param dom_vals   Values for sym.vDom domain variables (length = vDom.size())
  //! @param input_vals Values for _vInp inputs (length = _vInp.size()).
  //!                   Distributed inputs are treated as spatially constant.
  std::vector<arma::mat> _eval_symbol
    ( t_Symbol             const& sym,
      std::vector<double>  const& state_vals,
      std::vector<double>  const& cst_vals,
      std::vector<double>  const& dom_vals,
      std::vector<double>  const& input_vals = std::vector<double>(),
      arma::mat*                  A0 = nullptr )   //!< on request: the zeroth-order block (sym.vCoeff0), evaluated
    const;


protected:

  //! @brief Auxiliary-state definition generated by reduce_order().
  struct t_AuxDef
  {
    FFVar                            aux;          //!< Auxiliary state variable
    FFVar                            expr;         //!< Collocated expression defining aux
    std::set<FFVar,lt_FFVar>         dom;          //!< Auxiliary domain
    std::set<FFVar,lt_FFVar>         diff_dom;     //!< Domain directions differentiated in the defining link
    FFVar                            parent;       //!< Immediate parent state (primitive or auxiliary from a prior pass)
    int                              block_id = 0;
  };

  //! @brief Auxiliary definitions generated by reduce_order(), in creation order.
  std::vector<t_AuxDef>              _auxDef;

  //! @brief True when a collocated state is an auxiliary generated by reduce_order().
  bool _is_auxiliary_state
    ( FFVar const& state )
    const;

  //! @brief Immediate parent state of an auxiliary (primitive or auxiliary from a prior pass).
  //! Returns nullptr if the state is not an auxiliary or has no recorded parent.
  FFVar const* _auxiliary_parent
    ( FFVar const& state )
    const;

  //! @brief Ultimate primitive root of an auxiliary state.
  //! Walks the parent chain through multi-pass auxiliaries (e.g. Dz_Drd_Cd → Drd_Cd → Cd)
  //! until a non-auxiliary state is reached.  Returns nullptr if the state is not an auxiliary.
  //! Returns the state itself if it is not auxiliary (i.e. already primitive).
  FFVar const* _auxiliary_primitive_root
    ( FFVar const& state )
    const;


public:

  //! @brief Reference-free structural DAE decomposition of one equation block,
  //! produced by _structural_dae_decomposition().  Splits the block's states by
  //! time-differentiation and identifies the algebraic constraints.  All reads
  //! are structural (DAG FFPartial incidence) -- no reference evaluation.
  struct t_StructuralDecomp
  {
    //! States carrying a time-derivative d_t (dynamic / differential states).
    std::vector<FFVar> dyn;
    //! States appearing in the block but with NO d_t (algebraic / auxiliary).
    std::vector<FFVar> alg;
    //! Equations with NO d_t of any state (algebraic constraints).
    std::vector<FFVar> constraints;
    //! Per-constraint flag: 1 if it carries a spatial derivative (order-raising,
    //! PDE-auxiliary character), 0 if purely algebraic (DAE-index character).
    std::vector<char>  constraint_has_spatial;
    //! Per-constraint flag: 1 if EqnRole::LINK (an order-reduction auxiliary
    //! definition, excluded from the differential-index matching).
    std::vector<char>  constraint_is_link;
  };

  //! @brief Result of _structural_index_analysis(): the differential index of a
  //! block from the constraint <-> algebraic-variable matching (structural,
  //! reference-free).  index = 0 (no algebraic variables, or all of them are
  //! order-reduction auxiliaries), 1 (every algebraic variable is pinned by one
  //! algebraic constraint), 2/3/... (exact high index resolved by Stage-2
  //! structural-differentiation reachability), or -1 (a high-index witness no
  //! differentiation can expose -- structurally singular).
  struct t_IndexResult
  {
    int                index = 0;
    //! Index-relevant algebraic variables (order-reduction auxiliaries removed).
    std::vector<FFVar> index_alg;
    //! High-index witnesses: algebraic variables left unpinned by the matching.
    std::vector<FFVar> unmatched;
    //! A pinning constraint carries a spatial derivative -> parabolic character.
    bool               parabolic_character = false;
  };

  //! @brief Result returned by classify()
  struct t_Classify
  {
    //! PDE type classification.
    EqnType                 type = UNDETERMINED;
    //! Numerically evaluated coefficient matrices A_i at the query point
    std::vector<arma::mat>  Ai;
    //! Characteristic speeds/eigenvalues for the last sampled direction
    arma::cx_vec            eigenvalues;
    //! Per-direction eigenvalue data: {transverse/principal direction \xi, eigenvalues}
    std::vector< std::pair< std::vector<double>, arma::cx_vec > > eigendata;
    //! True when the chosen evolution coefficient matrix is numerically singular
    bool                    At_singular = false;
    //! True when the system has a singular evolution matrix and should be
    //! regarded as a descriptor first-order system unless further structure is found
    bool                    descriptor = false;
    //! True when descriptor structure plus LINK/interior metadata indicate a parabolic closure
    bool                    parabolic_structure_detected = false;
    //! True when an evolution direction is known/detected and all sampled speeds are real
    bool                    evolution_hyperbolic = false;
    //! True when real characteristic directions are detected but no evolution direction is available
    bool                    spatially_characteristic = false;
    //! True when real speeds are found but the eigenbasis is numerically defective/ill-conditioned
    bool                    weak_hyperbolic = false;
    //! True when the principal symbol is singular in all sampled directions
    bool                    degenerate = false;
    //! @brief ELLIPTIC established on the Douglis-Nirenberg WEIGHTED symbol of an order-reduced system (2026-10-04):
    //! the LINK rows' zeroth-order entries are principal there, and the first-derivative symbol alone is singular
    //! for every direction (Laplace reduced to (u, Du): the two LINK rows are parallel).  Report-only: the
    //! redundant-claim gate keeps treating such a block as it did a DEGENERATE one.
    bool                    dn_elliptic = false;
    //! Index of the chosen evolution domain in the local symbol ordering, or npos
    size_t                  evolution_dom_idx = std::numeric_limits<size_t>::max();
    //! Whether the evolution direction was supplied explicitly by the caller
    //! @brief the differential core's principal symbol is RECTANGULAR
    //! (retained differential rows differentiate more states than there are rows).
    //! A SHAPE fact about the symbol, ORTHOGONAL to the type: a block can be rectangular AND
    //! parabolic in character (measured: 22 of 22 rectangular probe lines corpus-wide report
    //! parab_char=y).  Recorded so the type axis can eventually carry character alone; the
    //! DIFFERENTIAL_RECTANGULAR enum value stays for now and its three consumers still read
    //! it (they migrate in STAGE B).
    bool                    symbol_rectangular = false;

    //! @brief TRUE when the principal-symbol analysis actually RAN for this
    //! block (the direction sweep / eigen analysis), FALSE on every early return.
    //! WHY IT EXISTS.  Retiring DIFFERENTIAL_RECTANGULAR (D5) sends the At-non-singular
    //! rectangular class to UNDETERMINED, where it joins genuinely inconclusive blocks.  Those
    //! are NOT the same statement -- the distinction is deliberate -- and
    //! symbol_rectangular records the shape, not whether anything was measured.  With this
    //! field the pair (symbol_analysed, symbol_rectangular) says exactly which case a block is:
    //!   analysed=1            -> the verdict is a measurement
    //!   analysed=0 rect=1     -> no square test applies (the old DIFFERENTIAL_RECTANGULAR)
    //!   analysed=0 rect=0     -> an earlier return: algebraic-lumped/field
    bool                    symbol_analysed = false;
    bool                    evolution_user_supplied = false;
    //! Whether a steady/marching evolution direction was auto-detected
    bool                    evolution_auto_detected = false;
    //! Numerical rank and conditioning diagnostics for the evolution matrix
    arma::uword             rank_evolution = 0;
    double                  sigma_min_evolution = 0.;
    double                  sigma_max_evolution = 0.;
    double                  cond_evolution = std::numeric_limits<double>::infinity();
    //! Largest imaginary part and eigenvector condition observed in the sweep
    double                  max_imag_eig = 0.;
    double                  max_eigvec_cond = 0.;
    //! Imaginary-part/rank tolerance used during classification
    double                  imag_tol = 1e-8;
    //! ---- (character, index) representation: structural differential index ----
    //! Differential (DAE) index from the reference-free structural analysis
    //! (_structural_index_analysis): 0 = regular (no algebraic states beyond
    //! order-reduction auxiliaries), 1 = index-1 (every algebraic state pinned
    //! by a single elimination), 2/3/... = exact high index from Stage-2
    //! structural differentiation, -1 = unresolved (structurally singular
    //! witness).  The `type` enum above remains the derived
    //! character view; this is the orthogonal index axis.
    int                     differential_index = 0;
    //! True when the index-1 elimination is closed by a SPATIAL (order-raising)
    //! constraint, i.e. the descriptor closes parabolically.  Reference-free.
    bool                    structural_parabolic = false;
    //! True when the block is a PURE differential-algebraic system: algebraic
    //! states pinned by ALGEBRAIC (non-spatial) constraints with NO transverse
    //! spatial domain.  The principal symbol is then rectangular (algebraic
    //! states are not derivative columns) so classify() returns UNDETERMINED;
    //! this structural verdict names it DAE.  Reference-free.
    bool                    structural_dae = false;
  };

  //! @brief Precomputed characteristic decomposition for one face direction.
  //! Built by _build_face_data() from the cached classification result.
  struct t_FaceData
  {
    //! Domain variable for this face direction
    FFVar const*  pdom_var;
    //! Index of this domain variable in _symbol.vDom / _classification.Ai
    size_t        dom_idx;
    //! Positive-eigenvalue part of A_i (right-going characteristics)
    arma::mat     Aplus;
    //! Negative-eigenvalue part of A_i (left-going characteristics)
    arma::mat     Aminus;
  };


protected:

  //! @brief Return the cached face decompositions associated with a block id.
  std::vector<t_FaceData> const& _face_data_for_block
    ( int block_id ) const;

  //! @brief Find the characteristic face decomposition for a domain variable.
  static t_FaceData const* _face_data_for_dom
    ( std::vector<t_FaceData> const& fdata, FFVar const* pvar );

  //! @brief Cached PDE classification result (populated by classify_pde).  Held on the HEAP, not by value:
  //! t_Classify carries arma::cx_vec, whose 16-byte alignment would otherwise make alignof(FFModel) 16, and
  //! GCC mis-lays an over-aligned VIRTUAL base (every solver virtually inherits FFModel) -- ODESLV_BASE's
  //! constructor then does a 16-byte-aligned store to an 8-aligned address and faults at -O2.  Measured
  //! 2026-09-26 on GCC 13.3; clang is unaffected.  The reference keeps every use site as it was.
  std::unique_ptr<t_Classify>               _pClassification;
  t_Classify&                               _classification;

  //! @brief Monotone generation for successful PDE classifications.
  //! Used to invalidate caches whose values depend on the principal symbol,
  //! even when _classified remains true across an in-place reclassification.
  size_t                                    _classificationSerial = 0;

  //! @brief Evolution-domain variable used by classify_pde() for IC_AUTO dispatch.
  FFVar                                     _evolution_dom_var;
  bool                                      _evolution_dom_user = false;
  bool                                      _evolution_dom_set  = false;

  //! @brief (row_id, dom_id) -> does the residual apply a Partial in that
  //! direction.  Model-structural; built once per setup by _eqn_differentiates_in(),
  //! invalidated with the other structural caches.  Consumed through PlanInput.
  mutable std::map<std::pair<size_t,size_t>,bool> _eqnDiffIn;
  mutable bool                              _eqnDiffInReady = false;

  std::map<std::pair<size_t,size_t>,bool> const& _eqn_differentiates_in() const;

  int         _evoInferProvenance  = 0;    //!< 0=none 1=INITIAL-pinned 2=name-based (weak)

  //! @brief Optional extra STATE references sampled by the reference-robustness
  //! probe to check that the principal-symbol type and characteristic split are
  //! invariant across the operating envelope.  Each entry is a partial map
  //! {state FFVar -> value}; unlisted states keep their base reference value.
  //! Populated via add_classification_reference_sample(); consumed only when
  //! MC__OCFESLV_REFERENCE_ROBUSTNESS_PROBE is defined (print-only).
  std::vector<std::map<FFVar,double,lt_FFVar>>  _classRefSamples;

  //! @brief Cached verdict of the most recent reference-robustness probe (true if
  //! no sampled state changed the PDE type or the per-face characteristic split).
  //! Defaults true; remains true when the probe is not compiled in.
  mutable bool                              _refRobustOk = true;

  //! @brief Per-face-direction characteristic decompositions (one entry per domain variable)
  std::vector<t_FaceData>                   _face_data;

  std::map<int,t_Classify>                  _blockClassification;
  std::map<int,std::vector<t_FaceData>>     _blockFaceData;

  //! @brief Per-block structural-index-analysis cache (memoisation): the classify
  //! pass, the structural probe, and the reduction-plan builder share ONE analysis
  //! per block.  Cleared at setup() start and after _reduce_high_index mutates _mEqn.
  mutable std::map<int,t_IndexResult>       _idxCache;

  //! @brief Build characteristic decomposition for each face direction.
  //! Called by classify_pde(); results are stored in _face_data.
  void _build_face_data
    ();

  //! @brief Build characteristic decomposition for an arbitrary block symbol/classification.
  std::vector<t_FaceData> _build_face_data
    ( t_Symbol const& sym, t_Classify const& cls )
    const;


  //! @brief Infer the evolution/marching domain for automatic classification.
  FFVar const* _infer_evolution_domain
    ();

  //! @brief Re-evaluate the principal symbol at several STATE references and check
  //! that the PDE type and per-face characteristic split are invariant.  Setup-time
  //! structural analysis linearizes ONCE at _classification_reference; for a
  //! nonlinear A_z(u) that point can get the structure wrong and the solver's
  //! per-iteration re-linearization cannot rescue it (the structure is baked into
  //! the assembled system).  Print-only; returns the all-invariant verdict.  Opt-in
  //! via MC__OCFESLV_REFERENCE_ROBUSTNESS_PROBE.
  bool _probe_reference_robustness
    ( std::vector<double> const& state_ref,
      std::vector<double> const& input_ref,
      std::vector<double> const& cst_ref,
      std::vector<double> const& dom_ref,
      FFVar const*               time_dom,
      unsigned const             n_sample,
      double   const             imag_tol ) const;


public:

  //! @brief Return cached PDE classification (valid after setup())
  t_Classify const& pde_type
    ()
    const
    { return _classification; }

  //! @brief The structural (differential) index of block @p block_id WITH RESPECT TO @p dir -- any declared domain,
  //! not only the evolution direction.  Index 0 means no algebraic state needs exposing; 1 that every algebraic state
  //! is pinned directly; k > 1 that some witness is reached only after k-1 differentiations of a constraint.
  //! ANALYSIS ONLY: FFModel reduces a high index in the EVOLUTION direction alone, so a high index reported here for
  //! a spatial direction is information for the consumer, not something setup() acts on.  What it means depends on
  //! that consumer: a collocation solver differentiates the interpolant and is unaffected (at a cost in conditioning),
  //! while an integrator marching in that direction, or a shooting or BVP solver, faces the usual high-index system.
  //! Uncached, unlike the evolution-direction analysis.
  t_IndexResult structural_index
    ( int const block_id, FFVar const& dir )
    const
    { return _structural_index_analysis_uncached( block_id, &dir ); }

  //! @brief What one FACE of one spatial direction of one block asks for, and what the model puts there.
  //! @a outgoing is the number of characteristic directions LEAVING the domain at this face -- the number of
  //! conditions it requires; @a covered how many of those the block's own rows already supply; @a appended how many
  //! the hyperbolic closure had to add; and @a rows_at_face how many rows the model pins at this face in total.
  //! @a rows_at_face exceeding @a incoming means data the face cannot take: at a face, the characteristics LEAVING
  //! are closed by rows derived from the PDE (what @a appended supplies) and only those ENTERING may be given data,
  //! so an excess over @a incoming is the classical ill-posedness of conditions at the wrong end -- invisible to any
  //! count, since the system stays square.  OCFESLV enforces this for its own solves (its hyperbolic BC guard, which
  //! refuses the model); the account is here so a consumer that is NOT OCFESLV can make the same check.
  struct t_FaceConditions
  {
    int         block_id     = 0;
    std::string direction;
    int         face         = 0;      //!< FFDom::LB or FFDom::UB
    size_t      outgoing     = 0;
    size_t      incoming     = 0;      //!< characteristics ENTERING here: the conditions this face needs
    size_t      covered      = 0;
    size_t      appended     = 0;
    size_t      rows_at_face = 0;
  };

  //! @brief The per-face account of the last setup() (see t_FaceConditions), one entry per block, spatial direction
  //! and face for which a characteristic split was computed.
  std::vector<t_FaceConditions> const& face_conditions
    ()
    const
    { return _faceConditions; }

  //! @brief One thing that makes this model ill-posed, or not obviously well-posed, for SOME consumer.  The model
  //! states it and does not act on it: a collocation solver may handle what an integrator cannot, so the reading
  //! belongs to whoever solves.  @a detail is a sentence a modeller can act on.
  struct t_Finding
  {
    enum Kind { COMPLEX_CHARACTERISTICS,   //!< speeds are not real: no continuous dependence on the data
                SINGULAR_INDEX,            //!< a witness no differentiation exposes (structurally singular)
                DOF_IMBALANCE,             //!< rows - unknowns is not zero
                FACE_EXCESS,               //!< more condition rows at a face than characteristics entering it
                REDUNDANT_INITIAL_DATA };  //!< more initial rows than the reduced system lets one choose
    Kind        kind = COMPLEX_CHARACTERISTICS;
    int         block_id = 0;
    std::string detail;
  };

  //! @brief Everything the model's own analyses found against it, collected by the last setup().  EMPTY is the
  //! good case.  Nothing here is fatal: which findings matter depends on the consumer -- collocation differentiates
  //! an interpolant and survives a spatial chain that would stop an integrator, while complex characteristics stop
  //! everyone.  Read it with report(), which prints the same list.
  std::vector<t_Finding> const& wellposedness
    ()
    const
    { return _findings; }

  //! @brief This model read as an ODE/DAE in the evolution direction -- what a time integrator would need to know
  //! before it could take it.  @a lumped is the precondition (every state on the evolution direction alone; a
  //! distributed model needs its spatial directions discretised first).  @a blocker is EMPTY when an integrator could
  //! take the model, and otherwise says in one sentence why not.
  struct t_DynamicForm
  {
    bool        lumped         = false;
    bool        decoupled      = true; //!< every row differentiates ONE state and every state's derivative appears
                                       //!< in ONE row, so the system is x' = f after dividing by each row's own
                                       //!< coefficient.  False means the derivatives are coupled through a MASS
                                       //!< MATRIX: still an ODE (index 0 when that matrix is nonsingular, which is
                                       //!< what the classification reports), but an explicit integrator needs it
                                       //!< inverted, while the (coefficient, remainder) pair suits IDAS directly.
    size_t      differential   = 0;   //!< states carrying a derivative in the evolution direction
    size_t      algebraic      = 0;   //!< the rest: an integrator needs a DAE solver for these
    int         declared_index = 0;   //!< the index the reduction found (0 when none ran)
    bool        reduced        = false;
    bool        resolved       = true;//!< the reduction completed, so the system handed over is index 0 or 1
    bool        initial_value  = true;//!< nothing is imposed at the FAR END of the evolution domain
    size_t      free_initial   = 0;   //!< initial values that may be chosen (see t_InitialData)
    std::string blocker;
  };

  //! @brief This model as an ODE/DAE (see t_DynamicForm), as read by the last setup().
  t_DynamicForm const& dynamic_form
    ()
    const
    { return _dynamicForm; }

  //! @brief Per-block PDE classification (block_id -> t_Classify), populated by
  //! setup()/classification when CLASSIFY.MODE is enabled.  Empty when only the
  //! aggregate (whole-system) classification was produced.  Diagnostic use.
  std::map<int,t_Classify> const& block_classification
    ()
    const
    { return _blockClassification; }


protected:

  //! @brief Canonical rule for "differential state": if op is OpP(state, evolution_var) for a state
  //! in varmap, return that state; else a null FFVar.  Single definition shared by the DAE
  //! decomposition (dyn set) and marching-transfer auto-detection.
  FFVar _time_derivative_state
    ( FFOp const* op, t_Var const& varmap )
    const;

  //! @brief Set of states carrying a d/d(evolution-var) partial anywhere in the given equations,
  //! via _time_derivative_state.  This IS the differential-state ("dyn") identification, evaluated
  //! wherever needed (globally for marching transfers, per-block in the decomposition).
  std::set<FFVar,lt_FFVar> _states_with_evolution_derivative
    ( t_Eqns const& eqns, t_Var const& varmap )
    const;

  //! @brief Structural (reference-free) DAE decomposition of a block: splits
  //! states into time-differentiated (dyn) vs algebraic/auxiliary (alg) and
  //! identifies algebraic constraints (no d_t), flagging those carrying a
  //! spatial derivative (order-raising) vs purely algebraic (DAE index).  Reads
  //! only DAG FFPartial incidence; no reference evaluation.  (Stage 1a)
  t_StructuralDecomp _structural_dae_decomposition
    ( int const block_id )
    const
    { return _structural_dae_decomposition( block_id, _evolution_dom_set? &_evolution_dom_var: nullptr ); }

  //! @brief As above, but splitting by derivatives in @p dir rather than in the evolution direction: a state is
  //! DIFFERENTIAL in @p dir when it appears under a derivative with respect to it.  This is what makes an index
  //! reading in a non-evolution direction meaningful -- with the evolution split, every state of a spatial chain
  //! looks algebraic and pinned, and the chain reads index 1.
  t_StructuralDecomp _structural_dae_decomposition
    ( int const block_id, FFVar const* dir )
    const;

  //! @brief Differential index of a block via constraint <-> algebraic-variable
  //! matching (structural, reference-free).  index-1 <=> every algebraic
  //! variable is pinned by an algebraic constraint (one elimination closes it);
  //! high-index <=> some algebraic variable is unpinned (Stage-2 differentiation
  //! needed).  Order-reduction (LINK) auxiliaries are excluded.  (Stage 1a)
  t_IndexResult _structural_index_analysis
    ( int const block_id )
    const;

  //! @brief Uncached body of _structural_index_analysis (see _idxCache), with respect to @p dir: the direction
  //! whose derivatives make a state differential, so that one differentiation of a constraint containing that state
  //! exposes its defining row.  Pass the evolution direction for the classical index; pass any other declared domain
  //! for that direction's index (the paper's nu_x).  A null @p dir leaves every state algebraic.
  t_IndexResult _structural_index_analysis_uncached
    ( int const block_id, FFVar const* dir )
    const;

  //! @brief Classify the PDE system from numerically evaluated coefficient matrices.
  //! @param Ai        One coefficient matrix per domain variable (from _eval_symbol())
  //! @param vDom      Domain variable ordering matching Ai (from t_Symbol::vDom)
  //! @param time_dom  Pointer to the time-like domain variable, or nullptr
  //! @param n_sample  Number of spatial directions sampled on the unit sphere
  //! @param imag_tol  Threshold for declaring an eigenvalue real
  //! @param has_algebraic_rows  TRUE when the block carries equations EXCLUDED from the
  //!        differential symbol (t_Symbol::vAlgEqn non-empty).  Without it the
  //!        classifier CANNOT SEE A DAE.  t_Symbol's own documentation says algebraic
  //!        equations are "excluded from the differential symbol (they would be zero
  //!        rows)", so for a genuine DAE the zero row never reaches A_e, ri.singular is
  //!        false, and the verdict came out DIFFERENTIAL_ORDINARY -- measured on
  //!        OCFE_DAE0's index-1 pendulum: 5 states, 5 interior equations, yet sym=4x4
  //!        with Ae_singular=n and sigma_min=sigma_max=1.  Defaulted false so every
  //!        existing caller keeps its behaviour unchanged.
  t_Classify classify
    ( std::vector<arma::mat> const& Ai,
      std::vector<FFVar>     const& vDom,
      FFVar const*                  time_dom = nullptr,
      unsigned const                n_sample  = 32,
      double   const                imag_tol  = 1e-8,
      bool     const                has_algebraic_rows = false,
      std::vector<arma::mat> const* Ai_dn = nullptr,   //!< weighted first-order symbol (see t_Classify::dn_elliptic)
      arma::mat const*              E0_dn = nullptr )  //!< its zeroth-order part
    const;


protected:

  //! @brief Compute and cache the PDE classification with input references.
  //! @param input_ref  Reference values for input FFVars (length = _vInp.size()).
  //!                   Distributed inputs are evaluated as spatially constant
  //!                   at these reference values during classification.
  bool _classify_pde
    ( std::vector<double> const& state_ref,
      std::vector<double> const& input_ref,
      std::vector<double> const& cst_ref,
      std::vector<double> const& dom_ref,
      FFVar const*               time_dom  = nullptr,
      unsigned const             n_sample  = 32,
      double   const             imag_tol  = 1e-8 );


private:

  typedef SMon<FFVar,lt_FFVar> t_SMon;


public:

  //! @brief Reason for the most recent setup() outcome.  setup() still returns a
  //! bool (false on any failure); setup_status() exposes WHY, since a bare false
  //! cannot distinguish a malformed model from a mis-placed boundary condition.
  //! Defined here (ahead of the _setupStatus member) so the member's type is
  //! complete at its point of declaration.
  enum class SetupStatus {
    OK = 0,                  //!< setup() completed successfully
    INCONSISTENT_MODEL,      //!< _check_consistency failed (malformed model)
    NODE_SETUP_FAILED,       //!< collocation-node generation failed for a domain
    CLASSIFICATION_FAILED,   //!< _classify_pde failed under CLASS_STRICT
    HYP_INCOMING_BC,         //!< hyperbolic incoming-BC count mismatch at a face
    HYP_BC_MISDIRECTED,      //!< hyperbolic incoming-BC pins no incoming characteristic
                             //!< in value OR derivative (criterion #3 direction reject)
    // HYP_CLOSURE_REPREP removed at rev269: declared but never set -- the auto-closure re-prep path reports
    // through HYP_INCOMING_BC / HYP_BC_MISDIRECTED.  The other 14 values are all reachable (checked by
    // counting SetupStatus:: references across both headers).
    INTERFACE_PLAN_INVALID,  //!< interface-plan draft failed to build/validate
    AUDIT_NONSQUARE,         //!< structural FFDep audit found a non-square system
    DERIV_CACHE_FAILED,      //!< eager derivative-cache setup failed
    LINEAR_CACHE_FAILED,     //!< eager linear-eval-cache setup failed
    REDUCED_DOF_INCONSISTENT,//!< post-index-reduction DOF audit found a structurally rank-deficient
                             //!< (or non-square) reduced system
    EVOLUTION_DOMAIN_UNSUITABLE, //!< _validate_evolution_domain(): unsuitable evolution domain or an IC
                             //!< that cannot be set consistently (fatal only under FATAL.EVOLUTION_DOMAIN)
    DETERMINACY_UNDETERMINED,//!< numerical determinacy audit found state content in the null
                             //!< space of the assembled Jacobian (fatal only under FATAL.DETERMINACY)
    MODEL_COPY_FAILED,       //!< _copy_usr_to_local failed (could not build working DAG)
    CAPTURE_NESTED_REFUSED,  //!< an evolution-direction reduction is nested inside another evolution
                             //!< reduction: the inner value is a post-solve captured input, so the
                             //!< outer's in-solve materialised operand would read a stale (window-
                             //!< lagged) value.  Refused rather than returned silently wrong.
      HYP_CLOSURE_MISSING,     //!< AUTO.HYP_CLOSURE off and the model leaves outgoing-characteristic rows missing at
                             //!< an outflow face of a hyperbolic block (2026-10-07; appended last: values unchanged)
      BACKEND_UNAVAILABLE      //!< an option explicitly requests a linear-algebra backend this build does not include
                             //!< (SOLVE_SPQR, DET_SPQR without SuiteSparseQR; DET_EIGEN without Eigen).  Default and
                             //!< AUTO choices fall back instead (2026-10-09; appended last: values unchanged)
  };

  //! @brief Process-wide lock for every operation that reads or modifies a USER DAG: setup() and fdiff() here, and
  //! the ODE solvers' setup().  Solver copies (setup( IVP ), and the COPY policy of the DAG operations) share their
  //! source's user DAG, and MC++ DAG traversals (subgraph) are not safe on a shared graph -- concurrent setups
  //! raced there (ODESLV_LVthreads Part B, 2026-09-28: failed copies and segfaults on a multi-core machine).
  //! Recursive, because setups nest (a SYMDIFF product is set up inside deriv).  Solves are not locked: each
  //! solver evaluates on its own working DAG.
  static std::recursive_mutex& dag_mutex
    ()
    { static std::recursive_mutex m;  return m; }

  //! @brief The status of the last setup() -- shared by OCFESLV and ODESLV.
  SetupStatus setup_status() const { return _setupStatus; }

  //! @brief Human-readable description of a SetupStatus value.
  static char const* setup_status_str( SetupStatus s )
  {
    switch( s ){
    case SetupStatus::OK:                     return "ok";
    case SetupStatus::INCONSISTENT_MODEL:     return "inconsistent model (consistency check failed)";
    case SetupStatus::NODE_SETUP_FAILED:      return "collocation-node setup failed";
    case SetupStatus::CLASSIFICATION_FAILED:  return "PDE classification failed (CLASS_STRICT)";
    case SetupStatus::HYP_INCOMING_BC:        return "hyperbolic incoming-BC count mismatch (wrong number of BCs at a face)";
    case SetupStatus::HYP_BC_MISDIRECTED:     return "hyperbolic incoming-BC mis-directed (pins no incoming characteristic in value or derivative; prescribes an outgoing characteristic)";
    case SetupStatus::INTERFACE_PLAN_INVALID: return "interface plan build/validate failed";
    case SetupStatus::AUDIT_NONSQUARE:        return "structural audit: collocation system not square";
    case SetupStatus::DERIV_CACHE_FAILED:     return "derivative-cache setup failed";
    case SetupStatus::LINEAR_CACHE_FAILED:    return "linear-evaluation-cache setup failed";
    case SetupStatus::REDUCED_DOF_INCONSISTENT: return "post-index-reduction DOF audit: IC/BC count "
        "inconsistent with the reduced degrees of freedom (structurally rank-deficient or non-square)";
    case SetupStatus::EVOLUTION_DOMAIN_UNSUITABLE: return "evolution-domain validation: unsuitable "
        "evolution domain or an initial condition that cannot be set consistently";
    case SetupStatus::DETERMINACY_UNDETERMINED: return "determinacy audit: the assembled Jacobian has "
        "a null space with state content (primal undetermined)";
    case SetupStatus::MODEL_COPY_FAILED:      return "could not copy user model into working DAG";
    case SetupStatus::HYP_CLOSURE_MISSING:    return "hyperbolic block: AUTO.HYP_CLOSURE is off and outgoing-characteristic "
        "rows are missing at an outflow face (close it in the model, or turn the option on)";
    case SetupStatus::BACKEND_UNAVAILABLE:    return "an option requests a linear-algebra backend this build does not "
        "include (see build_info()['backends']; the default or AUTO choice falls back instead)";
    case SetupStatus::CAPTURE_NESTED_REFUSED: return "nested evolution-direction reduction: an inner "
      "reduction's captured (post-solve) value feeds an outer reduction's in-solve operand, which "
      "would read a window-lagged value -- refused (consume the inner reduction in an output/objective, "
      "not inside another evolution reduction)";
    }
    return "unknown";
  }



protected:

  //! @brief Reason for the most recent setup() outcome; see setup_status().
  SetupStatus                        _setupStatus = SetupStatus::OK;

  //! @brief State-ids of GENUINE (non-auxiliary) algebraic states that are
  //! "value-slaved": pinned by a derivative-free algebraic constraint
  //! (e.g. P - Rg c = 0).  Such a state is determined POINTWISE at every node
  //! including interfaces, so its value-continuity is implied by the continuity
  //! of the states it depends on and an interface C0-continuity claim on it is
  //! redundant (the interface analogue of _generate_algebraic_boundary_closure's
  //! "closed by its constraint, never by continuity" rule).
  std::set<size_t> _value_slaved_algebraic_state_ids
    ()
    const;

  //! @brief Order>2 boundary closure: retire LINK_j traces at boundary faces
  //! over-determined by a co-located value+derivative condition on one
  //! primitive's chain.  Inert (byte-identical) unless such a co-occurrence
  //! exists.  See plan §13.
  void _displace_overdetermined_boundary_links
    ( std::vector<t_Eqn>& wEqn )
    const;

  //! @brief Reduce high-order PDEs to first-order in the local working model.
  //! Requires _issetup == true and is only called after _copy_usr_to_local().
  bool _reduce_order
    ( bool const reuse_aux = true );


public:

  //! @brief Persisted high-index reduction plan (Pantelides + dummy-derivative
  //! matching), decoupled from t_IndexResult and from any classify pass so a
  //! marching/MOL mode can query the reduction independently.  Filled by
  //! _build_reduction_plan() in setup(); consumed by _reduce_high_index().
  struct t_ReductionPlan
  {
    //! One differentiated-constraint assignment: differentiate `constraint`
    //! `n_diff` times to pin `pinned_var` (the dummy-derivative match).
    struct t_Assign
    {
      int   block_id = 0;
      FFVar constraint;      //!< original interior algebraic constraint residual
      int   n_diff   = 0;    //!< d/dt rounds applied (= local index - 1)
      FFVar pinned_var;      //!< the witness pinned by the differentiated form
    };
    std::vector<t_Assign> assigns;
    int  max_index = 0;      //!< highest block differential index encountered
    bool resolved  = true;   //!< every high-index witness was structurally exposable
    bool empty() const { return assigns.empty(); }
    void clear() { assigns.clear(); max_index = 0; resolved = true; }
  };


protected:

  //! @brief Add (or redeclare) a state of the WORKING model on @p dom -- an auxiliary minted during setup().  The
  //! classification is invalidated and the derivative caches are dropped, as for any model change.
  void _add_working_state
    ( FFVar const& var, std::vector<FFVar> const& dom )
    {
      auto [it,ins] = _mVar.insert( {var,{}} );
      if( !ins ) it->second.clear();
      for( auto const& d : dom ) it->second.insert( d );
      _classified = false;
      _on_model_changed( ModelChange::DERIVATIVES );
    }

  //! @brief Add (or redeclare) an input of the WORKING model on @p dom -- e.g. the input that holds a deferred value.
  void _add_working_input
    ( FFVar const& var, std::vector<FFVar> const& dom )
    {
      auto [it,ins] = _mInp.insert( {var,{}} );
      if( !ins ) it->second.clear();
      for( auto const& d : dom ) it->second.insert( d );
      _classified = false;
      _on_model_changed( ModelChange::DERIVATIVES );
    }

  //! @brief Record that @p aux is defined by @p expr on @p dom (differentiated in @p diff_dom), in block @p block_id.
  void _register_aux_def
    ( FFVar const& aux, FFVar const& expr, std::vector<FFVar> const& dom,
      std::set<FFVar,lt_FFVar> const& diff_dom, int block_id, FFVar const& parent = FFVar() )
    {
    t_AuxDef def;
    def.aux         = aux;
    def.expr        = expr;
    def.dom.insert( dom.begin(), dom.end() );
    def.diff_dom    = diff_dom;
    def.parent      = parent;
    def.block_id    = block_id;
    _auxVarID.insert( aux.id() );
    _auxDef.push_back( std::move( def ) );
  }


  //! @brief This model read as an ODE/DAE (see t_DynamicForm; public through dynamic_form()).
  t_DynamicForm                              _dynamicForm;

  //! @brief Fill _dynamicForm from the analyses already run; no new mathematics.
  void _build_dynamic_form
    ();

  //! @brief What the model's analyses found against it (see t_Finding; public through wellposedness()).
  std::vector<t_Finding>                     _findings;

  //! @brief Collect the findings of the other analyses into _findings; no new mathematics.
  void _collect_findings
    ();

  //! @brief The per-face condition account of the last setup() (see t_FaceConditions; public through face_conditions()).
  std::vector<t_FaceConditions>              _faceConditions;

  //! @brief The account of the initial data of the last setup() (see t_InitialData; public through initial_data()).
  t_InitialData                              _initialData;

  //! @brief The degree-of-freedom balance of the last setup() (see t_DofBalance; public through dof_balance()).
  t_DofBalance                               _dofBalance;

  //! @brief Count the working model's rows and unknowns symbolically and fill _dofBalance; report a model that is
  //! not balanced, and say nothing when it is.
  void _audit_dof_balance
    ();

  //! @brief The deferred values of this model (see t_DeferredValue; public through var_deferred()).
  std::vector<t_DeferredValue>               _deferredValue;

  //! @brief Is this reduction a DEFERRED VALUE?  A reduction consuming exactly the evolution direction -- a definite
  //! integral over it, or an evaluation at one point of it -- is CAPTURED: the node becomes an input written after the
  //! solve (LATCH: R(tau); ACCUM: the integral, summed over windows), feeding outputs and objectives only.  Every
  //! other reduction (spatial, or consuming several directions) keeps its in-solve quadrature or point form.
  //! This is the DECISION; _reduce_order() applies it, at three sites (an equation, an output integral, an output
  //! point), because the application needs that pass's local state.
  bool _is_deferred_reduction
    ( std::set<FFVar,lt_FFVar> const& consumed )
    const
    { return _evolution_dom_set && _evolution_dom_var.dag()
          && consumed.size() == 1 && consumed.count( _evolution_dom_var ); }

  //! @brief defining constraints kept (not value-slaved) because a bare term is an alias.
  mutable size_t _reuseAwareKept = 0;

  //! @brief Resolve WHICH INITIAL-role equation CLOSES each candidate state in @p cand, as opposed
  //! to merely REFERENCING it, returning an INJECTIVE map state -> index into @p eqns.
  //!
  //! An IC residual may legitimately carry other states as ARGUMENTS: a moving-mesh driver's
  //!     U_BC := u - 0.5*( 1 - tanh( ( x - s(t) ) / delta ) )
  //! CLOSES u but merely EVALUATES the (also differential) mesh state x.  The historical test
  //! ("the subgraph references st") attributed that single equation to BOTH u and x, so both
  //! claimed the SAME residual rows; the k>0 marching override then wrote two different value-pins
  //! into one row set, silently dropping one state's continuity and leaving the other UNCLOSED at
  //! the window LB -> singular/inconsistent augmented interface system at every window after the
  //! first.  Attribution must therefore be exclusive.
  //!
  //! Rule (a strict NARROWING of the old behaviour -- identical whenever an INITIAL residual
  //! references exactly one candidate, which is the entire existing corpus):
  //!   (a) refs(E) = candidates referenced by E.  |refs| == 1  ->  E closes that state.
  //!   (b) |refs| > 1 -> the closure is the candidate E is AFFINE in (FFDep type L).  A Dirichlet-
  //!       form IC enters its own state linearly, while an argument state enters through a
  //!       nonlinear wrapper (tanh/exp/product).  A unique affine candidate resolves E.
  //!   (c) still ambiguous (or FFDep unusable because E carries an external OpP/OpI/OpEval op) ->
  //!       first unclaimed candidate in refs order, reported at DISPLAY_LEVEL>=1.
  //! Equations are visited in insertion order, so a user INITIAL equation always takes precedence
  //! over a synthesized (consistent-IC/re-pivot) one.  A candidate absent from the returned map has
  //! NO INITIAL closure of its own -- its evolution-LB value is fixed by the algebraic closure
  //! (the index-2 equidistribution case) and it must not be given an IC row or a transfer.
  std::map<FFVar,size_t,lt_FFVar> _match_initial_closures
    ( t_Eqns const& eqns, t_Var const& varmap, std::set<FFVar,lt_FFVar> const& cand,
      char const* who ) const;

  //! @brief Persisted high-index reduction plan and post-reduction DOF audit
  //! (filled during setup(); queryable via reduction_plan()/reduced_dof_audit()).
  t_ReductionPlan                           _reductionPlan;

  //! @brief Auto-generate outgoing-characteristic boundary closure rows for
  //! EVOL_HYPERBOLIC blocks.  Appends collocated INTERIOR/classify=false rows to
  //! _mEqn and returns the count.  Must run after _classify_pde() (consumes
  //! _blockFaceData) and before the interface-plan build.  Gated at the call site
  //! by CRONOS_AUTO_HYP_CLOSURE (environment-only).  See OCEnv_interface_refactor_plan.md.
  size_t _generate_hyperbolic_boundary_closure( bool append = true );   //!< append=false: count the missing rows only
  //! outgoing rows the model already supplied, and at how many faces, in the last generator call
  size_t _hypClosureSkipped = 0, _hypClosureFaces = 0;

  //! @brief Auto-generate boundary-collocation rows for genuinely algebraic
  //! states in non-hyperbolic blocks.  An algebraic state takes no BC, so a
  //! constraint collocated on the spatial interior leaves its boundary DOFs
  //! unequationed; this appends the constraint at those faces.  Runs in the same
  //! single-pass slot as the hyperbolic closure (after _classify_pde, before the
  //! consistency/cache build).  Gated at the call site by CRONOS_AUTO_ALG_CLOSURE (environment-only).
  //! See OCEnv_interface_refactor_plan.md.
  size_t _generate_algebraic_boundary_closure();


public:

  //! @brief The persisted high-index reduction plan built during setup()
  //! (Pantelides + dummy-derivative matching).  Empty for index-0/1 systems.
  t_ReductionPlan const& reduction_plan() const { return _reductionPlan; }


protected:

  //! @brief Auto differential-elimination (Pantelides increment 1; gated by
  //! options.AUTO.DIFF_ELIM, default OFF).  For each rectangular block, an ELIMINABLE
  //! dangling differential column d/d(dir)s -- one referenced by a differential CONSUMER
  //! equation but determined by an algebraic SOURCE equation G that is LINEAR in it --
  //! is removed by solving G=0 for d/d(dir)s and substituting into the consumer, then
  //! RELOCATING the consumer onto G's transverse face (grafting G's face directions into
  //! the consumer's t_Eqn.dom) so the block becomes square.  Solution-exact (uses G=0,
  //! which holds at the root).  Mutates _mEqn in place; iterated to a fixpoint.  Runs
  //! post-reduce_order / pre-collocation-build.  Returns the number of eliminations.
  int _auto_diff_eliminate();

  //! @brief Stage-2 decision layer: assemble the persisted Pantelides reduction
  //! plan from the per-block classification (no equation mutation).
  void _build_reduction_plan();

  //! @brief Stage-2 execution layer: realise _reductionPlan by differentiating
  //! the flagged interior constraints in place (FAD chain-rule + pure-RHS subst).
  bool _reduce_high_index();


protected:

  //! @brief optional per-input interface-continuity declarations. Keys are input variables, values are maps from physical domain variables to the declared InputContinuity level along that domain. Missing input or direction means DISCONTINUOUS (the default). Asserted against supplied data by _check_input_continuity() at eval.
  t_Cont                             _mInpCont;

  //! @brief Import user-facing model data into a fresh private working DAG
  bool _copy_usr_to_local
    ();

  //! @brief Register DAG variables from another already set-up collocation environment
  bool _register
    ( std::map<FFVar::pt_idVar,size_t>& ndxAllVars, std::vector<FFVar>& vAllVarsLoc,
      FFGraph* dagSrc, std::vector<FFVar> const& vCstSrc, t_Dom const& mDomSrc,
      t_Var const& mVarSrc, t_Var const& mInpSrc, t_Eqns const& mEqnSrc, t_Fcts const& mFctSrc,
      std::map<int,t_Symbol> const* bSymSrc = nullptr,
      std::vector<t_AuxDef> const* auxDefSrc = nullptr,
      t_Transitions const* trnSrc = nullptr );

  //! @brief Register DAG variable from another already set-up collocation environment
  bool _register
    ( FFVar const& var, std::vector<FFVar>& all,  std::map<FFVar::pt_idVar,size_t>& ndx );

  //! @brief Helper for retrieving local var
  static FFVar _local
    ( FFVar const& var, std::vector<FFVar> const& all, 
      std::map<FFVar::pt_idVar,size_t> const& ndx );

  //! @brief Helper for copying FFVar vector
  static std::vector<FFVar> _copy
    ( std::vector<FFVar> const& in, std::vector<FFVar> const& all, 
      std::map<FFVar::pt_idVar,size_t> const& ndx );

  //! @brief Helper for copying FFVar set
  static std::set<FFVar,lt_FFVar> _copy
    ( std::set<FFVar,lt_FFVar> const& in, std::vector<FFVar> const& all,
      std::map<FFVar::pt_idVar,size_t> const& ndx );

  //! @brief Helper for copying FFVar map
  template <typename VAL>
  static std::map<FFVar,VAL,lt_FFVar> _copy
    ( std::map<FFVar,VAL,lt_FFVar> const& in, std::vector<FFVar> const& all,
      std::map<FFVar::pt_idVar,size_t> const& ndx );

  //! @brief Helper for copying optional input continuity-declaration maps (remaps the
  //! input and direction FFVar keys to the local DAG; levels copied verbatim).
  static t_Cont _copy_inp_cont
    ( t_Cont const& in, std::vector<FFVar> const& all,
      std::map<FFVar::pt_idVar,size_t> const& ndx );

  //! @brief Helper for copying reference-function maps and wrapping domain-coordinate keys
  static std::map<FFVar,t_Fun,lt_FFVar> _copy
    ( std::map<FFVar,t_Fun,lt_FFVar> const& in, t_Dom const& domSrc,
      std::vector<FFVar> const& all, std::map<FFVar::pt_idVar,size_t> const& ndx );

  //! @brief Helper for copying FFVar multimap
  template <typename VAL>
  static std::multimap<FFVar,VAL,lt_FFVar> _copy
    ( std::multimap<FFVar,VAL,lt_FFVar> const& in, std::vector<FFVar> const& all,
      std::map<FFVar::pt_idVar,size_t> const& ndx );

  //! @brief Helper for copying symbol map
  t_Symbol _copy
    ( t_Symbol const& in, std::vector<FFVar> const& all,
      std::map<FFVar::pt_idVar,size_t> const& ndx );


protected:

  //! @brief resolved addresses of the runtime thread-count setters/getters.
  //! Looked up by NAME (dlsym), never declared, so the header gains no link dependency
  //! and cannot collide with <omp.h> or a vendor BLAS header.
  struct t_ThreadHooks
  {
    void (*omp_set) ( int )         = nullptr;
    int  (*omp_get) ()              = nullptr;
    void (*blas_set)( int )         = nullptr;
    int  (*blas_get)()              = nullptr;
    void (*mkl_set) ( int )         = nullptr;
    int  (*mkl_get) ()              = nullptr;
    void (*blis_set)( int64_t )     = nullptr;
    int64_t (*blis_get)()           = nullptr;
    bool none() const { return !omp_set && !blas_set && !mkl_set && !blis_set; }
  };
  static t_ThreadHooks const& _thread_hooks();

  //! @brief RAII thread cap for Options::MAXTHREAD.  Sets the process-global
  //! thread count on entry to setup()/solve() and restores it on EVERY exit path.
  //! Depth-counted: the w-decide restore re-enters setup() from inside setup(), and only
  //! the outermost guard may act.  MAXTHREAD=0 constructs a guard that does nothing at
  //! all -- no call is made, so the default path is bit-identical to a build without the cap.
  class t_ThreadCap
  {
  public:
    explicit t_ThreadCap( size_t maxthread );
    ~t_ThreadCap();
    t_ThreadCap( t_ThreadCap const& ) = delete;
    t_ThreadCap& operator=( t_ThreadCap const& ) = delete;
    //! @brief the cap in force right now, 0 = none.  Read by the cholmod_common sites,
    //! which are static/const and cannot see options.
    static int active(){ return _active; }
  private:
    bool _counted   = false;      //!< this guard incremented the depth
    bool _outermost = false;      //!< ... and was the one that applied the cap
    inline static int _depth  = 0;
    inline static int _active = 0;
    //! the thread counts are PROCESS-wide, and solves may run concurrently (FFGraph::veval threads, each on its own
    //! copy of an embedded solver): the depth, the capture and the restore are serialised by this mutex, and the
    //! LAST guard out restores (2026-10-03; before, the outermost guard restored, possibly while others still ran)
    static std::mutex& _mutex(){ static std::mutex m; return m; }
    inline static int     _omp0  = 0;   //!< pristine values, captured by the outermost guard
    inline static int     _blas0 = 0;
    inline static int     _mkl0  = 0;
    inline static int64_t _blis0 = 0;
  };


protected:



  //! @brief Evaluate the default reference coordinate for a domain variable.
  double _default_dom_ref
    ( FFVar const& dvar ) const;


public:

  //! @brief Declared continuity of a distributed input across element interfaces, per
  //! domain direction.  Distributed inputs are DISCONTINUOUS by default (independent
  //! per-element DOFs, free to jump at interfaces).  A stronger declaration is an
  //! ASSERTION about the supplied data: over-declaring (claiming smoother than the data
  //! is) is the dangerous direction and is caught by _check_input_continuity() at eval.
  enum InputContinuity
  {
    DISCONTINUOUS = 0, //!< default: input may jump across element interfaces
    CONTINUOUS_C0 = 1, //!< input value matches across interfaces (C0)
    SMOOTH_C1     = 2  //!< input value and first derivative match across interfaces (C1)
  };

  //! @brief Add distributed input, optionally with a reference value used by automatic PDE classification or initialisation
  void add_input
    ( FFVar const& Var, std::vector<FFVar> const& vDom={},
      std::optional<double> ref=std::nullopt )
    {
      auto& mInp           = _usr._mInpUsr;
      auto& mInpDisc       = _usr._mInpDiscUsr;
      auto& classInpRef    = _usr._classInpRefUsr;
      auto& classInpRefFun = _usr._classInpRefFunUsr;
      auto [it,ins]        = mInp.insert( {Var,{}} );
      if( !ins ) it->second.clear();
      for( auto const& d : vDom )
        it->second.insert( d );

      // No explicit input discretisation: use the physical state/domain
      // discretisation in _mDomUsr for this input.
      mInpDisc.erase( Var );

      if( ref ){
        classInpRef[Var] = *ref;
        classInpRefFun.erase( Var );
      }
      _issetup    = false;
      _classified = false;
      _on_model_changed( ModelChange::DERIVATIVES );
    }

  //! @brief Add distributed input with its own local collocation type/order. The physical element partition matches the corresponding state/domain partition, but the input can use a different collocation type and/or number of nodes on each matching element. Passing n_node=0 keeps the node count of the physical domain; the type is set to the value supplied here.
  void add_input
    ( FFVar const& Var, std::vector<FFVar> const& vDom,
      FFDom::TYPE const& type, size_t const n_node=0,
      std::optional<double> ref=std::nullopt )
    {
      add_input( Var, vDom, ref );
      if( vDom.empty() ) return;

      std::map<FFVar,InpColloc,lt_FFVar> disc;
      for( auto const& dvar : vDom ){
        auto itdom = _usr._mDomUsr.find( dvar );
        if( itdom == _usr._mDomUsr.end() )
          throw Exceptions( Exceptions::INDEX );

        // The effective order is validated now, as before, against the domain as currently declared.
        if( ( n_node ? n_node : itdom->second.n_node ) == 0 )
          throw FFDom::Exceptions( FFDom::Exceptions::NODES );
        disc[dvar] = InpColloc{ type, n_node };
      }
      _usr._mInpDiscUsr[Var] = std::move( disc );
    }

  //! @brief Add distributed input with its own local collocation type/order. The physical element partition matches the corresponding state/domain partition, but the input can use a different collocation type and/or number of nodes on each matching element. Passing n_node=0 keeps the node count of the physical domain; the type is set to the value supplied here.
  void add_input
    ( FFVar const& Var, std::vector<FFVar> const& vDom,
      FFDom::TYPE const& type, size_t const n_node, t_Fun const& ref )
    { add_input( Var, vDom, type, n_node ); update_ref( Var, ref ); }

  //! @brief Add distributed input with a per-direction interface-continuity declaration.
  //! vCnt must match vDom in size; entry i declares the continuity of Var across element
  //! interfaces in direction vDom[i] (DISCONTINUOUS default / CONTINUOUS_C0 / SMOOTH_C1).
  //! The input uses the matching state/domain grid; the declaration is an ASSERTION about
  //! the supplied data, checked by _check_input_continuity() at eval (over-declaration is
  //! rejected).  Only non-default (>= CONTINUOUS_C0) entries are recorded.
  void add_input
    ( FFVar const& Var, std::vector<FFVar> const& vDom,
      std::vector<InputContinuity> const& vCnt,
      std::optional<double> ref=std::nullopt )
    {
      add_input( Var, vDom, ref );
      if( vDom.size() != vCnt.size() )
        throw Exceptions( Exceptions::INDEX );
      auto& mInpCont = _usr._mInpContUsr;
      mInpCont.erase( Var );
      std::map<FFVar,int,lt_FFVar> dirs;
      for( size_t i=0; i<vDom.size(); ++i )
        if( vCnt[i] != DISCONTINUOUS )
          dirs[ vDom[i] ] = static_cast<int>( vCnt[i] );
      if( !dirs.empty() )
        mInpCont[Var] = std::move( dirs );
    }

  //! @brief Add lumped input with a reference value used by automatic PDE classification.
  //! With @p is_decision=true the input is also flagged as a reduced-space decision variable (auto-registered as a
  //! control by setup() in declaration order), and @p ref serves as its nominal value.
  void add_input
    ( FFVar const& Var, double const ref, bool is_decision=false )
    { add_input( Var, std::vector<FFVar>{}, std::optional<double>( ref ) ); if( is_decision ) _flag_decision( Var ); }

  //! @brief Add distributed input with a domain-dependent reference function.  With
  //! @p is_decision=true the input is also flagged as a reduced-space decision variable (auto-registered as a
  //! control by setup() in declaration order), and @p ref serves as its nominal profile.
  void add_input
    ( FFVar const& Var, std::vector<FFVar> const& vDom, t_Fun const& ref, bool is_decision=false )
    { add_input( Var, vDom ); update_ref( Var, ref ); if( is_decision ) _flag_decision( Var ); }

  //! @brief One registered control in the canonical control vector: its DOF count in the declared
  //! function space, and its offset (block start) in that vector.
  struct ControlSpec
  {
    size_t ndof   = 0;
    size_t offset = 0;
  };

  //! @brief Register a declared input as a reduced-space control (decision variable).  Equivalent to
  //! add_input( ..., is_decision=true ), for a control decided after the input was declared.  Keyed by
  //! FFVar in a map, so the order is canonical (FFVar order, NOT registration order) and duplicates are
  //! impossible.  The DOF count comes from the DECLARATION, not from any solver -- see control_ndof().
  void register_control
    ( FFVar const& Var )
    { _flag_decision( Var ); _reindex_controls(); _on_controls_changed(); }

  //! @brief Clear the control registry.
  void clear_controls
    ()
    { _decisionFlag.clear(); _controls.clear(); _nControlDof = 0; _on_controls_changed(); }

  //! @brief The registered controls, keyed by input FFVar in canonical FFVar order, each carrying its
  //! DOF count and offset in the control vector.  Access by key: controls().at( myInput ).offset.
  std::map<FFVar,ControlSpec,lt_FFVar> const& controls
    () const
    { return _controls; }

  //! @brief Total control DOFs -- the dimension of the canonical control vector p.
  size_t n_control_dof
    () const
    { return _nControlDof; }

  //! @brief The {offset,ndof} block of @p Var in the control vector; ndof==0 if it is not registered.
  ControlSpec control_block
    ( FFVar const& Var )
    const
    { auto const it = _controls.find( Var ); return it==_controls.cend()? ControlSpec(): it->second; }

  //! @brief The parameter indices a solver's sensitivity runs over, in control order (see ODESLV_BASE);
  //! empty at model level.  fdiff_seed() reads it on the sensitivity model.
  virtual std::vector<size_t> const& sensitivity_index
    () const
    { static std::vector<size_t> const none; return none; }

  //! @brief The parameter indices of the DOFs of declared input @p w in a solver's extracted parameter
  //! vector (see ODESLV_BASE); empty at model level.
  virtual std::vector<size_t> parameter_index
    ( FFVar const& w ) const
    { return std::vector<size_t>(); }

  //! @brief The DOF count @p Var has as a control, from its DECLARED function space alone: the product,
  //! over the domains it is distributed over, of n_elem * n_node -- n_node taken from the input's own
  //! InpColloc where it declares one and from the domain otherwise.  1 for a time-invariant input, 0 if
  //! @p Var is not a declared input.  Registration is not required, and no solver is consulted.
  size_t control_ndof
    ( FFVar const& Var ) const;

  //! @brief The components of control-space VECTOR @p g (size n_control_dof()) belonging to @p Var -- an
  //! adjoint gradient, a step, a bound.  Empty if @p Var is not registered or @p g is not in control space.
  std::vector<double> control_slice
    ( FFVar const& Var, std::vector<double> const& g )
    const
    {
      ControlSpec const s = control_block( Var );
      if( !s.ndof || g.size() != _nControlDof ) return std::vector<double>();
      return std::vector<double>( g.cbegin()+s.offset, g.cbegin()+s.offset+s.ndof );
    }

  //! @brief The columns of ROW-MAJOR ( @p nrow x n_control_dof() ) Jacobian @p Jrm belonging to @p Var,
  //! returned row-major as nrow x ndof.  The strided gather that OCFESLV::sens_jacobian() and ODESLV's
  //! function-gradient block both need; empty on a size or registration mismatch.
  std::vector<double> control_columns
    ( FFVar const& Var, std::vector<double> const& Jrm, size_t const nrow )
    const
    {
      ControlSpec const s = control_block( Var );
      if( !s.ndof || !nrow || Jrm.size() != nrow*_nControlDof ) return std::vector<double>();
      std::vector<double> blk( nrow*s.ndof, 0. );
      for( size_t i=0; i<nrow; ++i )
        for( size_t j=0; j<s.ndof; ++j )
          blk[i*s.ndof+j] = Jrm[i*_nControlDof + s.offset + j];
      return blk;
    }

  //! @brief Hold the declared input @p w at KNOWN values: one value for a time-invariant input, its nodal DOFs
  //! (element-major, as control_ndof() counts them) for a distributed one.  A fixed input is not a parameter of
  //! the extracted model -- the solver substitutes the constants where it would mint levels -- and it cannot be
  //! a control: registering it is undone here.  Returns false if @p w is not a declared input or the count is
  //! wrong.
  //! @brief Whether a solver applies fixed-input values at every call (true: OCFESLV, so fix_input works at any
  //! time) or substitutes them at setup() (false: ODESLV, so fix_input after setup() would silently do nothing).
  virtual bool _fixed_inputs_applied_per_call
    () const
    { return true; }

  bool fix_input
    ( FFVar const& w, std::vector<double> const& values )
    {
      if( _issetup && !_fixed_inputs_applied_per_call() ){
        std::cerr << "FFModel::fix_input ** " << w.name() << ": this solver substitutes fixed values at setup(), so a"
                     " fix_input after setup() would have no effect -- REFUSED.  Call it before setup(), or run setup()"
                     " again; to vary a value between solves, give it as an input instead (map 2, or solve by name).\n";
        return false;
      }
      if( _usr._mInpUsr.find( w ) == _usr._mInpUsr.cend() ) return false;
      size_t const ndof = control_ndof( w );
      if( values.size() != ndof ) return false;
      _usr._mInpFixUsr[ w ] = values;
      // a fixed input is never a control: clear the FLAG too, since _reindex_controls() rebuilds the
      // registry from _decisionFlag at setup()
      bool changed = _controls.erase( w ) > 0;
      for( auto it = _decisionFlag.begin(); it != _decisionFlag.end(); )
        if( it->id().second == w.id().second ){ it = _decisionFlag.erase( it ); changed = true; } else ++it;
      if( changed ){ _reindex_controls(); _on_controls_changed(); }
      return true;
    }
  //! @brief The fixed values of input @p w, or nullptr if it is not fixed
  std::vector<double> const* fixed_input
    ( FFVar const& w ) const
    { auto const it = _usr._mInpFixUsr.find( w ); return it == _usr._mInpFixUsr.cend()? nullptr: &it->second; }

  //! @brief One DOF of an input: its element index and node index on each of the input's domains (both maps
  //! empty for a time-invariant input).  The element map has the shape of pos_input's ndx_el.
  struct DofIndex
  {
    std::map<FFVar,size_t,lt_FFVar> element;
    std::map<FFVar,size_t,lt_FFVar> node;
  };

  //! @brief A declared input and the variables for its DOFs, as a FLAT VECTOR in control_dofs() order, or a
  //! GENERATOR called once per DofIndex -- the form that never requires knowing the order, for multi-domain
  //! inputs in particular.
  struct InputArg
  {
    FFVar                                      input;
    std::vector<FFVar>                         dofs;
    std::function<FFVar(DofIndex const&)>      gen;
    InputArg( FFVar const& u, std::vector<FFVar> const& v ) : input( u ), dofs( v ) {}
    InputArg( FFVar const& u, std::function<FFVar(DofIndex const&)> g ) : input( u ), gen( std::move( g ) ) {}
  };

  //! @brief A declared input or constant and its VALUES: a flat list in control_dofs() order, or a generator called
  //! once per DofIndex.  A constant takes one value.  The value counterpart of InputArg, for the solvers' own calls.
  struct InputVal
  {
    FFVar                                      input;
    std::vector<double>                        vals;
    std::function<double(DofIndex const&)>     gen;
    InputVal( FFVar const& u, std::vector<double> const& v ) : input( u ), vals( v ) {}
    InputVal( FFVar const& u, std::function<double(DofIndex const&)> g ) : input( u ), gen( std::move( g ) ) {}
  };

  //! @brief The DOFs of declared input @p u, ENUMERATED in the order the solvers lay out its control block:
  //! element blocks outer, nodes inner, and at each level the first domain (canonical FFVar order) varying
  //! fastest -- as OCFESLV::get_input_values / sample_input traverse it, and as ODESLV mints its levels.
  //! control_dofs( u ).size() == control_ndof( u ).  Empty if @p u is not a declared input.
  std::vector<DofIndex> control_dofs
    ( FFVar const& u )
    const;

  //! @brief Resolve an InputArg into the DOF variables of its input, in control_dofs() order: the flat vector
  //! as given, or the generator called on each DofIndex.  Refuses (false, reason in @p err) a count that does
  //! not match control_ndof(), an undeclared input, or a generator that throws on some DofIndex.
  bool resolve_input_arg
    ( InputArg const& a, std::vector<FFVar>& dofs, std::string& err )
    const;

  //! @brief Declare into @p sens the FORWARD SENSITIVITY MODEL of this model along the declared input @p u,
  //! in Gateaux form: a new distributed input du carrying u's domains and collocation, one sensitivity state
  //! s_x per state x (same domains), every equation E joined by its directional derivative along (s, du),
  //! every output F by dF -- with domains, masks, options, output kind and coordinates carried over
  //! verbatim.  Outputs are laid out as the nf originals followed by the nf derivatives, so output nf+j is dF_j.  All other inputs, constants and domain variables are held fixed.  @p sens must be freshly
  //! constructed on the SAME user DAG as this model; it may be any FFModel-derived solver, and its own
  //! add_* overrides fire for every declaration.  Returns false, with a message in @p err, if @p u is not a
  //! declared input or a derivative could not be formed.  The new variables are returned through @p du and
  //! @p vS so the caller can set du's values (a unit vector per DOF) and read s.
  bool fdiff
    ( FFModel& sens, FFVar const& u, FFVar& du, std::vector<FFVar>& vS, std::string& err )
    const
    {
      std::vector<std::vector<FFVar>> vDU, vSS;
      if( !fdiff( sens, std::vector<FFVar>{ u }, 1, vDU, vSS, err ) ) return false;
      du = vDU[0][0];  vS = vSS[0];  return true;
    }

  //! @brief Multi-direction form.  Direction k perturbs EVERY control in @p vU at once through its own input
  //! du_c^(k) (declared in control c's function space and REGISTERED as a control of @p sens), with its own
  //! sensitivity states s^(k) and outputs dF^(k).  @p sens gets nx*(1+nDir) states and nf*(1+nDir) outputs,
  //! laid out [F | dF^(1) | ... | dF^(nDir)] (the F block omitted when @p keep_originals is false, giving the
  //! legacy layout f[k*nf+j]).  @p vDU[k][c] and @p vS[k][i] return the new variables; seed a
  //! direction with fdiff_seed().  nDir == 0 means one direction per DOF of the listed controls (the full
  //! Jacobian in one solve).
  bool fdiff
    ( FFModel& sens, std::vector<FFVar> const& vU, size_t nDir,
      std::vector<std::vector<FFVar>>& vDU, std::vector<std::vector<FFVar>>& vS, std::string& err,
      bool const keep_originals = true,
      std::vector<FFVar> const& vCdir = std::vector<FFVar>(),
      std::vector<std::vector<double>> const& seedC = std::vector<std::vector<double>>() )
    const;

  //! @brief Build the parameter vector of the sensitivity model @p sens (after its setup()) from a parameter
  //! vector @p P0 of the ORIGINAL model @p orig (after its setup()), with every direction k seeded by its k-th
  //! DOF unit vector -- the full Jacobian in one solve.  Each original input's DOFs are placed by the two
  //! models' parameter_index() layouts, so no ordering assumption is made: a time-invariant direction input
  //! sorts BEFORE the distributed originals in the extracted order, and this handles that.
  static std::vector<double> fdiff_seed_all
    ( FFModel const& sens, std::vector<std::vector<FFVar>> const& vDU, FFModel const& orig, std::vector<double> const& P0 )
    {
      size_t np = 0;                                             // the product's parameter count, from its layout
      for( auto const& [w, dom] : sens.var_declared_input() ) for( size_t i : sens.parameter_index( w ) ) np = std::max( np, i+1 );
      std::vector<double> P( np, 0. );
      for( auto const& [w, dom] : orig.var_declared_input() ){
        auto const io = orig.parameter_index( w ), is = sens.parameter_index( w );
        for( size_t j = 0; j < io.size() && j < is.size(); ++j ) if( io[j] < P0.size() ) P[ is[j] ] = P0[ io[j] ];
      }
      for( size_t k = 0; k < vDU.size(); ++k ){ std::vector<double> Pk( P.size(), 0. ); fdiff_seed( sens, vDU, k, Pk );
        for( size_t i = 0; i < P.size(); ++i ) if( Pk[i] != 0. ) P[i] = Pk[i]; }
      return P;
    }

  //! @brief In a parameter vector @p P of the sensitivity model @p sens (after its setup()), zero every
  //! direction input and set direction @p k to the k-th DOF unit vector, DOFs enumerated control by control
  //! in the order of @p vDU[k].  Returns false if @p k has no DOF.
  static bool fdiff_seed
    ( FFModel const& sens, std::vector<std::vector<FFVar>> const& vDU, size_t const k, std::vector<double>& P )
    {
      auto const& ndx = sens.sensitivity_index();
      for( auto const& dir : vDU )
        for( auto const& du : dir ){
          ControlSpec const b = sens.control_block( du );
          for( size_t j = 0; j < b.ndof; ++j ) if( b.offset+j < ndx.size() ) P[ ndx[b.offset+j] ] = 0.;
        }
      if( k >= vDU.size() ) return false;
      size_t dof = k;                                          // the k-th DOF across the controls of direction k
      for( auto const& du : vDU[k] ){
        ControlSpec const b = sens.control_block( du );
        if( dof < b.ndof ){ if( b.offset+dof < ndx.size() ){ P[ ndx[b.offset+dof] ] = 1.; return true; } return false; }
        dof -= b.ndof;
      }
      return false;
    }

protected:
  //! @brief Hook in the _on_* family: the control registry or its layout has changed.  A solver overrides
  //! it to drop whatever it caches against the control vector -- OCFESLV invalidates its reduced-space
  //! dependency pattern here.  NOT called from setup(), where solvers already invalidate wholesale.
  virtual void _on_controls_changed
    ()
    {}

  //! @brief Recompute every registered control's DOF count from the declaration, and its offset from the
  //! canonical (map) order.  Called by register_control() and by setup(), where the domains a control is
  //! distributed over are certainly declared -- an add_input( ..., is_decision=true ) may well precede
  //! its own add_domain(), and control_ndof() would report 0 at that point.
  void _reindex_controls
    ()
    {
      _controls.clear();
      for( auto const& u : _decisionFlag ) _controls[u].ndof = control_ndof( u );
      size_t off = 0;
      for( auto& [u,spec] : _controls ){ spec.offset = off; off += spec.ndof; }
      _nControlDof = off;
    }

  //! @brief Flag an input for automatic control registration at setup() (declaration order,
  //! deduplicated).  Used by the is_decision add_input() overloads.
  //! PROTECTED -- no corpus driver calls it and add_input() sets it; a PUBLIC member named with a
  //! leading underscore contradicted the convention either way.
  void _flag_decision
    ( FFVar const& Var )
    { for( auto const& v : _decisionFlag ) if( v.id() == Var.id() ) return; _decisionFlag.push_back( Var ); }

public:   // rev269: the DECLARATION API resumes here -- add_input/reset_input are what a driver calls.

  //! @brief Add distributed inputs, optionally with reference values used by automatic PDE classification
  void add_input
    ( std::vector<FFVar> const& vVar, std::vector<FFVar> const& vDom={},
      std::vector<double> const& ref=std::vector<double>() )
    {
      assert( ref.empty() || ref.size() == vVar.size() );
      for( size_t i=0; i<vVar.size(); ++i ){
        std::optional<double> optref;
        if( !ref.empty() ) optref = ref[i];
        add_input( vVar[i], vDom, optref );
      }
    }

  //! @brief Reset distributed inputs
  void reset_input
    ()
    {
      auto& mInp           = _usr._mInpUsr;
      auto& mInpDisc       = _usr._mInpDiscUsr;
      auto& classInpRef    = _usr._classInpRefUsr;
      auto& classInpRefFun = _usr._classInpRefFunUsr;
      mInp.clear();
      mInpDisc.clear();
      classInpRef.clear();
      classInpRefFun.clear();
      _issetup    = false;
      _classified = false;
      _on_model_changed( ModelChange::DERIVATIVES );
    }


protected:

  //! @brief Inputs flagged is_decision at add_input() time (in declaration order).  setup()
  //! auto-registers these as controls, so a driver can declare "this input is a decision
  //! variable" alongside its nominal and skip the separate register_control() call.
  std::vector<FFVar>                        _decisionFlag;

  //! @brief Canonical control layout: the registered controls in FFVar order, each with its DOF count
  //! and offset.  Rebuilt by _reindex_controls(); the single definition every solver reads.
  std::map<FFVar,ControlSpec,lt_FFVar>      _controls;

  //! @brief Dimension of the canonical control vector (sum of the controls' ndof).
  size_t                                    _nControlDof = 0;

};


inline size_t
FFModel::control_ndof
( FFVar const& Var )
const
{
  // The DECLARATION side throughout: _mInpUsr / _mDomUsr / _mInpDiscUsr are keyed by the modeller's own
  // FFVars, exactly as _decisionFlag is.  The post-setup maps hold DAG COPIES, so a lookup there would
  // compare mismatched DAG pointers and miss -- the same trap OCFESLV::_is_registered_control documents.
  auto const iti = _usr._mInpUsr.find( Var );
  if( iti == _usr._mInpUsr.cend() ) return 0;           // not a declared input
  if( iti->second.empty() )         return 1;           // time-invariant: a single value
  auto const itdisc = _usr._mInpDiscUsr.find( Var );
  size_t ndof = 1;
  for( auto const& dvar : iti->second ){
    auto const itd = _usr._mDomUsr.find( dvar );
    if( itd == _usr._mDomUsr.cend() ) return 0;         // domain not declared (yet)
    size_t n_node = itd->second.n_node;                 // InpColloc n_node==0 means "as the domain"
    if( itdisc != _usr._mInpDiscUsr.cend() ){
      auto const itc = itdisc->second.find( dvar );
      if( itc != itdisc->second.cend() && itc->second.n_node ) n_node = itc->second.n_node;
    }
    ndof *= itd->second.n_elem * n_node;                // no continuity: every element node is a DOF
  }
  return ndof;
}

inline bool
FFModel::fdiff
( FFModel& sens, std::vector<FFVar> const& vU, size_t nDir,
  std::vector<std::vector<FFVar>>& vDU, std::vector<std::vector<FFVar>>& vS, std::string& err,
  bool const keep_originals, std::vector<FFVar> const& vCdir, std::vector<std::vector<double>> const& seedC )
const
{
  std::lock_guard<std::recursive_mutex> dag_lock_( dag_mutex() );   // see dag_mutex()
  if( !_usr._mTrnUsr.empty() ){                                    // see add_transition
    err = "fdiff: the model has transitions, which the differentiated model does not carry yet -- refused";
    return false;
  }
  // CONSTANT directions (vCdir, seeds seedC[k][c]): a constant needs no direction input -- it enters the
  // directional derivative with a literal seed.  The declared constants stay constants of the product.
  for( auto const& c : vCdir ){
    bool found = false;
    for( auto const& cd : _usr._vCstUsr ) if( cd.id().second == c.id().second ){ found = true; break; }
    if( !found ){ err = "fdiff: " + c.name() + " is not a declared constant"; return false; }
  }
  if( !vCdir.empty() && seedC.size() != nDir ){ err = "fdiff: one constant-seed row per direction is required"; return false; }
  for( auto const& row : seedC ) if( row.size() != vCdir.size() ){ err = "fdiff: a constant-seed row has the wrong size"; return false; }
  err.clear();  vDU.clear();  vS.clear();
  FFGraph* const dag = _usr._dagUsr;
  if( !dag ){ err = "fdiff: this model has no user DAG (call set() first)"; return false; }
  if( sens._usr._dagUsr && sens._usr._dagUsr != dag ){ err = "fdiff: the target model is on a different DAG"; return false; }
  if( !sens._usr._dagUsr ) sens.set( dag );
  if( vU.empty() && vCdir.empty() ){ err = "fdiff: no direction given"; return false; }
  size_t ndof_total = 0;
  for( auto const& u : vU ){
    if( _usr._mInpUsr.find( u ) == _usr._mInpUsr.cend() ){ err = "fdiff: " + u.name() + " is not a declared input"; return false; }
    ndof_total += control_ndof( u );
  }
  if( !nDir ) nDir = ndof_total;                          // one direction per DOF: the full Jacobian in one solve
  if( !nDir ){ err = "fdiff: zero directions"; return false; }

  // ---- 1. constants and domains, verbatim -----------------------------------------------------------
  if( !_usr._vCstUsr.empty() ) sens.set_constant( _usr._vCstUsr, _usr._vCstValUsr );
  for( auto const& [dv, dom] : _usr._mDomUsr ){
    auto const itr = _usr._classDomRefUsr.find( dv );
    if( itr != _usr._classDomRefUsr.cend() ) sens.add_domain( std::vector<FFVar>{ dv }, dom, std::vector<double>{ itr->second } );
    else                                     sens.add_domain( std::vector<FFVar>{ dv }, dom );
  }
  if( _usr._evolution_dom_setUsr ) sens.set_evolution_domain( _usr._evolution_dom_varUsr );

  // ---- 2. states, then nDir sensitivity copies ------------------------------------------------------
  std::vector<FFVar> vX;  std::vector<std::vector<FFVar>> vXdom;
  for( auto const& [x, sdom] : _usr._mVarUsr ){
    std::vector<FFVar> vDom( sdom.cbegin(), sdom.cend() );
    auto const itr = _usr._classVarRefUsr.find( x );
    if( itr != _usr._classVarRefUsr.cend() ) sens.add_state( x, vDom, itr->second );
    else                                     sens.add_state( x, vDom );
    auto const itf = _usr._classVarRefFunUsr.find( x );      // function references too (OCFESLV starts Newton there)
    if( itf != _usr._classVarRefFunUsr.cend() ) sens.update_ref( x, itf->second );
    vX.push_back( x );  vXdom.push_back( vDom );
  }
  vS.assign( nDir, std::vector<FFVar>() );
  for( size_t k = 0; k < nDir; ++k )
    for( size_t i = 0; i < vX.size(); ++i ){
      FFVar s = dag->add_var( "s" + std::to_string( k+1 ) + "[" + vX[i].name() + "]" );
      sens.add_state( s, vXdom[i], 0. );
      vS[k].push_back( s );
    }

  // ---- 3. inputs verbatim, then nDir direction inputs per control, each registered as a control -----
  auto replay_input = [&]( FFVar const& w, FFVar const& as, bool zero_ref ){
    auto const& sdom = _usr._mInpUsr.at( w );
    std::vector<FFVar> vDom( sdom.cbegin(), sdom.cend() );
    auto const itd = _usr._mInpDiscUsr.find( w );
    if( vDom.empty() || itd == _usr._mInpDiscUsr.cend() || itd->second.empty() ) sens.add_input( as, vDom );
    else{ auto const& col = itd->second.cbegin()->second; sens.add_input( as, vDom, col.type, col.n_node ); }
    auto const itc = _usr._mInpContUsr.find( w );
    if( itc != _usr._mInpContUsr.cend() ) sens._usr._mInpContUsr[ as ] = itc->second;
    if( zero_ref ) sens.update_ref( as, 0. );
    else{ auto const itr = _usr._classInpRefUsr.find( w ); if( itr != _usr._classInpRefUsr.cend() ) sens.update_ref( as, itr->second );
          auto const itf = _usr._classInpRefFunUsr.find( w ); if( itf != _usr._classInpRefFunUsr.cend() ) sens.update_ref( as, itf->second ); }
  };
  for( auto const& [w, sdom] : _usr._mInpUsr ) replay_input( w, w, false );
  vDU.assign( nDir, std::vector<FFVar>() );
  for( size_t k = 0; k < nDir; ++k )
    for( auto const& u : vU ){
      FFVar du = dag->add_var( "d" + std::to_string( k+1 ) + "[" + u.name() + "]" );
      replay_input( u, du, true );
      sens.register_control( du );                             // so its block is discoverable for seeding
      vDU[k].push_back( du );
    }

  // ---- 4. the directional derivative of direction k: independents {x..., u...}, direction {s^(k)..., du^(k)...}
  auto ddir = [&]( FFVar const& f, size_t const k, FFVar& df ) -> bool {
    std::vector<FFVar> vIndep( vX ), vDir( vS[k] );
    for( size_t c = 0; c < vU.size(); ++c ){ vIndep.push_back( vU[c] ); vDir.push_back( vDU[k][c] ); }
    for( size_t c = 0; c < vCdir.size(); ++c ){ vIndep.push_back( vCdir[c] ); vDir.push_back( FFVar( seedC[k][c] ) ); }
    try{ FFVar* pd = dag->DFAD( 1, &f, (unsigned)vIndep.size(), vIndep.data(), vDir.data() ); df = pd[0]; delete[] pd; return true; }
    catch( FFBase::Exceptions& e ){ err = std::string("FFBase::Exceptions: ") + e.what(); return false; }
    catch( std::exception& e ){ err = std::string("exception: ") + e.what(); return false; }
    catch( ... ){ err = "unknown exception"; return false; }
  };

  // ---- 5. equations, each joined by its nDir derivatives -- same domains, masks and options -----------
  for( auto const& eq : _usr._mEqnUsr ){
    std::vector<FFVar> vDom;  std::vector<int> vLim;
    for( auto const& [dv, lim] : eq.dom ){ vDom.push_back( dv ); vLim.push_back( lim ); }
    sens.add_equation( eq.var, vDom, vLim, *eq.opt );
    for( size_t k = 0; k < nDir; ++k ){
      FFVar deq;
      if( !ddir( eq.var, k, deq ) ){
        FFSubgraph sg = dag->subgraph( 1, &eq.var );  std::ostringstream os;  os << FFExpr::subgraph( dag, sg )[0];
        err = "fdiff: cannot differentiate equation `" + os.str() + "`: " + err;  return false; }
      sens.add_equation( deq, vDom, vLim, *eq.opt );
    }
  }

  // ---- 6. outputs: the nf originals, then dF^(1), ..., dF^(nDir) -- each block nf long -----------------
  auto replay_output = [&]( FFVar const& f, t_Fct const& F ){
    if( F.kind == FctKind::DISTRIBUTED ){
      std::vector<FFVar> vDom;  std::vector<int> vLim;
      for( auto const& [dv, lim] : F.grid ){ vDom.push_back( dv ); vLim.push_back( lim ); }
      sens.add_output( f, vDom, vLim );
    }
    else if( !F.point.empty() ){                            // the SIDE travels with the point
      std::vector<FFVar> vDom;  std::vector<double> vVal;  std::vector<int> vSide;
      for( auto const& [dv, c] : F.point ){
        vDom.push_back( dv ); vVal.push_back( c );
        auto const is = F.side.find( dv ); vSide.push_back( is == F.side.end()? (int)FFDom::MINUS: is->second );
      }
      sens.add_output( f, vDom, vVal, vSide );
    }
    else sens.add_output( f );
  };
  if( keep_originals ) for( auto const& Fo : _usr._mFctUsr ) replay_output( Fo.var, Fo );
  for( size_t k = 0; k < nDir; ++k )
    for( auto const& Fo : _usr._mFctUsr ){
      FFVar dF;
      if( !ddir( Fo.var, k, dF ) ){ err = "fdiff: cannot differentiate an output: " + err; return false; }
      replay_output( dF, Fo );
    }
  return true;
}

inline std::vector<FFModel::DofIndex>
FFModel::control_dofs
( FFVar const& u )
const
{
  std::vector<DofIndex> out;
  auto const iti = _usr._mInpUsr.find( u );
  if( iti == _usr._mInpUsr.cend() ) return out;
  if( iti->second.empty() ){ out.emplace_back(); return out; }        // time-invariant: one DOF
  auto const itdisc = _usr._mInpDiscUsr.find( u );
  std::vector<FFVar> doms;  std::vector<size_t> nel, nno;
  for( auto const& d : iti->second ){                                  // canonical order: first varies fastest
    auto const itd = _usr._mDomUsr.find( d );
    if( itd == _usr._mDomUsr.cend() ) return std::vector<DofIndex>();
    size_t nn = itd->second.n_node;                                    // as control_ndof resolves it
    if( itdisc != _usr._mInpDiscUsr.cend() ){
      auto const itc = itdisc->second.find( d );
      if( itc != itdisc->second.cend() && itc->second.n_node ) nn = itc->second.n_node;
    }
    doms.push_back( d );  nel.push_back( itd->second.n_elem );  nno.push_back( nn );
  }
  size_t nblk = 1;  for( size_t e : nel ) nblk *= e;
  size_t blk  = 1;  for( size_t n : nno ) blk  *= n;
  out.reserve( nblk*blk );
  for( size_t lin = 0; lin < nblk; ++lin )                             // element blocks outer
    for( size_t nk = 0; nk < blk; ++nk ){                              // nodes inner
      DofIndex x;
      size_t r = lin, q = nk;
      for( size_t j = 0; j < doms.size(); ++j ){
        x.element[ doms[j] ] = nel[j]? r % nel[j]: 0;  if( nel[j] ) r /= nel[j];
        x.node   [ doms[j] ] = nno[j]? q % nno[j]: 0;  if( nno[j] ) q /= nno[j];
      }
      out.push_back( x );
    }
  return out;
}

inline bool
FFModel::resolve_input_arg
( InputArg const& a, std::vector<FFVar>& dofs, std::string& err )
const
{
  err.clear();  dofs.clear();
  size_t const nd = control_ndof( a.input );
  if( !nd ){ err = a.input.name() + " is not a declared input of the model (or one of its domains is not declared)"; return false; }
  if( a.gen ){
    for( auto const& d : control_dofs( a.input ) ){
      try{ dofs.push_back( a.gen( d ) ); }
      catch( std::exception const& e ){ err = "the generator for " + a.input.name() + " failed at one of its DOFs: " + e.what(); return false; }
    }
  }
  else dofs = a.dofs;
  if( dofs.size() != nd ){
    err = a.input.name() + " has " + std::to_string( nd ) + " DOFs (the product of n_elem x n_node over its domains), but "
        + std::to_string( dofs.size() ) + " variables were given";
    return false;
  }
  return true;
}

inline void
FFModel::_reset_model
()
{
  _vCstVal.clear();
  _vCst.clear();
  _mInp.clear();
  _mVar.clear();
  _mEqn.clear();
  _mFct.clear();
  _mTrn.clear();
  _trnLift.clear();
  _mDom.clear();
  _classVarRef.clear();
  _classInpRef.clear();
  _classDomRef.clear();
  _classVarRefFun.clear();
  _classInpRefFun.clear();
  _issetup    = false;
  _classified = false;
}

inline bool
FFModel::_deep_copy_model
( FFModel const& src, std::vector<FFVar> const& extraRoots,
  std::vector<FFVar>& vAllVarsLoc, std::map<FFVar::pt_idVar,size_t>& ndxAllVars )
{
  _reset_model();
  _usr.clear();
  _release_working_dag();

  options   = src.options;
  _usr      = src._usr;       // shallow copy: the user DAG is not owned

  _dag      = new FFGraph;
  _dagOwned = true;

  // Registration order (the solver's former _register): constants, domains, states, inputs, equations,
  // outputs, then the solver's extra roots; a variable is registered once, keyed by its source id.
  std::vector<FFVar> vAllVarsSrc;
  ndxAllVars.clear();
  vAllVarsLoc.clear();
  auto reg = [&]( FFVar const& v ){
    if( !v.dag() ) return;
    auto const id = v.id();
    if( ndxAllVars.find( id ) != ndxAllVars.end() ) return;
    ndxAllVars[id] = vAllVarsSrc.size();
    vAllVarsSrc.push_back( v );
  };
  for( auto const& cst : src._vCst )        reg( cst );
  for( auto const& [var, dom] : src._mDom ) reg( var );
  for( auto const& [var, dom] : src._mVar ) reg( var );
  for( auto const& [var, dom] : src._mInp ) reg( var );
  for( auto const& eqn : src._mEqn )        reg( eqn.var );
  for( auto const& fct : src._mFct )        reg( fct.var );
  for( auto const& tr : src._mTrn ){ for( auto const& v : tr.left ) reg( v ); for( auto const& v : tr.right ) reg( v ); reg( tr.dom ); }
  for( auto const& v : extraRoots )         reg( v );
  try{
    _dag->insert( src._dag, vAllVarsSrc, vAllVarsLoc );
  }
  catch(...){
    return false;
  }

  auto loc  = [&]( FFVar const& v ) -> FFVar {
    if( !v.dag() ) return v;
    auto it = ndxAllVars.find( v.id() );
    if( it == ndxAllVars.end() )
      throw Exceptions( Exceptions::INDEX );
    return vAllVarsLoc.at( it->second );
  };
  auto cset = [&]( std::set<FFVar,lt_FFVar> const& in ){
    std::set<FFVar,lt_FFVar> out;
    for( auto const& v : in ) out.insert( loc( v ) );
    return out;
  };
  auto cmap = [&]( auto const& in ){
    typename std::decay<decltype(in)>::type out;
    for( auto const& kv : in ) out.insert( { loc( kv.first ), kv.second } );
    return out;
  };
  // Reference functions: wrapped so that coordinates keyed by the source domain variables are supplied
  // (identical to the solver's historical wrapper).
  auto cfun = [&]( std::map<FFVar,t_Fun,lt_FFVar> const& in ){
    std::vector<std::pair<FFVar,FFVar>> dom_src_dst;
    dom_src_dst.reserve( src._mDom.size() );
    for( auto const& [dsrc, dom] : src._mDom ){
      (void)dom;
      dom_src_dst.emplace_back( dsrc, loc( dsrc ) );
    }
    std::map<FFVar,t_Fun,lt_FFVar> out;
    for( auto const& kv : in ){
      FFVar const& vsrc = kv.first;
      t_Fun const fun = kv.second;
      FFVar const vdst = loc( vsrc );
      out[vdst] = [fun, dom_src_dst]( t_Coord const& coord_dst ) -> double
      {
        t_Coord coord = coord_dst;
        for( auto const& [dsrc,ddst] : dom_src_dst ){
          auto it = coord_dst.find( ddst );
          if( it != coord_dst.end() ) coord[dsrc] = it->second;
        }
        return fun( coord );
      };
    }
    return out;
  };

  _vCstVal = src._vCstVal;
  _vCst.clear();
  _vCst.reserve( src._vCst.size() );
  for( auto const& v : src._vCst ) _vCst.push_back( loc( v ) );
  _mInp.clear();
  for( auto const& [var, dom] : src._mInp )
    _mInp.insert( { loc( var ), cset( dom ) } );
  _mVar.clear();
  for( auto const& [var, dom] : src._mVar )
    _mVar.insert( { loc( var ), cset( dom ) } );
  _mEqn.clear();
  _mEqn.reserve( src._mEqn.size() );
  for( auto const& eqn : src._mEqn )
    // rev311: CLONE, do not share.  Copying the shared_ptr made the source and the copy point at one options
    // object -- a latent coupling between environments that are meant to be independent.
    _mEqn.push_back( { loc( eqn.var ), cmap( eqn.dom ), eqn.opt->clone() } );
  _mTrn.clear();
  for( auto const& tr : src._mTrn ){
    t_Transition ltr;  ltr.dom = loc( tr.dom );  ltr.tau = tr.tau;
    for( auto const& v : tr.left ) ltr.left.push_back( loc( v ) );
    for( auto const& v : tr.right ) ltr.right.push_back( loc( v ) );
    _mTrn.push_back( ltr );
  }
  _mFct.clear();
  _mFct.reserve( src._mFct.size() );
  for( auto const& fct : src._mFct ){
    t_Fct lfct;
    lfct.var   = loc( fct.var );
    lfct.kind  = fct.kind;
    lfct.point = cmap( fct.point );
    lfct.grid  = cmap( fct.grid );
    lfct.row0  = fct.row0;
    lfct.nrow  = fct.nrow;
    for( auto const& [dv, sd] : fct.side ) lfct.side[ loc( dv ) ] = sd;   // the SIDE of a point output
    _mFct.push_back( std::move( lfct ) );
  }
  _mDom.clear();
  for( auto const& [var, dom] : src._mDom )
    _mDom.insert( { loc( var ), dom } );
  _classVarRef    = cmap( src._classVarRef );
  _classInpRef    = cmap( src._classInpRef );
  _classDomRef    = cmap( src._classDomRef );
  _classVarRefFun = cfun( src._classVarRefFun );
  _classInpRefFun = cfun( src._classInpRefFun );
  return true;
}

inline void
FFModel::_normalise_equation_options
( EqnOptions& opt, std::map<FFVar,int,lt_FFVar> const& eqndom )
const
{
  if( opt.role == EqnRole::AUTO ){
    opt.role_auto = true;
    opt.role = _auto_role( eqndom, _usr._evolution_dom_setUsr? &_usr._evolution_dom_varUsr: nullptr );
  }

  // Trace equations are not rows of the volume principal symbol.  (Solver: under weak
  // imposition, OCFESLV::_resolve_interface_type() gives an IC_AUTO trace equation IC_VALUE.)
  if( _ordinary_trace_role( opt.role ) )
    opt.participate_in_classification = false;

  if( opt.role == EqnRole::INTERFACE ){
    // Covered by the trace rule above; receive_sat and donation stay as supplied.
    opt.participate_in_classification = false;
  }

  if( opt.role == EqnRole::DIAGNOSTIC ){
    opt.participate_in_classification = false;
    // rev312: DIAGNOSTIC also receives no SAT and is never donated -- the SOLVER's rules, applied by its
    // override of this function.
  }

  // LINK equations (pointwise auxiliary definitions from reduce_order()) may enter
  // classification and keep receive_sat as supplied; they are never donated.  (Solver:
  // OCFESLV::_eval_eqn() routes weak interface jumps through LINK rows.)
  // rev312: LINK equations are never donated -- a SOLVER rule, applied by its override.

  // (Solver: INITIAL, BOUNDARY and INTERFACE equations enter the exact-imposition donor graph
  // only when donate_for_state_continuity permits.)

  // SURFACE is intentionally not treated as an ordinary BOUNDARY role: it is
  // a PDE on a lower-dimensional manifold and may have its own block.
}

inline void
FFModel::update_ref
( FFVar const& var, double const& val )
{
  auto same_var = []( FFVar const& a, FFVar const& b ) -> bool
  { return a.id() == b.id(); };

  auto& mVar         = _usr._mVarUsr;
  auto& mInp         = _usr._mInpUsr;
  auto& mDom         = _usr._mDomUsr;
  auto& vCst         = _usr._vCstUsr;
  auto& vCstVal      = _usr._vCstValUsr;
  auto& classStateRef= _usr._classVarRefUsr;
  auto& classInpRef  = _usr._classInpRefUsr;
  auto& classDomRef  = _usr._classDomRefUsr;
  auto& classStateRefFun = _usr._classVarRefFunUsr;
  auto& classInpRefFun   = _usr._classInpRefFunUsr;

  bool found = false;

  for( auto const& [v,dom] : mVar ){
    (void)dom;
    if( same_var( v, var ) ){
      classStateRef[v] = val;
      classStateRefFun.erase( v );
      found = true;
      break;
    }
  }

  if( !found ){
    for( auto const& [v,dom] : mInp ){
      (void)dom;
      if( same_var( v, var ) ){
        classInpRef[v] = val;
        classInpRefFun.erase( v );
        found = true;
        break;
      }
    }
  }

  if( !found ){
    for( size_t i = 0; i < vCst.size(); ++i ){
      if( same_var( vCst[i], var ) ){
        if( vCstVal.empty() )
          vCstVal.assign( vCst.size(), 0. );
        if( vCstVal.size() != vCst.size() )
          throw Exceptions( Exceptions::CSTVAL );
        vCstVal[i] = val;
        found = true;
        break;
      }
    }
  }

  if( !found ){
    for( auto const& [v,dom] : mDom ){
      (void)dom;
      if( same_var( v, var ) ){
        classDomRef[v] = val;
        found = true;
        break;
      }
    }
  }

  if( !found ) throw Exceptions( Exceptions::INDEX );

  _issetup    = false;
  _classified = false;
  _on_model_changed( ModelChange::DERIVATIVES );
}

inline void
FFModel::update_ref
( FFVar const& var, t_Fun const& fun )
{
  auto same_var = []( FFVar const& a, FFVar const& b ) -> bool
  { return a.id() == b.id(); };

  bool found = false;

  for( auto const& [v,dom] : _usr._mVarUsr ){
    (void)dom;
    if( same_var( v, var ) ){
      _usr._classVarRefFunUsr[v] = fun;
      _usr._classVarRefUsr.erase( v );
      found = true;
      break;
    }
  }

  if( !found ){
    for( auto const& [v,dom] : _usr._mInpUsr ){
      (void)dom;
      if( same_var( v, var ) ){
        _usr._classInpRefFunUsr[v] = fun;
        _usr._classInpRefUsr.erase( v );
        found = true;
        break;
      }
    }
  }

  if( !found ) throw Exceptions( Exceptions::INDEX );

  _issetup    = false;
  _classified = false;
  _on_model_changed( ModelChange::DERIVATIVES );
}

inline double
FFModel::ref
( FFVar const& var )
const
{
  auto same_var = []( FFVar const& a, FFVar const& b ) -> bool
  { return a.id() == b.id(); };

  auto make_default_coord = [&]( t_Dom const& mDom,
                                 std::map<FFVar,double,lt_FFVar> const& classDomRef )
    -> t_Coord
  {
    t_Coord coord;
    for( auto const& [d,dom] : mDom ){
      double val = 0.5 * ( dom.lo_dom + dom.up_dom );
      for( auto const& [dr,x] : classDomRef )
        if( dr.id() == d.id() ){ val = x; break; }
      coord[d] = val;
    }
    return coord;
  };

  auto find_scalar = []( FFVar const& v,
                         std::map<FFVar,double,lt_FFVar> const& ref,
                         double& val ) -> bool
  {
    for( auto const& [vr,x] : ref )
      if( vr.id() == v.id() ){ val = x; return true; }
    return false;
  };

  auto find_fun = []( FFVar const& v,
                      std::map<FFVar,t_Fun,lt_FFVar> const& ref,
                      t_Fun const*& fun ) -> bool
  {
    for( auto const& [vr,f] : ref )
      if( vr.id() == v.id() ){ fun = &f; return true; }
    fun = nullptr;
    return false;
  };

  auto lookup = [&]( t_Var const& mVar, t_Var const& mInp,
                     std::vector<FFVar> const& vCst, std::vector<double> const& vCstVal,
                     t_Dom const& mDom,
                     std::map<FFVar,double,lt_FFVar> const& classStateRef,
                     std::map<FFVar,double,lt_FFVar> const& classInpRef,
                     std::map<FFVar,double,lt_FFVar> const& classDomRef,
                     std::map<FFVar,t_Fun,lt_FFVar> const& classStateRefFun,
                     std::map<FFVar,t_Fun,lt_FFVar> const& classInpRefFun,
                     double& val ) -> bool
  {
    for( auto const& [v,dom] : mVar ){
      (void)dom;
      if( same_var( v, var ) ){
        t_Fun const* fun = nullptr;
        if( find_fun( v, classStateRefFun, fun ) ){
          auto coord = make_default_coord( mDom, classDomRef );
          val = (*fun)( coord );
        }
        else if( !find_scalar( v, classStateRef, val ) ) val = 0.;
        return true;
      }
    }

    for( auto const& [v,dom] : mInp ){
      (void)dom;
      if( same_var( v, var ) ){
        t_Fun const* fun = nullptr;
        if( find_fun( v, classInpRefFun, fun ) ){
          auto coord = make_default_coord( mDom, classDomRef );
          val = (*fun)( coord );
        }
        else if( !find_scalar( v, classInpRef, val ) ) val = 0.;
        return true;
      }
    }

    for( size_t i = 0; i < vCst.size(); ++i ){
      if( same_var( vCst[i], var ) ){
        if( vCstVal.empty() ){
          val = 0.;
          return true;
        }
        if( vCstVal.size() != vCst.size() )
          throw Exceptions( Exceptions::CSTVAL );
        val = vCstVal[i];
        return true;
      }
    }

    for( auto const& [v,dom] : mDom ){
      if( same_var( v, var ) ){
        auto it = classDomRef.find( v );
        val = it != classDomRef.end()? it->second
                                     : 0.5 * ( dom.lo_dom + dom.up_dom );
        return true;
      }
    }

    return false;
  };

  double val = 0.;
  if( _issetup ){
    if( lookup( _mVar, _mInp, _vCst, _vCstVal, _mDom,
                _classVarRef, _classInpRef, _classDomRef,
                _classVarRefFun, _classInpRefFun, val ) )
      return val;
    if( lookup( _usr._mVarUsr, _usr._mInpUsr, _usr._vCstUsr, _usr._vCstValUsr, _usr._mDomUsr,
                _usr._classVarRefUsr, _usr._classInpRefUsr, _usr._classDomRefUsr,
                _usr._classVarRefFunUsr, _usr._classInpRefFunUsr, val ) )
      return val;
  }
  else{
    if( lookup( _usr._mVarUsr, _usr._mInpUsr, _usr._vCstUsr, _usr._vCstValUsr, _usr._mDomUsr,
                _usr._classVarRefUsr, _usr._classInpRefUsr, _usr._classDomRefUsr,
                _usr._classVarRefFunUsr, _usr._classInpRefFunUsr, val ) )
      return val;
    if( lookup( _mVar, _mInp, _vCst, _vCstVal, _mDom,
                _classVarRef, _classInpRef, _classDomRef,
                _classVarRefFun, _classInpRefFun, val ) )
      return val;
  }

  throw Exceptions( Exceptions::INDEX );
}

inline const char*
FFModel::pde_type_name
( EqnType type )
{
  switch( type ){
    case DIFFERENTIAL_ORDINARY:    return "DIFFERENTIAL_ORDINARY";
    case DIFFERENTIAL_ALGEBRAIC:   return "DIFFERENTIAL_ALGEBRAIC";
    case DIFFERENTIAL_IMPLICIT:    return "DIFFERENTIAL_IMPLICIT";
    case ALGEBRAIC_FIELD:          return "ALGEBRAIC_FIELD";
    case RECTANGULAR_RETIRED_:     return "INVALID(retired DIFFERENTIAL_RECTANGULAR)";
    case ALGEBRAIC_LUMPED:         return "ALGEBRAIC_LUMPED";
    case EVOL_HYPERBOLIC:          return "EVOL_HYPERBOLIC";
    case PARABOLIC:                return "PARABOLIC";
    case ELLIPTIC:                 return "ELLIPTIC";
    case COMPLEX_CHARACTERISTIC:   return "COMPLEX_CHARACTERISTIC";
    case UNDETERMINED:             return "UNDETERMINED";
    case SPATIALLY_CHARACTERISTIC: return "SPATIALLY_CHARACTERISTIC";
    case DESCRIPTOR:               return "DESCRIPTOR";
    case DEGENERATE:               return "DEGENERATE";
    case DEGENERATE_LAST_:         return "INVALID(retired WEAK_HYPERBOLIC)";
    default:                       return "UNKNOWN";
  }
}

// ======================================================================
// OCFESLV::_strip_nonlocal_terms_for_symbol
// ======================================================================
inline std::vector<FFVar>
FFModel::_strip_nonlocal_terms_for_symbol
( std::vector<FFVar> const& eqns )
const
{
  std::vector<FFVar> targ;
  std::vector<FFVar> repl;

  auto sg = _dag->subgraph( eqns );
  std::set< FFVar::pt_idVar > seen;

  for( auto const& op : sg.l_op ){
    // Integrals and point evaluations are both nonlocal (a value over a domain, or at a fixed coordinate, not at
    // the collocation point): neither belongs in the principal symbol.  Differentiating an FFEval copied into the
    // working DAG would also mix DAGs (its point is keyed by the user's variable) and throw.
    if( !op->sameid( typeid(FFIntegral) ) && !op->sameid( typeid(FFEval) ) ) continue;

    for( auto const* pout : op->varout ){
      if( !pout ) continue;
      if( !seen.insert( pout->id() ).second ) continue;

      targ.push_back( *pout );
      repl.push_back( FFVar( 0. ) );
    }
  }

  if( targ.empty() ) return eqns;
  return _dag->substitute( eqns, targ, repl );
}

// ======================================================================
// OCFESLV::_principal_symbol
// ======================================================================
inline FFModel::t_Symbol
FFModel::_principal_symbol()
const
{
  return _principal_symbol( 0 );
}

inline FFModel::t_Symbol
FFModel::_principal_symbol
( int const block_id )
const
{
  t_Symbol sym;

  // 1. Identify equations that participate in this block's principal symbol.
  //    Equation metadata provides the third level of control: ordinary
  //    INITIAL/BOUNDARY trace equations are excluded, LINK/INTERIOR
  //    equations are included by default, and SURFACE equations may form their
  //    own lower-dimensional block.
  std::set<FFVar,lt_FFVar> dom_seen;
  std::map<size_t,EqnRole> role_of_eqn;                // rev283: per participating equation
  std::map<size_t,FFVar>   inlined_eqn;                // rev283: eqn id -> chain-inlined form (symbol only)
  for( auto const& eqn : _mEqn ){
    auto const& eqnvar = eqn.var;
    auto const& eqndom = eqn.dom;
    auto const& opt    = *eqn.opt;
    //EqnOptions const opt = _equation_options( eqnvar, eqndom );
    if( !opt.participate_in_classification || opt.block_id != block_id ) continue;

    sym.vEqn.push_back( eqnvar );
    role_of_eqn[ eqnvar.id().second ] = opt.role;   // rev283

    // Domain variables pinned to LB/UB are normal coordinates of a trace
    // equation and are not free manifold directions for the symbol.
    for( auto const& [v,lim] : eqndom ){
      if( lim == FFDom::LB || lim == FFDom::UB ) continue;
      if( _mDom.find(v) != _mDom.end() ) dom_seen.insert( v );
    }
  }

  // Preserve global map ordering for block-free domains.  For block-level
  // classification, keep only states that actually occur differentiated in
  // this block.  This prevents an unrelated state in another block from
  // making the block symbol artificially rectangular.
  for( auto const& [v, dom]  : _mDom )
    if( dom_seen.find(v) != dom_seen.end() ) sym.vDom.push_back( v );

  std::set<FFVar,lt_FFVar> state_seen;
  std::set<FFVar,lt_FFVar> all_block_states;        // every state that appears (differentiated or not)
  std::vector<FFVar> diff_eqn;                      // participating equations with >=1 differentiated state
  std::map<size_t, std::set<FFVar,lt_FFVar>> eqn_diff;  // eqn id -> states it differentiates
  std::map<size_t, int> diff_count;                     // state id -> # differential equations differentiating it
  for( auto const& eqnvar : sym.vEqn ){
    bool eqn_has_deriv = false;
    std::set<FFVar,lt_FFVar> this_diff;
    auto sg = _dag->subgraph( 1, &eqnvar );
    for( auto const& op : sg.l_op ){
      if( op->type == FFOp::VAR && !op->varout.empty() && op->varout[0]
          && _mVar.find( *op->varout[0] ) != _mVar.end() )
        all_block_states.insert( *op->varout[0] );
      if( !op->sameid( typeid(FFPartial) ) ) continue;
      auto const* pop = mc::type_cast<FFPartial const>( op );
      for( size_t jj = 0; jj < op->varin.size(); ++jj ){
        FFVar const* operand = op->varin[jj];
        if( _mVar.find( *operand ) == _mVar.end() ) continue;
        bool in_free_domain = false;
        for( auto const& [indep_var, ord] : pop->Indep().expr ){
          for( auto const& d : sym.vDom )
            if( d.id() == indep_var.id() ){ in_free_domain = true; break; }
          if( in_free_domain ) break;
        }
        if( in_free_domain ){ state_seen.insert( *operand ); this_diff.insert( *operand ); eqn_has_deriv = true; }
      }
    }
    // rev283 (T-P1): a balance row that reduce_order left ALGEBRAIC in a bare auxiliary is still a balance row.
    // Inline each such auxiliary one level (aux -> its defining derivative) and re-scan; the inlined form is used
    // for the symbol only.  INTERIOR rows only: LINK rows ARE the definitions.
    // 4.0b (2026-10-08): ALSO a differential row, for an auxiliary that appears in it only BARE (never under a
    // derivative in the row): u_t + Dz_u = 0 -- the advection PDE once a second-derivative BOUNDARY row has made
    // RED_FULL substitute Dz_u -- is read u_t + u_z = 0.  Without this its symbol column u_z was missing and the
    // block came out rectangular (UNDETERMINED).  An auxiliary differentiated in the row (u_t - D w_z, the reduced
    // diffusion form) is left alone: there the first-order system in (u, w) is the intended reading.
    bool const had_deriv = eqn_has_deriv;
    if( _knob_CHAIN_SYMBOL() ){
      auto itr = role_of_eqn.find( eqnvar.id().second );
      if( itr != role_of_eqn.end() && itr->second == EqnRole::INTERIOR ){
        std::vector<FFVar> targ, repl;
        for( auto const& aux : _auxDef ){
          if( all_block_states.find( aux.aux ) == all_block_states.end() ) continue;
          bool in_row = false;
          for( auto const& op : sg.l_op )
            if( op->type == FFOp::VAR && !op->varout.empty() && op->varout[0] && op->varout[0]->id() == aux.aux.id() ){ in_row = true; break; }
          if( !in_row ) continue;
          if( had_deriv ){                                     // a differential row: bare occurrences only
            bool aux_differentiated = false;
            for( auto const& op : sg.l_op ){
              if( !op->sameid( typeid(FFPartial) ) ) continue;
              for( auto const* vin : op->varin ) if( vin && vin->id() == aux.aux.id() ){ aux_differentiated = true; break; }
              if( aux_differentiated ) break;
            }
            if( aux_differentiated ) continue;
          }
          bool expr_has_deriv = false;
          auto sge = _dag->subgraph( 1, &aux.expr );
          for( auto const& ope : sge.l_op ) if( ope->sameid( typeid(FFPartial) ) ){ expr_has_deriv = true; break; }
          if( !expr_has_deriv ) continue;
          targ.push_back( aux.aux ); repl.push_back( aux.expr );
        }
        if( !targ.empty() ){
          FFVar const inl = _dag->substitute( std::vector<FFVar>{ eqnvar }, targ, repl )[0];
          auto sgi = _dag->subgraph( 1, &inl );
          for( auto const& op : sgi.l_op ){
            if( !op->sameid( typeid(FFPartial) ) ) continue;
            auto const* pop = mc::type_cast<FFPartial const>( op );
            for( size_t jj = 0; jj < op->varin.size(); ++jj ){
              FFVar const* operand = op->varin[jj];
              if( _mVar.find( *operand ) == _mVar.end() ) continue;
              bool in_free_domain = false;
              for( auto const& [indep_var, ord] : pop->Indep().expr ){
                for( auto const& d : sym.vDom ) if( d.id() == indep_var.id() ){ in_free_domain = true; break; }
                if( in_free_domain ) break; }
              if( in_free_domain ){ state_seen.insert( *operand ); this_diff.insert( *operand ); eqn_has_deriv = true; }
            }
          }
          if( eqn_has_deriv ){
            inlined_eqn[ eqnvar.id().second ] = inl;
            if( !had_deriv && options.DISPLAY_LEVEL >= 1 )
              std::cerr << "OCFESLV::_principal_symbol ** chain-aware: row " << eqnvar << " is algebraic in "
                        << targ.size() << " reduction auxiliar" << ( targ.size() == 1 ? "y" : "ies" )
                        << " -- read with the definition(s) inlined; it is a balance row" << std::endl;
          }
        }
      }
    }
    if( eqn_has_deriv ){
      diff_eqn.push_back( eqnvar );
      eqn_diff[ eqnvar.id().second ] = this_diff;
      for( auto const& sdiff : this_diff ) ++diff_count[ sdiff.id().second ];
    }
    else sym.vAlgEqn.push_back( eqnvar );
  }

  // Value-slaved / derived-chain refinement.  A PURE value-slaved state -- one that appears in the
  // block but is NEVER differentiated (so it has no principal-symbol column) -- is defined by an
  // ALGEBRAIC equation even when that equation carries a derivative of ANOTHER state.  Example: a
  // distributed OpP-chain field  w = d/dz( u * du/dz ).  Under whole-operand materialisation this is a
  // DERIVED CHAIN  u -> aux(=u*du/dz, LINK-defined) -> w(=d(aux)/dz):
  //     w  - d(aux)/dz          = 0   (defines w,   carries d(aux)/dz)
  //     u*du/dz - aux           = 0   (defines aux, carries du/dz)
  // Both defining equations belong in the ALGEBRAIC part; left in the differential symbol they are
  // spurious extra rows (rectangular / rank-deficient).  Route a defining equation to vAlgEqn iff
  // (a) it contains a value-slaved state, and (b) EVERY state it differentiates is COVERED --
  //   (b1) differentiated by another differential equation (diff_count > 1, so no symbol column loses
  //        its provider), OR
  //   (b2) a materialised auxiliary (in _auxVarID): it appears undifferentiated as the output of its
  //        own defining LINK, hence is fully determined and its derivative is a derivative of a known
  //        quantity -- exactly the LINK-defined coverability the planned Pantelides reduction formalises.
  // Routing one derived state's definition exposes the next in the chain (routing w's equation leaves
  // aux undifferentiated -> value-slaved -> route aux's LINK), so ITERATE to a fixpoint, recomputing
  // state_seen / diff_count from the shrinking differential set each pass.  Condition (b1)/(a) still
  // prevent misrouting a genuine PDE that merely USES a value-slaved coefficient: its own field's
  // derivative is neither supplied elsewhere nor an aux, so (b) fails -- an EOS/wetting-law coefficient
  // in a transport equation is untouched (and PDE22's flux divergence has NO value-slaved state, so
  // (a) fails and it never enters this path).
  {
    // 4.0b (2026-10-08): (b1) is judged per (state, DIRECTION), not per state.  A row is the provider of the
    // derivatives it takes; "covered elsewhere" must mean the SAME derivative is supplied by another row.  Counted
    // per state, a LINK Dz_u - u_z (a z-derivative) "covered" the PDE u_t + Dz_u (the t-derivative) once a
    // second-derivative BOUNDARY row had made reduction substitute Dz_u into it: the PDE was routed to vAlgEqn,
    // then u, then everything -- an EMPTY symbol, classified ALGEBRAIC_FIELD (measured, OCFE: u_zz = 0 at an outflow).
    bool routed_any = true;
    while( routed_any ){
      routed_any = false;
      std::set<FFVar,lt_FFVar> value_slaved;
      for( auto const& sv : all_block_states )
        if( state_seen.find( sv ) == state_seen.end() ) value_slaved.insert( sv );
      {
        static bool const kSymRow = Options::_env_flag( "CRONOS_AUDIT_SYMROW", false );
        if( kSymRow ){
          std::cerr << "  [symrow] pass: value-slaved states (appear under NO derivative"
                       " anywhere, so they hold no symbol column) = " << value_slaved.size();
          for( auto const& v : value_slaved ) std::cerr << " " << v.name();
          std::cerr << "\n";
        }
      }
      if( value_slaved.empty() ) break;
      std::map<size_t, std::set<std::pair<size_t,size_t>>> eqn_diff_dir;   // eqn id -> {(state id, direction id)}
      std::map<std::pair<size_t,size_t>, int> diff_count_dir;              // (state, direction) -> # rows taking it
      for( auto const& eqnvar : diff_eqn ){
        auto const it_i = inlined_eqn.find( eqnvar.id().second );            // rev283: the inlined row, as below
        auto sgd = _dag->subgraph( 1, it_i != inlined_eqn.end() ? &it_i->second : &eqnvar );
        auto& pr = eqn_diff_dir[ eqnvar.id().second ];
        for( auto const& op : sgd.l_op ){
          if( !op->sameid( typeid(FFPartial) ) ) continue;
          auto const* pop = mc::type_cast<FFPartial const>( op );
          for( size_t jj = 0; jj < op->varin.size(); ++jj ){
            FFVar const* operand = op->varin[jj];
            if( _mVar.find( *operand ) == _mVar.end() ) continue;
            for( auto const& [indep_var, ord] : pop->Indep().expr )
              for( auto const& d : sym.vDom )
                if( d.id() == indep_var.id() ) pr.insert( { operand->id().second, indep_var.id().second } );
          }
        }
        for( auto const& pd : pr ) ++diff_count_dir[ pd ];
      }
      // (b1) every direction in which THIS row differentiates sd is also taken by another row; or (b2)
      auto covered = [&]( FFVar const& sd, size_t eqn_id )->bool {
        if( _auxVarID.find( sd.id() ) != _auxVarID.end() ) return true;        // (b2) LINK-defined materialised aux
        for( auto const& pd : eqn_diff_dir[ eqn_id ] )
          if( pd.first == sd.id().second && diff_count_dir[ pd ] <= 1 ) return false;
        return true;
      };
      std::vector<FFVar> keep_diff;
      bool any_moved = false;
      for( auto const& eqnvar : diff_eqn ){
        bool has_vs = false;
        auto const it_inl = inlined_eqn.find( eqnvar.id().second );   // rev283: judge the inlined row
        FFVar const& row_for_scan = ( it_inl != inlined_eqn.end() ) ? it_inl->second : eqnvar;
        auto sg = _dag->subgraph( 1, &row_for_scan );
        for( auto const& op : sg.l_op )
          if( op->type == FFOp::VAR && !op->varout.empty() && op->varout[0]
              && value_slaved.find( *op->varout[0] ) != value_slaved.end() ){ has_vs = true; break; }
        bool all_covered = has_vs;
        if( has_vs )
          for( auto const& sdiff : eqn_diff[ eqnvar.id().second ] )
            if( !covered( sdiff, eqnvar.id().second ) ){ all_covered = false; break; }
        // rev69 INSTRUMENT [symrow].  This routing rule exists, by its own comment, to
        // "avoid spurious extra rows (rectangular / rank-deficient)" -- and on MMPDE27 it
        // PRODUCES a rectangular symbol: sym=3x4, hence type=UNDETERMINED with
        // evol_dom_idx=none and no evolution matrix formed at all.
        //
        // The asymmetry to expose: ROWS leave here when every state they differentiate is
        // covered elsewhere, but COLUMNS are fixed earlier by "appears under any OpP" and
        // never shrink to match.  If a state's every differentiating row is routed away,
        // its column survives with no provider -- three rows removed, zero columns removed.
        //
        // Printed unconditionally under CRONOS_AUDIT_SYMROW so the decision can be read
        // per equation instead of inferred from the final dimensions.  Inference already
        // failed once here: vEqn.push_back is UNCONDITIONAL for every classify=y row, so
        // the reduction happens ONLY in this loop.
        static bool const kSymRow = Options::_env_flag( "CRONOS_AUDIT_SYMROW", false );
        if( kSymRow ){
          std::cerr << "  [symrow] eqn " << eqnvar.name()
                    << "  has_value_slaved=" << ( has_vs ? "y" : "n" )
                    << "  all_diffd_states_covered=" << ( all_covered ? "y" : "n" )
                    << "  -> " << ( has_vs && all_covered ? "ROUTED to vAlgEqn"
                                                          : "kept in vEqn" ) << "\n";
          if( has_vs )
            for( auto const& sdiff : eqn_diff[ eqnvar.id().second ] ){
              auto it = diff_count.find( sdiff.id().second );
              std::cerr << "  [symrow]     differentiates " << sdiff.name()
                        << "  diff_count=" << ( it != diff_count.end() ? it->second : 0 )
                        << "  is_aux=" << ( _auxVarID.find( sdiff.id() ) != _auxVarID.end()
                                              ? "y" : "n" )
                        << "  covered=" << ( covered( sdiff, eqnvar.id().second ) ? "y" : "n" ) << "\n";
            }
        }
        if( has_vs && all_covered ){ sym.vAlgEqn.push_back( eqnvar ); any_moved = true; }
        else                        keep_diff.push_back( eqnvar );
      }
      if( !any_moved ) break;
      diff_eqn.swap( keep_diff );
      routed_any = true;
      // Recompute the differentiated-state bookkeeping from the reduced differential set, so the next
      // pass sees states that just became undifferentiated (the next link in the derived chain).
      state_seen.clear();
      diff_count.clear();
      eqn_diff.clear();
      for( auto const& eqnvar : diff_eqn ){
        std::set<FFVar,lt_FFVar> this_diff;
        auto const it_inl2 = inlined_eqn.find( eqnvar.id().second );   // rev283: rescan the inlined row
        auto sg = _dag->subgraph( 1, it_inl2 != inlined_eqn.end() ? &it_inl2->second : &eqnvar );
        for( auto const& op : sg.l_op ){
          if( !op->sameid( typeid(FFPartial) ) ) continue;
          auto const* pop = mc::type_cast<FFPartial const>( op );
          for( size_t jj = 0; jj < op->varin.size(); ++jj ){
            FFVar const* operand = op->varin[jj];
            if( _mVar.find( *operand ) == _mVar.end() ) continue;
            bool in_free_domain = false;
            for( auto const& [indep_var, ord] : pop->Indep().expr ){
              for( auto const& d : sym.vDom )
                if( d.id() == indep_var.id() ){ in_free_domain = true; break; }
              if( in_free_domain ) break;
            }
            if( in_free_domain ){ state_seen.insert( *operand ); this_diff.insert( *operand ); }
          }
        }
        eqn_diff[ eqnvar.id().second ] = this_diff;
        for( auto const& sdiff : this_diff ) ++diff_count[ sdiff.id().second ];
      }
    }
  }
  // A participating equation with NO differentiated state (an algebraic closure such as an
  // equilibrium or a wetting law omega = g(...)) contributes only a structurally-zero row to the
  // DIFFERENTIAL principal symbol, making it rectangular/rank-deficient and corrupting the block
  // classification.  Algebraic constraints are the province of the index analysis / DAE reduction,
  // not the differential symbol, so the symbol rows are the differential equations only; the
  // algebraic ones are recorded in sym.vAlgEqn (their bare states, e.g. omega, form the algebraic
  // block that the index analysis pins via its own INTERIOR closure).
  {
    // rev126: the FIXPOINT line MOVED.  It used to sit here and read sym.vState.size()
    // TWELVE LINES BEFORE sym.vState is filled, so it printed "columns (vState)=0"
    // unconditionally and its SQUARE/RECTANGULAR verdict compared rows against zero --
    // i.e. it announced *** RECTANGULAR *** on every model, including ones the classify
    // line then reports as sym=2x2.  Measured on OCFE_PDE6_blk0: "rows kept=2 ...
    // columns (vState)=0 *** RECTANGULAR ***" on BOTH formulations, while [classify]
    // says PARABOLIC sym=2x2 -- square.  The warning was pure artefact of placement.
    //
    // See the relocated line below, after sym.vState is populated.
  }
  sym.vEqn = diff_eqn;
  for( auto const& [v, sdom] : _mVar )
    if( state_seen.find(v) != state_seen.end() ) sym.vState.push_back( v );

  {
    // rev126: the FIXPOINT report, now AFTER sym.vState exists so its numbers are real.
    static bool const kSymRow = Options::_env_flag( "CRONOS_AUDIT_SYMROW", false );
    if( kSymRow )
      std::cerr << "  [symrow] FIXPOINT: rows kept=" << diff_eqn.size()
                << " routed to vAlgEqn=" << sym.vAlgEqn.size()
                << "  columns (vState)=" << sym.vState.size()
                << ( diff_eqn.size() == sym.vState.size()
                       ? "   SQUARE"
                       : "   *** RECTANGULAR: rows shrank, columns did not.  A column whose"
                         " every differentiating row was routed away has no provider left,"
                         " and classification returns UNDETERMINED ***" )
                << "\n";
  }

  size_t const nDom   = sym.vDom.size();
  size_t const nState = sym.vState.size();
  size_t const nEqn   = sym.vEqn.size();
  if( !nEqn || !nState || !nDom ) return sym;
  // rev283: the rows as the SYMBOL reads them (chain-inlined where flagged); sym.vEqn keeps the original vars.
  std::vector<FFVar> vEqnSym = sym.vEqn;
  for( size_t k = 0; k < vEqnSym.size(); ++k ){
    auto it = inlined_eqn.find( vEqnSym[k].id().second );
    if( it != inlined_eqn.end() ) vEqnSym[k] = it->second;
  }

  // 3. Locate FFPartial AUX output nodes across interior equations
  // deriv_node[i][j] = DAG AUX node for du_j/dx_i, nullptr if absent
  std::vector< std::vector<FFVar const*> > deriv_node(
      nDom, std::vector<FFVar const*>( nState, nullptr ) );

  for( auto const& eqnvar : vEqnSym ){   // rev283: inlined rows carry their derivative node here
    auto sg = _dag->subgraph( 1, &eqnvar );
    for( auto const& op : sg.l_op ){
      if( !op->sameid( typeid(FFPartial) ) ) continue;
      auto const* pop = mc::type_cast<FFPartial const>( op );
      for( size_t jj = 0; jj < op->varin.size(); ++jj ){
        FFVar const* operand = op->varin[jj];
        FFVar const* dag_out = op->varout[jj];
        size_t j = nState;
        for( size_t jj2 = 0; jj2 < nState; ++jj2 )
          if( sym.vState[jj2].id() == operand->id() ){ j = jj2; break; }
        if( j == nState ) continue;
        for( auto const& [indep_var, ord] : pop->Indep().expr ){
          size_t i = nDom;
          for( size_t ii = 0; ii < nDom; ++ii )
            if( sym.vDom[ii].id() == indep_var.id() ){ i = ii; break; }
          if( i == nDom ) continue;
          deriv_node[i][j] = dag_out;
        }
      }
    }
  }

  // 4. Create proxy VAR nodes for each (domain, state) derivative slot
  std::vector<FFVar> sub_targ, sub_repl;
  std::vector< std::vector<FFVar> > proxy( nDom, std::vector<FFVar>( nState ) );

  for( size_t i = 0; i < nDom; ++i )
    for( size_t j = 0; j < nState; ++j ){
      proxy[i][j] = _dag->add_var( "v_" + sym.vState[j].name()
                                        + "_" + sym.vDom[i].name() );
      if( deriv_node[i][j] ){
        sub_targ.push_back( *deriv_node[i][j] );
        sub_repl.push_back(  proxy[i][j] );
      }
    }

  // 5. Substitute derivative AUX nodes -> proxy VARs in interior equations.
  auto sub_eqns = _dag->substitute( vEqnSym, sub_targ, sub_repl );   // rev283

  // Flat proxy order: proxy[0][0], proxy[0][1], ..., proxy[1][0], ...
  // J[ k*(nDom*nState) + i*nState + j ] = dF_k/d(du_j/dx_i)
  std::vector<FFVar> all_proxies;
  all_proxies.reserve( nDom * nState );
  for( size_t i = 0; i < nDom; ++i )
    for( size_t j = 0; j < nState; ++j )
      all_proxies.push_back( proxy[i][j] );

  // Retain the proxy VARs so vCoeff0 (which may reference them when a derivative
  // coefficient is state-dependent) can be evaluated later (item 11 Part-2).
  sym.vDerivProxy = all_proxies;

  // 6. Nonlocal terms such as FFIntegral are lower order for the ordinary
  // local principal symbol, so they should not be passed directly to DAG FAD
  // because FFIntegral::deriv is intentionally undefined.  Before stripping
  // them, however, check whether a nonlocal term appears in a coefficient of
  // a derivative proxy, e.g. (Integral(u))*u_x.  Such a problem has a
  // nonlocal principal coefficient and cannot be classified by this local
  // symbol extractor.
  {
    std::vector<FFVar> nl_targ, nl_repl;
    std::set< FFVar::pt_idVar > nl_seen;

    auto sg = _dag->subgraph( sub_eqns );
    for( auto const& op : sg.l_op ){
      if( !op->sameid( typeid(FFIntegral) ) ) continue;
      for( auto const* pout : op->varout ){
        if( !pout ) continue;
        if( !nl_seen.insert( pout->id() ).second ) continue;
        nl_targ.push_back( *pout );
        nl_repl.push_back( _dag->add_var( "nonlocal_symbol_"
                         + std::to_string( nl_repl.size() ) ) );
      }
    }

    if( !nl_targ.empty() ){
      auto probe_eqns = _dag->substitute( sub_eqns, nl_targ, nl_repl );
      auto probe_jac  = _dag->FAD( probe_eqns, all_proxies );

      std::set< FFVar::pt_idVar > nl_var_id;
      for( auto const& z : nl_repl ) nl_var_id.insert( z.id() );

      bool nonlocal_principal_coeff = false;
      for( auto const& coeff : probe_jac ){
        auto csg = _dag->subgraph( 1, &coeff );
        for( auto const& op : csg.l_op ){
          if( op->type != FFOp::VAR || op->varout.empty() || !op->varout[0] )
            continue;
          if( nl_var_id.find( op->varout[0]->id() ) != nl_var_id.end() ){
            nonlocal_principal_coeff = true;
            break;
          }
        }
        if( nonlocal_principal_coeff ) break;
      }

      if( nonlocal_principal_coeff ){
        std::cerr << "OCFESLV::_principal_symbol ** block " << block_id
                  << ": nonlocal FFIntegral term appears in a principal "
                  << "coefficient; PDE classification set to "
                  << pde_type_name( UNDETERMINED ) << std::endl;
        return t_Symbol();
      }

      sub_eqns = _strip_nonlocal_terms_for_symbol( sub_eqns );
    }
  }

  // 7. Differentiate w.r.t. all proxy VARs to obtain the Jacobian.
  auto jac = _dag->FAD( sub_eqns, all_proxies );

  // 7. Distribute into per-domain coefficient matrices
  sym.vCoeff.assign( nDom, std::vector<FFVar>( nEqn * nState ) );
  for( size_t k = 0; k < nEqn; ++k )
    for( size_t i = 0; i < nDom; ++i )
      for( size_t j = 0; j < nState; ++j )
        sym.vCoeff[i][ k * nState + j ] = jac[ k * (nDom * nState) + i * nState + j ];

#ifdef MC__OCFESLV_SYMBOL_BALANCE_PROBE
  // Item 11 Part-2: zeroth-order state Jacobian from the SAME proxied system.
#endif
  // In sub_eqns every derivative term is an independent proxy VAR, so FAD w.r.t.
  // the state VALUES captures only the algebraic/lower-order coupling (e.g. the
  // advective Dz_Cl*(...) term), which the principal symbol omits.  Used by
  // _compute_interface_reads() for the symbol-AUTO magnitude balance.
  {
    auto jac0 = _dag->FAD( sub_eqns, sym.vState );
    sym.vCoeff0.assign( nEqn * nState, FFVar( 0. ) );
    for( size_t k = 0; k < nEqn; ++k )
      for( size_t j = 0; j < nState; ++j )
        sym.vCoeff0[ k * nState + j ] = jac0[ k * nState + j ];
  }

  // rev179 S0'': the DIFFERENTIAL rows' zeroth-order Jacobian over ALL states (blk0's rescued
  // edges receive on differential rows where the auxiliary is undifferentiated and outside
  // vState -- vCoeff0 cannot see them).  Same sub_eqns, wider column set, guarded.
  try{
    std::vector<FFVar> all_states;
    for( auto const& vv : _mVar ) all_states.push_back( vv.first );
    auto jacAll = _dag->FAD( sub_eqns, all_states );
    size_t const nS = all_states.size();
    sym.vAllState = all_states;
    sym.vCoeff0All.assign( nEqn * nS, FFVar( 0. ) );
    for( size_t k = 0; k < nEqn; ++k )
      for( size_t j = 0; j < nS; ++j )
        sym.vCoeff0All[ k * nS + j ] = jacAll[ k * nS + j ];
  }
  catch( ... ){ sym.vAllState.clear(); sym.vCoeff0All.clear();
    if( options.DISPLAY_LEVEL >= 1 ) std::cerr << "OCFESLV::_principal_symbol ** vCoeff0All not computed" << std::endl; }

  // rev296: C2's own table.  Built from the rows AS WRITTEN (only FFPartial outputs proxied, no stripping) so a
  // bare auxiliary is differentiated as itself; see the header doc for the measurement that motivated it.
  try{
    std::vector<FFVar> c2_rows( sym.vEqn );
    c2_rows.insert( c2_rows.end(), sym.vAlgEqn.begin(), sym.vAlgEqn.end() );
    if( !c2_rows.empty() ){
      std::vector<FFVar> c2_targ, c2_repl;
      { std::set<size_t> seen; size_t np = 0;
        auto sg = _dag->subgraph( c2_rows.size(), c2_rows.data() );
        for( auto const& op : sg.l_op ){
          // rev300: EVERY external operation, not just FFPartial.  FAD cannot differentiate through an
          // external (FFGraph::Exceptions ierr=-2), and a row carrying an FFIntegral or FFEval made the whole
          // table fail to build -- silently, before rev299 named it.  An external's OUTPUT is never a state,
          // so proxying it cannot hide the bare state this table exists to see.
          if( op->type != FFOp::EXTERN || op->varout.empty() || !op->varout[0] ) continue;
          FFVar const& out = *op->varout[0];
          if( !seen.insert( (size_t)out.id().second ).second ) continue;
          c2_targ.push_back( out );
          c2_repl.push_back( _dag->add_var( "c2coef_proxy_" + std::to_string( np++ ) ) ); } }
      auto c2_sub = c2_targ.empty() ? c2_rows : _dag->substitute( c2_rows, c2_targ, c2_repl );
      std::vector<FFVar> all_states;
      for( auto const& vv : _mVar ) all_states.push_back( vv.first );
      auto jacC2 = _dag->FAD( c2_sub, all_states );
      size_t const nR = c2_rows.size(), nS2 = all_states.size();
      sym.vC2Eqn = c2_rows; sym.vC2State = all_states;
      sym.vC2Coef0.assign( nR * nS2, FFVar( 0. ) );
      for( size_t k = 0; k < nR; ++k )
        for( size_t j = 0; j < nS2; ++j )
          sym.vC2Coef0[ k * nS2 + j ] = jacC2[ k * nS2 + j ];
    }
  }
  // rev299 (D2): say WHY.  A bare catch left 18 lines in 10 programs unexplained, and their models fall back to
  // the two symbol tables -- so the rev296 fix is unvalidated exactly there.
  catch( FFGraph::Exceptions& ex ){   // rev299: MC++'s own type -- the bare catch hid it
    sym.vC2Eqn.clear(); sym.vC2State.clear(); sym.vC2Coef0.clear();
    std::cerr << "OCFESLV::_principal_symbol ** vC2Coef0 not computed: FFGraph::Exceptions ierr=" << ex.ierr()
              << " (" << ex.what() << ")  [rows=" << ( sym.vEqn.size() + sym.vAlgEqn.size() )
              << " states=" << _mVar.size() << "]" << std::endl; }
  catch( std::exception const& ex ){
    sym.vC2Eqn.clear(); sym.vC2State.clear(); sym.vC2Coef0.clear();
    std::cerr << "OCFESLV::_principal_symbol ** vC2Coef0 not computed: " << ex.what()
              << "  [rows=" << ( sym.vEqn.size() + sym.vAlgEqn.size() ) << " states=" << _mVar.size() << "]"
              << std::endl; }
  catch( ... ){
    sym.vC2Eqn.clear(); sym.vC2State.clear(); sym.vC2Coef0.clear();
    std::cerr << "OCFESLV::_principal_symbol ** vC2Coef0 not computed: non-std exception"
              << "  [rows=" << ( sym.vEqn.size() + sym.vAlgEqn.size() ) << " states=" << _mVar.size() << "]"
              << std::endl; }

  // rev177 S0': the algebraic rows' zeroth-order Jacobian, over ALL states.  WHY: a claim on
  // an auxiliary that the symbol excludes (it appears only undifferentiated) receives on an
  // algebraic row, gets a ZERO natural coupling, and is rescued with a FABRICATED +1 -- which
  // on OCFE_PDE20f model C cancels the row's own unit coefficient of that auxiliary exactly
  // (1.69e+00 at +1.0, ~3.9e-06 at every other value).  The coefficient the symbol WOULD have
  // supplied is this one.  Same proxies, same stripping, same discipline as vCoeff0.  Guarded:
  // any failure leaves the arrays EMPTY and changes nothing else.
  if( !sym.vAlgEqn.empty() ){
    try{
      std::vector<FFVar> all_states;
      for( auto const& vv : _mVar ) all_states.push_back( vv.first );
      // The differential rows' proxies (sub_targ/sub_repl) cover only the derivative slots
      // FOUND IN vEqn.  An algebraic row can carry an FFPartial node with no such proxy --
      // under sigma reuse the PDE has no du/dz node at all, yet LINK and ALG both do -- and a
      // raw FFPartial under FAD throws.  So: collect every FFPartial OUTPUT in the algebraic
      // rows and give each its own proxy VAR before differentiating.
      std::vector<FFVar> alg_targ( sub_targ ), alg_repl( sub_repl );
      {
        std::set<size_t> seen;
        for( auto const& t : sub_targ ) seen.insert( (size_t)t.id().second );
        auto sg = _dag->subgraph( sym.vAlgEqn.size(), sym.vAlgEqn.data() );
        size_t np = 0;
        for( auto const& op : sg.l_op ){
          if( !op->sameid( typeid(FFPartial) ) || op->varout.empty() || !op->varout[0] ) continue;
          FFVar const& out = *op->varout[0];
          if( seen.count( (size_t)out.id().second ) ) continue;
          seen.insert( (size_t)out.id().second );
          alg_targ.push_back( out );
          alg_repl.push_back( _dag->add_var( "algcoef_proxy_" + std::to_string( np++ ) ) );
        }
      }
      auto alg_sub = _dag->substitute( sym.vAlgEqn, alg_targ, alg_repl );
      alg_sub = _strip_nonlocal_terms_for_symbol( alg_sub );
      auto jacA = _dag->FAD( alg_sub, all_states );
      size_t const nA = sym.vAlgEqn.size(), nS = all_states.size();
      sym.vAlgState = all_states;
      sym.vAlgCoeff0.assign( nA * nS, FFVar( 0. ) );
      for( size_t g = 0; g < nA; ++g )
        for( size_t j = 0; j < nS; ++j )
          sym.vAlgCoeff0[ g * nS + j ] = jacA[ g * nS + j ];
    }
    catch( std::exception const& e ){
      sym.vAlgState.clear(); sym.vAlgCoeff0.clear();
      if( options.DISPLAY_LEVEL >= 1 )
        std::cerr << "OCFESLV::_principal_symbol ** vAlgCoeff0 not computed: " << e.what() << std::endl;
    }
    catch( ... ){
      sym.vAlgState.clear(); sym.vAlgCoeff0.clear();
      if( options.DISPLAY_LEVEL >= 1 )
        std::cerr << "OCFESLV::_principal_symbol ** vAlgCoeff0 not computed (non-std exception)" << std::endl;
    }
  }

  return sym;
}

// ======================================================================
// OCFESLV::_eval_symbol
// ======================================================================
inline std::vector<arma::mat>
FFModel::_eval_symbol
( t_Symbol             const& sym,
  std::vector<double>  const& state_vals,
  std::vector<double>  const& cst_vals,
  std::vector<double>  const& dom_vals,
  std::vector<double>  const& input_vals,
  arma::mat*                  A0 )
const
{
  size_t const nDom   = sym.vDom.size();
  size_t const nState = sym.vState.size();
  size_t const nEqn   = sym.vEqn.size();

  std::vector<FFVar>  eval_vars;
  std::vector<double> eval_vals;
  eval_vars.reserve( nState + nDom + _vInp.size() + _vCst.size() );
  eval_vals.reserve( nState + nDom + _vInp.size() + _vCst.size() );

  // The block symbol may contain a subset of the global state/domain
  // variables.  Map reference values by FFVar identity, not by the local
  // block ordering, so multi-block/SURFACE classifications are evaluated at
  // the intended reference point.  The positional fallback preserves the
  // legacy behaviour for callers that pass block-local vectors.
  auto value_by_id = []( FFVar const& v, auto const& vars,
                         std::vector<double> const& vals, size_t fallback )
    -> double
  {
    for( size_t k = 0; k < vars.size() && k < vals.size(); ++k )
      if( vars[k].id() == v.id() ) return vals[k];
    return fallback < vals.size() ? vals[fallback] : 0.;
  };

  for( size_t j = 0; j < nState; ++j ){
    eval_vars.push_back( sym.vState[j] );
    eval_vals.push_back( value_by_id( sym.vState[j], _vVar, state_vals, j ) );
  }
  // Coefficients (sym.vCoeff) may reference NON-principal states absent from sym.vState -- an
  // algebraic coefficient-only state such as membrane wetting omega(z), which is never
  // differentiated and so is correctly excluded from the principal symbol, yet appears in the
  // interface/reaction coefficients.  Provision every remaining global state (value keyed to
  // _vVar via value_by_id) so the coefficient eval resolves it; the symbol math (FAD w.r.t.
  // sym.vState) is unaffected.
  for( size_t p = 0; p < _vVar.size(); ++p ){
    bool principal = false;
    for( size_t j = 0; j < nState; ++j )
      if( sym.vState[j].id() == _vVar[p].id() ){ principal = true; break; }
    if( principal ) continue;
    eval_vars.push_back( _vVar[p] );
    eval_vals.push_back( value_by_id( _vVar[p], _vVar, state_vals, p ) );
  }
  for( size_t i = 0; i < nDom; ++i ){
    eval_vars.push_back( sym.vDom[i] );
    eval_vals.push_back( value_by_id( sym.vDom[i], _vDom, dom_vals, i ) );
  }
  // Principal coefficients may depend on collocated inputs.  The symbol is
  // evaluated in FFVar space, so a distributed input is represented by one
  // scalar reference value and is therefore classified as spatially constant.
  for( size_t p = 0; p < _vInp.size(); ++p ){
    eval_vars.push_back( _vInp[p] );
    eval_vals.push_back( value_by_id( _vInp[p], _vInp, input_vals, p ) );
  }
  for( size_t c = 0; c < _vCst.size(); ++c ){
    eval_vars.push_back( _vCst[c] );
    eval_vals.push_back( c < cst_vals.size() ? cst_vals[c] : 0. );
  }

  // GUARD (rev67): nDom is sym.vDom.size(), but the loop below indexes sym.vCoeff.
  // Those are DIFFERENT vectors and nothing enforced that they had the same length.
  // A DISTRIBUTED PURELY ALGEBRAIC field -- states on a collocation mesh with NO
  // derivative of any state in any equation -- populates vDom from the environment's
  // domains while vCoeff stays EMPTY, and sym.vCoeff[0] then walked off the end,
  // aborting inside setup() with _GLIBCXX_ASSERTIONS on (rev66:14391).  Measured on
  // OCFE_ALGFIELD1.  The corpus never hit it because OCFE_AE0/AE1/AE2 are ZERO-domain,
  // so this loop does not execute for them.
  //
  // _eval_symbol returns the per-domain coefficient matrices and has no classification
  // result to fill, so it returns an EMPTY vector here; the CALLER (_classify_pde)
  // distinguishes "no principal symbol" from a populated one and sets ALGEBRAIC_FIELD /
  // ALGEBRAIC_LUMPED accordingly.  Returning empty rather than nDom zero matrices is
  // deliberate: a zero matrix is indistinguishable from a genuinely zero coefficient.
  if( A0 && sym.vCoeff0.size() == nEqn * nState ){          // the zeroth-order block, on request
    // sym.vCoeff0 is the state Jacobian of the PROXIED equations, so it may contain the derivative proxies
    // (v_Cg_z for dCg/dz under an advective term): at the reference point -- a constant field -- the derivatives
    // vanish, so the proxies are supplied at 0 (2026-10-05)
    std::vector<FFVar>  v0 = eval_vars;  std::vector<double> x0 = eval_vals;
    for( auto const& pr : sym.vDerivProxy ){
      bool have = false;
      for( auto const& w : v0 ) if( w.id() == pr.id() ){ have = true; break; }
      if( !have ){ v0.push_back( pr ); x0.push_back( 0. ); }
    }
    // evaluate only if every variable is supplied; otherwise no zeroth-order block, QUIETLY (FFGraph::eval would
    // print "Subgraph evaluation failed -- missing variable ..." before throwing)
    bool complete = true;
    { FFSubgraph sg = _dag->subgraph( sym.vCoeff0.size(), sym.vCoeff0.data() );
      for( auto const* op : sg.l_op ){
        if( op->type != FFOp::VAR || !op->varout[0] ) continue;
        bool have = false;
        for( auto const& w : v0 ) if( w.id() == op->varout[0]->id() ){ have = true; break; }
        if( !have ){ complete = false; break; }
      } }
    if( complete ){
      try{
        std::vector<double> c0( nEqn * nState, 0. );
        _dag->eval( sym.vCoeff0, c0, v0, x0 );
        *A0 = arma::mat( (arma::uword)nEqn, (arma::uword)nState, arma::fill::zeros );
        for( size_t k = 0; k < nEqn; ++k )
          for( size_t j = 0; j < nState; ++j ) (*A0)( (arma::uword)k, (arma::uword)j ) = c0[ k * nState + j ];
      }
      catch( ... ){ A0->reset(); }   // a coefficient this context cannot evaluate: no zeroth-order block
    }
    else A0->reset();
  }
  if( sym.vCoeff.size() < nDom ) return std::vector<arma::mat>();

  std::vector<arma::mat> Ai( nDom, arma::mat( (arma::uword)nEqn, (arma::uword)nState, arma::fill::zeros ) );

  for( size_t i = 0; i < nDom; ++i ){
    std::vector<double> coeff_vals( nEqn * nState, 0. );
    _dag->eval( sym.vCoeff[i], coeff_vals, eval_vars, eval_vals );
    for( size_t k = 0; k < nEqn; ++k )
      for( size_t j = 0; j < nState; ++j )
        Ai[i]( (arma::uword)k, (arma::uword)j ) = coeff_vals[ k * nState + j ];
  }

  return Ai;
}

inline FFModel::t_Symbol const&
FFModel::_symbol_for_block
( int block_id )
const
{
  auto it = _blockSymbol.find( block_id );
  return it != _blockSymbol.end() ? it->second : _symbol;
}

inline bool
FFModel::_is_auxiliary_state
( FFVar const& state )
const
{
  return _auxVarID.find( state.id() ) != _auxVarID.end();
}

inline FFVar const*
FFModel::_auxiliary_parent
( FFVar const& state )
const
{
  for( auto const& aux : _auxDef ){
    if( aux.aux.id() != state.id() ) continue;
    // A default-constructed FFVar has id().second < 0 and no DAG binding.
    // Treat it as "no parent recorded" (should not happen for well-formed
    // auxiliary definitions, but guard against incomplete legacy data).
    if( aux.parent.id().second < 0 ) return nullptr;
    return &aux.parent;
  }
  return nullptr;
}

inline FFVar const*
FFModel::_auxiliary_primitive_root
( FFVar const& state )
const
{
  if( !_is_auxiliary_state( state ) ) return nullptr;
  FFVar const* cur = &state;
  // Walk the parent chain.  The chain length is bounded by the number of
  // reduce_order() passes (at most total derivative order minus one), so
  // this terminates quickly.  A cycle guard limits iterations to _auxDef.size().
  for( size_t guard = 0; guard <= _auxDef.size(); ++guard ){
    FFVar const* p = _auxiliary_parent( *cur );
    if( !p ) return cur;                 // parent not recorded — return cur as best guess
    if( !_is_auxiliary_state( *p ) ) return p;  // reached a primitive state
    cur = p;
  }
  return cur;  // safety fallback (should not be reached)
}

// ======================================================================
// OCFESLV::_time_derivative_state / _states_with_evolution_derivative
//   The single definition of "differential state" (OpP incidence in the evolution direction),
//   shared by the DAE decomposition's dyn set and marching-transfer auto-detection.
// ======================================================================
inline FFVar
FFModel::_time_derivative_state
( FFOp const* op, t_Var const& varmap )
const
{
  if( !op || !op->sameid( typeid( FFPartial ) ) || !_evolution_dom_set || !_evolution_dom_var.dag() ) return FFVar();
  auto const* pop = mc::type_cast<FFPartial const>( op );
  if( !pop ) return FFVar();
  auto const evo_id2 = _evolution_dom_var.id().second;   // match id().second's type (signed) exactly
  for( size_t jj = 0; jj < op->varin.size(); ++jj ){
    FFVar const* operand = op->varin[jj];
    if( !operand || varmap.find( *operand ) == varmap.end() ) continue;   // state operand only
    for( auto const& [indep_var, ord] : pop->Indep().expr ){ (void)ord;
      if( indep_var.id().second == evo_id2 ) return *operand;             // OpP(state, evolution_var)
    }
  }
  return FFVar();
}

inline std::set<FFVar,lt_FFVar>
FFModel::_states_with_evolution_derivative
( t_Eqns const& eqns, t_Var const& varmap )
const
{
  std::set<FFVar,lt_FFVar> out;
  if( !_evolution_dom_var.dag() ) return out;
  for( auto const& eqn : eqns ){
    FFGraph* dag = static_cast<FFGraph*>( eqn.var.dag() );
    if( !dag ) continue;
    FFVar expr = eqn.var;
    FFSubgraph sg = dag->subgraph( 1, &expr );
    for( auto const* op : sg.l_op ){
      FFVar const st = _time_derivative_state( op, varmap );
      if( st.dag() ) out.insert( st );
    }
  }
  return out;
}

// ======================================================================
// OCFESLV::_structural_dae_decomposition   (Stage 1a, structural/reference-free)
// ======================================================================
inline FFModel::t_StructuralDecomp
FFModel::_structural_dae_decomposition
( int const block_id, FFVar const* dir )
const
{
  t_StructuralDecomp d;

  std::set<FFVar,lt_FFVar> states_any, states_dt;

  for( auto const& eqn : _mEqn ){
    if( !eqn.opt->participate_in_classification || eqn.opt->block_id != block_id )
      continue;
    FFVar const ev = eqn.var;

    bool has_dt = false, has_spatial = false;
    auto sg = _dag->subgraph( 1, &ev );
    for( auto const& op : sg.l_op ){
      // Any state operand (bare or differentiated) counts as "appearing".
      for( size_t jj = 0; jj < op->varin.size(); ++jj ){
        FFVar const* v = op->varin[jj];
        if( _mVar.find( *v ) != _mVar.end() ) states_any.insert( *v );
      }
      // FFPartial: classify the derivative direction (time vs spatial).  After
      // _reduce_order's conservative-form normalisation every partial is
      // OpP(state,dir), so the differentiated operand is a bare state -- no
      // operand-subgraph walk is needed.  Time-derivative detection uses the shared
      // _time_derivative_state rule (same definition as marching auto-detection).
      if( !op->sameid( typeid(FFPartial) ) ) continue;
      auto const* pop = mc::type_cast<FFPartial const>( op );
      { bool const dir_is_evol = ( dir && dir->dag() && _evolution_dom_set && _evolution_dom_var.dag()
                                   && dir->id().second == _evolution_dom_var.id().second );
        FFVar st_dt;
        if( dir_is_evol || !dir ) st_dt = _time_derivative_state( op, _mVar );
        else                                          // the same question, asked of another direction
          for( size_t jj = 0; jj < op->varin.size() && !st_dt.dag(); ++jj ){
            FFVar const* operand = op->varin[jj];
            if( !operand || _mVar.find( *operand ) == _mVar.end() ) continue;
            for( auto const& [iv, ord] : pop->Indep().expr ){ (void)ord;
              if( iv.id().second == dir->id().second ){ st_dt = *operand; break; } }
          }
        if( st_dt.dag() ){ has_dt = true; states_dt.insert( st_dt ); } }
      for( size_t jj = 0; jj < op->varin.size(); ++jj ){
        FFVar const* operand = op->varin[jj];
        if( _mVar.find( *operand ) == _mVar.end() ) continue;       // state operand
        for( auto const& [indep_var, ord] : pop->Indep().expr ){ (void)ord;
          bool const is_time = ( dir && dir->dag() && indep_var.id().second == dir->id().second );
          if( !is_time && _mDom.find( indep_var ) != _mDom.end() ) has_spatial = true;
        }
      }
    }

    if( !has_dt ){
      d.constraints.push_back( ev );
      d.constraint_has_spatial.push_back( has_spatial ? 1 : 0 );
      d.constraint_is_link.push_back( eqn.opt->role == EqnRole::LINK ? 1 : 0 );
    }
  }

  // dyn = states with d_t; alg = appearing states with no d_t.  Map order.
  for( auto const& [v, sd] : _mVar ){
    if( states_any.find( v ) == states_any.end() ) continue;
    if( states_dt.find( v ) != states_dt.end() ) d.dyn.push_back( v );
    else                                         d.alg.push_back( v );
  }

  return d;
}

// ======================================================================
// OCFESLV::_structural_index_analysis   (Stage 1a, cached wrapper)
// Memoise the per-block analysis so the classify pass, the structural probe,
// and the Stage-2 reduction-plan builder share ONE computation.  See _idxCache.
// ======================================================================
inline FFModel::t_IndexResult
FFModel::_structural_index_analysis
( int const block_id )
const
{
  auto it = _idxCache.find( block_id );
  if( it != _idxCache.end() ) return it->second;
  t_IndexResult R = _structural_index_analysis_uncached( block_id, _evolution_dom_set? &_evolution_dom_var: nullptr );
  _idxCache[ block_id ] = R;
  return R;
}

// ======================================================================
// OCFESLV::_structural_index_analysis_uncached   (Stage 1a, structural matching)
// ======================================================================
inline FFModel::t_IndexResult
FFModel::_structural_index_analysis_uncached
( int const block_id, FFVar const* dir )
const
{
  t_IndexResult R;
  t_StructuralDecomp const d = _structural_dae_decomposition( block_id, dir );

  // Parabolic character, source 2 (index-independent): a SPATIAL LINK
  //   Daux - d_x(state) = 0
  // is the reduce_order signature of a >=2nd-order spatial (diffusion) term.
  // It makes the block parabolic regardless of DAE index, and -- unlike the
  // index-1 spatial-constraint closure detected later (source 1) -- it is
  // present even when the block is structurally index-0 (the aux is a LINK).
  // This is what the old has_link && has_volume heuristic was really catching
  // (e.g. heat d_t u = D d_zz u: index 0 via the Dz_u LINK, yet PARABOLIC),
  // refined here to require the LINK be SPATIAL (a 2nd-order TIME reduction is
  // not parabolic).
  for( size_t i = 0; i < d.constraints.size(); ++i )
    if( d.constraint_is_link[i] && d.constraint_has_spatial[i] ){
      R.parabolic_character = true;                       // source 2
      break;
    }

  if( d.alg.empty() ){ R.index = 0; return R; }       // regular (no algebraic states)

  // Split constraints: INTERIOR algebraic (count toward index) vs LINK
  // (order-reduction auxiliary definitions, excluded from the index match).
  std::vector<FFVar> icon, lcon;
  std::vector<char>  icon_spatial;
  for( size_t i = 0; i < d.constraints.size(); ++i ){
    if( d.constraint_is_link[i] ) lcon.push_back( d.constraints[i] );
    else { icon.push_back( d.constraints[i] ); icon_spatial.push_back( d.constraint_has_spatial[i] ); }
  }

  // Structural "bare" (algebraic, non-derivative) coupling.  A state couples
  // algebraically to an equation iff it is a direct operand of a NON-external
  // op; a state appearing only under d/dx(.) or integral(.) is a derivative /
  // nonlocal coupling, not algebraic.  After _reduce_order's conservative-form
  // normalisation every FFPartial operand is a bare state and every op in the
  // equation subgraph feeds the root, so this flat incidence scan is exact --
  // no reachability pruning is needed.  No differentiation (FAD throws on
  // external ops); this is pure incidence.
  auto bare_states = [&]( FFVar const& eqn ) -> std::set<FFVar,lt_FFVar> {
    std::set<FFVar,lt_FFVar> bare;
    auto sg = _dag->subgraph( 1, &eqn );
    for( auto const& op : sg.l_op ){
      if( op->sameid( typeid(FFPartial) ) || op->sameid( typeid(FFIntegral) ) )
        continue;                                         // operands are under a derivative/integral
      for( auto const* in : op->varin )
        if( in && _mVar.find( *in ) != _mVar.end() ) bare.insert( *in );
    }
    return bare;
  };

  // Algebraic variables coupled (bare) by a LINK are order-reduction
  // auxiliaries: drop them from the index-relevant set.  A GENUINE algebraic
  // state (e.g. a wetting law omega=g(...)) can ALSO appear bare in a LINK -- as
  // a coefficient of an interface-flux constraint -- yet it is not an auxiliary
  // and IS index-relevant (pinned by its own INTERIOR closure).  Gate on the
  // auxiliary registry so such value-slaved algebraic states are retained.
  std::set<FFVar,lt_FFVar> link_aux;
  for( auto const& lc : lcon ){
    auto const bs = bare_states( lc );
    for( auto const& a : d.alg )
      if( bs.find( a ) != bs.end() && _is_auxiliary_state( a ) ) link_aux.insert( a );
  }
  for( auto const& a : d.alg )
    if( link_aux.find( a ) == link_aux.end() ) R.index_alg.push_back( a );

  if( R.index_alg.empty() ){ R.index = 0; return R; }  // all algebraic vars were order-reduction aux

  size_t const nC = icon.size(), nA = R.index_alg.size();

  // Structural algebraic incidence: interior constraint c pins algebraic var a
  // iff a appears bare (non-derivative) in c.
  std::vector< std::vector<char> > inc( nC, std::vector<char>( nA, 0 ) );
  for( size_t c = 0; c < nC; ++c ){
    auto const bs = bare_states( icon[c] );
    for( size_t a = 0; a < nA; ++a )
      inc[c][a] = ( bs.find( R.index_alg[a] ) != bs.end() ) ? 1 : 0;
  }

  // Maximum bipartite matching (constraints -> algebraic variables), Kuhn.
  std::vector<int> matchA( nA, -1 );
  std::function<bool(size_t,std::vector<char>&)> aug =
    [&]( size_t c, std::vector<char>& seen ) -> bool {
      for( size_t a = 0; a < nA; ++a ){
        if( !inc[c][a] || seen[a] ) continue;
        seen[a] = 1;
        if( matchA[a] < 0 || aug( (size_t)matchA[a], seen ) ){ matchA[a] = (int)c; return true; }
      }
      return false;
    };
  for( size_t c = 0; c < nC; ++c ){ std::vector<char> seen( nA, 0 ); aug( c, seen ); }

  size_t nmatched = 0;
  for( size_t a = 0; a < nA; ++a ) if( matchA[a] >= 0 ) ++nmatched;

  if( nmatched == nA ){
    R.index = 1;                                       // one elimination closes every algebraic var
    for( size_t a = 0; a < nA; ++a )
      if( matchA[a] >= 0 && icon_spatial[ matchA[a] ] ) R.parabolic_character = true;
  }
  else {
    // ---- Stage 2: resolve the EXACT index by structural-differentiation
    // reachability (incidence-only; FAD throws on the external d_t ops).  An
    // unmatched algebraic var a_u is exposed only by differentiating a
    // constraint until a_u appears in it: differentiating a constraint that
    // contains a dynamic state s introduces s', and s' is fixed by s's ODE, so
    // the differentiated constraint structurally acquires the bare states of
    // that ODE.  Chaining through the dynamic states, index = 1 + (rounds of
    // differentiation to first expose a_u).  Exact on the index ladder
    // (M2 x'=y,0=x-a -> 1 diff -> 2; M3 x'=v,v'=L,0=x-a -> 2 diffs -> 3).  A
    // witness no differentiation can expose (structurally singular) leaves -1.
    // (Single-deficiency chains are exact; arbitrary coupled high-index DAEs
    // would need full Pantelides -- deferred, flagged where it matters.)
    for( size_t a = 0; a < nA; ++a )
      if( matchA[a] < 0 ) R.unmatched.push_back( R.index_alg[a] );

    FFVar const* tdom = dir;                       // rev339: the direction is a parameter

    // dyn_rhs[s] = bare states of s's defining (d_t s) equation -- what ONE
    // differentiation of a constraint containing s structurally exposes.
    std::map<FFVar,std::set<FFVar,lt_FFVar>,lt_FFVar> dyn_rhs;
    for( auto const& eqn : _mEqn ){
      if( !eqn.opt->participate_in_classification || eqn.opt->block_id != block_id )
        continue;
      std::set<FFVar,lt_FFVar> dt_here;
      auto sg = _dag->subgraph( 1, &eqn.var );
      for( auto const& op : sg.l_op ){
        if( !op->sameid( typeid(FFPartial) ) ) continue;
        auto const* pop = mc::type_cast<FFPartial const>( op );
        for( auto const* operand : op->varin ){
          if( !operand || _mVar.find( *operand ) == _mVar.end() ) continue;
          for( auto const& [iv, ord] : pop->Indep().expr )
            if( tdom && iv.id() == tdom->id() ) dt_here.insert( *operand );
        }
      }
      if( dt_here.empty() ) continue;                 // a constraint, not an ODE
      auto const bs = bare_states( eqn.var );
      for( auto const& s : dt_here )
        dyn_rhs[s].insert( bs.begin(), bs.end() );
    }

    std::set<FFVar,lt_FFVar> const dyn_set( d.dyn.begin(), d.dyn.end() );

    // Round-0 exposed set: bare states of all interior constraints.
    std::set<FFVar,lt_FFVar> exposed0;
    for( size_t c = 0; c < nC; ++c ){
      auto const bs = bare_states( icon[c] );
      exposed0.insert( bs.begin(), bs.end() );
    }
    size_t const cap = nA + d.dyn.size() + 2;         // termination bound
    size_t max_diffs = 0;
    bool   resolved  = true;
    for( size_t a = 0; a < nA; ++a ){
      if( matchA[a] >= 0 ) continue;                  // pinned at index 1
      FFVar const& au = R.index_alg[a];
      std::set<FFVar,lt_FFVar> exposed = exposed0;
      size_t diffs = 0;
      bool found = ( exposed.find( au ) != exposed.end() );
      while( !found && diffs < cap ){
        std::set<FFVar,lt_FFVar> next = exposed;
        for( auto const& s : exposed ){
          if( dyn_set.find( s ) == dyn_set.end() ) continue;
          auto it = dyn_rhs.find( s );
          if( it != dyn_rhs.end() ) next.insert( it->second.begin(), it->second.end() );
        }
        if( next.size() == exposed.size() ) break;    // no growth -> unexposable
        exposed.swap( next );
        ++diffs;
        found = ( exposed.find( au ) != exposed.end() );
      }
      if( !found ){ resolved = false; break; }
      max_diffs = std::max( max_diffs, diffs );
    }
    R.index = resolved ? (int)( max_diffs + 1 ) : -1;
  }
  return R;
}

// ======================================================================
// OCFESLV::classify
// ======================================================================
inline FFModel::t_Classify
FFModel::classify
( std::vector<arma::mat> const& Ai,
  std::vector<FFVar>     const& vDom,
  FFVar const*                  time_dom,
  unsigned const                n_sample,
  double   const                imag_tol,
  bool     const                has_algebraic_rows,
  std::vector<arma::mat> const* Ai_dn,
  arma::mat const*              E0_dn )
const
{
  t_Classify res;
  res.Ai       = Ai;
  res.imag_tol = imag_tol;

  size_t const nDom = vDom.size();
  // rev67: separate "no principal symbol AT ALL" from "analysis inconclusive".
  //   nDom == 0                      -> zero collocation domains: ALGEBRAIC_LUMPED
  //   Ai.empty() with nDom > 0       -> distributed states, no derivative anywhere:
  //                                     ALGEBRAIC_FIELD
  // Both are CONCLUSIVE.  UNDETERMINED is reserved for its documented meaning,
  // "non-square system or analysis inconclusive", and the two carry different
  // downstream consequences: an algebraic system needs no initial condition, has
  // nothing to march, and has no characteristic treatment.
  if( nDom == 0 ){
    res.type = ALGEBRAIC_LUMPED;
    return res;
  }
  if( Ai.empty() ){
    res.type = ALGEBRAIC_FIELD;
    return res;
  }

  size_t const nRows = Ai[0].n_rows;
  size_t const nCols = Ai[0].n_cols;
  // rev70: name the rectangular case instead of leaving it as UNDETERMINED.  The
  // recovery path below (drop all-zero rows if that squares the system) is tried FIRST
  // and unchanged; only if it does not apply do we record the shape.
  if( nRows != nCols ){
    // A RECTANGULAR principal symbol means the block carries more interior
    // equations than derivative-bearing states.  The excess equations are
    // algebraic constraints (e.g. a value-slaved EOS) whose principal rows are
    // identically zero (no highest-order derivative term) and carry no spectral
    // character.  Drop the rows that are all-zero across EVERY domain; if that
    // squares the symbol, classify the reduced (differential) subsystem so the
    // eigenanalysis decides the character (parabolic / dispersive / hyperbolic)
    // instead of bailing to UNDETERMINED.  This only ever affects a block that
    // would otherwise be UNDETERMINED -- square symbols are untouched -- and the
    // verdict comes from the spectrum, never a heuristic, so it cannot mislabel
    // (a dispersive chain still resolves to an imaginary-spectrum verdict).
    if( nRows > nCols ){
      std::vector<arma::uword> keep;
      keep.reserve( nRows );
      for( arma::uword r = 0; r < nRows; ++r ){
        bool all_zero = true;
        for( auto const& A : Ai ){
          for( arma::uword c = 0; c < A.n_cols; ++c )
            if( A(r,c) != 0. ){ all_zero = false; break; }
          if( !all_zero ) break;
        }
        if( !all_zero ) keep.push_back( r );
      }
      if( keep.size() == nCols && keep.size() < nRows ){
        std::vector<arma::mat> Ai_sq;
        Ai_sq.reserve( Ai.size() );
        for( auto const& A : Ai ){
          arma::mat B( keep.size(), A.n_cols );
          for( size_t i = 0; i < keep.size(); ++i ) B.row( i ) = A.row( keep[i] );
          Ai_sq.push_back( std::move( B ) );
        }
        return classify( Ai_sq, vDom, time_dom, n_sample, imag_tol,
                         has_algebraic_rows );
      }
    }
    // rev70: a rectangular differential core is a KNOWN shape, not an inconclusive
    // analysis.  Naming it lets the evolution check distinguish "we could not tell" from
    // "the symbol cannot be square here", which carry different consequences for marching.
    // rev162 STAGE A (knob CRONOS_RECT_RECORD=1, default 0 = rev161a).
    // This return SKIPS everything below it, including the evolution-domain recording at the
    // end of classify().  `evolution_user_supplied` is assigned ONLY there, i.e. after this
    // point, so on this path it is ALWAYS false -- which makes the rev70 suitability clause
    // `type == DIFFERENTIAL_RECTANGULAR && evolution_user_supplied` UNSATISFIABLE.  MEASURED
    // corpus-wide (rev160 sweep 20260907_202623): 136 of 136 rectangular classify lines report
    // evol_dom_idx=none, zero report [user-supplied], and the "AUTO-DETECTED" note fires 20
    // times -- including on drivers that DO call set_evolution_domain().
    // The evolution recording is SVD-based and VALID on a rectangular matrix; the
    // eigen/character sweep below is not (it needs a square pencil).  STAGE A performs exactly
    // the valid half here and leaves the type untouched.  It also PRINTS the type this block
    // would receive if the type were not overwritten, applying nothing, so STAGE C's blast
    // radius is readable before STAGE C exists.
    res.symbol_rectangular = true;
    {
      // rev168: unconditional.  STAGE A's recording was accepted on sweep 20260908_132337
      // (plan-neutral, 1126/1126) and STAGE C -- now the only path -- requires it.
      {
        size_t i_time_r = nDom;
        if( time_dom )
          for( size_t i = 0; i < nDom; ++i )
            if( vDom[i].id() == time_dom->id() ){ i_time_r = i; break; }
        bool at_sing = true;
        if( i_time_r < nDom ){
          res.evolution_dom_idx       = i_time_r;
          res.evolution_user_supplied = true;
          arma::mat const& At = Ai[ i_time_r ];
          if( At.n_rows && At.n_cols ){
            arma::vec sv;
            arma::svd( sv, At );                      // valid for a rectangular At
            if( sv.n_elem ){
              double const smax = sv(0), smin = sv(sv.n_elem-1);
              double const tol  = imag_tol * std::max( 1.0, smax );
              arma::uword rk = 0;
              for( arma::uword i = 0; i < sv.n_elem; ++i ) if( sv(i) > tol ) ++rk;
              res.rank_evolution      = rk;
              res.sigma_min_evolution = smin;
              res.sigma_max_evolution = smax;
              res.cond_evolution      = ( smin > 0. ) ? smax / smin
                                        : std::numeric_limits<double>::infinity();
              at_sing = ( rk < std::min( At.n_rows, At.n_cols ) );  // mirrors rank_info()
              res.At_singular = at_sing;
            }
          }
        }
        std::cerr << "OCFESLV::classify ** RECT(A): core " << nRows << "x" << nCols
                  << " time_dom=" << ( time_dom ? "y" : "n" )
                  << " evol_dom_idx=" << ( i_time_r < nDom ? std::to_string( i_time_r )
                                                           : std::string( "none" ) )
                  << " At_singular=" << ( i_time_r < nDom ? ( at_sing ? "y" : "n" ) : "?" )
                  << " rank_evol=" << res.rank_evolution
                  << " | would-be type="
                  << ( i_time_r >= nDom ? "UNDETERMINED (no evolution domain)"
                     : at_sing         ? "DESCRIPTOR (then PARABOLIC iff parabolic_character)"
                                       : "UNDETERMINED (At non-singular; needs a square test)" )
                  << "  [recorded and applied]" << std::endl;
      }
    }
    // rev164 STAGE C (knob CRONOS_RECT_TYPE_FROM_CHARACTER=1, default 0 = rev163).
    // Stop writing a SHAPE fact into the CHARACTER axis.  `symbol_rectangular` (STAGE A) now
    // carries the shape and the three consumers read it (STAGE B), so the type is free to say
    // what the block IS:
    //   At_singular  -> DESCRIPTOR: a rectangular differential core whose evolution block is
    //                   rank-deficient IS a descriptor system, and the existing structural
    //                   path then upgrades it to PARABOLIC exactly when parabolic_character
    //                   holds -- setting parabolic_structure_detected as a by-product, which
    //                   is what _block_drop_eligible() reads.
    //   otherwise    -> UNDETERMINED: honest.  The square-symbol test does not apply and there
    //                   is no square test to give; suitability is preserved through
    //                   symbol_rectangular, not through the type.
    // CORPUS CENSUS (rev162 sweep 20260908_132337, 148 RECT(A) lines / 18 drivers):
    //   86 lines At_singular=y  -> DESCRIPTOR (13 drivers)
    //   61 lines At_singular=n  -> UNDETERMINED (PDE39, PDE39b, DAE5; all core 1x2)
    //    2 lines no evolution domain -> UNDETERMINED (autoelim)
    // rev168 STAGE C DEFAULTED ON, and STAGE D5: the type is set from structure
    // unconditionally and DIFFERENTIAL_RECTANGULAR is no longer assigned anywhere.
    // ACCEPTED on corpus sweep 20260908_234634: plan-neutral (1199/1199 determinacy lines
    // identical), the drop stayed locked, and PDE39/PDE39b/DAE5 took the new UNDETERMINED
    // path with no consequence.  CRONOS_RECT_TYPE_FROM_CHARACTER is retired with it.
    {
      bool at_sing_c = false;
      if( res.evolution_dom_idx != std::numeric_limits<size_t>::max() ) at_sing_c = res.At_singular;
      res.type       = at_sing_c ? DESCRIPTOR : UNDETERMINED;
      res.descriptor = at_sing_c;              // mirrors the square DESCRIPTOR branch below
      std::cerr << "OCFESLV::classify ** RECT(C): type set from structure -> "
                << pde_type_name( res.type )
                << " (symbol_rectangular=1, At_singular=" << ( at_sing_c ? "y" : "n" ) << ")"
                << std::endl;
      return res;
    }
  }
  // rev168 STAGE D3: past every early return (algebraic-lumped, algebraic-field, rectangular),
  // so the principal symbol IS formed and examined for this block from here on -- including the
  // DESCRIPTOR branch, whose singular-A_t verdict is a measurement even though the direction
  // sweep is not reached.  MEASURED while placing this: setting the flag at the direction sweep
  // instead left PDE31b's blk0 reading analysed=0 while carrying type=PARABOLIC from the
  // structural upgrade -- accurate about the sweep, wrong about the documented meaning.
  res.symbol_analysed = true;

  size_t const n = nRows;

  struct RankInfo {
    arma::uword rank = 0;
    double smin = 0.;
    double smax = 0.;
    double cond = std::numeric_limits<double>::infinity();
    bool singular = true;
  };

  auto rank_info = [&]( arma::mat const& A ) -> RankInfo
  {
    RankInfo r;
    if( A.n_rows == 0 || A.n_cols == 0 ) return r;
    arma::vec s;
    arma::svd( s, A );
    if( s.n_elem ){
      r.smax = s(0);
      r.smin = s(s.n_elem-1);
      double const tol = imag_tol * std::max( 1.0, r.smax );
      for( arma::uword i = 0; i < s.n_elem; ++i )
        if( s(i) > tol ) ++r.rank;
      r.singular = ( r.rank < std::min( A.n_rows, A.n_cols ) );
      r.cond = ( r.smin > 0. ) ? r.smax / r.smin
                               : std::numeric_limits<double>::infinity();
    }
    return r;
  };

  auto make_dirs = [&]( size_t dim ) -> std::vector< std::vector<double> >
  {
    std::vector< std::vector<double> > dirs;
    if( dim == 0 ){
      dirs.push_back( {} );
    }
    else if( dim == 1 ){
      dirs.push_back( {1.0} );
    }
    else if( dim == 2 ){
      for( unsigned s = 0; s < n_sample; ++s ){
        double theta = s * 3.14159265358979323846 / static_cast<double>( n_sample );   // pi (M_PI is not standard C++)
        dirs.push_back( { std::cos(theta), std::sin(theta) } );
      }
    }
    else {
      arma::arma_rng::set_seed( 42u );
      for( unsigned s = 0; s < n_sample; ++s ){
        arma::vec v = arma::randn<arma::vec>( (arma::uword)dim );
        double nv = arma::norm(v);
        if( nv <= 0. ) continue;
        v /= nv;
        std::vector<double> d( dim );
        for( size_t k = 0; k < dim; ++k ) d[k] = v(k);
        dirs.push_back( d );
      }
    }
    return dirs;
  };

  auto sweep_evolution = [&]( size_t i_evo,
                              std::vector<size_t> const& i_space,
                              bool record )
    -> std::pair<bool,bool>
  {
    bool all_real = true;
    bool any_real = false;
    bool weak     = false;
    auto dirs = make_dirs( i_space.size() );

    for( auto const& xi : dirs ){
      arma::mat S( (arma::uword)n, (arma::uword)n, arma::fill::zeros );
      for( size_t k = 0; k < i_space.size(); ++k )
        S += xi[k] * Ai[ i_space[k] ];

      arma::mat M = arma::solve( Ai[i_evo], S );
      arma::cx_vec eigvals;
      arma::cx_mat eigvec;
      bool ok = arma::eig_gen( eigvals, eigvec, M );
      if( !ok ){
        all_real = false;
        weak = true;
      }
      else {
        for( arma::uword e = 0; e < eigvals.n_elem; ++e ){
          double const im = std::abs( eigvals(e).imag() );
          res.max_imag_eig = std::max( res.max_imag_eig, im );
          if( im > imag_tol ) all_real = false;
          else                any_real = true;
        }
        double cnd = arma::cond( eigvec );
        if( !std::isfinite(cnd) ) cnd = std::numeric_limits<double>::infinity();
        res.max_eigvec_cond = std::max( res.max_eigvec_cond, cnd );
        if( !std::isfinite(cnd) || cnd > 1.0 / std::max( imag_tol, 1e-15 ) )
          weak = true;
      }

      if( record ){
        res.eigenvalues = eigvals;
        res.eigendata.push_back( { xi, eigvals } );
      }
    }

    return { all_real && any_real, weak };
  };

  // Locate the caller-supplied evolution/time variable, if any.
  size_t i_time = nDom;
  if( time_dom )
    for( size_t i = 0; i < nDom; ++i )
      if( vDom[i].id() == time_dom->id() ){ i_time = i; break; }

  if( i_time < nDom ){
    res.evolution_dom_idx       = i_time;
    res.evolution_user_supplied = true;

    RankInfo const ri = rank_info( Ai[i_time] );
    res.rank_evolution      = ri.rank;
    res.sigma_min_evolution = ri.smin;
    res.sigma_max_evolution = ri.smax;
    res.cond_evolution      = ri.cond;
    res.At_singular         = ri.singular;

    std::vector<size_t> i_space;
    for( size_t i = 0; i < nDom; ++i )
      if( i != i_time ) i_space.push_back( i );

    if( i_space.empty() ){
      // rev68: A_e's rank is NOT the whole test.
      //
      // The principal symbol is built over the DIFFERENTIAL SUBSYSTEM ONLY -- this
      // header's own t_Symbol documentation says algebraic equations are "excluded from
      // the differential symbol (they would be zero rows)" and vState holds "state
      // variables appearing in at least one derivative term".  So for a genuine DAE the
      // zero row never reaches A_e, ri.singular is FALSE, and the verdict came out
      // DIFFERENTIAL_ORDINARY.
      //
      // MEASURED on OCFE_DAE0 (index-3 pendulum, index-1 reduced): five states
      // (x,y,u,v,lam) and five interior equations -- four evolution plus one algebraic
      // constraint with no OpP -- yet sym=4x4, Ae_singular=n, sigma_min=sigma_max=1.
      // A_e is the 4x4 identity on the differential states, which is CORRECT for the
      // reduced symbol and wrong as a statement about the model: the reduction discarded
      // exactly the feature that makes it a DAE.  OCFE_ODE0 is the control at sym=3x3
      // with 3 states, where the identity is the right answer.
      //
      // vAlgEqn is already computed and holds those excluded rows, so the test costs
      // nothing: a non-empty algebraic block means the system is differential-algebraic
      // regardless of what the differential part's rank says.  The singular-A_e route is
      // kept because a rank-deficient differential block is a DAE too.
      // rev344: a NONSINGULAR but coupled evolution matrix is an IMPLICIT ODE, not an ordinary one.  The test is
      // structural -- one nonzero per row AND per column, i.e. a scaled permutation -- because a non-unit
      // coefficient (2*x' = f) and a permuted assignment are both explicit and neither is diagonal.
      bool decoupled = true;
      if( time_dom ){
        size_t evo = vDom.size();
        for( size_t d = 0; d < vDom.size(); ++d ) if( vDom[d].id() == time_dom->id() ) evo = d;
        if( evo < Ai.size() ){
          arma::mat const& Ae = Ai[evo];
          double const scale = Ae.n_elem? arma::abs( Ae ).max(): 0.;
          double const tol   = 1e-10 * std::max( 1.0, scale );
          for( arma::uword i = 0; i < Ae.n_rows; ++i ){
            size_t nz = 0;
            for( arma::uword j = 0; j < Ae.n_cols; ++j ) if( std::fabs( Ae(i,j) ) > tol ) ++nz;
            if( nz > 1 ) decoupled = false;
          }
          for( arma::uword j = 0; j < Ae.n_cols; ++j ){
            size_t nz = 0;
            for( arma::uword i = 0; i < Ae.n_rows; ++i ) if( std::fabs( Ae(i,j) ) > tol ) ++nz;
            if( nz > 1 ) decoupled = false;
          }
        }
      }
      res.type = ( ri.singular || has_algebraic_rows ) ? DIFFERENTIAL_ALGEBRAIC
               : ( decoupled? DIFFERENTIAL_ORDINARY : DIFFERENTIAL_IMPLICIT );
      res.descriptor = ri.singular;
      return res;
    }

    if( ri.singular ){
      // Singular A_e is a descriptor system.  classify_pde() may upgrade this
      // to PARABOLIC when equation metadata shows a diffusion/link closure.
      res.type       = DESCRIPTOR;
      res.descriptor = true;
      return res;
    }

    auto sw = sweep_evolution( i_time, i_space, true );
    if( sw.first ){
      res.evolution_hyperbolic = true;
      res.weak_hyperbolic     = sw.second;
      // rev167 STAGE D1: defectiveness is a QUALIFIER on hyperbolicity, not a character in its
      // own right, and every consumer that cares already reads the FLAG -- the upwind decision
      // tests `evolution_hyperbolic && !weak_hyperbolic` at three sites.  Writing it into the
      // type displaced EVOL_HYPERBOLIC and dropped the block out of the suitability whitelist,
      // so the encoding could only ever do harm.  MEASURED: WEAK_HYPERBOLIC occurs 0 times in
      // 1410 corpus classification lines, and OCFE_SYMBOLCASES1 (Jordan block) confirms the
      // detection itself works and sets the flag.
      res.type = EVOL_HYPERBOLIC;
    }
    else if( !res.eigendata.empty() && res.max_imag_eig <= imag_tol ){
      res.type = EVOL_HYPERBOLIC;
      res.evolution_hyperbolic = true;
    }
    else {
      res.type = res.max_imag_eig > imag_tol ? COMPLEX_CHARACTERISTIC : UNDETERMINED;
    }
    return res;
  }

  // No evolution direction was supplied.  First try to detect a steady
  // marching/evolution direction from an invertible coefficient matrix.  This
  // enables steady hyperbolic first-order systems without treating every real
  // characteristic variety as causal.
  if( nDom > 1 ){
    for( size_t i_evo = 0; i_evo < nDom; ++i_evo ){
      RankInfo const ri = rank_info( Ai[i_evo] );
      if( ri.singular ) continue;
      std::vector<size_t> i_space;
      for( size_t i = 0; i < nDom; ++i )
        if( i != i_evo ) i_space.push_back( i );
      if( i_space.empty() ) continue;

      t_Classify trial = res;
      res.eigendata.clear();
      res.max_imag_eig = 0.;
      res.max_eigvec_cond = 0.;
      auto sw = sweep_evolution( i_evo, i_space, true );
      if( sw.first ){
        res.evolution_dom_idx       = i_evo;
        res.evolution_auto_detected = true;
        res.evolution_hyperbolic    = true;
        res.weak_hyperbolic         = sw.second;   // rev167 D1: flag only, see the site above
        res.rank_evolution          = ri.rank;
        res.sigma_min_evolution     = ri.smin;
        res.sigma_max_evolution     = ri.smax;
        res.cond_evolution          = ri.cond;
        res.type = EVOL_HYPERBOLIC;                // rev167 D1
        return res;
      }
      res = trial;
    }
  }

  // No causal/evolution direction is available: classify the spatial
  // principal symbol P(\xi)=sum_i xi_i A_i by rank.  Real characteristic
  // directions are recorded as SPATIALLY_CHARACTERISTIC, not as upwindable
  // EVOL_HYPERBOLIC.
  bool any_characteristic = false;
  bool any_nonsingular    = false;
  for( auto const& xi : make_dirs( nDom ) ){
    arma::mat P( (arma::uword)n, (arma::uword)n, arma::fill::zeros );
    for( size_t i = 0; i < nDom; ++i ) P += xi[i] * Ai[i];

    arma::cx_vec eigvals;
    arma::eig_gen( eigvals, P );
    RankInfo const ri = rank_info( P );
    if( ri.singular ) any_characteristic = true;
    else              any_nonsingular    = true;
    res.eigenvalues = eigvals;
    res.eigendata.push_back( { xi, eigvals } );
  }

  // An order-reduced system: the first-derivative symbol is singular for every direction (the LINK rows of one
  // primitive are parallel), so test the Douglis-Nirenberg WEIGHTED symbol -- the LINK rows' zeroth-order entries
  // are principal there.  Nonsingular in every sampled direction: ELLIPTIC (report-only, see t_Classify::dn_elliptic).
  if( any_characteristic && Ai_dn && E0_dn && Ai_dn->size() == nDom && E0_dn->n_rows == (arma::uword)n ){
    bool dn_nonsingular = true;
    for( auto const& xi : make_dirs( nDom ) ){
      arma::mat S = *E0_dn;
      for( size_t i = 0; i < nDom; ++i ) S += xi[i] * (*Ai_dn)[i];
      if( rank_info( S ).singular ){ dn_nonsingular = false; break; }
    }
    if( dn_nonsingular ){
      res.type = ELLIPTIC;
      res.dn_elliptic = true;
      return res;
    }
  }

  if( !any_characteristic ){
    res.type = ELLIPTIC;
  }
  else if( !any_nonsingular ){
    res.type = DEGENERATE;
    res.degenerate = true;
  }
  else {
    res.type = SPATIALLY_CHARACTERISTIC;
    res.spatially_characteristic = true;
  }

  return res;
}

// ======================================================================
// OCFESLV automatic classification helpers
// ======================================================================
inline FFVar const*
FFModel::_infer_evolution_domain
()
{
  if( _evolution_dom_user && _evolution_dom_set )
    return &_evolution_dom_var;

  // Prefer domains pinned by INITIAL equations.  This catches the common
  // transient case t = LB/UB without requiring the user to name the domain.
  for( auto const& eqn : _mEqn ){
    auto const& eqndom = eqn.dom;
    auto const& opt    = *eqn.opt;
    if( opt.role != EqnRole::INITIAL ) continue;
    for( auto const& [v,lim] : eqndom ){
      if( lim == FFDom::LB || lim == FFDom::UB ){
        _evolution_dom_var = v;
        _evolution_dom_set = true;
        _evoInferProvenance = 1;   // rev99: STRUCTURAL -- an INITIAL row pins this domain
        return &_evolution_dom_var;
      }
    }
  }

  // Fall back to common names.  If no evolution domain is known, nullptr is
  // intentional: classify() then performs the no-time/spatial analysis.
  for( auto const& [v,dom] : _mDom ){
    (void)dom;
    std::string const nm = v.name();
    if( nm == "t" || nm == "T" || nm == "time" || nm == "Time" ){
      _evolution_dom_var = v;
      _evolution_dom_set = true;
      // ===================================================================
      // rev99: WEAK provenance -- a NAMING CONVENTION, not a structural fact
      // ===================================================================
      // The rule above finds the domain an INITIAL equation actually pins.  This one
      // matches a domain because of what it is CALLED.  Both set the same member and
      // rev98 reported them identically, so a plan built on a naming coincidence was
      // indistinguishable from one built on the model's own structure.
      //
      // That distinction did not matter while CRONOS_HOIST_EVODETECT defaulted to 0,
      // because nothing acted on the inference.  It matters now: at the promoted default
      // the inference decides whether the solve MARCHES, and marching in a wrong causal
      // direction is SILENTLY wrong where monolithic is merely slow.  No driver in the
      // 139-target corpus reaches this branch with a domain that is not genuinely the
      // evolution direction -- which is not the same as the branch being safe.
      //
      // So it is recorded and reported rather than trusted.  A model with a spatial
      // domain named `t`, or a transient whose IC is imposed some other way, lands here.
      _evoInferProvenance = 2;   // NAME-BASED
      return &_evolution_dom_var;
    }
  }

  _evoInferProvenance = 0;

  if( !_evolution_dom_user )
    _evolution_dom_set = false;
  return nullptr;
}

// ======================================================================
// OCFESLV::_probe_reference_robustness   (opt-in, print-only)
// ======================================================================
// All setup-time structural analysis (principal symbol / type, characteristic
// directions, auto-closure, incoming-BC guard, Peclet corrector, SAT scaling)
// linearizes ONCE at _classification_reference.  The symbol is sampled over many
// directions xi but at a SINGLE state -- directional robustness, no STATE
// robustness.  For a nonlinear/quasilinear A_z(u) the reference-linearized
// character can differ from the operating regime; and unlike the per-iteration
// Newton/LM re-linearization (which keeps the SOLVE correct), a wrong STRUCTURE
// is baked into the assembled system at setup and cannot be recovered.
//
// This probe re-evaluates the SAME _principal_symbol at several state references
// and reports, in priority order:
//   GATE 1 (hard) : PDE type invariant                  -> closure strategy stable
//   GATE 2 (hard) : per-face incoming/outgoing split invariant
//                   (a characteristic sign flip moves a row A+ <-> A-) -> the
//                   incoming/outgoing assignment the closure+guard rely on is stable
//   advisory      : the incoming-characteristic subspace does not rotate beyond tol
// It is PRINT-ONLY and returns the all-invariant verdict (a caller may later
// promote GATEs 1/2 to a hard setup failure).  The base symbol/classification are
// reused from the cache; only _eval_symbol/classify/_build_face_data re-run per
// sample (no DAG growth).  Samples = base + user-supplied envelope states
// (add_classification_reference_sample) or, if none, a LOCAL magnitude sweep of
// the base (sensitivity only -- a genuine sign flip needs envelope states).
inline bool
FFModel::_probe_reference_robustness
( std::vector<double> const& state_ref,
  std::vector<double> const& input_ref,
  std::vector<double> const& cst_ref,
  std::vector<double> const& dom_ref,
  FFVar const*               time_dom,
  unsigned const             n_sample,
  double   const             imag_tol ) const
{
  // -- assemble the sample set -----------------------------------------------
  // Each user sample is a partial map {state FFVar -> value}; resolve each key to
  // its _vVar index and override that entry of a copy of the base reference.
  std::vector<std::vector<double>> samples;
  std::vector<std::string>         labels;
  samples.push_back( state_ref );  labels.push_back( "base" );
  size_t n_user = 0;
  for( auto const& m : _classRefSamples ){
    std::vector<double> s = state_ref;
    size_t n_resolved = 0;
    for( auto const& kv : m ){
      size_t idx = _vVar.size();
      for( size_t j = 0; j < _vVar.size(); ++j )
        if( _vVar[j].id() == kv.first.id() ){ idx = j; break; }
      if( idx < s.size() ){ s[idx] = kv.second; ++n_resolved; }
      else std::cerr << "  WARN: sample key '" << kv.first.name()
                     << "' is not a state -- ignored\n";
    }
    if( n_resolved ){
      samples.push_back( s );
      labels.push_back( "user#" + std::to_string( n_user++ ) );
    }
  }
  if( n_user == 0 ){
    static double const facs[] = { 0.75, 0.90, 1.10, 1.25 };
    for( double f : facs ){
      std::vector<double> s = state_ref;
      for( auto& x : s ) x *= f;
      samples.push_back( s );
      std::ostringstream os; os << "scale" << std::setprecision(2) << f;
      labels.push_back( os.str() );
    }
  }
  // A sample whose state values are not finite cannot be classified on: +/-inf or NaN in a coefficient makes
  // the symbol degenerate for a reason that has nothing to do with the model.  Drop it and say so.
  {
    std::vector<std::vector<double>> keep;  std::vector<std::string> keep_lab;
    for( size_t s = 0; s < samples.size(); ++s ){
      size_t bad = 0;
      for( double x : samples[s] ) if( !std::isfinite( x ) ) ++bad;
      if( !bad ){ keep.push_back( samples[s] ); keep_lab.push_back( labels[s] ); continue; }
      std::cerr << "  DROPPED sample '" << labels[s] << "': " << bad
                << " non-finite state value(s); not classified on\n";
    }
    samples.swap( keep );  labels.swap( keep_lab );
    if( samples.empty() ){
      std::cerr << "  NOTE: every sample was dropped as non-finite -- the classification cannot be checked"
                   " against any reference point.\n";
      _refRobustOk = false;
      return false;
    }
  }

  size_t const nS = samples.size();

  // A multiplicative sweep of an all-zero reference is degenerate (every sample
  // equals the base): flag it so a vacuous PASS is not mistaken for coverage.
  bool base_nonzero = false;
  for( double x : state_ref ) if( std::fabs(x) > 1e-300 ){ base_nonzero = true; break; }
  bool const vacuous_sweep = ( n_user == 0 && !base_nonzero );

  std::cerr << "OCFESLV::reference_robustness ** sampling principal symbol at "
            << nS << " state(s) (1 base + " << (nS-1)
            << (n_user ? " user-supplied envelope" : " default LOCAL magnitude sweep")
            << ")\n";
  if( vacuous_sweep )
    std::cerr << "  NOTE: reference state is all-zero -> the default magnitude sweep is"
              << " VACUOUS (samples == base).\n"
              << "        Supply envelope states via add_classification_reference_sample()"
              << " for genuine coverage.\n";

  auto rank_of = []( arma::mat const& M ) -> arma::uword {
    if( M.n_rows == 0 || M.n_cols == 0 ) return 0;
    arma::vec s;
    if( !arma::svd( s, M ) || s.is_empty() ) return 0;
    double const tol = 1e-8 * std::max( 1.0, s(0) );
    arma::uword r = 0;
    for( arma::uword i = 0; i < s.n_elem; ++i ) if( s(i) > tol ) ++r;
    return r;
  };
  // Max principal angle (deg) between two equal-dimension column spaces; -1 if
  // the spaces are empty or of unequal dimension (not comparable).
  auto max_angle_deg = []( arma::mat const& A, arma::mat const& B ) -> double {
    arma::mat Qa = arma::orth( A ), Qb = arma::orth( B );
    if( Qa.n_cols == 0 || Qb.n_cols == 0 || Qa.n_cols != Qb.n_cols ) return -1.0;
    arma::vec sv;
    if( !arma::svd( sv, Qa.t() * Qb ) || sv.is_empty() ) return -1.0;
    double c = std::min( 1.0, std::max( -1.0, sv.min() ) );
    return std::acos( c ) * 180.0 / 3.14159265358979323846;
  };

  bool all_ok = true;

  for( auto const& bcpair : _blockClassification ){
    int        const  bid      = bcpair.first;
    t_Classify const& base_cls = bcpair.second;
    auto itsym = _blockSymbol.find( bid );
    auto itfd  = _blockFaceData.find( bid );
    if( itsym == _blockSymbol.end() || itfd == _blockFaceData.end() ) continue;
    t_Symbol const& sym       = itsym->second;
    auto     const& base_face = itfd->second;

    // The DESCRIPTOR->PARABOLIC upgrade is sample-independent: the structural
    // index is reference-free, so it is identical at every sample.  Reuse the
    // verdict the block loop already stored on the base classification and
    // apply it in eff_type to BOTH the stored base (already upgraded) and each
    // re-classified sample, so the comparison is apples-to-apples.
    bool const struct_parab = base_cls.structural_parabolic;
    bool const struct_dae   = base_cls.structural_dae;
    auto eff_type = [&]( t_Classify const& c ) -> EqnType {
      if( c.type == DESCRIPTOR   && time_dom && struct_parab ) return PARABOLIC;
      if( c.type == UNDETERMINED && time_dom && struct_dae )   return DIFFERENTIAL_ALGEBRAIC;
      return c.type;
    };

    EqnType const base_type = eff_type( base_cls );
    bool   block_type_ok = true, block_sign_ok = true;
    double worst_angle   = 0.0;
    std::ostringstream type_detail, sign_detail;

    for( size_t s = 1; s < nS; ++s ){               // sample 0 == base (cached)
      auto Ai = _eval_symbol( sym, samples[s], cst_ref, dom_ref, input_ref );
      t_Classify cls = classify( Ai, sym.vDom, time_dom, n_sample, imag_tol,
                                 !sym.vAlgEqn.empty() );
      EqnType const t = eff_type( cls );
      if( t != base_type ){
        block_type_ok = false;  all_ok = false;
        type_detail << "      ** type CHANGES at " << labels[s] << ": "
                    << pde_type_name( base_type ) << " -> " << pde_type_name( t ) << "\n";
      }
      if( base_cls.evolution_hyperbolic ){
        auto face = _build_face_data( sym, cls );
        for( auto const& bf : base_face ){
          t_FaceData const* sf = nullptr;
          for( auto const& f : face ) if( f.dom_idx == bf.dom_idx ){ sf = &f; break; }
          if( !sf ) continue;
          arma::uword const bIn = rank_of( bf.Aplus ), bOut = rank_of( bf.Aminus );
          arma::uword const sIn = rank_of( sf->Aplus ), sOut = rank_of( sf->Aminus );
          // The incoming/outgoing COUNT alone is insufficient: in a +-eigenvalue
          // system a sign flip SWAPS which characteristic is incoming while the
          // count is unchanged (the PSA flow-reversal case).  So GATE 2 also tests
          // the incoming SUBSPACE identity, threshold-free: a flip makes the
          // sample's right-going subspace (Aplus) align BETTER with the base's
          // LEFT-going subspace (Aminus) than with the base's right-going one.
          bool const count_change = ( sIn != bIn || sOut != bOut );
          double const a_keep = max_angle_deg( bf.Aplus,  sf->Aplus ); // base+ vs sample+
          double const a_swap = max_angle_deg( bf.Aminus, sf->Aplus ); // base- vs sample+
          bool swap = false;
          if( a_keep >= 0.0 && a_swap >= 0.0 ) swap = ( a_swap + 1.0 < a_keep );
          if( a_keep > worst_angle ) worst_angle = a_keep;   // residual drift (>=0)
          if( count_change || swap ){
            block_sign_ok = false;  all_ok = false;
            sign_detail << "      ** face " << ( bf.pdom_var ? bf.pdom_var->name() : "?" )
                        << " in/out " << bIn << "/" << bOut << " -> " << sIn << "/" << sOut
                        << " at " << labels[s];
            if( count_change ) sign_detail << " (count change)";
            if( swap ) sign_detail << " (Vin swap: sample-incoming aligns with base-"
                                   << "outgoing, a_keep=" << std::fixed << std::setprecision(1)
                                   << a_keep << " > a_swap=" << a_swap << " deg)"
                                   << std::defaultfloat;
            sign_detail << "\n";
          }
        }
      }
    }

    std::cerr << "  block " << bid << " [base " << pde_type_name( base_type ) << "]\n"
              << "    GATE 1 type      : "
              << ( block_type_ok ? "invariant  [PASS]" : "CHANGES  [FAIL]" ) << "\n";
    if( !block_type_ok ) std::cerr << type_detail.str();
    if( base_cls.evolution_hyperbolic ){
      std::cerr << "    GATE 2 Vin split : "
                << ( block_sign_ok ? "invariant  [PASS]" : "SIGN FLIP / SWAP  [FAIL]" ) << "\n";
      if( !block_sign_ok ) std::cerr << sign_detail.str();
      std::cerr << "    advisory eigvec  : max incoming-subspace drift "
                << std::fixed << std::setprecision(2)
                << ( worst_angle < 0.0 ? 0.0 : worst_angle ) << " deg"
                << ( !block_sign_ok      ? "  [see GATE 2 swap]"
                   : worst_angle > 5.0   ? "  [WARN >5 deg rotation]" : "  [ok]" )
                << std::defaultfloat << "\n";
    }
  }

  _refRobustOk = all_ok;
  std::cerr << "OCFESLV::reference_robustness ** verdict: "
            << ( all_ok
                 ? "PASS (type + characteristic split invariant across samples)"
                 : "FAIL (structural classification is reference-dependent -- see gates above)" )
            << "\n";
  return all_ok;
}

// ======================================================================
// OCFESLV::_build_face_data
// ======================================================================
inline std::vector<FFModel::t_FaceData>
FFModel::_build_face_data
( t_Symbol const& sym, t_Classify const& cls )
const
{
  std::vector<t_FaceData> face_data;
  size_t const nDom   = sym.vDom.size();
  size_t const nState = sym.vState.size();
  face_data.resize( nDom );

  for( size_t i = 0; i < nDom; ++i ){
    face_data[i].pdom_var = &sym.vDom[i];
    face_data[i].dom_idx  = i;
    face_data[i].Aplus    = arma::zeros<arma::mat>( (arma::uword)nState,
                                                    (arma::uword)nState );
    face_data[i].Aminus   = arma::zeros<arma::mat>( (arma::uword)nState,
                                                    (arma::uword)nState );
  }

  if( cls.Ai.empty() || !cls.evolution_hyperbolic ) return face_data;
  if( cls.evolution_dom_idx >= nDom || cls.evolution_dom_idx >= cls.Ai.size() )
    return face_data;

  arma::mat const& Ae = cls.Ai[ cls.evolution_dom_idx ];
  if( Ae.n_rows != nState || Ae.n_cols != nState ) return face_data;

  for( size_t i = 0; i < nDom; ++i ){
    if( i == cls.evolution_dom_idx || i >= cls.Ai.size() ) continue;

    arma::mat const& An = cls.Ai[i];
    if( An.n_rows != nState || An.n_cols != nState ) continue;

    // Split the normal operator through the evolution-normal pencil:
    //     B_n = A_e^{-1} A_n .
    // The reconstructed matrices A_e B_n^+ and A_e B_n^- retain the
    // physical scaling of the original normal coefficient A_n.  For the
    // scalar advection equation u_t + c u_x = 0, this gives A^+ = c, A^- = 0
    // when c > 0, preserving the previous residual scaling.
    arma::mat Bn = arma::solve( Ae, An );

    arma::cx_vec eigval;
    arma::cx_mat eigvec;
    bool ok = arma::eig_gen( eigval, eigvec, Bn );
    if( !ok || eigvec.n_rows != nState || eigvec.n_cols != nState ) continue;

    arma::cx_mat Lplus ( (arma::uword)nState, (arma::uword)nState, arma::fill::zeros );
    arma::cx_mat Lminus( (arma::uword)nState, (arma::uword)nState, arma::fill::zeros );
    for( arma::uword k = 0; k < (arma::uword)nState; ++k ){
      double const re = eigval(k).real();
      if( re > 0. ) Lplus (k,k) = re;
      else          Lminus(k,k) = re;
    }

    arma::cx_mat const Rinv = arma::inv( eigvec );
    face_data[i].Aplus  = arma::real( Ae * eigvec * Lplus  * Rinv );
    face_data[i].Aminus = arma::real( Ae * eigvec * Lminus * Rinv );
  }

  return face_data;
}

inline void
FFModel::_build_face_data
()
{
  _face_data = _build_face_data( _symbol, _classification );
}

inline std::vector<FFModel::t_FaceData> const&
FFModel::_face_data_for_block
( int block_id )
const
{
  auto it = _blockFaceData.find( block_id );
  return it != _blockFaceData.end() ? it->second : _face_data;
}

inline FFModel::t_FaceData const*
FFModel::_face_data_for_dom
( std::vector<t_FaceData> const& fdata, FFVar const* pvar )
{
  if( !pvar ) return nullptr;
  for( auto const& fd : fdata )
    if( fd.pdom_var && fd.pdom_var->id() == pvar->id() ) return &fd;
  return nullptr;
}

// rev152c: "Does residual row_id apply a Partial in direction dom?"  The scan the
// legacy decision layer performed inline (FFPartial ops in the residual subgraph, read
// for their differentiation multi-index), done once per setup for every row x direction.
inline std::map<std::pair<size_t,size_t>,bool> const&
FFModel::_eqn_differentiates_in() const
{
  if( _eqnDiffInReady ) return _eqnDiffIn;
  _eqnDiffIn.clear();
  if( _dag ){
    for( size_t row_id = 0; row_id < _mEqn.size(); ++row_id ){
      FFVar const& e = _mEqn[ row_id ].var;
      FFSubgraph sg = _dag->subgraph( 1, &e );
      for( auto const& kd : _mDom ){
        size_t const did = static_cast<size_t>( kd.first.id().second );
        bool found = false;
        for( auto const& op : sg.l_op ){
          if( !op || !op->sameid( typeid(FFPartial) ) ) continue;
          auto const* pop = mc::type_cast<FFPartial const>( op );
          if( !pop ) continue;
          for( auto const& [dv,ord] : pop->Indep().expr )
            if( ord && dv.id() == kd.first.id() ){ found = true; break; }
          if( found ) break;
        }
        _eqnDiffIn[ { row_id, did } ] = found;
      }
    }
  }
  _eqnDiffInReady = true;
  return _eqnDiffIn;
}


inline bool
FFModel::_classify_pde
( std::vector<double> const& state_ref,
  std::vector<double> const& input_ref,
  std::vector<double> const& cst_ref,
  std::vector<double> const& dom_ref,
  FFVar const*               time_dom,
  unsigned const             n_sample,
  double   const             imag_tol )
{
  // Store the effective evolution domain for per-direction IC_AUTO resolution.
  if( time_dom ){
    _evolution_dom_var = *time_dom;
    _evolution_dom_set = true;
    _reresolve_auto_roles( _mEqn, &_evolution_dom_var );    // the evolution domain may only now be known
  }
  else {
    _evolution_dom_user = false;
    _evolution_dom_set  = false;
  }

  _blockSymbol.clear();
  _blockClassification.clear();
  _blockFaceData.clear();

  // The classification is about to be recomputed; a derived solver drops what that invalidates.
  _on_classify_begin();

  // Identify all blocks that contain equations participating in classification.
  std::set<int> block_ids;
  for( auto const& eqn : _mEqn ){
    auto const& opt    = *eqn.opt;
    if( opt.participate_in_classification ) block_ids.insert( opt.block_id );
  }
  if( block_ids.empty() ) block_ids.insert( 0 );

  bool first_valid = true;
  for( int bid : block_ids ){
    t_Symbol sym_b = _principal_symbol( bid );
    arma::mat A0_b;
    auto Ai_b = _eval_symbol( sym_b, state_ref, cst_ref, dom_ref, input_ref, time_dom? nullptr: &A0_b );
    // The Douglis-Nirenberg weighted symbol of an order-reduced block (weights: 1 for a reduced primitive, 0 for
    // the other columns; 0 for a LINK row, 1 for the other rows): a non-LINK row keeps only its first-order
    // entries in the non-primitive columns; a LINK row keeps its first-order entries in the primitive columns
    // and its zeroth-order entries in the others.  Only for a block without an evolution direction.
    std::vector<arma::mat> Ai_dn;  arma::mat E0_dn;  bool have_dn = false;
    if( !time_dom && !Ai_b.empty() && A0_b.n_rows == (arma::uword)sym_b.vEqn.size() ){
      size_t const nE = sym_b.vEqn.size(), nS = sym_b.vState.size();
      std::vector<bool> is_link( nE, false ), reduced_prim( nS, false );
      for( size_t k = 0; k < nE; ++k )
        for( auto const& eqn : _mEqn )
          if( eqn.opt && eqn.var.id() == sym_b.vEqn[k].id() && eqn.opt->role == EqnRole::LINK ){ is_link[k] = true; break; }
      for( size_t k = 0; k < nE; ++k ) if( is_link[k] ){
        have_dn = true;
        for( size_t i = 0; i < Ai_b.size(); ++i )
          for( size_t j = 0; j < nS; ++j ) if( Ai_b[i]( k, j ) != 0. ) reduced_prim[j] = true;
      }
      if( have_dn ){
        Ai_dn = Ai_b;  E0_dn = arma::mat( nE, nS, arma::fill::zeros );
        for( size_t k = 0; k < nE; ++k )
          for( size_t j = 0; j < nS; ++j ){
            if( !is_link[k] && reduced_prim[j] ) for( auto& A : Ai_dn ) A( k, j ) = 0.;   // lower order
            if( is_link[k] && !reduced_prim[j] ) E0_dn( k, j ) = A0_b( k, j );              // principal
          }
      }
    }
    t_Classify cls_b = classify( Ai_b, sym_b.vDom, time_dom, n_sample, imag_tol,
                                 !sym_b.vAlgEqn.empty(), have_dn? &Ai_dn: nullptr, have_dn? &E0_dn: nullptr );

    // ---- (character, index) wiring -------------------------------------- //
    // Reference-free structural index (Pantelides-style incidence matching).
    // Recorded on EVERY block as the orthogonal index axis.  For a singular-
    // evolution descriptor with a real time direction, upgrade to PARABOLIC
    // exactly when the index-1 elimination is closed by a SPATIAL constraint.
    // This replaces the old has_link/has_volume role heuristic, which both
    // (a) MISSED first-order parabolic closures carrying no reduce_order LINK
    // (M4: d_z c + R u pins u with no LINK -> was left DESCRIPTOR), and
    // (b) could upgrade genuine high-index descriptors purely because a LINK
    // and an interior equation coexisted.  parabolic_character now unifies both
    // sources -- the index-1 spatial-constraint closure (M4/M7) AND the
    // spatial-LINK 2nd-order reduction (M5/heat, index 0) -- so the upgrade
    // fires on the character, not on the index value.
    t_IndexResult const ir = _structural_index_analysis( bid );
    cls_b.differential_index   = ir.index;
    cls_b.structural_parabolic = ir.parabolic_character;
    if( cls_b.type == DESCRIPTOR && time_dom && ir.parabolic_character ){
      cls_b.type = PARABOLIC;
      cls_b.parabolic_structure_detected = true;
    }
    // rev161 F1 (knob CRONOS_PARAB_STRUCT_CHAR=1, default 0 = rev160a).
    // parabolic_structure_detected is assigned NOWHERE ELSE in this header: it is a side
    // effect of the DESCRIPTOR -> PARABOLIC upgrade above, so a block that is parabolic in
    // CHARACTER but never passes through DESCRIPTOR cannot acquire it.  Its only consumer is
    // _block_drop_eligible(), which returns true on this flag alone -- hence such a block is
    // drop-ineligible and the redundancy/value-slaved drops are vetoed there.
    // MEASURED (OCFE_STRONGCORNER1 ladder 3, 2026-09-08): type=DIFFERENTIAL_RECTANGULAR with
    // parab_char=y, parab_struct=0, degen=0, evol_hyp=0 -> vetoed "3 of 3 column-dependent
    // flag(s)"; the SAME three claims one rung down (ladder 2, type upgraded to PARABOLIC,
    // parab_struct=1) are dropped as [w-decide] dropped=6 def0=0 ALL-IMPLIED.  Kept, they
    // become explicit tau columns that the assembled IC_STRONG system reports as ORPHANS
    // ("touched by 0 rows -- multiplier enforces nothing"), k' = explicit-tau count.
    // This knob sets the flag from the CHARACTER test instead, leaving cls_b.type untouched:
    // eligibility follows "is this block parabolic in character", which is what the
    // predicate's own comment says it means ("parabolic / parabolic-closure -> a flux
    // continuity is implied by value continuity + the LINK relation -> eligible").
    // NOT a substitute for the claim-local test (F2): this widens eligibility to every
    // parabolic-character block, including classes no measurement has covered, which is why
    // it ships dormant and why the corpus read must list every model whose plan moves.
    {
      static int const kParabChar = 0;   // CRONOS_PARAB_STRUCT_CHAR (retired 2026-10-07, WORKPLAN 3.B batch 2b)
      // rev161a: MEASURED (OCFE_nlpartial, 2026-09-08) that the unguarded form fires on an
      // ELLIPTIC block (parab_char=y, no evolution direction).  Eligibility means "a flux
      // continuity is implied by value continuity + the LINK relation", which presupposes an
      // evolution direction; an elliptic block has none.  So mirror the upgrade branch's
      // guards exactly, dropping only its `type == DESCRIPTOR` test:
      //   =1  character + time_dom   (tightened; the shipped candidate)
      //   =2  character alone        (the unguarded form, kept so the ELLIPTIC widening stays
      //                              reproducible rather than only recorded)
      bool const f1_ok = ( kParabChar == 2 ) || ( kParabChar == 1 && time_dom != nullptr );
      if( kParabChar && f1_ok && !cls_b.parabolic_structure_detected && ir.parabolic_character ){
        cls_b.parabolic_structure_detected = true;
        std::cerr << "OCFESLV::classify ** F1(" << kParabChar << "): blk" << bid
                  << " parab_struct set from CHARACTER (type=" << pde_type_name( cls_b.type )
                  << " unchanged, parab_char=y, time_dom=" << ( time_dom ? "y" : "n" ) << ")"
                  << std::endl;
      }
    }

    // Pure differential-algebraic block (no transverse spatial domain): the
    // algebraic states carry no time-derivative, so the principal symbol is
    // RECTANGULAR (they are not derivative columns) and classify() bailed to
    // UNDETERMINED at the squareness gate -- never reaching its own pure-time
    // ri.singular ? DAE : ODE branch.  The structural index confirms genuine
    // algebraic structure (index != 0) closed by ALGEBRAIC (non-spatial)
    // constraints (no spatial domain present), which is exactly the DAE enum
    // ("singular evolution, no transverse domain").  Name it DAE; the index
    // axis (differential_index 1 / high) already records its reducibility.
    bool block_has_spatial = false;
    for( auto const& dv : sym_b.vDom )
      if( !time_dom || dv.id() != time_dom->id() ){ block_has_spatial = true; break; }
    if( cls_b.type == UNDETERMINED && time_dom && !block_has_spatial
        && ir.index != 0 ){
      cls_b.type           = DIFFERENTIAL_ALGEBRAIC;
      cls_b.structural_dae = true;
    }

    // classify-probe (DISPLAY_LEVEL>=1): report the deciding facts so an
    // UNDETERMINED verdict can be traced to its exact cause -- rectangular
    // principal symbol (squareness-gate bail) vs a square-but-inconclusive
    // eigenanalysis (A_e nonsingular, complex-small eigenvalues) vs a singular-A_e
    // descriptor that failed the parabolic upgrade.  Pure diagnostic, no logic.
    if( options.DISPLAY_LEVEL >= 1 ){
      size_t const sR = sym_b.vEqn.size(), sC = sym_b.vState.size();
      size_t zero_rows = 0;
      if( !Ai_b.empty() ){
        for( arma::uword r = 0; r < Ai_b[0].n_rows; ++r ){
          bool az = true;
          for( auto const& A : Ai_b ){
            for( arma::uword c = 0; c < A.n_cols; ++c )
              if( A(r,c) != 0. ){ az = false; break; }
            if( !az ) break;
          }
          if( az ) ++zero_rows;
        }
      }
      _disp( 3 ) << "FFModel::classify [probe] blk" << bid
                << " type=" << pde_type_name( cls_b.type )
                << " sym=" << sR << "x" << sC
                << ( sR != sC ? " RECTANGULAR" : " square" )
                << " zerorows=" << zero_rows
                << " droppable=" << ( sR > sC && zero_rows >= sR - sC ? "y" : "n" )
                << " Ae_rank=" << cls_b.rank_evolution
                << " Ae_sing=" << ( cls_b.At_singular ? "y" : "n" )
                << " max_imag=" << cls_b.max_imag_eig
                << " idx=" << ir.index
                << " parab_char=" << ( ir.parabolic_character ? "y" : "n" )
                << " struct_dae=" << ( cls_b.structural_dae ? "y" : "n" )
                << "\n";

      // [DM] Dulmage-Mendelsohn diagnostic (RECTANGULAR blocks only): a maximum
      // bipartite matching on the block's PRINCIPAL incidence names the exact
      // dangling columns/rows behind a rectangular (hence UNDETERMINED) symbol.
      // vState columns are DIFFERENTIAL states ONLY (algebraic states are
      // partitioned into vAlgEqn), so an unmatched COLUMN is a differential state
      // whose determining derivative-row is absent from THIS block -- either
      // imported from a SIBLING face-block isomorphic to this one (e.g. {z} vs
      // {z,rd}@rd=LB: a cross-domain coupling a domain-overlap merge would
      // absorb) or genuinely under-determined.  Naming the columns makes that
      // call by inspection; a sibling-block auto-verdict is a follow-up.  An
      // unmatched ROW is an equation with no assignable principal here.
      // Read-only, print-only: no classification logic depends on this.
      if( sR != sC && !Ai_b.empty() ){
        std::vector< std::vector<char> > adj( sR, std::vector<char>( sC, 0 ) );
        for( auto const& A : Ai_b )
          for( arma::uword r = 0; r < A.n_rows && (size_t)r < sR; ++r )
            for( arma::uword c = 0; c < A.n_cols && (size_t)c < sC; ++c )
              if( A(r,c) != 0. ) adj[(size_t)r][(size_t)c] = 1;
        std::vector<int> colRow( sC, -1 ), rowCol( sR, -1 );
        auto aug = [&]( auto&& self, size_t c, std::vector<char>& seen )->bool{
          for( size_t r = 0; r < sR; ++r ){
            if( !adj[r][c] || seen[r] ) continue;
            seen[r] = 1;
            if( rowCol[r] < 0 || self( self, (size_t)rowCol[r], seen ) ){
              rowCol[r] = (int)c; colRow[c] = (int)r; return true;
            }
          }
          return false;
        };
        for( size_t c = 0; c < sC; ++c ){ std::vector<char> seen( sR, 0 ); aug( aug, c, seen ); }
        std::cerr << "OCFESLV::classify [DM] blk" << bid
                  << " defect=" << ( sC > sR ? sC - sR : sR - sC )
                  << " unmatched_cols={";
        { bool f = true; for( size_t c = 0; c < sC; ++c ) if( colRow[c] < 0 ){
            std::cerr << ( f ? "" : ", " ) << sym_b.vState[c].name(); f = false; } }
        std::cerr << "} unmatched_rows={";
        { bool f = true; for( size_t r = 0; r < sR; ++r ) if( rowCol[r] < 0 ){
            std::cerr << ( f ? "" : ", " ) << sym_b.vEqn[r].name(); f = false; } }
        std::cerr << "}\n";
      }
    }

    _blockSymbol[bid]         = sym_b;
    _blockClassification[bid] = cls_b;
    _blockFaceData[bid]       = _build_face_data( _blockSymbol[bid], cls_b );

    if( bid == 0 || first_valid ){
      _symbol         = _blockSymbol[bid];
      _classification = cls_b;
      _face_data      = _build_face_data( _symbol, _classification );
      first_valid     = false;
    }
  }

  // [DM-global] (v2) cross-block IMPORT vs BENIGN classification.  The per-block
  // [DM] pass names the dangling column/row behind a rectangular block but cannot
  // tell WHY it dangles.  A second maximum bipartite matching on the UNION of every
  // block's principal incidence separates the two causes, using a single robust
  // rule -- the DIFFERENT-BLOCK guard:
  //   * IMPORT -- a column unmatched WITHIN its block but matched GLOBALLY to a
  //     determining differential row in ANOTHER block (an isomorphic sibling face,
  //     e.g. {z} vs {z,rd}@rd=LB).  A domain-overlap merge would absorb it.  The
  //     only way a locally-unmatched column acquires a different-block global match
  //     is a genuine cross-block determining relation, so this cannot false-fire.
  //   * BENIGN -- no different-block determining row exists: a genuine algebraic
  //     closure (index-1/2 DAE derivative-aux, e.g. Dz_c) or a displaced order>2
  //     boundary LINK over-row.  Legitimately rectangular; the DOF audit certifies it.
  // Single-block problems have no sibling, so NOTHING can flip (max-matching
  // non-uniqueness re-matches stay same-block and are ignored) -- the pass cannot
  // false-positive a genuine DAE.  Read-only, print-only: no logic consumes it.
  if( options.DISPLAY_LEVEL >= 1 ){
    std::vector<int> rect_blocks;
    for( int bid : block_ids )
      if( _blockSymbol[bid].vEqn.size() != _blockSymbol[bid].vState.size() )
        rect_blocks.push_back( bid );
    if( !rect_blocks.empty() ){
      // ---- union incidence: unique columns by state id, rows tagged by owner block ----
      std::vector<FFVar> gCol;                             // unique global columns
      auto gcol_id = [&]( FFVar const& v )->int{
        for( size_t i = 0; i < gCol.size(); ++i ) if( gCol[i].id() == v.id() ) return (int)i;
        gCol.push_back( v ); return (int)gCol.size() - 1;
      };
      struct GRow { int bid; int lr; std::vector<int> cols; };
      std::vector<GRow>   gRow;
      std::map<int,int>   rowBase;                          // bid -> first global-row index
      for( int bid : block_ids ){
        rowBase[bid] = (int)gRow.size();
        t_Symbol const& s = _blockSymbol[bid];
        for( size_t c = 0; c < s.vState.size(); ++c ) (void)gcol_id( s.vState[c] );  // pre-register ALL columns
        auto Ai = _eval_symbol( s, state_ref, cst_ref, dom_ref, input_ref );          // (incl. numerically all-zero
        for( size_t r = 0; r < s.vEqn.size(); ++r ){                                 //  columns) so the frozen GNC
          GRow gr; gr.bid = bid; gr.lr = (int)r;                                     //  covers every per-block lookup
          for( size_t c = 0; c < s.vState.size(); ++c ){                             //  and gcol_id never grows gCol
            bool nz = false;                                                         //  after gColRow/gRowCol are sized.
            for( auto const& A : Ai )
              if( (arma::uword)r < A.n_rows && (arma::uword)c < A.n_cols && A(r,c) != 0. ){ nz = true; break; }
            if( nz ) gr.cols.push_back( gcol_id( s.vState[c] ) );
          }
          gRow.push_back( std::move( gr ) );
        }
      }
      size_t const GNR = gRow.size(), GNC = gCol.size();
      // ---- global Kuhn matching (row -> col), then the inverse (col -> row) ----
      std::vector<int> gRowCol( GNR, -1 ), gColRow( GNC, -1 );
      auto gaug = [&]( auto&& self, size_t rr, std::vector<char>& seen )->bool{
        for( int c : gRow[rr].cols ){
          if( seen[c] ) continue; seen[c] = 1;
          if( gColRow[c] < 0 || self( self, (size_t)gColRow[c], seen ) ){
            gColRow[c] = (int)rr; gRowCol[rr] = c; return true; }
        }
        return false;
      };
      for( size_t rr = 0; rr < GNR; ++rr ){ std::vector<char> seen( GNC, 0 ); gaug( gaug, rr, seen ); }
      // ---- per rectangular block: local matching, then DIFFERENT-BLOCK diff ----
      for( int bid : rect_blocks ){
        t_Symbol const& s = _blockSymbol[bid];
        size_t const nR = s.vEqn.size(), nC = s.vState.size();
        auto Ai = _eval_symbol( s, state_ref, cst_ref, dom_ref, input_ref );
        std::vector< std::vector<char> > adj( nR, std::vector<char>( nC, 0 ) );
        for( auto const& A : Ai )
          for( arma::uword r = 0; r < A.n_rows && (size_t)r < nR; ++r )
            for( arma::uword c = 0; c < A.n_cols && (size_t)c < nC; ++c )
              if( A(r,c) != 0. ) adj[(size_t)r][(size_t)c] = 1;
        std::vector<int> colRow( nC, -1 ), rowCol( nR, -1 );
        auto laug = [&]( auto&& self, size_t c, std::vector<char>& seen )->bool{
          for( size_t r = 0; r < nR; ++r ){
            if( !adj[r][c] || seen[r] ) continue; seen[r] = 1;
            if( rowCol[r] < 0 || self( self, (size_t)rowCol[r], seen ) ){ rowCol[r] = (int)c; colRow[c] = (int)r; return true; }
          }
          return false;
        };
        for( size_t c = 0; c < nC; ++c ){ std::vector<char> seen( nR, 0 ); laug( laug, c, seen ); }
        std::set<int> localCols; for( size_t c = 0; c < nC; ++c ) localCols.insert( gcol_id( s.vState[c] ) );
        std::cerr << "OCFESLV::classify [DM-global] blk" << bid << " :";
        bool any = false;
        // [DM-elim] (v3) For a benign (algebraically-closed) dangling column, name the
        // differential CONSUMER of d/d(dir) st -- found via the SYMBOLIC vCoeff (literal-zero
        // test, so it sees through continuation-gated coefficients like -epsE*kT*Vg that the
        // numeric Ai zeroes at setup) -- and the SOURCE, a vAlgEqn equation carrying d/d(dir) st
        // (the advection/EOS row the auto-elimination pass would solve for that derivative and
        // substitute).  Purely structural; no DAG mutation.  Consumer+source is exactly the
        // (E, G) pair the differential-elimination pass consumes.
        auto lit_zero = []( FFVar const& v ){ return v.cst() && v.num().val() == 0.; };
        auto elim_note = [&]( size_t c )->std::string{
          if( s.vCoeff.size() != s.vDom.size() ) return "benign:algebraic/DAE-closure";
          FFVar const& st = s.vState[c];
          int cons_row = -1, cons_dir = -1; std::set<int> dirs;
          for( size_t i = 0; i < s.vDom.size(); ++i ){
            if( s.vCoeff[i].size() != nR * nC ) continue;
            for( size_t k = 0; k < nR; ++k )
              if( !lit_zero( s.vCoeff[i][ k * nC + c ] ) ){
                if( cons_row < 0 ){ cons_row = (int)k; cons_dir = (int)i; }
                dirs.insert( (int)i );
              }
          }
          if( cons_row < 0 ) return "benign:algebraic/DAE-closure";      // no differential consumer
          int src_eq = -1;
          for( int i : dirs ){
            for( size_t g = 0; g < s.vAlgEqn.size(); ++g ){
              auto sg = _dag->subgraph( 1, &s.vAlgEqn[g] );
              bool carries = false;
              for( auto const& op : sg.l_op ){
                if( !op->sameid( typeid(FFPartial) ) ) continue;
                auto const* pop = mc::type_cast<FFPartial const>( op );
                for( size_t jj = 0; jj < op->varin.size() && !carries; ++jj ){
                  FFVar const* operand = op->varin[jj];
                  if( !operand || operand->id() != st.id() ) continue;
                  for( auto const& [iv, ord] : pop->Indep().expr ){ (void)ord;
                    if( iv.id() == s.vDom[i].id() ){ carries = true; break; } }
                }
                if( carries ) break;
              }
              if( carries ){ src_eq = (int)g; break; }
            }
            if( src_eq >= 0 ) break;
          }
          std::string note = "algebraic-closure; d/d" + s.vDom[cons_dir].name() + " "
                           + st.name() + " consumed by " + s.vEqn[cons_row].name();
          note += ( src_eq >= 0 )
                ? ( ", source " + s.vAlgEqn[src_eq].name() + " -> ELIMINABLE (differential substitution)" )
                : std::string( ", no in-block source -> genuine DAE closure" );
          return note;
        };
        // locally-unmatched COLUMN: import iff globally matched to a row in ANOTHER block
        for( size_t c = 0; c < nC; ++c ) if( colRow[c] < 0 ){
          any = true;
          int const grow = gColRow[ gcol_id( s.vState[c] ) ];
          if( grow >= 0 && gRow[grow].bid != bid )
            std::cerr << " col " << s.vState[c].name() << "[IMPORT<-"
                      << _blockSymbol[gRow[grow].bid].vEqn[gRow[grow].lr].name()
                      << "@blk" << gRow[grow].bid << "]";
          else
            std::cerr << " col " << s.vState[c].name() << "[" << elim_note( c ) << "]";
        }
        // locally-unmatched ROW: import iff globally matched to a column absent from this block.  rev324: otherwise it
        // is BENIGN only if it is a DUPLICATE -- the same equation (same residual) entered twice on different node sets,
        // as a displaced LINK and its corner restore arm are (or a driver's hand-split LINK, e.g. OCFE_PDE12m).  The
        // symbol cannot see node sets, so each duplicate adds one structurally empty row; at most multiplicity-1 rows
        // per residual are benign.  Any other unmatched row is a genuine surplus equation: reported as UNMATCHED.
        std::map<std::string,size_t> dup_mult, dup_used;
        for( auto const& v : s.vEqn ) ++dup_mult[ v.name() ];
        for( size_t r = 0; r < nR; ++r ) if( rowCol[r] < 0 ){
          any = true;
          int const gc = gRowCol[ rowBase[bid] + (int)r ];
          if( gc >= 0 && !localCols.count( gc ) )
            std::cerr << " row " << s.vEqn[r].name() << "[ROW-IMPORT->" << gCol[gc].name() << "]";
          else
          {
            std::string const rn = s.vEqn[r].name();
            if( dup_used[ rn ] + 1 < dup_mult[ rn ] ){
              ++dup_used[ rn ];
              std::cerr << " row " << rn << "[benign:duplicate-row(same equation on another node set)]";
            }
            else
              std::cerr << " row " << rn << "[UNMATCHED:surplus-row]";
          }
        }
        if( !any ) std::cerr << " (none)";
        std::cerr << "\n";
      }
    }
  }

#ifdef MC__OCFESLV_STRUCTURAL_PROBE
  // Stage-1a structural decomposition probe (print-only, reference-free).
  // Reports per block the time-differentiated (dyn) vs algebraic/auxiliary (alg)
  // states and the algebraic constraints, flagging spatial (order-raising) vs
  // purely algebraic (DAE-index) ones.  Validates algebraic-variable recovery
  // against the PDE20 corpus before any classification logic consumes it.
  for( int bid : block_ids ){
    t_StructuralDecomp const sd = _structural_dae_decomposition( bid );
    auto emit = [&]( char const* tag, std::vector<FFVar> const& vs ){
      std::cerr << "    " << std::left << std::setw(10) << tag << "{";
      for( size_t i = 0; i < vs.size(); ++i )
        std::cerr << ( i ? ", " : "" ) << vs[i].name();
      std::cerr << "}\n";
    };
    std::cerr << "OCFESLV::_structural_probe ** blk" << bid << "\n";
    emit( "dynamic",   sd.dyn );
    emit( "algebraic", sd.alg );
    std::cerr << "    constraint {";
    for( size_t i = 0; i < sd.constraints.size(); ++i )
      std::cerr << ( i ? ", " : "" ) << sd.constraints[i].name()
                << ( sd.constraint_has_spatial[i] ? "[spatial]" : "[algebraic]" );
    std::cerr << "}\n";
    std::cerr << "    shape: " << sd.dyn.size() << " dyn / " << sd.alg.size()
              << " alg / " << sd.constraints.size() << " constraint"
              << ( sd.alg.empty()
                   ? "   -> regular (index 0)"
                   : "   -> algebraic/auxiliary states present; index & character"
                     " resolved by the matching step (next)" )
              << "\n";
    t_IndexResult const ir = _structural_index_analysis( bid );
    std::cerr << "    INDEX: "
              << ( ir.index == 0 ? std::string( "0 (regular / order-reduction only)" )
                 : ir.index == 1 ? std::string( ir.parabolic_character
                                                ? "1  [parabolic character]" : "1" )
                 : ir.index >= 2 ? ( std::to_string( ir.index ) + " (Stage-2: "
                                     + std::to_string( ir.index - 1 ) + " differentiation"
                                     + ( ir.index - 1 == 1 ? "" : "s" ) + ")" )
                 : std::string( "high (>= 2; unresolved by Stage-2 differentiation)" ) );
    if( ir.index != 0 && ir.index != 1 && !ir.unmatched.empty() ){
      std::cerr << ( ir.index >= 2 ? "  witnesses {" : "  unpinned {" );
      for( size_t i = 0; i < ir.unmatched.size(); ++i )
        std::cerr << ( i ? ", " : "" ) << ir.unmatched[i].name();
      std::cerr << "}";
    }
    std::cerr << "\n";
  }
#endif

#ifdef MC__OCFESLV_SYMBOL_BALANCE_PROBE
  // Item 11 Part-3 (boundary characteristic split) -- READ-ONLY diagnostic, the
  // first step toward auto-generated hyperbolic boundary closure.  The cached
  // projectors A+ (right-going) and A- (left-going) already fix how many
  // conditions each domain end requires: at z=LB (outward normal -z) the
  // INCOMING characteristics are the right-going ones (A+); at z=UB they are the
  // left-going ones (A-).  The eventual generator consumes exactly this; here we
  // only REPORT it so the eigenstructure->boundary mapping can be checked
  // against the PDE14 manual oracle and its sign flip before any row is emitted.
  {
    auto rank_of = []( arma::mat const& M ) -> arma::uword {
      if( M.n_rows == 0 || M.n_cols == 0 ) return 0;
      arma::vec s;
      if( !arma::svd( s, M ) || s.is_empty() ) return 0;
      double const tol = 1e-8 * std::max( 1.0, s(0) );
      arma::uword r = 0;
      for( arma::uword i = 0; i < s.n_elem; ++i ) if( s(i) > tol ) ++r;
      return r;
    };
    for( auto const& [bid, cls] : _blockClassification ){
      if( !cls.evolution_hyperbolic ) continue;
      auto const& fdata = _blockFaceData[bid];
      auto const& symb  = _blockSymbol[bid];
      std::cerr << "OCFESLV::_boundary_characteristic_split ** blk" << bid
                << " (EVOL_HYPERBOLIC, " << symb.vState.size()
                << " states): incoming-condition demand from A+/A-\n";
      for( auto const& fd : fdata ){
        if( !fd.pdom_var ) continue;
        arma::uword const rp = rank_of( fd.Aplus );    // right-going
        arma::uword const rm = rank_of( fd.Aminus );   // left-going
        std::cerr << "    dir " << std::left << std::setw(6) << fd.pdom_var->name()
                  << std::right << " A+rank=" << rp << " (right-going)  A-rank="
                  << rm << " (left-going)\n"
                  << "        LB: " << rp << " incoming BC(s) [A+] + " << rm
                  << " outgoing closure [A-]   |   UB: " << rm
                  << " incoming BC(s) [A-] + " << rp << " outgoing closure [A+]\n";
        // Small blocks: print the projectors so the characteristic combination
        // (and its sign flip, A+ swapping w+ <-> w-) is visible against the oracle.
        if( symb.vState.size() <= 6 ){
          auto show = [&]( char const* tag, arma::mat const& P ){
            std::cerr << "        " << tag << ":";
            for( arma::uword r = 0; r < P.n_rows; ++r ){
              std::cerr << " [";
              for( arma::uword cc = 0; cc < P.n_cols; ++cc )
                std::cerr << ' ' << std::showpos << std::fixed
                          << std::setprecision(3) << P(r,cc);
              std::cerr << " ]";
            }
            std::cerr << std::noshowpos << "\n";
          };
          show( "A+", fd.Aplus );
          show( "A-", fd.Aminus );
        }
      }
    }
  }
#endif

  // The classification is determined but not yet published: a derived solver may derive from it now.
  _on_classify_ready( state_ref, input_ref, cst_ref, dom_ref );

  _classified = true;
  ++_classificationSerial;

#ifdef MC__OCFESLV_REFERENCE_ROBUSTNESS_PROBE
  bool const probe_robust = true;   // the historical compile-time opt-in still forces it on
#else
  bool const probe_robust = ( options.CLASSIFY.ROBUST != Options::ROBUST_OFF );
#endif
  if( probe_robust ){
    _probe_reference_robustness( state_ref, input_ref, cst_ref, dom_ref,
                                 time_dom, n_sample, imag_tol );
    if( !_refRobustOk && options.CLASSIFY.ROBUST == Options::ROBUST_STRICT )
      std::cerr << "OCFESLV::setup [classify] ** the classification is REFERENCE-SENSITIVE: it changes across the"
                   " sampled reference points (see the report above).  A claim resting on it is not safe.\n";
  }

  // Published: a derived solver invalidates what the new classification makes stale.
  if( !_on_classify_published() ) return false;

  return true;
}

inline bool
FFModel::_reduce_order
( bool const reuse_aux )
{
  if( !_issetup ) throw Exceptions( Exceptions::SETUP );

  FFPartial OpP;
  FFIntegral OpILoc;   // local instances used to RE-CREATE a reduction around a folded operand
  FFEval     OpEvalLoc; //   (reduce_order cannot fold inside an external op; it rebuilds, like OpP)

  _auxDef.clear();
  _auxVarID.clear();
  _deferredValue.clear();

  auto add_local_state = [&]( FFVar const& var, std::vector<FFVar> const& dom )
    { _add_working_state( var, dom ); };

  auto add_local_input = [&]( FFVar const& var, std::vector<FFVar> const& dom )
    { _add_working_input( var, dom ); };

  auto register_aux_def = [&]( FFVar const& aux, FFVar const& expr,
                               std::vector<FFVar> const& dom,
                               std::set<FFVar,lt_FFVar> const& diff_dom,
                               int block_id,
                               FFVar const& parent = FFVar() )
    { _register_aux_def( aux, expr, dom, diff_dom, block_id, parent ); };

  std::vector<t_Eqn> wEqn;
  wEqn.reserve( _mEqn.size() );
  for( auto const& eqn : _mEqn )
    wEqn.push_back( eqn );

  std::vector<t_Fct> wFct;
  wFct.reserve( _mFct.size() );
  for( auto const& fct : _mFct )
    wFct.push_back( fct );
  
  bool overall_changed = false;
  bool pass_changed;

  // Reduction-auxiliary name disambiguation.  The base name is built from the consumed (integrated /
  // evaluated) directions and the operand's PRIMARY state, so distinct reductions over the same state
  // (INT x dt vs INT x^2 dt -> both "Intt_x") would otherwise collide by NAME.  Append an ordinal to
  // the 2nd+ use of a base name so every auxiliary has a unique, readable label ("Intt_x", "Intt_x#2").
  // (Identity is always by DAG id; this only fixes the display/label -- and any downstream keying that
  // a caller might do by name, e.g. the accumulator carry input.)
  std::map<std::string,int> aux_name_seq;

  // Reusable probe leaves for the non-affine-partial test (REDUCE.NONLINEAR_PARTIALS).  The test
  // substitutes ONE first-order partial node by a fresh leaf and reads its mc::FFDep dependency
  // type; the substituted expression is discarded.  The pool is grown on demand and reused across
  // equations and passes so the DAG gains at most (max partials in any one residual) inert leaves,
  // never registered as states or inputs.
  std::vector<FFVar> lin_probe;
  auto probe_var = [&]( size_t k ) -> FFVar const& {
    while( lin_probe.size() <= k )
      lin_probe.push_back( _dag->add_var( "_roprobe" + std::to_string( lin_probe.size() ) ) );
    return lin_probe[k];
  };

  // Ids of the first-order bare-state partials in `res` that `res` is NOT affine in.
  // Affinity is tested ONE partial at a time -- that is what makes the criterion minimal:
  // in (c-x_t)*u_xi/x_xi only x_xi comes back non-affine (rational), while x_t and u_xi are
  // each affine on their own and are left as ordinary variable-coefficient partials.
  auto nonaffine_state_partials = [&]( FFVar const& res ) -> std::set<size_t> {
    std::set<size_t> out;
    if( !options.REDUCE.NONLINEAR_PARTIALS ) return out;
    FFVar e = res;
    FFSubgraph sg = _dag->subgraph( 1, &e );
    std::vector<FFVar> dnode;
    for( auto const* op : sg.l_op ){
      if( !op || !op->sameid( typeid(FFPartial) ) ) continue;
      auto const* pop2 = mc::type_cast<FFPartial const>( op );
      if( !pop2 || pop2->Indep().tord != 1 ) continue;                 // high order: peel path
      if( op->varin.empty() || op->varout.empty() ) continue;
      if( !op->varin[0] || !op->varout[0] ) continue;
      if( _mVar.find( *op->varin[0] ) == _mVar.end() ) continue;       // bare state operand only
      dnode.push_back( *op->varout[0] );
    }
    for( size_t k = 0; k < dnode.size(); ++k ){
      FFVar const& pr = probe_var( k );
      FFVar sub;
      try {
        sub = _dag->substitute( std::vector<FFVar>{ e },
                                std::vector<FFVar>{ dnode[k] },
                                std::vector<FFVar>{ pr } )[0];
      }
      catch( ... ){ continue; }
      FFSubgraph sg2 = _dag->subgraph( 1, &sub );
      std::vector<FFVar> leaf;
      for( auto const* op : sg2.l_op )
        if( op && op->type == FFOp::VAR && op->varout[0] ) leaf.push_back( *op->varout[0] );
      if( leaf.empty() ) continue;
      // Seed ONLY the probe; every other leaf carries no dependency, so the result's dependency
      // on the probe is exactly the residual's dependency on that one partial.  FFPartial,
      // FFIntegral and FFEval all propagate FFDep (differentiation/integration are linear), so a
      // residual carrying other external ops is handled without special-casing.
      std::vector<FFDep> dleaf( leaf.size() ), ddep( 1 );
      for( size_t i = 0; i < leaf.size(); ++i )
        if( leaf[i].id() == pr.id() ) dleaf[i].indep( (int)pr.id().second );
      try { _dag->eval( std::vector<FFVar>{ sub }, ddep, leaf, dleaf ); }
      catch( ... ){ continue; }
      if( ddep.size() != 1 ) continue;
      auto const d = ddep[0].dep( (int)pr.id().second );
      if( d.first && d.second != FFDep::TYPE::L ) out.insert( dnode[k].id().second );
    }
    return out;
  };

  do {
    pass_changed = false;
    size_t const n_eqn_this_pass = wEqn.size();

    // -- Substitution maps --------------------------------------------- //
    // Main:  dag_out  ->  first-order replacement  (applied to ALL wEqn)
    // Extra: link_node -> auxiliary wk             (applied to EXISTING
    //        wEqn[0..n_eqn_this_pass] only, before linking eqns are added)
    std::map<size_t, std::pair<FFVar,FFVar>> subst_map;
    std::map<size_t, std::pair<FFVar,FFVar>> extra_subst_map;

    // Linking equations deferred until after the extra substitution
    std::vector<t_Eqn> new_link_eqns;
    // Reduction (OpI/OpEval) LINKs.  Their LHS references the reduction node
    // itself (which IS the main-substitution target dag_out), so unlike the OpP
    // LINKs they must be appended AFTER Phase 4 -- otherwise the subst folds
    // w - reduction into w - w and leaves the auxiliary undefined.
    std::vector<t_Eqn> reduction_link_eqns;

    // Ensure a captured reduction's integrand/operand is a bare COLLOCATED STATE, so the march loop can
    // eval_colloc it (eval_colloc resolves registered states/inputs only -- not a raw reduction node or
    // a compound expression such as x*x).  A bare state is returned as-is; otherwise the operand is
    // materialised as an aux state (LINK operand - w = 0, pointwise/in-solve) and that state returned.

    // -- Helpers for first-order reduction --------------------------- //
    // Mixed partials are reduced by peeling off one first-order derivative
    // at a time.  For example,
    //   D_x D_y u  ->  D_y w,       w - D_x u = 0
    // and, on a later pass if needed,
    //   D_x D_y w  ->  D_y z,       z - D_x w = 0.
    // This preserves the full multi-index derivative structure while ensuring
    // that every final FFPartial operation has total order at most one.
    auto derivative_tag = []( t_SMon const& mon ) -> std::string {
      std::string tag = "D";
      for( auto const& [var,ord] : mon.expr ){
        tag += var.name();
        if( ord > 1 ) tag += std::to_string( ord );
      }
      return tag;
    };

    auto first_direction = []( t_SMon const& mon ) -> FFVar {
      return mon.expr.begin()->first;
    };

    auto monomial_minus_one =
      []( t_SMon const& mon, FFVar const& var ) -> t_SMon {
        auto expr = mon.expr;
        auto it = expr.find( var );
        assert( it != expr.end() && it->second > 0 );
        if( --(it->second) == 0 ) expr.erase( it );
        return t_SMon( expr );
      };

    auto defining_partial_multiindex =
      [&]( FFVar const* var, t_SMon& mon ) -> bool {
        auto sg = _dag->subgraph( 1, var );
        for( auto const& sop : sg.l_op ){
          if( !sop->sameid( typeid(FFPartial) ) ) continue;
          auto const* spop = mc::type_cast<FFPartial const>( sop );
          for( size_t j = 0; j < sop->varout.size(); ++j )
            if( sop->varout[j]->id() == var->id() ){
              mon = spop->Indep();
              return true;
            }
        }
        return false;
      };

    auto collect_operand_state_data =
      [&]( FFVar const* operand,
           std::vector<FFVar>& aux_dom_vec,
           std::string& state_name,
           FFVar& parent_state ) {
        auto op_sg = _dag->subgraph( 1, operand );
        std::set<FFVar, lt_FFVar> aux_dom_seen;
        aux_dom_vec.clear();
        state_name.clear();
        parent_state = FFVar();
        for( auto const& sop : op_sg.l_op ){
          if( sop->type != FFOp::VAR ) continue;
          auto it = _mVar.find( *sop->varout[0] );
          if( it == _mVar.end() ){
            // a TIME-VARYING INPUT in the operand also spans its directions (2026-09-28): without them the
            // materialised auxiliary of a pure-input operand, e.g. OpEval( u, t, tau ) or OpIntegral( u, t ), was
            // domain-less and its defining row u(t) - w kept a free direction -- refused by the model check.
            auto ii = _mInp.find( *sop->varout[0] );
            if( ii != _mInp.end() )
              for( auto const& d : ii->second )
                if( aux_dom_seen.insert(d).second ) aux_dom_vec.push_back(d);
            continue;
          }
          for( auto const& d : it->second )
            if( aux_dom_seen.insert(d).second ) aux_dom_vec.push_back(d);
          if( state_name.empty() ){
            state_name = sop->varout[0]->name();
            parent_state = *sop->varout[0];
          }
        }
        if( state_name.empty() )
          state_name = "expr" + std::to_string( operand->id().second );
      };

    // Nesting helper: does the operand's expression contain a REDUCTION node
    // (FFIntegral / FFEval)?  Used to DEFER a reduction whose operand is itself a
    // reduction -- the inner one is extracted first (subgraph is topological) and
    // its Phase-4 substitution folds it to its aux, so on the NEXT pass the outer
    // reduction sees a clean operand and extracts with the correct post-inner
    // domain and a single-reduction LINK.  (A `tainted` operand alone is too
    // coarse: OpI(D_z u, t) has a tainted operand but no nested reduction and must
    // extract now.)  Mirrors FFPartial's multi-pass resolution of a nested operand.
    auto operand_has_reduction =
      [&]( FFVar const* operand ) -> bool {
        auto osg = _dag->subgraph( 1, operand );
        for( auto const* sop : osg.l_op )
          if( sop && ( sop->sameid( typeid( FFIntegral ) )
                    || sop->sameid( typeid( FFEval ) ) ) )
            return true;
        return false;
      };

    // True if @p operand references (transitively) any already-captured reduction INPUT.  Used to
    // REFUSE nesting an evolution reduction inside another: the inner value is a post-solve captured
    // input, so the outer's in-solve materialised operand would read a window-lagged value.  A nested
    // reduction whose inner is SPATIAL is fine -- that inner is an in-solve aux STATE, not a capture --
    // so this keys strictly on captured inputs, letting spatial-in-evolution nesting through.
    auto operand_depends_on_capture =
      [&]( FFVar const* operand ) -> bool {
        if( _deferredValue.empty() ) return false;
        // FULL (type, index) ids: MC++ numbers variables and intermediate nodes separately, so comparing the index
        // alone matched an intermediate node (Z69) to a captured VARIABLE of index 69 -- a false "nested reduction"
        // refusal whose occurrence depended on how many nodes earlier fdiff calls had added to the shared user DAG
        // (ODESLV_NIFTE: second SYMDIFF on one solver, 2026-09-27).
        std::set<std::decay_t<decltype( std::declval<FFVar>().id() )>> capids;
        for( auto const& C : _deferredValue )
          if( C.input.dag() ) capids.insert( C.input.id() );
        if( capids.empty() ) return false;
        if( operand && capids.count( operand->id() ) ) return true;
        auto osg = _dag->subgraph( 1, operand );
        for( auto const* sop : osg.l_op ){
          if( !sop ) continue;
          for( auto const* v : sop->varin  ) if( v && capids.count( v->id() ) ) return true;
          for( auto const* v : sop->varout ) if( v && capids.count( v->id() ) ) return true;
        }
        return false;
      };

    // Loud refusal for a nested evolution-in-evolution reduction (see operand_depends_on_capture).
    auto refuse_nested_capture = [&]() {
      std::cerr << "OCFESLV::setup ** REFUSED: nested evolution-direction reduction -- an inner "
                   "reduction's captured (post-solve) value feeds an outer evolution reduction's "
                   "in-solve operand, which would be window-lagged under marching.  Consume the inner "
                   "reduction in an output/objective, not inside another evolution reduction.\n";
      _setupStatus = SetupStatus::CAPTURE_NESTED_REFUSED;
      throw OCBase::Exceptions( OCBase::Exceptions::DOMAIN );
    };

    // True if the reduction node @p dag_out is consumed NONLINEARLY within output subgraph @p sg --
    // i.e. it (or anything computed linearly from it) feeds a non-affine operation.  Only such outputs
    // need the captured form; a bare/linear evolution-integral output (dag_out feeds only PLUS/MINUS/
    // NEG/SHIFT/SCALE, or is the output root) stays on the proven fctacc/fctsum path, which the marched
    // accumulation, forward/reverse sensitivities, and SOLVE.REUSE all key on via
    // _fct_is_evolution_integral.  Capturing a linear one would drop it out of that path and zero its
    // sensitivity (STEP 4 will extend the sensitivity machinery to captured outputs).
    auto reduction_nonlinear_in_output =
      []( auto const& sg, FFVar const* dag_out ) -> bool {
        // FULL (type, index) ids -- see operand_depends_on_capture: an index-only match could flag a LINEAR
        // output as nonlinear, capture it, and so zero its sensitivity.
        std::set<std::decay_t<decltype( std::declval<FFVar>().id() )>> D; D.insert( dag_out->id() );
        for( auto const* op : sg.l_op ){
          if( !op ) continue;
          bool any = false;
          for( auto const* v : op->varin ) if( v && D.count( v->id() ) ){ any = true; break; }
          if( !any ) continue;
          switch( op->type ){
            case FFOp::PLUS: case FFOp::MINUS: case FFOp::NEG:
            case FFOp::SHIFT: case FFOp::SCALE:                 // affine in the reduction: still linear
              for( auto const* v : op->varout ) if( v ) D.insert( v->id() );
              break;
            default:                                            // TIMES/DIV/INV/IPOW/nonlinear unary
              return true;
          }
        }
        return false;
      };

    // Collect every domain direction differentiated anywhere inside an
    // operand's expression subgraph.  Unlike defining_partial_multiindex(),
    // which only succeeds when `operand` is itself the output node of an
    // FFPartial, this walks the whole subgraph and therefore also handles
    // composite operands such as  a * D_x(u)  arising from
    // variable-coefficient high-order derivatives  D_x( a * D_x(u) ).
    // Used to populate t_AuxDef::diff_dom in the operand_tainted branches,
    // which previously left diff_dom empty for such operands and so caused
    // _block_has_aux_link_in_direction() to miss the reduced direction,
    // dropping the right-element interface trace of the LINK equation and
    // leaving a C1 discontinuity at finite-element interfaces.
    auto collect_operand_diff_dom =
      [&]( FFVar const* operand ) -> std::set<FFVar,lt_FFVar> {
        std::set<FFVar,lt_FFVar> diff_dom;
        auto op_sg = _dag->subgraph( 1, operand );
        for( auto const& sop : op_sg.l_op ){
          if( !sop || !sop->sameid( typeid(FFPartial) ) ) continue;
          auto const* spop = mc::type_cast<FFPartial const>( sop );
          if( !spop ) continue;
          for( auto const& [dv,ord] : spop->Indep().expr )
            if( ord ) diff_dom.insert( dv );
        }
        return diff_dom;
      };

    // Materialise the first-order state-partials INSIDE a (possibly compound) operand as derivative
    // auxiliary STATES, so the operand becomes STATE-BASED before an OpP/OpI/OpEval consumes or
    // differentiates it.  This generalises the FFPartial derivative-aux trick to ANY operand:
    //   OpP(a,z) -> Dp_z_a (a real state), REUSED with RED_FULL via extra_subst_map keyed by
    //   OpP(...).id() (no duplicate).  Then OpI(OpP(a,z)^2,z) -> OpI(Dp^2,z), and
    //   OpP(a*OpP(u,x),x) -> OpP(a*Dp,x) (handled by the chain-rule expansion) -- no raw OpP survives
    //   in any LINK (which would leave the principal-symbol proxy unresolvable).  Higher-order OpPs
    //   are peeled from the inside out across passes.
    // include_self=false skips the operand's OWN bare-OpP output, leaving a bare high-order OpP to
    // the existing FFPartial peel path (so ordinary PDE reduction is undisturbed); OpI/OpEval pass
    // include_self=true so a bare OpP reduction operand (OpEval(OpP(a,z),z,z0)) also becomes Dp.
    // Returns true if any qualifying inner state-partial was found (caller DEFERS: subst folds it).
    // Materialise the first-order state-partials inside an operand as derivative auxes (reusing
    // RED_FULL's via extra_subst_map keyed by OpP(...).id()) and COLLECT the (OpP -> Dp) pairs so the
    // caller can fold the operand EXPRESSION and rebuild the enclosing op.  A derivative-aux LINK
    // OpP(a,z)-Dp goes in new_link_eqns (seen by the proxy pass, so the OpP proxy resolves).
    // include_self=false skips the operand's own bare-OpP output (leaves bare high-order OpP to the
    // peel path).  Returns true if any qualifying inner state-partial was found.
    auto materialize_operand_partials =
      [&]( FFVar const* operand, int const block_id, bool const include_self,
           std::vector<FFVar>& fold_targ, std::vector<FFVar>& fold_repl ) -> bool {
        bool found = false;
        auto osg = _dag->subgraph( 1, operand );
        for( auto const* sop : osg.l_op ){
          if( !sop || !sop->sameid( typeid(FFPartial) ) ) continue;
          auto const* spop = mc::type_cast<FFPartial const>( sop );
          if( !spop || spop->Indep().tord != 1 ) continue;            // higher-order: peeled from inside out
          if( sop->varin.empty() || sop->varout.empty() ) continue;
          FFVar const* pin  = sop->varin[0];
          FFVar const* pout = sop->varout[0];
          if( !include_self && pout->id().second == operand->id().second ) continue;  // leave bare OpP to peel path
          if( _mVar.find( *pin ) == _mVar.end() ) continue;           // OpP operand not a bare state
          found = true;
          FFVar Dp;
          auto it = extra_subst_map.find( pout->id().second );
          if( it != extra_subst_map.end() ){
            Dp = it->second.second;                                    // reuse RED_FULL's / prior mint
          }
          else{
            std::set<FFVar,lt_FFVar> const& pdom_set = _mVar.find( *pin )->second;
          std::vector<FFVar> pdom( pdom_set.begin(), pdom_set.end() );
            std::set<FFVar,lt_FFVar> pdiff;  std::string dtag;
            for( auto const& [dv,ord] : spop->Indep().expr ){ (void)ord; pdiff.insert( dv ); dtag += dv.name(); }
            Dp = _dag->add_var( "Dp" + dtag + "_" + pin->name() );
            add_local_state( Dp, pdom );
            register_aux_def( Dp, *pout, pdom, pdiff, block_id, *pin );
            std::map<FFVar,int,lt_FFVar> ldom;
            for( auto const& d : pdom ) ldom[d] = FFDom::ALL;
            new_link_eqns.push_back( { *pout - Dp, ldom,
              _make_equation_options( _link_eqn_options( block_id, true ) ) } );
            extra_subst_map[ pout->id().second ] = { *pout, Dp };
          }
          fold_targ.push_back( *pout );
          fold_repl.push_back( Dp );
        }
        return found;
      };

    // -- DETECTION ----------------------------------------------------- //
    for( size_t ieqn = 0; ieqn < n_eqn_this_pass; ++ieqn ){
      auto const& eqnvar = wEqn[ieqn].var;
      auto const& eqndom = wEqn[ieqn].dom;
      auto sg = _dag->subgraph( 1, &eqnvar );
      std::set<FFVar const*, lt_FFVar> tainted;
      // Pre-pass: which first-order bare-state partials does THIS residual use non-affinely?
      // Computed up front because the op walk below is topological -- when the FFPartial branch
      // reaches a node, its downstream consumers have not been visited yet.  Empty (and free)
      // unless REDUCE.NONLINEAR_PARTIALS is set.
      std::set<size_t> const nonaffine = nonaffine_state_partials( eqnvar );

      for( auto const& op : sg.l_op ){

        // Propagate taint through standard arithmetic
        if( op->type != FFOp::EXTERN ){
          bool any_t = false;
          for( auto const& vin : op->varin )
            if( tainted.count(vin) ){ any_t = true; break; }
          if( any_t )
            for( auto const& vout : op->varout ) tainted.insert( vout );
          continue;
        }

        // -- REDUCTION EXTRACTION (OpI / OpEval) -------------------------
        // An OpI/OpEval CONSUMES its direction(s).  Uniformly (like high-order
        // OpP) extract every such reduction into an auxiliary state w carrying
        // the REMAINING directions, with a single-reduction defining LINK
        //     w - OpI/OpEval(operand) = 0 ,
        // and substitute w for the reduction in the parent residual.  This makes
        // every reduction linear in its own LINK (so the section-4a g-correction
        // is always exact) and lets mixed / nonlinear-in-reduction residuals
        // reduce to reduction-free parents.  A LINK's OWN defining reduction is
        // left in place -- it is already the single-reduction form, and skipping
        // it terminates the multi-pass loop.
        if( op->sameid( typeid(FFIntegral) ) || op->sameid( typeid(FFEval) ) ){
          if( wEqn[ieqn].opt->role == EqnRole::LINK ){
            for( auto const& vout : op->varout ) tainted.insert( vout );
            continue;
          }

          bool const is_eval = op->sameid( typeid(FFEval) );
          std::set<FFVar,lt_FFVar> consumed;
          if( is_eval ){
            auto const* eop = mc::type_cast<FFEval const>( op );
            for( auto const& [dv,z0] : eop->Coord() ){ (void)z0; consumed.insert( dv ); }
          }
          else{
            auto const* iop = mc::type_cast<FFIntegral const>( op );
            for( auto const& [dv,ord] : iop->Indep().expr ){ (void)ord; consumed.insert( dv ); }
          }

          for( size_t j = 0; j < op->varin.size(); ++j ){
            FFVar const* operand = op->varin[j];
            FFVar const* dag_out = op->varout[j];
            if( subst_map.count( dag_out->id().second ) ){ tainted.insert( dag_out ); continue; }
            if( operand_has_reduction( operand ) ) continue;   // defer nested reduction (resolved next pass)
            // FFIntegral/FFEval have NO deriv() (by design, like FFPartial).  If the operand is otherwise
            // tainted -- it contains a non-differentiable op such as a first-order OpP that the FFPartial
            // branch accepts as-is (e.g. OpEval(OpP(a,z),z,z0)) -- FAD would reach the reduction.  Mirror
            // FFPartial's tainted-operand path: materialise the operand as an aux STATE (whose LINK is an
            // ordinary OpP(state,dir), FAD-safe), substitute it in, and DEFER; next pass the reduction sees
            // the aux (a clean state) and extracts with no deriv ever needed on a reduction operand.
            // A tainted operand carries state-partials (a first-order OpP the FFPartial branch
            // accepts as-is, possibly wrapped in a coefficient / nonlinear function).  Materialise
            // those partials as derivative STATES (include_self=true so a bare OpP operand becomes a
            // state too), so the reduction extracts against a STATE-based operand and FAD never
            // reaches a deriv().  DEFER: the subst folds the operand next pass.
            // Fold a tainted (state-partial-carrying) operand and REBUILD the reduction around it.
            // reduce_order cannot fold inside an external reduction node via subst (subst replaces a
            // DIRECT varin only, not a compound one) -- so, exactly as it rebuilds OpP with a
            // materialised operand, it materialises the inner OpPs (-> Dz), folds the operand
            // EXPRESSION (ordinary arithmetic, subst recurses), and re-creates OpI/OpEval around it.
            FFVar red_operand = *operand;
            FFVar red_out     = *dag_out;
            {
              std::vector<FFVar> ftarg, frepl;
              if( materialize_operand_partials( operand, wEqn[ieqn].opt->block_id, true, ftarg, frepl ) ){
                red_operand = _dag->substitute( std::vector<FFVar>{ *operand }, ftarg, frepl )[0];
                if( is_eval ){
                  auto const* eop = mc::type_cast<FFEval const>( op );
                  red_out = OpEvalLoc( red_operand, eop->Indep(), eop->Coord(), eop->Side() );   // keep the SIDE
                }
                else{
                  auto const* iop = mc::type_cast<FFIntegral const>( op );
                  red_out = OpILoc( red_operand, iop->Indep() );
                }
              }
              else if( tainted.count( operand ) ) continue;   // still tainted, no foldable partials -> defer
            }

            std::vector<FFVar> aux_dom_vec;
            std::string        state_name;
            FFVar              parent_state;
            collect_operand_state_data( &red_operand, aux_dom_vec, state_name, parent_state );

            // Auxiliary domain = operand's distributed directions MINUS the
            // consumed (integrated / point-evaluated) directions.
            std::vector<FFVar> w_dom;
            for( auto const& d : aux_dom_vec )
              if( !consumed.count( d ) ) w_dom.push_back( d );

            std::string cons_tag;
            for( auto const& d : consumed ) cons_tag += d.name();
            std::string const aux_base =
              std::string( is_eval ? "Ev" : "Int" ) + cons_tag + "_" + state_name;
            int const aux_seq = aux_name_seq[ aux_base ]++;
            std::string const aux_name =
              aux_seq ? aux_base + "#" + std::to_string( aux_seq + 1 ) : aux_base;

            // An evolution-direction reduction (definite integral OR interior point) is CAPTURED: the
            // reduction node becomes a captured scalar/profile INPUT (init 0), written POST-SOLVE from
            // the collocated states -- no in-block state, ODE, IC, or transfer.  Generic to monolithic
            // (one post-solve eval) and marching (per-window).  The value feeds outputs/objectives only
            // (in-solve state consumption is the causality refusal).  Spatial reductions keep the
            // in-solve quadrature/point form below.
            if( _is_deferred_reduction( consumed ) ){
              // Both modes read the integrand/operand as a COLLOCATED STATE via eval_colloc:
              //   LATCH: value = R(tau) at t=tau.
              //   ACCUM: value += INT R dt over the window, formed as an explicit element quadrature of
              //          that state (the reduction node itself is not eval_colloc-able).
              FFVar cap = _dag->add_var( aux_name + ( is_eval ? "_L" : "_Q" ) );
              add_local_input( cap, w_dom );
              t_DeferredValue C;
              C.input      = cap;
              C.block_id   = wEqn[ieqn].opt->block_id;
              C.accumulate = !is_eval;
              if( operand_depends_on_capture( &red_operand ) ) refuse_nested_capture();
              C.source     = red_operand;                                  // the expression itself: see t_DeferredValue
              C.source_dom.insert( aux_dom_vec.begin(), aux_dom_vec.end() );
              if( is_eval )
                for( auto const& [dv,z0] : mc::type_cast<FFEval const>( op )->Coord() )
                  if( dv.id().second == _evolution_dom_var.id().second ){
                    C.tau  = z0;
                    C.side = mc::type_cast<FFEval const>( op )->side( dv );
                  }
              _deferredValue.push_back( C );
              subst_map[ dag_out->id().second ] = { *dag_out, cap };
              tainted.insert( dag_out );
              pass_changed = overall_changed = true;
              continue;
            }

            FFVar w = _dag->add_var( aux_name );
            add_local_state( w, w_dom );
            // A reduction auxiliary is NOT a derivative: empty diff_dom (it claims
            // no differentiated interface direction).
            register_aux_def( w, red_out, w_dom, std::set<FFVar,lt_FFVar>(),
                              wEqn[ieqn].opt->block_id, parent_state );

            std::map<FFVar,int,lt_FFVar> link_dom;
            for( auto const& d : w_dom ) link_dom[d] = FFDom::ALL;
            // Defining LINK  w - reduction = 0.  new_link_eqns is appended AFTER
            // the substitution below, so the LINK keeps the real reduction node
            // (the subst does not fold it to w - w).
            reduction_link_eqns.push_back( { w - red_out, link_dom,
              _make_equation_options( _link_eqn_options( wEqn[ieqn].opt->block_id, true ) ) } );

            // Cross-equation reuse is automatic: subst_map is applied to ALL
            // wEqn in Phase 4, so a repeated reduction node folds to the same w.
            subst_map[ dag_out->id().second ] = { *dag_out, w };
            tainted.insert( dag_out );
            pass_changed = overall_changed = true;
          }
          continue;
        }

        if( !op->sameid( typeid(FFPartial) ) ) continue;
        auto const* pop = mc::type_cast<FFPartial const>( op );
        t_SMon const full_indep = pop->Indep();

        for( size_t j = 0; j < op->varin.size(); ++j ){
          FFVar const* operand = op->varin[j];
          FFVar const* dag_out = op->varout[j];

          bool const operand_tainted = tainted.count( operand ) > 0;
          bool const high_order      = full_indep.tord > 1;

          // Nesting: a partial derivative OVER a reduction, OpP(OpI/OpEval(...),dir), is
          // deferred so the inner reduction extracts first; next pass this OpP sees the
          // reduction's aux (a clean state) and reduces normally.  Without this the tainted-
          // operand path below would materialise a redundant duplicate of that aux.
          if( operand_has_reduction( operand ) ) continue;

          // An OpP of a COMPOUND operand carrying inner state-partials -- a conservative flux
          // OpP(D*du/dz,z) or a chain OpP(u*du/dz,z) -- is left to the whole-operand materialisation /
          // chain-rule path below.  That path materialises the WHOLE flux (e.g. D*du/dz) as a single
          // auxiliary, which is exactly what conservative flux continuity across a discontinuous
          // coefficient requires (materialising only the bare du/dz would force du/dz continuous
          // instead of the flux).  The raw OpP in such an aux's LINK resolves via the principal-symbol
          // proxy because an OpP LINK is a DIFFERENTIAL equation -- unlike a nonlocal OpI/OpEval
          // reduction LINK, which cannot expose an inner OpP, and is why ONLY the reduction operand
          // case (OpI/OpEval branch above) needs the inner-OpP fold + rebuild.

          // A first-order derivative of a clean operand.  Accept it AS-IS only
          // when the operand is a BARE STATE: the principal-symbol proxy pass
          // (which requires OpP(state,dir)) and the collocation Jacobian both
          // assume every partial differentiates a state.  A first-order partial
          // of a COMPOUND operand -- a conservative flux d_z(u*c) -- is instead
          // normalised by the chain rule
          //     d_dir f  ->  (df/d dir)|_explicit  +  sum_s (df/ds) d_dir s ,
          // so every surviving partial d_dir s is an ordinary state-partial.
          // The operand is clean (!operand_tainted) so it carries no nested
          // FFPartial/FFIntegral; the FAD forming df/ds and df/d(dir) therefore
          // never differentiates an external node and is safe.  No new
          // collocated state and no LINK are introduced -- d_dir s is evaluated
          // by the spatial differentiation operator exactly like any first-order
          // state partial.  NOTE: this replaces the conservative discrete
          // operator D_z(u*c) with the split form u*(D_z c)+c*(D_z u); they are
          // analytically identical but differ by spectral aliasing.
          if( !operand_tainted && !high_order ){
            if( _mVar.find( *operand ) != _mVar.end() ){
              // A first-order partial of a bare state used NON-AFFINELY -- divided by, or wrapped
              // in a transcendental.  The proxy pass and the collocation Jacobian both assume a
              // partial enters its residual affinely, which is why a raw 1/x_xi previously had to
              // be hand-rolled as a state by the driver.  Materialise it as a derivative auxiliary
              // (mint + defining LINK + extra_subst_map reuse -- the SAME protocol the OpI/OpEval
              // operand path uses), so the nonlinearity closes over a STATE and the raw OpP
              // survives only in its LINK, which is an ordinary differential equation.
              //   Reuse is by DAG node id, so the same partial appearing in several residuals -- or
              // several times in one -- yields exactly ONE auxiliary: this is what removes the
              // multiply-defined flux states in the ALE formulations.
              if( nonaffine.count( dag_out->id().second ) ){
                std::vector<FFVar> ftarg, frepl;
                if( materialize_operand_partials( dag_out, wEqn[ieqn].opt->block_id,
                                                  true, ftarg, frepl ) ){
                  if( options.DISPLAY_LEVEL >= 1 && !frepl.empty() )
                    std::cerr << "OCFESLV::setup ** reduce_order: non-affine first-order partial "
                              << *dag_out << " materialised as " << frepl.back()
                              << " (block " << wEqn[ieqn].opt->block_id << ")\n";
                  tainted.insert( dag_out );
                  pass_changed = overall_changed = true;
                  continue;
                }
              }
              tainted.insert( dag_out );             // operand is a bare state
              continue;
            }
            if( subst_map.count( dag_out->id().second ) ){
              tainted.insert( dag_out );
              continue;
            }

            // ---- conservative-form chain-rule expansion -------------------
            FFVar const dir_var = full_indep.expr.begin()->first; // tord==1

            std::vector<FFVar> fstates;                // states inside operand
            {
              auto fsg = _dag->subgraph( 1, operand );
              std::set<FFVar::pt_idVar> seen;
              for( auto const& sop : fsg.l_op ){
                if( sop->type != FFOp::VAR || sop->varout.empty() || !sop->varout[0] )
                  continue;
                FFVar const* v = sop->varout[0];
                if( _mVar.find( *v ) == _mVar.end() ) continue;   // states only
                if( !seen.insert( v->id() ).second ) continue;
                fstates.push_back( *v );
              }
            }

            // Jacobian of the algebraic operand wrt [dir_var, states...].
            // dir_var first captures explicit coordinate dependence, e.g.
            // d_z(z*u) = u + z*d_z u.
            std::vector<FFVar> wrt;
            wrt.reserve( 1 + fstates.size() );
            wrt.push_back( dir_var );
            for( auto const& s : fstates ) wrt.push_back( s );
            auto dfd = _dag->FAD( std::vector<FFVar>{ *operand }, wrt );  // 1 x (1+nS)

            FFVar expansion = dfd[0];                  // explicit (df/d dir)
            for( size_t s = 0; s < fstates.size(); ++s )
              expansion += dfd[1+s] * OpP( fstates[s], full_indep );

            subst_map[ dag_out->id().second ] = { *dag_out, expansion };
            tainted.insert( dag_out );
            pass_changed = overall_changed = true;
            continue;
          }

          if( subst_map.count( dag_out->id().second ) ){
            tainted.insert( dag_out );
            continue;
          }

          pass_changed = overall_changed = true;

          std::vector<FFVar> aux_dom_vec;
          std::string        state_name;
          FFVar              parent_state;
          collect_operand_state_data( operand, aux_dom_vec, state_name, parent_state );

          // Interior domain map for linking equations
          std::map<FFVar,int,lt_FFVar> link_dom;
          for( auto const& [var,lim] : eqndom ) link_dom[var] = FFDom::ALL;

          if( operand_tainted ){
            // Materialise the already-differentiated operand as an auxiliary
            // state.  This handles nested first-order partials such as
            // D_y(D_x u), and also nested high-order partials; any remaining
            // high total order is reduced in subsequent passes.
            t_SMon operand_indep;
            std::string aux_name;
            if( defining_partial_multiindex( operand, operand_indep ) && operand_indep.tord )
              aux_name = derivative_tag( operand_indep ) + "_" + state_name;
            else
              aux_name = "Daux" + std::to_string( operand->id().second )
                       + "_" + state_name;

            FFVar w = _dag->add_var( aux_name );
            add_local_state( w, aux_dom_vec );
            // Record the differentiated domain directions of the operand.
            // Scanning the whole operand subgraph (rather than requiring the
            // operand to be a bare FFPartial output) ensures composite
            // operands such as  a * D_x(u)  from variable-coefficient
            // high-order derivatives still yield a non-empty diff_dom, so the
            // reduced direction is recognised by
            // _block_has_aux_link_in_direction() and both one-sided LINK
            // interface traces are retained by _keep_eqn_node().
            std::set<FFVar,lt_FFVar> diff_dom =
              collect_operand_diff_dom( operand );
            register_aux_def( w, *operand, aux_dom_vec, diff_dom,
                              wEqn[ieqn].opt->block_id, parent_state );

            // Linking: operand_expr - w = 0.  LINK rows are pointwise
            // identities; they are kept at strong reduced-order interfaces
            // and receive weak SAT through explicit receive_sat metadata.
            new_link_eqns.push_back( { *operand - w, link_dom,
              _make_equation_options( _link_eqn_options( wEqn[ieqn].opt->block_id, true ) ) } );

            // Optionally reuse w wherever the same operand derivative appears
            // in equations already present at the start of this pass.
            // rev306: sharing is unconditional -- see the header note.  One auxiliary per distinct
            // derivative, whichever reduction is selected.
            extra_subst_map[ operand->id().second ] = { *operand, w };

            FFVar repl = full_indep.tord ? OpP( w, full_indep ) : w;
            subst_map[ dag_out->id().second ] = { *dag_out, repl };
          }
          else {
            // Clean high-order or mixed derivative.  Peel off a single
            // first-order derivative in one direction:
            //   D^alpha operand -> D^(alpha-e_d) w,  w = D_d operand.
            // The remaining multi-index alpha-e_d may still have total order
            // greater than one; it will be processed on the next pass.
            FFVar const peel_var = first_direction( full_indep );
            t_SMon const peel_smon( peel_var, 1u );
            t_SMon const rest_smon = monomial_minus_one( full_indep, peel_var );

            FFVar link_node = OpP( *operand, peel_smon );

            // Reuse an auxiliary already created (this pass) for this exact
            // peeled derivative, so a high-order BOUNDARY derivative such as
            // OpP(u,{x,2}) shares the PDE chain's D1 = OpP(u,x) instead of
            // spawning a duplicate aux of the same name.  Both peels occur in
            // the same pass and extra_subst_map is keyed by link_node id.
            FFVar w;
            auto const itreuse = extra_subst_map.find( link_node.id().second );
            if( itreuse != extra_subst_map.end() ){   // rev306: shared whichever reduction is selected
              w = itreuse->second.second;          // existing aux for link_node
            }
            else {
              std::string aux_name = derivative_tag( peel_smon )
                                   + "_" + state_name;
              w = _dag->add_var( aux_name );
              add_local_state( w, aux_dom_vec );
              std::set<FFVar,lt_FFVar> diff_dom{ peel_var };
              register_aux_def( w, link_node, aux_dom_vec, diff_dom,
                                wEqn[ieqn].opt->block_id, parent_state );

              // Linking: D_d(operand) - w = 0.  LINK rows are pointwise
              // identities; they are kept at strong reduced-order interfaces
              // and receive weak SAT through explicit receive_sat metadata.
              new_link_eqns.push_back( { link_node - w, link_dom,
                _make_equation_options( _link_eqn_options( wEqn[ieqn].opt->block_id, true ) ) } );

              // Optionally replace matching first-order derivative nodes in the
              // existing equations, e.g. boundary conditions containing D_d u.
              extra_subst_map[ link_node.id().second ] = { link_node, w };   // rev306
            }

            FFVar repl = rest_smon.tord ? OpP( w, rest_smon ) : w;
            subst_map[ dag_out->id().second ] = { *dag_out, repl };
          }

          tainted.insert( dag_out );
        }
      } // end op loop
    } // end equation loop

    // -- Optional RED_FULL output scan ---------------------------------
    // Scalar outputs do not participate in classification and should not
    // drive the main PDE block symbol.  With RED_FULL/reuse_aux enabled,
    // however, an output containing a high-order derivative can introduce
    // the same first-order auxiliary/linking equation machinery used by
    // equations, so that the output is evaluated from collocated auxiliary
    // states rather than directly from a high-order OCVar derivative.
    {   // rev306: outputs share the equations' auxiliaries whichever reduction is selected
      for( size_t ifct = 0; ifct < wFct.size(); ++ifct ){
        auto const& fctvar = wFct[ifct].var;

        auto sg = _dag->subgraph( 1, &fctvar );
        std::set<FFVar const*, lt_FFVar> tainted;

        for( auto const& op : sg.l_op ){

          // Propagate taint through standard arithmetic.  This mirrors the
          // equation scan: once an integral/high-order derivative has been
          // formed, nested partials above it are materialised as auxiliaries.
          if( op->type != FFOp::EXTERN ){
            bool any_t = false;
            for( auto const& vin : op->varin )
              if( tainted.count(vin) ){ any_t = true; break; }
            if( any_t )
              for( auto const& vout : op->varout ) tainted.insert( vout );
            continue;
          }

          if( op->sameid( typeid(FFIntegral) ) ){
            auto const* iop = mc::type_cast<FFIntegral const>( op );
            std::set<FFVar,lt_FFVar> consumed;
            for( auto const& [dv,ord] : iop->Indep().expr ){ (void)ord; consumed.insert( dv ); }
            bool const evo = _is_deferred_reduction( consumed );
            if( evo ){
              // Evolution-direction OUTPUT integral: CAPTURE it as an accumulating input (cap = INT R dt,
              // summed per window / one full-domain eval), so the enclosing output f(INT R) becomes
              // f(cap) -- ordinary algebra over a captured input, evaluated once at completion (the
              // last-window post-capture eval).  No aux state, so no causality-guard trip; and correct
              // for a NONLINEAR wrapper (unlike the per-window fctacc sum, which gives sum f(INT_win)).
              // The integrand R may be arbitrarily nonlinear in states (INT is additive over windows);
              // only a nested inner reduction defers one pass.  Spatial output integrals stay below.
              for( size_t j = 0; j < op->varin.size(); ++j ){
                FFVar const* operand = op->varin[j];
                FFVar const* dag_out = op->varout[j];
                if( subst_map.count( dag_out->id().second ) ){ tainted.insert( dag_out ); continue; }
                if( operand_has_reduction( operand ) ) continue;   // inner reduction first
                // Stage 2 (unified path): capture EVERY evolution-direction output integral -- linear
                // ones as f=identity over cap.  The capture value + forward/adjoint sensitivity + reuse
                // now match the (retiring) fctacc/fctsum path, so there is a single code path.
                (void) reduction_nonlinear_in_output;   // linearity gate retired; helper kept for now
                std::vector<FFVar> aux_dom_vec; std::string state_name; FFVar parent_state;
                collect_operand_state_data( operand, aux_dom_vec, state_name, parent_state );
                std::vector<FFVar> w_dom;
                for( auto const& d : aux_dom_vec ) if( !consumed.count( d ) ) w_dom.push_back( d );
                std::string cons_tag; for( auto const& d : consumed ) cons_tag += d.name();
                std::string const aux_base = std::string( "Int" ) + cons_tag + "_" + state_name;
                int const aux_seq = aux_name_seq[ aux_base ]++;
                std::string const aux_name =
                  aux_seq ? aux_base + "#" + std::to_string( aux_seq + 1 ) : aux_base;
                FFVar cap = _dag->add_var( aux_name + "_Q" );
                add_local_input( cap, w_dom );
                t_DeferredValue C;
                C.input = cap;
                if( operand_depends_on_capture( operand ) ) refuse_nested_capture();
                C.source = *operand;                                       // the expression itself: see t_DeferredValue
                C.source_dom.insert( aux_dom_vec.begin(), aux_dom_vec.end() );
                C.accumulate = true; C.block_id = 0;
                _deferredValue.push_back( C );
                subst_map[ dag_out->id().second ] = { *dag_out, cap };
                tainted.insert( dag_out );
                pass_changed = overall_changed = true;
              }
              continue;
            }
            // Spatial output integrals stay on the _eval_fct free-accumulate (in-solve quadrature).
            for( auto const& vout : op->varout ) tainted.insert( vout );
            continue;
          }

          // -- OUTPUT POINT REDUCTION (OpEval) -----------------------------
          // add_output( OpEval(c,z,z0) ) : the FFEval consumes z at z0.  Extract
          // it into an auxiliary state + defining LINK exactly as in the equation
          // scan, so the output becomes reduction-free (F = w) and _eval_fct reads
          // the collocated w -- no ref_el/eval_lagrange path needed for it.  This
          // is the output analogue of the equation OpEval extraction; fctrec.point
          // / fctrec.grid (output-node PRESENTATION) are untouched.
          if( op->sameid( typeid(FFEval) ) ){
            auto const* eop = mc::type_cast<FFEval const>( op );
            std::set<FFVar,lt_FFVar> consumed;
            for( auto const& [dv,z0] : eop->Coord() ){ (void)z0; consumed.insert( dv ); }

            // Evolution-direction OUTPUT point reduction: LATCH-capture it as an input (cap = R(tau),
            // written once at tau's window, held), so the enclosing output f(x(tau)) becomes f(cap) --
            // algebra over the captured input, correct under marching for an INTERIOR tau (an aux state
            // w = x(tau) could not be satisfied in windows whose range excludes tau).  Spatial OpEval
            // keeps the in-solve aux-state extraction below.
            //   2026-09-30: an output evaluation consuming the evolution direction is deferred WHATEVER else it consumes
            // (e.g. OpE(OpE(u,t,T),z,z0), folded onto {t,z}); its spatial coordinates stay an in-solve evaluation inside
            // the capture's source (below).  The in-solve extraction made a domain-less aux w = u(T,z0) that every
            // marching window carries but only the window holding T can satisfy -- the marching setup's structural
            // audit failed ("Index mismatch", "not square").  (_is_deferred_reduction, shared with equations and
            // output integrals, is unchanged.)
            bool const evo_point = _evolution_dom_set && _evolution_dom_var.dag() && consumed.count( _evolution_dom_var );
            if( evo_point ){
              double tau = 0.;  int tside = FFDom::MINUS;
              for( auto const& [dv,z0] : eop->Coord() )
                if( dv.id().second == _evolution_dom_var.id().second ){ tau = z0; tside = eop->side( dv ); }
              for( size_t j = 0; j < op->varin.size(); ++j ){
                FFVar const* operand = op->varin[j];
                FFVar const* dag_out = op->varout[j];
                if( subst_map.count( dag_out->id().second ) ){ tainted.insert( dag_out ); continue; }
                if( operand_has_reduction( operand ) ) continue;   // inner reduction first
                std::vector<FFVar> aux_dom_vec; std::string state_name; FFVar parent_state;
                collect_operand_state_data( operand, aux_dom_vec, state_name, parent_state );
                std::vector<FFVar> w_dom;
                for( auto const& d : aux_dom_vec ) if( !consumed.count( d ) ) w_dom.push_back( d );
                std::string cons_tag; for( auto const& d : consumed ) cons_tag += d.name();
                std::string const aux_base = std::string( "Ev" ) + cons_tag + "_" + state_name;
                int const aux_seq = aux_name_seq[ aux_base ]++;
                std::string const aux_name =
                  aux_seq ? aux_base + "#" + std::to_string( aux_seq + 1 ) : aux_base;
                FFVar cap = _dag->add_var( aux_name + "_L" );
                add_local_input( cap, w_dom );
                t_DeferredValue C;
                C.input = cap;
                if( operand_depends_on_capture( operand ) ) refuse_nested_capture();
                // An evaluation that ALSO consumes spatial directions (e.g. OpE(OpE(u,t,T),z,z0), folded into one node on
                // {t,z}): only the EVOLUTION coordinate is captured; the spatial ones stay an in-solve evaluation INSIDE
                // the source, OpE(u,z,z0) -- so the source is a quantity on the remaining directions, which the solver
                // materialises as usual.  Capturing the bare operand dropped z0: a profile u(tau,.) behind a scalar
                // capture, which broke the marching setup's structural audit (2026-09-30).
                FFVar src = *operand;
                std::set<FFVar,lt_FFVar> spatial;
                for( auto const& [dv,z0] : eop->Coord() )
                  if( dv.id().second != _evolution_dom_var.id().second ){
                    FFEval OpS;  src = OpS( src, dv, z0, eop->side( dv ) );  spatial.insert( dv ); }
                C.source = src;                                            // the expression itself: see t_DeferredValue
                for( auto const& d : aux_dom_vec ) if( !spatial.count( d ) ) C.source_dom.insert( d );
                C.accumulate = false; C.tau = tau; C.side = tside; C.block_id = 0;
                _deferredValue.push_back( C );
                subst_map[ dag_out->id().second ] = { *dag_out, cap };
                tainted.insert( dag_out );
                pass_changed = overall_changed = true;
              }
              continue;
            }

            for( size_t j = 0; j < op->varin.size(); ++j ){
              FFVar const* operand = op->varin[j];
              FFVar const* dag_out = op->varout[j];
              if( subst_map.count( dag_out->id().second ) ){ tainted.insert( dag_out ); continue; }
              if( operand_has_reduction( operand ) ) continue;   // defer nested reduction (resolved next pass)
              // FFIntegral/FFEval have NO deriv() (by design, like FFPartial).  If the operand is otherwise
              // tainted -- it contains a non-differentiable op such as a first-order OpP that the FFPartial
              // branch accepts as-is (e.g. OpEval(OpP(a,z),z,z0)) -- FAD would reach the reduction.  Mirror
              // FFPartial's tainted-operand path: materialise the operand as an aux STATE (whose LINK is an
              // ordinary OpP(state,dir), FAD-safe), substitute it in, and DEFER; next pass the reduction sees
              // the aux (a clean state) and extracts with no deriv ever needed on a reduction operand.
              if( tainted.count( operand ) ){
                std::vector<FFVar> mat_dom; std::string mat_name; FFVar mat_parent;
                collect_operand_state_data( operand, mat_dom, mat_name, mat_parent );
                // Empty diff_dom: the materialised operand is CONSUMED by a reduction (point/integral), NOT
  //   differentiated further, so it carries no derivative-continuity claim (unlike FFPartial's
  //   materialisation, whose aux feeds another OpP).  A derivative-continuity claim here is
  //   symbol-ineligible in the degenerate reduction block and breaks IC_STRONG squareness.
  std::set<FFVar,lt_FFVar> mat_diff;
                FFVar w_op = _dag->add_var( "Rmat" + std::to_string( operand->id().second ) + "_" + mat_name );
                add_local_state( w_op, mat_dom );
                register_aux_def( w_op, *operand, mat_dom, mat_diff, 0, mat_parent );
                std::map<FFVar,int,lt_FFVar> mat_link_dom;
                for( auto const& d : mat_dom ) mat_link_dom[d] = FFDom::ALL;
                reduction_link_eqns.push_back( { *operand - w_op, mat_link_dom,
                  _make_equation_options( _link_eqn_options( 0, true ) ) } );
                subst_map[ operand->id().second ] = { *operand, w_op };
                tainted.insert( operand );
                pass_changed = overall_changed = true;
                continue;
              }

              std::vector<FFVar> aux_dom_vec;
              std::string        state_name;
              FFVar              parent_state;
              collect_operand_state_data( operand, aux_dom_vec, state_name, parent_state );

              std::vector<FFVar> w_dom;
              for( auto const& d : aux_dom_vec )
                if( !consumed.count( d ) ) w_dom.push_back( d );

              std::string cons_tag;
              for( auto const& d : consumed ) cons_tag += d.name();
              std::string const aux_name = std::string( "Ev" ) + cons_tag + "_" + state_name;

              FFVar w = _dag->add_var( aux_name );
              add_local_state( w, w_dom );
              register_aux_def( w, *dag_out, w_dom, std::set<FFVar,lt_FFVar>(), 0, parent_state );

              std::map<FFVar,int,lt_FFVar> link_dom;
              for( auto const& d : w_dom ) link_dom[d] = FFDom::ALL;
              reduction_link_eqns.push_back( { w - *dag_out, link_dom,
                _make_equation_options( _link_eqn_options( 0, true ) ) } );

              subst_map[ dag_out->id().second ] = { *dag_out, w };
              tainted.insert( dag_out );
              pass_changed = overall_changed = true;
            }
            continue;
          }

          if( !op->sameid( typeid(FFPartial) ) ) continue;
          auto const* pop = mc::type_cast<FFPartial const>( op );
          t_SMon const full_indep = pop->Indep();

          for( size_t j = 0; j < op->varin.size(); ++j ){
            FFVar const* operand = op->varin[j];
            FFVar const* dag_out = op->varout[j];

            bool const operand_tainted = tainted.count( operand ) > 0;
            bool const high_order      = full_indep.tord > 1;

            if( !operand_tainted && !high_order ){
              tainted.insert( dag_out );
              continue;
            }

            if( subst_map.count( dag_out->id().second ) ){
              tainted.insert( dag_out );
              continue;
            }

            std::vector<FFVar> aux_dom_vec;
            std::string        state_name;
            FFVar              parent_state;
            collect_operand_state_data( operand, aux_dom_vec, state_name, parent_state );

            // Only materialise output-induced auxiliaries when the operand is
            // genuinely distributed.  Pure domain/constant high-order
            // derivatives can be evaluated directly and should not create new
            // collocated states or equations.
            if( aux_dom_vec.empty() ){
              tainted.insert( dag_out );
              continue;
            }

            pass_changed = overall_changed = true;

            std::map<FFVar,int,lt_FFVar> link_dom;
            for( auto const& d : aux_dom_vec ) link_dom[d] = FFDom::ALL;

            if( operand_tainted ){
              t_SMon operand_indep;
              std::string aux_name;
              if( defining_partial_multiindex( operand, operand_indep ) && operand_indep.tord )
                aux_name = derivative_tag( operand_indep ) + "_" + state_name;
              else
                aux_name = "Daux" + std::to_string( operand->id().second )
                         + "_" + state_name;

              FFVar w = _dag->add_var( aux_name );
              add_local_state( w, aux_dom_vec );
              // See the matching equation-scan branch: scan the whole operand
              // subgraph so composite operands (e.g. a * D_x(u)) still yield a
              // non-empty diff_dom and the reduced interface direction is
              // recognised downstream.
              std::set<FFVar,lt_FFVar> diff_dom =
                collect_operand_diff_dom( operand );
              register_aux_def( w, *operand, aux_dom_vec, diff_dom, 0, parent_state );

              new_link_eqns.push_back( { *operand - w, link_dom,
                _make_equation_options( _link_eqn_options( 0, false ) ) } );

              extra_subst_map[ operand->id().second ] = { *operand, w };

              FFVar repl = full_indep.tord ? OpP( w, full_indep ) : w;
              subst_map[ dag_out->id().second ] = { *dag_out, repl };
            }
            else {
              FFVar const peel_var = first_direction( full_indep );
              t_SMon const peel_smon( peel_var, 1u );
              t_SMon const rest_smon = monomial_minus_one( full_indep, peel_var );

              FFVar link_node = OpP( *operand, peel_smon );

              std::string aux_name = derivative_tag( peel_smon )
                                   + "_" + state_name;
              FFVar w = _dag->add_var( aux_name );
              add_local_state( w, aux_dom_vec );
              std::set<FFVar,lt_FFVar> diff_dom{ peel_var };
              register_aux_def( w, link_node, aux_dom_vec, diff_dom, 0, parent_state );

              new_link_eqns.push_back( { link_node - w, link_dom,
                _make_equation_options( _link_eqn_options( 0, false ) ) } );

              extra_subst_map[ link_node.id().second ] = { link_node, w };

              FFVar repl = rest_smon.tord ? OpP( w, rest_smon ) : w;
              subst_map[ dag_out->id().second ] = { *dag_out, repl };
            }

            tainted.insert( dag_out );
          }
        } // end output op loop
      } // end output loop
    }

    if( subst_map.empty() ) break;

    // -- PHASE 2: extra substitution on EXISTING equations only -------- //
    // Apply link_node -> wk substitutions to wEqn[0..n_eqn_this_pass]
    // BEFORE adding linking equations, so the linking equations themselves
    // (which define the relationship) are never rewritten.
    // rev306: THE axis the mode selects.  RED_FULL rewrites rows processed earlier to use the auxiliary;
    // RED_MAIN does not, which is what leaves a balance row its explicit derivative (PDE25: 5.4x better in the
    // exact modes, which read that row's coefficient; G4: the documented RED_FULL/RED_MAIN discrepancy).
    if( reuse_aux && !extra_subst_map.empty() ){
      std::vector<FFVar> extra_targ, extra_repl;
      extra_targ.reserve( extra_subst_map.size() );
      extra_repl.reserve( extra_subst_map.size() );
      for( auto const& [id, tr] : extra_subst_map ){
        extra_targ.push_back( tr.first );
        extra_repl.push_back( tr.second );
      }
      std::vector<FFVar> dep_vec;
      dep_vec.reserve( n_eqn_this_pass + wFct.size() );
      for( size_t i = 0; i < n_eqn_this_pass; ++i )
        dep_vec.push_back( wEqn[i].var );
      for( auto const& fct : wFct )
        dep_vec.push_back( fct.var );
      auto new_dep = _dag->substitute( dep_vec, extra_targ, extra_repl );
      size_t idep = 0;
      for( size_t i = 0; i < n_eqn_this_pass; ++i )
        wEqn[i].var = new_dep[idep++];
      for( size_t i = 0; i < wFct.size(); ++i )
        wFct[i].var = new_dep[idep++];
    }

    // -- PHASE 3: append linking equations ----------------------------- //
    // Added AFTER Phase 2 so they are not touched by the extra substitution.
    // Their left-hand sides contain the original first-order nodes that
    // define the auxiliary variables (e.g. Partial[x](T) - Dx_T = 0).
    for( auto& link_eqn : new_link_eqns )
      wEqn.push_back( std::move( link_eqn ) );

    // -- PHASE 4: main substitution on ALL equations -------------------- //
    // Replace dag_out (high-order nodes) with first-order replacements.
    std::vector<FFVar> targ_vec, repl_vec;
    targ_vec.reserve( subst_map.size() );
    repl_vec.reserve( subst_map.size() );
    for( auto const& [id, tr] : subst_map ){
      targ_vec.push_back( tr.first );
      repl_vec.push_back( tr.second );
    }
    std::vector<FFVar> dep_vec;
    dep_vec.reserve( wEqn.size() + wFct.size() );
    for( auto const& eqn : wEqn )
      dep_vec.push_back( eqn.var );
    for( auto const& fct : wFct )
      dep_vec.push_back( fct.var );
    auto new_dep = _dag->substitute( dep_vec, targ_vec, repl_vec );
    size_t idep = 0;
    for( size_t i = 0; i < wEqn.size(); ++i )
      wEqn[i].var = new_dep[idep++];
    for( size_t i = 0; i < wFct.size(); ++i )
      wFct[i].var = new_dep[idep++];

    // -- PHASE 5: append reduction (OpI/OpEval) LINKs, AFTER the main subst --
    // Their LHS carries the reduction node dag_out, which is a subst target;
    // appending here keeps  w - reduction  intact (Phase 4 would fold it to
    // w - w).  On the next pass these are EqnRole::LINK and are skipped.
    for( auto& link_eqn : reduction_link_eqns )
      wEqn.push_back( std::move( link_eqn ) );

  } while( pass_changed );

  // -- ORDER>2 BOUNDARY CLOSURE (plan §13) ------------------------------- //
  // Retire LINK_j traces at boundary faces that are over-determined by a
  // co-located value+derivative condition on one primitive's chain.  Inert
  // (byte-identical) unless such a co-occurrence exists.
  _displace_overdetermined_boundary_links( wEqn );

  // -- COMMIT ------------------------------------------------------------ //
  _mEqn.clear();
  _mEqn.reserve( wEqn.size() );
  for( auto eqn : wEqn ){
    eqn.opt = _normalised_options( *eqn.opt, eqn.dom );
    _mEqn.push_back( std::move( eqn ) );
  }
  
  _mFct.clear();
  _mFct.reserve( wFct.size() );
  for( auto const& fct : wFct )
    _mFct.push_back( fct );

  return overall_changed;
}

// ======================================================================
// OCFESLV::_displace_overdetermined_boundary_links
//
// Order>2 boundary closure (plan §13).  After reduce_order, a domain-boundary
// face may carry BOTH a value condition on a reduced primitive p AND a
// derivative condition on one of p's auxiliaries D_j (the BC pins D_j, and so
// does LINK_j: D_j = d_dir D_{j-1}).  That node is over-determined.  Retire
// LINK_j's trace at that face by restricting its collocation domain to
// ALL-face (handled downstream by the existing _keep_eqn_node lim switch), and
// normalise the derivative BC to the full face so D_j stays pinned at every
// dropped node (including the initial-time node).
//
// SAFETY.  The value / derivative classification requires the BC to reference
// exactly one state-of-interest (an aux, or a chain-root primitive), so
// multi-state interface conditions (e.g. flux continuity across two blocks)
// never trigger.  Displacement fires only when a value AND a derivative
// condition share the SAME root primitive at the SAME (block,dir,face); an
// order-2 flux wall (aux-defining BC with NO co-located value, e.g. PDE3/PDE4)
// therefore leaves every LINK in place.  Empty result => byte-identical, so the
// order<=2 suite is unaffected.
// ======================================================================
inline void
FFModel::_displace_overdetermined_boundary_links
( std::vector<t_Eqn>& wEqn )
const
{
  if( _auxDef.empty() ) return;                 // no reduction -> nothing to do
  bool const dbg = ( options.DISPLAY_LEVEL >= 2 );   // verbose detection trace (debug only)

  // Aux: primitive root + chain depth, and a (root,dir,depth) lookup so a
  // derivative BC of order j maps to the depth-j auxiliary D_j.
  std::map<size_t,size_t> root_of;              // aux id -> root primitive id
  std::map<size_t,FFVar>  prim_by_id;           // root id -> root var
  struct ADKey { size_t root, dir, depth;
    bool operator<( ADKey const& o ) const {
      if( root != o.root ) return root < o.root;
      if( dir  != o.dir  ) return dir  < o.dir;
      return depth < o.depth; } };
  std::map<ADKey,size_t> aux_lookup;
  std::map<size_t,size_t> depth_of;             // rev322: aux id -> chain depth from its root            // (root,dir,depth) -> aux id
  for( auto const& a : _auxDef ){
    size_t const aid = static_cast<size_t>( a.aux.id().second );
    FFVar const* root = _auxiliary_primitive_root( a.aux );
    if( !root ) continue;
    size_t const rid = static_cast<size_t>( root->id().second );
    root_of[aid] = rid;
    prim_by_id[rid] = *root;
    // chain depth = number of differentiation steps from the primitive root.
    size_t depth = 0;
    FFVar const* cur = &a.aux;
    for( size_t g = 0; g <= _auxDef.size(); ++g ){
      if( !_is_auxiliary_state( *cur ) ) break;
      FFVar const* p = _auxiliary_parent( *cur );
      ++depth;
      if( !p ) break;
      cur = p;
    }
    depth_of[aid] = depth;
    for( auto const& d : a.diff_dom )
      aux_lookup[ ADKey{ rid, static_cast<size_t>(d.id().second), depth } ] = aid;
  }
  if( dbg ){
    std::cout << "OCFESLV::reduce_order [displace-probe] auxes=" << _auxDef.size() << ":";
    for( auto const& a : _auxDef ){
      size_t const aid = static_cast<size_t>( a.aux.id().second );
      std::cout << " aux#" << aid << "(root#"
                << ( root_of.count(aid)? (long)root_of[aid] : -1 ) << ",dir";
      for( auto const& d : a.diff_dom ) std::cout << "#" << (long)d.id().second;
      std::cout << ")";
    }
    std::cout << "\n";
  }

  // States-of-interest (auxes or chain-root primitives) referenced by an eqn.
  auto refs = [&]( FFVar const& eqnvar,
                   std::set<size_t>& aux_hits,
                   std::set<size_t>& prim_hits )
  {
    FFSubgraph sg = _dag->subgraph( 1, &eqnvar );
    auto scan = [&]( FFVar const* v ){
      if( !v ) return;
      // Compare the FULL FFVar identity (type + index).  Operation-result nodes
      // draw indices from a separate sequence that can collide with a state
      // VARIABLE index, so an index-only test spuriously matches intermediate
      // nodes.  This mirrors _equation_has_aux_link_in_direction.
      for( auto const& a : _auxDef )
        if( v->id() == a.aux.id() )
          aux_hits.insert( static_cast<size_t>( a.aux.id().second ) );
      for( auto const& pr : prim_by_id )
        if( v->id() == pr.second.id() )
          prim_hits.insert( pr.first );
    };
    for( auto const& op : sg.l_op ){
      for( auto const* vi : op->varin )  scan( vi );
      for( auto const* vo : op->varout ) scan( vo );
    }
  };

  // Single boundary face (FFDom::LB / FFDom::UB) of a mask in direction dir.
  auto face_in = []( t_EqnDom const& dom, FFVar const& dir ) -> int {
    auto it = dom.find( dir );
    if( it == dom.end() ) return 0;
    if( it->second == FFDom::LB ) return FFDom::LB;
    if( it->second == FFDom::UB ) return FFDom::UB;
    return 0;
  };

  // Highest-order differentiation of each root primitive in each direction,
  // read directly from FFPartial nodes.  This recognises a derivative BC that
  // keeps the partial node  OpP(p,{dir,j})  rather than the substituted aux D_j
  // (the driver's OpP and reduce_order's internal OpP do not share a DAG node,
  // so the substitution never lands on the boundary equation).
  auto partials = [&]( FFVar const& eqnvar,
                       std::map<std::pair<size_t,size_t>,size_t>& ord )
  {
    FFSubgraph sg = _dag->subgraph( 1, &eqnvar );
    for( auto const& sop : sg.l_op ){
      if( !sop || !sop->sameid( typeid(FFPartial) ) ) continue;
      auto const* pop = mc::type_cast<FFPartial const>( sop );
      if( !pop ) continue;
      t_SMon const idx = pop->Indep();
      for( size_t j = 0; j < sop->varin.size(); ++j ){
        FFVar const* operand = sop->varin[j];
        if( !operand ) continue;
        size_t rootid; bool ok = false; size_t base = 0;   // rev322: base = depth of an auxiliary operand
        size_t const oid = static_cast<size_t>( operand->id().second );
        if( prim_by_id.count( oid ) ){ rootid = oid; ok = true; }
        else if( _is_auxiliary_state( *operand ) ){
          FFVar const* r = _auxiliary_primitive_root( *operand );
          if( r ){ rootid = static_cast<size_t>( r->id().second ); ok = true;
                 auto itd = depth_of.find( oid ); if( itd != depth_of.end() ) base = itd->second; }
        }
        if( !ok ) continue;
        for( auto const& [dv, o] : idx.expr ){
          if( o < 1 ) continue;
          auto const key = std::make_pair( rootid, static_cast<size_t>( dv.id().second ) );
          size_t const ov = static_cast<size_t>( o ) + base;   // rev322: OpP(D_k,dir^o) is order k+o on the root
          auto it = ord.find( key );
          if( it == ord.end() || ov > it->second ) ord[key] = ov;
        }
      }
    }
  };

  struct Key { int blk; size_t dir; int face;
    bool operator<( Key const& o ) const {
      if( blk  != o.blk  ) return blk  < o.blk;
      if( dir  != o.dir  ) return dir  < o.dir;
      return face < o.face; } };
  std::map<Key,std::set<size_t>> value_prims;   // (blk,dir,face) -> root prims with value cond
  std::map<Key,std::set<size_t>> deriv_auxes;   // (blk,dir,face) -> auxes with deriv cond
  std::map<size_t,size_t>        derivbc_eqn;   // aux id -> wEqn index of its derivative BC

  for( size_t ie=0; ie<wEqn.size(); ++ie ){
    auto const& eqn = wEqn[ie];
    if( eqn.opt->role != EqnRole::BOUNDARY ) continue;
    std::set<size_t> aux_hits, prim_hits;
    refs( eqn.var, aux_hits, prim_hits );
    std::map<std::pair<size_t,size_t>,size_t> pord;   // (root,dir) -> max partial order
    partials( eqn.var, pord );

    bool has_deriv = false;

    // (1) Derivative conditions from partial structure: OpP(p,{dir,j}) -> D_j.
    for( auto const& [rd, j] : pord ){
      size_t const rootid = rd.first, dirid = rd.second;
      int face = 0;
      for( auto const& [dv,lim] : eqn.dom )
        if( static_cast<size_t>( dv.id().second ) == dirid ){
          face = ( lim==FFDom::LB || lim==FFDom::UB ) ? lim : 0; break;
        }
      if( !face ) continue;
      auto itx = aux_lookup.find( ADKey{ rootid, dirid, j } );
      if( itx == aux_lookup.end() ) continue;        // no matching aux of that order
      deriv_auxes[ Key{ eqn.opt->block_id, dirid, face } ].insert( itx->second );
      derivbc_eqn[ itx->second ] = ie;
      has_deriv = true;
    }

    // (2) Derivative conditions from a direct aux reference (substituted form).
    for( size_t aid : aux_hits ){
      for( auto const& a : _auxDef ){
        if( static_cast<size_t>( a.aux.id().second ) != aid ) continue;
        for( auto const& d : a.diff_dom ){
          int const face = face_in( eqn.dom, d );
          if( !face ) continue;
          deriv_auxes[ Key{ eqn.opt->block_id, static_cast<size_t>(d.id().second), face } ].insert( aid );
          derivbc_eqn[aid] = ie;
          has_deriv = true;
        }
      }
    }

    // (3) Value condition: a referenced root primitive with NO derivative on
    // this equation, so OpP(u,x)-... is never misread as a value condition.
    if( !has_deriv ){
      for( size_t pid : prim_hits )
        for( auto const& [d,lim] : eqn.dom ){
          int const face = ( lim==FFDom::LB || lim==FFDom::UB ) ? lim : 0;
          if( !face ) continue;
          value_prims[ Key{ eqn.opt->block_id, static_cast<size_t>(d.id().second), face } ].insert( pid );
        }
    }

    if( dbg ){
      std::cout << "OCFESLV::reduce_order [displace-probe] BC#" << ie
                << " blk=" << eqn.opt->block_id << " deriv=" << (has_deriv?1:0)
                << " aux_hits={";
      for( size_t h : aux_hits ) std::cout << h << " ";
      std::cout << "} prim_hits={";
      for( size_t h : prim_hits ) std::cout << h << " ";
      std::cout << "} partials={";
      for( auto const& [rd,j] : pord ) std::cout << "p#"<<rd.first<<"/d#"<<rd.second<<"^"<<j<<" ";
      std::cout << "} faces={";
      for( auto const& [d,lim] : eqn.dom )
        std::cout << "#" << (long)d.id().second << ":" << lim << " ";
      std::cout << "}\n";
    }
  }

  // Co-occurrence: value AND derivative on the SAME root primitive at the SAME
  // (blk,dir,face) -> displace that aux's LINK there.
  std::map<size_t,std::set<std::pair<size_t,int>>> drop;   // aux id -> {(dir,face)}
  for( auto const& [k, auxes] : deriv_auxes ){
    auto itv = value_prims.find( k );
    if( itv == value_prims.end() ) continue;
    for( size_t aid : auxes ){
      auto itr = root_of.find( aid );
      if( itr == root_of.end() ) continue;
      if( !itv->second.count( itr->second ) ) continue;   // value on a different primitive
      drop[aid].insert( { k.dir, k.face } );
    }
  }
  if( dbg )
    std::cout << "OCFESLV::reduce_order [displace-probe] value_keys=" << value_prims.size()
              << " deriv_keys=" << deriv_auxes.size() << " drop=" << drop.size() << "\n";
  if( drop.empty() ) return;                    // nothing over-determined

  std::vector<t_Eqn> restore_arms;              // LINK arms to re-add at corners

  for( auto const& [aid, dirfaces] : drop ){
    // Locate LINK_j: the role-LINK equation defining this aux (var == expr-aux).
    FFVar link_target; bool have_target=false;
    for( auto const& a : _auxDef )
      if( static_cast<size_t>( a.aux.id().second ) == aid ){
        link_target = a.expr - a.aux; have_target = true; break;
      }
    if( !have_target ) continue;
    auto itb = derivbc_eqn.find( aid );

    for( auto& eqn : wEqn ){
      if( eqn.opt->role != EqnRole::LINK ) continue;
      if( !( eqn.var.id() == link_target.id() ) ) continue;

      // Re-add LINK_j on the part of the face the derivative BC does NOT cover,
      // so D_j stays tied to the SPECTRAL derivative there (e.g. the initial-
      // time corner) rather than being pinned to the exact boundary value --
      // which would over-constrain the polynomial solution and leave a residual
      // floor.  This reproduces the single-element oracle's closure: LINK_j is
      // present everywhere except where the derivative BC is imposed.
      if( itb != derivbc_eqn.end() ){
        auto const bcdom = wEqn[ itb->second ].dom;   // copy (eqn.dom edited below)
        for( auto const& [dirid, face] : dirfaces )
          for( auto const& [dv, m] : bcdom ){
            if( static_cast<size_t>( dv.id().second ) == dirid ) continue;  // face dir
            int excl = 0;
            if(      m == FFDom::ALL - FFDom::LB ) excl = FFDom::LB;
            else if( m == FFDom::ALL - FFDom::UB ) excl = FFDom::UB;
            if( !excl ) continue;               // derivative BC covers full extent here
            t_Eqn arm = eqn;                    // copy var/role/block
            for( auto& kv : arm.dom ){
              size_t const did = static_cast<size_t>( kv.first.id().second );
              if(      did == dirid )            kv.second = face;
              else if( kv.first.id() == dv.id() ) kv.second = excl;
              else                                kv.second = FFDom::ALL;
            }
            restore_arms.push_back( arm );
          }
      }

      // Drop LINK_j across the full face.
      for( auto const& [dirid, face] : dirfaces )
        for( auto& dv : eqn.dom )
          if( static_cast<size_t>( dv.first.id().second ) == dirid )
            dv.second = dv.second - face;       // ALL -> ALL-face
      break;
    }
    if( options.DISPLAY_LEVEL >= 1 )
      std::cout << "OCFESLV::setup ** order reduction: order>2 boundary-closure LINK displaced for aux id "
                << aid << " at " << dirfaces.size() << " boundary face(s)\n";
  }

  for( auto& arm : restore_arms ) wEqn.push_back( arm );
}

// ======================================================================
// OCFESLV::_auto_diff_eliminate  (Pantelides increment 1; gated by AUTO.DIFF_ELIM)
// ======================================================================
inline int
FFModel::_auto_diff_eliminate()
{
  // Block ids present in the (post-reduce_order) equation set.
  std::set<int> bids;
  for( auto const& e : _mEqn ) bids.insert( e.opt->block_id );

  auto lit_zero  = []( FFVar const& v ){ return v.cst() && v.num().val() == 0.; };
  auto refs_var  = [&]( FFVar const& expr, FFVar const& tgt )->bool{
    auto sg = _dag->subgraph( 1, &expr );
    for( auto const& op : sg.l_op )
      if( op->type == FFOp::VAR && !op->varout.empty() && op->varout[0]
          && op->varout[0]->id() == tgt.id() ) return true;
    return false;
  };
  auto eqn_index = [&]( FFVar const& target )->int{
    for( size_t k = 0; k < _mEqn.size(); ++k )
      if( _mEqn[k].var.id() == target.id() ) return (int)k;
    return -1;
  };

  int n_elim = 0;
  bool changed = true;
  int  guard   = 0;
  // Re-classify after every single elimination so the symbol (and the equation roots
  // it references) is always fresh -- the consumer's t_Eqn.var is rewritten in place.
  while( changed && guard++ < 256 ){
    changed = false;
    for( int bid : bids ){
      t_Symbol sym = _principal_symbol( bid );
      size_t const nR = sym.vEqn.size(), nC = sym.vState.size();
      if( nR == nC || sym.vCoeff.size() != sym.vDom.size() ) continue;   // square / malformed

      for( size_t c = 0; c < nC && !changed; ++c ){
        FFVar const st = sym.vState[c];

        // CONSUMER: a differential row whose SYMBOLIC coefficient of d/d(dir) st is not a
        // literal zero (structural -- sees through continuation-gated coefficients).
        int cons_row = -1, cons_dir = -1;
        for( size_t i = 0; i < sym.vDom.size() && cons_row < 0; ++i ){
          if( sym.vCoeff[i].size() != nR * nC ) continue;
          for( size_t k = 0; k < nR; ++k )
            if( !lit_zero( sym.vCoeff[i][ k * nC + c ] ) ){ cons_row = (int)k; cons_dir = (int)i; break; }
        }
        if( cons_row < 0 ) continue;
        FFVar const dvar = sym.vDom[cons_dir];

        // SOURCE: an algebraic (vAlgEqn) equation carrying FFPartial(st, dir); capture the node.
        FFVar const* p_out = nullptr;
        FFVar        src_eqn;
        for( auto const& G : sym.vAlgEqn ){
          auto sg = _dag->subgraph( 1, &G );
          for( auto const& op : sg.l_op ){
            if( !op->sameid( typeid(FFPartial) ) ) continue;
            auto const* pop = mc::type_cast<FFPartial const>( op );
            for( size_t jj = 0; jj < op->varin.size() && !p_out; ++jj ){
              FFVar const* operand = op->varin[jj];
              if( !operand || operand->id() != st.id() ) continue;
              bool dirmatch = false;
              for( auto const& [iv, ord] : pop->Indep().expr ){ (void)ord;
                if( iv.id() == dvar.id() ){ dirmatch = true; break; } }
              if( dirmatch && !op->varout.empty() && op->varout[0] ) p_out = op->varout[0];
            }
            if( p_out ) break;
          }
          if( p_out ){ src_eqn = G; break; }
        }
        if( !p_out ) continue;   // no algebraic source -> genuine DAE closure, leave alone

        // Minimality / safety guard (matching-INDEPENDENT).  Only eliminate if the consumer
        // retains a principal it OWNS -- i.e. it differentiates at least one OTHER state that
        // has NO algebraic source.  This (a) prevents dissolving a consumer whose every
        // differentiated state is source-backed (it would be left with no principal), and
        // (b) skips a column that has its own DEDICATED differential equation (that equation
        // differentiates only that column -> no owned principal remains -> not eliminated),
        // so a genuinely-provided column is never over-reduced.
        {
          bool retains = false;
          for( size_t cc = 0; cc < nC && !retains; ++cc ){
            if( cc == c ) continue;
            bool diffs = false;
            for( size_t i = 0; i < sym.vDom.size() && !diffs; ++i )
              if( sym.vCoeff[i].size() == nR * nC && !lit_zero( sym.vCoeff[i][ cons_row * nC + cc ] ) )
                diffs = true;
            if( !diffs ) continue;                 // cons_row does not differentiate cc
            FFVar const stcc = sym.vState[cc];
            bool cc_has_src = false;
            for( auto const& G : sym.vAlgEqn ){
              auto sg = _dag->subgraph( 1, &G );
              for( auto const& op : sg.l_op ){
                if( !op->sameid( typeid(FFPartial) ) ) continue;
                for( size_t jj = 0; jj < op->varin.size() && !cc_has_src; ++jj ){
                  FFVar const* operand = op->varin[jj];
                  if( operand && operand->id() == stcc.id() ) cc_has_src = true;
                }
                if( cc_has_src ) break;
              }
              if( cc_has_src ) break;
            }
            if( !cc_has_src ) retains = true;      // cc is a sourceless principal the consumer owns
          }
          if( !retains ) continue;                 // consumer would lose all owned principals -> skip
        }

        // Solve  G = alpha * (d/d(dir) st) + G_rest = 0  ->  d/d(dir) st = -G_rest / alpha.
        FFVar proxy = _dag->add_var( "ade_probe_" + std::to_string( guard ) );
        FFVar Gp    = _dag->substitute( std::vector<FFVar>{ src_eqn },
                                        std::vector<FFVar>{ *p_out },
                                        std::vector<FFVar>{ proxy } )[0];
        FFVar alpha = _dag->FAD( std::vector<FFVar>{ Gp }, std::vector<FFVar>{ proxy } )[0];
        if( lit_zero( alpha ) )        continue;   // not invertible in this derivative
        if( refs_var( alpha, proxy ) ) continue;   // source nonlinear in the derivative -> skip

        FFVar Grest = _dag->substitute( std::vector<FFVar>{ src_eqn },
                                        std::vector<FFVar>{ *p_out },
                                        std::vector<FFVar>{ FFVar( 0. ) } )[0];
        FFVar subst = -Grest / alpha;

        int const ke = eqn_index( sym.vEqn[cons_row] );
        int const kg = eqn_index( src_eqn );
        if( ke < 0 ) continue;

        // (a) substitute the derivative in the consumer's residual expression
        _mEqn[ke].var = _dag->substitute( std::vector<FFVar>{ _mEqn[ke].var },
                                          std::vector<FFVar>{ *p_out },
                                          std::vector<FFVar>{ subst } )[0];
        // (b) relocate: graft the source's transverse-face directions the consumer lacks
        if( kg >= 0 )
          for( auto const& [dv, lim] : _mEqn[kg].dom )
            if( _mEqn[ke].dom.find( dv ) == _mEqn[ke].dom.end() )
              _mEqn[ke].dom[ dv ] = lim;

        if( options.DISPLAY_LEVEL >= 1 )
          std::cerr << "OCFESLV::setup ** auto_diff_elim: eliminated d/d" << dvar.name()
                    << " " << st.name() << " in " << sym.vEqn[cons_row].name()
                    << " via " << src_eqn.name()
                    << " (substitute + relocate onto source face)\n";
        ++n_elim;
        changed = true;   // restart with a fresh _principal_symbol
      }
      if( changed ) break;
    }
  }
  if( options.DISPLAY_LEVEL >= 1 && n_elim )
    std::cerr << "OCFESLV::setup ** auto_diff_elim: " << n_elim
              << " differential-elimination substitution(s) applied\n";
  return n_elim;
}

// ======================================================================
// OCFESLV::_build_reduction_plan   (Stage 2, decision layer)
// Recovered 2026-07-02 from the 2026-06-30 transcripts (lost to a fork).
// Rebound to t_IndexResult::unmatched (the current header's witness list).
// ======================================================================
inline void
FFModel::_build_reduction_plan()
{
  _reductionPlan.clear();
  FFVar const* tdom = _evolution_dom_set ? &_evolution_dom_var : nullptr;

  // local replica of the bare-state collector used by _structural_index_analysis
  auto bare_states = [&]( FFVar const& expr ) -> std::set<FFVar,lt_FFVar> {
    std::set<FFVar,lt_FFVar> bare;
    auto sg = _dag->subgraph( 1, &expr );
    for( auto const& op : sg.l_op ){
      if( op->sameid( typeid(FFPartial) ) || op->sameid( typeid(FFIntegral) ) ) continue;
      for( auto const* in : op->varin )
        if( in && _mVar.find( *in ) != _mVar.end() ) bare.insert( *in );
    }
    return bare;
  };

  for( auto const& [bid, cls] : _blockClassification ){
    if( cls.differential_index < 2 ) continue;
    _reductionPlan.max_index = std::max( _reductionPlan.max_index, cls.differential_index );

    t_IndexResult const ir = _structural_index_analysis( bid );
    if( ir.index < 0 || ir.unmatched.empty() ){ _reductionPlan.resolved = false; continue; }

    // dynamic states of this block = those carrying a d/dt-defining equation
    std::set<FFVar,lt_FFVar> dyn_states;
    auto eqn_is_ode = [&]( FFVar const& var ) -> bool {
      if( !tdom ) return false;
      bool is_ode = false;
      auto sg = _dag->subgraph( 1, &var );
      for( auto const& op : sg.l_op ){
        if( !op->sameid( typeid(FFPartial) ) ) continue;
        auto const* pop = mc::type_cast<FFPartial const>( op );
        for( auto const* operand : op->varin ){
          if( !operand || _mVar.find( *operand ) == _mVar.end() ) continue;
          for( auto const& [iv, ord] : pop->Indep().expr )
            if( iv.id() == tdom->id() ){ dyn_states.insert( *operand ); is_ode = true; }
        }
      }
      return is_ode;
    };
    for( auto const& eqn : _mEqn )
      if( eqn.opt->block_id == bid && eqn.opt->participate_in_classification )
        eqn_is_ode( eqn.var );

    // ODE residual per dynamic state -> dyn_rhs (bare states of the RHS).  One
    // structural differentiation of a constraint exposes the time-derivatives of
    // its dynamic bare states; each such derivative is that state's ODE RHS, so
    // the reachable variables after k differentiations follow the dyn_rhs chain.
    std::map<FFVar,FFVar,lt_FFVar> ode_res;
    for( auto const& eqn : _mEqn ){
      if( eqn.opt->block_id != bid || !eqn.opt->participate_in_classification ) continue;
      auto sg = _dag->subgraph( 1, &eqn.var );
      for( auto const& op : sg.l_op ){
        if( !op->sameid( typeid(FFPartial) ) ) continue;
        auto const* pop = mc::type_cast<FFPartial const>( op );
        for( auto const* operand : op->varin ){
          if( !operand || _mVar.find( *operand ) == _mVar.end() ) continue;
          for( auto const& [iv, ord] : pop->Indep().expr )
            if( iv.id() == tdom->id() ) ode_res.emplace( *operand, eqn.var );
        }
      }
    }
    std::map<FFVar,std::set<FFVar,lt_FFVar>,lt_FFVar> dyn_rhs;
    for( auto const& sr : ode_res ) dyn_rhs[ sr.first ] = bare_states( sr.second );

    // interior algebraic constraints carrying a dynamic state (matching candidates)
    std::vector<FFVar> cons;
    for( auto const& eqn : _mEqn ){
      if( eqn.opt->block_id != bid || !eqn.opt->participate_in_classification ) continue;
      if( eqn.opt->role != EqnRole::INTERIOR ) continue;
      if( eqn_is_ode( eqn.var ) ) continue;
      auto const bs = bare_states( eqn.var );
      bool carries = false;
      for( auto const& s : bs ) if( dyn_states.count( s ) ){ carries = true; break; }
      if( carries ) cons.push_back( eqn.var );
    }

    // witnesses to pin (the high-index algebraic variables)
    std::vector<FFVar> wit;
    for( auto const& w : ir.unmatched ) wit.push_back( w );
    size_t const nW = wit.size(), nC = cons.size();

    // Per-constraint structural reachability (Pantelides differentiation order):
    // expose[ci][wi] = fewest d/dt rounds for constraint ci's chain to reach
    // witness wi.  M2: {x}->y at 1.  M3: {x}->v->L at 2.  Coupled: each
    // constraint reaches its own subset.
    std::vector<std::map<size_t,int>> expose( nC );
    for( size_t ci = 0; ci < nC; ++ci ){
      std::set<FFVar,lt_FFVar> seen, frontier;
      for( auto const& s : bare_states( cons[ci] ) )
        if( ode_res.count( s ) ){ frontier.insert( s ); seen.insert( s ); }
      int order = 0;
      int const bound = (int)ode_res.size() + 2;
      while( !frontier.empty() && order < bound ){
        ++order;
        std::set<FFVar,lt_FFVar> nf;
        for( auto const& s : frontier ){
          auto it = dyn_rhs.find( s );
          if( it == dyn_rhs.end() ) continue;
          for( auto const& r : it->second ){
            for( size_t wi = 0; wi < nW; ++wi )
              if( r.id() == wit[wi].id() && !expose[ci].count( wi ) )
                expose[ci][wi] = order;
            if( ode_res.count( r ) && !seen.count( r ) ){ nf.insert( r ); seen.insert( r ); }
          }
        }
        frontier.swap( nf );
      }
    }

    // Bipartite augmenting-path matching: each witness <- a DISTINCT constraint
    // whose chain exposes it.  This is the Pantelides matching that M2/M3 (single
    // constraint) satisfy trivially and that coupled high-index genuinely needs;
    // it replaces the earlier constraints.front() pick.  A witness with no
    // exposing constraint, or no perfect matching, leaves resolved=false and is
    // the genuinely-coupled subset-differentiation extension (not in the corpus).
    std::vector<int> matchC( nC, -1 ), matchW( nW, -1 );
    std::function<bool(size_t,std::vector<char>&)> aug =
      [&]( size_t wi, std::vector<char>& vis ) -> bool {
        for( size_t ci = 0; ci < nC; ++ci ){
          if( !expose[ci].count( wi ) || vis[ci] ) continue;
          vis[ci] = 1;
          if( matchC[ci] < 0 || aug( (size_t)matchC[ci], vis ) ){
            matchC[ci] = (int)wi; matchW[wi] = (int)ci; return true;
          }
        }
        return false;
      };
    for( size_t wi = 0; wi < nW; ++wi ){
      std::vector<char> vis( nC, 0 );
      if( !aug( wi, vis ) ) _reductionPlan.resolved = false;
    }

    // emit one assignment per matched witness
    for( size_t wi = 0; wi < nW; ++wi ){
      if( matchW[wi] < 0 ) continue;
      size_t const ci = (size_t)matchW[wi];
      t_ReductionPlan::t_Assign a;
      a.block_id   = bid;
      a.constraint = cons[ci];
      a.n_diff     = expose[ci].at( wi );
      a.pinned_var = wit[wi];
      _reductionPlan.assigns.push_back( a );
    }
  }
}

// ======================================================================
// OCFESLV::_reduce_high_index   (Stage 2, execution layer)
// Recovered 2026-07-02 from the 2026-06-30 transcripts (lost to a fork).
// This is the in-place corpus path (M2/M3/M8); the consistent-IC/re-pivot
// PSA extension is reapplied separately (Stage 3b).
// ----------------------------------------------------------------------
// Realise _reductionPlan by differentiating the flagged interior constraints
// in place.  The differentiated constraint is the total time derivative built
// by the FAD chain rule, with the PURE ODE right-hand side substituted for
// each state rate (so no d/dt-operator node leaks into the result -- the
// reduced constraint then re-classifies as an algebraic pin and survives any
// further FAD rounds).  Original constraints are retained at their IC traces.
// The detected differential index in _blockClassification is left untouched.
// ======================================================================
inline bool
FFModel::_reduce_high_index()
{
  if( _reductionPlan.empty() ) return true;
  FFVar const* tdom = _evolution_dom_set ? &_evolution_dom_var : nullptr;
  if( !tdom ){
    if( options.DISPLAY_LEVEL >= 1 )
      std::cerr << "OCFESLV::setup ** high-index reduction: no evolution direction; cannot reduce\n";
    return false;
  }
  FFPartial OpPt;
  bool all_ok = true;
  std::vector<t_Eqn> new_ics;   // synthesized INITIAL traces (consistent-IC + re-pivot)

  auto bare_states = [&]( FFVar const& expr ) -> std::set<FFVar,lt_FFVar> {
    std::set<FFVar,lt_FFVar> bare;
    auto sg = _dag->subgraph( 1, &expr );
    for( auto const& op : sg.l_op ){
      if( op->sameid( typeid(FFPartial) ) || op->sameid( typeid(FFIntegral) ) ) continue;
      for( auto const* in : op->varin )
        if( in && _mVar.find( *in ) != _mVar.end() ) bare.insert( *in );
    }
    return bare;
  };

  for( auto const& a : _reductionPlan.assigns ){
    // ODE residual per dynamic state in this block
    std::map<FFVar,FFVar,lt_FFVar> ode_res;
    for( auto const& eqn : _mEqn ){
      if( eqn.opt->block_id != a.block_id || !eqn.opt->participate_in_classification ) continue;
      auto sg = _dag->subgraph( 1, &eqn.var );
      for( auto const& op : sg.l_op ){
        if( !op->sameid( typeid(FFPartial) ) ) continue;
        auto const* pop = mc::type_cast<FFPartial const>( op );
        for( auto const* operand : op->varin ){
          if( !operand || _mVar.find( *operand ) == _mVar.end() ) continue;
          for( auto const& [iv, ord] : pop->Indep().expr )
            if( iv.id() == tdom->id() ) ode_res.emplace( *operand, eqn.var );
        }
      }
    }
    // pure ODE RHS per state: substitute OpP(s,t) -> 0 in the residual, negate.
    // residual r_s = OpP(s,t) - RHS_s  =>  r_s|_{OpP=0} = -RHS_s.
    std::map<FFVar,FFVar,lt_FFVar> pure_rhs;
    for( auto const& sr : ode_res ){
      FFVar dsdt = OpPt( sr.first, *tdom );
      std::vector<FFVar> dep{ sr.second }, targ{ dsdt }, repl{ FFVar( 0.0 ) };
      std::vector<FFVar> sub = _dag->substitute( dep, targ, repl );
      pure_rhs.emplace( sr.first, - sub[0] );
    }

    // total time derivative via FAD chain rule:
    //   dG/dt = (dG/dt)|_explicit + sum_s (dG/ds) * RHS_s
    auto total_ddt = [&]( FFVar const& G ) -> FFVar {
      auto const bs = bare_states( G );
      std::vector<FFVar> wrt;  wrt.push_back( *tdom );
      std::vector<FFVar> sts;
      for( auto const& s : bs )
        if( pure_rhs.count( s ) ){ wrt.push_back( s ); sts.push_back( s ); }
      auto dfd = _dag->FAD( std::vector<FFVar>{ G }, wrt );   // 1 x (1 + sts)
      FFVar out = dfd[0];                                     // explicit dG/dt
      for( size_t k = 0; k < sts.size(); ++k )
        out += dfd[1+k] * pure_rhs.at( sts[k] );
      return out;
    };

    FFVar G = a.constraint;
    std::vector<FFVar> hidden_levels;                  // levels 0..n_diff-1: they hold at the initial point
    for( int r = 0; r < a.n_diff; ++r ){ hidden_levels.push_back( G ); G = total_ddt( G ); }

    // Locate the INTERIOR equation carrying the consumed constraint.
    t_Eqn* eos_eqn = nullptr;
    for( auto& eqn : _mEqn ){
      if( eqn.opt->block_id != a.block_id ) continue;
      if( eqn.opt->role != EqnRole::INTERIOR ) continue;
      if( eqn.var.id() != a.constraint.id() ) continue;
      eos_eqn = &eqn; break;
    }
    if( !eos_eqn ){
      all_ok = false;
      if( options.DISPLAY_LEVEL >= 1 )
        std::cerr << "OCFESLV::setup ** high-index reduction: interior constraint not found (block "
                  << a.block_id << ")\n";
      continue;
    }

    // The hidden constraints: every level below the one the reduction pins with.  They hold at the initial point,
    // and materialising them there is what lets the model declare the FREE initial data only.
    if( options.REDUCE.HIDDEN_IC && !hidden_levels.empty() ){
      t_EqnDom ic_dom = eos_eqn->dom;
      ic_dom[ *tdom ] = FFDom::LB;
      for( auto const& g : hidden_levels )
        new_ics.push_back( t_Eqn{ g, ic_dom,
                                  _make_equation_options( EqnOptions( EqnRole::INITIAL, a.block_id ) ) } );
      if( options.DISPLAY_LEVEL >= 1 )
        std::cerr << "FFModel::setup ** high-index reduction: " << hidden_levels.size()
                  << " hidden constraint level(s) materialised at the initial point (block " << a.block_id
                  << "); the model declares the free initial data only" << std::endl;
    }

    // The dynamic state pinned by this constraint = the bare state of the
    // constraint that also carries a d/dt-defining ODE residual.
    FFVar dyn; bool have_dyn = false;
    for( auto const& s : bare_states( a.constraint ) )
      if( ode_res.count( s ) ){ dyn = s; have_dyn = true; break; }

    // Does the consumed constraint span the evolution-LB, and is there NO
    // separate INITIAL condition already pinning `dyn` there?  The abstract
    // t-ODE corpus (M2/M3/M8) supplies an explicit IC trace -> have_ic == true
    // -> the in-place replacement (else branch) is used unchanged.  A spatial
    // high-index block whose all-nodes constraint WAS the t=0 pin (e.g. PSA
    // total-continuity + EOS-density, prescribed P) has no such trace: consuming
    // the constraint to pin the witness would leave the dynamic state unpinned
    // at t=0 (singular root).  Synthesize the IC and re-pivot coverage.
    auto evol_mask = [&]( t_EqnDom const& d ) -> int {
      auto it = d.find( *tdom );
      return ( it == d.end() ? (int)FFDom::ALL : it->second );
    };
    int const eos_em = evol_mask( eos_eqn->dom );
    bool const spans_evol_lb = ( eos_em == FFDom::ALL || eos_em == FFDom::LB );
    bool have_ic = false;
    if( have_dyn ){
      for( auto const& eqn : _mEqn ){
        if( eqn.opt->block_id != a.block_id ) continue;
        if( eqn.opt->role != EqnRole::INITIAL ) continue;
        if( !bare_states( eqn.var ).count( dyn ) ) continue;
        int const m = evol_mask( eqn.dom );
        if( m == FFDom::ALL || m == FFDom::LB ){ have_ic = true; break; }
      }
    }

    if( !options.REDUCE.HIDDEN_IC && have_dyn && spans_evol_lb && !have_ic ){
      // Locate the dynamic equation (the ODE residual for `dyn`).
      FFVar const cont_res = ode_res.at( dyn );
      t_Eqn* cont_eqn = nullptr;
      for( auto& eqn : _mEqn ){
        if( eqn.opt->block_id != a.block_id ) continue;
        if( eqn.var.id() != cont_res.id() ) continue;
        cont_eqn = &eqn; break;
      }
      if( !cont_eqn ){
        eos_eqn->var = G;          // no dynamic equation found: fall back to in-place
      }
      else {
        // Consistent-IC + re-pivot (validated by the PSA index-2 reproducer):
        //   - G_reduced inherits the witness coverage (CONT's old domain);
        //   - CONT moves to the dynamic state's coverage = EOS's old domain with
        //     the evolution component reduced to NO_LB (the synthesized IC pins
        //     the evolution-LB);
        //   - the original constraint is retained at the evolution-LB as the IC.
        t_EqnDom const old_eos_dom  = eos_eqn->dom;
        t_EqnDom const old_cont_dom = cont_eqn->dom;
        eos_eqn->var  = G;
        eos_eqn->dom  = old_cont_dom;
        cont_eqn->dom = old_eos_dom;
        cont_eqn->dom[ *tdom ] = FFDom::ALL - FFDom::LB;
        t_EqnDom ic_dom = old_eos_dom;
        ic_dom[ *tdom ] = FFDom::LB;
        new_ics.push_back( t_Eqn{ a.constraint, ic_dom,
                                  _make_equation_options( EqnOptions( EqnRole::INITIAL, a.block_id ) ) } );
        if( options.DISPLAY_LEVEL >= 1 )
          std::cerr << "OCFESLV::setup ** high-index reduction: consistent-IC + re-pivot:"
                       " synthesized INITIAL trace at the evolution-LB and"
                       " re-pivoted CONT/G_reduced coverage (block "
                    << a.block_id << ")\n";
      }
    }
    else {
      // Abstract-corpus / existing-IC path: in-place replacement (unchanged).
      eos_eqn->var = G;
    }
  }

  for( auto& ic : new_ics ) _mEqn.push_back( std::move( ic ) );
  _idxCache.clear();   // _mEqn mutated by the reduction -> cached analysis is stale
  _on_model_changed( ModelChange::DERIVATIVES );
  return all_ok;
}

//-----------------------------------------------------------------------------
// OCFESLV::_generate_hyperbolic_boundary_closure
//
// First-order hyperbolic blocks carry no reduce_order LINK to close their
// boundaries, so a per-state PDE collocation leaves the OUTGOING characteristic
// at each domain end unequationed (under-determined by rank(A-) at LB / rank(A+)
// at UB).  This emits exactly those missing rows.
//
// Per EVOL_HYPERBOLIC block, per spatial (non-evolution) face direction:
//   * LB outgoing = left-going characteristics  -> rowspace of A- (negative part)
//   * UB outgoing = right-going characteristics -> rowspace of A+ (positive part)
// For each basis vector ell of that row space, emit  sum_k ell_k * vEqn[k]
// (the block's full PDE residuals) collocated on {face=end, evolution=ALL-LB,
// other spatial=ALL-LB-UB}, tagged INTERIOR / classify=false / sat=true so it
// closes the count without entering the principal symbol.  The sign of the wave
// speed is carried by A+/A- themselves (they swap when the speed flips), so this
// is sign-blind by construction -- matching the PDE14 manual FORWARD/REVERSE
// placement oracle.
//
// NON-SYMMETRIC symbols (validated by PDE15): A+/A- are the eigenvalue-weighted
// spectral projectors lambda_+- * P_+- onto the outgoing eigenspaces, so their
// ROW space (SVD V columns) is the LEFT characteristic combination (left
// eigenvector) in GENERAL -- not only for a symmetric symbol.  (Column space =
// right eigenvector; the incoming-BC guard correspondingly projects onto
// rowspace = V columns, see _validate_hyperbolic_incoming_bcs.)  PDE15's
// off-diagonal a*[[0,1],[4,0]] block -- where left (2,+-1) != right (1,+-2) --
// recovers the manufactured solution to ~5e-8 across all three modes and both
// wave-speed signs, confirming the rowspace closure is correct for the PSA
// (P,u) block; no left-vs-right correction is needed.
//-----------------------------------------------------------------------------
inline size_t
FFModel::_generate_hyperbolic_boundary_closure( bool const append )
{
  _faceConditions.clear();
  size_t added = 0;
  _hypClosureSkipped = 0; _hypClosureFaces = 0;
  for( auto const& bc : _blockClassification ){
    int             const  bid = bc.first;
    t_Classify      const& cls = bc.second;
    if( !cls.evolution_hyperbolic ) continue;

    auto itsym = _blockSymbol.find( bid );
    auto itfd  = _blockFaceData.find( bid );
    if( itsym == _blockSymbol.end() || itfd == _blockFaceData.end() ) continue;
    t_Symbol const& sym   = itsym->second;
    auto const&     fdata = itfd->second;

    size_t const nState = sym.vState.size();
    if( nState == 0 || sym.vEqn.size() != nState ) continue;   // need a square block
    size_t const evo_idx = cls.evolution_dom_idx;

    for( auto const& fd : fdata ){
      if( !fd.pdom_var ) continue;
      if( fd.dom_idx == evo_idx ) continue;                    // evolution dir: IC, not char BC
      if( fd.Aplus.n_rows  < (arma::uword)nState ||
          fd.Aminus.n_rows < (arma::uword)nState ) continue;

      struct EndSpec { int face; arma::mat const* A; };
      EndSpec const ends[2] = {
        { FFDom::LB, &fd.Aminus },   // LB outgoing = left-going  (A-)
        { FFDom::UB, &fd.Aplus  } }; // UB outgoing = right-going (A+)

      for( auto const& end : ends ){
        arma::mat U, V; arma::vec s;
        if( !arma::svd( U, s, V, *end.A ) || s.is_empty() ) continue;
        double const tol = 1e-8 * std::max( 1.0, s(0) );
        // rev317: which of the block's OWN equations already cover this face?  An equation covers it when its stored
        // mask (pooled over every stored copy of the same expression; DIAGNOSTIC rows excluded) contains the node set
        // the closure row would occupy: the face in this direction, ALL-LB in the evolution direction, the interior
        // in every other one.
        // rev317 (fix): FFDom masks are NOT bit flags -- ALL = 0, LB = -1, UB = -2, and `ALL - X` means "all nodes
        // except X".  Map a mask to the node classes it covers (1 = LB face, 2 = interior, 4 = UB face), mirroring
        // mask_has_face in _generate_algebraic_boundary_closure; any other value covers nothing (-> not covering,
        // the previous behaviour).  The closure row's own mask goes through the same map.
        auto nodes = []( int m ) -> int {
          int n = 0;
          if( m == FFDom::ALL || m == FFDom::LB || m == FFDom::ALL - FFDom::UB ) n |= 1;
          if( m == FFDom::ALL || m == FFDom::UB || m == FFDom::ALL - FFDom::LB ) n |= 4;
          if( m == FFDom::ALL || m == FFDom::ALL - FFDom::LB || m == FFDom::ALL - FFDom::UB
                              || m == FFDom::ALL - FFDom::LB - FFDom::UB )       n |= 2;
          return n;
        };
        auto covers = [&]( FFVar const& eqv ) -> bool {
          for( size_t dd = 0; dd < sym.vDom.size(); ++dd ){
            FFVar const& dv = sym.vDom[dd];
            int const need = nodes( ( dv.id() == fd.pdom_var->id() ) ? end.face
                                  : ( dd == evo_idx ? FFDom::ALL - FFDom::LB : FFDom::ALL - FFDom::LB - FFDom::UB ) );
            int have = 0;
            for( auto const& e : _mEqn ){
              if( e.var.id() != eqv.id() || !e.opt || e.opt->role == EqnRole::DIAGNOSTIC ) continue;
              auto it = e.dom.find( dv ); if( it != e.dom.end() ) have |= nodes( it->second );
            }
            if( !need || ( have & need ) != need ) return false;
          }
          return true;
        };
        std::vector<arma::vec> ells;                            // outgoing rowspace basis, in SVD order
        for( arma::uword js = 0; js < s.n_elem; ++js ) if( s(js) > tol ) ells.push_back( V.col(js) );
        std::vector<size_t> K;
        for( size_t k = 0; k < nState; ++k ) if( covers( sym.vEqn[k] ) ) K.push_back( k );
        std::vector<arma::vec> wts;                             // the weight vectors to emit
        if( K.empty() ){
          wts = ells;                                           // the previous behaviour, row for row
        }
        else{
          arma::mat C( ells.size(), K.size() );
          for( size_t j = 0; j < ells.size(); ++j ) for( size_t i = 0; i < K.size(); ++i ) C( j, i ) = ells[j]( K[i] );
          arma::vec sc = arma::svd( C );
          double const ctol = 1e-8 * std::max( 1.0, sc.is_empty() ? 0. : sc(0) );
          size_t rank = 0; for( arma::uword i = 0; i < sc.n_elem; ++i ) if( sc(i) > ctol ) ++rank;
          _hypClosureSkipped += rank;
          if( rank ) ++_hypClosureFaces;
          if( rank < ells.size() ){                            // partial: only the uncovered outgoing directions
            arma::mat Z = arma::null( C.t() );
            for( arma::uword iz = 0; iz < Z.n_cols; ++iz ){
              arma::vec w( nState, arma::fill::zeros );
              for( size_t j = 0; j < ells.size(); ++j ) w += Z( j, iz ) * ells[j];
              wts.push_back( w );
            }
            if( options.DISPLAY_LEVEL >= 2 ) std::cerr << "FFModel::setup ** auto_hyp_closure: block " << bid << " face " << ( end.face == FFDom::LB ? "LB" : "UB" )
                      << " in " << fd.pdom_var->name() << ": " << rank << " of " << ells.size()
                      << " outgoing row(s) supplied by the model -- appending the remaining " << Z.n_cols << "\n";
          }
        }
        {   // rev340: the per-face account (see t_FaceConditions)
          t_FaceConditions fc;
          fc.block_id  = bid;
          fc.direction = fd.pdom_var->name();
          fc.face      = end.face;
          fc.outgoing  = ells.size();
          fc.incoming  = (size_t)nState >= ells.size()? (size_t)nState - ells.size(): 0;
          fc.appended  = wts.size();
          fc.covered   = ells.size() >= wts.size()? ells.size() - wts.size(): 0;
          for( auto const& e : _mEqn ){
            if( !e.opt || e.opt->block_id != bid || e.opt->role == EqnRole::DIAGNOSTIC ) continue;
            auto it = e.dom.find( *fd.pdom_var );
            if( it != e.dom.end() && it->second == end.face ) ++fc.rows_at_face;
          }
          _faceConditions.push_back( fc );
          if( append && fc.rows_at_face > fc.incoming && options.DISPLAY_LEVEL >= 1 )   // not in the counting mode
            std::cerr << "FFModel::setup ** face " << ( end.face == FFDom::LB? "LB": "UB" ) << " in " << fc.direction
                      << " (block " << bid << "): " << fc.rows_at_face << " condition row(s) where only "
                      << fc.incoming << " characteristic(s) enter the domain -- data at the wrong end" << std::endl;
        }

        for( auto const& ell : wts ){

          FFVar comb; bool first = true;
          for( size_t k = 0; k < nState; ++k ){
            double const w = ell((arma::uword)k);
            if( std::fabs(w) < 1e-12 ) continue;
            comb  = first ? ( w * sym.vEqn[k] ) : ( comb + w * sym.vEqn[k] );
            first = false;
          }
          if( first ) continue;                                // degenerate all-zero row

          t_EqnDom md;
          for( size_t d = 0; d < sym.vDom.size(); ++d ){
            FFVar const& dv = sym.vDom[d];
            int mask;
            if( dv.id() == fd.pdom_var->id() ) mask = end.face;                        // LB / UB face
            else if( d == evo_idx )            mask = FFDom::ALL - FFDom::LB;           // drop IC node
            else                               mask = FFDom::ALL - FFDom::LB - FFDom::UB; // interior
            md.insert( { dv, mask } );
          }

          EqnOptions opt( EqnRole::INTERIOR, bid ); opt.participate_in_classification = false;   // rev312: sat/IC_AUTO are the solver defaults
          if( append ) _mEqn.push_back( { comb, md, _normalised_options( opt, md ) } );
          ++added;
        }
      }
    }
  }
  return added;
}

//-----------------------------------------------------------------------------
// OCFESLV::_generate_algebraic_boundary_closure
//
// A genuinely ALGEBRAIC state (no d_t, hence no characteristic and no boundary
// condition) is pinned ONLY by its defining algebraic constraint.  If a modeller
// collocates that constraint on the spatial INTERIOR (the natural analogue of a
// differential equation, whose boundaries are covered by BCs), the algebraic
// state's boundary DOFs are left unequationed -> non-square, silently, with the
// manufactured seed in the null space (M4/PDE19 = 36 short = u at z=LB,UB).
// Unlike a differential state, the constraint HOLDS POINTWISE at the boundary
// and the algebraic state is in general DISCONTINUOUS there (it inherits the
// interface jump of whatever derivative defines it), so the fix is to collocate
// the CONSTRAINT at the boundary -- NOT to add a C0-continuity row.  This is the
// non-hyperbolic analogue of _generate_hyperbolic_boundary_closure.
//
// Per non-hyperbolic block, for each interior non-LINK constraint (no d_t, so
// not a dynamic PDE; LINK = order-reduction aux, closed by its own continuity)
// pinning >=1 GENUINE (non-auxiliary) algebraic state via a BARE coupling: for
// each spatial face the constraint's mask EXCLUDES, append the constraint
// collocated at that face (classify=false, sat=true) -- but ONLY when no other
// block equation already collocates a pinned state there.  Conservative by
// construction: it adds a row only where an algebraic DOF is provably
// unequationed, so it can close an under-determined block but never over-
// determine one (a residual deficit, if any, still surfaces in the audit).
//-----------------------------------------------------------------------------
inline size_t
FFModel::_generate_algebraic_boundary_closure()
{
  FFVar const* tdom = _evolution_dom_set ? &_evolution_dom_var : nullptr;
  auto mask_has_face = []( int m, int face ) -> bool {
    return face == FFDom::LB
      ? ( m == FFDom::ALL || m == FFDom::LB || m == FFDom::ALL - FFDom::UB )
      : ( m == FFDom::ALL || m == FFDom::UB || m == FFDom::ALL - FFDom::LB );
  };
  auto bare_states = [&]( FFVar const& eqn ) -> std::set<FFVar,lt_FFVar> {
    std::set<FFVar,lt_FFVar> bare;
    auto sg = _dag->subgraph( 1, &eqn );
    for( auto const& op : sg.l_op ){
      if( op->sameid(typeid(FFPartial)) || op->sameid(typeid(FFIntegral)) ) continue;
      for( auto const* in : op->varin )
        if( in && _mVar.find(*in) != _mVar.end() ) bare.insert(*in);
    }
    return bare;
  };
  auto all_states = [&]( FFVar const& eqn ) -> std::set<FFVar,lt_FFVar> {
    std::set<FFVar,lt_FFVar> st;
    auto sg = _dag->subgraph( 1, &eqn );
    for( auto const& op : sg.l_op )
      for( auto const* in : op->varin )
        if( in && _mVar.find(*in) != _mVar.end() ) st.insert(*in);
    return st;
  };
  auto has_dt = [&]( FFVar const& eqn ) -> bool {
    auto sg = _dag->subgraph( 1, &eqn );
    for( auto const& op : sg.l_op ){
      if( !op->sameid(typeid(FFPartial)) ) continue;
      auto const* pop = mc::type_cast<FFPartial const>(op);
      for( auto const* operand : op->varin ){
        if( !operand || _mVar.find(*operand) == _mVar.end() ) continue;
        for( auto const& [iv,ord] : pop->Indep().expr )
          if( tdom && iv.id() == tdom->id() ) return true;
      }
    }
    return false;
  };

  // Accumulate new rows in a LOCAL buffer and splice them in only AFTER every
  // _mEqn scan completes.  Holding a reference into _mEqn (E / F below) while
  // push_back'ing to _mEqn would reallocate the vector and leave those
  // references dangling -> use-after-realloc segfault.  (The hyperbolic closure
  // sidesteps this because it builds `comb` from _blockSymbol, never an _mEqn
  // reference.)  _normalise_equation_options is pure (no _mEqn access), so it is
  // safe to call inside the loop.
  std::vector<t_Eqn> to_add;
#ifdef MC__OCFESLV_ALG_CLOSURE_PROBE
  size_t probe_total = 0, probe_uncovered = 0;
#endif
  for( auto const& bcpair : _blockClassification ){
    int const bid = bcpair.first;
    if( bcpair.second.evolution_hyperbolic ) continue;     // hyp closure handles these
    // DEGENERATE gate retired (2026-06-26).  Earlier this pass skipped every
    // non-PARABOLIC/DAE/ELLIPTIC block, on the theory that a constraint's
    // "excluded" boundary face might be an internal region interface pinned by
    // continuity built later.  Measurement (PDE2/PDE3 via the alg-closure probe)
    // showed the real over-fire was a DIMENSION mismatch -- the gas states Cg,Vg
    // are z-only but their equation spans {z,rd}, so the closure proposed a
    // spurious rd-face row -- which the dimension guard below now drops directly.
    // With that guard plus the coverage scan self-limiting the closure on every
    // non-hyperbolic block, the blanket skip is redundant (and would block a
    // multi-region PSA bed), so it is gone.  The genuinely reference-dependent
    // "state varies ALONG a multi-region interface" case is tracked separately
    // (latent; not exercised by the current suite).
#ifdef MC__OCFESLV_ALG_CLOSURE_PROBE
    EqnType const bt = bcpair.second.type;   // retained for the probe report only
#endif
    t_StructuralDecomp const d = _structural_dae_decomposition( bid );
    if( d.alg.empty() ) continue;
    std::set<FFVar,lt_FFVar> const alg_set( d.alg.begin(), d.alg.end() );

    size_t const n0 = _mEqn.size();                        // _mEqn is NOT mutated in this loop
    for( size_t ei = 0; ei < n0; ++ei ){
      t_Eqn const& E = _mEqn[ei];
      if( !E.opt->participate_in_classification || E.opt->block_id != bid ) continue;
      if( E.opt->role == EqnRole::LINK ) continue;          // order-reduction aux
      if( has_dt( E.var ) ) continue;                      // dynamic eqn, not a constraint
      auto const bs = bare_states( E.var );
      std::vector<FFVar> pinned;                           // GENUINE algebraic states pinned here
      for( auto const& s : bs )
        if( alg_set.count(s) && !_is_auxiliary_state(s) ) pinned.push_back(s);
      if( pinned.empty() ) continue;

      FFVar    const Evar = E.var;                          // copy: do not touch E after this
      t_EqnDom const Edom = E.dom;
      for( auto const& dm : Edom ){
        FFVar const& dv = dm.first; int const mask = dm.second;
        if( tdom && dv.id() == tdom->id() ) continue;      // evolution direction
        // A genuine-algebraic state has a boundary DOF only along the directions
        // it actually varies in.  If no pinned state depends on dv, an excluded
        // dv-face is not a real boundary for that state: collocating the
        // constraint there merely duplicates its existing pinning and over-
        // determines the block.  (PDE2/3: gas states Cg,Vg are z-only, but their
        // equation spans {z,rd} collocated at rd=LB to couple into the membrane,
        // so its mask excludes rd=UB -- without this guard the closure proposes a
        // spurious rd=UB row.  M4/PDE19: u varies along z and the face is z, so
        // the guard passes and the legitimate closure still fires.)
        bool dv_relevant = false;
        for( auto const& p : pinned ){
          auto itv = _mVar.find( p );
          if( itv != _mVar.end() && itv->second.count( dv ) ){ dv_relevant = true; break; }
        }
        if( !dv_relevant ){
#ifdef MC__OCFESLV_ALG_CLOSURE_PROBE
          std::cout << "  [alg-probe] blk=" << bid << " eqn=" << Evar.name()
                    << " dom=" << dv.name()
                    << "  [SKIP: no pinned state varies along this domain]\n";
#endif
          continue;
        }
        for( int face : { FFDom::LB, FFDom::UB } ){
          if( mask_has_face( mask, face ) ) continue;      // already collocated at this face
          bool covered = false;                            // any eqn pins a target state here?
          for( size_t ej = 0; ej < n0 && !covered; ++ej ){
            t_Eqn const& F = _mEqn[ej];
            if( F.opt->block_id != bid ) continue;
            auto itf = F.dom.find( dv );
            if( itf == F.dom.end() || !mask_has_face( itf->second, face ) ) continue;
            auto const fst = all_states( F.var );
            for( auto const& p : pinned ) if( fst.count(p) ){ covered = true; break; }
          }
          if( covered ) continue;
#ifdef MC__OCFESLV_ALG_CLOSURE_PROBE
          {
            ++probe_total;
            std::cout << "  [alg-probe] blk=" << bid
                      << " type=" << pde_type_name( bt )
                      << " eqn=" << Evar.name() << " state{";
            for( auto const& p : pinned ) std::cout << ' ' << p.name();
            std::cout << " } dom=" << dv.name()
                      << " face=" << ( face == FFDom::LB ? "LB" : "UB" )
                      << "  [would-ADD]"
                      << "\n";
            // GLOBAL (cross-block) coverage: does ANY equation anywhere collocate
            // at this (dom,face) AND reference a pinned target state?  A hit with
            // block!=bid is exactly what the per-block scan above misses
            // (hypothesis 1: pre-existing interface equation -> cheap fix = widen
            // the scan).  No hit at all means the cover, if any, is auto-generated
            // continuity built later in the interface plan (hypothesis 2 -> need
            // the continuity-policy predicate).
            size_t xcov = 0;
            for( size_t ej = 0; ej < n0; ++ej ){
              t_Eqn const& F = _mEqn[ej];
              auto itf = F.dom.find( dv );
              if( itf == F.dom.end() || !mask_has_face( itf->second, face ) ) continue;
              auto const fst = all_states( F.var );
              bool hits = false;
              for( auto const& p : pinned ) if( fst.count(p) ){ hits = true; break; }
              if( !hits ) continue;
              ++xcov;
              std::cout << "        cover: eqn=" << F.var.name()
                        << " blk=" << F.opt->block_id
                        << " role=" << (int)F.opt->role
                        << ( F.opt->block_id == bid ? "  (same block)"
                                                   : "  (CROSS-BLOCK)" ) << "\n";
            }
            if( !xcov ){
              ++probe_uncovered;
              std::cout << "        NO pre-existing eqn covers this face (any block)"
                           " -> cover, if any, is auto-generated continuity\n";
            }
          }
#endif
          t_EqnDom md = Edom;
          md[dv] = face;
          EqnOptions opt( EqnRole::INTERIOR, bid ); opt.participate_in_classification = false;   // rev312: sat/IC_AUTO are the solver defaults
          to_add.push_back( { Evar, md, _normalised_options( opt, md ) } );
        }
      }
    }
  }
#ifdef MC__OCFESLV_ALG_CLOSURE_PROBE
  std::cout << "  [alg-probe] SUMMARY: " << probe_total
            << " candidate boundary row(s); " << probe_uncovered
            << " with NO pre-existing cover (auto-continuity hypothesis), "
            << ( probe_total - probe_uncovered )
            << " already covered by a pre-existing equation the per-block scan missed.\n";
#endif
  _mEqn.reserve( _mEqn.size() + to_add.size() );
  for( auto& e : to_add ) _mEqn.push_back( std::move(e) );
  return to_add.size();
}

inline std::set<size_t>
FFModel::_value_slaved_algebraic_state_ids
() const
{
  // A GENUINE algebraic state pinned by a DERIVATIVE-FREE algebraic constraint
  // (no FFPartial/FFIntegral on any state, e.g. P - Rg c = 0) is determined
  // pointwise at every node; it is continuous exactly when the states it depends
  // on are, so an interface C0-continuity claim on it is redundant.  The
  // receiver-coupling redundancy proxy cannot see this (the implied state's claim
  // and the source state's claim sit on disjoint DOFs -- only the constraint links
  // them), so it must be detected structurally here.  Mirrors the lambdas of
  // _generate_algebraic_boundary_closure (bare/all/has_dt) so the same notion of
  // "algebraic constraint pinning a genuine state" is used at interfaces as at
  // boundaries.  Derivative-DEFINED algebraic states (u = -K dP/dz: all != bare)
  // are excluded -- they are in general discontinuous and correctly carry no claim.
  std::set<size_t> ids;
  FFVar const* tdom = _evolution_dom_set ? &_evolution_dom_var : nullptr;
  auto bare_states = [&]( FFVar const& eqn ) -> std::set<FFVar,lt_FFVar> {
    std::set<FFVar,lt_FFVar> bare;
    auto sg = _dag->subgraph( 1, &eqn );
    for( auto const& op : sg.l_op ){
      if( op->sameid(typeid(FFPartial)) || op->sameid(typeid(FFIntegral)) ) continue;
      for( auto const* in : op->varin )
        if( in && _mVar.find(*in) != _mVar.end() ) bare.insert(*in);
    }
    return bare;
  };
  auto all_states = [&]( FFVar const& eqn ) -> std::set<FFVar,lt_FFVar> {
    std::set<FFVar,lt_FFVar> st;
    auto sg = _dag->subgraph( 1, &eqn );
    for( auto const& op : sg.l_op )
      for( auto const* in : op->varin )
        if( in && _mVar.find(*in) != _mVar.end() ) st.insert(*in);
    return st;
  };
  auto has_dt = [&]( FFVar const& eqn ) -> bool {
    auto sg = _dag->subgraph( 1, &eqn );
    for( auto const& op : sg.l_op ){
      if( !op->sameid(typeid(FFPartial)) ) continue;
      auto const* pop = mc::type_cast<FFPartial const>(op);
      for( auto const* operand : op->varin ){
        if( !operand || _mVar.find(*operand) == _mVar.end() ) continue;
        for( auto const& [iv,ord] : pop->Indep().expr )
          if( tdom && iv.id() == tdom->id() ) return true;
      }
    }
    return false;
  };
  for( auto const& bcpair : _blockClassification ){
    int const bid = bcpair.first;
    t_StructuralDecomp const d = _structural_dae_decomposition( bid );
    if( d.alg.empty() ) continue;
    std::set<FFVar,lt_FFVar> const alg_set( d.alg.begin(), d.alg.end() );
    for( auto const& E : _mEqn ){
      if( !E.opt->participate_in_classification || E.opt->block_id != bid ) continue;
      if( E.opt->role == EqnRole::LINK ) continue;            // order-reduction aux
      if( has_dt( E.var ) ) continue;                        // dynamic eqn, not a constraint
      if( all_states( E.var ) != bare_states( E.var ) ) continue;  // derivative-DEFINED -> keep
      // rev182 R2: a BARE order-reduction auxiliary is a lifted derivative.  Under sigma reuse
      // the literal test above sees  w - cf*Dz_u  as derivative-free and classifies w value-
      // slaved; its seam claims are then dropped as implied by Dz_u's -- which are fabricated
      // and degenerate.  See through the alias: if any bare term is an auxiliary, the row is
      // derivative-DEFINED and its state is kept.  The EOS reference case (P = Rg c) has no
      // auxiliary term and is untouched.
      if( Options::_reuse_aware() ){
        bool aliased = false;
        for( auto const& s : bare_states( E.var ) ) if( _is_auxiliary_state(s) ){ aliased = true; break; }
        if( aliased ){ ++_reuseAwareKept; continue; }
      }
      for( auto const& s : bare_states( E.var ) )
        if( alg_set.count(s) && !_is_auxiliary_state(s) )
          ids.insert( (size_t)s.id().second );
    }
  }
  return ids;
}

// ======================================================================
// OCFESLV::_match_initial_closures
//   Exclusive attribution of INITIAL-role equations to the states they CLOSE.  See the declaration
//   for the rule and for the marching-override defect that motivated it.  The FFDep affinity test
//   is only consulted when an INITIAL residual references MORE THAN ONE candidate, so every model
//   whose ICs each name a single differential state is bit-for-bit unaffected.
// ======================================================================
inline std::map<FFVar,size_t,lt_FFVar>
FFModel::_match_initial_closures
( t_Eqns const& eqns, t_Var const& varmap, std::set<FFVar,lt_FFVar> const& cand,
  char const* who )
const
{
  std::map<FFVar,size_t,lt_FFVar> closure;
  if( cand.empty() ) return closure;
  (void)varmap;   // candidacy is carried by `cand`; varmap kept for signature symmetry/diagnostics

  // ---- Pass 0: per-INITIAL-equation candidate set + its affine (FFDep type L) subset ----
  struct t_Rec { size_t ndx = 0; std::vector<FFVar> refs, affine; };
  std::vector<t_Rec> rec;

  for( size_t ndx = 0; ndx < eqns.size(); ++ndx ){
    if( eqns[ndx].opt->role != EqnRole::INITIAL ) continue;
    FFVar expr = eqns[ndx].var;
    FFGraph* dag = static_cast<FFGraph*>( expr.dag() );
    if( !dag ) continue;
    FFSubgraph sg = dag->subgraph( 1, &expr );

    t_Rec c; c.ndx = ndx;
    std::vector<FFVar> leaf;          // every variable leaf: the FFDep independent set
    bool external = false;            // FFDep does not propagate through these external ops
    for( auto const* op : sg.l_op ){
      if( !op ) continue;
      if( op->type != FFOp::VAR ){
        if( op->sameid( typeid( FFPartial ) ) || op->sameid( typeid( FFIntegral ) )
         || op->sameid( typeid( FFEval    ) ) ) external = true;
        continue;
      }
      FFVar* pv = op->varout[0];
      if( !pv ) continue;
      leaf.push_back( *pv );
      if( cand.find( *pv ) != cand.end() ) c.refs.push_back( *pv );
    }
    if( c.refs.empty() ) continue;

    if( c.refs.size() == 1 || external || leaf.empty() ){
      c.affine = c.refs;              // nothing to disambiguate (or no usable dependency arithmetic)
    }
    else {
      // A Dirichlet-form IC is AFFINE in the state it pins; a state that only supplies an argument
      // reaches the residual through a nonlinear wrapper and comes back Q/P/R/N.
      std::vector<FFDep> dleaf( leaf.size() ), ddep( 1 );
      for( size_t i = 0; i < leaf.size(); ++i ) dleaf[i].indep( (int)leaf[i].id().second );
      bool ok = true;
      try {
        std::vector<FFVar> const vdep{ expr };
        dag->eval( vdep, ddep, leaf, dleaf );
      }
      catch( ... ){ ok = false; }
      if( ddep.size() != 1 ) ok = false;
      if( !ok ) c.affine = c.refs;
      else
        for( auto const& s : c.refs ){
          auto const d = ddep[0].dep( (int)s.id().second );
          if( d.first && d.second == FFDep::TYPE::L ) c.affine.push_back( s );
        }
      if( c.affine.empty() ) c.affine = c.refs;   // no affine candidate -> fall back to refs
    }
    rec.push_back( std::move( c ) );
  }

  // ---- Pass 1: unambiguous closure (exactly one still-unclaimed affine candidate) ----
  std::set<size_t> eqn_used;
  for( auto const& c : rec ){
    std::vector<FFVar> freev;
    for( auto const& s : c.affine ) if( !closure.count( s ) ) freev.push_back( s );
    if( freev.size() != 1 ) continue;
    closure[ freev[0] ] = c.ndx;
    eqn_used.insert( c.ndx );
  }

  // ---- Pass 2: residual ambiguity -- first unclaimed candidate, reported ----
  for( auto const& c : rec ){
    if( eqn_used.count( c.ndx ) ) continue;
    FFVar pick;
    for( auto const& s : c.affine ) if( !closure.count( s ) ){ pick = s; break; }
    if( !pick.dag() )
      for( auto const& s : c.refs ) if( !closure.count( s ) ){ pick = s; break; }
    if( !pick.dag() ) continue;       // every candidate this equation names is already closed
    closure[ pick ] = c.ndx;
    eqn_used.insert( c.ndx );
    if( options.DISPLAY_LEVEL >= 1 ){
      std::cerr << "OCFESLV::setup ** " << ( who ? who : "initial-closure" )
                << ": ambiguous INITIAL closure (equation " << c.ndx << " is affine in";
      for( auto const& s : c.affine ) std::cerr << " " << s;
      std::cerr << "); attributed to " << pick << ".\n";
    }
  }

  return closure;
}

inline bool
FFModel::_register
( FFVar const& var, std::vector<FFVar>& all, 
  std::map<FFVar::pt_idVar,size_t>& ndx )
{
  if( !var.dag() ) return false;
  auto const id = var.id();
  if( ndx.find( id ) != ndx.end() ) return false;
  ndx[id] = all.size();
  all.push_back( var );
  return true;
}

inline bool
FFModel::_register
( std::map<FFVar::pt_idVar,size_t>& ndxAllVars,
  std::vector<FFVar>& vAllVarsLoc,
  FFGraph* dagSrc,
  std::vector<FFVar> const& vCstSrc,
  t_Dom const& mDomSrc,
  t_Var const& mVarSrc,
  t_Var const& mInpSrc,
  t_Eqns const& mEqnSrc,
  t_Fcts const& mFctSrc,
  std::map<int,t_Symbol> const* bSymSrc,
  std::vector<t_AuxDef> const* auxDefSrc,
  t_Transitions const* trnSrc )
{
  std::vector<FFVar> vAllVarsSrc;
  ndxAllVars.clear();
  vAllVarsLoc. clear();

  for( auto const& cst : vCstSrc )        _register( cst, vAllVarsSrc, ndxAllVars );
  for( auto const& [var, dom] : mDomSrc ) _register( var, vAllVarsSrc, ndxAllVars );
  for( auto const& [var, dom] : mVarSrc ) _register( var, vAllVarsSrc, ndxAllVars );
  for( auto const& [var, dom] : mInpSrc ) _register( var, vAllVarsSrc, ndxAllVars );
  for( auto const& eqn : mEqnSrc ) _register( eqn.var, vAllVarsSrc, ndxAllVars );
  for( auto const& fct : mFctSrc ) _register( fct.var, vAllVarsSrc, ndxAllVars );
  if( trnSrc )                                              // transitions: their expressions are roots too
    for( auto const& tr : *trnSrc ){
      for( auto const& v : tr.left )  _register( v, vAllVarsSrc, ndxAllVars );
      for( auto const& v : tr.right ) _register( v, vAllVarsSrc, ndxAllVars );
      _register( tr.dom, vAllVarsSrc, ndxAllVars );
    }
  
  // Register variables in symbolic coefficient matrices of each block
  if( bSymSrc ){
    for( auto const& [blk, sym] : *bSymSrc )
      for( auto const& row : sym.vCoeff )
        for( auto const& var : row )      _register( var, vAllVarsSrc, ndxAllVars );
  }

  // Register auxiliary-definition roots.  The defining expression can be an
  // internal DAG node that is not otherwise exposed as a container key.
  if( auxDefSrc ){
    for( auto const& aux : *auxDefSrc ){
      _register( aux.aux,  vAllVarsSrc, ndxAllVars );
      _register( aux.expr, vAllVarsSrc, ndxAllVars );
      for( auto const& d : aux.dom )      _register( d, vAllVarsSrc, ndxAllVars );
      for( auto const& d : aux.diff_dom ) _register( d, vAllVarsSrc, ndxAllVars );
      _register( aux.parent, vAllVarsSrc, ndxAllVars );
    }
  }

  try{
    _dag->insert( dagSrc, vAllVarsSrc, vAllVarsLoc );
  }
  catch(...){
    return false;
  }

  return true;
}

inline FFVar
FFModel::_local
( FFVar const& var, std::vector<FFVar> const& all, 
  std::map<FFVar::pt_idVar,size_t> const& ndx )
{
  // Plain numeric FFVar constants need no DAG rebinding.  They can appear in
  // cached symbolic objects such as principal-symbol coefficients.
  if( !var.dag() ) return var;

  auto it = ndx.find( var.id() );
  if( it == ndx.end() )
    throw Exceptions( Exceptions::INDEX );
  return all.at( it->second );
};

inline std::vector<FFVar>
FFModel::_copy
( std::vector<FFVar> const& in, std::vector<FFVar> const& all, 
  std::map<FFVar::pt_idVar,size_t> const& ndx )
{
  std::vector<FFVar> out;
  out.reserve( in.size() );
  for( auto const& var : in )
    out.push_back( _local( var, all, ndx ) );
  return out;
}

inline std::set<FFVar,lt_FFVar>
FFModel::_copy
( std::set<FFVar,lt_FFVar> const& in, std::vector<FFVar> const& all,
  std::map<FFVar::pt_idVar,size_t> const& ndx )
{
  std::set<FFVar,lt_FFVar> out;
  for( auto const& var : in )
    out.insert( _local( var, all, ndx ) );
  return  out;
}

template <typename VAL>
inline std::map<FFVar,VAL,lt_FFVar>
FFModel::_copy
( std::map<FFVar,VAL,lt_FFVar> const& in, std::vector<FFVar> const& all,
  std::map<FFVar::pt_idVar,size_t> const& ndx )
{
  std::map<FFVar,VAL,lt_FFVar> out;
  for( auto const& [var, val] : in )
    out.insert( {_local( var, all, ndx ), val} );
  return out;
}

inline FFModel::t_Cont
FFModel::_copy_inp_cont
( t_Cont const& in, std::vector<FFVar> const& all,
  std::map<FFVar::pt_idVar,size_t> const& ndx )
{
  t_Cont out;
  for( auto const& [inp, dirs] : in ){
    std::map<FFVar,int,lt_FFVar> dirs_local;
    for( auto const& [dom, level] : dirs )
      dirs_local[ _local( dom, all, ndx ) ] = level;
    out[ _local( inp, all, ndx ) ] = std::move( dirs_local );
  }
  return out;
}

inline std::map<FFVar,FFModel::t_Fun,lt_FFVar>
FFModel::_copy
( std::map<FFVar,t_Fun,lt_FFVar> const& in, t_Dom const& domSrc,
  std::vector<FFVar> const& all, std::map<FFVar::pt_idVar,size_t> const& ndx )
{
  std::vector<std::pair<FFVar,FFVar>> dom_src_dst;
  dom_src_dst.reserve( domSrc.size() );
  for( auto const& [dsrc, dom] : domSrc ){
    (void)dom;
    dom_src_dst.emplace_back( dsrc, _local( dsrc, all, ndx ) );
  }

  std::map<FFVar,t_Fun,lt_FFVar> out;
  for( auto const& kv : in ){
    FFVar const& vsrc = kv.first;
    t_Fun const fun = kv.second;
    FFVar const vdst = _local( vsrc, all, ndx );
    out[vdst] = [fun, dom_src_dst]( t_Coord const& coord_dst ) -> double
    {
      t_Coord coord = coord_dst;
      for( auto const& [dsrc,ddst] : dom_src_dst ){
        auto it = coord_dst.find( ddst );
        if( it != coord_dst.end() ) coord[dsrc] = it->second;
      }
      return fun( coord );
    };
  }
  return out;
}

template <typename VAL>
inline std::multimap<FFVar,VAL,lt_FFVar>
FFModel::_copy
( std::multimap<FFVar,VAL,lt_FFVar> const& in, std::vector<FFVar> const& all,
  std::map<FFVar::pt_idVar,size_t> const& ndx )
{
  std::multimap<FFVar,VAL,lt_FFVar> out;
  for( auto const& [var, val] : in )
    out.insert( {_local( var, all, ndx ), val} );
  return out;
}

inline FFModel::t_Symbol
FFModel::_copy
( t_Symbol const& in, std::vector<FFVar> const& all,
  std::map<FFVar::pt_idVar,size_t> const& ndx )
{
  t_Symbol out;
  out.vEqn   = _copy( in.vEqn,   all, ndx );
  out.vState = _copy( in.vState, all, ndx );
  out.vDom   = _copy( in.vDom,   all, ndx );
  out.vCoeff.reserve( in.vCoeff.size() );
  for( auto const& row : in.vCoeff )
    out.vCoeff.push_back( _copy( row, all, ndx ) );
  return out;
};

inline bool
FFModel::_copy_usr_to_local
()
{
  _reset();
  _release_working_dag();
  
  _dag      = new FFGraph;
  _dagOwned = true;

  std::map<FFVar::pt_idVar,size_t> ndxAllVars;
  std::vector<FFVar> vAllVarsLoc;
  if( !_register( ndxAllVars, vAllVarsLoc, _usr._dagUsr, _usr._vCstUsr, _usr._mDomUsr,
                  _usr._mVarUsr, _usr._mInpUsr, _usr._mEqnUsr, _usr._mFctUsr, nullptr, nullptr, &_usr._mTrnUsr ) ){
    _reset();
    _release_working_dag();
    return false;
  }

  _vCstVal = _usr._vCstValUsr;
  _vCst = _copy( _usr._vCstUsr, vAllVarsLoc, ndxAllVars );
  _mDom = _copy( _usr._mDomUsr, vAllVarsLoc, ndxAllVars );

  _mInp.clear();
  for( auto const& [var, dom] : _usr._mInpUsr )
    _mInp.insert( { _local( var, vAllVarsLoc, ndxAllVars ),
                    _copy( dom, vAllVarsLoc, ndxAllVars ) } );
  // The working model now exists on the new DAG; a derived solver rebuilds what it derives from it.
  _on_model_imported( vAllVarsLoc, ndxAllVars );
  _mInpCont = _copy_inp_cont( _usr._mInpContUsr, vAllVarsLoc, ndxAllVars );

  _mVar.clear();
  for( auto const& [var, dom] : _usr._mVarUsr )
    _mVar.insert( { _local( var, vAllVarsLoc, ndxAllVars ),
                    _copy( dom, vAllVarsLoc, ndxAllVars ) } );

  _mEqn.clear();
  _mEqn.reserve( _usr._mEqnUsr.size() );
  for( auto const& eqn : _usr._mEqnUsr )
    _mEqn.push_back( { _local( eqn.var, vAllVarsLoc, ndxAllVars ),
                       _copy( eqn.dom, vAllVarsLoc, ndxAllVars ), eqn.opt } );

  _mFct.clear();
  _mFct.reserve( _usr._mFctUsr.size() );
  for( auto const& fct : _usr._mFctUsr ){
    t_Fct lfct;
    lfct.var   = _local( fct.var, vAllVarsLoc, ndxAllVars );
    lfct.kind  = fct.kind;
    lfct.point = _copy( fct.point, vAllVarsLoc, ndxAllVars );
    for( auto const& [dv, sd] : fct.side )                  // the SIDE of a point output (2026-09-29: was dropped here)
      lfct.side[ _local( dv, vAllVarsLoc, ndxAllVars ) ] = sd;
    lfct.grid  = _copy( fct.grid,  vAllVarsLoc, ndxAllVars );
    lfct.row0  = fct.row0;
    lfct.nrow  = fct.nrow;
    _mFct.push_back( std::move( lfct ) );
  }
  _mTrn.clear();                                            // transitions (add_transition) on the working DAG
  for( auto const& tr : _usr._mTrnUsr ){
    t_Transition ltr{ _copy( tr.left, vAllVarsLoc, ndxAllVars ), _copy( tr.right, vAllVarsLoc, ndxAllVars ),
                      _local( tr.dom, vAllVarsLoc, ndxAllVars ), tr.tau };
    if( tr.evaluated )                                        // each evaluation -> its operand: the explicit form, whose
      for( auto* side : { &ltr.left, &ltr.right } ){          // left is read at tau^- and right at tau^+ (checked)
        // symbolic replay with every leaf its own value, evaluations in STRIP mode (compose cannot: it substitutes
        // variables only, and an evaluation's result is an intermediate node)
        FFSubgraph sg = _dag->subgraph( side->size(), side->data() );
        std::vector<FFVar> vv, wk, res( side->size() );
        for( auto const* op : sg.l_op ) if( op && op->type == FFOp::VAR && op->varout[0] ) vv.push_back( *op->varout[0] );
        { FFEval::StripGuard const strip;
          _dag->eval( sg, wk, side->size(), side->data(), res.data(), vv.size(), vv.data(), vv.data() ); }
        *side = res;
      }
    _mTrn.push_back( ltr );
  }
  
  if( _usr._evolution_dom_setUsr ){
    _evolution_dom_var  = _local( _usr._evolution_dom_varUsr, vAllVarsLoc, ndxAllVars );
    _evolution_dom_user = _usr._evolution_dom_userUsr;
    _evolution_dom_set  = true;
  }
  else{
    _evolution_dom_var  = FFVar();
    _evolution_dom_user = false;
    _evolution_dom_set  = false;
  }

  _classVarRef = _copy( _usr._classVarRefUsr, vAllVarsLoc, ndxAllVars );
  _classInpRef = _copy( _usr._classInpRefUsr, vAllVarsLoc, ndxAllVars );
  _classDomRef = _copy( _usr._classDomRefUsr, vAllVarsLoc, ndxAllVars );

  // Copy reference functions and wrap them so that callers may key the
  // coordinate map with either the original user-domain variables or their
  // local DAG copies.
  _classVarRefFun = _copy( _usr._classVarRefFunUsr, _usr._mDomUsr, vAllVarsLoc, ndxAllVars );
  _classInpRefFun = _copy( _usr._classInpRefFunUsr, _usr._mDomUsr, vAllVarsLoc, ndxAllVars );

  return true;
}

inline bool
FFModel::set_model
( FFModel const& src )
{
  if( &src == this ) return true;
  if( !_usr._dagUsr ) return false;          // set() the user DAG first: there is nowhere to import into

  FFGraph* const dag = _usr._dagUsr;
  _reset();                                  // virtual: a derived solver clears its working state too
  _usr.clear_model();
  _release_working_dag();
  _on_model_replaced();
  _usr._dagUsr = dag;

  // Registration order, as everywhere else: constants, domains, states, inputs, equations, outputs.
  std::vector<FFVar> rootsSrc, rootsLoc;
  std::map<FFVar::pt_idVar,size_t> ndx;
  auto reg = [&]( FFVar const& v ){
    if( !v.dag() ) return;
    if( ndx.find( v.id() ) != ndx.end() ) return;
    ndx[v.id()] = rootsSrc.size();
    rootsSrc.push_back( v );
  };
  for( auto const& cst : src._usr._vCstUsr )        reg( cst );
  for( auto const& [var, dom] : src._usr._mDomUsr ) reg( var );
  for( auto const& [var, dom] : src._usr._mVarUsr ) reg( var );
  for( auto const& [var, dom] : src._usr._mInpUsr ) reg( var );
  for( auto const& eqn : src._usr._mEqnUsr )        reg( eqn.var );
  for( auto const& fct : src._usr._mFctUsr )        reg( fct.var );
  if( src._usr._evolution_dom_setUsr )              reg( src._usr._evolution_dom_varUsr );
  try{
    if( src._usr._dagUsr && src._usr._dagUsr != dag )
      dag->insert( src._usr._dagUsr, rootsSrc, rootsLoc );
    else
      rootsLoc = rootsSrc;                   // same user DAG: the variables are already the right ones
  }
  catch(...){
    return false;
  }

  _usr._vCstUsr    = _copy( src._usr._vCstUsr, rootsLoc, ndx );
  _usr._vCstValUsr = src._usr._vCstValUsr;
  _usr._mDomUsr.clear();
  for( auto const& [var, dom] : src._usr._mDomUsr )
    _usr._mDomUsr.insert( { _local( var, rootsLoc, ndx ), dom } );
  _usr._mVarUsr.clear();
  for( auto const& [var, dom] : src._usr._mVarUsr )
    _usr._mVarUsr.insert( { _local( var, rootsLoc, ndx ), _copy( dom, rootsLoc, ndx ) } );
  _usr._mInpUsr.clear();
  for( auto const& [var, dom] : src._usr._mInpUsr )
    _usr._mInpUsr.insert( { _local( var, rootsLoc, ndx ), _copy( dom, rootsLoc, ndx ) } );
  _usr._mInpDiscUsr.clear();
  for( auto const& [inp, decl] : src._usr._mInpDiscUsr ){
    std::map<FFVar,InpColloc,lt_FFVar> d;
    for( auto const& [dvar, ic] : decl ) d[ _local( dvar, rootsLoc, ndx ) ] = ic;
    _usr._mInpDiscUsr[ _local( inp, rootsLoc, ndx ) ] = std::move( d );
  }
  _usr._mInpContUsr = _copy_inp_cont( src._usr._mInpContUsr, rootsLoc, ndx );
  _usr._mEqnUsr.clear();
  _usr._mEqnUsr.reserve( src._usr._mEqnUsr.size() );
  for( auto const& eqn : src._usr._mEqnUsr )
    _usr._mEqnUsr.push_back( { _local( eqn.var, rootsLoc, ndx ), _copy( eqn.dom, rootsLoc, ndx ), eqn.opt } );
  _usr._mFctUsr.clear();
  _usr._mFctUsr.reserve( src._usr._mFctUsr.size() );
  for( auto const& fct : src._usr._mFctUsr ){
    t_Fct f;
    f.var   = _local( fct.var, rootsLoc, ndx );
    f.kind  = fct.kind;
    f.point = _copy( fct.point, rootsLoc, ndx );
    f.grid  = _copy( fct.grid,  rootsLoc, ndx );
    f.row0  = fct.row0;
    f.nrow  = fct.nrow;
    _usr._mFctUsr.push_back( std::move( f ) );
  }
  _usr._evolution_dom_setUsr  = src._usr._evolution_dom_setUsr;
  _usr._evolution_dom_userUsr = src._usr._evolution_dom_userUsr;
  if( src._usr._evolution_dom_setUsr )
    _usr._evolution_dom_varUsr = _local( src._usr._evolution_dom_varUsr, rootsLoc, ndx );
  _usr._classVarRefUsr = _copy( src._usr._classVarRefUsr, rootsLoc, ndx );
  _usr._classInpRefUsr = _copy( src._usr._classInpRefUsr, rootsLoc, ndx );
  _usr._classDomRefUsr = _copy( src._usr._classDomRefUsr, rootsLoc, ndx );
  _usr._classVarRefFunUsr = _copy( src._usr._classVarRefFunUsr, src._usr._mDomUsr, rootsLoc, ndx );
  _usr._classInpRefFunUsr = _copy( src._usr._classInpRefFunUsr, src._usr._mDomUsr, rootsLoc, ndx );

  _issetup    = false;
  _classified = false;
  _on_model_changed( ModelChange::DERIVATIVES );
  return true;
}

inline bool
FFModel::_check_model
( FFGraph* dag, std::vector<FFVar> const& vCst, std::vector<double> const* vCstVal,
  t_Dom const& mDom, t_Var const& mVar, t_Var const& mInp,
  t_Eqns const& mEqn, t_Fcts const& mFct,
  std::vector<FFSubgraph>* sgEqn, std::vector<FFSubgraph>* sgFct )
const
{
  if( !dag ){
    std::cerr << "  **ERROR: UNDEFINED DAG" << std::endl;
    return false;
  }

  if( vCstVal && !vCstVal->empty() && vCstVal->size() != vCst.size() ){
    std::cerr << "  **ERROR: CONSTANT VALUE VECTOR SIZE " << vCstVal->size()
              << " INCONSISTENT WITH NUMBER OF CONSTANTS " << vCst.size()
              << std::endl;
    return false;
  }

  auto same_dag = [dag]( FFVar const& v ) -> bool
  { return !v.dag() || v.dag() == dag; };

  auto is_constant = [&vCst]( FFVar const& var ) -> bool
  {
    if( !var.dag() ) return false;
    for( auto const& cst : vCst )
      if( var.id() == cst.id() ) return true;
    return false;
  };

  for( auto const& cst : vCst ){
    if( !cst.dag() || cst.dag() != dag ){
      std::cerr << "  **ERROR: CONSTANT " << cst
                << " DOES NOT BELONG TO THE MODEL DAG" << std::endl;
      return false;
    }
  }

  for( auto const& [var, dom] : mDom ){
    (void)dom;
    if( !var.dag() || var.dag() != dag ){
      std::cerr << "  **ERROR: DOMAIN " << var
                << " DOES NOT BELONG TO THE MODEL DAG" << std::endl;
      return false;
    }
  }

  auto check_var_domains = [&]( t_Var const& m, char const* label ) -> bool
  {
    for( auto const& [var, dom] : m ){
      if( !var.dag() || var.dag() != dag ){
        std::cerr << "  **ERROR: " << label << " " << var
                  << " DOES NOT BELONG TO THE MODEL DAG" << std::endl;
        return false;
      }
      for( auto const& dvar : dom ){
        if( mDom.find( dvar ) == mDom.end() ){
          std::cerr << "  **ERROR: DOMAIN " << dvar << " MISSING FOR "
                    << label << " " << var << std::endl;
          return false;
        }
        if( !same_dag( dvar ) ){
          std::cerr << "  **ERROR: DOMAIN " << dvar << " FOR " << label
                    << " " << var << " DOES NOT BELONG TO THE MODEL DAG"
                    << std::endl;
          return false;
        }
      }
    }
    return true;
  };

  if( !check_var_domains( mVar, "STATE" ) ) return false;
  if( !check_var_domains( mInp, "INPUT" ) ) return false;

  // Consistency check: equation domains should be consistent with variable domains.
  if( sgEqn ) sgEqn->resize( mEqn.size() );
  size_t ndx = 0;
  for( auto const& eqn : mEqn ){
    auto const& eqnvar = eqn.var;
    auto const& eqndom = eqn.dom;
    auto const& opt    = *eqn.opt;
    if( !eqnvar.dag() || eqnvar.dag() != dag ){
      std::cerr << "  **ERROR: EQUATION " << eqnvar
                << " DOES NOT BELONG TO THE MODEL DAG" << std::endl;
      return false;
    }

    (void)opt;
    std::set<FFVar,lt_FFVar> mindom;

    for( auto const& [var, lim] : eqndom ){
      auto itdom = mDom.find( var );

      // Check whether all required domains were defined.
      if( itdom == mDom.end() ){
        std::cerr << "  **ERROR: DOMAIN " << var << " MISSING" << std::endl;
        return false;
      }

      switch( lim ){
        case FFDom::ALL:
        case FFDom::ALL-FFDom::LB:
        case FFDom::ALL-FFDom::UB:
        case FFDom::ALL-FFDom::LB-FFDom::UB:
        case FFDom::LB:
        case FFDom::UB:                                                      break;
        default: std::cerr << "  **ERROR: MISSPECIFIED BOUNDARY/DOMAIN " << var
                           << " IN EQUATION " << eqnvar << std::endl;        return false;
      }
    }

    FFSubgraph sg = dag->subgraph( 1, &eqnvar );

    for( auto itop = sg.l_op.cbegin(); itop != sg.l_op.cend(); ++itop ){
      auto const& op = *itop;

      if( op->sameid( typeid( FFPartial ) )
       || op->sameid( typeid( FFIntegral ) ) )
        continue;
      if( op->type != FFOp::VAR || is_constant( *op->varout[0] ) )
        continue;

      // Check all participating variables in equation were defined.
      // Keep state and input iterators separate; comparing iterators from
      // different containers is undefined even if their types coincide.
      auto it_state = mVar.find( *op->varout[0] );
      auto it_input = mInp.find( *op->varout[0] );

      std::set<FFVar,lt_FFVar> dom;
      FFVar const* pvar = nullptr;
      if( it_state != mVar.end() ){
        pvar = &it_state->first;
        dom  = it_state->second;
      }
      else if( it_input != mInp.end() ){
        pvar = &it_input->first;
        dom  = it_input->second;
      }
      else{
        auto itdom = mDom.find( *op->varout[0] );
        if( itdom == mDom.end() ){
          std::cerr << "  **ERROR: VARIABLE " << *op->varout[0] << " MISSING" << std::endl;
          return false;
        }
        mindom.insert( itdom->first );
        continue;
      }

      if( dom.empty() ) continue;
      FFVar const& var = *pvar;

      // Check whether variable domain is a subset of equation domain.
      if( _subset( dom, eqndom ) )
        mindom.insert( dom.cbegin(), dom.cend() );

      // Identify extra independents.
      auto domfree = dom;

      // Search for quadrature terms.
      auto jtop = itop;
      for( ++jtop; jtop != sg.l_op.cend(); ++jtop ){
        if( (*jtop)->sameid( typeid( FFIntegral ) ) ){
          for( auto const& [v,e] : mc::type_cast<FFIntegral const>(*jtop)->Indep().expr ){
            (void)e;
            domfree.erase( v );
          }
        }
        else if( (*jtop)->sameid( typeid( FFEval ) ) ){
          // OpEval consumes its direction at a fixed coordinate, exactly like integration.
          for( auto const& [dv,z0] : mc::type_cast<FFEval const>(*jtop)->Coord() ){
            (void)z0;
            domfree.erase( dv );
          }
        }
      }

      // Still extra independents not participating in integrals?
      if( !_subset( domfree, eqndom ) ){
        std::cerr << "  **ERROR: VARIABLE DOMAIN " << var
                  << " INCONSISTENT WITH DOMAIN (" << _strset(eqndom,",")
                  << ") OF EQUATION " << eqnvar << std::endl;
        return false;
      }

      // Insert variable domain in mindom.
      mindom.insert( domfree.cbegin(), domfree.cend() );
    }

    // Check whether equation domain is tight/exact.
    auto unuseddom = _setminus( eqndom, mindom );
    if( unuseddom.size() == 1 ){
      std::cerr << "  **ERROR: DOMAIN VARIABLE " << _strset(unuseddom,",")
                << " DECLARED BUT NOT USED IN EQUATION " << eqnvar << std::endl;
      return false;
    }
    else if( unuseddom.size() > 1 ){
      std::cerr << "  **ERROR: DOMAIN VARIABLES (" << _strset(unuseddom,",")
                << ") DECLARED BUT NOT USED IN EQUATION " << eqnvar << std::endl;
      return false;
    }

    if( sgEqn ) (*sgEqn)[ndx] = std::move( sg );
    ++ndx;

  }

  // Consistency check and row count for output functions. Outputs are
  // evaluated after all equation residuals and do not participate in
  // principal-symbol classification or SAT/interface residual assembly.

  if( sgFct ) sgFct->resize( mFct.size() );
  size_t ndxf = 0;
  for( auto const& fct : mFct ){
    auto const& fctvar = fct.var;
    auto const& fctpoint = fct.point;
    auto const& fctgrid  = fct.grid;
    if( !fctvar.dag() || fctvar.dag() != dag ){
      std::cerr << "  **ERROR: OUTPUT " << fctvar
                << " DOES NOT BELONG TO THE MODEL DAG" << std::endl;
      return false;
    }

    if( fct.kind == FctKind::DISTRIBUTED && fctgrid.empty() ){
      std::cerr << "  **ERROR: DISTRIBUTED OUTPUT " << fctvar
                << " HAS NO DISTRIBUTED DOMAIN" << std::endl;
      return false;
    }

    std::set<FFVar,lt_FFVar> availabledom;

    for( auto const& [var, val] : fctpoint ){
      auto itdom = mDom.find( var );
      if( itdom == mDom.end() ){
        std::cerr << "  **ERROR: DOMAIN " << var << " MISSING IN OUTPUT "
                  << fctvar << std::endl;
        return false;
      }
      double lo = itdom->second.lo_dom, up = itdom->second.up_dom;
      _output_point_bounds( var, lo, up );   // a derived solver may widen the span (OCFESLV: the march)
      double const tol = 64. * DBL_EPSILON * std::max( 1., std::max( std::fabs(lo), std::fabs(up) ) );
      if( val < lo - tol || val > up + tol ){
        std::cerr << "  **ERROR: OUTPUT POINT " << var << "=" << val
                  << " OUTSIDE DOMAIN [" << lo << "," << up << "] IN "
                  << fctvar << std::endl;
        return false;
      }
      availabledom.insert( var );
    }

    for( auto const& [var, lim] : fctgrid ){
      auto itdom = mDom.find( var );
      if( itdom == mDom.end() ){
        std::cerr << "  **ERROR: DOMAIN " << var << " MISSING IN OUTPUT "
                  << fctvar << std::endl;
        return false;
      }
      switch( lim ){
        case FFDom::ALL:
        case FFDom::ALL-FFDom::LB:
        case FFDom::ALL-FFDom::UB:
        case FFDom::ALL-FFDom::LB-FFDom::UB:
        case FFDom::LB:
        case FFDom::UB: break;
        default:
          std::cerr << "  **ERROR: MISSPECIFIED BOUNDARY/DOMAIN " << var
                    << " IN OUTPUT " << fctvar << std::endl;
          return false;
      }
      availabledom.insert( var );
    }

    FFSubgraph sg = dag->subgraph( 1, &fctvar );
    std::set<FFVar,lt_FFVar> fctmindom;

    for( auto itop = sg.l_op.cbegin(); itop != sg.l_op.cend(); ++itop ){
      auto const& op = *itop;
      if( op->type != FFOp::VAR || is_constant( *op->varout[0] ) )
        continue;

      auto it_state = mVar.find( *op->varout[0] );
      auto it_input = mInp.find( *op->varout[0] );

      std::set<FFVar,lt_FFVar> dom;
      FFVar const* pvar = nullptr;
      if( it_state != mVar.end() ){
        pvar = &it_state->first;
        dom  = it_state->second;
      }
      else if( it_input != mInp.end() ){
        pvar = &it_input->first;
        dom  = it_input->second;
      }
      else{
        auto itdom = mDom.find( *op->varout[0] );
        if( itdom == mDom.end() ){
          std::cerr << "  **ERROR: VARIABLE " << *op->varout[0]
                    << " MISSING IN OUTPUT " << fctvar << std::endl;
          return false;
        }
        if( availabledom.find( itdom->first ) == availabledom.end() ){
          std::cerr << "  **ERROR: DOMAIN VARIABLE " << itdom->first
                    << " USED IN OUTPUT " << fctvar
                    << " BUT NO EVALUATION VALUE OR DISTRIBUTED LIMIT WAS PROVIDED" << std::endl;
          return false;
        }
        fctmindom.insert( itdom->first );
        continue;
      }

      if( dom.empty() ) continue;
      FFVar const& var = *pvar;

      auto domfree = dom;
      auto jtop = itop;
      for( ++jtop; jtop != sg.l_op.cend(); ++jtop ){
        if( (*jtop)->sameid( typeid( FFIntegral ) ) ){
          for( auto const& [v,e] : mc::type_cast<FFIntegral const>(*jtop)->Indep().expr ){
            (void)e;
            domfree.erase( v );
          }
        }
        else if( (*jtop)->sameid( typeid( FFEval ) ) ){
          // OpEval consumes its direction at a fixed coordinate, like integration.
          for( auto const& [dv,z0] : mc::type_cast<FFEval const>(*jtop)->Coord() ){
            (void)z0;
            domfree.erase( dv );
          }
        }
      }

      if( !_subset( domfree, availabledom ) ){
        std::cerr << "  **ERROR: VARIABLE DOMAIN " << var
                  << " INCONSISTENT WITH OUTPUT DOMAINS (" << _strset(availabledom,",")
                  << ") OF OUTPUT " << fctvar << std::endl;
        return false;
      }

      fctmindom.insert( domfree.cbegin(), domfree.cend() );
    }

    for( auto const& var : fctmindom ){
      if( availabledom.find( var ) == availabledom.end() ){
        std::cerr << "  **ERROR: DOMAIN VARIABLE " << var
                  << " REQUIRED BY OUTPUT " << fctvar
                  << " BUT NO EVALUATION VALUE OR DISTRIBUTED LIMIT WAS PROVIDED" << std::endl;
        return false;
      }
    }

    if( fct.kind == FctKind::DISTRIBUTED ){
      for( auto const& [var, lim] : fctgrid ){
        (void)lim;
        if( fctmindom.find( var ) == fctmindom.end() ){
          std::cerr << "  **ERROR: DISTRIBUTED DOMAIN VARIABLE " << var
                    << " DECLARED BUT NOT USED IN OUTPUT " << fctvar << std::endl;
          return false;
        }
      }
    }

    if( sgFct ) (*sgFct)[ndxf] = std::move( sg );
    ++ndxf;

  }

  return true;
}


inline bool
FFModel::_subset
( std::set<FFVar,lt_FFVar> const& set1, std::set<FFVar,lt_FFVar> const& set2 )
{
  for( auto const& el1 : set1 )
    if( set2.find( el1 ) == set2.cend() ) return false;
  return true;
}

inline std::set<FFVar,lt_FFVar>
FFModel::_setminus
( std::set<FFVar,lt_FFVar> set1, std::set<FFVar,lt_FFVar> const& set2 )
{
  // set1 passed as a copy - can be modified and returned
  for( auto const& el2 : set2 )
    if( !set1.erase( el2 ) )
      throw std::runtime_error("OCFESLV::setminus **error not a subset");     
  return set1;
}

inline std::string
FFModel::_strset
( std::set<FFVar,lt_FFVar> const& set1, std::string const& sep )
{
  std::ostringstream oss;
  bool first = true;
  for( auto const& el1 : set1 ){
    if( !first ) oss << sep; 
    oss << el1;
    first = false;
  }
  return oss.str();
}

template <typename U>
inline bool
FFModel::_subset
( std::set<FFVar,lt_FFVar> const& set1, std::map<FFVar,U,lt_FFVar> const& set2 )
{
  for( auto const& el1 : set1 )
    if( set2.find( el1 ) == set2.cend() ) return false;
  return true;
}

template <typename U>
inline std::map<FFVar,U,lt_FFVar>
FFModel::_setminus
( std::map<FFVar,U,lt_FFVar> set1, std::set<FFVar,lt_FFVar> const& set2 )
{
  // set1 passed as a copy - can be modified and returned
  for( auto const& el2 : set2 )
    if( !set1.erase( el2 ) )
      throw std::runtime_error("OCFESLV::setminus **error not a subset");     
  return set1;
}

template <typename U>
inline std::string
FFModel::_strset
( std::map<FFVar,U,lt_FFVar> const& set1, std::string const& sep )
{
  std::ostringstream oss;
  bool first = true;
  for( auto const& [el1,dum] : set1 ){
    if( !first ) oss << sep; 
    oss << el1;
    first = false;
  }
  return oss.str();
}



inline void
FFModel::report
( std::ostream& os )
const
{
  // A number as the report prints it (shortest of %g)
  auto num_str = []( double x ) -> std::string
  { std::ostringstream o; o << x; return o.str(); };
  // The bounds of domain @p n (by name: the working and declared DAGs hold different FFVar copies)
  auto dom_bounds = [&]( std::string const& n, double& lo, double& up ) -> bool
  {
    for( auto const* m : { &_mDom, &_usr._mDomUsr } )
      for( auto const& [v,d] : *m ) if( v.name() == n ){ lo = d.lo_dom; up = d.up_dom; return true; }
    return false;
  };
  // A region of domain @p n, explicitly: "t in (0, 10]", "t = 0"
  auto mask_str = [&]( std::string const& n, int mask ) -> std::string
  {
    double lo = 0., up = 0.;
    bool const b = dom_bounds( n, lo, up );
    std::string const L = b? num_str( lo ): "lo", U = b? num_str( up ): "up";
    switch( mask ){
      case FFDom::ALL:                     return n + " in [" + L + ", " + U + "]";
      case FFDom::ALL-FFDom::LB:           return n + " in (" + L + ", " + U + "]";
      case FFDom::ALL-FFDom::UB:           return n + " in [" + L + ", " + U + ")";
      case FFDom::ALL-FFDom::LB-FFDom::UB: return n + " in (" + L + ", " + U + ")";
      case FFDom::LB:                      return n + " = " + L;
      case FFDom::UB:                      return n + " = " + U;
      default:                             return n + " (mask " + std::to_string( mask ) + ")";
    }
  };
  // An expression of the DAG, written out by FFExpr (truncated beyond 160 characters)
  auto expr_str = []( FFGraph* g, FFVar const& v ) -> std::string
  {
    std::string e = v.name();
    if( g ) try{
      FFSubgraph sg = g->subgraph( 1, &v );
      std::vector<FFExpr> const ex = FFExpr::subgraph( g, sg );
      if( !ex.empty() ){ std::ostringstream o; o << ex[0]; e = o.str(); }
    }
    catch( ... ){}
    return e.size() > 160? e.substr( 0, 157 ) + "...": e;
  };
  auto role_str = []( EqnRole r ) -> char const*
  {
    switch( r ){
      case EqnRole::AUTO:       return "AUTO";
      case EqnRole::INTERIOR:   return "INTERIOR";
      case EqnRole::INITIAL:    return "INITIAL";
      case EqnRole::BOUNDARY:   return "BOUNDARY";
      case EqnRole::INTERFACE:  return "INTERFACE";
      case EqnRole::LINK:       return "LINK";
      case EqnRole::SURFACE:    return "SURFACE";
      case EqnRole::DIAGNOSTIC: return "DIAGNOSTIC";
    }
    return "?";
  };
  auto dom_list = []( std::set<FFVar,lt_FFVar> const& d ) -> std::string
  {
    std::string s;
    for( auto const& v : d ) s += ( s.empty()? "": "," ) + v.name();
    return s.empty()? std::string("-"): s;
  };

  // the revision is not printed (2026-10-04): revision() returns it on request
  os << "\nFFModel ** "
     << ( _issetup? "the WORKING model (after setup)": "the DECLARED model (setup has not run)" ) << "\n";
  if( !_dag ){ os << "  NO DAG\n"; return; }

  os << "\nDOMAINS (" << var_domain().size() << ")\n";
  for( auto const& [var,dom] : var_domain() )
    os << "  " << std::left << std::setw(14) << var.name() << std::right
       << " [" << dom.lo_dom << " : " << dom.up_dom << "]  elements=" << dom.n_elem << " nodes/element=" << dom.n_node
       << ( _evolution_dom_set && _evolution_dom_var.dag() && var.id() == _evolution_dom_var.id()
            ? ( _evolution_dom_user? "   <= EVOLUTION (set)": "   <= EVOLUTION (found)" ): "" ) << "\n";

  // The deferred-value plumbing -- the INPUT setup adds to hold each deferred value -- is an implementation detail:
  // the report shows the inputs the USER declared (2026-10-04).  Only those holding inputs are hidden.  (A deferred
  // value's SOURCE is any expression of the model -- for a value at a point of the evolution direction, a declared
  // state itself -- so states and LINK rows are never hidden on account of it: an earlier version did, and the
  // CSTR's C_B vanished from its report.)
  auto deferred_holder = [&]( FFVar const& v ){
    for( auto const& C : var_deferred() ) if( C.input.dag() && C.input.id() == v.id() ) return true;
    return false; };
  os << "\nSTATES (" << var_state().size() << ")\n";
  for( auto const& [var,dom] : var_state() ){
    os << "  " << std::left << std::setw(20) << var.name() << std::right << " on " << dom_list( dom ) << "\n";
  }

  size_t n_inputs = 0;
  for( auto const& [var,dom] : var_input() ) if( !deferred_holder( var ) ) ++n_inputs;
  os << "\nINPUTS (" << n_inputs << ")\n";
  for( auto const& [var,dom] : var_input() ){
    if( deferred_holder( var ) ) continue;
    os << "  " << std::left << std::setw(20) << var.name() << std::right << " on " << dom_list( dom ) << "\n";
  }

  os << "\nCONSTANTS (" << var_constant().size() << ")\n";
  for( size_t i = 0; i < var_constant().size(); ++i )
    os << "  " << std::left << std::setw(20) << var_constant()[i].name() << std::right
       << ( i < val_constant().size()? " = " + std::to_string( val_constant()[i] ): std::string() ) << "\n";

  os << "\nEQUATIONS (" << var_equation().size() << ")\n";
  for( auto const& eqn : var_equation() ){
    os << "  " << std::left << std::setw(10) << ( eqn.opt? role_str( eqn.opt->role ): "?" ) << std::right
       << " block=" << ( eqn.opt? eqn.opt->block_id: 0 ) << "   ";
    std::string d;
    for( auto const& [var,mask] : eqn.dom ) d += ( d.empty()? "": ", " ) + mask_str( var.name(), mask );
    os << "0 = " << expr_str( _dag, eqn.var ) << ( d.empty()? std::string(): "   on " + d ) << "\n";
  }

  bool const usr_out = _issetup && _usr._dagUsr && !_usr._mFctUsr.empty();
  t_Fcts const& fcts = usr_out? _usr._mFctUsr: var_output();
  FFGraph* fdag = usr_out? _usr._dagUsr: _dag;
  os << "\nOUTPUTS (" << fcts.size() << ")\n";
  for( auto const& fct : fcts ){
    os << "  " << expr_str( fdag, fct.var );
    std::string w;
    for( auto const& [v,c] : fct.point ){
      auto const is = fct.side.find( v );
      w += ( w.empty()? "": ", " ) + v.name() + " = " + num_str( c ) + ( is != fct.side.end() && is->second == FFDom::PLUS? "+": "" );
    }
    std::string g;
    for( auto const& [v,mask] : fct.grid ) g += ( g.empty()? "": ", " ) + mask_str( v.name(), mask );
    if( fct.kind == FctKind::DISTRIBUTED ) os << "   distributed on " << ( g.empty()? std::string("-"): g ) << ( w.empty()? "": "   at " + w );
    else if( !w.empty() )                  os << "   at " << w;
    os << "\n";
  }

  if( !_usr._mTrnUsr.empty() ){
    os << "\nTRANSITIONS (" << _usr._mTrnUsr.size() << ") -- left(tau^-) = right(tau^+), componentwise"
       << ( _issetup? "; validated at setup": "" ) << "\n";
    for( auto const& tr : _usr._mTrnUsr ){
      os << "  " << ( tr.evaluated? std::string("through evaluations"):
                      "at " + ( tr.dom.dag()? tr.dom.name(): std::string("t") ) + " = " + num_str( tr.tau ) ) << ":";
      for( size_t k = 0; k < tr.left.size() && k < tr.right.size(); ++k )
        os << ( k? ";": "" ) << "  " << expr_str( _usr._dagUsr, tr.left[k] ) << " = " << expr_str( _usr._dagUsr, tr.right[k] );
      os << "\n";
    }
  }


  os << "\nCLASSIFICATION\n";
  if( !_issetup || block_classification().empty() )
    os << "  (not classified)\n";
  else
    for( auto const& [blk,cls] : block_classification() ){
      os << "  block " << blk << ": " << pde_type_name( cls.type )
         << ( cls.symbol_rectangular? "  symbol RECTANGULAR": "  symbol square" )
         << ( cls.descriptor? "  descriptor": "" )
         // the evolution coefficient of the DIFFERENTIAL rows only (the symbol leaves the algebraic rows out): it
         // says whether those rows can be solved for their time derivatives, NOT the model's differential index,
         // which the INDEX REDUCTION and STRUCTURAL INDEX sections give (2026-10-04)
         << ( cls.At_singular? "  differential rows: evolution coefficient singular": "  differential rows: evolution coefficient regular" )
         << ( cls.evolution_hyperbolic? "  evolution-hyperbolic": "" )
         << ( cls.spatially_characteristic? "  spatially characteristic": "" )
         << ( cls.degenerate? "  DEGENERATE": "" ) << "\n";
    }
  if( _evolution_dom_set && _evolution_dom_var.dag() )
    os << "  evolution direction: " << _evolution_dom_var.name()
       << ( _evolution_dom_user? " (set by the model)": " (found during setup)" ) << "\n";
  else
    os << "  no evolution direction\n";

  if( !wellposedness().empty() ){
    os << "\nWELL-POSEDNESS: " << wellposedness().size() << " finding(s)\n";
    for( auto const& f : wellposedness() )
      os << "  block " << f.block_id << ": " << f.detail << "\n";
  }
  else if( _issetup )
    os << "\nWELL-POSEDNESS: nothing found against this model\n";

  if( !face_conditions().empty() ){
    os << "\nCONDITIONS PER FACE (characteristic count)\n";
    for( auto const& fc : face_conditions() ){
      os << "  block " << fc.block_id << "  " << fc.direction << " " << ( fc.face == FFDom::LB? "LB": "UB" )
         << ":  " << fc.outgoing << " leaving (" << fc.covered << " covered, " << fc.appended << " appended), "
         << fc.incoming << " entering; " << fc.rows_at_face << " condition row(s) declared here";
      if( fc.rows_at_face > fc.incoming ) os << "   <= EXCESS: data at the wrong end";
      os << "\n";
    }
  }

  if( _issetup && !block_classification().empty() && var_domain().size() > 1 ){
    os << "\nSTRUCTURAL INDEX, per direction\n";
    for( auto const& [blk,cls] : block_classification() ){
      os << "  block " << blk << ":";
      for( auto const& [var,dom] : var_domain() ){
        t_IndexResult const ir = structural_index( blk, var );
        os << "  " << var.name() << "=" << ir.index
           << ( _evolution_dom_set && _evolution_dom_var.dag() && var.id() == _evolution_dom_var.id()
                ? " (evolution)": "" );
      }
      os << "\n";
    }
    os << "  -- reduced in the evolution direction only; elsewhere this is information for the consumer\n";
  }

  if( !reduction_plan().empty() ){
    auto const& P = reduction_plan();
    os << "\nINDEX REDUCTION (dummy derivatives) -- highest differential index " << P.max_index
       << ( P.resolved? "": ", NOT fully resolved" ) << "\n";
    for( auto const& a : P.assigns )
      os << "  block " << a.block_id << ": constraint 0 = " << expr_str( _dag, a.constraint )
         << " differentiated " << a.n_diff << " time(s) along the evolution direction, which determines "
         << a.pinned_var.name() << "\n";
  }

  if( initial_data().reduced ){
    auto const& I = initial_data();
    os << "\nINITIAL DATA\n  " << I.declared << " INITIAL row(s) declared; " << I.differential
       << " state(s) differentiated in the evolution direction; the reduction implies " << I.hidden
       << " hidden\n  constraint level(s) at the initial point, so only " << I.free
       << " of the initial values may be CHOSEN.\n";
    if( I.redundant && I.materialised )
      os << "  the " << I.hidden << " hidden constraint level(s) were MATERIALISED as INITIAL rows at the initial"
            " point (REDUCE.HIDDEN_IC);\n  " << I.redundant << " declared row(s) are REDUNDANT with them: the initial"
            " point is OVER-DETERMINED (a least-squares\n  surplus) -- declare only the free data, or set"
            " REDUCE.HIDDEN_IC = false and make the declared data consistent by hand.\n";
    else if( I.redundant )
      os << "  " << I.redundant << " declared row(s) are REDUNDANT: they must be consistent with constraints this"
            " model does not enforce,\n  and nothing here checks that they are.\n";
    else if( I.materialised )
      os << "  the " << I.hidden << " hidden constraint level(s) were MATERIALISED as INITIAL rows at the initial"
            " point (REDUCE.HIDDEN_IC).\n";
    else if( I.declared < I.differential )
      os << "  " << ( I.differential - I.declared ) << " row(s) short of what collocation needs: the hidden"
            " constraints are not rows of this model (REDUCE.HIDDEN_IC materialises them).\n";
  }

  os << "\nDEGREES OF FREEDOM\n  rows - unknowns = " << dof_balance().str
     << ( dof_balance().balanced? "   (balanced)": "   (NOT balanced)" ) << "\n";
  os << std::endl;
}


inline void
FFModel::_build_dynamic_form
()
{
  _dynamicForm = t_DynamicForm();
  if( !_evolution_dom_set || !_evolution_dom_var.dag() ){
    _dynamicForm.blocker = "no evolution direction is set, so there is nothing for an integrator to advance";
    return;
  }

  _dynamicForm.lumped = true;                       // every state on the evolution direction alone?
  for( auto const& [var,dom] : _mVar )
    for( auto const& d : dom )
      if( d.id() != _evolution_dom_var.id() ) _dynamicForm.lumped = false;

  _dynamicForm.differential   = _initialData.differential;
  _dynamicForm.algebraic      = _mVar.size() >= _initialData.differential
                              ? _mVar.size() - _initialData.differential: 0;
  _dynamicForm.declared_index = _reductionPlan.max_index;
  _dynamicForm.reduced        = !_reductionPlan.assigns.empty();
  _dynamicForm.resolved       = _reductionPlan.assigns.empty() || _reductionPlan.resolved;
  _dynamicForm.free_initial   = _initialData.reduced? _initialData.free: _initialData.differential;

  {   // are the derivatives structurally decoupled?  Not "is the matrix diagonal": a non-unit coefficient and a
      // permuted assignment are both explicit, and both would fail a diagonality test.
    std::map<std::string,size_t> rows_per_state;
    for( auto const& eqn : _mEqn ){
      if( !eqn.opt || eqn.opt->role != EqnRole::INTERIOR ) continue;
      std::set<std::string> here;
      auto sg = _dag->subgraph( 1, &eqn.var );
      for( auto const& op : sg.l_op ){
        if( !op->sameid( typeid(FFPartial) ) ) continue;
        auto const* pop = mc::type_cast<FFPartial const>( op );
        if( !pop ) continue;
        for( size_t jj = 0; jj < op->varin.size(); ++jj ){
          FFVar const* operand = op->varin[jj];
          if( !operand || _mVar.find( *operand ) == _mVar.end() ) continue;
          for( auto const& [indep_var, ord] : pop->Indep().expr ){ (void)ord;
            if( indep_var.id() == _evolution_dom_var.id() ) here.insert( operand->name() );
          }
        }
      }
      if( here.size() > 1 ) _dynamicForm.decoupled = false;      // this row couples several derivatives
      for( auto const& n : here ) ++rows_per_state[n];
    }
    for( auto const& [n,c] : rows_per_state )
      if( c > 1 ) _dynamicForm.decoupled = false;                // this derivative is spread over several rows
  }

  for( auto const& eqn : _mEqn ){                   // anything imposed at the FAR END of the evolution domain?
    auto it = eqn.dom.find( _evolution_dom_var );
    if( it != eqn.dom.end() && it->second == FFDom::UB ) _dynamicForm.initial_value = false;
  }

  for( auto const& f : _findings )                  // the blocker: the first thing that stops an integrator
    if( f.kind == t_Finding::COMPLEX_CHARACTERISTICS ){
      _dynamicForm.blocker = "the characteristic speeds are not real: no method can advance this model";
      return;
    }
  if( !_dynamicForm.lumped )
    _dynamicForm.blocker = "distributed: the spatial directions must be discretised before an integrator sees it"
                           " (method of lines)";
  else if( !_dynamicForm.initial_value )
    _dynamicForm.blocker = "something is imposed at the far end of the evolution domain, so this is a boundary-value"
                           " problem in that direction: collocation solves it, an integrator cannot advance it";
  else if( _dynamicForm.reduced && !_dynamicForm.resolved )
    _dynamicForm.blocker = "a high index was found but not resolved, so the system handed over would still be"
                           " high index";
  else if( !_dynamicForm.decoupled )
    _dynamicForm.blocker = "the time derivatives are coupled through a mass matrix -- still an ODE, and collocation"
                           " solves it as it stands, but an explicit integrator needs dy/dt = f, so either that"
                           " matrix is inverted or the model goes to IDAS as a residual";
}
inline void
FFModel::_collect_findings
()
{
  _findings.clear();

  for( auto const& [blk,cls] : _blockClassification ){
    if( cls.type == COMPLEX_CHARACTERISTIC )
      _findings.push_back( { t_Finding::COMPLEX_CHARACTERISTICS, blk,
        "characteristic speeds are not real, so the solution does not depend continuously on its data:"
        " no method will give meaningful results until the model is corrected" } );
    for( auto const& [var,dom] : _mDom ){
      t_IndexResult const ir = _structural_index_analysis_uncached( blk, &var );
      if( ir.index < 0 )
        _findings.push_back( { t_Finding::SINGULAR_INDEX, blk,
          "structurally singular with respect to " + var.name()
          + ": a variable no differentiation in that direction exposes, so nothing marching or shooting along it"
            " can determine the solution" } );
    }
  }

  if( !_dofBalance.balanced )
    _findings.push_back( { t_Finding::DOF_IMBALANCE, 0,
      "rows - unknowns = " + _dofBalance.str
      + ": a surplus is least-squared and a deficit leaves the system underdetermined" } );

  for( auto const& fc : _faceConditions )
    if( fc.rows_at_face > fc.incoming )
      _findings.push_back( { t_Finding::FACE_EXCESS, fc.block_id,
        std::to_string( fc.rows_at_face ) + " condition row(s) at the " + ( fc.face == FFDom::LB? "LB": "UB" )
        + " face in " + fc.direction + " where only " + std::to_string( fc.incoming )
        + " characteristic(s) enter: the excess is data at the wrong end" } );

  if( _initialData.redundant )
    _findings.push_back( { t_Finding::REDUNDANT_INITIAL_DATA, 0,
      std::to_string( _initialData.redundant ) + " of the " + std::to_string( _initialData.declared )
      + " initial row(s) cannot be chosen -- the reduction implies them -- "
      + ( _initialData.materialised
          ? std::string( "and the hidden constraints are rows (REDUCE.HIDDEN_IC), so the initial point is over-determined" )
          : std::string( "so they must be consistent with constraints this model does not enforce" ) ) } );
}
inline void
FFModel::_audit_dof_balance
()
{
  _initialData = t_InitialData();
  for( auto const& eqn : _mEqn )
    if( eqn.opt && eqn.opt->role == EqnRole::INITIAL ) ++_initialData.declared;
  for( auto const& a : _reductionPlan.assigns )
    if( a.n_diff > 0 ) _initialData.hidden += (size_t)a.n_diff;
  _initialData.materialised = options.REDUCE.HIDDEN_IC;
  if( _initialData.materialised && _initialData.declared >= _initialData.hidden )
    _initialData.declared -= _initialData.hidden;   // the materialised levels are rows, but not the modeller's
  _initialData.reduced = !_reductionPlan.assigns.empty();
  if( _evolution_dom_set && _evolution_dom_var.dag() ){
    std::set<FFVar,lt_FFVar> diff_states;                  // states differentiated in the evolution direction
    for( auto const& eqn : _mEqn ){
      auto sg = _dag->subgraph( 1, &eqn.var );
      for( auto const& op : sg.l_op ){
        if( !op->sameid( typeid(FFPartial) ) ) continue;
        auto const* pop = mc::type_cast<FFPartial const>( op );
        if( !pop ) continue;
        for( size_t jj = 0; jj < op->varin.size(); ++jj ){
          FFVar const* operand = op->varin[jj];
          if( !operand || _mVar.find( *operand ) == _mVar.end() ) continue;
          for( auto const& [indep_var, ord] : pop->Indep().expr )
            if( indep_var.id() == _evolution_dom_var.id() ){ diff_states.insert( *operand ); break; }
        }
      }
    }
    _initialData.differential = diff_states.size();
  }
  _initialData.free      = _initialData.differential > _initialData.hidden
                         ? _initialData.differential - _initialData.hidden : 0;
  _initialData.redundant = _initialData.declared > _initialData.free
                         ? _initialData.declared - _initialData.free : 0;

  _dofBalance = t_DofBalance();
  typedef std::set<std::string>    t_Mono;      // a multilinear monomial: the domains whose N it multiplies
  typedef std::map<t_Mono,long>    t_Poly;

  auto mul = []( t_Poly const& p, t_Poly const& f ) -> t_Poly
  {
    t_Poly out;
    for( auto const& [m1,c1] : p )
      for( auto const& [m2,c2] : f ){
        t_Mono m = m1; m.insert( m2.cbegin(), m2.cend() );
        out[m] += c1 * c2;
      }
    return out;
  };
  auto factor = []( std::string const& name, int mask ) -> t_Poly
  {
    switch( mask ){
      case FFDom::ALL:                     return { { {name}, 1 } };                        // N
      case FFDom::ALL-FFDom::LB:
      case FFDom::ALL-FFDom::UB:           return { { {name}, 1 }, { t_Mono(), -1 } };      // N-1
      case FFDom::ALL-FFDom::LB-FFDom::UB: return { { {name}, 1 }, { t_Mono(), -2 } };      // N-2
      default:                             return { { t_Mono(), 1 } };                      // a single face
    }
  };

  t_Poly diff;
  for( auto const& eqn : _mEqn ){                       // rows
    t_Poly p; p[ t_Mono() ] = 1;
    for( auto const& [var,mask] : eqn.dom ) p = mul( p, factor( var.name(), mask ) );
    for( auto const& [m,c] : p ) diff[m] += c;
  }
  for( auto const& [var,dom] : _mVar ){                 // unknowns
    t_Mono m;
    for( auto const& d : dom ) m.insert( d.name() );
    diff[m] -= 1;
  }

  std::ostringstream os;
  for( auto it = diff.cbegin(); it != diff.cend(); ){    // drop the zeros, then write it out, highest degree first
    if( !it->second ) it = diff.erase( it );
    else ++it;
  }
  std::vector<std::pair<t_Mono,long>> terms( diff.cbegin(), diff.cend() );
  std::sort( terms.begin(), terms.end(),
             []( auto const& a, auto const& b ){ return a.first.size() != b.first.size()? a.first.size() > b.first.size()
                                                                                        : a.first < b.first; } );
  for( auto const& [m,c] : terms ){
    std::string t;
    for( auto const& n : m ) t += ( t.empty()? "": "*" ) + ( "N_" + n );
    if( t.empty() ) t = "1";
    if( !os.str().empty() ) os << ( c < 0 ? " - ": " + " );
    else if( c < 0 )        os << "-";
    long const a = c < 0 ? -c : c;
    if( a != 1 || t == "1" ) os << a << ( t == "1"? "": "*" );
    if( t != "1" ) os << t;
  }
  _dofBalance.diff     = diff;
  _dofBalance.balanced = diff.empty();
  _dofBalance.str      = diff.empty()? "0": os.str();

  if( !_dofBalance.balanced && options.DISPLAY_LEVEL >= 1 )
    std::cerr << "FFModel::setup ** DOF: the model is NOT balanced -- rows - unknowns = " << _dofBalance.str
              << " (a surplus is least-squared, a deficit leaves the system underdetermined)" << std::endl;
}
inline bool
FFModel::setup
()
{
  std::lock_guard<std::recursive_mutex> dag_lock_( dag_mutex() );   // see dag_mutex()
  // rev142b: cap the numerical backends for the duration of setup() (MAXTHREAD, 0 = no
  // cap).  RAII because setup() has many early returns AND the w-decide restore
  // re-enters setup() from inside setup(); the guard's depth counter handles both.
  // rev308: the caller's configuration lands here, BEFORE the thread cap reads MAXTHREAD and before any
  // phase reads an option.
  _apply_options();
  if( !_check_transitions() ){ _setupStatus = SetupStatus::INCONSISTENT_MODEL;  return false; }   // see add_transition
  if( !_check_equation_reductions() ){ _setupStatus = SetupStatus::INCONSISTENT_MODEL;  return false; }   // the causality refusal

  t_ThreadCap const _tcap( options.MAXTHREAD );

  // rev150 = rev144 (revision of record) + this override ONLY.  The rev145-149
  // faithful-projection lineage is retired to the diagnostic record: it proved
  // the forked IC_STRONG deficiency is a def0-class implication over the
  // dropped pivot rows (pattern-complete columns, rank invariant across three
  // independent structural repairs) and that the rev146 -tau row is unsound as
  // a default (turns tau into jump slack; PDE5 legacy STRONG dup-spread
  // regressed 3e-16 -> 5e-02, caught by the audit AND the driver).  Under
  // CRONOS_STRONG_VIA_TRACE=1 an IC_STRONG setup builds and solves the full

  _issetup = false;
  _classified = false;

  // A derived solver clears what a re-setup invalidates on its side.
  if( !_on_setup_begin() ) return false;


  // --- optional setup phase timing -------------------------------------------
  // When DISPLAY_LEVEL>=2, prints the wall-clock spent in each major setup phase
  // so the dominant cost on large problems can be localized.  Zero overhead when
  // DISPLAY_LEVEL==0 (a clock read per phase, no printing).
  _t_phase  = std::chrono::steady_clock::now();
  _t_setup0 = _t_phase;

  // Preflight the user-facing model before copying it into the private
  // working DAG. This catches malformed user containers before order
  // reduction or classification can operate on them.
  if( !_validate_model() )
    { _setupStatus = SetupStatus::INCONSISTENT_MODEL; return _issetup; }

  if( !_copy_usr_to_local() )
    { _setupStatus = SetupStatus::MODEL_COPY_FAILED; return _issetup; }


  // Solver: the evolution hoist and, when marching, the collapse to the first window -- order reduction then
  // runs against the horizon actually being solved.
  if( !_on_before_reduction() ) return false;

  _reindex_controls();   // domains are declared by now, so the DOF counts are final

  _issetup = true;

  // Optionally rewrite high-order PDEs to first-order before any setup data
  // structures, classification, or derivative sparsity caches are built.
  if( options.REDUCE.ORDER != Options::RED_NONE ){
    bool const reuse_aux = ( options.REDUCE.ORDER == Options::RED_FULL );
    try{
      _reduce_order( reuse_aux );
    }
    catch( OCBase::Exceptions const& ){
      // A nested evolution-in-evolution reduction is refused loudly (reduce_order set the status and
      // printed the reason).  Surface it as a failed setup rather than propagating; any other
      // reduce_order exception is re-thrown.
      if( _setupStatus == SetupStatus::CAPTURE_NESTED_REFUSED ){ _issetup = false; return false; }
      throw;
    }
  }
  _phase( "copy + reduce_order" );
  // Transitions (see _transitions_lifted) are lowered AFTER order reduction: their rows read x(tau^-) / x(tau^+)
  // IN-SOLVE through pinned evaluations, which the reduction would capture as post-solve latches (and its
  // substitution map is applied to every equation, so a node shared with a captured output would be replaced too).
  if( !_lower_transitions() ){ _issetup = false;  _setupStatus = SetupStatus::INCONSISTENT_MODEL;  return false; }

  // Auto differential-elimination (Pantelides increment 1; OFF by default).  Rewrites
  // ELIMINABLE algebraically-determined derivative couplings (substitute + relocate) on
  // _mEqn BEFORE the collocation system is built, so an affected block classifies square.
  if( options.AUTO.DIFF_ELIM ) _auto_diff_eliminate();
  _phase( "auto_diff_elim" );

  // Solver: the collocation arrays and the node cache, which the classification's mesh-wavenumber probe reads.
  if( !_on_model_discretise() ) return false;

  // Automatic classification is a setup-time discretisation decision:
  // IC_AUTO may select value matching or characteristic/upwind interfaces,
  // which in turn changes residual evaluation and Jacobian sparsity.
  if( options.CLASSIFY.MODE != Options::CLASS_NONE ){
    if( !_classify_pde() ){
      if( options.CLASSIFY.MODE == Options::CLASS_STRICT ){ _setupStatus = SetupStatus::CLASSIFICATION_FAILED; return false; }
      _classified = false;
      _blockSymbol.clear();
      _blockClassification.clear();
      _blockFaceData.clear();
      _face_data.clear();
    }
  }
  _phase( "classify_pde" );


  // Stage-2 high-index reduction: assemble the persisted Pantelides reduction
  // plan (decision layer; no equation mutation here -- _reduce_high_index executes it).
  if( _classified && options.CLASSIFY.MODE != Options::CLASS_NONE ){
    _build_reduction_plan();     // decision layer (Pantelides matching)
    _reduce_high_index();        // execution layer (FAD reduction + consistent-IC/re-pivot)
    _phase( "reduce_high_index" );
  }

  // Solver: a validation (unconditional since rev319) of user-supplied incoming BCs at hyperbolic boundary faces, before the
  // outgoing closure is generated.
  if( !_on_before_closures() ) return false;

  // Auto-generate the outgoing-characteristic boundary closure for hyperbolic
  // blocks.  Gated by CRONOS_AUTO_HYP_CLOSURE (environment-only since rev318).  Classification, the incoming-BC
  // guard and this generator all run BEFORE the consistency/count pass, the
  // dependency map and the eval-plans, so the rows appended here flow through all
  // three caches in a SINGLE pass -- no append-then-reprep is needed.
  // 2026-10-07 (WORKPLAN 3.C): an option again, AUTO.HYP_CLOSURE.  Off, the model must close every outflow face
  // itself: if rows are missing the system would be UNDER-determined and solved silently wrong (measured: 10% on
  // scalar advection), so setup refuses -- the same rule as consistent initial data for a high-index DAE.
  if( !options.AUTO.HYP_CLOSURE && _classified ){
    // The generator's coverage test recognises a face closed by the block's own PDE EXTENDED to it (rev317); a model
    // may also close it with rows written for the purpose (OCFE_PDE14/15's manual oracle).  So a face is short only
    // if BOTH say so: directions uncovered by the PDE (fc.appended), and fewer boundary/closure rows reaching the face
    // than the block has states (initial and diagnostic rows excluded).
    _generate_hyperbolic_boundary_closure( false );
    auto reaches = []( int mask, int face ){
      return mask == FFDom::ALL || mask == face
          || ( face == FFDom::UB && mask == FFDom::ALL - FFDom::LB )
          || ( face == FFDom::LB && mask == FFDom::ALL - FFDom::UB ); };
    size_t n_missing = 0;  std::ostringstream where;
    for( auto const& fc : _faceConditions ){
      if( !fc.appended ) continue;
      size_t rows = 0;
      for( auto const& e : _mEqn ){
        if( !e.opt || e.opt->block_id != fc.block_id ) continue;
        if( e.opt->role == EqnRole::INITIAL || e.opt->role == EqnRole::DIAGNOSTIC ) continue;
        for( auto const& [dv, mask] : e.dom )
          if( dv.name() == fc.direction && reaches( mask, fc.face ) ){ ++rows; break; }
      }
      size_t const need = fc.incoming + fc.outgoing;
      size_t const short_by = std::min( fc.appended, need > rows? need - rows: size_t(0) );
      if( !short_by ) continue;
      n_missing += short_by;
      where << " [block " << fc.block_id << ", " << fc.direction << " " << ( fc.face == FFDom::LB? "LB": "UB" ) << ": " << short_by << "]";
    }
    if( n_missing ){
      std::cerr << "FFModel::setup ** AUTO.HYP_CLOSURE is off and " << n_missing << " outgoing-characteristic row(s) are"
                   " MISSING at outflow face(s):" << where.str();
      std::cerr << " -- close them in the model (e.g. the PDE extended to the face) or turn AUTO.HYP_CLOSURE on"
                << std::endl;
      _issetup = false;  _setupStatus = SetupStatus::HYP_CLOSURE_MISSING;  return false;
    }
  }
  if( options.AUTO.HYP_CLOSURE && _classified ){
    size_t const n_clo = _generate_hyperbolic_boundary_closure();
  if( _hypClosureSkipped )
    if( options.DISPLAY_LEVEL >= 2 ) std::cerr << "FFModel::setup ** auto_hyp_closure: " << _hypClosureSkipped << " outgoing-characteristic row(s) at "
              << _hypClosureFaces << " face(s) already supplied by the model's own equations -- not duplicated (rev317)\n";
    if( n_clo ){
      if( options.DISPLAY_LEVEL >= 1 )
        if( options.DISPLAY_LEVEL >= 2 ) std::cerr << "FFModel::setup ** auto_hyp_closure: appended " << n_clo
                  << " outgoing-characteristic boundary row(s) before the "
                     "consistency/cache build (single pass)\n";
      _phase( "auto_hyp_closure (gen)" );
    }
  }

  // Auto-close the boundary DOFs of genuinely algebraic states (no BC) by
  // collocating their defining constraint at the spatial faces the modeller's
  // interior collocation excludes.  Gated by CRONOS_AUTO_ALG_CLOSURE (environment-only since rev316).  Same
  // single-pass slot as the hyperbolic closure: appended rows flow through the
  // consistency/dependency/eval caches built below without an append-then-reprep
  // cycle.
#ifdef MC__OCFESLV_ALG_CLOSURE_PROBE
  bool const run_alg_closure = _classified;   // probe: force the scan regardless of the option
#else
  bool const run_alg_closure = ( _knob_AUTO_ALG_CLOSURE() && _classified );   // rev316: environment-only
#endif
  if( run_alg_closure ){
    size_t const n_alg = _generate_algebraic_boundary_closure();
    if( n_alg ){
      if( options.DISPLAY_LEVEL >= 1 )
        std::cerr << "OCFESLV::setup ** auto_alg_closure: appended " << n_alg
                  << " algebraic-constraint boundary row(s) before the "
                     "consistency/cache build (single pass)\n";
      _phase( "auto_alg_closure (gen)" );
    }
  }


  // Solver: counts, dependency map, interface plan, collocated caches and audits.
  _audit_dof_balance();

  _collect_findings();
  _build_dynamic_form();

  if( !_on_model_ready() ) return false;

  return _issetup;
}

// rev142b: MAXTHREAD -- one cap for every numerical backend.
//
// WHY ONE CAP IS ENOUGH.  SPQR has no thread pool of its own: the TBB path was a
// non-default option upstream and was removed when TBB dropped the feature SPQR used,
// so cc.SPQR_grain and cc.SPQR_nthreads are inert and setting them would be a no-op
// that LOOKS like it works.  What actually threads is the dense BLAS/LAPACK inside the
// frontal factorizations plus CHOLMOD's own OpenMP regions -- and SuperLU, UMFPACK,
// KLU, min2norm and every arma::svd/arma::qr in the audit sit on that same layer.
// Capping it once covers all of them.
//
// WHY dlsym AND NOT A DECLARATION.  MEASURED on the reference build: libgomp resolves
// out of a Python venv's torch/lib, not /usr/lib.  Which runtime is loaded is a
// property of the process, not of the link line, so the setters are found by name at
// run time.  That also means no new link dependency and no clash with <omp.h>.
//
// WHY IT REPORTS.  A cap that silently resolves nothing would let a run be labelled
// single-threaded when it is not, and every timing taken after that would be
// unfalsifiable.  The first application prints what bound; if nothing bound it says so
// loudly.  Nothing is printed at MAXTHREAD=0, so the default path stays silent.
inline FFModel::t_ThreadHooks const&
FFModel::_thread_hooks()
{
  static t_ThreadHooks const H = [](){
    t_ThreadHooks h;
#if defined(__unix__) || defined(__APPLE__)
    // RTLD_DEFAULT searches every object already loaded, however it got there.
    h.omp_set  = reinterpret_cast<void(*)(int)>    ( dlsym( RTLD_DEFAULT, "omp_set_num_threads" ) );
    h.omp_get  = reinterpret_cast<int(*)()>        ( dlsym( RTLD_DEFAULT, "omp_get_max_threads" ) );
    h.blas_set = reinterpret_cast<void(*)(int)>    ( dlsym( RTLD_DEFAULT, "openblas_set_num_threads" ) );
    h.blas_get = reinterpret_cast<int(*)()>        ( dlsym( RTLD_DEFAULT, "openblas_get_num_threads" ) );
    h.mkl_set  = reinterpret_cast<void(*)(int)>    ( dlsym( RTLD_DEFAULT, "mkl_set_num_threads" ) );
    h.mkl_get  = reinterpret_cast<int(*)()>        ( dlsym( RTLD_DEFAULT, "mkl_get_max_threads" ) );
    h.blis_set = reinterpret_cast<void(*)(int64_t)>( dlsym( RTLD_DEFAULT, "bli_thread_set_num_threads" ) );
    h.blis_get = reinterpret_cast<int64_t(*)()>    ( dlsym( RTLD_DEFAULT, "bli_thread_get_num_threads" ) );
#endif
    return h;
  }();
  return H;
}

inline
FFModel::t_ThreadCap::t_ThreadCap
( size_t maxthread )
{
  if( !maxthread ) return;                  // 0 = no cap imposed: make NO call at all
  std::lock_guard<std::mutex> lock( _mutex() );
  _counted = true;
  if( _depth++ ) return;                    // an inner (w-decide re-entry) or concurrent guard: the cap is on
  _outermost = true;

  t_ThreadHooks const& H = _thread_hooks();
  int const cap = (int)maxthread;
  // Capture the pristine values BEFORE anything is changed; restore to these rather
  // than re-reading later, because a getter may be absent while its setter is present.
  _omp0  = H.omp_get  ? H.omp_get()  : 0;
  _blas0 = H.blas_get ? H.blas_get() : 0;
  _mkl0  = H.mkl_get  ? H.mkl_get()  : 0;
  _blis0 = H.blis_get ? H.blis_get() : 0;
  if( H.omp_set  ) H.omp_set ( cap );
  if( H.blas_set ) H.blas_set( cap );
  if( H.mkl_set  ) H.mkl_set ( cap );
  if( H.blis_set ) H.blis_set( (int64_t)cap );
  _active = cap;                            // picked up by the cholmod_common sites

  // Announce once per DISTINCT cap: a run that caps setup() and solve() differently
  // must not have the second one go unreported.
  static int announced = -1;
  if( announced != cap ){
    announced = cap;
    if( H.none() )
      std::cerr << "  [maxthread] ** NO THREAD SETTER RESOLVED ** cap=" << cap
                << " requested and NOT applied (cholmod nthreads_max only)."
                << "  Timings from this build must not be labelled thread-capped."
                << std::endl;
    else{
      std::cerr << "  [maxthread] cap=" << cap << " via";
      if( H.omp_set  ) std::cerr << " omp(was " << _omp0  << ")";
      if( H.blas_set ) std::cerr << " openblas(was " << _blas0 << ")";
      if( H.mkl_set  ) std::cerr << " mkl(was "  << _mkl0  << ")";
      if( H.blis_set ) std::cerr << " blis(was " << _blis0 << ")";
      std::cerr << " + cholmod.nthreads_max; restored on exit" << std::endl;
    }
  }
}

inline
FFModel::t_ThreadCap::~t_ThreadCap()
{
  if( !_counted ) return;
  std::lock_guard<std::mutex> lock( _mutex() );
  if( --_depth > 0 ) return;                // the LAST guard out restores, whichever applied the cap
  t_ThreadHooks const& H = _thread_hooks();
  if( H.omp_set  && _omp0  > 0 ) H.omp_set ( _omp0 );
  if( H.blas_set && _blas0 > 0 ) H.blas_set( _blas0 );
  if( H.mkl_set  && _mkl0  > 0 ) H.mkl_set ( _mkl0 );
  if( H.blis_set && _blis0 > 0 ) H.blis_set( _blis0 );
  _active = 0;
}

inline double
FFModel::_default_dom_ref
( FFVar const& dvar )
const
{
  auto itdom = _mDom.find( dvar );
  if( itdom == _mDom.end() ) throw Exceptions( Exceptions::INDEX );

  double val = 0.5 * ( itdom->second.lo_dom + itdom->second.up_dom );
  for( auto const& [v,x] : _classDomRef )
    if( v.id() == dvar.id() ){ val = x; break; }
  return val;
}

inline double
FFModel::_model_reference_value
( FFVar const& var, t_Coord const& coord, unsigned const depth )
const
{
  // An auxiliary has no declared reference of its own: its value is the derivative of its parent's value in the
  // direction it differentiates.  Taken by central difference at the reference coordinate -- no mesh involved.
  // (The collocated path instead evaluates the defining expression on the collocated arrays, which IS the mesh.)
  if( depth < 8 )
    for( auto const& aux : _auxDef ){
      if( aux.aux.id() != var.id() ) continue;
      double val = _model_reference_value( aux.parent, coord, depth+1 );
      for( auto const& d : aux.diff_dom ){
        auto itd = _mDom.find( d );
        double const span = ( itd != _mDom.end() )? ( itd->second.up_dom - itd->second.lo_dom ): 1.;
        double const h    = 1e-4 * ( std::fabs( span ) > 0.? std::fabs( span ): 1. );
        t_Coord cp = coord, cm = coord;
        cp[d] = coord.count( d )? coord.at( d ) + h: h;
        cm[d] = coord.count( d )? coord.at( d ) - h: -h;
        val = ( _model_reference_value( aux.parent, cp, depth+1 )
              - _model_reference_value( aux.parent, cm, depth+1 ) ) / ( 2.*h );
      }
      return val;
    }

  for( auto const& [v,f] : _classVarRefFun ) if( v.id() == var.id() ) return f( coord );
  for( auto const& [v,f] : _classInpRefFun ) if( v.id() == var.id() ) return f( coord );
  for( auto const& [v,x] : _classVarRef )    if( v.id() == var.id() ) return x;
  for( auto const& [v,x] : _classInpRef )    if( v.id() == var.id() ) return x;
  for( size_t i = 0; i < _vCst.size() && i < _vCstVal.size(); ++i )
    if( _vCst[i].id() == var.id() ) return _vCstVal[i];
  return 0.;
}

inline void
FFModel::_classification_reference
( std::vector<double>& state_ref, std::vector<double>& input_ref,
  std::vector<double>& cst_ref,   std::vector<double>& dom_ref )
{
  // The coordinates at which a reference FUNCTION is evaluated: each domain's reference coordinate.
  t_Coord coord;
  for( auto const& [dvar,dom] : _mDom ){
    (void)dom;
    coord[dvar] = _default_dom_ref( dvar );
  }

  size_t nonfinite = 0;
  auto report = [&nonfinite]( FFVar const& var, double v ) -> double
  {
    if( std::isfinite( v ) ) return v;
    ++nonfinite;
    std::cerr << "OCFESLV::setup [classify] ** reference value for " << var
              << " is " << ( std::isnan( v )? "NaN": "infinite" )
              << "; the classification at this reference point is not trustworthy\n";
    return v;
  };

  state_ref.assign( _vVar.size(), 0. );
  for( size_t i = 0; i < _vVar.size(); ++i )
    state_ref[i] = report( _vVar[i], _model_reference_value( _vVar[i], coord, 0 ) );

  input_ref.assign( _vInp.size(), 0. );
  for( size_t i = 0; i < _vInp.size(); ++i )
    input_ref[i] = report( _vInp[i], _model_reference_value( _vInp[i], coord, 0 ) );

  cst_ref.assign( _vCst.size(), 0. );
  for( size_t i = 0; i < _vCst.size() && i < _vCstVal.size(); ++i )
    cst_ref[i] = _vCstVal[i];

  dom_ref.assign( _vDom.size(), 0. );
  for( size_t i = 0; i < _vDom.size(); ++i )
    dom_ref[i] = _default_dom_ref( _vDom[i] );
}

inline bool
FFModel::_classify_pde
()
{
  std::vector<double> state_ref, input_ref, cst_ref, dom_ref;
  _classification_reference( state_ref, input_ref, cst_ref, dom_ref );
  FFVar const* evolution_dom = _infer_evolution_domain();
  return _classify_pde( state_ref, input_ref, cst_ref, dom_ref, evolution_dom,
                        options.CLASSIFY.NSAMPLE, options.CLASSIFY.IMAG_TOL );
}



inline std::set<FFVar,lt_FFVar>
FFModel::_row_bare_states
( size_t const row )
const
{
  std::set<FFVar,lt_FFVar> out;
  if( row >= _mEqn.size() || !_dag ) return out;
  FFVar e = _mEqn[row].var;
  // Substitute EVERY external-operation node (FFPartial/FFIntegral/FFEval) by a fresh probe leaf; whatever
  // remains a VAR leaf appears in the row outside all of them.
  FFSubgraph sg = _dag->subgraph( 1, &e );
  std::vector<FFVar> nodes;
  for( auto const* op : sg.l_op )
    if( op && ( op->sameid( typeid(FFPartial) ) || op->sameid( typeid(FFIntegral) ) || op->sameid( typeid(FFEval) ) )
        && !op->varout.empty() && op->varout[0] )
      nodes.push_back( *op->varout[0] );
  FFVar red = e;
  for( size_t k = 0; k < nodes.size(); ++k ){
    FFVar pr = _dag->add_var( "_c0bare" + std::to_string( k ) );
    try{ red = _dag->substitute( std::vector<FFVar>{ red }, std::vector<FFVar>{ nodes[k] },
                                 std::vector<FFVar>{ pr } )[0]; }
    catch( ... ){ return out; }
  }
  FFSubgraph sg2 = _dag->subgraph( 1, &red );
  for( auto const* op : sg2.l_op )
    if( op && op->type == FFOp::VAR && op->varout[0] && _mVar.find( *op->varout[0] ) != _mVar.end() )
      out.insert( *op->varout[0] );
  return out;
}

inline std::map<FFVar::pt_idVar,char>
FFModel::_row_invertibility
( size_t const row )
const
{
  std::map<FFVar::pt_idVar,char> out;
  if( row >= _mEqn.size() || !_dag ) return out;
  FFVar const ev = _mEqn[row].var;

  auto letter = []( FFInv::TYPE t ) -> char {
    switch( t ){ case FFInv::TYPE::L: return 'L'; case FFInv::TYPE::S: return 'S';
                 case FFInv::TYPE::N: return 'N'; default: return 'U'; }
  };
  auto leaves_of = [&]( FFVar const& e ){
    std::vector<FFVar> lv;
    FFVar x = e;
    FFSubgraph sg = _dag->subgraph( 1, &x );
    for( auto const* op : sg.l_op )
      if( op && op->type == FFOp::VAR && op->varout[0] ) lv.push_back( *op->varout[0] );
    return lv;
  };

  // (1) FFInv with the STATES seeded: answers L / S / N outright, U when the state reaches the row only
  //     through a derivative (FFPartial/FFIntegral/FFEval are not invertible operations).
  {
    auto leaf = leaves_of( ev );
    std::vector<FFInv> ileaf( leaf.size() ), iout( 1 );
    for( size_t i = 0; i < leaf.size(); ++i )
      if( _mVar.find( leaf[i] ) != _mVar.end() ) ileaf[i].indep( static_cast<int>( leaf[i].id().second ) );
    try{ _dag->eval( std::vector<FFVar>{ ev }, iout, leaf, ileaf ); }
    catch( ... ){ return out; }
    for( auto const& [ind, ty] : iout[0].inv() )
      out[ { FFVar::VAR, static_cast<long>( ind ) } ] = letter( ty );
  }

  // A state FFInv left as U (or did not mention) may still be genuinely nonlinear in the row: FFDep decides,
  // and it propagates through the external operations that make FFInv undetermined.
  {
    auto leaf = leaves_of( ev );
    std::vector<FFDep> dleaf( leaf.size() ), dout( 1 );
    for( size_t i = 0; i < leaf.size(); ++i )
      if( _mVar.find( leaf[i] ) != _mVar.end() ) dleaf[i].indep( static_cast<int>( leaf[i].id().second ) );
    try{
      _dag->eval( std::vector<FFVar>{ ev }, dout, leaf, dleaf );
      for( auto& [vid, ty] : out ){
        if( ty == 'L' || ty == 'S' ) continue;
        auto const dd = dout[0].dep( static_cast<int>( vid.second ) );
        if( dd.first && dd.second != FFDep::TYPE::L ) ty = 'N';
      }
    }
    catch( ... ){}
  }

  // (2) For every first-order bare-state partial in the row, substitute the partial NODE by a fresh probe leaf
  //     and read the coefficient of the DERIVATIVE in FFInv; FFDep on the state says whether the row is linear
  //     in it overall.  The result is marked with a prime: it describes the derivative's coefficient.
  {
    FFVar e = ev;
    FFSubgraph sg = _dag->subgraph( 1, &e );
    std::vector<std::pair<FFVar,FFVar>> dnode;          // (partial output node, its state operand)
    for( auto const* op : sg.l_op ){
      if( !op || !op->sameid( typeid(FFPartial) ) ) continue;
      if( op->varin.empty() || op->varout.empty() || !op->varin[0] || !op->varout[0] ) continue;
      if( _mVar.find( *op->varin[0] ) == _mVar.end() ) continue;
      dnode.emplace_back( *op->varout[0], *op->varin[0] );
    }
    for( size_t k = 0; k < dnode.size(); ++k ){
      auto it = out.find( dnode[k].second.id() );
      if( it != out.end() && it->second != 'U' ) continue;   // already answered directly
      FFVar pr = _dag->add_var( "_c0probe" + std::to_string( k ) );
      FFVar subst;
      try{ subst = _dag->substitute( std::vector<FFVar>{ e },
                                     std::vector<FFVar>{ dnode[k].first },
                                     std::vector<FFVar>{ pr } )[0]; }
      catch( ... ){ continue; }
      auto leaf = leaves_of( subst );
      if( leaf.empty() ) continue;
      // FFInv: seed the probe AND the other states, so a state-valued coefficient shows up as S on the probe.
      std::vector<FFInv> ileaf( leaf.size() ), iout( 1 );
      for( size_t i = 0; i < leaf.size(); ++i )
        if( leaf[i].id() == pr.id() || _mVar.find( leaf[i] ) != _mVar.end() )
          ileaf[i].indep( static_cast<int>( leaf[i].id().second ) );
      char c = 'U';
      try{
        _dag->eval( std::vector<FFVar>{ subst }, iout, leaf, ileaf );
        auto const pd = iout[0].inv( static_cast<int>( pr.id().second ) );
        if( pd.first ) c = letter( pd.second );
      }
      catch( ... ){}
      // FFDep on the state itself: it propagates through the partial, so it answers "linear in s overall?"
      // Run it whenever FFInv did not answer L or S -- in particular for U, which can mean "not invertible
      // under the current option set" (IPOW is not allowed by default, so q*q reads U) rather than nonlinear.
      if( true ){
        auto leaf2 = leaves_of( ev );
        std::vector<FFDep> dleaf( leaf2.size() ), dout( 1 );
        for( size_t i = 0; i < leaf2.size(); ++i )
          if( _mVar.find( leaf2[i] ) != _mVar.end() )
            dleaf[i].indep( static_cast<int>( leaf2[i].id().second ) );
        try{
          _dag->eval( std::vector<FFVar>{ ev }, dout, leaf2, dleaf );
          auto const dd = dout[0].dep( static_cast<int>( dnode[k].second.id().second ) );
          if( dd.first && dd.second != FFDep::TYPE::L ) c = 'N';
        }
        catch( ... ){}
      }
      if( c == 'L' ) c = 'l';        // lower case = answered through the derivative, printed as L'
      else if( c == 'S' ) c = 's';
      out[ dnode[k].second.id() ] = c;
    }
  }
  return out;
}

inline FFModel::t_ClaimMap
FFModel::claim_classification
( int const block_id )
const
{
  t_ClaimMap out;
  if( !_dag || _mEqn.empty() ) return out;   // needs the working model, which exists from the import on

  // ---- per interior row: the states it references, and the (state, direction, order) derivatives it holds ----
  struct Row { size_t idx; int block; std::set<FFVar,lt_FFVar> states;
               std::vector<std::tuple<FFVar,FFVar,unsigned>> derivs; };
  std::vector<Row> rows;
  for( size_t i = 0; i < _mEqn.size(); ++i ){
    auto const& eqn = _mEqn[i];
    if( !eqn.opt->participate_in_classification ) continue;
    if( _ordinary_trace_role( eqn.opt->role ) ) continue;      // face conditions are not pointwise relations
    if( block_id >= 0 && eqn.opt->block_id != block_id ) continue;
    Row r; r.idx = i; r.block = eqn.opt->block_id;
    FFVar const ev = eqn.var;
    FFSubgraph sg = _dag->subgraph( 1, &ev );
    for( auto const* op : sg.l_op ){
      if( !op ) continue;
      if( op->type == FFOp::VAR && op->varout[0] && _mVar.find( *op->varout[0] ) != _mVar.end() )
        r.states.insert( *op->varout[0] );
      if( !op->sameid( typeid(FFPartial) ) ) continue;
      auto const* pop = mc::type_cast<FFPartial const>( op );
      if( !pop ) continue;
      for( auto const* operand : op->varin ){
        if( !operand || _mVar.find( *operand ) == _mVar.end() ) continue;
        for( auto const& [dv, ord] : pop->Indep().expr )
          if( ord && _mDom.find( dv ) != _mDom.end() )
            r.derivs.emplace_back( *operand, dv, static_cast<unsigned>( ord ) );
      }
    }
    rows.push_back( std::move( r ) );
  }

  std::set<int> blocks;
  for( auto const& r : rows ) blocks.insert( r.block );

  for( int b : blocks ){
    std::set<FFVar,lt_FFVar> bstates;
    for( auto const& r : rows ) if( r.block == b ) bstates.insert( r.states.begin(), r.states.end() );

    for( auto const& [dvar, dom] : _mDom ){ (void)dom;
      // only states that DEPEND on this direction take part: a claim in a direction the state does not vary
      // in is meaningless
      std::set<FFVar,lt_FFVar> dstates;
      for( auto const& s : bstates ){
        auto iv = _mVar.find( s );
        if( iv != _mVar.end() && iv->second.find( dvar ) != iv->second.end() ) dstates.insert( s );
      }
      // Orders and implications depend on each other: a state is C1 in d if the state standing for its
      // d-derivative is CONTINUOUS in d -- natural, or implied through an admissible relation -- and which
      // relations are admissible depends on the orders.  Iterate to a fixpoint (monotone: the continuous set and
      // the admissible set only grow).
      std::set<FFVar,lt_FFVar> continuous;
      std::map<FFVar,t_Claim,lt_FFVar> claims;
      for( unsigned pass = 0; pass < 8; ++pass ){
      claims.clear();
        // NATURAL(d): differentiated in d in some interior row of the block; record the justifying row.
        std::map<FFVar,size_t,lt_FFVar> natural;
        for( auto const& r : rows ){
          if( r.block != b ) continue;
          for( auto const& [s, dd, k] : r.derivs )
            if( dd.id() == dvar.id() && natural.find( s ) == natural.end() ) natural[s] = r.idx;
        }
        // continuity order: C1 iff a lifted d-derivative of s is itself NATURAL(d)
        // continuity order: C1 iff some state a standing for ds/dd is itself NATURAL(d).  Two ways a can arise:
        // a reduction-minted auxiliary (recorded in _auxDef), or a USER-declared flux state defined by a row of
        // the form a - ds/dd (exactly one derivative, exactly one other state) -- the mixed-form heat equation
        // q = T_x is the plain case.
        auto order_of = [&]( FFVar const& s ) -> int {
          if( natural.find( s ) == natural.end() ) return -1;
          for( auto const& aux : _auxDef ){
            if( aux.parent.id() != s.id() || natural.find( aux.aux ) == natural.end() ) continue;
            for( auto const& d2 : aux.diff_dom ) if( d2.id() == dvar.id() ) return 1;
          }
          // a row with exactly one derivative, ds/dd, standing for "a = (state-free factor) ds/dd": in FFInv
          // with only the states seeded, a is L and NO state other than s and a carries a type -- a state
          // factor M would show up as S (it multiplies the partial).
          for( auto const& r : rows ){
            if( r.block != b || r.derivs.size() != 1 ) continue;
            auto const& [ds, dd, k] = r.derivs[0];
            if( ds.id() != s.id() || dd.id() != dvar.id() || k != 1 ) continue;
            auto const inv = _row_invertibility( r.idx );
            for( auto const& a : r.states ){
              if( a.id() == s.id() ) continue;
              if( natural.find( a ) == natural.end() && continuous.find( a ) == continuous.end() ) continue;
              auto ia = inv.find( a.id() );
              if( ia == inv.end() || ( ia->second != 'L' && ia->second != 'l' ) ) continue;
              bool clean = true;
              for( auto const& [vid, ty] : inv )
                if( vid != a.id() && vid != s.id() && ty != 'U' ){ clean = false; break; }
              if( clean ) return 1;
            }
          }
          return 0;
        };
        // ADMISSIBLE in d: every d-derivative is of a state continuous to that order
        auto admissible = [&]( Row const& r ) -> bool {
          for( auto const& [s, dd, k] : r.derivs )
            if( dd.id() == dvar.id() && order_of( s ) < static_cast<int>( k ) ) return false;
          return true;
        };

        std::set<FFVar,lt_FFVar> known;          // NATURAL states; inputs are known by construction
        for( auto const& [s,i] : natural ) known.insert( s );
        // `claims` is declared outside the fixpoint loop and cleared per pass
        for( auto const& s : dstates ){
          t_Claim c; c.algebraic = true;
          for( auto const& r : rows ) if( r.block == b )
            for( auto const& [s2, dd, k] : r.derivs )
              if( s2.id() == s.id() && dd.id() == _evolution_dom_var.id() && _evolution_dom_set ) c.algebraic = false;
          auto it = natural.find( s );
          if( it != natural.end() ){ c.tag = ClaimTag::NATURAL; c.order = order_of( s ); c.via = it->second; }
          claims[s] = c;
        }

        // IMPLIED(d): fixpoint.  Singleton relations first (each determines exactly one remaining unknown) ...
        std::vector<size_t> adm;
        for( size_t ri = 0; ri < rows.size(); ++ri )
          if( rows[ri].block == b && admissible( rows[ri] ) ) adm.push_back( ri );
        std::set<size_t> used;
        auto unknowns_of = [&]( Row const& r ){
          // only a BARE occurrence can be determined pointwise by this row
          auto const bare = _row_bare_states( r.idx );
          std::vector<FFVar> u;
          for( auto const& s : r.states )
            if( known.find( s ) == known.end() && bare.find( s ) != bare.end() ) u.push_back( s );
          return u;
        };
        auto imply = [&]( FFVar const& u, Row const& r ){
          t_Claim& c = claims[u];
          c.tag = ClaimTag::IMPLIED; c.order = 0; c.via = r.idx; c.rests_on.clear();
          for( auto const& s : r.states ) if( s.id() != u.id() ) c.rests_on.push_back( s );
          known.insert( u ); used.insert( r.idx );
        };
        for( bool changed = true; changed; ){
          changed = false;
          for( size_t ri : adm ){
            Row const& r = rows[ri];
            if( used.count( r.idx ) ) continue;
            auto u = unknowns_of( r );
            if( u.size() == 1 ){ imply( u[0], r ); changed = true; }
          }
        }
        // ... then a maximum matching over what is left, for COUPLED closure sets: a matched unknown whose
        // relation's every unknown is matched too belongs to a square, structurally nonsingular subsystem
        // resting only on known states.
        {
          std::vector<size_t> R; std::vector<FFVar> U;
          for( size_t ri : adm ) if( !used.count( rows[ri].idx ) && !unknowns_of( rows[ri] ).empty() ) R.push_back( ri );
          std::set<FFVar,lt_FFVar> Uset;
          for( size_t ri : R ) for( auto const& s : unknowns_of( rows[ri] ) ) Uset.insert( s );
          U.assign( Uset.begin(), Uset.end() );
          if( !R.empty() && !U.empty() ){
            std::map<FFVar::pt_idVar,size_t> uidx; for( size_t j = 0; j < U.size(); ++j ) uidx[U[j].id()] = j;
            std::vector<std::vector<size_t>> adj( R.size() );
            for( size_t i = 0; i < R.size(); ++i )
              for( auto const& s : unknowns_of( rows[R[i]] ) ) adj[i].push_back( uidx[s.id()] );
            std::vector<long> matchU( U.size(), -1 ), matchR( R.size(), -1 );
            std::function<bool(size_t,std::vector<char>&)> aug = [&]( size_t i, std::vector<char>& seen ){
              for( size_t j : adj[i] ){
                if( seen[j] ) continue; seen[j] = 1;
                if( matchU[j] < 0 || aug( static_cast<size_t>( matchU[j] ), seen ) ){ matchU[j] = i; matchR[i] = j; return true; }
              }
              return false;
            };
            for( size_t i = 0; i < R.size(); ++i ){ std::vector<char> seen( U.size(), 0 ); aug( i, seen ); }
            // closure: drop matched pairs whose relation still has an unmatched unknown, until stable
            std::vector<char> ok( R.size(), 0 );
            for( size_t i = 0; i < R.size(); ++i ) ok[i] = ( matchR[i] >= 0 );
            for( bool changed = true; changed; ){
              changed = false;
              for( size_t i = 0; i < R.size(); ++i ){
                if( !ok[i] ) continue;
                for( size_t j : adj[i] )
                  if( matchU[j] < 0 || !ok[static_cast<size_t>( matchU[j] )] ){ ok[i] = 0; changed = true; break; }
              }
            }
            for( size_t i = 0; i < R.size(); ++i )
              if( ok[i] ) imply( U[static_cast<size_t>( matchR[i] )], rows[R[i]] );
          }
        }

        // CLOSURE(d): a sub-case of IMPLIED, keyed on the RELATION -- it is CLOSURE when the determining relation
        // is not affine in the states with constant coefficients (any letter other than L on any state or probed
        // derivative).  The verdict is the same (no claim); the tag tells C2/C3 that no constant coupling exists.
        // The letter of s itself is kept: C1 separately needs "can the relation be solved for s linearly".
        for( auto& [s, c] : claims ){
          if( c.tag != ClaimTag::IMPLIED ) continue;
          auto const inv = _row_invertibility( c.via );
          auto it = inv.find( s.id() );
          c.coeff = ( it == inv.end() )? 'U': it->second;
          bool affine_constant = true;
          for( auto const& [vid, ty] : inv )
            if( ty != 'L' && ty != 'l' ){ affine_constant = false; break; }
          if( !affine_constant ) c.tag = ClaimTag::CLOSURE;
        }
        std::set<FFVar,lt_FFVar> now;
        for( auto const& [s, c] : claims )
          if( c.tag != ClaimTag::NONE ) now.insert( s );
        if( now == continuous ) break;
        continuous = now;
      }

      for( auto const& [s, c] : claims ) out[s][dvar] = c;
    }
  }
  return out;
}

inline void
FFModel::claim_report
( std::ostream& os, int const block_id )
const
{
  auto const cm = claim_classification( block_id );
  auto name = [&]( FFVar const& v ){ std::ostringstream o; o << v; return o.str(); };
  auto tag  = []( ClaimTag t ){ switch( t ){ case ClaimTag::NATURAL: return "NATURAL"; case ClaimTag::IMPLIED: return "IMPLIED";
                                             case ClaimTag::CLOSURE: return "CLOSURE"; default: return "NONE"; } };
  os << "OCFESLV::claim ** continuity-claim classification (C0, report only): " << cm.size() << " state(s)\n";
  for( auto const& [s, byd] : cm ){
    for( auto const& [d, c] : byd ){
      os << "  [claim] " << std::left << std::setw(12) << name( s ) << " in " << std::setw(6) << name( d )
         << ( c.algebraic? " alg " : " dyn " ) << std::setw(8) << tag( c.tag ) << std::right
         << " order=" << std::setw(2) << c.order;
      if( c.via != size_t(-1) ) os << "  via row " << c.via;
      if( c.coeff != '-' ){
        // lower case marks an answer obtained through the DERIVATIVE: the letter is the coefficient of ds/dd
        char const up = ( c.coeff >= 'a' && c.coeff <= 'z' )? static_cast<char>( c.coeff - 32 ): c.coeff;
        os << " (" << up << ( up != c.coeff? "'": "" ) << ")";
      }
      if( !c.rests_on.empty() ){
        os << "  rests on {";
        for( size_t i = 0; i < c.rests_on.size(); ++i ) os << ( i? ",": "" ) << name( c.rests_on[i] );
        os << "}";
      }
      os << "\n";
    }
  }
}

// 2026-09-28: over-alignment guard (GCC wrong-code with over-aligned virtual bases)
static_assert( alignof( FFModel ) <= alignof( void* ), "FFModel is a VIRTUAL base of the CRONOS solvers and must not be over-aligned (alignof > 8): GCC (6 to at least 16) emits aligned vector stores in base-object constructors assuming the full alignment, but a virtual-base subobject is only placed at its non-virtual alignment -> SIGSEGV at -O2/-O3 (see gccvb/pr_vbase_align.cpp). Keep over-aligned members (Armadillo/Eigen fixed-size types, alignas) behind a pointer, as FFModel::_pClassification does." );

} // namespace mc

#endif