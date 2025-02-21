#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>

#include "ffode.hpp"
//#include "odeslvs_cvodes.hpp" 

namespace py = pybind11;

void mc_odeslvs( py::module_ &m )
{

typedef mc::FFGraph FFGraph;
typedef mc::ODESLVS_CVODES ODESLVS;

py::class_<ODESLVS> pyODESLV( m, "ODESLV" );
pyODESLV
 .def(
   py::init<>()
 )
 .def_readwrite( 
   "options",
   &ODESLVS::options
 )
 .def( 
   "set",
   []( ODESLVS& self,  ODESLVS& ivp ){ self.set( ivp ); },
   py::keep_alive<1,2>(),
   "copy initial value problem"
 )
 .def( 
   "set_dag",
   []( ODESLVS& self,  FFGraph* dag ){ self.set_dag( dag ); },
   py::keep_alive<1,2>(),
   "set DAG"
 )
 .def(
   "dag",
   []( ODESLVS& self ){ return self.dag(); },
   py::return_value_policy::reference_internal,
   "get DAG"
 )
 .def_property_readonly(
   "ns",
   []( ODESLVS const& self ){ return self.nsmax(); },
   "number of stages"
 )
 .def(
   "set_time",
   []( ODESLVS& self, std::vector<double> const& ts, mc::FFVar const* t ){ self.set_time( ts, t ); },
   py::arg("stages"),
   py::arg("time")=nullptr,
   "set parameters"
 )
 .def_property_readonly(
   "var_time",
   []( ODESLVS const& self ){ return self.var_time()? std::vector<mc::FFVar>({*self.var_time()}): std::vector<mc::FFVar>(); },
   py::return_value_policy::reference_internal,
   "time variable"
 )
 .def_property_readonly(
   "val_stage",
   []( ODESLVS const& self ){ return self.val_stage(); },
   py::return_value_policy::reference_internal,
   "time stages"
 )
 .def_property_readonly(
   "nc",
   []( ODESLVS const& self ){ return self.nc(); },
   "number of constants"
 )
 .def(
   "set_constant",
   []( ODESLVS& self, std::vector<mc::FFVar> const& c ){ self.set_constant( c ); },
   "set constants"
 )
 .def_property_readonly(
   "var_constant",
   []( ODESLVS const& self ){ return self.var_constant(); },
   py::return_value_policy::reference_internal,
   "constants"
 )
 .def_property_readonly(
   "np",
   []( ODESLVS const& self ){ return self.np(); },
   "number of parameters"
 )
 .def(
   "set_parameter",
   []( ODESLVS& self, std::vector<mc::FFVar> const& p ){ self.set_parameter( p ); },
   "set parameters"
 )
 .def_property_readonly(
   "var_parameter",
   []( ODESLVS const& self ){ return self.var_parameter(); },
   py::return_value_policy::reference_internal,
   "parameters"
 )
 .def_property_readonly(
   "nx",
   []( ODESLVS const& self ){ return self.nx(); },
   "number of states"
 )
 .def(
   "set_state",
   []( ODESLVS& self, std::vector<mc::FFVar> const& x ){ self.set_state( x ); },
   "set state variables"
 )
 .def_property_readonly(
   "var_state",
   []( ODESLVS const& self ){ return self.var_state(); },
   py::return_value_policy::reference_internal,
   "state variables"
 )
 .def(
   "set_differential",
   []( ODESLVS& self, std::vector<mc::FFVar> const& rhs ){ self.set_differential( rhs ); },
   "set differential equations"
 )
 .def(
   "set_differential",
   []( ODESLVS& self, std::vector<std::vector<mc::FFVar>> const& rhs ){ self.set_differential( rhs ); },
   "set differential equations"
 )
 .def_property_readonly(
   "eqn_differential",
   []( ODESLVS const& self ){ return self.eqn_differential(); },
   py::return_value_policy::reference_internal,
   "differential equations"
 )
 .def_property_readonly(
   "nq",
   []( ODESLVS const& self ){ return self.nq(); },
   "number of quadratures"
 )
 .def(
   "set_quadrature",
   []( ODESLVS& self, std::vector<mc::FFVar> const& quad, std::vector<mc::FFVar> const& q ){ self.set_quadrature( quad, q ); },
   "set quadrature variables and equations"
 )
 .def(
   "set_quadrature",
   []( ODESLVS& self, std::vector<std::vector<mc::FFVar>> const& quad, std::vector<mc::FFVar> const& q ){ self.set_quadrature( quad, q ); },
   "set quadrature variables and equations"
 )
 .def_property_readonly(
   "var_quadrature",
   []( ODESLVS const& self ){ return self.var_quadrature(); },
   py::return_value_policy::reference_internal,
   "quadrature variables"
 )
 .def_property_readonly(
   "eqn_quadrature",
   []( ODESLVS const& self ){ return self.eqn_quadrature(); },
   py::return_value_policy::reference_internal,
   "quadrature equations"
 )
 .def(
   "set_initial",
   []( ODESLVS& self, std::vector<mc::FFVar> const& ic ){ self.set_initial( ic ); },
   "set initial conditions"
 )
 .def(
   "set_initial",
   []( ODESLVS& self, std::vector<std::vector<mc::FFVar>> const& ic ){ self.set_initial( ic ); },
   "set initial conditions"
 )
 .def_property_readonly(
   "eqn_initial",
   []( ODESLVS const& self ){ return self.eqn_initial(); },
   py::return_value_policy::reference_internal,
   "initial conditions"
 )
 .def_property_readonly(
   "nf",
   []( ODESLVS const& self ){ return self.nf(); },
   "number of functions"
 )
 .def(
   "set_function",
   []( ODESLVS& self, std::vector<mc::FFVar> const& fct ){ self.set_function( fct ); },
   "set state functions"
 )
 .def(
   "set_function",
   []( ODESLVS& self, std::vector<std::vector<mc::FFVar>> const& fct ){ self.set_function( fct ); },
   "set state functions"
 )
 .def_property_readonly(
   "eqn_function",
   []( ODESLVS const& self ){ return self.eqn_function(); },
   py::return_value_policy::reference_internal,
   "state functions"
 )
 .def(
   "setup",
   []( ODESLVS& self ){ self.setup(); },
   "setup initial value problem before solve"
 )
 .def(
   "fdiff",
   []( ODESLVS& self, std::vector<mc::FFVar> const& p ){ return self.fdiff( p ); },
   py::return_value_policy::take_ownership,
   "differentiate initial value problem with respect to selected parameters"
 )
 .def(
   "solve_state",
   []( ODESLVS& self, std::vector<double> const& p, std::vector<double> const& c ){ return self.solve_state( p, c ); },
   py::arg("p"), py::arg("c") = std::vector<double>(),
   "solve initial value problem"
 )
 .def(
   "solve_sensitivity",
   []( ODESLVS& self, std::vector<double> const& p, std::vector<double> const& c ){ return self.solve_sensitivity( p, c ); },
   py::arg("p"), py::arg("c") = std::vector<double>(),
   "solve initial value problem with forward sensitivity"
 )
 .def(
   "solve_adjoint",
   []( ODESLVS& self, std::vector<double> const& p, std::vector<double> const& c ){ return self.solve_adjoint( p, c ); },
   py::arg("p"), py::arg("c") = std::vector<double>(),
   "solve initial value problem with adjoint sensitivity"
 )
 .def_property_readonly(
   "final_time",
   []( ODESLVS const& self ){ return self.final_time(); },
   py::return_value_policy::reference_internal,
   "final time reached during solve"
 )
 .def_property_readonly(
   "final_stage",
   []( ODESLVS const& self ){ return self.final_stage(); },
   py::return_value_policy::reference_internal,
   "final stage reached during solve"
 )
 .def_property_readonly(
   "val_state",
   []( ODESLVS const& self ){ return self.val_state(); },
   py::return_value_policy::reference_internal,
   "state values"
 )
 .def_property_readonly(
   "val_quadrature",
   []( ODESLVS const& self ){ return self.val_quadrature(); },
   py::return_value_policy::reference_internal,
   "quadrature values"
 )
 .def_property_readonly(
   "val_function",
   []( ODESLVS const& self ){ return self.val_function(); },
   py::return_value_policy::reference_internal,
   "function values"
 )
 .def_property_readonly(
   "val_state_sensitivity",
   []( ODESLVS const& self ){ return self.val_state_sensitivity(); },
   py::return_value_policy::reference_internal,
   "state sensitivities"
 )
 .def_property_readonly(
   "val_state_adjoint",
   []( ODESLVS const& self ){ return self.val_state_adjoint(); },
   py::return_value_policy::reference_internal,
   "state adjoints"
 )
 .def_property_readonly(
   "val_quadrature_sensitivity",
   []( ODESLVS const& self ){ return self.val_quadrature_sensitivity(); },
   py::return_value_policy::reference_internal,
   "quadrature sensitivities"
 )
 .def_property_readonly(
   "val_function_gradient",
   []( ODESLVS const& self ){ return self.val_function_gradient(); },
   py::return_value_policy::reference_internal,
   "function gradient"
 )
 .def_readonly( 
   "results_state",
   &ODESLVS::results_state
 )
 .def_readonly( 
   "results_sensitivity",
   &ODESLVS::results_sensitivity
 )
;

py::enum_<ODESLVS::STATUS>(pyODESLV, "Status", py::module_local())
 .value("Normal",  ODESLVS::STATUS::NORMAL )
 .value("Failure", ODESLVS::STATUS::FAILURE)
 .value("Fatal",   ODESLVS::STATUS::FATAL  )
 .export_values()
;

py::class_<ODESLVS::Results> pyODESLVResults( pyODESLV, "Results" );

pyODESLVResults
 .def_readonly( "t", &ODESLVS::Results::t, "time" )
 .def_readonly( "x", &ODESLVS::Results::x, "values" )
 .def( "__str__",
       []( ODESLVS::Results& self ){
         std::ostringstream ss;
         ss << self.t << ": [";
         for( auto const& xi : self.x ) ss << " " << xi;
         ss << "]";
         return ss.str();
       }
 )
 .def( "__repr__",
       []( ODESLVS::Results& self ){
         std::ostringstream ss;
         ss << self.t << ": [";
         for( auto const& xi : self.x ) ss << " " << xi;
         ss << "]";
         return ss.str();
       }
 )
;

py::class_<ODESLVS::Options> pyODESLVOptions( pyODESLV, "Options" );//, py::module_local() );

pyODESLVOptions
 .def( py::init<>() )
 .def( py::init<ODESLVS::Options const&>() )
 .def_readwrite( "DISPLEVEL", &ODESLVS::Options::DISPLAY,   "display level [Default: 1]" )
 .def_readwrite( "RESRECORD", &ODESLVS::Options::RESRECORD, "record level for simulation results [Default: 0]" )
 .def_readwrite( "INTMETH",   &ODESLVS::Options::INTMETH,   "numerical integration method [Default: MSBDF]" )
 .def_readwrite( "NLINSOL",   &ODESLVS::Options::NLINSOL,   "Nonlinear solver method [Default: NEWTON]" )
 .def_readwrite( "LINSOL",    &ODESLVS::Options::LINSOL,    "Linear solver method and Jacobian approximation [Default: DIAG]" )
 .def_readwrite( "H0",        &ODESLVS::Options::H0,        "Initial step-size [Default: 0e0 (auto)]" )
 .def_readwrite( "HMIN",      &ODESLVS::Options::HMIN,      "Minimum step-size [Default: 0e0]" )
 .def_readwrite( "HMAX",      &ODESLVS::Options::HMAX,      "Maximum step-size [Default: 0e0 (+inf)]" )
 .def_readwrite( "NMAX",      &ODESLVS::Options::NMAX,      "Maximum number of steps in a time stage [Default: 2000]" )
 .def_readwrite( "MAXFAIL",   &ODESLVS::Options::MAXFAIL,   "Maximum number of error test failures per step [Default: 10]" )
 .def_readwrite( "MAXCORR",   &ODESLVS::Options::MAXCORR,   "Maximum  number of nonlinear solver iterations per step [Default: 5]" )
 .def_readwrite( "RTOL",      &ODESLVS::Options::RTOL,      "Relative integration tolerance [Default: 1e-7]" )
 .def_readwrite( "ATOL",      &ODESLVS::Options::ATOL,      "Absolute integration tolerance [Default: 1e-9]" )
 .def_readwrite( "QERR",      &ODESLVS::Options::QERR,      "Whether to perform error control on quadrature variables [Default: 1]" )
 .def_readwrite( "RTOLS",     &ODESLVS::Options::RTOLS,     "Relative integration tolerance for FSA [Default: 1e-7]" )
 .def_readwrite( "ATOLS",     &ODESLVS::Options::ATOLS,     "Absolute integration tolerance for FSA [Default: 1e-9]" )
 .def_readwrite( "AUTOTOLS",  &ODESLVS::Options::AUTOTOLS,  "Whether to set automatic integration tolerances for FSA [Default: 0]" )
 .def_readwrite( "QERRS",     &ODESLVS::Options::QERRS,     "Whether to perform error control on quadrature sensitivities [Default: 1]" )
 .def_readwrite( "FSAERR",    &ODESLVS::Options::FSAERR,    "Whether to perform error control on sensitivities [Default: 1]" )
 .def_readwrite( "FSACORR",   &ODESLVS::Options::FSACORR,   "FSA correction strategy [Default: STAGGERED]" )
 .def_readwrite( "RTOLB",     &ODESLVS::Options::RTOLB,     "Relative integration tolerance for ASA [Default: 1e-7]" )
 .def_readwrite( "ATOLB",     &ODESLVS::Options::ATOLB,     "Absolute integration tolerance for ASA [Default: 1e-9]" )
 .def_readwrite( "QERRB",     &ODESLVS::Options::QERRB,     "Whether to perform error control on quadrature adjoints [Default: 1]" )
 .def_readwrite( "ASAINTERP", &ODESLVS::Options::ASAINTERP, "ASA interpolation strategy [Default: HERMITE]" )
 .def_readwrite( "ASACHKPT",  &ODESLVS::Options::ASACHKPT,  "Number of steps between each check point for ASA [Default: 2000]" )
;

py::enum_<ODESLVS::Options::INTEGRATION_METHOD>(pyODESLVOptions, "INTEGRATION_METHOD")
 .value("MSADAMS", ODESLVS::Options::INTEGRATION_METHOD::MSADAMS, "Variable-coefficient linear multistep Adams method (non-stiff systems)")
 .value("MSBDF",   ODESLVS::Options::INTEGRATION_METHOD::MSBDF,   "Variable-coefficient linear multistep backward differentiation formula (BDF) method (stiff systems)")
 .export_values()
;

py::enum_<ODESLVS::Options::FSA_STRATEGY>(pyODESLVOptions, "FSA_STRATEGY")
 .value("SIMULTANEOUS", ODESLVS::Options::FSA_STRATEGY::SIMULTANEOUS, "Simultaneous state/sensitivity correction")
 .value("STAGGERED",    ODESLVS::Options::FSA_STRATEGY::STAGGERED,    "Simultaneous sensitivity corrections after state corrections")
 .value("STAGGERED1",   ODESLVS::Options::FSA_STRATEGY::STAGGERED1,   "Sequential sensitivity corrections after state corrections")
 .export_values()
;

py::enum_<ODESLVS::Options::ASA_STRATEGY>(pyODESLVOptions, "ASA_STRATEGY")
 .value("HERMITE",    ODESLVS::Options::ASA_STRATEGY::HERMITE,    "Cubic Hermite interpolation")
 .value("POLYNOMIAL", ODESLVS::Options::ASA_STRATEGY::POLYNOMIAL, "Variable degree polynomial interpolation")
 .export_values()
;

py::enum_<ODESLVS::Options::NONLINEAR_SOLVER>(pyODESLVOptions, "NONLINEAR_SOLVER")
 .value("FIXEDPOINT", ODESLVS::Options::NONLINEAR_SOLVER::FIXEDPOINT, "Fixed point nonlinear solver")
 .value("NEWTON",     ODESLVS::Options::NONLINEAR_SOLVER::NEWTON,     "Newton nonlinear solver")
 .export_values()
;

py::enum_<ODESLVS::Options::LINEAR_SOLVER>(pyODESLVOptions, "LINEAR_SOLVER")
 .value("DIAG",    ODESLVS::Options::LINEAR_SOLVER::DIAG,    "Approximate diagonal Jacobian formed by way of a difference quotient")
 .value("DENSE",   ODESLVS::Options::LINEAR_SOLVER::DENSE,   "Use analytic dense Jacobian and internal direct dense linear algebra functions")
 .value("DENSEDQ", ODESLVS::Options::LINEAR_SOLVER::DENSEDQ, "Use approximate dense Jacobian by way of a difference quotient and internal direct dense linear algebra functions")
#if defined( CRONOS__WITH_KLU )
 .value("SPARSE", ODESLVS::Options::LINEAR_SOLVER::SPARSE, "Use analytic sparse Jacobian and KLU for the direct solution of sparse nonsymmetric linear systems of equations")
#endif
 .export_values()
;
}

void mc_ffode( py::module_ &m )
{

typedef mc::ODESLVS_CVODES ODESLVS;

py::class_<mc::FFODE, mc::FFOp> pyFFODE( m, "FFODE" );//, py::module_local() );

pyFFODE
 .def(
   py::init<>(),
   "default constructor"
 )
 .def_readwrite_static(
   "options",
   &mc::FFODE::options
 )
 .def(
   "__call__",
   []( mc::FFODE& self, std::vector<mc::FFVar> const& vVar, ODESLVS* pODE )
   {
     auto pDep = self( vVar.size(), vVar.data(), pODE );
     return std::vector<mc::FFVar*>( pDep, pDep+pODE->nf() );
   },
   py::return_value_policy::reference_internal,
   "define ODE operation in DAG"
 )
 .def(
   "__call__",
   []( mc::FFODE& self, std::vector<mc::FFVar> const& vVar, std::vector<mc::FFVar> vCst, ODESLVS* pODE )
   {
     auto pDep = self( vVar.size(), vVar.data(), vCst.size(), vCst.data(), pODE );
     return std::vector<mc::FFVar*>( pDep, pDep+pODE->nf() );
   },
   py::return_value_policy::reference_internal,
   "define ODE operation in DAG"
 )
 .def(
   "__call__",
   []( mc::FFODE& self, unsigned const idep, std::vector<mc::FFVar> const& vVar, ODESLVS* pODE )
   {
     return self( idep, vVar.size(), vVar.data(), pODE );
   },
   py::return_value_policy::reference_internal,
   "define ODE operation in DAG"
 )
 .def(
   "__call__",
   []( mc::FFODE& self, unsigned const idep, std::vector<mc::FFVar> const& vVar, std::vector<mc::FFVar> vCst, ODESLVS* pODE )
   {
     return self( idep, vVar.size(), vVar.data(), vCst.size(), vCst.data(), pODE );
   },
   py::return_value_policy::reference_internal,
   "define ODE operation in DAG"
 )
// .def_property(
//   "options",
//   []( mc::FFODE& self )
//   {
//     return self.pODESLV()->options;
//   },
//   []( mc::FFODE& self, ODESLVS::Options& options )
//   {
//     self.pODESLV()->options = options;
//   },
//   py::return_value_policy::reference_internal,
//   "ODE options"
// )
// .def(
//   "setup",
//   []( mc::FFODE& self )
//   {
//     ODESLVS* data = static_cast<ODESLVS*>( self.data );
//     data->setup();
//   },
//   "ODE setup"
// )
// .def_property_readonly(
//   "ODESLV",
//   []( mc::FFODE& self ) -> mc::ODESLVS_CVODES*
//   {
//     return self.pODESLV();
//   },
//   "ODE object"
// )
 .def_readwrite(
   "type",
   &mc::FFODE::type,
   "retreive operation type"
 )
 .def_readwrite(
   "info",
   &mc::FFODE::info,
   "retreive operation id"
 )
 .def_readwrite(
   "varin",
   &mc::FFODE::varin
 )
 .def_readwrite(
   "varout",
   &mc::FFODE::varout
 )
 .def( "name",
   &mc::FFODE::name,
   "retreive operation name"
 )
 .def(
   "__str__",
   []( mc::FFODE const& O )
   {
     std::ostringstream Oss;
     Oss << O;
     return Oss.str();
   }
 )
 .def(
   "__repr__",
   []( mc::FFODE const& O )
   {
     std::ostringstream Oss;
     Oss << O;
     return Oss.str();
   }
 )
;

py::class_<mc::FFODE::Options>  pyFFODEOptions( pyFFODE, "Options" );

pyFFODEOptions
 .def( py::init<>() )
 .def( py::init<mc::FFODE::Options const&>() )
 .def_readwrite( "DIFF", &mc::FFODE::Options::DIFF,   "method of ODE differentiation [Default: NUM_P]" )
 .def_readwrite( "NP2NF", &mc::FFODE::Options::NP2NF, "parameter-to-function-size ratio above which adjoint sensitivity is applied instead of forward sensitivity [Default: 3]" )
;

py::enum_<mc::FFODE::Options::DERIV_TYPE>( pyFFODEOptions, "DERIV_TYPE" )
 .value("NUM_P", mc::FFODE::Options::DERIV_TYPE::NUM_P, "Derivatives w.r.t. parameters only through forward or adjoint sensitivity integration")
 .value("SYM_P", mc::FFODE::Options::DERIV_TYPE::SYM_P, "Derivatives w.r.t. parameters only through symbolic differentiation of ODEs")
 .value("SYM_C", mc::FFODE::Options::DERIV_TYPE::SYM_C, "Derivatives w.r.t. constants only through symbolic differentiation of ODEs")
 .value("SYM_PC", mc::FFODE::Options::DERIV_TYPE::SYM_PC, "Derivatives w.r.t. parameters and constants jointly through symbolic differentiation of ODEs")
 .export_values()
;
}

