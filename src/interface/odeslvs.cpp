#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>

#include "odeslvs_cvodes.hpp" 

namespace py = pybind11;

void mc_odeslvs( py::module_ &m )
{

typedef mc::FFGraph<> FFGraph;
typedef mc::ODESLVS_CVODES<> ODESLVS;

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
   "solve_state",
   []( ODESLVS& self, std::vector<double> const& p ){ return self.solve_state( p ); },
   "solve initial value problem"
 )
 .def(
   "solve_sensitivity",
   []( ODESLVS& self, std::vector<double> const& p ){ return self.solve_sensitivity( p ); },
   "solve initial value problem with forward sensitivity"
 )
 .def(
   "solve_adjoint",
   []( ODESLVS& self, std::vector<double> const& p ){ return self.solve_adjoint( p ); },
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

py::enum_<ODESLVS::STATUS>(pyODESLV, "ODESLV.STATUS", py::module_local())
 .value("NORMAL",  ODESLVS::STATUS::NORMAL )
 .value("FAILURE", ODESLVS::STATUS::FAILURE)
 .value("FATAL",   ODESLVS::STATUS::FATAL  )
 .export_values()
;

py::class_<ODESLVS::Results> pyODESLVResults( pyODESLV, "ODESLV.Results" );
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

py::class_<ODESLVS::Options> pyODESLVOptions( pyODESLV, "ODESLV.Options", py::module_local() );
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

py::enum_<ODESLVS::Options::INTEGRATION_METHOD>(pyODESLVOptions, "ODESLV.INTEGRATION_METHOD")
 .value("MSADAMS", ODESLVS::Options::INTEGRATION_METHOD::MSADAMS, "Variable-coefficient linear multistep Adams method (non-stiff systems)")
 .value("MSBDF",   ODESLVS::Options::INTEGRATION_METHOD::MSBDF,   "Variable-coefficient linear multistep backward differentiation formula (BDF) method (stiff systems)")
 .export_values()
;

py::enum_<ODESLVS::Options::FSA_STRATEGY>(pyODESLVOptions, "ODESLV.FSA_STRATEGY")
 .value("SIMULTANEOUS", ODESLVS::Options::FSA_STRATEGY::SIMULTANEOUS, "Simultaneous state/sensitivity correction")
 .value("STAGGERED",    ODESLVS::Options::FSA_STRATEGY::STAGGERED,    "Simultaneous sensitivity corrections after state corrections")
 .value("STAGGERED1",   ODESLVS::Options::FSA_STRATEGY::STAGGERED1,   "Sequential sensitivity corrections after state corrections")
 .export_values()
;

py::enum_<ODESLVS::Options::ASA_STRATEGY>(pyODESLVOptions, "ODESLV.ASA_STRATEGY")
 .value("HERMITE",    ODESLVS::Options::ASA_STRATEGY::HERMITE,    "Cubic Hermite interpolation")
 .value("POLYNOMIAL", ODESLVS::Options::ASA_STRATEGY::POLYNOMIAL, "Variable degree polynomial interpolation")
 .export_values()
;

py::enum_<ODESLVS::Options::NONLINEAR_SOLVER>(pyODESLVOptions, "ODESLV.NONLINEAR_SOLVER")
 .value("FIXEDPOINT", ODESLVS::Options::NONLINEAR_SOLVER::FIXEDPOINT, "Fixed point nonlinear solver")
 .value("NEWTON",     ODESLVS::Options::NONLINEAR_SOLVER::NEWTON,     "Newton nonlinear solver")
 .export_values()
;

py::enum_<ODESLVS::Options::LINEAR_SOLVER>(pyODESLVOptions, "ODESLV.LINEAR_SOLVER")
 .value("DIAG",    ODESLVS::Options::LINEAR_SOLVER::DIAG,    "Approximate diagonal Jacobian formed by way of a difference quotient")
 .value("DENSE",   ODESLVS::Options::LINEAR_SOLVER::DENSE,   "Use analytic dense Jacobian and internal direct dense linear algebra functions")
 .value("DENSEDQ", ODESLVS::Options::LINEAR_SOLVER::DENSEDQ, "Use approximate dense Jacobian by way of a difference quotient and internal direct dense linear algebra functions")
 .export_values()
;
}

