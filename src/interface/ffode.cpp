#include <pybind11/pybind11.h>
//#include <pybind11/operators.h>
#include <pybind11/stl.h>

#include <fstream>
#include "ffode.hpp" 

namespace py = pybind11;

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
 .def_readwrite( "SYMDIFF", &mc::FFODE::Options::SYMDIFF,   "method of ODE differentiation [Default: -]" )
 .def_readwrite( "NP2NF", &mc::FFODE::Options::NP2NF, "parameter-to-function-size ratio above which adjoint sensitivity is applied instead of forward sensitivity [Default: 3]" )
;

//py::enum_<mc::FFODE::Options::DERIV_TYPE>( pyFFODEOptions, "DERIV_TYPE" )
// .value("NUM_P", mc::FFODE::Options::DERIV_TYPE::NUM_P, "Derivatives w.r.t. parameters only through forward or adjoint sensitivity integration")
// .value("SYM_P", mc::FFODE::Options::DERIV_TYPE::SYM_P, "Derivatives w.r.t. parameters only through symbolic differentiation of ODEs")
// .value("SYM_C", mc::FFODE::Options::DERIV_TYPE::SYM_C, "Derivatives w.r.t. constants only through symbolic differentiation of ODEs")
// .value("SYM_PC", mc::FFODE::Options::DERIV_TYPE::SYM_PC, "Derivatives w.r.t. parameters and constants jointly through symbolic differentiation of ODEs")
// .export_values()
//;
}

