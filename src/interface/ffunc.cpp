#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>

#include <fstream>
#include "ffunc.hpp" 
#include "odeslvs_cvodes.hpp" 

namespace py = pybind11;

void mc_ffunc( py::module_ &m )
{

typedef mc::FFGraph< mc::FFODE<0>, mc::FFGRADODE<0> > FFGraphExt;
typedef mc::ODESLVS_CVODES<> ODESLVS;

py::class_<mc::FFODE<0>, mc::FFOp> pyFFODE( m, "FFODE" );//, py::module_local() );

pyFFODE
 .def(
   py::init<>(),
   "default constructor"
 )
 .def(
   "__call__",
   []( mc::FFODE<0>& self, std::vector<mc::FFVar> const& vVar, ODESLVS* pODE )
   {
     auto pDep = self( vVar.size(), vVar.data(), pODE );
     return std::vector<mc::FFVar*>( pDep, pDep+pODE->nf() );
   },
   py::return_value_policy::reference_internal,
   "define ODE operation in DAG"
 )
 .def_property_readonly(
   "options",
   []( mc::FFODE<0>& self )
   {
     ODESLVS* data = static_cast<ODESLVS*>( self.data );
     return data->options;
   },
   py::return_value_policy::reference_internal,
   "ODE options"
 )
 .def(
   "setup",
   []( mc::FFODE<0>& self )
   {
     ODESLVS* data = static_cast<ODESLVS*>( self.data );
     data->setup();
   },
   "ODE setup"
 )
 .def_property_readonly(
   "data",
   []( mc::FFODE<0>& self ) -> mc::ODESLVS_CVODES<>*
   {
     ODESLVS* data = static_cast<ODESLVS*>( self.data );
     return data;
   },
   "ODE object"
 )
 .def_readwrite(
   "type",
   &mc::FFODE<0>::type,
   "retreive operation type"
 )
 .def_readwrite(
   "varin",
   &mc::FFODE<0>::varin
 )
 .def_readwrite(
   "varout",
   &mc::FFODE<0>::varout
 )
 .def( "name",
   &mc::FFODE<0>::name,
   "retreive operation name"
 )
 .def(
   "__str__",
   []( mc::FFODE<0> const& O )
   {
     std::ostringstream Oss;
     Oss << O;
     return Oss.str();
   }
 )
 .def(
   "__repr__",
   []( mc::FFODE<0> const& O )
   {
     std::ostringstream Oss;
     Oss << O;
     return Oss.str();
   }
 )
;

py::class_<FFGraphExt, mc::FFBase> pyFFGraph( m, "FFGraphExt" );//, py::module_local() );

pyFFGraph
 .def(
   py::init<>(),
   "default constructor"
 )
 .def(
   "fdiff",
   []( FFGraphExt& G, std::vector<mc::FFVar const*> const& vDep, std::vector<mc::FFVar const*> const& vIndep )
   {
     return G.SFAD( vDep, vIndep );
   },
   py::return_value_policy::reference_internal,
   "apply forward differentiation"
 )
 .def(
   "fdiff",
   []( FFGraphExt& G, std::vector<mc::FFVar const*> const& vDep, std::vector<mc::FFVar const*> const& vIndep, std::vector<mc::FFVar const*> const& vDir )
   {
     return G.SFAD( vDep, vIndep, vDir );
   },
   py::return_value_policy::reference_internal,
   "apply directional forward differentiation"
 )
 .def(
   "bdiff",
   []( FFGraphExt& G, std::vector<mc::FFVar const*> const& vDep, std::vector<mc::FFVar const*> const& vIndep )
   {
     return G.SBAD( vDep, std::vector<mc::FFVar const*>(), vIndep );
   },
   py::return_value_policy::reference_internal,
   "apply backward differentiation"
 )
 .def(
   "bdiff",
   []( FFGraphExt& G, std::vector<mc::FFVar const*> const& vDep, std::vector<mc::FFVar const*> const& vDir,
       std::vector<mc::FFVar const*> const& vIndep )
   { 
     return G.SBAD( vDep, vDir, vIndep );
   },
   py::return_value_policy::reference_internal,
   "apply directional backward differentiation"
 )
 .def(
   "tdiff",
   []( FFGraphExt& G, unsigned int const ordermax, std::vector<mc::FFVar const*> const& vDep,
       std::vector<mc::FFVar const*> const& vVar, mc::FFVar const* const pIndep )
   {
     return G.TAD( ordermax, vDep, vVar, pIndep );
   },
   py::return_value_policy::reference_internal,
   "apply Taylor expansion"
 )
 .def(
   "compose",
   []( FFGraphExt& G, std::vector<mc::FFVar const*> const& vDepOut, std::vector< std::pair<mc::FFVar const*,
       mc::FFVar const*> > const& vDepIn )
   {
     return G.compose( vDepOut, vDepIn );
   },
   "apply compostion"
 )
 .def(
   "eval",
   []( FFGraphExt& G, mc::FFSubgraph& SgDep, std::vector<mc::FFVar> const& vDep, std::vector<mc::FFVar> const& vVar,
       std::vector<double> const& DVar )
   {
     size_t const nDep = vDep.size();
     std::vector<double> DDep( nDep );
     G.eval( SgDep, nDep, vDep.data(), DDep.data(), vVar.size(), vVar.data(), DVar.data() );
     return DDep;
   },
   py::return_value_policy::take_ownership,
   "evaluate subgraph in double arithmetic"
 )
 .def(
   "eval",
   []( FFGraphExt& G, std::vector<mc::FFVar> const& vDep, std::vector<mc::FFVar> const& vVar, std::vector<double> const& DVar )
   {
     size_t const nDep = vDep.size();
     std::vector<double> DDep( nDep );
     G.eval( nDep, vDep.data(), DDep.data(), vVar.size(), vVar.data(), DVar.data() );
     return DDep;
   },
   py::return_value_policy::take_ownership,
   "evaluate subgraph in double arithmetic"
 )
 .def(
   "__str__",
   []( FFGraphExt const& G )
   {
     std::ostringstream Gss;
     Gss << G;
     return Gss.str();
   }
 )
 .def(
   "__repr__",
   []( FFGraphExt const& G )
   {
     std::ostringstream Gss;
     Gss << G;
     return Gss.str();
   }
 )
;
}

