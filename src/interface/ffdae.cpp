// Copyright (C) Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the EPL-2.0 with GPL-2.0-or-later as a Secondary License; see the LICENSE file.
// SPDX-License-Identifier: EPL-2.0 OR GPL-2.0-or-later

// Python binding of FFDAESLV: the DAESLV_IDAS solver embedded as a DAG operation.
// FFDAESLV and FFDAESLV are instantiations of ONE template (FFDE<SOLVER>), so this binder is
// ffode.cpp with the solver type changed -- deliberately, so the two stay comparable line by line.
#include "ffde.hpp"

#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/typing.h>

#include <functional>

namespace py = pybind11;

namespace
{
typedef mc::FFModel M;
typedef mc::DAESLV_IDAS S;
typedef std::function<mc::FFVar(M::DofIndex const&)> t_Gen;
typedef py::typing::Dict<
    mc::FFVar,
    py::typing::Union<std::vector<mc::FFVar>,
                      py::typing::Callable<mc::FFVar(M::DofIndex const&)>>>
    t_Map;

// A dict {input: list of FFVar | gen(DofIndex) -> FFVar} -> the C++ InputArg
// list, in the dict's order (map 1's order is the order of the gradient
// columns)
std::vector<M::InputArg>
input_args(py::dict const& d)
{
  std::vector<M::InputArg> v;
  for (auto const& [k, val] : d)
  {
    mc::FFVar const u = k.cast<mc::FFVar>();
    if (PyCallable_Check(val.ptr()))
      v.emplace_back(u, val.cast<t_Gen>());
    else
      v.emplace_back(u, val.cast<std::vector<mc::FFVar>>());
  }
  return v;
}

// Python cannot hand over ownership of a solver it owns: TRANSFER is refused;
// under SHALLOW the solver is kept alive as long as its DAG's Python object
void
check_policy(int policy, S* solver, std::vector<M::InputArg> const& args,
             py::handle hsolver)
{
  if (policy == mc::FFDAESLV::TRANSFER)
    throw std::invalid_argument(
        "FFDAESLV: policy TRANSFER is not available from Python (the solver is "
        "owned by Python); use COPY (default) or SHALLOW");
  if (policy != mc::FFDAESLV::SHALLOW) return;
  mc::FFGraph* g = nullptr;
  for (auto const& a : args)
    if ((g = dynamic_cast<mc::FFGraph*>(a.input.dag()))) break;
  if (!g) g = solver->dag();
  py::handle hg = py::cast(g, py::return_value_policy::reference);
  if (hg && !hg.is_none()) py::detail::keep_alive_impl(hg, hsolver);
}

char const* DOC_MAPS = R"doc(
Parameters
----------
diff : dict
    Map 1: the inputs the numerical gradient is taken with respect to, in the
    order of the gradient columns. Each maps to its DOFs as DAG variables: a
    list in ``control_dofs`` order, or a generator ``gen(dof) -> FFVar``. It
    resets the solver's control registry to exactly these inputs.
rest : dict
    Map 2: every other declared input and every constant (one variable each),
    whose values pass through the DAG with zero numerical derivative.
solver : DAESLV
    The solver (set up). Under COPY (default) the operation holds its own
    copy; under SHALLOW it refers to ``solver``, which is then kept alive as
    long as the DAG.
policy : FFDAESLV.POLICY_TYPE, optional
    COPY (default) or SHALLOW.
name : str, optional
    Name of the operation in the DAG.
)doc";
}  // namespace

void
mc_ffdaeslv(py::module_& m)
{
  py::class_<mc::FFDAESLV, mc::FFOp> pyFFDAESLV(m, "FFDAESLV", R"doc(
A DAESLV solve embedded as an operation of the DAG: its outputs become DAG
variables, functions of the DAG variables given for the solver's inputs and
constants. Evaluating them integrates the model; differentiating them (FFGraph
``fdiff`` / ``bdiff``) uses forward or adjoint sensitivities
(``options.GRADIENT``), through FFGradDAESLV.

Example -- y(1) of the index-1 DAE  (2+p) y' + y = 0,  0 = z - p y,  y(0) = 1,
as a DAG variable of P:

>>> OpDAE = FFDAESLV()
>>> P = G.add_var("P")
>>> y = OpDAE({p: [P]}, solver)      # every input differentiated, no constant
>>> G.eval(y, [P], [0.7])             # [exp(-0.7)]
)doc");

  // --- FFDAESLV enums ---
  py::enum_<mc::FFDAESLV::POLICY_TYPE>(pyFFDAESLV, "POLICY_TYPE",
                                       "Copy policy of the embedded solver.")
      .value("SHALLOW", mc::FFDAESLV::SHALLOW,
             "Refer to the solver (kept alive as long as the DAG)")
      .value("COPY", mc::FFDAESLV::COPY,
             "Hold a deep copy of the solver (default)")
      .value("TRANSFER", mc::FFDAESLV::TRANSFER,
             "Not available from Python (ownership transfer)")
      .export_values();
  py::enum_<mc::FFDAESLV::FFBaseDE_t::GRADIENT_TYPE>(pyFFDAESLV, "GRADIENT_TYPE",
                                             "Sensitivity mode of derivatives.")
      .value("FORWARD", mc::FFDAESLV::FFBaseDE_t::FORWARD,
             "Forward sensitivities (one direction per parameter)")
      .value("ADJOINT", mc::FFDAESLV::FFBaseDE_t::ADJOINT,
             "Adjoint sensitivities (one backward sweep per output)")
      .value("AUTO", mc::FFDAESLV::FFBaseDE_t::AUTO,
             "ADJOINT when np > NP2NF * nf, else FORWARD (default)")
      .export_values();

  // --- FFDAESLV Options (static: shared by every FFDAESLV and FFGradDAESLV)
  // ---
  py::class_<mc::FFDAESLV::FFBaseDE_t::Options>(pyFFDAESLV, "Options", R"doc(
Options of FFDAESLV and FFGradDAESLV, shared by all operations (class-level:
``FFDAESLV.options``).
)doc")
      .def_readwrite("SYMDIFF", &mc::FFDAESLV::FFBaseDE_t::Options::SYMDIFF, R"doc(
Inputs or constants of the operation differentiated SYMBOLICALLY; the others
numerically (sensitivities).
)doc")
      .def_readwrite("NP2NF", &mc::FFDAESLV::FFBaseDE_t::Options::NP2NF, R"doc(
Parameter-to-output ratio above which AUTO uses adjoint sensitivities.
)doc")
      .def_readwrite("GRADIENT", &mc::FFDAESLV::FFBaseDE_t::Options::GRADIENT, R"doc(
Sensitivity mode: FFDAESLV.FORWARD, ADJOINT or AUTO.
)doc");
  pyFFDAESLV.def_property_static(
      "options", [](py::object) -> mc::FFDAESLV::FFBaseDE_t::Options&
      { return mc::FFDAESLV::FFBaseDE_t::options; },
      [](py::object, mc::FFDAESLV::FFBaseDE_t::Options const& o)
      { mc::FFDAESLV::FFBaseDE_t::options = o; }, py::return_value_policy::reference,
      "Options shared by every FFDAESLV and FFGradDAESLV operation.");

  // --- FFDAESLV Main Class ---
  pyFFDAESLV.def(py::init<>(), "Default constructor.")
      .def(
          "__call__",
          [](mc::FFDAESLV& self, t_Map const& diff, t_Map const& rest,
             S* solver, int policy, std::string const& name)
          {
            auto const a1 = input_args(diff), a2 = input_args(rest);
            auto a = a1;
            a.insert(a.end(), a2.begin(), a2.end());
            check_policy(policy, solver, a, py::cast(solver));
            return self(a1, a2, solver, policy, name);
          },
          py::arg("diff"), py::arg("rest"), py::arg("solver"),
          py::arg_v("policy", (int)mc::FFDAESLV::COPY, "FFDAESLV.COPY"),
          py::arg("name") = "",
          (std::string("TWO-MAP form: the outputs of the solve.\n") + DOC_MAPS)
              .c_str())
      .def(
          "__call__",
          [](mc::FFDAESLV& self, t_Map const& diff, S* solver, int policy,
             std::string const& name)
          {
            auto const a1 = input_args(diff);
            check_policy(policy, solver, a1, py::cast(solver));
            return self(a1, solver, policy, name);
          },
          py::arg("diff"), py::arg("solver"),
          py::arg_v("policy", (int)mc::FFDAESLV::COPY, "FFDAESLV.COPY"),
          py::arg("name") = "", R"doc(
ONE-MAP form: every declared input is differentiated and the model has no
constant (refused, naming what is missing, otherwise). Arguments as the
two-map form.
)doc")
      .def(
          "__call__",
          [](mc::FFDAESLV& self, unsigned idep, t_Map const& diff,
             t_Map const& rest, S* solver, int policy, std::string const& name)
          {
            auto const a1 = input_args(diff), a2 = input_args(rest);
            auto a = a1;
            a.insert(a.end(), a2.begin(), a2.end());
            check_policy(policy, solver, a, py::cast(solver));
            return self(idep, a1, a2, solver, policy, name);
          },
          py::arg("idep"), py::arg("diff"), py::arg("rest"), py::arg("solver"),
          py::arg_v("policy", (int)mc::FFDAESLV::COPY, "FFDAESLV.COPY"),
          py::arg("name") = "", "TWO-MAP form, output ``idep`` only.")
      .def(
          "__call__",
          [](mc::FFDAESLV& self, unsigned idep, t_Map const& diff, S* solver,
             int policy, std::string const& name)
          {
            auto const a1 = input_args(diff);
            check_policy(policy, solver, a1, py::cast(solver));
            return self(idep, a1, solver, policy, name);
          },
          py::arg("idep"), py::arg("diff"), py::arg("solver"),
          py::arg_v("policy", (int)mc::FFDAESLV::COPY, "FFDAESLV.COPY"),
          py::arg("name") = "", "ONE-MAP form, output ``idep`` only.");

  // --- FFGradDAESLV ---
  py::class_<mc::FFGradDAESLV, mc::FFOp>(m, "FFGradDAESLV", R"doc(
The derivative operation of FFDAESLV, created when an FFDAESLV output is
differentiated (FFGraph ``fdiff`` / ``bdiff``); not constructed directly.
)doc");
}
