// Copyright (C) Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

// Python binding of FFOCFESLV (an OCFESLV solve as a DAG operation) and
// FFOCFERES (the collocation residuals and outputs as a DAG operation), with
// their derivative operations.
#include "ffocfe.hpp"

#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/typing.h>

#include <functional>

namespace py = pybind11;

namespace
{
typedef mc::FFModel M;
typedef mc::OCFESLV O;
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
// under SHALLOW the solver is kept alive as long as the DAG's Python object
void
check_policy(char const* op, int policy, O* solver, mc::FFGraph* g,
             py::handle hsolver)
{
  if (policy == mc::FFOCFESLV::TRANSFER)
    throw std::invalid_argument(
        std::string(op) +
        ": policy TRANSFER is not available from Python (the solver is owned "
        "by Python); use COPY or SHALLOW");
  if (policy != mc::FFOCFESLV::SHALLOW) return;
  if (!g) g = solver->dag();
  py::handle hg = py::cast(g, py::return_value_policy::reference);
  if (hg && !hg.is_none()) py::detail::keep_alive_impl(hg, hsolver);
}

mc::FFGraph*
dag_of(std::vector<M::InputArg> const& args)
{
  for (auto const& a : args)
    if (auto* g = dynamic_cast<mc::FFGraph*>(a.input.dag())) return g;
  return nullptr;
}

mc::FFGraph*
dag_of(std::vector<mc::FFVar> const& v)
{
  for (auto const& x : v)
    if (auto* g = dynamic_cast<mc::FFGraph*>(x.dag())) return g;
  return nullptr;
}

char const* DOC_MAPS = R"doc(
Parameters
----------
diff : dict
    Map 1: the inputs the numerical gradient is taken with respect to (the
    controls), in the order of the gradient columns. Each maps to its DOFs as
    DAG variables: a list in ``control_dofs`` order, or a generator ``gen(dof)
    -> FFVar``. It resets the solver's control registry to exactly these
    inputs.
rest : dict
    Map 2: every other declared input and every constant (one variable each),
    whose values pass through the DAG with zero numerical derivative.
solver : OCFESLV
    The solver (set up). Under SHALLOW (default) the operation refers to
    ``solver``, which is then kept alive as long as the DAG; under COPY it
    holds its own copy.
policy : FFOCFESLV.POLICY_TYPE, optional
    SHALLOW (default) or COPY.
name : str, optional
    Name of the operation in the DAG.
)doc";
}  // namespace

void
mc_ffocfe(py::module_& m)
{
  // --- FFOCFESLV ---
  py::class_<mc::FFOCFESLV, mc::FFOp> pyFFOCFESLV(m, "FFOCFESLV", R"doc(
An OCFESLV solve embedded as an operation of the DAG: its outputs become DAG
variables, functions of the DAG variables given for the solver's inputs and
constants. Evaluating them solves the collocation system (monolithic or
marching); differentiating them (FFGraph ``fdiff`` / ``bdiff``) uses the
solver's forward or adjoint sensitivities (``options.GRADIENT``), through
FFGradOCFESLV.

Example -- the outputs of a solver whose inputs are its control ``q`` (one DOF
per time element) and a constant ``k``:

>>> OpOCFE = FFOCFESLV()
>>> Q = [G.add_var("Q%d" % i) for i in range(len(S.control_dofs(q)))]
>>> K = G.add_var("K")
>>> y = OpOCFE({q: Q}, {k: [K]}, S)
)doc");

  py::enum_<mc::FFOCFESLV::POLICY_TYPE>(pyFFOCFESLV, "POLICY_TYPE",
                                        "Copy policy of the embedded solver.")
      .value("SHALLOW", mc::FFOCFESLV::SHALLOW,
             "Refer to the solver (kept alive as long as the DAG; default)")
      .value("COPY", mc::FFOCFESLV::COPY, "Hold a deep copy of the solver")
      .value("TRANSFER", mc::FFOCFESLV::TRANSFER,
             "Not available from Python (ownership transfer)")
      .export_values();
  py::enum_<mc::FFBaseOCFE::GRADIENT_TYPE>(pyFFOCFESLV, "GRADIENT_TYPE",
                                           "Reduced-space gradient mode.")
      .value("FORWARD", mc::FFBaseOCFE::FORWARD,
             "One forward-sensitivity solve per control DOF")
      .value("ADJOINT", mc::FFBaseOCFE::ADJOINT,
             "One adjoint solve per output function")
      .value("AUTO", mc::FFBaseOCFE::AUTO,
             "ADJOINT when n_control_dof > NP2NF * n_colloc_fct, else FORWARD "
             "(default)")
      .export_values();

  // --- options (static: shared by every FFOCFESLV, FFGradOCFESLV, FFOCFERES)
  py::class_<mc::FFBaseOCFE::Options>(pyFFOCFESLV, "Options", R"doc(
Options of FFOCFESLV, FFGradOCFESLV and FFOCFERES, shared by all operations
(class-level: ``FFOCFESLV.options``).
)doc")
      .def_readwrite("SYMDIFF", &mc::FFBaseOCFE::Options::SYMDIFF, R"doc(
Inputs or constants of the operation differentiated SYMBOLICALLY; the others
numerically (sensitivities).
)doc")
      .def_readwrite("NP2NF", &mc::FFBaseOCFE::Options::NP2NF, R"doc(
Control-to-output ratio above which AUTO uses adjoint sensitivities.
)doc")
      .def_readwrite("GRADIENT", &mc::FFBaseOCFE::Options::GRADIENT, R"doc(
Reduced-space gradient mode: FFOCFESLV.FORWARD, ADJOINT or AUTO.
)doc");
  pyFFOCFESLV.def_property_static(
      "options", [](py::object) -> mc::FFBaseOCFE::Options&
      { return mc::FFBaseOCFE::options; },
      [](py::object, mc::FFBaseOCFE::Options const& o)
      { mc::FFBaseOCFE::options = o; }, py::return_value_policy::reference,
      "Options shared by every FFOCFESLV, FFGradOCFESLV and FFOCFERES "
      "operation.");

  pyFFOCFESLV.def(py::init<>(), "Default constructor.")
      .def(
          "__call__",
          [](mc::FFOCFESLV& self, t_Map const& diff, t_Map const& rest,
             O* solver, int policy, std::string const& name)
          {
            auto const a1 = input_args(diff), a2 = input_args(rest);
            mc::FFGraph* g = dag_of(a1);
            if (!g) g = dag_of(a2);
            check_policy("FFOCFESLV", policy, solver, g, py::cast(solver));
            return self(a1, a2, solver, policy, name);
          },
          py::arg("diff"), py::arg("rest"), py::arg("solver"),
          py::arg_v("policy", (int)mc::FFOCFESLV::SHALLOW, "FFOCFESLV.SHALLOW"),
          py::arg("name") = "",
          (std::string("TWO-MAP form: the outputs of the solve.\n") + DOC_MAPS)
              .c_str())
      .def(
          "__call__",
          [](mc::FFOCFESLV& self, t_Map const& diff, O* solver, int policy,
             std::string const& name)
          {
            auto const a1 = input_args(diff);
            check_policy("FFOCFESLV", policy, solver, dag_of(a1),
                         py::cast(solver));
            return self(a1, solver, policy, name);
          },
          py::arg("diff"), py::arg("solver"),
          py::arg_v("policy", (int)mc::FFOCFESLV::SHALLOW, "FFOCFESLV.SHALLOW"),
          py::arg("name") = "", R"doc(
ONE-MAP form: every declared input is differentiated and the model has no
constant (refused, naming what is missing, otherwise). Arguments as the
two-map form.
)doc");

  py::class_<mc::FFGradOCFESLV, mc::FFOp>(m, "FFGradOCFESLV", R"doc(
The derivative operation of FFOCFESLV, created when an FFOCFESLV output is
differentiated (FFGraph ``fdiff`` / ``bdiff``); not constructed directly.
)doc");

  // --- FFOCFERES ---
  py::class_<mc::FFOCFERES, mc::FFOp> pyFFOCFERES(m, "FFOCFERES", R"doc(
The collocation system of an OCFESLV as an operation of the DAG: given DAG
variables for the collocation state coefficients (``n_colloc_sta()``), the
collocation input coefficients (``n_colloc_inp()``) and the constants, its
values are the collocation residuals (``n_colloc_eqn()``) followed by the
outputs (``n_colloc_fct()``) -- the constraints and objective terms of a
full-discretisation formulation. Its derivatives are the sparse collocation
Jacobian, through FFGradOCFERES.
)doc");

  py::enum_<mc::FFOCFERES::POLICY_TYPE>(pyFFOCFERES, "POLICY_TYPE",
                                        "Copy policy of the embedded solver.")
      .value("SHALLOW", mc::FFOCFERES::SHALLOW,
             "Refer to the solver (kept alive as long as the DAG)")
      .value("COPY", mc::FFOCFERES::COPY,
             "Hold a deep copy of the solver (default)")
      .value("TRANSFER", mc::FFOCFERES::TRANSFER,
             "Not available from Python (ownership transfer)")
      .export_values();

  auto res_call = [](mc::FFOCFERES& self, std::vector<mc::FFVar> const& sta,
                     std::vector<mc::FFVar> const& inp,
                     std::vector<mc::FFVar> const& cst, O* solver, int policy,
                     std::string const& name)
  {
    if (sta.size() != solver->n_colloc_sta())
      throw std::invalid_argument(
          "FFOCFERES: " + std::to_string(solver->n_colloc_sta()) +
          " collocation state variables expected, " +
          std::to_string(sta.size()) + " given");
    if (inp.size() != solver->n_colloc_inp())
      throw std::invalid_argument(
          "FFOCFERES: " + std::to_string(solver->n_colloc_inp()) +
          " collocation input variables expected, " +
          std::to_string(inp.size()) + " given");
    if (cst.size() != solver->var_constant().size())
      throw std::invalid_argument(
          "FFOCFERES: " + std::to_string(solver->var_constant().size()) +
          " constant variables expected, " + std::to_string(cst.size()) +
          " given");
    mc::FFGraph* g = dag_of(sta);
    if (!g) g = dag_of(inp);
    check_policy("FFOCFERES", policy, solver, g, py::cast(solver));
    return self(sta, inp, cst, solver, policy, name);
  };
  pyFFOCFERES.def(py::init<>(), "Default constructor.")
      .def("__call__", res_call, py::arg("sta"), py::arg("inp"), py::arg("cst"),
           py::arg("solver"),
           py::arg_v("policy", (int)mc::FFOCFERES::COPY, "FFOCFERES.COPY"),
           py::arg("name") = "", R"doc(
The collocation residuals then the outputs, as DAG variables.

Parameters
----------
sta : list of FFVar
    DAG variables for the collocation state coefficients (``n_colloc_sta()``).
inp : list of FFVar
    DAG variables for the collocation input coefficients (``n_colloc_inp()``).
cst : list of FFVar
    DAG variables for the constants, in the order of ``var_constant``.
solver : OCFESLV
    The solver (set up, monolithic).
policy : FFOCFERES.POLICY_TYPE, optional
    COPY (default) or SHALLOW.
name : str, optional
    Name of the operation in the DAG.
)doc")
      .def(
          "__call__",
          [res_call](mc::FFOCFERES& self, std::vector<mc::FFVar> const& sta,
                     std::vector<mc::FFVar> const& inp, O* solver, int policy,
                     std::string const& name)
          {
            return res_call(self, sta, inp, std::vector<mc::FFVar>(), solver,
                            policy, name);
          },
          py::arg("sta"), py::arg("inp"), py::arg("solver"),
          py::arg_v("policy", (int)mc::FFOCFERES::COPY, "FFOCFERES.COPY"),
          py::arg("name") = "", "As above, for a model without constants.");

  py::class_<mc::FFGradOCFERES, mc::FFOp>(m, "FFGradOCFERES", R"doc(
The derivative operation of FFOCFERES, created when an FFOCFERES output is
differentiated; not constructed directly.
)doc");
}
