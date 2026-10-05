// Copyright (C) Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

// Python binding of FFModel: the model of a differential-algebraic system on
// the MC++ DAG (domains, states, inputs, constants, equations, outputs,
// transitions, controls), shared by every CRONOS solver.
//
// Python-specific forms, where C++ overloads cannot be told apart safely from
// Python values:
//  * add_output( fct, doms, point=..., side=... | masks=..., at=... ): a POINT
//  output takes coordinates, a
//    DISTRIBUTED one takes masks -- [1.] and [1] must not select different
//    overloads silently;
//  * add_input( var, doms, ref=..., type=..., n_node=..., continuity=...,
//  is_decision=... ): one entry for the
//    seven C++ overloads, dispatched on the keywords given (a combination C++
//    has no overload for is refused).
// Reference values may be a float or a callable taking {FFVar: coordinate} and
// returning a float.
#include "ffmodel.hpp"

#include <pybind11/functional.h>
#include <pybind11/iostream.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <optional>
#include <sstream>
#include <variant>

#include "ffexpr.hpp"
#include "gen_ffmodel_options.hpp"

namespace py = pybind11;

namespace
{
typedef mc::FFModel M;
typedef M::t_Fun t_Fun;
typedef std::variant<double, t_Fun> t_Ref;

// A region of domain n, written explicitly with symbolic bounds (a record does
// not know them; the model report writes the numbers): "t in (lo, up]", "t =
// lo"
std::string
region_str(mc::FFVar const& v, int mask)
{
  std::string const n = v.name();
  switch (mask)
  {
    case mc::FFDom::ALL:
      return n + " in [lo, up]";
    case mc::FFDom::ALL - mc::FFDom::LB:
      return n + " in (lo, up]";
    case mc::FFDom::ALL - mc::FFDom::UB:
      return n + " in [lo, up)";
    case mc::FFDom::ALL - mc::FFDom::LB - mc::FFDom::UB:
      return n + " in (lo, up)";
    case mc::FFDom::LB:
      return n + " = lo";
    case mc::FFDom::UB:
      return n + " = up";
    default:
      return n + " (mask " + std::to_string(mask) + ")";
  }
}

// An expression of the DAG written out by FFExpr (its name if it cannot be)
std::string
expr_str(mc::FFVar const& v)
{
  auto* g = dynamic_cast<mc::FFGraph*>(v.dag());
  if (!g) return v.name();
  try
  {
    mc::FFSubgraph sg                = g->subgraph(1, &v);
    std::vector<mc::FFExpr> const ex = mc::FFExpr::subgraph(g, sg);
    if (ex.empty()) return v.name();
    std::ostringstream o;
    o << ex[0];
    return o.str();
  }
  catch (...)
  {
    return v.name();
  }
}

std::string
eqn_str(M::t_Eqn const& e)
{
  std::string d;
  for (auto const& [v, mask] : e.dom)
    d += (d.empty() ? "" : ", ") + region_str(v, mask);
  return "0 = " + expr_str(e.var) + (d.empty() ? std::string() : "   on " + d);
}

std::string
fct_str(M::t_Fct const& f)
{
  std::string w, g;
  for (auto const& [v, c] : f.point)
  {
    auto const is = f.side.find(v);
    std::ostringstream o;
    o << c;
    w += (w.empty() ? "" : ", ") + v.name() + " = " + o.str() +
         (is != f.side.end() && is->second == mc::FFDom::PLUS ? "+" : "");
  }
  for (auto const& [v, mask] : f.grid)
    g += (g.empty() ? "" : ", ") + region_str(v, mask);
  std::string const x = expr_str(f.var);
  if (f.kind == M::FctKind::DISTRIBUTED)
    return x + "   distributed on " + (g.empty() ? std::string("-") : g) +
           (w.empty() ? "" : "   at " + w);
  return x + (w.empty() ? "" : "   at " + w);
}

void
add_input(M& self, mc::FFVar const& var, std::vector<mc::FFVar> const& doms,
          std::optional<t_Ref> const& ref,
          std::optional<mc::FFDom::TYPE> const& type, size_t n_node,
          std::optional<std::vector<M::InputContinuity>> const& cnt,
          bool is_decision)
{
  bool const fun = ref && std::holds_alternative<t_Fun>(*ref);
  std::optional<double> const val =
      (ref && !fun) ? std::optional<double>(std::get<double>(*ref))
                    : std::nullopt;
  if (cnt)
  {
    if (type || fun || is_decision)
      throw std::invalid_argument(
          "add_input: 'continuity' combines with a float 'ref' only");
    return self.add_input(var, doms, *cnt, val);
  }
  if (type)
  {
    if (is_decision)
      throw std::invalid_argument(
          "add_input: 'type' does not combine with 'is_decision'");
    if (fun)
      return self.add_input(var, doms, *type, n_node, std::get<t_Fun>(*ref));
    return self.add_input(var, doms, *type, n_node, val);
  }
  if (n_node)
    throw std::invalid_argument("add_input: 'n_node' requires 'type'");
  if (fun) return self.add_input(var, doms, std::get<t_Fun>(*ref), is_decision);
  if (is_decision)
  {
    if (!doms.empty() || !val)
      throw std::invalid_argument(
          "add_input: a decision input without domains needs a float 'ref'");
    return self.add_input(var, *val, true);
  }
  return self.add_input(var, doms, val);
}

void
add_output(M& self, mc::FFVar const& fct, std::vector<mc::FFVar> const& doms,
           std::optional<std::vector<double>> const& point,
           std::optional<std::vector<int>> const& side,
           std::optional<std::vector<int>> const& masks,
           std::optional<std::map<mc::FFVar, double, mc::lt_FFVar>> const& at)
{
  if (point && masks)
    throw std::invalid_argument(
        "add_output: give 'point' (a point output) or 'masks' (a distributed "
        "one), not both");
  if (masks)
  {
    if (side)
      throw std::invalid_argument(
          "add_output: 'side' applies to a point output");
    if (!at) return self.add_output(fct, doms, *masks);
    std::vector<mc::FFVar> pd;
    std::vector<double> pv;
    for (auto const& [d, v] : *at)
    {
      pd.push_back(d);
      pv.push_back(v);
    }
    return self.add_output(fct, doms, *masks, pd, pv);
  }
  if (at)
    throw std::invalid_argument(
        "add_output: 'at' applies to a distributed output (give 'masks')");
  if (side)
    return self.add_output(fct, doms, point ? *point : std::vector<double>(),
                           *side);
  return self.add_output(fct, doms, point ? *point : std::vector<double>());
}
}  // namespace

void
mc_ffmodel(py::module_& m)
{
  py::class_<M> pyM(m, "FFModel", R"doc(
Model of a differential-algebraic system on the DAG of MC++, shared by the
CRONOS solvers.

A model is declared on an ``FFGraph``: domains of the independent variables
(``add_domain``, with an ``FFDom``), states and inputs over subsets of these
domains, constants, equations with a region (masks per domain) and a role,
outputs (point values, integrals, or distributed values), and transitions
(jumps of the states at an interior point of the evolution domain).
``setup()`` validates and analyses the model; the solvers (``ODESLV``,
``OCFESLV``) derive from ``FFModel`` and call it themselves.

Example -- the scalar ODE x' = -k x, x(0) = 1, output x(1):

>>> from pymcpp import FFGraph, FFPartial, FFEval
>>> from cronos import FFDom, FFModel, ODESLV
>>> G = FFGraph(); t, x = G.add_var("t"), G.add_var("x")
>>> OpP, OpE = FFPartial(), FFEval()
>>> M = ODESLV(G)
>>> M.add_domain(t, FFDom(0., 1., 1)); M.set_evolution_domain(t)
>>> M.add_state(x, [t])
>>> M.add_equation(OpP(x, t) + .7*x, [t], [FFDom.ALL - FFDom.LB])
>>> M.add_equation(x - 1., [t], [FFDom.LB])
>>> M.add_output(OpE(x, t, 1.))
)doc");

  // --- FFModel DofIndex ---
  py::class_<M::DofIndex>(pyM, "DofIndex", R"doc(
Position of one degree of freedom of a distributed input: its element and its
node within the element, per domain.
)doc")
      .def(py::init<>())
      .def_readwrite("element", &M::DofIndex::element,
                     "Element index per domain, ``{domain: element}``.")
      .def_readwrite(
          "node", &M::DofIndex::node,
          "Node index within the element per domain, ``{domain: node}``.");

  // --- FFModel Exceptions ---
  py::class_<M::Exceptions> pyFFModelExceptions(pyM, "Exceptions", R"doc(
Exception information for errors raised by a model or its solvers (the same
type as OCFESLV.Exceptions).

Instances carry an error code (ierr) and a human-readable description (what).
)doc");
  py::enum_<M::Exceptions::TYPE>(pyFFModelExceptions, "TYPE",
                                 "Error codes for FFModel exceptions")
      .value("SETUP", M::Exceptions::SETUP,
             "Incomplete setup before evaluation")
      .value("INDEX", M::Exceptions::INDEX, "Index mismatch in data access")
      .value("DOMAIN", M::Exceptions::DOMAIN,
             "Variable-domain misspecification in expression")
      .value("ENV", M::Exceptions::ENV,
             "Environment mismatch in collocation operation")
      .value("CSTVAL", M::Exceptions::CSTVAL, "Undefined constant values")
      .value("UNDEF", M::Exceptions::UNDEF, "Undefined collocation operation")
      .value("INTERNAL", M::Exceptions::INTERNAL, "Internal error")
      .value("NOSTORE", M::Exceptions::NOSTORE,
             "No stored solution available for a buffer-free read")
      .export_values();
  pyFFModelExceptions
      .def("ierr", &M::Exceptions::ierr, R"doc(
Return the error code as an integer (see FFModel.Exceptions.TYPE).
)doc")
      .def(
          "what", [](M::Exceptions const& e) { return std::string(e.what()); },
          R"doc(
Return a string describing the error.
)doc");

  // --- enums ---
  py::enum_<M::EqnRole>(pyM, "EqnRole", "Role of an equation.")
      .value("AUTO", M::EqnRole::AUTO,
             "resolved on add: BOUNDARY if some domain mask is LB or UB alone, "
             "else INTERIOR")
      .value("INTERIOR", M::EqnRole::INTERIOR,
             "volume/interior PDE or algebraic residual")
      .value("INITIAL", M::EqnRole::INITIAL, "initial-condition trace residual")
      .value("BOUNDARY", M::EqnRole::BOUNDARY,
             "exterior boundary-condition trace residual")
      .value("INTERFACE", M::EqnRole::INTERFACE,
             "user-supplied transmission residual between physical subdomains")
      .value("LINK", M::EqnRole::LINK,
             "auxiliary first-order link equation generated by order reduction")
      .value("SURFACE", M::EqnRole::SURFACE,
             "residual posed on a lower-dimensional boundary or interface")
      .value("DIAGNOSTIC", M::EqnRole::DIAGNOSTIC,
             "evaluated only: not classified, receives no SAT, never donated");
  py::enum_<M::InputContinuity>(
      pyM, "InputContinuity",
      "Continuity of a distributed input across element interfaces.")
      .value("DISCONTINUOUS", M::DISCONTINUOUS,
             "default: the input may jump across element interfaces")
      .value("CONTINUOUS_C0", M::CONTINUOUS_C0,
             "the input value matches across interfaces (C0)")
      .value("SMOOTH_C1", M::SMOOTH_C1,
             "value and first derivative match across interfaces (C1)")
      .export_values();
  py::enum_<M::SetupStatus>(pyM, "SetupStatus", "Outcome of setup().")
      .value("OK", M::SetupStatus::OK, "setup succeeded")
      .value("INCONSISTENT_MODEL", M::SetupStatus::INCONSISTENT_MODEL,
             "the declared model is inconsistent")
      .value("NODE_SETUP_FAILED", M::SetupStatus::NODE_SETUP_FAILED,
             "collocation node setup failed")
      .value("CLASSIFICATION_FAILED", M::SetupStatus::CLASSIFICATION_FAILED,
             "PDE classification failed")
      .value("HYP_INCOMING_BC", M::SetupStatus::HYP_INCOMING_BC,
             "hyperbolic block: incoming boundary condition missing")
      .value(
          "HYP_BC_MISDIRECTED", M::SetupStatus::HYP_BC_MISDIRECTED,
          "hyperbolic block: boundary condition on an outgoing characteristic")
      .value("INTERFACE_PLAN_INVALID", M::SetupStatus::INTERFACE_PLAN_INVALID,
             "the interface plan is invalid")
      .value("AUDIT_NONSQUARE", M::SetupStatus::AUDIT_NONSQUARE,
             "the discretised system is not square")
      .value("DERIV_CACHE_FAILED", M::SetupStatus::DERIV_CACHE_FAILED,
             "derivative cache construction failed")
      .value("LINEAR_CACHE_FAILED", M::SetupStatus::LINEAR_CACHE_FAILED,
             "linear-evaluation cache construction failed");

  // --- options ---
  py::class_<M::EqnOptions>(
      pyM, "EqnOptions",
      "Options of one equation: its role and classification block.")
      .def(py::init<M::EqnRole, int>(),
           py::arg_v("role", M::EqnRole::AUTO, "FFModel.EqnRole.AUTO"),
           py::arg("block_id") = 0,
           "Options with role-derived defaults (INITIAL, BOUNDARY, INTERFACE "
           "and DIAGNOSTIC rows do not enter classification).")
      .def_readwrite("role", &M::EqnOptions::role,
                     "equation role (AUTO resolved on add)")
      .def_readwrite("participate_in_classification",
                     &M::EqnOptions::participate_in_classification,
                     "enters block principal-symbol classification")
      .def_readonly("role_auto", &M::EqnOptions::role_auto,
                    "True if the role was resolved from AUTO (re-resolved once "
                    "the evolution domain is known).")
      .def_readwrite("block_id", &M::EqnOptions::block_id,
                     "classification block");

  py::class_<M::Options> pyOpt(
      pyM, "Options",
      "Model options: display, order reduction, classification, tolerances.");
  pyOpt.def(py::init<>()).def(py::init<M::Options const&>());
  bind_ffmodel_options(pyOpt);

  // --- FFModel records ---
  py::enum_<M::FctKind>(pyM, "FctKind", "Kind of an output.")
      .value("POINT", M::FctKind::POINT,
             "a scalar output, or a value at a point of its domains")
      .value("DISTRIBUTED", M::FctKind::DISTRIBUTED,
             "an output distributed over a region of its domains")
      .export_values();

  py::class_<M::t_Eqn>(pyM, "t_Eqn", R"doc(
A declared equation: its residual, its region (mask per domain), and its
options.
)doc")
      .def_readonly("var", &M::t_Eqn::var,
                    "Residual expression (equation: var = 0).")
      .def_readonly("dom", &M::t_Eqn::dom,
                    "Domain mask per domain variable, ``{domain: mask}``.")
      .def_property_readonly(
          "opt", [](M::t_Eqn const& e) { return M::EqnOptions(*e.opt); },
          "Normalised equation options (a copy).")
      .def("__str__", &eqn_str,
           "The equation written out: ``0 = <expression>   on <region>``.")
      .def("__repr__",
           [](M::t_Eqn const& e) { return "t_Eqn(" + eqn_str(e) + ")"; });

  py::class_<M::t_Fct>(pyM, "t_Fct", R"doc(
A declared output: its expression, kind, fixed coordinates and sides, and the
masks of its distributed dimensions.
)doc")
      .def_readonly("var", &M::t_Fct::var, "Output expression.")
      .def_readonly("kind", &M::t_Fct::kind, "Point or distributed output.")
      .def_readonly("point", &M::t_Fct::point,
                    "Fixed coordinates, ``{domain: coordinate}``.")
      .def_readonly("side", &M::t_Fct::side,
                    "Sides given as PLUS, ``{domain: 1}`` (absent: MINUS).")
      .def_readonly("grid", &M::t_Fct::grid,
                    "Masks of the distributed dimensions, ``{domain: mask}``.")
      .def_readonly("row0", &M::t_Fct::row0,
                    "First row in the output array (set by the solver).")
      .def_readonly("nrow", &M::t_Fct::nrow,
                    "Number of output rows (set by the solver).")
      .def("__str__", &fct_str,
           "The output written out, with its point or its region.")
      .def("__repr__",
           [](M::t_Fct const& f) { return "t_Fct(" + fct_str(f) + ")"; });

  // --- model ---
  pyM.def_static("revision", []() { return std::string(M::revision()); },
                 "The revision of the model layer (ffmodel.hpp), e.g. for a bug report.  The model report\n"
                 "does not print it.");
  pyM.def(py::init<mc::FFGraph*>(), py::arg("dag"), py::keep_alive<1, 2>(),
          "Empty model declared on the DAG ``dag`` (kept alive by the model).")
      .def_readwrite("options", &M::options,
                     "Model options (the solvers extend them).")
      .def(
          "set_constant",
          [](M& self, std::vector<mc::FFVar> const& c,
             std::vector<double> const& v) { self.set_constant(c, v); },
          py::arg("vars"), py::arg("values") = std::vector<double>(), R"doc(
Declare constants: parameters whose values are given at solve time rather than
differentiated.

Parameters
----------
vars : list of FFVar
    The constants.
values : list of float, optional
    Their values, one per constant.
)doc")
      .def(
          "set_evolution_domain", [](M& self, mc::FFVar const& t)
          { self.set_evolution_domain(t); }, py::arg("dom"),
          "Declare which domain variable is the evolution direction (time); it "
          "must have been added with ``add_domain``.")
      .def(
          "add_domain",
          [](M& self, mc::FFVar const& v, mc::FFDom const& d,
             std::optional<double> ref) { self.add_domain(v, d, ref); },
          py::arg("var"), py::arg("dom"), py::arg("ref") = std::nullopt, R"doc(
Declare the domain of an independent variable.

Parameters
----------
var : FFVar
    The independent variable.
dom : FFDom
    Its bounds, finite elements and collocation nodes.
ref : float, optional
    Reference coordinate for the model analyses.
)doc")
      .def(
          "add_domain",
          [](M& self, std::vector<mc::FFVar> const& v, mc::FFDom const& d,
             std::vector<double> const& ref) { self.add_domain(v, d, ref); },
          py::arg("vars"), py::arg("dom"),
          py::arg("ref") = std::vector<double>(),
          "Declare several independent variables on the same domain.")
      .def(
          "add_state",
          [](M& self, mc::FFVar const& v, std::vector<mc::FFVar> const& d,
             std::optional<double> ref) { self.add_state(v, d, ref); },
          py::arg("var"), py::arg("doms"), py::arg("ref") = std::nullopt, R"doc(
Declare a state over the domains ``doms`` (empty for a scalar state).

Parameters
----------
var : FFVar
    The state.
doms : list of FFVar
    Its independent variables.
ref : float, optional
    Reference value (initial guess and analyses).
)doc")
      .def(
          "add_state",
          [](M& self, mc::FFVar const& v, std::vector<mc::FFVar> const& d,
             t_Fun const& ref) { self.add_state(v, d, ref); },
          py::arg("var"), py::arg("doms"), py::arg("ref"),
          "Declare a state with a reference function ``ref({dom: coordinate}) "
          "-> float``.")
      .def(
          "add_state",
          [](M& self, std::vector<mc::FFVar> const& v,
             std::vector<mc::FFVar> const& d, std::vector<double> const& ref)
          { self.add_state(v, d, ref); },
          py::arg("vars"), py::arg("doms"),
          py::arg("ref") = std::vector<double>(),
          "Declare several states over the same domains.")
      .def("add_input", &add_input, py::arg("var"),
           py::arg("doms") = std::vector<mc::FFVar>(), py::kw_only(),
           py::arg("ref") = std::nullopt, py::arg("type") = std::nullopt,
           py::arg("n_node") = 0, py::arg("continuity") = std::nullopt,
           py::arg("is_decision") = false, R"doc(
Declare an input: a parameter (no domains) or a distributed input over
``doms``.

Parameters
----------
var : FFVar
    The input.
doms : list of FFVar, optional
    Its independent variables (empty: a scalar parameter).
ref : float or callable, optional
    Reference value, or a function ``ref({dom: coordinate}) -> float``.
type : FFDom.TYPE, optional
    The input's own collocation node family (with ``n_node``); default: that
    of its domains.
n_node : int, optional
    Number of nodes per element for ``type``.
continuity : list of FFModel.InputContinuity, optional
    Continuity per domain across element interfaces (with a float ``ref``
    only).
is_decision : bool, optional
    Register the input as a control (decision) at setup.
)doc")
      .def(
          "add_input",
          [](M& self, std::vector<mc::FFVar> const& v,
             std::vector<mc::FFVar> const& d, std::vector<double> const& ref)
          { self.add_input(v, d, ref); },
          py::arg("vars"), py::arg("doms") = std::vector<mc::FFVar>(),
          py::arg("ref") = std::vector<double>(),
          "Declare several inputs over the same domains.")
      .def(
          "add_equation",
          [](M& self, mc::FFVar const& e, std::vector<mc::FFVar> const& d,
             std::vector<int> const& l, M::EqnOptions const& o)
          { self.add_equation(e, d, l, o); },
          py::arg("eqn"), py::arg("doms") = std::vector<mc::FFVar>(),
          py::arg("masks") = std::vector<int>(),
          py::arg_v("options", M::EqnOptions(), "FFModel.EqnOptions()"), R"doc(
Declare the equation ``eqn = 0`` on a region of its domains.

Parameters
----------
eqn : FFVar
    The residual.
doms : list of FFVar
    Its domains (empty for a scalar equation).
masks : list of int
    The region per domain: ``FFDom.ALL``, ``LB``, ``UB``, or combinations such
    as ``FFDom.ALL - FFDom.LB``.
options : FFModel.EqnOptions, optional
    Role and classification block (default: role resolved from the masks).
)doc")
      .def(
          "add_equation",
          [](M& self, std::vector<mc::FFVar> const& e,
             std::vector<mc::FFVar> const& d, std::vector<int> const& l,
             M::EqnOptions const& o) { self.add_equation(e, d, l, o); },
          py::arg("eqns"), py::arg("doms") = std::vector<mc::FFVar>(),
          py::arg("masks") = std::vector<int>(),
          py::arg_v("options", M::EqnOptions(), "FFModel.EqnOptions()"),
          "Declare several equations on the same region.")
      .def("add_output", &add_output, py::arg("fct"),
           py::arg("doms")  = std::vector<mc::FFVar>(), py::kw_only(),
           py::arg("point") = std::nullopt, py::arg("side") = std::nullopt,
           py::arg("masks") = std::nullopt, py::arg("at") = std::nullopt, R"doc(
Declare an output.

- ``add_output(f)``: a scalar output (values and integrals built with
  ``FFEval`` / ``FFIntegral``);
- ``add_output(f, doms, point=[...], side=[...])``: ``f`` at a point of
  ``doms``, from the given sides;
- ``add_output(f, doms, masks=[...])``: ``f`` distributed over a region of
  ``doms``;
- ``add_output(f, doms, masks=[...], at={z: z0})``: distributed over ``doms``
  at fixed coordinates of others.

Parameters
----------
fct : FFVar
    The output expression.
doms : list of FFVar, optional
    Domains of the point or of the distribution.
point : list of float, optional
    Coordinates, one per domain (point output).
side : list of int, optional
    Side per domain for a point output: ``FFDom.MINUS`` (default) or
    ``FFDom.PLUS``.
masks : list of int, optional
    Region per domain (distributed output).
at : dict of FFVar to float, optional
    Fixed coordinates of further domains (distributed output).
)doc")
      .def(
          "add_transition",
          [](M& self, std::vector<mc::FFVar> const& l,
             std::vector<mc::FFVar> const& r, mc::FFVar const& d, double tau)
          { self.add_transition(l, r, d, tau); },
          py::arg("left"), py::arg("right"), py::arg("dom"), py::arg("tau"),
          R"doc(
Declare a transition at ``dom = tau``: the equations ``left(tau^-) =
right(tau^+)``, which determine the post-jump values of the states they
involve. States not mentioned are continuous.

Parameters
----------
left, right : list of FFVar
    Expressions evaluated before and after the jump.
dom : FFVar
    The evolution domain.
tau : float
    Location of the jump (interior to the domain).
)doc")
      .def(
          "add_transition",
          [](M& self, mc::FFVar const& l, mc::FFVar const& r,
             mc::FFVar const& d, double tau)
          { self.add_transition(l, r, d, tau); },
          py::arg("left"), py::arg("right"), py::arg("dom"), py::arg("tau"),
          "Scalar form of ``add_transition``.")
      .def(
          "add_transition",
          [](M& self, std::vector<mc::FFVar> const& l,
             std::vector<mc::FFVar> const& r) { self.add_transition(l, r); },
          py::arg("left"), py::arg("right"),
          "Transition given through one-sided evaluations: ``left`` evaluated "
          "at tau^- (``FFDom.MINUS``), ``right`` at tau^+ (``FFDom.PLUS``).")
      .def(
          "add_transition", [](M& self, mc::FFVar const& l, mc::FFVar const& r)
          { self.add_transition(l, r); }, py::arg("left"), py::arg("right"),
          "Scalar form of the one-sided ``add_transition``.")
      .def(
          "update_ref", [](M& self, mc::FFVar const& v, double val)
          { self.update_ref(v, val); }, py::arg("var"), py::arg("ref"),
          "Set the reference value of a state or input.")
      .def(
          "update_ref", [](M& self, mc::FFVar const& v, t_Fun const& f)
          { self.update_ref(v, f); }, py::arg("var"), py::arg("ref"),
          "Set the reference function ``ref({dom: coordinate}) -> float`` of a "
          "state or input.")
      .def(
          "register_control", [](M& self, mc::FFVar const& v)
          { self.register_control(v); }, py::arg("var"),
          "Register a declared input as a control: its DOFs become the "
          "sensitivity directions, in registration order.")
      .def(
          "clear_controls", [](M& self) { self.clear_controls(); },
          "Clear the control registry.")
      .def(
          "n_control_dof", [](M const& self) { return self.n_control_dof(); },
          "Total number of control DOFs.")
      .def(
          "control_dofs", [](M const& self, mc::FFVar const& u)
          { return self.control_dofs(u); }, py::arg("var"),
          "The DOF indices of control ``var``, in the order of its values.")
      .def(
          "fix_input",
          [](M& self, mc::FFVar const& w, std::vector<double> const& v)
          { self.fix_input(w, v); }, py::arg("var"), py::arg("values"),
          "Hold a declared input at known values: substituted as constants, "
          "never a parameter or a control.")
      .def(
          "setup", [](M& self) { return self.setup(); },
          py::call_guard<py::scoped_ostream_redirect,
                         py::scoped_estream_redirect>(),
          R"doc(
Validate and analyse the model; returns False on failure (see
``setup_status``).

The setup report (``options.DISPLAY_LEVEL`` >= 1) and any diagnostics, written
by C++ to std::cout / std::cerr, are redirected to Python's ``sys.stdout`` /
``sys.stderr`` for the duration of the call, so they appear in notebooks.
)doc")
      .def(
          "is_setup", [](M const& self) { return self.is_setup(); },
          "True once ``setup()`` has succeeded.")
      .def_property_readonly(
          "setup_status", [](M const& self) { return self.setup_status(); },
          "Outcome of the last ``setup()``.")
      .def_property_readonly(
          "var_domain", [](M const& self) { return self.var_domain(); },
          "Declared domain variables.")
      .def_property_readonly(
          "var_state", [](M const& self) { return self.var_state(); },
          "Declared states.")
      .def_property_readonly(
          "var_input", [](M const& self) { return self.var_input(); },
          "Inputs of the model (declared ones and those created by setup).")
      .def_property_readonly(
          "var_declared_input",
          [](M const& self) { return self.var_declared_input(); },
          "Inputs declared by the modeller.")
      .def_property_readonly(
          "var_constant", [](M const& self) { return self.var_constant(); },
          "Declared constants.")
      .def_property_readonly(
          "var_equation", [](M const& self) { return self.var_equation(); },
          "Declared equations, as ``FFModel.t_Eqn`` records.")
      .def_property_readonly(
          "var_output", [](M const& self) { return self.var_output(); },
          "Declared outputs, as ``FFModel.t_Fct`` records.")
      .def(
          "report",
          [](M const& self)
          {
            std::ostringstream os;
            self.report(os);
            py::print(os.str(), py::arg("end") = "");
          },
          "Print the model report (the same text as ``str(model)``).")
      .def(
          "__str__",
          [](M const& self)
          {
            std::ostringstream os;
            self.report(os);
            return os.str();
          },
          "Report of the model.");
}
